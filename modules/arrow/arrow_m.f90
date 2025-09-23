module arrow_m
    use koch_m
    ! use jsonx
    implicit none

    ! BUILD GLOBAL VARIABLES FOR THE MODULE
    real :: mass
    real, dimension(3,3) :: inertia, inertia_inv
    real, dimension(6) :: FM
    
    ! DEFINE ARROW PROPERTIES
    real, parameter :: weight = 0.0697 !lbf
    real, parameter :: length = 2.3 !ft
    real, parameter :: S = 0.000218 !ft^2
    logical :: straight = .true.

    contains
  !=========================
  ! RK4 Integrator

    function rk4(t0, initial_state, delta_t) result(state)
        implicit none
        real, intent(in) :: t0, delta_t, initial_state(13)
        real, dimension(13) :: state, k1, k2, k3, k4

        ! DEFINE THE K TERMS FOR RK4 METHOD
        k1 = differential_equations(t0, initial_state)
        k2 = differential_equations(t0 + delta_t*0.5, initial_state + k1 * delta_t*0.5)
        k3 = differential_equations(t0 + delta_t*0.5, initial_state + k2 * delta_t*0.5)
        k4 = differential_equations(t0 + delta_t, initial_state + k3 * delta_t)

        ! DEFINE THE RESULT FROM RK4
        state = initial_state + (delta_t/6) * (k1 + 2*k2 + 2*k3 + k4)

    end function rk4

  !=========================
  ! Equations of Motion: (/u,v,w, p,q,r, x,y,z, e0,ex,ey,ez/)

    function differential_equations(t, state) result(dstate_dt)
      implicit none
      real, intent(in) :: t, state(13)
      real :: dstate_dt(13)
      real :: u, v, w, p, q, r, x, y, z, e0, ex, ey, ez
      real :: acceleration(3), angular_accelerations(3), velocity(3), quat_change(4)
      real :: quaternion(4), velocity_quat(4), quat_inv(4)
      real :: orientation_effect(3), angular_v_effect(3), gyroscopic_change(3), &
              inertia_effects(3), wind_velocity(3), gyroscopic_effects(3,3)
      real :: Ixx, Iyy, Izz, Ixy, Ixz, Iyz
      real :: quat_matrix(4,3)
      real :: hxb, hyb, hzb, hxb_dot, hyb_dot, hzb_dot, Vxw, Vyw, Vzw, &
              gravity_ft_per_sec2
      real :: avoid_warning


        call pseudo_aero(state)
        gravity_ft_per_sec2 = gravity_English(-state(9))

        avoid_warning = t

        ! UNPACK STATE VARIABLES
        u  = state(1)
        v  = state(2)
        w  = state(3)
        p  = state(4)
        q  = state(5)
        r  = state(6)
        x  = state(7)
        y  = state(8)   
        z  = state(9)
        e0 = state(10)
        ex = state(11)
        ey = state(12)
        ez = state(13)

        ! UNPACK INERTIA
        Ixx = inertia(1, 1)
        Iyy = inertia(2, 2)
        Izz = inertia(3, 3)
        Ixy = inertia(1, 2)
        Ixz = inertia(1, 3)
        Iyz = inertia(2, 3)

      ! SET GYROSCOPIC EFFECTS AND WIND VELOCITY TO ZERO
      hxb = 0.0
      hyb = 0.0
      hzb = 0.0
      hxb_dot = 0.0
      hyb_dot = 0.0
      hzb_dot = 0.0
      Vxw = 0.0
      Vyw = 0.0
      Vzw = 0.0

        ! BUILD MATRICES/ VECTORS USED IN DIFFERENTIAL EQUATION
        orientation_effect =  (/2 * (ex*ez - ey*e0), &
                               2 * (ey*ez - ex*e0), & 
                               ez**2 + e0**2 - ex**2 - ey**2/)

        angular_v_effect =    (/q*w - r*v, r*u - p*w, p*v - q*u/)

        gyroscopic_effects = reshape((/ 0.0, -hzb,  hyb, &
                                      hzb,  0.0, -hxb, &
                                      -hyb,  hxb,  0.0 /), (/3,3/))

        gyroscopic_change =   (/hxb_dot, hyb_dot, hzb_dot/)
        
        inertia_effects =     (/(Iyy - Izz)*p*r + Iyz*(q**2 -r**2) + Ixz*p*q - Ixy*p*r, &
                               (Izz - Ixx)*p*r + Ixz*(r**2 -p**2) + Ixy*q*r - Iyz*p*q, &
                               (Ixx - Iyy)*p*r + Ixy*(p**2 -q**2) + Iyz*p*r - Ixz*q*r/)
        wind_velocity = (/Vxw, Vyw, Vzw/)

        quaternion = (/e0, ex, ey, ez/)
        velocity_quat = (/0.0, u, v, w/)
        quat_inv = (/e0, -ex, -ey, -ez/)
        quat_matrix = reshape((/ -ex, -ey, -ez, &
                         e0, -ez,  ey, &
                         ez,  e0, -ex, &
                        -ey,  ex,  e0 /), (/4,3/))

        ! BUILD THE DIFFERENTIAL EQUATIONS
        ! ACCELERATION IN BODY FRAME
        acceleration = 1 / mass * FM(1:3) + gravity_ft_per_sec2 * orientation_effect + angular_v_effect

        ! ROLL, PITCH, YAW ACCELERATIONS
        angular_accelerations = matmul(inertia_inv , FM(4:6) + &
                                matmul(gyroscopic_effects, (/p, q, r/)) + inertia_effects - gyroscopic_change)

        ! VELOCITY IN THE INERITAL FRAME
        velocity = quat_base_to_dependent((/u, v, w/), quaternion) + wind_velocity

        ! AIRCRAFT ORIENTATION RATE OF CHANGE
        quat_change = 0.5 * matmul(quat_matrix, (/p, q, r/))

        ! RETURN THE STATE DERIVATIVE
        dstate_dt(1:3)  = acceleration
        dstate_dt(4:6)  = angular_accelerations
        dstate_dt(7:9)  = velocity
        dstate_dt(10:13) = quat_change
    end function differential_equations

  !=========================
  ! Aerodynamic Forces and Moments for straight fletchings

    subroutine pseudo_aero(state)
      implicit none
      real, intent(in) :: state(13)
      real :: Re, geometric_altitude_ft, geopotential_altitude_ft,     & 
        temp_R, pressure_lbf_per_ft2, density_slugs_per_ft3, & 
        dyn_viscosity_slug_per_ft_sec, sos_ft_per_sec
      real :: V, Uc(3)
      real, parameter :: CL_alpha = 4.929, CD0 = 5.096, CD2 = 48.138
      real, parameter :: Cm_alpha = -2.605, Cm_qbar = -9.06, Cl_pitch_pbar = -5.378
      real :: Cl_pitch0, alpha, beta, beta_f, pbar, qbar, rbar, angular_rates(3)
      real :: CL, CD, CS, L, D, S_force, Cl_pitch, Cm, Cn
      real :: ca, cb, sa, sb
        
      ! BUILD THE ATMOSPHERE 
      geometric_altitude_ft = -state(9)
      call std_atm_English(&
        geometric_altitude_ft, geopotential_altitude_ft,     & 
        temp_R, pressure_lbf_per_ft2, density_slugs_per_ft3, & 
        dyn_viscosity_slug_per_ft_sec, sos_ft_per_sec)

      ! CALCULATE VELOCITY UNIT VECTOR
      V =  (state(1)**2 + state(2)**2 + state(3)**2)**0.5
      Uc = (/state(1:3)/) / V

      ! CALCULATE ALPHA AND BETA 3.4.4 and 3.4.5
      alpha =  atan2(state(3) , state(1))
      beta =   asin(state(2) / V)
      beta_f = atan2(state(2) , state(1))


      ! CALCULATE PBAR, QBAR, AND RBAR from eq 3.4.23
      angular_rates = 1 / (2*V) * (/state(4) * length, state(5) * length, state(6) * length/)
      pbar = angular_rates(1)
      qbar = angular_rates(2)
      rbar = angular_rates(3)

      ! CALCULATE THE REYNOLDS NUMBER
      Re = density_slugs_per_ft3 * V * 2 * length / dyn_viscosity_slug_per_ft_sec

      ! CALCULATE THE LIFT, DRAG, AND SIDE COEFFICIENTS
      CL =      CL_alpha * alpha
      CS =     -CL_alpha * beta_f
      CD =      CD0 + CD2 * CL**2 + CD2 * CS**2
      L =       CL * (density_slugs_per_ft3 * V **2 * S)
      S_force = CS * (density_slugs_per_ft3 * V **2 * S)
      D =       CD * (density_slugs_per_ft3 * V **2 * S)

      ca = cos(alpha)
      cb = cos(beta)
      sa = sin(alpha)
      sb = sin(beta)

      FM(3) = (D*ca*cb + S*ca*sb - L*sa) * (-1.0)
      FM(2) = (S*cb - D*sb)
      FM(1) = (D*sa*cb + S*sa*sb + L*ca) * (-1.0)

      ! CALCULATE THE ROLL, PITCH, AND YAW COEFFICIENTS
      Cl_pitch0 = 0.0

      Cl_pitch = Cl_pitch0 + Cl_pitch_pbar * pbar  ! roll moment
      Cm =       Cm_alpha * alpha + Cm_qbar * qbar ! pitch moment
      Cn =      -Cm_alpha * beta_f + Cm_qbar * rbar ! yaw moment

      FM(4) = Cl_pitch * 0.5 * density_slugs_per_ft3 * V **2 * S * length
      FM(5) = Cm       * 0.5 * density_slugs_per_ft3 * V **2 * S * length
      FM(6) = Cn       * 0.5 * density_slugs_per_ft3 * V **2 * S * length

    end subroutine pseudo_aero

  !=========================
  ! Mass and Inertia

    subroutine mass_inertia(inertia)
        implicit none
        real, intent(inout) :: inertia(3,3)
        real :: I(3,3)

        ! DEFINE IDENTITY MATRIX 
        I = reshape((/ 0.0, 0.0, 0.0, &
               0.0, 0.0, 0.0, &
               0.0, 0.0, 0.0 /), (/3,3/))

        ! CALCULATE MASS AND INERTIA
        mass = weight / 32.17404855643
        inertia(1,1) = 0.0000194
        inertia(2,2) = 0.00097
        inertia(3,3) = 0.00097
        inertia_inv = matrix_inv(inertia)
    end subroutine mass_inertia

  !=========================
  ! Matrix Inverse

    function matrix_inv(A) result(A_inv)
      implicit none
      real, dimension(3,3), intent(in) :: A
      real, dimension(3,3) :: A_inv, A_adj
      real :: det

      ! CALCULATE THE DETERMINANT
      det = A(1,1) * A(2,2) * A(3,3) + A(1,2) * A(2,3) * A(3,1) + A(1,3) * A(2,1) * A(3,2) &
          - A(1,3) * A(2,2) * A(3,1) - A(1,2) * A(2,1) * A(3,3) - A(1,1) * A(2,3) * A(3,2)

      ! CALCULATE THE ADJUGATE MATRIX
      A_adj(1,1) = A(2,2) * A(3,3) - A(2,3) * A(3,2)
      A_adj(1,2) = A(1,3) * A(3,2) - A(1,2) * A(3,3)
      A_adj(1,3) = A(1,2) * A(2,3) - A(1,3) * A(2,2)

      A_adj(2,1) = A(2,3) * A(3,1) - A(2,1) * A(3,3)
      A_adj(2,2) = A(1,1) * A(3,3) - A(1,3) * A(3,1)
      A_adj(2,3) = A(1,3) * A(2,1) - A(1,1) * A(2,3)

      A_adj(3,1) = A(2,1) * A(3,2) - A(2,2) * A(3,1)
      A_adj(3,2) = A(1,2) * A(3,1) - A(1,1) * A(3,2)
      A_adj(3,3) = A(1,1) * A(2,2) - A(1,2) * A(2,1)

      ! CALCULATE THE INVERSE MATRIX
      A_inv = 1/det * A_adj
    end function matrix_inv

  !=========================
  ! Run Subroutine

    subroutine run()
      implicit none
      real :: t, dt, tf, initial_state(13), new_state(13), eul(3)
      integer :: io_unit
      ! character(len=*) :: filename
      ! type(json_type) :: j_main
      real :: V, h, elevation_angle_deg

      ! OPEN A FILE TO WRITE TO 
      open(newunit=io_unit, file='arrow_output.txt', status='replace', action='write')
      write(io_unit,*) 'Arrow Test'

      ! ! UNPACK VALUES FROM THE JSON FILE
      ! call get_command_argument(1, filename)

      ! call jsonx_load(filename, j_main)

      ! call jsonx_get(j_main, 'simulation.time_step[s]', dt)
      ! call jsonx_get(j_main, 'initial.airspeed[ft/s]', V)
      ! call jsonx_get(j_main, 'initial.altitude[ft]', h)
      ! h = h * -1
      ! call jsonx_get(j_main, 'initial.elevation_angle[deg]', elevation_angle_deg)

      dt = 0.01
      V = 210.0
      h = -5.0
      elevation_angle_deg = 5.0

      ! INITIALIZE TIME
      t = 0.0
      tf = 10.0

      ! BUILD INITIAL CONDITIONS
      initial_state = 0.0
      initial_state(1) = V !ft/sec
      initial_state(9) = h !altitude in ft
      eul = 0.0 ! zero deg orientation
      eul(2) = elevation_angle_deg * pi / 180.0
      initial_state(10:13) = euler_to_quat(eul)


      ! CALCULATE MASS AND INERTIA
      call mass_inertia(inertia)     

      ! BUILD THE LOOP AND WRITE THE OUTPUT
      write(io_unit,*) &
     "           t(sec)              u(ft/sec)               v(ft/sec)               w(ft/sec)              ", &
     "p(rad/sec)            q(rad/sec)            r(rad/sec)               x(ft)                   y(ft)                   ", &
     "z(ft)                 e0                     ex                    ey                    ez"
      write(io_unit,'(14ES23.12)') t,initial_state(:)
    
      do while(initial_state(9) <= 0)
        ! CALCULATE THE NEW STATE
        new_state = rk4(t, initial_state, dt)

        ! NORMALIZE THE QUATERNION
        call quat_norm(new_state(10:13))

        ! UPDATE THE STATE AND TIME
        initial_state = new_state
        t = t + dt
      write(io_unit,'(14ES23.12)') t,initial_state(:)
      end do 
      close(io_unit)

    end subroutine run
end module arrow_m