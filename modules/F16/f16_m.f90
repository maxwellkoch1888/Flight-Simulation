module f16_m
    use koch_m
    use json_m
    use jsonx_m
    implicit none

    ! BUILD GLOBAL VARIABLES FOR THE MODULE
    real :: mass
    real, dimension(3,3) :: inertia, inertia_inv
    real, dimension(6) :: FM
    integer :: io_unit

    
    ! DEFINE ARROW PROPERTIES
    real, parameter :: weight = 0.0697 !lbf
    real, parameter :: length = 2.3 !ft
    real, parameter :: planform_area = 0.000218 !ft^2
    logical :: print_steps = .true.

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

      ! PRINT STEPS IF SPECIFIED
      if(print_steps) then 
        write(io_unit,*) "diff_eq function called: "
        write(io_unit,*) "              RK4 call number =  1"
        write(io_unit,*) "diff_eq function called: "
        write(io_unit,*) "              RK4 call number =  2"
        write(io_unit,*) "diff_eq function called: "
        write(io_unit,*) "              RK4 call number =  3"
        write(io_unit,*) "diff_eq function called: "
        write(io_unit,*) "              RK4 call number =  4"
      end if

    end function rk4

  !=========================
  ! Equations of Motion: (/u,v,w, p,q,r, x,y,z, e0,ex,ey,ez/)

    function differential_equations(t, state) result(dstate_dt)
      implicit none 
      real, intent(in) :: t, state(13) 
      real :: dstate_dt(13) 
      real :: u, v, w, p, q, r, x, y, z, e0, ex, ey, ez 
      real :: acceleration(3), angular_accelerations(3), velocity(3), quat_change(4) 
      real :: velocity_quat(4), quat_inv(4) 
      real :: orientation_effect(3), angular_v_effect(3), gyroscopic_change(3), & 
              inertia_effects(3), wind_velocity(3), gyroscopic_effects(3,3) 
              real :: Ixx, Iyy, Izz, Ixy, Ixz, Iyz 
              real :: quat_matrix(4,3) 
      real :: hxb, hyb, hzb, hxb_dot, hyb_dot, hzb_dot, Vxw, Vyw, Vzw, & 
              gravity_ft_per_sec2 
      real :: avoid_warning, v1(4), v2(4)

      ! UNPACK STATES
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
      Ixx = inertia(1,1); Iyy = inertia(2,2); Izz = inertia(3,3)
      Ixy = inertia(1,2); Ixz = inertia(1,3); Iyz = inertia(2,3)
    
      ! CALCULATE FORCES AND MOMENTS
      call pseudo_aero(state)
      gravity_ft_per_sec2 = gravity_English(-state(9))

      avoid_warning = t

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
                              2 * (ey*ez + ex*e0), & 
                              ez**2 + e0**2 - ex**2 - ey**2/)

      angular_v_effect =    (/r*v - q*w, p*w - r*u, q*u - p*v/)

      gyroscopic_effects = reshape((/ 0.0, -hzb,  hyb, &
                                      hzb,  0.0, -hxb, &
                                     -hyb,  hxb,  0.0 /), (/3,3/))

      gyroscopic_change =   (/hxb_dot, hyb_dot, hzb_dot/)
      
      inertia_effects =     (/(Iyy - Izz)*q*r + Iyz*(q**2 -r**2) + Ixz*p*q - Ixy*p*r, &
                              (Izz - Ixx)*p*r + Ixz*(r**2 -p**2) + Ixy*q*r - Iyz*p*q, &
                              (Ixx - Iyy)*p*q + Ixy*(p**2 -q**2) + Iyz*p*r - Ixz*q*r/)
      wind_velocity = (/Vxw, Vyw, Vzw/)

      velocity_quat = (/0.0, u, v, w/)
      quat_inv = (/e0, -ex, -ey, -ez/)
      quat_matrix = reshape((/ -ex,  e0,  ez, -ey, &
                               -ey, -ez,  e0,  ex, &
                               -ez,  ey, -ex,  e0  /), (/4,3/))             

      ! BUILD THE DIFFERENTIAL EQUATIONS
      ! ACCELERATION IN BODY FRAME
      acceleration = 1.0 / mass * FM(1:3) + gravity_ft_per_sec2 * orientation_effect + angular_v_effect

      ! ROLL, PITCH, YAW ACCELERATIONS
      angular_accelerations = matmul(inertia_inv , FM(4:6) + &
                              matmul(gyroscopic_effects, (/p, q, r/)) + inertia_effects - gyroscopic_change)

      ! VELOCITY IN THE INERITAL FRAME
      v1 = quat_mult(state(10:13), (/0.0, state(1:3)/))
      v2 = quat_mult(v1, quat_inv)
      velocity = v2(2:4) + wind_velocity

      ! AIRCRAFT ORIENTATION RATE OF CHANGE
      quat_change = 0.5 * matmul(quat_matrix, (/p, q, r/))


      ! RETURN THE STATE DERIVATIVE
      dstate_dt(1:3)  = acceleration
      dstate_dt(4:6)  = angular_accelerations
      dstate_dt(7:9)  = velocity
      dstate_dt(10:13) = quat_change

      ! PRINT STEPS IF SPECIFIED
      if(print_steps) then
        write(io_unit,'(A,13(1X,ES20.12))') "                     time [s] =", t
        write(io_unit,'(A,13(1X,ES20.12))') '    State vector coming in    =', state
        write(io_unit,'(A,6 (1X,ES20.12))') "    Pseudo aerodynamics (F,M) =", FM
        write(io_unit,*) "v1", v1
        write(io_unit,*) "v2", v2
        write(io_unit,*) "velocity", velocity
        write(io_unit,*) "quat change", quat_change
        write(io_unit,'(A,13(1X,ES20.12))') "    Diff. Eq. results         =", dstate_dt
        write(io_unit,*) ''
      end if

      end function differential_equations

  !=========================
  ! Aerodynamic Forces and Moments for f16

    subroutine pseudo_aero(state)
      implicit none
      real, intent(in) :: state(13)
      real             :: Re, geometric_altitude_ft, geopotential_altitude_ft,     & 
                          temp_R, pressure_lbf_per_ft2, density_slugs_per_ft3, & 
                          dyn_viscosity_slug_per_ft_sec, sos_ft_per_sec
      real             :: V, Uc(3)
      real, parameter  :: CL_alpha = 4.929, CD0 = 5.096, CD2 = 48.138
      real, parameter  :: Cm_alpha = -2.605, Cm_qbar = -9.06, Cl_pitch_pbar = -5.378
      real             :: Cl_pitch0, alpha, beta, beta_f, pbar, qbar, rbar, angular_rates(3)
      real             :: CL, CD, CS, L, D, S, Cl_pitch, Cm, Cn
      real             :: ca, cb, sa, sb
        
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

      ! CALCULATE ALPHA_HAT USING EQN 3.4.20
      alpha_hat = dalpha_dtime * longitudinal_length / (2 * V)
      beta_ hat = dbeta_dtime * lateral_length / (2 * V)

      ! CALCULATE PBAR, QBAR, AND RBAR from eq 3.4.23
      angular_rates = 1 / (2*V) * (/state(4) * length, state(5) * length, state(6) * length/)
      pbar = angular_rates(1)
      qbar = angular_rates(2)
      rbar = angular_rates(3)

      ! CALCULATE THE REYNOLDS NUMBER
      Re = density_slugs_per_ft3 * V * 2 * length / dyn_viscosity_slug_per_ft_sec

      ! CALCULATE THE LIFT, DRAG, AND SIDE COEFFICIENTS

      CL1 =  CL0 + CL_alpha * alpha
      CL  = CL1 + CL_qbar * qbar + CL_alphahat * alpha_hat + CL_elevator * delta_e
      CS  = CS_beta * beta + (CS_pbar + CS_alpha_pbar * alpha) * pbar + CS_rbar * rbar &
           + CS_aileron * delta_a + CS_rudder * delta_r
      CD  =  CD_L0 + CD_L1 * CL1 + CD_L1_L1 * CL1 **2 + CD_CS_CS * CS **2 &
            + (CD_qbar + CD_alpha_qbar * alpha) * qbar + (CD_elevator + CD_alpha_elevator * alpha) &
            * delta_e + CD_elevator_elevator * delta_e ** 2

      L =   CL * (0.5 * density_slugs_per_ft3 * V **2 * planform_area)
      S =   CS * (0.5 * density_slugs_per_ft3 * V **2 * planform_area)
      D =   CD * (0.5 * density_slugs_per_ft3 * V **2 * planform_area)

      ca = cos(alpha)
      cb = cos(beta)
      sa = sin(alpha)
      sb = sin(beta)

      FM(1) = (D*ca*cb + S*ca*sb - L*sa) * (-1.0)
      FM(2) = (S*cb - D*sb)
      FM(3) = (D*sa*cb + S*sa*sb + L*ca) * (-1.0)

      ! CALCULATE THE ROLL, PITCH, AND YAW COEFFICIENTS
      Cl_pitch = Cl_beta * beta + Cl_pbar * pbar + (Cl_rbar + Cl_alpha_rbar * alpha) * rbar &
                 + Cl_aileron * delta_a + Cl_rudder * delta_r  ! roll moment
      Cm =       Cm_0 + Cm_alpha * alpha + Cm_qbar * qbar + Cm_alphahat * alpha_hat + Cm_elevator * delta_e ! pitch moment
      Cn =       Cn_beta * beta + (Cn_pbar + Cn_alpha_pbar * alpha) * pbar + Cn_rbar * rbar &
                 + (Cn_aileron + Cn_alpha_aileron * alpha) * delta_a + Cn_rudder * delta_r ! yaw moment

      FM(4) = Cl_pitch * (0.5 * density_slugs_per_ft3 * V **2 * planform_area * length)
      FM(5) = Cm       * (0.5 * density_slugs_per_ft3 * V **2 * planform_area * length)
      FM(6) = Cn       * (0.5 * density_slugs_per_ft3 * V **2 * planform_area * length)

      ! PRINT STEPS IF SPECIFIED
      if(print_steps) then
        write(io_unit,*) "alpha = ", alpha
        write(io_unit,*) "beta = ", beta
        write(io_unit,*) "beta_f = ", beta_f
        write(io_unit,*) "CL", CL
        write(io_unit,*) "CD", CD
        write(io_unit,*) "CS", CS
        write(io_unit,*) "L", L
        write(io_unit,*) "D", D
        write(io_unit,*) "S", S
        write(io_unit,*) "ca", ca
        write(io_unit,*) "cb", cb
        write(io_unit,*) "sa", sa
        write(io_unit,*) "sb", sb
      end if 

    end subroutine pseudo_aero

  !=========================
  ! Mass and Inertia

    subroutine mass_inertia(inertia)
      implicit none
      real, intent(inout) :: inertia(3,3)
      real :: I(3,3)

      ! DEFINE IDENTITY MATRIX 
      I = reshape((/ Ixx, Ixy, Ixz, &
             Ixy, Iyy, Iyz, &
             Ixz, Iyz, Izz /), (/3,3/))

      ! CALCULATE MASS AND INERTIA
      mass = weight / 32.17404855643
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
      real :: V, h, elevation_angle_deg
      character(len=40) :: filename
      type(json_value), pointer :: j_main

      ! UNPACK VALUES FROM THE JSON FILE
      call get_command_argument(1, filename)

      call jsonx_load(filename, j_main)

      call jsonx_get(j_main, 'simulation.time_step[s]',  dt)
      call jsonx_get(j_main, 'simulation.total_time[s]', tf)
      
      call jsonx_get(j_main, 'vehicle.CG_shift[ft]', CG_shift)

      call jsonx_get(j_main, 'vehicle.mass.weight[lbf]',     weight)
      call jsonx_get(j_main, 'vehicle.mass.Ixx[slug-ft^2]',  Ixx)
      call jsonx_get(j_main, 'vehicle.mass.Iyy[slug-ft^2]',  Iyy)
      call jsonx_get(j_main, 'vehicle.mass.Izz[slug-ft^2]',  Izz)
      call jsonx_get(j_main, 'vehicle.mass.Ixy[slug-ft^2]',  Ixy)
      call jsonx_get(j_main, 'vehicle.mass.Ixz[slug-ft^2]',  Ixz)
      call jsonx_get(j_main, 'vehicle.mass.Iyz[slug-ft^2]',  Iyz)
      call jsonx_get(j_main, 'vehicle.mass.hx[slug-ft^2/s]', hx)
      call jsonx_get(j_main, 'vehicle.mass.hy[slug-ft^2/s]', hy)
      call jsonx_get(j_main, 'vehicle.mass.hz[slug-ft^2/s]', hz)

      call jsonx_get(j_main, 'vehicle.thrust.T0[lbf]', T0)
      call jsonx_get(j_main, 'vehicle.thrust.Ta',      Ta)

      call jsonx_get(j_main, 'vehicle.aerodynamics.reference.area[ft^2]',              area)
      call jsonx_get(j_main, 'vehicle.aerodynamics.reference.longitudinal_length[ft]', longitudinal_length)
      call jsonx_get(j_main, 'vehicle.aerodynamics.reference.lateral_length[ft]',      lateral_length)

      call jsonx_get(j_main, 'vehicle.aerodynamics.coefficients.CL.0',        CL0)
      call jsonx_get(j_main, 'vehicle.aerodynamics.coefficients.CL.alpha',    CL_alpha)
      call jsonx_get(j_main, 'vehicle.aerodynamics.coefficients.CL.alphahat', CL_alphahat)
      call jsonx_get(j_main, 'vehicle.aerodynamics.coefficients.CL.qbar',     CL_qbar)
      call jsonx_get(j_main, 'vehicle.aerodynamics.coefficients.CL.elevator', CL_elevator)


      call jsonx_get(j_main, 'vehicle.aerodynamics.coefficients.CS.beta',       CS_beta)
      call jsonx_get(j_main, 'vehicle.aerodynamics.coefficients.CS.pbar',       CS_pbar)
      call jsonx_get(j_main, 'vehicle.aerodynamics.coefficients.CS.alpha_pbar', CS_alpha_pbar)
      call jsonx_get(j_main, 'vehicle.aerodynamics.coefficients.CS.rbar',       CS_rbar)
      call jsonx_get(j_main, 'vehicle.aerodynamics.coefficients.CS.aileron',    CS_aileron)
      call jsonx_get(j_main, 'vehicle.aerodynamics.coefficients.CS.rudder',     CS_rudder)


      call jsonx_get(j_main, 'vehicle.aerodynamics.coefficients.CD.L0',                CD_L0)
      call jsonx_get(j_main, 'vehicle.aerodynamics.coefficients.CD.CL1',               CD_L1)
      call jsonx_get(j_main, 'vehicle.aerodynamics.coefficients.CD.CL1_CL1',           CD_CL1_CL1)
      call jsonx_get(j_main, 'vehicle.aerodynamics.coefficients.CD.CS_CS',             CD_CS_CS)
      call jsonx_get(j_main, 'vehicle.aerodynamics.coefficients.CD.qbar',              CD_qbar)
      call jsonx_get(j_main, 'vehicle.aerodynamics.coefficients.CD.alpha_qbar',        CD_alpha_qbar)
      call jsonx_get(j_main, 'vehicle.aerodynamics.coefficients.CD.elevator',          CD_elevator)
      call jsonx_get(j_main, 'vehicle.aerodynamics.coefficients.CD.alpha_elevator',    CD_alpha_elevator)
      call jsonx_get(j_main, 'vehicle.aerodynamics.coefficients.CD.elevator_elevator', CD_elevator_elevator)

      call jsonx_get(j_main, 'vehicle.aerodynamics.coefficients.Cl.beta',       Cl_beta)
      call jsonx_get(j_main, 'vehicle.aerodynamics.coefficients.Cl.pbar',       Cl_pbar)
      call jsonx_get(j_main, 'vehicle.aerodynamics.coefficients.Cl.rbar',       Cl_rbar)
      call jsonx_get(j_main, 'vehicle.aerodynamics.coefficients.Cl.alpha_rbar', Cl_alpha_rbar)
      call jsonx_get(j_main, 'vehicle.aerodynamics.coefficients.Cl.aileron',    Cl_aileron)
      call jsonx_get(j_main, 'vehicle.aerodynamics.coefficients.Cl.rudder',     Cl_rudder)

      call jsonx_get(j_main, 'vehicle.aerodynamics.coefficients.Cm.0',        Cm_0)
      call jsonx_get(j_main, 'vehicle.aerodynamics.coefficients.Cm.alpha',    Cm_alpha)
      call jsonx_get(j_main, 'vehicle.aerodynamics.coefficients.Cm.qbar',     Cm_qbar)
      call jsonx_get(j_main, 'vehicle.aerodynamics.coefficients.Cm.alphahat', Cm_alphahat)
      call jsonx_get(j_main, 'vehicle.aerodynamics.coefficients.Cm.elevator', Cm_elevator)
    
      call jsonx_get(j_main, 'vehicle.aerodynamics.coefficients.Cn.beta',          Cn_beta)
      call jsonx_get(j_main, 'vehicle.aerodynamics.coefficients.Cn.pbar',          Cn_pbar)
      call jsonx_get(j_main, 'vehicle.aerodynamics.coefficients.Cn.alpha_pbar',    Cn_alpha_pbar)
      call jsonx_get(j_main, 'vehicle.aerodynamics.coefficients.Cn.rbar',          Cn_rbar)
      call jsonx_get(j_main, 'vehicle.aerodynamics.coefficients.Cn.aileron',       Cn_aileron)
      call jsonx_get(j_main, 'vehicle.aerodynamics.coefficients.Cn.alpha_aileron', Cn_alpha_aileron)
      call jsonx_get(j_main, 'vehicle.aerodynamics.coefficients.Cn.rudder',        Cn_rudder)
      ! (/1,2,3, 4,5,6, 7,8,9, 10,11,12,13/)
      ! (/u,v,w, p,q,r, x,y,z, e0,ex,ey,ez/)
      call jsonx_get(j_main, 'initial.airspeed[ft/s]',       airspeed)
      call jsonx_get(j_main, 'initial.altitude[ft]',         initial_state(9))
      call jsonx_get(j_main, 'initial.elevation_angle[deg]', eul(2))
      call jsonx_get(j_main, 'initial.bank_angle[deg]',      eul(1))
      call jsonx_get(j_main, 'initial.heading_angle[deg]',   eul(3))
      call jsonx_get(j_main, 'initial.alpha[deg]',           alpha_deg)
      call jsonx_get(j_main, 'initial.beta[deg]',            beta_deg)
      call jsonx_get(j_main, 'initial.p[deg/s]',             initial_state(4))
      call jsonx_get(j_main, 'initial.q[deg/s]',             initial_state(5))
      call jsonx_get(j_main, 'initial.r[deg/s]',             initial_state(6))
      call jsonx_get(j_main, 'initial.aileron[deg]',         aileron_deg)
      call jsonx_get(j_main, 'initial.elevator[deg]',        elevator_deg)
      call jsonx_get(j_main, 'initial.rudder[deg]',          rudder_deg)
      call jsonx_get(j_main, 'initial.throttle',             throttle_deg)

      ! OPEN A FILE TO WRITE TO 
      open(newunit=io_unit, file='arrow_output_straight.txt', status='replace', action='write')

      ! INITIALIZE TIME
      t = 0.0

      ! BUILD INITIAL CONDITIONS
      initial_state = 0.0
      initial_state(1) = V !ft/sec
      eul = eul * pi / 180.0
      initial_state(10:13) = euler_to_quat(eul)

      ! CALCULATE MASS AND INERTIA
      call mass_inertia(inertia)     

      ! BUIL.0 ! zero deg orientationD THE LOOP AND WRITE THE OUTPUT
      write(io_unit,*) " time[s]             u[ft/s]             v[ft/s]&
                   w[ft/s]             p[rad/s]            q[rad/s]      &
            r[rad/s]            xf[ft]              yf[ft]              &
      zf[ft]              e0                  ex                  ey     &
                   ez                  "
      write(io_unit,'(14ES20.12)') t,initial_state(:)
    
      do while(t < tf)
        if (print_steps) write(io_unit,*) ''

        ! CALCULATE THE NEW STATE
        new_state = rk4(t, initial_state, dt)

        ! NORMALIZE THE QUATERNION
        call quat_norm(new_state(10:13))

        ! UPDATE THE STATE AND TIME
        initial_state = new_state
        t = t + dt

        if (print_steps) write(io_unit,*) ''

        write(io_unit,'(14ES20.12)') t,initial_state(:)
  
        if (print_steps) write(io_unit,*) ''

      end do 

      close(io_unit)

    end subroutine run
end module f16_m