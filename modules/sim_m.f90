module sim_m
    use koch_m
    implicit none

    ! BUILD GLOBAL VARIABLES FOR THE MODULE
    real :: mass
    real, dimension(3,3) :: inertia, inertia_inv
    real, dimension(6) :: FM

    contains
  !=========================
  ! RK4 Integrator
  !=========================
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
  !=========================
    function differential_equations(t, state) result(dstate_dt)
      implicit none
      real, intent(in) :: t
      real, intent(in), dimension(13) :: state
      real :: dstate_dt(13)
      real :: acceleration(3), angular_accelerations(3), velocity(3), quat_change(4)
      real :: quaternion(4), velocity_quat(4), quat_inv(4)
      real :: orientation_effect(3), angular_v_effect(3), gyroscopic_change(3), inertia_effects(3)
      real :: wind_velocity(3)
      real :: u, v, w, p, q, r, x_pos, y_pos, z, e0, ex, ey, ez
      real :: Ixx, Iyy, Izz, Ixy, Ixz, Iyz
      real :: quat_matrix(4,3)
      real :: temp_vec3(3)
      real :: gyroscopic_effects(3,3)

        call pseudo_aero

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

        ! BUILD MATRICES/ VECTORS USED IN DIFFERENTIAL EQUATION
        orientation_effect =  (/2 * (ex*ez - ey*e0), &
                               2 * (ey*ez - ex*e0), & 
                               ez**2 + e0**2 - ex**2 - ey**2/)

        angular_v_effect =    (/q*w - r*v, r*u - p*w, p*v - q*u/)

        gyroscopic_effects = reshape((/ 0.0, -hzb,  hyb, &
                                      hzb,  0.0, -hxb, &
                                      -hyb,  hxb,  0.0 /), (/3,3/))

        gyroscopic_change =   (/hxb_dot, hyb_dot, dzb_dot/)
        
        inertia_effects =     (/(Iyy - Izz)*p*r + Iyz*(q**2 -r**2) + Ixz*p*q - Ixy*p*r, &
                               (Izz - Ixx)*p*r + Ixz*(r**2 -p**2) + Ixy*q*r - Iyz*p*q, &
                               (Ixx - Iyy)*p*r + Ixy*(p**2 -q**2) + Iyz*p*r - Ixz*q*r/)
        wind_velocity = (/Vx, Vy, Vz/)

        quaternion = (/e0, ex, ey, ez/)
        velocity_quat = (/0, u, v, w/)
        quat_inv = (/e0, -ex, -ey, -ez/)
        quat_matrix = reshape((/ -ex, -ey, -ez, &
                         e0, -ez,  ey, &
                         ez,  e0, -ex, &
                        -ey,  ex,  e0 /), (/4,3/))

        ! BUILD THE DIFFERENTIAL EQUATIONS
        ! ACCELERATION IN BODY FRAME
        acceleration = 1 / mass * FM(1:3) + g * orientation_effect + angular_v_effect

        ! ROLL, PITCH, YAW ACCELERATIONS
        angular_accelerations = matmul(inertia_inv , FM(4:6) + matmul(gyroscopic_efffects, (/p, q, r/)) + inertia_effects - gyroscopic_change)

        ! VELOCITY IN THE INERITAL FRAME
        velocity = quat_mult(quaterion, quat_mult(velocity_quat, quat_inv)) + wind_velocity

        ! AIRCRAFT ORIENTATION RATE OF CHANGE
        quat_change = 0.5 * matmul(quat_matrix, (/p, q, r/))

        ! RETURN THE STATE DERIVATIVE
        dstate_dt(1:3)  = acceleration
        dstate_dt(4:6)  = angular_accelerations
        dstate_dt(7:9)  = velocity
        dstate_dt(10:13) = quat_change
    end function differential_equations

  !=========================
  ! Aerodynamic Forces and Moments
  !=========================
    subroutine pseudo_aero()
        implicit none

        ! CALCULATE THE REYNOLDS NUMBER
        Re = density * V * 2 * radius / dyn_viscosity_pa_sec

        ! COMPUTE THE DRAG COEFFICIENT
        if(0.0 < Re .and. Re <= 450000.0) then
          CD = 24.0/Re + 6.0/(1.0 + Re**0.5) + 0.4
        else if(450000.0 < Re .and. Re <= 560000.0) then 
          CD = 1.0E29 * Re **(-5.211)
        else if(560000 < Re .and. Re <= 14000000) then
          CD = -2.0E-23 * Re**3 - 1.0E-16*Re**2 + 9.0E-09 * Re + 0.069
        else if(Re > 1400000) then
          CD = 0.12
        end if

        FM = 0
        FM(1) = -0.5 * density * V**2 * pi * radius**2 * CD * Uc
    end subroutine pseudo_aero

  !=========================
  ! Mass and Inertia
  !=========================
    subroutine mass_inertia()
        implicit none
    end subroutine mass_inertia

  !=========================
  ! Matrix Inverse
  !=========================
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
  !=========================
    subroutine run()
      implicit none
      real :: t, dt, tf, initial_state(13), new_state(13), eul(3)

      ! INITIALIZE TIME
      t = 0.0
      dt = 0.01
      tf = 10.0

      ! BUILD INITIAL CONDITIONS
      initial_state = 0.0
      initial_state(1) = 50.0 !ft/sec
      initial_state(9) = -200 !altitude in ft
      eul = 0.0 ! zero deg orientation
      initial_state(10:13) = euler_to_quat(eul)


      ! CALCULATE MASS AND INERTIA
      call mass_inertia

      ! BUILD THE LOOP AND WRITE THE OUTPUT
      write(*,*) " t(sec)              u(ft/sec)           v(ft/sec)       w(ft/sec)         "
      write(*,'(14ES17.9)') t,initial_state(:)
    
      do while(t<tf)
        ! CALCULATE THE NEW STATE
        new_state = rk4(t, initial_state, dt)

        ! NORMALIZE THE QUATERNION
        call quat_norm(new_state(10:13))

        ! UPDATE THE STATE AND TIME
        initial_state = new_state
        t = t + dt
      write(*,'(14ES17.9)') t,initial_state(:)
      end do 
    end subroutine run