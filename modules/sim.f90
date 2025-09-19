module sim_m
    implicit none
    use koch_m

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
        real, intent(in) :: t0, delta_t
        real, intent(in), dimension(:) :: initial_state
        real, dimension(13) :: state, k1, k2, k3, k4

        ! DEFINE THE K TERMS FOR RK4 METHOD
        k1 = differential_equations_vector(t0, initial_state)
        k2 = differential_equations_vector(t0 + delta_t*0.5, initial_state + k1 * delta_t*0.5)
        k3 = differential_equations_vector(t0 + delta_t*0.5, initial_state + k2 * delta_t*0.5)
        k4 = differential_equations_vector(t0 + delta_t, initial_state + k3 * delta_t)

        ! DEFINE THE RESULT FROM RK4
        state = initial_state + delta_t/6 * (k1 + 2*k2 + 2*k3 + k4)

    end function rk4

  !=========================
  ! Equations of Motion: [u,v,w, p,q,r, x,y,z, e0,ex,ey,ez]
  !=========================
    function differential_equations(t, state) result(dy_dt)
        implicit none
        real, intent(in) :: t
        real, intent(in), dimension(:) :: state
        real, dimension(13) :: dy_dt
        real :: u, v, w, p, q, r, x, y_pos, z, e0, ex, ey, ez

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
        orientation_effect =  [2 * (ex*ez - ey*e0), &
                               2 * (ey*ez - ex*e0), & 
                               ez**2 + e0**2 - ex**2 - ey**2]
        angular_v_effect =    [q*w - r*v, r*u - p*w, p*v - q*u]
        gyroscopic_efffects = [0, -hzb, hyb, hzb, 0, -hxb, -hyb, hxb, 0]
        gyroscopic_change =   [hxb_dot, hyb_dot, dzb_dot]
        inertia_effects =     [(Iyy - Izz)*p*r + Iyz*(q**2 -r**2) + Ixz*p*q - Ixy*p*r, &
                               (Izz - Ixx)*p*r + Ixz*(r**2 -p**2) + Ixy*q*r - Iyz*p*q, &
                               (Ixx - Iyy)*p*r + Ixy*(p**2 -q**2) + Iyz*p*r - Ixz*q*r]
        wind_velocity = [Vx, Vy, Vz]
        quaternion = [e0, ex, ey, ez]
        velocity_quat = [0, u, v, w]
        quat_inv = [e0, -ex, -ey, -ez]
        quat_matrix = [-ex, -ey, -ez, e0, -ez, ey, ez, e0, -ex, -ey, ex, e0]

        ! DEFINE THE DIFFERENTIAL EQUATIONS
        ! ACCELERATION IN BODY FRAME
        acceleration = 1 / mass * FM(1:3) + g * orientation_effect + angular_v_effect

        ! ROLL, PITCH, YAW ACCELERATIONS
        angular_accelerations = matmul(inertia_inv , FM(4:6) + matmul(gyroscopic_efffects, [p, q, r]) + inertia_effects - gyroscopic_change)

        ! VELOCITY IN THE INERITAL FRAME
        velocity_quat_res = quat_mult(velocity_quat, quat_inv)
        velocity = quat_mult(quaterion, velocity_quat_res) + wind_velocity

        ! AIRCRAFT ORIENTATION RATE OF CHANGE
        quat_change = 0.5 * matmul(quat_matrix, [p, q, r])

        ! RETURN THE STATE DERIVATIVE
        dy_dt(1:3)  = acceleration
        dy_dt(4:6)  = angular_accelerations
        dy_dt(7:9)  = velocity(1:3)
        dy_dt(10:13) = quat_change
    end function differential_equations

  !=========================
  ! Aerodynamic Forces and Moments
  !=========================
    subroutine pseudo_aero()
        implicit none
    end subroutine pseudo_aero

  !=========================
  ! Mass and Inertia
  !=========================
    subroutine mass_inertia()
        implicit none
    end subroutine mass_inertia

    subroutine run 
        y = 0.0
