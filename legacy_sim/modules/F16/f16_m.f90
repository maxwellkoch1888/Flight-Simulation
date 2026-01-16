module f16_m
  use koch_m
  use jsonx_m
  use micro_time_m
  use linalg_mod
  use connection_m

  ! GLOBAL VARIABLES
    implicit none
  
    ! JSON POINTER
    type(json_value), pointer :: j_main

    ! UDP POINTERS
    type(connection) :: graphics
    type(connection) :: user_controls

    ! BUILD GLOBAL VARIABLES FOR THE MODULE
    real :: mass
    real :: inertia(3,3)
    real :: inertia_inv(3,3)
    real, target :: h(3)
    real :: FM(6)
    real :: initial_state(13)
    real :: controls(4)
    real :: rho0
    real :: T0, Ta
    real, allocatable :: aero_ref_location(:)
    integer :: io_unit
    real :: init_airspeed 

    ! CONTROLLER VARIABLES 
    real :: Amat(6,6), Bmat(6,4), xd(13), ud(4) 

    ! DEFINE AERODYNAMIC PROPERTIES
    real :: planform_area, longitudinal_length, lateral_length, sweep
    real :: CL0, CL_alpha, CL_alphahat, CL_qbar, CL_elevator
    real :: CS_beta, CS_pbar, CS_alpha_pbar, CS_rbar, CS_aileron, CS_rudder
    real :: CD_L0, CD_L1, CD_L1_L1, CD_CS_CS, CD_qbar, CD_alpha_qbar, CD_elevator, CD_alpha_elevator, CD_elevator_elevator
    real :: Cl_beta, Cl_pbar, Cl_rbar, Cl_alpha_rbar, Cl_aileron, Cl_rudder
    real :: Cm_0, Cm_alpha, Cm_qbar, Cm_alphahat, Cm_elevator
    real :: Cn_beta, Cn_pbar, Cn_alpha_pbar, Cn_rbar, Cn_aileron, Cn_alpha_aileron, Cn_rudder
    real :: CL_lambda_b, CL_alpha_0, CL_alpha_s, CD_lambda_b, CD_alpha_0, CD_alpha_s, Cm_lambda_b
    real :: Cm_alpha_0, Cm_alpha_s, Cm_min
    logical :: compressibility, rk4_verbose, print_states, stall, test_stall, test_compressibility
    logical :: print_stall

    ! ADD VARIABLES FOR TRIM ALGORITHM
    character(:), allocatable :: sim_type
    character(:), allocatable :: trim_type
    real :: relaxation_factor, tolerance, max_iterations, finite_difference_step_size
    real :: bank_angle0, sideslip_angle0, climb_angle0, elevation_angle0
    logical :: trim_verbose, exam_answers

  
  contains
  ! INTEGRATOR AND EQN OF MOTION
    !=========================
    ! RK4 Integrator
      function rk4(t0, y1, delta_t) result(state)
        implicit none
        real, intent(in) :: t0, delta_t, y1(13)
        real, dimension(13) :: state, k1, k2, k3, k4

        ! DEFINE THE K TERMS FOR RK4 METHOD
        k1 = differential_equations(t0, y1)
        k2 = differential_equations(t0 + delta_t*0.5, y1 + k1 * delta_t*0.5)
        k3 = differential_equations(t0 + delta_t*0.5, y1 + k2 * delta_t*0.5)
        k4 = differential_equations(t0 + delta_t, y1 + k3 * delta_t)

        ! DEFINE THE RESULT FROM RK4
        state = y1 + (delta_t/6) * (k1 + 2*k2 + 2*k3 + k4)

      end function rk4

    !=========================
    ! Equations of Motion: (/u,v,w, p,q,r, x,y,z, e0,ex,ey,ez/)
      function differential_equations(t, state) result(dstate_dt)
        implicit none 
        real, intent(in) :: t
        real, target :: state(13) 
        real :: dstate_dt(13) 
        real :: acceleration(3), angular_accelerations(3), rhs(3), velocity(3), quat_change(4) 
        real :: quat_inv(4) 
        real :: orientation_effect(3), angular_v_effect(3), gyroscopic_change(3)
        real :: inertia_effects(3), wind_velocity(3), gyroscopic_effects(3,3) 
        real :: quat_matrix(4,3) 
        real :: hxb_dot, hyb_dot, hzb_dot
        real :: Vxw, Vyw, Vzw, gravity_ft_per_sec2 
        real :: angular_inertia(3,3), angular_inertia_inv(3,3)
        real :: avoid_warning
        real :: Ixx, Iyy, Izz, Ixy, Ixz, Iyz 
        real, pointer :: u, v, w, p, q, r, e0, ex, ey, ez
        real, pointer :: hxb, hyb, hzb

        avoid_warning = t
        
        ! UNPACK STATES
        u  => state(1)
        v  => state(2)
        w  => state(3)
        p  => state(4)
        q  => state(5)
        r  => state(6)
        e0 => state(10)
        ex => state(11)
        ey => state(12)
        ez => state(13) 


        ! UNPACK INERTIA
        Ixx = inertia(1,1)
        Iyy = inertia(2,2)
        Izz = inertia(3,3)
        Ixy = inertia(1,2)
        Ixz = inertia(1,3)
        Iyz = inertia(2,3)
      
        ! CALCULATE FORCES AND MOMENTS
        call pseudo_aero(state)
        gravity_ft_per_sec2 = gravity_English(-state(9))

        ! SET GYROSCOPIC EFFECTS
        hxb => h(1)
        hyb => h(2)
        hzb => h(3)

        ! SET GYROSCOPIC CHANGE AND WIND VELOCITY TO ZERO
        hxb_dot = 0.0
        hyb_dot = 0.0
        hzb_dot = 0.0
        Vxw = 0.0
        Vyw = 0.0
        Vzw = 0.0

        ! BUILD MATRICES/ VECTORS USED IN DIFFERENTIAL EQUATION
        orientation_effect(1) = 2 * (ex*ez - ey*e0)
        orientation_effect(2) = 2 * (ey*ez + ex*e0)
        orientation_effect(3) = ez**2 + e0**2 - ex**2 - ey**2

        angular_v_effect(1) = r*v - q*w 
        angular_v_effect(2) = p*w - r*u
        angular_v_effect(3) = q*u - p*v

        gyroscopic_effects(1,1) =  0.0
        gyroscopic_effects(1,2) =  hzb
        gyroscopic_effects(1,3) = -hyb

        gyroscopic_effects(2,1) = -hzb
        gyroscopic_effects(2,2) =  0.0
        gyroscopic_effects(2,3) =  hxb

        gyroscopic_effects(3,1) =  hyb
        gyroscopic_effects(3,2) = -hxb
        gyroscopic_effects(3,3) =  0.0

        angular_inertia(1,1) =  Ixx
        angular_inertia(1,2) = -Ixy
        angular_inertia(1,3) = -Ixz

        angular_inertia(2,1) = -Ixy
        angular_inertia(2,2) =  Iyy
        angular_inertia(2,3) = -Iyz

        angular_inertia(3,1) = -Ixz
        angular_inertia(3,2) = -Iyz
        angular_inertia(3,3) =  Izz

        angular_inertia_inv = matrix_inv(angular_inertia)

        gyroscopic_change(1) = hxb_dot
        gyroscopic_change(2) = hyb_dot
        gyroscopic_change(3) = hzb_dot

        inertia_effects(1) = (Iyy - Izz)*q*r + Iyz*(q**2 -r**2) + Ixz*p*q - Ixy*p*r
        inertia_effects(2) = (Izz - Ixx)*p*r + Ixz*(r**2 -p**2) + Ixy*q*r - Iyz*p*q
        inertia_effects(3) = (Ixx - Iyy)*p*q + Ixy*(p**2 -q**2) + Iyz*p*r - Ixz*q*r

        wind_velocity(1) = Vxw
        wind_velocity(2) = Vyw
        wind_velocity(3) = Vzw

        quat_inv(1) =  e0
        quat_inv(2) = -ex
        quat_inv(3) = -ey
        quat_inv(4) = -ez

        quat_matrix(1,1) = -ex
        quat_matrix(1,2) =  e0
        quat_matrix(1,3) =  ez

        quat_matrix(2,1) = -ey
        quat_matrix(2,2) = -ez
        quat_matrix(2,3) =  e0

        quat_matrix(3,1) = -ez
        quat_matrix(3,2) =  ey
        quat_matrix(3,3) = -ex

        quat_matrix(4,1) = -ey
        quat_matrix(4,2) =  ex
        quat_matrix(4,3) =  e0            

        ! BUILD THE DIFFERENTIAL EQUATIONS
        ! ACCELERATION IN BODY FRAME
        acceleration(1) = FM(1)/mass + gravity_ft_per_sec2*orientation_effect(1) + angular_v_effect(1)
        acceleration(2) = FM(2)/mass + gravity_ft_per_sec2*orientation_effect(2) + angular_v_effect(2)
        acceleration(3) = FM(3)/mass + gravity_ft_per_sec2*orientation_effect(3) + angular_v_effect(3)

        ! ROLL, PITCH, YAW ACCELERATIONS
        rhs(1) = FM(4) + gyroscopic_effects(1,1)*p + gyroscopic_effects(1,2)*q + gyroscopic_effects(1,3)*r + inertia_effects(1) - gyroscopic_change(1)
        rhs(2) = FM(5) + gyroscopic_effects(2,1)*p + gyroscopic_effects(2,2)*q + gyroscopic_effects(2,3)*r + inertia_effects(2) - gyroscopic_change(2)
        rhs(3) = FM(6) + gyroscopic_effects(3,1)*p + gyroscopic_effects(3,2)*q + gyroscopic_effects(3,3)*r + inertia_effects(3) - gyroscopic_change(3)

        angular_accelerations(1) = angular_inertia_inv(1,1)*rhs(1) + angular_inertia_inv(1,2)*rhs(2) + angular_inertia_inv(1,3)*rhs(3)
        angular_accelerations(2) = angular_inertia_inv(2,1)*rhs(1) + angular_inertia_inv(2,2)*rhs(2) + angular_inertia_inv(2,3)*rhs(3)
        angular_accelerations(3) = angular_inertia_inv(3,1)*rhs(1) + angular_inertia_inv(3,2)*rhs(2) + angular_inertia_inv(3,3)*rhs(3)

        ! VELOCITY IN THE INERITAL FRAME
        velocity = quat_base_to_dependent(state(1:3), state(10:13)) + wind_velocity

        ! AIRCRAFT ORIENTATION RATE OF CHANGE
        quat_change(1) = 0.5 * (quat_matrix(1,1)*p + quat_matrix(1,2)*q + quat_matrix(1,3)*r)
        quat_change(2) = 0.5 * (quat_matrix(2,1)*p + quat_matrix(2,2)*q + quat_matrix(2,3)*r)
        quat_change(3) = 0.5 * (quat_matrix(3,1)*p + quat_matrix(3,2)*q + quat_matrix(3,3)*r)
        quat_change(4) = 0.5 * (quat_matrix(4,1)*p + quat_matrix(4,2)*q + quat_matrix(4,3)*r)

        ! RETURN THE STATE DERIVATIVE
        dstate_dt(1:3)  = acceleration
        dstate_dt(4:6)  = angular_accelerations
        dstate_dt(7:9)  = velocity
        dstate_dt(10:13) = quat_change
      end function differential_equations
  ! 
  ! AERODYNAMICS AND FORCES
    !=========================
    ! Aerodynamic Forces and Moments for f16
      subroutine pseudo_aero(state)
        implicit none
        real, intent(in) :: state(13)
        real :: Re, geometric_altitude_ft, geopotential_altitude_ft
        real :: temp_R, pressure_lbf_per_ft2, density_slugs_per_ft3
        real :: dyn_viscosity_slug_per_ft_sec, sos_ft_per_sec
        real :: V, dyn_pressure
        real :: Cl_roll0, alpha, beta, beta_f, pbar, qbar, rbar, angular_rates(3)
        real :: CL, CL1, CD, CS, L, D, S, Cl_roll, Cm, Cn
        real :: CM1, CM2, mach_num, gamma, R
        real :: ca, cb, sa, sb
        real :: alpha_hat, beta_hat
        real :: delta_a, delta_e, delta_r
        real :: T, throttle
        real :: CL_ss, CD_ss, Cm_ss, CL_s, CD_s, Cm_s 
        real :: sigma_D, sigma_L, sigma_m, sign_a
        
        ! COMPRESSIBILITY
        real :: blend, drag_factor, lift_factor, moment_factor, drag_rise
        real :: pg_factor, kt_factor
        real :: m_high, m_low, m_trans_start, m_trans_end
        real :: k, var, sqrt_term
        real :: tran, max_CL_factor, max_CD_factor, max_Cl_roll_factor, transonic_blend
          
        ! BUILD THE ATMOSPHERE 
        geometric_altitude_ft = -state(9)
        call std_atm_English(&
          geometric_altitude_ft, geopotential_altitude_ft,     & 
          temp_R, pressure_lbf_per_ft2, density_slugs_per_ft3, & 
          dyn_viscosity_slug_per_ft_sec, sos_ft_per_sec)

        ! CALCULATE VELOCITY UNIT VECTOR
        V =  (state(1)**2 + state(2)**2 + state(3)**2)**0.5
        dyn_pressure = 0.5 * density_slugs_per_ft3 * V **2 * planform_area
        ! CALCULATE ALPHA AND BETA 3.4.4 and 3.4.5
        alpha =  atan2(state(3) , state(1))
        beta =   asin(state(2) / V)
        beta_f = atan2(state(2) , state(1))

        ca = cos(alpha)
        cb = cos(beta)
        sa = sin(alpha)
        sb = sin(beta)

        ! CALCULATE ALPHA_HAT USING EQN 3.4.20
        alpha_hat = 0.0
        beta_hat = 0.0

        ! CALCULATE PBAR, QBAR, AND RBAR from eq 3.4.23
        angular_rates(1) = 1 / (2*V) * state(4) * lateral_length
        angular_rates(2) = 1 / (2*V) * state(5) * longitudinal_length
        angular_rates(3) = 1 / (2*V) * state(6) * lateral_length

        pbar = angular_rates(1)
        qbar = angular_rates(2)
        rbar = angular_rates(3)

        ! CALCULATE THE REYNOLDS NUMBER
        ! Re = density_slugs_per_ft3 * V * 2 * longitudinal_length / dyn_viscosity_slug_per_ft_sec

        ! PULL OUT CONTROLS
        delta_a = controls(1)
        delta_e = controls(2)
        delta_r = controls(3)
        throttle = controls(4)

        ! CALCULATE THE LIFT, DRAG, AND SIDE FORCE COEFFICIENTS
        CL1 =  CL0 + CL_alpha * alpha
        CL_ss  = CL1 + CL_qbar * qbar + CL_alphahat * alpha_hat + CL_elevator * delta_e
        CS  = CS_beta * beta + (CS_pbar + CS_alpha_pbar * alpha) * pbar + CS_rbar * rbar &
            + CS_aileron * delta_a + CS_rudder * delta_r
        CD_ss  =  CD_L0 + CD_L1 * CL1 + CD_L1_L1 * CL1 **2 + CD_CS_CS * CS **2 &
              + (CD_qbar + CD_alpha_qbar * alpha) * qbar + (CD_elevator + CD_alpha_elevator * alpha) &
              * delta_e + CD_elevator_elevator * delta_e ** 2

        ! CALCULATE THE ROLL, PITCH, AND YAW COEFFICIENTS
        Cl_roll = Cl_beta * beta + Cl_pbar * pbar + (Cl_rbar + Cl_alpha_rbar * alpha) * rbar &
                  + Cl_aileron * delta_a + Cl_rudder * delta_r  ! roll moment
        Cm_ss =    Cm_0 + Cm_alpha * alpha + Cm_qbar * qbar + Cm_alphahat * alpha_hat + Cm_elevator * delta_e ! pitch moment
        Cn =       Cn_beta * beta + (Cn_pbar + Cn_alpha_pbar * alpha) * pbar + Cn_rbar * rbar &
                  + (Cn_aileron + Cn_alpha_aileron * alpha) * delta_a + Cn_rudder * delta_r ! yaw moment

        ! STALL MODEL FOR FORCES
        if (stall) then 
          sign_a = sign(1.0,alpha)
          CL_s = 2 * sign_a * sa**2 * ca 
          CD_s = 2 * (sin(abs(alpha)))**3
          sigma_L = calc_sigma(CL_lambda_b, CL_alpha_0, CL_alpha_s, alpha)
          sigma_D = calc_sigma(CD_lambda_b, CD_alpha_0, CD_alpha_s, alpha)
          
          CL = CL_ss * (1 - sigma_L) + CL_s * sigma_L 
          CD = CD_ss * (1 - sigma_D) + CD_s * sigma_D 

          Cm_s = CM_min * sa**2 * sign_a
          sigma_m = calc_sigma(Cm_lambda_b, Cm_alpha_0, Cm_alpha_s, alpha)
          Cm = Cm_ss * (1 - sigma_m) + Cm_s * sigma_m 
     
        end if 

        ! COMPRESSIBILIITY MODEL, USING PRANDTL-GLAUERT CORRECTION
        if (compressibility) then 
          CM1 = 2.13/ (sweep + 0.15)**2
          CM2 = 15.35*sweep**2 - 19.64*sweep +16.86
          mach_num = V / sos_ft_per_sec

          ! MACH BREAKPOINTS
          m_low         = 0.60    ! only use prandtl-glauert
          m_high        = 0.92    ! only use karman-tsien
          m_trans_start = 0.88    ! start transsonic region

          ! Prandtl-Glauert factor
          sqrt_term = sqrt(max(1.0 - mach_num**2, tol))
          pg_factor = 1.0 / sqrt_term

          ! Karman-Tsien factor
          kt_factor = pg_factor * (1.0 + mach_num**2 / (1.0 + sqrt_term))

          ! BLEND BETWEEN PG AND KT
          if (mach_num <= m_low) then
            blend = 0.0
          else if (mach_num >= m_high) then
            blend = 1.0
          else
            blend = (mach_num - m_low) / (m_high - m_low)
          end if

          ! FINAL FACTOR USING BLEND
          lift_factor = (1.0 - blend) * pg_factor + blend * kt_factor
          moment_factor = lift_factor
          drag_factor = 1.0 + CM1 * mach_num**CM2

          ! APPLY THE FACTORS
          max_CL_factor = 6.5
          max_CD_factor = 4
          max_Cl_roll_factor = 2.5

          if (mach_num < 0.92) then ! accurate range
            CL = CL * lift_factor
            CS = CS * lift_factor
            Cl_roll = Cl_roll * moment_factor
            Cm       = Cm * moment_factor
            Cn       = Cn * moment_factor
            CD       = CD * drag_factor
          else
            CL       = CL * min(lift_factor, max_CL_factor)
            CS       = CS * min(lift_factor, max_CL_factor)
            Cl_roll = Cl_roll * min(moment_factor, max_Cl_roll_factor)
            Cm       = Cm * min(moment_factor, max_Cl_roll_factor)
            Cn       = Cn * min(moment_factor, max_Cl_roll_factor)
          end if

          ! apply drag multiplier
          CD = CD * min(drag_factor, max_CD_factor)
        end if 

        L =   CL * dyn_pressure
        S =   CS * dyn_pressure
        D =   CD * dyn_pressure

        ! TABLE 3.4.4
        FM(1) = - (ca*(D*cb + S*sb) - L*sa)
        FM(2) = (S*cb - D*sb)
        FM(3) = - (sa*(D*cb + S*sb) + L*ca) * (-1.0)

        FM(4) = Cl_roll * dyn_pressure * lateral_length
        FM(5) = Cm       * dyn_pressure * longitudinal_length
        FM(6) = Cn       * dyn_pressure * lateral_length

        ! ADD THE ENGINE THRUST
        if (throttle < 0.0) then
            throttle = 0.0
          else if (throttle > 1.0) then
            throttle = 1.0
          else 
            throttle = throttle 
        end if 

        T = throttle * T0 * (density_slugs_per_ft3/rho0) ** Ta
        FM(1) = FM(1) + T

        ! SHIFT CG LOCATION
        FM(4:6) = FM(4:6) + cross_product(aero_ref_location, FM(1:3))
      end subroutine pseudo_aero

    !=========================
    ! Check stall blending funtion
      subroutine check_stall(state) 
        implicit none 
        integer :: i 
        real :: state(13)
        real :: alpha, beta, states(13), N, Y, A
        real :: ca, cb, sa, sb
        real :: CL, CD, Cm
        real :: const, geometric_altitude_ft, geopotential_altitude_ft, temp_R, pressure_lbf_per_ft2, density_slugs_per_ft3, dyn_viscosity_slug_per_ft_sec, sos_ft_per_sec

        geometric_altitude_ft = -initial_state(9)
        call std_atm_English(&
          geometric_altitude_ft, geopotential_altitude_ft,     & 
          temp_R, pressure_lbf_per_ft2, density_slugs_per_ft3, & 
          dyn_viscosity_slug_per_ft_sec, sos_ft_per_sec)

        const = (0.5 * density_slugs_per_ft3 * init_airspeed**2 * planform_area)
        write(io_unit,*) 'altitude', geometric_altitude_ft
        write(io_unit,*) 'rho', density_slugs_per_ft3 
        write(io_unit,*) 'airspeed', init_airspeed 
        write(io_unit,*) 'planform area', planform_area

        controls = 0.0 
        states = 0.0 

        write(io_unit,*) '  alpha[deg]                CL                        CD                        Cm'  
        do i=-180, 180, 1 
          alpha = real(i) * pi / 180.0
          beta = 0.0

          ca = cos(alpha) 
          cb = cos(beta) 
          sa = sin(alpha) 
          sb = sin(beta) 

          states(1) = init_airspeed *ca * cb 
          states(2) = init_airspeed * sb 
          states(3) = init_airspeed * sa * cb 
          states(9) = initial_state(9)

          call pseudo_aero(states) 
          A = -FM(1) 
          Y =  FM(2)
          N = -FM(3) 
          
          CL = N * ca - A* sa
          CD = A * ca * cb - Y * sb + N * sa * cb 
          Cm = FM(5) 

          CL = CL / const 
          CD = CD / const 
          Cm = Cm / const / longitudinal_length
          write(io_unit,*) alpha*180.0/pi, CL, CD, Cm
        end do 
      end subroutine

    !=========================
    ! Check Compressiblity model
      subroutine check_compressibility(state)
        implicit none

        integer :: i, j
        real :: state(13)
        real :: alpha, beta, states(13)
        real :: N, Y, A
        real :: ca, cb, sa, sb
        real :: CL, CD, Cm
        real :: const
        real :: geometric_altitude_ft, geopotential_altitude_ft
        real :: temp_R, pressure_lbf_per_ft2, density_slugs_per_ft3
        real :: dyn_viscosity_slug_per_ft_sec, sos_ft_per_sec
        real :: airspeed, mach_num
        real, dimension(9) :: mach_num_list

        ! Mach numbers to sweep
        mach_num_list = (/ 0.3, 0.6, 0.7, 0.8, 0.9, 0.95, 1.0, 1.05, 1.1 /)

        ! Get altitude for density
        geometric_altitude_ft = -initial_state(9)
        call std_atm_English( &
          geometric_altitude_ft, geopotential_altitude_ft, &
          temp_R, pressure_lbf_per_ft2, density_slugs_per_ft3, &
          dyn_viscosity_slug_per_ft_sec, sos_ft_per_sec)

        controls = 0.0
        states   = 0.0

        ! Main Mach loop -----------------------------------------------
        do j = 1, size(mach_num_list)
          mach_num = mach_num_list(j)

          write(io_unit,*) '--------------------------------------------------'
          write(io_unit,*) 'Mach Number = ', mach_num
          write(io_unit,*) '--------------------------------------------------'

          airspeed = mach_num * sos_ft_per_sec

          const = 0.5 * density_slugs_per_ft3 * airspeed**2 * planform_area

          write(io_unit,*) '   alpha(deg)                CL                        CD                        Cm'

          do i = -180, 180
            alpha = real(i) * pi / 180.0
            beta  = 0.0

            ca = cos(alpha)
            cb = cos(beta)
            sa = sin(alpha)
            sb = sin(beta)

            states = 0.0
            states(1) = airspeed * ca * cb
            states(2) = airspeed * sb
            states(3) = airspeed * sa * cb
            states(9) = initial_state(9)

            call pseudo_aero(states)

            A = -FM(1)
            Y =  FM(2)
            N = -FM(3)

            ! Aerodynamic coefficients
            CL = (N * ca - A * sa) / const
            CD = (A * ca * cb - Y * sb + N * sa * cb) / const
            Cm = FM(5) / (const * longitudinal_length)

            write(io_unit,*) alpha*180.0/pi, CL, CD, Cm
          end do

        end do

      end subroutine
       
    !=========================
    ! Calculate Sigma
      function calc_sigma(lambda_b, alpha_0, alpha_s, alpha) result(sigma)
        implicit none
        real, intent(in) :: lambda_b, alpha_0, alpha_s, alpha
        real :: sigma, neg, pos

        pos = exp( lambda_b * (alpha - alpha_0 + alpha_s))
        neg = exp(-lambda_b * (alpha - alpha_0 - alpha_s))
        sigma = (1.0 + neg + pos) / ((1.0 + neg)*(1.0 + pos))

      end function calc_sigma

    !=========================
    ! Mass and Inertia
      subroutine mass_inertia()
        implicit none
        real :: weight
        real :: Ixx, Iyy, Izz, Ixy, Ixz, Iyz

        ! READ CG SHIFT IF APPLICABLE
        call jsonx_get(j_main, 'vehicle.aero_ref_location[ft]', aero_ref_location, 0.0, 3)

        ! READ MASS PROPERTIES
        call jsonx_get(j_main, 'vehicle.mass.weight[lbf]',     weight)
        call jsonx_get(j_main, 'vehicle.mass.Ixx[slug-ft^2]',  Ixx)
        call jsonx_get(j_main, 'vehicle.mass.Iyy[slug-ft^2]',  Iyy)
        call jsonx_get(j_main, 'vehicle.mass.Izz[slug-ft^2]',  Izz)
        call jsonx_get(j_main, 'vehicle.mass.Ixy[slug-ft^2]',  Ixy)
        call jsonx_get(j_main, 'vehicle.mass.Ixz[slug-ft^2]',  Ixz)
        call jsonx_get(j_main, 'vehicle.mass.Iyz[slug-ft^2]',  Iyz)
        call jsonx_get(j_main, 'vehicle.mass.hx[slug-ft^2/s]', h(1))
        call jsonx_get(j_main, 'vehicle.mass.hy[slug-ft^2/s]', h(2))
        call jsonx_get(j_main, 'vehicle.mass.hz[slug-ft^2/s]', h(3))

        ! DEFINE IDENTITY MATRIX 
        inertia(1,1) = Ixx
        inertia(1,2) = Ixy
        inertia(1,3) = Ixz
        inertia(2,1) = Ixy
        inertia(2,2) = Iyy
        inertia(2,3) = Iyz
        inertia(3,1) = Ixz
        inertia(3,2) = Iyz
        inertia(3,3) = Izz

        ! CALCULATE MASS AND INERTIA
        mass = weight / 32.17404855643
        inertia_inv = matrix_inv(inertia)
      end subroutine mass_inertia
  ! 
  ! CALCULATE TRIM STATE
    !=========================
    ! Trim Algorithm
      function trim_algorithm(V_mag, height, euler, tolerance, trim_type) result(G)
        implicit none 
        real :: V_mag, height, tolerance
        real :: alpha, beta, p, q, r, delta_a, delta_e, delta_r, throttle 
        real :: elevation_angle, azimuth_angle
        real :: c_bank, c_elev, s_bank, s_elev, ca, cb, sa, sb, error, pw
        real :: u, v, w, velocities(3), gravity
        real :: G(6), res(6), iteration_residual(6)
        real :: angular_rates(3), euler(3), print_statement(13)
        real :: cgamma, sgamma, climb_angle, solution, theta1, theta2, solution1
        integer :: k, iteration, case_number
        character(*), intent(in) :: trim_type
        
        ! DETERMINE THE TRIM TYPE BEING SPECIFIED
        if (trim_type == 'sct') then 
          if (elevation_angle0 /= -999.0) then 
            case_number = 1 ! sct, elevation_angle, bank_angle
          else
            case_number = 2 ! sct, climb_angle, bank_angle
          end if 
        else 
          if (elevation_angle0 /= -999.0) then 
            if (sideslip_angle0 == -999.0) then 
              case_number = 3 ! shss, elevation_angle, bank_angle 
            else 
              case_number = 4 ! shss, elevation_angle, sideslip_angle
            end if 
          else 
            if (sideslip_angle0 == -999.0) then 
              case_number = 5 ! shss, climb_angle, bank_angle
            else 
              case_number = 6 ! shss, climb_angle, sideslip_angle
            end if 
          end if 
        end if 
        
        if (trim_verbose) then 
          write(io_unit,*) 'Case Number:', case_number
          write(io_unit,*) 'Initial Airspeed:', V_mag 
          write(io_unit,*) 'Initial Height', height 
        end if 

        ! CALCULATE GRAVITY
        gravity = gravity_English(-height)

        ! SET INITIAL GUESSSES TO ZERO
        G        = 0.0
        p        = 0.0
        q        = 0.0
        r        = 0.0
        alpha    = G(1)
        beta     = G(2)
        delta_a  = G(3)
        delta_e  = G(4)
        delta_r  = G(5)
        throttle = G(6)

        if (trim_verbose) then
          write(io_unit,*) 'Trimming Aircraft for ', trim_type
          write(io_unit,'(A,ES20.12)') '  --> Azimuth angle set to psi [deg] =', euler(3) * 180 / pi
          if (elevation_angle0 /= -999.0) then 
            write(io_unit,'(A,ES20.12)') '  --> Elevation angle set to theta [deg] =', elevation_angle0 * 180 / pi
          else 
            write(io_unit,'(A,ES20.12)') '  --> Elevation angle set to theta [deg] =', elevation_angle0
          end if 
          if (sideslip_angle0 /= -999.0) then 
            write(io_unit,'(A,ES20.12)') '  --> Sideslip angle set to beta [deg] =', sideslip_angle0 * 180 / pi 
          else 
            write(io_unit,'(A,ES20.12)') '  --> Bank angle set to phi [deg] =', bank_angle0 * 180 / pi
          end if
          write(io_unit,'(A)') ''
          if (elevation_angle0 /= -999.0) then 
            write(io_unit,'(A,ES20.12)') 'Initial theta [deg] =', elevation_angle0 * 180 / pi 
            write(io_unit,'(A,ES20.12)') 'Initial gamma [deg] =', climb_angle0 
          else 
            write(io_unit,'(A,ES20.12)') 'Initial theta [deg] =', elevation_angle0
            write(io_unit,'(A,ES20.12)') 'Initial gamma [deg] =', climb_angle0 * 180 / pi 
          end if 
          write(io_unit,'(A,ES20.12)') 'Initial phi [deg]   =', bank_angle0 * 180 / pi
          write(io_unit,'(A,ES20.12)') 'Initial beta [deg]  =', beta * 180 / pi 
          write(io_unit,'(A)') ''
          write(io_unit,'(A)') 'Newton Solver Settings:'
          write(io_unit,'(A,ES20.12)') 'Finite Difference Step Size =', finite_difference_step_size
          write(io_unit,'(A,ES20.12)') '          Relaxation Factor =', relaxation_factor
          write(io_unit,'(A,ES20.12)') '                  Tolerance =', tolerance
          write(io_unit,'(A)') ''
          write(io_unit,'(A)') ''
        end if

        ! PULL OUT AZIMUTH, BANK, AND ELEVATION ANGLES
        azimuth_angle = euler(3)
        
        if (any(case_number == [4,6])) then 
          beta = sideslip_angle0
          euler(1) = 0.0
        else 
          euler(1) = bank_angle0
        end if 

        if (any(case_number == [1,3,4])) then 
          euler(2) = elevation_angle0
        end if 

        ! LOOP FOR FINDING TRIM STATE
        error = 1.0
        iteration = 1
        do while (error > tolerance)
          if (iteration > max_iterations) then 
            write(*,*) 'No trim solutions... Reached max iterations...'
          end if     
          ! DEFINE COS AND SIN TERMS TO SAVE TIME
          alpha = G(1)
          if (any(case_number == [1,2,3,5])) then
            beta = G(2)
          else 
            euler(1) = G(2)
          end if 

          ca = cos(alpha)
          cb = cos(beta) 
          sa = sin(alpha) 
          sb = sin(beta)
          c_bank = cos(euler(1))
          s_bank = sin(euler(1))


          ! CALCULATE VELOCITIES FROM 3.4.12
          velocities = V_mag * (/ca*cb, sb, sa*cb/) 

          u = velocities(1)
          v = velocities(2)
          w = velocities(3)

          ! CALCULATE ELEVATION ANGLE IF CLIMB ANGLE SPECIFIED
          if (any(case_number == [2,5,6])) then 
            if (trim_verbose) then 
              write(io_unit,*) ' '            
              write(io_unit,*) 'Solving for elevation angle given a climb angle:'
            end if 
            climb_angle = climb_angle0
            cgamma = cos(climb_angle)
            sgamma = sin(climb_angle)

              theta1 = asin((u*V_mag*sgamma + (v*s_bank + w*c_bank) * & 
                sqrt(u**2 + (v*s_bank + w*c_bank)**2 - V_mag**2*sgamma**2)) &
                / (u**2 + (v*s_bank + w*c_bank)**2))

              theta2 = asin((u*V_mag*sgamma - (v*s_bank + w*c_bank) * & 
                        sqrt(u**2 + (v*s_bank + w*c_bank)**2 - V_mag**2*sgamma**2)) &
                        / (u**2 + (v*s_bank + w*c_bank)**2))      
              
              solution1 = u*sin(theta1) - (v*s_bank + w*c_bank)*cos(theta1)
              solution = V_mag * sgamma 

              if (abs(solution1-solution) < tol) then 
                euler(2) = theta1 
              else 
                euler(2) = theta2 
              end if 

              if (trim_verbose) then
                write(io_unit,*) '        theta 1 [deg] =', theta1 * 180 / pi 
                write(io_unit,*) '        theta 2 [deg] =', theta2 * 180 / pi 
                write(io_unit,*) '  Correct theta [deg] =', euler(2) * 180 / pi 
                write(io_unit,*) '  Correct theta [rad] =', euler(2) 
                write(io_unit,*) ' '              
              end if 
          end if 

          s_elev = sin(euler(2))
          c_elev = cos(euler(2)) 

          ! CALCULATE THE ANGULAR RATES
          if (any(case_number == [1,2])) then
            angular_rates = (gravity * s_bank * c_elev) / (u*c_elev*c_bank + w*s_elev) &
                            * (/-s_elev, s_bank*c_elev, c_bank*c_elev/)
          else if (any(case_number == [3,4,5,6])) then
            angular_rates = 0.0

          end if 
          p = angular_rates(1)
          q = angular_rates(2)
          r = angular_rates(3)

          if (trim_verbose) then 
            write(io_unit,*) 'Updating rotation rates for ', trim_type
            write(io_unit,'(A,ES20.12)') 'p [deg/s] = ', p * 180 / pi 
            write(io_unit,'(A,ES20.12)') 'q [deg/s] = ', q * 180 / pi 
            write(io_unit,'(A,ES20.12)') 'r [deg/s] = ', r * 180 / pi
            write(io_unit, '(A)') ''
          end if 
          res = calc_r(V_mag, height, euler, angular_rates, G)
          
          ! STATE THE DEFINITION OF G
          if (any(case_number == [4,6])) then 
            if (trim_verbose) then 
              write(io_unit,'(A)') 'G defined as G = [alpha, bank_angle, aileron, elevator, rudder, throttle]'
              write(io_unit, '(A,6(1X,ES20.12))') '      G =', (G(k), k=1,6)
              write(io_unit, '(A,6(1X,ES20.12))') '      R =', (res(k), k=1,6)
              write(io_unit, '(A)') ''
            end if
          else 
            if (trim_verbose) then 
              write(io_unit,'(A)') 'G defined as G = [alpha, beta, aileron, elevator, rudder, throttle]'
              write(io_unit, '(A,6(1X,ES20.12))') '      G =', (G(k), k=1,6)
              write(io_unit, '(A,6(1X,ES20.12))') '      R =', (res(k), k=1,6)
              write(io_unit, '(A)') ''
            end if
          end if 
          
          ! SOLVE FOR G 
          call newtons_method(V_mag, height, euler, angular_rates, G)

          ! CALCULATE THE ITERATION ERROR
          iteration_residual = calc_r(V_mag, height, euler, angular_rates, G)
          error = maxval(abs(iteration_residual))

          if (trim_verbose) then 
            write(io_unit,'(A)') 'New G:'
            write(io_unit,'(A,6(1X,ES20.12))') '      G =', (G(k),   k=1,6)
            write(io_unit,'(A,6(1X,ES20.12))') '      R =', (iteration_residual(k), k=1,6) 
            write(io_unit,*) ''
            write(io_unit, '(A)') &
              'Iteration   Residual           alpha[deg]           beta[deg]            '// & 
              'p[deg/s]             q[deg/s]             r[deg/s]   ' // &
              '          phi[deg]             theta[deg]           aileron[deg]         elevator[deg]   ' // &
              '     rudder[deg]          throttle[]'
            if (sideslip_angle0 /= -999.0) then 
              write(io_unit,'(I6,1X,12(1X,ES20.12))') iteration, error, G(1)*180/pi, sideslip_angle0 * 180 / pi, &
                p * 180 / pi, q * 180 / pi, r * 180 / pi, euler(1) * 180 / pi, euler(2) * 180 / pi, G(3)*180/pi, &
                G(4)*180/pi, G(5)*180/pi, G(6)
            else 
              write(io_unit,'(I6,1X,12(1X,ES20.12))') iteration, error, G(1)*180/pi, G(2)*180/pi, &
                p * 180 / pi, q * 180 / pi, r * 180 / pi, euler(1) * 180 / pi, euler(2) * 180 / pi, G(3)*180/pi, &
                G(4)*180/pi, G(5)*180/pi, G(6)
            end if 
          end if 

          ! ENSURE THROTTLE IS IN BOUNDS
          if (G(6) > 1.0) then 
            if (trim_verbose) then 
              write(io_unit,*) 'Overwriting throttle > 1.'
            G(6) = 1.0
            end if 
          else if (G(6) < 0.0) then 
            if (trim_verbose) then 
              write(io_unit,*) 'Overwriting throttle < 0.'
              write(io_unit,*) ''
            end if 
            G(6) = 0.0
          end if 

          iteration = iteration + 1
        end do

        ! SAVE THE TRIM STATE AS THE NEW INITIAL CONDITIONS
        ! SET CONTROLS
        controls(1:4) = G(3:6)

        ! PULL OUT ALPHA
        alpha = G(1)
        ca = cos(alpha)
        sa = sin(alpha)     

        ! PULL OUT BETA
        if (any(case_number == [4,6])) then 
          cb = cos(sideslip_angle0)
          sb = sin(sideslip_angle0)
          euler(1) = G(2)
        else 
          cb = cos(G(2))
          sb = sin(G(2))
        end if

        ! CALCUALTE INITIAL STATES
        initial_state(1:3)   = V_mag * (/ca*cb, sb, sa*cb/) 
        initial_state(4:6)   = angular_rates
        initial_state(9)     = height
        initial_state(10:13) = euler_to_quat(euler)    

        ! UPDATE DESIRED STATE FOR CONTROLLER
        xd = initial_state
        ud = controls

        if (trim_verbose) then 
          write(io_unit,*) '---------------------- Trim Solution ----------------------'
          write(io_unit,'(A30,F20.13)') '       azimuth_angle[deg]  :', euler(3)
          write(io_unit,'(A30,F20.13)') '       elevation_angle[deg]:', euler(2)
          write(io_unit,'(A30,F20.13)') '       bank_angle[deg]     :', euler(1) 
          write(io_unit,'(A30,F20.13)') '       alpha[deg]          :', alpha 
          write(io_unit,'(A30,F20.13)') '       beta[deg]           :', beta 
          write(io_unit,'(A30,F20.13)') '       p[deg]              :', p 
          write(io_unit,'(A30,F20.13)') '       q[deg]              :', q
          write(io_unit,'(A30,F20.13)') '       r[deg]              :', r
          write(io_unit,'(A30,F20.13)') '       p_w[deg]            :', initial_state(4)
          write(io_unit,'(A30,F20.13)') '       q_w[deg]            :', initial_state(5)
          write(io_unit,'(A30,F20.13)') '       r_w[deg]            :', initial_state(6)
          write(io_unit,'(A30,F20.13)') '       aileron[deg]        :', controls(1)
          write(io_unit,'(A30,F20.13)') '       elevator[deg]       :', controls(2)
          write(io_unit,'(A30,F20.13)') '       rudder[deg]         :', controls(3)
          write(io_unit,'(A30,F20.13)') '       throttle[deg]       :', controls(4)
          write(io_unit,'(A30,F20.13)') '       Climb Angle[deg]    :', climb_angle
        end if 

        if (exam_answers) then 
          write(io_unit,*) '---------------------- Exam Answers ----------------------'
          write(io_unit,'(A30,ES25.13E3)') '       elevation_angle[deg]:', euler(2)
          write(io_unit,'(A30,ES25.13E3)') '       bank_angle[deg]     :', euler(1)
          write(io_unit,'(A30,ES25.13E3)') '       alpha[deg]          :', alpha
          write(io_unit,'(A30,ES25.13E3)') '       beta[deg]           :', beta
          write(io_unit,'(A30,ES25.13E3)') '       p[deg]              :', p
          write(io_unit,'(A30,ES25.13E3)') '       q[deg]              :', q
          write(io_unit,'(A30,ES25.13E3)') '       r[deg]              :', r
          write(io_unit,'(A30,ES25.13E3)') '       aileron[deg]        :', controls(1)
          write(io_unit,'(A30,ES25.13E3)') '       elevator[deg]       :', controls(2)
          write(io_unit,'(A30,ES25.13E3)') '       rudder[deg]         :', controls(3)
          write(io_unit,'(A30,ES25.13E3)') '       throttle[deg]       :', controls(4)

        end if 
      end function trim_algorithm

    !=========================
    ! Newton's Method Solver to find G (alpha, beta, delta_a, delta_e, delta_r, throttle)
      subroutine newtons_method(V_mag, height, euler, angular_rates, G)
        implicit none
        real , allocatable :: res(:), delta_G(:)
        real, intent(inout) :: G(6)
        real :: R_pos(6), R_neg(6), step_size
        real :: angular_rates(3), jac(6,6), euler(3), height, dummy_res(13), V_mag
        integer :: k, i, j

        allocate(res(6), delta_G(6))

        ! CALCULATE JACOBIAN AND RESIDUAL
        if (trim_verbose) then
          write(io_unit,'(A)') 'Building Jacobian Matrix:'
          write(io_unit,'(A)') ''
        end if 

        ! USE CENTRAL DIFFERENCE METHOD TO FIND JACOBIAN
        step_size = finite_difference_step_size
        do j = 1,6        
          G(j) = G(j) + step_size
          R_pos = calc_r(V_mag, height, euler, angular_rates, G)
          
          if (trim_verbose) then
            write(io_unit, '(A,I0,A)') 'Computing gradient relative to G[', j-1, ']'
            write(io_unit, '(A)') '   Positive Finite-Difference Step '
            write(io_unit, '(A,6(1X,ES20.12))') '      G =', (G(i), i=1,6)
            write(io_unit, '(A,6(1X,ES20.12))') '      R =', (R_pos(i), i=1,6)
          end if 
          
          G(j) = G(j) - 2 * step_size
          R_neg = calc_r(V_mag, height, euler, angular_rates, G)
          
          if (trim_verbose) then
            write(io_unit, '(A)') '   Negative Finite-Difference Step'
            write(io_unit, '(A,6(1X,ES20.12))') '      G =', (G(i), i=1,6)
            write(io_unit, '(A,6(1X,ES20.12))') '      R =', (R_neg(i), i=1,6)
            write(io_unit,'(A)') ''
          end if

          do i = 1,6
            jac(i, j) = (R_pos(i) - R_neg(i)) / (2*step_size)

          end do 
          G(j) = G(j) + step_size
        end do 

        if (trim_verbose) then
          write(io_unit, '(A)') 'Jacobian Matrix ='
          do i = 1, size(jac,1)
              write(io_unit,'(*(1X,ES20.12))') (jac(i,j), j=1,size(jac,2))
          end do
        end if
        res = calc_r(V_mag, height, euler, angular_rates, G)
        res = -1* res

        ! CALCUALTE DELTA G AND ADD RELAXATION FACTOR
        call lu_solve(6, jac, res, delta_G)
        G = G + relaxation_factor * delta_G

        if (trim_verbose) then
          write(io_unit,*)
          write(io_unit,'(A,6(1X,ES20.12))') 'Delta G =', (delta_G(k), k=1,6)
        end if 

      end subroutine newtons_method

    !=========================
    ! Calculate Residual
      function calc_r(V_mag, height, euler, angular_rates, G) result(R)
        implicit none
        real, intent(in) :: V_mag, height, G(6), angular_rates(3)
        real, intent(inout) :: euler(3)
        real :: alpha, ca, cb, sa, sb
        real :: R(6), dummy_R(13), temp_state(13)

        temp_state = 0.0

        ! PULL OUT CONTROLS
        controls(1:4) = G(3:6)

        ! PULL OUT ALPHA
        alpha = G(1)
        ca = cos(alpha)
        sa = sin(alpha)     

        ! PULL OUT BETA
        if (sideslip_angle0 /= -999.0) then 
          cb = cos(sideslip_angle0)
          sb = sin(sideslip_angle0)
          euler(1) = G(2)
        else 
          cb = cos(G(2))
          sb = sin(G(2))
        end if

        ! CALCUALTE INITIAL STATES
        temp_state(1:3)   = V_mag * (/ca*cb, sb, sa*cb/) 
        temp_state(4:6)   = angular_rates
        temp_state(9)     = height
        temp_state(7:8)   = 0
        temp_state(10:13) = euler_to_quat(euler)

        ! CALCULATE RESIDUAL
        dummy_R = differential_equations(0.0, temp_state)
        R = dummy_R(1:6)

      end function calc_r
  ! 
  ! MISCELLANEOUS FUNCTIONS
    !=========================
    ! Matrix Inverse
      function matrix_inv(A) result(A_inv)
          implicit none
          real, dimension(3,3), intent(in) :: A
          real, dimension(3,3) :: A_inv
          real :: det
          real :: a11, a12, a13, a21, a22, a23, a31, a32, a33

          ! UNPACK ELEMENTS
          a11 = A(1,1); a12 = A(1,2); a13 = A(1,3)
          a21 = A(2,1); a22 = A(2,2); a23 = A(2,3)
          a31 = A(3,1); a32 = A(3,2); a33 = A(3,3)

          ! CALCULATE DETERMINANT
          det = a11*(a22*a33 - a23*a32) - a12*(a21*a33 - a23*a31) + a13*(a21*a32 - a22*a31)

          ! CALCULATE INVERSE ELEMENT-WISE
          A_inv(1,1) =  (a22*a33 - a23*a32)/det
          A_inv(1,2) = -(a12*a33 - a13*a32)/det
          A_inv(1,3) =  (a12*a23 - a13*a22)/det

          A_inv(2,1) = -(a21*a33 - a23*a31)/det
          A_inv(2,2) =  (a11*a33 - a13*a31)/det
          A_inv(2,3) = -(a11*a23 - a13*a21)/det

          A_inv(3,1) =  (a21*a32 - a22*a31)/det
          A_inv(3,2) = -(a11*a32 - a12*a31)/det
          A_inv(3,3) =  (a11*a22 - a12*a21)/det

      end function matrix_inv


    !=========================
    ! Cross product
      function cross_product(a, b) result(c)
          implicit none
          real, intent(in) :: a(3), b(3)
          real :: c(3)

          ! DIRECTLY COMPUTE CROSS PRODUCT
          c(1) = a(2)*b(3) - a(3)*b(2)
          c(2) = a(3)*b(1) - a(1)*b(3)
          c(3) = a(1)*b(2) - a(2)*b(1)
      end function cross_product

    !=========================
    ! CONVERT DEG TO RAD IF NEEDED
      function to_radians_if_valid(angle) result(new_angle)
        implicit none
        real, intent(in) :: angle
        real :: new_angle
        if (angle == -999.0) then
          new_angle = angle
        else
          new_angle = angle * pi / 180.0
        end if
      end function to_radians_if_valid  
  ! 
  ! INITIALIZE AND RUN
    !=========================
    ! Init Subroutine
      subroutine init(filename)
        implicit none 
        real :: alpha_deg, beta_deg
        real, allocatable :: eul(:)
        real :: alpha, beta, trim_state(6)
        character(100), intent(in) :: filename
        type(json_value), pointer :: j_connections, j_graphics, j_user_controls, j_test_controls
        integer :: i

        ! OPEN A FILE TO WRITE TO 
        open(newunit=io_unit, file='f16_output.txt', status='replace', action='write')

        ! OPEN THE SPECIFIED JSON FILE
        call jsonx_load(filename, j_main) 

        ! DETERMINE IF RK4_VERBOSE
        call jsonx_get(j_main, 'simulation.rk4_verbose', rk4_verbose, .false.)
        call jsonx_get(j_main, 'simulation.print_states', print_states, .false.)

        ! DEFINE THRUST COEFFICIENTS FOR THRUST MODEL
        call jsonx_get(j_main, 'vehicle.thrust.T0[lbf]', T0)
        call jsonx_get(j_main, 'vehicle.thrust.Ta',      Ta)

        ! READ IN ALL AERODYNAMIC DATA
        call jsonx_get(j_main, 'vehicle.aerodynamics.compressibility',                   compressibility, .false.)
        call jsonx_get(j_main, 'vehicle.aerodynamics.test_compressibility',              test_compressibility, .false.)
        call jsonx_get(j_main, 'vehicle.aerodynamics.sweep[deg]',                        sweep, 0.0)
        call jsonx_get(j_main, 'vehicle.aerodynamics.reference.area[ft^2]',              planform_area)
        call jsonx_get(j_main, 'vehicle.aerodynamics.reference.longitudinal_length[ft]', longitudinal_length)
        call jsonx_get(j_main, 'vehicle.aerodynamics.reference.lateral_length[ft]',      lateral_length)

        call jsonx_get(j_main, 'vehicle.aerodynamics.stall.stall_flag',  stall)
        call jsonx_get(j_main, 'vehicle.aerodynamics.stall.test_stall', test_stall)
        call jsonx_get(j_main, 'vehicle.aerodynamics.stall.CL.lambda_b', CL_lambda_b)
        call jsonx_get(j_main, 'vehicle.aerodynamics.stall.CL.alpha_0',  CL_alpha_0)
        call jsonx_get(j_main, 'vehicle.aerodynamics.stall.CL.alpha_s',  CL_alpha_s)
        call jsonx_get(j_main, 'vehicle.aerodynamics.stall.CD.lambda_b', CD_lambda_b)
        call jsonx_get(j_main, 'vehicle.aerodynamics.stall.CD.alpha_0',  CD_alpha_0)
        call jsonx_get(j_main, 'vehicle.aerodynamics.stall.CD.alpha_s',  CD_alpha_s)
        call jsonx_get(j_main, 'vehicle.aerodynamics.stall.Cm.lambda_b', Cm_lambda_b)
        call jsonx_get(j_main, 'vehicle.aerodynamics.stall.Cm.alpha_0',  Cm_alpha_0)
        call jsonx_get(j_main, 'vehicle.aerodynamics.stall.Cm.alpha_s',  Cm_alpha_s)
        call jsonx_get(j_main, 'vehicle.aerodynamics.stall.Cm.min',      Cm_min)
        print_stall = .false.

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
        call jsonx_get(j_main, 'vehicle.aerodynamics.coefficients.CD.CL1_CL1',           CD_L1_L1)
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

        ! TRIM VARIABLES
        call jsonx_get(j_main, 'initial.type',                                    sim_type)
        call jsonx_get(j_main, 'initial.init_airspeed[ft/s]',                     init_airspeed)
        call jsonx_get(j_main, 'initial.altitude[ft]',                            initial_state(9))
        call jsonx_get(j_main, 'initial.Euler_angles[deg]',                       eul, 0.0, 3)

        call jsonx_get(j_main, 'initial.trim.exam_answers',                       exam_answers, .false.)      
        call jsonx_get(j_main, 'initial.trim.bank_angle[deg]',                    bank_angle0, -999.0)
        call jsonx_get(j_main, 'initial.trim.sideslip[deg]',                      sideslip_angle0, -999.0)
        call jsonx_get(j_main, 'initial.trim.elevation_angle[deg]',               elevation_angle0, -999.0)
        call jsonx_get(j_main, 'initial.trim.climb_angle[deg]',                   climb_angle0, -999.0)
        call jsonx_get(j_main, 'initial.trim.type',                               trim_type)
        call jsonx_get(j_main, 'initial.trim.verbose',                            trim_verbose, .false.)
        call jsonx_get(j_main, 'initial.trim.solver.relaxation_factor',           relaxation_factor)
        call jsonx_get(j_main, 'initial.trim.solver.tolerance',                   tolerance)
        call jsonx_get(j_main, 'initial.trim.solver.max_iterations',              max_iterations)
        call jsonx_get(j_main, 'initial.trim.solver.finite_difference_step_size', finite_difference_step_size)

        call jsonx_get(j_main, 'initial.state.alpha[deg]', alpha_deg)
        call jsonx_get(j_main, 'initial.state.beta[deg]', beta_deg)
        call jsonx_get(j_main, 'initial.state.p[deg/s]',  initial_state(4))
        call jsonx_get(j_main, 'initial.state.q[deg/s]',  initial_state(5))
        call jsonx_get(j_main, 'initial.state.r[deg/s]',  initial_state(6))

        ! CALCULATE MASS AND INERTIA
        call mass_inertia()         

        ! CONVERT ALTITUDE TO CORRECT DIRECTION
        initial_state(9) = initial_state(9) * (-1.0)

        ! CONVERT DEGREE VALUES TO RADIANS
        alpha = alpha_deg * pi / 180.0
        beta = beta_deg * pi / 180.0
        initial_state(4:6) = initial_state(4:6) * pi / 180.0
        eul = eul * pi / 180.0

        elevation_angle0 = to_radians_if_valid(elevation_angle0)
        climb_angle0     = to_radians_if_valid(climb_angle0)
        sideslip_angle0  = to_radians_if_valid(sideslip_angle0)
        bank_angle0      = to_radians_if_valid(bank_angle0)

        CL_alpha_0 = CL_alpha_0 * pi / 180.0
        CL_alpha_s = CL_alpha_s * pi / 180.0
        CD_alpha_0 = CD_alpha_0 * pi / 180.0
        CD_alpha_s = CD_alpha_s * pi / 180.0
        Cm_alpha_0 = Cm_alpha_0 * pi / 180.0
        Cm_alpha_s = Cm_alpha_s * pi / 180.0

        ! CALCULATE INITIAL SPEED IN U, V, W DIRECTIONS
        initial_state(1) = init_airspeed * cos(alpha) * cos(beta)
        initial_state(2) = init_airspeed * sin(beta)
        initial_state(3) = init_airspeed * sin(alpha) * cos(beta)

        ! CALCULATE THE INITIAL ORIENTATION
        initial_state(10:13) = euler_to_quat(eul)
        
        ! BUILD THE INITIAL CONTROL VECTOR
        call jsonx_get(j_main, 'initial.aileron[deg]',         controls(1))
        call jsonx_get(j_main, 'initial.elevator[deg]',        controls(2))
        call jsonx_get(j_main, 'initial.rudder[deg]',          controls(3))
        call jsonx_get(j_main, 'initial.throttle',             controls(4))

        controls(1:3) = controls(1:3) * pi / 180.0

        ! CONVERT SWEEP TO RADIANS
        sweep = sweep * pi / 180.0

        ! STORE THE DENSITY AT SEA LEVEL
        rho0 = 2.3768921839070335E-03

        ! CALCULATE THE SPECIFIED TRIM CONDITION IF GIVEN 
        if (sim_type == 'trim') then 
          trim_state = trim_algorithm(init_airspeed, initial_state(9), eul, tolerance, trim_type)
          call compute_A()
          call compute_B()
          ! write(io_unit,*) 'Amat:'
          ! do i = 1, 6
          !     write(io_unit,'(6F12.6)') Amat(i,1:6)
          ! end do
          ! write(io_unit,*) 'Bmat:'
          ! do i = 1, 6
          !     write(io_unit,'(4F12.6)') Bmat(i,1:4)
          ! end do     
        end if 

        call jsonx_get(j_main, 'initial.x_pos_ft', initial_state(7))
        
        ! INITIALIZE CONNECTIONS
        call jsonx_get(j_main, 'connections', j_connections)
        call jsonx_get(j_connections, 'graphics', j_graphics)
        write(*,*) 'call graphics%init(j_graphics)'
        call graphics%init(j_graphics)

        ! ADD CONTROLS
        call jsonx_get(j_connections, 'user_controls', j_user_controls)
        call user_controls%init(j_user_controls)

        ! TEST STALL IF SPECIFIED 
        if (test_stall) then 
          call check_stall(initial_state)
        end if 

        ! TEST COMPRESSIBILITY IF SPECIFIED
        if (test_compressibility) then 
          call check_compressibility(initial_state)
        end if 

      end subroutine init

    !=========================
    ! Run Subroutine
      subroutine run()
        implicit none
        real :: t, dt, tf, y_new(13), s(14) 
        real :: cpu_start_time, cpu_end_time, actual_time, integrated_time, time_error
        real :: time_1, time_2, y_init(13)
        real :: controls_input(6)
        logical :: real_time = .false.

        call jsonx_get(j_main, 'simulation.time_step[s]',  dt, 0.0)
        call jsonx_get(j_main, 'simulation.total_time[s]', tf)

        ! INITIALIZE TIME AND STATE
        t = 0.0 
        integrated_time = 0.0
        y_init = initial_state

        ! SWITCH TO REAL TIME SIMULATION IF SPECIFIED
        if(abs(dt) < tol) then
          real_time = .true.
          time_1 = get_time()
          y_init = rk4(t,initial_state,dt)
          call quat_norm(y_init)
          time_2 = get_time()
          dt = time_2 - time_1
          y_init = initial_state
        end if

        ! BUILD THE LOOP AND WRITE THE OUTPUT
        if (print_states) then 
          write(io_unit,*) " time[s]             u[ft/s]             &
          v[ft/s]             w[ft/s]             p[rad/s]            &
          q[rad/s]            r[rad/s]            xf[ft]              &
          yf[ft]              zf[ft]              e0                  &
          ex                  ey                  ez"
          write(io_unit,'(14ES20.12)') t,y_init(:)
        end if 

        ! SAVE THE TIMESTAMP WHEN THE SIMULATION BEGINS
        cpu_start_time = get_time()

        ! START THE SIMULATION
        do while(t < tf)
          ! CALCULATE THE NEW STATE
          y_new = rk4(t, y_init, dt)

          ! NORMALIZE THE QUATERNION
          call quat_norm(y_new(10:13))
          
          if (print_states) then 
            write(io_unit,'(14ES20.12)') t,y_init(:)
          end if 

          ! SEND GRAPHICS OVER CONNECTION 
          s(1) = t 
          s(2:14) = y_new(1:13)
          call graphics%send(s)

          ! RECEIVE USER CONTROLS OVER CONNECTION
          controls_input = user_controls%recv()
          controls_input(2:4) = controls_input(2:4) * pi / 180
          controls = controls_input(2:5)

          ! UPDATE THE STATE AND TIME        
          if(real_time) then
            time_2 = get_time()
            dt = time_2 - time_1
            time_1 = time_2
          end if 

          y_init = y_new
          t = t + dt
          integrated_time = integrated_time + dt
          write(*,*) t, dt
          
        end do 

        ! SAVE THE TIMESTAMP FOR WHEN THE SIMULATION STOPPED
        cpu_end_time = get_time()
        actual_time = cpu_end_time - cpu_start_time
        time_error = actual_time - integrated_time

        ! WRITE OUT THE TIME DIFFERENCES
        write(*,*) "Actual Time Elapsed", actual_time
        write(*,*) "Integrated Time", integrated_time
        write(*,*) "Time Error", time_error

        close(io_unit)

      end subroutine run
  ! 
  ! LQR CONTROLLER 
    !=========================
    ! ESTIMATE A MATRIX  
      subroutine compute_A()
        implicit none
        real :: A(13,6)
        real :: f_plus(13), f_minus(13)
        real :: x_plus(13), x_minus(13)
        real :: eps
        integer :: i

        ! FINITE DIFFERENCE STEP
        eps = 1.0e-6   
        controls = ud  

        ! ESTIMATE A WITH CENTRAL DIFFERENCE METHOD
        do i = 1, 6
          ! DESIRED STATE PERTURBATIONS
          x_plus  = xd
          x_minus = xd

          x_plus(i)  = x_plus(i)  + eps
          x_minus(i) = x_minus(i) - eps

          ! Evaluate the system derivatives
          f_plus = differential_equations(0.0, x_plus)
          f_minus =  differential_equations(0.0, x_minus)

          ! Central difference for column i of A
          A(:, i) = (f_plus - f_minus) / (2.0 * eps)
        end do

        Amat = A(1:6, :)

    end subroutine compute_A

    !=========================    
    ! ESTIMATE B MATRIX 
    subroutine compute_B()
        implicit none
        real :: f_plus(13), f_minus(13)
        real :: u_plus(4), u_minus(4)
        real :: eps
        integer :: i

        ! FINITE DIFFERENCE STEP
        eps = 1.0e-6
        ! ESTIMATE B WITH CENTRAL DIFFERENCE METHOD
        do i = 1, 4
            ! Perturb control i
            u_plus  = ud
            u_minus = ud

            u_plus(i)  = u_plus(i)  + eps
            u_minus(i) = u_minus(i) - eps

            ! Evaluate system derivatives with perturbed controls
            controls = u_plus
            f_plus  = differential_equations(0.0, xd)
            controls = u_minus
            f_minus = differential_equations(0.0, xd)

            ! Central difference for column i of B
            Bmat(:,i) = (f_plus(1:6) - f_minus(1:6)) / (2.0 * eps)
        end do

    end subroutine compute_B

    !=========================
    ! GET CONTROL INPUT
    function get_control(state) result(u)
      real :: Q(6,6), R(4,4), K(4,6), matlab_K(6,4)
      real :: u(4), state(13), state6(6), xd6(6)

      ! Balanced Q
      Q(1,1) = 0.25    ! u
      Q(2,2) = 0.25    ! v
      Q(3,3) = 0.25    ! w
      Q(4,4) = 100.0   ! p
      Q(5,5) = 100.0   ! q
      Q(6,6) = 100.0   ! r

      ! Balanced R (throttle, elevator, aileron, rudder)
      R(1,1) = 11.111111
      R(2,2) = 100.0
      R(3,3) = 100.0
      R(4,4) = 100.0

      matlab_K = reshape( (/ &
      -0.0000,   -0.0126,    0.0000,   -0.9212,    0.0001,    0.3333, &
       0.0040,    0.0000,   -0.0455,    0.0001,   -2.5648,   -0.0006, &
      -0.0000,    0.0399,    0.0000,    0.3865,   -0.0004,   -3.9043, &
       0.1493,   -0.0000,    0.0094,   -0.0000,   -0.3653,    0.0002  /), (/6,4/) )
      K = transpose(matlab_K)

      state6 = state(1:6)
      xd6 = xd(1:6)
      u = -matmul(K, state6 - xd6)

    end function 


end module f16_m