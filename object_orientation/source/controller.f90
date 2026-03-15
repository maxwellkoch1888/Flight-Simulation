module controller_m
  use koch_m 
  use jsonx_m 
  use linalg_mod
  use connection_m 
  implicit none 
  
  !==================================================
  ! TYPE DECLARATIONS
  !==================================================
    !----------------------------------------
    ! PID type 
      type pid_t 
        character(len=:), allocatable :: name 
        real :: KP, KI, KD 
        real :: error, prev_error, error_int, error_deriv 
        real :: prev_time, prev_ans, update_rate 
        real, allocatable :: limit(:)
        character(len=:), allocatable :: units 
        real :: display_units = 1.0 
        logical :: dyp_schedule
      end type pid_t
    !----------------------------------------
    ! Controller type
      type controller_t 
        type (pid_t) :: p_da, q_de, r_dr, bank_p, gamma_q, V_tau
        integer :: num_pid 
        type(connection) :: pilot_conn 
        logical :: running = .false. 
      end type controller_t 
    !----------------------------------------
    ! Trim solver type
      type trim_solver_t
        real :: step_size, relaxation_factor, tolerance, max_iterations
      end type trim_solver_t
    !----------------------------------------
    ! Trim settings type
      type trim_settings_t
        character(len=:), allocatable :: type 

        logical, allocatable :: free_vars(:) 
        real :: climb_angle, sideslip_angle
        real :: load_factor
        logical :: verbose, solve_relative_climb_angle, solve_load_factor

        type(trim_solver_t) :: solver 
      end type trim_settings_t
    !----------------------------------------
    ! Control type
      type control_t 
        integer :: dynamics_order, state_ID 
        real :: commanded_value
        real, allocatable :: mag_limit(:), rate_limit(:), accel_limit(:) 
        real :: time_constant, natural_frequency, damping_ratio
      end type control_t       
      
    !----------------------------------------        
    ! Vehicle type
      
      type vehicle_t
        type(json_value), pointer :: j_vehicle
            
        character(len=:), allocatable :: name
        character(len=:), allocatable :: type
        character(100) :: states_filename, rk4_filename, trim_filename, latlong_filename

        logical :: run_physics
        logical :: save_states, limit_controls = .true. 
        integer :: iunit_states, iunit_rk4, iunit_trim, iunit_latlong

        ! Location 
        real :: latitude, longitude, prev_latitude, prev_longitude
        real :: course_angle

        ! Mass
        real :: mass
        real :: inertia(3,3)
        real :: inertia_inv(3,3)
        real, allocatable :: h(:)

        ! Aerodynamics
        real, allocatable :: aero_ref_location(:)
        real :: planform_area, longitudinal_length, lateral_length, sweep
        real :: CL0, CL_alpha, CL_alphahat, CL_qbar, CL_elevator
        real :: CS_beta, CS_pbar, CS_alpha_pbar, CS_rbar, CS_aileron, CS_rudder
        real :: CD_L0, CD_L1, CD_L1_L1, CD_CS_CS, CD_qbar, CD_alpha_qbar, CD_elevator, CD_alpha_elevator, CD_elevator_elevator
        real :: Cl_l0, Cl_beta, Cl_pbar, Cl_rbar, Cl_alpha_rbar, Cl_aileron, Cl_rudder
        real :: Cm_0, Cm_alpha, Cm_qbar, Cm_alphahat, Cm_elevator
        real :: Cn_beta, Cn_pbar, Cn_alpha_pbar, Cn_rbar, Cn_aileron, Cn_alpha_aileron, Cn_rudder
        real :: Cm_alpha_0, Cm_alpha_s, Cm_min

        ! Thrust 
        real :: T0, Ta, thrust_quat(4), rho0
        real, allocatable :: thrust_location(:), thrust_orientation(:)

        ! Debugging
        logical :: compressibility = .false., rk4_verbose, print_states, test_compressibility

        ! Stall model 
        logical :: stall = .false., test_stall
        real :: CL_lambda_b, CL_alpha_0, CL_alpha_s, CD_lambda_b, CD_alpha_0, CD_alpha_s, Cm_lambda_b
        
        ! Initialization constants
        real :: init_airspeed, init_alt, init_state(21)
        real, allocatable :: init_eul(:) ! has to be allocatable because will be read from json object

        ! States/controls
        real :: state(21)
        type(control_t) :: controls(4)
        type(controller_t) :: controller 

        type(trim_settings_t) :: trim
      end type vehicle_t
    !----------------------------------------   
  contains 
  !==================================================
  ! INITIALIZATION FUNCTIONS
  !==================================================    
    !----------------------------------------
    ! Controller initialization  
      subroutine controller_init(t, j_main) 
        implicit none 
        type(controller_t), intent(inout) :: t 
        type(json_value), pointer :: j_main, j_connections, j_pid, j_temp 

        write(*,*) 'Initializing controller...'
        t%running = .true. 

        write(*,*) 'Initializing connections...' 
        call jsonx_get(j_main, 'connections', j_connections) 
        call jsonx_get(j_connections, 'receive_pilot_commands', j_temp)
        call t%pilot_conn%init(j_temp,time = 0.0) 

        write(*,*) 'Initializing PID Controllers...'
        call jsonx_get(j_main, 'PID', j_pid)

        call json_value_get(j_pid, 'p->aileron', j_temp) 
        call pid_init(t%p_da, j_temp) 

        call json_value_get(j_pid, 'q->elevator', j_temp) 
        call pid_init(t%q_de, j_temp) 
        
        call json_value_get(j_pid, 'r->rudder', j_temp) 
        call pid_init(t%r_dr, j_temp)

        call json_value_get(j_pid, 'bank->p', j_temp) 
        call pid_init(t%bank_p, j_temp)

        call json_value_get(j_pid, 'gamma->q', j_temp) 
        call pid_init(t%gamma_q, j_temp)

        call json_value_get(j_pid, 'V->throttle', j_temp) 
        call pid_init(t%V_tau, j_temp)

        write(*,*) 'Controller Initialization Complete.'
      end subroutine controller_init 
    !----------------------------------------
    ! PID Initialization
      subroutine pid_init(t, j_pid) 
        implicit none 
        type(pid_t), intent(inout) :: t 
        type(json_value), pointer :: j_pid 

        t%name = j_pid%name 
        write(*,*) 'Initializing PID for ', t%name 
        call jsonx_get(j_pid, 'update_rate[hz]', t%update_rate) 
        call jsonx_get(j_pid, 'KP', t%KP)
        call jsonx_get(j_pid, 'KI', t%KI)
        call jsonx_get(j_pid, 'KD', t%KD)

        t%prev_error =  0.0 
        t%error_int  =  0.0 
        t%prev_ans   =  0.0 
        t%prev_time  = -1.0 

        call jsonx_get (j_pid, 'units', t%units, 'none') 
        if(t%units == 'deg') then 
          t%display_units = 180.0/pi 
        else 
          t%display_units = 1.0 
        end if 
        call jsonx_get(j_pid, 'limits', t%limit, 0.0, 2) 
        call jsonx_get(j_pid, 'dyp_schedule', t%dyp_schedule, .false.)
        t%limit(:) = t%limit(:)/t%display_units 

      end subroutine pid_init 
  !==================================================
  ! UPDATE FUNCTIONS
  !==================================================         
    !----------------------------------------
    ! Update a controller
      function controller_update(t, states, time) result(ans) 
        implicit none 
        type(controller_t), intent(inout) :: t 
        real, intent(in) :: states(21), time 
        real :: ans(4) 
        real :: u, v, w, p, q, r
        real :: eul(3), sp, st, cp, ct 
        real :: Vmag, gamma, g 
        real :: Z, temp, Pressure, rho, a, mu, dyp 
        real :: pilot_command(3), bank_sp, gamma_sp, V_sp, p_sp, q_sp, r_sp

        u = states(1)
        v = states(2)
        w = states(3)
        p = states(4)
        q = states(5)
        r = states(6)

        eul = quat_to_euler(states(10:13))
        sp  = sin(eul(1))
        st  = sin(eul(2))
        cp  = cos(eul(1))
        ct  = cos(eul(2))

        Vmag = sqrt(u**2 + v**2 + w**2) 

        gamma = asin((u*st - (v*sp + w*cp) * ct) / Vmag) 
        call std_atm_English(-states(9), Z, temp, Pressure, rho, a, mu) 
        dyp = 0.5 * rho * Vmag**2
        g = gravity_English(-states(9))

        ! Command orientation, climb, velocity
        pilot_command = t%pilot_conn%recv([time], time)
        bank_sp  = pilot_command(1) * pi / 180.0 
        gamma_sp = pilot_command(2) * pi /180.0 
        V_sp     = pilot_command(3) 

        p_sp = pid_get_command(t%bank_p, bank_sp, eul(1), time, dyp)
        q_sp = pid_get_command(t%gamma_q, gamma_sp, gamma, time, dyp)
        r_sp = (g*sp*ct + p*w)/u ! neglect gravity relief

        !12.6.2 example
        ! p_sp = 0.0
        ! q_sp = 0.0
        ! r_sp = 0.0

        ans(1) = pid_get_command(t%p_da,  p_sp, p,    time, dyp)
        ans(2) = pid_get_command(t%q_de,  q_sp, q,    time, dyp)
        ans(3) = pid_get_command(t%r_dr,  r_sp, r,    time, dyp)
        ans(4) = pid_get_command(t%V_tau, V_sp, Vmag, time, dyp)

        ! 12.6.2 example
        ! ans(4) = 0.76052985298660364E-01

      end function controller_update
    !----------------------------------------
    ! Get commanded values for a PID controller
      function pid_get_command(t, commanded, actual, time, dyp) result(ans) 
        implicit none 
        type(pid_t), intent(inout) :: t 
        real, intent(in) :: commanded, actual, time, dyp 
        real :: dt, ans 

        if(t%prev_time < 0.0) then ! first iteration update
          t%prev_time   = time 
          t%prev_error  = commanded - actual 
          t%error_int   = 0.0 
          t%error_deriv = 0.0 
          ans           = t%KP * t%prev_error 
          if(t%dyp_schedule) ans = ans/dyp 
          t%prev_ans    = ans 
          return 
        end if 

        if (time - t%prev_time >= 1.0/t%update_rate-1.0e-12) then ! only update if correct frequency 
          dt = time - t%prev_time 

          ! Proportional error
          t%error = commanded - actual 

          ! Integral error
          if ((t%prev_ans > t%limit(1)) .and. (t%prev_ans < t%limit(2))) then ! integrator clamping
            t%error_int = t%error_int + 0.5 * (t%prev_error + t%error) * dt 
          else 
            write(*,*) 'PID controller saturated at ', t%prev_ans * t%display_units, '. Using integrator clamping...' 
          end if 
          
          ! Derivative error
          if (dt>tol) t%error_deriv = (t%error  - t%prev_error) / dt 

          ! PID equation 
          ans = (t%KP*t%error + t%KI*t%error_int + t%KD*t%error_deriv) 
          if(t%dyp_schedule) ans = ans/dyp 

          t%prev_error = t%error 
          t%prev_time  = time 
          t%prev_ans   = ans 
        else ! don't update if not correct frequency 
          ans = t%prev_ans 
        end if 
      end function pid_get_command
    !----------------------------------------
    ! Dynamic Inversion Controller
      function dynamic_inversion(t, time, state) result(u)
        implicit none 
        type(vehicle_t) :: t
        type(pid_t) :: tau_controller
        real, intent(in) :: time, state(21)
        real :: u(4)
        real :: pilot_command(3)
        real :: omega(3), omegadot_des(3), Kp(3,3)
        real :: Vmag, alpha, beta, dyp, altitude
        real :: geometric_altitude_ft, geopotential_altitude_ft
        real :: temp_R, pressure_lbf_per_ft2, density_slugs_per_ft3
        real :: dyn_viscosity_slug_per_ft_sec, sos_ft_per_sec      
        real :: Ixx, Iyy, Izz, Ixy, Ixz, Iyz, Imat(3,3), Imat_inv(3,3)
        real :: hxb, hyb, hzb, hmat(3,3)
        real :: sign_a, sa, Cm_s, sigma_m, pos, neg 
        real :: angular_rates(3), pbar, qbar, rbar, p, q, r 
        real :: dyn_pressure_mat(3,3), C_states(3), C_control(3,3), pqr_term(3)
        real :: f(3), G(3,3), G_inv(3,3), delta(3), V_sp

        omega = state(4:6)

        ! Define desired angular acceleration with proportional controller
        Kp = 0.0
        Kp(1,1) = 2.0
        Kp(2,2) = 2.0
        Kp(3,3) = 2.0
        omegadot_des = matmul(-Kp, omega); ! desired is zero     

        ! ! Define reference signal
        ! omega_ref = pi/180 * [5*sin(2*t); -2*sin(t); sin(t)]; 
        ! omega_ref_dot = pi/180 * [10*cos(2*t); -2*cos(t); cos(t)]; 
        ! error = omega - omega_ref;

        ! ! Define error
        ! int_error = state(13:15); % integral error

        ! ! Define desired angular acceleration with proportional controller
        ! Kp = diag([4,2,4]);
        ! Ki = diag([1,1,1]); 
        ! omegadot_des = -Kp*error - Ki*int_error + omega_ref_dot;         

        ! Aircraft model
          ! Velocity and angles
          Vmag = sqrt(state(1)**2 + state(2)**2 + state(3)**2)
          alpha  = atan2(state(3) , state(1))
          beta   = asin(state(2) / Vmag)

          ! Density
            geometric_altitude_ft = -state(9)
            call std_atm_English(&
              geometric_altitude_ft, geopotential_altitude_ft,     & 
              temp_R, pressure_lbf_per_ft2, density_slugs_per_ft3, & 
              dyn_viscosity_slug_per_ft_sec, sos_ft_per_sec)      

          ! Inertia
          Ixx = t%inertia(1,1)
          Iyy = t%inertia(2,2)
          Izz = t%inertia(3,3)
          Ixy = t%inertia(1,2)
          Ixz = t%inertia(1,3)
          Iyz = t%inertia(2,3)    
          
          Imat(1,1) =  Ixx
          Imat(1,2) = -Ixy
          Imat(1,3) = -Ixz

          Imat(2,1) = -Ixy
          Imat(2,2) =  Iyy
          Imat(2,3) = -Iyz

          Imat(3,1) = -Ixz
          Imat(3,2) = -Iyz
          Imat(3,3) =  Izz

          Imat_inv = matrix_inv(Imat)   
          
          ! Angular momentum
          hxb = t%h(1)
          hyb = t%h(2)
          hzb = t%h(3)

          hmat(1,1) =  0.0
          hmat(1,2) = -hzb
          hmat(1,3) =  hyb

          hmat(2,1) =  hzb
          hmat(2,2) =  0.0
          hmat(2,3) = -hxb

          hmat(3,1) = -hyb
          hmat(3,2) =  hxb
          hmat(3,3) =  0.0

          ! Angular rates
          p = state(4)
          q = state(5)
          r = state(6)
          angular_rates(1) = 1 / (2*Vmag) * state(4) * t%lateral_length
          angular_rates(2) = 1 / (2*Vmag) * state(5) * t%longitudinal_length
          angular_rates(3) = 1 / (2*Vmag) * state(6) * t%lateral_length

          pbar = angular_rates(1)
          qbar = angular_rates(2)
          rbar = angular_rates(3)

          ! Dynamic Pressure
          dyn_pressure_mat = 0.0
          dyn_pressure_mat(1,1) = 0.5 * density_slugs_per_ft3 * Vmag **2 * t%planform_area * t%lateral_length
          dyn_pressure_mat(2,2) = 0.5 * density_slugs_per_ft3 * Vmag **2 * t%planform_area * t%longitudinal_length
          dyn_pressure_mat(3,3) = 0.5 * density_slugs_per_ft3 * Vmag **2 * t%planform_area * t%lateral_length

          ! pqr_term
          pqr_term(1) = (Iyy - Izz)*q*r + Iyz*(q**2 -r**2) + Ixz*p*q - Ixy*p*r
          pqr_term(2) = (Izz - Ixx)*p*r + Ixz*(r**2 -p**2) + Ixy*q*r - Iyz*p*q
          pqr_term(3) = (Ixx - Iyy)*p*q + Ixy*(p**2 -q**2) + Iyz*p*r - Ixz*q*r        

          ! Initial Moment Coefficient
          C_states(1) = t%Cl_beta * beta + t%Cl_pbar * pbar + (t%Cl_rbar + t%Cl_alpha_rbar * alpha) * rbar
          C_states(2) = t%Cm_0 + t%Cm_alpha * alpha + t%Cm_qbar * qbar
          C_states(3) = t%Cn_beta * beta + (t%Cn_pbar + t%Cn_alpha_pbar * alpha) * pbar + t%Cn_rbar * rbar

          C_control = 0.0 
          C_control(1,1) = t%Cl_aileron
          C_control(1,3) = t%Cl_rudder
          C_control(2,2) = t%Cm_elevator
          C_control(3,1) = t%Cn_aileron + t%Cn_alpha_aileron * alpha
          C_control(3,3) = t%Cn_rudder

          ! Compressibility model
          ! if (t%compressibility) then 
          !   CM1      = 2.13/ (t%sweep + 0.15)**2
          !   CM2      = 15.35*t%sweep**2 - 19.64*t%sweep +16.86
          !   mach_num = Vmag / sos_ft_per_sec

          !   ! Mach breakpoints
          !   m_low         = 0.60    ! only use prandtl-glauert
          !   m_high        = 0.92    ! only use karman-tsien
          !   m_trans_start = 0.88    ! start transsonic region

          !   ! Prandtl-Glauert factor
          !   sqrt_term = sqrt(max(1.0 - mach_num**2, tol))
          !   pg_factor = 1.0 / sqrt_term

          !   ! Karman-Tsien factor
          !   kt_factor = pg_factor * (1.0 + mach_num**2 / (1.0 + sqrt_term))

          !   ! Blending 
          !   if (mach_num <= m_low) then
          !     blend = 0.0
          !   else if (mach_num >= m_high) then
          !     blend = 1.0
          !   else
          !     blend = (mach_num - m_low) / (m_high - m_low)
          !   end if

          !   ! Final compressibility factor
          !   moment_factor = lift_factor

          !   ! Apply factors
          !   max_Cl_roll_factor = 2.5

          !   if (mach_num < 0.92) then ! accurate range
          !     C_states  = C_states  * moment_factor
          !     C_control = C_control * moment_factor
          !   else
          !     C_states  = C_states  * min(moment_factor, max_Cl_roll_factor)
          !     C_control = C_control * min(moment_factor, max_Cl_roll_factor)          
          !   end if

          ! end if 

          ! Stall model 
          if (t%stall) then 
            sign_a  = sign(1.0,alpha)
            sa = sin(alpha)

            Cm_s    = t%Cm_min * sa**2 * sign_a

            pos = exp( t%Cm_lambda_b * (alpha - t%Cm_alpha_0 + t%Cm_alpha_s))
            neg = exp(-t%Cm_lambda_b * (alpha - t%Cm_alpha_0 - t%Cm_alpha_s))
            sigma_m = (1.0 + neg + pos) / ((1.0 + neg)*(1.0 + pos))       

            C_control(2, 2) = C_control(2, 2) * (1 - sigma_m) + Cm_s * sigma_m 
            C_states(2)     = C_states(2)     * (1 - sigma_m) + Cm_s * sigma_m 
          end if 

        ! f(x)
        f = matmul(Imat_inv, matmul(dyn_pressure_mat,C_states) + pqr_term + matmul(hmat,omega))
        ! G(x)
        G = matmul(Imat_inv, matmul(dyn_pressure_mat, C_control))      
        G_inv = matrix_inv(G)

        ! Dynamic Inversion
        delta = matmul(G_inv, (omegadot_des - f))
        
        ! PID Terms
        dyp = 0.5 * density_slugs_per_ft3 * Vmag**2
        V_sp = pilot_command(3) 
        
        ! Control Vector 
        u(1:3) = delta;
        u(4) = pid_get_command(tau_controller, V_sp, Vmag, time, dyp) ! PID for throttle
      end function dynamic_inversion
    !----------------------------------------
      
end module controller_m 