module controller_m
  use koch_m 
  use jsonx_m 
  use linalg_mod
  use connection_m 
  implicit none 
  
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
      function controller_update(t, vehicle, states, time) result(ans) 
        implicit none 
        type(controller_t), intent(inout) :: t 
        type(vehicle_t) :: vehicle 
        real, intent(in) :: states(24), time 
        real :: ans(4) 
        real :: u, v, w, p, q, r
        real :: eul(3), sp, st, cp, ct 
        real :: Vmag, gamma, g 
        real :: Z, temp, Pressure, rho, a, mu, dyp 
        real :: bank_sp, gamma_sp, V_sp, p_sp, q_sp, r_sp, setpoint(4)
        ! Command bank and climb angles
        real :: pilot_command(3)
        ! Command angular rates
        ! real :: pilot_command(4)

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

        pilot_command = t%pilot_conn%recv([time], time)
        
        ! Command bank angle, climb angle, velocity
        bank_sp  = pilot_command(1) * pi / 180.0 
        gamma_sp = pilot_command(2) * pi /180.0 
        V_sp     = pilot_command(3) 

        p_sp = pid_get_command(t%bank_p, bank_sp, eul(1), time, dyp)
        q_sp = pid_get_command(t%gamma_q, gamma_sp, gamma, time, dyp)
        r_sp = (g*sp*ct + p*w)/u ! neglect gravity relief

        ! Command a roll, pitch, yaw rate
        ! p_sp = pilot_command(1) * pi / 180.0 
        ! q_sp = pilot_command(2) * pi / 180.0 
        ! r_sp = pilot_command(3) * pi / 180.0 
        ! V_sp = pilot_command(4) 

        ! PID Controller
        ! ans(1) = pid_get_command(t%p_da,  p_sp, p,    time, dyp)
        ! ans(2) = pid_get_command(t%q_de,  q_sp, q,    time, dyp)
        ! ans(3) = pid_get_command(t%r_dr,  r_sp, r,    time, dyp)
        ! ans(4) = pid_get_command(t%V_tau, V_sp, Vmag, time, dyp)

        ! Dynamic Inversion Controller
        setpoint(1) = p_sp
        setpoint(2) = q_sp 
        setpoint(3) = r_sp 
        setpoint(4) = V_sp 
        ans = dynamic_inversion(vehicle, time, states, setpoint)

      end function controller_update
    !----------------------------------------
    ! Get commanded values for a controller
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
            ! write(*,*) 'PID controller saturated at ', t%prev_ans * t%display_units, '. Using integrator clamping...' 
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
      function dynamic_inversion(t, time, state, sp) result(u)
        implicit none 
        type(vehicle_t) :: t
        real, intent(in) :: time, state(24), sp(4)
        real :: u(4)
        real :: omega(3), omegadot_des(3), omega_ref(3), omega_ref_dot(3)
        real :: Kp(3,3), Ki(3,3)
        real :: Vmag, alpha, beta, dyp
        real :: altitude, gh_dum, T_dum, P_dum, rho, mu_dum, a_dum    
        real :: Ixx, Iyy, Izz, Ixy, Ixz, Iyz, Imat(3,3), Imat_inv(3,3)
        real :: hxb, hyb, hzb, hmat(3,3)
        real :: sign_a, sa, Cm_s, sigma_m, pos, neg 
        real :: angular_rates(3), pbar, qbar, rbar, p, q, r 
        real :: dyn_pressure_mat(3,3), C_states(3), C_control(3,3), pqr_term(3)
        real :: f(3), G(3,3), G_inv(3,3), delta(3), V_sp
        real :: error(3), int_error(3)

        omega = state(4:6)

        ! ! TRACK DESIRED EULER ANGLES
        ! ! Define desired angular acceleration with proportional controller
        ! error = omega - sp(1:3)
        ! int_error = state(22:24)
        ! omega_ref_dot = 0.0

        ! TRACK A REFERENCE SIGNAL
        ! Define reference signal
        omega_ref = pi/180.0 * [5.0*sin(2.0*time), -2.0*sin(time), sin(time)]
        omega_ref_dot = pi/180.0 * [10.0*cos(2.0*time), -2.0*cos(time), cos(time)] 
        error = omega - omega_ref
        int_error = state(22:24) ! integral error

        ! Define desired angular acceleration with proportional controller
        Kp = 0.0 
        Ki = 0.0
        Kp(1,1) = 4.0
        Kp(2,2) = 2.0
        Kp(3,3) = 4.0 
        Ki(1,1) = 1.0
        Ki(2,2) = 1.0
        Ki(3,3) = 1.0

        omegadot_des = -matmul(Kp, error) - matmul(Ki, int_error) + omega_ref_dot  

        ! Aircraft model start
          ! Velocity and angles
          Vmag = sqrt(state(1)**2 + state(2)**2 + state(3)**2)
          alpha  = atan2(state(3) , state(1))
          beta   = asin(state(2) / Vmag)

          ! Density
          altitude = -state(9)
          call std_atm_English(altitude, gh_dum, T_dum, P_dum, rho, mu_dum, a_dum )      

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

          pbar = 1 / (2*Vmag) * state(4) * t%lateral_length
          qbar = 1 / (2*Vmag) * state(5) * t%longitudinal_length
          rbar = 1 / (2*Vmag) * state(6) * t%lateral_length

          ! Dynamic Pressure
          dyn_pressure_mat = 0.0
          dyp = 0.5 * rho * Vmag**2
          dyn_pressure_mat(1,1) = dyp * t%planform_area * t%lateral_length
          dyn_pressure_mat(2,2) = dyp * t%planform_area * t%longitudinal_length
          dyn_pressure_mat(3,3) = dyp * t%planform_area * t%lateral_length

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
        ! Aircraft model end 

        ! f(x)
        f = matmul(Imat_inv, matmul(dyn_pressure_mat,C_states) + pqr_term + matmul(hmat,omega))

        ! G(x)
        G = matmul(Imat_inv, matmul(dyn_pressure_mat, C_control))      
        G_inv = matrix_inv(G)

        ! Dynamic Inversion
        delta = matmul(G_inv, (omegadot_des - f))
                
        ! Control Vector 
        u(1:3) = delta;
        u(4) = pid_get_command(t%controller%V_tau, sp(4), Vmag, time, dyp) ! PID for throttle

        ! Update integral error
        t%zdot = error

      end function dynamic_inversion
    !----------------------------------------
      
end module controller_m 