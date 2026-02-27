module controller_m
  use koch_m 
  use jsonx_m 
  use linalg_mod
  use connection_m 
  implicit none 

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

  type controller_t 
    type (pid_t) :: p_da, q_de, r_dr, bank_p, gamma_q, V_tau
    integer :: num_pid 
    type(connection) :: pilot_conn 
    logical :: running = .false. 
  end type controller_t 

  contains 
    
    subroutine controller_init(t, j_main) 
      implicit none 
      type(controller_t), intent(inout) :: t 
      type(json_value), pointer :: j_main, j_connections, j_pid, j_temp 

      write(*,*) 'Initializing controller...'
      t%running = .true. 

      write(*,*) 'Initializing connections...' 
      call jsonx_get(j_main, 'connections', j_connections) 
      call jsonx_get(j_connections, 'receive_pilot_commends', j_temp)
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

      pilot_command = t%pilot_conn%recv([time], time)
      bank_sp  = pilot_command(1) * pi / 180.0 
      gamma_sp = pilot_command(2) * pi /180.0 
      V_sp     = pilot_command(3) 

      p_sp = pid_get_command(t%bank_p, bank_sp, eul(1), time, dyp)
      q_sp = pid_get_command(t%gamma_q, gamma_sp, gamma, time, dyp)
      r_sp = (g*sp*ct + p*w)/u ! neglect gravity relief

      ans(1) = pid_get_command(t%p_da,  p_sp, p,    time, dyp)
      ans(2) = pid_get_command(t%q_de,  q_sp, q,    time, dyp)
      ans(3) = pid_get_command(t%r_dr,  r_sp, r,    time, dyp)
      ans(4) = pid_get_command(t%V_tau, V_sp, Vmag, time, dyp)

    end function controller_update

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
        if (dt>tol) t%error_deriv = (t%error  - t%prev_error) * dt 

        ! PID equation 
        ans = (t%KP*t%error + t%KI*t%error_int + t%KD*t%error_deriv) / dyp 
        if(t%dyp_schedule) ans = ans/dyp 

        t%prev_error = t%error 
        t%prev_time  = time 
        t%prev_ans   = ans 
      else ! don't update if not correct frequency 
        ans = t%prev_ans 
      end if 
    end function pid_get_command

  end module controller_m 