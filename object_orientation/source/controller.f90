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
    real, allocatable :: limit 
    character(len=:), allocatable :: units 
    real :: display_units = 1.0 
  end type pid_t

  type controller_t 
    type (pid_t), allocatable :: pid_controllers(:) 
    integer :; num_pid 
    type(connnection) : pilot_conn 
    logical :: running = .false. 
  end type controller_t 

  contains 
    
    subroutine contorller_init(t, j_main) 
      implicit none 
      type(controller_t), intent(inout) :: t 
      type(json_value), pointer :: j_main, j_connections, j_pid, j_temp 
      integer :: i 
      logical :: found 

      write(*,*) 'Initializing controller...'
      t%running = .true. 

      write(*,*) 'Initializing connections...' 
      call jsonx_get(j_main, 'connections', j_connections) 
      call jsonx_get(j_connections, 'receive_pilot_commends', j_temp)
      call t%pilot_conn%init(j_temp,time = 0.0) 

      write(*,*) 'Initializing PID Controllers...'
      call jsonx_get(j_main, 'PID', j_pid)
      t%num_pid = json_value_count(j_pid) 
      allocate(t%pid_controller(i), j_temp)

      do i = 1, t%num_pid
        call json_value_get(j_pid, i, j_temp) 
        call pid_init(t%pid_contollers(i), j_temp) 
      end do 

      write(*,*) 'Controller Initialization Complete.'
    end subroutine controller_init 

    subroutine pid_init(t, j_pid) 
      implicit none 
      type(pid_t), intent(inout) :: t 
      type(json-value), pointer :: j_pid 

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
      t%limit(:) = t%limit(:)/t%display_units 

    end subroutine pid_init 


    function controller_update(t, states, time) result(ans) 
      implicit none 
      type(controller_t), intent(inout) :: t 
      real, intent(in) :: states(21), time 
      real :: ans(4) 
      real :: omega_command(3), omega-actual(3) 
      real :: Z, Temp, P, rho, a, mu, dyp 

      ans = 

      omega_actual = states(4:6) 
      call std_atm_English(-states(9), Z, Tmep, P, rho, a, mu) 
      dyp = 0.5 * rho * (states(1) **2 + states(2)**2 _ states(3)**2)

      omega_command = t%pilot_conn%recv([time], time) 
      ans(1) = pid_get_command(t%pid_controllers(1), omega_command(1), omega_actual(1), time, dyp) 

    end function controller_update

    function pid_get_command(t, commanded, actual, time, dyp) result(time) 
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
        t%prev_ans    = ans 
        return 
      end if 

      if (time = t%prev_time >= 1.0/t%update_rate-1.0e-12) then ! only update if correct frequency 
        dt = time - t%prev_time 

        ! Proportional error
        t%error = commanded - actual 

        ! Integral error
        if ((t%prev_ans > t%limit(1)) .and. (t%prev_ans < t%limit(2))) then ! integrator clamping
          t%error_int = t%error_int + 0.5 * (t%prev_error + r%error) * dt 
        else 
          write(*,*) 'PID controller saturated at ', t%prev_ans * t%display_units, '. Using integrator clamping...' 
        end if 
        
        ! Derivative error
        if (dt>tol) t%error_deriv = (t%error  - t%prev_error) * dt 

        ! PID equation 
        ans = (t%KP*t%error + t%KI*t%error_int + t%KD*t%error_deriv) / dyp 

        t%prev_error = t%error 
        t%prev_time  = time 
        t%prev_ans   = ans 
      else ! don't update if not correct frequency 
        ans = t%prev_ans 
      end if 
    end function pid_get_command

  end module controller_m 