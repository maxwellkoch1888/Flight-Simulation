module sim_m 
    use vehicle_m
    use jsonx_m 

    implicit none 
    integer :: io_unit
    type(json_value), pointer :: j_main
    type(vehicle_t), allocatable :: vehicles(:)
    integer :: num_vehicles

    contains 
    !=========================
    ! Init Subroutine
        subroutine init(filename)
            implicit none 
            character(100), intent(in) :: filename
            type(json_value), pointer :: j_vehicles, j_temp
            integer :: i

            ! Open file to write to  
            open(newunit=io_unit, file='f16_output.txt', status='replace', action='write')

            ! Begin simulation
            write(*,*) 'Initializing Simulation...'            
            call jsonx_load(filename, j_main) 

            ! ! Debugging flags
            ! call jsonx_get(j_main, 'simulation.rk4_verbose', rk4_verbose, .false.)
            ! call jsonx_get(j_main, 'simulation.save_states', save_states, .false.)

            ! Initialize vehicles
            write(*,*) 'Initializing vehicles...'
            call jsonx_get(j_main, 'vehicles', j_vehicles)
            num_vehicles = json_value_count(j_vehicles)
            allocate(vehicles(num_vehicles))

            do i = 1,num_vehicles
                call json_value_get(j_vehicles, i, j_temp)
                call vehicle_init(vehicles(i), j_temp)
            end do 
        end subroutine init

    !=========================
    ! Run Subroutine
        subroutine run()
            implicit none
            real :: time, dt, tf
            real :: cpu_start_time, cpu_end_time, actual_time, integrated_time, time_error
            real :: time_1, time_2
            integer :: i
            logical :: real_time = .false.

            call jsonx_get(j_main, 'simulation.time_step[s]',  dt, 0.0)
            call jsonx_get(j_main, 'simulation.end_time[s]', tf)

            ! INITIALIZE TIME AND STATE
            time = 0.0 
            integrated_time = 0.0

            ! SWITCH TO REAL TIME SIMULATION IF SPECIFIED
            if(abs(dt) < tol) then
                real_time = .true.
                time_1 = get_time()
                do i=1,num_vehicles
                    if (vehicles(i)%run_physics) call vehicle_tick_state(vehicles(i), time, dt) 
                end do 

                time_2 = get_time()
                dt = time_2 - time_1
                do i=1,num_vehicles
                    if (vehicles(i)%run_physics) vehicles(i)%state = vehicles(i)%init_state 
                end do 
            end if

            ! BUILD THE LOOP AND WRITE THE OUTPUT
            ! if (save_states) then 
            !     write(io_unit,*) " time[s]             u[ft/s]             &
            !     v[ft/s]             w[ft/s]             p[rad/s]            &
            !     q[rad/s]            r[rad/s]            xf[ft]              &
            !     yf[ft]              zf[ft]              e0                  &
            !     ex                  ey                  ez"
            !     write(io_unit,'(14ES20.12)') time,y_init(:)
            ! end if 

            ! SAVE THE TIMESTAMP WHEN THE SIMULATION BEGINS
            cpu_start_time = get_time()
            ! START THE SIMULATION
            do while(time < tf)
                ! CALCULATE THE NEW STATEs FOR EACH VEHICLE
                do i=1,num_vehicles
                    if(vehicles(i)%run_physics) call vehicle_tick_state(vehicles(i), time, dt)
                end do 

                ! UPDATE THE STATE AND TIME        
                if(real_time) then
                    time_2 = get_time()
                    dt = time_2 - time_1
                    time_1 = time_2
                end if 

                time = time + dt
                integrated_time = integrated_time + dt
                write(*,*) time, dt
                
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

end module sim_m 