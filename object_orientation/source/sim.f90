module sim_m 
    use vehicle_m
    use jsonx_m 

    implicit none 
    integer :: io_unit
    type(json_value), pointer :: j_main
    type(vehicle_t), allocatable :: vehicles(:)
    integer :: num_vehicles

    contains 
    !----------------------------------------
    ! Init Subroutine
        subroutine init(filename)
            implicit none 
            character(100), intent(in) :: filename
            type(json_value), pointer :: j_vehicles, j_temp
            integer :: i
            logical :: save_states, rk4_verbose

            ! Begin simulation
            write(*,*) 'Initializing Simulation...'            
            call jsonx_load(filename, j_main)              

            ! Initialize vehicles
            write(*,*) 'Initializing vehicles...'
            call jsonx_get(j_main, 'vehicles', j_vehicles)
            num_vehicles = json_value_count(j_vehicles)
            allocate(vehicles(num_vehicles))

            ! Write statement flags
            call jsonx_get(j_main, 'simulation.rk4_verbose',   rk4_verbose, .false.)
            call jsonx_get(j_main, 'simulation.save_states',   save_states, .false.)
            call jsonx_get(j_main, 'simulation.save_lat_long', save_lat_long, .false.)

            ! Geographic model
            call jsonx_get(j_main, 'simulation.geographic_model', geographic_model, 'none')
            geographic_model_ID = 0 ! flat earth default
            if (geographic_model == 'sphere') geographic_model_ID = 1
            if (geographic_model == 'ellipse') geographic_model_ID = 2

            ! Initialize vehicles
            do i=1,num_vehicles
                vehicles(i)%save_states = save_states
                vehicles(i)%rk4_verbose = rk4_verbose
            end do    

            do i = 1,num_vehicles
                call json_value_get(j_vehicles, i, j_temp)
                call vehicle_init(vehicles(i), j_temp)
            end do 
        end subroutine init

    !----------------------------------------
    ! Run Subroutine
        subroutine run()
            implicit none
            real :: time, dt, tf
            real :: cpu_start_time, cpu_end_time, actual_time, integrated_time, time_error
            real :: time_1 = 0.0, time_2 = 0.0 
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

            ! SAVE THE TIMESTAMP WHEN THE SIMULATION BEGINS
            cpu_start_time = get_time()

            ! START THE SIMULATION
            do i=1,num_vehicles
                if(vehicles(i)%run_physics) then 
                    if(vehicles(i)%save_states) call vehicle_write_state(vehicles(i), time)
                    if(save_lat_long) call vehicle_write_lat_long(vehicles(i), time)
                end if 
            end do     
                    
            do while (time < tf - tol)
                ! CALCULATE THE NEW STATES FOR EACH VEHICLE
                do i=1,num_vehicles
                    if(vehicles(i)%run_physics) then 
                        call vehicle_tick_state(vehicles(i), time, dt)     
                    end if 
                end do 

                ! UPDATE THE TIME        
                if(real_time) then
                    time_2 = get_time()
                    dt = time_2 - time_1
                    time_1 = time_2
                end if 

                time = time + dt
                integrated_time = integrated_time + dt
                write(*,*) time, dt
                
                do i=1,num_vehicles
                    if(vehicles(i)%run_physics) then 
                        if(vehicles(i)%save_states) call vehicle_write_state(vehicles(i), time)
                        if(save_lat_long) call vehicle_write_lat_long(vehicles(i), time)
                    end if 
                end do                 
                
                ! if (abs(time - 100.0*anint(time/100.0)) < 1e-2) then
                !     do i = 1, num_vehicles
                !         if (vehicles(i)%run_physics) then 
                !             if (vehicles(i)%save_states) call vehicle_write_state(vehicles(i), time)
                !             if (save_lat_long) call vehicle_write_lat_long(vehicles(i), time)
                !         end if
                !     end do
                ! end if
                
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

    !----------------------------------------


end module sim_m 