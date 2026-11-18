module sim_m 
    use vehicle_m
    use jsonx_m 

    implicit none 
    integer :: io_unit
    type(json_value), pointer :: j_main
    type(vehicle_t), allocatable :: vehicles(:)
    type(atmosphere_t) :: atm
    type(connection) :: graphics, user_controls
    integer :: num_vehicles
    logical :: input_controls = .false. 

    contains 
    !----------------------------------------
    ! Init Subroutine
        subroutine init(filename)
            implicit none 
            character(100), intent(in) :: filename
            type(json_value), pointer :: j_vehicles, j_temp, j_atmosphere, j_connections, j_graphics, j_user_controls
            integer :: i
            logical :: save_states, rk4_verbose, sigsev, found 

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
            call jsonx_get(j_main, 'simulator.sigsev_error', sigsev, .false.)

            ! Geographic model
            call jsonx_get(j_main, 'simulation.geographic_model', geographic_model, 'none')
            geographic_model_ID = 0 ! flat earth default
            if (geographic_model == 'sphere') geographic_model_ID = 1
            if (geographic_model == 'ellipse') geographic_model_ID = 2

            ! Initialize atmosphere
            write(*,*) 'Reading Atmospheric JSON Object...'
            call jsonx_get(j_main, 'atmosphere', j_atmosphere)

            ! Initialize connections 
            write(*,*) 'Reading connections json object...'
            call jsonx_get(j_main, 'connections', j_connections) 

            call jsonx_get(j_connections, 'graphics', j_graphics)
            call graphics%init(j_graphics)

            call json_get(j_connections, 'user_controls', j_user_controls, found)
            if (found) then 
                input_controls = .true. 
                call user_controls%init(j_user_controls)
            end if 


            ! Initialize vehicles
            do i=1,num_vehicles
                vehicles(i)%save_states = save_states
                vehicles(i)%rk4_verbose = rk4_verbose
                vehicles(i)%sigsev = sigsev
            end do    

            do i = 1,num_vehicles
                call json_value_get(j_vehicles, i, j_temp)
                call atmosphere_init(vehicles(i)%atm, j_atmosphere)
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
            real :: unreal_engine_state(14)
            real :: controls_input(6), commands(3)
            integer :: i
            logical :: real_time = .false. 

            call jsonx_get(j_main, 'simulation.time_step[s]',  dt, 0.0)
            call jsonx_get(j_main, 'simulation.end_time[s]', tf)

            ! Initialize time and state
            time = 0.0 
            integrated_time = 0.0

            ! Switch to real time sim if specified
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

            ! Save timestamp when sim begins
            cpu_start_time = get_time()

            ! Set velocity and controller gains if user controls turned on 
            if(input_controls) then 
                vehicles(1)%controller%gamma_q%KP = 0.4
                vehicles(1)%controller%gamma_q%KI = 0.05
                vehicles(1)%controller%gamma_q%KD = 0.0
                vehicles(1)%controller%gamma_q%update_rate = 1000
                vehicles(1)%controller%gamma_q%limit = [-1000.0, 1000.0]
                vehicles(1)%controller%gamma_q%display_units = 180.0/ pi 

                vehicles(1)%controller%V_tau%KP = 0.01 
                vehicles(1)%controller%V_tau%KI = 0.002
                vehicles(1)%controller%V_tau%KD = 0.0     
                vehicles(1)%controller%V_tau%update_rate = 1000
                vehicles(1)%controller%V_tau%limit = [0.0, 1.0]
                vehicles(1)%controller%V_tau%display_units = 1.0                
                commands(3) = sqrt(vehicles(1)%state(1)**2 + vehicles(1)%state(2)**2 + vehicles(1)%state(3)**2)
            end if 

            ! Start sim
            do i=1,num_vehicles
                if(vehicles(i)%run_physics) then 
                    if(vehicles(i)%save_states) call vehicle_write_state(vehicles(i), time)
                    if(save_lat_long) call vehicle_write_lat_long(vehicles(i), time)
                end if 
            end do     
            
            write(*,*) 'Simulation Running...'
            if(abs(dt) < tol) write(*,*) 'Real time Simulation, check Unreal Engine...'
            do while (time < tf - tol)
                ! New state for each vehicle
                do i=1,num_vehicles
                    if(vehicles(i)%run_physics) then 
                        call vehicle_tick_state(vehicles(i), time, dt)     
                        if(vehicles(i)%name == 'F16') unreal_engine_state(2:14) = vehicles(i)%state(1:13)
                    end if 
                end do 

                ! Send graphics over connection
                unreal_engine_state(1) = time
                call graphics%send(unreal_engine_state)

                ! Receive user controls 
                if(input_controls) then 
                    controls_input = user_controls%recv()

                    commands(1) = -60.0 * controls_input(2)
                    commands(2) = -60.0 * controls_input(3)
                    commands(3) = commands(3) + controls_input(5) * 0.001
                    vehicles(1)%state(14:17) = user_controller(vehicles(1)%controller, commands, vehicles(1), vehicles(1)%state, time)
                end if           

                ! Update time      
                if(real_time) then
                    time_2 = get_time()
                    dt = time_2 - time_1
                    time_1 = time_2
                end if 

                time = time + dt
                integrated_time = integrated_time + dt
                ! write(*,*) time, dt
                
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

            ! Save timestamp
            cpu_end_time = get_time()
            actual_time = cpu_end_time - cpu_start_time
            time_error = actual_time - integrated_time

            ! Write time differences
            write(*,*) "Actual Time Elapsed", actual_time
            write(*,*) "Integrated Time", integrated_time
            write(*,*) "Time Error", time_error

            do i=1,num_vehicles
                close(vehicles(i)%iunit_states)
            end do 

            close(io_unit)

        end subroutine run

    !----------------------------------------


end module sim_m 