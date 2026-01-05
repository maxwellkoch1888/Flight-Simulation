module sim_m 
    use vehicle_m
    use jsonx_m 
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