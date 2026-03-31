module vehicle_m
  use koch_m
  use jsonx_m
  use micro_time_m
  use linalg_mod
  use controller_m 
  use database_m
  use propulsion_m 

  implicit none
  real :: rho0
  real, parameter :: earth_radius_ft = 6366707.01949371/0.3048
  integer :: geographic_model_ID
  character(len=:), allocatable :: geographic_model 
  logical :: save_lat_long

  contains
  !==================================================
  ! INITIALIZATION FUNCTIONS
  !==================================================
    !----------------------------------------
    ! Vehicle initialization
      subroutine vehicle_init(t, j_vehicle_input)
        implicit none 
        type(vehicle_t), intent(inout) :: t
        type(json_value), pointer :: j_vehicle_input, j_controller, j_control, j_control_temp, j_propulsion, j_temp 
        character(len=:), allocatable :: init_type 
        real :: geopotential_altitude_ft,temp_R, pressure_lbf_per_ft2, density_slugs_per_ft3, dyn_viscosity_slug_per_ft_sec, sos_ft_per_sec
        integer :: i 
        logical :: found, allow_saturation

        t%j_vehicle => j_vehicle_input 
        t%name = t%j_vehicle%name 

        write(*,*) 'Initializing ', t%name 
        call jsonx_get(t%j_vehicle, 'type',t%type)
        write(*,*) '- type = ', t%type 
        call jsonx_get(t%j_vehicle, 'run_physics', t%run_physics, .true.)
        write(*,*) '- run physics = ', t%run_physics

        if(t%run_physics) then 
          if(t%save_states) then 
            t%states_filename = 'output_files/' // trim(t%name) // '_states.csv'
            open(newunit=t%iunit_states, file=t%states_filename, status='REPLACE')
            write(t%iunit_states,'(*(A,:,","))') &
                'time[s]', 'u[ft/s]', 'v[ft/s]', 'w[ft/s]', &
                'p[rad/s]', 'q[rad/s]', 'r[rad/s]', &
                'xf[ft]', 'yf[ft]', 'zf[ft]', &
                'e0', 'ex', 'ey', 'ez', &
                'aileron[deg]', 'elevator[deg]', 'rudder[deg]', 'throttle', &
                'daileron[deg]', 'delevator[deg]', 'drudder[deg]', 'dthrottle', &
                'p_integral_error', 'q_integral_error', 'r_integral_error'
            write(*,*) '- saving states to ', t%states_filename
          end if 

          if(t%rk4_verbose) then 
            t%rk4_filename = 'output_files/' // trim(t%name)//'_RK4.txt'
            open(newunit=t%iunit_rk4, file=t%rk4_filename, status='REPLACE')
            write(*,*) '- saving RK4 results to ', t%rk4_filename
          end if 

          if(save_lat_long) then 
            t%latlong_filename = 'output_files/' // trim(t%name)//'_latlong.txt'
            open(newunit=t%iunit_latlong, file=t%latlong_filename, status='REPLACE')
            write(t%iunit_latlong,*) "time[s]              longitude[deg]       &
            &latitude[deg]        azimuth[deg]"              
            write(*,*) '- saving latlong results to ', t%latlong_filename
          end if 

          write(*,*) '- mass' 
          call jsonx_get(t%j_vehicle, 'mass.weight[lbf]', t%mass)
          t%mass = t%mass/gravity_English(0.0)

          ! Read mass properties
          call jsonx_get(t%j_vehicle, 'mass.Ixx[slug-ft^2]',  t%inertia(1,1))
          call jsonx_get(t%j_vehicle, 'mass.Iyy[slug-ft^2]',  t%inertia(2,2))
          call jsonx_get(t%j_vehicle, 'mass.Izz[slug-ft^2]',  t%inertia(3,3))
          call jsonx_get(t%j_vehicle, 'mass.Ixy[slug-ft^2]',  t%inertia(1,2), 0.0)
          call jsonx_get(t%j_vehicle, 'mass.Ixz[slug-ft^2]',  t%inertia(1,3), 0.0)
          call jsonx_get(t%j_vehicle, 'mass.Iyz[slug-ft^2]',  t%inertia(2,3), 0.0)
          call jsonx_get(t%j_vehicle, 'mass.h[slug-ft^2/s]',  t%h, 0.0, 3)

          ! Read aerodynamic data      
          write(*,*) '- aerodynamics'
          call jsonx_get(t%j_vehicle, 'aerodynamics.compressibility',                   t%compressibility, .false.)
          call jsonx_get(t%j_vehicle, 'aerodynamics.test_compressibility',              t%test_compressibility, .false.)
          call jsonx_get(t%j_vehicle, 'aerodynamics.sweep[deg]',                        t%sweep, 0.0)
          call jsonx_get(t%j_vehicle, 'aerodynamics.reference.area[ft^2]',              t%planform_area, 0.0)
          call jsonx_get(t%j_vehicle, 'aerodynamics.reference.longitudinal_length[ft]', t%longitudinal_length, 0.0)
          call jsonx_get(t%j_vehicle, 'aerodynamics.reference.lateral_length[ft]',      t%lateral_length, 0.0 )
          call jsonx_get(t%j_vehicle, 'aerodynamics.reference.location[ft]',            t%aero_ref_location, 0.0)

          if(t%type == 'arrow') then 
            call jsonx_get(t%j_vehicle, 'aerodynamics.coefficients.CL.alpha',    t%CL_alpha)
            call jsonx_get(t%j_vehicle, 'aerodynamics.coefficients.CD.L0',       t%CD_L0)
            call jsonx_get(t%j_vehicle, 'aerodynamics.coefficients.CD.CL1_CL1',  t%CD_L1_L1)
            call jsonx_get(t%j_vehicle, 'aerodynamics.coefficients.Cl.0',        t%Cl_l0)
            call jsonx_get(t%j_vehicle, 'aerodynamics.coefficients.Cl.pbar',     t%Cl_pbar)
            call jsonx_get(t%j_vehicle, 'aerodynamics.coefficients.Cm.alpha',    t%Cm_alpha)
            call jsonx_get(t%j_vehicle, 'aerodynamics.coefficients.Cm.qbar',     t%Cm_qbar)
          end if 

          if(t%type == 'aircraft') then 
            call jsonx_get(t%j_vehicle, 'aerodynamics.coefficients.CL.0',        t%CL0)
            call jsonx_get(t%j_vehicle, 'aerodynamics.coefficients.CL.alpha',    t%CL_alpha)
            call jsonx_get(t%j_vehicle, 'aerodynamics.coefficients.CL.alphahat', t%CL_alphahat)
            call jsonx_get(t%j_vehicle, 'aerodynamics.coefficients.CL.qbar',     t%CL_qbar)
            call jsonx_get(t%j_vehicle, 'aerodynamics.coefficients.CL.elevator', t%CL_elevator)

            call jsonx_get(t%j_vehicle, 'aerodynamics.coefficients.CS.beta',       t%CS_beta)
            call jsonx_get(t%j_vehicle, 'aerodynamics.coefficients.CS.pbar',       t%CS_pbar)
            call jsonx_get(t%j_vehicle, 'aerodynamics.coefficients.CS.alpha_pbar', t%CS_alpha_pbar)
            call jsonx_get(t%j_vehicle, 'aerodynamics.coefficients.CS.rbar',       t%CS_rbar)
            call jsonx_get(t%j_vehicle, 'aerodynamics.coefficients.CS.aileron',    t%CS_aileron)
            call jsonx_get(t%j_vehicle, 'aerodynamics.coefficients.CS.rudder',     t%CS_rudder)

            call jsonx_get(t%j_vehicle, 'aerodynamics.coefficients.CD.L0',                t%CD_L0)
            call jsonx_get(t%j_vehicle, 'aerodynamics.coefficients.CD.CL1',               t%CD_L1)
            call jsonx_get(t%j_vehicle, 'aerodynamics.coefficients.CD.CL1_CL1',           t%CD_L1_L1)
            call jsonx_get(t%j_vehicle, 'aerodynamics.coefficients.CD.CS_CS',             t%CD_CS_CS)
            call jsonx_get(t%j_vehicle, 'aerodynamics.coefficients.CD.qbar',              t%CD_qbar)
            call jsonx_get(t%j_vehicle, 'aerodynamics.coefficients.CD.alpha_qbar',        t%CD_alpha_qbar)
            call jsonx_get(t%j_vehicle, 'aerodynamics.coefficients.CD.elevator',          t%CD_elevator)
            call jsonx_get(t%j_vehicle, 'aerodynamics.coefficients.CD.alpha_elevator',    t%CD_alpha_elevator)
            call jsonx_get(t%j_vehicle, 'aerodynamics.coefficients.CD.elevator_elevator', t%CD_elevator_elevator)

            call jsonx_get(t%j_vehicle, 'aerodynamics.coefficients.Cl.beta',       t%Cl_beta)
            call jsonx_get(t%j_vehicle, 'aerodynamics.coefficients.Cl.pbar',       t%Cl_pbar)
            call jsonx_get(t%j_vehicle, 'aerodynamics.coefficients.Cl.rbar',       t%Cl_rbar)
            call jsonx_get(t%j_vehicle, 'aerodynamics.coefficients.Cl.alpha_rbar', t%Cl_alpha_rbar)
            call jsonx_get(t%j_vehicle, 'aerodynamics.coefficients.Cl.aileron',    t%Cl_aileron)
            call jsonx_get(t%j_vehicle, 'aerodynamics.coefficients.Cl.rudder',     t%Cl_rudder)

            call jsonx_get(t%j_vehicle, 'aerodynamics.coefficients.Cm.0',        t%Cm_0)
            call jsonx_get(t%j_vehicle, 'aerodynamics.coefficients.Cm.alpha',    t%Cm_alpha)
            call jsonx_get(t%j_vehicle, 'aerodynamics.coefficients.Cm.qbar',     t%Cm_qbar)
            call jsonx_get(t%j_vehicle, 'aerodynamics.coefficients.Cm.alphahat', t%Cm_alphahat)
            call jsonx_get(t%j_vehicle, 'aerodynamics.coefficients.Cm.elevator', t%Cm_elevator)

            call jsonx_get(t%j_vehicle, 'aerodynamics.coefficients.Cn.beta',          t%Cn_beta)
            call jsonx_get(t%j_vehicle, 'aerodynamics.coefficients.Cn.pbar',          t%Cn_pbar)
            call jsonx_get(t%j_vehicle, 'aerodynamics.coefficients.Cn.alpha_pbar',    t%Cn_alpha_pbar)
            call jsonx_get(t%j_vehicle, 'aerodynamics.coefficients.Cn.rbar',          t%Cn_rbar)
            call jsonx_get(t%j_vehicle, 'aerodynamics.coefficients.Cn.aileron',       t%Cn_aileron)
            call jsonx_get(t%j_vehicle, 'aerodynamics.coefficients.Cn.alpha_aileron', t%Cn_alpha_aileron)
            call jsonx_get(t%j_vehicle, 'aerodynamics.coefficients.Cn.rudder',        t%Cn_rudder)

            ! Aerodynamic Databases
            call jsonx_get(t%j_vehicle, 'aerodynamics.rectilinear_databases', t%db_fn, 'none') 
            if(trim(t%db_fn(1)) /= 'none') then 
              t%n_db = size(t%db_fn) 
              allocate(t%db(t%n_db))
              call jsonx_get(t%j_vehicle, 'aerodynamics.use_database',                    t%use_database)
              call jsonx_get(t%j_vehicle, 'aerodynamics.database_directory',              t%db_path, '')
              call jsonx_get(t%j_vehicle, 'aerodynamics.database_allow_past_saturation',  allow_saturation)
              call jsonx_get(t%j_vehicle, 'aerodynamics.speed_brake[deg]',                t%speed_brake)
              call jsonx_get(t%j_vehicle, 'aerodynamics.leading_edge_flap[deg]',          t%le_flap) 
              t%speed_brake = t%speed_brake * pi / 180.0 
              t%le_flap     = t%le_flap * pi / 180.0 
              do i = 1,t%n_db 
                call t%db(i)%init(t%db_fn(i), pn=t%db_path, verbose=.true., presorted=.true.)
                t%db(i)%saturate = allow_saturation 
              end do 
            end if 

            ! Controller
            call json_get(t%j_vehicle, 'controller', j_controller, found)
            if (found) then 
              write(*,*) '   -controller'
              call controller_init(t%controller, j_controller)
            end if 

            ! Initialize integral error 
            t%zdot = 0.0 

          end if 

          ! Stall conditions
          if(t%type == 'arrow' .or. t%type == 'aircraft') then 
            call jsonx_get(t%j_vehicle, 'aerodynamics.stall.include_stall', t%stall)
            if(t%stall) then 
              call jsonx_get(t%j_vehicle, 'aerodynamics.stall.CL.lambda_b',      t%CL_lambda_b)
              call jsonx_get(t%j_vehicle, 'aerodynamics.stall.CL.alpha_0[deg]',  t%CL_alpha_0)
              call jsonx_get(t%j_vehicle, 'aerodynamics.stall.CL.alpha_s[deg]',  t%CL_alpha_s)
              call jsonx_get(t%j_vehicle, 'aerodynamics.stall.CD.lambda_b',      t%CD_lambda_b)
              call jsonx_get(t%j_vehicle, 'aerodynamics.stall.CD.alpha_0[deg]',  t%CD_alpha_0)
              call jsonx_get(t%j_vehicle, 'aerodynamics.stall.CD.alpha_s[deg]',  t%CD_alpha_s)
              call jsonx_get(t%j_vehicle, 'aerodynamics.stall.Cm.lambda_b',      t%Cm_lambda_b)
              call jsonx_get(t%j_vehicle, 'aerodynamics.stall.Cm.alpha_0[deg]',  t%Cm_alpha_0)
              call jsonx_get(t%j_vehicle, 'aerodynamics.stall.Cm.alpha_s[deg]',  t%Cm_alpha_s)
              call jsonx_get(t%j_vehicle, 'aerodynamics.stall.Cm.min',           t%Cm_min)

              t%CL_alpha_0 = t%CL_alpha_0 * pi / 180.0
              t%CL_alpha_s = t%CL_alpha_s * pi / 180.0
              t%CD_alpha_0 = t%CD_alpha_0 * pi / 180.0
              t%CD_alpha_s = t%CD_alpha_s * pi / 180.0
              t%Cm_alpha_0 = t%Cm_alpha_0 * pi / 180.0
              t%Cm_alpha_s = t%Cm_alpha_s * pi / 180.0
            end if 
          end if 

          ! Thrust coefficients
          write(*,*) '- propulsion'
          call json_get(t%j_vehicle, 'propulsion', j_propulsion) 
          t%num_props = json_value_count(j_propulsion) 
          allocate(t%props(t%num_props)) 

          do i=1,t%num_props 
            call json_value_get(j_propulsion, i, j_temp) 
            call propulsion_init(t%props(i), j_temp) 
          end do 

          ! Control Effectors
          write(*,*) '  -control effectors'
          call jsonx_get(t%j_vehicle, 'control_effectors', j_control) 
          call jsonx_get(j_control, '1', j_control_temp)
          call init_control(t, j_control_temp, 1) 
          call jsonx_get(j_control, '2', j_control_temp)
          call init_control(t, j_control_temp, 2) 
          call jsonx_get(j_control, '3', j_control_temp)
          call init_control(t, j_control_temp, 3) 
          call jsonx_get(j_control, '4', j_control_temp)
          call init_control(t, j_control_temp, 4)           

          ! Initial conditions
          write(*,*) '- Initial Conditions'
          t%init_state = 0.0 

          call jsonx_get(t%j_vehicle, 'initial.airspeed[ft/s]', t%init_airspeed)
          call jsonx_get(t%j_vehicle, 'initial.altitude[ft]',   t%init_alt)
          t%init_state(9) = -t%init_alt 
          call jsonx_get(t%j_vehicle, 'initial.latitude[deg]',  t%latitude)
          call jsonx_get(t%j_vehicle, 'initial.longitude[deg]', t%longitude)
          t%latitude  = t%latitude  * pi / 180.0
          t%longitude = t%longitude * pi / 180.0

          call jsonx_get(t%j_vehicle, 'initial.Euler_angles[deg]', t%init_eul, 0.0, 3)
          t%init_eul = t%init_eul * pi / 180.0 

          call jsonx_get(t%j_vehicle, 'initial.type', init_type)

          if(init_type == 'state') then 
            write(*,*) 'initializing state...'
            call init_to_state(t)
          else 
            write(*,*) 'initializing trim state...'
            call init_to_trim(t) 
          end if 
          write(*,*) 'init euler:'
          write(*,*) t%init_eul * 180.0 / pi 
          t%init_state(10:13) = euler_to_quat(t%init_eul) 
          t%state = t%init_state 
          write(*,*) 'initializing course angle...'
          t%course_angle   = t%init_eul(3)
          t%prev_latitude  = t%latitude - cos(t%course_angle)
          t%prev_longitude = t% longitude - sin(t%course_angle) 
          t%atm%prev_xyz = t%init_state(7:9)      

          ! call get_controller_input(t, 0.0)

        end if 
        write(*,*) 'Finished vehicle initialization.'
        
      end subroutine vehicle_init

    !----------------------------------------
    ! Control effector initialization
      subroutine init_control(t, j_control, ID) 
        implicit none 
        type(vehicle_t) :: t 
        type(json_value), pointer :: j_control 
        integer :: ID, i 

        call jsonx_get(j_control, 'name', t%controls(ID)%name) 
        write(*,*) '      -reading control effector: ', t%controls(ID)%name 
        
        if(t%controls(ID)%name == 'aileron') t%aileron_ID = ID + 13
        if(t%controls(ID)%name == 'elevator') t%elevator_ID = ID + 13
        if(t%controls(ID)%name == 'rudder') t%rudder_ID = ID + 13

        call jsonx_get(j_control, 'dynamics_order', t%controls(ID)%dynamics_order, 0) 
        call jsonx_get(j_control, 'units', t%controls(ID)%units, 'none') 
        if(t%controls(ID)%units == 'deg') then 
          t%controls(ID)%display_units = 180.0/pi 
        else 
          t%controls(ID)%display_units = 1.0 
        end if 
        call jsonx_get(j_control, 'magnitude_limits', t%controls(ID)%mag_limit, 0.0, 2) 
        call jsonx_get(j_control, 'rate_limits[/s]', t%controls(ID)%rate_limit, 0.0, 2) 
        t%controls(ID)%mag_limit(:) = t%controls(ID)%mag_limit(:) / t%controls(ID)%display_units 

        ! First order dynamics 
        if(t%controls(ID)%dynamics_order == 1) then 
          call jsonx_get(j_control, 'rate_limits[/s]', t%controls(ID)%rate_limit, 0.0, 2) 
          t%controls(ID)%rate_limit(:) = t%controls(ID)%rate_limit(:) / t%controls(ID)%display_units 
          call jsonx_get(j_control, 'time_constant[s]', t%controls(ID)%time_constant) 
        end if 

        ! Second order dynamics 
        if(t%controls(ID)%dynamics_order == 2) then 
          call jsonx_get(j_control, 'rate_limits[/s]', t%controls(ID)%rate_limit, 0.0, 2) 
          t%controls(ID)%rate_limit(:) = t%controls(ID)%rate_limit(:) / t%controls(ID)%display_units 

          call jsonx_get(j_control, 'acceleration_limits[/s^2]', t%controls(ID)%accel_limit, 0.0, 2) 
          t%controls(ID)%accel_limit(:) = t%controls(ID)%accel_limit(:) / t%controls(ID)%display_units 

          call jsonx_get(j_control, 'natural_frequency[rad/s]', t%controls(ID)%natural_frequency)
          call jsonx_get(j_control, 'damping_ratio', t%controls(ID)%damping_ratio) 
        end if 

        do i=1,t%num_props 
          if(t%controls(ID)%name == t%props(i)%name) then 
            t%props(i)%control_ID = ID 
            t%props(i)%units = t%controls(ID)%units 
          end if 
        end do         

        t%controls(ID)%state_ID = 13 + ID 

        t%controls(ID)%commanded_value = 0.0 
      end subroutine init_control 
    !----------------------------------------
    ! Write states to a file
      subroutine vehicle_write_state(t, time)
        implicit none 
        type(vehicle_t) :: t 
        real, intent(in) :: time 

        write(t%iunit_states,'(*(g0,","))') &
            time, t%state(1), t%state(2), t%state(3), t%state(4), t%state(5), t%state(6), &
            t%state(7), t%state(8), t%state(9), t%state(10), t%state(11), t%state(12), t%state(13), &
            t%state(14)*180.0/pi, t%state(15)*180.0/pi, t%state(16)*180.0/pi, &
            t%state(17), t%state(18)*180.0/pi, t%state(19)*180.0/pi, &
            t%state(20)*180/pi, t%state(21), t%state(22)*180.0/pi, t%state(23)*180.0/pi, t%state(24)*180.0/pi

      end subroutine 
    !----------------------------------------
    ! Write lat/long to a file
      subroutine vehicle_write_lat_long(t, time)
        implicit none 
        type(vehicle_t) :: t 
        real, intent(in) :: time 
        real :: eul(3)

        eul = quat_to_euler(t%state(10:13))
        write(t%iunit_latlong,'(14(ES20.13,1X))') time, t%latitude * 180.0 / pi , t%longitude * 180.0 / pi, eul(3) * 180.0 / pi 

      end subroutine         
    !----------------------------------------
    ! State initial condition
      subroutine init_to_state(t)
        implicit none 
        type(vehicle_t) :: t
        type(json_value), pointer :: j_initial, j_state 
        real :: alpha, beta
        real, allocatable :: temp_controls(:)
        integer :: i  

        write(*,*) '    - setting prescribed state'
        ! Get json state object 
        call jsonx_get(t%j_vehicle, 'initial', j_initial)
        call jsonx_get(j_initial, 'state', j_state)

        call jsonx_get(j_state, 'angle_of_attack[deg]', alpha) 
        call jsonx_get(j_state, 'sideslip_angle[deg]',  beta) 
        alpha = alpha * pi / 180.0
        beta  = beta  * pi / 180.0

        ! Calculate initial velocities
        t%init_state(1) = t%init_airspeed * cos(alpha) * cos(beta)
        t%init_state(2) = t%init_airspeed * sin(beta)
        t%init_state(3) = t%init_airspeed * sin(alpha) * cos(beta)

        ! Read initial angular velocities 
        call jsonx_get(j_state, 'p[deg/s]', t%init_state(4))
        call jsonx_get(j_state, 'q[deg/s]', t%init_state(5))
        call jsonx_get(j_state, 'r[deg/s]', t%init_state(6))
        t%init_state(4:6) = t%init_state(4:6) * pi / 180.0 

        ! Read initial controls
        if(t%type == 'aircraft' .or. t%type == 'quadrotor') then 
          call jsonx_get(j_state, 'controls', temp_controls, 0.0, 4) 
          do i=1,4 
            t%init_state(i+13) = temp_controls(i)/t%controls(i)%display_units 
            t%controls(i)%commanded_value = t%init_state(i+13) 
          end do 
        end if     
      end subroutine 
    !----------------------------------------
    ! Trim initial condition
      subroutine init_to_trim(t) 
        implicit none 
        type(vehicle_t) :: t
        type(json_value), pointer :: j_initial, j_trim, j_solver 
        integer, parameter :: n_vars = 9
        integer :: n_free, iter, i, j, k 
        real :: x(n_vars) 
        real, allocatable :: res(:), R_pos(:), R_neg(:), jac(:,:), delta_x(:) 
        integer, allocatable, dimension(:) :: idx_free 
        real :: error
        logical :: temp_turb 

        t%limit_controls = .false. 
        temp_turb = t%atm%use_turb
        t%atm%use_turb = .false. 

        allocate(t%trim%free_vars(n_vars))

        write(*,*) '  -trimming'
        ! Get json objects 
        call jsonx_get(t%j_vehicle, 'initial', j_initial) 
        call jsonx_get(j_initial, 'trim',      j_trim)
        call jsonx_get(j_trim, 'solver',       j_solver)
        call jsonx_get(j_trim, 'type',         t%trim%type)

        call jsonx_get(j_trim, 'sideslip_angle[deg]',  t%trim%sideslip_angle, -999.0) 
        call jsonx_get(j_trim, 'ref_climb_angle[deg]', t%trim%climb_angle, -999.0) 
        call jsonx_get(j_trim, 'load_factor',          t%trim%load_factor, -999.0) 

        write(*,*) '   -trimming vehicle for ', t%trim%type 

        call jsonx_get(j_solver, 'finite_difference_step_size', t%trim%solver%step_size)
        call jsonx_get(j_solver, 'relaxation_factor',           t%trim%solver%relaxation_factor)
        call jsonx_get(j_solver, 'tolerance',                   t%trim%solver%tolerance)
        call jsonx_get(j_solver, 'max_iterations',              t%trim%solver%max_iterations)
        call jsonx_get(j_solver, 'verbose',                     t%trim%verbose)

        if(t%trim%verbose) then 
          t%trim_filename = 'output_files/' // trim(t%name)//'_trim.txt'
          open(newunit=t%iunit_trim, file=t%trim_filename, status='REPLACE')
          write(t%iunit_trim, *) 'Trimming for type: ', t%trim%type
          write(t%iunit_trim, *) 
          write(t%iunit_trim, *) 'Newton Solver Settings: '
          write(t%iunit_trim, *) 'Finite Difference Step Size = ', t%trim%solver%step_size
          write(t%iunit_trim, *) '          Relaxation Factor = ', t%trim%solver%relaxation_factor 
          write(t%iunit_trim, *) '                  Tolerance = ', t%trim%solver%tolerance
        end if 

        write(*,*) 
        write(*,*) 'Newton Solver Settings: '
        write(*,*) 'Finite Difference Step Size = ', t%trim%solver%step_size
        write(*,*) '          Relaxation Factor = ', t%trim%solver%relaxation_factor 
        write(*,*) '                  Tolerance = ', t%trim%solver%tolerance

        t%trim%solve_relative_climb_angle = .false. 
        t%trim%solve_load_factor =          .false. 

        x = 0.0 
        x(7:9) = t%init_eul(:) 
        t%trim%free_vars(:) =   .true. 
        t%trim%free_vars(7:9) = .false. 

        if(t%trim%type == 'shss' .and. abs(t%trim%sideslip_angle+999.0) > tol) then 
          x(2) = t%trim%sideslip_angle
          x(2) = x(2) * pi / 180.0
          t%trim%free_vars(2) = .false. 
          t%trim%free_vars(7) = .true.
        end if 

        if(abs(t%trim%climb_angle+999.0) > tol) then 
          t%trim%climb_angle = t%trim%climb_angle * pi / 180.0
          x(8) = t%trim%climb_angle
          write(*,*) ' Reference Climb Angle[deg] = ', t%trim%climb_angle * 180.0 / pi 
          write(t%iunit_trim,*) ' Reference Climb Angle[deg] = ', t%trim%climb_angle * 180.0 / pi 
          t%trim%free_vars(8) = .true. 
          t%trim%solve_relative_climb_angle = .true. 
        end if 

        if(t%trim%type == 'sct' .and. abs(t%trim%load_factor+999.0) > tol) then 
          t%trim%solve_load_factor = .true. 
          x(7) = acos(cos(x(8))/t%trim%load_factor) ! load factor approximation
          t%trim%free_vars(7) = .true. 
          write(*,*) '                Load Factor =', t%trim%load_factor
          write(t%iunit_trim,*) '                Load Factor =', t%trim%load_factor
        end if 

        n_free = count(t%trim%free_vars)

        allocate(res(n_free))
        allocate(R_pos(n_free))
        allocate(R_neg(n_free))
        allocate(jac(n_free, n_free))
        
        idx_free = pack([(i,i=1,n_vars)], t%trim%free_vars) ! yields an array of indices of all free variables

        if(t%trim%verbose) then 
          write(t%iunit_trim, *)
          write(t%iunit_trim, *) ' x array defined as [alpha, beta, aileron, elevator, rudder, throttle, phi, theta, psi]'
          write(t%iunit_trim, *) 'n_vars = ', n_vars 
          write(t%iunit_trim, *) 'n_free = ', n_free 
          write(t%iunit_trim, *) 'free_vars = ', t%trim%free_vars
          write(t%iunit_trim, *) 'idx_free = ', idx_free(:)
        end if 
        write(*,*) '    n_vars = ', n_vars 
        write(*,*) '    n_free = ', n_free 
        write(*,*) '    free_vars = ', t%trim%free_vars
        write(*,*) '    idx_free = ', idx_free(:)

        iter = 0

        ! Set initial error
        if(t%trim%verbose) then 
          write(t%iunit_trim,*)
          write(t%iunit_trim,*) 'Initial Guess'
        end if         

        res = calc_r(t, x, n_free)
        error = maxval(abs(res))

        if(t%trim%verbose) then 
          write(t%iunit_trim,*)
          write(t%iunit_trim,*) 'Beginning Solution Process...'
        end if 

        ! BEGIN NEWTONS METHOD
        ! Use central difference method to find jacobian
        do while(abs(error) > t%trim%solver%tolerance)
          iter = iter + 1

          if(t%trim%verbose) then 
            write(t%iunit_trim,*)
            write(t%iunit_trim,*) '------------------------------ '
            write(t%iunit_trim,*) 'Beginning iteration ', iter 
            write(t%iunit_trim,*) '------------------------------ '
            write(t%iunit_trim,*) 
            write(t%iunit_trim,*) 'Building Jacobian'
            write(t%iunit_trim,*)
          end if           

          do i = 1,n_free      
            k = idx_free(i) 

            if(t%trim%verbose) then 
              write(t%iunit_trim, '(A,I0,A)') 'Computing gradient relative to x[', k, ']'
              write(t%iunit_trim, '(A)') '   Positive Finite-Difference Step '
            end if 

            x(k) = x(k) + t%trim%solver%step_size
            R_pos = calc_r(t, x, n_free)
            
            if (t%trim%verbose) then 
              write(t%iunit_trim, '(A)') '   Negative Finite-Difference Step'
            end if 
                          
            x(k) = x(k) - 2.0 * t%trim%solver%step_size
            R_neg = calc_r(t, x, n_free)
            
            if(t%trim%verbose) then 
              write(t%iunit_trim,*) 
            end if

            x(k) = x(k) + t%trim%solver%step_size
            jac(:,i) = (R_pos - R_neg) / 2.0 / t%trim%solver%step_size
          end do 

          if (t%trim%verbose) then
            write(t%iunit_trim, '(A)') 'Jacobian Matrix ='
            do i = 1, size(jac,1)
                write(t%iunit_trim,'(*(1X,G25.17))') (jac(i,j), j=1,size(jac,2))
            end do
          end if

          res = calc_r(t, x, n_free)

          ! Calculate delta x and add relaxation factor
          call lu_solve(n_free, jac, -res, delta_x)

          i = 1
          do i=1,n_free 
            k = idx_free(i) 
            x(k) = x(k) + t%trim%solver%relaxation_factor * delta_x(i)
          end do 

          res = calc_r(t,x, n_free) 
          error = maxval(abs(res))
          
          if (t%trim%verbose) then
            write(t%iunit_trim,*)
            write(t%iunit_trim, '(A,*(1X,G25.17))') ' Delta X =', (delta_x(k), k=1,n_free)
            write(t%iunit_trim, '(A,*(1X,G25.17))') '   New X = ', x(:)
            write(t%iunit_trim, '(A,*(1X,G25.17))') 'Residual = ', res(:)
            write(t%iunit_trim, '(A,*(1X,G25.17))') '   Error = ', error 
          end if 

        end do 

        if (t%trim%verbose) then
          write(t%iunit_trim,*)
          write(t%iunit_trim,*)
          write(t%iunit_trim,*) 'Trim Solution:'
          write(t%iunit_trim,*) '--------------------------'
          write(t%iunit_trim,*) 'max error = ', error
          write(t%iunit_trim, '(A,*(1X,G25.17))') ' bank angle[deg]          = ', x(7) * 180.0 / pi 
          write(t%iunit_trim, '(A,*(1X,G25.17))') ' elevation angle[deg]     = ', x(8) * 180.0 / pi 
          write(t%iunit_trim, '(A,*(1X,G25.17))') ' angle of attack[deg]     = ', x(1) * 180.0 / pi 
          write(t%iunit_trim, '(A,*(1X,G25.17))') ' sideslip angle[deg]      = ', x(2) * 180 / pi      

          write(t%iunit_trim, '(A,*(1X,G25.17))') ' p[deg/sec]               = ', t%init_state(4) * 180.0 / pi 
          write(t%iunit_trim, '(A,*(1X,G25.17))') ' q[deg/sec]               = ', t%init_state(5) * 180.0 / pi 
          write(t%iunit_trim, '(A,*(1X,G25.17))') ' r[deg/sec]               = ', t%init_state(6) * 180.0 / pi          

          write(t%iunit_trim, '(A,*(1X,G25.17))') ' aileron deflection[deg]  = ', x(3) * 180.0 / pi 
          write(t%iunit_trim, '(A,*(1X,G25.17))') ' elevator deflection[deg] = ', x(4) * 180.0 / pi 
          write(t%iunit_trim, '(A,*(1X,G25.17))') ' rudder deflection[deg]   = ', x(5) * 180.0 / pi 
          write(t%iunit_trim, '(A,*(1X,G25.17))') ' throttle[none]           = ', x(6)      

          write(t%iunit_trim, '(A,*(1X,G25.17))') ' azimuth   angle[deg]     = ', x(9) * 180.0 / pi
          write(t%iunit_trim,*)

        end if
        
        t%init_eul(:) = x(7:9)
        t%init_state(14:17) = x(3:6)
        t%init_state(18:24) = 0.0 

        t%controls(1)%commanded_value = t%state(14)
        t%controls(2)%commanded_value = t%state(15)
        t%controls(3)%commanded_value = t%state(16)
        t%controls(4)%commanded_value = t%state(17)

        t%limit_controls = .true. 
        t%atm%use_turb = temp_turb  

      end subroutine

    !----------------------------------------
    ! Calculate Residual
      function calc_r(t, x, n_free) result(ans)
        implicit none
        type(vehicle_t) :: t
        real, intent(in) :: x(9)
        integer, intent(in) :: n_free
        integer :: last
        real :: FM(9) 
        real :: g, ac, xyzdot(3)
        real :: ca, cb, sa, sb, cp, sp, ct, st
        real :: u, v, w, euler(3)
        real :: angular_rates(3)
        real :: ans(n_free), temp_state(24), dummy_res(24)

        ans = 0.0 

        ! Pull out controls
        t%state(14:17) = x(3:6)
        
        ! Pull out euler angles 
        euler = x(7:9)

        ! Pull out alpha, beta, phi, and theta
        ca = cos(x(1))
        sa = sin(x(1))     
        cb = cos(x(2))
        sb = sin(x(2))
        cp = cos(x(7))
        sp = sin(x(7))
        ct = cos(x(8))
        st = sin(x(8))

        ! Pull out velocities
        u = t%init_airspeed * ca * cb
        v = t%init_airspeed * sb 
        w = t%init_airspeed * sa * cb 
        
        ! Set init_states of vehicle based on passed in values
        t%init_state(1)  = u
        t%init_state(2)  = v
        t%init_state(3)  = w
        t%init_state(10:13) = euler_to_quat(euler)

        ! Calculate gravity and gravity relief
        g = gravity_English(t%init_alt)
        xyzdot = quat_dependent_to_base(t%init_state(1:3), t%init_state(10:13))
        ac = (xyzdot(1)**2 + xyzdot(2)**2) / (earth_radius_ft - t%init_state(9))

        ! Caculate angular rates
        angular_rates = 0.0 
        if (t%trim%type == 'sct') then 
          angular_rates = (g-ac)*sp*ct / (u*ct*cp + w*st) * (/-st, sp*ct, cp*ct/)
        end if 

        ! Set states
        temp_state = 0.0
        temp_state(1:3)   = t%init_airspeed * (/ca*cb, sb, sa*cb/) 
        temp_state(4:6)   = angular_rates(:)
        temp_state(9)     = -t%init_alt
        temp_state(7:8)   = 0.0
        temp_state(10:13) = euler_to_quat(euler)
        
        ! Set controls 
        temp_state(14:17) = x(3:6) 

        ! Calculate residual
        dummy_res = diff_eq(t, 0.0, temp_state)
        ans = dummy_res(1:6)
        last = 6

        if(t%trim%solve_relative_climb_angle) then
          last = last + 1
          ans(last) = t%trim%climb_angle - calc_relative_climb_angle(temp_state)
        end if 

        if(t%trim%solve_load_factor) then 
          last = last + 1
          FM = pseudo_aero(t, temp_state)  
          ans(last) = t%trim%load_factor - ((FM(1)*sa - FM(3)*ca) / (t%mass*(g-ac)))
        end if 

        if (last /= n_free) then
          write(*,*) 'ERROR: residual size mismatch, last = ', last, ', n_free = ', n_free
          stop
        end if        

        if (t%trim%verbose) then 
          write(t%iunit_trim, '(A,*(1X,G25.17))') '       x =', x
          if(t%trim%type == 'sct') then 
            write(t%iunit_trim, *) '         p[deg/s] = ', angular_rates(1) * 180.0 / pi 
            write(t%iunit_trim, *) '         q[deg/s] = ', angular_rates(2) * 180.0 / pi 
            write(t%iunit_trim, *) '         r[deg/s] = ', angular_rates(3) * 180.0 / pi 
          end if 
            write(t%iunit_trim, '(A,*(1X,G25.17))') '       R =', ans
        end if 
        ! write(*,*) temp_state
        t%init_state(4:6) = temp_state(4:6)     
      end function calc_r

      function calc_relative_climb_angle(y) result(ans)
        implicit none 
        real :: y(24), ans 
        real :: xyzdot(3), Vmag 

        xyzdot = quat_dependent_to_base(y(1:3), y(10:13))
        Vmag = sqrt(y(1)**2 + y(2)**2 + y(3)**2)
        ans = asin(-xyzdot(3) / Vmag) 
      end function      
    !-------------------------
    ! Tick a vehicle forward in time 
      subroutine vehicle_tick_state(t, time, dt)
        implicit none 
        type(vehicle_t) :: t 
        real, intent(in) :: time, dt 
        real :: x1(24), x(24), out(11)
        real :: r(3), q(3), c(3), position(3)
        integer :: N
        real, allocatable :: waypoints(:,:)
        real :: rho, lambda, radius
        real :: hc, chic, cmd(2)      
        integer :: i, flag 
        
        ! allocate(waypoints(3,N))

        x = t%state
        ! Step vehicle forward in time
        x1 = rk4(t, time, t%state, dt) 

        do i = 1,4
          x1(13+i) = max(min(x1(13+i), t%controls(i)%mag_limit(2)),  t%controls(i)%mag_limit(1))
          x1(17+i) = max(min(x1(17+i), t%controls(i)%rate_limit(2)), t%controls(i)%rate_limit(1))
        end do 

        t%state = x1 

        ! Call geographic model 
        if(geographic_model_ID > 0) call spherical_earth(t, x, x1)

        ! Normalize quaternion
        call quat_norm(t%state(10:13)) 

        ! ! Path follower
        ! out = follow_wpp_fillet(waypoints, n, position, radius)

        ! flag   = int(out(1))
        ! r      = out(2:4)
        ! q      = out(5:7)
        ! c      = out(8:10)
        ! rho    = out(11)

        ! if (flag == 1) then
        !     cmd = straight_line_follow(t, r, q)
        ! else
        !     cmd = follow_orbit(t, c,rho,lambda)
        ! end if

        ! hc   = cmd(1)
        ! chic = cmd(2)
        
        call get_controller_input(t, time+dt) 
        
      end subroutine

    !----------------------------------------   
    ! Get a controller input
      subroutine get_controller_input(t,time) 
        implicit none 
        type(vehicle_t) :: t 
        real, intent(in) :: time 
        real :: controls_setpoint(4) 
        integer :: i 

        if(t%controller%running) then 
          controls_setpoint(:) = controller_update(t%controller, t, t%state, time) 
          ! controls_setpoint(:) = dynamic_inversion(t, time, t%state, [0.0, 0.0, 0.0, 350.0])
          ! write(*,*) 'controls_setpoint = ', controls_setpoint(1:3)*180.0/pi, controls_setpoint(4)
          ! write(*,*)
          do i = 1,4 
            t%controls(i)%commanded_value = controls_setpoint(i)
            if(t%controls(i)%dynamics_order == 0) t%state(13+i) = max(min(t%controls(i)%commanded_value, t%controls(i)%mag_limit(2)), t%controls(i)%mag_limit(1))
          end do 
        end if 

      end subroutine get_controller_input
    !----------------------------------------   
    ! Spherical Earth model 
      subroutine spherical_earth(t, y1, y2)
        implicit none 
        type(vehicle_t) :: t 
        real, intent(in) :: y1(13)
        real, intent(inout) :: y2(13) 
        real :: d, dxf, dyf, dzf 
        real :: phi1, psi1 
        real :: dpsi_g, psi_g1 
        real :: theta, H1
        real :: xhat, yhat, zhat, xhat_prime, yhat_prime, zhat_prime
        real :: rhat, Chat, Shat
        real :: cphi1, sphi1, ct, st, cg1, sg1
        real :: cg, sg, quat(4)        
        write(*,*) 'spherical earth'

        ! Pull out change in coordinate
        dxf = y2(7) - y1(7)
        dyf = y2(8) - y1(8)
        dzf = y2(9) - y1(9)

        d = sqrt(dxf**2 + dyf**2)
        if (d < tol) then 
          ! nothing 
        else 
          ! Pull out lat and long 
          H1   = -y1(9)
          phi1 = t%latitude
          psi1 = t%longitude
          
          cphi1 = cos(phi1)
          sphi1 = sin(phi1) 

          theta = d / (earth_radius_ft + H1 - dzf*0.5) 

          ct    = cos(theta) 
          st    = sin(theta) 
          
          psi_g1 = atan2(dyf,dxf)
          cg1    = cos(psi_g1)
          sg1    = sin(psi_g1)

          xhat = cphi1*ct - sphi1*st*cg1
          yhat = st*sg1
          zhat = sphi1*ct + cphi1*st*cg1

          xhat_prime = -cphi1*st - sphi1*ct*cg1
          yhat_prime = ct * sg1
          zhat_prime = -sphi1 * st + cphi1 * ct * cg1
          rhat       = sqrt(xhat**2 + yhat**2) 

          t%latitude  = atan2(zhat, rhat) 
          t%longitude = psi1 + atan2(yhat, xhat) 

          Chat   = xhat**2 * zhat_prime 
          Shat   = (xhat*yhat_prime - yhat * xhat_prime) * (cos(t%latitude))**2 * (cos(t%longitude-psi1))**2 
          dpsi_g = atan2(Shat, Chat) - psi_g1
          
          ! Limit  geographic coordinates
          if(t%longitude >  pi) t%longitude = t%longitude - 2.0*pi 
          if(t%longitude < -pi) t%longitude = t%longitude + 2.0*pi 

          ! Rotate flat earth quat according to delta bearing (7.5.17)
          cg = cos(0.5 * dpsi_g)
          sg = sin(0.5 * dpsi_g)
          quat(1) = -t%state(13)
          quat(2) = -t%state(12)
          quat(3) =  t%state(11)
          quat(4) =  t%state(10)
          t%state(10:13) = cg * t%state(10:13) + sg*quat(:)              
        end if 

      end subroutine 

    !----------------------------------------         
  !==================================================
  ! INTEGRATOR AND DIFFERENTIAL EQUATIONS
  !==================================================
    !----------------------------------------
    ! RK4 Integrator
      function rk4(t, t0, x1, delta_t) result(state)
        implicit none
        type(vehicle_t) :: t
        real, intent(in) :: t0, delta_t, x1(24)
        real, dimension(24) :: state, k1, k2, k3, k4

        if(t%rk4_verbose) then 
          write(t%iunit_rk4,*)
          write(t%iunit_rk4,*) '  state of the vehicle at the beginning of this RK4 integration step:'
          write(t%iunit_rk4,*) "time[s]             u[ft/s]              &
            &v[ft/s]              w[ft/s]              p[rad/s]             &
            &q[rad/s]             r[rad/s]             xf[ft]               &
            &yf[ft]               zf[ft]               e0                   &
            &ex                   ey                   ez"

          write(t%iunit_rk4,'(X,22(ES19.12,1X))') t0, t%state
          write(t%iunit_rk4,*)
          write(t%iunit_rk4,*) '  rk4 function called...'
          write(t%iunit_rk4,*)
        end if 

        ! K terms for RK4
        if (t%rk4_verbose) write(t%iunit_rk4,'(A,X,ES19.12,1X)')'    |           RK4 call number =             1'
        k1 = diff_eq(t, t0, x1)
        if (t%rk4_verbose) write(t%iunit_rk4,'(A,X,ES19.12,1X)')'    |           RK4 call number =             2'
        k2 = diff_eq(t, t0 + delta_t*0.5, x1 + k1 * delta_t*0.5)
        if (t%rk4_verbose) write(t%iunit_rk4,'(A,X,ES19.12,1X)')'    |           RK4 call number =             3'
        k3 = diff_eq(t, t0 + delta_t*0.5, x1 + k2 * delta_t*0.5)
        if (t%rk4_verbose) write(t%iunit_rk4,'(A,X,ES19.12,1X)')'    |           RK4 call number =             4'
        k4 = diff_eq(t, t0 + delta_t, x1 + k3 * delta_t)

        ! RK4 result
        state = x1 + (delta_t/6) * (k1 + 2*k2 + 2*k3 + k4)

        if (t%rk4_verbose) then 
          write(t%iunit_rk4,*) '  state of the vehicle after running RK4:'
          write(t%iunit_rk4,*) 'time[s]             u[ft/s]              &
            &v[ft/s]              w[ft/s]              p[rad/s]             &
            &q[rad/s]             r[rad/s]             xf[ft]               &
            &yf[ft]               zf[ft]               e0                   &
            &ex                   ey                   ez'
          write(t%iunit_rk4,'(X,22(ES19.12,1X))') t0+delta_t, state
          write(t%iunit_rk4,*) ' --------------------------- End of single RK4 integration step. ---------------------------'
        end if 

      end function rk4
    
    !----------------------------------------
    ! Differential Equations
      function diff_eq(t, time, state) result(dstate_dt)
        implicit none 
        type(vehicle_t), intent(inout) :: t
        real :: time
        real, target :: state(24) 
        real :: FMh(9) 
        real :: dstate_dt(24) 
        real :: acceleration(3), angular_accelerations(3), rhs(3), velocity(3), quat_change(4) 
        real :: quat_inv(4) 
        real :: orientation_effect(3), angular_v_effect(3), hmat_dot(3)
        real :: pqr_term(3), wind_velocity(3), hmat(3,3) 
        real :: quat_matrix(4,3) 
        real :: hxb_dot, hyb_dot, hzb_dot
        real :: Vxw, Vyw, Vzw, gravity_ft_per_sec2, ac
        real :: Imat(3,3), Imat_inv(3,3)
        real :: Ixx, Iyy, Izz, Ixy, Ixz, Iyz 
        real :: hxb, hyb, hzb        
        real, pointer :: u, v, w, p, q, r, e0, ex, ey, ez

        integer :: i 
        real :: wn, zeta
        real :: delta, d_delta, dd_delta

        if (t%rk4_verbose) then 
          write(t%iunit_rk4,'(A,X,ES19.12,1X)') '    |                  time [s] = ', time 
          write(t%iunit_rk4,'(A,X,24(ES19.12,1X))') '    |    state vector coming in = ', state 
        end if 
        ! Unpack states
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

        ! Unpack inertia
        Ixx = t%inertia(1,1)
        Iyy = t%inertia(2,2)
        Izz = t%inertia(3,3)
        Ixy = t%inertia(1,2)
        Ixz = t%inertia(1,3)
        Iyz = t%inertia(2,3)
      
        ! Calculate forces and moments
        FMh = pseudo_aero(t, state)
        hxb = t%h(1) + FMh(7) 
        hyb = t%h(2) + FMh(8) 
        hzb = t%h(3) + FMh(9) 

        gravity_ft_per_sec2 = gravity_English(-state(9))

        ! Set gyroscopic change and wind velocity to zero
        hxb_dot = 0.0
        hyb_dot = 0.0
        hzb_dot = 0.0
        Vxw = 0.0
        Vyw = 0.0
        Vzw = 0.0

        ! Build matrices/vectors in diff eqs
        orientation_effect(1) = 2 * (ex*ez - ey*e0)
        orientation_effect(2) = 2 * (ey*ez + ex*e0)
        orientation_effect(3) = ez**2 + e0**2 - ex**2 - ey**2

        angular_v_effect(1) = r*v - q*w 
        angular_v_effect(2) = p*w - r*u
        angular_v_effect(3) = q*u - p*v

        hmat(1,1) =  0.0
        hmat(1,2) = -hzb
        hmat(1,3) =  hyb

        hmat(2,1) =  hzb
        hmat(2,2) =  0.0
        hmat(2,3) = -hxb

        hmat(3,1) = -hyb
        hmat(3,2) =  hxb
        hmat(3,3) =  0.0

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

        hmat_dot(1) = hxb_dot
        hmat_dot(2) = hyb_dot
        hmat_dot(3) = hzb_dot

        pqr_term(1) = (Iyy - Izz)*q*r + Iyz*(q**2 -r**2) + Ixz*p*q - Ixy*p*r
        pqr_term(2) = (Izz - Ixx)*p*r + Ixz*(r**2 -p**2) + Ixy*q*r - Iyz*p*q
        pqr_term(3) = (Ixx - Iyy)*p*q + Ixy*(p**2 -q**2) + Iyz*p*r - Ixz*q*r

        wind_velocity(1) = Vxw
        wind_velocity(2) = Vyw
        wind_velocity(3) = Vzw

        quat_inv(1) =  e0
        quat_inv(2) = -ex
        quat_inv(3) = -ey
        quat_inv(4) = -ez

        quat_matrix(1,1) = -ex
        quat_matrix(1,2) = -ey
        quat_matrix(1,3) = -ez

        quat_matrix(2,1) =  e0
        quat_matrix(2,2) = -ez
        quat_matrix(2,3) =  ey

        quat_matrix(3,1) =  ez
        quat_matrix(3,2) =  e0
        quat_matrix(3,3) = -ex

        quat_matrix(4,1) = -ey
        quat_matrix(4,2) =  ex
        quat_matrix(4,3) =  e0           

        ! Differential equations
        ! Roll, pitch, yaw accel
        rhs(1) = FMh(4) + hmat(1,1)*p + hmat(1,2)*q + hmat(1,3)*r + pqr_term(1) - hmat_dot(1)
        rhs(2) = FMh(5) + hmat(2,1)*p + hmat(2,2)*q + hmat(2,3)*r + pqr_term(2) - hmat_dot(2)
        rhs(3) = FMh(6) + hmat(3,1)*p + hmat(3,2)*q + hmat(3,3)*r + pqr_term(3) - hmat_dot(3)

        angular_accelerations(1) = Imat_inv(1,1)*rhs(1) + Imat_inv(1,2)*rhs(2) + Imat_inv(1,3)*rhs(3)
        angular_accelerations(2) = Imat_inv(2,1)*rhs(1) + Imat_inv(2,2)*rhs(2) + Imat_inv(2,3)*rhs(3)
        angular_accelerations(3) = Imat_inv(3,1)*rhs(1) + Imat_inv(3,2)*rhs(2) + Imat_inv(3,3)*rhs(3)

        ! Velocity in inertial
        velocity = quat_base_to_dependent(state(1:3), state(10:13)) + wind_velocity

        ! Gravity relief
        ac = (velocity(1)**2 + velocity(2)**2) / (earth_radius_ft - state(9))

        ! Aacceleration in body
        acceleration(1) = FMh(1)/t%mass + (gravity_ft_per_sec2 - ac)*orientation_effect(1) + angular_v_effect(1)
        acceleration(2) = FMh(2)/t%mass + (gravity_ft_per_sec2 - ac)*orientation_effect(2) + angular_v_effect(2)
        acceleration(3) = FMh(3)/t%mass + (gravity_ft_per_sec2 - ac)*orientation_effect(3) + angular_v_effect(3)

        ! Orientation rate of change
        quat_change(1) = 0.5 * (quat_matrix(1,1)*p + quat_matrix(1,2)*q + quat_matrix(1,3)*r)
        quat_change(2) = 0.5 * (quat_matrix(2,1)*p + quat_matrix(2,2)*q + quat_matrix(2,3)*r)
        quat_change(3) = 0.5 * (quat_matrix(3,1)*p + quat_matrix(3,2)*q + quat_matrix(3,3)*r)
        quat_change(4) = 0.5 * (quat_matrix(4,1)*p + quat_matrix(4,2)*q + quat_matrix(4,3)*r)

        ! State derivative
        dstate_dt(1:3)  = acceleration
        dstate_dt(4:6)  = angular_accelerations
        dstate_dt(7:9)  = velocity
        dstate_dt(10:13) = quat_change

        ! Actuator dynamics
          do i = 1,4
            ! Read in state and limit magnitude 
            delta = state(13 + i)
            d_delta = state(17 + i) 
            delta = max(min(delta, t%controls(i)%mag_limit(2)),  t%controls(i)%mag_limit(1))

            ! Compute dynamics
            select case (t%controls(i)%dynamics_order)

            case(0) ! No dynamics
              d_delta = 0.0 
              dd_delta = 0.0 

            case(1) ! First order dynamics             
              d_delta = (t%controls(i)%commanded_value - delta) / t%controls(i)%time_constant
              d_delta = max(min(d_delta, t%controls(i)%rate_limit(2)),  t%controls(i)%rate_limit(1))

              ! Position limit 
              if (delta <= t%controls(i)%mag_limit(1) + tol .and. d_delta < 0.0) d_delta = 0.0 
              if (delta >= t%controls(i)%mag_limit(2) - tol .and. d_delta > 0.0) d_delta = 0.0 

              dd_delta = 0.0 
            
            case(2) ! Second order dynamics 
              wn = t%controls(i)%natural_frequency 
              zeta = t%controls(i)%damping_ratio
              d_delta = max(min(d_delta, t%controls(i)%rate_limit(2)),  t%controls(i)%rate_limit(1))

              if (delta <= t%controls(i)%mag_limit(1) + tol .and. d_delta < 0.0) d_delta = 0.0 
              if (delta >= t%controls(i)%mag_limit(2) - tol .and. d_delta > 0.0) d_delta = 0.0 

              dd_delta = wn**2 * (t%controls(i)%commanded_value - delta) - 2.0*zeta*wn*d_delta
              dd_delta = max(min(dd_delta, t%controls(i)%accel_limit(2)),  t%controls(i)%accel_limit(1))
              
              ! Rate limit 
              if (d_delta <= t%controls(i)%rate_limit(1) + tol .and. dd_delta < 0.0) dd_delta = 0.0 
              if (d_delta >= t%controls(i)%rate_limit(2) - tol .and. dd_delta > 0.0) dd_delta = 0.0 

              ! Position limit
              if (delta <= t%controls(i)%mag_limit(1) + tol .and. t%controls(i)%commanded_value <= t%controls(i)%mag_limit(1) + tol) dd_delta = 0.0 
              if (delta >= t%controls(i)%mag_limit(2) - tol .and. t%controls(i)%commanded_value >= t%controls(i)%mag_limit(2) - tol) dd_delta = 0.0 
            
            end select 

            dstate_dt(13+i) = d_delta 
            dstate_dt(17+i) = dd_delta 

          end do 

          ! Integral error
          dstate_dt(22:24) = t%zdot 

        if (t%rk4_verbose) then 
          write(t%iunit_rk4,'(A,X,6(ES19.12,1X))')  '    | pseudo aerodynamics (F,M) = ', FMh
          write(t%iunit_rk4,'(A,X,24(ES19.12,1X))') '    |           diff_eq results = ', dstate_dt
        end if 

      end function diff_eq

    !----------------------------------------
  !==================================================
  ! AERODYNAMICS AND FORCES
  !==================================================
    !----------------------------------------
    ! Aerodynamic Forces and Moments for f16
      function pseudo_aero(t, state) result(FMh)
        implicit none
        type(vehicle_t) :: t
        real, intent(in) :: state(24)
        real :: local_state(24) 
        real :: FMh(9) 
        real :: Re, geometric_altitude_ft, geopotential_altitude_ft
        real :: temp_R, pressure_lbf_per_ft2, density_slugs_per_ft3
        real :: dyn_viscosity_slug_per_ft_sec, sos_ft_per_sec
        real :: V, dyn_pressure
        real :: alpha, beta, beta_f, pbar, qbar, rbar
        real :: CL, CL1, CD, CS, L, D, S, Cl_roll, Cm, Cn
        real :: CM1, CM2, mach_num
        real :: ca, cb, sa, sb
        real :: alpha_hat, beta_hat
        real :: delta_a, delta_e, delta_r, tau 
        real :: CL_s, CD_s, Cm_s 
        real :: sigma_D, sigma_L, sigma_m, sign_a
        real :: turbulence(6)
        real :: blend, drag_factor, lift_factor, moment_factor
        real :: pg_factor, kt_factor
        real :: m_high, m_low, m_trans_start
        real :: sqrt_term
        real :: max_CL_factor, max_CD_factor, max_Cl_roll_factor
        real :: alphad, betad, ded, speedbrake, lef, Cxyzlmn(6) 
        real :: db1(1), db2(2), db3(3), db6(6)
        integer :: i, throttle_ID

        local_state = state 
        
        ! Build atmosphere
        geometric_altitude_ft = -state(9)
        call std_atm_English(&
          geometric_altitude_ft, geopotential_altitude_ft,     & 
          temp_R, pressure_lbf_per_ft2, density_slugs_per_ft3, & 
          dyn_viscosity_slug_per_ft_sec, sos_ft_per_sec)

        ! Add turbulence
        if (t%atm%use_turb) then
          turbulence = atmosphere_get_turbulence(t%atm, state)
          local_state(1:6) = state(1:6) + turbulence
        end if 

        ! Calculate velocity unit vector
        V            =  sqrt(local_state(1)**2 + local_state(2)**2 + local_state(3)**2)
        dyn_pressure = 0.5 * density_slugs_per_ft3 * V **2 * t%planform_area

        ! Calculate alpha and beta 3.4.4 and 3.4.5
        alpha  = atan2(local_state(3) , local_state(1))
        beta   = asin(local_state(2) / V)
        beta_f = atan2(local_state(2) , local_state(1))

        ca = cos(alpha)
        cb = cos(beta)
        sa = sin(alpha)
        sb = sin(beta)

        ! Define alpha_hat and beta_hat
        alpha_hat = 0.0
        beta_hat  = 0.0

        ! Calculate rotation rates from eq 3.4.23
        pbar = 1 / (2*V) * local_state(4) * t%lateral_length
        qbar = 1 / (2*V) * local_state(5) * t%longitudinal_length
        rbar = 1 / (2*V) * local_state(6) * t%lateral_length

        if(t%type == 'aircraft') then 
          ! Pull out controls
          if(t%limit_controls) then 
            delta_a  = max(t%controls(1)%mag_limit(1), min(t%controls(1)%mag_limit(2), state(t%aileron_ID)))
            delta_e  = max(t%controls(2)%mag_limit(1), min(t%controls(2)%mag_limit(2), state(t%elevator_ID)))
            delta_r  = max(t%controls(3)%mag_limit(1), min(t%controls(3)%mag_limit(2), state(t%rudder_ID)))
          else 
            delta_a  = state(t%aileron_ID)
            delta_e  = state(t%elevator_ID)
            delta_r  = state(t%rudder_ID)
          end if 

          ! Database aerodynamics
          Cxyzlmn = 0.0 
          if(t%use_database) then 
            if(allocated(t%db)) then 
              alphad = alpha * 180.0 / pi 
              betad = beta * 180.0 / pi 
              ded = delta_e * 180.0 / pi 
              speedbrake = t%speed_brake 
              lef = t%le_flap 

              ! C(x,y,z,l,m,n)o(elevator,alpha,beta)
              db6 = t%db(1)%interpolate([ded,alphad,betad])
              Cxyzlmn = Cxyzlmn + db6 

              ! DC(m)(alpha)
              db1 = t%db(2)%interpolate([alphad])
              Cxyzlmn(5) = Cxyzlmn(5) + db1(1) 

              ! DC(m)(elevator,alpha) 
              db1 = t%db(3)%interpolate([ded,alphad])
              Cxyzlmn(5) = Cxyzlmn(5) + db1(1) 

              ! DC(l,n),beta(alpha) 
              db2(:) = t%db(4)%interpolate([alphad])
              Cxyzlmn(4) = Cxyzlmn(4) + db2(1)*beta
              Cxyzlmn(6) = Cxyzlmn(6) + db2(2)*beta 

              ! DC(x,y,z,l,m,n,leg(alpha,beta) 
              db6 = t%db(5)%interpolate([alphad,betad]) 
              Cxyzlmn = Cxyzlmn + db6*lef 

              ! DC(x,z,m),qbar_lef(alpha) 
              db3 = t%db(6)%interpolate([alphad])
              Cxyzlmn(1) = Cxyzlmn(1) + db3(1)*qbar*lef 
              Cxyzlmn(3) = Cxyzlmn(3) + db3(2)*qbar*lef 
              Cxyzlmn(5) = Cxyzlmn(5) + db3(3)*qbar*lef 

              ! DC(s,z,m),qbar(alpha) 
              db3 = t%db(7)%interpolate([alphad])
              Cxyzlmn(1) = Cxyzlmn(1) + db3(1)*qbar
              Cxyzlmn(3) = Cxyzlmn(3) + db3(2)*qbar
              Cxyzlmn(5) = Cxyzlmn(5) + db3(3)*qbar

              ! DC(x,z,m),sb(alpha)
              db3 = t%db(8)%interpolate([alphad])
              Cxyzlmn(1) = Cxyzlmn(1) + db3(1)*speedbrake 
              Cxyzlmn(3) = Cxyzlmn(3) + db3(2)*speedbrake 
              Cxyzlmn(5) = Cxyzlmn(5) + db3(3)*speedbrake 

              ! DC(y,l,n),aileron_lef(alpha,beta)
              db3 = t%db(9)%interpolate([alphad,betad])
              Cxyzlmn(2) = Cxyzlmn(2) + db3(1)*delta_a*lef 
              Cxyzlmn(4) = Cxyzlmn(4) + db3(2)*delta_a*lef 
              Cxyzlmn(6) = Cxyzlmn(6) + db3(3)*delta_a*lef 

              ! DC(y,l,n),aileron(alpha,beta)
              db3 = t%db(10)%interpolate([alphad,betad])
              Cxyzlmn(2) = Cxyzlmn(2) + db3(1)*delta_a
              Cxyzlmn(4) = Cxyzlmn(4) + db3(2)*delta_a
              Cxyzlmn(6) = Cxyzlmn(6) + db3(3)*delta_a
              
              ! DC(y,l,n),pbar_lef(alpha)
              db3 = t%db(11)%interpolate([alphad])
              Cxyzlmn(2) = Cxyzlmn(2) + db3(1)*pbar*lef 
              Cxyzlmn(4) = Cxyzlmn(4) + db3(2)*pbar*lef 
              Cxyzlmn(6) = Cxyzlmn(6) + db3(3)*pbar*lef 
              
              ! DC(y,l,n),pbar(alpha)
              db3 = t%db(12)%interpolate([alphad])
              Cxyzlmn(2) = Cxyzlmn(2) + db3(1)*pbar
              Cxyzlmn(4) = Cxyzlmn(4) + db3(2)*pbar
              Cxyzlmn(6) = Cxyzlmn(6) + db3(3)*pbar
              
              ! DC(y,l,n),rbar_lef(alpha)
              db3 = t%db(13)%interpolate([alphad])
              Cxyzlmn(2) = Cxyzlmn(2) + db3(1)*rbar*lef 
              Cxyzlmn(4) = Cxyzlmn(4) + db3(2)*rbar*lef 
              Cxyzlmn(6) = Cxyzlmn(6) + db3(3)*rbar*lef 
              
              ! DC(y,l,n),rbar(alpha)
              db3 = t%db(14)%interpolate([alphad])
              Cxyzlmn(2) = Cxyzlmn(2) + db3(1)*rbar 
              Cxyzlmn(4) = Cxyzlmn(4) + db3(2)*rbar 
              Cxyzlmn(6) = Cxyzlmn(6) + db3(3)*rbar 
              
              ! DC(y,l,n),rudder(alpha,beta)
              db3 = t%db(15)%interpolate([alphad,betad])
              Cxyzlmn(2) = Cxyzlmn(2) + db3(1)*delta_r 
              Cxyzlmn(4) = Cxyzlmn(4) + db3(2)*delta_r 
              Cxyzlmn(6) = Cxyzlmn(6) + db3(3)*delta_r 
            end if  
          else 
            ! Force and moment coefficients 
            CL1 =  t%CL0 + t%CL_alpha * alpha
            CL  = CL1 + t%CL_qbar * qbar + t%CL_alphahat * alpha_hat + t%CL_elevator * delta_e
            CS  = t%CS_beta * beta + (t%CS_pbar + t%CS_alpha_pbar * alpha) * pbar + t%CS_rbar * rbar + t%CS_aileron * delta_a + t%CS_rudder * delta_r
            CD  =  t%CD_L0 + t%CD_L1 * CL1 + t%CD_L1_L1 * CL1 **2 + t%CD_CS_CS * CS **2 + (t%CD_qbar + t%CD_alpha_qbar * alpha) * qbar + (t%CD_elevator + t%CD_alpha_elevator * alpha) * delta_e + t%CD_elevator_elevator * delta_e ** 2
            Cl_roll = t%Cl_beta * beta + t%Cl_pbar * pbar + (t%Cl_rbar + t%Cl_alpha_rbar * alpha) * rbar + t%Cl_aileron * delta_a + t%Cl_rudder * delta_r  ! roll moment
            Cm      = t%Cm_0 + t%Cm_alpha * alpha + t%Cm_qbar * qbar + t%Cm_alphahat * alpha_hat + t%Cm_elevator * delta_e ! pitch moment
            Cn      = t%Cn_beta * beta + (t%Cn_pbar + t%Cn_alpha_pbar * alpha) * pbar + t%Cn_rbar * rbar + (t%Cn_aileron + t%Cn_alpha_aileron * alpha) * delta_a + t%Cn_rudder * delta_r ! yaw moment         
 
          end if                

        else if (t%type == 'arrow') then 
          CL = t%CL_alpha * alpha 
          CS = -t%CL_alpha * beta_f 
          CD = t%CD_L0 + t%CD_L1_L1*CL**2 + t%CD_L1_L1*CS**2 
          
          Cl_roll  = t%Cl_l0 + t%Cl_pbar*pbar
          Cm       = t%Cm_alpha*alpha + t%Cm_qbar*qbar 
          Cn       = -t%Cm_alpha*beta_f + t%Cm_qbar*rbar

        else if (t%type == 'sphere') then 
          CL      = 0.0 
          CS      = 0.0 
          Cl_roll = 0.0 
          Cm      = 0.0 
          Cn      = 0.0
          
          ! Calculate reynolds number
          Re = density_slugs_per_ft3 * V * 2 * t%longitudinal_length / dyn_viscosity_slug_per_ft_sec

          ! Compute drag coef
          if(0.0 < Re .and. Re <= 450000.0) then
            CD = 24.0/Re + 6.0/(1.0 + Re**0.5) + 0.4
          else if(450000.0 < Re .and. Re <= 560000.0) then 
            CD = 1.0E29 * Re **(-5.211)
          else if(560000 < Re .and. Re <= 14000000) then
            CD = -2.0E-23 * Re**3 - 1.0E-16*Re**2 + 9.0E-09 * Re + 0.069
          else if(Re > 1400000) then
            CD = 0.12
          end if

        end if 

        ! Stall model 
        if (t%stall) then 
          sign_a  = sign(1.0,alpha)
          CL_s    = 2 * sign_a * sa**2 * ca 
          CD_s    = 2 * (sin(abs(alpha)))**3
          sigma_L = calc_sigma(t%CL_lambda_b, t%CL_alpha_0, t%CL_alpha_s, alpha)
          sigma_D = calc_sigma(t%CD_lambda_b, t%CD_alpha_0, t%CD_alpha_s, alpha)
          
          CL = CL * (1 - sigma_L) + CL_s * sigma_L 
          CD = CD * (1 - sigma_D) + CD_s * sigma_D 

          Cm_s    = t%Cm_min * sa**2 * sign_a
          sigma_m = calc_sigma(t%Cm_lambda_b, t%Cm_alpha_0, t%Cm_alpha_s, alpha)
          Cm      = Cm * (1 - sigma_m) + Cm_s * sigma_m 
    
        end if 

        ! Compressibility model
        if (t%compressibility) then 
          CM1      = 2.13/ (t%sweep + 0.15)**2
          CM2      = 15.35*t%sweep**2 - 19.64*t%sweep +16.86
          mach_num = V / sos_ft_per_sec

          ! Mach breakpoints
          m_low         = 0.60    ! only use prandtl-glauert
          m_high        = 0.92    ! only use karman-tsien
          m_trans_start = 0.88    ! start transsonic region

          ! Prandtl-Glauert factor
          sqrt_term = sqrt(max(1.0 - mach_num**2, tol))
          pg_factor = 1.0 / sqrt_term

          ! Karman-Tsien factor
          kt_factor = pg_factor * (1.0 + mach_num**2 / (1.0 + sqrt_term))

          ! Blending 
          if (mach_num <= m_low) then
            blend = 0.0
          else if (mach_num >= m_high) then
            blend = 1.0
          else
            blend = (mach_num - m_low) / (m_high - m_low)
          end if

          ! Final compressibility factor
          lift_factor   = (1.0 - blend) * pg_factor + blend * kt_factor
          moment_factor = lift_factor
          drag_factor   = 1.0 + CM1 * mach_num**CM2

          ! Apply factors
          max_CL_factor      = 6.5
          max_CD_factor      = 4
          max_Cl_roll_factor = 2.5

          if (mach_num < 0.92) then ! accurate range
            CL       = CL * lift_factor
            CS       = CS * lift_factor
            Cl_roll  = Cl_roll * moment_factor
            Cm       = Cm * moment_factor
            Cn       = Cn * moment_factor
            CD       = CD * drag_factor
          else
            CL       = CL * min(lift_factor, max_CL_factor)
            CS       = CS * min(lift_factor, max_CL_factor)
            Cl_roll  = Cl_roll * min(moment_factor, max_Cl_roll_factor)
            Cm       = Cm * min(moment_factor, max_Cl_roll_factor)
            Cn       = Cn * min(moment_factor, max_Cl_roll_factor)
          end if

          ! Apply drag multiplier
          CD = CD * min(drag_factor, max_CD_factor)
        end if 

        if(t%use_database) then 
          FMh(1:3) = Cxyzlmn(1:3) 
          FMh(4)   = Cxyzlmn(4) * t%lateral_length 
          FMh(5)   = Cxyzlmn(5) * t%longitudinal_length
          FMh(6)   = Cxyzlmn(6) * t%lateral_length

          FMh(1:6) = FMh(1:6) * dyn_pressure
        else if(t%type .ne. 'quadrotor') then 
          L =   CL * dyn_pressure
          S =   CS * dyn_pressure
          D =   CD * dyn_pressure

          ! Table 3.4.4
          FMh(1) = - (ca*(D*cb + S*sb) - L*sa)
          FMh(2) = (S*cb - D*sb)
          FMh(3) = - (sa*(D*cb + S*sb) + L*ca)

          FMh(4) = Cl_roll  * dyn_pressure * t%lateral_length
          FMh(5) = Cm       * dyn_pressure * t%longitudinal_length
          FMh(6) = Cn       * dyn_pressure * t%lateral_length       
        end if 
        
        ! Add propulsion
        do i=1,t%num_props 
          throttle_ID = t%props(i)%control_ID 
          if(t%limit_controls) then 
            tau = max(t%controls(throttle_ID)%mag_limit(1), min(t%controls(throttle_ID)%mag_limit(2), state(throttle_ID+13)))
          else 
            tau = state(throttle_ID + 13) 
          end if 
          FMh(:) = FMh(:) + propulsion_get_FMh(t%props(i), state, tau)
        end do 

        ! Shift CG location 
        FMh(4:6) = FMh(4:6) + cross_product(t%aero_ref_location, FMh(1:3))
      end function pseudo_aero


    !----------------------------------------
    ! Calculate Sigma
      function calc_sigma(lambda_b, alpha_0, alpha_s, alpha) result(sigma)
        implicit none
        real, intent(in) :: lambda_b, alpha_0, alpha_s, alpha
        real :: sigma, neg, pos

        pos = exp( lambda_b * (alpha - alpha_0 + alpha_s))
        neg = exp(-lambda_b * (alpha - alpha_0 - alpha_s))
        sigma = (1.0 + neg + pos) / ((1.0 + neg)*(1.0 + pos))

      end function calc_sigma

    !----------------------------------------
  !==================================================
  ! WAYPOINT FINDING
  !==================================================
    !----------------------------------------
    ! Follow a straight line 
      function straight_line_follow(t, r, q) result(command)
      implicit none
      type(vehicle_t) :: t
      real, intent(in) :: r(3), q(3)
      real :: command(2)
      real :: chi_c, h_c
      real :: chi_q
      real :: chi_inf, k_path
      real :: pn, pe, pd
      real :: e_py
      real :: s

      chi_inf = 1.047      ! 60 deg
      k_path  = 0.05

      ! MAV position
      pn = t%state(1)
      pe = t%state(2)
      pd = t%state(3)

      ! Desired path course
      chi_q = atan2(q(2), q(1))

      ! Cross-track error
      e_py = -(pn-r(1))*sin(chi_q) + (pe-r(2))*cos(chi_q)

      ! Course command
      chi_c = chi_q - chi_inf*(2.0/pi)*atan(k_path*e_py)

      ! Along-track horizontal distance
      s = sqrt((pn-r(1))**2 + (pe-r(2))**2)

      ! Altitude command
      h_c = -r(3) - s * q(3)/sqrt(q(1)**2 + q(2)**2)

      command(1) = h_c
      command(2) = chi_c

      end function straight_line_follow
    !----------------------------------------
    ! Follow a circular orbit 
        function follow_orbit(t, c, rho, lambda) result(command)
        implicit none
        type(vehicle_t) :: t
        real, intent(in) :: c(3)
        real, intent(in) :: rho
        real, intent(in) :: lambda
        real :: command(2)
        real :: h_c, chi_c
        real :: pn, pe
        real :: d
        real :: phi
        real :: k_orbit

        ! Guidance gain
        k_orbit = 4.0

        ! Vehicle position
        pn = t%state(1)
        pe = t%state(2)

        ! Commanded altitude
        h_c = -c(3)

        ! Distance from orbit center
        d = sqrt((pn-c(1))**2 + (pe-c(2))**2)

        ! Angle from center to aircraft
        phi = atan2(pe - c(2), pn - c(1))

        ! Commanded course
        chi_c = phi + lambda*(pi/2.0 + atan(k_orbit*(d-rho)/rho))

        command(1) = h_c
        command(2) = chi_c

        end function follow_orbit 
    !----------------------------------------
    ! Follow filleted path 
      function follow_wpp_fillet(w, N, p, radius) result(out)
      implicit none
      integer, intent(in) :: N
      real, intent(in) :: w(3,N)
      real, intent(in) :: p(3)
      real, intent(in) :: radius
      real :: out(11)
      integer, save :: i = 2
      integer, save :: state = 1
      integer :: flag
      real :: r(3), q(3), c(3)
      real :: rho
      real :: lambda
      real :: qi_prev(3), qi(3)
      real :: z(3)
      real :: varrho
      real :: normv

      c = 0.0
      rho = 0.0

      ! unit vectors
      qi_prev = (w(:,i) - w(:,i-1))
      qi_prev = qi_prev / sqrt(sum(qi_prev**2))

      qi = (w(:,i+1) - w(:,i))
      qi = qi / sqrt(sum(qi**2))

      ! angle between segments
      varrho = acos(-dot_product(qi_prev,qi))

      if (state == 1) then
          flag = 1
          r = w(:,i-1)
          q = qi_prev
          z = w(:,i) - (radius/tan(varrho/2.0))*qi_prev
          if (dot_product(p-z, qi_prev) > 0.0) then
              state = 2
          end if
      else if (state == 2) then
          flag = 2
          normv = sqrt(sum((qi_prev-qi)**2))
          c = w(:,i) - (radius/sin(varrho/2.0)) * (qi_prev-qi)/normv
          rho = radius
          lambda = sign(1.0, qi_prev(1)*qi(2) - qi_prev(2)*qi(1))
          z = w(:,i) + (radius/tan(varrho/2.0))*qi
          if (dot_product(p-z, qi) > 0.0) then
              if (i < N-1) then
                  i = i + 1
              end if
              state = 1
          end if
      end if

      out(1) = flag
      out(2:4) = r
      out(5:7) = q
      out(8:10) = c
      out(11) = rho

      end function follow_wpp_fillet
    !----------------------------------------

end module vehicle_m