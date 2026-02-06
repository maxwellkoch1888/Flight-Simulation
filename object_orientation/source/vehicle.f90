module vehicle_m
    use koch_m
    use jsonx_m
    use micro_time_m
    use linalg_mod

    implicit none
    real :: rho0
    real, parameter :: earth_radius_ft = 6366707.01949371/0.3048
    integer :: geographic_model_ID
    character(len=:), allocatable :: geographic_model 
    logical :: save_lat_long

    !==================================================
    ! TYPE DECLARATIONS
    !==================================================
      !----------------------------------------
      ! Trim solver type
        type trim_solver_t
          real :: step_size, relaxation_factor, tolerance, max_iterations
        end type trim_solver_t
      !----------------------------------------
      ! Trim solver type
        type trim_settings_t
          character(len=:), allocatable :: type 

          logical, allocatable :: free_vars(:) 
          real :: climb_angle, sideslip_angle
          real :: load_factor
          logical :: verbose, solve_relative_climb_angle, solve_load_factor

          type(trim_solver_t) :: solver 
        end type trim_settings_t
      !----------------------------------------
      ! Vehicle type
      type vehicle_t
        type(json_value), pointer :: j_vehicle
            
        character(len=:), allocatable :: name
        character(len=:), allocatable :: type
        character(100) :: states_filename, rk4_filename, trim_filename, latlong_filename

        logical :: run_physics
        logical :: save_states
        integer :: iunit_states, iunit_rk4, iunit_trim, iunit_latlong

        ! Location 
        real :: latitude, longitude

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
        real :: init_airspeed, init_alt, init_state(13)
        real, allocatable :: init_eul(:) ! has to be allocatable because will be read from json object

        ! States/controls
        real :: initial_state(13), state(13)
        real :: controls(4)

        type(trim_settings_t) :: trim
      end type vehicle_t
      !----------------------------------------
    contains
    !==================================================
    ! INITIALIZATION FUNCTIONS
    !==================================================
      !----------------------------------------
      ! Vehicle initialization
        subroutine vehicle_init(t, j_vehicle_input)
          implicit none 
          type(vehicle_t), intent(inout) :: t
          type(json_value), pointer :: j_vehicle_input
          character(len=:), allocatable :: init_type 
          real :: geopotential_altitude_ft,temp_R, pressure_lbf_per_ft2, density_slugs_per_ft3, dyn_viscosity_slug_per_ft_sec, sos_ft_per_sec

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
              write(t%iunit_states,*) " time[s]             u[ft/s]              &
              &v[ft/s]              w[ft/s]              p[rad/s]             &
              &q[rad/s]             r[rad/s]             xf[ft]               &
              &yf[ft]               zf[ft]               e0                   &
              &ex                   ey                   ez                        &
              &aileron[deg]              elevator[deg]             rudder[deg]               throttle"
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
              write(*,*) '- saving RK4 results to ', t%latlong_filename
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
            call jsonx_get(t%j_vehicle, 'aerodynamics.reference.area[ft^2]',              t%planform_area)
            call jsonx_get(t%j_vehicle, 'aerodynamics.reference.longitudinal_length[ft]', t%longitudinal_length)
            call jsonx_get(t%j_vehicle, 'aerodynamics.reference.lateral_length[ft]',      t%lateral_length)
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
            write(*,*) '- thrust'
            call jsonx_get(t%j_vehicle, 'thrust.T0[lbf]',          t%T0, 0.0)
            call jsonx_get(t%j_vehicle, 'thrust.Ta',               t%Ta, 0.0)
            call jsonx_get(t%j_vehicle, 'thrust.location[ft]',     t%thrust_location, 0.0, 3)
            call jsonx_get(t%j_vehicle, 'thrust.orientation[deg]', t%thrust_orientation, 0.0, 3)
            
            t%thrust_orientation = t%thrust_orientation * pi / 180.0
            t%thrust_quat = euler_to_quat(t%thrust_orientation)

            ! Calculate rho0 for thrust
            call std_atm_English(0.0, geopotential_altitude_ft,     & 
              temp_R, pressure_lbf_per_ft2, density_slugs_per_ft3, & 
              dyn_viscosity_slug_per_ft_sec, sos_ft_per_sec)
            t%rho0 = density_slugs_per_ft3

            ! Initial conditions
            write(*,*) '- Initial Conditions'
            t%init_state = 0.0 
            t%controls = 0.0

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

          end if 
          write(*,*) 'Finished vehicle initialization.'
        end subroutine

      !----------------------------------------
      ! Write states to a file
        subroutine vehicle_write_state(t, time)
          implicit none 
          type(vehicle_t) :: t 
          real, intent(in) :: time 

          ! write(t%iunit_states,*) time, t%state, t%controls 
          write(t%iunit_states,'(14(ES20.13,1X))') time, t%state

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

          ! Calculate initial angular velocities 
          call jsonx_get(j_state, 'p[deg/s]', t%init_state(4))
          call jsonx_get(j_state, 'q[deg/s]', t%init_state(5))
          call jsonx_get(j_state, 'r[deg/s]', t%init_state(6))
          t%init_state(4:6) = t%init_state(4:6) * pi / 180.0 

          ! Calculate initial controls
          t%controls = 0.0 
          if(t%type == 'aircraft') then 
            call jsonx_get(j_state, 'aileron[deg]' , t%controls(1))
            call jsonx_get(j_state, 'elevator[deg]', t%controls(2))
            call jsonx_get(j_state, 'rudder[deg]'  , t%controls(3))
            call jsonx_get(j_state, 'throttle'     , t%controls(4))
            t%controls(1:3) = t%controls(1:3) * pi / 180.0 
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

          ! t%trim%solve_fixed_climb_angle =    .false. 
          ! t%trim%solve_relative_climb_angle = .false. 
          ! t%trim%solve_load_factor =          .false. 

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
            write(t%iunit_trim, '(A,*(1X,G25.17))') 'elevation angle[deg] = ', x(8) * 180.0 / pi 
            write(t%iunit_trim, '(A,*(1X,G25.17))') 'bank angle[deg]      = ', x(7) * 180.0 / pi 
            write(t%iunit_trim, '(A,*(1X,G25.17))') 'angle of attack[deg] = ', x(1) * 180.0 / pi 
            write(t%iunit_trim, '(A,*(1X,G25.17))') 'sideslip angle[deg]  = ', x(2) * 180 / pi 
            write(t%iunit_trim, '(A,*(1X,G25.17))') 'p[deg/sec]  = ', t%state(4) * 180.0 / pi 
            write(t%iunit_trim, '(A,*(1X,G25.17))') 'q[deg/sec]  = ', t%state(5) * 180.0 / pi 
            write(t%iunit_trim, '(A,*(1X,G25.17))') 'r[deg/sec]  = ', t%state(6) * 180.0 / pi 
            write(t%iunit_trim, '(A,*(1X,G25.17))') 'aileron deflection[deg]  = ', x(3) * 180.0 / pi 
            write(t%iunit_trim, '(A,*(1X,G25.17))') 'elevator deflection[deg] = ', x(4) * 180.0 / pi 
            write(t%iunit_trim, '(A,*(1X,G25.17))') 'rudder deflection[deg]   = ', x(5) * 180.0 / pi 
            write(t%iunit_trim, '(A,*(1X,G25.17))') 'throttle[none]           = ', x(6)

          end if           
        
        end subroutine

      !----------------------------------------
      ! Calculate Residual
        function calc_r(t, x, n_free) result(ans)
          implicit none
          type(vehicle_t) :: t
          real, intent(in) :: x(9)
          integer, intent(in) :: n_free
          integer :: last
          real :: FM(6) 
          real :: g, ac, xyzdot(3)
          real :: ca, cb, sa, sb, cp, sp, ct, st
          real :: u, v, w, euler(3)
          real :: angular_rates(3)
          real :: ans(n_free), temp_state(13), dummy_res(13)

          ! Pull out controls
          t%controls(1:4) = x(3:6)
          
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
          temp_state(1:3)   = t%init_airspeed * (/ca*cb, sb, sa*cb/) 
          temp_state(4:6)   = angular_rates(:)
          temp_state(9)     = -t%init_alt
          temp_state(7:8)   = 0.0
          temp_state(10:13) = euler_to_quat(euler)
          
          ! Set controls 
          t%controls = x(3:6) 

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

          if (t%trim%verbose) then 
            write(t%iunit_trim, '(A,*(1X,G25.17))') '       x =', x
            if(t%trim%type == 'sct') then 
              write(t%iunit_trim, *) '         p[deg/s] = ', angular_rates(1) * 180.0 / pi 
              write(t%iunit_trim, *) '         q[deg/s] = ', angular_rates(2) * 180.0 / pi 
              write(t%iunit_trim, *) '         r[deg/s] = ', angular_rates(3) * 180.0 / pi 
            end if 
              write(t%iunit_trim, '(A,*(1X,G25.17))') '       R =', ans
          end if         

        end function calc_r

        function calc_relative_climb_angle(y) result(ans)
          implicit none 
          real :: y(13), ans 
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
          real :: x1(13), x(13)

          x = t%state

          ! Step vehicle forward in time
          x1 = rk4(t, time, t%state, dt) 

          ! Call geographic model 
          if(geographic_model_ID > 0) call spherical_earth(t, x, x1)

          ! Normalize quaternion
          call quat_norm(x1(10:13)) 

          ! Update states
          t%state = x1 

        end subroutine

      !----------------------------------------   
      ! Spherical Earth model 
        subroutine spherical_earth(t, y1, y2)
          implicit none 
          type(vehicle_t) :: t 
          real, intent(in) :: y1(13), y2(13) 
          real :: d, dxf, dyf, dzf 
          real :: phi1, psi1 
          real :: dpsi_g, psi_g1 
          real :: theta, H1
          real :: xhat, yhat, zhat, xhat_prime, yhat_prime, zhat_prime
          real :: rhat, Chat, Shat
          real :: cphi1, sphi1, ct, st, cg1, sg1
          real :: cg, sg, quat(4)

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

            theta = d/ (earth_radius_ft + H1 - dzf*0.5) 
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
          end if 

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
        end subroutine 

      !----------------------------------------         
    !==================================================
    ! INTEGRATOR AND DIFFERENTIAL EQUATIONS
    !==================================================
      !----------------------------------------
      ! RK4 Integrator
        function rk4(t, t0, y1, delta_t) result(state)
          implicit none
          type(vehicle_t) :: t
          real, intent(in) :: t0, delta_t, y1(13)
          real, dimension(13) :: state, k1, k2, k3, k4

          if(t%rk4_verbose) then 
            write(t%iunit_rk4,*) '  state of the vehicle at the beginning of this RK4 integration step:'
            write(t%iunit_rk4,*) "time[s]             u[ft/s]              &
              &v[ft/s]              w[ft/s]              p[rad/s]             &
              &q[rad/s]             r[rad/s]             xf[ft]               &
              &yf[ft]               zf[ft]               e0                   &
              &ex                   ey                   ez"

            write(t%iunit_rk4,'(X,14(ES19.12,1X))') t0, t%state

            write(t%iunit_rk4,*) '  rk4 function called...'
          end if 

          ! K terms for RK4
          if (t%rk4_verbose) write(t%iunit_rk4,'(A,X,ES19.12,1X)')'    |           RK4 call number =             1'
          k1 = diff_eq(t, t0, y1)
          if (t%rk4_verbose) write(t%iunit_rk4,'(A,X,ES19.12,1X)')'    |           RK4 call number =             2'
          k2 = diff_eq(t, t0 + delta_t*0.5, y1 + k1 * delta_t*0.5)
          if (t%rk4_verbose) write(t%iunit_rk4,'(A,X,ES19.12,1X)')'    |           RK4 call number =             3'
          k3 = diff_eq(t, t0 + delta_t*0.5, y1 + k2 * delta_t*0.5)
          if (t%rk4_verbose) write(t%iunit_rk4,'(A,X,ES19.12,1X)')'    |           RK4 call number =             4'
          k4 = diff_eq(t, t0 + delta_t, y1 + k3 * delta_t)

          ! RK4 result
          state = y1 + (delta_t/6) * (k1 + 2*k2 + 2*k3 + k4)

          if (t%rk4_verbose) then 
            write(t%iunit_rk4,*) '  state of the vehicle after running RK4:'
            write(t%iunit_rk4,*) 'time[s]             u[ft/s]              &
              &v[ft/s]              w[ft/s]              p[rad/s]             &
              &q[rad/s]             r[rad/s]             xf[ft]               &
              &yf[ft]               zf[ft]               e0                   &
              &ex                   ey                   ez'
            write(t%iunit_rk4,'(X,14(ES19.12,1X))') t0+delta_t, state
            write(t%iunit_rk4,*) ' --------------------------- End of single RK4 integration step. ---------------------------'
          end if 

        end function rk4
      
      !----------------------------------------
      ! Differential Equations
        function diff_eq(t, time, state) result(dstate_dt)
          implicit none 
          type(vehicle_t), intent(inout) :: t
          real :: time
          real, target :: state(13) 
          real :: FM(6) 
          real :: dstate_dt(13) 
          real :: acceleration(3), angular_accelerations(3), rhs(3), velocity(3), quat_change(4) 
          real :: quat_inv(4) 
          real :: orientation_effect(3), angular_v_effect(3), gyroscopic_change(3)
          real :: inertia_effects(3), wind_velocity(3), gyroscopic_effects(3,3) 
          real :: quat_matrix(4,3) 
          real :: hxb_dot, hyb_dot, hzb_dot
          real :: Vxw, Vyw, Vzw, gravity_ft_per_sec2, ac
          real :: angular_inertia(3,3), angular_inertia_inv(3,3)
          real :: Ixx, Iyy, Izz, Ixy, Ixz, Iyz 
          real :: hxb, hyb, hzb        
          real, pointer :: u, v, w, p, q, r, e0, ex, ey, ez

          if (t%rk4_verbose) then 
            write(t%iunit_rk4,'(A,X,ES19.12,1X)') '    |                  time [s] = ', time 
            write(t%iunit_rk4,'(A,X,13(ES19.12,1X))') '    |    state vector coming in = ', state 
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
          FM = pseudo_aero(t, state)
          gravity_ft_per_sec2 = gravity_English(-state(9))

          ! Set gyroscopic effects
          hxb = t%h(1)
          hyb = t%h(2)
          hzb = t%h(3)

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

          gyroscopic_effects(1,1) =  0.0
          gyroscopic_effects(1,2) = -hzb
          gyroscopic_effects(1,3) =  hyb

          gyroscopic_effects(2,1) =  hzb
          gyroscopic_effects(2,2) =  0.0
          gyroscopic_effects(2,3) = -hxb

          gyroscopic_effects(3,1) = -hyb
          gyroscopic_effects(3,2) =  hxb
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
          rhs(1) = FM(4) + gyroscopic_effects(1,1)*p + gyroscopic_effects(1,2)*q + gyroscopic_effects(1,3)*r + inertia_effects(1) - gyroscopic_change(1)
          rhs(2) = FM(5) + gyroscopic_effects(2,1)*p + gyroscopic_effects(2,2)*q + gyroscopic_effects(2,3)*r + inertia_effects(2) - gyroscopic_change(2)
          rhs(3) = FM(6) + gyroscopic_effects(3,1)*p + gyroscopic_effects(3,2)*q + gyroscopic_effects(3,3)*r + inertia_effects(3) - gyroscopic_change(3)

          angular_accelerations(1) = angular_inertia_inv(1,1)*rhs(1) + angular_inertia_inv(1,2)*rhs(2) + angular_inertia_inv(1,3)*rhs(3)
          angular_accelerations(2) = angular_inertia_inv(2,1)*rhs(1) + angular_inertia_inv(2,2)*rhs(2) + angular_inertia_inv(2,3)*rhs(3)
          angular_accelerations(3) = angular_inertia_inv(3,1)*rhs(1) + angular_inertia_inv(3,2)*rhs(2) + angular_inertia_inv(3,3)*rhs(3)

          ! Velocity in inertial
          velocity = quat_base_to_dependent(state(1:3), state(10:13)) + wind_velocity

          ! Gravity relief
          ac = (velocity(1)**2 + velocity(2)**2) / (earth_radius_ft - state(9))

          ! Aacceleration in body
          acceleration(1) = FM(1)/t%mass + (gravity_ft_per_sec2 - ac)*orientation_effect(1) + angular_v_effect(1)
          acceleration(2) = FM(2)/t%mass + (gravity_ft_per_sec2 - ac)*orientation_effect(2) + angular_v_effect(2)
          acceleration(3) = FM(3)/t%mass + (gravity_ft_per_sec2 - ac)*orientation_effect(3) + angular_v_effect(3)

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

          if (t%rk4_verbose) then 
            write(t%iunit_rk4,'(A,X,6(ES19.12,1X))') '    | pseudo aerodynamics (F,M) = ', FM
            write(t%iunit_rk4,'(A,X,13(ES19.12,1X))') '    |           diff_eq results = ', dstate_dt
          end if 
        end function diff_eq

      !----------------------------------------
    !==================================================
    ! AERODYNAMICS AND FORCES
    !==================================================
      !----------------------------------------
      ! Aerodynamic Forces and Moments for f16
        function pseudo_aero(t, state) result(FM)
          implicit none
          type(vehicle_t) :: t
          real, intent(in) :: state(13)
          real :: FM(6) 
          real :: Re, geometric_altitude_ft, geopotential_altitude_ft
          real :: temp_R, pressure_lbf_per_ft2, density_slugs_per_ft3
          real :: dyn_viscosity_slug_per_ft_sec, sos_ft_per_sec
          real :: V, dyn_pressure
          real :: alpha, beta, beta_f, pbar, qbar, rbar, angular_rates(3)
          real :: CL, CL1, CD, CS, L, D, S, Cl_roll, Cm, Cn
          real :: CM1, CM2, mach_num
          real :: ca, cb, sa, sb
          real :: alpha_hat, beta_hat
          real :: delta_a, delta_e, delta_r
          real :: thrust, throttle
          real :: CL_s, CD_s, Cm_s 
          real :: sigma_D, sigma_L, sigma_m, sign_a
          
          ! Ccmpressibility
          real :: blend, drag_factor, lift_factor, moment_factor
          real :: pg_factor, kt_factor
          real :: m_high, m_low, m_trans_start
          real :: sqrt_term
          real :: max_CL_factor, max_CD_factor, max_Cl_roll_factor
          
          ! Build atmosphere
          geometric_altitude_ft = -state(9)
          call std_atm_English(&
            geometric_altitude_ft, geopotential_altitude_ft,     & 
            temp_R, pressure_lbf_per_ft2, density_slugs_per_ft3, & 
            dyn_viscosity_slug_per_ft_sec, sos_ft_per_sec)

          ! Calculate velocity unit vector
          V            =  (state(1)**2 + state(2)**2 + state(3)**2)**0.5
          dyn_pressure = 0.5 * density_slugs_per_ft3 * V **2 * t%planform_area

          ! Calculate alpha and beta 3.4.4 and 3.4.5
          alpha  =  atan2(state(3) , state(1))
          beta   =   asin(state(2) / V)
          beta_f = atan2(state(2) , state(1))

          ca = cos(alpha)
          cb = cos(beta)
          sa = sin(alpha)
          sb = sin(beta)

          ! Define alpha_hat and beta_hat
          alpha_hat = 0.0
          beta_hat  = 0.0

          ! Calculate rotation rates from eq 3.4.23
          angular_rates(1) = 1 / (2*V) * state(4) * t%lateral_length
          angular_rates(2) = 1 / (2*V) * state(5) * t%longitudinal_length
          angular_rates(3) = 1 / (2*V) * state(6) * t%lateral_length

          pbar = angular_rates(1)
          qbar = angular_rates(2)
          rbar = angular_rates(3)

          if(t%type == 'aircraft') then 
            ! Pull out controls
            delta_a  =  t%controls(1)
            delta_e  =  t%controls(2)
            delta_r  =  t%controls(3)
            throttle = t%controls(4)

            ! Calculate lift, drag, side force coef
            CL1 =  t%CL0 + t%CL_alpha * alpha
            CL  = CL1 + t%CL_qbar * qbar + t%CL_alphahat * alpha_hat + t%CL_elevator * delta_e
            CS  = t%CS_beta * beta + (t%CS_pbar + t%CS_alpha_pbar * alpha) * pbar + t%CS_rbar * rbar &
                  + t%CS_aileron * delta_a + t%CS_rudder * delta_r
            CD  =  t%CD_L0 + t%CD_L1 * CL1 + t%CD_L1_L1 * CL1 **2 + t%CD_CS_CS * CS **2 &
                  + (t%CD_qbar + t%CD_alpha_qbar * alpha) * qbar + (t%CD_elevator + t%CD_alpha_elevator * alpha) &
                  * delta_e + t%CD_elevator_elevator * delta_e ** 2

            ! Calculate roll, pitch, yaw coef
            Cl_roll = t%Cl_beta * beta + t%Cl_pbar * pbar + (t%Cl_rbar + t%Cl_alpha_rbar * alpha) * rbar &
                      + t%Cl_aileron * delta_a + t%Cl_rudder * delta_r  ! roll moment
            Cm      = t%Cm_0 + t%Cm_alpha * alpha + t%Cm_qbar * qbar + t%Cm_alphahat * alpha_hat + t%Cm_elevator * delta_e ! pitch moment
            Cn      = t%Cn_beta * beta + (t%Cn_pbar + t%Cn_alpha_pbar * alpha) * pbar + t%Cn_rbar * rbar &
                      + (t%Cn_aileron + t%Cn_alpha_aileron * alpha) * delta_a + t%Cn_rudder * delta_r ! yaw moment
          
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

            Cm_s    = t%CM_min * sa**2 * sign_a
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

          L =   CL * dyn_pressure
          S =   CS * dyn_pressure
          D =   CD * dyn_pressure

          ! Table 3.4.4
          FM(1) = - (ca*(D*cb + S*sb) - L*sa)
          FM(2) = (S*cb - D*sb)
          FM(3) = - (sa*(D*cb + S*sb) + L*ca)

          FM(4) = Cl_roll * dyn_pressure * t%lateral_length
          FM(5) = Cm       * dyn_pressure * t%longitudinal_length
          FM(6) = Cn       * dyn_pressure * t%lateral_length

          ! Add thrust
          thrust = throttle * t%T0 * (density_slugs_per_ft3/t%rho0) ** t%Ta
          FM(1)  = FM(1) + thrust

          ! Shift CG location 
          FM(4:6) = FM(4:6) + cross_product(t%aero_ref_location, FM(1:3))

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

    end module vehicle_m