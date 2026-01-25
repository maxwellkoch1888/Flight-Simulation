module vehicle_m
    use koch_m
    use jsonx_m
    use micro_time_m
    use linalg_mod

    implicit none

    ! Trim solver type
    type trim_solver_t
      real :: step_size, relaxation_factor, tolerance, max_iterations
      logical :: verbose 
    end type trim_solver_t

    ! Trim solver type
    type trim_settings_t
      character(len=:), allocatable :: type 

      logical, allocatable :: free_vars(:) 
      real :: climb_angle, sideslip_angle
      real :: load_factor

      type(trim_solver_t) :: solver 
    end type trim_settings_t

    ! Vehicle type
    type vehicle_t
      type(json_value), pointer :: j_vehicle
          
      character(len=:), allocatable :: name
      character(len=:), allocatable :: type
      character(100) :: states_filename, rk4_filename, trim_filename

      logical :: run_physics
      logical :: save_states
      integer :: iunit_states, iunit_rk4, iunit_trim

      ! Location 
      real :: lat, long 

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

    ! BUILD GLOBAL VARIABLES FOR THE MODULE
    real :: FM(6)
    real :: rho0

    ! ! ADD VARIABLES FOR TRIM ALGORITHM
    ! character(:), allocatable :: sim_type
    ! character(:), allocatable :: trim_type
    ! real :: relaxation_factor, tolerance, max_iterations, finite_difference_step_size
    ! real :: bank_angle0, sideslip_angle0, climb_angle0, elevation_angle0
    ! logical :: trim_verbose, exam_answers

  ! 
  contains
  ! INITIALIZATION FUNCTIONS
    !=========================
    ! VEHICLE INITIALIZATION
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

          write(*,*) '- mass' 
          call jsonx_get(t%j_vehicle, 'mass.weight[lbf]', t%mass)
          t%mass = t%mass/gravity_English(0.0)

          ! READ MASS PROPERTIES
          call jsonx_get(t%j_vehicle, 'mass.Ixx[slug-ft^2]',  t%inertia(1,1))
          call jsonx_get(t%j_vehicle, 'mass.Iyy[slug-ft^2]',  t%inertia(2,2))
          call jsonx_get(t%j_vehicle, 'mass.Izz[slug-ft^2]',  t%inertia(3,3))
          call jsonx_get(t%j_vehicle, 'mass.Ixy[slug-ft^2]',  t%inertia(1,2), 0.0)
          call jsonx_get(t%j_vehicle, 'mass.Ixz[slug-ft^2]',  t%inertia(1,3), 0.0)
          call jsonx_get(t%j_vehicle, 'mass.Iyz[slug-ft^2]',  t%inertia(2,3), 0.0)
          call jsonx_get(t%j_vehicle, 'mass.h[slug-ft^2/s]',  t%h, 0.0, 3)

          ! READ IN ALL AERODYNAMIC DATA      
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

          ! STALL CONDITIONS
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

          ! DEFINE THRUST COEFFICIENTS FOR THRUST MODEL
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
          call jsonx_get(t%j_vehicle, 'initial.latitude[deg]',  t%lat)
          call jsonx_get(t%j_vehicle, 'initial.longitude[deg]', t%long)
          t%lat  = t%lat  * pi / 180.0
          t%long = t%long * pi / 180.0

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

    !=========================
    ! WRITE STATES TO A FILE
      subroutine vehicle_write_state(t, time)
        implicit none 
        type(vehicle_t) :: t 
        real, intent(in) :: time 

        ! write(t%iunit_states,*) time, t%state, t%controls 
        write(t%iunit_states,'(14(ES20.13,1X))') time, t%state

      end subroutine 

    !=========================
    ! STATE INITIAL CONDITION
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
    !=========================
    ! TRIM INITIAL
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
        real :: sa, sb, ca, cb

        allocate(t%trim%free_vars(n_vars))

        write(*,*) '  -trimming'
        ! Get json objects 
        call jsonx_get(t%j_vehicle, 'initial', j_initial) 
        call jsonx_get(j_initial, 'trim',      j_trim)
        call jsonx_get(j_trim, 'solver',       j_solver)
        call jsonx_get(j_trim, 'type',         t%trim%type) 

        write(*,*)
        write(*,*) '    Trimming vehicle for ', t%trim%type 

        call jsonx_get(j_solver, 'finite_difference_step_size', t%trim%solver%step_size)
        call jsonx_get(j_solver, 'relaxation_factor',           t%trim%solver%relaxation_factor)
        call jsonx_get(j_solver, 'tolerance',                   t%trim%solver%tolerance)
        call jsonx_get(j_solver, 'max_iterations',              t%trim%solver%max_iterations)
        call jsonx_get(j_solver, 'verbose',                     t%trim%solver%verbose)

        if(t%trim%solver%verbose) then 
          t%trim_filename = 'output_files/' // trim(t%name)//'_trim.txt'
          open(newunit=t%iunit_trim, file=t%trim_filename, status='REPLACE')
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
        x(7:9) = t%init_eul 
        t%trim%free_vars(:) =   .true. 
        t%trim%free_vars(7:9) = .false. 
        t%trim%climb_angle =    0.0
        t%trim%load_factor =    1.0

        n_free = count(t%trim%free_vars)

        allocate(res(n_free))
        allocate(R_pos(n_free))
        allocate(R_neg(n_free))
        allocate(jac(n_free, n_free))

        idx_free = pack([(i,i=1,n_vars)], t%trim%free_vars) ! yields an array of indices of all free variables

        if(t%trim%solver%verbose) then 
          write(t%iunit_trim, *)
          write(t%iunit_trim, *) 'n_vars = ', n_vars 
          write(t%iunit_trim, *) 'n_free = ', n_free 
          write(t%iunit_trim, *) 'free_vars = ', t%trim%free_vars
          write(t%iunit_trim, *) 'idx_free = ', idx_free(:)
        end if 
        write(*,*) 'n_vars = ', n_vars 
        write(*,*) 'n_free = ', n_free 
        write(*,*) 'free_vars = ', t%trim%free_vars
        write(*,*) 'idx_free = ', idx_free(:)
        write(*,*) 
        write(*,*) 'iteration Res alphs deg beta deg' 

        iter = 0

        if(t%trim%solver%verbose) then 
          write(t%iunit_trim,*)
          write(t%iunit_trim,*) 'Initial Guess: '
        end if 

        ! Calculate residual and initial error
        res = calc_r(t, x) 
        error = maxval(abs(res))

        if(t%trim%solver%verbose) then 
          write(t%iunit_trim,*)
          write(t%iunit_trim,*) 'Beginning Solution Process...'
        end if 

        ! BEGIN NEWTONS METHOD
        ! Use central difference method to find jacobian
        do while(error > t%trim%solver%tolerance)
          iter = iter + 1

          if(t%trim%solver%verbose) then 
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
            x(k) = x(k) + t%trim%solver%step_size
            R_pos = calc_r(t, x)
            
            if(t%trim%solver%verbose) then 
              write(t%iunit_trim, '(A,I0,A)') 'Computing gradient relative to x[', k, ']'
              write(t%iunit_trim, '(A)') '   Positive Finite-Difference Step '
              write(t%iunit_trim, '(A,6(1X,ES20.12))') '      x =', (x(j), j=1,6)
              write(t%iunit_trim, '(A,6(1X,ES20.12))') '      R =', (R_pos(j), j=1,6)
            end if 
            
            x(k) = x(k) - 2.0 * t%trim%solver%step_size
            R_neg = calc_r(t, x)
            
            if(t%trim%solver%verbose) then 
              write(t%iunit_trim, '(A)') '   Negative Finite-Difference Step'
              write(t%iunit_trim, '(A,6(1X,ES20.12))') '      x =', (x(j), j=1,6)
              write(t%iunit_trim, '(A,6(1X,ES20.12))') '      R =', (R_neg(j), j=1,6)
              write(t%iunit_trim,'(A)') ''
            end if

            x(k) = x(k) + t%trim%solver%step_size
            jac(:,i) = (R_pos - R_neg) / 2.0 / t%trim%solver%step_size
          end do 

          if (t%trim%solver%verbose) then
            write(t%iunit_trim, '(A)') 'Jacobian Matrix ='
            do i = 1, size(jac,1)
                write(t%iunit_trim,'(*(1X,G25.17))') (jac(i,j), j=1,size(jac,2))
            end do
          end if

          res = calc_r(t, x)

          ! Calculate delta x and add relaxation factor
          call lu_solve(6, jac, -res, delta_x)
          x = x + t%trim%solver%relaxation_factor * delta_x

          res = calc_r(t,x) 
          error = maxval(abs(res))

          if (t%trim%solver%verbose) then
            write(t%iunit_trim,*)
            write(t%iunit_trim, '(A,*(1X,G25.17))') ' Delta X =', (delta_x(k), k=1,6)
            write(t%iunit_trim, '(A,*(1X,G25.17))') '   New X = ', x(:)
            write(t%iunit_trim, '(A,*(1X,G25.17))') 'Residual = ', res(:)
            write(t%iunit_trim, '(A,*(1X,G25.17))') '   Error = ', error 
          end if 

        end do 

        sa = sin(x(1))
        ca = cos(x(1))
        sb = sin(x(2))
        cb = cos(x(2)) 
      
      end subroutine
    !=========================
    ! TICK A VEHICLE FORWARD IN TIME 
      subroutine vehicle_tick_state(t, time, dt)
        implicit none 
        type(vehicle_t) :: t 
        real, intent(in) :: time, dt 
        real :: x(13) 

        ! STEP VEHICLE FORWARD IN TIME 
        x = rk4(t, time, t%state, dt) 

        ! NORMALIZE THE QUATERNION 
        call quat_norm(x(10:13)) 

        ! CHECK FOR HIGH ROTATION RATES 
        ! if (sqrt(y1(4)**2 + y1(5)**2 + y1(6)**2)/2.0/pi/dt > 0.1 ) write(*,*) 'Warning: High vehicle rotation rate relative to timestep.'

        ! UPDATE STATES
        t%state = x 

      end subroutine
  ! 
  ! INTEGRATOR AND EQN OF MOTION
    !=========================
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

        ! DEFINE THE K TERMS FOR RK4 METHOD
        if (t%rk4_verbose) write(t%iunit_rk4,'(A,X,ES19.12,1X)')'    |           RK4 call number =             1'
        k1 = diff_eq(t, t0, y1)
        if (t%rk4_verbose) write(t%iunit_rk4,'(A,X,ES19.12,1X)')'    |           RK4 call number =             2'
        k2 = diff_eq(t, t0 + delta_t*0.5, y1 + k1 * delta_t*0.5)
        if (t%rk4_verbose) write(t%iunit_rk4,'(A,X,ES19.12,1X)')'    |           RK4 call number =             3'
        k3 = diff_eq(t, t0 + delta_t*0.5, y1 + k2 * delta_t*0.5)
        if (t%rk4_verbose) write(t%iunit_rk4,'(A,X,ES19.12,1X)')'    |           RK4 call number =             4'
        k4 = diff_eq(t, t0 + delta_t, y1 + k3 * delta_t)

        ! DEFINE THE RESULT FROM RK4
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
    
    !=========================
    ! Equations of Motion: (/u,v,w, p,q,r, x,y,z, e0,ex,ey,ez/)
      function diff_eq(t, time, state) result(dstate_dt)
        implicit none 
        type(vehicle_t), intent(inout) :: t
        real :: time
        real, target :: state(13) 
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
        
        ! UNPACK STATES
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

        ! UNPACK INERTIA
        Ixx = t%inertia(1,1)
        Iyy = t%inertia(2,2)
        Izz = t%inertia(3,3)
        Ixy = t%inertia(1,2)
        Ixz = t%inertia(1,3)
        Iyz = t%inertia(2,3)
      
        ! CALCULATE FORCES AND MOMENTS
        call pseudo_aero(t, state)
        gravity_ft_per_sec2 = gravity_English(-state(9))

        ! SET GYROSCOPIC EFFECTS
        hxb = t%h(1)
        hyb = t%h(2)
        hzb = t%h(3)

        ! SET GYROSCOPIC CHANGE AND WIND VELOCITY TO ZERO
        hxb_dot = 0.0
        hyb_dot = 0.0
        hzb_dot = 0.0
        Vxw = 0.0
        Vyw = 0.0
        Vzw = 0.0

        ! BUILD MATRICES/ VECTORS USED IN DIFFERENTIAL EQUATION
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

        ! BUILD THE DIFFERENTIAL EQUATIONS
        ! ROLL, PITCH, YAW ACCELERATIONS
        rhs(1) = FM(4) + gyroscopic_effects(1,1)*p + gyroscopic_effects(1,2)*q + gyroscopic_effects(1,3)*r + inertia_effects(1) - gyroscopic_change(1)
        rhs(2) = FM(5) + gyroscopic_effects(2,1)*p + gyroscopic_effects(2,2)*q + gyroscopic_effects(2,3)*r + inertia_effects(2) - gyroscopic_change(2)
        rhs(3) = FM(6) + gyroscopic_effects(3,1)*p + gyroscopic_effects(3,2)*q + gyroscopic_effects(3,3)*r + inertia_effects(3) - gyroscopic_change(3)

        angular_accelerations(1) = angular_inertia_inv(1,1)*rhs(1) + angular_inertia_inv(1,2)*rhs(2) + angular_inertia_inv(1,3)*rhs(3)
        angular_accelerations(2) = angular_inertia_inv(2,1)*rhs(1) + angular_inertia_inv(2,2)*rhs(2) + angular_inertia_inv(2,3)*rhs(3)
        angular_accelerations(3) = angular_inertia_inv(3,1)*rhs(1) + angular_inertia_inv(3,2)*rhs(2) + angular_inertia_inv(3,3)*rhs(3)

        ! VELOCITY IN THE INERITAL FRAME
        velocity = quat_base_to_dependent(state(1:3), state(10:13)) + wind_velocity

        ! ADD GRAVITY RELIEF FOR EARTH'S CURVATURE
        ac = (velocity(1)**2 + velocity(2)**2) / (6366707.01949371/0.3048 - state(9))

        ! ACCELERATION IN BODY FRAME
        acceleration(1) = FM(1)/t%mass + (gravity_ft_per_sec2 - ac)*orientation_effect(1) + angular_v_effect(1)
        acceleration(2) = FM(2)/t%mass + (gravity_ft_per_sec2 - ac)*orientation_effect(2) + angular_v_effect(2)
        acceleration(3) = FM(3)/t%mass + (gravity_ft_per_sec2 - ac)*orientation_effect(3) + angular_v_effect(3)

        ! AIRCRAFT ORIENTATION RATE OF CHANGE
        quat_change(1) = 0.5 * (quat_matrix(1,1)*p + quat_matrix(1,2)*q + quat_matrix(1,3)*r)
        quat_change(2) = 0.5 * (quat_matrix(2,1)*p + quat_matrix(2,2)*q + quat_matrix(2,3)*r)
        quat_change(3) = 0.5 * (quat_matrix(3,1)*p + quat_matrix(3,2)*q + quat_matrix(3,3)*r)
        quat_change(4) = 0.5 * (quat_matrix(4,1)*p + quat_matrix(4,2)*q + quat_matrix(4,3)*r)

        ! RETURN THE STATE DERIVATIVE
        dstate_dt(1:3)  = acceleration
        dstate_dt(4:6)  = angular_accelerations
        dstate_dt(7:9)  = velocity
        dstate_dt(10:13) = quat_change

        if (t%rk4_verbose) then 
          write(t%iunit_rk4,'(A,X,6(ES19.12,1X))') '    | pseudo aerodynamics (F,M) = ', FM
          write(t%iunit_rk4,'(A,X,13(ES19.12,1X))') '    |           diff_eq results = ', dstate_dt
        end if 
      end function diff_eq
  ! 
  ! AERODYNAMICS AND FORCES
    !=========================
    ! Aerodynamic Forces and Moments for f16
      subroutine pseudo_aero(t, state)
        implicit none
        type(vehicle_t) :: t
        real, intent(in) :: state(13)
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
        
        ! COMPRESSIBILITY
        real :: blend, drag_factor, lift_factor, moment_factor
        real :: pg_factor, kt_factor
        real :: m_high, m_low, m_trans_start
        real :: sqrt_term
        real :: max_CL_factor, max_CD_factor, max_Cl_roll_factor
        
        ! BUILD THE ATMOSPHERE 
        geometric_altitude_ft = -state(9)
        call std_atm_English(&
          geometric_altitude_ft, geopotential_altitude_ft,     & 
          temp_R, pressure_lbf_per_ft2, density_slugs_per_ft3, & 
          dyn_viscosity_slug_per_ft_sec, sos_ft_per_sec)

        ! CALCULATE VELOCITY UNIT VECTOR
        V =  (state(1)**2 + state(2)**2 + state(3)**2)**0.5
        dyn_pressure = 0.5 * density_slugs_per_ft3 * V **2 * t%planform_area

        ! CALCULATE ALPHA AND BETA 3.4.4 and 3.4.5
        alpha =  atan2(state(3) , state(1))
        beta =   asin(state(2) / V)
        beta_f = atan2(state(2) , state(1))

        ca = cos(alpha)
        cb = cos(beta)
        sa = sin(alpha)
        sb = sin(beta)

        ! CALCULATE ALPHA_HAT USING EQN 3.4.20
        alpha_hat = 0.0
        beta_hat = 0.0

        ! CALCULATE PBAR, QBAR, AND RBAR from eq 3.4.23
        angular_rates(1) = 1 / (2*V) * state(4) * t%lateral_length
        angular_rates(2) = 1 / (2*V) * state(5) * t%longitudinal_length
        angular_rates(3) = 1 / (2*V) * state(6) * t%lateral_length

        pbar = angular_rates(1)
        qbar = angular_rates(2)
        rbar = angular_rates(3)

        if(t%type == 'aircraft') then 
          ! PULL OUT CONTROLS
          delta_a =  t%controls(1)
          delta_e =  t%controls(2)
          delta_r =  t%controls(3)
          throttle = t%controls(4)

          ! CALCULATE THE LIFT, DRAG, AND SIDE FORCE COEFFICIENTS
          CL1 =  t%CL0 + t%CL_alpha * alpha
          CL  = CL1 + t%CL_qbar * qbar + t%CL_alphahat * alpha_hat + t%CL_elevator * delta_e
          CS  = t%CS_beta * beta + (t%CS_pbar + t%CS_alpha_pbar * alpha) * pbar + t%CS_rbar * rbar &
              + t%CS_aileron * delta_a + t%CS_rudder * delta_r
          CD  =  t%CD_L0 + t%CD_L1 * CL1 + t%CD_L1_L1 * CL1 **2 + t%CD_CS_CS * CS **2 &
                + (t%CD_qbar + t%CD_alpha_qbar * alpha) * qbar + (t%CD_elevator + t%CD_alpha_elevator * alpha) &
                * delta_e + t%CD_elevator_elevator * delta_e ** 2

          ! CALCULATE THE ROLL, PITCH, AND YAW COEFFICIENTS
          Cl_roll = t%Cl_beta * beta + t%Cl_pbar * pbar + (t%Cl_rbar + t%Cl_alpha_rbar * alpha) * rbar &
                     + t%Cl_aileron * delta_a + t%Cl_rudder * delta_r  ! roll moment
          Cm    = t%Cm_0 + t%Cm_alpha * alpha + t%Cm_qbar * qbar + t%Cm_alphahat * alpha_hat + t%Cm_elevator * delta_e ! pitch moment
          Cn       = t%Cn_beta * beta + (t%Cn_pbar + t%Cn_alpha_pbar * alpha) * pbar + t%Cn_rbar * rbar &
                     + (t%Cn_aileron + t%Cn_alpha_aileron * alpha) * delta_a + t%Cn_rudder * delta_r ! yaw moment
        
        else if (t%type == 'arrow') then 
          CL = t%CL_alpha * alpha 
          CS = -t%CL_alpha * beta_f 
          CD = t%CD_L0 + t%CD_L1_L1*CL**2 + t%CD_L1_L1*CS**2 
          
          Cl_roll  = t%Cl_l0 + t%Cl_pbar*pbar
          Cm       = t%Cm_alpha*alpha + t%Cm_qbar*qbar 
          Cn       = -t%Cm_alpha*beta_f + t%Cm_qbar*rbar

        else if (t%type == 'sphere') then 
          CL       = 0.0 
          CS       = 0.0 
          Cl_roll = 0.0 
          Cm       = 0.0 
          Cn       = 0.0
          
          ! CALCULATE THE REYNOLDS NUMBER
          Re = density_slugs_per_ft3 * V * 2 * t%longitudinal_length / dyn_viscosity_slug_per_ft_sec

          ! COMPUTE THE DRAG COEFFICIENT
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

        ! STALL MODEL FOR FORCES
        if (t%stall) then 
          sign_a = sign(1.0,alpha)
          CL_s = 2 * sign_a * sa**2 * ca 
          CD_s = 2 * (sin(abs(alpha)))**3
          sigma_L = calc_sigma(t%CL_lambda_b, t%CL_alpha_0, t%CL_alpha_s, alpha)
          sigma_D = calc_sigma(t%CD_lambda_b, t%CD_alpha_0, t%CD_alpha_s, alpha)
          
          CL = CL * (1 - sigma_L) + CL_s * sigma_L 
          CD = CD * (1 - sigma_D) + CD_s * sigma_D 

          Cm_s = t%CM_min * sa**2 * sign_a
          sigma_m = calc_sigma(t%Cm_lambda_b, t%Cm_alpha_0, t%Cm_alpha_s, alpha)
          Cm = Cm * (1 - sigma_m) + Cm_s * sigma_m 
     
        end if 

        ! COMPRESSIBILIITY MODEL, USING PRANDTL-GLAUERT CORRECTION
        if (t%compressibility) then 
          CM1 = 2.13/ (t%sweep + 0.15)**2
          CM2 = 15.35*t%sweep**2 - 19.64*t%sweep +16.86
          mach_num = V / sos_ft_per_sec

          ! MACH BREAKPOINTS
          m_low         = 0.60    ! only use prandtl-glauert
          m_high        = 0.92    ! only use karman-tsien
          m_trans_start = 0.88    ! start transsonic region

          ! Prandtl-Glauert factor
          sqrt_term = sqrt(max(1.0 - mach_num**2, tol))
          pg_factor = 1.0 / sqrt_term

          ! Karman-Tsien factor
          kt_factor = pg_factor * (1.0 + mach_num**2 / (1.0 + sqrt_term))

          ! BLEND BETWEEN PG AND KT
          if (mach_num <= m_low) then
            blend = 0.0
          else if (mach_num >= m_high) then
            blend = 1.0
          else
            blend = (mach_num - m_low) / (m_high - m_low)
          end if

          ! FINAL FACTOR USING BLEND
          lift_factor = (1.0 - blend) * pg_factor + blend * kt_factor
          moment_factor = lift_factor
          drag_factor = 1.0 + CM1 * mach_num**CM2

          ! APPLY THE FACTORS
          max_CL_factor = 6.5
          max_CD_factor = 4
          max_Cl_roll_factor = 2.5

          if (mach_num < 0.92) then ! accurate range
            CL = CL * lift_factor
            CS = CS * lift_factor
            Cl_roll = Cl_roll * moment_factor
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

          ! apply drag multiplier
          CD = CD * min(drag_factor, max_CD_factor)
        end if 

        L =   CL * dyn_pressure
        S =   CS * dyn_pressure
        D =   CD * dyn_pressure

        ! TABLE 3.4.4
        FM(1) = - (ca*(D*cb + S*sb) - L*sa)
        FM(2) = (S*cb - D*sb)
        FM(3) = - (sa*(D*cb + S*sb) + L*ca)

        FM(4) = Cl_roll * dyn_pressure * t%lateral_length
        FM(5) = Cm       * dyn_pressure * t%longitudinal_length
        FM(6) = Cn       * dyn_pressure * t%lateral_length

        ! ADD THE ENGINE THRUST
        thrust = throttle * t%T0 * (density_slugs_per_ft3/t%rho0) ** t%Ta
        FM(1) = FM(1) + thrust

        ! SHIFT CG LOCATION
        FM(4:6) = FM(4:6) + cross_product(t%aero_ref_location, FM(1:3))

      end subroutine pseudo_aero

    !=========================
    ! Check Compressiblity model
      ! subroutine check_compressibility(t, state)
      !   implicit none
      !   type(vehicle_t) :: t
      !   integer :: i, j
      !   real :: state(13)
      !   real :: alpha, beta, states(13)
      !   real :: N, Y, A
      !   real :: ca, cb, sa, sb
      !   real :: CL, CD, Cm
      !   real :: const
      !   real :: geometric_altitude_ft, geopotential_altitude_ft
      !   real :: temp_R, pressure_lbf_per_ft2, density_slugs_per_ft3
      !   real :: dyn_viscosity_slug_per_ft_sec, sos_ft_per_sec
      !   real :: airspeed, mach_num
      !   real, dimension(9) :: mach_num_list

      !   ! Mach numbers to sweep
      !   mach_num_list = (/ 0.3, 0.6, 0.7, 0.8, 0.9, 0.95, 1.0, 1.05, 1.1 /)

      !   ! Get altitude for density
      !   geometric_altitude_ft = -t%initial_state(9)
      !   call std_atm_English( &
      !     geometric_altitude_ft, geopotential_altitude_ft, &
      !     temp_R, pressure_lbf_per_ft2, density_slugs_per_ft3, &
      !     dyn_viscosity_slug_per_ft_sec, sos_ft_per_sec)

      !   t%controls = 0.0
      !   states   = 0.0

      !   ! Main Mach loop -----------------------------------------------
      !   do j = 1, size(mach_num_list)
      !     mach_num = mach_num_list(j)

      !     ! write(io_unit,*) '--------------------------------------------------'
      !     ! write(io_unit,*) 'Mach Number = ', mach_num
      !     ! write(io_unit,*) '--------------------------------------------------'

      !     airspeed = mach_num * sos_ft_per_sec

      !     const = 0.5 * density_slugs_per_ft3 * airspeed**2 * t%planform_area

      !     ! write(io_unit,*) '   alpha(deg)                CL                        CD                        Cm'

      !     do i = -180, 180
      !       alpha = real(i) * pi / 180.0
      !       beta  = 0.0

      !       ca = cos(alpha)
      !       cb = cos(beta)
      !       sa = sin(alpha)
      !       sb = sin(beta)

      !       states = 0.0
      !       states(1) = airspeed * ca * cb
      !       states(2) = airspeed * sb
      !       states(3) = airspeed * sa * cb
      !       states(9) = t%initial_state(9)

      !       call pseudo_aero(t, states)

      !       A = -FM(1)
      !       Y =  FM(2)
      !       N = -FM(3)

      !       ! Aerodynamic coefficients
      !       CL = (N * ca - A * sa) / const
      !       CD = (A * ca * cb - Y * sb + N * sa * cb) / const
      !       Cm = FM(5) / (const * t%longitudinal_length)

      !       ! write(io_unit,*) alpha*180.0/pi, CL, CD, Cm
      !     end do

      !   end do

      ! end subroutine
       
    !=========================
    ! Calculate Sigma
      function calc_sigma(lambda_b, alpha_0, alpha_s, alpha) result(sigma)
        implicit none
        real, intent(in) :: lambda_b, alpha_0, alpha_s, alpha
        real :: sigma, neg, pos

        pos = exp( lambda_b * (alpha - alpha_0 + alpha_s))
        neg = exp(-lambda_b * (alpha - alpha_0 - alpha_s))
        sigma = (1.0 + neg + pos) / ((1.0 + neg)*(1.0 + pos))

      end function calc_sigma

    !=========================
    ! Mass and Inertia
      subroutine mass_inertia(t)
        implicit none
        type(vehicle_t) :: t
        real :: weight
        real :: Ixx, Iyy, Izz, Ixy, Ixz, Iyz

        ! DEFINE IDENTITY MATRIX 
        t%inertia(1,1) = Ixx
        t%inertia(1,2) = Ixy
        t%inertia(1,3) = Ixz
        t%inertia(2,1) = Ixy
        t%inertia(2,2) = Iyy
        t%inertia(2,3) = Iyz
        t%inertia(3,1) = Ixz
        t%inertia(3,2) = Iyz
        t%inertia(3,3) = Izz

        ! CALCULATE MASS AND INERTIA
        t%mass = weight / 32.17404855643
        t%inertia_inv = matrix_inv(t%inertia)
      end subroutine mass_inertia
  ! 
  ! CALCULATE TRIM STATE
    !=========================
    ! Old Trim Algorithm
      ! function trim_algorithm(V_mag, height, euler, tolerance, trim_type) result(G)
      !   implicit none 
      !   real :: V_mag, height, tolerance
      !   real :: alpha, beta, p, q, r, delta_a, delta_e, delta_r, throttle 
      !   real :: elevation_angle, azimuth_angle
      !   real :: c_bank, c_elev, s_bank, s_elev, ca, cb, sa, sb, error, pw
      !   real :: u, v, w, velocities(3), gravity
      !   real :: G(6), res(6), iteration_residual(6)
      !   real :: angular_rates(3), euler(3), print_statement(13)
      !   real :: cgamma, sgamma, climb_angle, solution, theta1, theta2, solution1
      !   integer :: k, iteration, case_number
      !   character(*), intent(in) :: trim_type
        
      !   ! DETERMINE THE TRIM TYPE BEING SPECIFIED
      !   if (trim_type == 'sct') then 
      !     if (elevation_angle0 /= -999.0) then 
      !       case_number = 1 ! sct, elevation_angle, bank_angle
      !     else
      !       case_number = 2 ! sct, climb_angle, bank_angle
      !     end if 
      !   else 
      !     if (elevation_angle0 /= -999.0) then 
      !       if (t%sideslip_angle0 == -999.0) then 
      !         case_number = 3 ! shss, elevation_angle, bank_angle 
      !       else 
      !         case_number = 4 ! shss, elevation_angle, sideslip_angle
      !       end if 
      !     else 
      !       if (sideslip_angle0 == -999.0) then 
      !         case_number = 5 ! shss, climb_angle, bank_angle
      !       else 
      !         case_number = 6 ! shss, climb_angle, sideslip_angle
      !       end if 
      !     end if 
      !   end if 
        
      !   if (trim_verbose) then 
      !     write(io_unit,*) 'Case Number:', case_number
      !     write(io_unit,*) 'Initial Airspeed:', V_mag 
      !     write(io_unit,*) 'Initial Height', height 
      !   end if 

      !   ! CALCULATE GRAVITY
      !   gravity = gravity_English(-height)

      !   ! SET INITIAL GUESSSES TO ZERO
      !   G        = 0.0
      !   p        = 0.0
      !   q        = 0.0
      !   r        = 0.0
      !   alpha    = G(1)
      !   beta     = G(2)
      !   delta_a  = G(3)
      !   delta_e  = G(4)
      !   delta_r  = G(5)
      !   throttle = G(6)

      !   if (trim_verbose) then
      !     write(io_unit,*) 'Trimming Aircraft for ', trim_type
      !     write(io_unit,'(A,ES20.12)') '  --> Azimuth angle set to psi [deg] =', euler(3) * 180 / pi
      !     if (elevation_angle0 /= -999.0) then 
      !       write(io_unit,'(A,ES20.12)') '  --> Elevation angle set to theta [deg] =', elevation_angle0 * 180 / pi
      !     else 
      !       write(io_unit,'(A,ES20.12)') '  --> Elevation angle set to theta [deg] =', elevation_angle0
      !     end if 
      !     if (sideslip_angle0 /= -999.0) then 
      !       write(io_unit,'(A,ES20.12)') '  --> Sideslip angle set to beta [deg] =', sideslip_angle0 * 180 / pi 
      !     else 
      !       write(io_unit,'(A,ES20.12)') '  --> Bank angle set to phi [deg] =', bank_angle0 * 180 / pi
      !     end if
      !     write(io_unit,'(A)') ''
      !     if (elevation_angle0 /= -999.0) then 
      !       write(io_unit,'(A,ES20.12)') 'Initial theta [deg] =', elevation_angle0 * 180 / pi 
      !       write(io_unit,'(A,ES20.12)') 'Initial gamma [deg] =', climb_angle0 
      !     else 
      !       write(io_unit,'(A,ES20.12)') 'Initial theta [deg] =', elevation_angle0
      !       write(io_unit,'(A,ES20.12)') 'Initial gamma [deg] =', climb_angle0 * 180 / pi 
      !     end if 
      !     write(io_unit,'(A,ES20.12)') 'Initial phi [deg]   =', bank_angle0 * 180 / pi
      !     write(io_unit,'(A,ES20.12)') 'Initial beta [deg]  =', beta * 180 / pi 
      !     write(io_unit,'(A)') ''
      !     write(io_unit,'(A)') 'Newton Solver Settings:'
      !     write(io_unit,'(A,ES20.12)') 'Finite Difference Step Size =', finite_difference_step_size
      !     write(io_unit,'(A,ES20.12)') '          Relaxation Factor =', relaxation_factor
      !     write(io_unit,'(A,ES20.12)') '                  Tolerance =', tolerance
      !     write(io_unit,'(A)') ''
      !     write(io_unit,'(A)') ''
      !   end if

      !   ! PULL OUT AZIMUTH, BANK, AND ELEVATION ANGLES
      !   azimuth_angle = euler(3)
        
      !   if (any(case_number == [4,6])) then 
      !     beta = sideslip_angle0
      !     euler(1) = 0.0
      !   else 
      !     euler(1) = bank_angle0
      !   end if 

      !   if (any(case_number == [1,3,4])) then 
      !     euler(2) = elevation_angle0
      !   end if 

      !   ! LOOP FOR FINDING TRIM STATE
      !   error = 1.0
      !   iteration = 1
      !   do while (error > tolerance)
      !     if (iteration > max_iterations) then 
      !       write(*,*) 'No trim solutions... Reached max iterations...'
      !     end if     
      !     ! DEFINE COS AND SIN TERMS TO SAVE TIME
      !     alpha = G(1)
      !     if (any(case_number == [1,2,3,5])) then
      !       beta = G(2)
      !     else 
      !       euler(1) = G(2)
      !     end if 

      !     ca = cos(alpha)
      !     cb = cos(beta) 
      !     sa = sin(alpha) 
      !     sb = sin(beta)
      !     c_bank = cos(euler(1))
      !     s_bank = sin(euler(1))


      !     ! CALCULATE VELOCITIES FROM 3.4.12
      !     velocities = V_mag * (/ca*cb, sb, sa*cb/) 

      !     u = velocities(1)
      !     v = velocities(2)
      !     w = velocities(3)

      !     ! CALCULATE ELEVATION ANGLE IF CLIMB ANGLE SPECIFIED
      !     if (any(case_number == [2,5,6])) then 
      !       if (trim_verbose) then 
      !         write(io_unit,*) ' '            
      !         write(io_unit,*) 'Solving for elevation angle given a climb angle:'
      !       end if 
      !       climb_angle = climb_angle0
      !       cgamma = cos(climb_angle)
      !       sgamma = sin(climb_angle)

      !         theta1 = asin((u*V_mag*sgamma + (v*s_bank + w*c_bank) * & 
      !           sqrt(u**2 + (v*s_bank + w*c_bank)**2 - V_mag**2*sgamma**2)) &
      !           / (u**2 + (v*s_bank + w*c_bank)**2))

      !         theta2 = asin((u*V_mag*sgamma - (v*s_bank + w*c_bank) * & 
      !                   sqrt(u**2 + (v*s_bank + w*c_bank)**2 - V_mag**2*sgamma**2)) &
      !                   / (u**2 + (v*s_bank + w*c_bank)**2))      
              
      !         solution1 = u*sin(theta1) - (v*s_bank + w*c_bank)*cos(theta1)
      !         solution = V_mag * sgamma 

      !         if (abs(solution1-solution) < tol) then 
      !           euler(2) = theta1 
      !         else 
      !           euler(2) = theta2 
      !         end if 

      !         if (trim_verbose) then
      !           write(io_unit,*) '        theta 1 [deg] =', theta1 * 180 / pi 
      !           write(io_unit,*) '        theta 2 [deg] =', theta2 * 180 / pi 
      !           write(io_unit,*) '  Correct theta [deg] =', euler(2) * 180 / pi 
      !           write(io_unit,*) '  Correct theta [rad] =', euler(2) 
      !           write(io_unit,*) ' '              
      !         end if 
      !     end if 

      !     s_elev = sin(euler(2))
      !     c_elev = cos(euler(2)) 

      !     ! CALCULATE THE ANGULAR RATES
      !     if (any(case_number == [1,2])) then
      !       angular_rates = (gravity * s_bank * c_elev) / (u*c_elev*c_bank + w*s_elev) &
      !                       * (/-s_elev, s_bank*c_elev, c_bank*c_elev/)
      !     else if (any(case_number == [3,4,5,6])) then
      !       angular_rates = 0.0

      !     end if 
      !     p = angular_rates(1)
      !     q = angular_rates(2)
      !     r = angular_rates(3)

      !     if (trim_verbose) then 
      !       write(io_unit,*) 'Updating rotation rates for ', trim_type
      !       write(io_unit,'(A,ES20.12)') 'p [deg/s] = ', p * 180 / pi 
      !       write(io_unit,'(A,ES20.12)') 'q [deg/s] = ', q * 180 / pi 
      !       write(io_unit,'(A,ES20.12)') 'r [deg/s] = ', r * 180 / pi
      !       write(io_unit, '(A)') ''
      !     end if 
      !     res = calc_r(V_mag, height, euler, angular_rates, G)
          
      !     ! STATE THE DEFINITION OF G
      !     if (any(case_number == [4,6])) then 
      !       if (trim_verbose) then 
      !         write(io_unit,'(A)') 'G defined as G = [alpha, bank_angle, aileron, elevator, rudder, throttle]'
      !         write(io_unit, '(A,6(1X,ES20.12))') '      G =', (G(k), k=1,6)
      !         write(io_unit, '(A,6(1X,ES20.12))') '      R =', (res(k), k=1,6)
      !         write(io_unit, '(A)') ''
      !       end if
      !     else 
      !       if (trim_verbose) then 
      !         write(io_unit,'(A)') 'G defined as G = [alpha, beta, aileron, elevator, rudder, throttle]'
      !         write(io_unit, '(A,6(1X,ES20.12))') '      G =', (G(k), k=1,6)
      !         write(io_unit, '(A,6(1X,ES20.12))') '      R =', (res(k), k=1,6)
      !         write(io_unit, '(A)') ''
      !       end if
      !     end if 
          
      !     ! SOLVE FOR G 
      !     call newtons_method(V_mag, height, euler, angular_rates, G)

      !     ! CALCULATE THE ITERATION ERROR
      !     iteration_residual = calc_r(V_mag, height, euler, angular_rates, G)
      !     error = maxval(abs(iteration_residual))

      !     if (trim_verbose) then 
      !       write(io_unit,'(A)') 'New G:'
      !       write(io_unit,'(A,6(1X,ES20.12))') '      G =', (G(k),   k=1,6)
      !       write(io_unit,'(A,6(1X,ES20.12))') '      R =', (iteration_residual(k), k=1,6) 
      !       write(io_unit,*) ''
      !       write(io_unit, '(A)') &
      !         'Iteration   Residual           alpha[deg]           beta[deg]            '// & 
      !         'p[deg/s]             q[deg/s]             r[deg/s]   ' // &
      !         '          phi[deg]             theta[deg]           aileron[deg]         elevator[deg]   ' // &
      !         '     rudder[deg]          throttle[]'
      !       if (sideslip_angle0 /= -999.0) then 
      !         write(io_unit,'(I6,1X,12(1X,ES20.12))') iteration, error, G(1)*180/pi, sideslip_angle0 * 180 / pi, &
      !           p * 180 / pi, q * 180 / pi, r * 180 / pi, euler(1) * 180 / pi, euler(2) * 180 / pi, G(3)*180/pi, &
      !           G(4)*180/pi, G(5)*180/pi, G(6)
      !       else 
      !         write(io_unit,'(I6,1X,12(1X,ES20.12))') iteration, error, G(1)*180/pi, G(2)*180/pi, &
      !           p * 180 / pi, q * 180 / pi, r * 180 / pi, euler(1) * 180 / pi, euler(2) * 180 / pi, G(3)*180/pi, &
      !           G(4)*180/pi, G(5)*180/pi, G(6)
      !       end if 
      !     end if 

      !     ! ENSURE THROTTLE IS IN BOUNDS
      !     if (G(6) > 1.0) then 
      !       if (trim_verbose) then 
      !         write(io_unit,*) 'Overwriting throttle > 1.'
      !       G(6) = 1.0
      !       end if 
      !     else if (G(6) < 0.0) then 
      !       if (trim_verbose) then 
      !         write(io_unit,*) 'Overwriting throttle < 0.'
      !         write(io_unit,*) ''
      !       end if 
      !       G(6) = 0.0
      !     end if 

      !     iteration = iteration + 1
      !   end do

      !   ! SAVE THE TRIM STATE AS THE NEW INITIAL CONDITIONS
      !   ! SET CONTROLS
      !   controls(1:4) = G(3:6)

      !   ! PULL OUT ALPHA
      !   alpha = G(1)
      !   ca = cos(alpha)
      !   sa = sin(alpha)     

      !   ! PULL OUT BETA
      !   if (any(case_number == [4,6])) then 
      !     cb = cos(sideslip_angle0)
      !     sb = sin(sideslip_angle0)
      !     euler(1) = G(2)
      !   else 
      !     cb = cos(G(2))
      !     sb = sin(G(2))
      !   end if

      !   ! CALCUALTE INITIAL STATES
      !   initial_state(1:3)   = V_mag * (/ca*cb, sb, sa*cb/) 
      !   initial_state(4:6)   = angular_rates
      !   initial_state(9)     = height
      !   initial_state(10:13) = euler_to_quat(euler)    

      !   ! UPDATE DESIRED STATE FOR CONTROLLER
      !   xd = initial_state
      !   ud = controls

      !   if (trim_verbose) then 
      !     write(io_unit,*) '---------------------- Trim Solution ----------------------'
      !     write(io_unit,'(A30,F20.13)') '       azimuth_angle[deg]  :', euler(3)
      !     write(io_unit,'(A30,F20.13)') '       elevation_angle[deg]:', euler(2)
      !     write(io_unit,'(A30,F20.13)') '       bank_angle[deg]     :', euler(1) 
      !     write(io_unit,'(A30,F20.13)') '       alpha[deg]          :', alpha 
      !     write(io_unit,'(A30,F20.13)') '       beta[deg]           :', beta 
      !     write(io_unit,'(A30,F20.13)') '       p[deg]              :', p 
      !     write(io_unit,'(A30,F20.13)') '       q[deg]              :', q
      !     write(io_unit,'(A30,F20.13)') '       r[deg]              :', r
      !     write(io_unit,'(A30,F20.13)') '       p_w[deg]            :', initial_state(4)
      !     write(io_unit,'(A30,F20.13)') '       q_w[deg]            :', initial_state(5)
      !     write(io_unit,'(A30,F20.13)') '       r_w[deg]            :', initial_state(6)
      !     write(io_unit,'(A30,F20.13)') '       aileron[deg]        :', controls(1)
      !     write(io_unit,'(A30,F20.13)') '       elevator[deg]       :', controls(2)
      !     write(io_unit,'(A30,F20.13)') '       rudder[deg]         :', controls(3)
      !     write(io_unit,'(A30,F20.13)') '       throttle[deg]       :', controls(4)
      !     write(io_unit,'(A30,F20.13)') '       Climb Angle[deg]    :', climb_angle
      !   end if 

      !   if (exam_answers) then 
      !     write(io_unit,*) '---------------------- Exam Answers ----------------------'
      !     write(io_unit,'(A30,ES25.13E3)') '       elevation_angle[deg]:', euler(2)
      !     write(io_unit,'(A30,ES25.13E3)') '       bank_angle[deg]     :', euler(1)
      !     write(io_unit,'(A30,ES25.13E3)') '       alpha[deg]          :', alpha
      !     write(io_unit,'(A30,ES25.13E3)') '       beta[deg]           :', beta
      !     write(io_unit,'(A30,ES25.13E3)') '       p[deg]              :', p
      !     write(io_unit,'(A30,ES25.13E3)') '       q[deg]              :', q
      !     write(io_unit,'(A30,ES25.13E3)') '       r[deg]              :', r
      !     write(io_unit,'(A30,ES25.13E3)') '       aileron[deg]        :', controls(1)
      !     write(io_unit,'(A30,ES25.13E3)') '       elevator[deg]       :', controls(2)
      !     write(io_unit,'(A30,ES25.13E3)') '       rudder[deg]         :', controls(3)
      !     write(io_unit,'(A30,ES25.13E3)') '       throttle[deg]       :', controls(4)

      !   end if 
      ! end function trim_algorithm

    

    !=========================
    ! Calculate Residual
      function calc_r(t, G) result(R)
        implicit none
        type(vehicle_t) :: t
        real, intent(in) :: G(6)
        real :: ca, cb, sa, sb
        real :: R(6), temp_state(13), dummy_res(13)

        ! Pull out controls
        t%controls(1:4) = G(3:6)

        ! Pull out alpha and beta
        ca = cos(G(1))
        sa = sin(G(1))     
        cb = cos(G(2))
        sb = sin(G(2))

        ! Set states
        temp_state(1:3)   = t%init_airspeed * (/ca*cb, sb, sa*cb/) 
        temp_state(4:6)   = 0.0
        temp_state(9)     = -t%init_alt
        temp_state(7:8)   = 0.0
        temp_state(10:13) = euler_to_quat(t%init_eul)

        ! Set controls 
        t%controls = G(3:6) 

        ! Calculate residual
        dummy_res = diff_eq(t, 0.0, temp_state)
        R = dummy_res(1:6)

      end function calc_r

end module vehicle_m