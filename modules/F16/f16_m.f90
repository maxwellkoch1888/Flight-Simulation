module f16_m
    use koch_m
    use jsonx_m
    use micro_time_m
    use linalg_mod
    use connection_m
    implicit none
  
    ! JSON POINTER
    type(json_value), pointer :: j_main

    ! UDP POINTERS
    type(connection) :: graphics
    type(connection) :: user_controls

    ! BUILD GLOBAL VARIABLES FOR THE MODULE
    real :: mass
    real :: inertia(3,3)
    real :: inertia_inv(3,3)
    real :: h(3)
    real :: FM(6)
    real :: initial_state(13)
    real :: controls(4)
    real :: rho0
    real :: T0, Ta
    real, allocatable :: aero_ref_location(:)
    integer :: io_unit

    ! DEFINE AERODYNAMIC PROPERTIES
    real:: planform_area, longitudinal_length, lateral_length, sweep
    real:: CL0, CL_alpha, CL_alphahat, CL_qbar, CL_elevator
    real:: CS_beta, CS_pbar, CS_alpha_pbar, CS_rbar, CS_aileron, CS_rudder
    real:: CD_L0, CD_L1, CD_L1_L1, CD_CS_CS, CD_qbar, CD_alpha_qbar, CD_elevator, CD_alpha_elevator, CD_elevator_elevator
    real:: Cl_beta, Cl_pbar, Cl_rbar, Cl_alpha_rbar, Cl_aileron, Cl_rudder
    real:: Cm_0, Cm_alpha, Cm_qbar, Cm_alphahat, Cm_elevator
    real:: Cn_beta, Cn_pbar, Cn_alpha_pbar, Cn_rbar, Cn_aileron, Cn_alpha_aileron, Cn_rudder
    logical :: compressibility, rk4_verbose, print_states

    ! ADD VARIABLES FOR TRIM ALGORITHM
    character(:), allocatable :: sim_type
    character(:), allocatable :: trim_type
    real :: relaxation_factor, tolerance, max_iterations, finite_difference_step_size
    real :: bank_angle0, sideslip_angle0, climb_angle0, elevation_angle0
    logical :: trim_verbose, exam_answers


    contains
  !=========================
  ! RK4 Integrator
    function rk4(t0, y1, delta_t) result(state)
      implicit none
      real, intent(in) :: t0, delta_t, y1(13)
      real, dimension(13) :: state, k1, k2, k3, k4

      ! DEFINE THE K TERMS FOR RK4 METHOD
      if(rk4_verbose) then 
        write(io_unit,*) "diff_eq function called: "
        write(io_unit,*) "              RK4 call number =  1"
      end if
      k1 = differential_equations(t0, y1)
      if(rk4_verbose) then 
        write(io_unit,*) "diff_eq function called: "
        write(io_unit,*) "              RK4 call number =  2"
      end if
      k2 = differential_equations(t0 + delta_t*0.5, y1 + k1 * delta_t*0.5)
      if(rk4_verbose) then 
        write(io_unit,*) "diff_eq function called: "
        write(io_unit,*) "              RK4 call number =  3"
        write(io_unit,*) " initial state", y1
        write(io_unit,*) " k2", k2   
      end if
      k3 = differential_equations(t0 + delta_t*0.5, y1 + k2 * delta_t*0.5)
      if(rk4_verbose) then 
        write(io_unit,*) "diff_eq function called: "
        write(io_unit,*) "              RK4 call number =  4"
      end if
      k4 = differential_equations(t0 + delta_t, y1 + k3 * delta_t)


      ! DEFINE THE RESULT FROM RK4
      state = y1 + (delta_t/6) * (k1 + 2*k2 + 2*k3 + k4)

    end function rk4

  !=========================
  ! Equations of Motion: (/u,v,w, p,q,r, x,y,z, e0,ex,ey,ez/)
    function differential_equations(t, state) result(dstate_dt)
      implicit none 
      real, intent(in) :: t, state(13) 
      real :: dstate_dt(13) 
      real :: u, v, w, p, q, r, x, y, z, e0, ex, ey, ez 
      real :: acceleration(3), angular_accelerations(3), velocity(3), quat_change(4) 
      real :: velocity_quat(4), quat_inv(4) 
      real :: orientation_effect(3), angular_v_effect(3), gyroscopic_change(3)
      real :: inertia_effects(3), wind_velocity(3), gyroscopic_effects(3,3) 
      real :: Ixx, Iyy, Izz, Ixy, Ixz, Iyz 
      real :: quat_matrix(4,3) 
      real :: hxb, hyb, hzb, hxb_dot, hyb_dot, hzb_dot
      real :: Vxw, Vyw, Vzw, gravity_ft_per_sec2 
      real :: v1(4), v2(4), angular_inertia(3,3)
      real :: avoid_warning

      avoid_warning = t

      ! UNPACK STATES
      u  = state(1)  
      v  = state(2)  
      w  = state(3)
      p  = state(4)  
      q  = state(5)  
      r  = state(6)
      x  = state(7)  
      y  = state(8)  
      z  = state(9)
      e0 = state(10) 
      ex = state(11) 
      ey = state(12) 
      ez = state(13) 

      ! UNPACK INERTIA
      Ixx = inertia(1,1); Iyy = inertia(2,2); Izz = inertia(3,3)
      Ixy = inertia(1,2); Ixz = inertia(1,3); Iyz = inertia(2,3)
    
      ! CALCULATE FORCES AND MOMENTS
      call pseudo_aero(state)
      gravity_ft_per_sec2 = gravity_English(-state(9))

      ! SET GYROSCOPIC EFFECTS
      hxb = h(1)
      hyb = h(2)
      hzb = h(3)

      ! SET GYROSCOPIC CHANGE AND WIND VELOCITY TO ZERO
      hxb_dot = 0.0
      hyb_dot = 0.0
      hzb_dot = 0.0
      Vxw = 0.0
      Vyw = 0.0
      Vzw = 0.0

      ! BUILD MATRICES/ VECTORS USED IN DIFFERENTIAL EQUATION
      orientation_effect =  (/2 * (ex*ez - ey*e0), &
                              2 * (ey*ez + ex*e0), & 
                              ez**2 + e0**2 - ex**2 - ey**2/)

      angular_v_effect =    (/r*v - q*w, p*w - r*u, q*u - p*v/)

      gyroscopic_effects = reshape((/ 0.0,  hzb,  -hyb, &
                                     -hzb,  0.0,   hxb, &
                                      hyb, -hxb,    0.0 /), (/3,3/))

      angular_inertia = reshape((/  Ixx, -Ixy, -Ixz, &
                                   -Ixy,  Iyy, -Iyz, &
                                   -Ixz, -Iyz,  Izz /), (/3,3/))

      gyroscopic_change =   (/hxb_dot, hyb_dot, hzb_dot/)
      
      inertia_effects =     (/(Iyy - Izz)*q*r + Iyz*(q**2 -r**2) + Ixz*p*q - Ixy*p*r, &
                              (Izz - Ixx)*p*r + Ixz*(r**2 -p**2) + Ixy*q*r - Iyz*p*q, &
                              (Ixx - Iyy)*p*q + Ixy*(p**2 -q**2) + Iyz*p*r - Ixz*q*r/)
      wind_velocity = (/Vxw, Vyw, Vzw/)

      velocity_quat = (/0.0, u, v, w/)
      quat_inv = (/e0, -ex, -ey, -ez/)
      quat_matrix = reshape((/ -ex,  e0,  ez, -ey, &
                               -ey, -ez,  e0,  ex, &
                               -ez,  ey, -ex,  e0  /), (/4,3/))             

      ! BUILD THE DIFFERENTIAL EQUATIONS
      ! ACCELERATION IN BODY FRAME
      acceleration = 1.0 / mass * FM(1:3) + gravity_ft_per_sec2 * orientation_effect + angular_v_effect

      ! ROLL, PITCH, YAW ACCELERATIONS
      angular_accelerations = matmul(matrix_inv(angular_inertia) , FM(4:6) + &
                              matmul(gyroscopic_effects, (/p, q, r/)) + inertia_effects - gyroscopic_change)

      ! VELOCITY IN THE INERITAL FRAME
      v1 = quat_mult(state(10:13), (/0.0, state(1:3)/))
      v2 = quat_mult(v1, quat_inv)
      velocity = v2(2:4) + wind_velocity

      ! AIRCRAFT ORIENTATION RATE OF CHANGE
      quat_change = 0.5 * matmul(quat_matrix, (/p, q, r/))


      ! RETURN THE STATE DERIVATIVE
      dstate_dt(1:3)  = acceleration
      dstate_dt(4:6)  = angular_accelerations
      dstate_dt(7:9)  = velocity
      dstate_dt(10:13) = quat_change

      ! PRINT STEPS IF SPECIFIED
      if(rk4_verbose) then
        write(io_unit,'(A,13(1X,ES20.12))') "                     time [s] =", t
        write(io_unit,'(A,*(1X,ES20.12))') '    State vector coming in    =', state, controls

        write(io_unit,'(A,6 (1X,ES20.12))') "    Pseudo aerodynamics (F,M) =", FM
        ! write(io_unit,*) "v1", v1
        ! write(io_unit,*) "v2", v2
        ! write(io_unit,*) "velocity", velocity
        ! write(io_unit,*) "quat change", quat_change
        write(io_unit,'(A,13(1X,ES20.12))') "    Diff. Eq. results         =", dstate_dt
        write(io_unit,*) ''
      end if

      end function differential_equations

  !=========================
  ! Aerodynamic Forces and Moments for f16
    subroutine pseudo_aero(state)
      implicit none
      real, intent(in) :: state(13)
      real :: Re, geometric_altitude_ft, geopotential_altitude_ft
      real :: temp_R, pressure_lbf_per_ft2, density_slugs_per_ft3
      real :: dyn_viscosity_slug_per_ft_sec, sos_ft_per_sec
      real :: V, Uc(3)
      real :: Cl_pitch0, alpha, beta, beta_f, pbar, qbar, rbar, angular_rates(3)
      real :: CL, CL1, CD, CS, L, D, S, Cl_pitch, Cm, Cn
      real :: CM1, CM2, mach_num, gamma, R
      real :: ca, cb, sa, sb
      real :: alpha_hat, beta_hat
      real :: delta_a, delta_e, delta_r
      real :: T, throttle
        
      ! BUILD THE ATMOSPHERE 
      geometric_altitude_ft = -state(9)
      call std_atm_English(&
        geometric_altitude_ft, geopotential_altitude_ft,     & 
        temp_R, pressure_lbf_per_ft2, density_slugs_per_ft3, & 
        dyn_viscosity_slug_per_ft_sec, sos_ft_per_sec)

      ! CALCULATE VELOCITY UNIT VECTOR
      V =  (state(1)**2 + state(2)**2 + state(3)**2)**0.5
      Uc = (/state(1:3)/) / V

      ! CALCULATE ALPHA AND BETA 3.4.4 and 3.4.5
      alpha =  atan2(state(3) , state(1))
      beta =   asin(state(2) / V)
      beta_f = atan2(state(2) , state(1))

      ! CALCULATE ALPHA_HAT USING EQN 3.4.20
      ! alpha_hat = dalpha_dtime * longitudinal_length / (2 * V)
      ! beta_ hat = dbeta_dtime * lateral_length / (2 * V)
      alpha_hat = 0.0
      beta_hat = 0.0

      ! CALCULATE PBAR, QBAR, AND RBAR from eq 3.4.23
      angular_rates = 1 / (2*V) * (/state(4) * lateral_length, state(5) * longitudinal_length, state(6) * lateral_length/)
      pbar = angular_rates(1)
      qbar = angular_rates(2)
      rbar = angular_rates(3)

      ! CALCULATE THE REYNOLDS NUMBER
      Re = density_slugs_per_ft3 * V * 2 * longitudinal_length / dyn_viscosity_slug_per_ft_sec

      ! PULL OUT CONTROLS
      delta_a = controls(1)
      delta_e = controls(2)
      delta_r = controls(3)
      throttle = controls(4)

      ! CALCULATE THE LIFT, DRAG, AND SIDE COEFFICIENTS
      CL1 =  CL0 + CL_alpha * alpha
      CL  = CL1 + CL_qbar * qbar + CL_alphahat * alpha_hat + CL_elevator * delta_e
      CS  = CS_beta * beta + (CS_pbar + CS_alpha_pbar * alpha) * pbar + CS_rbar * rbar &
           + CS_aileron * delta_a + CS_rudder * delta_r
      CD  =  CD_L0 + CD_L1 * CL1 + CD_L1_L1 * CL1 **2 + CD_CS_CS * CS **2 &
            + (CD_qbar + CD_alpha_qbar * alpha) * qbar + (CD_elevator + CD_alpha_elevator * alpha) &
            * delta_e + CD_elevator_elevator * delta_e ** 2

      ! ACCOUNT FOR COMPRESSIBILITY IF SPECIFIED
      ! CHAPTER 3.10 IN BOOK, PRANDTL-GLAUERT CORRECTION
      if (compressibility == .true.) then 
        CM1 = 2.13/ (sweep + 0.15)**2
        CM2 = 15.35*sweep**2 - 19.64*sweep +16.86

        CL = CL / (sqrt(1-sos_ft_per_sec**2))
        CS = CS / (sqrt(1-sos_ft_per_sec**2))
        CD = CD * (1.0 + CM1 * sos_ft_per_sec**CM2)
      end if 

      L =   CL * (0.5 * density_slugs_per_ft3 * V **2 * planform_area)
      S =   CS * (0.5 * density_slugs_per_ft3 * V **2 * planform_area)
      D =   CD * (0.5 * density_slugs_per_ft3 * V **2 * planform_area)

      ca = cos(alpha)
      cb = cos(beta)
      sa = sin(alpha)
      sb = sin(beta)

      ! TABLE 3.4.4
      FM(1) = (D*ca*cb + S*ca*sb - L*sa) * (-1.0)
      FM(2) = (S*cb - D*sb)
      FM(3) = (D*sa*cb + S*sa*sb + L*ca) * (-1.0)

      ! ADD THE ENGINE THRUST
      if (throttle < 0.0) then
        throttle = 0.0
      else if (throttle > 1.0) then
        throttle = 1.0
      else 
        throttle = throttle 
      end if 

      T = throttle * T0 * (density_slugs_per_ft3/rho0) ** Ta
      FM(1) = FM(1) + T

      ! CALCULATE THE ROLL, PITCH, AND YAW COEFFICIENTS
      Cl_pitch = Cl_beta * beta + Cl_pbar * pbar + (Cl_rbar + Cl_alpha_rbar * alpha) * rbar &
                 + Cl_aileron * delta_a + Cl_rudder * delta_r  ! roll moment
      Cm =       Cm_0 + Cm_alpha * alpha + Cm_qbar * qbar + Cm_alphahat * alpha_hat + Cm_elevator * delta_e ! pitch moment
      Cn =       Cn_beta * beta + (Cn_pbar + Cn_alpha_pbar * alpha) * pbar + Cn_rbar * rbar &
                 + (Cn_aileron + Cn_alpha_aileron * alpha) * delta_a + Cn_rudder * delta_r ! yaw moment

      ! ACCOUNT FOR COMPRESSIBILITY IF SPECIFIED
      if (compressibility == .true.) then 
        cl_pitch = cl_pitch / (sqrt(1-sos_ft_per_sec**2))
        Cm       = Cm       / (sqrt(1-sos_ft_per_sec**2))
        Cn       = Cn       / (sqrt(1-sos_ft_per_sec**2))
      end if 
      
      FM(4) = Cl_pitch * (0.5 * density_slugs_per_ft3 * V **2 * planform_area * lateral_length)
      FM(5) = Cm       * (0.5 * density_slugs_per_ft3 * V **2 * planform_area * longitudinal_length)
      FM(6) = Cn       * (0.5 * density_slugs_per_ft3 * V **2 * planform_area * lateral_length)

      ! SHIFT CG LOCATION
      FM(4:6) = FM(4:6) + cross_product(aero_ref_location, FM(1:3))
      
      ! PRINT STEPS IF SPECIFIED
      if(rk4_verbose) then
        write(io_unit,*) "V_mag = ", V
        write(io_unit,*) "pbar = ", pbar
        write(io_unit,*) "qbar = ", qbar
        write(io_unit,*) "rbar = ", rbar
        write(io_unit,*) "cos(alpha) = ", cos(alpha)
        write(io_unit,*) "sin(alpha) = ", sin(alpha)
        write(io_unit,*) "cos(beta) = ", cos(beta)
        write(io_unit,*) "sin(beta) = ", sin(beta)
        write(io_unit,*) "alpha = ", alpha
        write(io_unit,*) "beta = ", beta
        write(io_unit,*) "beta_flank = ", beta_f
        write(io_unit,*) "CL1", CL1
        write(io_unit,*) "CL", CL
        write(io_unit,*) "CS", CS
        write(io_unit,*) "CD", CD
        write(io_unit,*) "C_l", Cl_pitch
        write(io_unit,*) "Cm", Cm
        write(io_unit,*) "Cn", Cn
        write(io_unit,*) "delta_a", delta_a * 180.0 / pi
        write(io_unit,*) "delta_e", delta_e * 180.0 / pi
        write(io_unit,*) "delta_r", delta_r * 180.0 / pi
        write(io_unit,*) "T", T
      end if 

    end subroutine pseudo_aero

  !=========================
  ! Mass and Inertia
    subroutine mass_inertia()
      implicit none
      real :: weight
      real :: Ixx, Iyy, Izz, Ixy, Ixz, Iyz

      ! READ CG SHIFT IF APPLICABLE
      call jsonx_get(j_main, 'vehicle.aero_ref_location[ft]', aero_ref_location, 0.0, 3)

      ! READ MASS PROPERTIES
      call jsonx_get(j_main, 'vehicle.mass.weight[lbf]',     weight)
      call jsonx_get(j_main, 'vehicle.mass.Ixx[slug-ft^2]',  Ixx)
      call jsonx_get(j_main, 'vehicle.mass.Iyy[slug-ft^2]',  Iyy)
      call jsonx_get(j_main, 'vehicle.mass.Izz[slug-ft^2]',  Izz)
      call jsonx_get(j_main, 'vehicle.mass.Ixy[slug-ft^2]',  Ixy)
      call jsonx_get(j_main, 'vehicle.mass.Ixz[slug-ft^2]',  Ixz)
      call jsonx_get(j_main, 'vehicle.mass.Iyz[slug-ft^2]',  Iyz)
      call jsonx_get(j_main, 'vehicle.mass.hx[slug-ft^2/s]', h(1))
      call jsonx_get(j_main, 'vehicle.mass.hy[slug-ft^2/s]', h(2))
      call jsonx_get(j_main, 'vehicle.mass.hz[slug-ft^2/s]', h(3))

      ! DEFINE IDENTITY MATRIX 
      inertia = reshape((/ Ixx, Ixy, Ixz, &
                           Ixy, Iyy, Iyz, &
                           Ixz, Iyz, Izz /), (/3,3/))

      ! CALCULATE MASS AND INERTIA
      mass = weight / 32.17404855643
      inertia_inv = matrix_inv(inertia)
    end subroutine mass_inertia

  !=========================
  ! Trim Algorithm
    function trim_algorithm(V_mag, height, euler, tolerance, trim_type) result(G)
      implicit none 
      real :: V_mag, height, tolerance
      real :: alpha, beta, p, q, r, delta_a, delta_e, delta_r, throttle 
      real :: elevation_angle, azimuth_angle
      real :: c_bank, c_elev, s_bank, s_elev, ca, cb, sa, sb, error, pw
      real :: u, v, w, velocities(3), gravity
      real :: G(6), res(6), iteration_residual(6)
      real :: angular_rates(3), euler(3), print_statement(13)
      real :: cgamma, sgamma, climb_angle, solution, theta1, theta2, solution1
      integer :: k, iteration, case_number
      character(*), intent(in) :: trim_type
      
      ! DETERMINE THE TRIM TYPE BEING SPECIFIED
      if (trim_type == 'sct') then 
        if (elevation_angle0 /= -999.0) then 
          case_number = 1 ! sct, elevation_angle, bank_angle
        else
          case_number = 2 ! sct, climb_angle, bank_angle
        end if 
      else 
        if (elevation_angle0 /= -999.0) then 
          if (sideslip_angle0 == -999.0) then 
            case_number = 3 ! shss, elevation_angle, bank_angle 
          else 
            case_number = 4 ! shss, elevation_angle, sideslip_angle
          end if 
        else 
          if (sideslip_angle0 == -999.0) then 
            case_number = 5 ! shss, climb_angle, bank_angle
          else 
            case_number = 6 ! shss, climb_angle, sideslip_angle
          end if 
        end if 
      end if 
      
      if (trim_verbose) then 
        write(io_unit,*) 'Case Number:', case_number
        write(io_unit,*) 'Initial Airspeed:', V_mag 
        write(io_unit,*) 'Initial Height', height 
      end if 

      ! CALCULATE GRAVITY
      gravity = gravity_English(-height)

      ! SET INITIAL GUESSSES TO ZERO
      G        = 0.0
      p        = 0.0
      q        = 0.0
      r        = 0.0
      alpha    = G(1)
      beta     = G(2)
      delta_a  = G(3)
      delta_e  = G(4)
      delta_r  = G(5)
      throttle = G(6)

      if (trim_verbose) then
        write(io_unit,*) 'Trimming Aircraft for ', trim_type
        write(io_unit,'(A,ES20.12)') '  --> Azimuth angle set to psi [deg] =', euler(3) * 180 / pi
        if (elevation_angle0 /= -999.0) then 
          write(io_unit,'(A,ES20.12)') '  --> Elevation angle set to theta [deg] =', elevation_angle0 * 180 / pi
        else 
          write(io_unit,'(A,ES20.12)') '  --> Elevation angle set to theta [deg] =', elevation_angle0
        end if 
        if (sideslip_angle0 /= -999.0) then 
          write(io_unit,'(A,ES20.12)') '  --> Sideslip angle set to beta [deg] =', sideslip_angle0 * 180 / pi 
        else 
          write(io_unit,'(A,ES20.12)') '  --> Bank angle set to phi [deg] =', bank_angle0 * 180 / pi
        end if
        write(io_unit,'(A)') ''
        if (elevation_angle0 /= -999.0) then 
          write(io_unit,'(A,ES20.12)') 'Initial theta [deg] =', elevation_angle0 * 180 / pi 
          write(io_unit,'(A,ES20.12)') 'Initial gamma [deg] =', climb_angle0 
        else 
          write(io_unit,'(A,ES20.12)') 'Initial theta [deg] =', elevation_angle0
          write(io_unit,'(A,ES20.12)') 'Initial gamma [deg] =', climb_angle0 * 180 / pi 
        end if 
        write(io_unit,'(A,ES20.12)') 'Initial phi [deg]   =', bank_angle0 * 180 / pi
        write(io_unit,'(A,ES20.12)') 'Initial beta [deg]  =', beta * 180 / pi 
        write(io_unit,'(A)') ''
        write(io_unit,'(A)') 'Newton Solver Settings:'
        write(io_unit,'(A,ES20.12)') 'Finite Difference Step Size =', finite_difference_step_size
        write(io_unit,'(A,ES20.12)') '          Relaxation Factor =', relaxation_factor
        write(io_unit,'(A,ES20.12)') '                  Tolerance =', tolerance
        write(io_unit,'(A)') ''
        write(io_unit,'(A)') ''
      end if

      ! PULL OUT AZIMUTH, BANK, AND ELEVATION ANGLES
      azimuth_angle = euler(3)
      
      if (any(case_number == [4,6])) then 
        beta = sideslip_angle0
        euler(1) = 0.0
      else 
        euler(1) = bank_angle0
      end if 

      if (any(case_number == [1,3,4])) then 
        euler(2) = elevation_angle0
      end if 

      ! LOOP FOR FINDING TRIM STATE
      error = 1.0
      iteration = 1
      do while (error > tolerance)
        if (iteration > max_iterations) then 
          write(*,*) 'No trim solutions... Reached max iterations...'
        end if     
        ! DEFINE COS AND SIN TERMS TO SAVE TIME
        alpha = G(1)
        if (any(case_number == [1,2,3,5])) then
          beta = G(2)
        else 
          euler(1) = G(2)
        end if 

        ca = cos(alpha)
        cb = cos(beta) 
        sa = sin(alpha) 
        sb = sin(beta)
        c_bank = cos(euler(1))
        s_bank = sin(euler(1))


        ! CALCULATE VELOCITIES FROM 3.4.12
        velocities = V_mag * (/ca*cb, sb, sa*cb/) 

        u = velocities(1)
        v = velocities(2)
        w = velocities(3)

        ! CALCULATE ELEVATION ANGLE IF CLIMB ANGLE SPECIFIED
        if (any(case_number == [2,5,6])) then 
          if (trim_verbose) then 
            write(io_unit,*) ' '            
            write(io_unit,*) 'Solving for elevation angle given a climb angle:'
          end if 
          climb_angle = climb_angle0
          cgamma = cos(climb_angle)
          sgamma = sin(climb_angle)

            theta1 = asin((u*V_mag*sgamma + (v*s_bank + w*c_bank) * & 
              sqrt(u**2 + (v*s_bank + w*c_bank)**2 - V_mag**2*sgamma**2)) &
              / (u**2 + (v*s_bank + w*c_bank)**2))

            theta2 = asin((u*V_mag*sgamma - (v*s_bank + w*c_bank) * & 
                      sqrt(u**2 + (v*s_bank + w*c_bank)**2 - V_mag**2*sgamma**2)) &
                      / (u**2 + (v*s_bank + w*c_bank)**2))      
            
            solution1 = u*sin(theta1) - (v*s_bank + w*c_bank)*cos(theta1)
            solution = V_mag * sgamma 

            if (abs(solution1-solution) < tol) then 
              euler(2) = theta1 
            else 
              euler(2) = theta2 
            end if 

            if (trim_verbose) then
              write(io_unit,*) '        theta 1 [deg] =', theta1 * 180 / pi 
              write(io_unit,*) '        theta 2 [deg] =', theta2 * 180 / pi 
              write(io_unit,*) '  Correct theta [deg] =', euler(2) * 180 / pi 
              write(io_unit,*) '  Correct theta [rad] =', euler(2) 
              write(io_unit,*) ' '              
            end if 
        end if 

        s_elev = sin(euler(2))
        c_elev = cos(euler(2)) 


        ! CALCULATE THE ANGULAR RATES
        if (any(case_number == [1,2])) then
          angular_rates = (gravity * s_bank * c_elev) / (u*c_elev*c_bank + w*s_elev) &
                          * (/-s_elev, s_bank*c_elev, c_bank*c_elev/)
        else if (any(case_number == [3,4,5,6])) then
          angular_rates = 0.0

        end if 
        p = angular_rates(1)
        q = angular_rates(2)
        r = angular_rates(3)

        if (trim_verbose) then 
          write(io_unit,*) 'Updating rotation rates for ', trim_type
          write(io_unit,'(A,ES20.12)') 'p [deg/s] = ', p * 180 / pi 
          write(io_unit,'(A,ES20.12)') 'q [deg/s] = ', q * 180 / pi 
          write(io_unit,'(A,ES20.12)') 'r [deg/s] = ', r * 180 / pi
          write(io_unit, '(A)') ''
        end if 
        res = calc_r(V_mag, height, euler, angular_rates, G)
        
        ! STATE THE DEFINITION OF G
        if (any(case_number == [4,6])) then 
          if (trim_verbose) then 
            write(io_unit,'(A)') 'G defined as G = [alpha, bank_angle, aileron, elevator, rudder, throttle]'
            write(io_unit, '(A,6(1X,ES20.12))') '      G =', (G(k), k=1,6)
            write(io_unit, '(A,6(1X,ES20.12))') '      R =', (res(k), k=1,6)
            write(io_unit, '(A)') ''
          end if
        else 
          if (trim_verbose) then 
            write(io_unit,'(A)') 'G defined as G = [alpha, beta, aileron, elevator, rudder, throttle]'
            write(io_unit, '(A,6(1X,ES20.12))') '      G =', (G(k), k=1,6)
            write(io_unit, '(A,6(1X,ES20.12))') '      R =', (res(k), k=1,6)
            write(io_unit, '(A)') ''
          end if
        end if 
        
        ! SOLVE FOR G 
        call newtons_method(V_mag, height, euler, angular_rates, G)

        ! CALCULATE THE ITERATION ERROR
        iteration_residual = calc_r(V_mag, height, euler, angular_rates, G)
        error = maxval(abs(iteration_residual))

        if (trim_verbose) then 
          write(io_unit,'(A)') 'New G:'
          write(io_unit,'(A,6(1X,ES20.12))') '      G =', (G(k),   k=1,6)
          write(io_unit,'(A,6(1X,ES20.12))') '      R =', (iteration_residual(k), k=1,6) 
          write(io_unit,*) ''
          write(io_unit, '(A)') &
            'Iteration   Residual           alpha[deg]           beta[deg]            '// & 
            'p[deg/s]             q[deg/s]             r[deg/s]   ' // &
            '          phi[deg]             theta[deg]           aileron[deg]         elevator[deg]   ' // &
            '     rudder[deg]          throttle[]'
          if (sideslip_angle0 /= -999.0) then 
            write(io_unit,'(I6,1X,12(1X,ES20.12))') iteration, error, G(1)*180/pi, sideslip_angle0 * 180 / pi, &
              p * 180 / pi, q * 180 / pi, r * 180 / pi, euler(1) * 180 / pi, euler(2) * 180 / pi, G(3)*180/pi, &
              G(4)*180/pi, G(5)*180/pi, G(6)
          else 
            write(io_unit,'(I6,1X,12(1X,ES20.12))') iteration, error, G(1)*180/pi, G(2)*180/pi, &
              p * 180 / pi, q * 180 / pi, r * 180 / pi, euler(1) * 180 / pi, euler(2) * 180 / pi, G(3)*180/pi, &
              G(4)*180/pi, G(5)*180/pi, G(6)
          end if 
        end if 

        ! ENSURE THROTTLE IS IN BOUNDS
        if (G(6) > 1.0) then 
          if (trim_verbose) then 
            write(io_unit,*) 'Overwriting throttle > 1.'
          G(6) = 1.0
          end if 
        else if (G(6) < 0.0) then 
          if (trim_verbose) then 
            write(io_unit,*) 'Overwriting throttle < 0.'
            write(io_unit,*) ''
          end if 
          G(6) = 0.0
        end if 

        iteration = iteration + 1
      end do

      ! SAVE THE TRIM STATE AS THE NEW INITIAL CONDITIONS
      ! SET CONTROLS
      controls(1:4) = G(3:6)

      ! PULL OUT ALPHA
      alpha = G(1)
      ca = cos(alpha)
      sa = sin(alpha)     

      ! PULL OUT BETA
      if (any(case_number == [4,6])) then 
        cb = cos(sideslip_angle0)
        sb = sin(sideslip_angle0)
        euler(1) = G(2)
      else 
        cb = cos(G(2))
        sb = sin(G(2))
      end if

      ! CALCUALTE INITIAL STATES
      initial_state(1:3)   = V_mag * (/ca*cb, sb, sa*cb/) 
      initial_state(4:6)   = angular_rates
      initial_state(9)     = height
      initial_state(10:13) = euler_to_quat(euler)    
      if (trim_verbose) then 
        write(io_unit,*) '---------------------- Trim Solution ----------------------'
        write(io_unit,'(A30,F20.13)') '       azimuth_angle[deg]  :', euler(3)
        write(io_unit,'(A30,F20.13)') '       elevation_angle[deg]:', euler(2)
        write(io_unit,'(A30,F20.13)') '       bank_angle[deg]     :', euler(1) 
        write(io_unit,'(A30,F20.13)') '       alpha[deg]          :', alpha 
        write(io_unit,'(A30,F20.13)') '       beta[deg]           :', beta 
        write(io_unit,'(A30,F20.13)') '       p[deg]              :', p 
        write(io_unit,'(A30,F20.13)') '       q[deg]              :', q
        write(io_unit,'(A30,F20.13)') '       r[deg]              :', r
        write(io_unit,'(A30,F20.13)') '       p_w[deg]            :', initial_state(4)
        write(io_unit,'(A30,F20.13)') '       q_w[deg]            :', initial_state(5)
        write(io_unit,'(A30,F20.13)') '       r_w[deg]            :', initial_state(6)
        write(io_unit,'(A30,F20.13)') '       aileron[deg]        :', controls(1)
        write(io_unit,'(A30,F20.13)') '       elevator[deg]       :', controls(2)
        write(io_unit,'(A30,F20.13)') '       rudder[deg]         :', controls(3)
        write(io_unit,'(A30,F20.13)') '       throttle[deg]       :', controls(4)
        write(io_unit,'(A30,F20.13)') '       Climb Angle[deg]    :', climb_angle
      end if 

      if (exam_answers) then 
        write(io_unit,*) '---------------------- Exam Answers ----------------------'
        write(io_unit,'(A30,ES25.13E3)') '       elevation_angle[deg]:', euler(2)
        write(io_unit,'(A30,ES25.13E3)') '       bank_angle[deg]     :', euler(1)
        write(io_unit,'(A30,ES25.13E3)') '       alpha[deg]          :', alpha
        write(io_unit,'(A30,ES25.13E3)') '       beta[deg]           :', beta
        write(io_unit,'(A30,ES25.13E3)') '       p[deg]              :', p
        write(io_unit,'(A30,ES25.13E3)') '       q[deg]              :', q
        write(io_unit,'(A30,ES25.13E3)') '       r[deg]              :', r
        write(io_unit,'(A30,ES25.13E3)') '       aileron[deg]        :', controls(1)
        write(io_unit,'(A30,ES25.13E3)') '       elevator[deg]       :', controls(2)
        write(io_unit,'(A30,ES25.13E3)') '       rudder[deg]         :', controls(3)
        write(io_unit,'(A30,ES25.13E3)') '       throttle[deg]       :', controls(4)

      end if 
      end function trim_algorithm
  !=========================
  ! Newton's Method Solver to find G (alpha, beta, delta_a, delta_e, delta_r, throttle)
    subroutine newtons_method(V_mag, height, euler, angular_rates, G)
      implicit none
      real , allocatable :: res(:), delta_G(:)
      real, intent(inout) :: G(6)
      real :: R_pos(6), R_neg(6), step_size
      real :: angular_rates(3), jac(6,6), euler(3), height, dummy_res(13), V_mag
      integer :: k, i, j

      allocate(res(6), delta_G(6))

      ! CALCULATE JACOBIAN AND RESIDUAL
      if (trim_verbose) then
        write(io_unit,'(A)') 'Building Jacobian Matrix:'
        write(io_unit,'(A)') ''
      end if 

      ! USE CENTRAL DIFFERENCE METHOD TO FIND JACOBIAN
      step_size = finite_difference_step_size
      do j = 1,6        
        G(j) = G(j) + step_size
        R_pos = calc_r(V_mag, height, euler, angular_rates, G)
        
        if (trim_verbose) then
          write(io_unit, '(A,I0,A)') 'Computing gradient relative to G[', j-1, ']'
          write(io_unit, '(A)') '   Positive Finite-Difference Step '
          write(io_unit, '(A,6(1X,ES20.12))') '      G =', (G(i), i=1,6)
          write(io_unit, '(A,6(1X,ES20.12))') '      R =', (R_pos(i), i=1,6)
        end if 
        
        G(j) = G(j) - 2 * step_size
        R_neg = calc_r(V_mag, height, euler, angular_rates, G)
        
        if (trim_verbose) then
          write(io_unit, '(A)') '   Negative Finite-Difference Step'
          write(io_unit, '(A,6(1X,ES20.12))') '      G =', (G(i), i=1,6)
          write(io_unit, '(A,6(1X,ES20.12))') '      R =', (R_neg(i), i=1,6)
          write(io_unit,'(A)') ''
        end if

        do i = 1,6
          jac(i, j) = (R_pos(i) - R_neg(i)) / (2*step_size)

        end do 
        G(j) = G(j) + step_size
      end do 

      if (trim_verbose) then
        write(io_unit, '(A)') 'Jacobian Matrix ='
        do i = 1, size(jac,1)
            write(io_unit,'(*(1X,ES20.12))') (jac(i,j), j=1,size(jac,2))
        end do
      end if
      res = calc_r(V_mag, height, euler, angular_rates, G)
      res = -1* res

      ! CALCUALTE DELTA G AND ADD RELAXATION FACTOR
      call lu_solve(6, jac, res, delta_G)
      G = G + relaxation_factor * delta_G

      if (trim_verbose) then
        write(io_unit,*)
        write(io_unit,'(A,6(1X,ES20.12))') 'Delta G =', (delta_G(k), k=1,6)
      end if 

    end subroutine newtons_method

  !=========================
  ! Calculate Residual
    function calc_r(V_mag, height, euler, angular_rates, G) result(R)
      implicit none
      real, intent(in) :: V_mag, height, G(6), angular_rates(3)
      real, intent(inout) :: euler(3)
      real :: alpha, ca, cb, sa, sb
      real :: R(6), dummy_R(13), temp_state(13)

      temp_state = 0.0

      ! PULL OUT CONTROLS
      controls(1:4) = G(3:6)

      ! PULL OUT ALPHA
      alpha = G(1)
      ca = cos(alpha)
      sa = sin(alpha)     

      ! PULL OUT BETA
      if (sideslip_angle0 /= -999.0) then 
        cb = cos(sideslip_angle0)
        sb = sin(sideslip_angle0)
        euler(1) = G(2)
      else 
        cb = cos(G(2))
        sb = sin(G(2))
      end if

      ! CALCUALTE INITIAL STATES
      temp_state(1:3)   = V_mag * (/ca*cb, sb, sa*cb/) 
      temp_state(4:6)   = angular_rates
      temp_state(9)     = height
      temp_state(7:8)   = 0
      temp_state(10:13) = euler_to_quat(euler)


      ! CALCULATE RESIDUAL
      dummy_R = differential_equations(0.0, temp_state)
      R = dummy_R(1:6)

    end function calc_r

  !=========================
  ! Matrix Inverse
    function matrix_inv(A) result(A_inv)
      implicit none
      real, dimension(3,3), intent(in) :: A
      real, dimension(3,3) :: A_inv, A_adj
      real :: det

      ! CALCULATE THE DETERMINANT
      det = A(1,1) * A(2,2) * A(3,3) + A(1,2) * A(2,3) * A(3,1) + A(1,3) * A(2,1) * A(3,2) &
          - A(1,3) * A(2,2) * A(3,1) - A(1,2) * A(2,1) * A(3,3) - A(1,1) * A(2,3) * A(3,2)

      ! CALCULATE THE ADJUGATE MATRIX
      A_adj(1,1) = A(2,2) * A(3,3) - A(2,3) * A(3,2)
      A_adj(1,2) = A(1,3) * A(3,2) - A(1,2) * A(3,3)
      A_adj(1,3) = A(1,2) * A(2,3) - A(1,3) * A(2,2)

      A_adj(2,1) = A(2,3) * A(3,1) - A(2,1) * A(3,3)
      A_adj(2,2) = A(1,1) * A(3,3) - A(1,3) * A(3,1)
      A_adj(2,3) = A(1,3) * A(2,1) - A(1,1) * A(2,3)

      A_adj(3,1) = A(2,1) * A(3,2) - A(2,2) * A(3,1)
      A_adj(3,2) = A(1,2) * A(3,1) - A(1,1) * A(3,2)
      A_adj(3,3) = A(1,1) * A(2,2) - A(1,2) * A(2,1)

      ! CALCULATE THE INVERSE MATRIX
      A_inv = 1/det * A_adj
    end function matrix_inv

  !=========================
  ! CROSS PRODUCT
    function cross_product(vector_a, vector_b) result(vector_c)
      implicit none 
      real, intent(in) :: vector_a(3), vector_b(3)
      real :: vector_c(3)
      real :: a1, a2, a3, b1, b2, b3

      ! BREAK DOWN VECTOR 1 AND 2
      a1 = vector_a(1)
      a2 = vector_a(2)
      a3 = vector_a(3)

      b1 = vector_b(1)
      b2 = vector_b(2)
      b3 = vector_b(3)

      ! CALCULATE ORTHOGONAL VECTOR
      vector_c(1) = a2*b3 - a3*b2
      vector_c(2) = a3*b1 - a1*b3
      vector_c(3) = a1*b2 - a2*b1
    end function cross_product

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

  !=========================
  ! Init Subroutine
    subroutine init(filename)
      implicit none 
      real :: airspeed, alpha_deg, beta_deg
      real, allocatable :: eul(:)
      real :: alpha, beta, trim_state(6)
      character(100), intent(in) :: filename
      type(json_value), pointer :: j_connections, j_graphics, j_user_controls, j_test_controls

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
      call jsonx_get(j_main, 'vehicle.aerodynamics.sweep[deg]',                        sweep, 0.0)
      call jsonx_get(j_main, 'vehicle.aerodynamics.reference.area[ft^2]',              planform_area)
      call jsonx_get(j_main, 'vehicle.aerodynamics.reference.longitudinal_length[ft]', longitudinal_length)
      call jsonx_get(j_main, 'vehicle.aerodynamics.reference.lateral_length[ft]',      lateral_length)

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
      call jsonx_get(j_main, 'initial.airspeed[ft/s]',                          airspeed)
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

      ! CALCULATE INITIAL SPEED IN U, V, W DIRECTIONS
      initial_state(1) = airspeed * cos(alpha) * cos(beta)
      initial_state(2) = airspeed * sin(beta)
      initial_state(3) = airspeed * sin(alpha) * cos(beta)

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

      ! CALCULATE MASS AND INERTIA
      call mass_inertia() 

      ! CALCULATE THE SPECIFIED TRIM CONDITION IF GIVEN 
      if (sim_type == 'trim') then 
        trim_state = trim_algorithm(airspeed, initial_state(9), eul, tolerance, trim_type)
      end if 

      call jsonx_get(j_main, 'initial.x_pos_ft', initial_state(7))

      ! INITIALIZE CONNECTIONS
      call jsonx_get(j_main, 'connections', j_connections)
      call jsonx_get(j_connections, 'graphics', j_graphics)
      call graphics%init(j_graphics)

      ! TEST SECTION
      call jsonx_get(j_connections, 'user_controls', j_user_controls)
      call user_controls%init(j_user_controls)

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
        write(io_unit,*) controls_input

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
end module f16_m