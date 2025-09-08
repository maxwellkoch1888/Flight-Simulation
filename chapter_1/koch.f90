module koch_m
    implicit none
    REAL, PARAMETER :: PI = 3.14159265358979323846264338327950288419716939937510
    real, parameter :: tol = 1.0e-12

    contains
! -------------------------------------------------------------------------
! PROBLEM 1.8.3
    function quat_mult(quat_a, quat_b) result(quat_c)
        implicit none
        real, dimension(4), intent(in) :: quat_a, quat_b
        real, dimension(4) :: quat_c
        real :: a0, ax, ay, az
        real :: b0, bx, by, bz
 
        ! EXTRACT THE VALUES FROM THE TWO QUATERNIONS
        a0 = quat_a(1)
        ax = quat_a(2)
        ay = quat_a(3)
        az = quat_a(4)

        b0 = quat_b(1)
        bx = quat_b(2)
        by = quat_b(3)
        bz = quat_b(4)

        ! BUILD THE 0, X, Y, AND Z COMPONENTS OF THE NEW QUATERNION
        quat_c(1) = a0*b0 - ax*bx - ay*by - az*bz
        quat_c(2) = a0*bx + ax*b0 + ay*bz - az*by
        quat_c(3) = a0*by - ax*bz + ay*b0 + az*bx
        quat_c(4) = a0*bz + ax*by - ay*bx + az*b0

    end function quat_mult
! -------------------------------------------------------------------------
! PROBLEM 1.8.4
    function quat_base_to_dependent(base_vec, quat) result(dependent_vec)
        implicit none
        real, dimension(4), intent(in) :: quat
        real, dimension(3), intent(in) :: base_vec
        real, dimension(3) :: dependent_vec
        real, dimension(4) :: vec_quat, quat_conj, temp_quat, quat_sol
        
        ! BUILD VARIABLES FROM EQUATION 1.5.4 
        vec_quat = (/0.0, base_vec(1), base_vec(2), base_vec(3)/)
        quat_conj = (/quat(1), -quat(2), -quat(3), -quat(4)/)

        ! FIRST ROTATION FROM 1.5.4
        temp_quat = quat_mult(vec_quat, quat)

        ! SECOND ROTATION FROM 1.5.4
        quat_sol = quat_mult(quat_conj, temp_quat)
        dependent_vec = (/quat_sol(2), quat_sol(3), quat_sol(4)/)

    end function quat_base_to_dependent
! -------------------------------------------------------------------------
! PROBLEM 1.8.5
    function quat_dependent_to_base(dependent_vec, quat) result(base_vec)
        implicit none
        real, dimension(4), intent(in) :: quat
        real, dimension(3), intent(in) :: dependent_vec
        real, dimension(3) :: base_vec
        real, dimension(4) :: vec_quat, quat_conj, temp_quat, quat_sol

        ! BUILD VARIABLES FROM EQUATION 1.5.4 
        vec_quat = (/0.0, dependent_vec(1), dependent_vec(2), dependent_vec(3)/)
        quat_conj = (/quat(1), -quat(2), -quat(3), -quat(4)/)
        
        ! FIRST ROTATION FROM 1.5.4
        temp_quat = quat_mult(vec_quat, quat_conj)

        ! SECOND ROTATION FROM 1.5.4
        quat_sol = quat_mult(quat, temp_quat)
        base_vec = (/quat_sol(2), quat_sol(3), quat_sol(4)/)
    end function quat_dependent_to_base
! -------------------------------------------------------------------------
! PROBLEM 1.8.6
    subroutine quat_norm (quat)
        implicit none
        real, dimension(4), intent(inout) :: quat
        real :: quat_mag

        ! FIND QUATERNION MAGNITUDE
        quat_mag = (quat(1)**2 + quat(2)**2 + quat(3)**2 + quat(4)**2)**0.5

        ! NORMALIZE EACH COMPONENT OF THE QUATERNION
        quat(1) = quat(1) / quat_mag
        quat(2) = quat(2) / quat_mag
        quat(3) = quat(3) / quat_mag
        quat(4) = quat(4) / quat_mag
        
    end subroutine quat_norm
! -------------------------------------------------------------------------
! PROBLEM 1.8.7
    function euler_to_quat(euler_angles) result(quat)
        implicit none
        real, dimension(3), intent(in) :: euler_angles
        real, dimension(4) :: quat
        real :: pitch_angle, bank_angle, azimuth_angle
        real :: e0, ex, ey, ez, cb, sb, cp, sp, ca, sa

        ! EXTRACT THE PITCH, BANK, AND AZIMUTH ANGLES
        pitch_angle = euler_angles(2)
        bank_angle = euler_angles(1)
        azimuth_angle = euler_angles(3)

        ! BUILD THE SIN AND COSINE VARIABLES
        cb = COS(bank_angle/2.0)
        sb = SIN(bank_angle/2.0)
        cp = COS(pitch_angle/2.0)
        sp = SIN(pitch_angle/2.0)
        ca = COS(azimuth_angle/2.0)
        sa = SIN(azimuth_angle/2.0)

        ! CALCULATE e0, ex, ey, ez
        e0 = cb*cp*ca + sb*sp*sa
        ex = sb*cp*ca - cb*sp*sa
        ey = cb*sp*ca + sb*cp*sa
        ez = cb*cp*sa - sb*sp*ca

        ! BUILD THE QUATERNION
        quat = (/e0, ex, ey, ez/)
    end function euler_to_quat 
! -------------------------------------------------------------------------
! PROBLEM 1.8.8
    function quat_to_euler(quat) result(euler_angles)
        implicit none
        real, dimension(4), intent(in) :: quat
        real, dimension(3) :: euler_angles
        real :: pitch_angle, bank_angle, azimuth_angle
        real :: e0, ex, ey, ez, res

        ! EXTRACT e0, ex, ey, AND ez FROM THE quat
        e0 = quat(1)
        ex = quat(2)
        ey = quat(3)
        ez = quat(4)

        ! CALCULATE RESULT FOR GIMBLE LOCK
        res = e0*ey - ex*ez

        ! CALCULATE THE EQUIVALENT EULER ANGLES
        if(abs(res - 0.5) < tol) then
            azimuth_angle = 0.0
            bank_angle = 2*ASIN(ex/COS(pi/4.0) + azimuth_angle)
            pitch_angle = pi/2.0
        else if(abs(res + 0.5) < tol) then
            azimuth_angle = 0.0
            bank_angle = 2.0*ASIN(ex/COS(pi/4.0) - azimuth_angle)
            pitch_angle = -pi/2
        else
            bank_angle = ATAN2(2.0*(e0*ex + ey*ez), (e0**2 + ez**2 - ex**2 - ey**2))
            pitch_angle = ASIN(2.0*(e0*ey - ex*ez))
            azimuth_angle = ATAN2(2.0*(e0*ez + ex*ey), (e0**2 + ex**2 - ey**2 - ez**2))
        end if

        ! RETURN THE EULER ANGLE
        euler_angles = (/bank_angle, pitch_angle, azimuth_angle/)
    end function quat_to_euler
! -------------------------------------------------------------------------
! PROBLEM 1.8.9
! MODIFY THE CODE TO MAKE IT RUN AS QUICKLY AS POSSIBLE
! -------------------------------------------------------------------------
!PROBLEM 3.13.1
! WRITE A FUNCTION TO COMPUTE THE GRAVITY AS A FUNCTION OF ALTITUDE USING EQ 3.2.1 IN SI UNITS
    function gravity_SI(geometric_altitude_m) result(gravity_m_per_sec2)
        implicit none
        real, intent(in) :: geometric_altitude_m
        real :: gravity_m_per_sec2, gssl, Rez
        
        ! DEFINE GRAVITY AT STANDARD SEA LEVEL AND RADIUS OF EARTH IN THE US
        gssl = 9.80665
        Rez = 6356766.0

        ! CALCULATE GRAVITY AT ALTITUDE
        gravity_m_per_sec2 = gssl * (Rez / (Rez + geometric_altitude_m)) ** 2
    end function gravity_SI
! -------------------------------------------------------------------------

! PROBLEM 3.13.2 WRITE A FUNCTION TO COMPUTE THE GRAVITY AS A FUNCTION OF GEOMETRIC ALTITUDE IN ENGLISH UNITS
    function gravity_English(geometric_altitude_ft) result(gravity_ft_per_sec2)
        implicit none
        real, intent(in) :: geometric_altitude_ft
        real :: gravity_ft_per_sec2, gssl, Rez 
        ! DEFINE GRAVITY AT STANDARD SEA LEVEL AND RADIUS OF EARTH IN THE US
        gssl = 9.80665 / 0.3048
        Rez = 6356766.0 / 0.3048

        ! CALCULATE GRAVITY AT ALTITUDE
        gravity_ft_per_sec2 = gssl * (Rez / (Rez + geometric_altitude_ft)) ** 2
    end function gravity_English
! -------------------------------------------------------------------------

! PROBLEM 3.13.3 WRITE A FUNCTION TO COMPUTE STANDARD ATMOSPHERIC PROPERTIES IN SI UNITS
    subroutine std_atm_SI(&
        geometric_altitude_m, geopotential_altitude_m, & 
        temp_k, pressure_N_per_m2, density_kg_per_m3, sos_m_per_sec)
        implicit none
        real, intent(in) :: geometric_altitude_m
        real, intent(out) :: geopotential_altitude_m, temp_k, pressure_N_per_m2, density_kg_per_m3, sos_m_per_sec
        real, dimension(8) :: ti, ti_prime
        real, dimension(8) :: p_i
        real, dimension(9) :: zi
        real :: Rez, R, gssl, power, gamma
        logical :: altitude_condition = .false.
        integer :: i, j

        ! DEFINE RADIUS OF EARTH IN US, GAS CONSTANT, INITIAL PRESSURE, AND THE TEMPERATURE VARIATION TABLE
        Rez = 6356766.0
        gssl = 9.80665
        R = 287.0528
        gamma = 1.4
        p_i(1) = 101325.0
        zi = (/0.0, 11000.0, 20000.0, 32000.0, 47000.0, 52000.0, 61000.0, 79000.0, 90000.0/)
        ti = (/288.150, 216.650, 216.650, 228.650, 270.650, 270.650, 252.650, 180.650/)
        ti_prime = (/-0.0065, 0.0, 0.001, 0.0028, 0.0, -0.0020, -0.004, 0.0/)

        ! CALCULATE THE GEOPOTENTIAL ALTITUDE IN km
        geopotential_altitude_m = Rez * geometric_altitude_m / (Rez + geometric_altitude_m)

        ! print*, "Geometric Altitude"
        ! print*, geometric_altitude_m
        ! print*, ""
        ! print*, "Geopotential Altitude" 
        ! print*, geopotential_altitude_m
        ! print*, ""

        ! CALCULATE THE TEMPERATURE
        i = 1
        do while (zi(i) <= geopotential_altitude_m)
            i = i + 1
        end do 
        i = i - 1

        temp_k = ti(i) + ti_prime(i) * (geopotential_altitude_m - zi(i))

        ! print*, "Base temp range"
        ! print*, ti(i)
        ! print*, ""
        ! print*, "Base T prime"
        ! print*, ti_prime(i)
        ! print*, "" 

        ! COMPUTE THE PRESSURE RANGES
        j = 1
        do while (j < 8)
            if (ti_prime(j) == 0.0) then
                power = -(gssl * (zi(j+1) - zi(j))) / (R * ti(j))
                p_i(j+1) = p_i(j) * exp(power)

            else
                power = -gssl / (R * ti_prime(j))
                p_i(j+1) = p_i(j) * ((ti(j) + ti_prime(j) * (zi(j+1) - zi(j))) / (ti(j))) ** (power)
            end if
            altitude_condition = geopotential_altitude_m <= zi(j)
            j = j + 1
        end do 

        if (ti_prime(i) == 0.0) then
            power = -(gssl * (geopotential_altitude_m - zi(i))) / (R * ti(i))
            pressure_N_per_m2 = p_i(i) * exp(power)

        else
            power = -gssl / (R * ti_prime(i))
            pressure_N_per_m2 = p_i(i) * ((ti(i) + ti_prime(i) * (geopotential_altitude_m - zi(i))) / (ti(i))) ** (power)
        end if

        ! print*, "Base pressure range"
        ! print*, p_i(i)
        ! print*, ""

        ! CALCULATE THE DENSITY 3.2.8
        density_kg_per_m3 = pressure_N_per_m2 / (R * temp_k)

        ! CALCULATE THE SPEED OF SOUND 3.2.9
        sos_m_per_sec = (gamma * R * temp_k) ** 0.5
    end subroutine std_atm_SI


! -------------------------------------------------------------------------

! PROBLEM 3.13.4 WRITE A FUNCTION TO COMPUTE STANDARD ATMOSPHERIC PROPERTIES IN IMPERIAL UNITS
    subroutine std_atm_English(&
        geometric_altitude_ft, geopotential_altitude_ft, & 
        temp_R, pressure_lbf_per_ft2, density_slugs_per_ft3, sos_ft_per_sec)
        
        implicit none
        real, intent(in) :: geometric_altitude_ft
        real, intent(out) :: geopotential_altitude_ft, temp_R, pressure_lbf_per_ft2, density_slugs_per_ft3, sos_ft_per_sec
        real :: geometric_altitude_m, geopotential_altitude_m, temp_k, pressure_N_per_m2, density_kg_per_m3, sos_m_per_sec

        ! CONVERT THE ALTITUDE TO FT
        geometric_altitude_m = geometric_altitude_ft * 0.3048

        ! USE THE SI VERSION OF ATMOSPHERE PROPERTIES
        call std_atm_SI(geometric_altitude_m, geopotential_altitude_m, & 
        temp_k, pressure_N_per_m2, density_kg_per_m3, sos_m_per_sec)

        ! CONVERT THE UNITS BACK TO IMPERIAL
        geopotential_altitude_ft = geopotential_altitude_m / 0.3048
        temp_R = (temp_k * 1.8) 
        pressure_lbf_per_ft2 = pressure_N_per_m2 / 47.880258
        density_slugs_per_ft3 = density_kg_per_m3 / 515.379
        sos_ft_per_sec = sos_m_per_sec / 0.3048
    
    end subroutine std_atm_English
        

end module koch_m