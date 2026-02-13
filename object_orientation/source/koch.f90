module koch_m
    implicit none
    REAL, PARAMETER :: PI = 3.14159265358979323846264338327950288419716939937510
    real, parameter :: tol = 1.0e-12

    contains
!=========================
! Multiply two quaternions
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
!=========================
! Convert a quaternion from a base coordinate to a dependent coordinate
    function quat_base_to_dependent(base_vec, quat) result(dependent_vec)
        implicit none
        real, dimension(4), intent(in) :: quat
        real, dimension(3), intent(in) :: base_vec
        real, dimension(3) :: dependent_vec
        real, dimension(4) :: vec_quat, quat_conj, temp_quat, quat_sol
        
        ! BUILD VARIABLES FROM EQUATION 1.5.4 
        vec_quat = (/0.0, base_vec(1), base_vec(2), base_vec(3)/)
        quat_conj = (/quat(1), -quat(2), -quat(3), -quat(4)/)

        temp_quat = quat_mult(quat, vec_quat)
        quat_sol  = quat_mult(temp_quat, quat_conj)

        dependent_vec = quat_sol(2:4)

    end function quat_base_to_dependent
!=========================
! Convert a quaternion from a dependent coordinate to a base coordinate
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
!=========================
! Normalize a quaternion
    subroutine quat_norm (quat)
        implicit none
        real, dimension(4), intent(inout) :: quat
        real :: quat_mag

        ! FIND QUATERNION MAGNITUDE
        quat_mag = sqrt(sum(quat**2))

        ! NORMALIZE EACH COMPONENT OF THE QUATERNION
        quat(1) = quat(1) / quat_mag
        quat(2) = quat(2) / quat_mag
        quat(3) = quat(3) / quat_mag
        quat(4) = quat(4) / quat_mag
        
    end subroutine quat_norm
!=========================
! Convert orientation angles to a quaternion
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
!=========================
! Convert a quaternion to orientation angles
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
            if (azimuth_angle < 0.0) azimuth_angle = azimuth_angle + 2*pi  

        end if

        ! RETURN THE EULER ANGLE
        euler_angles = (/bank_angle, pitch_angle, azimuth_angle/)
    end function quat_to_euler
!=========================
! Calculate gravity in SI units
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
!=========================
! Calculate gravity in english units
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

!=========================
! Compute standard atm in SI units
    subroutine std_atm_SI( &
        geometric_altitude_m, geopotential_altitude_m, &
        temp_k, pressure_N_per_m2, density_kg_per_m3, &
        dyn_viscosity_pa_sec, sos_m_per_sec)

        implicit none
        real, intent(in)  :: geometric_altitude_m
        real, intent(out) :: geopotential_altitude_m, temp_k, pressure_N_per_m2, &
                            density_kg_per_m3, dyn_viscosity_pa_sec, sos_m_per_sec

        real, dimension(8) :: ti = (/288.150, 216.650, 216.650, 228.650, &
                                     270.650, 270.650, 252.650, 180.650/)
        real, dimension(8) :: ti_prime = (/-0.0065, 0.0, 0.001, 0.0028, &
                                        0.0, -0.0020, -0.004, 0.0/)
        real, dimension(8) :: p_i = (/1.0132500000000000E+05, 2.2632031822221168E+04, &
                                      5.4748735282708267E+03, 8.6801476908672271E+02, &
                                      1.1090558898922531E+02, 5.9000524278924367E+01, &
                                      1.8209924905017658E+01, 1.0377004548920223E+00/)
        real, dimension(9) :: zi = (/0.0, 11000.0, 20000.0, 32000.0, &
                                     47000.0, 52000.0, 61000.0, 79000.0, 90000.0/)
        real:: Rez = 6356766.0, R = 287.0528, gssl = 9.80665, gamma = 1.4, Ts = 273.15, mu_s = 1.716E-05, Ks = 110.4        
        integer :: i

        ! CALCULATE THE GEOPOTENTIAL ALTITUDE IN km
        geopotential_altitude_m = Rez * geometric_altitude_m / (Rez + geometric_altitude_m)

        ! CALCULATE THE TEMPERATURE
        i = 1
        do while (zi(i) <= geopotential_altitude_m)
            i = i + 1
        end do 
        i = i - 1

        temp_k = ti(i) + ti_prime(i) * (geopotential_altitude_m - zi(i))

        ! COMPUTE THE PRESSURE
        if(abs(ti_prime(i)) <= tol) then
            pressure_N_per_m2 = p_i(i) * exp(-(gssl * (geopotential_altitude_m - zi(i))) / (R * ti(i)))

        else
            pressure_N_per_m2 = p_i(i) * ((ti(i) + ti_prime(i) *  &
            (geopotential_altitude_m - zi(i))) / (ti(i))) ** (-gssl / (R * ti_prime(i)))
        end if

        ! CALCULATE THE DENSITY 3.2.8
        density_kg_per_m3 = pressure_N_per_m2 / (R * temp_k)

        ! COMPUTE VISCOSITY
        dyn_viscosity_pa_sec = mu_s * ( (Ts + Ks) / (temp_k + Ks)) * (temp_k / Ts)**1.5

        ! CALCULATE THE SPEED OF SOUND 3.2.9
        sos_m_per_sec = (gamma * R * temp_k) ** 0.5
    end subroutine std_atm_SI

!=========================
! Compute standard atm in engligh units
    subroutine std_atm_English(&
        geometric_altitude_ft, geopotential_altitude_ft,     & 
        temp_R, pressure_lbf_per_ft2, density_slugs_per_ft3, & 
        dyn_viscosity_slug_per_ft_sec, sos_ft_per_sec)
        
        implicit none
        real, intent(in) :: geometric_altitude_ft
        real, intent(out) :: geopotential_altitude_ft, temp_R, pressure_lbf_per_ft2, &
                            density_slugs_per_ft3, sos_ft_per_sec, dyn_viscosity_slug_per_ft_sec
        real :: geometric_altitude_m, geopotential_altitude_m, temp_k, pressure_N_per_m2, &
                density_kg_per_m3, sos_m_per_sec, dyn_viscosity_pa_sec

        ! CONVERT THE ALTITUDE TO FT
        geometric_altitude_m = geometric_altitude_ft * 0.3048

        ! USE THE SI VERSION OF ATMOSPHERE PROPERTIES
        call std_atm_SI(geometric_altitude_m, geopotential_altitude_m, & 
        temp_k, pressure_N_per_m2, density_kg_per_m3, dyn_viscosity_pa_sec, sos_m_per_sec)

        ! CONVERT THE UNITS BACK TO IMPERIAL
        geopotential_altitude_ft = geopotential_altitude_m / 0.3048
        temp_R = (temp_k * 1.8) 
        pressure_lbf_per_ft2 = pressure_N_per_m2 / 47.880258
        density_slugs_per_ft3 = density_kg_per_m3 / 515.379
        dyn_viscosity_slug_per_ft_sec = dyn_viscosity_pa_sec / 47.880258
        sos_ft_per_sec = sos_m_per_sec / 0.3048
    
    end subroutine std_atm_English

!=========================
! Matrix Inverse
    function matrix_inv(A) result(A_inv)
        implicit none
        real, dimension(3,3), intent(in) :: A
        real, dimension(3,3) :: A_inv
        real :: det
        real :: a11, a12, a13, a21, a22, a23, a31, a32, a33

        ! UNPACK ELEMENTS
        a11 = A(1,1); a12 = A(1,2); a13 = A(1,3)
        a21 = A(2,1); a22 = A(2,2); a23 = A(2,3)
        a31 = A(3,1); a32 = A(3,2); a33 = A(3,3)

        ! CALCULATE DETERMINANT
        det = a11*(a22*a33 - a23*a32) - a12*(a21*a33 - a23*a31) + a13*(a21*a32 - a22*a31)

        ! CALCULATE INVERSE ELEMENT-WISE
        A_inv(1,1) =  (a22*a33 - a23*a32)/det
        A_inv(1,2) = -(a12*a33 - a13*a32)/det
        A_inv(1,3) =  (a12*a23 - a13*a22)/det

        A_inv(2,1) = -(a21*a33 - a23*a31)/det
        A_inv(2,2) =  (a11*a33 - a13*a31)/det
        A_inv(2,3) = -(a11*a23 - a13*a21)/det

        A_inv(3,1) =  (a21*a32 - a22*a31)/det
        A_inv(3,2) = -(a11*a32 - a12*a31)/det
        A_inv(3,3) =  (a11*a22 - a12*a21)/det

    end function matrix_inv


!=========================
! Cross product
    function cross_product(a, b) result(c)
        implicit none
        real, intent(in) :: a(3), b(3)
        real :: c(3)
        ! DIRECTLY COMPUTE CROSS PRODUCT
        c(1) = a(2)*b(3) - a(3)*b(2)
        c(2) = a(3)*b(1) - a(1)*b(3)
        c(3) = a(1)*b(2) - a(2)*b(1)
    end function cross_product 
end module koch_m