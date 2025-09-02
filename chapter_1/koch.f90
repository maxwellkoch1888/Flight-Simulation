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

        ! BUILD A LOCAL FUNCTION TO CORRECT ROUNDING ERRORS FOR PI
        ! contains
        !     function clean_value(x) result(y)
        !         real, intent(in) :: x
        !         real :: y
        !         if(abs(x) < tol) then
        !             y = 0
        !         else if(abs(x - 1) < tol) then
        !             y = 1
        !         else if(abs(x+1) < tol) then
        !             y = -1
        !         else
        !             y = x 
        !         end if 
        !     end function clean_value
    end function euler_to_quat 
! -------------------------------------------------------------------------
! PROBLEM 1.8.8
    function quat_to_euler(quat) result(euler_angles)
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
            azimuth_angle = 0
            bank_angle = 2*ASIN(ex/COS(pi/4) + azimuth_angle)
            pitch_angle = pi/2
        else if(abs(res + 0.5) < tol) then
            azimuth_angle = 0
            bank_angle = 2*ASIN(ex/COS(pi/4) - azimuth_angle)
            pitch_angle = -pi/2
        else
            bank_angle = ATAN2(2*(e0*ex + ey*ez), (e0**2 + ez**2 - ex**2 - ey**2))
            pitch_angle = ASIN(2*(e0*ey - ex*ez))
            azimuth_angle = ATAN2(2*(e0*ez + ex*ey), (e0**2 + ex**2 - ey**2 - ez**2))
        end if

        ! RETURN THE EULER ANGLE
        euler_angles = (/bank_angle, pitch_angle, azimuth_angle/)
    end function quat_to_euler
! -------------------------------------------------------------------------
! PROBLEM 1.8.9
! MODIFY THE CODE TO MAKE IT RUN AS QUICKLY AS POSSIBLE
end module koch_m