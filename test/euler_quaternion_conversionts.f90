program euler_quaternion_conversions
    ! CONVERTS EULER ANGLES IN DEGREES TO QUATERNIONS
    real(8), dimension(3) :: euler_angle, rebuilt_euler_angle
    real(8), dimension(4) :: quaternion

    ! EULER ANGLES (BANK, PITCH, AZIMUTH)
    euler_angle = (/40, 20, 70/)

    quaternion = eul_to_quat(euler_angle)
    WRITE(*,*) quaternion

    rebuilt_euler_angle = quat_to_eul(quaternion)
    WRITE(*,*) rebuilt_euler_angle

    contains 
    
    function eul_to_quat(euler_angles) result(quaternion)
        real(8), dimension(3), intent(in) :: euler_angles
        real(8), dimension(4) :: quaternion
        real(8) :: pitch_angle, bank_angle, azimuth_angle
        real(8) :: e0, ex, ey, ez
        real(8) :: pi 
        pi = 4.D0 * DATAN(1.D0)

        ! EXTRACT THE PITCH, BANK, AND AZIMUTH ANGLES
        pitch_angle = euler_angles(2) * pi / 180
        bank_angle = euler_angles(1) * pi / 180
        azimuth_angle = euler_angles(3) * pi / 180

        ! CALCULATE e0, ex, ey, ez
        e0 = COS(bank_angle/2) * COS(pitch_angle/2) * &
             COS(azimuth_angle/2) + SIN(bank_angle/2) * &
             SIN(pitch_angle/2) * SIN(azimuth_angle/2)
        ex = SIN(bank_angle/2) * COS(pitch_angle/2) * &
             COS(azimuth_angle/2) - COS(bank_angle/2) * &
             SIN(pitch_angle/2) * SIN(azimuth_angle/2)
        ey = COS(bank_angle/2) * SIN(pitch_angle/2) * &
             COS(azimuth_angle/2) + SIN(bank_angle/2) * &
             COS(pitch_angle/2) * SIN(azimuth_angle/2)
        ez = COS(bank_angle/2) * COS(pitch_angle/2) * &
             SIN(azimuth_angle/2) - SIN(bank_angle/2) * &
             SIN(pitch_angle/2) * COS(azimuth_angle/2)

        ! BUILD THE QUATERNION
        quaternion = (/e0, ex, ey, ez/)
    end function eul_to_quat 

    function quat_to_eul(quaternion) result(euler_angles)
        real(8), dimension(4), intent(in) :: quaternion
        real(8), dimension(3) :: euler_angles
        real(8) :: pitch_angle, bank_angle, azimuth_angle
        real(8) :: e0, ex, ey, ez
        real(8) :: pi 
        pi = 4.D0 * DATAN(1.D0)

        ! EXTRACT e0, ex, ey, AND ez FROM THE QUATERNION
        e0 = quaternion(1)
        ex = quaternion(2)
        ey = quaternion(3)
        ez = quaternion(4)

        ! CALCULATE THE EQUIVALENT EULER ANGLES IN DEGREES
        bank_angle = ATAN2(2*(e0*ex + ey*ez), (e0**2 + ez**2 - ex**2 - ey**2)) * 180 / pi
        pitch_angle = ASIN(2*(e0*ey - ex*ez)) * 180 / pi
        azimuth_angle = ATAN2(2*(e0*ez + ex*ey), (e0**2 + ex**2 - ey**2 - ez**2)) * 180 / pi

        ! RETURN THE EULER ANGLE
        euler_angles = (/bank_angle, pitch_angle, azimuth_angle/)
    end function quat_to_eul
end program euler_quaternion_conversions