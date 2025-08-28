program quaternion_multiplication
    ! WRITE A FUNCTION NAMED quat_mult THAT ACCEPTS TWO ARRAYS OF LENGTH FOUR.
    ! REPRESENTS TWO QUATERNIONS AND RETURNS THE QUATERNION PRODUCT. 
    real(8), dimension(4) :: quaternion, velocity, inertial_velocity, quaternion_conj, temporary_quat

    ! TEST USING EXAMPLE 11.7.1 MECHANICS OF FLIGHT
    quaternion = (/0.79212226029964539, 0.18231628330847405, 0.32686024176107303, 0.48214674108020061/)
    quaternion_conj = (/0.79212226029964539, -0.18231628330847405, -0.32686024176107303, -0.48214674108020061/)
    velocity = (/ 0.0, 825.96, 300.63, 476.87/)

    temporary_quat = quat_mult(velocity, quaternion_conj)
    WRITE(*,*) temporary_quat

    inertial_velocity = quat_mult(quaternion, temporary_quat)

    WRITE(*,*) inertial_velocity

    
    contains
    
    function quat_mult(quat_a, quat_b) result(quat_c)
        real(8), dimension(4), intent(in) :: quat_a, quat_b
        real(8), dimension(4) :: quat_c
        real(8) :: a0, ax, ay, az
        real(8) :: b0, bx, by, bz
        real(8) :: c0, cx, cy, cz
 
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
        c0 = a0*b0 - ax*bx - ay*by - az*bz
        cx = a0*bx + ax*b0 + ay*bz - az*by
        cy = a0*by - ax*bz + ay*b0 + az*bx
        cz = a0*bz + ax*by - ay*bx + az*b0

        ! BUILD THE RESULT OF THE QUATERNION PRODUCT
        quat_c = (/c0, cx, cy, cz/)
    end function quat_mult

end program quaternion_multiplication



        
            
