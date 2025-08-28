program quat_mult
    ! WRITE A FUNCTION NAMED quat_mult THAT ACCEPTS TWO ARRAYS OF LENGTH FOUR.
    ! REPRESENTS TWO QUATERNIONS AND RETURNS THE QUATERNION PRODUCT. 

    contains
    
    function quaternion_product(quat_a, quat_b) result(quat_c)
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
        c0 = a0*b0 - ax*bx - ay*by - az*bz
        cx = a0*bx + ax*b0 + ay*bz - az*by
        cy = a0*by - ax*bz + ay*b0 + az*bx
        cz = a0*bz + ax*by - ay*bx + az*b0

        ! BUILD THE RESULT OF THE QUATERNION PRODUCT
        quat_c = (/c0, cx, cy, cz/)
    end function quaternion_product

end program quat_mult



        
            
