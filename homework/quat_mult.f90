program quat_mult
    ! WRITE A FUNCTION NAMED quat_mult THAT ACCEPTS TWO ARRAYS OF LENGTH FOUR.
    ! REPRESENTS TWO QUATERNIONS AND RETURNS THE QUATERNION PRODUCT. 

    contains
    
    function quaternion_product(quatA, quatB) result(quatC)
        real, dimension(4), intent(in) :: quatA, quatB
        real, dimension(4) :: quatC
        real :: A0, AX, AY, AZ
        real :: B0, BX, BY, BZ
 
        ! EXTRACT THE VALUES FROM THE TWO QUATERNIONS
        A0 = quatA(1)
        AX = quatA(2)
        AY = quatA(3)
        AZ = quatA(4)

        B0 = quatB(1)
        BX = quatB(2)
        BY = quatB(3)
        BZ = quatB(4)

        ! BUILD THE 0, X, Y, AND Z COMPONENTS OF THE NEW QUATERNION
        C0 = 2
        CX = 1
        CY = 3
        CZ = 3


        
            
