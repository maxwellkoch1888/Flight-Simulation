program quaternion_base_to_dependent
    ! Comptes a vector transformation using the euler rodriguez quaternion
    ! Accepts an arroay of length three containing the three components of a vector in the base coordinate system
    ! Also accepts an array of length four containing the four components of the quaternion that represent the orientation

    contains

    function quat_base_to_dependent(quaternion, base_vector) result(dependent_coordinates_vector)
        implicit none
        real(8), dimension(4), intent(in) :: quaternion
        real(8), dimension(3), intent(in) :: base_vector
        real(8), dimension(3) :: dependent_coordinates_vector
        real(8), dimension(4) :: vector, quat_conj
        
        ! BUILD VARIABLES FROM EQUATION 1.5.4 
        vector = (/0.0d0, base_vector(1), base_vector(2), base_vector(3)/)
        quat_conj = (/quaternion(1), -quaternion(2), -quaternion(3), -quaternion(4)/)

        ! FIRST ROTATION FROM 1.5.4
        temp_quat = quat_mult(vector, quaternion)

        ! SECOND ROTATION FROM 1.5.4
        dependent_coordinates_vector = quat_mult(quat_conj, temp_quat)

    end function quat_base_to_dependent

    function quat_dependent_to_base(quaternion, vector) result(dbase_coordinates_vector)
        implicit none
        real(8), dimension(4), intent(in) :: quaternion, vector
        real(8), dimension(3), intent(in) :: base_vector
        real(8), dimension(3) :: dependent_coordinates_vector


end program quaternion_base_to_dependent