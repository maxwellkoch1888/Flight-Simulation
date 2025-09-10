program rk4_test
    use koch_m
    implicit none
    real:: t0, tf, y0, delta_t
    real, dimension(2) :: y0_vector
    t0 = 0.0 
    tf = 0.25 
    y0 = 0.0
    delta_t = 0.025

    ! TEST SCALAR FUNCTION
    call test_main(t0, tf, y0, delta_t)

    ! RESET STARTING VARIABLES
    t0 = 0.0 
    tf = 0.25 
    y0 = 0.0
    delta_t = 0.025
    

    ! TEST VECTOR_VERSION
    call test_main_vector(t0, tf, y0_vector, delta_t)

end program rk4_test