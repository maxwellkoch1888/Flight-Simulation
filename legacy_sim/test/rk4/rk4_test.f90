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
    
    subroutine test_main_vector(t0, tf, y0, delta_t)
        implicit none
        real, intent(in) :: delta_t, tf
        real, intent(inout) :: t0
        real, intent(inout), dimension(:) :: y0
        real, dimension(size(y0)) :: y1, k1, k2, k3, k4
        integer :: io_unit

        ! OPEN AN OUTPUT FILE
        open(newunit=io_unit, file='5.9.2_output.txt', status='replace', action='write')

        ! LOOP THE FUNCTION
        write(io_unit, '(A8, A18, A24)') 't', 'x', 'z'


        do while(t0 < tf)
            ! CALCULATE K VALUES FOR TABLE
            k1 = differential_equations_vector(t0, y0)
            k2 = differential_equations_vector(t0 + delta_t*0.5, y0 + k1 * delta_t*0.5)
            k3 = differential_equations_vector(t0 + delta_t*0.5, y0 + k2 * delta_t*0.5)
            k4 = differential_equations_vector(t0 + delta_t, y0 + k3 * delta_t)

            ! CALCULATE NEW Y VALUE
            y1 = runge_kutta_vector(t0, y0, delta_t)

            ! WRITE VALUES TO THE TABLE
            write(io_unit,'(1F10.3, 2F24.14)') t0, y0

            ! UPDATE VALUES FOR t0 AND y0
            t0 = t0 + delta_t
            y0 = y1

        end do
            
    end subroutine test_main_vector
    ! RESET STARTING VARIABLES
    t0 = 0.0 
    tf = 0.25 
    y0 = 0.0
    delta_t = 0.025
    

    ! TEST VECTOR_VERSION
    call test_main_vector(t0, tf, y0_vector, delta_t)

end program rk4_test