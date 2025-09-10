program rk4_test
    use koch_m
    implicit none
    real:: t0, tf, y0, delta_t
    t0 = 0.0 
    tf = 0.25 
    y0 = 0.0
    delta_t = 0.025

    call test_main(t0, tf, y0, delta_t)

end program rk4_test