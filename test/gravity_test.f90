program gravity_test
    use koch_m
    implicit none
    real :: gravity
    integer :: i

    do i = 0, 100000, 5000
        gravity = gravity_at_altitude_SI(real(i))
        write(*,'(I8, 2X, F16.12)') i, gravity
    end do
    
    do i = 0, 200000, 10000
        gravity = gravity_at_altitude_imperial(real(i))
        write(*,'(I8, 2X, F15.11)') i, gravity
    end do

    ! COMPILE STEPS
    ! 1: IN Flight_Simulation RUN
    ! gfortran -fdefault-real-8 -Wall chapter_1/koch.f90 test/gravity_test.f90 -o test/gravity_test.exe
    ! .\gravity_test.exe
end program gravity_test