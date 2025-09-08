program main
    use koch_m
    implicit none
    integer :: i,j
    real :: Z,T,P,rho,a,g
    real :: start_time, end_time
    integer :: io_unit

    open(newunit=io_unit, file='3.13.7_output.txt', status='replace', action='write')

    write(io_unit,*) 'SI atmospheric test'
    call cpu_time(start_time)
    do i=0,90000,100
        do j = 0,10000
            call std_atm_SI(real(i),Z,T,P,rho,a)
            g = gravity_SI(real(i))
        end do
        write(io_unit,'(7ES25.11)') real(i),Z,T,P,rho,a,g
    end do
    call cpu_time(end_time)
    print *, 'SI time total [sec]:          ', end_time - start_time

    write(io_unit,*) 'English atmospheric test'
    call cpu_time(start_time)
    do i=0,200000,200
        do j=1,10000
            call std_atm_English(real(i),Z,T,P,rho,a)
            g = gravity_English(real(i))
        end do
        write(io_unit,'(7ES25.11)') real(i),Z,T,P,rho,a,g
    end do
    call cpu_time(end_time)
    print *, 'English time total [sec]:     ', end_time - start_time

    close(io_unit)
end program main
