program main
    use koch_m
    implicit none
    integer :: i, j, k, m, n
    real :: q1(4), q2(4), q3(4)
    real :: euler(3)
    real :: vec1(3), vec2(3)
    real :: start_time, end_time
    integer :: io_unit

    open(newunit=io_unit, file='1.8.9_output.txt', status='replace', action='write')

    ! ----------------- quat_mult test --------------------
    write(io_unit,*) 'quat_mult test'
    call cpu_time(start_time)
    q1 = [0.5, 0.5, 0.5, 0.5]
    q3 = q1
    do i = 1, 10
        do j = 1, 100000
            q2 = q3
            do k = 1, 100
                q3 = quat_mult(q1, q2)
            end do
        end do
        write(io_unit, '(I5, 4ES25.11)') i, q3
    end do
    call cpu_time(end_time)
    print *, 'quat_mult time total [sec]:              ', end_time - start_time

    ! ----------------- quat_base_to_dependent test --------------------
    write(io_unit,*) 'quat_base_to_dependent test'
    call cpu_time(start_time)
    vec1 = [0.0, 0.0, 1.0]
    do i = 1, 101, 4
        do j = 1, 101
            do k = 1, 101
                do m = 1, 101
                    q1 = [0.01*i, 0.01*j, 0.01*k, 0.01*m]
                    call quat_norm(q1)
                    vec2 = quat_base_to_dependent(vec1, q1)
                end do
            end do
        end do
        write(io_unit, '(I5, 3ES25.11)') i, vec2
    end do
    call cpu_time(end_time)
    print *, 'quat_base_to_dependent time total [sec]: ', end_time - start_time

    ! ----------------- quat_dependent_to_base test --------------------
    write(io_unit,*) 'quat_dependent_to_base test'
    call cpu_time(start_time)
    vec1 = [0.0, 0.0, 1.0]
    do i = 1, 101, 4
        do j = 1, 101
            do k = 1, 101
                do m = 1, 101
                    q1 = [0.01*i, 0.01*j, 0.01*k, 0.01*m]
                    call quat_norm(q1)
                    vec2 = quat_dependent_to_base(vec1, q1)
                end do
            end do
        end do
        write(io_unit, '(I5, 3ES25.11)') i, vec2
    end do
    call cpu_time(end_time)
    print *, 'quat_dependent_to_base time total [sec]: ', end_time - start_time

    ! ----------------- euler_to_quat test --------------------
    write(io_unit,*) 'euler_to_quat test'
    call cpu_time(start_time)
    do i = -180, 180, 2
        euler(1) = real(i) * PI / 180.0
        do j = -90, 90, 2
            euler(2) = real(j) * PI / 180.0
            do k = 0, 360, 2
                euler(3) = real(k) * PI / 180.0
                do m = 1, 10
                    q1 = euler_to_quat(euler)
                end do
            end do
        end do
        write(io_unit, '(I5, 4ES25.11)') i, q1
    end do
    call cpu_time(end_time)
    print *, 'euler_to_quat time total [sec]:          ', end_time - start_time

    ! ----------------- quat_to_euler test --------------------
    write(io_unit,*) 'quat_to_euler test'
    call cpu_time(start_time)
    do i = 0, 100, 4
        do j = 0, 100, 4
            do k = 0, 100, 2
                do m = 0, 100, 2
                    q1 = [0.01*i, 0.01*j, 0.01*k, 0.01*m]
                    call quat_norm(q1)
                    do n = 1, 10
                        euler = quat_to_euler(q1)
                    end do
                end do
            end do
        end do
        write(io_unit, '(I5, 3ES25.11)') i, euler * 180.0 / PI
    end do
    call cpu_time(end_time)
    print *, 'quat_to_euler time total [sec]:          ', end_time - start_time

    close(io_unit)
end program main
