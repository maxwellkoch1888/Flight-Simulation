program variables
    implicit none

    integer :: amount
    real :: pi, e 
    complex :: frequency
    character :: initial
    logical :: is_okay

    amount = 10
    pi = 3.1415927
    e = 2.71
    frequency = (1.0, -0.5)
    initial = 'A'
    is_okay = .false.

    print *, 'The value of amount is ', amount
    print *, 'The value of pi is ', pi
    print *, 'The value of frequency is ', frequency
    print *, 'The value of initial is ', initial
    print *, 'The value of is_okay is ', is_okay

end program variables