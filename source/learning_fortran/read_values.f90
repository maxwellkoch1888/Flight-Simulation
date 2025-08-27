program read_values
    implicit none
    real :: x, y

    print *, "Please enter two numbers."
    read(*,*) x, y

    print *, "The sum and product of the numbers are ", x+y, x*y

end program read_values