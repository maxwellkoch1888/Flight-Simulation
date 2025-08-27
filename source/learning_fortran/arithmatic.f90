program arithmetic
    implicit none

    real :: pi, radius, height, cylinder_area, volume

    pi = 3.1415927

    print *, "Enter cylinder radius."
    read(*,*) radius

    print*, "Enter cylinder height."
    read(*,*) height

    cylinder_area = pi * radius **2 
    volume = cylinder_area * height

    print *, 'Cylinder Volume:', volume

end program arithmetic