program atmosphere_test_imperial 
    use koch_m
    implicit none
    real :: geometric_altitude_ft, geopotential_altitude_ft
    real :: temp_R, pressure_lbf_per_ft2, density_slugs_per_ft3, sos_ft_per_sec
    real, dimension(6) :: res
    integer :: i, io_unit

    ! OPEN AN OUTPUT FILE
    open(newunit=io_unit, file='3.13.6_output.txt', status='replace', action='write')

    ! BUILD THE HEADER
    write(io_unit,'(8X,A22,4X,A25,1X,A14,12X,A18,1X,A26,7X,A20)') &
        "Geometric_Altitude[ft]", "Geopotential_Altitude[ft]", "Temperature[R]", &
        "Pressure[lbf/ft^2]", "Density[slugs/ft^3]","Speed_of_Sound[ft/s]"
    
    ! BUILD THE RESULTS 
    do i = 0, 200000, 10000
        geometric_altitude_ft = real(i)

        call std_atm_English(geometric_altitude_ft, geopotential_altitude_ft, & 
        temp_R, pressure_lbf_per_ft2, density_slugs_per_ft3, sos_ft_per_sec)
    
        res = (/geometric_altitude_ft, geopotential_altitude_ft, temp_R, &
        pressure_lbf_per_ft2, density_slugs_per_ft3, sos_ft_per_sec/)

        write(io_unit,'(6ES26.12)') res

    end do

end program atmosphere_test_imperial