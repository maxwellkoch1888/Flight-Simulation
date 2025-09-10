program atmosphere_test 
    use koch_m
    implicit none
    real :: geometric_altitude_m, geopotential_altitude_m
    real :: temp_k, pressure_N_per_m2, density_kg_per_m3, sos_m_per_sec
    real, dimension(6) :: res
    integer :: i, io_unit

    ! OPEN A OUTPUT FILE
    open(newunit=io_unit, file='3.13.4_output.txt', status='replace', action='write')

    ! BUILD THE HEADER
    write(io_unit,'(8X,A21,5X,A24,2X,A14,12X,A15,A26,11X,A19)') &
        "Geometric_Altitude[m]", "Geopotential_Altitude[m]", "Temperature[K]", &
        "Pressure[N/m^2]", "Density[kg/m^3]","Speed_of_Sound[m/s]"
    
    ! BUILD THE RESULTS 
    do i = 0, 90000, 5000
        geometric_altitude_m = real(i)

        call std_atm_SI(geometric_altitude_m, geopotential_altitude_m, & 
            temp_k, pressure_N_per_m2, density_kg_per_m3, sos_m_per_sec)

        res = (/geometric_altitude_m, geopotential_altitude_m, temp_k, pressure_N_per_m2, density_kg_per_m3, sos_m_per_sec/)

        write(io_unit,'(6ES26.12)') res

    end do

end program atmosphere_test