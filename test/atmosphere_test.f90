program atmosphere_test 
    use koch_m
    implicit none
    real :: geometric_altitude_m, geopotential_altitude_m
    real :: temp_k, pressure_N_per_m2, density_kg_per_m2, sos_m_per_sec

    geometric_altitude_m = 65000.0
    print*, "Building atmosphere"

    call atmospheric_properties_SI(&
        geometric_altitude_m, geopotential_altitude_m, & 
        temp_k, pressure_N_per_m2, density_kg_per_m2, sos_m_per_sec)

end program atmosphere_test