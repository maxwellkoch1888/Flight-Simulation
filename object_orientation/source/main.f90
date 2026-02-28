program main 
    use sim_m
    use udp_m
    implicit none
    character(100) :: filename
    call udp_initialize()

    ! GET VALUES FROM SPECIFIED JSON
    call get_command_argument(1, filename)
    call init(filename)

    ! RUN SIMULATION
    call run()
    call udp_finalize()
    
end program main