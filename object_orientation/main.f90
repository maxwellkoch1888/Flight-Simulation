program main 
    use sim_m
    ! use udp_m, only: udp_initialize, udp_finalize 
    use udp_m 
    implicit none
    character(100) :: filename

    ! INITIALIZE UDP CONNECTION
    call udp_initialize()

    ! GET VALUES FROM SPECIFIED JSON
    call get_command_argument(1, filename)
    call init(filename)

    ! RUN SIMULATION
    call run()

    ! END UDP CONNECTION
    call udp_finalize()
    
end program main