program main 
    use sim_m
    implicit none
    character(100) :: filename

    ! GET VALUES FROM SPECIFIED JSON
    call get_command_argument(1, filename)
    call init(filename)

    ! RUN SIMULATION
    call run()
    
end program main