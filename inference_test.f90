! MLP inference driver that evaluates a droplet's parameters and prints out the
! approximation produced by the droplet_model module.  This serves to verify
! that serialized models behave as expected before launching a full fledged
! NTLP run.

program inference_test

    use droplet_model, only:    initialize_model, estimate

    implicit none

    real*4, dimension(7) :: input
    real*4, dimension(2) :: output
    integer              :: number_arguments
    integer              :: argument_index
    character(len=100)   :: argument
    integer              :: ios

    ! Check if we have the correct number of command line arguments
    number_arguments = command_argument_count()

    if (number_arguments /= 7) then
        write( *, "(A)" ) "Error: Expected 7 command line arguments"
        write( *, "(A)" ) "Usage: ./inference_test <radius> <temperature> <salt_solute> <air_temperature> <rh> <air_density> <dt>"
        write( *, "(A)" ) "Example: ./inference_test 1.09402e-07 295.765 1.75657e-20 295.203 1.08863 1.28343 0.1"
        stop 1
    end if

    ! Read command line arguments and convert to real numbers
    do argument_index = 1, 7
        call get_command_argument( argument_index, argument )
        read( argument, *, iostat=ios ) input(argument_index)

        if (ios /= 0) then
            write( *, "(A,I1,A)" ) "Error: Could not parse argument ", argument_index, " as a number"
            write( *, "(A,A)" ) "Argument value: ", trim( argument )
            stop 1
        end if
    end do

    ! Initialize the weights and biases.
    call initialize_model()

    call estimate( input, output )

    ! Report our inputs and the estimate.
    write( *, "(A,ES15.7,A)" ) "Input Radius:       ", input(1), " m"
    write( *, "(A,F12.7,A)" )  "Input Temperature:   ", input(2), " K"
    write( *, "(A,ES15.7,A)" ) "Salt Solute:        ", input(3), " kg"
    write( *, "(A,F12.7,A)" )  "Air Temperature:     ", input(4), " K"
    write( *, "(A,F12.7,A)" )  "Relative Humidity:   ", input(5) * 100.0, " %"
    write( *, "(A,F12.7,A)" )  "Air Density:       ", input(6), " kg/m^3"
    write( *, "(A,F6.3,A)" )   "Integration Time:    ", input(7), " s"

    write( *, * )
    write( *, "(A,ES15.7,A)" ) "Output Radius:      ", output(1), " m"
    write( *, "(A,F12.7,A)" )  "Output Temperature:  ", output(2), " K"

end program inference_test
