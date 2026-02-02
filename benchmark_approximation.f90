! Benchmarking driver that measures and compares the execution time required for
! estimating a droplet's characteristics at some point in the future.  This also
! serves as a simple, manual tool for verifying that the Fortran implementations
! reproduce an expected value as a (very narrowly focused) sanity check.
!
! Currently this only benchmarks the implementation in the droplet_model module.
! Preliminary timings on Dr. Richter's Ivy Bridge and Epyc 7002 compute nodes
! result in the following observations:
!
!   1. Naive, looped implementation is significantly faster than Gauss-Newton
!      iterations with speedups of 98x and 34x on Ivy Bridge and Epyc nodes,
!      respectively.
!

module benchmarking

    integer :: measurement_id_hardcoded

end module benchmarking

subroutine time_hardcoded( number_iterations )

    use benchmarking, only:     measurement_id_hardcoded
    use droplet_model, only:    initialize_model, estimate
    use measure, only:          REDUCTION_TYPE_MAX, &
                                end_phase, &
                                get_duration, &
                                start_phase

    implicit none

    integer, intent(in)  :: number_iterations

    real*4, dimension(7) :: input
    real*4, dimension(2) :: output

    integer              :: number_arguments
    integer              :: argument_index
    character(len=100)   :: argument
    integer              :: ios

    real*4               :: dummy_sum
    real*4               :: duration
    integer              :: iteration_index

    ! Check if we have the correct number of command line arguments
    number_arguments = command_argument_count()

    if (number_arguments /= 0 .and. number_arguments /= 7) then
        write( *, "(A)" ) "Error: Expected either one or 7 command line arguments"
        write( *, "(A)" ) "Usage: ./benchmark_approximation.x [<radius> <temperature> <salt_solute> <air_temperature> <rh> <air_density> <dt>]"
        write( *, "(A)" ) "Example: ./benchmark_approximation.x 1.09402e-07 295.765 1.75657e-20 295.203 1.08863 1.28343 0.1"
        stop 1
    end if

    ! Read in the user's droplet parameters if provided, otherwise default to a
    ! random one.
    if (number_arguments == 7) then
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
    else
        input = [3.59981e-5, 290.672, 2.06007e-12, 309.598, 1.08544, 1.18962, 1.0]
    end if

    ! We accumulate the outputs to prevent the compiler from optimizing away
    ! our calls.
    dummy_sum = 0.0

    ! Initialize the weights and biases.
    call initialize_model()

    ! Report our inputs.
    write( *, "(A,ES15.7,A)" ) "Input Radius:       ", input(1), " m"
    write( *, "(A,F12.7,A)" )  "Input Temperature:   ", input(2), " K"
    write( *, "(A,ES15.7,A)" ) "Salt Solute:        ", input(3), " kg"
    write( *, "(A,F12.7,A)" )  "Air Temperature:     ", input(4), " K"
    write( *, "(A,F12.7,A)" )  "Relative Humidity:   ", input(5) * 100.0, " %"
    write( *, "(A,F12.7,A)" )  "Air Density:       ", input(6), " kg/m^3"
    write( *, "(A,F6.3,A)" )   "Integration Time:    ", input(7), " s"

    do iteration_index = 1, number_iterations
        call start_phase( measurement_id_hardcoded )
        call estimate( input, output )
        call end_phase( measurement_id_hardcoded )

        ! Ensure that the compiler keeps this loop by using the output.
        dummy_sum = dummy_sum + sum( output )
    end do

    duration = get_duration( measurement_id_hardcoded, REDUCTION_TYPE_MAX )

    write( *, * )
    write( *, "(A,ES15.7,A)" ) "Output Radius:      ", output(1), " m"
    write( *, "(A,F12.7,A)" )  "Output Temperature:  ", output(2), " K"

    write( *, * )
    write( *, "(I0,A,ES23.16,A,ES23.16,A)" ) number_iterations, " hardcoded iterations took ", duration, " seconds at ", &
         duration / number_iterations, " seconds per iteration."

end subroutine time_hardcoded

subroutine setup_measurements()

    use benchmarking, only:    measurement_id_hardcoded
    use measure, only:         create_phase, initialize_measurements

    implicit none

    call initialize_measurements()

    measurement_id_hardcoded = create_phase( "hardcoded implementation times" )

end subroutine setup_measurements

program benchmark_approximation

    implicit none

    integer, parameter :: NUMBER_HARDCODED_ITERATIONS = 10000000

    call setup_measurements()

    call time_hardcoded( NUMBER_HARDCODED_ITERATIONS )

end program benchmark_approximation
