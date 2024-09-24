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

    real*4               :: dummy_sum
    real*4               :: duration
    integer              :: iteration_index

    ! Benchmark estimating a droplet with fixed parameters so we can report the
    ! output radius and temperature for verification purposes.
    input     = [3.59981e-5, 290.672, 2.06007e-12, 309.598, 1.08544, 1.18962, 1.0]
    dummy_sum = 0.0

    ! Initialize the weights and biases.
    call initialize_model()

    do iteration_index = 1, number_iterations
        call start_phase( measurement_id_hardcoded )
        call estimate( input, output )
        call end_phase( measurement_id_hardcoded )

        ! Ensure that the compiler keeps this loop by using the output.
        dummy_sum = dummy_sum + sum( output )
    end do

    duration = get_duration( measurement_id_hardcoded, REDUCTION_TYPE_MAX )

    write( *, "(I0,A,ES23.16,A,ES23.16,A)" ) number_iterations, " hardcoded iterations took ", duration, " seconds at ", &
         duration / number_iterations, " seconds per iteration."

    write( *, * ) output

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
