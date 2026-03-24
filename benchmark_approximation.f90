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

    integer :: measurement_id_hardcoded_tvm, measurement_id_hardcoded_fortran

end module benchmarking

subroutine time_hardcoded_fortran( number_iterations, input )

    use benchmarking, only:     measurement_id_hardcoded_fortran
    use droplet_model, only:    initialize_model, estimate
    use measure, only:          REDUCTION_TYPE_MAX, &
                                end_phase, &
                                get_duration, &
                                start_phase

    implicit none

    integer, intent(in)              :: number_iterations
    real*4, dimension(7), intent(in) :: input
    real*4, dimension(2)             :: output

    real*4               :: dummy_sum
    real*4               :: duration
    integer              :: iteration_index

    ! Benchmark estimating a droplet with fixed parameters so we can report the
    ! output radius and temperature for verification purposes.
    dummy_sum = 0.0

    ! Initialize the weights and biases.
    call initialize_model()

    do iteration_index = 1, number_iterations
        call start_phase( measurement_id_hardcoded_fortran )
        call estimate( input, output )
        call end_phase( measurement_id_hardcoded_fortran )
        ! Ensure that the compiler keeps this loop by using the output.
        dummy_sum = dummy_sum + sum( output )
    end do

    duration = get_duration( measurement_id_hardcoded_fortran, REDUCTION_TYPE_MAX )

    write( *, "(I0,A,ES23.16,A,ES23.16,A)" ) number_iterations, " hardcoded fortran iterations took ", duration, " seconds at ", &
         duration / number_iterations, " seconds per iteration."

    write( *, * ) "Fortran Output: ", output
    write( *,* ) ""

end subroutine time_hardcoded_fortran


subroutine setup_measurements()

    use benchmarking, only:    measurement_id_hardcoded_fortran, &
                               measurement_id_hardcoded_tvm
    use measure, only:         create_phase, initialize_measurements

    implicit none

    call initialize_measurements()

    measurement_id_hardcoded_tvm     = create_phase( "hardcoded implementation times tvm" )
    measurement_id_hardcoded_fortran = create_phase( "hardcoded implementation times fortran" )

end subroutine setup_measurements

subroutine time_hardcoded_tvm( number_iterations, input )
    use benchmarking, only:     measurement_id_hardcoded_tvm
    use measure, only:          REDUCTION_TYPE_MAX, &
                                get_duration, &
                                start_phase, &
                                end_phase

    use inference_wrapper

    implicit none

    integer, intent(in)                     :: number_iterations
    real(c_float), dimension(7), intent(in) :: input
    integer                                 :: padded_number_iterations

    real*4  :: dummy_sum
    real*4  :: duration
    integer :: iteration_index, i, stat
    integer :: batch_count

    real(c_float), dimension(INPUT_COUNT, BATCH_SIZE)  :: tvm_inputs
    real(c_float), dimension(OUTPUT_COUNT, BATCH_SIZE) :: tvm_outputs

    batch_count              = ( number_iterations + BATCH_SIZE - 1 ) / BATCH_SIZE
    padded_number_iterations = batch_count * BATCH_SIZE

    if (padded_number_iterations /= number_iterations) then
        write(*,*) "Batch size ", BATCH_SIZE, " does not divide iteration count ", number_iterations, ". Padding to ", &
        padded_number_iterations, " iterations."
    endif

    write( *,* ) "Batch count: ", batch_count

    do i=1,BATCH_SIZE
       tvm_inputs(:, i) = input
    end do

    stat = init_inference()
    if (stat /= 0) then
        print *, "Failed to initialize TVM function."
        stop
    end if

    dummy_sum = 0.0
    do iteration_index = 1, batch_count
        call start_phase( measurement_id_hardcoded_tvm )
        call run_inference(tvm_inputs, tvm_outputs)
        call end_phase( measurement_id_hardcoded_tvm )

        dummy_sum = dummy_sum + sum( tvm_outputs )
    end do


   duration = get_duration( measurement_id_hardcoded_tvm, REDUCTION_TYPE_MAX )

   write( *, "(I0,A,ES23.16,A,ES23.16,A)" ) padded_number_iterations, " hardcoded tvm iterations took ", duration, " seconds at ", &
       duration / padded_number_iterations, " seconds per iteration."
   write(*,*) "TVM Outputs: "
   write(*,*) tvm_outputs(1, :)
   write(*,*) tvm_outputs(2, :)
   write(*,*) ""

end subroutine time_hardcoded_tvm

program benchmark_approximation
    use, intrinsic :: iso_c_binding

    implicit none

    integer, parameter :: NUMBER_HARDCODED_ITERATIONS = 64*102400
    real(c_float), parameter    :: DROPLET_INPUT(7) = [2.0e-7, 285.0, 10.0**-17.885, 285.0, 1.02, 1.00, 0.1]

    call setup_measurements()

    call time_hardcoded_fortran( NUMBER_HARDCODED_ITERATIONS, DROPLET_INPUT )

    call time_hardcoded_tvm( NUMBER_HARDCODED_ITERATIONS, DROPLET_INPUT )


end program benchmark_approximation

