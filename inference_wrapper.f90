module inference_wrapper
    use, intrinsic :: iso_c_binding
    implicit none

    integer(c_int), parameter :: INPUT_COUNT  = 7
    integer(c_int), parameter :: OUTPUT_COUNT = 2
    integer(c_int), parameter :: BATCH_SIZE   = 4096 ! This must match the IR model batch size 

    interface
        ! Setup function
        integer(c_int) function raw_init_inference( input_count, output_count, batch_size ) bind(c, name="init_inference")
            import :: c_int
            integer(c_int), value, intent(in) :: input_count
            integer(c_int), value, intent(in) :: output_count
            integer(c_int), value, intent(in) :: batch_size
        end function raw_init_inference

        ! Hot loop function
        subroutine run_inference(in_data, out_data) bind(c, name="run_inference")
            import :: INPUT_COUNT, OUTPUT_COUNT, BATCH_SIZE
            import :: c_float, c_int
            real(c_float), intent(in)  :: in_data(INPUT_COUNT, BATCH_SIZE)
            real(c_float), intent(out) :: out_data(OUTPUT_COUNT, BATCH_SIZE)
        end subroutine run_inference

    end interface

CONTAINS
    ! Wrapper function to pass model shape constants
    integer(c_int) function init_inference() result(tvm_status)
        import :: c_int
        write(*,*) INPUT_COUNT
        tvm_status = raw_init_inference( INPUT_COUNT, OUTPUT_COUNT, BATCH_SIZE )
    end function init_inference

end module inference_wrapper
