module derived_types
    use iso_c_binding
    use activation_functions

    implicit none

    abstract interface
        function func (z) result(output)
            REAL, intent(in) :: z(:,:)
            real :: output(size(z,1), size(z,2))
        end function func
    end interface
    
    TYPE linLayer
        REAL (c_float), ALLOCATABLE, DIMENSION(:,:) :: weights
        REAL (c_float), ALLOCATABLE, DIMENSION(:) :: biases
    ENDTYPE linLayer

    TYPE lstmLayer
        REAL (c_float), ALLOCATABLE, DIMENSION(:,:,:) :: whh
        REAL (c_float), ALLOCATABLE, DIMENSION(:,:,:) :: wih
        REAL (c_float), ALLOCATABLE, DIMENSION(:) :: bhh
        REAL (c_float), ALLOCATABLE, DIMENSION(:) :: bih
        REAL (c_float), ALLOCATABLE, DIMENSION(:,:,:) :: hid
        REAL (c_float), ALLOCATABLE, DIMENSION(:,:,:) :: cell
    ENDTYPE lstmLayer

    TYPE convLayer
        REAL (c_float), ALLOCATABLE, DIMENSION(:,:,:,:) :: weights
        REAL (c_float), ALLOCATABLE, DIMENSION(:) :: biases
        !==stride
    ENDTYPE convLayer

    TYPE maxpoolLayer
        INTEGER :: kernel_size
    ENDTYPE maxpoolLayer

    TYPE avgpoolLayer
        INTEGER :: kernel_size
    ENDTYPE avgpoolLayer

    TYPE addLayer
        REAL (c_float), ALLOCATABLE, DIMENSION(:,:,:,:) :: adder
    ENDTYPE addLayer

    TYPE reshapeLayer
        REAL (c_float), ALLOCATABLE, DIMENSION(:,:) :: reshape2d
        REAL (c_float), ALLOCATABLE, DIMENSION(:,:,:) :: reshape3d
        REAL (c_float), ALLOCATABLE, DIMENSION(:,:,:,:) :: reshape4d
    ENDTYPE

    
end module derived_types