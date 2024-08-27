module model
    !---------------
    ! adding any number of layers to our neural network
    !----------------


    ! ===============================================================
    ! USE filereader !<loading in weights, biases
    USE activation_functions !<getting activation functions
    USE model_layers
    USE reader
    USE iso_c_binding
    ! ===============================================================


    IMPLICIT NONE
    contains
    !===============================================THIS BLOCK COMES FROM FYPP OUTPUT "openNP.fpp"=======================================================

        SUBROUTINE use_model(        i0,         o0) bind(c,name="use_model")
            IMPLICIT NONE
            !INPUTS CORRESPONDING TO C
            REAL (c_float), INTENT(INOUT), DIMENSION(        1, 6) :: i0
            !===========================

            !INPUTS CORRESPONDING TO INTERMEDIARY PROCESSING
            !================================================

            !OUTPUTS CORRESPONDING TO C
            REAL (c_float), INTENT(OUT), DIMENSION(        1, 2) :: o0
            !===========================

            REAL (c_float), ALLOCATABLE, DIMENSION(:,:) :: Gemm_0

            REAL :: T1, T2

            ALLOCATE(Gemm_0(        1, 6))

            Gemm_0 = i0


            CALL CPU_TIME(T1)

            !========Gemm Layer============
            CALL linear_layer(Gemm_0, linLayers(1),0)

            Gemm_0 = relu2d(Gemm_0)

            !========Gemm Layer============
            CALL linear_layer(Gemm_0, linLayers(2),0)

            Gemm_0 = relu2d(Gemm_0)

            !========Gemm Layer============
            CALL linear_layer(Gemm_0, linLayers(3),0)

            Gemm_0 = relu2d(Gemm_0)

            !========Gemm Layer============
            CALL linear_layer(Gemm_0, linLayers(4),0)

            call CPU_TIME(T2)
            o0 = Gemm_0
        end SUBROUTINE
    !===================================================================================================================================================


END module model

!=======================================================================================
