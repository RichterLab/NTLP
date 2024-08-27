#:set architecture = [('Gemm', ['Gemm_0', 1], None), ('Relu', ['Gemm_0'], [2]), ('Gemm', ['Gemm_0', 1], None), ('Relu', ['Gemm_0'], [2]), ('Gemm', ['Gemm_0', 1], None), ('Relu', ['Gemm_0'], [2]), ('Gemm', ['Gemm_0', 1], None)]
#:set inputs = []
#:set trueInputs = [['Gemm_0', [1, 6]]]
#:set outShape = [['15', [1, 2]]]
#:set outputs = {'15': 'Gemm_0'}
