#include <tvm/runtime/module.h>
#include <tvm/ffi/function.h> // Modern replacement for PackedFunc
#include <tvm/ffi/extra/c_env_api.h>
#include <tvm/ffi/container/tensor.h>
#include <tvm/runtime/memory/memory_manager.h>
#include <tvm/runtime/object.h>
#include <cstring>

// Cache the function so we don't reload it

static tvm::ffi::Function g_inference_func;

void nop(DLManagedTensor* self) { }

extern "C" {
    int INPUT_COUNT;
    int OUTPUT_COUNT;
    int BATCH_SIZE;

    int OUTPUT_BUFFER_SIZE;

    static DLManagedTensor input_tensor;
    static tvm::ffi::Tensor input_view;
    static int64_t input_shape[2];

    int init_inference( int input_count, int output_count, int batch_size ) {
        INPUT_COUNT  = input_count;
        OUTPUT_COUNT = output_count;
        BATCH_SIZE   = batch_size;

        INPUT_BUFFER_SIZE  = 4 * BATCH_SIZE * INPUT_COUNT;  // Allocate for 32 bit floats
        OUTPUT_BUFFER_SIZE = 4 * BATCH_SIZE * OUTPUT_COUNT; // Allocate for 32 bit floats

        // Access the global registry using tvm::ffi::GetGlobal. This gets the static droplet model library.
        auto system_lib_getter = tvm::ffi::Function::GetGlobal("ffi.SystemLib");
        
        if (!system_lib_getter) {
            std::cerr << "C++ Error: SystemLib not found. "
                      << "Ensure you linked with -Wl,--whole-archive your_model.a" << std::endl;
            return -1;
        }

        // Invoke the getter to obtain the System Module
        tvm::ffi::Module lib = system_lib_getter.value()().cast<tvm::ffi::Module>();

        // Get the VM loader from the module
        auto vm_load_executable = lib->GetFunction("vm_load_executable");
        if (!vm_load_executable) {
            std::cerr << "C++ Error: vm_load_executable not found in static library." << std::endl;
            return -1;
        }

        // 4. Instantiate the Virtual Machine
        tvm::ffi::Module relax_vm = vm_load_executable.value()().cast<tvm::ffi::Module>();

        // 5. Initialize the VM
        auto vm_initialization = relax_vm->GetFunction("vm_initialization");
        if (!vm_initialization) {
            std::cerr << "C++ Error: VM initialization function not found." << std::endl;
            return -1;
        }

        // Arguments: (device_type, device_id, allocator_type)
        // kDLCPU is typically 1. kPooled is typically 1.
        vm_initialization.value()(1, 0, 1);
        g_inference_func = relax_vm->GetFunction("main").value();

        input_shape[0]         		    = BATCH_SIZE;
        input_shape[1]         		    = INPUT_COUNT;
        input_tensor.dl_tensor.device       = {kDLCPU, 0};
        input_tensor.dl_tensor.ndim         = 2;
        input_tensor.dl_tensor.dtype        = {kDLFloat, 32, 1};
        input_tensor.dl_tensor.shape        = input_shape;
        input_tensor.dl_tensor.strides      = nullptr;
        input_tensor.dl_tensor.byte_offset  = 0;
        input_tensor.deleter                = nop;
        input_tensor.manager_ctx            = nullptr;

        return 0;
    }

    void run_inference(float* in_data, float* out_data) {
        // TVM Runtime expects a tensor, so wrap our data accordingly
        input_tensor.dl_tensor.data  = in_data;

        auto input_view  = tvm::ffi::Tensor::FromDLPack( &input_tensor );

        float* output = (float *)g_inference_func(input_view).cast<tvm::ffi::Tensor>().data_ptr();

        std::memcpy( out_data, output, OUTPUT_BUFFER_SIZE );
    }
}
