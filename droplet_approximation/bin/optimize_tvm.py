#!/usr/bin/env python3

import os

import torch
from torch.export import export

import tvm
from tvm import relax
from tvm.contrib import cc
from tvm.relax.frontend.torch import from_exported_program

import droplet_approximation

# This number is a guess and certainly not optimal
BATCH_SIZE = 4096

# Set number of threads to exactly one for single-threaded optimizations
# This might not be necessary
os.environ["TVM_NUM_THREADS"] = "1"
tvm.runtime.num_threads = 1

# Load checkpoint
torch_model = droplet_approximation.load_model_checkpoint( "../models/devilish_taker.pt" )[0]

example_args   = (torch.randn( BATCH_SIZE, 7, dtype=torch.float32),)

with torch.no_grad():
    exported_program = export( torch_model.eval(), example_args )
    mod_from_torch = from_exported_program(
        exported_program, keep_params_as_input=False, unwrap_unit_return_tuple=True
    )


target = tvm.target.Target({
    "kind": "llvm",
    "mcpu": "znver3",          # Target Zen-3
    "mattr": ["avx2","fma"],   # Target AVX2/Fused Multiply-Add
    "num-cores": 1,            # Target single core
    "libs": ["libxsmm"]        # Should be use libxsmm. Currently a NOP
})

# opt_level=3 = '-O3' flag
with tvm.transform.PassContext(opt_level=3):
    # Lock in model weights/params
    # Note - Testing is required to see if this sequence is necessary.
    #        Some preliminary testing suggests leaving the fusion to
    #        static_shape_tuning is preferable.

    seq = tvm.transform.Sequential([
        relax.transform.FoldConstant(),
        relax.transform.FuseOps(), # Fuse MatMul + Bias + Activation while they are still Relax ops
        relax.transform.FuseTIR(),
    ])
    mod_from_torch, params_from_torch = relax.frontend.detach_params(mod_from_torch)
    mod_from_torch = relax.transform.BindParams(
        "main", params_from_torch
    )(mod_from_torch)

    mod_from_torch = seq( mod_from_torch )

    # Show unoptimized IR
    mod_from_torch.show()

    # Tune model
    optimized_mod = relax.get_pipeline(
        "static_shape_tuning",
        target=target,
        work_dir="./tuning_logs",
        total_trials=10000,       #10000
        max_trials_per_task=256,  #256
        cpu_weight_prepack=True
    )(mod_from_torch)

    # Save IR for later use
    mod_json = tvm.ir.save_json(optimized_mod)
    with open("model_ir.json", "w") as f:
        f.write(mod_json)

    ex = relax.build(optimized_mod, target, system_lib=True, exec_mode="compiled")
    ex.export_library("tvm_model.a", fcompile=cc.create_staticlib)

