# LES-GPU

## Compilation
To build on the CRC machines you must run the following commands:
```
module load cmake/3.2.2
module load mvapich2/2.1-intel-15.0-mlx
module load intel/15.0

mkdir build
cd build 
cmake ..
make
```
