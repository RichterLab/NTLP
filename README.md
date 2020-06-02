# NTLP (NCAR Turbulence with Lagrangian Particles)

## Compilation
To build on the CRC machines you must run the following commands:
```
module load cmake
module load mvapich2
module load intel

mkdir build
cd build 
cmake ..
make
```

## SETUP AND RUNNING
To run, make a directory ("case1" or something) where les.run and params.in will go
(i.e., not out of the same directory as les.F)
Make sure all paths in these directories point to the proper locations

