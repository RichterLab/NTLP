# NTLP (NCAR Turbulence with Lagrangian Particles)

## Compilation
To build on the CRC machines you must run the following commands:
```
module load mvapich2
module load intel
module load netcdf

make clean
make
```

## SETUP AND RUNNING
To run, make a directory ("case1" or something) where les.run and params.in will go (i.e., not out of the same directory as the code files)

Set up the I/O directory on scratch: /scratch365/netID/tutorial/case1

Make sure all paths in params.in and les.run point to the proper locations


