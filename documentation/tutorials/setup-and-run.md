# Tutorial Overview
This tutorial demonstrates how to check out a copy of NTLP from Github, build an
executable, and then launch said executable to run the pi chamber simulation on
one of the group's clusters.

***NOTE:*** This tutorial is reliant on Notre Dame resources.  External
users need to adapt to their compilation environment (e.g. installing
dependencies and setting paths correctly), manually launching the simulation
instead of submitting it to a job scheduler, and replace the paths used when
configuring the pi chamber test case.

# Instructions

## Connect to Notre Dame
Connect to the Notre Dame network via VPN if you're off campus.

***NOTE:*** Skip this step if you're already on Notre Dame's network.

Login to the cluster's cswarmfe head node with your netID:

```shell
$ ssh -C <netID>@cswarmfe.crc.nd.edu
```

Enter your password when prompted and read the system's message of the day
(MOTD) for any important information from CRC.

## Clone the NTLP Repository
Clone the NTLP repository from Github so you have a local working copy:

```shell
$ git clone https://github.com/RichterLab/NTLP.git ~/NTLP
```

The NTLP code hosted on Github at https://github.com/RichterLab/NTLP is
downloaded a sub-directory of your home directory named NTLP.

## List the Repository Contents
Change into the working copy and list its contents.

```shell
$ cd NTLP
$ ls
defs.F  documentation/  fft.f  kdtree.f90  les.F  les.run  Makefile  netcdf_io.f90 particles.f90  postprocessing  README.md  tec_io.f90  test_cases
```

The following files and directories are of note:

| File or Directory Name | Notes |
| --- | --- |
| `defs.F` | Contains important modules |
| `documentation/` | Documentation directory tree |
| `fft.f` | Fourier transform package |
| `les.F` | This is the main code |
| `les.run` | Templated job submission script |
| `kdtree.f90` | Package which supports nearest neighbor search using a K-tree algorithm |
| `Makefile` | Rules and commands for compiling NTLP |
| `netcdf_io.f90` | Input and output routines for working with netCDF-formated files |
| `postprocessing/` | Directory containing code for post-processing NTLP outputs |
| `README.md` | Summary and initial introduction to NTLP documentation |
| `tec_io.f90` | Output routines for writing TecPlot-compatible files |
| `test_cases/` | Directory containing parameters for previously-validated test cases |

## Compiling NTLP
Now you must compile the code.  The compiler and MPI and NetCDF libraries must be
loaded into your terminal's environment:

```shell
$ module load mvapich2 intel netcdf
Loading mvapich2/2.3.1/intel/19.0
  Loading requirement: intel/19.0
```

Next, compile the code. This will take roughly a minute:

```shell
$ make -j
mpif90 -i4 -r8 -assume byterecl -xHost -fpp -traceback -O2  -c defs.F -I/opt/crc/n/netcdf/4.7.4/parallel/include -L/opt/crc/n/netcdf/4.7.4/parallel/lib
mpif90 -i4 -r8 -assume byterecl -xHost -fpp -traceback -O2  -c fft.f -I/opt/crc/n/netcdf/4.7.4/parallel/include -L/opt/crc/n/netcdf/4.7.4/parallel/lib
mpif90 -i4 -r8 -assume byterecl -xHost -fpp -traceback -O2  -c kdtree.f90 -I/opt/crc/n/netcdf/4.7.4/parallel/include -L/opt/crc/n/netcdf/4.7.4/parallel/lib
mpif90 -i4 -r8 -assume byterecl -xHost -fpp -traceback -O2  -c particles.f90 -I/opt/crc/n/netcdf/4.7.4/parallel/include -L/opt/crc/n/netcdf/4.7.4/parallel/lib
mpif90 -i4 -r8 -assume byterecl -xHost -fpp -traceback -O2  -c netcdf_io.f90 -I/opt/crc/n/netcdf/4.7.4/parallel/include -L/opt/crc/n/netcdf/4.7.4/parallel/lib
mpif90 -i4 -r8 -assume byterecl -xHost -fpp -traceback -O2  -c tec_io.f90 -I/opt/crc/n/netcdf/4.7.4/parallel/include -L/opt/crc/n/netcdf/4.7.4/parallel/lib
mpif90 -i4 -r8 -assume byterecl -xHost -fpp -traceback -O2  -c les.F -I/opt/crc/n/netcdf/4.7.4/parallel/include -L/opt/crc/n/netcdf/4.7.4/parallel/lib
mpif90 defs.o fft.o kdtree.o les.o netcdf_io.o particles.o tec_io.o -o lesmpi.a  -i4 -r8 -assume byterecl -xHost -fpp -traceback -O2  -I/opt/crc/n/netcdf/4.7.4/parallel/include -L/opt/crc/n/netcdf/4.7.4/parallel/lib -lnetcdf -lnetcdff
```

By default NTLP is compiled with optimizations and targets the Ivy Bridge
cluster.  Each of the source files is compiled individually before linking their
compiled outputs into an executable named `lesmpi.a`.

## Prepare the pi Chamber Test Case
Change into the pi chamber's test directory and update the job submission script
so it is customized to you and your environment:

```shell
$ cd test_cases/pi_chamber
$ sed -i -e "s#PATH_TO_SCRATCH#/scratch365/${USER}#" params.in
$ sed -i -e "s#PATH_TO_SCRATCH#/scratch365/${USER}#" \
         -e "s#PATH_TO_CODE#${HOME}/NTLP#" pi_chamber.run
```

Both the parameters file (params.in) and the job submission script
(`pi_chamber.run`) are generic templates. These need to be instantiated so they
refer to paths that actually exist.

## Submit the Test Case to the Cluster
Now you are ready to run the simulation. Submit the job to the scheduler:

```shell
$ qsub pi_chamber.run
Your job 211545 ("test_case_pi_chamber") has been submitted
```

Make a note of the job identifier (`211545`) as that is used to name some of the
files generated when running the simulation.  Job identifiers are also required
when interacting with the scheduling system should you need to stop a running
job or query for its current state.

## Monitor the Test Case's Output
The pi chamber test case runs 10,000 time steps and takes roughly 40 minutes on
the Ivy Bridge cluster.  You can watch the simulation's execution by examining
the contents of the log file as it is updated:

```shell
$ tail -F /scratch365/${USER}/pi_chamber/pi_chamber.out.0000000
 Starting time loop
 it,time =           20  0.623764565020764
 time,tnumpart:  0.623764565020764                0
 radmin,radmax:   0.623765E+00   0.100000E+04   0.000000E+00
tempmin,tempmax   0.623765E+00   0.100000E+04   0.000000E+00
     qmin,qmax:   0.623765E+00   0.100000E+04   0.000000E+00
   time,radavg:   0.623765E+00   0.000000E+00
time,100,1000,i   0.623765E+00           0           0           0
```

On successful completion of the simulation will indicate the job has completed
by writing lines like the following, one per processor used (64 in this
tutorial):

```shell
Job Execution Time XXX
```

Due to how MPI applications execute concurrently the lines above may not be the
last lines in the file, but rather some of the last lines.

Kill the `tail` command by pressing Ctrl-c.

## Test Case Cleanup
Cleanup the job-related files generated so they do not clutter your working
copy.

```shell
$ rm *.out test_case_pi_chamber.*
```

The per-processor files (`cou.mp.NNNNN.out`) contain information about the
sub-problem each processor handles, while the job output files
(`test_case_pi_chamber.*`) contain any errors generated when launching and running
the simulation.

Note that the simulation outputs, stored beneath
`/scratch365/<netID>/pi_chamber/`, are left on disk so they may be inspected and
visualized.

# Summary
In this tutorial you downloaded the NTLP code base to the cluster's login node,
compiled it, and ran the pi chamber test case for 10,000 steps using the Richter
Lab's Ivy Bridge cluster.  This covers the general workflow for working with
NTLP and running simulations via a cluster that is shared by multiple users.
