FORTRAN=mpiifort
#FORTRAN=mpif90

# Specify the architecture to optimize for.  Should be one of the following:
#
#   host    Optimize for the current system's architecture.  This cannot
#           be used when running on AMD processors as the resulting executable
#           will crash with an illegal instruction error.
#   avx     Optimize for Intel's AVX instruction set.  This is the most advanced
#           optimization for Intel Ivy Bridge processors and AMD's Bulldozer
#           processors.
#   avx2    Optimize for Intel's AVX2 instruction set.  This is the most advanced
#           optimization for Intel Haswell processors, or newer, and AMD's Epyc
#           7001, 7002, and 7003 (Zen 1-3) processors.
#
# NOTE: Specifying this incorrectly will, at best, result in degraded execution
#       times and, at worst, result in crashes due to illegal instructions.
#
ARCH ?= host

ifeq ($(ARCH), host)
ARCH_FLAGS = -xHost
endif
ifeq ($(ARCH), avx)
ARCH_FLAGS = -march=corei7-avx
endif
ifeq ($(ARCH), avx2)
ARCH_FLAGS = -march=core-avx2
endif

# Always used compilation flags:
#
#   -i4               Use 32-bit integers
#   -r2               Use 64-bit floating point values
#   -assume byterecl  Interpret open()'s recl option in bytes to eliminate
#                     Fortran record headers being used in binary files.
#   $(ARCH_FLAGS)     Architecture-specific flags, from above
#   -fpp              Enable pre-processing of source files regardless of extension
#   -traceback        Generate detailed stack traces when run-time errors occur.
#                     This doesn't have a performance cost as it only modifies
#                     the error code path.  Note that it does not depend on the
#                     level of debugging support.
#
FLAGS=-i4 -r8 -assume byterecl $(ARCH_FLAGS) -fpp -traceback -fc=ifx

# Are we building a debug build?  This enables options useful for debugging
# the solver's behavior but are not desirable to unconditionally enable.
# See below for details.
ifeq ($(DEBUG), yes)

# Generate debugging information.  This is necessary for running the solver
# under a debugger and can be useful for other tools (e.g. profilers).
#
# NOTE: This disables optimizations and enables run-time checks!
#
DEBUG_FLAGS = -g

# Enable all compilation warnings except for when temporary arrays are created
# when passing to a subroutine or Fortran.  Temporary arrays occur throughout
# the code base and will be addressed at a later date.
DEBUG_FLAGS += -check all,noarg_temp_created

# Exit with a SIGFPE whenever a floating point exception (FPE) is detected.
# This is useful for identifying precisely where an invalid value (infinities
# and NaNs) are introduced.
DEBUG_FLAGS += -fpe0

# Initialize floating point values, both scalars and arrays, with signalling
# NaNs.  Combined with exiting on FPEs this makes it trivial to identify the use
# of uninitialized floating point values.
DEBUG_FLAGS += -init=arrays -init=snan

else # DEBUG == no

# Enable optimizations that are almost always beneficial.
FLAGS += -O2

endif

# NetCDF output is always enabled.  Make sure we can find its headers and
# libraries.
ifeq ($(strip $(NETCDFBASE)),)

# Don't enforce the variable being set if we just want to clean the directory.
ifneq ($(strip $(MAKECMDGOALS)), clean)
$(info NETCDFBASE is not set!  Did you forget to load a module (e.g. mvapich2)? )
$(info )
$(error "Please set NETCDFBASE to the location of the netCDF install." )
endif


endif

OUTPUTINC = -I$(NETCDFBASE)/include
OUTPUTLIB = -L$(NETCDFBASE)/lib
LINKOPTS  = -lnetcdf -lnetcdff

# TecPlot output is only used when requested.
TECPLOT ?= no
TECFLAGS = -DTECIO
TECLIB   = ~/Research/tecio/libteciompi.a
TECLINK  = -lm -lstdc++ -lgcc_eh

SRC = data_structures.f90 \
      defs.F \
      fft.f \
      kdtree.f90 \
      les.F \
      measurement.f90 \
      netcdf_io.f90 \
      particles.f90 \
      tec_io.f90

OBJS = $(addsuffix .o, $(basename $(SRC)))

ifeq ($(TECPLOT), yes)
FLAGS    += $(TECFLAGS)
LINKOPTS += $(TECLINK) $(TECLIB)
endif

lesmpi.a: $(OBJS)
	$(FORTRAN) $^ -o $@  $(FLAGS) $(DEBUG_FLAGS) $(OUTPUTINC) $(OUTPUTLIB) $(LINKOPTS)

%.o: %.f
	$(FORTRAN) $(FLAGS) $(DEBUG_FLAGS) -c $< $(OUTPUTINC) $(OUTPUTLIB)

%.o: %.f90
	$(FORTRAN) $(FLAGS) $(DEBUG_FLAGS) -c $< $(OUTPUTINC) $(OUTPUTLIB)

%.o: %.F
	$(FORTRAN) $(FLAGS) $(DEBUG_FLAGS) -c $< $(OUTPUTINC) $(OUTPUTLIB)


clean:
	rm -f *.o *.mod lesmpi.a mach.file

# Dependencies between the individual objects.
les.o: defs.o measurement.o netcdf_io.o particles.o tec_io.o
measurement.o: data_structures.o
particles.o: defs.o measurement.o
netcdf_io.o: particles.o
tec_io.o: particles.o
