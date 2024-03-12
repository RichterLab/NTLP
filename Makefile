FORTRAN=mpif90
F90=ifort

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

FLAGS=-i4 -r8 -O2 -assume byterecl $(ARCH_FLAGS) -fpp

## UNCOMMENT TO RUN IN DEBUG MODE
#DEBUG_FLAGS=-g -traceback
#DEBUG_FLAGS+=-check all,noarg_temp_created
#DEBUG_FLAGS+=-fpe0
#DEBUG_FLAGS+=-init=arrays -init=snan

# NetCDF output is always enabled.
OUTPUTINC = -I$(NETCDFBASE)/include
OUTPUTLIB = -L$(NETCDFBASE)/lib
LINKOPTS  = -lnetcdf -lnetcdff

# TecPlot output is only used when requested.
TECPLOT ?= no
TECFLAGS = -DTECIO
TECLIB   = ~/Research/tecio/libteciompi.a
TECLINK  = -lm -lstdc++ -lgcc_eh

SRC = defs.F \
      fft.f \
      kdtree.f90 \
      les.F \
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
les.o: defs.o netcdf_io.o particles.o tec_io.o
particles.o: defs.o
netcdf_io.o: particles.o
tec_io.o: particles.o
