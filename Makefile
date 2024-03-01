FORTRAN=mpif90
F90=ifort

#FLAGS=-i4 -r8 -O2 -assume byterecl -xHost -fpp
FLAGS=-i4 -r8 -O2 -assume byterecl -march=core-avx2 -fpp

## UNCOMMENT TO RUN IN DEBUG MODE
DEBUG_FLAGS=-g -traceback
#DEBUG_FLAGS+=-check bounds

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
