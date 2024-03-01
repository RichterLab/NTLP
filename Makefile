FORTRAN=mpif90
F90=ifort

#FLAGS=-i4 -r8 -O2 -assume byterecl -xHost -fpp
FLAGS=-i4 -r8 -O2 -assume byterecl -march=core-avx2 -fpp

## UNCOMMENT TO RUN IN DEBUG MODE
DEBUG_FLAGS=-g -traceback
#DEBUG_FLAGS+=-check bounds

##UNCOMMENT TO RUN WITH TECPLOT I/O
#Provide location of the mpi-enabled tecio library
#TECINCLUDE=~/Research/tecio/libteciompi.a
#TECLINK=-lm -lstdc++ -lgcc_eh -DTECIO

OUTPUTINC = -I$(NETCDFBASE)/include
OUTPUTLIB = -L$(NETCDFBASE)/lib
LINKOPTS  = -lnetcdf -lnetcdff


SRC = defs.F \
      fft.f \
      kdtree.f90 \
      les.F \
      netcdf_io.f90 \
      particles.f90 \
      tec_io.f90

OBJS = $(addsuffix .o, $(basename $(SRC)))


lesmpi.a: $(OBJS)
	$(FORTRAN) $^ -o $@  $(FLAGS) $(DEBUG_FLAGS) $(OUTPUTINC) $(OUTPUTLIB) $(LINKOPTS) $(TECLINK) $(TECINCLUDE)

%.o: %.f
	$(FORTRAN) $(FLAGS) $(DEBUG_FLAGS) -c $< $(OUTPUTINC) $(OUTPUTLIB) $(TECLINK)

%.o: %.f90
	$(FORTRAN) $(FLAGS) $(DEBUG_FLAGS) -c $< $(OUTPUTINC) $(OUTPUTLIB) $(TECLINK)

%.o: %.F
	$(FORTRAN) $(FLAGS) $(DEBUG_FLAGS) -c $< $(OUTPUTINC) $(OUTPUTLIB) $(TECLINK)


clean:
	rm -f *.o *.mod lesmpi.a mach.file

# Dependencies between the individual objects.
les.o: defs.o netcdf_io.o particles.o
particles.o: defs.o
netcdf_io.o: particles.o
tec_io.o: particles.o
