FORTRAN=mpif90
F90=ifort

## UNCOMMENT ONLY IF RUNNING IN HPC OR LONG QUEUE
FLAGS=-i4 -r8 -O2 -assume byterecl -fpp

## UNCOMMENT IF RUNNING IN RICHTER QUEUE
#FLAGS=-i4 -r8 -O2 -assume byterecl -xHost -fpp

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


SRC = 	fft.f \
	kdtree.f90 \
        defs.F \
        particles.f90\
        netcdf_io.f90\
        tec_io.f90

OBJS = $(addsuffix .o, $(basename $(SRC)))


lesmpi.a: $(OBJS) les.F
	$(FORTRAN) $^ -o $@  $(FLAGS) $(DEBUG_FLAGS) $(OUTPUTINC) $(OUTPUTLIB) $(LINKOPTS) $(TECLINK) $(TECINCLUDE)

%.o: %.f
	$(FORTRAN) $(FLAGS) $(DEBUG_FLAGS) -c $< $(OUTPUTINC) $(OUTPUTLIB) $(TECLINK)

%.o: %.f90
	$(FORTRAN) $(FLAGS) $(DEBUG_FLAGS) -c $< $(OUTPUTINC) $(OUTPUTLIB) $(TECLINK)

%.o: %.F
	$(FORTRAN) $(FLAGS) $(DEBUG_FLAGS) -c $< $(OUTPUTINC) $(OUTPUTLIB) $(TECLINK)


clean:
	rm -f *.o *.mod lesmpi.a mach.file

