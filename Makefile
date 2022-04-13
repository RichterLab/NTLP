FORTRAN=mpif90
F90=ifort

FLAGS=-i4 -r8 -O2 -assume byterecl -xHost

## UNCOMMENT TO RUN IN DEBUG MODE
DEBUG_FLAGS=-g -traceback 

OUTPUTINC = -I$(NETCDFBASE)/include
OUTPUTLIB = -L$(NETCDFBASE)/lib
OUTPUTOPT = -DNETCDF -DNCFPLUS
LINKOPTS  = -lnetcdf -lnetcdff

TECINCLUDE=~/Research/tecio/libteciompi.a
TECLINK=-lm -lstdc++ -lgcc_eh

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
	$(FORTRAN) $(FLAGS) $(DEBUG_FLAGS) -c $< $(OUTPUTINC) $(OUTPUTLIB) $(TECLINK) $(TECINCLUDE)

%.o: %.f90
	$(FORTRAN) $(FLAGS) $(DEBUG_FLAGS) -c $< $(OUTPUTINC) $(OUTPUTLIB) $(TECLINK) $(TECINCLUDE)

%.o: %.F
	$(FORTRAN) $(FLAGS) $(DEBUG_FLAGS) -c $< $(OUTPUTINC) $(OUTPUTLIB) $(TECLINK) $(TECINCLUDE)


clean:
	rm -f *.o *.mod lesmpi.a mach.file

