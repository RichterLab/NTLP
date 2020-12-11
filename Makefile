FORTRAN=mpif90
F90=ifort

FLAGS=-i4 -r8 -O2 -assume byterecl

## UNCOMMENT TO RUN IN DEBUG MODE
# DEBUG_FLAGS=-g -traceback -check bounds

OUTPUTINC = -I$(NETCDFBASE)/include
OUTPUTLIB = -L$(NETCDFBASE)/lib
OUTPUTOPT = -DNETCDF -DNCFPLUS
LINKOPTS  = -lnetcdf -lnetcdff

SRC = 	fft.f \
	kdtree.f90 \
        defs.F \
        netcdf_io.f90

OBJS = $(addsuffix .o, $(basename $(SRC)))


lesmpi.a: $(OBJS) les.F
	$(FORTRAN) les.F -o $@ $^ $(FLAGS) $(DEBUG_FLAGS) $(OUTPUTINC) $(OUTPUTLIB) $(LINKOPTS)

%.o: %.f
	$(FORTRAN) $(FLAGS) $(DEBUG_FLAGS) -c $< $(OUTPUTINC) $(OUTPUTLIB)

%.o: %.f90
	$(FORTRAN) $(FLAGS) $(DEBUG_FLAGS) -c $< $(OUTPUTINC) $(OUTPUTLIB)

%.o: %.F
	$(FORTRAN) $(FLAGS) $(DEBUG_FLAGS) -c $< $(OUTPUTINC) $(OUTPUTLIB)


clean:
	rm -f *.o *.mod lesmpi.a mach.file

