FORTRAN=mpif90
F90=ifort

FLAGS=-i4 -r8 -O2 -assume byterecl
DEBUG_FLAGS=-g -traceback -check all

OUTPUTINC = -I$(NETCDFBASE)/include
OUTPUTLIB = -L$(NETCDFBASE)/lib
OUTPUTOPT = -DNETCDF -DNCFPLUS
LINKOPTS  = -lnetcdf -lnetcdff

SRC = 	fft.f \
	kdtree.f90 \
        defs.F \
        netcdf_io.f90

OBJS = $(addsuffix .o, $(basename $(SRC)))


les.F:	$(OBJS)
	$(FORTRAN) les.F -o lesmpi.a  $(OBJS) $(FLAGS) $(OUTPUTINC) $(OUTPUTLIB) $(LINKOPTS)

%.o:	
	$(FORTRAN) $(FLAGS) $(SRC) -c $(OUTPUTINC) $(OUTPUTLIB)

clean:
	rm -f *.o *.mod lesmpi.a mach.file *.*~

