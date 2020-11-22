FORTRAN=mpif90
F90=ifort

FLAGS=-i4 -r8 -O2 -assume byterecl
DEBUG_FLAGS=-g -traceback -check all

OUTPUTINC = -I$(NETCDFBASE)/include
OUTPUTLIB = -L$(NETCDFBASE)/lib
OUTPUTOPT = -DNETCDF -DNCFPLUS
LINKOPTS  = -lnetcdf -lnetcdff

SRC = fft.f \
      kdtree.f90

OBJS = $(addsuffix .o, $(basename $(SRC)))


all:	$(OBJS)
	$(FORTRAN) les.F -o lesmpi.a  $(OBJS) $(FLAGS) $(OUTPUTINC) $(OUTPUTLIB) $(LINKOPTS)

%.o:	
	$(F90) $(FLAGS) $(SRC) -c

debug:  $(OBJS)
	$(FORTRAN) $(FLAGS) $(DEBUG_FLAGS) $(INCLUDE) les.F $(OFILES)

clean:
	rm -f *.o *.mod lesmpi.a mach.file *.*~


