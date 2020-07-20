SHELL:=/bin/bash

FC = mpif90
O_DIR =  
LIBHOME    = $(PWD)/libs/
HDF5HOME   = /afs/crc.nd.edu/x86_64_linux/hdf/hdf5-1.8.6-linux-x86_64-static/lib
NETCDFHOME = /afs/crc.nd.edu/x86_64_linux/n/netcdf/4.7.0/intel/18.0/
GSWHOME       := $(LIBHOME)GSW-Fortran/build/
DATETIMEHOME  := $(LIBHOME)datetime-fortran/build/
incdir = -I$(DATETIMEHOME)include -I$(NETCDFHOME)include -I$(GSWHOME)gsw/
libdir = -L$(NETCDFHOME)lib -L$(HDF5HOME) -L$(GSWHOME) -L$(DATETIMEHOME)lib/
FFLAGS = -O2 -c $(incdir) -traceback #-g -check bounds 
LFLAGS = -O2 $(libdir) -lnetcdf -lnetcdff -lgsw -ldatetime
LINK = $(FC)
TARGET = OGCM_DL.a

$(O_DIR):
	mkdir -p $@

SRC  =  OGCM_DL.f90
	
OBJ:=   $(patsubst %.f90, $(O_DIR)%.o, $(SRC) )

all: $(TARGET)

clean: 
	-rm -rf $(O_DIR)*.o

$(TARGET): $(OBJ)
	$(LINK) -o $@ $(OBJ) $(LFLAGS) 

$(O_DIR)%.o : %.f90
	$(FC) $(FFLAGS) $< -o $@
