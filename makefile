FC = mpif90
O_DIR = #odir1/ 
GSWHOME       := ../lib/GSW-Fortran-master/
DATETIMEHOME  := ../lib/datetime-fortran-master/build/
HDF5HOME=/opt/crc/h/hdf5/intel/17.1/build/lib/ 
NETCDFHOME=/afs/crc.nd.edu/x86_64_linux/n/netcdf/4.4.1/intel-17.1/build/
incdir = -I$(DATETIMEHOME)include -I$(NETCDFHOME)include -I$(GSWHOME)modules
libdir = -L$(NETCDFHOME)/lib -L$(HDF5HOME) -L$(GSWHOME) -L$(DATETIMEHOME)lib/
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
