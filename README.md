# OC-GTSM_utility

This utility downloads GOFS 3.1 Ocean Global Circulation Model (OGCM) salinity and temperature (and velocity if required) data, and computes 2D variables required for a ocean-circulation coupled Global Tide and Surge Model (OC-GTSM).

See https://www.hycom.org/dataserver/gofs-3pt1/analysis for details on GOFS 3.1 (HYCOM model)

This utility requires the following libraries:
1) GSW Fortran toolbox: https://github.com/TEOS-10/GSW-Fortran
2) datetime-fortran module:  https://github.com/wavebitscientific/datetime-fortran
3) Netcdf/hdf5
4) MPI

--------------
Example compilation instructions on Linux system:

---
-> load netcdf, MPI compatible compiler, and cmake  modules
mkdir libs
cd libs
git clone https://github.com/TEOS-10/GSW-Fortran.git
cd GSW-Fortran
-> follow cmake compilation instructions on GSW-Fortran README (Make sure -DCMAKE_Fortran_COMPILER="compiler" that you will use for OC-GTSM_utility)
cd ../
git clone https://github.com/wavebitscientific/datetime-fortran.git
cd datetime-fortran
-> follow cmake compilation instructions on datetime-fortran README (Make sure -DCMAKE_Fortran_COMPILER="compiler" that you will use for OC-GTSM_utility)
cd ../../
make
---

Note: I used MVAPICH2 Intel compiler,and set the -DCMAKE_Fortran_COMPILER=ifort, i.e., use:
"cmake .. -DCMAKE_Fortran_COMPILER=ifort"
when compiling GSW and datetime
