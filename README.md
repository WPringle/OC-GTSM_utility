# OC-GTSM_utility

This utility downloads GOFS 3.1 Ocean Global Circulation Model (OGCM) salinity and temperature (and velocity if required) data, and computes 2D variables required for a ocean-circulation coupled Global Tide and Surge Model (OC-GTSM).

See https://www.hycom.org/dataserver/gofs-3pt1/analysis for details on GOFS 3.1 (HYCOM model)

This utility requires the following libraries:
1) GSW Fortran toolbox: https://github.com/TEOS-10/GSW-Fortran
2) datetime-fortran module:  https://github.com/wavebitscientific/datetime-fortran
3) Netcdf/hdf5
4) MPI
