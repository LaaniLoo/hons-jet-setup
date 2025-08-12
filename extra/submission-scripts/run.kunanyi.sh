#!/bin/bash -l
#PBS -lselect=100:ncpus=28:mpiprocs=28
#PBS -lwalltime=48:00:00
#PBS -N 3D-off1r-30
#PBS -m abe
#PBS -M patrick.yates@utas.edu.au

cd $PBS_O_WORKDIR
module load intel impi hdf5 libpng

mkdir -p /u/pmyates/pluto-simulations/asymmetric-jets/3D/off1r/data-off1r-30
mkdir -p /scratch/pmyates/simulations/log/3D/off1r-30

# export HDF5_USE_FILE_LOCKING="FALSE"

# HDF5_USE_FILE_LOCKING="FALSE" mpiexec ./pluto -i pluto.off1r.ini -h5restart 19

mpiexec ./pluto -i
