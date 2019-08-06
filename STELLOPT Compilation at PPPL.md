STELLOPT Compilation at the PPPL cluster
========================================

![](images/PPPL-LOGO-FNLWH-GRADIENT_300px_WEB.jpg)

This page details how to compile the STELLOPT family of codes at the
[PPPL](@http://www.pppl.gov/) cluster. The proper modules need to be
loaded then compilation can begin. Please note that if you require
additional module to be loaded you should do this before loading the
compiler. This will prevent the \$PATH variable from searching the wrong
directories.

------------------------------------------------------------------------

### GENERAL Instructions

On the PPPL cluster, STELLOPT is now maintained as an installed module.
Use the following command to load stellopt

    module load stellopt

To load the most current version of the code. This will load all the
necessary modules. The /bin/ directory will be added to your path
variable as well so codes can be called without specifying the full
path. Such as:

    mpirun -np 64 xstelloptv2 input.test

------------------------------------------------------------------------

### GNU

Load the appropriate module files for the [GNU](@https://gcc.gnu.org/)
compiler.

    module purge
    module load git
    module load gcc/6.1.0
    module load gsl/1.16 silo nag/fll6a24dfl pgplot szip/2.1 ntcc/gcc6.1.0
    module load acml/6.1.0/gfortran64 hdf/4.2.12
    module load openmpi/3.0.0
    module load blacs fftw hdf5-parallel/1.10.1 scalapack
    module load petsc_complex/3.8.3
    module load slepc_complex/3.8.3
    module load curl/7.59.0 netcdf-c/4.4.1
    module load netcdf-fortran/4.4.4 netcdf-cxx4/4.3.0
    module load qt/4.8.3
    module load python/3.6.4

You will now be able to follow the instructions on the main compilation
page to compile the code. Please be sure to point your the make.inc at
SHARE/make\_pppl.inc.