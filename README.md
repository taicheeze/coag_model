#Introduction
This is the repository for the 2D model for coagulation used by Petzold Group researchers at UCSB. The code was developed by Dr. Sheng Wu and Dr. Matthew Buoni during their time at UCSB.

##download and install petsc
This model uses PETSc as its solver. This code was originally run with petsc version 3.4.0, and we cannot guarantee that it will work with the same efficiency with current or future versions of petsc. PETSc can be downloaded from the [official website](https://www.mcs.anl.gov/petsc/download/index.html)
./configure PETSC_ARCH=arch-linux2-c-opt --download-superlu_dist --download-parmetis --download-metis --with-debugging=no

# DO NOT USE THIS ONE./configure PETSC_ARCH=arch-linux2-c-opt --with-cc=mpicc --with-fc=mpif90 --download-f-blas-lapack --download-mpich --download-superlu_dist --download-parmetis --download-metis --with-debugging=no

make all test

#add $PETSC_DIR and $PETSC_ARCH into ~/.bash_profile or ~/.bashrc

#make coagulation solver
make coagulationo3
