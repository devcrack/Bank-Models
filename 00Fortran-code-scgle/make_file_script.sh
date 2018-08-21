gfortran variables_modules.f90 -c
gfortran matrix_module.f90 -c
gfortran integration_module.f90 -c
gfortran writting_module.f90 -c
gfortran structure_module.f90 -c
gfortran SCGLE_module.f90 -c
gfortran NESCGLE_module.f90 -c
gfortran program.f90 variables_modules.o structure_module.o matrix_module.o integration_module.o SCGLE_module.o writting_module.o NESCGLE_module.o -o program.out -fopenmp
