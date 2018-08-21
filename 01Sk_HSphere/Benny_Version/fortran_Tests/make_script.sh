gfortran variables_modules.f90 -c
gfortran structure_module.f90 -c
gfortran writting_module.f90 -c
gfortran -fopenmp main_program.f90 variables_modules.o structure_module.o writting_module.o -o program.out 
