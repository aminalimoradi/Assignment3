mpiifort main.f90 worker.f90 manager.f90 -o amin.x -qmkl
mpirun -np 4 ./amin.x

