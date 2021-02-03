# PA2: MPI Monte Carlo and Cannon's algorithm

* Team 4: Geoffrey Meier, Floriana Ciaglia, Ravi Shankar, Dylan Leman
* Class: CS 430/530 - Parallel Computing
* Date: Fall '20

## Overview

Parallelization:
- PI: Monte Carlo - splits between processors to calculate pi
- Matrix-Matrix Multiply - Multiplies square matrices using Cannon's Algorithm

## Manifest

```
├── README.md
├── Makefile
├── slurm_MPI.bash
├── bin
├── lib
├── obj
├── include
│   ├── mmio.c
│   └── mmio.h
├── src
│   ├── mm
│   │   ├── mm.c
│   │   ├── mm.h
│   │   └── parallel-mm.c
│   └── pi-mc
│       ├── pi-mc.c
│       ├── pi-mc.h
│       └── pi-mc-parallel-main.c
├── input
│   ├── expected2x2.m
│   ├── expected32x32.m
│   ├── expected9x9.m
│   ├── generateMatrix
│   ├── matrix1000x1000.m
│   ├── matrix100x100.m
│   ├── matrix2x2.m
│   ├── matrix32x32.m
│   ├── matrix4x4.m
│   ├── matrix500x500.m
│   ├── matrix50x50.m
│   ├── matrix9x9.m
└── test
    ├── canons-test.cpp
    ├── matrix-test.sh
    ├── pi-mc-test.cpp
    └── test-output
```

## Building and Running

* Modules required:
  * `gcc/7.2.0`
  * `slurm/17.11.12`
* After cloning you will need to initialize submodules for googletest:
  * `git submodule init`
  * `git submodle update`
* To build all targets, use: `make all`
* To build in DEBUG mode, use: `make debug`
* To build individual targets:
  * `make pi-mc`
  * `make mm`
* Run as follows: `bin/<target-name>`


## Testing

* A combination of GoogleTest and a bash script is used for testing.
  * Monte Carlo Pi estimation program has tests written in googletest
  * Matrix-Matrix multiplication program is tested with a bash script
* All tests are located in the `test/` folder in the root project directory.
* To make and run all tests, use: `make test`

### Test 4 Cores

Run the following:

`$ mpirun -np 4 bin/mm input/matrix2x2.m input/matrix2x2.m output2x2.m`

Perform `diff` between this output file and `input/expected2x2.m`

### Test 9 Cores

Run the following:

`$ mpirun -np 9 bin/mm input/matrix9x9.m input/matrix9x9.m output9x9.m`

Perform `diff` between this output file and `input/expected9x9.m`

### Test 16 Cores

Run the following:

`$ mpirun -np 16 bin/mm input/matrix32x32.m input/matrix32x32.m output32x32.m`

Perform `diff` between this output file and `input/expected32x32.m`


## Known Bugs

Having other modules loaded may cause build issues. It is recommended that only `gcc/7.2.0` and `mpich/ge/gcc/64/3.2.1` be loaded. 

## Sources

* [Cannon's Algorithm Distributed Matrix Multiplication](https://iq.opengenus.org/cannon-algorithm-distributed-matrix-multiplication/)

* [Cannon's Algorithm Matrix Multiplication with MPI Topologies](http://boron.physics.metu.edu.tr/ozdogan/GraduateParallelComputing.old/week11/node3.html)

* [Cannon's Algorithm](http://cseweb.ucsd.edu/classes/fa12/cse260-b/Lectures/Lec13.pdf)

* Dr. Olschanowsky's Cannon's code

