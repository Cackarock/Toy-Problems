#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>
#include "pi-mc.h"


int main(int argc, char** argv){
  int numprocs, id, iterations;

  MPI_Init(&argc, &argv);

  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &id);
		
  if(argc != 2){
    printf("Usage: ./pi <iterations>\n");
    return(-1);
  }

  iterations = atoi(argv[1]);

  if(iterations == 0){
    printf("Provide a non-0 integer interval\n");
    exit(1);
  }

  pi_proc(iterations, id, numprocs);
	
	MPI_Finalize();

  return 0;
}

