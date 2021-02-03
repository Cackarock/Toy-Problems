#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "mm.h"
#include "../../include/mmio.h"

int main(int argc, char *argv[]){
  
  if(argc != 4){
    print_usage();
    return -1;
  }
  
  //Initialize MPI
  int my_rank, world_size; 
  
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  
  MPI_Barrier(MPI_COMM_WORLD);


  //Call some function that reads and validates 2 matrices A and B using the matrix factory
  //The code is currently in a single file "mm.c" in the current folder
  //I need another function, lets call it "mat-mult" that takes 4 arguments: int dimension, double* a, double* b, double* c 
  //and multiples a and b and stores it in c. 
  int mm_return = mpi_matrix_mult(argv, my_rank, world_size);
  if(mm_return != 0){
    MPI_Finalize();
    return mm_return;
  }


  MPI_Finalize();
  return 0;
}
