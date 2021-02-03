//////////////////////////////////////////////////////////////////////////////////////
//Authors: Kenslee Moy, Ravi Shankar, Dylan Leman, Floriana Ciaglia, Geoffrey Meier //
//////////////////////////////////////////////////////////////////////////////////////
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <string.h>
#include <sys/time.h>
#include "mm.h"
#include "../../include/mmio.h"

/** 
 * Prints program usage
 */
void print_usage(void)
{
  printf("\nUsage:");
  printf("\n      ./bin/mm <matrix file 1> <matrix file 2> <output file>\n\n");
  return;
}

// Simple Serial
int mm(double *A, int c, int rc, int r, double *B, double *C)
{

  for (int i = 0; i < c; i++)
  {
    for (int j = 0; j < r; j++)
    {
      for (int k = 0; k < rc; k++)
      {
        ARRAY(C, c, r, i, j) += ARRAY(A, c, rc, i, k) * ARRAY(B, rc, r, k, j);
      }
    }
  }
  return 0;
}

/**
 * Calls all methods reqd for cannon's multiplication of
 * two matrices
 *
 */
int mpi_matrix_mult(char *argv[], int my_rank, int world_size)
{
  return read_matrix(argv, my_rank, world_size);
}

int read_matrix(char *argv[], int my_rank, int world_size)
{

  FILE *m1, *m2, *output;
  MM_typecode m1matcode;
  MM_typecode m2matcode;
  int m1M, m2M, m2N, m1Size, m2Size, outputMSize;

  m1 = fopen(argv[1], "r");
  m2 = fopen(argv[2], "r");
  output = fopen(argv[3], "w");

  if (m1 == NULL || m2 == NULL)
  {
    fprintf(stderr, "Usage: ./bin/mm <matrix file 1> <matrix file 2> <output file>\n");
    return (-1);
  }

  if ((mm_read_banner(m1, &m1matcode) != 0) || (mm_read_banner(m2, &m2matcode) != 0))
  {
    fprintf(stderr, "Usage: ./bin/mm <matrix file 1> <matrix file 2> <output file>\n");
    printf("Could not process Matrix Market banner from files.\n");
    return (-2);
  }

  if (mm_is_complex(m1matcode) || mm_is_coordinate(m1matcode))
  {
    printf("Sorry, this application does not support ");
    printf("Market Market type: [%s]\n", mm_typecode_to_str(m1matcode));
    return (-3);
  }

  if (mm_is_complex(m2matcode) || mm_is_coordinate(m2matcode))
  {
    printf("Sorry, this application does not support ");
    printf("Market Market type: [%s]\n", mm_typecode_to_str(m2matcode));
    return (-4);
  }

  if ((mm_read_mtx_array_size(m1, &m1M, &m1M) != 0) || (mm_read_mtx_array_size(m2, &m2M, &m2N) != 0))
    return (-5);

  m1Size = m1M * m1M;
  m2Size = m2M * m2N;
  outputMSize = m1M * m2N;

  if (m1M != m2M)
  {
    fclose(m1);
    fclose(m2);
    fprintf(stderr, "Cannot perform matrix multiply on given matrices\n");
    return (-6);
  }

  double matrix1[m1Size], matrix2[m2Size], C[outputMSize];

  //Read A - matrix1
  for (int i = 0; i < m1M; i++)
  {
    for (int j = 0; j < m1M; j++)
    {
      fscanf(m1, "%lf", &matrix1[(m1M * i) + j]);
    }
  }

  //Read B - matrix2
  for (int i = 0; i < m2N; i++)
  {
    for (int j = 0; j < m2M; j++)
    {
      fscanf(m2, "%lf", &matrix2[(m2N * i) + j]);
    }
  }

  fclose(m1);
  fclose(m2);

  cannon(matrix1, matrix2, C, my_rank, world_size, m1M);

  MPI_Barrier(MPI_COMM_WORLD);
  if (my_rank==0) {
    write_matrix(C, output, m1M, outputMSize, m1matcode);
  }

  return 0;
}

//Matrix A and B are now in row major form and ready to be split into blocks
//This portion will create blocks for A and B so that each process will have its own block
int cannon(double *matrix1, double *matrix2, double *C, int my_rank, int world_size, int m1M)
{

  int block_size = m1M / sqrt(world_size);
  int num_blocks_row_col = m1M / block_size; //number of blocks to skip over
  int dim = sqrt(world_size);
  int myi = my_rank / dim;
  int myj = my_rank % dim;
  MPI_Status status;

  // How much space do I need?

  int data_dim = m1M / dim;
  if (dim * dim != world_size)
  {
    fprintf(stderr, "ERROR: Number of processes must be a perfect square!");
    return 1;
  }
  if (m1M % dim)
  {
    fprintf(stderr, "ERROR: Matrix size not compatible with number of processes specified!\n");
    return 2;
  }

  int mySize = data_dim * data_dim;

  //TODO: Malloc space for final result C

  double *block_a, *block_b, *block_c, *tmp;
  block_a = (double *)malloc(block_size * block_size * sizeof(double));
  block_b = (double *)malloc(block_size * block_size * sizeof(double));
  block_c = (double *)malloc(block_size * block_size * sizeof(double));
  tmp = (double *)malloc(sizeof(double) * block_size * block_size);

  //Check that malloc didn't fail
  if (block_a == NULL || block_b == NULL || block_c == NULL || tmp == NULL)
  {
    fprintf(stderr, "malloc failed\n");
    return 3;
  }

  //Starting indices for each block depends on world_rank
  int block_col_index = my_rank % num_blocks_row_col;
  int block_row_index = (my_rank - block_col_index) / num_blocks_row_col;

  int index = 0;
  int row, col;

  //Populate blocks a, b and c
  for (row = block_row_index * block_size; row < block_row_index * block_size + block_size; row++)
  {
    for (col = block_col_index * block_size; col < block_col_index * block_size + block_size; col++)
    {
      block_a[index] = matrix1[row * m1M + col];
      block_b[index] = matrix2[row * m1M + col];
      block_c[index] = 0;

#ifdef DEBUG
      fprintf(stderr, "[Rank %d] block_a[%d]: matrix1[%d * %d + %d = %d] = %f\n", my_rank, index, row, m1M, col, row * m1M + col, matrix1[row * m1M + col]);
      fprintf(stderr, "[Rank %d] block_a[%d]: matrix2[%d * %d + %d = %d] = %f\n", my_rank, index, row, m1M, col, row * m1M + col, matrix2[row * m1M + col]);
#endif

      index++;
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);

  struct timeval tval_before, tval_after, tval_result;
  if (my_rank == 0)
  {
    gettimeofday(&tval_before, NULL);
  }

  // Beginning of Cannon's

  // 1. Skew A (to the left)
  // (i,j) --> (i,(i-j)%dim)
  if (myi)
  {

    MPI_Request send_request, recv_request;

    int dest = myj - myi;
    if (dest < 0)
    {
      dest = dim + dest;
    }
    dest = NAME_TO_RANK(dim, dim, myi, dest);

#ifdef DEBUG
    fprintf(stderr, "[Skew A]: Block %d sending data to block %d\n", my_rank, dest);
#endif

    MPI_Isend(block_a, mySize, MPI_DOUBLE, dest, 0, MPI_COMM_WORLD, &send_request);
    MPI_Irecv(tmp, mySize, MPI_DOUBLE, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &recv_request);
    MPI_Wait(&send_request, &status);
    MPI_Wait(&recv_request, &status);
    memcpy(block_a, tmp, sizeof(double) * mySize);
  }

  // 2. Skew B (up)
  //  (i,j) --> (i-j,j);
  if (myj)
  {
    int dest = myi - myj;
    if (dest < 0)
    {
      dest = dim + dest;
    }
    dest = NAME_TO_RANK(dim, dim, dest, myj);
    MPI_Request send_request, recv_request;

#ifdef DEBUG
    fprintf(stderr, "[Skew B]: Block %d sending data to block %d\n", my_rank, dest);
#endif

    MPI_Isend(block_b, mySize, MPI_DOUBLE, dest, 2, MPI_COMM_WORLD, &send_request);
    MPI_Irecv(tmp, mySize, MPI_DOUBLE, MPI_ANY_SOURCE, 2, MPI_COMM_WORLD, &recv_request);
    MPI_Wait(&send_request, &status);
    MPI_Wait(&recv_request, &status);
    memcpy(block_b, tmp, sizeof(double) * mySize);
  }

  // In a Loop (multiply -> shift)
  for (int j = 1; j < dim; j++)
  {
    // do the multiplication
    int retval = mm(block_a, data_dim, data_dim, data_dim, block_b, block_c);
    if (retval)
    {
      fprintf(stderr, "There was an error in MM\n");
      return 4;
    }
    MPI_Request send_request, recv_request;

    //shift  (i,j) --> (i,(j-1)%dim)
    int dest = myj - 1;
    if (dest < 0)
    {
      dest = dim + dest;
    }
    dest = NAME_TO_RANK(dim, dim, myi, dest);

    MPI_Isend(block_a, mySize, MPI_DOUBLE, dest, j, MPI_COMM_WORLD, &send_request);
    MPI_Irecv(tmp, mySize, MPI_DOUBLE, MPI_ANY_SOURCE, j, MPI_COMM_WORLD, &recv_request);
    MPI_Wait(&send_request, &status);
    MPI_Wait(&recv_request, &status);
    memcpy(block_a, tmp, sizeof(double) * mySize);

    //shift  (i,j) --> ((i-1)%dim,j);
    dest = myi - 1;
    if (dest < 0)
    {
      dest = dim + dest;
    }

    dest = NAME_TO_RANK(dim, dim, dest, myj);
    MPI_Isend(block_b, mySize, MPI_DOUBLE, dest, j + world_size, MPI_COMM_WORLD, &send_request);
    MPI_Irecv(tmp, mySize, MPI_DOUBLE, MPI_ANY_SOURCE, j + world_size, MPI_COMM_WORLD, &recv_request);
    MPI_Wait(&send_request, &status);
    MPI_Wait(&recv_request, &status);
    memcpy(block_b, tmp, sizeof(double) * mySize);
  }

  // do the multiplication
  int retval = mm(block_a, data_dim, data_dim, data_dim, block_b, block_c);
  if (retval)
  {
    fprintf(stderr, "There was an error in MM\n");
    return 5;
  }

#ifdef DEBUG
  for (int i = 0; i < mySize; i++)
  {
    fprintf(stderr, "Cannon's done: [Rank %d] block_a[%d] =  %f, block_b[%d] = %f, block_c = %f\n", my_rank, i, block_a[i], i, block_b[i], block_c[i]);
  }
#endif
  // End of Cannon's

  // need to collect all of the data back to rank0 and put them together
  if (my_rank == 0)
  {
    for (int i = 0; i < data_dim; i++)
    {
      for (int j = 0; j < data_dim; j++)
      {
        ARRAY(C, m1M, m1M, i, j) = ARRAY(block_c, data_dim, data_dim, i, j);
      }
    }

    // I get to collect the answer
    for (int r = 1; r < world_size; r++)
    {
      MPI_Recv(block_c, mySize, MPI_DOUBLE, r, 0, MPI_COMM_WORLD, &status);
      for (int i = 0; i < data_dim; i++)
      {
        for (int j = 0; j < data_dim; j++)
        {
          ARRAY(C, m1M, m1M, ((r / dim) * data_dim + i), (r % dim) * data_dim + j) = ARRAY(block_c, data_dim, data_dim, i, j);
        }
      }
    }
  }
  else
  {
    MPI_Send(block_c, mySize, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
  }

  MPI_Barrier(MPI_COMM_WORLD);
  if (my_rank == 0)
  {
    gettimeofday(&tval_after, NULL);
    timersub(&tval_after, &tval_before, &tval_result);
    printf("Time elapsed: %ld.%06ld\n", (long int)tval_result.tv_sec, (long int)tval_result.tv_usec);
  }

  return 0;
}


// This function will write matrix C to the output file
int write_matrix(double C[], FILE *output, int m1M, int outputMSize, MM_typecode mm_matcode)
{

  mm_write_banner(output, mm_matcode);
  mm_write_mtx_array_size(output, m1M, m1M);
  for (int i = 0; i < outputMSize; i++) {
    fprintf(output, "%lf\n", C[i]);
  }

  fclose(output);

  return 0;
}


// This function performs serial multiplication of two row-major matrices
double *matrix_mul(int dim, double *a, double *b, double *c)
{

  for (int i = 0; i < dim; i++)
  {
    for (int j = 0; j < dim; j++)
    {
      for (int k = 0; k < dim; k++)
      {

        c[j + i * dim] += a[k + i * dim] * b[j + k * dim];
      }
    }
  }

  return c;
}
