#ifndef MATRIX_MULT_H 
#define MATRIX_MULT_H
#include "../../include/mmio.h"

void print_usage(void);
int mpi_matrix_mult(char *argv[], int, int);
int read_matrix(char* argv[], int, int);
int cannon(double*, double*, double*, int, int, int);
double* matrix_mul(int dimensions, double* a, double* b, double* c);
int mm(double* A,int c,int rc,int r,double* B,double* C);
int write_matrix(double C[], FILE *output, int m1M, int outputMSize, MM_typecode mm_typcode);

#define EPS .0000001
#define COMPARE(V1,V2) V1 < V2 + EPS && V1 > V2 - EPS ? 0 : 1 
#define ARRAY(D,c,r,i,j) (D)[(i)*(r) + (j)]
#define ARRAYl(D,c,r,i,j) &((D)[(i)*(r) + (j)])
#define NAME_TO_RANK(c,r,i,j) (c)*(i) + (j) 

#endif
