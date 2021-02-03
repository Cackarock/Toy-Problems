//////////////////////////////////////////////////////////////////////
//Authors: Kenslee Moy, Ravi Shankar, Dylan Leman, Floriana Ciaglia //
//////////////////////////////////////////////////////////////////////

#include"pi-mc.h"
#include<stdio.h>
#include<ctype.h>
#include<string.h>
#include<math.h>
#include<stdlib.h>
#include<errno.h>
#include<time.h>
#include<mpi.h>

/** 
 * Prints program usage
 */ 
void print_usage(void){
  printf("\nUsage:");
  printf("\n      ./pi-mc <positve integer>\n\n");
  return;
}

/** 
 * Check that there is only 1 command line argument
 */ 
int check_num_args(int argc){
  if(argc != 2){ 
    return -1; 
  }
  return 0;
}

/** 
 * Check that the command line arg was an integer
 */ 
int check_args(char *argv[]){
  printf("Parsing '%s':\n", argv[1]);
  // errno can be set to any non-zero value by a library function call
  // regardless of whether there was an error, so it needs to be cleared
  // in order to check the error set by strtol
  errno = 0;
  char *end;
   
  int i = strtol(argv[1], &end, 10);
  if (argv[1] != end && *end =='\0' && i >= 0 ){
    return i;
  }
  return -1; 
}


/** 
 * Estimate pi to n iterations using monte carlo method
 */ 
double s_monte_carlo(int n_iter){

  double x,y,z,pi;
  int i;
  int count = 0;
  srand(time(NULL));
  
  for (i=0; i<n_iter; ++i){
    //get random points
    x = (double)random()/RAND_MAX;
    y = (double)random()/RAND_MAX;
    z = sqrt((x*x)+(y*y));
    //check to see if point is in unit circle
    if (z<=1){
        ++count;
    }
  }
  pi = ((double)count/(double)n_iter)*4.0;          //p = 4(m/n)
  return pi;
}

//distribute tasks to processes
void pi_proc(int n, int id, int numprocs){
    int partition = n / numprocs, remainder = n % numprocs;
    long double pi = 0;
    if(id != numprocs -1){
        pi = p_monte_carlo(partition, id);
    }else{
        pi = p_monte_carlo(partition + remainder, id);
    }   

    //printf("%d: %Lf\n", id, pi);
    
    long double gpi;
    MPI_Reduce(&pi, &gpi, 1, MPI_LONG_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    
    if(id == 0){ 
        gpi = (long double) gpi / numprocs;
        printf("%Lf\n", gpi);
    }   
} 

//monte carlo method of calculating pi, id must be specified for each procs seed to be different
long double p_monte_carlo(int n, int id){
    int circle_points = 0, square_points = 0, i;
        long double x, y, d, pi; 

        time_t t;

        srand( (unsigned) time(&t) + id);
    
        for(i = 0; i < n; i++){
                x = (long double) (rand() % (n +1)) / n;
                y = (long double) (rand() % (n +1)) / n;
                d = (x * x) + (y * y); 

                if(d <= 1) circle_points++;

                square_points++;
        }   

        pi = (long double) 4 * circle_points / square_points;
    
    return pi;  
}
