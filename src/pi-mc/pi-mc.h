#ifndef MONTE_CARLO_H 
#define MONTE_CARLO_H


void print_usage(void);
int check_args(char *argv[]);
int check_num_args(int);

//Serial
double s_monte_carlo(int);

//Parallel
void pi_proc(int n, int id, int numprocs);
long double p_monte_carlo(int n, int id);

#endif
