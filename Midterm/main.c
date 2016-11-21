#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <omp.h>

#define Tolerance 0.00001
#define TRUE 1
#define FALSE 0

#define N 5000

double ** serial_A;
double ** openMP_A;
double **temp;

int num_threads;

extern void serialSolver(double **A, double** temp, int n);
extern void openMPSolver(double **A, double **temp, int n, int num_threads);

int initialize (double **A, int n)
{
    int i,j;

    for ( j = 0 ; j < n+1 ; j++ )
    {
        A[0][j] = 1.0;
    }
    
    for ( i = 1 ; i < n+1 ; i++ )
    {
        A[i][0] = 1.0;
        for ( j = 1 ; j < n+1 ; j++ ) A[i][j]=0.0;
    }
    
    return 0;

}

long usecs (void)
{
    struct timeval t;

    gettimeofday(&t,NULL);
    return t.tv_sec*1000000+t.tv_usec;
}

int check(double **A, double **B , int n) {
    
    int i,j;
    int correct = 1;
    for( i = 1 ; i < n ; i++ ) 
    {
        for( j = 1 ; j < n ; j++ ) 
        {
            if( A[i][j] != B[i][j] ) {
                correct = 0;
                break;
            }
        }
    }
    return correct;
}

int main(int argc, char * argv[])
{
    int i;
    long serial_t_start,serial_t_end;
    long openMP_t_start,openMP_t_end;
    double serial_time,openMP_time;

    if( argc < 2 ) {
        printf("Give number of threads arguments");
        return 0;
    }
    num_threads = atoi(argv[1]);
    
    serial_A = malloc((N+2) * sizeof(double *));
    openMP_A = malloc((N+2) * sizeof(double *));
    temp = malloc((N+2) * sizeof(double *));
    for (i = 0 ; i < N+2 ; i++)
    {
        serial_A[i] = malloc((N+2) * sizeof(double));
        openMP_A[i] = malloc((N+2) * sizeof(double));
        temp[i] = malloc((N+2) * sizeof(double));    
    }

    initialize(serial_A, N);
    initialize(openMP_A,N);
    printf("\nInitialization Done");
    
    serial_t_start=usecs();
    serialSolver(serial_A,temp,N);
    serial_t_end=usecs();
    
    openMP_t_start=usecs();
    openMPSolver(openMP_A,temp,N,num_threads);
    openMP_t_end=usecs();

    printf("\nSolving Done..");

    serial_time = ((double)(serial_t_end-serial_t_start))/1000000;
    openMP_time = ((double)(openMP_t_end-openMP_t_start))/1000000;
    
    printf("\nSerial Computation time = %fs\n", serial_time);
    printf("\nOpenMP Computation time = %fs\n", openMP_time);
    
    if( check(serial_A,openMP_A,N) != 1 ) {
        printf("\n The output doesn't match with serial execution ");
    }
    else {
        printf("\n The serial output matches with the openMP execution");
    }
    
    printf("\n Speedup Achieved is : %fx with %d threads",(serial_time/openMP_time),num_threads);
    return 0;

}