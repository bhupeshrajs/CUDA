#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>

#define Tolerance 0.00001
#define TRUE 1
#define FALSE 0

void serialSolver(double **A, double **temp, int n)
{
    int convergence = FALSE;
    double diff;
    int i,j, iters=0;

    do
    {
        diff = 0.0;

        for (i = 1 ; i < n ; i++ )
        {
            for (j = 1 ; j < n ; j++ )
            {
                temp[i][j] = 0.2*(A[i][j] + A[i][j-1] + A[i-1][j] + A[i][j+1] + A[i+1][j]);
                diff += fabs(temp[i][j] - A[i][j]);
            }
        }
        iters++;

        for(i = 1 ; i < n ; i++ )
        {
            for(j = 1 ; j < n ; j++ )
            {
                A[i][j] = temp[i][j];
            }
        }
       
        if (diff/((double)n*(double)n) < Tolerance) {
            convergence=TRUE;
            printf("\n Convergence at iteration : %d and difference is : %f Serial Execution",iters,diff);
        }
       
        //printf("\nDifference is %f at iteration %d",diff,iters);
    } while(convergence == FALSE ); /*for*/
    
}