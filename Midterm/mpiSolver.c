#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <mpi.h>

#define Tolerance 0.00001
#define TRUE 1
#define FALSE 0

#define N 5000

double **A;
double **temp;

double **serial_A;
double **serial_temp;

extern void serialSolver(double **A, double** temp, int n);

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

void assignRow(double **A, double *temp_row , int row , int n) {
    
    int i;
    for( i = 1 ; i < n ; i++ ) {
        temp_row[i] = A[row][i];
    }
    
}

void copyData(double **A, double *result , int startRow, int endRow , int n ) {
    
    int i,j;
     for ( i = startRow ; i <= endRow ; i++ )
    {
        for ( j = 1 ; j < n ; j++ )
        {
            result[ (i-startRow)*(n+2) + j ] = A[i][j];
        }
    }
    
}

void copyFromData(double **A, double *result , int startRow, int endRow , int n ) {
    
    int i,j;
     for ( i = startRow ; i <= endRow ; i++ )
    {
        for ( j = 1 ; j < n ; j++ )
        {
            A[i][j] = result[ (i-startRow)*(n+2) +j];
        }
    }
    
}

void assignFromRow(double **A, double *temp_row , int row , int n) {
    
    int i;
    for( i = 1 ; i < n ; i++ ) {
        A[row][i] = temp_row[i];
    }
    
}

long usecs (void)
{
    struct timeval t;

    gettimeofday(&t,NULL);
    return t.tv_sec*1000000+t.tv_usec;
}

void print(double**A , int n ) {
    int i,j;
    for( i = 0 ; i < n+2 ; i++ ) {
        for( j = 0 ; j < n+2 ; j++ ) {
            printf(" %f",A[i][j]);
        }
        printf("\n");
    }
}
int check(double **A, double **B , int n) {
    
    int i,j;
    int correct = 1;
    for( i = 1 ; i <= n ; i++ ) 
    {
        for( j = 1 ; j <= n ; j++ ) 
        {
            if( A[i][j] != B[i][j] ) {
                correct = 0;
                break;
            }
        }
    }
    return correct;
}

int main()
{
    // Initialize the MPI environment
    MPI_Init(NULL, NULL);
    // Get the number of processes
    int totalRank;
    MPI_Comm_size(MPI_COMM_WORLD, &totalRank);
    // Get the rank of the process
    int myRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);    
    
    
    int n = N;
    
    A = malloc((n+2) * sizeof(double *));
    temp = malloc((n+2) * sizeof(double *));
    
    serial_A = malloc((n+2) * sizeof(double *));
    serial_temp = malloc((n+2) * sizeof(double *));
    for (int i = 0 ; i < n+2 ; i++)
    {
        A[i] = malloc((n+2) * sizeof(double));
        temp[i] = malloc((n+2) * sizeof(double));    
        serial_A[i] = malloc((n+2) * sizeof(double));
        serial_temp[i] = malloc((n+2) * sizeof(double)); 
    }
    
    initialize(A,n);
    initialize(serial_A,n);
    
    int workTag = 1;
    int dataTag = 2;
    int stopTag = 3;
    int workDoneTag = 4;
    
    double diff = 0.0;
    double global_diff = 0.0;
    
    double serial_start_time,serial_end_time,serial_time;
    double mpi_start_time,mpi_end_time,mpi_time;
    double communication_time = 0.0;
    double communication_start_time,communication_end_time;
    double global_communication_time = 0.0;
    
    
    if( myRank == 0 ) {
        
        mpi_start_time=usecs();
        
        int iters = 0;
        int reply;
        int i;
        int convergence = FALSE;
        while (convergence == FALSE) {
            
            iters++;
            global_diff = 0.0;
            
            
            communication_start_time = usecs();
            for( i = 1 ; i < totalRank ; i++ ) {
                MPI_Send(&reply,1,MPI_INT,i,workTag,MPI_COMM_WORLD);    
            }
            communication_end_time = usecs();
            communication_time += (communication_end_time - communication_start_time);
            
            
            double difference = 0.0;
            MPI_Reduce(&difference,&global_diff,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
            //printf("\nThe sum diff at iteration %d is : %f",iters,global_diff);
            
            if (global_diff/((double)n*(double)n) < Tolerance) {
                convergence=TRUE;
                printf("\n Convergence at iteration : %d and difference is : %f Parallel Execution",iters,global_diff);
            }
            
        }
        
        
        communication_start_time = usecs();
        for( i = 1 ; i < totalRank ; i++ ) {
                MPI_Send(&reply,1,MPI_INT,i,stopTag,MPI_COMM_WORLD);    
        }
        communication_end_time = usecs();
        communication_time += (communication_end_time - communication_start_time);
        //printf("\nConvergence at iteration : %d and difference is : %f",iters,global_diff);
        
        for( i = 1 ; i < totalRank ; i++ ) {
            
            int success;
            int startRow;
            int endRow;
            
            MPI_Status status;
           
            communication_start_time = usecs();  
            
            MPI_Recv(&success,1,MPI_INT,MPI_ANY_SOURCE,workDoneTag,MPI_COMM_WORLD,&status);
            int id = status.MPI_SOURCE;
            MPI_Recv(&startRow,1,MPI_INT,id,dataTag,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            MPI_Recv(&endRow,1,MPI_INT,id,dataTag,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            
            communication_end_time = usecs();
            communication_time += (communication_end_time - communication_start_time);
            
            
            int sendRows = endRow - startRow + 1;
            
            double *result =  malloc((sendRows) *(n+2) * sizeof(double));
    
            communication_start_time = usecs();
            MPI_Recv(&(result[0]),(n+2)*sendRows,MPI_DOUBLE,id,dataTag,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            communication_end_time = usecs();
            communication_time += (communication_end_time - communication_start_time);
            
            
            copyFromData(A,result,startRow,endRow,n);
            
        }
        
        mpi_end_time = usecs();
        
        serial_start_time=usecs();
        serialSolver(serial_A,serial_temp,n);
        serial_end_time=usecs();
        
        serial_time = ((double)(serial_end_time-serial_start_time))/1000000;
        mpi_time = ((double)(mpi_end_time-mpi_start_time))/1000000;
        
        printf("\n Serial Computation time = %fs\n", serial_time);
        printf("\n MPI Computation time = %fs\n", mpi_time);
    
        MPI_Reduce(&communication_time,&global_communication_time,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
        
        printf("\n MPI Communication time = %fs\n", global_communication_time/1000000);
        
        
        //print(serial_A,n);
        //printf("\n The parallel execution \n");
        //print(A,n);
        if( check(serial_A,A,n) != 1 ) {
            printf("\n The output doesn't match with serial execution ");
        }
        else {
            printf("\n The serial output matches with the MPI execution");
        }
    
        printf("\n Speedup Achieved is : %fx with %d threads\n\n",(serial_time/mpi_time),totalRank-1);

        
        
    }
    else {
        
        int numRowsPerProcess = n/(totalRank-1);
        
        int startRow = (myRank-1)*numRowsPerProcess + 1;
        
        int startGhostRow = (myRank-1)*numRowsPerProcess;
        int endGhostRow;
        
        if( myRank == totalRank - 1 ) {
            endGhostRow = n;
        }
        else 
        {
            endGhostRow = (myRank)*numRowsPerProcess +1;
        }
        
        int endRow = endGhostRow-1;
        
        int totalRows = endGhostRow - startGhostRow + 1;
        
        int iter = 0;
        
        //printf("\nHello from process %d with startRow %d and endRow %d",myRank,startRow,endRow);
        
        int reply;
        MPI_Status status;
        MPI_Request request;
        
        double *temp_row = malloc((n+2) * sizeof(double));
        
        while ( ((MPI_Recv(&reply, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD,&status)) == MPI_SUCCESS) &&
            (status.MPI_TAG == workTag) ) 
        {
            
            communication_start_time = usecs();
            
            
            if( myRank != 1 && myRank != totalRank - 1 ) {
                assignRow(A,temp_row,startRow,n);
                //printf("Sending message from process %d to %d and %d",myRank,myRank-1,myRank+1);
                MPI_Isend(temp_row,n+2,MPI_DOUBLE,myRank-1,dataTag,MPI_COMM_WORLD,&request);
                
                assignRow(A,temp_row,endRow,n);
                MPI_Isend(temp_row,n+2,MPI_DOUBLE,myRank+1,dataTag,MPI_COMM_WORLD,&request);
            }
            else if (myRank == 1 ) {
                //printf("Sending message from process %d to %d",myRank,myRank+1);
                assignRow(A,temp_row,endRow,n);
                MPI_Isend(temp_row,n+2,MPI_DOUBLE,myRank+1,dataTag,MPI_COMM_WORLD,&request);
            }
            else if( myRank == totalRank - 1) {
                //printf("Sending message from process %d to %d",myRank,myRank-1);
                assignRow(A,temp_row,startRow,n);
                MPI_Isend(temp_row,n+2,MPI_DOUBLE,myRank-1,dataTag,MPI_COMM_WORLD,&request);
            }
            
            
            if( myRank != 1 && myRank != totalRank - 1 ) {
                
                MPI_Recv(temp_row,n+2,MPI_DOUBLE,myRank-1,dataTag,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
                assignFromRow(A,temp_row,startGhostRow,n);
                
                MPI_Recv(temp_row,n+2,MPI_DOUBLE,myRank+1,dataTag,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
                assignFromRow(A,temp_row,endGhostRow,n);
            }
            else if (myRank == 1 ) {
                MPI_Recv(temp_row,n+2,MPI_DOUBLE,myRank+1,dataTag,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
                assignFromRow(A,temp_row,endGhostRow,n);
            }
            else if( myRank == totalRank - 1) {
                MPI_Recv(temp_row,n+2,MPI_DOUBLE,myRank-1,dataTag,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
                assignFromRow(A,temp_row,startGhostRow,n);
            }
            
            communication_end_time = usecs();
            communication_time += (communication_end_time - communication_start_time);
            
            
            iter++;
            diff = 0.0;
        
            for (int i = startRow ; i <= endRow ; i++ )
            {
                for (int j = 1 ; j < n ; j++ )
                {
                    temp[i][j] = 0.2*(A[i][j] + A[i][j-1] + A[i-1][j] + A[i][j+1] + A[i+1][j]);
                    diff += fabs(temp[i][j] - A[i][j]);
                }
            }
            
            MPI_Reduce(&diff,&global_diff,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
            
            for (int i = startRow ; i <= endRow ; i++ )
            {
                for (int j = 1 ; j < n ; j++ )
                {
                    A[i][j] = temp[i][j];
                }
            }
        
            //printf("\nThe sum diff at iteration %d is : %f",iter,diff);
            
        }
        
        
        int success;
        
        communication_start_time = usecs();
        
        MPI_Send(&success,1,MPI_INT,0,workDoneTag,MPI_COMM_WORLD);
        MPI_Send(&startRow,1,MPI_INT,0,dataTag,MPI_COMM_WORLD);
        MPI_Send(&endRow,1,MPI_INT,0,dataTag,MPI_COMM_WORLD);
        
        communication_end_time = usecs();
        communication_time += (communication_end_time - communication_start_time);
            
        int sendRows = endRow - startRow + 1;
        
        double* result;
        result =  malloc( (sendRows) * (n+2) * sizeof(double));
    
        copyData(A,result,startRow,endRow,n);
        
        communication_start_time = usecs();
        
        MPI_Send(&(result[0]),sendRows*(n+2),MPI_DOUBLE,0,dataTag,MPI_COMM_WORLD);
        communication_end_time = usecs();
        communication_time += (communication_end_time - communication_start_time);
            
        MPI_Reduce(&communication_time,&global_communication_time,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
        
    }
    
    MPI_Finalize();
    
    
}