#include <mpi.h>
#include <stdio.h>
#include <math.h>

/* define f(x) */
double f(double);
double f(double x){
	return ( 4.0/ (1.0+x*x));
}

int main(int argc, char *argv[]){
	
    int done=0, myid, numprocs, i, namelen;
	/* define a const PI for test */
	double PI = 3.141592653589793238462643;
	double mypi, pi, h, sum, x;
	double startwtime=0.0, endwtime;

    char processor_name[MPI_MAX_PROCESSOR_NAME];
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD,&numprocs); 
    MPI_Comm_rank(MPI_COMM_WORLD,&myid); 
    MPI_Get_processor_name(processor_name,&namelen);
    fprintf(stdout, "Process %d of %d on %s\n", myid, numprocs, processor_name);
    int n;
	if (myid==0){
		printf("Please give N=\n"); 
        scanf("%d", &n);
        startwtime = MPI_Wtime(); 
    }
    /* broadcast the n*/
	MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    /* Get the width of the rectangle. All rectangles have the same width */
    h = 1.0 / (double) n;
	/*Assign initial value to rectangle area*/
    sum = 0.0; 
	for (i = myid + 1; i <= n; i += numprocs){
    /* Each process calculates the area of a part of the rectangle. If the total number of processes numprocs is 4, the 0-1 interval is divided into 100 rectangles
    0 proc 1 5 9 13 ... 97 
    1 proc 2 6 10 14 ... 98 
    2 proc 3 7 11 15 ... 99 
    3 proc 4 8 12 16 ... 100 */

		x = h * ((double)i - 0.5); 
        sum += f(x);
        }
    mypi = h * sum;/* The partial sum of each process is calculated in parallel */
    MPI_Reduce(&mypi, &pi, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    /* The area of all rectangles is obtained by adding the partial sum, which is the approximate π value */ 
    if (myid == 0) {
    /* The approximate value is printed out by executing the accumulated process No.0 */
    printf("pi is approximately %.16f, Error is %.16f\n", pi, fabs(pi - PI));
    endwtime = MPI_Wtime();
    printf("wall clock time = %f\n", endwtime-startwtime);
    fflush( stdout );
    }

    MPI_Finalize();

    }
