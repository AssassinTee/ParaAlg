/*  Kompilieren mit icc -qopenmp testdgemv_neu.c -o testdgemv_neu_c -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -liomp5 -lpthread -lm
*/

#include <stdlib.h>
#include <stdio.h>
#include <mkl.h>
#include <mpi.h>
#include <math.h>

int main(int argc, char** argv) {
	MPI_Init(&argc, &argv);
    int n, i, j, k;
	double *x, *y, *A, *yn;
    double alpha=1.0;
	double beta = 1.0;
    double one = 1.0;
    double errmax=0.0, tempa;
    double time;
	FILE *fp, *aus;
    int errcnt, incx = 1, incy = 1, lda;
	int nstart =500, nend = 4000, nincr = 500;
	int maxiter=1;
	int size, my_rank;
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	
	if(my_rank == 0){
		fp = fopen ("ingemv", "r");
		fscanf(fp," %d %d %d \n",&nstart, &nend, &nincr);
		fscanf(fp," %d \n",&maxiter);
		fscanf(fp," %lf %lf \n",&alpha,&beta);
		fclose(fp);
		/*	int nstart = 500, nend = 4000, nincr = 500; */
		
		
		aus = fopen("out_dgemv_c_mpi","aw");	
	}
	
	MPI_Bcast(&nstart, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&nend, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&nincr, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&maxiter, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&alpha, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&beta, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	/* Loop over different problem sizes */

	
	
	for(n = nstart; n <= nend; n+= nincr) {
	  x = malloc(sizeof(double) * n);
	  y = malloc(sizeof(double) * n/size);
	  yn = malloc(sizeof(double) * n);
	  A = malloc(sizeof(double) * n*n/size);
	  lda=n;
	   
	   
	  /*  Put values to vectors which allow to control results */
	   
	  for(i = 0; i < n; ++i) {
	     x[i] = 1.0;
	  }
	  for(i = 0; i < n/size; ++i){
		 y[i] = (n-1)*1.0;
	  }
	  for(i = 0; i < n/size; ++i) {
	    for(j = 0; j < n; ++j) {
	      A[i*n + j] = (my_rank * n * n/size + i*1.0*n + (j+1.0)) ;
	    }
	  }
		if(my_rank == 0){
			time = MPI_Wtime();
		}
		
		
		
		cblas_dgemv(CblasRowMajor, CblasNoTrans, n/size, n, alpha, A, lda, x, incx, beta, y, incy);	
		
		MPI_Gather(y, n/size, MPI_DOUBLE, yn, n/size, MPI_DOUBLE, 0, MPI_COMM_WORLD); 
		 
		if(my_rank == 0){ 
		  time = (MPI_Wtime()-time)/maxiter;
		  fprintf(aus,"Benoetigte Zeit fuer (%d X %d)-Matrix: %f\n", n, n, time);
	 
		  errmax=0.0;
		  errcnt=0;
		  if (maxiter==1){
			for (i = 0; i < n; i++){
				//fprintf(aus, "yn %f\n", yn[i]);
				tempa = yn[i]-(alpha*(2.0*i+1.0)*n*n+(alpha+2.0*beta)*n-2.0*beta)/2.0;
				/*fprintf(aus,"y[%d] = %f, sollte sein= %f \n",i,y[i],(alpha*(2.0*i+1.0)*n*n+(alpha+2.0*beta)*n-2.0*beta)/2.0); */
				if (fabs(tempa)>errmax){
				  errmax=fabs(tempa);
				  errcnt=errcnt+1;
				}
			 }
			fprintf(aus,"Maximaler absoluter Fehler: %f \n",errmax);
		  }
	
		fflush(aus);
		}
		
		free(A); free(x); free(y); free(yn);
	 	  
	}
	
	if(my_rank == 0){
		fclose(aus);
	}
	MPI_Finalize();
	return 0;
}

