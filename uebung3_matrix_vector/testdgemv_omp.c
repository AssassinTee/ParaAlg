/*  Kompilieren mit icc -qopenmp testdgemv_neu.c -o testdgemv_neu_c -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -liomp5 -lpthread -lm
*/

#include <stdlib.h>
#include <stdio.h>
#include <mkl.h>
#include <omp.h>
#include <math.h>

int main(int argc, char** argv) {
    int n, i, j, k;
	double *x, *y, *A;
    double alpha=1.0;
	double beta = 1.0;
    double one = 1.0;
    double errmax=0.0, tempa;
    double time;
	FILE *fp, *aus;
    int errcnt, incx = 1, incy = 1, lda;
	int nstart =500, nend = 4000, nincr = 500;
	int maxiter=1;
	
	fp = fopen ("ingemv", "r");
	fscanf(fp," %d %d %d \n",&nstart, &nend, &nincr);
	fscanf(fp," %d \n",&maxiter);
	fscanf(fp," %lf %lf \n",&alpha,&beta);
        fclose(fp);
	/*	int nstart = 500, nend = 4000, nincr = 500; */
	
	/* Loop over different problem sizes */

	aus = fopen("out_dgemv_c_omp","aw");
	
	for(n = nstart; n <= nend; n+= nincr) {
	  fprintf(aus,"(%d X %d)-Matrix \n", n, n);
	  x = malloc(sizeof(double) * n);
	  y = malloc(sizeof(double) * n);
	  A = malloc(sizeof(double) * n*n);
	  lda=n;
	   
	   
	  /*  Put values to vectors which allow to control results */
	   
	  for(i = 0; i < n; ++i) {
	    x[i] = 1.0;
	    y[i] = (n-1)*1.0;
	  }
	  for(i = 0; i < n; ++i) {
	    for(j = 0; j < n; ++j) {
	      A[i*n + j] = (i*1.0*n + (j+1.0)) ;
	    }
	  }
		
	  

	int thr[6] = {1, 2, 4, 8, 16, 32};  

	for(int t = 0; t < 6; t++){
		omp_set_num_threads(thr[t]);
		/*for(i = 0; i < n; ++i) {
			y[i] = (n-1)*1.0;
		}*/
		
		time = omp_get_wtime();
		
		#pragma omp parallel
		{
			int numThreads = omp_get_num_threads();
			int num = omp_get_thread_num();
			if(num == 0){
				fprintf(aus, "numThreads %d \n", omp_get_num_threads());
			}
			double * newA = A + num * n * n/numThreads;
			double * newY = y + num * n/numThreads;
			cblas_dgemv(CblasRowMajor, CblasNoTrans, n/numThreads, n, alpha, newA, lda, x, incx, beta, newY, incy);	
		}
		 
		  time = (omp_get_wtime()-time)/maxiter;
		  fprintf(aus,"Benoetigte Zeit fuer (%d X %d)-Matrix: %f\n", n, n, time);
	 
		  errmax=0.0;
		  errcnt=0;
		  if (maxiter==1){
			for (i = 0; i < n; i++){
				tempa = y[i]-(alpha*(2.0*i+1.0)*n*n+(alpha+2.0*beta)*n-2.0*beta)/2.0;
				/*fprintf(aus,"y[%d] = %f, sollte sein= %f \n",i,y[i],(alpha*(2.0*i+1.0)*n*n+(alpha+2.0*beta)*n-2.0*beta)/2.0); */
				if (fabs(tempa)>errmax){
				  errmax=fabs(tempa);
				  errcnt=errcnt+1;
				}
			 }
			fprintf(aus,"Maximaler absoluter Fehler: %f \n",errmax);
		  }
	}
	 	fflush(aus);
		free(A); free(x); free(y);
	 	  
	}
	
	
	fclose(aus);
	return 0;
}

