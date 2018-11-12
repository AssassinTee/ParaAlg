/* Sequential program */
/* Compilieren mit icc -o mergeseq mergeseq.c */
/* Ausfuehren mit ./mergeseq  */
/* Sortiert n integer Werte, n<100000 */

#include <stdio.h>
#include <stdlib.h>

#define MAXNUM 100000

  void mergesort(int* iein, int istart, int iend);
  void merge(int* iein, int istart, int iend);

int main(int argc, char** argv)
{
  int n;   /* Anzahl der zu sortierenden Zahlen */
  int ieing[MAXNUM];
  int i,j,k, istart, iend, ihalf;

  FILE *fp, *aus;


  n=13;
  for (i=0;i<n;++i) {
    ieing[i]=rand();
    printf("%d ",ieing[i]);
  }
    printf("\n ");
  istart=0;
  iend=n-1;
  mergesort(ieing,istart,iend);
  /* aus = fopen ("Ausgabe", "w");
  fprintf(aus,"Es waren %d Elemente zu sortieren, Ergebnis:\n", n); */
  printf("Es waren %d Elemente zu sortieren, Ergebnis:\n", n);
  for(i=0;i<n;++i) {
    /* fprintf(aus,"%d \n",ieing[i]); */
    printf("%d ",ieing[i]);
  }
    printf("\n ");
  /* fclose(aus); */
}

void mergesort(int* iein, int istart, int iend){
  if(iend>istart){
    if(iend>istart+1){
      mergesort(iein,istart,(istart+iend-1)/2);
      mergesort(iein,(istart+iend-1)/2+1,iend);
    }
    merge(iein,istart,iend);
  }
}


void merge(int* iein, int istart, int iend){
  int ihalf, i ,j , k, l;
  int *tempa;

  tempa = (int *) malloc((iend-istart+1) * sizeof(int));
  ihalf=(iend+istart-1)/2;
  i=istart;
  j=ihalf+1;
  k=0;
  do{
    if(iein[i]<iein[j]){
      tempa[k]=iein[i];
      i=i+1;
    }
    else{
      tempa[k]=iein[j];
      j=j+1;
    }
    k=k+1;
  } while ( (i<=ihalf) && (j<=iend));
  if(i>ihalf){
    for(l=0;l<=iend-j;++l){
      tempa[k+l]=iein[j+l];
    }
  }
  else{
    for(l=0;l<=ihalf-i;++l){
      tempa[k+l]=iein[i+l];
    }
  }
  for(l=0;l<=iend-istart;++l){
    iein[l+istart]=tempa[l];
  }
  free(tempa);
}
