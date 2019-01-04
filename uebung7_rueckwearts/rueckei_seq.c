/***************************************************/
/* Kompilieren: icc -o rueckei_seq rueckei_seq.c           */
/* Ausfuehren:  ./rueckei_seq                          */
/* Das Ergebnis lautet:                            */
/*   1.000000    2.000000    3.000000    4.000000  */
/***************************************************/
#include <stdlib.h>
#include <stdio.h>

int main(int argc, char *argv[])
{
  int i, j;
  
  double matrix[4][4]={{4,3,2,1},{0,3,2,1},{0,0,2,1},{0,0,0,1}};
  double vectorB[4]={20,16,10,4};
  double vectorX[4];


  for(i=3; i>=0; --i)
    {
      /* Berechnung Element des Ergebnisvektors*/
      vectorX[i] = vectorB[i] / matrix[i][i];
      
      /* Update mit gerade berechnetem Ergebnis */
      for(j=i-1; j>=0; --j)
      vectorB[j] -= matrix[j][i]*vectorX[i];
    }


  /* Ausgabe */
      printf("Das Ergebnis lautet:\n");
      for(i=0;i<4;++i)
	printf("  %f  ",vectorX[i]);
      printf("\n");

  return(0);
}
