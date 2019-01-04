/***************************************************/
/* Kompilieren: icc -o rueckei_seq rueckei_seq.c           */
/* Ausfuehren:  ./rueckei_seq                          */
/* Das Ergebnis lautet:                            */
/*   1.000000    2.000000    3.000000    4.000000  */
/***************************************************/
#include <stdlib.h>
#include <stdio.h>
#define mat_size 20
int main(int argc, char *argv[])
{
  int i, j;
    int world_rank, world_size;


  double matrix[mat_size][mat_size];
  double vectorB[mat_size];
  double vectorX[mat_size];

    for(int i = 0; i < mat_size; ++i)
    {
        for(int j = 0; j < mat_size; ++j)
        {
            if(j < i)
                matrix[i][j] = 0;
            else
                matrix[i][j] = mat_size-j;
        }
        b[i] = i+1;
    }

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);


  for(i=mat_size-1; i>=0; --i)
    {
        int root = i%world_size;
        if(root == world_rank) {
            /* Berechnung Element des Ergebnisvektors*/
            vectorX[i] = vectorB[i] / matrix[i][i];
            //TODO send
        }
        MPI_Bcast(&vectorX[i], 1, MPI_DOUBLE, root, MPI_COMM_WORLD);
          /* Update mit gerade berechnetem Ergebnis */
          /* Update every worldrank+k*worldsize spalte*/
          for(j=world_rank; j<mat_size; j+=world_size)
            vectorB[j] -= matrix[j][i]*vectorX[i];
    }


  /* Ausgabe */
    if(world_rank == 0) {
        printf("Das Ergebnis lautet:\n");
        for(i=0;i<4;++i)
            printf("  %f  ",vectorX[i]);
        printf("\n");
    }
    /********** Finalize MPI **********/
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
    return(0);
}
