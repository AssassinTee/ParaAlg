/***************************************************/
/* Kompilieren: icc -o rueckei_seq rueckei_seq.c           */
/* Ausfuehren:  ./rueckei_seq                          */
/* Das Ergebnis lautet:                            */
/*   1.000000    2.000000    3.000000    4.000000  */
/***************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#define mat_size 17
int main(int argc, char *argv[])
{
  int i, j;
    int world_rank, world_size;

  //double matrix[mat_size][mat_size];
  //double vectorB[mat_size];
  //double vectorX[mat_size];
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  if(mat_size%world_size != 0) {
    MPI_Abort(MPI_COMM_WORLD, 42);
    return 1;
  }

  int np = mat_size/world_size;
  double matrix[np][mat_size];
  double vectorX[mat_size];
  double vectorB[np];

    for(int i = 0; i < np; ++i)
    {
	int i2 = i*world_size+world_rank;
        for(int j = 0; j < mat_size; ++j)
        {
            if(j < i2)
                matrix[i][j] = 0;
            else
                matrix[i][j] = mat_size-j;
        }
        vectorB[i] = i2+1;
    }

    //MPI_Init(&argc, &argv);
    //MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    //MPI_Comm_size(MPI_COMM_WORLD, &world_size);


  for(i=mat_size-1; i>=0; --i)
    {
        int root = i%world_size;
	int i2 = i/world_size;
        if(root == world_rank) {
            /* Berechnung Element des Ergebnisvektors*/
            vectorX[i] = vectorB[i2] / matrix[i2][i];
            //TODO send
        }
        MPI_Bcast(&vectorX[i], 1, MPI_DOUBLE, root, MPI_COMM_WORLD);
          /* Update mit gerade berechnetem Ergebnis */
          /* Update every worldrank+k*worldsize spalte*/
          for(j=0; j < mat_size/world_size; ++j)
            vectorB[j] -= matrix[j][i]*vectorX[i];
    }


  /* Ausgabe */
    if(world_rank == 0) {
        printf("Das Ergebnis lautet:\n");
        for(i=0;i<mat_size;++i)
            printf("  %f  ",vectorX[i]);
        printf("\n");
    }
    /********** Finalize MPI **********/
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
    return(0);
}
