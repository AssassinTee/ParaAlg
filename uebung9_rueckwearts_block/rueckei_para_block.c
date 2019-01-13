/***************************************************/
/* Kompilieren: icc -o rueckei_seq rueckei_seq.c           */
/* Ausfuehren:  ./rueckei_seq                          */
/* Das Ergebnis lautet:                            */
/*   1.000000    2.000000    3.000000    4.000000  */
/***************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#define mat_size 16
#define nb 2
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

    if(mat_size%world_size != 0)
    {
        MPI_Abort(MPI_COMM_WORLD, 42);
        return 1;
    }

    int blocks = mat_size/nb;
    int local_blocks = blocks/world_size;
    double matrix[local_blocks*nb][mat_size];
    double vectorX[mat_size];
    double vectorB[local_blocks*nb];

    for(int i = 0; i < local_blocks; ++i)//iterate over blocks
    {
        for(int k = 0; k < nb; ++k)  //iteriere über zeilen im block
        {
            int i2 = i*blocks+k;//globaler index i2, i-ter block, k-te zeile
            int i3 = i*nb+k;//localer index i3, i-ter block, k-te spalte
            for(int j = 0; j < mat_size; ++j)//iteriere über spalten
            {
                if(j < i2)
                    matrix[i3][j] = 0;
                else
                    matrix[i3][j] = mat_size-j;
            }
            vectorB[i3] = i2+1;
        }
    }

    //MPI_Init(&argc, &argv);
    //MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    //MPI_Comm_size(MPI_COMM_WORLD, &world_size);


    for(i=blocks-1; i>=0; --i)
    {
        int root = i%world_size;//root ist block%world_size
        int i2 = i/world_size;//i2 ist der i-te block vom root prozess
        if(root == world_rank)
        {
            /* Berechnung Block des Ergebnisvektors*/
            for(int k=nb-1; k>= 0; --k)
            {
                vectorX[i+k] = vectorB[i2] / matrix[i2][i+k];
                //Update
                for(j = k-1; j >= 0; --j)
                    vectorB[j] -= matrix[j][i+k]*vectorX[i+k];
            }
            //TODO send
        }
        MPI_Bcast(&vectorX[i], nb, MPI_DOUBLE, root, MPI_COMM_WORLD);
        /* Update mit gerade berechnetem Ergebnis */
        for(j=0; j < mat_size/world_size; ++j)//matrix vector mult
        {
            for(int k = nb-1; k >= 0; --k)
            {
                vectorB[j+k] -= matrix[j][i]*vectorX[i];
            }
        }
    }


    /* Ausgabe */
    if(world_rank == 0)
    {
        printf("Das Ergebnis lautet:\n");
        for(i=0; i<mat_size; ++i)
            printf("  %f  ",vectorX[i]);
        printf("\n");
    }
    /********** Finalize MPI **********/
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    return(0);
}
