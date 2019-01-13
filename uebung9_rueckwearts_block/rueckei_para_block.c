/***************************************************/
/* Kompilieren: icc -o rueckei_seq rueckei_seq.c           */
/* Ausfuehren:  ./rueckei_seq                          */
/* Das Ergebnis lautet:                            */
/*   1.000000    2.000000    3.000000    4.000000  */
/***************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#define glob_mat_size 16
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

    if(glob_mat_size%world_size != 0)
    {
        MPI_Abort(MPI_COMM_WORLD, 42);
        return 1;
    }

    int glob_blocks = glob_mat_size/nb;
    int local_blocks = glob_blocks/world_size;
    int local_row_num = local_blocks*nb;
    double matrix[local_row_num][glob_mat_size];
    double vectorX[glob_mat_size];
    double vectorB[local_row_num];

    for(int i = 0; i < local_blocks; ++i)//iterate over blocks
    {
        for(int k = 0; k < nb; ++k)  //iteriere über zeilen im block
        {
            int glob_block = i*world_size+world_rank;

            int glob_index_row = nb*glob_block+k;//globaler zeilenindex, i-ter block, k-te zeile
            int loc_index_row = i*nb+k;//localer index i3, i-ter block, k-te spalte

            if(world_rank == 0) {
                printf("glob block: %d glob_index_row: %d loc_index_row: %d\n", glob_block, glob_index_row, loc_index_row);
            }
            for(int col = 0; col < glob_mat_size; ++col)//iteriere über spalten
            {
                if(col < glob_index_row)
                    matrix[loc_index_row][col] = 0;
                else
                    matrix[loc_index_row][col] = glob_mat_size-col;
            }
            vectorB[loc_index_row] = glob_index_row+1;
        }
    }

    for(int z = 0; z < glob_mat_size; ++z)
    {
        MPI_Barrier(MPI_COMM_WORLD);
        if((z/nb)%world_size == world_rank)
        {
            int offset = z%nb;
            int block = (z/nb)/world_size;
            int row = offset;
            printf("%d (B%d,Z%d):", world_rank, block, row);
            for(int v = 0; v < glob_mat_size; ++v)
            {
                printf("  %f  ",matrix[block*nb+row][v]);
            }
            printf("\n");
        }
    }

    for(i=glob_blocks-1; i>=0; --i)
    {
        int root = i%world_size;//root ist block%world_size
        int loc_block_index = i/world_size;//i2 ist der i-te block vom root prozess
        if(root == world_rank)
        {
            /* Berechnung Block des Ergebnisvektors*/
            for(int off=nb-1; off>= 0; --off)//offset ist zeilenindex des blocks
            {
                //Berechne indizies
                int glob_index_row = nb*i+off;
                int loc_index_row = loc_block_index+off;

                //Berechne Matrix vector
                vectorX[glob_index_row] = vectorB[loc_index_row] / matrix[loc_index_row][glob_index_row];

                //Update B Vektor NUR am lokalen Block
                for(int uoff = off-1; uoff >= 0; --uoff) {//uoff = update offset
                    int loc_index_row_update = loc_index_row+uoff;
                    vectorB[loc_index_row_update] -= matrix[loc_index_row_update][glob_index_row]*vectorX[glob_index_row];
                }
            }
        }
        /*Broadcast Block i (size=nb) an alle*/
        MPI_Bcast(&vectorX[i], nb, MPI_DOUBLE, root, MPI_COMM_WORLD);
        /* Update mit Block i */
        if(root != world_rank) {//Root hat schon geupdatet!!!
            for(j=0; j < local_row_num; ++j)//matrix vector mult
            {
                for(int off = nb-1; off >= 0; --off)//offset in block i
                {
                    int glob_index_row = nb*i+off;
                    int loc_index_row = j;
                    vectorB[loc_index_row] -= matrix[loc_index_row][glob_index_row]*vectorX[glob_index_row];
                }
            }
        }
    }


    /* Ausgabe */
    if(world_rank == 0)
    {
        printf("Das Ergebnis lautet:\n");
        for(i=0; i<glob_mat_size; ++i)
            printf("  %f  ",vectorX[i]);
        printf("\n");
    }
    /********** Finalize MPI **********/
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    return(0);
}
