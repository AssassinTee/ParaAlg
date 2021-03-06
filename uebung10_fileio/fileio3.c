/**********************************************************************************************/
/* Kompilieren: mpicc -o fileio3 fileio3.c                                                    */
/* Ausfuehren mit sbatch batch_fileio3                                                        */
/* Mit dem Befehl 'od -e output3' l�sst sich der output �berpr�fen!                           */
/**********************************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>

#define STR_HELPER(x) #x
#define STR(x) STR_HELPER(x)

#ifndef len_vector
#define len_vector 500
#endif

#ifndef len_block
#define len_block 25
#endif

#pragma message "vectorl�nge: " STR(len_vector)
#pragma message "blockl�nge: " STR(len_block)

int main(int argc, char *argv[])
{
    int world_rank, world_size;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    //Pr�fe ob die bl�cke sich gleichm��ig blockzyklisch aufteilen lassen
    if((len_vector/len_block)%world_size != 0)
    {
        MPI_Abort(MPI_COMM_WORLD, 42);
        return 1;
    }

    //Berechne lokale gr��en
    int loc_num_blocks = (len_vector/len_block)/world_size;
    int loc_vec_size = len_block*loc_num_blocks;
    double loc_vec[loc_vec_size];

    //Iteriere �ber lokale bl�cke
    for(int loc_block = 0; loc_block < loc_num_blocks; ++loc_block){
        //Berechne globalen Blockindex
        int glob_block_index = loc_block*world_size+world_rank;
        for(int off = 0; off < len_block; ++off)//offset im lokalen Block
        {
            //Berechne lokalen vectorreihen index
            int loc_row_index = loc_block*len_block+off;
            //Berechne globalen vector reihen index
            int glob_row_index = glob_block_index*len_block+off;
            loc_vec[loc_row_index] = glob_row_index+1;
        }
    }

    //Init
    MPI_File myfile;
    MPI_Status status;

    //�ffne u/o erzeuge Outputdatei
    MPI_File_open(MPI_COMM_WORLD, "output3", MPI_MODE_WRONLY|MPI_MODE_CREATE, MPI_INFO_NULL, &myfile);

    //Initialisiere Werte f�r subarray
    int ndims = 1;//Vector ist 1D
    int sizes[1], subsizes[1], starts[1];
    sizes[0] = world_size   * len_block;
    subsizes[0] = 1         * len_block;
    starts[0] = world_rank  * len_block;//start ist am world_rank'ten block

    //Create new datatype
    MPI_Datatype mydatatype;
    MPI_Type_create_subarray(ndims, sizes, subsizes, starts, MPI_ORDER_C, MPI_DOUBLE, &mydatatype);
    MPI_Type_commit(&mydatatype);

    //Create Fileview//Remember to use OLD TYPE
    MPI_File_set_view(myfile, 0, MPI_DOUBLE, mydatatype, "native", MPI_INFO_NULL);


    //write//Remember to use OLD TYPE
    MPI_File_write(myfile, &loc_vec, loc_vec_size, MPI_DOUBLE, &status);

    MPI_Finalize();
}

