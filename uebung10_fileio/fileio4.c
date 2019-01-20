/**********************************************************************************************/
/* Kompilieren: mpicc -o fileio4 fileio4.c                                                    */
/* Ausfuehren mit sbatch batch_fileio4                                                        */
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
#define len_block 50
#endif

#pragma message "content of glob_mat_size: " STR(len_vector)
#pragma message "content of nb: " STR(len_block)

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

    int loc_num_blocks = (len_vector/len_block)/world_size;
    int loc_vec_size = len_block*loc_num_blocks;
    double loc_vec[loc_vec_size];

    //Init vars
    MPI_File myfile;
    MPI_Status status;

    //�ffne u/o erzeuge outputdatei
    MPI_File_open(MPI_COMM_WORLD, "output3", MPI_MODE_RDONLY|MPI_MODE_CREATE, MPI_INFO_NULL, &myfile);

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


    //read//Remember to use OLD TYPE
    MPI_File_read(myfile, &loc_vec, loc_vec_size, MPI_DOUBLE, &status);

    //Init buffer
    double *glob_vec_pointer;
    if(!world_rank)
    {
        double glob_vec[len_vector];
        glob_vec_pointer = &glob_vec[0];//why use malloc
    }

    //Sammel Daten
    MPI_Gather(&loc_vec, loc_vec_size, MPI_DOUBLE, &glob_vec_pointer, len_vector, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    //Lokale test routine
    if(!world_rank)
    {
        for(int i = 0; i < len_vector; ++i)
        {
            //�berpr�fe ganz sauber mit schwellwert
            if(abs(glob_vec_pointer[i]-i-1) > 10e-10)
                printf("line %d is wrong: %f != %d", i, glob_vec_pointer[i], i+1);
        }
    }

    MPI_Finalize();
}

