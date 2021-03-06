/**********************************************************************************************/
/* Kompilieren: mpicc -o fileio4 fileio4.c -Dlen_block=25                                     */
/* Ausfuehren mit sbatch batch_fileio4                                                        */
/* Mit dem Befehl 'od -e output3' l�sst sich der Input �berpr�fen!                            */
/* Getestet mit:                                                                              */
/*      - Dlen_block=25, 4 procs -> 5 lokale Bl�cke                                           */
/*      - Dlen_block=125, 4 procs -> 1 lokaler Block                                          */
/*      - Dlen_block=10, 10 procs -> 5 lokale Bl�cke                                          */
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

#pragma message "vectorl�nge: " STR(len_vector)
#pragma message "blockl�nge: " STR(len_block)

static void handle_error(int errcode, char *str)
{
    char msg[MPI_MAX_ERROR_STRING];
    int resultlen;
    MPI_Error_string(errcode, msg, &resultlen);
    printf("%s: %s\n", str, msg);
    MPI_Abort(MPI_COMM_WORLD, 1);
}

#define MPI_CHECK(fn) { int errcode; errcode = (fn);\
     if (errcode != MPI_SUCCESS) handle_error  (errcode, #fn ); }

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

    //�ffne output datei
    MPI_CHECK(MPI_File_open(MPI_COMM_WORLD, "output3", MPI_MODE_RDONLY, MPI_INFO_NULL, &myfile));

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
    MPI_CHECK(MPI_File_read(myfile, &loc_vec, loc_vec_size, MPI_DOUBLE, &status));

    //Debug
    /*if(!world_rank)
    {
        for(int i = 0; i < loc_vec_size; ++i)
            printf("wr %d, line %d: %f\n",world_rank, i , loc_vec[i]);
    }*/

    //Sammel Daten
    //Init buffer
    double *glob_vec;
    if(!world_rank)
        glob_vec = malloc(sizeof(double)*len_vector);//why use malloc

    int* displs = malloc(sizeof(int)*world_size);
    int* rcounts = malloc(sizeof(int)*world_size);
    for(int i=0; i < world_size; ++i)
    {
        displs[i] = i*len_block;
        rcounts[i] = len_block;
    }

    //Sammel world_size bl�cke pro gather, loc_num_blocks mal, geht bestimmt besser mit eigenen datentypen
    //MPI_Gather(loc_vec, 1, mydatatype, glob_vec, len_vector, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    for(int i = 0; i < loc_num_blocks; ++i)
    {
        MPI_Gatherv(loc_vec+i*len_block, len_block, MPI_DOUBLE, glob_vec+i*world_size*len_block, rcounts, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }

    free(displs);
    free(rcounts);

    //Lokale test routine
    if(!world_rank)
    {
        int err = 0;
        int correct = 0;
        for(int i = 0; i < len_vector; ++i)
        {
            //�berpr�fe ganz sauber mit schwellwert
            if(abs(glob_vec[i]-i-1) > 10e-10) {
                printf("line %d is wrong: %f != %d\n", i, glob_vec[i], i+1);
                err++;
            }
            else
                correct++;
        }
        if(err)
            printf("Vektor is not okay! %d correct, %d false\n", correct, err);
        else
            printf("Vector is okay! %d correct\n", correct);
        free(glob_vec);
    }

    MPI_Finalize();
}

