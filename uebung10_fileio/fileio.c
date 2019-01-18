/**********************************************************************************************/
/* Kompilieren: mpicc -o fileio fileio.c                                                      */
/* Ausfuehren mit sbatch                                                                      */
/**********************************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>

int main(int argc, char *argv[])
{
    int world_rank, world_size;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    MPI_File myfile;
    MPI_Status status;
    MPI_File_open(MPI_COMM_WORLD, "output", MPI_MODE_WRONLY|MPI_MODE_CREATE, MPI_INFO_NULL, &myfile);

    //char buf[world_size];
    int ndims = 1;
    int sizes[1], subsizes[1], starts[1];
    sizes[0] = world_size;
    subsizes[0] = 1;
    starts[0] = world_rank;
    MPI_Datatype mydatatype;
    MPI_Type_create_subarray(ndims, sizes, subsizes, starts, MPI_ORDER_C, MPI_CHAR, &mydatatype);
    MPI_Type_commit(&mydatatype);
    MPI_File_set_view(myfile, 0, MPI_CHAR, mydatatype, "native", MPI_INFO_NULL);

    //create text
    int loop = 10;
    char text[loop];
    for(int i = 0; i < loop; ++i)
        text[i] = '0'+world_rank;

    //write
    MPI_File_write(myfile, &text, loop, MPI_CHAR, &status);

    MPI_Finalize();
}
