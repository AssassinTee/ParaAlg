/**********************************************************************************************/
/* Kompilieren: mpicc -o fileio fileio.c                                                      */
/* Ausfuehren mit sbatch                                                                      */ */
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
    char text = '0'+world_rank;
    MPI_File_write_at(myfile, world_rank, &text, 1, MPI_CHAR, &status);

    MPI_Finalize();
}
