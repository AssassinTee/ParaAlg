#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <mpi.h>
#include <math.h>
#define MAXNUM 10000

void merge(int *, int *, int, int, int);
void mergeSort(int *, int *, int, int);

int main(int argc, char** argv) {

	/********** Create and populate the array **********/
    int n;   /* Anzahl der zu sortierenden Zahlen */
    int ieing[MAXNUM];
    int i,j,k, c;

    n=13;
    for (i=0;i<n;++i) {
        ieing[i]=rand();
        printf("%d ",ieing[i]);
    }
    printf("\n ");

	/********** Initialize MPI **********/
	int world_rank;
	int world_size;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);

	MPI_Status status;
	MPI_Request req;

	/********** Divide the array in equal-sized chunks **********/
	int size = n/world_size;

	/********** Send each subarray to each process **********/
	int *sub_array = malloc(n * sizeof(int));
	//Todo: Change with Scatter V
	MPI_Scatter(&ieing[0], size, MPI_INT, sub_array, size, MPI_INT, 0, MPI_COMM_WORLD);
	printf("scattered; size: %d\n", size);

	/********** Perform the mergesort on each process **********/
	int *tmp_array = malloc(n * sizeof(int));
	mergeSort(sub_array, tmp_array, 0, (size - 1));
	printf("rank: %d; merge sorted sub array", world_rank);

	/********** Gather the sorted subarrays into one **********/
	/********** Fan In ***********/
	int fans = log(n);
	int selected = 2;
	for(int i = 0; i < fans; ++i)
	{
        if(world_rank%selected == 0)
        {
            //recieve
            int recv_from = world_rank+selected/2;

            if(recv_from < n)//scatter_v
            {
                //recv blocked
                printf("%d recieving from %d\n", world_rank, recv_from);
                MPI_Recv(tmp_array+size, size, MPI_INT, recv_from, 0, MPI_COMM_WORLD, &status);
                merge(sub_array, tmp_array, 0, size, 2*size);
                memcpy(tmp_array, sub_array, size);
                size*=2;
            }


        }
        else if(world_rank%selected/2 == 0)
        {
            int send_to = world_rank-world_rank%selected;
            printf("%d sending to %d\n", world_rank, send_to);
            MPI_Isend(tmp_array, size, MPI_INT, send_to, 0, MPI_COMM_WORLD, &req);
            break;
        }
        selected*=2;
	}

	/*int *sorted = NULL;
	if(world_rank == 0) {

		sorted = malloc(n * sizeof(int));

		}

	MPI_Gather(sub_array, size, MPI_INT, sorted, size, MPI_INT, 0, MPI_COMM_WORLD);
*/
	/********** Make the final mergeSort call **********/
	if(world_rank == 0) {
/*
		int *other_array = malloc(n * sizeof(int));
		mergeSort(sorted, other_array, 0, (n - 1));
*/
		/********** Display the sorted array **********/
		printf("This is the sorted array: ");
		for(c = 0; c < n; c++) {

			printf("%d ", tmp_array[c]);

			}

		printf("\n");
		printf("\n");

		/********** Clean up root **********/
		//free(tmp_array);
		//free(other_array);

		}

	/********** Clean up rest **********/
	//free(original_array);
	free(sub_array);
	free(tmp_array);

	/********** Finalize MPI **********/
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();

	}

/********** Merge Function **********/
void merge(int *a, int *b, int l, int m, int r) {

	int h, i, j, k;
	h = l;
	i = l;
	j = m + 1;

	while((h <= m) && (j <= r)) {

		if(a[h] <= a[j]) {

			b[i] = a[h];
			h++;

			}

		else {

			b[i] = a[j];
			j++;

			}

		i++;

		}

	if(m < h) {

		for(k = j; k <= r; k++) {

			b[i] = a[k];
			i++;

			}

		}

	else {

		for(k = h; k <= m; k++) {

			b[i] = a[k];
			i++;

			}

		}

	for(k = l; k <= r; k++) {

		a[k] = b[k];

		}

	}

/********** Recursive Merge Function **********/
void mergeSort(int *a, int *b, int l, int r) {

	int m;

	if(l < r) {

		m = (l + r)/2;

		mergeSort(a, b, l, m);
		mergeSort(a, b, (m + 1), r);
		merge(a, b, l, m, r);

		}

	}
