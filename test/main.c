#include "stdio.h"
#include "mpi.h"

#include "amracut.h"
#include <stdlib.h>


int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);
    MPI_Comm comm = MPI_COMM_WORLD;
    int procs_n, my_rank;
    MPI_Comm_size(comm, &procs_n);
    MPI_Comm_rank(comm, &my_rank);


    // printf("rank: %d, PID: %d\n", my_rank, getpid());

    amracut_ctrl ctrl;

    amracut_uint_t vtx_dist[5] = {0, 5, 10, 15};

    amracut_uint_t xadj_0[6] = {0, 2, 5, 8, 11, 13};
    amracut_uint_t adjncy_0[13] = {1, 5, 0, 2, 6, 1, 3, 7, 2, 4, 8, 3, 9};

    amracut_uint_t xadj_1[6] = {0, 3, 7, 11, 15, 18};
    amracut_uint_t adjncy_1[18] = {0, 6, 10, 1, 5, 7, 11, 2, 6, 8, 12, 3, 7, 9, 13, 4, 8, 14};

    amracut_uint_t xadj_2[6] = {0, 2, 5, 8, 11, 13};
    amracut_uint_t adjncy_2[13] = {5, 11, 6, 10, 12, 7, 11, 13, 8, 12, 14, 9, 13};

    amracut_uint_t *xadj;
    amracut_uint_t *adjncy;
    // amracut_uint_t *parts;

    switch (my_rank)
    {
    case 0:
        xadj = xadj_0;
        adjncy = adjncy_0;
        break;
    case 1:
        xadj = xadj_1;
        adjncy = adjncy_1;
        break;
    case 2:
        xadj = xadj_2;
        adjncy = adjncy_2;
        break;
    
    default:
        break;
    }

    amracut_setup(&ctrl, vtx_dist, xadj, adjncy, NULL, NULL, AMRACUT_UNWEIGHTED, &comm);

    amracut_uint_t *parts = (amracut_uint_t *)malloc(10*sizeof(amracut_uint_t));
    amracut_partgraph(&ctrl,parts, true, &comm, 1);

    // printf("%d\n", rest);
    amracut_destroy(&ctrl);

    MPI_Finalize();
    free(parts);

    return 0;
}
