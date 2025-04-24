#include "stdio.h"
#include "mpi.h"

#include "amracut.h"
#include <stdlib.h>
#include <limits.h>
#include <string.h> 

int main(int argc, char *argv[])
{
  MPI_Init(&argc, &argv);
  MPI_Comm comm = MPI_COMM_WORLD;
  int procs_n, my_rank;
  MPI_Comm_size(comm, &procs_n);
  MPI_Comm_rank(comm, &my_rank);

  /**
   * example octree structure
   * -----------------------------------
   * |     | 1 | 2 |         |  8 |  9 | 
   * |  0  |---|---|    7    |---------|
   * |     | 3 | 4 |         | 10 | 11 |
   * -----------------------------------
   * |     |       |         |         |
   * |  5  |   6   |   12    |   13    |
   * |     |       |         |         |
   * -----------------------------------
   * |     |       | 18 | 19 |         |
   * | 14  |  15   |---------|   22    |
   * |     |       | 20 | 21 |         |
   * -----------------------------------
   * |     |       |         |         |
   * |  16 |   17  |   23    |   24    |
   * |     |       |         |         |
   * -----------------------------------
   * 
   * 
   * RANK0 has elements = [0, 6]
   * RANK1 has elements = [7, 13]
   * RANK2 has elements = [14, 17]
   * RANK3 has elements = [18, 24]
   * 
   * 
   * 
   * 
   */

  amracut_uint_t vtx_dist[5] = {0, 7, 14, 18, 25};

  oct_element *local_elements;

  switch (my_rank)
  {
  case 0:
    local_elements = malloc(sizeof(*local_elements) * 7);
    for (size_t i = 0; i < 7; i++)
    {
      for (size_t j = 0; j < 6; j++)
      {
        local_elements[i].e2e[j] = UINT_MAX;
      }
    }

    local_elements[0].eid = 0;
    local_elements[0].e2e[0] = 5;

    local_elements[1].eid = 1;
    local_elements[1].e2e[0] = 0;
    local_elements[1].e2e[1] = 3;
    local_elements[1].e2e[2] = 2;

    local_elements[2].eid = 2;
    local_elements[2].e2e[0] = 1;
    local_elements[2].e2e[1] = 4;
    local_elements[2].e2e[2] = 7;

    local_elements[3].eid = 3;
    local_elements[3].e2e[0] = 0;
    local_elements[3].e2e[1] = 6;
    local_elements[3].e2e[2] = 4;
    local_elements[3].e2e[3] = 1;

    local_elements[4].eid = 4;
    local_elements[4].e2e[0] = 3;
    local_elements[4].e2e[1] = 6;
    local_elements[4].e2e[2] = 7;
    local_elements[4].e2e[3] = 2;

    local_elements[5].eid = 5;
    local_elements[5].e2e[0] = 14;
    local_elements[5].e2e[1] = 6;
    local_elements[5].e2e[2] = 0;

    local_elements[6].eid = 6;
    local_elements[6].e2e[0] = 5;
    local_elements[6].e2e[1] = 15;
    local_elements[6].e2e[2] = 12;
    break;

  case 1:
    local_elements = malloc(sizeof(*local_elements) * 7);
    for (size_t i = 0; i < 7; i++)
    {
      for (size_t j = 0; j < 6; j++)
      {
        local_elements[i].e2e[j] = UINT_MAX;
      }
    }

    local_elements[0].eid = 7;
    local_elements[0].e2e[0] = 12;

    local_elements[1].eid = 8;
    local_elements[1].e2e[0] = 7; // Connected to larger element 7
    local_elements[1].e2e[1] = 10;
    local_elements[1].e2e[2] = 9;

    local_elements[2].eid = 9;
    local_elements[2].e2e[0] = 8;
    local_elements[2].e2e[1] = 11;

    local_elements[3].eid = 10;
    local_elements[3].e2e[0] = 8;
    local_elements[3].e2e[1] = 13;
    local_elements[3].e2e[2] = 11;
    local_elements[3].e2e[3] = 7;


    local_elements[4].eid = 11;
    local_elements[4].e2e[0] = 10;
    local_elements[4].e2e[1] = 13;
    local_elements[4].e2e[2] = 9;

    local_elements[5].eid = 12;
    local_elements[5].e2e[0] = 6;
    local_elements[5].e2e[1] = 13;
    local_elements[5].e2e[3] = 7;

    local_elements[6].eid = 13;
    local_elements[6].e2e[0] = 12;
    local_elements[6].e2e[1] = 22;


    break;

  case 2:
    local_elements = malloc(sizeof(*local_elements) * 4);
    for (size_t i = 0; i < 4; i++)
    {
      for (size_t j = 0; j < 6; j++)
      {
        local_elements[i].e2e[j] = UINT_MAX;
      }
    }

    local_elements[0].eid = 14;
    local_elements[0].e2e[0] = 16;
    local_elements[0].e2e[1] = 15;
    local_elements[0].e2e[2] = 5;

    local_elements[1].eid = 15;
    local_elements[1].e2e[0] = 14;
    local_elements[1].e2e[1] = 17;
    local_elements[1].e2e[2] = 6;

    local_elements[2].eid = 16;
    local_elements[2].e2e[0] = 14;
    local_elements[2].e2e[1] = 17;

    local_elements[3].eid = 17;
    local_elements[3].e2e[0] = 16;
    local_elements[3].e2e[1] = 15;
    local_elements[3].e2e[2] = 23;
    break;

  case 3:
    local_elements = malloc(sizeof(*local_elements) * 7);
    for (size_t i = 0; i < 7; i++)
    {
      for (size_t j = 0; j < 6; j++)
      {
        local_elements[i].e2e[j] = UINT_MAX;
      }
    }

    local_elements[0].eid = 18;
    local_elements[0].e2e[0] = 15;
    local_elements[0].e2e[1] = 20;
    local_elements[0].e2e[2] = 19;
    local_elements[0].e2e[3] = 12;

    local_elements[1].eid = 19;
    local_elements[1].e2e[0] = 18;
    local_elements[1].e2e[1] = 21;
    local_elements[1].e2e[2] = 22;
    local_elements[1].e2e[3] = 12;


    local_elements[2].eid = 20;
    local_elements[2].e2e[0] = 18;
    local_elements[2].e2e[1] = 21;
    local_elements[2].e2e[2] = 23;
    local_elements[2].e2e[3] = 15;


    local_elements[3].eid = 21;
    local_elements[3].e2e[0] = 20;
    local_elements[3].e2e[1] = 23;
    local_elements[3].e2e[2] = 19;
    local_elements[3].e2e[3] = 22;

    local_elements[4].eid = 22;
    local_elements[4].e2e[0] = 24;
    local_elements[4].e2e[1] = 13; 

    local_elements[5].eid = 23;
    local_elements[5].e2e[0] = 17; 
    local_elements[5].e2e[1] = 24;

    local_elements[6].eid = 24;
    local_elements[6].e2e[0] = 23;
    local_elements[6].e2e[1] = 22;
    break;

  default:
    break;
  }

  amracut_uint_t *parts = (amracut_uint_t *)malloc(10 * sizeof(amracut_uint_t));


  amracut_partgraph_octree(vtx_dist, local_elements, parts, &comm);
  char output[500] = {'\0'};
  char tempbuf[10];
  snprintf(tempbuf, 10, "%d", my_rank);
  strcat(output, tempbuf);
  strcat(output, " : [");

  for (amracut_uint_t i = 0; i < vtx_dist[my_rank + 1] - vtx_dist[my_rank]; i++)
  {
    snprintf(tempbuf, 10, "%d", parts[i]);
    strcat(output, tempbuf);
    strcat(output, ", ");
  }
  strcat(output, "]\n");

  printf("%s", output);
  


  MPI_Finalize();
  free(parts);

  return 0;
}
