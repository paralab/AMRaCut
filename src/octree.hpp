#ifndef AMRACUT_OCTREE_H__
#define AMRACUT_OCTREE_H__

#include "mpi.h"
#include "amracut.h"
#include "dgraph.hpp"

namespace amracut
{
  /**
   * @brief Creates a graph from octree connectivity
   * 
   * @param vtx_dist__        process P_i stores the vertices from vtxdist[i] up to (but not including) vertex vtxdist[i + 1]
   * @param local_elements__  array of type `oct_element` containing element connectivity
   * @param comm              MPI communicator
   * @param dgraph            output dgraph
   * 
   * @note The caller is responsible to call `delete` on `dgraph`
   */
  void dgraph_from_octree(const amracut_uint_t *vtx_dist__,
                          const oct_element *local_elements__, MPI_Comm *comm,
                          DGraph** dgraph);

}
#endif