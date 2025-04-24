#include "mpi.h"
#include <vector>
#include <algorithm>
#include <unordered_set>
#include <numeric>
#include <climits>

#include "octree.hpp"
#include "utils.hpp"

#define OCTREE_ELEMENT_FACES_N 6 
#define OCTREE_BOUNDARY UINT_MAX


// Overloading the << operator for oct_element
std::ostream& operator<<(std::ostream& os, const oct_element& obj){
    os << "(" << obj.eid << ", [" << obj.e2e[0];
    for (int i = 1; i < 6; i++)
    {
      os << "," << obj.e2e[i];
    }
    os << "] )";
    
    return os;
}

namespace amracut
{

  void dgraph_from_octree(const amracut_uint_t *vtx_dist__,
                          const oct_element *local_elements__, MPI_Comm *comm,
                          DGraph** dgraph)
  {
    int procs_n, my_rank;
    MPI_Comm_size(*comm, &procs_n);
    MPI_Comm_rank(*comm, &my_rank);

    const amracut_uint_t local_element_count = vtx_dist__[my_rank + 1] - vtx_dist__[my_rank];

    std::vector<oct_element> local_elements(local_elements__, local_elements__ + local_element_count);
    std::vector<amracut_uint_t> vtx_dist(vtx_dist__, vtx_dist__ + (procs_n + 1));


    /**
     * sort local elements by global element id
     */
    std::sort(local_elements.begin(), local_elements.end(),
              [](const oct_element &a, const oct_element &b)
              { 
                return a.eid < b.eid; 
              });

    std::vector<std::unordered_set<amracut_uint_t>> e2e_local(local_element_count);
    std::vector<std::pair<amracut_uint_t, amracut_uint_t>> e2e_boundary_pairs;  // a pair will be formed as (bdry element, inside element)

    // for (amracut_uint_t elem_i = 0; elem_i < local_element_count; elem_i++)
    // {
    //   e2e_local[elem_i] = {};
    // }

    for (const oct_element &e : local_elements)
    {
      amracut_uint_t local_id = e.eid - vtx_dist[my_rank];
      // e2e_local[local_id] = {};
      for (int f = 0; f < OCTREE_ELEMENT_FACES_N; f++)
      {
        amracut_uint_t neighbor_id = e.e2e[f];
        if (neighbor_id != OCTREE_BOUNDARY)
        {
          if (neighbor_id < vtx_dist__[my_rank] || vtx_dist__[my_rank + 1] <= neighbor_id)
          { 
            // boundary edge

            e2e_boundary_pairs.push_back({neighbor_id, e.eid});
            e2e_local[local_id].insert(neighbor_id);
          }else
          {
            // local edge, add edge to both elements

            e2e_local[local_id].insert(neighbor_id);
            amracut_uint_t neighbor_id_local = neighbor_id - vtx_dist[my_rank];
            e2e_local[neighbor_id_local].insert(e.eid);
          }
          
        }
      }
    }

    /**
     * sort bdry pairs by bdry element id
     */
    std::sort(e2e_boundary_pairs.begin(), e2e_boundary_pairs.end(),
              [](const auto &a, const auto &b)
              { 
                return a.first < b.first; 
              });

    std::vector<int> pair_send_counts(procs_n, 0);
    {
      int current_proc = 0;
      for (amracut_uint_t i = 0; i < e2e_boundary_pairs.size(); i++)
      {
        if (e2e_boundary_pairs[i].first < vtx_dist[current_proc + 1])
        {
          pair_send_counts[current_proc]++;
        }
        else
        {
          while (1)
          {
            current_proc++;
            if (e2e_boundary_pairs[i].first < vtx_dist[current_proc + 1])
            {
              pair_send_counts[current_proc]++;
              break;
            }
          }
        }
      }
    }
    
    std::vector<int> pair_send_displs(procs_n);
    std::exclusive_scan(pair_send_counts.begin(), pair_send_counts.end(), pair_send_displs.begin(), 0);

    std::vector<int> pair_recv_counts(procs_n, 0);
    MPI_Alltoall(pair_send_counts.data(), 1, MPI_INT, pair_recv_counts.data(), 1, MPI_INT, *comm);

    std::vector<int> pair_recv_displs(procs_n);
    std::exclusive_scan(pair_recv_counts.begin(), pair_recv_counts.end(), pair_recv_displs.begin(), 0);


    int pair_recv_total_count = pair_recv_counts[procs_n-1] + pair_recv_displs[procs_n-1];
    std::vector<std::pair<amracut_uint_t, amracut_uint_t>> e2e_boundary_pairs_recevied(pair_recv_total_count);

    MPI_Alltoallv(e2e_boundary_pairs.data(),
                  pair_send_counts.data(), pair_send_displs.data(), Mpi_pairtype<amracut_uint_t, amracut_uint_t>::value(), 
                  e2e_boundary_pairs_recevied.data(), pair_recv_counts.data(), pair_recv_displs.data(),
                  Mpi_pairtype<amracut_uint_t, amracut_uint_t>::value(), *comm);
    // print_log_mpi(my_rank, "sent: ", VectorToString(e2e_boundary_pairs));
    // print_log_mpi(my_rank, "received: ", VectorToString(e2e_boundary_pairs_recevied));
    for (const std::pair<amracut_uint_t, amracut_uint_t> &recev_bdry_pair : e2e_boundary_pairs_recevied)
    {
      amracut_uint_t local_id = recev_bdry_pair.first - vtx_dist[my_rank];
      e2e_local[local_id].insert(recev_bdry_pair.second);
    }



    /**
     * now we are creating distributed CSR encoding for the graph
     * 
     */
    std::vector<amracut_uint_t> xadj(local_element_count+1);
    xadj[0] = 0;

    for (amracut_uint_t elem_i = 0; elem_i < local_element_count; elem_i++)
    {
      xadj[elem_i + 1] = xadj[elem_i] + e2e_local[elem_i].size();
    }

    std::vector<amracut_uint_t> adjncy(xadj[local_element_count]);

    amracut_uint_t running_idx = 0;

    for (amracut_uint_t elem_i = 0; elem_i < local_element_count; elem_i++)
    {
      for (auto neighbor_id : e2e_local[elem_i])
      {
        adjncy[running_idx++] = neighbor_id;
      }
    }

    *dgraph = new DGraph(vtx_dist.data(), xadj.data(), adjncy.data(), NULL, NULL, AMRACUT_UNWEIGHTED, comm);
  }

}