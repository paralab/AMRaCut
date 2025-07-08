#include "dgraph.hpp"
#include <algorithm>
#include <numeric>
#include <string>
#include <cstring>
#include <chrono>
#include <unordered_map>
#include <unordered_set>


#define PRINT_RANK "[",this->my_rank,"] "

#if 0
#define TIMED_COMM(...) \
{ \
  auto com_start_ = std::chrono::high_resolution_clock::now(); \
  __VA_ARGS__ \
  auto com_end_ = std::chrono::high_resolution_clock::now(); \
  this->com_duration += std::chrono::duration_cast<std::chrono::microseconds>(com_end_ - com_start_); \
}
#endif

namespace amracut
{
  // Overloading the << operator for BFSValue
  std::ostream &operator<<(std::ostream &os, const amracut::BFSValue &obj)
  {
    os << "(d:" << obj.distance << ",l:" << obj.label << ")";

    return os;
  }

  // Overloading the << operator for DiffusionValue
  std::ostream &operator<<(std::ostream &os, const amracut::DiffusionValue &obj)
  {
    os << "(v:" << obj.value << ",l:" << obj.label << ")";

    return os;
  }

  // Overloading the << operator for GhostBFSValue
  std::ostream &operator<<(std::ostream &os, const amracut::GhostBFSValue &obj)
  {
    os << "(d:" << obj.distance << ",l:" << obj.label << ",o:" << obj.offset << ")";
    // os << obj.global_idx;

    return os;
  }
}

amracut::DGraph::DGraph(const amracut_uint_t *vtx_dist_, const amracut_uint_t *xadj_, const amracut_uint_t *adjncy_, 
                         const amracut_uint_t *vwgt_, const amracut_uint_t *adjwgt_, const amracut_uint_t wgtflag_, 
                         MPI_Comm *comm_)
{
  this->comm = comm_;

  MPI_Comm_size(*comm, &(this->procs_n));
  MPI_Comm_rank(*comm, &(this->my_rank));

  this->wgtflag = wgtflag_;

  this->global_n  = 0;
  this->local_n   = 0;
  this->ghost_n   = 0;
  this->send_n    = 0;

  this->ghost_procs_n   = 0;
  this->send_procs_n    = 0;

  this->has_ghost = false;

  this->local_n = vtx_dist_[my_rank + 1] - vtx_dist_[my_rank];
  this->global_n = vtx_dist_[this->procs_n];


  // copy adjacency structure into local buffers
  this->local_xadj.resize(this->local_n + 1);
  this->local_adjncy.resize(xadj_[this->local_n]);
  std::memcpy(this->local_xadj.data(), xadj_, sizeof(amracut_uint_t) * this->local_xadj.size());
  std::memcpy(this->local_adjncy.data(), adjncy_, sizeof(amracut_uint_t) * this->local_adjncy.size());     // local_adjncy will be modified later with local indexing

  this->vtx_dist.resize(this->procs_n + 1);
  std::memcpy(this->vtx_dist.data(), vtx_dist_, sizeof(amracut_uint_t) * this->vtx_dist.size()); 

  if(this->wgtflag == AMRACUT_VTX_WEIGHTED || this->wgtflag == AMRACUT_VTX_EDGE_WEIGHTED)
  {
    this->vwgt.resize(this->local_n);
    std::memcpy(this->vwgt.data(), vwgt_, sizeof(amracut_uint_t) * this->vwgt.size());
  }

  if(this->wgtflag == AMRACUT_EDGE_WEIGHTED || this->wgtflag == AMRACUT_VTX_EDGE_WEIGHTED)
  {
    this->adjwgt.resize(this->local_adjncy.size());
    std::memcpy(this->adjwgt.data(), adjwgt_, sizeof(amracut_uint_t) * this->adjwgt.size());
  }

  amracut_uint_t ghost_edges_n = 0;

  // TODO: can be a parallel reduction
  for (size_t i = 0; i < this->local_adjncy.size(); i++)
  {
    if (this->local_adjncy[i] < vtx_dist[my_rank] || this->local_adjncy[i] >= vtx_dist[my_rank + 1])  // check if a ghost edge
    {
      ghost_edges_n++;
    }
  }
  if (ghost_edges_n == 0)
  {
    amracut::print_log("no ghost");
    return ;
  }

  this->has_ghost = true;

  std::vector<GhostEdge> ghost_edges(ghost_edges_n);
  std::vector<int> tmp_ghost_counts_per_proc(this->procs_n, 0);
  std::vector<int> tmp_ghost_edge_counts_per_proc(this->procs_n, 0);
  std::vector<int> tmp_ghost_edge_counts_per_proc_scanned(this->procs_n, 0);
  std::vector<int> tmp_send_counts_per_proc(this->procs_n, 0);


  amracut_uint_t idx = 0;

  for (amracut_uint_t local_v = 0; local_v < this->local_n; local_v++)    // for each local vertex v
  {
    for (amracut_uint_t adjncy_i = this->local_xadj[local_v]; adjncy_i < this->local_xadj[local_v+1]; adjncy_i++)   // for each edge in vertex v
    {
      amracut_uint_t u = this->local_adjncy[adjncy_i];
      if (u >= vtx_dist[my_rank] && u < vtx_dist[my_rank + 1])    // local edge
      {
        this->local_adjncy[adjncy_i] -= vtx_dist[my_rank];    // change this->local_adjncy into local indexing 
        continue;
      }

      // Else, we have a ghost edge now, add it to buffer
      ghost_edges[idx].local_vertex   = local_v;
      ghost_edges[idx].ghost_vertex   = u;
      ghost_edges[idx].adjncy_idx     = adjncy_i;
      if (this->wgtflag == AMRACUT_EDGE_WEIGHTED || this->wgtflag == AMRACUT_VTX_EDGE_WEIGHTED)
      {
        ghost_edges[idx].wgt          = this->adjwgt[adjncy_i];
      }
      
      idx++;
      
    }
    
  }

    
  // sort by the global index of the ghost element
  // TODO: use a parallel sort?
  std::sort(ghost_edges.begin(), ghost_edges.end(),
            [](const GhostEdge &a, const GhostEdge &b)
            {
              return a.ghost_vertex < b.ghost_vertex;
            });



  // now we assign local indices for ghost vertices.
  // we do a contiguous assignment s.t. local indices of ghost vertices are contiguous for each ghost proc and can be directly used with MPIAllToAllV




  amracut_uint_t ghost_local_idx = this->local_n;     // local indices for ghosts start from the end of local vertex indices.

  ghost_edges[0].ghost_local_idx = ghost_local_idx;
  this->local_adjncy[ghost_edges[0].adjncy_idx] = ghost_local_idx;

  int ghost_proc_i = 0;
  while (!(ghost_edges[0].ghost_vertex >= this->vtx_dist[ghost_proc_i] && ghost_edges[0].ghost_vertex < this->vtx_dist[ghost_proc_i+1]))  // identifying which ghost proc this ghost vertex belongs to
  {
    ghost_proc_i++;
  }
  tmp_ghost_counts_per_proc[ghost_proc_i]++;
  tmp_ghost_edge_counts_per_proc[ghost_proc_i]++;
  ghost_edges[0].ghost_proc = ghost_proc_i;
  

  for (amracut_uint_t i = 1; i < ghost_edges_n; i++)
  {

    if (ghost_edges[i-1].ghost_vertex == ghost_edges[i].ghost_vertex)   // ghost duplicate detected
    {
      ghost_edges[i].ghost_local_idx = ghost_local_idx;
      this->local_adjncy[ghost_edges[i].adjncy_idx] = ghost_local_idx;
      ghost_edges[i].ghost_proc = ghost_proc_i;
      tmp_ghost_edge_counts_per_proc[ghost_proc_i]++;
      continue;
    }
    
    ghost_local_idx++;    // not a duplicate
    
    ghost_edges[i].ghost_local_idx = ghost_local_idx;
    this->local_adjncy[ghost_edges[i].adjncy_idx] = ghost_local_idx;

    while (!(ghost_edges[i].ghost_vertex >= this->vtx_dist[ghost_proc_i] && ghost_edges[i].ghost_vertex < this->vtx_dist[ghost_proc_i+1]))    // update proc ownership of ghost
    {
      ghost_proc_i++;
    }
    tmp_ghost_counts_per_proc[ghost_proc_i]++;
    tmp_ghost_edge_counts_per_proc[ghost_proc_i]++;
    ghost_edges[i].ghost_proc = ghost_proc_i;
       
  }

  this->ghost_n = ghost_local_idx + 1 - this->local_n;

  // now we add back-edges from ghosts to local vertices, into the adjacency structure
  // this is useful when we have to continue BFS from ghosts to inside

  this->local_xadj.resize(this->local_n + this->ghost_n + 1);
  this->local_adjncy.resize(this->local_xadj[local_n] + ghost_edges_n);
  if (this->wgtflag == AMRACUT_EDGE_WEIGHTED || this->wgtflag == AMRACUT_VTX_EDGE_WEIGHTED)
  {
    this->adjwgt.resize(this->local_xadj[local_n] + ghost_edges_n);
  }
  

  std::vector<amracut_uint_t> ghost_degrees(this->ghost_n, 0);

 
  for (amracut_uint_t i = 0; i < ghost_edges_n; i++)
  {
    ghost_degrees[ghost_edges[i].ghost_local_idx - this->local_n]++;
  }

  {
    amracut_uint_t ghost_edge_idx = 0;
    for (amracut_uint_t g_i = 0; g_i < this->ghost_n; g_i++)
    {
      amracut_uint_t xadj_idx = this->local_n + g_i + 1;
      this->local_xadj[xadj_idx] = this->local_xadj[xadj_idx - 1] + ghost_degrees[g_i];
      for (amracut_uint_t adjncy_i = this->local_xadj[xadj_idx - 1]; adjncy_i < this->local_xadj[xadj_idx]; adjncy_i++)
      {
        this->local_adjncy[adjncy_i] = ghost_edges[ghost_edge_idx].local_vertex;
        if (this->wgtflag == AMRACUT_EDGE_WEIGHTED || this->wgtflag == AMRACUT_VTX_EDGE_WEIGHTED)
        {
          this->adjwgt[adjncy_i] = ghost_edges[ghost_edge_idx].wgt;
        }
        ghost_edge_idx++;
      }
      
    }
  }

  // preprocess adjwgt into normalized values in range (epsilon, 1]. this is useful when performing diffusion
  // use use a small positive constant `epsilon` to prevent edge weights normalizing to 0
  if (this->wgtflag == AMRACUT_EDGE_WEIGHTED || this->wgtflag == AMRACUT_VTX_EDGE_WEIGHTED)
  {
    auto wgt_minmax  = std::minmax_element(this->adjwgt.begin(), this->adjwgt.end());
    amracut_uint_t local_edge_wgt_min = *wgt_minmax.first;
    amracut_uint_t local_edge_wgt_max = *wgt_minmax.second;

    amracut_uint_t global_edge_wgt_min = std::numeric_limits<amracut_uint_t>::max();
    amracut_uint_t global_edge_wgt_max = 0;

    MPI_Allreduce(&local_edge_wgt_min, &global_edge_wgt_min, 1, Mpi_datatype<amracut_uint_t>::value(), MPI_MIN, *(this->comm));
    MPI_Allreduce(&local_edge_wgt_max, &global_edge_wgt_max, 1, Mpi_datatype<amracut_uint_t>::value(), MPI_MAX, *(this->comm));

    this->adjwgt_normalized.resize(this->adjwgt.size());

    if (global_edge_wgt_min == global_edge_wgt_max)     // if all edge weights are equal
    {
      std::fill(this->adjwgt_normalized.begin(), this->adjwgt_normalized.end(), static_cast<amracut_real_t>(1));
    } else
    {
      std::transform(this->adjwgt.begin(), this->adjwgt.end(), this->adjwgt_normalized.begin(),
                    [global_edge_wgt_min, global_edge_wgt_max](auto w)
                    { 
                      amracut_real_t epsilon = 0.8;
                                            
                      return epsilon + 
                             (1-epsilon) * static_cast<amracut_real_t>(w - global_edge_wgt_min) / (global_edge_wgt_max - global_edge_wgt_min); 
                    });
    }
    
  }

  // now populate ghost procs information. This is useful when doing an optimized sparse mpi all to all with neighbor procs only.
  for (int p_i = 0; p_i < procs_n; p_i++)
  {
    this->ghost_procs_n += (int)(tmp_ghost_counts_per_proc[p_i]>0);
    
  }
  this->ghost_procs.resize(this->ghost_procs_n);
  this->ghost_counts_per_proc.resize(this->ghost_procs_n);
  this->ghost_counts_per_proc_scanned.resize(this->ghost_procs_n);
  this->updated_ghost_counts_per_proc.resize(this->ghost_procs_n);
  this->updated_ghost_counts_per_proc_scanned.resize(this->ghost_procs_n);
  this->updated_ghost_buffer.resize(this->ghost_n);


  std::fill(this->ghost_counts_per_proc.begin(), this->ghost_counts_per_proc.end(),0);

  ghost_proc_i = 0;
  for (int p_i = 0; p_i < procs_n; p_i++)
  {
    if (tmp_ghost_counts_per_proc[p_i] > 0)
    {
      this->ghost_procs[ghost_proc_i] = p_i;
      this->ghost_counts_per_proc[ghost_proc_i] = tmp_ghost_counts_per_proc[p_i];
      ghost_proc_i++;

    }
    
  }
  std::exclusive_scan(this->ghost_counts_per_proc.begin(), this->ghost_counts_per_proc.end(),this->ghost_counts_per_proc_scanned.begin(),0);


  // now we build the sending scatter map

  // TODO: parallel scan
  std::exclusive_scan(tmp_ghost_edge_counts_per_proc.begin(), tmp_ghost_edge_counts_per_proc.end(),tmp_ghost_edge_counts_per_proc_scanned.begin(),0);

  for (int proc_i = 0; proc_i < procs_n; proc_i++)
  {
    int start = tmp_ghost_edge_counts_per_proc_scanned[proc_i];
    int count = tmp_ghost_edge_counts_per_proc[proc_i];
    if (count > 0)
    {
      std::sort(ghost_edges.begin() + start, ghost_edges.begin() + start + count,
            [](const GhostEdge &a, const GhostEdge &b)
            {
              return a.local_vertex < b.local_vertex;
            });
      tmp_send_counts_per_proc[proc_i]++;
      this->send_n++;
      for (int i = start+1; i < start+count; i++)
      {
        if (ghost_edges[i-1].local_vertex != ghost_edges[i].local_vertex)   // filtering out duplicate send vertices
        {
          tmp_send_counts_per_proc[proc_i]++;
          this->send_n++;
        }
      }
    }  
  }


  this->send_scatter_map.resize(this->send_n);
  idx = -1;

  for (int proc_i = 0; proc_i < procs_n; proc_i++)
  {
    int start = tmp_ghost_edge_counts_per_proc_scanned[proc_i];
    int count = tmp_ghost_edge_counts_per_proc[proc_i];
    if (count > 0)
    {
      idx++;
      this->send_scatter_map[idx] = ghost_edges[start].local_vertex;
      for (int i = start+1; i < start+count; i++)
      {
        if (ghost_edges[i-1].local_vertex != ghost_edges[i].local_vertex)   // filtering out duplicate send vertices
        {
          idx++;
          this->send_scatter_map[idx] = ghost_edges[i].local_vertex;
        }
      }
    }  
  }

  // now populate send procs information. This is useful when doing an optimized sparse mpi all to all with neighbor procs only.
  for (int p_i = 0; p_i < procs_n; p_i++)
  {
    this->send_procs_n += (int)(tmp_send_counts_per_proc[p_i]>0);
    
  }
  this->send_procs.resize(this->send_procs_n);
  this->send_counts_per_proc.resize(this->send_procs_n);
  this->send_counts_per_proc_scanned.resize(this->send_procs_n);
  this->updated_send_counts_per_proc.resize(this->send_procs_n);
  this->updated_send_counts_per_proc_scanned.resize(this->send_procs_n);
  this->updated_count_mpi_requests.resize(this->ghost_procs_n + this->send_procs_n);
  this->updated_count_mpi_statuses.resize(this->ghost_procs_n + this->send_procs_n);
  this->updated_send_buffer.resize(this->send_n);

  std::fill(this->send_counts_per_proc.begin(), this->send_counts_per_proc.end(),0);

  int send_proc_i = 0;
  for (int p_i = 0; p_i < procs_n; p_i++)
  {
    if (tmp_send_counts_per_proc[p_i] > 0)
    {
      this->send_procs[send_proc_i] = p_i;
      this->send_counts_per_proc[send_proc_i] = tmp_send_counts_per_proc[p_i];
      send_proc_i++;

    }
    
  }
  std::exclusive_scan(this->send_counts_per_proc.begin(), this->send_counts_per_proc.end(),this->send_counts_per_proc_scanned.begin(),0);

 
}


std::string amracut::DGraph::PrintLocal(){
    std::ostringstream output;
    for (size_t vertex_i = 0; vertex_i < this->local_n + this->ghost_n; vertex_i++)
    {
        output << vertex_i << "\t->";
        for (size_t neigh_i = this->local_xadj[vertex_i]; neigh_i < this->local_xadj[vertex_i+1]; neigh_i++)
        {
          if (this->wgtflag == AMRACUT_EDGE_WEIGHTED || this->wgtflag == AMRACUT_VTX_EDGE_WEIGHTED)
          {
            output << this->local_adjncy[neigh_i] << "(" << this->adjwgt[neigh_i] << "),";
          }
          else
          {
            output << this->local_adjncy[neigh_i] << ",";
          }
          
            
        }
        output << "\n";
        
    }

    return output.str();
    
}

amracut_uint_t amracut::DGraph::PartGraph(amracut_uint_t* partition_labels_out, bool use_diffusion, int verbose_)
{
  auto total_time_start = std::chrono::high_resolution_clock::now();
  this->com_duration    = std::chrono::microseconds(0);
  this->recv_size = 0;
  this->send_size = 0;

  this->verbose = verbose_;


  int round_counter = 0;
  int stop_guess = std::min(this->procs_n, 12);
  int guess_counter = 0;
  bool is_not_stable_global = true;      // global BFS stability


  // initialize temporary buffers for BFS
  this->frontier_buffer.resize(this->local_n + this->ghost_n);
  this->is_candidate.resize(this->local_n);
  this->candidates.resize(this->local_n);


  BFSValue bfs_init_value = {.label = AMRACUT_BFS_NO_LABEL, .distance =  AMRACUT_BFS_INFINITY};
  std::vector<BFSValue> bfs_vector(this->local_n + this->ghost_n, bfs_init_value);

  DiffusionValue diffusion_init_value = {.label = 0, .value = 1};
  std::vector<DiffusionValue> diffusion_vector(bfs_vector.size(), diffusion_init_value);


  this->RunFirstBFSIteration(bfs_vector);


  if (this->procs_n > 1)
  {
    std::vector<BFSValue> send_buffer(this->send_n);
    std::vector<BFSValue> send_buffer_prev(this->send_n);
    std::vector<BFSValue> recv_buffer(this->ghost_n);

    std::vector<graph_indexing_t> updated_ghosts(this->ghost_n);
    amracut_uint_t updated_ghost_n = 0;


    while (is_not_stable_global)
    {

      round_counter++;
      guess_counter++;
      is_not_stable_global = false;
      bool is_not_stable_local = true;



      if (round_counter > 1)
      {
        this->StartReceivingUpdatedOnlyGhostCounts();
        is_not_stable_local = this->RunBFSFromGhostUpdates(bfs_vector,updated_ghosts, updated_ghost_n);
      }

      

      for (size_t send_i = 0; send_i < this->send_n; send_i++)
      {
          send_buffer[send_i] = bfs_vector[this->send_scatter_map[send_i]];
      }

      if (round_counter > 1) // after the first round, exchanging updated only ghosts is better for communication
      {
        this->ExchangeUpdatedOnlyGhost(send_buffer, send_buffer_prev, recv_buffer);
      }
      else
      {
        this->ExchangeAllGhosts(send_buffer.data(), recv_buffer.data(), BFS_GHOST_EXCHANGE_TAG);
      }

      std::copy(send_buffer.begin(), send_buffer.end(), send_buffer_prev.begin());      // for the next iteration

      updated_ghost_n = 0;

      // Now update bfs_vector using received ghost values
      bool ghost_any_is_not_stable = false;

      for (size_t recv_i = 0; recv_i < this->ghost_n; recv_i++)
      {
        auto offset = this->local_n; // ghost elements are in the last section of the vector, in sorted order
        if (bfs_vector[offset + recv_i].distance != recv_buffer[recv_i].distance ||
            bfs_vector[offset + recv_i].label != recv_buffer[recv_i].label)
        {
          bfs_vector[offset + recv_i].distance = recv_buffer[recv_i].distance;
          bfs_vector[offset + recv_i].label = recv_buffer[recv_i].label;

          ghost_any_is_not_stable = true;
          updated_ghosts[updated_ghost_n++] = recv_i;
        }
      }
      is_not_stable_local = is_not_stable_local || ghost_any_is_not_stable;

      if (guess_counter >= stop_guess)
      {
        MPI_Allreduce(&is_not_stable_local, &is_not_stable_global, 1, MPI_CXX_BOOL, MPI_LOR, *(this->comm));
      }
      else
      {
        is_not_stable_global = true;
      }

    }   // end while (is_not_stable_global)


    // if(!this->my_rank) print_log("n:",this->global_n, "\tp:", this->procs_n, "\tBFS_iter_count:", round_counter); 
    
  }   // end if (this->procs_n > 1)


  for (size_t i = 0; i < bfs_vector.size(); i++)
  {
    diffusion_vector[i].label = bfs_vector[i].label;
    diffusion_vector[i].value = 1;
  }

  // now we run refinement using diffusion
  if(use_diffusion && (this->procs_n > 1))
  {
    this->RefineByDiffusion(diffusion_vector);
  }  

  for (amracut_uint_t i = 0; i < this->local_n; i++)      // copying labels to output
  {
    if (diffusion_vector[i].label == AMRACUT_BFS_NO_LABEL)
    {
      throw std::runtime_error("amracut ERROR: vertex "+ std::to_string(i+vtx_dist[this->my_rank]) + " was not assigned a partition label\n");
    }
    partition_labels_out[i] = diffusion_vector[i].label;

  }
  auto total_time_end = std::chrono::high_resolution_clock::now();
  int64_t total_time = std::chrono::duration_cast<std::chrono::microseconds>(total_time_end - total_time_start).count();
  // int64_t com_time = this->com_duration.count();

  if (this->verbose > 0)
  {
    if (!my_rank)
    {
      print_log("\n========================================");
    }
    
    int64_t total_time_min, total_time_max, total_time_avg;

    total_time_min = total_time_max = total_time_avg = 0;


    MPI_Reduce(&total_time, &total_time_min, 1, MPI_INT64_T, MPI_MIN, 0, *this->comm);
    MPI_Reduce(&total_time, &total_time_max, 1, MPI_INT64_T, MPI_MAX, 0, *this->comm);
    MPI_Reduce(&total_time, &total_time_avg, 1, MPI_INT64_T, MPI_SUM, 0, *this->comm);
    total_time_avg = total_time_avg / this->procs_n;

    if (!my_rank)
    {
      print_log("label propagation sync rounds : ", round_counter, "\n");
      print_log(amracut::FormatStatsV1(total_time_min, total_time_max, total_time_avg));
      print_log("========================================\n");
    }
  }
 

  bfs_vector.clear();

  this->frontier_buffer.clear();
  this->is_candidate.clear();
  this->candidates.clear();



  
  return 0;
}



void amracut::DGraph::RunFirstBFSIteration(std::vector<BFSValue> & bfs_vector)
{
  const bfs_label_t init_label = this->my_rank;
  const graph_indexing_t seed = this->local_n / 2;
  bfs_vector[seed].distance = 0;
  bfs_vector[seed].label = init_label;


  graph_indexing_t curr_frontier_start = 0;
  graph_indexing_t curr_frontier_size = 1;

  this->frontier_buffer[curr_frontier_start] = seed; // first frontier has only the seed

  bfs_distance_t curr_distance = 1;



  while (curr_frontier_size > 0)
  {

    graph_indexing_t next_frontier_size = 0;
    graph_indexing_t next_frontier_slot = curr_frontier_start + curr_frontier_size;
    for (graph_indexing_t frontier_i = curr_frontier_start;
         frontier_i < curr_frontier_start + curr_frontier_size;
         frontier_i++)
    {
      graph_indexing_t frontier_vertex = this->frontier_buffer[frontier_i];

      // looping neighbors
      for (graph_indexing_t neighbor_i = this->local_xadj[frontier_vertex];
           neighbor_i < this->local_xadj[frontier_vertex + 1]; neighbor_i++)
      {
        graph_indexing_t neighbor = this->local_adjncy[neighbor_i];

        if (neighbor >= this->local_n)
          continue; // we dont update the ghost vertices

        if (bfs_vector[neighbor].label == AMRACUT_BFS_NO_LABEL)
        {
          // neighbor is not visited, visit it now
          bfs_vector[neighbor].label = init_label;
          bfs_vector[neighbor].distance = curr_distance;

          // add to next frontier
          this->frontier_buffer[next_frontier_slot++] = neighbor;
          next_frontier_size++;
        }
      }
    }
    curr_distance++;
    curr_frontier_start += curr_frontier_size;
    curr_frontier_size = next_frontier_size;
  }
}

bool amracut::DGraph::RunBFSFromGhostUpdates(std::vector<BFSValue> &bfs_vector,
                                              std::vector<graph_indexing_t> &updated_ghosts, amracut_uint_t updated_ghosts_n)
{

  bool changed = false; // to detect if the BFS continued
  amracut_uint_t frontier_size = updated_ghosts_n;

  for (amracut_uint_t i = 0; i < updated_ghosts_n; i++)
  {
    this->frontier_buffer[i] = (updated_ghosts[i] + this->local_n); // offest for correct local indexing of ghosts
  }

  while (frontier_size > 0)
  {
    amracut_uint_t candidates_n = 0;
    std::fill(this->is_candidate.begin(), this->is_candidate.end(), 0);
    for (amracut_uint_t i = 0; i < frontier_size; i++)
    {
      graph_indexing_t frontier_vertex = frontier_buffer[i];

      for (graph_indexing_t neighbor_i = this->local_xadj[frontier_vertex];
           neighbor_i < this->local_xadj[frontier_vertex + 1]; neighbor_i++) // looping neighbors, add to next_updated
      {
        auto neighbor = this->local_adjncy[neighbor_i];
        if (neighbor >= this->local_n)
        {
          continue; // we dont update the ghost vertices
        }
        if (this->is_candidate[neighbor] == 0) // avoiding uplicate candidates
        {
          this->candidates[candidates_n++] = {.local_idx = neighbor, .bfs_value = bfs_vector[neighbor]};
          this->is_candidate[neighbor] = 1;
        }
      }
    }


    frontier_size = 0;

    // now considering only the neighborhoods of candidates
    for (amracut_uint_t c_i = 0; c_i < candidates_n; c_i++)
    {
      BFSCandidate &c = this->candidates[c_i];
      bool value_changed = false;

      for (graph_indexing_t neighbor_i = this->local_xadj[c.local_idx];
           neighbor_i < this->local_xadj[c.local_idx + 1]; neighbor_i++)
      {
        auto neighbor = this->local_adjncy[neighbor_i];

        if (bfs_vector[neighbor].label != AMRACUT_BFS_NO_LABEL &&
            c.bfs_value.distance > (bfs_vector[neighbor].distance + 1))
        {
          c.bfs_value.label = bfs_vector[neighbor].label;
          c.bfs_value.distance = bfs_vector[neighbor].distance + 1;

          if (!value_changed) // to prevent duplicates in next frontier
          {
            this->frontier_buffer[frontier_size++] = c.local_idx;
          }
          value_changed = true;
          changed = true;
        }
      }
    }

    // writing updated values to BFS vector
    for (amracut_uint_t c_i = 0; c_i < candidates_n; c_i++)
    {
      const BFSCandidate &c = this->candidates[c_i];
      bfs_vector[c.local_idx].label = c.bfs_value.label;
      bfs_vector[c.local_idx].distance = c.bfs_value.distance;
    }
  }

  return changed;
}


void amracut::DGraph::RefineByDiffusion(std::vector<DiffusionValue> &diffusion_vector)
{
  const int max_diffusion_iterations = 5;     // TODO: make this configurable

  this->diffusion_send_buffer.resize(this->send_n);
  this->diffusion_recv_buffer.resize(this->ghost_n);


  std::vector<int> level_1_2_neigh_procs;

  //   this->DiscoverLevel1and2NeighProcs(level_1_2_neigh_procs);



  std::vector<amracut_uint_t> partition_sizes(this->procs_n, 0);


  this->SyncGlobalPartSizes(diffusion_vector, partition_sizes);


  amracut_uint_t total_global_size = std::accumulate(partition_sizes.begin(), partition_sizes.end(), 0);

  std::unordered_map<amracut_uint_t, amracut_sint_t> local_parts_size_deltas;
  std::vector<amracut_real_t> diffusion_rates(this->procs_n);

  this->UpdateDiffusionRates(partition_sizes, total_global_size, diffusion_rates);

  std::vector<uint8_t> compact_diffusion_level(this->local_n);

  DiffusionGhostAsyncInfo diffusion_asyc_comm;
  diffusion_asyc_comm.requests = (MPI_Request*)malloc((this->ghost_procs_n + this->send_procs_n) * (sizeof(MPI_Request)));
  diffusion_asyc_comm.statuses = (MPI_Status*)malloc((this->ghost_procs_n + this->send_procs_n) * (sizeof(MPI_Status)));



  for (int d_i = 0; d_i < max_diffusion_iterations; d_i++)
  {
    // starting ghost sync
    this->StartReceivingDiffusionGhost(&diffusion_asyc_comm);
    this->DiffuseOnce(diffusion_vector, diffusion_rates);

    this->UpdateDiffusionLabels(diffusion_vector, partition_sizes, diffusion_rates,
                                local_parts_size_deltas, compact_diffusion_level);

    // finishing ghost sync
    for (size_t send_i = 0; send_i < this->send_n; send_i++)
    {
      this->diffusion_send_buffer[send_i] = diffusion_vector[this->send_scatter_map[send_i]];
    }
    this->EndReceivingDiffusionGhost(&diffusion_asyc_comm);
      // this->ExchangeAllGhosts(this->diffusion_send_buffer.data(), &diffusion_vector[this->local_n], DIFFUSION_GHOST_EXCHANGE_TAG);

    for (size_t recv_i = 0; recv_i < this->ghost_n; recv_i++)
    {
      diffusion_vector[this->local_n + recv_i] = this->diffusion_recv_buffer[recv_i];
    }

    // this->FillDiffusionCompactLevel2(compact_diffusion_level);

    // this->DiffuseCompact(diffusion_vector, diffusion_rates, compact_diffusion_level);


      // this->SyncNeighborPartSizes(partition_sizes, local_parts_size_deltas, level_1_2_neigh_procs);
      // this->SyncNeighborPartSizes(partition_sizes, local_parts_size_deltas, this->ghost_procs);

    this->SyncGlobalPartSizes(diffusion_vector, partition_sizes);

    this->UpdateDiffusionRates(partition_sizes, total_global_size,diffusion_rates );
  }

  this->diffusion_send_buffer.clear();
  this->diffusion_recv_buffer.clear();

  free(diffusion_asyc_comm.requests);
  free(diffusion_asyc_comm.statuses);


  
  
}


void amracut::DGraph::DiffuseOnce(std::vector<DiffusionValue> &diffusion_vector,
                                   const std::vector<amracut_real_t> &partition_diffusion_rates)
{
  std::vector<DiffusionValue> local_diffusion_vector_new(this->local_n);
  for (amracut_uint_t vertex = 0; vertex < this->local_n; vertex++)
  {
    amracut_uint_t part_v = diffusion_vector[vertex].label;

    amracut_real_t outgoing_flux = diffusion_vector[vertex].value * partition_diffusion_rates[part_v];

    amracut_real_t incoming_flux = 0;
    amracut_real_t edge_wgt_total = 0;


    for (graph_indexing_t neighbor_i = this->local_xadj[vertex];
          neighbor_i < this->local_xadj[vertex + 1]; neighbor_i++) // looping neighbors
    {
      auto neighbor = this->local_adjncy[neighbor_i];
      amracut_real_t edge_wgt = (this->wgtflag == AMRACUT_EDGE_WEIGHTED || this->wgtflag == AMRACUT_VTX_EDGE_WEIGHTED) ? 
                                  this->adjwgt_normalized[neighbor_i] : 1.0;

      amracut_uint_t part_neigh = diffusion_vector[neighbor].label;
      incoming_flux += diffusion_vector[neighbor].value * partition_diffusion_rates[part_neigh] * edge_wgt;
      edge_wgt_total+= edge_wgt;

    }
    outgoing_flux*= edge_wgt_total;
    local_diffusion_vector_new[vertex].label = part_v;    // no change to labels
    local_diffusion_vector_new[vertex].value = std::clamp(diffusion_vector[vertex].value + incoming_flux - outgoing_flux, 
                                                          (amracut_real_t)0.0, (amracut_real_t)1.0);
  }
  // copy the new updated values
  std::copy(local_diffusion_vector_new.begin(), local_diffusion_vector_new.end(), diffusion_vector.begin());

}                                   

void amracut::DGraph::UpdateDiffusionRates(const std::vector<amracut_uint_t> &est_part_sizes, 
                                            amracut_uint_t total_global_size,
                                            std::vector<amracut_real_t> &partition_diffusion_rates)
{ 
  auto ideal_part_size = total_global_size / this->procs_n;
  const float MAX_IMBALANCE = 1.1;

  std::transform(est_part_sizes.begin(), est_part_sizes.end(), partition_diffusion_rates.begin(),
                 [ideal_part_size, MAX_IMBALANCE](auto part_size)
                 {
                   if (part_size < 0.5 * ideal_part_size)
                   {
                     return static_cast<amracut_real_t>(0.01);
                   }
                   if (part_size < MAX_IMBALANCE * ideal_part_size)
                   {
                     return static_cast<amracut_real_t>(0.2);
                   }
                   return static_cast<amracut_real_t>(0.8);
                 });
}


void amracut::DGraph::UpdateDiffusionLabels(std::vector<DiffusionValue> &diffusion_vector,
                                             const std::vector<amracut_uint_t> &est_part_sizes,
                                             const std::vector<amracut_real_t> &partition_diffusion_rates,
                                             std::unordered_map<amracut_uint_t, amracut_sint_t> &local_p_size_deltas,
                                             std::vector<uint8_t> &compact_diffusion_level)
{

  // local_p_size_deltas.clear();
  // std::fill(compact_diffusion_level.begin(), compact_diffusion_level.end(), COMPACT_LEVEL_NONE);

  std::vector<DiffusionValue> local_diffusion_vector_new(this->local_n);


  for (amracut_uint_t vertex = 0; vertex < this->local_n; vertex++)
  {
    auto value_v = diffusion_vector[vertex].value;
    auto part_v = diffusion_vector[vertex].label;

    

    if (value_v < 0.5)    // TODO: make 0.5 configurable
    {
      // compact_diffusion_level[vertex] = COMPACT_LEVEL_0;

      amracut_uint_t best_part      = AMRACUT_BFS_NO_LABEL;
      amracut_real_t best_part_rate = std::numeric_limits<amracut_real_t>::max();

      amracut_real_t incoming_flux = 0;
      amracut_real_t edge_wgt_total = 0;

      for (graph_indexing_t neighbor_i = this->local_xadj[vertex];
            neighbor_i < this->local_xadj[vertex + 1]; neighbor_i++) // looping neighbors
      {
        auto neighbor   = this->local_adjncy[neighbor_i];
        auto part_neigh = diffusion_vector[neighbor].label;
        amracut_real_t edge_wgt = (this->wgtflag == AMRACUT_EDGE_WEIGHTED || this->wgtflag == AMRACUT_VTX_EDGE_WEIGHTED) ? 
                                    this->adjwgt_normalized[neighbor_i] : 1.0;

        incoming_flux += diffusion_vector[neighbor].value * partition_diffusion_rates[part_neigh] * edge_wgt;
        edge_wgt_total+= edge_wgt;

        if (partition_diffusion_rates[part_neigh] < best_part_rate)
        {
          best_part      = part_neigh;     // in the neighborhood, we assign to the partition with lowest diffusion rate
          best_part_rate = partition_diffusion_rates[part_neigh];
        } 

        // if (neighbor < this->local_n && compact_diffusion_level[neighbor] > COMPACT_LEVEL_1)
        // {
        //   compact_diffusion_level[neighbor] = COMPACT_LEVEL_1;
        // }
               
      }
      local_diffusion_vector_new[vertex].value = std::clamp(incoming_flux / (edge_wgt_total * best_part_rate), 
                                                          (amracut_real_t)0.0, (amracut_real_t)1.0);  
      local_diffusion_vector_new[vertex].label = best_part;  

      // updating partition size deltas
      // if (part_v != best_part)
      // {

      //   local_p_size_deltas.emplace(part_v, 0);
      //   local_p_size_deltas.emplace(best_part, 0);

      //   amracut_sint_t w = (this->wgtflag == AMRACUT_VTX_WEIGHTED || this->wgtflag == AMRACUT_VTX_EDGE_WEIGHTED)? 
      //                        this->vwgt[vertex] : 1;

      //   local_p_size_deltas[part_v]    -= w; 
      //   local_p_size_deltas[best_part] += w;
      // }

    } else
    {
      local_diffusion_vector_new[vertex].value = value_v;  
      local_diffusion_vector_new[vertex].label = part_v;  
    } 
  }

  // copy the new updated values
  std::copy(local_diffusion_vector_new.begin(), local_diffusion_vector_new.end(), diffusion_vector.begin());

}    

#if 0
void amracut::DGraph::FillDiffusionCompactLevel2(std::vector<uint8_t> &compact_level)
{
  for (amracut_uint_t vertex = 0; vertex < this->local_n; vertex++)
  {
    if (compact_level[vertex] != COMPACT_LEVEL_1)
    {
      continue;
    }
    

    for (graph_indexing_t neighbor_i = this->local_xadj[vertex];
          neighbor_i < this->local_xadj[vertex + 1]; neighbor_i++) // looping neighbors
    {
      auto neighbor = this->local_adjncy[neighbor_i];

      if (neighbor < this->local_n && compact_level[neighbor] > COMPACT_LEVEL_2)
      {
        compact_level[neighbor] = COMPACT_LEVEL_2;
      }
      
    }
  }
}
#endif

void amracut::DGraph::SyncGlobalPartSizes(const std::vector<DiffusionValue> &diffusion_vector,
                                           std::vector<amracut_uint_t> &global_partition_sizes)
{
  std::vector<amracut_uint_t> local_partition_sizes(this->procs_n, 0);

  for (size_t i = 0; i < this->local_n; i++)
  {
    auto label = diffusion_vector[i].label;
    if(this->wgtflag == AMRACUT_VTX_WEIGHTED || this->wgtflag == AMRACUT_VTX_EDGE_WEIGHTED)
    {
      local_partition_sizes[label]+= this->vwgt[i];
    } else
    {
      local_partition_sizes[label]++;
    }
    
  }


  MPI_Allreduce(local_partition_sizes.data(), global_partition_sizes.data(),
                procs_n, Mpi_datatype<amracut_uint_t>::value(), MPI_SUM, *this->comm);

}

#if 0
void amracut::DGraph::DiscoverLevel1and2NeighProcs(std::vector<int> &neighbor_procs)
{
  // first we get the count of ghost procs of each of my ghost procs (i.e. the count of level 2 neighbors)
  std::vector<int> recv_cnts(this->ghost_procs_n, 0);

  std::vector<MPI_Request> requests(2*this->ghost_procs_n);
  std::vector<MPI_Status> statuses(2*this->ghost_procs_n);

  for (int p_i = 0; p_i < this->ghost_procs_n; p_i++)
  {
    MPI_Irecv(&recv_cnts[p_i], 1, Mpi_datatype<int>::value(), this->ghost_procs[p_i], 
              LEVEL2_GHOST_PROCS_COUNTS_DISCOVER_TAG, *(this->comm), &(requests[p_i]));
  }

  for (int p_i = 0; p_i < this->ghost_procs_n; p_i++)
  {
    MPI_Isend(&(this->ghost_procs_n), 1, Mpi_datatype<int>::value(), this->ghost_procs[p_i], 
              LEVEL2_GHOST_PROCS_COUNTS_DISCOVER_TAG, *(this->comm), &(requests[p_i + this->ghost_procs_n]));
  }

  MPI_Waitall(2 * this->ghost_procs_n, requests.data(), statuses.data());

  // now we exchange the level 2 ghost proc ids

  std::vector<int> send_cnts(this->ghost_procs_n, this->ghost_procs_n);
  std::vector<int> send_displs(this->ghost_procs_n, 0);     // because we are sending the same buffer `this->ghost_procs` to all

  std::vector<int> recv_displs(this->ghost_procs_n, 0);
  
  std::exclusive_scan(recv_cnts.begin(), recv_cnts.end(), recv_displs.begin(), 0);

  std::vector<int> level2_procs(recv_cnts[ghost_procs_n-1] + recv_displs[ghost_procs_n-1]);

  this->MpiNeighborAllToAllVSparse(this->ghost_procs.data(), this->ghost_procs_n, this->ghost_procs.data(),
                                   send_cnts.data(), send_displs.data(), 
                                   level2_procs.data(), this->ghost_procs_n, this->ghost_procs.data(), 
                                   recv_cnts.data(), recv_displs.data(), LEVEL2_GHOST_PROCS_DISCOVER_TAG);

  // level2_procs includes my rank and potential duplicates, filter them
  std::unordered_set<int> level2_procs_filtered(level2_procs.begin(), level2_procs.end());
  level2_procs_filtered.erase(this->my_rank);

  // copy to output vector
  neighbor_procs.resize(level2_procs_filtered.size());
  std::copy(level2_procs_filtered.begin(), level2_procs_filtered.end(), neighbor_procs.begin());
}


void amracut::DGraph::SyncNeighborPartSizes(std::vector<amracut_uint_t> &est_part_sizes,
                                             const std::unordered_map<amracut_uint_t, amracut_sint_t> &local_p_size_deltas,
                                             std::vector<int> &neighbor_procs)
{
  int neighbor_procs_n = neighbor_procs.size();
  std::vector<PartAndDelta> local_nnz_parts_delta;

  for (const auto [ p, d ] : local_p_size_deltas)
  {
    if (d != 0)
    {
      local_nnz_parts_delta.push_back({.partition = p, .delta = d});  
    }  
  }
  int snd_cnt = local_nnz_parts_delta.size();

  std::vector<int> rcv_cnts(neighbor_procs_n, 0);

  // first we exchange counts
  {
    std::vector<MPI_Request> requests(2 * neighbor_procs_n);
    std::vector<MPI_Status> statuses(2 * neighbor_procs_n);

    int i_= 0;

    for (int p_i = 0; p_i < neighbor_procs_n; p_i++)
    {
      MPI_Irecv(&rcv_cnts[p_i], 1, Mpi_datatype<int>::value(), neighbor_procs[p_i], 
                PART_DELTA_COUNT_EXCHANGE_TAG, *(this->comm), &(requests[i_++]));
    }

    for (int p_i = 0; p_i < neighbor_procs_n; p_i++)
    {
      MPI_Isend(&snd_cnt, 1, Mpi_datatype<int>::value(), neighbor_procs[p_i], 
                PART_DELTA_COUNT_EXCHANGE_TAG, *(this->comm), &(requests[i_++]));
    }

    MPI_Waitall(2 * neighbor_procs_n, requests.data(), statuses.data());
  } // count exchange done 


  // now we exchnage delta values

  std::vector<int> snd_cnts(neighbor_procs_n, snd_cnt);
  std::vector<int> snd_displs(neighbor_procs_n, 0);   // since we are sending the same thing for all neighbor_procs

  std::vector<int> recv_displs(neighbor_procs_n);



  std::exclusive_scan(rcv_cnts.begin(), rcv_cnts.end(), recv_displs.begin(), 0);

  std::vector<PartAndDelta> recv_parts_delta(rcv_cnts[neighbor_procs_n-1] + recv_displs[neighbor_procs_n-1]);


  this->MpiNeighborAllToAllVSparse<PartAndDelta>(local_nnz_parts_delta.data(), neighbor_procs_n, neighbor_procs.data(),
                                   snd_cnts.data(), snd_displs.data(), 
                                   recv_parts_delta.data(), neighbor_procs_n, neighbor_procs.data(), 
                                   rcv_cnts.data(), recv_displs.data(), PART_DELTA_EXCHANGE_TAG);


  // now update partition size estimates
  for (const auto& part_delta : local_nnz_parts_delta)
  {
    est_part_sizes[part_delta.partition] = static_cast<amracut_uint_t>(
                                                                        static_cast<amracut_sint_t>(est_part_sizes[part_delta.partition]) 
                                                                        + part_delta.delta
                                                                       );
  }
  
  for (const auto& part_delta : recv_parts_delta)
  {
    est_part_sizes[part_delta.partition] = static_cast<amracut_uint_t>(
                                                                        static_cast<amracut_sint_t>(est_part_sizes[part_delta.partition]) 
                                                                        + part_delta.delta
                                                                       );
  }
}                                             
#endif

void amracut::DGraph::StartReceivingDiffusionGhost(DiffusionGhostAsyncInfo* comm_info)
{
  for (int i = 0; i < this->ghost_procs_n; i++)
  {
    MPI_Irecv(&(this->diffusion_recv_buffer[this->ghost_counts_per_proc_scanned[i]]), 
                this->ghost_counts_per_proc[i], Mpi_datatype<DiffusionValue>::value(), this->ghost_procs[i],
                DIFFUSION_GHOST_EXCHANGE_TAG, *(this->comm), &(comm_info->requests[i]));
  }
}

void amracut::DGraph::EndReceivingDiffusionGhost(DiffusionGhostAsyncInfo* comm_info)
{
  // since the first part of MPI req and status buffers contain the recev requests (from ghost procs), we have to offset that
  int mpi_idx_offest = this->ghost_procs_n;

  for (int i = 0; i < this->send_procs_n; i++)
  {
    MPI_Isend(&(this->diffusion_send_buffer[this->send_counts_per_proc_scanned[i]]), 
                this->send_counts_per_proc[i], Mpi_datatype<DiffusionValue>::value(), this->send_procs[i],
                DIFFUSION_GHOST_EXCHANGE_TAG, *(this->comm), &(comm_info->requests[mpi_idx_offest + i]));
  }

  MPI_Waitall(this->ghost_procs_n + this->send_procs_n, comm_info->requests, comm_info->statuses);
}

#if 0
void amracut::DGraph::DiffuseCompact(std::vector<DiffusionValue> &diffusion_vector,
                                      const std::vector<amracut_real_t> &partition_diffusion_rates,
                                      const std::vector<uint8_t> &compact_diffusion_level)
{
  std::vector<uint8_t> levels = {this->COMPACT_LEVEL_0, this->COMPACT_LEVEL_1, this->COMPACT_LEVEL_2};
  

  for (int level_i = 0; level_i < 3; level_i++)
  {
    uint8_t level = levels[level_i];
    std::vector<DiffusionValue> local_diffusion_vector_new(diffusion_vector.begin(), diffusion_vector.begin() + this->local_n);
    
    for (amracut_uint_t vertex = 0; vertex < this->local_n; vertex++)
    {
      if (compact_diffusion_level[vertex] > level)
      {
        continue;
      }
      

      amracut_uint_t degree_v = this->local_xadj[vertex + 1] - this->local_xadj[vertex];
      amracut_uint_t part_v = diffusion_vector[vertex].label;

      amracut_real_t dv_dt = -diffusion_vector[vertex].value * degree_v * partition_diffusion_rates[part_v];

      for (graph_indexing_t neighbor_i = this->local_xadj[vertex];
            neighbor_i < this->local_xadj[vertex + 1]; neighbor_i++) // looping neighbors
      {
        auto neighbor = this->local_adjncy[neighbor_i];
        amracut_uint_t part_neigh = diffusion_vector[neighbor].label;
        dv_dt += diffusion_vector[neighbor].value * partition_diffusion_rates[part_neigh];

      }
      local_diffusion_vector_new[vertex].label = part_v;    // no change to labels
      local_diffusion_vector_new[vertex].value = std::clamp(diffusion_vector[vertex].value + dv_dt, 
                                                            (amracut_real_t)0.0, (amracut_real_t)1.0);
    }
    // copy the new updated values
    std::copy(local_diffusion_vector_new.begin(), local_diffusion_vector_new.end(), diffusion_vector.begin());

    // sync ghost
    for (size_t send_i = 0; send_i < this->send_n; send_i++)
    {
      this->diffusion_send_buffer[send_i] = diffusion_vector[this->send_scatter_map[send_i]];
    }
    this->ExchangeAllGhosts(this->diffusion_send_buffer.data(), &diffusion_vector[this->local_n], COMPACT_DIFFUSION_GHOST_EXCHANGE_TAG);

  }
}     
#endif




void amracut::DGraph::StartReceivingUpdatedOnlyGhostCounts()
{

  std::fill(this->updated_ghost_counts_per_proc.begin(), this->updated_ghost_counts_per_proc.end(), 0);
  for (int i = 0; i < this->ghost_procs_n; i++)
  {
    MPI_Irecv(&(this->updated_ghost_counts_per_proc[i]), 1, MPI_INT, this->ghost_procs[i],
              this->GHOST_COUNT_EXCHANGE_TAG, *(this->comm), &(this->updated_count_mpi_requests[i]));
    this->recv_size += sizeof(int);
  }
}

void amracut::DGraph::EndReceivingUpdatedOnlyGhostCounts()
{

  // since the first part of MPI req and status buffers contain the recev requests (from ghost procs), we have to offset that
  int mpi_idx_offest = this->ghost_procs.size();

  for (int i = 0; i < this->send_procs_n; i++)
  {
    MPI_Isend(&(this->updated_send_counts_per_proc[i]), 1, MPI_INT, this->send_procs[i],
              this->GHOST_COUNT_EXCHANGE_TAG, *(this->comm), &(this->updated_count_mpi_requests[i + mpi_idx_offest]));
    this->send_size += sizeof(int);
  }


  MPI_Waitall(this->ghost_procs_n + this->send_procs_n, this->updated_count_mpi_requests.data(), this->updated_count_mpi_statuses.data());

}


void amracut::DGraph::ExchangeUpdatedOnlyGhost(const std::vector<BFSValue> &ghost_send_buffer,
                                                const std::vector<BFSValue> &ghost_send_buffer_prev,
                                                std::vector<BFSValue> &ghost_recv_buffer)
{
  std::fill(this->updated_send_counts_per_proc.begin(), this->updated_send_counts_per_proc.end(), 0);

  int slot = 0;
  for (int i = 0; i < this->send_procs_n; i++)
  {
    for (int send_offset = 0; send_offset < this->send_counts_per_proc[i]; send_offset++)
    {
      int idx = send_offset + this->send_counts_per_proc_scanned[i];

      if (ghost_send_buffer[idx].distance != ghost_send_buffer_prev[idx].distance ||
          ghost_send_buffer[idx].label != ghost_send_buffer_prev[idx].label) // updated only check
      {
        this->updated_send_counts_per_proc[i]++;
        this->updated_send_buffer[slot++] = {.label = ghost_send_buffer[idx].label,
                                             .distance = ghost_send_buffer[idx].distance,
                                             .offset = static_cast<graph_indexing_t>(send_offset)};
      }
    }
  }

  std::exclusive_scan(this->updated_send_counts_per_proc.begin(), this->updated_send_counts_per_proc.end(),
                      this->updated_send_counts_per_proc_scanned.begin(), 0);

  this->EndReceivingUpdatedOnlyGhostCounts();

  std::exclusive_scan(this->updated_ghost_counts_per_proc.begin(), this->updated_ghost_counts_per_proc.end(),
                      this->updated_ghost_counts_per_proc_scanned.begin(), 0);

  this->MpiNeighborAllToAllVSparse(updated_send_buffer.data(), this->send_procs_n, this->send_procs.data(),
                                   updated_send_counts_per_proc.data(), updated_send_counts_per_proc_scanned.data(),
                                   updated_ghost_buffer.data(), this->ghost_procs_n, this->ghost_procs.data(),
                                   updated_ghost_counts_per_proc.data(), updated_ghost_counts_per_proc_scanned.data(),
                                   BFS_GHOST_EXCHANGE_TAG);

  // now place the received updated values in correct places in the original receive buffer

  for (int i = 0; i < this->ghost_procs_n; i++)
  {

    for (int received_i = 0; received_i < updated_ghost_counts_per_proc[i]; received_i++)
    {
      int idx_in_updated_ghost_buffer = updated_ghost_counts_per_proc_scanned[i] + received_i;
      GhostBFSValue received_value = updated_ghost_buffer[idx_in_updated_ghost_buffer];

      int idx_in_original_recv_buffer = static_cast<int>(received_value.offset) + this->ghost_counts_per_proc_scanned[i];

      ghost_recv_buffer[idx_in_original_recv_buffer].label = received_value.label;
      ghost_recv_buffer[idx_in_original_recv_buffer].distance = received_value.distance;
    }
  }
}

template <typename T>
void amracut::DGraph::MpiNeighborAllToAllVSparse(T *sendbuf, int send_p_cnt, int *send_p, int *sendcnts, int *sdispls,
                                                  T *recvbuf, int recv_p_cnt, int *recv_p, int *recvcnts, int *rdispls,
                                                  int tag)
{

  int nnz_p_count = 0;

  for (int i = 0; i < send_p_cnt; i++)
  {
    if (sendcnts[i] > 0)
      nnz_p_count++;
  }

  for (int i = 0; i < recv_p_cnt; i++)
  {
    if (recvcnts[i] > 0)
      nnz_p_count++;
  }

  std::vector<MPI_Request> requests(nnz_p_count);
  std::vector<MPI_Status> statuses(nnz_p_count);

  int proc_i = 0;

  // First place all recv requests.
  for (int recv_i = 0; recv_i < recv_p_cnt; recv_i++)
  {
    if (recvcnts[recv_i] == 0)
      continue;

    MPI_Irecv(&(recvbuf[rdispls[recv_i]]), recvcnts[recv_i], Mpi_datatype<T>::value(), recv_p[recv_i], tag,
              *this->comm, &(requests[proc_i]));
    proc_i++;
    this->recv_size += sizeof(T) * recvcnts[recv_i];
  }

  // Next send the messages.
  for (int send_i = 0; send_i < send_p_cnt; send_i++)
  {
    if (sendcnts[send_i] == 0)
      continue;

    MPI_Isend(&(sendbuf[sdispls[send_i]]), sendcnts[send_i], Mpi_datatype<T>::value(), send_p[send_i], tag,
              *this->comm, &(requests[proc_i]));
    proc_i++;
    this->send_size += sizeof(T) * sendcnts[send_i];
  }

  MPI_Waitall(nnz_p_count, requests.data(), statuses.data());

  return;
}

template <typename T>
void amracut::DGraph::ExchangeAllGhosts(T* ghost_send_buffer, T* ghost_recv_buffer, int tag)
{
  this->MpiNeighborAllToAllVSparse(ghost_send_buffer, this->send_procs_n, this->send_procs.data(),
                                   this->send_counts_per_proc.data(), this->send_counts_per_proc_scanned.data(), 
                                   ghost_recv_buffer, this->ghost_procs_n, this->ghost_procs.data(), 
                                   this->ghost_counts_per_proc.data(), this->ghost_counts_per_proc_scanned.data(),
                                   tag);
  
}