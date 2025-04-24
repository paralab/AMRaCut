/**
 * @file    dgraph.hpp
 * @author  
 * @brief   Interface for distributed graph data structure class
 * @version 0.1
 * @date 2024-10-25
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef AMRACUT_DGRAPH_H__
#define AMRACUT_DGRAPH_H__

#include <vector>
#include <array>
#include "utils.hpp" 
#include <mpi.h>
#include <chrono>
#include <unordered_map>


// #if BFS_DISTANCE_TYPE == 64
// using bfs_distance_t = uint64_t;
// #define DIST_GRAPH_BFS_INFINITY UINT64_MAX

// #elif BFS_DISTANCE_TYPE == 32
// using bfs_distance_t = uint32_t;
// #define DIST_GRAPH_BFS_INFINITY UINT32_MAX

// #elif BFS_DISTANCE_TYPE == 16
// using bfs_distance_t = uint16_t;
// #define DIST_GRAPH_BFS_INFINITY UINT16_MAX

// #else
// #error "Invalid BFS_DISTANCE_TYPE specified. Allowed values: 16, 32, 64"
// #endif

// #if BFS_LABEL_INTEGER_WIDTH == 64
// using bfs_label_t = uint64_t;
// #define DIST_GRAPH_BFS_NO_LABEL UINT64_MAX

// #elif BFS_LABEL_INTEGER_WIDTH == 32
// using bfs_label_t = uint32_t;
// #define DIST_GRAPH_BFS_NO_LABEL UINT32_MAX

// #elif BFS_LABEL_INTEGER_WIDTH == 16
// using bfs_label_t = uint16_t;
// #define DIST_GRAPH_BFS_NO_LABEL UINT16_MAX

// #else
// #error "Invalid BFS_LABEL_INTEGER_WIDTH specified. Allowed values: 16, 32, 64"
// #endif

namespace amracut
{
  struct BFSValue
  {
    bfs_label_t     label;
    bfs_distance_t  distance;
  };

  

  template <>
  class Mpi_datatype<BFSValue>    // MPI_Datatype for the C++ datatype "BFSValue"
  {
  public:
    static MPI_Datatype value()
    {
      static bool first = true;
      static MPI_Datatype custom_mpi_type;

      if (first)
      {
        first = false;
        int block_lengths[2] = {1, 1};
        MPI_Datatype types[2] = {Mpi_datatype<bfs_label_t>::value(), Mpi_datatype<bfs_distance_t>::value()};
        MPI_Aint offsets[2];
        offsets[0] = offsetof(BFSValue, label);
        offsets[1] = offsetof(BFSValue, distance);

        MPI_Type_create_struct(2, block_lengths, offsets, types, &custom_mpi_type);
        MPI_Type_commit(&custom_mpi_type);
      }

      return custom_mpi_type;
    }
  };

  struct DiffusionValue
  {
    bfs_label_t         label;
    amracut_real_t     value;
  };

  template <>
  class Mpi_datatype<DiffusionValue>    // MPI_Datatype for the C++ datatype "DiffusionValue"
  {
  public:
    static MPI_Datatype value()
    {
      static bool first = true;
      static MPI_Datatype custom_mpi_type;

      if (first)
      {
        first = false;
        int block_lengths[2] = {1, 1};
        MPI_Datatype types[2] = {Mpi_datatype<bfs_label_t>::value(), Mpi_datatype<amracut_real_t>::value()};
        MPI_Aint offsets[2];
        offsets[0] = offsetof(DiffusionValue, label);
        offsets[1] = offsetof(DiffusionValue, value);

        MPI_Type_create_struct(2, block_lengths, offsets, types, &custom_mpi_type);
        MPI_Type_commit(&custom_mpi_type);
      }

      return custom_mpi_type;
    }
  };

  /**
   * To be used for BFS value exchange for updated only ghost vertices
   */
  struct GhostBFSValue
  {
    bfs_label_t       label;
    bfs_distance_t    distance;
    graph_indexing_t  offset;       // offset w.r.t original send buffer
  };

  template <>
  class Mpi_datatype<GhostBFSValue>     // MPI_Datatype for the C++ datatype "GhostBFSValue"
  {
  public:
    static MPI_Datatype value()
    {
      static bool first = true;
      static MPI_Datatype custom_mpi_type;

      if (first)
      {
        first = false;
        int block_lengths[3] = {1, 1, 1};
        MPI_Datatype types[3] = {Mpi_datatype<bfs_label_t>::value(), Mpi_datatype<bfs_distance_t>::value(), Mpi_datatype<graph_indexing_t>::value()};
        MPI_Aint offsets[3];
        offsets[0] = offsetof(GhostBFSValue, label);
        offsets[1] = offsetof(GhostBFSValue, distance);
        offsets[2] = offsetof(GhostBFSValue, offset);

        MPI_Type_create_struct(3, block_lengths, offsets, types, &custom_mpi_type);
        MPI_Type_commit(&custom_mpi_type);
      }

      return custom_mpi_type;
    }
  };

  struct PartAndDelta
  {
    bfs_label_t         partition;
    amracut_sint_t     delta;
  };

  template <>
  class Mpi_datatype<PartAndDelta>    // MPI_Datatype for the C++ datatype "PartAndDelta"
  {
  public:
    static MPI_Datatype value()
    {
      static bool first = true;
      static MPI_Datatype custom_mpi_type;

      if (first)
      {
        first = false;
        int block_lengths[2] = {1, 1};
        MPI_Datatype types[2] = {Mpi_datatype<bfs_label_t>::value(), Mpi_datatype<amracut_sint_t>::value()};
        MPI_Aint offsets[2];
        offsets[0] = offsetof(PartAndDelta, partition);
        offsets[1] = offsetof(PartAndDelta, delta);

        MPI_Type_create_struct(2, block_lengths, offsets, types, &custom_mpi_type);
        MPI_Type_commit(&custom_mpi_type);
      }

      return custom_mpi_type;
    }
  };

  struct DiffusionGhostAsyncInfo
  {
    MPI_Request*  requests;
    MPI_Status*   statuses;
  };
  



  /**
   * Used when initializing DGraph
   */
  struct GhostEdge
  {
    amracut_uint_t  local_vertex;             // in local indexing
    amracut_uint_t  ghost_vertex;             // in global indexing
    amracut_uint_t  ghost_local_idx;          
    amracut_uint_t  wgt;
    amracut_uint_t  adjncy_idx;               // index of this edge in adjncy
    int              ghost_proc;               // process which the ghost vertex belongs to

  };


  enum MPI_COMM_TAGS {
    PART_DELTA_EXCHANGE_TAG               = 100,
    PART_DELTA_COUNT_EXCHANGE_TAG         = 101,
    BFS_GHOST_EXCHANGE_TAG                = 102,
    DIFFUSION_GHOST_EXCHANGE_TAG          = 103,
    COMPACT_DIFFUSION_GHOST_EXCHANGE_TAG  = 104,
    LEVEL2_GHOST_PROCS_COUNTS_DISCOVER_TAG= 105,
    LEVEL2_GHOST_PROCS_DISCOVER_TAG       = 106

    



  }; 

  /**
   * @brief Distributed graph data structure
   * 
   */
  class DGraph
  {
  private:
    MPI_Comm*         comm;
    int               procs_n, my_rank;
    amracut_uint_t    global_n;
    amracut_uint_t    local_n;
    amracut_uint_t    wgtflag;
  
    
    std::vector<graph_indexing_t>    local_xadj;     // CSR format pointers
    std::vector<graph_indexing_t>    local_adjncy;   // CSR format values
    std::vector<amracut_uint_t>     vtx_dist;
    std::vector<amracut_uint_t>     vwgt;
    std::vector<amracut_uint_t>     adjwgt;
    std::vector<amracut_real_t>     adjwgt_normalized;



    bool has_ghost;

    amracut_uint_t        ghost_n;
    int                   ghost_procs_n;
    std::vector<int>      ghost_procs;
    std::vector<int>      ghost_counts_per_proc;
    std::vector<int>      ghost_counts_per_proc_scanned;
    // std::vector<BFSValue> recev_buffer;

    std::vector<int>            updated_ghost_counts_per_proc;
    std::vector<int>            updated_ghost_counts_per_proc_scanned;
    std::vector<MPI_Request>    updated_count_mpi_requests;      // recev + send
    std::vector<MPI_Status>     updated_count_mpi_statuses;      // recev + send
    std::vector<GhostBFSValue>  updated_ghost_buffer;

    amracut_uint_t              send_n;
    int                         send_procs_n;
    std::vector<int>            send_procs;
    std::vector<int>            send_counts_per_proc;
    std::vector<int>            send_counts_per_proc_scanned;
    // std::vector<BFSValue>       send_buffer;
    std::vector<amracut_uint_t> send_scatter_map;

    std::vector<DiffusionValue> diffusion_send_buffer;
    std::vector<DiffusionValue> diffusion_recv_buffer;


    std::vector<int>            updated_send_counts_per_proc;
    std::vector<int>            updated_send_counts_per_proc_scanned;
    std::vector<GhostBFSValue>  updated_send_buffer;

    //temp buffers for bfs
    struct BFSCandidate
    {
      graph_indexing_t  local_idx;
      BFSValue          bfs_value;
    };
    
    std::vector<graph_indexing_t>   frontier_buffer;
    std::vector<uint8_t>            is_candidate;
    std::vector<BFSCandidate>       candidates;



    // for optimized async operations of exchanging updated only ghost counts
    static constexpr int GHOST_VALUE_EXCHANGE_TAG = 2;
    static constexpr int GHOST_COUNT_EXCHANGE_TAG = 3;

    std::vector<MPI_Request> comm_requests; // recev + send
    std::vector<MPI_Status> comm_statuses;  // recev + send

    std::vector<int> updated_only_recv_counts;

    static constexpr uint8_t COMPACT_LEVEL_0    = (uint8_t)0;
    static constexpr uint8_t COMPACT_LEVEL_1    = (uint8_t)1;
    static constexpr uint8_t COMPACT_LEVEL_2    = (uint8_t)2;
    static constexpr uint8_t COMPACT_LEVEL_NONE = std::numeric_limits<uint8_t>::max();




    // Adjust related constants
    static constexpr int alpha_count = 6;
    static constexpr std::array<float,  alpha_count>  alpha_vals    = {3.2, 2.5, 2.0, 1.7, 1.5, 1.2};
    static constexpr std::array<int,    alpha_count>  d_plus_vals   = {6, 5, 4, 3, 2, 1};

    static constexpr int beta_count = 4;
    static constexpr std::array<float,  beta_count>   beta_vals     = {0.2, 0.4, 0.6, 0.8};
    static constexpr std::array<int,    beta_count>   d_minus_vals  = {-4, -3, -2, -1};

    // stats related
    int verbose;
    std::chrono::microseconds com_duration;
    size_t recv_size;
    size_t send_size;



    /**
     * @brief Runs the first iteration of local BFS. 
     * 
     * Seed will be calculated as the middle vertex and label will MPI_rank.
     * 
     * Does not update the BFS status of ghost vertices.
     * 
     * @param[in,out] bfs_vector BFS status vector of local and ghost vertices
     */
    void RunFirstBFSIteration(std::vector<BFSValue> & bfs_vector);

    /**
     * @brief Continues local BFS iterations from last updated ghost vertices.
     * 
     * Does not update the BFS status of ghost vertices.
     * 
     * @param bfs_vector            BFS state vector of local and ghost vertices
     * @param updated_ghosts        Indices of updated ghost vertices
     * @param updated_ghosts_n      Number of updated ghost vertices
     * @return true                 if BFS continued
     * @return false                if BFS is already stable
     */
    bool RunBFSFromGhostUpdates(std::vector<BFSValue> &bfs_vector,
                                std::vector<graph_indexing_t> &updated_ghosts, amracut_uint_t updated_ghosts_n);

    /**
     * @brief Continues local BFS iteration from a given BFS state
     * 
     * Does not update the BFS status of ghost vertices.
     * 
     * @param bfs_vector            BFS state vector of local and ghost vertices
     * @return true                 if BFS continued
     * @return false                if BFS is already stable
     */
    bool RunBFS(std::vector<BFSValue> &bfs_vector);

    /**
     * @brief Apply the BFS distance perturbation by comparing the current partition sizes.
     * 
     * Can be called only after the BFS is in global stable state.
     * 
     * @param bfs_vector 
     */
    void Adjust(std::vector<BFSValue> &bfs_vector);

    /**
     * @brief Refines a given partition labelling by performing graph laplacian type diffusion on `diffusion_vector`
     * 
     * @param diffusion_vector    Should be an already initialized diffusion vector
     */
    void RefineByDiffusion(std::vector<DiffusionValue> &diffusion_vector);


    void DiffuseOnce(std::vector<DiffusionValue> &diffusion_vector,
                     const std::vector<amracut_real_t> &partition_diffusion_rates);

    void DiffuseCompact(std::vector<DiffusionValue> &diffusion_vector,
                        const std::vector<amracut_real_t> &partition_diffusion_rates,
                        const std::vector<uint8_t> &compact_diffusion_level);
    void FillDiffusionCompactLevel2(std::vector<uint8_t> &compact_level);

    void SyncNeighborPartSizes(std::vector<amracut_uint_t> &est_part_sizes,
                               const std::unordered_map<amracut_uint_t, amracut_sint_t> &local_p_size_deltas,
                               std::vector<int> &neighbor_procs);
    void DiscoverLevel1and2NeighProcs(std::vector<int> &neighbor_procs);
    void SyncGlobalPartSizes(const std::vector<DiffusionValue> &diffusion_vector,
                             std::vector<amracut_uint_t> &global_partition_sizes);                            

    /**
     * @brief Peform a MPI_Alltoallv sparse operation only for a given set of processes.
     * 
     * @tparam T          Data type to send. `amracut::Mpi_datatype` should be instantiated with `T`.
     * @param sendbuf     Send buffer
     * @param send_p_cnt  No. of destinations
     * @param send_p      Destination processes
     * @param sendcnts    Send counts
     * @param sdispls     Send displacements in `sendbuf`
     * @param recvbuf     Receiving buffer
     * @param recv_p_cnt  No. of sources
     * @param recv_p      Source processes
     * @param recvcnts    Receving counts
     * @param rdispls     Receiving displacements in `recvbuf`
     */
    template <typename T>
    void MpiNeighborAllToAllVSparse(T *sendbuf, int send_p_cnt, int *send_p, int *sendcnts, int *sdispls,
                                    T *recvbuf, int recv_p_cnt, int *recv_p, int *recvcnts, int *rdispls,
                                    int tag);
    /**
     * @brief Exchange status of all ghost verices
     * 
     * @param ghost_send_buffer   send buffer
     * @param ghost_recv_buffer   receive buffer 
     */
    template <typename T>
    void ExchangeAllGhosts(T* ghost_send_buffer,
                           T* ghost_recv_buffer,
                           int tag);

    /**
     * @brief Start async receiving updated ghost counts.
     * 
     * Starts updating `this->updated_ghost_counts_per_proc`
     * 
     * @return * void 
     */
    void StartReceivingUpdatedOnlyGhostCounts();

    /**
     * @brief End async receiving updated ghost counts.
     * 
     * Counterpart (i.e. the matching end) procedure of StartReceivingUpdatedOnlyGhostCounts
     * 
     * Sends counts from `this->updated_send_counts_per_proc`.
     * Ends updating `this->updated_ghost_counts_per_proc`
     * 
     */
    void EndReceivingUpdatedOnlyGhostCounts();

    /**
     * @brief Exchange BFS status of updated only ghost vertices. 
     * 
     * Checks the previous sent buffer and current sending buffer, and sends only updated values.
     * `ghost_recv_buffer` is updated only for the received values, others are unchanged
     *
     * 
     * `this->StartReceivingUpdatedOnlyGhostCounts` should be called before calling this function.
     * Internally calls `this->EndReceivingUpdatedOnlyGhostCounts`.
     * 
     * @param ghost_send_buffer 
     * @param ghost_send_buffer_prev 
     * @param ghost_recv_buffer 
     */
    void ExchangeUpdatedOnlyGhost(const std::vector<BFSValue> &ghost_send_buffer,
                                  const std::vector<BFSValue> &ghost_send_buffer_prev,
                                  std::vector<BFSValue> &ghost_recv_buffer);

    void StartReceivingDiffusionGhost(DiffusionGhostAsyncInfo* comm_info);
    void EndReceivingDiffusionGhost(DiffusionGhostAsyncInfo* comm_info);


    void UpdateDiffusionRates(const std::vector<amracut_uint_t> &est_part_sizes, 
                              amracut_uint_t total_global_size,
                              std::vector<amracut_real_t> &partition_diffusion_rates);
    void UpdateDiffusionLabels(std::vector<DiffusionValue> &diffusion_vector,
                               const std::vector<amracut_uint_t> &est_part_sizes,
                               const std::vector<amracut_real_t> &partition_diffusion_rates,
                               std::unordered_map<amracut_uint_t, amracut_sint_t> &local_p_size_deltas,
                               std::vector<uint8_t> &compact_diffusion_level);

  public:
    /**
     * @brief Construct a new DGraph object for distributed graph representation
     * 
     * @param vtx_dist  process P_i stores the vertices from vtxdist[i] up to (but not including) vertex vtxdist[i + 1]
     * @param xadj      local graph structure CSR encoding pointers
     * @param adjncy    local graph structure CSR encoding values
     * @param vwgt      Local vertex weights
     * @param adjwgt_   Edge weights
     * @param wgtflag   Flag to indicate if the graph is weighted.
     * @param comm      MPI communicator
     */
    DGraph(const amracut_uint_t *vtx_dist_, const amracut_uint_t *xadj_, const amracut_uint_t *adjncy_, 
           const amracut_uint_t *vwgt_, const amracut_uint_t *adjwgt_, const amracut_uint_t wgtflag_, 
           MPI_Comm *comm_);
    std::string PrintLocal();

    /**
     * @brief Compute graph partitioning
     * 
     * @param[out]  partition_labels_out  Computed partition labels in range [0, p -1]
     * @param[in]   verbose_               verbosity level for performance stats
     * @return      amracut_uint_t        Partiitoning status. 0 if suceess.
     */
    amracut_uint_t PartGraph(amracut_uint_t *partition_labels_out, bool use_diffusion, int verbose_);
    // ~DistGraph();
  };

} // namespace amracut


std::ostream &operator<<(std::ostream &os, const amracut::BFSValue &obj);
std::ostream &operator<<(std::ostream &os, const amracut::GhostBFSValue &obj);
std::ostream &operator<<(std::ostream &os, const amracut::DiffusionValue &obj);

#endif /* AMRACUT_DGRAPH_H__ */