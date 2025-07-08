#include "amracut.h"
#include "dgraph.hpp"
#include "utils.hpp"
#include "octree.hpp"
#include <cassert>
#include <stdexcept>

amracut_uint_t amracut_setup(amracut_ctrl *ctrl, const amracut_uint_t *vtx_dist,
                             const amracut_uint_t *xadj, const amracut_uint_t *adjncy,
                             const amracut_uint_t *vwgt, const amracut_uint_t *adjwgt, 
                             const amracut_uint_t wgtflag, MPI_Comm *comm)
{          
  if (!(wgtflag == AMRACUT_UNWEIGHTED || wgtflag == AMRACUT_VTX_WEIGHTED || 
        wgtflag == AMRACUT_EDGE_WEIGHTED|| wgtflag == AMRACUT_VTX_EDGE_WEIGHTED))
  {
    throw std::runtime_error("amracut ERROR: wgtflag should be one of 0, 1, 2, 3");
  }

  if ((wgtflag == AMRACUT_VTX_WEIGHTED || wgtflag == AMRACUT_VTX_EDGE_WEIGHTED) && vwgt == NULL)
  {
    throw std::runtime_error("amracut ERROR: vwgt can not be NULL when wgtflag is set to 1 or 3\n");
  }
  if ((wgtflag == AMRACUT_EDGE_WEIGHTED || wgtflag == AMRACUT_VTX_EDGE_WEIGHTED) && adjwgt == NULL)
  {
    throw std::runtime_error("amracut ERROR: adjwgt can not be NULL when wgtflag is set to 2 or 3\n");
  }
                        
  ctrl->internal_ = new amracut::DGraph(vtx_dist, xadj, adjncy, vwgt, adjwgt, wgtflag, comm);

  return 0;
}

amracut_uint_t amracut_partgraph(amracut_ctrl *ctrl, amracut_uint_t *parts, bool use_diffusion, int verbose)
{

  amracut::DGraph *dgraph = static_cast<amracut::DGraph *>(ctrl->internal_);
  dgraph->PartGraph(parts, use_diffusion, verbose);

  return 0;
}

amracut_uint_t amracut_destroy(amracut_ctrl *ctrl)
{
  assert(ctrl->internal_ != nullptr);
  amracut::DGraph *dgraph = static_cast<amracut::DGraph *>(ctrl->internal_);
  delete dgraph;

  return 0;
}


amracut_uint_t amracut_partgraph_octree( const amracut_uint_t *vtx_dist,
                                          const oct_element* local_elements, amracut_uint_t *parts, 
                                          MPI_Comm *comm)
{
  amracut::DGraph* dgraph;
  amracut::dgraph_from_octree(vtx_dist, local_elements, comm, &dgraph);
  dgraph->PartGraph(parts, true, 1);
  delete dgraph;
  return 0;
}                                          