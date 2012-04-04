#include "DiLeptonAnalysis/NTupleProducer/interface/ETHVertexInfo.h"

ETHVertexInfo::ETHVertexInfo(int nvtx, 
                             std::vector<float> vtxx, std::vector<float> vtxy, std::vector<float> vtxz, 
                             int ntracks, float * tkpx, float * tkpy, float * tkpz,
                             float * tkPtErr, int * tkVtxId,
                             float * tkd0, float * tkd0Err, float * tkdz, float * tkdzErr,
			     bool * tkIsHighPurity, std::vector<std::vector<unsigned short> > vtx_std_tkind, std::vector<std::vector<float> > vtx_std_tkweight, int * vtx_std_ntks
                             ) :
  nvtx_(nvtx),
  ntracks_(ntracks),
  tkpx_(tkpx),
  tkpy_(tkpy),
  tkpz_(tkpz),
  tkPtErr_(tkPtErr),
  tkVtxId_(tkVtxId),
  tkd0_(tkd0),
  tkd0Err_(tkd0Err),
  tkdz_(tkdz),
  tkdzErr_(tkdzErr),
  tkIsHighPurity_(tkIsHighPurity),
  vtx_std_tkind_(vtx_std_tkind),
  vtx_std_tkweight_(vtx_std_tkweight),
  vtx_std_ntks_(vtx_std_ntks)
{
  // Copy input vectors into local arrays
  assert(vtxx.size()==vtxy.size());
  assert(vtxx.size()==vtxz.size());
  for (size_t i=0; i<vtxx.size(); ++i) {
    vtxx_[i] = vtxx.at(i);
    vtxy_[i] = vtxy.at(i);
    vtxz_[i] = vtxz.at(i);
  }

    
  for (unsigned int i=0; i<vtx_std_tkind_.size();i++){
    for (unsigned int j=0; j<vtx_std_tkind_.at(i).size();j++){
      vtx_std_tkind_helper_[i][j]=vtx_std_tkind_.at(i).at(j);
    }
  }
    
  for (unsigned int i=0; i<vtx_std_tkweight_.size();i++){
    for (unsigned int j=0; j<vtx_std_tkweight_.at(i).size();j++){
      vtx_std_tkweight_helper_[i][j]=vtx_std_tkweight_.at(i).at(j);
    }
  }
}
