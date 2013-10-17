#include "DiLeptonAnalysis/NTupleProducer/interface/ETHVertexInfo.h"

ETHVertexInfo::ETHVertexInfo(int nvtx, std::vector<TVector3> * vtxes,
			     int ntracks, std::vector<TVector3> * tracks_p3,
			     std::vector<float> * tkPtErr, std::vector<int> * tkVtxId, 
			     std::vector<float> * tkd0, std::vector<float> * tkd0Err, std::vector<float> * tkdz, std::vector<float> * tkdzErr,
			     std::vector<bool> * tkIsHighPurity, std::vector<std::vector<unsigned short> > * vtx_std_tkind,
			     std::vector<std::vector<float> > * vtx_std_tkweight, std::vector<int> * vtx_std_ntks
			     ) :
  nvtx_(nvtx),
  vtxes_(vtxes),
  ntracks_(ntracks),
  tracks_p3_(tracks_p3),
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
  assert (vtxes_->size()<__VTX_AUX_ARRAYS_DIM__);
  assert (tracks_p3_->size()<__TRK_AUX_ARRAYS_DIM__);
  
  for (unsigned int i=0; i<vtx_std_tkind_->size();i++){
    for (unsigned int j=0; j<vtx_std_tkind_->at(i).size();j++){
      vtx_std_tkind_helper_[i][j]=vtx_std_tkind_->at(i).at(j);
    }
  }
  
  for (unsigned int i=0; i<vtx_std_tkweight_->size();i++){
    for (unsigned int j=0; j<vtx_std_tkweight_->at(i).size();j++){
      vtx_std_tkweight_helper_[i][j]=vtx_std_tkweight_->at(i).at(j);
    }
  }
}
