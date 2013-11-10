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

  vtx_std_tkind_helper_.resize(vtx_std_tkind_->size(),NULL);
  vtx_std_tkweight_helper_.resize(vtx_std_tkweight_->size(),NULL);
  
  for (size_t i=0; i<vtx_std_tkind_->size();i++){
    vtx_std_tkind_helper_.at(i) = new unsigned short[vtx_std_tkind_->at(i).size()];
    for (size_t j=0; j<vtx_std_tkind_->at(i).size();j++){
      (vtx_std_tkind_helper_.at(i))[j] = vtx_std_tkind_->at(i).at(j);
    }
  }
  
  for (size_t i=0; i<vtx_std_tkweight_->size();i++){
    vtx_std_tkweight_helper_.at(i) = new float[vtx_std_tkweight_->at(i).size()];
    for (size_t j=0; j<vtx_std_tkweight_->at(i).size();j++){
      (vtx_std_tkweight_helper_.at(i))[j]=vtx_std_tkweight_->at(i).at(j);
    }
  }

}

ETHVertexInfo::~ETHVertexInfo(){

  for (size_t i=0; i<vtx_std_tkind_helper_.size();i++){
    if (vtx_std_tkind_helper_.at(i)) delete[] vtx_std_tkind_helper_.at(i);
  }

  for (size_t i=0; i<vtx_std_tkweight_helper_.size();i++){
    if (vtx_std_tkweight_helper_.at(i)) delete[] vtx_std_tkweight_helper_.at(i);
  }

}
