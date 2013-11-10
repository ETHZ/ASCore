#ifndef __DiLeptonAnalysis_NTupleProducer_ETHVertexInfo_H__
#define __DiLeptonAnalysis_NTupleProducer_ETHVertexInfo_H__


#include "TVector3.h"
#include <vector>
#include "h2gglobe/VertexAnalysis/interface/HggVertexAnalyzer.h"

class ETHVertexInfo : public VertexInfoAdapter
{
public:

  ETHVertexInfo(int nvtx, std::vector<TVector3> * vtxes,
		int ntracks, std::vector<TVector3> * tracks_p3,
		std::vector<float> * tkPtErr, std::vector<int> * tkVtxId, 
		std::vector<float> * tkd0, std::vector<float> * tkd0Err, std::vector<float> * tkdz, std::vector<float> * tkdzErr,
		std::vector<bool> * tkIsHighPurity, std::vector<std::vector<unsigned short> > * vtx_std_tkind,
		std::vector<std::vector<float> > * vtx_std_tkweight, std::vector<int> * vtx_std_ntks
		);

  virtual int nvtx() const    { return nvtx_; };
  virtual int ntracks() const { return ntracks_; };

  virtual bool hasVtxTracks()  const { return true; };
  virtual const unsigned short * vtxTracks(int ii) const { return vtx_std_tkind_helper_.at(ii); };
  virtual int vtxNTracks(int ii) const { return vtx_std_ntks_->at(ii); };
  virtual const float * vtxTkWeights(int ii) const { return vtx_std_tkweight_helper_.at(ii); };

  virtual float tkpx(int ii) const { return tracks_p3_->at(ii).X(); };
  virtual float tkpy(int ii) const { return tracks_p3_->at(ii).Y(); };
  virtual float tkpz(int ii) const { return tracks_p3_->at(ii).Z(); };
	
  virtual float tkPtErr(int ii) const { return tkPtErr_->at(ii); };
  virtual int   tkVtxId(int ii) const { return -1; };

  virtual float tkWeight(int ii, int jj) const { return (vtx_std_tkweight_helper_.at(jj))[ii]; };
	
  virtual float vtxx(int ii) const { return vtxes_->at(ii).X(); };
  virtual float vtxy(int ii) const { return vtxes_->at(ii).Y(); };
  virtual float vtxz(int ii) const { return vtxes_->at(ii).Z(); };

  virtual float tkd0(int ii, int jj) const { return 0; } // CHECK
  virtual float tkd0Err(int ii, int jj) const { return 1; } // CHECK

  virtual float tkdz(int ii, int jj) const { return 0; } // CHECK
  virtual float tkdzErr(int ii, int jj) const { return 1; } // CHECK

  virtual bool tkIsHighPurity(int ii) const { return tkIsHighPurity_->at(ii); };

  ~ETHVertexInfo();
       	
private:

  int nvtx_;
  std::vector<TVector3> * vtxes_;
  int ntracks_;
  std::vector<TVector3> * tracks_p3_;
  std::vector<float> * tkPtErr_;
  std::vector<int> * tkVtxId_;
  std::vector<float> * tkd0_;
  std::vector<float> * tkd0Err_;
  std::vector<float> * tkdz_;
  std::vector<float> * tkdzErr_;
  std::vector<bool> * tkIsHighPurity_;
  std::vector<std::vector<unsigned short> > * vtx_std_tkind_;
  std::vector<std::vector<float> > * vtx_std_tkweight_;
  std::vector<int> * vtx_std_ntks_;

  std::vector<unsigned short *> vtx_std_tkind_helper_;
  std::vector<float *> vtx_std_tkweight_helper_;


};

#endif
