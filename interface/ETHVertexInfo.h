#ifndef __DiLeptonAnalysis_NTupleProducer_ETHVertexInfo_H__
#define __DiLeptonAnalysis_NTupleProducer_ETHVertexInfo_H__


#include "h2gglobe/VertexAnalysis/interface/HggVertexAnalyzer.h"

class ETHVertexInfo : public VertexInfoAdapter
{
public:

  ETHVertexInfo(int nvtx, 
                std::vector<float> vtxx, std::vector<float> vtxy, std::vector<float> vtxz, 
		int ntracks, float * tkpx, float * tkpy, float * tkpz,
		float * tkPtErr, int * tkVtxId, 
		float * tkd0, float * tkd0Err, float * tkdz, float * tkdzErr,
		bool * tkIsHighPurity, std::vector<std::vector<unsigned short> > vtx_std_tkind, std::vector<std::vector<float> > vtx_std_tkweight, int * vtx_std_ntks
		);

  virtual int nvtx() const    { return nvtx_; };
  virtual int ntracks() const { return ntracks_; };

  virtual bool hasVtxTracks()  const { return true; };
  virtual const unsigned short * vtxTracks(int ii) const { return vtx_std_tkind_helper_[ii]+0; };
  virtual int vtxNTracks(int ii) const { return vtx_std_ntks_[ii]; };
  virtual const float * vtxTkWeights(int ii) const { return vtx_std_tkweight_helper_[ii]+0; };

  virtual float tkpx(int ii) const { return tkpx_ != 0 ? tkpx_[ii] : 0.; };
  virtual float tkpy(int ii) const { return tkpx_ != 0 ? tkpy_[ii] : 0.; };
  virtual float tkpz(int ii) const { return tkpx_ != 0 ? tkpz_[ii] : 0.; };
	
  virtual float tkPtErr(int ii) const { return tkPtErr_  != 0 ? tkPtErr_[ii] : 999.; };
  virtual int   tkVtxId(int ii) const { return tkVtxId_  != 0 ? tkVtxId_[ii] : 999; };

  //	virtual float tkWeight(int ii, int jj) const { return tkWeight_ != 0 ? tkWeight_[ii]*(float)( tkVtxId(ii) == jj) : 0.; };
  virtual float tkWeight(int ii, int jj) const { 
    if (jj>=(int)(vtx_std_tkind_.size())) {std::cout << "wrong vertex index call" << std::endl; return 0;}

    int trkid=-1;
    int n=0;
    for (std::vector<unsigned short>::const_iterator it=vtx_std_tkind_.at(jj).begin(); it!=vtx_std_tkind_.at(jj).end(); it++){
      if (*it==ii) trkid=n;
      n++;
    }

    float out=0;
    if (trkid!=-1) out=vtx_std_tkweight_.at(jj).at(trkid);
    return out;

  };
  
	
  virtual float vtxx(int ii) const { return vtxx_ != 0 ? vtxx_[ii] : 0.; };
  virtual float vtxy(int ii) const { return vtxy_ != 0 ? vtxy_[ii] : 0.; };
  virtual float vtxz(int ii) const { return vtxz_ != 0 ? vtxz_[ii] : 0.; };

  virtual float tkd0(int ii, int jj) const { assert(tkVtxId(ii) == jj); return tkd0_ != 0 ? tkd0_[ii] : 0.; };
  virtual float tkd0Err(int ii, int jj) const { assert(tkVtxId(ii) == jj); return tkd0Err_ != 0 ? tkd0Err_[ii] : 0.; };

  virtual float tkdz(int ii, int jj) const { assert(tkVtxId(ii) == jj); return tkdz_ != 0 ? tkdz_[ii] : 0.; };
  virtual float tkdzErr(int ii, int jj) const { assert(tkVtxId(ii) == jj); return tkdzErr_ != 0 ? tkdzErr_[ii] : 0.; };

  virtual bool tkIsHighPurity(int ii) const { return tkIsHighPurity_ != 0 ? tkIsHighPurity_[ii] : 0.; };

  virtual ~ETHVertexInfo() {}
       	
private:

  int nvtx_;
  float * vtxx_;
  float * vtxy_;
  float * vtxz_;

  static const unsigned int __TRK_AUX_ARRAYS_DIM__ = 2000;
  static const unsigned int __VTX_AUX_ARRAYS_DIM__ = 100;

  float vtxx_array[__VTX_AUX_ARRAYS_DIM__];
  float vtxy_array[__VTX_AUX_ARRAYS_DIM__];
  float vtxz_array[__VTX_AUX_ARRAYS_DIM__];

  int ntracks_;
  float * tkpx_;
  float * tkpy_;
  float * tkpz_;
  float * tkPtErr_;
  int * tkVtxId_;	

  float * tkd0_;
  float * tkd0Err_;
  float * tkdz_;
  float * tkdzErr_;

  bool * tkIsHighPurity_;
  
  std::vector<std::vector<unsigned short> > vtx_std_tkind_;
  std::vector<std::vector<float> > vtx_std_tkweight_;
  int * vtx_std_ntks_;

  unsigned short vtx_std_tkind_helper_[__VTX_AUX_ARRAYS_DIM__][__TRK_AUX_ARRAYS_DIM__];
  float vtx_std_tkweight_helper_[__VTX_AUX_ARRAYS_DIM__][__TRK_AUX_ARRAYS_DIM__];


};

#endif
