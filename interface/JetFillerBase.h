#ifndef __DiLeptonAnalysis_NTupleProducer_JetFillerBase_H__
#define __DiLeptonAnalysis_NTupleProducer_JetFillerBase_H__
// 
// Package: NTupleProducer
// Class:   JetFillerBase
//
/* class JetFillerBase
   JetFillerBase.h
   Description:  generic class for basic jet information dumper

*/
//
// $Id: JetFillerBase.h,v 1.5 2011/03/29 15:49:11 pnef Exp $
//
//

#include <string>
#include <vector>

#include "TTree.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "DiLeptonAnalysis/NTupleProducer/interface/FillerBase.h"

class JetFillerBase : public FillerBase {
public:
  // Enumerate jet types
  enum JetType {
    CALO, PF, JPT, unknown
  };

  /// Constructor: set pointer to tree
  JetFillerBase( const edm::ParameterSet& cfg, TTree* tree, const bool& isRealData );
  virtual ~JetFillerBase(void);

  /// Define all branches
  void createBranches(void);
  /// Reset all branch containers
  void reset(void);

  /// Fill all branches (needs to be implemented in specialised classes)
  virtual const int fillBranches(const edm::Event&, const edm::EventSetup& ) = 0;


protected:

  /// Set and get jet type
  void setJetType( const JetType& type ) { fJetType = type; }
  const JetType jetType(void) const { return fJetType; }

  // To order indices by pt
  typedef std::pair<unsigned int,double> OrderPair;
  struct IndexByPt {
    const bool operator()(const OrderPair& j1, const OrderPair& j2 ) const {
      return j1.second > j2.second;
    }
  };

  JetType fJetType;	

  size_t gMaxnobjs;

  // Tree variables
  int     fTnobj;
  double* fTpx;
  double* fTpy;
  double* fTpz;
  double* fTpt;
  double* fTe;
  double* fTet;
  double* fTeta;
  double* fTphi;
  double* fTscale;
  double* fTarea;
  // b-taggers
  double* fTjbTagProbTkCntHighEff ;
  double* fTjbTagProbTkCntHighPur ;
  double* fTjbTagProbSimpSVHighEff;
  double* fTjbTagProbSimpSVHighPur;
  // jet algo dependant info
  int*    fTNConstituents;
  double* fTjChfrac;
  double* fTEMfrac;
  int*    fTNeuMult;
  int*    fTChMult;
  double* fTID_HPD; 
  double* fTID_RBX;    
  double* fTID_n90Hits;
  double* fTID_resEMF; 
  int*    fTjnAssoTracks;
  double* fTChHadFrac;
  double* fTNeuHadFrac;
  double* fTChEmFrac;
  double* fTNeuEmFrac;

};


#endif
