#ifndef __DiLeptonAnalysis_NTupleProducer_JetFillerReco_H__
#define __DiLeptonAnalysis_NTupleProducer_JetFillerReco_H__
// 
// Package: NTupleProducer
// Class:   JetFillerReco
//
/* class JetFillerReco
   JetFillerReco.h
   Description:  generic class for basic jet dumper

*/
//
// $Id: $
//
//

#include <string>
#include <vector>

#include "TTree.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "DiLeptonAnalysis/NTupleProducer/interface/JetFillerBase.h"


class JetFillerReco : public JetFillerBase {
public:
  /// Constructor: set pointer to tree
  JetFillerReco( const edm::ParameterSet&, TTree* tree, 
                 const bool& isPat, const bool& isRealData );
  virtual ~JetFillerReco(void) {}

  /// Fill all branches
  virtual const int fillBranches(const edm::Event&, const edm::EventSetup& );

private:


  //- Configuration parameters
  edm::InputTag fTag; 
  edm::InputTag fJetID;	
  edm::InputTag fJetTracksTag;
  std::string fJetCorrs; 
  bool fIsRealData;         /// Global switch

  // Pre-selection
  double fMinpt;
  double fMaxeta;
	
};


#endif
