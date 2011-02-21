#ifndef __DiLeptonAnalysis_NTupleProducer_JetFillerPat_H__
#define __DiLeptonAnalysis_NTupleProducer_JetFillerPat_H__
// 
// Package: NTupleProducer
// Class:   JetFillerPat
//
/* class JetFillerPat
   JetFillerPat.h
   Description:  generic class for basic jet dumper

*/
//
// $Id: JetFillerPat.h,v 1.2 2010/11/13 14:31:06 pnef Exp $
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


class JetFillerPat : public JetFillerBase {
public:
  /// Constructor: set pointer to tree
  JetFillerPat( const edm::ParameterSet&, TTree* tree, 
                 const bool& isPat, const bool& isRealData );
  virtual ~JetFillerPat(void) {}

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
  double fBtagMatchdeltaR;
	
};


#endif
