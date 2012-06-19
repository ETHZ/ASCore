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
// $Id: JetFillerPat.h,v 1.3.2.2 2012/04/04 12:21:32 fronga Exp $
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
  JetFillerPat( const edm::ParameterSet&, const bool& isRealData );
  virtual ~JetFillerPat(void) {}

  /// Fill all branches
  virtual void fillProducts(edm::Event&, const edm::EventSetup& );

private:


  //- Configuration parameters
  edm::InputTag fTag; 

  // Pre-selection
  double fMinpt;
  double fMaxeta;
	
};


#endif
