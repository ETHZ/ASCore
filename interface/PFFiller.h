#ifndef __DiLeptonAnalysis_NTupleProducer_PFFiller_H__
#define __DiLeptonAnalysis_NTupleProducer_PFFiller_H__
// 
// Package: NTupleProducer
// Class:   PFFiller
//
/* class PFFiller
   PFFiller.h
   Description:  generic class for basic jet dumper

*/
//
// $Id: PFFiller.h,v 1.3.2.2 2012/04/04 12:21:32 fronga Exp $
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


class PFFiller : public FillerBase {
public:
  /// Constructor: set pointer to tree
  PFFiller( const edm::ParameterSet&, const bool& isRealData );
  virtual ~PFFiller(void) {}

  /// Define all branches
  virtual const std::vector<filler::PPair> declareProducts(void);
  /// Reset all branch containers
  virtual void resetProducts(void);
  /// Fill variables (needs to be implemented in specialised classes)
  virtual void fillProducts(edm::Event&, const edm::EventSetup& );
  /// Put products in the event data
  virtual void putProducts( edm::Event& );

private:

  size_t gMaxnobjs;

  //- Configuration parameters
  edm::InputTag fTag; 

  // Pre-selection
  double fMinpt;
  double fMaxeta;

    // Stored variables
  std::auto_ptr<int>     fTNObjs;
  std::auto_ptr<std::vector<float> > fTPx;
  std::auto_ptr<std::vector<float> > fTPy;
  std::auto_ptr<std::vector<float> > fTPz;
  std::auto_ptr<std::vector<float> > fTPt;
  std::auto_ptr<std::vector<float> > fTE;
  std::auto_ptr<std::vector<float> > fTEt;
  std::auto_ptr<std::vector<float> > fTEta;
  std::auto_ptr<std::vector<float> > fTPhi;
  std::auto_ptr<std::vector<int> >   fTType;
  std::auto_ptr<std::vector<float> > fTVx;
  std::auto_ptr<std::vector<float> > fTVy;
  std::auto_ptr<std::vector<float> > fTVz;

	
};


#endif
