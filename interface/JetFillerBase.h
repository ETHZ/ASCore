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
// $Id: JetFillerBase.h,v 1.10.2.2 2012/04/04 12:21:31 fronga Exp $
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
  JetFillerBase( const edm::ParameterSet& cfg, const bool& isRealData );
  virtual ~JetFillerBase(void) {}

  /// Define all branches
  virtual const std::vector<filler::PPair> declareProducts(void);
  /// Reset all branch containers
  virtual void resetProducts(void);
  /// Fill variables (needs to be implemented in specialised classes)
  virtual void fillProducts(edm::Event&, const edm::EventSetup& ) = 0;
  /// Put products in the event data
  virtual void putProducts( edm::Event& );



protected:

  /// Set and get jet type
  void setJetType( const JetType& type ) { fJetType = type; }
  const JetType jetType(void) const { return fJetType; }

  static const unsigned int gMaxNBtags      = 10;

  // To order indices by pt
  typedef std::pair<unsigned int,double> OrderPair;
  struct IndexByPt {
    const bool operator()(const OrderPair& j1, const OrderPair& j2 ) const {
      return j1.second > j2.second;
    }
  };

  JetType fJetType;	

  size_t gMaxnobjs;

  //- Configuration parameters
  std::vector<std::string> fBtagNames;

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
  std::auto_ptr<std::vector<float> > fTScale;
  std::auto_ptr<std::vector<float> > fTL1FastJetScale;
  std::auto_ptr<std::vector<float> > fTArea;
  std::auto_ptr<std::vector<int> >   fTFlavour;
  std::auto_ptr<std::vector<int> >   fTIDLoose;
  // b-taggers
  std::auto_ptr<std::vector<float> >  fTJbTagProb[gMaxNBtags];
  // jet algo dependant info
  std::auto_ptr<std::vector<int> >   fTNConstituents;
  std::auto_ptr<std::vector<int> >   fTNAssoTracks;
  std::auto_ptr<std::vector<float> > fTChfrac;
  std::auto_ptr<std::vector<float> > fTEMfrac;
  std::auto_ptr<std::vector<float> > fTID_HPD; 
  std::auto_ptr<std::vector<float> > fTID_RBX;    
  std::auto_ptr<std::vector<float> > fTID_n90Hits;
  std::auto_ptr<std::vector<int> >   fTn90;
  std::auto_ptr<std::vector<float> > fTID_resEMF; 
  std::auto_ptr<std::vector<int> >   fTChMult;
  std::auto_ptr<std::vector<int> >   fTNeuMult;
  std::auto_ptr<std::vector<float> > fTChHadfrac;
  std::auto_ptr<std::vector<float> > fTNeuHadfrac;
  std::auto_ptr<std::vector<float> > fTChEmfrac;
  std::auto_ptr<std::vector<float> > fTNeuEmfrac;
  std::auto_ptr<std::vector<float> > fTChMufrac;
  std::auto_ptr<std::vector<float> > fTPhofrac;
  std::auto_ptr<std::vector<float> > fTHFHadfrac;
  std::auto_ptr<std::vector<float> > fTHFEMfrac;

};


#endif
