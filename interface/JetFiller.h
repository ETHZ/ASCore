#ifndef __DiLeptonAnalysis_NTupleProducer_JetFiller_H__
#define __DiLeptonAnalysis_NTupleProducer_JetFiller_H__
// 
// Package: NTupleProducer
// Class:   JetFiller
//
/* class JetFiller
   JetFiller.h
   Description:  generic class for basic jet dumper

*/
//
// $Id: JetFiller.h,v 1.1 2010/05/17 09:54:20 fronga Exp $
//
//

#include <string>
#include <vector>

#include "TTree.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/VertexReco/interface/Vertex.h"


class JetFiller {
public:
  /// Constructor: set pointer to tree
  JetFiller( const edm::ParameterSet&, TTree* tree, 
                 const bool& isPat, const bool& isRealData );
  virtual ~JetFiller(void);

  /// Define all branches
  virtual void createBranches(void);
  /// Fill all branches
  virtual const int fillBranches(const edm::Event&, const edm::EventSetup& );
  /// Reset all branch containers
  virtual void reset(void);

protected:

  /// Add a branch to the tree (takes care of prefixing branch names)
  const bool addBranch(const char* name, const char* type, 
                       void* address, const char* size = 0);
  void resetDouble(double* v, size_t size = 1);
  void resetInt(int* v, size_t size = 1);

  // To order indices by pt
  typedef std::pair<unsigned int,double> OrderPair;
  struct IndexByPt {
    const bool operator()(const OrderPair& j1, const OrderPair& j2 ) const {
      return j1.second > j2.second;
    }
  };

  
  std::string fPrefix;        /// Prefix for branches
  TTree* fTree;               /// Pointer to tree to fill
  bool   fIsPat, fIsRealData; /// Global switches

  size_t gMaxnobjs;

  //- Configuration parameters
  edm::InputTag fTag;
  std::string fJetCorrs; 

  // Pre-selection
  double fMinpt;
  double fMaxeta;

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
  int*    fTNConstituents;
  // Calo jets
  double* fTEMfrac;
  double* fTHadFrac;
  // PF jets
  double* fTChHadFrac;
  double* fTNeuHadFrac;
  double* fTChEmFrac;
  double* fTNeuEmFrac;
  int*    fTMuonMultiplicity;

};


#endif
