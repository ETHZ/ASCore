#ifndef __DiLeptonAnalysis_NTupleProducer_FillerBase_H__
#define __DiLeptonAnalysis_NTupleProducer_FillerBase_H__
// 
// Package: NTupleProducer
// Class:   FillerBase
//
/* class FillerBase
   FillerBase.h
   Description:  generic class for basic information dumper

*/
//
// $Id: FillerBase.h,v 1.4 2011/02/23 19:34:29 stiegerb Exp $
//
//

#include <string>
#include <vector>

#include "TTree.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"


class FillerBase {
public:
  /// Constructor: set pointer to tree
  FillerBase( const edm::ParameterSet& cfg, TTree* tree, const bool& isRealData );
  virtual ~FillerBase(void) {}

  // Interface: needs to be implemented in specialised classes
  /// Define all branches
  virtual void createBranches(void) = 0;
  /// Reset all branch containers
  virtual void reset(void) = 0;
  /// Fill all branches
  virtual const int fillBranches(const edm::Event&, const edm::EventSetup& ) = 0;


protected:

  /// Add a branch to the tree (takes care of prefixing branch names)
  const bool addBranch(const char* name, const char* type, 
                       void* address, const char* size = 0);
  /// Resetting
  void resetDouble(double* v, size_t size = 1);
  void resetInt(int* v, size_t size = 1);

  std::string fPrefix;        /// Prefix for branches
  TTree* fTree;               /// Pointer to tree to fill
  bool   fIsRealData;         /// Global switch

};


#endif
