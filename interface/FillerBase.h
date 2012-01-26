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
// $Id: FillerBase.h,v 1.2 2012/01/06 09:08:00 pnef Exp $
//
//

#include <string>
#include <vector>
#include <typeinfo>

#include "TTree.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/TypeID.h"

namespace filler {
  // Gory details of production (need to hand over to producer class)
  // This is because of how the EDFilter::produce() method works.
  // So we make a list of products and type names, 
  // and use a specialized EDFilter::produce() method.
  typedef std::pair<edm::TypeID,std::string> PPair;
}


class FillerBase {

public:
  /// Constructor: set pointer to tree
  FillerBase( const edm::ParameterSet& cfg, const bool& isRealData );
  virtual ~FillerBase(void) {}

  // Interface: needs to be implemented in specialised classes
  /// Define all variables
  virtual const std::vector<filler::PPair> declareProducts(void) = 0;
  /// Reset all variable containers
  virtual void resetProducts(void) = 0;
  /// Fill variables (called for each event)
  virtual void putProducts(edm::Event&, const edm::EventSetup&) = 0;

protected:

  /// Add a product to the list
  const bool addProduct(const char* name, const type_info& address);

  std::string fPrefix;        /// Prefix for branches
  bool   fIsRealData;         /// Global switch
  std::vector<filler::PPair> productList;
  

};


#endif
