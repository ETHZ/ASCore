#include "DiLeptonAnalysis/NTupleProducer/interface/FillerBase.h"


//________________________________________________________________________________________
FillerBase::FillerBase( const edm::ParameterSet& cfg, const bool& isRealData )
  : fIsRealData(isRealData){
	
  // Retrieve configuration parameters
  fPrefix = cfg.getParameter<std::string>("prefix");

}


//________________________________________________________________________________________
void FillerBase::addProduct( const char* name, const type_info& type)
{
  
  // Form input
  edm::TypeID tid(type);

  typeList.push_back( std::make_pair(tid,fullName(name)) );

}
