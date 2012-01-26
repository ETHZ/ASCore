#include "DiLeptonAnalysis/NTupleProducer/interface/FillerBase.h"


//________________________________________________________________________________________
FillerBase::FillerBase( const edm::ParameterSet& cfg, const bool& isRealData )
  : fIsRealData(isRealData){
	
  // Retrieve configuration parameters
  fPrefix = cfg.getUntrackedParameter<std::string>("prefix");

}


//________________________________________________________________________________________
const bool FillerBase::addProduct( const char* name,const type_info& type )
{
  
  // Form input
  std::string fullname(fPrefix+name);
  edm::TypeID tid(type);

  //std::cout << "Adding product " << fullname << " of type " << tid  << std::endl;
  productList.push_back( std::make_pair(tid,fullname) );

}
