#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DiLeptonAnalysis/NTupleProducer/interface/FillerBase.h"


//________________________________________________________________________________________
FillerBase::FillerBase( const edm::ParameterSet& cfg, TTree* tree, const bool& isRealData )
  : fTree(tree),fIsRealData(isRealData){
	
  // Retrieve configuration parameters
  fPrefix = cfg.getUntrackedParameter<std::string>("prefix");

}


//________________________________________________________________________________________
void FillerBase::resetDouble(double* v, size_t size)
{
  for(size_t i = 0; i < size; ++i) v[i] = -999.9;
}

//________________________________________________________________________________________
void FillerBase::resetInt(int* v, size_t size)
{
  for(size_t i = 0; i < size; ++i) v[i] = -999;
}


//________________________________________________________________________________________
const bool FillerBase::addBranch( const char* name, const char* type, 
                                  void* address, const char* size )
{
  
  // Form input
  std::string fullname(fPrefix+name);
  
  std::string branchType(fullname);
  if ( size ) // Size needs to be pre-fixed
    branchType += "[" + std::string(fPrefix+size) + "]";
  branchType += "/"+std::string(type);

  // Declare branch
  TBranch* b = fTree->Branch(fullname.c_str(),address,branchType.c_str());

  return !(b==0); // return 1 if branch was successfully created, 0 otherwise

}
