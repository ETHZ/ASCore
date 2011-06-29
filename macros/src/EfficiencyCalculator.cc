#include "EfficiencyCalculator.hh"
#include "RunEfficiency.hh"

#include "base/TreeAnalyzerBase.hh"
#include "base/TreeReader.hh"

using namespace std;

EfficiencyCalculator::EfficiencyCalculator(TTree *tree, std::string dataType, bool fullCleaning) 
  : TreeAnalyzerBase(tree) {
  runEfficiency = new RunEfficiency(fTR,dataType,fullCleaning);
}

EfficiencyCalculator::~EfficiencyCalculator(){
	delete runEfficiency;
	if(!fTR->fChain) cout << "EfficiencyCalculator ==> No chain!" << endl;
}

// Method for looping over the tree
void EfficiencyCalculator::Loop(){
  Long64_t nentries = fTR->GetEntries();
  cout << " total events in ntuples = " << fTR->GetEntries() << endl;
  
  // loop over all ntuple entries
  if(fMaxEvents==-1)nentries=fTR->GetEntries();
  if(fMaxEvents>0)nentries=fMaxEvents;
  
  for( Long64_t jentry = 0; jentry < nentries; jentry++ ){
    PrintProgress(jentry);
    fTR->GetEntry(jentry);
    if ( fCurRun != fTR->Run ) {
      fCurRun = fTR->Run;
      runEfficiency->BeginRun(fCurRun);
      skipRun = false;
      if ( !CheckRun() ) skipRun = true;
    }
    // Check if new lumi is in JSON file
    if ( !skipRun && fCurLumi != fTR->LumiSection ) {
      fCurLumi = fTR->LumiSection;
      skipLumi = false; // Re-initialise
      if ( !CheckRunLumi() ) skipLumi = true;
    }
    if ( !(skipRun || skipLumi) ) runEfficiency->Analyze();
  }
}

// Method called before starting the event loop
void EfficiencyCalculator::BeginJob(string fdata_PileUp, string fmc_PileUp){
	runEfficiency->outputFileName_ = outputFileName_;
	runEfficiency->fVerbose = fVerbose;
	runEfficiency->SetPileUpSrc(fdata_PileUp, fmc_PileUp);
	runEfficiency->Begin();

}

// Method called after finishing the event loop
void EfficiencyCalculator::EndJob(){
	runEfficiency->End();
}
