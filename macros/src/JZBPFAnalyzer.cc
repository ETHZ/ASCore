#include "JZBPFAnalyzer.hh"
#include "JZBAnalysis.hh"

#include "base/TreeAnalyzerBase.hh"
#include "base/TreeReader.hh"

using namespace std;

JZBPFAnalyzer::JZBPFAnalyzer(TTree *tree, std::string dataType, bool fullCleaning) 
  : TreeAnalyzerBase(tree) {
  fJZBPFAnalysis = new JZBPFAnalysis(fTR,dataType,fullCleaning);
}

JZBPFAnalyzer::~JZBPFAnalyzer(){
	delete fJZBPFAnalysis;
	if(!fTR->fChain) cout << "JZBPFAnalyzer ==> No chain!" << endl;
}

// Method for looping over the tree
void JZBPFAnalyzer::Loop(){
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
      fJZBPFAnalysis->BeginRun(fCurRun);
      skipRun = false;
      if ( !CheckRun() ) skipRun = true;
    }
    // Check if new lumi is in JSON file
    if ( !skipRun && fCurLumi != fTR->LumiSection ) {
      fCurLumi = fTR->LumiSection;
      skipLumi = false; // Re-initialise
      if ( !CheckRunLumi() ) skipLumi = true;
    }
    if ( !(skipRun || skipLumi) ) fJZBPFAnalysis->Analyze();
  }
}

// Method called before starting the event loop
void JZBPFAnalyzer::BeginJob(string fdata_PileUp, string fmc_PileUp){
	fJZBPFAnalysis->outputFileName_ = outputFileName_;
	fJZBPFAnalysis->fVerbose = fVerbose;
	fJZBPFAnalysis->SetPileUpSrc(fdata_PileUp, fmc_PileUp);
	fJZBPFAnalysis->Begin();

}

// Method called after finishing the event loop
void JZBPFAnalyzer::EndJob(){
	fJZBPFAnalysis->End();
}
