#include "JZBAnalyzer.hh"
#include "JZBAnalysis.hh"

#include "base/TreeAnalyzerBase.hh"
#include "base/TreeReader.hh"

using namespace std;

JZBAnalyzer::JZBAnalyzer(TTree *tree, std::string dataType, bool fullCleaning) 
  : TreeAnalyzerBase(tree) {
  fJZBAnalysis = new JZBAnalysis(fTR,dataType,fullCleaning);
}

JZBAnalyzer::~JZBAnalyzer(){
	delete fJZBAnalysis;
	delete fJZBPFAnalysis;
	if(!fTR->fChain) cout << "JZBAnalyzer ==> No chain!" << endl;
}

// Method for looping over the tree
void JZBAnalyzer::Loop(){
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
      fJZBAnalysis->BeginRun(fCurRun);
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
    if ( !(skipRun || skipLumi) ) {
	fJZBAnalysis->Analyze();
	fJZBPFAnalysis->Analyze();
    }
  }
}

// Method called before starting the event loop
void JZBAnalyzer::BeginJob(string fdata_PileUp, string fmc_PileUp){
	fJZBAnalysis->outputFileName_ = outputFileName_;
	fJZBPFAnalysis->outputFileName_ = outputFileName_;
	fJZBAnalysis->fVerbose = fVerbose;
	fJZBPFAnalysis->fVerbose = fVerbose;
	fJZBAnalysis->SetPileUpSrc(fdata_PileUp, fmc_PileUp);
	fJZBPFAnalysis->SetPileUpSrc(fdata_PileUp, fmc_PileUp);
	fJZBAnalysis->Begin();
	fJZBPFAnalysis->Begin();
}

// Method called after finishing the event loop
void JZBAnalyzer::EndJob(){
	fJZBAnalysis->End();
}
