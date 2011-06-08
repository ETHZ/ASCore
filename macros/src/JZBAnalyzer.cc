#include "JZBAnalyzer.hh"
#include "JZBAnalysis.hh"

#include "base/TreeAnalyzerBase.hh"
#include "base/TreeReader.hh"

using namespace std;

JZBAnalyzer::JZBAnalyzer(TTree *tree, std::string dataType, bool fullCleaning) 
  : TreeAnalyzerBase(tree) {
  fJZBAnalysis = new JZBAnalysis(fTR,dataType,fullCleaning);
  fJZBPFAnalysis = new JZBPFAnalysis(fTR,dataType,fullCleaning);
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

TFile *fHistFile;

// Method called before starting the event loop
void JZBAnalyzer::BeginJob(string fdata_PileUp, string fmc_PileUp){
	fHistFile = new TFile(outputFileName_.c_str(), "RECREATE");
//	fJZBAnalysis->outputFileName_ = outputFileName_;
//	fJZBPFAnalysis->outputFileName_ = outputFileName_;
	fJZBAnalysis->fVerbose = fVerbose;
	fJZBPFAnalysis->fVerbose = fVerbose;
	fJZBAnalysis->SetPileUpSrc(fdata_PileUp, fmc_PileUp);
	fJZBPFAnalysis->SetPileUpSrc(fdata_PileUp, fmc_PileUp);
cout << "Starting the regular analysis " << endl;
	fJZBAnalysis->Begin(fHistFile);
cout << "Starting the new analysis " << endl;
	fJZBPFAnalysis->Begin(fHistFile);
cout << "Done with beginning ... " << endl;
}

// Method called after finishing the event loop
void JZBAnalyzer::EndJob(){
	cout << "Ok gotten here ... " << endl;
	fJZBAnalysis->End(fHistFile);
	cout << "Survived the first gracefully ... " << endl;
	fJZBPFAnalysis->End(fHistFile);
	cout << "Survived the second gracefully ... " << endl;

fHistFile->Close();
}
