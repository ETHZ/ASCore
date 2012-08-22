#include "JZBAnalyzer.hh"
#include "JZBAnalysis.hh"

#include "base/TreeAnalyzerBase.hh"
#include "base/TreeReader.hh"

using namespace std;

JZBAnalyzer::JZBAnalyzer(std::vector<std::string>& fileList, std::string dataType, bool fullCleaning, bool isModelScan, bool makeSmall, bool doGenInfo)  
   : TreeAnalyzerBase(fileList) {
  f_isModelScan=isModelScan;
  f_doGenInfo=doGenInfo;
  fJZBAnalysis = new JZBAnalysis(fTR,dataType,fullCleaning,isModelScan,makeSmall,doGenInfo);
}

JZBAnalyzer::~JZBAnalyzer() {
  delete fJZBAnalysis;
}

// Method for looping over the tree
void JZBAnalyzer::Loop(){
  Long64_t nentries = fTR->GetEntries();
  cout << " total events in ntuples = " << fTR->GetEntries() << endl;
  if ( fTR->GetEntries() < 1 ) {
    cerr << " No entry found: stopping here" << endl;
    return;
  }
  
  // loop over all ntuple entries
  if(fMaxEvents==-1)nentries=fTR->GetEntries();
  if(fMaxEvents>0)nentries=fMaxEvents;
  
  Long64_t jentry=0;
  for ( fTR->ToBegin();!(fTR->AtEnd()) && (jentry<fMaxEvents || fMaxEvents<0);++(*fTR) ) {
    PrintProgress(jentry++);
    if ( fCurRun != fTR->Run ) {
      fCurRun = fTR->Run;
      fJZBAnalysis->BeginRun(fCurRun);
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
    }
  }
}

TFile *fHistFile;

// Method called before starting the event loop
void JZBAnalyzer::BeginJob(string fdata_PileUp, string fmc_PileUp){
  fHistFile = new TFile(outputFileName_.c_str(), "RECREATE");
  //	Note: the next line is commented out because we are now not saving in the analysis routine anymore but at the "analyzer level"
  //	fJZBAnalysis->outputFileName_ = outputFileName_;
  fJZBAnalysis->SetVerbose(fVerbose);
  fJZBAnalysis->SetPileUpSrc(fdata_PileUp, fmc_PileUp);
  fJZBAnalysis->Begin(fHistFile);
}

// Method called after finishing the event loop
void JZBAnalyzer::EndJob(){
  fJZBAnalysis->End(fHistFile);
  fHistFile->Close();
}
