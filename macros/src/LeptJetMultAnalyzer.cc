#include "LeptJetMultAnalyzer.hh"

#include "base/TreeAnalyzerBase.hh"
#include "base/TreeReader.hh"
#include "MassAnalysis.hh"

using namespace std;

LeptJetMultAnalyzer::LeptJetMultAnalyzer(TTree *tree) : TreeAnalyzerBase(tree) {
	fMassAnalysis             = new MassAnalysis(fTR);
	Util::SetStyle();
}

LeptJetMultAnalyzer::~LeptJetMultAnalyzer(){
	delete fMassAnalysis;
	if(!fTR->fChain) cout << "LeptJetMultAnalyzer ==> No chain!" << endl;
}

// Method for looping over the tree
void LeptJetMultAnalyzer::Loop(){
	Long64_t nentries = fTR->GetEntries();
	cout << " total events in ntuples = " << fTR->GetEntries() << endl;
	
	if(fMaxEvents==-1)nentries=fTR->GetEntries();
	if(fMaxEvents > 0){
		nentries=min((Long64_t)fMaxEvents, fTR->GetEntries());
		cout << " only running on first " << nentries << " events" << endl;
	}
	// loop over all ntuple entries
	for( Long64_t jentry = 0; jentry < nentries; jentry++ ){
		PrintProgress(jentry);
		fTR->GetEntry(jentry);
                if ( fCurRun != fTR->Run ) {
        		fCurRun = fTR->Run;
			fMassAnalysis        ->BeginRun(fCurRun);
                  	skipRun = false; // re-initialize
                  	if ( !CheckRun() ) skipRun = true;
                }
                // Check if new lumi is in JSON file
                if ( !skipRun && fCurLumi != fTR->LumiSection ) {
                  	fCurLumi = fTR->LumiSection;
                  	skipLumi = false; // Re-initialise
                  	if ( !CheckRunLumi() ) skipLumi = true;
                }
		if ( !(skipRun || skipLumi) ) {
			fMassAnalysis        ->Analyze();
		}
	}
}

// Method called before starting the event loop
void LeptJetMultAnalyzer::BeginJob(TString filename, TString setofcuts, float lumi, bool isData, string data_PileUp, string mc_PileUp){
	fMassAnalysis             ->ReadCuts(setofcuts);
	fMassAnalysis             ->SetType(isData);
	fMassAnalysis             ->SetPileUpSrc(data_PileUp, mc_PileUp);
	fMassAnalysis             ->SetOutputDir(fOutputDir);
	fMassAnalysis             ->fVerbose        = fVerbose;
	fMassAnalysis             ->Begin(filename);

        fMassAnalysis             ->isS3        = isS3;
	fMassAnalysis             ->noPU        = noPU;
	//if(is)
}

// Method called after finishing the event loop
void LeptJetMultAnalyzer::EndJob(){
  fMassAnalysis         ->End();
  cout << " LeptJetMultAnalyzer::End()                                             " << endl;
  
}
