#include "TTWAcceptanceAnalyzer.hh"
#include "TTWAcceptanceAnalysis.hh"

#include "base/TreeAnalyzerBase.hh"
#include "base/TreeReader.hh"

using namespace std;

TTWAcceptanceAnalyzer::TTWAcceptanceAnalyzer(TTree *tree) : TreeAnalyzerBase(tree) {
	fTTWAcceptanceAnalysis = new TTWAcceptanceAnalysis(fTR);
}

TTWAcceptanceAnalyzer::~TTWAcceptanceAnalyzer(){
	delete fTTWAcceptanceAnalysis;
	if(!fTR->fChain) cout << "TTWAcceptanceAnalyzer ==> No chain!" << endl;
}

// Method for looping over the tree
void TTWAcceptanceAnalyzer::Loop(){
	Long64_t nentries = fTR->GetEntries();
	cout << " total events in ntuples = " << fTR->GetEntries() << endl;
	
	// loop over all ntuple entries
	// nentries = 200;
	for( Long64_t jentry = 0; jentry < nentries; jentry++ ){
		PrintProgress(jentry);
		fTR->GetEntry(jentry);
		if ( fCurRun != fTR->Run ) {
			fCurRun = fTR->Run;
			fTTWAcceptanceAnalysis->BeginRun(fCurRun);
		}
		
		fTTWAcceptanceAnalysis->Analyze();
	}
}

// Method called before starting the event loop
void TTWAcceptanceAnalyzer::BeginJob(){
	fTTWAcceptanceAnalysis->SetOutputDir(fOutputDir);
	fTTWAcceptanceAnalysis->fVerbose = fVerbose;

	fTTWAcceptanceAnalysis->Begin();

}

// Method called after finishing the event loop
void TTWAcceptanceAnalyzer::EndJob(){
	fTTWAcceptanceAnalysis->End();
}
