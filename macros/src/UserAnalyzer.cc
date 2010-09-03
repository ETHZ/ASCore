#include "UserAnalyzer.hh"
#include "UserAnalysis.hh"

#include "base/TreeAnalyzerBase.hh"
#include "base/TreeReader.hh"

using namespace std;

UserAnalyzer::UserAnalyzer(TTree *tree) : TreeAnalyzerBase(tree) {
	fUserAnalysis = new UserAnalysis(fTR);
}

UserAnalyzer::~UserAnalyzer(){
	delete fUserAnalysis;
	if(!fTR->fChain) cout << "UserAnalyzer ==> No chain!" << endl;
}

// Method for looping over the tree
void UserAnalyzer::Loop(){
	Long64_t nentries = fTR->GetEntries();
	cout << " total events in ntuples = " << fTR->GetEntries() << endl;
	
	// loop over all ntuple entries
	// nentries = 200;
	for( Long64_t jentry = 0; jentry < nentries; jentry++ ){
		PrintProgress(jentry);
		fTR->GetEntry(jentry);
                if ( fCurRun != fTR->Run ) {
                  fCurRun = fTR->Run;
                  fUserAnalysis->BeginRun(fCurRun);
                 }

		fUserAnalysis->Analyze();
	}
}

// Method called before starting the event loop
void UserAnalyzer::BeginJob(){
	fUserAnalysis->SetOutputDir(fOutputDir);
	fUserAnalysis->fVerbose = fVerbose;

	fUserAnalysis->Begin();

}

// Method called after finishing the event loop
void UserAnalyzer::EndJob(){
	fUserAnalysis->End();
}
