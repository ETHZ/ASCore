#include "QuickAnalyzer.hh"
#include "QuickAnalysis.hh"

#include "base/TreeAnalyzerBase.hh"
#include "base/TreeReader.hh"

using namespace std;

QuickAnalyzer::QuickAnalyzer(TTree *tree) : TreeAnalyzerBase(tree) {
	fQuickAnalysis = new QuickAnalysis(fTR);
}

QuickAnalyzer::~QuickAnalyzer(){
	delete fQuickAnalysis;
	if(!fTR->fChain) cout << "QuickAnalyzer ==> No chain!" << endl;
}

// Method for looping over the tree
void QuickAnalyzer::Loop(){
	Long64_t nentries = fTR->GetEntries();
	cout << " total events in ntuples = " << fTR->GetEntries() << endl;
	
	// loop over all ntuple entries
	// nentries = 200;
	for( Long64_t jentry = 0; jentry < nentries; jentry++ ){
		PrintProgress(jentry);
		fTR->GetEntry(jentry);
                if ( fCurRun != fTR->Run ) {
                  fCurRun = fTR->Run;
                  fQuickAnalysis->BeginRun(fCurRun);
                 }

		fQuickAnalysis->Analyze();
	}
}

// Method called before starting the event loop
void QuickAnalyzer::BeginJob(TString filename){
	fQuickAnalysis->SetOutputDir(fOutputDir);
	fQuickAnalysis->fVerbose = fVerbose;

	fQuickAnalysis->Begin(filename);

}

// Method called after finishing the event loop
void QuickAnalyzer::EndJob(){
	fQuickAnalysis->End();
}
