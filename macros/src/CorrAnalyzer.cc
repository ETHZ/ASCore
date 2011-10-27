#include "CorrAnalyzer.hh"
#include "ZeeAnalysis.hh"
#include "HggAnalysis.hh"

#include "base/TreeAnalyzerBase.hh"
#include "base/TreeReader.hh"

using namespace std;

CorrAnalyzer::CorrAnalyzer(TTree *tree) : TreeAnalyzerBase(tree) {
	fZeeAnalysis = new ZeeAnalysis(fTR);
	fHggAnalysis = new HggAnalysis(fTR);
}

CorrAnalyzer::~CorrAnalyzer(){
	delete fZeeAnalysis;
	delete fHggAnalysis;
	if(!fTR->fChain) cout << "CorrAnalyzer ==> No chain!" << endl;
}

// Method for looping over the tree
void CorrAnalyzer::Loop(){
	Long64_t nentries = fTR->GetEntries();
	cout << " total events in ntuples = " << fTR->GetEntries() << endl;
	
	// loop over all ntuple entries
	// nentries = 200;
	for( Long64_t jentry = 0; jentry < nentries; jentry++ ){
		PrintProgress(jentry);
		fTR->GetEntry(jentry);
                if ( fCurRun != fTR->Run ) {
                  fCurRun = fTR->Run;
                  fZeeAnalysis->BeginRun(fCurRun);
                  fHggAnalysis->BeginRun(fCurRun);
                 }

		fZeeAnalysis->Analyze();
		fHggAnalysis->Analyze();
	}
}

// Method called before starting the event loop
void CorrAnalyzer::BeginJob(){
	fZeeAnalysis->SetOutputDir(fOutputDir);
	fZeeAnalysis->fVerbose = fVerbose;

	fZeeAnalysis->Begin();

	fHggAnalysis->SetOutputDir(fOutputDir);
	fHggAnalysis->fVerbose = fVerbose;

	fHggAnalysis->Begin();

}

// Method called after finishing the event loop
void CorrAnalyzer::EndJob(){
	fZeeAnalysis->End();
	fHggAnalysis->End();
}
