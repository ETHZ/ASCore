#include "QuickAnalyzer.hh"
#include "QuickAnalysis.hh"

#include "base/TreeAnalyzerBase.hh"
#include "base/TreeReader.hh"

using namespace std;

QuickAnalyzer::QuickAnalyzer(std::vector<std::string>& fileList) : TreeAnalyzerBase(fileList) {
	fQuickAnalysis = new QuickAnalysis(fTR);
}

QuickAnalyzer::~QuickAnalyzer(){
	delete fQuickAnalysis;
}

// Method for looping over the tree
void QuickAnalyzer::Loop(){
	Long64_t nentries = fTR->GetEntries();
	cout << " total events in ntuples = " << fTR->GetEntries() << endl;
	
	// loop over all ntuple entries
	// nentries = 200;
	Long64_t jentry=0;
	for ( fTR->ToBegin();!(fTR->AtEnd()) && (jentry<fMaxEvents || fMaxEvents<0);++(*fTR) ) {
		PrintProgress(jentry++);
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
	fQuickAnalysis->SetVerbose(fVerbose);

	fQuickAnalysis->Begin(filename);

}

// Method called after finishing the event loop
void QuickAnalyzer::EndJob(){
	fQuickAnalysis->End();
}
