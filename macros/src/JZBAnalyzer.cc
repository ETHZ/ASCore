#include "JZBAnalyzer.hh"
#include "JZBAnalysis.hh"

#include "base/TreeAnalyzerBase.hh"
#include "base/TreeReader.hh"

using namespace std;

JZBAnalyzer::JZBAnalyzer(TTree *tree) : TreeAnalyzerBase(tree) {
	fJZBAnalysis = new JZBAnalysis(fTR);
}

JZBAnalyzer::~JZBAnalyzer(){
	delete fJZBAnalysis;
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
                 }
		fJZBAnalysis->Analyze();
	}
}

// Method called before starting the event loop
void JZBAnalyzer::BeginJob(){
	//fJZBAnalysis->SetOutputFile(fOutputFile);
	fJZBAnalysis->outputFileName_ = outputFileName_;
	fJZBAnalysis->fVerbose = fVerbose;

	fJZBAnalysis->Begin();

}

// Method called after finishing the event loop
void JZBAnalyzer::EndJob(){
	fJZBAnalysis->End();
}
