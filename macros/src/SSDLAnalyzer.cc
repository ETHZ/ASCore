#include "SSDLAnalyzer.hh"
#include "base/TreeReader.hh"

using namespace std;

SSDLAnalyzer::SSDLAnalyzer(TTree *tree) : TreeAnalyzerBase(tree) {
	fSSDLAnalysis = new SSDLAnalysis(fTR);
}

SSDLAnalyzer::~SSDLAnalyzer(){
	delete fSSDLAnalysis;
	if(!fTR->fChain) cout << "SSDLAnalyzer ==> No chain!" << endl;
}

// Method for looping over the tree
void SSDLAnalyzer::Loop(Long64_t maxEvents, Int_t prescale){
	Long64_t nentries = fTR->GetEntries();
	cout << " total events in ntuples = " << nentries << endl;
	if ( maxEvents>=0 ) { 
		nentries = maxEvents;
		cout << " processing only " << nentries << " events out of it" << endl;
	}
	if ( prescale>1 ) cout << " processing only every " << prescale << " events" << endl;
	
	//nentries = 1000;
	for( Long64_t jentry = 0; jentry < nentries; jentry++ ){
		PrintProgress(jentry);
		if ( prescale>1 && jentry%prescale ) continue; // Prescale processing...
		fTR->GetEntry(jentry);

		if ( fCurRun != fTR->Run ) {
			fCurRun = fTR->Run;
			fSSDLAnalysis->BeginRun(fCurRun);
		}
		fSSDLAnalysis->Analyze();
	}
}

// Method called before starting the event loop
void SSDLAnalyzer::BeginJob(){
	fSSDLAnalysis->SetOutputDir(fOutputDir);
	fSSDLAnalysis->SetOutputFile("SSDLTree");
	fSSDLAnalysis->SetVerbose(fVerbose);
	fSSDLAnalysis->Begin();
}

// Method called after finishing the event loop
void SSDLAnalyzer::EndJob(){
	fSSDLAnalysis->End();
}
