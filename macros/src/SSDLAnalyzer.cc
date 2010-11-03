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
void SSDLAnalyzer::Loop(Int_t prescale){
	Long64_t nentries = fTR->GetEntries();
	cout << " total events in ntuples = " << nentries << endl;

	if( fMaxEvents > -1 ){
		cout << " will run on " << fMaxEvents << " events..." << endl;
		nentries = fMaxEvents;
	}

	if( prescale>1 ) cout << " processing only every " << prescale << " events" << endl;
	
	for( Long64_t jentry = 0; jentry < nentries; jentry++ ){
		PrintProgress(jentry);

		// Prescale processing...
		if ( prescale>1 && jentry%prescale ) continue;

		fTR->GetEntry(jentry);

		// Upper Pt Hat cut
		if( fPtHatCut > -1.0 )
			if( fTR->PtHat > fPtHatCut ) continue;

		if( fCurRun != fTR->Run ) {
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