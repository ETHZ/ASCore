#include <stdlib.h>
#include "base/TreeAnalyzerBase.hh"
#include <TTree.h>
#include <TString.h>

using namespace std;

TreeAnalyzerBase::TreeAnalyzerBase(TTree *tree) {
	fTR = new TreeReader(tree);
	fNEntries = fTR->GetEntries();
	fVerbose = 0;
	fMaxEvents = -1;
	fCurRun = -1; // Initialise to dummy value
}

TreeAnalyzerBase::~TreeAnalyzerBase(){
	if(!fTR->fChain) cout << "TreeAnalyzerBase ==> No chain!" << endl;
	delete fTR;
}

// Method for looping over the tree
void TreeAnalyzerBase::Loop(){
	Long64_t nentries = fTR->GetEntries();
	cout << " total events in ntuples = " << fTR->GetEntries() << endl;
	if(fMaxEvents > -1){
		cout << " will run on " << fMaxEvents << " events..." << endl;
		nentries = fMaxEvents;
	}
	for( Long64_t jentry = 0; jentry < nentries; jentry++ ){
		PrintProgress(jentry);
		fTR->GetEntry(jentry);
	}

}

// Method called before starting the event loop
void TreeAnalyzerBase::BeginJob(){}

// Method called after finishing the event loop
void TreeAnalyzerBase::EndJob(){}

// Method that prints the progress in reasonable frequency
void TreeAnalyzerBase::PrintProgress(Long64_t entry){
	bool print = false;
	int step = 0;
	step = (fNEntries/20)/100;
	step *= 100;
	if( step < 200 ) step = 200;
	if( step > 10000 ) step = 10000;
	print = (entry%step == 0);
	if( print ) cout << ">>> Processing event # " << entry << endl;
}
