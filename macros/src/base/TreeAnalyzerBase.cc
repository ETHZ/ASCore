#include <stdlib.h>
#include "base/TreeAnalyzerBase.hh"
#include <TTree.h>
#include <TString.h>

using namespace std;

TreeAnalyzerBase::TreeAnalyzerBase(TTree *tree) {
	fTR = new TreeReader(tree);
	fVerbose = false;
}

TreeAnalyzerBase::~TreeAnalyzerBase(){
	if(!fTR->fChain) cout << "TreeAnalyzerBase ==> No chain!" << endl;
}

// Method for looping over the tree
void TreeAnalyzerBase::Loop(){
	Long64_t nentries = fTR->GetEntries();
	cout << " total events in ntuples = " << fTR->GetEntries() << endl;
	// nentries = 10;
	for( Long64_t jentry = 0; jentry < nentries; jentry++ ){
		if( jentry%200 == 0 ) cout << ">>> Processing event # " << jentry << endl;
		fTR->GetEntry(jentry);

	}
}

// Method called before starting the event loop
void TreeAnalyzerBase::BeginJob(){}

// Method called after finishing the event loop
void TreeAnalyzerBase::EndJob(){}
