#include "TreeSkimmer.hh"

#include "base/TreeAnalyzerBase.hh"
#include "base/TreeReader.hh"

#include <TTree.h>
#include <TFile.h>

using namespace std;

TreeSkimmer::TreeSkimmer(TTree *tree) : TreeAnalyzerBase(tree){
	fTree = tree;
}

TreeSkimmer::~TreeSkimmer(){
}

// Method for looping over the tree
void TreeSkimmer::Loop(){
	Long64_t nentries = fTree->GetEntries();
	cout << " total events in ntuples = " << fTree->GetEntries() << endl;
	// nentries = 10;
	for( Long64_t jentry = 0; jentry < nentries; jentry++ ){
		PrintProgress(jentry);

		fTree->GetEntry(jentry);

		fNtot++;
		if(!EventSelection()) continue;
		fNsel++;
		if(fVerbose > 1) cout << " TreeSkimmer ==> Run/Event " << fTR->Run << "/" << fTR->Event << " selected, " << fNsel << " so far" << endl;
		fSkimmedTree->Fill();
	}
}

bool TreeSkimmer::EventSelection(){
	///////////////////////////////
	// YOUR EVENT SELECTION HERE //
	///////////////////////////////
	if(fTR->NJets < 10) return false;
	return true;
}

// Method called before starting the event loop
void TreeSkimmer::BeginJob(){
	fSkimmedTreeFile = new TFile(fOutputDir + "SkimmedTree.root", "RECREATE");
	fSkimmedTreeFile->mkdir("analyze", "analyze");
	fSkimmedTreeFile->cd("analyze");
	fSkimmedTree = fTree->CloneTree(0);
	fSkimmedTree->CopyAddresses(fTree);
	
	fNsel = 0;
	fNtot = 0;
}

// Method called after finishing the event loop
void TreeSkimmer::EndJob(){
	if(fVerbose > 0) cout << "--------------------------------------------------------------------" << endl;
	if(fVerbose > 0) cout << " TreeSkimmer > Events selected / Total Events: " << fNsel << "/" << fNtot << endl;
	if(fVerbose > 0) cout << "--------------------------------------------------------------------" << endl;
	
	fSkimmedTreeFile->cd("analyze");
	fSkimmedTree->Write();
	fSkimmedTreeFile->Close();
}
