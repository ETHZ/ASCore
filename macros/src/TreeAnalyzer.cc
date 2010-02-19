#include "TreeAnalyzer.hh"

#include "TreeAnalyzerBase.hh"
#include "TreeReader.hh"
#include "TreeCleaner.hh"
#include "DiLeptonAnalysis.hh"
#include "MultiplicityAnalysis.hh"

using namespace std;

TreeAnalyzer::TreeAnalyzer(TTree *tree) : TreeAnalyzerBase(tree) {
	fTR = new TreeReader(tree);
	fTree = tree;
	fVerbose = false;

	SetStyle();

	// Initialize UserAnalyses here:
	fTreeCleaner          = new TreeCleaner(fTR);
	fDiLeptonAnalysis     = new DiLeptonAnalysis(fTR);
	fMultiplicityAnalysis = new MultiplicityAnalysis(fTR);
	fSignificanceAnalysis = new SignificanceAnalysis(fTR);
}

TreeAnalyzer::~TreeAnalyzer(){
	delete fTreeCleaner;
	delete fDiLeptonAnalysis;
	delete fMultiplicityAnalysis;
	delete fSignificanceAnalysis;
	delete fTR;
	if(!fTree) cout << "TreeAnalyzer ==> No chain!" << endl;
}

// Method for looping over the tree
void TreeAnalyzer::Loop(){
	Long64_t nentries = fTree->GetEntries();
	cout << " total events in ntuples = " << fTree->GetEntries() << endl;
	// nentries = 10;
	for( Long64_t jentry = 0; jentry < nentries; jentry++ ){
		if( jentry%200 == 0 ) cout << ">>> Processing event # " << jentry << endl;
		fTree->GetEntry(jentry);

		// cout << &(fTR->TCMET) << endl;

		fTreeCleaner         ->Analyze();
		fDiLeptonAnalysis    ->Analyze();
		fMultiplicityAnalysis->Analyze();
		fSignificanceAnalysis->Analyze();
	}
}

// Method called before starting the event loop
void TreeAnalyzer::BeginJob(){
	fTreeCleaner->SetOutputDir(fOutputDir);
	fTreeCleaner->fClean = true;
	fTreeCleaner->fVerbose = fVerbose;
	
	fDiLeptonAnalysis->SetOutputDir(fOutputDir);
	fDiLeptonAnalysis->fVerbose = fVerbose;
	
	fMultiplicityAnalysis->SetOutputDir(fOutputDir);
	fMultiplicityAnalysis->fVerbose = fVerbose;
	
	fSignificanceAnalysis->SetOutputDir(fOutputDir);
	fSignificanceAnalysis->fVerbose = fVerbose;
	
	fTreeCleaner         ->Begin();
	fDiLeptonAnalysis    ->Begin();
	fMultiplicityAnalysis->Begin();
	fSignificanceAnalysis->Begin();
}

// Method called after finishing the event loop
void TreeAnalyzer::EndJob(){
	fTreeCleaner         ->End();
	fDiLeptonAnalysis    ->End();
	fMultiplicityAnalysis->End();
	fSignificanceAnalysis->End();
}
