#include "PhysQCAnalyzer.hh"

#include "TreeAnalyzerBase.hh"
#include "TreeReader.hh"
#include "TreeCleaner.hh"
#include "DiLeptonAnalysis.hh"
#include "MultiplicityAnalysis.hh"

using namespace std;

PhysQCAnalyzer::PhysQCAnalyzer(TTree *tree) : TreeAnalyzerBase(tree) {
	fTR = new TreeReader(tree);
	fTree = tree;
	fVerbose = 0;

	fTreeCleaner = new TreeCleaner(fTR);
	fPhysQCAnalysis = new PhysQCAnalysis(fTR);

	SetStyle();
}

PhysQCAnalyzer::~PhysQCAnalyzer(){
	delete fPhysQCAnalysis;
	delete fTreeCleaner;
	delete fTR;
	if(!fTree) cout << "PhysQCAnalyzer ==> No chain!" << endl;
}

// Method for looping over the tree
void PhysQCAnalyzer::Loop(){
	Long64_t nentries = fTree->GetEntries();
	cout << " total events in ntuples = " << fTree->GetEntries() << endl;
	// nentries = 10;
	for( Long64_t jentry = 0; jentry < nentries; jentry++ ){
		if( jentry%200 == 0 ) cout << ">>> Processing event # " << jentry << endl;
		fTree->GetEntry(jentry);

		fTreeCleaner   ->Analyze();
		fPhysQCAnalysis->Analyze();
	}
}

// Method called before starting the event loop
void PhysQCAnalyzer::BeginJob(){	
	// Initialize UserAnalyses here:
	fTreeCleaner->SetOutputDir(fOutputDir);
	fTreeCleaner->fClean = true;
	fTreeCleaner->fVerbose = fVerbose;

	fPhysQCAnalysis->SetOutputDir(fOutputDir);
	fPhysQCAnalysis->fVerbose = fVerbose;
	
	fTreeCleaner   ->Begin();
	fPhysQCAnalysis->Begin();
	fPhysQCAnalysis->PlotTriggerStats();
	fPhysQCAnalysis->MakePlots("plots_uncleaned.dat", fTree);
}

// Method called after finishing the event loop
void PhysQCAnalyzer::EndJob(){
	fTreeCleaner   ->End();
	// This is a stupid workaround, there must be a more clever way using fTreeCleaner->fCleanTree
	TFile *f = TFile::Open(fOutputDir + "CleanTree.root");
	f->cd("analyze");
	TTree *cleantree = (TTree*)gDirectory->Get("Analysis");
	fPhysQCAnalysis->MakePlots("plots_cleaned.dat", cleantree);
	fPhysQCAnalysis->End();
}
