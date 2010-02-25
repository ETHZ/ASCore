#include "PhysQCAnalyzer.hh"

#include "TreeAnalyzerBase.hh"
#include "TreeReader.hh"
#include "TreeCleaner.hh"
#include "DiLeptonAnalysis.hh"
#include "MultiplicityAnalysis.hh"

using namespace std;

PhysQCAnalyzer::PhysQCAnalyzer(TTree *tree) : TreeAnalyzerBase(tree) {
	fTR = new TreeReader(tree);
	fVerbose = 0;
	fTreeCleaner          = new TreeCleaner(fTR);
	fPhysQCAnalysis       = new PhysQCAnalysis(fTR, fTreeCleaner);
	fMultiplicityAnalysis = new MultiplicityAnalysis(fTR);
}

PhysQCAnalyzer::~PhysQCAnalyzer(){
	delete fPhysQCAnalysis;
	delete fTreeCleaner;
	delete fMultiplicityAnalysis;
	delete fTR;
	if(!fTR->fChain) cout << "PhysQCAnalyzer ==> No chain!" << endl;
}

// Method for looping over the tree
void PhysQCAnalyzer::Loop(){
	Long64_t nentries = fTR->GetEntries();
	cout << " total events in ntuples = " << fTR->GetEntries() << endl;
	// nentries = 10;
	for( Long64_t jentry = 0; jentry < nentries; jentry++ ){
		if( jentry%200 == 0 ) cout << ">>> Processing event # " << jentry << endl;
		fTR->GetEntry(jentry);
		fTreeCleaner         ->Analyze();
		fPhysQCAnalysis      ->Analyze();
		fMultiplicityAnalysis->Analyze();
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

	fMultiplicityAnalysis->SetOutputDir(fOutputDir);
	fMultiplicityAnalysis->fVerbose = fVerbose;
	
	fTreeCleaner         ->Begin();
	fPhysQCAnalysis      ->Begin();
	fMultiplicityAnalysis->Begin();

	fPhysQCAnalysis->PlotTriggerStats();
	fPhysQCAnalysis->MakePlots("plots_uncleaned.dat", fTR->fChain);
}

// Method called after finishing the event loop
void PhysQCAnalyzer::EndJob(){
	fPhysQCAnalysis->MakePlots("plots_cleaned.dat", fTreeCleaner->fCleanTree);

	fTreeCleaner         ->End();
	fPhysQCAnalysis      ->End();
	fMultiplicityAnalysis->End();
}
