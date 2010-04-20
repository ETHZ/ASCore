#include "PhysQCAnalyzer.hh"

#include "base/TreeAnalyzerBase.hh"
#include "base/TreeReader.hh"
#include "TreeCleaner.hh"
#include "DiLeptonAnalysis.hh"
#include "MultiplicityAnalysis.hh"

using namespace std;

PhysQCAnalyzer::PhysQCAnalyzer(TTree *tree) : TreeAnalyzerBase(tree) {
	fTreeCleaner          = new TreeCleaner(fTR);
	fPhysQCAnalysis       = new PhysQCAnalysis(fTR, fTreeCleaner);
	fMultiplicityAnalysis = new MultiplicityAnalysis(fTR);
}

PhysQCAnalyzer::~PhysQCAnalyzer(){
	delete fPhysQCAnalysis;
	delete fTreeCleaner;
	delete fMultiplicityAnalysis;
	if(!fTR->fChain) cout << "PhysQCAnalyzer ==> No chain!" << endl;
}

// Method for looping over the tree
void PhysQCAnalyzer::Loop(){
	Long64_t nentries = fTR->GetEntries();
	cout << " total events in ntuples = " << fTR->GetEntries() << endl;
	// nentries = 10;
	for( Long64_t jentry = 0; jentry < nentries; jentry++ ){
		PrintProgress(jentry);
		fTR->GetEntry(jentry);

		fPhysQCAnalysis      ->Analyze1();
		fTreeCleaner         ->Analyze();
		fPhysQCAnalysis      ->Analyze2();
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

	TCut select = "NMus+NEles+NPhotons+NJets>0";
	fPhysQCAnalysis->PlotTriggerStats();
	fPhysQCAnalysis->MakePlots("plots_uncleaned.dat", fTR->fChain, select);
	fPhysQCAnalysis->MakeElIDPlots(fTR->fChain, select);
}

// Method called after finishing the event loop
void PhysQCAnalyzer::EndJob(){
	TCut select = "GoodEvent==0 || (NMus+NEles+NPhotons+NJets)>0";
	fPhysQCAnalysis->MakePlots("plots_cleaned.dat", fTreeCleaner->fCleanTree, select);

	fTreeCleaner         ->End();
	fPhysQCAnalysis      ->End();
	fMultiplicityAnalysis->End();
}
