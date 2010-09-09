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
void PhysQCAnalyzer::Loop(Long64_t maxEvents, Int_t prescale){
	Long64_t nentries = fTR->GetEntries();
	cout << " total events in ntuples = " << nentries << endl;
        if ( maxEvents>=0 ) { 
          nentries = maxEvents;
          cout << " processing only " << nentries << " events out of it" << endl;
        }
        if ( prescale>1 ) {
          cout << " processing only every " << prescale << " events" << endl;
        }
	for( Long64_t jentry = 0; jentry < nentries; jentry++ ){
		PrintProgress(jentry);
                if ( prescale>1 && jentry%prescale ) continue; // Prescale processing...
		fTR->GetEntry(jentry);

                if ( fCurRun != fTR->Run ) {
                  fCurRun = fTR->Run;
                  fPhysQCAnalysis      ->BeginRun(fCurRun);
                  fTreeCleaner         ->BeginRun(fCurRun);
                  fMultiplicityAnalysis->BeginRun(fCurRun);
                }

		fPhysQCAnalysis      ->Analyze();
		fTreeCleaner         ->Analyze();
		fPhysQCAnalysis      ->Analyze2();
		fMultiplicityAnalysis->Analyze();
	}
}

// Method called before starting the event loop
void PhysQCAnalyzer::BeginJob(){
	// Initialize UserAnalyses here:
	fTreeCleaner->SetOutputDir(fOutputDir);
	fTreeCleaner->SetClean(true);
	fTreeCleaner->SetSkim(true);
	fTreeCleaner->SetVerbose(fVerbose);

	fPhysQCAnalysis->SetOutputDir(fOutputDir);
	fPhysQCAnalysis->SetVerbose(fVerbose);

	fMultiplicityAnalysis->SetOutputDir(fOutputDir);
	fMultiplicityAnalysis->SetVerbose(fVerbose);
	fMultiplicityAnalysis->ReadCuts("multiplicity_cuts/default.dat");
	
	fTreeCleaner         ->Begin();
	fPhysQCAnalysis      ->Begin();
	fMultiplicityAnalysis->Begin();
	

	fPhysQCAnalysis->PlotTriggerStats();
	TCut select = "(NMus + NEles + NJets) > 0";
	fPhysQCAnalysis->MakePlots("plots_uncleaned.dat", select, fTR->fChain);
	fPhysQCAnalysis->MakeElIDPlots(select, fTR->fChain);
}

// Method called after finishing the event loop
void PhysQCAnalyzer::EndJob(){
	TCut select = "GoodEvent == 0 || (NMus + NEles + NJets) > 0";
	fPhysQCAnalysis->MakePlots("plots_cleaned.dat", select, fTreeCleaner->fCleanTree);

	fTreeCleaner         ->End();
	fPhysQCAnalysis      ->End();
	fMultiplicityAnalysis->End();
}
