#include "LeptJetMultAnalyzer.hh"

#include "base/TreeAnalyzerBase.hh"
#include "base/TreeReader.hh"
#include "MultiplicityAnalysis.hh"
#include "MassAnalysis.hh"
#include "TreeCleaner.hh"


using namespace std;

LeptJetMultAnalyzer::LeptJetMultAnalyzer(TTree *tree) : TreeAnalyzerBase(tree) {
	fTreeCleaner          = new TreeCleaner(fTR);
	fMultiplicityAnalysis = new MultiplicityAnalysis(fTR);
	fMassAnalysis         = new MassAnalysis(fTR);
	Util::SetStyle();
}

LeptJetMultAnalyzer::~LeptJetMultAnalyzer(){
	delete fTreeCleaner;
	delete fMultiplicityAnalysis;
	delete fMassAnalysis;
	if(!fTR->fChain) cout << "LeptJetMultAnalyzer ==> No chain!" << endl;
}

// Method for looping over the tree
void LeptJetMultAnalyzer::Loop(){
	Long64_t nentries = fTR->GetEntries();
	cout << " total events in ntuples = " << fTR->GetEntries() << endl;


	
	// loop over all ntuple entries
	for( Long64_t jentry = 0; jentry < nentries; jentry++ ){
		PrintProgress(jentry);
		fTR->GetEntry(jentry);
        if ( fCurRun != fTR->Run ) {
        	fCurRun = fTR->Run;
			fTreeCleaner         ->BeginRun(fCurRun);
			fMultiplicityAnalysis->BeginRun(fCurRun);
			fMassAnalysis        ->BeginRun(fCurRun);
		}
		fTreeCleaner         ->Analyze();
		fMassAnalysis        ->Analyze();
		fMultiplicityAnalysis->Analyze();		
	}
}

// Method called before starting the event loop
void LeptJetMultAnalyzer::BeginJob(TString filename, TString setofcuts, float lumi, std::vector<std::string>* requiredHLT, std::vector<std::string>* vetoedHLT){
	
	fTreeCleaner         ->SetOutputDir(fOutputDir);
	fTreeCleaner         ->fClean = true;       // if not specified: true, set to false if not-cleaned objects wanted. 
	fTreeCleaner         ->SetSkim(false);      // do not skim the tree! 
	fTreeCleaner         ->fVerbose = fVerbose;
	fTreeCleaner         ->Begin();
	
	fMultiplicityAnalysis->SetOutputDir(fOutputDir);
	fMultiplicityAnalysis->fVerbose        =fVerbose;	
	fMultiplicityAnalysis->fRequiredHLT    =requiredHLT; 
	fMultiplicityAnalysis->fVetoedHLT      =vetoedHLT; 
	fMultiplicityAnalysis->fSetofCuts      =setofcuts;
	fMultiplicityAnalysis->fLumi           =lumi;
	fMultiplicityAnalysis->Begin(filename);
	
	fMassAnalysis        ->SetOutputDir(fOutputDir);
	fMassAnalysis        ->fVerbose        = fVerbose;
	fMassAnalysis        ->Begin();


}

// Method called after finishing the event loop
void LeptJetMultAnalyzer::EndJob(){
	fTreeCleaner          ->End();
	fMassAnalysis         ->End();
	fMultiplicityAnalysis ->End();
}
