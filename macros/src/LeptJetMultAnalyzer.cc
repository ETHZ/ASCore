#include "LeptJetMultAnalyzer.hh"

#include "base/TreeAnalyzerBase.hh"
#include "base/TreeReader.hh"
#include "MultiplicityAnalysis.hh"
#include "MassAnalysis.hh"


using namespace std;

LeptJetMultAnalyzer::LeptJetMultAnalyzer(TTree *tree) : TreeAnalyzerBase(tree) {
	fMultiplicityAnalysis     = new MultiplicityAnalysis(fTR);
	fMassAnalysis             = new MassAnalysis(fTR);
	Util::SetStyle();
}

LeptJetMultAnalyzer::~LeptJetMultAnalyzer(){
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
			fMultiplicityAnalysis->BeginRun(fCurRun);
			fMassAnalysis        ->BeginRun(fCurRun);
		}
		fMassAnalysis        ->Analyze();
		fMultiplicityAnalysis->Analyze();		
	}
}

// Method called before starting the event loop
void LeptJetMultAnalyzer::BeginJob(TString filename, TString setofcuts, float lumi){
	
	fMultiplicityAnalysis     ->ReadCuts(setofcuts);
	fMultiplicityAnalysis     ->SetOutputDir(fOutputDir);
	fMultiplicityAnalysis     ->fVerbose        =fVerbose;	
	fMultiplicityAnalysis     ->fLumi           =lumi;
	fMultiplicityAnalysis     ->Begin(filename);
	
	fMassAnalysis             ->ReadCuts(setofcuts);
	fMassAnalysis             ->SetOutputDir(fOutputDir);
	fMassAnalysis             ->fVerbose        = fVerbose;
	fMassAnalysis             ->Begin();


}

// Method called after finishing the event loop
void LeptJetMultAnalyzer::EndJob(){
	fMassAnalysis         ->End();
	fMultiplicityAnalysis ->End();
}
