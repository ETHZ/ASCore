#include "LeptJetMultAnalyzer.hh"

#include "base/TreeAnalyzerBase.hh"
#include "base/TreeReader.hh"
#include "MultiplicityAnalysis.hh"
#include "TreeCleaner.hh"


using namespace std;

LeptJetMultAnalyzer::LeptJetMultAnalyzer(TTree *tree) : TreeAnalyzerBase(tree) {
	fTreeCleaner          = new TreeCleaner(fTR);
	fMultiplicityAnalysis = new MultiplicityAnalysis(fTR);
	Util::SetStyle();
}

LeptJetMultAnalyzer::~LeptJetMultAnalyzer(){
	delete fTreeCleaner;
	delete fMultiplicityAnalysis;
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
                }
		fTreeCleaner         ->Analyze();
		fMultiplicityAnalysis->Analyze();		
	}
}

// Method called before starting the event loop
void LeptJetMultAnalyzer::BeginJob(TString filename, TString setofcuts, float lumi){
	
	fTreeCleaner->SetOutputDir(fOutputDir);
	fTreeCleaner->fClean = true;       // is not specified: true, set to false if not-cleaned objects wanted. 
	fTreeCleaner->fVerbose = fVerbose;
	
	fMultiplicityAnalysis->SetOutputDir(fOutputDir);
	fMultiplicityAnalysis->fVerbose = fVerbose;	
	

	fTreeCleaner         ->Begin();
	fMultiplicityAnalysis->fSetofCuts      =setofcuts;
	fMultiplicityAnalysis->fLumi           =lumi;
	fMultiplicityAnalysis->Begin(filename);

}

// Method called after finishing the event loop
void LeptJetMultAnalyzer::EndJob(){
	fTreeCleaner          ->End();
	fMultiplicityAnalysis ->End();
}
