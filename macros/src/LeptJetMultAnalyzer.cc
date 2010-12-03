#include "LeptJetMultAnalyzer.hh"

#include "base/TreeAnalyzerBase.hh"
#include "base/TreeReader.hh"
#include "MultiplicityAnalysis.hh"
#include "MassAnalysis.hh"
#include "RatioAnalysis.hh"



using namespace std;

LeptJetMultAnalyzer::LeptJetMultAnalyzer(TTree *tree) : TreeAnalyzerBase(tree) {
	fMultiplicityAnalysis     = new MultiplicityAnalysis(fTR);
	fMassAnalysis             = new MassAnalysis(fTR);
	fRatioAnalysis            = new RatioAnalysis(fTR);
	Util::SetStyle();
}

LeptJetMultAnalyzer::~LeptJetMultAnalyzer(){
	delete fMultiplicityAnalysis;
	delete fMassAnalysis;
	delete fRatioAnalysis;
	if(!fTR->fChain) cout << "LeptJetMultAnalyzer ==> No chain!" << endl;
}

// Method for looping over the tree
void LeptJetMultAnalyzer::Loop(){
	Long64_t nentries = fTR->GetEntries();
	cout << " total events in ntuples = " << fTR->GetEntries() << endl;
	
	if(fMaxEvents==-1)nentries=fTR->GetEntries();
	if(fMaxEvents > 0){
		nentries=fMaxEvents;
		cout << " only running on first " << fMaxEvents << " events" << endl;
	}
	// loop over all ntuple entries
	for( Long64_t jentry = 0; jentry < nentries; jentry++ ){
		PrintProgress(jentry);
		fTR->GetEntry(jentry);
        	if ( fCurRun != fTR->Run ) {
        		fCurRun = fTR->Run;
			fMultiplicityAnalysis->BeginRun(fCurRun);
			fMassAnalysis        ->BeginRun(fCurRun);
			fRatioAnalysis       ->BeginRun(fCurRun);
		}
		fMassAnalysis        ->Analyze();
		fMultiplicityAnalysis->Analyze();	
		fRatioAnalysis       ->Analyze();	
	}
}

// Method called before starting the event loop
void LeptJetMultAnalyzer::BeginJob(TString filename, TString setofcuts, float lumi){
	
	fMultiplicityAnalysis     ->ReadCuts(setofcuts);
	fMultiplicityAnalysis     ->SetOutputDir(fOutputDir);
	fMultiplicityAnalysis     ->fVerbose        =fVerbose;	
	fMultiplicityAnalysis     ->fLumi           =lumi;
	fMultiplicityAnalysis     ->Begin();
	
	fMassAnalysis             ->ReadCuts(setofcuts);
	fMassAnalysis             ->SetOutputDir(fOutputDir);
	fMassAnalysis             ->fVerbose        = fVerbose;
	fMassAnalysis             ->Begin(filename);
	
	fRatioAnalysis            ->ReadCuts(setofcuts);
	fRatioAnalysis            ->SetOutputDir(fOutputDir);
	fRatioAnalysis            ->fVerbose        = fVerbose;
	fRatioAnalysis            ->Begin();


}

// Method called after finishing the event loop
void LeptJetMultAnalyzer::EndJob(){
	fMassAnalysis         ->End();
	fMultiplicityAnalysis ->End();
//	fRatioAnalysis        ->End();
}
