#include "DiPhotonJetsAnalyzer.hh"
#include "DiPhotonPurity.hh"
#include "DiPhotonMiniTree.hh"


#include "base/TreeAnalyzerBase.hh"
#include "base/TreeReader.hh"

using namespace std;

DiPhotonJetsAnalyzer::DiPhotonJetsAnalyzer(TTree *tree, std::string dataType, Float_t aw, Float_t* _kfac) : TreeAnalyzerBase(tree), AddWeight(aw), kfactors(_kfac) {
  //  fDiPhotonPurity = new DiPhotonPurity(fTR,dataType,AddWeight);
  fDiPhotonMiniTree = new DiPhotonMiniTree(fTR,dataType,AddWeight,kfactors);
}

DiPhotonJetsAnalyzer::~DiPhotonJetsAnalyzer(){
  //	delete fDiPhotonPurity;
	delete fDiPhotonMiniTree;
	if(!fTR->fChain) cout << "DiPhotonJetsAnalyzer ==> No chain!" << endl;
}

// Method for looping over the tree
void DiPhotonJetsAnalyzer::Loop(){
	Long64_t nentries = fTR->GetEntries();
	cout << " total events in ntuples = " << fTR->GetEntries() << endl;

	if(fMaxEvents==-1)nentries=fTR->GetEntries();
	if(fMaxEvents>0)nentries=fMaxEvents;
	
	// loop over all ntuple entries
	// nentries = 200;
	for( Long64_t jentry = 0; jentry < nentries; jentry++ ){
		PrintProgress(jentry);
		fTR->GetEntry(jentry);
		if ( fCurRun != fTR->Run ) {
		  fCurRun = fTR->Run;
		  //		  fDiPhotonPurity->BeginRun(fCurRun);
		  fDiPhotonMiniTree->BeginRun(fCurRun);
		  skipRun = false;
		  if ( !CheckRun() ) skipRun = true;
		}
		// Check if new lumi is in JSON file
		if ( !skipRun && fCurLumi != fTR->LumiSection ) {
		  fCurLumi = fTR->LumiSection;
		  skipLumi = false; // Re-initialise
		  if ( !CheckRunLumi() ) skipLumi = true;
		}
		if ( !(skipRun || skipLumi) ) {
		  //		  fDiPhotonPurity->Analyze();
		  fDiPhotonMiniTree->Analyze();
		}
	}
}

// Method called before starting the event loop
void DiPhotonJetsAnalyzer::BeginJob(string fdata_PileUp, string fmc_PileUp){
// 	fDiPhotonPurity->SetOutputDir(fOutputDir);
// 	fDiPhotonPurity->SetOutputFile(fOutputFile);
// 	fDiPhotonPurity->fVerbose = fVerbose;
// 	fDiPhotonPurity->SetPileUpSrc(fdata_PileUp, fmc_PileUp);
// 	fDiPhotonPurity->Begin();

	fDiPhotonMiniTree->SetOutputDir(fOutputDir);
	fDiPhotonMiniTree->SetOutputFile(fOutputFile);
	fDiPhotonMiniTree->fVerbose = fVerbose;
	fDiPhotonMiniTree->SetPileUpSrc(fdata_PileUp, fmc_PileUp);
	fDiPhotonMiniTree->Begin();

}

// Method called after finishing the event loop
void DiPhotonJetsAnalyzer::EndJob(){
  //	fDiPhotonPurity->End();
	fDiPhotonMiniTree->End();
}

