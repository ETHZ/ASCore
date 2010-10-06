#include "MuonAnalyzer.hh"
#include "base/TreeAnalyzerBase.hh"
#include "base/TreeReader.hh"
#include "MuonAnalysis.hh"
#include "helper/Utilities.hh"

using namespace std;

MuonAnalyzer::MuonAnalyzer(TTree *tree) : TreeAnalyzerBase(tree) {
	Util::SetStyle();

	// Initialize UserAnalyses here:
	fMuonAnalysis = new MuonAnalysis(fTR);
	fMuonAnalysisDi = new MuonAnalysis(fTR);
	fMuonAnalysisSS = new MuonAnalysis(fTR);
}

MuonAnalyzer::~MuonAnalyzer(){
	delete fMuonAnalysis;
	delete fMuonAnalysisDi;
	delete fMuonAnalysisSS;
}

// Method for looping over the tree
void MuonAnalyzer::Loop(){
	Long64_t nentries = fTR->GetEntries();
	cout << " total events in ntuples = " << fTR->GetEntries() << endl;
	if(fMaxEvents > -1){
		cout << " will run on " << fMaxEvents << " events..." << endl;
		nentries = fMaxEvents;
	}
	for( Long64_t jentry = 0; jentry < nentries; jentry++ ){
		PrintProgress(jentry);
		fTR->GetEntry(jentry);

		if( fCurRun != fTR->Run ){
			fCurRun = fTR->Run;
			fMuonAnalysis->BeginRun(fCurRun);
			fMuonAnalysisDi->BeginRun(fCurRun);
			fMuonAnalysisSS->BeginRun(fCurRun);
		}

		fMuonAnalysis->Analyze();
		fMuonAnalysisDi->AnalyzeDi();
		fMuonAnalysisSS->AnalyzeSS();
	}
}

// Method called before starting the event loop
void MuonAnalyzer::BeginJob(){
	fMuonAnalysis->SetOutputDir(fOutputDir);
	fMuonAnalysis->SetOutputFile("SingleMuonSelection");
	fMuonAnalysis->SetVerbose(fVerbose);
	
	fMuonAnalysisDi->SetOutputDir(fOutputDir);
	fMuonAnalysisDi->SetOutputFile("DiMuonSelection");
	fMuonAnalysisDi->SetVerbose(fVerbose);
	
	fMuonAnalysisSS->SetOutputDir(fOutputDir);
	fMuonAnalysisSS->SetOutputFile("SSMuonSelection");
	fMuonAnalysisSS->SetVerbose(fVerbose);

	fMuonAnalysis->Begin();
	fMuonAnalysisDi->Begin();
	fMuonAnalysisSS->Begin();
}

// Method called after finishing the event loop
void MuonAnalyzer::EndJob(){
	fMuonAnalysis->End();
	fMuonAnalysisDi->End();
	fMuonAnalysisSS->End();
}
