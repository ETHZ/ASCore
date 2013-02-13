#include "SSDLAnalyzer.hh"
#include "base/TreeReader.hh"

using namespace std;

SSDLAnalyzer::SSDLAnalyzer(std::vector<std::string>& fileList) : TreeAnalyzerBase(fileList) {
	fSSDLAnalysis = new SSDLAnalysis(fTR);
	fDoFillEffTree = false;
}

// MARC SSDLAnalyzer::SSDLAnalyzer(TTree *tree) : TreeAnalyzerBase(tree) {
// MARC 	fSSDLAnalysis = new SSDLAnalysis(fTR);
// MARC 	fDoFillEffTree = false;
// MARC }

SSDLAnalyzer::~SSDLAnalyzer(){
	delete fSSDLAnalysis;
	// MARC if(!fTR->fChain) cout << "SSDLAnalyzer ==> No chain!" << endl;
}

// Method for looping over the tree
void SSDLAnalyzer::Loop(){
	Long64_t nentries = fTR->GetEntries();
	cout << " total events in ntuples = " << nentries << endl;

	if( fMaxEvents > -1 ){
		cout << " will run on " << fMaxEvents << " events..." << endl;
		nentries = fMaxEvents;
	}

	// MARC if( prescale>1 ) cout << " processing only every " << prescale << " events" << endl;
	// MARC for( Long64_t jentry = 0; jentry < nentries; jentry++ ){

    Long64_t jentry=0;
	for ( fTR->ToBegin(); 
			!(fTR->AtEnd()) && (jentry<fMaxEvents || fMaxEvents<0); 
			++(*fTR) ) 
		{
			PrintProgress(jentry++);

			// MARC PrintProgress(jentry);
			// MARC // Prescale processing...
			// MARC if ( prescale>1 && jentry%prescale ) continue;

			fTR->GetEntry(jentry);

			// Upper Pt Hat cut
			if( (fPtHatCut > -1.0) && (fTR->PtHat > fPtHatCut) ) continue;

			// Run processing
			if( fCurRun != fTR->Run ) { // new run
				fCurRun = fTR->Run;
				skipRun = false;
				if ( CheckRun() == false ) skipRun = true;
				else fSSDLAnalysis->BeginRun(fCurRun);
				fSSDLAnalysis->BeginRun(fCurRun);
			}
			
			// Check if new lumi is in JSON file
			if( fCurLumi != fTR->LumiSection ) { // new lumisection
				fCurLumi = fTR->LumiSection;
				skipLumi = false;
				if ( CheckRunLumi() == false ) skipLumi = true;
			}
			// disable lumi checking in SSDLAnalysis if(skipRun || skipLumi) continue;
			fSSDLAnalysis->Analyze();
	}
	cout << endl;
}

// Method called before starting the event loop
// there are two function for with/without PU reweighting
void SSDLAnalyzer::BeginJob(){
	// fSSDLAnalysis->SetOutputDir(fOutputDir);
	fSSDLAnalysis->SetOutputFile(fOutputFile);
	fSSDLAnalysis->SetVerbose(fVerbose);
	fSSDLAnalysis->SetData(fIsData);
	fSSDLAnalysis->DoFillEffTree(fDoFillEffTree);
	fSSDLAnalysis->Begin();
}

void SSDLAnalyzer::BeginJob(std::string dataPuFile, std::string mcPuFile){
	// fSSDLAnalysis->SetOutputDir(fOutputDir);
	fSSDLAnalysis->SSDLAnalysis::SetPileUpSrc( dataPuFile , mcPuFile); //HERE
	fSSDLAnalysis->SetOutputFile(fOutputFile);
	fSSDLAnalysis->SetVerbose(fVerbose);
	fSSDLAnalysis->SetData(fIsData);
	fSSDLAnalysis->DoFillEffTree(fDoFillEffTree);
	fSSDLAnalysis->Begin();
}

// Method called after finishing the event loop
void SSDLAnalyzer::EndJob(){
	fSSDLAnalysis->End();
}
