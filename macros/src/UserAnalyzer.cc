#include "UserAnalyzer.hh"
#include "UserAnalysis.hh"

#include "base/TreeAnalyzerBase.hh"
#include "base/TreeReader.hh"

using namespace std;

UserAnalyzer::UserAnalyzer(std::vector<std::string>& fileList) 
    : TreeAnalyzerBase(fileList) {
    fUserAnalysis = new UserAnalysis(fTR);
}

UserAnalyzer::~UserAnalyzer(){
    delete fUserAnalysis;
}

// Method for looping over the tree
void UserAnalyzer::Loop(){
    Long64_t nentries = fTR->GetEntries();
    cout << " total events in ntuples = " << fTR->GetEntries() << endl;
	
    // loop over all ntuple entries
    // nentries = 200;
    Long64_t jentry=0;
    for ( fTR->ToBegin(); 
          !(fTR->AtEnd()) && (jentry<fMaxEvents || fMaxEvents<0); 
          ++(*fTR) ) 
        {
            PrintProgress(jentry++);
            if ( fCurRun != fTR->Run ) {
                fCurRun = fTR->Run;
                fUserAnalysis->BeginRun(fCurRun);
            }
            fUserAnalysis->Analyze();
        }
}

// Method called before starting the event loop
void UserAnalyzer::BeginJob(){
    fUserAnalysis->SetOutputDir(fOutputDir);
    fUserAnalysis->SetVerbose(fVerbose);

    fUserAnalysis->Begin();

}

// Method called after finishing the event loop
void UserAnalyzer::EndJob(){
    fUserAnalysis->End();
}
