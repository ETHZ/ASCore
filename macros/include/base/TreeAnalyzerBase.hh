#ifndef TreeAnalyzerBase_hh
#define TreeAnalyzerBase_hh

#include <TString.h>
#include "TreeReader.hh"
#include "helper/Utilities.hh"
#include <vector>
#include <string>

// FWLite includes
#include "DataFormats/FWLite/interface/Event.h"

class TreeAnalyzerBase {
public:
    TreeAnalyzerBase(std::vector<std::string>& fileList);
    virtual ~TreeAnalyzerBase();
    virtual void BeginJob() {} // Method called before starting the event loop
    virtual void EndJob() {} // Method called after finishing the event loop
    virtual void Loop() = 0;
    virtual void PrintProgress(Long64_t);

    inline virtual void SetVerbose(int verbose){fVerbose = verbose;};
    inline virtual void SetMaxEvents(Long64_t maxevents){fMaxEvents = maxevents;};
    inline virtual void SetData(bool isdata){fIsData = isdata;};

    inline virtual void SetOutputDir(TString dir){ fOutputDir = Util::MakeOutputDir(dir); };
    inline virtual void SetOutputFile(TString file){ fOutputFile = file; };
	
    bool fIsData;
    TString fOutputDir;
    TString fOutputFile;
    int fVerbose;
    int fNEntries;
    Long64_t fMaxEvents;
  

    Int_t fCurRun;
    Int_t fCurLumi;
    bool skipLumi;
    bool skipRun;

    // stuff for JSON reading
    virtual void ReadJSON(const char* JSONpath);
    virtual const bool CheckRunLumi(void) const;
    virtual const bool CheckRun(void) const;
    struct RunLumi{
        int run;
        std::vector<int> lumi_min;
        std::vector<int> lumi_max;
    } fCurRunLumi; 
    std::vector<RunLumi> fRunLumis;


    TreeReader *fTR;


};
#endif
