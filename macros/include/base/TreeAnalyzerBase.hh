#ifndef TreeAnalyzerBase_hh
#define TreeAnalyzerBase_hh

#include <TTree.h>
#include <TString.h>
#include "TreeReader.hh"
#include "helper/Utilities.hh"
#include <vector>

class TreeAnalyzerBase {
public:
	TreeAnalyzerBase(TTree *tree = 0);
	virtual ~TreeAnalyzerBase();
	virtual void BeginJob();
	virtual void EndJob();
	virtual void Loop();
	virtual void PrintProgress(Long64_t);

	inline virtual void SetVerbose(int verbose){fVerbose = verbose;};
	inline virtual void SetMaxEvents(Long64_t maxevents){fMaxEvents = maxevents;};
	inline virtual void SetData(bool isdata){fIsData = isdata;};

	inline void SetOutputDir(TString dir){ fOutputDir = Util::MakeOutputDir(dir); };
	inline void SetOutputFile(TString file){ fOutputFile = file; };
	
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

private:	
};
#endif
