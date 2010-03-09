#ifndef TreeAnalyzerBase_hh
#define TreeAnalyzerBase_hh

#include <TTree.h>
#include <TString.h>
#include "TreeReader.hh"
#include "helper/Utilities.hh"

class TreeAnalyzerBase {
public:
	TreeAnalyzerBase(TTree *tree = 0);
	virtual ~TreeAnalyzerBase();
	virtual void BeginJob();
	virtual void EndJob();
	virtual void Loop();
	virtual void PrintProgress(Long64_t);

	inline virtual void SetVerbose(int verbose){fVerbose = verbose;};
	inline void SetOutputDir(TString dir){ fOutputDir = Util::MakeOutputDir(dir); };
	
	TString fOutputDir;
	int fVerbose;
	int fNEntries;

	TreeReader *fTR;

private:	
};
#endif
