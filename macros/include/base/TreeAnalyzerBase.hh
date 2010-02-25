#ifndef TreeAnalyzerBase_hh
#define TreeAnalyzerBase_hh

#include <TTree.h>
#include <TString.h>
#include "TreeReader.hh"

class TreeAnalyzerBase {
public:
	TreeAnalyzerBase(TTree *tree = 0);
	virtual ~TreeAnalyzerBase();
	virtual void BeginJob();
	virtual void EndJob();
	virtual void Loop();

	virtual void SetStyle();
	virtual void SetOutputDir(TString);
	inline virtual void SetVerbose(int verbose){fVerbose = verbose;};

	TString fOutputDir;
	int fVerbose;

	TStyle *fStyle;
	TreeReader *fTR;

private:	
};
#endif
