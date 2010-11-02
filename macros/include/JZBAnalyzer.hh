#ifndef JZBAnalyzer_hh
#define JZBAnalyzer_hh

#include <TTree.h>
#include <TStyle.h>

#include "base/TreeAnalyzerBase.hh"
#include "base/TreeReader.hh"
#include "JZBAnalysis.hh"
#include <string>

class JZBAnalyzer : public TreeAnalyzerBase {
public:
	JZBAnalyzer(TTree *tree = 0);
	virtual ~JZBAnalyzer();
	void BeginJob();
	void EndJob();
	void Loop();
	void SetMaxEvents(int a){fMaxEvents=a;}
	void SetOutputFile(TString a){fOutputFile=a;}
	void SetOutputFileName(string a){outputFileName_=a;}

private:
	JZBAnalysis *fJZBAnalysis;
        int fMaxEvents;
	TString fOutputFile;
	string outputFileName_;

};
#endif
