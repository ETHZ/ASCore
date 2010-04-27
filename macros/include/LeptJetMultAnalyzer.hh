#ifndef LeptJetMultAnalyzer_hh
#define LeptJetMultAnalyzer_hh

#include <TTree.h>
#include <TStyle.h>

#include "base/TreeAnalyzerBase.hh"
#include "base/TreeReader.hh"
#include "TreeCleaner.hh"
#include "MultiplicityAnalysis.hh"


class LeptJetMultAnalyzer : public TreeAnalyzerBase {
public:
	LeptJetMultAnalyzer(TTree *tree = 0);
	virtual ~LeptJetMultAnalyzer();
	void BeginJob(TString filename="Multiplicity.root" , TString setofcuts="default", float lumi=-999.99);
	void EndJob();
	void Loop();

private:
	TreeCleaner *fTreeCleaner;
	MultiplicityAnalysis *fMultiplicityAnalysis;	
};
#endif
