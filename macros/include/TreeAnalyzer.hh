#ifndef TreeAnalyzer_hh
#define TreeAnalyzer_hh

#include <TTree.h>
#include <TStyle.h>

#include "base/TreeAnalyzerBase.hh"
#include "base/TreeReader.hh"
#include "TreeCleaner.hh"
#include "DiLeptonAnalysis.hh"
#include "MultiplicityAnalysis.hh"
#include "SignificanceAnalysis.hh"

class TreeAnalyzer : public TreeAnalyzerBase {
public:
	TreeAnalyzer(TTree *tree = NULL);
	virtual ~TreeAnalyzer();
	void BeginJob();
	void EndJob();
	void Loop();

private:

	TTree *fTree;
	TreeCleaner *fTreeCleaner;
	DiLeptonAnalysis *fDiLeptonAnalysis;
	MultiplicityAnalysis *fMultiplicityAnalysis;
	SignificanceAnalysis *fSignificanceAnalysis;
};
#endif
