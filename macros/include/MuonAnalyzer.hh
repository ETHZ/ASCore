#ifndef MuonAnalyzer_hh
#define MuonAnalyzer_hh

#include <TTree.h>
#include <TStyle.h>

#include "base/TreeAnalyzerBase.hh"
#include "base/TreeReader.hh"
#include "TreeCleaner.hh"
// #include "MuonFakeAnalysis.hh"
#include "MuonAnalysis.hh"

class MuonAnalyzer : public TreeAnalyzerBase {
public:
	MuonAnalyzer(TTree *tree = NULL);
	virtual ~MuonAnalyzer();
	void BeginJob();
	void EndJob();
	void Loop();

private:

	TTree *fTree;
	TreeCleaner *fTreeCleaner;
	MuonAnalysis *fMuonAnalysis;
	MuonAnalysis *fMuonAnalysisDi;
	MuonAnalysis *fMuonAnalysisSS;
};
#endif
