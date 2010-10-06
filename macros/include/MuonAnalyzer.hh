#ifndef MuonAnalyzer_hh
#define MuonAnalyzer_hh

#include "base/TreeAnalyzerBase.hh"
#include "MuonAnalysis.hh"

class MuonAnalyzer : public TreeAnalyzerBase {
public:
	MuonAnalyzer(TTree *tree = NULL);
	virtual ~MuonAnalyzer();
	void BeginJob();
	void EndJob();
	void Loop();

private:

	MuonAnalysis *fMuonAnalysis;
	MuonAnalysis *fMuonAnalysisDi;
	MuonAnalysis *fMuonAnalysisSS;
};
#endif
