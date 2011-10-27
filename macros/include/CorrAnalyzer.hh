#ifndef CorrAnalyzer_hh
#define CorrAnalyzer_hh

#include <TTree.h>
#include <TStyle.h>

#include "base/TreeAnalyzerBase.hh"
#include "base/TreeReader.hh"
#include "ZeeAnalysis.hh"
#include "HggAnalysis.hh"

class CorrAnalyzer : public TreeAnalyzerBase {
public:
	CorrAnalyzer(TTree *tree = 0);
	virtual ~CorrAnalyzer();
	void BeginJob();
	void EndJob();
	void Loop();

private:
  ZeeAnalysis *fZeeAnalysis;
  HggAnalysis *fHggAnalysis;
  

};
#endif
