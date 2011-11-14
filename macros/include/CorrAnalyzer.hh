#ifndef CorrAnalyzer_hh
#define CorrAnalyzer_hh

#include <TTree.h>
#include <TStyle.h>

#include "base/TreeAnalyzerBase.hh"
#include "base/TreeReader.hh"
#include "ZeeAnalysis.hh"
#include "ZeeMiniTree.hh"
#include "HggAnalysis.hh"
#include <string>

class CorrAnalyzer : public TreeAnalyzerBase {
public:
  CorrAnalyzer(TTree *tree = 0, std::string dataType="data");
	virtual ~CorrAnalyzer();
	void BeginJob(string data_PileUp, string mc_PileUp);
	void EndJob();
	void Loop();
  void SetMaxEvents(int a){fMaxEvents=a;}

private:
  ZeeAnalysis *fZeeAnalysis;
  ZeeMiniTree *fZeeMiniTree;
  HggAnalysis *fHggAnalysis;
  int fMaxEvents;

};
#endif
