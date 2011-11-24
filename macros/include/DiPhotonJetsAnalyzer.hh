#ifndef DiPhotonJetsAnalyzer_hh
#define DiPhotonJetsAnalyzer_hh

#include <TTree.h>
#include <TStyle.h>

#include "base/TreeAnalyzerBase.hh"
#include "base/TreeReader.hh"
#include "DiPhotonPurity.hh"
#include "DiPhotonMiniTree.hh"
#include <string>

class DiPhotonJetsAnalyzer : public TreeAnalyzerBase {
public:
  DiPhotonJetsAnalyzer(TTree *tree = 0, std::string dataType="data", double aw=-999, double* _kfac=NULL);
	virtual ~DiPhotonJetsAnalyzer();
	void BeginJob(string data_PileUp, string mc_PileUp);
	void EndJob();
	void Loop();
  void SetMaxEvents(int a){fMaxEvents=a;}

private:
  DiPhotonPurity *fDiPhotonPurity;
  DiPhotonMiniTree *fDiPhotonMiniTree;
  int fMaxEvents;
  double AddWeight;
  double* kfactors;

};
#endif
