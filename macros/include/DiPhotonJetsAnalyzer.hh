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
  DiPhotonJetsAnalyzer(TTree *tree = 0, std::string dataType="data", Float_t aw=-999, Float_t* _kfac=NULL, Float_t minthrpfphotoncandEE=0);
	virtual ~DiPhotonJetsAnalyzer();
	void BeginJob(string data_PileUp, string mc_PileUp);
	void EndJob();
	void Loop();
  void SetMaxEvents(int a){fMaxEvents=a;}

private:
  DiPhotonPurity *fDiPhotonPurity;
  DiPhotonMiniTree *fDiPhotonMiniTree;
  int fMaxEvents;
  Float_t AddWeight;
  Float_t* kfactors;
  std::string templateChoice;
  Float_t minthrpfphotoncandEE;

};
#endif
