#ifndef LeptJetMultAnalyzer_hh
#define LeptJetMultAnalyzer_hh

#include <TTree.h>
#include <TStyle.h>

#include "base/TreeAnalyzerBase.hh"
#include "base/TreeReader.hh"
#include "MassAnalysis.hh"


class LeptJetMultAnalyzer : public TreeAnalyzerBase {
public:
	LeptJetMultAnalyzer(TTree *tree = 0);
	virtual ~LeptJetMultAnalyzer();
	void BeginJob(TString filename="Multiplicity.root" , TString setofcuts="default", float lumi=-999.99, 
	              bool isData=false, string data_PileUp="", string mc_PileUp="");
	void EndJob();
	void Loop();
	void SetMaxEvents(int a){fMaxEvents=a;}
  	bool isS3;
  bool noPU;  
private:
	MassAnalysis             *fMassAnalysis;
  	int fMaxEvents;
};
#endif
