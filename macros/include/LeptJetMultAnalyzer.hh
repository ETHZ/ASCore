#ifndef LeptJetMultAnalyzer_hh
#define LeptJetMultAnalyzer_hh

#include <TTree.h>
#include <TStyle.h>

#include "base/TreeAnalyzerBase.hh"
#include "base/TreeReader.hh"
#include "MultiplicityAnalysis.hh"
#include "MassAnalysis.hh"
#include "RatioAnalysis.hh"



class LeptJetMultAnalyzer : public TreeAnalyzerBase {
public:
	LeptJetMultAnalyzer(TTree *tree = 0);
	virtual ~LeptJetMultAnalyzer();
	void BeginJob(TString filename="Multiplicity.root" , TString setofcuts="default", float lumi=-999.99, float weight=1.);
	void EndJob();
	void Loop();

private:
	MultiplicityAnalysis     *fMultiplicityAnalysis;
	MassAnalysis             *fMassAnalysis;
	RatioAnalysis            *fRatioAnalysis;		
		
};
#endif
