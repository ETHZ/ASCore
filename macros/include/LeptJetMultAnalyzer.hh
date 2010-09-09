#ifndef LeptJetMultAnalyzer_hh
#define LeptJetMultAnalyzer_hh

#include <TTree.h>
#include <TStyle.h>

#include "base/TreeAnalyzerBase.hh"
#include "base/TreeReader.hh"
#include "TreeCleaner.hh"
#include "MultiplicityAnalysis.hh"
#include "MassAnalysis.hh"


class LeptJetMultAnalyzer : public TreeAnalyzerBase {
public:
	LeptJetMultAnalyzer(TTree *tree = 0);
	virtual ~LeptJetMultAnalyzer();
	void BeginJob(TString filename="Multiplicity.root" , TString setofcuts="default", float lumi=-999.99, std::vector<std::string>* requiredHLT=NULL, std::vector<std::string>* vetoedHLT=NULL);
	void EndJob();
	void Loop();
	std::vector<std::string>* fRequiredHLT; 
	std::vector<std::string>* fVetoedHLT;

private:
	TreeCleaner              *fTreeCleaner;
	MultiplicityAnalysis     *fMultiplicityAnalysis;
	MassAnalysis             *fMassAnalysis;	
		
};
#endif
