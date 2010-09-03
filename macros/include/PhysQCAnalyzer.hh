#ifndef PhysQCAnalyzer_hh
#define PhysQCAnalyzer_hh

#include <TTree.h>
#include <TStyle.h>

#include "base/TreeAnalyzerBase.hh"
#include "base/TreeReader.hh"
#include "TreeCleaner.hh"
#include "PhysQCAnalysis.hh"
#include "helper/AnaClass.hh"
#include "MultiplicityAnalysis.hh"

class PhysQCAnalyzer : public TreeAnalyzerBase {
public:
	PhysQCAnalyzer(TTree *tree = 0);
	virtual ~PhysQCAnalyzer();
	void BeginJob();
	void EndJob();
	void Loop(Long64_t maxEvents=-1, Int_t prescale=1);

private:
	AnaClass *fAnaClass;
	TreeCleaner *fTreeCleaner;
	PhysQCAnalysis *fPhysQCAnalysis;
	MultiplicityAnalysis *fMultiplicityAnalysis;

};
#endif
