#ifndef PhysQCAnalyzer_hh
#define PhysQCAnalyzer_hh

#include <TTree.h>
#include <TStyle.h>

#include "TreeAnalyzerBase.hh"
#include "TreeReader.hh"
#include "TreeCleaner.hh"
#include "PhysQCAnalysis.hh"
#include "MultiplicityAnalysis.hh"
#include "AnaClass.hh"

class PhysQCAnalyzer : public TreeAnalyzerBase {
public:
	PhysQCAnalyzer(TTree *tree = 0);
	virtual ~PhysQCAnalyzer();
	void BeginJob();
	void EndJob();
	void Loop();

private:
	AnaClass *fAnaClass;
	TreeCleaner *fTreeCleaner;
	PhysQCAnalysis *fPhysQCAnalysis;
	MultiplicityAnalysis *fMultiplicityAnalysis;

};
#endif
