#ifndef SSDLAnalyzer_hh
#define SSDLAnalyzer_hh

#include <TTree.h>
#include <TStyle.h>

#include "base/TreeAnalyzerBase.hh"
#include "base/TreeReader.hh"
#include "TreeCleaner.hh"
#include "SSDLAnalysis.hh"
#include "helper/AnaClass.hh"
#include "MultiplicityAnalysis.hh"

class SSDLAnalyzer : public TreeAnalyzerBase {
public:
	SSDLAnalyzer(TTree *tree = 0);
	virtual ~SSDLAnalyzer();
	void BeginJob();
	void EndJob();
	
	void Loop(Long64_t maxEvents=-1, Int_t prescale=1);

private:
	AnaClass				*fAnaClass;
	TreeCleaner				*fTreeCleaner;
	SSDLAnalysis			*fSSDLAnalysis;
	MultiplicityAnalysis	*fMultiplicityAnalysis;

};
#endif
