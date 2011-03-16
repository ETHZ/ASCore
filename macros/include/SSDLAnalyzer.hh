#ifndef SSDLAnalyzer_hh
#define SSDLAnalyzer_hh

#include <TTree.h>
#include <TStyle.h>

#include "base/TreeAnalyzerBase.hh"
#include "SSDLAnalysis.hh"

class SSDLAnalyzer : public TreeAnalyzerBase {
public:
	SSDLAnalyzer(TTree *tree = 0);
	virtual ~SSDLAnalyzer();
	void BeginJob	();
	void EndJob		();
	void Loop		(Int_t prescale=1);
	inline void SetPtHatCut(float cut){fPtHatCut = cut;};

private:
	SSDLAnalysis	*fSSDLAnalysis;
	float fPtHatCut;
};
#endif
