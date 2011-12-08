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
	void BeginJob	( std::string, std::string );
	void EndJob		();
	void Loop		(Int_t prescale=1);
	inline void SetPtHatCut(float cut){fPtHatCut = cut;};
	inline void DoFillEffTree(bool fill){fDoFillEffTree = fill;};

private:
	SSDLAnalysis	*fSSDLAnalysis;
	float fPtHatCut;
	bool fDoFillEffTree;
};
#endif
