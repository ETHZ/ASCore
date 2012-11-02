#ifndef TTWAcceptanceAnalyzer_hh
#define TTWAcceptanceAnalyzer_hh

#include <TTree.h>
#include <TStyle.h>

#include "base/TreeAnalyzerBase.hh"
#include "base/TreeReader.hh"
#include "TTWAcceptanceAnalysis.hh"

class TTWAcceptanceAnalyzer : public TreeAnalyzerBase {
public:
	TTWAcceptanceAnalyzer(TTree *tree = 0);
	virtual ~TTWAcceptanceAnalyzer();
	void BeginJob();
	void EndJob();
	void Loop();

private:
	TTWAcceptanceAnalysis *fTTWAcceptanceAnalysis;

};
#endif
