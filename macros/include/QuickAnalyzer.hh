#ifndef QuickAnalyzer_hh
#define QuickAnalyzer_hh

#include <TTree.h>
#include <TStyle.h>

#include "base/TreeAnalyzerBase.hh"
#include "base/TreeReader.hh"
#include "QuickAnalysis.hh"

class QuickAnalyzer : public TreeAnalyzerBase {
public:
	QuickAnalyzer(std::vector<std::string>& fileList);
	virtual ~QuickAnalyzer();
	void BeginJob(TString filename="QuickHistos.root");
	void EndJob();
	void Loop();

private:
	QuickAnalysis *fQuickAnalysis;

};
#endif
