#ifndef UserAnalyzer_hh
#define UserAnalyzer_hh

#include <TTree.h>
#include <TStyle.h>

#include "base/TreeAnalyzerBase.hh"
#include "base/TreeReader.hh"
#include "UserAnalysis.hh"

class UserAnalyzer : public TreeAnalyzerBase {
public:
	UserAnalyzer(TTree *tree = 0);
	virtual ~UserAnalyzer();
	void BeginJob();
	void EndJob();
	void Loop();

private:
	UserAnalysis *fUserAnalysis;

};
#endif
