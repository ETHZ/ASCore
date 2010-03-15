#ifndef TreeSkimmer_hh
#define TreeSkimmer_hh

#include <TTree.h>
#include <TFile.h>

#include "base/TreeAnalyzerBase.hh"
#include "base/TreeReader.hh"

class TreeSkimmer : public TreeAnalyzerBase {
public:
	TreeSkimmer(TTree *tree = NULL);
	virtual ~TreeSkimmer();
	void BeginJob();
	void EndJob();
	void Loop();

	bool EventSelection();

private:
	unsigned int fNsel;
	unsigned int fNtot;

	TTree *fTree;
	TTree *fSkimmedTree;
	TFile *fSkimmedTreeFile;
};
#endif
