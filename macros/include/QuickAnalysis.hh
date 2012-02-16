#ifndef QuickAnalysis_hh
#define QuickAnalysis_hh


#include <cmath>
#include <string>
#include <iostream>
#include <fstream>

#include <TH1D.h>

#include "base/TreeReader.hh"
#include "base/UserAnalysisBase.hh"

class QuickAnalysis : public UserAnalysisBase{
public:
	QuickAnalysis(TreeReader *tr = NULL);
	virtual ~QuickAnalysis();

	void Begin(const char* filename = "QuickHistos.root");
	void Analyze();
	void End();

private:
	
	// file for histograms:
	TFile *fHistFile;
	// histo	
	TH1D *fHpileup;
};
#endif
