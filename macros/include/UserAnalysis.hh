#ifndef UserAnalysis_hh
#define UserAnalysis_hh


#include <cmath>
#include <string>
#include <iostream>
#include <fstream>

#include <TH1D.h>

#include "base/TreeReader.hh"
#include "base/UserAnalysisBase.hh"

class UserAnalysis : public UserAnalysisBase{
public:
	UserAnalysis(TreeReader *tr = NULL);
	virtual ~UserAnalysis();

	void Begin();
	void Analyze();
	void End();

private:
	
	// file for histograms:
	TFile *fHistFile;
	
	TH1D *fHjetMult;
	TH1D *fHjetPt;
	TH1D *fHmuMult;
	TH1D *fHmuPt;
};
#endif
