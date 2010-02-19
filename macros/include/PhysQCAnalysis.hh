#ifndef PhysQCAnalysis_hh
#define PhysQCAnalysis_hh


#include <vector>
#include <cmath>
#include <string>
#include <iostream>
#include <fstream>

#include <TH2D.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TStyle.h>
#include <TLatex.h>
#include <TLeaf.h>
#include <TBranch.h>

#include "TreeReader.hh"
#include "AnaClass.hh"
#include "UserAnalysisBase.hh"

class PhysQCAnalysis : public UserAnalysisBase{
public:
	PhysQCAnalysis(TreeReader *tr = 0);
	virtual ~PhysQCAnalysis();

	void Begin(const char* filename = "PhysQC.root");
	void Analyze();
	void End();
	void MakePlots(TString, TTree*);
	void PlotTriggerStats();
	
	TreeReader *fTR;
	AnaClass *fAC;

private:

};
#endif
