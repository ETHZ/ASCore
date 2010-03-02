#ifndef MultiplicityAnalysis_hh
#define MultiplicityAnalysis_hh


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

#include "base/TreeReader.hh"
#include "base/UserAnalysisBase.hh"
#include "helper/LeptJetStat.h"

class MultiplicityAnalysis : public UserAnalysisBase{
public:
	MultiplicityAnalysis(TreeReader *tr = 0);
	virtual ~MultiplicityAnalysis();

	void Begin(const char* filename = "MultiplicityPlots.root");
	void Analyze();
	void End();

	void PlotMPSummary();
	void PlotMPEffic();

private:
	// Multiplicity Plots Variables:
	TFile *fMPHistFile;
	LeptJetStat *fMyLeptJetStat;
	TH2D *fHljMult;
	TH2D *fHemuMult;
	TH1F *fHemuEff;



};
#endif
