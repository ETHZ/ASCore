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
	const char* fSetofCuts;  
	float fLumi; 	// lumi to scale plots with xsection


private:
	void ReadCuts();
	
	// Multiplicity Plots Variables:
	TFile *fMPHistFile;
	LeptJetStat *fMyLeptJetStat;
	TH2D *fHljMult;
	TH2D *fHemuMult;
	TH1F *fHemuEff;
	
	//  ---- set of cuts ---
	TString fSetName;
	// ---- Electrons
	float fCut_ElPt;
	float fCut_ElEta;
	float fCut_ElRelIso04;
	int   fCut_ElIDRobustTight;  
	// ---- Muons
	float fCut_MuPt;
	float fCut_MuEta;
	float fCut_MuRelIso03;
	float fCut_MuNTkHits;
	int   fCut_MuTrackerMu;
	// ---- Jets
	float fCut_JPt;
	float fCut_JEta;



};
#endif
