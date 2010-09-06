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

	void PlotMPSummary(LeptJetStat *fLeptJetStat, TH2D *fHljMult);
	void PlotMPEffic(LeptJetStat *fLeptJetStat, TH2D *fHemuMult, TH1F *fHemuEff);
	const char* fSetofCuts;  
	float fLumi; 	// lumi to scale plots with xsection
	std::vector<std::string>* fRequiredHLT; 
	std::vector<std::string>* fVetoedHLT;
	

private:
	void ReadCuts();
	void Reset();
	
	std::vector<int> elecs;
	std::vector<int> muons;
	
	
	// Multiplicity Plots Variables:
	TFile *fMPHistFile;
	LeptJetStat *fMyLeptJetStat;
	LeptJetStat *fMyLeptJetStat_bjets;
	TH2D *fHljMult_alljets;
	TH2D *fHemuMult_alljets;
	TH1F *fHemuEff;
	TH2D *fHljMult_bjets;
	TH2D *fHemuMult_bjets;
	
	//  ---- set of cuts ---
	TString fSetName;
	float fCut_PFMET;
	
};
#endif
