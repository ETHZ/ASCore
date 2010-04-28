#ifndef PhysQCAnalysis_hh
#define PhysQCAnalysis_hh


#include <vector>
#include <cmath>
#include <string>
#include <iostream>
#include <fstream>
#include <time.h>

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
#include "helper/AnaClass.hh"
#include "base/UserAnalysisBase.hh"
#include "TreeCleaner.hh"

class PhysQCAnalysis : public UserAnalysisBase{
public:
	PhysQCAnalysis(TreeReader *tr = NULL, TreeCleaner *tc = NULL);
	virtual ~PhysQCAnalysis();

	void Begin();
	void Analyze1();
	void Analyze2();
	void End();
	void MakePlots(TString, TCut, TTree*);
	void MakeElIDPlots(TCut, TTree*);
	void PlotTriggerStats();
	void PrintHisto(TH1D*, TString, TString, TString, 
		double x1=999., double x2=999., int logy=0);
	void PrintInfoStart(int nEntries);
	double invMass(double p1[], double p2[]);
	double DeltaPhi(double v1, double v2);
	void GetEvtEmChFrac(double & fracEm, double & fracCh);
	
	TreeCleaner *fTC;
	AnaClass *fAC;

private:
	
	TH1D *fMuHistos[7];
	TH1D *fElHistos[20];
	TH1D *fPhHistos[6];
	TH1D *fJHistos[4];
	TH1D *fMETHistos[7];
	TH2D *fMETDphi12;
	TH2D *fMETR12R21;
	TH2D *fMETR12Dphij12;
	TH1D *fMuCIHistos[4];
	TH1D *fElCIHistos[5];
	TH1D *fPhCIHistos[4];
	TH1D *fJCIHistos[4];
	TH1D *fInvMHistos[10];
	
};
#endif
