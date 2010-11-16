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

	void Begin(const char* filename = "PhysQCAnalysis.root" );
	void Analyze1();
	void Analyze2();
	void End();
	void BookHistos();
	void MakePlots(TString, TCut, TTree*);
	void MakeElIDPlots(TCut, TTree*);
	void PlotTriggerStats();
	void PrintHisto(TH1D*, TString, TString, TString, 
		double x1=999., double x2=999., int logy=0, int fract=0);
	void PrintInfoStart(int nEntries);
	
	TreeCleaner *fTC;
	AnaClass *fAC;

private:

	TFile* fHistFile; // Where all histograms will be saved

	TH1D *fPvxHistos[4];
	TH1D *fMuHistos[9];
	TH1D *fElHistos[20];
	TH1D *fPhHistos[7];
	TH1D *fJHistos[8];
	TH2D *fJChEMfrac;
	TH2D *fJEMfracEta;
	TH2D *fJChfracEta;
	TH1D *fMETHistos[10];
	TH2D *fMETDphi12;
	TH2D *fMETR12R21;
	TH2D *fMETR12Etaj12;
	TH2D *fMETR12j12PtRat;
	TH2D *fMETR12Dphij12;
	TH2D *fMETR12dRj12;
	TH1D *fEvtHistos[2];
	TH1D *fMuCIHistos[4];
	TH1D *fElCIHistos[5];
	TH1D *fPhCIHistos[4];
	TH1D *fJCIHistos[14];
	TH1D *fInvMHistos[12];
	
};
#endif
