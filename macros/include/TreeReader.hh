#ifndef TreeReader_hh
#define TreeReader_hh


#include <vector>
#include <cmath>
#include <string>
#include <iostream>

#include <TH2D.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TStyle.h>
#include <TLatex.h>

#include "TreeClassBase.h"
#include "LeptJetStat.h"
#include "Davismt2.h"
// #include "ETHStyle.h"

class TreeReader : public TreeClassBase{
public:
	TreeReader(TTree *tree=0, int flag = 111);
	virtual ~TreeReader();
	void DefStyle();
	void BeginJob();
	void EndJob();
	void Loop();

	void SetOutputDir(TString dir);
	inline void SetTag(TString tag){fTag = tag;};

	void BookSignHists(const char* filename = "SignificancePlots.root");
	void FillSignHists(Int_t part);
	void WriteSignHists();

	void BookMPHistos(const char* filename = "MultiplicityPlots.root");
	void FillMPHistos();
	void PrintMPOutput();
	void PlotMPSummary();
	void PlotMPEffic();

	void InitDiLepTree(const char *filename = "DiLepTree.root");
	void FillDiLepTree();
	void ResetDiLepTree();
	void WriteDiLepTree();

	void PlotMultiplicity();

	void printPNG(TCanvas*, TString, TString);
	void printEPS(TCanvas*, TString, TString);
	
	double getEta(double, double, double);
	
private:
	// Global parameters:
	TString fOutputDir;
	TString fTag;
	bool fDiLep;
	bool fMPHist;
	bool fSignHist;
	
	TStyle *fStyle;
	Davismt2 *fMT2;
	TLatex *fTlat;

	// Significance Plots:
	TFile *fSignHistsFile;
	int fNBinsEta[5];
	int fNBinsPhi;
	TH2D *fH_ptdev[5];
	TH2D *fH_ptsum[5];
	TH2D *fH_pt2sum[5];
	TH2I *fH_ptevt[5];
	TH2D *fH_ptavg[5];
	TH1D *fH_ptsumeta[5];
	TH1I *fH_ptevteta[5];

	// Multiplicity Plots Variables:
	TFile *fMPHistFile;
	LeptJetStat *fMyLeptJetStat;
	TH2D *fHljMult;
	TH2D *fHemuMult;
	TH1F *fHemuEff;

	// DiLepton Tree Variables:
	TFile *fDiLepTreeFile;
	TTree *fDiLepTree;
	int    fTRunNumber;
	int    fTEventNumber;
	int    fTLumiSection;

	int fTMu1charge;
	int fTMu2charge;
	int fTNqualmu;
	double fTMu1pt;
	double fTMu2pt;
	double fTMu1eta;
	double fTMu2eta;
	double fTMu1iso;
	double fTMu2iso;
	double fTMu1d0;
	double fTMu2d0;
	double fTMu1ntkhits;
	double fTMu2ntkhits;
	double fTMuminv;
	double fTMumt2_50;
	double fTMumt2_100;

	int fTNqualel;
	int fTEl1charge;
	int fTEl2charge;
	double fTEl1pt;
	double fTEl2pt;
	double fTEl1eta;
	double fTEl2eta;
	double fTEl1iso;
	double fTEl2iso;
	double fTEl1d0;
	double fTEl2d0;
	double fTElminv;
	double fTElmt2_50;
	double fTElmt2_100;
	
};
#endif
