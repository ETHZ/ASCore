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

#include "TreeClassBase.h"
#include "LeptJetStat.h"
#include "Davismt2.h"

class TreeReader : public TreeClassBase{
public:
	TreeReader(TTree *tree=0, int flag = 111);
	virtual ~TreeReader();
	void BeginJob();
	void EndJob();
	void Loop();

	void setOutputDir(TString dir);

	void BookSignHists(const char* filename = "SignificancePlots.root");
	void FillSignHists(Int_t part);
	void WriteSignHists();

	void BookMPHistos(const char* filename = "MultiplicityPlots.root");
	void FillMPHistos();
	void PrintMPOutput();
	void PlotMPSummary(int it);
	void PlotMPEffic(int it);

	void InitDiLepTree(const char *filename = "DiLepTree.root");
	void FillDiLepTree();
	void ResetDiLepTree();
	void WriteDiLepTree();

	double getEta(double, double, double);
	
private:
	// Global parameters:
	TString fOutputDir;
	bool fDiLep;
	bool fMPHist;
	bool fSignHist;
	
	Davismt2 *fMT2;

	// Significance Plots:
	TFile *fSignHistsFile;
	int fNBinsEta;
	int fNBinsPhi;
	TH2D *fH_ptdev;
	TH2D *fH_ptsum;
	TH2D *fH_pt2sum;
	TH2I *fH_ptevt;
	TH2D *fH_ptavg;
	TH1D *fH_ptsumeta;
	TH1I *fH_ptevteta;

	// Multiplicity Plots Variables:
	TFile *fMPHistFile;
	std::vector<LeptJetStat *> fMyLeptJetStat;
	std::vector<TH2D *> fMyhljMult;
	std::vector<TH2D *> fMyhemuMult;
	std::vector<TH1F *> fMyhemuEff;
	int fNcuts;

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
