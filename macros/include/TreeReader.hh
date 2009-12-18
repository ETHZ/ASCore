#ifndef TreeReader_hh
#define TreeReader_hh


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
	
	void ReadObjCuts(const char* = "objsel.dat");
	void ReadEvtSel(const char* = "evtsel.dat");
	
private:
	
	// Functions performing the cleaning and duplicate tagging 
	virtual void TagCleanObjects(void);
	virtual int CleanPrimaryVertex(void);
	virtual int IsFromPrimaryVx(int ipart, int ichk);
	virtual int CleanMuon(int ichk);
	virtual bool DuplicateMuon(int ichk);
	virtual int CleanElectron(int ichk);
	virtual bool DuplicateElectron(int ichk);
	virtual int CleanPhoton(int ichk);
	virtual int CleanJet(int ichk);
	virtual bool ElectronJet(int ichk);
	virtual int FindNearestJet(double eta, double phi);
	virtual double DeltaPhi(double v1, double v2);
	virtual double GetDeltaR(double eta1, double eta2, double phi1, double phi2);
	
	// Functions to actually perform the cleaning
	virtual void DecideIso(void);
	virtual void InitCleaning(int flag);
	virtual void DoCleanObjects(void);
	virtual void AddToJet(int ipart, int ichk, int iJet);
	virtual void SubtrFromJet(int ipart, int ichk, int iJet);
	virtual int CleanEvent(void);
	virtual int CleanMET(double met, double metphi);
	virtual int NextMuClean(void);
	virtual int NextElClean(void);
	virtual int NextJClean(void);
	
	// Global parameters:
	TString fOutputDir;
	TString fTag;
	bool fClean;
	bool fDiLep;
	bool fMPHist;
	bool fSignHist;
	
	TStyle *fStyle;
	Davismt2 *fMT2;
	TLatex *fTlat;

	// Cleaning variables
	int fDoClean;
	int fEvtClean;
	int fNMuClean;
	int fMuClean[20];
	int fNElClean;
	int fElClean[20];
	int fNJClean;
	int fJClean[50];
	int iMuNext;
	int iElNext;
	int iJNext;

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
	
	// Object quality cuts:
	struct Cut{
		TBranch *branch;
		double upperbound;
		double lowerbound;
	};
	bool IsGoodObj(int, std::vector<Cut>*);
	bool IsGoodEvt(std::vector<Cut>*);
	std::vector<Cut> fMuCuts;
	std::vector<Cut> fElCuts;
	std::vector<Cut> fJetCuts;
	std::vector<Cut> fEvtSelCuts;

};
#endif
