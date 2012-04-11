/*****************************************************************************
*   Collection of tools for producing plots for same-sign dilepton analysis  *
*                                                                            *
*                                                  (c) 2010 Benjamin Stieger *
*****************************************************************************/
#ifndef SSDLPlotter_HH
#define SSDLPlotter_HH

#include "SSDLDumper.hh"
#include "helper/Monitor.hh"

#include "TLorentzVector.h"



struct SSDLPrediction {

	float bg_mm;
	float bg_em;
	float bg_ee;

	float bg_mm_err;
	float bg_em_err;
	float bg_ee_err;

	float s_mm;
	float s_em;
	float s_ee;

};


class SSDLPlotter : public SSDLDumper{

public:
	static TString gHiLoLabel[3];

	static double gEChMisIDB;
	static double gEChMisIDB_E;
	static double gEChMisIDE;
	static double gEChMisIDE_E;

	static float gMMTrigScale;
	static float gEMTrigScale;
	static float gEETrigScale;

	SSDLPlotter();
	SSDLPlotter(TString);
	SSDLPlotter(TString, TString);
	virtual ~SSDLPlotter();

	virtual void init(TString filename = "samples.dat");

	virtual void doAnalysis();
	virtual void sandBox();
	virtual void msugraKfacs(TFile *);
	virtual void msugraLOxsecs(TFile *);
	virtual void msugraNLOxsecs(TFile *);
	virtual void scanMSUGRA(const char * filestring);
	virtual void scanSMS(const char * filestring);
	virtual void plotWeightedHT();

	//////////////////////////////
	// Plots
	void makeMuIsolationPlots(bool = false);
	void makeElIsolationPlots(bool = false);
	void makeElIdPlots();
	
	void makeNT2KinPlots(gHiLoSwitch = HighPt);
	void makeMETvsHTPlot(vector<int>, vector<int>, vector<int>, gHiLoSwitch = HighPt);
	void makeMETvsHTPlotPRL();
	void makeMETvsHTPlot0HT();
	void makeMETvsHTPlotTau();
	
	void makeFRvsPtPlots(gChannel, gFPSwitch);
	void makeFRvsEtaPlots(gChannel);
	void makeFRvsPtPlotsForPAS(gChannel);
	void makeFRvsEtaPlotsForPAS(gChannel);
	void makeRatioPlots(gChannel);
	void make2DRatioPlots(gChannel);
	void makeNTightLoosePlots(gChannel);

	void makeOriginPlots(gRegion);
	
	void makeIsoVsMETPlot(gSample);
	void makePileUpPlots(bool write = true);
	
	void makePRLPlot1();
	
	void makeNT012Plots(vector<int>, gChannel, gRegion = Baseline, gHiLoSwitch = HighPt);
	void makeNT012Plots(gChannel, vector<int>, bool(SSDLPlotter::*)(int&, int&), TString = "");

	void makeAllIntPredictions();
	void makeIntPrediction(TString, gRegion, gHiLoSwitch = HighPt);
	void makeIntPredictionTTW(TString, gRegion, gHiLoSwitch = HighPt);
	SSDLPrediction makePredictionSignalEvents(float, float, float, float, int, int, int, float=20., float=10., bool=false);
	void makeDiffPrediction();

	void makeAllClosureTests();
	void makeIntMCClosure(vector<int>, TString, gRegion = Baseline, gHiLoSwitch = HighPt);

	void makeTTbarClosure();
	void makeRelIsoTTSigPlots();
	
	void storeWeightedPred();
	float getFRatio(gChannel, float, int = 0);        // diff in pt only
	float getFRatio(gChannel, float, float, int = 0); // diff in pt/eta
	float getPRatio(gChannel, float, int = 0);
	
	//////////////////////////////
	// Fake ratios
	// Calculate from pre stored numbers, with fixed selections:
	void fillRatios(vector<int>, vector<int>, int = 0);
	TH1D* fillRatioPt(gChannel, int, gFPSwitch, bool = false);
	TH1D* fillRatioPt(gChannel, vector<int>, gFPSwitch, bool = false);
	TH2D* fillRatio(gChannel, int, gFPSwitch, bool = false);
	TH2D* fillRatio(gChannel, vector<int>, gFPSwitch, bool = false);

	void calculateRatio(vector<int>, gChannel, gFPSwitch, TH2D*&, bool = false);
	void calculateRatio(vector<int>, gChannel, gFPSwitch, TH2D*&, TH1D*&, TH1D*&, bool = false);
	void calculateRatio(vector<int>, gChannel, gFPSwitch, float&, float&);
	void calculateRatio(vector<int>, gChannel, gFPSwitch, float&, float&, float&);
	
	void getPassedTotal(vector<int>,  gChannel, gFPSwitch, TH2D*&, TH2D*&, bool = false, gHiLoSwitch = HighPt);
	TH1D* getFRatio(vector<int>, gChannel, int = 0, bool = false);

	void ratioWithBinomErrors(float, float, float&, float&);
	void ratioWithPoissErrors(float, float, float&, float&);
	void ratioWithAsymmCPErrors(int, int, float&, float&, float&);

	void printYields(gChannel, float = -1.0);
	void printYieldsShort(float = -1);
	
	void printAllYieldTables();
	void printMCYieldTable(TString, gRegion = Baseline);

	TGraph* getSigEventGraph(gChannel, gRegion = Baseline);
	TGraph* getSigEventGraph(gChannel, float, float, float, float);

	void fixPRatios();
	
	// Geninfo stuff:
	int muIndexToBin(int);
	int elIndexToBin(int);
	TString muBinToLabel(int);
	TString elBinToLabel(int);
	void labelOriginAxis(TAxis*, gChannel);
	void label2OriginAxes(TAxis*, TAxis*, gChannel);
	
	void printOrigins(gRegion = Baseline);
	void printMuOriginTable(gRegion = Baseline);
	void printMuOriginHeader(TString);
	void printMuOriginFromSample(Sample*, int, gRegion = Baseline, gHiLoSwitch = HighPt);
	void print2MuOriginsFromSample(Sample*, int, gRegion = Baseline, gHiLoSwitch = HighPt);

	void printElOriginTable(gRegion = Baseline);
	void printElOriginHeader(TString);
	void printElOriginFromSample(Sample*, int, gRegion = Baseline, gHiLoSwitch = HighPt);
	void print2ElOriginsFromSample(Sample*, int, gRegion = Baseline, gHiLoSwitch = HighPt);

	void printEMuOriginTable(gRegion = Baseline);
	void printEMuOriginHeader(TString);
	void printEMuOriginsFromSample(Sample*, int, gRegion = Baseline, gHiLoSwitch = HighPt);

	void printOriginSummary(vector<int>, int, gChannel, gRegion = Baseline, gHiLoSwitch = HighPt);
	void printOriginSummary2L(vector<int>, int, gChannel, gRegion = Baseline, gHiLoSwitch = HighPt);
	
	virtual void drawTopLine(float = 0.70);
	virtual void drawRegionSel(gRegion);
	virtual void drawDiffCuts(int);
	
	const int     getNPtBins (gChannel);
	const double *getPtBins  (gChannel);
	const int     getNPt2Bins(gChannel);
	const double *getPt2Bins (gChannel);
	const int     getNEtaBins(gChannel);
	const double *getEtaBins (gChannel);
	
	vector<int> fMCBG;    // SM background MC samples
	vector<int> fMCBGNoQCDNoGJets;    // SM background MC samples without QCD or GJets
	vector<int> fMCBGNoQCDNoGJetsSig;    // SM background MC samples without QCD or GJets, with signal
	vector<int> fMCBGSig; // SM background + LM0 signal samples
	vector<int> fMCBGMuEnr;    // SM background MC samples with Muon enriched QCD
	vector<int> fMCBGMuEnrSig; // SM background + LM0 signal samples with Muon enriched QCD
	vector<int> fMCRareSM; // Rare SM backgrounds
	vector<int> fMuData;  // Muon data samples
	vector<int> fEGData;  // EG data samples
	vector<int> fMuEGData;  // MuEG dataset
	vector<int> fMuHadData;  // Muon data samples
	vector<int> fEleHadData;  // EG data samples
	vector<int> fHighPtData;  // All high pt triggered data
	vector<int> fLowPtData;   // All lepton cross HT triggered data

	TLatex *fLatex;
	
	float fLumiNorm;      // Normalize everything to this luminosity
	float fBinWidthScale; // Normalize bin contents to this width

	TH1D *fH1D_MufRatio;
	TH1D *fH1D_MupRatio;
	TH1D *fH1D_ElfRatio;
	TH1D *fH1D_ElpRatio;

	TH1D *fH1D_MufRatio_MC;
	TH1D *fH1D_MupRatio_MC;
	TH1D *fH1D_ElfRatio_MC;
	TH1D *fH1D_ElpRatio_MC;

	TH2D *fH2D_MufRatio;
	TH2D *fH2D_MupRatio;
	TH2D *fH2D_ElfRatio;
	TH2D *fH2D_ElpRatio;

	TH2D *fH2D_MufRatio_MC;
	TH2D *fH2D_MupRatio_MC;
	TH2D *fH2D_ElfRatio_MC;
	TH2D *fH2D_ElpRatio_MC;

	private:
	
};

#endif
