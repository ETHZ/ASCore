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


class SSDLPlotter : public SSDLDumper{

public:
	static TString gHiLoLabel[3];

	static double gEChMisIDB;
	static double gEChMisIDB_E;
	static double gEChMisIDE;
	static double gEChMisIDE_E;

	SSDLPlotter();
	SSDLPlotter(TString);
	SSDLPlotter(TString, TString);
	virtual ~SSDLPlotter();

	virtual void init(TString filename = "samples.dat");
	virtual void InitMC(TTree*); // remove a few branches
	virtual void readSamples(const char* filename = "samples.dat");

	virtual void doAnalysis();
	virtual void sandBox();
	virtual void load_kfacs(TFile *);
	virtual void load_loxsecs(TFile *);
	virtual void load_msugraInfo(const char * filestring);

	//////////////////////////////
	// Plots
	void makeMuIsolationPlots();
	void makeElIsolationPlots();
	void makeElIdPlots();
	
	void makeNT2KinPlots(gHiLoSwitch = HighPt);
	void makeMETvsHTPlot(vector<int>, vector<int>, vector<int>, gHiLoSwitch = HighPt);
	void makeMETvsHTPlotPRL();
	void makeMETvsHTPlotTau();
	
	void makeFRvsPtPlots(gChannel, gFPSwitch);
	void makeFRvsEtaPlots(gChannel);
	void makeFRvsPtPlotsForPAS(gChannel);
	void makeFRvsEtaPlotsForPAS(gChannel);
	void makeRatioPlots(gChannel);
	void makeNTightLoosePlots(gChannel);
	
	void makeIsoVsMETPlot(gSample);
	void makePileUpPlots(bool write = true);
	
	void makePRLPlot1();
	
	void makeMCClosurePlots(vector<int>);
	void makeDataClosurePlots();
	void makeNT012Plots(vector<int>, gChannel, gRegion = Baseline, gHiLoSwitch = HighPt);
	void makeNT012Plots(gChannel, vector<int>, bool(SSDLPlotter::*)(int&, int&), TString = "");

	void makeAllIntPredictions();
	void makeIntPrediction(TString, gRegion, gHiLoSwitch = HighPt);
	void makeDiffPrediction();
	void makeIntMCClosure(TString, gHiLoSwitch = HighPt);	
	void makeTTbarClosure();
	void makeRelIsoTTSigPlots();
	
	//////////////////////////////
	// Fake ratios
	// Produce from tree, with given selections:
	void produceRatio(gChannel, int, int, bool(SSDLPlotter::*)(), bool(SSDLPlotter::*)(int), TH2D*&, TH1D*&, TH1D*&, bool = false);
	void produceRatio(gChannel, vector<int>, int, bool(SSDLPlotter::*)(), bool(SSDLPlotter::*)(int), TH2D*&, TH1D*&, TH1D*&, bool = false);

	TH1D* fillMuRatioPt(int, int, bool(SSDLPlotter::*)(), bool(SSDLPlotter::*)(int), bool = false);
	TH1D* fillMuRatioPt(vector<int>, int, bool(SSDLPlotter::*)(), bool(SSDLPlotter::*)(int), bool = false);
	TH1D* fillMuRatioPt(vector<int>, int, bool(SSDLPlotter::*)(), bool(SSDLPlotter::*)(int), const int, const double*, const int, const double*, bool = false);

	// Calculate from pre stored numbers, with fixed selections:
	void fillMuElRatios(vector<int>);

	TH1D* fillMuRatioPt(int, gFPSwitch, bool = false);
	TH1D* fillMuRatioPt(vector<int>, gFPSwitch, bool = false);
	TH1D* fillElRatioPt(int, gFPSwitch, bool = false);
	TH1D* fillElRatioPt(vector<int>, gFPSwitch, bool = false);

	void calculateRatio(vector<int>, gChannel, gFPSwitch, TH2D*&, bool = false);
	void calculateRatio(vector<int>, gChannel, gFPSwitch, TH2D*&, TH1D*&, TH1D*&, bool = false);
	void calculateRatio(vector<int>, gChannel, gFPSwitch, float&, float&);
	void calculateRatio(vector<int>, gChannel, gFPSwitch, float&, float&, float&);
	
	void getPassedTotal(vector<int>,  gChannel, gFPSwitch, TH2D*&, TH2D*&, bool = false, gHiLoSwitch = HighPt);
	TH1D* getFRatio(vector<int>, gChannel, int = 0, bool = false);

	void ratioWithBinomErrors(float, float, float&, float&);
	void ratioWithPoissErrors(float, float, float&, float&);
	void ratioWithAsymmCPErrors(int, int, float&, float&, float&);

	void storeNumbers(Sample*, gChannel, gRegion);
	
	void printYields(gChannel, float = -1.0);
	void printYieldsShort(float = -1);


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
	
	virtual void drawTopLine();
	virtual void drawDiffCuts(int);
	
	const int     getNPtBins (gChannel);
	const double *getPtBins  (gChannel);
	const int     getNPt2Bins(gChannel);
	const double *getPt2Bins (gChannel);
	const int     getNEtaBins(gChannel);
	const double *getEtaBins (gChannel);
	
	vector<int> fMCBG;    // SM background MC samples
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

	private:
	
};

#endif
