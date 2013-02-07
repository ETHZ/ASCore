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

	float bg;
	float bg_err;

	float bg_mm;
	float bg_em;
	float bg_ee;

	float bg_mm_err;
	float bg_em_err;
	float bg_ee_err;

	float s_ttwz;
	float s_ttwz_mm;
	float s_ttwz_em;
	float s_ttwz_ee;

	float s_ttw;
	float s_ttw_mm;
	float s_ttw_em;
	float s_ttw_ee;

	float s_ttz;
	float s_ttz_mm;
	float s_ttz_em;
	float s_ttz_ee;

	int ns_ttw_mm;
	int ns_ttw_em;
	int ns_ttw_ee;

	int ns_ttz_mm;
	int ns_ttz_em;
	int ns_ttz_ee;

	float s_mm;
	float s_em;
	float s_ee;

	int obs_mm;
	int obs_em;
	int obs_ee;
	int obs;

};

struct TTWZPrediction {
	int obs;
	int obs_mm;
	int obs_ee;
	int obs_em;
	
	float ttw;
	float ttw_mm;
	float ttw_ee;
	float ttw_em;
	
	float ttz;
	float ttz_mm;
	float ttz_ee;
	float ttz_em;
	
	float ttwz;
	float ttwz_mm;
	float ttwz_ee;
	float ttwz_em;
	
	float fake;
	float fake_mm;
	float fake_ee;
	float fake_em;
	float fake_err;
	float fake_err_mm;
	float fake_err_ee;
	float fake_err_em;
	
	float cmid;
	float cmid_ee;
	float cmid_em;
	float cmid_err;
	float cmid_err_ee;
	float cmid_err_em;
	
	float wz;
	float wz_mm;
	float wz_ee;
	float wz_em;
	float wz_err;
	float wz_err_mm;
	float wz_err_ee;
	float wz_err_em;
	
	float rare;
	float rare_mm;
	float rare_ee;
	float rare_em;
	float rare_err;
	float rare_err_mm;
	float rare_err_ee;
	float rare_err_em;
};
struct SSPrediction {
	int obs;
        int obs_mm;
        int obs_ee;
        int obs_em;

        float s;         
        float s_mm; 
        float s_ee; 
        float s_em; 
  
        float bkg;
        float bkg_err;
        float bkg_stat;                  

	float fake;
	float fake_mm;
	float fake_ee;
	float fake_em;
	float fake_err;
	float fake_err_mm;
	float fake_err_ee;
	float fake_err_em;
	float fake_stat;
	float fake_stat_mm;
	float fake_stat_ee;
	float fake_stat_em;

	float sf;
	float sf_mm;
	float sf_ee;
	float sf_em;
	float sf_err;
	float sf_err_mm;
	float sf_err_ee;
	float sf_err_em;
	float sf_stat;
	float sf_stat_mm;
	float sf_stat_ee;
	float sf_stat_em;

        float df;
	float df_mm;
	float df_ee;
	float df_em;
	float df_err;
	float df_err_mm;
	float df_err_ee;
	float df_err_em;
	float df_stat;
	float df_stat_mm;
	float df_stat_ee;
	float df_stat_em;

	float cmid;
	float cmid_mm;
	float cmid_ee;
	float cmid_em;
	float cmid_err;
	float cmid_err_mm;
	float cmid_err_ee;
	float cmid_err_em;
	float cmid_stat;
	float cmid_stat_mm;
	float cmid_stat_ee;
	float cmid_stat_em;

	float wz;
	float wz_mm;
	float wz_ee;
	float wz_em;
	float wz_err;
	float wz_err_mm;
	float wz_err_ee;
	float wz_err_em;
	float wz_stat;
	float wz_stat_mm;
	float wz_stat_ee;
	float wz_stat_em;

	float rare;
	float rare_mm;
	float rare_ee;
	float rare_em;
	float rare_err;
	float rare_err_mm;
	float rare_err_ee;
	float rare_err_em;
	float rare_stat;
	float rare_stat_mm;
	float rare_stat_ee;
	float rare_stat_em;
};


class SSDLPlotter : public SSDLDumper{

public:
	static TString gHiLoLabel[3];

	static double gEChMisIDE;
	static double gEChMisIDE_E;
	static double gEChMisIDB;
	static double gEChMisIDB_E;

	static float gMMTrigScale;
	static float gEMTrigScale;
	static float gEETrigScale;

	bool fDO_OPT;

	SSDLPlotter(TString);
	SSDLPlotter(TString, TString);
	virtual ~SSDLPlotter();

	virtual void init(TString filename = "samples.dat");
        virtual void doSMSscans(TString);
	virtual void doAnalysis();
	virtual void sandBox();
	virtual void pythiaMadgraph(bool pyhia);
	virtual void makeROCCurve();
	virtual void msugraKfacs(TFile *);
	virtual void msugraLOxsecs(TFile *);
	virtual void msugraNLOxsecs(TFile *);
	virtual void scanMSUGRA(const char * filestring);
	virtual void scanSMS( const char * filestring, int);
	// old virtual void scanSMS( const char * filestring, gRegion);
	// virtual void scanSMSEWKpaper( const char * filestring, float minHT, float maxHT, float minMET, float maxMET, float pt1, float pt2, bool pass3rdVeto);
	virtual void plotWeightedHT();

	//////////////////////////////
	// Plots
	void makeMuIsolationPlots(bool = false);
	void makeElIsolationPlots(bool = false);
	void makeElIdPlots();
	
	void makeNT2KinPlots(bool=false);
	void makeMETvsHTPlot(vector<int>, vector<int>, vector<int>, gHiLoSwitch = HighPt);
	void makeMETvsHTPlotPRL();
	void makeMETvsHTPlot0HT();
	void makeMETvsHTPlotTau();
	
	void makeFRvsPtPlots(gChannel, gFPSwitch);
	void makeFRvsNVPlots(gChannel, gFPSwitch);
	void makeFRvsEtaPlots(gChannel);
	void makeRatioPlots(gChannel);
	void make2DRatioPlots(gChannel);
	void makeNTightLoosePlots(gChannel);
	void makeFakeGenIDTables();
	void makeFakeGenIDTable(vector<int>, gChannel, gFPSwitch, TString);

	void makeOriginPlots(int);
	// old void makeOriginPlots(gRegion);
	
	void makeIsoVsMETPlot(gSample);
	void makePileUpPlots(bool write = true);
	
	void makePRLPlot1();
	
	void makeNT012Plots(vector<int>, gChannel, int region = 0, gHiLoSwitch = HighPt);
	// old void makeNT012Plots(vector<int>, gChannel, gRegion = Baseline, gHiLoSwitch = HighPt);
	void makeNT012Plots(gChannel, vector<int>, bool(SSDLPlotter::*)(int&, int&), TString = "");

	void makeAllIntPredictions();
        void makeWZValidation(SSPrediction);
        SSPrediction makeIntPrediction(TString, int);
        // old SSPrediction makeIntPrediction(TString, gRegion);
	void makeTTWIntPredictions();
	TTWZPrediction makeIntPredictionTTW(TString, int);
	// old TTWZPrediction makeIntPredictionTTW(TString, gRegion);
	void makeSystPlot(TString outputname, TString label, TH1D *nom, TH1D *plus, TH1D *minus=NULL);

	void makeTTWIntPredictionsSigEvent(float minHT, float maxHT, float minMET, float maxMET, int minNjets, int minNbjetsL, int minNbjetsM, float pT1=20., float pT2=10., int chVeto = 0, bool ttw    =false, int flag=0);
	TTWZPrediction makePredictionSignalEvents(float minHT, float maxHT, float minMET, float maxMET, int minNjets, int minNbjetsL, int minNbjetsM, float pT1=20., float pT2=10., int chVeto = 0, bool ttw=false, int flag=0);
	void makeDiffPrediction();
	void makeTTWDiffPredictions();
	void makeDiffPredictionTTW(int);

        void makeAllClosureTests();
	void makeIntMCClosure(vector<int>, TString, int region = 0);
	void makeAllClosureTestsTTW();
	void makeIntMCClosureTTW(vector<int> samples, TString filename, int reg);
	// void makeIntMCClosure(vector<int>, TString, gRegion = Baseline);

	void makeTTbarClosure();
        void makeRelIsoTTSigPlots();
	
	void storeWeightedPred(int reg);
	// void storeWeightedPred(gRegion reg);
	float getFRatio(gChannel, float, int = 0);        // diff in pt only
	float getFRatio(gChannel, float, float, int = 0); // diff in pt/eta
	float getPRatio(gChannel, float, int = 0);
	
	//////////////////////////////
	// Fake ratios
	// Calculate from pre stored numbers, with fixed selections:
	void fillRatios(vector<int>, vector<int>, int = 0);
	void fillRatios(vector<int>, vector<int>, vector<int>, vector<int>, int = 0);
	TH1D* fillRatioPt(gChannel, int, gFPSwitch, bool = false);
	TH1D* fillRatioPt(gChannel, vector<int>, gFPSwitch, bool = false);
	TH2D* fillRatio(gChannel, int, gFPSwitch, bool = false);
	TH2D* fillRatio(gChannel, vector<int>, gFPSwitch, bool = false);

	void calculateRatio(vector<int>, gChannel, gFPSwitch, TH2D*&, bool = false);
	void calculateRatio(vector<int>, gChannel, gFPSwitch, TH2D*&, TH1D*&, TH1D*&, bool = false);
	void calculateRatio(vector<int>, gChannel, gFPSwitch, TH2D*&, TH1D*&, TH1D*&, TH1D*&, bool = false);
	void calculateRatio(vector<int>, gChannel, gFPSwitch, float&, float&);
	void calculateRatio(vector<int>, gChannel, gFPSwitch, float&, float&, float&);
        
	void calculateChMisIdProb(vector<int>, gChMisIdReg, float&, float&);

        TEfficiency *getMergedEfficiency(vector<int> samples, gChannel chan, gFPSwitch fp, int pteta=0);
	TGraphAsymmErrors *getCombEfficiency(vector<int> samples, gChannel chan, gFPSwitch fp, int pteta=0);
	
	void getPassedTotal(vector<int>,  gChannel, gFPSwitch, TH2D*&, TH2D*&, bool = false);
	void getPassedTotal(vector<int>,  gChannel, gFPSwitch, TH2D*&, TH2D*&, TH1D*&, TH1D*&, bool = false);
	void getPassedTotalTTbar(vector<int>,  gChannel, gFPSwitch, TH2D*&, TH2D*&, bool = false);
	TH1D* getFRatio(vector<int>, gChannel, int = 0, bool = false);

	void ratioWithBinomErrors(float, float, float&, float&);
	void ratioWithPoissErrors(float, float, float&, float&);
	void ratioWithAsymmCPErrors(int, int, float&, float&, float&);

	void printYields(gChannel, float = -1.0);
	void printYieldsShort(float = -1);
	
	void printAllYieldTables();
	void printMCYieldTable(TString, int region = 0);
	// old void printMCYieldTable(TString, gRegion = Baseline);

	TGraph* getSigEventGraph(gChannel, int region = 0);
	// old TGraph* getSigEventGraph(gChannel, gRegion = Baseline);
	// TGraph* getSigEventGraph(gChannel, float, float, float, float);
	TGraph* getSigEventGraph(gChannel, float, float, float, float, int = 0);

	void fixPRatios();
	
	// Geninfo stuff:
	int muIndexToBin(int);
	int elIndexToBin(int);
	TString muBinToLabel(int);
	TString elBinToLabel(int);
	void labelOriginAxis(TAxis*, gChannel);
	void label2OriginAxes(TAxis*, TAxis*, gChannel);
	
	void printOrigins(int region = 0);
	void printMuOriginTable(int region = 0);
	// old void printOrigins(gRegion = Baseline);
	// old void printMuOriginTable(gRegion = Baseline);
	void printMuOriginHeader(TString);
	void printMuOriginFromSample(Sample*, int,  int region = 0, gHiLoSwitch = HighPt);
	void print2MuOriginsFromSample(Sample*, int,  int region = 0, gHiLoSwitch = HighPt);
	// old void printMuOriginFromSample(Sample*, int, gRegion = Baseline, gHiLoSwitch = HighPt);
	// old void print2MuOriginsFromSample(Sample*, int, gRegion = Baseline, gHiLoSwitch = HighPt);

	void printElOriginTable( int region = 0);
	// old void printElOriginTable(gRegion = Baseline);
	void printElOriginHeader(TString);
	void printElOriginFromSample(Sample*, int,  int region = 0, gHiLoSwitch = HighPt);
	void print2ElOriginsFromSample(Sample*, int,  int region = 0, gHiLoSwitch = HighPt);
	// old void printElOriginFromSample(Sample*, int, gRegion = Baseline, gHiLoSwitch = HighPt);
	// old void print2ElOriginsFromSample(Sample*, int, gRegion = Baseline, gHiLoSwitch = HighPt);

	void printEMuOriginTable( int region = 0);
	// old void printEMuOriginTable(gRegion = Baseline);
	void printEMuOriginHeader(TString);
	void printEMuOriginsFromSample(Sample*, int,  int region = 0, gHiLoSwitch = HighPt);
	// old void printEMuOriginsFromSample(Sample*, int, gRegion = Baseline, gHiLoSwitch = HighPt);

	void printOriginSummary(vector<int>, int, gChannel,  int region = 0, gHiLoSwitch = HighPt);
	void printOriginSummary2L(vector<int>, int, gChannel,  int region = 0, gHiLoSwitch = HighPt);
	// old void printOriginSummary(vector<int>, int, gChannel, gRegion = Baseline, gHiLoSwitch = HighPt);
	// old void printOriginSummary2L(vector<int>, int, gChannel, gRegion = Baseline, gHiLoSwitch = HighPt);
	
	virtual void drawTopLine(float = 0.60, float = 1.0, float = 0.13);
	virtual void drawTopLineSim(float = 0.60, float = 1.0, float = 0.13);
	virtual void drawRegionSel(int);
	// old virtual void drawRegionSel(gRegion);
	virtual void drawDiffCuts(int);
	
	vector<int> fMCBG;    // SM background MC samples
	vector<int> fMCBGNoQCDNoGJets;    // SM background MC samples without QCD or GJets
	vector<int> fMCBGNoQCDNoGJetsSig;    // SM background MC samples without QCD or GJets, with signal
	vector<int> fMCBGSig; // SM background + LM0 signal samples
	vector<int> fMCBGMuEnr;    // SM background MC samples with Muon enriched QCD
	vector<int> fMCBGEMEnr;    // SM background MC samples with Elec enriched QCD
	vector<int> fMCBGMuEnrSig; // SM background + LM0 signal samples with Muon enriched QCD
	vector<int> fMCRareSM; // Rare SM backgrounds
	vector<int> fMuData;  // Muon data samples
	vector<int> fEGData;  // EG data samples
	vector<int> fMuEGData;  // MuEG dataset
	vector<int> fMuHadData;  // Muon data samples
	vector<int> fEleHadData;  // EG data samples
	vector<int> fHighPtData;  // All high pt triggered data
	vector<int> fLowPtData;   // All lepton cross HT triggered data
	vector<int> fClosureSamples;
	vector<int> fTTJets;
	vector<int> fWJets;
	vector<int> fMuEnr;
	vector<int> fEMEnr;
	vector<int> fDYJets;

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
