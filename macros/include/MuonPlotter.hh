/*****************************************************************************
*   Collection of tools for producing plots for Muon Fake Rate Analysis      *
*                                                                            *
*                                                  (c) 2010 Benjamin Stieger *
*****************************************************************************/
#ifndef MuonPlotter_HH
#define MuonPlotter_HH

#include "helper/AnaClass.hh"
#include "helper/FPRatios.hh"
#include "helper/Monitor.hh"
#include "helper/Davismt2.h"

#include "TLorentzVector.h"
#include "TEfficiency.h"
#include "TGraphAsymmErrors.h"


class MuonPlotter : public AnaClass{

public:
	// Binning
	static const int gNMuPtbins  = 5;
	static const int gNMuPt2bins = 5;
	static const int gNMuEtabins = 4;
	// static const int gNMuEtabins = 8;
	static const int gNElPtbins  = 5;
	static const int gNElPt2bins = 5;
	static const int gNElEtabins = 3;
	// static const int gNElEtabins = 6;
	
	static const int gNNVrtxBins = 9;
	static double gNVrtxBins[gNNVrtxBins+1];

	static const int gNJetPtBins = 10;
	static double gJetPtBins[gNJetPtBins+1];
	
	static double gMuPtbins [gNMuPtbins+1];
	static double gMuPt2bins[gNMuPt2bins+1];
	static double gMuEtabins[gNMuEtabins+1];

	static double gElPtbins [gNElPtbins+1];
	static double gElPt2bins[gNElPt2bins+1];
	static double gElEtabins[gNElEtabins+1];

	static double gEChMisIDB;
	static double gEChMisIDB_E;
	static double gEChMisIDE;
	static double gEChMisIDE_E;

	// This enum has to correspond to the content of the samples.dat file
	enum gSample {
		sample_begin,
		DoubleMu1 = sample_begin, DoubleMu2, DoubleEle1, DoubleEle2, MuEG1, MuEG2,
		MuHad1, MuHad2, EleHad1, EleHad2,
		TTJetsSync,
		TTJets, TJets_t, TJets_tW, TJets_s, WJets, DYJets, DYJets50, DYJets10to50,
		GJets40, GJets100, GJets200, GVJets, WW, WZ, ZZ, VVTo4L,
		LM0, LM1, LM2, LM3, LM4, LM5, LM6, LM7, LM8, LM9, LM11, LM12, LM13, 
		QCDMuEnr10,
		QCD5,
		QCD15,
		QCD30,
		QCD50,
		QCD80,
		QCD120,
		QCD170,
		QCD300,
		QCD470,
		QCD600,
		QCD800,
		QCD1000,
		QCD1400,
		QCD1800,
		QCD50MG,
		QCD100MG,
		QCD250MG,
		QCD500MG,
		QCD1000MG,
		gNSAMPLES
	};
	enum gFPSwitch{
		SigSup,
		ZDecay
	};
	enum gHiLoSwitch{
		HighPt,
		LowPt,
		TauChan
	};
	enum gRegion {
		region_begin,
		Baseline = region_begin,
		Control,
		HT80MET100,
		HT200MET30,
		HT200MET120,
		HT400MET50,
		HT400MET120,
		gNREGIONS
	};
	enum gChannel {
		channels_begin,
		Muon = channels_begin,
		EMu,
		Electron,
		gNCHANNELS
	};
	struct lepton{
		lepton(){};
		lepton(TLorentzVector vec, int ch, int ty, int ind){
			p = vec;
			charge = ch;
			type = ty;
			index = ind;
		};
		TLorentzVector p;
		int charge;
		int type; // -1(unknown), 0(mu), 1(ele)
		int index;
	};
	
	struct NumberSet{
		long nt2;
		long nt10;
		long nt01;
		long nt0;
		long nsst;
		long nssl;
		long nzl;
		long nzt;
	};

	struct Channel{ // like mumu, emu, ee
		TString name;
		TString sname;
		TH2D *nt20_pt; // pt1 vs pt2
		TH2D *nt10_pt;
		TH2D *nt01_pt; // only filled for e/mu
		TH2D *nt00_pt;
		TH2D *nt20_eta; // eta1 vs eta2
		TH2D *nt10_eta;
		TH2D *nt01_eta;
		TH2D *nt00_eta;

		TH2D *fntight; // pt vs eta
		TH2D *fnloose; 
		TH2D *pntight; // pt vs eta
		TH2D *pnloose;

		// Gen matched yields: t = tight, p = prompt, etc.
		TH2D *npp_pt; // overall pp/fp/.., binned in pt1 vs pt2
		TH2D *npp_cm_pt; // charge misid
		TH2D *nfp_pt;
		TH2D *npf_pt; // only filled for e/mu
		TH2D *nff_pt;
		TH2D *nt2pp_pt; // pp/fp/.. in tt window, binned in pt1 vs pt2
		TH2D *nt2pp_cm_pt; // charge misid
		TH2D *nt2fp_pt;
		TH2D *nt2pf_pt; // only filled for e/mu
		TH2D *nt2ff_pt;

		// Origin histos
		TH2D *nt11_origin;
		TH2D *nt10_origin;
		TH2D *nt01_origin;
		TH2D *nt00_origin;
		TH1D *sst_origin;
		TH1D *ssl_origin;
		TH1D *zt_origin;
		TH1D *zl_origin;

		// OS Yields
		// Only filled for electrons
		// For e/mu channel, use only BB and EE to count e's in barrel and endcaps
		TH2D *nt20_OS_BB_pt; // binned in pt1 vs pt2
		TH2D *nt20_OS_EE_pt; // BB = barrel/barrel, EE = endcap/endcap
		TH2D *nt20_OS_EB_pt; // EB = barrel/endcap
	};
	
	struct Region{
		static TString sname[gNREGIONS];
		// Two different pt cuts
		static float minMu1pt[2];
		static float minMu2pt[2];
		static float minEl1pt[2];
		static float minEl2pt[2];
		// Custom selections for every region
		static float minHT   [gNREGIONS];
		static float maxHT   [gNREGIONS];
		static float minMet  [gNREGIONS];
		static float maxMet  [gNREGIONS];
		static int   minNjets[gNREGIONS];
		
		Channel mm;
		Channel em;
		Channel ee;
	};
	
	// static const int gNRatioVars = 8;
	static const int gNRatioVars = 9;
	struct FRatioPlots{
		static TString var_name[gNRatioVars];
		static int nbins[gNRatioVars];
		static float xmin[gNRatioVars];
		static float xmax[gNRatioVars];
		TH1D *ntight[gNRatioVars];
		TH1D *nloose[gNRatioVars];
	};
	
	static const int gNKinVars = 10;
	struct KinPlots{
		static TString var_name[gNKinVars];
		static TString axis_label[gNKinVars];
		static int nbins[gNKinVars];
		static float xmin[gNKinVars];
		static float xmax[gNKinVars];
		TH1D *hvar[gNKinVars];
		TH2D *hmetvsht;
		static const int nMETBins = 100;
		static const int nHTBins = 200;
		static const float METmin = 0.;
		static const float METmax = 400.;
		static const float HTmin = 0.;
		static const float HTmax = 1000.;
	};
	
	static const int gNHWWVars = 10;
	struct HWWPlots{
		static TString var_name[gNHWWVars];
		static TString axis_label[gNHWWVars];
		static int nbins[gNHWWVars];
		static float xmin[gNHWWVars];
		static float xmax[gNHWWVars];
		TH1D *hvar[gNHWWVars];
	};
	
	static const int gNSels = 2;
	struct IsoPlots{
		static TString sel_name[gNSels];
		static int nbins[gNSels];
		TH1D *hiso[gNSels];
		TH1D *hiso_pt[gNSels][gNMuPt2bins];
		TH1D *hiso_nv[gNSels][gNNVrtxBins];
	};
	
	static const int gNHWWSels = 3;
	static TString gHWWSelNames[gNHWWSels];
	static const int gNKinSels = 3;
	static TString gKinSelNames[gNKinSels];
	static TString gEMULabel[2];
	static TString gHiLoLabel[3];

	class Sample{
	public:
		Sample(){};
		
		TString name;
		TString sname;
		TFile *file;
		TTree *tree;
		float lumi;
		int color;
		int datamc; // 0: Data, 1: SM MC, 2: Signal MC
		Region region[gNREGIONS][2];
		NumberSet numbers[gNREGIONS][gNCHANNELS]; // summary of integrated numbers
		KinPlots    kinplots[gNKinSels][2]; // tt and ll and signal for both low and high pt analysis
		HWWPlots    hwwplots[gNHWWSels]; // 0: no event sel, 1: N-1 sel
		IsoPlots    isoplots[2]; // e and mu
		FRatioPlots ratioplots[2]; // e and mu
		TGraph *sigevents[gNCHANNELS][2];
	};
	
	MuonPlotter();
	MuonPlotter(TString);
	MuonPlotter(TString, TString);
	virtual ~MuonPlotter();

	inline void setCharge(int charge){ fChargeSwitch = charge; };

	void init(TString filename = "samples.dat");
	void InitMC(TTree*); // remove a few branches
	void readSamples(const char* filename = "samples.dat");
	void loadSamples();
	void setBinning();

	void doAnalysis();
	void doLoop();
	
	void sandBox();
	
	//////////////////////////////
	// Plots
	void makeMufRatioPlots(bool = false, gHiLoSwitch = HighPt);
	void makeMupRatioPlots(bool = false, gHiLoSwitch = HighPt);
	void makeElfRatioPlots(bool = false, gHiLoSwitch = HighPt);
	void makeElpRatioPlots(bool = false, gHiLoSwitch = HighPt);

	void makeMuIsolationPlots();
	void makeElIsolationPlots();
	
	void makeNT2KinPlots(gHiLoSwitch = HighPt);
	void makeMETvsHTPlot(vector<int>, vector<int>, vector<int>, gHiLoSwitch = HighPt);
	void makeMETvsHTPlotCustom();
	void makeMETvsHTPlotTau();
	
	void makeFRvsPtPlots(gChannel, gFPSwitch);
	void makeFRvsEtaPlots(gChannel);
	void makeFRvsPtPlotsForPAS(gChannel, gFPSwitch);
	void makeFRvsEtaPlotsForPAS(gChannel);
	void makeRatioPlots(gChannel);
	
	void makeHWWPlots();
	
	void makeIsoVsMETPlot(gSample);
	
	void makeMCClosurePlots(vector<int>);
	void makeDataClosurePlots();
	void makeNT012Plots(vector<int>, gChannel, gRegion = Baseline, gHiLoSwitch = HighPt);
	void makeNT012Plots(gChannel, vector<int>, bool(MuonPlotter::*)(int&, int&), TString = "");

	void makeIntPrediction(TString, gRegion, gHiLoSwitch = HighPt);
	void makeIntPredictionOLD(TString, gRegion, gHiLoSwitch = HighPt);
	void makeIntMCClosure(TString, gHiLoSwitch = HighPt);
	
	void makeTTbarClosure();
	
	void printSyncExercise();
	
	//////////////////////////////
	// Fake ratios
	// Produce from tree, with given selections:
	void produceRatio(gChannel, int, int, bool(MuonPlotter::*)(), bool(MuonPlotter::*)(int), TH2D*&, TH1D*&, TH1D*&, bool = false);
	void produceRatio(gChannel, vector<int>, int, bool(MuonPlotter::*)(), bool(MuonPlotter::*)(int), TH2D*&, TH1D*&, TH1D*&, bool = false);

	TH1D* fillMuRatioPt(int, int, bool(MuonPlotter::*)(), bool(MuonPlotter::*)(int), bool = false);
	TH1D* fillMuRatioPt(vector<int>, int, bool(MuonPlotter::*)(), bool(MuonPlotter::*)(int), bool = false);
	TH1D* fillMuRatioPt(vector<int>, int, bool(MuonPlotter::*)(), bool(MuonPlotter::*)(int), const int, const double*, const int, const double*, bool = false);

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

	//////////////////////////////
	// Predictions
	void makeSSMuMuPredictionPlots(vector<int>, bool = false, gHiLoSwitch = HighPt);
	void makeSSElElPredictionPlots(vector<int>, bool = false, gHiLoSwitch = HighPt);
	void makeSSElMuPredictionPlots(vector<int>, bool = false, gHiLoSwitch = HighPt);
	void NObs(gChannel, TH1D *&, vector<int>, bool(MuonPlotter::*)());
	void NObs(gChannel, TH1D *&, vector<int>, gRegion = Baseline, gHiLoSwitch = HighPt);
	void NObs(gChannel, THStack *&, vector<int>, gRegion = Baseline, gHiLoSwitch = HighPt);
	vector<TH1D*> MuMuFPPrediction(TH2D* fratio, TH2D* pratio, TH2D* nt2, TH2D* nt1, TH2D* nt0, bool output = false);
	vector<TH1D*> ElElFPPrediction(TH2D* fratio, TH2D* pratio, TH2D* nt2, TH2D* nt1, TH2D* nt0, bool output = false);
	vector<TH1D*> ElMuFPPrediction(TH2D* mufratio, TH2D* mupratio, TH2D* elfratio, TH2D* elpratio,  TH2D* nt2, TH2D* nt10, TH2D* nt01, TH2D* nt0, bool output = false);
	
	void initCounters(gSample);
	void storeNumbers(Sample*, gChannel, gRegion);
	void printCutFlows(TString);
	void printCutFlow(gChannel, int, int);
	void printCutFlowsOld(TString);
	
	void printYields(gChannel, float = -1.0);
	void printYieldsShort(float = -1);

	//////////////////////////////
	// Fillers
	void fillYields(Sample*, gRegion, gHiLoSwitch);
	void fillRatioPlots(Sample*);
	void fillHWWPlots(Sample*);

	void fillMuIsoPlots(gSample);
	void fillElIsoPlots(gSample);
	void fillKinPlots(gSample, gHiLoSwitch);
	
	//////////////////////////////
	// I/O
	void bookHistos(Sample*);
	void deleteHistos(Sample*);
	void writeHistos(Sample*, TFile*);
	void writeSigGraphs(Sample*, gChannel, gHiLoSwitch, TFile*);
	int readHistos(TString);
	int readSigGraphs(TString);
	
	void bookRatioHistos();
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

	// Trigger selections:
	bool  singleMuTrigger();
	float singleMuPrescale();
	bool  singleElTrigger();
	float singleElPrescale();

	bool mumuSignalTrigger();
	bool elelSignalTrigger();
	bool elmuSignalTrigger();

	bool doubleMuTrigger();
	bool doubleElTrigger();
	bool eMuTrigger();

	bool doubleMuHTTrigger();
	bool doubleElHTTrigger();
	bool eMuHTTrigger();
	
	bool isGoodRun(gSample);

	// Event and Object selectors:
	int getNJets();
	int getNBTags();
	float getHT();
	float getMT2(int, int, int);
	float getClosestJetPt(int, gChannel);
	float getClosestJetDR(int, gChannel);
	float getSecondClosestJetDR(int, gChannel);
	float getAwayJetPt(int, gChannel);
	float getMaxJPt();
	
	int isSSLLEvent(int&, int&);
	int isOSLLEvent(int&, int&);
	int isSSEvent(int&, bool(MuonPlotter::*)(int), int&, bool(MuonPlotter::*)(int));
	int isOSEvent(int&, bool(MuonPlotter::*)(int), int&, bool(MuonPlotter::*)(int));
	vector<lepton> sortLeptonsByPt(vector<lepton> &leptons);
	vector<lepton> sortLeptonsByTypeAndPt(vector<lepton> &leptons);

	bool isGoodEvent();
	bool isGoodMuEvent();
	int hasLooseMuons(int&, int&);
	int hasLooseMuons();
	int hasLooseElectrons(int&, int&);
	int hasLooseElectrons();
	bool passesNJetCut(int=2);
	bool passesNJetCut_LooseLep(int=2);
	bool passesJet50Cut();
	
	bool passesHTCut(float, float = 7000.);
	bool passesMETCut(float = -1., float = 7000.);
	bool passesZVeto(int, int, gChannel, float = 15.); // cut with mZ +/- cut value
	bool passesZVeto(float = 15.); // cut with mZ +/- cut value
	bool passesMllEventVeto(float = 5.);
	bool passesMllEventVeto(int, int, int, float = 5.);

	// HWW Stuff
	bool passesAddLepVeto(int, int, int);

	bool isSigSupMuEvent();
	bool isZMuMuEvent();

	bool isSigSupElEvent();
	bool isZElElEvent(int&);

	bool isGenMatchedSUSYDiLepEvent();
	bool isGenMatchedSUSYDiLepEvent(int&, int&);

	bool isGenMatchedSUSYEEEvent();
	bool isGenMatchedSUSYEMuEvent();

	bool isSSLLMuEvent(int&, int&);
	bool isSSTTMuEvent(int&, int&);

	bool isSSLLElEvent(int&, int&);
	bool isSSTTElEvent(int&, int&);

	bool isSSLLElMuEvent(int&, int&);
	bool isSSTTElMuEvent(int&, int&);

	bool isGoodMuon(int, float = -1.);
	bool isLooseMuon(int);
	bool isTightMuon(int);
	bool isLooseNoTightMuon(int);
	bool isGoodPrimMuon(int, float = -1.);
	bool isGoodSecMuon(int, float = -1.);

	bool isFakeMuon(int);
	bool isPromptMuon(int);
	bool isChargeMatchedMuon(int);

	bool isGoodElectron(int, float = -1.);
	bool isLooseElectron(int);
	bool isTightElectron(int);
	bool isGoodPrimElectron(int, float = -1.);
	bool isGoodSecElectron(int, float = -1.);

	bool isFakeElectron(int);
	bool isPromptElectron(int);
	bool isChargeMatchedElectron(int);

	bool isBarrelElectron(int);

	bool isGoodJet(int, float = 30.);
	bool isGoodJet_LooseLep(int);

private:
	float fC_minHT;
	float fC_minMet;
	float fC_maxHT;
	float fC_maxMet;
	int   fC_minNjets;
	float fC_minMu1pt;
	float fC_minMu2pt;
	float fC_minEl1pt;
	float fC_minEl2pt;
	float fC_maxMet_Control;
	float fC_maxMt_Control;
	
	void resetHypLeptons();
	void setHypLepton1(int, gChannel);
	void setHypLepton2(int, gChannel);
	lepton fHypLepton1;
	lepton fHypLepton2;
	
	const int     getNPtBins (gChannel);
	const double *getPtBins  (gChannel);
	const int     getNPt2Bins(gChannel);
	const double *getPt2Bins (gChannel);
	const int     getNEtaBins(gChannel);
	const double *getEtaBins (gChannel);
	
	Monitor fCounters[gNSAMPLES][3];
	bool fDoCounting;
	gSample fCurrentSample;
	gChannel fCurrentChannel;
	int fCurrentRun;
	ofstream fOUTSTREAM, fOUTSTREAM2;

	int fChargeSwitch;    // 0 for SS, 1 for OS

	vector<int> fMCBG;    // SM background MC samples
	vector<int> fMCBGSig; // SM background + LM0 signal samples
	vector<int> fMCBGMuEnr;    // SM background MC samples with Muon enriched QCD
	vector<int> fMCBGMuEnrSig; // SM background + LM0 signal samples with Muon enriched QCD
	vector<int> fMuData;  // Muon data samples
	vector<int> fEGData;  // EG data samples
	vector<int> fMuEGData;  // MuEG dataset
	vector<int> fMuHadData;  // Muon data samples
	vector<int> fEleHadData;  // EG data samples
	vector<int> fHighPtData;  // All high pt triggered data
	vector<int> fLowPtData;   // All lepton cross HT triggered data
	
	float fLumiNorm;      // Normalize everything to this luminosity
	float fBinWidthScale; // Normalize bin contents to this width

	vector<Sample*>::iterator fS;
	vector<Sample*> fSamples;
	vector<Sample*> fMCSamples;
	map<TString, Sample*> fSampleMap;	// Mapping of sample to name
	
	vector<float> fSigEv_HI_MM_HT;
	vector<float> fSigEv_HI_MM_MET;
	vector<float> fSigEv_LO_MM_HT;
	vector<float> fSigEv_LO_MM_MET;
	vector<float> fSigEv_HI_EE_HT;
	vector<float> fSigEv_HI_EE_MET;
	vector<float> fSigEv_LO_EE_HT;
	vector<float> fSigEv_LO_EE_MET;
	vector<float> fSigEv_HI_EM_HT;
	vector<float> fSigEv_HI_EM_MET;
	vector<float> fSigEv_LO_EM_HT;
	vector<float> fSigEv_LO_EM_MET;
	
	
	TFile *fStorageFile;
	TString fOutputFileName;
	TLatex *fLatex;
	
	
	TH2D *fH2D_MufRatio;
	TH1D *fH1D_MufRatioPt;
	TH1D *fH1D_MufRatioEta;
	TH2D *fH2D_MupRatio;
	TH1D *fH1D_MupRatioPt;
	TH1D *fH1D_MupRatioEta;

	TH2D *fH2D_ElfRatio;
	TH1D *fH1D_ElfRatioPt;
	TH1D *fH1D_ElfRatioEta;
	TH2D *fH2D_ElpRatio;
	TH1D *fH1D_ElpRatioPt;
	TH1D *fH1D_ElpRatioEta;
};

#endif
