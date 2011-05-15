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
	static const int gNMuEtabins = 1;
	static const int gNElPtbins  = 5;
	static const int gNElPt2bins = 5;
	static const int gNElEtabins = 1;
	
	static const int gNNVrtxBins = 10;
	static int gNVrtxBins[gNNVrtxBins+1];
	
	static double gMuPtbins [gNMuPtbins+1];
	static double gMuPt2bins[gNMuPt2bins+1];
	static double gMuEtabins[gNMuEtabins+1];

	static double gElPtbins [gNElPtbins+1];
	static double gElPt2bins[gNElPt2bins+1];
	static double gElEtabins[gNElEtabins+1];


	// This enum has to correspond to the content of the samples.dat file
	enum gSample {
		sample_begin,
		// DoubleMu = sample_begin, DoubleEle, MuEG, MuHad, ElectronHad,
		DoubleMu1 = sample_begin, DoubleMu2, DoubleEle1, DoubleEle2, MuEG1, MuEG2,
		TTJets, WJets, DYJets, WW, WZ, ZZ,
		LM0, LM1, LM2, LM3, LM4, LM5, LM6, LM7, LM8, LM9, LM11, LM12, LM13, 
		// LM0, LM9, LM11,
		QCDMuEnr10,
		QCD5to15,
		QCD15to30,
		QCD30to50,
		QCD50to80,
		QCD80to120,
		QCD120to170,
		QCD170to300,
		QCD300to470,
		QCD470to600,
		QCD600to800,
		QCD800to1000,
		QCD1000to1400,
		QCD1400to1800,
		QCD1800,
		// QCD50to100MG,
		// QCD100to250MG,
		// QCD250to500MG,
		// QCD500to1000MG,
		// QCD1000MG,
		gNSAMPLES
	};
	enum gOldSample { // temporary
		oldsample_begin,
		MuA = sample_begin, MuB, EGA, EGB, JMA, JMB, MultiJet,
		TTbar, ZJets, AstarJets, VVJets, QCD15, QCD30, QCD80, QCD170,
		SSWWDPS, SSWWSPSPos, SSWWSPSNeg,
		InclMu,
		gNOLDSAMPLES
	};
	enum gRegion {
		region_begin,
		Signal = region_begin,
		SigSup,
		ZDecay,
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

		TH2D *nt11_origin;
		TH2D *nt10_origin;
		TH2D *nt01_origin;
		TH2D *nt00_origin;
		TH1D *sst_origin;
		TH1D *ssl_origin;
		TH1D *zt_origin;
		TH1D *zl_origin;
	};
	
	struct Region{ // different binnings and or selections cuts, e.g. florida vs surfturf
		TString name;
		TString sname;
		Channel mm;
		Channel em;
		Channel ee;
	};
	
	static const int gNKinVars = 10;
	struct KinPlots{
		static TString var_name[gNKinVars];
		static int nbins[gNKinVars];
		static float xmin[gNKinVars];
		static float xmax[gNKinVars];
		TH1D *hvar[gNKinVars];
	};
	
	static const int gNSels = 2;
	struct IsoPlots{
		static TString sel_name[gNSels];
		static int nbins[gNSels];
		TH1D *hiso[gNSels];
		TH1D *hiso_pt[gNSels][gNMuPt2bins];
		TH1D *hiso_nv[gNSels][gNNVrtxBins];
	};
	
	struct Sample{
		Sample(){};
		TString name;
		TString sname;
		TFile *file;
		TTree *tree;
		float lumi;
		int color;
		int datamc; // 0: Data, 1: SM MC, 2: Signal MC
		Region region[gNREGIONS];
		NumberSet numbers[gNCHANNELS]; // summary of integrated numbers
		KinPlots kinplots[3]; // tt and ll and signal
		IsoPlots isoplots[2]; // e and mu
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
	void makeMufRatioPlots(bool = false);
	void makeMupRatioPlots(bool = false);
	void makeElfRatioPlots(bool = false);
	void makeElpRatioPlots(bool = false);

	void makeMufEffPlots(bool = false);
	void makeElfEffPlots(bool = false);
	void makeMupEffPlots(bool = false);
	void makeElpEffPlots(bool = false);

	void makeMuIsolationPlotsOld();
	
	void makeMuIsolationPlots();
	void makeElIsolationPlots();
	void fillMuIsoHistos(gSample);
	void fillElIsoHistos(gSample);
	
	void makeNT2KinPlots();
	void fillKinematicHistos(gSample);
	
	void makeMCClosurePlots(vector<int>);
	void makeDataClosurePlots();
	void makeNT012Plots(vector<int>, gChannel, gRegion = Signal);
	void makeNT012Plots(gChannel, vector<int>, bool(MuonPlotter::*)(int&, int&), TString = "");

	void makeIntPrediction(TString);
	
	void makeMuIsoVsPtPlot(int, int, TCut, int, int, TCut, TString = "IsovsPt", bool = false);
	void makeMuIsoVsPtPlot(int, int, bool(MuonPlotter::*)(), bool(MuonPlotter::*)(int), int, int, bool(MuonPlotter::*)(), bool(MuonPlotter::*)(int), TString = "IsovsPt", bool = false);
	void makeMuIsoVsPtPlot(vector<int>, int, bool(MuonPlotter::*)(), bool(MuonPlotter::*)(int), vector<int>, int, bool(MuonPlotter::*)(), bool(MuonPlotter::*)(int), TString = "IsovsPt", bool = false);
	void makeMuIsoVsNJetsPlot(int, int, TCut, int, int, TCut, TString = "IsovsNJets", bool = false);
	
	//////////////////////////////
	// Fake ratios
	// Produce from tree, with given selections:
	void produceRatio(gChannel, int, int, bool(MuonPlotter::*)(), bool(MuonPlotter::*)(int), TH2D*&, TH1D*&, TH1D*&, bool = false);
	void produceRatio(gChannel, vector<int>, int, bool(MuonPlotter::*)(), bool(MuonPlotter::*)(int), TH2D*&, TH1D*&, TH1D*&, bool = false);

	TH1D* fillMuRatioPt(int, int, bool(MuonPlotter::*)(), bool(MuonPlotter::*)(int), bool = false);
	TH1D* fillMuRatioPt(vector<int>, int, bool(MuonPlotter::*)(), bool(MuonPlotter::*)(int), bool = false);
	TH1D* fillMuRatioPt(vector<int>, int, bool(MuonPlotter::*)(), bool(MuonPlotter::*)(int), const int, const double*, const int, const double*, bool = false);
	TH2D* fillMuRatio(int, int, bool(MuonPlotter::*)(), bool(MuonPlotter::*)(int), const int, const double*, const int, const double*);
	TH2D* fillMuRatio(vector<int>, int, bool(MuonPlotter::*)(), bool(MuonPlotter::*)(int), const int, const double*, const int, const double*);

	// Calculate from pre stored numbers, with fixed selections:
	void fillMuElRatios(vector<int>);

	TH1D* fillMuRatioPt(int, gRegion, bool = false);
	TH1D* fillMuRatioPt(vector<int>, gRegion, bool = false);
	TH1D* fillElRatioPt(int, gRegion, bool = false);
	TH1D* fillElRatioPt(vector<int>, gRegion, bool = false);

	void calculateRatio(vector<int>, gChannel, gRegion, TH2D*&, bool = false);
	void calculateRatio(vector<int>, gChannel, gRegion, TH2D*&, TH1D*&, TH1D*&, bool = false);
	void calculateRatio(vector<int>, gChannel, gRegion, float&, float&);
	void calculateRatio(vector<int>, gChannel, gRegion, float&, float&, float&);
	
	TEfficiency *mergeDataEfficiencies(vector<int>, gChannel, gRegion, bool = false, TEfficiency::EStatOption = TEfficiency::kBUniform, double beta = 1., double alpha = 1.);
	TEfficiency *getEfficiency(Sample*, gChannel, gRegion, int = 0, bool = false);
	TGraphAsymmErrors *combineMCEfficiencies(vector<int>, gChannel, gRegion, bool = false, TEfficiency::EStatOption = TEfficiency::kBUniform, double beta = 1., double alpha = 1.);

	void getPassedTotal(vector<int>, gChannel, gRegion, TH2D*&, TH2D*&, bool = false);

	void ratioWithBinomErrors(float, float, float&, float&);
	void ratioWithPoissErrors(float, float, float&, float&);
	void ratioWithAsymmCPErrors(int, int, float&, float&, float&);

	//////////////////////////////
	// Predictions
	void makeSSMuMuPredictionPlots(vector<int>, bool = false);
	void makeSSElElPredictionPlots(vector<int>, bool = false);
	void makeSSElMuPredictionPlots(vector<int>, bool = false);
	void NObs(gChannel, TH1D *&, vector<int>, bool(MuonPlotter::*)());
	void NObs(gChannel, TH1D *&, vector<int>, gRegion = Signal);
	void NObs(gChannel, THStack *&, vector<int>, gRegion = Signal);
	vector<TH1D*> MuMuFPPrediction(TH2D* fratio, TH2D* pratio, TH2D* nt2, TH2D* nt1, TH2D* nt0, bool output = false);
	vector<TH1D*> ElElFPPrediction(TH2D* fratio, TH2D* pratio, TH2D* nt2, TH2D* nt1, TH2D* nt0, bool output = false);
	vector<TH1D*> ElMuFPPrediction(TH2D* mufratio, TH2D* mupratio, TH2D* elfratio, TH2D* elpratio,  TH2D* nt2, TH2D* nt10, TH2D* nt01, TH2D* nt0, bool output = false);
	
	void fillYields(Sample*);
	
	void initCounters(gSample);
	void storeNumbers(Sample*, gChannel);
	void printCutFlows(TString);
	
	void printYields(gChannel, float = -1.0);
	void printYieldsShort(float = -1);

	//////////////////////////////
	// I/O
	void bookHistos();
	void writeHistos();
	int readHistos(TString);
	
	void bookRatioHistos();
	void fixPRatios();
	
	// Geninfo stuff:
	int muIndexToBin(int);
	int elIndexToBin(int);
	TString muBinToLabel(int);
	TString elBinToLabel(int);
	void labelMuOriginAxis(TAxis *);
	void labelElOriginAxis(TAxis *);
	
	void printOrigins();
	void printMuOriginTable();
	void printMuOriginHeader(TString);
	void printMuOriginFromSample(Sample*, int);
	void print2MuOriginsFromSample(Sample*, int);

	void printElOriginTable();
	void printElOriginHeader(TString);
	void printElOriginFromSample(Sample*, int);
	void print2ElOriginsFromSample(Sample*, int);

	void printEMuOriginTable();
	void printEMuOriginHeader(TString);
	void printEMuOriginsFromSample(Sample*, int);

	void printOriginSummary(int, gChannel);
	void printOriginSummary2L(int, gChannel);

	// Event and Object selectors:
	int getNJets();
	float getHT();
	float getMT2(int, int, int);
	
	bool isGoodEvent();
	bool isGoodMuEvent();
	int hasLooseMuons(int&, int&);
	int hasLooseMuons();
	int hasLooseElectrons(int&, int&);
	int hasLooseElectrons();
	bool passesNJetCut(int=2);
	bool passesNJetCut_LooseLep(int=2);
	bool passesJet50Cut();
	
	bool passesHTCut(float);
	bool passesMETCut(float = -1.);
	bool passesZVeto(float = 15.); // cut with mZ +/- cut value
	bool passesMllEventVeto(float = 12.);

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

	bool isSigSupMuEvent();
	bool isZMuMuEvent();

	bool isSigSupElEvent();
	bool isZElElEvent(int&);

	bool isGenMatchedSUSYDiLepEvent();
	bool isGenMatchedSUSYDiLepEvent(int&, int&);

	bool isGenMatchedSUSYEEEvent();
	bool isGenMatchedSUSYEMuEvent();

	int isSSLLEvent(int&, int&);
	int isOSLLEvent(int&, int&);
	vector<lepton> sortLeptonsByPt(vector<lepton> &leptons);

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
	bool isFakeTTbarMuon(int);
	bool isPromptTTbarMuon(int);
	bool isPromptSUSYMuon(int);

	bool isGoodElectron(int, float = -1.);
	bool isLooseElectron(int);
	bool isTightElectron(int);
	bool isGoodPrimElectron(int, float = -1.);
	bool isGoodSecElectron(int, float = -1.);
	bool isPromptSUSYElectron(int);

	bool isGoodJet(int, float = 30.);
	bool isGoodJet_LooseLep(int);

private:
	float fC_minHT;
	float fC_minMet;
	int   fC_minNjets;
	float fC_minMu1pt;
	float fC_minMu2pt;
	float fC_minEl1pt;
	float fC_minEl2pt;
	
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
	ofstream fOUTSTREAM;
	map<string, int> fHLTLabelMap; // Mapping of HLT trigger bit names

	int fChargeSwitch;    // 0 for SS, 1 for OS

	vector<int> fMCBG;    // SM background MC samples
	vector<int> fMCBGSig; // SM background + LM0 signal samples
	vector<int> fMCBGMuEnr;    // SM background MC samples with Muon enriched QCD
	vector<int> fMCBGMuEnrSig; // SM background + LM0 signal samples with Muon enriched QCD
	vector<int> fMuData;  // Muon data samples
	vector<int> fEGData;  // EG data samples
	vector<int> fMuEGData;  // MuEG dataset
	vector<int> fJMData;  // JetMET data samples
	
	float fLumiNorm;      // Normalize everything to this luminosity
	float fBinWidthScale; // Normalize bin contents to this width

	vector<Sample*> fSamples;
	map<TString, int> fSampleMap;	// Mapping of sample number to name
	
	TFile *fStorageFile;
	TString fOutputFileName;
	
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
