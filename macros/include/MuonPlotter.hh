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
#include "TLorentzVector.h"

#include "TEfficiency.h"
#include "TGraphAsymmErrors.h"


class MuonPlotter : public AnaClass{

public:
	
	// This enum has to correspond to the content of the samples.dat file
	enum gSample {
		sample_begin,
		MuA = sample_begin, MuB, EGA, EGB, JMA, JMB, MultiJet,
		TTbar, WJets, ZJets, AstarJets, VVJets, QCD15, QCD30, QCD80, QCD170,
		SSWWDPS, SSWWSPSPos, SSWWSPSNeg,
		LM0, InclMu,
		gNSAMPLES
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
	};
	
	struct Region{ // different binnings and or selections cuts, e.g. florida vs surfturf
		TString name;
		TString sname;
		Channel mm;
		Channel em;
		Channel ee;
	};
	
	struct Sample{
		TString name;
		TString sname;
		TFile *file;
		TTree *tree;
		float lumi;
		int color;
		bool isdata;
		bool isfilled; // check if this has been read in -- not implemented yet
		Region region[gNREGIONS];
		NumberSet numbers[gNCHANNELS]; // summary of integrated numbers
	};
	
	MuonPlotter();
	MuonPlotter(TString);
	MuonPlotter(TString, TString);
	virtual ~MuonPlotter();

	inline void setSelection(int sel){ fSelectionSwitch = sel; };
	inline void setCharge(int charge){ fChargeSwitch = charge; };

	void init(TString filename = "samples.dat");
	void loadSamples(const char* filename = "samples.dat");
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

	void makeMuIsolationPlots();
	void makeMuIsolationPlot();
	void makeMuPtPlots();
	
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

	void plotMuRatio(int, int, bool(MuonPlotter::*)(), bool(MuonPlotter::*)(int), TString = "");
	void plotMuRatio(vector<int>, int, bool(MuonPlotter::*)(), bool(MuonPlotter::*)(int), TString = "");

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
	
	void initCounters(int = -1);
	void storeNumbers(Sample*, gChannel);
	void printCutFlows(TString);
	
	void printYields(gChannel = Muon);
	void printYields(gChannel, float);
	void printYields(gChannel, int, float = -1.0);
	void printYields(gChannel, vector<int>, float = -1.0);
	void printYieldsShort(float = -1);

	//////////////////////////////
	// I/O
	void bookHistos();
	void writeHistos();
	int readHistos(TString);
	
	void bookRatioHistos();
	void fixPRatios();

	// Event and Object selectors:
	bool isGoodEvent();
	bool isGoodMuEvent();
	bool isGoodElEvent();
	bool isGoodElMuEvent();
	bool passesNJetCut(int=2);
	bool passesNJetCut_LooseLep(int=2);
	bool passesHTCut(float);
	bool passesMETCut(float = -1.);
	bool passesZVeto(float = 15.); // cut with mZ +/- cut value
	bool passesMllEventVeto(float = 12.);

	bool isMuTriggeredEvent();
	bool isElTriggeredEvent();
	bool isJetTriggeredEvent();
	bool isHTTriggeredEvent();
	bool isGoodRun(int sample);

	bool isSigSupMuEvent();
	bool isSigSupMuEventTRG();
	bool isSigSupOSMuMuEvent(int&, int&);
	bool isSigSupSSMuMuEvent(int&, int&);
	bool isZMuMuEvent();
	bool isZMuMuEventTRG();

	bool isSigSupElEvent();
	bool isSigSupElEventTRG();
	bool isZElElEvent(int&);
	bool isZElElEventTRG();
	bool isZElElEventTRG(int&);

	bool isGenMatchedSUSYDiLepEvent();
	bool isGenMatchedSUSYDiLepEvent(int&, int&);

	bool isGenMatchedSUSYEEEvent();
	bool isGenMatchedSUSYEMuEvent();

	int isSSLLEvent(int&, int&);
	int isOSLLEvent(int&, int&);
	vector<lepton> sortLeptonsByPt(vector<lepton> &leptons);

	bool isSSLLMuEvent(   int&, int&);
	bool isSSLLMuEventInvMETTRG(int&, int&);
	bool isSSLLMuEventHTControlTRG(int&, int&);
	bool isSSLLMuEventTRG(int&, int&);
	bool isSSTTMuEvent(   int&, int&);
	bool isSSTTMuEventTRG(int&, int&);

	bool isSSLLElEvent(   int&, int&);
	bool isSSLLElEventTRG(int&, int&);
	bool isSSTTElEvent(   int&, int&);
	bool isSSTTElEventTRG(int&, int&);

	bool isSSLLElMuEvent(   int&, int&);
	bool isSSLLElMuEventInvMETTRG(int&, int&);
	bool isOSLLElMuEventTRG(int&, int&);
	bool isSSLLElMuEventTRG(int&, int&);
	bool isSSTTElMuEvent(   int&, int&);
	bool isSSTTElMuEventTRG(int&, int&);

	bool isGoodMuon(int);
	bool isLooseMuon(int);
	bool isTightMuon(int);
	bool isLooseNoTightMuon(int);
	bool isGoodPrimMuon(int);
	bool isGoodSecMuon(int);
	bool isFakeTTbarMuon(int);
	bool isPromptTTbarMuon(int);
	bool isPromptSUSYMuon(int);

	bool isGoodElectron(int);
	bool isLooseElectron(int);
	bool isTightElectron(int);
	bool isGoodPrimElectron(int);
	bool isGoodSecElectron(int);
	bool isPromptSUSYElectron(int);

	bool isGoodJet(int);
	bool isGoodJet_LooseLep(int);

private:
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
	ofstream fOUTSTREAM;

	int fSelectionSwitch; // 0 for UCSD, 1 for UFlorida
	int fChargeSwitch;    // 0 for SS, 1 for OS

	vector<int> fAllSamples;
	vector<int> fMCBG;    // SM background MC samples
	vector<int> fMCBGSig; // SM background + LM0 signal samples
	vector<int> fMCBGMuEnr;    // SM background MC samples with Muon enriched QCD
	vector<int> fMCBGMuEnrSig; // SM background + LM0 signal samples with Muon enriched QCD
	vector<int> fMuData;  // Muon data samples
	vector<int> fEGData;  // EG data samples
	vector<int> fJMData;  // JetMET data samples
	
	int fNJetsMin; // Cut on minimal number of jets
	float fMinPt1; // Cut on pt of harder muon
	float fMinPt2; // Cut on pt of softer muon

	float fLumiNorm;      // Normalize everything to this luminosity
	float fBinWidthScale; // Normalize bin contents to this width

	vector<Sample> fSamples;
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
