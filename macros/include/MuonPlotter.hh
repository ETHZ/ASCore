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


class MuonPlotter : public AnaClass{

public:
	
	// This enum has to correspond to the content of the samples.dat file
	enum gSample {
		sample_begin,
		MuA = sample_begin, MuB, EGA, EGB, JMA, JMB, MultiJet,
		TTbar, WJets, ZJets, VVJets, QCD15, QCD30, QCD80, QCD170,
		SSWWDPS, SSWWSPSPos, SSWWSPSNeg,
		LM0, InclMu,
		gNSAMPLES
	};
	enum gChannel {
		Muon,
		EMu,
		Electron
	};
	struct lepton{
		TLorentzVector p;
		int charge;
		int type; // -1(unknown), 0(mu), 1(ele)
		int index;
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
	
	//////////////////////////////
	// Plots
	void makefRatioPlots();
	void makepRatioPlots();
	void makeIsolationPlots();
	void makeIsolationPlot();
	void makePtPlots();
	
	void makeDiffPredictionPlots();
	void makeIntPrediction(TString);
	void makeIntPredictionMuMu(vector<int>);
	
	void makeIsoVsPtPlot(int, int, TCut, int, int, TCut, TString = "IsovsPt", bool = false);
	void makeIsoVsPtPlot(int, int, bool(MuonPlotter::*)(), bool(MuonPlotter::*)(int), int, int, bool(MuonPlotter::*)(), bool(MuonPlotter::*)(int), TString = "IsovsPt", bool = false);
	void makeIsoVsPtPlot(vector<int>, int, bool(MuonPlotter::*)(), bool(MuonPlotter::*)(int), vector<int>, int, bool(MuonPlotter::*)(), bool(MuonPlotter::*)(int), TString = "IsovsPt", bool = false);
	void makeIsoVsNJetsPlot(int, int, TCut, int, int, TCut, TString = "IsovsNJets", bool = false);
	
	//////////////////////////////
	// Fake ratios
	// Produce from tree, with given selections:
	void produceRatio(int, int, bool(MuonPlotter::*)(), bool(MuonPlotter::*)(int), TH2D*&, TH1D*&, TH1D*&, bool = false);
	void produceRatio(vector<int>, int, bool(MuonPlotter::*)(), bool(MuonPlotter::*)(int), TH2D*&, TH1D*&, TH1D*&, bool = false);
	vector<double> produceRatio(vector<int>, int, bool(MuonPlotter::*)(), bool(MuonPlotter::*)(int));

	void fillfRatio(int, int);
	void fillfRatio(vector<int>, int);
	void fillfRatio(int, int, const int, const double*, const int, const double*);
	void fillfRatio(vector<int>, int, const int, const double*, const int, const double*);

	void fillpRatio(int, int);
	void fillpRatio(vector<int>, int);
	void fillpRatio(int, int, const int, const double*, const int, const double*);
	void fillpRatio(vector<int>, int, const int, const double*, const int, const double*);
	
	TH1D* fillRatioPt(int, int, bool(MuonPlotter::*)(), bool(MuonPlotter::*)(int), bool = false);
	TH1D* fillRatioPt(vector<int>, int, bool(MuonPlotter::*)(), bool(MuonPlotter::*)(int), bool = false);
	TH1D* fillRatioPt(vector<int>, int, bool(MuonPlotter::*)(), bool(MuonPlotter::*)(int), const int, const double*, const int, const double*, bool = false);
	TH2D* fillRatio(int, int, bool(MuonPlotter::*)(), bool(MuonPlotter::*)(int), const int, const double*, const int, const double*);
	TH2D* fillRatio(vector<int>, int, bool(MuonPlotter::*)(), bool(MuonPlotter::*)(int), const int, const double*, const int, const double*);

	void plotRatio(int, int, bool(MuonPlotter::*)(), bool(MuonPlotter::*)(int), TString = "");
	void plotRatio(vector<int>, int, bool(MuonPlotter::*)(), bool(MuonPlotter::*)(int), TString = "");

	// Calculate from pre stored numbers, with fixed selections:
	TH1D* fillRatioPt(int, int, bool = false);
	TH1D* fillRatioPt(vector<int>, int, bool = false);

	void calculateRatio(vector<int>, gChannel, int, TH2D*&, bool = false);
	void calculateRatio(vector<int>, gChannel, int, TH2D*&, TH1D*&, TH1D*&, bool = false);
	void calculateRatio(vector<int>, gChannel, int, float&, float&, bool = false);

	//////////////////////////////
	// Predictions
	void makeSSPredictionPlots(vector<int>);
	void NObs(TH1D *&, vector<int>, bool(MuonPlotter::*)());
	void NObs(TH1D *&, vector<int>);
	vector<TH1D*> NsigPredFromFPRatios(const int, bool = false);
	
	void fillYields();                 // All samples
	void fillYields(int);              // One sample
	void fillYields(vector<int>); // List of samples
	
	void initCounters(int = -1);
	void storeNumbers();
	void storeNumbers(gChannel);
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
	bool isZMuMuEvent();
	bool isZMuMuEventTRG();

	bool isSigSupElEvent();
	bool isSigSupElEventTRG();
	bool isZElElEvent(int&);
	bool isZElElEventTRG(int&);

	bool isGenMatchedSUSYDiLepEvent();
	bool isGenMatchedSUSYDiLepEvent(int&, int&);

	int isSSLLEvent(int&, int&);
	int isOSLLEvent(int&, int&);
	vector<lepton> sortLeptonsByPt(vector<lepton> &leptons);

	bool isSSLLMuEvent(   int&, int&);
	bool isSSLLMuEventTRG(int&, int&);
	bool isSSTTMuEvent(   int&, int&);
	bool isSSTTMuEventTRG(int&, int&);

	bool isSSLLElEvent(   int&, int&);
	bool isSSLLElEventTRG(int&, int&);
	bool isSSTTElEvent(   int&, int&);
	bool isSSTTElEventTRG(int&, int&);

	bool isSSLLElMuEvent(   int&, int&);
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

	bool isGoodJet(int);
	bool isGoodJet_LooseLep(int);

private:
	const int     getNPtBins();
	const double *getPtBins();
	const int     getNPt2Bins();
	const double *getPt2Bins();
	const int     getNEtaBins();
	const double *getEtaBins();
	
	Monitor fCounters[gNSAMPLES][3];
	bool fDoCounting;
	gSample fCurrentSample;
	gChannel fCurrentChannel;
	ofstream fOUTSTREAM;

	int fSelectionSwitch; // 0 for UCSD, 1 for UFlorida
	int fChargeSwitch;    // 0 for SS, 1 for OS

	vector<int> fAllSamples;
	vector<int> fMCBG;    // SM background MC samples
	vector<int> fMCBGMuEnr;    // SM background MC samples
	vector<int> fMCBGSig; // SM background + LM0 signal samples
	vector<int> fMuData;  // Muon data samples
	vector<int> fEGData;  // EG data samples
	vector<int> fJMData;  // JetMET data samples
	
	struct numberset{
		long nt2;
		long nt10;
		long nt01;
		long nt0;
		long nsst;
		long nssl;
		long nzl;
		long nzt;
	};
	
	struct lthistos{
		TH2D *h_ntight;
		TH2D *h_nloose;
		TH1D *h_ntight_pt;
		TH1D *h_nloose_pt;
		TH1D *h_ntight_eta;
		TH1D *h_nloose_eta;
	};

	struct NThistos{
		TH2D *h_nt2;
		TH1D *h_nt2_pt;
		TH1D *h_nt2_eta;
		TH2D *h_nt10;    // Only this is used for EE/MuMu
		TH1D *h_nt10_pt;
		TH1D *h_nt10_eta;
		TH2D *h_nt01;    // See notation in FR Note
		TH1D *h_nt01_pt;
		TH1D *h_nt01_eta;
		TH2D *h_nt0;
		TH1D *h_nt0_pt;
		TH1D *h_nt0_eta;
	};
	
	struct channel{
		NThistos nthistos;
		lthistos fhistos;
		lthistos phistos;
		numberset numbers;
	};
	
	struct sample{
		TString name;
		TString sname;
		TFile *file;
		TTree *tree;
		float lumi;
		int color;
		bool isdata;
		channel mumu;
		channel emu;
		channel ee;
	};
	
	int fNJetsMin; // Cut on minimal number of jets
	float fMinPt1; // Cut on pt of harder muon
	float fMinPt2; // Cut on pt of softer muon

	float fLumiNorm;      // Normalize everything to this luminosity
	float fBinWidthScale; // Normalize bin contents to this width

	vector<sample> fSamples;
	map<TString, int> fSampleMap;	// Mapping of sample number to name
	
	TFile *fStorageFile;
	TString fOutputFileName;
	
	TH2D *fH2D_fRatio;
	TH1D *fH1D_fRatioPt;
	TH1D *fH1D_fRatioEta;
	TH2D *fH2D_pRatio;
	TH1D *fH1D_pRatioPt;
	TH1D *fH1D_pRatioEta;
};

#endif
