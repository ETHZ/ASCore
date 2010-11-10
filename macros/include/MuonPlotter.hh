/*****************************************************************************
*   Collection of tools for producing plots for Muon Fake Rate Analysis      *
*                                                                            *
*                                                  (c) 2010 Benjamin Stieger *
*****************************************************************************/
#ifndef MuonPlotter_HH
#define MuonPlotter_HH

#include "helper/AnaClass.hh"
#include "helper/FPRatios.hh"

class MuonPlotter : public AnaClass{

public:
	MuonPlotter();
	MuonPlotter(TString);
	MuonPlotter(TString, TString);
	virtual ~MuonPlotter();

	void init(TString filename = "samples.dat");
	void loadSamples(const char* filename = "samples.dat");

	void makePlots();

	//////////////////////////////
	// Plots
	void makefRatioPlots();
	void makepRatioPlots();
	void makeIsolationPlots();
	void makePtPlots();
	
	void makeDiffPredictionPlots();
	void makeIntPrediction();
	
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

	void calculateRatio(vector<int>, int, TH2D*&, bool = false);
	void calculateRatio(vector<int>, int, TH2D*&, TH1D*&, TH1D*&, bool = false);
	void calculateRatio(vector<int>, int, float&, float&, bool = false);

	//////////////////////////////
	// Predictions
	void makeSSPredictionPlots(vector<int>);
	void NObs(TH1D *&, vector<int>, bool(MuonPlotter::*)());
	void NObs(TH1D *&, vector<int>);
	vector<TH1D*> NsigPredFromFPRatios(const int, bool = false);
	
	void fillYields();                 // All samples
	void fillYields(int);              // One sample
	void fillYields(vector<int>); // List of samples
	void printYields();
	void printYields(float);
	void printYields(int);
	void printYields(vector<int>, float = -1.0);

	//////////////////////////////
	// I/O
	void bookHistos();
	void writeHistos();
	int readHistos(TString);


	// Event and Object selectors:
	bool isGoodEvent();
	bool passesNJetCut(int=2);
	bool passesHTCut(float);
	bool passesMETCut(float = -1.);
	bool passesZVeto(float = 15.); // cut with mZ +/- cut value
	bool passesZVeto(float, float); // cut in this window
	bool passesMllEventVeto(float = 12.);

	bool isMuTriggeredEvent();
	bool isJetTriggeredEvent();
	bool isHTTriggeredEvent();

	bool isSignalSuppressedEvent();
	bool isSignalSuppressedEventTRG();
	bool isZMuMuEvent();
	bool isZMuMuEventTRG();
	bool isGenMatchedSUSYDiLepEvent();
	bool isSSLLEvent();
	bool isSSLLEventTRG();
	bool isSSTTEvent();
	bool isSSTTEventTRG();

	bool isGoodMuon(int);
	bool isLooseMuon(int);
	bool isTightMuon(int);
	bool isLooseNoTightMuon(int);
	bool isGoodPrimMuon(int);
	bool isGoodSecMuon(int);
	bool isFakeTTbarMuon(int);
	bool isPromptTTbarMuon(int);
	bool isPromptSUSYMuon(int);

	bool isTightElectron(int);
	bool isGoodJet(int);

private:

	FPRatios *fFPRatios;

	vector<int> fAllSamples;
	vector<int> fMCBG;    // SM background MC samples
	vector<int> fMCBGSig; // SM background + LM0 signal samples
	vector<int> fMuData;  // Muon data samples
	vector<int> fEGData;  // EG data samples
	vector<int> fJMData;  // JetMET data samples

	struct numberset{
		long nt2;
		long nt1;
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

		TH2D *h_ratio;
		TH1D *h_ratio_pt;
		TH1D *h_ratio_eta;
	};

	struct NThistos{
		TH2D *h_nt2;
		TH1D *h_nt2_pt;
		TH1D *h_nt2_eta;
		TH2D *h_nt1;
		TH1D *h_nt1_pt;
		TH1D *h_nt1_eta;
		TH2D *h_nt0;
		TH1D *h_nt0_pt;
		TH1D *h_nt0_eta;
	};
	
	struct sample{
		TString name;
		TString sname;
		TFile *file;
		TTree *tree;
		float lumi;
		bool isdata;
		NThistos nthistos;
		numberset numbers;
		lthistos fhistos;
		lthistos phistos;
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
