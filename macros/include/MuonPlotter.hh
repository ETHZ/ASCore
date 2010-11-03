/*****************************************************************************
*   Collection of tools for producing plots for Muon Fake Rate Analysis      *
*                                                                            *
*                                                  (c) 2010 Benjamin Stieger *
*****************************************************************************/
#ifndef MuonPlotter_HH
#define MuonPlotter_HH

#include "helper/AnaClass.hh"
#include "helper/FPRatios.hh"

// Binning ///////////////////////////////////////////////////////////////////////
// For ttbar Analysis
// static const int gNPtbins = 1;
// static const double gPtbins[gNPtbins+1] = {10., 150.};
// static const int gNPtbins = 5;
// static const double gPtbins[gNPtbins+1] = {10., 15., 20., 35., 65., 150.};
// static const int gNEtabins = 1;
// static const double gEtabins[gNEtabins+1] = {-2.4, 2.4};
// static const int gNEtabins = 5;
// static const double gEtabins[gNEtabins+1] = {-2.4, -1.4, -0.5, 0.5, 1.4, 2.4};

// For SS Analysis
// static const int gNPtbins = 5;
// static const double gPtbins[gNPtbins+1] = {10., 20., 30., 40., 50., 100.};
// static const int gNPtbins = 5;
// static const double gPtbins[gNPtbins+1] = {10., 30., 45., 60., 100., 200.};

// For Analysis on Data
// static const int gNPtbins = 4;
// static const double gPtbins[gNPtbins+1] = {10., 20., 30., 40., 60.};
static const int gNPtbins = 5;
static const double gPtbins[gNPtbins+1] = {20., 30., 40., 50., 75., 100.};
static const int gNPt2bins = 6;
static const double gPt2bins[gNPt2bins+1] = {10., 20., 30., 40., 50., 75., 100.};
static const int gNEtabins = 1;
static const double gEtabins[gNEtabins+1] = {-2.4, 2.4};

//////////////////////////////////////////////////////////////////////////////////

class MuonPlotter : public AnaClass{

public:
	MuonPlotter();
	MuonPlotter(TString);
	MuonPlotter(TString, TString);
	virtual ~MuonPlotter();

	void init(TString filename = "samples.dat");
	void loadSamples(const char* filename = "samples.dat");

	void makePlots();

	void makefRatioPlots();
	void makepRatioPlots();
	void makeIsolationPlots();
	void makePtPlots();
	
	void make100pbMCPredictionPlots();
	void makePrediction();
	
	void makeIsoVsPtPlot(int, int, TCut, int, int, TCut, TString = "IsovsPt", bool = false);
	void makeIsoVsPtPlot(int, int, bool(MuonPlotter::*)(), bool(MuonPlotter::*)(int), int, int, bool(MuonPlotter::*)(), bool(MuonPlotter::*)(int), TString = "IsovsPt", bool = false);
	void makeIsoVsNJetsPlot(int, int, TCut, int, int, TCut, TString = "IsovsNJets", bool = false);
	
	// Fake ratios
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

	void makeSSNsigPredictionPlots();
	void makeSSPredictionPlots(std::vector<int>);
	
	void NObs(TH1D *&, std::vector<int>, bool(MuonPlotter::*)());
	
	std::vector<TH1D*> NsigPredFromFPRatios(const int, bool = false);

	// Event and Object selectors:
	bool isGoodEvent();
	bool isMuTriggeredEvent();
	bool isJetTriggeredEvent();
	bool isSignalSuppressedEvent();
	bool isSignalSuppressedEventTRG();
	bool isZEvent();
	bool isZEventTRG();
	bool isGenMatchedSUSYDiLepEvent();
	bool isSSLLEvent();
	bool isSSLLEventTRG();
	bool isSSTTEvent();

	bool isGoodMuon(int);
	bool isLooseMuon(int);
	bool isTightMuon(int);
	bool isLooseNoTightMuon(int);
	bool isGoodPrimMuon(int);
	bool isGoodSecMuon(int);
	bool isFakeTTbarMuon(int);
	bool isPromptTTbarMuon(int);
	bool isPromptSUSYMuon(int);

private:

	FPRatios *fFPRatios;

	vector<int> fMCBG;    // SM background MC samples
	vector<int> fMCBGSig; // SM background + LM0 signal samples
	vector<int> fMuData;  // Muon data samples
	vector<int> fEGData;  // EG data samples
	vector<int> fJMData;  // JetMET data samples

	struct sample{
		TString name;
		TString sname;
		TFile *file;
		TTree *tree;
		float lumi;
		int color;
	};
	
	int fNJetsMin; // Cut on minimal number of jets
	float fMinPt1; // Cut on pt of harder muon
	float fMinPt2; // Cut on pt of softer muon

	float fLumiNorm;      // Normalize everything to this luminosity
	float fBinWidthScale; // Normalize bin contents to this width

	std::vector<sample> fSamples;
	map<TString, int> fSampleMap;	// Mapping of sample number to name
	
	TH2D *fH2D_fRatio;
	TH1D *fH1D_fRatioPt;
	TH1D *fH1D_fRatioEta;
	TH2D *fH2D_pRatio;
	TH1D *fH1D_pRatioPt;
	TH1D *fH1D_pRatioEta;
};

#endif
