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
static const int gNPtbins = 5;
static const double gPtbins[gNPtbins+1] = {10., 20., 30., 40., 50., 100.};
// static const int gNPtbins = 5;
// static const double gPtbins[gNPtbins+1] = {10., 30., 45., 60., 100., 200.};
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
	void makeIsolationPlots();
	void makePtPlots();
	void makeIsoVsPtPlot(int, int, TCut, int, int, TCut, TString = "IsovsPt", bool = false);
	void makeIsoVsNJetsPlot(int, int, TCut, int, int, TCut, TString = "IsovsNJets", bool = false);
	
	// Fake ratios
	void producefRatioFromTree(int, int, const int, const double*, const int, const double*);
	void producepRatioFromTree(int, int, const int, const double*, const int, const double*);
	void producefRatioFromTreeGenMatch(int, int, const int, const double*, const int, const double*);
	void producepRatioFromTreeGenMatch(int, int, const int, const double*, const int, const double*);
	double getFakeRatio(double, double);
	double getFakeRatioPt(double);
	double getFakeRatioEta(double);

	void makeWJetsNt2PredFromTree();

	void makeTTbarNsigPredictionPlots();
	void makeSSNsigPredictionPlots();
	void makeSSPredictionPlots(std::vector<int>);
	
	void NsigObsTTbar(TH1D*&);
	void NsigObsSSLM0(TH1D*&, std::vector<int>);
	void NT2ObsSS(TH1D*&, std::vector<int>);
	
	std::vector<TH1D*> NsigPredFromFPRatios(const int, bool = false);

	// Old
	void loadFakeRatio(int sample, int muon = 1);
	void makeWJetsNt2PredFromHistos();
	TH1D *NsigPredFromTree(const int, const int, const double*, bool = false);
	void makeWJetsFRPredictionPlotsOld();
	void makeSSFRPredictionPlots();
	void plotOrigin(int, int, int);
	void fillMuControlPlots(TCut, TCut, int = 7);

	// Utilities
	// void plotVar(const char* var, const int sample, int nbins, double xmin, double xmax, const TCut reqs, TString ofilename="ofilename", bool logy = false, Option_t *drawopt = "", double line1x = -999.99, double line2x = -999.99, TFile* file = 0 );
private:

	FPRatios *fFPRatios;
	
	struct ratio_histos{
		TH2D *munt2;
		TH2D *munt1;
		TH2D *munt0;
		TH2D *mu1tight;
		TH2D *mu1loose;
		TH2D *mu1loosenotight;
		TH2D *mu1ratio1;
		TH1D *mu1ratio1pt;
		TH1D *mu1ratio1eta;
		TH2D *mu1ratio2;
		TH1D *mu1ratio2pt;
		TH1D *mu1ratio2eta;
		TH2D *mu2tight;
		TH2D *mu2loose;
		TH2D *mu2loosenotight;
		TH2D *mu2ratio1;
		TH1D *mu2ratio1pt;
		TH1D *mu2ratio1eta;
		TH2D *mu2ratio2;
		TH1D *mu2ratio2pt;
		TH1D *mu2ratio2eta;
	};
	ratio_histos getRatios(TFile *file);

	struct sample{
		TString name;
		TString sname;
		TFile *file;
		TTree *tree;
		ratio_histos histos;
		float lumi;
		int color;
	};
	
	int fNJetsMin; // Cut on minimal number of jets
	float fLumiNorm; // Normalize everything to this luminosity
	float fBinWidthScale; // Normalize bin contents to this width

	std::vector<sample> fSamples;
	TH2D *fH2D_fRatio;
	TH1D *fH1D_fRatioPt;
	TH1D *fH1D_fRatioEta;
	TH2D *fH2D_pRatio;
	TH1D *fH1D_pRatioPt;
	TH1D *fH1D_pRatioEta;
	double fD_fRatio, fD_fRatioE;
	double fD_pRatio, fD_pRatioE;
};

#endif
