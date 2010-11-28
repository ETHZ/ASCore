/*****************************************************************************
*   Small Class to make plots for MassAnalysis                               *
*****************************************************************************/

#ifndef MassPlotter_HH
#define MassPlotter_HH

#include "MT2tree.hh"
#include "helper/Utilities.hh"
#include "THStack.h"
#include "TTree.h"

static const int gNJPtbins =99;
static const double gJPtbins[gNJPtbins+1] = {0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 210, 220, 230, 240, 250, 260, 270, 280, 290, 300, 310, 320, 330, 340, 350, 360, 370, 380, 390, 400, 410, 420, 430, 440, 450, 460, 470, 480, 490, 500, 510, 520, 530, 540, 550, 560, 570, 580, 590, 600, 610, 620, 630, 640, 650, 660, 670, 680, 690, 700, 710, 720, 730, 740, 750, 760, 770, 780, 790, 800, 810, 820, 830, 840, 850, 860, 870, 880, 890, 900, 910, 920, 930, 940, 950, 960, 970, 980, 990};

static const int gNMT2bins                   = 19;
static const double  gMT2bins[gNMT2bins+1]   = {0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 60, 70, 85, 105, 135, 180, 260, 360, 500}; 	

static const int gNMT2Leptonbins                   = 11;
static const double  gMT2Leptonbins[gNMT2Leptonbins+1]   = {0, 10, 20, 30, 45,  65,  100, 140, 180, 260, 360, 500}; 	

static const int gNMCTLeptonbins                   = 14;
static const double  gMCTLeptonbins[gNMCTLeptonbins+1]   = {0, 50, 100, 150,  200,  250,  350, 450, 600, 750, 1000, 1250, 1500, 1750, 2000 }; 	

static const int gNMT2predbins                   = 9;
static const double  gMT2predbins[gNMT2predbins+1]   = {50, 60, 70, 85, 105, 135, 180, 260, 360, 500}; 	

static const int gNMT2Normbins                       = 16;
static const double  gMT2Normbins[gNMT2Normbins+1]   = {0, 0.033, 0.067, 0.1, 0.133, 0.167, 0.2, 0.25, 0.3, 0.35, 0.4, 0.5, 0.6, 0.8, 1.00, 1.5, 2.0}; 	

static const int gNDPhibins                  = 20;
static const double  gDPhibins[gNDPhibins+1] = {0, 0.16, 0.32, 0.48, 0.64, 0.8, 0.96, 1.12, 1.28, 1.44, 1.6, 1.76, 1.92, 2.08, 2.24, 
	                                        2.4, 2.56, 2.72, 2.88, 3.04, 3.2};
static const int gNHTbins                    = 8;
static const double  gHTbins[gNHTbins+1]     = {200., 250., 300., 350., 400., 450., 500., 550., 600.}; 


static const double VectorSumPtCut=50;
static const int    LeptConfigCut =9;
static const int    NJetsCut      =3;
static const int    MaxNJetsCut   =100;
static const double dPhiJetsMET   =0.0;
static const double DPhiJetsMetMinJpt=50;

static const int    PrintRunLumiEvt=1;

//________________________________________________________________________________
class MassPlotter  {

public:
	MassPlotter();
	MassPlotter(TString);
	MassPlotter(TString, TString);
	virtual ~MassPlotter();

	void init(TString filename = "samples.dat");
	void loadSamples(const char* filename = "samples.dat");

	void makePlots();

	void setVerbose(int v){ fVerbose = v;};
	void setOutputDir(TString dir){ fOutputDir = Util::MakeOutputDir(dir); };
	void setOutputFile(TString filename){ fOutputFile = Util::MakeOutputFile(fOutputDir + filename); };

private:

	TString fOutputDir;
	TFile *fOutputFile;
	int fVerbose;
	TString fPath;

	struct sample{
		TString name;
		TString sname;
		TString type;
		TFile *file;
		TTree *tree;
		float xsection;
		float nevents;
		float kfact;
		float lumi;
		int color;
	};

	std::vector<sample>  fSamples;
	MT2tree* fMT2tree;
	TTree*   fTree;

	void ControlPlot();
        void MakeMT2PredictionAndPlots(bool cleaned , double dPhisplit[], double fudgefactor);
        void MakePlot(std::vector<sample> Samples, TString var="misc.PseudoJetMT2", TString cuts="NJets>=2", 
		      TString xtitle="MT2 [GeV]", const int nbins=50, const double min=0, const double max=1, 
		      bool cleaned=false, bool logflag=true, bool composited=false, bool ratio=false, 
		      bool stacked=true, bool overlaySUSY=false, float overlayScale = 0);
        void MakePlot(std::vector<sample> Samples, TString var="misc.PseudoJetMT2", TString cuts="NJets>=2", 
		      TString xtitle="MT2 [GeV]", const int nbins=gNMT2bins, const double *bins=gMT2bins, 
		      bool cleaned=false, bool logflag=true, bool composited=false, bool ratio=false, 
		      bool stacked=true, bool overlaySUSY=false, float overlayScale = 0);
	void MakePlot(std::vector<sample> Samples, TString branch_name, const int nbins, const double bins[], bool cleaned, bool logflag, 
		      TString cut_branch_name, double ysplit[], TString option, TString version, TString prediction, double factor );
	void printHisto(THStack h, TString canvname, Option_t *drawopt,  bool logflag);
	void printHisto(THStack h, TH1* h_data, TLegend* leg,  TString canvname, Option_t *drawopt,  bool logflag, TString xtitle, TString ytitle);
        void printHisto(THStack h, TH1* h_data, TH1* h_prediction, TLegend* leg,  TString canvname, Option_t *drawopt, bool logflag, TString xtitle, TString ytitle, bool stacked=true);
	void printHisto(THStack h, TH1* h_data, TH1* h_prediction, TH1* h_susy, TLegend* leg,  TString canvname, Option_t *drawopt, bool logflag, TString xtitle, TString ytitle, float overlayScale=0);
	void printHisto(TH1* h, TString canvname, Option_t *drawopt, bool logflag);
	void ABCD_MT2(TString branch_name, double ysplit[], TString option, const int nbins, const double bins[], bool cleaned, TString sname, TString type);
	void ABCD_MT2(TString branch_name, double ysplit[], TString option, const int nbins, const double bins[], TString version, bool cleaned, TString sname, TString type);
	void plotRatio(TH1* h1, TH1* h2, bool logflag, bool normalize, TString name, TLegend* leg, TString xtitle, TString ytitle);
	void plotRatioStack(THStack hstack, TH1* h1, TH1* h2, bool logflag, bool normalize, TString name, TLegend* leg, TString xtitle, TString ytitle);
	void plotRatioStack(THStack hstack, TH1* h1, TH1* h2, TH1* h3, bool logflag, bool normalize, TString name, TLegend* leg, TString xtitle, TString ytitle, float overlayScale=0);
	void FillRatioHistdPhi(TString cut_branch_name, double lower_cut, double middle_cut, double upper_cut, TString version, bool cleaned);
	void FillRatioHistHT(TString version, bool cleaned);

	double ExpoFitTesting(double *x, double *par);
	TH1D* FillRatioHist(TString branch_name, float MT2plit[], float ysplit[], TString option, const int nbins, const double bins[], TString version,  bool cleaned);
	TH1D*  GetPrediction(TString branch_name, const int nbins, const double bins[], bool cleaned, TString cut_branch_name, double lower_cut, double upper_cut, bool data, double factor );
};

#endif
