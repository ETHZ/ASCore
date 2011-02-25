/*****************************************************************************
*   Small Class to make plots for MassAnalysis                               *
*****************************************************************************/

#ifndef MassPlotter_HH
#define MassPlotter_HH

#include "MT2tree.hh"
#include "helper/Utilities.hh"
#include "helper/Monitor.hh"
#include "THStack.h"
#include "TTree.h"
#include <map>

static const int gNMT2bins                   = 19;
static const double  gMT2bins[gNMT2bins+1]   = {0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 60, 70, 85, 105, 135, 180, 260, 360, 500}; 	

static const int gNMT2dijetM                    = 20;
static const double  gMT2dijetM[gNMT2dijetM+1]  = {0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 60, 70, 85, 105, 135, 180, 260, 360, 500, 800}; 	

static const int gNMT2Massivebins                          = 16;
static const double  gMT2Massivebins[gNMT2Massivebins+1]   = {0, 10, 20, 30, 40, 50, 60, 70, 85, 105, 135, 180, 260, 360, 500, 700, 900}; 	

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
	void makeZnunu();
	void makeSmallCopy(int nevents, int sample);
	
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
 
	void setVerbose(int v){ fVerbose = v;};
	void setOutputDir(TString dir){ fOutputDir = Util::MakeOutputDir(dir); };
	void setOutputFile(TString filename){ fOutputFile = Util::MakeOutputFile(fOutputDir + filename); };
        void makePlot(TString var="misc.PseudoJetMT2", TString cuts="misc.HBHENoiseFlag == 1", 
		      int njets=-2, int nleps=0, //njets: 2 -> njets==2, -2 -> njets>=2
		      TString xtitle="MT2 [GeV]", const int nbins=50, const double min=0., const double max=1., 
		      bool cleaned=false, bool logflag=true, bool composited=false, bool ratio=false, 
		      bool stacked=true, bool overlaySUSY=false, float overlayScale = 0);
        void MakePlot(TString var="misc.PseudoJetMT2", TString cuts="misc.HBHENoiseFlag == 1", 
		      int njets=-2, int nleps=0, //njets: 2 -> njets==2, -2 -> njets>=2
		      TString xtitle="MT2 [GeV]", const int nbins=gNMT2bins, const double *bins=gMT2bins, 
		      bool cleaned=false, bool logflag=true, bool composited=false, bool ratio=false, 
		      bool stacked=true, bool overlaySUSY=false, float overlayScale = 0);
        void plotSig(TString var="misc.PseudoJetMT2", TString cuts="misc.HBHENoiseFlag == 1", TString xtitle="MT2 [GeV]", 
		     int nbins=50, double min=0., double max=1., bool cleaned=false, int type=0 ); // 0: s/sqrt(b), 1: s/sqrt(s+b), 3:s/b
	void PrintCutFlow(int njets=-2, int nleps=0);
        void FillMonitor(Monitor *count, TString sname, TString type, TString cut, float weight);
        void abcd_MT2(TString var="misc.MinMetJetDPhi", TString basecut="misc.HBHENoiseFlag == 1", 
		      TString upper_cut="misc.MinMetJetDPhi<0.2", TString lower_cut="misc.MinMetJetDPhi>0.3", 
		      const int nbins=100, const double min=0., const double max=380., double fit_min=40., double fit_max=100.);
        void ABCD_MT2(TString var="misc.MinMetJetDPhi", TString basecut="misc.HBHENoiseFlag == 1", 
		      TString upper_cut="misc.MinMetJetDPhi<0.2", TString lower_cut="misc.MinMetJetDPhi>0.3", 
		      const int nbins=gNMT2bins, const double *bins=gMT2bins, double fit_min=40., double fit_max=100.);

	void PrintZllEfficiency(sample Sample, bool data, std::string lept );
	void CompSamples(TString var, TString cuts, TString optcut, bool RemoveLepts, TString xtitle, 
			  const int nbins, const double *bins, bool add_underflow, bool logflag, double scale_factor, bool normalize);
	void compSamples(TString var, TString cuts, TString optcut,  bool RemoveLepts, TString xtitle, 
			  const int nbins, const double min, const double max, bool add_underflow, bool logflag, double scale_factor, bool normalize);

private:

	TString fOutputDir;
	TFile *fOutputFile;
	int fVerbose;
	TString fPath;
	
	typedef std::map <TString, TString> MapType;
	MapType RemoveLeptMap;

	MT2tree* fMT2tree;
	TTree*   fTree;


        void MakeMT2PredictionAndPlots(bool cleaned , double dPhisplit[], double fudgefactor);
        void PrintABCDPredictions(TString var, TString basecut, TString upper_cut, TString lower_cut, TF1* func_qcd, TF1* func_sub);
        void MakePlot(std::vector<sample> Samples, TString var="misc.PseudoJetMT2", TString cuts="misc.HBHENoiseFlag == 1", 
		      int njets=-2, int nleps=0, //njets: 2 -> njets==2, -2 -> njets>=2
		      TString xtitle="MT2 [GeV]", const int nbins=50, const double min=0, const double max=1, 
		      bool cleaned=false, bool logflag=true, bool composited=false, bool ratio=false, 
		      bool stacked=true, bool overlaySUSY=false, float overlayScale = 0);
        void MakePlot(std::vector<sample> Samples, TString var="misc.PseudoJetMT2", TString cuts="misc.HBHENoiseFlag == 1", 
		      int njets=-2, int nleps=0, //njets: 2 -> njets==2, -2 -> njets>=2
		      TString xtitle="MT2 [GeV]", const int nbins=gNMT2bins, const double *bins=gMT2bins, 
		      bool cleaned=false, bool logflag=true, bool composited=false, bool ratio=false, 
		      bool stacked=true, bool overlaySUSY=false, float overlayScale = 0);
	void MakePlot(std::vector<sample> Samples, TString branch_name, const int nbins, const double bins[], bool cleaned, bool logflag, 
		      TString cut_branch_name, double ysplit[], TString option, TString version, TString prediction, double factor );
	void printHisto(THStack* h, TString canvname, Option_t *drawopt,  bool logflag);
	void printHisto(THStack* h, TH1* h_data, TLegend* leg,  TString canvname, Option_t *drawopt,  bool logflag, TString xtitle, TString ytitle);
        void printHisto(THStack* h, TH1* h_data, TH1* h_prediction, TLegend* leg,  TString canvname, Option_t *drawopt, bool logflag, TString xtitle, TString ytitle,int njets=-2, int nleps=0, bool stacked=true);
	void printHisto(THStack* h, TH1* h_data, TH1* h_prediction, TH1* h_susy, TLegend* leg,  TString canvname, Option_t *drawopt, bool logflag, TString xtitle, TString ytitle,int njets=-2, int nleps=0, float overlayScale=0);
	void printHisto(TH1* h, TString canvname, Option_t *drawopt, bool logflag);
	void plotRatio(TH1* h1, TH1* h2, bool logflag, bool normalize, TString name, TLegend* leg, TString xtitle, TString ytitle);
	void plotRatioStack(THStack* hstack, TH1* h1, TH1* h2, bool logflag, bool normalize, TString name, TLegend* leg, TString xtitle, TString ytitle,int njets=-2, int nleps=0);
	void plotRatioStack(THStack* hstack, TH1* h1, TH1* h2, TH1* h3, bool logflag, bool normalize, TString name, TLegend* leg, TString xtitle, TString ytitle,int njets=-2, int nleps=0, float overlayScale=0);
	void FillRatioHistdPhi(TString cut_branch_name, double lower_cut, double middle_cut, double upper_cut, TString version, bool cleaned);
	void FillRatioHistHT(TString version, bool cleaned);

	double ExpoFitTesting(double *x, double *par);
	TH1D* FillRatioHist(TString branch_name, float MT2plit[], float ysplit[], TString option, const int nbins, const double bins[], TString version,  bool cleaned);
	TH1D* GetPrediction(TString branch_name, const int nbins, const double bins[], bool cleaned, TString cut_branch_name, double lower_cut, double upper_cut, bool data, double factor );
	void  CompSamples(std::vector<sample> Samples, TString var, TString cuts, TString optcut, bool RemoveLepts, TString xtitle, 
			  const int nbins, const double *bins, bool add_underflow, bool logflag, double scale_factor, bool normalize);
	void  CompSamples(std::vector<sample> Samples, TString var, TString cuts, TString optcut,  bool RemoveLepts, TString xtitle, 
			  const int nbins, const double min, const double max, bool add_underflow, bool logflag, double scale_factor, bool normalize);
};

#endif
