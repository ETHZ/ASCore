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
	

	// probability for W->lnu event to be reconstructed and identified. 
	struct Wprediction{
		// efficiencies for W->lnu
		double Wenu_acc;
		double Wenu_acc_err;
		double Wenu_rec;
		double Wenu_rec_err;
		double Wenu_prob;
		double Wenu_prob_err;
		double Wmunu_acc;
		double Wmunu_acc_err;
		double Wmunu_rec;
		double Wmunu_rec_err;
		double Wmunu_prob;
		double Wmunu_prob_err;
		double W_prob(std::string lept){
			if     (lept == "ele" && Wenu_acc  >0 && Wenu_rec  >0){return Wenu_acc *Wenu_rec;}
			else if(lept == "muo" && Wmunu_acc >0 && Wmunu_rec >0){return Wmunu_acc*Wmunu_rec;}
			else return -1;
		}
		double W_prob_err(std::string lept){
			if     (lept == "ele" && W_prob("ele")>0){return sqrt(pow(Wenu_acc *Wenu_rec_err ,2) + pow(Wenu_rec *Wenu_acc_err ,2));}
			else if(lept == "muo" && W_prob("muo")>0){return sqrt(pow(Wmunu_acc*Wmunu_rec_err,2) + pow(Wmunu_rec*Wmunu_acc_err,2));}
			else return -1;
		}
		// efficiencies for Top W->lnu
		double TopWenu_acc;
		double TopWenu_acc_err;
		double TopWenu_rec;
		double TopWenu_rec_err;
		double TopWenu_prob;
		double TopWenu_prob_err;
		double TopWmunu_acc;
		double TopWmunu_acc_err;
		double TopWmunu_rec;
		double TopWmunu_rec_err;
		double TopWmunu_prob;
		double TopWmunu_prob_err;
		double TopW_prob(std::string lept){
			if     (lept == "ele" && TopWenu_acc  >0 && TopWenu_rec  >0){return TopWenu_acc *TopWenu_rec;}
			else if(lept == "muo" && TopWmunu_acc >0 && TopWmunu_rec >0){return TopWmunu_acc*TopWmunu_rec;}
			else return -1;
		}
		double TopW_prob_err(std::string lept){
			if     (lept == "ele" && TopW_prob("ele")>0){return sqrt(pow(TopWenu_acc *TopWenu_rec_err ,2) + pow(TopWenu_rec *TopWenu_acc_err ,2));}
			else if(lept == "muo" && TopW_prob("muo")>0){return sqrt(pow(TopWmunu_acc*TopWmunu_rec_err,2) + pow(TopWmunu_rec*TopWmunu_acc_err,2));}
			else return -1;
		}
		double QCD_bg_e;
		double QCD_bg_mu;
		double Top_bg_e;
		double Top_bg_mu;
		double W_bg_e;
		double W_bg_mu;
		double Z_bg_e;
		double Z_bg_mu;
		double Other_bg_e;
		double Other_bg_mu;
	} fWpred;

	struct Zprediction{
		double ele_reco;
		double ele_reco_err;
		double muo_reco;
		double muo_reco_err;
		double nu_acc;
		double nu_acc_err;
		double R(std::string lept){
			if     (lept=="ele" && ele_reco>0 && nu_acc > 0){return 1./(ele_reco*nu_acc);}
			else if(lept=="muo" && muo_reco>0 && nu_acc > 0){return 1./(muo_reco*nu_acc);}
			else return -1;
		}
		double R_err(std::string lept){
			// error propagation
			if     (lept=="ele" && R("ele")!=-1){return R("ele")*sqrt(pow(ele_reco_err/ele_reco,2)+pow(nu_acc_err/nu_acc,2));}
			else if(lept=="muo" && R("muo")!=-1){return R("muo")*sqrt(pow(muo_reco_err/muo_reco,2)+pow(nu_acc_err/nu_acc,2));}
			else return -1;
		}
		double QCD_bg_e;
		double QCD_bg_mu;
		double Top_bg_e;
		double Top_bg_mu;
		double W_bg_e;
		double W_bg_mu;
		double Other_bg_e;
		double Other_bg_mu;
	} fZpred;
	// -----------------
	

	typedef std::map <TString, TString> MapType;
	MapType RemoveLeptMap;

	void setVerbose(int v){ fVerbose = v;};
	void setOutputDir(TString dir){ fOutputDir = Util::MakeOutputDir(dir); };
	void setOutputFile(TString filename){ fOutputFile = Util::MakeOutputFile(fOutputDir + filename); };
        void makePlot(TString var="misc.PseudoJetMT2", TString cuts="misc.HBHENoiseFlag == 1", 
		      int njets=-2, int nleps=0, TString HLT="", //njets: 2 -> njets==2, -2 -> njets>=2
		      TString xtitle="MT2 [GeV]", const int nbins=50, const double min=0., const double max=1., 
		      bool cleaned=false, bool logflag=true, bool composited=false, bool ratio=false, 
		      bool stacked=true, bool overlaySUSY=false, float overlayScale = 0);
        void MakePlot(TString var="misc.PseudoJetMT2", TString cuts="misc.HBHENoiseFlag == 1", 
		      int njets=-2, int nleps=0, TString HLT="",//njets: 2 -> njets==2, -2 -> njets>=2
		      TString xtitle="MT2 [GeV]", const int nbins=gNMT2bins, const double *bins=gMT2bins, 
		      bool cleaned=false, bool logflag=true, bool composited=false, bool ratio=false, 
		      bool stacked=true, bool overlaySUSY=false, float overlayScale = 0);
        void plotSig(TString var="misc.PseudoJetMT2", TString cuts="misc.HBHENoiseFlag == 1", TString xtitle="MT2 [GeV]", 
		     int nbins=50, double min=0., double max=1., bool cleaned=false, int type=0 ); // 0: s/sqrt(b), 1: s/sqrt(s+b), 3:s/b
  	void PrintCutFlow(int njets=-2, int nleps=0, TString trigger="", TString cuts="");
        void FillMonitor(Monitor *count, TString sname, TString type, TString cut, double weight);
	void PrintZllEfficiency(int sample_index, bool data, std::string lept, Long64_t nevents, double lower_mass, double upper_mass, bool pileup_weight);
	void PrintWEfficiency(int sample_index ,TString process, std::string lept, Long64_t nevents, bool includeTaus);
        void abcd_MT2(TString var="misc.MinMetJetDPhi", TString basecut="misc.HBHENoiseFlag == 1", 
		      TString upper_cut="misc.MinMetJetDPhi<0.2", TString lower_cut="misc.MinMetJetDPhi>0.3", 
		      const int nbins=100, const double min=0., const double max=380., double fit_min=40., double fit_max=100.);
        void ABCD_MT2(TString var="misc.MinMetJetDPhi", TString basecut="misc.HBHENoiseFlag == 1", 
		      TString upper_cut="misc.MinMetJetDPhi<0.2", TString lower_cut="misc.MinMetJetDPhi>0.3", 
		      const int nbins=gNMT2bins, const double *bins=gMT2bins, double fit_min=40., double fit_max=100.);

	void CompSamples(TString var, TString cuts, TString optcut, bool RemoveLepts, TString xtitle, 
			  const int nbins, const double *bins, bool add_underflow, bool logflag, double scale_factor, bool normalize);
	void compSamples(TString var, TString cuts, TString optcut,  bool RemoveLepts, TString xtitle, 
			  const int nbins, const double min, const double max, bool add_underflow, bool logflag, double scale_factor, bool normalize);

private:

	TString fOutputDir;
	TFile *fOutputFile;
	int fVerbose;
	TString fPath;

	MT2tree* fMT2tree;
	TTree*   fTree;


        void MakeMT2PredictionAndPlots(bool cleaned , double dPhisplit[], double fudgefactor);
        void PrintABCDPredictions(TString var, TString basecut, TString upper_cut, TString lower_cut, TF1* func_qcd, TF1* func_sub, TF1* func_qcd_model);
        void printEstimation(TH1D* h_pred, TH1D* h_pred_c, int nbins, float min, float max);
        void MakePlot(std::vector<sample> Samples, TString var="misc.PseudoJetMT2", TString cuts="misc.HBHENoiseFlag == 1", 
		      int njets=-2, int nleps=0, TString HLT="", //njets: 2 -> njets==2, -2 -> njets>=2
		      TString xtitle="MT2 [GeV]", const int nbins=50, const double min=0, const double max=1, 
		      bool cleaned=false, bool logflag=true, bool composited=false, bool ratio=false, 
		      bool stacked=true, bool overlaySUSY=false, float overlayScale = 0);
        void MakePlot(std::vector<sample> Samples, TString var="misc.PseudoJetMT2", TString cuts="misc.HBHENoiseFlag == 1", 
		      int njets=-2, int nleps=0, TString HLT="", //njets: 2 -> njets==2, -2 -> njets>=2
		      TString xtitle="MT2 [GeV]", const int nbins=gNMT2bins, const double *bins=gMT2bins, 
		      bool cleaned=false, bool logflag=true, bool composited=false, bool ratio=false, 
		      bool stacked=true, bool overlaySUSY=false, float overlayScale = 0);
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
	void  CompSamples(std::vector<sample> Samples, TString var, TString cuts, TString optcut, bool RemoveLepts, TString xtitle, 
			  const int nbins, const double *bins, bool add_underflow, bool logflag, double scale_factor, bool normalize);
	void  CompSamples(std::vector<sample> Samples, TString var, TString cuts, TString optcut,  bool RemoveLepts, TString xtitle, 
			  const int nbins, const double min, const double max, bool add_underflow, bool logflag, double scale_factor, bool normalize);
};

#endif
