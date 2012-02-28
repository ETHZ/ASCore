/*****************************************************************************
*   Collection of tools for producing plots for same-sign dilepton analysis  *
*                                                                            *
*                                                  (c) 2010 Benjamin Stieger *
*****************************************************************************/
#ifndef SSDLDumper_HH
#define SSDLDumper_HH

#include "helper/AnaClass.hh"
#include "helper/Monitor.hh"

#include "TLorentzVector.h"

using namespace std;

class SSDLDumper : public AnaClass{

public:
	// Binning
	static const int gNMuPtbins  = 8;
	static const int gNMuPt2bins = 8;
	static const int gNMuEtabins = 5;
	static const int gNElPtbins  = 7;
	static const int gNElPt2bins = 7;
	static const int gNElEtabins = 5;
	
	static const int gNNVrtxBins = 9;

	static double gNVrtxBins[gNNVrtxBins+1];

	static double gMuPtbins [gNMuPtbins+1];
	static double gMuPt2bins[gNMuPt2bins+1];
	static double gMuEtabins[gNMuEtabins+1];

	static double gElPtbins [gNElPtbins+1];
	static double gElPt2bins[gNElPt2bins+1];
	static double gElEtabins[gNElEtabins+1];

	static const int gNDiffHTBins   = 5;
	static const int gNDiffMETBins  = 4;
	static const int gNDiffMET3Bins = 9;
	static const int gNDiffNJBins   = 5;
	static const int gNDiffMT2Bins  = 3;
	static const int gNDiffPT1Bins  = 5;
	static const int gNDiffPT2Bins  = 5;
	static const int gNDiffNBJBins  = 4;

	static double gDiffHTBins  [gNDiffHTBins+1];
	static double gDiffMETBins [gNDiffMETBins+1];
	static double gDiffMET3Bins[gNDiffMET3Bins+1];
	static double gDiffNJBins  [gNDiffNJBins+1];
	static double gDiffMT2Bins [gNDiffMT2Bins+1];
	static double gDiffPT1Bins [gNDiffPT1Bins+1];
	static double gDiffPT2Bins [gNDiffPT2Bins+1];
	static double gDiffNBJBins [gNDiffNBJBins+1];

	static const int gM0bins  = 150  ; // these values
	static const int gM0min   = 0    ; // are 
	static const int gM0max   = 3000 ; // the same
	static const int gM12bins = 50   ; // as in
	static const int gM12min  = 0    ; // SSDLAnalysis
	static const int gM12max  = 1000 ; //

	// This enum has to correspond to the content of the samples.dat file
	enum gSample {
		sample_begin,
		DoubleMu1 = sample_begin, DoubleMu2, DoubleMu3, DoubleMu4, DoubleMu5,
		DoubleEle1, DoubleEle2, DoubleEle3, DoubleEle4, DoubleEle5,
		MuEG1, MuEG2, MuEG3, MuEG4, MuEG5,
		MuHad1, MuHad2, EleHad1, EleHad2,
		TTJets, TJets_t, TbarJets_t, TJets_tW, TbarJets_tW, TJets_s, TbarJets_s, WJets, DYJets,
		GJets40, GJets100, GJets200,
		GVJets, WGstarE, WGstarMu, WGstarTau, WW, WZ, ZZ, 
		TTbarW, TTbarZ, TTbarG, DPSWW, WWZ, WZZ, WWG, ZZZ, WWW, WpWp, WmWm,
		LM0, LM1, LM2, LM3, LM4, LM5, LM6, LM7, LM8, LM9, LM11, LM12, LM13, 
		QCDMuEnr10,
		QCD15, QCD30, QCD50, QCD80, QCD120, QCD170, QCD300, QCD470, QCD600, QCD800,
		QCD1000, QCD1400, QCD1800,
		gNSAMPLES
	};
	enum gHiLoSwitch{
		HighPt,
		LowPt,
		TauChan
	};
	enum gFPSwitch{
		SigSup,
		ZDecay
	};
	enum gRegion {
		region_begin,
		Baseline = region_begin,
		Control,
		HT80MET120,
		HT80MET120x,  // exclusive (HT < 200)
		HT200MET30,
		HT200MET120,
		HT200MET120x, // exclusive (HT < 450)
		HT450MET50,
		HT450MET50x,  // exclusive (MET < 120)
		HT450MET120,
		HT450MET0,
		HT0MET120,
		HT0MET200,
		HT0MET120JV,  // Jet veto
		HT0MET200JV,
		HT80MET302b,  // > 1 bjet
		HT200MET302b,
		HT80MET1202b,
		gNREGIONS
	};
	enum gChannel {
		channels_begin,
		Muon = channels_begin,
		ElMu,
		Elec,
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
		long nsst;
		long nssl;
		long nzl;
		long nzt;

		long nt2; // get these from signal events tree (overwrite the old ones, should agree)
		long nt10;
		long nt01;
		long nt0;
		
		float npp; // sum of event by event weights
		float npf;
		float nfp;
		float nff;
	};

	struct Channel{ // MM/EE/EM
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
		static float minHT    [gNREGIONS];
		static float maxHT    [gNREGIONS];
		static float minMet   [gNREGIONS];
		static float maxMet   [gNREGIONS];
		static int   minNjets [gNREGIONS];
		static int   minNbjets[gNREGIONS];
		
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
	
	static const int gNKinVars = 11;
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
	
	static const int gNSels = 2;
	struct IsoPlots{
		static TString sel_name[gNSels];
		static int nbins[gNSels];
		TH1D *hiso[gNSels];
		TH1D *hiso_pt[gNSels][gNMuPt2bins];
		TH1D *hiso_nv[gNSels][gNNVrtxBins];
	};

	struct IdPlots{
		static TString sel_name[gNSels];
		static int nbins[gNSels];
		TH1D *hhoe[gNSels];
		TH1D *hsiesie[gNSels];
		TH1D *hdphi[gNSels];
		TH1D *hdeta[gNSels];
		//TH1D *hid_pt[gNSels][gNMuPt2bins];
		//TH1D *hid_nv[gNSels][gNNVrtxBins];
	};
	
	static const int gNDiffVars = 10;
	struct DiffPredYields{
		static TString var_name[gNDiffVars];
		static TString axis_label[gNDiffVars];
		static int nbins[gNDiffVars];
		static double* bins[gNDiffVars];

		TH1D *hnt11[gNDiffVars]; // SS yields
		TH1D *hnt10[gNDiffVars];
		TH1D *hnt01[gNDiffVars];
		TH1D *hnt00[gNDiffVars];

		TH1D *hnpp[gNDiffVars]; // EbE predictions
		TH1D *hnpf[gNDiffVars];
		TH1D *hnfp[gNDiffVars];
		TH1D *hnff[gNDiffVars];

		TH1D *hnt2_os_BB[gNDiffVars]; // OS yields
		TH1D *hnt2_os_EB[gNDiffVars];
		TH1D *hnt2_os_EE[gNDiffVars];
	};
	
	static const int gNKinSels = 3;
	static TString gKinSelNames[gNKinSels];
	static TString gEMULabel[2];
	static TString gChanLabel[gNCHANNELS];
	static TString gHiLoLabel[3];

	class Sample{
	public:
		Sample(){};
		Sample(TString loc, TString tag, int dm, int cs = -1, float xs = 1.0, long ng = 1, int col = 1){
			location = loc;
			sname    = tag;
			datamc   = dm;
			chansel  = cs;
			xsec     = xs;
			ngen     = ng;
			color    = col;
		};
		~Sample(){};
		
		TString name;
		TString sname;
		TString location;
		TFile *file;
		TTree *tree;
		float lumi; // simulated lumi = ngen/xsec
		float xsec; // cross-section
		int ngen;   // generated number of events
		int color;
		int datamc;  // 0: Data, 1: SM MC, 2: Signal MC, 3: rare MC, 4: rare MC (no pileup)
		int chansel; // -1: Ignore, 0: mumu, 1: elel, 2: elmu
		int proc;    // process type, to group binned samples like QCD and give latex namees
		TH1D *cutFlowHisto[gNCHANNELS];
		Region region[gNREGIONS][2];
		DiffPredYields diffyields[gNCHANNELS];
		NumberSet numbers[gNREGIONS][gNCHANNELS]; // summary of integrated numbers
		KinPlots kinplots[gNKinSels][2]; // tt and ll and signal for both low and high pt analysis
		IsoPlots isoplots[2]; // e and mu
		IdPlots  idplots; // only for electrons
		FRatioPlots ratioplots[2]; // e and mu
		TGraph *sigevents[gNCHANNELS][2];

		float getLumi(){
			if(ngen > 0 && xsec > 0) return float(ngen)/xsec;
			else return -1.;
		}

		float getError(int n){
			// If n passed of ngen generated, what is upper limit
			// on number of events passing?
			if(ngen <= 0) return -1.;
			TEfficiency *eff = new TEfficiency();
			float upper = eff->ClopperPearson(ngen, n, 0.68, true);
			float delta = upper - float(n)/float(ngen);
			delete eff;
			return delta * float(ngen);
		}

		int getType(){ // -1: undef, 0: data, 1: QCD, 2: top, 3: EWK, 4: Rare SM, 5: diboson
			if(datamc == 0) return 0;
			if( (sname.Contains("QCD")) ||
			    (sname) == "MuEnr10" ) return 1;
			if( (sname.Contains("SingleT")) ||
			    (sname) == "TTJets" )  return 2;
			if( (sname.Contains("DYJets")) ||
			    (sname.Contains("GJets"))  ||
			    (sname) == "WJets" )   return 3;
			if( (sname) == "TTbarW"    ||
			    (sname) == "TTbarZ"    ||
			    (sname) == "TTbarG"    ||
			    (sname) == "DPSWW"     ||
			    (sname) == "WWZ"       ||
			    (sname) == "WZZ"       ||
			    (sname) == "WZZ"       ||
			    (sname) == "WWG"       ||
			    (sname) == "ZZZ"       ||
			    (sname) == "WWW"       ||
			    (sname) == "W+W+"      ||
			    (sname) == "W-W-")     return 4;
			if( (sname.Contains("GVJets"))    ||
			    (sname.Contains("WWTo2L2Nu")) ||
			    (sname.Contains("WZTo3LNu"))  ||
			    (sname.Contains("ZZTo4L")) ) return 5;
			else {
				cout << "SSDLDumper::Sample::getType() ==> ERROR: "<< sname << " has no defined type!" << endl;
				return -1;
			}
		}
		int getProc(){ // used for binned samples
			if(datamc == 0)                                 return 0;
			if(sname == "TTJets" )                          return 1;
			if(sname.Contains("SingleT"))                   return 2;
			if(sname == "WJets" )                           return 3;
			if(sname.Contains("DYJets"))                    return 4;
			if(sname.Contains("GJets"))                     return 5;
			if(sname.Contains("WWTo2L2Nu"))                 return 6; 
			if(sname.Contains("WZTo3LNu"))                  return 7; 
			if(sname.Contains("ZZTo4L"))                    return 8; 
			if(sname.Contains("GVJets"))                    return 9; 
			if(sname == "TTbarW")                           return 10;
			if(sname == "TTbarZ")                           return 11;
			if(sname == "TTbarG")                           return 12;
			if(sname == "W+W+" || sname == "W-W-")          return 13;
			if(sname == "WWZ" ||                                     
			   sname == "WZZ" ||                                     
			   sname == "WZZ" ||                                     
			   sname == "WWG" ||                                     
			   sname == "ZZZ" ||                                     
			   sname == "WWW" )                             return 14;
			if(sname == "DPSWW")                            return 15;
			if(sname.Contains("QCD") || sname == "MuEnr10") return 16;
			else {
				cout << "SSDLDumper::Sample::getProc() ==> ERROR: "<< sname << " has no defined process!" << endl;
				return -1;
			}
		}
		TString getProcName(int proc){ // used for binned samples
			if(proc ==  0) return "Data";
			if(proc ==  1) return "t\\bar{t}";
			if(proc ==  2) return "Single t";
			if(proc ==  3) return "W+jets";
			if(proc ==  4) return "Z+jets";
			if(proc ==  5) return "\\gamma+jets";
			if(proc ==  6) return "WW";
			if(proc ==  7) return "WZ";
			if(proc ==  8) return "ZZ";
			if(proc ==  9) return "V\\gamma+jets";
			if(proc == 10) return "t\\bar{t}W";
			if(proc == 11) return "t\\bar{t}Z";
			if(proc == 12) return "t\\bar{t}\\gamma";
			if(proc == 13) return "W^{\\pm}W^{\\pm}";
			if(proc == 14) return "Tri-Boson";
			if(proc == 15) return "DPS (2\\times W+jets)";
			if(proc == 16) return "QCD";
			else {
				cout << "SSDLDumper::Sample::getProcName() ==> ERROR: "<< proc << " has no defined process name!" << endl;
				return "";
			}
		}
		inline int getNProcs(){return 17;} // make sure this number corresponds to the number of
		                                   // processes define in the previous method
		TTree* getTree(){
			file = TFile::Open(location);
			if(file->IsZombie()){
				cout << "SSDLDumper::Sample::getTree() ==> Error opening file " << location << endl;
				exit(1);
			}
			tree = (TTree*)file->Get("Analysis");
			return tree;
		}
		
		void cleanUp(){
			tree = NULL;
			file->Close();
		}
	};
	
	SSDLDumper();
	virtual ~SSDLDumper();

	virtual void init();
	virtual void init(TString); // samples read from a datacard
	virtual void init(TString, TString, int, int = -1, float = 1.0); // running on single sample
	virtual void InitMC(TTree*); // removed a few branches (triggers)

	virtual void readDatacard(TString); // read in a datacard

	virtual void loop();              // loop on all samples if several
	virtual void loopEvents(Sample*); // perform loop on single sample
	
	void storeNumbers(Sample*, gChannel, gRegion);
	
	// Cutflow
	void initCutNames();
	void initCounters();
	void fillCutFlowHistos(Sample*);
	void printCutFlow(gChannel, gSample, gSample);
	void printCutFlows(TString);
	
	//////////////////////////////
	// Fillers
	void fillSigEventTree(Sample*);
	void resetSigEventTree();
	void fillYields(Sample*, gRegion, gHiLoSwitch);
	void fillDiffYields(Sample*);
	void fillRatioPlots(Sample*);
	void fillMuIsoPlots(Sample*);
	void fillElIsoPlots(Sample*);
	void fillElIdPlots (Sample*);
	void fillKinPlots(Sample*, gHiLoSwitch);
	
	//////////////////////////////
	// I/O
	void bookSigEvTree();
	void writeSigEvTree(TFile*);
	int getSampleType(Sample*);

	void bookHistos(Sample*);
	void deleteHistos(Sample*);
	void writeHistos(Sample*, TFile*);
	void writeSigGraphs(Sample*, gChannel, gHiLoSwitch, TFile*);
	int readHistos(TString);
	int readSigGraphs(TString);
	
	void bookRatioHistos();

	// Geninfo stuff:
	int muIndexToBin(int);
	int elIndexToBin(int);
	TString muBinToLabel(int);
	TString elBinToLabel(int);
	void labelOriginAxis(TAxis*, gChannel);
	void label2OriginAxes(TAxis*, TAxis*, gChannel);
	
	// Trigger selections:
	virtual bool  singleMuTrigger();
	virtual float singleMuPrescale();
	virtual bool  singleElTrigger();
	virtual float singleElPrescale();

	virtual bool mumuSignalTrigger();
	virtual bool elelSignalTrigger();
	virtual bool elmuSignalTrigger();

	virtual bool doubleMuTrigger();
	virtual bool doubleElTrigger();
	virtual bool eMuTrigger();

	virtual bool doubleMuHTTrigger();
	virtual bool doubleElHTTrigger();
	virtual bool eMuHTTrigger();
	
	virtual bool isGoodRun(Sample*);

	// Event and Object selectors:
	virtual int getNJets();
	virtual int getNBTags();
	virtual float getHT();
	virtual float getMT2(int, int, gChannel);
	virtual float getMll(int, int, gChannel);
	virtual int   getClosestJet(int, gChannel);
	virtual float getClosestJetPt(int, gChannel);
	virtual float getClosestJetDR(int, gChannel);
	virtual float getSecondClosestJetDR(int, gChannel);
	virtual float getAwayJetPt(int, gChannel);
	virtual float getMaxJPt();
	
	virtual int isSSLLEvent(int&, int&);
	virtual int isOSLLEvent(int&, int&);
	virtual int isSSEvent(int&, bool(SSDLDumper::*)(int), int&, bool(SSDLDumper::*)(int));
	virtual int isOSEvent(int&, bool(SSDLDumper::*)(int), int&, bool(SSDLDumper::*)(int));
	std::vector<lepton> sortLeptonsByPt(std::vector<lepton> &leptons);
	vector<lepton> sortLeptonsByTypeAndPt(vector<lepton> &leptons);

	virtual bool isGoodEvent();
	virtual bool isGoodMuEvent();
	virtual int hasLooseMuons(int&, int&);
	virtual int hasLooseMuons();
	virtual int hasLooseElectrons(int&, int&);
	virtual int hasLooseElectrons();
	virtual bool passesNJetCut(int=2);
	virtual bool passesJet50Cut();
	virtual bool passesNbJetCut(int=2);
	
	virtual bool passesHTCut(float, float = 7000.);
	virtual bool passesMETCut(float = -1., float = 7000.);
	virtual bool passesZVeto(bool(SSDLDumper::*)(int), bool(SSDLDumper::*)(int), float = 15.); // cut with mZ +/- cut value and specified obj selectors
	virtual bool passesZVeto(float = 15.); // cut with mZ +/- cut value
	virtual bool passesMllEventVeto(float = 5.);
	virtual bool passesMllEventVeto(int, int, int, float = 5.);

	virtual bool isSigSupMuEvent();
	virtual bool isZMuMuEvent(int&, int&);

	virtual bool isSigSupElEvent();
	virtual bool isZElElEvent(int&, int&);

	virtual bool isGenMatchedSUSYDiLepEvent();
	virtual bool isGenMatchedSUSYDiLepEvent(int&, int&);

	virtual bool isGenMatchedSUSYEEEvent();
	virtual bool isGenMatchedSUSYEMuEvent();

	virtual bool isSSLLMuEvent(int&, int&);
	virtual bool isSSLLElEvent(int&, int&);
	virtual bool isSSLLElMuEvent(int&, int&);

	virtual bool isGoodMuon(int, float = -1.);
	virtual bool isLooseMuon(int);
	virtual bool isTightMuon(int);
	virtual bool isGoodPrimMuon(int, float = -1.);
	virtual bool isGoodSecMuon(int, float = -1.);

	virtual bool isFakeMuon(int);
	virtual bool isPromptMuon(int);
	virtual bool isChargeMatchedMuon(int);

	virtual bool isGoodElectron(int, float = -1.);
	virtual bool isLooseElectron(int);
	virtual bool isTightElectron(int);
	virtual bool isGoodPrimElectron(int, float = -1.);
	virtual bool isGoodSecElectron(int, float = -1.);

	virtual bool isFakeElectron(int);
	virtual bool isPromptElectron(int);
	virtual bool isChargeMatchedElectron(int);

	virtual bool isBarrelElectron(int);

	virtual bool isGoodJet(int, float = 40.);

	float fC_minHT;
	float fC_minMet;
	float fC_maxHT;
	float fC_maxMet;
	int   fC_minNjets;
	int   fC_minNbjets;
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
	
	const double *getDiffPredBins(int);
	
	Monitor fCounters[gNSAMPLES][3];
	vector<string> fMMCutNames;
	vector<string> fEMCutNames;
	vector<string> fEECutNames;
	bool fDoCounting;
	gSample fCurrentSample;
	gChannel fCurrentChannel;
	ofstream fOUTSTREAM, fOUTSTREAM2, fOUTSTREAM3;

	int fChargeSwitch;    // 0 for SS, 1 for OS

	float fLumiNorm;      // Normalize everything to this luminosity

	vector<Sample*>::iterator fS;
	vector<Sample*> fSamples;
	vector<Sample*> fMCSamples;
	map<TString, Sample*> fSampleMap;	// Mapping of sample to name
	// map<TString, gSample> fSampleMap;	// Mapping of sample number to name
	
	TTree *fSigEv_Tree;
	float       fSETree_PUWeight;
	float       fSETree_SLumi;
	std::string fSETree_SName;
	int         fSETree_SType; // DoubleMu(0), DoubleEle(1), MuEG(2), MC(10)
	int         fSETree_Run;
	int         fSETree_LS;
	long        fSETree_Event;
	int         fSETree_Flavor; // mm(0), em(1), ee(2)
	int         fSETree_Charge;
	int         fSETree_TLCat; // TL category: TT(0), TL(1), LT(2), LL(3)
	float       fSETree_HT;
	float       fSETree_MET;
	int         fSETree_NJ;
	int         fSETree_NbJ;
	float       fSETree_MT2;
	float       fSETree_Mll;
	float       fSETree_pT1;
	float       fSETree_pT2;

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

	Sample *fSample;

	private:
	
	Monitor fCounter[3];	
};

#endif
