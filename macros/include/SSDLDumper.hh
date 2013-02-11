/*****************************************************************************
*   Collection of tools for producing plots for same-sign dilepton analysis  *
*                                                                            *
*                                                  (c) 2010 Benjamin Stieger *
*****************************************************************************/
#ifndef SSDLDumper_HH
#define SSDLDumper_HH

#include "helper/AnaClass.hh"
#include "helper/Monitor.hh"
//#include "helper/BTagSFUtil/BTagSFUtil.h"
#include "helper/BTagSF.hh"
#include "helper/GoodRunList.h"

#include "TRandom3.h"
#include "TLorentzVector.h"

using namespace std;

class SSDLDumper : public AnaClass{

public:
	// Binning
	static const int gNMuFPtBins = 6;
	static const int gNMuPPtbins = 10;
	static const int gNMuEtabins = 5;
	static const int gNElFPtBins = 8;
	static const int gNElPPtbins = 10;
	static const int gNElEtabins = 5;
        static const int gNElCMIdbins = 2;

        static const int gNNVrtxBins = 10;

	static double gNVrtxBins[gNNVrtxBins+1];

	static double gMuPPtbins[gNMuPPtbins+1];
	static double gMuFPtBins[gNMuFPtBins+1];
	static double gMuEtabins[gNMuEtabins+1];

	static double gElPPtbins[gNElPPtbins+1];
	static double gElFPtBins[gNElFPtBins+1];
	static double gElEtabins[gNElEtabins+1];
        static double gElCMIdbins[gNElCMIdbins+1];

	static const int gNDiffHTBins   = 6;
	static const int gNDiffMETBins  = 6;
	static const int gNDiffMET3Bins = 9;
	static const int gNDiffNJBins   = 6;
	static const int gNDiffMT2Bins  = 3;
	static const int gNDiffPT1Bins  = 9;
	static const int gNDiffPT2Bins  = 8;
	static const int gNDiffNBJBins  = 4;
	static const int gNDiffNBJMBins = 3;

	static double gDiffHTBins  [gNDiffHTBins+1];
	static double gDiffMETBins [gNDiffMETBins+1];
	static double gDiffMET3Bins[gNDiffMET3Bins+1];
	static double gDiffNJBins  [gNDiffNJBins+1];
	static double gDiffMT2Bins [gNDiffMT2Bins+1];
	static double gDiffPT1Bins [gNDiffPT1Bins+1];
	static double gDiffPT2Bins [gNDiffPT2Bins+1];
	static double gDiffNBJBins [gNDiffNBJBins+1];
	static double gDiffNBJMBins[gNDiffNBJMBins+1];

	static const int gM0bins  = 150  ; // these values
	static const int gM0min   = 0    ; // are
	static const int gM0max   = 3000 ; // the same
	static const int gM12bins = 50   ; // as in
	static const int gM12min  = 0    ; // SSDLAnalysis
	static const int gM12max  = 1000 ; //

	float gHLTSF;
	float gBtagSF;
	float gBtagSF1;
	float gBtagSF2;
	float gEventWeight;
	float gEventWeight0;
	float gEventWeight1;
	float gEventWeight2;

	// This enum has to correspond to the content of the samples.dat file
	enum gSample {
		sample_begin,
		// data samples
		DoubleMu1 = sample_begin, DoubleMu1a , DoubleMu2 , DoubleMu3 , DoubleMu4 , DoubleMu5 , DoubleMu5a ,
		DoubleEle1              , DoubleEle1a, DoubleEle2, DoubleEle3, DoubleEle4, DoubleEle5, DoubleEle5a,
		MuEG1                   , MuEG1a     , MuEG2     , MuEG3     , MuEG4     , MuEG5     , MuEG5a     ,
		// fake samples
		TTJets, SingleT_t, SingleTbar_t, SingleT_tW, SingleTbar_tW, SingleT_s, SingleTbar_s,
		WJets,
		DYJets,
		GJets200, GJets400, WW,
		// start of the rares
		WZ,ZZ,
		TTbarH, TTbarW, TTbarZ, TTbarG, TbZ, DPSWW,
		WWZ, WZZ, 
		WWG, ZZZ, WWW,
		TTbarWW,
		WpWp, WmWm,
		// QCD samples
		QCDMuEnr15,
		QCD80, QCD120, QCD170, QCD300, QCD470, QCD600, QCD800,
		QCDEM20, QCDEM30, QCDEM80, QCDEM170, QCDEM250,
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
	// enum gRegion {
	// 	region_begin,
	// 	Baseline = region_begin,
	// 	HT80MET30, HT80MET30b, HT80MET30bpp,
	// 	HT0MET120,HT0MET80NJ2,HT0MET80NJ2bV,
	// 	HT0MET120V, HT0MET120NJ2, HT0MET120NJ2bV,	// EWino regions MET > 120
	// 	HT0MET200, HT0MET120NJ2bVlV, HT0MET200lV,   
	// 	TTbarWPresel, TTbarWSel, //TTbarV selections
	// 	HT0MET50, HT0MET80, HT0MET50lV, HT0MET80lV, HT0MET120lV,
	// 	gNREGIONS
	// };
	enum gChannel {
		channels_begin,
		Muon = channels_begin,
		ElMu,
		Elec,
		gNCHANNELS
	};
        enum gChMisIdReg {
	  chmisidreg_begin,
	  BB = chmisidreg_begin,
	  EB,
	  EE
	};
	struct ValueAndError {
		float val;
		float err;
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
		
		float tt_avweight;
		float tl_avweight;
		float lt_avweight;
		float ll_avweight;
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

		TH1D *fntight_nv; // nvertices
		TH1D *fnloose_nv;
		TH1D *pntight_nv;
		TH1D *pnloose_nv;

		// duplicate for only ttbar use
		TH2D *fntight_ttbar; // pt vs eta
		TH2D *fnloose_ttbar; 
		TH2D *pntight_ttbar; // pt vs eta
		TH2D *pnloose_ttbar;
		
		// gen ID
		TH1D *fntight_genID;
		TH1D *fnloose_genID;
		TH1D *pntight_genID;
		TH1D *pnloose_genID;
		
		TEfficiency *fratio_pt;
		TEfficiency *pratio_pt;
		TEfficiency *fratio_eta;
		TEfficiency *pratio_eta;
		TEfficiency *fratio_nv;
		TEfficiency *pratio_nv;
		// Charged miss-id calculation
		TH2D *ospairs;
		TH2D *sspairs;

                // Charged miss-id calculation
	        TEfficiency *chmid_BB_pt;
	        TEfficiency *chmid_EE_pt;
    	        TEfficiency *chmid_BE_pt;

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
		TH2D *nt10_OS_BB_pt;
		TH2D *nt10_OS_EE_pt;
		TH2D *nt10_OS_EB_pt;
		TH2D *nt01_OS_BB_pt;
		TH2D *nt01_OS_EE_pt;
		TH2D *nt01_OS_EB_pt;
		TH2D *nt00_OS_BB_pt;
		TH2D *nt00_OS_EE_pt;
		TH2D *nt00_OS_EB_pt;
	};
	
	struct Region{
		TString sname;
		// Two different pt cuts
		float minMu1pt;
		float minMu2pt;
		float minEl1pt;
		float minEl2pt;
		// Custom selections for every region
		float minHT    ;
		float maxHT    ;
		float minMet   ;
		float maxMet   ;
		float minJetPt ;
		int   minNjets ;
		int   maxNjets ;
		int   minNbjets;
		int   maxNbjets;
		int   minNbjmed;
		int   maxNbjmed;
		int   app3rdVet;
		int   vetoTTZSel;
		int   chargeVeto;
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
	
	static const int gNKinVars = 12;
	struct KinPlots{
		static TString var_name[gNKinVars];
		static TString axis_label[gNKinVars];
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
		TH1D *hiso_pt[gNSels][gNMuFPtBins];
		TH1D *hiso_nv[gNSels][gNNVrtxBins];
	};

	struct IdPlots{
		static TString sel_name[gNSels];
		static int nbins[gNSels];
		TH1D *hhoe[gNSels];
		TH1D *hsiesie[gNSels];
		TH1D *hdphi[gNSels];
		TH1D *hdeta[gNSels];
		TH1D *hmvaid[gNSels];
		TH1D *hmedwp[gNSels];
		//TH1D *hid_pt[gNSels][gNMuFPtBins];
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
	std::vector< SSDLDumper::Region* > gRegions;
	std::vector< SSDLDumper::Region* >::iterator regIt;
	std::map<TString , int> gRegion;
	int gNREGIONS;
	TString gBaseRegion;
	std::map<TString , int> gSystematics;
	std::map<TString , int>::const_iterator gsystIt;
	bool gApplyZVeto;
    bool    gDoWZValidation;

	class Sample{
	public:
		Sample(){};
		Sample(TString loc, TString tag, int dm, float xs, int nregions, int cs = -1, int col = 1){
			location = loc;
			sname    = tag;
			datamc   = dm;
			chansel  = cs;
			xsec     = xs;
			color    = col;
			numbers = new NumberSet*[nregions];
			region  = new Region   *[nregions];
            for ( size_t i = 0; i<nregions; ++i ) region[i]  = new Region   [2];
            for ( size_t i = 0; i<nregions; ++i ) numbers[i] = new NumberSet[gNCHANNELS];
            //region[0] = new Region[gNREGIONS*gNCHANNELS];
            file = NULL;
            tree = NULL;
		};
		~Sample(){
			delete[] region[0];
			delete[] region;
			delete[] numbers[0];
			delete[] numbers;
        };
		
		TString name;
		TString sname;
		TString location;
		TFile *file;
		TTree *tree;
		float lumi; // simulated lumi = ngen/xsec
		float xsec; // cross-section
		int ngen;   // number of generated events
		TH1F *evcount; // count number of generated events
		int color;
		int datamc;  // 0: Data, 1: SM MC, 2: Signal MC, 3: rare MC, 4: rare MC (no pileup)
		int chansel; // -1: Ignore, 0: mumu, 1: elel, 2: elmu
		int proc;    // process type, to group binned samples like QCD and give latex namees
		TH1D *cutFlowHisto[gNCHANNELS];
//		Region region[gNREGIONS][2];
		Region **region;
		DiffPredYields diffyields[gNCHANNELS];
		NumberSet **numbers; // summary of integrated numbers
//		NumberSet numbers[gNREGIONS][gNCHANNELS]; // summary of integrated numbers
		KinPlots kinplots[gNKinSels][2]; // tt and ll and signal for both low and high pt analysis
	        KinPlots kinplots_wz[gNKinSels];
	        IsoPlots isoplots[2]; // e and mu
 	        IdPlots  idplots; // only for electrons
		FRatioPlots ratioplots[2]; // e and mu
		TGraph *sigevents[gNCHANNELS][2];

		float getLumi(){
			if(datamc == 0) return -1;
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
		float getError2(int n){
			float err = getError(n);
			return err*err;
		}

		int getType(){ // -1: undef, 0: data, 1: QCD, 2: top, 3: EWK, 4: Rare SM, 5: diboson
			if(datamc == 0) return 0;
			if( (sname.Contains("QCD")) ||
			    (sname) == "MuEnr15"    ||
			    (sname) == "MuEnr10"    ||
			    (sname) == "EMEnr20"    ||
			    (sname) == "EMEnr30"    ) return 1;
			if( (sname.Contains("SingleT")) ||
			    (sname.Contains("TTJets"))  ) return 2;
			if( (sname.Contains("DYJets")) ||
			    (sname.Contains("GJets"))  ||
			    (sname) == "WJets" )   return 3;
			if( (sname) == "TTbarH"    ||
			    (sname) == "TTbarW"    ||
			    (sname) == "TTbarZ"    ||
			    (sname) == "TTbarG"    ||
			    (sname) == "TbZ"       ||
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
			if(sname == "TTbarH")                           return 10;
			if(sname == "TTbarW")                           return 11;
			if(sname == "TTbarZ")                           return 12;
			if(sname == "TTbarG")                           return 13;
			if(sname == "W+W+" || sname == "W-W-")          return 14;
			if(sname == "WWZ" ||                                     
			   sname == "WZZ" ||                                     
			   sname == "WWG" ||                                     
			   sname == "ZZZ" ||                                     
			   sname == "WWW" )                             return 15;
			if(sname == "DPSWW")                            return 16;
			if(sname.Contains("QCD") || 
			   sname == "MuEnr10"    ||
			   sname == "MuEnr15"    ||
			   sname == "EMEnr20"    ||
			   sname == "EMEnr30")                          return 17;
			if(sname == "TbZ")                              return 18;
			if(sname == "TTJets_matchingdown")				return 19;
			if(sname == "TTJets_matchingup")				return 20;
			if(sname == "TTJets_scaledown")					return 21;
			if(sname == "TTJets_scaleup")					return 22;
			else {
				cout << "SSDLDumper::Sample::getProc() ==> ERROR: "<< sname << " has no defined process!" << endl;
				return -1;
			}
		}
		TString getProcName(int proc){ // used for binned samples
			if(proc ==  0) return "Data";
			if(proc ==  1) return "$t\\bar{t}$";
			if(proc ==  2) return "Single t";
			if(proc ==  3) return "W+jets";
			if(proc ==  4) return "Z+jets";
			if(proc ==  5) return "$\\gamma$+jets";
			if(proc ==  6) return "WW";
			if(proc ==  7) return "WZ";
			if(proc ==  8) return "ZZ";
			if(proc ==  9) return "V$\\gamma$+jets";
			if(proc == 10) return "$t\\bar{t}$H";
			if(proc == 11) return "$t\\bar{t}$W";
			if(proc == 12) return "$t\\bar{t}$Z";
			if(proc == 13) return "$t\\bar{t}\\gamma$";
			if(proc == 14) return "W$^{\\pm}$W$^{\\pm}$";
			if(proc == 15) return "Tri-Boson";
			if(proc == 16) return "DPS (2$\\times$ W+jets)";
			if(proc == 17) return "QCD";
			if(proc == 18) return "TbZ";
			if(proc == 19) return "$t\\bar{t}$matchingdown";
			if(proc == 20) return "$t\\bar{t}$matchingup";
			if(proc == 21) return "$t\\bar{t}$scaledown";
			if(proc == 22) return "$t\\bar{t}$scaleup";
			else {
				cout << "SSDLDumper::Sample::getProcName() ==> ERROR: "<< proc << " has no defined process name!" << endl;
				return "";
			}
		}
		inline int getNProcs(){return 17;} // make sure this number corresponds to the number of
		                                   // processes define in the previous method
		TTree* getTree(){
			if (!file || !file->IsOpen()) file = TFile::Open(location);
			if(file->IsZombie()){
				cout << "SSDLDumper::Sample::getTree() ==> Error opening file " << location << endl;
				exit(1);
			}
			tree = (TTree*)file->Get("Analysis");
			return tree;
		}
		TH1F* getEvCount(){
			if (!file || !file->IsOpen()) file = TFile::Open(location);
			if(file->IsZombie()){
				cout << "SSDLDumper::Sample::getEvCount() ==> Error opening file " << location << endl;
				exit(1);
			}
			return (TH1F*)file->Get("EventCount");
		}
		
		void cleanUp(){
			tree = NULL;
			file->Close();
		}
	};

	SSDLDumper(TString configfile = "dumperconfig.cfg");
	virtual ~SSDLDumper();

	virtual void init();
	virtual void init(TString); // samples read from a datacard
	virtual void init(TString, TString, int, float, int = -1); // running on single sample

	virtual void readDatacard(TString); // read in a datacard

	virtual void loop();              // loop on all samples if several
	virtual void loopEvents(Sample*); // perform loop on single sample
	
        /////////////////////////
        bool IsInJSON();
	// old void storeNumbers(Sample*, gChannel, gRegion);
	void storeNumbers(Sample*, gChannel, int);
	
	// Cutflow
	void initCutNames();
	void initCounters();
	void fillCutFlowHistos(Sample*);
	void printCutFlow(gChannel, gSample, gSample);
        void printCutFlow(gChannel);
	void printCutFlows(TString);
	
	//////////////////////////////
	// Fillers
	void fillSigEventTree(Sample*, int);
	void resetSigEventTree();
	// old void fillYields(Sample*, gRegion);
	void fillYields(Sample*, int);
	void fillDiffYields(Sample*);
	void fillDiffVar(Sample* S, int lep1, int lep2, float val, int bin, gChannel chan);
	void fillDiffVarOS(Sample* S, int lep1, int lep2, float val, int bin, gChannel chan);
	void fillRatioPlots(Sample*);
	void fillMuIsoPlots(Sample*);
	void fillElIsoPlots(Sample*);
	void fillElIdPlots (Sample*);
        void fillKinPlots(Sample*, int);
        void fillSyncCounters(Sample*); 
	
	//////////////////////////////
	// I/O
	void bookSigEvTree();
	void writeSigEvTree(TFile*);
	int getSampleType(Sample*);

	void bookHistos(Sample*);
	void deleteHistos(Sample*);
	void writeHistos(Sample*, TFile*);
	void writeSigGraphs(Sample*, gChannel, TFile*);
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

	virtual bool isGoodRun(Sample*);

	// Event and Object selectors:
	virtual void scaleBTags(Sample *S, int flag = 0);
	virtual void saveBTags();
	virtual void resetBTags();
	virtual void smearJetPts(Sample *S, int flag = 0);
	virtual void scaleLeptons(Sample *S, int flag = 0);
	virtual void smearMET(Sample *S);
	virtual void propagateMET(TLorentzVector vec1, TLorentzVector vec2);
	virtual void scaleMET(Sample *S, int flag = 0);
	virtual float getJetPt(int); // for shifting and smearing
	virtual float getM3();
	virtual float getMET();
	virtual float getMETPhi();
	virtual int getNJets();
	virtual int getNBTags();
	virtual int getNBTagsMed();
	virtual std::vector< int > getNBTagsMedIndices();
	virtual float getHT();
	virtual float getWeightedHT();
	virtual float getMT(int, gChannel);
	virtual float getMT2(int, int, gChannel);
	virtual float getMll(int, int, gChannel);
  //        virtual float getDPhi(int, int, gChannel);
	virtual int   getClosestJet(int, gChannel);
	virtual float getClosestJetPt(int, gChannel);
	virtual float getClosestJetDR(int, gChannel);
	virtual float getClosestJetDPhi(int, gChannel);
	virtual float getSecondClosestJetDR(int, gChannel);
	virtual float getAwayJetPt(int, gChannel);
	virtual float getMaxJPt();
	
	virtual int getNTightMuons();
	virtual int getNTightElectrons();
	
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
	virtual bool passesJet50Cut();
	
	virtual bool passesHTCut(float, float = 8000.);
	virtual bool passesMETCut(float = -1., float = 8000.);
	virtual bool passesZVeto(bool(SSDLDumper::*)(int), bool(SSDLDumper::*)(int), float = 15.); // cut with mZ +/- cut value and specified obj selectors
	virtual bool passesZVeto(float = 15.); // cut with mZ +/- cut value
	virtual bool passesZVetoNew(int l1, int l2, int toggle, float dm = 15.);
	virtual bool passesGammaStarVeto(int l1, int l2, int toggle, float mass = 12.);
	virtual bool passesChVeto(int = 1);
	virtual bool passesMllEventVeto(float = 5.);
	virtual bool passesMllEventVeto(int, int, int, float = 5.);
	virtual bool passes3rdLepVeto();
	virtual bool passesTauVeto();
	virtual bool passesTTZSel();
        //virtual bool passesGammaStarVeto();

	virtual int isSigSupMuEvent();
	virtual bool isZMuMuEvent(int&, int&);

	virtual int isSigSupElEvent();
	virtual bool isZElElEvent(int&, int&);
	virtual bool isZElElChMisIdEvent(int&, int&);

	// virtual bool isGenMatchedSUSYDiLepEvent();
	// virtual bool isGenMatchedSUSYDiLepEvent(int&, int&);

	virtual bool AvoidDoubleCountingOfFakes(Sample*);
	virtual bool isGenMatchedSUSYMuMuEvent();
	virtual bool isGenMatchedSUSYEEEvent();
	virtual bool isGenMatchedSUSYEMuEvent();

	virtual bool isSSLLMuEvent(int&, int&);
	virtual bool isSSLLElEvent(int&, int&);
	virtual bool isSSLLElMuEvent(int&, int&);

	// old virtual bool passesPtCuts(float pT1, float pT2, gRegion reg, gChannel chan);
	virtual bool passesPtCuts(float pT1, float pT2, int reg, gChannel chan);
	
	virtual bool isGoodMuon(int, float = -1.);
	virtual bool isGoodMuonForZVeto(int);
	virtual bool isGoodMuonForGammaStarVeto(int);  
	virtual bool isGoodMuonFor3rdLepVeto(int);
	virtual bool isGoodMuonForTTZ(int, float = 20.);
	virtual bool isLooseMuon(int);
	virtual bool isTightMuon(int);
	virtual bool isGoodPrimMuon(int, float = -1.);
	virtual bool isGoodSecMuon(int, float = -1.);

	virtual bool isFakeMuon(int);
	virtual bool isPromptMuon(int);
	virtual bool isChargeMatchedMuon(int);

	virtual bool isGoodElectron(int, float = -1.);
	virtual bool isGoodEleForZVeto(int);
	virtual bool isGoodEleForGammaStarVeto(int);  
	virtual bool isGoodEleFor3rdLepVeto(int);
	virtual bool isGoodEleForTTZ(int, float = 20.);
	virtual bool isGoodEleForChMId(int, float = 20.);
	virtual bool isLooseElectron(int);
	virtual bool isTightElectron(int);
	virtual bool isGoodPrimElectron(int, float = -1.);
	virtual bool isGoodSecElectron(int, float = -1.);

	virtual bool isFakeElectron(int);
	virtual bool isPromptElectron(int);
	virtual bool isChargeMatchedElectron(int);

	virtual bool isBarrelElectron(int);

	virtual bool isGoodTau(int);

	virtual bool isGoodJet(int, float = 20.);
	virtual float getJERScale(int);
	

	float getHLTSF_DoubleElectron( float pt, float eta, const std::string& runPeriod="" );
	float getHLTSF_MuEG(           float pt, float eta, const std::string& runPeriod="" );
	float diMuonHLTSF2012();
	float muEleHLTSF2012();
	float diEleHLTSF2012();

	ValueAndError getMuMuHLTSF( float pt, float eta, const std::string& runPeriod );
	ValueAndError getMuonRecoSF(              float pt, float eta );
	ValueAndError getMuonIsoSF(               float pt, float eta );
	ValueAndError getElectronRecoSF(          float pt, float eta );
	ValueAndError getElectronIsoSF(           float pt, float eta );

	float getMuScale(float pt, float eta);
	float getElScale(float pt, float eta);

	float fC_minHT;
	float fC_minMet;
	float fC_maxHT;
	float fC_maxMet;
	float fC_minJetPt;
	int   fC_minNjets;
	int   fC_maxNjets;
	int   fC_minNbjets;
	int   fC_maxNbjets;
	int   fC_minNbjmed;
	int   fC_maxNbjmed;
	float fC_minMu1pt;
	float fC_minMu2pt;
	float fC_minEl1pt;
	float fC_minEl2pt;
	float fC_maxMet_Control;
	float fC_maxMt_Control;
	int   fC_app3rdVet; // 3rd lepton veto
	int   fC_vetoTTZSel; // ttZ veto
        int   fC_chargeVeto;
        int   fC_GStarVeto;
	void resetHypLeptons();
	void setHypLepton1(int, gChannel);
	void setHypLepton2(int, gChannel);
	void setHypLepton3(int, gChannel);
	lepton fHypLepton1;
	lepton fHypLepton2;
	lepton fHypLepton3;
	
	// old void setRegionCuts(gRegion reg = Baseline);
	void setRegionCuts(int reg);
	
	const int     getNFPtBins(gChannel); // fake ratios
	const double *getFPtBins (gChannel);
	const int     getNPPtBins (gChannel); // prompt ratios
	const double *getPPtBins  (gChannel);
	const int     getNEtaBins(gChannel);
	const double *getEtaBins (gChannel);
        const int     getNCMidbins() { return gNElCMIdbins; };
        const double *getCMIdbins()  { return gElCMIdbins;  };
	const double *getDiffPredBins(int);
	
	Monitor fCounters[gNSAMPLES][3];
	vector<string> fMMCutNames;
	vector<string> fEMCutNames;
	vector<string> fEECutNames;
	bool fDoCounting;
	gSample fCurrentSample;
	gChannel fCurrentChannel;
	ofstream fOUTSTREAM, fOUTSTREAM2, fOUTSTREAM3, fOUTSTREAM4;

	int fChargeSwitch;    // 0 for SS, 1 for OS

	float fLumiNorm;      // Normalize everything to this luminosity

	vector<Sample*>::iterator fS;
	vector<Sample*> fSamples;
	vector<Sample*> fMCSamples;
	map<TString, Sample*> fSampleMap;	// Mapping of sample to name
	// map<TString, gSample> fSampleMap;	// Mapping of sample number to name
	
	TTree *fSigEv_Tree;
	int         fSETree_SystFlag; // 0 nominal, 1 jets up, 2 jets dn, 3 jets smeared, 4 btag up, 5 btag dn, 6 lep up, 7 lep dn
	float       fSETree_PUWeight;
	float       fSETree_HLTSF;
	float       fSETree_BtagSF1;
	float       fSETree_BtagSF2;
	float       fSETree_SLumi;
	std::string fSETree_SName;
	int         fSETree_SType; // DoubleMu(0), DoubleEle(1), MuEG(2), MC(10)
	int         fSETree_Run;
	int         fSETree_LS;
	long        fSETree_Event;
	int         fSETree_Flavor; // mm(0), em(1), ee(2)
	int         fSETree_Charge;
	int         fSETree_TLCat; // TL category: TT(0), TL(1), LT(2), LL(3)
	int         fSETree_ZVeto;   // passes Z veto
	int         fSETree_3rdVeto; // passes 3rd lepton veto
	int         fSETree_ttZSel;  // passes ttz sel
	float       fSETree_HT;
	float       fSETree_MET;
	int         fSETree_NM; // number of tight muons
	int         fSETree_NE; // number of tight electrons
	int         fSETree_NJ;
	int         fSETree_NbJ;
	int         fSETree_NbJmed;
	float       fSETree_M3;
	float       fSETree_MT2;
	float       fSETree_Mll;
	float       fSETree_pT1;
	float       fSETree_pT2;
	float       fSETree_eta1;
	float       fSETree_eta2;
	float		fSETree_PFIso1;
	float		fSETree_PFIso2;
	float		fSETree_MVAID1;
	float		fSETree_MVAID2;
	float		fSETree_medWP1;
	float		fSETree_medWP2;

	vector<float> fSigEv_HI_MM_HT;
	vector<float> fSigEv_HI_MM_MET;
	vector<float> fSigEv_HI_EE_HT;
	vector<float> fSigEv_HI_EE_MET;
	vector<float> fSigEv_HI_EM_HT;
	vector<float> fSigEv_HI_EM_MET;
	
	vector<float> fSaved_Tags;
	
	TFile *fStorageFile;
	TString fOutputFileName;
	Sample *fSample;
	
	// BTagSFUtil *fBTagSFUtil;
	BTagSF *fBTagSF;
        GoodRunList *fGoodRunList;
	TRandom3 *fRand3;

        Int_t fCurRun;
        Int_t fCurLumi;
        bool skipLumi;
        bool skipRun;
	bool isTChiSlepSnu;
	private:
	
	Monitor fCounter[3];	
        Monitor fCounterSync[3];
	vector<string> fSyncCutNames;
};

#endif
