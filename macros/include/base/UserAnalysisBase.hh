#ifndef UserAnalysisBase_hh
#define UserAnalysisBase_hh

#include "TreeReader.hh"
#include "helper/pdgparticle.hh"
#include "helper/Utilities.hh"
#include "helper/PUWeight.h"
#include "helper/Lumi3DReWeighting_standalone.hh"
#include <map>
#include <string>

#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

class UserAnalysisBase{
public:
	UserAnalysisBase(TreeReader *tr = 0);
	virtual ~UserAnalysisBase();
  
	virtual void Begin() {}
    	virtual void BeginRun(Int_t& run);
    	virtual void Analyze() {}
    	virtual void End() {}
	inline virtual void SetTag(TString tag){fTag = tag;};
	inline virtual void SetVerbose(int verbose){fVerbose = verbose;};
	inline virtual void SetData(bool isdata){fIsData = isdata;};

	inline void SetOutputDir(TString dir){ fOutputDir = Util::MakeOutputDir(dir); };
	inline void SetOutputFile(TString file){ fOutputFile = Util::MakeOutputFile(fOutputDir + file); };

	virtual void ReadPDGTable(const char* filename = "pdgtable.txt");
	virtual int GetPDGParticle(pdgparticle&, int);
	virtual void GetHLTNames(Int_t& run);
	virtual int GetHLTBit(string);
	virtual bool GetHLTResult(string);
	virtual int GetHLTPrescale(string);

	// PileUp reweighting;
	virtual void  SetPileUpSrc(string, string = "");
  	virtual void  SetPileUp3DSrc(string, string);

	virtual float GetPUWeight(int);
	virtual float GetPUWeight(int, int);
  virtual float GetPUWeight3D( int , int , int );


	TreeReader *fTR;
	TString fOutputDir;
	TFile *fOutputFile;
	TString fTag;
	TLatex *fTlat;

	bool fIsData;

	int fVerbose;
	map<int, pdgparticle> fPDGMap; // Mapping of PDG ID names
	map<string, int> fHLTLabelMap; // Mapping of HLT trigger bit names
	vector<string>   fHLTLabels;   // Vector with current HLT names
	

	// Jet Selectors
	virtual bool IsGoodBasicJet(int);
	virtual bool IsGoodBasicPFJet   ( int, double = 30., double = 2.5);
	virtual bool IsGoodPFJetMedium  ( int, double = 30., double = 2.5);
	virtual bool IsGoodPFJetTight   ( int, double = 30., double = 2.5);
	virtual bool IsGoodBasicPFJetPAT( int, double = 30., double = 2.5);
	virtual bool IsGoodPFJetMediumPAT(int, double = 30., double = 2.5);
	virtual bool IsGoodPFJetTightPAT( int, double = 30., double = 2.5);

	// Muon Selectors
	virtual bool IsGoodBasicMu(int);
	virtual bool IsTightMu(int);
	virtual bool IsLooseMu(int);
	virtual bool IsLooseNoTightMu(int);

	// Electron Selectors
	virtual bool IsGoodBasicEl(int);
	virtual bool IsGoodElId_WP80(int);
	virtual bool IsGoodElId_WP90(int);
	virtual bool ElPassesWP80_ConvRej(int);
	virtual bool IsLooseEl(int);
	virtual float relElIso(int);

	// Photon Selectors
	virtual bool IsGoodBasicPho(int);

	// Event Selectors
	virtual bool IsGoodEvent();
	virtual bool IsGoodMuEvent();
	virtual bool IsGoodElEvent();
	virtual bool IsGoodElFakesEvent();
	virtual bool IsGoodHadronicEvent();
	virtual vector<int> MuonSelection(    bool(UserAnalysisBase::*muonSelector)(int) = NULL);
	virtual vector<int> ElectronSelection(bool(UserAnalysisBase::*eleSelector)(int) = NULL);
	virtual vector<int> PhotonSelection(  bool(UserAnalysisBase::*phoSelector)(int) = NULL);
	virtual vector<int> JetSelection(     bool(UserAnalysisBase::*jetSelector)(int) = NULL);
	virtual vector<int> PFJetSelection(double ptcut = 30., double absetacut = 2.5,   bool(UserAnalysisBase::*pfjetSelector)(int, double, double) = NULL);
	virtual bool SingleMuonSelection(int&);
	virtual bool DiMuonSelection(int&, int&, int = 0);
	virtual bool SSDiMuonSelection(int &, int&);	
	virtual bool SingleElectronSelection(int&, bool(UserAnalysisBase::*eleSelector)(int) = NULL);
	virtual bool DiElectronSelection(int&, int&, int = 0, bool(UserAnalysisBase::*eleSelector)(int) = NULL);
	virtual bool SSDiElectronSelection(int&, int&, bool(UserAnalysisBase::*eleSelector)(int) = NULL);
	virtual bool OSDiElectronSelection(int&, int&, bool(UserAnalysisBase::*eleSelector)(int) = NULL);

	// TDL & RA5
	virtual bool IsGoodElEvent_TDL();
	virtual bool IsGoodElEvent_RA5();

	// Print interesting event
	virtual void EventPrint();
	virtual void GetEvtEmChFrac(double & fracEm, double & fracCh);

        // Put all JES-related stuff between precompiler flags
#ifdef DOJES 
        FactorizedJetCorrector *fJetCorrector;
        JetCorrectionUncertainty *jecUnc;
	virtual float GetJetPtNoResidual(int);
#endif

private:

	// Object quality cuts:
	struct Cut{
		TBranch *branch;
		double upperbound;
		double lowerbound;
	};
	void ReadObjCuts(const char* = "objsel.dat");
	void ReadEvtSel(const char* = "evtsel.dat");
	bool IsGoodObj(int, std::vector<Cut>*);
	bool IsGoodEvt(std::vector<Cut>*);
	std::vector<Cut> fMuCuts;
	std::vector<Cut> fElCuts;
	std::vector<Cut> fJetCuts;
	std::vector<Cut> fEvtSelCuts;
	
	// Pile UP reweighting
  bool fDoPileUpReweight;
  bool fDoPileUpReweight3D;

	PUWeight  *fPUWeight;
  Lumi3DReWeighting   *fPUWeight3D;


};

#endif
