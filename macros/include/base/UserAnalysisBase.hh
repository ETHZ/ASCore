#ifndef UserAnalysisBase_hh
#define UserAnalysisBase_hh

#include "TreeReader.hh"
#include "helper/pdgparticle.hh"
#include "helper/Utilities.hh"
#include "helper/PUWeight.h"
#include <map>
#include <string>

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

	inline void SetOutputDir(TString dir){ fOutputDir = Util::MakeOutputDir(dir); };
	inline void SetOutputFile(TString file){ fOutputFile = Util::MakeOutputFile(fOutputDir + file); };

	virtual void ReadPDGTable(const char* filename = "pdgtable.txt");
	virtual int GetPDGParticle(pdgparticle&, int);
	virtual void GetHLTNames(Int_t& run);
	virtual int GetHLTBit(string);
	virtual bool GetHLTResult(string);
	virtual int GetHLTPrescale(string);

	// PileUp reweighting;
	virtual void  SetPileUpSrc(string, string);
	virtual float GetPUWeight(int);

	TreeReader *fTR;
	TString fOutputDir;
	TFile *fOutputFile;
	TString fTag;
	TLatex *fTlat;

	
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
	virtual bool IsGoodBasicEl		(int);
	virtual bool IsElInGap(int);
	virtual bool IsElFromPrimaryVx	(int);

	virtual bool IsGoodElId_WP80	(int);
	virtual bool IsGoodElId_WP90	(int);

	virtual bool IsTightEl			(int);
	virtual bool IsLooseEl			(int);
	virtual bool IsLooseNoTightEl	(int);
	
	virtual double hybRelElIso		(int);
	virtual bool IsIsolatedEl		(int, double, double);

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
	PUWeight  *fPUWeight;

};

#endif
