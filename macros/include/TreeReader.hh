#ifndef TreeReader_hh
#define TreeReader_hh


#include <vector>
#include <cmath>
#include <string>
#include <iostream>
#include <fstream>

#include <TH2D.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TStyle.h>
#include <TLatex.h>
#include <TLeaf.h>
#include <TBranch.h>

#include "TreeClassBase.h"
#include "LeptJetStat.h"
#include "Davismt2.h"
// #include "ETHStyle.h"

class TreeReader : public TreeClassBase{
public:
	TreeReader(TTree *tree=0, int flag = 111);
	virtual ~TreeReader();
	void DefStyle();
	void BeginJob();
	void EndJob();
	void Loop();

	void SetOutputDir(TString dir);
	inline void SetTag(TString tag){fTag = tag;};

	void BookSignHists(const char* filename = "SignificancePlots.root");
	void FillSignHists(Int_t part);
	void WriteSignHists();

	void BookMPHistos(const char* filename = "MultiplicityPlots.root");
	void FillMPHistos();
	void PrintMPOutput();
	void PlotMPSummary();
	void PlotMPEffic();

	void InitDiLepTree(const char *filename = "DiLepTree.root");
	void FillDiLepTree();
	void ResetDiLepTree();
	void WriteDiLepTree();

	void PlotMultiplicity();

	void printPNG(TCanvas*, TString, TString);
	void printEPS(TCanvas*, TString, TString);
	
	double getEta(double, double, double);
	
	void ReadObjCuts(const char* = "objsel.dat");
	void ReadEvtSel(const char* = "evtsel.dat");
	
	void WriteTree(TString);
	
private:
	
	// Functions performing the cleaning and duplicate tagging 
	virtual void ReadCleaningParameters(const char* filename = "cleaningparms.dat");
	virtual void TagCleanObjects(void);
	virtual int CleanPrimaryVertex(void);
	virtual int IsFromPrimaryVx(int ipart, int ichk);
	virtual int CleanMuon(int ichk);
	virtual bool DuplicateMuon(int ichk);
	virtual int CleanElectron(int ichk);
	virtual bool DuplicateElectron(int ichk);
	virtual int CleanPhoton(int ichk);
	virtual int CleanJet(int ichk);
	virtual bool ElectronJet(int ichk);
	virtual int FindNearestJet(double eta, double phi);
	virtual double DeltaPhi(double v1, double v2);
	virtual double GetDeltaR(double eta1, double eta2, double phi1, double phi2);
	
	// Functions to actually perform the cleaning
	virtual void DecideIso(void);
	virtual void InitCleaning();
	virtual void DoCleanObjects(void);
	virtual void AddToJet(int ipart, int ichk, int iJet);
	virtual void SubtrFromJet(int ipart, int ichk, int iJet);
	virtual int CleanEvent(void);
	virtual int CleanMET(double met, double metphi);

	virtual void PutMuon(int, int);
	virtual void PutElectron(int, int);
	virtual void PutPhoton(int, int);
	virtual void PutJet(int, int);


	// functions to make cleaning statistics
	virtual void StatInit(const char *filename = "CleaningStats.root");
	virtual void StatFill(void);
	virtual void StatPrint(void);
	virtual void StatHistos(void);

// variables for cleaning statistics
	int fNumTotEvt;
	int fNumTotEvtReject;
	int fNumTotEvtEmpty;
	int fNumTotEvtCleanEmpty;
	int fNumTotEvtLtFem;
	int fNumTotEvtLtFch;
	int fNumTotEvtPfMETJet;
	int fNumTotEvtPfMETRij;
	int fNumTotEvtCaMETJet;
	int fNumTotEvtCaMETRij;
	int fNumTotEvtTcMETJet;
	int fNumTotEvtTcMETRij;
	int fNumTotEvtBadHardJet;

	int fNumTotMuons;  
	int fNumTotMuonGoodIso;  
	int fNumTotMuonGoodNonIso;  
	int fNumTotMuonBadIso;  
	int fNumTotMuonBadNonIso;  
	int fNumTotMuonDupl;
	int fNumTotMuonNotPrimaryTrk;
	int fNumTotMuonNotClean;
	int fNumTotMuonBadDpop;
	int fNumTotMuonBadChi2;  
	int fNumTotMuonBadNhit;  

	int fNumTotElectrons;
	int fNumTotElecGoodIso;
	int fNumTotElecGoodNonIso;
	int fNumTotElecBadIso;
	int fNumTotElecBadNonIso;
	int fNumTotElecDupl;
	int fNumTotElecNotPrimaryTrk;
	int fNumTotElecNotClean;
	int fNumTotElecBadHoE;
	int fNumTotElecBadShsh;
	int fNumTotElecBadTmat;

	int fNumTotPhotons;
	int fNumTotPhotGoodIso;
	int fNumTotPhotGoodNonIso;
	int fNumTotPhotBadIso;
	int fNumTotPhotBadNonIso;
	int fNumTotPhotDupl;
	int fNumTotPhotNotClean;
	int fNumTotPhotBadHoE;
	int fNumTotPhotBadShsh;

	int fNumTotJets;  
	int fNumTotJetGood;  
	int fNumTotJetBad;  
	int fNumTotJetDuplElJet;
	int fNumTotJetNotPrimaryTrk;
	int fNumTotJetNotClean;
	int fNumTotJetPgtE;
	int fNumTotJetGtFem;
	int fNumTotJetLtFem;
	int fNumTotJetLtFch;
	int fNumTotBJets;  

	// Global parameters:
	TString fOutputDir;
	TString fTag;
	bool fClean;
	bool fDiLep;
	bool fMPHist;
	bool fSignHist;
	
	TStyle *fStyle;
	Davismt2 *fMT2;
	TLatex *fTlat;

	// Cleaning variables:
	int fNMuClean;
	int fNElClean;
	int fNPhClean;
	int fNJClean;

	// Cleaning parameters:
	// Number in comments are initial values
	// -- Primary vertex:
	float fClean_chisqVxmax;            // = 5.0   // Max nchi2, nchi2 is also cut at > 0
	float fClean_dRVxmax;               // = 0.25  // Max transverse distance to beamspot
	float fClean_dzVxmax;               // = 20.0  // Max longitudinal distance to beamspot
	float fClean_sumPtTkfromVxmin;      // = 0.0   // Min summed pt of tracks assoc. with vtx

	// -- Vertex compatibility (of muons, electrons, photons and jets)
	float fClean_distVxmax;             // = 5.0   // Max deviation (in sigmas) for d0 and dz
	
	// -- Muons:
	float fClean_MuonDPbyPmax;          // = 0.5   // Max PtError/Pt
	float fClean_MuonChi2max;           // = 10.0  // Max nchi2
	float fClean_MuonNHitsmin;          // = 11.0  // Min number of tracker hits
	float fClean_dRSSmuonmax;           // = 0.1   // Max delta R of same sign muon duplicate check
	
	// -- Electrons:
	float fClean_ElecHoverEBarmax;      // = 0.045 // Max ElHcalOverEcal (barrel)
	float fClean_ElecHoverEEndmax;      // = 0.05  // Max ElHcalOverEcal (endcap)
	float fClean_ElecSigmaEtaEtaBarmax; // = 0.011 // Max ElSigmaIetaIeta (barrel)
	float fClean_ElecSigmaEtaEtaEndmax; // = 0.025 // Max ElSigmaIetaIeta (endcap)
	float fClean_ElecEoverPInBarmin;    // = 0.3   // Min ElESuperClusterOverP (barrel)
	float fClean_ElecEoverPInEndmin;    // = 0.4   // Min ElESuperClusterOverP (endcap)
	float fClean_ElecDeltaEtaInBarmax;  // = 0.007 // Max ElDeltaEtaSuperClusterAtVtx (barrel)
	float fClean_ElecDeltaEtaInEndmax;  // = 0.007 // Max ElDeltaEtaSuperClusterAtVtx (endcap)
	float fClean_ElecDeltaPhiInBarmax;  // = 0.06  // Max ElDeltaPhiSuperClusterAtVtx (barrel)
	float fClean_ElecDeltaPhiInEndmax;  // = 0.06  // Max ElDeltaPhiSuperClusterAtVtx (endcap)
	float fClean_ElecDeltaPhiOutBarmax; // = 999.0 // Max ElDeltaPhiSeedClusterAtCalo (barrel)
	float fClean_ElecDeltaPhiOutEndmax; // = 999.0 // Max ElDeltaPhiSeedClusterAtCalo (endcap)
	float fClean_dRSSelecmax;           // = 10.   // Max delta R of same sign elec. duplicate check

	// -- Photons:
	float fClean_PhotHoverEBarmax;      // = 0.2   // Max PhotonHcalOverEcal (barrel)
	float fClean_PhotHoverEEndmax;      // = 0.2   // Max PhotonHcalOverEcal (endcap)
	float fClean_PhotSigmaEtaEtaBarmax; // = 0.011 // Max PhotonSigmaIetaIeta (barrel)
	float fClean_PhotSigmaEtaEtaEndmax; // = 0.025 // Max PhotonSigmaIetaIeta (endcap)
	
	// -- Jets:
	float fClean_FracEmmaxJet;          // = 1.0   // Max JEMfrac
	float fClean_FracEmminJet;          // = 0.01  // Min JEMfrac
	float fClean_FracChminJet;          // = 0.05  // Min JChfrac

	float fClean_deltaRElecJetmax;      // = 0.5   // Max delta R of e and j for electron jet check
	float fClean_elecbyJetEratio;       // = 0.7   // Min eE/jE to be considered electron jet
	
	// -- Isolation:
	float fClean_MuonIsomax;            // = 1.    // Max relative iso cut (muons)
	float fClean_ElecIsomax;            // = 1.    // Max relative iso cut (electrons)
	float fClean_PhotIsomax;            // = 1.    // Max relative iso cut (photons)
	
	// -- Event cleaning:
	float fClean_FracChmin;             // = 0.1   // Min charge fraction in event
	float fClean_FracEmmin;             // = 0.175 // Min EM fraction in event
	
	// -- MET:
	float fClean_METmin;                // = 50.0  // Min MET to be considered
	float fClean_dPhiJetMETmin;         // = 0.0   // Min phi distance of MET to closest jet
	float fClean_dR12min;               // = 0.5   // Min R12 = sqrt(dPhi1^2 + (PI-dPhi2)^2)
	float fClean_dR21min;               // = 0.5   // Min R21 = sqrt(dPhi2^2 + (PI-dPhi1)^2)
	
	TTree *fCleanTree;
	TFile *fCleanTreeFile;
	
	// Significance Plots:
	TFile *fSignHistsFile;
	int fNBinsEta[5];
	int fNBinsPhi;
	TH2D *fH_ptdev[5];
	TH2D *fH_ptsum[5];
	TH2D *fH_pt2sum[5];
	TH2I *fH_ptevt[5];
	TH2D *fH_ptavg[5];
	TH1D *fH_ptsumeta[5];
	TH1I *fH_ptevteta[5];

	// cleaning statistica histogram
	TFile *fHstatFile;
	TH1D *fHstatHistos;
	
	// Multiplicity Plots Variables:
	TFile *fMPHistFile;
	LeptJetStat *fMyLeptJetStat;
	TH2D *fHljMult;
	TH2D *fHemuMult;
	TH1F *fHemuEff;

	// DiLepton Tree Variables:
	TFile *fDiLepTreeFile;
	TTree *fDiLepTree;
	int    fTRunNumber;
	int    fTEventNumber;
	int    fTLumiSection;

	int fTMu1charge;
	int fTMu2charge;
	int fTNqualmu;
	double fTMu1pt;
	double fTMu2pt;
	double fTMu1eta;
	double fTMu2eta;
	double fTMu1iso;
	double fTMu2iso;
	double fTMu1d0;
	double fTMu2d0;
	double fTMu1ntkhits;
	double fTMu2ntkhits;
	double fTMuminv;
	double fTMumt2_50;
	double fTMumt2_100;

	int fTNqualel;
	int fTEl1charge;
	int fTEl2charge;
	double fTEl1pt;
	double fTEl2pt;
	double fTEl1eta;
	double fTEl2eta;
	double fTEl1iso;
	double fTEl2iso;
	double fTEl1d0;
	double fTEl2d0;
	double fTElminv;
	double fTElmt2_50;
	double fTElmt2_100;
	
	// Object quality cuts:
	struct Cut{
		TBranch *branch;
		double upperbound;
		double lowerbound;
	};
	bool IsGoodObj(int, std::vector<Cut>*);
	bool IsGoodEvt(std::vector<Cut>*);
	std::vector<Cut> fMuCuts;
	std::vector<Cut> fElCuts;
	std::vector<Cut> fJetCuts;
	std::vector<Cut> fEvtSelCuts;

};
#endif
