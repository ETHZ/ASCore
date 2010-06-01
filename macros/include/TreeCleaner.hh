#ifndef TreeCleaner_hh
#define TreeCleaner_hh


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

#include "base/TreeReader.hh"
#include "base/UserAnalysisBase.hh"
#include "helper/Davismt2.h"

class TreeCleaner : public UserAnalysisBase{
public:
	TreeCleaner(TreeReader *tr = 0);
	virtual ~TreeCleaner();

	void Begin();
	void Analyze();
	void End();
	void DoTagging();
	void DoCleaning();
	void DoSkimTree();
	void Reset();
	void StatWrite(TString fCleanerStats);

	bool fClean;
	bool fSkim;
	inline void SetClean(bool clean){fClean = clean;};
	inline void SetSkim(bool skim){fSkim = skim;};

	TTree *fCleanTree;
	TFile *fCleanTreeFile;
	
	TFile *fHstatFile;
	TH1D *fHstatHistos;
	
	double fR12;
	double fR21;

	// Cleaning parameters:
	// Number in comments are initial values
        // -- Cut values
        float fMinJetPt;                    // = 20.   // Min Pt of jets
	// -- Primary vertex:
	float fClean_chisqVxmax;            // = 5.0   // Max nchi2, nchi2 is also cut at > 0
	float fClean_dRVxmax;               // = 0.25  // Max transverse distance to beamspot
	float fClean_dzVxmax;               // = 20.0  // Max longitudinal distance to beamspot
	float fClean_sumPtTkfromVxmin;      // = 0.0   // Min summed pt of tracks assoc. with vtx
	int fClean_PrimVtxNdofmin;          // = 5     // Min Ndof for Primary Vertex fit

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
	float fClean_ElecSigmaEtaEtaBarmin; // = 0.002 // Min ElecSigmaEtaEta (barrel), against spikes
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
        float fClean_ElecConvPartTrackDistmax; // = 0.    // Max distance for conversion
        float fClean_ElecConvPartTrackDCotmax; // = 0.    // Max dCosTheta for conversion
        float fClean_ElecNMissHitsmax;      // = 1.    // Max number of missing hits
	float fClean_dRSSelecmax;           // = 10.   // Max delta R of same sign elec. duplicate check

	// -- Photons:
	float fClean_PhotHoverEBarmax;      // = 0.2   // Max PhotonHcalOverEcal (barrel)
	float fClean_PhotHoverEEndmax;      // = 0.2   // Max PhotonHcalOverEcal (endcap)
	float fClean_PhotSigmaEtaEtaBarmin; // = 0.002 // Min PhotonSigmaEtaEta (barrel), against spikes
	float fClean_PhotSigmaEtaEtaBarmax; // = 0.011 // Max PhotonSigmaIetaIeta (barrel)
	float fClean_PhotSigmaEtaEtaEndmax; // = 0.025 // Max PhotonSigmaIetaIeta (endcap)
	
	// -- Jets:
	float fClean_FracEmmaxJet;          // = 1.0   // Max JEMfrac
	float fClean_FracEmminJet;          // = 0.01  // Min JEMfrac
	float fClean_FracChminJet;          // = 0.05  // Min JChfrac
	float fClean_JID_n90Hitsmin;        // = 2     // Min JID_n90Hits
	float fClean_JID_HPDmax;            // = 0.98  // Max JID_HPD
	float fClean_JID_RBXmax;            // = 0.95  // Max JID_RBX

	float fClean_deltaRElecJetmax;      // = 0.5   // Max delta R of e and j for electron jet check
	float fClean_elecbyJetEratio;       // = 0.7   // Min eE/jE to be considered electron jet
	
	// -- Isolation:
	float fClean_MuonIsomax;            // = 1.    // Max relative iso cut (muons)
	float fClean_ElecIsomax;            // = 1.    // Max relative iso cut (electrons)
	float fClean_PhotIsomax;            // = 1.    // Max relative iso cut (photons)
	
	// -- Event cleaning:
	float fClean_FracChmin;             // = 0.1   // Min charge fraction in event
	float fClean_FracEmmin;             // = 0.175 // Min EM fraction in event
	float fClean_JetBadHardPtmin;       // = 50.   // Min Pt for jet to trigger bad event
	
	// -- MET:
	float fClean_METmin;                // = 50.0  // Min MET to be considered
	float fClean_dPhiJetMETmin;         // = 0.0   // Min phi distance of MET to closest jet
	float fClean_dR12min;               // = 0.5   // Min R12 = sqrt(dPhi1^2 + (PI-dPhi2)^2)
	float fClean_dR21min;               // = 0.5   // Min R21 = sqrt(dPhi2^2 + (PI-dPhi1)^2)

private:

	// Functions performing the cleaning and duplicate tagging 
	void ReadCleaningParameters(const char* filename = "cleaningparms.dat");
	void TagCleanObjects(void);
	void TagDuplObjects(void);
	int CleanPrimaryVertex(void);
	int IsFromPrimaryVx(int ipart, int ichk);
	int CleanMuon(int ichk);
	bool DuplicateMuon(int ichk);
	int CleanElectron(int ichk);
        bool DuplPhotonElectron(int ichk);
	bool DuplicateElectron(int ichk);
	int CleanPhoton(int ichk);
        int CleanJet(int ichk, int itype);
	bool ElectronJet(int ichk);
	bool PhotonJet(int ichk);
	int FindNearestJet(double eta, double phi);
	double DeltaPhi(double v1, double v2);
	double GetDeltaR(double eta1, double eta2, double phi1, double phi2);
	
	// Functions to actually perform the cleaning
	void DecideIso(void);
	void InitCleaning();
	void DoCleanObjects(void);
	void AddToJet(int ipart, int ichk, int iJet);
	void SubtrFromJet(int ipart, int ichk, int iJet);
	int CleanEvent(void);
	int CleanMET(double met, double metphi);

	// functions to make cleaning statistics
	void StatInit(const char *filename = "CleaningStats.root");
	void StatFill(void);
	void StatPrint(void);
	void StatHistos(void);

	void PutMuon(int, int);
	void PutElectron(int, int);
	void PutPhoton(int, int);
	void PutJet(int, int);

	// Cleaning Statistics Counters
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
	int fNumTotMuonBadGlTr;  

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
	int fNumTotElecBadSpik;
	int fNumTotElecBadHits;
	int fNumTotElecBadConv;

	int fNumTotPhotons;
	int fNumTotPhotGoodIso;
	int fNumTotPhotGoodNonIso;
	int fNumTotPhotBadIso;
	int fNumTotPhotBadNonIso;
	int fNumTotPhotDupl;
	int fNumTotPhotNotClean;
	int fNumTotPhotBadHoE;
	int fNumTotPhotBadShsh;
	int fNumTotPhotBadSpik;

	int fNumTotJets;  
	int fNumTotJetGood;  
	int fNumTotJetBad;  
	int fNumTotJetDuplElJet;
	int fNumTotJetDuplPhoJet;
	int fNumTotJetNotPrimaryTrk;
	int fNumTotJetNotClean;
	int fNumTotJetPgtE;
	int fNumTotJetGtFem;
	int fNumTotJetLtFem;
	int fNumTotJetLtFch;
	int fNumTotJetLtn90hits;
	int fNumTotJetGtfHPD;
	int fNumTotJetGtfRBX;
	int fNumTotBJets;  

	// Cleaning variables:
	int fNMuClean;
	int fNElClean;
	int fNPhClean;
	int fNJClean;
	
};
#endif
