// -*- C++ -*-
//
// Package:    NTupleProducer
// Class:      NTupleProducer
//
/* class NTupleProducer
NTupleProducer.h
AnalysisExamples/NTupleProducer/src/NTupleProducer.h
Description: Template to produce NTuples for ETH SUSY Analysis

Implementation:

*/
//
// Original Author:  Benjamin Stieger
//         Created:  Wed Sep  2 16:43:05 CET 2009
// $Id: NTupleProducer.h,v 1.22 2009/11/27 15:48:38 stiegerb Exp $
//
//


// system include files
#include <vector>

// ROOT includes
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"

// Framework include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "PhysicsTools/UtilAlgos/interface/TFileService.h"

// Data formats
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/JetReco/interface/Jet.h"

// Helpers
#include "RecoJets/JetAlgorithms/interface/JetIDHelper.h"



class NTupleProducer : public edm::EDAnalyzer {
public:
	explicit NTupleProducer(const edm::ParameterSet&);
	~NTupleProducer();
	vector<double> calcMuIso(const reco::Muon *Mu, const edm::Event& iEvent);
	vector<double> calcMuIso2(const reco::Muon *Mu, const edm::Event& iEvent);
	vector<double> calcMuIso3(const reco::Muon *Mu, const edm::Event& iEvent, edm::Ref<reco::MuonCollection> muonRef);
	vector<double> calcElIso(const reco::GsfElectron *El, const edm::Event& iEvent);
	vector<int> matchMuCand(const reco::Muon *Mu, const edm::Event& iEvent);
	double DeltaPhi(double, double);
	double GetDeltaR(double, double, double, double);
	vector<const reco::Muon*> sortMus(vector<const reco::Muon*>);
	void switchDouble(double &, double &);
	void switchInt(int &, int &);
	void resetDouble(double *v, unsigned int size = 20);
	void resetInt(int *v, unsigned int size = 20);
	void resetTree();
	vector<const reco::Track*> FindAssociatedTracks(const reco::Jet *jet, const reco::TrackCollection *tracks);
private:
	virtual void beginJob(const edm::EventSetup&);
	virtual void beginRun(const edm::Run&, const edm::EventSetup&);
	virtual void analyze(const edm::Event&, const edm::EventSetup&);
	virtual void endJob();
	virtual void endRun(const edm::Run&, const edm::EventSetup&);

	typedef pair<int,double> OrderPair;
	struct IndexByPt {
		const bool operator()(const OrderPair& j1, 
		const OrderPair& j2 ) const {
			return j1.second > j2.second;
		}
	};

// ----------member data ---------------------------
	edm::Service<TFileService> fTFileService;
	reco::helper::JetIDHelper jetIDHelper;

	static const int fMaxnmus = 20;
	static const int fMaxneles = 20;
	static const int fMaxnjets = 100;
	static const int fMaxntrks = 500;

	bool fIsRealData;
	bool fIsPat;
	int fNTotEvents;
	int fNFillTree;

	edm::InputTag fMuonTag;
	edm::InputTag fElectronTag;
	edm::InputTag fEleIsoTkTag;
	edm::InputTag fEleIsoECTag;
	edm::InputTag fEleIsoHCTag;
	edm::InputTag fEleIsoDepTkTag;
	edm::InputTag fEleIsoDepECTag;
	edm::InputTag fEleIsoDepHCTag;
	edm::InputTag fMuIsoDepTkTag;
	edm::InputTag fMuIsoDepECTag;
	edm::InputTag fMuIsoDepHCTag;
	edm::InputTag fSCTag;
	edm::InputTag fJetTag;
	edm::InputTag fBtagTag;
	edm::InputTag fMET1Tag;
	edm::InputTag fMET2Tag;
	edm::InputTag fMET3Tag;
	edm::InputTag fMET4Tag;
	edm::InputTag fMET5Tag;
	edm::InputTag fVertexTag;
	edm::InputTag fTrackTag;
	edm::InputTag fCalTowTag;
	edm::InputTag fGenPartTag;
	edm::InputTag fTriggerTag;

	double fMinmupt;
	double fMaxmueta;
	double fMinelpt;
	double fMaxeleta;
	double fMaxeliso;
	double fMaxeld0;
	double fMinjpt;
	double fMaxjeta;
	double fMinjemfrac;
	double fMintrkpt;
	double fMaxtrketa;
	double fMaxtrknchi2;
	int fMintrknhits;
	double fIso_MuTkDRin;
	double fIso_MuTkDRout;
	double fIso_MuTkSeed;
	double fIso_MuCalDRin;
	double fIso_MuCalDRout;
	double fIso_MuCalSeed;

	double fJUNC_px_match[fMaxnjets];
	double fJUNC_py_match[fMaxnjets];
	double fJUNC_pz_match[fMaxnjets];

	TH1I *fHtrigstat; // Added to keep track of trigger names
	bool fFirstevent;

////////////////////////////////////////////////////////
// Trees:
	TTree *fRunTree;

	int fRTrunnumber;
	double fRTextxslo;
	double fRTextxsnlo;
	double fRTintxs;

	double fRTMinmupt;
	double fRTMaxmueta;
	double fRTMinelpt;
	double fRTMaxeleta;
	double fRTMaxeliso;
	double fRTMaxeld0;
	double fRTMinjpt;
	double fRTMaxjeta;
	double fRTMinjemfrac;

	double fRTMintrkpt;
	double fRTMaxtrketa;
	double fRTMaxtrknchi2;
	int fRTMintrknhits;

	double fRTIsoMuTkDRin;
	double fRTIsoMuTkDRout;
	double fRTIsoMuTkSeed;
	double fRTIsoMuCalDRin;
	double fRTIsoMuCalDRout;
	double fRTIsoMuCalSeed;

	TTree *fEventTree;

// General event information
	int fTrunnumber;
	int fTeventnumber;
	int fTlumisection;
	int fTsigprocid;
	double fTextxslo;
	double fTintxs;
	double fTweight;

	double fTprimvtxx;
	double fTprimvtxy;
	double fTprimvtxz;
	double fTprimvtxxE;
	double fTprimvtxyE;
	double fTprimvtxzE;
	double fTpvtxznchi2;
	int fTpvtxntracks;

	double fTbeamspotx;
	double fTbeamspoty;
	double fTbeamspotz;

	int fTNCaloTowers;

// Trigger
	int fTtrigres[200];

// Flags
	int fTgoodevent;      // 1 for good events, 0 for bad events
	int fTflagmaxmuexc;   // Found more than 20 muons in event (0 is good, 1 is bad)
	int fTflagmaxelexc;   // Found more than 20 electrons in event
	int fTflagmaxujetexc; // Found more than 50 jets in event
	int fTflagmaxjetexc;  // Found more than 50 uncorrected jets in event
	int fTflagmaxtrkexc;  // Found more than 500 tracks in event

// Muons:
	unsigned int fTnmu;
	int fTgoodmu[fMaxnmus];
	double fTmupx[fMaxnmus];
	double fTmupy[fMaxnmus];
	double fTmupz[fMaxnmus];
	double fTmue[fMaxnmus];
	double fTmuet[fMaxnmus];
	double fTmupt[fMaxnmus];
	double fTmueta[fMaxnmus];
	double fTmuphi[fMaxnmus];
	int fTmucharge[fMaxnmus];

// - Isolation Variables
	double fTmuetsum[fMaxnmus];
	double fTmuptsum[fMaxnmus];
	double fTmuiso[fMaxnmus];
	double fTmueecal[fMaxnmus];
	double fTmuehcal[fMaxnmus];

// - Impact Parameters
	double fTmud0bs[fMaxnmus];
	double fTmud0pv[fMaxnmus];
	double fTmud0E[fMaxnmus];
	double fTmudzbs[fMaxnmus];
	double fTmudzpv[fMaxnmus];
	double fTmudzE[fMaxnmus];

// - MuID Variables
	double fTmunchi2[fMaxnmus];
	int fTmunglhits[fMaxnmus];
	int fTmunmuhits[fMaxnmus];
	int fTmuntkhits[fMaxnmus];
	int fTmunmatches[fMaxnmus];
	int fTmunchambers[fMaxnmus];
	double fTmucalocomp[fMaxnmus];
	double fTmusegmcomp[fMaxnmus];
	int fTmutrackermu[fMaxnmus];
	int fTmuisGMPT[fMaxnmus];

// - Gen Info:
	int fTmuid[fMaxnmus];
	int fTmumid[fMaxnmus];

// Electrons:
	int fTneles;
	int fTgoodel[fMaxneles];
	double fTepx[fMaxneles];
	double fTepy[fMaxneles];
	double fTepz[fMaxneles];
	double fTept[fMaxneles];
	double fTee[fMaxneles];
	double fTeet[fMaxneles];
	double fTeeta[fMaxneles];
	double fTephi[fMaxneles];
	double fTed0bs[fMaxneles];
	double fTed0pv[fMaxneles];
	double fTed0E[fMaxneles];
	double fTedzbs[fMaxneles];
	double fTedzpv[fMaxneles];
	double fTedzE[fMaxneles];
	double fTeiso[fMaxneles];
	double fTeptsum[fMaxneles];
	double fTeetsum[fMaxneles];
	double fTenchi2[fMaxneles];
	int fTeIDTight[fMaxneles];
	int fTeIDLoose[fMaxneles];
	int fTeIDRobustTight[fMaxneles];
	int fTeIDRobustLoose[fMaxneles];
	int fTecharge[fMaxneles];
	int fTeInGap[fMaxneles];  // seed crystal next to a gap
	int fTeEcalDriven[fMaxneles];
	int fTeTrackerDriven[fMaxneles];
	int fTeBasicClustersSize[fMaxneles];
	double fTefbrem[fMaxneles];
	double fTeHcalOverEcal[fMaxneles];
	double fTeE5x5[fMaxneles];                      // 5x5 arround seed
	double fTeE2x5Max[fMaxneles];                   // 2x5 arround seed
	double fTeSigmaIetaIeta[fMaxneles];             // shower shape covariance
	double fTeDeltaPhiSeedClusterAtCalo[fMaxneles]; // Dphi (seed-track) at calo from p_out
	double fTeDeltaEtaSeedClusterAtCalo[fMaxneles]; // outermost track state extrapolated at calo
	double fTeDeltaPhiSuperClusterAtVtx[fMaxneles]; // Dphi (sc-track) at calo extrapolated from p_in
	double fTeDeltaEtaSuperClusterAtVtx[fMaxneles]; // Deta (sc-track) at calo extrapolated from p_in
	double fTecaloenergy[fMaxneles];                // caloEnergy() = supercluster energy 99.9% of the time
	double fTtrkmomatvtx[fMaxneles];                // trackMomentumAtVtx().R()
	double fTeESuperClusterOverP[fMaxneles];        // Esc/Pin


// Jets:
	int fTnjets;
	int fTgoodjet[fMaxnjets];
	double fTjpx[fMaxnjets];
	double fTjpy[fMaxnjets];
	double fTjpz[fMaxnjets];
	double fTje[fMaxnjets];
	double fTjet[fMaxnjets];
	double fTjpt[fMaxnjets];
	double fTjeta[fMaxnjets];
	double fTjphi[fMaxnjets];
	double fTjemfrac[fMaxnjets];
	double fTjID_HPD[fMaxnjets];
	double fTjID_RBX[fMaxnjets];
	double fTjID_n90Hits[fMaxnjets];
	double fTjID_SubDet1[fMaxnjets];
	double fTjID_SubDet2[fMaxnjets];
	double fTjID_SubDet3[fMaxnjets];
	double fTjID_SubDet4[fMaxnjets];
	double fTjID_resEMF[fMaxnjets];
	double fTjID_HCALTow[fMaxnjets];
	double fTjID_ECALTow[fMaxnjets];
	double fTJEtaEMrms[fMaxnjets];
	double fTJEtaHADrms[fMaxnjets];
	double fTJPhiEMrms[fMaxnjets];
	double fTJPhiHADrms[fMaxnjets];
	double fTbTagProb[fMaxnjets];
	double fTChfrac[fMaxnjets];
	int fTnAssoTracks[fMaxnjets];
	double fTtrk1px[fMaxnjets];
	double fTtrk1py[fMaxnjets];
	double fTtrk1pz[fMaxnjets];
	double fTtrk2px[fMaxnjets];
	double fTtrk2py[fMaxnjets];
	double fTtrk2pz[fMaxnjets];
	double fTtrk3px[fMaxnjets];
	double fTtrk3py[fMaxnjets];
	double fTtrk3pz[fMaxnjets];
	double fTjEcorr[fMaxnjets];
	double fTjeMinDR[fMaxnjets];
	double fTjetVtxx[fMaxnjets];
	double fTjetVtxy[fMaxnjets];
	double fTjetVtxz[fMaxnjets];
	double fTjetVtxExx[fMaxnjets];
	double fTjetVtxEyx[fMaxnjets];
	double fTjetVtxEyy[fMaxnjets];
	double fTjetVtxEzy[fMaxnjets];
	double fTjetVtxEzz[fMaxnjets];
	double fTjetVtxEzx[fMaxnjets];
	double fTjetVtxNChi2[fMaxnjets];

// Tracks:
	int fTntracks;
	int fTgoodtrk[fMaxntrks];
	double fTtrkpt[fMaxntrks];
	double fTtrketa[fMaxntrks];
	double fTtrkphi[fMaxntrks];
	double fTtrknchi2[fMaxntrks];
	double fTtrknhits[fMaxntrks];

// (M)E(T):
	double fTTrkPtSumx;
	double fTTrkPtSumy;
	double fTTrkPtSumphi;
	double fTTrkPtSum;
	double fTSumEt;
	double fTECALSumEt;
	double fTHCALSumEt;
	double fTECALEsumx;
	double fTECALEsumy;
	double fTECALEsumz;
	double fTECALMETphi;
	double fTECALMETeta;
	double fTECALMET;
	double fTHCALEsumx;
	double fTHCALEsumy;
	double fTHCALEsumz;
	double fTHCALMETphi;
	double fTHCALMETeta;
	double fTHCALMET;
	double fTRawMET;
	double fTRawMETpx;
	double fTRawMETpy;
	double fTRawMETphi;
	double fTMuCorrMET;
	double fTMuCorrMETpx;
	double fTMuCorrMETpy;
	double fTMuCorrMETphi;
	double fTTCMET;
	double fTTCMETpx;
	double fTTCMETpy;
	double fTTCMETphi;
	double fTMuJESCorrMET;
	double fTMuJESCorrMETpx;
	double fTMuJESCorrMETpy;
	double fTMuJESCorrMETphi;
	double fTPFMET;
	double fTPFMETpx;
	double fTPFMETpy;
	double fTPFMETphi;
////////////////////////////////////////////////////////
};
