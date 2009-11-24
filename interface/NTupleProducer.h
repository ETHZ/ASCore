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
// $Id: NTupleProducer.h,v 1.20 2009/11/19 17:19:34 sordini Exp $
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
	double fMintrknhits;
	double fIso_MuTkDRin;
	double fIso_MuTkDRout;
	double fIso_MuTkSeed;
	double fIso_MuCalDRin;
	double fIso_MuCalDRout;
	double fIso_MuCalSeed;

	double fJUNC_px_match[50];
	double fJUNC_py_match[50];
	double fJUNC_pz_match[50];

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
	int fTgoodmu[20];
	double fTmupx[20];
	double fTmupy[20];
	double fTmupz[20];
	double fTmue[20];
	double fTmuet[20];
	double fTmupt[20];
	double fTmueta[20];
	double fTmuphi[20];
	int fTmucharge[20];

// - Isolation Variables
	double fTmuetsum[20];
	double fTmuptsum[20];
	double fTmuiso[20];
	double fTmueecal[20];
	double fTmuehcal[20];

// - Impact Parameters
	double fTmud0bs[20];
	double fTmud0pv[20];
	double fTmud0E[20];
	double fTmudzbs[20];
	double fTmudzpv[20];
	double fTmudzE[20];

// - MuID Variables
	double fTmunchi2[20];
	int fTmunglhits[20];
	int fTmunmuhits[20];
	int fTmuntkhits[20];
	int fTmunmatches[20];
	int fTmunchambers[20];
	double fTmucalocomp[20];
	double fTmusegmcomp[20];
	int fTmutrackermu[20];
	int fTmuisGMPT[20];

// - Gen Info:
	int fTmuid[20];
	int fTmumid[20];

// Electrons:
	int fTneles;
	int fTgoodel[20];
	double fTepx[20];
	double fTepy[20];
	double fTepz[20];
	double fTept[20];
	double fTee[20];
	double fTeet[20];
	double fTeeta[20];
	double fTephi[20];
	double fTed0bs[20];
	double fTed0pv[20];
	double fTed0E[20];
	double fTedzbs[20];
	double fTedzpv[20];
	double fTedzE[20];
	double fTeiso[20];
	double fTeptsum[20];
	double fTeetsum[20];
	double fTenchi2[20];
	int fTeIDTight[20];
	int fTeIDLoose[20];
	int fTeIDRobustTight[20];
	int fTeIDRobustLoose[20];
	int fTecharge[20];
	int fTeInGap[20];  // seed crystal next to a gap
	int fTeEcalDriven[20];
	int fTeTrackerDriven[20];
	int fTeBasicClustersSize[20];
	double fTefbrem[20];
	double fTeHcalOverEcal[20];
	double fTeE5x5[20];                      // 5x5 arround seed
	double fTeE2x5Max[20];                   // 2x5 arround seed
	double fTeSigmaIetaIeta[20];             // shower shape covariance
	double fTeDeltaPhiSeedClusterAtCalo[20]; // Dphi (seed-track) at calo from p_out
	double fTeDeltaEtaSeedClusterAtCalo[20]; // outermost track state extrapolated at calo
	double fTeDeltaPhiSuperClusterAtVtx[20]; // Dphi (sc-track) at calo extrapolated from p_in
	double fTeDeltaEtaSuperClusterAtVtx[20]; // Deta (sc-track) at calo extrapolated from p_in
	double fTecaloenergy[20];                // caloEnergy() = supercluster energy 99.9% of the time
	double fTtrkmomatvtx[20];                // trackMomentumAtVtx().R()
	double fTeESuperClusterOverP[20];        // Esc/Pin


// Jets:
	int fTnjets;
	int fTgoodjet[50];
	double fTjpx[50];
	double fTjpy[50];
	double fTjpz[50];
	double fTje[50];
	double fTjet[50];
	double fTjpt[50];
	double fTjeta[50];
	double fTjphi[50];
	double fTjemfrac[50];
	double fTjID_HPD[50];
	double fTjID_RBX[50];
	double fTjID_n90Hits[50];
	double fTjID_SubDet1[50];
	double fTjID_SubDet2[50];
	double fTjID_SubDet3[50];
	double fTjID_SubDet4[50];
	double fTjID_resEMF[50];
	double fTjID_HCALTow[50];
	double fTjID_ECALTow[50];
	double fTJEtaEMrms[50];
	double fTJEtaHADrms[50];
	double fTJPhiEMrms[50];
	double fTJPhiHADrms[50];
	double fTbTagProb[50];
	double fTChfrac[50];
	int fTnAssoTracks[50];
	double fTtrk1px[50];
	double fTtrk1py[50];
	double fTtrk1pz[50];
	double fTtrk2px[50];
	double fTtrk2py[50];
	double fTtrk2pz[50];
	double fTtrk3px[50];
	double fTtrk3py[50];
	double fTtrk3pz[50];
	double fTjEcorr[50];
	double fTjeMinDR[50];
	double fTjetVtxx[50];
	double fTjetVtxy[50];
	double fTjetVtxz[50];
	double fTjetVtxExx[50];
	double fTjetVtxEyx[50];
	double fTjetVtxEyy[50];
	double fTjetVtxEzy[50];
	double fTjetVtxEzz[50];
	double fTjetVtxEzx[50];
	double fTjetVtxNChi2[50];

// Tracks:
	int fTntracks;
	int fTgoodtrk[500];
	double fTtrkpt[500];
	double fTtrketa[500];
	double fTtrkphi[500];
	double fTtrknchi2[500];
	double fTtrknhits[500];

// MET:
	double fTTrkPtSumx;
	double fTTrkPtSumy;
	double fTTrkPtSumphi;
	double fTTrkPtSum;
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
