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
// $Id: NTupleProducer.h,v 1.27 2009/12/11 20:02:23 stiegerb Exp $
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

	virtual void beginJob(const edm::EventSetup&);
	virtual void beginRun(const edm::Run&, const edm::EventSetup&);
	virtual void analyze(const edm::Event&, const edm::EventSetup&);
	virtual void endJob();
	virtual void endRun(const edm::Run&, const edm::EventSetup&);

	vector<double> calcMuIso(const reco::Muon *Mu, const edm::Event& iEvent);
	vector<double> calcMuIso2(const reco::Muon *Mu, const edm::Event& iEvent);
	vector<double> calcMuIso3(const reco::Muon *Mu, const edm::Event& iEvent, edm::Ref<reco::MuonCollection> muonRef);
	vector<double> calcElIso(const reco::GsfElectron *El, const edm::Event& iEvent);
	vector<int> matchMuCand(const reco::Muon *Mu, const edm::Event& iEvent);
	vector<const reco::Muon*> sortMus(vector<const reco::Muon*>);
	void switchDouble(double &, double &);
	void switchInt(int &, int &);
	void resetDouble(double *v, unsigned int size = 20);
	void resetInt(int *v, unsigned int size = 20);
	void resetTree();
	vector<const reco::Track*> FindAssociatedTracks(const reco::Jet *jet, const reco::TrackCollection *tracks);
private:

	virtual void ElectronDuplicate(vector<const SuperCluster*> elecPtr, vector<const GsfTrack*> trckPtr);
	virtual void ElJetOverlap(vector<const Jet*> jets, vector<const SuperCluster*> electrons, edm::Handle<CaloTowerCollection> calotowers);
	virtual bool IsEMObjectInJet(const SuperCluster* theElecSC, const CaloJet* theJet, edm::Handle<CaloTowerCollection> calotowers, math::XYZVector* sharedMomentum);
	virtual bool EMCaloTowerWindow(const SuperCluster* superCluster, float & phimin, float & phimax, float & etamin, float & etamax);
	virtual float CaloTowerSizePhi(float eta);
	virtual float CaloTowerSizeEta(float eta);
	virtual bool IsInPhiWindow(float phi, float phimin, float phimax);
	virtual float DeltaPhiSigned(float v1, float v2);
	virtual float GetPhiMin(float phi1, float phi2);
	virtual float GetPhiMax(float phi1, float phi2);

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

	static const int gMaxnmus  = 20;
	static const int gMaxneles = 20;
	static const int gMaxnjets = 100;
	static const int gMaxntrks = 500;
	static const int gMaxnphos = 500;

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
	edm::InputTag fPhotonTag;
	edm::InputTag fCalTowTag;
	edm::InputTag fGenPartTag;
	edm::InputTag fL1TriggerTag;
	edm::InputTag fHLTTriggerTag;

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

	double fMinphopt;
	double fMaxphoeta;

	double fJUNC_px_match[gMaxnjets];
	double fJUNC_py_match[gMaxnjets];
	double fJUNC_pz_match[gMaxnjets];

	TH1I *fHhltstat; // Added to keep track of trigger names
	TH1I *fHl1physstat;
	TH1I *fHl1techstat;
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

	double fRTMinphopt;
	double fRTMaxphoeta;

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
	double fTpvtxptsum;

	double fTbeamspotx;
	double fTbeamspoty;
	double fTbeamspotz;

	int fTNCaloTowers;

// Trigger
	static const unsigned int gMaxhltbits = 200;
	static const unsigned int gMaxl1physbits = 128;
	static const unsigned int gMaxl1techbits = 64;
	int fTHLTres[200];
	int fTL1physres[gMaxl1physbits];
	int fTL1techres[gMaxl1techbits];

// Flags
	int fTgoodevent;      // 1 for good events, 0 for bad events
	int fTflagmaxmuexc;   // Found more than 20 muons in event (0 is good, 1 is bad)
	int fTflagmaxelexc;   // Found more than 20 electrons in event
	int fTflagmaxujetexc; // Found more than 50 jets in event
	int fTflagmaxjetexc;  // Found more than 50 uncorrected jets in event
	int fTflagmaxtrkexc;  // Found more than 500 tracks in event
	int fTflagmaxphoexc;  // Found more than 500 photons in event

// Muons:
	unsigned int fTnmu;
	int fTgoodmu[gMaxnmus];
	double fTmupx[gMaxnmus];
	double fTmupy[gMaxnmus];
	double fTmupz[gMaxnmus];
	double fTmue[gMaxnmus];
	double fTmuet[gMaxnmus];
	double fTmupt[gMaxnmus];
	double fTmuptE[gMaxnmus];
	double fTmueta[gMaxnmus];
	double fTmuphi[gMaxnmus];
	int fTmucharge[gMaxnmus];

// - Isolation Variables
	double fTmuetsum[gMaxnmus];
	double fTmuptsum[gMaxnmus];
	double fTmuiso[gMaxnmus];
	double fTmueecal[gMaxnmus];
	double fTmuehcal[gMaxnmus];

// - Impact Parameters
	double fTmud0bs[gMaxnmus];
	double fTmud0pv[gMaxnmus];
	double fTmud0E[gMaxnmus];
	double fTmudzbs[gMaxnmus];
	double fTmudzpv[gMaxnmus];
	double fTmudzE[gMaxnmus];

// - MuID Variables
	double fTmunchi2[gMaxnmus];
	int fTmunglhits[gMaxnmus];
	int fTmunmuhits[gMaxnmus];
	int fTmuntkhits[gMaxnmus];
	int fTmunmatches[gMaxnmus];
	int fTmunchambers[gMaxnmus];
	double fTmucalocomp[gMaxnmus];
	double fTmusegmcomp[gMaxnmus];
	int fTmutrackermu[gMaxnmus];
	int fTmuisGMPT[gMaxnmus];

// - Gen Info:
	int fTmuid[gMaxnmus];
	int fTmumid[gMaxnmus];

// Electrons:
	int fTneles;
	int fTgoodel[gMaxneles];
	double fTepx[gMaxneles];
	double fTepy[gMaxneles];
	double fTepz[gMaxneles];
	double fTept[gMaxneles];
	double fTeptE[gMaxneles];
	double fTee[gMaxneles];
	double fTeet[gMaxneles];
	double fTeeta[gMaxneles];
	double fTephi[gMaxneles];
	double fTed0bs[gMaxneles];
	double fTed0pv[gMaxneles];
	double fTed0E[gMaxneles];
	double fTedzbs[gMaxneles];
	double fTedzpv[gMaxneles];
	double fTedzE[gMaxneles];
	double fTeiso[gMaxneles];
	double fTeptsum[gMaxneles];
	double fTeetsum[gMaxneles];
	double fTenchi2[gMaxneles];
	int fTeIDTight[gMaxneles];
	int fTeIDLoose[gMaxneles];
	int fTeIDRobustTight[gMaxneles];
	int fTeIDRobustLoose[gMaxneles];
	int fTecharge[gMaxneles];
	int fTeInGap[gMaxneles];  // seed crystal next to a gap
	int fTeEcalDriven[gMaxneles];
	int fTeTrackerDriven[gMaxneles];
	int fTeBasicClustersSize[gMaxneles];
	double fTefbrem[gMaxneles];
	double fTeHcalOverEcal[gMaxneles];
	double fTeE5x5[gMaxneles];                      // 5x5 arround seed
	double fTeE2x5Max[gMaxneles];                   // 2x5 arround seed
	double fTeSigmaIetaIeta[gMaxneles];             // shower shape covariance
	double fTeDeltaPhiSeedClusterAtCalo[gMaxneles]; // Dphi (seed-track) at calo from p_out
	double fTeDeltaEtaSeedClusterAtCalo[gMaxneles]; // outermost track state extrapolated at calo
	double fTeDeltaPhiSuperClusterAtVtx[gMaxneles]; // Dphi (sc-track) at calo extrapolated from p_in
	double fTeDeltaEtaSuperClusterAtVtx[gMaxneles]; // Deta (sc-track) at calo extrapolated from p_in
	double fTecaloenergy[gMaxneles];                // caloEnergy() = supercluster energy 99.9% of the time
	double fTetrkmomatvtx[gMaxneles];               // trackMomentumAtVtx().R()
	double fTeESuperClusterOverP[gMaxneles];        // Esc/Pin
	int fTeIsInJet[gMaxneles];
	double fTeSharedPx[gMaxneles];
	double fTeSharedPy[gMaxneles];
	double fTeSharedPz[gMaxneles];
	double fTeSharedEnergy[gMaxneles];
	int fTeDupEl[gMaxneles];

// Jets:
	int fTnjets;
	int fTgoodjet[gMaxnjets];
	double fTjpx[gMaxnjets];
	double fTjpy[gMaxnjets];
	double fTjpz[gMaxnjets];
	double fTje[gMaxnjets];
	double fTjet[gMaxnjets];
	double fTjpt[gMaxnjets];
	double fTjeta[gMaxnjets];
	double fTjphi[gMaxnjets];
	double fTjemfrac[gMaxnjets];
	int fTjNconstituents[gMaxnjets];
	double fTjID_HPD[gMaxnjets];
	double fTjID_RBX[gMaxnjets];
	double fTjID_n90Hits[gMaxnjets];
	double fTjID_SubDet1[gMaxnjets];
	double fTjID_SubDet2[gMaxnjets];
	double fTjID_SubDet3[gMaxnjets];
	double fTjID_SubDet4[gMaxnjets];
	double fTjID_resEMF[gMaxnjets];
	double fTjID_HCALTow[gMaxnjets];
	double fTjID_ECALTow[gMaxnjets];
	double fTJEtaEMrms[gMaxnjets];
	double fTJEtaHADrms[gMaxnjets];
	double fTJPhiEMrms[gMaxnjets];
	double fTJPhiHADrms[gMaxnjets];
	double fTbTagProb[gMaxnjets];
	double fTChfrac[gMaxnjets];
	int fTnAssoTracks[gMaxnjets];
	double fTtrk1px[gMaxnjets];
	double fTtrk1py[gMaxnjets];
	double fTtrk1pz[gMaxnjets];
	double fTtrk2px[gMaxnjets];
	double fTtrk2py[gMaxnjets];
	double fTtrk2pz[gMaxnjets];
	double fTtrk3px[gMaxnjets];
	double fTtrk3py[gMaxnjets];
	double fTtrk3pz[gMaxnjets];
	double fTjEcorr[gMaxnjets];
	double fTjeMinDR[gMaxnjets];
	double fTjetVtxx[gMaxnjets];
	double fTjetVtxy[gMaxnjets];
	double fTjetVtxz[gMaxnjets];
	double fTjetVtxExx[gMaxnjets];
	double fTjetVtxEyx[gMaxnjets];
	double fTjetVtxEyy[gMaxnjets];
	double fTjetVtxEzy[gMaxnjets];
	double fTjetVtxEzz[gMaxnjets];
	double fTjetVtxEzx[gMaxnjets];
	double fTjetVtxNChi2[gMaxnjets];

// Tracks:
	int fTntracks;
	int fTgoodtrk[gMaxntrks];
	double fTtrkpt[gMaxntrks]; // this is actually charge*pt
	double fTtrketa[gMaxntrks];
	double fTtrkphi[gMaxntrks];
	double fTtrknchi2[gMaxntrks];
	double fTtrknhits[gMaxntrks];

  //Photons
  int    fTnphotons;
  double fTPhotonPt[gMaxnphos];  
  double fTPhotonPx[gMaxnphos];  
  double fTPhotonPy[gMaxnphos];  
  double fTPhotonPz[gMaxnphos];  
  double fTPhotonEta[gMaxnphos];
  double fTPhotonPhi[gMaxnphos];
  double fTPhotonEnergy[gMaxnphos];
  double fTPhotoncaloPositionX[gMaxnphos];
  double fTPhotoncaloPositionY[gMaxnphos];
  double fTPhotoncaloPositionZ[gMaxnphos];
  double fTPhotonHoverE[gMaxnphos];
  double fTPhotonH1overE[gMaxnphos];
  double fTPhotonH2overE[gMaxnphos];
  int    fTPhotonHasPixSeed[gMaxnphos];
  int    fTPhotonHasConvTrks[gMaxnphos];
  
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
