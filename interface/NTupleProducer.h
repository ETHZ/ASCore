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
// $Id: NTupleProducer.h,v 1.92 2011/05/17 14:04:04 fronga Exp $
//
//


// system include files
#include <vector>
#include <string>

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
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

// Data formats
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"


// Helpers
#include "Math/VectorUtil.h"

// Local classes
#include "DiLeptonAnalysis/NTupleProducer/interface/JetFillerReco.h"
#include "DiLeptonAnalysis/NTupleProducer/interface/JetFillerPat.h"
#include "DiLeptonAnalysis/NTupleProducer/interface/LeptonFillerPat.h"

typedef math::XYZTLorentzVector LorentzVector;
using namespace reco;

class NTupleProducer : public edm::EDAnalyzer {
public:
	explicit NTupleProducer(const edm::ParameterSet&);
	~NTupleProducer();

	virtual void beginJob(); //336  beginJob(const edm::EventSetup&)
	virtual void beginRun(const edm::Run&, const edm::EventSetup&);
	virtual void analyze(const edm::Event&, const edm::EventSetup&);
	virtual void endJob();
	virtual void endRun(const edm::Run&, const edm::EventSetup&);
	std::vector<const reco::GenParticle*> matchRecoCand(const reco::RecoCandidate *Cand, const edm::Event& iEvent);
	const int matchJet(const reco::Jet *jet);
	void resetDouble(double *v, unsigned int size = 20);
	void resetFloat(float *v, unsigned int size = 20);
	void resetInt(int *v, unsigned int size = 20);
	void resetTree();

private:

	virtual void ElectronDuplicate(std::vector<const SuperCluster*> elecPtr, std::vector<const GsfTrack*> trckPtr);
	virtual void PhotonElectronDuplicate(std::vector<const SuperCluster*>, std::vector<const SuperCluster*>);
	virtual void ElJetOverlap(std::vector<const Jet*> jets, std::vector<const SuperCluster*> electrons, edm::Handle<CaloTowerCollection> calotowers);
	virtual void PhotonJetOverlap(std::vector<const Jet*> jets, std::vector<const SuperCluster*> electrons, edm::Handle<CaloTowerCollection> calotowers);
	virtual bool IsEMObjectInJet(const SuperCluster* theElecSC, std::vector<CaloTowerPtr> jetCaloRefs, edm::Handle<CaloTowerCollection> calotowers, math::XYZVector* sharedMomentum);
	virtual bool EMCaloTowerWindow(const SuperCluster* superCluster, float & phimin, float & phimax, float & etamin, float & etamax);
	virtual float CaloTowerSizePhi(float eta);
	virtual float CaloTowerSizeEta(float eta);
	virtual bool IsInPhiWindow(float phi, float phimin, float phimax);
	virtual float DeltaPhiSigned(float v1, float v2);
	virtual float GetPhiMin(float phi1, float phi2);
	virtual float GetPhiMax(float phi1, float phi2);

        typedef std::pair<int,double> OrderPair;
	struct IndexByPt {
		const bool operator()(const OrderPair& j1, const OrderPair& j2 ) const {
			return j1.second > j2.second;
		}
	};


// ----------member data ---------------------------
	edm::Service<TFileService> fTFileService;
        AdaptiveVertexFitter avFitter;


	std::vector<JetFillerBase*>     jetFillers;
        std::vector<PatMuonFiller*>     muonFillers;
        std::vector<PatElectronFiller*> electronFillers;
        std::vector<PatTauFiller*>      tauFillers;

	bool fIsRealData;
	int fNTotEvents;
	int fNFillTree;

	static const int gMaxnmus     = 30;
	static const int gMaxneles    = 20;
	static const int gMaxntaus    = 20;
	static const int gMaxnjets    = 100;
	static const int gMaxntrks    = 500;
	static const int gMaxnphos    = 50;
	static const int gMaxngenlept = 100;
	static const int gMaxngenjets = 100;
	static const int gMaxnvrtx    = 25;
	static const int gMaxnpileup  = 50;
	static const int gMaxnEBhits  = 20;

	edm::InputTag fMuonTag;
	edm::InputTag fElectronTag;
        std::string fEleIdWP;
	edm::InputTag fMuIsoDepTkTag;
	edm::InputTag fMuIsoDepECTag;
	edm::InputTag fMuIsoDepHCTag;
	edm::InputTag fJetTag;
        std::string fJetCorrs;
	edm::InputTag fBtag1Tag;
	edm::InputTag fBtag2Tag;
	edm::InputTag fBtag3Tag;
	edm::InputTag fBtag4Tag;
	edm::InputTag fRawCaloMETTag;
	edm::InputTag fTCMETTag;
	edm::InputTag fPFMETTag;
	edm::InputTag fPFMETPATTag;
	edm::InputTag fCorrCaloMETTag;
	edm::InputTag fGenMETTag;
	edm::InputTag fVertexTag;
	edm::InputTag fTrackTag;
	edm::InputTag fPhotonTag;
	edm::InputTag fCalTowTag;
	edm::InputTag fEBRecHitsTag;
	edm::InputTag fEERecHitsTag;
	edm::InputTag fGenPartTag;
	edm::InputTag fGenJetTag;
	edm::InputTag fL1TriggerTag;
	edm::InputTag fHLTTrigEventTag;
	edm::InputTag fHBHENoiseResultTag;
	edm::InputTag fSrcRho;

	float fMinmupt;
	float fMaxmueta;
	float fMinelpt;
	float fMaxeleta;
	float fMincorjpt;
	float fMinrawjpt;
	float fMaxjeta;
	float fMinjemfrac;
	float fMintrkpt;
	float fMaxtrketa;
	float fMaxtrknchi2;
	int	fMintrknhits;
	float fMinphopt;
	float fMaxphoeta;

	float fMingenleptpt; 
	float fMaxgenlepteta;
	float fMingenjetpt;
	float fMaxgenjeteta;
	
	float fBtagMatchdeltaR;

	TH1I *fHhltstat; // Added to keep track of trigger names
	TH1I *fHl1physstat;
	TH1I *fHl1techstat;
	bool fFirstevent;

        // Trigger stuff
        std::string fProcessName; // process name of (HLT) process for which to get HLT configuration
	HLTConfigProvider fHltConfig;

        ////////////////////////////////////////////////////////
        // Trees:
	TTree *fRunTree;

	int fRTrunnumber;
	float fRTextxslo;
	float fRTextxsnlo;
	float fRTintxs;

	float fRTMinmupt;
	float fRTMaxmueta;
	float fRTMinelpt;
	float fRTMaxeleta;
	float fRTMinjpt;
	float fRTMinrawjpt;
	float fRTMaxjeta;
	float fRTMinjemfrac;

	float fRTMintrkpt;
	float fRTMaxtrketa;
	float fRTMaxtrknchi2;
	int fRTMintrknhits;

	float fRTMinphopt;
	float fRTMaxphoeta;


	int fRTmaxnmu;
	int fRTmaxnel;
	int fRTmaxnjet;
	int fRTmaxntrk;
	int fRTmaxnphot;

	float fMinebrechitE; 

	TTree *fEventTree;

	// General event information
	int fTrunnumber;
	int fTeventnumber;
	int fTlumisection;
	float fTpthat;
	int fTsigprocid;
	float fTpdfscalePDF;
	int fTpdfid1;
	int fTpdfid2;
	float fTpdfx1;
	float fTpdfx2;
	float fTpdfxPDF1;
	float fTpdfxPDF2;
	float fTextxslo;
	float fTintxs;
	float fTweight;

	// Pile-up
	int fTpuNumInteractions;
	float fTpuZpositions[gMaxnpileup];
	float fTpuSumpT_lowpT[gMaxnpileup];
	float fTpuSumpT_highpT[gMaxnpileup];
	float fTpuNtrks_lowpT[gMaxnpileup];
	float fTpuNtrks_highpT[gMaxnpileup];
	float fTrho; // rho from L1FastJetCorrection
	// float fTpuInstLumi[gMaxnpileup];
        TH1I* fHpileupstat;

	// ECAL & HCAL Noise
	int fTHBHENoiseFlag;

	int fTnEBhits;
	float fTEBrechitE[gMaxnEBhits];
	float fTEBrechitPt[gMaxnEBhits];
	float fTEBrechitEta[gMaxnEBhits];
	float fTEBrechitPhi[gMaxnEBhits];
	float fTEBrechitChi2[gMaxnEBhits];
	float fTEBrechitTime[gMaxnEBhits];
	float fTEBrechitE4oE1[gMaxnEBhits];
	float fTEBrechitE2oE9[gMaxnEBhits];
	
	// int fTEcalDeadCellBEFlag;
	// static const unsigned int gMaxnECALGapClusters = 50;
	// unsigned int fTnECALGapClusters;
	// float fTEcalGapBE[gMaxnECALGapClusters];
	// int fTEcalGapClusterSize[gMaxnECALGapClusters];

	int fTgoodvtx;
	float fTprimvtxx;
	float fTprimvtxy;
	float fTprimvtxz;
	float fTprimvtxrho;
	float fTprimvtxxE;
	float fTprimvtxyE;
	float fTprimvtxzE;
	float fTpvtxznchi2;
	float fTpvtxndof;
	int   fTpvtxisfake;
	float fTpvtxptsum;

	float fTbeamspotx;
	float fTbeamspoty;
	float fTbeamspotz;

	int fTnvrtx;
	float fTvrtxx[gMaxnvrtx];
	float fTvrtxy[gMaxnvrtx];
	float fTvrtxz[gMaxnvrtx];
	float fTvrtxxE[gMaxnvrtx];
	float fTvrtxyE[gMaxnvrtx];
	float fTvrtxzE[gMaxnvrtx];
	float fTvrtxndof[gMaxnvrtx];
	float fTvrtxchi2[gMaxnvrtx];
	float fTvrtxntrks[gMaxnvrtx];
	float fTvrtxsumpt[gMaxnvrtx];
	int fTvrtxisfake[gMaxnvrtx];

	int fTNCaloTowers;

// Trigger
	static const unsigned int gMaxhltbits = 400;
	static const unsigned int gMaxl1physbits = 128;
	static const unsigned int gMaxl1techbits = 64;
	int fTHLTres[gMaxhltbits];
	int fTL1physres[gMaxl1physbits];
	int fTL1techres[gMaxl1techbits];
	int fTHLTprescale[gMaxhltbits];
	std::vector<std::string> fTHLTmenu;
	std::vector<std::string> fTL1physmenu;

	static const unsigned int gMaxhltnobjs  = 10;
	std::vector<std::string> fTHltLabels; // HLT Paths to store the triggering objects of
	unsigned int fTNpaths;
	unsigned int fTNHLTobjects;
	int**    fTHLTObjectID;
	float** fTHLTObjectPt;
	float** fTHLTObjectEta;
	float** fTHLTObjectPhi;

// Flags
	int fTgoodevent;         // 0 for good events, 1 for bad events
	int fTflagmaxmuexc;      // Found more than 20 muons in event (0 is good, 1 is bad)
	int fTflagmaxelexc;      // Found more than 20 electrons in event
	int fTflagmaxujetexc;    // Found more than 50 jets in event
	int fTflagmaxjetexc;     // Found more than 50 uncorrected jets in event
	int fTflagmaxtrkexc;     // Found more than 500 tracks in event
	int fTflagmaxphoexc;     // Found more than 500 photons in event
	int fTflagmaxgenleptexc; // Found more than 100 genleptons in event
	int fTflagmaxgenjetexc;  // Found more than 100 genjets in event
	int fTflagmaxvrtxexc;    // Found more than 25 vertices in event

// GenLeptons
	int   fTngenleptons;      
	int   fTGenLeptonId[gMaxngenlept];      
	float fTGenLeptonPt[gMaxngenlept];      
	float fTGenLeptonEta[gMaxngenlept];     
	float fTGenLeptonPhi[gMaxngenlept];     
	int   fTGenLeptonMId[gMaxngenlept];     
	int   fTGenLeptonMStatus[gMaxngenlept]; 
	float fTGenLeptonMPt[gMaxngenlept];     
	float fTGenLeptonMEta[gMaxngenlept];    
	float fTGenLeptonMPhi[gMaxngenlept];    
	int   fTGenLeptonGMId[gMaxngenlept];    
	int   fTGenLeptonGMStatus[gMaxngenlept];
	float fTGenLeptonGMPt[gMaxngenlept];    
	float fTGenLeptonGMEta[gMaxngenlept];   
	float fTGenLeptonGMPhi[gMaxngenlept]; 

// GenJets
	int   fTNGenJets;
	float fTGenJetPt  [gMaxngenjets];
	float fTGenJetEta [gMaxngenjets];
	float fTGenJetPhi [gMaxngenjets];
	float fTGenJetE   [gMaxngenjets];
	float fTGenJetemE [gMaxngenjets];
	float fTGenJethadE[gMaxngenjets];
	float fTGenJetinvE[gMaxngenjets];


// Muons:
	int fTnmu;
	int fTnmutot; // before preselection
	int fTnglobalmu;
	int fTntrackermu;
	int fTgoodmu     [gMaxnmus];
	int fTmuIsIso    [gMaxnmus];
	int fTmuIsGM     [gMaxnmus];
	int fTmuIsTM     [gMaxnmus];
	float fTmupx     [gMaxnmus];
	float fTmupy     [gMaxnmus];
	float fTmupz     [gMaxnmus];
	float fTmue      [gMaxnmus];
	float fTmuet     [gMaxnmus];
	float fTmupt     [gMaxnmus];
	float fTmuinnerpt[gMaxnmus];
	float fTmuptE    [gMaxnmus];
	float fTmueta    [gMaxnmus];
	float fTmuphi    [gMaxnmus];
	int fTmucharge   [gMaxnmus];

// - Isolation Variables
	float fTmuiso           [gMaxnmus];
	float fTmuIso03sumPt    [gMaxnmus];
	float fTmuIso03emEt     [gMaxnmus];
	float fTmuIso03hadEt    [gMaxnmus];
	float fTmuIso03emVetoEt [gMaxnmus];
	float fTmuIso03hadVetoEt[gMaxnmus];
	float fTmuIso05sumPt    [gMaxnmus];
	float fTmuIso05emEt     [gMaxnmus];
	float fTmuIso05hadEt    [gMaxnmus];
	float fTmueecal         [gMaxnmus];
	float fTmuehcal         [gMaxnmus];

// - Impact Parameters
	float fTmud0bs[gMaxnmus];
	float fTmud0pv[gMaxnmus];
	float fTmud0E [gMaxnmus];
	float fTmudzbs[gMaxnmus];
	float fTmudzpv[gMaxnmus];
	float fTmudzE [gMaxnmus];

// - MuID Variables
	float fTmunchi2     [gMaxnmus];
	int fTmunglhits     [gMaxnmus];
	int fTmunmuhits     [gMaxnmus];
	int fTmuntkhits     [gMaxnmus];
	int fTmunpxhits     [gMaxnmus];
	float fTmuinntknchi2[gMaxnmus];
	int fTmunmatches    [gMaxnmus];
	int fTmunchambers   [gMaxnmus];
	float fTmucalocomp  [gMaxnmus];
	float fTmusegmcomp  [gMaxnmus];

	int fTmuIsGMPT                 [gMaxnmus];
	int fTmuIsGMTkChiComp          [gMaxnmus];
	int fTmuIsGMStaChiComp         [gMaxnmus];
	int fTmuIsGMTkKinkTight        [gMaxnmus];
	int fTmuIsAllStaMuons          [gMaxnmus];
	int fTmuIsAllTrkMuons          [gMaxnmus];
	int fTmuIsTrkMuArb             [gMaxnmus];
	int fTmuIsAllArb               [gMaxnmus];
	int fTmuIsTMLastStationLoose   [gMaxnmus];
	int fTmuIsTMLastStationTight   [gMaxnmus];
	int fTmuIsTM2DCompLoose        [gMaxnmus];
	int fTmuIsTM2DCompTight        [gMaxnmus];
	int fTmuIsTMOneStationLoose    [gMaxnmus];
	int fTmuIsTMOneStationTight    [gMaxnmus];
	int fTmuIsTMLSOPL              [gMaxnmus];
	int fTmuIsTMLastStationAngLoose[gMaxnmus];
	int fTmuIsTMLastStationAngTight[gMaxnmus];
	int fTmuIsTMOneStationAngTight [gMaxnmus];
	int fTmuIsTMOneStationAngLoose [gMaxnmus];

// - Gen Info:
	int   fTGenMuId      [gMaxnmus];
	int   fTGenMuStatus  [gMaxnmus];
	float fTGenMuPt      [gMaxnmus];
	float fTGenMuEta     [gMaxnmus];
	float fTGenMuPhi     [gMaxnmus];
	float fTGenMuE       [gMaxnmus];
	int   fTGenMuMId     [gMaxnmus];
	int   fTGenMuMStatus [gMaxnmus];
	float fTGenMuMPt     [gMaxnmus];
	float fTGenMuMEta    [gMaxnmus];
	float fTGenMuMPhi    [gMaxnmus];
	float fTGenMuME      [gMaxnmus];
	int   fTGenMuGMId    [gMaxnmus];
	int   fTGenMuGMStatus[gMaxnmus];
	float fTGenMuGMPt    [gMaxnmus];
	float fTGenMuGMEta   [gMaxnmus];
	float fTGenMuGMPhi   [gMaxnmus];
	float fTGenMuGME     [gMaxnmus];


// Electrons:
	int fTneles;
	int fTnelestot; // before preselection
	int fTgoodel[gMaxneles];
	int fTeIsIso[gMaxneles];
	int fTeChargeMisIDProb[gMaxneles];
// Kinematics
	float fTepx[gMaxneles];
	float fTepy[gMaxneles];
	float fTepz[gMaxneles];
	float fTept[gMaxneles];
	float fTeptE[gMaxneles];
	float fTee[gMaxneles];
	float fTeet[gMaxneles];
	float fTeeta[gMaxneles];
	float fTephi[gMaxneles];
	float fTetheta[gMaxneles];
	float fTesceta[gMaxneles];
	float fTegsfpt[gMaxneles];
	float fTegsfeta[gMaxneles];
	float fTegsfphi[gMaxneles];
	float fTetrkmomerror[gMaxneles];
	float fTeecalergerror[gMaxneles];
	float fTeelemomerror[gMaxneles];
	int   fTenbrems[gMaxneles];

// Impact parameter
	float fTed0bs[gMaxneles];
	float fTed0pv[gMaxneles];
	float fTed0E[gMaxneles];
	float fTedzbs[gMaxneles];
	float fTedzpv[gMaxneles];
	float fTedzE[gMaxneles];
// Isolation
	float fTeiso03[gMaxneles];
	float fTeiso04[gMaxneles];
	float fTdr03tksumpt[gMaxneles];
	float fTdr04tksumpt[gMaxneles];
	float fTdr03ecalrechitsumet[gMaxneles];
	float fTdr04ecalrechitsumet[gMaxneles];
	float fTdr03hcaltowersumet[gMaxneles];
	float fTdr04hcaltowersumet[gMaxneles];
	float fTenchi2[gMaxneles];
// Electron ID
	float fTeIDMva[gMaxneles];
	int fTeIDTight[gMaxneles];
	int fTeIDLoose[gMaxneles];
	int fTeIDRobustTight[gMaxneles];
	int fTeIDRobustLoose[gMaxneles];
	int fTeIDsimpleWPrelIso[gMaxneles];
	int fTeIDsimpleWP95relIso[gMaxneles];
	int fTeIDsimpleWP90relIso[gMaxneles];
	int fTeIDsimpleWP85relIso[gMaxneles];
	int fTeIDsimpleWP80relIso[gMaxneles];
	int fTecharge[gMaxneles];
	int fTeCInfoIsGsfCtfCons[gMaxneles];
	int fTeCInfoIsGsfCtfScPixCons[gMaxneles];
	int fTeCInfoIsGsfScPixCons[gMaxneles];
	int fTeCInfoScPixCharge[gMaxneles];
	float fTeClosestCtfTrackpt[gMaxneles];
	float fTeClosestCtfTracketa[gMaxneles];
	float fTeClosestCtfTrackphi[gMaxneles];
	int fTeClosestCtfTrackcharge[gMaxneles];
	int fTeInGap[gMaxneles];  // seed crystal next to a gap
	int fTeEcalDriven[gMaxneles];
	int fTeTrackerDriven[gMaxneles];
	int fTeBasicClustersSize[gMaxneles];
	float fTefbrem[gMaxneles];
	float fTeHcalOverEcal[gMaxneles];
	float fTeE1x5[gMaxneles];                      // 5x5 arround seed
	float fTeE5x5[gMaxneles];                      // 5x5 arround seed
	float fTeE2x5Max[gMaxneles];                   // 2x5 arround seed
	float fTeSigmaIetaIeta[gMaxneles];             // shower shape covariance
	float fTeDeltaPhiSeedClusterAtCalo[gMaxneles]; // Dphi (seed-track) at calo from p_out
	float fTeDeltaEtaSeedClusterAtCalo[gMaxneles]; // outermost track state extrapolated at calo
	float fTeDeltaPhiSuperClusterAtVtx[gMaxneles]; // Dphi (sc-track) at calo extrapolated from p_in
	float fTeDeltaEtaSuperClusterAtVtx[gMaxneles]; // Deta (sc-track) at calo extrapolated from p_in
	float fTecaloenergy[gMaxneles];                // caloEnergy() = supercluster energy 99.9% of the time
	float fTetrkmomatvtx[gMaxneles];               // trackMomentumAtVtx().R()
	float fTeESuperClusterOverP[gMaxneles];        // Esc/Pin
	int fTeNumberOfMissingInnerHits[gMaxneles];
	// int fTeIsInJet[gMaxneles];
	// float fTeSharedPx[gMaxneles];
	// float fTeSharedPy[gMaxneles];
	// float fTeSharedPz[gMaxneles];
	// float fTeSharedEnergy[gMaxneles];
	// int fTeDupEl[gMaxneles];
	float fTeConvPartTrackDist[gMaxneles];
	float fTeConvPartTrackDCot[gMaxneles];
	float fTeConvPartTrackPt[gMaxneles];
	float fTeConvPartTrackEta[gMaxneles];
	float fTeConvPartTrackPhi[gMaxneles];
	float fTeConvPartTrackCharge[gMaxneles];
// Spike cleaning
	int fTeScSeedSeverity[gMaxneles];
	float fTeE1OverE9[gMaxneles];
	float fTeS4OverS1[gMaxneles];


// - Gen Info:
	int   fTGenElId[gMaxneles];
	int   fTGenElStatus[gMaxneles];
	float fTGenElPt[gMaxneles];
	float fTGenElEta[gMaxneles];
	float fTGenElPhi[gMaxneles];
	float fTGenElE[gMaxneles];
	int   fTGenElMId[gMaxneles];
	int   fTGenElMStatus[gMaxneles];
	float fTGenElMPt[gMaxneles];
	float fTGenElMEta[gMaxneles];
	float fTGenElMPhi[gMaxneles];
	float fTGenElME[gMaxneles];
	int   fTGenElGMId[gMaxneles];
	int   fTGenElGMStatus[gMaxneles];
	float fTGenElGMPt[gMaxneles];
	float fTGenElGMEta[gMaxneles];
	float fTGenElGMPhi[gMaxneles];
	float fTGenElGME[gMaxneles];

// Photons:
	int fTnphotons;
	int fTnphotonstot; // before preselection
	int fTgoodphoton[gMaxnphos];
	int fTPhotIsIso[gMaxnphos];
	float fTPhotPt[gMaxnphos];  
	float fTPhotPx[gMaxnphos];  
	float fTPhotPy[gMaxnphos];  
	float fTPhotPz[gMaxnphos];  
	float fTPhotEta[gMaxnphos];
	float fTPhotPhi[gMaxnphos];
	float fTPhotEnergy[gMaxnphos];
	float fTPhotIso03Ecal[gMaxnphos];
	float fTPhotIso03Hcal[gMaxnphos];
	float fTPhotIso03TrkSolid[gMaxnphos];
	float fTPhotIso03TrkHollow[gMaxnphos];
	float fTPhotIso03[gMaxnphos];
	float fTPhotIso04Ecal[gMaxnphos];
	float fTPhotIso04Hcal[gMaxnphos];
	float fTPhotIso04TrkSolid[gMaxnphos];
	float fTPhotIso04TrkHollow[gMaxnphos];
	float fTPhotIso04[gMaxnphos];
	float fTPhotR9[gMaxnphos];
	float fTPhotSCEnergy[gMaxnphos];
	float fTPhotcaloPosX[gMaxnphos];
	float fTPhotcaloPosY[gMaxnphos];
	float fTPhotcaloPosZ[gMaxnphos];
	float fTPhotHoverE[gMaxnphos];
	float fTPhotH1overE[gMaxnphos];
	float fTPhotH2overE[gMaxnphos];
	float fTPhotSigmaIetaIeta[gMaxnphos];
	float fTPhotSCEtaWidth[gMaxnphos];
	float fTPhotSCSigmaPhiPhi[gMaxnphos];
	int   fTPhotHasPixSeed[gMaxnphos];
	int   fTPhotHasConvTrks[gMaxnphos];
	// int   fTPhotIsInJet[gMaxnphos];
	// int   fTPhotDupEl[gMaxnphos];
	// float fTPhotSharedPx[gMaxnphos];
	// float fTPhotSharedPy[gMaxnphos];
	// float fTPhotSharedPz[gMaxnphos];
	// float fTPhotSharedEnergy[gMaxnphos];
// Spike cleaning
	int   fTPhotScSeedSeverity[gMaxnphos];
	float fTPhotE1OverE9[gMaxnphos];
	float fTPhotS4OverS1[gMaxnphos];


// Jets:
	int fTnjets;
	int fTnjetstot; // before preselection
	int fTgoodjet[gMaxnjets];
	float fTjpx[gMaxnjets];
	float fTjpy[gMaxnjets];
	float fTjpz[gMaxnjets];
	float fTje[gMaxnjets];
	float fTjet[gMaxnjets];
	float fTjpt[gMaxnjets];
	float fTjeta[gMaxnjets];
	float fTjphi[gMaxnjets];
	float fTjEcorr[gMaxnjets];
	float fTjArea[gMaxnjets];

	float fTJEtaRms[gMaxnjets];
	float fTJPhiRms[gMaxnjets];

	int fTjNconstituents[gMaxnjets];
	int fTjChMult[gMaxnjets];
	int fTjNeuMult[gMaxnjets];
	float fTjChHadFrac[gMaxnjets];
	float fTjNeuHadFrac[gMaxnjets];
	float fTjChEmFrac[gMaxnjets];
	float fTjNeuEmFrac[gMaxnjets];
	float fTjChMuEFrac[gMaxnjets];

	float fTjbTagProbTkCntHighEff[gMaxnjets];
	float fTjbTagProbTkCntHighPur[gMaxnjets];
	float fTjbTagProbSimpSVHighEff[gMaxnjets];
	float fTjbTagProbSimpSVHighPur[gMaxnjets];
	float fTjMass[gMaxnjets];
	float fTjtrk1px[gMaxnjets];
	float fTjtrk1py[gMaxnjets];
	float fTjtrk1pz[gMaxnjets];
	float fTjtrk2px[gMaxnjets];
	float fTjtrk2py[gMaxnjets];
	float fTjtrk2pz[gMaxnjets];
	float fTjtrk3px[gMaxnjets];
	float fTjtrk3py[gMaxnjets];
	float fTjtrk3pz[gMaxnjets];
	float fTjeMinDR[gMaxnjets];
	float fTjetVtxx[gMaxnjets];
	float fTjetVtxy[gMaxnjets];
	float fTjetVtxz[gMaxnjets];
	float fTjetVtxExx[gMaxnjets];
	float fTjetVtxEyx[gMaxnjets];
	float fTjetVtxEyy[gMaxnjets];
	float fTjetVtxEzy[gMaxnjets];
	float fTjetVtxEzz[gMaxnjets];
	float fTjetVtxEzx[gMaxnjets];
	float fTjetVtxNChi2[gMaxnjets];
	
	int   fTjetGenJetIndex[gMaxnjets];

// Tracks:
	int fTntracks;
	int fTntrackstot; // before preselection
	int fTgoodtrk[gMaxntrks];
	float fTtrkpt[gMaxntrks]; // this is actually charge*pt
	float fTtrketa[gMaxntrks];
	float fTtrkphi[gMaxntrks];
	float fTtrknchi2[gMaxntrks];
	float fTtrknhits[gMaxntrks];

// (M)E(T):
	float fTTrkPtSumx;
	float fTTrkPtSumy;
	float fTTrkPtSumphi;
	float fTTrkPtSum;
	float fTSumEt;
	float fTECALSumEt;
	float fTHCALSumEt;
	float fTECALEsumx;
	float fTECALEsumy;
	float fTECALEsumz;
	float fTECALMETphi;
	float fTECALMETeta;
	float fTECALMET;
	float fTHCALEsumx;
	float fTHCALEsumy;
	float fTHCALEsumz;
	float fTHCALMETphi;
	float fTHCALMETeta;
	float fTHCALMET;
	float fTRawMET;
	float fTRawMETpx;
	float fTRawMETpy;
	float fTRawMETphi;
	float fTRawMETemEtFrac;
	float fTRawMETemEtInEB;
	float fTRawMETemEtInEE;
	float fTRawMETemEtInHF;
	float fTRawMEThadEtFrac;
	float fTRawMEThadEtInHB;
	float fTRawMEThadEtInHE;
	float fTRawMEThadEtInHF;
	float fTRawMETSignificance;
	float fTGenMET;
	float fTGenMETpx;
	float fTGenMETpy;
	float fTGenMETphi;
	float fTTCMET;
	float fTTCMETpx;
	float fTTCMETpy;
	float fTTCMETphi;
	float fTTCMETSignificance;
	float fTMuJESCorrMET;
	float fTMuJESCorrMETpx;
	float fTMuJESCorrMETpy;
	float fTMuJESCorrMETphi;
	float fTPFMET;
	float fTPFMETpx;
	float fTPFMETpy;
	float fTPFMETphi;
	float fTPFMETSignificance;
	float fTPFMETPAT;
	float fTPFMETPATpx;
	float fTPFMETPATpy;
	float fTPFMETPATphi;
	float fTPFMETPATSignificance;
	float fTMETR12;
	float fTMETR21;
////////////////////////////////////////////////////////
};
