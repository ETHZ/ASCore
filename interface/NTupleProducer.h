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
// $Id: NTupleProducer.h,v 1.125 2012/02/27 14:09:41 peruzzi Exp $
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

// PF stuff
#include "DataFormats/ParticleFlowReco/interface/PFBlock.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlockElement.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "RecoParticleFlow/PFClusterTools/interface/ClusterClusterMapping.h"
#include "CommonTools/ParticleFlow/plugins/PFPileUp.h"
#include "DataFormats/ParticleFlowCandidate/interface/PileUpPFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PileUpPFCandidateFwd.h"
#include "RecoParticleFlow/PFClusterTools/interface/ClusterClusterMapping.h"

#include "RecoEcal/EgammaCoreTools/interface/EcalClusterFunctionBaseClass.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterFunctionFactory.h"

// Helpers
#include "Math/VectorUtil.h"
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"

// Local classes
#include "DiLeptonAnalysis/NTupleProducer/interface/JetFillerReco.h"
#include "DiLeptonAnalysis/NTupleProducer/interface/JetFillerPat.h"
#include "DiLeptonAnalysis/NTupleProducer/interface/LeptonFillerPat.h"

#include "h2gglobe/VertexAnalysis/interface/HggVertexAnalyzer.h"
#include "h2gglobe/VertexAnalysis/interface/HggVertexFromConversions.h"
#include "h2gglobe/VertexAnalysis/interface/PhotonInfo.h"
#include "h2gglobe/VertexAnalysis/interface/VertexAlgoParameters.h"
#include "TMVA/Reader.h"

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

  double DeltaR(double phi1, double phi2, double eta1, double eta2);
  void FillPhotonIsoVariables(double photonEta, double photonPhi, double photonVz, int type, bool isPU, edm::Handle<reco::PFCandidateCollection>& pfCandidates, int ipf, int phoqi);
//   void FillPhotonIsoVariables_Frixione_Neutrals(int type, int ipf, int phoqi);
//   void FillPhotonIsoVariables_Frixione_ChHad(int type, bool isPU, int ipf, int phoqi);
  reco::VertexRef chargedHadronVertex( const edm::Handle<reco::VertexCollection>& vertices, const reco::PFCandidate& pfcand ) const ;
  int FindPFCandType(int id);
  bool isInPhiCracks(double phi, double eta);
  bool isInEtaCracks(double eta);
  bool CheckPhotonPFCandOverlap(reco::SuperClusterRef scRef, edm::Handle<reco::PFCandidateCollection>& pfCandidates, int i);
  double DeltaPhi(double phi1, double phi2);
  double phiNorm(float &phi);
  double etaTransformation(float EtaParticle , float Zvertex);

  EcalClusterFunctionBaseClass *CrackCorrFunc;
  EcalClusterFunctionBaseClass *LocalCorrFunc;

  std::string perVtxMvaWeights, perVtxMvaMethod;
  std::string perEvtMvaWeights, perEvtMvaMethod;
  VertexAlgoParameters vtxAlgoParams;
  std::vector<std::string> rankVariables;

  PhotonInfo fillPhotonInfos(int p1, bool useAllConvs);
  std::vector<int> HggVertexSelection(HggVertexAnalyzer & vtxAna, HggVertexFromConversions & vtxAnaFromConv, 
				      PhotonInfo & pho1, PhotonInfo & pho2, std::vector<std::string> & vtxVarNames, 
				      bool useMva, TMVA::Reader * tmvaReader, std::string tmvaMethod);
  int  matchPhotonToConversion(int lpho);
  bool tkIsHighPurity(reco::TrackRef tk) const;
  bool TrackCut(reco::TrackRef tk) const;
  bool ConversionsCut(const reco::Conversion &conv);

// ----------member data ---------------------------
	edm::Service<TFileService> fTFileService;
        AdaptiveVertexFitter avFitter;

  //for OOT reweighting in Summer11_S3 samples
  edm::LumiReWeighting LumiWeights_;


	std::vector<JetFillerBase*>     jetFillers;
        std::vector<PatMuonFiller*>     muonFillers;
        std::vector<PatElectronFiller*> electronFillers;
        std::vector<PatTauFiller*>      tauFillers;

	bool fIsRealData;
	bool fIsModelScan;
	int fNTotEvents;
	int fNFillTree;

        bool doVertexingFlag;

	static const int gMaxnmus     = 30;
	static const int gMaxneles    = 20;
	static const int gMaxntaus    = 20;
	static const int gMaxnjets    = 100;
	static const int gMaxntrks    = 800;
	static const int gMaxnphos    = 50;
        static const int gMaxnconv    = 50;
        static const int gMaxnSC      = 100;
	static const int gMaxngenlept = 100;
	static const int gMaxngenphot = 100;
	static const int gMaxngenjets = 100;
	static const int gMaxnvrtx    = 25;
	static const int gMaxnpileup  = 50;
	static const int gMaxnEBhits  = 20;
        static const int gMaxngenvtx = 60;
        static const int nStoredGenParticles = 2000;

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
	edm::InputTag fHBHENoiseResultTagIso;
	edm::InputTag fSrcRho;
	edm::InputTag fSrcRhoPFnoPU;
	edm::InputTag fpdfWeightTag;
        edm::InputTag pfphotonsProducerTag;
        edm::InputTag pfProducerTag;
        edm::InputTag fSCTagBarrel;
        edm::InputTag fSCTagEndcap;
        edm::InputTag fTrackCollForVertexing;
        edm::InputTag fallConversionsCollForVertexing;

	int NPdfs;
	float fTpdfW[100];
	float fTpdfWsum;
	int fTprocess;

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
        float fMinSCraw;

	float fMingenleptpt; 
	float fMaxgenlepteta;
	float fMingenphotpt; 
	float fMaxgenphoteta;
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
        float fRTMinSCraw;


	int fRTmaxnmu;
	int fRTmaxnel;
	int fRTmaxnjet;
	int fRTmaxntrk;
	int fRTmaxnphot;

	float fMinebrechitE; 

	TTree *fEventTree;

	// General event information
	int fTrunnumber;
	unsigned int fTeventnumber;
	int fTlumisection;
	float fTpthat;
	float fTqcdPartonicHT;
	int fTsigprocid;
	float fTpdfscalePDF;
	int fTpdfid1;
	int fTpdfid2;
	float fTpdfx1;
	float fTpdfx2;
	float fTpdfxPDF1;
	float fTpdfxPDF2;
	float fTgenweight;
	float fTextxslo;
	float fTintxs;
	float fTweight;
        float fTMassGlu;
        float fTMassChi;
        float fTMassLSP;
        float fTxSMS;
        float fTxbarSMS;
        float fTSUSYScanM0;
        float fTSUSYScanM12;
        float fTSUSYScanA0;
        float fTSUSYScanMu;
        float fTSUSYScanCrossSection;
        float fTSUSYScanTanBeta;


// Generator information
	int fTnGenParticles;
        int fTgenInfoId[nStoredGenParticles];
        int fTgenInfoStatus[nStoredGenParticles];
        float fTgenInfoMass[nStoredGenParticles];
        int fTgenInfoNMo[nStoredGenParticles];
	float fTgenInfoMo1Pt[nStoredGenParticles];
	float fTgenInfoMo2Pt[nStoredGenParticles];
        int fTgenInfoNDa[nStoredGenParticles];
        int fTgenInfoMo1[nStoredGenParticles];
        int fTgenInfoMo2[nStoredGenParticles];
        int fTgenInfoDa1[nStoredGenParticles];
        int fTgenInfoDa2[nStoredGenParticles];
        float fTgenInfoPt[nStoredGenParticles];
        float fTgenInfoEta[nStoredGenParticles];
        float fTgenInfoPhi[nStoredGenParticles];
        float fTgenInfoPx[nStoredGenParticles];
        float fTgenInfoPy[nStoredGenParticles];
        float fTgenInfoPz[nStoredGenParticles];
        float fTgenInfoM[nStoredGenParticles];
	float fTgenInfoPromptFlag[nStoredGenParticles];
	int fTgenInfoMoIndex[nStoredGenParticles];
	int fTPromptnessLevel[nStoredGenParticles];


	// Pile-up
	int fTpuNumInteractions;
	int fTpuNumTrueInteractions;
        int fTpuNumFilled;
        int fTpuOOTNumInteractionsLate;
        int fTpuOOTNumInteractionsEarly;

	float fTpuZpositions[gMaxnpileup];
	float fTpuSumpT_lowpT[gMaxnpileup];
	float fTpuSumpT_highpT[gMaxnpileup];
	float fTpuNtrks_lowpT[gMaxnpileup];
	float fTpuNtrks_highpT[gMaxnpileup];
	float fTrho; // rho from L1FastJetCorrection
	float fTrhoPFnoPU; // rho from L1FastJetCorrection running PFnoPU
	// float fTpuInstLumi[gMaxnpileup];
        TH1I* fHpileupstat;
        TH1I* fHtruepileupstat;
  	float fTpuWeightTotal;
  	float fTpuWeightInTime;
	std::vector<std::string> fTPileUpHistoData;
	std::vector<std::string> fTPileUpHistoMC;

	// ECAL & HCAL Noise
	int fTHBHENoiseFlag;
	int fTHBHENoiseFlagIso;
	int fRecovRecHitFilterFlag;
	int fTra2TrackingFailureFilterFlag;
	//FR int fPBNRFlag;

	int fTnEBhits;
	float fTEBrechitE[gMaxnEBhits];
	float fTEBrechitPt[gMaxnEBhits];
	float fTEBrechitEta[gMaxnEBhits];
	float fTEBrechitPhi[gMaxnEBhits];
	float fTEBrechitChi2[gMaxnEBhits];
	float fTEBrechitTime[gMaxnEBhits];
	float fTEBrechitE4oE1[gMaxnEBhits];
	float fTEBrechitE2oE9[gMaxnEBhits];
	
	int fTecalDeadTPFilterFlag;
	// int fTEcalDeadCellBEFlag;
	// static const unsigned int gMaxnECALGapClusters = 50;
	// unsigned int fTnECALGapClusters;
	// float fTEcalGapBE[gMaxnECALGapClusters];
	// int fTEcalGapClusterSize[gMaxnECALGapClusters];
	//
	
	// CSCBeamHalo
	int fTcscTightHaloID;

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
	static const unsigned int gMaxhltbits = 500;
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
	int fTflagmaxgenphotexc; // Found more than 100 genphotons in event
	int fTflagmaxgenjetexc;  // Found more than 100 genjets in event
	int fTflagmaxvrtxexc;    // Found more than 25 vertices in event
	int fTflagmaxgenpartexc; // Found more than nStoredGenParticles in event

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

  // GenPhotons
  int fTngenphotons;
  float fTGenPhotonPt[gMaxngenphot];      
  float fTGenPhotonEta[gMaxngenphot];     
  float fTGenPhotonPhi[gMaxngenphot];    
  float fTGenPhotonPartonMindR[gMaxngenphot];    
  int fTGenPhotonMotherID[gMaxngenphot];
  int fTGenPhotonMotherStatus[gMaxngenphot];

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
        int fTElSCindex[gMaxneles];
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
  float fTPhotSigmaEtaEta[gMaxnphos];
  float fTPhote1x5[gMaxnphos];
  float fTPhote2x5[gMaxnphos];
  float fTPhote3x3[gMaxnphos];
  float fTPhote5x5[gMaxnphos];
  float fTPhotmaxEnergyXtal[gMaxnphos];
  float fTPhotIso03HcalDepth1[gMaxnphos];
  float fTPhotIso03HcalDepth2[gMaxnphos];
  float fTPhotIso04HcalDepth1[gMaxnphos];
  float fTPhotIso04HcalDepth2[gMaxnphos];
  int fTPhotIso03nTrksSolid[gMaxnphos];
  int fTPhotIso03nTrksHollow[gMaxnphos];
  int fTPhotIso04nTrksSolid[gMaxnphos];
  int fTPhotIso04nTrksHollow[gMaxnphos];
  int fTPhotisEB[gMaxnphos];
  int fTPhotisEE[gMaxnphos];
  int fTPhotisEBEtaGap[gMaxnphos];
  int fTPhotisEBPhiGap[gMaxnphos];
  int fTPhotisEERingGap[gMaxnphos];
  int fTPhotisEEDeeGap[gMaxnphos];
  int fTPhotisEBEEGap[gMaxnphos];
  int fTPhotisPFlowPhoton[gMaxnphos];
  int fTPhotisStandardPhoton[gMaxnphos];
  int fTPhotMCmatchindex[gMaxnphos];
  int fTPhotMCmatchexitcode[gMaxnphos];
 
  float fT_pho_ChargedHadronIso[gMaxnphos];
  float fT_pho_NeutralHadronIso[gMaxnphos];
  float fT_pho_PhotonIso[gMaxnphos];
  int fT_pho_isPFPhoton[gMaxnphos];
  int fT_pho_isPFElectron[gMaxnphos];
  int fTPhotSCindex[gMaxnphos];

float fT_pho_Cone04PhotonIso_dR0_dEta0_pt0[gMaxnphos];
float fT_pho_Cone04PhotonIso_dR0_dEta0_pt5[gMaxnphos];
float fT_pho_Cone04PhotonIso_dR8_dEta0_pt0[gMaxnphos];
float fT_pho_Cone04PhotonIso_dR8_dEta0_pt5[gMaxnphos];
float fT_pho_Cone01PhotonIso_dR045EB070EE_dEta015_pt08EB1EE_mvVtx[gMaxnphos];
float fT_pho_Cone02PhotonIso_dR045EB070EE_dEta015_pt08EB1EE_mvVtx[gMaxnphos];
float fT_pho_Cone03PhotonIso_dR045EB070EE_dEta015_pt08EB1EE_mvVtx[gMaxnphos];
float fT_pho_Cone04PhotonIso_dR045EB070EE_dEta015_pt08EB1EE_mvVtx[gMaxnphos];
float fT_pho_Cone04NeutralHadronIso_dR0_dEta0_pt0[gMaxnphos];
float fT_pho_Cone04NeutralHadronIso_dR0_dEta0_pt5[gMaxnphos];
float fT_pho_Cone04NeutralHadronIso_dR0_dEta0_pt0_nocracks[gMaxnphos];
float fT_pho_Cone04NeutralHadronIso_dR0_dEta0_pt5_nocracks[gMaxnphos];
float fT_pho_Cone04NeutralHadronIso_dR7_dEta0_pt0[gMaxnphos];
float fT_pho_Cone04NeutralHadronIso_dR7_dEta0_pt5[gMaxnphos];
float fT_pho_Cone01NeutralHadronIso_dR0_dEta0_pt0_mvVtx[gMaxnphos];
float fT_pho_Cone02NeutralHadronIso_dR0_dEta0_pt0_mvVtx[gMaxnphos];
float fT_pho_Cone03NeutralHadronIso_dR0_dEta0_pt0_mvVtx[gMaxnphos];
float fT_pho_Cone04NeutralHadronIso_dR0_dEta0_pt0_mvVtx[gMaxnphos];
float fT_pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_dz0_old[gMaxnphos];
float fT_pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_PFnoPU_old[gMaxnphos];
float fT_pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_dz0_old[gMaxnphos];
float fT_pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_PFnoPU_old[gMaxnphos];
float fT_pho_Cone01ChargedHadronIso_dR0_dEta0_pt0_dz0[gMaxnphos];
float fT_pho_Cone01ChargedHadronIso_dR0_dEta0_pt0_dz1_dxy01[gMaxnphos];
float fT_pho_Cone01ChargedHadronIso_dR0_dEta0_pt0_PFnoPU[gMaxnphos];
float fT_pho_Cone01ChargedHadronIso_dR015_dEta0_pt0_dz0[gMaxnphos];
float fT_pho_Cone01ChargedHadronIso_dR015_dEta0_pt0_dz1_dxy01[gMaxnphos];
float fT_pho_Cone01ChargedHadronIso_dR015_dEta0_pt0_PFnoPU[gMaxnphos];
float fT_pho_Cone02ChargedHadronIso_dR0_dEta0_pt0_dz0[gMaxnphos];
float fT_pho_Cone02ChargedHadronIso_dR0_dEta0_pt0_dz1_dxy01[gMaxnphos];
float fT_pho_Cone02ChargedHadronIso_dR0_dEta0_pt0_PFnoPU[gMaxnphos];
float fT_pho_Cone02ChargedHadronIso_dR015_dEta0_pt0_dz0[gMaxnphos];
float fT_pho_Cone02ChargedHadronIso_dR015_dEta0_pt0_dz1_dxy01[gMaxnphos];
float fT_pho_Cone02ChargedHadronIso_dR015_dEta0_pt0_PFnoPU[gMaxnphos];
float fT_pho_Cone03ChargedHadronIso_dR0_dEta0_pt0_dz0[gMaxnphos];
float fT_pho_Cone03ChargedHadronIso_dR0_dEta0_pt0_dz1_dxy01[gMaxnphos];
float fT_pho_Cone03ChargedHadronIso_dR0_dEta0_pt0_PFnoPU[gMaxnphos];
float fT_pho_Cone03ChargedHadronIso_dR015_dEta0_pt0_dz0[gMaxnphos];
float fT_pho_Cone03ChargedHadronIso_dR015_dEta0_pt0_dz1_dxy01[gMaxnphos];
float fT_pho_Cone03ChargedHadronIso_dR015_dEta0_pt0_PFnoPU[gMaxnphos];
float fT_pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_dz0[gMaxnphos];
float fT_pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_dz1_dxy01[gMaxnphos];
float fT_pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_PFnoPU[gMaxnphos];
float fT_pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_dz0[gMaxnphos];
float fT_pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_dz1_dxy01[gMaxnphos];
float fT_pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_PFnoPU[gMaxnphos];

  std::vector<TVector3> pho_conv_vtx;
  std::vector<TVector3> pho_conv_refitted_momentum;
  bool pho_conv_validvtx[gMaxnphos];
  int pho_conv_ntracks[gMaxnphos];
  float pho_conv_chi2_probability[gMaxnphos];
  float pho_conv_eoverp[gMaxnphos];
  
  int conv_n;
  std::vector<TVector3> conv_vtx;
  std::vector<TVector3> conv_refitted_momentum;
  bool conv_validvtx[gMaxnconv];
  int conv_ntracks[gMaxnconv];
  float conv_chi2_probability[gMaxnconv];
  float conv_eoverp[gMaxnconv];
  float conv_zofprimvtxfromtrks[gMaxnconv];

  Int_t gv_n;
  std::vector<TVector3> gv_pos;
  std::vector<TVector3> gv_p3;
  Float_t gv_sumPtHi[gMaxngenvtx];
  Float_t gv_sumPtLo[gMaxngenvtx];
  Int_t gv_nTkHi[gMaxngenvtx];
  Int_t gv_nTkLo[gMaxngenvtx];

  // SC
  int fTnSC;
  float fTSCraw[gMaxnSC];
  float fTSCpre[gMaxnSC];
  float fTSCenergy[gMaxnSC];
  float fTSCeta[gMaxnSC];
  float fTSCphi[gMaxnSC];
  float fTSCsigmaPhi[gMaxnSC];
  float fTSCsigmaEta[gMaxnSC];
  float fTSCbrem[gMaxnSC];
  float fTSCR9[gMaxnSC];
  float fTSCcrackcorrseed[gMaxnSC];
  float fTSCcrackcorr[gMaxnSC];
  float fTSClocalcorrseed[gMaxnSC];
  float fTSClocalcorr[gMaxnSC];
  float fTSCcrackcorrseedfactor[gMaxnSC];
  float fTSClocalcorrseedfactor[gMaxnSC];


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
	float fTjPhoFrac[gMaxnjets];
	float fTjHFHadFrac[gMaxnjets];
	float fTjHFEMFrac[gMaxnjets];

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
	float fTtrkVtxDz[gMaxntrks];
	float fTtrkVtxDxy[gMaxntrks];

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
        float fTPFSumEt;
	float fTPFMETPAT;
	float fTPFMETPATpx;
	float fTPFMETPATpy;
	float fTPFMETPATphi;
	float fTPFMETPATSignificance;
	float fTMETR12;
	float fTMETR21;
////////////////////////////////////////////////////////
};



class ETHVertexInfo : public VertexInfoAdapter
{
public:

  ETHVertexInfo(int nvtx, float * vtxx, float * vtxy, float * vtxz, 
		int ntracks, float * tkpx, float * tkpy, float * tkpz,
		float * tkPtErr, int * tkVtxId, 
		float * tkd0, float * tkd0Err, float * tkdz, float * tkdzErr,
		bool * tkIsHighPurity, std::vector<unsigned short> * vtx_std_tkind, std::vector<float> * vtx_std_tkweight, int * vtx_std_ntks
		);

  virtual int nvtx() const    { return nvtx_; };
  virtual int ntracks() const { return ntracks_; };

  virtual bool hasVtxTracks()  const { return true; };
  virtual const unsigned short * vtxTracks(int ii) const { return &(vtx_std_tkind_[ii][0]); };
  virtual int vtxNTracks(int ii) const { return vtx_std_ntks_[ii]; };
  virtual const float * vtxTkWeights(int ii) const { return &(vtx_std_tkweight_[ii][0]); };

  virtual float tkpx(int ii) const { return tkpx_ != 0 ? tkpx_[ii] : 0.; };
  virtual float tkpy(int ii) const { return tkpx_ != 0 ? tkpy_[ii] : 0.; };
  virtual float tkpz(int ii) const { return tkpx_ != 0 ? tkpz_[ii] : 0.; };
	
  virtual float tkPtErr(int ii) const { return tkPtErr_  != 0 ? tkPtErr_[ii] : 999.; };
  virtual int   tkVtxId(int ii) const { return tkVtxId_  != 0 ? tkVtxId_[ii] : 999; };

  //	virtual float tkWeight(int ii, int jj) const { return tkWeight_ != 0 ? tkWeight_[ii]*(float)( tkVtxId(ii) == jj) : 0.; };
  virtual float tkWeight(int ii, int jj) const { return vtx_std_tkweight_[jj][ii]; };

	
  virtual float vtxx(int ii) const { return vtxx_ != 0 ? vtxx_[ii] : 0.; };
  virtual float vtxy(int ii) const { return vtxy_ != 0 ? vtxy_[ii] : 0.; };
  virtual float vtxz(int ii) const { return vtxz_ != 0 ? vtxz_[ii] : 0.; };

  virtual float tkd0(int ii, int jj) const { assert(tkVtxId(ii) == jj); return tkd0_ != 0 ? tkd0_[ii] : 0.; };
  virtual float tkd0Err(int ii, int jj) const { assert(tkVtxId(ii) == jj); return tkd0Err_ != 0 ? tkd0Err_[ii] : 0.; };

  virtual float tkdz(int ii, int jj) const { assert(tkVtxId(ii) == jj); return tkdz_ != 0 ? tkdz_[ii] : 0.; };
  virtual float tkdzErr(int ii, int jj) const { assert(tkVtxId(ii) == jj); return tkdzErr_ != 0 ? tkdzErr_[ii] : 0.; };

  virtual bool tkIsHighPurity(int ii) const { return tkIsHighPurity_ != 0 ? tkIsHighPurity_[ii] : 0.; };

  virtual ~ETHVertexInfo();
       	
private:

  int nvtx_;
  float * vtxx_;
  float * vtxy_;
  float * vtxz_;

  int ntracks_;
  float * tkpx_;
  float * tkpy_;
  float * tkpz_;
  float * tkPtErr_;
  int * tkVtxId_;	

  float * tkd0_;
  float * tkd0Err_;
  float * tkdz_;
  float * tkdzErr_;

  bool * tkIsHighPurity_;
  
  std::vector<unsigned short> * vtx_std_tkind_;
  std::vector<float> * vtx_std_tkweight_;
  int * vtx_std_ntks_;

};
