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
// $Id: NTupleProducer.h,v 1.114.2.31 2012/08/15 11:02:36 pandolf Exp $
//
//


// system include files
#include <vector>
#include <string>

// ROOT includes
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TVector3.h"

// Framework include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
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

#include "EGamma/EGammaAnalysisTools/interface/EGammaMvaEleEstimator.h"

#include "Muon/MuonAnalysisTools/interface/MuonMVAEstimator.h"


// Helpers
#include "Math/VectorUtil.h"
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"

// Local classes
#include "DiLeptonAnalysis/NTupleProducer/interface/JetFillerReco.h"
#include "DiLeptonAnalysis/NTupleProducer/interface/JetFillerPat.h"
#include "DiLeptonAnalysis/NTupleProducer/interface/LeptonFillerPat.h"
#include "DiLeptonAnalysis/NTupleProducer/interface/PFFiller.h"

#include "h2gglobe/VertexAnalysis/interface/HggVertexFromConversions.h"
#include "h2gglobe/VertexAnalysis/interface/PhotonInfo.h"
#include "h2gglobe/VertexAnalysis/interface/VertexAlgoParameters.h"
#include "TMVA/Reader.h"

typedef math::XYZTLorentzVector LorentzVector;
using namespace reco;

class NTupleProducer : public edm::EDFilter {
public:
  explicit NTupleProducer(const edm::ParameterSet&);
  ~NTupleProducer() {}
  
  virtual void beginJob(void);
  virtual bool beginRun(edm::Run&, const edm::EventSetup&);
  virtual bool filter(edm::Event&, const edm::EventSetup&);
  virtual void endJob(void);
  virtual bool endRun(edm::Run&, const edm::EventSetup&);
  std::vector<const reco::GenParticle*> matchRecoCand(const reco::RecoCandidate *Cand, const edm::Event& iEvent);
  const int matchJet(const reco::Jet *jet);
  
  void declareProducts(void);
  void resetProducts(void);        // Called for each event
  void resetRunProducts(void);     // Called in beginRun
  void putProducts( edm::Event& ); // Called for each event
  void putRunProducts( edm::Event& ); // Called in endRun
  

private:

  typedef std::pair<int,double> OrderPair;
  struct IndexByPt {
    const bool operator()(const OrderPair& j1, const OrderPair& j2 ) const {
      return j1.second > j2.second;
    }

  };
  struct JetRefCompare :
    public std::binary_function<edm::RefToBase<reco::Jet>, edm::RefToBase<reco::Jet>, bool> {
    inline bool operator () (const edm::RefToBase<reco::Jet> &j1,
                             const edm::RefToBase<reco::Jet> &j2) const
    { return j1.id() < j2.id() || (j1.id() == j2.id() && j1.key() < j2.key()); }
  };

  typedef std::map<edm::RefToBase<reco::Jet>, unsigned int, JetRefCompare> FlavourMap;
  
  void FillPhotonIsoVariables(double photonEta, double photonPhi, double photonVz, int type, bool isPU, edm::Handle<reco::PFCandidateCollection>& pfCandidates, int ipf, int phoqi);
  //   void FillPhotonIsoVariables_Frixione_Neutrals(int type, int ipf, int phoqi);
  //   void FillPhotonIsoVariables_Frixione_ChHad(int type, bool isPU, int ipf, int phoqi);
  reco::VertexRef chargedHadronVertex( const edm::Handle<reco::VertexCollection>& vertices, const reco::PFCandidate& pfcand ) const ;
  int FindPFCandType(int id);
  bool isInPhiCracks(double phi, double eta);
  bool isInEtaCracks(double eta);
  bool CheckPhotonPFCandOverlap(reco::SuperClusterRef scRef, edm::Handle<reco::PFCandidateCollection>& pfCandidates, int i);
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

  double DeltaPhi( double phi1, double phi2);

  // ----------member data ---------------------------
  AdaptiveVertexFitter avFitter;

  //for OOT reweighting in Summer11_S3 samples
  edm::LumiReWeighting LumiWeights_;

  std::vector<JetFillerBase*>     jetFillers;
  std::vector<PatMuonFiller*>     muonFillers;
  std::vector<PatElectronFiller*> electronFillers;
  std::vector<PatTauFiller*>      tauFillers;
  std::vector<PFFiller*>          pfFillers;

  bool fIsRealData;
  bool fIsModelScan;
  bool fIsFastSim;
  int fNTotEvents;
  int fNFillTree;

  bool doVertexingFlag;

  static const int gMaxNMus     = 30;
  static const int gMaxNEles    = 20;
  static const int gMaxNJets    = 100;
  static const int gMaxNTrks    = 800;
  static const int gMaxNPhotons = 50;
  static const int gMaxNConv    = 50;
  static const int gMaxNSC      = 100;
  static const int gMaxNGenLept = 100;
  static const int gMaxNGenPhot = 100;
  static const int gMaxNGenJets = 100;
  static const int gMaxNVrtx    = 45;
  static const int gMaxNPileup  = 60;
  static const int gMaxNEBhits  = 20;
  static const int gMaxNGenVtx = 60;
  static const int gMaxNGenParticles = 2000;

  // Maximum configurable number of tags (b-tagging and PF iso)
  static const unsigned int gMaxNPfIsoTags  = 20;
  static const unsigned int gMaxNBtags      = 10;

  static const unsigned int __TRK_AUX_ARRAYS_DIM__ = 2000;
  static const unsigned int __VTX_AUX_ARRAYS_DIM__ = 100;


  edm::InputTag fMuonTag;
  std::vector<edm::InputTag> fMuonPfIsoTagsCustom;
  edm::InputTag fElectronTag;
  std::vector<edm::InputTag> fElePfIsoTagsCustom;
  std::vector<edm::InputTag> fElePfIsoTagsEvent;
  std::string fEleIdWP;
  edm::InputTag fMuIsoDepTkTag;
  edm::InputTag fMuIsoDepECTag;
  edm::InputTag fMuIsoDepHCTag;
  edm::InputTag fJetTag;
  std::string fJetCorrs;
  std::vector<edm::InputTag> fBtagTags;
  edm::InputTag fPartonMatch;
  edm::InputTag fRawCaloMETTag;
  edm::InputTag fTCMETTag;
  edm::InputTag fPFMETTag;
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
  edm::InputTag fSrcRho;
  edm::InputTag fSrcRhoForIso;
  edm::InputTag fpdfWeightTag;
  edm::InputTag pfphotonsProducerTag;
  edm::InputTag pfProducerTag;
  edm::InputTag fSCTagBarrel;
  edm::InputTag fSCTagEndcap;
  edm::InputTag fTrackCollForVertexing;
  edm::InputTag fAllConversionsCollForVertexing;
    
  PFIsolationEstimator isolator;

  EGammaMvaEleEstimator* electronIDMVANonTrig_;
  EGammaMvaEleEstimator* electronIDMVATrig_;

  MuonMVAEstimator* fMuonIsoMVA;

  // Selection cuts
  float fMinMuPt;
  float fMaxMuEta;
  float fMinElPt;
  float fMaxElEta;
  float fMinCorJPt;
  float fMinRawJPt;
  float fMaxJEta;
  float fMinJEMFrac;

  float fMinTrkPt;
  float fMaxTrkEta;
  float fMaxTrkNChi2;
  int	fMinTrkNHits;

  float fMinPhotonPt;
  float fMaxPhotonEta;
  float fMinSCraw;
  float fMinEBRechitE; 
  
  float fMinGenLeptPt; 
  float fMaxGenLeptEta;
  float fMinGenPhotPt; 
  float fMaxGenPhotEta;
  float fMinGenJetPt;
  float fMaxGenJetEta;
	
  bool fFirstevent;

  // Trigger stuff
  HLTConfigProvider fHltConfig;

  /////////////////////////////////////////////////////////////////////////
  // Products stored in the output file

  ////////////////////////////////////////////////////////
  // Run information:
  std::auto_ptr<float> fRExtXSecLO;
  std::auto_ptr<float> fRExtXSecNLO;
  std::auto_ptr<float> fRIntXSec;

  std::auto_ptr<float> fRMinMuPt;
  std::auto_ptr<float> fRMaxMuEta;
  std::auto_ptr<float> fRMinElPt;
  std::auto_ptr<float> fRMaxElEta;
  std::auto_ptr<float> fRMinJPt;
  std::auto_ptr<float> fRMinRawJPt;
  std::auto_ptr<float> fRMaxJEta;
  std::auto_ptr<float> fRMinJEMFrac;

  std::auto_ptr<float> fRMinTrkPt;
  std::auto_ptr<float> fRMaxTrkEta;
  std::auto_ptr<float> fRMaxTrkNChi2;
  std::auto_ptr<int>   fRMinTrkNHits;

  std::auto_ptr<float> fRMinPhotonPt;
  std::auto_ptr<float> fRMaxPhotonEta;
  std::auto_ptr<float> fRMinSCraw;
  std::auto_ptr<float> fRMinEBRechitE;

  std::auto_ptr<float> fRMinGenLeptPt;
  std::auto_ptr<float> fRMaxGenLeptEta;
  std::auto_ptr<float> fRMinGenPhotPt; 
  std::auto_ptr<float> fRMaxGenPhotEta;
  std::auto_ptr<float> fRMinGenJetPt;
  std::auto_ptr<float> fRMaxGenJetEta;
	
  std::auto_ptr<int>   fRMaxNMus;
  std::auto_ptr<int>   fRMaxNEles;
  std::auto_ptr<int>   fRMaxNJets;
  std::auto_ptr<int>   fRMaxNTrks;
  std::auto_ptr<int>   fRMaxNPhotons;    
  std::auto_ptr<int>   fRMaxNSC;
  std::auto_ptr<int>   fRMaxNGenLept;
  std::auto_ptr<int>   fRMaxNGenPhot;
  std::auto_ptr<int>   fRMaxNGenJets;
  std::auto_ptr<int>   fRMaxNVrtx;
  std::auto_ptr<int>   fRMaxNPileup;
  std::auto_ptr<int>   fRMaxNEBhits; 

  std::auto_ptr<std::vector<std::string> > fRHLTNames;  // Full HLT menu
  std::auto_ptr<std::vector<std::string> > fRL1PhysMenu;

  // These vectors are stored in the run info, but they are set in the constructor
  // Since the "Run::put" operation deletes the pointers, we need to reset them
  // at each run, so we have to duplicate the vectors to keep the original info.
  std::auto_ptr<std::vector<std::string> > fRHLTLabels; // HLT Paths to store the triggering objects of
  std::auto_ptr<std::vector<std::string> > fRPileUpData;
  std::auto_ptr<std::vector<std::string> > fRPileUpMC;
  std::vector<std::string> fHLTLabels;
  std::vector<std::string> fPileUpData;
  std::vector<std::string> fPileUpMC;


  ////////////////////////////////////////////////////////
  // Event information:
  // General event information
  std::auto_ptr<int>   fTRun;
  std::auto_ptr<int>   fTEvent;
  std::auto_ptr<int>   fTLumiSection;

  std::auto_ptr<std::vector<float> >  fTpdfW;
  std::auto_ptr<float>  fTpdfWsum;
  std::auto_ptr<int>    fTNPdfs;
  std::auto_ptr<int>    fTprocess;


  // Generator event information
  std::auto_ptr<float> fTPtHat;
  std::auto_ptr<float> fTQCDPartonicHT;
  std::auto_ptr<int>   fTSigProcID;
  std::auto_ptr<float> fTPDFScalePDF;
  std::auto_ptr<int>   fTPDFID1;
  std::auto_ptr<int>   fTPDFID2;
  std::auto_ptr<float> fTPDFx1;
  std::auto_ptr<float> fTPDFx2;
  std::auto_ptr<float> fTPDFxPDF1;
  std::auto_ptr<float> fTPDFxPDF2;
  std::auto_ptr<float> fTGenWeight;
  std::auto_ptr<float> fTWeight;
  std::auto_ptr<float> fTMassGlu;
  std::auto_ptr<float> fTMassChi;
  std::auto_ptr<float> fTMassLSP;
  std::auto_ptr<float> fTxSMS;
  std::auto_ptr<float> fTxbarSMS;
  std::auto_ptr<float> fTM0;
  std::auto_ptr<float> fTM12;
  std::auto_ptr<float> fTA0;
  std::auto_ptr<float> fTsignMu;
  std::auto_ptr<float> fTtanBeta;
  std::auto_ptr<float> fTCrossSection;

  // Generator information
  std::auto_ptr<int> fTnGenParticles;
  std::auto_ptr<std::vector<int> > fTgenInfoId;
  std::auto_ptr<std::vector<int> > fTgenInfoStatus;
  std::auto_ptr<std::vector<int> > fTgenInfoNMo;
  std::auto_ptr<std::vector<int> > fTgenInfoMo1;
  std::auto_ptr<std::vector<int> > fTgenInfoMo2;
  std::auto_ptr<std::vector<int> > fTPromptnessLevel;
  std::auto_ptr<std::vector<float> > fTgenInfoPt;
  std::auto_ptr<std::vector<float> > fTgenInfoEta;
  std::auto_ptr<std::vector<float> > fTgenInfoPhi;
  std::auto_ptr<std::vector<float> > fTgenInfoM;
  std::auto_ptr<std::vector<float> > fTgenInfoPromptFlag;


  // Pile-up event info
  std::auto_ptr<int>  fTPUnumInteractions;
  std::auto_ptr<int>  fTPUnumTrueInteractions;
  std::auto_ptr<int>  fTPUnumFilled;
  std::auto_ptr<int>  fTPUOOTnumInteractionsEarly;
  std::auto_ptr<int>  fTPUOOTnumInteractionsLate;

  std::auto_ptr<std::vector<float> >  fTPUzPositions;
  std::auto_ptr<std::vector<float> >  fTPUsumPtLowPt;
  std::auto_ptr<std::vector<float> >  fTPUsumPtHighPt;
  std::auto_ptr<std::vector<float> >  fTPUnTrksLowPt;
  std::auto_ptr<std::vector<float> >  fTPUnTrksHighPt;
  std::auto_ptr<float>  fTRho; // rho from L1FastJetCorrection
  std::auto_ptr<float>  fTRhoForIso; // rho computed up to eta=2.5
 
  std::auto_ptr<float>  fTPUWeightTotal;
  std::auto_ptr<float>  fTPUWeightInTime;

  //FR std::auto_ptr<int> fPBNRFlag;
  std::auto_ptr<float> fTPFType1MET;
  std::auto_ptr<float> fTPFType1METpx;
  std::auto_ptr<float> fTPFType1METpy;
  std::auto_ptr<float> fTPFType1METphi;
  std::auto_ptr<float> fTPFType1METSignificance;
  std::auto_ptr<float> fTPFType1SumEt;
  
  std::auto_ptr<int>  fTNEBhits;
  std::auto_ptr<std::vector<float> >  fTEBrechitE;
  std::auto_ptr<std::vector<float> >  fTEBrechitPt;
  std::auto_ptr<std::vector<float> >  fTEBrechitEta;
  std::auto_ptr<std::vector<float> >  fTEBrechitPhi;
  std::auto_ptr<std::vector<float> >  fTEBrechitChi2;
  std::auto_ptr<std::vector<float> >  fTEBrechitTime;
  std::auto_ptr<std::vector<float> >  fTEBrechitE4oE1;
  std::auto_ptr<std::vector<float> >  fTEBrechitE2oE9;

  // CSCBeamHalo 
  std::auto_ptr<int>  fTCSCTightHaloID;

  // Vertex information
  std::auto_ptr<int>    fTPrimVtxGood;
  std::auto_ptr<float>  fTPrimVtxx;
  std::auto_ptr<float>  fTPrimVtxy;
  std::auto_ptr<float>  fTPrimVtxz;
  std::auto_ptr<float>  fTPrimVtxRho;
  std::auto_ptr<float>  fTPrimVtxxE;
  std::auto_ptr<float>  fTPrimVtxyE;
  std::auto_ptr<float>  fTPrimVtxzE;
  std::auto_ptr<float>  fTPrimVtxNChi2;
  std::auto_ptr<float>  fTPrimVtxNdof;
  std::auto_ptr<int>    fTPrimVtxIsFake;
  std::auto_ptr<float>  fTPrimVtxPtSum;

  std::auto_ptr<float>  fTBeamspotx;
  std::auto_ptr<float>  fTBeamspoty;
  std::auto_ptr<float>  fTBeamspotz;

  std::auto_ptr<int>  fTNVrtx;
  std::auto_ptr<std::vector<float> >  fTVrtxX;
  std::auto_ptr<std::vector<float> >  fTVrtxY;
  std::auto_ptr<std::vector<float> >  fTVrtxZ;
  std::auto_ptr<std::vector<float> >  fTVrtxXE;
  std::auto_ptr<std::vector<float> >  fTVrtxYE;
  std::auto_ptr<std::vector<float> >  fTVrtxZE;
  std::auto_ptr<std::vector<float> >  fTVrtxNdof;
  std::auto_ptr<std::vector<float> >  fTVrtxChi2;
  std::auto_ptr<std::vector<float> >  fTVrtxNtrks;
  std::auto_ptr<std::vector<float> >  fTVrtxSumPt;
  std::auto_ptr<std::vector<int> >  fTVrtxIsFake;
  
  std::auto_ptr<int>  fTNCaloTowers;

  // Trigger
  static const unsigned int gMaxHltBits = 400;
  static const unsigned int gMaxL1PhysBits = 128;
  static const unsigned int gMaxL1TechBits = 64;

  std::auto_ptr<std::vector<int> >  fTHLTResults;
  std::auto_ptr<std::vector<int> >  fTHLTPrescale;
  std::auto_ptr<std::vector<int> >  fTL1PhysResults;
  std::auto_ptr<std::vector<int> >  fTL1TechResults;

  static const unsigned int gMaxHltNObjs  = 10;
  std::auto_ptr<std::vector<int> >  fTNHLTObjs;
  std::auto_ptr<std::vector<int> >  fTHLTObjectID[10];
  std::auto_ptr<std::vector<float> >  fTHLTObjectPt[10];
  std::auto_ptr<std::vector<float> >  fTHLTObjectEta[10];
  std::auto_ptr<std::vector<float> >  fTHLTObjectPhi[10];

  unsigned int fTNpaths;

  // Flags
  std::auto_ptr<int>  fTGoodEvent;         // 0 for good events, 1 for bad events                     
  std::auto_ptr<int>  fTMaxMuExceed;       // Found more than 20 muons in event (0 is good, 1 is bad) 
  std::auto_ptr<int>  fTMaxElExceed;       // Found more than 20 electrons in event                   
  std::auto_ptr<int>  fTMaxJetExceed;      // Found more than 50 jets in event                        
  std::auto_ptr<int>  fTMaxUncJetExceed;   // Found more than 50 uncorrected jets in event            
  std::auto_ptr<int>  fTMaxTrkExceed;      // Found more than 800 tracks in event                     
  std::auto_ptr<int>  fTMaxPhotonsExceed;  // Found more than 500 photons in event                    
  std::auto_ptr<int>  fTMaxGenLepExceed;   // Found more than 100 genleptons in event                 
  std::auto_ptr<int>  fTMaxGenPhoExceed;   // Found more than 100 genphotons in event                 
  std::auto_ptr<int>  fTMaxGenJetExceed;   // Found more than 100 genjets in event                    
  std::auto_ptr<int>  fTMaxVerticesExceed; // Found more than 25 vertices in event                    
  std::auto_ptr<int>  fTMaxGenPartExceed;  // Found more than 2000 gen particles in event

  // GenLeptons
  std::auto_ptr<int>  fTNGenLeptons;
  std::auto_ptr<std::vector<int> >  fTGenLeptonID;
  std::auto_ptr<std::vector<float> >  fTGenLeptonPt;
  std::auto_ptr<std::vector<float> >  fTGenLeptonEta;
  std::auto_ptr<std::vector<float> >  fTGenLeptonPhi;
  std::auto_ptr<std::vector<int> >  fTGenLeptonMID;
  std::auto_ptr<std::vector<int> >  fTGenLeptonMStatus;
  std::auto_ptr<std::vector<float> >  fTGenLeptonMPt;
  std::auto_ptr<std::vector<float> >  fTGenLeptonMEta;
  std::auto_ptr<std::vector<float> >  fTGenLeptonMPhi;
  std::auto_ptr<std::vector<int> >  fTGenLeptonGMID;
  std::auto_ptr<std::vector<int> >  fTGenLeptonGMStatus;
  std::auto_ptr<std::vector<float> >  fTGenLeptonGMPt;
  std::auto_ptr<std::vector<float> >  fTGenLeptonGMEta;
  std::auto_ptr<std::vector<float> >  fTGenLeptonGMPhi;

  // GenPhotons
  std::auto_ptr<int>  fTNGenPhotons;
  std::auto_ptr<std::vector<float> >  fTGenPhotonPt;
  std::auto_ptr<std::vector<float> >  fTGenPhotonEta;
  std::auto_ptr<std::vector<float> >  fTGenPhotonPhi;
  std::auto_ptr<std::vector<float> >  fTGenPhotonPartonMindR;
  std::auto_ptr<std::vector<int> >  fTGenPhotonMotherID;
  std::auto_ptr<std::vector<int> >  fTGenPhotonMotherStatus;

  // GenJets
  std::auto_ptr<int>  fTNGenJets;
  std::auto_ptr<std::vector<float> >  fTGenJetPt;
  std::auto_ptr<std::vector<float> >  fTGenJetEta;
  std::auto_ptr<std::vector<float> >  fTGenJetPhi;
  std::auto_ptr<std::vector<float> >  fTGenJetE;
  std::auto_ptr<std::vector<float> >  fTGenJetEmE;
  std::auto_ptr<std::vector<float> >  fTGenJetHadE;
  std::auto_ptr<std::vector<float> >  fTGenJetInvE;


  // Muons:
  std::auto_ptr<int>  fTNMus;
  std::auto_ptr<int>  fTNMusTot; // before preselection
  std::auto_ptr<int>  fTNGMus;
  std::auto_ptr<int>  fTNTMus;
  std::auto_ptr<std::vector<int> >  fTMuGood;
  std::auto_ptr<std::vector<int> >  fTMuIsIso;
  std::auto_ptr<std::vector<int> >  fTMuIsGlobalMuon;
  std::auto_ptr<std::vector<int> >  fTMuIsTrackerMuon;
  std::auto_ptr<std::vector<int> >  fTMuIsPFMuon;
  std::auto_ptr<std::vector<float> >  fTMuPx;
  std::auto_ptr<std::vector<float> >  fTMuPy;
  std::auto_ptr<std::vector<float> >  fTMuPz;
  std::auto_ptr<std::vector<float> >  fTMuPt;
  std::auto_ptr<std::vector<float> >  fTMuInnerTkPt;
  std::auto_ptr<std::vector<float> >  fTMuPtE;
  std::auto_ptr<std::vector<float> >  fTMuTkPtE;
  std::auto_ptr<std::vector<float> >  fTMuTkD0E;
  std::auto_ptr<std::vector<float> >  fTMuTkDzE;  
  std::auto_ptr<std::vector<float> >  fTMuE;
  std::auto_ptr<std::vector<float> >  fTMuEt;
  std::auto_ptr<std::vector<float> >  fTMuEta;
  std::auto_ptr<std::vector<float> >  fTMuPhi;
  std::auto_ptr<std::vector<int> >  fTMuCharge;
  // Isolation variables
  std::auto_ptr<std::vector<float> >  fTMuRelIso03;
  std::auto_ptr<std::vector<float> >  fTMuIso03SumPt;
  std::auto_ptr<std::vector<float> >  fTMuIso03EmEt;
  std::auto_ptr<std::vector<float> >  fTMuIso03HadEt;
  std::auto_ptr<std::vector<float> >  fTMuIso03EMVetoEt;
  std::auto_ptr<std::vector<float> >  fTMuIso03HadVetoEt;
  std::auto_ptr<std::vector<float> >  fTMuIso05SumPt;
  std::auto_ptr<std::vector<float> >  fTMuIso05EmEt;
  std::auto_ptr<std::vector<float> >  fTMuIso05HadEt;
  std::auto_ptr<std::vector<float> >  fTMuPfIsoR03ChHad;
  std::auto_ptr<std::vector<float> >  fTMuPfIsoR03NeHad;
  std::auto_ptr<std::vector<float> >  fTMuPfIsoR03Photon;
  std::auto_ptr<std::vector<float> >  fTMuPfIsoR03NeHadHighThresh  ;
  std::auto_ptr<std::vector<float> >  fTMuPfIsoR03PhotonHighThresh ;
  std::auto_ptr<std::vector<float> >  fTMuPfIsoR03SumPUPt;
  std::auto_ptr<std::vector<float> >  fTMuPfIsoR04ChHad;
  std::auto_ptr<std::vector<float> >  fTMuPfIsoR04NeHad;
  std::auto_ptr<std::vector<float> >  fTMuPfIsoR04Photon;
  std::auto_ptr<std::vector<float> >  fTMuPfIsoR04NeHadHighThresh  ;
  std::auto_ptr<std::vector<float> >  fTMuPfIsoR04PhotonHighThresh ;
  std::auto_ptr<std::vector<float> >  fTMuPfIsoR04SumPUPt;
  std::auto_ptr<std::vector<float> >  fTMuEem;
  std::auto_ptr<std::vector<float> >  fTMuEhad;
  // PF isolation variables
  std::auto_ptr<std::vector<float> >  fTMuPfIsosCustom[gMaxNPfIsoTags];
  // Impact parameters
  std::auto_ptr<std::vector<float> >  fTMuD0BS;
  std::auto_ptr<std::vector<float> >  fTMuD0PV;
  std::auto_ptr<std::vector<float> >  fTMuD0E;
  std::auto_ptr<std::vector<float> >  fTMuDzBS;
  std::auto_ptr<std::vector<float> >  fTMuDzPV;
  std::auto_ptr<std::vector<float> >  fTMuDzE;
  // Mu ID variables
  std::auto_ptr<std::vector<float> >  fTMuNChi2;
  std::auto_ptr<std::vector<int> >  fTMuNGlHits;
  std::auto_ptr<std::vector<int> >  fTMuNMuHits;
  std::auto_ptr<std::vector<int> >  fTMuNTkHits;
  std::auto_ptr<std::vector<int> >  fTMuNPxHits;
  std::auto_ptr<std::vector<float> >  fTMuInnerTkNChi2;
  std::auto_ptr<std::vector<int> >  fTMuNSiLayers;
  std::auto_ptr<std::vector<int> >  fTMuNMatches;
  std::auto_ptr<std::vector<int> >  fTMuNMatchedStations;
  std::auto_ptr<std::vector<int> >  fTMuNChambers;
  std::auto_ptr<std::vector<float> >  fTMuIsoMVA;
  std::auto_ptr<std::vector<float> >  fTMuCaloComp;
  std::auto_ptr<std::vector<float> >  fTMuSegmComp;
  std::auto_ptr<std::vector<int> >  fTMuIsGMPT;
  std::auto_ptr<std::vector<int> >  fTMuIsGMTkChiComp;
  std::auto_ptr<std::vector<int> >  fTMuIsGMStaChiComp;
  std::auto_ptr<std::vector<int> >  fTMuIsGMTkKinkTight;
  std::auto_ptr<std::vector<int> >  fTMuIsAllStaMuons;
  std::auto_ptr<std::vector<int> >  fTMuIsAllTrkMuons;
  std::auto_ptr<std::vector<int> >  fTMuIsTrkMuonArbitrated;
  std::auto_ptr<std::vector<int> >  fTMuIsAllArbitrated;
  std::auto_ptr<std::vector<int> >  fTMuIsTMLSLoose;
  std::auto_ptr<std::vector<int> >  fTMuIsTMLSTight;
  std::auto_ptr<std::vector<int> >  fTMuIsTM2DCompLoose;
  std::auto_ptr<std::vector<int> >  fTMuIsTM2DCompTight;
  std::auto_ptr<std::vector<int> >  fTMuIsTMOneStationLoose;
  std::auto_ptr<std::vector<int> >  fTMuIsTMOneStationTight;
  std::auto_ptr<std::vector<int> >  fTMuIsTMLSOptLowPtLoose;
  std::auto_ptr<std::vector<int> >  fTMuIsTMLSAngLoose;
  std::auto_ptr<std::vector<int> >  fTMuIsTMLSAngTight;
  std::auto_ptr<std::vector<int> >  fTMuIsTMOneStationAngTight;
  std::auto_ptr<std::vector<int> >  fTMuIsTMOneStationAngLoose;
  // Generator information
  std::auto_ptr<std::vector<int> >  fTMuGenID;
  std::auto_ptr<std::vector<int> >  fTMuGenStatus;
  std::auto_ptr<std::vector<float> >  fTMuGenPt;
  std::auto_ptr<std::vector<float> >  fTMuGenEta;
  std::auto_ptr<std::vector<float> >  fTMuGenPhi;
  std::auto_ptr<std::vector<float> >  fTMuGenE;
  std::auto_ptr<std::vector<int> >  fTMuGenMID;
  std::auto_ptr<std::vector<int> >  fTMuGenMStatus;
  std::auto_ptr<std::vector<float> >  fTMuGenMPt;
  std::auto_ptr<std::vector<float> >  fTMuGenMEta;
  std::auto_ptr<std::vector<float> >  fTMuGenMPhi;
  std::auto_ptr<std::vector<float> >  fTMuGenME;
  std::auto_ptr<std::vector<int> >  fTMuGenGMID;
  std::auto_ptr<std::vector<int> >  fTMuGenGMStatus;
  std::auto_ptr<std::vector<float> >  fTMuGenGMPt;
  std::auto_ptr<std::vector<float> >  fTMuGenGMEta;
  std::auto_ptr<std::vector<float> >  fTMuGenGMPhi;
  std::auto_ptr<std::vector<float> >  fTMuGenGME;

  //- Electrons:
  std::auto_ptr<int>  fTNEles;
  std::auto_ptr<int>  fTNElesTot;
  std::auto_ptr<std::vector<int> >  fTElGood;
  std::auto_ptr<std::vector<int> >  fTElIsIso;
  std::auto_ptr<std::vector<int> >  fTElChargeMisIDProb;
  // Kinematics
  std::auto_ptr<std::vector<float> >  fTElPx;
  std::auto_ptr<std::vector<float> >  fTElPy;
  std::auto_ptr<std::vector<float> >  fTElPz;
  std::auto_ptr<std::vector<float> >  fTElPt;
  std::auto_ptr<std::vector<float> >  fTElPtE;
  std::auto_ptr<std::vector<float> >  fTElE;
  std::auto_ptr<std::vector<float> >  fTElEt;
  std::auto_ptr<std::vector<float> >  fTElEta;
  std::auto_ptr<std::vector<float> >  fTElTheta;
  std::auto_ptr<std::vector<float> >  fTElSCEta;
  std::auto_ptr<std::vector<float> >  fTElPhi;
  std::auto_ptr<std::vector<float> >  fTElGsfTkPt;
  std::auto_ptr<std::vector<float> >  fTElGsfTkEta;
  std::auto_ptr<std::vector<float> >  fTElGsfTkPhi;
  std::auto_ptr<std::vector<float> >  fTElTrkMomentumError;
  std::auto_ptr<std::vector<float> >  fTElEcalEnergyError;
  std::auto_ptr<std::vector<float> >  fTElEleMomentumError;
  std::auto_ptr<std::vector<int> >  fTElNBrems;
  // Impact Parameter
  std::auto_ptr<std::vector<float> >  fTElD0BS;
  std::auto_ptr<std::vector<float> >  fTElD0PV;
  std::auto_ptr<std::vector<float> >  fTElD0E;
  std::auto_ptr<std::vector<float> >  fTElDzBS;
  std::auto_ptr<std::vector<float> >  fTElDzPV;
  std::auto_ptr<std::vector<float> >  fTElDzE;
  // Isolation variables
  std::auto_ptr<std::vector<float> >  fTElRelIso03;
  std::auto_ptr<std::vector<float> >  fTElRelIso04;
  std::auto_ptr<std::vector<float> >  fTElPfIsoChHad03;
  std::auto_ptr<std::vector<float> >  fTElPfIsoNeHad03;
  std::auto_ptr<std::vector<float> >  fTElPfIsoPhoton03;
  std::auto_ptr<std::vector<float> >  fTElDR03TkSumPt;
  std::auto_ptr<std::vector<float> >  fTElDR04TkSumPt;
  std::auto_ptr<std::vector<float> >  fTElDR03EcalRecHitSumEt;
  std::auto_ptr<std::vector<float> >  fTElDR04EcalRecHitSumEt;
  std::auto_ptr<std::vector<float> >  fTElDR03HcalTowerSumEt;
  std::auto_ptr<std::vector<float> >  fTElDR04HcalTowerSumEt;
  std::auto_ptr<std::vector<float> >  fTElNChi2;
  std::auto_ptr<std::vector<float> >  fTElKfTrkchi2;
  std::auto_ptr<std::vector<float> >  fTElKfTrkhits;
  // PF Isolation Variables
  std::auto_ptr<std::vector<float> >   fTElPfIsosCustom[gMaxNPfIsoTags];
  std::auto_ptr<std::vector<float> >  fTElPfIsosEvent[gMaxNPfIsoTags];
  // Electron ID variables
  std::auto_ptr<std::vector<int> >  fTElCharge;
  std::auto_ptr<std::vector<int> >  fTElCInfoIsGsfCtfCons;
  std::auto_ptr<std::vector<int> >  fTElCInfoIsGsfCtfScPixCons;
  std::auto_ptr<std::vector<int> >  fTElCInfoIsGsfScPixCons;
  std::auto_ptr<std::vector<int> >  fTElScPixCharge;
  std::auto_ptr<std::vector<float> >  fTElClosestCtfTrackPt;
  std::auto_ptr<std::vector<float> >  fTElClosestCtfTrackEta;
  std::auto_ptr<std::vector<float> >  fTElClosestCtfTrackPhi;
  std::auto_ptr<std::vector<int> >  fTElClosestCtfTrackCharge;
  std::auto_ptr<std::vector<float> >  fTElIDMva;
  std::auto_ptr<std::vector<float> >  fTElIDMVATrig;
  std::auto_ptr<std::vector<float> >  fTElIDMVANoTrig;
  std::auto_ptr<std::vector<int> >  fTElIDTight;
  std::auto_ptr<std::vector<int> >  fTElIDLoose;
  std::auto_ptr<std::vector<int> >  fTElIDRobustTight;
  std::auto_ptr<std::vector<int> >  fTElIDRobustLoose;
  std::auto_ptr<std::vector<int> >  fTElIDsimpleWPrelIso;
  std::auto_ptr<std::vector<int> >  fTElIDsimpleWP80relIso;
  std::auto_ptr<std::vector<int> >  fTElIDsimpleWP85relIso;
  std::auto_ptr<std::vector<int> >  fTElIDsimpleWP90relIso;
  std::auto_ptr<std::vector<int> >  fTElIDsimpleWP95relIso;
  std::auto_ptr<std::vector<int> >  fTElInGap; // seed crystal next to a gap
  std::auto_ptr<std::vector<int> >  fTElEcalDriven;
  std::auto_ptr<std::vector<int> >  fTElTrackerDriven;
  std::auto_ptr<std::vector<int> >  fTElBasicClustersSize;
  std::auto_ptr<std::vector<float> >  fTElfbrem;
  std::auto_ptr<std::vector<float> >  fTElEoPout;
  std::auto_ptr<std::vector<float> >  fTElIoEmIoP;
  std::auto_ptr<std::vector<float> >  fTElHcalOverEcal;
  std::auto_ptr<std::vector<float> >  fTElE1x5;                       // 5x5 arround seed                                     
  std::auto_ptr<std::vector<float> >  fTElE5x5;                       // 5x5 arround seed                                     
  std::auto_ptr<std::vector<float> >  fTElE2x5Max;                    // 2x5 arround seed                                     
  std::auto_ptr<std::vector<float> >  fTElR9;                  
  std::auto_ptr<std::vector<float> >  fTElPreShowerOverRaw;                  
  std::auto_ptr<std::vector<float> >  fTElScEtaWidth;                  
  std::auto_ptr<std::vector<float> >  fTElScPhiWidth;                  
  std::auto_ptr<std::vector<float> >  fTElSigmaIetaIeta;              // shower shape covariance                              
  std::auto_ptr<std::vector<float> >  fTElSigmaIphiIphi;              // shower shape covariance                              
  std::auto_ptr<std::vector<float> >  fTElDeltaPhiSeedClusterAtCalo;  // Dphi (seed-track) at calo from p_out                 
  std::auto_ptr<std::vector<float> >  fTElDeltaEtaSeedClusterAtCalo;  // outermost track state extrapolated at calo           
  std::auto_ptr<std::vector<float> >  fTElDeltaPhiSuperClusterAtVtx;  // Dphi (sc-track) at calo extrapolated from p_in       
  std::auto_ptr<std::vector<float> >  fTElDeltaEtaSuperClusterAtVtx;  // Deta (sc-track) at calo extrapolated from p_in       
  std::auto_ptr<std::vector<float> >  fTElCaloEnergy;                 // caloEnergy() = supercluster energy 99.9% of the time 
  std::auto_ptr<std::vector<float> >  fTElTrkMomAtVtx;                // trackMomentumAtVtx().R()                             
  std::auto_ptr<std::vector<float> >  fTElESuperClusterOverP;         // Esc/Pin                                              
  std::auto_ptr<std::vector<int> >  fTElNumberOfMissingInnerHits;
  std::auto_ptr<std::vector<int> >  fTElSCindex;
  std::auto_ptr<std::vector<bool> >  fTElPassConversionVeto;
  std::auto_ptr<std::vector<float> >  fTElConvPartnerTrkDist;
  std::auto_ptr<std::vector<float> >  fTElConvPartnerTrkDCot;
  std::auto_ptr<std::vector<float> >  fTElConvPartnerTrkPt;
  std::auto_ptr<std::vector<float> >  fTElConvPartnerTrkEta;
  std::auto_ptr<std::vector<float> >  fTElConvPartnerTrkPhi;
  std::auto_ptr<std::vector<float> >  fTElConvPartnerTrkCharge;
  // Spike cleaning
  std::auto_ptr<std::vector<int> >  fTElScSeedSeverity;
  std::auto_ptr<std::vector<float> >  fTElE1OverE9;
  std::auto_ptr<std::vector<float> >  fTElS4OverS1;
  // Gen info
  std::auto_ptr<std::vector<int> >  fTElGenID;
  std::auto_ptr<std::vector<int> >  fTElGenStatus;
  std::auto_ptr<std::vector<float> >  fTElGenPt;
  std::auto_ptr<std::vector<float> >  fTElGenEta;
  std::auto_ptr<std::vector<float> >  fTElGenPhi;
  std::auto_ptr<std::vector<float> >  fTElGenE;
  std::auto_ptr<std::vector<int> >  fTElGenMID;
  std::auto_ptr<std::vector<int> >  fTElGenMStatus;
  std::auto_ptr<std::vector<float> >  fTElGenMPt;
  std::auto_ptr<std::vector<float> >  fTElGenMEta;
  std::auto_ptr<std::vector<float> >  fTElGenMPhi;
  std::auto_ptr<std::vector<float> >  fTElGenME;
  std::auto_ptr<std::vector<int> >  fTElGenGMID;
  std::auto_ptr<std::vector<int> >  fTElGenGMStatus;
  std::auto_ptr<std::vector<float> >  fTElGenGMPt;
  std::auto_ptr<std::vector<float> >  fTElGenGMEta;
  std::auto_ptr<std::vector<float> >  fTElGenGMPhi;
  std::auto_ptr<std::vector<float> >  fTElGenGME;

  //- Photons:
  std::auto_ptr<int>  fTNPhotons;
  std::auto_ptr<int>  fTNPhotonsTot;
  std::auto_ptr<std::vector<bool> >  fTPhoPassConversionVeto;
  std::auto_ptr<std::vector<int> >  fTPhoGood;
  std::auto_ptr<std::vector<int> >  fTPhoIsIso;
  std::auto_ptr<std::vector<float> >  fTPhoPt;
  std::auto_ptr<std::vector<float> >  fTPhoPx;
  std::auto_ptr<std::vector<float> >  fTPhoPy;
  std::auto_ptr<std::vector<float> >  fTPhoPz;
  std::auto_ptr<std::vector<float> >  fTPhoEta;
  std::auto_ptr<std::vector<float> >  fTPhoPhi;
  std::auto_ptr<std::vector<float> >  fTPhoEnergy;
  std::auto_ptr<std::vector<float> >  fTPhoIso03Ecal;
  std::auto_ptr<std::vector<float> >  fTPhoIso03Hcal;
  std::auto_ptr<std::vector<float> >  fTPhoIso03TrkSolid;
  std::auto_ptr<std::vector<float> >  fTPhoIso03TrkHollow;
  std::auto_ptr<std::vector<float> >  fTPhoIso03;
  std::auto_ptr<std::vector<float> >  fTPhoIso04Ecal;
  std::auto_ptr<std::vector<float> >  fTPhoIso04Hcal;
  std::auto_ptr<std::vector<float> >  fTPhoIso04TrkSolid;
  std::auto_ptr<std::vector<float> >  fTPhoIso04TrkHollow;
  std::auto_ptr<std::vector<float> >  fTPhoIso04;
  std::auto_ptr<std::vector<float> >  fTPhoR9;
  std::auto_ptr<std::vector<float> >  fTPhoCaloPositionX;
  std::auto_ptr<std::vector<float> >  fTPhoCaloPositionY;
  std::auto_ptr<std::vector<float> >  fTPhoCaloPositionZ;
  std::auto_ptr<std::vector<float> >  fTPhoHoverE;
  std::auto_ptr<std::vector<float> >  fTPhoH1overE;
  std::auto_ptr<std::vector<float> >  fTPhoH2overE;
  std::auto_ptr<std::vector<float> >  fTPhoSigmaIetaIeta;
  std::auto_ptr<std::vector<float> >  fTPhoSCRawEnergy;
  std::auto_ptr<std::vector<float> >  fTPhoSCEtaWidth;
  std::auto_ptr<std::vector<float> >  fTPhoSCSigmaPhiPhi;
  std::auto_ptr<std::vector<int> >  fTPhoHasPixSeed;
  std::auto_ptr<std::vector<int> >  fTPhoHasConvTrks;
  // Spike cleaning
  std::auto_ptr<std::vector<int> >  fTPhoScSeedSeverity;
  std::auto_ptr<std::vector<float> >  fTPhoE1OverE9;
  std::auto_ptr<std::vector<float> >  fTPhoS4OverS1;
  // ID
  std::auto_ptr<std::vector<float> >  fTPhoSigmaEtaEta;
  std::auto_ptr<std::vector<float> >  fTPhoHCalIso2012ConeDR03;
  std::auto_ptr<std::vector<float> >  fTPhoNewIsoPFCharged;
  std::auto_ptr<std::vector<float> >  fTPhoNewIsoPFPhoton;
  std::auto_ptr<std::vector<float> >  fTPhoNewIsoPFNeutral;
  std::auto_ptr<std::vector<float> >  fTPhoE1x5;
  std::auto_ptr<std::vector<float> >  fTPhoE2x5;
  std::auto_ptr<std::vector<float> >  fTPhoE3x3;
  std::auto_ptr<std::vector<float> >  fTPhoE5x5;
  std::auto_ptr<std::vector<float> >  fTPhomaxEnergyXtal;
  std::auto_ptr<std::vector<float> >  fTPhoIso03HcalDepth1;
  std::auto_ptr<std::vector<float> >  fTPhoIso03HcalDepth2;
  std::auto_ptr<std::vector<float> >  fTPhoIso04HcalDepth1;
  std::auto_ptr<std::vector<float> >  fTPhoIso04HcalDepth2;
  std::auto_ptr<std::vector<int> >  fTPhoIso03nTrksSolid;
  std::auto_ptr<std::vector<int> >  fTPhoIso03nTrksHollow;
  std::auto_ptr<std::vector<int> >  fTPhoIso04nTrksSolid;
  std::auto_ptr<std::vector<int> >  fTPhoIso04nTrksHollow;
  std::auto_ptr<std::vector<int> >  fTPhoisEB;
  std::auto_ptr<std::vector<int> >  fTPhoisEE;
  std::auto_ptr<std::vector<int> >  fTPhoisEBEtaGap;
  std::auto_ptr<std::vector<int> >  fTPhoisEBPhiGap;
  std::auto_ptr<std::vector<int> >  fTPhoisEERingGap;
  std::auto_ptr<std::vector<int> >  fTPhoisEEDeeGap;
  std::auto_ptr<std::vector<int> >  fTPhoisEBEEGap;
  std::auto_ptr<std::vector<int> >  fTPhoisPFlowPhoton;
  std::auto_ptr<std::vector<int> >  fTPhoisStandardPhoton;
  std::auto_ptr<std::vector<int> >  fTPhoMCmatchindex;
  std::auto_ptr<std::vector<int> >  fTPhoMCmatchexitcode;
  std::auto_ptr<std::vector<float> >  fTPhoChargedHadronIso;
  std::auto_ptr<std::vector<float> >  fTPhoNeutralHadronIso;
  std::auto_ptr<std::vector<float> >  fTPhoPhotonIso;
  std::auto_ptr<std::vector<int> >  fTPhoisPFPhoton;
  std::auto_ptr<std::vector<int> >  fTPhoisPFElectron;
  std::auto_ptr<std::vector<int> >  fTPhotSCindex;
  std::auto_ptr<std::vector<float> >  fTPhoCone04PhotonIsodR0dEta0pt0;
  std::auto_ptr<std::vector<float> >  fTPhoCone04PhotonIsodR0dEta0pt5;
  std::auto_ptr<std::vector<float> >  fTPhoCone04PhotonIsodR8dEta0pt0;
  std::auto_ptr<std::vector<float> >  fTPhoCone04PhotonIsodR8dEta0pt5;
  std::auto_ptr<std::vector<float> >  fTPhoCone01PhotonIsodR045EB070EEdEta015pt08EB1EEmvVtx;
  std::auto_ptr<std::vector<float> >  fTPhoCone02PhotonIsodR045EB070EEdEta015pt08EB1EEmvVtx;
  std::auto_ptr<std::vector<float> >  fTPhoCone03PhotonIsodR045EB070EEdEta015pt08EB1EEmvVtx;
  std::auto_ptr<std::vector<float> >  fTPhoCone04PhotonIsodR045EB070EEdEta015pt08EB1EEmvVtx;
  std::auto_ptr<std::vector<float> >  fTPhoCone04NeutralHadronIsodR0dEta0pt0;
  std::auto_ptr<std::vector<float> >  fTPhoCone04NeutralHadronIsodR0dEta0pt5;
  std::auto_ptr<std::vector<float> >  fTPhoCone04NeutralHadronIsodR0dEta0pt0nocracks;
  std::auto_ptr<std::vector<float> >  fTPhoCone04NeutralHadronIsodR0dEta0pt5nocracks;
  std::auto_ptr<std::vector<float> >  fTPhoCone04NeutralHadronIsodR7dEta0pt0;
  std::auto_ptr<std::vector<float> >  fTPhoCone04NeutralHadronIsodR7dEta0pt5;
  std::auto_ptr<std::vector<float> >  fTPhoCone01NeutralHadronIsodR0dEta0pt0mvVtx;
  std::auto_ptr<std::vector<float> >  fTPhoCone02NeutralHadronIsodR0dEta0pt0mvVtx;
  std::auto_ptr<std::vector<float> >  fTPhoCone03NeutralHadronIsodR0dEta0pt0mvVtx;
  std::auto_ptr<std::vector<float> >  fTPhoCone04NeutralHadronIsodR0dEta0pt0mvVtx;
  std::auto_ptr<std::vector<float> >  fTPhoCone04ChargedHadronIsodR0dEta0pt0dz0old;
  std::auto_ptr<std::vector<float> >  fTPhoCone04ChargedHadronIsodR0dEta0pt0PFnoPUold;
  std::auto_ptr<std::vector<float> >  fTPhoCone04ChargedHadronIsodR015dEta0pt0dz0old;
  std::auto_ptr<std::vector<float> >  fTPhoCone04ChargedHadronIsodR015dEta0pt0PFnoPUold;
  std::auto_ptr<std::vector<float> >  fTPhoCone01ChargedHadronIsodR0dEta0pt0dz0;
  std::auto_ptr<std::vector<float> >  fTPhoCone01ChargedHadronIsodR0dEta0pt0dz1dxy01;
  std::auto_ptr<std::vector<float> >  fTPhoCone01ChargedHadronIsodR0dEta0pt0PFnoPU;
  std::auto_ptr<std::vector<float> >  fTPhoCone01ChargedHadronIsodR015dEta0pt0dz0;
  std::auto_ptr<std::vector<float> >  fTPhoCone01ChargedHadronIsodR015dEta0pt0dz1dxy01;
  std::auto_ptr<std::vector<float> >  fTPhoCone01ChargedHadronIsodR015dEta0pt0PFnoPU;
  std::auto_ptr<std::vector<float> >  fTPhoCone02ChargedHadronIsodR0dEta0pt0dz0;
  std::auto_ptr<std::vector<float> >  fTPhoCone02ChargedHadronIsodR0dEta0pt0dz1dxy01;
  std::auto_ptr<std::vector<float> >  fTPhoCone02ChargedHadronIsodR0dEta0pt0PFnoPU;
  std::auto_ptr<std::vector<float> >  fTPhoCone02ChargedHadronIsodR015dEta0pt0dz0;
  std::auto_ptr<std::vector<float> >  fTPhoCone02ChargedHadronIsodR015dEta0pt0dz1dxy01;
  std::auto_ptr<std::vector<float> >  fTPhoCone02ChargedHadronIsodR015dEta0pt0PFnoPU;
  std::auto_ptr<std::vector<float> >  fTPhoCone03ChargedHadronIsodR0dEta0pt0dz0;
  std::auto_ptr<std::vector<float> >  fTPhoCone03ChargedHadronIsodR0dEta0pt0dz1dxy01;
  std::auto_ptr<std::vector<float> >  fTPhoCone03ChargedHadronIsodR0dEta0pt0PFnoPU;
  std::auto_ptr<std::vector<float> >  fTPhoCone03ChargedHadronIsodR015dEta0pt0dz0;
  std::auto_ptr<std::vector<float> >  fTPhoCone03ChargedHadronIsodR015dEta0pt0dz1dxy01;
  std::auto_ptr<std::vector<float> >  fTPhoCone03ChargedHadronIsodR015dEta0pt0PFnoPU;
  std::auto_ptr<std::vector<float> >  fTPhoCone04ChargedHadronIsodR0dEta0pt0dz0;
  std::auto_ptr<std::vector<float> >  fTPhoCone04ChargedHadronIsodR0dEta0pt0dz1dxy01;
  std::auto_ptr<std::vector<float> >  fTPhoCone04ChargedHadronIsodR0dEta0pt0PFnoPU;
  std::auto_ptr<std::vector<float> >  fTPhoCone04ChargedHadronIsodR015dEta0pt0dz0;
  std::auto_ptr<std::vector<float> >  fTPhoCone04ChargedHadronIsodR015dEta0pt0dz1dxy01;
  std::auto_ptr<std::vector<float> >  fTPhoCone04ChargedHadronIsodR015dEta0pt0PFnoPU;

  TVector3 pho_conv_vtx[gMaxNPhotons];
  TVector3 pho_conv_refitted_momentum[gMaxNPhotons];
  TVector3 conv_vtx[gMaxNPhotons];
  TVector3 conv_refitted_momentum[gMaxNPhotons];

  std::auto_ptr<std::vector<bool> >  fTPhoConvValidVtx;
  std::auto_ptr<std::vector<int> >   fTPhoConvNtracks;
  std::auto_ptr<std::vector<float> > fTPhoConvChi2Probability;
  std::auto_ptr<std::vector<float> > fTPhoConvEoverP;
  
  std::auto_ptr<int> fTNconv;
  std::auto_ptr<std::vector<bool> >  fTConvValidVtx;
  std::auto_ptr<std::vector<int> >   fTConvNtracks;
  std::auto_ptr<std::vector<float> > fTConvChi2Probability;
  std::auto_ptr<std::vector<float> > fTConvEoverP;
  std::auto_ptr<std::vector<float> > fTConvZofPrimVtxFromTrks;

  TVector3 gv_pos[gMaxNGenVtx];
  TVector3 gv_p3[gMaxNGenVtx];
  std::auto_ptr<int> fTNgv;
  std::auto_ptr<std::vector<float> > fTgvSumPtHi;
  std::auto_ptr<std::vector<float> > fTgvSumPtLo;
  std::auto_ptr<std::vector<int> >   fTgvNTkHi;
  std::auto_ptr<std::vector<int> >   fTgvNTkLo;


  //- Superclusters:
  std::auto_ptr<int> fTNGoodSuperClusters;
  std::auto_ptr<std::vector<float> >  fTGoodSCEnergy;
  std::auto_ptr<std::vector<float> >  fTGoodSCEta;
  std::auto_ptr<std::vector<float> >  fTGoodSCPhi;

  std::auto_ptr<int>  fTNSuperClusters;
  std::auto_ptr<std::vector<float> >  fTSCRaw;
  std::auto_ptr<std::vector<float> >  fTSCPre;
  std::auto_ptr<std::vector<float> >  fTSCEnergy;
  std::auto_ptr<std::vector<float> >  fTSCEta;
  std::auto_ptr<std::vector<float> >  fTSCPhi;
  std::auto_ptr<std::vector<float> >  fTSCPhiWidth;
  std::auto_ptr<std::vector<float> >  fTSCEtaWidth;
  std::auto_ptr<std::vector<float> >  fTSCBrem;
  std::auto_ptr<std::vector<float> >  fTSCR9;
  std::auto_ptr<std::vector<float> >  fTSCcrackcorrseed;
  std::auto_ptr<std::vector<float> >  fTSCcrackcorr;
  std::auto_ptr<std::vector<float> >  fTSClocalcorrseed;
  std::auto_ptr<std::vector<float> >  fTSClocalcorr;
  std::auto_ptr<std::vector<float> >  fTSCcrackcorrseedfactor;
  std::auto_ptr<std::vector<float> >  fTSClocalcorrseedfactor;

  //- Jets:
  std::auto_ptr<int>  fTNJets;
  std::auto_ptr<int>  fTNJetsTot;
  std::auto_ptr<std::vector<int> >  fTJGood;
  std::auto_ptr<std::vector<float> >  fTJPx;
  std::auto_ptr<std::vector<float> >  fTJPy;
  std::auto_ptr<std::vector<float> >  fTJPz;
  std::auto_ptr<std::vector<float> >  fTJPt;
  std::auto_ptr<std::vector<float> >  fTJE;
  std::auto_ptr<std::vector<float> >  fTJEt;
  std::auto_ptr<std::vector<float> >  fTJEta;
  std::auto_ptr<std::vector<float> >  fTJPhi;
  std::auto_ptr<std::vector<float> >  fTJEcorr;
  std::auto_ptr<std::vector<float> >  fTJArea;
  
  std::auto_ptr<std::vector<float> >  fTJEtaRms;
  std::auto_ptr<std::vector<float> >  fTJPhiRms;
  
  std::auto_ptr<std::vector<int> >  fTJNConstituents;
  std::auto_ptr<std::vector<int> >  fTJNAssoTracks;
  std::auto_ptr<std::vector<int> >  fTJNNeutrals;
  std::auto_ptr<std::vector<float> > fTJChargedEmFrac;
  std::auto_ptr<std::vector<float> > fTJNeutralEmFrac;
  std::auto_ptr<std::vector<float> > fTJChargedHadFrac;
  std::auto_ptr<std::vector<float> > fTJNeutralHadFrac;
  std::auto_ptr<std::vector<float> > fTJChargedMuEnergyFrac;
  std::auto_ptr<std::vector<float> > fTJPhoFrac;
  std::auto_ptr<std::vector<float> > fTJHFHadFrac;
  std::auto_ptr<std::vector<float> > fTJHFEMFrac;
  std::auto_ptr<std::vector<float> > fTJPtD;
  std::auto_ptr<std::vector<float> > fTJRMSCand;
  std::auto_ptr<std::vector<float> >  fTJeMinDR;

  std::auto_ptr<std::vector<float> >  fTJbTagProb[gMaxNBtags];
  std::auto_ptr<std::vector<int>   >  fTJPartonFlavour;

  std::auto_ptr<std::vector<float> >  fTJMass;
  std::auto_ptr<std::vector<float> >  fTJBetaStar;
  std::auto_ptr<std::vector<float> >  fTJtrk1px;
  std::auto_ptr<std::vector<float> >  fTJtrk1py;
  std::auto_ptr<std::vector<float> >  fTJtrk1pz;
  std::auto_ptr<std::vector<float> >  fTJtrk2px;
  std::auto_ptr<std::vector<float> >  fTJtrk2py;
  std::auto_ptr<std::vector<float> >  fTJtrk2pz;
  std::auto_ptr<std::vector<float> >  fTJtrk3px;
  std::auto_ptr<std::vector<float> >  fTJtrk3py;
  std::auto_ptr<std::vector<float> >  fTJtrk3pz;

  std::auto_ptr<std::vector<float> >  fTJVtxx;
  std::auto_ptr<std::vector<float> >  fTJVtxy;
  std::auto_ptr<std::vector<float> >  fTJVtxz;
  std::auto_ptr<std::vector<float> >  fTJVtxExx;
  std::auto_ptr<std::vector<float> >  fTJVtxEyx;
  std::auto_ptr<std::vector<float> >  fTJVtxEyy;
  std::auto_ptr<std::vector<float> >  fTJVtxEzy;
  std::auto_ptr<std::vector<float> >  fTJVtxEzz;
  std::auto_ptr<std::vector<float> >  fTJVtxEzx;
  std::auto_ptr<std::vector<float> >  fTJVtxNChi2;

  std::auto_ptr<std::vector<int> >  fTJGenJetIndex;

  //- Tracks:
  std::auto_ptr<int>  fTNTracks;
  std::auto_ptr<int>  fTNTracksTot;
  std::auto_ptr<std::vector<int> >  fTTrkGood;
  std::auto_ptr<std::vector<float> >  fTTrkPt;
  std::auto_ptr<std::vector<float> >  fTTrkEta;
  std::auto_ptr<std::vector<float> >  fTTrkPhi;
  std::auto_ptr<std::vector<float> >  fTTrkNChi2;
  std::auto_ptr<std::vector<float> >  fTTrkNHits;
  std::auto_ptr<std::vector<float> >  fTTrkVtxDz;
  std::auto_ptr<std::vector<float> >  fTTrkVtxDxy;

  //- (M)E(T):
  std::auto_ptr<float>  fTTrkPtSumx;
  std::auto_ptr<float>  fTTrkPtSumy;
  std::auto_ptr<float>  fTTrkPtSum;
  std::auto_ptr<float>  fTTrkPtSumPhi;
  std::auto_ptr<float>  fTSumEt;
  std::auto_ptr<float>  fTECALSumEt;
  std::auto_ptr<float>  fTHCALSumEt;
  std::auto_ptr<float>  fTECALEsumx;
  std::auto_ptr<float>  fTECALEsumy;
  std::auto_ptr<float>  fTECALEsumz;
  std::auto_ptr<float>  fTECALMET;
  std::auto_ptr<float>  fTECALMETPhi;
  std::auto_ptr<float>  fTECALMETEta;
  std::auto_ptr<float>  fTHCALEsumx;
  std::auto_ptr<float>  fTHCALEsumy;
  std::auto_ptr<float>  fTHCALEsumz;
  std::auto_ptr<float>  fTHCALMET;
  std::auto_ptr<float>  fTHCALMETPhi;
  std::auto_ptr<float>  fTHCALMETeta;
  std::auto_ptr<float>  fTRawMET;
  std::auto_ptr<float>  fTRawMETpx;
  std::auto_ptr<float>  fTRawMETpy;
  std::auto_ptr<float>  fTRawMETphi;
  std::auto_ptr<float>  fTRawMETemEtFrac;
  std::auto_ptr<float>  fTRawMETemEtInEB;
  std::auto_ptr<float>  fTRawMETemEtInEE;
  std::auto_ptr<float>  fTRawMETemEtInHF;
  std::auto_ptr<float>  fTRawMEThadEtFrac;
  std::auto_ptr<float>  fTRawMEThadEtInHB;
  std::auto_ptr<float>  fTRawMEThadEtInHE;
  std::auto_ptr<float>  fTRawMEThadEtInHF;
  std::auto_ptr<float>  fTRawMETSignificance;
  std::auto_ptr<float>  fTGenMET;
  std::auto_ptr<float>  fTGenMETpx;
  std::auto_ptr<float>  fTGenMETpy;
  std::auto_ptr<float>  fTGenMETphi;
  std::auto_ptr<float>  fTTCMET;
  std::auto_ptr<float>  fTTCMETpx;
  std::auto_ptr<float>  fTTCMETpy;
  std::auto_ptr<float>  fTTCMETphi;
  std::auto_ptr<float>  fTTCMETSignificance;
  std::auto_ptr<float>  fTMuJESCorrMET;
  std::auto_ptr<float>  fTMuJESCorrMETpx;
  std::auto_ptr<float>  fTMuJESCorrMETpy;
  std::auto_ptr<float>  fTMuJESCorrMETphi;
  std::auto_ptr<float>  fTPFMET;
  std::auto_ptr<float>  fTPFMETpx;
  std::auto_ptr<float>  fTPFMETpy;
  std::auto_ptr<float>  fTPFMETphi;
  std::auto_ptr<float>  fTPFMETSignificance;
  std::auto_ptr<float>  fTPFSumEt;
  std::auto_ptr<float>  fTMETR12;
  std::auto_ptr<float>  fTMETR21;
};


