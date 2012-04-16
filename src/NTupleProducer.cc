// -*- C++ -*-
//
// Package:    NTupleProducer
// Class:      NTupleProducer
//
/* class NTupleProducer
   NTupleProducer.cc
   AnalysisExamples/NTupleProducer/src/NTupleProducer.cc
   Description: Template to produce NTuples for ETH SUSY Analysis

   Implementation:

*/
//
// Original Author:  Benjamin Stieger
//         Created:  Wed Sep  2 16:43:05 CET 2009
// $Id: NTupleProducer.cc,v 1.146.2.4 2012/04/04 12:21:36 fronga Exp $
//
//


// system include files
#include <memory>
#include <iostream>
#include <cmath>
#include <sys/types.h>
#include <dirent.h>
#include <errno.h>
#include <sstream>
#include <map>

// ROOT includes
#include "TLorentzVector.h"

// Framework include files
#include "FWCore/Framework/interface/TriggerNamesService.h"
#include "FWCore/Framework/interface/MakerMacros.h"

// Utilities
#include "CommonTools/Utils/interface/PtComparator.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"

#include "RecoLocalCalo/EcalRecAlgos/interface/EcalSeverityLevelAlgo.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
#include "CondFormats/DataRecord/interface/EcalChannelStatusRcd.h"
#include "CondFormats/EcalObjects/interface/EcalChannelStatus.h"

// Data formats
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "DataFormats/BTauReco/interface/JetTag.h"
#include "DataFormats/JetReco/interface/JetID.h"
#include "DataFormats/JetReco/interface/JetCollection.h"

#include "DataFormats/METReco/interface/METCollection.h"
#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/METReco/interface/CaloMETCollection.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/GenMET.h"
#include "DataFormats/PatCandidates/interface/MET.h"


#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "CondFormats/L1TObjects/interface/L1GtTriggerMenu.h"
#include "CondFormats/DataRecord/interface/L1GtTriggerMenuRcd.h"

#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/Candidate/interface/Candidate.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"

#include "DataFormats/METReco/interface/BeamHaloSummary.h"

#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

/*
  #include "DataFormats/AnomalousEcalDataFormats/interface/AnomalousECALVariables.h"
  #include "PhysicsTools/EcalAnomalousEventFilter/interface/EcalBoundaryInfoCalculator.h"
*/

#include "MagneticField/Engine/interface/MagneticField.h"

#include "Geometry/Records/interface/CaloTopologyRecord.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h"
#include "RecoParticleFlow/PFClusterTools/interface/ClusterClusterMapping.h"

#include "DiLeptonAnalysis/NTupleProducer/interface/ETHVertexInfo.h"

// Interface
#include "DiLeptonAnalysis/NTupleProducer/interface/NTupleProducer.h"

namespace LHAPDF {
  void initPDFSet(const std::string& filename, int member=0);
  void initPDF(int member=0);
  int numberPDF(int nset);
  int numberPDF();
  void usePDFMember(int nset, int member);
  double xfx( double x, double Q, int fl);
  double getXmin(int nset, int member);
  double getXmax(int nset, int member);
  double getQ2min(int nset, int member);
  double getQ2max(int nset, int member);
  void extrapolate(bool extrapolate=true);
}

void FlipGenStoreFlag(int index, int Promptness[], int genMo1Index[], int genMo2Index[], bool StoreFlag[]) {
  if(StoreFlag[index]) return;//particle has already been marked. 
  StoreFlag[index]=true;
  if(genMo1Index[index]>0) FlipGenStoreFlag(genMo1Index[index], Promptness, genMo1Index, genMo2Index, StoreFlag);
  if(genMo2Index[index]>0) FlipGenStoreFlag(genMo2Index[index], Promptness, genMo1Index, genMo2Index, StoreFlag);
}
  

NTupleProducer::NTupleProducer(const edm::ParameterSet& iConfig){
  // Main settings
  fIsRealData = iConfig.getUntrackedParameter<bool>("isRealData");
  fIsModelScan = iConfig.getUntrackedParameter<bool>("isModelScan");
  if(fIsRealData&&fIsModelScan) fIsModelScan=false; // avoiding possible mistakes

  // InputTags
  fMuonTag            = iConfig.getUntrackedParameter<edm::InputTag>("tag_muons");
  fElectronTag        = iConfig.getUntrackedParameter<edm::InputTag>("tag_electrons");
  fEleIdWP            = iConfig.getUntrackedParameter<std::string>("tag_elidWP");
  fMuIsoDepTkTag      = iConfig.getUntrackedParameter<edm::InputTag>("tag_muisodeptk");
  fMuIsoDepECTag      = iConfig.getUntrackedParameter<edm::InputTag>("tag_muisodepec");
  fMuIsoDepHCTag      = iConfig.getUntrackedParameter<edm::InputTag>("tag_muisodephc");
  fJetTag             = iConfig.getUntrackedParameter<edm::InputTag>("tag_jets");
  fJetCorrs           = iConfig.getUntrackedParameter<std::string>("jetCorrs");
  fBtag1Tag           = iConfig.getUntrackedParameter<edm::InputTag>("tag_btag1");
  fBtag2Tag           = iConfig.getUntrackedParameter<edm::InputTag>("tag_btag2");
  fBtag3Tag           = iConfig.getUntrackedParameter<edm::InputTag>("tag_btag3");
  fBtag4Tag           = iConfig.getUntrackedParameter<edm::InputTag>("tag_btag4");
  fRawCaloMETTag      = iConfig.getUntrackedParameter<edm::InputTag>("tag_rawcalomet");
  fTCMETTag           = iConfig.getUntrackedParameter<edm::InputTag>("tag_tcmet");
  fPFMETTag           = iConfig.getUntrackedParameter<edm::InputTag>("tag_pfmet");
  fPFMETPATTag        = iConfig.getUntrackedParameter<edm::InputTag>("tag_pfmetPAT");
  fCorrCaloMETTag     = iConfig.getUntrackedParameter<edm::InputTag>("tag_corrcalomet");
  fGenMETTag          = iConfig.getUntrackedParameter<edm::InputTag>("tag_genmet");
  fVertexTag          = iConfig.getUntrackedParameter<edm::InputTag>("tag_vertex");
  fTrackTag           = iConfig.getUntrackedParameter<edm::InputTag>("tag_tracks");
  fPhotonTag          = iConfig.getUntrackedParameter<edm::InputTag>("tag_photons");
  fCalTowTag          = iConfig.getUntrackedParameter<edm::InputTag>("tag_caltow");
  fEBRecHitsTag       = iConfig.getUntrackedParameter<edm::InputTag>("tag_EBrechits");
  fEERecHitsTag       = iConfig.getUntrackedParameter<edm::InputTag>("tag_EErechits");
  fGenPartTag         = iConfig.getUntrackedParameter<edm::InputTag>("tag_genpart");
  fGenJetTag          = iConfig.getUntrackedParameter<edm::InputTag>("tag_genjets");
  fL1TriggerTag       = iConfig.getUntrackedParameter<edm::InputTag>("tag_l1trig");
  fHLTTrigEventTag    = iConfig.getUntrackedParameter<edm::InputTag>("tag_hlttrigevent");
  if(!fIsModelScan) fHBHENoiseResultTag    = iConfig.getUntrackedParameter<edm::InputTag>("tag_hcalnoise");
  if(!fIsModelScan) fHBHENoiseResultTagIso = iConfig.getUntrackedParameter<edm::InputTag>("tag_hcalnoiseIso");
  fSrcRho             = iConfig.getUntrackedParameter<edm::InputTag>("tag_srcRho");
  fSrcRhoPFnoPU       = iConfig.getUntrackedParameter<edm::InputTag>("tag_srcRhoPFnoPU");
  pfphotonsProducerTag = iConfig.getUntrackedParameter<edm::InputTag>("tag_pfphotonsProducer");
  pfProducerTag = iConfig.getUntrackedParameter<edm::InputTag>("tag_pfProducer");
  fSCTagBarrel = iConfig.getUntrackedParameter<edm::InputTag>("tag_SC_barrel");
  fSCTagEndcap = iConfig.getUntrackedParameter<edm::InputTag>("tag_SC_endcap");
  fTrackCollForVertexing = iConfig.getUntrackedParameter<edm::InputTag>("tag_fTrackCollForVertexing");
  fAllConversionsCollForVertexing = iConfig.getUntrackedParameter<edm::InputTag>("tag_fallConversionsCollForVertexing");
  perVtxMvaWeights = iConfig.getUntrackedParameter<std::string>("tag_perVtxMvaWeights");
  perVtxMvaMethod = iConfig.getUntrackedParameter<std::string>("tag_perVtxMvaMethod");
  perEvtMvaWeights = iConfig.getUntrackedParameter<std::string>("tag_perEvtMvaWeights");
  perEvtMvaMethod = iConfig.getUntrackedParameter<std::string>("tag_perEvtMvaMethod");

  doVertexingFlag = iConfig.getUntrackedParameter<bool>("tag_doVertexing");
  if (fIsModelScan) doVertexingFlag=false;

  // Event Selection
  fMinMuPt        = iConfig.getParameter<double>("sel_minmupt");
  fMaxMuEta       = iConfig.getParameter<double>("sel_maxmueta");
  fMinElPt        = iConfig.getParameter<double>("sel_minelpt");
  fMaxElEta       = iConfig.getParameter<double>("sel_maxeleta");
  fMinCorJPt      = iConfig.getParameter<double>("sel_mincorjpt");
  fMinRawJPt      = iConfig.getParameter<double>("sel_minrawjpt");
  fMaxJEta        = iConfig.getParameter<double>("sel_maxjeta");
  fMinJEMFrac     = iConfig.getParameter<double>("sel_minjemfrac");

  fMinTrkPt       = iConfig.getParameter<double>("sel_mintrkpt");
  fMaxTrkEta      = iConfig.getParameter<double>("sel_maxtrketa");
  fMaxTrkNChi2    = iConfig.getParameter<double>("sel_maxtrknchi2");
  fMinTrkNHits    = iConfig.getParameter<int>("sel_mintrknhits");

  fMinPhotonPt    = iConfig.getParameter<double>("sel_minphopt");
  fMaxPhotonEta   = iConfig.getParameter<double>("sel_maxphoeta");
  fMinSCraw       = iConfig.getParameter<double>("sel_minSCraw");
  fMinEBRechitE   = iConfig.getParameter<double>("sel_fminebrechitE");

  fMinGenLeptPt   = iConfig.getParameter<double>("sel_mingenleptpt");
  fMaxGenLeptEta  = iConfig.getParameter<double>("sel_maxgenlepteta");
  fMinGenPhotPt   = iConfig.getParameter<double>("sel_mingenphotpt");
  fMaxGenPhotEta  = iConfig.getParameter<double>("sel_maxgenphoteta");
  fMinGenJetPt    = iConfig.getParameter<double>("sel_mingenjetpt");
  fMaxGenJetEta   = iConfig.getParameter<double>("sel_maxgenjeteta");

  if(fIsModelScan) {
    LHAPDF::initPDFSet("cteq66.LHgrid",1);
    *fTNPdfs = LHAPDF::numberPDF();
  }

  CrackCorrFunc    = EcalClusterFunctionFactory::get()->create("EcalClusterCrackCorrection", iConfig);
  LocalCorrFunc    = EcalClusterFunctionFactory::get()->create("EcalClusterLocalContCorrection",iConfig);

  // Dump the full configuration
  edm::LogVerbatim("NTP") << "---------------------------------";
  edm::LogVerbatim("NTP") << " ==> NTupleProducer Constructor ...";
  edm::LogVerbatim("NTP") << iConfig;

  // Create additional jet fillers
  std::vector<edm::ParameterSet> jConfigs = iConfig.getParameter<std::vector<edm::ParameterSet> >("jets");
  for (size_t i=0; i<jConfigs.size(); ++i)
    if ( jConfigs[i].getUntrackedParameter<bool>("isPat") ) jetFillers.push_back( new JetFillerPat(jConfigs[i], fIsRealData) );
    else jetFillers.push_back( new JetFillerReco(jConfigs[i], fIsRealData) );

  // Create additional lepton fillers
  std::vector<edm::ParameterSet> lConfigs = iConfig.getParameter<std::vector<edm::ParameterSet> >("leptons");
  for (size_t i=0; i<lConfigs.size(); ++i) {
    std::string type(lConfigs[i].getUntrackedParameter<std::string>("type"));
    if ( type == "electron" ) 
      electronFillers.push_back( new PatElectronFiller(lConfigs[i], fIsRealData) );
    else if ( type == "muon" ) 
      muonFillers.push_back( new PatMuonFiller(lConfigs[i], fIsRealData) );
    else if ( type == "tau" ) 
      tauFillers.push_back( new PatTauFiller(lConfigs[i], fIsRealData) );
  }
        

  // Get list of trigger paths to store the triggering object info. of
  //std::vector<std::string> v(iConfig.getUntrackedParameter<std::vector<std::string> >("hlt_labels"));
  fRHLTLabels.reset( new std::vector<std::string>(iConfig.getUntrackedParameter<std::vector<std::string> >("hlt_labels")) );
  fTNpaths = fRHLTLabels->size();
	
  //OOT pu reweighting
  if(!fIsRealData){
    *fRPileUpData = iConfig.getUntrackedParameter<std::vector<std::string> >("pu_data");
    *fRPileUpMC   = iConfig.getUntrackedParameter<std::vector<std::string> >("pu_mc");
    if(! (*fRPileUpData)[0].empty() && !(*fRPileUpMC)[0].empty() ){
      LumiWeights_      = edm::LumiReWeighting((*fRPileUpMC)[0], (*fRPileUpData)[0], (*fRPileUpMC)[1], (*fRPileUpData)[1]);
    }
  }

  // Declare all products to be stored (needs to be done at construction time)
  declareProducts();
  std::vector<filler::PPair > list;
  typedef std::vector<filler::PPair>::iterator PPI;
  for ( std::vector<JetFillerBase*>::iterator it = jetFillers.begin();
        it != jetFillers.end(); ++it ) {
    list = (*it)->declareProducts();
    for ( PPI ip = list.begin(); ip != list.end(); ++ip )
      produces<edm::InEvent>( ip->first, ip->second );
  }
  for ( std::vector<PatMuonFiller*>::iterator it = muonFillers.begin(); 
        it != muonFillers.end(); ++it ) {
    list = (*it)->declareProducts();
    for ( PPI ip = list.begin(); ip != list.end(); ++ip )
      produces<edm::InEvent>( ip->first, ip->second );
  }
  for ( std::vector<PatElectronFiller*>::iterator it = electronFillers.begin(); 
        it != electronFillers.end(); ++it ) {
    list = (*it)->declareProducts();
    for (PPI ip = list.begin(); ip != list.end(); ++ip)
      produces<edm::InEvent>( ip->first, ip->second );
  }
  for ( std::vector<PatTauFiller*>::iterator it = tauFillers.begin(); 
        it != tauFillers.end(); ++it ) {
    list = (*it)->declareProducts();
    for ( PPI ip = list.begin();  ip != list.end(); ++ip )
      produces<edm::InEvent>( ip->first, ip->second );
  }

}

//________________________________________________________________________________________
// Method called once per job
void NTupleProducer::beginJob(void) {

  fNTotEvents = 0;
  fNFillTree  = 0;
  fFirstevent = true;

}


//________________________________________________________________________________________
// Method called once for each event
bool NTupleProducer::filter(edm::Event& iEvent, const edm::EventSetup& iSetup){

  ++fNTotEvents;

  using namespace edm;
  using namespace std;
  using namespace reco;
  using reco::MuonCollection;
  using reco::JetTagCollection;

  // Reset all the variables
  resetProducts();
  for ( std::vector<JetFillerBase*>::iterator it = jetFillers.begin(); 
        it != jetFillers.end(); ++it ) 
    (*it)->resetProducts();
  for ( std::vector<PatMuonFiller*>::iterator it = muonFillers.begin(); 
        it != muonFillers.end(); ++it ) 
    (*it)->resetProducts();
  for ( std::vector<PatElectronFiller*>::iterator it = electronFillers.begin(); 
        it != electronFillers.end(); ++it ) 
    (*it)->resetProducts();
  for ( std::vector<PatTauFiller*>::iterator it = tauFillers.begin(); 
        it != tauFillers.end(); ++it ) 
    (*it)->resetProducts();


  ////////////////////////////////////////////////////////////////////////////////
  // Get the collections /////////////////////////////////////////////////////////
  Handle<View<Muon> > muons;
  iEvent.getByLabel(fMuonTag,muons); // 'muons'

  Handle<View<GsfElectron> > electrons;
  iEvent.getByLabel(fElectronTag, electrons); // 'gsfElectrons'

  // Jets and Jet Correctors
  Handle<View<Jet> > jets;
  iEvent.getByLabel(fJetTag,jets);

  // rho for L1FastJet
  edm::Handle<double> rho;
  iEvent.getByLabel(fSrcRho,rho);
  *fTRho = *rho;
	
  // rho for L1FastJet running PFnoPU
  edm::Handle<double> rhoNoPU;
  iEvent.getByLabel(fSrcRhoPFnoPU,rhoNoPU);
  *fTRhoPFnoPU = *rhoNoPU;

  // beam halo
  edm::Handle<BeamHaloSummary> TheBeamHaloSummary;
  iEvent.getByLabel("BeamHaloSummary",TheBeamHaloSummary);
  const BeamHaloSummary TheSummary = (*TheBeamHaloSummary.product());
  *fTCSCTightHaloID = static_cast<int>  (TheSummary.CSCTightHaloId());


  // collect information for b-tagging (4 tags)
  Handle<JetTagCollection> jetsAndProbsTkCntHighEff;
  iEvent.getByLabel(fBtag1Tag,jetsAndProbsTkCntHighEff);

  Handle<JetTagCollection> jetsAndProbsTkCntHighPur;
  iEvent.getByLabel(fBtag2Tag,jetsAndProbsTkCntHighPur);

  Handle<JetTagCollection> jetsAndProbsSimpSVHighEff;
  iEvent.getByLabel(fBtag3Tag,jetsAndProbsSimpSVHighEff);

  Handle<JetTagCollection> jetsAndProbsSimpSVHighPur;
  iEvent.getByLabel(fBtag4Tag,jetsAndProbsSimpSVHighPur);

  //Get Tracks collection
  Handle<TrackCollection> tracks;
  iEvent.getByLabel(fTrackTag, tracks);

  //Get Photon collection
  Handle<View<Photon> > photons;
  iEvent.getByLabel(fPhotonTag, photons);

  //Get SC collections
  Handle<SuperClusterCollection> BarrelSuperClusters;
  Handle<SuperClusterCollection> EndcapSuperClusters;
  iEvent.getByLabel(fSCTagBarrel,BarrelSuperClusters);
  iEvent.getByLabel(fSCTagEndcap,EndcapSuperClusters);
	
  //PFcandidates
  edm::Handle<reco::PFCandidateCollection> pfCandidates;
  iEvent.getByLabel(pfProducerTag, pfCandidates);
	
  //Electron collection
  edm::Handle<reco::GsfElectronCollection> electronHandle;
  iEvent.getByLabel(fElectronTag, electronHandle);
	
  //PF Photon collection
  edm::Handle<reco::PhotonCollection> pfPhotonHandle;
  iEvent.getByLabel(pfphotonsProducerTag,pfPhotonHandle);
	
  edm::Handle<reco::VertexCollection> alternativeVertexHandle;
  iEvent.getByLabel(fVertexTag, alternativeVertexHandle);
	
  // MET
  Handle<CaloMETCollection> calomet;
  iEvent.getByLabel(fRawCaloMETTag, calomet);

  Handle<METCollection> tcmet;
  iEvent.getByLabel(fTCMETTag, tcmet);

  Handle<View<PFMET> > pfmet;
  iEvent.getByLabel(fPFMETTag, pfmet);
	
  Handle<View<pat::MET> > pfMETpat;
  iEvent.getByLabel(fPFMETPATTag,pfMETpat);  //'pfMET PAT'


  Handle<CaloMETCollection> corrmujesmet;
  iEvent.getByLabel(fCorrCaloMETTag, corrmujesmet);

  // Get beamspot for d0 determination
  BeamSpot beamSpot;
  Handle<BeamSpot> beamSpotHandle;
  iEvent.getByLabel("offlineBeamSpot", beamSpotHandle);
  beamSpot = *beamSpotHandle;

  // Primary vertex
  edm::Handle<VertexCollection> vertices;
  iEvent.getByLabel(fVertexTag, vertices);
  const reco::Vertex *primVtx = &(*(vertices.product()))[0]; // Just take first vertex ...

  // Get Muon IsoDeposits
  // ECAL:
  edm::Handle<edm::ValueMap<reco::IsoDeposit> > IsoDepECValueMap;
  iEvent.getByLabel(fMuIsoDepECTag, IsoDepECValueMap);
  const edm::ValueMap<reco::IsoDeposit> &ECDepMap = *IsoDepECValueMap.product();
  // HCAL:
  edm::Handle<edm::ValueMap<reco::IsoDeposit> > IsoDepHCValueMap;
  iEvent.getByLabel(fMuIsoDepHCTag, IsoDepHCValueMap);
  const edm::ValueMap<reco::IsoDeposit> &HCDepMap = *IsoDepHCValueMap.product();

  // Get CaloTowers
  edm::Handle<CaloTowerCollection> calotowers;
  iEvent.getByLabel(fCalTowTag, calotowers);

  // For ECAL cleaning: rechit and channel status
  edm::Handle<EcalRecHitCollection> ebRecHits;
  edm::Handle<EcalRecHitCollection> eeRecHits;
  iEvent.getByLabel(fEBRecHitsTag,ebRecHits);
  iEvent.getByLabel(fEERecHitsTag,eeRecHits);
  // edm::ESHandle<EcalChannelStatus> chStatus;
  // iSetup.get<EcalChannelStatusRcd>().get(chStatus);
  // const EcalChannelStatus * channelStatus = chStatus.product();

  edm::ESHandle<CaloGeometry> geometry ;
  iSetup.get<CaloGeometryRecord>().get(geometry);

  edm::ESHandle<CaloTopology> theCaloTopo;
  iSetup.get<CaloTopologyRecord>().get(theCaloTopo);
  const CaloTopology *topology = theCaloTopo.product();

  CrackCorrFunc->init(iSetup);
  LocalCorrFunc->init(iSetup);

  // ECAL dead cell Trigger Primitive filter
  edm::Handle<bool> EcalDeadTPFilterFlag;
  iEvent.getByLabel("ecalDeadCellTPfilter",EcalDeadTPFilterFlag);
  *fTEcalDeadTPFilterFlag = (*EcalDeadTPFilterFlag) ? 1 : 0;
	
  // Stevens recovRecHitFilter
  // http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/UserCode/lowette/SandBox/Skims/python/recovRecHitFilter_cfi.py?sortby=date&view=log
  // https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFilters#RecovRecHitFilter
  edm::Handle<bool> RecovRecHitFilterFlag;
  iEvent.getByLabel("recovRecHitFilter","Result",RecovRecHitFilterFlag);
  *fTRecovRecHitFilterFlag = (*RecovRecHitFilterFlag) ? 1 : 0;

  // RA2 tracking tailure filter
  // https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFilters
  // http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/UserCode/seema/SandBox/Skims/python/trackingFailureFilter_cfi.py?hideattic=0&revision=1.3&view=markup&pathrev=MAIN
  edm::Handle<bool> RA2TrackingFailureFlag;
  iEvent.getByLabel("trackingFailureFilter",RA2TrackingFailureFlag);
  *fTRA2TrackingFailureFilterFlag = (int) *RA2TrackingFailureFlag;

  // Colin's PBNR filter
  edm::Handle<bool> ParticleBasedNoiseRejectionFlag;
  //FR: this requires modifications in the CVS code. Disabled.
  //iEvent.getByLabel("jetIDFailure",ParticleBasedNoiseRejectionFlag);
  //fPBNRFlag = (int) *ParticleBasedNoiseRejectionFlag;

  /*
  // TEMPORARILY DISABLED FOR RUNNING ON CMSSW_3_9_X

  // Ecal dead cells: boundary energy check:
  // fTEcalDeadCellBEFlag ==0 if >24 dead cells with >10 GeV boundary energy
  edm::InputTag ecalAnomalousFilterTag("EcalAnomalousEventFilter","anomalousECALVariables");
  Handle<AnomalousECALVariables> anomalousECALvarsHandle;
  iEvent.getByLabel(ecalAnomalousFilterTag, anomalousECALvarsHandle);
  AnomalousECALVariables anomalousECALvars;
  if (anomalousECALvarsHandle.isValid()) {
  anomalousECALvars = *anomalousECALvarsHandle;
  } else {
  edm::LogError("NTP") << " anomalous ECAL Vars not valid/found ";
  }
  // default call to isEcalNoise() returns true if
  //      1. highestEnergyDepositAroundDeadCell > 10 (maxBoundaryEnergy)
  //   && 2. DeadClusterSize > 24 (minDeadClusterSize)
  //   to change default values: call isEcalNoise(minDeadClusterSize, maxBoundaryEnergy)
  if(anomalousECALvars.isEcalNoise()==1){
  fTEcalDeadCellBEFlag = 0;
  }else{
  fTEcalDeadCellBEFlag = 1;
  }
  // ------------------------------
  // ECAL GAP energy
  // minimal BE spefified as cutBoundEnergyGapEE/cutBoundEnergyGapEB in cfi.py
  fTnECALGapClusters=0;
  // EB
  for (int i = 0; i < (int) anomalousECALvars.v_enNeighboursGap_EB.size(); ++i) {
  BoundaryInformation bInfo =anomalousECALvars.v_enNeighboursGap_EB[i];
  fTEcalGapBE[fTnECALGapClusters]           = bInfo.boundaryEnergy;
  fTEcalGapClusterSize[fTnECALGapClusters]  = bInfo.neighboursWithSameFlag.size()+1;
  fTnECALGapClusters++;
  if(fTnECALGapClusters >= gMaxnECALGapClusters){
  edm::LogWarning("NTP") << "@SUB=analyze()"
  << "More than " << static_cast<int>(gMaxnECALGapClusters)
  << " ECAL GAP Clusters!!";
  break;
  }
  }
  // EE
  for (int i = 0; i < (int) anomalousECALvars.v_enNeighboursGap_EE.size(); ++i) {
  BoundaryInformation bInfo = anomalousECALvars.v_enNeighboursGap_EE[i];
  fTEcalGapBE[fTnECALGapClusters]           = bInfo.boundaryEnergy;
  fTEcalGapClusterSize[fTnECALGapClusters]  = bInfo.neighboursWithSameFlag.size()+1;
  fTnECALGapClusters++;
  if(fTnECALGapClusters >= gMaxnECALGapClusters){
  edm::LogWarning("NTP") << "@SUB=analyze()"
  << "More than " << static_cast<int>(gMaxnECALGapClusters)
  << " ECAL GAP Clusters!!";
  break;
  }
  }
  */

  // Retrieve HB/HE noise flag
	
  if(!fIsModelScan) {
    edm::Handle<bool> hbHeNoiseFlag;
    iEvent.getByLabel(fHBHENoiseResultTag,hbHeNoiseFlag);
    *fTHBHENoiseFlag    = static_cast<int>(*hbHeNoiseFlag);
	   
    edm::Handle<bool> hbHeNoiseFlagIso;
    iEvent.getByLabel(fHBHENoiseResultTagIso,hbHeNoiseFlagIso);
    *fTHBHENoiseFlagIso = static_cast<int>(*hbHeNoiseFlagIso);
  }

  // Get Transient Track Builder
  ESHandle<TransientTrackBuilder> theB;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB);

  // Get GenEventInfoProduct
  edm::Handle<GenEventInfoProduct> genEvtInfo;
  edm::Handle<std::vector<PileupSummaryInfo> > pileupInfo;

	
  double MyWeightTotal   =-999;
  double MyWeightInTime  =-999;
  // If need to run on MC prior to Spring11 (CMSSW_3_9_X series) change code accordingly here:
  // edm::Handle<PileupSummaryInfo> pileupInfo; 
  if(!fIsRealData){
    // Get LHEEventProduct with partonic momenta. 	
    Handle<LHEEventProduct> evt;
    bool LHEEventProduct_found= iEvent.getByType( evt );
    if(LHEEventProduct_found){ 
      const lhef::HEPEUP hepeup_ = evt->hepeup();
      const std::vector<lhef::HEPEUP::FiveVector> pup_ = hepeup_.PUP; // px, py, pz, E, M
      float partonicHT=0;
      for(unsigned int i=0; i<pup_.size(); ++i){
        if(hepeup_.ISTUP[i]!=1) continue; // quarks and gluons come out of MG/ME as stable
        int idabs=abs(hepeup_.IDUP[i]); 
        if( idabs != 21 && (idabs<1 || idabs>6) ) continue;  // gluons and quarks
        float ptPart = sqrt( pow(hepeup_.PUP[i][0],2) + pow(hepeup_.PUP[i][1],2) ); // first entry is px, second py
        partonicHT +=ptPart; // sum QCD partonic HT
      }
      *fTQCDPartonicHT = partonicHT;
    }

    iEvent.getByLabel("generator", genEvtInfo);
    *fTPtHat       = genEvtInfo->hasBinningValues() ? (genEvtInfo->binningValues())[0] : 0.0;
    *fTSigProcID   = genEvtInfo->signalProcessID();
    *fTPDFScalePDF = genEvtInfo->pdf()->scalePDF;
    *fTPDFID1      = genEvtInfo->pdf()->id.first;
    *fTPDFID2      = genEvtInfo->pdf()->id.second;
    *fTPDFx1       = genEvtInfo->pdf()->x.first;
    *fTPDFx2       = genEvtInfo->pdf()->x.second;
    *fTPDFxPDF1    = genEvtInfo->pdf()->xPDF.first;
    *fTPDFxPDF2    = genEvtInfo->pdf()->xPDF.second;
    *fTGenWeight   = genEvtInfo->weight();
     
    iEvent.getByLabel("addPileupInfo", pileupInfo);
    std::vector<PileupSummaryInfo>::const_iterator PVI;

    for (PVI = pileupInfo->begin(); PVI !=pileupInfo->end(); ++PVI) {
      if( PVI->getBunchCrossing() == 0 ) { // in-time PU
        *fTPUnumInteractions     = PVI->getPU_NumInteractions();
        *fTPUnumTrueInteractions = PVI->getTrueNumInteractions();
		    
        if(*fTPUnumInteractions > gMaxNPileup){
          edm::LogWarning("NTP") << "@SUB=analyze()"
                                 << "More than " << static_cast<int>(gMaxNPileup)
                                 << " generated Pileup events found, increase size!";
          *fTGoodEvent = 1;
        }
		    
        *fTPUnumFilled = (int)PVI->getPU_zpositions().size();
        for( int i = 0; i < *fTPUnumFilled; i++) {
          if(i >= gMaxNPileup) break; // hard protection
          fTPUzPositions  ->push_back( PVI->getPU_zpositions()[i]   );
          fTPUsumPtLowPt  ->push_back( PVI->getPU_sumpT_lowpT()[i]  );
          fTPUsumPtHighPt ->push_back( PVI->getPU_sumpT_highpT()[i] );
          fTPUnTrksLowPt  ->push_back( PVI->getPU_ntrks_lowpT()[i]  );
          fTPUnTrksHighPt ->push_back( PVI->getPU_ntrks_highpT()[i] );
        }
      } else if ( PVI->getBunchCrossing() == 1 ) { // OOT pile-Up: this is the 50ns late Bunch
        *fTPUOOTnumInteractionsLate = PVI->getPU_NumInteractions();
      } else if ( PVI->getBunchCrossing() == -1 ) { // OOT pile-Up: this is the 50ns early Bunch
        *fTPUOOTnumInteractionsEarly = PVI->getPU_NumInteractions();
      }
    }
    //see https://twiki.cern.ch/twiki/bin/view/CMS/PileupMCReweightingUtilities 
    // as well as http://cmslxr.fnal.gov/lxr/source/PhysicsTools/Utilities/src/LumiReWeighting.cc
    if(!(*fRPileUpData)[0].empty() && !(*fRPileUpMC)[0].empty() ){
      //const EventBase* iEventB = dynamic_cast<const EventBase*>(&iEvent);
      MyWeightTotal  = LumiWeights_.weightOOT( iEvent ); // this is the total weight inTimeWeight * WeightOOTPU * Correct_Weights2011
      MyWeightInTime = LumiWeights_.weight   ( iEvent ); // this is the inTimeWeight only
    }
    if(fIsModelScan) {
      edm::Handle<GenEventInfoProduct> pdfstuff;
      if (!iEvent.getByLabel("generator", pdfstuff)) {
        edm::LogError("PDFWeightProducer") << ">>> PdfInfo not found !!!";
        return false;
      }
		
      float Q = pdfstuff->pdf()->scalePDF;
      int id1 = pdfstuff->pdf()->id.first;
      double x1 = pdfstuff->pdf()->x.first;
      //       double pdf1 = pdfstuff->pdf()->xPDF.first;
      int id2 = pdfstuff->pdf()->id.second;
      double x2 = pdfstuff->pdf()->x.second;
      //       double pdf2 = pdfstuff->pdf()->xPDF.second;

      fTpdfW->push_back(1);
      LHAPDF::initPDF(0);

      double newpdf1_0 = LHAPDF::xfx(x1, Q, id1)/x1;
      double newpdf2_0 = LHAPDF::xfx(x2, Q, id2)/x2;
			
      float pdfWsum=0;
      for(int pdf=1; pdf <= *fTNPdfs; pdf++){
        LHAPDF::initPDF(pdf);
        double newpdf1 = LHAPDF::xfx(x1, Q, id1)/x1;
        double newpdf2 = LHAPDF::xfx(x2, Q, id2)/x2;
        fTpdfW->push_back( newpdf1/newpdf1_0*newpdf2/newpdf2_0 );
        pdfWsum += fTpdfW->back();
      }
      *fTpdfWsum = pdfWsum;

      int process = 0;
      if (*fTSigProcID>=237 && *fTSigProcID<=242) process=1;//"ng";
      else if(*fTSigProcID>=246 && *fTSigProcID<=256) process=2;//"ns";
      else if(*fTSigProcID>=216 && *fTSigProcID<=236) process=3;//"nn";
      else if(*fTSigProcID>=201 && *fTSigProcID<=214) process=4;//"ll";
      else if( (*fTSigProcID>=274 && *fTSigProcID<=280) || (*fTSigProcID>=284 && *fTSigProcID<=286)) process=5;//"sb";
      else if( (*fTSigProcID>=271 && *fTSigProcID<=273) || (*fTSigProcID>=281 && *fTSigProcID<=283) 
               || (*fTSigProcID>=291 && *fTSigProcID<=293)) process=6;//"ss";
      else if(*fTSigProcID>=261 && *fTSigProcID<=265) process=7;//"tb";
      else if( (*fTSigProcID>=287 && *fTSigProcID<=290) || *fTSigProcID==296) process=8;//"bb";
      else if(*fTSigProcID>=243 && *fTSigProcID<=244) process=9;//"gg";
      else if( (*fTSigProcID>=258 && *fTSigProcID<=259) || (*fTSigProcID>=294 && *fTSigProcID<=295) ) process=10;//"sg";
      *fTprocess = process;
    }
  }
  *fTPUWeightTotal  = MyWeightTotal;
  *fTPUWeightInTime = MyWeightInTime;
	
	
  //////////////////////////////////////////////////////////////////////////////
  // Trigger information
  Handle<L1GlobalTriggerReadoutRecord> l1GtReadoutRecord;
  iEvent.getByLabel(fL1TriggerTag, l1GtReadoutRecord);

  Handle<trigger::TriggerEvent> triggerEventHLT;
  iEvent.getByLabel("hltTriggerSummaryAOD", triggerEventHLT);
  
  // Retrieve trigger results, with process name auto-discovered in beginRun()
  Handle<TriggerResults> triggers;
  iEvent.getByLabel(InputTag("TriggerResults","",fHltConfig.processName()), triggers);
  const TriggerResults& tr = *triggers;

  // Get trigger results and prescale
  for(unsigned int i = 0; i < tr.size(); i++ ){
    bool fired = tr[i].accept();
    fTHLTResults->push_back( fired ? 1:0 );
    fTHLTPrescale->push_back( fHltConfig.prescaleValue(iEvent, iSetup, (*fRHLTNames)[i]) );
    //FIXME: TOO MANY ERRORS FOR THE MOMENT.
    //           // Check that there is only 1 L1 seed (otherwise we can't get the combined prescale)
    //           if ( fHltConfig.hltL1GTSeeds(fTHLTmenu[i]).size() == 1 ) {
    //             const std::string l1tname(fHltConfig.hltL1GTSeeds(fTHLTmenu[i]).at(0).second);
    //             if ( l1tname.find(" OR ") == string::npos ) { // Poor man's check...
    //               std::pair<int,int> ps = fHltConfig.prescaleValues(iEvent, iSetup, fTHLTmenu[i]);
    //               fTHLTprescale[i] = ps.first*ps.second; // Multiply L1 and HLT prescales
    //             }
    //           }
  }
  for( unsigned int i = 0; i < gMaxL1PhysBits; ++i ){
    bool fired = l1GtReadoutRecord->decisionWord()[i];
    fTL1PhysResults->push_back( fired ? 1:0 );
  }
  for( unsigned int i = 0; i < gMaxL1TechBits; ++i){
    bool fired = l1GtReadoutRecord->technicalTriggerWord()[i];
    fTL1TechResults->push_back( fired ? 1:0 );
  }

  // Store information for some trigger paths
  edm::Handle<trigger::TriggerEvent> trgEvent;
  iEvent.getByLabel(fHLTTrigEventTag, trgEvent);

  InputTag collectionTag;
  // Loop over path names and get related objects
  for (size_t i=0; i<fTNpaths; ++i) {
    collectionTag = edm::InputTag((*fRHLTLabels)[i],"",fHltConfig.processName());
    size_t  filterIndex_ = trgEvent->filterIndex(collectionTag);
    if (filterIndex_<trgEvent->sizeFilters()) {
      const trigger::TriggerObjectCollection& TOC(trgEvent->getObjects());
      const trigger::Keys& keys = trgEvent->filterKeys(filterIndex_);
      // Loop over objects
      for ( size_t hlto = 0; hlto<keys.size(); ++hlto ) {
        if (hlto>=gMaxHltNObjs) {
          edm::LogWarning("NTP") << "@SUB=analyze()"
                                 << "Maximum number of triggering objects exceeded"
                                 << " for filter " << (*fRHLTLabels)[i];
          break;
        }
        // Update number of objects stored
        if ( *fTNHLTObjs<(int)hlto ) *fTNHLTObjs = hlto+1; // Not an index...
        const trigger::TriggerObject& TO(TOC[keys[hlto]]);
        fTHLTObjectID[i]->push_back( TO.id() );
        fTHLTObjectPt[i]->push_back( TO.pt() );
        fTHLTObjectEta[i]->push_back( TO.eta() );
        fTHLTObjectPhi[i]->push_back( TO.phi() );
      }
    }
  }

  ////////////////////////////////////////////////////////////////////////////////
  // Dump tree variables /////////////////////////////////////////////////////////
  *fTRun   = iEvent.id().run();
  *fTEvent = iEvent.id().event();
  *fTLumiSection = iEvent.luminosityBlock();

  *fTWeight = 1.0; // To be filled at some point?

  // Save position of primary vertex
  *fTPrimVtxx   = primVtx->x();
  *fTPrimVtxy   = primVtx->y();
  *fTPrimVtxz   = primVtx->z();
  *fTPrimVtxRho = primVtx->position().rho();
  *fTPrimVtxxE  = primVtx->xError();
  *fTPrimVtxyE  = primVtx->yError();
  *fTPrimVtxzE  = primVtx->zError();
  *fTPrimVtxNChi2 = primVtx->normalizedChi2();
  *fTPrimVtxNdof  = primVtx->ndof();
  *fTPrimVtxIsFake= primVtx->isFake();

  *fTPrimVtxPtSum = 0.;
  for(std::vector<TrackBaseRef>::const_iterator trackit = primVtx->tracks_begin(); trackit != primVtx->tracks_end(); ++trackit){
    *fTPrimVtxPtSum += (*trackit)->pt();
  }
  *fTPrimVtxGood = 0;

  // get all vertices
  int countVrtx = -1;
  for(VertexCollection::const_iterator vertexit = vertices->begin(); vertexit != vertices->end(); ++vertexit) {
    countVrtx++;
    if(countVrtx >= gMaxNVrtx){
      edm::LogWarning("NTP") << "@SUB=analyze()"
                             << "Maximum number of vertices exceeded";
      *fTGoodEvent = 1;
      *fTMaxVerticesExceed = 1;
      break;
    }

    fTVrtxX     ->push_back( vertexit->x() );
    fTVrtxY     ->push_back( vertexit->y() );
    fTVrtxZ     ->push_back( vertexit->z() );
    fTVrtxXE    ->push_back( vertexit->xError() );
    fTVrtxYE    ->push_back( vertexit->yError() );
    fTVrtxZE    ->push_back( vertexit->zError() );
    fTVrtxNdof  ->push_back( vertexit->ndof() );
    fTVrtxChi2  ->push_back( vertexit->normalizedChi2() );
    fTVrtxNtrks ->push_back( vertexit->nTracks() );
    fTVrtxSumPt ->push_back( vertexit->p4().pt() );
    fTVrtxIsFake->push_back( vertexit->isFake() );
  }
  *fTNVrtx = countVrtx+1;


  // Save position of beamspot
  *fTBeamspotx = (beamSpot.position()).x();
  *fTBeamspoty = (beamSpot.position()).y();
  *fTBeamspotz = (beamSpot.position()).z();

  IndexByPt indexComparator; // Need this to sort collections

  /////////////////////////////////////////
  /// GenVertices 
  if (!fIsRealData && doVertexingFlag){

    edm::Handle<reco::GenParticleCollection> gpH;
    iEvent.getByLabel(fGenPartTag, gpH);

    const float lowPtThrGenVtx = 0.1;
    const float highPtThrGenVtx = 0.5;

    *fTNgv = 0;
    for(reco::GenParticleCollection::const_iterator it_gen = gpH->begin(); 
        it_gen!= gpH->end(); ++it_gen){   

      if (*fTNgv>=gMaxNGenVtx){
        edm::LogWarning("NTP") << "@SUB=analyze"
                               << "Maximum number of gen-vertices exceeded..";
        *fTGoodEvent = 1;
        break;
      }

      if( it_gen->status() != 3 || !(it_gen->vx()!=0. || it_gen->vy()!=0. || it_gen->vz()!=0.)  ) continue; 

      // check for duplicate vertex
      bool duplicate = false;
      for(Int_t itv = 0; itv < *fTNgv; itv++) {
        TVector3 checkVtx = gv_pos[itv];
        if( (fabs(it_gen->vx()-checkVtx.X())<1e-5) &&  (fabs(it_gen->vy()-checkVtx.Y())<1e-5) && (fabs(it_gen->vz()-checkVtx.Z())<1e-5)) {
          duplicate = true;
          break;
        }
      }

      if (duplicate) continue;
    
      gv_pos[*fTNgv].SetXYZ(it_gen->vx(), it_gen->vy(), it_gen->vz());
    
      TVector3  this_gv_pos = gv_pos[*fTNgv];
      TVector3 p3(0,0,0);
    
      fTgvSumPtLo->push_back(0);
      fTgvNTkLo->push_back(0);
      fTgvSumPtHi->push_back(0);
      fTgvNTkHi->push_back(0);

      for(reco::GenParticleCollection::const_iterator part = gpH->begin(); 
          part!= gpH->end(); ++part){   
        if (part->pt()==0) continue;
        if( part->status() == 1 && part->charge() != 0 && fabs(part->eta())<2.5 &&
            ( fabs(part->vx()-this_gv_pos.X())<1.e-5 && fabs(part->vy()-this_gv_pos.Y())<1.e-5 && fabs(part->vz()-this_gv_pos.Z())<1.e-5 ) )  {
	
          TVector3 m(part->px(),part->py(),part->pz());
          p3 += m;
          if( m.Pt() > lowPtThrGenVtx ) {
            (*fTgvSumPtLo)[*fTNgv] += m.Pt();
            (*fTgvNTkLo)[*fTNgv] += 1;
            if( m.Pt() > highPtThrGenVtx ) {
              (*fTgvSumPtHi)[*fTNgv] += m.Pt();
              (*fTgvNTkHi)[*fTNgv] += 1;
            }
          }
        }
      }

      gv_p3[*fTNgv].SetXYZ(p3.X(),p3.Y(),p3.Z());

      (*fTNgv)++;
    }

  } // end gen vertices



  ////////////////////////////////////////////////////////////////////////////////
  // Get GenLeptons (+ Mother and GMother)
  if(!fIsRealData){
    edm::Handle<GenParticleCollection> gen;
    iEvent.getByLabel(fGenPartTag, gen);
    GenParticleCollection::const_iterator g_part;
    GenParticleCollection::const_iterator g_end = gen->end();

    std::vector<const GenParticle*> gen_lepts;
    std::vector<const GenParticle*> gen_moms;
    std::vector<const GenParticle*> gen_gmoms;


    // loop over genparticles to get gen_els and gen_mus
    for( g_part = gen->begin(); g_part != g_end; ++g_part ){
      // select stable leptons
      if( abs(g_part->pdgId()) !=  5  // b jets
          && abs(g_part->pdgId()) != 11
          && abs(g_part->pdgId()) != 12
          && abs(g_part->pdgId()) != 13
          && abs(g_part->pdgId()) != 14
          && abs(g_part->pdgId()) != 15
          && abs(g_part->pdgId()) != 16 ) continue;

      if(!( g_part->status() ==1 || (g_part->status() ==2 && (abs(g_part->pdgId())==5 || abs(g_part->pdgId())==15)))) continue;

      bool GenMomExists  (true);
      bool GenGrMomExists(true);

      if( g_part->pt()        < fMinGenLeptPt )  continue;
      if( fabs(g_part->eta()) > fMaxGenLeptEta ) continue;

      int gen_id= g_part->pdgId();
      const GenParticle* gen_lept = &(*g_part);

      // get mother of gen_lept
      const GenParticle* gen_mom = static_cast<const GenParticle*> (gen_lept->mother());
      if(gen_mom==NULL){
        edm::LogWarning("NTP") << "@SUB=analyze"
                               << " WARNING: GenParticle does not have a mother ";
        GenMomExists=false;
      }
      int m_id=-999;
      if(GenMomExists) m_id = gen_mom -> pdgId();

      if(m_id != gen_id || !GenMomExists);
      else{
        int id= m_id;
        while(id == gen_id && GenMomExists){
          gen_mom = static_cast<const GenParticle*> (gen_mom->mother());
          if(gen_mom==NULL){
            edm::LogWarning("NTP") << "@SUB=analyze"
                                   << " WARNING: GenParticle does not have a mother ";
            GenMomExists=false;
          }
          if(GenMomExists) id=gen_mom->pdgId();
        }
      }
      if(GenMomExists) m_id = gen_mom->pdgId();

      // get grand mother of gen_lept
      const GenParticle* gen_gmom=NULL;
      if(GenMomExists) gen_gmom  = static_cast<const GenParticle*>(gen_mom->mother());
      if(gen_gmom==NULL){
        edm::LogWarning("NTP") << "@SUB=analyze"
                               << " WARNING: GenParticle (" << m_id << ") does not have a GrandMother ";
        GenGrMomExists=false;
      }
      int gm_id=-999;
      if(GenGrMomExists) gm_id = gen_gmom->pdgId();
      if (m_id != gm_id || !GenGrMomExists);
      else{
        int id=gm_id;
        while(id == m_id && GenGrMomExists){
          gen_gmom  = static_cast<const GenParticle*>(gen_gmom->mother());
          if(gen_gmom==NULL){
            edm::LogWarning("NTP") << "@SUB=analyze"
                                   << " WARNING: GenParticle does not have a GrandMother ";
            GenGrMomExists=false;
          }
          if(GenGrMomExists) id = gen_gmom->pdgId();
        }
      }

      gen_lepts.push_back(gen_lept);
      gen_moms.push_back(gen_mom);
      gen_gmoms.push_back(gen_gmom);
    }

    // set variables
    if(gen_lepts.size()!=gen_moms.size() || gen_lepts.size()!=gen_gmoms.size()){
      edm::LogWarning("NTP") << "@SUB=analyze"
                             << "ERROR in filling of GenLeptons!! ";
    }else{
      *fTNGenLeptons = gen_lepts.size();
      for(int i=0; i<*fTNGenLeptons; ++i){
        if( i >= gMaxNGenLept) {
          edm::LogWarning("NTP") << "@SUB=analyze"
                                 << "Maximum number of gen-leptons exceeded..";
          *fTMaxGenLepExceed = 1;
          *fTGoodEvent = 1;
          break;
        }

        fTGenLeptonID->push_back( gen_lepts[i]->pdgId() );
        fTGenLeptonPt->push_back( gen_lepts[i]->pt() );
        fTGenLeptonEta->push_back( gen_lepts[i]->eta() );
        fTGenLeptonPhi->push_back( gen_lepts[i]->phi() );

        fTGenLeptonMID->push_back( (gen_moms[i]!=NULL )  ? gen_moms[i]->pdgId()  : -999 );
        fTGenLeptonMStatus->push_back( (gen_moms[i]!=NULL )  ? gen_moms[i]->status() : -999 );
        fTGenLeptonMPt->push_back( (gen_moms[i]!=NULL )  ? gen_moms[i]->pt()     : -999 );
        fTGenLeptonMEta->push_back( (gen_moms[i]!=NULL )  ? gen_moms[i]->eta()    : -999 );
        fTGenLeptonMPhi->push_back( (gen_moms[i]!=NULL )  ? gen_moms[i]->phi()    : -999 );

        fTGenLeptonGMID->push_back( (gen_gmoms[i]!=NULL ) ? gen_gmoms[i]->pdgId() : -999 );
        fTGenLeptonGMStatus->push_back( (gen_gmoms[i]!=NULL ) ? gen_gmoms[i]->status(): -999 );
        fTGenLeptonGMPt->push_back( (gen_gmoms[i]!=NULL ) ? gen_gmoms[i]->pt()    : -999 );
        fTGenLeptonGMEta->push_back( (gen_gmoms[i]!=NULL ) ? gen_gmoms[i]->eta()   : -999 );
        fTGenLeptonGMPhi->push_back( (gen_gmoms[i]!=NULL ) ? gen_gmoms[i]->phi()   : -999 );
      }
    }
  }

  ////////////////////////////////////////////////////////////////////////////////
  // Gen GenPhotons
  if(!fIsRealData){
    edm::Handle<GenParticleCollection> gen;
    iEvent.getByLabel(fGenPartTag, gen);
    GenParticleCollection::const_iterator g_part;
	
    // Steve Mrenna's status 2 parton jets
    edm::Handle<GenJetCollection> partonGenJets;
    bool partonGenJets_found = iEvent.getByLabel("partonGenJets", partonGenJets);
    GenJetCollection::const_iterator pGenJet;

    edm::Handle<View<Candidate> > partons;
    bool partons_found = iEvent.getByLabel("partons", partons);
	 
    std::vector<const GenParticle*> gen_photons;
    std::vector<const GenParticle*> gen_photons_mothers;

    for(g_part = gen->begin(); g_part != gen->end(); g_part++){

      if( g_part->pdgId() != 22 ) continue;
      if( g_part->status()!=1 ) continue;
      if( g_part->pt() < fMinGenPhotPt )  continue;
      if( fabs(g_part->eta()) > fMaxGenPhotEta ) continue;

      const GenParticle* gen_phot = &(*g_part);
      const GenParticle* gen_phot_mom = static_cast<const GenParticle*> (g_part->mother());

      if(gen_phot_mom==NULL){
        edm::LogWarning("NTP") << "@SUB=analyze" << " WARNING: GenPhoton does not have a mother ";
      }

      gen_photons.push_back(gen_phot);
      gen_photons_mothers.push_back(gen_phot_mom);

    }

    *fTNGenPhotons = gen_photons.size();

    for(int i=0; i<*fTNGenPhotons; ++i){
      if( i >= gMaxNGenPhot){
        edm::LogWarning("NTP") << "@SUB=analyze" << "Maximum number of gen-photons exceeded..";
        *fTMaxPhotonsExceed = 1;
        *fTGoodEvent = 1;
        break;
      }

      fTGenPhotonPt->push_back( gen_photons[i]->pt() );
      fTGenPhotonEta->push_back( gen_photons[i]->eta() );
      fTGenPhotonPhi->push_back( gen_photons[i]->phi() );
      fTGenPhotonMotherID->push_back( gen_photons_mothers[i]!=NULL ? gen_photons_mothers[i]->pdgId() : -999 );
      fTGenPhotonMotherStatus->push_back( gen_photons_mothers[i]!=NULL ? gen_photons_mothers[i]->status() : -999 );

      // use Steve Mrenna's status 2 parton jets to compute dR to closest jet of prompt photon
      if((*fTGenPhotonMotherStatus)[i]!=3) continue;
      TLorentzVector photon(0,0,0,0);
      photon.SetPtEtaPhiM((*fTGenPhotonPt)[i],(*fTGenPhotonEta)[i],(*fTGenPhotonPhi)[i],0);
      double minDR=10;
      for(pGenJet = partonGenJets->begin(); pGenJet != partonGenJets->end(); pGenJet++){
        TLorentzVector pJ(pGenJet->px(), pGenJet->py(), pGenJet->pz(), pGenJet->energy());
        float dR = photon.DeltaR(pJ);
        if(dR < minDR) minDR=dR;
      }
      fTGenPhotonPartonMindR->push_back( minDR );
    }
  }


  ////////////////////////////////////////////////////////////////////////////////
  // Get GenJets
  if(!fIsRealData){
    edm::Handle<GenJetCollection> genjets;
    iEvent.getByLabel(fGenJetTag, genjets);
    GenJetCollection::const_iterator gjet;
		
    int jqi=-1;
    for(gjet = genjets->begin(); gjet != genjets->end(); gjet++) {
      // Preselection
      if(gjet->pt() < fMinGenJetPt) continue;
      if(fabs(gjet->eta()) > fMaxGenJetEta) continue;
      jqi++;
      if( jqi >= gMaxNGenJets){
        edm::LogWarning("NTP") << "@SUB=analyze"
                               << "Maximum number of gen-jets exceeded..";
        *fTMaxGenJetExceed = 1;
        *fTGoodEvent = 1;
        break;
      }
			
      fTGenJetPt  ->push_back( gjet->pt() );
      fTGenJetEta ->push_back( gjet->eta() );
      fTGenJetPhi ->push_back( gjet->phi() );
      fTGenJetE   ->push_back( gjet->energy() );
      fTGenJetEmE ->push_back( gjet->emEnergy() );
      fTGenJetHadE->push_back( gjet->hadEnergy() );
      fTGenJetInvE->push_back( gjet->invisibleEnergy() );
    }
    *fTNGenJets = jqi+1;
  }


  ////////////////////////////////////////////////////////
  // Muon Variables:
  int mqi(0);  // Index of qualified muons
  *fTNMusTot = 0; // Total number of tracker&&global muons

  // Get muons, order them by pt and apply selection
  std::vector<OrderPair> muOrdered;
  int muIndex(0);
  for ( View<Muon>::const_iterator Mit = muons->begin(); Mit != muons->end();
        ++Mit,++muIndex ) {
    // Check if maximum number of muons is exceeded already:
    if(mqi >= gMaxNMus){
      edm::LogWarning("NTP") << "@SUB=analyze()"
                             << "Maximum number of muons exceeded";
      *fTMaxMuExceed = 1;
      *fTGoodEvent = 1;
      break;
    }
    // Muon preselection:
    // Only consider global and trackermuons:
    if(!(Mit->isGlobalMuon()) && !(Mit->isTrackerMuon())) continue;
    (*fTNMusTot)++;     // Count all
    if(Mit->pt() < fMinMuPt) continue;
    if(fabs(Mit->eta()) > fMaxMuEta) continue;
    ++mqi;          // Count how many we'll eventually store
    muOrdered.push_back(make_pair(muIndex,Mit->pt()));
  }
  std::sort(muOrdered.begin(),muOrdered.end(),indexComparator);
  *fTNMus = muOrdered.size();
  mqi = 0;

  // Dump muon properties in tree variables
  for (std::vector<OrderPair>::const_iterator it = muOrdered.begin();
       it != muOrdered.end(); ++it, ++mqi ) {
    int index = it->first;
    const Muon& muon = (*muons)[index];

    fTMuIsGlobalMuon->push_back( muon.isGlobalMuon() ? 1:0 );
    fTMuIsTrackerMuon->push_back( muon.isTrackerMuon() ? 1:0 );

    // Combined methods for Global and Tracker muons:
    fTMuPx->push_back( muon.px() );
    fTMuPy->push_back( muon.py() );
    fTMuPz->push_back( muon.pz() );
    fTMuPt->push_back( muon.pt() );
    fTMuInnerTkPt->push_back( muon.innerTrack()->pt() );
    fTMuEta->push_back( muon.eta() );
    fTMuPhi->push_back( muon.phi() );
    fTMuE->push_back( muon.energy() );
    fTMuEt->push_back( muon.et() );
    fTMuCharge->push_back( muon.charge() );

    fTMuRelIso03->push_back( (muon.isolationR03().sumPt + muon.isolationR03().emEt + muon.isolationR03().hadEt) / muon.pt() );
    fTMuIso03SumPt->push_back( muon.isolationR03().sumPt );
    fTMuIso03EmEt->push_back( muon.isolationR03().emEt );
    fTMuIso03HadEt->push_back( muon.isolationR03().hadEt );
    fTMuIso03EMVetoEt->push_back( muon.isolationR03().emVetoEt );
    fTMuIso03HadVetoEt->push_back( muon.isolationR03().hadVetoEt );
    fTMuIso05SumPt->push_back( muon.isolationR05().sumPt );
    fTMuIso05EmEt->push_back( muon.isolationR05().emEt );
    fTMuIso05HadEt->push_back( muon.isolationR05().hadEt );

    fTMuCaloComp->push_back( muon.caloCompatibility() );
    fTMuSegmComp->push_back( muon::segmentCompatibility(muon) );

    // MuID Flags:
    fTMuIsGMPT->push_back( muon::isGoodMuon(muon, muon::GlobalMuonPromptTight) ? 1:0 );
    fTMuIsGMTkChiComp->push_back( muon::isGoodMuon(muon, muon::GMTkChiCompatibility) ? 1:0 );
    fTMuIsGMStaChiComp->push_back( muon::isGoodMuon(muon, muon::GMStaChiCompatibility) ? 1:0 );
    fTMuIsGMTkKinkTight->push_back( muon::isGoodMuon(muon, muon::GMTkKinkTight) ? 1:0 );
    fTMuIsAllStaMuons->push_back( muon::isGoodMuon(muon, muon::AllStandAloneMuons) ? 1:0 );
    fTMuIsAllTrkMuons->push_back( muon::isGoodMuon(muon, muon::AllTrackerMuons) ? 1:0 );
    fTMuIsTrkMuonArbitrated->push_back( muon::isGoodMuon(muon, muon::TrackerMuonArbitrated) ? 1:0 );
    fTMuIsAllArbitrated->push_back( muon::isGoodMuon(muon, muon::AllArbitrated) ? 1:0 );
    fTMuIsTMLSLoose->push_back( muon::isGoodMuon(muon, muon::TMLastStationLoose) ? 1:0 );
    fTMuIsTMLSTight->push_back( muon::isGoodMuon(muon, muon::TMLastStationTight) ? 1:0 );
    fTMuIsTM2DCompLoose->push_back( muon::isGoodMuon(muon, muon::TM2DCompatibilityLoose) ? 1:0 );
    fTMuIsTM2DCompTight->push_back( muon::isGoodMuon(muon, muon::TM2DCompatibilityTight) ? 1:0 );
    fTMuIsTMOneStationLoose->push_back( muon::isGoodMuon(muon, muon::TMOneStationLoose) ? 1:0 );
    fTMuIsTMOneStationTight->push_back( muon::isGoodMuon(muon, muon::TMOneStationTight) ? 1:0 );
    fTMuIsTMLSOptLowPtLoose->push_back( muon::isGoodMuon(muon, muon::TMLastStationOptimizedLowPtLoose) ? 1:0 );
    fTMuIsTMLSAngLoose->push_back( muon::isGoodMuon(muon, muon::TMLastStationAngLoose) ? 1:0 );
    fTMuIsTMLSAngTight->push_back( muon::isGoodMuon(muon, muon::TMLastStationAngTight) ? 1:0 );
    fTMuIsTMOneStationAngLoose->push_back( muon::isGoodMuon(muon, muon::TMOneStationAngLoose) ? 1:0 );
    fTMuIsTMOneStationAngTight->push_back( muon::isGoodMuon(muon, muon::TMOneStationAngTight) ? 1:0 );

    Ref<View<Muon> > muonRef(muons,index);
    const reco::IsoDeposit ECDep = ECDepMap[muonRef];
    const reco::IsoDeposit HCDep = HCDepMap[muonRef];
    fTMuEem->push_back( ECDep.candEnergy() );
    fTMuEhad->push_back( HCDep.candEnergy() );

    fTMuD0BS->push_back( -1.0*muon.innerTrack()->dxy(beamSpot.position()) );
    fTMuD0PV->push_back( -1.0*muon.innerTrack()->dxy(primVtx->position()) );
    fTMuDzBS->push_back( muon.innerTrack()->dz(beamSpot.position()) );
    fTMuDzPV->push_back( muon.innerTrack()->dz(primVtx->position()) );
    fTMuInnerTkNChi2->push_back( muon.innerTrack()->normalizedChi2() );

    // Separate methods:
    if((*fTMuIsTrackerMuon)[mqi]){ // Tracker Muons
      (*fTNTMus)++;
      fTMuPtE->push_back( muon.innerTrack()->ptError() );
      fTMuD0E->push_back( muon.innerTrack()->dxyError() );
      fTMuDzE->push_back( muon.innerTrack()->dzError() );

      fTMuNChi2->push_back( (*fTMuInnerTkNChi2)[mqi] ); // No difference for TM
      fTMuNGlHits->push_back( 0 );
      fTMuNTkHits->push_back( muon.innerTrack()->hitPattern().numberOfValidHits() );
      fTMuNPxHits->push_back( muon.innerTrack()->hitPattern().numberOfValidPixelHits() );
      fTMuNMuHits->push_back( 0 );
      fTMuNMatches->push_back( 0 );
      fTMuNChambers->push_back( 0 );
    }
    if((*fTMuIsGlobalMuon)[mqi]){ // Global Muons
      (*fTNGMus)++;
      fTMuPtE->push_back( muon.globalTrack()->ptError() );
      fTMuD0E->push_back( muon.globalTrack()->dxyError() );
      fTMuDzE->push_back( muon.globalTrack()->dzError() );

      fTMuNChi2->push_back( muon.globalTrack()->normalizedChi2() );
      fTMuNGlHits->push_back( muon.globalTrack()->hitPattern().numberOfValidHits() );
      fTMuNTkHits->push_back( muon.innerTrack()->hitPattern().numberOfValidHits() );
      fTMuNPxHits->push_back( muon.innerTrack()->hitPattern().numberOfValidPixelHits() );
      fTMuNMuHits->push_back( muon.outerTrack()->hitPattern().numberOfValidHits() );
      fTMuNMatches->push_back( muon.numberOfMatches() );
      fTMuNChambers->push_back( muon.numberOfChambers() );
    }

    // MC Matching
    if(!fIsRealData){
      std::vector<const GenParticle*> MuMatch = matchRecoCand(&muon, iEvent);
      if(MuMatch[0] != NULL){
        fTMuGenID->push_back( MuMatch[0]->pdgId() );
        fTMuGenStatus->push_back( MuMatch[0]->status() );
        fTMuGenPt->push_back( MuMatch[0]->pt() );
        fTMuGenEta->push_back( MuMatch[0]->eta() );
        fTMuGenPhi->push_back( MuMatch[0]->phi() );
        fTMuGenE->push_back( MuMatch[0]->energy() );

        fTMuGenMID->push_back( MuMatch[1]->pdgId() );
        fTMuGenMStatus->push_back( MuMatch[1]->status() );
        fTMuGenMPt->push_back( MuMatch[1]->pt() );
        fTMuGenMEta->push_back( MuMatch[1]->eta() );
        fTMuGenMPhi->push_back( MuMatch[1]->phi() );
        fTMuGenME->push_back( MuMatch[1]->energy() );

        fTMuGenGMID->push_back( MuMatch[2]->pdgId() );
        fTMuGenGMStatus->push_back( MuMatch[2]->status() );
        fTMuGenGMPt->push_back( MuMatch[2]->pt() );
        fTMuGenGMEta->push_back( MuMatch[2]->eta() );
        fTMuGenGMPhi->push_back( MuMatch[2]->phi() );
        fTMuGenGME->push_back( MuMatch[2]->energy() );
      }
      MuMatch.clear();
    }
    fTMuGood->push_back( 0 );
    fTMuIsIso->push_back( 1 );
  }


  // SC variables
  (*fTNSuperClusters)=0;
  for (SuperClusterCollection::const_iterator sc = BarrelSuperClusters->begin(); sc!=BarrelSuperClusters->end(); ++sc){

    if (sc->rawEnergy()<fMinSCraw) continue;

    if (*fTNSuperClusters>gMaxNSC) {
      edm::LogWarning("NTP") << "@SUB=analyze" << "Maximum number of Super Clusters exceeded"; 
      *fTGoodEvent = 1; 
      break;
    }

    fTSCRaw->push_back(sc->rawEnergy());
    fTSCPre->push_back(sc->preshowerEnergy());
    fTSCEnergy->push_back(sc->energy());
    fTSCEta->push_back(sc->eta());
    fTSCPhi->push_back(sc->phi());
    fTSCPhiWidth->push_back(sc->phiWidth());
    fTSCEtaWidth->push_back(sc->etaWidth());
    fTSCBrem->push_back((sc->etaWidth()!=0) ? sc->phiWidth()/sc->etaWidth() : -1);
    fTSCR9->push_back(sc->rawEnergy()!=0 ? EcalClusterTools::e3x3(  *(sc->seed()), ebRecHits.product(), &(*topology)) / sc->rawEnergy() : -1);
    {
      float crackcorrseedenergy = sc->rawEnergy();
      float localcorrseedenergy = sc->rawEnergy();
      float crackcorrenergy = sc->rawEnergy();
      float localcorrenergy = sc->rawEnergy();
      int index=0;
      for(reco::CaloCluster_iterator itClus = sc->clustersBegin(); itClus != sc->clustersEnd(); ++itClus) {
        const reco::CaloClusterPtr cc = *itClus;
        if (&(**itClus)==&(*sc->seed())){
          crackcorrseedenergy += (*itClus)->energy()*(CrackCorrFunc->getValue(*cc)-1);
          localcorrseedenergy += (*itClus)->energy()*(LocalCorrFunc->getValue(*cc)-1);
          fTSCcrackcorrseedfactor->push_back(CrackCorrFunc->getValue(*cc));
          fTSClocalcorrseedfactor->push_back(LocalCorrFunc->getValue(*cc));		
        }
        crackcorrenergy += (*itClus)->energy()*(CrackCorrFunc->getValue(*cc)-1);
        localcorrenergy += (*itClus)->energy()*(LocalCorrFunc->getValue(*cc)-1);
        index++;
      }
      fTSCcrackcorrseed->push_back(crackcorrseedenergy/sc->rawEnergy());
      fTSCcrackcorr->push_back(crackcorrenergy/sc->rawEnergy());
      fTSClocalcorrseed->push_back(localcorrseedenergy/sc->rawEnergy());
      fTSClocalcorr->push_back(localcorrenergy/sc->rawEnergy());
    }
    (*fTNSuperClusters)++;
  }

  for (SuperClusterCollection::const_iterator sc = EndcapSuperClusters->begin(); sc!=EndcapSuperClusters->end(); ++sc){

    if (sc->rawEnergy()<fMinSCraw) continue;

    if (*fTNSuperClusters>gMaxNSC) {
      edm::LogWarning("NTP") << "@SUB=analyze" << "Maximum number of Super Clusters exceeded"; 
      *fTGoodEvent = 1; 
      break;
    }

    fTSCRaw->push_back(sc->rawEnergy());
    fTSCPre->push_back(sc->preshowerEnergy());
    fTSCEnergy->push_back(sc->energy());
    fTSCEta->push_back(sc->eta());
    fTSCPhi->push_back(sc->phi());
    fTSCPhiWidth->push_back(sc->phiWidth());
    fTSCEtaWidth->push_back(sc->etaWidth());
    fTSCBrem->push_back((sc->etaWidth()!=0) ? sc->phiWidth()/sc->etaWidth() : -1);
    fTSCR9->push_back(sc->rawEnergy()!=0 ? EcalClusterTools::e3x3(  *(sc->seed()), eeRecHits.product(), &(*topology)) / sc->rawEnergy() : -1);
    {
      float crackcorrseedenergy = sc->rawEnergy();
      float localcorrseedenergy = sc->rawEnergy();
      float crackcorrenergy = sc->rawEnergy();
      float localcorrenergy = sc->rawEnergy();
      int index=0;
      for(reco::CaloCluster_iterator itClus = sc->clustersBegin(); itClus != sc->clustersEnd(); ++itClus) {
        const reco::CaloClusterPtr cc = *itClus;
        if (&(**itClus)==&(*sc->seed())){
          crackcorrseedenergy += (*itClus)->energy()*(CrackCorrFunc->getValue(*cc)-1);
          localcorrseedenergy += (*itClus)->energy()*(LocalCorrFunc->getValue(*cc)-1);
          fTSCcrackcorrseedfactor->push_back(CrackCorrFunc->getValue(*cc));
          fTSClocalcorrseedfactor->push_back(LocalCorrFunc->getValue(*cc));		
        }
        crackcorrenergy += (*itClus)->energy()*(CrackCorrFunc->getValue(*cc)-1);
        localcorrenergy += (*itClus)->energy()*(LocalCorrFunc->getValue(*cc)-1);
        index++;
      }
      fTSCcrackcorrseed->push_back(crackcorrseedenergy/sc->rawEnergy());
      fTSCcrackcorr->push_back(crackcorrenergy/sc->rawEnergy());
      fTSClocalcorrseed->push_back(localcorrseedenergy/sc->rawEnergy());
      fTSClocalcorr->push_back(localcorrenergy/sc->rawEnergy());
    }
    (*fTNSuperClusters)++;
  }




  ////////////////////////////////////////////////////////
  // Electron variables:
  // Keep pointers to electron superCluster in original collections
  std::vector<const SuperCluster*> elecPtr;
  std::vector<const GsfTrack*> trckPtr;
  int eqi(0);                    // Index of qualified electrons
  (*fTNElesTot) = electrons->size(); // Total number of electrons

  if (electrons->size() > 0) {
    // Get electrons, order them by pt and apply selection
    std::vector<OrderPair> elOrdered;
    int elIndex(0);
    for( View<GsfElectron>::const_iterator El = electrons->begin();
         El != electrons->end(); ++El, ++elIndex ) {
      // Check if maximum number of electrons is exceeded already:
      if(eqi >= gMaxNEles) {
        edm::LogWarning("NTP") << "@SUB=analyze"
                               << "Maximum number of electrons exceeded..";
        *fTMaxElExceed = 1;
        *fTGoodEvent = 1;
        break;
      }
      // Electron preselection:
      if(El->pt() < fMinElPt) continue;
      if(fabs(El->eta()) > fMaxElEta) continue;

      eqi++; // Count how many we'll eventually store
      elOrdered.push_back(make_pair(elIndex,El->pt()));
    }
    std::sort(elOrdered.begin(),elOrdered.end(),indexComparator);
    (*fTNEles) = elOrdered.size();
    eqi = 0;

    // Read eID results
    std::vector<Handle<ValueMap<float> > > eIDValueMap(9);
    iEvent.getByLabel( "eidRobustLoose",      eIDValueMap[0]);
    iEvent.getByLabel( "eidRobustTight",      eIDValueMap[1]);
    iEvent.getByLabel( "eidLoose",            eIDValueMap[2]);
    iEvent.getByLabel( "eidTight",            eIDValueMap[3]);
    iEvent.getByLabel( fEleIdWP,              eIDValueMap[4]);
    iEvent.getByLabel( "simpleEleId80relIso", eIDValueMap[5]);
    iEvent.getByLabel( "simpleEleId85relIso", eIDValueMap[6]);
    iEvent.getByLabel( "simpleEleId90relIso", eIDValueMap[7]);
    iEvent.getByLabel( "simpleEleId95relIso", eIDValueMap[8]);
    const ValueMap<float> &eIDmapRL         = *eIDValueMap[0];  // Robust-Loose
    const ValueMap<float> &eIDmapRT         = *eIDValueMap[1];  // Robust-Tight
    const ValueMap<float> &eIDmapL          = *eIDValueMap[2];  // Loose
    const ValueMap<float> &eIDmapT          = *eIDValueMap[3];  // Tight
    const ValueMap<float> &eIDmapsimpleWP   = *eIDValueMap[4];  // WP from config
    const ValueMap<float> &eIDmapsimpleWP80 = *eIDValueMap[5];  // WP80
    const ValueMap<float> &eIDmapsimpleWP85 = *eIDValueMap[6];  // WP85
    const ValueMap<float> &eIDmapsimpleWP90 = *eIDValueMap[7];  // WP90
    const ValueMap<float> &eIDmapsimpleWP95 = *eIDValueMap[8];  // WP95
    eIDValueMap.clear();


    // Dump electron properties in tree variables
    for( std::vector<OrderPair>::const_iterator it = elOrdered.begin();
         it != elOrdered.end(); ++it, ++eqi ) {

      int index = it->first;
      const GsfElectron& electron = (*electrons)[index];

      // Save the electron SuperCluster pointer
      elecPtr.push_back(&(*electron.superCluster()));
      trckPtr.push_back(&(*electron.gsfTrack()));

      fTElPx                     ->push_back(electron.px());
      fTElPy                     ->push_back(electron.py());
      fTElPz                     ->push_back(electron.pz());
      fTElPt                     ->push_back(electron.pt());
      fTElPtE                    ->push_back(electron.gsfTrack()->ptError());
      fTElEta                    ->push_back(electron.eta());
      fTElPhi                    ->push_back(electron.phi());
      fTElGsfTkPt                ->push_back(electron.gsfTrack()->pt());
      fTElGsfTkEta               ->push_back(electron.gsfTrack()->eta());
      fTElGsfTkPhi               ->push_back(electron.gsfTrack()->phi());
      fTElTrkMomentumError       ->push_back(electron.trackMomentumError());
      fTElEcalEnergyError        ->push_back(electron.ecalEnergyError());
      // 4_2: take the error on the default momentum
      fTElEleMomentumError       ->push_back(electron.p4Error( electron.candidateP4Kind() ));
      fTElNBrems                 ->push_back(electron.numberOfBrems());
      fTElE                      ->push_back(electron.energy());
      fTElEt                     ->push_back(electron.et());
      fTElD0BS                   ->push_back(-1.0*electron.gsfTrack()->dxy(beamSpot.position()));
      fTElD0PV                   ->push_back(-1.0*electron.gsfTrack()->dxy(primVtx->position()));
      fTElD0E                    ->push_back(electron.gsfTrack()->dxyError());
      fTElDzBS                   ->push_back(electron.gsfTrack()->dz(beamSpot.position()));
      fTElDzPV                   ->push_back(electron.gsfTrack()->dz(primVtx->position()));
      fTElDzE                    ->push_back(electron.gsfTrack()->dzError());
      fTElNChi2                  ->push_back(electron.gsfTrack()->normalizedChi2());
      fTElDR03TkSumPt            ->push_back(electron.dr03TkSumPt());
      fTElDR03EcalRecHitSumEt    ->push_back(electron.dr03EcalRecHitSumEt());
      fTElDR03HcalTowerSumEt     ->push_back(electron.dr03HcalTowerSumEt());
      fTElDR04TkSumPt            ->push_back(electron.dr04TkSumPt());
      fTElDR04EcalRecHitSumEt    ->push_back(electron.dr04EcalRecHitSumEt());
      fTElDR04HcalTowerSumEt     ->push_back(electron.dr04HcalTowerSumEt());
      fTElRelIso03               ->push_back(((*fTElDR03TkSumPt)[eqi] + (*fTElDR03EcalRecHitSumEt)[eqi] + (*fTElDR03HcalTowerSumEt)[eqi]) / (*fTElPt)[eqi]);
      fTElRelIso04               ->push_back(((*fTElDR04TkSumPt)[eqi] + (*fTElDR04EcalRecHitSumEt)[eqi] + (*fTElDR04HcalTowerSumEt)[eqi]) / (*fTElPt)[eqi]);
      fTElCharge                 ->push_back(electron.charge());
      fTElInGap                  ->push_back(electron.isGap() ? 1:0);
      fTElEcalDriven             ->push_back(electron.ecalDrivenSeed() ? 1:0);
      fTElTrackerDriven          ->push_back(electron.trackerDrivenSeed() ? 1:0);
      fTElCInfoIsGsfCtfCons      ->push_back(electron.chargeInfo().isGsfCtfConsistent ? 1:0);
      fTElCInfoIsGsfCtfScPixCons ->push_back(electron.chargeInfo().isGsfCtfScPixConsistent ? 1:0);
      fTElCInfoIsGsfScPixCons    ->push_back(electron.chargeInfo().isGsfScPixConsistent ? 1:0);
      fTElScPixCharge            ->push_back(electron.chargeInfo().scPixCharge);
      if( electron.closestCtfTrackRef().isNonnull() ) {
        fTElClosestCtfTrackPt      ->push_back((electron.closestCtfTrack().ctfTrack)->pt());
        fTElClosestCtfTrackEta     ->push_back((electron.closestCtfTrack().ctfTrack)->eta());
        fTElClosestCtfTrackPhi     ->push_back((electron.closestCtfTrack().ctfTrack)->phi());
        fTElClosestCtfTrackCharge  ->push_back((electron.closestCtfTrack().ctfTrack)->charge());
      }
      fTElBasicClustersSize         ->push_back(electron.basicClustersSize());
      fTElfbrem                     ->push_back(electron.fbrem());
      fTElHcalOverEcal              ->push_back(electron.hcalOverEcal());
      fTElE1x5                      ->push_back(electron.e1x5());
      fTElE5x5                      ->push_back(electron.e5x5());
      fTElE2x5Max                   ->push_back(electron.e2x5Max());
      fTElSigmaIetaIeta             ->push_back(electron.sigmaIetaIeta());
      fTElDeltaEtaSeedClusterAtCalo ->push_back(electron.deltaEtaSeedClusterTrackAtCalo());
      fTElDeltaPhiSeedClusterAtCalo ->push_back(electron.deltaPhiSeedClusterTrackAtCalo());
      fTElDeltaPhiSuperClusterAtVtx ->push_back(electron.deltaPhiSuperClusterTrackAtVtx());
      fTElDeltaEtaSuperClusterAtVtx ->push_back(electron.deltaEtaSuperClusterTrackAtVtx());
      fTElCaloEnergy                ->push_back(electron.caloEnergy());
      fTElTrkMomAtVtx               ->push_back(electron.trackMomentumAtVtx().R());
      fTElESuperClusterOverP        ->push_back(electron.eSuperClusterOverP());
      fTElNumberOfMissingInnerHits  ->push_back(electron.gsfTrack()->trackerExpectedHitsInner().numberOfHits());
      fTElTheta                     ->push_back(electron.superCluster()->position().theta());
      fTElSCEta                     ->push_back(electron.superCluster()->eta());



      // DISABLED: NO SEED IN AOD (UPDATE IT IN 4_2)
      // 			if ( electron.superCluster()->seed()->caloID().detector( reco::CaloID::DET_ECAL_BARREL ) ) {
      //                           fTeScSeedSeverity ->push_back(EcalSeverityLevelAlgo::severityLevel( electron.superCluster()->seed()->seed(), *ebRecHits, *channelStatus ));
      //                           fTeE1OverE9       ->push_back(EcalSeverityLevelAlgo::E1OverE9(   electron.superCluster()->seed()->seed(), *ebRecHits ));
      //                           fTeS4OverS1       ->push_back(EcalSeverityLevelAlgo::swissCross( electron.superCluster()->seed()->seed(), *ebRecHits ));
      // 			} else if ( electron.superCluster()->seed()->caloID().detector( reco::CaloID::DET_ECAL_ENDCAP ) ) {
      //                           fTeScSeedSeverity ->push_back(EcalSeverityLevelAlgo::severityLevel( electron.superCluster()->seed()->seed(), *eeRecHits, *channelStatus ));
      //                           fTeE1OverE9       ->push_back(EcalSeverityLevelAlgo::E1OverE9(   electron.superCluster()->seed()->seed(), *eeRecHits ));
      //                           fTeS4OverS1       ->push_back(EcalSeverityLevelAlgo::swissCross( electron.superCluster()->seed()->seed(), *eeRecHits ));
      // 			} else {
      // 				edm::LogWarning("NTP") << "Electron supercluster seed crystal neither in EB nor in EE!";
      // 			}

      // Read in Electron ID
      fTElIDMva ->push_back(electron.mva());
      Ref<View<GsfElectron> > electronRef(electrons,index);
      fTElIDTight            ->push_back(eIDmapT[electronRef]  ? 1:0);
      fTElIDLoose            ->push_back(eIDmapL[electronRef]  ? 1:0);
      fTElIDRobustTight      ->push_back(eIDmapRT[electronRef] ? 1:0);
      fTElIDRobustLoose      ->push_back(eIDmapRL[electronRef] ? 1:0);
      fTElIDsimpleWPrelIso   ->push_back(eIDmapsimpleWP[electronRef]);
      fTElIDsimpleWP95relIso ->push_back(eIDmapsimpleWP95[electronRef]);
      fTElIDsimpleWP90relIso ->push_back(eIDmapsimpleWP90[electronRef]);
      fTElIDsimpleWP85relIso ->push_back(eIDmapsimpleWP85[electronRef]);
      fTElIDsimpleWP80relIso ->push_back(eIDmapsimpleWP80[electronRef]);

      {
        (*fTElSCindex)[eqi] = -1;
        float diff=1e+4;
        for (int scind=0; scind<*fTNSuperClusters; scind++){
          if (fabs((*fTSCRaw)[scind]-electron.superCluster()->rawEnergy())<diff) {
            (*fTElSCindex)[eqi]=scind;
            diff=fabs((*fTSCRaw)[scind]-electron.superCluster()->rawEnergy());
          }
        }

        if ((*fTElSCindex)[eqi]!=-1)
          if (fabs((*fTSCEta)[(*fTElSCindex)[eqi]]-electron.superCluster()->eta())>0.1 || 
              DeltaPhi((*fTSCPhi)[(*fTElSCindex)[eqi]],electron.superCluster()->phi())>0.1){
            // std::cout << fTSCraw[fTElSCindex[eqi]] << " " << electron.superCluster()->rawEnergy() << std::endl;
            // std::cout << fTSCeta[fTElSCindex[eqi]] << " " << electron.superCluster()->eta() << std::endl;
            // std::cout << fTSCphi[fTElSCindex[eqi]] << " " << electron.superCluster()->phi() << std::endl;
            (*fTElSCindex)[eqi] = -1;			    
          }

        if ((*fTElSCindex)[eqi]==-1) {
          //edm::LogWarning("NTP") << "@SUB=analyze" << "No matching SC found for electron"; 
          //	*fTGoodEvent = 1; 
          //	break;
        }
      }


      // FIXME: TO BE REMOVED ONCE WE ARE HAPPY WITH THE FULL GEN. INFO
      // MC Matching
      if(!fIsRealData){
        std::vector<const GenParticle*> ElMatch = matchRecoCand(&electron, iEvent);
        if(ElMatch[0] != NULL){
          fTElGenID       ->push_back(ElMatch[0]->pdgId());
          fTElGenStatus   ->push_back(ElMatch[0]->status());
          fTElGenPt       ->push_back(ElMatch[0]->pt());
          fTElGenEta      ->push_back(ElMatch[0]->eta());
          fTElGenPhi      ->push_back(ElMatch[0]->phi());
          fTElGenE        ->push_back(ElMatch[0]->energy());

          fTElGenMID      ->push_back(ElMatch[1]->pdgId());
          fTElGenMStatus  ->push_back(ElMatch[1]->status());
          fTElGenMPt      ->push_back(ElMatch[1]->pt());
          fTElGenMEta     ->push_back(ElMatch[1]->eta());
          fTElGenMPhi     ->push_back(ElMatch[1]->phi());
          fTElGenME       ->push_back(ElMatch[1]->energy());

          fTElGenGMID     ->push_back(ElMatch[2]->pdgId());
          fTElGenGMStatus ->push_back(ElMatch[2]->status());
          fTElGenGMPt     ->push_back(ElMatch[2]->pt());
          fTElGenGMEta    ->push_back(ElMatch[2]->eta());
          fTElGenGMPhi    ->push_back(ElMatch[2]->phi());
          fTElGenGME      ->push_back(ElMatch[2]->energy());
        }
        ElMatch.clear();
      }


      // Conversion Information
      reco::GsfElectron::ConversionRejection ConvRejVars = electron.conversionRejectionVariables();
      reco::TrackBaseRef ConvPartnerTrack = ConvRejVars.partner;
      if( ConvPartnerTrack.isNonnull() ){
        fTElConvPartnerTrkDist   ->push_back(ConvRejVars.dist);
        fTElConvPartnerTrkDCot   ->push_back(ConvRejVars.dcot);
        fTElConvPartnerTrkPt     ->push_back(ConvPartnerTrack->pt());
        fTElConvPartnerTrkEta    ->push_back(ConvPartnerTrack->eta());
        fTElConvPartnerTrkPhi    ->push_back(ConvPartnerTrack->phi());
        fTElConvPartnerTrkCharge ->push_back(ConvPartnerTrack->charge());
      }


      fTElGood            ->push_back(0);
      fTElIsIso           ->push_back(1);
      fTElChargeMisIDProb ->push_back(0);
    }
  }

  (*fTNEBhits) = 0;
  for(EcalRecHitCollection::const_iterator ecalrechit = ebRecHits->begin(); ecalrechit!=ebRecHits->end() ; ++ecalrechit)
    {
      double energy = ecalrechit->energy();
      if(energy<fMinEBRechitE)continue; 
      if((*fTNEBhits)>gMaxNEBhits)
        {
          edm::LogWarning("NTP") << "@SUB=analyze" << "Maximum number of EB rechits exceeded"; 
          *fTGoodEvent = 1; 
          break;
        }

      double time = ecalrechit->time();
      double chi2 = ecalrechit->chi2();

      EBDetId ebDetId = (*ecalrechit).id();	
      const GlobalPoint p ( geometry->getPosition( ebDetId ) ) ;
      TVector3 hitPos(p.x(),p.y(),p.z());
      hitPos *= 1.0/hitPos.Mag();
      hitPos *= energy;	
      //4_2 		float e4oe1 = EcalSeverityLevelAlgo::swissCross( ebDetId , *ebRecHits );
      //4_2 		float e2oe9 =EcalSeverityLevelAlgo::E2overE9(ebDetId , *ebRecHits );

      fTEBrechitE    ->push_back( hitPos.Mag());
      fTEBrechitPt   ->push_back( hitPos.Pt());
      fTEBrechitEta  ->push_back( hitPos.Eta());
      fTEBrechitPhi  ->push_back(hitPos.Phi() );
      fTEBrechitTime ->push_back(time);
      fTEBrechitChi2 ->push_back(chi2);
      //4_2 		fTEBrechitE4oE1[fTnEBhits] =e4oe1 ;
      //4_2 		fTEBrechitE2oE9[fTnEBhits] = e2oe9;

      //	cout << "ebrechit P =" << fTEBrechitE[fTnEBhits] << " Pt = " << fTEBrechitPt[fTnEBhits] ;
      //	cout << " eta = " << fTEBrechitEta[fTnEBhits] << "phi = " << fTEBrechitPhi[fTnEBhits] << " time = "<< fTEBrechitTime[fTnEBhits] ;
      //	cout << " chi2 = "<< fTEBrechitChi2[fTnEBhits] << " e4oe1 = "<< fTEBrechitE4oE1[fTnEBhits] << " e2oe9 = " << e2oe9 << endl;

      (*fTNEBhits)++;
    }

	


  ////////////////////////////////////////////////////////
  // Photon Variables:
  // Keep pointers to superclusters for cross cleaning
  std::vector<const SuperCluster*> photSCs;
  int phoqi(0); // Index of qualified photons
  (*fTNPhotonsTot) = photons->size();



  // Get photonss, order them by pt and apply selection
  std::vector<OrderPair> phoOrdered;
  int phoIndex(0);
  for( View<Photon>::const_iterator ip = photons->begin();
       ip != photons->end(); ++ip, ++phoIndex ){
    // Check if maximum number of photons exceeded
    if(phoqi >= gMaxNPhotons){
      edm::LogWarning("NTP") << "@SUB=analyze"
                             << "Maximum number of photons exceeded";
      *fTMaxPhotonsExceed = 1;
      *fTGoodEvent = 1;
      break;
    }
    // Preselection
    if(ip->pt() < fMinPhotonPt) continue;
    if(fabs(ip->eta()) > fMaxPhotonEta) continue;

    phoqi++; // Count how many we'll eventually store
    phoOrdered.push_back(make_pair(phoIndex,ip->pt()));
  }
  std::sort(phoOrdered.begin(),phoOrdered.end(),indexComparator);
  (*fTNPhotons) = phoOrdered.size();
  phoqi = 0;

  for (std::vector<OrderPair>::const_iterator it = phoOrdered.begin();
       it != phoOrdered.end(); ++it, ++phoqi ) {

    int index = it->first;
    const Photon& photon = (*photons)[index];

    // Save photon supercluster position
    photSCs.push_back(&(*photon.superCluster()));

    fTPhoPt              ->push_back(photon.pt());
    fTPhoPx              ->push_back(photon.px());
    fTPhoPy              ->push_back(photon.py());
    fTPhoPz              ->push_back(photon.pz());
    fTPhoEta             ->push_back(photon.eta());
    fTPhoPhi             ->push_back(photon.phi());
    fTPhoEnergy          ->push_back(photon.energy());
    fTPhoIso03Ecal       ->push_back(photon.ecalRecHitSumEtConeDR03());
    fTPhoIso03Hcal       ->push_back(photon.hcalTowerSumEtConeDR03());
    fTPhoIso03TrkSolid   ->push_back(photon.trkSumPtSolidConeDR03());
    fTPhoIso03TrkHollow  ->push_back(photon.trkSumPtHollowConeDR03());
    fTPhoIso03           ->push_back(((*fTPhoIso03TrkHollow)[phoqi] + (*fTPhoIso03Ecal)[phoqi] + (*fTPhoIso03Hcal)[phoqi]) / (*fTPhoPt)[phoqi]);
    fTPhoIso04Ecal       ->push_back(photon.ecalRecHitSumEtConeDR04());
    fTPhoIso04Hcal       ->push_back(photon.hcalTowerSumEtConeDR04());
    fTPhoIso04TrkSolid   ->push_back(photon.trkSumPtSolidConeDR04());
    fTPhoIso04TrkHollow  ->push_back(photon.trkSumPtHollowConeDR04());
    fTPhoIso04           ->push_back(((*fTPhoIso04TrkHollow)[phoqi] + (*fTPhoIso04Ecal)[phoqi] + (*fTPhoIso04Hcal)[phoqi]) / (*fTPhoPt)[phoqi]);
    fTPhoR9              ->push_back(photon.r9());
    fTPhoCaloPositionX   ->push_back(photon.caloPosition().X());
    fTPhoCaloPositionY   ->push_back(photon.caloPosition().Y());
    fTPhoCaloPositionZ   ->push_back(photon.caloPosition().Z());
    fTPhoHoverE          ->push_back(photon.hadronicOverEm());
    fTPhoH1overE         ->push_back(photon.hadronicDepth1OverEm());
    fTPhoH2overE         ->push_back(photon.hadronicDepth2OverEm());
    fTPhoSigmaIetaIeta   ->push_back(photon.sigmaIetaIeta());
    fTPhoSigmaEtaEta     ->push_back(photon.sigmaEtaEta());

    fTPhoE1x5 ->push_back(photon.e1x5());
    fTPhoE2x5 ->push_back(photon.e2x5());
    fTPhoE3x3 ->push_back(photon.e3x3());
    fTPhoE5x5 ->push_back(photon.e5x5());
    fTPhomaxEnergyXtal   ->push_back(photon.maxEnergyXtal());
    fTPhoIso03HcalDepth1 ->push_back(photon.hcalDepth1TowerSumEtConeDR03());
    fTPhoIso03HcalDepth2 ->push_back(photon.hcalDepth2TowerSumEtConeDR03());
    fTPhoIso04HcalDepth1 ->push_back(photon.hcalDepth1TowerSumEtConeDR04());
    fTPhoIso04HcalDepth2 ->push_back(photon.hcalDepth2TowerSumEtConeDR04());
    fTPhoIso03nTrksSolid ->push_back(photon.nTrkSolidConeDR03());
    fTPhoIso03nTrksHollow->push_back(photon.nTrkHollowConeDR03());
    fTPhoIso04nTrksSolid ->push_back(photon.nTrkSolidConeDR04());
    fTPhoIso04nTrksHollow->push_back(photon.nTrkHollowConeDR04());

    fTPhoSCRawEnergy       ->push_back(photon.superCluster()->rawEnergy());
    fTPhoSCEtaWidth     ->push_back(photon.superCluster()->etaWidth());
    fTPhoHasPixSeed     ->push_back(photon.hasPixelSeed() ? 1:0);
    fTPhoHasConvTrks    ->push_back(photon.hasConversionTracks() ? 1:0);
    // DISABLED: NO SEED IN AOD
    // EcalClusterLazyTools *lazyTools = new EcalClusterLazyTools(iEvent, iSetup, edm::InputTag("reducedEcalRecHitsEB"), edm::InputTag("reducedEcalRecHitsEE"));
    // 	fTPhoSCSigmaPhiPhi  ->push_back(lazyTools->covariances(*(photon.superCluster()->seed())).at(2));

    fTPhoGood     ->push_back(0);
    fTPhoIsIso    ->push_back(1);

    fTPhoisEB->push_back(photon.isEB());
    fTPhoisEE->push_back(photon.isEE());
    fTPhoisEBEtaGap->push_back(photon.isEBEtaGap());
    fTPhoisEBPhiGap->push_back(photon.isEBPhiGap());
    fTPhoisEERingGap->push_back(photon.isEERingGap());
    fTPhoisEEDeeGap->push_back(photon.isEEDeeGap());
    fTPhoisEBEEGap->push_back(photon.isEBEEGap());
    fTPhoisPFlowPhoton->push_back(photon.isPFlowPhoton());
    fTPhoisStandardPhoton->push_back(photon.isStandardPhoton());

    if (doVertexingFlag && photon.hasConversionTracks()) { // photon conversions

      reco::ConversionRefVector conversions = photon.conversions();
      if (conversions.size()<1) { std::cout << "something wrong here" << std::endl; }
      reco::ConversionRef conv = conversions[0];
      if (!conv->conversionVertex().isValid()) continue;

      fTPhoConvValidVtx->push_back(conv->conversionVertex().isValid());
      fTPhoConvChi2Probability->push_back(-999.);
      fTPhoConvNtracks->push_back(-999.);
      fTPhoConvEoverP->push_back(-999.);
      for (unsigned int i=0; i<conversions.size(); i++) {
        conv=conversions[i];
        if(ConversionsCut(*conv)) continue;
        reco::Vertex vtx=conv->conversionVertex();
        pho_conv_vtx[phoqi].SetXYZ(vtx.x(), vtx.y(), vtx.z());
        (*fTPhoConvChi2Probability)[phoqi] = ChiSquaredProbability(vtx.chi2(), vtx.ndof());
        (*fTPhoConvNtracks)[phoqi]         = conv->nTracks();
        (*fTPhoConvEoverP)[phoqi]          = conv->EoverPrefittedTracks();
        pho_conv_refitted_momentum[phoqi].SetXYZ(conv->refittedPairMomentum().x(), conv->refittedPairMomentum().y(), conv->refittedPairMomentum().z());
      }

    }
    

    if (!fIsRealData){

      edm::Handle<GenParticleCollection> gen;
      iEvent.getByLabel(fGenPartTag, gen);
      GenParticleCollection::const_iterator g_part;
      std::vector<const reco::GenParticle*> matched = matchRecoCand(&photon,iEvent);
      if (matched[0]==NULL) {
        fTPhoMCmatchexitcode->push_back(-1);
        fTPhoMCmatchindex->push_back(-999);
      }
      else if (matched[0]->pdgId()!=22) {
        fTPhoMCmatchexitcode->push_back(0);
        fTPhoMCmatchindex->push_back(-999);
      }
      else {
        fTPhoMCmatchindex->push_back(-999);
        fTPhoMCmatchexitcode->push_back(-1); // Initialize
        for(int i=0; i<*fTNGenPhotons; ++i){
          if ( (fabs((*fTGenPhotonPt)[i]-matched[0]->pt())<0.01*matched[0]->pt()) 
               && (fabs((*fTGenPhotonEta)[i]-matched[0]->eta())<0.01) 
               && ( DeltaPhi((*fTGenPhotonPhi)[i],matched[0]->phi())<0.01 ) ) {
            (*fTPhoMCmatchindex)[phoqi] = i;
          }
        }

        if ((*fTPhoMCmatchindex)[phoqi] != -999){
          if ( (*fTGenPhotonMotherID)[(*fTPhoMCmatchindex)[phoqi]]>=-6 
               && (*fTGenPhotonMotherID)[(*fTPhoMCmatchindex)[phoqi]]<=6) (*fTPhoMCmatchexitcode)[phoqi]=1;
          else if ((*fTGenPhotonMotherID)[phoqi]==21 ) (*fTPhoMCmatchexitcode)[phoqi]=1;
          else if ((*fTGenPhotonMotherID)[phoqi]==22 
                   && (*fTGenPhotonMotherStatus)[(*fTPhoMCmatchindex)[phoqi]]==3) (*fTPhoMCmatchexitcode)[phoqi]=2;
          else (*fTPhoMCmatchexitcode)[phoqi] = 3;
        }
        else (*fTPhoMCmatchexitcode)[phoqi] = 2;
      }

    }

    {
      fTPhotSCindex ->push_back(-1); // Initialize
      float diff=1e+4;
      for (int scind=0; scind<*fTNSuperClusters; scind++){
        if (fabs((*fTSCRaw)[scind]-photon.superCluster()->rawEnergy())<diff) {
          (*fTPhotSCindex)[phoqi] = scind;
          diff=fabs((*fTSCRaw)[scind]-photon.superCluster()->rawEnergy());
        }
      }

      if ((*fTPhotSCindex)[phoqi]!=-1)
        if (fabs((*fTSCEta)[(*fTPhotSCindex)[phoqi]]-photon.superCluster()->eta())>0.1 ||
            DeltaPhi((*fTSCPhi)[(*fTPhotSCindex)[phoqi]],photon.superCluster()->phi())>0.1){
          (*fTPhotSCindex)[phoqi] = -1;	    
        }

      if ( (*fTPhotSCindex)[phoqi]==-1) {
        edm::LogWarning("NTP") << "@SUB=analyze" << "No matching SC found for photon"; 
        // *fTGoodEvent = 1; 
        // break;
      }
    }


    /*
    { // start PF stuff from Nicholas

      reco::PhotonCollection::const_iterator gamIterSl;

      const Photon* gamIter = &photon;
	 
      // Initialization
      fTPhoCone04PhotonIsodR0dEta0pt0 ->push_back(0);
      fTPhoCone04PhotonIsodR0dEta0pt5 ->push_back(0);
      fTPhoCone04PhotonIsodR8dEta0pt0 ->push_back(0);
      fTPhoCone04PhotonIsodR8dEta0pt5 ->push_back(0);

      fTPhoCone01PhotonIsodR045EB070EEdEta015pt08EB1EEmvVtx ->push_back(0);
      fTPhoCone02PhotonIsodR045EB070EEdEta015pt08EB1EEmvVtx ->push_back(0);
      fTPhoCone03PhotonIsodR045EB070EEdEta015pt08EB1EEmvVtx ->push_back(0);  
      fTPhoCone04PhotonIsodR045EB070EEdEta015pt08EB1EEmvVtx ->push_back(0);

      fTPhoCone04NeutralHadronIsodR0dEta0pt0 ->push_back(0);
      fTPhoCone04NeutralHadronIsodR0dEta0pt5 ->push_back(0);
      fTPhoCone04NeutralHadronIsodR0dEta0pt0nocracks ->push_back(0);
      fTPhoCone04NeutralHadronIsodR0dEta0pt5nocracks ->push_back(0);
      fTPhoCone04NeutralHadronIsodR7dEta0pt0 ->push_back(0);
      fTPhoCone04NeutralHadronIsodR7dEta0pt5 ->push_back(0);

      fTPhoCone01NeutralHadronIsodR0dEta0pt0mvVtx ->push_back(0);
      fTPhoCone02NeutralHadronIsodR0dEta0pt0mvVtx ->push_back(0);
      fTPhoCone03NeutralHadronIsodR0dEta0pt0mvVtx ->push_back(0);
      fTPhoCone04NeutralHadronIsodR0dEta0pt0mvVtx ->push_back(0);

      fTPhoCone04ChargedHadronIsodR0dEta0pt0dz0old ->push_back(0);
      fTPhoCone04ChargedHadronIsodR0dEta0pt0PFnoPUold ->push_back(0);
      fTPhoCone04ChargedHadronIsodR015dEta0pt0dz0old ->push_back(0);
      fTPhoCone04ChargedHadronIsodR015dEta0pt0PFnoPUold ->push_back(0);

      fTPhoCone01ChargedHadronIsodR0dEta0pt0dz0 ->push_back(0);
      fTPhoCone01ChargedHadronIsodR0dEta0pt0dz1dxy01 ->push_back(0);
      fTPhoCone01ChargedHadronIsodR0dEta0pt0PFnoPU ->push_back(0);
      fTPhoCone01ChargedHadronIsodR015dEta0pt0dz0 ->push_back(0);
      fTPhoCone01ChargedHadronIsodR015dEta0pt0dz1dxy01 ->push_back(0);
      fTPhoCone01ChargedHadronIsodR015dEta0pt0PFnoPU ->push_back(0);

      fTPhoCone02ChargedHadronIsodR0dEta0pt0dz0 ->push_back(0);
      fTPhoCone02ChargedHadronIsodR0dEta0pt0dz1dxy01 ->push_back(0);
      fTPhoCone02ChargedHadronIsodR0dEta0pt0PFnoPU ->push_back(0);
      fTPhoCone02ChargedHadronIsodR015dEta0pt0dz0 ->push_back(0);
      fTPhoCone02ChargedHadronIsodR015dEta0pt0dz1dxy01 ->push_back(0);
      fTPhoCone02ChargedHadronIsodR015dEta0pt0PFnoPU ->push_back(0);

      fTPhoCone03ChargedHadronIsodR0dEta0pt0dz0 ->push_back(0);
      fTPhoCone03ChargedHadronIsodR0dEta0pt0dz1dxy01 ->push_back(0);
      fTPhoCone03ChargedHadronIsodR0dEta0pt0PFnoPU ->push_back(0);
      fTPhoCone03ChargedHadronIsodR015dEta0pt0dz0 ->push_back(0);
      fTPhoCone03ChargedHadronIsodR015dEta0pt0dz1dxy01 ->push_back(0);
      fTPhoCone03ChargedHadronIsodR015dEta0pt0PFnoPU ->push_back(0);

      fTPhoCone04ChargedHadronIsodR0dEta0pt0dz0 ->push_back(0);
      fTPhoCone04ChargedHadronIsodR0dEta0pt0dz1dxy01 ->push_back(0);
      fTPhoCone04ChargedHadronIsodR0dEta0pt0PFnoPU ->push_back(0); //pour reference
      fTPhoCone04ChargedHadronIsodR015dEta0pt0dz0 ->push_back(0);
      fTPhoCone04ChargedHadronIsodR015dEta0pt0dz1dxy01 ->push_back(0);
      fTPhoCone04ChargedHadronIsodR015dEta0pt0PFnoPU ->push_back(0);

      //Look for associated PF objects
      bool FoundPFPhoton   = false;
      bool FoundPFElectron = false;

      int iPFduplicata=-1;

      //Find PFPhoton
      int iphot=-1;
      int ncand=pfCandidates->size();
      for( int i=0; i<ncand; ++i ) {
        if ((*pfCandidates)[i].particleId()==reco::PFCandidate::gamma){
          if ((*pfCandidates)[i].mva_nothing_gamma()>0){
            if( (*pfCandidates)[i].superClusterRef()==gamIter->superCluster()) {
              iphot = i;
            }
          }
        }
      }    

      if (iphot!=-1) FoundPFPhoton=true;

      //Find PFElectron
      bool foundEgSC = false;
      int iel = -1;
      reco::GsfElectronCollection::const_iterator elIterSl;
      for (reco::GsfElectronCollection::const_iterator elIter = electronHandle->begin(); elIter != electronHandle->end(); ++elIter){
        if (gamIter->superCluster()==elIter->superCluster()) {
          elIterSl = elIter;
          foundEgSC = true;
        }

        if (foundEgSC){
          int iid = 0;
          double MVACut_ = -0.1; //42X
          //double MVACut_ = -1.; //44X
          int ncand=pfCandidates->size();
          for( int i=0; i<ncand; ++i ) {
            if ((*pfCandidates)[i].particleId()==reco::PFCandidate::e && (*pfCandidates)[i].gsfTrackRef().isNull()==false && (*pfCandidates)[i].mva_e_pi()>MVACut_ && (*pfCandidates)[i].gsfTrackRef()==elIterSl->gsfTrack()){
              iel = i;
              iid++;
            }
          }
        }
     
      }


      if (iel!=-1) FoundPFElectron=true;
   
      double photonPhi;
      double photonEta;
      double photonVz;

      double phoSC_GeomEta = gamIter->superCluster()->eta();
      double phoSC_GeomPhi = gamIter->superCluster()->phi();

      if (FoundPFPhoton){ //use PFPhoton

        fTPhoisPFPhoton   ->push_back(1);
        fTPhoisPFElectron ->push_back(0);

        if (pfPhotonHandle->size()>0){
          double Emin=4000;
          reco::PhotonCollection::const_iterator pfgamIterSl;
          for (reco::PhotonCollection::const_iterator pfgamIter = pfPhotonHandle->begin(); pfgamIter != pfPhotonHandle->end(); ++pfgamIter){
            if (fabs((*pfCandidates)[iphot].energy()-pfgamIter->energy())<Emin) {
              Emin = fabs((*pfCandidates)[iphot].energy()-pfgamIter->energy());
              pfgamIterSl=pfgamIter;       
            }
          }
          if (Emin<4000){
            fTPhoChargedHadronIso ->push_back(pfgamIterSl->chargedHadronIso());
            fTPhoNeutralHadronIso ->push_back(pfgamIterSl->neutralHadronIso());
            fTPhoPhotonIso        ->push_back(pfgamIterSl->photonIso());
          }
        }

        iPFduplicata = iphot;
        photonPhi = (*pfCandidates)[iphot].phi();
        photonEta = (*pfCandidates)[iphot].eta();
        //photonVz = (*pfCandidates)[iphot].vz();
        photonVz = gamIter->vz();
      }

      else if (!FoundPFPhoton && FoundPFElectron){ //use PFElectron


        fTPhoisPFPhoton ->push_back(0);
        fTPhoisPFElectron ->push_back(1);

        fTPhoChargedHadronIso ->push_back(elIterSl->pfIsolationVariables().chargedHadronIso);
        fTPhoNeutralHadronIso ->push_back(elIterSl->pfIsolationVariables().neutralHadronIso);
        fTPhoPhotonIso ->push_back(elIterSl->pfIsolationVariables().photonIso);

        iPFduplicata = iel;
        photonPhi = (*pfCandidates)[iel].phi();
        photonEta = (*pfCandidates)[iel].eta();
        //photonVz = (*pfCandidates)[iphot].vz();
        photonVz = gamIter->vz();
 
      }

      else if (!FoundPFPhoton && !FoundPFElectron){ //e/g only => use cones

        fTPhoisPFPhoton ->push_back(0);
        fTPhoisPFElectron ->push_back(0);


        photonPhi = phoSC_GeomPhi;
        photonEta = phoSC_GeomEta;
        //photonVz = (*pfCandidates)[iel].vz();
        photonVz = gamIter->vz();
      }



      if (FoundPFPhoton && FoundPFElectron){ //both => we have used PFPhoton

        fTPhoisPFPhoton ->push_back(1);
        fTPhoisPFElectron ->push_back(1);

      }

 
      int pho_Cone06NbPfCand;
      int pho_Cone06PfCandType[500];
      int pho_Cone06PfCandOverlap[500];
      float pho_Cone06PfCandEta[500];
      float pho_Cone06PfCandPhi[500];
      float pho_Cone06PfCandDeltaR[500];
      float pho_Cone06PfCandDeltaEta[500];
      float pho_Cone06PfCandDeltaPhi[500];
      float pho_Cone06PfCandPt[500];
      float pho_Cone06PfCandDz[500]; 
      int pho_Cone06PfCandIsFromPU[500]; 
      float pho_Cone06PfCandDxy[500];
      float pho_Cone06PfCandDeltaRrecomputed[500];
      float pho_Cone06PfCandDeltaEtarecomputed[500];
      float pho_Cone06PfCandDeltaPhirecomputed[500];
      float pho_Cone06PfCandPtrecomputed[500];

      //Recompute pflow isolation keeping all the pfcandidates in a 0.6 cone
      pho_Cone06NbPfCand = 0;
      int ipf=0;

      //if (iphot!=-1){
    
      double dR;
      int type = -1;

      for( int i=0; i<ncand; ++i ) {

        if (FoundPFPhoton && i==iphot) continue;
        if (FoundPFElectron && i==iel) continue;

        dR = reco::deltaR(photonEta,photonPhi,(*pfCandidates)[i].eta(),(*pfCandidates)[i].phi());
        if (dR<0.6 && dR>1e-05){
	   
          type = FindPFCandType((*pfCandidates)[i].pdgId());
          //cout << "type=" << type<<endl;
	   
          if (type==0 || type==1 || type==2){
	   
            if (dR<0.6){
	       
              bool isOverlapping = false;
              //isOverlapping = CheckPhotonPFCandOverlap(gamIter->superCluster(), pfCandidates, i);
	       
              pho_Cone06PfCandType ->push_back(type);
              pho_Cone06PfCandOverlap ->push_back(isOverlapping);
	       
              pho_Cone06PfCandEta ->push_back((*pfCandidates)[i].eta());
              pho_Cone06PfCandPhi ->push_back((*pfCandidates)[i].phi());
              pho_Cone06PfCandDeltaR ->push_back(reco::deltaR( photonEta,photonPhi,(*pfCandidates)[i].eta(),(*pfCandidates)[i].phi()));
              pho_Cone06PfCandDeltaEta ->push_back((*pfCandidates)[i].eta()-photonEta);
              pho_Cone06PfCandDeltaPhi ->push_back(DeltaPhi(photonPhi,(*pfCandidates)[i].phi()));
              pho_Cone06PfCandPt ->push_back((*pfCandidates)[i].pt());
              pho_Cone06PfCandDz ->push_back(fabs((*pfCandidates)[i].vz()-photonVz));
              pho_Cone06PfCandDxy ->push_back(( -((*pfCandidates)[i].vx() - gamIter->vx())*(*pfCandidates)[i].py() + ((*pfCandidates)[i].vy() - gamIter->vy())*(*pfCandidates)[i].px()) / (*pfCandidates)[i].pt());

              math::XYZVector vCand = math::XYZVector(gamIter->superCluster()->x(), gamIter->superCluster()->y(), gamIter->superCluster()->z());
              float r = vCand.R();
              math::XYZVector pfvtx((*pfCandidates)[i].vx(), (*pfCandidates)[i].vy(), (*pfCandidates)[i].vz());
              math::XYZVector pvm(((*pfCandidates)[i].momentum()*r/(*pfCandidates)[i].momentum().R()) + pfvtx);

              //float dR = deltaR(vCand.Eta(), vCand.Phi(), pvm.Eta(), pvm.Phi());
              //float dEta = fabs(vCand.Eta() - pvm.Eta());
              //double dPhi = fabs(vCand.Phi() - pvm.Phi());

              pho_Cone06PfCandDeltaRrecomputed ->push_back(reco::deltaR( vCand.Eta(),vCand.Phi(),pvm.Eta(),pvm.Phi()));
              pho_Cone06PfCandDeltaEtarecomputed ->push_back(pvm.Eta() - vCand.Eta());
              pho_Cone06PfCandDeltaPhirecomputed ->push_back(DeltaPhi(vCand.Phi(),pvm.Phi()));
              pho_Cone06PfCandPtrecomputed ->push_back((*pfCandidates)[i].pt());	       


              pho_Cone06PfCandIsFromPU ->push_back(-1);
              if (type==1){
                //cout << "A" << endl;
                reco::VertexRef chvtx = chargedHadronVertex(alternativeVertexHandle, (*pfCandidates)[i]);
                //cout << "B"<<endl;
                if (chvtx.isNull() || chvtx.key()==0) pho_Cone06PfCandIsFromPU ->push_back(0);
                else pho_Cone06PfCandIsFromPU ->push_back(1);
                //cout << "C"<<endl;
              }
              ipf++;
	       
              if (dR<0.4 && isOverlapping==false) {
                FillPhotonIsoVariables(photonEta, photonPhi, photonVz, type, pho_Cone06PfCandIsFromPU[ipf], pfCandidates, i, phoqi);
                { //		 FillPhotonIsoVariables_Frixione_Neutrals(type, ipf, phoqi);


                  double pt = pho_Cone06PfCandPt[ipf];
                  double dEta = pho_Cone06PfCandDeltaEtarecomputed[ipf];
                  double dPhi = pho_Cone06PfCandDeltaPhirecomputed[ipf];
                  double dR = sqrt(dEta*dEta+dPhi*dPhi);
  


                  if (type==0){ //Neutral Hadron
    
                    if (dR<0.1) fTPhoCone01NeutralHadronIso_dR0_dEta0_pt0_mvVtx[phoqi] += pt;
                    if (dR<0.2) fTPhoCone02NeutralHadronIso_dR0_dEta0_pt0_mvVtx[phoqi] += pt;
                    if (dR<0.3) fTPhoCone03NeutralHadronIso_dR0_dEta0_pt0_mvVtx[phoqi] += pt;
                    if (dR<0.4) fTPhoCone04NeutralHadronIso_dR0_dEta0_pt0_mvVtx[phoqi] += pt;
                  }

                  if (type==2) { //Photon

                    if ((dR>0.045 && gamIter->isEB() && pt>0.08)||(dR>0.07 && gamIter->isEE() && pt>0.1)){
                      if (dR<0.1) fTPhoCone01PhotonIso_dR045EB070EE_dEta015_pt08EB1EE_mvVtx[phoqi] += pt;
                      if (dR<0.2) fTPhoCone02PhotonIso_dR045EB070EE_dEta015_pt08EB1EE_mvVtx[phoqi] += pt;
                      if (dR<0.3) fTPhoCone03PhotonIso_dR045EB070EE_dEta015_pt08EB1EE_mvVtx[phoqi] += pt;  
                      if (dR<0.4) fTPhoCone04PhotonIso_dR045EB070EE_dEta015_pt08EB1EE_mvVtx[phoqi] += pt;
                    }
                  }
                }

                { // FillPhotonIsoVariables_Frixione_ChHad(type, pho_Cone06PfCandIsFromPU[ipf], ipf, phoqi);

                  bool isPU=pho_Cone06PfCandIsFromPU[ipf];

                  double pt = pho_Cone06PfCandPt[ipf];
                  double dEta = pho_Cone06PfCandDeltaEtarecomputed[ipf];
                  double dPhi = pho_Cone06PfCandDeltaPhirecomputed[ipf];
                  double dR = sqrt(dEta*dEta+dPhi*dPhi);
                  double dz = pho_Cone06PfCandDz[ipf];
                  double dxy = pho_Cone06PfCandDxy[ipf];
  
                  if (type==1){ //Charged Hadron
    
                    //no cut
                    if (dR<0.1) fTPhoCone01ChargedHadronIso_dR0_dEta0_pt0_dz0[phoqi] += pt;
                    if (dR<0.2) fTPhoCone02ChargedHadronIso_dR0_dEta0_pt0_dz0[phoqi] += pt;
                    if (dR<0.3) fTPhoCone03ChargedHadronIso_dR0_dEta0_pt0_dz0[phoqi] += pt;
                    if (dR<0.4) fTPhoCone04ChargedHadronIso_dR0_dEta0_pt0_dz0[phoqi] += pt;

                    if (dR>0.015){
                      if (dR<0.1) fTPhoCone01ChargedHadronIso_dR015_dEta0_pt0_dz0[phoqi] += pt;
                      if (dR<0.2) fTPhoCone02ChargedHadronIso_dR015_dEta0_pt0_dz0[phoqi] += pt;
                      if (dR<0.3) fTPhoCone03ChargedHadronIso_dR015_dEta0_pt0_dz0[phoqi] += pt;
                      if (dR<0.4) fTPhoCone04ChargedHadronIso_dR015_dEta0_pt0_dz0[phoqi] += pt;
                    }
    
                    //dz/dxy
                    if (dz<1. && dxy<0.1){

                      if (dR<0.1) fTPhoCone01ChargedHadronIso_dR0_dEta0_pt0_dz1_dxy01[phoqi] += pt;
                      if (dR<0.2) fTPhoCone02ChargedHadronIso_dR0_dEta0_pt0_dz1_dxy01[phoqi] += pt;
                      if (dR<0.3) fTPhoCone03ChargedHadronIso_dR0_dEta0_pt0_dz1_dxy01[phoqi] += pt;
                      if (dR<0.4) fTPhoCone04ChargedHadronIso_dR0_dEta0_pt0_dz1_dxy01[phoqi] += pt;

                      if (dR>0.015){
                        if (dR<0.1) fTPhoCone01ChargedHadronIso_dR015_dEta0_pt0_dz1_dxy01[phoqi] += pt;
                        if (dR<0.2) fTPhoCone02ChargedHadronIso_dR015_dEta0_pt0_dz1_dxy01[phoqi] += pt;
                        if (dR<0.3) fTPhoCone03ChargedHadronIso_dR015_dEta0_pt0_dz1_dxy01[phoqi] += pt;
                        if (dR<0.4) fTPhoCone04ChargedHadronIso_dR015_dEta0_pt0_dz1_dxy01[phoqi] += pt;
                      }

                    }

                    //pfNoPU
                    if (isPU==false){

                      if (dR<0.1) fTPhoCone01ChargedHadronIso_dR0_dEta0_pt0_PFnoPU[phoqi] += pt;
                      if (dR<0.2) fTPhoCone02ChargedHadronIso_dR0_dEta0_pt0_PFnoPU[phoqi] += pt;
                      if (dR<0.3) fTPhoCone03ChargedHadronIso_dR0_dEta0_pt0_PFnoPU[phoqi] += pt;
                      if (dR<0.4) fTPhoCone04ChargedHadronIso_dR0_dEta0_pt0_PFnoPU[phoqi] += pt;

                      if (dR>0.015){
                        if (dR<0.1) fTPhoCone01ChargedHadronIso_dR015_dEta0_pt0_PFnoPU[phoqi] += pt;
                        if (dR<0.2) fTPhoCone02ChargedHadronIso_dR015_dEta0_pt0_PFnoPU[phoqi] += pt;
                        if (dR<0.3) fTPhoCone03ChargedHadronIso_dR015_dEta0_pt0_PFnoPU[phoqi] += pt;
                        if (dR<0.4) fTPhoCone04ChargedHadronIso_dR015_dEta0_pt0_PFnoPU[phoqi] += pt;
                      }
                    }

                  }


                }



		

              }
            }
          }
	   
        }

      }
      pho_Cone06NbPfCand = ipf;



    }     // end PF stuff from Nicholas
    */


		


    // DISABLED: NO SEED IN AOD (UPDATE IT IN 4_2)
    // 	// Spike removal information
    // 	if ( photon.superCluster()->seed()->caloID().detector( reco::CaloID::DET_ECAL_BARREL ) ) {
    // 		fTPhotScSeedSeverity ->push_back(EcalSeverityLevelAlgo::severityLevel( photon.superCluster()->seed()->seed(), *ebRecHits, *channelStatus ));
    // 		fTPhotE1OverE9 ->push_back(EcalSeverityLevelAlgo::E1OverE9(   photon.superCluster()->seed()->seed(), *ebRecHits ));
    // 		fTPhotS4OverS1 ->push_back(EcalSeverityLevelAlgo::swissCross( photon.superCluster()->seed()->seed(), *ebRecHits ));
    // 	} else if ( photon.superCluster()->seed()->caloID().detector( reco::CaloID::DET_ECAL_ENDCAP ) ) {
    // 		fTPhotScSeedSeverity ->push_back(EcalSeverityLevelAlgo::severityLevel( photon.superCluster()->seed()->seed(), *eeRecHits, *channelStatus ));
    // 		fTPhotE1OverE9 ->push_back(EcalSeverityLevelAlgo::E1OverE9(   photon.superCluster()->seed()->seed(), *eeRecHits ));
    // 		fTPhotS4OverS1 ->push_back(1.0-EcalSeverityLevelAlgo::swissCross( photon.superCluster()->seed()->seed(), *eeRecHits ));
    //      } else
    // 			edm::LogWarning("NTP") << "Photon supercluster seed crystal neither in EB nor in EE!";
  } // end photons

  /*
    USAGE OF VERTEX CHOICE FOR DIPHOTON EVENTS:
    
    diphotons_{first,second} are vectors of {photon_1_index,photon_2_index}
    
    vtx_dipho_??? are, for each diphoton pair, vectors of vertex indices (as ranked by the different algos)
    
    For example: best vertex for diphoton pair 3, with photon_1_index=diphotons_first[3] and photon_2_index=diphotons_second[3]: vtx_dipho_bla[3].at(0), second choice vtx_dipho_bla[3].at(1) ...
    
  */

  std::vector<int> diphotons_first;
  std::vector<int> diphotons_second;
  std::vector<std::vector<int> > vtx_dipho_h2gglobe;
  std::vector<std::vector<int> > vtx_dipho_mva;
  std::vector<std::vector<int> > vtx_dipho_productrank;
  
  if (doVertexingFlag) { // start vertex selection stuff with MVA from Hgg (Musella) UserCode/HiggsAnalysis/HiggsTo2photons/h2gglobe/VertexAnalysis tag vertex_mva_v4
    
    bool VTX_MVA_DEBUG = false;
    
    edm::Handle<reco::TrackCollection> tkH;
    iEvent.getByLabel(fTrackCollForVertexing, tkH);
    
    edm::Handle<VertexCollection> vtxH = vertices;
    
    edm::Handle<reco::ConversionCollection> convH;
    iEvent.getByLabel(fAllConversionsCollForVertexing, convH);
    
    int tk_n = 0; 
    
    float tk_px[__TRK_AUX_ARRAYS_DIM__];
    float tk_py[__TRK_AUX_ARRAYS_DIM__];
    float tk_pz[__TRK_AUX_ARRAYS_DIM__];
    float tk_d0[__TRK_AUX_ARRAYS_DIM__];
    float tk_dz[__TRK_AUX_ARRAYS_DIM__];
    float tk_d0err[__TRK_AUX_ARRAYS_DIM__];
    float tk_dzerr[__TRK_AUX_ARRAYS_DIM__];
    float tk_pterr[__TRK_AUX_ARRAYS_DIM__];
    bool tk_ishighpurity[__TRK_AUX_ARRAYS_DIM__];
    std::vector<std::vector<unsigned short> > vtx_std_tkind;
    std::vector<std::vector<float> > vtx_std_tkweight;
    int vtx_std_ntks[__VTX_AUX_ARRAYS_DIM__];
    int tkVtxId[__TRK_AUX_ARRAYS_DIM__];
    
    
    { // tracks
      if (VTX_MVA_DEBUG)	   	   cout << "tracks begin" << endl;
      std::vector<reco::TrackBaseRef>::const_iterator tk;
          
      for(unsigned int i=0; i<vtxH->size(); i++) {

        if (VTX_MVA_DEBUG) cout << "working on vtx " << i << endl;

        if (vtxH->size()>__VTX_AUX_ARRAYS_DIM__) std::cout << "Too many vertices in the event; please enlarge the value of __VTX_AUX_ARRAYS_DIM__" << std::endl;

        reco::VertexRef vtx(vtxH, i);
	  
        vtx_std_ntks[i]=vtx->tracksSize();
        if (VTX_MVA_DEBUG)	     	     cout << "vtx tracks " << vtx->tracksSize() << endl;
	     
        std::vector<unsigned short> temp;
        std::vector<float> temp_float;

        if (vtx->tracksSize()>0){
          for(tk=vtx->tracks_begin();tk!=vtx->tracks_end();++tk) {
            if (VTX_MVA_DEBUG)		 		 cout << "processing vtx track (out of " << vtx->tracksSize() << ")" << endl;
            int index = 0;
            bool ismatched = false; 
            for(reco::TrackCollection::size_type j = 0; j<tkH->size(); ++j) {
              if (VTX_MVA_DEBUG)		   		   std::cout << j << std::endl;
              reco::TrackRef track(tkH, j);
              if(TrackCut(track)) continue; 
              if (&(**tk) == &(*track)) {
                temp.push_back(index);
                temp_float.push_back(vtx->trackWeight(track));
                ismatched = true;
                if (VTX_MVA_DEBUG)		     		     cout << "matching found index" << index << " weight " << vtx->trackWeight(track) << endl;
                break;
              }
              index++;
            }
            if(!ismatched) {
              temp.push_back(-9999);
              temp_float.push_back(-9999);
            }
          }
        }
        else {
          if (VTX_MVA_DEBUG)	       	       cout << "no vertex tracks found" << endl;
          temp = std::vector<unsigned short>(0);
          temp_float = std::vector<float>(0);
        }
      
        //	      for (std::vector<unsigned short>::const_iterator it=temp.begin(); it!=temp.end(); it++) {int k=0; vtx_std_tkind[i][k]=*it; k++;}
        //	      for (std::vector<float>::const_iterator it=temp_float.begin(); it!=temp_float.end(); it++) {int k=0; vtx_std_tkweight[i][k]=*it; k++;}	    
        vtx_std_tkind.push_back(temp);
        vtx_std_tkweight.push_back(temp_float);
        if (VTX_MVA_DEBUG)	     	     	     std::cout << "tracks: " <<  temp.size() << std::endl;
      }	  

      if (VTX_MVA_DEBUG){
        std::cout << "tkWeight is " << std::endl;
        for (int a=0; a<(int)(vtx_std_tkind.size()); a++) std::cout << a << ":" << vtx_std_tkind.at(a).size() << " " ;
        std::cout << std::endl;
      }


      for(unsigned int i=0; i<tkH->size(); i++) {

        if (tkH->size()>__TRK_AUX_ARRAYS_DIM__) std::cout << "Too many tracks in the event; please enlarge the value of __TRK_AUX_ARRAYS_DIM__" << std::endl;

        reco::TrackRef tk(tkH, i);

        if(TrackCut(tk))continue; 
	
        tk_px[tk_n] = tk->px();
        tk_py[tk_n] = tk->py();
        tk_pz[tk_n] = tk->pz();
        tk_d0[tk_n] = tk->d0();
        tk_dz[tk_n] = tk->dz();
        tk_dzerr[tk_n] = tk->dzError();
        tk_d0err[tk_n] = tk->d0Error();   
        tk_pterr[tk_n] = tk->ptError();
        tk_ishighpurity[tk_n] = tkIsHighPurity(tk);
        tkVtxId[tk_n] = -1; 
	 

        tk_n++;
      } // for i (loop over all tracks)

      if (VTX_MVA_DEBUG)	  cout << "done tracks" << endl;
    }

    (*fTNconv)=0;

    { // all conversions

	  
      for( reco::ConversionCollection::const_iterator  iConv = convH->begin(); iConv != convH->end(); iConv++) {

        reco::Conversion localConv = reco::Conversion(*iConv);
  
        if(ConversionsCut(localConv)) continue;
  
        if ( !localConv.conversionVertex().isValid() ) { (*fTNconv)++; continue; } // FR: TO BE CHECKED!
  
        reco::Vertex vtx=localConv.conversionVertex();
        conv_vtx[*fTNconv].SetXYZ(vtx.x(), vtx.y(), vtx.z());
        fTConvValidVtx          ->push_back( localConv.conversionVertex().isValid() );
        fTConvNtracks           ->push_back(localConv.nTracks());
        fTConvChi2Probability   ->push_back(ChiSquaredProbability(vtx.chi2(), vtx.ndof()));
        fTConvEoverP            ->push_back(localConv.EoverPrefittedTracks());
        fTConvZofPrimVtxFromTrks->push_back(localConv.zOfPrimaryVertexFromTracks());

        conv_refitted_momentum[*fTNconv].SetXYZ(localConv.refittedPairMomentum().x(), localConv.refittedPairMomentum().y(), localConv.refittedPairMomentum().z());
        (*fTNconv)++;

      }

    }

    if (VTX_MVA_DEBUG)	  cout << "done convs" << endl;

    VertexAlgoParameters vtxAlgoParams;
    vector<string> rankVariables;
    HggVertexAnalyzer vAna(vtxAlgoParams,(*fTNVrtx)); 
    HggVertexFromConversions vConv(vtxAlgoParams);
	 
	
      

    //	 std::string perVtxMvaWeights, perVtxMvaMethod;
	 
    //	 std::string perEvtMvaWeights, perEvtMvaMethod;
    TMVA::Reader * perVtxReader;	 
    TMVA::Reader * perEvtReader;
    vAna.setupWithDefaultOptions(perVtxMvaWeights, perEvtMvaWeights, rankVariables, perVtxReader, perVtxMvaMethod, perEvtReader, perEvtMvaMethod);
    std::vector<std::string> vtxVarNames;
    vtxVarNames.push_back("ptbal"), vtxVarNames.push_back("ptasym"), vtxVarNames.push_back("logsumpt2");

    if (VTX_MVA_DEBUG)	 cout << "ready: remember delete readers" << endl;	 


    if( (*fTNPhotons) < 2 ) {

      vtx_dipho_h2gglobe.push_back( std::vector<int>() );
      vtx_dipho_mva.push_back( std::vector<int>() );
      vtx_dipho_productrank.push_back( std::vector<int>() );
	   
      for(int ii=0;ii<(*fTNVrtx); ++ii) {vtx_dipho_h2gglobe.back().push_back(ii); }
      for(int ii=0;ii<(*fTNVrtx); ++ii) {vtx_dipho_mva.back().push_back(ii); }
      for(int ii=0;ii<(*fTNVrtx); ++ii) {vtx_dipho_productrank.back().push_back(ii); }

    } else {

      if (VTX_MVA_DEBUG)	   cout << "temp" << endl;

      // fully combinatorial vertex selection
      for(int ip=0; ip<(*fTNPhotons); ++ip) {
        for(int jp=ip+1; jp<(*fTNPhotons); ++jp) {
          diphotons_first.push_back(ip);
          diphotons_second.push_back(jp);
        }
      }


      for(unsigned int id=0; id<diphotons_first.size(); ++id ) {
			
        if (VTX_MVA_DEBUG)	     cout << "processing diphoton pair " << id << endl;

        int ipho1 = diphotons_first[id];
        int ipho2 = diphotons_second[id];

		  
        PhotonInfo pho1=fillPhotonInfos(ipho1,vtxAlgoParams.useAllConversions);
        PhotonInfo pho2=fillPhotonInfos(ipho2,vtxAlgoParams.useAllConversions);
	     
        ETHVertexInfo vinfo(int((*fTNVrtx)),
                            (*fTVrtxX),
                            (*fTVrtxY),
                            (*fTVrtxZ),
                            int(tk_n),
                            tk_px+0,
                            tk_py+0,
                            tk_pz+0,
                            tk_pterr+0,
                            tkVtxId+0,
                            tk_d0+0,
                            tk_d0err+0,
                            tk_dz+0,
                            tk_dzerr+0,
                            tk_ishighpurity+0,
                            vtx_std_tkind,
                            vtx_std_tkweight,
                            vtx_std_ntks+0
                            );


        if (VTX_MVA_DEBUG)	     cout << "filled photon/tuplevertex info" << endl;
        if (VTX_MVA_DEBUG)	     cout << vinfo.nvtx() << " vertices" << endl;

        vAna.analyze(vinfo,pho1,pho2);

        if (VTX_MVA_DEBUG)	     cout << "initialized vAna" << endl;

        // make sure that vertex analysis indexes are in synch 
        assert( int(id) == vAna.pairID(ipho1,ipho2) );

        if (VTX_MVA_DEBUG)	     cout << "starting rankings" << endl;

        if (VTX_MVA_DEBUG)	     cout << "rankprod" << endl;
        /// rank product vertex selection. Including pre-selection based on conversions information.
        vtx_dipho_productrank.push_back(vAna.rankprod(rankVariables));


        if (VTX_MVA_DEBUG)	     cout << "mva pasquale" << endl;
        /// MVA vertex selection
        vtx_dipho_mva.push_back(vAna.rank(*perVtxReader,perVtxMvaMethod));


        // vertex probability through per-event MVA (not used so far)
        // float vtxEvtMva = vAna.perEventMva( *perEvtReader,  perEvtMvaMethod, vtx_dipho_mva->back() );
        // float vtxProb = vAna.vertexProbability( vtxEvtMva );
	     
        if (VTX_MVA_DEBUG)	     cout << "mva hgg globe" << endl;
        // Globe vertex selection with conversions
        vtx_dipho_h2gglobe.push_back(HggVertexSelection(vAna, vConv, pho1, pho2, vtxVarNames,false,0,""));


        if (VTX_MVA_DEBUG){
          cout << "ranking : ";
          for (std::vector<int>::const_iterator it=vtx_dipho_productrank.at(id).begin(); it!=vtx_dipho_productrank.at(id).end(); it++) cout << *it;
          cout << endl;
          cout << "pasquale : ";
          for (std::vector<int>::const_iterator it=vtx_dipho_mva.at(id).begin(); it!=vtx_dipho_mva.at(id).end(); it++) cout << *it;
          cout << endl;
          cout << "globe : ";
          for (std::vector<int>::const_iterator it=vtx_dipho_h2gglobe.at(id).begin(); it!=vtx_dipho_h2gglobe.at(id).end(); it++) cout << *it;
          cout << endl;
        }
	     

      } // end diphoton loop
 

    } // end else
	

    if (VTX_MVA_DEBUG)	 cout << "deleting reader" << endl;
    delete perVtxReader;
    delete perEvtReader;

  } // end vertex selection for diphoton events



  //       cout << "end vertex selection MVA" << endl;



      //   ////////////////////////////////////////////////////////
      //   // Jet Variables:
      //   const JetCorrector* jetCorr = JetCorrector::getJetCorrector(fJetCorrs, iSetup);
      //   std::vector<OrderPair> corrIndices;  // Vector of indices and pt of corr. jets to re-order them
      //   int iraw(0);
      //   fTnjetstot = jets->size();
      //   // Loop over uncorr. jets
      //   for(View<Jet>::const_iterator Jit = jets->begin(); Jit != jets->end(); ++Jit, ++iraw){
      //     // Cut on uncorrected pT (for startup)
      //     if(Jit->pt() < fMinrawjpt) continue;
      //     // Save only the gMaxNJets first uncorrected jets
      //     if(iraw >= gMaxNJets){
      //       edm::LogWarning("NTP") << "@SUB=analyze"
      //                              << "Found more than " << static_cast<int>(gMaxNJets) << " uncorrected jets, I'm scared ...";
      //       fTflagmaxujetexc = 1;
      //       *fTGoodEvent = 1;
      //       break;
      //     }
                
      //     JetBaseRef jetRef(edm::Ref<JetView>(jets,iraw));
      //     double scale = jetCorr->correction(*Jit,jetRef,iEvent,iSetup);
      //     corrIndices.push_back(make_pair(iraw, scale*Jit->pt()));
      //   }
	
      //   // Sort corrected jet collection by decreasing pt
      //   std::sort(corrIndices.begin(), corrIndices.end(), indexComparator);
	
      //   // Determine corrected jets
      //   int jqi(-1); // counts # of qualified jets
      //   // Loop over corr. jet indices
      //   for(std::vector<OrderPair>::const_iterator it = corrIndices.begin(); it != corrIndices.end(); ++it ) {
      //     // Check if maximum number of jets is exceeded already
      //     if(jqi >= gMaxNJets-1) {
      //       edm::LogWarning("NTP") << "@SUB=analyze"
      //                              << "Maximum number of jets exceeded";
      //       fTflagmaxjetexc = 1;
      //       *fTGoodEvent = 1;
      //       break;
      //     }
      //     int index = it->first;
      //     const PFJet* cojet = static_cast<const PFJet*>( &((*jets)[index]) ); // look away...
      //     std::auto_ptr<PFJet> jet(new PFJet(*cojet));

      //     // The correction was calculated above: use it
      //     double scale = it->second/jet->pt();
      //     jet->scaleEnergy(scale);
	
      //     // Jet preselection
      //     if(jet->pt() < fMincorjpt) continue; 
      //     if(fabs(jet->eta()) > fMaxjeta) continue;
      //     jqi++;

      //     // Dump jet properties into tree variables
      //     fTjpx    [jqi] = jet->px();
      //     fTjpy    [jqi] = jet->py();
      //     fTjpz    [jqi] = jet->pz();
      //     fTjpt    [jqi] = jet->pt();
      //     fTjeta   [jqi] = jet->eta();
      //     fTjphi   [jqi] = jet->phi();
      //     fTje     [jqi] = jet->energy();
      //     fTjet    [jqi] = jet->et();
      //     fTjEcorr [jqi] = scale;
      //     fTJEtaRms[jqi] = sqrt(jet->etaetaMoment());
      //     fTJPhiRms[jqi] = sqrt(jet->etaphiMoment());
      //     fTjArea  [jqi] = jet->jetArea();

      //     fTjNconstituents[jqi] = jet->nConstituents();
      //     fTjChMult       [jqi] = jet->chargedMultiplicity(); // do it the pf way...
      //     fTjNeuMult      [jqi] = jet->neutralMultiplicity(); 
		
      //     // energy fractions for JID need to be computed w.r.t. uncorrected jet energy!!
      //     // see for instance https://twiki.cern.ch/twiki/bin/view/CMS/JetID
      //     // or http://cmssdt.cern.ch/SDT/doxygen/CMSSW_4_1_3/doc/html/dc/dd5/classPFJetIDSelectionFunctor.html  
      //     double uncorr_energy  = jet->energy()/scale;
      //     fTjChHadFrac    [jqi] = jet->chargedHadronEnergy()/uncorr_energy;
      //     fTjNeuHadFrac   [jqi] = jet->neutralHadronEnergy()/uncorr_energy + jet->HFHadronEnergy()/uncorr_energy;
      //     fTjChEmFrac     [jqi] = jet->chargedEmEnergy()/uncorr_energy;
      //     fTjNeuEmFrac    [jqi] = jet->neutralEmEnergy()/uncorr_energy;
      //     fTjChMuEFrac    [jqi] = jet->chargedMuEnergy()/uncorr_energy;
      // 		fTjPhoFrac       [jqi] = jet->photonEnergy()/uncorr_energy; // photons also count for neutralEmEnergy
      // 		fTjHFHadFrac     [jqi] = jet->HFHadronEnergy()/uncorr_energy;
      // 		fTjHFEMFrac      [jqi] = jet->HFEMEnergy()/uncorr_energy;   // also contained in neutralEmEnergy
      // 		// see CMSSW/RecoJets/JetProducers/src/JetSpecific.cc

      // 		vector<PFCandidatePtr> pfCandidates = jet->getPFConstituents();
		
      // 		float sumPt_cands=0.;
      // 		float sumPt2_cands=0.;
      // 		float rms_cands=0.;
		
      // 		TLorentzVector jetp4;
      // 		jetp4.SetPtEtaPhiE(jet->pt(), jet->eta(), jet->phi(), jet->energy());

      // 		for (vector<PFCandidatePtr>::const_iterator jCand = pfCandidates.begin(); jCand != pfCandidates.end(); ++jCand) {
		  
      // 		  math::XYZTLorentzVectorD const& pCand_t = (*jCand)->p4();
      // 		  TLorentzVector pCand(pCand_t.px(), pCand_t.py(), pCand_t.pz(), pCand_t.energy());
      // 		  if(pCand.Pt()>0.){
      // 		    sumPt_cands += pCand.Pt();
      // 		    sumPt2_cands += (pCand.Pt()*pCand.Pt());

      // 		    float deltaR = pCand.DeltaR(jetp4);
      // 		    rms_cands += (pCand.Pt()*pCand.Pt()*deltaR*deltaR);

      // 		  }

      // 		} //for PFCandidates

      // 		fTjPtD      [jqi] = sqrt( sumPt2_cands )/sumPt_cands;
      // 		fTjRMSCand  [jqi] = rms_cands/sumPt2_cands;



      //     // Calculate the DR wrt the closest electron
      //     float ejDRmin = 10.; // Default when no electrons previously selected
      //     for( int j = 0; j < fTneles; j++ ){
      //       float ejDR = reco::deltaR(jet->eta(), jet->phi(), fTeeta[j], fTephi[j]);
      //       if(ejDR<ejDRmin) ejDRmin = ejDR;
      //     }
      //     fTjeMinDR[jqi] = ejDRmin;

      // 		// B-tagging probability (for 4 b-taggings)
      // 		// remember: 'index' is the index of the uncorrected jet, as saved in the event
      // 		fTjbTagProbTkCntHighEff[jqi]  = (*jetsAndProbsTkCntHighEff) [index].second;
      // 		fTjbTagProbTkCntHighPur[jqi]  = (*jetsAndProbsTkCntHighPur) [index].second;
      // 		fTjbTagProbSimpSVHighEff[jqi] = (*jetsAndProbsSimpSVHighEff)[index].second;
      // 		fTjbTagProbSimpSVHighPur[jqi] = (*jetsAndProbsSimpSVHighPur)[index].second;

      //     // Jet-track association: get associated tracks
      //     const reco::TrackRefVector& tracks = jet->getTrackRefs();
      //     std::vector<const reco::Track*> AssociatedTracks;
      //     for( TrackRefVector::iterator it = tracks.begin(); it != tracks.end(); ++it ) AssociatedTracks.push_back( it->get() );
			
      //     // Below save the momenta of the three leading tracks associated to the jet
      //     float pT1(0.), pT2(0.), pT3(0.);
      //     int idx1(-1), idx2(-1), idx3(-1);
			
      //     // Jet-track association: make transient tracks and store information
      //     std::vector<TransientTrack> AssociatedTTracks;
      //     fTjMass[jqi] = 0.;
      //     if(fabs(jet->eta())<2.9){ // when the cone of dR=0.5 around the jet is (at least partially) inside the tracker acceptance
      //       // Tmp variables for vectorial sum of pt of tracks
      //       double pXtmp(0.), pYtmp(0.), pZtmp(0.), E2tmp(0.);
      //       const double trkmass = 0.; // Assumed mass for tracks

      //       // Loop over associated tracks:
      //       for(size_t t = 0; t < AssociatedTracks.size(); ++t){
      //         AssociatedTTracks.push_back(theB->build(AssociatedTracks[t])); // build transient tracks for vertex fitting below
      //         if(AssociatedTracks[t]->normalizedChi2()<10. && AssociatedTracks[t]->numberOfValidHits()>10 && AssociatedTracks[t]->pt()>1.){
      //           pXtmp += AssociatedTracks[t]->px();
      //           pYtmp += AssociatedTracks[t]->py();
      //           pZtmp += AssociatedTracks[t]->pz();
      //           E2tmp += trkmass*trkmass + pXtmp*pXtmp + pYtmp*pYtmp + pZtmp*pZtmp;
      //         }
      //         // Find the three highest pT tracks
      //         if(AssociatedTracks[t]->pt() > pT1 && AssociatedTracks.size() >= 1){
      //           pT1=AssociatedTracks[t]->pt();
      //           idx3=idx2;
      //           idx2=idx1;
      //           idx1=t;
      //         } else if (AssociatedTracks[t]->pt() < pT1 && AssociatedTracks[t]->pt() > pT2 && AssociatedTracks.size() >= 2) {
      //           pT2=AssociatedTracks[t]->pt();
      //           idx3=idx2;
      //           idx2=t;
      //         } else if (AssociatedTracks[t]->pt() < pT2 && AssociatedTracks[t]->pt() > pT3 && AssociatedTracks.size() >= 3){
      //           pT3=AssociatedTracks[t]->pt();
      //           idx3=t;
      //         }
      //       }
      //       // Fill the momenta
      //       if(AssociatedTracks.size()>=1){
      //         fTjtrk1px[jqi] = AssociatedTracks[idx1]->px();
      //         fTjtrk1py[jqi] = AssociatedTracks[idx1]->py();
      //         fTjtrk1pz[jqi] = AssociatedTracks[idx1]->pz();
      //       }
      //       if(AssociatedTracks.size()>=2){
      //         fTjtrk2px[jqi] = AssociatedTracks[idx2]->px();
      //         fTjtrk2py[jqi] = AssociatedTracks[idx2]->py();
      //         fTjtrk2pz[jqi] = AssociatedTracks[idx2]->pz();
      //       }
      //       if(AssociatedTracks.size()>=3){
      //         fTjtrk3px[jqi] = AssociatedTracks[idx3]->px();
      //         fTjtrk3py[jqi] = AssociatedTracks[idx3]->py();
      //         fTjtrk3pz[jqi] = AssociatedTracks[idx3]->pz();
      //       }
			
      //       fTjMass[jqi]   = sqrt(E2tmp - pXtmp*pXtmp - pYtmp*pYtmp - pZtmp*pZtmp);
      //       if(fTjChMult[jqi] > 0) fTjMass[jqi] *= fTjNconstituents[jqi]/fTjChMult[jqi]; // apparantly there ARE cases where ChMult is 0, but we still end up here...
      //       else fTjMass[jqi] = 0.;
      //     } else { // The whole cone used for jet-tracks association is outside of the tracker acceptance
      //       fTjMass[jqi] = -888.88;
      //     }

      //     // Do a vertex fitting with the tracks
      //     if(AssociatedTTracks.size() > 1) {
      //       TransientVertex jetVtx = avFitter.vertex(AssociatedTTracks);
      //       if(jetVtx.isValid()){
      //         fTjetVtxx    [jqi] = jetVtx.position().x();
      //         fTjetVtxy    [jqi] = jetVtx.position().y();
      //         fTjetVtxz    [jqi] = jetVtx.position().z();
      //         fTjetVtxExx  [jqi] = jetVtx.positionError().cxx();
      //         fTjetVtxEyx  [jqi] = jetVtx.positionError().cyx();
      //         fTjetVtxEyy  [jqi] = jetVtx.positionError().cyy();
      //         fTjetVtxEzy  [jqi] = jetVtx.positionError().czy();
      //         fTjetVtxEzz  [jqi] = jetVtx.positionError().czz();
      //         fTjetVtxEzx  [jqi] = jetVtx.positionError().czx();
      //         fTjetVtxNChi2[jqi] = jetVtx.normalisedChiSquared();
      //       }else{
      //         fTjetVtxx    [jqi] = -777.77;
      //         fTjetVtxy    [jqi] = -777.77;
      //         fTjetVtxz    [jqi] = -777.77;
      //         fTjetVtxExx  [jqi] = -777.77;
      //         fTjetVtxEyx  [jqi] = -777.77;
      //         fTjetVtxEyy  [jqi] = -777.77;
      //         fTjetVtxEzy  [jqi] = -777.77;
      //         fTjetVtxEzz  [jqi] = -777.77;
      //         fTjetVtxEzx  [jqi] = -777.77;
      //         fTjetVtxNChi2[jqi] = -777.77;
      //       }
      //     }else{
      //       fTjetVtxx    [jqi] = -888.88;
      //       fTjetVtxy    [jqi] = -888.88;
      //       fTjetVtxz    [jqi] = -888.88;
      //       fTjetVtxExx  [jqi] = -888.88;
      //       fTjetVtxEyx  [jqi] = -888.88;
      //       fTjetVtxEyy  [jqi] = -888.88;
      //       fTjetVtxEzy  [jqi] = -888.88;
      //       fTjetVtxEzz  [jqi] = -888.88;
      //       fTjetVtxEzx  [jqi] = -888.88;
      //       fTjetVtxNChi2[jqi] = -888.88;
      //     }
      //     AssociatedTracks.clear();
      //     AssociatedTTracks.clear();
	
      //     // GenJet matching
      //     if (!fIsRealData && fTNGenJets > 0) fTjetGenJetIndex[jqi] = matchJet(&(*jet));
      //     fTgoodjet[jqi] = 0;
      //   }
      //   fTnjets = jqi+1;
      //   // corrJets.clear();
      //   corrIndices.clear();

      //   // Check electron duplication
      //   ElectronDuplicate(elecPtr, trckPtr);
      //   // Check photon/electron duplication
      //   PhotonElectronDuplicate(elecPtr, photSCs);

      //   // These don't work for pf jets yet
      //   // Check electron/jet duplication
      //   // ElJetOverlap(jetPtr, elecPtr, calotowers);
      //   // Check photon/jet duplication
      //   // PhotonJetOverlap(jetPtr, photSCs, calotowers);

      //   ////////////////////////////////////////////////////////
      //   // Get and Dump (M)E(T) Variables:
      //   // Tracks:
      //   int nqtrk(-1);
      //   fTntrackstot = tracks->size();
      //   fTTrkPtSumx = 0.; fTTrkPtSumy = 0.;
      //   for( TrackCollection::const_iterator it = tracks->begin(); it != tracks->end() ; ++it ){
      //     fTTrkPtSumx += it->px(); // Calculated for ALL tracks
      //     fTTrkPtSumy += it->py();
      //     if(it->pt() < fMintrkpt) continue;
      //     if(fabs(it->eta()) > fMaxtrketa) continue;
      //     if(it->normalizedChi2() > fMaxtrknchi2) continue;
      //     if(it->numberOfValidHits() < fMintrknhits) continue;
      //     nqtrk++; // starts at 0
      //     // Check if maximum number of tracks is exceeded already
      //     if(nqtrk >= gMaxNTrks) {
      //       edm::LogWarning("NTP") << "@SUB=analyze"
      //                              << "Maximum number of tracks exceeded";
      //       fTflagmaxtrkexc = 1;
      //       *fTGoodEvent = 1;
      //       break;
      //     }
      //     fTtrkpt[nqtrk]    = it->pt()*it->charge();
      //     fTtrketa[nqtrk]   = it->eta();
      //     fTtrkphi[nqtrk]   = it->phi();
      //     fTtrknchi2[nqtrk] = it->normalizedChi2();
      //     fTtrknhits[nqtrk] = it->numberOfValidHits();
      //     fTtrkVtxDz[nqtrk] = it->dz(primVtx->position());
      //     fTtrkVtxDxy[nqtrk]= it->dxy(primVtx->position());	
      //     fTgoodtrk[nqtrk]  = 0;
      //   }
      //   fTntracks = nqtrk+1;

      //   fTTrkPtSum = sqrt(fTTrkPtSumx*fTTrkPtSumx + fTTrkPtSumy*fTTrkPtSumy);
      //   TVector3 trkPtSum(fTTrkPtSumx, fTTrkPtSumy, 0.);
      //   fTTrkPtSumphi = trkPtSum.Phi();

      //   // Calotowers:
      //   fTNCaloTowers = calotowers->size();
      //   fTECALEsumx = 0.; fTECALEsumy = 0.; fTECALEsumz = 0.;
      //   fTHCALEsumx = 0.; fTHCALEsumy = 0.; fTHCALEsumz = 0.;
      //   fTSumEt = 0.; fTECALSumEt = 0.; fTHCALSumEt = 0.;
      //   for(CaloTowerCollection::const_iterator itow = calotowers->begin();
      //       itow!=calotowers->end(); ++itow ){
      //     if(itow->energy() == 0.) continue; // Check against zero energy towers
      //     fTSumEt += itow->et();
      //     fTECALSumEt += itow->emEt();
      //     fTHCALSumEt += itow->hadEt();
      //     double emFrac = itow->emEnergy()/itow->energy();
      //     double hadFrac = itow->hadEnergy()/itow->energy();
      //     fTECALEsumx += itow->px()*emFrac;
      //     fTECALEsumy += itow->py()*emFrac;
      //     fTECALEsumz += itow->pz()*emFrac;
      //     fTHCALEsumx += itow->px()*hadFrac;
      //     fTHCALEsumy += itow->py()*hadFrac;
      //     fTHCALEsumz += itow->pz()*hadFrac;
      //   }
      //   TVector3 ecalMET(fTECALEsumx, fTECALEsumy, fTECALEsumz);
      //   TVector3 hcalMET(fTHCALEsumx, fTHCALEsumy, fTHCALEsumz);
      //   fTECALMET    = ecalMET.Mag();
      //   fTECALMETphi = ecalMET.Phi();
      //   fTHCALMET    = hcalMET.Mag();
      //   fTHCALMETphi = hcalMET.Phi();
      //   if(fTECALEsumz != 0.) fTECALMETeta = ecalMET.Eta();
      //   else fTECALMETeta = 0.;
      //   if(fTHCALEsumz != 0.) fTHCALMETeta = hcalMET.Eta();
      //   else fTHCALMETeta = 0.;

      //   // MET Collections:
      //   fTRawMET             = (calomet->at(0)).pt();
      //   fTRawMETpx           = (calomet->at(0)).px();
      //   fTRawMETpy           = (calomet->at(0)).py();
      //   fTRawMETphi          = (calomet->at(0)).phi();
      //   fTRawMETemEtFrac     = (calomet->at(0)).emEtFraction();
      //   fTRawMETemEtInEB     = (calomet->at(0)).emEtInEB();
      //   fTRawMETemEtInEE     = (calomet->at(0)).emEtInEE();
      //   fTRawMETemEtInHF     = (calomet->at(0)).emEtInHF();
      //   fTRawMEThadEtFrac    = (calomet->at(0)).etFractionHadronic();
      //   fTRawMEThadEtInHB    = (calomet->at(0)).hadEtInHB();
      //   fTRawMEThadEtInHE    = (calomet->at(0)).hadEtInHE();
      //   fTRawMEThadEtInHF    = (calomet->at(0)).hadEtInHF();
      //   fTRawMETSignificance = (calomet->at(0)).significance();

      //   if (!fIsRealData) {
      //     Handle<View<GenMET> > GenMET;
      //     iEvent.getByLabel(fGenMETTag, GenMET);
      //     fTGenMET    = (GenMET->front()).pt();
      //     fTGenMETpx  = (GenMET->front()).px();
      //     fTGenMETpy  = (GenMET->front()).py();
      //     fTGenMETphi = (GenMET->front()).phi();
      //   }

      //   fTTCMET    = (tcmet->at(0)).pt();
      //   fTTCMETpx  = (tcmet->at(0)).px();
      //   fTTCMETpy  = (tcmet->at(0)).py();
      //   fTTCMETphi = (tcmet->at(0)).phi();
      //   fTTCMETSignificance = (tcmet->at(0)).significance();

      //   fTPFMET    = (pfmet->front()).pt();
      //   fTPFMETpx  = (pfmet->front()).px();
      //   fTPFMETpy  = (pfmet->front()).py();
      //   fTPFMETphi = (pfmet->front()).phi();
      //   fTPFMETSignificance = (pfmet->front()).significance();
      //   fTPFSumEt  = (pfmet->front()).sumEt();

      //   fTMuJESCorrMET    = (corrmujesmet->at(0)).pt();
      //   fTMuJESCorrMETpx  = (corrmujesmet->at(0)).px();
      //   fTMuJESCorrMETpy  = (corrmujesmet->at(0)).py();
      //   fTMuJESCorrMETphi = (corrmujesmet->at(0)).phi();

      //   if(fTnjets > 1){
      //     double dPhiMJ1 = TMath::Abs(reco::deltaPhi(fTjphi[0], fTMuJESCorrMETphi));
      //     double dPhiMJ2 = TMath::Abs(reco::deltaPhi(fTjphi[1], fTMuJESCorrMETphi));
      //     fTMETR12 = TMath::Sqrt(dPhiMJ1*dPhiMJ1 + (TMath::Pi()-dPhiMJ2)*(TMath::Pi()-dPhiMJ2) );
      //     fTMETR21 = TMath::Sqrt(dPhiMJ2*dPhiMJ2 + (TMath::Pi()-dPhiMJ1)*(TMath::Pi()-dPhiMJ1) );
      //   }

      //   fTPFMETPAT     = (pfMETpat->front()).pt();
      //   fTPFMETPATpx   = (pfMETpat->front()).px();
      //   fTPFMETPATpy   = (pfMETpat->front()).py();
      //   fTPFMETPATphi  = (pfMETpat->front()).phi();
      //   fTPFMETPATSignificance = (pfMETpat->at(0)).significance();



      //   ////////////////////////////////////////////////////////////////////////////////
      //   // Special stuff for Model Scans ///////////////////////////////////////////////

      // 	fTxSMS=-1;
      // 	fTxbarSMS=-1;

      //   if(!fIsRealData&&fIsModelScan) {
      //     Handle<LHEEventProduct> product;
      //     bool LHEEventProduct_found=iEvent.getByLabel("source", product);
      //     LHEEventProduct::comments_const_iterator c_begin = product->comments_begin();
      //     LHEEventProduct::comments_const_iterator c_end = product->comments_end();

      //     float mGL;
      //     float mLSP;
      //     float xCHI;

      //     for(LHEEventProduct::comments_const_iterator cit=c_begin; cit!=c_end; ++cit) {
      //       size_t found = (*cit).find("model");

      //       //# model T5zz_0.5_925.0_400.0to1000.0_450.0.lhe
      //       if( found != std::string::npos)   {
      //         size_t foundLength = (*cit).size();
      //         found = (*cit).find("T5zz");
      //         std::string smaller = (*cit).substr(found+1,foundLength);
      //         found = smaller.find("_");
      //         smaller = smaller.substr(found+1,smaller.size());
      //         std::istringstream iss(smaller);
      //         iss >> xCHI;
      //         iss.clear();
      //         found = smaller.find("_");
      //         smaller = smaller.substr(found+1,smaller.size());
      //         iss.str(smaller);
      //         iss >> mGL;
      //         iss.clear();
      //         found = smaller.find("_");
      //         smaller = smaller.substr(found+1,smaller.size());
      //         iss.str(smaller);
      //         iss >> mLSP;
      //         iss.clear();

      //         //model msugra_1900_100_10_0_1
      //         float m0,m12,tanb,A0,signMu;
      //         foundLength = (*cit).size();
      //         found = (*cit).find("=");
      //         smaller = (*cit).substr(found+1,foundLength);
      //         found = smaller.find("_");
      //         smaller = smaller.substr(found+1,smaller.size());
      //         //
      //         iss.clear();
      //         iss.str(smaller);
      //         iss >> m0;
      //         iss.clear();
      //         //
      //         found = smaller.find("_");
      //         smaller = smaller.substr(found+1,smaller.size());
      //         iss.str(smaller);
      //         iss >> m12;
      //         iss.clear();
      //         //
      //         found = smaller.find("_");
      //         smaller = smaller.substr(found+1,smaller.size());
      //         iss.str(smaller);
      //         iss >> tanb;
      //         iss.clear();
      //         //
      //         found = smaller.find("_");
      //         smaller = smaller.substr(found+1,smaller.size());
      //         iss.str(smaller);
      //         iss >> A0;
      //         iss.clear();
      //         //
      //         found = smaller.find("_");
      //         smaller = smaller.substr(found+1,smaller.size());
      //         iss.str(smaller);
      //         iss >> signMu;
      //         iss.clear();

      //         // mSUGRA scan
      //         fTSUSYScanM0=m0;
      //         fTSUSYScanM12=m12;
      //         fTSUSYScanMu=signMu;
      //         fTSUSYScanA0=A0;

      //       }
      //     }
      //     fTMassGlu = mGL;
      //     fTMassChi = xCHI;
      //     fTMassLSP = mLSP;
      //   }


      // 	bool blabalot=false;
      // 	bool BloatWithGenInfo=true;

	
      // 	if(BloatWithGenInfo && !fIsRealData) {
      // 	  int genIndex[1200];
      // 	  int genNIndex[1200];
      // 	  int genID[1200];
      // 	  float genPx[1200];
      // 	  float genPy[1200];
      // 	  float genPz[1200];
      // 	  float genPt[1200];
      // 	  float genEta[1200];
      // 	  float genPhi[1200];
      // 	  float genM[1200];
      // 	  int genMo1Index[1200];
      // 	  int genMo2Index[1200];
      // 	  int genNMo[1200];
      // 	  int genStatus[1200];
      // 	  int Promptness[1200];
      // 	  bool StoreFlag[1200];
	  
      // 	  int nGenParticles=0;
	  
      // 	  Handle<GenParticleCollection> genParticles;
      // 	  iEvent.getByLabel("genParticles", genParticles);
	  
      // 	  // STEP 1: Loop over all particles and store the information.
      // 	  for(size_t i = 0; i < genParticles->size(); ++ i) {
      // 	    nGenParticles++;
      // 	    const GenParticle & p = (*genParticles)[i];
      // 	    genIndex[i]=i;
      // 	    genNIndex[i]=0;
      // 	    genID[i]=p.pdgId();
      // 	    genPx[i]=p.px();
      // 	    genPy[i]=p.py();
      // 	    genPz[i]=p.pz();
      // 	    genPt[i]=p.pt();
      // 	    genPhi[i]=p.phi();
      // 	    genEta[i]=p.eta();
      // 	    genM[i]=p.mass();
      // 	    genStatus[i]=p.status();
      // 	    genMo1Index[i]=-1;
      // 	    genMo2Index[i]=-1;
      // 	    genNMo[i]=p.numberOfMothers();
      // 	    StoreFlag[i]=false;
      // 	    if(genID[i]==2212) Promptness[i]=0;
      // 	    else Promptness[i]=-1;
	    

      // 	    if(blabalot) cout << "Reading particle " << i << " (is pdgid=" << genID[i] << ") with " << genNMo[i] << " mothers" << endl;
      // 	    for(int j=0;j<genNMo[i]&&j<2;j++) {
      // 	      int idx = -1;
      // 	      const GenParticle* mom = static_cast<const GenParticle*> (p.mother(j));
      // 	      const GenParticle* gmom = static_cast<const GenParticle*> (mom->mother());
      // 	      if (mom == NULL) break;
      // 	      for (unsigned int k = 0; k < genParticles->size() && k<i; ++k) {
      // 		const reco::GenParticle& testm = (*genParticles)[k];
      // 		const GenParticle* testgm = static_cast<const GenParticle*> (testm.mother());
      // 		if (testm.pt() == mom->pt()) {
      // 		  if(gmom == NULL || (testgm!=NULL) && (gmom->pt() == testgm->pt())) {
      // 		    idx = k;
      // 		    if(blabalot) cout << "     Found a hit for a mother for index " << i << "! It's index " << idx << " (pdgid " << genID[idx] << ")" << endl;
      // 		    break;
      // 		  }
      // 		}
      // 	      }
      // 	      if(j==0) genMo1Index[i]=idx;
      // 	      if(j==1) genMo2Index[i]=idx;
      // 	    }
      // 	    if(Promptness[i]==-1&&genMo1Index[i]>=0&&Promptness[genMo1Index[i]]>-1) Promptness[i]=Promptness[genMo1Index[i]]+1;
      // 	    if(blabalot) cout << "                    Promtpness: " << Promptness[i] << endl;
      // 	    if(blabalot) cout << "          Mother 1: " << genMo1Index[i] << " with Promptness " << Promptness[genMo1Index[i]] << endl;
      // 	    if(blabalot) cout << "          Mother 2: " << genMo2Index[i] << endl;
      // 	  }
	  
      // 	  // STEP 2: Loop from end to start flipping storeflag when necessary
      // 	  float genPtThreshold=5.0; // any particle with less pt than this will not be stored.
      // 	  for(int i=nGenParticles-1;i>=0;i--) {
      // 	    if(genStatus[i]!=1) continue;
      // 	    if(genPt[i]<genPtThreshold) continue;
      // 	    FlipGenStoreFlag(i,Promptness,genMo1Index,genMo2Index,StoreFlag);
      // 	  }
	  
      // 	  // Intermediate step: Make sure that all first particles are stored, and that the earliest particles are stored (i.e. promptness criteria are met)
      // 	  for(int i=nGenParticles-1;i>=0;i--) {
      // 	    if(Promptness[i]<4&&Promptness[i]>=0) StoreFlag[i]=true;
      // 	    if(i<20) StoreFlag[i]=true;
      // 	  }
	  
      // 	  // STEP 3: Loop again, setting new index ("nindex") and replacing Mo1 by if(Mo1>-1) Mo1=nindex[index];, same for Mo2. If storeflag, store it directly.
      // 	  fTnGenParticles=-1;
      // 	  for(int i=0;i<nGenParticles;i++) {
      // 	    if(!StoreFlag[i]) continue;
      // 	    fTnGenParticles++;
      // 	    if(fTnGenParticles>nStoredGenParticles) {
      // 	      edm::LogWarning("NTP") << "@SUB=analyze()"
      // 	      << "Maximum number of gen particles exceeded";
      // 	      fTflagmaxgenpartexc = 1;
      // 	      break;
      // 	    }
	    
      // 	    genNIndex[i]=fTnGenParticles;
	    
      // 	    //store everything
      // 	    fTgenInfoId[fTnGenParticles] = genID[i];
      // 	    fTgenInfoStatus[fTnGenParticles] = genStatus[i];
      // 	    fTgenInfoNMo[fTnGenParticles] = genNMo[i];
      // 	    fTgenInfoMo1Pt[fTnGenParticles]=0;
      // 	    fTgenInfoMo2Pt[fTnGenParticles]=0;
      // 	    if(genMo1Index[i]>=0) {
      // 	      fTgenInfoMo1Pt[fTnGenParticles]=genPt[genMo1Index[i]];
      // 	      fTgenInfoMo1[fTnGenParticles]=genNIndex[genMo1Index[i]];
      // 	    } else {
      // 	      fTgenInfoMo1[fTnGenParticles]=-1;
      // 	    }
      // 	    if(genMo2Index[i]>=0) {
      // 	      fTgenInfoMo2Pt[fTnGenParticles]=genPt[genMo2Index[i]];
      // 	      fTgenInfoMo2[fTnGenParticles]=genNIndex[genMo2Index[i]];
      // 	    } else {
      // 	      fTgenInfoMo2[fTnGenParticles]=-1;
      // 	    }
	    
      // 	    fTPromptnessLevel[fTnGenParticles]=Promptness[i];
      // 	    fTgenInfoPt[fTnGenParticles] = genPt[i];
      // 	    fTgenInfoEta[fTnGenParticles] = genEta[i];
      // 	    fTgenInfoPhi[fTnGenParticles] = genPhi[i];
      // 	    fTgenInfoPx[fTnGenParticles] = genPx[i];
      // 	    fTgenInfoPy[fTnGenParticles] = genPy[i];
      // 	    fTgenInfoPz[fTnGenParticles] = genPz[i];
      // 	    fTgenInfoM[fTnGenParticles] = genM[i];

      // 	    if(blabalot) {
      // 	      cout << "Working particle " << i << " with Promptness: " << Promptness[i] << "  (storage number: " << fTnGenParticles << ")" << endl;
      // 	      cout << "    Mother is Particle " << genMo1Index[i] << " with Promptness " << Promptness[genMo1Index[i]] << endl;
      // 	      cout << "Particle " << fTnGenParticles << "  (" << i << "): The particle has ID = " << genID[i]                     << " and its mother has index " << genNIndex[genMo1Index[i]]     << " mo pt : " << genPt[genMo1Index[i]] << endl;
      // 	      cout << "stored:  " << fTnGenParticles << "  (" << i << ")                      " << fTgenInfoId[fTnGenParticles] << "                          " << fTgenInfoMo1[fTnGenParticles] << "         " << fTgenInfoMo1Pt[fTnGenParticles] << endl;
      // 	    }
	    
      // 	  }
	  
      // 	  fTnGenParticles++;
      // 	  if(blabalot) cout << "A total of " << fTnGenParticles << " Particles  have been stored out of " << nGenParticles << " ( " << 100*fTnGenParticles/(float)nGenParticles << " %)" << endl;
	  
	  
      // 	  float gluinomass=0;
      // 	  float ngluinomass=0;
      // 	  float chi2mass=0;
      // 	  float nchi2mass=0;
      // 	  float chi1mass=0;
      // 	  float nchi1mass=0;
	  
      // 	  for(size_t i = 0; i < genParticles->size(); ++ i) {
      // 	    const GenParticle & p = (*genParticles)[i];
      // 	    int id = p.pdgId();
      // 	    double mass = p.mass();
      // 	    if(id==1000021) {
      // 	      gluinomass+=mass;
      // 	      ngluinomass++;
      // 	    }
      // 	    if(id==1000023) {
      // 	      chi2mass+=mass;
      // 	      nchi2mass++;
      // 	    }
      // 	    if(id==1000022) {
      // 	      chi1mass+=mass;
      // 	      nchi1mass++;
      // 	    }
      // 	  }
	  
      // 	  if(ngluinomass>0&&nchi2mass>0&&nchi1mass>0) {
      // 	    fTxSMS = (chi2mass/nchi1mass - chi1mass/nchi1mass) / (gluinomass/ngluinomass - chi1mass/nchi1mass);
      // 	    fTxbarSMS = 1 - fTxSMS; // Mariarosaria's definition of x
      // 	  }
      // 	}// end of bloat with gen information


      ////////////////////////////////////////////////////////
      // Process other jet collections, as configured
      for ( std::vector<JetFillerBase*>::iterator it = jetFillers.begin();
            it != jetFillers.end(); ++it )
        (*it)->fillProducts(iEvent,iSetup);
      for ( std::vector<PatMuonFiller*>::iterator it = muonFillers.begin(); 
            it != muonFillers.end(); ++it ) 
        (*it)->fillProducts(iEvent,iSetup);
      for ( std::vector<PatElectronFiller*>::iterator it = electronFillers.begin(); 
            it != electronFillers.end(); ++it ) 
        (*it)->fillProducts(iEvent,iSetup);
      for ( std::vector<PatTauFiller*>::iterator it = tauFillers.begin(); 
            it != tauFillers.end(); ++it ) 
        (*it)->fillProducts(iEvent,iSetup);

      ///////////////////////////////////////////////////////////////////////////////
      // Fill Tree //////////////////////////////////////////////////////////////////
      putProducts( iEvent );
      for ( std::vector<JetFillerBase*>::iterator it = jetFillers.begin();
            it != jetFillers.end(); ++it )
        (*it)->putProducts(iEvent);
      for ( std::vector<PatMuonFiller*>::iterator it = muonFillers.begin(); 
            it != muonFillers.end(); ++it ) 
        (*it)->putProducts(iEvent);
      for ( std::vector<PatElectronFiller*>::iterator it = electronFillers.begin(); 
            it != electronFillers.end(); ++it ) 
        (*it)->putProducts(iEvent);
      for ( std::vector<PatTauFiller*>::iterator it = tauFillers.begin(); 
            it != tauFillers.end(); ++it ) 
        (*it)->putProducts(iEvent);

      fNFillTree++;

      // Not used as a filter right now
      return true;
    }

    //________________________________________________________________________________________
    // Declare all products
    void NTupleProducer::declareProducts(void) { 

      // Run products
      produces<float,edm::InRun>("ExtXSecLO"     );
      produces<float,edm::InRun>("ExtXSecNLO"    );
      produces<float,edm::InRun>("IntXSec"       );
      produces<float,edm::InRun>("MinMuPt"       );
      produces<float,edm::InRun>("MaxMuEta"      );
      produces<float,edm::InRun>("MinElPt"       );
      produces<float,edm::InRun>("MaxElEta"      );
      produces<float,edm::InRun>("MinJPt"        );
      produces<float,edm::InRun>("MinRawJPt"     );
      produces<float,edm::InRun>("MaxJEta"       );
      produces<float,edm::InRun>("MinJEMfrac"    );
                                                                
      produces<float,edm::InRun>("MinTrkPt"      );
      produces<float,edm::InRun>("MaxTrkEta"     );
      produces<float,edm::InRun>("MaxTrkNChi2"   );
      produces<int,edm::InRun>  ("MinTrkNHits"   );
                                                                
      produces<float,edm::InRun>("MinPhotonPt"   );
      produces<float,edm::InRun>("MaxPhotonEta"  );
      produces<float,edm::InRun>("MinSCraw"      );
      produces<float,edm::InRun>("MinEBRechitE"  );
  
      produces<float,edm::InRun>("MinGenLeptPt"  );
      produces<float,edm::InRun>("MaxGenLeptEta" );
      produces<float,edm::InRun>("MinGenPhotPt"  );
      produces<float,edm::InRun>("MaxGenPhotEta" );
      produces<float,edm::InRun>("MinGenJetPt"   );
      produces<float,edm::InRun>("MaxGenJetEta"  );

      produces<float,edm::InRun>("BtagMatchDeltaR"  );
                                                                
      produces<int,edm::InRun>("MaxNMus"       );
      produces<int,edm::InRun>("MaxNEles"      );
      produces<int,edm::InRun>("MaxNJets"      );
      produces<int,edm::InRun>("MaxNTrks"      );
      produces<int,edm::InRun>("MaxNPhotons"   );
      produces<int,edm::InRun>("MaxNSC"        );
      produces<int,edm::InRun>("MaxNGenLep"    );
      produces<int,edm::InRun>("MaxNGenPho"    );
      produces<int,edm::InRun>("MaxNGenJet"    );
      produces<int,edm::InRun>("MaxNVrtx"      );
      produces<int,edm::InRun>("MaxNPileup"    );
      produces<int,edm::InRun>("MaxNEBhits"    );

      produces<std::vector<std::string>,edm::InRun>("HLTNames");
      produces<std::vector<std::string>,edm::InRun>("L1PhysMenu");
      produces<std::vector<std::string>,edm::InRun>("HLTLabels");

      produces<std::vector<std::string>,edm::InRun>("PileUpData");
      produces<std::vector<std::string>,edm::InRun>("PileUpMC");

      // Event products
      produces<int>("Run");
      produces<int>("Event");
      produces<int>("LumiSection");
      produces<float>("PtHat");
      produces<float>("QCDPartonicHT");
      produces<int>("SigProcID");
      produces<float>("PDFScalePDF");
      produces<int>("PDFID1");
      produces<int>("PDFID2");
      produces<float>("PDFx1");
      produces<float>("PDFx2");
      produces<float>("PDFxPDF1");
      produces<float>("PDFxPDF2");
      produces<float>("GenWeight");
      produces<std::vector<float> >("pdfW");
      produces<float>("pdfWsum");
      produces<int>("NPdfs");
      produces<int>("PUnumInteractions");
      produces<int>("PUnumTrueInteractions");
      produces<int>("PUnumFilled");
      produces<int>("PUOOTnumInteractionsEarly");
      produces<int>("PUOOTnumInteractionsLate");
      produces<std::vector<float> >("PUzPositions");
      produces<std::vector<float> >("PUsumPtLowPt");
      produces<std::vector<float> >("PUsumPtHighPt");
      produces<std::vector<float> >("PUnTrksLowPt");
      produces<std::vector<float> >("PUnTrksHighPt");
      produces<float>("Rho");
      produces<float>("RhoPFnoPU");
      produces<float>("Weight");
      produces<std::vector<int> >("HLTResults");
      produces<std::vector<int> >("HLTPrescale");
      produces<std::vector<int> >("L1PhysResults");
      produces<std::vector<int> >("L1TechResults");
      produces<int>("NHLTObjs");
      for ( size_t i=0; i<gMaxHltNObjs; ++i ) {
        std::ostringstream s;
        s << i;
        produces<std::vector<int> >(("HLTObjectID"+s.str()).c_str());
        produces<std::vector<float> >(("HLTObjectPt"+s.str()).c_str());
        produces<std::vector<float> >(("HLTObjectEta"+s.str()).c_str());
        produces<std::vector<float> >(("HLTObjectPhi"+s.str()).c_str());
      }
      produces<float>("PUWeightTotal");
      produces<float>("PUWeightInTime");
      produces<float>("MassGlu");
      produces<float>("MassChi");
      produces<float>("MassLSP");
      produces<float>("xSMS");
      produces<float>("xbarSMS");
      produces<float>("M0");
      produces<float>("M12");
      produces<float>("signMu");
      produces<float>("A0");
      produces<int>("process");

      produces<int>("MaxGenPartExceed");
      produces<int>("nGenParticles");
      produces<std::vector<int> >("genInfoId");
      produces<std::vector<int> >("genInfoStatus");
      produces<std::vector<int> >("genInfoNMo");
      produces<std::vector<int> >("genInfoNDa");
      produces<std::vector<int> >("genInfoMo1");
      produces<std::vector<int> >("genInfoMo2");
      produces<std::vector<int> >("genInfoDa1");
      produces<std::vector<int> >("genInfoDa2");
      produces<std::vector<int> >("genInfoMoIndex");
      produces<std::vector<int> >("PromptnessLevel");
      produces<std::vector<float> >("genInfoMass");
      produces<std::vector<float> >("genInfoMo1Pt");
      produces<std::vector<float> >("genInfoMo2Pt");
      produces<std::vector<float> >("genInfoPt");
      produces<std::vector<float> >("genInfoEta");
      produces<std::vector<float> >("genInfoPhi");
      produces<std::vector<float> >("genInfoPx");
      produces<std::vector<float> >("genInfoPy");
      produces<std::vector<float> >("genInfoPz");
      produces<std::vector<float> >("genInfoM");
      produces<std::vector<float> >("genInfoPromptFlag");

      produces<int>("PrimVtxGood");
      produces<float>("PrimVtxx");
      produces<float>("PrimVtxy");
      produces<float>("PrimVtxz");
      produces<float>("PrimVtxRho");
      produces<float>("PrimVtxxE");
      produces<float>("PrimVtxyE");
      produces<float>("PrimVtxzE");
      produces<float>("PrimVtxNChi2");
      produces<float>("PrimVtxNdof");
      produces<int>("PrimVtxIsFake");
      produces<float>("PrimVtxPtSum");
      produces<float>("Beamspotx");
      produces<float>("Beamspoty");
      produces<float>("Beamspotz");
      produces<int>("NCaloTowers");
      produces<int>("GoodEvent");
      produces<int>("MaxMuExceed");
      produces<int>("MaxElExceed");
      produces<int>("MaxJetExceed");
      produces<int>("MaxUncJetExceed");
      produces<int>("MaxTrkExceed");
      produces<int>("MaxPhotonsExceed");
      produces<int>("MaxGenLepExceed");
      produces<int>("MaxGenPhoExceed");
      produces<int>("MaxGenJetExceed");
      produces<int>("MaxVerticesExceed");
      produces<int>("HBHENoiseFlag");
      produces<int>("HBHENoiseFlagIso");
      produces<int>("CSCTightHaloID");
      produces<int>("EcalDeadTPFilterFlag");
      produces<int>("RecovRecHitFilterFlag");
      produces<int>("RA2TrackingFailureFilterFlag");
      //FR produces<int>("PBNRFlag");
      produces<int>("NGenLeptons");
      produces<std::vector<int> >("GenLeptonID");
      produces<std::vector<float> >("GenLeptonPt");
      produces<std::vector<float> >("GenLeptonEta");
      produces<std::vector<float> >("GenLeptonPhi");
      produces<std::vector<int> >("GenLeptonMID");
      produces<std::vector<int> >("GenLeptonMStatus");
      produces<std::vector<float> >("GenLeptonMPt");
      produces<std::vector<float> >("GenLeptonMEta");
      produces<std::vector<float> >("GenLeptonMPhi");
      produces<std::vector<int> >("GenLeptonGMID");
      produces<std::vector<int> >("GenLeptonGMStatus");
      produces<std::vector<float> >("GenLeptonGMPt");
      produces<std::vector<float> >("GenLeptonGMEta");
      produces<std::vector<float> >("GenLeptonGMPhi");
      produces<int>("NGenPhotons");
      produces<std::vector<float> >("GenPhotonPt");
      produces<std::vector<float> >("GenPhotonEta");
      produces<std::vector<float> >("GenPhotonPhi");
      produces<std::vector<float> >("GenPhotonPartonMindR");
      produces<std::vector<int> >("GenPhotonMotherID");
      produces<std::vector<int> >("GenPhotonMotherStatus");
      produces<int>("NGenJets");
      produces<std::vector<float> >("GenJetPt");
      produces<std::vector<float> >("GenJetEta");
      produces<std::vector<float> >("GenJetPhi");
      produces<std::vector<float> >("GenJetE");
      produces<std::vector<float> >("GenJetEmE");
      produces<std::vector<float> >("GenJetHadE");
      produces<std::vector<float> >("GenJetInvE");
      produces<int>("NVrtx");
      produces<std::vector<float> >("VrtxX");
      produces<std::vector<float> >("VrtxY");
      produces<std::vector<float> >("VrtxZ");
      produces<std::vector<float> >("VrtxXE");
      produces<std::vector<float> >("VrtxYE");
      produces<std::vector<float> >("VrtxZE");
      produces<std::vector<float> >("VrtxNdof");
      produces<std::vector<float> >("VrtxChi2");
      produces<std::vector<float> >("VrtxNtrks");
      produces<std::vector<float> >("VrtxSumPt");
      produces<std::vector<int> >("VrtxIsFake");
      produces<int>("NMus");
      produces<int>("NMusTot");
      produces<int>("NGMus");
      produces<int>("NTMus");
      produces<std::vector<int> >("MuGood");
      produces<std::vector<int> >("MuIsIso");
      produces<std::vector<int> >("MuIsGlobalMuon");
      produces<std::vector<int> >("MuIsTrackerMuon");
      produces<std::vector<float> >("MuPx");
      produces<std::vector<float> >("MuPy");
      produces<std::vector<float> >("MuPz");
      produces<std::vector<float> >("MuPt");
      produces<std::vector<float> >("MuInnerTkPt");
      produces<std::vector<float> >("MuPtE");
      produces<std::vector<float> >("MuE");
      produces<std::vector<float> >("MuEt");
      produces<std::vector<float> >("MuEta");
      produces<std::vector<float> >("MuPhi");
      produces<std::vector<int> >("MuCharge");
      produces<std::vector<float> >("MuRelIso03");
      produces<std::vector<float> >("MuIso03SumPt");
      produces<std::vector<float> >("MuIso03EmEt");
      produces<std::vector<float> >("MuIso03HadEt");
      produces<std::vector<float> >("MuIso03EMVetoEt");
      produces<std::vector<float> >("MuIso03HadVetoEt");
      produces<std::vector<float> >("MuIso05SumPt");
      produces<std::vector<float> >("MuIso05EmEt");
      produces<std::vector<float> >("MuIso05HadEt");
      produces<std::vector<float> >("MuEem");
      produces<std::vector<float> >("MuEhad");
      produces<std::vector<float> >("MuD0BS");
      produces<std::vector<float> >("MuD0PV");
      produces<std::vector<float> >("MuD0E");
      produces<std::vector<float> >("MuDzBS");
      produces<std::vector<float> >("MuDzPV");
      produces<std::vector<float> >("MuDzE");
      produces<std::vector<float> >("MuNChi2");
      produces<std::vector<int> >("MuNGlHits");
      produces<std::vector<int> >("MuNMuHits");
      produces<std::vector<int> >("MuNTkHits");
      produces<std::vector<int> >("MuNPxHits");
      produces<std::vector<float> >("MuInnerTkNChi2");
      produces<std::vector<int> >("MuNMatches");
      produces<std::vector<int> >("MuNChambers");
      produces<std::vector<float> >("MuCaloComp");
      produces<std::vector<float> >("MuSegmComp");
      produces<std::vector<int> >("MuIsGMPT");
      produces<std::vector<int> >("MuIsGMTkChiComp");
      produces<std::vector<int> >("MuIsGMStaChiComp");
      produces<std::vector<int> >("MuIsGMTkKinkTight");
      produces<std::vector<int> >("MuIsAllStaMuons");
      produces<std::vector<int> >("MuIsAllTrkMuons");
      produces<std::vector<int> >("MuIsTrkMuonArbitrated");
      produces<std::vector<int> >("MuIsAllArbitrated");
      produces<std::vector<int> >("MuIsTMLSLoose");
      produces<std::vector<int> >("MuIsTMLSTight");
      produces<std::vector<int> >("MuIsTM2DCompLoose");
      produces<std::vector<int> >("MuIsTM2DCompTight");
      produces<std::vector<int> >("MuIsTMOneStationLoose");
      produces<std::vector<int> >("MuIsTMOneStationTight");
      produces<std::vector<int> >("MuIsTMLSOptLowPtLoose");
      produces<std::vector<int> >("MuIsTMLSAngLoose");
      produces<std::vector<int> >("MuIsTMLSAngTight");
      produces<std::vector<int> >("MuIsTMOneStationAngTight");
      produces<std::vector<int> >("MuIsTMOneStationAngLoose");
      produces<std::vector<int> >("MuGenID");
      produces<std::vector<int> >("MuGenStatus");
      produces<std::vector<float> >("MuGenPt");
      produces<std::vector<float> >("MuGenEta");
      produces<std::vector<float> >("MuGenPhi");
      produces<std::vector<float> >("MuGenE");
      produces<std::vector<int> >("MuGenMID");
      produces<std::vector<int> >("MuGenMStatus");
      produces<std::vector<float> >("MuGenMPt");
      produces<std::vector<float> >("MuGenMEta");
      produces<std::vector<float> >("MuGenMPhi");
      produces<std::vector<float> >("MuGenME");
      produces<std::vector<int> >("MuGenGMID");
      produces<std::vector<int> >("MuGenGMStatus");
      produces<std::vector<float> >("MuGenGMPt");
      produces<std::vector<float> >("MuGenGMEta");
      produces<std::vector<float> >("MuGenGMPhi");
      produces<std::vector<float> >("MuGenGME");
      produces<int>("NEBhits");
      produces<std::vector<float> >("EBrechitE");
      produces<std::vector<float> >("EBrechitPt");
      produces<std::vector<float> >("EBrechitEta");
      produces<std::vector<float> >("EBrechitPhi");
      produces<std::vector<float> >("EBrechitChi2");
      produces<std::vector<float> >("EBrechitTime");
      produces<std::vector<float> >("EBrechitE4oE1");
      produces<std::vector<float> >("EBrechitE2oE9");
      produces<int>("NEles");
      produces<int>("NElesTot");
      produces<std::vector<int> >("ElGood");
      produces<std::vector<int> >("ElIsIso");
      produces<std::vector<int> >("ElChargeMisIDProb");
      produces<std::vector<float> >("ElPx");
      produces<std::vector<float> >("ElPy");
      produces<std::vector<float> >("ElPz");
      produces<std::vector<float> >("ElPt");
      produces<std::vector<float> >("ElPtE");
      produces<std::vector<float> >("ElE");
      produces<std::vector<float> >("ElEt");
      produces<std::vector<float> >("ElEta");
      produces<std::vector<float> >("ElTheta");
      produces<std::vector<float> >("ElSCEta");
      produces<std::vector<float> >("ElPhi");
      produces<std::vector<float> >("ElGsfTkPt");
      produces<std::vector<float> >("ElGsfTkEta");
      produces<std::vector<float> >("ElGsfTkPhi");
      produces<std::vector<float> >("ElTrkMomentumError");
      produces<std::vector<float> >("ElEcalEnergyError");
      produces<std::vector<float> >("ElEleMomentumError");
      produces<std::vector<int> >("ElNBrems");
      produces<std::vector<float> >("ElD0BS");
      produces<std::vector<float> >("ElD0PV");
      produces<std::vector<float> >("ElD0E");
      produces<std::vector<float> >("ElDzBS");
      produces<std::vector<float> >("ElDzPV");
      produces<std::vector<float> >("ElDzE");
      produces<std::vector<float> >("ElRelIso03");
      produces<std::vector<float> >("ElRelIso04");
      produces<std::vector<float> >("ElDR03TkSumPt");
      produces<std::vector<float> >("ElDR04TkSumPt");
      produces<std::vector<float> >("ElDR03EcalRecHitSumEt");
      produces<std::vector<float> >("ElDR04EcalRecHitSumEt");
      produces<std::vector<float> >("ElDR03HcalTowerSumEt");
      produces<std::vector<float> >("ElDR04HcalTowerSumEt");
      produces<std::vector<float> >("ElNChi2");
      produces<std::vector<int> >("ElCharge");
      produces<std::vector<int> >("ElCInfoIsGsfCtfCons");
      produces<std::vector<int> >("ElCInfoIsGsfCtfScPixCons");
      produces<std::vector<int> >("ElCInfoIsGsfScPixCons");
      produces<std::vector<int> >("ElScPixCharge");
      produces<std::vector<float> >("ElClosestCtfTrackPt");
      produces<std::vector<float> >("ElClosestCtfTrackEta");
      produces<std::vector<float> >("ElClosestCtfTrackPhi");
      produces<std::vector<int> >("ElClosestCtfTrackCharge");
      produces<std::vector<float> >("ElIDMva");
      produces<std::vector<int> >("ElIDTight");
      produces<std::vector<int> >("ElIDLoose");
      produces<std::vector<int> >("ElIDRobustTight");
      produces<std::vector<int> >("ElIDRobustLoose");
      produces<std::vector<int> >("ElIDsimpleWPrelIso");
      produces<std::vector<int> >("ElIDsimpleWP80relIso");
      produces<std::vector<int> >("ElIDsimpleWP85relIso");
      produces<std::vector<int> >("ElIDsimpleWP90relIso");
      produces<std::vector<int> >("ElIDsimpleWP95relIso");
      produces<std::vector<int> >("ElInGap");
      produces<std::vector<int> >("ElEcalDriven");
      produces<std::vector<int> >("ElTrackerDriven");
      produces<std::vector<int> >("ElBasicClustersSize");
      produces<std::vector<float> >("Elfbrem");
      produces<std::vector<float> >("ElHcalOverEcal");
      produces<std::vector<float> >("ElE1x5");
      produces<std::vector<float> >("ElE5x5");
      produces<std::vector<float> >("ElE2x5Max");
      produces<std::vector<float> >("ElSigmaIetaIeta");
      produces<std::vector<float> >("ElDeltaPhiSeedClusterAtCalo");
      produces<std::vector<float> >("ElDeltaEtaSeedClusterAtCalo");
      produces<std::vector<float> >("ElDeltaPhiSuperClusterAtVtx");
      produces<std::vector<float> >("ElDeltaEtaSuperClusterAtVtx");
      produces<std::vector<float> >("ElCaloEnergy");
      produces<std::vector<float> >("ElTrkMomAtVtx");
      produces<std::vector<float> >("ElESuperClusterOverP");
      produces<std::vector<int> >("ElNumberOfMissingInnerHits");
      produces<std::vector<int> >("ElSCindex");
      produces<std::vector<float> >("ElConvPartnerTrkDist");
      produces<std::vector<float> >("ElConvPartnerTrkDCot");
      produces<std::vector<float> >("ElConvPartnerTrkPt");
      produces<std::vector<float> >("ElConvPartnerTrkEta");
      produces<std::vector<float> >("ElConvPartnerTrkPhi");
      produces<std::vector<float> >("ElConvPartnerTrkCharge");
      produces<std::vector<int> >("ElScSeedSeverity");
      produces<std::vector<float> >("ElE1OverE9");
      produces<std::vector<float> >("ElS4OverS1");
      produces<std::vector<int> >("ElGenID");
      produces<std::vector<int> >("ElGenStatus");
      produces<std::vector<float> >("ElGenPt");
      produces<std::vector<float> >("ElGenEta");
      produces<std::vector<float> >("ElGenPhi");
      produces<std::vector<float> >("ElGenE");
      produces<std::vector<int> >("ElGenMID");
      produces<std::vector<int> >("ElGenMStatus");
      produces<std::vector<float> >("ElGenMPt");
      produces<std::vector<float> >("ElGenMEta");
      produces<std::vector<float> >("ElGenMPhi");
      produces<std::vector<float> >("ElGenME");
      produces<std::vector<int> >("ElGenGMID");
      produces<std::vector<int> >("ElGenGMStatus");
      produces<std::vector<float> >("ElGenGMPt");
      produces<std::vector<float> >("ElGenGMEta");
      produces<std::vector<float> >("ElGenGMPhi");
      produces<std::vector<float> >("ElGenGME");
      produces<int>("NPhotons");
      produces<int>("NPhotonsTot");
      produces<std::vector<int> >("PhoGood");
      produces<std::vector<int> >("PhoIsIso");
      produces<std::vector<float> >("PhoPt");
      produces<std::vector<float> >("PhoPx");
      produces<std::vector<float> >("PhoPy");
      produces<std::vector<float> >("PhoPz");
      produces<std::vector<float> >("PhoEta");
      produces<std::vector<float> >("PhoPhi");
      produces<std::vector<float> >("PhoEnergy");
      produces<std::vector<float> >("PhoIso03Ecal");
      produces<std::vector<float> >("PhoIso03Hcal");
      produces<std::vector<float> >("PhoIso03TrkSolid");
      produces<std::vector<float> >("PhoIso03TrkHollow");
      produces<std::vector<float> >("PhoIso03");
      produces<std::vector<float> >("PhoIso04Ecal");
      produces<std::vector<float> >("PhoIso04Hcal");
      produces<std::vector<float> >("PhoIso04TrkSolid");
      produces<std::vector<float> >("PhoIso04TrkHollow");
      produces<std::vector<float> >("PhoIso04");
      produces<std::vector<float> >("PhoR9");
      produces<std::vector<float> >("PhoCaloPositionX");
      produces<std::vector<float> >("PhoCaloPositionY");
      produces<std::vector<float> >("PhoCaloPositionZ");
      produces<std::vector<float> >("PhoHoverE");
      produces<std::vector<float> >("PhoH1overE");
      produces<std::vector<float> >("PhoH2overE");
      produces<std::vector<float> >("PhoSigmaIetaIeta");
      produces<std::vector<float> >("PhoSCRawEnergy");
      produces<std::vector<float> >("PhoSCEtaWidth");
      produces<std::vector<float> >("PhoSCSigmaPhiPhi");
      produces<std::vector<int> >("PhoHasPixSeed");
      produces<std::vector<int> >("PhoHasConvTrks");
      produces<std::vector<int> >("PhoScSeedSeverity");
      produces<std::vector<float> >("PhoE1OverE9");
      produces<std::vector<float> >("PhoS4OverS1");
      produces<std::vector<float> >("PhoSigmaEtaEta");
      produces<std::vector<float> >("PhoE1x5");
      produces<std::vector<float> >("PhoE2x5");
      produces<std::vector<float> >("PhoE3x3");
      produces<std::vector<float> >("PhoE5x5");
      produces<std::vector<float> >("PhomaxEnergyXtal");
      produces<std::vector<float> >("PhoIso03HcalDepth1");
      produces<std::vector<float> >("PhoIso03HcalDepth2");
      produces<std::vector<float> >("PhoIso04HcalDepth1");
      produces<std::vector<float> >("PhoIso04HcalDepth2");
      produces<std::vector<int> >("PhoIso03nTrksSolid");
      produces<std::vector<int> >("PhoIso03nTrksHollow");
      produces<std::vector<int> >("PhoIso04nTrksSolid");
      produces<std::vector<int> >("PhoIso04nTrksHollow");
      produces<std::vector<int> >("PhoisEB");
      produces<std::vector<int> >("PhoisEE");
      produces<std::vector<int> >("PhoisEBEtaGap");
      produces<std::vector<int> >("PhoisEBPhiGap");
      produces<std::vector<int> >("PhoisEERingGap");
      produces<std::vector<int> >("PhoisEEDeeGap");
      produces<std::vector<int> >("PhoisEBEEGap");
      produces<std::vector<int> >("PhoisPFlowPhoton");
      produces<std::vector<int> >("PhoisStandardPhoton");
      produces<std::vector<int> >("PhoMCmatchindex");
      produces<std::vector<int> >("PhoMCmatchexitcode");
      produces<std::vector<float> >("PhoChargedHadronIso");
      produces<std::vector<float> >("PhoNeutralHadronIso");
      produces<std::vector<float> >("PhoPhotonIso");
      produces<std::vector<int> >("PhoisPFPhoton");
      produces<std::vector<int> >("PhoisPFElectron");
      produces<std::vector<int> >("PhotSCindex");
      produces<std::vector<float> >("PhoCone04PhotonIsodR0dEta0pt0");
      produces<std::vector<float> >("PhoCone04PhotonIsodR0dEta0pt5");
      produces<std::vector<float> >("PhoCone04PhotonIsodR8dEta0pt0");
      produces<std::vector<float> >("PhoCone04PhotonIsodR8dEta0pt5");
      produces<std::vector<float> >("PhoCone01PhotonIsodR045EB070EEdEta015pt08EB1EEmvVtx");
      produces<std::vector<float> >("PhoCone02PhotonIsodR045EB070EEdEta015pt08EB1EEmvVtx");
      produces<std::vector<float> >("PhoCone03PhotonIsodR045EB070EEdEta015pt08EB1EEmvVtx");
      produces<std::vector<float> >("PhoCone04PhotonIsodR045EB070EEdEta015pt08EB1EEmvVtx");
      produces<std::vector<float> >("PhoCone04NeutralHadronIsodR0dEta0pt0");
      produces<std::vector<float> >("PhoCone04NeutralHadronIsodR0dEta0pt5");
      produces<std::vector<float> >("PhoCone04NeutralHadronIsodR0dEta0pt0nocracks");
      produces<std::vector<float> >("PhoCone04NeutralHadronIsodR0dEta0pt5nocracks");
      produces<std::vector<float> >("PhoCone04NeutralHadronIsodR7dEta0pt0");
      produces<std::vector<float> >("PhoCone04NeutralHadronIsodR7dEta0pt5");
      produces<std::vector<float> >("PhoCone01NeutralHadronIsodR0dEta0pt0mvVtx");
      produces<std::vector<float> >("PhoCone02NeutralHadronIsodR0dEta0pt0mvVtx");
      produces<std::vector<float> >("PhoCone03NeutralHadronIsodR0dEta0pt0mvVtx");
      produces<std::vector<float> >("PhoCone04NeutralHadronIsodR0dEta0pt0mvVtx");
      produces<std::vector<float> >("PhoCone04ChargedHadronIsodR0dEta0pt0dz0old");
      produces<std::vector<float> >("PhoCone04ChargedHadronIsodR0dEta0pt0PFnoPUold");
      produces<std::vector<float> >("PhoCone04ChargedHadronIsodR015dEta0pt0dz0old");
      produces<std::vector<float> >("PhoCone04ChargedHadronIsodR015dEta0pt0PFnoPUold");
      produces<std::vector<float> >("PhoCone01ChargedHadronIsodR0dEta0pt0dz0");
      produces<std::vector<float> >("PhoCone01ChargedHadronIsodR0dEta0pt0dz1dxy01");
      produces<std::vector<float> >("PhoCone01ChargedHadronIsodR0dEta0pt0PFnoPU");
      produces<std::vector<float> >("PhoCone01ChargedHadronIsodR015dEta0pt0dz0");
      produces<std::vector<float> >("PhoCone01ChargedHadronIsodR015dEta0pt0dz1dxy01");
      produces<std::vector<float> >("PhoCone01ChargedHadronIsodR015dEta0pt0PFnoPU");
      produces<std::vector<float> >("PhoCone02ChargedHadronIsodR0dEta0pt0dz0");
      produces<std::vector<float> >("PhoCone02ChargedHadronIsodR0dEta0pt0dz1dxy01");
      produces<std::vector<float> >("PhoCone02ChargedHadronIsodR0dEta0pt0PFnoPU");
      produces<std::vector<float> >("PhoCone02ChargedHadronIsodR015dEta0pt0dz0");
      produces<std::vector<float> >("PhoCone02ChargedHadronIsodR015dEta0pt0dz1dxy01");
      produces<std::vector<float> >("PhoCone02ChargedHadronIsodR015dEta0pt0PFnoPU");
      produces<std::vector<float> >("PhoCone03ChargedHadronIsodR0dEta0pt0dz0");
      produces<std::vector<float> >("PhoCone03ChargedHadronIsodR0dEta0pt0dz1dxy01");
      produces<std::vector<float> >("PhoCone03ChargedHadronIsodR0dEta0pt0PFnoPU");
      produces<std::vector<float> >("PhoCone03ChargedHadronIsodR015dEta0pt0dz0");
      produces<std::vector<float> >("PhoCone03ChargedHadronIsodR015dEta0pt0dz1dxy01");
      produces<std::vector<float> >("PhoCone03ChargedHadronIsodR015dEta0pt0PFnoPU");
      produces<std::vector<float> >("PhoCone04ChargedHadronIsodR0dEta0pt0dz0");
      produces<std::vector<float> >("PhoCone04ChargedHadronIsodR0dEta0pt0dz1dxy01");
      produces<std::vector<float> >("PhoCone04ChargedHadronIsodR0dEta0pt0PFnoPU");
      produces<std::vector<float> >("PhoCone04ChargedHadronIsodR015dEta0pt0dz0");
      produces<std::vector<float> >("PhoCone04ChargedHadronIsodR015dEta0pt0dz1dxy01");
      produces<std::vector<float> >("PhoCone04ChargedHadronIsodR015dEta0pt0PFnoPU");
      produces<std::vector<bool> > ("PhoConvValidVtx");
      produces<std::vector<int> >  ("PhoConvNtracks");
      produces<std::vector<float> >("PhoConvChi2Probability");
      produces<std::vector<float> >("PhoConvEoverP");
      produces<int>("Nconv");
      produces<std::vector<bool> >("ConvValidVtx");
      produces<std::vector<int> >("ConvNtracks");
      produces<std::vector<float> >("ConvChi2Probability");
      produces<std::vector<float> >("ConvEoverP");
      produces<std::vector<float> >("ConvZofPrimVtxFromTrks");
      //     produces<int>("Ngv");
      //     produces<std::vector<float> >("gvSumPtHi");
      //     produces<std::vector<float> >("gvSumPtLo");
      //     produces<std::vector<int> >("gvNTkHi");
      //     produces<std::vector<int> >("gvNTkLo");
      produces<int>("NSuperClusters");
      produces<std::vector<float> >("SCRaw");
      produces<std::vector<float> >("SCPre");
      produces<std::vector<float> >("SCEnergy");
      produces<std::vector<float> >("SCEta");
      produces<std::vector<float> >("SCPhi");
      produces<std::vector<float> >("SCPhiWidth");
      produces<std::vector<float> >("SCEtaWidth");
      produces<std::vector<float> >("SCBrem");
      produces<std::vector<float> >("SCR9");
      produces<std::vector<float> >("SCcrackcorrseed");
      produces<std::vector<float> >("SCcrackcorr");
      produces<std::vector<float> >("SClocalcorrseed");
      produces<std::vector<float> >("SClocalcorr");
      produces<std::vector<float> >("SCcrackcorrseedfactor");
      produces<std::vector<float> >("SClocalcorrseedfactor");
      produces<int>("NJets");
      produces<int>("NJetsTot");
      produces<std::vector<int> >("JGood");
      produces<std::vector<float> >("JPx");
      produces<std::vector<float> >("JPy");
      produces<std::vector<float> >("JPz");
      produces<std::vector<float> >("JPt");
      produces<std::vector<float> >("JE");
      produces<std::vector<float> >("JEt");
      produces<std::vector<float> >("JEta");
      produces<std::vector<float> >("JPhi");
      produces<std::vector<float> >("JEcorr");
      produces<std::vector<float> >("JArea");
      produces<std::vector<float> >("JEtaRms");
      produces<std::vector<float> >("JPhiRms");
      produces<std::vector<int> >("JNConstituents");
      produces<std::vector<int> >("JNAssoTracks");
      produces<std::vector<int> >("JNNeutrals");
      produces<std::vector<float> >("JChargedEmFrac");
      produces<std::vector<float> >("JNeutralEmFrac");
      produces<std::vector<float> >("JChargedHadFrac");
      produces<std::vector<float> >("JNeutralHadFrac");
      produces<std::vector<float> >("JChargedMuEnergyFrac");
      produces<std::vector<float> >("JPhoFrac");
      produces<std::vector<float> >("JHFHadFrac");
      produces<std::vector<float> >("JHFEMFrac");
      produces<std::vector<float> >("JPtD");
      produces<std::vector<float> >("JRMSCand");
      produces<std::vector<float> >("JeMinDR");
      produces<std::vector<float> >("JbTagProbTkCntHighEff");
      produces<std::vector<float> >("JbTagProbTkCntHighPur");
      produces<std::vector<float> >("JbTagProbSimpSVHighEff");
      produces<std::vector<float> >("JbTagProbSimpSVHighPur");
      produces<std::vector<float> >("JMass");
      produces<std::vector<float> >("Jtrk1px");
      produces<std::vector<float> >("Jtrk1py");
      produces<std::vector<float> >("Jtrk1pz");
      produces<std::vector<float> >("Jtrk2px");
      produces<std::vector<float> >("Jtrk2py");
      produces<std::vector<float> >("Jtrk2pz");
      produces<std::vector<float> >("Jtrk3px");
      produces<std::vector<float> >("Jtrk3py");
      produces<std::vector<float> >("Jtrk3pz");
      produces<std::vector<float> >("JVtxx");
      produces<std::vector<float> >("JVtxy");
      produces<std::vector<float> >("JVtxz");
      produces<std::vector<float> >("JVtxExx");
      produces<std::vector<float> >("JVtxEyx");
      produces<std::vector<float> >("JVtxEyy");
      produces<std::vector<float> >("JVtxEzy");
      produces<std::vector<float> >("JVtxEzz");
      produces<std::vector<float> >("JVtxEzx");
      produces<std::vector<float> >("JVtxNChi2");
      produces<std::vector<int> >("JGenJetIndex");
      produces<int>("NTracks");
      produces<int>("NTracksTot");
      produces<std::vector<int> >("TrkGood");
      produces<std::vector<float> >("TrkPt");
      produces<std::vector<float> >("TrkEta");
      produces<std::vector<float> >("TrkPhi");
      produces<std::vector<float> >("TrkNChi2");
      produces<std::vector<float> >("TrkNHits");
      produces<std::vector<float> >("TrkVtxDz");
      produces<std::vector<float> >("TrkVtxDxy");
      produces<float>("TrkPtSumx");
      produces<float>("TrkPtSumy");
      produces<float>("TrkPtSum");
      produces<float>("TrkPtSumPhi");
      produces<float>("SumEt");
      produces<float>("ECALSumEt");
      produces<float>("HCALSumEt");
      produces<float>("ECALEsumx");
      produces<float>("ECALEsumy");
      produces<float>("ECALEsumz");
      produces<float>("ECALMET");
      produces<float>("ECALMETPhi");
      produces<float>("ECALMETEta");
      produces<float>("HCALEsumx");
      produces<float>("HCALEsumy");
      produces<float>("HCALEsumz");
      produces<float>("HCALMET");
      produces<float>("HCALMETPhi");
      produces<float>("HCALMETeta");
      produces<float>("RawMET");
      produces<float>("RawMETpx");
      produces<float>("RawMETpy");
      produces<float>("RawMETphi");
      produces<float>("RawMETemEtFrac");
      produces<float>("RawMETemEtInEB");
      produces<float>("RawMETemEtInEE");
      produces<float>("RawMETemEtInHF");
      produces<float>("RawMEThadEtFrac");
      produces<float>("RawMEThadEtInHB");
      produces<float>("RawMEThadEtInHE");
      produces<float>("RawMEThadEtInHF");
      produces<float>("RawMETSignificance");
      produces<float>("GenMET");
      produces<float>("GenMETpx");
      produces<float>("GenMETpy");
      produces<float>("GenMETphi");
      produces<float>("TCMET");
      produces<float>("TCMETpx");
      produces<float>("TCMETpy");
      produces<float>("TCMETphi");
      produces<float>("TCMETSignificance");
      produces<float>("MuJESCorrMET");
      produces<float>("MuJESCorrMETpx");
      produces<float>("MuJESCorrMETpy");
      produces<float>("MuJESCorrMETphi");
      produces<float>("PFMET");
      produces<float>("PFMETpx");
      produces<float>("PFMETpy");
      produces<float>("PFMETphi");
      produces<float>("PFMETSignificance");
      produces<float>("PFSumEt");
      produces<float>("PFMETPAT");
      produces<float>("PFMETPATpx");
      produces<float>("PFMETPATpy");
      produces<float>("PFMETPATphi");
      produces<float>("PFMETPATSignificance");
      produces<float>("METR12");
      produces<float>("METR21");

    }

    //________________________________________________________________________________________
    // Reset all event variables
    void NTupleProducer::resetProducts( void ) {
  
      fTRun.reset(new int(-999));
      fTEvent.reset(new int(-999));
      fTLumiSection.reset(new int(-999));
      fTPtHat.reset(new float(-999.99));
      fTQCDPartonicHT.reset(new float(-999.99));
      fTSigProcID.reset(new int(-999));
      fTPDFScalePDF.reset(new float(-999.99));
      fTPDFID1.reset(new int(-999));
      fTPDFID2.reset(new int(-999));
      fTPDFx1.reset(new float(-999.99));
      fTPDFx2.reset(new float(-999.99));
      fTPDFxPDF1.reset(new float(-999.99));
      fTPDFxPDF2.reset(new float(-999.99));
      fTGenWeight.reset(new float(-999.99));
      fTpdfW.reset(new std::vector<float> );
      fTpdfWsum.reset(new float(-999.99));
      fTNPdfs.reset(new int(-999));
      fTPUnumInteractions.reset(new int(-999));
      fTPUnumTrueInteractions.reset(new int(-999));
      fTPUnumFilled.reset(new int(-999));
      fTPUOOTnumInteractionsEarly.reset(new int(-999));
      fTPUOOTnumInteractionsLate.reset(new int(-999));
      fTPUzPositions.reset(new std::vector<float> );
      fTPUsumPtLowPt.reset(new std::vector<float> );
      fTPUsumPtHighPt.reset(new std::vector<float> );
      fTPUnTrksLowPt.reset(new std::vector<float> );
      fTPUnTrksHighPt.reset(new std::vector<float> );
      fTRho.reset(new float(-999.99));
      fTRhoPFnoPU.reset(new float(-999.99));
      fTWeight.reset(new float(-999.99));
      fTHLTResults.reset(new std::vector<int> );
      fTHLTPrescale.reset(new std::vector<int> );
      fTL1PhysResults.reset(new std::vector<int> );
      fTL1TechResults.reset(new std::vector<int> );
      fTNHLTObjs.reset(new int(-999));
      for ( size_t i=0; i<gMaxHltNObjs; ++i ) {
        fTHLTObjectID[i].reset(new std::vector<int> );
        fTHLTObjectPt[i].reset(new std::vector<float> );
        fTHLTObjectEta[i].reset(new std::vector<float> );
        fTHLTObjectPhi[i].reset(new std::vector<float> );
      }
      fTPUWeightTotal.reset(new float(-999.99));
      fTPUWeightInTime.reset(new float(-999.99));
      fTMassGlu.reset(new float(-999.99));
      fTMassChi.reset(new float(-999.99));
      fTMassLSP.reset(new float(-999.99));
      fTxSMS.reset(new float(-999.99));
      fTxbarSMS.reset(new float(-999.99));
      fTM0.reset(new float(-999.99));
      fTM12.reset(new float(-999.99));
      fTsignMu.reset(new float(-999.99));
      fTA0.reset(new float(-999.99));
      fTprocess.reset(new int(-999));

      fTMaxGenPartExceed.reset(new int(-999));
      fTnGenParticles.reset(new int(-999));
      fTgenInfoId.reset(new std::vector<int>);
      fTgenInfoStatus.reset(new std::vector<int>);
      fTgenInfoNMo.reset(new std::vector<int>);
      fTgenInfoNDa.reset(new std::vector<int>);
      fTgenInfoMo1.reset(new std::vector<int>);
      fTgenInfoMo2.reset(new std::vector<int>);
      fTgenInfoDa1.reset(new std::vector<int>);
      fTgenInfoDa2.reset(new std::vector<int>);
      fTgenInfoMoIndex.reset(new std::vector<int>);
      fTPromptnessLevel.reset(new std::vector<int>);
      fTgenInfoMass.reset(new std::vector<float>);
      fTgenInfoMo1Pt.reset(new std::vector<float>);
      fTgenInfoMo2Pt.reset(new std::vector<float>);
      fTgenInfoPt.reset(new std::vector<float>);
      fTgenInfoEta.reset(new std::vector<float>);
      fTgenInfoPhi.reset(new std::vector<float>);
      fTgenInfoPx.reset(new std::vector<float>);
      fTgenInfoPy.reset(new std::vector<float>);
      fTgenInfoPz.reset(new std::vector<float>);
      fTgenInfoM.reset(new std::vector<float>);
      fTgenInfoPromptFlag.reset(new std::vector<float>);

      fTPrimVtxGood.reset(new int(-999));
      fTPrimVtxx.reset(new float(-999.99));
      fTPrimVtxy.reset(new float(-999.99));
      fTPrimVtxz.reset(new float(-999.99));
      fTPrimVtxRho.reset(new float(-999.99));
      fTPrimVtxxE.reset(new float(-999.99));
      fTPrimVtxyE.reset(new float(-999.99));
      fTPrimVtxzE.reset(new float(-999.99));
      fTPrimVtxNChi2.reset(new float(-999.99));
      fTPrimVtxNdof.reset(new float(-999.99));
      fTPrimVtxIsFake.reset(new int(-999));
      fTPrimVtxPtSum.reset(new float(-999.99));
      fTBeamspotx.reset(new float(-999.99));
      fTBeamspoty.reset(new float(-999.99));
      fTBeamspotz.reset(new float(-999.99));
      fTNCaloTowers.reset(new int(-999));
      fTGoodEvent.reset(new int(-999));
      fTMaxMuExceed.reset(new int(-999));
      fTMaxElExceed.reset(new int(-999));
      fTMaxJetExceed.reset(new int(-999));
      fTMaxUncJetExceed.reset(new int(-999));
      fTMaxTrkExceed.reset(new int(-999));
      fTMaxPhotonsExceed.reset(new int(-999));
      fTMaxGenLepExceed.reset(new int(-999));
      fTMaxGenPhoExceed.reset(new int(-999));
      fTMaxGenJetExceed.reset(new int(-999));
      fTMaxVerticesExceed.reset(new int(-999));
      fTHBHENoiseFlag.reset(new int(-999));
      fTHBHENoiseFlagIso.reset(new int(-999));
      fTCSCTightHaloID.reset(new int(-999));
      fTEcalDeadTPFilterFlag.reset(new int(-999));
      fTRecovRecHitFilterFlag.reset(new int(-999));
      fTRA2TrackingFailureFilterFlag.reset(new int(-999));
      //FR fPBNRFlag.reset(new int(-999));
      fTNGenLeptons.reset(new int(-999));
      fTGenLeptonID.reset(new std::vector<int> );
      fTGenLeptonPt.reset(new std::vector<float> );
      fTGenLeptonEta.reset(new std::vector<float> );
      fTGenLeptonPhi.reset(new std::vector<float> );
      fTGenLeptonMID.reset(new std::vector<int> );
      fTGenLeptonMStatus.reset(new std::vector<int> );
      fTGenLeptonMPt.reset(new std::vector<float> );
      fTGenLeptonMEta.reset(new std::vector<float> );
      fTGenLeptonMPhi.reset(new std::vector<float> );
      fTGenLeptonGMID.reset(new std::vector<int> );
      fTGenLeptonGMStatus.reset(new std::vector<int> );
      fTGenLeptonGMPt.reset(new std::vector<float> );
      fTGenLeptonGMEta.reset(new std::vector<float> );
      fTGenLeptonGMPhi.reset(new std::vector<float> );
      fTNGenPhotons.reset(new int(-999));
      fTGenPhotonPt.reset(new std::vector<float> );
      fTGenPhotonEta.reset(new std::vector<float> );
      fTGenPhotonPhi.reset(new std::vector<float> );
      fTGenPhotonPartonMindR.reset(new std::vector<float> );
      fTGenPhotonMotherID.reset(new std::vector<int> );
      fTGenPhotonMotherStatus.reset(new std::vector<int> );
      fTNGenJets.reset(new int(-999));
      fTGenJetPt.reset(new std::vector<float> );
      fTGenJetEta.reset(new std::vector<float> );
      fTGenJetPhi.reset(new std::vector<float> );
      fTGenJetE.reset(new std::vector<float> );
      fTGenJetEmE.reset(new std::vector<float> );
      fTGenJetHadE.reset(new std::vector<float> );
      fTGenJetInvE.reset(new std::vector<float> );
      fTNVrtx.reset(new int(-999));
      fTVrtxX.reset(new std::vector<float> );
      fTVrtxY.reset(new std::vector<float> );
      fTVrtxZ.reset(new std::vector<float> );
      fTVrtxXE.reset(new std::vector<float> );
      fTVrtxYE.reset(new std::vector<float> );
      fTVrtxZE.reset(new std::vector<float> );
      fTVrtxNdof.reset(new std::vector<float> );
      fTVrtxChi2.reset(new std::vector<float> );
      fTVrtxNtrks.reset(new std::vector<float> );
      fTVrtxSumPt.reset(new std::vector<float> );
      fTVrtxIsFake.reset(new std::vector<int> );
      fTNMus.reset(new int(-999));
      fTNMusTot.reset(new int(-999));
      fTNGMus.reset(new int(-999));
      fTNTMus.reset(new int(-999));
      fTMuGood.reset(new std::vector<int> );
      fTMuIsIso.reset(new std::vector<int> );
      fTMuIsGlobalMuon.reset(new std::vector<int> );
      fTMuIsTrackerMuon.reset(new std::vector<int> );
      fTMuPx.reset(new std::vector<float> );
      fTMuPy.reset(new std::vector<float> );
      fTMuPz.reset(new std::vector<float> );
      fTMuPt.reset(new std::vector<float> );
      fTMuInnerTkPt.reset(new std::vector<float> );
      fTMuPtE.reset(new std::vector<float> );
      fTMuE.reset(new std::vector<float> );
      fTMuEt.reset(new std::vector<float> );
      fTMuEta.reset(new std::vector<float> );
      fTMuPhi.reset(new std::vector<float> );
      fTMuCharge.reset(new std::vector<int> );
      fTMuRelIso03.reset(new std::vector<float> );
      fTMuIso03SumPt.reset(new std::vector<float> );
      fTMuIso03EmEt.reset(new std::vector<float> );
      fTMuIso03HadEt.reset(new std::vector<float> );
      fTMuIso03EMVetoEt.reset(new std::vector<float> );
      fTMuIso03HadVetoEt.reset(new std::vector<float> );
      fTMuIso05SumPt.reset(new std::vector<float> );
      fTMuIso05EmEt.reset(new std::vector<float> );
      fTMuIso05HadEt.reset(new std::vector<float> );
      fTMuEem.reset(new std::vector<float> );
      fTMuEhad.reset(new std::vector<float> );
      fTMuD0BS.reset(new std::vector<float> );
      fTMuD0PV.reset(new std::vector<float> );
      fTMuD0E.reset(new std::vector<float> );
      fTMuDzBS.reset(new std::vector<float> );
      fTMuDzPV.reset(new std::vector<float> );
      fTMuDzE.reset(new std::vector<float> );
      fTMuNChi2.reset(new std::vector<float> );
      fTMuNGlHits.reset(new std::vector<int> );
      fTMuNMuHits.reset(new std::vector<int> );
      fTMuNTkHits.reset(new std::vector<int> );
      fTMuNPxHits.reset(new std::vector<int> );
      fTMuInnerTkNChi2.reset(new std::vector<float> );
      fTMuNMatches.reset(new std::vector<int> );
      fTMuNChambers.reset(new std::vector<int> );
      fTMuCaloComp.reset(new std::vector<float> );
      fTMuSegmComp.reset(new std::vector<float> );
      fTMuIsGMPT.reset(new std::vector<int> );
      fTMuIsGMTkChiComp.reset(new std::vector<int> );
      fTMuIsGMStaChiComp.reset(new std::vector<int> );
      fTMuIsGMTkKinkTight.reset(new std::vector<int> );
      fTMuIsAllStaMuons.reset(new std::vector<int> );
      fTMuIsAllTrkMuons.reset(new std::vector<int> );
      fTMuIsTrkMuonArbitrated.reset(new std::vector<int> );
      fTMuIsAllArbitrated.reset(new std::vector<int> );
      fTMuIsTMLSLoose.reset(new std::vector<int> );
      fTMuIsTMLSTight.reset(new std::vector<int> );
      fTMuIsTM2DCompLoose.reset(new std::vector<int> );
      fTMuIsTM2DCompTight.reset(new std::vector<int> );
      fTMuIsTMOneStationLoose.reset(new std::vector<int> );
      fTMuIsTMOneStationTight.reset(new std::vector<int> );
      fTMuIsTMLSOptLowPtLoose.reset(new std::vector<int> );
      fTMuIsTMLSAngLoose.reset(new std::vector<int> );
      fTMuIsTMLSAngTight.reset(new std::vector<int> );
      fTMuIsTMOneStationAngTight.reset(new std::vector<int> );
      fTMuIsTMOneStationAngLoose.reset(new std::vector<int> );
      fTMuGenID.reset(new std::vector<int> );
      fTMuGenStatus.reset(new std::vector<int> );
      fTMuGenPt.reset(new std::vector<float> );
      fTMuGenEta.reset(new std::vector<float> );
      fTMuGenPhi.reset(new std::vector<float> );
      fTMuGenE.reset(new std::vector<float> );
      fTMuGenMID.reset(new std::vector<int> );
      fTMuGenMStatus.reset(new std::vector<int> );
      fTMuGenMPt.reset(new std::vector<float> );
      fTMuGenMEta.reset(new std::vector<float> );
      fTMuGenMPhi.reset(new std::vector<float> );
      fTMuGenME.reset(new std::vector<float> );
      fTMuGenGMID.reset(new std::vector<int> );
      fTMuGenGMStatus.reset(new std::vector<int> );
      fTMuGenGMPt.reset(new std::vector<float> );
      fTMuGenGMEta.reset(new std::vector<float> );
      fTMuGenGMPhi.reset(new std::vector<float> );
      fTMuGenGME.reset(new std::vector<float> );
      fTNEBhits.reset(new int(-999));
      fTEBrechitE.reset(new std::vector<float> );
      fTEBrechitPt.reset(new std::vector<float> );
      fTEBrechitEta.reset(new std::vector<float> );
      fTEBrechitPhi.reset(new std::vector<float> );
      fTEBrechitChi2.reset(new std::vector<float> );
      fTEBrechitTime.reset(new std::vector<float> );
      fTEBrechitE4oE1.reset(new std::vector<float> );
      fTEBrechitE2oE9.reset(new std::vector<float> );
      fTNEles.reset(new int(-999));
      fTNElesTot.reset(new int(-999));
      fTElGood.reset(new std::vector<int> );
      fTElIsIso.reset(new std::vector<int> );
      fTElChargeMisIDProb.reset(new std::vector<int> );
      fTElPx.reset(new std::vector<float> );
      fTElPy.reset(new std::vector<float> );
      fTElPz.reset(new std::vector<float> );
      fTElPt.reset(new std::vector<float> );
      fTElPtE.reset(new std::vector<float> );
      fTElE.reset(new std::vector<float> );
      fTElEt.reset(new std::vector<float> );
      fTElEta.reset(new std::vector<float> );
      fTElTheta.reset(new std::vector<float> );
      fTElSCEta.reset(new std::vector<float> );
      fTElPhi.reset(new std::vector<float> );
      fTElGsfTkPt.reset(new std::vector<float> );
      fTElGsfTkEta.reset(new std::vector<float> );
      fTElGsfTkPhi.reset(new std::vector<float> );
      fTElTrkMomentumError.reset(new std::vector<float> );
      fTElEcalEnergyError.reset(new std::vector<float> );
      fTElEleMomentumError.reset(new std::vector<float> );
      fTElNBrems.reset(new std::vector<int> );
      fTElD0BS.reset(new std::vector<float> );
      fTElD0PV.reset(new std::vector<float> );
      fTElD0E.reset(new std::vector<float> );
      fTElDzBS.reset(new std::vector<float> );
      fTElDzPV.reset(new std::vector<float> );
      fTElDzE.reset(new std::vector<float> );
      fTElRelIso03.reset(new std::vector<float> );
      fTElRelIso04.reset(new std::vector<float> );
      fTElDR03TkSumPt.reset(new std::vector<float> );
      fTElDR04TkSumPt.reset(new std::vector<float> );
      fTElDR03EcalRecHitSumEt.reset(new std::vector<float> );
      fTElDR04EcalRecHitSumEt.reset(new std::vector<float> );
      fTElDR03HcalTowerSumEt.reset(new std::vector<float> );
      fTElDR04HcalTowerSumEt.reset(new std::vector<float> );
      fTElNChi2.reset(new std::vector<float> );
      fTElCharge.reset(new std::vector<int> );
      fTElCInfoIsGsfCtfCons.reset(new std::vector<int> );
      fTElCInfoIsGsfCtfScPixCons.reset(new std::vector<int> );
      fTElCInfoIsGsfScPixCons.reset(new std::vector<int> );
      fTElScPixCharge.reset(new std::vector<int> );
      fTElClosestCtfTrackPt.reset(new std::vector<float> );
      fTElClosestCtfTrackEta.reset(new std::vector<float> );
      fTElClosestCtfTrackPhi.reset(new std::vector<float> );
      fTElClosestCtfTrackCharge.reset(new std::vector<int> );
      fTElIDMva.reset(new std::vector<float> );
      fTElIDTight.reset(new std::vector<int> );
      fTElIDLoose.reset(new std::vector<int> );
      fTElIDRobustTight.reset(new std::vector<int> );
      fTElIDRobustLoose.reset(new std::vector<int> );
      fTElIDsimpleWPrelIso.reset(new std::vector<int> );
      fTElIDsimpleWP80relIso.reset(new std::vector<int> );
      fTElIDsimpleWP85relIso.reset(new std::vector<int> );
      fTElIDsimpleWP90relIso.reset(new std::vector<int> );
      fTElIDsimpleWP95relIso.reset(new std::vector<int> );
      fTElInGap.reset(new std::vector<int> );
      fTElEcalDriven.reset(new std::vector<int> );
      fTElTrackerDriven.reset(new std::vector<int> );
      fTElBasicClustersSize.reset(new std::vector<int> );
      fTElfbrem.reset(new std::vector<float> );
      fTElHcalOverEcal.reset(new std::vector<float> );
      fTElE1x5.reset(new std::vector<float> );
      fTElE5x5.reset(new std::vector<float> );
      fTElE2x5Max.reset(new std::vector<float> );
      fTElSigmaIetaIeta.reset(new std::vector<float> );
      fTElDeltaPhiSeedClusterAtCalo.reset(new std::vector<float> );
      fTElDeltaEtaSeedClusterAtCalo.reset(new std::vector<float> );
      fTElDeltaPhiSuperClusterAtVtx.reset(new std::vector<float> );
      fTElDeltaEtaSuperClusterAtVtx.reset(new std::vector<float> );
      fTElCaloEnergy.reset(new std::vector<float> );
      fTElTrkMomAtVtx.reset(new std::vector<float> );
      fTElESuperClusterOverP.reset(new std::vector<float> );
      fTElNumberOfMissingInnerHits.reset(new std::vector<int> );
      fTElSCindex.reset(new std::vector<int> );
      fTElConvPartnerTrkDist.reset(new std::vector<float> );
      fTElConvPartnerTrkDCot.reset(new std::vector<float> );
      fTElConvPartnerTrkPt.reset(new std::vector<float> );
      fTElConvPartnerTrkEta.reset(new std::vector<float> );
      fTElConvPartnerTrkPhi.reset(new std::vector<float> );
      fTElConvPartnerTrkCharge.reset(new std::vector<float> );
      fTElScSeedSeverity.reset(new std::vector<int> );
      fTElE1OverE9.reset(new std::vector<float> );
      fTElS4OverS1.reset(new std::vector<float> );
      fTElGenID.reset(new std::vector<int> );
      fTElGenStatus.reset(new std::vector<int> );
      fTElGenPt.reset(new std::vector<float> );
      fTElGenEta.reset(new std::vector<float> );
      fTElGenPhi.reset(new std::vector<float> );
      fTElGenE.reset(new std::vector<float> );
      fTElGenMID.reset(new std::vector<int> );
      fTElGenMStatus.reset(new std::vector<int> );
      fTElGenMPt.reset(new std::vector<float> );
      fTElGenMEta.reset(new std::vector<float> );
      fTElGenMPhi.reset(new std::vector<float> );
      fTElGenME.reset(new std::vector<float> );
      fTElGenGMID.reset(new std::vector<int> );
      fTElGenGMStatus.reset(new std::vector<int> );
      fTElGenGMPt.reset(new std::vector<float> );
      fTElGenGMEta.reset(new std::vector<float> );
      fTElGenGMPhi.reset(new std::vector<float> );
      fTElGenGME.reset(new std::vector<float> );
      fTNPhotons.reset(new int(-999));
      fTNPhotonsTot.reset(new int(-999));
      fTPhoGood.reset(new std::vector<int> );
      fTPhoIsIso.reset(new std::vector<int> );
      fTPhoPt.reset(new std::vector<float> );
      fTPhoPx.reset(new std::vector<float> );
      fTPhoPy.reset(new std::vector<float> );
      fTPhoPz.reset(new std::vector<float> );
      fTPhoEta.reset(new std::vector<float> );
      fTPhoPhi.reset(new std::vector<float> );
      fTPhoEnergy.reset(new std::vector<float> );
      fTPhoIso03Ecal.reset(new std::vector<float> );
      fTPhoIso03Hcal.reset(new std::vector<float> );
      fTPhoIso03TrkSolid.reset(new std::vector<float> );
      fTPhoIso03TrkHollow.reset(new std::vector<float> );
      fTPhoIso03.reset(new std::vector<float> );
      fTPhoIso04Ecal.reset(new std::vector<float> );
      fTPhoIso04Hcal.reset(new std::vector<float> );
      fTPhoIso04TrkSolid.reset(new std::vector<float> );
      fTPhoIso04TrkHollow.reset(new std::vector<float> );
      fTPhoIso04.reset(new std::vector<float> );
      fTPhoR9.reset(new std::vector<float> );
      fTPhoCaloPositionX.reset(new std::vector<float> );
      fTPhoCaloPositionY.reset(new std::vector<float> );
      fTPhoCaloPositionZ.reset(new std::vector<float> );
      fTPhoHoverE.reset(new std::vector<float> );
      fTPhoH1overE.reset(new std::vector<float> );
      fTPhoH2overE.reset(new std::vector<float> );
      fTPhoSigmaIetaIeta.reset(new std::vector<float> );
      fTPhoSCRawEnergy.reset(new std::vector<float> );
      fTPhoSCEtaWidth.reset(new std::vector<float> );
      fTPhoSCSigmaPhiPhi.reset(new std::vector<float> );
      fTPhoHasPixSeed.reset(new std::vector<int> );
      fTPhoHasConvTrks.reset(new std::vector<int> );
      fTPhoScSeedSeverity.reset(new std::vector<int> );
      fTPhoE1OverE9.reset(new std::vector<float> );
      fTPhoS4OverS1.reset(new std::vector<float> );
      fTPhoSigmaEtaEta.reset(new std::vector<float> );
      fTPhoE1x5.reset(new std::vector<float> );
      fTPhoE2x5.reset(new std::vector<float> );
      fTPhoE3x3.reset(new std::vector<float> );
      fTPhoE5x5.reset(new std::vector<float> );
      fTPhomaxEnergyXtal.reset(new std::vector<float> );
      fTPhoIso03HcalDepth1.reset(new std::vector<float> );
      fTPhoIso03HcalDepth2.reset(new std::vector<float> );
      fTPhoIso04HcalDepth1.reset(new std::vector<float> );
      fTPhoIso04HcalDepth2.reset(new std::vector<float> );
      fTPhoIso03nTrksSolid.reset(new std::vector<int> );
      fTPhoIso03nTrksHollow.reset(new std::vector<int> );
      fTPhoIso04nTrksSolid.reset(new std::vector<int> );
      fTPhoIso04nTrksHollow.reset(new std::vector<int> );
      fTPhoisEB.reset(new std::vector<int> );
      fTPhoisEE.reset(new std::vector<int> );
      fTPhoisEBEtaGap.reset(new std::vector<int> );
      fTPhoisEBPhiGap.reset(new std::vector<int> );
      fTPhoisEERingGap.reset(new std::vector<int> );
      fTPhoisEEDeeGap.reset(new std::vector<int> );
      fTPhoisEBEEGap.reset(new std::vector<int> );
      fTPhoisPFlowPhoton.reset(new std::vector<int> );
      fTPhoisStandardPhoton.reset(new std::vector<int> );
      fTPhoMCmatchindex.reset(new std::vector<int> );
      fTPhoMCmatchexitcode.reset(new std::vector<int> );
      fTPhoChargedHadronIso.reset(new std::vector<float> );
      fTPhoNeutralHadronIso.reset(new std::vector<float> );
      fTPhoPhotonIso.reset(new std::vector<float> );
      fTPhoisPFPhoton.reset(new std::vector<int> );
      fTPhoisPFElectron.reset(new std::vector<int> );
      fTPhotSCindex.reset(new std::vector<int> );
      fTPhoCone04PhotonIsodR0dEta0pt0.reset(new std::vector<float> );
      fTPhoCone04PhotonIsodR0dEta0pt5.reset(new std::vector<float> );
      fTPhoCone04PhotonIsodR8dEta0pt0.reset(new std::vector<float> );
      fTPhoCone04PhotonIsodR8dEta0pt5.reset(new std::vector<float> );
      fTPhoCone01PhotonIsodR045EB070EEdEta015pt08EB1EEmvVtx.reset(new std::vector<float> );
      fTPhoCone02PhotonIsodR045EB070EEdEta015pt08EB1EEmvVtx.reset(new std::vector<float> );
      fTPhoCone03PhotonIsodR045EB070EEdEta015pt08EB1EEmvVtx.reset(new std::vector<float> );
      fTPhoCone04PhotonIsodR045EB070EEdEta015pt08EB1EEmvVtx.reset(new std::vector<float> );
      fTPhoCone04NeutralHadronIsodR0dEta0pt0.reset(new std::vector<float> );
      fTPhoCone04NeutralHadronIsodR0dEta0pt5.reset(new std::vector<float> );
      fTPhoCone04NeutralHadronIsodR0dEta0pt0nocracks.reset(new std::vector<float> );
      fTPhoCone04NeutralHadronIsodR0dEta0pt5nocracks.reset(new std::vector<float> );
      fTPhoCone04NeutralHadronIsodR7dEta0pt0.reset(new std::vector<float> );
      fTPhoCone04NeutralHadronIsodR7dEta0pt5.reset(new std::vector<float> );
      fTPhoCone01NeutralHadronIsodR0dEta0pt0mvVtx.reset(new std::vector<float> );
      fTPhoCone02NeutralHadronIsodR0dEta0pt0mvVtx.reset(new std::vector<float> );
      fTPhoCone03NeutralHadronIsodR0dEta0pt0mvVtx.reset(new std::vector<float> );
      fTPhoCone04NeutralHadronIsodR0dEta0pt0mvVtx.reset(new std::vector<float> );
      fTPhoCone04ChargedHadronIsodR0dEta0pt0dz0old.reset(new std::vector<float> );
      fTPhoCone04ChargedHadronIsodR0dEta0pt0PFnoPUold.reset(new std::vector<float> );
      fTPhoCone04ChargedHadronIsodR015dEta0pt0dz0old.reset(new std::vector<float> );
      fTPhoCone04ChargedHadronIsodR015dEta0pt0PFnoPUold.reset(new std::vector<float> );
      fTPhoCone01ChargedHadronIsodR0dEta0pt0dz0.reset(new std::vector<float> );
      fTPhoCone01ChargedHadronIsodR0dEta0pt0dz1dxy01.reset(new std::vector<float> );
      fTPhoCone01ChargedHadronIsodR0dEta0pt0PFnoPU.reset(new std::vector<float> );
      fTPhoCone01ChargedHadronIsodR015dEta0pt0dz0.reset(new std::vector<float> );
      fTPhoCone01ChargedHadronIsodR015dEta0pt0dz1dxy01.reset(new std::vector<float> );
      fTPhoCone01ChargedHadronIsodR015dEta0pt0PFnoPU.reset(new std::vector<float> );
      fTPhoCone02ChargedHadronIsodR0dEta0pt0dz0.reset(new std::vector<float> );
      fTPhoCone02ChargedHadronIsodR0dEta0pt0dz1dxy01.reset(new std::vector<float> );
      fTPhoCone02ChargedHadronIsodR0dEta0pt0PFnoPU.reset(new std::vector<float> );
      fTPhoCone02ChargedHadronIsodR015dEta0pt0dz0.reset(new std::vector<float> );
      fTPhoCone02ChargedHadronIsodR015dEta0pt0dz1dxy01.reset(new std::vector<float> );
      fTPhoCone02ChargedHadronIsodR015dEta0pt0PFnoPU.reset(new std::vector<float> );
      fTPhoCone03ChargedHadronIsodR0dEta0pt0dz0.reset(new std::vector<float> );
      fTPhoCone03ChargedHadronIsodR0dEta0pt0dz1dxy01.reset(new std::vector<float> );
      fTPhoCone03ChargedHadronIsodR0dEta0pt0PFnoPU.reset(new std::vector<float> );
      fTPhoCone03ChargedHadronIsodR015dEta0pt0dz0.reset(new std::vector<float> );
      fTPhoCone03ChargedHadronIsodR015dEta0pt0dz1dxy01.reset(new std::vector<float> );
      fTPhoCone03ChargedHadronIsodR015dEta0pt0PFnoPU.reset(new std::vector<float> );
      fTPhoCone04ChargedHadronIsodR0dEta0pt0dz0.reset(new std::vector<float> );
      fTPhoCone04ChargedHadronIsodR0dEta0pt0dz1dxy01.reset(new std::vector<float> );
      fTPhoCone04ChargedHadronIsodR0dEta0pt0PFnoPU.reset(new std::vector<float> );
      fTPhoCone04ChargedHadronIsodR015dEta0pt0dz0.reset(new std::vector<float> );
      fTPhoCone04ChargedHadronIsodR015dEta0pt0dz1dxy01.reset(new std::vector<float> );
      fTPhoCone04ChargedHadronIsodR015dEta0pt0PFnoPU.reset(new std::vector<float> );
      fTPhoConvValidVtx.reset(new std::vector<bool>);
      fTPhoConvNtracks.reset(new std::vector<int>);
      fTPhoConvChi2Probability.reset(new std::vector<float>);
      fTPhoConvEoverP.reset(new std::vector<float>);
      fTNconv.reset(new int(-999));
      fTConvValidVtx.reset(new std::vector<bool>);
      fTConvNtracks.reset(new std::vector<int>);
      fTConvChi2Probability.reset(new std::vector<float>);
      fTConvEoverP.reset(new std::vector<float>);
      fTConvZofPrimVtxFromTrks.reset(new std::vector<float>);
      //     fTNgv.reset(new int(-999));
      //     fTgvSumPtHi.reset(new std::vector<float>);
      //     fTgvSumPtLo.reset(new std::vector<float>);
      //     fTgvNTkHi.reset(new std::vector<int>);
      //     fTgvNTkLo.reset(new std::vector<int>);
      for (int i=0; i<gMaxNPhotons; i++) {
        pho_conv_vtx[i]=TVector3();
        pho_conv_refitted_momentum[i]=TVector3();
        conv_vtx[i]=TVector3();
        conv_refitted_momentum[i]=TVector3();
      }
      for (int i=0; i<gMaxNGenVtx; i++) {
        gv_pos[i]=TVector3();
        gv_p3[i]=TVector3();
      }
      fTNSuperClusters.reset(new int(-999));
      fTSCRaw.reset(new std::vector<float> );
      fTSCPre.reset(new std::vector<float> );
      fTSCEnergy.reset(new std::vector<float> );
      fTSCEta.reset(new std::vector<float> );
      fTSCPhi.reset(new std::vector<float> );
      fTSCPhiWidth.reset(new std::vector<float> );
      fTSCEtaWidth.reset(new std::vector<float> );
      fTSCBrem.reset(new std::vector<float> );
      fTSCR9.reset(new std::vector<float> );
      fTSCcrackcorrseed.reset(new std::vector<float> );
      fTSCcrackcorr.reset(new std::vector<float> );
      fTSClocalcorrseed.reset(new std::vector<float> );
      fTSClocalcorr.reset(new std::vector<float> );
      fTSCcrackcorrseedfactor.reset(new std::vector<float> );
      fTSClocalcorrseedfactor.reset(new std::vector<float> );
      fTNJets.reset(new int(-999));
      fTNJetsTot.reset(new int(-999));
      fTJGood.reset(new std::vector<int> );
      fTJPx.reset(new std::vector<float> );
      fTJPy.reset(new std::vector<float> );
      fTJPz.reset(new std::vector<float> );
      fTJPt.reset(new std::vector<float> );
      fTJE.reset(new std::vector<float> );
      fTJEt.reset(new std::vector<float> );
      fTJEta.reset(new std::vector<float> );
      fTJPhi.reset(new std::vector<float> );
      fTJEcorr.reset(new std::vector<float> );
      fTJArea.reset(new std::vector<float> );
      fTJEtaRms.reset(new std::vector<float> );
      fTJPhiRms.reset(new std::vector<float> );
      fTJNConstituents.reset(new std::vector<int> );
      fTJNAssoTracks.reset(new std::vector<int> );
      fTJNNeutrals.reset(new std::vector<int> );
      fTJChargedEmFrac.reset(new std::vector<float> );
      fTJNeutralEmFrac.reset(new std::vector<float> );
      fTJChargedHadFrac.reset(new std::vector<float> );
      fTJNeutralHadFrac.reset(new std::vector<float> );
      fTJChargedMuEnergyFrac.reset(new std::vector<float> );
      fTJPhoFrac.reset(new std::vector<float>);
      fTJHFHadFrac.reset(new std::vector<float>);
      fTJHFEMFrac.reset(new std::vector<float>);
      fTJPtD.reset(new std::vector<float>);
      fTJRMSCand.reset(new std::vector<float>);
      fTJeMinDR.reset(new std::vector<float> );
      fTJbTagProbTkCntHighEff.reset(new std::vector<float> );
      fTJbTagProbTkCntHighPur.reset(new std::vector<float> );
      fTJbTagProbSimpSVHighEff.reset(new std::vector<float> );
      fTJbTagProbSimpSVHighPur.reset(new std::vector<float> );
      fTJMass.reset(new std::vector<float> );
      fTJtrk1px.reset(new std::vector<float> );
      fTJtrk1py.reset(new std::vector<float> );
      fTJtrk1pz.reset(new std::vector<float> );
      fTJtrk2px.reset(new std::vector<float> );
      fTJtrk2py.reset(new std::vector<float> );
      fTJtrk2pz.reset(new std::vector<float> );
      fTJtrk3px.reset(new std::vector<float> );
      fTJtrk3py.reset(new std::vector<float> );
      fTJtrk3pz.reset(new std::vector<float> );
      fTJVtxx.reset(new std::vector<float> );
      fTJVtxy.reset(new std::vector<float> );
      fTJVtxz.reset(new std::vector<float> );
      fTJVtxExx.reset(new std::vector<float> );
      fTJVtxEyx.reset(new std::vector<float> );
      fTJVtxEyy.reset(new std::vector<float> );
      fTJVtxEzy.reset(new std::vector<float> );
      fTJVtxEzz.reset(new std::vector<float> );
      fTJVtxEzx.reset(new std::vector<float> );
      fTJVtxNChi2.reset(new std::vector<float> );
      fTJGenJetIndex.reset(new std::vector<int> );
      fTNTracks.reset(new int(-999));
      fTNTracksTot.reset(new int(-999));
      fTTrkGood.reset(new std::vector<int> );
      fTTrkPt.reset(new std::vector<float> );
      fTTrkEta.reset(new std::vector<float> );
      fTTrkPhi.reset(new std::vector<float> );
      fTTrkNChi2.reset(new std::vector<float> );
      fTTrkNHits.reset(new std::vector<float> );
      fTTrkVtxDz.reset(new std::vector<float> );
      fTTrkVtxDxy.reset(new std::vector<float> );
      fTTrkPtSumx.reset(new float(-999.99));
      fTTrkPtSumy.reset(new float(-999.99));
      fTTrkPtSum.reset(new float(-999.99));
      fTTrkPtSumPhi.reset(new float(-999.99));
      fTSumEt.reset(new float(-999.99));
      fTECALSumEt.reset(new float(-999.99));
      fTHCALSumEt.reset(new float(-999.99));
      fTECALEsumx.reset(new float(-999.99));
      fTECALEsumy.reset(new float(-999.99));
      fTECALEsumz.reset(new float(-999.99));
      fTECALMET.reset(new float(-999.99));
      fTECALMETPhi.reset(new float(-999.99));
      fTECALMETEta.reset(new float(-999.99));
      fTHCALEsumx.reset(new float(-999.99));
      fTHCALEsumy.reset(new float(-999.99));
      fTHCALEsumz.reset(new float(-999.99));
      fTHCALMET.reset(new float(-999.99));
      fTHCALMETPhi.reset(new float(-999.99));
      fTHCALMETeta.reset(new float(-999.99));
      fTRawMET.reset(new float(-999.99));
      fTRawMETpx.reset(new float(-999.99));
      fTRawMETpy.reset(new float(-999.99));
      fTRawMETphi.reset(new float(-999.99));
      fTRawMETemEtFrac.reset(new float(-999.99));
      fTRawMETemEtInEB.reset(new float(-999.99));
      fTRawMETemEtInEE.reset(new float(-999.99));
      fTRawMETemEtInHF.reset(new float(-999.99));
      fTRawMEThadEtFrac.reset(new float(-999.99));
      fTRawMEThadEtInHB.reset(new float(-999.99));
      fTRawMEThadEtInHE.reset(new float(-999.99));
      fTRawMEThadEtInHF.reset(new float(-999.99));
      fTRawMETSignificance.reset(new float(-999.99));
      fTGenMET.reset(new float(-999.99));
      fTGenMETpx.reset(new float(-999.99));
      fTGenMETpy.reset(new float(-999.99));
      fTGenMETphi.reset(new float(-999.99));
      fTTCMET.reset(new float(-999.99));
      fTTCMETpx.reset(new float(-999.99));
      fTTCMETpy.reset(new float(-999.99));
      fTTCMETphi.reset(new float(-999.99));
      fTTCMETSignificance.reset(new float(-999.99));
      fTMuJESCorrMET.reset(new float(-999.99));
      fTMuJESCorrMETpx.reset(new float(-999.99));
      fTMuJESCorrMETpy.reset(new float(-999.99));
      fTMuJESCorrMETphi.reset(new float(-999.99));
      fTPFMET.reset(new float(-999.99));
      fTPFMETpx.reset(new float(-999.99));
      fTPFMETpy.reset(new float(-999.99));
      fTPFMETphi.reset(new float(-999.99));
      fTPFMETSignificance.reset(new float(-999.99));
      fTPFSumEt.reset(new float(-999.99));
      fTPFMETPAT.reset(new float(-999.99));
      fTPFMETPATpx.reset(new float(-999.99));
      fTPFMETPATpy.reset(new float(-999.99));
      fTPFMETPATphi.reset(new float(-999.99));
      fTPFMETPATSignificance.reset(new float(-999.99));
      fTMETR12.reset(new float(-999.99));
      fTMETR21.reset(new float(-999.99));

    }

    //____________________________________________________________________
    void NTupleProducer::resetRunProducts( void ) {
  
      // Reset run variables
      fRExtXSecLO   .reset(new float(-999.99));
      fRExtXSecNLO  .reset(new float(-999.99));
      fRIntXSec     .reset(new float(-999.99));
                                         
      fRMinMuPt     .reset(new float(-999.99));
      fRMaxMuEta    .reset(new float(-999.99));
      fRMinElPt     .reset(new float(-999.99));
      fRMaxElEta    .reset(new float(-999.99));
      fRMinJPt      .reset(new float(-999.99));
      fRMinRawJPt   .reset(new float(-999.99));
      fRMaxJEta     .reset(new float(-999.99));
      fRMinJEMFrac  .reset(new float(-999.99));
                                         
      fRMinTrkPt    .reset(new float(-999.99));
      fRMaxTrkEta   .reset(new float(-999.99));
      fRMaxTrkNChi2 .reset(new float(-999.99));
      fRMinTrkNHits .reset(new int(-999));
                                         
      fRMinPhotonPt .reset(new float(-999.99));
      fRMaxPhotonEta.reset(new float(-999.99));
      fRMinSCraw    .reset(new float(-999.99));                                         
      fRMinEBRechitE.reset(new float(-999.99)); 

      fRMinGenLeptPt .reset(new float(-999.99)); 
      fRMaxGenLeptEta.reset(new float(-999.99)); 
      fRMinGenPhotPt .reset(new float(-999.99)); 
      fRMaxGenPhotEta.reset(new float(-999.99)); 
      fRMinGenJetPt  .reset(new float(-999.99)); 
      fRMaxGenJetEta .reset(new float(-999.99)); 
               
      fRMaxNMus     .reset(new int(-999));
      fRMaxNEles    .reset(new int(-999));
      fRMaxNJets    .reset(new int(-999));
      fRMaxNTrks    .reset(new int(-999));
      fRMaxNPhotons .reset(new int(-999));
      fRMaxNSC      .reset(new int(-999));
      fRMaxNGenLept .reset(new int(-999));
      fRMaxNGenPhot .reset(new int(-999));
      fRMaxNGenJets .reset(new int(-999));
      fRMaxNVrtx    .reset(new int(-999));
      fRMaxNPileup  .reset(new int(-999));
      fRMaxNEBhits  .reset(new int(-999));

      fRL1PhysMenu  .reset(new std::vector<std::string>);

      fRPileUpData  .reset(new std::vector<std::string>);
      fRPileUpMC    .reset(new std::vector<std::string>);

    }

    //________________________________________________________________________________________
    void NTupleProducer::putProducts( edm::Event& event ) {
  
      event.put(fTEvent, "Event");
      event.put(fTLumiSection, "LumiSection");
      event.put(fTPtHat, "PtHat");
      event.put(fTQCDPartonicHT, "QCDPartonicHT");
      event.put(fTSigProcID, "SigProcID");
      event.put(fTPDFScalePDF, "PDFScalePDF");
      event.put(fTPDFID1, "PDFID1");
      event.put(fTPDFID2, "PDFID2");
      event.put(fTPDFx1, "PDFx1");
      event.put(fTPDFx2, "PDFx2");
      event.put(fTPDFxPDF1, "PDFxPDF1");
      event.put(fTPDFxPDF2, "PDFxPDF2");
      event.put(fTGenWeight, "GenWeight");
      event.put(fTpdfW, "pdfW");
      event.put(fTpdfWsum, "pdfWsum");
      event.put(fTNPdfs, "NPdfs");
      event.put(fTPUnumInteractions, "PUnumInteractions");
      event.put(fTPUnumTrueInteractions, "PUnumTrueInteractions");
      event.put(fTPUnumFilled, "PUnumFilled");
      event.put(fTPUOOTnumInteractionsEarly, "PUOOTnumInteractionsEarly");
      event.put(fTPUOOTnumInteractionsLate, "PUOOTnumInteractionsLate");
      event.put(fTPUzPositions, "PUzPositions");
      event.put(fTPUsumPtLowPt, "PUsumPtLowPt");
      event.put(fTPUsumPtHighPt, "PUsumPtHighPt");
      event.put(fTPUnTrksLowPt, "PUnTrksLowPt");
      event.put(fTPUnTrksHighPt, "PUnTrksHighPt");
      event.put(fTRho, "Rho");
      event.put(fTRhoPFnoPU, "RhoPFnoPU");
      event.put(fTWeight, "Weight");
      event.put(fTHLTResults, "HLTResults");
      event.put(fTHLTPrescale, "HLTPrescale");
      event.put(fTL1PhysResults, "L1PhysResults");
      event.put(fTL1TechResults, "L1TechResults");
      event.put(fTNHLTObjs, "NHLTObjs");
      for ( size_t i=0; i<gMaxHltNObjs; ++i ) {
        std::ostringstream s;
        s << i;
        event.put(fTHLTObjectID[i], ("HLTObjectID"+s.str()).c_str());
        event.put(fTHLTObjectPt[i], ("HLTObjectPt"+s.str()).c_str());
        event.put(fTHLTObjectEta[i], ("HLTObjectEta"+s.str()).c_str());
        event.put(fTHLTObjectPhi[i], ("HLTObjectPhi"+s.str()).c_str());
      }
      event.put(fTPUWeightTotal, "PUWeightTotal");
      event.put(fTPUWeightInTime, "PUWeightInTime");
      event.put(fTMassGlu, "MassGlu");
      event.put(fTMassChi, "MassChi");
      event.put(fTMassLSP, "MassLSP");
      event.put(fTxSMS, "xSMS");
      event.put(fTxbarSMS, "xbarSMS");
      event.put(fTM0, "M0");
      event.put(fTM12, "M12");
      event.put(fTsignMu, "signMu");
      event.put(fTA0, "A0");
      event.put(fTprocess, "process");
      event.put(fTMaxGenPartExceed,"MaxGenPartExceed");
      event.put(fTnGenParticles,"nGenParticles");
      event.put(fTgenInfoId,"genInfoId");
      event.put(fTgenInfoStatus,"genInfoStatus");
      event.put(fTgenInfoNMo,"genInfoNMo");
      event.put(fTgenInfoNDa,"genInfoNDa");
      event.put(fTgenInfoMo1,"genInfoMo1");
      event.put(fTgenInfoMo2,"genInfoMo2");
      event.put(fTgenInfoDa1,"genInfoDa1");
      event.put(fTgenInfoDa2,"genInfoDa2");
      event.put(fTgenInfoMoIndex,"genInfoMoIndex");
      event.put(fTPromptnessLevel,"PromptnessLevel");
      event.put(fTgenInfoMass,"genInfoMass");
      event.put(fTgenInfoMo1Pt,"genInfoMo1Pt");
      event.put(fTgenInfoMo2Pt,"genInfoMo2Pt");
      event.put(fTgenInfoPt,"genInfoPt");
      event.put(fTgenInfoEta,"genInfoEta");
      event.put(fTgenInfoPhi,"genInfoPhi");
      event.put(fTgenInfoPx,"genInfoPx");
      event.put(fTgenInfoPy,"genInfoPy");
      event.put(fTgenInfoPz,"genInfoPz");
      event.put(fTgenInfoM,"genInfoM");
      event.put(fTgenInfoPromptFlag,"genInfoPromptFlag");
      event.put(fTPrimVtxGood, "PrimVtxGood");
      event.put(fTPrimVtxx, "PrimVtxx");
      event.put(fTPrimVtxy, "PrimVtxy");
      event.put(fTPrimVtxz, "PrimVtxz");
      event.put(fTPrimVtxRho, "PrimVtxRho");
      event.put(fTPrimVtxxE, "PrimVtxxE");
      event.put(fTPrimVtxyE, "PrimVtxyE");
      event.put(fTPrimVtxzE, "PrimVtxzE");
      event.put(fTPrimVtxNChi2, "PrimVtxNChi2");
      event.put(fTPrimVtxNdof, "PrimVtxNdof");
      event.put(fTPrimVtxIsFake, "PrimVtxIsFake");
      event.put(fTPrimVtxPtSum, "PrimVtxPtSum");
      event.put(fTBeamspotx, "Beamspotx");
      event.put(fTBeamspoty, "Beamspoty");
      event.put(fTBeamspotz, "Beamspotz");
      event.put(fTNCaloTowers, "NCaloTowers");
      event.put(fTGoodEvent, "GoodEvent");
      event.put(fTMaxMuExceed, "MaxMuExceed");
      event.put(fTMaxElExceed, "MaxElExceed");
      event.put(fTMaxJetExceed, "MaxJetExceed");
      event.put(fTMaxUncJetExceed, "MaxUncJetExceed");
      event.put(fTMaxTrkExceed, "MaxTrkExceed");
      event.put(fTMaxPhotonsExceed, "MaxPhotonsExceed");
      event.put(fTMaxGenLepExceed, "MaxGenLepExceed");
      event.put(fTMaxGenPhoExceed, "MaxGenPhoExceed");
      event.put(fTMaxGenJetExceed, "MaxGenJetExceed");
      event.put(fTMaxVerticesExceed, "MaxVerticesExceed");
      event.put(fTHBHENoiseFlag, "HBHENoiseFlag");
      event.put(fTHBHENoiseFlagIso, "HBHENoiseFlagIso");
      event.put(fTCSCTightHaloID, "CSCTightHaloID");
      event.put(fTEcalDeadTPFilterFlag, "EcalDeadTPFilterFlag");
      event.put(fTRecovRecHitFilterFlag, "RecovRecHitFilterFlag");
      event.put(fTRA2TrackingFailureFilterFlag,"RA2TrackingFailureFilterFlag");
      //FR event.put(fPBNRFlag,"PBNRFlag");
      event.put(fTNGenLeptons, "NGenLeptons");
      event.put(fTGenLeptonID, "GenLeptonID");
      event.put(fTGenLeptonPt, "GenLeptonPt");
      event.put(fTGenLeptonEta, "GenLeptonEta");
      event.put(fTGenLeptonPhi, "GenLeptonPhi");
      event.put(fTGenLeptonMID, "GenLeptonMID");
      event.put(fTGenLeptonMStatus, "GenLeptonMStatus");
      event.put(fTGenLeptonMPt, "GenLeptonMPt");
      event.put(fTGenLeptonMEta, "GenLeptonMEta");
      event.put(fTGenLeptonMPhi, "GenLeptonMPhi");
      event.put(fTGenLeptonGMID, "GenLeptonGMID");
      event.put(fTGenLeptonGMStatus, "GenLeptonGMStatus");
      event.put(fTGenLeptonGMPt, "GenLeptonGMPt");
      event.put(fTGenLeptonGMEta, "GenLeptonGMEta");
      event.put(fTGenLeptonGMPhi, "GenLeptonGMPhi");
      event.put(fTNGenPhotons, "NGenPhotons");
      event.put(fTGenPhotonPt, "GenPhotonPt");
      event.put(fTGenPhotonEta, "GenPhotonEta");
      event.put(fTGenPhotonPhi, "GenPhotonPhi");
      event.put(fTGenPhotonPartonMindR, "GenPhotonPartonMindR");
      event.put(fTGenPhotonMotherID, "GenPhotonMotherID");
      event.put(fTGenPhotonMotherStatus, "GenPhotonMotherStatus");
      event.put(fTNGenJets, "NGenJets");
      event.put(fTGenJetPt, "GenJetPt");
      event.put(fTGenJetEta, "GenJetEta");
      event.put(fTGenJetPhi, "GenJetPhi");
      event.put(fTGenJetE, "GenJetE");
      event.put(fTGenJetEmE, "GenJetEmE");
      event.put(fTGenJetHadE, "GenJetHadE");
      event.put(fTGenJetInvE, "GenJetInvE");
      event.put(fTNVrtx, "NVrtx");
      event.put(fTVrtxX, "VrtxX");
      event.put(fTVrtxY, "VrtxY");
      event.put(fTVrtxZ, "VrtxZ");
      event.put(fTVrtxXE, "VrtxXE");
      event.put(fTVrtxYE, "VrtxYE");
      event.put(fTVrtxZE, "VrtxZE");
      event.put(fTVrtxNdof, "VrtxNdof");
      event.put(fTVrtxChi2, "VrtxChi2");
      event.put(fTVrtxNtrks, "VrtxNtrks");
      event.put(fTVrtxSumPt, "VrtxSumPt");
      event.put(fTVrtxIsFake, "VrtxIsFake");
      event.put(fTNMus, "NMus");
      event.put(fTNMusTot, "NMusTot");
      event.put(fTNGMus, "NGMus");
      event.put(fTNTMus, "NTMus");
      event.put(fTMuGood, "MuGood");
      event.put(fTMuIsIso, "MuIsIso");
      event.put(fTMuIsGlobalMuon, "MuIsGlobalMuon");
      event.put(fTMuIsTrackerMuon, "MuIsTrackerMuon");
      event.put(fTMuPx, "MuPx");
      event.put(fTMuPy, "MuPy");
      event.put(fTMuPz, "MuPz");
      event.put(fTMuPt, "MuPt");
      event.put(fTMuInnerTkPt, "MuInnerTkPt");
      event.put(fTMuPtE, "MuPtE");
      event.put(fTMuE, "MuE");
      event.put(fTMuEt, "MuEt");
      event.put(fTMuEta, "MuEta");
      event.put(fTMuPhi, "MuPhi");
      event.put(fTMuCharge, "MuCharge");
      event.put(fTMuRelIso03, "MuRelIso03");
      event.put(fTMuIso03SumPt, "MuIso03SumPt");
      event.put(fTMuIso03EmEt, "MuIso03EmEt");
      event.put(fTMuIso03HadEt, "MuIso03HadEt");
      event.put(fTMuIso03EMVetoEt, "MuIso03EMVetoEt");
      event.put(fTMuIso03HadVetoEt, "MuIso03HadVetoEt");
      event.put(fTMuIso05SumPt, "MuIso05SumPt");
      event.put(fTMuIso05EmEt, "MuIso05EmEt");
      event.put(fTMuIso05HadEt, "MuIso05HadEt");
      event.put(fTMuEem, "MuEem");
      event.put(fTMuEhad, "MuEhad");
      event.put(fTMuD0BS, "MuD0BS");
      event.put(fTMuD0PV, "MuD0PV");
      event.put(fTMuD0E, "MuD0E");
      event.put(fTMuDzBS, "MuDzBS");
      event.put(fTMuDzPV, "MuDzPV");
      event.put(fTMuDzE, "MuDzE");
      event.put(fTMuNChi2, "MuNChi2");
      event.put(fTMuNGlHits, "MuNGlHits");
      event.put(fTMuNMuHits, "MuNMuHits");
      event.put(fTMuNTkHits, "MuNTkHits");
      event.put(fTMuNPxHits, "MuNPxHits");
      event.put(fTMuInnerTkNChi2, "MuInnerTkNChi2");
      event.put(fTMuNMatches, "MuNMatches");
      event.put(fTMuNChambers, "MuNChambers");
      event.put(fTMuCaloComp, "MuCaloComp");
      event.put(fTMuSegmComp, "MuSegmComp");
      event.put(fTMuIsGMPT, "MuIsGMPT");
      event.put(fTMuIsGMTkChiComp, "MuIsGMTkChiComp");
      event.put(fTMuIsGMStaChiComp, "MuIsGMStaChiComp");
      event.put(fTMuIsGMTkKinkTight, "MuIsGMTkKinkTight");
      event.put(fTMuIsAllStaMuons, "MuIsAllStaMuons");
      event.put(fTMuIsAllTrkMuons, "MuIsAllTrkMuons");
      event.put(fTMuIsTrkMuonArbitrated, "MuIsTrkMuonArbitrated");
      event.put(fTMuIsAllArbitrated, "MuIsAllArbitrated");
      event.put(fTMuIsTMLSLoose, "MuIsTMLSLoose");
      event.put(fTMuIsTMLSTight, "MuIsTMLSTight");
      event.put(fTMuIsTM2DCompLoose, "MuIsTM2DCompLoose");
      event.put(fTMuIsTM2DCompTight, "MuIsTM2DCompTight");
      event.put(fTMuIsTMOneStationLoose, "MuIsTMOneStationLoose");
      event.put(fTMuIsTMOneStationTight, "MuIsTMOneStationTight");
      event.put(fTMuIsTMLSOptLowPtLoose, "MuIsTMLSOptLowPtLoose");
      event.put(fTMuIsTMLSAngLoose, "MuIsTMLSAngLoose");
      event.put(fTMuIsTMLSAngTight, "MuIsTMLSAngTight");
      event.put(fTMuIsTMOneStationAngTight, "MuIsTMOneStationAngTight");
      event.put(fTMuIsTMOneStationAngLoose, "MuIsTMOneStationAngLoose");
      event.put(fTMuGenID, "MuGenID");
      event.put(fTMuGenStatus, "MuGenStatus");
      event.put(fTMuGenPt, "MuGenPt");
      event.put(fTMuGenEta, "MuGenEta");
      event.put(fTMuGenPhi, "MuGenPhi");
      event.put(fTMuGenE, "MuGenE");
      event.put(fTMuGenMID, "MuGenMID");
      event.put(fTMuGenMStatus, "MuGenMStatus");
      event.put(fTMuGenMPt, "MuGenMPt");
      event.put(fTMuGenMEta, "MuGenMEta");
      event.put(fTMuGenMPhi, "MuGenMPhi");
      event.put(fTMuGenME, "MuGenME");
      event.put(fTMuGenGMID, "MuGenGMID");
      event.put(fTMuGenGMStatus, "MuGenGMStatus");
      event.put(fTMuGenGMPt, "MuGenGMPt");
      event.put(fTMuGenGMEta, "MuGenGMEta");
      event.put(fTMuGenGMPhi, "MuGenGMPhi");
      event.put(fTMuGenGME, "MuGenGME");
      event.put(fTNEBhits, "NEBhits");
      event.put(fTEBrechitE, "EBrechitE");
      event.put(fTEBrechitPt, "EBrechitPt");
      event.put(fTEBrechitEta, "EBrechitEta");
      event.put(fTEBrechitPhi, "EBrechitPhi");
      event.put(fTEBrechitChi2, "EBrechitChi2");
      event.put(fTEBrechitTime, "EBrechitTime");
      event.put(fTEBrechitE4oE1, "EBrechitE4oE1");
      event.put(fTEBrechitE2oE9, "EBrechitE2oE9");
      event.put(fTNEles, "NEles");
      event.put(fTNElesTot, "NElesTot");
      event.put(fTElGood, "ElGood");
      event.put(fTElIsIso, "ElIsIso");
      event.put(fTElChargeMisIDProb, "ElChargeMisIDProb");
      event.put(fTElPx, "ElPx");
      event.put(fTElPy, "ElPy");
      event.put(fTElPz, "ElPz");
      event.put(fTElPt, "ElPt");
      event.put(fTElPtE, "ElPtE");
      event.put(fTElE, "ElE");
      event.put(fTElEt, "ElEt");
      event.put(fTElEta, "ElEta");
      event.put(fTElTheta, "ElTheta");
      event.put(fTElSCEta, "ElSCEta");
      event.put(fTElPhi, "ElPhi");
      event.put(fTElGsfTkPt, "ElGsfTkPt");
      event.put(fTElGsfTkEta, "ElGsfTkEta");
      event.put(fTElGsfTkPhi, "ElGsfTkPhi");
      event.put(fTElTrkMomentumError, "ElTrkMomentumError");
      event.put(fTElEcalEnergyError, "ElEcalEnergyError");
      event.put(fTElEleMomentumError, "ElEleMomentumError");
      event.put(fTElNBrems, "ElNBrems");
      event.put(fTElD0BS, "ElD0BS");
      event.put(fTElD0PV, "ElD0PV");
      event.put(fTElD0E, "ElD0E");
      event.put(fTElDzBS, "ElDzBS");
      event.put(fTElDzPV, "ElDzPV");
      event.put(fTElDzE, "ElDzE");
      event.put(fTElRelIso03, "ElRelIso03");
      event.put(fTElRelIso04, "ElRelIso04");
      event.put(fTElDR03TkSumPt, "ElDR03TkSumPt");
      event.put(fTElDR04TkSumPt, "ElDR04TkSumPt");
      event.put(fTElDR03EcalRecHitSumEt, "ElDR03EcalRecHitSumEt");
      event.put(fTElDR04EcalRecHitSumEt, "ElDR04EcalRecHitSumEt");
      event.put(fTElDR03HcalTowerSumEt, "ElDR03HcalTowerSumEt");
      event.put(fTElDR04HcalTowerSumEt, "ElDR04HcalTowerSumEt");
      event.put(fTElNChi2, "ElNChi2");
      event.put(fTElCharge, "ElCharge");
      event.put(fTElCInfoIsGsfCtfCons, "ElCInfoIsGsfCtfCons");
      event.put(fTElCInfoIsGsfCtfScPixCons, "ElCInfoIsGsfCtfScPixCons");
      event.put(fTElCInfoIsGsfScPixCons, "ElCInfoIsGsfScPixCons");
      event.put(fTElScPixCharge, "ElScPixCharge");
      event.put(fTElClosestCtfTrackPt, "ElClosestCtfTrackPt");
      event.put(fTElClosestCtfTrackEta, "ElClosestCtfTrackEta");
      event.put(fTElClosestCtfTrackPhi, "ElClosestCtfTrackPhi");
      event.put(fTElClosestCtfTrackCharge, "ElClosestCtfTrackCharge");
      event.put(fTElIDMva, "ElIDMva");
      event.put(fTElIDTight, "ElIDTight");
      event.put(fTElIDLoose, "ElIDLoose");
      event.put(fTElIDRobustTight, "ElIDRobustTight");
      event.put(fTElIDRobustLoose, "ElIDRobustLoose");
      event.put(fTElIDsimpleWPrelIso, "ElIDsimpleWPrelIso");
      event.put(fTElIDsimpleWP80relIso, "ElIDsimpleWP80relIso");
      event.put(fTElIDsimpleWP85relIso, "ElIDsimpleWP85relIso");
      event.put(fTElIDsimpleWP90relIso, "ElIDsimpleWP90relIso");
      event.put(fTElIDsimpleWP95relIso, "ElIDsimpleWP95relIso");
      event.put(fTElInGap, "ElInGap");
      event.put(fTElEcalDriven, "ElEcalDriven");
      event.put(fTElTrackerDriven, "ElTrackerDriven");
      event.put(fTElBasicClustersSize, "ElBasicClustersSize");
      event.put(fTElfbrem, "Elfbrem");
      event.put(fTElHcalOverEcal, "ElHcalOverEcal");
      event.put(fTElE1x5, "ElE1x5");
      event.put(fTElE5x5, "ElE5x5");
      event.put(fTElE2x5Max, "ElE2x5Max");
      event.put(fTElSigmaIetaIeta, "ElSigmaIetaIeta");
      event.put(fTElDeltaPhiSeedClusterAtCalo, "ElDeltaPhiSeedClusterAtCalo");
      event.put(fTElDeltaEtaSeedClusterAtCalo, "ElDeltaEtaSeedClusterAtCalo");
      event.put(fTElDeltaPhiSuperClusterAtVtx, "ElDeltaPhiSuperClusterAtVtx");
      event.put(fTElDeltaEtaSuperClusterAtVtx, "ElDeltaEtaSuperClusterAtVtx");
      event.put(fTElCaloEnergy, "ElCaloEnergy");
      event.put(fTElTrkMomAtVtx, "ElTrkMomAtVtx");
      event.put(fTElESuperClusterOverP, "ElESuperClusterOverP");
      event.put(fTElNumberOfMissingInnerHits, "ElNumberOfMissingInnerHits");
      event.put(fTElSCindex, "ElSCindex");
      event.put(fTElConvPartnerTrkDist, "ElConvPartnerTrkDist");
      event.put(fTElConvPartnerTrkDCot, "ElConvPartnerTrkDCot");
      event.put(fTElConvPartnerTrkPt, "ElConvPartnerTrkPt");
      event.put(fTElConvPartnerTrkEta, "ElConvPartnerTrkEta");
      event.put(fTElConvPartnerTrkPhi, "ElConvPartnerTrkPhi");
      event.put(fTElConvPartnerTrkCharge, "ElConvPartnerTrkCharge");
      event.put(fTElScSeedSeverity, "ElScSeedSeverity");
      event.put(fTElE1OverE9, "ElE1OverE9");
      event.put(fTElS4OverS1, "ElS4OverS1");
      event.put(fTElGenID, "ElGenID");
      event.put(fTElGenStatus, "ElGenStatus");
      event.put(fTElGenPt, "ElGenPt");
      event.put(fTElGenEta, "ElGenEta");
      event.put(fTElGenPhi, "ElGenPhi");
      event.put(fTElGenE, "ElGenE");
      event.put(fTElGenMID, "ElGenMID");
      event.put(fTElGenMStatus, "ElGenMStatus");
      event.put(fTElGenMPt, "ElGenMPt");
      event.put(fTElGenMEta, "ElGenMEta");
      event.put(fTElGenMPhi, "ElGenMPhi");
      event.put(fTElGenME, "ElGenME");
      event.put(fTElGenGMID, "ElGenGMID");
      event.put(fTElGenGMStatus, "ElGenGMStatus");
      event.put(fTElGenGMPt, "ElGenGMPt");
      event.put(fTElGenGMEta, "ElGenGMEta");
      event.put(fTElGenGMPhi, "ElGenGMPhi");
      event.put(fTElGenGME, "ElGenGME");
      event.put(fTNPhotons, "NPhotons");
      event.put(fTNPhotonsTot, "NPhotonsTot");
      event.put(fTPhoGood, "PhoGood");
      event.put(fTPhoIsIso, "PhoIsIso");
      event.put(fTPhoPt, "PhoPt");
      event.put(fTPhoPx, "PhoPx");
      event.put(fTPhoPy, "PhoPy");
      event.put(fTPhoPz, "PhoPz");
      event.put(fTPhoEta, "PhoEta");
      event.put(fTPhoPhi, "PhoPhi");
      event.put(fTPhoEnergy, "PhoEnergy");
      event.put(fTPhoIso03Ecal, "PhoIso03Ecal");
      event.put(fTPhoIso03Hcal, "PhoIso03Hcal");
      event.put(fTPhoIso03TrkSolid, "PhoIso03TrkSolid");
      event.put(fTPhoIso03TrkHollow, "PhoIso03TrkHollow");
      event.put(fTPhoIso03, "PhoIso03");
      event.put(fTPhoIso04Ecal, "PhoIso04Ecal");
      event.put(fTPhoIso04Hcal, "PhoIso04Hcal");
      event.put(fTPhoIso04TrkSolid, "PhoIso04TrkSolid");
      event.put(fTPhoIso04TrkHollow, "PhoIso04TrkHollow");
      event.put(fTPhoIso04, "PhoIso04");
      event.put(fTPhoR9, "PhoR9");
      event.put(fTPhoCaloPositionX, "PhoCaloPositionX");
      event.put(fTPhoCaloPositionY, "PhoCaloPositionY");
      event.put(fTPhoCaloPositionZ, "PhoCaloPositionZ");
      event.put(fTPhoHoverE, "PhoHoverE");
      event.put(fTPhoH1overE, "PhoH1overE");
      event.put(fTPhoH2overE, "PhoH2overE");
      event.put(fTPhoSigmaIetaIeta, "PhoSigmaIetaIeta");
      event.put(fTPhoSCRawEnergy, "PhoSCRawEnergy");
      event.put(fTPhoSCEtaWidth, "PhoSCEtaWidth");
      event.put(fTPhoSCSigmaPhiPhi, "PhoSCSigmaPhiPhi");
      event.put(fTPhoHasPixSeed, "PhoHasPixSeed");
      event.put(fTPhoHasConvTrks, "PhoHasConvTrks");
      event.put(fTPhoScSeedSeverity, "PhoScSeedSeverity");
      event.put(fTPhoE1OverE9, "PhoE1OverE9");
      event.put(fTPhoS4OverS1, "PhoS4OverS1");
      event.put(fTPhoSigmaEtaEta, "PhoSigmaEtaEta");
      event.put(fTPhoE1x5, "PhoE1x5");
      event.put(fTPhoE2x5, "PhoE2x5");
      event.put(fTPhoE3x3, "PhoE3x3");
      event.put(fTPhoE5x5, "PhoE5x5");
      event.put(fTPhomaxEnergyXtal, "PhomaxEnergyXtal");
      event.put(fTPhoIso03HcalDepth1, "PhoIso03HcalDepth1");
      event.put(fTPhoIso03HcalDepth2, "PhoIso03HcalDepth2");
      event.put(fTPhoIso04HcalDepth1, "PhoIso04HcalDepth1");
      event.put(fTPhoIso04HcalDepth2, "PhoIso04HcalDepth2");
      event.put(fTPhoIso03nTrksSolid, "PhoIso03nTrksSolid");
      event.put(fTPhoIso03nTrksHollow, "PhoIso03nTrksHollow");
      event.put(fTPhoIso04nTrksSolid, "PhoIso04nTrksSolid");
      event.put(fTPhoIso04nTrksHollow, "PhoIso04nTrksHollow");
      event.put(fTPhoisEB, "PhoisEB");
      event.put(fTPhoisEE, "PhoisEE");
      event.put(fTPhoisEBEtaGap, "PhoisEBEtaGap");
      event.put(fTPhoisEBPhiGap, "PhoisEBPhiGap");
      event.put(fTPhoisEERingGap, "PhoisEERingGap");
      event.put(fTPhoisEEDeeGap, "PhoisEEDeeGap");
      event.put(fTPhoisEBEEGap, "PhoisEBEEGap");
      event.put(fTPhoisPFlowPhoton, "PhoisPFlowPhoton");
      event.put(fTPhoisStandardPhoton, "PhoisStandardPhoton");
      event.put(fTPhoMCmatchindex, "PhoMCmatchindex");
      event.put(fTPhoMCmatchexitcode, "PhoMCmatchexitcode");
      event.put(fTPhoChargedHadronIso, "PhoChargedHadronIso");
      event.put(fTPhoNeutralHadronIso, "PhoNeutralHadronIso");
      event.put(fTPhoPhotonIso, "PhoPhotonIso");
      event.put(fTPhoisPFPhoton, "PhoisPFPhoton");
      event.put(fTPhoisPFElectron, "PhoisPFElectron");
      event.put(fTPhotSCindex, "PhotSCindex");
      event.put(fTPhoCone04PhotonIsodR0dEta0pt0, "PhoCone04PhotonIsodR0dEta0pt0");
      event.put(fTPhoCone04PhotonIsodR0dEta0pt5, "PhoCone04PhotonIsodR0dEta0pt5");
      event.put(fTPhoCone04PhotonIsodR8dEta0pt0, "PhoCone04PhotonIsodR8dEta0pt0");
      event.put(fTPhoCone04PhotonIsodR8dEta0pt5, "PhoCone04PhotonIsodR8dEta0pt5");
      event.put(fTPhoCone01PhotonIsodR045EB070EEdEta015pt08EB1EEmvVtx, "PhoCone01PhotonIsodR045EB070EEdEta015pt08EB1EEmvVtx");
      event.put(fTPhoCone02PhotonIsodR045EB070EEdEta015pt08EB1EEmvVtx, "PhoCone02PhotonIsodR045EB070EEdEta015pt08EB1EEmvVtx");
      event.put(fTPhoCone03PhotonIsodR045EB070EEdEta015pt08EB1EEmvVtx, "PhoCone03PhotonIsodR045EB070EEdEta015pt08EB1EEmvVtx");
      event.put(fTPhoCone04PhotonIsodR045EB070EEdEta015pt08EB1EEmvVtx, "PhoCone04PhotonIsodR045EB070EEdEta015pt08EB1EEmvVtx");
      event.put(fTPhoCone04NeutralHadronIsodR0dEta0pt0, "PhoCone04NeutralHadronIsodR0dEta0pt0");
      event.put(fTPhoCone04NeutralHadronIsodR0dEta0pt5, "PhoCone04NeutralHadronIsodR0dEta0pt5");
      event.put(fTPhoCone04NeutralHadronIsodR0dEta0pt0nocracks, "PhoCone04NeutralHadronIsodR0dEta0pt0nocracks");
      event.put(fTPhoCone04NeutralHadronIsodR0dEta0pt5nocracks, "PhoCone04NeutralHadronIsodR0dEta0pt5nocracks");
      event.put(fTPhoCone04NeutralHadronIsodR7dEta0pt0, "PhoCone04NeutralHadronIsodR7dEta0pt0");
      event.put(fTPhoCone04NeutralHadronIsodR7dEta0pt5, "PhoCone04NeutralHadronIsodR7dEta0pt5");
      event.put(fTPhoCone01NeutralHadronIsodR0dEta0pt0mvVtx, "PhoCone01NeutralHadronIsodR0dEta0pt0mvVtx");
      event.put(fTPhoCone02NeutralHadronIsodR0dEta0pt0mvVtx, "PhoCone02NeutralHadronIsodR0dEta0pt0mvVtx");
      event.put(fTPhoCone03NeutralHadronIsodR0dEta0pt0mvVtx, "PhoCone03NeutralHadronIsodR0dEta0pt0mvVtx");
      event.put(fTPhoCone04NeutralHadronIsodR0dEta0pt0mvVtx, "PhoCone04NeutralHadronIsodR0dEta0pt0mvVtx");
      event.put(fTPhoCone04ChargedHadronIsodR0dEta0pt0dz0old, "PhoCone04ChargedHadronIsodR0dEta0pt0dz0old");
      event.put(fTPhoCone04ChargedHadronIsodR0dEta0pt0PFnoPUold, "PhoCone04ChargedHadronIsodR0dEta0pt0PFnoPUold");
      event.put(fTPhoCone04ChargedHadronIsodR015dEta0pt0dz0old, "PhoCone04ChargedHadronIsodR015dEta0pt0dz0old");
      event.put(fTPhoCone04ChargedHadronIsodR015dEta0pt0PFnoPUold, "PhoCone04ChargedHadronIsodR015dEta0pt0PFnoPUold");
      event.put(fTPhoCone01ChargedHadronIsodR0dEta0pt0dz0, "PhoCone01ChargedHadronIsodR0dEta0pt0dz0");
      event.put(fTPhoCone01ChargedHadronIsodR0dEta0pt0dz1dxy01, "PhoCone01ChargedHadronIsodR0dEta0pt0dz1dxy01");
      event.put(fTPhoCone01ChargedHadronIsodR0dEta0pt0PFnoPU, "PhoCone01ChargedHadronIsodR0dEta0pt0PFnoPU");
      event.put(fTPhoCone01ChargedHadronIsodR015dEta0pt0dz0, "PhoCone01ChargedHadronIsodR015dEta0pt0dz0");
      event.put(fTPhoCone01ChargedHadronIsodR015dEta0pt0dz1dxy01, "PhoCone01ChargedHadronIsodR015dEta0pt0dz1dxy01");
      event.put(fTPhoCone01ChargedHadronIsodR015dEta0pt0PFnoPU, "PhoCone01ChargedHadronIsodR015dEta0pt0PFnoPU");
      event.put(fTPhoCone02ChargedHadronIsodR0dEta0pt0dz0, "PhoCone02ChargedHadronIsodR0dEta0pt0dz0");
      event.put(fTPhoCone02ChargedHadronIsodR0dEta0pt0dz1dxy01, "PhoCone02ChargedHadronIsodR0dEta0pt0dz1dxy01");
      event.put(fTPhoCone02ChargedHadronIsodR0dEta0pt0PFnoPU, "PhoCone02ChargedHadronIsodR0dEta0pt0PFnoPU");
      event.put(fTPhoCone02ChargedHadronIsodR015dEta0pt0dz0, "PhoCone02ChargedHadronIsodR015dEta0pt0dz0");
      event.put(fTPhoCone02ChargedHadronIsodR015dEta0pt0dz1dxy01, "PhoCone02ChargedHadronIsodR015dEta0pt0dz1dxy01");
      event.put(fTPhoCone02ChargedHadronIsodR015dEta0pt0PFnoPU, "PhoCone02ChargedHadronIsodR015dEta0pt0PFnoPU");
      event.put(fTPhoCone03ChargedHadronIsodR0dEta0pt0dz0, "PhoCone03ChargedHadronIsodR0dEta0pt0dz0");
      event.put(fTPhoCone03ChargedHadronIsodR0dEta0pt0dz1dxy01, "PhoCone03ChargedHadronIsodR0dEta0pt0dz1dxy01");
      event.put(fTPhoCone03ChargedHadronIsodR0dEta0pt0PFnoPU, "PhoCone03ChargedHadronIsodR0dEta0pt0PFnoPU");
      event.put(fTPhoCone03ChargedHadronIsodR015dEta0pt0dz0, "PhoCone03ChargedHadronIsodR015dEta0pt0dz0");
      event.put(fTPhoCone03ChargedHadronIsodR015dEta0pt0dz1dxy01, "PhoCone03ChargedHadronIsodR015dEta0pt0dz1dxy01");
      event.put(fTPhoCone03ChargedHadronIsodR015dEta0pt0PFnoPU, "PhoCone03ChargedHadronIsodR015dEta0pt0PFnoPU");
      event.put(fTPhoCone04ChargedHadronIsodR0dEta0pt0dz0, "PhoCone04ChargedHadronIsodR0dEta0pt0dz0");
      event.put(fTPhoCone04ChargedHadronIsodR0dEta0pt0dz1dxy01, "PhoCone04ChargedHadronIsodR0dEta0pt0dz1dxy01");
      event.put(fTPhoCone04ChargedHadronIsodR0dEta0pt0PFnoPU, "PhoCone04ChargedHadronIsodR0dEta0pt0PFnoPU");
      event.put(fTPhoCone04ChargedHadronIsodR015dEta0pt0dz0, "PhoCone04ChargedHadronIsodR015dEta0pt0dz0");
      event.put(fTPhoCone04ChargedHadronIsodR015dEta0pt0dz1dxy01, "PhoCone04ChargedHadronIsodR015dEta0pt0dz1dxy01");
      event.put(fTPhoCone04ChargedHadronIsodR015dEta0pt0PFnoPU, "PhoCone04ChargedHadronIsodR015dEta0pt0PFnoPU");
      event.put(fTPhoConvValidVtx, "PhoConvValidVtx");
      event.put(fTPhoConvNtracks, "PhoConvNtracks");
      event.put(fTPhoConvChi2Probability, "PhoConvChi2Probability");
      event.put(fTPhoConvEoverP, "PhoConvEoverP");
      event.put(fTNconv, "Nconv");
      event.put(fTConvValidVtx, "ConvValidVtx");
      event.put(fTConvNtracks, "ConvNtracks");
      event.put(fTConvChi2Probability, "ConvChi2Probability");
      event.put(fTConvEoverP, "ConvEoverP");
      event.put(fTConvZofPrimVtxFromTrks, "ConvZofPrimVtxFromTrks");
      //     event.put(fTNgv, "Ngv");
      //     event.put(fTgvSumPtHi, "gvSumPtHi");
      //     event.put(fTgvSumPtLo, "gvSumPtLo");
      //     event.put(fTgvNTkHi, "gvNTkHi");
      //     event.put(fTgvNTkLo, "gvNTkLo");
      event.put(fTNSuperClusters, "NSuperClusters");
      event.put(fTSCRaw, "SCRaw");
      event.put(fTSCPre, "SCPre");
      event.put(fTSCEnergy, "SCEnergy");
      event.put(fTSCEta, "SCEta");
      event.put(fTSCPhi, "SCPhi");
      event.put(fTSCPhiWidth, "SCPhiWidth");
      event.put(fTSCEtaWidth, "SCEtaWidth");
      event.put(fTSCBrem, "SCBrem");
      event.put(fTSCR9, "SCR9");
      event.put(fTSCcrackcorrseed, "SCcrackcorrseed");
      event.put(fTSCcrackcorr, "SCcrackcorr");
      event.put(fTSClocalcorrseed, "SClocalcorrseed");
      event.put(fTSClocalcorr, "SClocalcorr");
      event.put(fTSCcrackcorrseedfactor, "SCcrackcorrseedfactor");
      event.put(fTSClocalcorrseedfactor, "SClocalcorrseedfactor");
      event.put(fTNJets, "NJets");
      event.put(fTNJetsTot, "NJetsTot");
      event.put(fTJGood, "JGood");
      event.put(fTJPx, "JPx");
      event.put(fTJPy, "JPy");
      event.put(fTJPz, "JPz");
      event.put(fTJPt, "JPt");
      event.put(fTJE, "JE");
      event.put(fTJEt, "JEt");
      event.put(fTJEta, "JEta");
      event.put(fTJPhi, "JPhi");
      event.put(fTJEcorr, "JEcorr");
      event.put(fTJArea, "JArea");
      event.put(fTJEtaRms, "JEtaRms");
      event.put(fTJPhiRms, "JPhiRms");
      event.put(fTJNConstituents, "JNConstituents");
      event.put(fTJNAssoTracks, "JNAssoTracks");
      event.put(fTJNNeutrals, "JNNeutrals");
      event.put(fTJChargedEmFrac, "JChargedEmFrac");
      event.put(fTJNeutralEmFrac, "JNeutralEmFrac");
      event.put(fTJChargedHadFrac, "JChargedHadFrac");
      event.put(fTJNeutralHadFrac, "JNeutralHadFrac");
      event.put(fTJChargedMuEnergyFrac, "JChargedMuEnergyFrac");
      event.put(fTJPhoFrac, "JPhoFrac");
      event.put(fTJHFHadFrac, "JHFHadFrac");
      event.put(fTJHFEMFrac, "JHFEMFrac");
      event.put(fTJPtD, "JPtD");
      event.put(fTJRMSCand, "JRMSCand");
      event.put(fTJeMinDR, "JeMinDR");
      event.put(fTJbTagProbTkCntHighEff, "JbTagProbTkCntHighEff");
      event.put(fTJbTagProbTkCntHighPur, "JbTagProbTkCntHighPur");
      event.put(fTJbTagProbSimpSVHighEff, "JbTagProbSimpSVHighEff");
      event.put(fTJbTagProbSimpSVHighPur, "JbTagProbSimpSVHighPur");
      event.put(fTJMass, "JMass");
      event.put(fTJtrk1px, "Jtrk1px");
      event.put(fTJtrk1py, "Jtrk1py");
      event.put(fTJtrk1pz, "Jtrk1pz");
      event.put(fTJtrk2px, "Jtrk2px");
      event.put(fTJtrk2py, "Jtrk2py");
      event.put(fTJtrk2pz, "Jtrk2pz");
      event.put(fTJtrk3px, "Jtrk3px");
      event.put(fTJtrk3py, "Jtrk3py");
      event.put(fTJtrk3pz, "Jtrk3pz");
      event.put(fTJVtxx, "JVtxx");
      event.put(fTJVtxy, "JVtxy");
      event.put(fTJVtxz, "JVtxz");
      event.put(fTJVtxExx, "JVtxExx");
      event.put(fTJVtxEyx, "JVtxEyx");
      event.put(fTJVtxEyy, "JVtxEyy");
      event.put(fTJVtxEzy, "JVtxEzy");
      event.put(fTJVtxEzz, "JVtxEzz");
      event.put(fTJVtxEzx, "JVtxEzx");
      event.put(fTJVtxNChi2, "JVtxNChi2");
      event.put(fTJGenJetIndex, "JGenJetIndex");
      event.put(fTNTracks, "NTracks");
      event.put(fTNTracksTot, "NTracksTot");
      event.put(fTTrkGood, "TrkGood");
      event.put(fTTrkPt, "TrkPt");
      event.put(fTTrkEta, "TrkEta");
      event.put(fTTrkPhi, "TrkPhi");
      event.put(fTTrkNChi2, "TrkNChi2");
      event.put(fTTrkNHits, "TrkNHits");
      event.put(fTTrkVtxDz, "TrkVtxDz");
      event.put(fTTrkVtxDxy, "TrkVtxDxy");
      event.put(fTTrkPtSumx, "TrkPtSumx");
      event.put(fTTrkPtSumy, "TrkPtSumy");
      event.put(fTTrkPtSum, "TrkPtSum");
      event.put(fTTrkPtSumPhi, "TrkPtSumPhi");
      event.put(fTSumEt, "SumEt");
      event.put(fTECALSumEt, "ECALSumEt");
      event.put(fTHCALSumEt, "HCALSumEt");
      event.put(fTECALEsumx, "ECALEsumx");
      event.put(fTECALEsumy, "ECALEsumy");
      event.put(fTECALEsumz, "ECALEsumz");
      event.put(fTECALMET, "ECALMET");
      event.put(fTECALMETPhi, "ECALMETPhi");
      event.put(fTECALMETEta, "ECALMETEta");
      event.put(fTHCALEsumx, "HCALEsumx");
      event.put(fTHCALEsumy, "HCALEsumy");
      event.put(fTHCALEsumz, "HCALEsumz");
      event.put(fTHCALMET, "HCALMET");
      event.put(fTHCALMETPhi, "HCALMETPhi");
      event.put(fTHCALMETeta, "HCALMETeta");
      event.put(fTRawMET, "RawMET");
      event.put(fTRawMETpx, "RawMETpx");
      event.put(fTRawMETpy, "RawMETpy");
      event.put(fTRawMETphi, "RawMETphi");
      event.put(fTRawMETemEtFrac, "RawMETemEtFrac");
      event.put(fTRawMETemEtInEB, "RawMETemEtInEB");
      event.put(fTRawMETemEtInEE, "RawMETemEtInEE");
      event.put(fTRawMETemEtInHF, "RawMETemEtInHF");
      event.put(fTRawMEThadEtFrac, "RawMEThadEtFrac");
      event.put(fTRawMEThadEtInHB, "RawMEThadEtInHB");
      event.put(fTRawMEThadEtInHE, "RawMEThadEtInHE");
      event.put(fTRawMEThadEtInHF, "RawMEThadEtInHF");
      event.put(fTRawMETSignificance, "RawMETSignificance");
      event.put(fTGenMET, "GenMET");
      event.put(fTGenMETpx, "GenMETpx");
      event.put(fTGenMETpy, "GenMETpy");
      event.put(fTGenMETphi, "GenMETphi");
      event.put(fTTCMET, "TCMET");
      event.put(fTTCMETpx, "TCMETpx");
      event.put(fTTCMETpy, "TCMETpy");
      event.put(fTTCMETphi, "TCMETphi");
      event.put(fTTCMETSignificance, "TCMETSignificance");
      event.put(fTMuJESCorrMET, "MuJESCorrMET");
      event.put(fTMuJESCorrMETpx, "MuJESCorrMETpx");
      event.put(fTMuJESCorrMETpy, "MuJESCorrMETpy");
      event.put(fTMuJESCorrMETphi, "MuJESCorrMETphi");
      event.put(fTPFMET, "PFMET");
      event.put(fTPFMETpx, "PFMETpx");
      event.put(fTPFMETpy, "PFMETpy");
      event.put(fTPFMETphi, "PFMETphi");
      event.put(fTPFMETSignificance, "PFMETSignificance");
      event.put(fTPFSumEt, "PFSumEt");
      event.put(fTPFMETPAT, "PFMETPAT");
      event.put(fTPFMETPATpx, "PFMETPATpx");
      event.put(fTPFMETPATpy, "PFMETPATpy");
      event.put(fTPFMETPATphi, "PFMETPATphi");
      event.put(fTPFMETPATSignificance, "PFMETPATSignificance");
      event.put(fTMETR12, "METR12");
      event.put(fTMETR21, "METR21");


    }

    // Method called once before each run
    bool NTupleProducer::beginRun(edm::Run& r, const edm::EventSetup& es){

      resetRunProducts();

      // Retrieve and fill RunTree information

      // Retrieve HLT trigger menu 
      bool changed;
      if ( !(fHltConfig.init(r,es,"*",changed)) ) { // "*": auto-discovery
        edm::LogError("NTP") << " hlt config extraction failure with process name "
                             << fHltConfig.processName();
      }
      fRHLTNames.reset(new std::vector<std::string>(fHltConfig.triggerNames()));

      // Store L1 trigger names
      edm::ESHandle<L1GtTriggerMenu> menuRcd;
      es.get<L1GtTriggerMenuRcd>().get(menuRcd);
      const L1GtTriggerMenu *menu = menuRcd.product();
      const AlgorithmMap& algoMap = menu->gtAlgorithmMap();
      (*fRL1PhysMenu).resize( gMaxL1PhysBits );
      for( AlgorithmMap::const_iterator it = algoMap.begin(); it != algoMap.end(); ++it ){
        (*fRL1PhysMenu)[(*it).second.algoBitNumber()] = std::string((*it).first); 
      }

      if(!fIsRealData) {
        edm::Handle<GenRunInfoProduct> genRunInfo;
        r.getByLabel("generator", genRunInfo);
        *fRExtXSecLO  = genRunInfo->externalXSecLO().value();
        *fRIntXSec    = genRunInfo->internalXSec().value();
        *fRExtXSecNLO = genRunInfo->externalXSecNLO().value();
      }
  
      *fRMinMuPt       = fMinMuPt;
      *fRMaxMuEta      = fMaxMuEta;
      *fRMinElPt       = fMinElPt;
      *fRMaxElEta      = fMaxElEta;
      *fRMinJPt        = fMinCorJPt;
      *fRMinRawJPt     = fMinRawJPt;
      *fRMaxJEta       = fMaxJEta;
      *fRMinJEMFrac    = fMinJEMFrac;

      *fRMinTrkPt      = fMinTrkPt;
      *fRMaxTrkEta     = fMaxTrkEta;
      *fRMaxTrkNChi2   = fMaxTrkNChi2;
      *fRMinTrkNHits   = fMinTrkNHits;

      *fRMinPhotonPt   = fMinPhotonPt;
      *fRMaxPhotonEta  = fMaxPhotonEta;
      *fRMinSCraw      = fMinSCraw;
      *fRMinEBRechitE  = fMinEBRechitE;
  
      *fRMinGenLeptPt  = fMinGenLeptPt;
      *fRMaxGenLeptEta = fMaxGenLeptEta;
      *fRMinGenPhotPt  = fMinGenPhotPt;
      *fRMaxGenPhotEta = fMaxGenPhotEta;
      *fRMinGenJetPt   = fMinGenJetPt;
      *fRMaxGenJetEta  = fMaxGenJetEta;

      *fRMaxNMus       = gMaxNMus;
      *fRMaxNEles      = gMaxNEles;
      *fRMaxNJets      = gMaxNJets;
      *fRMaxNTrks      = gMaxNTrks;
      *fRMaxNPhotons   = gMaxNPhotons;
      *fRMaxNSC        = gMaxNSC;
      *fRMaxNGenLept   = gMaxNGenLept;
      *fRMaxNGenPhot   = gMaxNGenPhot;
      *fRMaxNGenJets   = gMaxNGenJets;
      *fRMaxNVrtx      = gMaxNVrtx;
      *fRMaxNPileup    = gMaxNPileup;
      *fRMaxNEBhits    = gMaxNEBhits;

      return true; // Not an actual filter
    }

    // Method called once after each run
    bool NTupleProducer::endRun(edm::Run& r, const edm::EventSetup&){

      // Store run information
      r.put(fRExtXSecLO     ,"ExtXSecLO"     );
      r.put(fRExtXSecNLO    ,"ExtXSecNLO"    );
      r.put(fRIntXSec       ,"IntXSec"       );
                         
      r.put(fRMinMuPt       ,"MinMuPt"       );
      r.put(fRMaxMuEta      ,"MaxMuEta"      );
      r.put(fRMinElPt       ,"MinElPt"       );
      r.put(fRMaxElEta      ,"MaxElEta"      );
      r.put(fRMinJPt        ,"MinJPt"        );
      r.put(fRMinRawJPt     ,"MinRawJPt"     );
      r.put(fRMaxJEta       ,"MaxJEta"       );
      r.put(fRMinJEMFrac    ,"MinJEMfrac"    );
                                           
      r.put(fRMinTrkPt      ,"MinTrkPt"      );
      r.put(fRMaxTrkEta     ,"MaxTrkEta"     );
      r.put(fRMaxTrkNChi2   ,"MaxTrkNChi2"   );
      r.put(fRMinTrkNHits   ,"MinTrkNHits"   );
                                           
      r.put(fRMinPhotonPt   ,"MinPhotonPt"   );
      r.put(fRMaxPhotonEta  ,"MaxPhotonEta"  );
      r.put(fRMinSCraw      ,"MinSCraw"      );
      r.put(fRMinEBRechitE  ,"MinEBRechitE"  );
                                                                        
      r.put(fRMinGenLeptPt  ,"MinGenLeptPt"  );
      r.put(fRMaxGenLeptEta ,"MaxGenLeptEta" );
      r.put(fRMinGenPhotPt  ,"MinGenPhotPt"  );
      r.put(fRMaxGenPhotEta ,"MaxGenPhotEta" );
      r.put(fRMinGenJetPt   ,"MinGenJetPt"   );
      r.put(fRMaxGenJetEta  ,"MaxGenJetEta"  );
                                                                        
      r.put(fRMaxNMus    ,"MaxNMus"       );
      r.put(fRMaxNEles   ,"MaxNEles"      );
      r.put(fRMaxNJets   ,"MaxNJets"      );
      r.put(fRMaxNTrks   ,"MaxNTrks"      );
      r.put(fRMaxNPhotons,"MaxNPhotons"   );
      r.put(fRMaxNSC     ,"MaxNSC"        );
      r.put(fRMaxNGenLept,"MaxNGenLep"    );
      r.put(fRMaxNGenPhot,"MaxNGenPho"    );
      r.put(fRMaxNGenJets,"MaxNGenJet"    );
      r.put(fRMaxNVrtx  , "MaxNVrtx"      );
      r.put(fRMaxNPileup, "MaxNPileup"    );
      r.put(fRMaxNEBhits, "MaxNEBhits"    );

      r.put(fRHLTNames   ,"HLTNames"      );
      r.put(fRL1PhysMenu ,"L1PhysMenu"    );
      r.put(fRHLTLabels  ,"HLTLabels"     );

      r.put(fRPileUpData ,"PileUpData"    );
      r.put(fRPileUpMC   ,"PileUpMC"      );

      return true;
    }

    //________________________________________________________________________________________
    // Method called once each job just after ending the event loop
    void NTupleProducer::endJob( void ){

      edm::LogVerbatim("NTP") << " ---------------------------------------------------";
      edm::LogVerbatim("NTP") << " ==> NTupleProducer::endJob() ...";
      edm::LogVerbatim("NTP") << "  Total number of processed Events: " << fNTotEvents;
      edm::LogVerbatim("NTP") << "  Number of times Tree was filled:  " << fNFillTree;
      edm::LogVerbatim("NTP") << " ---------------------------------------------------";

    }

    //________________________________________________________________________________________
    // Method for matching of reco candidates
    std::vector<const reco::GenParticle*> NTupleProducer::matchRecoCand(const reco::RecoCandidate *Cand, const edm::Event& iEvent){
      using namespace std;
      const GenParticle *GenCand = NULL;
      const GenParticle *GenMom  = NULL;
      const GenParticle *GenGMom = NULL;
      std::vector<const GenParticle*> res;
      if(fIsRealData){
        edm::LogWarning("NTP") << "@SUB=matchRecoCand"
                               << "Trying to access generator info on real data...";
        res.push_back(GenCand);
        res.push_back(GenMom);
        res.push_back(GenGMom);
        return res;
      }

      int id(0), mid(0), gmid(0);

      edm::Handle<GenParticleCollection> genparts;
      iEvent.getByLabel(fGenPartTag, genparts);
      GenParticleCollection::const_iterator gpart;
      GenCand = new GenParticle();
      GenMom  = new GenParticle();
      GenGMom = new GenParticle();
      // initialize the results
      // btw, this is a little memory leak
      res.push_back(GenCand);
      res.push_back(GenMom);
      res.push_back(GenGMom);

      bool matched = false;

      // Try to match the reco candidate to a generator object
      double mindr(999.99);
      for(gpart = genparts->begin(); gpart != genparts->end(); gpart++){
        if( gpart->status() != 1 ) continue;

        // Restrict to cone of 0.1 in DR around candidate
        double dr = reco::deltaR(gpart->eta(), gpart->phi(), Cand->eta(), Cand->phi());
        if(dr > 0.1) continue;

        // Restrict to pt match within a factor of 2
        double ndpt = fabs(gpart->pt() - Cand->pt())/gpart->pt();
        if(ndpt > 2.) continue;

        // Minimize DeltaR
        if(dr > mindr) continue;
        mindr = dr;

        matched = true;
        GenCand = &(*gpart);
      }


      // Fill generator information in case match was successful
      if(matched){
        // update the results
        res[0]=GenCand;
        id  = GenCand->pdgId();
        if ( !GenCand->mother() ) return res;

        // Determine mother object of matched gen object:
        // (Make sure that the mother has a different PDG ID)
        GenMom = static_cast<const GenParticle*>(GenCand->mother());
        mid = GenMom->pdgId();
        if(mid != id); // do nothing
        else{
          int tempid = mid;
          while(tempid == id){
            GenMom = static_cast<const GenParticle*>(GenMom->mother());
            tempid = GenMom->pdgId();
          }
        }
        mid = GenMom->pdgId();
        // update the results
        res[1] = GenMom;

        if ( !GenMom->mother() ) return res;
        // Determine grand-mother object of matched gen object:
        // (Make sure that the grand-mother has a different PDG ID than the mother)
        GenGMom = static_cast<const GenParticle*>(GenMom->mother());
        gmid = GenGMom->pdgId();
        if(gmid != mid); // do nothing
        else{
          int tempid = gmid;
          while(tempid == mid){
            GenGMom = static_cast<const GenParticle*>(GenGMom->mother());
            tempid = GenGMom->pdgId();
          }
        }

        res[2] = GenGMom;
      }
      return res;
    }

    //________________________________________________________________________________________
    // Method for matching of jets
    const int NTupleProducer::matchJet(const reco::Jet* jet){
      // match to already filled genjets in ntuple:
      // this can only be called AFTER the genjets have been stored!
      if(fIsRealData || *fTNGenJets < 1){
        edm::LogWarning("NTP") << "@SUB=matchJet"
                               << "Trying to access generator info on real data...";
        return -1;
      }

      // Try to match the reco jet to a stored generator jet
      double mindr(999.99);
      int matchedindex = -1;
      for(int i = 0; i < *fTNGenJets; i++){

        // Restrict to cone of 0.1 in DR around candidate
        double dr = reco::deltaR((*fTGenJetEta)[i], (*fTGenJetPhi)[i], jet->eta(), jet->phi());
        if(dr > 0.3) continue;

        // Restrict to pt match within a factor of 2
        double ndpt = fabs((*fTGenJetPt)[i] - jet->pt())/(*fTGenJetPt)[i];
        if(ndpt > 2.) continue;

        // Minimize DeltaR
        if(dr > mindr) continue;
        mindr = dr;
        matchedindex = i;
      }
      return matchedindex;
    }

    /////////////////////////////////////////////////////////////////////////////////////////////
    // Cleaning methods
    //________________________________________________________________________________________
    void NTupleProducer::ElectronDuplicate(std::vector<const SuperCluster*> elecPtr, std::vector<const GsfTrack*> trckPtr) {
      // Looks for duplication among electrons
      if( *fTNEles <= 0 ) return;

      // loop over the electrons
      for (int i = 0; i < *fTNEles; ++i) {
        const SuperCluster* supercluster = elecPtr[i];
        const GsfTrack* eletrack = trckPtr[i];

        // loop over the electrons again
        for (int j = i+1; j < *fTNEles; ++j) {
          const SuperCluster* newsuper = elecPtr[j];
          const GsfTrack* newtrack = trckPtr[j];

          // check if duplicate
          if (newsuper == supercluster || newtrack == eletrack){
            // fTeDupEl[i] = j;
            // fTeDupEl[j] = i;
            break;
          }
        }
      }

      return;
    }

    //________________________________________________________________________________________
    void NTupleProducer::PhotonElectronDuplicate(std::vector<const SuperCluster*> elecPtr, std::vector<const SuperCluster*> phoPtr) {
      // Looks for duplication between photons and electrons
      if( *fTNEles <= 0 ) return;
      if( *fTNPhotons <= 0 ) return;

      // loop over the photons
      for( int i = 0; i < *fTNPhotons; ++i ){
        const SuperCluster* phoSC = phoPtr[i];

        // loop over the electrons again
        for( int j = 0; j < *fTNEles; ++j ){
          const SuperCluster* elSC = elecPtr[j];

          // check if duplicate
          if( elSC == phoSC ){
            // fTPhotDupEl[i] = j;
            break;
          }
        }
      }

      return;
    }

    //________________________________________________________________________________________
    void NTupleProducer::ElJetOverlap(std::vector<const Jet*> jets, std::vector<const SuperCluster*> electrons, edm::Handle<CaloTowerCollection> calotowers){
      // checks for jets made from electrons
      // jetIndex and elecIndex contain the indices of the selected jets and
      //   electrons in the Views
      // (electrons and jets should be filled in the ntuple before checking)
      if (*fTNJets <= 0) return;
      if (*fTNEles <= 0) return;

      std::vector<CaloTowerPtr> jetCaloRefs;

      // loop over the jets
      for (int i = 0; i < *fTNJets; ++i) {

        // Collect the CaloTowers detIds for the jet
        const CaloJet* theJet = static_cast<const CaloJet*>(&(*jets[i]));
        jetCaloRefs = theJet->getCaloConstituents();

        // loop over the electrons
        for( int j = 0; j < *fTNEles; ++j ){
          const SuperCluster* theElecSC = electrons[j];

          math::XYZVector sharedP(0., 0., 0.);
          bool isInJet = IsEMObjectInJet(theElecSC, jetCaloRefs, calotowers, &sharedP);
          // float sharedE = sqrt(sharedP.X()*sharedP.X() + sharedP.Y()*sharedP.Y() + sharedP.Z()*sharedP.Z() );
          if (isInJet) {
            // fTeIsInJet[j] = i;
            // fTeSharedPx[j] = sharedP.X();
            // fTeSharedPy[j] = sharedP.Y();
            // fTeSharedPz[j] = sharedP.Z();
            // fTeSharedEnergy[j] = sharedE;
            break;
          }
        }
      }

      return;
    }

    //________________________________________________________________________________________
    void NTupleProducer::PhotonJetOverlap(std::vector<const Jet*> jets, std::vector<const SuperCluster*> superclusters, edm::Handle<CaloTowerCollection> calotowers){
      // checks for jets made from photons
      // (photons and jets should be filled in the ntuple before checking)
      if( *fTNJets <= 0 ) return;
      if( *fTNPhotons <= 0 ) return;

      std::vector<CaloTowerPtr> jetCaloRefs;

      // loop over the jets
      for( int i = 0; i < *fTNJets; ++i ){

        // Collect the CaloTowers detIds for the jet
        const CaloJet* theJet = static_cast<const CaloJet*>(&(*jets[i]));
        jetCaloRefs = theJet->getCaloConstituents();

        // loop over the photons
        for( int j = 0; j < *fTNPhotons; ++j ){
          const SuperCluster* theSC = superclusters[j];

          math::XYZVector sharedP(0., 0., 0.);
          bool isInJet = IsEMObjectInJet(theSC, jetCaloRefs, calotowers, &sharedP);
          // float sharedE = sqrt(sharedP.X()*sharedP.X() + sharedP.Y()*sharedP.Y() + sharedP.Z()*sharedP.Z() );
          if( isInJet ){
            // ftPhoisInJet[j] = i;
            // fTPhotSharedPx[j] = sharedP.X();
            // fTPhotSharedPy[j] = sharedP.Y();
            // fTPhotSharedPz[j] = sharedP.Z();
            // fTPhotSharedEnergy[j] = sharedE;
            break;
          }
        }
      }

      return;
    }

    //________________________________________________________________________________________
    bool NTupleProducer::IsEMObjectInJet(const SuperCluster* elecSC, std::vector<CaloTowerPtr> jetCaloRefs, edm::Handle<CaloTowerCollection> calotowers, math::XYZVector* sharedMomentum){
      // Checks whether an electron or photon is included in the jet energy
      // and if true, it returns the momentum vector shared by the two

      // Define a window in eta,phi for the SuperCluster
      float phimin=0., phimax=0., etamin=0., etamax=0.;
      bool window = EMCaloTowerWindow(elecSC, phimin, phimax, etamin, etamax);
      if (!window){ return false;}

      // Collect the CaloTowers inside this window, save their detId in a vector
      std::vector<CaloTowerDetId> eleDetId;
      std::vector<float> eleTowerEnergy;
      std::vector<float> eleTowerEta;
      std::vector<float> eleTowerPhi;
      float eleEnergySum = 0.;
      CaloTowerCollection::const_iterator calo;
      for (calo = calotowers->begin(); calo != calotowers->end(); ++calo ){
        float tow_eta = calo->eta();
        float tow_phi = calo->phi();
        if (IsInPhiWindow (tow_phi, phimin, phimax)
            && (tow_eta-etamin)*(tow_eta-etamax) <= 0.){
          eleDetId.push_back(calo->id());
          eleTowerEnergy.push_back(calo->emEnergy());
          eleTowerEta.push_back(calo->eta());
          eleTowerPhi.push_back(calo->phi());
          eleEnergySum += calo->emEnergy();
        }
      }

      // Loop over the detIds vector of the electron
      float sharedEnergy = 0.;
      float sharedPx = 0.;
      float sharedPy = 0.;
      float sharedPz = 0.;
      for (unsigned int i = 0; i < eleDetId.size(); ++i){
        // find whether this detId is in the jet detId list
        for (unsigned int j = 0; j < jetCaloRefs.size(); ++j){
          // if yes, add its energy to the sum
          if (eleDetId[i] == jetCaloRefs[j]->id() ){
            sharedEnergy += eleTowerEnergy[i];
            float eleTowerTheta = 2. * atan(exp(-eleTowerEta[i]));
            if (eleTowerTheta < 0.) {eleTowerTheta += 3.141592654;}
            float sintheta = sin(eleTowerTheta);
            sharedPx += eleTowerEnergy[i]*sintheta*cos(eleTowerPhi[i]);
            sharedPy += eleTowerEnergy[i]*sintheta*sin(eleTowerPhi[i]);
            sharedPz += eleTowerEnergy[i]*cos(eleTowerTheta);
          }
        }
      }

      eleDetId.clear();
      eleTowerEnergy.clear();
      eleTowerEta.clear();
      eleTowerPhi.clear();
      jetCaloRefs.clear();

      sharedMomentum->SetXYZ(sharedPx,sharedPy,sharedPz);

      if (sharedEnergy > 0.) {return true;}
      else {return false;}

    }

    //________________________________________________________________________________________
    bool NTupleProducer::EMCaloTowerWindow(const SuperCluster* superCluster, float & phimin, float & phimax, float & etamin, float & etamax){
      // Define a window for CaloTowers around an electron or photon
      // First, find the extremes from the basicCluster positions
      phimin=0.;
      phimax=0.;
      etamin=0.;
      etamax=0.;
      float pi    = 3.141592654;
      float twopi = 6.283185307;
      float clusterEsum = 0.;
      reco::CaloCluster_iterator iCluster = superCluster->clustersBegin();
      for (; iCluster != superCluster->clustersEnd(); ++iCluster){
        math::XYZPoint clusterXYZ = (*iCluster)->position();
        float clusterE = (*iCluster)->energy();
        clusterEsum += clusterE;
        float clusterphi = atan2(clusterXYZ.Y(), clusterXYZ.X());
        float clustertheta = acos(clusterXYZ.Z() /
                                  sqrt(clusterXYZ.X()*clusterXYZ.X()+clusterXYZ.Y()*clusterXYZ.Y()
                                       +clusterXYZ.Z()*clusterXYZ.Z()) );
        if (clustertheta < 0.) {clustertheta = clustertheta + 3.141592654;}
        float clustereta = -log(tan(0.5*clustertheta));
        if (iCluster == superCluster->clustersBegin() ){
          etamin = clustereta;
          etamax = clustereta;
          phimin = clusterphi;
          phimax = clusterphi;
        } else {
          if (etamin > clustereta){etamin = clustereta;}
          if (etamax < clustereta){etamax = clustereta;}
          phimin = GetPhiMin(phimin, clusterphi);
          phimax = GetPhiMax(phimax, clusterphi);
        }
      }

      // Then put a tolerance for the bigger caloTowers
      // (adding 1/2 of the caloTower size)
      float etamean = fabs(0.5*(etamin+etamax) );
      phimin -= CaloTowerSizePhi(etamean);
      phimax += CaloTowerSizePhi(etamean);
      etamin -= CaloTowerSizeEta(etamin);
      etamax += CaloTowerSizeEta(etamax);
      // Check that they are still correctly normalized (-pi < phi < pi)
      if (phimin < -pi) {phimin += twopi;}
      if (phimax >  pi) {phimax -= twopi;}

      return true;

    }

    float NTupleProducer::CaloTowerSizePhi(float eta){
      // Returns the half size of a CaloTower in phi
      // numbers from ptdr1, p.201

      float sizePhi = 0.;
      if (fabs(eta) <= 1.74){
        sizePhi = 0.0435+0.0174;
      } else {
        sizePhi = 0.087+0.0174;
      }

      return sizePhi;
    }

    float NTupleProducer::CaloTowerSizeEta(float eta){
      // Returns the half size of a CaloTower in eta
      // numbers from ptdr1, p.201

      float abseta = fabs(eta);
      float sizeEta = 0.;
      if (abseta <= 1.74){
        sizeEta = 0.0435+0.0174;
      } else if (abseta <= 2.5){
        sizeEta = 0.0435 + 0.0678*(abseta-1.74)+0.0174;
      } else {
        sizeEta = 0.0875+0.0174;
      }

      return sizeEta;
    }

    bool NTupleProducer::IsInPhiWindow(float phi, float phimin, float phimax){
      // Checks whether phi is inside a given window

      float dphimin = reco::deltaPhi(phi, phimin);
      float dphimax = reco::deltaPhi(phi, phimax);

      float pi    = 3.141592654;
      if (dphimin*dphimax <= 0. && fabs(dphimin-dphimax) < pi){return true;}
      else {return false;}

    }

    float NTupleProducer::GetPhiMin(float phi1, float phi2){
      // Computes the minimum of two phi values

      float pi    = 3.141592654;
      float phimin = phi1;
      if ((phimin-phi2)>0. && (phimin-phi2)< pi){phimin = phi2;}
      else if ((phimin-phi2)<-pi){phimin = phi2;}

      return phimin;

    }

    float NTupleProducer::GetPhiMax(float phi1, float phi2){
      // Computes the minimum of two phi values

      float pi    = 3.141592654;
      float phimax = phi1;
      if ((phimax-phi2)<0. && (phimax-phi2)>-pi){phimax = phi2;}
      else if ((phimax-phi2)> pi){phimax = phi2;}

      return phimax;

    }

    //________________________________________________________________________________________
    double NTupleProducer::DeltaPhi(double phi1, double phi2){

      double dphi=phi1-phi2;
      if (dphi>TMath::Pi()) dphi=2*TMath::Pi()-dphi;
      if (dphi<-TMath::Pi()) dphi=-2*TMath::Pi()-dphi;

      return dphi;
    }


    //________________________________________________________________________________________
    // void NTupleProducer::FillPhotonIsoVariables(double photonEta, double photonPhi, double photonVz, int type, bool isPU, edm::Handle<reco::PFCandidateCollection>& pfCandidates, int ipf, int phoqi){

    //   double pt = (*pfCandidates)[ipf].pt();
    //   double dEta = fabs(photonEta - (*pfCandidates)[ipf].eta());
    //   double dPhi = DeltaPhi(photonPhi,(*pfCandidates)[ipf].phi());
    //   double dR = sqrt(dEta*dEta+dPhi*dPhi);
    //   //  double dz = fabs(photonVz - (*pfCandidates)[ipf].vz());


    //   //cout << "FillPhotonIsoVariables pt="<<pt<<endl;
    //   if (type==0){

    //     fTPhoCone04NeutralHadronIsodR0dEta0pt0[phoqi] += pt;
    //     if (pt>0.5) fTPhoCone04NeutralHadronIsodR0dEta0pt5[phoqi] += pt;
  
    //     if (dR>0.07) {
    //       fTPhoCone04NeutralHadronIsodR7dEta0pt0[phoqi] += pt;
    //       if (pt>0.5) fTPhoCone04NeutralHadronIsodR7dEta0pt5[phoqi] += pt;
    //     }
    //     if (isInEtaCracks((*pfCandidates)[ipf].eta())==false && isInPhiCracks((*pfCandidates)[ipf].phi(),(*pfCandidates)[ipf].eta())==false){
    //       fTPhoCone04NeutralHadronIsodR0dEta0pt0nocracks[phoqi] += pt;
    //       if (pt>0.5) fTPhoCone04NeutralHadronIsodR0dEta0pt5nocracks[phoqi] += pt;
    //     }
    //   }
  
    //   if (type==1) { //Charged hadron
    //     fTPhoCone04ChargedHadronIsodR0dEta0pt0dz0old[phoqi] += pt;

    //     if (isPU==0) fTPhoCone04ChargedHadronIsodR0dEta0pt0PFnoPUold[phoqi] += pt;
    //     if (dR>0.015) {
    //       fTPhoCone04ChargedHadronIsodR015dEta0pt0dz0old[phoqi] += pt;
    //       if (isPU==0) fTPhoCone04ChargedHadronIsodR015dEta0pt0PFnoPUold[phoqi] += pt;
    //     }
    //   }

    //   if (type==2) { //Photon
    //     fTPhoCone04PhotonIsodR0dEta0pt0[phoqi] += pt;
    //     if (pt>0.5) fTPhoCone04PhotonIsodR0dEta0pt5[phoqi] += pt;
    //     if (dR>0.08) {
    //       fTPhoCone04PhotonIsodR8dEta0pt0[phoqi] += pt;
    //       if (pt>0.5) fTPhoCone04PhotonIsodR8dEta0pt5[phoqi] += pt;
    //     }
    //   }
  
    //   return;
    // }

    //________________________________________________________________________________________
    reco::VertexRef NTupleProducer::chargedHadronVertex( const edm::Handle<reco::VertexCollection>& vertices, const reco::PFCandidate& pfcand ) const {

      //PfPileUp candidates!

      //  cout << "chargedHadronVertex finding" << endl;

      reco::TrackBaseRef trackBaseRef( pfcand.trackRef() );
  
      size_t  iVertex = 0;
      unsigned index=0;
      unsigned nFoundVertex = 0;
      typedef reco::VertexCollection::const_iterator IV;
      float bestweight=0;
      for(IV iv=vertices->begin(); iv!=vertices->end(); ++iv, ++index) {

        const reco::Vertex& vtx = *iv;
    
        typedef reco::Vertex::trackRef_iterator IT;
    
        // loop on tracks in vertices
        for(IT iTrack=vtx.tracks_begin(); 
            iTrack!=vtx.tracks_end(); ++iTrack) {
	 
          const reco::TrackBaseRef& baseRef = *iTrack;

          // one of the tracks in the vertex is the same as 
          // the track considered in the function
          float w = vtx.trackWeight(baseRef);
          if(baseRef == trackBaseRef ) {
            //select the vertex for which the track has the highest weight
            if (w > bestweight){
              bestweight=w;
              iVertex=index;
              nFoundVertex++;
            }	 	
          }
        }
      }

      if (nFoundVertex>0){
        if (nFoundVertex!=1)
          edm::LogWarning("TrackOnTwoVertex")<<"a track is shared by at least two verteces. Used to be an assert";
        return reco::VertexRef( vertices, iVertex);
      }
      // no vertex found with this track. 

      bool checkClosestZVertex_ = true;

      // optional: as a secondary solution, associate the closest vertex in z
      if ( checkClosestZVertex_ ) {

        double dzmin = 10000;
        double ztrack = pfcand.vertex().z();
        bool foundVertex = false;
        index = 0;
        for(IV iv=vertices->begin(); iv!=vertices->end(); ++iv, ++index) {

          double dz = fabs(ztrack - iv->z());
          if(dz<dzmin) {
            dzmin = dz; 
            iVertex = index;
            foundVertex = true;
          }
        }

        if( foundVertex ) 
          return reco::VertexRef( vertices, iVertex);  

      }


      return reco::VertexRef();
    }

    //________________________________________________________________________________________
    int NTupleProducer::FindPFCandType(int id){

      int type = -1;

      if (id==111 || id==130 || id==310 || id==2112) type=0; //neutral hadrons
      if (fabs(id)==211 || fabs(id)==321 || id==999211 || fabs(id)==2212) type=1; //charged hadrons
      if (id==22) type=2; //photons
      if (fabs(id)==11) type=3; //electrons
      if (fabs(id)==13) type=4; //muons

      return type;
    }

    //________________________________________________________________________________________
    bool NTupleProducer::isInPhiCracks(double phi, double eta){


      // tranform radiants [-pi,pi] in degrees [0,360]
      phi = (phi+TMath::Pi()) *180/TMath::Pi();

      // each supermodule is 20 degrees wide
      Double_t moduleWidth = 20;

      // the first module is centered at phi=0, so the first cracks are +10 and -10
      Double_t phi0 = 10.;

      // set a fiducial cut around the crack of +2-2 degrees
      Double_t fiducialCut = 2.;

      bool OK = false;
      if (fabs(eta)<1.44){
        for (Int_t i = 0 ; i < 18; ++i){
          if ((phi0 + moduleWidth*i -fiducialCut) <= phi && phi <= (phi0 + moduleWidth*i + fiducialCut)) OK = true;
          //        cout << " PHI " << (phi0 + moduleWidth*i -fiducialCut) << " " << phi << " " <<  (phi0 + moduleWidth*i + fiducialCut)  << " " << OK << endl ;
        }
      }

      //  cout << "is in phi crack ? " << OK << endl;
      return OK;
    }

    //________________________________________________________________________________________
    bool NTupleProducer::isInEtaCracks(double eta){

      /*
        Yuri Maravin eta cracks def :
        double emAbsEta = fabs(phoSC_GeomEta);
        pho_isInCrack = 0;
        if ( emAbsEta < 0.018 ||
        (emAbsEta > 0.423 && emAbsEta < 0.461) || 
        (emAbsEta > 0.770 && emAbsEta < 0.806) || 
        (emAbsEta > 1.127 && emAbsEta < 1.163) || 
        (emAbsEta > 1.460 && emAbsEta < 1.558) )
        pho_isInCrack = 1;
      */

      const Int_t nBinsEta = 5;
      Double_t leftEta [nBinsEta]       = {0.00, 0.42, 0.77, 1.13, 1.46};
      Double_t rightEta[nBinsEta]       = {0.02, 0.46, 0.81, 1.16, 9999.};

      bool OK = false;
      if (TMath::Abs(eta)<1.44) {
        for (Int_t i = 0; i< nBinsEta; ++i){
          if (leftEta[i] < TMath::Abs(eta) && TMath::Abs(eta) < rightEta[i] ) OK = true;
        }
      }
      else if (TMath::Abs(eta)>1.44 && TMath::Abs(eta)<1.56) OK = true;
      else if (TMath::Abs(eta)>1.56) OK = false;
      //    cout << leftEta[i] << " " << TMath::Abs(eta) << " " << rightEta[i] <<  " " << OK << endl;

      //  cout << "IS IN CRACK ? " << OK << endl;
      return OK;
    }

    //________________________________________________________________________________________
    bool NTupleProducer::CheckPhotonPFCandOverlap(reco::SuperClusterRef scRef, edm::Handle<reco::PFCandidateCollection>& pfCandidates, int i){

      //This tool can be used in RECO only

      bool isOverlapping = false;
  
      const reco::SuperCluster* supercluster = new reco::SuperCluster(*(scRef));
      //const reco::SuperCluster* supercluster = (const reco::SuperCluster *)gamIter->superCluster());
      std::vector<const reco::SuperCluster*> sc;
      sc.push_back(supercluster);
		 
      for(unsigned iele=0; iele<(*pfCandidates)[i].elementsInBlocks().size(); ++iele) {
        // first get the block 
        reco::PFBlockRef blockRef = (*pfCandidates)[i].elementsInBlocks()[iele].first;
        //
        unsigned elementIndex = (*pfCandidates)[i].elementsInBlocks()[iele].second;
        // check it actually exists 
        if(!blockRef.isNull()){
		 
          // then get the elements of the block
          const edm::OwnVector< reco::PFBlockElement >&  elements = (*blockRef).elements();
		 
          const reco::PFBlockElement & pfbe (elements[elementIndex]); 
          // The first ECAL element should be the cluster associated to the GSF; defined as the seed
          //if(pfbe.type()==reco::PFBlockElement::ECAL || pfbe.type()==reco::PFBlockElement::PS1 || pfbe.type()==reco::PFBlockElement::PS2 || pfbe.type()==reco::PFBlockElement::HCAL)
          //{	  
		 
          reco::PFClusterRef myPFClusterRef = pfbe.clusterRef();
          if(!myPFClusterRef.isNull()){  
		 
            const reco::PFCluster & myPFCluster (*myPFClusterRef);
            //reco::PFCluster & myPFCluster (*myPFClusterRef);
	
            //	cout << "PFcand has a ClusterRef"<<endl;
            int hasOverlap = ClusterClusterMapping::checkOverlap(myPFCluster,sc);
            //	if (hasOverlap==-1) cout << "NO Overlap with E/gamma SC"<<endl;
            //	else cout << "Overlap with E/gamma SC"<<endl;
            if (hasOverlap!=-1) isOverlapping = true;

          }
        }
      }

      return isOverlapping;
	
    }

    // ---------------------------------------------------------------------------------------------------------------------------------------------
    PhotonInfo NTupleProducer::fillPhotonInfos(int p1, bool useAllConvs) 
    {
	
      int iConv1 = useAllConvs ? matchPhotonToConversion(p1) : -1;
	
      if ( iConv1 >= 0) {
        // conversions infos
        return PhotonInfo(p1,
                          TVector3((*fTPhoCaloPositionX)[p1],(*fTPhoCaloPositionY)[p1],(*fTPhoCaloPositionZ)[p1]),
                          TVector3((*fTBeamspotx),(*fTBeamspoty),(*fTBeamspotz)),
                          conv_vtx[iConv1],
                          conv_refitted_momentum[iConv1],
                          (*fTPhoEnergy)[p1],
                          (*fTPhoisEB)[p1],
                          (*fTConvNtracks)[iConv1],
                          (*fTPhoConvValidVtx)[iConv1],
                          (*fTConvChi2Probability)[iConv1],
                          (*fTConvEoverP)[iConv1]
                          );
      } 

	
      return PhotonInfo(p1, 
                        TVector3((*fTPhoCaloPositionX)[p1],(*fTPhoCaloPositionY)[p1],(*fTPhoCaloPositionZ)[p1]),
                        TVector3((*fTBeamspotx),(*fTBeamspoty),(*fTBeamspotz)),
                        pho_conv_vtx[p1],
                        pho_conv_refitted_momentum[p1],
                        (*fTPhoEnergy)[p1],
                        (*fTPhoisEB)[p1],
                        (*fTPhoConvNtracks)[p1],                                                                                                                             
                        (*fTPhoConvValidVtx)[p1],                                                                                                                            
                        (*fTPhoConvChi2Probability)[p1] ,                                                                                                                   
                        (*fTPhoConvEoverP)[p1]                                                                                                                               
                        );
    }

    // ---------------------------------------------------------------------------------------------------------------------------------------------
    std::vector<int> NTupleProducer::HggVertexSelection(HggVertexAnalyzer & vtxAna, HggVertexFromConversions & vtxAnaFromConv, 
                                                        PhotonInfo & pho1, PhotonInfo & pho2, std::vector<std::string> & vtxVarNames, 
                                                        bool useMva=false, TMVA::Reader * tmvaReader=0, std::string tmvaMethod="")
    {
      int p1 = pho1.id(), p2 = pho2.id();
      // assert( p1 == vtxAna.pho1() && p2 == vtxAna.pho2() );
      vtxAna.setPairID(p1,p2);
	
      // preselect vertices : all vertices
      std::vector<int> preselAll;
      for(int i=0; i<(*fTNVrtx) ; i++) {
        preselAll.push_back(i); 
      }

      float zconv = 0; 
      float dzconv = 0;
      std::vector<int> preselConv;

      if ( (pho1.isAConversion() || pho2.isAConversion() ) )  {
	  
        if (pho1.isAConversion()  && !pho2.isAConversion() ){
          zconv  = vtxAnaFromConv.vtxZ(pho1);
          dzconv = vtxAnaFromConv.vtxdZ(pho1);
        }
	  
        if (pho2.isAConversion() && !pho1.isAConversion()){
          zconv  = vtxAnaFromConv.vtxZ(pho2);
          dzconv = vtxAnaFromConv.vtxdZ(pho2);
        }
	  
        if ( pho1.isAConversion() && pho2.isAConversion()){
          float z1  = vtxAnaFromConv.vtxZ(pho1);
          float dz1 = vtxAnaFromConv.vtxdZ(pho1);
            
          float z2  = vtxAnaFromConv.vtxZ(pho2);
          float dz2 = vtxAnaFromConv.vtxdZ(pho2);
            
          zconv  = (z1/dz1/dz1 + z2/dz2/dz2)/(1./dz1/dz1 + 1./dz2/dz2 );  // weighted average
          dzconv = sqrt( 1./(1./dz1/dz1 + 1./dz2/dz2)) ;
        }
	  
        // preselect vertices : only vertices in a window zconv +/- dzconv
        for(int i=0; i < (*fTNVrtx); i++) {
          TVector3 * vtxpos= new TVector3((*fTVrtxX)[i],(*fTVrtxY)[i],(*fTVrtxZ)[i]);
          if ( fabs(zconv - vtxpos->Z() ) < dzconv ) 
            preselConv.push_back(i); 
        }
	  
      }
	
      // ---- METHOD 1 	
      // preselection 
      //	if ( preselConv.size()==0 )
      //  vtxAna.preselection(preselAll);
      //else
      //  vtxAna.preselection(preselConv);
	
      //std::vector<int> rankprod = vtxAna.rankprod(vtxVarNames);
	
      // ---- METHOD 1 

	
      // ---- NEW METHOD 2 (suggested by MarcoP) : first use ranking, then conversions info, e.g. on the N vtxs with best rank
      // preselection  : all vtxs
      std::vector<int> rankprodAll = useMva ? vtxAna.rank(*tmvaReader,tmvaMethod) : vtxAna.rankprod(vtxVarNames);
      int iClosestConv = -1;
      float dminconv = 9999999;
	
      TLorentzVector dipho = TLorentzVector((*fTPhoPx)[p1],(*fTPhoPy)[p1],(*fTPhoPz)[p1],(*fTPhoEnergy)[p1]) 
        + TLorentzVector((*fTPhoPx)[p2],(*fTPhoPy)[p2],(*fTPhoPz)[p2],(*fTPhoEnergy)[p2]);
	
      unsigned int nbest ;
      if (  dipho.Pt() < 30 ) nbest = 5;
      else nbest = 3; 
      if (rankprodAll.size() < nbest ) nbest = rankprodAll.size();

      for (unsigned int ii = 0; ii < nbest; ii++ ){
        TVector3 * vtxpos= new TVector3((*fTVrtxX)[rankprodAll[ii]],(*fTVrtxY)[rankprodAll[ii]],(*fTVrtxZ)[rankprodAll[ii]]);
        if ( fabs( vtxpos->Z()-zconv ) < dzconv && fabs(vtxpos->Z() - zconv ) < dminconv){
          iClosestConv = rankprodAll[ii];
          dminconv = fabs(vtxpos->Z()-zconv );
        }
      }
      std::vector<int> rankprod;
      rankprod.clear();
      if (iClosestConv!=-1 ) rankprod.push_back(iClosestConv);

      //for (int kk = 0; kk  < nbest; kk++ ){
      for (unsigned int kk = 0; kk  < rankprodAll.size(); kk++ ){
        if ( iClosestConv == rankprodAll[kk] ) continue;
        else rankprod.push_back(rankprodAll[kk]);
      }
      // ---- METHOD 2


 
      return rankprod;
    }

    int  NTupleProducer::matchPhotonToConversion( int lpho) {

      int result=-99;
      double conv_eta=-999.;
      double conv_phi=-999.;
             
      int p1=lpho;
  
      float sc_eta  = TVector3((*fTPhoCaloPositionX)[p1],(*fTPhoCaloPositionY)[p1],(*fTPhoCaloPositionZ)[p1]).Eta();
      float  phi  = TVector3((*fTPhoCaloPositionX)[p1],(*fTPhoCaloPositionY)[p1],(*fTPhoCaloPositionZ)[p1]).Phi();
      double sc_phi = phiNorm(phi);

      // FIXME: Not used anywhere
      //   TLorentzVector * p4 = new TLorentzVector((*fTPhoPx)[p1],(*fTPhoPy)[p1],(*fTPhoPz)[p1],(*fTPhoEnergy)[p1]);
      //   float et = (*fTPhoPt)[p1];
  
      float detaMin=999.;
      float dphiMin=999.;   
      float dRMin = 999.;

      float mconv_pt=-999999;
      int iMatch=-1;     
      float conv_pt = -9999;

      // if(LDEBUG)  cout << "   LoopAll::matchPhotonToConversion conv_n " << conv_n << endl; 
      for(int iconv=0; iconv<(*fTNconv); iconv++) {
        TVector3 refittedPairMomentum= conv_refitted_momentum[iconv];
        conv_pt =  refittedPairMomentum.Pt();
        if (conv_pt < 1 ) continue;    
        if ( !(*fTPhoConvValidVtx)[iconv] || (*fTConvNtracks)[iconv]!=2 || (*fTConvChi2Probability)[iconv]<0.000001) continue;

        phi  = conv_refitted_momentum[iconv].Phi();
        conv_phi  = phiNorm(phi);
        float eta  = conv_refitted_momentum[iconv].Eta();
        conv_eta = etaTransformation(eta, (*fTConvZofPrimVtxFromTrks)[iconv] );

        double delta_phi = acos( cos(conv_phi - sc_phi) );       
        double delta_eta = conv_eta - sc_eta;

        delta_phi*=delta_phi;
        delta_eta*=delta_eta;
        float dR = sqrt( delta_phi + delta_eta ); 
    
        if ( fabs(delta_eta) < detaMin && fabs(delta_phi) < dphiMin ) {
          // if ( dR < dRMin ) {
          detaMin=  fabs(delta_eta);
          dphiMin=  fabs(delta_phi);
          dRMin=dR;
          iMatch=iconv;
          mconv_pt = conv_pt;
        }
    
      }
  
      if ( detaMin < 0.1 && dphiMin < 0.1 ) {
        result = iMatch;
      } else {
        result = -1;
      }
  
      return result;
    }

    //________________________________________________________________________________________
    bool NTupleProducer::tkIsHighPurity(reco::TrackRef tk) const { return ( tk->qualityMask() & (1<<2) ) >> 2; }
    bool NTupleProducer::TrackCut(reco::TrackRef tk) const { return false; }
    bool NTupleProducer::ConversionsCut(const reco::Conversion &conv) { 
      //  return (sqrt(conv.refittedPairMomentum().perp2()) <  0); 
      return false;
    }

    //________________________________________________________________________________________
    double NTupleProducer::phiNorm(float &phi) {

      const float pi = 3.1415927;
      const float twopi = 2.0*pi;

      if(phi >  pi) {phi = phi - twopi;}
      if(phi < -pi) {phi = phi + twopi;}

      return phi;
    }

    //________________________________________________________________________________________
    double NTupleProducer::etaTransformation(  float EtaParticle , float Zvertex)  {

      //---Definitions
      const float pi = 3.1415927;

      //---Definitions for ECAL
      const float R_ECAL           = 136.5;
      const float Z_Endcap         = 328.0;
      const float etaBarrelEndcap  = 1.479; 
   
      //---ETA correction

      float Theta = 0.0  ; 
      float ZEcal = R_ECAL*sinh(EtaParticle)+Zvertex;

      if(ZEcal != 0.0) Theta = atan(R_ECAL/ZEcal);
      if(Theta<0.0) Theta = Theta+pi ;
      double ETA = - log(tan(0.5*Theta));
         
      if( fabs(ETA) > etaBarrelEndcap )
        {
          float Zend = Z_Endcap ;
          if(EtaParticle<0.0 )  Zend = -Zend ;
          float Zlen = Zend - Zvertex ;
          float RR = Zlen/sinh(EtaParticle); 
          Theta = atan(RR/Zend);
          if(Theta<0.0) Theta = Theta+pi ;
          ETA = - log(tan(0.5*Theta));		      
        } 
      //---Return the result
      return ETA;
      //---end
    }


    //define this as a plug-in
    DEFINE_FWK_MODULE(NTupleProducer);
