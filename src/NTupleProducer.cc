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
// $Id: NTupleProducer.cc,v 1.171.2.28 2012/12/15 18:30:08 peruzzi Exp $
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

#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "PFIsolation/SuperClusterFootprintRemoval/interface/SuperClusterFootprintRemoval.h"

// Interface
#include "DiLeptonAnalysis/NTupleProducer/interface/NTupleProducer.h"

#include <map>
//#include <pair>
#include <sstream>


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
	fSrcSigma             = iConfig.getUntrackedParameter<edm::InputTag>("tag_srcSigma");
	fSrcRhoPFnoPU       = iConfig.getUntrackedParameter<edm::InputTag>("tag_srcRhoPFnoPU");
	pfphotonsProducerTag = iConfig.getUntrackedParameter<edm::InputTag>("tag_pfphotonsProducer");
	pfProducerTag = iConfig.getUntrackedParameter<edm::InputTag>("tag_pfProducer");
	fSCTagBarrel = iConfig.getUntrackedParameter<edm::InputTag>("tag_SC_barrel");
        fSCTagEndcap = iConfig.getUntrackedParameter<edm::InputTag>("tag_SC_endcap");
	fTrackCollForVertexing = iConfig.getUntrackedParameter<edm::InputTag>("tag_fTrackCollForVertexing");
	fallConversionsCollForVertexing = iConfig.getUntrackedParameter<edm::InputTag>("tag_fallConversionsCollForVertexing");
	perVtxMvaWeights = iConfig.getUntrackedParameter<std::string>("tag_perVtxMvaWeights");
	perVtxMvaMethod = iConfig.getUntrackedParameter<std::string>("tag_perVtxMvaMethod");
	perEvtMvaWeights = iConfig.getUntrackedParameter<std::string>("tag_perEvtMvaWeights");
       	perEvtMvaMethod = iConfig.getUntrackedParameter<std::string>("tag_perEvtMvaMethod");

	doVertexingFlag = iConfig.getUntrackedParameter<bool>("tag_doVertexing");
	if (fIsModelScan) doVertexingFlag=false;

	// Event Selection
	fMinmupt        = iConfig.getParameter<double>("sel_minmupt");
	fMaxmueta       = iConfig.getParameter<double>("sel_maxmueta");
	fMinelpt        = iConfig.getParameter<double>("sel_minelpt");
	fMaxeleta       = iConfig.getParameter<double>("sel_maxeleta");
	fMincorjpt      = iConfig.getParameter<double>("sel_mincorjpt");
	fMinrawjpt      = iConfig.getParameter<double>("sel_minrawjpt");
	fMaxjeta        = iConfig.getParameter<double>("sel_maxjeta");
	fMinjemfrac     = iConfig.getParameter<double>("sel_minjemfrac");
	fMintrkpt       = iConfig.getParameter<double>("sel_mintrkpt");
	fMaxtrketa      = iConfig.getParameter<double>("sel_maxtrketa");
	fMaxtrknchi2    = iConfig.getParameter<double>("sel_maxtrknchi2");
	fMintrknhits    = iConfig.getParameter<int>("sel_mintrknhits");
	fMinphopt       = iConfig.getParameter<double>("sel_minphopt");
	fMaxphoeta      = iConfig.getParameter<double>("sel_maxphoeta");
	fMinSCraw       = iConfig.getParameter<double>("sel_minSCraw");
	fMinSCrawPt     = iConfig.getParameter<double>("sel_minSCrawPt");
	fMingenleptpt   = iConfig.getParameter<double>("sel_mingenleptpt");
	fMaxgenlepteta  = iConfig.getParameter<double>("sel_maxgenlepteta");
	fMingenjetpt    = iConfig.getParameter<double>("sel_mingenjetpt");
	fMaxgenjeteta   = iConfig.getParameter<double>("sel_maxgenjeteta");
	fMinebrechitE   = iConfig.getParameter<double>("sel_fminebrechitE");
	fMingenphotpt   = iConfig.getParameter<double>("sel_mingenphotpt");
	fMaxgenphoteta  = iConfig.getParameter<double>("sel_maxgenphoteta");

	
	if(fIsModelScan) {
		LHAPDF::initPDFSet("cteq66.LHgrid",1);
		NPdfs = LHAPDF::numberPDF();
	}

	CrackCorrFunc    = EcalClusterFunctionFactory::get()->create("EcalClusterCrackCorrection", iConfig);
	LocalCorrFunc    = EcalClusterFunctionFactory::get()->create("EcalClusterLocalContCorrection",iConfig);


	// Create histograms and trees
	fHhltstat        = fTFileService->make<TH1I>("HLTTriggerStats",    "HLTTriggerStatistics",    gMaxhltbits+2,    0, gMaxhltbits+2);
	fHl1physstat     = fTFileService->make<TH1I>("L1PhysTriggerStats", "L1PhysTriggerStatistics", gMaxl1physbits+2, 0, gMaxl1physbits+2);
	fHl1techstat     = fTFileService->make<TH1I>("L1TechTriggerStats", "L1TechTriggerStatistics", gMaxl1techbits+2, 0, gMaxl1techbits+2);
        fHpileupstat     = fTFileService->make<TH1I>("PileUpStats", "PileUpStats", 100, 0, 100 ); // Keep track of pileup distribution
        fHtruepileupstat = fTFileService->make<TH1I>("TruePileUpStats", "TruePileUpStats", 100, 0, 100 ); // Keep track of pileup distribution

	fRunTree         = fTFileService->make<TTree>("RunInfo", "ETHZRunAnalysisTree");
	fEventTree       = fTFileService->make<TTree>("Analysis", "ETHZAnalysisTree");

	// Dump the full configuration
	edm::LogVerbatim("NTP") << "---------------------------------";
	edm::LogVerbatim("NTP") << " ==> NTupleProducer Constructor ...";
	edm::LogVerbatim("NTP") << iConfig;

	// Create additional jet fillers
	std::vector<edm::ParameterSet> jConfigs = iConfig.getParameter<std::vector<edm::ParameterSet> >("jets");
	for (size_t i=0; i<jConfigs.size(); ++i)
          if ( jConfigs[i].getUntrackedParameter<bool>("isPat") ) jetFillers.push_back( new JetFillerPat(jConfigs[i], fEventTree, fIsRealData) );
          else jetFillers.push_back( new JetFillerReco(jConfigs[i], fEventTree, fIsRealData) );

	// Create additional lepton fillers
	std::vector<edm::ParameterSet> lConfigs = iConfig.getParameter<std::vector<edm::ParameterSet> >("leptons");
	for (size_t i=0; i<lConfigs.size(); ++i) {
          std::string type(lConfigs[i].getUntrackedParameter<std::string>("type"));
          if ( type == "electron" ) 
            electronFillers.push_back( new PatElectronFiller(lConfigs[i], fEventTree, fIsRealData) );
          else if ( type == "muon" ) 
            muonFillers.push_back( new PatMuonFiller(lConfigs[i], fEventTree, fIsRealData) );
          else if ( type == "tau" ) 
            tauFillers.push_back( new PatTauFiller(lConfigs[i], fEventTree, fIsRealData) );
        }
        

	// Get list of trigger paths to store the triggering object info. of
	fTHltLabels = iConfig.getUntrackedParameter<std::vector<std::string> >("hlt_labels");
	fTNpaths = fTHltLabels.size();

	// Close your eyes...
	fTHLTObjectID    = new int*[fTNpaths];
	fTHLTObjectID[0] = new int[fTNpaths*gMaxhltnobjs];
	fTHLTObjectPt    = new float*[fTNpaths];
	fTHLTObjectPt[0] = new float[fTNpaths*gMaxhltnobjs];
	fTHLTObjectEta   = new float*[fTNpaths];
	fTHLTObjectEta[0]= new float[fTNpaths*gMaxhltnobjs];
	fTHLTObjectPhi   = new float*[fTNpaths];
	fTHLTObjectPhi[0]= new float[fTNpaths*gMaxhltnobjs];
	for ( size_t i=1; i<fTNpaths; ++i ) {
		fTHLTObjectID[i]  = fTHLTObjectID[i-1]+gMaxhltnobjs;
		fTHLTObjectPt[i]  = fTHLTObjectPt[i-1]+gMaxhltnobjs;
		fTHLTObjectEta[i] = fTHLTObjectEta[i-1]+gMaxhltnobjs;
		fTHLTObjectPhi[i] = fTHLTObjectPhi[i-1]+gMaxhltnobjs;
	}
	
	
	//OOT pu reweighting
        if(!fIsRealData){
		fTPileUpHistoData = iConfig.getUntrackedParameter<std::vector<std::string> >("pu_data");
		fTPileUpHistoMC   = iConfig.getUntrackedParameter<std::vector<std::string> >("pu_mc");
	 	if(! fTPileUpHistoData[0].empty() && !fTPileUpHistoMC[0].empty() ){
			LumiWeights_      = edm::LumiReWeighting(fTPileUpHistoMC[0], fTPileUpHistoData[0], fTPileUpHistoMC[1], fTPileUpHistoData[1]);
		}
	}

//	vAna = new HggVertexAnalyzer(vtxAlgoParams);
//	vConv= new HggVertexFromConversions(vtxAlgoParams);
//	if (doVertexingFlag) vAna->setupWithDefaultOptions(perVtxMvaWeights, perEvtMvaWeights, rankVariables, perVtxReader, perVtxMvaMethod, perEvtReader, perEvtMvaMethod);

}

//________________________________________________________________________________________
NTupleProducer::~NTupleProducer(){
	delete [] fTHLTObjectID[0];
	delete [] fTHLTObjectPt[0];
	delete [] fTHLTObjectEta[0];
	delete [] fTHLTObjectPhi[0];
	delete [] fTHLTObjectID;
	delete [] fTHLTObjectPt;
	delete [] fTHLTObjectEta;
	delete [] fTHLTObjectPhi;
}

//________________________________________________________________________________________
// Method called once for each event
void NTupleProducer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){

	++fNTotEvents;

	using namespace edm;
	using namespace std;
	using namespace reco;
	using reco::MuonCollection;
	using reco::JetTagCollection;

	// Reset all the tree variables
	resetTree();
	for ( std::vector<JetFillerBase*>::iterator it = jetFillers.begin(); 
              it != jetFillers.end(); ++it ) 
          (*it)->reset();
	for ( std::vector<PatMuonFiller*>::iterator it = muonFillers.begin(); 
              it != muonFillers.end(); ++it ) 
          (*it)->reset();
	for ( std::vector<PatElectronFiller*>::iterator it = electronFillers.begin(); 
              it != electronFillers.end(); ++it ) 
          (*it)->reset();
	for ( std::vector<PatTauFiller*>::iterator it = tauFillers.begin(); 
              it != tauFillers.end(); ++it ) 
          (*it)->reset();


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
	fTrho = *rho;

	// sigma for L1FastJet
	edm::Handle<double> sigma;
	iEvent.getByLabel(fSrcSigma,sigma);
	fTsigma = *sigma;
	
	// rho for L1FastJet running PFnoPU
	edm::Handle<double> rhoNoPU;
	iEvent.getByLabel(fSrcRhoPFnoPU,rhoNoPU);
	fTrhoPFnoPU = *rhoNoPU;

	// beam halo
	edm::Handle<BeamHaloSummary> TheBeamHaloSummary;
	iEvent.getByLabel("BeamHaloSummary",TheBeamHaloSummary);
	const BeamHaloSummary TheSummary = (*TheBeamHaloSummary.product());
	fTcscTightHaloID = static_cast<int>  (TheSummary.CSCTightHaloId());

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

	// All conversions
	edm::Handle<reco::ConversionCollection> convH;
	iEvent.getByLabel(fallConversionsCollForVertexing, convH);

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
	const reco::BeamSpot &beamspot = *beamSpotHandle.product();

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
	const CaloSubdetectorGeometry *barrelGeometry = geometry->getSubdetectorGeometry(DetId::Ecal, EcalBarrel);
	const CaloSubdetectorGeometry *endcapGeometry = geometry->getSubdetectorGeometry(DetId::Ecal, EcalEndcap);

	edm::ESHandle<CaloTopology> theCaloTopo;
	iSetup.get<CaloTopologyRecord>().get(theCaloTopo);
	const CaloTopology *topology = theCaloTopo.product();

	CrackCorrFunc->init(iSetup);
	LocalCorrFunc->init(iSetup);

	// ECAL dead cell Trigger Primitive filter
	edm::Handle<bool> EcalDeadTPFilterFlag;
	iEvent.getByLabel("ecalDeadCellTPfilter",EcalDeadTPFilterFlag);
	fTecalDeadTPFilterFlag = (int) *EcalDeadTPFilterFlag;
	
	// Stevens recovRecHitFilter
	// http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/UserCode/lowette/SandBox/Skims/python/recovRecHitFilter_cfi.py?sortby=date&view=log
	// https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFilters#RecovRecHitFilter
	edm::Handle<bool> RecovRecHitFilterFlag;
	iEvent.getByLabel("recovRecHitFilter","Result",RecovRecHitFilterFlag);
	fRecovRecHitFilterFlag =(int) *RecovRecHitFilterFlag;

	// RA2 tracking tailure filter
	// https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFilters
	// http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/UserCode/seema/SandBox/Skims/python/trackingFailureFilter_cfi.py?hideattic=0&revision=1.3&view=markup&pathrev=MAIN
	edm::Handle<bool> RA2TrackingFailureFlag;
	iEvent.getByLabel("trackingFailureFilter",RA2TrackingFailureFlag);
	fTra2TrackingFailureFilterFlag = (int) *RA2TrackingFailureFlag;

	// Colin's PBNR filter
	edm::Handle<bool> ParticleBasedNoiseRejectionFlag;
	//FR: this requires modifications in the CVS code.
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
	   fTHBHENoiseFlag    = static_cast<int>(*hbHeNoiseFlag);
	   
	   edm::Handle<bool> hbHeNoiseFlagIso;
	   iEvent.getByLabel(fHBHENoiseResultTagIso,hbHeNoiseFlagIso);
	   fTHBHENoiseFlagIso = static_cast<int>(*hbHeNoiseFlagIso);
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
			fTqcdPartonicHT=partonicHT;
		}

		iEvent.getByLabel("generator", genEvtInfo);
		fTpthat       = genEvtInfo->hasBinningValues() ? (genEvtInfo->binningValues())[0] : 0.0;
		fTsigprocid   = genEvtInfo->signalProcessID();
		fTpdfscalePDF = genEvtInfo->pdf()->scalePDF;
		fTpdfid1      = genEvtInfo->pdf()->id.first;
		fTpdfid2      = genEvtInfo->pdf()->id.second;
		fTpdfx1       = genEvtInfo->pdf()->x.first;
		fTpdfx2       = genEvtInfo->pdf()->x.second;
		fTpdfxPDF1    = genEvtInfo->pdf()->xPDF.first;
		fTpdfxPDF2    = genEvtInfo->pdf()->xPDF.second;
		fTgenweight   = genEvtInfo->weight();
     
		iEvent.getByLabel("addPileupInfo", pileupInfo);
		std::vector<PileupSummaryInfo>::const_iterator PVI;

		for (PVI = pileupInfo->begin(); PVI !=pileupInfo->end(); ++PVI){
		  if( PVI->getBunchCrossing() == 0 ){ // in-time PU
		    fTpuNumInteractions  = PVI->getPU_NumInteractions();
		    fHpileupstat->Fill( fTpuNumInteractions );
		    fTpuNumTrueInteractions  = PVI->getTrueNumInteractions();
		    fHtruepileupstat->Fill( fTpuNumTrueInteractions );
		    
		    if(fTpuNumInteractions > gMaxnpileup){
		      edm::LogWarning("NTP") << "@SUB=analyze()"
					     << "More than " << static_cast<int>(gMaxnpileup)
					     << " generated Pileup events found, increase size!";
		      fTgoodevent = 1;
		    }
		    
		    fTpuNumFilled = (int)PVI->getPU_zpositions().size();
		    for( int i = 0; i < fTpuNumFilled; i++){
		      if(i >= gMaxnpileup) break; // hard protection
		      fTpuZpositions[i]   = PVI->getPU_zpositions()[i];
		      fTpuSumpT_lowpT[i]  = PVI->getPU_sumpT_lowpT()[i];
		      fTpuSumpT_highpT[i] = PVI->getPU_sumpT_highpT()[i];
		      fTpuNtrks_lowpT[i]  = PVI->getPU_ntrks_lowpT()[i];
		      fTpuNtrks_highpT[i] = PVI->getPU_ntrks_highpT()[i];
		      
		      // fTpuInstLumi[i]     = PVI->getPU_instLumi()[i];
		    }
		    
		  }
		  else if( PVI->getBunchCrossing() == 1 ){ // OOT pile-Up: this is the 50ns late Bunch
		    fTpuOOTNumInteractionsLate = PVI->getPU_NumInteractions();
		  }
		  else if( PVI->getBunchCrossing() == -1 ){ // OOT pile-Up: this is the 50ns early Bunch
		    fTpuOOTNumInteractionsEarly = PVI->getPU_NumInteractions();
		  }
		  
		  
		}
		//see https://twiki.cern.ch/twiki/bin/view/CMS/PileupMCReweightingUtilities 
		// as well as http://cmslxr.fnal.gov/lxr/source/PhysicsTools/Utilities/src/LumiReWeighting.cc
		if(!fTPileUpHistoData[0].empty() && !fTPileUpHistoMC[0].empty() ){
		  //const EventBase* iEventB = dynamic_cast<const EventBase*>(&iEvent);
		  MyWeightTotal  = LumiWeights_.weightOOT( iEvent ); // this is the total weight inTimeWeight * WeightOOTPU * Correct_Weights2011
		  MyWeightInTime = LumiWeights_.weight   ( iEvent ); // this is the inTimeWeight only
		}
		if(fIsModelScan) {
			edm::Handle<GenEventInfoProduct> pdfstuff;
			if (!iEvent.getByLabel("generator", pdfstuff)) {
				edm::LogError("PDFWeightProducer") << ">>> PdfInfo not found !!!";
				return;
			}
		
			float Q = pdfstuff->pdf()->scalePDF;
			int id1 = pdfstuff->pdf()->id.first;
			double x1 = pdfstuff->pdf()->x.first;
			double pdf1 = pdfstuff->pdf()->xPDF.first;
			int id2 = pdfstuff->pdf()->id.second;
			double x2 = pdfstuff->pdf()->x.second;
			double pdf2 = pdfstuff->pdf()->xPDF.second;

			fTpdfW[0]=1;
			LHAPDF::initPDF(0);

			double newpdf1_0 = LHAPDF::xfx(x1, Q, id1)/x1;
			double newpdf2_0 = LHAPDF::xfx(x2, Q, id2)/x2;
			
			float pdfWsum=0;
			for(int pdf=1; pdf<= NPdfs; pdf++){
				LHAPDF::initPDF(pdf);
				double newpdf1 = LHAPDF::xfx(x1, Q, id1)/x1;
				double newpdf2 = LHAPDF::xfx(x2, Q, id2)/x2;
				fTpdfW[pdf] = newpdf1/newpdf1_0*newpdf2/newpdf2_0;
				pdfWsum+=fTpdfW[pdf];
			}
			fTpdfWsum=pdfWsum;

			int process = 0;
			if(fTsigprocid>=237 && fTsigprocid<=242) process=1;//"ng";
			else if(fTsigprocid>=246 && fTsigprocid<=256) process=2;//"ns";
			else if(fTsigprocid>=216 && fTsigprocid<=236) process=3;//"nn";
			else if(fTsigprocid>=201 && fTsigprocid<=214) process=4;//"ll";
			else if( (fTsigprocid>=274 && fTsigprocid<=280) || (fTsigprocid>=284 && fTsigprocid<=286)) process=5;//"sb";
			else if( (fTsigprocid>=271 && fTsigprocid<=273) || (fTsigprocid>=281 && fTsigprocid<=283) || (fTsigprocid>=291 && fTsigprocid<=293)) process=6;//"ss";
			else if(fTsigprocid>=261 && fTsigprocid<=265) process=7;//"tb";
			else if( (fTsigprocid>=287 && fTsigprocid<=290) || fTsigprocid==296) process=8;//"bb";
			else if(fTsigprocid>=243 && fTsigprocid<=244) process=9;//"gg";
			else if( (fTsigprocid>=258 && fTsigprocid<=259) || (fTsigprocid>=294 && fTsigprocid<=295) ) process=10;//"sg";
			fTprocess=process;
		}
	} else {
          // Just store the number of primary vertices
          fHpileupstat->Fill(vertices->size());
        }
	fTpuWeightTotal  = MyWeightTotal;
	fTpuWeightInTime = MyWeightInTime;
	
	
	//////////////////////////////////////////////////////////////////////////////
	// Trigger information
	Handle<L1GlobalTriggerReadoutRecord> l1GtReadoutRecord;
	iEvent.getByLabel(fL1TriggerTag, l1GtReadoutRecord);

	Handle<trigger::TriggerEvent> triggerEventHLT;
	iEvent.getByLabel("hltTriggerSummaryAOD", triggerEventHLT);

	fProcessName = triggerEventHLT.provenance()->processName();
	InputTag trigResultsTag("TriggerResults","",fProcessName);

	Handle<TriggerResults> triggers;
	iEvent.getByLabel(trigResultsTag, triggers);
	const TriggerResults& tr = *triggers;

	bool changed(true);
	// Trigger menu could have changed
	if (fHltConfig.init(iEvent.getRun(),iSetup,fProcessName,changed)) {
		if ( changed ) { //fixme: could book trigger histos here
		}
	} else {
		// problem: init failed
		edm::LogError("NTP") << " hlt config extraction failure with process name "
			<< fProcessName;
	}

	// Get trigger menus
	if(tr.size() >= gMaxhltbits){
		edm::LogWarning("NTP") << "@SUB=analyze()"
			<< "More than " << static_cast<int>(gMaxhltbits)
			<< " HLT trigger bits, increase length!";
		fTgoodevent = 1;
	}
	std::vector<string> triggernames;
	triggernames.reserve(tr.size());
	Service<service::TriggerNamesService> tns;
	tns->getTrigPaths(*triggers, triggernames);
	fTHLTmenu = triggernames; // Used below and stored in run tree at end of run

	if( fFirstevent ){ // do these things only once
		for( unsigned int i = 0; i < tr.size(); i++ ) fHhltstat->GetXaxis()->SetBinLabel(i+1, TString(triggernames[i]));

		// Add a warning about the shift between trigger bits and bin numbers:
		fHhltstat->GetXaxis()->SetBinLabel(gMaxhltbits+2, "Bin#=Bit#+1");

		// Add algorithm names to L1 physical bits in trigger histogram
		ESHandle<L1GtTriggerMenu> menuRcd;
		iSetup.get<L1GtTriggerMenuRcd>().get(menuRcd);
		const L1GtTriggerMenu *menu = menuRcd.product();
		const AlgorithmMap& algoMap = menu->gtAlgorithmMap();
		for( AlgorithmMap::const_iterator it = algoMap.begin(); it != algoMap.end(); ++it ){
			fHl1physstat->GetXaxis()->SetBinLabel((*it).second.algoBitNumber() + 1, TString((*it).first));
			fTL1physmenu[(*it).second.algoBitNumber()] = std::string((*it).first); // Actually stored in run tree (but easier so)
		}
	}

	fFirstevent = false;

	// Get trigger results and prescale
	for(unsigned int i = 0; i < tr.size(); i++ ){
          bool fired = tr[i].accept();
          if(fired) fHhltstat->Fill(i);
          fTHLTres[i] = fired ? 1:0;
          fTHLTprescale[i] = fHltConfig.prescaleValue(iEvent, iSetup, fTHLTmenu[i]);
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
	for( unsigned int i = 0; i < gMaxl1physbits; ++i ){
		bool fired = l1GtReadoutRecord->decisionWord()[i];
		if(fired) fHl1physstat->Fill(i);
		fTL1physres[i] = fired ? 1:0;
	}
	for( unsigned int i = 0; i < gMaxl1techbits; ++i){
		bool fired = l1GtReadoutRecord->technicalTriggerWord()[i];
		if(fired) fHl1techstat->Fill(i);
		fTL1techres[i] = fired ? 1:0;
	}

	// Store information for some trigger paths
	edm::Handle<trigger::TriggerEvent> trgEvent;
	iEvent.getByLabel(fHLTTrigEventTag, trgEvent);

	InputTag collectionTag;
	// Loop over path names and get related objects
	for (size_t i=0; i<fTNpaths; ++i) {
		collectionTag = edm::InputTag(fTHltLabels[i],"",fProcessName);
		size_t  filterIndex_ = trgEvent->filterIndex(collectionTag);
		if (filterIndex_<trgEvent->sizeFilters()) {
			const trigger::TriggerObjectCollection& TOC(trgEvent->getObjects());
			const trigger::Keys& keys = trgEvent->filterKeys(filterIndex_);
		// Loop over objects
			for ( size_t hlto = 0; hlto<keys.size(); ++hlto ) {
				if (hlto>=gMaxhltnobjs) {
					edm::LogWarning("NTP") << "@SUB=analyze()"
						<< "Maximum number of triggering objects exceeded"
						<< " for filter " << fTHltLabels[i];
					break;
				}
			// Update number of objects stored
				if ( hlto>=fTNHLTobjects ) fTNHLTobjects = hlto+1; // Not an index...
				const trigger::TriggerObject& TO(TOC[keys[hlto]]);
				fTHLTObjectID[i][hlto]  = TO.id();
				fTHLTObjectPt[i][hlto]  = TO.pt();
				fTHLTObjectEta[i][hlto] = TO.eta();
				fTHLTObjectPhi[i][hlto] = TO.phi();
			}
		}
	}

	////////////////////////////////////////////////////////////////////////////////
	// Dump tree variables /////////////////////////////////////////////////////////
	fTrunnumber   = iEvent.id().run();
	fTeventnumber = iEvent.id().event();
	fTlumisection = iEvent.luminosityBlock();

	fTweight = 1.0; // To be filled at some point?

	// Save position of primary vertex
	fTprimvtxx   = primVtx->x();
	fTprimvtxy   = primVtx->y();
	fTprimvtxz   = primVtx->z();
	fTprimvtxrho = primVtx->position().rho();
	fTprimvtxxE  = primVtx->xError();
	fTprimvtxyE  = primVtx->yError();
	fTprimvtxzE  = primVtx->zError();
	fTpvtxznchi2 = primVtx->normalizedChi2();
	fTpvtxndof   = primVtx->ndof();
	fTpvtxisfake = primVtx->isFake();

	fTpvtxptsum = 0.;
	for(std::vector<TrackBaseRef>::const_iterator trackit = primVtx->tracks_begin(); trackit != primVtx->tracks_end(); ++trackit){
		fTpvtxptsum += (*trackit)->pt();
	}
	fTgoodvtx = 0;

	// get all vertices
	int countVrtx = -1;
	for(VertexCollection::const_iterator vertexit = vertices->begin(); vertexit != vertices->end(); ++vertexit) {
		countVrtx++;
		if(countVrtx >= gMaxnvrtx){
			edm::LogWarning("NTP") << "@SUB=analyze()"
				<< "Maximum number of vertices exceeded";
			fTgoodevent = 1;
			fTflagmaxvrtxexc = 1;
			break;
		}

		fTvrtxx     [countVrtx] = vertexit->x();
		fTvrtxy     [countVrtx] = vertexit->y();
		fTvrtxz     [countVrtx] = vertexit->z();
		fTvrtxxE    [countVrtx] = vertexit->xError();
		fTvrtxyE    [countVrtx] = vertexit->yError();
		fTvrtxzE    [countVrtx] = vertexit->zError();
		fTvrtxndof  [countVrtx] = vertexit->ndof();
		fTvrtxchi2  [countVrtx] = vertexit->normalizedChi2();
		fTvrtxntrks [countVrtx] = vertexit->nTracks();
		fTvrtxsumpt [countVrtx] = vertexit->p4().pt();
		fTvrtxisfake[countVrtx] = vertexit->isFake();
	}
	fTnvrtx = countVrtx+1;


	// Save position of beamspot
	fTbeamspotx = (beamSpot.position()).x();
	fTbeamspoty = (beamSpot.position()).y();
	fTbeamspotz = (beamSpot.position()).z();

	IndexByPt indexComparator; // Need this to sort collections

	/////////////////////////////////////////
	/// GenVertices 

	if (!fIsRealData && doVertexingFlag){

	  edm::Handle<reco::GenParticleCollection> gpH;
	  iEvent.getByLabel(fGenPartTag, gpH);

	  const float lowPtThrGenVtx = 0.1;
	  const float highPtThrGenVtx = 0.5;

	  gv_n = 0;

	  for(reco::GenParticleCollection::const_iterator it_gen = 
		gpH->begin(); it_gen!= gpH->end(); it_gen++){   

	    if (gv_n>=gMaxngenvtx){
	      edm::LogWarning("NTP") << "@SUB=analyze"
				     << "Maximum number of gen-vertices exceeded..";
	      fTgoodevent = 1;
	      break;
	    }

	    if( it_gen->status() != 3 || !(it_gen->vx()!=0. || it_gen->vy()!=0. || it_gen->vz()!=0.)  ) { continue; }

	    // check for duplicate vertex
	    bool duplicate = false;
	    for(Int_t itv = 0; itv < gv_n; itv++) {
	      TVector3 checkVtx = gv_pos[itv];
	      if( (fabs(it_gen->vx()-checkVtx.X())<1e-5) &&  (fabs(it_gen->vy()-checkVtx.Y())<1e-5) && (fabs(it_gen->vz()-checkVtx.Z())<1e-5)) {
		duplicate = true;
		break;
	      }
	    }

	    if (duplicate) continue;
    
	    gv_pos[gv_n].SetXYZ(it_gen->vx(), it_gen->vy(), it_gen->vz());
    
	    TVector3  this_gv_pos = gv_pos[gv_n];
	    TVector3 p3(0,0,0);
    
	    gv_sumPtLo[gv_n] = 0;
	    gv_nTkLo[gv_n] = 0;
	    gv_sumPtHi[gv_n] = 0;
	    gv_nTkHi[gv_n] = 0;

	    for(reco::GenParticleCollection::const_iterator part = gpH->begin(); part!= gpH->end(); part++){   
	      if (part->pt()==0) continue;
	      if( part->status() == 1 && part->charge() != 0 && fabs(part->eta())<2.5 &&
		  ( fabs(part->vx()-this_gv_pos.X())<1.e-5 && fabs(part->vy()-this_gv_pos.Y())<1.e-5 && fabs(part->vz()-this_gv_pos.Z())<1.e-5 ) )  {
	
		TVector3 m(part->px(),part->py(),part->pz());
		p3 += m;
		if( m.Pt() > lowPtThrGenVtx ) {
		  gv_sumPtLo[gv_n] += m.Pt();
		  gv_nTkLo[gv_n] += 1;
		  if( m.Pt() > highPtThrGenVtx ) {
		    gv_sumPtHi[gv_n] += m.Pt();
		    gv_nTkHi[gv_n] += 1;
		  }
		}
	      }
	    }

	    gv_p3[gv_n].SetXYZ(p3.X(),p3.Y(),p3.Z());

	    gv_n++;
	  }

	} // end gen vertices



	////////////////////////////////////////////////////////////////////////////////
	// Get GenLeptons (+ Mother and GMother)
	if(!fIsRealData){
		edm::Handle<GenParticleCollection> gen;
		iEvent.getByLabel(fGenPartTag, gen);
		GenParticleCollection::const_iterator g_part;

		std::vector<const GenParticle*> gen_lepts;
		std::vector<const GenParticle*> gen_moms;
		std::vector<const GenParticle*> gen_gmoms;


		// loop over genparticles to get gen_els and gen_mus
		for(g_part = gen->begin(); g_part != gen->end(); g_part++){
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

			if( g_part->pt()        < fMingenleptpt )  continue;
			if( fabs(g_part->eta()) > fMaxgenlepteta ) continue;

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
				<< " WARNING: GenParticle does not have a GrandMother ";
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
			fTngenleptons = gen_lepts.size();
			for(int i=0; i<fTngenleptons; ++i){
				if( i >= gMaxngenlept){
					edm::LogWarning("NTP") << "@SUB=analyze"
						<< "Maximum number of gen-leptons exceeded...";
					fTflagmaxgenleptexc = 1;
					fTgoodevent = 1;
					break;
				}

				fTGenLeptonId[i]       =   gen_lepts[i]->pdgId();
				fTGenLeptonPt[i]       =   gen_lepts[i]->pt();
				fTGenLeptonEta[i]      =   gen_lepts[i]->eta();
				fTGenLeptonPhi[i]      =   gen_lepts[i]->phi();

				fTGenLeptonMId[i]      =   (gen_moms[i]!=NULL )  ? gen_moms[i]->pdgId()  : -999;
				fTGenLeptonMStatus[i]  =   (gen_moms[i]!=NULL )  ? gen_moms[i]->status() : -999;
				fTGenLeptonMPt[i]      =   (gen_moms[i]!=NULL )  ? gen_moms[i]->pt()     : -999;
				fTGenLeptonMEta[i]     =   (gen_moms[i]!=NULL )  ? gen_moms[i]->eta()    : -999;
				fTGenLeptonMPhi[i]     =   (gen_moms[i]!=NULL )  ? gen_moms[i]->phi()    : -999;

				fTGenLeptonGMId[i]     =   (gen_gmoms[i]!=NULL ) ? gen_gmoms[i]->pdgId() : -999;
				fTGenLeptonGMStatus[i] =   (gen_gmoms[i]!=NULL ) ? gen_gmoms[i]->status(): -999;
				fTGenLeptonGMPt[i]     =   (gen_gmoms[i]!=NULL ) ? gen_gmoms[i]->pt()    : -999;
				fTGenLeptonGMEta[i]    =   (gen_gmoms[i]!=NULL ) ? gen_gmoms[i]->eta()   : -999;
				fTGenLeptonGMPhi[i]    =   (gen_gmoms[i]!=NULL ) ? gen_gmoms[i]->phi()   : -999;
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
           if( g_part->pt() < fMingenphotpt )  continue;
           if( fabs(g_part->eta()) > fMaxgenphoteta ) continue;

           const GenParticle* gen_phot = &(*g_part);
           const GenParticle* gen_phot_mom = static_cast<const GenParticle*> (g_part->mother());

           if(gen_phot_mom==NULL){
             edm::LogWarning("NTP") << "@SUB=analyze" << " WARNING: GenPhoton does not have a mother ";
           }

           gen_photons.push_back(gen_phot);
           gen_photons_mothers.push_back(gen_phot_mom);

         }

         fTngenphotons = gen_photons.size();

         for(int i=0; i<fTngenphotons; ++i){
           if( i >= gMaxngenphot){
             edm::LogWarning("NTP") << "@SUB=analyze" << "Maximum number of gen-photons exceeded..";
             fTflagmaxgenphotexc = 1;
             fTgoodevent = 1;
             break;
           }

           fTGenPhotonPt[i]       =   gen_photons[i]->pt();
           fTGenPhotonEta[i]      =   gen_photons[i]->eta();
           fTGenPhotonPhi[i]      =   gen_photons[i]->phi();
           fTGenPhotonMotherID[i] =   gen_photons_mothers[i]!=NULL ? gen_photons_mothers[i]->pdgId() : -999;
           fTGenPhotonMotherStatus[i] = gen_photons_mothers[i]!=NULL ? gen_photons_mothers[i]->status() : -999;

	   fTGenPhotonIsoDR03[i] = GenPartonicIso_allpart(*(gen_photons[i]),gen,0.3);
	   fTGenPhotonIsoDR04[i] = GenPartonicIso_allpart(*(gen_photons[i]),gen,0.4);


	   // use Steve Mrenna's status 2 parton jets to compute dR to closest jet of prompt photon
	   if(fTGenPhotonMotherStatus[i]!=3) continue;
	   TLorentzVector photon(0,0,0,0);
	   photon.SetPtEtaPhiM(fTGenPhotonPt[i],fTGenPhotonEta[i],fTGenPhotonPhi[i],0);
	   double minDR=10;
	   for(pGenJet = partonGenJets->begin(); pGenJet != partonGenJets->end(); pGenJet++){
		TLorentzVector pJ(pGenJet->px(), pGenJet->py(), pGenJet->pz(), pGenJet->energy());
		float dR = photon.DeltaR(pJ);
		if(dR < minDR) minDR=dR;
	   }
	   fTGenPhotonPartonMindR[i] = minDR;

	 
         }
	 
       }


	////////////////////////////////////////////////////////////////////////////////
	// Get GenJets
	if(!fIsRealData){
		edm::Handle<GenJetCollection> genjets;
		iEvent.getByLabel(fGenJetTag, genjets);
		GenJetCollection::const_iterator gjet;
		
		int jqi=-1;
		for(gjet = genjets->begin(); gjet != genjets->end(); gjet++){
			// Preselection
			if(gjet->pt() < fMingenjetpt) continue;
			if(fabs(gjet->eta()) > fMaxgenjeteta) continue;
			jqi++;
			if( jqi >= gMaxngenjets){
				edm::LogWarning("NTP") << "@SUB=analyze"
					<< "Maximum number of gen-jets exceeded..";
				fTflagmaxgenjetexc = 1;
				fTgoodevent = 1;
				break;
			}
			
			fTGenJetPt  [jqi] = gjet->pt();
			fTGenJetEta [jqi] = gjet->eta();
			fTGenJetPhi [jqi] = gjet->phi();
			fTGenJetE   [jqi] = gjet->energy();
			fTGenJetemE [jqi] = gjet->emEnergy();
			fTGenJethadE[jqi] = gjet->hadEnergy();
			fTGenJetinvE[jqi] = gjet->invisibleEnergy();
		}
		fTNGenJets = jqi+1;
	}


	////////////////////////////////////////////////////////
	// Muon Variables:
	int mqi(0);  // Index of qualified muons
	fTnmutot = 0; // Total number of tracker&&global muons

	// Get muons, order them by pt and apply selection
	std::vector<OrderPair> muOrdered;
	int muIndex(0);
	for ( View<Muon>::const_iterator Mit = muons->begin(); Mit != muons->end();
	++Mit,++muIndex ) {
		// Check if maximum number of muons is exceeded already:
		if(mqi >= gMaxnmus){
			edm::LogWarning("NTP") << "@SUB=analyze()"
				<< "Maximum number of muons exceeded";
			fTflagmaxmuexc = 1;
			fTgoodevent = 1;
			break;
		}
		// Muon preselection:
		// Only consider global and trackermuons:
		if(!(Mit->isGlobalMuon()) && !(Mit->isTrackerMuon())) continue;
		fTnmutot++;     // Count all
		if(Mit->pt() < fMinmupt) continue;
		if(fabs(Mit->eta()) > fMaxmueta) continue;
		++mqi;          // Count how many we'll eventually store
		muOrdered.push_back(make_pair(muIndex,Mit->pt()));
	}
	std::sort(muOrdered.begin(),muOrdered.end(),indexComparator);
	fTnmu = muOrdered.size();
	mqi = 0;

	// Dump muon properties in tree variables
	for (std::vector<OrderPair>::const_iterator it = muOrdered.begin();
	it != muOrdered.end(); ++it, ++mqi ) {
		int index = it->first;
		const Muon& muon = (*muons)[index];

		fTmuIsGM[mqi]   = muon.isGlobalMuon() ? 1:0;
		fTmuIsTM[mqi]   = muon.isTrackerMuon() ? 1:0;

		// Combined methods for Global and Tracker muons:
		fTmupx[mqi]      = muon.px();
		fTmupy[mqi]      = muon.py();
		fTmupz[mqi]      = muon.pz();
		fTmupt[mqi]      = muon.pt();
		fTmuinnerpt[mqi] = muon.innerTrack()->pt();
		fTmueta[mqi]     = muon.eta();
		fTmuphi[mqi]     = muon.phi();
		fTmue[mqi]       = muon.energy();
		fTmuet[mqi]      = muon.et();
		fTmucharge[mqi]  = muon.charge();

		fTmuiso[mqi]            = (muon.isolationR03().sumPt + muon.isolationR03().emEt + muon.isolationR03().hadEt) / muon.pt();
		fTmuIso03sumPt[mqi]     = muon.isolationR03().sumPt;
		fTmuIso03emEt[mqi]      = muon.isolationR03().emEt;
		fTmuIso03hadEt[mqi]     = muon.isolationR03().hadEt;
		fTmuIso03emVetoEt[mqi]  = muon.isolationR03().emVetoEt;
		fTmuIso03hadVetoEt[mqi] = muon.isolationR03().hadVetoEt;
		fTmuIso05sumPt[mqi]     = muon.isolationR05().sumPt;
		fTmuIso05emEt[mqi]      = muon.isolationR05().emEt;
		fTmuIso05hadEt[mqi]     = muon.isolationR05().hadEt;

		fTmucalocomp[mqi] = muon.caloCompatibility();
		fTmusegmcomp[mqi] = muon::segmentCompatibility(muon);

		// MuID Flags:
		fTmuIsGMPT[mqi]                  = muon::isGoodMuon(muon, muon::GlobalMuonPromptTight) ? 1:0;
		fTmuIsGMTkChiComp[mqi]           = muon::isGoodMuon(muon, muon::GMTkChiCompatibility) ? 1:0;
		fTmuIsGMStaChiComp[mqi]          = muon::isGoodMuon(muon, muon::GMStaChiCompatibility) ? 1:0;
		fTmuIsGMTkKinkTight[mqi]         = muon::isGoodMuon(muon, muon::GMTkKinkTight) ? 1:0;
		fTmuIsAllStaMuons[mqi]           = muon::isGoodMuon(muon, muon::AllStandAloneMuons) ? 1:0;
		fTmuIsAllTrkMuons[mqi]           = muon::isGoodMuon(muon, muon::AllTrackerMuons) ? 1:0;
		fTmuIsTrkMuArb[mqi]              = muon::isGoodMuon(muon, muon::TrackerMuonArbitrated) ? 1:0;
		fTmuIsAllArb[mqi]                = muon::isGoodMuon(muon, muon::AllArbitrated) ? 1:0;
		fTmuIsTMLastStationLoose[mqi]    = muon::isGoodMuon(muon, muon::TMLastStationLoose) ? 1:0;
		fTmuIsTMLastStationTight[mqi]    = muon::isGoodMuon(muon, muon::TMLastStationTight) ? 1:0;
		fTmuIsTM2DCompLoose[mqi]         = muon::isGoodMuon(muon, muon::TM2DCompatibilityLoose) ? 1:0;
		fTmuIsTM2DCompTight[mqi]         = muon::isGoodMuon(muon, muon::TM2DCompatibilityTight) ? 1:0;
		fTmuIsTMOneStationLoose[mqi]     = muon::isGoodMuon(muon, muon::TMOneStationLoose) ? 1:0;
		fTmuIsTMOneStationTight[mqi]     = muon::isGoodMuon(muon, muon::TMOneStationTight) ? 1:0;
		fTmuIsTMLSOPL[mqi]               = muon::isGoodMuon(muon, muon::TMLastStationOptimizedLowPtLoose) ? 1:0;
		fTmuIsTMLastStationAngLoose[mqi] = muon::isGoodMuon(muon, muon::TMLastStationAngLoose) ? 1:0;
		fTmuIsTMLastStationAngTight[mqi] = muon::isGoodMuon(muon, muon::TMLastStationAngTight) ? 1:0;
		fTmuIsTMOneStationAngLoose[mqi]  = muon::isGoodMuon(muon, muon::TMOneStationAngLoose) ? 1:0;
		fTmuIsTMOneStationAngTight[mqi]  = muon::isGoodMuon(muon, muon::TMOneStationAngTight) ? 1:0;

		Ref<View<Muon> > muonRef(muons,index);
		const reco::IsoDeposit ECDep = ECDepMap[muonRef];
		const reco::IsoDeposit HCDep = HCDepMap[muonRef];
		fTmueecal[mqi] = ECDep.candEnergy();
		fTmuehcal[mqi] = HCDep.candEnergy();

		fTmud0bs[mqi] = -1.0*muon.innerTrack()->dxy(beamSpot.position());
		fTmud0pv[mqi] = -1.0*muon.innerTrack()->dxy(primVtx->position());
		fTmudzbs[mqi] = muon.innerTrack()->dz(beamSpot.position());
		fTmudzpv[mqi] = muon.innerTrack()->dz(primVtx->position());
		fTmuinntknchi2[mqi] = muon.innerTrack()->normalizedChi2();

		// Separate methods:
		if(fTmuIsTM[mqi]){ // Tracker Muons
			fTntrackermu++;
			fTmuptE[mqi]  = muon.innerTrack()->ptError();
			fTmud0E[mqi]  = muon.innerTrack()->dxyError();
			fTmudzE[mqi]  = muon.innerTrack()->dzError();

			fTmunchi2[mqi]      = fTmuinntknchi2[mqi]; // No difference for TM
			fTmunglhits[mqi]    = 0;
			fTmuntkhits[mqi]    = muon.innerTrack()->hitPattern().numberOfValidHits();
			fTmunpxhits[mqi]    = muon.innerTrack()->hitPattern().numberOfValidPixelHits();
			fTmunmuhits[mqi]    = 0;
			fTmunmatches[mqi]   = 0;
			fTmunchambers[mqi]  = 0;
		}
		if(fTmuIsGM[mqi]){ // Global Muons
			fTnglobalmu++;
			fTmuptE[mqi]  = muon.globalTrack()->ptError();
			fTmud0E[mqi]  = muon.globalTrack()->dxyError();
			fTmudzE[mqi]  = muon.globalTrack()->dzError();

			fTmunchi2[mqi]      = muon.globalTrack()->normalizedChi2();
			fTmunglhits[mqi]    = muon.globalTrack()->hitPattern().numberOfValidHits();
			fTmuntkhits[mqi]    = muon.innerTrack()->hitPattern().numberOfValidHits();
			fTmunpxhits[mqi]    = muon.innerTrack()->hitPattern().numberOfValidPixelHits();
			fTmunmuhits[mqi]    = muon.outerTrack()->hitPattern().numberOfValidHits();
			fTmunmatches[mqi]   = muon.numberOfMatches();
			fTmunchambers[mqi]  = muon.numberOfChambers();
		}

		// MC Matching
		if(!fIsRealData){
			std::vector<const GenParticle*> MuMatch = matchRecoCand(&muon, iEvent);
			if(MuMatch[0] != NULL){
				fTGenMuId[mqi]       = MuMatch[0]->pdgId();
				fTGenMuStatus[mqi]   = MuMatch[0]->status();
				fTGenMuPt[mqi]       = MuMatch[0]->pt();
				fTGenMuEta[mqi]      = MuMatch[0]->eta();
				fTGenMuPhi[mqi]      = MuMatch[0]->phi();
				fTGenMuE[mqi]        = MuMatch[0]->energy();

				fTGenMuMId[mqi]      = MuMatch[1]->pdgId();
				fTGenMuMStatus[mqi]  = MuMatch[1]->status();
				fTGenMuMPt[mqi]      = MuMatch[1]->pt();
				fTGenMuMEta[mqi]     = MuMatch[1]->eta();
				fTGenMuMPhi[mqi]     = MuMatch[1]->phi();
				fTGenMuME[mqi]       = MuMatch[1]->energy();

				fTGenMuGMId[mqi]     = MuMatch[2]->pdgId();
				fTGenMuGMStatus[mqi] = MuMatch[2]->status();
				fTGenMuGMPt[mqi]     = MuMatch[2]->pt();
				fTGenMuGMEta[mqi]    = MuMatch[2]->eta();
				fTGenMuGMPhi[mqi]    = MuMatch[2]->phi();
				fTGenMuGME[mqi]      = MuMatch[2]->energy();
			}
			MuMatch.clear();
		}
		fTgoodmu[mqi]  = 0;
		fTmuIsIso[mqi] = 1;
	}


	// SC variables
	fTnSC=0;
	for (SuperClusterCollection::const_iterator sc = BarrelSuperClusters->begin(); sc!=BarrelSuperClusters->end(); ++sc){

	  if (sc->rawEnergy()<fMinSCraw) continue;
	  if (sc->rawEnergy()/TMath::CosH(sc->eta())<fMinSCrawPt) continue;

	  if (fTnSC>=gMaxnSC) {
	    edm::LogWarning("NTP") << "@SUB=analyze" << "Maximum number of Super Clusters exceeded"; 
	    fTgoodevent = 1; 
	    break;
          }
	  
	  fTSCx[fTnSC] = sc->x();
	  fTSCy[fTnSC] = sc->y();
	  fTSCz[fTnSC] = sc->z();
	  fTSCraw[fTnSC] = sc->rawEnergy();
	  fTSCpre[fTnSC] = sc->preshowerEnergy();
	  fTSCenergy[fTnSC] = sc->energy();
	  fTSCeta[fTnSC] = sc->eta();
	  fTSCphi[fTnSC] = sc->phi();
	  fTSCsigmaPhi[fTnSC] = sc->phiWidth();
	  fTSCsigmaEta[fTnSC] = sc->etaWidth();
	  fTSCbrem[fTnSC] = (sc->etaWidth()!=0) ? sc->phiWidth()/sc->etaWidth() : -1;
	  fTSCR9[fTnSC] = sc->rawEnergy()!=0 ? EcalClusterTools::e3x3(  *(sc->seed()), ebRecHits.product(), &(*topology)) / sc->rawEnergy() : -1;
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
		fTSCcrackcorrseedfactor[fTnSC] = CrackCorrFunc->getValue(*cc);
		fTSClocalcorrseedfactor[fTnSC] = LocalCorrFunc->getValue(*cc);		
	      }
	      crackcorrenergy += (*itClus)->energy()*(CrackCorrFunc->getValue(*cc)-1);
	      localcorrenergy += (*itClus)->energy()*(LocalCorrFunc->getValue(*cc)-1);
	      index++;
	    }
	    fTSCcrackcorrseed[fTnSC] = crackcorrseedenergy/sc->rawEnergy();
	    fTSCcrackcorr[fTnSC] = crackcorrenergy/sc->rawEnergy();
	    fTSClocalcorrseed[fTnSC] = localcorrseedenergy/sc->rawEnergy();
	    fTSClocalcorr[fTnSC] = localcorrenergy/sc->rawEnergy();
	  }

	  {
	    std::vector<DetId> cristalli;	  
	    for (reco::CaloCluster_iterator bc=sc->clustersBegin(); bc!=sc->clustersEnd(); ++bc){
	      const std::vector< std::pair<DetId, float> > & seedrechits = (*bc)->hitsAndFractions();
	      for (uint i=0; i<seedrechits.size(); i++) cristalli.push_back(seedrechits[i].first);
	      sort(cristalli.begin(),cristalli.end());
	      std::vector<DetId>::iterator it;
	      it = unique(cristalli.begin(),cristalli.end());
	      cristalli.resize(it-cristalli.begin());
	    }
	  
	    uint i=0;
	    for (i=0; i<cristalli.size(); i++){

	      if ((int)i>=gMaxnSCxtals){
		edm::LogWarning("NTP") << "@SUB=analyze" << "Maximum number of SC xtals exceeded!";
		fTgoodevent = 1;
		break;
	      }

	      CaloCellGeometry *cellGeometry = NULL;
	      if (cristalli.at(i).subdetId()!=EcalBarrel) {
		edm::LogWarning("NTP") << "@SUB=analyze" << "Problem with xtals subdetId()";
		continue;
	      } 
	      EBDetId ebDetId  = cristalli.at(i);
	      cellGeometry = (CaloCellGeometry*)(barrelGeometry->getGeometry(ebDetId));
	      TVector3 xtal_position(cellGeometry->getPosition().x(),cellGeometry->getPosition().y(),cellGeometry->getPosition().z());
	      float dphi=(dynamic_cast<const EcalBarrelGeometry*>(barrelGeometry))->deltaPhi(ebDetId);
	      float deta=(dynamic_cast<const EcalBarrelGeometry*>(barrelGeometry))->deltaEta(ebDetId);
	      fTSCxtalX[fTnSC][i]=xtal_position.x();
	      fTSCxtalY[fTnSC][i]=xtal_position.y();
	      fTSCxtalZ[fTnSC][i]=xtal_position.z();
	      fTSCxtalEtaWidth[fTnSC][i]=deta;
	      fTSCxtalPhiWidth[fTnSC][i]=dphi;
	      const CaloCellGeometry::CornersVec& cellCorners (cellGeometry->getCorners());
	      for (int k=0; k<4; k++){
		fTSCxtalfrontX[fTnSC][i][k]=(float)(cellCorners[k].x());
		fTSCxtalfrontY[fTnSC][i][k]=(float)(cellCorners[k].y());
		fTSCxtalfrontZ[fTnSC][i][k]=(float)(cellCorners[k].z());
	      }
	    }
	    fTSCNXtals[fTnSC]=i;
	  }

	  fTnSC++;
	}

	for (SuperClusterCollection::const_iterator sc = EndcapSuperClusters->begin(); sc!=EndcapSuperClusters->end(); ++sc){

	  if (sc->rawEnergy()<fMinSCraw) continue;
	  if (sc->rawEnergy()/TMath::CosH(sc->eta())<fMinSCrawPt) continue;

	  if (fTnSC>=gMaxnSC) {
	    edm::LogWarning("NTP") << "@SUB=analyze" << "Maximum number of Super Clusters exceeded"; 
	    fTgoodevent = 1; 
	    break;
          }

	  fTSCx[fTnSC] = sc->x();
	  fTSCy[fTnSC] = sc->y();
	  fTSCz[fTnSC] = sc->z();
	  fTSCraw[fTnSC] = sc->rawEnergy();
	  fTSCpre[fTnSC] = sc->preshowerEnergy();
	  fTSCenergy[fTnSC] = sc->energy();
	  fTSCeta[fTnSC] = sc->eta();
	  fTSCphi[fTnSC] = sc->phi();
	  fTSCsigmaPhi[fTnSC] = sc->phiWidth();
	  fTSCsigmaEta[fTnSC] = sc->etaWidth();
	  fTSCbrem[fTnSC] = (sc->etaWidth()!=0) ? sc->phiWidth()/sc->etaWidth() : -1;
	  fTSCR9[fTnSC] = sc->rawEnergy()!=0 ? EcalClusterTools::e3x3(  *(sc->seed()), eeRecHits.product(), &(*topology)) / sc->rawEnergy() : -1;
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
		fTSCcrackcorrseedfactor[fTnSC] = CrackCorrFunc->getValue(*cc);
		fTSClocalcorrseedfactor[fTnSC] = LocalCorrFunc->getValue(*cc);		
	      }
	      crackcorrenergy += (*itClus)->energy()*(CrackCorrFunc->getValue(*cc)-1);
	      localcorrenergy += (*itClus)->energy()*(LocalCorrFunc->getValue(*cc)-1);
	      index++;
	    }
	    fTSCcrackcorrseed[fTnSC] = crackcorrseedenergy/sc->rawEnergy();
	    fTSCcrackcorr[fTnSC] = crackcorrenergy/sc->rawEnergy();
	    fTSClocalcorrseed[fTnSC] = localcorrseedenergy/sc->rawEnergy();
	    fTSClocalcorr[fTnSC] = localcorrenergy/sc->rawEnergy();
	  }

	  {
	    std::vector<DetId> cristalli;	  
	    for (reco::CaloCluster_iterator bc=sc->clustersBegin(); bc!=sc->clustersEnd(); ++bc){
	      const std::vector< std::pair<DetId, float> > & seedrechits = (*bc)->hitsAndFractions();
	      for (uint i=0; i<seedrechits.size(); i++) cristalli.push_back(seedrechits[i].first);
	      sort(cristalli.begin(),cristalli.end());
	      std::vector<DetId>::iterator it;
	      it = unique(cristalli.begin(),cristalli.end());
	      cristalli.resize(it-cristalli.begin());
	    }
	  
	    uint i=0;
	    for (i=0; i<cristalli.size(); i++){

	      if ((int)i>=gMaxnSCxtals){
		edm::LogWarning("NTP") << "@SUB=analyze" << "Maximum number of SC xtals exceeded!";
		fTgoodevent = 1;
		break;
	      }

	      CaloCellGeometry *cellGeometry = NULL;
	      if (cristalli.at(i).subdetId()!=EcalEndcap) {
		edm::LogWarning("NTP") << "@SUB=analyze" << "Problem with xtals subdetId()";
		continue;
	      } 
	      EEDetId eeDetId  = cristalli.at(i);
	      cellGeometry = (CaloCellGeometry*)(endcapGeometry->getGeometry(eeDetId));
	      TVector3 xtal_position(cellGeometry->getPosition().x(),cellGeometry->getPosition().y(),cellGeometry->getPosition().z());
	      float dphi=(dynamic_cast<const EcalEndcapGeometry*>(endcapGeometry))->deltaPhi(eeDetId);
	      float deta=(dynamic_cast<const EcalEndcapGeometry*>(endcapGeometry))->deltaEta(eeDetId);
	      fTSCxtalX[fTnSC][i]=xtal_position.x();
	      fTSCxtalY[fTnSC][i]=xtal_position.y();
	      fTSCxtalZ[fTnSC][i]=xtal_position.z();
	      fTSCxtalEtaWidth[fTnSC][i]=deta;
	      fTSCxtalPhiWidth[fTnSC][i]=dphi;
	      const CaloCellGeometry::CornersVec& cellCorners (cellGeometry->getCorners());
	      for (int k=0; k<4; k++){
		fTSCxtalfrontX[fTnSC][i][k]=(float)(cellCorners[k].x());
		fTSCxtalfrontY[fTnSC][i][k]=(float)(cellCorners[k].y());
		fTSCxtalfrontZ[fTnSC][i][k]=(float)(cellCorners[k].z());
	      }
	    }
	    fTSCNXtals[fTnSC]=i;
	  }

	  fTnSC++;
	}




	////////////////////////////////////////////////////////
	// Electron variables:
	// Keep pointers to electron superCluster in original collections
	std::vector<const SuperCluster*> elecPtr;
	std::vector<const GsfTrack*> trckPtr;
	int eqi(0);                    // Index of qualified electrons
	fTnelestot = electrons->size(); // Total number of electrons

	if(electrons->size() > 0){
		// Get electrons, order them by pt and apply selection
		std::vector<OrderPair> elOrdered;
		int elIndex(0);
		for( View<GsfElectron>::const_iterator El = electrons->begin();
		El != electrons->end(); ++El, ++elIndex ) {
		// Check if maximum number of electrons is exceeded already:
			if(eqi >= gMaxneles) {
				edm::LogWarning("NTP") << "@SUB=analyze"
					<< "Maximum number of electrons exceeded..";
				fTflagmaxelexc = 1;
				fTgoodevent = 1;
				break;
			}
		// Electron preselection:
			if(El->pt() < fMinelpt) continue;
			if(fabs(El->eta()) > fMaxeleta) continue;

			eqi++; // Count how many we'll eventually store
			elOrdered.push_back(make_pair(elIndex,El->pt()));
		}
		std::sort(elOrdered.begin(),elOrdered.end(),indexComparator);
		fTneles = elOrdered.size();
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

			fTepx[eqi]                     = electron.px();
			fTepy[eqi]                     = electron.py();
			fTepz[eqi]                     = electron.pz();
			fTept[eqi]                     = electron.pt();
			fTeptE[eqi]                    = electron.gsfTrack()->ptError();
			fTeeta[eqi]                    = electron.eta();
			fTephi[eqi]                    = electron.phi();
			fTegsfpt [eqi]                 = electron.gsfTrack()->pt();
			fTegsfeta[eqi]                 = electron.gsfTrack()->eta();
			fTegsfphi[eqi]                 = electron.gsfTrack()->phi();
			fTetrkmomerror[eqi]            = electron.trackMomentumError();
			fTeecalergerror[eqi]           = electron.ecalEnergyError();
                        // 4_2: take the error on the default momentum
			fTeelemomerror[eqi]            = electron.p4Error( electron.candidateP4Kind() );
			fTenbrems[eqi]                 = electron.numberOfBrems();
			fTee[eqi]                      = electron.energy();
			fTeet[eqi]                     = electron.et();
			fTed0bs[eqi]                   = -1.0*electron.gsfTrack()->dxy(beamSpot.position());
			fTed0pv[eqi]                   = -1.0*electron.gsfTrack()->dxy(primVtx->position());
			fTed0E[eqi]                    = electron.gsfTrack()->dxyError();
			fTedzbs[eqi]                   = electron.gsfTrack()->dz(beamSpot.position());
			fTedzpv[eqi]                   = electron.gsfTrack()->dz(primVtx->position());
			fTedzE[eqi]                    = electron.gsfTrack()->dzError();
			fTenchi2[eqi]                  = electron.gsfTrack()->normalizedChi2();
			fTdr03tksumpt[eqi]             = electron.dr03TkSumPt();
			fTdr03ecalrechitsumet[eqi]     = electron.dr03EcalRecHitSumEt();
			fTdr03hcaltowersumet[eqi]      = electron.dr03HcalTowerSumEt();
			fTdr04tksumpt[eqi]             = electron.dr04TkSumPt();
			fTdr04ecalrechitsumet[eqi]     = electron.dr04EcalRecHitSumEt();
			fTdr04hcaltowersumet[eqi]      = electron.dr04HcalTowerSumEt();
			fTeiso03[eqi]                  = (fTdr03tksumpt[eqi] + fTdr03ecalrechitsumet[eqi] + fTdr03hcaltowersumet[eqi]) / fTept[eqi];
			fTeiso04[eqi]                  = (fTdr04tksumpt[eqi] + fTdr04ecalrechitsumet[eqi] + fTdr04hcaltowersumet[eqi]) / fTept[eqi];
			fTecharge[eqi]                 = electron.charge();
			fTeInGap[eqi]                  = electron.isGap() ? 1:0;
			fTeEcalDriven[eqi]             = electron.ecalDrivenSeed() ? 1:0;
			fTeTrackerDriven[eqi]          = electron.trackerDrivenSeed() ? 1:0;
			fTeCInfoIsGsfCtfCons[eqi]      = electron.chargeInfo().isGsfCtfConsistent ? 1:0;
			fTeCInfoIsGsfCtfScPixCons[eqi] = electron.chargeInfo().isGsfCtfScPixConsistent ? 1:0;
			fTeCInfoIsGsfScPixCons[eqi]    = electron.chargeInfo().isGsfScPixConsistent ? 1:0;
			fTeCInfoScPixCharge[eqi]       = electron.chargeInfo().scPixCharge;
			if( electron.closestCtfTrackRef().isNonnull() ){
				fTeClosestCtfTrackpt[eqi]      = (electron.closestCtfTrack().ctfTrack)->pt();
				fTeClosestCtfTracketa[eqi]     = (electron.closestCtfTrack().ctfTrack)->eta();
				fTeClosestCtfTrackphi[eqi]     = (electron.closestCtfTrack().ctfTrack)->phi();
				fTeClosestCtfTrackcharge[eqi]  = (electron.closestCtfTrack().ctfTrack)->charge();
			}

			fTeBasicClustersSize[eqi]         = electron.basicClustersSize();
			fTefbrem[eqi]                     = electron.fbrem();
			fTeHcalOverEcal[eqi]              = electron.hcalOverEcal();
			fTeE1x5[eqi]                      = electron.e1x5();
			fTeE5x5[eqi]                      = electron.e5x5();
			fTeE2x5Max[eqi]                   = electron.e2x5Max();
			fTeSigmaIetaIeta[eqi]             = electron.sigmaIetaIeta();
			fTeDeltaEtaSeedClusterAtCalo[eqi] = electron.deltaEtaSeedClusterTrackAtCalo();
			fTeDeltaPhiSeedClusterAtCalo[eqi] = electron.deltaPhiSeedClusterTrackAtCalo();
			fTeDeltaPhiSuperClusterAtVtx[eqi] = electron.deltaPhiSuperClusterTrackAtVtx();
			fTeDeltaEtaSuperClusterAtVtx[eqi] = electron.deltaEtaSuperClusterTrackAtVtx();
			fTecaloenergy[eqi]                = electron.caloEnergy();
			fTetrkmomatvtx[eqi]               = electron.trackMomentumAtVtx().R();
			fTeESuperClusterOverP[eqi]        = electron.eSuperClusterOverP();
			fTeNumberOfMissingInnerHits[eqi]  = electron.gsfTrack()->trackerExpectedHitsInner().numberOfHits();
			fTetheta[eqi]                     = electron.superCluster()->position().theta();
			fTesceta[eqi]                     = electron.superCluster()->eta();



// DISABLED: NO SEED IN AOD (UPDATE IT IN 4_2)
// 			if ( electron.superCluster()->seed()->caloID().detector( reco::CaloID::DET_ECAL_BARREL ) ) {
//                           fTeScSeedSeverity[eqi] = EcalSeverityLevelAlgo::severityLevel( electron.superCluster()->seed()->seed(), *ebRecHits, *channelStatus );
//                           fTeE1OverE9[eqi]       = EcalSeverityLevelAlgo::E1OverE9(   electron.superCluster()->seed()->seed(), *ebRecHits );
//                           fTeS4OverS1[eqi]       = EcalSeverityLevelAlgo::swissCross( electron.superCluster()->seed()->seed(), *ebRecHits );
// 			} else if ( electron.superCluster()->seed()->caloID().detector( reco::CaloID::DET_ECAL_ENDCAP ) ) {
//                           fTeScSeedSeverity[eqi] = EcalSeverityLevelAlgo::severityLevel( electron.superCluster()->seed()->seed(), *eeRecHits, *channelStatus );
//                           fTeE1OverE9[eqi]       = EcalSeverityLevelAlgo::E1OverE9(   electron.superCluster()->seed()->seed(), *eeRecHits );
//                           fTeS4OverS1[eqi]       = EcalSeverityLevelAlgo::swissCross( electron.superCluster()->seed()->seed(), *eeRecHits );
// 			} else {
// 				edm::LogWarning("NTP") << "Electron supercluster seed crystal neither in EB nor in EE!";
// 			}

			// Read in Electron ID
			fTeIDMva[eqi] = electron.mva();
			Ref<View<GsfElectron> > electronRef(electrons,index);
			fTeIDTight[eqi]            = eIDmapT[electronRef]  ? 1:0;
			fTeIDLoose[eqi]            = eIDmapL[electronRef]  ? 1:0;
			fTeIDRobustTight[eqi]      = eIDmapRT[electronRef] ? 1:0;
			fTeIDRobustLoose[eqi]      = eIDmapRL[electronRef] ? 1:0;
			fTeIDsimpleWPrelIso[eqi]   = eIDmapsimpleWP[electronRef];
			fTeIDsimpleWP95relIso[eqi] = eIDmapsimpleWP95[electronRef];
			fTeIDsimpleWP90relIso[eqi] = eIDmapsimpleWP90[electronRef];
			fTeIDsimpleWP85relIso[eqi] = eIDmapsimpleWP85[electronRef];
			fTeIDsimpleWP80relIso[eqi] = eIDmapsimpleWP80[electronRef];

			{
			  fTElSCindex[eqi] = -1;
			  float diff=1e+4;
			  for (int scind=0; scind<fTnSC; scind++){
			        if (fabs(fTSCraw[scind]-electron.superCluster()->rawEnergy())<diff) {
				  fTElSCindex[eqi]=scind;
				  diff=fabs(fTSCraw[scind]-electron.superCluster()->rawEnergy());
				}
			    }

			  if (fTElSCindex[eqi]!=-1)
			  if (fabs(fTSCeta[fTElSCindex[eqi]]-electron.superCluster()->eta())>0.1 || \
			      fabs(reco::deltaPhi(fTSCphi[fTElSCindex[eqi]],electron.superCluster()->phi()))>0.1){
			    // 			    std::cout << fTSCraw[fTElSCindex[eqi]] << " " << electron.superCluster()->rawEnergy() << std::endl;
// 			    			    std::cout << fTSCeta[fTElSCindex[eqi]] << " " << electron.superCluster()->eta() << std::endl;
// 			    			    std::cout << fTSCphi[fTElSCindex[eqi]] << " " << electron.superCluster()->phi() << std::endl;
			    fTElSCindex[eqi] = -1;			    
			  }

			  if (fTElSCindex[eqi]==-1) {
			    //edm::LogWarning("NTP") << "@SUB=analyze" << "No matching SC found for electron"; 
			    //			    fTgoodevent = 1; 
			    //			    break;
			  }
			}


			// MC Matching
			if(!fIsRealData){
				std::vector<const GenParticle*> ElMatch = matchRecoCand(&electron, iEvent);
				if(ElMatch[0] != NULL){
					fTGenElId[eqi]       = ElMatch[0]->pdgId();
					fTGenElStatus[eqi]   = ElMatch[0]->status();
					fTGenElPt[eqi]       = ElMatch[0]->pt();
					fTGenElEta[eqi]      = ElMatch[0]->eta();
					fTGenElPhi[eqi]      = ElMatch[0]->phi();
					fTGenElE[eqi]        = ElMatch[0]->energy();

					fTGenElMId[eqi]      = ElMatch[1]->pdgId();
					fTGenElMStatus[eqi]  = ElMatch[1]->status();
					fTGenElMPt[eqi]      = ElMatch[1]->pt();
					fTGenElMEta[eqi]     = ElMatch[1]->eta();
					fTGenElMPhi[eqi]     = ElMatch[1]->phi();
					fTGenElME[eqi]       = ElMatch[1]->energy();

					fTGenElGMId[eqi]     = ElMatch[2]->pdgId();
					fTGenElGMStatus[eqi] = ElMatch[2]->status();
					fTGenElGMPt[eqi]     = ElMatch[2]->pt();
					fTGenElGMEta[eqi]    = ElMatch[2]->eta();
					fTGenElGMPhi[eqi]    = ElMatch[2]->phi();
					fTGenElGME[eqi]      = ElMatch[2]->energy();
				}
				ElMatch.clear();
			}


			// Conversion Information
			reco::GsfElectron::ConversionRejection ConvRejVars = electron.conversionRejectionVariables();
			reco::TrackBaseRef ConvPartnerTrack = ConvRejVars.partner;
			if( ConvPartnerTrack.isNonnull() ){
				fTeConvPartTrackDist[eqi]   = ConvRejVars.dist;
				fTeConvPartTrackDCot[eqi]   = ConvRejVars.dcot;
				fTeConvPartTrackPt[eqi]     = ConvPartnerTrack->pt();
				fTeConvPartTrackEta[eqi]    = ConvPartnerTrack->eta();
				fTeConvPartTrackPhi[eqi]    = ConvPartnerTrack->phi();
				fTeConvPartTrackCharge[eqi] = ConvPartnerTrack->charge();
			}


			// fTeIsInJet[eqi] = -1;
			// fTeSharedPx[eqi] = 0.;
			// fTeSharedPy[eqi] = 0.;
			// fTeSharedPz[eqi] = 0.;
			// fTeSharedEnergy[eqi] = 0.;

			fTgoodel[eqi] = 0;
			fTeIsIso[eqi] = 1;
			fTeChargeMisIDProb[eqi] = 0;
			// fTeDupEl[eqi] = -1;
		}
	}

	fTnEBhits=0;
	for(EcalRecHitCollection::const_iterator ecalrechit = ebRecHits->begin(); ecalrechit!=ebRecHits->end() ; ++ecalrechit)
	{
		double energy = ecalrechit->energy();
		if(energy<fMinebrechitE)continue; 
		if(fTnEBhits>=gMaxnEBhits)
		{
                  edm::LogWarning("NTP") << "@SUB=analyze" << "Maximum number of EB rechits exceeded"; 
                  fTgoodevent = 1; 
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

		fTEBrechitE[fTnEBhits] =  hitPos.Mag();
		fTEBrechitPt[fTnEBhits] =  hitPos.Pt();
		fTEBrechitEta[fTnEBhits] =  hitPos.Eta();
		fTEBrechitPhi[fTnEBhits] = hitPos.Phi() ;
		fTEBrechitTime[fTnEBhits] = time;
		fTEBrechitChi2[fTnEBhits] = chi2;
//4_2 		fTEBrechitE4oE1[fTnEBhits] =e4oe1 ;
//4_2 		fTEBrechitE2oE9[fTnEBhits] = e2oe9;

	//	cout << "ebrechit P =" << fTEBrechitE[fTnEBhits] << " Pt = " << fTEBrechitPt[fTnEBhits] ;
	//	cout << " eta = " << fTEBrechitEta[fTnEBhits] << "phi = " << fTEBrechitPhi[fTnEBhits] << " time = "<< fTEBrechitTime[fTnEBhits] ;
	//	cout << " chi2 = "<< fTEBrechitChi2[fTnEBhits] << " e4oe1 = "<< fTEBrechitE4oE1[fTnEBhits] << " e2oe9 = " << e2oe9 << endl;

		fTnEBhits++;
	}








	////////////////////////////////////////////////////////
	// Photon Variables:
	// Keep pointers to superclusters for cross cleaning
	std::vector<const SuperCluster*> photSCs;
	int phoqi(0); // Index of qualified photons
	fTnphotonstot = photons->size();



	// Get photonss, order them by pt and apply selection
	std::vector<OrderPair> phoOrdered;
	int phoIndex(0);
	for( View<Photon>::const_iterator ip = photons->begin();
	ip != photons->end(); ++ip, ++phoIndex ){
		// Check if maximum number of photons exceeded
		if(phoqi >= gMaxnphos){
			edm::LogWarning("NTP") << "@SUB=analyze"
				<< "Maximum number of photons exceeded";
			fTflagmaxphoexc = 1;
			fTgoodevent = 1;
			break;
		}
		// Preselection
		if(ip->superCluster()->rawEnergy()/TMath::CosH(ip->superCluster()->eta()) < fMinphopt) continue;
		if(fabs(ip->eta()) > fMaxphoeta) continue;

		phoqi++; // Count how many we'll eventually store
		phoOrdered.push_back(make_pair(phoIndex,ip->pt()));
	}
	std::sort(phoOrdered.begin(),phoOrdered.end(),indexComparator);
	fTnphotons = phoOrdered.size();
	phoqi = 0;

	std::vector<int> storethispfcand(pfCandidates->size(),0);
	std::vector<int> PhotonToPFPhotonMatchingArray(gMaxnphos,-999);
	std::vector<int> PhotonToPFElectronMatchingArray(gMaxnphos,-999);

	for (std::vector<OrderPair>::const_iterator it = phoOrdered.begin();
	it != phoOrdered.end(); ++it, ++phoqi ) {

	int index = it->first;
	const Photon& photon = (*photons)[index];

	// Save photon supercluster position
	photSCs.push_back(&(*photon.superCluster()));

	fTPhotPt[phoqi]             = photon.pt();
	fTPhotPx[phoqi]             = photon.px();
	fTPhotPy[phoqi]             = photon.py();
	fTPhotPz[phoqi]             = photon.pz();
	fTPhotEta[phoqi]            = photon.eta();
	fTPhotPhi[phoqi]            = photon.phi();
	fTPhotEnergy[phoqi]         = photon.energy();
	fTPhotIso03Ecal[phoqi]      = photon.ecalRecHitSumEtConeDR03();
	fTPhotIso03Hcal[phoqi]      = photon.hcalTowerSumEtConeDR03();
	fTPhotIso03TrkSolid[phoqi]  = photon.trkSumPtSolidConeDR03();
	fTPhotIso03TrkHollow[phoqi] = photon.trkSumPtHollowConeDR03();
	fTPhotIso03[phoqi]          = (fTPhotIso03TrkHollow[phoqi] + fTPhotIso03Ecal[phoqi] + fTPhotIso03Hcal[phoqi]) / fTPhotPt[phoqi];
	fTPhotIso04Ecal[phoqi]      = photon.ecalRecHitSumEtConeDR04();
	fTPhotIso04Hcal[phoqi]      = photon.hcalTowerSumEtConeDR04();
	fTPhotIso04TrkSolid[phoqi]  = photon.trkSumPtSolidConeDR04();
	fTPhotIso04TrkHollow[phoqi] = photon.trkSumPtHollowConeDR04();
	fTPhotIso04[phoqi]          = (fTPhotIso04TrkHollow[phoqi] + fTPhotIso04Ecal[phoqi] + fTPhotIso04Hcal[phoqi]) / fTPhotPt[phoqi];
	fTPhotR9[phoqi]             = photon.r9();
	fTPhotcaloPosX[phoqi]       = photon.caloPosition().X();
	fTPhotcaloPosY[phoqi]       = photon.caloPosition().Y();
	fTPhotcaloPosZ[phoqi]       = photon.caloPosition().Z();
	fTPhotHoverE[phoqi]         = photon.hadronicOverEm();
	fTPhotH1overE[phoqi]        = photon.hadronicDepth1OverEm();
	fTPhotH2overE[phoqi]        = photon.hadronicDepth2OverEm();
	fTPhotSigmaIetaIeta[phoqi]  = photon.sigmaIetaIeta();

	fTPhotVx[phoqi] = photon.vx();
	fTPhotVy[phoqi] = photon.vy();
	fTPhotVz[phoqi] = photon.vz();

       fTPhotSigmaEtaEta[phoqi] = photon.sigmaEtaEta();
       fTPhote1x5[phoqi]= photon.e1x5();
       fTPhote2x5[phoqi]= photon.e2x5();
       fTPhote3x3[phoqi]= photon.e3x3();
       fTPhote5x5[phoqi]= photon.e5x5();
       fTPhotmaxEnergyXtal[phoqi]= photon.maxEnergyXtal();
       fTPhotIso03HcalDepth1[phoqi]= photon.hcalDepth1TowerSumEtConeDR03();
       fTPhotIso03HcalDepth2[phoqi]= photon.hcalDepth2TowerSumEtConeDR03();
       fTPhotIso04HcalDepth1[phoqi]= photon.hcalDepth1TowerSumEtConeDR04();
       fTPhotIso04HcalDepth2[phoqi]= photon.hcalDepth2TowerSumEtConeDR04();
       fTPhotIso03nTrksSolid[phoqi]= photon.nTrkSolidConeDR03();
       fTPhotIso03nTrksHollow[phoqi]= photon.nTrkHollowConeDR03();
       fTPhotIso04nTrksSolid[phoqi]= photon.nTrkSolidConeDR04();
       fTPhotIso04nTrksHollow[phoqi]= photon.nTrkHollowConeDR04();

	// EcalClusterLazyTools *lazyTools = new EcalClusterLazyTools(iEvent, iSetup, edm::InputTag("reducedEcalRecHitsEB"), edm::InputTag("reducedEcalRecHitsEE"));

	fTPhotSCEnergy[phoqi]       = photon.superCluster()->rawEnergy();
	fTPhotSCEtaWidth[phoqi]     = photon.superCluster()->etaWidth();
// DISABLED: NO SEED IN AOD (UPDATE IT IN 4_2)
// 	fTPhotSCSigmaPhiPhi[phoqi]  = lazyTools->covariances(*(photon.superCluster()->seed())).at(2);
	fTPhotHasPixSeed[phoqi]     = photon.hasPixelSeed() ? 1:0;
	fTPhotPassConvSafeElectronVeto[phoqi] = !ConversionTools::hasMatchedPromptElectron(photon.superCluster(), electronHandle, convH, beamspot.position());
	fTPhotHasConvTrks[phoqi]    = photon.hasConversionTracks() ? 1:0;

	// fTPhotIsInJet[phoqi]      = -1;
	// fTPhotDupEl[phoqi]        = -1;
	// fTPhotSharedPx[phoqi]     = 0.;
	// fTPhotSharedPy[phoqi]     = 0.;
	// fTPhotSharedPz[phoqi]     = 0.;
	// fTPhotSharedEnergy[phoqi] = 0.;
	fTgoodphoton[phoqi]       = 0;
	fTPhotIsIso[phoqi]        = 1;

       fTPhotisEB[phoqi]= photon.isEB();
       fTPhotisEE[phoqi]= photon.isEE();
       fTPhotisEBEtaGap[phoqi]= photon.isEBEtaGap();
       fTPhotisEBPhiGap[phoqi]= photon.isEBPhiGap();
       fTPhotisEERingGap[phoqi]= photon.isEERingGap();
       fTPhotisEEDeeGap[phoqi]= photon.isEEDeeGap();
       fTPhotisEBEEGap[phoqi]= photon.isEBEEGap();
       fTPhotisPFlowPhoton[phoqi]= photon.isPFlowPhoton();
       fTPhotisStandardPhoton[phoqi]= photon.isStandardPhoton();

       {
	 SuperClusterFootprintRemoval remover(iEvent,iSetup);
	 fTPhoSCRemovalPFIsoCharged[phoqi] = (remover.PFIsolation("charged",photon.superCluster(),-1));
	 fTPhoSCRemovalPFIsoChargedPrimVtx[phoqi] = (vertices->size()>0) ? (remover.PFIsolation("charged",photon.superCluster(),0)) : fTPhoSCRemovalPFIsoCharged[phoqi];
	 fTPhoSCRemovalPFIsoNeutral[phoqi] = (remover.PFIsolation("neutral",photon.superCluster()));
	 fTPhoSCRemovalPFIsoPhoton[phoqi] = (remover.PFIsolation("photon",photon.superCluster()));
	 PFIsolation_RandomCone_struct risos = remover.RandomConeIsolation(photon.superCluster(),-1);
	 fTPhoSCRemovalPFIsoCharged_RCone[phoqi] = risos.chargediso;
	 fTPhoSCRemovalPFIsoChargedPrimVtx_RCone[phoqi] = risos.chargediso_primvtx;
	 fTPhoSCRemovalPFIsoNeutral_RCone[phoqi] = risos.neutraliso;
	 fTPhoSCRemovalPFIsoPhoton_RCone[phoqi] = risos.photoniso;
	 fTPhoSCRemoval_RCone_Eta[phoqi]=risos.randomcone_eta;
	 fTPhoSCRemoval_RCone_Phi[phoqi]=risos.randomcone_phi;
       }

       if (doVertexingFlag && photon.hasConversionTracks()) { // photon conversions

	 reco::ConversionRefVector conversions = photon.conversions();
	 if (conversions.size()<1) { std::cout << "something wrong here" << std::endl; }
	 reco::ConversionRef conv = conversions[0];
	 pho_conv_validvtx[phoqi]=conv->conversionVertex().isValid();

	 if(pho_conv_validvtx[phoqi] && ConversionsCut(*conv)) {
	   reco::Vertex vtx=conv->conversionVertex();
	   pho_conv_vtx[phoqi].SetXYZ(vtx.x(), vtx.y(), vtx.z());
	   pho_conv_chi2_probability[phoqi]=ChiSquaredProbability(vtx.chi2(), vtx.ndof());
	   pho_conv_ntracks[phoqi]=conv->nTracks();
	   pho_conv_eoverp[phoqi]=conv->EoverPrefittedTracks();
	   pho_conv_refitted_momentum[phoqi].SetXYZ(conv->refittedPairMomentum().x(), conv->refittedPairMomentum().y(), conv->refittedPairMomentum().z());
	 }

       }



       if (!fIsRealData){

         edm::Handle<GenParticleCollection> gen;
         iEvent.getByLabel(fGenPartTag, gen);
         GenParticleCollection::const_iterator g_part;

         std::vector<const reco::GenParticle*> matched = matchRecoCand(&photon,iEvent);

         if (matched[0]==NULL) {
	   fTPhotMCmatchexitcode[phoqi]=-1;
	   fTPhotMCmatchindex[phoqi]=-999;
	 }
         else if (matched[0]->pdgId()!=22) {
	   fTPhotMCmatchexitcode[phoqi]=0;
	   fTPhotMCmatchindex[phoqi]=-999;
	 }
         else {

           fTPhotMCmatchindex[phoqi]=-999;

           for(int i=0; i<fTngenphotons; ++i){
             if ( (fabs(fTGenPhotonPt[i]-matched[0]->pt())<0.01*matched[0]->pt()) \
                  && (fabs(fTGenPhotonEta[i]-matched[0]->eta())<0.01) \
                  && ( fabs(reco::deltaPhi(fTGenPhotonPhi[i],matched[0]->phi()))<0.01 ) ) {
               fTPhotMCmatchindex[phoqi]=i;
             }
           }

           if (fTPhotMCmatchindex[phoqi]!=-999){
	     if (fTGenPhotonMotherID[fTPhotMCmatchindex[phoqi]]>=-6 && fTGenPhotonMotherID[fTPhotMCmatchindex[phoqi]]<=6) fTPhotMCmatchexitcode[phoqi]=1;
	     else if (fTGenPhotonMotherID[fTPhotMCmatchindex[phoqi]]==21) fTPhotMCmatchexitcode[phoqi]=1;
             else if (fTGenPhotonMotherID[fTPhotMCmatchindex[phoqi]]==22 && fTGenPhotonMotherStatus[fTPhotMCmatchindex[phoqi]]==3) fTPhotMCmatchexitcode[phoqi]=2;
             else fTPhotMCmatchexitcode[phoqi]=3;
           }
	   else fTPhotMCmatchexitcode[phoqi]=-2;

         }

       }

       {
	 fTPhotSCindex[phoqi] = -1;
	 float diff=1e+4;
	 for (int scind=0; scind<fTnSC; scind++){
	   if (fabs(fTSCraw[scind]-photon.superCluster()->rawEnergy())<diff) {
	     fTPhotSCindex[phoqi]=scind;
	     diff=fabs(fTSCraw[scind]-photon.superCluster()->rawEnergy());
	   }
	 }

	 if (fTPhotSCindex[phoqi]!=-1)
	   if (fabs(fTSCeta[fTPhotSCindex[phoqi]]-photon.superCluster()->eta())>0.1 || \
	       fabs(reco::deltaPhi(fTSCphi[fTPhotSCindex[phoqi]],photon.superCluster()->phi()))>0.1){
	     //			    std::cout << fTSCraw[fTPhotSCindex[phoqi]] << " " << photon.superCluster()->rawEnergy() << std::endl;
	     //			    std::cout << fTSCeta[fTPhotSCindex[phoqi]] << " " << photon.superCluster()->eta() << std::endl;
	     //			    std::cout << fTSCphi[fTPhotSCindex[phoqi]] << " " << photon.superCluster()->phi() << std::endl;
	     fTPhotSCindex[phoqi] = -1;			    
	   }

	 if (fTPhotSCindex[phoqi]==-1) {
	   edm::LogWarning("NTP") << "@SUB=analyze" << "No matching SC found for photon"; 
	   cout << photon.superCluster()->rawEnergy() << " " << photon.superCluster()->eta() << " " << photon.superCluster()->rawEnergy()/TMath::CosH(photon.superCluster()->eta()) << endl;
	   //			    fTgoodevent = 1; 
	   //			    break;
	 }
       }


       { //Look for associated PF objects
	 
	 bool FoundPFPhoton=false;
	 bool FoundPFElectron=false;

	 fT_pho_isPFPhoton[phoqi] = 0;
	 fT_pho_isPFElectron[phoqi] = 0;

	 //Find PFPhoton
	 int iphot=-1;
	 int ncand=pfCandidates->size();
	 for( int i=0; i<ncand; ++i ) {
	   if ((*pfCandidates)[i].particleId()==reco::PFCandidate::gamma){
	     if ((*pfCandidates)[i].mva_nothing_gamma()>0){
	       if( (*pfCandidates)[i].superClusterRef()==photon.superCluster()) {
		 iphot = i;
	       }
	     }
	   }
	 }    
	 if (iphot!=-1) {
	   FoundPFPhoton=true;
	   fT_pho_isPFPhoton[phoqi] = 1;
	   PhotonToPFPhotonMatchingArray[phoqi] = iphot;
	 }

	 //Find PFElectron
	 bool foundEgSC = false;
	 int iel = -1;
	 reco::GsfElectronCollection::const_iterator elIterSl;
	 for (reco::GsfElectronCollection::const_iterator elIter = electronHandle->begin(); elIter != electronHandle->end(); ++elIter){
	   if (photon.superCluster()==elIter->superCluster()) {
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
	 if (iel!=-1) {
	   FoundPFElectron=true;
	   fT_pho_isPFElectron[phoqi] = 1;
	   PhotonToPFElectronMatchingArray[phoqi] = iel;
	 }

       }


       /*       
       { // start PF stuff from Nicolas

	 reco::PhotonCollection::const_iterator gamIterSl;

	 const Photon* gamIter = &photon;
	 
	 //e/gammma agreed recommandation
	 
	 fT_pho_Cone01PhotonIso_dEta015EB_dR070EE_mvVtx[phoqi] = 0;
	 fT_pho_Cone02PhotonIso_dEta015EB_dR070EE_mvVtx[phoqi] = 0;
	 fT_pho_Cone03PhotonIso_dEta015EB_dR070EE_mvVtx[phoqi] = 0;
	 fT_pho_Cone04PhotonIso_dEta015EB_dR070EE_mvVtx[phoqi] = 0;

	 fT_pho_Cone01NeutralHadronIso_mvVtx[phoqi] = 0;
	 fT_pho_Cone02NeutralHadronIso_mvVtx[phoqi] = 0;
	 fT_pho_Cone03NeutralHadronIso_mvVtx[phoqi] = 0;
	 fT_pho_Cone04NeutralHadronIso_mvVtx[phoqi] = 0;

	 fT_pho_Cone01ChargedHadronIso_dR02_dz02_dxy01[phoqi] = 0;
	 fT_pho_Cone02ChargedHadronIso_dR02_dz02_dxy01[phoqi] = 0;
	 fT_pho_Cone03ChargedHadronIso_dR02_dz02_dxy01[phoqi] = 0;
	 fT_pho_Cone04ChargedHadronIso_dR02_dz02_dxy01[phoqi] = 0; 

	 int PfCandType[10000];
	 float PfCandPt[10000];
	 float PfCandPx[10000];
	 float PfCandPy[10000];
	 float PfCandPz[10000];
	 float PfCandPtAtVtx[10000];
	 float PfCandPxAtVtx[10000];
	 float PfCandPyAtVtx[10000];
	 float PfCandPzAtVtx[10000];
	 float PfCandVx[10000];
	 float PfCandVy[10000];
	 float PfCandVz[10000];
	 float PfCandDxy[10000];
	 float PfCandDz[10000];
	 int PfCandIsFromPU[10000]; 
	 float PfCandDeltaRrecomputed[10000];
	 float PfCandDeltaEtarecomputed[10000];
	 float PfCandDeltaPhirecomputed[10000];
	 

	 //Recompute pflow isolation keeping all the pfcandidates in a 0.4 cone

	     
	 float photonVx = gamIter->vx();
	 float photonVy = gamIter->vy();
	 float photonVz = gamIter->vz();
	     
	 for( int i=0; i<ncand; ++i ) {
	   
	   int type = FindPFCandType((*pfCandidates)[i].pdgId());

	   if ((FoundPFPhoton && i==iphot) || (FoundPFElectron && i==iel)) {
	     storethispfcand[i]=1;
	     continue;
	   }
	   
	   if (type==0 || type==1 || type==2){
	   
	     PfCandType[i] = type;

	     PfCandPt[i] = (*pfCandidates)[i].pt();
	     PfCandPx[i] = (*pfCandidates)[i].px();
	     PfCandPy[i] = (*pfCandidates)[i].py();
	     PfCandPz[i] = (*pfCandidates)[i].pz();
	     PfCandVx[i] = (*pfCandidates)[i].vx();
	     PfCandVy[i] = (*pfCandidates)[i].vy();
	     PfCandVz[i] = (*pfCandidates)[i].vz();

	     math::XYZVector photon_direction = math::XYZVector(gamIter->superCluster()->x()-PfCandVx[i], \
								gamIter->superCluster()->y()-PfCandVy[i], \
								gamIter->superCluster()->z()-PfCandVz[i]);
	     float photon_eta = photon_direction.eta();
	     float photon_phi = photon_direction.phi();
	 
	     float pfcand_eta = (*pfCandidates)[i].eta();
	     float pfcand_phi = (*pfCandidates)[i].phi();
	     
	     PfCandDeltaRrecomputed[i] = reco::deltaR(photon_eta,photon_phi,pfcand_eta,pfcand_phi);
	     PfCandDeltaEtarecomputed[i] = photon_eta-pfcand_eta;
	     PfCandDeltaPhirecomputed[i] = reco::deltaPhi(photon_phi,pfcand_phi);

	     bool usetrackref=false;
	     if (type==1 && (*pfCandidates)[i].trackRef().isNonnull()) usetrackref=true;

	     if (usetrackref){ // use momentum (and vtx) of the track for dxy and dz calculation
	       PfCandPtAtVtx[i] = (*pfCandidates)[i].trackRef()->pt();
	       PfCandPxAtVtx[i] = (*pfCandidates)[i].trackRef()->px();
	       PfCandPyAtVtx[i] = (*pfCandidates)[i].trackRef()->py();
	       PfCandPzAtVtx[i] = (*pfCandidates)[i].trackRef()->pz();

	       bool bad=false;
	       if ((*pfCandidates)[i].vx()!=(*pfCandidates)[i].trackRef()->vx()) bad=true;
	       if ((*pfCandidates)[i].vy()!=(*pfCandidates)[i].trackRef()->vy()) bad=true;
	       if ((*pfCandidates)[i].vz()!=(*pfCandidates)[i].trackRef()->vz()) bad=true;
	       if (bad) edm::LogWarning("NTP") << "@SUB=analyze"
					       << "Something wrong with trackRef vertex for charged hadron pfcandidate";
	     }
	     else {
	       PfCandPtAtVtx[i] = PfCandPt[i];
               PfCandPxAtVtx[i] = PfCandPx[i];
               PfCandPyAtVtx[i] = PfCandPy[i];
               PfCandPzAtVtx[i] = PfCandPz[i];
	     }
	  
	     PfCandDxy[i] =  ( -(PfCandVx[i]-photonVx)*PfCandPyAtVtx[i] +(PfCandVy[i]-photonVy)*PfCandPxAtVtx[i] ) / PfCandPtAtVtx[i];
	     PfCandDz[i] = (PfCandVz[i]-photonVz) - ( (PfCandVx[i]-photonVx)*PfCandPxAtVtx[i]+(PfCandVy[i]-photonVy)*PfCandPyAtVtx[i] )/PfCandPtAtVtx[i] * PfCandPzAtVtx[i]/PfCandPtAtVtx[i];
	     PfCandDxy[i] = fabs(PfCandDxy[i]);
	     PfCandDz[i] = fabs(PfCandDz[i]);
 

	     PfCandIsFromPU[i] = -1;
	     if (type==1){
	       reco::VertexRef chvtx = chargedHadronVertex(alternativeVertexHandle, (*pfCandidates)[i]);
	       if (chvtx.isNull() || chvtx.key()==0) PfCandIsFromPU[i] = 0;
	       else PfCandIsFromPU[i] = 1;
	     }

	     double pt = PfCandPt[i];
	     double dEta = PfCandDeltaEtarecomputed[i];
	     double dPhi = PfCandDeltaPhirecomputed[i];
	     double dR = PfCandDeltaRrecomputed[i];
	     double dz = PfCandDz[i];
	     double dxy = PfCandDxy[i];

	     if (type==1 && dR<0.4) storethispfcand[i]=true;
	     if (fabs(dEta)<0.4) storethispfcand[i]=true;

	     { // determination of distance for footprint removal method
	       TVector3 photon_scposition(gamIter->superCluster()->x(),gamIter->superCluster()->y(),gamIter->superCluster()->z());
	       bool isbarrel = gamIter->isEB();
	       TVector3 pfvertex(PfCandVx[i],PfCandVy[i],PfCandVz[i]);
	       TVector3 pfmomentum(PfCandPx[i],PfCandPy[i],PfCandPz[i]);
	       pfmomentum = pfmomentum.Unit();
	       TVector3 ecalpfhit(0,0,0);
	       bool good=false;
	       if (isbarrel){
		 TGeoTube ebgeom(0,photon_scposition.Perp(),1e+10);
 		 double p[3] = {pfvertex.x(),pfvertex.y(),pfvertex.z()};
		 double d[3] = {pfmomentum.x(),pfmomentum.y(),pfmomentum.z()};
		 if (ebgeom.Contains(p)){
		   double dist = ebgeom.DistFromInside(p,d);
		   ecalpfhit = pfvertex + dist*pfmomentum;
		   good=true;
		 }
	       }
	       else { // EE
		 TGeoPara eegeom(1e+10,1e+10,fabs(photon_scposition.z()),0,0,0);
 		 double p[3] = {pfvertex.x(),pfvertex.y(),pfvertex.z()};
		 double d[3] = {pfmomentum.x(),pfmomentum.y(),pfmomentum.z()};
		 if (eegeom.Contains(p)){
		   double dist = eegeom.DistFromInside(p,d);
		   ecalpfhit = pfvertex + dist*pfmomentum;
		   good=true;
		 }
	       }
	       if (good && ecalpfhit.Perp()!=0 && photon_scposition.Perp()!=0){
		 if (fabs(ecalpfhit.Eta()-photon_scposition.Eta())<0.4) storethispfcand[i]=true;
	       }
	     }

	       
	     if (type==0){ //Neutral Hadron
	       if (dR<0.1) fT_pho_Cone01NeutralHadronIso_mvVtx[phoqi] += pt;
	       if (dR<0.2) fT_pho_Cone02NeutralHadronIso_mvVtx[phoqi] += pt;
	       if (dR<0.3) fT_pho_Cone03NeutralHadronIso_mvVtx[phoqi] += pt;
	       if (dR<0.4) fT_pho_Cone04NeutralHadronIso_mvVtx[phoqi] += pt;
	     }

	     if (type==2) { //Photon
	       bool vetoed=false;
	       if (fabs(dEta)<0.015 && gamIter->isEB()) vetoed=true;
	       else if (gamIter->isEE()){
		 float sceta = gamIter->superCluster()->eta();
		 float limit_dR = 0.00864*fabs(sinh(sceta))*4;
		 if (dR<limit_dR) vetoed=true;
	       }
	       if (!vetoed){
		 if (dR<0.1) fT_pho_Cone01PhotonIso_dEta015EB_dR070EE_mvVtx[phoqi] += pt;
		 if (dR<0.2) fT_pho_Cone02PhotonIso_dEta015EB_dR070EE_mvVtx[phoqi] += pt;
		 if (dR<0.3) fT_pho_Cone03PhotonIso_dEta015EB_dR070EE_mvVtx[phoqi] += pt;  
		 if (dR<0.4) fT_pho_Cone04PhotonIso_dEta015EB_dR070EE_mvVtx[phoqi] += pt;
	       }
	     }
	   
	     if (type==1){ //Charged Hadron
	       //dz/dxy
	       if (fabs(dz)<0.2 && fabs(dxy)<0.1 && dR>0.02){
		 if (dR<0.1) fT_pho_Cone01ChargedHadronIso_dR02_dz02_dxy01[phoqi] += pt;
		 if (dR<0.2) fT_pho_Cone02ChargedHadronIso_dR02_dz02_dxy01[phoqi] += pt;
		 if (dR<0.3) fT_pho_Cone03ChargedHadronIso_dR02_dz02_dxy01[phoqi] += pt;
		 if (dR<0.4) fT_pho_Cone04ChargedHadronIso_dR02_dz02_dxy01[phoqi] += pt;
	       }
		 
	     }
	       
	       
	   }

	 }
	     
	 
	 fT_pho_Cone03PFCombinedIso[phoqi] = (fT_pho_Cone03ChargedHadronIso_dR02_dz02_dxy01[phoqi]+fT_pho_Cone03NeutralHadronIso_mvVtx[phoqi]+fT_pho_Cone03PhotonIso_dEta015EB_dR070EE_mvVtx[phoqi]) / fTPhotPt[phoqi];
	 fT_pho_Cone04PFCombinedIso[phoqi] = (fT_pho_Cone04ChargedHadronIso_dR02_dz02_dxy01[phoqi]+fT_pho_Cone04NeutralHadronIso_mvVtx[phoqi]+fT_pho_Cone04PhotonIso_dEta015EB_dR070EE_mvVtx[phoqi]) / fTPhotPt[phoqi];


       }     // end PF stuff from Nicolas
       */


		


// DISABLED: NO SEED IN AOD (UPDATE IT IN 4_2)
// 	// Spike removal information
// 	if ( photon.superCluster()->seed()->caloID().detector( reco::CaloID::DET_ECAL_BARREL ) ) {
// 		fTPhotScSeedSeverity[phoqi] = EcalSeverityLevelAlgo::severityLevel( photon.superCluster()->seed()->seed(), *ebRecHits, *channelStatus );
// 		fTPhotE1OverE9[phoqi] = EcalSeverityLevelAlgo::E1OverE9(   photon.superCluster()->seed()->seed(), *ebRecHits );
// 		fTPhotS4OverS1[phoqi] = EcalSeverityLevelAlgo::swissCross( photon.superCluster()->seed()->seed(), *ebRecHits );
// 	} else if ( photon.superCluster()->seed()->caloID().detector( reco::CaloID::DET_ECAL_ENDCAP ) ) {
// 		fTPhotScSeedSeverity[phoqi] = EcalSeverityLevelAlgo::severityLevel( photon.superCluster()->seed()->seed(), *eeRecHits, *channelStatus );
// 		fTPhotE1OverE9[phoqi] = EcalSeverityLevelAlgo::E1OverE9(   photon.superCluster()->seed()->seed(), *eeRecHits );
// 		fTPhotS4OverS1[phoqi] = 1.0-EcalSeverityLevelAlgo::swissCross( photon.superCluster()->seed()->seed(), *eeRecHits );
//      } else
// 			edm::LogWarning("NTP") << "Photon supercluster seed crystal neither in EB nor in EE!";



 	} // end photon loop


	
       /*

	 USAGE OF VERTEX CHOICE FOR DIPHOTON EVENTS:

	 diphotons_{first,second} are vectors of {photon_1_index,photon_2_index}

	 vtx_dipho_??? are, for each diphoton pair, vectors of vertex indices (as ranked by the different algos)

	 For example: best vertex for diphoton pair 3, with photon_1_index=diphotons_first[3] and photon_2_index=diphotons_second[3]: vtx_dipho_bla[3].at(0), second choice vtx_dipho_bla[3].at(1) ...

       */

	/*
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

if (VTX_MVA_DEBUG)	     	     cout << "working on vtx " << i << endl;

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

	 conv_n=0;

	 { // all conversions

	     
	   for( reco::ConversionCollection::const_iterator  iConv = convH->begin(); iConv != convH->end(); iConv++) {
	          
	     reco::Conversion localConv = reco::Conversion(*iConv);
  
	     if(ConversionsCut(localConv)) continue;

	     if (conv_n >= gMaxnconv){
	       edm::LogWarning("NTP") << "@SUB=analyze"
				      << "Maximum number of conversions exceeded";
	       fTgoodevent = 1;
	       break;
	     }

	     conv_validvtx[conv_n]=localConv.conversionVertex().isValid();
  
	     if ( localConv.conversionVertex().isValid() ){
	       reco::Vertex vtx=localConv.conversionVertex();
	       conv_vtx[conv_n].SetXYZ(vtx.x(), vtx.y(), vtx.z());
	       conv_ntracks[conv_n]=localConv.nTracks();
	       conv_chi2_probability[conv_n]=ChiSquaredProbability(vtx.chi2(), vtx.ndof());
	       conv_eoverp[conv_n]=localConv.EoverPrefittedTracks();
	       conv_zofprimvtxfromtrks[conv_n]=localConv.zOfPrimaryVertexFromTracks();
	       conv_refitted_momentum[conv_n].SetXYZ(localConv.refittedPairMomentum().x(), localConv.refittedPairMomentum().y(), localConv.refittedPairMomentum().z());
	     }

	     conv_n++;

	   }

	 }


	 if (VTX_MVA_DEBUG)	  cout << "done convs" << endl;

	 std::vector<std::string> vtxVarNames;
	 vtxVarNames.push_back("ptbal"), vtxVarNames.push_back("ptasym"), vtxVarNames.push_back("logsumpt2");

	 if (VTX_MVA_DEBUG)	 cout << "ready: remember delete readers" << endl;	 


	 if( fTnphotons < 2 ) {

	   vtx_dipho_h2gglobe.push_back( std::vector<int>() );
	   vtx_dipho_mva.push_back( std::vector<int>() );
	   vtx_dipho_productrank.push_back( std::vector<int>() );
	   
	   for(int ii=0;ii<fTnvrtx; ++ii) {vtx_dipho_h2gglobe.back().push_back(ii); }
	   for(int ii=0;ii<fTnvrtx; ++ii) {vtx_dipho_mva.back().push_back(ii); }
	   for(int ii=0;ii<fTnvrtx; ++ii) {vtx_dipho_productrank.back().push_back(ii); }

	 } else {

	   if (VTX_MVA_DEBUG)	   cout << "temp" << endl;

	   // fully combinatorial vertex selection
	   for(int ip=0; ip<fTnphotons; ++ip) {
	     for(int jp=ip+1; jp<fTnphotons; ++jp) {
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
	     
	     ETHVertexInfo vinfo(int(fTnvrtx),
				 fTvrtxx+0,
				 fTvrtxy+0,
				 fTvrtxz+0,
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

	     vAna->analyze(vinfo,pho1,pho2);

	     if (VTX_MVA_DEBUG)	     cout << "initialized vAna" << endl;

	     // make sure that vertex analysis indexes are in synch 
	     assert( int(id) == vAna->pairID(ipho1,ipho2) );

	     if (VTX_MVA_DEBUG)	     cout << "starting rankings" << endl;

	     if (VTX_MVA_DEBUG)	     cout << "rankprod" << endl;
	     /// rank product vertex selection. Including pre-selection based on conversions information.
	     vtx_dipho_productrank.push_back(vAna->rankprod(rankVariables));


	     if (VTX_MVA_DEBUG)	     cout << "mva pasquale" << endl;
	     /// MVA vertex selection
	     vtx_dipho_mva.push_back(vAna->rank(*perVtxReader,perVtxMvaMethod));


	     // vertex probability through per-event MVA (not used so far)
	     // float vtxEvtMva = vAna->perEventMva( *perEvtReader,  perEvtMvaMethod, vtx_dipho_mva->back() );
	     // float vtxProb = vAna->vertexProbability( vtxEvtMva );
	     
	     if (VTX_MVA_DEBUG)	     cout << "mva hgg globe" << endl;
	     // Globe vertex selection with conversions
	     vtx_dipho_h2gglobe.push_back(HggVertexSelection(*vAna, *vConv, pho1, pho2, vtxVarNames,false,0,""));


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
	

         { // write output of vertexing

           int i1;
           int i2;

           // [diphoton_pair][ranking_of_vertices]

           i1=0;
           for (std::vector<int>::const_iterator it=diphotons_first.begin(); it!=diphotons_first.end() && i1<gMax_vertexing_diphoton_pairs; it++) {f_diphotons_first[i1]=*it; i1++;}
           i1=0;
           for (std::vector<int>::const_iterator it=diphotons_second.begin(); it!=diphotons_second.end() && i1<gMax_vertexing_diphoton_pairs; it++) {f_diphotons_second[i1]=*it; i1++;}

           for (i1=0; i1<(int)vtx_dipho_h2gglobe.size() && i1<gMax_vertexing_diphoton_pairs; i1++)
             for (i2=0; i2<(int)vtx_dipho_h2gglobe.at(i1).size() && i2<gMax_vertexing_vtxes; i2++) {
       if (!(i1<gMax_vertexing_diphoton_pairs && i2<gMax_vertexing_vtxes)) cout<<"wrong!!!"<<endl;
               f_vtx_dipho_h2gglobe[i1][i2]=vtx_dipho_h2gglobe.at(i1).at(i2); }

           for (i1=0; i1<(int)vtx_dipho_mva.size() && i1<gMax_vertexing_diphoton_pairs; i1++)
             for (i2=0; i2<(int)vtx_dipho_mva.at(i1).size() && i2<gMax_vertexing_vtxes; i2++) {
       if (!(i1<gMax_vertexing_diphoton_pairs && i2<gMax_vertexing_vtxes)) cout<<"wrong!!!"<<endl;
               f_vtx_dipho_mva[i1][i2]=vtx_dipho_mva.at(i1).at(i2); }

           for (i1=0; i1<(int)vtx_dipho_productrank.size() && i1<gMax_vertexing_diphoton_pairs; i1++)
             for (i2=0; i2<(int)vtx_dipho_productrank.at(i1).size() && i2<gMax_vertexing_vtxes; i2++) {
       if (!(i1<gMax_vertexing_diphoton_pairs && i2<gMax_vertexing_vtxes)) cout<<"wrong!!!"<<endl;
               f_vtx_dipho_productrank[i1][i2]=vtx_dipho_productrank.at(i1).at(i2); }

         }


       } // end vertex selection for diphoton events



       //       cout << "end vertex selection MVA" << endl;

       */


	std::vector<std::vector<int> > Jets_PfCand_content;

	////////////////////////////////////////////////////////
	// Jet Variables:
	const JetCorrector* jetCorr = JetCorrector::getJetCorrector(fJetCorrs, iSetup);
	std::vector<OrderPair> corrIndices;  // Vector of indices and pt of corr. jets to re-order them
	int iraw(0);
	fTnjetstot = jets->size();
	// Loop over uncorr. jets
	for(View<Jet>::const_iterator Jit = jets->begin(); Jit != jets->end(); ++Jit, ++iraw){
		// Cut on uncorrected pT (for startup)
		if(Jit->pt() < fMinrawjpt) continue;
		// Save only the gMaxnjets first uncorrected jets
		if(iraw >= gMaxnjets){
			edm::LogWarning("NTP") << "@SUB=analyze"
                                               << "Found more than " << static_cast<int>(gMaxnjets) << " uncorrected jets, I'm scared ...";
			fTflagmaxujetexc = 1;
			fTgoodevent = 1;
			break;
		}
                
                JetBaseRef jetRef(edm::Ref<JetView>(jets,iraw));
		double scale = jetCorr->correction(*Jit,jetRef,iEvent,iSetup);
		corrIndices.push_back(make_pair(iraw, scale*Jit->pt()));
	}
	
	// Sort corrected jet collection by decreasing pt
	std::sort(corrIndices.begin(), corrIndices.end(), indexComparator);
	
	// Determine corrected jets
	int jqi(-1); // counts # of qualified jets
	// Loop over corr. jet indices
	for(std::vector<OrderPair>::const_iterator it = corrIndices.begin(); it != corrIndices.end(); ++it ) {
		// Check if maximum number of jets is exceeded already
		if(jqi >= gMaxnjets-1) {
                  edm::LogWarning("NTP") << "@SUB=analyze"
				<< "Maximum number of jets exceeded";
			fTflagmaxjetexc = 1;
			fTgoodevent = 1;
			break;
		}
		int index = it->first;
		const PFJet* cojet = static_cast<const PFJet*>( &((*jets)[index]) ); // look away...
                std::auto_ptr<PFJet> jet(new PFJet(*cojet));

		// The correction was calculated above: use it
		double scale = it->second/jet->pt();
		jet->scaleEnergy(scale);
	
		// Jet preselection
		if(jet->pt() < fMincorjpt) continue; 
		if(fabs(jet->eta()) > fMaxjeta) continue;
		jqi++;

		// Dump jet properties into tree variables
		fTjpx    [jqi] = jet->px();
		fTjpy    [jqi] = jet->py();
		fTjpz    [jqi] = jet->pz();
		fTjpt    [jqi] = jet->pt();
		fTjeta   [jqi] = jet->eta();
		fTjphi   [jqi] = jet->phi();
		fTje     [jqi] = jet->energy();
		fTjet    [jqi] = jet->et();
		fTjEcorr [jqi] = scale;
		fTJEtaRms[jqi] = sqrt(jet->etaetaMoment());
		fTJPhiRms[jqi] = sqrt(jet->etaphiMoment());
		fTjArea  [jqi] = jet->jetArea();

		fTjNconstituents[jqi] = jet->nConstituents();
		fTjChMult       [jqi] = jet->chargedMultiplicity(); // do it the pf way...
		fTjNeuMult      [jqi] = jet->neutralMultiplicity(); 
		
		// energy fractions for JID need to be computed w.r.t. uncorrected jet energy!!
		// see for instance https://twiki.cern.ch/twiki/bin/view/CMS/JetID
		// or http://cmssdt.cern.ch/SDT/doxygen/CMSSW_4_1_3/doc/html/dc/dd5/classPFJetIDSelectionFunctor.html  
		double uncorr_energy  = jet->energy()/scale;
		fTjChHadFrac     [jqi] = jet->chargedHadronEnergy()/uncorr_energy;
		fTjNeuHadFrac    [jqi] = jet->neutralHadronEnergy()/uncorr_energy + jet->HFHadronEnergy()/uncorr_energy;
		fTjChEmFrac      [jqi] = jet->chargedEmEnergy()/uncorr_energy;
		fTjNeuEmFrac     [jqi] = jet->neutralEmEnergy()/uncorr_energy;
		fTjChMuEFrac     [jqi] = jet->chargedMuEnergy()/uncorr_energy;
		fTjPhoFrac       [jqi] = jet->photonEnergy()/uncorr_energy; // photons also count for neutralEmEnergy
		fTjHFHadFrac     [jqi] = jet->HFHadronEnergy()/uncorr_energy;
		fTjHFEMFrac      [jqi] = jet->HFEMEnergy()/uncorr_energy;   // also contained in neutralEmEnergy
		// see CMSSW/RecoJets/JetProducers/src/JetSpecific.cc

		vector<PFCandidatePtr> JetpfCandidates = jet->getPFConstituents();
		assert((int)(Jets_PfCand_content.size())==jqi);
		Jets_PfCand_content.push_back(std::vector<int>());
		
		float sumPt_cands=0.;
		float sumPt2_cands=0.;
		float rms_cands=0.;
		
		TLorentzVector jetp4;
		jetp4.SetPtEtaPhiE(jet->pt(), jet->eta(), jet->phi(), jet->energy());

		for (vector<PFCandidatePtr>::const_iterator jCand = JetpfCandidates.begin(); jCand != JetpfCandidates.end(); ++jCand) {
		  
		  math::XYZTLorentzVectorD const& pCand_t = (*jCand)->p4();
		  TLorentzVector pCand(pCand_t.px(), pCand_t.py(), pCand_t.pz(), pCand_t.energy());
		  if(pCand.Pt()>0.){
		    sumPt_cands += pCand.Pt();
		    sumPt2_cands += (pCand.Pt()*pCand.Pt());

		    float deltaR = pCand.DeltaR(jetp4);
		    rms_cands += (pCand.Pt()*pCand.Pt()*deltaR*deltaR);

		  }

		  for (int j=0; j<(int)(pfCandidates->size()); j++){
		    if ((*jCand) == edm::Ptr<reco::PFCandidate>(pfCandidates,j)) Jets_PfCand_content.back().push_back(j);
		  }

		} //for PFCandidates

		fTjPtD      [jqi] = sqrt( sumPt2_cands )/sumPt_cands;
		fTjRMSCand  [jqi] = rms_cands/sumPt2_cands;



		// Calculate the DR wrt the closest electron
		float ejDRmin = 10.; // Default when no electrons previously selected
		for( int j = 0; j < fTneles; j++ ){
			float ejDR = reco::deltaR(jet->eta(), jet->phi(), fTeeta[j], fTephi[j]);
			if(ejDR<ejDRmin) ejDRmin = ejDR;
		}
		fTjeMinDR[jqi] = ejDRmin;

		// B-tagging probability (for 4 b-taggings)
		// remember: 'index' is the index of the uncorrected jet, as saved in the event
		fTjbTagProbTkCntHighEff[jqi]  = (*jetsAndProbsTkCntHighEff) [index].second;
		fTjbTagProbTkCntHighPur[jqi]  = (*jetsAndProbsTkCntHighPur) [index].second;
		fTjbTagProbSimpSVHighEff[jqi] = (*jetsAndProbsSimpSVHighEff)[index].second;
		fTjbTagProbSimpSVHighPur[jqi] = (*jetsAndProbsSimpSVHighPur)[index].second;

		// Jet-track association: get associated tracks
		const reco::TrackRefVector& tracks = jet->getTrackRefs();
		std::vector<const reco::Track*> AssociatedTracks;
		for( TrackRefVector::iterator it = tracks.begin(); it != tracks.end(); ++it ) AssociatedTracks.push_back( it->get() );
			
		// Below save the momenta of the three leading tracks associated to the jet
		float pT1(0.), pT2(0.), pT3(0.);
		int idx1(-1), idx2(-1), idx3(-1);
			
		// Jet-track association: make transient tracks and store information
		std::vector<TransientTrack> AssociatedTTracks;
		fTjMass[jqi] = 0.;
		if(fabs(jet->eta())<2.9){ // when the cone of dR=0.5 around the jet is (at least partially) inside the tracker acceptance
			// Tmp variables for vectorial sum of pt of tracks
			double pXtmp(0.), pYtmp(0.), pZtmp(0.), E2tmp(0.);
			const double trkmass = 0.; // Assumed mass for tracks

			// Loop over associated tracks:
			for(size_t t = 0; t < AssociatedTracks.size(); ++t){
				AssociatedTTracks.push_back(theB->build(AssociatedTracks[t])); // build transient tracks for vertex fitting below
				if(AssociatedTracks[t]->normalizedChi2()<10. && AssociatedTracks[t]->numberOfValidHits()>10 && AssociatedTracks[t]->pt()>1.){
					pXtmp += AssociatedTracks[t]->px();
					pYtmp += AssociatedTracks[t]->py();
					pZtmp += AssociatedTracks[t]->pz();
					E2tmp += trkmass*trkmass + pXtmp*pXtmp + pYtmp*pYtmp + pZtmp*pZtmp;
				}
				// Find the three highest pT tracks
				if(AssociatedTracks[t]->pt() > pT1 && AssociatedTracks.size() >= 1){
					pT1=AssociatedTracks[t]->pt();
					idx3=idx2;
					idx2=idx1;
					idx1=t;
				} else if (AssociatedTracks[t]->pt() < pT1 && AssociatedTracks[t]->pt() > pT2 && AssociatedTracks.size() >= 2) {
					pT2=AssociatedTracks[t]->pt();
					idx3=idx2;
					idx2=t;
				} else if (AssociatedTracks[t]->pt() < pT2 && AssociatedTracks[t]->pt() > pT3 && AssociatedTracks.size() >= 3){
					pT3=AssociatedTracks[t]->pt();
					idx3=t;
				}
			}
			// Fill the momenta
			if(AssociatedTracks.size()>=1){
				fTjtrk1px[jqi] = AssociatedTracks[idx1]->px();
				fTjtrk1py[jqi] = AssociatedTracks[idx1]->py();
				fTjtrk1pz[jqi] = AssociatedTracks[idx1]->pz();
			}
			if(AssociatedTracks.size()>=2){
				fTjtrk2px[jqi] = AssociatedTracks[idx2]->px();
				fTjtrk2py[jqi] = AssociatedTracks[idx2]->py();
				fTjtrk2pz[jqi] = AssociatedTracks[idx2]->pz();
			}
			if(AssociatedTracks.size()>=3){
				fTjtrk3px[jqi] = AssociatedTracks[idx3]->px();
				fTjtrk3py[jqi] = AssociatedTracks[idx3]->py();
				fTjtrk3pz[jqi] = AssociatedTracks[idx3]->pz();
			}
			
			fTjMass[jqi]   = sqrt(E2tmp - pXtmp*pXtmp - pYtmp*pYtmp - pZtmp*pZtmp);
			if(fTjChMult[jqi] > 0) fTjMass[jqi] *= fTjNconstituents[jqi]/fTjChMult[jqi]; // apparantly there ARE cases where ChMult is 0, but we still end up here...
			else fTjMass[jqi] = 0.;
		} else { // The whole cone used for jet-tracks association is outside of the tracker acceptance
			fTjMass[jqi] = -888.88;
		}

		// Do a vertex fitting with the tracks
		if(AssociatedTTracks.size() > 1) {
			TransientVertex jetVtx = avFitter.vertex(AssociatedTTracks);
			if(jetVtx.isValid()){
				fTjetVtxx    [jqi] = jetVtx.position().x();
				fTjetVtxy    [jqi] = jetVtx.position().y();
				fTjetVtxz    [jqi] = jetVtx.position().z();
				fTjetVtxExx  [jqi] = jetVtx.positionError().cxx();
				fTjetVtxEyx  [jqi] = jetVtx.positionError().cyx();
				fTjetVtxEyy  [jqi] = jetVtx.positionError().cyy();
				fTjetVtxEzy  [jqi] = jetVtx.positionError().czy();
				fTjetVtxEzz  [jqi] = jetVtx.positionError().czz();
				fTjetVtxEzx  [jqi] = jetVtx.positionError().czx();
				fTjetVtxNChi2[jqi] = jetVtx.normalisedChiSquared();
			}else{
				fTjetVtxx    [jqi] = -777.77;
				fTjetVtxy    [jqi] = -777.77;
				fTjetVtxz    [jqi] = -777.77;
				fTjetVtxExx  [jqi] = -777.77;
				fTjetVtxEyx  [jqi] = -777.77;
				fTjetVtxEyy  [jqi] = -777.77;
				fTjetVtxEzy  [jqi] = -777.77;
				fTjetVtxEzz  [jqi] = -777.77;
				fTjetVtxEzx  [jqi] = -777.77;
				fTjetVtxNChi2[jqi] = -777.77;
			}
		}else{
			fTjetVtxx    [jqi] = -888.88;
			fTjetVtxy    [jqi] = -888.88;
			fTjetVtxz    [jqi] = -888.88;
			fTjetVtxExx  [jqi] = -888.88;
			fTjetVtxEyx  [jqi] = -888.88;
			fTjetVtxEyy  [jqi] = -888.88;
			fTjetVtxEzy  [jqi] = -888.88;
			fTjetVtxEzz  [jqi] = -888.88;
			fTjetVtxEzx  [jqi] = -888.88;
			fTjetVtxNChi2[jqi] = -888.88;
		}
		AssociatedTracks.clear();
		AssociatedTTracks.clear();
	
		// GenJet matching
		if (!fIsRealData && fTNGenJets > 0) fTjetGenJetIndex[jqi] = matchJet(&(*jet));
		fTgoodjet[jqi] = 0;
	}
	fTnjets = jqi+1;
	// corrJets.clear();
	corrIndices.clear();

	// Check electron duplication
	ElectronDuplicate(elecPtr, trckPtr);
	// Check photon/electron duplication
	PhotonElectronDuplicate(elecPtr, photSCs);

	// These don't work for pf jets yet
	// Check electron/jet duplication
	// ElJetOverlap(jetPtr, elecPtr, calotowers);
	// Check photon/jet duplication
	// PhotonJetOverlap(jetPtr, photSCs, calotowers);


	////////////////////////////////////////////////////////
	// PfCandidates Variables:

	int pfcandIndex(0);
	for (unsigned int i=0; i<pfCandidates->size(); i++){

	  int type = FindPFCandType((*pfCandidates)[i].pdgId());
	  if (type==2) storethispfcand[i]=true;

	  if (storethispfcand[i]==0) continue;

	  if (pfcandIndex >= gMaxnpfcand){
	    edm::LogWarning("NTP") << "@SUB=analyze"
				   << "Maximum number of pf candidates exceeded";
	    fTgoodevent = 1;
	    break;
	  }

	  fTPfCandPdgId[pfcandIndex] = (*pfCandidates)[i].pdgId();
	  fTPfCandPt[pfcandIndex] = (*pfCandidates)[i].pt();
	  fTPfCandEta[pfcandIndex] = (*pfCandidates)[i].eta();
	  fTPfCandPhi[pfcandIndex] = (*pfCandidates)[i].phi();
//	  fTPfCandPx[pfcandIndex] = (*pfCandidates)[i].px();
//	  fTPfCandPy[pfcandIndex] = (*pfCandidates)[i].py();
//	  fTPfCandPz[pfcandIndex] = (*pfCandidates)[i].pz();	  
	  fTPfCandEnergy[pfcandIndex] = (*pfCandidates)[i].energy();
	  fTPfCandVx[pfcandIndex] = (*pfCandidates)[i].vx();
	  fTPfCandVy[pfcandIndex] = (*pfCandidates)[i].vy();
	  fTPfCandVz[pfcandIndex] = (*pfCandidates)[i].vz();
//	  fTPfCandMomX[pfcandIndex] = (*pfCandidates)[i].momentum().x();
//	  fTPfCandMomY[pfcandIndex] = (*pfCandidates)[i].momentum().y();
//	  fTPfCandMomZ[pfcandIndex] = (*pfCandidates)[i].momentum().z();

	  for (int j=0; j<fTnphotons; j++){
	    if (PhotonToPFPhotonMatchingArray[j]==(int)i) fT_pho_matchedPFPhotonCand[j]=pfcandIndex;
	    if (PhotonToPFElectronMatchingArray[j]==(int)i) fT_pho_matchedPFElectronCand[j]=pfcandIndex;
	  }

	  for (int jet=0; jet<(int)(Jets_PfCand_content.size()); jet++){
	    for (int k=0; k<(int)(Jets_PfCand_content.at(jet).size()); k++){
	      if (Jets_PfCand_content.at(jet).at(k)==(int)i){
		if (fTPfCandBelongsToJet[pfcandIndex]>=0) {cout << "WRONG: this PFCand is already in another jet!" << endl; continue;}
		fTPfCandBelongsToJet[pfcandIndex]=jet;
	      }
	    }
	  }


	  /*
	  if ( (type==1) && ( !((*pfCandidates)[i].trackRef()) ) ) type=-1;
	  reco::HitPattern pattern; 
	  if (type==1) pattern=(*pfCandidates)[i].trackRef()->hitPattern(); 
	  fTPfCandHasHitInFirstPixelLayer[pfcandIndex] = (type==1) ? (pattern.hasValidHitInFirstPixelBarrel() || pattern.hasValidHitInFirstPixelEndcap()) : -999;
	  fTPfCandTrackRefPx[pfcandIndex] = (type==1) ? (*pfCandidates)[i].trackRef()->px() : -999;
	  fTPfCandTrackRefPy[pfcandIndex] = (type==1) ? (*pfCandidates)[i].trackRef()->py() : -999;
	  fTPfCandTrackRefPz[pfcandIndex] = (type==1) ? (*pfCandidates)[i].trackRef()->pz() : -999;
	  fTPfCandTrackRefVx[pfcandIndex] = (type==1) ? (*pfCandidates)[i].trackRef()->vx() : -999;
	  fTPfCandTrackRefVy[pfcandIndex] = (type==1) ? (*pfCandidates)[i].trackRef()->vy() : -999;
	  fTPfCandTrackRefVz[pfcandIndex] = (type==1) ? (*pfCandidates)[i].trackRef()->vz() : -999;
	  */

	  pfcandIndex++;

	}
	fTNPfCand=pfcandIndex;







	////////////////////////////////////////////////////////
	// Process other jet collections, as configured
	for ( std::vector<JetFillerBase*>::iterator it = jetFillers.begin();
              it != jetFillers.end(); ++it )
          (*it)->fillBranches(iEvent,iSetup);
	for ( std::vector<PatMuonFiller*>::iterator it = muonFillers.begin(); 
             it != muonFillers.end(); ++it ) 
          (*it)->fillBranches(iEvent,iSetup);
        for ( std::vector<PatElectronFiller*>::iterator it = electronFillers.begin(); 
              it != electronFillers.end(); ++it ) 
          (*it)->fillBranches(iEvent,iSetup);
	for ( std::vector<PatTauFiller*>::iterator it = tauFillers.begin(); 
              it != tauFillers.end(); ++it ) 
          (*it)->fillBranches(iEvent,iSetup);

	////////////////////////////////////////////////////////
	// Get and Dump (M)E(T) Variables:
	// Tracks:
	int nqtrk(-1);
	fTntrackstot = tracks->size();
	fTTrkPtSumx = 0.; fTTrkPtSumy = 0.;
	for( TrackCollection::const_iterator it = tracks->begin(); it != tracks->end() ; ++it ){
		fTTrkPtSumx += it->px(); // Calculated for ALL tracks
		fTTrkPtSumy += it->py();
		if(it->pt() < fMintrkpt) continue;
		if(fabs(it->eta()) > fMaxtrketa) continue;
		if(it->normalizedChi2() > fMaxtrknchi2) continue;
		if(it->numberOfValidHits() < fMintrknhits) continue;
		nqtrk++; // starts at 0
		// Check if maximum number of tracks is exceeded already
		if(nqtrk >= gMaxntrks) {
			edm::LogWarning("NTP") << "@SUB=analyze"
				<< "Maximum number of tracks exceeded";
			fTflagmaxtrkexc = 1;
			fTgoodevent = 1;
			break;
		}
		fTtrkpt[nqtrk]    = it->pt()*it->charge();
		fTtrketa[nqtrk]   = it->eta();
		fTtrkphi[nqtrk]   = it->phi();
		fTtrknchi2[nqtrk] = it->normalizedChi2();
		fTtrknhits[nqtrk] = it->numberOfValidHits();
		fTtrkVtxDz[nqtrk] = it->dz(primVtx->position());
		fTtrkVtxDxy[nqtrk]= it->dxy(primVtx->position());	
		fTgoodtrk[nqtrk]  = 0;
	}
	fTntracks = nqtrk+1;

	fTTrkPtSum = sqrt(fTTrkPtSumx*fTTrkPtSumx + fTTrkPtSumy*fTTrkPtSumy);
	TVector3 trkPtSum(fTTrkPtSumx, fTTrkPtSumy, 0.);
	fTTrkPtSumphi = trkPtSum.Phi();

	// Calotowers:
	fTNCaloTowers = calotowers->size();
	fTECALEsumx = 0.; fTECALEsumy = 0.; fTECALEsumz = 0.;
	fTHCALEsumx = 0.; fTHCALEsumy = 0.; fTHCALEsumz = 0.;
	fTSumEt = 0.; fTECALSumEt = 0.; fTHCALSumEt = 0.;
	for(CaloTowerCollection::const_iterator itow = calotowers->begin();
	itow!=calotowers->end(); ++itow ){
		if(itow->energy() == 0.) continue; // Check against zero energy towers
		fTSumEt += itow->et();
		fTECALSumEt += itow->emEt();
		fTHCALSumEt += itow->hadEt();
		double emFrac = itow->emEnergy()/itow->energy();
		double hadFrac = itow->hadEnergy()/itow->energy();
		fTECALEsumx += itow->px()*emFrac;
		fTECALEsumy += itow->py()*emFrac;
		fTECALEsumz += itow->pz()*emFrac;
		fTHCALEsumx += itow->px()*hadFrac;
		fTHCALEsumy += itow->py()*hadFrac;
		fTHCALEsumz += itow->pz()*hadFrac;
	}
	TVector3 ecalMET(fTECALEsumx, fTECALEsumy, fTECALEsumz);
	TVector3 hcalMET(fTHCALEsumx, fTHCALEsumy, fTHCALEsumz);
	fTECALMET    = ecalMET.Mag();
	fTECALMETphi = ecalMET.Phi();
	fTHCALMET    = hcalMET.Mag();
	fTHCALMETphi = hcalMET.Phi();
	if(fTECALEsumz != 0.) fTECALMETeta = ecalMET.Eta();
	else fTECALMETeta = 0.;
	if(fTHCALEsumz != 0.) fTHCALMETeta = hcalMET.Eta();
	else fTHCALMETeta = 0.;

	// MET Collections:
	fTRawMET             = (calomet->at(0)).pt();
	fTRawMETpx           = (calomet->at(0)).px();
	fTRawMETpy           = (calomet->at(0)).py();
	fTRawMETphi          = (calomet->at(0)).phi();
	fTRawMETemEtFrac     = (calomet->at(0)).emEtFraction();
	fTRawMETemEtInEB     = (calomet->at(0)).emEtInEB();
	fTRawMETemEtInEE     = (calomet->at(0)).emEtInEE();
	fTRawMETemEtInHF     = (calomet->at(0)).emEtInHF();
	fTRawMEThadEtFrac    = (calomet->at(0)).etFractionHadronic();
	fTRawMEThadEtInHB    = (calomet->at(0)).hadEtInHB();
	fTRawMEThadEtInHE    = (calomet->at(0)).hadEtInHE();
	fTRawMEThadEtInHF    = (calomet->at(0)).hadEtInHF();
	fTRawMETSignificance = (calomet->at(0)).significance();

	if (!fIsRealData) {
		Handle<View<GenMET> > GenMET;
		iEvent.getByLabel(fGenMETTag, GenMET);
		fTGenMET    = (GenMET->front()).pt();
		fTGenMETpx  = (GenMET->front()).px();
		fTGenMETpy  = (GenMET->front()).py();
		fTGenMETphi = (GenMET->front()).phi();
	}

	fTTCMET    = (tcmet->at(0)).pt();
	fTTCMETpx  = (tcmet->at(0)).px();
	fTTCMETpy  = (tcmet->at(0)).py();
	fTTCMETphi = (tcmet->at(0)).phi();
	fTTCMETSignificance = (tcmet->at(0)).significance();

	fTPFMET    = (pfmet->front()).pt();
	fTPFMETpx  = (pfmet->front()).px();
	fTPFMETpy  = (pfmet->front()).py();
	fTPFMETphi = (pfmet->front()).phi();
	fTPFMETSignificance = (pfmet->front()).significance();
        fTPFSumEt  = (pfmet->front()).sumEt();

	fTMuJESCorrMET    = (corrmujesmet->at(0)).pt();
	fTMuJESCorrMETpx  = (corrmujesmet->at(0)).px();
	fTMuJESCorrMETpy  = (corrmujesmet->at(0)).py();
	fTMuJESCorrMETphi = (corrmujesmet->at(0)).phi();

	if(fTnjets > 1){
		double dPhiMJ1 = TMath::Abs(reco::deltaPhi(fTjphi[0], fTMuJESCorrMETphi));
		double dPhiMJ2 = TMath::Abs(reco::deltaPhi(fTjphi[1], fTMuJESCorrMETphi));
		fTMETR12 = TMath::Sqrt(dPhiMJ1*dPhiMJ1 + (TMath::Pi()-dPhiMJ2)*(TMath::Pi()-dPhiMJ2) );
		fTMETR21 = TMath::Sqrt(dPhiMJ2*dPhiMJ2 + (TMath::Pi()-dPhiMJ1)*(TMath::Pi()-dPhiMJ1) );
	}

	fTPFMETPAT     = (pfMETpat->front()).pt();
	fTPFMETPATpx   = (pfMETpat->front()).px();
	fTPFMETPATpy   = (pfMETpat->front()).py();
	fTPFMETPATphi  = (pfMETpat->front()).phi();
	fTPFMETPATSignificance = (pfMETpat->at(0)).significance();



        ////////////////////////////////////////////////////////////////////////////////
	// Special stuff for Model Scans ///////////////////////////////////////////////

	fTxSMS=-1;
	fTxbarSMS=-1;

	if(!fIsRealData&&fIsModelScan) {
		Handle<LHEEventProduct> product;
		bool LHEEventProduct_found=iEvent.getByLabel("source", product);
		LHEEventProduct::comments_const_iterator c_begin = product->comments_begin();
		LHEEventProduct::comments_const_iterator c_end = product->comments_end();

		float mGL;
		float mLSP;
	        float xCHI;

	        for(LHEEventProduct::comments_const_iterator cit=c_begin; cit!=c_end; ++cit) {
	          size_t found = (*cit).find("model");

	          //# model T5zz_0.5_925.0_400.0to1000.0_450.0.lhe
	          if( found != std::string::npos)   {
	            size_t foundLength = (*cit).size();
	            found = (*cit).find("T5zz");
	            std::string smaller = (*cit).substr(found+1,foundLength);
	            found = smaller.find("_");
	            smaller = smaller.substr(found+1,smaller.size());
	            std::istringstream iss(smaller);
	            iss >> xCHI;
	            iss.clear();
	            found = smaller.find("_");
	            smaller = smaller.substr(found+1,smaller.size());
	            iss.str(smaller);
	            iss >> mGL;
	            iss.clear();
	            found = smaller.find("_");
	            smaller = smaller.substr(found+1,smaller.size());
	            iss.str(smaller);
	            iss >> mLSP;
	            iss.clear();

                     //model msugra_1900_100_10_0_1
	             float m0,m12,tanb,A0,signMu;
	             foundLength = (*cit).size();
	             found = (*cit).find("=");
	             smaller = (*cit).substr(found+1,foundLength);
	             found = smaller.find("_");
	             smaller = smaller.substr(found+1,smaller.size());
	             //
	             iss.clear();
	             iss.str(smaller);
	             iss >> m0;
	             iss.clear();
	             //
	             found = smaller.find("_");
	             smaller = smaller.substr(found+1,smaller.size());
	             iss.str(smaller);
	             iss >> m12;
	             iss.clear();
	             // 	 
	             found = smaller.find("_");
	             smaller = smaller.substr(found+1,smaller.size());
	             iss.str(smaller);
	             iss >> tanb;
	             iss.clear();
	             // 	 
	             found = smaller.find("_");
	             smaller = smaller.substr(found+1,smaller.size());
	             iss.str(smaller);
	             iss >> A0;
	             iss.clear();
	             // 	 
	             found = smaller.find("_");
	             smaller = smaller.substr(found+1,smaller.size());
	             iss.str(smaller);
	             iss >> signMu;
	             iss.clear();

	             // mSUGRA scan
	             fTSUSYScanM0=m0;
	             fTSUSYScanM12=m12;
	             fTSUSYScanMu=signMu;
	             fTSUSYScanA0=A0;

	          }
	        }
	        fTMassGlu = mGL;
	        fTMassChi = xCHI;
	        fTMassLSP = mLSP;
	}
	

	bool blabalot=false;
	bool BloatWithGenInfo=false;

	
	if(BloatWithGenInfo && !fIsRealData) {
	  int genIndex[1200];
	  int genNIndex[1200];
	  int genID[1200];
	  float genPx[1200];
	  float genPy[1200];
	  float genPz[1200];
	  float genPt[1200];
	  float genEta[1200];
	  float genPhi[1200];
	  float genM[1200];
	  int genMo1Index[1200];
	  int genMo2Index[1200];
	  int genNMo[1200];
	  int genStatus[1200];
	  int Promptness[1200];
	  bool StoreFlag[1200];
	  
	  int nGenParticles=0;
	  
	  Handle<GenParticleCollection> genParticles;
	  iEvent.getByLabel("genParticles", genParticles);
	  
	  // STEP 1: Loop over all particles and store the information.
	  for(size_t i = 0; i < genParticles->size(); ++ i) {
	    nGenParticles++;
	    const GenParticle & p = (*genParticles)[i];
	    genIndex[i]=i;
	    genNIndex[i]=0;
	    genID[i]=p.pdgId();
	    genPx[i]=p.px();
	    genPy[i]=p.py();
	    genPz[i]=p.pz();
	    genPt[i]=p.pt();
	    genPhi[i]=p.phi();
	    genEta[i]=p.eta();
	    genM[i]=p.mass();
	    genStatus[i]=p.status();
	    genMo1Index[i]=-1;
	    genMo2Index[i]=-1;
	    genNMo[i]=p.numberOfMothers();
	    StoreFlag[i]=false;
	    if(genID[i]==2212) Promptness[i]=0;
	    else Promptness[i]=-1;
	    

	    if(blabalot) cout << "Reading particle " << i << " (is pdgid=" << genID[i] << ") with " << genNMo[i] << " mothers" << endl;
	    for(int j=0;j<genNMo[i]&&j<2;j++) {
	      int idx = -1;
	      const GenParticle* mom = static_cast<const GenParticle*> (p.mother(j));
	      const GenParticle* gmom = static_cast<const GenParticle*> (mom->mother());
	      if (mom == NULL) break;
	      for (unsigned int k = 0; k < genParticles->size() && k<i; ++k) {
		const reco::GenParticle& testm = (*genParticles)[k];
		const GenParticle* testgm = static_cast<const GenParticle*> (testm.mother());
		if (testm.pt() == mom->pt()) {
		  if(gmom == NULL || (testgm!=NULL) && (gmom->pt() == testgm->pt())) {
		    idx = k;
		    if(blabalot) cout << "     Found a hit for a mother for index " << i << "! It's index " << idx << " (pdgid " << genID[idx] << ")" << endl;
		    break;
		  }
		}
	      }
	      if(j==0) genMo1Index[i]=idx;
	      if(j==1) genMo2Index[i]=idx;
	    }
	    if(Promptness[i]==-1&&genMo1Index[i]>=0&&Promptness[genMo1Index[i]]>-1) Promptness[i]=Promptness[genMo1Index[i]]+1;
	    if(blabalot) cout << "                    Promtpness: " << Promptness[i] << endl;
	    if(blabalot) cout << "          Mother 1: " << genMo1Index[i] << " with Promptness " << Promptness[genMo1Index[i]] << endl;
	    if(blabalot) cout << "          Mother 2: " << genMo2Index[i] << endl;
	  }
	  
	  // STEP 2: Loop from end to start flipping storeflag when necessary
	  float genPtThreshold=5.0; // any particle with less pt than this will not be stored.
	  for(int i=nGenParticles-1;i>=0;i--) {
	    if(genStatus[i]!=1) continue;
	    if(genPt[i]<genPtThreshold) continue;
	    FlipGenStoreFlag(i,Promptness,genMo1Index,genMo2Index,StoreFlag);
	  }
	  
	  // Intermediate step: Make sure that all first particles are stored, and that the earliest particles are stored (i.e. promptness criteria are met)
	  for(int i=nGenParticles-1;i>=0;i--) {
	    if(Promptness[i]<4&&Promptness[i]>=0) StoreFlag[i]=true;
	    if(i<20) StoreFlag[i]=true;
	  }
	  
	  // STEP 3: Loop again, setting new index ("nindex") and replacing Mo1 by if(Mo1>-1) Mo1=nindex[index];, same for Mo2. If storeflag, store it directly.
	  fTnGenParticles=-1;
	  for(int i=0;i<nGenParticles;i++) {
	    if(!StoreFlag[i]) continue;
	    fTnGenParticles++;
	    if(fTnGenParticles>nStoredGenParticles) {
	      edm::LogWarning("NTP") << "@SUB=analyze()"
	      << "Maximum number of gen particles exceeded";
	      fTflagmaxgenpartexc = 1;
	      break;
	    }
	    
	    genNIndex[i]=fTnGenParticles;
	    
	    //store everything
	    fTgenInfoId[fTnGenParticles] = genID[i];
	    fTgenInfoStatus[fTnGenParticles] = genStatus[i];
	    fTgenInfoNMo[fTnGenParticles] = genNMo[i];
	    fTgenInfoMo1Pt[fTnGenParticles]=0;
	    fTgenInfoMo2Pt[fTnGenParticles]=0;
	    if(genMo1Index[i]>=0) {
	      fTgenInfoMo1Pt[fTnGenParticles]=genPt[genMo1Index[i]];
	      fTgenInfoMo1[fTnGenParticles]=genNIndex[genMo1Index[i]];
	    } else {
	      fTgenInfoMo1[fTnGenParticles]=-1;
	    }
	    if(genMo2Index[i]>=0) {
	      fTgenInfoMo2Pt[fTnGenParticles]=genPt[genMo2Index[i]];
	      fTgenInfoMo2[fTnGenParticles]=genNIndex[genMo2Index[i]];
	    } else {
	      fTgenInfoMo2[fTnGenParticles]=-1;
	    }
	    
	    fTPromptnessLevel[fTnGenParticles]=Promptness[i];
	    fTgenInfoPt[fTnGenParticles] = genPt[i];
	    fTgenInfoEta[fTnGenParticles] = genEta[i];
	    fTgenInfoPhi[fTnGenParticles] = genPhi[i];
	    fTgenInfoPx[fTnGenParticles] = genPx[i];
	    fTgenInfoPy[fTnGenParticles] = genPy[i];
	    fTgenInfoPz[fTnGenParticles] = genPz[i];
	    fTgenInfoM[fTnGenParticles] = genM[i];

	    if(blabalot) {
	      cout << "Working particle " << i << " with Promptness: " << Promptness[i] << "  (storage number: " << fTnGenParticles << ")" << endl;
	      cout << "    Mother is Particle " << genMo1Index[i] << " with Promptness " << Promptness[genMo1Index[i]] << endl;
	      cout << "Particle " << fTnGenParticles << "  (" << i << "): The particle has ID = " << genID[i]                     << " and its mother has index " << genNIndex[genMo1Index[i]]     << " mo pt : " << genPt[genMo1Index[i]] << endl;
	      cout << "stored:  " << fTnGenParticles << "  (" << i << ")                      " << fTgenInfoId[fTnGenParticles] << "                          " << fTgenInfoMo1[fTnGenParticles] << "         " << fTgenInfoMo1Pt[fTnGenParticles] << endl;
	    }
	    
	  }
	  
	  fTnGenParticles++;
	  if(blabalot) cout << "A total of " << fTnGenParticles << " Particles  have been stored out of " << nGenParticles << " ( " << 100*fTnGenParticles/(float)nGenParticles << " %)" << endl;
	  
	  
	  float gluinomass=0;
	  float ngluinomass=0;
	  float chi2mass=0;
	  float nchi2mass=0;
	  float chi1mass=0;
	  float nchi1mass=0;
	  
	  for(size_t i = 0; i < genParticles->size(); ++ i) {
	    const GenParticle & p = (*genParticles)[i];
	    int id = p.pdgId();
	    double mass = p.mass();
	    if(id==1000021) {
	      gluinomass+=mass;
	      ngluinomass++;
	    }
	    if(id==1000023) {
	      chi2mass+=mass;
	      nchi2mass++;
	    }
	    if(id==1000022) {
	      chi1mass+=mass;
	      nchi1mass++;
	    }
	  }
	  
	  if(ngluinomass>0&&nchi2mass>0&&nchi1mass>0) {
	    fTxSMS = (chi2mass/nchi1mass - chi1mass/nchi1mass) / (gluinomass/ngluinomass - chi1mass/nchi1mass);
	    fTxbarSMS = 1 - fTxSMS; // Mariarosaria's definition of x
	  }
	}// end of bloat with gen information



	////////////////////////////////////////////////////////////////////////////////
	// Fill Tree ///////////////////////////////////////////////////////////////////
	fEventTree->Fill();
	fNFillTree++;
}

// Method called once each job just before starting event loop
void NTupleProducer::beginJob(){ //336 beginJob(const edm::EventSetup&)
	fNTotEvents = 0;
	fNFillTree  = 0;
	fFirstevent = true;

	fRunTree->Branch("Run"            ,&fRTrunnumber,      "Run/I");
	fRunTree->Branch("ExtXSecLO"      ,&fRTextxslo,        "ExtXSecLO/F");
	fRunTree->Branch("ExtXSecNLO"     ,&fRTextxsnlo,       "ExtXSecNLO/F");
	fRunTree->Branch("IntXSec"        ,&fRTintxs,          "IntXSec/F");
	fRunTree->Branch("MinMuPt"        ,&fRTMinmupt,        "MinMuPt/F");
	fRunTree->Branch("MaxMuEta"       ,&fRTMaxmueta,       "MaxMuEta/F");
	fRunTree->Branch("MinElPt"        ,&fRTMinelpt,        "MinElPt/F");
	fRunTree->Branch("MaxElEta"       ,&fRTMaxeleta,       "MaxElEta/F");
	fRunTree->Branch("MinJPt"         ,&fRTMinjpt,         "MinJPt/F");
	fRunTree->Branch("MinRawJPt"      ,&fRTMinrawjpt,      "MinRawJPt/F");
	fRunTree->Branch("MaxJEta"        ,&fRTMaxjeta,        "MaxJEta/F");
	fRunTree->Branch("MinJEMfrac"     ,&fRTMinjemfrac,     "MinJEMfrac/F");

	fRunTree->Branch("MinTrkPt"       ,&fRTMintrkpt,       "MinTrkPt/F");
	fRunTree->Branch("MaxTrkEta"      ,&fRTMaxtrketa,      "MaxTrkEta/F");
	fRunTree->Branch("MaxTrkNChi2"    ,&fRTMaxtrknchi2,    "MaxTrkNChi2/F");
	fRunTree->Branch("MinTrkNHits"    ,&fRTMintrknhits,    "MinTrkNHits/I");

	fRunTree->Branch("MinPhotonPt"    ,&fRTMinphopt,       "MinPhotonPt/F");
	fRunTree->Branch("MaxPhotonEta"   ,&fRTMaxphoeta,      "MaxPhotonEta/F");

	fRunTree->Branch("MaxNMus"        ,&fRTmaxnmu,         "MaxTmu/I");
	fRunTree->Branch("MaxNEles"       ,&fRTmaxnel,         "MaxTel/I");
	fRunTree->Branch("MaxNJets"       ,&fRTmaxnjet,        "MaxTjet/I");
	fRunTree->Branch("MaxNTrks"       ,&fRTmaxntrk,        "MaxTtrk/I");
	fRunTree->Branch("MaxNPhotons"    ,&fRTmaxnphot,       "MaxTphot/I");
	fRunTree->Branch("HLTLabels"      ,&fTHltLabels );
	fRunTree->Branch("HLTNames"       ,&fTHLTmenu );
	fRunTree->Branch("L1PhysMenu"     ,&fTL1physmenu );
	
	fRunTree->Branch("PileUpData"     ,&fTPileUpHistoData);
	fRunTree->Branch("PileUpMC"       ,&fTPileUpHistoMC);

	// Event information:
	fEventTree->Branch("Run"              ,&fTrunnumber       ,"Run/I");
	fEventTree->Branch("Event"            ,&fTeventnumber     ,"Event/L");
	fEventTree->Branch("LumiSection"      ,&fTlumisection     ,"LumiSection/I");
	fEventTree->Branch("PtHat"            ,&fTpthat           ,"PtHat/F");
	fEventTree->Branch("QCDPartonicHT"    ,&fTqcdPartonicHT   ,"QCDPartonicHT/F");
	fEventTree->Branch("SigProcID"        ,&fTsigprocid       ,"SigProcID/I");
	fEventTree->Branch("PDFScalePDF"      ,&fTpdfscalePDF     ,"PDFScalePDF/F");
	fEventTree->Branch("PDFID1"           ,&fTpdfid1          ,"PDFID1/I");
	fEventTree->Branch("PDFID2"           ,&fTpdfid2          ,"PDFID2/I");
	fEventTree->Branch("PDFx1"            ,&fTpdfx1           ,"PDFx1/F");
	fEventTree->Branch("PDFx2"            ,&fTpdfx2           ,"PDFx2/F");
	fEventTree->Branch("PDFxPDF1"         ,&fTpdfxPDF1        ,"PDFxPDF1/F");
	fEventTree->Branch("PDFxPDF2"         ,&fTpdfxPDF2        ,"PDFxPDF2/F");
	fEventTree->Branch("GenWeight"        ,&fTgenweight       ,"GenWeight/F");
	fEventTree->Branch("ExtXSecLO"        ,&fTextxslo         ,"ExtXSecLO/F");
	fEventTree->Branch("IntXSec"          ,&fTintxs           ,"IntXSec/F");
	fEventTree->Branch("pdfW"             ,&fTpdfW             ,"pdfW[100]/F");
	fEventTree->Branch("pdfWsum"          ,&fTpdfWsum          ,"pdfWsum/F");
	fEventTree->Branch("NPdfs"            ,&NPdfs              ,"NPdfs/I");

	// Pile-Up information:
	fEventTree->Branch("PUnumInteractions",   &fTpuNumInteractions   ,"PUnumInteractions/I");
	fEventTree->Branch("PUnumTrueInteractions",   &fTpuNumTrueInteractions   ,"PUnumTrueInteractions/I");
	fEventTree->Branch("PUnumFilled",&fTpuNumFilled,"PUnumFilled/I");
	fEventTree->Branch("PUOOTnumInteractionsEarly",&fTpuOOTNumInteractionsEarly,"PUOOTnumInteractionsEarly/I");
	fEventTree->Branch("PUOOTnumInteractionsLate",&fTpuOOTNumInteractionsLate,"PUOOTnumInteractionsLate/I");
	fEventTree->Branch("PUzPositions"     ,&fTpuZpositions     ,"PUzPositions[PUnumFilled]/F");
	fEventTree->Branch("PUsumPtLowPt"     ,&fTpuSumpT_lowpT    ,"PUsumPtLowPt[PUnumFilled]/F");
	fEventTree->Branch("PUsumPtHighPt"    ,&fTpuSumpT_highpT   ,"PUsumPtHighPt[PUnumFilled]/F");
	fEventTree->Branch("PUnTrksLowPt"     ,&fTpuNtrks_lowpT    ,"PUnTrksLowPt[PUnumFilled]/F");
	fEventTree->Branch("PUnTrksHighPt"    ,&fTpuNtrks_highpT   ,"PUnTrksHighPt[PUnumFilled]/F");
	fEventTree->Branch("Rho"              ,&fTrho              ,"Rho/F");
	fEventTree->Branch("Sigma"              ,&fTsigma              ,"Sigma/F");
	fEventTree->Branch("RhoPFnoPU"        ,&fTrhoPFnoPU        ,"RhoPFnoPU/F");
	// fEventTree->Branch("PUinstLumi"       ,&fTpuInstLumi       ,"PUinstLumi[PUnumFilled]/F");
	fEventTree->Branch("Weight"           ,&fTweight          ,"Weight/F");
	fEventTree->Branch("HLTResults"       ,&fTHLTres          , Form("HLTResults[%d]/I", gMaxhltbits));
	fEventTree->Branch("HLTPrescale"      ,&fTHLTprescale     , Form("HLTPrescale[%d]/I", gMaxhltbits));
	fEventTree->Branch("L1PhysResults"    ,&fTL1physres       , Form("L1PhysResults[%d]/I", gMaxl1physbits));
	fEventTree->Branch("L1TechResults"    ,&fTL1techres       , Form("L1TechResults[%d]/I", gMaxl1techbits));
	fEventTree->Branch("NHLTObjs"         ,&fTNHLTobjects     ,"NHLTObjs/I");
	fEventTree->Branch("HLTObjectID"      ,fTHLTObjectID[0]   , Form("HLTObjectID[%d][%d]/I",  fTNpaths, gMaxhltnobjs));
	fEventTree->Branch("HLTObjectPt"      ,fTHLTObjectPt[0]   , Form("HLTObjectPt[%d][%d]/F",  fTNpaths, gMaxhltnobjs));
	fEventTree->Branch("HLTObjectEta"     ,fTHLTObjectEta[0]  , Form("HLTObjectEta[%d][%d]/F", fTNpaths, gMaxhltnobjs));
	fEventTree->Branch("HLTObjectPhi"     ,fTHLTObjectPhi[0]  , Form("HLTObjectPhi[%d][%d]/F", fTNpaths, gMaxhltnobjs));
	fEventTree->Branch("PUWeightTotal"    ,&fTpuWeightTotal   , "PUWeightTotal/F");
	fEventTree->Branch("PUWeightInTime"   ,&fTpuWeightInTime  , "PUWeightInTime/F");
        fEventTree->Branch("MassGlu"          ,&fTMassGlu         ,"MassGlu/F");
        fEventTree->Branch("MassChi"          ,&fTMassChi         ,"MassChi/F");
        fEventTree->Branch("MassLSP"          ,&fTMassLSP         ,"MassLSP/F");
        fEventTree->Branch("xSMS"             ,&fTxSMS            ,"xSMS/F");
        fEventTree->Branch("xbarSMS"          ,&fTxbarSMS         ,"xbarSMS/F");
        fEventTree->Branch("M0"               ,&fTSUSYScanM0      ,"M0/F");
        fEventTree->Branch("M12"              ,&fTSUSYScanM12     ,"M12/F");
        fEventTree->Branch("signMu"           ,&fTSUSYScanMu      ,"signMu/F");
        fEventTree->Branch("A0"               ,&fTSUSYScanA0      ,"A0/F");
        fEventTree->Branch("process"          ,&fTprocess         ,"process/I");
//        fEventTree->Branch("nStoredGeneratorParticles",(int *) &nStoredGenParticles,"nStoredGeneratorParticles/I");


//        fEventTree->Branch("FlagMaxGenPartExceeded",&fTflagmaxgenpartexc,"FlagMaxGenPartExceeded/I");
//        fEventTree->Branch("nGenParticles",&fTnGenParticles,"nGenParticles/I");
//        fEventTree->Branch("genInfoId"        ,&fTgenInfoId         ,"genInfoId[nGenParticles]/I");
//        fEventTree->Branch("genInfoStatus"    ,&fTgenInfoStatus     ,"genInfoStatus[nGenParticles]/I");
//        fEventTree->Branch("genInfoMass"      ,&fTgenInfoMass       ,"genInfoMass[nGenParticles]/F");
//        fEventTree->Branch("genInfoNMo"       ,&fTgenInfoNMo        ,"genInfoNMo[nGenParticles]/I");
//        fEventTree->Branch("genInfoMo1Pt"     ,&fTgenInfoMo1Pt      ,"genInfoMo1Pt[nGenParticles]/F");
//        fEventTree->Branch("genInfoMo2Pt"     ,&fTgenInfoMo2Pt      ,"genInfoMo2Pt[nGenParticles]/F");
//        fEventTree->Branch("genInfoNDa"       ,&fTgenInfoNDa        ,"genInfoNDa[nGenParticles]/I");
//        fEventTree->Branch("genInfoMo1"       ,&fTgenInfoMo1        ,"genInfoMo1[nGenParticles]/I");
//        fEventTree->Branch("genInfoMo2"       ,&fTgenInfoMo2        ,"genInfoMo2[nGenParticles]/I");
//        fEventTree->Branch("genInfoDa1"       ,&fTgenInfoDa1        ,"genInfoDa1[nGenParticles]/I");
//        fEventTree->Branch("genInfoDa2"       ,&fTgenInfoDa2        ,"genInfoDa2[nGenParticles]/I");
//        fEventTree->Branch("genInfoPt"        ,&fTgenInfoPt         ,"genInfoPt[nGenParticles]/F");
//        fEventTree->Branch("genInfoEta"       ,&fTgenInfoEta        ,"genInfoEta[nGenParticles]/F");
//        fEventTree->Branch("genInfoPhi"       ,&fTgenInfoPhi        ,"genInfoPhi[nGenParticles]/F");
//        fEventTree->Branch("genInfoPx"        ,&fTgenInfoPx         ,"genInfoPx[nGenParticles]/F");
//        fEventTree->Branch("genInfoPy"        ,&fTgenInfoPy         ,"genInfoPy[nGenParticles]/F");
//        fEventTree->Branch("genInfoPz"        ,&fTgenInfoPz         ,"genInfoPz[nGenParticles]/F");
//        fEventTree->Branch("genInfoM"         ,&fTgenInfoM          ,"genInfoM[nGenParticles]/F");
//	fEventTree->Branch("genInfoMoIndex"   ,&fTgenInfoMoIndex    ,"genInfoMoIndex[nGenParticles]/I");
//	fEventTree->Branch("PromptnessLevel"  ,&fTPromptnessLevel   ,"PromptnessLevel[nGenParticles]/I");


	fEventTree->Branch("PrimVtxGood"      ,&fTgoodvtx           ,"PrimVtxGood/I");
	fEventTree->Branch("PrimVtxx"         ,&fTprimvtxx          ,"PrimVtxx/F");
	fEventTree->Branch("PrimVtxy"         ,&fTprimvtxy          ,"PrimVtxy/F");
	fEventTree->Branch("PrimVtxz"         ,&fTprimvtxz          ,"PrimVtxz/F");
	fEventTree->Branch("PrimVtxRho"       ,&fTprimvtxrho        ,"PrimVtxRho/F");
	fEventTree->Branch("PrimVtxxE"        ,&fTprimvtxxE         ,"PrimVtxxE/F");
	fEventTree->Branch("PrimVtxyE"        ,&fTprimvtxyE         ,"PrimVtxyE/F");
	fEventTree->Branch("PrimVtxzE"        ,&fTprimvtxzE         ,"PrimVtxzE/F");
	fEventTree->Branch("PrimVtxNChi2"     ,&fTpvtxznchi2        ,"PrimVtxNChi2/F");
	fEventTree->Branch("PrimVtxNdof"      ,&fTpvtxndof          ,"PrimVtxNdof/F");
	fEventTree->Branch("PrimVtxIsFake"    ,&fTpvtxisfake        ,"PrimVtxIsFake/I");
	fEventTree->Branch("PrimVtxPtSum"     ,&fTpvtxptsum         ,"PrimVtxPtSum/F");
	fEventTree->Branch("Beamspotx"        ,&fTbeamspotx         ,"Beamspotx/F");
	fEventTree->Branch("Beamspoty"        ,&fTbeamspoty         ,"Beamspoty/F");
	fEventTree->Branch("Beamspotz"        ,&fTbeamspotz         ,"Beamspotz/F");
	fEventTree->Branch("NCaloTowers"      ,&fTNCaloTowers       ,"NCaloTowers/I");
	fEventTree->Branch("GoodEvent"        ,&fTgoodevent         ,"GoodEvent/I");
	fEventTree->Branch("MaxMuExceed"      ,&fTflagmaxmuexc      ,"MaxMuExceed/I");
	fEventTree->Branch("MaxElExceed"      ,&fTflagmaxelexc      ,"MaxElExceed/I");
	fEventTree->Branch("MaxJetExceed"     ,&fTflagmaxjetexc     ,"MaxJetExceed/I");
	fEventTree->Branch("MaxUncJetExceed"  ,&fTflagmaxujetexc    ,"MaxUncJetExceed/I");
	fEventTree->Branch("MaxTrkExceed"     ,&fTflagmaxtrkexc     ,"MaxTrkExceed/I");
	fEventTree->Branch("MaxPhotonsExceed" ,&fTflagmaxphoexc     ,"MaxPhotonsExceed/I");
	fEventTree->Branch("MaxGenLepExceed"  ,&fTflagmaxgenleptexc ,"MaxGenLepExceed/I");
	fEventTree->Branch("MaxGenPhoExceed"  ,&fTflagmaxgenphotexc ,"MaxGenPhoExceed/I");
	fEventTree->Branch("MaxGenJetExceed"  ,&fTflagmaxgenjetexc  ,"MaxGenJetExceed/I");
	fEventTree->Branch("MaxVerticesExceed",&fTflagmaxvrtxexc    ,"MaxVerticesExceed/I");
	fEventTree->Branch("HBHENoiseFlag"    ,&fTHBHENoiseFlag     ,"HBHENoiseFlag/I");
	fEventTree->Branch("HBHENoiseFlagIso" ,&fTHBHENoiseFlagIso  ,"HBHENoiseFlagIso/I");
	fEventTree->Branch("CSCTightHaloID"   ,&fTcscTightHaloID    ,"CSCTightHaloID/I");
	fEventTree->Branch("EcalDeadTPFilterFlag"        ,&fTecalDeadTPFilterFlag        ,"EcalDeadTPFilterFlag/I");
	fEventTree->Branch("RecovRecHitFilterFlag"       ,&fRecovRecHitFilterFlag        ,"RecovRecHitFilterFlag/I");
	fEventTree->Branch("RA2TrackingFailureFilterFlag",&fTra2TrackingFailureFilterFlag,"RA2TrackingFailureFilterFlag/I");
	//FR fEventTree->Branch("PBNRFlag"         ,&fPBNRFlag           ,"PBNRFlag/I");
	// fEventTree->Branch("EcalDeadCellBEFlag",&fTEcalDeadCellBEFlag,"EcalDeadCellBEFlag/I");
	// fEventTree->Branch("NECALGapClusters"  ,&fTnECALGapClusters  ,"NECALGapClusters/I");
	// fEventTree->Branch("EcalGapBE"         ,&fTEcalGapBE         ,"EcalGapBE[NECALGapClusters]/F");
	// fEventTree->Branch("EcalGapClusterSize",&fTEcalGapClusterSize,"EcalGapClusterSize[NECALGapClusters]/I");

	// Gen-Leptons
	fEventTree->Branch("NGenLeptons"      ,&fTngenleptons         ,"NGenLeptons/I");
	fEventTree->Branch("GenLeptonID"      ,&fTGenLeptonId         ,"GenLeptonID[NGenLeptons]/I");
	fEventTree->Branch("GenLeptonPt"      ,&fTGenLeptonPt         ,"GenLeptonPt[NGenLeptons]/F");
	fEventTree->Branch("GenLeptonEta"     ,&fTGenLeptonEta        ,"GenLeptonEta[NGenLeptons]/F");
	fEventTree->Branch("GenLeptonPhi"     ,&fTGenLeptonPhi        ,"GenLeptonPhi[NGenLeptons]/F");
	fEventTree->Branch("GenLeptonMID"     ,&fTGenLeptonMId        ,"GenLeptonMID[NGenLeptons]/I");
	fEventTree->Branch("GenLeptonMStatus" ,&fTGenLeptonMStatus    ,"GenLeptonMStatus[NGenLeptons]/I");
	fEventTree->Branch("GenLeptonMPt"     ,&fTGenLeptonMPt        ,"GenLeptonMPt[NGenLeptons]/F");
	fEventTree->Branch("GenLeptonMEta"    ,&fTGenLeptonMEta       ,"GenLeptonMEta[NGenLeptons]/F");
	fEventTree->Branch("GenLeptonMPhi"    ,&fTGenLeptonMPhi       ,"GenLeptonMPhi[NGenLeptons]/F");
	fEventTree->Branch("GenLeptonGMID"    ,&fTGenLeptonGMId       ,"GenLeptonGMID[NGenLeptons]/I");
	fEventTree->Branch("GenLeptonGMStatus",&fTGenLeptonGMStatus   ,"GenLeptonGMStatus[NGenLeptons]/I");
	fEventTree->Branch("GenLeptonGMPt"    ,&fTGenLeptonGMPt       ,"GenLeptonGMPt[NGenLeptons]/F");
	fEventTree->Branch("GenLeptonGMEta"   ,&fTGenLeptonGMEta      ,"GenLeptonGMEta[NGenLeptons]/F");
	fEventTree->Branch("GenLeptonGMPhi"   ,&fTGenLeptonGMPhi      ,"GenLeptonGMPhi[NGenLeptons]/F");

	// Gen-Photons
	fEventTree->Branch("NGenPhotons"      ,&fTngenphotons         ,"NGenPhotons/I");
	fEventTree->Branch("GenPhotonPt"      ,&fTGenPhotonPt         ,"GenPhotonPt[NGenPhotons]/F");
	fEventTree->Branch("GenPhotonEta"     ,&fTGenPhotonEta        ,"GenPhotonEta[NGenPhotons]/F");
	fEventTree->Branch("GenPhotonPhi"     ,&fTGenPhotonPhi        ,"GenPhotonPhi[NGenPhotons]/F");
	fEventTree->Branch("GenPhotonPartonMindR",&fTGenPhotonPartonMindR ,"GenPhotonPartonMindR[NGenPhotons]/F");
	fEventTree->Branch("GenPhotonMotherID"     ,&fTGenPhotonMotherID        ,"GenPhotonMotherID[NGenPhotons]/I");
	fEventTree->Branch("GenPhotonMotherStatus"     ,&fTGenPhotonMotherStatus        ,"GenPhotonMotherStatus[NGenPhotons]/I");
	fEventTree->Branch("GenPhotonIsoDR03"     ,&fTGenPhotonIsoDR03        ,"GenPhotonIsoDR03[NGenPhotons]/F");
	fEventTree->Branch("GenPhotonIsoDR04"     ,&fTGenPhotonIsoDR04        ,"GenPhotonIsoDR04[NGenPhotons]/F");

	// Gen-Jets
	fEventTree->Branch("NGenJets"    ,&fTNGenJets   ,"NGenJets/I");
	fEventTree->Branch("GenJetPt"    ,&fTGenJetPt   ,"GenJetPt[NGenJets]/F");
	fEventTree->Branch("GenJetEta"   ,&fTGenJetEta  ,"GenJetEta[NGenJets]/F");
	fEventTree->Branch("GenJetPhi"   ,&fTGenJetPhi  ,"GenJetPhi[NGenJets]/F");
	fEventTree->Branch("GenJetE"     ,&fTGenJetE    ,"GenJetE[NGenJets]/F");
	fEventTree->Branch("GenJetEmE"   ,&fTGenJetemE  ,"GenJetEmE[NGenJets]/F");
	fEventTree->Branch("GenJetHadE"  ,&fTGenJethadE ,"GenJetHadE[NGenJets]/F");
	fEventTree->Branch("GenJetInvE"  ,&fTGenJetinvE ,"GenJetInvE[NGenJets]/F");

	// Vertices:
	fEventTree->Branch("NVrtx",            &fTnvrtx           ,"NVrtx/I");
	fEventTree->Branch("VrtxX",            &fTvrtxx           ,"VrtxX[NVrtx]/F");
	fEventTree->Branch("VrtxY",            &fTvrtxy           ,"VrtxY[NVrtx]/F");
	fEventTree->Branch("VrtxZ",            &fTvrtxz           ,"VrtxZ[NVrtx]/F");
	fEventTree->Branch("VrtxXE",           &fTvrtxxE          ,"VrtxXE[NVrtx]/F");
	fEventTree->Branch("VrtxYE",           &fTvrtxyE          ,"VrtxYE[NVrtx]/F");
	fEventTree->Branch("VrtxZE",           &fTvrtxzE          ,"VrtxZE[NVrtx]/F");
	fEventTree->Branch("VrtxNdof",         &fTvrtxndof        ,"VrtxNdof[NVrtx]/F");
	fEventTree->Branch("VrtxChi2",         &fTvrtxchi2        ,"VrtxChi2[NVrtx]/F");
	fEventTree->Branch("VrtxNtrks",        &fTvrtxntrks       ,"VrtxNtrks[NVrtx]/F");
	fEventTree->Branch("VrtxSumPt",        &fTvrtxsumpt       ,"VrtxSumPt[NVrtx]/F");
	fEventTree->Branch("VrtxIsFake",       &fTvrtxisfake      ,"VrtxIsFake[NVrtx]/I");

	// Muons:
	fEventTree->Branch("NMus"             ,&fTnmu              ,"NMus/I");
	fEventTree->Branch("NMusTot"          ,&fTnmutot           ,"NMusTot/I");
	fEventTree->Branch("NGMus"            ,&fTnglobalmu        ,"NGMus/I");
	fEventTree->Branch("NTMus"            ,&fTntrackermu       ,"NTMus/I");
	fEventTree->Branch("MuGood"           ,&fTgoodmu           ,"MuGood[NMus]/I");
	fEventTree->Branch("MuIsIso"          ,&fTmuIsIso          ,"MuIsIso[NMus]/I");
	fEventTree->Branch("MuIsGlobalMuon"   ,&fTmuIsGM           ,"MuIsGlobalMuon[NMus]/I");
	fEventTree->Branch("MuIsTrackerMuon"  ,&fTmuIsTM           ,"MuIsTrackerMuon[NMus]/I");
	fEventTree->Branch("MuPx"             ,&fTmupx             ,"MuPx[NMus]/F");
	fEventTree->Branch("MuPy"             ,&fTmupy             ,"MuPy[NMus]/F");
	fEventTree->Branch("MuPz"             ,&fTmupz             ,"MuPz[NMus]/F");
	fEventTree->Branch("MuPt"             ,&fTmupt             ,"MuPt[NMus]/F");
	fEventTree->Branch("MuInnerTkPt"      ,&fTmuinnerpt        ,"MuInnerTkPt[NMus]/F");
	fEventTree->Branch("MuPtE"            ,&fTmuptE            ,"MuPtE[NMus]/F");
	fEventTree->Branch("MuE"              ,&fTmue              ,"MuE[NMus]/F");
	fEventTree->Branch("MuEt"             ,&fTmuet             ,"MuEt[NMus]/F");
	fEventTree->Branch("MuEta"            ,&fTmueta            ,"MuEta[NMus]/F");
	fEventTree->Branch("MuPhi"            ,&fTmuphi            ,"MuPhi[NMus]/F");
	fEventTree->Branch("MuCharge"         ,&fTmucharge         ,"MuCharge[NMus]/I");
	fEventTree->Branch("MuRelIso03"       ,&fTmuiso            ,"MuRelIso03[NMus]/F");
	fEventTree->Branch("MuIso03SumPt"     ,&fTmuIso03sumPt     ,"MuIso03SumPt[NMus]/F");
	fEventTree->Branch("MuIso03EmEt"      ,&fTmuIso03emEt      ,"MuIso03EmEt[NMus]/F");
	fEventTree->Branch("MuIso03HadEt"     ,&fTmuIso03hadEt     ,"MuIso03HadEt[NMus]/F");
	fEventTree->Branch("MuIso03EMVetoEt"  ,&fTmuIso03emVetoEt  ,"MuIso03EMVetoEt[NMus]/F");
	fEventTree->Branch("MuIso03HadVetoEt" ,&fTmuIso03hadVetoEt ,"MuIso03HadVetoEt[NMus]/F");
	fEventTree->Branch("MuIso05SumPt"     ,&fTmuIso05sumPt     ,"MuIso05SumPt[NMus]/F");
	fEventTree->Branch("MuIso05EmEt"      ,&fTmuIso05emEt      ,"MuIso05EmEt[NMus]/F");
	fEventTree->Branch("MuIso05HadEt"     ,&fTmuIso05hadEt     ,"MuIso05HadEt[NMus]/F");
	fEventTree->Branch("MuEem"            ,&fTmueecal          ,"MuEem[NMus]/F");
	fEventTree->Branch("MuEhad"           ,&fTmuehcal          ,"MuEhad[NMus]/F");
	fEventTree->Branch("MuD0BS"           ,&fTmud0bs           ,"MuD0BS[NMus]/F");
	fEventTree->Branch("MuD0PV"           ,&fTmud0pv           ,"MuD0PV[NMus]/F");
	fEventTree->Branch("MuD0E"            ,&fTmud0E            ,"MuD0E[NMus]/F");
	fEventTree->Branch("MuDzBS"           ,&fTmudzbs           ,"MuDzBS[NMus]/F");
	fEventTree->Branch("MuDzPV"           ,&fTmudzpv           ,"MuDzPV[NMus]/F");
	fEventTree->Branch("MuDzE"            ,&fTmudzE            ,"MuDzE[NMus]/F");
	fEventTree->Branch("MuNChi2"          ,&fTmunchi2          ,"MuNChi2[NMus]/F");
	fEventTree->Branch("MuNGlHits"        ,&fTmunglhits        ,"MuNGlHits[NMus]/I");
	fEventTree->Branch("MuNMuHits"        ,&fTmunmuhits        ,"MuNMuHits[NMus]/I");
	fEventTree->Branch("MuNTkHits"        ,&fTmuntkhits        ,"MuNTkHits[NMus]/I");
	fEventTree->Branch("MuNPxHits"        ,&fTmunpxhits        ,"MuNPxHits[NMus]/I");
	fEventTree->Branch("MuInnerTkNChi2"   ,&fTmuinntknchi2     ,"MuInnerTkNChi2[NMus]/F");
	fEventTree->Branch("MuNMatches"       ,&fTmunmatches       ,"MuNMatches[NMus]/I");
	fEventTree->Branch("MuNChambers"      ,&fTmunchambers      ,"MuNChambers[NMus]/I");
	fEventTree->Branch("MuCaloComp"       ,&fTmucalocomp       ,"MuCaloComp[NMus]/F");
	fEventTree->Branch("MuSegmComp"       ,&fTmusegmcomp       ,"MuSegmComp[NMus]/F");

	fEventTree->Branch("MuIsGMPT"                  ,&fTmuIsGMPT                  ,"MuIsGMPT[NMus]/I");
	fEventTree->Branch("MuIsGMTkChiComp"           ,&fTmuIsGMTkChiComp           ,"MuIsGMTkChiComp[NMus]/I");
	fEventTree->Branch("MuIsGMStaChiComp"          ,&fTmuIsGMStaChiComp          ,"MuIsGMStaChiComp[NMus]/I");
	fEventTree->Branch("MuIsGMTkKinkTight"         ,&fTmuIsGMTkKinkTight         ,"MuIsGMTkKinkTight[NMus]/I");
	fEventTree->Branch("MuIsAllStaMuons"           ,&fTmuIsAllStaMuons           ,"MuIsAllStaMuons[NMus]/I");
	fEventTree->Branch("MuIsAllTrkMuons"           ,&fTmuIsAllTrkMuons           ,"MuIsAllTrkMuons[NMus]/I");
	fEventTree->Branch("MuIsTrkMuonArbitrated"     ,&fTmuIsTrkMuArb              ,"MuIsTrkMuonArbitrated[NMus]/I");
	fEventTree->Branch("MuIsAllArbitrated"         ,&fTmuIsAllArb                ,"MuIsAllArbitrated[NMus]/I");
	fEventTree->Branch("MuIsTMLSLoose"             ,&fTmuIsTMLastStationLoose    ,"MuIsTMLSLoose[NMus]/I");
	fEventTree->Branch("MuIsTMLSTight"             ,&fTmuIsTMLastStationTight    ,"MuIsTMLSTight[NMus]/I");
	fEventTree->Branch("MuIsTM2DCompLoose"         ,&fTmuIsTM2DCompLoose         ,"MuIsTM2DCompLoose[NMus]/I");
	fEventTree->Branch("MuIsTM2DCompTight"         ,&fTmuIsTM2DCompTight         ,"MuIsTM2DCompTight[NMus]/I");
	fEventTree->Branch("MuIsTMOneStationLoose"     ,&fTmuIsTMOneStationLoose     ,"MuIsTMOneStationLoose[NMus]/I");
	fEventTree->Branch("MuIsTMOneStationTight"     ,&fTmuIsTMOneStationTight     ,"MuIsTMOneStationTight[NMus]/I");
	fEventTree->Branch("MuIsTMLSOptLowPtLoose"     ,&fTmuIsTMLSOPL               ,"MuIsTMLSOptLowPtLoose[NMus]/I");
	fEventTree->Branch("MuIsTMLSAngLoose"          ,&fTmuIsTMLastStationAngLoose ,"MuIsTMLSAngLoose[NMus]/I");
	fEventTree->Branch("MuIsTMLastStationAngTight" ,&fTmuIsTMLastStationAngTight ,"MuIsTMLastStationAngTight[NMus]/I");
	fEventTree->Branch("MuIsTMOneStationAngTight"  ,&fTmuIsTMOneStationAngTight  ,"MuIsTMOneStationAngTight[NMus]/I");
	fEventTree->Branch("MuIsTMOneStationAngLoose"  ,&fTmuIsTMOneStationAngLoose  ,"MuIsTMOneStationAngLoose[NMus]/I");

	fEventTree->Branch("MuGenID"          ,&fTGenMuId         ,"MuGenID[NMus]/I");
	fEventTree->Branch("MuGenStatus"      ,&fTGenMuStatus     ,"MuGenStatus[NMus]/I");
	fEventTree->Branch("MuGenPt"          ,&fTGenMuPt         ,"MuGenPt[NMus]/F");
	fEventTree->Branch("MuGenEta"         ,&fTGenMuEta        ,"MuGenEta[NMus]/F");
	fEventTree->Branch("MuGenPhi"         ,&fTGenMuPhi        ,"MuGenPhi[NMus]/F");
	fEventTree->Branch("MuGenE"           ,&fTGenMuE          ,"MuGenE[NMus]/F");
	fEventTree->Branch("MuGenMID"         ,&fTGenMuMId        ,"MuGenMID[NMus]/I");
	fEventTree->Branch("MuGenMStatus"     ,&fTGenMuMStatus    ,"MuGenMStatus[NMus]/I");
	fEventTree->Branch("MuGenMPt"         ,&fTGenMuMPt        ,"MuGenMPt[NMus]/F");
	fEventTree->Branch("MuGenMEta"        ,&fTGenMuMEta       ,"MuGenMEta[NMus]/F");
	fEventTree->Branch("MuGenMPhi"        ,&fTGenMuMPhi       ,"MuGenMPhi[NMus]/F");
	fEventTree->Branch("MuGenME"          ,&fTGenMuME         ,"MuGenME[NMus]/F");
	fEventTree->Branch("MuGenGMID"        ,&fTGenMuGMId       ,"MuGenGMID[NMus]/I");
	fEventTree->Branch("MuGenGMStatus"    ,&fTGenMuGMStatus   ,"MuGenGMStatus[NMus]/I");
	fEventTree->Branch("MuGenGMPt"        ,&fTGenMuGMPt       ,"MuGenGMPt[NMus]/F");
	fEventTree->Branch("MuGenGMEta"       ,&fTGenMuGMEta      ,"MuGenGMEta[NMus]/F");
	fEventTree->Branch("MuGenGMPhi"       ,&fTGenMuGMPhi      ,"MuGenGMPhi[NMus]/F");
	fEventTree->Branch("MuGenGME"         ,&fTGenMuGME        ,"MuGenGME[NMus]/F");

	fEventTree->Branch("NEBhits"              ,&fTnEBhits              ,"NEBhits/I");
	fEventTree->Branch("EBrechitE"            ,&fTEBrechitE            ,"EBrechitE[NEBhits]/F");
	fEventTree->Branch("EBrechitPt"           ,&fTEBrechitPt           ,"EBrechitPt[NEBhits]/F");
	fEventTree->Branch("EBrechitEta"          ,&fTEBrechitEta          ,"EBrechitEta[NEBhits]/F");
	fEventTree->Branch("EBrechitPhi"          ,&fTEBrechitPhi          ,"EBrechitPhi[NEBhits]/F");
	fEventTree->Branch("EBrechitChi2"         ,&fTEBrechitChi2         ,"EBrechitChi2[NEBhits]/F");
	fEventTree->Branch("EBrechitTime"         ,&fTEBrechitTime         ,"EBrechitTime[NEBhits]/F");
	fEventTree->Branch("EBrechitE4oE1"        ,&fTEBrechitE4oE1        ,"EBrechitE4oE1[NEBhits]/F");
	fEventTree->Branch("EBrechitE2oE9"        ,&fTEBrechitE2oE9        ,"EBrechitE2oE9[NEBhits]/F");

	// Electrons:
	fEventTree->Branch("NEles"                       ,&fTneles            ,"NEles/I");
	fEventTree->Branch("NElesTot"                    ,&fTnelestot         ,"NElesTot/I");
	fEventTree->Branch("ElGood"                      ,&fTgoodel           ,"ElGood[NEles]/I");
	fEventTree->Branch("ElIsIso"                     ,&fTeIsIso           ,"ElIsIso[NEles]/I");
	fEventTree->Branch("ElChargeMisIDProb"           ,&fTeChargeMisIDProb ,"ElChargeMisIDProb[NEles]/I");
	fEventTree->Branch("ElPx"                        ,&fTepx              ,"ElPx[NEles]/F");
	fEventTree->Branch("ElPy"                        ,&fTepy              ,"ElPy[NEles]/F");
	fEventTree->Branch("ElPz"                        ,&fTepz              ,"ElPz[NEles]/F");
	fEventTree->Branch("ElPt"                        ,&fTept              ,"ElPt[NEles]/F");
	fEventTree->Branch("ElPtE"                       ,&fTeptE             ,"ElPtE[NEles]/F");
	fEventTree->Branch("ElE"                         ,&fTee               ,"ElE[NEles]/F");
	fEventTree->Branch("ElEt"                        ,&fTeet              ,"ElEt[NEles]/F");
	fEventTree->Branch("ElEta"                       ,&fTeeta             ,"ElEta[NEles]/F");
	fEventTree->Branch("ElTheta"                     ,&fTetheta           ,"ElTheta[NEles]/F");
	fEventTree->Branch("ElSCEta"                     ,&fTesceta           ,"ElSCEta[NEles]/F");
	fEventTree->Branch("ElPhi"                       ,&fTephi             ,"ElPhi[NEles]/F");
	fEventTree->Branch("ElGsfTkPt"                   ,&fTegsfpt           ,"ElGsfTkPt[NEles]/F");
	fEventTree->Branch("ElGsfTkEta"                  ,&fTegsfeta          ,"ElGsfTkEta[NEles]/F");
	fEventTree->Branch("ElGsfTkPhi"                  ,&fTegsfphi          ,"ElGsfTkPhi[NEles]/F");
	fEventTree->Branch("ElTrkMomentumError"          ,&fTetrkmomerror     ,"ElTrkMomentumError[NEles]/F");
	fEventTree->Branch("ElEcalEnergyError"           ,&fTeecalergerror    ,"ElEcalEnergyError[NEles]/F");
	fEventTree->Branch("ElEleMomentumError"          ,&fTeelemomerror     ,"ElEleMomentumError[NEles]/F");
	fEventTree->Branch("ElNBrems"                    ,&fTenbrems          ,"ElNBrems[NEles]/I");
	fEventTree->Branch("ElD0BS"                      ,&fTed0bs            ,"ElD0BS[NEles]/F");
	fEventTree->Branch("ElD0PV"                      ,&fTed0pv            ,"ElD0PV[NEles]/F");
	fEventTree->Branch("ElD0E"                       ,&fTed0E             ,"ElD0E[NEles]/F");
	fEventTree->Branch("ElDzBS"                      ,&fTedzbs            ,"ElDzBS[NEles]/F");
	fEventTree->Branch("ElDzPV"                      ,&fTedzpv            ,"ElDzPV[NEles]/F");
	fEventTree->Branch("ElDzE"                       ,&fTedzE             ,"ElDzE[NEles]/F");
	fEventTree->Branch("ElRelIso03"                  ,&fTeiso03           ,"ElRelIso03[NEles]/F");
	fEventTree->Branch("ElRelIso04"                  ,&fTeiso04           ,"ElRelIso04[NEles]/F");
	fEventTree->Branch("ElDR03TkSumPt"               ,&fTdr03tksumpt      ,"ElDR03TkSumPt[NEles]/F");
	fEventTree->Branch("ElDR04TkSumPt"               ,&fTdr04tksumpt      ,"ElDR04TkSumPt[NEles]/F");
	fEventTree->Branch("ElDR03EcalRecHitSumEt"       ,&fTdr03ecalrechitsumet     ,"ElDR03EcalRecHitSumEt[NEles]/F");
	fEventTree->Branch("ElDR04EcalRecHitSumEt"       ,&fTdr04ecalrechitsumet     ,"ElDR04EcalRecHitSumEt[NEles]/F");
	fEventTree->Branch("ElDR03HcalTowerSumEt"        ,&fTdr03hcaltowersumet      ,"ElDR03HcalTowerSumEt[NEles]/F");
	fEventTree->Branch("ElDR04HcalTowerSumEt"        ,&fTdr04hcaltowersumet      ,"ElDR04HcalTowerSumEt[NEles]/F");
	fEventTree->Branch("ElNChi2"                     ,&fTenchi2                  ,"ElNChi2[NEles]/F");
	fEventTree->Branch("ElCharge"                    ,&fTecharge                 ,"ElCharge[NEles]/I");
	fEventTree->Branch("ElCInfoIsGsfCtfCons"         ,&fTeCInfoIsGsfCtfCons      ,"ElCInfoIsGsfCtfCons[NEles]/I");
	fEventTree->Branch("ElCInfoIsGsfCtfScPixCons"    ,&fTeCInfoIsGsfCtfScPixCons ,"ElCInfoIsGsfCtfScPixCons[NEles]/I");
	fEventTree->Branch("ElCInfoIsGsfScPixCons"       ,&fTeCInfoIsGsfScPixCons    ,"ElCInfoIsGsfScPixCons[NEles]/I");
	fEventTree->Branch("ElScPixCharge"               ,&fTeCInfoScPixCharge       ,"ElScPixCharge[NEles]/I");
	fEventTree->Branch("ElClosestCtfTrackPt"         ,&fTeClosestCtfTrackpt      ,"ElClosestCtfTrackPt[NEles]/F");
	fEventTree->Branch("ElClosestCtfTrackEta"        ,&fTeClosestCtfTracketa     ,"ElClosestCtfTrackEta[NEles]/F");
	fEventTree->Branch("ElClosestCtfTrackPhi"        ,&fTeClosestCtfTrackphi     ,"ElClosestCtfTrackPhi[NEles]/F");
	fEventTree->Branch("ElClosestCtfTrackCharge"     ,&fTeClosestCtfTrackcharge  ,"ElClosestCtfTrackCharge[NEles]/I");
	fEventTree->Branch("ElIDMva"                     ,&fTeIDMva             ,"ElIDMva[NEles]/F");
	fEventTree->Branch("ElIDTight"                   ,&fTeIDTight           ,"ElIDTight[NEles]/I");
	fEventTree->Branch("ElIDLoose"                   ,&fTeIDLoose           ,"ElIDLoose[NEles]/I");
	fEventTree->Branch("ElIDRobustTight"             ,&fTeIDRobustTight     ,"ElIDRobustTight[NEles]/I");
	fEventTree->Branch("ElIDRobustLoose"             ,&fTeIDRobustLoose     ,"ElIDRobustLoose[NEles]/I");
	fEventTree->Branch("ElIDsimpleWPrelIso"          ,&fTeIDsimpleWPrelIso    ,"ElIDsimpleWPrelIso[NEles]/I");
	fEventTree->Branch("ElIDsimpleWP80relIso"        ,&fTeIDsimpleWP80relIso  ,"ElIDsimpleWP80relIso[NEles]/I");
	fEventTree->Branch("ElIDsimpleWP85relIso"        ,&fTeIDsimpleWP85relIso  ,"ElIDsimpleWP85relIso[NEles]/I");
	fEventTree->Branch("ElIDsimpleWP90relIso"        ,&fTeIDsimpleWP90relIso  ,"ElIDsimpleWP90relIso[NEles]/I");
	fEventTree->Branch("ElIDsimpleWP95relIso"        ,&fTeIDsimpleWP95relIso  ,"ElIDsimpleWP95relIso[NEles]/I");
	fEventTree->Branch("ElInGap"                     ,&fTeInGap             ,"ElInGap[NEles]/I");
	fEventTree->Branch("ElEcalDriven"                ,&fTeEcalDriven        ,"ElEcalDriven[NEles]/I");
	fEventTree->Branch("ElTrackerDriven"             ,&fTeTrackerDriven     ,"ElTrackerDriven[NEles]/I");
	fEventTree->Branch("ElBasicClustersSize"         ,&fTeBasicClustersSize ,"ElBasicClustersSize[NEles]/I");
	fEventTree->Branch("Elfbrem"                     ,&fTefbrem             ,"Elfbrem[NEles]/F");
	fEventTree->Branch("ElHcalOverEcal"              ,&fTeHcalOverEcal      ,"ElHcalOverEcal[NEles]/F");
	fEventTree->Branch("ElE1x5"                      ,&fTeE1x5              ,"ElE1x5[NEles]/F");
	fEventTree->Branch("ElE5x5"                      ,&fTeE5x5              ,"ElE5x5[NEles]/F");
	fEventTree->Branch("ElE2x5Max"                   ,&fTeE2x5Max           ,"ElE2x5Max[NEles]/F");
	fEventTree->Branch("ElSigmaIetaIeta"             ,&fTeSigmaIetaIeta     ,"ElSigmaIetaIeta[NEles]/F");
	fEventTree->Branch("ElDeltaPhiSeedClusterAtCalo" ,&fTeDeltaPhiSeedClusterAtCalo ,"ElDeltaPhiSeedClusterAtCalo[NEles]/F");
	fEventTree->Branch("ElDeltaEtaSeedClusterAtCalo" ,&fTeDeltaEtaSeedClusterAtCalo ,"ElDeltaEtaSeedClusterAtCalo[NEles]/F");
	fEventTree->Branch("ElDeltaPhiSuperClusterAtVtx" ,&fTeDeltaPhiSuperClusterAtVtx ,"ElDeltaPhiSuperClusterAtVtx[NEles]/F");
	fEventTree->Branch("ElDeltaEtaSuperClusterAtVtx" ,&fTeDeltaEtaSuperClusterAtVtx ,"ElDeltaEtaSuperClusterAtVtx[NEles]/F");
	fEventTree->Branch("ElCaloEnergy"                ,&fTecaloenergy        ,"ElCaloEnergy[NEles]/F");
	fEventTree->Branch("ElTrkMomAtVtx"               ,&fTetrkmomatvtx        ,"ElTrkMomAtVtx[NEles]/F");
	fEventTree->Branch("ElESuperClusterOverP"        ,&fTeESuperClusterOverP        ,"ElESuperClusterOverP[NEles]/F");
	fEventTree->Branch("ElNumberOfMissingInnerHits"  ,&fTeNumberOfMissingInnerHits  ,"ElNumberOfMissingInnerHits[NEles]/I");
	fEventTree->Branch("ElSCindex",&fTElSCindex,"ElSCindex[NEles]/I");

	// fEventTree->Branch("ElIsInJet"                   ,&fTeIsInJet           ,"ElIsInJet[NEles]/I");
	// fEventTree->Branch("ElSharedPx"                  ,&fTeSharedPx          ,"ElSharedPx[NEles]/F");
	// fEventTree->Branch("ElSharedPy"                  ,&fTeSharedPy          ,"ElSharedPy[NEles]/F");
	// fEventTree->Branch("ElSharedPz"                  ,&fTeSharedPz          ,"ElSharedPz[NEles]/F");
	// fEventTree->Branch("ElSharedEnergy"              ,&fTeSharedEnergy      ,"ElSharedEnergy[NEles]/F");
	// fEventTree->Branch("ElDuplicateEl"               ,&fTeDupEl             ,"ElDuplicateEl[NEles]/I");
	fEventTree->Branch("ElConvPartnerTrkDist"        ,&fTeConvPartTrackDist   ,"ElConvPartnerTrkDist[NEles]/F");
	fEventTree->Branch("ElConvPartnerTrkDCot"        ,&fTeConvPartTrackDCot   ,"ElConvPartnerTrkDCot[NEles]/F");
	fEventTree->Branch("ElConvPartnerTrkPt"          ,&fTeConvPartTrackPt     ,"ElConvPartnerTrkPt[NEles]/F");
	fEventTree->Branch("ElConvPartnerTrkEta"         ,&fTeConvPartTrackEta    ,"ElConvPartnerTrkEta[NEles]/F");
	fEventTree->Branch("ElConvPartnerTrkPhi"         ,&fTeConvPartTrackPhi    ,"ElConvPartnerTrkPhi[NEles]/F");
	fEventTree->Branch("ElConvPartnerTrkCharge"      ,&fTeConvPartTrackCharge ,"ElConvPartnerTrkCharge[NEles]/F");
	fEventTree->Branch("ElScSeedSeverity"            ,&fTeScSeedSeverity      ,"ElScSeedSeverity[NEles]/I");
	fEventTree->Branch("ElE1OverE9"                  ,&fTeE1OverE9            ,"ElE1OverE9[NEles]/F");
	fEventTree->Branch("ElS4OverS1"                  ,&fTeS4OverS1            ,"ElS4OverS1[NEles]/F");

	fEventTree->Branch("ElGenID"                     ,&fTGenElId         ,"ElGenID[NEles]/I");
	fEventTree->Branch("ElGenStatus"                 ,&fTGenElStatus     ,"ElGenStatus[NEles]/I");
	fEventTree->Branch("ElGenPt"                     ,&fTGenElPt         ,"ElGenPt[NEles]/F");
	fEventTree->Branch("ElGenEta"                    ,&fTGenElEta        ,"ElGenEta[NEles]/F");
	fEventTree->Branch("ElGenPhi"                    ,&fTGenElPhi        ,"ElGenPhi[NEles]/F");
	fEventTree->Branch("ElGenE"                      ,&fTGenElE          ,"ElGenE[NEles]/F");
	fEventTree->Branch("ElGenMID"                    ,&fTGenElMId        ,"ElGenMID[NEles]/I");
	fEventTree->Branch("ElGenMStatus"                ,&fTGenElMStatus    ,"ElGenMStatus[NEles]/I");
	fEventTree->Branch("ElGenMPt"                    ,&fTGenElMPt        ,"ElGenMPt[NEles]/F");
	fEventTree->Branch("ElGenMEta"                   ,&fTGenElMEta       ,"ElGenMEta[NEles]/F");
	fEventTree->Branch("ElGenMPhi"                   ,&fTGenElMPhi       ,"ElGenMPhi[NEles]/F");
	fEventTree->Branch("ElGenME"                     ,&fTGenElME         ,"ElGenME[NEles]/F");
	fEventTree->Branch("ElGenGMID"                   ,&fTGenElGMId       ,"ElGenGMID[NEles]/I");
	fEventTree->Branch("ElGenGMStatus"               ,&fTGenElGMStatus   ,"ElGenGMStatus[NEles]/I");
	fEventTree->Branch("ElGenGMPt"                   ,&fTGenElGMPt       ,"ElGenGMPt[NEles]/F");
	fEventTree->Branch("ElGenGMEta"                  ,&fTGenElGMEta      ,"ElGenGMEta[NEles]/F");
	fEventTree->Branch("ElGenGMPhi"                  ,&fTGenElGMPhi      ,"ElGenGMPhi[NEles]/F");
	fEventTree->Branch("ElGenGME"                    ,&fTGenElGME        ,"ElGenGME[NEles]/F");

	//PfCandidates:
	
	fEventTree->Branch("NPfCand"                     ,&fTNPfCand         ,"NPfCand/I");
	fEventTree->Branch("PfCandPdgId"                 ,&fTPfCandPdgId     ,"PfCandPdgId[NPfCand]/F");
	fEventTree->Branch("PfCandPt"                    ,&fTPfCandPt        ,"PfCandPt[NPfCand]/F");
	fEventTree->Branch("PfCandEta"                   ,&fTPfCandEta       ,"PfCandEta[NPfCand]/F");
	fEventTree->Branch("PfCandPhi"                   ,&fTPfCandPhi       ,"PfCandPhi[NPfCand]/F");
//	fEventTree->Branch("PfCandPx"                    ,&fTPfCandPx        ,"PfCandPx[NPfCand]/F");
//	fEventTree->Branch("PfCandPy"                    ,&fTPfCandPy        ,"PfCandPy[NPfCand]/F");
//	fEventTree->Branch("PfCandPz"                    ,&fTPfCandPz        ,"PfCandPz[NPfCand]/F");
	fEventTree->Branch("PfCandEnergy"                ,&fTPfCandEnergy    ,"PfCandEnergy[NPfCand]/F");
	fEventTree->Branch("PfCandVx"                    ,&fTPfCandVx        ,"PfCandVx[NPfCand]/F");
	fEventTree->Branch("PfCandVy"                    ,&fTPfCandVy        ,"PfCandVy[NPfCand]/F");
	fEventTree->Branch("PfCandVz"                    ,&fTPfCandVz        ,"PfCandVz[NPfCand]/F");
	fEventTree->Branch("PfCandBelongsToJet"          ,&fTPfCandBelongsToJet,"PfCandBelongsToJet[NPfCand]/I");
//	fEventTree->Branch("PfCandMomX"                  ,&fTPfCandMomX      ,"PfCandMomX[NPfCand]/F");
//	fEventTree->Branch("PfCandMomY"                  ,&fTPfCandMomY      ,"PfCandMomY[NPfCand]/F");
//	fEventTree->Branch("PfCandMomZ"                  ,&fTPfCandMomZ      ,"PfCandMomZ[NPfCand]/F");
//	fEventTree->Branch("PfCandHasHitInFirstPixelLayer", &fTPfCandHasHitInFirstPixelLayer, "PfCandHasHitInFirstPixelLayer[NPfCand]/I");
//	fEventTree->Branch("PfCandTrackRefPx", &fTPfCandTrackRefPx, "PfCandTrackRefPx[NPfCand]/F");
//	fEventTree->Branch("PfCandTrackRefPy", &fTPfCandTrackRefPy, "PfCandTrackRefPy[NPfCand]/F");
//	fEventTree->Branch("PfCandTrackRefPz", &fTPfCandTrackRefPz, "PfCandTrackRefPz[NPfCand]/F");
//	fEventTree->Branch("PfCandTrackRefVx", &fTPfCandTrackRefVx, "PfCandTrackRefVx[NPfCand]/F");
//	fEventTree->Branch("PfCandTrackRefVy", &fTPfCandTrackRefVy, "PfCandTrackRefVy[NPfCand]/F");
//	fEventTree->Branch("PfCandTrackRefVz", &fTPfCandTrackRefVz, "PfCandTrackRefVz[NPfCand]/F");

	// Photons:
	fEventTree->Branch("NPhotons"         ,&fTnphotons          ,"NPhotons/I");
	fEventTree->Branch("NPhotonsTot"      ,&fTnphotonstot       ,"NPhotonsTot/I");
	fEventTree->Branch("PhoGood"          ,&fTgoodphoton        ,"PhoGood[NPhotons]/I");
	fEventTree->Branch("PhoIsIso"         ,&fTPhotIsIso         ,"PhoIsIso[NPhotons]/I");
	fEventTree->Branch("PhoPt"            ,&fTPhotPt            ,"PhoPt[NPhotons]/F");
	fEventTree->Branch("PhoPx"            ,&fTPhotPx            ,"PhoPx[NPhotons]/F");
	fEventTree->Branch("PhoPy"            ,&fTPhotPy            ,"PhoPy[NPhotons]/F");
	fEventTree->Branch("PhoPz"            ,&fTPhotPz            ,"PhoPz[NPhotons]/F");
	fEventTree->Branch("PhoEta"           ,&fTPhotEta           ,"PhoEta[NPhotons]/F");
	fEventTree->Branch("PhoPhi"           ,&fTPhotPhi           ,"PhoPhi[NPhotons]/F");
	fEventTree->Branch("PhoEnergy"        ,&fTPhotEnergy        ,"PhoEnergy[NPhotons]/F");
	fEventTree->Branch("PhoIso03Ecal"     ,&fTPhotIso03Ecal     ,"PhoIso03Ecal[NPhotons]/F");
	fEventTree->Branch("PhoIso03Hcal"     ,&fTPhotIso03Hcal     ,"PhoIso03Hcal[NPhotons]/F");
	fEventTree->Branch("PhoIso03TrkSolid" ,&fTPhotIso03TrkSolid ,"PhoIso03TrkSolid[NPhotons]/F");
	fEventTree->Branch("PhoIso03TrkHollow",&fTPhotIso03TrkHollow,"PhoIso03TrkHollow[NPhotons]/F");
	fEventTree->Branch("PhoIso03"         ,&fTPhotIso03         ,"PhoIso03[NPhotons]/F");
	fEventTree->Branch("PhoIso04Ecal"     ,&fTPhotIso04Ecal     ,"PhoIso04Ecal[NPhotons]/F");
	fEventTree->Branch("PhoIso04Hcal"     ,&fTPhotIso04Hcal     ,"PhoIso04Hcal[NPhotons]/F");
	fEventTree->Branch("PhoIso04TrkSolid" ,&fTPhotIso04TrkSolid ,"PhoIso04TrkSolid[NPhotons]/F");
	fEventTree->Branch("PhoIso04TrkHollow",&fTPhotIso04TrkHollow,"PhoIso04TrkHollow[NPhotons]/F");
	fEventTree->Branch("PhoIso04"         ,&fTPhotIso04         ,"PhoIso04[NPhotons]/F");
	fEventTree->Branch("PhoR9"            ,&fTPhotR9            ,"PhoR9[NPhotons]/F");
	fEventTree->Branch("PhoCaloPositionX" ,&fTPhotcaloPosX      ,"PhoCaloPositionX[NPhotons]/F");
	fEventTree->Branch("PhoCaloPositionY" ,&fTPhotcaloPosY      ,"PhoCaloPositionY[NPhotons]/F");
	fEventTree->Branch("PhoCaloPositionZ" ,&fTPhotcaloPosZ      ,"PhoCaloPositionZ[NPhotons]/F");
	fEventTree->Branch("PhoHoverE"        ,&fTPhotHoverE        ,"PhoHoverE[NPhotons]/F");
	fEventTree->Branch("PhoH1overE"       ,&fTPhotH1overE       ,"PhoH1overE[NPhotons]/F");
	fEventTree->Branch("PhoH2overE"       ,&fTPhotH2overE       ,"PhoH2overE[NPhotons]/F");
	fEventTree->Branch("PhoSigmaIetaIeta" ,&fTPhotSigmaIetaIeta ,"PhoSigmaIetaIeta[NPhotons]/F");
	fEventTree->Branch("PhoSCRawEnergy"   ,&fTPhotSCEnergy      ,"PhoSCRawEnergy[NPhotons]/F");
	fEventTree->Branch("PhoSCEtaWidth"    ,&fTPhotSCEtaWidth    ,"PhoSCEtaWidth[NPhotons]/F");
	fEventTree->Branch("PhoSCSigmaPhiPhi" ,&fTPhotSCSigmaPhiPhi ,"PhoSCSigmaPhiPhi[NPhotons]/F");
	fEventTree->Branch("PhoHasPixSeed"    ,&fTPhotHasPixSeed    ,"PhoHasPixSeed[NPhotons]/I");
	fEventTree->Branch("PhoPassConvSafeElectronVeto"    ,&fTPhotPassConvSafeElectronVeto    ,"PhoPassConvSafeElectronVeto[NPhotons]/I");
	fEventTree->Branch("PhoHasConvTrks"   ,&fTPhotHasConvTrks   ,"PhoHasConvTrks[NPhotons]/I");
	fEventTree->Branch("PhoScSeedSeverity",&fTPhotScSeedSeverity,"PhoScSeedSeverity[NPhotons]/I");
	fEventTree->Branch("PhoE1OverE9"      ,&fTPhotE1OverE9      ,"PhoE1OverE9[NPhotons]/F");
	fEventTree->Branch("PhoS4OverS1"      ,&fTPhotS4OverS1      ,"PhoS4OverS1[NPhotons]/F");
       fEventTree->Branch("PhoSigmaEtaEta"   ,&fTPhotSigmaEtaEta   ,"PhoSigmaEtaEta[NPhotons]/F");
       fEventTree->Branch("PhoE1x5"   ,&fTPhote1x5   ,"PhoE1x5[NPhotons]/F");
       fEventTree->Branch("PhoE2x5"   ,&fTPhote2x5   ,"PhoE2x5[NPhotons]/F");
       fEventTree->Branch("PhoE3x3"   ,&fTPhote3x3   ,"PhoE3x3[NPhotons]/F");
       fEventTree->Branch("PhoE5x5"   ,&fTPhote5x5   ,"PhoE5x5[NPhotons]/F");
       fEventTree->Branch("PhomaxEnergyXtal"   ,&fTPhotmaxEnergyXtal   ,"PhomaxEnergyXtal[NPhotons]/F");
       fEventTree->Branch("PhoIso03HcalDepth1"   ,&fTPhotIso03HcalDepth1   ,"PhoIso03HcalDepth1[NPhotons]/F");
       fEventTree->Branch("PhoIso03HcalDepth2"   ,&fTPhotIso03HcalDepth2   ,"PhoIso03HcalDepth2[NPhotons]/F");
       fEventTree->Branch("PhoIso04HcalDepth1"   ,&fTPhotIso04HcalDepth1   ,"PhoIso04HcalDepth1[NPhotons]/F");
       fEventTree->Branch("PhoIso04HcalDepth2"   ,&fTPhotIso04HcalDepth2   ,"PhoIso04HcalDepth2[NPhotons]/F");
       fEventTree->Branch("PhoIso03nTrksSolid"   ,&fTPhotIso03nTrksSolid   ,"PhoIso03nTrksSolid[NPhotons]/I");
       fEventTree->Branch("PhoIso03nTrksHollow"   ,&fTPhotIso03nTrksHollow   ,"PhoIso03nTrksHollow[NPhotons]/I");
       fEventTree->Branch("PhoIso04nTrksSolid"   ,&fTPhotIso04nTrksSolid   ,"PhoIso04nTrksSolid[NPhotons]/I");
       fEventTree->Branch("PhoIso04nTrksHollow"   ,&fTPhotIso04nTrksHollow   ,"PhoIso04nTrksHollow[NPhotons]/I");
       fEventTree->Branch("PhoisEB"   ,&fTPhotisEB   ,"PhoisEB[NPhotons]/I");
       fEventTree->Branch("PhoisEE"   ,&fTPhotisEE   ,"PhoisEE[NPhotons]/I");
       fEventTree->Branch("PhoisEBEtaGap"   ,&fTPhotisEBEtaGap   ,"PhoisEBEtaGap[NPhotons]/I");
       fEventTree->Branch("PhoisEBPhiGap"   ,&fTPhotisEBPhiGap   ,"PhoisEBPhiGap[NPhotons]/I");
       fEventTree->Branch("PhoisEERingGap"   ,&fTPhotisEERingGap   ,"PhoisEERingGap[NPhotons]/I");
       fEventTree->Branch("PhoisEEDeeGap"   ,&fTPhotisEEDeeGap   ,"PhoisEEDeeGap[NPhotons]/I");
       fEventTree->Branch("PhoisEBEEGap"   ,&fTPhotisEBEEGap   ,"PhoisEBEEGap[NPhotons]/I");
       fEventTree->Branch("PhoisPFlowPhoton"   ,&fTPhotisPFlowPhoton   ,"PhoisPFlowPhoton[NPhotons]/I");
       fEventTree->Branch("PhoisStandardPhoton"   ,&fTPhotisStandardPhoton   ,"PhoisStandardPhoton[NPhotons]/I"); 
       fEventTree->Branch("PhoMCmatchindex"   ,&fTPhotMCmatchindex   ,"PhoMCmatchindex[NPhotons]/I");
       fEventTree->Branch("PhoMCmatchexitcode"   ,&fTPhotMCmatchexitcode   ,"PhoMCmatchexitcode[NPhotons]/I");
       fEventTree->Branch("Pho_ChargedHadronIso",&fT_pho_ChargedHadronIso,"Pho_ChargedHadronIso[NPhotons]/F");
       fEventTree->Branch("Pho_NeutralHadronIso",&fT_pho_NeutralHadronIso,"Pho_NeutralHadronIso[NPhotons]/F");
       fEventTree->Branch("Pho_PhotonIso",&fT_pho_PhotonIso,"Pho_PhotonIso[NPhotons]/F");
       fEventTree->Branch("PhoSCRemovalPFIsoCharged",&fTPhoSCRemovalPFIsoCharged,"PhoSCRemovalPFIsoCharged[NPhotons]/F");
       fEventTree->Branch("PhoSCRemovalPFIsoChargedPrimVtx",&fTPhoSCRemovalPFIsoChargedPrimVtx,"PhoSCRemovalPFIsoChargedPrimVtx[NPhotons]/F");
       fEventTree->Branch("PhoSCRemovalPFIsoNeutral",&fTPhoSCRemovalPFIsoNeutral,"PhoSCRemovalPFIsoNeutral[NPhotons]/F");
       fEventTree->Branch("PhoSCRemovalPFIsoPhoton",&fTPhoSCRemovalPFIsoPhoton,"PhoSCRemovalPFIsoPhoton[NPhotons]/F");
       fEventTree->Branch("PhoSCRemovalPFIsoCharged_RCone",&fTPhoSCRemovalPFIsoCharged_RCone,"PhoSCRemovalPFIsoCharged_RCone[NPhotons]/F");
       fEventTree->Branch("PhoSCRemovalPFIsoChargedPrimVtx_RCone",&fTPhoSCRemovalPFIsoChargedPrimVtx_RCone,"PhoSCRemovalPFIsoChargedPrimVtx_RCone[NPhotons]/F");
       fEventTree->Branch("PhoSCRemovalPFIsoNeutral_RCone",&fTPhoSCRemovalPFIsoNeutral_RCone,"PhoSCRemovalPFIsoNeutral_RCone[NPhotons]/F");
       fEventTree->Branch("PhoSCRemovalPFIsoPhoton_RCone",&fTPhoSCRemovalPFIsoPhoton_RCone,"PhoSCRemovalPFIsoPhoton_RCone[NPhotons]/F");
       fEventTree->Branch("fTPhoSCRemoval_RCone_Eta",&fTPhoSCRemoval_RCone_Eta,"fTPhoSCRemoval_RCone_Eta[NPhotons]/F");
       fEventTree->Branch("fTPhoSCRemoval_RCone_Phi",&fTPhoSCRemoval_RCone_Phi,"fTPhoSCRemoval_RCone_Phi[NPhotons]/F");
       fEventTree->Branch("Pho_isPFPhoton",&fT_pho_isPFPhoton,"Pho_isPFPhoton[NPhotons]/I");
       fEventTree->Branch("Pho_isPFElectron",&fT_pho_isPFElectron,"Pho_isPFElectron[NPhotons]/I");
       fEventTree->Branch("PhotSCindex",&fTPhotSCindex,"PhotSCindex[NPhotons]/I");
       fEventTree->Branch("pho_matchedPFPhotonCand",&fT_pho_matchedPFPhotonCand,"pho_matchedPFPhotonCand[NPhotons]/I");
       fEventTree->Branch("pho_matchedPFElectronCand",&fT_pho_matchedPFElectronCand,"pho_matchedPFElectronCand[NPhotons]/I");
       fEventTree->Branch("PhoVx",&fTPhotVx,"PhoVx[NPhotons]/F");
       fEventTree->Branch("PhoVy",&fTPhotVy,"PhoVy[NPhotons]/F");
       fEventTree->Branch("PhoVz",&fTPhotVz,"PhoVz[NPhotons]/F");

//       fEventTree->Branch("pho_Cone01PhotonIso_dEta015EB_dR070EE_mvVtx",&fT_pho_Cone01PhotonIso_dEta015EB_dR070EE_mvVtx,"pho_Cone01PhotonIso_dEta015EB_dR070EE_mvVtx[NPhotons]/F");
//       fEventTree->Branch("pho_Cone02PhotonIso_dEta015EB_dR070EE_mvVtx",&fT_pho_Cone02PhotonIso_dEta015EB_dR070EE_mvVtx,"pho_Cone02PhotonIso_dEta015EB_dR070EE_mvVtx[NPhotons]/F");
//       fEventTree->Branch("pho_Cone03PhotonIso_dEta015EB_dR070EE_mvVtx",&fT_pho_Cone03PhotonIso_dEta015EB_dR070EE_mvVtx,"pho_Cone03PhotonIso_dEta015EB_dR070EE_mvVtx[NPhotons]/F");
//       fEventTree->Branch("pho_Cone04PhotonIso_dEta015EB_dR070EE_mvVtx",&fT_pho_Cone04PhotonIso_dEta015EB_dR070EE_mvVtx,"pho_Cone04PhotonIso_dEta015EB_dR070EE_mvVtx[NPhotons]/F");
//
//       fEventTree->Branch("pho_Cone01NeutralHadronIso_mvVtx",&fT_pho_Cone01NeutralHadronIso_mvVtx,"pho_Cone01NeutralHadronIso_mvVtx[NPhotons]/F");
//       fEventTree->Branch("pho_Cone02NeutralHadronIso_mvVtx",&fT_pho_Cone02NeutralHadronIso_mvVtx,"pho_Cone02NeutralHadronIso_mvVtx[NPhotons]/F");
//       fEventTree->Branch("pho_Cone03NeutralHadronIso_mvVtx",&fT_pho_Cone03NeutralHadronIso_mvVtx,"pho_Cone03NeutralHadronIso_mvVtx[NPhotons]/F");
//       fEventTree->Branch("pho_Cone04NeutralHadronIso_mvVtx",&fT_pho_Cone04NeutralHadronIso_mvVtx,"pho_Cone04NeutralHadronIso_mvVtx[NPhotons]/F");
//
//       fEventTree->Branch("pho_Cone01ChargedHadronIso_dR02_dz02_dxy01",&fT_pho_Cone01ChargedHadronIso_dR02_dz02_dxy01,"pho_Cone01ChargedHadronIso_dR02_dz02_dxy01[NPhotons]/F");
//       fEventTree->Branch("pho_Cone02ChargedHadronIso_dR02_dz02_dxy01",&fT_pho_Cone02ChargedHadronIso_dR02_dz02_dxy01,"pho_Cone02ChargedHadronIso_dR02_dz02_dxy01[NPhotons]/F");
//       fEventTree->Branch("pho_Cone03ChargedHadronIso_dR02_dz02_dxy01",&fT_pho_Cone03ChargedHadronIso_dR02_dz02_dxy01,"pho_Cone03ChargedHadronIso_dR02_dz02_dxy01[NPhotons]/F");
//       fEventTree->Branch("pho_Cone04ChargedHadronIso_dR02_dz02_dxy01",&fT_pho_Cone04ChargedHadronIso_dR02_dz02_dxy01,"pho_Cone04ChargedHadronIso_dR02_dz02_dxy01[NPhotons]/F");
//
//       fEventTree->Branch("pho_Cone03PFCombinedIso",&fT_pho_Cone03PFCombinedIso,"pho_Cone03CombinedIso[NPhotons]/F");
//       fEventTree->Branch("pho_Cone04PFCombinedIso",&fT_pho_Cone04PFCombinedIso,"pho_Cone04CombinedIso[NPhotons]/F");


       /*
fEventTree->Branch("Pho_Cone04PhotonIso_dR0_dEta0_pt0",&fT_pho_Cone04PhotonIso_dR0_dEta0_pt0,"Pho_Cone04PhotonIso_dR0_dEta0_pt0[NPhotons]/F");
fEventTree->Branch("Pho_Cone04PhotonIso_dR0_dEta0_pt5",&fT_pho_Cone04PhotonIso_dR0_dEta0_pt5,"Pho_Cone04PhotonIso_dR0_dEta0_pt5[NPhotons]/F");
fEventTree->Branch("Pho_Cone04PhotonIso_dR8_dEta0_pt0",&fT_pho_Cone04PhotonIso_dR8_dEta0_pt0,"Pho_Cone04PhotonIso_dR8_dEta0_pt0[NPhotons]/F");
fEventTree->Branch("Pho_Cone04PhotonIso_dR8_dEta0_pt5",&fT_pho_Cone04PhotonIso_dR8_dEta0_pt5,"Pho_Cone04PhotonIso_dR8_dEta0_pt5[NPhotons]/F");
fEventTree->Branch("Pho_Cone01PhotonIso_dR045EB070EE_dEta015_pt08EB1EE_mvVtx",&fT_pho_Cone01PhotonIso_dR045EB070EE_dEta015_pt08EB1EE_mvVtx,"Pho_Cone01PhotonIso_dR045EB070EE_dEta015_pt08EB1EE_mvVtx[NPhotons]/F");
fEventTree->Branch("Pho_Cone02PhotonIso_dR045EB070EE_dEta015_pt08EB1EE_mvVtx",&fT_pho_Cone02PhotonIso_dR045EB070EE_dEta015_pt08EB1EE_mvVtx,"Pho_Cone02PhotonIso_dR045EB070EE_dEta015_pt08EB1EE_mvVtx[NPhotons]/F");
fEventTree->Branch("Pho_Cone03PhotonIso_dR045EB070EE_dEta015_pt08EB1EE_mvVtx",&fT_pho_Cone03PhotonIso_dR045EB070EE_dEta015_pt08EB1EE_mvVtx,"Pho_Cone03PhotonIso_dR045EB070EE_dEta015_pt08EB1EE_mvVtx[NPhotons]/F");
fEventTree->Branch("Pho_Cone04PhotonIso_dR045EB070EE_dEta015_pt08EB1EE_mvVtx",&fT_pho_Cone04PhotonIso_dR045EB070EE_dEta015_pt08EB1EE_mvVtx,"Pho_Cone04PhotonIso_dR045EB070EE_dEta015_pt08EB1EE_mvVtx[NPhotons]/F");
fEventTree->Branch("Pho_Cone04NeutralHadronIso_dR0_dEta0_pt0",&fT_pho_Cone04NeutralHadronIso_dR0_dEta0_pt0,"Pho_Cone04NeutralHadronIso_dR0_dEta0_pt0[NPhotons]/F");
fEventTree->Branch("Pho_Cone04NeutralHadronIso_dR0_dEta0_pt5",&fT_pho_Cone04NeutralHadronIso_dR0_dEta0_pt5,"Pho_Cone04NeutralHadronIso_dR0_dEta0_pt5[NPhotons]/F");
fEventTree->Branch("Pho_Cone04NeutralHadronIso_dR0_dEta0_pt0_nocracks",&fT_pho_Cone04NeutralHadronIso_dR0_dEta0_pt0_nocracks,"Pho_Cone04NeutralHadronIso_dR0_dEta0_pt0_nocracks[NPhotons]/F");
fEventTree->Branch("Pho_Cone04NeutralHadronIso_dR0_dEta0_pt5_nocracks",&fT_pho_Cone04NeutralHadronIso_dR0_dEta0_pt5_nocracks,"Pho_Cone04NeutralHadronIso_dR0_dEta0_pt5_nocracks[NPhotons]/F");
fEventTree->Branch("Pho_Cone04NeutralHadronIso_dR7_dEta0_pt0",&fT_pho_Cone04NeutralHadronIso_dR7_dEta0_pt0,"Pho_Cone04NeutralHadronIso_dR7_dEta0_pt0[NPhotons]/F");
fEventTree->Branch("Pho_Cone04NeutralHadronIso_dR7_dEta0_pt5",&fT_pho_Cone04NeutralHadronIso_dR7_dEta0_pt5,"Pho_Cone04NeutralHadronIso_dR7_dEta0_pt5[NPhotons]/F");
fEventTree->Branch("Pho_Cone01NeutralHadronIso_dR0_dEta0_pt0_mvVtx",&fT_pho_Cone01NeutralHadronIso_dR0_dEta0_pt0_mvVtx,"Pho_Cone01NeutralHadronIso_dR0_dEta0_pt0_mvVtx[NPhotons]/F");
fEventTree->Branch("Pho_Cone02NeutralHadronIso_dR0_dEta0_pt0_mvVtx",&fT_pho_Cone02NeutralHadronIso_dR0_dEta0_pt0_mvVtx,"Pho_Cone02NeutralHadronIso_dR0_dEta0_pt0_mvVtx[NPhotons]/F");
fEventTree->Branch("Pho_Cone03NeutralHadronIso_dR0_dEta0_pt0_mvVtx",&fT_pho_Cone03NeutralHadronIso_dR0_dEta0_pt0_mvVtx,"Pho_Cone03NeutralHadronIso_dR0_dEta0_pt0_mvVtx[NPhotons]/F");
fEventTree->Branch("Pho_Cone04NeutralHadronIso_dR0_dEta0_pt0_mvVtx",&fT_pho_Cone04NeutralHadronIso_dR0_dEta0_pt0_mvVtx,"Pho_Cone04NeutralHadronIso_dR0_dEta0_pt0_mvVtx[NPhotons]/F");
fEventTree->Branch("Pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_dz0_old",&fT_pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_dz0_old,"Pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_dz0_old[NPhotons]/F");
fEventTree->Branch("Pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_PFnoPU_old",&fT_pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_PFnoPU_old,"Pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_PFnoPU_old[NPhotons]/F");
fEventTree->Branch("Pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_dz0_old",&fT_pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_dz0_old,"Pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_dz0_old[NPhotons]/F");
fEventTree->Branch("Pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_PFnoPU_old",&fT_pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_PFnoPU_old,"Pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_PFnoPU_old[NPhotons]/F");
fEventTree->Branch("Pho_Cone01ChargedHadronIso_dR0_dEta0_pt0_dz0",&fT_pho_Cone01ChargedHadronIso_dR0_dEta0_pt0_dz0,"Pho_Cone01ChargedHadronIso_dR0_dEta0_pt0_dz0[NPhotons]/F");
fEventTree->Branch("Pho_Cone01ChargedHadronIso_dR0_dEta0_pt0_dz1_dxy01",&fT_pho_Cone01ChargedHadronIso_dR0_dEta0_pt0_dz1_dxy01,"Pho_Cone01ChargedHadronIso_dR0_dEta0_pt0_dz1_dxy01[NPhotons]/F");
fEventTree->Branch("Pho_Cone01ChargedHadronIso_dR0_dEta0_pt0_PFnoPU",&fT_pho_Cone01ChargedHadronIso_dR0_dEta0_pt0_PFnoPU,"Pho_Cone01ChargedHadronIso_dR0_dEta0_pt0_PFnoPU[NPhotons]/F");
fEventTree->Branch("Pho_Cone01ChargedHadronIso_dR015_dEta0_pt0_dz0",&fT_pho_Cone01ChargedHadronIso_dR015_dEta0_pt0_dz0,"Pho_Cone01ChargedHadronIso_dR015_dEta0_pt0_dz0[NPhotons]/F");
fEventTree->Branch("Pho_Cone01ChargedHadronIso_dR015_dEta0_pt0_dz1_dxy01",&fT_pho_Cone01ChargedHadronIso_dR015_dEta0_pt0_dz1_dxy01,"Pho_Cone01ChargedHadronIso_dR015_dEta0_pt0_dz1_dxy01[NPhotons]/F");
fEventTree->Branch("Pho_Cone01ChargedHadronIso_dR015_dEta0_pt0_PFnoPU",&fT_pho_Cone01ChargedHadronIso_dR015_dEta0_pt0_PFnoPU,"Pho_Cone01ChargedHadronIso_dR015_dEta0_pt0_PFnoPU[NPhotons]/F");
fEventTree->Branch("Pho_Cone02ChargedHadronIso_dR0_dEta0_pt0_dz0",&fT_pho_Cone02ChargedHadronIso_dR0_dEta0_pt0_dz0,"Pho_Cone02ChargedHadronIso_dR0_dEta0_pt0_dz0[NPhotons]/F");
fEventTree->Branch("Pho_Cone02ChargedHadronIso_dR0_dEta0_pt0_dz1_dxy01",&fT_pho_Cone02ChargedHadronIso_dR0_dEta0_pt0_dz1_dxy01,"Pho_Cone02ChargedHadronIso_dR0_dEta0_pt0_dz1_dxy01[NPhotons]/F");
fEventTree->Branch("Pho_Cone02ChargedHadronIso_dR0_dEta0_pt0_PFnoPU",&fT_pho_Cone02ChargedHadronIso_dR0_dEta0_pt0_PFnoPU,"Pho_Cone02ChargedHadronIso_dR0_dEta0_pt0_PFnoPU[NPhotons]/F");
fEventTree->Branch("Pho_Cone02ChargedHadronIso_dR015_dEta0_pt0_dz0",&fT_pho_Cone02ChargedHadronIso_dR015_dEta0_pt0_dz0,"Pho_Cone02ChargedHadronIso_dR015_dEta0_pt0_dz0[NPhotons]/F");
fEventTree->Branch("Pho_Cone02ChargedHadronIso_dR015_dEta0_pt0_dz1_dxy01",&fT_pho_Cone02ChargedHadronIso_dR015_dEta0_pt0_dz1_dxy01,"Pho_Cone02ChargedHadronIso_dR015_dEta0_pt0_dz1_dxy01[NPhotons]/F");
fEventTree->Branch("Pho_Cone02ChargedHadronIso_dR015_dEta0_pt0_PFnoPU",&fT_pho_Cone02ChargedHadronIso_dR015_dEta0_pt0_PFnoPU,"Pho_Cone02ChargedHadronIso_dR015_dEta0_pt0_PFnoPU[NPhotons]/F");
fEventTree->Branch("Pho_Cone03ChargedHadronIso_dR0_dEta0_pt0_dz0",&fT_pho_Cone03ChargedHadronIso_dR0_dEta0_pt0_dz0,"Pho_Cone03ChargedHadronIso_dR0_dEta0_pt0_dz0[NPhotons]/F");
fEventTree->Branch("Pho_Cone03ChargedHadronIso_dR0_dEta0_pt0_dz1_dxy01",&fT_pho_Cone03ChargedHadronIso_dR0_dEta0_pt0_dz1_dxy01,"Pho_Cone03ChargedHadronIso_dR0_dEta0_pt0_dz1_dxy01[NPhotons]/F");
fEventTree->Branch("Pho_Cone03ChargedHadronIso_dR0_dEta0_pt0_PFnoPU",&fT_pho_Cone03ChargedHadronIso_dR0_dEta0_pt0_PFnoPU,"Pho_Cone03ChargedHadronIso_dR0_dEta0_pt0_PFnoPU[NPhotons]/F");
fEventTree->Branch("Pho_Cone03ChargedHadronIso_dR015_dEta0_pt0_dz0",&fT_pho_Cone03ChargedHadronIso_dR015_dEta0_pt0_dz0,"Pho_Cone03ChargedHadronIso_dR015_dEta0_pt0_dz0[NPhotons]/F");
fEventTree->Branch("Pho_Cone03ChargedHadronIso_dR015_dEta0_pt0_dz1_dxy01",&fT_pho_Cone03ChargedHadronIso_dR015_dEta0_pt0_dz1_dxy01,"Pho_Cone03ChargedHadronIso_dR015_dEta0_pt0_dz1_dxy01[NPhotons]/F");
fEventTree->Branch("Pho_Cone03ChargedHadronIso_dR015_dEta0_pt0_PFnoPU",&fT_pho_Cone03ChargedHadronIso_dR015_dEta0_pt0_PFnoPU,"Pho_Cone03ChargedHadronIso_dR015_dEta0_pt0_PFnoPU[NPhotons]/F");
fEventTree->Branch("Pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_dz0",&fT_pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_dz0,"Pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_dz0[NPhotons]/F");
fEventTree->Branch("Pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_dz1_dxy01",&fT_pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_dz1_dxy01,"Pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_dz1_dxy01[NPhotons]/F");
fEventTree->Branch("Pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_PFnoPU",&fT_pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_PFnoPU,"Pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_PFnoPU[NPhotons]/F");
fEventTree->Branch("Pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_dz0",&fT_pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_dz0,"Pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_dz0[NPhotons]/F");
fEventTree->Branch("Pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_dz1_dxy01",&fT_pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_dz1_dxy01,"Pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_dz1_dxy01[NPhotons]/F");
fEventTree->Branch("Pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_PFnoPU",&fT_pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_PFnoPU,"Pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_PFnoPU[NPhotons]/F");
       */

 fEventTree->Branch("Pho_conv_validvtx",&pho_conv_validvtx,"Pho_conv_validvtx[NPhotons]/O");
 fEventTree->Branch("Pho_conv_ntracks",&pho_conv_ntracks,"Pho_conv_ntracks[NPhotons]/I");
 fEventTree->Branch("Pho_conv_chi2_probability",&pho_conv_chi2_probability,"Pho_conv_chi2_probability[NPhotons]/F");
 fEventTree->Branch("Pho_conv_eoverp",&pho_conv_eoverp,"Pho_conv_eoverp[NPhotons]/F");

 fEventTree->Branch("Conv_n",&conv_n,"Conv_n/I");
 fEventTree->Branch("Conv_validvtx",&conv_validvtx,"Conv_validvtx[Conv_n]/O");
 fEventTree->Branch("Conv_ntracks",&conv_ntracks,"Conv_ntracks[Conv_n]/I");
 fEventTree->Branch("Conv_chi2_probability",&conv_chi2_probability,"Conv_chi2_probability[Conv_n]/F");
 fEventTree->Branch("Conv_eoverp",&conv_eoverp,"Conv_eoverp[Conv_n]/F");
 fEventTree->Branch("Conv_zofprimvtxfromtrks",&conv_zofprimvtxfromtrks,"Conv_zofprimvtxfromtrks[Conv_n]/F");

// fEventTree->Branch("Gvn",&gv_n,"Gvn/I");
// fEventTree->Branch("Gv_sumPtHi",&gv_sumPtHi,"Gv_sumPtHi[Gvn]/F");
// fEventTree->Branch("Gv_sumPtLo",&gv_sumPtLo,"Gv_sumPtLo[Gvn]/F");
// fEventTree->Branch("Gv_nTkHi",&gv_nTkHi,"Gv_nTkHi[Gvn]/I");
// fEventTree->Branch("Gv_nTkLo",&gv_nTkLo,"Gv_nTkLo[Gvn]/I");

// fEventTree->Branch("diphotons_first", &f_diphotons_first, Form("diphotons_first[%d]/I",gMax_vertexing_diphoton_pairs));
// fEventTree->Branch("diphotons_second", &f_diphotons_second, Form("diphotons_second[%d]/I",gMax_vertexing_diphoton_pairs));

// fEventTree->Branch("vtx_dipho_h2gglobe", &f_vtx_dipho_h2gglobe, Form("vtx_dipho_h2gglobe[%d][%d]/I",gMax_vertexing_diphoton_pairs,gMax_vertexing_vtxes));
// fEventTree->Branch("vtx_dipho_mva", &f_vtx_dipho_mva, Form("vtx_dipho_mva[%d][%d]/I",gMax_vertexing_diphoton_pairs,gMax_vertexing_vtxes));
// fEventTree->Branch("vtx_dipho_productrank", &f_vtx_dipho_productrank, Form("vtx_dipho_productrank[%d][%d]/I",gMax_vertexing_diphoton_pairs,gMax_vertexing_vtxes));

	fEventTree->Branch("NSuperClusters",&fTnSC ,"NSuperClusters/I");
	fEventTree->Branch("SCRaw",&fTSCraw ,"SCRaw[NSuperClusters]/F");
	fEventTree->Branch("SCPre",&fTSCpre ,"SCPre[NSuperClusters]/F");
	fEventTree->Branch("SCEnergy",&fTSCenergy ,"SCEnergy[NSuperClusters]/F");
	fEventTree->Branch("SCEta",&fTSCeta ,"SCEta[NSuperClusters]/F");
	fEventTree->Branch("SCPhi",&fTSCphi ,"SCPhi[NSuperClusters]/F");
	fEventTree->Branch("SCPhiWidth",&fTSCsigmaPhi ,"SCPhiWidth[NSuperClusters]/F");
	fEventTree->Branch("SCEtaWidth",&fTSCsigmaEta ,"SCEtaWidth[NSuperClusters]/F");
	fEventTree->Branch("SCBrem",&fTSCbrem ,"SCBrem[NSuperClusters]/F");
	fEventTree->Branch("SCR9",&fTSCR9 ,"SCR9[NSuperClusters]/F");
	fEventTree->Branch("SCcrackcorrseed",&fTSCcrackcorrseed ,"SCcrackcorrseed[NSuperClusters]/F");
	fEventTree->Branch("SCcrackcorr",&fTSCcrackcorr ,"SCcrackcorr[NSuperClusters]/F");
	fEventTree->Branch("SClocalcorrseed",&fTSClocalcorrseed ,"SClocalcorrseed[NSuperClusters]/F");
	fEventTree->Branch("SClocalcorr",&fTSClocalcorr ,"SClocalcorr[NSuperClusters]/F");
	fEventTree->Branch("SCcrackcorrseedfactor",&fTSCcrackcorrseedfactor ,"SCcrackcorrseedfactor[NSuperClusters]/F");
	fEventTree->Branch("SClocalcorrseedfactor",&fTSClocalcorrseedfactor ,"SClocalcorrseedfactor[NSuperClusters]/F");
	fEventTree->Branch("SCx",&fTSCx,"SCx[NSuperClusters]/F");
	fEventTree->Branch("SCy",&fTSCy,"SCy[NSuperClusters]/F");
	fEventTree->Branch("SCz",&fTSCz,"SCz[NSuperClusters]/F");
	fEventTree->Branch("SCNXtals",&fTSCNXtals,"SCNXtals[NSuperClusters]/I");
	fEventTree->Branch("SCxtalX",&fTSCxtalX,Form("SCxtalX[%d][%d]/F",gMaxnSC,gMaxnSCxtals));
	fEventTree->Branch("SCxtalY",&fTSCxtalY,Form("SCxtalY[%d][%d]/F",gMaxnSC,gMaxnSCxtals));
	fEventTree->Branch("SCxtalZ",&fTSCxtalZ,Form("SCxtalZ[%d][%d]/F",gMaxnSC,gMaxnSCxtals));
	fEventTree->Branch("SCxtalEtaWidth",&fTSCxtalEtaWidth,Form("SCxtalEtaWidth[%d][%d]/F",gMaxnSC,gMaxnSCxtals));
	fEventTree->Branch("SCxtalPhiWidth",&fTSCxtalPhiWidth,Form("SCxtalPhiWidth[%d][%d]/F",gMaxnSC,gMaxnSCxtals));
	fEventTree->Branch("SCxtalfrontX",&fTSCxtalfrontX,Form("SCxtalfrontX[%d][%d][4]/F",gMaxnSC,gMaxnSCxtals));
	fEventTree->Branch("SCxtalfrontY",&fTSCxtalfrontY,Form("SCxtalfrontY[%d][%d][4]/F",gMaxnSC,gMaxnSCxtals));
	fEventTree->Branch("SCxtalfrontZ",&fTSCxtalfrontZ,Form("SCxtalfrontZ[%d][%d][4]/F",gMaxnSC,gMaxnSCxtals));

	// Jets:
	fEventTree->Branch("NJets"          ,&fTnjets          ,"NJets/I");
	fEventTree->Branch("NJetsTot"       ,&fTnjetstot       ,"NJetsTot/I");
	fEventTree->Branch("JGood"          ,&fTgoodjet        ,"JGood[NJets]/I");
	fEventTree->Branch("JPx"            ,&fTjpx            ,"JPx[NJets]/F");
	fEventTree->Branch("JPy"            ,&fTjpy            ,"JPy[NJets]/F");
	fEventTree->Branch("JPz"            ,&fTjpz            ,"JPz[NJets]/F");
	fEventTree->Branch("JPt"            ,&fTjpt            ,"JPt[NJets]/F");
	fEventTree->Branch("JE"             ,&fTje             ,"JE[NJets]/F");
	fEventTree->Branch("JEt"            ,&fTjet            ,"JEt[NJets]/F");
	fEventTree->Branch("JEta"           ,&fTjeta           ,"JEta[NJets]/F");
	fEventTree->Branch("JPhi"           ,&fTjphi           ,"JPhi[NJets]/F");
	fEventTree->Branch("JEcorr"         ,&fTjEcorr         ,"JEcorr[NJets]/F");
	fEventTree->Branch("JArea"          ,&fTjArea          ,"JArea[NJets]/F");

	fEventTree->Branch("JEtaRms"             ,&fTJEtaRms        ,"JEtaRms[NJets]/F");
	fEventTree->Branch("JPhiRms"             ,&fTJPhiRms        ,"JPhiRms[NJets]/F");
	fEventTree->Branch("JNConstituents"      ,&fTjNconstituents ,"JNConstituents[NJets]/I");
	fEventTree->Branch("JNAssoTracks"        ,&fTjChMult        ,"JNAssoTracks[NJets]/I");
	fEventTree->Branch("JNNeutrals"          ,&fTjNeuMult       ,"JNNeutrals[NJets]/I");
	fEventTree->Branch("JChargedEmFrac"      ,&fTjChEmFrac      ,"JChargedEmFrac[NJets]/F");
	fEventTree->Branch("JNeutralEmFrac"      ,&fTjNeuEmFrac     ,"JNeutralEmFrac[NJets]/F");
	fEventTree->Branch("JChargedHadFrac"     ,&fTjChHadFrac     ,"JChargedHadFrac[NJets]/F");
	fEventTree->Branch("JNeutralHadFrac"     ,&fTjNeuHadFrac    ,"JNeutralHadFrac[NJets]/F");
	fEventTree->Branch("JChargedMuEnergyFrac",&fTjChMuEFrac     ,"JChargedMuEnergyFrac[NJets]/F");
	fEventTree->Branch("JPhotonEnergyFrac"   ,&fTjPhoFrac       ,"JPhotonEnergyFrac[NJets]/F");
	fEventTree->Branch("JHFHadEnergyFrac"    ,&fTjHFHadFrac     ,"JHFHadEnergyFrac[NJets]/F");
	fEventTree->Branch("JHFEMEnergyFrac"     ,&fTjHFEMFrac      ,"JHFEMEnergyFrac[NJets]/F");
	fEventTree->Branch("JPtD"                ,&fTjPtD           ,"JPtD[NJets]/F");
	fEventTree->Branch("JRMSCand"            ,&fTjRMSCand       ,"JRMSCand[NJets]/F");

	fEventTree->Branch("JeMinDR"               ,&fTjeMinDR               ,"JeMinDR[NJets]/F");
	fEventTree->Branch("JbTagProbTkCntHighEff" ,&fTjbTagProbTkCntHighEff ,"JbTagProbTkCntHighEff[NJets]/F");
	fEventTree->Branch("JbTagProbTkCntHighPur" ,&fTjbTagProbTkCntHighPur ,"JbTagProbTkCntHighPur[NJets]/F");
	fEventTree->Branch("JbTagProbSimpSVHighEff",&fTjbTagProbSimpSVHighEff,"JbTagProbSimpSVHighEff[NJets]/F");
	fEventTree->Branch("JbTagProbSimpSVHighPur",&fTjbTagProbSimpSVHighPur,"JbTagProbSimpSVHighPur[NJets]/F");

	fEventTree->Branch("JMass"          ,&fTjMass          ,"JMass[NJets]/F");
	fEventTree->Branch("Jtrk1px"        ,&fTjtrk1px        ,"Jtrk1px[NJets]/F");
	fEventTree->Branch("Jtrk1py"        ,&fTjtrk1py        ,"Jtrk1py[NJets]/F");
	fEventTree->Branch("Jtrk1pz"        ,&fTjtrk1pz        ,"Jtrk1pz[NJets]/F");
	fEventTree->Branch("Jtrk2px"        ,&fTjtrk2px        ,"Jtrk2px[NJets]/F");
	fEventTree->Branch("Jtrk2py"        ,&fTjtrk2py        ,"Jtrk2py[NJets]/F");
	fEventTree->Branch("Jtrk2pz"        ,&fTjtrk2pz        ,"Jtrk2pz[NJets]/F");
	fEventTree->Branch("Jtrk3px"        ,&fTjtrk3px        ,"Jtrk3px[NJets]/F");
	fEventTree->Branch("Jtrk3py"        ,&fTjtrk3py        ,"Jtrk3py[NJets]/F");
	fEventTree->Branch("Jtrk3pz"        ,&fTjtrk3pz        ,"Jtrk3pz[NJets]/F");
	fEventTree->Branch("JVtxx"          ,&fTjetVtxx        ,"JVtxx[NJets]/F");
	fEventTree->Branch("JVtxy"          ,&fTjetVtxy        ,"JVtxy[NJets]/F");
	fEventTree->Branch("JVtxz"          ,&fTjetVtxz        ,"JVtxz[NJets]/F");
	fEventTree->Branch("JVtxExx"        ,&fTjetVtxExx      ,"JVtxExx[NJets]/F");
	fEventTree->Branch("JVtxEyx"        ,&fTjetVtxEyx      ,"JVtxEyx[NJets]/F");
	fEventTree->Branch("JVtxEyy"        ,&fTjetVtxEyy      ,"JVtxEyy[NJets]/F");
	fEventTree->Branch("JVtxEzy"        ,&fTjetVtxEzy      ,"JVtxEzy[NJets]/F");
	fEventTree->Branch("JVtxEzz"        ,&fTjetVtxEzz      ,"JVtxEzz[NJets]/F");
	fEventTree->Branch("JVtxEzx"        ,&fTjetVtxEzx      ,"JVtxEzx[NJets]/F");
	fEventTree->Branch("JVtxNChi2"      ,&fTjetVtxNChi2    ,"JVtxNChi2[NJets]/F");
	fEventTree->Branch("JGenJetIndex"   ,&fTjetGenJetIndex ,"JGenJetIndex[NJets]/I");

	// Additional Jet collections
	for ( std::vector<JetFillerBase*>::iterator it = jetFillers.begin();
              it != jetFillers.end(); ++it )
          (*it)->createBranches();
	for ( std::vector<PatMuonFiller*>::iterator it = muonFillers.begin(); 
              it != muonFillers.end(); ++it ) 
          (*it)->createBranches();
	for ( std::vector<PatElectronFiller*>::iterator it = electronFillers.begin(); 
              it != electronFillers.end(); ++it ) 
          (*it)->createBranches();
	for ( std::vector<PatTauFiller*>::iterator it = tauFillers.begin(); 
              it != tauFillers.end(); ++it ) 
          (*it)->createBranches();


	// Tracks:
	fEventTree->Branch("NTracks"        ,&fTntracks      ,"NTracks/I");
	fEventTree->Branch("NTracksTot"     ,&fTntrackstot   ,"NTracksTot/I");
	fEventTree->Branch("TrkGood"        ,&fTgoodtrk      ,"TrkGood[NTracks]/I");
	fEventTree->Branch("TrkPt"          ,&fTtrkpt        ,"TrkPt[NTracks]/F");
	fEventTree->Branch("TrkEta"         ,&fTtrketa       ,"TrkEta[NTracks]/F");
	fEventTree->Branch("TrkPhi"         ,&fTtrkphi       ,"TrkPhi[NTracks]/F");
	fEventTree->Branch("TrkNChi2"       ,&fTtrknchi2     ,"TrkNChi2[NTracks]/F");
	fEventTree->Branch("TrkNHits"       ,&fTtrknhits     ,"TrkNHits[NTracks]/F");
	fEventTree->Branch("TrkVtxDz"       ,&fTtrkVtxDz     ,"TrkVtxDz[NTracks]/F");
	fEventTree->Branch("TrkVtxDxy"      ,&fTtrkVtxDxy    ,"TrkVtxDxy[NTracks]/F");
	fEventTree->Branch("TrkPtSumx"      ,&fTTrkPtSumx      ,"TrkPtSumx/F");
	fEventTree->Branch("TrkPtSumy"      ,&fTTrkPtSumy      ,"TrkPtSumy/F");
	fEventTree->Branch("TrkPtSum"       ,&fTTrkPtSum       ,"TrkPtSum/F");
	fEventTree->Branch("TrkPtSumPhi"    ,&fTTrkPtSumphi    ,"TrkPtSumPhi/F");

	// MET:
	fEventTree->Branch("SumEt"                 ,&fTSumEt               ,"SumEt/F");
	fEventTree->Branch("ECALSumEt"             ,&fTECALSumEt           ,"ECALSumEt/F");
	fEventTree->Branch("HCALSumEt"             ,&fTHCALSumEt           ,"HCALSumEt/F");
	fEventTree->Branch("ECALEsumx"             ,&fTECALEsumx           ,"ECALEsumx/F");
	fEventTree->Branch("ECALEsumy"             ,&fTECALEsumy           ,"ECALEsumy/F");
	fEventTree->Branch("ECALEsumz"             ,&fTECALEsumz           ,"ECALEsumz/F");
	fEventTree->Branch("ECALMET"               ,&fTECALMET             ,"ECALMET/F");
	fEventTree->Branch("ECALMETPhi"            ,&fTECALMETphi          ,"ECALMETPhi/F");
	fEventTree->Branch("ECALMETEta"            ,&fTECALMETeta          ,"ECALMETEta/F");
	fEventTree->Branch("HCALEsumx"             ,&fTHCALEsumx           ,"HCALEsumx/F");
	fEventTree->Branch("HCALEsumy"             ,&fTHCALEsumy           ,"HCALEsumy/F");
	fEventTree->Branch("HCALEsumz"             ,&fTHCALEsumz           ,"HCALEsumz/F");
	fEventTree->Branch("HCALMET"               ,&fTHCALMET             ,"HCALMET/F");
	fEventTree->Branch("HCALMETPhi"            ,&fTHCALMETphi          ,"HCALMETPhi/F");
	fEventTree->Branch("HCALMETeta"            ,&fTHCALMETeta          ,"HCALMETEta/F");
	fEventTree->Branch("RawMET"                ,&fTRawMET              ,"RawMET/F");
	fEventTree->Branch("RawMETpx"              ,&fTRawMETpx            ,"RawMETpx/F");
	fEventTree->Branch("RawMETpy"              ,&fTRawMETpy            ,"RawMETpy/F");
	fEventTree->Branch("RawMETphi"             ,&fTRawMETphi           ,"RawMETphi/F");
	fEventTree->Branch("RawMETemEtFrac"        ,&fTRawMETemEtFrac      ,"RawMETemEtFrac/F");
	fEventTree->Branch("RawMETemEtInEB"        ,&fTRawMETemEtInEB      ,"RawMETemEtInEB/F");
	fEventTree->Branch("RawMETemEtInEE"        ,&fTRawMETemEtInEE      ,"RawMETemEtInEE/F");
	fEventTree->Branch("RawMETemEtInHF"        ,&fTRawMETemEtInHF      ,"RawMETemEtInHF/F");
	fEventTree->Branch("RawMEThadEtFrac"       ,&fTRawMEThadEtFrac     ,"RawMEThadEtFrac/F");
	fEventTree->Branch("RawMEThadEtInHB"       ,&fTRawMEThadEtInHB     ,"RawMEThadEtInHB/F");
	fEventTree->Branch("RawMEThadEtInHE"       ,&fTRawMEThadEtInHE     ,"RawMEThadEtInHE/F");
	fEventTree->Branch("RawMEThadEtInHF"       ,&fTRawMEThadEtInHF     ,"RawMEThadEtInHF/F");
	fEventTree->Branch("RawMETSignificance"    ,&fTRawMETSignificance  ,"RawMETSignificance/F");
	fEventTree->Branch("GenMET"                ,&fTGenMET              ,"GenMET/F");
	fEventTree->Branch("GenMETpx"              ,&fTGenMETpx            ,"GenMETpx/F");
	fEventTree->Branch("GenMETpy"              ,&fTGenMETpy            ,"GenMETpy/F");
	fEventTree->Branch("GenMETphi"             ,&fTGenMETphi           ,"GenMETphi/F");
	fEventTree->Branch("TCMET"                 ,&fTTCMET               ,"TCMET/F");
	fEventTree->Branch("TCMETpx"               ,&fTTCMETpx             ,"TCMETpx/F");
	fEventTree->Branch("TCMETpy"               ,&fTTCMETpy             ,"TCMETpy/F");
	fEventTree->Branch("TCMETphi"              ,&fTTCMETphi            ,"TCMETphi/F");
	fEventTree->Branch("TCMETSignificance"     ,&fTTCMETSignificance   ,"TCMETSignificance/F");
	fEventTree->Branch("MuJESCorrMET"          ,&fTMuJESCorrMET        ,"MuJESCorrMET/F");
	fEventTree->Branch("MuJESCorrMETpx"        ,&fTMuJESCorrMETpx      ,"MuJESCorrMETpx/F");
	fEventTree->Branch("MuJESCorrMETpy"        ,&fTMuJESCorrMETpy      ,"MuJESCorrMETpy/F");
	fEventTree->Branch("MuJESCorrMETphi"       ,&fTMuJESCorrMETphi     ,"MuJESCorrMETphi/F");
	fEventTree->Branch("PFMET"                 ,&fTPFMET               ,"PFMET/F");
	fEventTree->Branch("PFMETpx"               ,&fTPFMETpx             ,"PFMETpx/F");
	fEventTree->Branch("PFMETpy"               ,&fTPFMETpy             ,"PFMETpy/F");
	fEventTree->Branch("PFMETphi"              ,&fTPFMETphi            ,"PFMETphi/F");
	fEventTree->Branch("PFMETSignificance"     ,&fTPFMETSignificance   ,"PFMETSignificance/F");
        fEventTree->Branch("PFSumEt"               ,&fTPFSumEt             ,"PFSumEt/F");
	fEventTree->Branch("PFMETPAT"              ,&fTPFMETPAT            ,"PFMETPAT/F");
	fEventTree->Branch("PFMETPATpx"            ,&fTPFMETPATpx          ,"PFMETPATpx/F");
	fEventTree->Branch("PFMETPATpy"            ,&fTPFMETPATpy          ,"PFMETPATpy/F");
	fEventTree->Branch("PFMETPATphi"           ,&fTPFMETPATphi         ,"PFMETPATphi/F");
	fEventTree->Branch("PFMETPATSignificance"  ,&fTPFMETPATSignificance   ,"PFMETPATSignificance/F");
	fEventTree->Branch("METR12"                ,&fTMETR12              ,"METR12/F");
	fEventTree->Branch("METR21"                ,&fTMETR21              ,"METR21/F");
}

// Method called once before each run
void NTupleProducer::beginRun(const edm::Run& r, const edm::EventSetup& es){
	if(!fIsRealData){
		edm::Handle<GenRunInfoProduct> genRunInfo;
		r.getByLabel("generator", genRunInfo);
		// Fill RunTree information
		fTextxslo       = genRunInfo->externalXSecLO().value();
		fTintxs         = genRunInfo->internalXSec().value();
	}



}

// Method called once after each run
void NTupleProducer::endRun(const edm::Run& r, const edm::EventSetup&){
	// Reset RunTree information
	fRTrunnumber     = -999;
	fRTextxslo       = -999.99;
	fRTextxsnlo      = -999.99;
	fRTintxs         = -999.99;
	fRTMinmupt       = -999.99;
	fRTMaxmueta      = -999.99;
	fRTMinelpt       = -999.99;
	fRTMaxeleta      = -999.99;
	fRTMinjpt        = -999.99;
	fRTMinrawjpt     = -999.99;
	fRTMaxjeta       = -999.99;
	fRTMinjemfrac    = -999.99;
	fRTMintrkpt      = -999.99;
	fRTMaxtrketa     = -999.99;
	fRTMaxtrknchi2   = -999.99;
	fRTMintrknhits   = -999;
	fRTMinphopt      = -999.99;
	fRTMaxphoeta     = -999.99;
	fRTMinSCraw      = -999.99;
	fRTMinSCrawPt    = -999.99;

	fRTmaxnmu   = -999;
	fRTmaxnel   = -999;
	fRTmaxnjet  = -999;
	fRTmaxntrk  = -999;
	fRTmaxnphot = -999;


	// Fill RunTree information
	fRTrunnumber     = r.id().run();
	edm::Handle<GenRunInfoProduct> genRunInfo;
	if(!fIsRealData){
		r.getByLabel("generator", genRunInfo);
		fRTextxslo       = genRunInfo->externalXSecLO().value();
		fRTextxsnlo      = genRunInfo->externalXSecNLO().value();
		fRTintxs         = genRunInfo->internalXSec().value();
	}
	fRTMinmupt       = fMinmupt;
	fRTMaxmueta      = fMaxmueta;
	fRTMinelpt       = fMinelpt;
	fRTMaxeleta      = fMaxeleta;
	fRTMinjpt        = fMincorjpt;
	fRTMinrawjpt     = fMinrawjpt;
	fRTMaxjeta       = fMaxjeta;
	fRTMinjemfrac    = fMinjemfrac;
	fRTMintrkpt      = fMintrkpt;
	fRTMaxtrketa     = fMaxtrketa;
	fRTMaxtrknchi2   = fMaxtrknchi2;
	fRTMintrknhits   = fMintrknhits;

	fRTMinphopt      = fMinphopt;
	fRTMaxphoeta     = fMaxphoeta;

	fRTMinSCraw      = fMinSCraw;
	fRTMinSCrawPt    = fMinSCrawPt;

	fRTmaxnmu   = gMaxnmus;
	fRTmaxnel   = gMaxneles;
	fRTmaxnjet  = gMaxnjets;
	fRTmaxntrk  = gMaxntrks;
	fRTmaxnphot = gMaxnphos;

	fRunTree->Fill();
}

// Method called once each job just after ending the event loop
void NTupleProducer::endJob(){

	// Build index in tree for fast retrieval by run
	fRunTree->BuildIndex("Run","0");
	edm::LogVerbatim("NTP") << " ---------------------------------------------------";
	edm::LogVerbatim("NTP") << " ==> NTupleProducer::endJob() ...";
	edm::LogVerbatim("NTP") << "  Total number of processed Events: " << fNTotEvents;
	edm::LogVerbatim("NTP") << "  Number of times Tree was filled:  " << fNFillTree;
	edm::LogVerbatim("NTP") << " ---------------------------------------------------";

	if (fNTotEvents!=fNFillTree) {
	  edm::LogVerbatim("NTP") << " ---------------------------------------------------";
	  edm::LogVerbatim("NTP") << " ==> WARNING!";
	  edm::LogVerbatim("NTP") << "  Total number of processed Events is not the same as Number of times Tree was filled ";
	  edm::LogVerbatim("NTP") << " ---------------------------------------------------";
	  
	}


}

// Method to reset the TTree variables for each event
void NTupleProducer::resetTree(){
	resetInt(fTHLTres, gMaxhltbits);
	resetInt(fTL1physres, gMaxl1physbits);
	resetInt(fTL1techres, gMaxl1techbits);
	resetInt(fTHLTprescale, gMaxhltbits);

	fTNHLTobjects = 0;
	for ( size_t i=0; i<fTNpaths; ++i ) {
		resetInt(fTHLTObjectID[i],gMaxhltnobjs);
		resetFloat(fTHLTObjectPt[i],gMaxhltnobjs);
		resetFloat(fTHLTObjectEta[i],gMaxhltnobjs);
		resetFloat(fTHLTObjectPhi[i],gMaxhltnobjs);
	}
	fTHLTmenu.clear();    fTHLTmenu.resize(gMaxhltbits);
	fTL1physmenu.clear(); fTL1physmenu.resize(gMaxl1physbits);

	fTrunnumber         = -999;
	fTeventnumber       = -999;
	fTlumisection       = -999;
	fTpthat             = -999.99;
	fTqcdPartonicHT     = -999.99;
	fTsigprocid         = -999;
	fTpdfscalePDF       = -999.99;
	fTpdfid1            = -999;
	fTpdfid2            = -999;
	fTpdfx1             = -999.99;
	fTpdfx2             = -999.99;
	fTpdfxPDF1          = -999.99;
	fTpdfxPDF2          = -999.99;
	fTgenweight         = -999.99;

	fTweight            = -999.99;
        fTpuWeightTotal     = -999.99;
        fTpuWeightInTime    = -999.99;

	fTgoodvtx           = -999;
	fTprimvtxx          = -999.99;
	fTprimvtxy          = -999.99;
	fTprimvtxz          = -999.99;
	fTprimvtxrho        = -999.99;
	fTprimvtxxE         = -999.99;
	fTprimvtxyE         = -999.99;
	fTprimvtxzE         = -999.99;
	fTpvtxznchi2        = -999.99;
	fTpvtxisfake        = -999;
	fTpvtxndof          = -999.99;
	fTpvtxptsum         = -999.99;
	fTbeamspotx         = -999.99;
	fTbeamspoty         = -999.99;
	fTbeamspotz         = -999.99;
	fTNCaloTowers       = -999;
	fTHBHENoiseFlag     = -999;
	fTHBHENoiseFlagIso  = -999;
	fTcscTightHaloID    = -999;
	fTrho               = -999;
	fTsigma               = -999;
	fTrhoPFnoPU         = -999;

	fTecalDeadTPFilterFlag = -999;
	fRecovRecHitFilterFlag = -999;
	fTra2TrackingFailureFilterFlag = -999;
	//FR fPBNRFlag              = -999;
	// fTEcalDeadCellBEFlag= -999;
	// fTnECALGapClusters  = 0;
	// resetFloat(fTEcalGapBE, gMaxnECALGapClusters);
	// resetInt   (fTEcalGapClusterSize, gMaxnECALGapClusters);

	resetFloat(fTvrtxx,     gMaxnvrtx);
	resetFloat(fTvrtxy,     gMaxnvrtx);
	resetFloat(fTvrtxz,     gMaxnvrtx);
	resetFloat(fTvrtxxE,    gMaxnvrtx);
	resetFloat(fTvrtxyE,    gMaxnvrtx);
	resetFloat(fTvrtxzE,    gMaxnvrtx);
	resetFloat(fTvrtxndof,  gMaxnvrtx);
	resetFloat(fTvrtxchi2,  gMaxnvrtx);
	resetFloat(fTvrtxntrks, gMaxnvrtx);
	resetFloat(fTvrtxsumpt, gMaxnvrtx);
	resetInt(fTvrtxisfake,   gMaxnvrtx);

	fTnGenParticles = 0;
	resetInt(fTgenInfoStatus,nStoredGenParticles);
	resetInt(fTgenInfoId,nStoredGenParticles);
	resetFloat(fTgenInfoMass,nStoredGenParticles);
	resetInt(fTgenInfoNMo,nStoredGenParticles);
	resetFloat(fTgenInfoMo1Pt,nStoredGenParticles);
	resetFloat(fTgenInfoMo2Pt,nStoredGenParticles);
	resetInt(fTgenInfoNDa,nStoredGenParticles);
	resetInt(fTgenInfoMo1,nStoredGenParticles);
	resetInt(fTgenInfoMo2,nStoredGenParticles);
	resetInt(fTgenInfoDa1,nStoredGenParticles);
	resetInt(fTgenInfoDa2,nStoredGenParticles);
	resetFloat(fTgenInfoPt,nStoredGenParticles);
	resetFloat(fTgenInfoEta,nStoredGenParticles);
	resetFloat(fTgenInfoPhi,nStoredGenParticles);
	resetFloat(fTgenInfoPx,nStoredGenParticles);
	resetFloat(fTgenInfoPy,nStoredGenParticles);
	resetFloat(fTgenInfoPz,nStoredGenParticles);
	resetFloat(fTgenInfoM,nStoredGenParticles);
	resetInt(fTgenInfoMoIndex,nStoredGenParticles);
	resetInt(fTPromptnessLevel,nStoredGenParticles);

	// Pile-up
	fTpuNumInteractions    = -999;
	fTpuNumFilled = -999;
	fTpuOOTNumInteractionsEarly = -999;
	fTpuOOTNumInteractionsLate = -999;

	resetFloat(fTpuZpositions   ,gMaxnpileup);
	resetFloat(fTpuSumpT_lowpT  ,gMaxnpileup);
	resetFloat(fTpuSumpT_highpT ,gMaxnpileup);
	resetFloat(fTpuNtrks_lowpT  ,gMaxnpileup);
	resetFloat(fTpuNtrks_highpT ,gMaxnpileup);


	// resetFloat(fTpuInstLumi     ,gMaxnpileup);

	fTgoodevent         = 0;
	fTflagmaxmuexc      = 0;
	fTflagmaxelexc      = 0;
	fTflagmaxujetexc    = 0;
	fTflagmaxjetexc     = 0;
	fTflagmaxtrkexc     = 0;
	fTflagmaxphoexc     = 0;
	fTflagmaxgenleptexc = 0;
	fTflagmaxgenphotexc = 0;
	fTflagmaxgenjetexc = 0;
	fTflagmaxvrtxexc    = 0;
	fTflagmaxgenpartexc = 0;


	fTnmu         = 0;
	fTnmutot      = 0;
	fTnglobalmu   = 0;
	fTntrackermu  = 0;
	fTneles       = 0;
	fTnelestot    = 0;
	fTnjets       = 0;
	fTnjetstot    = 0;
	fTntracks     = 0;
	fTntrackstot  = 0;
	fTNPfCand     = 0;
	fTnphotons    = 0;
	fTnphotonstot = 0;
	fTngenleptons = 0;
	fTngenphotons = 0;
	fTNGenJets    = 0;
	fTnvrtx       = 0;
	fTnSC         = 0;
	conv_n        = 0;
	gv_n          = 0;

	fTnGenParticles = 0;
	resetInt(fTGenLeptonId       ,gMaxngenlept);
	resetFloat(fTGenLeptonPt    ,gMaxngenlept);
	resetFloat(fTGenLeptonEta   ,gMaxngenlept);
	resetFloat(fTGenLeptonPhi   ,gMaxngenlept);
	resetInt(fTGenLeptonMId      ,gMaxngenlept);
	resetInt(fTGenLeptonMStatus  ,gMaxngenlept);
	resetFloat(fTGenLeptonMPt   ,gMaxngenlept);
	resetFloat(fTGenLeptonMEta  ,gMaxngenlept);
	resetFloat(fTGenLeptonMPhi  ,gMaxngenlept);
	resetInt(fTGenLeptonGMId     ,gMaxngenlept);
	resetInt(fTGenLeptonGMStatus ,gMaxngenlept);
	resetFloat(fTGenLeptonGMPt  ,gMaxngenlept);
	resetFloat(fTGenLeptonGMEta ,gMaxngenlept);
	resetFloat(fTGenLeptonGMPhi ,gMaxngenlept);

       resetFloat(fTGenPhotonPt    ,gMaxngenphot);
       resetFloat(fTGenPhotonEta   ,gMaxngenphot);
       resetFloat(fTGenPhotonPhi   ,gMaxngenphot);
       resetFloat(fTGenPhotonPartonMindR   ,gMaxngenphot);
       resetInt(fTGenPhotonMotherID   ,gMaxngenphot);
       resetInt(fTGenPhotonMotherStatus ,gMaxngenphot);
       resetFloat(fTGenPhotonIsoDR03 ,gMaxngenphot);
       resetFloat(fTGenPhotonIsoDR04 ,gMaxngenphot);

	resetFloat(fTGenJetPt   ,gMaxngenjets);
	resetFloat(fTGenJetEta  ,gMaxngenjets);
	resetFloat(fTGenJetPhi  ,gMaxngenjets);
	resetFloat(fTGenJetE    ,gMaxngenjets);
	resetFloat(fTGenJetemE  ,gMaxngenjets);
	resetFloat(fTGenJethadE ,gMaxngenjets);
	resetFloat(fTGenJetinvE ,gMaxngenjets);

	resetInt(fTgoodmu, gMaxnmus);
	resetInt(fTmuIsIso, gMaxnmus);
	resetInt(fTmuIsGM, gMaxnmus);
	resetInt(fTmuIsTM, gMaxnmus);
	resetFloat(fTmupx, gMaxnmus);
	resetFloat(fTmupy, gMaxnmus);
	resetFloat(fTmupz, gMaxnmus);
	resetFloat(fTmupt, gMaxnmus);
	resetFloat(fTmuinnerpt, gMaxnmus);
	resetFloat(fTmuptE, gMaxnmus);
	resetFloat(fTmue, gMaxnmus);
	resetFloat(fTmuet, gMaxnmus);
	resetFloat(fTmueta, gMaxnmus);
	resetFloat(fTmuphi, gMaxnmus);
	resetInt(fTmucharge, gMaxnmus);
	resetFloat(fTmuiso, gMaxnmus);
	resetFloat(fTmuIso03sumPt, gMaxnmus);
	resetFloat(fTmuIso03emEt, gMaxnmus);
	resetFloat(fTmuIso03hadEt, gMaxnmus);
	resetFloat(fTmuIso03emVetoEt, gMaxnmus);
	resetFloat(fTmuIso03hadVetoEt, gMaxnmus);
	resetFloat(fTmuIso05sumPt, gMaxnmus);
	resetFloat(fTmuIso05emEt, gMaxnmus);
	resetFloat(fTmuIso05hadEt, gMaxnmus);
	resetFloat(fTmueecal, gMaxnmus);
	resetFloat(fTmuehcal, gMaxnmus);
	resetFloat(fTmud0bs, gMaxnmus);
	resetFloat(fTmud0pv, gMaxnmus);
	resetFloat(fTmud0E, gMaxnmus);
	resetFloat(fTmudzbs, gMaxnmus);
	resetFloat(fTmudzpv, gMaxnmus);
	resetFloat(fTmudzE, gMaxnmus);
	resetFloat(fTmunchi2, gMaxnmus);
	resetInt(fTmunglhits, gMaxnmus);
	resetInt(fTmunmuhits, gMaxnmus);
	resetInt(fTmuntkhits, gMaxnmus);
	resetInt(fTmunpxhits, gMaxnmus);
	resetFloat(fTmuinntknchi2, gMaxnmus);
	resetInt(fTmunmatches, gMaxnmus);
	resetInt(fTmunchambers, gMaxnmus);
	resetFloat(fTmucalocomp, gMaxnmus);
	resetFloat(fTmusegmcomp, gMaxnmus);

	resetInt(fTmuIsGMPT, gMaxnmus);
	resetInt(fTmuIsGMTkChiComp, gMaxnmus);
	resetInt(fTmuIsGMStaChiComp, gMaxnmus);
	resetInt(fTmuIsGMTkKinkTight, gMaxnmus);
	resetInt(fTmuIsAllStaMuons, gMaxnmus);
	resetInt(fTmuIsAllTrkMuons, gMaxnmus);
	resetInt(fTmuIsTrkMuArb, gMaxnmus);
	resetInt(fTmuIsAllArb, gMaxnmus);
	resetInt(fTmuIsTMLastStationLoose, gMaxnmus);
	resetInt(fTmuIsTMLastStationTight, gMaxnmus);
	resetInt(fTmuIsTM2DCompLoose, gMaxnmus);
	resetInt(fTmuIsTM2DCompTight, gMaxnmus);
	resetInt(fTmuIsTMOneStationLoose, gMaxnmus);
	resetInt(fTmuIsTMOneStationTight, gMaxnmus);
	resetInt(fTmuIsTMLSOPL, gMaxnmus);
	resetInt(fTmuIsTMLastStationAngLoose, gMaxnmus);
	resetInt(fTmuIsTMLastStationAngTight, gMaxnmus);
	resetInt(fTmuIsTMOneStationAngTight, gMaxnmus);
	resetInt(fTmuIsTMOneStationAngLoose, gMaxnmus);

	resetInt(fTGenMuId, gMaxnmus);
	resetInt(fTGenMuStatus, gMaxnmus);
	resetFloat(fTGenMuPt, gMaxnmus);
	resetFloat(fTGenMuEta, gMaxnmus);
	resetFloat(fTGenMuPhi, gMaxnmus);
	resetFloat(fTGenMuE, gMaxnmus);
	resetInt(fTGenMuMId, gMaxnmus);
	resetInt(fTGenMuMStatus, gMaxnmus);
	resetFloat(fTGenMuMPt, gMaxnmus);
	resetFloat(fTGenMuMEta, gMaxnmus);
	resetFloat(fTGenMuMPhi, gMaxnmus);
	resetFloat(fTGenMuME, gMaxnmus);
	resetInt(fTGenMuGMId, gMaxnmus);
	resetInt(fTGenMuGMStatus, gMaxnmus);
	resetFloat(fTGenMuGMPt, gMaxnmus);
	resetFloat(fTGenMuGMEta, gMaxnmus);
	resetFloat(fTGenMuGMPhi, gMaxnmus);
	resetFloat(fTGenMuGME, gMaxnmus);

	resetInt(fTgoodel, gMaxneles);
	resetInt(fTeIsIso, gMaxneles);
	resetInt(fTeChargeMisIDProb, gMaxneles);
	resetFloat(fTepx, gMaxneles);
	resetFloat(fTepy, gMaxneles);
	resetFloat(fTepz, gMaxneles);
	resetFloat(fTee, gMaxneles);
	resetFloat(fTeet, gMaxneles);
	resetFloat(fTept, gMaxneles);
	resetFloat(fTeptE, gMaxneles);
	resetFloat(fTeeta, gMaxneles);
	resetFloat(fTephi, gMaxneles);
	resetFloat(fTegsfpt , gMaxneles);
	resetFloat(fTegsfeta, gMaxneles);
	resetFloat(fTegsfphi, gMaxneles);
	resetFloat(fTetrkmomerror, gMaxneles);
	resetFloat(fTeecalergerror, gMaxneles);
	resetFloat(fTeelemomerror, gMaxneles);
	resetInt(  fTenbrems, gMaxneles);
	resetFloat(fTed0bs, gMaxneles);
	resetFloat(fTed0pv, gMaxneles);
	resetFloat(fTed0E, gMaxneles);
	resetFloat(fTedzbs, gMaxneles);
	resetFloat(fTedzpv, gMaxneles);
	resetFloat(fTedzE, gMaxneles);
	resetFloat(fTenchi2, gMaxneles);
	resetFloat(fTeiso03, gMaxneles);
	resetFloat(fTeiso04, gMaxneles);
	resetFloat(fTdr03tksumpt, gMaxneles);
	resetFloat(fTdr04tksumpt, gMaxneles);
	resetFloat(fTdr03ecalrechitsumet, gMaxneles);
	resetFloat(fTdr04ecalrechitsumet, gMaxneles);
	resetFloat(fTdr03hcaltowersumet, gMaxneles);
	resetFloat(fTdr04hcaltowersumet, gMaxneles);
	resetFloat(fTetheta, gMaxneles);
	resetFloat(fTesceta, gMaxneles);
	resetInt(fTecharge, gMaxneles);
	resetInt(fTeCInfoIsGsfCtfCons, gMaxneles);
	resetInt(fTeCInfoIsGsfCtfScPixCons, gMaxneles);
	resetInt(fTeCInfoIsGsfScPixCons, gMaxneles);
	resetInt(fTeCInfoScPixCharge, gMaxneles);
	resetFloat(fTeClosestCtfTrackpt, gMaxneles);
	resetFloat(fTeClosestCtfTracketa, gMaxneles);
	resetFloat(fTeClosestCtfTrackphi, gMaxneles);
	resetInt(fTeClosestCtfTrackcharge, gMaxneles);
	resetInt(fTeInGap, gMaxneles);
	resetInt(fTeEcalDriven, gMaxneles);
	resetInt(fTeTrackerDriven, gMaxneles);
	resetInt(fTeBasicClustersSize, gMaxneles);
	resetFloat(fTefbrem, gMaxneles);
	resetFloat(fTeHcalOverEcal, gMaxneles);
	resetFloat(fTeE1x5, gMaxneles);
	resetFloat(fTeE5x5, gMaxneles);
	resetFloat(fTeE2x5Max, gMaxneles);
	resetFloat(fTeSigmaIetaIeta, gMaxneles);
	resetFloat(fTeDeltaPhiSeedClusterAtCalo, gMaxneles);
	resetFloat(fTeDeltaEtaSeedClusterAtCalo, gMaxneles);
	resetFloat(fTeDeltaPhiSuperClusterAtVtx, gMaxneles);
	resetFloat(fTeDeltaEtaSuperClusterAtVtx, gMaxneles);
	resetFloat(fTecaloenergy, gMaxneles);
	resetFloat(fTetrkmomatvtx, gMaxneles);
	resetFloat(fTeESuperClusterOverP, gMaxneles);
	resetInt(fTeNumberOfMissingInnerHits, gMaxneles);
	// resetInt(fTeIsInJet, gMaxneles);
	// resetFloat(fTeSharedPx, gMaxneles);
	// resetFloat(fTeSharedPy, gMaxneles);
	// resetFloat(fTeSharedPz, gMaxneles);
	// resetFloat(fTeSharedEnergy, gMaxneles);
	// resetInt(fTeDupEl, gMaxneles);
	resetInt (fTElSCindex, gMaxnphos);

	resetFloat(fTeConvPartTrackDist, gMaxneles);
	resetFloat(fTeConvPartTrackDCot, gMaxneles);
	resetFloat(fTeConvPartTrackPt, gMaxneles);
	resetFloat(fTeConvPartTrackEta, gMaxneles);
	resetFloat(fTeConvPartTrackPhi, gMaxneles);
	resetFloat(fTeConvPartTrackCharge, gMaxneles);

	resetInt(fTeScSeedSeverity, gMaxneles);
	resetFloat(fTeS4OverS1, gMaxneles);
	resetFloat(fTeE1OverE9, gMaxneles);

	resetFloat(fTeIDMva, gMaxneles);
	resetInt(fTeIDTight, gMaxneles);
	resetInt(fTeIDLoose, gMaxneles);
	resetInt(fTeIDRobustTight, gMaxneles);
	resetInt(fTeIDRobustLoose, gMaxneles);
	resetInt(fTeIDsimpleWPrelIso, gMaxneles);
	resetInt(fTeIDsimpleWP95relIso, gMaxneles);
	resetInt(fTeIDsimpleWP90relIso, gMaxneles);
	resetInt(fTeIDsimpleWP85relIso, gMaxneles);
	resetInt(fTeIDsimpleWP80relIso, gMaxneles);
	resetInt(fTGenElId, gMaxneles);
	resetInt(fTGenElStatus, gMaxneles);
	resetFloat(fTGenElPt, gMaxneles);
	resetFloat(fTGenElEta, gMaxneles);
	resetFloat(fTGenElPhi, gMaxneles);
	resetFloat(fTGenElE, gMaxneles);
	resetInt(fTGenElMId, gMaxneles);
	resetInt(fTGenElMStatus, gMaxneles);
	resetFloat(fTGenElMPt, gMaxneles);
	resetFloat(fTGenElMEta, gMaxneles);
	resetFloat(fTGenElMPhi, gMaxneles);
	resetFloat(fTGenElME, gMaxneles);
	resetInt(fTGenElGMId, gMaxneles);
	resetInt(fTGenElGMStatus, gMaxneles);
	resetFloat(fTGenElGMPt, gMaxneles);
	resetFloat(fTGenElGMEta, gMaxneles);
	resetFloat(fTGenElGMPhi, gMaxneles);
	resetFloat(fTGenElGME, gMaxneles);

	resetInt(fTgoodjet, gMaxnjets);
	resetFloat(fTjpx,  gMaxnjets);
	resetFloat(fTjpy,  gMaxnjets);
	resetFloat(fTjpz,  gMaxnjets);
	resetFloat(fTje,   gMaxnjets);
	resetFloat(fTjet,  gMaxnjets);
	resetFloat(fTjpt,  gMaxnjets);
	resetFloat(fTjeta, gMaxnjets);
	resetFloat(fTjphi, gMaxnjets);
	resetFloat(fTJEtaRms,gMaxnjets );
	resetFloat(fTJPhiRms,gMaxnjets );
	resetFloat(fTjMass,   gMaxnjets);
	resetFloat(fTjArea,   gMaxnjets);
	resetFloat(fTjEcorr,   gMaxnjets);

	resetInt  (fTjNconstituents, gMaxnjets);
	resetInt  (fTjChMult, gMaxnjets);
	resetInt  (fTjNeuMult, gMaxnjets);
	resetFloat(fTjChHadFrac, gMaxnjets);
	resetFloat(fTjNeuHadFrac, gMaxnjets);
	resetFloat(fTjChEmFrac, gMaxnjets);
	resetFloat(fTjNeuEmFrac, gMaxnjets);
	resetFloat(fTjChMuEFrac, gMaxnjets);
	resetFloat(fTjPhoFrac, gMaxnjets);
	resetFloat(fTjHFHadFrac, gMaxnjets);
	resetFloat(fTjHFEMFrac, gMaxnjets);
	resetFloat(fTjPtD, gMaxnjets);
	resetFloat(fTjRMSCand, gMaxnjets);

	resetFloat(fTjbTagProbTkCntHighEff, gMaxnjets);
	resetFloat(fTjbTagProbTkCntHighPur, gMaxnjets);
	resetFloat(fTjbTagProbSimpSVHighEff, gMaxnjets);
	resetFloat(fTjbTagProbSimpSVHighPur, gMaxnjets);
	resetFloat(fTjtrk1px, gMaxnjets);
	resetFloat(fTjtrk1py, gMaxnjets);
	resetFloat(fTjtrk1pz, gMaxnjets);
	resetFloat(fTjtrk2px, gMaxnjets);
	resetFloat(fTjtrk2py, gMaxnjets);
	resetFloat(fTjtrk2pz, gMaxnjets);
	resetFloat(fTjtrk3px, gMaxnjets);
	resetFloat(fTjtrk3py, gMaxnjets);
	resetFloat(fTjtrk3pz, gMaxnjets);
	resetFloat(fTjeMinDR, gMaxnjets);
	resetFloat(fTjetVtxx, gMaxnjets);
	resetFloat(fTjetVtxy, gMaxnjets);
	resetFloat(fTjetVtxz, gMaxnjets);
	resetFloat(fTjetVtxExx, gMaxnjets);
	resetFloat(fTjetVtxEyx, gMaxnjets);
	resetFloat(fTjetVtxEyy, gMaxnjets);
	resetFloat(fTjetVtxEzy, gMaxnjets);
	resetFloat(fTjetVtxEzz, gMaxnjets);
	resetFloat(fTjetVtxEzx, gMaxnjets);
	resetFloat(fTjetVtxNChi2, gMaxnjets);

	resetInt(fTjetGenJetIndex, gMaxnjets);

	resetInt(fTgoodtrk,  gMaxntrks);
	resetFloat(fTtrkpt, gMaxntrks);
	resetFloat(fTtrketa, gMaxntrks);
	resetFloat(fTtrkphi, gMaxntrks);
	resetFloat(fTtrknchi2, gMaxntrks);
	resetFloat(fTtrknhits, gMaxntrks);
	resetFloat(fTtrkVtxDxy, gMaxntrks);
	resetFloat(fTtrkVtxDz, gMaxntrks);


	
	resetFloat(fTPfCandPdgId     ,gMaxnpfcand);
	resetFloat(fTPfCandEta       ,gMaxnpfcand);
	resetFloat(fTPfCandPhi       ,gMaxnpfcand);
	resetFloat(fTPfCandPx        ,gMaxnpfcand);
	resetFloat(fTPfCandPy        ,gMaxnpfcand);
	resetFloat(fTPfCandPz        ,gMaxnpfcand);
	resetFloat(fTPfCandEnergy    ,gMaxnpfcand);
	resetFloat(fTPfCandPt        ,gMaxnpfcand);
	resetFloat(fTPfCandVx        ,gMaxnpfcand);
	resetFloat(fTPfCandVy        ,gMaxnpfcand);
	resetFloat(fTPfCandVz        ,gMaxnpfcand);
	resetInt(fTPfCandBelongsToJet,gMaxnpfcand);
	resetFloat(fTPfCandMomX      ,gMaxnpfcand);
	resetFloat(fTPfCandMomY      ,gMaxnpfcand);
	resetFloat(fTPfCandMomZ      ,gMaxnpfcand);
	resetInt(fTPfCandHasHitInFirstPixelLayer, gMaxnpfcand);
	resetFloat(fTPfCandTrackRefPx, gMaxnpfcand);
	resetFloat(fTPfCandTrackRefPy, gMaxnpfcand);
	resetFloat(fTPfCandTrackRefPz, gMaxnpfcand);
	resetFloat(fTPfCandTrackRefVx, gMaxnpfcand);
	resetFloat(fTPfCandTrackRefVy, gMaxnpfcand);
	resetFloat(fTPfCandTrackRefVz, gMaxnpfcand);


	resetFloat(fTPhotVx,gMaxnphos);
	resetFloat(fTPhotVy,gMaxnphos);
	resetFloat(fTPhotVz,gMaxnphos);
	resetFloat(fTPhotPt,gMaxnphos);
	resetFloat(fTPhotPx,gMaxnphos);
	resetFloat(fTPhotPy,gMaxnphos);
	resetFloat(fTPhotPz,gMaxnphos);
	resetFloat(fTPhotEta,gMaxnphos);
	resetFloat(fTPhotPhi,gMaxnphos);
	resetFloat(fTPhotEnergy,gMaxnphos);
	resetFloat(fTPhotIso03Ecal,gMaxnphos);
	resetFloat(fTPhotIso03Hcal,gMaxnphos);
	resetFloat(fTPhotIso03TrkSolid,gMaxnphos);
	resetFloat(fTPhotIso03TrkHollow,gMaxnphos);
	resetFloat(fTPhotIso03,gMaxnphos);
	resetFloat(fTPhotIso04Ecal,gMaxnphos);
	resetFloat(fTPhotIso04Hcal,gMaxnphos);
	resetFloat(fTPhotIso04TrkSolid,gMaxnphos);
	resetFloat(fTPhotIso04TrkHollow,gMaxnphos);
	resetFloat(fTPhotIso04,gMaxnphos);
	resetFloat(fTPhotR9,gMaxnphos);
	resetFloat(fTPhotcaloPosX,gMaxnphos);
	resetFloat(fTPhotcaloPosY,gMaxnphos);
	resetFloat(fTPhotcaloPosZ,gMaxnphos);
	resetFloat(fTPhotHoverE,gMaxnphos);
	resetFloat(fTPhotH1overE,gMaxnphos);
	resetFloat(fTPhotH2overE,gMaxnphos);
	resetFloat(fTPhotSigmaIetaIeta,gMaxnphos);
	resetFloat(fTPhotSCEnergy,gMaxnphos);
	resetFloat(fTPhotSCEtaWidth,gMaxnphos);
	resetFloat(fTPhotSCSigmaPhiPhi,gMaxnphos);
	resetInt(fTPhotHasPixSeed,gMaxnphos);
	resetInt(fTPhotPassConvSafeElectronVeto,gMaxnphos);
	resetInt(fTPhotHasConvTrks,gMaxnphos);
	resetInt(fTgoodphoton,gMaxnphos);
	resetInt(fTPhotIsIso,gMaxnphos);
	// resetInt(fTPhotIsInJet,gMaxnphos);
	// resetInt(fTPhotDupEl,gMaxnphos);
	// resetFloat(fTPhotSharedPx, gMaxnphos);
	// resetFloat(fTPhotSharedPy, gMaxnphos);
	// resetFloat(fTPhotSharedPz, gMaxnphos);
	// resetFloat(fTPhotSharedEnergy, gMaxnphos);

	resetInt(fTPhotScSeedSeverity, gMaxnphos);
	resetFloat(fTPhotS4OverS1, gMaxnphos);
	resetFloat(fTPhotE1OverE9, gMaxnphos);

       resetFloat(fTPhotSigmaEtaEta,gMaxnphos);
       resetFloat(fTPhote1x5,gMaxnphos);
       resetFloat(fTPhote2x5,gMaxnphos);
       resetFloat(fTPhote3x3,gMaxnphos);
       resetFloat(fTPhote5x5,gMaxnphos);
       resetFloat(fTPhotmaxEnergyXtal,gMaxnphos);
       resetFloat(fTPhotIso03HcalDepth1,gMaxnphos);
       resetFloat(fTPhotIso03HcalDepth2,gMaxnphos);
       resetFloat(fTPhotIso04HcalDepth1,gMaxnphos);
       resetFloat(fTPhotIso04HcalDepth2,gMaxnphos);
       resetInt(fTPhotIso03nTrksSolid,gMaxnphos);
       resetInt(fTPhotIso03nTrksHollow,gMaxnphos);
       resetInt(fTPhotIso04nTrksSolid,gMaxnphos);
       resetInt(fTPhotIso04nTrksHollow,gMaxnphos);
       resetInt(fTPhotisEB,gMaxnphos);
       resetInt(fTPhotisEE,gMaxnphos);
       resetInt(fTPhotisEBEtaGap,gMaxnphos);
       resetInt(fTPhotisEBPhiGap,gMaxnphos);
       resetInt(fTPhotisEERingGap,gMaxnphos);
       resetInt(fTPhotisEEDeeGap,gMaxnphos);
       resetInt(fTPhotisEBEEGap,gMaxnphos);
       resetInt(fTPhotisPFlowPhoton,gMaxnphos);
       resetInt(fTPhotisStandardPhoton,gMaxnphos);
       resetInt(fTPhotMCmatchindex,gMaxnphos);
       resetInt(fTPhotMCmatchexitcode,gMaxnphos);
       resetFloat( fT_pho_ChargedHadronIso, gMaxnphos);
       resetFloat( fT_pho_NeutralHadronIso, gMaxnphos);
       resetFloat( fT_pho_PhotonIso, gMaxnphos);
       resetFloat( fTPhoSCRemovalPFIsoCharged, gMaxnphos);
       resetFloat( fTPhoSCRemovalPFIsoChargedPrimVtx, gMaxnphos);
       resetFloat( fTPhoSCRemovalPFIsoNeutral, gMaxnphos);
       resetFloat( fTPhoSCRemovalPFIsoPhoton, gMaxnphos);
       resetFloat( fTPhoSCRemovalPFIsoCharged_RCone, gMaxnphos);
       resetFloat( fTPhoSCRemovalPFIsoChargedPrimVtx_RCone, gMaxnphos);
       resetFloat( fTPhoSCRemovalPFIsoNeutral_RCone, gMaxnphos);
       resetFloat( fTPhoSCRemovalPFIsoPhoton_RCone, gMaxnphos);
       resetFloat( fTPhoSCRemoval_RCone_Eta, gMaxnphos);
       resetFloat( fTPhoSCRemoval_RCone_Phi, gMaxnphos);
       resetInt( fT_pho_isPFPhoton, gMaxnphos);
       resetInt( fT_pho_isPFElectron, gMaxnphos);
       resetInt (fTPhotSCindex, gMaxnphos);
       resetInt(fT_pho_matchedPFPhotonCand, gMaxnphos);
       resetInt(fT_pho_matchedPFElectronCand, gMaxnphos);

       for (int i=0; i<gMaxnphos; i++) pho_conv_validvtx[i]=false;
       resetFloat(pho_conv_chi2_probability,gMaxnphos);
       resetFloat(pho_conv_eoverp,gMaxnphos);
       resetInt(pho_conv_ntracks,gMaxnphos);

       for (int i=0; i<gMaxnconv; i++) conv_validvtx[i]=false;
       resetInt(conv_ntracks,gMaxnconv);
       resetFloat(conv_chi2_probability,gMaxnconv);
       resetFloat(conv_eoverp,gMaxnconv);
       resetFloat(conv_zofprimvtxfromtrks,gMaxnconv);

       resetFloat(gv_sumPtHi,gMaxngenvtx);
       resetFloat(gv_sumPtLo,gMaxngenvtx);
       resetInt(gv_nTkHi,gMaxngenvtx);
       resetInt(gv_nTkLo,gMaxngenvtx);

       resetInt(f_diphotons_first,gMax_vertexing_diphoton_pairs);
       resetInt(f_diphotons_second,gMax_vertexing_diphoton_pairs);

       for (int i=0; i<gMax_vertexing_diphoton_pairs; i++) {
	 for (int j=0; j<gMax_vertexing_vtxes; j++) {
           f_vtx_dipho_h2gglobe[i][j]=-999;
           f_vtx_dipho_mva[i][j]=-999;
           f_vtx_dipho_productrank[i][j]=-999;
         }
       }


       resetFloat(fT_pho_Cone01PhotonIso_dEta015EB_dR070EE_mvVtx, gMaxnphos);
       resetFloat(fT_pho_Cone02PhotonIso_dEta015EB_dR070EE_mvVtx, gMaxnphos);
       resetFloat(fT_pho_Cone03PhotonIso_dEta015EB_dR070EE_mvVtx, gMaxnphos);
       resetFloat(fT_pho_Cone04PhotonIso_dEta015EB_dR070EE_mvVtx, gMaxnphos);
       
       resetFloat(fT_pho_Cone01NeutralHadronIso_mvVtx, gMaxnphos);
       resetFloat(fT_pho_Cone02NeutralHadronIso_mvVtx, gMaxnphos);
       resetFloat(fT_pho_Cone03NeutralHadronIso_mvVtx, gMaxnphos);
       resetFloat(fT_pho_Cone04NeutralHadronIso_mvVtx, gMaxnphos);
       
       resetFloat(fT_pho_Cone01ChargedHadronIso_dR02_dz02_dxy01, gMaxnphos);
       resetFloat(fT_pho_Cone02ChargedHadronIso_dR02_dz02_dxy01, gMaxnphos);
       resetFloat(fT_pho_Cone03ChargedHadronIso_dR02_dz02_dxy01, gMaxnphos);
       resetFloat(fT_pho_Cone04ChargedHadronIso_dR02_dz02_dxy01, gMaxnphos);

       resetFloat(fT_pho_Cone03PFCombinedIso, gMaxnphos);
       resetFloat(fT_pho_Cone04PFCombinedIso, gMaxnphos);

       /*
       resetFloat(fT_pho_Cone04PhotonIso_dR0_dEta0_pt0,gMaxnphos);
       resetFloat(fT_pho_Cone04PhotonIso_dR0_dEta0_pt5,gMaxnphos);
       resetFloat(fT_pho_Cone04PhotonIso_dR8_dEta0_pt0,gMaxnphos);
       resetFloat(fT_pho_Cone04PhotonIso_dR8_dEta0_pt5,gMaxnphos);
       resetFloat(fT_pho_Cone01PhotonIso_dR045EB070EE_dEta015_pt08EB1EE_mvVtx,gMaxnphos);
       resetFloat(fT_pho_Cone02PhotonIso_dR045EB070EE_dEta015_pt08EB1EE_mvVtx,gMaxnphos);
       resetFloat(fT_pho_Cone03PhotonIso_dR045EB070EE_dEta015_pt08EB1EE_mvVtx,gMaxnphos);
       resetFloat(fT_pho_Cone04PhotonIso_dR045EB070EE_dEta015_pt08EB1EE_mvVtx,gMaxnphos);
       resetFloat(fT_pho_Cone04NeutralHadronIso_dR0_dEta0_pt0,gMaxnphos);
       resetFloat(fT_pho_Cone04NeutralHadronIso_dR0_dEta0_pt5,gMaxnphos);
       resetFloat(fT_pho_Cone04NeutralHadronIso_dR0_dEta0_pt0_nocracks,gMaxnphos);
       resetFloat(fT_pho_Cone04NeutralHadronIso_dR0_dEta0_pt5_nocracks,gMaxnphos);
       resetFloat(fT_pho_Cone04NeutralHadronIso_dR7_dEta0_pt0,gMaxnphos);
       resetFloat(fT_pho_Cone04NeutralHadronIso_dR7_dEta0_pt5,gMaxnphos);
       resetFloat(fT_pho_Cone01NeutralHadronIso_dR0_dEta0_pt0_mvVtx,gMaxnphos);
       resetFloat(fT_pho_Cone02NeutralHadronIso_dR0_dEta0_pt0_mvVtx,gMaxnphos);
       resetFloat(fT_pho_Cone03NeutralHadronIso_dR0_dEta0_pt0_mvVtx,gMaxnphos);
       resetFloat(fT_pho_Cone04NeutralHadronIso_dR0_dEta0_pt0_mvVtx,gMaxnphos);
       resetFloat(fT_pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_dz0_old,gMaxnphos);
       resetFloat(fT_pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_PFnoPU_old,gMaxnphos);
       resetFloat(fT_pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_dz0_old,gMaxnphos);
       resetFloat(fT_pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_PFnoPU_old,gMaxnphos);
       resetFloat(fT_pho_Cone01ChargedHadronIso_dR0_dEta0_pt0_dz0,gMaxnphos);
       resetFloat(fT_pho_Cone01ChargedHadronIso_dR0_dEta0_pt0_dz1_dxy01,gMaxnphos);
       resetFloat(fT_pho_Cone01ChargedHadronIso_dR0_dEta0_pt0_PFnoPU,gMaxnphos);
       resetFloat(fT_pho_Cone01ChargedHadronIso_dR015_dEta0_pt0_dz0,gMaxnphos);
       resetFloat(fT_pho_Cone01ChargedHadronIso_dR015_dEta0_pt0_dz1_dxy01,gMaxnphos);
       resetFloat(fT_pho_Cone01ChargedHadronIso_dR015_dEta0_pt0_PFnoPU,gMaxnphos);
       resetFloat(fT_pho_Cone02ChargedHadronIso_dR0_dEta0_pt0_dz0,gMaxnphos);
       resetFloat(fT_pho_Cone02ChargedHadronIso_dR0_dEta0_pt0_dz1_dxy01,gMaxnphos);
       resetFloat(fT_pho_Cone02ChargedHadronIso_dR0_dEta0_pt0_PFnoPU,gMaxnphos);
       resetFloat(fT_pho_Cone02ChargedHadronIso_dR015_dEta0_pt0_dz0,gMaxnphos);
       resetFloat(fT_pho_Cone02ChargedHadronIso_dR015_dEta0_pt0_dz1_dxy01,gMaxnphos);
       resetFloat(fT_pho_Cone02ChargedHadronIso_dR015_dEta0_pt0_PFnoPU,gMaxnphos);
       resetFloat(fT_pho_Cone03ChargedHadronIso_dR0_dEta0_pt0_dz0,gMaxnphos);
       resetFloat(fT_pho_Cone03ChargedHadronIso_dR0_dEta0_pt0_dz1_dxy01,gMaxnphos);
       resetFloat(fT_pho_Cone03ChargedHadronIso_dR0_dEta0_pt0_PFnoPU,gMaxnphos);
       resetFloat(fT_pho_Cone03ChargedHadronIso_dR015_dEta0_pt0_dz0,gMaxnphos);
       resetFloat(fT_pho_Cone03ChargedHadronIso_dR015_dEta0_pt0_dz1_dxy01,gMaxnphos);
       resetFloat(fT_pho_Cone03ChargedHadronIso_dR015_dEta0_pt0_PFnoPU,gMaxnphos);
       resetFloat(fT_pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_dz0,gMaxnphos);
       resetFloat(fT_pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_dz1_dxy01,gMaxnphos);
       resetFloat(fT_pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_PFnoPU,gMaxnphos);
       resetFloat(fT_pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_dz0,gMaxnphos);
       resetFloat(fT_pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_dz1_dxy01,gMaxnphos);
       resetFloat(fT_pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_PFnoPU,gMaxnphos);
       */

       for (int i=0; i<gMaxnphos; i++) {
	 pho_conv_vtx[i]=TVector3();
	 pho_conv_refitted_momentum[i]=TVector3();
       }

       for (int i=0; i<gMaxnconv; i++) {
	 conv_vtx[i]=TVector3();
	 conv_refitted_momentum[i]=TVector3();
       }
       
       for (int i=0; i<gMaxngenvtx; i++) {
	 gv_pos[i]=TVector3();
	 gv_p3[i]=TVector3();
       }
       
       //       vAna->clear();

	resetFloat(fTSCraw,gMaxnSC);
	resetFloat(fTSCpre,gMaxnSC);
	resetFloat(fTSCenergy,gMaxnSC);
	resetFloat(fTSCeta,gMaxnSC);
	resetFloat(fTSCphi,gMaxnSC);
	resetFloat(fTSCsigmaPhi,gMaxnSC);
	resetFloat(fTSCsigmaEta,gMaxnSC);
	resetFloat(fTSCbrem,gMaxnSC);
	resetFloat(fTSCR9,gMaxnSC);
	resetFloat(fTSCcrackcorrseed,gMaxnSC);
	resetFloat(fTSCcrackcorr,gMaxnSC);
	resetFloat(fTSClocalcorrseed,gMaxnSC);
	resetFloat(fTSClocalcorr,gMaxnSC);
	resetFloat(fTSCcrackcorrseedfactor,gMaxnSC);
	resetFloat(fTSClocalcorrseedfactor,gMaxnSC);
	resetFloat(fTSCx,gMaxnSC);
	resetFloat(fTSCy,gMaxnSC);
	resetFloat(fTSCz,gMaxnSC);
	resetInt(fTSCNXtals,gMaxnSC);
	
	for (int i=0; i<gMaxnSC; i++) {
	  for (int j=0; j<gMaxnSCxtals; j++) {
	    fTSCxtalX[i][j]=-999;
	    fTSCxtalY[i][j]=-999;
	    fTSCxtalZ[i][j]=-999;
	    fTSCxtalEtaWidth[i][j]=-999;
	    fTSCxtalPhiWidth[i][j]=-999;
	    for (int k=0; k<4; k++){                                                                                                                                            
	      fTSCxtalfrontX[i][j][k]=-999;
	      fTSCxtalfrontY[i][j][k]=-999;
	      fTSCxtalfrontZ[i][j][k]=-999;
	    } 
	  }
	}

	fTTrkPtSumx          = -999.99;
	fTTrkPtSumy          = -999.99;
	fTTrkPtSum           = -999.99;
	fTTrkPtSumphi        = -999.99;
	fTSumEt              = -999.99;
	fTECALSumEt          = -999.99;
	fTHCALSumEt          = -999.99;
	fTECALEsumx          = -999.99;
	fTECALEsumy          = -999.99;
	fTECALEsumz          = -999.99;
	fTECALMET            = -999.99;
	fTECALMETphi         = -999.99;
	fTECALMETeta         = -999.99;
	fTHCALEsumx          = -999.99;
	fTHCALEsumy          = -999.99;
	fTHCALEsumz          = -999.99;
	fTHCALMET            = -999.99;
	fTHCALMETphi         = -999.99;
	fTHCALMETeta         = -999.99;
	fTRawMET             = -999.99;
	fTRawMETpx           = -999.99;
	fTRawMETpy           = -999.99;
	fTRawMETphi          = -999.99;
	fTRawMETemEtFrac     = -999.99;
	fTRawMETemEtInEB     = -999.99;
	fTRawMETemEtInEE     = -999.99;
	fTRawMETemEtInHF     = -999.99;
	fTRawMEThadEtFrac    = -999.99;
	fTRawMEThadEtInHB    = -999.99;
	fTRawMEThadEtInHE    = -999.99;
	fTRawMEThadEtInHF    = -999.99;
	fTRawMETSignificance = -999.99;
	fTGenMET             = -999.99;
	fTGenMETpx           = -999.99;
	fTGenMETpy           = -999.99;
	fTGenMETphi          = -999.99;
	fTTCMET              = -999.99;
	fTTCMETpx            = -999.99;
	fTTCMETpy            = -999.99;
	fTTCMETphi           = -999.99;
	fTTCMETSignificance  = -999.99;
	fTMuJESCorrMET       = -999.99;
	fTMuJESCorrMETpx     = -999.99;
	fTMuJESCorrMETpy     = -999.99;
	fTMuJESCorrMETphi    = -999.99;
	fTPFMET              = -999.99;
	fTPFMETpx            = -999.99;
	fTPFMETpy            = -999.99;
	fTPFMETphi           = -999.99;
	fTPFMETSignificance  = -999.99;
        fTPFSumEt            = -999.99;
	fTPFMETPAT           = -999.99;
	fTPFMETPATpx         = -999.99;
	fTPFMETPATpy         = -999.99;
	fTPFMETPATphi        = -999.99;
	fTPFMETPATSignificance= -999.99;
	fTMETR12             = -999.99;
	fTMETR21             = -999.99;
}

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

// Method for matching of jets
const int NTupleProducer::matchJet(const reco::Jet* jet){
	// match to already filled genjets in ntuple:
	// this can only be called AFTER the genjets have been stored!
	if(fIsRealData || fTNGenJets < 1){
		edm::LogWarning("NTP") << "@SUB=matchJet"
			<< "Trying to access generator info on real data...";
		return -1;
	}

	// Try to match the reco jet to a stored generator jet
	double mindr(999.99);
	int matchedindex = -1;
	for(int i = 0; i < fTNGenJets; i++){

		// Restrict to cone of 0.1 in DR around candidate
		double dr = reco::deltaR(fTGenJetEta[i], fTGenJetPhi[i], jet->eta(), jet->phi());
		if(dr > 0.3) continue;

		// Restrict to pt match within a factor of 2
		double ndpt = fabs(fTGenJetPt[i] - jet->pt())/fTGenJetPt[i];
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
void NTupleProducer::ElectronDuplicate(std::vector<const SuperCluster*> elecPtr, std::vector<const GsfTrack*> trckPtr) {
	// Looks for duplication among electrons
	if( fTneles <= 0 ) return;

	// loop over the electrons
	for (int i = 0; i < fTneles; ++i) {
		const SuperCluster* supercluster = elecPtr[i];
		const GsfTrack* eletrack = trckPtr[i];

		// loop over the electrons again
		for (int j = i+1; j < fTneles; ++j) {
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

void NTupleProducer::PhotonElectronDuplicate(std::vector<const SuperCluster*> elecPtr, std::vector<const SuperCluster*> phoPtr) {
	// Looks for duplication between photons and electrons
	if( fTneles <= 0 ) return;
	if( fTnphotons <= 0 ) return;

	// loop over the photons
	for( int i = 0; i < fTnphotons; ++i ){
		const SuperCluster* phoSC = phoPtr[i];

		// loop over the electrons again
		for( int j = 0; j < fTneles; ++j ){
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

void NTupleProducer::ElJetOverlap(std::vector<const Jet*> jets, std::vector<const SuperCluster*> electrons, edm::Handle<CaloTowerCollection> calotowers){
	// checks for jets made from electrons
	// jetIndex and elecIndex contain the indices of the selected jets and
	//   electrons in the Views
	// (electrons and jets should be filled in the ntuple before checking)
	if (fTnjets <= 0) return;
	if (fTneles <= 0) return;

	std::vector<CaloTowerPtr> jetCaloRefs;

	// loop over the jets
	for (int i = 0; i < fTnjets; ++i) {

		// Collect the CaloTowers detIds for the jet
		const CaloJet* theJet = static_cast<const CaloJet*>(&(*jets[i]));
		jetCaloRefs = theJet->getCaloConstituents();

		// loop over the electrons
		for( int j = 0; j < fTneles; ++j ){
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

void NTupleProducer::PhotonJetOverlap(std::vector<const Jet*> jets, std::vector<const SuperCluster*> superclusters, edm::Handle<CaloTowerCollection> calotowers){
	// checks for jets made from photons
	// (photons and jets should be filled in the ntuple before checking)
	if( fTnjets <= 0 ) return;
	if( fTnphotons <= 0 ) return;

	std::vector<CaloTowerPtr> jetCaloRefs;

	// loop over the jets
	for( int i = 0; i < fTnjets; ++i ){

		// Collect the CaloTowers detIds for the jet
		const CaloJet* theJet = static_cast<const CaloJet*>(&(*jets[i]));
		jetCaloRefs = theJet->getCaloConstituents();

		// loop over the photons
		for( int j = 0; j < fTnphotons; ++j ){
			const SuperCluster* theSC = superclusters[j];

			math::XYZVector sharedP(0., 0., 0.);
			bool isInJet = IsEMObjectInJet(theSC, jetCaloRefs, calotowers, &sharedP);
			// float sharedE = sqrt(sharedP.X()*sharedP.X() + sharedP.Y()*sharedP.Y() + sharedP.Z()*sharedP.Z() );
			if( isInJet ){
				// fTPhotIsInJet[j] = i;
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

	float dphimin = DeltaPhiSigned(phi, phimin);
	float dphimax = DeltaPhiSigned(phi, phimax);

	float pi    = 3.141592654;
	if (dphimin*dphimax <= 0. && fabs(dphimin-dphimax) < pi){return true;}
	else {return false;}

}

float NTupleProducer::DeltaPhiSigned(float v1, float v2){
// Computes the clockwise phi difference v1-v2
// v1, v2 = phi of object 1 and 2

	float pi    = 3.141592654;
	float twopi = 6.283185307;

	float diff = v2 - v1;
	if (diff >  pi){ diff -= twopi;}
	else if (diff < -pi){ diff += twopi;}
	return diff;

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

void NTupleProducer::resetDouble(double *v, unsigned int size){
	for(size_t i = 0; i < size; ++i){
		v[i] = -999.99;
	}
}

void NTupleProducer::resetFloat(float *v, unsigned int size){
	for(size_t i = 0; i < size; ++i){
		v[i] = -999.99;
	}
}

void NTupleProducer::resetInt(int *v, unsigned int size){
	for(size_t i = 0; i < size; ++i){
		v[i] = -999;
	}
}

//double NTupleProducer::DeltaR(double phi1, double phi2, double eta1, double eta2){
//
//  double dphi = reco::deltaPhi(phi1,phi2);
//  double dR=sqrt(dphi*dphi+(eta2-eta1)*(eta2-eta1));
//
//  return dR;
//}
/*
void NTupleProducer::FillPhotonIsoVariables(double photonEta, double photonPhi, double photonVz, int type, bool isPU, edm::Handle<reco::PFCandidateCollection>& pfCandidates, int ipf, int phoqi){

  double pt = (*pfCandidates)[ipf].pt();
  double dEta = fabs(photonEta - (*pfCandidates)[ipf].eta());
  double dPhi = DeltaPhi(photonPhi,(*pfCandidates)[ipf].phi());
  double dR = sqrt(dEta*dEta+dPhi*dPhi);
  //  double dz = fabs(photonVz - (*pfCandidates)[ipf].vz());


  //cout << "FillPhotonIsoVariables pt="<<pt<<endl;
  if (type==0){

    fT_pho_Cone04NeutralHadronIso_dR0_dEta0_pt0[phoqi] += pt;
    if (pt>0.5) fT_pho_Cone04NeutralHadronIso_dR0_dEta0_pt5[phoqi] += pt;
  
    if (dR>0.07) {
      fT_pho_Cone04NeutralHadronIso_dR7_dEta0_pt0[phoqi] += pt;
      if (pt>0.5) fT_pho_Cone04NeutralHadronIso_dR7_dEta0_pt5[phoqi] += pt;
    }
    if (isInEtaCracks((*pfCandidates)[ipf].eta())==false && isInPhiCracks((*pfCandidates)[ipf].phi(),(*pfCandidates)[ipf].eta())==false){
      fT_pho_Cone04NeutralHadronIso_dR0_dEta0_pt0_nocracks[phoqi] += pt;
      if (pt>0.5) fT_pho_Cone04NeutralHadronIso_dR0_dEta0_pt5_nocracks[phoqi] += pt;
    }
  }
  
  if (type==1) { //Charged hadron
    fT_pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_dz0_old[phoqi] += pt;

    if (isPU==0) fT_pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_PFnoPU_old[phoqi] += pt;
    if (dR>0.015) {
      fT_pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_dz0_old[phoqi] += pt;
      if (isPU==0) fT_pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_PFnoPU_old[phoqi] += pt;
    }
  }

  if (type==2) { //Photon
    fT_pho_Cone04PhotonIso_dR0_dEta0_pt0[phoqi] += pt;
    if (pt>0.5) fT_pho_Cone04PhotonIso_dR0_dEta0_pt5[phoqi] += pt;
    if (dR>0.08) {
      fT_pho_Cone04PhotonIso_dR8_dEta0_pt0[phoqi] += pt;
      if (pt>0.5) fT_pho_Cone04PhotonIso_dR8_dEta0_pt5[phoqi] += pt;
    }
  }
  
  return;
}
*/
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

int NTupleProducer::FindPFCandType(int id){

  int type = -1;

  if (id==111 || id==130 || id==310 || id==2112) type=0; //neutral hadrons
  if (fabs(id)==211 || fabs(id)==321 || id==999211 || fabs(id)==2212) type=1; //charged hadrons
  if (id==22) type=2; //photons
  if (fabs(id)==11) type=3; //electrons
  if (fabs(id)==13) type=4; //muons

  return type;
}

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
				  TVector3(fTPhotcaloPosX[p1],fTPhotcaloPosY[p1],fTPhotcaloPosZ[p1]),
				  TVector3(fTbeamspotx,fTbeamspoty,fTbeamspotz),
				  conv_vtx[iConv1],
				  conv_refitted_momentum[iConv1],
				  fTPhotEnergy[p1],
				  fTPhotisEB[p1],
				  conv_ntracks[iConv1],
				  conv_validvtx[iConv1],
				  conv_chi2_probability[iConv1],
				  conv_eoverp[iConv1]
			);
	} 

	
	return PhotonInfo(p1, 
			  TVector3(fTPhotcaloPosX[p1],fTPhotcaloPosY[p1],fTPhotcaloPosZ[p1]),
			  TVector3(fTbeamspotx,fTbeamspoty,fTbeamspotz),
			  pho_conv_vtx[p1],
			  pho_conv_refitted_momentum[p1],
			  fTPhotEnergy[p1],
			  fTPhotisEB[p1],
			  pho_conv_ntracks[p1],                                                                                                                             
			  pho_conv_validvtx[p1],                                                                                                                            
			  pho_conv_chi2_probability[p1] ,                                                                                                                   
			  pho_conv_eoverp[p1]                                                                                                                               
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
        for(int i=0; i<fTnvrtx ; i++) {
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
	  for(int i=0; i < fTnvrtx; i++) {
	    TVector3 * vtxpos= new TVector3(fTvrtxx[i],fTvrtxy[i],fTvrtxz[i]);
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
	
	TLorentzVector dipho = TLorentzVector(fTPhotPx[p1],fTPhotPy[p1],fTPhotPz[p1],fTPhotEnergy[p1]) + TLorentzVector(fTPhotPx[p2],fTPhotPy[p2],fTPhotPz[p2],fTPhotEnergy[p2]);
	
	unsigned int nbest ;
	if (  dipho.Pt() < 30 ) nbest = 5;
	else nbest = 3; 
	if (rankprodAll.size() < nbest ) nbest = rankprodAll.size();

	for (unsigned int ii = 0; ii < nbest; ii++ ){
	  TVector3 * vtxpos= new TVector3(fTvrtxx[rankprodAll[ii]],fTvrtxy[rankprodAll[ii]],fTvrtxz[rankprodAll[ii]]);
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
  
  float sc_eta  = TVector3(fTPhotcaloPosX[p1],fTPhotcaloPosY[p1],fTPhotcaloPosZ[p1]).Eta();
  float  phi  = TVector3(fTPhotcaloPosX[p1],fTPhotcaloPosY[p1],fTPhotcaloPosZ[p1]).Phi();
  double sc_phi = phiNorm(phi);
  
  TLorentzVector * p4 = new TLorentzVector(fTPhotPx[p1],fTPhotPy[p1],fTPhotPz[p1],fTPhotEnergy[p1]);
  float et = fTPhotPt[p1];
  
  float detaMin=999.;
  float dphiMin=999.;   
  float dRMin = 999.;

  float mconv_pt=-999999;
  int iMatch=-1;     
  float conv_pt = -9999;

  // if(LDEBUG)  cout << "   LoopAll::matchPhotonToConversion conv_n " << conv_n << endl; 
  for(int iconv=0; iconv<conv_n; iconv++) {
    TVector3 refittedPairMomentum= conv_refitted_momentum[iconv];
    conv_pt =  refittedPairMomentum.Pt();
    if (conv_pt < 1 ) continue;    
    if ( !conv_validvtx[iconv] || conv_ntracks[iconv]!=2 || conv_chi2_probability[iconv]<0.000001) continue;

    phi  = conv_refitted_momentum[iconv].Phi();
    conv_phi  = phiNorm(phi);
    float eta  = conv_refitted_momentum[iconv].Eta();
    conv_eta = etaTransformation(eta, conv_zofprimvtxfromtrks[iconv] );

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

bool NTupleProducer::tkIsHighPurity(reco::TrackRef tk) const { return ( tk->qualityMask() & (1<<2) ) >> 2; }
bool NTupleProducer::TrackCut(reco::TrackRef tk) const { return false; }
bool NTupleProducer::ConversionsCut(const reco::Conversion &conv) { 
  //  return (sqrt(conv.refittedPairMomentum().perp2()) <  0); 
  return false;
}

double NTupleProducer::phiNorm(float &phi) {

  const float pi = 3.1415927;
  const float twopi = 2.0*pi;

  if(phi >  pi) {phi = phi - twopi;}
  if(phi < -pi) {phi = phi + twopi;}

  return phi;
}


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


double NTupleProducer::GenPartonicIso_allpart(const reco::GenParticle & photon,    edm::Handle <reco::GenParticleCollection> & genparticles, double dRcone){

  double etsum=0;
  double dR=0;

  for (unsigned int ipart=0; ipart< genparticles->size(); ipart++){
    const reco::GenParticle & mygenpart = (*genparticles)[ ipart ];

    if (mygenpart.status()==1){

//      bool isDoubleCounted = false;
//      for (uint i=0; i<genparticleToVeto.size(); i++) {
//        if (genparticleToVeto[i]->pt()==(*genparticles)[ ipart ].pt()) isDoubleCounted=true;
//      }
//
//      if (isDoubleCounted==true) continue;

      if (mygenpart.pdgId()!=22 || (fabs(mygenpart.pt()-photon.pt())>0.01 && mygenpart.pdgId()==22)){
        dR = reco::deltaR(photon.eta(), photon.phi(), mygenpart.eta(), mygenpart.phi());
        if (dR<dRcone && dR>1e-05){
          etsum += mygenpart.et();
        }
      }
    }
  }

  //if (etsum>0) cout << "GenPartonicIso etsum="<<etsum<<endl;

  return etsum;

}




ETHVertexInfo::ETHVertexInfo(int nvtx, float * vtxx, float * vtxy, float * vtxz, 
				 int ntracks, float * tkpx, float * tkpy, float * tkpz,
				 float * tkPtErr, int * tkVtxId,
				 float * tkd0, float * tkd0Err, float * tkdz, float * tkdzErr,
			     bool * tkIsHighPurity, std::vector<std::vector<unsigned short> > vtx_std_tkind, std::vector<std::vector<float> > vtx_std_tkweight, int * vtx_std_ntks
) :
	nvtx_(nvtx),
	vtxx_(vtxx),
	vtxy_(vtxy),
	vtxz_(vtxz),
	ntracks_(ntracks),
	tkpx_(tkpx),
	tkpy_(tkpy),
	tkpz_(tkpz),
	tkPtErr_(tkPtErr),
	tkVtxId_(tkVtxId),
	tkd0_(tkd0),
	tkd0Err_(tkd0Err),
	tkdz_(tkdz),
	tkdzErr_(tkdzErr),
	tkIsHighPurity_(tkIsHighPurity),
	vtx_std_tkind_(vtx_std_tkind),
	vtx_std_tkweight_(vtx_std_tkweight),
	vtx_std_ntks_(vtx_std_ntks)
{

  for (unsigned int i=0; i<vtx_std_tkind_.size();i++){
    for (unsigned int j=0; j<vtx_std_tkind_.at(i).size();j++){
      vtx_std_tkind_helper_[i][j]=vtx_std_tkind_.at(i).at(j);
    }
  }

  for (unsigned int i=0; i<vtx_std_tkweight_.size();i++){
    for (unsigned int j=0; j<vtx_std_tkweight_.at(i).size();j++){
      vtx_std_tkweight_helper_[i][j]=vtx_std_tkweight_.at(i).at(j);
    }
  }

}
  
ETHVertexInfo::~ETHVertexInfo() {}








//define this as a plug-in
DEFINE_FWK_MODULE(NTupleProducer);
