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
// $Id: NTupleProducer.cc,v 1.128 2011/08/29 13:43:37 buchmann Exp $
//
//


// system include files
#include <memory>
#include <iostream>
#include <cmath>

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

#include "FWCore/Framework/interface/ESHandle.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"

#include "DataFormats/METReco/interface/BeamHaloSummary.h"

#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"

/*
#include "DataFormats/AnomalousEcalDataFormats/interface/AnomalousECALVariables.h"
#include "PhysicsTools/EcalAnomalousEventFilter/interface/EcalBoundaryInfoCalculator.h"
*/

#include "MagneticField/Engine/interface/MagneticField.h"

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
	if(!fIsModelScan) fHBHENoiseResultTag = iConfig.getUntrackedParameter<edm::InputTag>("tag_hcalnoise");
	fSrcRho             = iConfig.getUntrackedParameter<edm::InputTag>("tag_srcRho");

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
	fMingenleptpt   = iConfig.getParameter<double>("sel_mingenleptpt");
	fMaxgenlepteta  = iConfig.getParameter<double>("sel_maxgenlepteta");
	fMingenjetpt    = iConfig.getParameter<double>("sel_mingenjetpt");
	fMaxgenjeteta   = iConfig.getParameter<double>("sel_maxgenjeteta");
	fMinebrechitE   = iConfig.getParameter<double>("sel_fminebrechitE");

	
	if(fIsModelScan) {
		LHAPDF::initPDFSet("cteq66.LHgrid",1);
		NPdfs = LHAPDF::numberPDF();
	}

	fBtagMatchdeltaR = iConfig.getParameter<double>("btag_matchdeltaR"); // 0.25

	// Create histograms and trees
	fHhltstat    = fTFileService->make<TH1I>("HLTTriggerStats",    "HLTTriggerStatistics",    gMaxhltbits+2,    0, gMaxhltbits+2);
	fHl1physstat = fTFileService->make<TH1I>("L1PhysTriggerStats", "L1PhysTriggerStatistics", gMaxl1physbits+2, 0, gMaxl1physbits+2);
	fHl1techstat = fTFileService->make<TH1I>("L1TechTriggerStats", "L1TechTriggerStatistics", gMaxl1techbits+2, 0, gMaxl1techbits+2);
        fHpileupstat = fTFileService->make<TH1I>("PileUpStats", "PileUpStats", 40, 0, 40 ); // Keep track of pileup distribution
	fRunTree     = fTFileService->make<TTree>("RunInfo", "ETHZRunAnalysisTree");
	fEventTree   = fTFileService->make<TTree>("Analysis", "ETHZAnalysisTree");

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

	// ECAL dead cell Trigger Primitive filter
	edm::Handle<bool> EcalDeadTPFilterFlag;
	iEvent.getByLabel("ecalDeadCellTPfilter",EcalDeadTPFilterFlag);
	fTecalDeadTPFilterFlag = (int) *EcalDeadTPFilterFlag;
	

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
	edm::Handle<bool> hbHeNoiseFlag;
	if(!fIsModelScan) {
	   iEvent.getByLabel(fHBHENoiseResultTag,hbHeNoiseFlag);
	   fTHBHENoiseFlag = static_cast<int>(*hbHeNoiseFlag);
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
     
		iEvent.getByLabel("addPileupInfo", pileupInfo);
		std::vector<PileupSummaryInfo>::const_iterator PVI;

		for (PVI = pileupInfo->begin(); PVI !=pileupInfo->end(); ++PVI){
		  if( PVI->getBunchCrossing() == 0 ){ // in-time PU
		    fTpuNumInteractions  = PVI->getPU_NumInteractions();
		    fHpileupstat->Fill( fTpuNumInteractions );
		    
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

			for(int pdf=1; pdf<= NPdfs; pdf++){
				LHAPDF::initPDF(pdf);
				double newpdf1 = LHAPDF::xfx(x1, Q, id1)/x1;
				double newpdf2 = LHAPDF::xfx(x2, Q, id2)/x2;
				fTpdfW[pdf] = newpdf1/newpdf1_0*newpdf2/newpdf2_0;
			}

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

			if(!( g_part->status() ==1 || (g_part->status() ==2 && abs(g_part->pdgId())==5))) continue;

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
						<< "Maximum number of gen-leptons exceeded..";
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
		if(fTnEBhits>gMaxnEBhits)
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

	// Get electrons, order them by pt and apply selection
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
		if(ip->pt() < fMinphopt) continue;
		if(fabs(ip->eta()) > fMaxphoeta) continue;

		phoqi++; // Count how many we'll eventually store
		phoOrdered.push_back(make_pair(phoIndex,ip->pt()));
	}
	std::sort(phoOrdered.begin(),phoOrdered.end(),indexComparator);
	fTnphotons = phoOrdered.size();
	phoqi = 0;

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
	
	// EcalClusterLazyTools *lazyTools = new EcalClusterLazyTools(iEvent, iSetup, edm::InputTag("reducedEcalRecHitsEB"), edm::InputTag("reducedEcalRecHitsEE"));

	fTPhotSCEnergy[phoqi]       = photon.superCluster()->rawEnergy();
	fTPhotSCEtaWidth[phoqi]     = photon.superCluster()->etaWidth();
// DISABLED: NO SEED IN AOD (UPDATE IT IN 4_2)
// 	fTPhotSCSigmaPhiPhi[phoqi]  = lazyTools->covariances(*(photon.superCluster()->seed())).at(2);
	fTPhotHasPixSeed[phoqi]     = photon.hasPixelSeed() ? 1:0;
	fTPhotHasConvTrks[phoqi]    = photon.hasConversionTracks() ? 1:0;

	// fTPhotIsInJet[phoqi]      = -1;
	// fTPhotDupEl[phoqi]        = -1;
	// fTPhotSharedPx[phoqi]     = 0.;
	// fTPhotSharedPy[phoqi]     = 0.;
	// fTPhotSharedPz[phoqi]     = 0.;
	// fTPhotSharedEnergy[phoqi] = 0.;
	fTgoodphoton[phoqi]       = 0;
	fTPhotIsIso[phoqi]        = 1;

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
 	}

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
		fTjChHadFrac    [jqi] = jet->chargedHadronEnergy()/uncorr_energy;
		fTjNeuHadFrac   [jqi] = jet->neutralHadronEnergy()/uncorr_energy + jet->HFHadronEnergy()/uncorr_energy;
		fTjChEmFrac     [jqi] = jet->chargedEmEnergy()/uncorr_energy;
		fTjNeuEmFrac    [jqi] = jet->neutralEmEnergy()/uncorr_energy;
		fTjChMuEFrac    [jqi] = jet->chargedMuEnergy()/uncorr_energy;

		// Calculate the DR wrt the closest electron
		float ejDRmin = 10.; // Default when no electrons previously selected
		for( int j = 0; j < fTneles; j++ ){
			float ejDR = reco::deltaR(jet->eta(), jet->phi(), fTeeta[j], fTephi[j]);
			if(ejDR<ejDRmin) ejDRmin = ejDR;
		}
		fTjeMinDR[jqi] = ejDRmin;

		// B-tagging probability (for 4 b-taggings), with delta R matching for now
		float mindr(999.99);
		for( size_t i = 0; i < jetsAndProbsTkCntHighEff->size(); i++ ){
			float deltar = reco::deltaR( jet->eta(), jet->phi(), (*jetsAndProbsTkCntHighEff)[i].first->eta(), (*jetsAndProbsTkCntHighEff)[i].first->phi());
			if( deltar <= fBtagMatchdeltaR && deltar < mindr)  {
				fTjbTagProbTkCntHighEff[jqi] = (*jetsAndProbsTkCntHighEff)[i].second;
				mindr = deltar;
			}
		}
		mindr = 999.99;
		for( size_t i = 0; i < jetsAndProbsTkCntHighPur->size(); i++ ){
			float deltar = reco::deltaR( jet->eta(), jet->phi(), (*jetsAndProbsTkCntHighPur)[i].first->eta(), (*jetsAndProbsTkCntHighPur)[i].first->phi());
			if( deltar <= fBtagMatchdeltaR && deltar < mindr)  {
				fTjbTagProbTkCntHighPur[jqi] = (*jetsAndProbsTkCntHighPur)[i].second;
				mindr = deltar;
			}
		}
		mindr = 999.99;
		for( size_t i = 0; i < jetsAndProbsSimpSVHighEff->size(); i++ ){
			float deltar = reco::deltaR( jet->eta(), jet->phi(), (*jetsAndProbsSimpSVHighEff)[i].first->eta(), (*jetsAndProbsSimpSVHighEff)[i].first->phi());
			if( deltar <= fBtagMatchdeltaR && deltar < mindr)  {
				fTjbTagProbSimpSVHighEff[jqi] = (*jetsAndProbsSimpSVHighEff)[i].second;
				mindr = deltar;
			}
		}
		mindr = 999.99;
		for( size_t i = 0; i < jetsAndProbsSimpSVHighPur->size(); i++ ){
			float deltar = reco::deltaR( jet->eta(), jet->phi(), (*jetsAndProbsSimpSVHighPur)[i].first->eta(), (*jetsAndProbsSimpSVHighPur)[i].first->phi());
			if( deltar <= fBtagMatchdeltaR && deltar < mindr)  {
				fTjbTagProbSimpSVHighPur[jqi] = (*jetsAndProbsSimpSVHighPur)[i].second;
				mindr = deltar;
			}
		}

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
	fTPFMETSignificance = (pfmet->at(0)).significance();
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

	if(!fIsRealData&&fIsModelScan) {
		Handle<LHEEventProduct> product;
		iEvent.getByLabel("source", product);

		LHEEventProduct::comments_const_iterator c_begin = product->comments_begin();
		LHEEventProduct::comments_const_iterator c_end = product->comments_end();

		float mGL;
		float mLSP;
	        float xCHI;
	        for(LHEEventProduct::comments_const_iterator cit=c_begin; cit!=c_end; ++cit) {
	          size_t found = (*cit).find("model");
	          //# model T5zz_0.5_925.0_400.0to1000.0_450.0.lhe
	          if( found != std::string::npos)   {
		    // Simplified models
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

//std::cout << "m0=" << m0 << " m12=" << m12 << " signMu=" << signMu << " and A0=" << A0 << std::endl;
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
	fEventTree->Branch("Event"            ,&fTeventnumber     ,"Event/I");
	fEventTree->Branch("LumiSection"      ,&fTlumisection     ,"LumiSection/I");
	fEventTree->Branch("PtHat"            ,&fTpthat           ,"PtHat/F");
	fEventTree->Branch("SigProcID"        ,&fTsigprocid       ,"SigProcID/I");
	fEventTree->Branch("PDFScalePDF"      ,&fTpdfscalePDF     ,"PDFScalePDF/F");
	fEventTree->Branch("PDFID1"           ,&fTpdfid1          ,"PDFID1/I");
	fEventTree->Branch("PDFID2"           ,&fTpdfid2          ,"PDFID2/I");
	fEventTree->Branch("PDFx1"            ,&fTpdfx1           ,"PDFx1/F");
	fEventTree->Branch("PDFx2"            ,&fTpdfx2           ,"PDFx2/F");
	fEventTree->Branch("PDFxPDF1"         ,&fTpdfxPDF1        ,"PDFxPDF1/F");
	fEventTree->Branch("PDFxPDF2"         ,&fTpdfxPDF2        ,"PDFxPDF2/F");
	fEventTree->Branch("ExtXSecLO"        ,&fTextxslo         ,"ExtXSecLO/F");
	fEventTree->Branch("IntXSec"          ,&fTintxs           ,"IntXSec/F");
	fEventTree->Branch("pdfW"             ,&fTpdfW             ,"pdfW[100]/F");
	fEventTree->Branch("NPdfs"            ,&NPdfs              ,"NPdfs/I");

	// Pile-Up information:
	fEventTree->Branch("PUnumInteractions",   &fTpuNumInteractions   ,"PUnumInteractions/I");
	fEventTree->Branch("PUnumFilled",&fTpuNumFilled,"PUnumFilled/I");
	fEventTree->Branch("PUOOTnumInteractionsEarly",&fTpuOOTNumInteractionsEarly,"PUOOTnumInteractionsEarly/I");
	fEventTree->Branch("PUOOTnumInteractionsLate",&fTpuOOTNumInteractionsLate,"PUOOTnumInteractionsLate/I");
	fEventTree->Branch("PUzPositions"     ,&fTpuZpositions     ,"PUzPositions[PUnumFilled]/F");
	fEventTree->Branch("PUsumPtLowPt"     ,&fTpuSumpT_lowpT    ,"PUsumPtLowPt[PUnumFilled]/F");
	fEventTree->Branch("PUsumPtHighPt"    ,&fTpuSumpT_highpT   ,"PUsumPtHighPt[PUnumFilled]/F");
	fEventTree->Branch("PUnTrksLowPt"     ,&fTpuNtrks_lowpT    ,"PUnTrksLowPt[PUnumFilled]/F");
	fEventTree->Branch("PUnTrksHighPt"    ,&fTpuNtrks_highpT   ,"PUnTrksHighPt[PUnumFilled]/F");
	fEventTree->Branch("Rho"              ,&fTrho              ,"Rho/F");
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
        fEventTree->Branch("M0"               ,&fTSUSYScanM0      ,"M0/F");
        fEventTree->Branch("M12"              ,&fTSUSYScanM12     ,"M12/F");
        fEventTree->Branch("signMu"           ,&fTSUSYScanMu      ,"signMu/F");
        fEventTree->Branch("A0"               ,&fTSUSYScanA0      ,"A0/F");

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
	fEventTree->Branch("MaxGenJetExceed"  ,&fTflagmaxgenjetexc  ,"MaxGenJetExceed/I");
	fEventTree->Branch("MaxVerticesExceed",&fTflagmaxvrtxexc    ,"MaxVerticesExceed/I");
	fEventTree->Branch("HBHENoiseFlag"    ,&fTHBHENoiseFlag     ,"HBHENoiseFlag/I");
	fEventTree->Branch("CSCTightHaloID"   ,&fTcscTightHaloID    ,"CSCTightHaloID/I");
	fEventTree->Branch("EcalDeadTPFilterFlag",&fTecalDeadTPFilterFlag,"EcalDeadTPFilterFlag/I");
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
	fEventTree->Branch("PhoHasConvTrks"   ,&fTPhotHasConvTrks   ,"PhoHasConvTrks[NPhotons]/I");
	// fEventTree->Branch("PhoIsInJet"       ,&fTPhotIsInJet       ,"PhoIsInJet[NPhotons]/I");
	// fEventTree->Branch("PhoIsElDupl"      ,&fTPhotDupEl         ,"PhoIsElDupl[NPhotons]/I");
	// fEventTree->Branch("PhoSharedPx"      ,&fTPhotSharedPx      ,"PhoSharedPx[NPhotons]/F");
	// fEventTree->Branch("PhoSharedPy"      ,&fTPhotSharedPy      ,"PhoSharedPy[NPhotons]/F");
	// fEventTree->Branch("PhoSharedPz"      ,&fTPhotSharedPz      ,"PhoSharedPz[NPhotons]/F");
	// fEventTree->Branch("PhoSharedEnergy"  ,&fTPhotSharedEnergy  ,"PhoSharedEnergy[NPhotons]/F");
	fEventTree->Branch("PhoScSeedSeverity",&fTPhotScSeedSeverity,"PhoScSeedSeverity[NPhotons]/I");
	fEventTree->Branch("PhoE1OverE9"      ,&fTPhotE1OverE9      ,"PhoE1OverE9[NPhotons]/F");
	fEventTree->Branch("PhoS4OverS1"      ,&fTPhotS4OverS1      ,"PhoS4OverS1[NPhotons]/F");

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
	fTsigprocid         = -999;
	fTpdfscalePDF       = -999.99;
	fTpdfid1            = -999;
	fTpdfid2            = -999;
	fTpdfx1             = -999.99;
	fTpdfx2             = -999.99;
	fTpdfxPDF1          = -999.99;
	fTpdfxPDF2          = -999.99;

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
	fTcscTightHaloID    = -999;
	fTrho               = -999;

	fTecalDeadTPFilterFlag = -999;
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
	fTflagmaxgenjetexc = 0;
	fTflagmaxvrtxexc    = 0;

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
	fTnphotons    = 0;
	fTnphotonstot = 0;
	fTngenleptons = 0;
	fTNGenJets    = 0;
	fTnvrtx       = 0;

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

//define this as a plug-in
DEFINE_FWK_MODULE(NTupleProducer);
