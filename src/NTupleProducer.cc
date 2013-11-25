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
// $Id: NTupleProducer.cc,v 1.146.2.51 2013/02/11 10:47:15 fronga Exp $
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
#include "FWCore/Utilities/interface/Exception.h"

// Utilities
#include "CommonTools/Utils/interface/PtComparator.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
// this is for 3d-ip testing at the moment. marc.
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "RecoVertex/GaussianSumVertexFit/interface/GsfVertexTrackCompatibilityEstimator.h"
// done test

#include "RecoLocalCalo/EcalRecAlgos/interface/EcalSeverityLevelAlgo.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
#include "CondFormats/DataRecord/interface/EcalChannelStatusRcd.h"
#include "CondFormats/EcalObjects/interface/EcalChannelStatus.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"

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
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
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
#include "EGamma/EGammaAnalysisTools/src/PFIsolationEstimator.cc"

/*
  #include "DataFormats/AnomalousEcalDataFormats/interface/AnomalousECALVariables.h"
  #include "PhysicsTools/EcalAnomalousEventFilter/interface/EcalBoundaryInfoCalculator.h"
*/

#include "MagneticField/Engine/interface/MagneticField.h"

#include "Geometry/Records/interface/CaloTopologyRecord.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h"
#include "RecoParticleFlow/PFClusterTools/interface/ClusterClusterMapping.h"

#include "SimDataFormats/JetMatching/interface/JetFlavour.h"
#include "SimDataFormats/JetMatching/interface/JetFlavourMatching.h"
#include "SimDataFormats/JetMatching/interface/MatchedPartons.h"
#include "SimDataFormats/JetMatching/interface/JetMatchedPartons.h"

#include "DiLeptonAnalysis/NTupleProducer/interface/ETHVertexInfo.h"
#include "PFIsolation/SuperClusterFootprintRemoval/interface/SuperClusterFootprintRemoval.h"

#include "CMGTools/External/interface/PileupJetIdentifier.h"

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
  fIsRealData = iConfig.getParameter<bool>("isRealData");
  fIsModelScan = iConfig.getParameter<bool>("isModelScan");
  fIsFastSim = iConfig.getParameter<bool>("isFastSim");

  if(fIsRealData&&fIsModelScan) fIsModelScan=false; // avoiding possible mistakes

  // InputTags
  fMuonTag             = iConfig.getParameter<edm::InputTag>("tag_muons");
  fMuonPfIsoTagsCustom = iConfig.getParameter<std::vector<edm::InputTag> >("tag_muonpfisosCustom");
  fElectronTag         = iConfig.getParameter<edm::InputTag>("tag_electrons");
  fElePfIsoTagsCustom  = iConfig.getParameter<std::vector<edm::InputTag> >("tag_elepfisosCustom");
  fElePfIsoTagsEvent   = iConfig.getParameter<std::vector<edm::InputTag> >("tag_elepfisosEvent");
  fMuIsoDepTkTag       = iConfig.getParameter<edm::InputTag>("tag_muisodeptk");
  fMuIsoDepECTag       = iConfig.getParameter<edm::InputTag>("tag_muisodepec");
  fMuIsoDepHCTag       = iConfig.getParameter<edm::InputTag>("tag_muisodephc");
  fJetTag              = iConfig.getParameter<edm::InputTag>("tag_jets");
  fJetCorrs            = iConfig.getParameter<std::string>("jetCorrs");
  fBtagTags            = iConfig.getParameter<std::vector<edm::InputTag> >("tag_btags");
  fPartonMatch         = iConfig.getParameter<edm::InputTag>("tag_partonmatch");
  fRawCaloMETTag       = iConfig.getParameter<edm::InputTag>("tag_rawcalomet");
  fTCMETTag            = iConfig.getParameter<edm::InputTag>("tag_tcmet");
  fPFMETTag            = iConfig.getParameter<edm::InputTag>("tag_pfmet");
  fCorrCaloMETTag      = iConfig.getParameter<edm::InputTag>("tag_corrcalomet");
  fGenMETTag           = iConfig.getParameter<edm::InputTag>("tag_genmet");
  fVertexTag           = iConfig.getParameter<edm::InputTag>("tag_vertex");
  fVertexTagWithBS     = iConfig.getParameter<edm::InputTag>("tag_vertex_withbs");
  fTrackTag            = iConfig.getParameter<edm::InputTag>("tag_tracks");
  fPhotonTag           = iConfig.getParameter<edm::InputTag>("tag_photons");
  fCalTowTag           = iConfig.getParameter<edm::InputTag>("tag_caltow");
  fEBRecHitsTag        = iConfig.getParameter<edm::InputTag>("tag_EBrechits");
  fEERecHitsTag        = iConfig.getParameter<edm::InputTag>("tag_EErechits");
  fGenPartTag          = iConfig.getParameter<edm::InputTag>("tag_genpart");
  fGenJetTag           = iConfig.getParameter<edm::InputTag>("tag_genjets");
  fL1TriggerTag        = iConfig.getParameter<edm::InputTag>("tag_l1trig");
  fHLTTrigEventTag     = iConfig.getParameter<edm::InputTag>("tag_hlttrigevent");
  fSrcRho              = iConfig.getParameter<edm::InputTag>("tag_srcRho");
  fSrcSigma            = iConfig.getParameter<edm::InputTag>("tag_srcSigma");
  fSrcRhoForIso        = iConfig.getParameter<edm::InputTag>("tag_srcRhoForIso");
  pfphotonsProducerTag = iConfig.getParameter<edm::InputTag>("tag_pfphotonsProducer");
  pfProducerTag        = iConfig.getParameter<edm::InputTag>("tag_pfProducer");
  fSCTagBarrel = iConfig.getParameter<edm::InputTag>("tag_SC_barrel");
  fSCTagEndcap = iConfig.getParameter<edm::InputTag>("tag_SC_endcap");
  doPhotonStuff   = iConfig.getParameter<bool>("tag_doPhotonStuff");
  if (fIsModelScan) doPhotonStuff=false;

  if (doPhotonStuff){
    fIsModelScan = false;
    fVertexTag = fVertexTagWithBS;
    fTrackCollForVertexing          = iConfig.getParameter<edm::InputTag>("tag_fTrackCollForVertexing");
    fAllConversionsCollForVertexing = iConfig.getParameter<edm::InputTag>("tag_fallConversionsCollForVertexing");
    regrVersion      = iConfig.getParameter<int>("tag_regressionVersion");
    jetMVAAlgos = iConfig.getParameter<std::vector<edm::ParameterSet> >("tag_puJetIDAlgos");
    QGSystString = iConfig.getParameter<std::string>("tag_QGSyst");
  }

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
  fMinSCrawPt     = iConfig.getParameter<double>("sel_minSCrawPt");
  fMaxPfCandEta   = iConfig.getParameter<double>("sel_maxpfcandeta");
  fMinEBRechitE   = iConfig.getParameter<double>("sel_fminebrechitE");

  fMinGenLeptPt   = iConfig.getParameter<double>("sel_mingenleptpt");
  fMaxGenLeptEta  = iConfig.getParameter<double>("sel_maxgenlepteta");
  fMinGenPhotPt   = iConfig.getParameter<double>("sel_mingenphotpt");
  fMaxGenPhotEta  = iConfig.getParameter<double>("sel_maxgenphoteta");
  fMinGenJetPt    = iConfig.getParameter<double>("sel_mingenjetpt");
  fMaxGenJetEta   = iConfig.getParameter<double>("sel_maxgenjeteta");

  ///commenting out, chrashes for unknow reasons
  /*
  if(fIsModelScan) {
    LHAPDF::initPDFSet("cteq66.LHgrid",1);
    std::cout << LHAPDF::numberPDF() << std::endl;
    *fTNPdfs.get() = 44;//(int)LHAPDF::numberPDF();
  }
  */
  CrackCorrFunc    = EcalClusterFunctionFactory::get()->create("EcalClusterCrackCorrection", iConfig);
  LocalCorrFunc    = EcalClusterFunctionFactory::get()->create("EcalClusterLocalContCorrection",iConfig);

  // Check on number of allowed tags
  if ( fBtagTags.size()>gMaxNBtags )
    throw cms::Exception("BadConfig") << "Too many (" << fBtagTags.size() << ") btagging algos requested. "
                                      << "Maximum is " << gMaxNBtags;
  if ( fMuonPfIsoTagsCustom.size()>gMaxNPfIsoTags )
    throw cms::Exception("BadConfig") << "Too many (" << fMuonPfIsoTagsCustom.size() << ") muon PF iso tags requested. "
                                      << "Maximum is " << gMaxNPfIsoTags;
  if ( fElePfIsoTagsCustom.size()>gMaxNPfIsoTags )
    throw cms::Exception("BadConfig") << "Too many (" << fElePfIsoTagsCustom.size() << ") electron PF iso tags requested. "
                                      << "Maximum is " << gMaxNPfIsoTags;
  if ( fElePfIsoTagsEvent.size()>gMaxNPfIsoTags )
    throw cms::Exception("BadConfig") << "Too many (" << fElePfIsoTagsEvent.size() << ") electron PF iso tags requested. "
                                      << "Maximum is " << gMaxNPfIsoTags;

  // Dump the full configuration FIXME: not needed with EDM and provenance...
  edm::LogVerbatim("NTP") << "---------------------------------";
  edm::LogVerbatim("NTP") << " ==> NTupleProducer Constructor ...";
  edm::LogVerbatim("NTP") << iConfig;

  if (!doPhotonStuff) {
  // Create additional jet fillers
  std::vector<edm::ParameterSet> jConfigs = iConfig.getParameter<std::vector<edm::ParameterSet> >("jets");
  for (size_t i=0; i<jConfigs.size(); ++i)
    if ( jConfigs[i].getParameter<bool>("isPat") ) jetFillers.push_back( new JetFillerPat(jConfigs[i], fIsRealData) );
    else jetFillers.push_back( new JetFillerReco(jConfigs[i], fIsRealData) );

  // Create additional lepton fillers
  std::vector<edm::ParameterSet> lConfigs = iConfig.getParameter<std::vector<edm::ParameterSet> >("leptons");
   for (size_t i=0; i<lConfigs.size(); ++i) {
    std::string type(lConfigs[i].getParameter<std::string>("type"));
    if ( type == "electron" ) 
      electronFillers.push_back( new PatElectronFiller(lConfigs[i], fIsRealData) );
    else if ( type == "muon" ) 
      muonFillers.push_back( new PatMuonFiller(lConfigs[i], fIsRealData) );
    else if ( type == "tau" ) 
      tauFillers.push_back( new PatTauFiller(lConfigs[i], fIsRealData) );
   }

  // Create PF candidate fillers
  std::vector<edm::ParameterSet> pfConfigs = iConfig.getParameter<std::vector<edm::ParameterSet> >("pfCandidates");
  for (size_t i=0; i<pfConfigs.size(); ++i) pfFillers.push_back( new PFFiller(pfConfigs[i], fIsRealData) );
  }

  // Get list of trigger paths to store the triggering object info. of
  //std::vector<std::string> v(iConfig.getParameter<std::vector<std::string> >("hlt_labels"));
  fHLTLabels = iConfig.getParameter<std::vector<std::string> >("hlt_labels");
  if (fHLTLabels.size()>gMaxHltNPaths){
    edm::LogWarning("NTP") << "@SUB=analyze()"
			   << "Maximum number of triggering paths exceeded";
    fHLTLabels.resize(gMaxHltNPaths);
  }
  fTNpaths = fHLTLabels.size();
  
  //OOT pu reweighting
  if( !fIsRealData ) {
    fPileUpData = iConfig.getParameter<std::vector<std::string> >("pu_data");
    fPileUpMC   = iConfig.getParameter<std::vector<std::string> >("pu_mc");
    if(!fPileUpData[0].empty() && !fPileUpMC[0].empty() ){
      LumiWeights_      = edm::LumiReWeighting(fPileUpMC[0], fPileUpData[0], fPileUpMC[1], fPileUpData[1]);
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
  for ( std::vector<PFFiller*>::iterator it = pfFillers.begin(); 
        it != pfFillers.end(); ++it ) {
    list = (*it)->declareProducts();
    for ( PPI ip = list.begin();  ip != list.end(); ++ip )
      produces<edm::InEvent>( ip->first, ip->second );
  }

  isolator.initializePhotonIsolation(kTRUE);
  isolator. setConeSize(0.3);


  if (!doPhotonStuff) {
  // instantiate and initialize the electron ID MVA classes:

  std::vector<std::string> myManualCatWeigths;
  myManualCatWeigths.push_back(edm::FileInPath("DiLeptonAnalysis/NTupleProducer/data/eleIDMVA_weightFiles/Electrons_BDTG_NonTrigV0_Cat1.weights.xml").fullPath());
  myManualCatWeigths.push_back(edm::FileInPath("DiLeptonAnalysis/NTupleProducer/data/eleIDMVA_weightFiles/Electrons_BDTG_NonTrigV0_Cat2.weights.xml").fullPath());
  myManualCatWeigths.push_back(edm::FileInPath("DiLeptonAnalysis/NTupleProducer/data/eleIDMVA_weightFiles/Electrons_BDTG_NonTrigV0_Cat3.weights.xml").fullPath());
  myManualCatWeigths.push_back(edm::FileInPath("DiLeptonAnalysis/NTupleProducer/data/eleIDMVA_weightFiles/Electrons_BDTG_NonTrigV0_Cat4.weights.xml").fullPath());
  myManualCatWeigths.push_back(edm::FileInPath("DiLeptonAnalysis/NTupleProducer/data/eleIDMVA_weightFiles/Electrons_BDTG_NonTrigV0_Cat5.weights.xml").fullPath());
  myManualCatWeigths.push_back(edm::FileInPath("DiLeptonAnalysis/NTupleProducer/data/eleIDMVA_weightFiles/Electrons_BDTG_NonTrigV0_Cat6.weights.xml").fullPath());

  Bool_t manualCat = true;
  
  electronIDMVANonTrig_ = new EGammaMvaEleEstimator();
  electronIDMVANonTrig_->initialize("BDT",
			   EGammaMvaEleEstimator::kNonTrig,
			   manualCat, 
			   myManualCatWeigths);
  

  std::vector<std::string> myManualCatWeigthsTrig;
  myManualCatWeigthsTrig.push_back(edm::FileInPath("DiLeptonAnalysis/NTupleProducer/data/eleIDMVA_weightFiles/Electrons_BDTG_TrigV0_Cat1.weights.xml").fullPath());
  myManualCatWeigthsTrig.push_back(edm::FileInPath("DiLeptonAnalysis/NTupleProducer/data/eleIDMVA_weightFiles/Electrons_BDTG_TrigV0_Cat2.weights.xml").fullPath());
  myManualCatWeigthsTrig.push_back(edm::FileInPath("DiLeptonAnalysis/NTupleProducer/data/eleIDMVA_weightFiles/Electrons_BDTG_TrigV0_Cat3.weights.xml").fullPath());
  myManualCatWeigthsTrig.push_back(edm::FileInPath("DiLeptonAnalysis/NTupleProducer/data/eleIDMVA_weightFiles/Electrons_BDTG_TrigV0_Cat4.weights.xml").fullPath());
  myManualCatWeigthsTrig.push_back(edm::FileInPath("DiLeptonAnalysis/NTupleProducer/data/eleIDMVA_weightFiles/Electrons_BDTG_TrigV0_Cat5.weights.xml").fullPath());
  myManualCatWeigthsTrig.push_back(edm::FileInPath("DiLeptonAnalysis/NTupleProducer/data/eleIDMVA_weightFiles/Electrons_BDTG_TrigV0_Cat6.weights.xml").fullPath());

  electronIDMVATrig_ = new EGammaMvaEleEstimator();
  electronIDMVATrig_->initialize("BDT",
			EGammaMvaEleEstimator::kTrig,
			manualCat,
			myManualCatWeigthsTrig); 


  // initialize the muon isolation MVA
  fMuonIsoMVA = new MuonMVAEstimator();
  vector<string> muoniso_weightfiles;
  muoniso_weightfiles.push_back(edm::FileInPath("Muon/MuonAnalysisTools/data/MuonIsoMVA_sixie-BarrelPt5To10_V0_BDTG.weights.xml").fullPath());
  muoniso_weightfiles.push_back(edm::FileInPath("Muon/MuonAnalysisTools/data/MuonIsoMVA_sixie-EndcapPt5To10_V0_BDTG.weights.xml").fullPath());
  muoniso_weightfiles.push_back(edm::FileInPath("Muon/MuonAnalysisTools/data/MuonIsoMVA_sixie-BarrelPt10ToInf_V0_BDTG.weights.xml").fullPath());
  muoniso_weightfiles.push_back(edm::FileInPath("Muon/MuonAnalysisTools/data/MuonIsoMVA_sixie-EndcapPt10ToInf_V0_BDTG.weights.xml").fullPath());
  muoniso_weightfiles.push_back(edm::FileInPath("Muon/MuonAnalysisTools/data/MuonIsoMVA_sixie-Tracker_V0_BDTG.weights.xml").fullPath());
  muoniso_weightfiles.push_back(edm::FileInPath("Muon/MuonAnalysisTools/data/MuonIsoMVA_sixie-Global_V0_BDTG.weights.xml").fullPath());
  fMuonIsoMVA->initialize("MuonIso_BDTG_IsoRings",
                          MuonMVAEstimator::kIsoRings,
                          true,
                          muoniso_weightfiles);
  fMuonIsoMVA->SetPrintMVADebug(false);
  }


  if (doPhotonStuff) {  // initialize diphoton vertex MVA

    SetupVtxAlgoParams2012(vtxAlgoParams);

    vAna  = new HggVertexAnalyzer(vtxAlgoParams);
    vConv = new HggVertexFromConversions(vtxAlgoParams);
    rankVariables.push_back("ptbal"); rankVariables.push_back("ptasym"); rankVariables.push_back("logsumpt2");

    perVtxVariables = rankVariables;
    if( addConversionToMva ) {
      perVtxVariables.push_back("limPullToConv");
      perVtxVariables.push_back("nConv");
    }
    
    perVtxReader = new TMVA::Reader( "!Color:!Silent" );
    HggVertexAnalyzer::bookVariables( *perVtxReader, perVtxVariables );
    perVtxReader->BookMVA( perVtxMvaMethod, perVtxMvaWeights );
    
    perEvtReader = new TMVA::Reader( "!Color:!Silent" );
    HggVertexAnalyzer::bookPerEventVariables( *perEvtReader );
    perEvtReader->BookMVA( perEvtMvaMethod, perEvtMvaWeights );
    
  }

  if (doPhotonStuff) {
  for (uint i=0; i<jetMVAAlgos.size(); i++) PileupJetIdAlgos.push_back(new PileupJetIdAlgo(jetMVAAlgos.at(i)));
  if (PileupJetIdAlgos.size()>gMaxNPileupJetIDAlgos) PileupJetIdAlgos.resize(gMaxNPileupJetIDAlgos);
  }

  if (doPhotonStuff){
  photonIDMVA_reader_EB = new TMVA::Reader("!Color:Silent");
  photonIDMVA_reader_EB->AddVariable("ph.scrawe",   &(photonIDMVA_variables.scrawe) );
  photonIDMVA_reader_EB->AddVariable("ph.r9",   &(photonIDMVA_variables.r9) );
  photonIDMVA_reader_EB->AddVariable("ph.sigietaieta",   &(photonIDMVA_variables.sieie) );
  photonIDMVA_reader_EB->AddVariable("ph.scetawidth",   &(photonIDMVA_variables.etawidth) );
  photonIDMVA_reader_EB->AddVariable("ph.scphiwidth",   &(photonIDMVA_variables.phiwidth) );
  photonIDMVA_reader_EB->AddVariable("ph.idmva_CoviEtaiPhi",   &(photonIDMVA_variables.sieip) );
  photonIDMVA_reader_EB->AddVariable("ph.idmva_s4ratio",   &(photonIDMVA_variables.s4ratio) );
  photonIDMVA_reader_EB->AddVariable("ph.idmva_GammaIso",   &(photonIDMVA_variables.pfphotoniso03) );
  photonIDMVA_reader_EB->AddVariable("ph.idmva_ChargedIso_selvtx",   &(photonIDMVA_variables.pfchargedisogood03) );
  photonIDMVA_reader_EB->AddVariable("ph.idmva_ChargedIso_worstvtx",   &(photonIDMVA_variables.pfchargedisobad03) );
  photonIDMVA_reader_EB->AddVariable("ph.sceta",   &(photonIDMVA_variables.sceta) );
  photonIDMVA_reader_EB->AddVariable("rho",   &(photonIDMVA_variables.eventrho) );

  photonIDMVA_reader_EE = new TMVA::Reader("!Color:Silent");
  photonIDMVA_reader_EE->AddVariable("ph.scrawe",   &(photonIDMVA_variables.scrawe) );
  photonIDMVA_reader_EE->AddVariable("ph.r9",   &(photonIDMVA_variables.r9) );
  photonIDMVA_reader_EE->AddVariable("ph.sigietaieta",   &(photonIDMVA_variables.sieie) );
  photonIDMVA_reader_EE->AddVariable("ph.scetawidth",   &(photonIDMVA_variables.etawidth) );
  photonIDMVA_reader_EE->AddVariable("ph.scphiwidth",   &(photonIDMVA_variables.phiwidth) );
  photonIDMVA_reader_EE->AddVariable("ph.idmva_CoviEtaiPhi",   &(photonIDMVA_variables.sieip) );
  photonIDMVA_reader_EE->AddVariable("ph.idmva_s4ratio",   &(photonIDMVA_variables.s4ratio) );
  photonIDMVA_reader_EE->AddVariable("ph.idmva_GammaIso",   &(photonIDMVA_variables.pfphotoniso03) );
  photonIDMVA_reader_EE->AddVariable("ph.idmva_ChargedIso_selvtx",   &(photonIDMVA_variables.pfchargedisogood03) );
  photonIDMVA_reader_EE->AddVariable("ph.idmva_ChargedIso_worstvtx",   &(photonIDMVA_variables.pfchargedisobad03) );
  photonIDMVA_reader_EE->AddVariable("ph.sceta",   &(photonIDMVA_variables.sceta) );
  photonIDMVA_reader_EE->AddVariable("rho",   &(photonIDMVA_variables.eventrho) );
  photonIDMVA_reader_EE->AddVariable("ph.idmva_PsEffWidthSigmaRR",   &(photonIDMVA_variables.ESEffSigmaRR) );

  {
    WeightsPhotonIDMVA_EB = iConfig.getParameter<std::string>("tag_WeightsPhotonIDMVA_EB");
    WeightsPhotonIDMVA_EE = iConfig.getParameter<std::string>("tag_WeightsPhotonIDMVA_EE");
    TString descr = getenv("CMSSW_BASE");
    photonIDMVA_reader_EB->BookMVA("AdaBoost",Form("%s/src/DiLeptonAnalysis/NTupleProducer/data/%s",descr.Data(),WeightsPhotonIDMVA_EB.c_str()));
    photonIDMVA_reader_EE->BookMVA("AdaBoost",Form("%s/src/DiLeptonAnalysis/NTupleProducer/data/%s",descr.Data(),WeightsPhotonIDMVA_EE.c_str()));
  }
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
  for ( std::vector<PFFiller*>::iterator it = pfFillers.begin(); 
        it != pfFillers.end(); ++it ) 
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
	
  // sigma for L1FastJet
  edm::Handle<double> sigma;
  iEvent.getByLabel(fSrcSigma,sigma);
  *fTSigma = *sigma;

  // rho for correcting isolation
  edm::Handle<double> rhoForIso;
  iEvent.getByLabel(fSrcRhoForIso,rhoForIso);
  *fTRhoForIso = *rhoForIso;

  // rho for QG tagger systematics 
  edm::Handle<double> rhoForQG;

  if (doPhotonStuff)  {
    iEvent.getByLabel("kt6PFJetsForQGSyst","rho",rhoForQG);
    TString descr = getenv("CMSSW_BASE");
    std::string systDB_fullPath;
    if (QGSystString=="pythia") systDB_fullPath = Form("%s/src/QuarkGluonTagger/EightTeV/data/SystDatabase.txt",descr.Data()); // for Pythia
    else if (QGSystString=="herwig++") systDB_fullPath = Form("%s/src/QuarkGluonTagger/EightTeV/data/SystDatabase_Hpp.txt",descr.Data()); // for Herwig++
    qgsyst.ReadDatabaseDoubleMin(systDB_fullPath);
  }

  edm::Handle<reco::GenParticleCollection> GlobalGenParticles;
  if (!fIsRealData) iEvent.getByLabel(fGenPartTag, GlobalGenParticles);

  // beam halo
  if(!fIsFastSim){
  edm::Handle<BeamHaloSummary> TheBeamHaloSummary;
  iEvent.getByLabel("BeamHaloSummary",TheBeamHaloSummary);
  const BeamHaloSummary TheSummary = (*TheBeamHaloSummary.product());
  *fTCSCTightHaloID = (TheSummary.CSCTightHaloId()) ? 0:1;
  }
  // collect information for b-tagging (4 tags)
  Handle<JetTagCollection> jetsBtag[gMaxNBtags];
  size_t ibtag = 0;
  for ( std::vector<edm::InputTag>::const_iterator it=fBtagTags.begin(); it!=fBtagTags.end(); ++it ) 
    iEvent.getByLabel((*it),jetsBtag[ibtag++]);

  FlavourMap flavours;
  if(!fIsRealData){
  // Get matching parton flavour for jets
  edm::Handle<reco::JetFlavourMatchingCollection> jetMC;
  iEvent.getByLabel(fPartonMatch, jetMC);
  for (reco::JetFlavourMatchingCollection::const_iterator iter = jetMC->begin();
       iter != jetMC->end(); iter++) {
    int fl = iter->second.getFlavour();
    //std::cout << "flavour " << fl << " ";
    flavours.insert(FlavourMap::value_type(iter->first, fl));
    }
  }

  //Get Tracks collection
  Handle<TrackCollection> tracks;
  iEvent.getByLabel(fTrackTag, tracks);

  //Get Photon collection
  Handle<View<Photon> > photons;
  iEvent.getByLabel(fPhotonTag, photons);

  //Get SC collections
  Handle<SuperClusterCollection> BarrelSuperClusters;
  Handle<SuperClusterCollection> EndcapSuperClusters;
  Handle<edm::View<reco::Candidate> > GoodSuperClusters;

  iEvent.getByLabel(fSCTagBarrel,BarrelSuperClusters);
  iEvent.getByLabel(fSCTagEndcap,EndcapSuperClusters);
  iEvent.getByLabel("goodSuperClustersClean", GoodSuperClusters);

	
  // PFcandidates
  edm::Handle<reco::PFCandidateCollection> pfCandidates;
  iEvent.getByLabel(pfProducerTag, pfCandidates);
  const  PFCandidateCollection thePfColl = *(pfCandidates.product());

  // PF candidate isolation
  Handle< edm::ValueMap<float> > muonPfIsoTagsCustom[gMaxNPfIsoTags];
  size_t ipfisotag = 0;
  for ( std::vector<edm::InputTag>::const_iterator it=fMuonPfIsoTagsCustom.begin(); 
        it!=fMuonPfIsoTagsCustom.end(); ++it ) 
    iEvent.getByLabel((*it),muonPfIsoTagsCustom[ipfisotag++]);
  Handle< edm::ValueMap<float> > elePfIsoTagsCustom[gMaxNPfIsoTags];
  ipfisotag = 0;
  for ( std::vector<edm::InputTag>::const_iterator it=fElePfIsoTagsCustom.begin(); 
        it!=fElePfIsoTagsCustom.end(); ++it ) 
    iEvent.getByLabel((*it),elePfIsoTagsCustom[ipfisotag++]);
  Handle< edm::ValueMap<double> > elePfIsoTagsEvent[gMaxNPfIsoTags];
  ipfisotag = 0;
  for ( std::vector<edm::InputTag>::const_iterator it=fElePfIsoTagsEvent.begin(); 
        it!=fElePfIsoTagsEvent.end(); ++it ) 
    iEvent.getByLabel((*it),elePfIsoTagsEvent[ipfisotag++]);
	
  //Electron collection
  edm::Handle<reco::GsfElectronCollection> electronHandle;
  iEvent.getByLabel(fElectronTag, electronHandle);
	
  //PF Photon collection
  edm::Handle<reco::PhotonCollection> pfPhotonHandle;
  iEvent.getByLabel(pfphotonsProducerTag,pfPhotonHandle);
	
  // MET
  Handle<CaloMETCollection> calomet;
  iEvent.getByLabel(fRawCaloMETTag, calomet);

  Handle<METCollection> tcmet;
  iEvent.getByLabel(fTCMETTag, tcmet);

  Handle<View<PFMET> > pfmet;
  iEvent.getByLabel(fPFMETTag, pfmet);
	
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
  const reco::Vertex *primVtx = (vertices->size()>0) ? &(*(vertices.product()))[0] : NULL; // Just take first vertex ...

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

  EcalClusterLazyTools lazyTools( iEvent, iSetup, edm::InputTag("reducedEcalRecHitsEB"), edm::InputTag("reducedEcalRecHitsEE") );

  if (doPhotonStuff) {
  if (!corSemiParm.IsInitialized() && (regrVersion==5 || regrVersion==8)) {
    char filename[500];
    char* descr = getenv("CMSSW_BASE");
    std::string energyRegFilename;
    if (regrVersion==8) {
      energyRegFilename="regweights_v8_7TeV_forest_ph.root";
      sprintf(filename, "%s/src/HiggsAnalysis/GBRLikelihoodEGTools/data/%s", descr, energyRegFilename.c_str());
      corSemiParm.Initialize(filename,8);
    }
    if (regrVersion==5) {
      energyRegFilename="regweights_v5_forest_ph.root";
      sprintf(filename, "%s/src/HiggsAnalysis/GBRLikelihoodEGTools/data/%s", descr, energyRegFilename.c_str());
      corSemiParm.Initialize(filename,5);
    }
  }
  }

  Handle<double> hRhoRegr;
  iEvent.getByLabel(edm::InputTag("kt6PFJets","rho"), hRhoRegr); 

  edm::Handle<edm::ValueMap<float> >  QGTagsHandleMLP;
  edm::Handle<edm::ValueMap<float> >  QGTagsHandleLikelihood;
  if (doPhotonStuff){
    iEvent.getByLabel("QGTagger","qgMLP", QGTagsHandleMLP);
    iEvent.getByLabel("QGTagger","qgLikelihood", QGTagsHandleLikelihood);
  }

  // type-I corrected MET for 2012 analyses
  edm::Handle<View<PFMET> > typeICorMET;
  iEvent.getByLabel("pfType1CorrectedMet",typeICorMET);
  // *fTtypeICorMET = (float) *typeICorMET;
  *fTPFType1MET             = (typeICorMET->front()).pt();
  *fTPFType1METpx           = (typeICorMET->front()).px();
  *fTPFType1METpy           = (typeICorMET->front()).py();
  *fTPFType1METphi          = (typeICorMET->front()).phi();
  double sigmaX2= (typeICorMET->front() ).getSignificanceMatrix()(0,0);
  double sigmaY2= (typeICorMET->front() ).getSignificanceMatrix()(1,1);
  double significance = 0;
  if(sigmaX2<1.e10 && sigmaY2<1.e10) significance = (typeICorMET->front() ).significance();
  *fTPFType1METSignificance = significance;
  *fTPFType1SumEt           = (typeICorMET->front()).sumEt();


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

  if(!fIsRealData && !fIsFastSim){
    // Get LHEEventProduct with partonic momenta. 	
    Handle<LHEEventProduct> evt;
    bool LHEEventProduct_found= iEvent.getByType( evt );
    if(LHEEventProduct_found){ 
      const lhef::HEPEUP hepeup_ = evt->hepeup();
      const std::vector<lhef::HEPEUP::FiveVector> pup_ = hepeup_.PUP; // px, py, pz, E, M
//      cout << "------------" << endl;
//      cout << pup_.size() << " NUP " << hepeup_.NUP << " IDPRUP " << hepeup_.IDPRUP <<  endl;
      float partonicHT=0;
      for(unsigned int i=0; i<pup_.size(); ++i){
//	cout << " ID " << hepeup_.IDUP[i] << " status " << hepeup_.ISTUP[i] << " mother " << hepeup_.MOTHUP[i].first << " " << hepeup_.MOTHUP[i].second << endl;
	fTLHEEventID                  ->push_back(hepeup_.IDUP[i]);
	fTLHEEventStatus              ->push_back(hepeup_.ISTUP[i]);
	fTLHEEventMotherFirst         ->push_back(hepeup_.MOTHUP[i].first);
	fTLHEEventMotherSecond        ->push_back(hepeup_.MOTHUP[i].second);
	fTLHEEventPx                  ->push_back(hepeup_.PUP[i][0]);
	fTLHEEventPy                  ->push_back(hepeup_.PUP[i][1]);
	fTLHEEventPz                  ->push_back(hepeup_.PUP[i][2]);
	fTLHEEventE                   ->push_back(hepeup_.PUP[i][3]);
	fTLHEEventM                   ->push_back(hepeup_.PUP[i][4]);
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
    if(!fPileUpData[0].empty() && !fPileUpMC[0].empty() ){
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
		
//      float Q = pdfstuff->pdf()->scalePDF;
//      int id1 = pdfstuff->pdf()->id.first;
//      double x1 = pdfstuff->pdf()->x.first;
      //       double pdf1 = pdfstuff->pdf()->xPDF.first;
//      int id2 = pdfstuff->pdf()->id.second;
//      double x2 = pdfstuff->pdf()->x.second;
      //       double pdf2 = pdfstuff->pdf()->xPDF.second;

      fTpdfW->push_back(1);
      //LHAPDF::initPDF(0);
      /*
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
      */
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
  const std::vector<std::string> allTrigNames = fHltConfig.triggerNames();
  for (size_t i=0; i<fTNpaths; ++i) {
    size_t j = 0;
    for(; j <allTrigNames.size();++j) {
      if(allTrigNames[j].find(fHLTLabels[i])!=string::npos) break;
    }
    if (j==allTrigNames.size()) continue;
    const std::vector<std::string> filtertags = fHltConfig.moduleLabels(allTrigNames[j]);
    if(filtertags.size() < 2) continue;
    collectionTag = edm::InputTag(filtertags[filtertags.size()-2],"",fHltConfig.processName());
    size_t  filterIndex_ = trgEvent->filterIndex(collectionTag);
    if (filterIndex_<trgEvent->sizeFilters()) {
      const trigger::TriggerObjectCollection& TOC(trgEvent->getObjects());
      const trigger::Keys& keys = trgEvent->filterKeys(filterIndex_);
      // Loop over objects
      int hlto;
      for (hlto = 0; hlto< (int) keys.size(); ++hlto ) {
        const trigger::TriggerObject& TO(TOC[keys[hlto]]);
        fTHLTObjectID[i]->push_back( TO.id() );
        fTHLTObjectPt[i]->push_back( TO.pt() );
        fTHLTObjectEta[i]->push_back( TO.eta() );
        fTHLTObjectPhi[i]->push_back( TO.phi() );
      }
      fTNHLTObjs->push_back(hlto);
    }
  }


  ////////////////////////////////////////////////////////////////////////////////
  // Dump tree variables /////////////////////////////////////////////////////////
  *fTRun   = iEvent.id().run();
  if (iEvent.id().event()>std::numeric_limits<unsigned int>::max()) {cout << "WARNING: DATA FORMAT LONG NOT ENOUGH TO CONTAIN EVENT NUMBER" << endl;}
  *fTEvent = iEvent.id().event();
  *fTLumiSection = iEvent.luminosityBlock();

  *fTWeight = 1.0; // To be filled at some point?

  if (primVtx){
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
  }

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
  *fTNVrtx = fTVrtxX->size();


  // Save position of beamspot
  *fTBeamspotx = (beamSpot.position()).x();
  *fTBeamspoty = (beamSpot.position()).y();
  *fTBeamspotz = (beamSpot.position()).z();

  IndexByPt indexComparator; // Need this to sort collections

  /////////////////////////////////////////
  /// GenVertices 
  if (!fIsRealData && doPhotonStuff){

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
        if (!(part->pt()>0.)) continue;
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
  // FIXME: TO BE REMOVED ONCE WE ARE HAPPY WITH THE FULL GEN. INFO
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
        if (abs(m_id)!=2212) edm::LogWarning("NTP") << "@SUB=analyze"
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
  // FIXME: TO BE REMOVED ONCE WE ARE HAPPY WITH THE FULL GEN. INFO
  if(!fIsRealData){
    edm::Handle<GenParticleCollection> gen;
    iEvent.getByLabel(fGenPartTag, gen);
    GenParticleCollection::const_iterator g_part;
	
    // Steve Mrenna's status 2 parton jets
    edm::Handle<GenJetCollection> partonGenJets;
    iEvent.getByLabel("partonGenJets", partonGenJets);
    GenJetCollection::const_iterator pGenJet;

    edm::Handle<View<Candidate> > partons;
    iEvent.getByLabel("partons", partons);
	 
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
      fTGenPhotonPartonMindR->push_back( -999.99 ); // Initialize

      fTGenPhotonIsoDR03->push_back(GenPartonicIso_allpart(*(gen_photons[i]),gen,0.3));
      fTGenPhotonIsoDR04->push_back(GenPartonicIso_allpart(*(gen_photons[i]),gen,0.4));

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
      (*fTGenPhotonPartonMindR)[i] = minDR;
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
    Ref<View<Muon> > muonRef(muons,index);

    fTMuIsGlobalMuon->push_back( muon.isGlobalMuon() ? 1:0 );
    fTMuIsTrackerMuon->push_back( muon.isTrackerMuon() ? 1:0 );
    fTMuIsPFMuon->push_back( muon.isPFMuon() ? 1:0);
    fTMuIsStandaloneMuon->push_back( muon.isStandAloneMuon() ? 1:0 );

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

    fTMuPfIsoR03ChHad  ->push_back( muon.pfIsolationR03().sumChargedHadronPt );
    fTMuPfIsoR03NeHad  ->push_back( muon.pfIsolationR03().sumNeutralHadronEt );
    fTMuPfIsoR03Photon ->push_back( muon.pfIsolationR03().sumPhotonEt );
    fTMuPfIsoR03NeHadHighThresh  ->push_back( muon.pfIsolationR03().sumNeutralHadronEtHighThreshold );
    fTMuPfIsoR03PhotonHighThresh ->push_back( muon.pfIsolationR03().sumPhotonEtHighThreshold );
    fTMuPfIsoR03SumPUPt->push_back( muon.pfIsolationR03().sumPUPt );

    fTMuPfIsoR04ChHad  ->push_back( muon.pfIsolationR04().sumChargedHadronPt );
    fTMuPfIsoR04NeHad  ->push_back( muon.pfIsolationR04().sumNeutralHadronEt );
    fTMuPfIsoR04Photon ->push_back( muon.pfIsolationR04().sumPhotonEt );
    fTMuPfIsoR04NeHadHighThresh  ->push_back( muon.pfIsolationR04().sumNeutralHadronEtHighThreshold );
    fTMuPfIsoR04PhotonHighThresh ->push_back( muon.pfIsolationR04().sumPhotonEtHighThreshold );
    fTMuPfIsoR04SumPUPt->push_back( muon.pfIsolationR04().sumPUPt );

    // PF isolations:
    ipfisotag = 0;
    for ( std::vector<edm::InputTag>::const_iterator it=fMuonPfIsoTagsCustom.begin(); 
          it!=fMuonPfIsoTagsCustom.end(); ++it ) {
      fTMuPfIsosCustom[ipfisotag]->push_back( (*muonPfIsoTagsCustom[ipfisotag])[muonRef] );
      ++ipfisotag;
    }
 
    // mva iso:
    const reco::GsfElectronCollection dummyIdentifiedEleCollection;
    const reco::MuonCollection dummyIdentifiedMuCollection;
    if (!doPhotonStuff){
    double isomva = fMuonIsoMVA->mvaValue( muon,
                                        vertices->front(),
                                        *pfCandidates,
                                        *rho,
                                        MuonEffectiveArea::kMuEAFall11MC,
                                        dummyIdentifiedEleCollection,
                                        dummyIdentifiedMuCollection);
    fTMuIsoMVA->push_back( isomva );
    }

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

    const reco::IsoDeposit ECDep = ECDepMap[muonRef];
    const reco::IsoDeposit HCDep = HCDepMap[muonRef];
    fTMuEem->push_back( ECDep.candEnergy() );
    fTMuEhad->push_back( HCDep.candEnergy() );

    // 3D impact parameter
    TransientTrack mutt = theB->build( muon.innerTrack() );
    
    if (primVtx){
      Measurement1D muip3dpv = IPTools::absoluteImpactParameter3D(mutt, *(primVtx)).second;      
      fTMuD03DPV->push_back( muip3dpv.value() );
      fTMuD03DE ->push_back( muip3dpv.error() );
    }
    else{
      fTMuD03DPV->push_back( 999 );
      fTMuD03DE ->push_back( 999 );
    }
    
    fTMuD0BS->push_back( -1.0*muon.innerTrack()->dxy(beamSpot.position()) );
    fTMuD0PV->push_back( (primVtx) ? -1.0*muon.innerTrack()->dxy(primVtx->position()) : 999 );
    fTMuDzPV->push_back( (primVtx) ? muon.innerTrack()->dz(primVtx->position()) : 999);
    fTMuInnerTkNChi2->push_back( muon.innerTrack()->normalizedChi2() );
    fTMuNSiLayers->push_back (muon.innerTrack()->hitPattern().trackerLayersWithMeasurement());

    // Separate methods:
    if((*fTMuIsTrackerMuon)[mqi]){ // Tracker Muons
      (*fTNTMus)++;
      fTMuTkPtE->push_back( muon.innerTrack()->ptError() );
      fTMuTkD0E->push_back( muon.innerTrack()->dxyError() );
      fTMuTkDzE->push_back( muon.innerTrack()->dzError() );

      fTMuNTkHits->push_back( muon.innerTrack()->hitPattern().numberOfValidHits() );
      fTMuNPxHits->push_back( muon.innerTrack()->hitPattern().numberOfValidPixelHits() );
    } else {
      fTMuTkPtE->push_back( 0.0 );
      fTMuTkD0E->push_back( 0.0 );
      fTMuTkDzE->push_back( 0.0 );
      fTMuNTkHits->push_back( 0.0 );
      fTMuNPxHits->push_back( 0.0 );
    }
    if((*fTMuIsGlobalMuon)[mqi]){ // Global Muons
      (*fTNGMus)++;
      fTMuPtE->push_back( muon.globalTrack()->ptError() );
      fTMuD0E->push_back( muon.globalTrack()->dxyError() );
      fTMuDzE->push_back( muon.globalTrack()->dzError() );

      fTMuNChi2->push_back( muon.globalTrack()->normalizedChi2() );
      fTMuNGlHits->push_back( muon.globalTrack()->hitPattern().numberOfValidHits() );
      fTMuNGlMuHits->push_back( muon.globalTrack()->hitPattern().numberOfValidMuonHits() );
      fTMuNMuHits->push_back( muon.outerTrack()->hitPattern().numberOfValidHits() );
      fTMuNMatches->push_back( muon.numberOfMatches() );
      fTMuNMatchedStations->push_back( muon.numberOfMatchedStations() );
      fTMuNChambers->push_back( muon.numberOfChambers() );
    } else {
      fTMuPtE->push_back( 0.0 );
      fTMuD0E->push_back( 0.0 );
      fTMuDzE->push_back( 0.0 );
      fTMuNChi2->push_back( 0.0 );
      fTMuNGlHits->push_back( 0.0 );
      fTMuNGlMuHits->push_back( 0.0 );
      fTMuNMuHits->push_back( 0.0 );
      fTMuNMatches->push_back( 0.0 );
      fTMuNMatchedStations->push_back( 0.0 );
      fTMuNChambers->push_back( 0.0 );
    }

    // MC Matching
    // FIXME: TO BE REMOVED ONCE WE ARE HAPPY WITH THE FULL GEN. INFO
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


  (*fTNGoodSuperClusters)=0;
  for (edm::View<reco::Candidate>::const_iterator sc = GoodSuperClusters->begin(); sc!=GoodSuperClusters->end(); ++sc){

    if (*fTNGoodSuperClusters>=gMaxNSC) {
      edm::LogWarning("NTP") << "@SUB=analyze" << "Maximum number of Super Clusters exceeded";
      *fTGoodEvent = 1;
      break;
    }

    fTGoodSCEnergy->push_back(sc->energy());
    fTGoodSCEta->push_back(sc->eta());
    fTGoodSCPhi->push_back(sc->phi());
    (*fTNGoodSuperClusters)++;
  }


  // SC variables
  (*fTNSuperClusters)=0;
  std::vector<DetId> cristalli_tokeep;
  for (SuperClusterCollection::const_iterator sc = BarrelSuperClusters->begin(); sc!=BarrelSuperClusters->end(); ++sc){

    if (sc->rawEnergy()<fMinSCraw) continue;
    if (sc->rawEnergy()/TMath::CosH(sc->eta())<fMinSCrawPt) continue;

    if (*fTNSuperClusters>=gMaxNSC) {
      edm::LogWarning("NTP") << "@SUB=analyze" << "Maximum number of Super Clusters exceeded"; 
      *fTGoodEvent = 1; 
      break;
    }

    fTSCX->push_back(sc->x());
    fTSCY->push_back(sc->y());
    fTSCZ->push_back(sc->z());
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

    if (doPhotonStuff) {
      std::vector<DetId> cristalli;
      for (reco::CaloCluster_iterator bc=sc->clustersBegin(); bc!=sc->clustersEnd(); ++bc){
	const std::vector< std::pair<DetId, float> > & seedrechits = (*bc)->hitsAndFractions();
	for (uint i=0; i<seedrechits.size(); i++) cristalli.push_back(seedrechits[i].first);
	sort(cristalli.begin(),cristalli.end());
	std::vector<DetId>::iterator it;
	it = unique(cristalli.begin(),cristalli.end());
	cristalli.resize(it-cristalli.begin());
      }

      fTSCXtalListStart->push_back(cristalli_tokeep.size());

      int added=0;      
      for (unsigned int i=0; i<cristalli.size(); i++){

	if (cristalli.at(i).subdetId()!=EcalBarrel) {
	  edm::LogWarning("NTP") << "@SUB=analyze" << "Problem with xtals subdetId()";
	  continue;
	}

	cristalli_tokeep.push_back(cristalli.at(i));	
	added++;

      }

      fTSCNXtals->push_back(added);

    }

    (*fTNSuperClusters)++;
  }

  for (SuperClusterCollection::const_iterator sc = EndcapSuperClusters->begin(); sc!=EndcapSuperClusters->end(); ++sc){

    if (sc->rawEnergy()<fMinSCraw) continue;
    if (sc->rawEnergy()/TMath::CosH(sc->eta())<fMinSCrawPt) continue;

    if (*fTNSuperClusters>=gMaxNSC) {
      edm::LogWarning("NTP") << "@SUB=analyze" << "Maximum number of Super Clusters exceeded"; 
      *fTGoodEvent = 1; 
      break;
    }

    fTSCX->push_back(sc->x());
    fTSCY->push_back(sc->y());
    fTSCZ->push_back(sc->z());
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

    if (doPhotonStuff) {
      std::vector<DetId> cristalli;
      for (reco::CaloCluster_iterator bc=sc->clustersBegin(); bc!=sc->clustersEnd(); ++bc){
	const std::vector< std::pair<DetId, float> > & seedrechits = (*bc)->hitsAndFractions();
	for (uint i=0; i<seedrechits.size(); i++) cristalli.push_back(seedrechits[i].first);
	sort(cristalli.begin(),cristalli.end());
	std::vector<DetId>::iterator it;
	it = unique(cristalli.begin(),cristalli.end());
	cristalli.resize(it-cristalli.begin());
      }
      
      fTSCXtalListStart->push_back(cristalli_tokeep.size());

      int added=0;
      for (unsigned int i=0; i<cristalli.size(); i++){

	if (cristalli.at(i).subdetId()!=EcalEndcap) {
	  edm::LogWarning("NTP") << "@SUB=analyze" << "Problem with xtals subdetId()";
	  continue;
	}

	cristalli_tokeep.push_back(cristalli.at(i));	
	added++;

      }

      fTSCNXtals->push_back(added);

    }

    (*fTNSuperClusters)++;
  }


  ////// Xtal position information

  if (doPhotonStuff) {
  (*fTNXtals)=cristalli_tokeep.size();  
  for (unsigned int i=0; i<cristalli_tokeep.size(); i++){
    bool isEB = (cristalli_tokeep.at(i).subdetId()==EcalBarrel);
    CaloCellGeometry *cellGeometry = NULL;
    TVector3 xtal_position;
    TVector3 cell_corners[4];
    float dphi;
    float deta;
    if (isEB){
      EBDetId ebDetId  = cristalli_tokeep.at(i);
      cellGeometry = (CaloCellGeometry*)(barrelGeometry->getGeometry(ebDetId));
      xtal_position = TVector3(cellGeometry->getPosition().x(),cellGeometry->getPosition().y(),cellGeometry->getPosition().z());
      dphi=(dynamic_cast<const EcalBarrelGeometry*>(barrelGeometry))->deltaPhi(ebDetId);
      deta=(dynamic_cast<const EcalBarrelGeometry*>(barrelGeometry))->deltaEta(ebDetId);
    }
    else {
      EEDetId eeDetId  = cristalli_tokeep.at(i);
      cellGeometry = (CaloCellGeometry*)(endcapGeometry->getGeometry(eeDetId));
      xtal_position = TVector3(cellGeometry->getPosition().x(),cellGeometry->getPosition().y(),cellGeometry->getPosition().z());
      dphi=(dynamic_cast<const EcalEndcapGeometry*>(endcapGeometry))->deltaPhi(eeDetId);
      deta=(dynamic_cast<const EcalEndcapGeometry*>(endcapGeometry))->deltaEta(eeDetId);
    }
    const CaloCellGeometry::CornersVec& cellCorners (cellGeometry->getCorners());
    for (int k=0; k<4; k++) cell_corners[k] = TVector3((float)(cellCorners[k].x()),(float)(cellCorners[k].y()),(float)(cellCorners[k].z()));
  
    fTXtalX->push_back(xtal_position.x());
    fTXtalY->push_back(xtal_position.y());
    fTXtalZ->push_back(xtal_position.z());
    fTXtalPhiWidth->push_back(dphi);
    fTXtalEtaWidth->push_back(deta);
    fTXtalFront1X->push_back(cell_corners[0].x());
    fTXtalFront1Y->push_back(cell_corners[0].y());
    fTXtalFront1Z->push_back(cell_corners[0].z());
    fTXtalFront2X->push_back(cell_corners[1].x());
    fTXtalFront2Y->push_back(cell_corners[1].y());
    fTXtalFront2Z->push_back(cell_corners[1].z());
    fTXtalFront3X->push_back(cell_corners[2].x());
    fTXtalFront3Y->push_back(cell_corners[2].y());
    fTXtalFront3Z->push_back(cell_corners[2].z());
    fTXtalFront4X->push_back(cell_corners[3].x());
    fTXtalFront4Y->push_back(cell_corners[3].y());
    fTXtalFront4Z->push_back(cell_corners[3].z());
 
  }
  }
   

  ////////////////////////////////////////////////////////
  // Electron variables:
  // Keep pointers to electron superCluster in original collections
  std::vector<const SuperCluster*> elecPtr;
  std::vector<const GsfTrack*> trckPtr;
  std::vector<edm::Ptr<GsfElectron> > elPtrVector;
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

    // Dump electron properties in tree variables
    for( std::vector<OrderPair>::const_iterator it = elOrdered.begin();
         it != elOrdered.end(); ++it, ++eqi ) {

      int index = it->first;
      const GsfElectron& electron = (*electrons)[index];
      Ref<View<GsfElectron> > electronRef(electrons,index);

      const EcalRecHitCollection *rechits = 0;

      float seedEta = electron.superCluster()->eta();
      if( fabs(seedEta) < 1.479 ) rechits = ebRecHits.product();
      else rechits = eeRecHits.product(); 

      // Save the electron SuperCluster pointer
      elecPtr.push_back(&(*electron.superCluster()));
      trckPtr.push_back(&(*electron.gsfTrack()));
      elPtrVector.push_back(edm::Ptr<GsfElectron>(electrons,index));

      fTElPx                     ->push_back(electron.px());
      fTElPy                     ->push_back(electron.py());
      fTElPz                     ->push_back(electron.pz());
      fTElPt                     ->push_back(electron.pt());
      fTElPtE                    ->push_back(electron.gsfTrack()->ptError());
      fTElEta                    ->push_back(electron.eta());
      fTElPhi                    ->push_back(electron.phi());
      fTElIsEB                   ->push_back(electron.isEB());
      fTElIsEE                   ->push_back(electron.isEE());
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
      fTElD0PV                   ->push_back((primVtx) ? -1.0*electron.gsfTrack()->dxy(primVtx->position()) : 999);
      fTElD0E                    ->push_back(electron.gsfTrack()->dxyError());
      fTElDzBS                   ->push_back(electron.gsfTrack()->dz(beamSpot.position()));
      fTElDzPV                   ->push_back((primVtx) ? electron.gsfTrack()->dz(primVtx->position()) : 999);
      fTElDzE                    ->push_back(electron.gsfTrack()->dzError());
      fTElNChi2                  ->push_back(electron.gsfTrack()->normalizedChi2());

      // 3D impact parameter
      TransientTrack eltt = theB->build( electron.gsfTrack() );

      if (primVtx){
	Measurement1D elip3dpv = IPTools::absoluteImpactParameter3D(eltt, *(primVtx)).second;
	fTElD03DPV->push_back( elip3dpv.value() );
	fTElD03DE ->push_back( elip3dpv.error() );
      }
      else{
	fTElD03DPV->push_back( 999 );
	fTElD03DE ->push_back( 999 );
      }

      // ctf track info:
      bool validKF= false; 
      reco::TrackRef myTrackRef = electron.closestCtfTrackRef();
      validKF = (myTrackRef.isAvailable() && myTrackRef.isNonnull()); 
      fTElKfTrkchi2              ->push_back( validKF ? myTrackRef->normalizedChi2() : 0 );
      fTElKfTrkhits              ->push_back( validKF ? myTrackRef->hitPattern().trackerLayersWithMeasurement() : -1. );


      // Isolation:
      fTElDR03TkSumPt            ->push_back(electron.dr03TkSumPt());
      fTElDR03EcalRecHitSumEt    ->push_back(electron.dr03EcalRecHitSumEt());
      fTElDR03HcalTowerSumEt     ->push_back(electron.dr03HcalTowerSumEt());
      fTElDR04TkSumPt            ->push_back(electron.dr04TkSumPt());
      fTElDR04EcalRecHitSumEt    ->push_back(electron.dr04EcalRecHitSumEt());
      fTElDR04HcalTowerSumEt     ->push_back(electron.dr04HcalTowerSumEt());
      fTElRelIso03               ->push_back(((*fTElDR03TkSumPt)[eqi] + (*fTElDR03EcalRecHitSumEt)[eqi] + (*fTElDR03HcalTowerSumEt)[eqi]) / (*fTElPt)[eqi]);
      fTElRelIso04               ->push_back(((*fTElDR04TkSumPt)[eqi] + (*fTElDR04EcalRecHitSumEt)[eqi] + (*fTElDR04HcalTowerSumEt)[eqi]) / (*fTElPt)[eqi]);

      // PF isolations:
      fTElPfIsoChHad03           ->push_back(electron.pfIsolationVariables().chargedHadronIso );
      fTElPfIsoNeHad03           ->push_back(electron.pfIsolationVariables().neutralHadronIso );
      fTElPfIsoPhoton03          ->push_back(electron.pfIsolationVariables().photonIso );

      ipfisotag = 0;
      for ( std::vector<edm::InputTag>::const_iterator it=fElePfIsoTagsCustom.begin(); 
            it!=fElePfIsoTagsCustom.end(); ++it ) {
        fTElPfIsosCustom[ipfisotag]->push_back( (*elePfIsoTagsCustom[ipfisotag])[electronRef] );
        ++ipfisotag;
      }
      ipfisotag = 0;
      for ( std::vector<edm::InputTag>::const_iterator it=fElePfIsoTagsEvent.begin(); 
            it!=fElePfIsoTagsEvent.end(); ++it ) {
        fTElPfIsosEvent[ipfisotag]->push_back( (*elePfIsoTagsEvent[ipfisotag])[electronRef] );
        ++ipfisotag;
      }

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
      } else {
        fTElClosestCtfTrackPt      ->push_back(-999.99);
        fTElClosestCtfTrackEta     ->push_back(-999.99);
        fTElClosestCtfTrackPhi     ->push_back(-999.99);
        fTElClosestCtfTrackCharge  ->push_back(-999.99);
      }
      fTElBasicClustersSize         ->push_back(electron.basicClustersSize());
      fTElfbrem                     ->push_back(electron.fbrem());
      fTElHcalOverEcal              ->push_back(electron.hcalOverEcal());
      fTElHcalOverEcalBc            ->push_back(electron.hcalOverEcalBc());
      fTElE1x5                      ->push_back(electron.e1x5());
      fTElE5x5                      ->push_back(electron.e5x5());
      fTElE2x5Max                   ->push_back(electron.e2x5Max());
      fTElR9                        ->push_back(electron.superCluster()->rawEnergy()!=0. ? EcalClusterTools::e3x3(*(electron.superCluster()->seed()),rechits, topology) / electron.superCluster()->rawEnergy() : 0.);
      fTElSigmaIetaIeta             ->push_back(electron.sigmaIetaIeta());
      std::vector<float> vCov = EcalClusterTools::localCovariances(*(electron.superCluster()->seed()), &(*rechits), topology ) ;
      fTElSigmaIphiIphi             ->push_back(!isnan(vCov[2]) ? sqrt(vCov[2]) : 0.);
      fTElScEtaWidth                ->push_back(electron.superCluster()->etaWidth());
      fTElScPhiWidth                ->push_back(electron.superCluster()->phiWidth());
      fTElDeltaEtaSeedClusterAtCalo ->push_back(electron.deltaEtaSeedClusterTrackAtCalo());
      fTElDeltaPhiSeedClusterAtCalo ->push_back(electron.deltaPhiSeedClusterTrackAtCalo());
      fTElDeltaPhiSuperClusterAtVtx ->push_back(electron.deltaPhiSuperClusterTrackAtVtx());
      fTElDeltaEtaSuperClusterAtVtx ->push_back(electron.deltaEtaSuperClusterTrackAtVtx());
      fTElCaloEnergy                ->push_back(electron.caloEnergy());
      fTElTrkMomAtVtx               ->push_back(electron.trackMomentumAtVtx().R());
      fTElESuperClusterOverP        ->push_back(electron.eSuperClusterOverP());
      fTElIoEmIoP                   ->push_back((1.0/electron.ecalEnergy()) - (1.0 / electron.p()));
      fTElEoPout                    ->push_back(electron.eEleClusterOverPout());
      fTElPreShowerOverRaw          ->push_back(electron.superCluster()->preshowerEnergy() / electron.superCluster()->rawEnergy());
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


      const TransientTrackBuilder thebuilder = *(theB.product());

      if (!doPhotonStuff) {
      fTElIDMVATrig          ->push_back( electronIDMVATrig_->mvaValue( electron, vertices->front(), thebuilder, lazyTools, false ) );
      fTElIDMVANoTrig        ->push_back( electronIDMVANonTrig_->mvaValue( electron, vertices->front(), thebuilder, lazyTools, false ) );
      }

      {
        fTElSCindex->push_back( -1 ); // Initialize
        float diff=1e+4;
        for (int scind=0; scind<*fTNSuperClusters; scind++){
	  if (fabs((*fTSCEta)[scind]-electron.superCluster()->eta())>0.1) continue;
	  if (fabs((*fTSCPhi)[scind]-electron.superCluster()->phi())>0.1) continue;
	  if (fabs((*fTSCRaw)[scind]-electron.superCluster()->rawEnergy())<diff && fabs((*fTSCRaw)[scind]/electron.superCluster()->rawEnergy()-1)<0.5) {
            (*fTElSCindex)[eqi]=scind;
            diff=fabs((*fTSCRaw)[scind]-electron.superCluster()->rawEnergy());
          }
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
      edm::Handle<reco::BeamSpot> beamspotHandle;
      iEvent.getByLabel("offlineBeamSpot", beamspotHandle);
      const reco::BeamSpot &beamspot = *beamspotHandle.product();
      
      edm::Handle<reco::ConversionCollection> hConversions;
      iEvent.getByLabel("allConversions", hConversions);
      
      bool passconversionveto = !ConversionTools::hasMatchedConversion(electron,hConversions, beamspot.position());
      fTElPassConversionVeto->push_back(passconversionveto);
        
      reco::GsfElectron::ConversionRejection ConvRejVars = electron.conversionRejectionVariables();
      reco::TrackBaseRef ConvPartnerTrack = ConvRejVars.partner;
      if( ConvPartnerTrack.isNonnull() ){
        fTElConvPartnerTrkDist   ->push_back(ConvRejVars.dist);
        fTElConvPartnerTrkDCot   ->push_back(ConvRejVars.dcot);
        fTElConvPartnerTrkPt     ->push_back(ConvPartnerTrack->pt());
        fTElConvPartnerTrkEta    ->push_back(ConvPartnerTrack->eta());
        fTElConvPartnerTrkPhi    ->push_back(ConvPartnerTrack->phi());
        fTElConvPartnerTrkCharge ->push_back(ConvPartnerTrack->charge());
      } else {
        fTElConvPartnerTrkDist   ->push_back(-999.99);
        fTElConvPartnerTrkDCot   ->push_back(-999.99);
        fTElConvPartnerTrkPt     ->push_back(-999.99);
        fTElConvPartnerTrkEta    ->push_back(-999.99);
        fTElConvPartnerTrkPhi    ->push_back(-999.99);
        fTElConvPartnerTrkCharge ->push_back(-999.99);
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
      if((*fTNEBhits)>=gMaxNEBhits)
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

  std::vector<bool> storethispfcand(pfCandidates->size(),false);
  std::vector<int> PhotonToPFPhotonOrElectronMatchingArray(gMaxNPhotons,-999);
  std::vector<int> PhotonToPFPhotonOrElectronMatchingArrayTranslator(gMaxNPhotons,-999);
  std::vector<std::vector<int> > list_pfcand_footprint(gMaxNPhotons,std::vector<int>());
  std::vector<std::vector<int> > list_pfcand_footprintTranslator(gMaxNPhotons,std::vector<int>());

  for (std::vector<OrderPair>::const_iterator it = phoOrdered.begin();
       it != phoOrdered.end(); ++it, ++phoqi ) {

    int index = it->first;
    const Photon& photon = (*photons)[index];

    fTPhoVrtxListStart->push_back((fTPhoVrtxListStart->size()==0) ? 0 : fTPhoVrtxListStart->back()+(*fTNVrtx));

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
    fTPhoHoverE2012      ->push_back(photon.hadTowOverEm());
    fTPhoSigmaIetaIeta   ->push_back(photon.sigmaIetaIeta());
    fTPhoSigmaEtaEta     ->push_back(photon.sigmaEtaEta());
    fTPhoSigmaRR         ->push_back(lazyTools.eseffsirir(*photon.superCluster()));
    fTPhoVx             ->push_back(photon.vx());
    fTPhoVy             ->push_back(photon.vy());
    fTPhoVz             ->push_back(photon.vz());

    if (doPhotonStuff) {
      double ecor, sigeovere, mean, sigma, alpha1, n1, alpha2, n2, pdfval;
      ecor=-999;
      sigeovere=-999;
      sigma=-999;
      if (regrVersion==8) corSemiParm.CorrectedEnergyWithErrorV8(photon, *(vertices.product()), *hRhoRegr, lazyTools, iSetup, ecor, sigeovere, mean, sigma, alpha1, n1, alpha2, n2, pdfval);
      else if (regrVersion==5) corSemiParm.CorrectedEnergyWithErrorV5(photon, *(vertices.product()), *hRhoRegr, lazyTools, iSetup, ecor, sigma, alpha1, n1, alpha2, n2, pdfval);

      fTPhoRegrEnergy     ->push_back(ecor);
      fTPhoRegrEnergyErr  ->push_back((regrVersion==8) ? sigeovere : sigma);
    }


  /*  
  the following lines are from:
  https://twiki.cern.ch/twiki/bin/view/CMS/HoverE2012
*/
    edm::Handle<reco::BeamSpot> beamspotHandle;
    iEvent.getByLabel("offlineBeamSpot", beamspotHandle);
    const reco::BeamSpot &beamspot = *beamspotHandle.product();

    edm::Handle<reco::ConversionCollection> hConversions;
    iEvent.getByLabel("allConversions", hConversions);

    bool passed_PhotonVeto = !ConversionTools::hasMatchedPromptElectron(photon.superCluster(), electronHandle, hConversions, beamspot.position());
    fTPhoPassConversionVeto->push_back(passed_PhotonVeto);

    if (vertices->size()>0){
      isolator.fGetIsolation(&photon,&thePfColl, VertexRef(vertices,0), vertices);
      fTPhoNewIsoPFCharged->push_back(isolator.getIsolationCharged());
      fTPhoNewIsoPFPhoton->push_back(isolator.getIsolationPhoton());
      fTPhoNewIsoPFNeutral->push_back(isolator.getIsolationNeutral());
    }
    else {
      fTPhoNewIsoPFCharged->push_back(999);
      fTPhoNewIsoPFPhoton->push_back(999);
      fTPhoNewIsoPFNeutral->push_back(999);
    }

  
    float PhoHCalIso2012ConeDR03 = photon.hcalTowerSumEtConeDR03() + (photon.hadronicOverEm() - photon.hadTowOverEm())*photon.superCluster()->energy()/cosh(photon.superCluster()->eta());
    fTPhoHCalIso2012ConeDR03->push_back(PhoHCalIso2012ConeDR03);

    {
      edm::ParameterSet myiConfig = edm::ParameterSet();
      myiConfig.insert(true,"tag_jets",edm::Entry("tag_jets",edm::InputTag("ak5PFJetsCorrected"),false));
      SuperClusterFootprintRemoval remover(iEvent,iSetup,myiConfig);
      PFIsolation_struct isos = remover.PFIsolation(photon.superCluster(),(vertices->size()>0) ? edm::Ptr<Vertex>(vertices,0) : edm::Ptr<Vertex>());
      fTPhoSCRemovalPFIsoCharged->push_back(isos.chargediso);
      fTPhoSCRemovalPFIsoChargedPrimVtx->push_back(isos.chargediso_primvtx);
      fTPhoSCRemovalPFIsoNeutral->push_back(isos.neutraliso);
      fTPhoSCRemovalPFIsoPhoton->push_back(isos.photoniso);
      fTPhoSCRemovalPFIsoChargedRCone->push_back(isos.chargediso_rcone);
      fTPhoSCRemovalPFIsoChargedPrimVtxRCone->push_back(isos.chargediso_primvtx_rcone);
      fTPhoSCRemovalPFIsoNeutralRCone->push_back(isos.neutraliso_rcone);
      fTPhoSCRemovalPFIsoPhotonRCone->push_back(isos.photoniso_rcone);
      fTPhoSCRemovalRConeEta->push_back(isos.eta_rcone);
      fTPhoSCRemovalRConePhi->push_back(isos.phi_rcone);
      for (size_t i=0; i<isos.pfcandindex_footprint.size(); i++) list_pfcand_footprint.at(phoqi).push_back(isos.pfcandindex_footprint.at(i));
    }

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
    fTPhoSCX            ->push_back(photon.superCluster()->x());
    fTPhoSCY            ->push_back(photon.superCluster()->y());
    fTPhoSCZ            ->push_back(photon.superCluster()->z());
    fTPhoSCEta          ->push_back(photon.superCluster()->eta());
    fTPhoSCEtaWidth     ->push_back(photon.superCluster()->etaWidth());
    fTPhoSCPhiWidth     ->push_back(photon.superCluster()->phiWidth());
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

    {
      std::vector<float> viCov;
      const reco::CaloClusterPtr  seed_clu = photon.superCluster()->seed();
      viCov = lazyTools.localCovariances(*seed_clu);
      fTPhoSigmaIetaIphi   ->push_back(viCov[1]);
      fTPhoSigmaIphiIphi   ->push_back(viCov[2]);
      float lambdaMinus         = (viCov[0] + viCov[2] - sqrt(pow(viCov[0] - viCov[2], 2) + 4*pow(viCov[1], 2)));
      float lambdaPlus          = (viCov[0] + viCov[2] + sqrt(pow(viCov[0] - viCov[2], 2) + 4*pow(viCov[1], 2)));
      fTPhoLambdaRatio     ->push_back(lambdaMinus/lambdaPlus);

      float e2x2 =  lazyTools.e2x2(*seed_clu);
      float bc_s25 = EcalClusterTools::e5x5(*(photon.superCluster()->seed()), (fTPhoisEB->at(phoqi)) ? &(*ebRecHits) : &(*eeRecHits), &(*topology));
      fTPhoS4Ratio         ->push_back(e2x2/bc_s25);
    }


    fTPhoConvValidVtx->push_back(false);
    fTPhoConvChi2Probability->push_back(-999.);
    fTPhoConvNtracks->push_back(-999.);
    fTPhoConvEoverP->push_back(-999.);

    if (doPhotonStuff && photon.hasConversionTracks()) { // photon conversions

      reco::ConversionRefVector conversions = photon.conversions();
      if (conversions.size()<1) { std::cout << "something wrong here" << std::endl; }
      reco::ConversionRef conv = conversions[0];
      (*fTPhoConvValidVtx)[phoqi]=conv->conversionVertex().isValid();

      if ((*fTPhoConvValidVtx)[phoqi] && ConversionsCut(*conv)) {
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
        fTPhoMCmatchexitcode->push_back(-999);
        fTPhoMCmatchindex->push_back(-999);
        for(int i=0; i<*fTNGenPhotons; ++i){
          if ( (fabs((*fTGenPhotonPt)[i]-matched[0]->pt())<0.01*matched[0]->pt()) 
               && (fabs((*fTGenPhotonEta)[i]-matched[0]->eta())<0.01) 
               && ( fabs(reco::deltaPhi((*fTGenPhotonPhi)[i],matched[0]->phi()))<0.01 ) ) {
            (*fTPhoMCmatchindex)[phoqi] = i;
          }
        }

        if ((*fTPhoMCmatchindex)[phoqi] != -999) {
          if ( (*fTGenPhotonMotherID)[(*fTPhoMCmatchindex)[phoqi]]>=-6 && (*fTGenPhotonMotherID)[(*fTPhoMCmatchindex)[phoqi]]<=6) (*fTPhoMCmatchexitcode)[phoqi]=1;
          else if ((*fTGenPhotonMotherID)[(*fTPhoMCmatchindex)[phoqi]]==21 ) (*fTPhoMCmatchexitcode)[phoqi]=1;
          else if ((*fTGenPhotonMotherID)[(*fTPhoMCmatchindex)[phoqi]]==22 && (*fTGenPhotonMotherStatus)[(*fTPhoMCmatchindex)[phoqi]]==3) (*fTPhoMCmatchexitcode)[phoqi]=2;
          else (*fTPhoMCmatchexitcode)[phoqi] = 3;
        }
        else (*fTPhoMCmatchexitcode)[phoqi] = -2;
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
            fabs(reco::deltaPhi((*fTSCPhi)[(*fTPhotSCindex)[phoqi]],photon.superCluster()->phi()))>0.1){
          (*fTPhotSCindex)[phoqi] = -1;	    
        }

//      if ( (*fTPhotSCindex)[phoqi]==-1) {
//        edm::LogWarning("NTP") << "@SUB=analyze" << "No matching SC found for photon"; 
//        // *fTGoodEvent = 1; 
//        // break;
//      }
    }

    { //Look for associated PF objects

      // Initialization
      fTPhoisPFPhoton  ->push_back(0);
      fTPhoisPFElectron->push_back(0);

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
	(*fTPhoisPFPhoton)[phoqi] = 1;
	PhotonToPFPhotonOrElectronMatchingArray[phoqi] = iphot;
      }
      

      //Find PFElectron
      bool foundEgSC = false;
      int iel = -1;

      edm::Ptr<GsfElectron> elPtrSl;
      for (int i=0; i<(*fTNEles); i++){
	if (fTElSCindex->at(i)>=0 && fTElSCindex->at(i)==fTPhotSCindex->at(phoqi)) {
	  elPtrSl = elPtrVector.at(i);
	  foundEgSC = true;
	  break; 
	}
      }

      if (foundEgSC){
	double MVACut_ = -0.1; //42X                                                                                                                                                                                                                             
	//double MVACut_ = -1.; //44X 
	for( int i=0; i<ncand; ++i ) {
	  if ((*pfCandidates)[i].particleId()==reco::PFCandidate::e && (*pfCandidates)[i].gsfTrackRef().isNull()==false && (*pfCandidates)[i].mva_e_pi()>MVACut_ && (*pfCandidates)[i].gsfTrackRef()==elPtrSl->gsfTrack())
	    iel = i;
	}
      }

      if (iel!=-1) {
	(*fTPhoisPFElectron)[phoqi] = 1;
	PhotonToPFPhotonOrElectronMatchingArray[phoqi] = iel;
      }

    }



    { // CiC isolation
      fTPhoCiCPFIsoNeutralDR03->push_back(pfEcalIsoCiC(phoqi,*pfCandidates,0,0.3,0.,0.,0.,0.,0.,0.));
      fTPhoCiCPFIsoPhotonDR03->push_back(pfEcalIsoCiC(phoqi,*pfCandidates,2,0.3,0.,0.070,0.015,0.,0.,0.));
      fTPhoCiCPFIsoNeutralDR04->push_back(pfEcalIsoCiC(phoqi,*pfCandidates,0,0.4,0.,0.,0.,0.,0.,0.));
      fTPhoCiCPFIsoPhotonDR04->push_back(pfEcalIsoCiC(phoqi,*pfCandidates,2,0.4,0.,0.070,0.015,0.,0.,0.));
      for (int i=0; i<(*fTNVrtx); i++) {
	fTPhoCiCPFIsoChargedDR03->push_back(pfTkIsoWithVertexCiC(phoqi,i,*pfCandidates,1,0.3,0.02,0.02,0.0,0.2,0.1));
	fTPhoCiCPFIsoChargedDR04->push_back(pfTkIsoWithVertexCiC(phoqi,i,*pfCandidates,1,0.4,0.02,0.02,0.0,0.2,0.1));
      }
    }

    if (doPhotonStuff) { // Photon ID MVA

      photonIDMVA_variables.isrescaled = false;

      photonIDMVA_variables.pfchargedisobad03=-999;
      for (int ivtx=0; ivtx<*fTNVrtx; ivtx++)
	if (fTPhoCiCPFIsoChargedDR03->at(fTPhoVrtxListStart->at(phoqi)+ivtx)>photonIDMVA_variables.pfchargedisobad03)
	  photonIDMVA_variables.pfchargedisobad03=fTPhoCiCPFIsoChargedDR03->at(fTPhoVrtxListStart->at(phoqi)+ivtx);
      photonIDMVA_variables.pfphotoniso03=fTPhoCiCPFIsoPhotonDR03->at(phoqi);
      photonIDMVA_variables.pfneutraliso03=fTPhoCiCPFIsoNeutralDR03->at(phoqi);
      photonIDMVA_variables.sieie=fTPhoSigmaIetaIeta->at(phoqi);
      photonIDMVA_variables.sieip=fTPhoSigmaIetaIphi->at(phoqi);
      photonIDMVA_variables.etawidth=fTPhoSCEtaWidth->at(phoqi);
      photonIDMVA_variables.scrawe=fTPhoSCRawEnergy->at(phoqi);
      photonIDMVA_variables.phiwidth=fTPhoSCPhiWidth->at(phoqi);
      photonIDMVA_variables.r9=fTPhoR9->at(phoqi);
      photonIDMVA_variables.lambdaratio=fTPhoLambdaRatio->at(phoqi);
      photonIDMVA_variables.s4ratio=fTPhoS4Ratio->at(phoqi);
      photonIDMVA_variables.eventrho=*fTRho;
      photonIDMVA_variables.sceta=fTPhoSCEta->at(phoqi);
      photonIDMVA_variables.ESEffSigmaRR=fTPhoSigmaRR->at(phoqi);

      rescaleClusterShapes(photonIDMVA_variables, fTPhoisEB->at(phoqi));

      for (int ivtx=0; ivtx<*fTNVrtx; ivtx++){
	photonIDMVA_variables.pfchargedisogood03=fTPhoCiCPFIsoChargedDR03->at(fTPhoVrtxListStart->at(phoqi)+ivtx);
	fTPhoIDMVA->push_back((fTPhoisEB->at(phoqi)) ? photonIDMVA_reader_EB->EvaluateMVA("AdaBoost") : photonIDMVA_reader_EE->EvaluateMVA("AdaBoost"));
      }

    }


  } // end photon loop


  

  ///////////////////////////////////////////////////////
  // USAGE OF VERTEX CHOICE FOR DIPHOTON EVENTS:
  //
  // diphotons_{first,second} are vectors of {photon_1_index,photon_2_index}
  // vtx_dipho_??? are, for each diphoton pair, vectors of vertex indices (as ranked by the different algos)
  //
  // For example: best vertex for diphoton pair 3, with photon_1_index=diphotons_first[3] and photon_2_index=diphotons_second[3]: vtx_dipho_bla[3].at(0), second choice vtx_dipho_bla[3].at(1) ...
  //
  std::vector<int> diphotons_first;
  std::vector<int> diphotons_second;
  std::vector<std::vector<int> > vtx_dipho_h2gglobe;
  std::vector<std::vector<int> > vtx_dipho_mva;
  std::vector<std::vector<int> > vtx_dipho_productrank;
  
  if (doPhotonStuff) { // start vertex selection stuff with MVA from Hgg (Musella) UserCode/HiggsAnalysis/HiggsTo2photons/h2gglobe/VertexAnalysis tag vertex_mva_v4
    
    bool VTX_MVA_DEBUG = false;
    
    edm::Handle<reco::TrackCollection> tkH;
    iEvent.getByLabel(fTrackCollForVertexing, tkH);
    
    edm::Handle<VertexCollection> vtxH = vertices;
    
    std::vector<TVector3> vtxes;
    std::vector<TVector3> tracks_p3;
    std::vector<int> tkVtxId;
    std::vector<float> tk_d0;
    std::vector<float> tk_d0Err;
    std::vector<float> tk_dz;
    std::vector<float> tk_dzErr;
    std::vector<float> tk_PtErr;
    std::vector<bool> tk_ishighpurity;
    std::vector<std::vector<unsigned short> > vtx_std_tkind;
    std::vector<std::vector<float> > vtx_std_tkweight; // remember: [vertex][track]
    std::vector<int> vtx_std_ntks;
    
    
    { // tracks
      if (VTX_MVA_DEBUG)	   	   cout << "tracks begin" << endl;
      std::vector<reco::TrackBaseRef>::const_iterator tk;
          
      for(unsigned int i=0; i<vtxH->size(); i++) {

        if (VTX_MVA_DEBUG) cout << "working on vtx " << i << endl;

        reco::VertexRef vtx(vtxH, i);
	  
        if (VTX_MVA_DEBUG)	     	     cout << "vtx tracks " << vtx->tracksSize() << endl;
	     
        std::vector<unsigned short> temp;
        std::vector<float> temp_float;

        if (vtx->tracksSize()>0){
          for(tk=vtx->tracks_begin();tk!=vtx->tracks_end();++tk) {
            if (VTX_MVA_DEBUG)		 		 cout << "processing vtx track (out of " << vtx->tracksSize() << ")" << endl;
            int index = 0;
            for(reco::TrackCollection::size_type j = 0; j<tkH->size(); ++j) {
              if (VTX_MVA_DEBUG)		   		   std::cout << j << std::endl;
              reco::TrackRef track(tkH, j);
              if(TrackCut(track)) continue; 
              if (&(**tk) == &(*track)) {
                temp.push_back(index);
                temp_float.push_back(vtx->trackWeight(track));
                if (VTX_MVA_DEBUG)		     		     cout << "matching found index" << index << " weight " << vtx->trackWeight(track) << endl;
                break;
              }
              index++;
            }
          }
        }
        else {
          if (VTX_MVA_DEBUG)	       	       cout << "no vertex tracks found" << endl;
        }

	vtxes.push_back(TVector3(vtx->x(),vtx->y(),vtx->z()));
	vtx_std_ntks.push_back(temp.size());
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

        reco::TrackRef tk(tkH, i);

        if(TrackCut(tk))continue; 
	
	tracks_p3.push_back(TVector3(tk->px(),tk->py(),tk->pz()));
        tk_d0.push_back(tk->d0());
        tk_dz.push_back(tk->dz());
        tk_dzErr.push_back(tk->dzError());
        tk_d0Err.push_back(tk->d0Error());   
        tk_PtErr.push_back(tk->ptError());
        tk_ishighpurity.push_back(tkIsHighPurity(tk));
        tkVtxId.push_back(-1); 
	 
      } // for i (loop over all tracks)

      if (VTX_MVA_DEBUG)	  cout << "done tracks" << endl;
    }

    (*fTNconv)=0;

    { // all conversions

      edm::Handle<reco::ConversionCollection> hConversions;
      iEvent.getByLabel("allConversions", hConversions);

      for( reco::ConversionCollection::const_iterator  iConv = hConversions->begin(); iConv != hConversions->end(); iConv++) {

        reco::Conversion localConv = reco::Conversion(*iConv);
  
        if(ConversionsCut(localConv)) continue;
  
	if (*fTNconv >= gMaxNConv){
	  edm::LogWarning("NTP") << "@SUB=analyze"
				 << "Maximum number of conversions exceeded";
	  *fTGoodEvent = 1;
	  break;
	}

	fTConvValidVtx->push_back(false);
	fTConvNtracks           ->push_back(-999);
	fTConvChi2Probability   ->push_back(-999);
	fTConvEoverP            ->push_back(-999);
	fTConvZofPrimVtxFromTrks->push_back(-999);
	conv_refitted_momentum[*fTNconv].SetXYZ(-999,-999,-999);
	conv_singleleg_momentum[*fTNconv].SetXYZ(-999,-999,-999);

	fTConvValidVtx->at(*fTNconv)=localConv.conversionVertex().isValid();
        if ( localConv.conversionVertex().isValid() ) {
	  reco::Vertex vtx=localConv.conversionVertex();
	  conv_vtx[*fTNconv].SetXYZ(vtx.x(), vtx.y(), vtx.z());
	  fTConvNtracks           ->at(*fTNconv)=(localConv.nTracks());	  
	  fTConvChi2Probability   ->at(*fTNconv)=(ChiSquaredProbability(vtx.chi2(), vtx.ndof()));
	  //	  fTConvEoverP            ->at(*fTNconv)=(localConv.EoverPrefittedTracks()); // it is commented out in globe
	  fTConvZofPrimVtxFromTrks->at(*fTNconv)=(localConv.zOfPrimaryVertexFromTracks());
	  conv_refitted_momentum[*fTNconv].SetXYZ(localConv.refittedPairMomentum().x(), localConv.refittedPairMomentum().y(), localConv.refittedPairMomentum().z());
	  conv_singleleg_momentum[*fTNconv].SetXYZ(-999,-999,-999);
	}

        (*fTNconv)++;  
      }


      for (edm::View<reco::Photon>::const_iterator  localPho=photons->begin(); localPho!=photons->end(); localPho++) {
	for (std::vector<reco::Photon>::const_iterator  iPfCand=pfPhotonHandle->begin(); iPfCand!=pfPhotonHandle->end(); iPfCand++) {
      if (localPho->superCluster()!=iPfCand->superCluster()) continue;
      reco::ConversionRefVector convsingleleg = iPfCand->conversionsOneLeg();
      for (unsigned int iconvoneleg=0; iconvoneleg<convsingleleg.size(); iconvoneleg++){
      
      reco::Conversion localConv = reco::Conversion(*convsingleleg[iconvoneleg]);

      if(ConversionsCut(localConv)) continue;

        if (*fTNconv >= gMaxNConv){
          edm::LogWarning("NTP") << "@SUB=analyze"
                                 << "Maximum number of conversions exceeded";
				 *fTGoodEvent = 1;
          break;
        }

	fTConvValidVtx->push_back(false);
	fTConvNtracks           ->push_back(-999);
	fTConvChi2Probability   ->push_back(-999);
	fTConvEoverP            ->push_back(-999);
	fTConvZofPrimVtxFromTrks->push_back(-999);
	conv_refitted_momentum[*fTNconv].SetXYZ(-999,-999,-999);
	conv_singleleg_momentum[*fTNconv].SetXYZ(-999,-999,-999);

	fTConvValidVtx->at(*fTNconv)=localConv.conversionVertex().isValid();
        if ( localConv.conversionVertex().isValid() ) {
          reco::Vertex vtx=localConv.conversionVertex();
          conv_vtx[*fTNconv].SetXYZ(vtx.x(), vtx.y(), vtx.z());
          fTConvNtracks           ->at(*fTNconv)=(localConv.nTracks());
          fTConvChi2Probability   ->at(*fTNconv)=(ChiSquaredProbability(vtx.chi2(), vtx.ndof()));
	  //          fTConvEoverP            ->at(*fTNconv)=(localConv.EoverPrefittedTracks()); // it is commented out in globe
          fTConvZofPrimVtxFromTrks->at(*fTNconv)=(localConv.zOfPrimaryVertexFromTracks());
          conv_refitted_momentum[*fTNconv].SetXYZ(localConv.refittedPairMomentum().x(), localConv.refittedPairMomentum().y(), localConv.refittedPairMomentum().z());
	  conv_singleleg_momentum[*fTNconv].SetXYZ(-999,-999,-999);
	  if( localConv.nTracks()) {
	    const std::vector<edm::RefToBase<reco::Track> > tracks = localConv.tracks();
	    conv_singleleg_momentum[*fTNconv].SetXYZ(tracks[0]->px(), tracks[0]->py(), tracks[0]->pz());
	  }
        }

        (*fTNconv)++;
	
	}
	
	}
      }
	

    }

    if (VTX_MVA_DEBUG)       cout << "done convs" << endl;

    ETHVertexInfo vinfo(vtxes.size(),&vtxes,tracks_p3.size(),&tracks_p3,&tk_PtErr,&tkVtxId,&tk_d0,&tk_d0Err,&tk_dz,&tk_dzErr,&tk_ishighpurity,&vtx_std_tkind,&vtx_std_tkweight,&vtx_std_ntks);
    
    if (VTX_MVA_DEBUG)	 cout << "ready" << endl;	 
    
    // fully combinatorial vertex selection
    for (int ip=0; ip<(*fTNPhotons); ++ip) {
    for (int jp=ip+1; jp<(*fTNPhotons); ++jp) {
    diphotons_first.push_back(ip);
    diphotons_second.push_back(jp);
    }
    }
    
    
    for(unsigned int id=0; id<diphotons_first.size(); ++id ) {
    
    if (VTX_MVA_DEBUG)	     cout << "processing diphoton pair " << id << endl;
    
    int ipho1 = diphotons_first[id];
    int ipho2 = diphotons_second[id];
    
    PhotonInfo pho1=fillPhotonInfos(ipho1,vtxAlgoParams.useAllConversions,(*fTPhoRegrEnergy)[ipho1]);
    PhotonInfo pho2=fillPhotonInfos(ipho2,vtxAlgoParams.useAllConversions,(*fTPhoRegrEnergy)[ipho2]);

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
    vtx_dipho_h2gglobe.push_back(HggVertexSelection(*vAna, *vConv, pho1, pho2, perVtxVariables, mvaVertexSelection, perVtxReader, perVtxMvaMethod));

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
 
    for (int i=0; i<(int)diphotons_first.size() && i<gMax_vertexing_diphoton_pairs; i++) { // write output of vertexing
      fTDiphotonsfirst->push_back(diphotons_first.at(i));
      fTDiphotonssecond->push_back(diphotons_second.at(i));
      fTVtxdiphoh2gglobe->push_back(vtx_dipho_h2gglobe.at(i).at(0));
      fTVtxdiphomva->push_back(vtx_dipho_mva.at(i).at(0));
      fTVtxdiphoproductrank->push_back(vtx_dipho_productrank.at(i).at(0));
    }

    vAna->clear();

  } // end vertex selection for diphoton events

  //       cout << "end vertex selection MVA" << endl;



  std::vector<std::vector<int> > Jets_PfCand_content;

  ////////////////////////////////////////////////////////
  // Jet Variables:
  const JetCorrector* jetCorr = JetCorrector::getJetCorrector(fJetCorrs, iSetup);
  std::vector<OrderPair> corrIndices;  // Vector of indices and pt of corr. jets to re-order them
  int iraw(0);
  *fTNJetsTot = jets->size();
  // Loop over uncorr. jets
  for(View<Jet>::const_iterator Jit = jets->begin(); Jit != jets->end(); ++Jit, ++iraw){
    // Cut on uncorrected pT (for startup)
    if(Jit->pt() < fMinRawJPt) continue;
    double scale = jetCorr->correction(*Jit,iEvent,iSetup);
    corrIndices.push_back(make_pair(iraw, scale*Jit->pt()));
  }
	
  // Sort corrected jet collection by decreasing pt
  std::sort(corrIndices.begin(), corrIndices.end(), indexComparator);
	
  // Determine corrected jets
  int jqi(-1); // counts # of qualified jets
  // Loop over corr. jet indices
  for(std::vector<OrderPair>::const_iterator it = corrIndices.begin(); 
      it != corrIndices.end(); ++it ) {
    // Check if maximum number of jets is exceeded already
    if(jqi >= gMaxNJets-1) {
      edm::LogWarning("NTP") << "@SUB=analyze"
                             << "Maximum number of jets exceeded";
      *fTMaxJetExceed = 1;
      *fTGoodEvent = 1;
      break;
    }
    int index = it->first;
    const PFJet* cojet = static_cast<const PFJet*>( &((*jets)[index]) ); // look away...
    std::auto_ptr<PFJet> jet(new PFJet(*cojet));
    
    // The correction was calculated above: use it
    double scale = it->second/jet->pt();
    jet->scaleEnergy(scale);
	
    // Jet preselection
    if(jet->pt() < fMinCorJPt) continue; 
    if(fabs(jet->eta()) > fMaxJEta) continue;
    jqi++;

    // Dump jet properties into tree variables
    fTJPx     ->push_back(jet->px());
    fTJPy     ->push_back(jet->py());
    fTJPz     ->push_back(jet->pz());
    fTJPt     ->push_back(jet->pt());
    fTJEta    ->push_back(jet->eta());
    fTJPhi    ->push_back(jet->phi());
    fTJE      ->push_back(jet->energy());
    fTJEt     ->push_back(jet->et());
    fTJEcorr  ->push_back(scale);
    fTJEtaRms ->push_back(sqrt(jet->etaetaMoment()));
    fTJPhiRms ->push_back(sqrt(jet->etaphiMoment()));
    fTJArea   ->push_back(jet->jetArea());

    fTJNConstituents->push_back(jet->nConstituents());
    fTJNAssoTracks  ->push_back(jet->chargedMultiplicity()); // do it the pf way...
    fTJNNeutrals    ->push_back(jet->neutralMultiplicity()); 
		
    // energy fractions for JID need to be computed w.r.t. uncorrected jet energy!!
    // see for instance https://twiki.cern.ch/twiki/bin/view/CMS/JetID
    // or http://cmssdt.cern.ch/SDT/doxygen/CMSSW_4_1_3/doc/html/dc/dd5/classPFJetIDSelectionFunctor.html  
    double uncorr_energy  = jet->energy()/scale;
    fTJChargedHadFrac     ->push_back(jet->chargedHadronEnergy()/uncorr_energy);
    fTJNeutralHadFrac     ->push_back(jet->neutralHadronEnergy()/uncorr_energy + jet->HFHadronEnergy()/uncorr_energy);
    fTJChargedEmFrac      ->push_back(jet->chargedEmEnergy()/uncorr_energy);
    fTJNeutralEmFrac      ->push_back(jet->neutralEmEnergy()/uncorr_energy);
    fTJChargedMuEnergyFrac->push_back(jet->chargedMuEnergy()/uncorr_energy);
    fTJPhoFrac            ->push_back(jet->photonEnergy()/uncorr_energy); // photons also count for neutralEmEnergy
    fTJHFHadFrac          ->push_back(jet->HFHadronEnergy()/uncorr_energy);
    fTJHFEMFrac           ->push_back(jet->HFEMEnergy()/uncorr_energy);   // also contained in neutralEmEnergy
    // see CMSSW/RecoJets/JetProducers/src/JetSpecific.cc

    vector<PFCandidatePtr> JetpfCandidates = jet->getPFConstituents();
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

    fTJPtD       ->push_back(sqrt( sumPt2_cands )/sumPt_cands);
    fTJRMSCand   ->push_back(rms_cands/sumPt2_cands);



    // Calculate the DR wrt the closest electron
    float ejDRmin = 10.; // Default when no electrons previously selected
    for( int j = 0; j < (*fTNEles); j++ ){
      float ejDR = reco::deltaR(jet->eta(), jet->phi(), (*fTElEta)[j], (*fTElPhi)[j]);
      if(ejDR<ejDRmin) ejDRmin = ejDR;
    }
    fTJeMinDR ->push_back(ejDRmin);

    // B-tagging probability
    // remember: 'index' is the index of the uncorrected jet, as saved in the event
    ibtag = 0;
    for ( std::vector<edm::InputTag>::const_iterator it = fBtagTags.begin();
          it != fBtagTags.end(); ++it ) {
      fTJbTagProb[ibtag]->push_back((*jetsBtag[ibtag])[index].second);
      ++ibtag;
    }

      int myFlavour=0;// unassigned == 0
      RefToBase<reco::Jet> aJet;
      aJet = (*jetsBtag[0])[index].first; // fill it from the collection you want to probe!
      if (flavours.find (aJet) == flavours.end()) {
          //std::cout <<" Cannot access flavour for this jet " << index << " - not in the Map"<<std::endl;
	  fTJPartonFlavour->push_back(0);
      } else {
          myFlavour = flavours[aJet];
          //std::cout << " Jet " << index << " has flavour = " << flavours[aJet] << " check " << myFlavour << std::endl;//comment after checks are done
	  fTJPartonFlavour->push_back(myFlavour);
      }


    // start computation of betaStar variable (pileUp ID) -- adding also beta variable, which cuts on the dz rather than the vertex association (marc feb5 2013)
    float sumTrkPt = 0.;
    float sumTrkPtBetaStar = 0.;
    float sumTrkPtBeta = 0.;
    float sumTrkPtSq = 0.;
    float sumTrkPtBetaSq = 0.;

    // Jet-track association: get associated tracks
    const reco::TrackRefVector& tracks = jet->getTrackRefs();
    std::vector<const reco::Track*> AssociatedTracks;

    for( TrackRefVector::iterator i_trk = tracks.begin(); i_trk != tracks.end(); ++i_trk )  { 

      // calculate first beta, then move on to beta*
      sumTrkPtSq += (*i_trk)->pt()*(*i_trk)->pt();
      AssociatedTracks.push_back( i_trk->get() );

      if ( vertices->size() == 0) continue;
      sumTrkPt += (*i_trk)->pt();
      
      // check if track is associated to primary vertex
      bool isFirstVtx=false;
      // loop over the tracks associated with the vertex
      if (!((*vertices)[0].isFake()) && (*vertices)[0].ndof() >= 4 && fabs((*vertices)[0].z()) <= 24.) {
        for(reco::Vertex::trackRef_iterator i_vtxTrk = (*vertices)[0].tracks_begin(); i_vtxTrk != (*vertices)[0].tracks_end(); ++i_vtxTrk) {
          // match the jet track to the track from the vertex 
          reco::TrackRef trkRef(i_vtxTrk->castTo<reco::TrackRef>());
          // check if the tracks match
          if (trkRef == (*i_trk)) {
            isFirstVtx=true; 
            // for the beta calculation. if the track is associated to the PV, cut on 0.5 cm of dz
            if ((*i_trk)->dz((*vertices)[0].position()) < 0.5) sumTrkPtBetaSq += (*i_trk)->pt()*(*i_trk)->pt();
            if ((*i_trk)->dz((*vertices)[0].position()) < 0.5) sumTrkPtBeta += (*i_trk)->pt();
            break;
          }
        }
      }
      
      // if not associated to primary vertex, check other vertices:
      bool isOtherVtx = false;
      if (!isFirstVtx) {   
        for(unsigned iotherVtx=1; iotherVtx<vertices->size();iotherVtx++) {
          if (!((*vertices)[iotherVtx].isFake()) && (*vertices)[iotherVtx].ndof() >= 4 && fabs((*vertices)[iotherVtx].z()) <= 24.) {
            for(reco::Vertex::trackRef_iterator i_vtxTrk = (*vertices)[iotherVtx].tracks_begin(); 
                i_vtxTrk != (*vertices)[iotherVtx].tracks_end(); ++i_vtxTrk) {
              reco::TrackRef trkRef(i_vtxTrk->castTo<reco::TrackRef>());
              if (trkRef == (*i_trk) ) {
                isOtherVtx=true;
                break;
              }
            } // for tracks
          } // if vtx is good
        } // for other vertices
      } // if !isFirstVt

      if(!isFirstVtx && isOtherVtx) { sumTrkPtBetaStar += (*i_trk)->pt(); } 

    } // for tracks
    
    
    float betaStar = -999.;
    float beta = -999.;
    float betaSq = -999.;
    if (sumTrkPt > 0.) 
      betaStar = sumTrkPtBetaStar/sumTrkPt;
    fTJBetaStar->push_back( betaStar ); 
    if (sumTrkPt   > 0.) beta = sumTrkPtBeta  /sumTrkPt  ; fTJBeta->push_back( beta ); 
    if (sumTrkPtSq > 0.) betaSq = sumTrkPtBetaSq/sumTrkPtSq; fTJBetaSq->push_back( betaSq ); 

			
    // Below save the momenta of the three leading tracks associated to the jet
    float pT1(0.), pT2(0.), pT3(0.);
    int idx1(-1), idx2(-1), idx3(-1);
			
    // Jet-track association: make transient tracks and store information
    std::vector<TransientTrack> AssociatedTTracks;
    // Initialization
    fTJMass ->push_back(0.);
    fTJtrk1px->push_back(-999.99 );
    fTJtrk1py->push_back(-999.99 );
    fTJtrk1pz->push_back(-999.99 );
    fTJtrk2px->push_back(-999.99 );
    fTJtrk2py->push_back(-999.99 );
    fTJtrk2pz->push_back(-999.99 );
    fTJtrk3px->push_back(-999.99 );
    fTJtrk3py->push_back(-999.99 );
    fTJtrk3pz->push_back(-999.99 );
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
        (*fTJtrk1px)[jqi] = AssociatedTracks[idx1]->px();
        (*fTJtrk1py)[jqi] = AssociatedTracks[idx1]->py();
        (*fTJtrk1pz)[jqi] = AssociatedTracks[idx1]->pz();
      }
      if(AssociatedTracks.size()>=2){
        (*fTJtrk2px)[jqi] = AssociatedTracks[idx2]->px();
        (*fTJtrk2py)[jqi] = AssociatedTracks[idx2]->py();
        (*fTJtrk2pz)[jqi] = AssociatedTracks[idx2]->pz();
      }
      if(AssociatedTracks.size()>=3){
        (*fTJtrk3px)[jqi] = AssociatedTracks[idx3]->px();
        (*fTJtrk3py)[jqi] = AssociatedTracks[idx3]->py();
        (*fTJtrk3pz)[jqi] = AssociatedTracks[idx3]->pz();
      }
			
      (*fTJMass)[jqi]   = sqrt(E2tmp - pXtmp*pXtmp - pYtmp*pYtmp - pZtmp*pZtmp);
      // apparently there ARE cases where ChMult is 0, but we still end up here...
      if((*fTJNAssoTracks)[jqi] > 0) (*fTJMass)[jqi] *= (*fTJNConstituents)[jqi]/(*fTJNAssoTracks)[jqi];
    } else { // The whole cone used for jet-tracks association is outside of the tracker acceptance
      (*fTJMass)[jqi] = -888.88;
    }

    // Do a vertex fitting with the tracks
    if(AssociatedTTracks.size() > 1) {
      TransientVertex jetVtx = avFitter.vertex(AssociatedTTracks);
      if(jetVtx.isValid()){
        fTJVtxx     ->push_back(jetVtx.position().x());
        fTJVtxy     ->push_back(jetVtx.position().y());
        fTJVtxz     ->push_back(jetVtx.position().z());
        fTJVtxExx   ->push_back(jetVtx.positionError().cxx());
        fTJVtxEyx   ->push_back(jetVtx.positionError().cyx());
        fTJVtxEyy   ->push_back(jetVtx.positionError().cyy());
        fTJVtxEzy   ->push_back(jetVtx.positionError().czy());
        fTJVtxEzz   ->push_back(jetVtx.positionError().czz());
        fTJVtxEzx   ->push_back(jetVtx.positionError().czx());
        fTJVtxNChi2 ->push_back(jetVtx.normalisedChiSquared());
      }else{
        fTJVtxx     ->push_back(-777.77);
        fTJVtxy     ->push_back(-777.77);
        fTJVtxz     ->push_back(-777.77);
        fTJVtxExx   ->push_back(-777.77);
        fTJVtxEyx   ->push_back(-777.77);
        fTJVtxEyy   ->push_back(-777.77);
        fTJVtxEzy   ->push_back(-777.77);
        fTJVtxEzz   ->push_back(-777.77);
        fTJVtxEzx   ->push_back(-777.77);
        fTJVtxNChi2 ->push_back(-777.77);
      }
    }else{
      fTJVtxx     ->push_back(-888.88);
      fTJVtxy     ->push_back(-888.88);
      fTJVtxz     ->push_back(-888.88);
      fTJVtxExx   ->push_back(-888.88);
      fTJVtxEyx   ->push_back(-888.88);
      fTJVtxEyy   ->push_back(-888.88);
      fTJVtxEzy   ->push_back(-888.88);
      fTJVtxEzz   ->push_back(-888.88);
      fTJVtxEzx   ->push_back(-888.88);
      fTJVtxNChi2 ->push_back(-888.88);
    }
    AssociatedTracks.clear();
    AssociatedTTracks.clear();
	
    // GenJet matching
    if (!fIsRealData && (*fTNGenJets) > 0) fTJGenJetIndex->push_back( matchJet(&(*jet)) );
    fTJGood->push_back( 0 );

    fTJVrtxListStart->push_back((fTJVrtxListStart->size()==0) ? 0 : fTJVrtxListStart->back()+(*fTNVrtx));

    if (doPhotonStuff){

      reco::VertexCollection vtxColl = *(vertices.product());

      if (PileupJetIdAlgos.size()>0) for (int ivtx=0; ivtx<(*fTNVrtx); ivtx++){
	  
	PileupJetIdentifier jetIdentifer_vars = PileupJetIdAlgos.at(0)->computeIdVariables( &(*jet), scale, Ptr<Vertex>(vertices,ivtx).get(), vtxColl );

	  for (uint i=0; i<PileupJetIdAlgos.size(); i++){
	    PileupJetIdAlgo* ialgo = PileupJetIdAlgos.at(i);
	    ialgo->set(jetIdentifer_vars);
	    PileupJetIdentifier id = ialgo->computeMva();
	    fTJPassPileupIDL[i]->push_back(false);
	    fTJPassPileupIDM[i]->push_back(false);
	    fTJPassPileupIDT[i]->push_back(false);
	    int idflag = id.idFlag();
	    if( PileupJetIdentifier::passJetId( idflag, PileupJetIdentifier::kLoose )) {
	      fTJPassPileupIDL[i]->back()=true;
	    }
	    if( PileupJetIdentifier::passJetId( idflag, PileupJetIdentifier::kMedium )) {
	      fTJPassPileupIDM[i]->back()=true;
	    }
	    if( PileupJetIdentifier::passJetId( idflag, PileupJetIdentifier::kTight )) {
	      fTJPassPileupIDT[i]->back()=true;
	    }
	    
	  } 
	}
    }
    
    if (doPhotonStuff) { // QG tagging
      edm::RefToBase<reco::Jet> jetRef(jets->refAt(index));
      fTJQGTagLD->push_back(-999);
      fTJQGTagMLP->push_back(-999);
      fTJSmearedQGL->push_back(-999);
      if (QGTagsHandleLikelihood.isValid()) {
	fTJQGTagLD->back()=(*QGTagsHandleLikelihood)[jetRef];

	// see https://twiki.cern.ch/twiki/bin/viewauth/CMS/GluonTag for further information
	// type is the jet flavor and must be either "quark" or "gluon" or "all"
	// (to be matched at parton level, i.e. status 3; if no match or pileup set to all)

	fTJSmearedQGL->back()=fTJQGTagLD->back();
	
	if (!fIsRealData){
	  std::string jet_type = "all";
	  for (reco::GenParticleCollection::const_iterator gpart = GlobalGenParticles->begin(); gpart != GlobalGenParticles->end(); gpart++){
	    if( gpart->status() != 3 ) continue;
	    if( (gpart->pdgId()<-4) || (gpart->pdgId()>4 && gpart->pdgId()!=21)) continue;
	    double dr = reco::deltaR(gpart->eta(), gpart->phi(), fTJEta->back(), fTJPhi->back());
	    if(dr > 0.3) continue;
	    double ndpt = fabs(gpart->pt() - fTJPt->back())/gpart->pt();
	    if(ndpt > 2.) continue;
	    if (gpart->pdgId()==21) jet_type="gluon";
	    else jet_type="quark";
	    break;
	  }
	  fTJSmearedQGL->back() = qgsyst.Smear(fTJPt->back(), fTJEta->back(), *rhoForQG, fTJQGTagLD->back(), jet_type);
	}

      }

      if (QGTagsHandleMLP.isValid()) {
	fTJQGTagMLP->back()=(*QGTagsHandleMLP)[jetRef];
      }

    }

  }
  (*fTNJets) = jqi+1;
  corrIndices.clear();




  ////////////////////////////////////////////////////////
  // PfCandidates Variables:

  int pfcandIndex(0);
  for (unsigned int i=0; i<pfCandidates->size() && doPhotonStuff; i++){

    if(fabs((*pfCandidates)[i].eta()) > fMaxPfCandEta) continue;

    int type = FindPFCandType((*pfCandidates)[i].pdgId());
    if (type==2) storethispfcand[i]=true;

    for (int j=0; j<(*fTNPhotons); j++){
      if (PhotonToPFPhotonOrElectronMatchingArray[j]==(int)i) {
	PhotonToPFPhotonOrElectronMatchingArrayTranslator[j]=pfcandIndex;
	storethispfcand[i]=true;
      }
    }

    for (int j=0; j<(*fTNPhotons); j++){
      for (size_t k=0; k<list_pfcand_footprint.at(j).size(); k++){
	if (list_pfcand_footprint.at(j).at(k)==(int)i){
	  list_pfcand_footprintTranslator.at(j).push_back(pfcandIndex);
	  storethispfcand[i]=true;
	}
      }
    }
  

    if (storethispfcand[i]==false) continue;

    if (pfcandIndex >= gMaxNPfCand){
      edm::LogWarning("NTP") << "@SUB=analyze"
			     << "Maximum number of pf candidates exceeded";
      *fTGoodEvent = 1;
      break;
    }

    fTPfCandPdgId->push_back( (*pfCandidates)[i].pdgId() );
    fTPfCandPt->push_back( (*pfCandidates)[i].pt() );
    fTPfCandEta->push_back( (*pfCandidates)[i].eta() );
    fTPfCandPhi->push_back( (*pfCandidates)[i].phi() );
    fTPfCandEnergy->push_back( (*pfCandidates)[i].energy() );
    fTPfCandEcalEnergy->push_back( (*pfCandidates)[i].ecalEnergy() );
    fTPfCandVx->push_back( (*pfCandidates)[i].vx() );
    fTPfCandVy->push_back( (*pfCandidates)[i].vy() );
    fTPfCandVz->push_back( (*pfCandidates)[i].vz() );

//    if ( (type==1) && ( !((*pfCandidates)[i].trackRef()) ) ) type=-1;
//    reco::HitPattern pattern; 
//    if (type==1) pattern=(*pfCandidates)[i].trackRef()->hitPattern(); 
//    fTPfCandHasHitInFirstPixelLayer->push_back( (type==1) ? (pattern.hasValidHitInFirstPixelBarrel() || pattern.hasValidHitInFirstPixelEndcap()) : -999 );
//    fTPfCandTrackRefPx->push_back( (type==1) ? (*pfCandidates)[i].trackRef()->px() : -999 );
//    fTPfCandTrackRefPy->push_back( (type==1) ? (*pfCandidates)[i].trackRef()->py() : -999 );
//    fTPfCandTrackRefPz->push_back( (type==1) ? (*pfCandidates)[i].trackRef()->pz() : -999 );

    fTPfCandBelongsToJet->push_back(-999);
    for (int jet=0; jet<(int)(Jets_PfCand_content.size()); jet++){
      for (int k=0; k<(int)(Jets_PfCand_content.at(jet).size()); k++){
	if (Jets_PfCand_content.at(jet).at(k)==(int)i) {
	  if (fTPfCandBelongsToJet->back()>=0) {cout << "WRONG: this PFCand is already in another jet!" << endl; continue;}
	  fTPfCandBelongsToJet->back()=jet;
	}
      }
    }

    pfcandIndex++;

  }
  *fTNPfCand=pfcandIndex;
  
  for (int j=0; j<(*fTNPhotons) && doPhotonStuff; j++){
    fTPhoMatchedPFPhotonOrElectronCand->push_back(PhotonToPFPhotonOrElectronMatchingArrayTranslator[j]);
    fTPhoFootprintPfCandsListStart->push_back(fTPhoFootprintPfCands->size());
    for (size_t k=0; k<list_pfcand_footprintTranslator.at(j).size(); k++) fTPhoFootprintPfCands->push_back(list_pfcand_footprintTranslator.at(j).at(k));
  }


  if (!doPhotonStuff){
  ////////////////////////////////////////////////////////
  // Get and Dump (M)E(T) Variables:
  // Tracks:
  int nqtrk(-1);
  *fTNTracksTot = tracks->size();
  *fTTrkPtSumx = 0.; *fTTrkPtSumy = 0.;
  for( TrackCollection::const_iterator it = tracks->begin(); it != tracks->end() ; ++it ){
    *fTTrkPtSumx += it->px(); // Calculated for ALL tracks
    *fTTrkPtSumy += it->py();
    if(it->pt() < fMinTrkPt) continue;
    if(fabs(it->eta()) > fMaxTrkEta) continue;
    if(it->normalizedChi2() > fMaxTrkNChi2) continue;
    if(it->numberOfValidHits() < fMinTrkNHits) continue;
    nqtrk++; // starts at 0
    // Check if maximum number of tracks is exceeded already
    if(nqtrk >= gMaxNTrks) {
      edm::LogWarning("NTP") << "@SUB=analyze"
                             << "Maximum number of tracks exceeded";
      *fTMaxTrkExceed = 1;
      *fTGoodEvent = 1;
      break;
    }
    fTTrkPt    ->push_back(it->pt()*it->charge());
    fTTrkEta   ->push_back(it->eta());
    fTTrkPhi   ->push_back(it->phi());
    fTTrkNChi2 ->push_back(it->normalizedChi2());
    fTTrkNHits ->push_back(it->numberOfValidHits());
    fTTrkVtxDz ->push_back((primVtx) ? it->dz(primVtx->position()) : 999);
    fTTrkVtxDxy->push_back((primVtx) ? it->dxy(primVtx->position()) : 999);	
    fTTrkGood  ->push_back(0);
  }
  *fTNTracks = nqtrk+1;

  *fTTrkPtSum = sqrt((*fTTrkPtSumx)*(*fTTrkPtSumx) + (*fTTrkPtSumy)*(*fTTrkPtSumy));
  TVector3 trkPtSum(*fTTrkPtSumx, *fTTrkPtSumy, 0.);
  *fTTrkPtSumPhi = trkPtSum.Phi();

  // Calotowers:
  (*fTNCaloTowers) = calotowers->size();
  *fTECALEsumx = 0.; *fTECALEsumy = 0.; *fTECALEsumz = 0.;
  *fTHCALEsumx = 0.; *fTHCALEsumy = 0.; *fTHCALEsumz = 0.;
  *fTSumEt = 0.; *fTECALSumEt = 0.; *fTHCALSumEt = 0.;
  for(CaloTowerCollection::const_iterator itow = calotowers->begin();
      itow!=calotowers->end(); ++itow ){
    if(!(itow->energy()>0.)) continue; // Check against zero energy towers
    *fTSumEt += itow->et();
    *fTECALSumEt += itow->emEt();
    *fTHCALSumEt += itow->hadEt();
    double emFrac = itow->emEnergy()/itow->energy();
    double hadFrac = itow->hadEnergy()/itow->energy();
    *fTECALEsumx += itow->px()*emFrac;
    *fTECALEsumy += itow->py()*emFrac;
    *fTECALEsumz += itow->pz()*emFrac;
    *fTHCALEsumx += itow->px()*hadFrac;
    *fTHCALEsumy += itow->py()*hadFrac;
    *fTHCALEsumz += itow->pz()*hadFrac;
  }
  TVector3 ecalMET(*fTECALEsumx, *fTECALEsumy, *fTECALEsumz);
  TVector3 hcalMET(*fTHCALEsumx, *fTHCALEsumy, *fTHCALEsumz);
  *fTECALMET    = ecalMET.Mag();
  *fTECALMETPhi = ecalMET.Phi();
  *fTHCALMET    = hcalMET.Mag();
  *fTHCALMETPhi = hcalMET.Phi();
  if(*fTECALEsumz != 0.) *fTECALMETEta = ecalMET.Eta();
  else *fTECALMETEta = 0.;
  if(*fTHCALEsumz != 0.) *fTHCALMETeta = hcalMET.Eta();
  else *fTHCALMETeta = 0.;

  // MET Collections:
  *fTRawMET             = (calomet->at(0)).pt();
  *fTRawMETpx           = (calomet->at(0)).px();
  *fTRawMETpy           = (calomet->at(0)).py();
  *fTRawMETphi          = (calomet->at(0)).phi();
  *fTRawMETemEtFrac     = (calomet->at(0)).emEtFraction();
  *fTRawMETemEtInEB     = (calomet->at(0)).emEtInEB();
  *fTRawMETemEtInEE     = (calomet->at(0)).emEtInEE();
  *fTRawMETemEtInHF     = (calomet->at(0)).emEtInHF();
  *fTRawMEThadEtFrac    = (calomet->at(0)).etFractionHadronic();
  *fTRawMEThadEtInHB    = (calomet->at(0)).hadEtInHB();
  *fTRawMEThadEtInHE    = (calomet->at(0)).hadEtInHE();
  *fTRawMEThadEtInHF    = (calomet->at(0)).hadEtInHF();
  sigmaX2= (calomet->front() ).getSignificanceMatrix()(0,0);
  sigmaY2= (calomet->front() ).getSignificanceMatrix()(1,1);
  significance = 0;
  if(sigmaX2<1.e10 && sigmaY2<1.e10) significance = (calomet->front() ).significance();
  *fTRawMETSignificance = significance;

  if (!fIsRealData) {
    Handle<View<GenMET> > GenMET;
    iEvent.getByLabel(fGenMETTag, GenMET);
    *fTGenMET    = (GenMET->front()).pt();
    *fTGenMETpx  = (GenMET->front()).px();
    *fTGenMETpy  = (GenMET->front()).py();
    *fTGenMETphi = (GenMET->front()).phi();
  }

  *fTTCMET    = (tcmet->at(0)).pt();
  *fTTCMETpx  = (tcmet->at(0)).px();
  *fTTCMETpy  = (tcmet->at(0)).py();
  *fTTCMETphi = (tcmet->at(0)).phi();
  sigmaX2= (tcmet->front() ).getSignificanceMatrix()(0,0);
  sigmaY2= (tcmet->front() ).getSignificanceMatrix()(1,1);
  significance = 0;
  if(sigmaX2<1.e10 && sigmaY2<1.e10) significance = (tcmet->front() ).significance();
  *fTTCMETSignificance = significance;

  *fTPFMET    = (pfmet->front()).pt();
  *fTPFMETpx  = (pfmet->front()).px();
  *fTPFMETpy  = (pfmet->front()).py();
  *fTPFMETphi = (pfmet->front()).phi();
  sigmaX2= (pfmet->front() ).getSignificanceMatrix()(0,0);
  sigmaY2= (pfmet->front() ).getSignificanceMatrix()(1,1);
  significance = 0;
  if(sigmaX2<1.e10 && sigmaY2<1.e10) significance = (pfmet->front() ).significance();
  *fTPFMETSignificance = significance;
  *fTPFSumEt  = (pfmet->front()).sumEt();

  *fTMuJESCorrMET    = (corrmujesmet->at(0)).pt();
  *fTMuJESCorrMETpx  = (corrmujesmet->at(0)).px();
  *fTMuJESCorrMETpy  = (corrmujesmet->at(0)).py();
  *fTMuJESCorrMETphi = (corrmujesmet->at(0)).phi();

  if(*fTNJets > 1){
    double dPhiMJ1 = TMath::Abs(reco::deltaPhi((*fTJPhi)[0], *fTMuJESCorrMETphi));
    double dPhiMJ2 = TMath::Abs(reco::deltaPhi((*fTJPhi)[1], *fTMuJESCorrMETphi));
    *fTMETR12 = TMath::Sqrt(dPhiMJ1*dPhiMJ1 + (TMath::Pi()-dPhiMJ2)*(TMath::Pi()-dPhiMJ2) );
    *fTMETR21 = TMath::Sqrt(dPhiMJ2*dPhiMJ2 + (TMath::Pi()-dPhiMJ1)*(TMath::Pi()-dPhiMJ1) );
  }

  // Information to be able to apply MET corrections on-the-fly
  // Need to store Pt, Eta, EMfraction and Area of raw Jet, after muon subtraction
  // See JetMETCorrections/Type1MET/interface/PFJetMETcorrInputProducerT.h
  for(View<Jet>::const_iterator Jit = jets->begin(); Jit != jets->end(); ++Jit, ++iraw) {
    const PFJet* rawJet = static_cast<const PFJet*>( &((*Jit)) );
    reco::Candidate::LorentzVector rawJetP4 = rawJet->p4();
    std::vector<reco::PFCandidatePtr> cands = rawJet->getPFConstituents();
    for ( std::vector<reco::PFCandidatePtr>::const_iterator cand = cands.begin();cand != cands.end(); ++cand ) {
      if ( (*cand)->muonRef().isNonnull() && 
           ((*cand)->muonRef()->isGlobalMuon() || (*cand)->muonRef()->isStandAloneMuon() ) ) { // Loose muon selection
        rawJetP4 -= (*cand)->p4();
      }
    }
    if ( rawJetP4.pt() < 3 ) continue; // Discard very low momentum jets (beware: JECs can be higher than 200%...)
    fTJMetCorrNoMuPt ->push_back( rawJetP4.pt() );
    fTJMetCorrPhi    ->push_back( rawJetP4.phi() );
    fTJMetCorrRawEta ->push_back( rawJet->eta() ); // Official tool uses original jet for corrections...
    fTJMetCorrRawPt  ->push_back( rawJet->pt() );
    fTJMetCorrEMF    ->push_back( rawJet->chargedEmEnergyFraction() + rawJet->neutralEmEnergyFraction() );
    fTJMetCorrArea   ->push_back( rawJet->jetArea() );
  }
  }

  ////////////////////////////////////////////////////////////////////////////////
  // Special stuff for Model Scans ///////////////////////////////////////////////
  *fTxSMS=-1;
  *fTxbarSMS=-1;

  if(!fIsRealData&&fIsModelScan && !fIsFastSim) {
    
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
        *fTM0      = m0;
        *fTM12     = m12;
        *fTsignMu  = signMu;
        *fTA0      = A0;
        *fTtanBeta = tanb;
      }
    }
    *fTMassGlu = mGL;
    *fTMassChi = xCHI;
    *fTMassLSP = mLSP;
  }


  ////////////////////////////////////////////////////////////////////////////////
  // Full generator information ///////////////////////////////////////////////
  bool blabalot=false;
  if(!fIsRealData && !doPhotonStuff) {
    static const size_t maxNGenLocal = 2000;
    int genNIndex[maxNGenLocal];
    int genID[maxNGenLocal];
    float genPt[maxNGenLocal];
    float genEta[maxNGenLocal];
    float genPhi[maxNGenLocal];
    float genM[maxNGenLocal];
    int genMo1Index[maxNGenLocal];
    int genMo2Index[maxNGenLocal];
    int genNMo[maxNGenLocal];
    int genStatus[maxNGenLocal];
    int Promptness[maxNGenLocal];
    bool StoreFlag[maxNGenLocal];

    int nGenParticles=0;
    
    Handle<GenParticleCollection> genParticles;
    iEvent.getByLabel("genParticles", genParticles);
	  
    // STEP 1: Loop over all particles and store the information.
    for(size_t i = 0; i < genParticles->size(); ++ i) {
      if(i>=maxNGenLocal) {
        edm::LogWarning("NTP") << "@SUB=analyze()"
                               << "Maximum number of gen particles for local array exceeded";
        break;
      }
      nGenParticles++;
      const GenParticle & p = (*genParticles)[i];
      genNIndex[i]=0;
      genID[i]=p.pdgId();
      genPt[i]=p.pt();
      genPhi[i]=p.phi();
      genEta[i]=p.eta();
      genM[i]=p.mass();
      genStatus[i]=p.status();
      genMo1Index[i]=-1;
      genMo2Index[i]=-1;
      genNMo[i]=p.numberOfMothers();
      StoreFlag[i]=false;
      if(genID[i]==2212 && p.numberOfMothers()==0) Promptness[i]=0;
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
            if( gmom == NULL || ((testgm!=NULL) && (gmom->pt() == testgm->pt())) ) {
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
    *fTnGenParticles =  -1;
    for(int i=0;i<nGenParticles;i++) {
      if(!StoreFlag[i]) continue;
      (*fTnGenParticles)++;
      if( *fTnGenParticles > gMaxNGenParticles ) {
        edm::LogWarning("NTP") << "@SUB=analyze()"
                               << "Maximum number of gen particles exceeded";
        *fTMaxGenPartExceed = 1;
        break;
      }
	    
      genNIndex[i] = *fTnGenParticles;
	    
      //store everything
      fTgenInfoId ->push_back( genID[i] );
      fTgenInfoStatus ->push_back( genStatus[i] );
      fTgenInfoNMo ->push_back( genNMo[i] );
      if(genMo1Index[i]>=0) {
        fTgenInfoMo1->push_back( genNIndex[genMo1Index[i]] );
      } else {
        fTgenInfoMo1->push_back( -1 );
      }
      if(genMo2Index[i]>=0) {
        fTgenInfoMo2->push_back( genNIndex[genMo2Index[i]] );
      } else {
        fTgenInfoMo2->push_back( -1 );
      }
          
      fTPromptnessLevel->push_back( Promptness[i] );
      fTgenInfoPt ->push_back( genPt[i] );
      fTgenInfoEta ->push_back( genEta[i] );
      fTgenInfoPhi ->push_back( genPhi[i] );
      fTgenInfoM ->push_back( genM[i] );

      if(blabalot) {
        cout << "Working particle " << i << " with Promptness: " << Promptness[i] << "  (storage number: " << *fTnGenParticles << ")" << endl;
        cout << "    Mother is Particle " << genMo1Index[i] << " with Promptness " << Promptness[genMo1Index[i]] << endl;
        cout << "Particle " << *fTnGenParticles << "  (" << i << "): The particle has ID = " << genID[i] 
             << " and its mother has index " << genNIndex[genMo1Index[i]]     << " mo pt : " << genPt[genMo1Index[i]] << endl;
        cout << "stored:  " << *fTnGenParticles << "  (" << i << ")                      " << (*fTgenInfoId)[*fTnGenParticles] 
             << "                          " << (*fTgenInfoMo1)[*fTnGenParticles] << endl;
      }
    }
	  
    (*fTnGenParticles)++;
    if(blabalot) cout << "A total of " << *fTnGenParticles << " Particles  have been stored out of " << nGenParticles << " ( " << 100*(*fTnGenParticles)/(float)nGenParticles << " %)" << endl;
	  
	  
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
      *fTxSMS = (chi2mass/nchi1mass - chi1mass/nchi1mass) / (gluinomass/ngluinomass - chi1mass/nchi1mass);
      *fTxbarSMS = 1 - *fTxSMS; // Mariarosaria's definition of x
    }
  }// end of bloat with gen information


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
  for ( std::vector<PFFiller*>::iterator it = pfFillers.begin(); 
        it != pfFillers.end(); ++it ) 
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
  for ( std::vector<PFFiller*>::iterator it = pfFillers.begin(); 
        it != pfFillers.end(); ++it ) 
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
  produces<float,edm::InRun>("MinSCrawPt"    );
  produces<float,edm::InRun>("MaxPfCandEta"  );
  produces<float,edm::InRun>("MinEBRechitE"  );
  
  produces<float,edm::InRun>("MinGenLeptPt"  );
  produces<float,edm::InRun>("MaxGenLeptEta" );
  produces<float,edm::InRun>("MinGenPhotPt"  );
  produces<float,edm::InRun>("MaxGenPhotEta" );
  produces<float,edm::InRun>("MinGenJetPt"   );
  produces<float,edm::InRun>("MaxGenJetEta"  );

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
  produces<int,edm::InRun>("MaxNConv"    );
  produces<int,edm::InRun>("MaxNPfCand"    );
  produces<int,edm::InRun>("MaxNXtals"    );

  produces<std::vector<std::string>,edm::InRun>("HLTNames");
  produces<std::vector<std::string>,edm::InRun>("L1PhysMenu");
  produces<std::vector<std::string>,edm::InRun>("HLTLabels");

  produces<std::vector<std::string>,edm::InRun>("PileUpData");
  produces<std::vector<std::string>,edm::InRun>("PileUpMC");

  // Event products
  produces<int>("Run");
  produces<unsigned int>("Event");
  produces<int>("LumiSection");
  produces<float>("PtHat");
  produces<float>("QCDPartonicHT");
  produces<std::vector<int> >("LHEEventID");
  produces<std::vector<int> >("LHEEventStatus");
  produces<std::vector<int> >("LHEEventMotherFirst");
  produces<std::vector<int> >("LHEEventMotherSecond");
  produces<std::vector<float> >("LHEEventPx");
  produces<std::vector<float> >("LHEEventPy");
  produces<std::vector<float> >("LHEEventPz");
  produces<std::vector<float> >("LHEEventE");
  produces<std::vector<float> >("LHEEventM");
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
  produces<float>("RhoForIso");
  produces<float>("Weight");
  produces<std::vector<int> >("HLTResults");
  produces<std::vector<int> >("HLTPrescale");
  produces<std::vector<int> >("L1PhysResults");
  produces<std::vector<int> >("L1TechResults");
  produces<std::vector<int> >("NHLTObjs");
  for ( size_t i=0; i<gMaxHltNPaths; ++i ) {
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

  produces<std::vector<int> >("PhoVrtxListStart");
  produces<std::vector<int> >("JVrtxListStart");

  produces<int>("MaxGenPartExceed");
  produces<int>("nGenParticles");
  produces<std::vector<int> >("genInfoId");
  produces<std::vector<int> >("genInfoStatus");
  produces<std::vector<int> >("genInfoNMo");
  produces<std::vector<int> >("genInfoMo1");
  produces<std::vector<int> >("genInfoMo2");
  produces<std::vector<int> >("PromptnessLevel");
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
  produces<int>("CSCTightHaloID");
  produces<float>("PFType1MET");
  produces<float>("PFType1METpx");
  produces<float>("PFType1METpy");
  produces<float>("PFType1METphi");
  produces<float>("PFType1METSignificance");
  produces<float>("PFType1SumEt");
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
  produces<std::vector<int> >("MuIsPFMuon");
  produces<std::vector<int> >("MuIsStandaloneMuon");
  produces<std::vector<float> >("MuPx");
  produces<std::vector<float> >("MuPy");
  produces<std::vector<float> >("MuPz");
  produces<std::vector<float> >("MuPt");
  produces<std::vector<float> >("MuInnerTkPt");
  produces<std::vector<float> >("MuTkPtE");
  produces<std::vector<float> >("MuTkD0E");
  produces<std::vector<float> >("MuTkDzE");
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
  produces<std::vector<float> >("MuPfIsoR03ChHad");
  produces<std::vector<float> >("MuPfIsoR03NeHad");
  produces<std::vector<float> >("MuPfIsoR03Photon");
  produces<std::vector<float> >("MuPfIsoR03NeHadHighThresh");
  produces<std::vector<float> >("MuPfIsoR03PhotonHighThresh");
  produces<std::vector<float> >("MuPfIsoR03SumPUPt");
  produces<std::vector<float> >("MuPfIsoR04ChHad");
  produces<std::vector<float> >("MuPfIsoR04NeHad");
  produces<std::vector<float> >("MuPfIsoR04Photon");
  produces<std::vector<float> >("MuPfIsoR04NeHadHighThresh");
  produces<std::vector<float> >("MuPfIsoR04PhotonHighThresh");
  produces<std::vector<float> >("MuPfIsoR04SumPUPt");
  for ( std::vector<edm::InputTag>::const_iterator it = fMuonPfIsoTagsCustom.begin();
        it != fMuonPfIsoTagsCustom.end(); ++it ) {
    produces<std::vector<float> >(("Mu"+(*it).label()).c_str());
  }
  produces<std::vector<float> >("MuEem");
  produces<std::vector<float> >("MuEhad");
  produces<std::vector<float> >("MuD0BS");
  produces<std::vector<float> >("MuD0PV");
  produces<std::vector<float> >("MuD03DPV");
  produces<std::vector<float> >("MuD03DE");
  produces<std::vector<float> >("MuD0E");
  produces<std::vector<float> >("MuDzBS");
  produces<std::vector<float> >("MuDzPV");
  produces<std::vector<float> >("MuDzE");
  produces<std::vector<float> >("MuNChi2");
  produces<std::vector<int> >("MuNGlHits");
  produces<std::vector<int> >("MuNGlMuHits");
  produces<std::vector<int> >("MuNMuHits");
  produces<std::vector<int> >("MuNTkHits");
  produces<std::vector<int> >("MuNPxHits");
  produces<std::vector<float> >("MuInnerTkNChi2");
  produces<std::vector<int> >("MuNSiLayers");
  produces<std::vector<int> >("MuNMatches");
  produces<std::vector<int> >("MuNMatchedStations");
  produces<std::vector<int> >("MuNChambers");
  produces<std::vector<float> >("MuIsoMVA");
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
  produces<std::vector<int> >("ElIsEB");
  produces<std::vector<int> >("ElIsEE");
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
  produces<std::vector<float> >("ElD03DPV");
  produces<std::vector<float> >("ElD03DE");
  produces<std::vector<float> >("ElDzBS");
  produces<std::vector<float> >("ElDzPV");
  produces<std::vector<float> >("ElDzE");
  produces<std::vector<float> >("ElRelIso03");
  produces<std::vector<float> >("ElRelIso04");
  produces<std::vector<float> >("ElPfIsoChHad03");
  produces<std::vector<float> >("ElPfIsoNeHad03");
  produces<std::vector<float> >("ElPfIsoPhoton03");
  produces<std::vector<float> >("ElDR03TkSumPt");
  produces<std::vector<float> >("ElDR04TkSumPt");
  produces<std::vector<float> >("ElDR03EcalRecHitSumEt");
  produces<std::vector<float> >("ElDR04EcalRecHitSumEt");
  produces<std::vector<float> >("ElDR03HcalTowerSumEt");
  produces<std::vector<float> >("ElDR04HcalTowerSumEt");
  for ( std::vector<edm::InputTag>::const_iterator it = fElePfIsoTagsCustom.begin();
        it != fElePfIsoTagsCustom.end(); ++it ) {
    produces<std::vector<float> >(("El"+(*it).label()).c_str());
  }
  for ( std::vector<edm::InputTag>::const_iterator it = fElePfIsoTagsEvent.begin();
        it != fElePfIsoTagsEvent.end(); ++it ) {
    produces<std::vector<float> >(("ElEvent"+(*it).label()).c_str());
  }
  produces<std::vector<float> >("ElNChi2");
  produces<std::vector<float> >("ElKfTrkchi2");
  produces<std::vector<float> >("ElKfTrkhits");
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
  produces<std::vector<float> >("ElIDMVATrig");
  produces<std::vector<float> >("ElIDMVANoTrig");
  produces<std::vector<int> >("ElInGap");
  produces<std::vector<int> >("ElEcalDriven");
  produces<std::vector<int> >("ElTrackerDriven");
  produces<std::vector<int> >("ElBasicClustersSize");
  produces<std::vector<float> >("Elfbrem");
  produces<std::vector<float> >("ElHcalOverEcal");
  produces<std::vector<float> >("ElHcalOverEcalBc");
  produces<std::vector<float> >("ElE1x5");
  produces<std::vector<float> >("ElE5x5");
  produces<std::vector<float> >("ElE2x5Max");
  produces<std::vector<float> >("ElR9");
  produces<std::vector<float> >("ElSigmaIetaIeta");
  produces<std::vector<float> >("ElSigmaIphiIphi");
  produces<std::vector<float> >("ElScEtaWidth");
  produces<std::vector<float> >("ElScPhiWidth");
  produces<std::vector<float> >("ElDeltaPhiSeedClusterAtCalo");
  produces<std::vector<float> >("ElDeltaEtaSeedClusterAtCalo");
  produces<std::vector<float> >("ElDeltaPhiSuperClusterAtVtx");
  produces<std::vector<float> >("ElDeltaEtaSuperClusterAtVtx");
  produces<std::vector<float> >("ElCaloEnergy");
  produces<std::vector<float> >("ElTrkMomAtVtx");
  produces<std::vector<float> >("ElESuperClusterOverP");
  produces<std::vector<float> >("ElPreShowerOverRaw");
  produces<std::vector<float> >("ElEoPout");
  produces<std::vector<float> >("ElIoEmIoP");
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
  produces<std::vector<float> >("PhoHoverE2012");
  produces<std::vector<float> >("PhoSigmaIetaIeta");
  produces<std::vector<float> >("PhoSigmaIetaIphi");
  produces<std::vector<float> >("PhoSigmaIphiIphi");
  produces<std::vector<float> >("PhoS4Ratio");
  produces<std::vector<float> >("PhoLambdaRatio");
  produces<std::vector<float> >("PhoSCRawEnergy");
  produces<std::vector<float> >("PhoSCEtaWidth");
  produces<std::vector<float> >("PhoSCSigmaPhiPhi");
  produces<std::vector<int> >("PhoHasPixSeed");
  produces<std::vector<int> >("PhoHasConvTrks");
  produces<std::vector<int> >("PhoScSeedSeverity");
  produces<std::vector<float> >("PhoE1OverE9");
  produces<std::vector<float> >("PhoS4OverS1");
  produces<std::vector<float> >("PhoSigmaEtaEta");
  produces<std::vector<float> >("PhoSigmaRR");
  produces<std::vector<float> >("PhoNewIsoPFCharged");
  produces<std::vector<float> >("PhoNewIsoPFPhoton");
  produces<std::vector<float> >("PhoNewIsoPFNeutral");
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
//  produces<std::vector<float> >("PhoCone04PhotonIsodR0dEta0pt0");
//  produces<std::vector<float> >("PhoCone04PhotonIsodR0dEta0pt5");
//  produces<std::vector<float> >("PhoCone04PhotonIsodR8dEta0pt0");
//  produces<std::vector<float> >("PhoCone04PhotonIsodR8dEta0pt5");
//  produces<std::vector<float> >("PhoCone01PhotonIsodR045EB070EEdEta015pt08EB1EEmvVtx");
//  produces<std::vector<float> >("PhoCone02PhotonIsodR045EB070EEdEta015pt08EB1EEmvVtx");
//  produces<std::vector<float> >("PhoCone03PhotonIsodR045EB070EEdEta015pt08EB1EEmvVtx");
//  produces<std::vector<float> >("PhoCone04PhotonIsodR045EB070EEdEta015pt08EB1EEmvVtx");
//  produces<std::vector<float> >("PhoCone04NeutralHadronIsodR0dEta0pt0");
//  produces<std::vector<float> >("PhoCone04NeutralHadronIsodR0dEta0pt5");
//  produces<std::vector<float> >("PhoCone04NeutralHadronIsodR0dEta0pt0nocracks");
//  produces<std::vector<float> >("PhoCone04NeutralHadronIsodR0dEta0pt5nocracks");
//  produces<std::vector<float> >("PhoCone04NeutralHadronIsodR7dEta0pt0");
//  produces<std::vector<float> >("PhoCone04NeutralHadronIsodR7dEta0pt5");
//  produces<std::vector<float> >("PhoCone01NeutralHadronIsodR0dEta0pt0mvVtx");
//  produces<std::vector<float> >("PhoCone02NeutralHadronIsodR0dEta0pt0mvVtx");
//  produces<std::vector<float> >("PhoCone03NeutralHadronIsodR0dEta0pt0mvVtx");
//  produces<std::vector<float> >("PhoCone04NeutralHadronIsodR0dEta0pt0mvVtx");
//  produces<std::vector<float> >("PhoCone04ChargedHadronIsodR0dEta0pt0dz0old");
//  produces<std::vector<float> >("PhoCone04ChargedHadronIsodR0dEta0pt0PFnoPUold");
//  produces<std::vector<float> >("PhoCone04ChargedHadronIsodR015dEta0pt0dz0old");
//  produces<std::vector<float> >("PhoCone04ChargedHadronIsodR015dEta0pt0PFnoPUold");
//  produces<std::vector<float> >("PhoCone01ChargedHadronIsodR0dEta0pt0dz0");
//  produces<std::vector<float> >("PhoCone01ChargedHadronIsodR0dEta0pt0dz1dxy01");
//  produces<std::vector<float> >("PhoCone01ChargedHadronIsodR0dEta0pt0PFnoPU");
//  produces<std::vector<float> >("PhoCone01ChargedHadronIsodR015dEta0pt0dz0");
//  produces<std::vector<float> >("PhoCone01ChargedHadronIsodR015dEta0pt0dz1dxy01");
//  produces<std::vector<float> >("PhoCone01ChargedHadronIsodR015dEta0pt0PFnoPU");
//  produces<std::vector<float> >("PhoCone02ChargedHadronIsodR0dEta0pt0dz0");
//  produces<std::vector<float> >("PhoCone02ChargedHadronIsodR0dEta0pt0dz1dxy01");
//  produces<std::vector<float> >("PhoCone02ChargedHadronIsodR0dEta0pt0PFnoPU");
//  produces<std::vector<float> >("PhoCone02ChargedHadronIsodR015dEta0pt0dz0");
//  produces<std::vector<float> >("PhoCone02ChargedHadronIsodR015dEta0pt0dz1dxy01");
//  produces<std::vector<float> >("PhoCone02ChargedHadronIsodR015dEta0pt0PFnoPU");
//  produces<std::vector<float> >("PhoCone03ChargedHadronIsodR0dEta0pt0dz0");
//  produces<std::vector<float> >("PhoCone03ChargedHadronIsodR0dEta0pt0dz1dxy01");
//  produces<std::vector<float> >("PhoCone03ChargedHadronIsodR0dEta0pt0PFnoPU");
//  produces<std::vector<float> >("PhoCone03ChargedHadronIsodR015dEta0pt0dz0");
//  produces<std::vector<float> >("PhoCone03ChargedHadronIsodR015dEta0pt0dz1dxy01");
//  produces<std::vector<float> >("PhoCone03ChargedHadronIsodR015dEta0pt0PFnoPU");
//  produces<std::vector<float> >("PhoCone04ChargedHadronIsodR0dEta0pt0dz0");
//  produces<std::vector<float> >("PhoCone04ChargedHadronIsodR0dEta0pt0dz1dxy01");
//  produces<std::vector<float> >("PhoCone04ChargedHadronIsodR0dEta0pt0PFnoPU");
//  produces<std::vector<float> >("PhoCone04ChargedHadronIsodR015dEta0pt0dz0");
//  produces<std::vector<float> >("PhoCone04ChargedHadronIsodR015dEta0pt0dz1dxy01");
//  produces<std::vector<float> >("PhoCone04ChargedHadronIsodR015dEta0pt0PFnoPU");
  produces<std::vector<float> >("PhoCiCPFIsoChargedDR03");
  produces<std::vector<float> >("PhoCiCPFIsoNeutralDR03");
  produces<std::vector<float> >("PhoCiCPFIsoPhotonDR03");
  produces<std::vector<float> >("PhoCiCPFIsoChargedDR04");
  produces<std::vector<float> >("PhoCiCPFIsoNeutralDR04");
  produces<std::vector<float> >("PhoCiCPFIsoPhotonDR04");
  produces<std::vector<float> >("PhoSCEta");
  produces<std::vector<float> >("PhoSCPhiWidth");
  produces<std::vector<float> >("PhoIDMVA");
  produces<std::vector<bool> > ("PhoConvValidVtx");
  produces<std::vector<bool> > ("ElPassConversionVeto");
  produces<std::vector<bool> > ("PhoPassConversionVeto");
  produces<std::vector<float> >("PhoHCalIso2012ConeDR03");
  produces<std::vector<int> >  ("PhoConvNtracks");
  produces<std::vector<float> >("PhoConvChi2Probability");
  produces<std::vector<float> >("PhoConvEoverP");
  produces<int>("Nconv");
  produces<std::vector<bool> >("ConvValidVtx");
  produces<std::vector<int> >("ConvNtracks");
  produces<std::vector<float> >("ConvChi2Probability");
  produces<std::vector<float> >("ConvEoverP");
  produces<std::vector<float> >("ConvZofPrimVtxFromTrks");
  produces<int>("Ngv");
  produces<std::vector<float> >("gvSumPtHi");
  produces<std::vector<float> >("gvSumPtLo");
  produces<std::vector<int> >("gvNTkHi");
  produces<std::vector<int> >("gvNTkLo");
  produces<int>("NGoodSuperClusters");
  produces<std::vector<float> >("GoodSCEnergy");
  produces<std::vector<float> >("GoodSCEta");
  produces<std::vector<float> >("GoodSCPhi");
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
  for ( std::vector<edm::InputTag>::const_iterator it = fBtagTags.begin();
	it != fBtagTags.end(); ++it ) {
    produces<std::vector<float> >(("J"+(*it).label()).c_str());
  }
  produces<std::vector<int>   >("JPartonFlavour");
  produces<std::vector<float> >("JMass");
  produces<std::vector<float> >("JBetaStar");
  produces<std::vector<float> >("JBeta");
  produces<std::vector<float> >("JBetaSq");
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
  produces<std::vector<float> >("JMetCorrRawEta"); 
  produces<std::vector<float> >("JMetCorrPhi");  
  produces<std::vector<float> >("JMetCorrNoMuPt");  
  produces<std::vector<float> >("JMetCorrRawPt");  
  produces<std::vector<float> >("JMetCorrEMF"); 
  produces<std::vector<float> >("JMetCorrArea");
  for ( size_t i=0; i<gMaxNPileupJetIDAlgos; ++i ) {
    std::ostringstream s;
    s << i;
    produces<std::vector<bool> >(("JPassPileupIDL"+s.str()).c_str());
    produces<std::vector<bool> >(("JPassPileupIDM"+s.str()).c_str());
    produces<std::vector<bool> >(("JPassPileupIDT"+s.str()).c_str());
  }
  produces<std::vector<float> >("JQGTagLD");
  produces<std::vector<float> >("JQGTagMLP");
  produces<std::vector<float> >("JSmearedQGL");
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
  produces<float>("METR12");
  produces<float>("METR21");

produces<float>("Sigma");
produces<std::vector<float> >("GenPhotonIsoDR03");
produces<std::vector<float> >("GenPhotonIsoDR04");
produces<std::vector<float> >("SCX");
produces<std::vector<float> >("SCY");
produces<std::vector<float> >("SCZ");
produces<std::vector<int> >("SCXtalListStart");
produces<std::vector<int> >("SCNXtals");
produces<int>("NXtals");
produces<std::vector<float> >("XtalX");
produces<std::vector<float> >("XtalY");
produces<std::vector<float> >("XtalZ");
produces<std::vector<float> >("XtalEtaWidth");
produces<std::vector<float> >("XtalPhiWidth");
produces<std::vector<float> >("XtalFront1X");
produces<std::vector<float> >("XtalFront1Y");
produces<std::vector<float> >("XtalFront1Z");
produces<std::vector<float> >("XtalFront2X");
produces<std::vector<float> >("XtalFront2Y");
produces<std::vector<float> >("XtalFront2Z");
produces<std::vector<float> >("XtalFront3X");
produces<std::vector<float> >("XtalFront3Y");
produces<std::vector<float> >("XtalFront3Z");
produces<std::vector<float> >("XtalFront4X");
produces<std::vector<float> >("XtalFront4Y");
produces<std::vector<float> >("XtalFront4Z");
produces<int>("NPfCand");
produces<std::vector<int> >("PfCandPdgId");
produces<std::vector<float> >("PfCandEta");
produces<std::vector<float> >("PfCandPhi");
produces<std::vector<float> >("PfCandEnergy");
produces<std::vector<float> >("PfCandEcalEnergy");
produces<std::vector<float> >("PfCandPt");
produces<std::vector<float> >("PfCandVx");
produces<std::vector<float> >("PfCandVy");
produces<std::vector<float> >("PfCandVz");
produces<std::vector<int> >("PfCandBelongsToJet");
//produces<std::vector<int> >("PfCandHasHitInFirstPixelLayer");
//produces<std::vector<float> >("PfCandTrackRefPx");
//produces<std::vector<float> >("PfCandTrackRefPy");
//produces<std::vector<float> >("PfCandTrackRefPz");
produces<std::vector<int> >("PhoMatchedPFPhotonOrElectronCand");
produces<std::vector<int> >("PhoFootprintPfCandsListStart");
produces<std::vector<int> >("PhoFootprintPfCands");
produces<std::vector<float> >("PhoVx");
produces<std::vector<float> >("PhoVy");
produces<std::vector<float> >("PhoVz");
produces<std::vector<float> >("PhoRegrEnergy");
produces<std::vector<float> >("PhoRegrEnergyErr");
//produces<std::vector<float> >("PhoCone01PhotonIsodEta015EBdR070EEmvVtx");
//produces<std::vector<float> >("PhoCone02PhotonIsodEta015EBdR070EEmvVtx");
//produces<std::vector<float> >("PhoCone03PhotonIsodEta015EBdR070EEmvVtx");
//produces<std::vector<float> >("PhoCone04PhotonIsodEta015EBdR070EEmvVtx");
//produces<std::vector<float> >("PhoCone01NeutralHadronIsomvVtx");
//produces<std::vector<float> >("PhoCone02NeutralHadronIsomvVtx");
//produces<std::vector<float> >("PhoCone03NeutralHadronIsomvVtx");
//produces<std::vector<float> >("PhoCone04NeutralHadronIsomvVtx");
//produces<std::vector<float> >("PhoCone01ChargedHadronIsodR02dz02dxy01");
//produces<std::vector<float> >("PhoCone02ChargedHadronIsodR02dz02dxy01");
//produces<std::vector<float> >("PhoCone03ChargedHadronIsodR02dz02dxy01");
//produces<std::vector<float> >("PhoCone04ChargedHadronIsodR02dz02dxy01");
//produces<std::vector<float> >("PhoCone03PFCombinedIso");
//produces<std::vector<float> >("PhoCone04PFCombinedIso");
produces<std::vector<int> >("Diphotonsfirst");
produces<std::vector<int> >("Diphotonssecond");
produces<std::vector<int> >("Vtxdiphoh2gglobe");
produces<std::vector<int> >("Vtxdiphomva");
produces<std::vector<int> >("Vtxdiphoproductrank");
produces<std::vector<float> >("PhoSCRemovalPFIsoCharged");
produces<std::vector<float> >("PhoSCRemovalPFIsoChargedPrimVtx");
produces<std::vector<float> >("PhoSCRemovalPFIsoNeutral");
produces<std::vector<float> >("PhoSCRemovalPFIsoPhoton");
produces<std::vector<float> >("PhoSCRemovalPFIsoChargedRCone");
produces<std::vector<float> >("PhoSCRemovalPFIsoChargedPrimVtxRCone");
produces<std::vector<float> >("PhoSCRemovalPFIsoNeutralRCone");
produces<std::vector<float> >("PhoSCRemovalPFIsoPhotonRCone");
produces<std::vector<float> >("PhoSCRemovalRConeEta");
produces<std::vector<float> >("PhoSCRemovalRConePhi");


}

//________________________________________________________________________________________
// Reset all event variables
void NTupleProducer::resetProducts( void ) {
  
  fTRun.reset(new int(-999));
  fTEvent.reset(new unsigned int(0));
  fTLumiSection.reset(new int(-999));
  fTPtHat.reset(new float(-999.99));
  fTQCDPartonicHT.reset(new float(-999.99));
  fTLHEEventID.reset(new std::vector<int> );
  fTLHEEventStatus.reset(new std::vector<int> );
  fTLHEEventMotherFirst.reset(new std::vector<int> );
  fTLHEEventMotherSecond.reset(new std::vector<int> );
  fTLHEEventPx.reset(new std::vector<float> );
  fTLHEEventPy.reset(new std::vector<float> );
  fTLHEEventPz.reset(new std::vector<float> );
  fTLHEEventE .reset(new std::vector<float> );
  fTLHEEventM .reset(new std::vector<float> );
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
  fTpdfWsum.reset(new float(0.0));
  fTNPdfs.reset(new int(1));
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
  fTRhoForIso.reset(new float(-999.99));
  fTWeight.reset(new float(-999.99));
  fTHLTResults.reset(new std::vector<int> );
  fTHLTPrescale.reset(new std::vector<int> );
  fTL1PhysResults.reset(new std::vector<int> );
  fTL1TechResults.reset(new std::vector<int> );
  fTNHLTObjs.reset(new std::vector<int> );
  for ( size_t i=0; i<gMaxHltNPaths; ++i ) {
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

  fTPhoVrtxListStart.reset(new std::vector<int>);
  fTJVrtxListStart.reset(new std::vector<int>);

  fTMaxGenPartExceed.reset(new int(-999));
  fTnGenParticles.reset(new int(0));
  fTgenInfoId.reset(new std::vector<int>);
  fTgenInfoStatus.reset(new std::vector<int>);
  fTgenInfoNMo.reset(new std::vector<int>);
  fTgenInfoMo1.reset(new std::vector<int>);
  fTgenInfoMo2.reset(new std::vector<int>);
  fTPromptnessLevel.reset(new std::vector<int>);
  fTgenInfoPt.reset(new std::vector<float>);
  fTgenInfoEta.reset(new std::vector<float>);
  fTgenInfoPhi.reset(new std::vector<float>);
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
  fTNCaloTowers.reset(new int(0));
  fTGoodEvent.reset(new int(0));
  fTMaxMuExceed.reset(new int(0));
  fTMaxElExceed.reset(new int(0));
  fTMaxJetExceed.reset(new int(0));
  fTMaxUncJetExceed.reset(new int(0));
  fTMaxTrkExceed.reset(new int(0));
  fTMaxPhotonsExceed.reset(new int(0));
  fTMaxGenLepExceed.reset(new int(0));
  fTMaxGenPhoExceed.reset(new int(0));
  fTMaxGenJetExceed.reset(new int(0));
  fTMaxVerticesExceed.reset(new int(0));
  fTCSCTightHaloID.reset(new int(-999));
  fTPFType1MET.reset(new float(-999.99));
  fTPFType1METpx.reset(new float(-999.99));
  fTPFType1METpy.reset(new float(-999.99));
  fTPFType1METphi.reset(new float(-999.99));
  fTPFType1METSignificance.reset(new float(-999.99));
  fTPFType1SumEt.reset(new float(-999.99));
  //FR fPBNRFlag.reset(new int(-999));
  fTNGenLeptons.reset(new int(0));
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
  fTNGenPhotons.reset(new int(0));
  fTGenPhotonPt.reset(new std::vector<float> );
  fTGenPhotonEta.reset(new std::vector<float> );
  fTGenPhotonPhi.reset(new std::vector<float> );
  fTGenPhotonPartonMindR.reset(new std::vector<float> );
  fTGenPhotonMotherID.reset(new std::vector<int> );
  fTGenPhotonMotherStatus.reset(new std::vector<int> );
  fTNGenJets.reset(new int(0));
  fTGenJetPt.reset(new std::vector<float> );
  fTGenJetEta.reset(new std::vector<float> );
  fTGenJetPhi.reset(new std::vector<float> );
  fTGenJetE.reset(new std::vector<float> );
  fTGenJetEmE.reset(new std::vector<float> );
  fTGenJetHadE.reset(new std::vector<float> );
  fTGenJetInvE.reset(new std::vector<float> );
  fTNVrtx.reset(new int(0));
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

  fTNMus.reset(new int(0));
  fTNMusTot.reset(new int(0));
  fTNGMus.reset(new int(0));
  fTNTMus.reset(new int(0));
  fTMuGood.reset(new std::vector<int> );
  fTMuIsIso.reset(new std::vector<int> );
  fTMuIsGlobalMuon.reset(new std::vector<int> );
  fTMuIsTrackerMuon.reset(new std::vector<int> );
  fTMuIsPFMuon.reset(new std::vector<int> );
  fTMuIsStandaloneMuon.reset(new std::vector<int> );
  fTMuPx.reset(new std::vector<float> );
  fTMuPy.reset(new std::vector<float> );
  fTMuPz.reset(new std::vector<float> );
  fTMuPt.reset(new std::vector<float> );
  fTMuInnerTkPt.reset(new std::vector<float> );
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
  fTMuPfIsoR03ChHad  .reset(new std::vector<float> );
  fTMuPfIsoR03NeHad  .reset(new std::vector<float> );
  fTMuPfIsoR03Photon .reset(new std::vector<float> );
  fTMuPfIsoR03NeHadHighThresh  .reset(new std::vector<float> );
  fTMuPfIsoR03PhotonHighThresh .reset(new std::vector<float> );
  fTMuPfIsoR03SumPUPt.reset(new std::vector<float> );
  fTMuPfIsoR04ChHad  .reset(new std::vector<float> );
  fTMuPfIsoR04NeHad  .reset(new std::vector<float> );
  fTMuPfIsoR04Photon .reset(new std::vector<float> );
  fTMuPfIsoR04NeHadHighThresh  .reset(new std::vector<float> );
  fTMuPfIsoR04PhotonHighThresh .reset(new std::vector<float> );
  fTMuPfIsoR04SumPUPt.reset(new std::vector<float> );
  size_t ipfisotag = 0;
  for ( std::vector<edm::InputTag>::const_iterator it = fMuonPfIsoTagsCustom.begin();
        it != fMuonPfIsoTagsCustom.end(); ++it ) {
    fTMuPfIsosCustom[ipfisotag++].reset(new std::vector<float> );
  }
  fTMuEem.reset(new std::vector<float> );
  fTMuEhad.reset(new std::vector<float> );
  fTMuD0BS.reset(new std::vector<float> );
  fTMuD0PV.reset(new std::vector<float> );
  fTMuD03DPV.reset(new std::vector<float> );
  fTMuD03DE.reset(new std::vector<float> );
  fTMuDzBS.reset(new std::vector<float> );
  fTMuDzPV.reset(new std::vector<float> );
  fTMuTkPtE.reset(new std::vector<float> );
  fTMuTkD0E.reset(new std::vector<float> );
  fTMuTkDzE.reset(new std::vector<float> );
  fTMuPtE.reset(new std::vector<float> );
  fTMuD0E.reset(new std::vector<float> );
  fTMuDzE.reset(new std::vector<float> );
  fTMuNChi2.reset(new std::vector<float> );
  fTMuNGlHits.reset(new std::vector<int> );
  fTMuNGlMuHits.reset(new std::vector<int> );
  fTMuNMuHits.reset(new std::vector<int> );
  fTMuNTkHits.reset(new std::vector<int> );
  fTMuNPxHits.reset(new std::vector<int> );
  fTMuInnerTkNChi2.reset(new std::vector<float> );
  fTMuNSiLayers.reset(new std::vector<int> );
  fTMuNMatches.reset(new std::vector<int> );
  fTMuNMatchedStations.reset(new std::vector<int> );
  fTMuNChambers.reset(new std::vector<int> );
  fTMuIsoMVA.reset(new std::vector<float> );
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
  fTNEBhits.reset(new int(0));
  fTEBrechitE.reset(new std::vector<float> );
  fTEBrechitPt.reset(new std::vector<float> );
  fTEBrechitEta.reset(new std::vector<float> );
  fTEBrechitPhi.reset(new std::vector<float> );
  fTEBrechitChi2.reset(new std::vector<float> );
  fTEBrechitTime.reset(new std::vector<float> );
  fTEBrechitE4oE1.reset(new std::vector<float> );
  fTEBrechitE2oE9.reset(new std::vector<float> );
  fTNEles.reset(new int(0));
  fTNElesTot.reset(new int(0));
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
  fTElIsEB.reset(new std::vector<int> );
  fTElIsEE.reset(new std::vector<int> );
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
  fTElD03DPV.reset(new std::vector<float> );
  fTElD03DE.reset(new std::vector<float> );
  fTElDzBS.reset(new std::vector<float> );
  fTElDzPV.reset(new std::vector<float> );
  fTElDzE.reset(new std::vector<float> );
  fTElRelIso03.reset(new std::vector<float> );
  fTElRelIso04.reset(new std::vector<float> );
  fTElPfIsoChHad03.reset(new std::vector<float> );
  fTElPfIsoNeHad03.reset(new std::vector<float> );
  fTElPfIsoPhoton03.reset(new std::vector<float> );
  fTElDR03TkSumPt.reset(new std::vector<float> );
  fTElDR04TkSumPt.reset(new std::vector<float> );
  fTElDR03EcalRecHitSumEt.reset(new std::vector<float> );
  fTElDR04EcalRecHitSumEt.reset(new std::vector<float> );
  fTElDR03HcalTowerSumEt.reset(new std::vector<float> );
  fTElDR04HcalTowerSumEt.reset(new std::vector<float> );
  ipfisotag = 0;
  for ( std::vector<edm::InputTag>::const_iterator it = fElePfIsoTagsCustom.begin();
        it != fElePfIsoTagsCustom.end(); ++it ) {
    fTElPfIsosCustom[ipfisotag++].reset(new std::vector<float> );
  }
  ipfisotag = 0;
  for ( std::vector<edm::InputTag>::const_iterator it = fElePfIsoTagsEvent.begin();
        it != fElePfIsoTagsEvent.end(); ++it ) {
    fTElPfIsosEvent[ipfisotag++].reset(new std::vector<float> );
  }
  fTElNChi2.reset(new std::vector<float> );
  fTElKfTrkchi2.reset(new std::vector<float> );
  fTElKfTrkhits.reset(new std::vector<float> );
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
  fTElIDMVATrig.reset(new std::vector<float> );
  fTElIDMVANoTrig.reset(new std::vector<float> );
  fTElInGap.reset(new std::vector<int> );
  fTElEcalDriven.reset(new std::vector<int> );
  fTElTrackerDriven.reset(new std::vector<int> );
  fTElBasicClustersSize.reset(new std::vector<int> );
  fTElfbrem.reset(new std::vector<float> );
  fTElHcalOverEcal.reset(new std::vector<float> );
  fTElHcalOverEcalBc.reset(new std::vector<float> );
  fTElE1x5.reset(new std::vector<float> );
  fTElE5x5.reset(new std::vector<float> );
  fTElE2x5Max.reset(new std::vector<float> );
  fTElR9.reset(new std::vector<float> );
  fTElSigmaIetaIeta.reset(new std::vector<float> );
  fTElSigmaIphiIphi.reset(new std::vector<float> );
  fTElScEtaWidth.reset(new std::vector<float> );
  fTElScPhiWidth.reset(new std::vector<float> );
  fTElDeltaPhiSeedClusterAtCalo.reset(new std::vector<float> );
  fTElDeltaEtaSeedClusterAtCalo.reset(new std::vector<float> );
  fTElDeltaPhiSuperClusterAtVtx.reset(new std::vector<float> );
  fTElDeltaEtaSuperClusterAtVtx.reset(new std::vector<float> );
  fTElCaloEnergy.reset(new std::vector<float> );
  fTElTrkMomAtVtx.reset(new std::vector<float> );
  fTElESuperClusterOverP.reset(new std::vector<float> );
  fTElIoEmIoP.reset(new std::vector<float> );
  fTElEoPout.reset(new std::vector<float> );
  fTElPreShowerOverRaw.reset(new std::vector<float> );
  fTElNumberOfMissingInnerHits.reset(new std::vector<int> );
  fTElSCindex.reset(new std::vector<int> );
  fTElConvPartnerTrkDist.reset(new std::vector<float> );
  fTElPassConversionVeto.reset(new std::vector<bool>);
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
  fTPhoPassConversionVeto.reset(new std::vector<bool>);
  fTNPhotons.reset(new int(0));
  fTNPhotonsTot.reset(new int(0));
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
  fTPhoHoverE2012.reset(new std::vector<float> );
  fTPhoSigmaIetaIeta.reset(new std::vector<float> );
  fTPhoSigmaIetaIphi.reset(new std::vector<float> );
  fTPhoSigmaIphiIphi.reset(new std::vector<float> );
  fTPhoS4Ratio.reset(new std::vector<float> );
  fTPhoLambdaRatio.reset(new std::vector<float> );
  fTPhoSCRawEnergy.reset(new std::vector<float> );
  fTPhoSCEtaWidth.reset(new std::vector<float> );
  fTPhoSCSigmaPhiPhi.reset(new std::vector<float> );
  fTPhoHasPixSeed.reset(new std::vector<int> );
  fTPhoHasConvTrks.reset(new std::vector<int> );
  fTPhoScSeedSeverity.reset(new std::vector<int> );
  fTPhoE1OverE9.reset(new std::vector<float> );
  fTPhoS4OverS1.reset(new std::vector<float> );
  fTPhoSigmaEtaEta.reset(new std::vector<float> );
  fTPhoSigmaRR.reset(new std::vector<float> );
  fTPhoHCalIso2012ConeDR03.reset(new std::vector<float> );
  fTPhoNewIsoPFCharged.reset(new std::vector<float> );
  fTPhoNewIsoPFPhoton.reset(new std::vector<float> );
  fTPhoNewIsoPFNeutral.reset(new std::vector<float> );
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
//  fTPhoCone04PhotonIsodR0dEta0pt0.reset(new std::vector<float> );
//  fTPhoCone04PhotonIsodR0dEta0pt5.reset(new std::vector<float> );
//  fTPhoCone04PhotonIsodR8dEta0pt0.reset(new std::vector<float> );
//  fTPhoCone04PhotonIsodR8dEta0pt5.reset(new std::vector<float> );
//  fTPhoCone01PhotonIsodR045EB070EEdEta015pt08EB1EEmvVtx.reset(new std::vector<float> );
//  fTPhoCone02PhotonIsodR045EB070EEdEta015pt08EB1EEmvVtx.reset(new std::vector<float> );
//  fTPhoCone03PhotonIsodR045EB070EEdEta015pt08EB1EEmvVtx.reset(new std::vector<float> );
//  fTPhoCone04PhotonIsodR045EB070EEdEta015pt08EB1EEmvVtx.reset(new std::vector<float> );
//  fTPhoCone04NeutralHadronIsodR0dEta0pt0.reset(new std::vector<float> );
//  fTPhoCone04NeutralHadronIsodR0dEta0pt5.reset(new std::vector<float> );
//  fTPhoCone04NeutralHadronIsodR0dEta0pt0nocracks.reset(new std::vector<float> );
//  fTPhoCone04NeutralHadronIsodR0dEta0pt5nocracks.reset(new std::vector<float> );
//  fTPhoCone04NeutralHadronIsodR7dEta0pt0.reset(new std::vector<float> );
//  fTPhoCone04NeutralHadronIsodR7dEta0pt5.reset(new std::vector<float> );
//  fTPhoCone01NeutralHadronIsodR0dEta0pt0mvVtx.reset(new std::vector<float> );
//  fTPhoCone02NeutralHadronIsodR0dEta0pt0mvVtx.reset(new std::vector<float> );
//  fTPhoCone03NeutralHadronIsodR0dEta0pt0mvVtx.reset(new std::vector<float> );
//  fTPhoCone04NeutralHadronIsodR0dEta0pt0mvVtx.reset(new std::vector<float> );
//  fTPhoCone04ChargedHadronIsodR0dEta0pt0dz0old.reset(new std::vector<float> );
//  fTPhoCone04ChargedHadronIsodR0dEta0pt0PFnoPUold.reset(new std::vector<float> );
//  fTPhoCone04ChargedHadronIsodR015dEta0pt0dz0old.reset(new std::vector<float> );
//  fTPhoCone04ChargedHadronIsodR015dEta0pt0PFnoPUold.reset(new std::vector<float> );
//  fTPhoCone01ChargedHadronIsodR0dEta0pt0dz0.reset(new std::vector<float> );
//  fTPhoCone01ChargedHadronIsodR0dEta0pt0dz1dxy01.reset(new std::vector<float> );
//  fTPhoCone01ChargedHadronIsodR0dEta0pt0PFnoPU.reset(new std::vector<float> );
//  fTPhoCone01ChargedHadronIsodR015dEta0pt0dz0.reset(new std::vector<float> );
//  fTPhoCone01ChargedHadronIsodR015dEta0pt0dz1dxy01.reset(new std::vector<float> );
//  fTPhoCone01ChargedHadronIsodR015dEta0pt0PFnoPU.reset(new std::vector<float> );
//  fTPhoCone02ChargedHadronIsodR0dEta0pt0dz0.reset(new std::vector<float> );
//  fTPhoCone02ChargedHadronIsodR0dEta0pt0dz1dxy01.reset(new std::vector<float> );
//  fTPhoCone02ChargedHadronIsodR0dEta0pt0PFnoPU.reset(new std::vector<float> );
//  fTPhoCone02ChargedHadronIsodR015dEta0pt0dz0.reset(new std::vector<float> );
//  fTPhoCone02ChargedHadronIsodR015dEta0pt0dz1dxy01.reset(new std::vector<float> );
//  fTPhoCone02ChargedHadronIsodR015dEta0pt0PFnoPU.reset(new std::vector<float> );
//  fTPhoCone03ChargedHadronIsodR0dEta0pt0dz0.reset(new std::vector<float> );
//  fTPhoCone03ChargedHadronIsodR0dEta0pt0dz1dxy01.reset(new std::vector<float> );
//  fTPhoCone03ChargedHadronIsodR0dEta0pt0PFnoPU.reset(new std::vector<float> );
//  fTPhoCone03ChargedHadronIsodR015dEta0pt0dz0.reset(new std::vector<float> );
//  fTPhoCone03ChargedHadronIsodR015dEta0pt0dz1dxy01.reset(new std::vector<float> );
//  fTPhoCone03ChargedHadronIsodR015dEta0pt0PFnoPU.reset(new std::vector<float> );
//  fTPhoCone04ChargedHadronIsodR0dEta0pt0dz0.reset(new std::vector<float> );
//  fTPhoCone04ChargedHadronIsodR0dEta0pt0dz1dxy01.reset(new std::vector<float> );
//  fTPhoCone04ChargedHadronIsodR0dEta0pt0PFnoPU.reset(new std::vector<float> );
//  fTPhoCone04ChargedHadronIsodR015dEta0pt0dz0.reset(new std::vector<float> );
//  fTPhoCone04ChargedHadronIsodR015dEta0pt0dz1dxy01.reset(new std::vector<float> );
//  fTPhoCone04ChargedHadronIsodR015dEta0pt0PFnoPU.reset(new std::vector<float> );
  fTPhoCiCPFIsoChargedDR03.reset(new std::vector<float> );
  fTPhoCiCPFIsoNeutralDR03.reset(new std::vector<float> );
  fTPhoCiCPFIsoPhotonDR03.reset(new std::vector<float> );
  fTPhoCiCPFIsoChargedDR04.reset(new std::vector<float> );
  fTPhoCiCPFIsoNeutralDR04.reset(new std::vector<float> );
  fTPhoCiCPFIsoPhotonDR04.reset(new std::vector<float> );
  fTPhoSCX.reset(new std::vector<float> );
  fTPhoSCY.reset(new std::vector<float> );
  fTPhoSCZ.reset(new std::vector<float> );
  fTPhoSCEta.reset(new std::vector<float> );
  fTPhoSCPhiWidth.reset(new std::vector<float> );
  fTPhoIDMVA.reset(new std::vector<float> );
  fTPhoConvValidVtx.reset(new std::vector<bool>);
  fTPhoConvNtracks.reset(new std::vector<int>);
  fTPhoConvChi2Probability.reset(new std::vector<float>);
  fTPhoConvEoverP.reset(new std::vector<float>);
  fTNconv.reset(new int(0));
  fTConvValidVtx.reset(new std::vector<bool>);
  fTConvNtracks.reset(new std::vector<int>);
  fTConvChi2Probability.reset(new std::vector<float>);
  fTConvEoverP.reset(new std::vector<float>);
  fTConvZofPrimVtxFromTrks.reset(new std::vector<float>);
  fTNgv.reset(new int(0));
  fTgvSumPtHi.reset(new std::vector<float>);
  fTgvSumPtLo.reset(new std::vector<float>);
  fTgvNTkHi.reset(new std::vector<int>);
  fTgvNTkLo.reset(new std::vector<int>);
  for (int i=0; i<gMaxNPhotons; i++) {
    pho_conv_vtx[i]=TVector3();
    pho_conv_refitted_momentum[i]=TVector3();
    conv_vtx[i]=TVector3();
    conv_refitted_momentum[i]=TVector3();
    conv_singleleg_momentum[i]=TVector3();
  }
  for (int i=0; i<gMaxNGenVtx; i++) {
    gv_pos[i]=TVector3();
    gv_p3[i]=TVector3();
  }
  fTNGoodSuperClusters.reset(new int(0));
  fTGoodSCEnergy.reset(new std::vector<float> );
  fTGoodSCEta.reset(new std::vector<float> );
  fTGoodSCPhi.reset(new std::vector<float> );
  fTNSuperClusters.reset(new int(0));
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
  fTNJets.reset(new int(0));
  fTNJetsTot.reset(new int(0));
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
  size_t ibtag = 0;
  for ( std::vector<edm::InputTag>::const_iterator it = fBtagTags.begin();
	it != fBtagTags.end(); ++it ) {
    fTJbTagProb[ibtag++].reset(new std::vector<float> );
  }
  fTJPartonFlavour.reset(new std::vector<int> );
  fTJMass.reset(new std::vector<float> );
  fTJBetaStar.reset(new std::vector<float> );
  fTJBeta.reset(new std::vector<float> );
  fTJBetaSq.reset(new std::vector<float> );
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
  fTJMetCorrRawEta.reset(new std::vector<float> );
  fTJMetCorrPhi.reset(new std::vector<float> );
  fTJMetCorrNoMuPt.reset(new std::vector<float> );
  fTJMetCorrRawPt.reset(new std::vector<float> );
  fTJMetCorrEMF.reset(new std::vector<float> );
  fTJMetCorrArea.reset(new std::vector<float> );
  for ( unsigned int i=0; i<gMaxNPileupJetIDAlgos; ++i ) {
    fTJPassPileupIDL[i].reset(new std::vector<bool> );
    fTJPassPileupIDM[i].reset(new std::vector<bool> );
    fTJPassPileupIDT[i].reset(new std::vector<bool> );
  }
  fTJQGTagLD.reset(new std::vector<float> );
  fTJQGTagMLP.reset(new std::vector<float> );
  fTJSmearedQGL.reset(new std::vector<float> );
  fTNTracks.reset(new int(0));
  fTNTracksTot.reset(new int(0));
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
  fTMETR12.reset(new float(-999.99));
  fTMETR21.reset(new float(-999.99));


fTSigma.reset(new float(-999.99) );
fTGenPhotonIsoDR03.reset(new std::vector<float>  );
fTGenPhotonIsoDR04.reset(new std::vector<float>  );
fTSCX.reset(new std::vector<float>  );
fTSCY.reset(new std::vector<float>  );
fTSCZ.reset(new std::vector<float>  );
fTSCXtalListStart.reset(new std::vector<int>  );
fTSCNXtals.reset(new std::vector<int>  );
fTNXtals.reset(new int(0) );
fTXtalX.reset(new std::vector<float> );
fTXtalY.reset(new std::vector<float> );
fTXtalZ.reset(new std::vector<float> );
fTXtalEtaWidth.reset(new std::vector<float> );
fTXtalPhiWidth.reset(new std::vector<float> );
fTXtalFront1X.reset(new std::vector<float>  );
fTXtalFront1Y.reset(new std::vector<float>  );
fTXtalFront1Z.reset(new std::vector<float>  );
fTXtalFront2X.reset(new std::vector<float>  );
fTXtalFront2Y.reset(new std::vector<float>  );
fTXtalFront2Z.reset(new std::vector<float>  );
fTXtalFront3X.reset(new std::vector<float>  );
fTXtalFront3Y.reset(new std::vector<float>  );
fTXtalFront3Z.reset(new std::vector<float>  );
fTXtalFront4X.reset(new std::vector<float>  );
fTXtalFront4Y.reset(new std::vector<float>  );
fTXtalFront4Z.reset(new std::vector<float>  );
fTNPfCand.reset(new int(0) );
fTPfCandPdgId.reset(new std::vector<int>  );
fTPfCandEta.reset(new std::vector<float>  );
fTPfCandPhi.reset(new std::vector<float>  );
fTPfCandEnergy.reset(new std::vector<float>  );
fTPfCandEcalEnergy.reset(new std::vector<float>  );
fTPfCandPt.reset(new std::vector<float>  );
fTPfCandVx.reset(new std::vector<float>  );
fTPfCandVy.reset(new std::vector<float>  );
fTPfCandVz.reset(new std::vector<float>  );
fTPfCandBelongsToJet.reset(new std::vector<int>  );
//fTPfCandHasHitInFirstPixelLayer.reset(new std::vector<int>  );
//fTPfCandTrackRefPx.reset(new std::vector<float>  );
//fTPfCandTrackRefPy.reset(new std::vector<float>  );
//fTPfCandTrackRefPz.reset(new std::vector<float>  );
fTPhoMatchedPFPhotonOrElectronCand.reset(new std::vector<int>  );
fTPhoFootprintPfCandsListStart.reset(new std::vector<int>  );
fTPhoFootprintPfCands.reset(new std::vector<int>  );
fTPhoVx.reset(new std::vector<float>  );
fTPhoVy.reset(new std::vector<float>  );
fTPhoVz.reset(new std::vector<float>  );
fTPhoRegrEnergy.reset(new std::vector<float>  );
fTPhoRegrEnergyErr.reset(new std::vector<float>  );
//fTPhoCone01PhotonIsodEta015EBdR070EEmvVtx.reset(new std::vector<float>  );
//fTPhoCone02PhotonIsodEta015EBdR070EEmvVtx.reset(new std::vector<float>  );
//fTPhoCone03PhotonIsodEta015EBdR070EEmvVtx.reset(new std::vector<float>  );
//fTPhoCone04PhotonIsodEta015EBdR070EEmvVtx.reset(new std::vector<float>  );
//fTPhoCone01NeutralHadronIsomvVtx.reset(new std::vector<float>  );
//fTPhoCone02NeutralHadronIsomvVtx.reset(new std::vector<float>  );
//fTPhoCone03NeutralHadronIsomvVtx.reset(new std::vector<float>  );
//fTPhoCone04NeutralHadronIsomvVtx.reset(new std::vector<float>  );
//fTPhoCone01ChargedHadronIsodR02dz02dxy01.reset(new std::vector<float>  );
//fTPhoCone02ChargedHadronIsodR02dz02dxy01.reset(new std::vector<float>  );
//fTPhoCone03ChargedHadronIsodR02dz02dxy01.reset(new std::vector<float>  );
//fTPhoCone04ChargedHadronIsodR02dz02dxy01.reset(new std::vector<float>  );
//fTPhoCone03PFCombinedIso.reset(new std::vector<float>  );
//fTPhoCone04PFCombinedIso.reset(new std::vector<float>  );
fTDiphotonsfirst.reset(new std::vector<int>  );
fTDiphotonssecond.reset(new std::vector<int>  );
fTVtxdiphoh2gglobe.reset(new std::vector<int>  );
fTVtxdiphomva.reset(new std::vector<int>  );
fTVtxdiphoproductrank.reset(new std::vector<int>  );
fTPhoSCRemovalPFIsoCharged.reset(new std::vector<float>  );
fTPhoSCRemovalPFIsoChargedPrimVtx.reset(new std::vector<float>  );
fTPhoSCRemovalPFIsoNeutral.reset(new std::vector<float>  );
fTPhoSCRemovalPFIsoPhoton.reset(new std::vector<float>  );
fTPhoSCRemovalPFIsoChargedRCone.reset(new std::vector<float>  );
fTPhoSCRemovalPFIsoChargedPrimVtxRCone.reset(new std::vector<float>  );
fTPhoSCRemovalPFIsoNeutralRCone.reset(new std::vector<float>  );
fTPhoSCRemovalPFIsoPhotonRCone.reset(new std::vector<float>  );
fTPhoSCRemovalRConeEta.reset(new std::vector<float>  );
fTPhoSCRemovalRConePhi.reset(new std::vector<float>  );

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
  fRMinSCrawPt  .reset(new float(-999.99));
  fRMaxPfCandEta.reset(new float(-999.99));                                         
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
  fRMaxNConv    .reset(new int(-999));
  fRMaxNPfCand  .reset(new int(-999));
  fRMaxNXtals .reset(new int(-999));

  fRL1PhysMenu  .reset(new std::vector<std::string>);

}

//________________________________________________________________________________________
void NTupleProducer::putProducts( edm::Event& event ) {
  
  event.put(fTRun,   "Run");
  event.put(fTEvent, "Event");
  event.put(fTLumiSection, "LumiSection");
  event.put(fTPtHat, "PtHat");
  event.put(fTQCDPartonicHT, "QCDPartonicHT");
  event.put(fTLHEEventID, "LHEEventID");
  event.put(fTLHEEventStatus, "LHEEventStatus");
  event.put(fTLHEEventMotherFirst, "LHEEventMotherFirst");
  event.put(fTLHEEventMotherSecond, "LHEEventMotherSecond");
  event.put(fTLHEEventPx, "LHEEventPx");
  event.put(fTLHEEventPy, "LHEEventPy");
  event.put(fTLHEEventPz, "LHEEventPz");
  event.put(fTLHEEventE,  "LHEEventE");
  event.put(fTLHEEventM,  "LHEEventM");
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
  event.put(fTRho,       "Rho");
  event.put(fTRhoForIso, "RhoForIso");
  event.put(fTWeight, "Weight");
  event.put(fTHLTResults, "HLTResults");
  event.put(fTHLTPrescale, "HLTPrescale");
  event.put(fTL1PhysResults, "L1PhysResults");
  event.put(fTL1TechResults, "L1TechResults");
  event.put(fTNHLTObjs, "NHLTObjs");
  for ( size_t i=0; i<fTNpaths; ++i ) {
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
  event.put(fTPhoVrtxListStart ,"PhoVrtxListStart");
  event.put(fTJVrtxListStart, "JVrtxListStart");
  event.put(fTMaxGenPartExceed,"MaxGenPartExceed");
  event.put(fTnGenParticles,"nGenParticles");
  event.put(fTgenInfoId,"genInfoId");
  event.put(fTgenInfoStatus,"genInfoStatus");
  event.put(fTgenInfoNMo,"genInfoNMo");
  event.put(fTgenInfoMo1,"genInfoMo1");
  event.put(fTgenInfoMo2,"genInfoMo2");
  event.put(fTPromptnessLevel,"PromptnessLevel");
  event.put(fTgenInfoPt,"genInfoPt");
  event.put(fTgenInfoEta,"genInfoEta");
  event.put(fTgenInfoPhi,"genInfoPhi");
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
  event.put(fTCSCTightHaloID, "CSCTightHaloID");
  event.put(fTPFType1MET, "PFType1MET");
  event.put(fTPFType1METpx, "PFType1METpx");
  event.put(fTPFType1METpy, "PFType1METpy");
  event.put(fTPFType1METphi, "PFType1METphi");
  event.put(fTPFType1METSignificance, "PFType1METSignificance");
  event.put(fTPFType1SumEt, "PFType1SumEt");
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
  event.put(fTMuIsPFMuon, "MuIsPFMuon");
  event.put(fTMuIsStandaloneMuon, "MuIsStandaloneMuon");
  event.put(fTMuPx, "MuPx");
  event.put(fTMuPy, "MuPy");
  event.put(fTMuPz, "MuPz");
  event.put(fTMuPt, "MuPt");
  event.put(fTMuInnerTkPt, "MuInnerTkPt");
  event.put(fTMuTkPtE,"MuTkPtE");
  event.put(fTMuTkD0E,"MuTkD0E");
  event.put(fTMuTkDzE,"MuTkDzE");
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
  event.put(fTMuPfIsoR03ChHad  , "MuPfIsoR03ChHad");
  event.put(fTMuPfIsoR03NeHad  , "MuPfIsoR03NeHad");
  event.put(fTMuPfIsoR03Photon , "MuPfIsoR03Photon");
  event.put(fTMuPfIsoR03NeHadHighThresh  , "MuPfIsoR03NeHadHighThresh");
  event.put(fTMuPfIsoR03PhotonHighThresh , "MuPfIsoR03PhotonHighThresh");
  event.put(fTMuPfIsoR03SumPUPt, "MuPfIsoR03SumPUPt");
  event.put(fTMuPfIsoR04ChHad  , "MuPfIsoR04ChHad");
  event.put(fTMuPfIsoR04NeHad  , "MuPfIsoR04NeHad");
  event.put(fTMuPfIsoR04Photon , "MuPfIsoR04Photon");
  event.put(fTMuPfIsoR04NeHadHighThresh  , "MuPfIsoR04NeHadHighThresh");
  event.put(fTMuPfIsoR04PhotonHighThresh , "MuPfIsoR04PhotonHighThresh");
  event.put(fTMuPfIsoR04SumPUPt, "MuPfIsoR04SumPUPt");
  size_t ipfisotag = 0;
  for ( std::vector<edm::InputTag>::const_iterator it = fMuonPfIsoTagsCustom.begin();
        it != fMuonPfIsoTagsCustom.end(); ++it ) {
    event.put(fTMuPfIsosCustom[ipfisotag++], ("Mu"+(*it).label()).c_str());
  }
  event.put(fTMuEem, "MuEem");
  event.put(fTMuEhad, "MuEhad");
  event.put(fTMuD0BS, "MuD0BS");
  event.put(fTMuD0PV, "MuD0PV");
  event.put(fTMuD03DPV, "MuD03DPV");
  event.put(fTMuD03DE, "MuD03DE");
  event.put(fTMuD0E, "MuD0E");
  event.put(fTMuDzBS, "MuDzBS");
  event.put(fTMuDzPV, "MuDzPV");
  event.put(fTMuDzE, "MuDzE");
  event.put(fTMuNChi2, "MuNChi2");
  event.put(fTMuNGlHits, "MuNGlHits");
  event.put(fTMuNGlMuHits, "MuNGlMuHits");
  event.put(fTMuNMuHits, "MuNMuHits");
  event.put(fTMuNTkHits, "MuNTkHits");
  event.put(fTMuNPxHits, "MuNPxHits");
  event.put(fTMuInnerTkNChi2, "MuInnerTkNChi2");
  event.put(fTMuNSiLayers, "MuNSiLayers");
  event.put(fTMuNMatches, "MuNMatches");
  event.put(fTMuNMatchedStations, "MuNMatchedStations");
  event.put(fTMuNChambers, "MuNChambers");
  event.put(fTMuIsoMVA, "MuIsoMVA");
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
  event.put(fTElIsEB, "ElIsEB");
  event.put(fTElIsEE, "ElIsEE");
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
  event.put(fTElD03DPV, "ElD03DPV");
  event.put(fTElD03DE, "ElD03DE");
  event.put(fTElDzBS, "ElDzBS");
  event.put(fTElDzPV, "ElDzPV");
  event.put(fTElDzE, "ElDzE");
  event.put(fTElRelIso03, "ElRelIso03");
  event.put(fTElRelIso04, "ElRelIso04");
  event.put(fTElPfIsoChHad03, "ElPfIsoChHad03");
  event.put(fTElPfIsoNeHad03, "ElPfIsoNeHad03");
  event.put(fTElPfIsoPhoton03, "ElPfIsoPhoton03");
  event.put(fTElDR03TkSumPt, "ElDR03TkSumPt");
  event.put(fTElDR04TkSumPt, "ElDR04TkSumPt");
  event.put(fTElDR03EcalRecHitSumEt, "ElDR03EcalRecHitSumEt");
  event.put(fTElDR04EcalRecHitSumEt, "ElDR04EcalRecHitSumEt");
  event.put(fTElDR03HcalTowerSumEt, "ElDR03HcalTowerSumEt");
  event.put(fTElDR04HcalTowerSumEt, "ElDR04HcalTowerSumEt");
  ipfisotag = 0;
  for ( std::vector<edm::InputTag>::const_iterator it = fElePfIsoTagsCustom.begin();
        it != fElePfIsoTagsCustom.end(); ++it ) {
    event.put(fTElPfIsosCustom[ipfisotag++], ("El"+(*it).label()).c_str());
  }
  ipfisotag = 0;
  for ( std::vector<edm::InputTag>::const_iterator it = fElePfIsoTagsEvent.begin();
        it != fElePfIsoTagsEvent.end(); ++it ) {
    event.put(fTElPfIsosEvent[ipfisotag++], ("ElEvent"+(*it).label()).c_str());
  }
  event.put(fTElNChi2, "ElNChi2");
  event.put(fTElKfTrkchi2, "ElKfTrkchi2");
  event.put(fTElKfTrkhits, "ElKfTrkhits");
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
  event.put(fTElIDMVATrig, "ElIDMVATrig");
  event.put(fTElIDMVANoTrig, "ElIDMVANoTrig");
  event.put(fTElInGap, "ElInGap");
  event.put(fTElEcalDriven, "ElEcalDriven");
  event.put(fTElTrackerDriven, "ElTrackerDriven");
  event.put(fTElBasicClustersSize, "ElBasicClustersSize");
  event.put(fTElfbrem, "Elfbrem");
  event.put(fTElHcalOverEcal, "ElHcalOverEcal");
  event.put(fTElHcalOverEcalBc, "ElHcalOverEcalBc");
  event.put(fTElE1x5, "ElE1x5");
  event.put(fTElE5x5, "ElE5x5");
  event.put(fTElE2x5Max, "ElE2x5Max");
  event.put(fTElR9, "ElR9");
  event.put(fTElSigmaIetaIeta, "ElSigmaIetaIeta");
  event.put(fTElSigmaIphiIphi, "ElSigmaIphiIphi");
  event.put(fTElScEtaWidth, "ElScEtaWidth");
  event.put(fTElScPhiWidth, "ElScPhiWidth");
  event.put(fTElDeltaPhiSeedClusterAtCalo, "ElDeltaPhiSeedClusterAtCalo");
  event.put(fTElDeltaEtaSeedClusterAtCalo, "ElDeltaEtaSeedClusterAtCalo");
  event.put(fTElDeltaPhiSuperClusterAtVtx, "ElDeltaPhiSuperClusterAtVtx");
  event.put(fTElDeltaEtaSuperClusterAtVtx, "ElDeltaEtaSuperClusterAtVtx");
  event.put(fTElCaloEnergy, "ElCaloEnergy");
  event.put(fTElTrkMomAtVtx, "ElTrkMomAtVtx");
  event.put(fTElESuperClusterOverP, "ElESuperClusterOverP");
  event.put(fTElIoEmIoP, "ElIoEmIoP");
  event.put(fTElEoPout, "ElEoPout");
  event.put(fTElPreShowerOverRaw, "ElPreShowerOverRaw");
  event.put(fTElNumberOfMissingInnerHits, "ElNumberOfMissingInnerHits");
  event.put(fTElSCindex, "ElSCindex");
  event.put(fTElConvPartnerTrkDist, "ElConvPartnerTrkDist");
  event.put(fTElPassConversionVeto,"ElPassConversionVeto");
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
  event.put(fTPhoPassConversionVeto,"PhoPassConversionVeto");
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
  event.put(fTPhoHoverE2012, "PhoHoverE2012");
  event.put(fTPhoSigmaIetaIeta, "PhoSigmaIetaIeta");
  event.put(fTPhoSigmaIetaIphi, "PhoSigmaIetaIphi");
  event.put(fTPhoSigmaIphiIphi, "PhoSigmaIphiIphi");
  event.put(fTPhoS4Ratio, "PhoS4Ratio");
  event.put(fTPhoLambdaRatio, "PhoLambdaRatio");
  event.put(fTPhoSCRawEnergy, "PhoSCRawEnergy");
  event.put(fTPhoSCEtaWidth, "PhoSCEtaWidth");
  event.put(fTPhoSCSigmaPhiPhi, "PhoSCSigmaPhiPhi");
  event.put(fTPhoHasPixSeed, "PhoHasPixSeed");
  event.put(fTPhoHasConvTrks, "PhoHasConvTrks");
  event.put(fTPhoScSeedSeverity, "PhoScSeedSeverity");
  event.put(fTPhoE1OverE9, "PhoE1OverE9");
  event.put(fTPhoS4OverS1, "PhoS4OverS1");
  event.put(fTPhoSigmaEtaEta, "PhoSigmaEtaEta");
  event.put(fTPhoSigmaRR, "PhoSigmaRR");
  event.put(fTPhoHCalIso2012ConeDR03, "PhoHCalIso2012ConeDR03");
  event.put(fTPhoNewIsoPFCharged, "PhoNewIsoPFCharged");
  event.put(fTPhoNewIsoPFPhoton, "PhoNewIsoPFPhoton");
  event.put(fTPhoNewIsoPFNeutral, "PhoNewIsoPFNeutral");
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
//  event.put(fTPhoCone04PhotonIsodR0dEta0pt0, "PhoCone04PhotonIsodR0dEta0pt0");
//  event.put(fTPhoCone04PhotonIsodR0dEta0pt5, "PhoCone04PhotonIsodR0dEta0pt5");
//  event.put(fTPhoCone04PhotonIsodR8dEta0pt0, "PhoCone04PhotonIsodR8dEta0pt0");
//  event.put(fTPhoCone04PhotonIsodR8dEta0pt5, "PhoCone04PhotonIsodR8dEta0pt5");
//  event.put(fTPhoCone01PhotonIsodR045EB070EEdEta015pt08EB1EEmvVtx, "PhoCone01PhotonIsodR045EB070EEdEta015pt08EB1EEmvVtx");
//  event.put(fTPhoCone02PhotonIsodR045EB070EEdEta015pt08EB1EEmvVtx, "PhoCone02PhotonIsodR045EB070EEdEta015pt08EB1EEmvVtx");
//  event.put(fTPhoCone03PhotonIsodR045EB070EEdEta015pt08EB1EEmvVtx, "PhoCone03PhotonIsodR045EB070EEdEta015pt08EB1EEmvVtx");
//  event.put(fTPhoCone04PhotonIsodR045EB070EEdEta015pt08EB1EEmvVtx, "PhoCone04PhotonIsodR045EB070EEdEta015pt08EB1EEmvVtx");
//  event.put(fTPhoCone04NeutralHadronIsodR0dEta0pt0, "PhoCone04NeutralHadronIsodR0dEta0pt0");
//  event.put(fTPhoCone04NeutralHadronIsodR0dEta0pt5, "PhoCone04NeutralHadronIsodR0dEta0pt5");
//  event.put(fTPhoCone04NeutralHadronIsodR0dEta0pt0nocracks, "PhoCone04NeutralHadronIsodR0dEta0pt0nocracks");
//  event.put(fTPhoCone04NeutralHadronIsodR0dEta0pt5nocracks, "PhoCone04NeutralHadronIsodR0dEta0pt5nocracks");
//  event.put(fTPhoCone04NeutralHadronIsodR7dEta0pt0, "PhoCone04NeutralHadronIsodR7dEta0pt0");
//  event.put(fTPhoCone04NeutralHadronIsodR7dEta0pt5, "PhoCone04NeutralHadronIsodR7dEta0pt5");
//  event.put(fTPhoCone01NeutralHadronIsodR0dEta0pt0mvVtx, "PhoCone01NeutralHadronIsodR0dEta0pt0mvVtx");
//  event.put(fTPhoCone02NeutralHadronIsodR0dEta0pt0mvVtx, "PhoCone02NeutralHadronIsodR0dEta0pt0mvVtx");
//  event.put(fTPhoCone03NeutralHadronIsodR0dEta0pt0mvVtx, "PhoCone03NeutralHadronIsodR0dEta0pt0mvVtx");
//  event.put(fTPhoCone04NeutralHadronIsodR0dEta0pt0mvVtx, "PhoCone04NeutralHadronIsodR0dEta0pt0mvVtx");
//  event.put(fTPhoCone04ChargedHadronIsodR0dEta0pt0dz0old, "PhoCone04ChargedHadronIsodR0dEta0pt0dz0old");
//  event.put(fTPhoCone04ChargedHadronIsodR0dEta0pt0PFnoPUold, "PhoCone04ChargedHadronIsodR0dEta0pt0PFnoPUold");
//  event.put(fTPhoCone04ChargedHadronIsodR015dEta0pt0dz0old, "PhoCone04ChargedHadronIsodR015dEta0pt0dz0old");
//  event.put(fTPhoCone04ChargedHadronIsodR015dEta0pt0PFnoPUold, "PhoCone04ChargedHadronIsodR015dEta0pt0PFnoPUold");
//  event.put(fTPhoCone01ChargedHadronIsodR0dEta0pt0dz0, "PhoCone01ChargedHadronIsodR0dEta0pt0dz0");
//  event.put(fTPhoCone01ChargedHadronIsodR0dEta0pt0dz1dxy01, "PhoCone01ChargedHadronIsodR0dEta0pt0dz1dxy01");
//  event.put(fTPhoCone01ChargedHadronIsodR0dEta0pt0PFnoPU, "PhoCone01ChargedHadronIsodR0dEta0pt0PFnoPU");
//  event.put(fTPhoCone01ChargedHadronIsodR015dEta0pt0dz0, "PhoCone01ChargedHadronIsodR015dEta0pt0dz0");
//  event.put(fTPhoCone01ChargedHadronIsodR015dEta0pt0dz1dxy01, "PhoCone01ChargedHadronIsodR015dEta0pt0dz1dxy01");
//  event.put(fTPhoCone01ChargedHadronIsodR015dEta0pt0PFnoPU, "PhoCone01ChargedHadronIsodR015dEta0pt0PFnoPU");
//  event.put(fTPhoCone02ChargedHadronIsodR0dEta0pt0dz0, "PhoCone02ChargedHadronIsodR0dEta0pt0dz0");
//  event.put(fTPhoCone02ChargedHadronIsodR0dEta0pt0dz1dxy01, "PhoCone02ChargedHadronIsodR0dEta0pt0dz1dxy01");
//  event.put(fTPhoCone02ChargedHadronIsodR0dEta0pt0PFnoPU, "PhoCone02ChargedHadronIsodR0dEta0pt0PFnoPU");
//  event.put(fTPhoCone02ChargedHadronIsodR015dEta0pt0dz0, "PhoCone02ChargedHadronIsodR015dEta0pt0dz0");
//  event.put(fTPhoCone02ChargedHadronIsodR015dEta0pt0dz1dxy01, "PhoCone02ChargedHadronIsodR015dEta0pt0dz1dxy01");
//  event.put(fTPhoCone02ChargedHadronIsodR015dEta0pt0PFnoPU, "PhoCone02ChargedHadronIsodR015dEta0pt0PFnoPU");
//  event.put(fTPhoCone03ChargedHadronIsodR0dEta0pt0dz0, "PhoCone03ChargedHadronIsodR0dEta0pt0dz0");
//  event.put(fTPhoCone03ChargedHadronIsodR0dEta0pt0dz1dxy01, "PhoCone03ChargedHadronIsodR0dEta0pt0dz1dxy01");
//  event.put(fTPhoCone03ChargedHadronIsodR0dEta0pt0PFnoPU, "PhoCone03ChargedHadronIsodR0dEta0pt0PFnoPU");
//  event.put(fTPhoCone03ChargedHadronIsodR015dEta0pt0dz0, "PhoCone03ChargedHadronIsodR015dEta0pt0dz0");
//  event.put(fTPhoCone03ChargedHadronIsodR015dEta0pt0dz1dxy01, "PhoCone03ChargedHadronIsodR015dEta0pt0dz1dxy01");
//  event.put(fTPhoCone03ChargedHadronIsodR015dEta0pt0PFnoPU, "PhoCone03ChargedHadronIsodR015dEta0pt0PFnoPU");
//  event.put(fTPhoCone04ChargedHadronIsodR0dEta0pt0dz0, "PhoCone04ChargedHadronIsodR0dEta0pt0dz0");
//  event.put(fTPhoCone04ChargedHadronIsodR0dEta0pt0dz1dxy01, "PhoCone04ChargedHadronIsodR0dEta0pt0dz1dxy01");
//  event.put(fTPhoCone04ChargedHadronIsodR0dEta0pt0PFnoPU, "PhoCone04ChargedHadronIsodR0dEta0pt0PFnoPU");
//  event.put(fTPhoCone04ChargedHadronIsodR015dEta0pt0dz0, "PhoCone04ChargedHadronIsodR015dEta0pt0dz0");
//  event.put(fTPhoCone04ChargedHadronIsodR015dEta0pt0dz1dxy01, "PhoCone04ChargedHadronIsodR015dEta0pt0dz1dxy01");
//  event.put(fTPhoCone04ChargedHadronIsodR015dEta0pt0PFnoPU, "PhoCone04ChargedHadronIsodR015dEta0pt0PFnoPU");
  event.put(fTPhoCiCPFIsoChargedDR03, "PhoCiCPFIsoChargedDR03");
  event.put(fTPhoCiCPFIsoNeutralDR03, "PhoCiCPFIsoNeutralDR03");
  event.put(fTPhoCiCPFIsoPhotonDR03, "PhoCiCPFIsoPhotonDR03");
  event.put(fTPhoCiCPFIsoChargedDR04, "PhoCiCPFIsoChargedDR04");
  event.put(fTPhoCiCPFIsoNeutralDR04, "PhoCiCPFIsoNeutralDR04");
  event.put(fTPhoCiCPFIsoPhotonDR04, "PhoCiCPFIsoPhotonDR04");
  event.put(fTPhoSCEta, "PhoSCEta");
  event.put(fTPhoSCPhiWidth, "PhoSCPhiWidth");
  event.put(fTPhoIDMVA, "PhoIDMVA");
//  event.put(fTPhoConvValidVtx, "PhoConvValidVtx");
//  event.put(fTPhoConvNtracks, "PhoConvNtracks");
//  event.put(fTPhoConvChi2Probability, "PhoConvChi2Probability");
//  event.put(fTPhoConvEoverP, "PhoConvEoverP");
//  event.put(fTNconv, "Nconv");
//  event.put(fTConvValidVtx, "ConvValidVtx");
//  event.put(fTConvNtracks, "ConvNtracks");
//  event.put(fTConvChi2Probability, "ConvChi2Probability");
//  event.put(fTConvEoverP, "ConvEoverP");
//  event.put(fTConvZofPrimVtxFromTrks, "ConvZofPrimVtxFromTrks");
  event.put(fTNgv, "Ngv");
  event.put(fTgvSumPtHi, "gvSumPtHi");
  event.put(fTgvSumPtLo, "gvSumPtLo");
  event.put(fTgvNTkHi, "gvNTkHi");
  event.put(fTgvNTkLo, "gvNTkLo");
  event.put(fTNGoodSuperClusters, "NGoodSuperClusters");
  event.put(fTGoodSCEnergy, "GoodSCEnergy");
  event.put(fTGoodSCEta, "GoodSCEta");
  event.put(fTGoodSCPhi, "GoodSCPhi");
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
  size_t ibtag = 0;
  for ( std::vector<edm::InputTag>::const_iterator it = fBtagTags.begin();
	it != fBtagTags.end(); ++it ) {
    event.put(fTJbTagProb[ibtag++], ("J"+(*it).label()).c_str());
  }
  event.put(fTJPartonFlavour, "JPartonFlavour");
  event.put(fTJMass, "JMass");
  event.put(fTJBetaStar, "JBetaStar");
  event.put(fTJBeta, "JBeta");
  event.put(fTJBetaSq, "JBetaSq");
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
  event.put(fTJMetCorrRawEta, "JMetCorrRawEta"); 
  event.put(fTJMetCorrNoMuPt,  "JMetCorrNoMuPt");  
  event.put(fTJMetCorrRawPt,  "JMetCorrRawPt");  
  event.put(fTJMetCorrPhi,  "JMetCorrPhi");  
  event.put(fTJMetCorrEMF, "JMetCorrEMF"); 
  event.put(fTJMetCorrArea,"JMetCorrArea");
  for ( unsigned int i=0; i<PileupJetIdAlgos.size(); ++i ) {
    std::ostringstream s;
    s << i;
    event.put(fTJPassPileupIDL[i], ("JPassPileupIDL"+s.str()).c_str());
    event.put(fTJPassPileupIDM[i], ("JPassPileupIDM"+s.str()).c_str());
    event.put(fTJPassPileupIDT[i], ("JPassPileupIDT"+s.str()).c_str());
  }
  event.put(fTJQGTagLD,"JQGTagLD");
  event.put(fTJQGTagMLP,"JQGTagMLP");
  event.put(fTJSmearedQGL,"JSmearedQGL");
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
  event.put(fTMETR12, "METR12");
  event.put(fTMETR21, "METR21");

event.put(fTSigma,"Sigma");
event.put(fTGenPhotonIsoDR03,"GenPhotonIsoDR03");
event.put(fTGenPhotonIsoDR04,"GenPhotonIsoDR04");
event.put(fTSCX,"SCX");
event.put(fTSCY,"SCY");
event.put(fTSCZ,"SCZ");
event.put(fTSCXtalListStart,"SCXtalListStart");
event.put(fTSCNXtals,"SCNXtals");
event.put(fTNXtals,"NXtals");
event.put(fTXtalX,"XtalX");
event.put(fTXtalY,"XtalY");
event.put(fTXtalZ,"XtalZ");
event.put(fTXtalEtaWidth,"XtalEtaWidth");
event.put(fTXtalPhiWidth,"XtalPhiWidth");
event.put(fTXtalFront1X,"XtalFront1X");
event.put(fTXtalFront1Y,"XtalFront1Y");
event.put(fTXtalFront1Z,"XtalFront1Z");
event.put(fTXtalFront2X,"XtalFront2X");
event.put(fTXtalFront2Y,"XtalFront2Y");
event.put(fTXtalFront2Z,"XtalFront2Z");
event.put(fTXtalFront3X,"XtalFront3X");
event.put(fTXtalFront3Y,"XtalFront3Y");
event.put(fTXtalFront3Z,"XtalFront3Z");
event.put(fTXtalFront4X,"XtalFront4X");
event.put(fTXtalFront4Y,"XtalFront4Y");
event.put(fTXtalFront4Z,"XtalFront4Z");
event.put(fTNPfCand,"NPfCand");
event.put(fTPfCandPdgId,"PfCandPdgId");
event.put(fTPfCandEta,"PfCandEta");
event.put(fTPfCandPhi,"PfCandPhi");
event.put(fTPfCandEnergy,"PfCandEnergy");
event.put(fTPfCandEcalEnergy,"PfCandEcalEnergy");
event.put(fTPfCandPt,"PfCandPt");
event.put(fTPfCandVx,"PfCandVx");
event.put(fTPfCandVy,"PfCandVy");
event.put(fTPfCandVz,"PfCandVz");
event.put(fTPfCandBelongsToJet,"PfCandBelongsToJet");
//event.put(fTPfCandHasHitInFirstPixelLayer,"PfCandHasHitInFirstPixelLayer");
//event.put(fTPfCandTrackRefPx,"PfCandTrackRefPx");
//event.put(fTPfCandTrackRefPy,"PfCandTrackRefPy");
//event.put(fTPfCandTrackRefPz,"PfCandTrackRefPz");
event.put(fTPhoMatchedPFPhotonOrElectronCand,"PhoMatchedPFPhotonOrElectronCand");
event.put(fTPhoFootprintPfCandsListStart,"PhoFootprintPfCandsListStart");
event.put(fTPhoFootprintPfCands,"PhoFootprintPfCands");
event.put(fTPhoVx,"PhoVx");
event.put(fTPhoVy,"PhoVy");
event.put(fTPhoVz,"PhoVz");
event.put(fTPhoRegrEnergy,"PhoRegrEnergy");
event.put(fTPhoRegrEnergyErr,"PhoRegrEnergyErr");
//event.put(fTPhoCone01PhotonIsodEta015EBdR070EEmvVtx,"PhoCone01PhotonIsodEta015EBdR070EEmvVtx");
//event.put(fTPhoCone02PhotonIsodEta015EBdR070EEmvVtx,"PhoCone02PhotonIsodEta015EBdR070EEmvVtx");
//event.put(fTPhoCone03PhotonIsodEta015EBdR070EEmvVtx,"PhoCone03PhotonIsodEta015EBdR070EEmvVtx");
//event.put(fTPhoCone04PhotonIsodEta015EBdR070EEmvVtx,"PhoCone04PhotonIsodEta015EBdR070EEmvVtx");
//event.put(fTPhoCone01NeutralHadronIsomvVtx,"PhoCone01NeutralHadronIsomvVtx");
//event.put(fTPhoCone02NeutralHadronIsomvVtx,"PhoCone02NeutralHadronIsomvVtx");
//event.put(fTPhoCone03NeutralHadronIsomvVtx,"PhoCone03NeutralHadronIsomvVtx");
//event.put(fTPhoCone04NeutralHadronIsomvVtx,"PhoCone04NeutralHadronIsomvVtx");
//event.put(fTPhoCone01ChargedHadronIsodR02dz02dxy01,"PhoCone01ChargedHadronIsodR02dz02dxy01");
//event.put(fTPhoCone02ChargedHadronIsodR02dz02dxy01,"PhoCone02ChargedHadronIsodR02dz02dxy01");
//event.put(fTPhoCone03ChargedHadronIsodR02dz02dxy01,"PhoCone03ChargedHadronIsodR02dz02dxy01");
//event.put(fTPhoCone04ChargedHadronIsodR02dz02dxy01,"PhoCone04ChargedHadronIsodR02dz02dxy01");
//event.put(fTPhoCone03PFCombinedIso,"PhoCone03PFCombinedIso");
//event.put(fTPhoCone04PFCombinedIso,"PhoCone04PFCombinedIso");
event.put(fTDiphotonsfirst,"Diphotonsfirst");
event.put(fTDiphotonssecond,"Diphotonssecond");
event.put(fTVtxdiphoh2gglobe,"Vtxdiphoh2gglobe");
event.put(fTVtxdiphomva,"Vtxdiphomva");
event.put(fTVtxdiphoproductrank,"Vtxdiphoproductrank");
event.put(fTPhoSCRemovalPFIsoCharged,"PhoSCRemovalPFIsoCharged");
event.put(fTPhoSCRemovalPFIsoChargedPrimVtx,"PhoSCRemovalPFIsoChargedPrimVtx");
event.put(fTPhoSCRemovalPFIsoNeutral,"PhoSCRemovalPFIsoNeutral");
event.put(fTPhoSCRemovalPFIsoPhoton,"PhoSCRemovalPFIsoPhoton");
event.put(fTPhoSCRemovalPFIsoChargedRCone,"PhoSCRemovalPFIsoChargedRCone");
event.put(fTPhoSCRemovalPFIsoChargedPrimVtxRCone,"PhoSCRemovalPFIsoChargedPrimVtxRCone");
event.put(fTPhoSCRemovalPFIsoNeutralRCone,"PhoSCRemovalPFIsoNeutralRCone");
event.put(fTPhoSCRemovalPFIsoPhotonRCone,"PhoSCRemovalPFIsoPhotonRCone");
event.put(fTPhoSCRemovalRConeEta, "PhoSCRemovalRConeEta");
event.put(fTPhoSCRemovalRConePhi, "PhoSCRemovalRConePhi");

}

// Method called once before each run
bool NTupleProducer::beginRun(edm::Run& r, const edm::EventSetup& es){

  resetRunProducts();

  // Retrieve and fill RunTree information
  
  // Need to reset at each run, since "put" deletes the pointers
  fRHLTLabels.reset( new std::vector<std::string>(fHLTLabels) );
  fRPileUpData.reset( new std::vector<std::string>(fPileUpData) );
  fRPileUpMC.reset( new std::vector<std::string>(fPileUpMC) );

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
  *fRMinSCrawPt    = fMinSCrawPt;
  *fRMaxPfCandEta  = fMaxPfCandEta;
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
  *fRMaxNConv      = gMaxNConv;
  *fRMaxNPfCand    = gMaxNPfCand;
  *fRMaxNXtals   = gMaxNXtals;

  ReadEnergyScale(r.run());

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
  r.put(fRMinSCrawPt    ,"MinSCrawPt"    );
  r.put(fRMaxPfCandEta  ,"MaxPfCandEta"  );
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
  r.put(fRMaxNConv,   "MaxNConv"      );
  r.put(fRMaxNPfCand, "MaxNPfCand"    );
  r.put(fRMaxNXtals,  "MaxNXtals"    );

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

//reco::VertexRef NTupleProducer::chargedHadronVertex( const edm::Handle<reco::VertexCollection>& vertices, const reco::PFCandidate& pfcand ) const {
//
//  //PfPileUp candidates!
//
//  //  cout << "chargedHadronVertex finding" << endl;
//
//  reco::TrackBaseRef trackBaseRef( pfcand.trackRef() );
//  
//  size_t  iVertex = 0;
//  unsigned index=0;
//  unsigned nFoundVertex = 0;
//  typedef reco::VertexCollection::const_iterator IV;
//  float bestweight=0;
//  for(IV iv=vertices->begin(); iv!=vertices->end(); ++iv, ++index) {
//
//    const reco::Vertex& vtx = *iv;
//    
//    typedef reco::Vertex::trackRef_iterator IT;
//    
//    // loop on tracks in vertices
//    for(IT iTrack=vtx.tracks_begin(); 
//	iTrack!=vtx.tracks_end(); ++iTrack) {
//       
//      const reco::TrackBaseRef& baseRef = *iTrack;
//
//      // one of the tracks in the vertex is the same as 
//      // the track considered in the function
//      float w = vtx.trackWeight(baseRef);
//      if(baseRef == trackBaseRef ) {
//	//select the vertex for which the track has the highest weight
//	if (w > bestweight){
//	  bestweight=w;
//	  iVertex=index;
//	  nFoundVertex++;
//	} 
//      }
//    }
//  }
//
//  if (nFoundVertex>0){
//    if (nFoundVertex!=1)
//      edm::LogWarning("TrackOnTwoVertex")<<"a track is shared by at least two verteces. Used to be an assert";
//    return reco::VertexRef( vertices, iVertex);
//  }
//  // no vertex found with this track. 
//
//  bool checkClosestZVertex_ = true;
//
//  // optional: as a secondary solution, associate the closest vertex in z
//  if ( checkClosestZVertex_ ) {
//
//    double dzmin = 10000;
//    double ztrack = pfcand.vertex().z();
//    bool foundVertex = false;
//    index = 0;
//    for(IV iv=vertices->begin(); iv!=vertices->end(); ++iv, ++index) {
//
//      double dz = fabs(ztrack - iv->z());
//      if(dz<dzmin) {
//	dzmin = dz; 
//	iVertex = index;
//	foundVertex = true;
//      }
//    }
//
//    if( foundVertex ) 
//      return reco::VertexRef( vertices, iVertex);  
//
//  }
//
//
//  return reco::VertexRef();
//}


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
PhotonInfo NTupleProducer::fillPhotonInfos(int p1, int useAllConvs, float correnergy) 
{
	
  int iConv1 = useAllConvs>0 ? matchPhotonToConversion(p1,useAllConvs) : -1;
	
  if ( iConv1 >= 0) {
    // conversions infos
    return PhotonInfo(p1,
		      TVector3(fTPhoSCX->at(p1),fTPhoSCY->at(p1),fTPhoSCZ->at(p1)),
                      TVector3((*fTBeamspotx),(*fTBeamspoty),(*fTBeamspotz)),
                      conv_vtx[iConv1],
		      (*fTConvNtracks)[iConv1]==1 ? conv_singleleg_momentum[iConv1] : conv_refitted_momentum[iConv1],
                      correnergy<=0 ? (*fTPhoEnergy)[p1] : correnergy,
                      (*fTPhoisEB)[p1],
                      (*fTConvNtracks)[iConv1],
                      (*fTConvValidVtx)[iConv1],
                      (*fTConvChi2Probability)[iConv1],
                      (*fTConvEoverP)[iConv1]
                      );
  } 

  return PhotonInfo(p1, 
		    TVector3(fTPhoSCX->at(p1),fTPhoSCY->at(p1),fTPhoSCZ->at(p1)),
                    TVector3((*fTBeamspotx),(*fTBeamspoty),(*fTBeamspotz)),
                    pho_conv_vtx[p1],
                    pho_conv_refitted_momentum[p1],
		    correnergy<=0 ? fTPhoEnergy->at(p1) : correnergy,
                    fTPhoisEB->at(p1),
		    0,
		    0,
		    0,
                    fTPhoConvEoverP->at(p1)                                                                                                                               
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
  if( useMva ) { return vtxAna.rank(*tmvaReader,tmvaMethod); }
	
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
	  
  } // end if at least one photon is a conversion
	
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
 
  return rankprod;
}

int NTupleProducer::matchPhotonToConversion( int lpho, int useAllConvs) {

  int iMatch=-1;
  TVector3 Photonxyz = TVector3(fTPhoSCX->at(lpho),fTPhoSCY->at(lpho),fTPhoSCZ->at(lpho));

//  float detaMin=999.;
//  float dphiMin=999.;   
  float dRMin = 999.;

  for(int iconv=0; iconv<(*fTNconv); iconv++) {

    TVector3 refittedPairMomentum= (*fTConvNtracks)[iconv]==1 ? conv_singleleg_momentum[iconv] : conv_refitted_momentum[iconv];

    //Conversion Selection
    if ( refittedPairMomentum.Pt() < 10 ) continue;
    if ( useAllConvs==1 && (*fTConvNtracks)[iconv]!=1 ) continue;
    if ( useAllConvs==2 && (*fTConvNtracks)[iconv]!=2 ) continue;
    if ( useAllConvs==3 && (*fTConvNtracks)[iconv]!=1 && (*fTConvNtracks)[iconv]!=2 ) continue;
    if ( (*fTConvNtracks)[iconv]==2 && (!(*fTConvValidVtx)[iconv] || (*fTConvChi2Probability)[iconv]<0.000001) ) continue; // Changed back based on meeting on 21.03.2012

    //New matching technique from meeting on 06.08.12
    TVector3 ConversionVertex = conv_vtx[iconv];
    TVector3 NewPhotonxyz = Photonxyz-ConversionVertex;
    double dR = NewPhotonxyz.DeltaR(refittedPairMomentum);
//    double delta_eta = NewPhotonxyz.Eta()-refittedPairMomentum.Eta();
//    double delta_phi = NewPhotonxyz.DeltaPhi(refittedPairMomentum);

    if ( dR < dRMin ) {
//      detaMin=fabs(delta_eta);
//      dphiMin=fabs(delta_phi);
      dRMin=dR;
      iMatch=iconv;
    }
    
  }

  if ( dRMin< 0.1 ) return iMatch;
  else return -1;
    
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
        dR = reco::deltaR( photon.eta(), photon.phi(), mygenpart.eta(), mygenpart.phi());
        if (dR<dRcone && dR>1e-05){
          etsum += mygenpart.et();
        }
      }
    }
  }

  //if (etsum>0) cout << "GenPartonicIso etsum="<<etsum<<endl;

  return etsum;

}


int NTupleProducer::fexist(char *filename) {
  struct stat buffer ;
  if (stat( filename, &buffer ))
    return 1;
  return 0;
}

void NTupleProducer::SetupVtxAlgoParams2012(VertexAlgoParameters &p){

  // Taken from h2gglobe upstream commit 0ab06402a23a441285c2404bac96963e6c4118c2, Oct 17th 2013
  // AnalysisScripts/common/photonanalysis.dat  and AnalysisScripts/common/vertex_selection.dat

  p.fixTkIndex=0;

  p.rescaleTkPtByError=0;
  p.trackCountThr=0;
  p.highPurityOnly=0;

  p.maxD0Signif=999.;
  p.maxDzSignif=999.;

  p.removeTracksInCone=1;
  p.coneSize=0.05;

  // Conversion Flag: 0 is no conversions. 1 is single leg conversions. 2 is double leg conversions. 3 is both single and double leg conversions.
  p.useAllConversions=2;

  p.sigma1Pix=0.011;
  p.sigma1Tib=0.491;
  p.sigma1Tob=4.196;
  p.sigma1PixFwd=0.057;
  p.sigma1Tid=0.316;
  p.sigma1Tec=1.064;

  p.sigma2Pix=0.023;
  p.sigma2Tib=0.301;
  p.sigma2Tob=1.603;
  p.sigma2PixFwd=0.158;
  p.sigma2Tid=0.413;
  p.sigma2Tec=1.026;

  p.singlelegsigma1Pix=0.007;
  p.singlelegsigma1Tib=0.946;
  p.singlelegsigma1Tob=2.495;
  p.singlelegsigma1PixFwd=0.039;
  p.singlelegsigma1Tid=0.349;
  p.singlelegsigma1Tec=1.581;

  p.singlelegsigma2Pix=0.017;
  p.singlelegsigma2Tib=0.472;
  p.singlelegsigma2Tob=0.540;
  p.singlelegsigma2PixFwd=0.119;
  p.singlelegsigma2Tid=0.392;
  p.singlelegsigma2Tec=1.022;

  p.vtxProbFormula="1.-0.49*(x+1)*(y>0.)";

  mvaVertexSelection=true;
  addConversionToMva=true;

  TString descr = getenv("CMSSW_BASE");
  perVtxMvaWeights=Form("%s/src/h2gglobe/VertexAnalysis/data/%s",descr.Data(),"TMVAClassification_BDTCat_conversions_tmva_407.weights.xml");
  perVtxMvaMethod="BDTCat";
  perEvtMvaWeights=Form("%s/src/h2gglobe/VertexAnalysis/data/%s",descr.Data(),"TMVAClassification_BDTvtxprob2012.weights.xml");
  perEvtMvaMethod="BDTvtxprob2012";

}

TLorentzVector NTupleProducer::get_pho_p4(int phoindex, int vtxInd, float energy){
  assert (vtxInd<(*fTNVrtx));
  PhotonInfo pho(phoindex, TVector3(fTPhoSCX->at(phoindex),fTPhoSCY->at(phoindex),fTPhoSCZ->at(phoindex)), (energy>0) ? energy : fTPhoEnergy->at(phoindex));
  return pho.p4(fTVrtxX->at(vtxInd),fTVrtxY->at(vtxInd),fTVrtxZ->at(vtxInd));
}

float NTupleProducer::pfTkIsoWithVertexCiC(int phoindex, int vtxInd, const reco::PFCandidateCollection &pfcands, int pfToUse,
					float dRmax, float dRvetoBarrel, float dRvetoEndcap, float ptMin, float dzMax, float dxyMax) {
  
  assert (pfToUse==1); // protection

  float dRveto;
  if ((*fTPhoisEB)[phoindex])
    dRveto = dRvetoBarrel;
  else
    dRveto = dRvetoEndcap;
  
  TLorentzVector photonDirectionWrtVtx = get_pho_p4(phoindex, vtxInd, 0);
  
  float sum = 0;
  // Loop over the PFCandidates
  for(reco::PFCandidateCollection::const_iterator pf=pfcands.begin(); pf!=pfcands.end(); pf++) {
    
    //require that PFCandidate is a charged hadron
    if (FindPFCandType(pf->pdgId()) == pfToUse) {
	
      if (!(pf->pt()>0)) continue;

      TLorentzVector pfc;
      pfc.SetPtEtaPhiE(pf->pt(),pf->eta(),pf->phi(),pf->energy());

      if (pfc.Pt() < ptMin)
	continue;

      TVector3 vtx(fTVrtxX->at(vtxInd),fTVrtxY->at(vtxInd),fTVrtxZ->at(vtxInd));
      TVector3 pfCandVtx(pf->vx(),pf->vy(),pf->vz());

      float dz = fabs(pfCandVtx.Z() - vtx.Z());
	
      if (dz > dzMax) 
	continue;
	
      double dxy = (-(pfCandVtx.X() - vtx.X())*pfc.Py() + (pfCandVtx.Y() - vtx.Y())*pfc.Px()) / pfc.Pt();
      if(fabs(dxy) > dxyMax) 
	continue;
	
      float dR = photonDirectionWrtVtx.DeltaR(pfc);
      if(dR > dRmax || dR < dRveto) 
	continue;
	    
      sum += pfc.Pt();
    }
  }
    
  return sum;
}

float NTupleProducer::pfEcalIsoCiC(int phoindex, const reco::PFCandidateCollection &pfcands, int pfToUse, float dRmax, float dRVetoBarrel, 
				   float dRVetoEndcap, float etaStripBarrel, float etaStripEndcap, float thrBarrel, float thrEndcaps) {
  
  assert (pfToUse==0 || pfToUse==2); // protection

  float dRVeto, etaStrip, thr;
  if ((*fTPhoisEB)[phoindex]) {
    dRVeto = dRVetoBarrel;
    etaStrip = etaStripBarrel;
    thr = thrBarrel;
  } else {
    dRVeto = dRVetoEndcap;
    etaStrip = etaStripEndcap;
    thr = thrEndcaps;
  }
  
  float sum = 0;
  for(reco::PFCandidateCollection::const_iterator pf=pfcands.begin(); pf!=pfcands.end(); pf++) {
    
    if (FindPFCandType(pf->pdgId()) == pfToUse) {
      
      if (!(pf->pt()>0)) continue;

      TVector3 pfvtx(pf->vx(),pf->vy(),pf->vz());
      TVector3 phoEcalPos(fTPhoSCX->at(phoindex),fTPhoSCY->at(phoindex),fTPhoSCZ->at(phoindex));
      
      TVector3 photonDirectionWrtVtx = TVector3(phoEcalPos.X() - pfvtx.X(),
						phoEcalPos.Y() - pfvtx.Y(),
						phoEcalPos.Z() - pfvtx.Z());
      
      TLorentzVector pfc;
      pfc.SetPtEtaPhiE(pf->pt(),pf->eta(),pf->phi(),pf->energy());
      
      if( pfc.Pt() < thr ) 
	continue;
      
      float dEta = fabs(photonDirectionWrtVtx.Eta() - pfc.Eta());
      float dR = photonDirectionWrtVtx.DeltaR(pfc.Vect());
      
      if (dEta < etaStrip)
	continue;
      
      if(dR > dRmax || dR < dRVeto)
	continue;
      
      sum += pfc.Pt();
    }
  }
  
  return sum;
}

void NTupleProducer::rescaleClusterShapes(struct_photonIDMVA_variables &str, bool isEB){

  return; // NOT APPLYING ANY SCALING FOR THE RUN-DEPENDENT MC

  //  cout << "WARNING: ARE YOU SURE THAT V7N HAS TO BE RESCALED?" << endl;

  if (str.isrescaled) {
    cout << "ERROR: rescaling same variables twice" << endl;
    return;
  }

  str.isrescaled = true;

  if (fIsRealData) return;

  if (isEB){
    str.r9 = 1.0045*str.r9 + 0.0010;
    str.s4ratio = 1.01894*str.s4ratio - 0.01034;
    str.sieie = 0.891832*str.sieie + 0.0009133;
    str.etawidth =  1.04302*str.etawidth - 0.000618;
    str.phiwidth =  1.00002*str.phiwidth - 0.000371;
  }
  else {
    str.r9 = 1.0086*str.r9 - 0.0007;
    str.s4ratio = 1.04969*str.s4ratio - 0.03642;
    str.sieie = 0.99470*str.sieie + 0.00003;
    str.etawidth =  0.903254*str.etawidth + 0.001346;
    str.phiwidth =  0.99992*str.phiwidth - 0.00000048;
  }

}

float NTupleProducer::GetEnergyScaleCorrection(int run, float eta, float r9){

  if (!fIsRealData) return 1; // apply scale correction on data only

  assert (run==energy_scales.run);

  eta = fabs(eta);
  bool highr9 = (r9>0.94);

  if (highr9){
    if (eta<1) return energy_scales.EBLowEtaGold;
    if (eta<1.5) return energy_scales.EBHighEtaGold;
    if (eta<2) return energy_scales.EELowEtaGold;
    if (eta<3) return energy_scales.EEHighEtaGold;
  }
  else {
    if (eta<1) return energy_scales.EBLowEtaBad;
    if (eta<1.5) return energy_scales.EBHighEtaBad;
    if (eta<2) return energy_scales.EELowEtaBad;
    if (eta<3) return energy_scales.EEHighEtaBad;
  }

  return 1; // protection

}

void NTupleProducer::ReadEnergyScale(int run){ // TO BE IMPLEMENTED

  energy_scales.run=run;
  energy_scales.EBLowEtaGold = 1;
  energy_scales.EBHighEtaGold = 1;
  energy_scales.EELowEtaGold = 1;
  energy_scales.EEHighEtaGold = 1;
  energy_scales.EBLowEtaBad = 1;
  energy_scales.EBHighEtaBad = 1;
  energy_scales.EELowEtaBad = 1;
  energy_scales.EEHighEtaBad = 1;

}

//define this as a plug-in
DEFINE_FWK_MODULE(NTupleProducer);
