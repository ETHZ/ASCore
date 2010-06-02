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
// $Id: NTupleProducer.cc,v 1.61 2010/06/01 10:06:11 predragm Exp $
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
#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"

#include "RecoLocalCalo/EcalRecAlgos/interface/EcalSeverityLevelAlgo.h"
#include "CondFormats/DataRecord/interface/EcalChannelStatusRcd.h"
#include "CondFormats/EcalObjects/interface/EcalChannelStatus.h"

// Data formats
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "DataFormats/BTauReco/interface/JetTag.h"
#include "DataFormats/JetReco/interface/JetTracksAssociation.h"
#include "DataFormats/JetReco/interface/JetID.h"

#include "DataFormats/METReco/interface/METCollection.h"
#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/METReco/interface/CaloMETCollection.h"
#include "DataFormats/METReco/interface/PFMET.h"

#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "CondFormats/L1TObjects/interface/L1GtTriggerMenu.h"
#include "CondFormats/DataRecord/interface/L1GtTriggerMenuRcd.h"

#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Jet.h"

#include "MagneticField/Engine/interface/MagneticField.h"

// Interface
#include "DiLeptonAnalysis/NTupleProducer/interface/NTupleProducer.h"


NTupleProducer::NTupleProducer(const edm::ParameterSet& iConfig){
  // Main settings
  fIsRealData = iConfig.getUntrackedParameter<bool>("isRealData");
  fIsPat      = iConfig.getUntrackedParameter<bool>("isPat");

  // InputTags
  fMuonTag        = iConfig.getUntrackedParameter<edm::InputTag>("tag_muons");
  fElectronTag    = iConfig.getUntrackedParameter<edm::InputTag>("tag_electrons");
  fMuIsoDepTkTag  = iConfig.getUntrackedParameter<edm::InputTag>("tag_muisodeptk");
  fMuIsoDepECTag  = iConfig.getUntrackedParameter<edm::InputTag>("tag_muisodepec");
  fMuIsoDepHCTag  = iConfig.getUntrackedParameter<edm::InputTag>("tag_muisodephc");
  fJetTag         = iConfig.getUntrackedParameter<edm::InputTag>("tag_jets");
  fJetCorrs       = iConfig.getUntrackedParameter<string>("jetCorrs");
  fBtagTag        = iConfig.getUntrackedParameter<edm::InputTag>("tag_btag");
  fJetTracksTag   = iConfig.getUntrackedParameter<edm::InputTag>("tag_jetTracks");
  fJetIDTag       = iConfig.getUntrackedParameter<edm::InputTag>("tag_jetID");
  fMET1Tag        = iConfig.getUntrackedParameter<edm::InputTag>("tag_met1");
  fMET2Tag        = iConfig.getUntrackedParameter<edm::InputTag>("tag_met2");
  fMET3Tag        = iConfig.getUntrackedParameter<edm::InputTag>("tag_met3");
  fMET4Tag        = iConfig.getUntrackedParameter<edm::InputTag>("tag_met4");
  fMET5Tag        = iConfig.getUntrackedParameter<edm::InputTag>("tag_met5");
  fVertexTag      = iConfig.getUntrackedParameter<edm::InputTag>("tag_vertex");
  fTrackTag       = iConfig.getUntrackedParameter<edm::InputTag>("tag_tracks");
  fPhotonTag      = iConfig.getUntrackedParameter<edm::InputTag>("tag_photons");
  fCalTowTag      = iConfig.getUntrackedParameter<edm::InputTag>("tag_caltow");
  fEBRecHitsTag   = iConfig.getUntrackedParameter<edm::InputTag>("tag_EBrechits");
  fEERecHitsTag   = iConfig.getUntrackedParameter<edm::InputTag>("tag_EErechits");
  fGenPartTag     = iConfig.getUntrackedParameter<edm::InputTag>("tag_genpart");
  fGenJetTag      = iConfig.getUntrackedParameter<edm::InputTag>("tag_genjets");
  fL1TriggerTag   = iConfig.getUntrackedParameter<edm::InputTag>("tag_l1trig");
  fHLTTriggerTag  = iConfig.getUntrackedParameter<edm::InputTag>("tag_hlttrig");
  fHBHENoiseResultTag = iConfig.getUntrackedParameter<edm::InputTag>("tag_hcalnoise");

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

  // Create histograms and trees
  // - Histograms with trigger information
  fHhltstat    = fTFileService->make<TH1I>("HLTTriggerStats", "HLTTriggerStatistics", gMaxhltbits+2, 0, gMaxhltbits+2);
  fHl1physstat = fTFileService->make<TH1I>("L1PhysTriggerStats", "L1PhysTriggerStatistics", gMaxl1physbits+2, 0, gMaxl1physbits+2);
  fHl1techstat = fTFileService->make<TH1I>("L1TechTriggerStats", "L1TechTriggerStatistics", gMaxl1techbits+2, 0, gMaxl1techbits+2);
  // - Tree with run information
  fRunTree = fTFileService->make<TTree>("RunInfo", "ETHZRunAnalysisTree");
  // Tree with event information
  fEventTree = fTFileService->make<TTree>("Analysis", "ETHZAnalysisTree");

  edm::LogVerbatim("NTP") << "---------------------------------";
  edm::LogVerbatim("NTP") << " ==> NTupleProducer Constructor ...";
  edm::LogVerbatim("NTP") << endl;
  edm::LogVerbatim("NTP") << "  Processing Real Data: " << (fIsRealData?"ON":"OFF"); 
  edm::LogVerbatim("NTP") << "  Processing PAT:       " << (fIsPat?"ON":"OFF"); 
  edm::LogVerbatim("NTP") << endl;
  edm::LogVerbatim("NTP") << "  Input Tags:";
  edm::LogVerbatim("NTP") << "    fMuonTag        = " << fMuonTag.label()        ;
  edm::LogVerbatim("NTP") << "    fElectronTag    = " << fElectronTag.label()    ;
  edm::LogVerbatim("NTP") << "    fMuIsoDepTkTag  = " << fMuIsoDepTkTag.label()  ;
  edm::LogVerbatim("NTP") << "    fMuIsoDepECTag  = " << fMuIsoDepECTag.label()  ;
  edm::LogVerbatim("NTP") << "    fMuIsoDepHCTag  = " << fMuIsoDepHCTag.label()  ;
  edm::LogVerbatim("NTP") << "    fJetTag         = " << fJetTag.label()         ;
  edm::LogVerbatim("NTP") << "    fJetCorrs       = " << fJetCorrs               ;
  edm::LogVerbatim("NTP") << "    fMET1Tag        = " << fMET1Tag.label()        ;
  edm::LogVerbatim("NTP") << "    fMET2Tag        = " << fMET2Tag.label()        ;
  edm::LogVerbatim("NTP") << "    fMET3Tag        = " << fMET3Tag.label()        ;
  edm::LogVerbatim("NTP") << "    fMET4Tag        = " << fMET4Tag.label()        ;
  edm::LogVerbatim("NTP") << "    fMET5Tag        = " << fMET5Tag.label()        ;
  edm::LogVerbatim("NTP") << "    fVertexTag      = " << fVertexTag.label()      ;
  edm::LogVerbatim("NTP") << "    fTrackTag       = " << fTrackTag.label()       ;
  edm::LogVerbatim("NTP") << "    fPhotonTag      = " << fPhotonTag.label()      ;
  edm::LogVerbatim("NTP") << "    fCalTowTag      = " << fCalTowTag.label()      ;
  edm::LogVerbatim("NTP") << "    fGenPartTag     = " << fGenPartTag.label()     ;
  edm::LogVerbatim("NTP") << "    fGenJetTag      = " << fGenJetTag.label()      ;
  edm::LogVerbatim("NTP") << endl;
  edm::LogVerbatim("NTP") << "  Event Selection Parameters:";
  edm::LogVerbatim("NTP") << "    fMinmupt        = " << fMinmupt    ;
  edm::LogVerbatim("NTP") << "    fMaxmueta       = " << fMaxmueta   ;
  edm::LogVerbatim("NTP") << "    fMinelpt        = " << fMinelpt    ;
  edm::LogVerbatim("NTP") << "    fMaxeleta       = " << fMaxeleta   ;
  edm::LogVerbatim("NTP") << "    fMincorjpt      = " << fMincorjpt  ;
  edm::LogVerbatim("NTP") << "    fMinrawjpt      = " << fMinrawjpt  ;
  edm::LogVerbatim("NTP") << "    fMaxjeta        = " << fMaxjeta    ;
  edm::LogVerbatim("NTP") << "    fMinjemfrac     = " << fMinjemfrac ;
  edm::LogVerbatim("NTP") << "    fMintrkpt       = " << fMintrkpt   ;
  edm::LogVerbatim("NTP") << "    fMaxtrketa      = " << fMaxtrketa  ;
  edm::LogVerbatim("NTP") << "    fMaxtrknchi2    = " << fMaxtrknchi2;
  edm::LogVerbatim("NTP") << "    fMintrknhits    = " << fMintrknhits;
  edm::LogVerbatim("NTP") << "    fMinphopt       = " << fMinphopt   ;
  edm::LogVerbatim("NTP") << "    fMaxphoeta      = " << fMaxphoeta  ;
  edm::LogVerbatim("NTP") << endl;
  edm::LogVerbatim("NTP") << "---------------------------------" ;

  // Create additional jet fillers
  std::vector<edm::ParameterSet> jConfigs = iConfig.getParameter<std::vector<edm::ParameterSet> >("jets");
  for (size_t i=0; i<jConfigs.size(); ++i)
    jetFillers.push_back( new JetFiller(jConfigs[i], fEventTree, fIsPat, fIsRealData) );

}

NTupleProducer::~NTupleProducer(){}

// Method called once for each event
void NTupleProducer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){

  fNTotEvents++;
  using namespace edm;
  using namespace std;
  using namespace reco;
  using reco::MuonCollection;
  using reco::JetTagCollection;

  // Reset all the tree variables
  resetTree();
  for ( std::vector<JetFiller*>::iterator it = jetFillers.begin(); 
        it != jetFillers.end(); ++it )
    (*it)->reset();


  ////////////////////////////////////////////////////////////////////////////////
  // Get the collections /////////////////////////////////////////////////////////
  Handle<View<Muon> > muons;
  iEvent.getByLabel(fMuonTag,muons); // 'muons'

  Handle<View<GsfElectron> > electrons;
  iEvent.getByLabel(fElectronTag, electrons); // 'gsfElectrons'

  // Jets and Jet Correctors
  Handle<View<Jet> > jets;
  iEvent.getByLabel(fJetTag,jets); // 'sisCone5CaloJets'

  // collect information for b-tagging
  Handle<JetTagCollection> jetsAndProbs;
  iEvent.getByLabel(fBtagTag,jetsAndProbs);

  // Jet tracks association (already done in PAT)
  Handle<reco::JetTracksAssociation::Container> jetTracksAssoc;
  if ( !fIsPat ) iEvent.getByLabel(fJetTracksTag,jetTracksAssoc);

  // Jet ID association (already done in PAT)
  edm::Handle<reco::JetIDValueMap> jetIDMap;
  if ( !fIsPat ) iEvent.getByLabel(fJetIDTag,jetIDMap);

  //Get Tracks collection
  Handle<TrackCollection> tracks;
  iEvent.getByLabel(fTrackTag, tracks);

  //Get Photon collection
  Handle<View<Photon> > photons;
  iEvent.getByLabel(fPhotonTag, photons);

  // MET
  Handle<CaloMETCollection> calomet;
  iEvent.getByLabel(fMET1Tag, calomet);

  Handle<CaloMETCollection> corrmumet;
  iEvent.getByLabel(fMET2Tag, corrmumet);

  Handle<METCollection> tcmet;
  iEvent.getByLabel(fMET3Tag, tcmet);

  Handle<View<PFMET> > pfmet;
  iEvent.getByLabel(fMET4Tag, pfmet);

  Handle<CaloMETCollection> corrmujesmet;
  iEvent.getByLabel(fMET5Tag, corrmujesmet);

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
  iEvent.getByLabel(fCalTowTag,calotowers);

  // For ECAL cleaning: rechit and channel status
  edm::Handle<EcalRecHitCollection> ebRecHits;
  edm::Handle<EcalRecHitCollection> eeRecHits;
  iEvent.getByLabel(fEBRecHitsTag,ebRecHits);
  iEvent.getByLabel(fEERecHitsTag,eeRecHits);
  edm::ESHandle<EcalChannelStatus> chStatus;
  iSetup.get<EcalChannelStatusRcd>().get(chStatus);
  const EcalChannelStatus * channelStatus = chStatus.product();

  // Retrieve HB/HE noise flag
  edm::Handle<bool> hbHeNoiseFlag;
  iEvent.getByLabel(fHBHENoiseResultTag,hbHeNoiseFlag);
  fTHBHENoiseFlag = static_cast<int>(*hbHeNoiseFlag);

  // Get Transient Track Builder
  ESHandle<TransientTrackBuilder> theB;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB);

  // Get GenEventInfoProduct
  edm::Handle<GenEventInfoProduct> genEvtInfo;
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
  }

  // Get Magnetic Field
  edm::ESHandle<MagneticField> magfield;
  iSetup.get<IdealMagneticFieldRecord>().get(magfield);


  // Dump trigger bits
  Handle<L1GlobalTriggerReadoutRecord> l1GtReadoutRecord;
  iEvent.getByLabel(fL1TriggerTag, l1GtReadoutRecord);

  Handle<TriggerResults> triggers;
  iEvent.getByLabel(fHLTTriggerTag, triggers);
  const TriggerResults& tr = *triggers;

  if(tr.size() >= gMaxhltbits){
    edm::LogWarning("NTP") << "@SUB=analyze()"
                           << "More than " << static_cast<int>(gMaxhltbits) << " HLT trigger bits, increase length!";
    fTgoodevent = 1;
  }

  if(fFirstevent){
    fFirstevent = false;
    vector<string> triggernames;
    triggernames.reserve(tr.size());
    Service<service::TriggerNamesService> tns;
    tns->getTrigPaths(*triggers, triggernames);
    for( unsigned int i = 0; i < tr.size(); i++ ){
      fHhltstat->GetXaxis()->SetBinLabel(i+1, TString(triggernames[i]));
    }
    // Add a warning about the shift between trigger bits and bin numbers:
    fHhltstat->GetXaxis()->SetBinLabel(gMaxhltbits+2, "Bin#=Bit#+1");
    triggernames.clear();

    ESHandle<L1GtTriggerMenu> menuRcd;
    iSetup.get<L1GtTriggerMenuRcd>().get(menuRcd);
    const L1GtTriggerMenu *menu = menuRcd.product();
    const AlgorithmMap& algoMap=menu->gtAlgorithmMap();
    for (AlgorithmMap::const_iterator it=algoMap.begin();it!=algoMap.end();++it) {
      fHl1physstat->GetXaxis()->SetBinLabel((*it).second.algoBitNumber() + 1, TString((*it).first));
    }
  }

  for(unsigned int i = 0; i < tr.size(); i++ ){
    if(tr[i].accept())	fHhltstat->Fill(i);
    fTHLTres[i] = tr[i].accept() ? 1:0;
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

  ////////////////////////////////////////////////////////////////////////////////
  // Dump tree variables /////////////////////////////////////////////////////////
  fTrunnumber = iEvent.id().run();
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
  fTpvtxndof   = static_cast<int>(primVtx->ndof());
  fTpvtxisfake = primVtx->isFake();
  // fTpvtxntracks = primVtx->tracksSize();

  fTpvtxptsum = 0.;
  for(vector<TrackBaseRef>::const_iterator trackit = primVtx->tracks_begin(); trackit != primVtx->tracks_end(); ++trackit){
    fTpvtxptsum += (*trackit)->pt();
  }
  fTgoodvtx     = 0;


  // Save position of beamspot
  fTbeamspotx = (beamSpot.position()).x();
  fTbeamspoty = (beamSpot.position()).y();
  fTbeamspotz = (beamSpot.position()).z();

	////////////////////////////////////////////////////////
	// Muon Variables:
	int mi(-1), mqi(-1); // index of all muons and qualified muons respectively
	fTnmutot = 0;
	for(View<Muon>::const_iterator Mit = muons->begin(); Mit != muons->end(); ++Mit){
		// Check if maximum number of muons is exceeded already:
		if(mqi >= gMaxnmus){
			edm::LogWarning("NTP") << "@SUB=analyze()"
				<< "Maximum number of muons exceeded";
			fTflagmaxmuexc = 1;
			fTgoodevent = 1;
			break;
		}
		mi++;
		// Muon preselection:
		// Only consider global and trackermuons:
		if(!(Mit->isGlobalMuon()) && !(Mit->isTrackerMuon())) continue;
		fTnmutot++;
		if(Mit->pt() < fMinmupt) continue;
		if(fabs(Mit->eta()) > fMaxmueta) continue;
		mqi++;

		fTmuIsGM[mqi]   = Mit->isGlobalMuon() ? 1:0;
		fTmuIsTM[mqi]   = Mit->isTrackerMuon() ? 1:0;

		// Combined methods for Global and Tracker muons:
		fTmupx[mqi]     = Mit->px();
		fTmupy[mqi]     = Mit->py();
		fTmupz[mqi]     = Mit->pz();
		fTmupt[mqi]     = Mit->pt();
		fTmueta[mqi]    = Mit->eta();
		fTmuphi[mqi]    = Mit->phi();
		fTmue[mqi]      = Mit->energy();
		fTmuet[mqi]     = Mit->et();
		fTmucharge[mqi] = Mit->charge();

		fTmuiso[mqi]            = (Mit->isolationR03().sumPt + Mit->isolationR03().emEt + Mit->isolationR03().hadEt) / Mit->pt();
		fTmuIso03sumPt[mqi]     = Mit->isolationR03().sumPt;
		fTmuIso03emEt[mqi]      = Mit->isolationR03().emEt;
		fTmuIso03hadEt[mqi]     = Mit->isolationR03().hadEt;
		fTmuIso03emVetoEt[mqi]  = Mit->isolationR03().emVetoEt;
		fTmuIso03hadVetoEt[mqi] = Mit->isolationR03().hadVetoEt;
		fTmuIso05sumPt[mqi]     = Mit->isolationR05().sumPt;
		fTmuIso05emEt[mqi]      = Mit->isolationR05().emEt;
		fTmuIso05hadEt[mqi]     = Mit->isolationR05().hadEt;

		fTmucalocomp[mqi] = Mit->caloCompatibility();
		fTmusegmcomp[mqi] = muon::segmentCompatibility(*Mit);

		// MuID Flags:
		fTmuIsGMPT[mqi]                  = muon::isGoodMuon(*Mit, muon::GlobalMuonPromptTight) ? 1:0;
		fTmuIsGMTkChiComp[mqi]           = muon::isGoodMuon(*Mit, muon::GMTkChiCompatibility) ? 1:0;
		fTmuIsGMStaChiComp[mqi]          = muon::isGoodMuon(*Mit, muon::GMStaChiCompatibility) ? 1:0;
		fTmuIsGMTkKinkTight[mqi]         = muon::isGoodMuon(*Mit, muon::GMTkKinkTight) ? 1:0;
		fTmuIsAllStaMuons[mqi]           = muon::isGoodMuon(*Mit, muon::AllStandAloneMuons) ? 1:0;
		fTmuIsAllTrkMuons[mqi]           = muon::isGoodMuon(*Mit, muon::AllTrackerMuons) ? 1:0;
		fTmuIsTrkMuArb[mqi]              = muon::isGoodMuon(*Mit, muon::TrackerMuonArbitrated) ? 1:0;
		fTmuIsAllArb[mqi]                = muon::isGoodMuon(*Mit, muon::AllArbitrated) ? 1:0;
		fTmuIsTMLastStationLoose[mqi]    = muon::isGoodMuon(*Mit, muon::TMLastStationLoose) ? 1:0;
		fTmuIsTMLastStationTight[mqi]    = muon::isGoodMuon(*Mit, muon::TMLastStationTight) ? 1:0;
		fTmuIsTM2DCompLoose[mqi]         = muon::isGoodMuon(*Mit, muon::TM2DCompatibilityLoose) ? 1:0;
		fTmuIsTM2DCompTight[mqi]         = muon::isGoodMuon(*Mit, muon::TM2DCompatibilityTight) ? 1:0;
		fTmuIsTMOneStationLoose[mqi]     = muon::isGoodMuon(*Mit, muon::TMOneStationLoose) ? 1:0;
		fTmuIsTMOneStationTight[mqi]     = muon::isGoodMuon(*Mit, muon::TMOneStationTight) ? 1:0;
		fTmuIsTMLSOPL[mqi]               = muon::isGoodMuon(*Mit, muon::TMLastStationOptimizedLowPtLoose) ? 1:0;
		fTmuIsTMLastStationAngLoose[mqi] = muon::isGoodMuon(*Mit, muon::TMLastStationAngLoose) ? 1:0;
		fTmuIsTMLastStationAngTight[mqi] = muon::isGoodMuon(*Mit, muon::TMLastStationAngTight) ? 1:0;
		fTmuIsTMOneStationAngLoose[mqi]  = muon::isGoodMuon(*Mit, muon::TMOneStationAngLoose) ? 1:0;
		fTmuIsTMOneStationAngTight[mqi]  = muon::isGoodMuon(*Mit, muon::TMOneStationAngTight) ? 1:0;

		// Isolation is embedded in PAT.
		if ( !fIsPat ) {
			Ref<View<Muon> > muonRef(muons,mi);
			const reco::IsoDeposit ECDep = ECDepMap[muonRef];
			const reco::IsoDeposit HCDep = HCDepMap[muonRef];
			fTmueecal[mqi] = ECDep.candEnergy();
			fTmuehcal[mqi] = HCDep.candEnergy();
		} else {
			const pat::Muon* pMuon =  static_cast<const pat::Muon*>(&(*Mit));
			fTmueecal[mqi] = pMuon->ecalIsoDeposit()->candEnergy();
			fTmuehcal[mqi] = pMuon->hcalIsoDeposit()->candEnergy();
		}

		fTmud0bs[mqi] = -1.0*Mit->innerTrack()->dxy(beamSpot.position());
		fTmud0pv[mqi] = -1.0*Mit->innerTrack()->dxy(primVtx->position());
		fTmudzbs[mqi] = Mit->innerTrack()->dz(beamSpot.position());
		fTmudzpv[mqi] = Mit->innerTrack()->dz(primVtx->position());
		fTmuinntknchi2[mqi] = Mit->innerTrack()->normalizedChi2();

		// Separate methods:
		if(fTmuIsTM[mqi]){ // Tracker Muons
			fTntrackermu++;
			fTmuptE[mqi]  = Mit->innerTrack()->ptError();
			fTmud0E[mqi]  = Mit->innerTrack()->dxyError();
			fTmudzE[mqi]  = Mit->innerTrack()->dzError();

			fTmunchi2[mqi]      = fTmuinntknchi2[mqi]; // No difference for TM
			fTmunglhits[mqi]    = 0;
			fTmuntkhits[mqi]    = Mit->innerTrack()->hitPattern().numberOfValidHits();
			fTmunmuhits[mqi]    = 0;
			fTmunmatches[mqi]   = 0;
			fTmunchambers[mqi]  = 0;

			fTmuoutposrad[mqi]   = Mit->innerTrack()->outerRadius();
			fTmuoutposx[mqi]     = Mit->innerTrack()->outerX();
			fTmuoutposy[mqi]     = Mit->innerTrack()->outerY();
			fTmuoutposz[mqi]     = Mit->innerTrack()->outerZ();

			fTmuoutmomx[mqi]     = Mit->innerTrack()->outerMomentum().x();
			fTmuoutmomy[mqi]     = Mit->innerTrack()->outerMomentum().z();
			fTmuoutmomz[mqi]     = Mit->innerTrack()->outerMomentum().y();
			fTmuoutmomphi[mqi]   = Mit->innerTrack()->outerPhi();
			fTmuoutmometa[mqi]   = Mit->innerTrack()->outerEta();
			fTmuoutmomtheta[mqi] = Mit->innerTrack()->outerTheta();

		}
		if(fTmuIsGM[mqi]){ // Global Muons
			fTnglobalmu++;
			fTmuptE[mqi]  = Mit->globalTrack()->ptError();
			fTmud0E[mqi]  = Mit->globalTrack()->dxyError();
			fTmudzE[mqi]  = Mit->globalTrack()->dzError();

			fTmunchi2[mqi]      = Mit->globalTrack()->normalizedChi2();
			fTmunglhits[mqi]    = Mit->globalTrack()->hitPattern().numberOfValidHits();
			fTmuntkhits[mqi]    = Mit->innerTrack()->hitPattern().numberOfValidHits();
			fTmunmuhits[mqi]    = Mit->outerTrack()->hitPattern().numberOfValidHits();
			fTmunmatches[mqi]   = Mit->numberOfMatches();
			fTmunchambers[mqi]  = Mit->numberOfChambers();

			fTmuoutposrad[mqi]   = Mit->globalTrack()->outerRadius();
			fTmuoutposx[mqi]     = Mit->globalTrack()->outerX();
			fTmuoutposy[mqi]     = Mit->globalTrack()->outerY();
			fTmuoutposz[mqi]     = Mit->globalTrack()->outerZ();

			fTmuoutmomx[mqi]     = Mit->globalTrack()->outerMomentum().x();
			fTmuoutmomy[mqi]     = Mit->globalTrack()->outerMomentum().z();
			fTmuoutmomz[mqi]     = Mit->globalTrack()->outerMomentum().y();
			fTmuoutmomphi[mqi]   = Mit->globalTrack()->outerPhi();
			fTmuoutmometa[mqi]   = Mit->globalTrack()->outerEta();
			fTmuoutmomtheta[mqi] = Mit->globalTrack()->outerTheta();

		}

		// MC Matching
		if(!fIsRealData){
			vector<const GenParticle*> MuMatch = matchRecoCand(&(*Mit), iEvent);
			if(MuMatch[0] != NULL){
				fTGenMuId[mqi]       = MuMatch[0]->pdgId();
				fTGenMuStatus[mqi]   = MuMatch[0]->status();
				fTGenMuCharge[mqi]   = MuMatch[0]->charge();
				fTGenMuPt[mqi]       = MuMatch[0]->pt();
				fTGenMuEta[mqi]      = MuMatch[0]->eta();
				fTGenMuPhi[mqi]      = MuMatch[0]->phi();
				fTGenMuE[mqi]        = MuMatch[0]->energy();

				fTGenMuMId[mqi]      = MuMatch[1]->pdgId();
				fTGenMuMStatus[mqi]  = MuMatch[1]->status();
				fTGenMuMCharge[mqi]  = MuMatch[1]->charge();
				fTGenMuMPt[mqi]      = MuMatch[1]->pt();
				fTGenMuMEta[mqi]     = MuMatch[1]->eta();
				fTGenMuMPhi[mqi]     = MuMatch[1]->phi();
				fTGenMuME[mqi]       = MuMatch[1]->energy();

				fTGenMuGMId[mqi]     = MuMatch[2]->pdgId();
				fTGenMuGMStatus[mqi] = MuMatch[2]->status();
				fTGenMuGMCharge[mqi] = MuMatch[2]->charge();
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
	fTnmu = mqi+1;

  ////////////////////////////////////////////////////////
  // Electron variables:
  // Keep pointers to electron superCluster in original collections
  vector<const SuperCluster*> elecPtr;
  vector<const GsfTrack*> trckPtr;
  int eqi(-1),ei(-1); // counts # of qualified electrons
  fTnelestot = 0;
  if(electrons->size() > 0){
    // Read eID results
    vector<Handle<ValueMap<float> > > eIDValueMap(4); 
    // Robust-Loose 
    iEvent.getByLabel( "eidRobustLoose", eIDValueMap[0] ); 
    const ValueMap<float> & eIDmapRL = *eIDValueMap[0] ;
    // Robust-Tight 
    iEvent.getByLabel( "eidRobustTight", eIDValueMap[1] ); 
    const ValueMap<float> & eIDmapRT = *eIDValueMap[1] ;
    // Loose 
    iEvent.getByLabel( "eidLoose", eIDValueMap[2] ); 
    const ValueMap<float> & eIDmapL  = *eIDValueMap[2] ;
    // Tight 
    iEvent.getByLabel( "eidTight", eIDValueMap[3] ); 
    const ValueMap<float> & eIDmapT  = *eIDValueMap[3] ;
    eIDValueMap.clear();

    // Loop over electrons
    for( View<GsfElectron>::const_iterator El = electrons->begin(); 
         El != electrons->end(); ++El ) {
      ++ei;
      // Check if maximum number of electrons is exceeded already:
      if(eqi >= gMaxneles){
        edm::LogWarning("NTP") << "@SUB=analyze"
                               << "Maximum number of electrons exceeded..";
        fTflagmaxelexc = 1;
        fTgoodevent = 1;
        break;
      }
      fTnelestot++;
      // Electron preselection:   
      if(El->pt() < fMinelpt) continue;
      if(fabs(El->eta()) > fMaxeleta) continue;

      // Save the electron SuperCluster pointer
      const SuperCluster* superCluster = &(*(El->superCluster()));
      elecPtr.push_back(superCluster);
      trckPtr.push_back(&(*(El->gsfTrack()) ) );

      // Dump electron properties in tree variables
      eqi++;

      fTepx[eqi]       = El->gsfTrack()->px();
      fTepy[eqi]       = El->gsfTrack()->py();
      fTepz[eqi]       = El->gsfTrack()->pz();
      fTept[eqi]       = El->pt();
      fTeptE[eqi]      = El->gsfTrack()->ptError();
      fTeeta[eqi]      = El->eta();
      fTephi[eqi]      = El->phi();
      fTee[eqi]        = El->energy();
      fTeet[eqi]       = El->et(); // i know: it's the same as pt, still...
      fTed0bs[eqi]     = -1.0*El->gsfTrack()->dxy(beamSpot.position());
      fTed0pv[eqi]     = -1.0*El->gsfTrack()->dxy(primVtx->position());
      fTed0E[eqi]      = El->gsfTrack()->dxyError();
      fTedzbs[eqi]     = El->gsfTrack()->dz(beamSpot.position());
      fTedzpv[eqi]     = El->gsfTrack()->dz(primVtx->position());
      fTedzE[eqi]      = El->gsfTrack()->dzError();
      fTenchi2[eqi]    = El->gsfTrack()->normalizedChi2();
      fTdr03tksumpt[eqi]         = El->dr03TkSumPt();
      fTdr03ecalrechitsumet[eqi] = El->dr03EcalRecHitSumEt();
      fTdr03hcaltowersumet[eqi]  = El->dr03HcalTowerSumEt();
      fTdr04tksumpt[eqi]         = El->dr04TkSumPt();
      fTdr04ecalrechitsumet[eqi] = El->dr04EcalRecHitSumEt();
      fTdr04hcaltowersumet[eqi]  = El->dr04HcalTowerSumEt();
      fTeiso03[eqi]              = (fTdr03tksumpt[eqi] + fTdr03ecalrechitsumet[eqi] + fTdr03hcaltowersumet[eqi]) / fTept[eqi];
	  fTeiso04[eqi]              = (fTdr04tksumpt[eqi] + fTdr04ecalrechitsumet[eqi] + fTdr04hcaltowersumet[eqi]) / fTept[eqi];
      fTecharge[eqi]   = El->charge();
      fTeCInfoIsGsfCtfCons[eqi]      = El->chargeInfo().isGsfCtfConsistent ? 1:0;
      fTeCInfoIsGsfCtfScPixCons[eqi] = El->chargeInfo().isGsfCtfScPixConsistent ? 1:0;
      fTeCInfoIsGsfScPixCons[eqi]    = El->chargeInfo().isGsfScPixConsistent ? 1:0;
      fTeCInfoScPixCharge[eqi]       = El->chargeInfo().scPixCharge;
      if( El->closestCtfTrackRef().isNonnull() ){
        fTeClosestCtfTrackpt[eqi]      = (El->closestCtfTrack().ctfTrack)->pt();
        fTeClosestCtfTracketa[eqi]     = (El->closestCtfTrack().ctfTrack)->eta();
        fTeClosestCtfTrackphi[eqi]     = (El->closestCtfTrack().ctfTrack)->phi();
        fTeClosestCtfTrackcharge[eqi]  = (El->closestCtfTrack().ctfTrack)->charge();
      }
      fTeInGap[eqi]    = El->isGap() ? 1:0;

      fTeEcalDriven[eqi]    = El->ecalDrivenSeed() ? 1:0;		
      fTeTrackerDriven[eqi] = El->trackerDrivenSeed() ? 1:0;


      fTeBasicClustersSize[eqi]         = El->basicClustersSize();
      fTefbrem[eqi]                     = El->fbrem();
      fTeHcalOverEcal[eqi]              = El->hcalOverEcal();                               
      fTeE1x5[eqi]                      = El->e1x5();                                           
      fTeE5x5[eqi]                      = El->e5x5();                                           
      fTeE2x5Max[eqi]                   = El->e2x5Max();                                           
      fTeSigmaIetaIeta[eqi]             = El->sigmaIetaIeta();                         
      fTeDeltaEtaSeedClusterAtCalo[eqi] = El->deltaEtaSeedClusterTrackAtCalo(); 
      fTeDeltaPhiSeedClusterAtCalo[eqi] = El->deltaPhiSeedClusterTrackAtCalo(); 
      fTeDeltaPhiSuperClusterAtVtx[eqi] = El->deltaPhiSuperClusterTrackAtVtx(); 
      fTeDeltaEtaSuperClusterAtVtx[eqi] = El->deltaEtaSuperClusterTrackAtVtx(); 
      fTecaloenergy[eqi]                = El->caloEnergy();
      fTetrkmomatvtx[eqi]               = El->trackMomentumAtVtx().R();
      fTeESuperClusterOverP[eqi]        = El->eSuperClusterOverP();
		fTeNumberOfMissingInnerHits[eqi]  = El->gsfTrack()->trackerExpectedHitsInner().numberOfHits();
      fTetheta[eqi]                     = El->superCluster()->position().theta();
      fTesceta[eqi]                     = El->superCluster()->eta();

      if ( El->superCluster()->seed()->caloID().detector( reco::CaloID::DET_ECAL_BARREL ) ) {
        fTeScSeedSeverity[eqi] = EcalSeverityLevelAlgo::severityLevel( El->superCluster()->seed()->seed(), *ebRecHits, *channelStatus );
        fTeE1OverE9[eqi] = EcalSeverityLevelAlgo::E1OverE9(   El->superCluster()->seed()->seed(), *ebRecHits );
        fTeS4OverS1[eqi] = EcalSeverityLevelAlgo::swissCross( El->superCluster()->seed()->seed(), *ebRecHits );
      } else if ( El->superCluster()->seed()->caloID().detector( reco::CaloID::DET_ECAL_ENDCAP ) ) {
        fTeScSeedSeverity[eqi] = EcalSeverityLevelAlgo::severityLevel( El->superCluster()->seed()->seed(), *eeRecHits, *channelStatus );
        fTeE1OverE9[eqi] = EcalSeverityLevelAlgo::E1OverE9(   El->superCluster()->seed()->seed(), *eeRecHits );
        fTeS4OverS1[eqi] = EcalSeverityLevelAlgo::swissCross( El->superCluster()->seed()->seed(), *eeRecHits );
      } else
        edm::LogWarning("NTP") << "Electron supercluster seed crystal neither in EB nor in EE!";
                                                

      // Read in Electron ID
      if ( !fIsPat ) { // Electron ID is embedded in PAT => switch
        Ref<View<GsfElectron> > electronRef(electrons,ei);
        fTeIDTight[eqi]       = eIDmapT[electronRef]  ? 1:0;
        fTeIDLoose[eqi]       = eIDmapL[electronRef]  ? 1:0;
        fTeIDRobustTight[eqi] = eIDmapRT[electronRef] ? 1:0;
        fTeIDRobustLoose[eqi] = eIDmapRL[electronRef] ? 1:0;
      } else {
        const pat::Electron* pEl = static_cast<const pat::Electron*>(&(*El));
        fTeIDTight[eqi]       = pEl->electronID("eidTight") > 0. ? 1:0;
        fTeIDLoose[eqi]       = pEl->electronID("eidLoose") > 0. ? 1:0;
        fTeIDRobustTight[eqi] = pEl->electronID("eidRobustTight") > 0. ? 1:0;
        fTeIDRobustLoose[eqi] = pEl->electronID("eidRobustLoose") > 0. ? 1:0;
      }
      // Matching
      if(!fIsRealData){
        vector<const GenParticle*> ElMatch = matchRecoCand(&(*El), iEvent);
        if(ElMatch[0] != NULL){
          fTGenElId[eqi]       = ElMatch[0]->pdgId();
          fTGenElStatus[eqi]   = ElMatch[0]->status();
          fTGenElCharge[eqi]   = ElMatch[0]->charge();
          fTGenElPt[eqi]       = ElMatch[0]->pt();
          fTGenElEta[eqi]      = ElMatch[0]->eta();
          fTGenElPhi[eqi]      = ElMatch[0]->phi();
          fTGenElE[eqi]        = ElMatch[0]->energy();

          fTGenElMId[eqi]      = ElMatch[1]->pdgId();
          fTGenElMStatus[eqi]  = ElMatch[1]->status();
          fTGenElMCharge[eqi]  = ElMatch[1]->charge();
          fTGenElMPt[eqi]      = ElMatch[1]->pt();
          fTGenElMEta[eqi]     = ElMatch[1]->eta();
          fTGenElMPhi[eqi]     = ElMatch[1]->phi();
          fTGenElME[eqi]       = ElMatch[1]->energy();

          fTGenElGMId[eqi]     = ElMatch[2]->pdgId();
          fTGenElGMStatus[eqi] = ElMatch[2]->status();
          fTGenElGMCharge[eqi] = ElMatch[2]->charge();
          fTGenElGMPt[eqi]     = ElMatch[2]->pt();
          fTGenElGMEta[eqi]    = ElMatch[2]->eta();
          fTGenElGMPhi[eqi]    = ElMatch[2]->phi();
          fTGenElGME[eqi]      = ElMatch[2]->energy();
        }
        ElMatch.clear();
      }
      // Conversion Information
      GlobalPoint origin_point(0,0,0);
      double Dist(0.), DCot(0.);
      const float bFieldAtOrigin = magfield.product()->inTesla(origin_point).mag();
      reco::TrackRef ConvPartnerTrack = getConversionPartnerTrack(*El, tracks, bFieldAtOrigin, Dist, DCot);
      if( ConvPartnerTrack.isNonnull() ){
        fTeConvPartTrackDist[eqi]   = Dist;
        fTeConvPartTrackDCot[eqi]   = DCot;
        fTeConvPartTrackPt[eqi]     = ConvPartnerTrack->pt();
        fTeConvPartTrackEta[eqi]    = ConvPartnerTrack->eta();
        fTeConvPartTrackPhi[eqi]    = ConvPartnerTrack->phi();
        fTeConvPartTrackCharge[eqi] = ConvPartnerTrack->charge();				
      }
			
      fTeIsInJet[eqi] = -1;
      fTeSharedPx[eqi] = 0.;
      fTeSharedPy[eqi] = 0.;
      fTeSharedPz[eqi] = 0.;
      fTeSharedEnergy[eqi] = 0.;
			
      fTgoodel[eqi] = 0;
      fTeIsIso[eqi] = 1;
      fTeChargeMisIDProb[eqi] = 0;
      fTeDupEl[eqi] = -1;
    }
  }
  fTneles = eqi+1;


  ////////////////////////////////////////////////////////
  // Photon Variables:
  int nqpho(-1);
  // Keep pointers to superclusters for cross cleaning
  vector<const SuperCluster*> photSCs;
  fTnphotonstot = photons->size();
  for( View<Photon>::const_iterator ip = photons->begin(); ip != photons->end(); ++ip ){	
    // Preselection
    if(ip->pt() < fMinphopt) continue;
    if(fabs(ip->eta()) > fMaxphoeta) continue;
    nqpho++;
    // Check if maximum number of photons exceeded
    if(nqpho>=gMaxnphos){
      edm::LogWarning("NTP") << "@SUB=analyze"
                             << "Maximum number of photons exceeded";
      fTflagmaxphoexc = 1;
      fTgoodevent = 1;
      break;
    }
    // Save photon supercluster position
    const SuperCluster* superCluster = &(*(ip->superCluster()));
    photSCs.push_back(superCluster);

    fTPhotPt[nqpho]     = ip->pt();
    fTPhotPx[nqpho]     = ip->px();
    fTPhotPy[nqpho]     = ip->py();
    fTPhotPz[nqpho]     = ip->pz();
    fTPhotEta[nqpho]    = ip->eta();
    fTPhotPhi[nqpho]    = ip->phi();
    fTPhotEnergy[nqpho] = ip->energy();
    fTPhotIso03Ecal[nqpho]      = ip->ecalRecHitSumEtConeDR03();
    fTPhotIso03Hcal[nqpho]      = ip->hcalTowerSumEtConeDR03();
    fTPhotIso03TrkSolid[nqpho]  = ip->trkSumPtSolidConeDR03();
    fTPhotIso03TrkHollow[nqpho] = ip->trkSumPtHollowConeDR03();
    fTPhotIso03[nqpho]          = (fTPhotIso03TrkSolid[nqpho] + fTPhotIso03Ecal[nqpho] + fTPhotIso03Hcal[nqpho]) / fTPhotPt[nqpho];
    fTPhotIso04Ecal[nqpho]      = ip->ecalRecHitSumEtConeDR04();
    fTPhotIso04Hcal[nqpho]      = ip->hcalTowerSumEtConeDR04();
    fTPhotIso04TrkSolid[nqpho]  = ip->trkSumPtSolidConeDR04();
    fTPhotIso04TrkHollow[nqpho] = ip->trkSumPtHollowConeDR04();
    fTPhotIso04[nqpho]          = (fTPhotIso04TrkSolid[nqpho] + fTPhotIso04Ecal[nqpho] + fTPhotIso04Hcal[nqpho]) / fTPhotPt[nqpho];
    fTPhotcaloPosX[nqpho]       = ip->caloPosition().X();
    fTPhotcaloPosY[nqpho]       = ip->caloPosition().Y();
    fTPhotcaloPosZ[nqpho]       = ip->caloPosition().Z();
    fTPhotHoverE[nqpho]         = ip->hadronicOverEm();
    fTPhotH1overE[nqpho]        = ip->hadronicDepth1OverEm();
    fTPhotH2overE[nqpho]        = ip->hadronicDepth2OverEm();
    fTPhotSigmaIetaIeta[nqpho]  = ip->sigmaIetaIeta();
    fTPhotHasPixSeed[nqpho]     = ip->hasPixelSeed() ? 1:0;
    fTPhotHasConvTrks[nqpho]    = ip->hasConversionTracks() ? 1:0; 
    fTPhotIsInJet[nqpho] = -1;
    fTPhotDupEl[nqpho] = -1;
    fTPhotSharedPx[nqpho] = 0.;
    fTPhotSharedPy[nqpho] = 0.;
    fTPhotSharedPz[nqpho] = 0.;
    fTPhotSharedEnergy[nqpho] = 0.;
    fTgoodphoton[nqpho]  = 0;
    fTPhotIsIso[nqpho] = 1;

    // Spike removal information
    if ( ip->superCluster()->seed()->caloID().detector( reco::CaloID::DET_ECAL_BARREL ) ) {
      fTPhotScSeedSeverity[nqpho] = EcalSeverityLevelAlgo::severityLevel( ip->superCluster()->seed()->seed(), *ebRecHits, *channelStatus );
      fTPhotE1OverE9[nqpho] = EcalSeverityLevelAlgo::E1OverE9(   ip->superCluster()->seed()->seed(), *ebRecHits );
      fTPhotS4OverS1[nqpho] = EcalSeverityLevelAlgo::swissCross( ip->superCluster()->seed()->seed(), *ebRecHits );
    } else if ( ip->superCluster()->seed()->caloID().detector( reco::CaloID::DET_ECAL_ENDCAP ) ) {
      fTPhotScSeedSeverity[nqpho] = EcalSeverityLevelAlgo::severityLevel( ip->superCluster()->seed()->seed(), *eeRecHits, *channelStatus );
      fTPhotE1OverE9[nqpho] = EcalSeverityLevelAlgo::E1OverE9(   ip->superCluster()->seed()->seed(), *eeRecHits );
      fTPhotS4OverS1[nqpho] = 1.0-EcalSeverityLevelAlgo::swissCross( ip->superCluster()->seed()->seed(), *eeRecHits );
    } else
      edm::LogWarning("NTP") << "Photon supercluster seed crystal neither in EB nor in EE!";
                                                
  }
  fTnphotons = nqpho+1;

  ////////////////////////////////////////////////////////
  // Jet Variables:
  // Keep pointers to jets in original collections
  vector<const Jet*> jetPtr;
  // Jet corrections
  vector<const Jet*> corrJets;
  corrJets.reserve(jets->size()); // Speed up push-backs
  vector<OrderPair> corrIndices;  // Vector of indices and pt of corr. jets to re-order them
  // Save here px, py, pz for uncorrected jets
  // used later for matching jets with b-tagging information
  int itag(-1);
  fTnjetstot = jets->size();
  // Loop over uncorr. jets
  for(View<Jet>::const_iterator Jit = jets->begin(); Jit != jets->end(); ++Jit){
    // Cut on uncorrected pT (for startup)
    if(Jit->pt() < fMinrawjpt) continue;
    itag++;
    // Save only the gMaxnjets first uncorrected jets
    if(itag >= gMaxnjets){
      edm::LogWarning("NTP") << "@SUB=analyze"
                             << "Found more than " << static_cast<int>(gMaxnjets) << " uncorrected jets, I'm scared ...";
      fTflagmaxujetexc = 1;
      fTgoodevent = 1;
      break;
    }
    if ( !fIsPat ) {
      // Cast and make a new (non-const) copy to apply corrections
      const CaloJet* jet = static_cast<const CaloJet*>(&(*Jit));
      CaloJet* j = new CaloJet(*jet);
      fJUNC_px_match[itag] = j->px();
      fJUNC_py_match[itag] = j->py();
      fJUNC_pz_match[itag] = j->pz();
      // Apply L2 and L3 JetCorrections
      const JetCorrector* jetCorr = JetCorrector::getJetCorrector(fJetCorrs,iSetup);
      j->scaleEnergy(jetCorr->correction(j->p4()));

      corrJets.push_back(j);
      corrIndices.push_back(make_pair(itag,j->pt()));
    } else { 
      // PAT Jets are already corrected: revert to raw
      pat::Jet pJunc = (static_cast<const pat::Jet*>(&(*Jit)))->correctedJet("RAW");
      fJUNC_px_match[itag] = pJunc.px();
      fJUNC_py_match[itag] = pJunc.py();
      fJUNC_pz_match[itag] = pJunc.pz();
      // Sorting not needed for PAT, but still need index
      corrJets.push_back( &(*Jit) );
      corrIndices.push_back(make_pair(itag,(*Jit).pt()));
    } 
  }

  // Sort corrected jet collection by decreasing pt
  IndexByPt indexComparator; // Need this to use the uncorrected jets below...
  std::sort(corrIndices.begin(), corrIndices.end(), indexComparator);

  // Determine corrected jets
  int jqi(-1); // counts # of qualified jets
  // Loop over corr. jet indices
  for(vector<OrderPair>::const_iterator it = corrIndices.begin(); it != corrIndices.end(); ++it ) {
    // Check if maximum number of jets is exceeded already
    if(jqi >= gMaxnjets) {
      edm::LogWarning("NTP") << "@SUB=analyze"
                             << "Maximum number of jets exceeded";
      fTflagmaxjetexc = 1;
      fTgoodevent = 1;
      break;
    }
    int index = it->first;

    // Jet preselection
    const Jet* jet = corrJets[index];
    if(jet->pt() < fMincorjpt) continue;
    if(fabs(jet->eta()) > fMaxjeta) continue;
    // if(jet->emEnergyFraction() < fMinjemfrac) continue;
    jqi++;

    // Save the jet pointer
    jetPtr.push_back(jet);
    // Dump jet properties into tree variables
    fTjpx[jqi]  = jet->px();
    fTjpy[jqi]  = jet->py();
    fTjpz[jqi]  = jet->pz();
    fTjpt[jqi]  = jet->pt();
    fTjeta[jqi] = jet->eta();
    fTjphi[jqi] = jet->phi();
    fTje[jqi]   = jet->energy();
    fTjet[jqi]  = jet->et(); // i know: it's the same as pt, still...
    fTjNconstituents[jqi]  = jet->getJetConstituents().size();
    fTjEcorr[jqi] = jet->px()/fJUNC_px_match[index];

    // Calculate the DR wrt the closest electron
    double ejDRmin = 10.; // Default when no electrons previously selected
    for( int j = 0; j < fTneles; j++ ){
      double ejDR = reco::deltaR(jet->eta(), jet->phi(), fTeeta[j], fTephi[j]);
      if(ejDR<ejDRmin) ejDRmin = ejDR;
    }
    fTjeMinDR[jqi] = ejDRmin;

    // B-tagging probability
    if( !fIsPat ){
      for (unsigned int i = 0; i < jetsAndProbs->size(); i++){
        // Match by pt between the two "collections"
        if( fabs( (fJUNC_px_match[index] - (*jetsAndProbs)[i].first->px())/fJUNC_px_match[index]) < 0.00001 && 
            fabs( (fJUNC_py_match[index] - (*jetsAndProbs)[i].first->py())/fJUNC_py_match[index]) < 0.00001 &&
            fabs( (fJUNC_pz_match[index] - (*jetsAndProbs)[i].first->pz())/fJUNC_pz_match[index]) < 0.00001 ){
          fTjbTagProb[jqi]=(*jetsAndProbs)[i].second;
          break;
        }
      }
    } else {
      const pat::Jet* pJ = static_cast<const pat::Jet*>(jet);      
      fTjbTagProb[jqi] = pJ->bDiscriminator(fBtagTag.label());
    }

    // Jet-track association: get associated tracks
    vector<const reco::Track*> AssociatedTracks;
    if( !fIsPat ){
      edm::RefToBase<reco::Jet> jetRef = jets->refAt(jqi);
      const reco::TrackRefVector& tracks 
        = JetTracksAssociation::getValue(*(jetTracksAssoc.product()),jetRef);
      for ( TrackRefVector::iterator it = tracks.begin();
            it != tracks.end(); ++it )
        AssociatedTracks.push_back( it->get() );
    } else {
      const pat::Jet* pJ = static_cast<const pat::Jet*>(jet);
      const reco::TrackRefVector& tracks = pJ->associatedTracks();
      for ( TrackRefVector::iterator it = tracks.begin();
            it != tracks.end(); ++it )
        AssociatedTracks.push_back( it->get() );
      // DO NOT use this to keep in sync. with non-PAT code
      // fTjChfrac[jqi] = pJ.jetCharge();
      // fTjnAssoTracks[jqi] = pJ.associatedTracks().size();
    }

    //Look inside the calotowers for each jet
    double sumEtaEM=0.;
    double sumEta2EM=0;
    double sumPhiEM=0.;
    double sumPhi2EM=0;
    double sumEtaHAD=0.;
    double sumEta2HAD=0;
    double sumPhiHAD=0.;
    double sumPhi2HAD=0;
    double sumEMWeights=0.;
    double sumHADWeights=0.;
    int NEM=0;
    int NHAD=0;
    double meanEtaEM=0.;
    double meanPhiEM=0.;
    double varianceEtaEM=0.;
    double variancePhiEM=0.;
    double meanEtaHAD=0.;
    double meanPhiHAD=0.;
    double varianceEtaHAD=0.;
    double variancePhiHAD=0.;
    double phioffset=0.; 

    if( !fIsPat ){

      // Cast and make a new copy of corrected calojet, to loop over the calotowers
      const CaloJet* calojet = static_cast<const CaloJet*>(&(*jet));
      vector<CaloTowerPtr> towers = calojet->getCaloConstituents();
      for(vector<CaloTowerPtr>::const_iterator iTower = towers.begin(); iTower!= towers.end(); ++iTower){

        if(iTower == towers.begin()) phioffset = (*iTower)->phi();//calculate all phi coordinates wrt this offset

        //calotowers EM energy
        if((*iTower)->emEnergy()>0.) {
          NEM+=1;
          //eta
          sumEtaEM  +=  (*iTower)->eta()*(*iTower)->emEnergy();
          sumEta2EM += (*iTower)->eta()*(*iTower)->eta()*(*iTower)->emEnergy();
          //phi
          sumPhiEM  += reco::deltaPhi((*iTower)->phi(),phioffset)*(*iTower)->emEnergy();
          sumPhi2EM += reco::deltaPhi((*iTower)->phi(),phioffset)*reco::deltaPhi((*iTower)->phi(),phioffset)*(*iTower)->emEnergy();
          //weights
          sumEMWeights += (*iTower)->emEnergy();
        }
        //calotowers HAD energy
        if((*iTower)->hadEnergy()>0.){
          NHAD += 1;
          //eta
          sumEtaHAD += (*iTower)->eta()*(*iTower)->hadEnergy();
          sumEta2HAD += (*iTower)->eta()*(*iTower)->eta()*(*iTower)->hadEnergy();
          //phi
          sumPhiHAD  += reco::deltaPhi((*iTower)->phi(),phioffset)*(*iTower)->hadEnergy();
          sumPhi2HAD += reco::deltaPhi((*iTower)->phi(),phioffset)*reco::deltaPhi((*iTower)->phi(),phioffset)*(*iTower)->hadEnergy();
          //weight
          sumHADWeights += (*iTower)->hadEnergy();
        }
      }
      if(sumEMWeights>0. && NEM>1){
        meanEtaEM = sumEtaEM/sumEMWeights;
        meanPhiEM = sumPhiEM/sumEMWeights;
        varianceEtaEM = (float(NEM)/float(NEM-1))*(sumEta2EM/sumEMWeights - meanEtaEM*meanEtaEM);
        variancePhiEM = (float(NEM)/float(NEM-1))*(sumPhi2EM/sumEMWeights - meanPhiEM*meanPhiEM);
      }
      if(sumHADWeights>0. && NHAD>1){
        meanEtaHAD = sumEtaHAD/sumHADWeights;
        meanPhiHAD = sumPhiHAD/sumHADWeights;
        varianceEtaHAD = (float(NHAD)/float(NHAD-1))*(sumEta2HAD/sumHADWeights - meanEtaHAD*meanEtaHAD);
        variancePhiHAD = (float(NHAD)/float(NHAD-1))*(sumPhi2HAD/sumHADWeights - meanPhiHAD*meanPhiHAD);
      }
      towers.clear();
    }else { 
      //for the moment do nothing, the variables variance* are initialized to zero
    }
    fTJEtaEMrms[jqi]  = sqrt(varianceEtaEM);
    fTJEtaHADrms[jqi] = sqrt(varianceEtaHAD);
    fTJPhiEMrms[jqi]  = sqrt(variancePhiEM);
    fTJPhiHADrms[jqi] = sqrt(variancePhiHAD);

    //Below save the momenta of the three leading tracks associated to the jet
    double pT1 =0.;
    double pT2 =0.;
    double pT3 =0.;
    int idx1=-1;
    int idx2=-1;
    int idx3=-1;
    // Jet-track association: make transient tracks and store information
    vector<TransientTrack> AssociatedTTracks;
    fTjnAssoTracks[jqi] = 0;
    fTjChfrac[jqi] = -1.; // Default (if jet-tracks association cone is outside tracker acceptance)
    fTjMass[jqi] = 0.;
    if(fabs(jet->eta())<2.9){ // when the cone of dR=0.5 around the jet is (at least partially) inside the tracker acceptance
      // Tmp variables for vectorial sum of pt of tracks
      double pXtmp(0.), pYtmp(0.), pZtmp(0.), E2tmp(0.);
      const double trkmass = 0.; // Assumed mass for tracks
      // Loop over associated tracks:
      for(size_t t = 0; t < AssociatedTracks.size(); ++t){
        AssociatedTTracks.push_back(theB->build(AssociatedTracks[t])); // build transient tracks
        if(AssociatedTracks[t]->normalizedChi2()<10. && AssociatedTracks[t]->numberOfValidHits()>10 && AssociatedTracks[t]->pt()>1.){
          pXtmp += AssociatedTracks[t]->px();
          pYtmp += AssociatedTracks[t]->py();
          pZtmp += AssociatedTracks[t]->pz();
          E2tmp += trkmass*trkmass + pXtmp*pXtmp + pYtmp*pYtmp + pZtmp*pZtmp;
          fTjnAssoTracks[jqi]++;
        }
        // Find the three highest pT tracks
        if(AssociatedTracks[t]->pt()>pT1 && AssociatedTracks.size()>=1){
          pT1=AssociatedTracks[t]->pt();
          idx3=idx2;
          idx2=idx1;
          idx1=t;
        } else if (AssociatedTracks[t]->pt()<pT1 && AssociatedTracks[t]->pt()>pT2 && AssociatedTracks.size()>=2) {
          pT2=AssociatedTracks[t]->pt();
          idx3=idx2;
          idx2=t;
        } else if (AssociatedTracks[t]->pt()<pT2 && AssociatedTracks[t]->pt()>pT3 && AssociatedTracks.size()>=3){
          pT3=AssociatedTracks[t]->pt();
          idx3=t;
        }
      }
      // Fill the momenta
      if(AssociatedTracks.size()>=1){
        fTjtrk1px[jqi] =AssociatedTracks[idx1]->px();
        fTjtrk1py[jqi] =AssociatedTracks[idx1]->py();
        fTjtrk1pz[jqi] =AssociatedTracks[idx1]->pz();
      }
      if(AssociatedTracks.size()>=2){
        fTjtrk2px[jqi] =AssociatedTracks[idx2]->px();
        fTjtrk2py[jqi] =AssociatedTracks[idx2]->py();
        fTjtrk2pz[jqi] =AssociatedTracks[idx2]->pz();
      }
      if(AssociatedTracks.size()>=3){
        fTjtrk3px[jqi] =AssociatedTracks[idx3]->px();
        fTjtrk3py[jqi] =AssociatedTracks[idx3]->py();
        fTjtrk3pz[jqi] =AssociatedTracks[idx3]->pz();
      }

      fTjChfrac[jqi] = sqrt(pXtmp*pXtmp + pYtmp*pYtmp) / jet->pt();
      fTjMass[jqi]   = sqrt(E2tmp - pXtmp*pXtmp - pYtmp*pYtmp - pZtmp*pZtmp);
      fTjMass[jqi]   *= 1/fTjChfrac[jqi];
    } else { // The whole cone used for jet-tracks association is outside of the tracker acceptance
      fTjChfrac[jqi] = -1.;
      fTjMass[jqi]   = -1.;
    }

    // Convert tracks to transient tracks for vertex fitting
    if(AssociatedTTracks.size() > 1) {
      AdaptiveVertexFitter *fitter = new AdaptiveVertexFitter();
      TransientVertex jetVtx = fitter->vertex(AssociatedTTracks);
      if(jetVtx.isValid()){
        fTjetVtxx[jqi]  = jetVtx.position().x();
        fTjetVtxy[jqi]  = jetVtx.position().y();
        fTjetVtxz[jqi]  = jetVtx.position().z();
        fTjetVtxExx[jqi] = jetVtx.positionError().cxx();
        fTjetVtxEyx[jqi] = jetVtx.positionError().cyx();
        fTjetVtxEyy[jqi] = jetVtx.positionError().cyy();
        fTjetVtxEzy[jqi] = jetVtx.positionError().czy();
        fTjetVtxEzz[jqi] = jetVtx.positionError().czz();
        fTjetVtxEzx[jqi] = jetVtx.positionError().czx();
        fTjetVtxNChi2[jqi] = jetVtx.normalisedChiSquared();
      }else{
        fTjetVtxx[jqi]     = -777.77;
        fTjetVtxy[jqi]     = -777.77;
        fTjetVtxz[jqi]     = -777.77;
        fTjetVtxExx[jqi]   = -777.77;
        fTjetVtxEyx[jqi]   = -777.77;
        fTjetVtxEyy[jqi]   = -777.77;
        fTjetVtxEzy[jqi]   = -777.77;
        fTjetVtxEzz[jqi]   = -777.77;
        fTjetVtxEzx[jqi]   = -777.77;
        fTjetVtxNChi2[jqi] = -777.77;
      }
    }else{
      fTjetVtxx[jqi]     = -888.88;
      fTjetVtxy[jqi]     = -888.88;
      fTjetVtxz[jqi]     = -888.88;
      fTjetVtxExx[jqi]   = -888.88;
      fTjetVtxEyx[jqi]   = -888.88;
      fTjetVtxEyy[jqi]   = -888.88;
      fTjetVtxEzy[jqi]   = -888.88;
      fTjetVtxEzz[jqi]   = -888.88;
      fTjetVtxEzx[jqi]   = -888.88;
      fTjetVtxNChi2[jqi] = -888.88;
    }
    AssociatedTracks.clear();
    AssociatedTTracks.clear();

    // CaloJet specific variables (embedded in PAT)
    reco::JetID jetID;
    if ( !fIsPat ) {
      const CaloJet* cJ = static_cast<const CaloJet*>(jet);
      fTjemfrac[jqi]    = cJ->emEnergyFraction();
      fTjEfracHadr[jqi] = cJ->energyFractionHadronic();
      edm::RefToBase<reco::Jet> jetRef = jets->refAt(jqi);
      jetID = (*jetIDMap)[ jetRef ];
    } else {
      const pat::Jet* pJ = static_cast<const pat::Jet*>(jet);
      fTjemfrac[jqi]      = pJ->emEnergyFraction();
      fTjEfracHadr[jqi]   = pJ->energyFractionHadronic();
      jetID = pJ->jetID();
    }
    fTjID_HPD[jqi]      = jetID.fHPD;
    fTjID_RBX[jqi]      = jetID.fRBX;
    fTjID_n90Hits[jqi]  = jetID.n90Hits;
    fTjID_resEMF[jqi]   = jetID.restrictedEMF;
    fTjID_HCALTow[jqi]  = jetID.nHCALTowers;
    fTjID_ECALTow[jqi]  = jetID.nECALTowers;

    // GenJet matching
    if (!fIsRealData) {
      const GenJet *matchedJet = matchJet(&(*jet), iEvent);
      if( matchedJet != NULL ){
        fTjetGenPt[jqi]   = matchedJet->pt();
        fTjetGenEta[jqi]  = matchedJet->eta();
        fTjetGenPhi[jqi]  = matchedJet->phi();
        fTjetGenE[jqi]    = matchedJet->energy();
        fTjetGenemE[jqi]  = matchedJet->emEnergy();
        fTjetGenhadE[jqi] = matchedJet->hadEnergy();
        fTjetGeninvE[jqi] = matchedJet->invisibleEnergy();
      }
    }
    fTgoodjet[jqi] = 0;
  }
  fTnjets = jqi+1;
  corrJets.clear();
  corrIndices.clear();

  // Check electron duplication
  ElectronDuplicate(elecPtr, trckPtr);
  // Check photon/electron duplication
  PhotonElectronDuplicate(elecPtr, photSCs);
  // Check electron/jet duplication
  ElJetOverlap(jetPtr, elecPtr, calotowers);
  // Check photon/jet duplication
  PhotonJetOverlap(jetPtr, photSCs, calotowers);

  // Process other jet collections, as configured
  for ( std::vector<JetFiller*>::iterator it = jetFillers.begin(); 
        it != jetFillers.end(); ++it )
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

  fTMuCorrMET    = (corrmumet->at(0)).pt();
  fTMuCorrMETpx  = (corrmumet->at(0)).px();
  fTMuCorrMETpy  = (corrmumet->at(0)).py();
  fTMuCorrMETphi = (corrmumet->at(0)).phi();

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
  fRunTree->Branch("ExtXSecLO"      ,&fRTextxslo,        "ExtXSecLO/D");
  fRunTree->Branch("ExtXSecNLO"     ,&fRTextxsnlo,       "ExtXSecNLO/D");
  fRunTree->Branch("IntXSec"        ,&fRTintxs,          "IntXSec/D");
  fRunTree->Branch("MinMuPt"        ,&fRTMinmupt,        "MinMuPt/D");
  fRunTree->Branch("MaxMuEta"       ,&fRTMaxmueta,       "MaxMuEta/D");
  fRunTree->Branch("MinElPt"        ,&fRTMinelpt,        "MinElPt/D");
  fRunTree->Branch("MaxElEta"       ,&fRTMaxeleta,       "MaxElEta/D");
  fRunTree->Branch("MinJPt"         ,&fRTMinjpt,         "MinJPt/D");
  fRunTree->Branch("MinRawJPt"      ,&fRTMinrawjpt,      "MinRawJPt/D");
  fRunTree->Branch("MaxJEta"        ,&fRTMaxjeta,        "MaxJEta/D");
  fRunTree->Branch("MinJEMfrac"     ,&fRTMinjemfrac,     "MinJEMfrac/D");

  fRunTree->Branch("MinTrkPt"       ,&fRTMintrkpt,       "MinTrkPt/D");
  fRunTree->Branch("MaxTrkEta"      ,&fRTMaxtrketa,      "MaxTrkEta/D");
  fRunTree->Branch("MaxTrkNChi2"    ,&fRTMaxtrknchi2,    "MaxTrkNChi2/D");
  fRunTree->Branch("MinTrkNHits"    ,&fRTMintrknhits,    "MinTrkNHits/I");

  fRunTree->Branch("MinPhotonPt"    ,&fRTMinphopt,       "MinPhotonPt/D");
  fRunTree->Branch("MaxPhotonEta"   ,&fRTMaxphoeta,      "MaxPhotonEta/D");

  fRunTree->Branch("MaxNMus"        ,&fRTmaxnmu,         "MaxTmu/I");
  fRunTree->Branch("MaxNEles"       ,&fRTmaxnel,         "MaxTel/I");
  fRunTree->Branch("MaxNJets"       ,&fRTmaxnjet,        "MaxTjet/I");
  fRunTree->Branch("MaxNTrks"       ,&fRTmaxntrk,        "MaxTtrk/I");
  fRunTree->Branch("MaxNPhotons"    ,&fRTmaxnphot,       "MaxTphot/I");

  // Event information:
  fEventTree->Branch("Run"              ,&fTrunnumber       ,"Run/I");
  fEventTree->Branch("Event"            ,&fTeventnumber     ,"Event/I");
  fEventTree->Branch("LumiSection"      ,&fTlumisection     ,"LumiSection/I");
  fEventTree->Branch("PtHat"            ,&fTpthat           ,"PtHat/D");
  fEventTree->Branch("SigProcID"        ,&fTsigprocid       ,"SigProcID/I");
  fEventTree->Branch("PDFScalePDF"      ,&fTpdfscalePDF     ,"PDFScalePDF/D");
  fEventTree->Branch("PDFID1"           ,&fTpdfid1          ,"PDFID1/I");
  fEventTree->Branch("PDFID2"           ,&fTpdfid2          ,"PDFID2/I");
  fEventTree->Branch("PDFx1"            ,&fTpdfx1           ,"PDFx1/D");
  fEventTree->Branch("PDFx2"            ,&fTpdfx2           ,"PDFx2/D");
  fEventTree->Branch("PDFxPDF1"         ,&fTpdfxPDF1        ,"PDFxPDF1/D");
  fEventTree->Branch("PDFxPDF2"         ,&fTpdfxPDF2        ,"PDFxPDF2/D");
  fEventTree->Branch("ExtXSecLO"        ,&fTextxslo         ,"ExtXSecLO/D");
  fEventTree->Branch("IntXSec"          ,&fTintxs           ,"IntXSec/D");
  fEventTree->Branch("Weight"           ,&fTweight          ,"Weight/D");
  fEventTree->Branch("HLTResults"       ,&fTHLTres          ,"HLTResults[200]/I");
  fEventTree->Branch("L1PhysResults"    ,&fTL1physres       ,"L1PhysResults[128]/I");
  fEventTree->Branch("L1TechResults"    ,&fTL1techres       ,"L1TechResults[64]/I");
  fEventTree->Branch("PrimVtxGood"      ,&fTgoodvtx         ,"PrimVtxGood/I");
  fEventTree->Branch("PrimVtxx"         ,&fTprimvtxx        ,"PrimVtxx/D");
  fEventTree->Branch("PrimVtxy"         ,&fTprimvtxy        ,"PrimVtxy/D");
  fEventTree->Branch("PrimVtxz"         ,&fTprimvtxz        ,"PrimVtxz/D");	
  fEventTree->Branch("PrimVtxRho"       ,&fTprimvtxrho      ,"PrimVtxRho/D");	
  fEventTree->Branch("PrimVtxxE"        ,&fTprimvtxxE       ,"PrimVtxxE/D");
  fEventTree->Branch("PrimVtxyE"        ,&fTprimvtxyE       ,"PrimVtxyE/D");
  fEventTree->Branch("PrimVtxzE"        ,&fTprimvtxzE       ,"PrimVtxzE/D");
  fEventTree->Branch("PrimVtxNChi2"     ,&fTpvtxznchi2      ,"PrimVtxNChi2/D");
  fEventTree->Branch("PrimVtxNdof"      ,&fTpvtxndof        ,"PrimVtxNdof/I");
  fEventTree->Branch("PrimVtxIsFake"    ,&fTpvtxisfake      ,"PrimVtxIsFake/I");
  fEventTree->Branch("PrimVtxPtSum"     ,&fTpvtxptsum       ,"PrimVtxPtSum/D");
  fEventTree->Branch("Beamspotx"        ,&fTbeamspotx       ,"Beamspotx/D");
  fEventTree->Branch("Beamspoty"        ,&fTbeamspoty       ,"Beamspoty/D");
  fEventTree->Branch("Beamspotz"        ,&fTbeamspotz       ,"Beamspotz/D");
  fEventTree->Branch("NCaloTowers"      ,&fTNCaloTowers     ,"NCaloTowers/I");
  fEventTree->Branch("GoodEvent"        ,&fTgoodevent       ,"GoodEvent/I");
  fEventTree->Branch("MaxMuExceed"      ,&fTflagmaxmuexc    ,"MaxMuExceed/I");
  fEventTree->Branch("MaxElExceed"      ,&fTflagmaxelexc    ,"MaxElExceed/I");
  fEventTree->Branch("MaxJetExceed"     ,&fTflagmaxjetexc   ,"MaxJetExceed/I");
  fEventTree->Branch("MaxUncJetExceed"  ,&fTflagmaxujetexc  ,"MaxUncJetExceed/I");
  fEventTree->Branch("MaxTrkExceed"     ,&fTflagmaxtrkexc   ,"MaxTrkExceed/I");
  fEventTree->Branch("MaxPhotonsExceed" ,&fTflagmaxphoexc   ,"MaxPhotonsExceed/I");
  fEventTree->Branch("HBHENoiseFlag",    &fTHBHENoiseFlag   ,"HBHENoiseFlag/I");

  // Muons:
  fEventTree->Branch("NMus"             ,&fTnmu              ,"NMus/I");
  fEventTree->Branch("NMusTot"          ,&fTnmutot           ,"NMusTot/I");
  fEventTree->Branch("NGMus"            ,&fTnglobalmu        ,"NGMus/I");
  fEventTree->Branch("NTMus"            ,&fTntrackermu       ,"NTMus/I");
  fEventTree->Branch("MuGood"           ,&fTgoodmu           ,"MuGood[NMus]/I");
  fEventTree->Branch("MuIsIso"          ,&fTmuIsIso          ,"MuIsIso[NMus]/I");
  fEventTree->Branch("MuIsGlobalMuon"   ,&fTmuIsGM           ,"MuIsGlobalMuon[NMus]/I");
  fEventTree->Branch("MuIsTrackerMuon"  ,&fTmuIsTM           ,"MuIsTrackerMuon[NMus]/I");
  fEventTree->Branch("MuPx"             ,&fTmupx             ,"MuPx[NMus]/D");
  fEventTree->Branch("MuPy"             ,&fTmupy             ,"MuPy[NMus]/D");
  fEventTree->Branch("MuPz"             ,&fTmupz             ,"MuPz[NMus]/D");
  fEventTree->Branch("MuPt"             ,&fTmupt             ,"MuPt[NMus]/D");
  fEventTree->Branch("MuPtE"            ,&fTmuptE            ,"MuPtE[NMus]/D");
  fEventTree->Branch("MuE"              ,&fTmue              ,"MuE[NMus]/D");
  fEventTree->Branch("MuEt"             ,&fTmuet             ,"MuEt[NMus]/D");
  fEventTree->Branch("MuEta"            ,&fTmueta            ,"MuEta[NMus]/D");
  fEventTree->Branch("MuPhi"            ,&fTmuphi            ,"MuPhi[NMus]/D");
  fEventTree->Branch("MuCharge"         ,&fTmucharge         ,"MuCharge[NMus]/I");
  fEventTree->Branch("MuRelIso03"       ,&fTmuiso            ,"MuRelIso03[NMus]/D");
  fEventTree->Branch("MuIso03SumPt"     ,&fTmuIso03sumPt     ,"MuIso03SumPt[NMus]/D");
  fEventTree->Branch("MuIso03EmEt"      ,&fTmuIso03emEt      ,"MuIso03EmEt[NMus]/D");
  fEventTree->Branch("MuIso03HadEt"     ,&fTmuIso03hadEt     ,"MuIso03HadEt[NMus]/D");
  fEventTree->Branch("MuIso03EMVetoEt"  ,&fTmuIso03emVetoEt  ,"MuIso03EMVetoEt[NMus]/D");
  fEventTree->Branch("MuIso03HadVetoEt" ,&fTmuIso03hadVetoEt ,"MuIso03HadVetoEt[NMus]/D");
  fEventTree->Branch("MuIso05SumPt"     ,&fTmuIso05sumPt     ,"MuIso05SumPt[NMus]/D");
  fEventTree->Branch("MuIso05EmEt"      ,&fTmuIso05emEt      ,"MuIso05EmEt[NMus]/D");
  fEventTree->Branch("MuIso05HadEt"     ,&fTmuIso05hadEt     ,"MuIso05HadEt[NMus]/D");
  fEventTree->Branch("MuEem"            ,&fTmueecal          ,"MuEem[NMus]/D");
  fEventTree->Branch("MuEhad"           ,&fTmuehcal          ,"MuEhad[NMus]/D");
  fEventTree->Branch("MuD0BS"           ,&fTmud0bs           ,"MuD0BS[NMus]/D");
  fEventTree->Branch("MuD0PV"           ,&fTmud0pv           ,"MuD0PV[NMus]/D");
  fEventTree->Branch("MuD0E"            ,&fTmud0E            ,"MuD0E[NMus]/D");
  fEventTree->Branch("MuDzBS"           ,&fTmudzbs           ,"MuDzBS[NMus]/D");
  fEventTree->Branch("MuDzPV"           ,&fTmudzpv           ,"MuDzPV[NMus]/D");
  fEventTree->Branch("MuDzE"            ,&fTmudzE            ,"MuDzE[NMus]/D");
  fEventTree->Branch("MuNChi2"          ,&fTmunchi2          ,"MuNChi2[NMus]/D");
  fEventTree->Branch("MuNGlHits"        ,&fTmunglhits        ,"MuNGlHits[NMus]/I");
  fEventTree->Branch("MuNMuHits"        ,&fTmunmuhits        ,"MuNMuHits[NMus]/I");
  fEventTree->Branch("MuNTkHits"        ,&fTmuntkhits        ,"MuNTkHits[NMus]/I");
  fEventTree->Branch("MuInnerTkNChi2"   ,&fTmuinntknchi2     ,"MuInnerTkNChi2[NMus]/D");
  fEventTree->Branch("MuNMatches"       ,&fTmunmatches       ,"MuNMatches[NMus]/I");
  fEventTree->Branch("MuNChambers"      ,&fTmunchambers      ,"MuNChambers[NMus]/I");
  fEventTree->Branch("MuCaloComp"       ,&fTmucalocomp       ,"MuCaloComp[NMus]/D");
  fEventTree->Branch("MuSegmComp"       ,&fTmusegmcomp       ,"MuSegmComp[NMus]/D");

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

  fEventTree->Branch("MuOutPosRadius" ,&fTmuoutposrad      ,"MuOutPosRadius[NMus]/D");
  fEventTree->Branch("MuOutPosX"      ,&fTmuoutposx        ,"MuOutPosX[NMus]/D");
  fEventTree->Branch("MuOutPosY"      ,&fTmuoutposy        ,"MuOutPosY[NMus]/D");
  fEventTree->Branch("MuOutPosZ"      ,&fTmuoutposz        ,"MuOutPosZ[NMus]/D");
  fEventTree->Branch("MuOutMomx"      ,&fTmuoutmomx        ,"MuOutMomx[NMus]/D");
  fEventTree->Branch("MuOutMomy"      ,&fTmuoutmomy        ,"MuOutMomy[NMus]/D");
  fEventTree->Branch("MuOutMomz"      ,&fTmuoutmomz        ,"MuOutMomz[NMus]/D");
  fEventTree->Branch("MuOutMomPhi"    ,&fTmuoutmomphi      ,"MuOutMomPhi[NMus]/D");
  fEventTree->Branch("MuOutMomEta"    ,&fTmuoutmometa      ,"MuOutMomEta[NMus]/D");
  fEventTree->Branch("MuOutMomTheta"  ,&fTmuoutmomtheta    ,"MuOutMomTheta[NMus]/D");

  fEventTree->Branch("MuGenID"          ,&fTGenMuId         ,"MuGenID[NMus]/I");
  fEventTree->Branch("MuGenStatus"      ,&fTGenMuStatus     ,"MuGenStatus[NMus]/I");
  fEventTree->Branch("MuGenCharge"      ,&fTGenMuCharge     ,"MuGenCharge[NMus]/I");
  fEventTree->Branch("MuGenPt"          ,&fTGenMuPt         ,"MuGenPt[NMus]/D");
  fEventTree->Branch("MuGenEta"         ,&fTGenMuEta        ,"MuGenEta[NMus]/D");
  fEventTree->Branch("MuGenPhi"         ,&fTGenMuPhi        ,"MuGenPhi[NMus]/D");
  fEventTree->Branch("MuGenE"           ,&fTGenMuE          ,"MuGenE[NMus]/D");
  fEventTree->Branch("MuGenMID"         ,&fTGenMuMId        ,"MuGenMID[NMus]/I");
  fEventTree->Branch("MuGenMStatus"     ,&fTGenMuMStatus    ,"MuGenMStatus[NMus]/I");
  fEventTree->Branch("MuGenMCharge"     ,&fTGenMuMCharge    ,"MuGenMCharge[NMus]/I");
  fEventTree->Branch("MuGenMPt"         ,&fTGenMuMPt        ,"MuGenMPt[NMus]/D");
  fEventTree->Branch("MuGenMEta"        ,&fTGenMuMEta       ,"MuGenMEta[NMus]/D");
  fEventTree->Branch("MuGenMPhi"        ,&fTGenMuMPhi       ,"MuGenMPhi[NMus]/D");
  fEventTree->Branch("MuGenME"          ,&fTGenMuME         ,"MuGenME[NMus]/D");
  fEventTree->Branch("MuGenGMID"        ,&fTGenMuGMId       ,"MuGenGMID[NMus]/I");
  fEventTree->Branch("MuGenGMStatus"    ,&fTGenMuGMStatus   ,"MuGenGMStatus[NMus]/I");
  fEventTree->Branch("MuGenGMCharge"    ,&fTGenMuGMCharge   ,"MuGenGMCharge[NMus]/I");
  fEventTree->Branch("MuGenGMPt"        ,&fTGenMuGMPt       ,"MuGenGMPt[NMus]/D");
  fEventTree->Branch("MuGenGMEta"       ,&fTGenMuGMEta      ,"MuGenGMEta[NMus]/D");
  fEventTree->Branch("MuGenGMPhi"       ,&fTGenMuGMPhi      ,"MuGenGMPhi[NMus]/D");
  fEventTree->Branch("MuGenGME"         ,&fTGenMuGME        ,"MuGenGME[NMus]/D");

  // Electrons:
  fEventTree->Branch("NEles"                       ,&fTneles       ,"NEles/I");
  fEventTree->Branch("NElesTot"                    ,&fTnelestot    ,"NElesTot/I");
  fEventTree->Branch("ElGood"                      ,&fTgoodel      ,"ElGood[NEles]/I");
  fEventTree->Branch("ElIsIso"                     ,&fTeIsIso      ,"ElIsIso[NEles]/I");
  fEventTree->Branch("ElChargeMisIDProb"           ,&fTeChargeMisIDProb ,"ElChargeMisIDProb[NEles]/I");
  fEventTree->Branch("ElPx"                        ,&fTepx         ,"ElPx[NEles]/D");
  fEventTree->Branch("ElPy"                        ,&fTepy         ,"ElPy[NEles]/D");
  fEventTree->Branch("ElPz"                        ,&fTepz         ,"ElPz[NEles]/D");
  fEventTree->Branch("ElPt"                        ,&fTept         ,"ElPt[NEles]/D");
  fEventTree->Branch("ElPtE"                       ,&fTeptE        ,"ElPtE[NEles]/D");
  fEventTree->Branch("ElE"                         ,&fTee          ,"ElE[NEles]/D");
  fEventTree->Branch("ElEt"                        ,&fTeet         ,"ElEt[NEles]/D");
  fEventTree->Branch("ElEta"                       ,&fTeeta        ,"ElEta[NEles]/D");
  fEventTree->Branch("ElTheta"                     ,&fTetheta      ,"ElTheta[NEles]/D");
  fEventTree->Branch("ElSCEta"                     ,&fTesceta      ,"ElSCEta[NEles]/D");
  fEventTree->Branch("ElPhi"                       ,&fTephi        ,"ElPhi[NEles]/D");
  fEventTree->Branch("ElD0BS"                      ,&fTed0bs       ,"ElD0BS[NEles]/D");
  fEventTree->Branch("ElD0PV"                      ,&fTed0pv       ,"ElD0PV[NEles]/D");
  fEventTree->Branch("ElD0E"                       ,&fTed0E        ,"ElD0E[NEles]/D");
  fEventTree->Branch("ElDzBS"                      ,&fTedzbs       ,"ElDzBS[NEles]/D");
  fEventTree->Branch("ElDzPV"                      ,&fTedzpv       ,"ElDzPV[NEles]/D");
  fEventTree->Branch("ElDzE"                       ,&fTedzE        ,"ElDzE[NEles]/D");
  fEventTree->Branch("ElRelIso03"                  ,&fTeiso03      ,"ElRelIso03[NEles]/D");
  fEventTree->Branch("ElRelIso04"                  ,&fTeiso04      ,"ElRelIso04[NEles]/D");
  fEventTree->Branch("ElDR03TkSumPt"               ,&fTdr03tksumpt ,"ElDR03TkSumPt[NEles]/D");
  fEventTree->Branch("ElDR04TkSumPt"               ,&fTdr04tksumpt ,"ElDR04TkSumPt[NEles]/D");
  fEventTree->Branch("ElDR03EcalRecHitSumEt"       ,&fTdr03ecalrechitsumet     ,"ElDR03EcalRecHitSumEt[NEles]/D");
  fEventTree->Branch("ElDR04EcalRecHitSumEt"       ,&fTdr04ecalrechitsumet     ,"ElDR04EcalRecHitSumEt[NEles]/D");
  fEventTree->Branch("ElDR03HcalTowerSumEt"        ,&fTdr03hcaltowersumet      ,"ElDR03HcalTowerSumEt[NEles]/D");
  fEventTree->Branch("ElDR04HcalTowerSumEt"        ,&fTdr04hcaltowersumet      ,"ElDR04HcalTowerSumEt[NEles]/D");
  fEventTree->Branch("ElNChi2"                     ,&fTenchi2                  ,"ElNChi2[NEles]/D");
  fEventTree->Branch("ElCharge"                    ,&fTecharge                 ,"ElCharge[NEles]/I");
  fEventTree->Branch("ElCInfoIsGsfCtfCons"         ,&fTeCInfoIsGsfCtfCons      ,"ElCInfoIsGsfCtfCons[NEles]/I");
  fEventTree->Branch("ElCInfoIsGsfCtfScPixCons"    ,&fTeCInfoIsGsfCtfScPixCons ,"ElCInfoIsGsfCtfScPixCons[NEles]/I");
  fEventTree->Branch("ElCInfoIsGsfScPixCons"       ,&fTeCInfoIsGsfScPixCons    ,"ElCInfoIsGsfScPixCons[NEles]/I");
  fEventTree->Branch("ElScPixCharge"               ,&fTeCInfoScPixCharge       ,"ElScPixCharge[NEles]/I");
  fEventTree->Branch("ElClosestCtfTrackPt"         ,&fTeClosestCtfTrackpt      ,"ElClosestCtfTrackPt[NEles]/D");
  fEventTree->Branch("ElClosestCtfTrackEta"        ,&fTeClosestCtfTracketa     ,"ElClosestCtfTrackEta[NEles]/D");
  fEventTree->Branch("ElClosestCtfTrackPhi"        ,&fTeClosestCtfTrackphi     ,"ElClosestCtfTrackPhi[NEles]/D");
  fEventTree->Branch("ElClosestCtfTrackCharge"     ,&fTeClosestCtfTrackcharge  ,"ElClosestCtfTrackCharge[NEles]/I");
  fEventTree->Branch("ElIDTight"                   ,&fTeIDTight           ,"ElIDTight[NEles]/I");
  fEventTree->Branch("ElIDLoose"                   ,&fTeIDLoose           ,"ElIDLoose[NEles]/I");
  fEventTree->Branch("ElIDRobustTight"             ,&fTeIDRobustTight     ,"ElIDRobustTight[NEles]/I");
  fEventTree->Branch("ElIDRobustLoose"             ,&fTeIDRobustLoose     ,"ElIDRobustLoose[NEles]/I");
  fEventTree->Branch("ElInGap"                     ,&fTeInGap             ,"ElInGap[NEles]/I");
  fEventTree->Branch("ElEcalDriven"                ,&fTeEcalDriven        ,"ElEcalDriven[NEles]/I");
  fEventTree->Branch("ElTrackerDriven"             ,&fTeTrackerDriven     ,"ElTrackerDriven[NEles]/I");
  fEventTree->Branch("ElBasicClustersSize"         ,&fTeBasicClustersSize ,"ElBasicClustersSize[NEles]/I");
  fEventTree->Branch("Elfbrem"                     ,&fTefbrem             ,"Elfbrem[NEles]/D");
  fEventTree->Branch("ElHcalOverEcal"              ,&fTeHcalOverEcal      ,"ElHcalOverEcal[NEles]/D");
  fEventTree->Branch("ElE1x5"                      ,&fTeE1x5              ,"ElE1x5[NEles]/D");
  fEventTree->Branch("ElE5x5"                      ,&fTeE5x5              ,"ElE5x5[NEles]/D");
  fEventTree->Branch("ElE2x5Max"                   ,&fTeE2x5Max           ,"ElE2x5Max[NEles]/D");
  fEventTree->Branch("ElSigmaIetaIeta"             ,&fTeSigmaIetaIeta     ,"ElSigmaIetaIeta[NEles]/D");
  fEventTree->Branch("ElDeltaPhiSeedClusterAtCalo" ,&fTeDeltaPhiSeedClusterAtCalo ,"ElDeltaPhiSeedClusterAtCalo[NEles]/D");
  fEventTree->Branch("ElDeltaEtaSeedClusterAtCalo" ,&fTeDeltaEtaSeedClusterAtCalo ,"ElDeltaEtaSeedClusterAtCalo[NEles]/D");
  fEventTree->Branch("ElDeltaPhiSuperClusterAtVtx" ,&fTeDeltaPhiSuperClusterAtVtx ,"ElDeltaPhiSuperClusterAtVtx[NEles]/D");
  fEventTree->Branch("ElDeltaEtaSuperClusterAtVtx" ,&fTeDeltaEtaSuperClusterAtVtx ,"ElDeltaEtaSuperClusterAtVtx[NEles]/D");
  fEventTree->Branch("ElCaloEnergy"                ,&fTecaloenergy        ,"ElCaloEnergy[NEles]/D");
  fEventTree->Branch("ElTrkMomAtVtx"               ,&fTetrkmomatvtx        ,"ElTrkMomAtVtx[NEles]/D");
  fEventTree->Branch("ElESuperClusterOverP"        ,&fTeESuperClusterOverP        ,"ElESuperClusterOverP[NEles]/D");
  fEventTree->Branch("ElNumberOfMissingInnerHits"  ,&fTeNumberOfMissingInnerHits  ,"ElNumberOfMissingInnerHits[NEles]/I");

  fEventTree->Branch("ElIsInJet"                   ,&fTeIsInJet           ,"ElIsInJet[NEles]/I");
  fEventTree->Branch("ElSharedPx"                  ,&fTeSharedPx          ,"ElSharedPx[NEles]/D");
  fEventTree->Branch("ElSharedPy"                  ,&fTeSharedPy          ,"ElSharedPy[NEles]/D");
  fEventTree->Branch("ElSharedPz"                  ,&fTeSharedPz          ,"ElSharedPz[NEles]/D");
  fEventTree->Branch("ElSharedEnergy"              ,&fTeSharedEnergy      ,"ElSharedEnergy[NEles]/D");
  fEventTree->Branch("ElDuplicateEl"               ,&fTeDupEl             ,"ElDuplicateEl[NEles]/I");
  fEventTree->Branch("ElConvPartnerTrkDist"        ,&fTeConvPartTrackDist   ,"ElConvPartnerTrkDist[NEles]/D");
  fEventTree->Branch("ElConvPartnerTrkDCot"        ,&fTeConvPartTrackDCot   ,"ElConvPartnerTrkDCot[NEles]/D");
  fEventTree->Branch("ElConvPartnerTrkPt"          ,&fTeConvPartTrackPt     ,"ElConvPartnerTrkPt[NEles]/D");
  fEventTree->Branch("ElConvPartnerTrkEta"         ,&fTeConvPartTrackEta    ,"ElConvPartnerTrkEta[NEles]/D");
  fEventTree->Branch("ElConvPartnerTrkPhi"         ,&fTeConvPartTrackPhi    ,"ElConvPartnerTrkPhi[NEles]/D");
  fEventTree->Branch("ElConvPartnerTrkCharge"      ,&fTeConvPartTrackCharge ,"ElConvPartnerTrkCharge[NEles]/D");
  fEventTree->Branch("ElScSeedSeverity"            ,&fTeScSeedSeverity      ,"ElScSeedSeverity[NEles]/I");
  fEventTree->Branch("ElE1OverE9"                  ,&fTeE1OverE9            ,"ElE1OverE9[NEles]/D");
  fEventTree->Branch("ElS4OverS1"                  ,&fTeS4OverS1            ,"ElS4OverS1[NEles]/D");
  fEventTree->Branch("ElGenID"                     ,&fTGenElId         ,"ElGenID[NEles]/I");
  fEventTree->Branch("ElGenStatus"                 ,&fTGenElStatus     ,"ElGenStatus[NEles]/I");
  fEventTree->Branch("ElGenCharge"                 ,&fTGenElCharge     ,"ElGenCharge[NEles]/I");
  fEventTree->Branch("ElGenPt"                     ,&fTGenElPt         ,"ElGenPt[NEles]/D");
  fEventTree->Branch("ElGenEta"                    ,&fTGenElEta        ,"ElGenEta[NEles]/D");
  fEventTree->Branch("ElGenPhi"                    ,&fTGenElPhi        ,"ElGenPhi[NEles]/D");
  fEventTree->Branch("ElGenE"                      ,&fTGenElE          ,"ElGenE[NEles]/D");
  fEventTree->Branch("ElGenMID"                    ,&fTGenElMId        ,"ElGenMID[NEles]/I");
  fEventTree->Branch("ElGenMStatus"                ,&fTGenElMStatus    ,"ElGenMStatus[NEles]/I");
  fEventTree->Branch("ElGenMCharge"                ,&fTGenElMCharge    ,"ElGenMCharge[NEles]/I");
  fEventTree->Branch("ElGenMPt"                    ,&fTGenElMPt        ,"ElGenMPt[NEles]/D");
  fEventTree->Branch("ElGenMEta"                   ,&fTGenElMEta       ,"ElGenMEta[NEles]/D");
  fEventTree->Branch("ElGenMPhi"                   ,&fTGenElMPhi       ,"ElGenMPhi[NEles]/D");
  fEventTree->Branch("ElGenME"                     ,&fTGenElME         ,"ElGenME[NEles]/D");
  fEventTree->Branch("ElGenGMID"                   ,&fTGenElGMId       ,"ElGenGMID[NEles]/I");
  fEventTree->Branch("ElGenGMStatus"               ,&fTGenElGMStatus   ,"ElGenGMStatus[NEles]/I");
  fEventTree->Branch("ElGenGMCharge"               ,&fTGenElGMCharge   ,"ElGenGMCharge[NEles]/I");
  fEventTree->Branch("ElGenGMPt"                   ,&fTGenElGMPt       ,"ElGenGMPt[NEles]/D");
  fEventTree->Branch("ElGenGMEta"                  ,&fTGenElGMEta      ,"ElGenGMEta[NEles]/D");
  fEventTree->Branch("ElGenGMPhi"                  ,&fTGenElGMPhi      ,"ElGenGMPhi[NEles]/D");
  fEventTree->Branch("ElGenGME"                    ,&fTGenElGME        ,"ElGenGME[NEles]/D");

  // Photons:
  fEventTree->Branch("NPhotons"         ,&fTnphotons          ,"NPhotons/I");
  fEventTree->Branch("NPhotonsTot"      ,&fTnphotonstot       ,"NPhotonsTot/I");
  fEventTree->Branch("PhoGood"          ,&fTgoodphoton        ,"PhoGood[NPhotons]/I");
  fEventTree->Branch("PhoIsIso"         ,&fTPhotIsIso         ,"PhoIsIso[NPhotons]/I");
  fEventTree->Branch("PhoPt"            ,&fTPhotPt            ,"PhoPt[NPhotons]/D");
  fEventTree->Branch("PhoPx"            ,&fTPhotPx            ,"PhoPx[NPhotons]/D");
  fEventTree->Branch("PhoPy"            ,&fTPhotPy            ,"PhoPy[NPhotons]/D");
  fEventTree->Branch("PhoPz"            ,&fTPhotPz            ,"PhoPz[NPhotons]/D");
  fEventTree->Branch("PhoEta"           ,&fTPhotEta           ,"PhoEta[NPhotons]/D");
  fEventTree->Branch("PhoPhi"           ,&fTPhotPhi           ,"PhoPhi[NPhotons]/D");
  fEventTree->Branch("PhoEnergy"        ,&fTPhotEnergy        ,"PhoEnergy[NPhotons]/D");
  fEventTree->Branch("PhoIso03Ecal"     ,&fTPhotIso03Ecal     ,"PhoIso03Ecal[NPhotons]/D");
  fEventTree->Branch("PhoIso03Hcal"     ,&fTPhotIso03Hcal     ,"PhoIso03Hcal[NPhotons]/D");
  fEventTree->Branch("PhoIso03TrkSolid" ,&fTPhotIso03TrkSolid ,"PhoIso03TrkSolid[NPhotons]/D");
  fEventTree->Branch("PhoIso03TrkHollow",&fTPhotIso03TrkHollow,"PhoIso03TrkHollow[NPhotons]/D");
  fEventTree->Branch("PhoIso03"         ,&fTPhotIso03         ,"PhoIso03[NPhotons]/D");	
  fEventTree->Branch("PhoIso04Ecal"     ,&fTPhotIso04Ecal     ,"PhoIso04Ecal[NPhotons]/D");
  fEventTree->Branch("PhoIso04Hcal"     ,&fTPhotIso04Hcal     ,"PhoIso04Hcal[NPhotons]/D");
  fEventTree->Branch("PhoIso04TrkSolid" ,&fTPhotIso04TrkSolid ,"PhoIso04TrkSolid[NPhotons]/D");
  fEventTree->Branch("PhoIso04TrkHollow",&fTPhotIso04TrkHollow,"PhoIso04TrkHollow[NPhotons]/D");
  fEventTree->Branch("PhoIso04"         ,&fTPhotIso04         ,"PhoIso04[NPhotons]/D");	
  fEventTree->Branch("PhoCaloPositionX" ,&fTPhotcaloPosX      ,"PhoCaloPositionX[NPhotons]/D");
  fEventTree->Branch("PhoCaloPositionY" ,&fTPhotcaloPosY      ,"PhoCaloPositionY[NPhotons]/D");
  fEventTree->Branch("PhoCaloPositionZ" ,&fTPhotcaloPosZ      ,"PhoCaloPositionZ[NPhotons]/D");
  fEventTree->Branch("PhoHoverE"        ,&fTPhotHoverE        ,"PhoHoverE[NPhotons]/D");
  fEventTree->Branch("PhoH1overE"       ,&fTPhotH1overE       ,"PhoH1overE[NPhotons]/D");
  fEventTree->Branch("PhoH2overE"       ,&fTPhotH2overE       ,"PhoH2overE[NPhotons]/D");
  fEventTree->Branch("PhoSigmaIetaIeta" ,&fTPhotSigmaIetaIeta ,"PhoSigmaIetaIeta[NPhotons]/D");
  fEventTree->Branch("PhoHasPixSeed"    ,&fTPhotHasPixSeed    ,"PhoHasPixSeed[NPhotons]/I");
  fEventTree->Branch("PhoHasConvTrks"   ,&fTPhotHasConvTrks   ,"PhoHasConvTrks[NPhotons]/I");
  fEventTree->Branch("PhoIsInJet"       ,&fTPhotIsInJet       ,"PhoIsInJet[NPhotons]/I");
  fEventTree->Branch("PhoIsElDupl"      ,&fTPhotDupEl         ,"PhoIsElDupl[NPhotons]/I");
  fEventTree->Branch("PhoSharedPx"      ,&fTPhotSharedPx      ,"PhoSharedPx[NPhotons]/D");
  fEventTree->Branch("PhoSharedPy"      ,&fTPhotSharedPy      ,"PhoSharedPy[NPhotons]/D");
  fEventTree->Branch("PhoSharedPz"      ,&fTPhotSharedPz      ,"PhoSharedPz[NPhotons]/D");
  fEventTree->Branch("PhoSharedEnergy"  ,&fTPhotSharedEnergy  ,"PhoSharedEnergy[NPhotons]/D");
  fEventTree->Branch("PhoScSeedSeverity",&fTPhotScSeedSeverity,"PhoScSeedSeverity[NPhotons]/I");
  fEventTree->Branch("PhoE1OverE9"      ,&fTPhotE1OverE9      ,"PhoE1OverE9[NPhotons]/D");
  fEventTree->Branch("PhoS4OverS1"      ,&fTPhotS4OverS1      ,"PhoS4OverS1[NPhotons]/D");


  // Jets:
  fEventTree->Branch("NJets"          ,&fTnjets          ,"NJets/I");
  fEventTree->Branch("NJetsTot"       ,&fTnjetstot       ,"NJetsTot/I");
  fEventTree->Branch("JGood"          ,&fTgoodjet        ,"JGood[NJets]/I");
  fEventTree->Branch("JPx"            ,&fTjpx            ,"JPx[NJets]/D");
  fEventTree->Branch("JPy"            ,&fTjpy            ,"JPy[NJets]/D");
  fEventTree->Branch("JPz"            ,&fTjpz            ,"JPz[NJets]/D");
  fEventTree->Branch("JPt"            ,&fTjpt            ,"JPt[NJets]/D");
  fEventTree->Branch("JE"             ,&fTje             ,"JE[NJets]/D");
  fEventTree->Branch("JEt"            ,&fTjet            ,"JEt[NJets]/D");
  fEventTree->Branch("JEta"           ,&fTjeta           ,"JEta[NJets]/D");
  fEventTree->Branch("JPhi"           ,&fTjphi           ,"JPhi[NJets]/D");
  fEventTree->Branch("JEMfrac"        ,&fTjemfrac        ,"JEMfrac[NJets]/D");
  fEventTree->Branch("JNConstituents" ,&fTjNconstituents ,"JNConstituents[NJets]/I");
  fEventTree->Branch("JID_HPD"        ,&fTjID_HPD        ,"JID_HPD[NJets]/D");
  fEventTree->Branch("JID_RBX"        ,&fTjID_RBX        ,"JID_RBX[NJets]/D");
  fEventTree->Branch("JID_n90Hits"    ,&fTjID_n90Hits    ,"JID_n90Hits[NJets]/D");
  fEventTree->Branch("JID_resEMF"     ,&fTjID_resEMF     ,"JID_resEMF[NJets]/D");
  fEventTree->Branch("JID_HCALTow"    ,&fTjID_HCALTow    ,"JID_HCALTow[NJets]/D");
  fEventTree->Branch("JID_ECALTow"    ,&fTjID_ECALTow    ,"JID_ECALTow[NJets]/D");
  fEventTree->Branch("JEtaEMrms"      ,&fTJEtaEMrms      ,"JEtaEMrms[NJets]/D");
  fEventTree->Branch("JEtaHADrms"     ,&fTJEtaHADrms     ,"JEtaHADrms[NJets]/D");
  fEventTree->Branch("JPhiEMrms"      ,&fTJPhiEMrms      ,"JPhiEMrms[NJets]/D");
  fEventTree->Branch("JPhiHADrms"     ,&fTJPhiHADrms     ,"JPhiHADrms[NJets]/D");
  fEventTree->Branch("JbTagProb"      ,&fTjbTagProb      ,"JbTagProb[NJets]/D");
  fEventTree->Branch("JChfrac"        ,&fTjChfrac        ,"JChfrac[NJets]/D");
  fEventTree->Branch("JEFracHadronic" ,&fTjEfracHadr     ,"JEFracHadronic[NJets]/D");
  fEventTree->Branch("JMass"          ,&fTjMass          ,"JMass[NJets]/D");
  fEventTree->Branch("JNAssoTracks"   ,&fTjnAssoTracks   ,"JNAssoTracks[NJets]/I");
  fEventTree->Branch("Jtrk1px"        ,&fTjtrk1px        ,"Jtrk1px[NJets]/D");
  fEventTree->Branch("Jtrk1py"        ,&fTjtrk1py        ,"Jtrk1py[NJets]/D");
  fEventTree->Branch("Jtrk1pz"        ,&fTjtrk1pz        ,"Jtrk1pz[NJets]/D");
  fEventTree->Branch("Jtrk2px"        ,&fTjtrk2px        ,"Jtrk2px[NJets]/D");
  fEventTree->Branch("Jtrk2py"        ,&fTjtrk2py        ,"Jtrk2py[NJets]/D");
  fEventTree->Branch("Jtrk2pz"        ,&fTjtrk2pz        ,"Jtrk2pz[NJets]/D");
  fEventTree->Branch("Jtrk3px"        ,&fTjtrk3px        ,"Jtrk3px[NJets]/D");
  fEventTree->Branch("Jtrk3py"        ,&fTjtrk3py        ,"Jtrk3py[NJets]/D");
  fEventTree->Branch("Jtrk3pz"        ,&fTjtrk3pz        ,"Jtrk3pz[NJets]/D");
  fEventTree->Branch("JEcorr"         ,&fTjEcorr         ,"JEcorr[NJets]/D");
  fEventTree->Branch("JeMinDR"        ,&fTjeMinDR        ,"JeMinDR[NJets]/D");
  fEventTree->Branch("JVtxx"          ,&fTjetVtxx        ,"JVtxx[NJets]/D");
  fEventTree->Branch("JVtxy"          ,&fTjetVtxy        ,"JVtxy[NJets]/D");
  fEventTree->Branch("JVtxz"          ,&fTjetVtxz        ,"JVtxz[NJets]/D");
  fEventTree->Branch("JVtxExx"        ,&fTjetVtxExx      ,"JVtxExx[NJets]/D");
  fEventTree->Branch("JVtxEyx"        ,&fTjetVtxEyx      ,"JVtxEyx[NJets]/D");
  fEventTree->Branch("JVtxEyy"        ,&fTjetVtxEyy      ,"JVtxEyy[NJets]/D");
  fEventTree->Branch("JVtxEzy"        ,&fTjetVtxEzy      ,"JVtxEzy[NJets]/D");
  fEventTree->Branch("JVtxEzz"        ,&fTjetVtxEzz      ,"JVtxEzz[NJets]/D");
  fEventTree->Branch("JVtxEzx"        ,&fTjetVtxEzx      ,"JVtxEzx[NJets]/D");
  fEventTree->Branch("JVtxNChi2"      ,&fTjetVtxNChi2    ,"JVtxNChi2[NJets]/D");
  fEventTree->Branch("JGenPt"         ,&fTjetGenPt       ,"JGenPt[NJets]/D");
  fEventTree->Branch("JGenEta"        ,&fTjetGenEta      ,"JGenEta[NJets]/D");
  fEventTree->Branch("JGenPhi"        ,&fTjetGenPhi      ,"JGenPhi[NJets]/D");
  fEventTree->Branch("JGenE"          ,&fTjetGenE        ,"JGenE[NJets]/D");
  fEventTree->Branch("JGenEmE"        ,&fTjetGenemE      ,"JGenEmE[NJets]/D");
  fEventTree->Branch("JGenHadE"       ,&fTjetGenhadE     ,"JGenHadE[NJets]/D");
  fEventTree->Branch("JGenInvE"       ,&fTjetGeninvE     ,"JGenInvE[NJets]/D");

  for ( std::vector<JetFiller*>::iterator it = jetFillers.begin(); 
        it != jetFillers.end(); ++it )
    (*it)->createBranches();

  // Tracks:
  fEventTree->Branch("NTracks"        ,&fTntracks      ,"NTracks/I");
  fEventTree->Branch("NTracksTot"     ,&fTntrackstot   ,"NTracksTot/I");
  fEventTree->Branch("TrkGood"        ,&fTgoodtrk      ,"TrkGood[NTracks]/I");
  fEventTree->Branch("TrkPt"          ,&fTtrkpt        ,"TrkPt[NTracks]/D");
  fEventTree->Branch("TrkEta"         ,&fTtrketa       ,"TrkEta[NTracks]/D");
  fEventTree->Branch("TrkPhi"         ,&fTtrkphi       ,"TrkPhi[NTracks]/D");
  fEventTree->Branch("TrkNChi2"       ,&fTtrknchi2     ,"TrkNChi2[NTracks]/D");
  fEventTree->Branch("TrkNHits"       ,&fTtrknhits     ,"TrkNHits[NTracks]/D");
  fEventTree->Branch("TrkPtSumx"      ,&fTTrkPtSumx      ,"TrkPtSumx/D");
  fEventTree->Branch("TrkPtSumy"      ,&fTTrkPtSumy      ,"TrkPtSumy/D");
  fEventTree->Branch("TrkPtSum"       ,&fTTrkPtSum       ,"TrkPtSum/D");
  fEventTree->Branch("TrkPtSumPhi"    ,&fTTrkPtSumphi    ,"TrkPtSumPhi/D");

  // MET:
  fEventTree->Branch("SumEt"              ,&fTSumEt               ,"SumEt/D");
  fEventTree->Branch("ECALSumEt"          ,&fTECALSumEt           ,"ECALSumEt/D");
  fEventTree->Branch("HCALSumEt"          ,&fTHCALSumEt           ,"HCALSumEt/D");
  fEventTree->Branch("ECALEsumx"          ,&fTECALEsumx           ,"ECALEsumx/D");
  fEventTree->Branch("ECALEsumy"          ,&fTECALEsumy           ,"ECALEsumy/D");
  fEventTree->Branch("ECALEsumz"          ,&fTECALEsumz           ,"ECALEsumz/D");
  fEventTree->Branch("ECALMET"            ,&fTECALMET             ,"ECALMET/D");
  fEventTree->Branch("ECALMETPhi"         ,&fTECALMETphi          ,"ECALMETPhi/D");
  fEventTree->Branch("ECALMETEta"         ,&fTECALMETeta          ,"ECALMETEta/D");
  fEventTree->Branch("HCALEsumx"          ,&fTHCALEsumx           ,"HCALEsumx/D");
  fEventTree->Branch("HCALEsumy"          ,&fTHCALEsumy           ,"HCALEsumy/D");
  fEventTree->Branch("HCALEsumz"          ,&fTHCALEsumz           ,"HCALEsumz/D");
  fEventTree->Branch("HCALMET"            ,&fTHCALMET             ,"HCALMET/D");
  fEventTree->Branch("HCALMETPhi"         ,&fTHCALMETphi          ,"HCALMETPhi/D");
  fEventTree->Branch("HCALMETeta"         ,&fTHCALMETeta          ,"HCALMETEta/D");
  fEventTree->Branch("RawMET"             ,&fTRawMET              ,"RawMET/D");
  fEventTree->Branch("RawMETpx"           ,&fTRawMETpx            ,"RawMETpx/D");
  fEventTree->Branch("RawMETpy"           ,&fTRawMETpy            ,"RawMETpy/D");
  fEventTree->Branch("RawMETphi"          ,&fTRawMETphi           ,"RawMETphi/D");
  fEventTree->Branch("RawMETemEtFrac"     ,&fTRawMETemEtFrac      ,"RawMETemEtFrac/D");
  fEventTree->Branch("RawMETemEtInEB"     ,&fTRawMETemEtInEB      ,"RawMETemEtInEB/D");
  fEventTree->Branch("RawMETemEtInEE"     ,&fTRawMETemEtInEE      ,"RawMETemEtInEE/D");
  fEventTree->Branch("RawMETemEtInHF"     ,&fTRawMETemEtInHF      ,"RawMETemEtInHF/D");
  fEventTree->Branch("RawMEThadEtFrac"    ,&fTRawMEThadEtFrac     ,"RawMEThadEtFrac/D");
  fEventTree->Branch("RawMEThadEtInHB"    ,&fTRawMEThadEtInHB     ,"RawMEThadEtInHB/D");
  fEventTree->Branch("RawMEThadEtInHE"    ,&fTRawMEThadEtInHE     ,"RawMEThadEtInHE/D");
  fEventTree->Branch("RawMEThadEtInHF"    ,&fTRawMEThadEtInHF     ,"RawMEThadEtInHF/D");
  fEventTree->Branch("RawMETSignificance" ,&fTRawMETSignificance  ,"RawMETSignificance/D");
  fEventTree->Branch("MuCorrMET"          ,&fTMuCorrMET           ,"MuCorrMET/D");
  fEventTree->Branch("MuCorrMETpx"        ,&fTMuCorrMETpx         ,"MuCorrMETpx/D");
  fEventTree->Branch("MuCorrMETpy"        ,&fTMuCorrMETpy         ,"MuCorrMETpy/D");
  fEventTree->Branch("MuCorrMETphi"       ,&fTMuCorrMETphi        ,"MuCorrMETphi/D");
  fEventTree->Branch("TCMET"              ,&fTTCMET               ,"TCMET/D");
  fEventTree->Branch("TCMETpx"            ,&fTTCMETpx             ,"TCMETpx/D");
  fEventTree->Branch("TCMETpy"            ,&fTTCMETpy             ,"TCMETpy/D");
  fEventTree->Branch("TCMETphi"           ,&fTTCMETphi            ,"TCMETphi/D");
  fEventTree->Branch("TCMETSignificance"  ,&fTTCMETSignificance   ,"TCMETSignificance/D");
  fEventTree->Branch("MuJESCorrMET"       ,&fTMuJESCorrMET        ,"MuJESCorrMET/D");
  fEventTree->Branch("MuJESCorrMETpx"     ,&fTMuJESCorrMETpx      ,"MuJESCorrMETpx/D");
  fEventTree->Branch("MuJESCorrMETpy"     ,&fTMuJESCorrMETpy      ,"MuJESCorrMETpy/D");
  fEventTree->Branch("MuJESCorrMETphi"    ,&fTMuJESCorrMETphi     ,"MuJESCorrMETphi/D");
  fEventTree->Branch("PFMET"              ,&fTPFMET               ,"PFMET/D");
  fEventTree->Branch("PFMETpx"            ,&fTPFMETpx             ,"PFMETpx/D");
  fEventTree->Branch("PFMETpy"            ,&fTPFMETpy             ,"PFMETpy/D");
  fEventTree->Branch("PFMETphi"           ,&fTPFMETphi            ,"PFMETphi/D");
  fEventTree->Branch("PFMETSignificance"  ,&fTPFMETSignificance   ,"PFMETSignificance/D");
  fEventTree->Branch("METR12"             ,&fTMETR12              ,"METR12/D");
  fEventTree->Branch("METR21"             ,&fTMETR21              ,"METR21/D");
}

// Method called once before each run
void NTupleProducer::beginRun(const edm::Run& r, const edm::EventSetup&){
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
  edm::LogVerbatim("NTP") << " ---------------------------------------------------";
  edm::LogVerbatim("NTP") << " ==> NTupleProducer::endJob() ...";
  edm::LogVerbatim("NTP") << "  Total number of processed Events: " << fNTotEvents;
  edm::LogVerbatim("NTP") << "  Number of times Tree was filled:  " << fNFillTree;
  edm::LogVerbatim("NTP") << " ---------------------------------------------------";
  //violatemp


}

// Method to reset the TTree variables for each event
void NTupleProducer::resetTree(){
  resetInt(fTHLTres, gMaxhltbits);
  resetInt(fTL1physres, gMaxl1physbits);
  resetInt(fTL1techres, gMaxl1techbits);
  fTrunnumber   = -999;
  fTeventnumber = -999;
  fTlumisection = -999;
  fTpthat       = -999.99;
  fTsigprocid   = -999;
  fTpdfscalePDF = -999.99;
  fTpdfid1      = -999;
  fTpdfid2      = -999;
  fTpdfx1       = -999.99;
  fTpdfx2       = -999.99;
  fTpdfxPDF1    = -999.99;
  fTpdfxPDF2    = -999.99;
	
  fTweight      = -999.99;
  fTgoodvtx     = -999;
  fTprimvtxx    = -999.99;
  fTprimvtxy    = -999.99;
  fTprimvtxz    = -999.99;
  fTprimvtxrho  = -999.99;
  fTprimvtxxE   = -999.99;
  fTprimvtxyE   = -999.99;
  fTprimvtxzE   = -999.99;
  fTpvtxznchi2  = -999.99;
  fTpvtxisfake  = -999;
  fTpvtxndof    = -999;
  fTpvtxptsum   = -999.99;
  fTbeamspotx   = -999.99;
  fTbeamspoty   = -999.99;
  fTbeamspotz   = -999.99;
  fTNCaloTowers = -999;
  fTHBHENoiseFlag = -999;

  fTgoodevent        = 0;
  fTflagmaxmuexc     = 0;
  fTflagmaxelexc     = 0;
  fTflagmaxujetexc = 0;
  fTflagmaxjetexc    = 0;
  fTflagmaxtrkexc    = 0;
  fTflagmaxphoexc    = 0;

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

  resetInt(fTgoodmu, gMaxnmus);
  resetInt(fTmuIsIso, gMaxnmus);
  resetInt(fTmuIsGM, gMaxnmus);
  resetInt(fTmuIsTM, gMaxnmus);
  resetDouble(fTmupx, gMaxnmus);
  resetDouble(fTmupy, gMaxnmus);
  resetDouble(fTmupz, gMaxnmus);
  resetDouble(fTmupt, gMaxnmus);
  resetDouble(fTmuptE, gMaxnmus);
  resetDouble(fTmue, gMaxnmus);
  resetDouble(fTmuet, gMaxnmus);
  resetDouble(fTmueta, gMaxnmus);
  resetDouble(fTmuphi, gMaxnmus);
  resetInt(fTmucharge, gMaxnmus);
  resetDouble(fTmuiso, gMaxnmus);
  resetDouble(fTmuIso03sumPt, gMaxnmus);
  resetDouble(fTmuIso03emEt, gMaxnmus);
  resetDouble(fTmuIso03hadEt, gMaxnmus);
  resetDouble(fTmuIso03emVetoEt, gMaxnmus);
  resetDouble(fTmuIso03hadVetoEt, gMaxnmus);
  resetDouble(fTmuIso05sumPt, gMaxnmus);
  resetDouble(fTmuIso05emEt, gMaxnmus);
  resetDouble(fTmuIso05hadEt, gMaxnmus);
  resetDouble(fTmueecal, gMaxnmus);
  resetDouble(fTmuehcal, gMaxnmus);
  resetDouble(fTmud0bs, gMaxnmus);
  resetDouble(fTmud0pv, gMaxnmus);
  resetDouble(fTmud0E, gMaxnmus);
  resetDouble(fTmudzbs, gMaxnmus);
  resetDouble(fTmudzpv, gMaxnmus);
  resetDouble(fTmudzE, gMaxnmus);
  resetDouble(fTmunchi2, gMaxnmus);
  resetInt(fTmunglhits, gMaxnmus);
  resetInt(fTmunmuhits, gMaxnmus);
  resetInt(fTmuntkhits, gMaxnmus);
  resetDouble(fTmuinntknchi2, gMaxnmus);
  resetInt(fTmunmatches, gMaxnmus);
  resetInt(fTmunchambers, gMaxnmus);
  resetDouble(fTmucalocomp, gMaxnmus);
  resetDouble(fTmusegmcomp, gMaxnmus);
  resetDouble(fTmuoutmomx, gMaxnmus);
  resetDouble(fTmuoutmomy, gMaxnmus);
  resetDouble(fTmuoutmomz, gMaxnmus);
  resetDouble(fTmuoutmomphi, gMaxnmus);
  resetDouble(fTmuoutmometa, gMaxnmus);
  resetDouble(fTmuoutmomtheta, gMaxnmus);
  resetDouble(fTmuoutposrad, gMaxnmus);
  resetDouble(fTmuoutposx, gMaxnmus);
  resetDouble(fTmuoutposy, gMaxnmus);
  resetDouble(fTmuoutposz, gMaxnmus);

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
  resetInt(fTGenMuCharge, gMaxnmus);
  resetDouble(fTGenMuPt, gMaxnmus);
  resetDouble(fTGenMuEta, gMaxnmus);
  resetDouble(fTGenMuPhi, gMaxnmus);
  resetDouble(fTGenMuE, gMaxnmus);
  resetInt(fTGenMuMId, gMaxnmus);
  resetInt(fTGenMuMStatus, gMaxnmus);
  resetInt(fTGenMuMCharge, gMaxnmus);
  resetDouble(fTGenMuMPt, gMaxnmus);
  resetDouble(fTGenMuMEta, gMaxnmus);
  resetDouble(fTGenMuMPhi, gMaxnmus);
  resetDouble(fTGenMuME, gMaxnmus);
  resetInt(fTGenMuGMId, gMaxnmus);
  resetInt(fTGenMuGMStatus, gMaxnmus);
  resetInt(fTGenMuGMCharge, gMaxnmus);
  resetDouble(fTGenMuGMPt, gMaxnmus);
  resetDouble(fTGenMuGMEta, gMaxnmus);
  resetDouble(fTGenMuGMPhi, gMaxnmus);
  resetDouble(fTGenMuGME, gMaxnmus);

  resetInt(fTgoodel, gMaxneles);
  resetInt(fTeIsIso, gMaxneles);
  resetInt(fTeChargeMisIDProb, gMaxneles);
  resetDouble(fTepx, gMaxneles);
  resetDouble(fTepy, gMaxneles);
  resetDouble(fTepz, gMaxneles);
  resetDouble(fTee, gMaxneles);
  resetDouble(fTeet, gMaxneles);
  resetDouble(fTept, gMaxneles);
  resetDouble(fTeptE, gMaxneles);
  resetDouble(fTeeta, gMaxneles);
  resetDouble(fTephi, gMaxneles);
  resetDouble(fTed0bs, gMaxneles);
  resetDouble(fTed0pv, gMaxneles);
  resetDouble(fTed0E, gMaxneles);
  resetDouble(fTedzbs, gMaxneles);
  resetDouble(fTedzpv, gMaxneles);
  resetDouble(fTedzE, gMaxneles);
  resetDouble(fTenchi2, gMaxneles);
  resetDouble(fTeiso03, gMaxneles);
  resetDouble(fTeiso04, gMaxneles);
  resetDouble(fTdr03tksumpt, gMaxneles);
  resetDouble(fTdr04tksumpt, gMaxneles);
  resetDouble(fTdr03ecalrechitsumet, gMaxneles);
  resetDouble(fTdr04ecalrechitsumet, gMaxneles);
  resetDouble(fTdr03hcaltowersumet, gMaxneles);
  resetDouble(fTdr04hcaltowersumet, gMaxneles);
  resetDouble(fTetheta, gMaxneles);
  resetDouble(fTesceta, gMaxneles);
  resetInt(fTecharge, gMaxneles);
  resetInt(fTeCInfoIsGsfCtfCons, gMaxneles);
  resetInt(fTeCInfoIsGsfCtfScPixCons, gMaxneles);
  resetInt(fTeCInfoIsGsfScPixCons, gMaxneles);
  resetInt(fTeCInfoScPixCharge, gMaxneles);
  resetDouble(fTeClosestCtfTrackpt, gMaxneles);
  resetDouble(fTeClosestCtfTracketa, gMaxneles);
  resetDouble(fTeClosestCtfTrackphi, gMaxneles);
  resetInt(fTeClosestCtfTrackcharge, gMaxneles);
  resetInt(fTeInGap, gMaxneles);
  resetInt(fTeEcalDriven, gMaxneles);
  resetInt(fTeTrackerDriven, gMaxneles);
  resetInt(fTeBasicClustersSize, gMaxneles);
  resetDouble(fTefbrem, gMaxneles);
  resetDouble(fTeHcalOverEcal, gMaxneles);
  resetDouble(fTeE1x5, gMaxneles);
  resetDouble(fTeE5x5, gMaxneles);
  resetDouble(fTeE2x5Max, gMaxneles);
  resetDouble(fTeSigmaIetaIeta, gMaxneles);
  resetDouble(fTeDeltaPhiSeedClusterAtCalo, gMaxneles);
  resetDouble(fTeDeltaEtaSeedClusterAtCalo, gMaxneles);
  resetDouble(fTeDeltaPhiSuperClusterAtVtx, gMaxneles);
  resetDouble(fTeDeltaEtaSuperClusterAtVtx, gMaxneles);
  resetDouble(fTecaloenergy, gMaxneles);
  resetDouble(fTetrkmomatvtx, gMaxneles);
  resetDouble(fTeESuperClusterOverP, gMaxneles);
  resetInt(fTeNumberOfMissingInnerHits, gMaxneles);
  resetInt(fTeIsInJet, gMaxneles);
  resetDouble(fTeSharedPx, gMaxneles);
  resetDouble(fTeSharedPy, gMaxneles);
  resetDouble(fTeSharedPz, gMaxneles);
  resetDouble(fTeSharedEnergy, gMaxneles);
  resetInt(fTeDupEl, gMaxneles);

  resetDouble(fTeConvPartTrackDist, gMaxneles);
  resetDouble(fTeConvPartTrackDCot, gMaxneles);
  resetDouble(fTeConvPartTrackPt, gMaxneles);
  resetDouble(fTeConvPartTrackEta, gMaxneles);
  resetDouble(fTeConvPartTrackPhi, gMaxneles);
  resetDouble(fTeConvPartTrackCharge, gMaxneles);

  resetInt(fTeScSeedSeverity, gMaxneles);
  resetDouble(fTeS4OverS1, gMaxneles);
  resetDouble(fTeE1OverE9, gMaxneles);

  resetInt(fTeIDTight, gMaxneles);
  resetInt(fTeIDLoose, gMaxneles);
  resetInt(fTeIDRobustTight, gMaxneles);
  resetInt(fTeIDRobustLoose, gMaxneles);

  resetInt(fTGenElId, gMaxneles);
  resetInt(fTGenElStatus, gMaxneles);
  resetInt(fTGenElCharge, gMaxneles);
  resetDouble(fTGenElPt, gMaxneles);
  resetDouble(fTGenElEta, gMaxneles);
  resetDouble(fTGenElPhi, gMaxneles);
  resetDouble(fTGenElE, gMaxneles);
  resetInt(fTGenElMId, gMaxneles);
  resetInt(fTGenElMStatus, gMaxneles);
  resetInt(fTGenElMCharge, gMaxneles);
  resetDouble(fTGenElMPt, gMaxneles);
  resetDouble(fTGenElMEta, gMaxneles);
  resetDouble(fTGenElMPhi, gMaxneles);
  resetDouble(fTGenElME, gMaxneles);
  resetInt(fTGenElGMId, gMaxneles);
  resetInt(fTGenElGMStatus, gMaxneles);
  resetInt(fTGenElGMCharge, gMaxneles);
  resetDouble(fTGenElGMPt, gMaxneles);
  resetDouble(fTGenElGMEta, gMaxneles);
  resetDouble(fTGenElGMPhi, gMaxneles);
  resetDouble(fTGenElGME, gMaxneles);

  resetInt(fTgoodjet, gMaxnjets);
  resetDouble(fTjpx,  gMaxnjets);
  resetDouble(fTjpy,  gMaxnjets);
  resetDouble(fTjpz,  gMaxnjets);
  resetDouble(fTje,   gMaxnjets);
  resetDouble(fTjet,  gMaxnjets);
  resetDouble(fTjpt,  gMaxnjets);
  resetDouble(fTjeta, gMaxnjets);
  resetDouble(fTjphi, gMaxnjets);
  resetDouble(fTjemfrac, gMaxnjets);
  resetDouble(fTjID_HPD, gMaxnjets);
  resetDouble(fTjID_RBX, gMaxnjets);
  resetDouble(fTjID_n90Hits, gMaxnjets);
  resetDouble(fTjID_resEMF,  gMaxnjets);
  resetDouble(fTjID_HCALTow, gMaxnjets);
  resetDouble(fTjID_ECALTow, gMaxnjets);
  resetDouble(fTjbTagProb, gMaxnjets);
  resetDouble(fTjChfrac,   gMaxnjets);
  resetDouble(fTjEfracHadr, gMaxnjets);
  resetDouble(fTjMass,   gMaxnjets);
  resetInt(fTjnAssoTracks, gMaxnjets);
  resetDouble(fTjtrk1px, gMaxnjets);
  resetDouble(fTjtrk1py, gMaxnjets);
  resetDouble(fTjtrk1pz, gMaxnjets);
  resetDouble(fTjtrk2px, gMaxnjets);
  resetDouble(fTjtrk2py, gMaxnjets);
  resetDouble(fTjtrk2pz, gMaxnjets);
  resetDouble(fTjtrk3px, gMaxnjets);
  resetDouble(fTjtrk3py, gMaxnjets);
  resetDouble(fTjtrk3pz, gMaxnjets);
  resetDouble(fTjeMinDR, gMaxnjets);
  resetDouble(fTjetVtxx, gMaxnjets);
  resetDouble(fTjetVtxy, gMaxnjets);
  resetDouble(fTjetVtxz, gMaxnjets);
  resetDouble(fTjetVtxExx, gMaxnjets);
  resetDouble(fTjetVtxEyx, gMaxnjets);
  resetDouble(fTjetVtxEyy, gMaxnjets);
  resetDouble(fTjetVtxEzy, gMaxnjets);
  resetDouble(fTjetVtxEzz, gMaxnjets);
  resetDouble(fTjetVtxEzx, gMaxnjets);
  resetDouble(fTjetVtxNChi2, gMaxnjets);
  resetInt(fTjNconstituents, gMaxnjets);
  resetDouble(fTjetGenPt, gMaxnjets);
  resetDouble(fTjetGenEta, gMaxnjets);
  resetDouble(fTjetGenPhi, gMaxnjets);
  resetDouble(fTjetGenE, gMaxnjets);
  resetDouble(fTjetGenemE, gMaxnjets);
  resetDouble(fTjetGenhadE, gMaxnjets);
  resetDouble(fTjetGeninvE, gMaxnjets);

  resetDouble(fJUNC_px_match, gMaxnjets);
  resetDouble(fJUNC_py_match, gMaxnjets);
  resetDouble(fJUNC_pz_match, gMaxnjets);

  resetInt(fTgoodtrk,  gMaxntrks);
  resetDouble(fTtrkpt, gMaxntrks);
  resetDouble(fTtrketa, gMaxntrks);
  resetDouble(fTtrkphi, gMaxntrks);
  resetDouble(fTtrknchi2, gMaxntrks);
  resetDouble(fTtrknhits, gMaxntrks);

  resetDouble(fTPhotPt,gMaxnphos);
  resetDouble(fTPhotPx,gMaxnphos);
  resetDouble(fTPhotPy,gMaxnphos);
  resetDouble(fTPhotPz,gMaxnphos);
  resetDouble(fTPhotEta,gMaxnphos);
  resetDouble(fTPhotPhi,gMaxnphos);
  resetDouble(fTPhotEnergy,gMaxnphos);
  resetDouble(fTPhotIso03Ecal);
  resetDouble(fTPhotIso03Hcal);
  resetDouble(fTPhotIso03TrkSolid);
  resetDouble(fTPhotIso03TrkHollow);
  resetDouble(fTPhotIso03);
  resetDouble(fTPhotIso04Ecal);
  resetDouble(fTPhotIso04Hcal);
  resetDouble(fTPhotIso04TrkSolid);
  resetDouble(fTPhotIso04TrkHollow);
  resetDouble(fTPhotIso04);
  resetDouble(fTPhotcaloPosX,gMaxnphos);
  resetDouble(fTPhotcaloPosY,gMaxnphos);
  resetDouble(fTPhotcaloPosZ,gMaxnphos);
  resetDouble(fTPhotHoverE,gMaxnphos);
  resetDouble(fTPhotH1overE,gMaxnphos);
  resetDouble(fTPhotH2overE,gMaxnphos);
  resetDouble(fTPhotSigmaIetaIeta,gMaxnphos);
  resetInt(fTPhotHasPixSeed,gMaxnphos);
  resetInt(fTPhotHasConvTrks,gMaxnphos);
  resetInt(fTgoodphoton,gMaxnphos);
  resetInt(fTPhotIsIso,gMaxnphos);
  resetInt(fTPhotIsInJet,gMaxnphos);
  resetInt(fTPhotDupEl,gMaxnphos);
  resetDouble(fTPhotSharedPx, gMaxnphos);
  resetDouble(fTPhotSharedPy, gMaxnphos);
  resetDouble(fTPhotSharedPz, gMaxnphos);
  resetDouble(fTPhotSharedEnergy, gMaxnphos);

  resetInt(fTPhotScSeedSeverity, gMaxnphos);
  resetDouble(fTPhotS4OverS1, gMaxnphos);
  resetDouble(fTPhotE1OverE9, gMaxnphos);

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
  fTMuCorrMET          = -999.99;
  fTMuCorrMETpx        = -999.99;
  fTMuCorrMETpy        = -999.99;
  fTMuCorrMETphi       = -999.99;
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
  fTMETR12             = -999.99;
  fTMETR21             = -999.99;
}

// Method for matching of reco candidates
vector<const reco::GenParticle*> NTupleProducer::matchRecoCand(const reco::RecoCandidate *Cand, const edm::Event& iEvent){
  const GenParticle *GenCand = NULL;
  const GenParticle *GenMom  = NULL;
  const GenParticle *GenGMom = NULL;
  vector<const GenParticle*> res;
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
  bool matched = false;

  // Try to match the reco candidate to a generator object
  double mindr(999.99);
  for(gpart = genparts->begin(); gpart != genparts->end(); gpart++){
    if( gpart->status() != 1 ) continue;

    // Restrict to cone of 0.1 in DR around candidate
    double dr = reco::deltaR(gpart->eta(), gpart->phi(), Cand->eta(), Cand->phi());
    if(dr > 0.1) continue;

    // Select same charge
    // if(gpart->charge() != Cand->charge()) continue;

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
    id  = GenCand->pdgId();

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

    res.push_back(GenCand);
    res.push_back(GenMom);
    res.push_back(GenGMom);
  }
  else{
    res.push_back(GenCand);
    res.push_back(GenMom);
    res.push_back(GenGMom);
  }
  return res;
}

// Method for matching of jets
const reco::GenJet* NTupleProducer::matchJet(const reco::Jet* jet, const edm::Event& iEvent){
  const reco::GenJet *genjet = NULL;
  if(fIsRealData){
    edm::LogWarning("NTP") << "@SUB=matchJet"
                           << "Trying to access generator info on real data...";
    return genjet;
  }

  edm::Handle<GenJetCollection> genjets;
  iEvent.getByLabel(fGenJetTag, genjets);
  GenJetCollection::const_iterator gjet;
  bool matched = false;

  // Try to match the reco jet to a generator jet
  double mindr(999.99);
  for(gjet = genjets->begin(); gjet != genjets->end(); gjet++){

    // Restrict to cone of 0.1 in DR around candidate
    double dr = reco::deltaR(gjet->eta(), gjet->phi(), jet->eta(), jet->phi());
    if(dr > 0.3) continue;

    // Restrict to pt match within a factor of 2
    double ndpt = fabs(gjet->pt() - jet->pt())/gjet->pt();
    if(ndpt > 2.) continue;

    // Minimize DeltaR
    if(dr > mindr) continue;
    mindr = dr;

    matched = true;
    genjet = &(*gjet);
  }
  return genjet;
}

/////////////////////////////////////////////////////////////////////////////////////////////
// Electron Conversion Information
reco::TrackRef NTupleProducer::getConversionPartnerTrack(const reco::GsfElectron& gsfElectron, const edm::Handle<reco::TrackCollection>& track_h, const float bFieldAtOrigin, double& Dist, double& DCot, const float maxAbsDist, const float maxAbsDCot, const float minFracSharedHits){
  using namespace edm;
  using namespace reco;
  const reco::TrackRef el_ctftrack = gsfElectron.closestCtfTrackRef();
  const TrackCollection *ctftracks = track_h.product();

  const reco::Track* el_track = getElectronTrack(gsfElectron, minFracSharedHits);
  int ctfidx = -999;
  int el_q   = el_track->charge();
  LorentzVector el_tk_p4(el_track->px(), el_track->py(), el_track->pz(), el_track->p());
  double el_d0 = el_track->d0();

  if(el_ctftrack.isNonnull() && gsfElectron.shFracInnerHits() > minFracSharedHits)
    ctfidx = static_cast<int>(el_ctftrack.key());

  int tk_i = 0;
  double mindR = 999;

  //make a null Track Ref
  TrackRef ctfTrackRef = TrackRef() ;

  for(TrackCollection::const_iterator tk = ctftracks->begin();
      tk != ctftracks->end(); tk++, tk_i++) {
    //if the general Track is the same one as made by the electron, skip it
    if((tk_i == ctfidx)  &&  (gsfElectron.shFracInnerHits() > minFracSharedHits))
      continue;


    LorentzVector tk_p4 = LorentzVector(tk->px(), tk->py(),
                                        tk->pz(), tk->p());

    //look only in a cone of 0.3
    double dR = deltaR(el_tk_p4, tk_p4);
    if(dR > 0.3)
      continue;

    int tk_q = tk->charge();
    double tk_d0 = tk->d0();

    //the electron and track must be opposite charge
    if(tk_q + el_q != 0)
      continue;

    std::pair<double, double> convInfo =  getConversionInfo(el_tk_p4, el_q, el_d0,
                                                            tk_p4, tk_q, tk_d0,
                                                            bFieldAtOrigin);

    double dist = convInfo.first;
    double dcot = convInfo.second;

    if(fabs(dist) < maxAbsDist && fabs(dcot) < maxAbsDCot && dR < mindR) {
      ctfTrackRef = reco::TrackRef(track_h, tk_i);
      mindR = dR;
      Dist = dist ;
      DCot = dcot ;
    }

  }//track loop

  return ctfTrackRef;
}

const reco::Track* NTupleProducer::getElectronTrack(const reco::GsfElectron& electron, const float minFracSharedHits) {
  if(electron.closestCtfTrackRef().isNonnull() &&
     electron.shFracInnerHits() > minFracSharedHits)
    return (const reco::Track*)electron.closestCtfTrackRef().get();
  return (const reco::Track*)(electron.gsfTrack().get());
}

std::pair<double, double> NTupleProducer::getConversionInfo(LorentzVector trk1_p4, int trk1_q, float trk1_d0, LorentzVector trk2_p4, int trk2_q, float trk2_d0, float bFieldAtOrigin) {

  double tk1Curvature = -0.3*bFieldAtOrigin*(trk1_q/trk1_p4.pt())/100.;
  double rTk1 = fabs(1./tk1Curvature);
  double xTk1 = (1./tk1Curvature - trk1_d0)*cos(trk1_p4.phi());
  double yTk1 = (1./tk1Curvature - trk1_d0)*sin(trk1_p4.phi());

  double tk2Curvature = -0.3*bFieldAtOrigin*(trk2_q/trk2_p4.pt())/100.;
  double rTk2 = fabs(1./tk2Curvature);
  double xTk2 = (1./tk2Curvature - trk2_d0)*cos(trk2_p4.phi());
  double yTk2 = (1./tk2Curvature - trk2_d0)*sin(trk2_p4.phi());

  double dist = sqrt(pow(xTk1-xTk2, 2) + pow(yTk1-yTk2 , 2));
  dist = dist - (rTk1 + rTk2);

  double dcot = 1/tan(trk1_p4.theta()) - 1/tan(trk2_p4.theta());

  return std::make_pair(dist, dcot);

}

/////////////////////////////////////////////////////////////////////////////////////////////
// Cleaning methods
void NTupleProducer::ElectronDuplicate(vector<const SuperCluster*> elecPtr, vector<const GsfTrack*> trckPtr) {
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
        fTeDupEl[i] = j;
        fTeDupEl[j] = i;
        break;
      }
    }
  }

  return;
}

void NTupleProducer::PhotonElectronDuplicate(vector<const SuperCluster*> elecPtr, vector<const SuperCluster*> phoPtr) {
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
        fTPhotDupEl[i] = j;
        break;
      }
    }
  }

  return;
}

void NTupleProducer::ElJetOverlap(vector<const Jet*> jets, vector<const SuperCluster*> electrons, edm::Handle<CaloTowerCollection> calotowers){
  // checks for jets made from electrons
  // jetIndex and elecIndex contain the indices of the selected jets and
  //   electrons in the Views
  // (electrons and jets should be filled in the ntuple before checking)   
  if (fTnjets <= 0) return;
  if (fTneles <= 0) return;

  vector<CaloTowerPtr> jetCaloRefs;

  // loop over the jets
  for (int i = 0; i < fTnjets; ++i) {

    // Collect the CaloTowers detIds for the jet 
    if (!fIsPat) {
      const CaloJet* theJet = static_cast<const CaloJet*>(&(*jets[i]));			
      jetCaloRefs = theJet->getCaloConstituents();
    } else {
      const pat::Jet* theJet = static_cast<const pat::Jet*>(&(*jets[i]));
      jetCaloRefs = theJet->getCaloConstituents();
    }

    // loop over the electrons
    for( int j = 0; j < fTneles; ++j ){
      const SuperCluster* theElecSC = electrons[j];

      math::XYZVector sharedP(0., 0., 0.);
      bool isInJet = IsEMObjectInJet(theElecSC, jetCaloRefs, calotowers, &sharedP);
      float sharedE = sqrt(sharedP.X()*sharedP.X() 
                           + sharedP.Y()*sharedP.Y()
                           + sharedP.Z()*sharedP.Z() );
      if (isInJet) {
        fTeIsInJet[j] = i;
        // was: fTeIsInJet[j] = isDuplicate ? 1:0;
        fTeSharedPx[j] = sharedP.X();
        fTeSharedPy[j] = sharedP.Y();
        fTeSharedPz[j] = sharedP.Z();
        fTeSharedEnergy[j] = sharedE;
        break;
      }
    }
  }

  return;
}

void NTupleProducer::PhotonJetOverlap(vector<const Jet*> jets, vector<const SuperCluster*> superclusters, edm::Handle<CaloTowerCollection> calotowers){
  // checks for jets made from photons
  // (photons and jets should be filled in the ntuple before checking)   
  if( fTnjets <= 0 ) return;
  if( fTnphotons <= 0 ) return;

  vector<CaloTowerPtr> jetCaloRefs;

  // loop over the jets
  for( int i = 0; i < fTnjets; ++i ){

    // Collect the CaloTowers detIds for the jet 
    if (!fIsPat) {
      const CaloJet* theJet = static_cast<const CaloJet*>(&(*jets[i]));			
      jetCaloRefs = theJet->getCaloConstituents();
    } else {
      const pat::Jet* theJet = static_cast<const pat::Jet*>(&(*jets[i]));
      jetCaloRefs = theJet->getCaloConstituents();
    }

    // loop over the photons
    for( int j = 0; j < fTnphotons; ++j ){
      const SuperCluster* theSC = superclusters[j];

      math::XYZVector sharedP(0., 0., 0.);
      bool isInJet = IsEMObjectInJet(theSC, jetCaloRefs, calotowers, &sharedP);
      float sharedE = sqrt(sharedP.X()*sharedP.X() 
                           + sharedP.Y()*sharedP.Y()
                           + sharedP.Z()*sharedP.Z() );
      if( isInJet ){
        fTPhotIsInJet[j] = i;
        fTPhotSharedPx[j] = sharedP.X();
        fTPhotSharedPy[j] = sharedP.Y();
        fTPhotSharedPz[j] = sharedP.Z();
        fTPhotSharedEnergy[j] = sharedE;
        break;
      }
    }
  }

  return;
}

bool NTupleProducer::IsEMObjectInJet(const SuperCluster* elecSC, vector<CaloTowerPtr> jetCaloRefs, edm::Handle<CaloTowerCollection> calotowers, math::XYZVector* sharedMomentum){
  // Checks whether an electron or photon is included in the jet energy
  // and if true, it returns the momentum vector shared by the two

  // Define a window in eta,phi for the SuperCluster
  float phimin=0., phimax=0., etamin=0., etamax=0.;
  bool window = EMCaloTowerWindow(elecSC, phimin, phimax, etamin, etamax); 
  if (!window){ return false;}

  // Collect the CaloTowers inside this window, save their detId in a vector
  vector<CaloTowerDetId> eleDetId;
  vector<float> eleTowerEnergy;
  vector<float> eleTowerEta;
  vector<float> eleTowerPhi;
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


  // Define a window in eta,phi for the SuperCluster
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

void NTupleProducer::switchDouble(double &d1, double &d2){
  double temp = d1;
  d1 = d2;
  d2 = temp;
}

void NTupleProducer::switchInt(int &i1, int &i2){
  int temp = i1;
  i1 = i2;
  i2 = temp;
}

void NTupleProducer::resetDouble(double *v, unsigned int size){
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
