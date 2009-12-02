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
// $Id: NTupleProducer.cc,v 1.28 2009/12/01 19:46:14 stiegerb Exp $
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

// Data formats
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "DataFormats/BTauReco/interface/JetTag.h"

#include "DataFormats/METReco/interface/METCollection.h"
#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/METReco/interface/CaloMETCollection.h"
#include "DataFormats/METReco/interface/PFMET.h"

#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/Common/interface/TriggerResults.h"

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Jet.h"

// Interface
#include "DiLeptonAnalysis/NTupleProducer/interface/NTupleProducer.h"


NTupleProducer::NTupleProducer(const edm::ParameterSet& iConfig){
  // Main settings
  fIsRealData = iConfig.getUntrackedParameter<bool>("isRealData");
  fIsPat      = iConfig.getUntrackedParameter<bool>("isPat");

  // InputTags
  fMuonTag        = iConfig.getUntrackedParameter<edm::InputTag>("tag_muons");
  fElectronTag    = iConfig.getUntrackedParameter<edm::InputTag>("tag_electrons");
  fEleIsoTkTag    = iConfig.getUntrackedParameter<edm::InputTag>("tag_elisotk");
  fEleIsoECTag    = iConfig.getUntrackedParameter<edm::InputTag>("tag_elisoec");
  fEleIsoHCTag    = iConfig.getUntrackedParameter<edm::InputTag>("tag_elisohc");
  fEleIsoDepTkTag = iConfig.getUntrackedParameter<edm::InputTag>("tag_elisodeptk");
  fEleIsoDepECTag = iConfig.getUntrackedParameter<edm::InputTag>("tag_elisodepec");
  fEleIsoDepHCTag = iConfig.getUntrackedParameter<edm::InputTag>("tag_elisodephc");
  fMuIsoDepTkTag  = iConfig.getUntrackedParameter<edm::InputTag>("tag_muisodeptk");
  fMuIsoDepECTag  = iConfig.getUntrackedParameter<edm::InputTag>("tag_muisodepec");
  fMuIsoDepHCTag  = iConfig.getUntrackedParameter<edm::InputTag>("tag_muisodephc");
  fSCTag          = iConfig.getUntrackedParameter<edm::InputTag>("tag_sc");
  fJetTag         = iConfig.getUntrackedParameter<edm::InputTag>("tag_jets");
  fBtagTag        = iConfig.getUntrackedParameter<edm::InputTag>("tag_btag");
  fMET1Tag        = iConfig.getUntrackedParameter<edm::InputTag>("tag_met1");
  fMET2Tag        = iConfig.getUntrackedParameter<edm::InputTag>("tag_met2");
  fMET3Tag        = iConfig.getUntrackedParameter<edm::InputTag>("tag_met3");
  fMET4Tag        = iConfig.getUntrackedParameter<edm::InputTag>("tag_met4");
  fMET5Tag        = iConfig.getUntrackedParameter<edm::InputTag>("tag_met5");
  fVertexTag      = iConfig.getUntrackedParameter<edm::InputTag>("tag_vertex");
  fTrackTag       = iConfig.getUntrackedParameter<edm::InputTag>("tag_tracks");
  fCalTowTag      = iConfig.getUntrackedParameter<edm::InputTag>("tag_caltow");
  fGenPartTag     = iConfig.getUntrackedParameter<edm::InputTag>("tag_genpart");
  fTriggerTag     = iConfig.getUntrackedParameter<edm::InputTag>("tag_triggers");

  // Jet ID helper
  jetIDHelper = reco::helper::JetIDHelper(iConfig.getParameter<edm::ParameterSet>("jetID")  );

  // Event Selection
  fMinmupt        = iConfig.getParameter<double>("sel_minmupt");
  fMaxmueta       = iConfig.getParameter<double>("sel_maxmueta");
  fMinelpt        = iConfig.getParameter<double>("sel_minelpt");
  fMaxeleta       = iConfig.getParameter<double>("sel_maxeleta");
  fMaxeliso       = iConfig.getParameter<double>("sel_maxeliso");
  fMaxeld0        = iConfig.getParameter<double>("sel_maxeld0");
  fMinjpt         = iConfig.getParameter<double>("sel_minjpt");
  fMaxjeta        = iConfig.getParameter<double>("sel_maxjeta");
  fMinjemfrac     = iConfig.getParameter<double>("sel_minjemfrac");
  fMintrkpt       = iConfig.getParameter<double>("sel_mintrkpt");
  fMaxtrketa      = iConfig.getParameter<double>("sel_maxtrketa");
  fMaxtrknchi2    = iConfig.getParameter<double>("sel_maxtrknchi2");
  fMintrknhits    = iConfig.getParameter<int>("sel_mintrknhits");

  // Isolation Parameters
  fIso_MuTkDRin   = iConfig.getParameter<double>("iso_MuTkDRin");
  fIso_MuTkDRout  = iConfig.getParameter<double>("iso_MuTkDRout");
  fIso_MuTkSeed   = iConfig.getParameter<double>("iso_MuTkSeed");
  fIso_MuCalDRin  = iConfig.getParameter<double>("iso_MuCalDRin");
  fIso_MuCalDRout = iConfig.getParameter<double>("iso_MuCalDRout");
  fIso_MuCalSeed  = iConfig.getParameter<double>("iso_MuCalSeed");

  edm::LogVerbatim("NTP") << "---------------------------------";
  edm::LogVerbatim("NTP") << " ==> NTupleProducer Constructor ...";
  edm::LogVerbatim("NTP") << endl;
  edm::LogVerbatim("NTP") << "  Processing Real Data: " << (fIsRealData?"ON":"OFF"); 
  edm::LogVerbatim("NTP") << "  Processing PAT:       " << (fIsPat?"ON":"OFF"); 
  edm::LogVerbatim("NTP") << endl;
  edm::LogVerbatim("NTP") << "  Input Tags:" ;
  edm::LogVerbatim("NTP") << "    fMuonTag        = " << fMuonTag.label()        ;
  edm::LogVerbatim("NTP") << "    fElectronTag    = " << fElectronTag.label()    ;
  edm::LogVerbatim("NTP") << "    fEleIsoTkTag    = " << fEleIsoTkTag.label()    ;
  edm::LogVerbatim("NTP") << "    fEleIsoECTag    = " << fEleIsoECTag.label()    ;
  edm::LogVerbatim("NTP") << "    fEleIsoHCTag    = " << fEleIsoHCTag.label()    ;
  edm::LogVerbatim("NTP") << "    fEleIsoDepTkTag = " << fEleIsoDepTkTag.label() ;
  edm::LogVerbatim("NTP") << "    fEleIsoDepECTag = " << fEleIsoDepECTag.label() ;
  edm::LogVerbatim("NTP") << "    fEleIsoDepHCTag = " << fEleIsoDepHCTag.label() ;
  edm::LogVerbatim("NTP") << "    fMuIsoDepTkTag  = " << fMuIsoDepTkTag.label()  ;
  edm::LogVerbatim("NTP") << "    fMuIsoDepECTag  = " << fMuIsoDepECTag.label()  ;
  edm::LogVerbatim("NTP") << "    fMuIsoDepHCTag  = " << fMuIsoDepHCTag.label()  ;
  edm::LogVerbatim("NTP") << "    fSCTag          = " << fSCTag.label()          ;
  edm::LogVerbatim("NTP") << "    fJetTag         = " << fJetTag.label()         ;
  edm::LogVerbatim("NTP") << "    fMET1Tag        = " << fMET1Tag.label()        ;
  edm::LogVerbatim("NTP") << "    fMET2Tag        = " << fMET2Tag.label()        ;
  edm::LogVerbatim("NTP") << "    fMET3Tag        = " << fMET3Tag.label()        ;
  edm::LogVerbatim("NTP") << "    fMET4Tag        = " << fMET4Tag.label()        ;
  edm::LogVerbatim("NTP") << "    fMET5Tag        = " << fMET5Tag.label()        ;
  edm::LogVerbatim("NTP") << "    fVertexTag      = " << fVertexTag.label()      ;
  edm::LogVerbatim("NTP") << "    fTrackTag       = " << fTrackTag.label()       ;
  edm::LogVerbatim("NTP") << "    fCalTowTag      = " << fCalTowTag.label()      ;
  edm::LogVerbatim("NTP") << "    fGenPartTag     = " << fGenPartTag.label()     ;
  edm::LogVerbatim("NTP") << endl;
  edm::LogVerbatim("NTP") << "  Event Selection Parameters:" ;
  edm::LogVerbatim("NTP") << "    fMinmupt        = " << fMinmupt    ;
  edm::LogVerbatim("NTP") << "    fMaxmueta       = " << fMaxmueta   ;
  edm::LogVerbatim("NTP") << "    fMinelpt        = " << fMinelpt    ;
  edm::LogVerbatim("NTP") << "    fMaxeleta       = " << fMaxeleta   ;
  edm::LogVerbatim("NTP") << "    fMaxeliso       = " << fMaxeliso   ;
  edm::LogVerbatim("NTP") << "    fMaxeld0        = " << fMaxeld0    ;
  edm::LogVerbatim("NTP") << "    fMinjpt         = " << fMinjpt     ;
  edm::LogVerbatim("NTP") << "    fMaxjeta        = " << fMaxjeta    ;
  edm::LogVerbatim("NTP") << "    fMinjemfrac     = " << fMinjemfrac ;
  edm::LogVerbatim("NTP") << "    fMintrkpt       = " << fMintrkpt   ;
  edm::LogVerbatim("NTP") << "    fMaxtrketa      = " << fMaxtrketa  ;
  edm::LogVerbatim("NTP") << "    fMaxtrknchi2    = " << fMaxtrknchi2 ;
  edm::LogVerbatim("NTP") << "    fMintrknhits    = " << fMintrknhits ;
  edm::LogVerbatim("NTP") << endl;
  edm::LogVerbatim("NTP") << "  Isolation Parameters:" ;
  edm::LogVerbatim("NTP") << "    fIso_MuTkDRin   = " << fIso_MuTkDRin   ;
  edm::LogVerbatim("NTP") << "    fIso_MuTkDRout  = " << fIso_MuTkDRout  ;
  edm::LogVerbatim("NTP") << "    fIso_MuTkSeed   = " << fIso_MuTkSeed   ;
  edm::LogVerbatim("NTP") << "    fIso_MuCalDRin  = " << fIso_MuCalDRin  ;
  edm::LogVerbatim("NTP") << "    fIso_MuCalDRout = " << fIso_MuCalDRout ;
  edm::LogVerbatim("NTP") << "    fIso_MuCalSeed  = " << fIso_MuCalSeed  ;
  edm::LogVerbatim("NTP") << endl;
  edm::LogVerbatim("NTP") << "---------------------------------" ;
}

NTupleProducer::~NTupleProducer() {}

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

  ////////////////////////////////////////////////////////////////////////////////
  // Get the collections /////////////////////////////////////////////////////////
  Handle<View<Muon> > muons;
  iEvent.getByLabel(fMuonTag,muons); // 'muons'

  Handle<View<GsfElectron> > electrons;
  iEvent.getByLabel(fElectronTag, electrons); // 'gsfElectrons'

  // Jets and Jet Correctors
  Handle<View<Jet> > jets;
  iEvent.getByLabel(fJetTag,jets); // 'sisCone5CaloJets'
  const JetCorrector* L2JetCorrector = JetCorrector::getJetCorrector ("L2RelativeJetCorrectorSC5Calo",iSetup);
  const JetCorrector* L3JetCorrector = JetCorrector::getJetCorrector ("L3AbsoluteJetCorrectorSC5Calo",iSetup);

  // collect information for b-tagging
  Handle<JetTagCollection> jetsAndProbs;
  iEvent.getByLabel(fBtagTag,jetsAndProbs);

  //Get Tracks collection
  Handle<TrackCollection> tracks;
  iEvent.getByLabel(fTrackTag, tracks);

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

  // Get Transient Track Builder
  ESHandle<TransientTrackBuilder> theB;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB);

  // Get GenEventInfoProduct
  edm::Handle<GenEventInfoProduct> genEvtInfo;
  if(!fIsRealData){
    iEvent.getByLabel("generator", genEvtInfo);
    fTsigprocid = genEvtInfo->signalProcessID();		
  }

  // Dump HLT trigger bits
  Handle<TriggerResults> triggers;
  iEvent.getByLabel(fTriggerTag, triggers);
  const TriggerResults& tr = *triggers;
  if(fFirstevent){
    fFirstevent = false;
    vector<string> triggernames;
    Service<service::TriggerNamesService> tns;
    tns->getTrigPaths(*triggers, triggernames);
    for( unsigned int i = 0; i < tr.size(); i++ ){
      //cout << "i " << i << " of " << tr.size() << " : " << tr[i].accept() << " : " << triggernames[i] << endl;
      fHtrigstat->GetXaxis()->SetBinLabel(i+1, TString(triggernames[i]));
    }
    // Add a warning about the shift between trigger bits and bin numbers:
    fHtrigstat->GetXaxis()->SetBinLabel(tr.size()+2, "Bin#=Bit#+1");
    triggernames.clear();
  }
  for( unsigned int i = 0; i < tr.size(); i++ ){
    if(tr[i].accept())	fHtrigstat->Fill(i);
    fTtrigres[i] = tr[i].accept() ? 1:0;
  }

  ////////////////////////////////////////////////////////////////////////////////
  // Dump tree variables /////////////////////////////////////////////////////////
  fTrunnumber = iEvent.id().run();
  fTeventnumber = iEvent.id().event();
  fTlumisection = iEvent.luminosityBlock();

  fTweight = 1.0; // To be filled at some point?

  // Save position of primary vertex
  fTprimvtxx    = primVtx->x();
  fTprimvtxy    = primVtx->y();
  fTprimvtxz    = primVtx->z();
  fTprimvtxxE   = primVtx->xError();
  fTprimvtxyE   = primVtx->yError();
  fTprimvtxzE   = primVtx->zError();
  fTpvtxznchi2  = primVtx->normalizedChi2();
  fTpvtxntracks = primVtx->tracksSize();

  // Save position of beamspot
  fTbeamspotx = (beamSpot.position()).x();
  fTbeamspoty = (beamSpot.position()).y();
  fTbeamspotz = (beamSpot.position()).z();

  ////////////////////////////////////////////////////////
  // Muon Variables:
  int mi(-1), mqi(-1); // index of all muons and qualified muons respectively
  for(View<Muon>::const_iterator Mit = muons->begin(); Mit != muons->end(); ++Mit){
    // Check if maximum number of electrons is exceeded already:
    if(mqi >= fMaxnmus){
      edm::LogWarning("NTP") << "@SUB=analyze()"
                             << "Maximum number of muons exceeded";
      fTflagmaxmuexc = 1;
      fTgoodevent = 0;
      break;
    }
    mi++;

    // Muon preselection:
    if(!(Mit->isGlobalMuon())) continue;
    if(Mit->globalTrack()->pt() < fMinmupt) continue;
    if(fabs(Mit->globalTrack()->eta()) > fMaxmueta) continue;
    mqi++;

    // Dump muon properties in tree variables
    fTmupx[mqi]     = Mit->globalTrack()->px();
    fTmupy[mqi]     = Mit->globalTrack()->py();
    fTmupz[mqi]     = Mit->globalTrack()->pz();
    fTmupt[mqi]     = Mit->globalTrack()->pt();
    fTmueta[mqi]    = Mit->globalTrack()->eta();
    fTmuphi[mqi]    = Mit->globalTrack()->phi();
    fTmue[mqi]      = Mit->energy();
    fTmuet[mqi]     = Mit->et();
    fTmucharge[mqi] = Mit->charge();


    vector<double> MuIso = calcMuIso(&(*Mit), iEvent);
    fTmuiso[mqi]   = MuIso[0];
    fTmuptsum[mqi] = MuIso[1];
    fTmuetsum[mqi] = MuIso[2];
    MuIso.clear();

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
    fTmud0E[mqi]  = Mit->globalTrack()->dxyError();
    fTmudzbs[mqi] = Mit->globalTrack()->dz(beamSpot.position());
    fTmudzpv[mqi] = Mit->globalTrack()->dz(primVtx->position());
    fTmudzE[mqi]  = Mit->globalTrack()->dzError();

    fTmunchi2[mqi]     = Mit->globalTrack()->normalizedChi2();
    fTmunglhits[mqi]   = Mit->globalTrack()->hitPattern().numberOfValidHits();
    fTmunmuhits[mqi]   = Mit->outerTrack()->hitPattern().numberOfValidHits(); // always 0
    fTmuntkhits[mqi]   = Mit->innerTrack()->hitPattern().numberOfValidHits();
    fTmunmatches[mqi]  = Mit->numberOfMatches();
    fTmunchambers[mqi] = Mit->numberOfChambers();

    fTmucalocomp[mqi] = Mit->caloCompatibility();
    fTmusegmcomp[mqi] = muon::segmentCompatibility(*Mit);
    fTmutrackermu[mqi] = Mit->isTrackerMuon() ? 1:0;
    fTmuisGMPT[mqi] = muon::isGoodMuon(*Mit, muon::GlobalMuonPromptTight) ? 1:0;

    // Matching
    if(!fIsRealData){
      vector<int> MuMatch = matchMuCand(&(*Mit), iEvent);
      if( MuMatch[0] ){
        fTmuid[mqi]  = MuMatch[1];
        fTmumid[mqi] = MuMatch[2];
      }
      MuMatch.clear();
    }
    fTgoodmu[mqi] = 1;
  }
  fTnmu = mqi+1;

  ////////////////////////////////////////////////////////
  // Electron variables:
  int eqi(-1),ei(-1); // counts # of qualified electrons
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
      if(eqi >= fMaxneles){
        edm::LogWarning("NTP") << "@SUB=analyze"
                               << "Maximum number of electrons exceeded..";
        fTflagmaxelexc = 1;
        fTgoodevent = 0;
        break;
      }

      // Electron preselection:   
      if(El->pt() < fMinelpt) continue;
      if(fabs(El->eta()) > fMaxeleta) continue;
      if(fabs(El->gsfTrack()->dxy(beamSpot.position())) > fMaxeld0) continue;

      // Read EleIsolation
      double eTkIso,eECIso,eHCIso;
      if ( !fIsPat ) { // Isolation is embedded in PAT => switch
        Handle<edm::ValueMap<double> > eIsoTkValueMap;
        iEvent.getByLabel(fEleIsoTkTag, eIsoTkValueMap);
        const ValueMap<double> &eTkIsoMap = *eIsoTkValueMap.product();
        Handle<edm::ValueMap<double> > eIsoECValueMap;
        iEvent.getByLabel(fEleIsoECTag, eIsoECValueMap);
        const ValueMap<double> &eECIsoMap = *eIsoECValueMap.product();
        Handle<edm::ValueMap<double> > eIsoHCValueMap;
        iEvent.getByLabel(fEleIsoHCTag, eIsoHCValueMap);
        const ValueMap<double> &eHCIsoMap = *eIsoHCValueMap.product();

        Ref<View<GsfElectron> > electronRef(electrons,ei);
        eTkIso = eTkIsoMap[electronRef];
        eECIso = eECIsoMap[electronRef];
        eHCIso = eHCIsoMap[electronRef];
      } else {
        const pat::Electron* pEl = static_cast<const pat::Electron*>(&(*El));
        eTkIso = pEl->trackIso();
        eECIso = pEl->ecalIso();
        eHCIso = pEl->hcalIso();
      }
      double eliso  = (eTkIso+eECIso+eHCIso)/El->pt();
      if(eliso > fMaxeliso) continue;

      // Dump electron properties in tree variables
      eqi++;

      fTepx[eqi]     = El->gsfTrack()->px();
      fTepy[eqi]     = El->gsfTrack()->py();
      fTepz[eqi]     = El->gsfTrack()->pz();
      fTept[eqi]     = El->pt();
      fTeeta[eqi]    = El->eta();
      fTephi[eqi]    = El->phi();
      fTee[eqi]      = El->energy();
      fTeet[eqi]     = El->et(); // i know: it's the same as pt, still...
      fTed0bs[eqi]   = -1.0*El->gsfTrack()->dxy(beamSpot.position());
      fTed0pv[eqi]   = -1.0*El->gsfTrack()->dxy(primVtx->position());
      fTed0E[eqi]    = El->gsfTrack()->dxyError();
      fTedzbs[eqi]   = El->gsfTrack()->dz(beamSpot.position());
      fTedzpv[eqi]   = El->gsfTrack()->dz(primVtx->position());
      fTedzE[eqi]    = El->gsfTrack()->dzError();
      fTenchi2[eqi]  = El->gsfTrack()->normalizedChi2();
      fTeiso[eqi]    = eliso;
      fTeptsum[eqi]  = eTkIso;
      fTeetsum[eqi]  = eECIso+eHCIso;
      fTecharge[eqi] = El->charge();
      fTeInGap[eqi]         = El->isGap() ? 1:0;  // no bug... //*** Beni
      fTeEcalDriven[eqi]    = El->isEcalDriven() ? 1:0;
      fTeTrackerDriven[eqi] = El->isTrackerDriven() ? 1:0;
      fTeBasicClustersSize[eqi] = El->basicClustersSize();
      fTefbrem[eqi]  = El->fbrem();
      fTeHcalOverEcal[eqi] = El->hcalOverEcal();                               
      fTeE5x5[eqi]   = El->e5x5();                                           
      fTeE2x5Max[eqi] = El->e2x5Max();                                           
      fTeSigmaIetaIeta[eqi] = El->sigmaIetaIeta();                         
      fTeDeltaEtaSeedClusterAtCalo[eqi] = El->deltaEtaSeedClusterTrackAtCalo(); 
      fTeDeltaPhiSeedClusterAtCalo[eqi] = El->deltaPhiSeedClusterTrackAtCalo(); 
      fTeDeltaPhiSuperClusterAtVtx[eqi] = El->deltaPhiSuperClusterTrackAtVtx(); 
      fTeDeltaEtaSuperClusterAtVtx[eqi] = El->deltaEtaSuperClusterTrackAtVtx(); 
      fTecaloenergy[eqi] = El->caloEnergy();
      fTtrkmomatvtx[eqi] = El->trackMomentumAtVtx().R();
      fTeESuperClusterOverP[eqi] = El->eSuperClusterOverP();
      // Read in Electron ID
      if ( !fIsPat ) { // Electron ID is embedded in PAT => switch
        Ref<View<GsfElectron> > electronRef(electrons,ei);
        fTeIDTight[eqi]       = eIDmapT[electronRef]  ? 1:0;
        fTeIDLoose[eqi]       = eIDmapL[electronRef]  ? 1:0;
        fTeIDRobustTight[eqi] = eIDmapRT[electronRef] ? 1:0;
        fTeIDRobustLoose[eqi] = eIDmapRL[electronRef] ? 1:0;
      } else {
        const pat::Electron* pEl = static_cast<const pat::Electron*>(&(*El));
        fTeIDTight[eqi]       = pEl->electronID("eidTight")>0. ? 1:0;
        fTeIDLoose[eqi]       = pEl->electronID("eidLoose")>0. ? 1:0;
        fTeIDRobustTight[eqi] = pEl->electronID("eidRobustTight")>0. ? 1:0;
        fTeIDRobustLoose[eqi] = pEl->electronID("eidRobustLoose")>0. ? 1:0;
      }

      fTgoodel[eqi] = 1;
    }
  }
  fTneles = eqi+1;

  ////////////////////////////////////////////////////////
  // Jet Variables:
  // Check Jet collection's size
  if( (int)jets->size() > fMaxnjets ){
    edm::LogWarning("NTP") << "@SUB=analyze"
                           << "Found more than " 
                           << fMaxnjets << " uncorrected jets, I'm scared ...";
    fTflagmaxujetexc = 1;
    fTgoodevent = 0;
    fEventTree->Fill();
    return;
  }

  // Jet corrections
  vector<const Jet*> corrJets;
  corrJets.reserve(jets->size()); // Speed up push-backs
  vector<OrderPair> corrIndices;  // Vector of indices and pt of corrected jets to re-order them
  // I save here px, py, pz for uncorrected jets
  // used later for matching jets with b-tagging information
  int itag(-1);
  for(View<Jet>::const_iterator Jit = jets->begin(); Jit != jets->end(); ++Jit){// Loop over uncorr. jets
    itag++;
    if ( !fIsPat ) {
      // Cast and make a new (non-const) copy to apply corrections
      const CaloJet* jet = static_cast<const CaloJet*>(&(*Jit));
      CaloJet* j = new CaloJet(*jet);
      fJUNC_px_match[itag] = j->px();
      fJUNC_py_match[itag] = j->py();
      fJUNC_pz_match[itag] = j->pz();
      // Apply L2 and L3 JetCorrections
      j->scaleEnergy(L2JetCorrector->correction(j->p4()));
      j->scaleEnergy(L3JetCorrector->correction(j->p4()));

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

  // Determine qualified jets
  int jqi(-1); // counts # of qualified jets
  for(vector<OrderPair>::const_iterator it = corrIndices.begin();
      it != corrIndices.end(); ++it ) { // Loop over corr. jet indices
    // Check if maximum number of jets is exceeded already
    if(jqi >= fMaxnjets) {
      edm::LogWarning("NTP") << "@SUB=analyze"
                             << "Maximum number of jets exceeded";
      fTflagmaxjetexc = 1;
      fTgoodevent = 0;
      break;
    }
    int index = it->first;

    // Jet preselection
    const Jet* jet = corrJets[index];
    if(jet->pt() < fMinjpt) continue;
    if(fabs(jet->eta()) > fMaxjeta) continue;
    // if(jet->emEnergyFraction() < fMinjemfrac) continue;
    jqi++;

    // Dump jet properties into tree variables
    fTjpx[jqi]  = jet->px();
    fTjpy[jqi]  = jet->py();
    fTjpz[jqi]  = jet->pz();
    fTjpt[jqi]  = jet->pt();
    fTjeta[jqi] = jet->eta();
    fTjphi[jqi] = jet->phi();
    fTje[jqi]   = jet->energy();
    fTjet[jqi]  = jet->et(); // i know: it's the same as pt, still...
    fTjEcorr[jqi] =  jet->px()/fJUNC_px_match[index];

    // Calculate the DR wrt the closest electron
    double ejDRmin=10.; // Default when no electrons previously selected
    for(int j = 0; j < fTneles; j++) {
      double ejDR = reco::deltaR(jet->eta(), jet->phi(), fTeeta[j], fTephi[j]);
      if(ejDR<ejDRmin) ejDRmin = ejDR;
    }
    fTjeMinDR[jqi]=ejDRmin;

    // B-tagging probability
    if ( !fIsPat ) {
      for (unsigned int i = 0; i < jetsAndProbs->size(); i++){
        // Match by pt between the two "collections"
        if (fabs( (fJUNC_px_match[index] - (*jetsAndProbs)[i].first->px())/fJUNC_px_match[index]) < 0.00001 && 
            fabs( (fJUNC_py_match[index] - (*jetsAndProbs)[i].first->py())/fJUNC_py_match[index]) < 0.00001 &&
            fabs( (fJUNC_pz_match[index] - (*jetsAndProbs)[i].first->pz())/fJUNC_pz_match[index]) < 0.00001)  {
          fTbTagProb[jqi]=(*jetsAndProbs)[i].second;
          break;
        }
      }
    } else {
      const pat::Jet* pJ = static_cast<const pat::Jet*>(jet);      
      fTbTagProb[jqi] = pJ->bDiscriminator(fBtagTag.label());
    }

    // Jet-track association: get associated tracks
    vector<const reco::Track*> AssociatedTracks;
    fTnAssoTracks[jqi]=0;
    fTChfrac[jqi]=-1.; // Default (if jet-tracks association cone is outside tracker acceptance)
    if ( !fIsPat ) {
      AssociatedTracks = FindAssociatedTracks(jet, tracks.product());
    } else {
      const pat::Jet* pJ = static_cast<const pat::Jet*>(jet);
      for ( TrackRefVector::iterator it = pJ->associatedTracks().begin();
            it != pJ->associatedTracks().end(); ++it ) {
        AssociatedTracks.push_back( it->get() );
      }
      //DO NOT use this to keep in sync. with non-PAT code
      //fTChfrac[jqi] = pJ.jetCharge();
      //fTnAssoTracks[jqi] = pJ.associatedTracks().size();
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

      // Cast and make a new (non-const) copy of corrected calojet, to loop over the calotowers
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
    fTnAssoTracks[jqi]=0;
    if(fabs(jet->eta())<2.9){ // when the cone of dR=0.5 around the jet is (at least partially) inside the tracker acceptance
      // tmp variables for vectorial sum of pt of tracks
      double pXtmp=0.;
      double pYtmp=0.;

      for(size_t t = 0; t < AssociatedTracks.size(); ++t){
        AssociatedTTracks.push_back(theB->build(AssociatedTracks[t]));
        if(AssociatedTracks[t]->normalizedChi2()<10. && AssociatedTracks[t]->numberOfValidHits()>10 && AssociatedTracks[t]->pt()>1.){
          pXtmp+=AssociatedTracks[t]->px();
          pYtmp+=AssociatedTracks[t]->py();
          fTnAssoTracks[jqi]+=1;
        }
	//loop and find the three higher pT tracks
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
      //fill the momenta
      if(AssociatedTracks.size()>=1){
        fTtrk1px[jqi] =AssociatedTracks[idx1]->px();
        fTtrk1py[jqi] =AssociatedTracks[idx1]->py();
        fTtrk1pz[jqi] =AssociatedTracks[idx1]->pz();
      }
      if(AssociatedTracks.size()>=2){
        fTtrk2px[jqi] =AssociatedTracks[idx2]->px();
        fTtrk2py[jqi] =AssociatedTracks[idx2]->py();
        fTtrk2pz[jqi] =AssociatedTracks[idx2]->pz();
      }
      if(AssociatedTracks.size()>=3){
        fTtrk3px[jqi] =AssociatedTracks[idx3]->px();
        fTtrk3py[jqi] =AssociatedTracks[idx3]->py();
        fTtrk3pz[jqi] =AssociatedTracks[idx3]->pz();
      }
      //
      fTChfrac[jqi]=sqrt(pow(pXtmp,2)+pow(pYtmp,2))/jet->pt();
    } else {//the whole cone used for jet-tracks association is outside of the tracker acceptance
      fTChfrac[jqi]=-1.;
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

    // CaloJet specific variables (embedded in PAT)
    if ( !fIsPat ) {
      const CaloJet* cJ = static_cast<const CaloJet*>(jet);
      fTjemfrac[jqi] =  cJ->emEnergyFraction();

      // jetID variables
      // FIXME: should be calculated before correction (for fractions)
      // => re-correct all fractions
      jetIDHelper.calculate(iEvent, *cJ);
      fTjID_HPD[jqi]      = jetIDHelper.fHPD()*fTjEcorr[jqi];
      fTjID_RBX[jqi]      = jetIDHelper.fRBX()*fTjEcorr[jqi];
      fTjID_n90Hits[jqi]  = jetIDHelper.n90Hits();
      fTjID_SubDet1[jqi]  = jetIDHelper.fSubDetector1()*fTjEcorr[jqi];
      fTjID_SubDet2[jqi]  = jetIDHelper.fSubDetector2()*fTjEcorr[jqi];
      fTjID_SubDet3[jqi]  = jetIDHelper.fSubDetector3()*fTjEcorr[jqi];
      fTjID_SubDet4[jqi]  = jetIDHelper.fSubDetector4()*fTjEcorr[jqi];
      fTjID_resEMF[jqi]   = jetIDHelper.restrictedEMF();
      fTjID_HCALTow[jqi]  = jetIDHelper.nHCALTowers();
      fTjID_ECALTow[jqi]  = jetIDHelper.nECALTowers();
    } else {
      const pat::Jet* pJ = static_cast<const pat::Jet*>(jet);
      fTjemfrac[jqi]      = pJ->emEnergyFraction();
      fTjID_HPD[jqi]      = pJ->jetID().fHPD;
      fTjID_RBX[jqi]      = pJ->jetID().fRBX;
      fTjID_n90Hits[jqi]  = pJ->jetID().n90Hits;
      fTjID_SubDet1[jqi]  = pJ->jetID().fSubDetector1;
      fTjID_SubDet2[jqi]  = pJ->jetID().fSubDetector2;
      fTjID_SubDet3[jqi]  = pJ->jetID().fSubDetector3;
      fTjID_SubDet4[jqi]  = pJ->jetID().fSubDetector4;
      fTjID_resEMF[jqi]   = pJ->jetID().restrictedEMF;
      fTjID_HCALTow[jqi]  = pJ->jetID().nHCALTowers;
      fTjID_ECALTow[jqi]  = pJ->jetID().nECALTowers;
    }
    AssociatedTracks.clear();
    AssociatedTTracks.clear();

    fTgoodjet[jqi] = 1;
  }
  fTnjets = jqi+1;
  corrJets.clear();
  corrIndices.clear();

  ////////////////////////////////////////////////////////
  // Get and Dump (M)E(T) Variables:
  // Tracks:
  int nqtrk(-1);
  fTTrkPtSumx = 0.; fTTrkPtSumy = 0.;
  for(TrackCollection::const_iterator it = tracks->begin(); 
      it != tracks->end() ; ++it ) {
    if(it->pt() < fMintrkpt) continue;
    if(fabs(it->eta()) > fMaxtrketa) continue;
    if(it->normalizedChi2() > fMaxtrknchi2) continue;
    if(it->numberOfValidHits() < fMintrknhits) continue;
    nqtrk++; // starts at 0
    // Check if maximum number of tracks is exceeded already
    if(nqtrk >= fMaxntrks) {
      edm::LogWarning("NTP") << "@SUB=analyze"
                             << "Maximum number of tracks exceeded";
      fTflagmaxtrkexc = 1;
      fTgoodevent = 0;
      break;
    }

    fTTrkPtSumx += it->px();
    fTTrkPtSumy += it->py();
    fTtrkpt[nqtrk]    = it->pt();
    fTtrketa[nqtrk]   = it->eta();
    fTtrkphi[nqtrk]   = it->phi();
    fTtrknchi2[nqtrk] = it->normalizedChi2();
    fTtrknhits[nqtrk] = it->numberOfValidHits();
    fTgoodtrk[nqtrk]  = 1;
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
  fTRawMET    = (calomet->at(0)).pt();
  fTRawMETpx  = (calomet->at(0)).px();
  fTRawMETpy  = (calomet->at(0)).py();
  fTRawMETphi = (calomet->at(0)).phi();

  fTMuCorrMET    = (corrmumet->at(0)).pt();
  fTMuCorrMETpx  = (corrmumet->at(0)).px();
  fTMuCorrMETpy  = (corrmumet->at(0)).py();
  fTMuCorrMETphi = (corrmumet->at(0)).phi();

  fTTCMET    = (tcmet->at(0)).pt();
  fTTCMETpx  = (tcmet->at(0)).px();
  fTTCMETpy  = (tcmet->at(0)).py();
  fTTCMETphi = (tcmet->at(0)).phi();

  fTPFMET    = (pfmet->front()).pt();
  fTPFMETpx  = (pfmet->front()).px();
  fTPFMETpy  = (pfmet->front()).py();
  fTPFMETphi = (pfmet->front()).phi();

  fTMuJESCorrMET    = (corrmujesmet->at(0)).pt();
  fTMuJESCorrMETpx  = (corrmujesmet->at(0)).px();
  fTMuJESCorrMETpy  = (corrmujesmet->at(0)).py();
  fTMuJESCorrMETphi = (corrmujesmet->at(0)).phi();

  ////////////////////////////////////////////////////////////////////////////////
  // Fill Tree ///////////////////////////////////////////////////////////////////
  fEventTree->Fill();
  fNFillTree++;
}

// Method called once each job just before starting event loop
void NTupleProducer::beginJob(const edm::EventSetup&){
  fNTotEvents = 0;
  fNFillTree  = 0;
  fFirstevent = true;
  fHtrigstat = fTFileService->make<TH1I>("TriggerStats", "TriggerStatistics", 200, 0, 200);

  // Tree with run information
  fRunTree = fTFileService->make<TTree>("RunInfo", "ETHZRunAnalysisTree");
  fRunTree->Branch("Run"            ,&fRTrunnumber,      "Run/I");
  fRunTree->Branch("ExtXSecLO"      ,&fRTextxslo,        "ExtXSecLO/D");
  fRunTree->Branch("ExtXSecNLO"     ,&fRTextxsnlo,       "ExtXSecNLO/D");
  fRunTree->Branch("IntXSec"        ,&fRTintxs,          "IntXSec/D");
  fRunTree->Branch("MinMuPt"        ,&fRTMinmupt,        "MinMuPt/D");
  fRunTree->Branch("MaxMuEta"       ,&fRTMaxmueta,       "MaxMuEta/D");
  fRunTree->Branch("MinElPt"        ,&fRTMinelpt,        "MinElPt/D");
  fRunTree->Branch("MaxElEta"       ,&fRTMaxeleta,       "MaxElEta/D");
  fRunTree->Branch("MaxElIso"       ,&fRTMaxeliso,       "MaxElIso/D");
  fRunTree->Branch("MaxElD0"        ,&fRTMaxeld0,        "MaxElD0/D");
  fRunTree->Branch("MinJPt"         ,&fRTMinjpt,         "MinJPt/D");
  fRunTree->Branch("MaxJEta"        ,&fRTMaxjeta,        "MaxJEta/D");
  fRunTree->Branch("MinJEMfrac"     ,&fRTMinjemfrac,     "MinJEMfrac/D");

  fRunTree->Branch("MinTrkPt"       ,&fRTMintrkpt,        "MinTrkPt/D");
  fRunTree->Branch("MaxTrkEta"      ,&fRTMaxtrketa,       "MaxTrkEta/D");
  fRunTree->Branch("MaxTrkNChi2"    ,&fRTMaxtrknchi2,     "MaxTrkNChi2/D");
  fRunTree->Branch("MinTrkNHits"    ,&fRTMintrknhits,     "MinTrkNHits/I");

  fRunTree->Branch("IsoMuTkDRin"    ,&fRTIsoMuTkDRin,    "IsoMuTkDRin/D");
  fRunTree->Branch("IsoMuTkDRout"   ,&fRTIsoMuTkDRout,   "IsoMuTkDRout/D");
  fRunTree->Branch("IsoMuTkSeed"    ,&fRTIsoMuTkSeed,    "IsoMuTkSeed/D");
  fRunTree->Branch("IsoMuCalDRin"   ,&fRTIsoMuCalDRin,   "IsoMuCalDRin/D");
  fRunTree->Branch("IsoMuCalDRout"  ,&fRTIsoMuCalDRout,  "IsoMuCalDRout/D");
  fRunTree->Branch("IsoMuCalSeed"   ,&fRTIsoMuCalSeed,   "IsoMuCalSeed/D");

  // Tree with event information
  fEventTree = fTFileService->make<TTree>("Analysis", "ETHZAnalysisTree");
  // Event information:
  fEventTree->Branch("Run"            ,&fTrunnumber      ,"Run/I");
  fEventTree->Branch("Event"          ,&fTeventnumber    ,"Event/I");
  fEventTree->Branch("LumiSection"    ,&fTlumisection    ,"LumiSection/I");
  fEventTree->Branch("SigProcID"      ,&fTsigprocid      ,"SigProcID/I");
  fEventTree->Branch("ExtXSecLO"      ,&fTextxslo        ,"ExtXSecLO/D");
  fEventTree->Branch("IntXSec"        ,&fTintxs          ,"IntXSec/D");
  fEventTree->Branch("Weight"         ,&fTweight         ,"Weight/D");
  fEventTree->Branch("TrigResults"    ,&fTtrigres        ,"TrigResults[200]/I");
  fEventTree->Branch("PrimVtxx"       ,&fTprimvtxx       ,"PrimVtxx/D");
  fEventTree->Branch("PrimVtxy"       ,&fTprimvtxy       ,"PrimVtxy/D");
  fEventTree->Branch("PrimVtxz"       ,&fTprimvtxz       ,"PrimVtxz/D");	
  fEventTree->Branch("PrimVtxxE"      ,&fTprimvtxxE      ,"PrimVtxxE/D");
  fEventTree->Branch("PrimVtxyE"      ,&fTprimvtxyE      ,"PrimVtxyE/D");
  fEventTree->Branch("PrimVtxzE"      ,&fTprimvtxzE      ,"PrimVtxzE/D");
  fEventTree->Branch("PrimVtxNChi2"   ,&fTpvtxznchi2     ,"PrimVtxNChi2/D");
  fEventTree->Branch("PrimVtxNTracks" ,&fTpvtxntracks    ,"PrimVtxNTracks/I");
  fEventTree->Branch("Beamspotx"      ,&fTbeamspotx      ,"Beamspotx/D");
  fEventTree->Branch("Beamspoty"      ,&fTbeamspoty      ,"Beamspoty/D");
  fEventTree->Branch("Beamspotz"      ,&fTbeamspotz      ,"Beamspotz/D");
  fEventTree->Branch("NCaloTowers"    ,&fTNCaloTowers    ,"NCaloTowers/I");
  fEventTree->Branch("GoodEvent"      ,&fTgoodevent      ,"GoodEvent/I");
  fEventTree->Branch("MaxMuExceed"    ,&fTflagmaxmuexc   ,"MaxMuExceed/I");
  fEventTree->Branch("MaxElExceed"    ,&fTflagmaxelexc   ,"MaxElExceed/I");
  fEventTree->Branch("MaxJetExceed"   ,&fTflagmaxjetexc  ,"MaxJetExceed/I");
  fEventTree->Branch("MaxUncJetExceed",&fTflagmaxujetexc ,"MaxUncJetExceed/I");

  // Muons:
  fEventTree->Branch("NMus"           ,&fTnmu            ,"NMus/I");
  fEventTree->Branch("MuGood"         ,&fTgoodmu         ,"MuGood[NMus]/I");
  fEventTree->Branch("MuPx"           ,&fTmupx           ,"MuPx[NMus]/D");
  fEventTree->Branch("MuPy"           ,&fTmupy           ,"MuPy[NMus]/D");
  fEventTree->Branch("MuPz"           ,&fTmupz           ,"MuPz[NMus]/D");
  fEventTree->Branch("MuPt"           ,&fTmupt           ,"MuPt[NMus]/D");
  fEventTree->Branch("MuE"            ,&fTmue            ,"MuE[NMus]/D");
  fEventTree->Branch("MuEt"           ,&fTmuet           ,"MuEt[NMus]/D");
  fEventTree->Branch("MuEta"          ,&fTmueta          ,"MuEta[NMus]/D");
  fEventTree->Branch("MuPhi"          ,&fTmuphi          ,"MuPhi[NMus]/D");
  fEventTree->Branch("MuCharge"       ,&fTmucharge       ,"MuCharge[NMus]/I");
  fEventTree->Branch("MuPtsum"        ,&fTmuptsum        ,"MuPtsum[NMus]/D");
  fEventTree->Branch("MuEtsum"        ,&fTmuetsum        ,"MuEtsum[NMus]/D");
  fEventTree->Branch("MuIso"          ,&fTmuiso          ,"MuIso[NMus]/D");
  fEventTree->Branch("MuEem"          ,&fTmueecal        ,"MuEem[NMus]/D");
  fEventTree->Branch("MuEhad"         ,&fTmuehcal        ,"MuEhad[NMus]/D");
  fEventTree->Branch("MuD0BS"         ,&fTmud0bs         ,"MuD0BS[NMus]/D");
  fEventTree->Branch("MuD0PV"         ,&fTmud0pv         ,"MuD0PV[NMus]/D");
  fEventTree->Branch("MuD0E"          ,&fTmud0E          ,"MuD0E[NMus]/D");
  fEventTree->Branch("MuDzBS"         ,&fTmudzbs         ,"MuDzBS[NMus]/D");
  fEventTree->Branch("MuDzPV"         ,&fTmudzpv         ,"MuDzPV[NMus]/D");
  fEventTree->Branch("MuDzE"          ,&fTmudzE          ,"MuDzE[NMus]/D");
  fEventTree->Branch("MuNChi2"        ,&fTmunchi2        ,"MuNChi2[NMus]/D");
  fEventTree->Branch("MuNGlHits"      ,&fTmunglhits      ,"MuNGlHits[NMus]/I");
  fEventTree->Branch("MuNMuHits"      ,&fTmunmuhits      ,"MuNMuHits[NMus]/I");
  fEventTree->Branch("MuNTkHits"      ,&fTmuntkhits      ,"MuNTkHits[NMus]/I");
  fEventTree->Branch("MuNMatches"     ,&fTmunmatches     ,"MuNMatches[NMus]/I");
  fEventTree->Branch("MuNChambers"    ,&fTmunchambers    ,"MuNChambers[NMus]/I");
  fEventTree->Branch("MuCaloComp"     ,&fTmucalocomp     ,"MuCaloComp[NMus]/D");
  fEventTree->Branch("MuSegmComp"     ,&fTmusegmcomp     ,"MuSegmComp[NMus]/D");
  fEventTree->Branch("MuTrackerMu"    ,&fTmutrackermu    ,"MuTrackerMu[NMus]/I");
  fEventTree->Branch("MuGMPT"         ,&fTmuisGMPT       ,"MuGMPT[NMus]/I");
  fEventTree->Branch("MuID"           ,&fTmuid           ,"MuID[NMus]/I");
  fEventTree->Branch("MuMID"          ,&fTmumid          ,"MuMID[NMus]/I");

  // Electrons:
  fEventTree->Branch("NEles"            ,&fTneles           ,"NEles/I");
  fEventTree->Branch("ElGood"           ,&fTgoodel          ,"ElGood[NEles]/I");
  fEventTree->Branch("ElPx"             ,&fTepx             ,"ElPx[NEles]/D");
  fEventTree->Branch("ElPy"             ,&fTepy             ,"ElPy[NEles]/D");
  fEventTree->Branch("ElPz"             ,&fTepz             ,"ElPz[NEles]/D");
  fEventTree->Branch("ElPt"             ,&fTept             ,"ElPt[NEles]/D");
  fEventTree->Branch("ElE"              ,&fTee              ,"ElE[NEles]/D");
  fEventTree->Branch("ElEt"             ,&fTeet             ,"ElEt[NEles]/D");
  fEventTree->Branch("ElEta"            ,&fTeeta            ,"ElEta[NEles]/D");
  fEventTree->Branch("ElPhi"            ,&fTephi            ,"ElPhi[NEles]/D");
  fEventTree->Branch("ElD0BS"           ,&fTed0bs           ,"ElD0BS[NEles]/D");
  fEventTree->Branch("ElD0PV"           ,&fTed0pv           ,"ElD0PV[NEles]/D");
  fEventTree->Branch("ElD0E"            ,&fTed0E            ,"ElD0E[NEles]/D");
  fEventTree->Branch("ElDzBS"           ,&fTedzbs           ,"ElDzBS[NEles]/D");
  fEventTree->Branch("ElDzPV"           ,&fTedzpv           ,"ElDzPV[NEles]/D");
  fEventTree->Branch("ElDzE"            ,&fTedzE            ,"ElDzE[NEles]/D");
  fEventTree->Branch("ElIso"            ,&fTeiso            ,"ElIso[NEles]/D");
  fEventTree->Branch("ElPtSum"          ,&fTeptsum          ,"ElPtSum[NEles]/D");
  fEventTree->Branch("ElEtSum"          ,&fTeetsum          ,"ElEtSum[NEles]/D");
  fEventTree->Branch("ElNChi2"          ,&fTenchi2          ,"ElNChi2[NEles]/D");
  fEventTree->Branch("ElCharge"         ,&fTecharge         ,"ElCharge[NEles]/I");
  fEventTree->Branch("ElIDTight"        ,&fTeIDTight        ,"ElIDTight[NEles]/I");
  fEventTree->Branch("ElIDLoose"        ,&fTeIDLoose        ,"ElIDLoose[NEles]/I");
  fEventTree->Branch("ElIDRobustTight"  ,&fTeIDRobustTight  ,"ElIDRobustTight[NEles]/I");
  fEventTree->Branch("ElIDRobustLoose"  ,&fTeIDRobustLoose  ,"ElIDRobustLoose[NEles]/I");
  fEventTree->Branch("ElInGap"             ,&fTeInGap             ,"ElInGap[NEles]/I");
  fEventTree->Branch("ElEcalDriven"        ,&fTeEcalDriven        ,"ElEcalDriven[NEles]/I");
  fEventTree->Branch("ElTrackerDriven"     ,&fTeTrackerDriven     ,"ElTrackerDriven[NEles]/I");
  fEventTree->Branch("ElBasicClustersSize" ,&fTeBasicClustersSize ,"ElBasicClustersSize[NEles]/I");
  fEventTree->Branch("Elfbrem"             ,&fTefbrem             ,"Elfbrem[NEles]/D");
  fEventTree->Branch("ElHcalOverEcal"      ,&fTeHcalOverEcal      ,"ElHcalOverEcal[NEles]/D");
  fEventTree->Branch("ElE5x5"              ,&fTeE5x5              ,"ElE5x5[NEles]/D");
  fEventTree->Branch("ElE2x5Max"           ,&fTeE2x5Max           ,"ElE2x5Max[NEles]/D");
  fEventTree->Branch("ElSigmaIetaIeta"     ,&fTeSigmaIetaIeta     ,"ElSigmaIetaIeta[NEles]/D");
  fEventTree->Branch("ElDeltaPhiSeedClusterAtCalo" ,&fTeDeltaPhiSeedClusterAtCalo ,"ElDeltaPhiSeedClusterAtCalo[NEles]/D");
  fEventTree->Branch("ElDeltaEtaSeedClusterAtCalo" ,&fTeDeltaEtaSeedClusterAtCalo ,"ElDeltaEtaSeedClusterAtCalo[NEles]/D");
  fEventTree->Branch("ElDeltaPhiSuperClusterAtVtx" ,&fTeDeltaPhiSuperClusterAtVtx ,"ElDeltaPhiSuperClusterAtVtx[NEles]/D");
  fEventTree->Branch("ElDeltaEtaSuperClusterAtVtx" ,&fTeDeltaEtaSuperClusterAtVtx ,"ElDeltaEtaSuperClusterAtVtx[NEles]/D");
  fEventTree->Branch("ElCaloEnergy"    ,&fTecaloenergy        ,"ElCaloEnergy[NEles]/D");
  fEventTree->Branch("ElTrkMomAtVtx"   ,&fTtrkmomatvtx        ,"ElTrkMomAtVtx[NEles]/D");
  fEventTree->Branch("ElESuperClusterOverP"        ,&fTeESuperClusterOverP        ,"ElESuperClusterOverP[NEles]/D");

  // Jets:
  fEventTree->Branch("NJets"          ,&fTnjets        ,"NJets/I");
  fEventTree->Branch("JGood"          ,&fTgoodjet      ,"JGood[NJets]/I");
  fEventTree->Branch("JPx"            ,&fTjpx          ,"JPx[NJets]/D");
  fEventTree->Branch("JPy"            ,&fTjpy          ,"JPy[NJets]/D");
  fEventTree->Branch("JPz"            ,&fTjpz          ,"JPz[NJets]/D");
  fEventTree->Branch("JPt"            ,&fTjpt          ,"JPt[NJets]/D");
  fEventTree->Branch("JE"             ,&fTje           ,"JE[NJets]/D");
  fEventTree->Branch("JEt"            ,&fTjet          ,"JEt[NJets]/D");
  fEventTree->Branch("JEta"           ,&fTjeta         ,"JEta[NJets]/D");
  fEventTree->Branch("JPhi"           ,&fTjphi         ,"JPhi[NJets]/D");
  fEventTree->Branch("JEMfrac"        ,&fTjemfrac      ,"JEMfrac[NJets]/D");
  fEventTree->Branch("JID_HPD"        ,&fTjID_HPD      ,"JID_HPD[NJets]/D");
  fEventTree->Branch("JID_RBX"        ,&fTjID_RBX      ,"JID_RBX[NJets]/D");
  fEventTree->Branch("JID_n90Hits"    ,&fTjID_n90Hits  ,"JID_n90Hits[NJets]/D");
  fEventTree->Branch("JID_SubDet1"    ,&fTjID_SubDet1  ,"JID_SubDet1[NJets]/D");
  fEventTree->Branch("JID_SubDet2"    ,&fTjID_SubDet2  ,"JID_SubDet2[NJets]/D");
  fEventTree->Branch("JID_SubDet3"    ,&fTjID_SubDet3  ,"JID_SubDet3[NJets]/D");
  fEventTree->Branch("JID_SubDet4"    ,&fTjID_SubDet4  ,"JID_SubDet4[NJets]/D");
  fEventTree->Branch("JID_resEMF"     ,&fTjID_resEMF   ,"JID_resEMF[NJets]/D");
  fEventTree->Branch("JID_HCALTow"    ,&fTjID_HCALTow  ,"JID_HCALTow[NJets]/D");
  fEventTree->Branch("JID_ECALTow"    ,&fTjID_ECALTow  ,"JID_ECALTow[NJets]/D");
  fEventTree->Branch("JEtaEMrms"      ,&fTJEtaEMrms    ,"JEtaEMrms[NJets]/D");
  fEventTree->Branch("JEtaHADrms"     ,&fTJEtaHADrms   ,"JEtaHADrms[NJets]/D");
  fEventTree->Branch("JPhiEMrms"      ,&fTJPhiEMrms    ,"JPhiEMrms[NJets]/D");
  fEventTree->Branch("JPhiHADrms"     ,&fTJPhiHADrms   ,"JPhiHADrms[NJets]/D");
  fEventTree->Branch("JbTagProb"      ,&fTbTagProb     ,"JbTagProb[NJets]/D");
  fEventTree->Branch("JChfrac"        ,&fTChfrac       ,"JChfrac[NJets]/D");
  fEventTree->Branch("JNAssoTracks"   ,&fTnAssoTracks  ,"JNAssoTracks[NJets]/I");
  fEventTree->Branch("Jtrk1px"        ,&fTtrk1px       ,"Jtrk1px[NJets]/D");
  fEventTree->Branch("Jtrk1py"        ,&fTtrk1py       ,"Jtrk1py[NJets]/D");
  fEventTree->Branch("Jtrk1pz"        ,&fTtrk1pz       ,"Jtrk1pz[NJets]/D");
  fEventTree->Branch("Jtrk2px"        ,&fTtrk2px       ,"Jtrk2px[NJets]/D");
  fEventTree->Branch("Jtrk2py"        ,&fTtrk2py       ,"Jtrk2py[NJets]/D");
  fEventTree->Branch("Jtrk2pz"        ,&fTtrk2pz       ,"Jtrk2pz[NJets]/D");
  fEventTree->Branch("Jtrk3px"        ,&fTtrk3px       ,"Jtrk3px[NJets]/D");
  fEventTree->Branch("Jtrk3py"        ,&fTtrk3py       ,"Jtrk3py[NJets]/D");
  fEventTree->Branch("Jtrk3pz"        ,&fTtrk3pz       ,"Jtrk3pz[NJets]/D");
  fEventTree->Branch("JEcorr"         ,&fTjEcorr       ,"JEcorr[NJets]/D");
  fEventTree->Branch("JeMinDR"        ,&fTjeMinDR      ,"JeMinDR[NJets]/D");
  fEventTree->Branch("JVtxx"          ,&fTjetVtxx      ,"JVtxx[NJets]/D");
  fEventTree->Branch("JVtxy"          ,&fTjetVtxy      ,"JVtxy[NJets]/D");
  fEventTree->Branch("JVtxz"          ,&fTjetVtxz      ,"JVtxz[NJets]/D");
  fEventTree->Branch("JVtxExx"        ,&fTjetVtxExx    ,"JVtxExx[NJets]/D");
  fEventTree->Branch("JVtxEyx"        ,&fTjetVtxEyx    ,"JVtxEyx[NJets]/D");
  fEventTree->Branch("JVtxEyy"        ,&fTjetVtxEyy    ,"JVtxEyy[NJets]/D");
  fEventTree->Branch("JVtxEzy"        ,&fTjetVtxEzy    ,"JVtxEzy[NJets]/D");
  fEventTree->Branch("JVtxEzz"        ,&fTjetVtxEzz    ,"JVtxEzz[NJets]/D");
  fEventTree->Branch("JVtxEzx"        ,&fTjetVtxEzx    ,"JVtxEzx[NJets]/D");
  fEventTree->Branch("JVtxNChi2"      ,&fTjetVtxNChi2  ,"JVtxNChi2[NJets]/D");

  // Tracks:
  fEventTree->Branch("NTracks"        ,&fTntracks      ,"NTracks/I");
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
  fEventTree->Branch("SumEt"          ,&fTSumEt          ,"SumEt/D");
  fEventTree->Branch("ECALSumEt"      ,&fTECALSumEt      ,"ECALSumEt/D");
  fEventTree->Branch("HCALSumEt"      ,&fTHCALSumEt      ,"HCALSumEt/D");
  fEventTree->Branch("ECALEsumx"      ,&fTECALEsumx      ,"ECALEsumx/D");
  fEventTree->Branch("ECALEsumy"      ,&fTECALEsumy      ,"ECALEsumy/D");
  fEventTree->Branch("ECALEsumz"      ,&fTECALEsumz      ,"ECALEsumz/D");
  fEventTree->Branch("ECALMET"        ,&fTECALMET        ,"ECALMET/D");
  fEventTree->Branch("ECALMETPhi"     ,&fTECALMETphi     ,"ECALMETPhi/D");
  fEventTree->Branch("ECALMETEta"     ,&fTECALMETeta     ,"ECALMETEta/D");
  fEventTree->Branch("HCALEsumx"      ,&fTHCALEsumx      ,"HCALEsumx/D");
  fEventTree->Branch("HCALEsumy"      ,&fTHCALEsumy      ,"HCALEsumy/D");
  fEventTree->Branch("HCALEsumz"      ,&fTHCALEsumz      ,"HCALEsumz/D");
  fEventTree->Branch("HCALMET"        ,&fTHCALMET        ,"HCALMET/D");
  fEventTree->Branch("HCALMETPhi"     ,&fTHCALMETphi     ,"HCALMETPhi/D");
  fEventTree->Branch("HCALMETeta"     ,&fTHCALMETeta     ,"HCALMETEta/D");
  fEventTree->Branch("RawMET"         ,&fTRawMET         ,"RawMET/D");
  fEventTree->Branch("RawMETpx"       ,&fTRawMETpx       ,"RawMETpx/D");
  fEventTree->Branch("RawMETpy"       ,&fTRawMETpy       ,"RawMETpy/D");
  fEventTree->Branch("RawMETphi"      ,&fTRawMETphi      ,"RawMETphi/D");
  fEventTree->Branch("MuCorrMET"      ,&fTMuCorrMET      ,"MuCorrMET/D");
  fEventTree->Branch("MuCorrMETpx"    ,&fTMuCorrMETpx    ,"MuCorrMETpx/D");
  fEventTree->Branch("MuCorrMETpy"    ,&fTMuCorrMETpy    ,"MuCorrMETpy/D");
  fEventTree->Branch("MuCorrMETphi"   ,&fTMuCorrMETphi   ,"MuCorrMETphi/D");
  fEventTree->Branch("TCMET"          ,&fTTCMET          ,"TCMET/D");
  fEventTree->Branch("TCMETpx"        ,&fTTCMETpx        ,"TCMETpx/D");
  fEventTree->Branch("TCMETpy"        ,&fTTCMETpy        ,"TCMETpy/D");
  fEventTree->Branch("TCMETphi"       ,&fTTCMETphi       ,"TCMETphi/D");
  fEventTree->Branch("MuJESCorrMET"   ,&fTMuJESCorrMET   ,"MuJESCorrMET/D");
  fEventTree->Branch("MuJESCorrMETpx" ,&fTMuJESCorrMETpx ,"MuJESCorrMETpx/D");
  fEventTree->Branch("MuJESCorrMETpy" ,&fTMuJESCorrMETpy ,"MuJESCorrMETpy/D");
  fEventTree->Branch("MuJESCorrMETphi",&fTMuJESCorrMETphi,"MuJESCorrMETphi/D");
  fEventTree->Branch("PFMET"          ,&fTPFMET          ,"PFMET/D");
  fEventTree->Branch("PFMETpx"        ,&fTPFMETpx        ,"PFMETpx/D");
  fEventTree->Branch("PFMETpy"        ,&fTPFMETpy        ,"PFMETpy/D");
  fEventTree->Branch("PFMETphi"       ,&fTPFMETphi       ,"PFMETphi/D");
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
  fRTMaxeliso      = -999.99;
  fRTMaxeld0       = -999.99;
  fRTMinjpt        = -999.99;
  fRTMaxjeta       = -999.99;
  fRTMinjemfrac    = -999.99;
  fRTMintrkpt      = -999.99;
  fRTMaxtrketa     = -999.99;
  fRTMaxtrknchi2   = -999.99;
  fRTMintrknhits   = -999;
  fRTIsoMuTkDRin   = -999.99;
  fRTIsoMuTkDRout  = -999.99;
  fRTIsoMuTkSeed   = -999.99;
  fRTIsoMuCalDRin  = -999.99;
  fRTIsoMuCalDRout = -999.99;
  fRTIsoMuCalSeed  = -999.99;


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
  fRTMaxeliso      = fMaxeliso;
  fRTMaxeld0       = fMaxeld0;
  fRTMinjpt        = fMinjpt;
  fRTMaxjeta       = fMaxjeta;
  fRTMinjemfrac    = fMinjemfrac;
  fRTMintrkpt      = fMintrkpt;
  fRTMaxtrketa     = fMaxtrketa;
  fRTMaxtrknchi2   = fMaxtrknchi2;
  fRTMintrknhits   = fMintrknhits;

  fRTIsoMuTkDRin   = fIso_MuTkDRin;
  fRTIsoMuTkDRout  = fIso_MuTkDRout;
  fRTIsoMuTkSeed   = fIso_MuTkSeed;
  fRTIsoMuCalDRin  = fIso_MuCalDRin;
  fRTIsoMuCalDRout = fIso_MuCalDRout;
  fRTIsoMuCalSeed  = fIso_MuCalSeed;

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
  resetInt(fTtrigres, 200);
  fTrunnumber   = -999;
  fTeventnumber = -999;
  fTlumisection = -999;
  fTsigprocid   = -999;
  fTweight      = -999.99;
  fTprimvtxx    = -999.99;
  fTprimvtxy    = -999.99;
  fTprimvtxz    = -999.99;
  fTprimvtxxE   = -999.99;
  fTprimvtxyE   = -999.99;
  fTprimvtxzE   = -999.99;
  fTpvtxznchi2  = -999.99;
  fTpvtxntracks = -999;
  fTbeamspotx   = -999.99;
  fTbeamspoty   = -999.99;
  fTbeamspotz   = -999.99;
  fTNCaloTowers = -999;

  fTgoodevent        = 1;
  fTflagmaxmuexc     = 0;
  fTflagmaxelexc     = 0;
  fTflagmaxujetexc = 0;
  fTflagmaxjetexc    = 0;
  fTflagmaxtrkexc    = 0;

  fTnmu   = 0;
  fTneles = 0;
  fTnjets = 0;
  resetInt(fTgoodmu);
  resetDouble(fTmupx);
  resetDouble(fTmupy);
  resetDouble(fTmupz);
  resetDouble(fTmupt);
  resetDouble(fTmue);
  resetDouble(fTmuet);
  resetDouble(fTmueta);
  resetDouble(fTmuphi);
  resetInt(fTmucharge);
  resetDouble(fTmuptsum);
  resetDouble(fTmuetsum);
  resetDouble(fTmuiso);
  resetDouble(fTmueecal);
  resetDouble(fTmuehcal);
  resetDouble(fTmud0bs);
  resetDouble(fTmud0pv);
  resetDouble(fTmud0E);
  resetDouble(fTmudzbs);
  resetDouble(fTmudzpv);
  resetDouble(fTmudzE);
  resetDouble(fTmunchi2);
  resetInt(fTmunglhits);
  resetInt(fTmunmuhits);
  resetInt(fTmuntkhits);
  resetInt(fTmunmatches);
  resetInt(fTmunchambers);
  resetDouble(fTmucalocomp);
  resetDouble(fTmusegmcomp);
  resetInt(fTmutrackermu);
  resetInt(fTmuisGMPT);
  resetInt(fTmuid);
  resetInt(fTmumid);

  resetInt(fTgoodel);
  resetDouble(fTepx);
  resetDouble(fTepy);
  resetDouble(fTepz);
  resetDouble(fTee);
  resetDouble(fTeet);
  resetDouble(fTept);
  resetDouble(fTeeta);
  resetDouble(fTephi);
  resetDouble(fTed0bs);
  resetDouble(fTed0pv);
  resetDouble(fTed0E);
  resetDouble(fTedzbs);
  resetDouble(fTedzpv);
  resetDouble(fTedzE);
  resetDouble(fTenchi2);
  resetDouble(fTeiso);
  resetDouble(fTeptsum);
  resetDouble(fTeetsum);
  resetInt(fTecharge);
  resetInt(fTeInGap);
  resetInt(fTeEcalDriven);
  resetInt(fTeTrackerDriven);
  resetInt(fTeBasicClustersSize);
  resetDouble(fTefbrem);
  resetDouble(fTeHcalOverEcal);
  resetDouble(fTeE5x5);
  resetDouble(fTeE2x5Max);
  resetDouble(fTeSigmaIetaIeta);
  resetDouble(fTeDeltaPhiSeedClusterAtCalo);
  resetDouble(fTeDeltaEtaSeedClusterAtCalo);
  resetDouble(fTeDeltaPhiSuperClusterAtVtx);
  resetDouble(fTeDeltaEtaSuperClusterAtVtx);
  resetDouble(fTecaloenergy);
  resetDouble(fTtrkmomatvtx);
  resetDouble(fTeESuperClusterOverP);

  resetInt(fTgoodjet, fMaxnjets);
  resetDouble(fTjpx,  fMaxnjets);
  resetDouble(fTjpy,  fMaxnjets);
  resetDouble(fTjpz,  fMaxnjets);
  resetDouble(fTje,   fMaxnjets);
  resetDouble(fTjet,  fMaxnjets);
  resetDouble(fTjpt,  fMaxnjets);
  resetDouble(fTjeta, fMaxnjets);
  resetDouble(fTjphi, fMaxnjets);
  resetDouble(fTjemfrac, fMaxnjets);
  resetDouble(fTjID_HPD, fMaxnjets);
  resetDouble(fTjID_RBX, fMaxnjets);
  resetDouble(fTjID_n90Hits, fMaxnjets);
  resetDouble(fTjID_SubDet1, fMaxnjets);
  resetDouble(fTjID_SubDet2, fMaxnjets);
  resetDouble(fTjID_SubDet3, fMaxnjets);
  resetDouble(fTjID_SubDet4, fMaxnjets);
  resetDouble(fTjID_resEMF,  fMaxnjets);
  resetDouble(fTjID_HCALTow, fMaxnjets);
  resetDouble(fTjID_ECALTow, fMaxnjets);
  resetDouble(fTbTagProb, fMaxnjets);
  resetDouble(fTChfrac,   fMaxnjets);
  resetInt(fTnAssoTracks, fMaxnjets);
  resetDouble(fTtrk1px, fMaxnjets);
  resetDouble(fTtrk1py, fMaxnjets);
  resetDouble(fTtrk1pz, fMaxnjets);
  resetDouble(fTtrk2px, fMaxnjets);
  resetDouble(fTtrk2py, fMaxnjets);
  resetDouble(fTtrk2pz, fMaxnjets);
  resetDouble(fTtrk3px, fMaxnjets);
  resetDouble(fTtrk3py, fMaxnjets);
  resetDouble(fTtrk3pz, fMaxnjets);
  resetDouble(fTjeMinDR, fMaxnjets);
  resetDouble(fTjetVtxx, fMaxnjets);
  resetDouble(fTjetVtxy, fMaxnjets);
  resetDouble(fTjetVtxz, fMaxnjets);
  resetDouble(fTjetVtxExx, fMaxnjets);
  resetDouble(fTjetVtxEyx, fMaxnjets);
  resetDouble(fTjetVtxEyy, fMaxnjets);
  resetDouble(fTjetVtxEzy, fMaxnjets);
  resetDouble(fTjetVtxEzz, fMaxnjets);
  resetDouble(fTjetVtxEzx, fMaxnjets);
  resetDouble(fTjetVtxNChi2, fMaxnjets);

  resetDouble(fJUNC_px_match, fMaxnjets);
  resetDouble(fJUNC_py_match, fMaxnjets);
  resetDouble(fJUNC_pz_match, fMaxnjets);

  resetInt(fTgoodtrk,  fMaxntrks);
  resetDouble(fTtrkpt, fMaxntrks);
  resetDouble(fTtrketa, fMaxntrks);
  resetDouble(fTtrkphi, fMaxntrks);
  resetDouble(fTtrknchi2, fMaxntrks);
  resetDouble(fTtrknhits, fMaxntrks);

  resetInt(fTeIDTight);
  resetInt(fTeIDLoose);
  resetInt(fTeIDRobustTight);
  resetInt(fTeIDRobustLoose);

  fTTrkPtSumx       = -999.99;
  fTTrkPtSumy       = -999.99;
  fTTrkPtSum        = -999.99;
  fTTrkPtSumphi     = -999.99;
  fTSumEt           = -999.99;
  fTECALSumEt       = -999.99;
  fTHCALSumEt       = -999.99;
  fTECALEsumx       = -999.99;
  fTECALEsumy       = -999.99;
  fTECALEsumz       = -999.99;
  fTECALMET         = -999.99;
  fTECALMETphi      = -999.99;
  fTECALMETeta      = -999.99;
  fTHCALEsumx       = -999.99;
  fTHCALEsumy       = -999.99;
  fTHCALEsumz       = -999.99;
  fTHCALMET         = -999.99;
  fTHCALMETphi      = -999.99;
  fTHCALMETeta      = -999.99;
  fTRawMET          = -999.99;
  fTRawMETpx        = -999.99;
  fTRawMETpy        = -999.99;
  fTRawMETphi       = -999.99;
  fTMuCorrMET       = -999.99;
  fTMuCorrMETpx     = -999.99;
  fTMuCorrMETpy     = -999.99;
  fTMuCorrMETphi    = -999.99;
  fTTCMET           = -999.99;
  fTTCMETpx         = -999.99;
  fTTCMETpy         = -999.99;
  fTTCMETphi        = -999.99;
  fTMuJESCorrMET    = -999.99;
  fTMuJESCorrMETpx  = -999.99;
  fTMuJESCorrMETpy  = -999.99;
  fTMuJESCorrMETphi = -999.99;
  fTPFMET           = -999.99;
  fTPFMETpx         = -999.99;
  fTPFMETpy         = -999.99;
  fTPFMETphi        = -999.99;
}

// Method for calculating the muon isolation
vector<double> NTupleProducer::calcMuIso(const reco::Muon *Mu, const edm::Event& iEvent){
  // Track pT sum:
  edm::Handle<VertexCollection> vertices;
  iEvent.getByLabel(fVertexTag, vertices);
  const reco::Vertex *primVtx = &(*(vertices.product()))[0]; // Just take first vertex ...
  // Compare against eta and phi of trackermuon:
  double mupt  = Mu->innerTrack()->pt();
  double mueta = Mu->innerTrack()->eta();
  double muphi = Mu->innerTrack()->phi();
  double ptsum(0.);
  for(vector<TrackBaseRef>::const_iterator itk = primVtx->tracks_begin(); itk!=primVtx->tracks_end(); ++itk){
    if((*itk)->pt() < fIso_MuTkSeed) continue;
    double eta = (*itk)->eta();
    double phi = (*itk)->phi();
    double DR = reco::deltaR(eta, phi, mueta, muphi);
    // if(DR <= 0.) DR = 0.001; // Why do I need this?
    if(DR < fIso_MuTkDRin) continue;
    if(DR > fIso_MuTkDRout) continue;
    ptsum += (*itk)->pt();
  }
  // Subtract Muon pT in case inner cone is set to 0.0 radius
  // This doesn't seem to work very well...
  // Muon track not part of primVtx tracks?1
  if(fIso_MuTkDRin == 0.0) ptsum -= mupt;

  // Calo ET sum:
  edm::Handle<CaloTowerCollection> calotowers;
  iEvent.getByLabel(fCalTowTag,calotowers);
  double etsum(0.);
  for(CaloTowerCollection::const_iterator itow = calotowers->begin(); itow!=calotowers->end(); ++itow){
    double eta = itow->eta();
    if(itow->energy()/cosh(eta) < fIso_MuCalSeed) continue;
    double phi = itow->phi();
    double DR = reco::deltaR(eta, phi, mueta, muphi);
    if(DR <= 0.) DR = 0.001;
    if(DR < fIso_MuCalDRin) continue;
    if(DR > fIso_MuCalDRout) continue;
    etsum += itow->energy()/cosh(eta);
  }
  vector<double> res;
  res.push_back((ptsum+etsum)/mupt);
  res.push_back(ptsum);
  res.push_back(etsum);
  return res;
}

// Method for calculating the muon isolation
vector<double> NTupleProducer::calcMuIso2(const reco::Muon *Mu, const edm::Event& iEvent){
  // Track pT sum:
  edm::Handle<TrackCollection> tracks;
  iEvent.getByLabel(fTrackTag, tracks);
  // Compare against eta and phi of trackermuon:
  double mupt  = Mu->innerTrack()->pt();
  double mueta = Mu->innerTrack()->eta();
  double muphi = Mu->innerTrack()->phi();
  double ptsum(0.);
  for(TrackCollection::const_iterator itk = tracks->begin(); itk!=tracks->end(); ++itk){
    if(itk->pt() < fIso_MuTkSeed) continue;
    double eta = itk->eta();
    double phi = itk->phi();
    double DR = reco::deltaR(eta, phi, mueta, muphi);
    // if(DR <= 0.) DR = 0.001; // Why do I need this?
    if(DR < fIso_MuTkDRin) continue;
    if(DR > fIso_MuTkDRout) continue;
    ptsum += itk->pt();
  }
  // Subtract Muon pT in case inner cone is set to 0.0 radius
  // This doesn't seem to work very well...
  // Muon track not part of primVtx tracks?
  if(fIso_MuTkDRin == 0.0) ptsum -= mupt;

  // Calo ET sum:
  edm::Handle<CaloTowerCollection> calotowers;
  iEvent.getByLabel(fCalTowTag,calotowers);
  double etsum = 0.;
  for(CaloTowerCollection::const_iterator itow = calotowers->begin(); itow!=calotowers->end(); ++itow){
    double eta = itow->eta();
    if(itow->energy()/cosh(eta) < fIso_MuCalSeed) continue;
    double phi = itow->phi();
    double DR = reco::deltaR(eta, phi, mueta, muphi);
    if(DR <= 0.) DR = 0.001;
    if(DR < fIso_MuCalDRin) continue;
    if(DR > fIso_MuCalDRout) continue;
    etsum += itow->energy()/cosh(eta);
  }
  vector<double> res;
  res.push_back((ptsum+etsum)/mupt);
  res.push_back(ptsum);
  res.push_back(etsum);
  return res;
}

// Method for calculating the muon isolation
vector<double> NTupleProducer::calcMuIso3(const reco::Muon *Mu, const edm::Event& iEvent, edm::Ref<reco::MuonCollection> muonRef){
  double mupt = Mu->innerTrack()->pt();

  // Read MuIsoDeposits
  // Tracker:
  edm::Handle<edm::ValueMap<reco::IsoDeposit> > IsoDepTkValueMap;
  iEvent.getByLabel(fMuIsoDepTkTag, IsoDepTkValueMap);
  const edm::ValueMap<reco::IsoDeposit> &TkDepMap = *IsoDepTkValueMap.product();
  // ECAL:
  edm::Handle<edm::ValueMap<reco::IsoDeposit> > IsoDepECValueMap;
  iEvent.getByLabel(fMuIsoDepECTag, IsoDepECValueMap);
  const edm::ValueMap<reco::IsoDeposit> &ECDepMap = *IsoDepECValueMap.product();
  // HCAL:
  edm::Handle<edm::ValueMap<reco::IsoDeposit> > IsoDepHCValueMap;
  iEvent.getByLabel(fMuIsoDepHCTag, IsoDepHCValueMap);
  const edm::ValueMap<reco::IsoDeposit> &HCDepMap = *IsoDepHCValueMap.product();

  const reco::IsoDeposit TkDep = TkDepMap[muonRef];
  const reco::IsoDeposit ECDep = ECDepMap[muonRef];
  const reco::IsoDeposit HCDep = HCDepMap[muonRef];
  double ptsum = TkDep.depositWithin(fIso_MuTkDRout);
  double ecsum = ECDep.depositWithin(fIso_MuCalDRout);
  double hcsum = HCDep.depositWithin(fIso_MuCalDRout);
  double etsum = ecsum + hcsum;
  double muiso = (ptsum + etsum)/mupt;
  vector<double> res;
  res.push_back(muiso);
  res.push_back(ptsum);
  res.push_back(etsum);
  return res;
}

// Method for calculating the electron isolation
vector<double> NTupleProducer::calcElIso(const reco::GsfElectron *El, const edm::Event& iEvent){
  // Track pT sum:
  edm::Handle<VertexCollection> vertices;
  iEvent.getByLabel(fVertexTag, vertices);
  const reco::Vertex *primVtx = &(*(vertices.product()))[0]; // Just take first vertex ...
  // Compare against eta and phi of electron:
  double elpt  = El->pt();
  double eleta = El->eta();
  double elphi = El->phi();
  double ptsum(0.);
  for(vector<TrackBaseRef>::const_iterator itk = primVtx->tracks_begin(); itk!=primVtx->tracks_end(); ++itk){
    if((*itk)->pt() < fIso_MuTkSeed) continue;
    double eta = (*itk)->eta();
    double phi = (*itk)->phi();
    double DR = reco::deltaR(eta, phi, eleta, elphi);
    if(DR <= 0.) DR = 0.001; // Why do I need this?
    if(DR < fIso_MuTkDRin) continue;
    if(DR > fIso_MuTkDRout) continue;
    ptsum += (*itk)->pt();
  }
  // Subtract electron pT in case inner cone is set to 0.0 radius
  if(fIso_MuTkDRin == 0.0) ptsum -= elpt;

  // SC ET sum:
  edm::Handle<SuperClusterCollection> superclusters;
  iEvent.getByLabel(fSCTag,superclusters);
  double etsum = 0.;
  for(SuperClusterCollection::const_iterator isc = superclusters->begin(); isc!=superclusters->end(); ++isc){
    double eta = isc->eta();
    if(isc->energy()/cosh(eta) < fIso_MuCalSeed) continue;
    double phi = isc->phi();
    double DR = reco::deltaR(eta, phi, eleta, elphi);
    if(DR <= 0.) DR = 0.001; //FIXME:???
    if(DR < fIso_MuCalDRin) continue;
    if(DR > fIso_MuCalDRout) continue;
    etsum += isc->energy()/cosh(eta);
  }
  if(fIso_MuCalDRin == 0.0) etsum -= El->caloEnergy();
  // // Calo ET sum:
  // edm::Handle<CaloTowerCollection> calotowers;
  // iEvent.getByLabel(fCalTowTag,calotowers);
  // double etsum = 0.;
  // for(CaloTowerCollection::const_iterator itow = calotowers->begin(); itow!=calotowers->end(); ++itow){
  //  double eta = itow->eta();
  //  if(itow->energy()/cosh(eta) < fIso_MuCalSeed) continue;
  //  double phi = itow->phi();
  //  double DR = reco::deltaR(eta, phi, mueta, muphi);
  //  if(DR <= 0.) DR = 0.001;
  //  if(DR < fIso_MuCalDRin) continue;
  //  if(DR > fIso_MuCalDRout) continue;
  //  etsum += itow->energy()/cosh(eta);
  // }
  // if(fIso_MuCalDRin == 0.0) etsum -= El->caloEnergy();

  vector<double> res;
  res.push_back((ptsum+etsum)/elpt);
  res.push_back(ptsum);
  res.push_back(etsum);
  return res;
}

// Method for matching of muon candidates
vector<int> NTupleProducer::matchMuCand(const reco::Muon *Mu, const edm::Event& iEvent){
  if(fIsRealData){
    edm::LogWarning("NTP") << "@SUB=matchMuCand"
                           << "Trying to access generator info on real data...";
    vector<int> fakeres;
    fakeres.push_back(0); fakeres.push_back(0); fakeres.push_back(0);
    return fakeres;
  }
  int muid(0);
  int mumid(0);
  int match(0);
  vector<int> res;
  double mupt = Mu->globalTrack()->pt();
  double mueta = Mu->globalTrack()->eta();
  double muphi = Mu->globalTrack()->phi();
  int mucharge = Mu->globalTrack()->charge();
  edm::Handle<GenParticleCollection> genparts;
  iEvent.getByLabel(fGenPartTag, genparts);
  GenParticleCollection::const_iterator gpart;
  const GenParticle *GenMu = new GenParticle();
  bool matched = false;
  // Try to match the muon candidate to a generator muon
  double mindr(999.99);
  for(gpart = genparts->begin(); gpart != genparts->end(); gpart++){
    // Select muons, pions and kaons
    int id = abs(gpart->pdgId());
    if(id!=13 && id!=211 && id!=321) continue;
    // Restrict to cone of 0.1 in DR around mu
    double dr = reco::deltaR(gpart->eta(), gpart->phi(), mueta, muphi);
    if(dr > 0.1) continue;
    // Select same charge
    int ch = gpart->charge();
    if(ch!=mucharge) continue;
    // Restrict to pt match within a factor of 2
    double ndpt = fabs(gpart->pt()-mupt)/gpart->pt();
    if(ndpt > 2.) continue;

    // Minimize DeltaR
    if(dr > mindr) continue;
    mindr = dr;
    matched = true;
    GenMu = &(*gpart);
  }

  // Fill generator information in case match was successful
  if(matched){
    match = 1;
    muid  = GenMu->pdgId();
    // Determine mother object of muon candidate:
    // (Make sure that the mother has a different PDG ID)
    mumid = GenMu->mother()->pdgId();
    // Distinguish the three cases: Mu, Pi, K
    if(abs(mumid) != abs(muid));
    else{
      int tempid = abs(muid);
      const Candidate *temppart = GenMu->mother()->mother();
      while(tempid == abs(muid)){
        tempid = temppart->pdgId();
        temppart = temppart->mother();
      }
      mumid = tempid;
    }
    res.push_back(match);
    res.push_back(muid);
    res.push_back(mumid);
  }
  else{
    match = 0;
    res.push_back(match);
    res.push_back(muid);
    res.push_back(mumid);
  }
  return res;
}

// Method for sorting the muons
vector<const reco::Muon*> NTupleProducer::sortMus(vector<const reco::Muon*> Mus){
  /* - This will arrange the vector of muons such that they are ordered in pT   */
  unsigned int size = Mus.size();
  vector<double> v_mupt;
  vector<const reco::Muon*> v_mu;
  // Fill mupt vector with dummies
  for(size_t i = 0; i < size; ++i){
    v_mupt.push_back(-888.88);
    const reco::Muon *mu = new reco::Muon();
    v_mu.push_back(mu);
  }
  double mupt(0.);
  for(size_t i = 0; i < size; ++i){ // Loop over muons
    mupt = Mus[i]->globalTrack()->pt();
    for(size_t j = 0; j < size; ++j){ // Loop over container slots
      if(mupt > v_mupt[j]){
        for(size_t k = size-1; k > j; k--){ // Shift lower cont. entries
          v_mupt[k] = v_mupt[k-1];
          v_mu[k] = v_mu[k-1];
        }
        v_mupt[j] = mupt;
        v_mu[j] = Mus[i];
        break;
      }
    }
  }
  return v_mu;
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

vector<const reco::Track*> NTupleProducer::FindAssociatedTracks(const reco::Jet *jet, const reco::TrackCollection *tracks){
  //returns a list of tracks associated to a given jet
  vector<const reco::Track*> AssociatedTracks;
  double jeteta=jet->eta();
  double jetphi=jet->phi();
  for(TrackCollection::const_iterator itk = tracks->begin(); itk!=tracks->end(); ++itk){
    // Use safe reco::deltaR function (phi should be between -pi and pi...)
    if ( reco::deltaR( jeteta, jetphi, itk->eta(), itk->phi() ) < 0.5 )
      AssociatedTracks.push_back(&(*itk));
  }
  return AssociatedTracks;
}

//define this as a plug-in
DEFINE_FWK_MODULE(NTupleProducer);
