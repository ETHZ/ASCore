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
// $Id: NTupleProducer.cc,v 1.8 2009/10/02 14:34:01 sordini Exp $
//
//

#include "DiLeptonAnalysis/NTupleProducer/interface/NTupleProducer.h"


NTupleProducer::NTupleProducer(const edm::ParameterSet& iConfig){
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
	fTriggerTag     = iConfig.getParameter<edm::InputTag>("tag_triggers");

	// Event Selection
	fMinmupt        = iConfig.getUntrackedParameter<double>("sel_minmupt");
	fMaxmueta       = iConfig.getUntrackedParameter<double>("sel_maxmueta");
	fMinelpt        = iConfig.getUntrackedParameter<double>("sel_minelpt");
	fMaxeleta       = iConfig.getUntrackedParameter<double>("sel_maxeleta");
	fMaxeliso       = iConfig.getUntrackedParameter<double>("sel_maxeliso");
	fMaxeld0        = iConfig.getUntrackedParameter<double>("sel_maxeld0");
	fMinjpt         = iConfig.getUntrackedParameter<double>("sel_minjpt");
	fMaxjeta        = iConfig.getUntrackedParameter<double>("sel_maxjeta");
	fMinjemfrac     = iConfig.getUntrackedParameter<double>("sel_minjemfrac");

	// Isolation Parameters
	fIso_MuTkDRin   = iConfig.getUntrackedParameter<double>("iso_MuTkDRin");
	fIso_MuTkDRout  = iConfig.getUntrackedParameter<double>("iso_MuTkDRout");
	fIso_MuTkSeed   = iConfig.getUntrackedParameter<double>("iso_MuTkSeed");
	fIso_MuCalDRin  = iConfig.getUntrackedParameter<double>("iso_MuCalDRin");
	fIso_MuCalDRout = iConfig.getUntrackedParameter<double>("iso_MuCalDRout");
	fIso_MuCalSeed  = iConfig.getUntrackedParameter<double>("iso_MuCalSeed");

	cout << "---------------------------------" << endl;
	cout << " ==> NTupleProducer Constructor ..." << endl;
	cout << endl;
	cout << "  Input Tags:" << endl;
	cout << "    fMuonTag        = " << fMuonTag.label()        << endl;
	cout << "    fElectronTag    = " << fElectronTag.label()    << endl;
	cout << "    fEleIsoTkTag    = " << fEleIsoTkTag.label()    << endl;
	cout << "    fEleIsoECTag    = " << fEleIsoECTag.label()    << endl;
	cout << "    fEleIsoHCTag    = " << fEleIsoHCTag.label()    << endl;
	cout << "    fEleIsoDepTkTag = " << fEleIsoDepTkTag.label() << endl;
	cout << "    fEleIsoDepECTag = " << fEleIsoDepECTag.label() << endl;
	cout << "    fEleIsoDepHCTag = " << fEleIsoDepHCTag.label() << endl;
	cout << "    fMuIsoDepTkTag  = " << fMuIsoDepTkTag.label()  << endl;
	cout << "    fMuIsoDepECTag  = " << fMuIsoDepECTag.label()  << endl;
	cout << "    fMuIsoDepHCTag  = " << fMuIsoDepHCTag.label()  << endl;
	cout << "    fSCTag          = " << fSCTag.label()          << endl;
	cout << "    fJetTag         = " << fJetTag.label()         << endl;
	cout << "    fMET1Tag        = " << fMET1Tag.label()        << endl;
	cout << "    fMET2Tag        = " << fMET2Tag.label()        << endl;
	cout << "    fMET3Tag        = " << fMET3Tag.label()        << endl;
	cout << "    fMET4Tag        = " << fMET4Tag.label()        << endl;
	cout << "    fMET5Tag        = " << fMET5Tag.label()        << endl;
	cout << "    fVertexTag      = " << fVertexTag.label()      << endl;
	cout << "    fTrackTag       = " << fTrackTag.label()       << endl;
	cout << "    fCalTowTag      = " << fCalTowTag.label()      << endl;
	cout << "    fGenPartTag     = " << fGenPartTag.label()     << endl;
	cout << endl;
	cout << "  Event Selection Parameters:" << endl;
	cout << "    fMinmupt        = " << fMinmupt    << endl;
	cout << "    fMaxmueta       = " << fMaxmueta   << endl;
	cout << "    fMinelpt        = " << fMinelpt    << endl;
	cout << "    fMaxeleta       = " << fMaxeleta   << endl;
	cout << "    fMaxeliso       = " << fMaxeliso   << endl;
	cout << "    fMaxeld0        = " << fMaxeld0    << endl;
	cout << "    fMinjpt         = " << fMinjpt     << endl;
	cout << "    fMaxjeta        = " << fMaxjeta    << endl;
	cout << "    fMinjemfrac     = " << fMinjemfrac << endl;
	cout << endl;
	cout << "  Isolation Parameters:" << endl;
	cout << "    fIso_MuTkDRin   = " << fIso_MuTkDRin   << endl;
	cout << "    fIso_MuTkDRout  = " << fIso_MuTkDRout  << endl;
	cout << "    fIso_MuTkSeed   = " << fIso_MuTkSeed   << endl;
	cout << "    fIso_MuCalDRin  = " << fIso_MuCalDRin  << endl;
	cout << "    fIso_MuCalDRout = " << fIso_MuCalDRout << endl;
	cout << "    fIso_MuCalSeed  = " << fIso_MuCalSeed  << endl;
	cout << endl;
	cout << "---------------------------------" << endl;
}


NTupleProducer::~NTupleProducer(){
}

// Method called once for each event
void NTupleProducer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){
	fNTotEvents++;
	using namespace edm;
	using namespace std;
	using namespace reco;
	using reco::MuonCollection;
	using reco::CaloJetCollection;
	using reco::JetTagCollection;

	// Reset all the tree variables
	resetTree();

	// Get the collections:
	Handle<MuonCollection> muons;
	iEvent.getByLabel(fMuonTag,muons); // 'muons'

	Handle<GsfElectronCollection> electrons;
	iEvent.getByLabel(fElectronTag, electrons); // 'gsfElectrons'

	// Jets and Jet Correctors
	Handle<CaloJetCollection> jets;
	iEvent.getByLabel(fJetTag,jets); // 'sisCone5CaloJets'
	const JetCorrector* L2JetCorrector = JetCorrector::getJetCorrector ("L2RelativeJetCorrectorSC5Calo",iSetup);
	const JetCorrector* L3JetCorrector = JetCorrector::getJetCorrector ("L3AbsoluteJetCorrectorSC5Calo",iSetup);
	//
	reco::helper::JetID jetID;

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

	// Dump HLT trigger bits
	Handle<TriggerResults> triggers;
	iEvent.getByLabel(fTriggerTag, triggers);
	TriggerResults tr = *triggers;
	if(fFirstevent){
		fFirstevent = false;
		vector<string> triggernames;
		Service<service::TriggerNamesService> tns;
		tns->getTrigPaths(*triggers, triggernames);
		for( unsigned int i = 0; i < tr.size(); i++ ){
			fHtrigstat->GetXaxis()->SetBinLabel(i+1, TString(triggernames[i]));
			// cout << "i " << i << " : " << tr[i].accept() << " : " << triggernames[i] << endl;
		}
		// Add a warning about the shift between trigger bits and bin numbers:
		fHtrigstat->GetXaxis()->SetBinLabel(tr.size()+2, "Bin#=Bit#+1");
	}
	for( unsigned int i = 0; i < tr.size(); i++ ){
		if(tr[i].accept())	fHtrigstat->Fill(i);
		fTtrigres[i] = tr[i].accept() ? 1:0;
	}

////////////////////////////////////////////////////////////////////////////////
// Event Selection /////////////////////////////////////////////////////////////
	bool acceptEvent = true;

	fTrunnumber = iEvent.id().run();
	fTeventnumber = iEvent.id().event();
	fTlumisection = iEvent.luminosityBlock();

	fTweight = 1.0; // To be filled at some point?

	// Save position of primary vertex
	fTprimvtxx   = primVtx->x();
	fTprimvtxy   = primVtx->y();
	fTprimvtxz   = primVtx->z();
	fTprimvtxxE  = primVtx->xError();
	fTprimvtxyE  = primVtx->yError();
	fTprimvtxzE  = primVtx->zError();
	fTpvtxznchi2 = primVtx->normalizedChi2();

	// Save position of beamspot
	fTbeamspotx = (beamSpot.position()).x();
	fTbeamspoty = (beamSpot.position()).y();
	fTbeamspotz = (beamSpot.position()).z();

	// Loop over MuonCollection to count muons
	int mi(-1), mqi(-1); // index of all muons and qualified muons respectively
	for(MuonCollection::const_iterator Mit = muons->begin(); Mit != muons->end(); ++Mit){
		// Check if maximum number of electrons is exceeded already:
		if(mqi > 19){
			cout << "NTupleSkeleton::analyze() ==> Warning: Maximum number of muons exceeded..." << endl;
			break;
		}
		mi++;
		// Muon preselection:
		if(!(Mit->isGlobalMuon())) continue;
		if(Mit->globalTrack()->pt() < fMinmupt) continue;
		if(TMath::Abs(Mit->globalTrack()->eta()) > fMaxmueta) continue;
		mqi++;
		const Muon *cm = Mit->clone();
		edm::Ref<reco::MuonCollection> muonRef(muons, mi);

		// Dump muon properties in tree variables
		fTmupx[mqi]     = Mit->globalTrack()->px();
		fTmupy[mqi]     = Mit->globalTrack()->py();
		fTmupz[mqi]     = Mit->globalTrack()->pz();
		fTmue[mqi]      = Mit->energy();
		fTmuet[mqi]     = Mit->et();
		fTmupt[mqi]     = Mit->globalTrack()->pt();
		fTmueta[mqi]    = Mit->globalTrack()->eta();
		fTmuphi[mqi]    = Mit->globalTrack()->phi();
		fTmucharge[mqi] = Mit->charge();

		vector<double> MuIso = calcMuIso(&(*Mit), iEvent);
		fTmuiso[mqi]   = MuIso[0];
		fTmuptsum[mqi] = MuIso[1];
		fTmuetsum[mqi] = MuIso[2];

		const reco::IsoDeposit ECDep = ECDepMap[muonRef];
		const reco::IsoDeposit HCDep = HCDepMap[muonRef];

		fTmueecal[mqi] = ECDep.candEnergy();
		fTmuehcal[mqi] = HCDep.candEnergy();

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
		fTmusegmcomp[mqi] = muon::segmentCompatibility(*cm);
		fTmutrackermu[mqi] = Mit->isTrackerMuon() ? 1:0;
		fTmuisGMPT[mqi] = muon::isGoodMuon(*cm, muon::GlobalMuonPromptTight) ? 1:0;

		// Matching
		vector<int> MuMatch = matchMuCand(&(*Mit), iEvent);
		if( MuMatch[0] ){
			fTmuid[mqi]  = MuMatch[1];
			fTmumid[mqi] = MuMatch[2];
		}
	}
	fTnmu = mqi+1;

	// Read Electrons
	int eqi(-1); // counts # of qualified electrons
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

		// Loop over electrons
		for( unsigned int i = 0; i < electrons->size(); i++){
			// Check if maximum number of electrons is exceeded already:
			if(eqi > 19){
				cout << "NTupleProducer::analyze() ==> Warning: Maximum number of electrons exceeded..." << endl;
				break;
			}
			Ref<GsfElectronCollection> electronRef(electrons,i);

			// Electron preselection:			
			const GsfElectron *El = &(*(electrons.product()))[i];
			if(El->pt() < fMinelpt) continue;
			if(fabs(El->eta()) > fMaxeleta) continue;
			if(fabs(El->gsfTrack()->dxy(beamSpot.position())) > fMaxeld0) continue;

			// Read EleIsolation
			Handle<edm::ValueMap<double> > eIsoTkValueMap;
			iEvent.getByLabel(fEleIsoTkTag, eIsoTkValueMap);
			const ValueMap<double> &eTkIsoMap = *eIsoTkValueMap.product();
			Handle<edm::ValueMap<double> > eIsoECValueMap;
			iEvent.getByLabel(fEleIsoECTag, eIsoECValueMap);
			const ValueMap<double> &eECIsoMap = *eIsoECValueMap.product();
			Handle<edm::ValueMap<double> > eIsoHCValueMap;
			iEvent.getByLabel(fEleIsoHCTag, eIsoHCValueMap);
			const ValueMap<double> &eHCIsoMap = *eIsoHCValueMap.product();

			double eTkIso = eTkIsoMap[electronRef];
			double eECIso = eECIsoMap[electronRef];
			double eHCIso = eHCIsoMap[electronRef];
			double eliso  = (eTkIso+eECIso+eHCIso)/El->pt();
			if(eliso > fMaxeliso) continue;

			bool isGap = El->isGap(); 
			bool EcalDriven = El-> isEcalDriven();	
			bool TrackerDriven = El->isTrackerDriven();



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
			fTecharge[eqi] = El->charge();
			fTeInGap[eqi]  = isGap ? 1:0;  // no bug... //*** Beni
			fTeEcalDriven[eqi] = EcalDriven ? 1:0;
			fTeTrackerDriven[eqi] = TrackerDriven ? 1:0;
			fTeBasicClustersSize[eqi] = El->basicClustersSize();
			fTefbrem[eqi] = El->fbrem();
			fTeHcalOverEcal[eqi] = El->hcalOverEcal();                               
			fTeE5x5[eqi] = El->e5x5();                                           
			fTeE2x5Max[eqi] = El->e2x5Max();                                           
			fTeSigmaIetaIeta[eqi] = El->sigmaIetaIeta();                         
			fTeDeltaPhiSeedClusterAtCalo[eqi] = El->deltaPhiSeedClusterTrackAtCalo(); 
			fTeDeltaPhiSuperClusterAtVtx[eqi] = El->deltaPhiSuperClusterTrackAtVtx(); 
			fTeESuperClusterOverP[eqi] = El-> eSuperClusterOverP();              
			fTeDeltaEtaSeedClusterAtCalo[eqi] = El->deltaEtaSeedClusterTrackAtCalo(); 

			fTeID[eqi][0] = eIDmapT[electronRef]  ? 1:0;
			fTeID[eqi][1] = eIDmapL[electronRef]  ? 1:0;
			fTeID[eqi][2] = eIDmapRT[electronRef] ? 1:0;
			fTeID[eqi][3] = eIDmapRL[electronRef] ? 1:0;
		}
	}
	fTneles = eqi+1;


	int itag(-1); 
	// Apply L2 and L3 JetCorrections
	vector<const reco::CaloJet*> CorrJets; // Contains all (corrected) jets
	CaloJetCollection::const_iterator Jit;
	for(Jit = jets->begin(); Jit != jets->end(); ++Jit){// Loop over uncorr. jets
		CaloJet *j = Jit->clone();
	//I save here px, py, pz for uncorrected jets
	//used later for matching jets with jets in the jet-tracks associator collections
		itag++;
		fTUNC_px_match[itag]=j->px();
		fTUNC_py_match[itag]=j->py();
		fTUNC_pz_match[itag]=j->pz();
	//energy corrections
		j->scaleEnergy(L2JetCorrector->correction(j->p4()));
		j->scaleEnergy(L3JetCorrector->correction(j->p4()));
		const CaloJet *cj = j->clone();
		CorrJets.push_back(cj);
	}

	// Determine qualified jets
	int jqi(-1); // counts # of qualified jets
	int itagcorr(-1); // counts # of jets for b-tagging (this index will talk to the uncorrected jets)
	for(size_t i = 0; i < CorrJets.size(); ++i){// Loop over corr. jets
		// Check if maximum number of jets is exceeded already:
		if(jqi > 19){
			cout << "NTupleSkeleton::analyze() ==> Warning: Maximum number of jets exceeded..." << endl;
			break;
		}
		itagcorr++;
		CaloJet *jet = CorrJets[i]->clone();
		// Jet preselection
		if(jet->pt() < fMinjpt) continue;
		if(TMath::Abs(jet->eta()) > fMaxjeta) continue;
		if(jet->emEnergyFraction() < fMinjemfrac) continue;
		jqi++;
		//calculate the DR wrt the closest electron
		double ejDRmin= 10.;
		//when no electrons in the event, DR will be 10.
		for(int j = 0; j < fTneles; j++){
			double ejDR = GetDeltaR(jet->eta(), fTeeta[j], jet->phi(), fTephi[j]);
			if(ejDR<ejDRmin) ejDRmin = ejDR;
		}
		fTjeMinDR[jqi]=ejDRmin;
		//jetID variables
		jetID.calculate(iEvent,*jet);
		//fill the b-tagging probability
		for (unsigned int i = 0; i < jetsAndProbs->size(); i++){
			if (fabs( (fTUNC_px_match[itagcorr] - (*jetsAndProbs)[i].first->px())/fTUNC_px_match[itagcorr]) < 0.00001 && 
				fabs( (fTUNC_py_match[itagcorr] - (*jetsAndProbs)[i].first->py())/fTUNC_py_match[itagcorr]) < 0.00001 &&
			fabs( (fTUNC_pz_match[itagcorr] - (*jetsAndProbs)[i].first->pz())/fTUNC_pz_match[itagcorr]) < 0.00001)  {
				fTbTagProb[jqi]=(*jetsAndProbs)[i].second;
				break;
			}
		}
		vector<const reco::Track*> AssociatedTracks = FindAssociatedTracks(&(*jet), tracks.product());
		vector<TransientTrack> AssociatedTTracks;
		if(fabs(jet->eta())<2.9){ //when the cone of dR=0.5 around the jet is (at least partially) inside the tracker acceptance
		//tmp variables for vectorial sum of pt of tracks
		double pXtmp=0.;
		double pYtmp=0.;
		for(size_t t = 0; t < AssociatedTracks.size(); ++t){
		  AssociatedTTracks.push_back(theB->build(AssociatedTracks[t]));
		  if(AssociatedTracks[t]->normalizedChi2()<10. && AssociatedTracks[t]->numberOfValidHits()>10 && AssociatedTracks[t]->pt()>1.){
		    pXtmp+=AssociatedTracks[t]->px();
		    pYtmp+=AssociatedTracks[t]->py();
		  }
		}
		fTChfrac[jqi]=sqrt(pow(pXtmp,2)+pow(pYtmp,2))/jet->pt();
		} else {//the whole cone used for jet-tracks association is outside of the tracker acceptance
		  fTChfrac[jqi]=-1.;
		}
		// Convert tracks to transient tracks for vertex fitting
		if(AssociatedTTracks.size() > 1){
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
		
		
		// Dump jet properties into tree variables
		fTjpx[jqi]  = jet->px();
		fTjpy[jqi]  = jet->py();
		fTjpz[jqi]  = jet->pz();
		fTjpt[jqi]  = jet->pt();
		fTjeta[jqi] = jet->eta();
		fTjphi[jqi] = jet->phi();
		fTje[jqi]   = jet->energy();
		fTjet[jqi]  = jet->et(); // i know: it's the same as pt, still...
		fTjemfrac[jqi] = jet->emEnergyFraction();
		fTjID_HPD[jqi]  = jetID.fHPD();
		fTjID_RBX[jqi]  = jetID.fRBX();
		fTjID_n90Hits[jqi]  = jetID.n90Hits();
		fTjID_SubDet1[jqi]  = jetID.fSubDetector1();
		fTjID_SubDet2[jqi]  = jetID.fSubDetector2();
		fTjID_SubDet3[jqi]  = jetID.fSubDetector3();
		fTjID_SubDet4[jqi]  = jetID.fSubDetector4();
		fTjID_resEMF[jqi]  = jetID.restrictedEMF();
		fTjID_HCALTow[jqi]  = jetID.nHCALTowers();
		fTjID_ECALTow[jqi]  = jetID.nECALTowers();
		fTjEcorr[jqi]  =  jet->px()/fTUNC_px_match[itagcorr];
	}
	fTnjets = jqi+1;

	// Get and Dump MET Variables:
	fTTrkPtSumx = 0.; fTTrkPtSumy = 0.;
	for(TrackCollection::const_iterator it = tracks->begin(); it != tracks->end() ; ++it){
		fTTrkPtSumx += it->px();
		fTTrkPtSumy += it->py();
	}
	fTTrkPtSum = sqrt(fTTrkPtSumx*fTTrkPtSumx + fTTrkPtSumy*fTTrkPtSumy);

	fTECALEsumx = 0.; fTECALEsumy = 0.; fTECALEsumz = 0.;
	fTHCALEsumx = 0.; fTHCALEsumy = 0.; fTHCALEsumz = 0.;
	for(CaloTowerCollection::const_iterator itow = calotowers->begin(); itow!=calotowers->end(); ++itow){
		double emFrac = itow->emEnergy()/itow->energy();
		double hadFrac = itow->hadEnergy()/itow->energy();
		fTECALEsumx += itow->px()*emFrac;
		fTECALEsumy += itow->py()*emFrac;
		fTECALEsumz += itow->pz()*emFrac;
		fTHCALEsumx += itow->px()*hadFrac;
		fTHCALEsumy += itow->py()*hadFrac;
		fTHCALEsumz += itow->pz()*hadFrac;
	}
	fTECALMET = sqrt(fTECALEsumx*fTECALEsumx + fTECALEsumy*fTECALEsumy);
	fTHCALMET = sqrt(fTHCALEsumx*fTHCALEsumx + fTHCALEsumy*fTHCALEsumy);

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

	if(acceptEvent){
		fTree->Fill();
		fNFillTree++;
	}
}

// Method called once each job just before starting event loop
void NTupleProducer::beginJob(const edm::EventSetup&){
	fNTotEvents = 0;
	fNFillTree  = 0;
	fFirstevent = true;
	fHtrigstat = fTFileService->make<TH1I>("TriggerStats", "TriggerStatistics", 200, 0, 200);
	fTree = fTFileService->make<TTree>("Analysis", "ETHZAnalysisTree");
	// Event information:
	fTree->Branch("Run"            ,&fTrunnumber      ,"Run/I");
	fTree->Branch("Event"          ,&fTeventnumber    ,"Event/I");
	fTree->Branch("LumiSection"    ,&fTlumisection    ,"LumiSection/I");
	fTree->Branch("Weight"         ,&fTweight         ,"Weight/D");
	fTree->Branch("TrigResults"    ,&fTtrigres        ,"TrigResults[200]/I");
	fTree->Branch("PrimVtxx"       ,&fTprimvtxx       ,"PrimVtxx/D");
	fTree->Branch("PrimVtxy"       ,&fTprimvtxy       ,"PrimVtxy/D");
	fTree->Branch("PrimVtxz"       ,&fTprimvtxz       ,"PrimVtxz/D");	
	fTree->Branch("PrimVtxxE"      ,&fTprimvtxxE      ,"PrimVtxxE/D");
	fTree->Branch("PrimVtxyE"      ,&fTprimvtxyE      ,"PrimVtxyE/D");
	fTree->Branch("PrimVtxzE"      ,&fTprimvtxzE      ,"PrimVtxzE/D");	
	fTree->Branch("PrimVtxzNChi2"  ,&fTpvtxznchi2     ,"PrimVtxzNChi2/D");	
	fTree->Branch("Beamspotx"      ,&fTbeamspotx      ,"Beamspotx/D");
	fTree->Branch("Beamspoty"      ,&fTbeamspoty      ,"Beamspoty/D");
	fTree->Branch("Beamspotz"      ,&fTbeamspotz      ,"Beamspotz/D");

	// Muons:
	fTree->Branch("NMus"           ,&fTnmu            ,"NMus/I");
	fTree->Branch("MuPx"           ,&fTmupx           ,"MuPx[NMus]/D");
	fTree->Branch("MuPy"           ,&fTmupy           ,"MuPy[NMus]/D");
	fTree->Branch("MuPz"           ,&fTmupz           ,"MuPz[NMus]/D");
	fTree->Branch("MuPt"           ,&fTmupt           ,"MuPt[NMus]/D");
	fTree->Branch("MuE"            ,&fTmue            ,"MuE[NMus]/D");
	fTree->Branch("MuEt"           ,&fTmuet           ,"MuEt[NMus]/D");
	fTree->Branch("MuEta"          ,&fTmueta          ,"MuEta[NMus]/D");
	fTree->Branch("MuPhi"          ,&fTmuphi          ,"MuPhi[NMus]/D");
	fTree->Branch("MuCharge"       ,&fTmucharge       ,"MuCharge[NMus]/I");
	fTree->Branch("MuPtsum"        ,&fTmuptsum        ,"MuPtsum[NMus]/D");
	fTree->Branch("MuEtsum"        ,&fTmuetsum        ,"MuEtsum[NMus]/D");
	fTree->Branch("MuIso"          ,&fTmuiso          ,"MuIso[NMus]/D");
	fTree->Branch("MuEem"          ,&fTmueecal        ,"MuEem[NMus]/D");
	fTree->Branch("MuEhad"         ,&fTmuehcal        ,"MuEhad[NMus]/D");
	fTree->Branch("MuD0BS"         ,&fTmud0bs         ,"MuD0BS[NMus]/D");
	fTree->Branch("MuD0PV"         ,&fTmud0pv         ,"MuD0PV[NMus]/D");
	fTree->Branch("MuD0E"          ,&fTmud0E          ,"MuD0E[NMus]/D");
	fTree->Branch("MuDzBS"         ,&fTmudzbs         ,"MuDzBS[NMus]/D");
	fTree->Branch("MuDzPV"         ,&fTmudzpv         ,"MuDzPV[NMus]/D");
	fTree->Branch("MuDzE"          ,&fTmudzE          ,"MuDzE[NMus]/D");
	fTree->Branch("MuNChi2"        ,&fTmunchi2        ,"MuNChi2[NMus]/D");
	fTree->Branch("MuNGlHits"      ,&fTmunglhits      ,"MuNGlHits[NMus]/I");
	fTree->Branch("MuNMuHits"      ,&fTmunmuhits      ,"MuNMuHits[NMus]/I");
	fTree->Branch("MuNTkHits"      ,&fTmuntkhits      ,"MuNTkHits[NMus]/I");
	fTree->Branch("MuNMatches"     ,&fTmunmatches     ,"MuNMatches[NMus]/I");
	fTree->Branch("MuNChambers"    ,&fTmunchambers    ,"MuNChambers[NMus]/I");
	fTree->Branch("MuCaloComp"     ,&fTmucalocomp     ,"MuCaloComp[NMus]/D");
	fTree->Branch("MuSegmComp"     ,&fTmusegmcomp     ,"MuSegmComp[NMus]/D");
	fTree->Branch("MuTrackerMu"    ,&fTmutrackermu    ,"MuTrackerMu[NMus]/I");
	fTree->Branch("MuGMPT"         ,&fTmuisGMPT       ,"MuGMPT[NMus]/I");
	fTree->Branch("MuID"           ,&fTmuid           ,"MuID[NMus]/I");
	fTree->Branch("MuMID"          ,&fTmumid          ,"MuMID[NMus]/I");

	// Electrons:
	fTree->Branch("NEles"           ,&fTneles          ,"NEles/I");
	fTree->Branch("ElPx"            ,&fTepx            ,"ElPx[NEles]/D");
	fTree->Branch("ElPy"            ,&fTepy            ,"ElPy[NEles]/D");
	fTree->Branch("ElPz"            ,&fTepz            ,"ElPz[NEles]/D");
	fTree->Branch("ElPt"            ,&fTept            ,"ElPt[NEles]/D");
	fTree->Branch("ElE"             ,&fTee             ,"ElE[NEles]/D");
	fTree->Branch("ElEt"            ,&fTeet            ,"ElEt[NEles]/D");
	fTree->Branch("ElEta"           ,&fTeeta           ,"ElEta[NEles]/D");
	fTree->Branch("ElPhi"           ,&fTephi           ,"ElPhi[NEles]/D");
	fTree->Branch("ElD0BS"          ,&fTed0bs          ,"ElD0BS[NEles]/D");
	fTree->Branch("ElD0PV"          ,&fTed0pv          ,"ElD0PV[NEles]/D");
	fTree->Branch("ElD0E"           ,&fTed0E           ,"ElD0E[NEles]/D");
	fTree->Branch("ElDzBS"          ,&fTedzbs          ,"ElDzBS[NEles]/D");
	fTree->Branch("ElDzPV"          ,&fTedzpv          ,"ElDzPV[NEles]/D");
	fTree->Branch("ElDzE"           ,&fTedzE           ,"ElDzE[NEles]/D");
	fTree->Branch("ElIso"           ,&fTeiso           ,"ElIso[NEles]/D");
	fTree->Branch("ElNChi2"         ,&fTenchi2         ,"ElNChi2[NEles]/D");
	fTree->Branch("ElCharge"        ,&fTecharge        ,"ElCharge[NEles]/I");
	fTree->Branch("ElID"            ,&fTeID            ,"ElID[NEles][4]/I");
	fTree->Branch("ElInGap"         ,&fTeInGap         ,"ElInGap[NEles]/I");
	fTree->Branch("ElEcalDriven"        ,&fTeEcalDriven        ,"ElEcalDriven[NEles]/I");
	fTree->Branch("ElTrackerDriven"     ,&fTeTrackerDriven     ,"ElTrackerDriven[NEles]/I");
	fTree->Branch("ElBasicClustersSize" ,&fTeBasicClustersSize ,"ElBasicClustersSize[NEles]/I");
	fTree->Branch("Elfbrem"             ,&fTefbrem             ,"Elfbrem[NEles]/D");
	fTree->Branch("ElHcalOverEcal"      ,&fTeHcalOverEcal      ,"ElHcalOverEcal[NEles]/D");
	fTree->Branch("ElE5x5"              ,&fTeE5x5              ,"ElE5x5[NEles]/D");
	fTree->Branch("ElE2x5Max"           ,&fTeE2x5Max           ,"ElE2x5Max[NEles]/D");
	fTree->Branch("ElSigmaIetaIeta"     ,&fTeSigmaIetaIeta     ,"ElSigmaIetaIeta[NEles]/D");
	fTree->Branch("ElDeltaPhiSeedClusterAtCalo" ,&fTeDeltaPhiSeedClusterAtCalo ,"ElDeltaPhiSeedClusterAtCalo[NEles]/D");
	fTree->Branch("ElDeltaPhiSuperClusterAtVtx" ,&fTeDeltaPhiSuperClusterAtVtx ,"ElDeltaPhiSuperClusterAtVtx[NEles]/D");
	fTree->Branch("ElESuperClusterOverP"        ,&fTeESuperClusterOverP        ,"ElESuperClusterOverP[NEles]/D");
	fTree->Branch("ElDeltaEtaSeedClusterAtCalo" ,&fTeDeltaEtaSeedClusterAtCalo ,"ElDeltaEtaSeedClusterAtCalo[NEles]/D");

	// Jets:
	fTree->Branch("NJets"          ,&fTnjets        ,"NJets/I");
	fTree->Branch("JPx"            ,&fTjpx          ,"JPx[NJets]/D");
	fTree->Branch("JPy"            ,&fTjpy          ,"JPy[NJets]/D");
	fTree->Branch("JPz"            ,&fTjpz          ,"JPz[NJets]/D");
	fTree->Branch("JPt"            ,&fTjpt          ,"JPt[NJets]/D");
	fTree->Branch("JE"             ,&fTje           ,"JE[NJets]/D");
	fTree->Branch("JEt"            ,&fTjet          ,"JEt[NJets]/D");
	fTree->Branch("JEta"           ,&fTjeta         ,"JEta[NJets]/D");
	fTree->Branch("JPhi"           ,&fTjphi         ,"JPhi[NJets]/D");
	fTree->Branch("JEMfrac"        ,&fTjemfrac      ,"JEMfrac[NJets]/D");
	fTree->Branch("JID_HPD"        ,&fTjID_HPD      ,"JID_HPD[NJets]/D");
	fTree->Branch("JID_RBX"        ,&fTjID_RBX      ,"JID_RBX[NJets]/D");
	fTree->Branch("JID_n90Hits"    ,&fTjID_n90Hits  ,"JID_n90Hits[NJets]/D");
	fTree->Branch("JID_SubDet1"    ,&fTjID_SubDet1  ,"JID_SubDet1[NJets]/D");
	fTree->Branch("JID_SubDet2"    ,&fTjID_SubDet2  ,"JID_SubDet2[NJets]/D");
	fTree->Branch("JID_SubDet3"    ,&fTjID_SubDet3  ,"JID_SubDet3[NJets]/D");
	fTree->Branch("JID_SubDet4"    ,&fTjID_SubDet4  ,"JID_SubDet4[NJets]/D");
	fTree->Branch("JID_resEMF"     ,&fTjID_resEMF   ,"JID_resEMF[NJets]/D");
	fTree->Branch("JID_HCALTow"    ,&fTjID_HCALTow  ,"JID_HCALTow[NJets]/D");
	fTree->Branch("JID_ECALTow"    ,&fTjID_ECALTow  ,"JID_ECALTow[NJets]/D");
	fTree->Branch("JbTagProb"      ,&fTbTagProb     ,"JbTagProb[NJets]/D");
	fTree->Branch("JChfrac"        ,&fTChfrac       ,"JChfrac[NJets]/D");
	fTree->Branch("JEcorr"         ,&fTjEcorr       ,"JEcorr[NJets]/D");
	fTree->Branch("JeMinDR"        ,&fTjeMinDR      ,"JeMinDR[NJets]/D");
	fTree->Branch("JVtxx"          ,&fTjetVtxx      ,"JVtxx[NJets]/D");
	fTree->Branch("JVtxy"          ,&fTjetVtxy      ,"JVtxy[NJets]/D");
	fTree->Branch("JVtxz"          ,&fTjetVtxz      ,"JVtxz[NJets]/D");
	fTree->Branch("JVtxExx"        ,&fTjetVtxExx    ,"JVtxExx[NJets]/D");
	fTree->Branch("JVtxEyx"        ,&fTjetVtxEyx    ,"JVtxEyx[NJets]/D");
	fTree->Branch("JVtxEyy"        ,&fTjetVtxEyy    ,"JVtxEyy[NJets]/D");
	fTree->Branch("JVtxEzy"        ,&fTjetVtxEzy    ,"JVtxEzy[NJets]/D");
	fTree->Branch("JVtxEzz"        ,&fTjetVtxEzz    ,"JVtxEzz[NJets]/D");
	fTree->Branch("JVtxEzx"        ,&fTjetVtxEzx    ,"JVtxEzx[NJets]/D");
	fTree->Branch("JVtxNChi2"      ,&fTjetVtxNChi2  ,"JVtxNChi2[NJets]/D");

	// MET:
	fTree->Branch("TrkPtSumx"      ,&fTTrkPtSumx      ,"TrkPtSumx/D");
	fTree->Branch("TrkPtSumy"      ,&fTTrkPtSumy      ,"TrkPtSumy/D");
	fTree->Branch("TrkPtSum"       ,&fTTrkPtSum       ,"TrkPtSum/D");
	fTree->Branch("ECALEsumx"      ,&fTECALEsumx      ,"ECALEsumx/D");
	fTree->Branch("ECALEsumy"      ,&fTECALEsumy      ,"ECALEsumy/D");
	fTree->Branch("ECALEsumz"      ,&fTECALEsumz      ,"ECALEsumz/D");
	fTree->Branch("ECALMET"        ,&fTECALMET        ,"ECALMET/D");
	fTree->Branch("HCALEsumx"      ,&fTHCALEsumx      ,"HCALEsumx/D");
	fTree->Branch("HCALEsumy"      ,&fTHCALEsumy      ,"HCALEsumy/D");
	fTree->Branch("HCALEsumz"      ,&fTHCALEsumz      ,"HCALEsumz/D");
	fTree->Branch("HCALMET"        ,&fTHCALMET        ,"HCALMET/D");
	fTree->Branch("RawMET"         ,&fTRawMET         ,"RawMET/D");
	fTree->Branch("RawMETpx"       ,&fTRawMETpx       ,"RawMETpx/D");
	fTree->Branch("RawMETpy"       ,&fTRawMETpy       ,"RawMETpy/D");
	fTree->Branch("RawMETphi"      ,&fTRawMETphi      ,"RawMETphi/D");
	fTree->Branch("MuCorrMET"      ,&fTMuCorrMET      ,"MuCorrMET/D");
	fTree->Branch("MuCorrMETpx"    ,&fTMuCorrMETpx    ,"MuCorrMETpx/D");
	fTree->Branch("MuCorrMETpy"    ,&fTMuCorrMETpy    ,"MuCorrMETpy/D");
	fTree->Branch("MuCorrMETphi"   ,&fTMuCorrMETphi   ,"MuCorrMETphi/D");
	fTree->Branch("TCMET"          ,&fTTCMET          ,"TCMET/D");
	fTree->Branch("TCMETpx"        ,&fTTCMETpx        ,"TCMETpx/D");
	fTree->Branch("TCMETpy"        ,&fTTCMETpy        ,"TCMETpy/D");
	fTree->Branch("TCMETphi"       ,&fTTCMETphi       ,"TCMETphi/D");
	fTree->Branch("MuJESCorrMET"   ,&fTMuJESCorrMET   ,"MuJESCorrMET/D");
	fTree->Branch("MuJESCorrMETpx" ,&fTMuJESCorrMETpx ,"MuJESCorrMETpx/D");
	fTree->Branch("MuJESCorrMETpy" ,&fTMuJESCorrMETpy ,"MuJESCorrMETpy/D");
	fTree->Branch("MuJESCorrMETphi",&fTMuJESCorrMETphi,"MuJESCorrMETphi/D");
	fTree->Branch("PFMET"          ,&fTPFMET          ,"PFMET/D");
	fTree->Branch("PFMETpx"        ,&fTPFMETpx        ,"PFMETpx/D");
	fTree->Branch("PFMETpy"        ,&fTPFMETpy        ,"PFMETpy/D");
	fTree->Branch("PFMETphi"       ,&fTPFMETphi       ,"PFMETphi/D");
}

// Method called once each job just after ending the event loop
void NTupleProducer::endJob(){
	cout << " ---------------------------------------------------" << endl;
	cout << " ==> NTupleProducer::endJob() ..." << endl;
	cout << "  Total number of processed Events: " << fNTotEvents << endl;
	cout << "  Number of times Tree was filled:  " << fNFillTree  << endl;
	cout << " ---------------------------------------------------" << endl;	
}

// Method to reset the TTree variables for each event
void NTupleProducer::resetTree(){
	resetInt(fTtrigres, 200);
	fTrunnumber = -999;
	fTeventnumber = -999;
	fTlumisection = -999;
	fTweight = -999.99;
	fTprimvtxx = -999.99;
	fTprimvtxy = -999.99;
	fTprimvtxz = -999.99;
	fTprimvtxxE = -999.99;
	fTprimvtxyE = -999.99;
	fTprimvtxzE = -999.99;
	fTpvtxznchi2 = -999.99;
	fTbeamspotx = -999.99;
	fTbeamspoty = -999.99;
	fTbeamspotz = -999.99;
	fTnmu = 0;
	fTneles = 0;
	fTnjets = 0;
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
	resetDouble(fTeDeltaPhiSuperClusterAtVtx);
	resetDouble(fTeESuperClusterOverP);
	resetDouble(fTeDeltaEtaSeedClusterAtCalo);

	resetDouble(fTjpx);
	resetDouble(fTjpy);
	resetDouble(fTjpz);
	resetDouble(fTje);
	resetDouble(fTjet);
	resetDouble(fTjpt);
	resetDouble(fTjeta);
	resetDouble(fTjphi);
	resetDouble(fTjemfrac);
	resetDouble(fTjID_HPD);
	resetDouble(fTjID_RBX);
	resetDouble(fTjID_n90Hits);
	resetDouble(fTjID_SubDet1);
	resetDouble(fTjID_SubDet2);
	resetDouble(fTjID_SubDet3);
	resetDouble(fTjID_SubDet4);
	resetDouble(fTjID_resEMF);
	resetDouble(fTjID_HCALTow);
	resetDouble(fTjID_ECALTow);
	resetDouble(fTbTagProb);
	resetDouble(fTChfrac);
	resetDouble(fTjeMinDR);
	resetDouble(fTjetVtxx);
	resetDouble(fTjetVtxy);
	resetDouble(fTjetVtxz);
	resetDouble(fTjetVtxExx);
	resetDouble(fTjetVtxEyx);
	resetDouble(fTjetVtxEyy);
	resetDouble(fTjetVtxEzy);
	resetDouble(fTjetVtxEzz);
	resetDouble(fTjetVtxEzx);
	resetDouble(fTjetVtxNChi2);
	resetDouble(fTUNC_px_match,50);
	resetDouble(fTUNC_py_match,50);
	resetDouble(fTUNC_pz_match,50);

	for(size_t i = 0; i < 20; ++i){
		for(size_t j = 0; j < 4; ++j){
			fTeID[i][j] = -999;
		}
	}

	fTTrkPtSumx       = -999.99;
	fTTrkPtSumy       = -999.99;
	fTTrkPtSum        = -999.99;
	fTECALEsumx       = -999.99;
	fTECALEsumy       = -999.99;
	fTECALEsumz       = -999.99;
	fTECALMET         = -999.99;
	fTHCALEsumx       = -999.99;
	fTHCALEsumy       = -999.99;
	fTHCALEsumz       = -999.99;
	fTHCALMET         = -999.99;
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
		double DR = GetDeltaR(eta, mueta, phi, muphi);
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
		double DR = GetDeltaR(mueta, eta, muphi, phi);
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
		double DR = GetDeltaR(eta, mueta, phi, muphi);
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
		double DR = GetDeltaR(mueta, eta, muphi, phi);
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
		double DR = GetDeltaR(eta, eleta, phi, elphi);
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
		double DR = GetDeltaR(eleta, eta, elphi, phi);
		if(DR <= 0.) DR = 0.001;
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
	// 	double eta = itow->eta();
	// 	if(itow->energy()/cosh(eta) < fIso_MuCalSeed) continue;
	// 	double phi = itow->phi();
	// 	double DR = GetDeltaR(eleta, eta, elphi, phi);
	// 	if(DR <= 0.) DR = 0.001;
	// 	if(DR < fIso_MuCalDRin) continue;
	// 	if(DR > fIso_MuCalDRout) continue;
	// 	etsum += itow->energy()/cosh(eta);
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
	int muid(0);
	int mumid(0);
	int match(0);
	vector<int> res;
	double mupt = Mu->globalTrack()->pt();
	double mueta = Mu->globalTrack()->eta();
	double muphi = Mu->globalTrack()->phi();
	int mucharge = Mu->globalTrack()->charge();
	edm::Handle<GenParticleCollection> genparts;
	iEvent.getByLabel(fGenPartTag,genparts);
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
		double dr = GetDeltaR(gpart->eta(), mueta, gpart->phi(), muphi);
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
/*	- This will arrange the vector of muons such that they are ordered in pT   */
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

double NTupleProducer::DeltaPhi(double v1, double v2){
	// Computes the correctly normalized phi difference
	// v1, v2 = phi of object 1 and 2
	double diff = fabs(v2 - v1);
	double corr = 2*acos(-1.) - diff;
	if(diff < acos(-1.)) return diff;
	else return corr;
}

double NTupleProducer::GetDeltaR(double eta1, double eta2, double phi1, double phi2){
	// Computes the DeltaR of two objects from their eta and phi values
	return sqrt( (eta1-eta2)*(eta1-eta2)
		+ DeltaPhi(phi1, phi2)*DeltaPhi(phi1, phi2) );
}

vector<const reco::Track*> NTupleProducer::FindAssociatedTracks(const reco::Jet *jet, const reco::TrackCollection *tracks){
//returns a list of tracks associated to a given jet
	vector<const reco::Track*> AssociatedTracks;
	double jeteta=jet->eta();
	double jetphi=jet->phi();
	for(TrackCollection::const_iterator itk = tracks->begin(); itk!=tracks->end(); ++itk){
		if(sqrt(pow(itk->eta()-jeteta,2) + pow(itk->phi()-jetphi,2))<0.5) {
			const Track *t = new Track(*itk);
			AssociatedTracks.push_back(t);
		}
	}
	return AssociatedTracks;
}

//define this as a plug-in
DEFINE_FWK_MODULE(NTupleProducer);
