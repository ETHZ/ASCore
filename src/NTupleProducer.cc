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
// $Id: NTupleProducer.cc,v 1.2 2009/09/18 12:33:08 stiegerb Exp $
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
	fMET1Tag        = iConfig.getUntrackedParameter<edm::InputTag>("tag_met1");
	fMET2Tag        = iConfig.getUntrackedParameter<edm::InputTag>("tag_met2");
	fMET3Tag        = iConfig.getUntrackedParameter<edm::InputTag>("tag_met3");
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

	// MET
	Handle<CaloMETCollection> calomet;
	iEvent.getByLabel(fMET1Tag, calomet);

	Handle<CaloMETCollection> corrmumet;
	iEvent.getByLabel(fMET2Tag, corrmumet);

	Handle<METCollection> tcmet;
	iEvent.getByLabel(fMET3Tag, tcmet);

	// Get beamspot for d0 determination
	BeamSpot beamSpot;
	Handle<BeamSpot> beamSpotHandle;
	iEvent.getByLabel("offlineBeamSpot", beamSpotHandle);
	beamSpot = *beamSpotHandle;

	// Get Muon IsoDeposits
	// ECAL:
	edm::Handle<edm::ValueMap<reco::IsoDeposit> > IsoDepECValueMap;
	iEvent.getByLabel(fMuIsoDepECTag, IsoDepECValueMap);
	const edm::ValueMap<reco::IsoDeposit> &ECDepMap = *IsoDepECValueMap.product();
	// HCAL:
	edm::Handle<edm::ValueMap<reco::IsoDeposit> > IsoDepHCValueMap;
	iEvent.getByLabel(fMuIsoDepHCTag, IsoDepHCValueMap);
	const edm::ValueMap<reco::IsoDeposit> &HCDepMap = *IsoDepHCValueMap.product();

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
			fHtrigstat->GetXaxis()->SetBinLabel(i+1,TString(triggernames[i]));
			// cout << "i " << i << " : " << tr[i].accept() << " : " << triggernames[i] << endl;
		}
		// Add a warning about the shift between trigger bits and bin numbers:
		fHtrigstat->GetXaxis()->SetBinLabel(tr.size()+2, "Bin#=Bit#+1");
	}
	for( unsigned int i = 0; i < tr.size(); i++ ){
		if(tr[i].accept()){
			fTtrigres[i] = 1;
			fHtrigstat->Fill(i);
		}
		if(!tr[i].accept())fTtrigres[i] = 0;
	}

////////////////////////////////////////////////////////////////////////////////
// Event Selection /////////////////////////////////////////////////////////////
	bool acceptEvent = true;

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

		fTmud0[mqi]    = -1.0*Mit->innerTrack()->dxy(beamSpot.position());
		fTmud0E[mqi]   = Mit->globalTrack()->dxyError();
		fTmud0Sig[mqi] = fabs(fTmud0[mqi])/fTmud0E[mqi];
		fTmudz[mqi]    = Mit->globalTrack()->dz(beamSpot.position());
		fTmudzE[mqi]   = Mit->globalTrack()->dzError();
		fTmudzSig[mqi] = fabs(fTmudz[mqi])/fTmudzE[mqi];

		fTmunchi2[mqi]     = Mit->globalTrack()->normalizedChi2();
		fTmunglhits[mqi]   = Mit->globalTrack()->hitPattern().numberOfValidHits();
		fTmunmuhits[mqi]   = Mit->outerTrack()->hitPattern().numberOfValidHits(); // always 0
		fTmuntkhits[mqi]   = Mit->innerTrack()->hitPattern().numberOfValidHits();
		fTmunmatches[mqi]  = Mit->numberOfMatches();
		fTmunchambers[mqi] = Mit->numberOfChambers();

		fTmucalocomp[mqi] = Mit->caloCompatibility();
		// fTmusegmcomp[mqi] = Mit->segmentCompatibility();
		if(Mit->isTrackerMuon()) fTmutrackermu[mqi] = 1;
		else fTmutrackermu[mqi] = 0;
		if(muon::isGoodMuon(*cm, muon::GlobalMuonPromptTight)) fTmuisGMPT[mqi] = 1;
		else fTmuisGMPT[mqi] = 0;

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
			fTed0[eqi]     = El->gsfTrack()->dxy(beamSpot.position());
			fTed0E[eqi]    = El->gsfTrack()->dxyError();
			fTedz[eqi]     = El->gsfTrack()->dz();
			fTedzE[eqi]    = El->gsfTrack()->dzError();
			fTenchi2[eqi]  = El->gsfTrack()->normalizedChi2();
			fTeiso[eqi]    = eliso;
			fTecharge[eqi] = El->charge();
			fTeInGap[eqi]  = isGap ? 1:0;  // my first contribution in beni's code (isBug? 1:1)   //*** kostas
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
			
			if(eIDmapT[electronRef])  fTeID[eqi][0] = 1;
			else fTeID[eqi][0]                      = 0;
			if(eIDmapL[electronRef])  fTeID[eqi][1] = 1;
			else fTeID[eqi][1]                      = 0;
			if(eIDmapRT[electronRef]) fTeID[eqi][2] = 1;
			else fTeID[eqi][2]                      = 0;
			if(eIDmapRL[electronRef]) fTeID[eqi][3] = 1;
			else fTeID[eqi][3]                      = 0;
		}
	}
	fTneles = eqi+1;

	// Apply L2 and L3 JetCorrections
	vector<const reco::CaloJet*> CorrJets; // Contains all (corrected) jets
	CaloJetCollection::const_iterator Jit;
	for(Jit = jets->begin(); Jit != jets->end(); ++Jit){// Loop over uncorr. jets
		CaloJet *j = Jit->clone();
		j->scaleEnergy(L2JetCorrector->correction(j->p4()));
		j->scaleEnergy(L3JetCorrector->correction(j->p4()));
		const CaloJet *cj = j->clone();
		CorrJets.push_back(cj);
	}

	// Determine qualified jets
	int jqi(-1); // counts # of qualified jets
	for(size_t i = 0; i < CorrJets.size(); ++i){// Loop over corr. jets
		// Check if maximum number of jets is exceeded already:
		if(jqi > 19){
			cout << "NTupleSkeleton::analyze() ==> Warning: Maximum number of jets exceeded..." << endl;
			break;
		}
		CaloJet *jet = CorrJets[i]->clone();
		// Jet preselection
		if(jet->pt() < fMinjpt) continue;
		if(TMath::Abs(jet->eta()) > fMaxjeta) continue;
		if(jet->emEnergyFraction() < fMinjemfrac) continue;
		// Ignore jets within 0.4 DR of a selected electron
		bool eleprox = false;
		for(int j = 0; j < fTneles; j++){
			double ejDR = GetDeltaR(jet->eta(), fTeeta[j], jet->phi(), fTephi[j]);
			if(ejDR < 0.4) eleprox = true;
		}
		if(eleprox) continue;
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
		fTjemfrac[jqi] = jet->emEnergyFraction();
	}
	fTnjets = jqi+1;

	// Get and Dump MET Variables:
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
	fHtrigstat = fTFileService->make<TH1I>("TriggerStats", "TriggerStatistics", 200,0,200);
	fTree = fTFileService->make<TTree>("Analysis", "MuonAnalysisTree");
	fTree->Branch("TrigResults"    ,&fTtrigres        ,"TrigResults[200]/I");
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
	fTree->Branch("MuD0"           ,&fTmud0           ,"MuD0[NMus]/D");
	fTree->Branch("MuD0E"          ,&fTmud0E          ,"MuD0E[NMus]/D");
	fTree->Branch("MuD0Sig"        ,&fTmud0Sig        ,"MuD0Sig[NMus]/D");
	fTree->Branch("MuDz"           ,&fTmudz           ,"MuDz[NMus]/D");
	fTree->Branch("MuDzE"          ,&fTmudzE          ,"MuDzE[NMus]/D");
	fTree->Branch("MuDzSig"        ,&fTmudzSig        ,"MuDzSig[NMus]/D");
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
	fTree->Branch("ElD0"            ,&fTed0            ,"ElD0[NEles]/D");
	fTree->Branch("ElD0E"           ,&fTed0E           ,"ElD0E[NEles]/D");
	fTree->Branch("ElDz"            ,&fTedz            ,"ElDz[NEles]/D");
	fTree->Branch("ElDzE"           ,&fTedzE           ,"ElDzE[NEles]/D");
	fTree->Branch("ElIso"           ,&fTeiso           ,"ElIso[NEles]/D");
	fTree->Branch("ElNChi2"         ,&fTenchi2         ,"ElNChi2[NEles]/D");
	fTree->Branch("ElCharge"        ,&fTecharge        ,"ElCharge[NEles]/I");
	fTree->Branch("ElID"            ,&fTeID            ,"ElID[NEles][4]/I");
	fTree->Branch("ElInGap"         ,&fTeInGap         ,"ElInGap[NEles]/I");
	fTree->Branch("ElEcalDriven"         ,&fTeEcalDriven         ,"ElEcalDriven[NEles]/I");
	fTree->Branch("ElTrackerDriven"         ,&fTeTrackerDriven         ,"ElTrackerDriven[NEles]/I");
	fTree->Branch("ElBasicClustersSize"         ,&fTeBasicClustersSize         ,"ElBasicClustersSize[NEles]/I");
	fTree->Branch("Elfbrem"            ,&fTefbrem            ,"Elfbrem[NEles]/D");
	fTree->Branch("ElHcalOverEcal"            ,&fTeHcalOverEcal            ,"ElHcalOverEcal[NEles]/D");
	fTree->Branch("ElE5x5"            ,&fTeE5x5            ,"ElE5x5[NEles]/D");
	fTree->Branch("ElE2x5Max"            ,&fTeE2x5Max            ,"ElE2x5Max[NEles]/D");
	fTree->Branch("ElSigmaIetaIeta"            ,&fTeSigmaIetaIeta            ,"ElSigmaIetaIeta[NEles]/D");
	fTree->Branch("ElDeltaPhiSeedClusterAtCalo"            ,&fTeDeltaPhiSeedClusterAtCalo            ,"ElDeltaPhiSeedClusterAtCalo[NEles]/D");
	fTree->Branch("ElDeltaPhiSuperClusterAtVtx"            ,&fTeDeltaPhiSuperClusterAtVtx            ,"ElDeltaPhiSuperClusterAtVtx[NEles]/D");
	fTree->Branch("ElESuperClusterOverP"            ,&fTeESuperClusterOverP            ,"ElESuperClusterOverP[NEles]/D");
	fTree->Branch("ElDeltaEtaSeedClusterAtCalo"            ,&fTeDeltaEtaSeedClusterAtCalo            ,"ElDeltaEtaSeedClusterAtCalo[NEles]/D");

	// Jets:
	fTree->Branch("NJets"          ,&fTnjets          ,"NJets/I");
	fTree->Branch("JPx"            ,&fTjpx            ,"JPx[NJets]/D");
	fTree->Branch("JPy"            ,&fTjpy            ,"JPy[NJets]/D");
	fTree->Branch("JPz"            ,&fTjpz            ,"JPz[NJets]/D");
	fTree->Branch("JPt"            ,&fTjpt            ,"JPt[NJets]/D");
	fTree->Branch("JE"             ,&fTje             ,"JE[NJets]/D");
	fTree->Branch("JEt"            ,&fTjet            ,"JEt[NJets]/D");
	fTree->Branch("JEta"           ,&fTjeta           ,"JEta[NJets]/D");
	fTree->Branch("JPhi"           ,&fTjphi           ,"JPhi[NJets]/D");
	fTree->Branch("JEMfrac"        ,&fTjemfrac        ,"JEMfrac[NJets]/D");
	
	// MET:
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
	resetDouble(fTmud0);
	resetDouble(fTmud0E);
	resetDouble(fTmud0Sig);
	resetDouble(fTmudz);
	resetDouble(fTmudzE);
	resetDouble(fTmudzSig);
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
	resetDouble(fTed0);
	resetDouble(fTed0E);
	resetDouble(fTedz);
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

	for(size_t i = 0; i < 20; ++i){
		for(size_t j = 0; j < 4; ++j){
			fTeID[i][j] = -999;
		}
	}

	fTRawMET         = -999.99;
	fTRawMETpx       = -999.99;
	fTRawMETpy       = -999.99;
	fTRawMETphi      = -999.99;
	fTMuCorrMET      = -999.99;
	fTMuCorrMETpx    = -999.99;
	fTMuCorrMETpy    = -999.99;
	fTMuCorrMETphi   = -999.99;
	fTTCMET          = -999.99;
	fTTCMETpx        = -999.99;
	fTTCMETpy        = -999.99;
	fTTCMETphi       = -999.99;

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

//define this as a plug-in
DEFINE_FWK_MODULE(NTupleProducer);
