// -*- C++ -*-
//
// Package:    NTupleSkeleton
// Class:      NTupleSkeleton
//
/* class NTupleSkeleton
   NTupleSkeleton.cc
   AnalysisExamples/NTupleSkeleton/src/NTupleSkeleton.cc
   Description: Template to produce NTuples for ETH SUSY Analysis

   Implementation: None
	
*/
//
// Original Author:  Benjamin Stieger
//         Created:  Wed Sep  2 16:43:05 CET 2009
// $Id: NTupleSkeleton.cc,v 1.1 2009/09/16 14:57:25 stiegerb Exp $
//
//

#include "DiLeptonAnalysis/NTupleProducer/interface/NTupleSkeleton.h"


NTupleSkeleton::NTupleSkeleton(const edm::ParameterSet& iConfig){
	// InputTags
	fMuonTag     = iConfig.getUntrackedParameter<edm::InputTag>("tag_muons");
	fElectronTag = iConfig.getUntrackedParameter<edm::InputTag>("tag_electrons");
	fJetTag      = iConfig.getUntrackedParameter<edm::InputTag>("tag_jets");
	fMETTag      = iConfig.getUntrackedParameter<edm::InputTag>("tag_met");

	// Event Selection
	fMinmupt     = iConfig.getUntrackedParameter<double>("sel_minmupt");
	fMaxmueta    = iConfig.getUntrackedParameter<double>("sel_maxmueta");

	fMinelpt     = iConfig.getUntrackedParameter<double>("sel_minelpt");
	fMaxeleta    = iConfig.getUntrackedParameter<double>("sel_maxeleta");

	fMinjpt      = iConfig.getUntrackedParameter<double>("sel_minjpt");
	fMaxjeta     = iConfig.getUntrackedParameter<double>("sel_maxjeta");

	cout << " --------------------------------" << endl;
	cout << " ==> NTupleSkeleton Constructor ..." << endl;
	cout << endl;
	cout << "  Input Tags:" << endl;
	cout << "    fMuonTag        = " << fMuonTag.label()        << endl;
	cout << "    fElectronTag    = " << fElectronTag.label()    << endl;
	cout << "    fJetTag         = " << fJetTag.label()         << endl;
	cout << "    fMETTag         = " << fMETTag.label()        << endl;
	cout << endl;
	cout << "  Event Selection Parameters:" << endl;
	cout << "    fMinmupt        = " << fMinmupt    << endl;
	cout << "    fMaxmueta       = " << fMaxmueta   << endl;
	cout << "    fMinelpt        = " << fMinelpt    << endl;
	cout << "    fMaxeleta       = " << fMaxeleta   << endl;
	cout << "    fMinjpt         = " << fMinjpt     << endl;
	cout << "    fMaxjeta        = " << fMaxjeta    << endl;
	cout << endl;
	cout << " --------------------------------" << endl;
}


NTupleSkeleton::~NTupleSkeleton(){
}

// Method called once for each event
void NTupleSkeleton::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){
	using namespace edm;
	using namespace std;
	using namespace reco;
	using reco::MuonCollection;
	using reco::CaloJetCollection;

	// Reset all the tree variables
	resetTree();

	// Get the collections:
	Handle<MuonCollection> muons;
	iEvent.getByLabel(fMuonTag, muons); // 'muons'

	Handle<GsfElectronCollection> electrons;
	iEvent.getByLabel(fElectronTag, electrons); // 'gsfElectrons'

	// Jets and Jet Correctors
	Handle<CaloJetCollection> jets;
	iEvent.getByLabel(fJetTag,jets); // 'sisCone5CaloJets'
	const JetCorrector* L2JetCorrector = JetCorrector::getJetCorrector ("L2RelativeJetCorrectorSC5Calo",iSetup);
	const JetCorrector* L3JetCorrector = JetCorrector::getJetCorrector ("L3AbsoluteJetCorrectorSC5Calo",iSetup);

	// MET
	Handle<CaloMETCollection> calomet;
	iEvent.getByLabel(fMETTag, calomet);

	// Get beamspot for d0 determination
	BeamSpot beamSpot;
	Handle<BeamSpot> beamSpotHandle;
	iEvent.getByLabel("offlineBeamSpot", beamSpotHandle);
	beamSpot = *beamSpotHandle;

	// Muon loop ////////////////////////////////////////////////////////////////
	vector<const reco::Muon*> Mus;
	int mqi(-1); // counts # of qualified muons
	for(MuonCollection::const_iterator Mit = muons->begin(); Mit != muons->end(); ++Mit){
		// Check if maximum number of muons is exceeded already:
		if(mqi > 19){
			cout << "NTupleSkeleton::analyze() ==> Warning: Maximum number of muons exceeded..." << endl;
			break;
		}
		// Muon preselection:
		if(!(Mit->isGlobalMuon())) continue;
		if(Mit->globalTrack()->pt() < fMinmupt) continue;
		if(TMath::Abs(Mit->globalTrack()->eta()) > fMaxmueta) continue;
		mqi++;
		const Muon *cm = Mit->clone();

		// Dump muon properties in tree variables
		fTmupt[mqi] = Mit->globalTrack()->pt();
		fTmud0[mqi] = -1.0*Mit->innerTrack()->dxy(beamSpot.position());
		if(muon::isGoodMuon(*cm, muon::GlobalMuonPromptTight)) fTmuisGMPT[mqi] = 1;
		else fTmuisGMPT[mqi] = 0;
	}
	fTnmu = mqi+1;

	// Electron loop ////////////////////////////////////////////////////////////
	int eqi(-1); // counts # of qualified electrons
	for(GsfElectronCollection::const_iterator El = electrons->begin(); El!=electrons->end(); ++El){
		// Check if maximum number of electrons is exceeded already:
		if(eqi > 19){
			cout << "NTupleSkeleton::analyze() ==> Warning: Maximum number of electrons exceeded..." << endl;
			break;
		}
		// Electron preselection:			
		if(El->pt() < fMinelpt) continue;
		if(fabs(El->eta()) > fMaxeleta) continue;

		// Dump electron properties in tree variables
		eqi++;
		fTept[eqi] = El->pt();
	}
	fTneles = eqi+1;

	// Jet loop /////////////////////////////////////////////////////////////////
	// Apply L2 and L3 JetCorrections
	vector<const reco::CaloJet*> CorrJets; // Contains all (corrected) jets
	for(CaloJetCollection::const_iterator Jit = jets->begin(); Jit != jets->end(); ++Jit){// Loop over uncorr. jets
		CaloJet *j = Jit->clone();
		j->scaleEnergy(L2JetCorrector->correction(j->p4()));
		j->scaleEnergy(L3JetCorrector->correction(j->p4()));
		const CaloJet *cj = j->clone();
		CorrJets.push_back(cj);
	}

	// Loop over corrected jets
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

		jqi++;
		
		// Dump jet properties into tree variables
		fTjpt[jqi]  = jet->pt();
	}
	fTnjets = jqi+1;

	// Dump MET
	fTRawMET = (calomet->at(0)).pt();

	// Fill TTree
	fTree->Fill();
	fNFillTree++;
}

// Method called once each job just before starting event loop
void NTupleSkeleton::beginJob(const edm::EventSetup&){
	fNFillTree = 0;
	fTree = fTFileService->make<TTree>("Analysis", "MuonAnalysisTree");
	fTree->Branch("NMus"           ,&fTnmu            ,"NMus/I");
	fTree->Branch("MuPt"           ,&fTmupt           ,"MuPt[NMus]/D");
	fTree->Branch("MuD0"           ,&fTmud0           ,"MuD0[NMus]/D");
	fTree->Branch("MuGMPT"         ,&fTmuisGMPT       ,"MuGMPT[NMus]/I");
	fTree->Branch("NEles"          ,&fTneles          ,"NEles/I");
	fTree->Branch("ElPt"           ,&fTept            ,"ElPt[NEles]/D");
	fTree->Branch("NJets"          ,&fTnjets          ,"NJets/I");
	fTree->Branch("JPt"            ,&fTjpt            ,"JPt[NJets]/D");
	fTree->Branch("RawMET"         ,&fTRawMET         ,"RawMET/D");
}

// Method called once each job just after ending the event loop
void NTupleSkeleton::endJob(){
	cout << " ---------------------------------------------" << endl;
	cout << " ==> NTupleSkeleton::endJob() ..." << endl;
	cout << "  Number of times Tree was filled: " << fNFillTree  << endl;
	cout << " ---------------------------------------------" << endl;	
}

// Method to reset the TTree variables for each event
void NTupleSkeleton::resetTree(){
	fTnmu   = 0;
	fTneles = 0;
	fTnjets = 0;
	fTRawMET = -999.99;
	resetDouble(fTmupt);
	resetDouble(fTmud0);
	resetInt(fTmuisGMPT);
	resetDouble(fTept);
	resetDouble(fTjpt);
}

void NTupleSkeleton::resetDouble(double (&v)[20], unsigned int size){
	for(size_t i = 0; i < size; ++i){
		v[i] = -999.99;
	}
}

void NTupleSkeleton::resetInt(int (&v)[20], unsigned int size){
	for(size_t i = 0; i < size; ++i){
		v[i] = -999;
	}
}

//define this as a plug-in
DEFINE_FWK_MODULE(NTupleSkeleton);
