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
// $Id: NTupleProducer.h,v 1.1 2009/09/16 14:57:24 stiegerb Exp $
//
//

// system include files
#include <memory>
#include <vector>
#include <iostream>
#include <cmath>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonEnergy.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/JetReco/interface/JetCollection.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/CaloTowers/interface/CaloTower.h"
#include "DataFormats/CaloTowers/interface/CaloTowerCollection.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"

#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/METReco/interface/METCollection.h"
#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/METReco/interface/CaloMETCollection.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "PhysicsTools/UtilAlgos/interface/TFileService.h"

#include "RecoEgamma/Examples/plugins/ElectronIDAnalyzer.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "DataFormats/RecoCandidate/interface/IsoDeposit.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Framework/interface/TriggerNamesService.h"

// #include "SSDiLeptonAnalysis/NTupleProducer/interface/MyMuon.h"

#include "TH1.h"
#include "TH2.h"
#include "TLorentzVector.h"
#include "TTree.h"

class NTupleProducer : public edm::EDAnalyzer {
public:
	explicit NTupleProducer(const edm::ParameterSet&);
	~NTupleProducer();
	vector<double> calcMuIso(const reco::Muon *Mu, const edm::Event& iEvent);
	vector<double> calcMuIso2(const reco::Muon *Mu, const edm::Event& iEvent);
	vector<double> calcMuIso3(const reco::Muon *Mu, const edm::Event& iEvent, edm::Ref<reco::MuonCollection> muonRef);
	vector<double> calcElIso(const reco::GsfElectron *El, const edm::Event& iEvent);
	vector<int> matchMuCand(const reco::Muon *Mu, const edm::Event& iEvent);
	double DeltaPhi(double, double);
	double GetDeltaR(double, double, double, double);
	vector<const reco::Muon*> sortMus(vector<const reco::Muon*>);
	void switchDouble(double &, double &);
	void switchInt(int &, int &);
	void resetDouble(double *v, unsigned int size = 20);
	void resetInt(int *v, unsigned int size = 20);
	void resetTree();
private:
	virtual void beginJob(const edm::EventSetup&) ;
	virtual void analyze(const edm::Event&, const edm::EventSetup&);
	virtual void endJob() ;
		
		// ----------member data ---------------------------
	int fNTotEvents;
	int fNFillTree;

	edm::InputTag fMuonTag;
	edm::InputTag fElectronTag;
	edm::InputTag fEleIsoTkTag;
	edm::InputTag fEleIsoECTag;
	edm::InputTag fEleIsoHCTag;
	edm::InputTag fEleIsoDepTkTag;
	edm::InputTag fEleIsoDepECTag;
	edm::InputTag fEleIsoDepHCTag;
	edm::InputTag fMuIsoDepTkTag;
	edm::InputTag fMuIsoDepECTag;
	edm::InputTag fMuIsoDepHCTag;
	edm::InputTag fSCTag;
	edm::InputTag fJetTag;
	edm::InputTag fMET1Tag;
	edm::InputTag fMET2Tag;
	edm::InputTag fMET3Tag;
	edm::InputTag fVertexTag;
	edm::InputTag fTrackTag;
	edm::InputTag fCalTowTag;
	edm::InputTag fGenPartTag;
	edm::InputTag fTriggerTag;
	edm::Service<TFileService> fTFileService;
	double fMinmupt;
	double fMaxmueta;
	double fMinelpt;
	double fMaxeleta;
	double fMaxeliso;
	double fMaxeld0;

	double fMinjpt;
	double fMaxjeta;
	double fMinjemfrac;
	
	double fIso_MuTkDRin;
	double fIso_MuTkDRout;
	double fIso_MuTkSeed;
	double fIso_MuCalDRin;
	double fIso_MuCalDRout;
	double fIso_MuCalSeed;
	
	////////////////////////////////////////////////////////
	// Tree:
	TTree *fTree;
	
	// Trigger
	int fTtrigres[200];
	TH1I *fHtrigstat; // Added to keep track of trigger names
	bool fFirstevent;
	// Muons:
	unsigned int fTnmu;
	double fTmupx[20];
	double fTmupy[20];
	double fTmupz[20];
	double fTmue[20];
	double fTmuet[20];
	double fTmupt[20];
	double fTmueta[20];
	double fTmuphi[20];
	int fTmucharge[20];

	// - Isolation Variables
	double fTmuetsum[20];
	double fTmuptsum[20];
	double fTmuiso[20];
	double fTmueecal[20];
	double fTmuehcal[20];

	// - Impact Parameters
	double fTmud0[20];
	double fTmud0E[20];
	double fTmud0Sig[20];
	double fTmudz[20];
	double fTmudzE[20];
	double fTmudzSig[20];

	// - MuID Variables
	double fTmunchi2[20];
	int fTmunglhits[20];
	int fTmunmuhits[20];
	int fTmuntkhits[20];
	int fTmunmatches[20];
	int fTmunchambers[20];
	double fTmucalocomp[20];
	double fTmusegmcomp[20];
	int fTmutrackermu[20];
	int fTmuisGMPT[20];

	// - Gen Info:
	int fTmuid[20];
	int fTmumid[20];

	// MET:
	double fTRawMET;
	double fTRawMETpx;
	double fTRawMETpy;
	double fTRawMETphi;
	double fTMuCorrMET;
	double fTMuCorrMETpx;
	double fTMuCorrMETpy;
	double fTMuCorrMETphi;
	double fTTCMET;
	double fTTCMETpx;
	double fTTCMETpy;
	double fTTCMETphi;

	// Electrons:
	int fTneles;
	double fTepx[20];
	double fTepy[20];
	double fTepz[20];
	double fTept[20];
	double fTee[20];
	double fTeet[20];
	double fTeeta[20];
	double fTephi[20];
	double fTed0[20];
	double fTed0E[20];
	double fTedz[20];
	double fTedzE[20];
	double fTeiso[20];
	double fTenchi2[20];
	int fTeID[20][4];    // eID flags: 0->Tight, 1->Loose, 2->RobustTight, 3->RobustLoose
	int fTecharge[20];

	// Jets:
	int fTnjets;
	double fTjpx[20];
	double fTjpy[20];
	double fTjpz[20];
	double fTje[20];
	double fTjet[20];
	double fTjpt[20];
	double fTjeta[20];
	double fTjphi[20];
	double fTjemfrac[20];
	////////////////////////////////////////////////////////
};
