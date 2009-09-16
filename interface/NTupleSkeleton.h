// -*- C++ -*-
//
// Package:    NTupleSkeleton
// Class:      NTupleSkeleton
//
/* class NTupleSkeleton
   NTupleSkeleton.h
   AnalysisExamples/NTupleSkeleton/src/NTupleSkeleton.h
   Description: Template to produce NTuples for ETH SUSY Analysis

   Implementation:
	
*/
//
// Original Author:  Benjamin Stieger
//         Created:  Wed Sep  2 16:43:05 CET 2009
// $Id$
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
#include "DataFormats/CaloTowers/interface/CaloTower.h"
#include "DataFormats/CaloTowers/interface/CaloTowerCollection.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"

#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/METReco/interface/METCollection.h"
#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/METReco/interface/CaloMETCollection.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "PhysicsTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/Common/interface/ValueMap.h"

#include "TH1.h"
#include "TH2.h"
#include "TLorentzVector.h"
#include "TTree.h"

class NTupleSkeleton : public edm::EDAnalyzer {
public:
	explicit NTupleSkeleton(const edm::ParameterSet&);
	~NTupleSkeleton();
	void resetDouble(double (&v)[20], unsigned int size = 20);
	void resetInt(int (&v)[20], unsigned int size = 20);
	void resetTree();
private:
	virtual void beginJob(const edm::EventSetup&) ;
	virtual void analyze(const edm::Event&, const edm::EventSetup&);
	virtual void endJob() ;
		
		// ----------member data ---------------------------
	int fNFillTree;

	edm::InputTag fMuonTag;
	edm::InputTag fElectronTag;
	edm::InputTag fJetTag;
	edm::InputTag fMETTag;
	edm::Service<TFileService> fTFileService;
	double fMinmupt;
	double fMaxmueta;
	double fMinelpt;
	double fMaxeleta;
	double fMinjpt;
	double fMaxjeta;
	
	////////////////////////////////////////////////////////
	// Tree:
	TTree *fTree;

	// Muons:
	unsigned int fTnmu;
	double fTmupt[20];
	double fTmud0[20];
	int fTmuisGMPT[20];

	// Electrons:
	int fTneles;
	double fTept[20];

	// Jets:
	int fTnjets;
	double fTjpt[20];

	// MET:
	double fTRawMET;
	////////////////////////////////////////////////////////
};
