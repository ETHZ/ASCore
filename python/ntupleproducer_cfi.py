import FWCore.ParameterSet.Config as cms


## !! mu, el photon eta changed from 2.4 to 2.5 !!
analyze = cms.EDAnalyzer('NTupleProducer',
	# Main settings
	isRealData = cms.untracked.bool(False),
	isPat      = cms.untracked.bool(False), 
	# Collections
	tag_muons      = cms.untracked.InputTag('muons'),
	tag_muisodeptk = cms.untracked.InputTag('muIsoDepositTk'),
	tag_muisodepec = cms.untracked.InputTag("muIsoDepositCalByAssociatorTowers","ecal"),
	tag_muisodephc = cms.untracked.InputTag("muIsoDepositCalByAssociatorTowers","hcal"),
	tag_electrons  = cms.untracked.InputTag('gsfElectrons'),	
	tag_elidWP     = cms.untracked.string('simpleEleId90relIso'),	
	tag_jets       = cms.untracked.InputTag('ak5CaloJets'),
	jetCorrs       = cms.untracked.string('ak5CaloL2L3'),
	tag_btag1       = cms.untracked.InputTag('trackCountingHighEffBJetTags'), 
	tag_btag2       = cms.untracked.InputTag('trackCountingHighPurBJetTags'), 
	tag_btag3       = cms.untracked.InputTag('simpleSecondaryVertexHighEffBJetTags'), 
	tag_btag4       = cms.untracked.InputTag('simpleSecondaryVertexHighPurBJetTags'), 		
	tag_jetTracks  = cms.untracked.InputTag('ak5JetTracksAssociatorAtVertex'),
	tag_jetID      = cms.untracked.InputTag('ak5JetID'),
	#trackCountingHighPurBJetTags #jetProbabilityBJetTags
	tag_met1       = cms.untracked.InputTag('met'),
	tag_met2       = cms.untracked.InputTag('corMetGlobalMuons'),
	tag_met3       = cms.untracked.InputTag('tcMet'),
	tag_met4       = cms.untracked.InputTag('pfMet'),
	tag_met5       = cms.untracked.InputTag('metMuonJESCorAK5'),
	tag_vertex     = cms.untracked.InputTag('offlinePrimaryVertices'),
	# tag_vertex   = cms.untracked.InputTag('offlinePrimaryVertices','','REVERTEX'),
	tag_tracks     = cms.untracked.InputTag('generalTracks'),
	# tag_tracks   = cms.untracked.InputTag('generalTracks','','RETRACK'),
	tag_photons    = cms.untracked.InputTag('photons'),
	tag_caltow     = cms.untracked.InputTag('towerMaker'),
	tag_EBrechits   = cms.untracked.InputTag('ecalRecHit:EcalRecHitsEB'),
	tag_EErechits   = cms.untracked.InputTag('ecalRecHit:EcalRecHitsEE'),
	tag_genpart    = cms.untracked.InputTag('genParticles'),
	tag_genjets    = cms.untracked.InputTag('ak5GenJets'),
	tag_l1trig     = cms.untracked.InputTag("gtDigis"),
	tag_hlttrig    = cms.untracked.InputTag("TriggerResults","","HLT"),
	tag_hcalnoise  = cms.untracked.InputTag('HBHENoiseFilterResultProducer','HBHENoiseFilterResult'),
        
	# Event Selection Criteria
	# Muons:
	sel_minmupt     = cms.double(5.0),
	sel_maxmueta    = cms.double(2.5),
	# Electrons: 
	sel_minelpt     = cms.double(5.0),
	sel_maxeleta    = cms.double(2.5),
	# Jets:
	sel_mincorjpt   = cms.double(20.0),
	sel_minrawjpt   = cms.double(0.0),
	sel_maxjeta     = cms.double(10.0),
	sel_minjemfrac  = cms.double(0.0),
	# Tracks:
	sel_mintrkpt    = cms.double(1.0),
	sel_maxtrketa   = cms.double(10.0),
	sel_maxtrknchi2 = cms.double(1e15),
	sel_mintrknhits = cms.int32(0),
	# Photons
	sel_minphopt    = cms.double(5.0),
	sel_maxphoeta   = cms.double(2.5),
	# GenLeptons
	sel_mingenleptpt    = cms.double(5.0),
	sel_maxgenlepteta   = cms.double(10),	

	# Additional jet collections
	jets = cms.VPSet()

)

