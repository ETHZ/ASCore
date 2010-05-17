import FWCore.ParameterSet.Config as cms

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
	tag_jets       = cms.untracked.InputTag('antikt5CaloJets'),
	jetCorrs       = cms.untracked.string('L2L3JetCorrectorAK5Calo'),
	tag_btag       = cms.untracked.InputTag('simpleSecondaryVertexBJetTags'),
        tag_jetTracks  = cms.untracked.InputTag('ak5JetTracksAssociatorAtVertex'),
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
	tag_genjets    = cms.untracked.InputTag('sisCone5GenJets'),
	tag_l1trig     = cms.untracked.InputTag("gtDigis"),
	tag_hlttrig    = cms.untracked.InputTag("TriggerResults","","HLT"),
        
	# Jet ID configuration
	jetID = cms.PSet(),

        # Event Selection Criteria
	# Muons:
	sel_minmupt     = cms.double(5.0),
	sel_maxmueta    = cms.double(2.4),
	# Electrons:
	sel_minelpt     = cms.double(5.0),
	sel_maxeleta    = cms.double(2.4),
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
	sel_maxphoeta   = cms.double(2.4),

        # Additional jet collections
        jets = cms.VPSet()

)
