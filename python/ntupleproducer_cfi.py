import FWCore.ParameterSet.Config as cms

analyze = cms.EDAnalyzer('NTupleProducer',
	# Main settings
	isRealData = cms.untracked.bool(True),
	isPat      = cms.untracked.bool(False),
	# Collections
	tag_muons      = cms.untracked.InputTag('muons'),
	tag_muisodeptk = cms.untracked.InputTag('muIsoDepositTk'),
	tag_muisodepec = cms.untracked.InputTag("muIsoDepositCalByAssociatorTowers","ecal"),
	tag_muisodephc = cms.untracked.InputTag("muIsoDepositCalByAssociatorTowers","hcal"),

	tag_electrons  = cms.untracked.InputTag('gsfElectrons'),	
	tag_elisotk    = cms.untracked.InputTag('eleIsoFromDepsTk::NTupleProducer'),
	tag_elisoec    = cms.untracked.InputTag('eleIsoFromDepsEcalFromHits::NTupleProducer'),
	tag_elisohc    = cms.untracked.InputTag('eleIsoFromDepsHcalFromTowers::NTupleProducer'),
	tag_elisodeptk = cms.untracked.InputTag('eleIsoDepositTk::NTupleProducer'),
	tag_elisodepec = cms.untracked.InputTag('eleIsoDepositEcalFromHits::NTupleProducer'),
	tag_elisodephc = cms.untracked.InputTag('eleIsoDepositHcalFromTowers::NTupleProducer'),
	
	tag_sc       = cms.untracked.InputTag('correctedHybridSuperClusters'),
	tag_jets     = cms.untracked.InputTag('antikt5CaloJets'),
	# tag_jets     = cms.untracked.InputTag('sisCone5CaloJets'),
	tag_btag     = cms.untracked.InputTag('simpleSecondaryVertexBJetTags'),  
	#trackCountingHighPurBJetTags #jetProbabilityBJetTags
	tag_met1     = cms.untracked.InputTag('met'),
	tag_met2     = cms.untracked.InputTag('corMetGlobalMuons'),
	tag_met3     = cms.untracked.InputTag('tcMet'),
	tag_met4     = cms.untracked.InputTag('pfMet'),
	tag_met5     = cms.untracked.InputTag('metMuonJESCorAK5'),
	# tag_met5     = cms.untracked.InputTag('metMuonJESCorSC5'),
	tag_vertex   = cms.untracked.InputTag('offlinePrimaryVertices'),
	# tag_vertex   = cms.untracked.InputTag('offlinePrimaryVertices','','REVERTEX'),
	tag_tracks   = cms.untracked.InputTag('generalTracks'),
        tag_photons   = cms.untracked.InputTag('photons'),
	# tag_tracks   = cms.untracked.InputTag('generalTracks','','RETRACK'),
	tag_caltow   = cms.untracked.InputTag('towerMaker'),
	tag_genpart  = cms.untracked.InputTag('genParticles'),
	tag_l1trig   = cms.untracked.InputTag("gtDigis"),
	tag_hlttrig  = cms.untracked.InputTag("TriggerResults","","HLT"),

	# Jet ID configuration
	jetID = cms.PSet(),

   # Event Selection Criteria
	# Muons:
	sel_minmupt     = cms.double(1.0),
	sel_maxmueta    = cms.double(2.4),
	# Electrons:
	sel_minelpt     = cms.double(1.0),
	sel_maxeleta    = cms.double(2.5),
	sel_maxeliso    = cms.double(1e15),
	sel_maxeld0     = cms.double(1e15),
	# Jets:
	sel_minjpt      = cms.double(1.0),
	sel_maxjeta     = cms.double(10.0),
	sel_minjemfrac  = cms.double(0.0),
	# Tracks:
	sel_mintrkpt    = cms.double(1.0),
	sel_maxtrketa   = cms.double(10.0),
	sel_maxtrknchi2 = cms.double(1e15),
	sel_mintrknhits = cms.int32(0),
        # Photons
	sel_minphopt    = cms.double(1.0),
	sel_maxphoeta   = cms.double(2.5),            
	# Isolation Parameters
	iso_MuTkDRin    = cms.double(0.015),
	iso_MuTkDRout   = cms.double(0.3),
	iso_MuTkSeed    = cms.double(0.1),
	iso_MuCalDRin   = cms.double(0.0),
	iso_MuCalDRout  = cms.double(0.3),
	iso_MuCalSeed   = cms.double(0.1)
)
