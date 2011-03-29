import FWCore.ParameterSet.Config as cms

analyze = cms.EDAnalyzer('NTupleProducer',
	# Main settings
	isRealData       = cms.untracked.bool(False),
	# Collections
	tag_muons        = cms.untracked.InputTag('muons'),
	tag_pfmuons        = cms.untracked.InputTag('patMuonsPF'),
	tag_muisodeptk   = cms.untracked.InputTag('muIsoDepositTk'),
	tag_muisodepec   = cms.untracked.InputTag("muIsoDepositCalByAssociatorTowers","ecal"),
	tag_muisodephc   = cms.untracked.InputTag("muIsoDepositCalByAssociatorTowers","hcal"),
	tag_electrons    = cms.untracked.InputTag('gsfElectrons'),
	tag_pfelectrons  = cms.untracked.InputTag('patElectronsPF'),
	tag_pftaus       = cms.untracked.InputTag('patTausPF'),
	tag_elidWP       = cms.untracked.string('simpleEleId90relIso'),
	tag_jets         = cms.untracked.InputTag('ak5PFJets'),
	jetCorrs         = cms.untracked.string('ak5PFL1FastL2L3'),
	tag_btag1        = cms.untracked.InputTag('trackCountingHighEffBJetTags'),
	tag_btag2        = cms.untracked.InputTag('trackCountingHighPurBJetTags'),
	tag_btag3        = cms.untracked.InputTag('simpleSecondaryVertexHighEffBJetTags'),
	tag_btag4        = cms.untracked.InputTag('simpleSecondaryVertexHighPurBJetTags'),
	tag_rawcalomet   = cms.untracked.InputTag('met'),
	tag_tcmet        = cms.untracked.InputTag('tcMet'),
	tag_pfmet        = cms.untracked.InputTag('pfMet'),
	tag_corrcalomet  = cms.untracked.InputTag('metMuonJESCorAK5'),
	tag_genmet       = cms.untracked.InputTag('genMetTrue'),
	tag_vertex       = cms.untracked.InputTag('offlinePrimaryVertices'),
	tag_tracks       = cms.untracked.InputTag('generalTracks'),
	tag_photons      = cms.untracked.InputTag('photons'),
	tag_caltow       = cms.untracked.InputTag('towerMaker'),
	tag_EBrechits    = cms.untracked.InputTag('reducedEcalRecHitsEB'),
	tag_EErechits    = cms.untracked.InputTag('reducedEcalRecHitsEE'),
	tag_genpart      = cms.untracked.InputTag('genParticles'),
	tag_genjets      = cms.untracked.InputTag('ak5GenJets'),
	tag_l1trig       = cms.untracked.InputTag("gtDigis"),
	tag_hlttrigevent = cms.untracked.InputTag("hltTriggerSummaryAOD"),
	tag_hcalnoise    = cms.untracked.InputTag('HBHENoiseFilterResultProducer','HBHENoiseFilterResult'),
	tag_srcRho       = cms.untracked.InputTag('kt6PFJets','rho'),

	# Trigger paths to store the triggering object information of
	hlt_labels = cms.untracked.vstring('hltSingleMu3L3Filtered3',
	                              'hltSingleMu5L3Filtered5',
	                              'hltSingleMu9L3Filtered9',
	                              'hltL1NonIsoHLTNonIsoSingleElectronLWEt10EleIdDphiFilteroHLTNonIsoSingleElectronLWEt10PixelMatchFilter',
	                              'hltL1NonIsoHLTNonIsoSingleElectronLWEt10EleIdDphiFilter',
	                              'hltL1NonIsoHLTNonIsoSingleElectronLWEt15PixelMatchFilter',
	                              'hltL1NonIsoHLTNonIsoSinglePhotonEt10HcalIsolFilter'),

	# Event Selection Criteria
	# Muons:
	sel_minmupt       = cms.double(5.0),
	sel_maxmueta      = cms.double(2.5),
	# Electrons:
	sel_minelpt       = cms.double(5.0),
	sel_maxeleta      = cms.double(2.5),
	# Jets:
	sel_mincorjpt     = cms.double(20.0),
	sel_minrawjpt     = cms.double(0.0),
	sel_maxjeta       = cms.double(10.0),
	sel_minjemfrac    = cms.double(0.0),
	# Tracks:
	sel_mintrkpt      = cms.double(1.0),
	sel_maxtrketa     = cms.double(10.0),
	sel_maxtrknchi2   = cms.double(1e15),
	sel_mintrknhits   = cms.int32(0),
	# Photons
	sel_minphopt      = cms.double(5.0),
	sel_maxphoeta     = cms.double(2.5),
	# GenLeptons
	sel_mingenleptpt  = cms.double(2.0),
	sel_maxgenlepteta = cms.double(10),
	# GenJets
	sel_mingenjetpt  = cms.double(10.0),
	sel_maxgenjeteta = cms.double(6.0),

	btag_matchdeltaR = cms.double(0.25),

	# EB rechits
        sel_fminebrechitE = cms.double(20.),

	# Additional jet collections
	jets = cms.VPSet()

)

