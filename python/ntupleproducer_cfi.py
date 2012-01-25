import FWCore.ParameterSet.Config as cms

analyze = cms.EDAnalyzer('NTupleProducer',
	# Main settings
	isRealData       = cms.untracked.bool(False),
	# Collections
	tag_muons        = cms.untracked.InputTag('muons'),
	tag_muisodeptk   = cms.untracked.InputTag('muIsoDepositTk'),
	tag_muisodepec   = cms.untracked.InputTag("muIsoDepositCalByAssociatorTowers","ecal"),
	tag_muisodephc   = cms.untracked.InputTag("muIsoDepositCalByAssociatorTowers","hcal"),
	tag_electrons    = cms.untracked.InputTag('gsfElectrons'),
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
	tag_pfmetPAT     = cms.untracked.InputTag('patMETsPF3'),
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
	tag_hcalnoise    = cms.untracked.InputTag('HBHENoiseFilterResultProducerStd','HBHENoiseFilterResult'),
	tag_hcalnoiseIso = cms.untracked.InputTag('HBHENoiseFilterResultProducerIso','HBHENoiseFilterResult'),
	tag_srcRho       = cms.untracked.InputTag('kt6PFJets','rho'),
	tag_srcRhoPFnoPU = cms.untracked.InputTag('kt6PFJetsPF3','rho'),
	tag_pfphotonsProducer  = cms.untracked.InputTag("pfPhotonTranslator:pfphot"),
	tag_pfProducer = cms.untracked.InputTag("particleFlow"),
        tag_SC_barrel    = cms.untracked.InputTag("correctedHybridSuperClusters"),
        tag_SC_endcap    = cms.untracked.InputTag("correctedMulti5x5SuperClustersWithPreshower"),

        tag_doVertexing = cms.untracked.bool(False),
        tag_fTrackCollForVertexing = cms.untracked.InputTag("generalTracks"),
        tag_fallConversionsCollForVertexing = cms.untracked.InputTag("allConversions"),
        tag_perVtxMvaWeights = cms.untracked.string("/shome/peruzzi/localafs/afs/cern.ch/user/m/musella/public/higgs/vertex_likelihoods/TMVAClassification_BDTCat_conversions_tmva_407.weights.xml"),
        tag_perVtxMvaMethod = cms.untracked.string("BDTCat"),
        tag_perEvtMvaWeights = cms.untracked.string("/shome/peruzzi/localafs/afs/cern.ch/user/m/musella/public/higgs/vertex_likelihoods/TMVAClassification_evtBDTG_conversions_tmva_407.weights.xml"),
        tag_perEvtMvaMethod = cms.untracked.string("evtBDTG"),

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
        # SC
        sel_minSCraw      = cms.double(5.0),
	# GenLeptons
	sel_mingenleptpt  = cms.double(0.0),
	sel_maxgenlepteta = cms.double(100),
	# GenPhotons
	sel_mingenphotpt = cms.double(5.0),
	sel_maxgenphoteta = cms.double(2.5),
	# GenJets
	sel_mingenjetpt  = cms.double(10.0),
	sel_maxgenjeteta = cms.double(6.0),

	btag_matchdeltaR = cms.double(0.25),

	# EB rechits
        sel_fminebrechitE = cms.double(20.),

	# Additional jet collections
	jets = cms.VPSet(),
        leptons = cms.VPSet(),

	# tag pile up distributions: replace empty strings in order to calculate in time and OOT pileup weights
	pu_data = cms.untracked.vstring('', ''), # replace this by cms.untracked.vstring('data_pileup.root', 'name_of_histo')
	pu_mc   = cms.untracked.vstring('', '')  # replace this by cms.untracked.vstring('mc_pileup.root'  , 'name_of_histo')

)

