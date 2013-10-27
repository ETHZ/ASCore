import FWCore.ParameterSet.Config as cms
from CMGTools.External.puJetIDAlgo_cff import full_53x,cutbased

analyze = cms.EDFilter('NTupleProducer',
	# Main settings
	isRealData       = cms.bool(False),
	isFastSim        = cms.bool(False),
        isModelScan      = cms.bool(False),
	# Collections
	tag_muons        = cms.InputTag('muons'),
        tag_muonpfisosCustom = cms.VInputTag(''), # No defaults: set in _cfg.py
	tag_muisodeptk   = cms.InputTag('muons','muIsoDepositTk'),
	tag_muisodepec   = cms.InputTag('muons','ecal'),
	tag_muisodephc   = cms.InputTag('muons','hcal'),
	tag_electrons    = cms.InputTag('gsfElectrons'),
        tag_elepfisosCustom = cms.VInputTag(''), # No defaults: set in _cfg.py
        tag_elepfisosEvent  = cms.VInputTag(''), # No defaults: set in _cfg.py
	tag_jets         = cms.InputTag('ak5PFJets'),
	jetCorrs         = cms.string('ak5PFL1FastL2L3'),
	tag_btags        = cms.VInputTag(''), # No defaults: set in _cfg.py
	tag_partonmatch  = cms.InputTag('AK5PFbyValAlgo'),
	tag_rawcalomet   = cms.InputTag('met'),
	tag_tcmet        = cms.InputTag('tcMet'),
	tag_pfmet        = cms.InputTag('pfMet'),
	tag_corrcalomet  = cms.InputTag('metMuonJESCorAK5'),
	tag_genmet       = cms.InputTag('genMetTrue'),
	tag_vertex       = cms.InputTag('offlinePrimaryVertices'), # WARNING: overwritten by default with 'goodVertices' in _cfg.py
	tag_vertex_withbs= cms.InputTag('offlinePrimaryVerticesWithBS'),
	tag_tracks       = cms.InputTag('generalTracks'),
	tag_photons      = cms.InputTag('photons'),
	tag_caltow       = cms.InputTag('towerMaker'),
	tag_EBrechits    = cms.InputTag('reducedEcalRecHitsEB'),
	tag_EErechits    = cms.InputTag('reducedEcalRecHitsEE'),
	tag_genpart      = cms.InputTag('genParticles'),
	tag_genjets      = cms.InputTag('ak5GenJets'),
	tag_l1trig       = cms.InputTag("gtDigis"),
	tag_hlttrigevent = cms.InputTag("hltTriggerSummaryAOD"),
	tag_srcRho       = cms.InputTag('kt6PFJets','rho'),
	tag_srcSigma       = cms.InputTag('kt6PFJets','sigma'),
        tag_srcRhoForIso = cms.InputTag('kt6PFJetsForIso','rho'),
	tag_pfphotonsProducer  = cms.InputTag("pfPhotonTranslator:pfphot"),
	tag_pfProducer = cms.InputTag("particleFlow"),
        tag_SC_barrel    = cms.InputTag("correctedHybridSuperClusters"),
        tag_SC_endcap    = cms.InputTag("correctedMulti5x5SuperClustersWithPreshower"),
                         
	# Trigger paths to store the triggering object information of
        hlt_labels = cms.vstring('HLT_IsoMu24_eta2p1_v',
                                 'HLT_Ele27_WP80_v',
                                 'HLT_Photon36_R9Id85_OR_CaloId10_Iso50_Photon22_R9Id85_OR_CaloId10_Iso50_v'),
	# Event Selection Criteria
	# Muons:
	sel_minmupt       = cms.double(5.0),
	sel_maxmueta      = cms.double(2.5),
	# Taus:
	sel_mintaupt      = cms.double(10.0),
	sel_maxtaueta     = cms.double(2.5),
	# Electrons:
	sel_minelpt       = cms.double(5.0),
	sel_maxeleta      = cms.double(2.5),
	# Jets:
	sel_mincorjpt     = cms.double(15.0),
	sel_minrawjpt     = cms.double(0.0),
	sel_maxjeta       = cms.double(10.0),
	sel_minjemfrac    = cms.double(0.0),
	# Tracks:
	sel_mintrkpt      = cms.double(1.0),
	sel_maxtrketa     = cms.double(10.0),
	sel_maxtrknchi2   = cms.double(1e15),
	sel_mintrknhits   = cms.int32(0),
	# Photons
	sel_minphopt      = cms.double(10.0),
	sel_maxphoeta     = cms.double(3.0),
        # SC
        sel_minSCraw      = cms.double(10.0),
        sel_minSCrawPt    = cms.double(10.0),
	# GenLeptons
	sel_mingenleptpt  = cms.double(0.0),
	sel_maxgenlepteta = cms.double(100),
	# GenPhotons
	sel_mingenphotpt = cms.double(5.0),
	sel_maxgenphoteta = cms.double(3.0),
	# GenJets
	sel_mingenjetpt  = cms.double(10.0),
	sel_maxgenjeteta = cms.double(6.0),
        # PFCandidates
        sel_maxpfcandeta = cms.double(3.0),

	# EB rechits
        sel_fminebrechitE = cms.double(20.),

	# Additional collections
	jets    = cms.VPSet(),
        leptons = cms.VPSet(),
        pfCandidates = cms.VPSet(),

	# tag pile up distributions: replace empty strings in order to calculate in time and OOT pileup weights
	pu_data = cms.vstring('', ''), # replace this by cms.vstring('data_pileup.root', 'name_of_histo')
	pu_mc   = cms.vstring('', ''), # replace this by cms.vstring('mc_pileup.root'  , 'name_of_histo')

        tag_doPhotonStuff = cms.bool(False), # overwritten from test/ntupleproducer_cfg.py

        tag_fTrackCollForVertexing = cms.InputTag("generalTracks"),
        tag_fallConversionsCollForVertexing = cms.InputTag("allConversions"),
        tag_regressionVersion = cms.int32(5), # turned off by default; use version number 5 for 2012 @ 8 TeV, number 8 for 2011 @ 7 TeV

        tag_QGSyst = cms.string("pythia"),
        tag_puJetIDAlgos = cms.VPSet(cutbased,full_53x),
        tag_WeightsPhotonIDMVA_EB = cms.string("2013FinalPaper_PhotonID_Barrel_BDT_TrainRangePT15.weights.xml"),
        tag_WeightsPhotonIDMVA_EE = cms.string("2013FinalPaper_PhotonID_Endcap_BDT_TrainRangePT15.weights.xml"),

)



