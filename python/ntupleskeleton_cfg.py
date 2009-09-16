import FWCore.ParameterSet.Config as cms

process = cms.Process("NTupleSkeleton")

process.MessageLogger = cms.Service("MessageLogger",
	# destinations = cms.untracked.vstring('output.txt')
)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
	# replace with the source file you want to use
	fileNames = cms.untracked.vstring(
		'file:scratch/ppEleX_Summer09-MC_31X_V3-v1_GEN-SIM-RECO.root'
	)
)

process.TFileService = cms.Service("TFileService",
	fileName = cms.string("NTupleSkeleton.root"),
	closeFileFast = cms.untracked.bool(True)
)

############# Jet Corrections ##################################################
process.load("JetMETCorrections.Configuration.L2L3Corrections_Winter09_cff")
process.L2JetCorrectorSC5Calo = cms.ESSource("L2RelativeCorrectionService", 
	tagName = cms.string('Summer08Redigi_L2Relative_SC5Calo'),
	label = cms.string('L2RelativeJetCorrectorSC5Calo')
)
process.L3JetCorrectorSC5Calo = cms.ESSource("L3AbsoluteCorrectionService", 
	tagName = cms.string('Summer08Redigi_L3Absolute_SC5Calo'),
	label = cms.string('L3AbsoluteJetCorrectorSC5Calo')
)
process.L2L3CorJetSC5Calo = cms.EDProducer("CaloJetCorrectionProducer",
    src = cms.InputTag("sisCone5CaloJets"),
    correctors = cms.vstring('L2L3JetCorrectorSC5Calo')
)
process.prefer("L2JetCorrectorSC5Calo")

############# User Analysis ####################################################
process.analyze = cms.EDAnalyzer('NTupleSkeleton',
	# Collections
	tag_muons   = cms.untracked.InputTag('muons'),
	tag_electrons = cms.untracked.InputTag('gsfElectrons'),
	tag_jets    = cms.untracked.InputTag('sisCone5CaloJets'),
	tag_met     = cms.untracked.InputTag('met'),
	# Event Selection Criteria
	# Muons:
	sel_minmupt    = cms.untracked.double(5.0),
	sel_maxmueta   = cms.untracked.double(2.4),
	# Electrons:
	sel_minelpt    = cms.untracked.double(5.0),
	sel_maxeleta   = cms.untracked.double(2.5),
	# Jets:
	sel_minjpt     = cms.untracked.double(30.0),
	sel_maxjeta    = cms.untracked.double(5.0),
)

############# Path #############################################################
process.p = cms.Path(process.L2L3CorJetSC5Calo)
process.o = cms.EndPath(process.analyze)
