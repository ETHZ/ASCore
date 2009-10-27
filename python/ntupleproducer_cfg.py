import FWCore.ParameterSet.Config as cms
import HLTrigger.HLTfilters.hltHighLevel_cfi as hlt

process = cms.Process("NTupleProducer")

# This stuff is for the MET correction...
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Geometry.CommonDetUnit.globalTrackingGeometry_cfi")
process.load("Geometry.CommonDetUnit.bareGlobalTrackingGeometry_cfi")
# process.GlobalTag.globaltag="IDEAL_V5::All"
process.GlobalTag.globaltag="MC_31X_V3::All"
process.load("TrackingTools.TrackAssociator.default_cfi")
process.load("TrackingTools.TrackAssociator.DetIdAssociatorESProducer_cff")
process.load("FWCore.MessageService.MessageLogger_cfi")
### b-tagging
process.load("RecoBTag.Configuration.RecoBTag_cff")

process.MessageLogger = cms.Service("MessageLogger",
)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
	# replace 'myfile.root' with the source file you want to use
	fileNames = cms.untracked.vstring(
		# This is from /TTJets-madgraph/Fall08_IDEAL_V11_redigi_v10/AODSIM:
		# 'file:/data/14/stiegerb/SUSY/31X/ppMuXLoose_Summer09-MC_31X_V3-v1_GEN-SIM-RECO.root'
		'file:/data/14/stiegerb/SUSY/31X/ppEleX_Summer09-MC_31X_V3-v1_GEN-SIM-RECO.root'
		# 'file:EleIsoTest.root'
		# This is from /SUSY_LM0-sftsht/Summer08_IDEAL_V11_v1/GEN-SIM-RECO:
		# 'file:/data/fronga/LM0-RECO-312.root'
	),
	duplicateCheckMode = cms.untracked.string("noDuplicateCheck")
)

process.TFileService = cms.Service("TFileService",
	fileName = cms.string("NTupleProducer.root"),
	closeFileFast = cms.untracked.bool(True)
)

############# Jet Corrections ##################################################
process.load("JetMETCorrections.Configuration.L2L3Corrections_Summer09_cff")
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
# process.prefer("L2JetCorrectorSC5Calo")

############# b-tagging for SC (antikt) ########################################
process.impactParameterTagInfos.jetTracks = cms.InputTag("sisCone5JetTracksAssociatorAtVertex")

### JES MET Corrections ########################################################
from JetMETCorrections.Configuration.L2L3Corrections_Summer09_cff import *
from JetMETCorrections.Type1MET.MetType1Corrections_cff import metJESCorSC5CaloJet

process.metMuonJESCorSC5 = metJESCorSC5CaloJet.clone()
process.metMuonJESCorSC5.inputUncorJetsLabel = "sisCone5CaloJets"
process.metMuonJESCorSC5.corrector = "L2L3JetCorrectorSC5Calo"
process.metMuonJESCorSC5.inputUncorMetLabel = "corMetGlobalMuons"

process.metCorSequence = cms.Sequence(process.metMuonJESCorSC5)

############# Egamma Isolation #################################################
# Produce eleIsoDeposits first!
from RecoEgamma.EgammaIsolationAlgos.eleTrackExtractorBlocks_cff import *
process.eleIsoDepositTk = cms.EDProducer("CandIsoDepositProducer",
    src = cms.InputTag("gsfElectrons"),
    trackType = cms.string('candidate'),
    MultipleDepositsFlag = cms.bool(False),
    ExtractorPSet = cms.PSet(EleIsoTrackExtractorBlock)
)
from RecoEgamma.EgammaIsolationAlgos.eleEcalExtractorBlocks_cff import *
process.eleIsoDepositEcalFromHits = cms.EDProducer("CandIsoDepositProducer",
    src = cms.InputTag("gsfElectrons"),
    trackType = cms.string('candidate'),
    MultipleDepositsFlag = cms.bool(False),
    ExtractorPSet = cms.PSet(EleIsoEcalFromHitsExtractorBlock)
)
from RecoEgamma.EgammaIsolationAlgos.eleHcalExtractorBlocks_cff import *
process.eleIsoDepositHcalFromTowers = cms.EDProducer("CandIsoDepositProducer",
    src = cms.InputTag("gsfElectrons"),
    trackType = cms.string('candidate'),
    MultipleDepositsFlag = cms.bool(False),
    ExtractorPSet = cms.PSet(EleIsoHcalFromTowersExtractorBlock)
)
# from RecoEgamma.EgammaIsolationAlgos.eleIsoDeposits_cff import *
# from RecoEgamma.EgammaIsolationAlgos.eleIsoFromDepsModules_cff  import *
process.eleIsoDeposits = cms.Sequence(
    process.eleIsoDepositTk +
    process.eleIsoDepositEcalFromHits +
    process.eleIsoDepositHcalFromTowers
)
#################################
process.eleIsoFromDepsTk = cms.EDProducer("CandIsolatorFromDeposits",
	deposits = cms.VPSet(
		cms.PSet(
			mode = cms.string('sum'),
			src = cms.InputTag("eleIsoDepositTk"),
			weight = cms.string('1'),
			deltaR = cms.double(0.3),
			vetos = cms.vstring('0.04','Threshold(0.7)'),
			skipDefaultVeto = cms.bool(True)
		)
	)
)

process.eleIsoFromDepsEcalFromHits= cms.EDProducer("CandIsolatorFromDeposits",
	deposits = cms.VPSet(
		cms.PSet(
			mode = cms.string('sum'),
			src = cms.InputTag("eleIsoDepositEcalFromHits"),
			weight = cms.string('1'),
			deltaR = cms.double(0.3),
			vetos = cms.vstring('EcalBarrel:0.045',
			            'EcalBarrel:RectangularEtaPhiVeto(-0.02,0.02,-0.5,0.5)',
			            'EcalBarrel:ThresholdFromTransverse(0.08)',
                     'EcalEndcaps:ThresholdFromTransverse(0.3)',
			            'EcalEndcaps:0.070',
			            'EcalEndcaps:RectangularEtaPhiVeto(-0.02,0.02,-0.5,0.5)'), 
			skipDefaultVeto = cms.bool(True)
		)
	)
)

process.eleIsoFromDepsHcalFromTowers = cms.EDProducer("CandIsolatorFromDeposits",
	deposits = cms.VPSet(
		cms.PSet(
	   	src             = cms.InputTag("eleIsoDepositHcalFromTowers"),
	   	deltaR          = cms.double(0.4),
	   	weight          = cms.string('1'),
	   	vetos           = cms.vstring('0.15'),
	   	skipDefaultVeto = cms.bool(True),
	   	mode            = cms.string('sum')
		)
	)
)

process.eleIsoFromDeposits = cms.Sequence(
	process.eleIsoFromDepsTk * 
	process.eleIsoFromDepsEcalFromHits * 
	process.eleIsoFromDepsHcalFromTowers
)
 

############# User Analysis ####################################################
process.analyze = cms.EDAnalyzer('NTupleProducer',
	# Collections
	tag_muons   = cms.untracked.InputTag('muons'),
	# tag_electrons = cms.untracked.InputTag('pixelMatchGsfElectrons'),
	tag_electrons = cms.untracked.InputTag('gsfElectrons'),
	
	tag_elisotk = cms.untracked.InputTag('eleIsoFromDepsTk::NTupleProducer'),
	tag_elisoec = cms.untracked.InputTag('eleIsoFromDepsEcalFromHits::NTupleProducer'),
	tag_elisohc = cms.untracked.InputTag('eleIsoFromDepsHcalFromTowers::NTupleProducer'),
	tag_elisodeptk = cms.untracked.InputTag('eleIsoDepositTk::NTupleProducer'),
	tag_elisodepec = cms.untracked.InputTag('eleIsoDepositEcalFromHits::NTupleProducer'),
	tag_elisodephc = cms.untracked.InputTag('eleIsoDepositHcalFromTowers::NTupleProducer'),
	
	tag_muisodeptk = cms.untracked.InputTag('muIsoDepositTk'),
	tag_muisodepec = cms.untracked.InputTag("muIsoDepositCalByAssociatorTowers","ecal"),
	tag_muisodephc = cms.untracked.InputTag("muIsoDepositCalByAssociatorTowers","hcal"),
	tag_sc      = cms.untracked.InputTag('correctedHybridSuperClusters'),
	tag_jets    = cms.untracked.InputTag('sisCone5CaloJets'),
        tag_btag    = cms.untracked.InputTag('simpleSecondaryVertexBJetTags'),  #trackCountingHighPurBJetTags #jetProbabilityBJetTags
	tag_met1    = cms.untracked.InputTag('met'),
	tag_met2    = cms.untracked.InputTag('corMetGlobalMuons'),
	tag_met3    = cms.untracked.InputTag('tcMet'),
	tag_met4    = cms.untracked.InputTag('pfMet'),
	tag_met5    = cms.untracked.InputTag('metMuonJESCorSC5'),
	tag_vertex  = cms.untracked.InputTag('offlinePrimaryVertices'),
	tag_tracks  = cms.untracked.InputTag('generalTracks'),
	tag_caltow  = cms.untracked.InputTag('towerMaker'),
	tag_genpart = cms.untracked.InputTag('genParticles'),
	tag_triggers   = cms.InputTag("TriggerResults","","HLT"),
	
	# Event Selection Criteria
	# Muons:
	sel_minmupt    = cms.untracked.double(5.0),
	sel_maxmueta   = cms.untracked.double(2.4),
	# Electrons:
	sel_minelpt    = cms.untracked.double(5.0),
	sel_maxeleta   = cms.untracked.double(2.5),
	sel_maxeliso   = cms.untracked.double(100.0),
	sel_maxeld0    = cms.untracked.double(100.0),
	# Jets:
	sel_minjpt     = cms.untracked.double(20.0),
	sel_maxjeta    = cms.untracked.double(5.0),
	sel_minjemfrac = cms.untracked.double(0.0),
	# Isolation Parameters
	iso_MuTkDRin   = cms.untracked.double(0.015),
	iso_MuTkDRout  = cms.untracked.double(0.3),
	iso_MuTkSeed   = cms.untracked.double(0.1),
	iso_MuCalDRin  = cms.untracked.double(0.0),
	iso_MuCalDRout = cms.untracked.double(0.3),
	iso_MuCalSeed  = cms.untracked.double(0.1)
)

############# Path #############################################################
process.p = cms.Path(process.L2L3CorJetSC5Calo + process.metCorSequence)
mybtag = cms.Sequence(process.impactParameterTagInfos*process.simpleSecondaryVertexBJetTags)
process.l = cms.Path(mybtag)
process.o = cms.EndPath(process.eleIsoDeposits + process.eleIsoFromDeposits + process.analyze)
