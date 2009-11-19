import FWCore.ParameterSet.Config as cms
import HLTrigger.HLTfilters.hltHighLevel_cfi as hlt

process = cms.Process("NTupleProducer")

### Message Logger #############################################################
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.categories.append('NTP')
process.MessageLogger.cerr.NTP = cms.untracked.PSet(
    limit = cms.untracked.int32(-1),
    reportEvery = cms.untracked.int32(1)
    )
process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

### Running conditions #########################################################
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag="MC_31X_V3::All"

### b-tagging ##################################################################
process.load("Geometry.CommonDetUnit.globalTrackingGeometry_cfi")
process.load("Geometry.CommonDetUnit.bareGlobalTrackingGeometry_cfi")
process.load("TrackingTools.TrackAssociator.default_cfi")
process.load("TrackingTools.TrackAssociator.DetIdAssociatorESProducer_cff")
process.load("RecoBTag.Configuration.RecoBTag_cff")
# b-tagging for SC5
process.impactParameterTagInfos.jetTracks = cms.InputTag("sisCone5JetTracksAssociatorAtVertex")

### Input/Output ###############################################################
# Input
process.source = cms.Source("PoolSource",
	# replace 'myfile.root' with the source file you want to use
	fileNames = cms.untracked.vstring(
		# This is from /TTJets-madgraph/Fall08_IDEAL_V11_redigi_v10/AODSIM:
		# 'file:/data/14/stiegerb/SUSY/31X/ppMuXLoose_Summer09-MC_31X_V3-v1_GEN-SIM-RECO.root'
		#'file:/data/14/stiegerb/SUSY/31X/ppEleX_Summer09-MC_31X_V3-v1_GEN-SIM-RECO.root'
		# 'file:EleIsoTest.root'
		# This is from /SUSY_LM0-sftsht/Summer08_IDEAL_V11_v1/GEN-SIM-RECO:
		 'file:/data/14/stiegerb/SUSY/31X/ppMuXLoose_Summer09-MC_31X_V3-v1_GEN-SIM-RECO.root'
	),
	duplicateCheckMode = cms.untracked.string("noDuplicateCheck")
)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )
# Output
process.TFileService = cms.Service("TFileService",
	fileName = cms.string("NTupleProducer.root"),
	closeFileFast = cms.untracked.bool(True)
)

#### Jet Corrections ###########################################################
process.load("JetMETCorrections.Configuration.L2L3Corrections_Summer09_cff")
process.L2JetCorrectorSC5Calo = cms.ESSource("L2RelativeCorrectionService", 
	tagName = cms.string('Summer09_L2Relative_SC5Calo'),
	label = cms.string('L2RelativeJetCorrectorSC5Calo')
)
process.L3JetCorrectorSC5Calo = cms.ESSource("L3AbsoluteCorrectionService", 
	tagName = cms.string('Summer09_L3Absolute_SC5Calo'),
	label = cms.string('L3AbsoluteJetCorrectorSC5Calo')
)
process.L2L3CorJetSC5Calo = cms.EDProducer("CaloJetCorrectionProducer",
    src = cms.InputTag("sisCone5CaloJets"),
    correctors = cms.vstring('L2L3JetCorrectorSC5Calo')
)
# process.prefer("L2JetCorrectorSC5Calo")

### JES MET Corrections ########################################################
from JetMETCorrections.Type1MET.MetType1Corrections_cff import metJESCorSC5CaloJet

process.metMuonJESCorSC5 = metJESCorSC5CaloJet.clone()
process.metMuonJESCorSC5.inputUncorJetsLabel = "sisCone5CaloJets"
process.metMuonJESCorSC5.corrector = "L2L3JetCorrectorSC5Calo"
process.metMuonJESCorSC5.inputUncorMetLabel = "corMetGlobalMuons"

process.metCorSequence = cms.Sequence(process.metMuonJESCorSC5)

### Egamma Isolation ###########################################################
# Produce eleIsoDeposits first!
process.load("RecoEgamma.EgammaIsolationAlgos.eleIsoDeposits_cff")
# Make isolation from deposits
process.load("RecoEgamma.EgammaIsolationAlgos.eleIsoFromDeposits_cff")
# Only keep the three main eleIsoDeposits
process.eleIsoDeposits.remove(process.eleIsoDepositHcalDepth1FromTowers)
process.eleIsoDeposits.remove(process.eleIsoDepositHcalDepth2FromTowers)
process.eleIsoFromDeposits.remove(process.eleIsoFromDepsHcalDepth1FromTowers)
process.eleIsoFromDeposits.remove(process.eleIsoFromDepsHcalDepth2FromTowers)

# Example configuration
process.eleIsoFromDepsEcalFromHits.deposits[0].deltaR = 0.3

### Analysis configuration #####################################################
process.load("DiLeptonAnalysis.NTupleProducer.ntupleproducer_cfi")

#### Path ######################################################################
process.mybtag = cms.Sequence(   process.impactParameterTagInfos
                               * process.simpleSecondaryVertexBJetTags )
process.p = cms.Path(
      process.L2L3CorJetSC5Calo
    + process.metCorSequence
    + process.mybtag
    + process.eleIsoDeposits
    + process.eleIsoFromDeposits
    + process.analyze
    )

