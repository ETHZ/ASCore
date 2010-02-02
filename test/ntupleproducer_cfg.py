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
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(False) )

### Running conditions #########################################################
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = "MC_31X_V3::All"

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
		'file:scratch/TTbar4Jets_40GeVthreshold-alpgen-GEN-SIM-RECO-MC_31X_V3_7TeV-v3.root'
		# 'file:/data26/papel/ttbar_3_1_2.root'
	),
	duplicateCheckMode = cms.untracked.string("noDuplicateCheck")
)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
# Output
process.TFileService = cms.Service("TFileService",
	fileName = cms.string("NTupleProducer.root"),
	closeFileFast = cms.untracked.bool(True)
)

#### Jet Corrections ###########################################################
# process.load("JetMETCorrections.Configuration.L2L3Corrections_900GeV_cff")
# process.load("JetMETCorrections.Configuration.L2L3Corrections_2360GeV_cff")
process.load("JetMETCorrections.Configuration.L2L3Corrections_Summer09_cff")
## antiKt5
process.L2JetCorrectorAK5Calo = cms.ESSource("L2RelativeCorrectionService", 
	tagName = cms.string('Summer09_L2Relative_AK5Calo'),
	# tagName = cms.string('900GeV_L2Relative_AK5Calo'),
	label = cms.string('L2RelativeJetCorrectorAK5Calo')
)
process.L3JetCorrectorAK5Calo = cms.ESSource("L3AbsoluteCorrectionService", 
	tagName = cms.string('Summer09_L3Absolute_AK5Calo'),
	# tagName = cms.string('900GeV_L3Absolute_AK5Calo'),
	label = cms.string('L3AbsoluteJetCorrectorAK5Calo')
)
process.L2L3CorJetAK5Calo = cms.EDProducer("CaloJetCorrectionProducer",
    src = cms.InputTag("antikt5CaloJets"),
    correctors = cms.vstring('L2L3JetCorrectorAK5Calo')
)

### JES MET Corrections ########################################################
from JetMETCorrections.Type1MET.MetType1Corrections_cff import metJESCorAK5CaloJet

process.metMuonJESCorAK5 = metJESCorAK5CaloJet.clone()
process.metMuonJESCorAK5.inputUncorJetsLabel = "antikt5CaloJets"
process.metMuonJESCorAK5.corrector = "L2L3JetCorrectorAK5Calo"
process.metMuonJESCorAK5.inputUncorMetLabel = "corMetGlobalMuons"
process.metCorSequence = cms.Sequence(process.metMuonJESCorAK5)

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

from RecoJets.JetProducers.JetIDParams_cfi import JetIDParams
process.analyze.jetID = JetIDParams

#### Path ######################################################################
process.mybtag = cms.Sequence(   process.impactParameterTagInfos
                               * process.simpleSecondaryVertexBJetTags )
process.p = cms.Path(
      process.L2L3CorJetAK5Calo
      # process.L2L3CorJetSC5Calo
    + process.metCorSequence
    + process.mybtag
    + process.eleIsoDeposits
    + process.eleIsoFromDeposits
    + process.analyze
    )
