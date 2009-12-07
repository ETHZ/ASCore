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
		'file:/home/xv/stiegerb/FirstData/Run123596/bit40-123596-stdReco-GR09_P_V7-lumi1-68.root'
		# 'file:/data/theofil08/Data/Beam09/bit40or41skim_expressPhysics_run123596_full.root'
		# 'file:/home/xv/stiegerb/Skims/python/CollisionSkim_123596_1.root'
		# 'file:/home/xv/stiegerb/FirstData/Run123592/241D28BF-1EE2-DE11-9CF3-001617E30CC8.root',
		# 'file:/home/xv/stiegerb/FirstData/Run123592/5865B623-1EE2-DE11-82E2-001D09F2545B.root',
		# 'file:/home/xv/stiegerb/FirstData/Run123592/5E32F723-1EE2-DE11-ACC5-001D09F252E9.root',
		# 'file:/home/xv/stiegerb/FirstData/Run123592/68BAA366-1FE2-DE11-A742-001617C3B66C.root',
		# 'file:/home/xv/stiegerb/FirstData/Run123592/A6C66722-1EE2-DE11-BCFC-001D09F25208.root',
		# 'file:/home/xv/stiegerb/FirstData/Run123592/D0ACE867-1FE2-DE11-AF2E-000423D94A20.root',
		# 'file:/home/xv/stiegerb/FirstData/Run123592/CA74FE6C-1FE2-DE11-9C9A-001617E30D4A.root',
		# 'file:/home/xv/stiegerb/FirstData/Run123592/761A0BBF-1EE2-DE11-8E1B-001617C3B76E.root',
		# 'file:/home/xv/stiegerb/FirstData/Run123592/863F0D1D-20E2-DE11-9D98-001617C3B706.root',
		# 'file:/home/xv/stiegerb/FirstData/Run123592/98B3021F-20E2-DE11-A8A9-0030487D0D3A.root',
		# 'file:/home/xv/stiegerb/FirstData/Run123592/F25C511D-20E2-DE11-82D0-001617C3B6DE.root',
		# 'file:/home/xv/stiegerb/FirstData/Run123592/ACA6611D-20E2-DE11-B804-001617E30CC8.root',
		# 'file:/home/xv/stiegerb/FirstData/Run123592/80B9DACB-20E2-DE11-9526-001617C3B706.root'

	),
	duplicateCheckMode = cms.untracked.string("noDuplicateCheck")
)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )
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

from RecoJets.JetProducers.JetIDParams_cfi import JetIDParams
process.analyze.jetID = JetIDParams

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
