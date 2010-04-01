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

### Switch for type of run (data, MC) and reconstruction (RECO, PAT, PF) #######
runon = 'data'
# runon = 'MC31x'
# runon = 'MC34x'
recoType = 'RECO'
# recoType = 'PAT'
# recoType = 'PF'

### Running conditions #########################################################
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
if runon=='data':
    process.GlobalTag.globaltag = "GR10_P_V4::All"
else:
    # CMSSW_3_4_X:
    # process.GlobalTag.globaltag = "MC_3XY_V18::All"
    # CMSSW_3_3_X:
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
      fileNames = cms.untracked.vstring(
     '/store/data/Commissioning10/MinimumBias/RAW-RECO/v7/000/132/440/0636C91D-4C3C-DF11-9A45-001A649747B0.root'
    ),
#Enable if you see duplicate error      duplicateCheckMode = cms.untracked.string("noDuplicateCheck")
)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(500) )
# Output
process.TFileService = cms.Service("TFileService",
# keep track of the type of data source and reco type in the ntuple file name
	fileName = cms.string("NTupleProducer_35X_"+runon+"_"+recoType+".root"),  
	closeFileFast = cms.untracked.bool(True)
)

#### Parameterisation for Jet Corrections and JES ME Corrections ###############
if runon=='MC31x':
    recoJet_src = "antikt5CaloJets"
    genJet_src = "antikt5GenJets"
else:
    recoJet_src = "ak5CaloJets"
    genJet_src = "ak5GenJets"

#### Jet Corrections ###########################################################
# See https://twiki.cern.ch/twiki/bin/view/CMS/WorkBookJetAnalysis#CorrInstructions
# to check what cff to use
if runon=='data':
    process.load("JetMETCorrections.Configuration.L2L3Corrections_900GeV_cff")
else:
    process.load("JetMETCorrections.Configuration.L2L3Corrections_Summer09_cff")
process.L2L3CorJetAK5Calo.src = recoJet_src
# NB: also check the analysis input below.

### JES MET Corrections ########################################################
from JetMETCorrections.Type1MET.MetType1Corrections_cff import metJESCorAK5CaloJet

process.metMuonJESCorAK5 = metJESCorAK5CaloJet.clone()
process.metMuonJESCorAK5.inputUncorJetsLabel = recoJet_src
process.metMuonJESCorAK5.corrector = "L2L3JetCorrectorAK5Calo"
process.metMuonJESCorAK5.inputUncorMetLabel = "corMetGlobalMuons"
process.metCorSequence = cms.Sequence(process.metMuonJESCorAK5)

### Analysis configuration #####################################################
process.load("DiLeptonAnalysis.NTupleProducer.ntupleproducer_cfi")
process.analyze.tag_jets   = recoJet_src
process.analyze.isRealData = cms.untracked.bool(runon=='data')
# Synchronise with Jet Corrections above
process.analyze.jetCorrs   = 'L2L3JetCorrectorAK5Calo' 

from RecoJets.JetProducers.JetIDParams_cfi import JetIDParams
process.analyze.jetID = JetIDParams

#### Path ######################################################################
process.mybtag = cms.Sequence(   process.impactParameterTagInfos
                               * process.simpleSecondaryVertexBJetTags )
process.p = cms.Path(
      process.L2L3CorJetAK5Calo
    + process.metCorSequence
    + process.mybtag
    + process.analyze
    )
