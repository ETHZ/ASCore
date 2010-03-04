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
# runon = 'data'
# runon = 'MC31x'
runon = 'MC34x'
recoType = 'RECO'
# recoType = 'PAT'
# recoType = 'PF'

### Running conditions #########################################################
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
if runon=='data':
    process.GlobalTag.globaltag = "GR09_R_34X_V3::All"
else:
    # CMSSW_3_4_X:
    process.GlobalTag.globaltag = "MC_3XY_V18::All"
    # CMSSW_3_3_X:
    #process.GlobalTag.globaltag = "MC_31X_V3::All"

### b-tagging ##################################################################
process.load("Geometry.CommonDetUnit.globalTrackingGeometry_cfi")
process.load("Geometry.CommonDetUnit.bareGlobalTrackingGeometry_cfi")
process.load("TrackingTools.TrackAssociator.default_cfi")
process.load("TrackingTools.TrackAssociator.DetIdAssociatorESProducer_cff")
process.load("RecoBTag.Configuration.RecoBTag_cff")
# b-tagging for SC5
process.impactParameterTagInfos.jetTracks = cms.InputTag("sisCone5JetTracksAssociatorAtVertex")

### Input/Output ###############################################################
# parameterization
if runon=='data':
    source_fileNames = cms.untracked.vstring(
        '/store/data/BeamCommissioning09/MinimumBias/RECO/Dec19thReReco_341_v1/0006/DA397007-23EE-DE11-A754-00261894388D.root'
        )
elif runon=='MC31x':
    source_fileNames = cms.untracked.vstring(
        # 'file:/data/14/stiegerb/SUSY/31X/TTbar4Jets_40GeVthreshold-alpgen-GEN-SIM-RECO-MC_31X_V3_7TeV-v3.root'   
        'file:/data/14/stiegerb/MC/TTbar-Summer09-MC_31X_V9_Feb15-v1-GEN-SIM-RECO/F08AC9AF-831A-DF11-9AE4-001E0BE922E2.root'
        )
elif runon=='MC34x':
    source_fileNames = cms.untracked.vstring(
        '/store/mc/Summer09/MinBias/GEN-SIM-RECO/V16D_900GeV-v1/0001/FCD70794-F216-DF11-931A-0015170AC494.root'
        #'file:/data/14/stiegerb/MC/TTbar-Summer09-MC_31X_V9_Feb15-v1-GEN-SIM-RECO/F08AC9AF-831A-DF11-9AE4-001E0BE922E2.root'
        )
# Input
process.source = cms.Source("PoolSource",
	fileNames = source_fileNames,
	duplicateCheckMode = cms.untracked.string("noDuplicateCheck")
)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )
# Output
process.TFileService = cms.Service("TFileService",
# keep track of the type of data source and reco type in the ntuple file name
	fileName = cms.string("NTupleProducer_34X_"+runon+"_"+recoType+".root"),  
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
