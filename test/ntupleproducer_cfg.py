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
runon = 'MC35x'
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
    process.GlobalTag.globaltag = "START3X_V25B::All"

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
#     '/store/data/Commissioning10/MinimumBias/RAW-RECO/v7/000/132/440/0636C91D-4C3C-DF11-9A45-001A649747B0.root'
#     '/store/mc/Spring10/MinBias/GEN-SIM-RECO/START3X_V25B_356ReReco-v1/0004/F44C3A4D-F73B-DF11-B42D-003048678B94.root'
      'file:////data/susy/data/MinBias-Spring10-START3X_V25B_356ReReco-v1-GEN-SIM-RECO.root'
    ),
#Enable if you see duplicate error      duplicateCheckMode = cms.untracked.string("noDuplicateCheck")
)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(500) )
# Output
process.TFileService = cms.Service("TFileService",
# keep track of the type of data source and reco type in the ntuple file name
	fileName = cms.string("NTupleProducer_36X_"+runon+"_"+recoType+".root"),  
	closeFileFast = cms.untracked.bool(True)
)

#### Parameterisation for Jet Corrections and JES ME Corrections ###############
recoJet_src = "ak5CaloJets"
genJet_src = "ak5GenJets"

# Jet ID: add the ones you want
process.load('RecoJets.Configuration.JetIDProducers_cff')
process.recoJetIdSequence = cms.Sequence( process.ak5JetID )

#### Jet Corrections ###########################################################
# See https://twiki.cern.ch/twiki/bin/view/CMS/WorkBookJetEnergyCorrections
# to check what cff to use
process.load('JetMETCorrections.Configuration.DefaultJEC_cff')
process.jecCorSequence = cms.Sequence(
    process.ak5CaloJetsL2L3*process.ak5PFJetsL2L3
    )
### NB: also check the analysis input below.

### JES MET Corrections ########################################################
from JetMETCorrections.Type1MET.MetType1Corrections_cff import metJESCorAK5CaloJet

process.metMuonJESCorAK5 = metJESCorAK5CaloJet.clone()
process.metMuonJESCorAK5.inputUncorJetsLabel = recoJet_src
#[FR: already default] process.metMuonJESCorAK5.corrector = "ak5CaloL2L3"
process.metMuonJESCorAK5.inputUncorMetLabel = "corMetGlobalMuons"
process.metCorSequence = cms.Sequence(process.metMuonJESCorAK5)

### Analysis configuration #####################################################
process.load("DiLeptonAnalysis.NTupleProducer.ntupleproducer_cfi")
process.analyze.isRealData = cms.untracked.bool(runon=='data')
# Synchronise with Jet configuration above (defaults to ak5)
process.analyze.tag_jets  = recoJet_src
process.analyze.jetCorrs  = 'ak5CaloL2L3'
process.analyze.tag_jetID = 'ak5JetID'
# Add some jet collections
process.analyze.jets = (
    # PF jets
    cms.PSet( prefix = cms.untracked.string('PF'),
              tag = cms.untracked.InputTag('ak5PFJets'),
              sel_minpt  = process.analyze.sel_mincorjpt,
              sel_maxeta = process.analyze.sel_maxjeta,
              corrections = cms.string('ak5PFL2L3'),
              ),
    # Calo jets (for cross-check)
    cms.PSet( prefix = cms.untracked.string('CA'),
              tag = cms.untracked.InputTag('ak5CaloJets'),
              sel_minpt  = process.analyze.sel_mincorjpt,
              sel_maxeta = process.analyze.sel_maxjeta,
              corrections = cms.string('ak5CaloL2L3'),
              ),
    )

#### Path ######################################################################
process.mybtag = cms.Sequence(   process.impactParameterTagInfos
                               * process.simpleSecondaryVertexBJetTags )
process.p = cms.Path(
      process.jecCorSequence
    + process.recoJetIdSequence
    + process.metCorSequence
    + process.mybtag
    + process.analyze
    )
