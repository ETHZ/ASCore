import FWCore.ParameterSet.Config as cms
import HLTrigger.HLTfilters.hltHighLevel_cfi as hlt
import FWCore.ParameterSet.VarParsing as VarParsing

process = cms.Process("NTupleProducer")

### Message Logger #############################################################
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.categories.append('NTP')
process.MessageLogger.cerr.NTP = cms.untracked.PSet(
    limit = cms.untracked.int32(-1),
    reportEvery = cms.untracked.int32(1)
    )
process.MessageLogger.categories.append('EcalSeverityLevelError')
process.MessageLogger.cerr.EcalSeverityLevelError = cms.untracked.PSet(
    limit = cms.untracked.int32(1),
    )
process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(False),
                                      fileMode = cms.untracked.string("NOMERGE")
                                    )

### Parsing of command line parameters #############################################
### (type of run: data, MC; reconstruction: RECO, PAT, PF) #####################
options = VarParsing.VarParsing ('standard') # set 'standard'  options
options.register ('runon', # register 'runon' option
                  'data',  # the default value
                  VarParsing.VarParsing.multiplicity.singleton, # singleton or list
                  VarParsing.VarParsing.varType.string,         # string, int, or float
                  "Type of sample to run on: data (default), MC")
# get and parse the command line arguments
# set NTupleProducer defaults (override the output, files and maxEvents parameter)

#options.files= '/store/data/Commissioning10/MinimumBias/RAW-RECO/v9/000/135/735/FAB17A5D-4465-DF11-8DBF-00E08178C031.root'
#options.files= '/store/mc/Spring10/TTbarJets-madgraph/GEN-SIM-RECO/START3X_V26_S09-v1/0011/A4121AB4-0747-DF11-8984-0030487F171B.root'
options.files= 'file:/scratch/buchmann/82F80C67-E111-E011-8089-003048D2910A.root'
options.maxEvents = -1 # If it is different from -1, string "_numEventXX" will be added to the output file name

# Now parse arguments from command line (might overwrite defaults)
options.parseArguments()
options.output='NTupleProducer_38X_'+options.runon+'.root'

### Running conditions #########################################################
# See https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideFrontierConditions
# to check what cff to use
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
if options.runon=='data':
    # CMSSW_3_8_X:
    process.GlobalTag.globaltag = "GR_R_38X_V14::All"
#     process.GlobalTag.globaltag = "GR10_P_V11::All"
else:
    # CMSSW_3_8_X:
    process.GlobalTag.globaltag = "START38_V13::All"


############ PF2PAT ##########################################
# load the standard PAT config
process.load("PhysicsTools.PatAlgos.patSequences_cff")
from PhysicsTools.PatAlgos.patEventContent_cff import patEventContent

# Add a pro forma output module because PF2PAT complains otherwise...
process.out = cms.OutputModule("PoolOutputModule",
      # save only events passing the full path
      SelectEvents   = cms.untracked.PSet( SelectEvents = cms.vstring('p') ),
      # save PAT Layer 1 output; you need a '*' to
      outputCommands = cms.untracked.vstring('drop *', *patEventContent )
      )

# Configure PAT to use PF2PAT instead of AOD sources
from PhysicsTools.PatAlgos.tools.pfTools import *

postfix = "PF"

if options.runon == 'data':
	usePF2PAT(process,runPF2PAT=True, jetAlgo='AK5', runOnMC=False, postfix=postfix)

else:
	usePF2PAT(process,runPF2PAT=True, jetAlgo='AK5', runOnMC=True, postfix=postfix)

# turn to false when running on data and MC (for the moment)
getattr(process, "patElectrons"+postfix).embedGenMatch = False
getattr(process, "patMuons"+postfix).embedGenMatch = False

process.pfIsolatedMuonsPF.combinedIsolationCut = cms.double(0.15)
process.pfIsolatedElectronsPF.combinedIsolationCut = cms.double(0.15)
process.patJetCorrFactorsPF.levels = cms.vstring()


### b-tagging ##################################################################
# Simple SV and TrackCounting based algos
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
process.load("RecoBTau.JetTagComputer.jetTagRecord_cfi")
process.load("RecoBTag.SecondaryVertex.simpleSecondaryVertex2TrkES_cfi")
process.load("RecoBTag.SecondaryVertex.simpleSecondaryVertex3TrkES_cfi")
process.load("RecoBTag.SecondaryVertex.simpleSecondaryVertexHighEffBJetTags_cfi")
process.load("RecoBTag.SecondaryVertex.simpleSecondaryVertexHighPurBJetTags_cfi")

process.mybtag = cms.Sequence( process.simpleSecondaryVertexHighPurBJetTags
                             * process.simpleSecondaryVertexHighEffBJetTags )

### Input/Output ###############################################################
# Input
process.source = cms.Source("PoolSource",
      fileNames = cms.untracked.vstring(options.files)
      # Enable if you see duplicate error:
      # duplicateCheckMode = cms.untracked.string("noDuplicateCheck")
)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents) )
# Output
process.TFileService = cms.Service("TFileService",
# Keep track of the type of data source and reco type in the ntuple file name
	fileName = cms.string(options.output),
	closeFileFast = cms.untracked.bool(True)
)

#### Electron ID ##############################################################
process.load("DiLeptonAnalysis.NTupleProducer.simpleEleIdSequence_cff")

#### Parameterisation for Jet Corrections and JES ME Corrections ###############
recoJet_src = "ak5CaloJets"
genJet_src = "ak5GenJets"

# Jet ID: add the ones you want
process.load('RecoJets.Configuration.JetIDProducers_cff')
process.recoJetIdSequence = cms.Sequence( process.ak5JetID )

#### Jet Corrections ###########################################################
# See https://twiki.cern.ch/twiki/bin/view/CMS/WorkBookJetEnergyCorrections
# and http://indico.cern.ch/getFile.py/access?contribId=38&sessionId=4&resId=0&materialId=slides&confId=110072 slide 13
# to check what cff to use
process.load('JetMETCorrections.Configuration.DefaultJEC_cff')
process.ak5CaloL2Relative.useCondDB = False
process.ak5CaloL3Absolute.useCondDB = False
process.ak5CaloResidual.useCondDB = False
process.ak5PFL2Relative.useCondDB = False
process.ak5PFL3Absolute.useCondDB = False
process.ak5PFResidual.useCondDB = False
if options.runon=='data':
	process.jecCorSequence = cms.Sequence(
		process.ak5CaloJetsL2L3Residual*process.ak5PFJetsL2L3Residual
	)
else:
	process.jecCorSequence = cms.Sequence(
		process.ak5CaloJetsL2L3*process.ak5PFJetsL2L3
	)
### NB: also check the analysis input below.

### JES MET Corrections ########################################################
from JetMETCorrections.Type1MET.MetType1Corrections_cff import metJESCorAK5CaloJet

process.metMuonJESCorAK5 = metJESCorAK5CaloJet.clone()
process.metMuonJESCorAK5.inputUncorJetsLabel = recoJet_src
process.metMuonJESCorAK5.inputUncorMetLabel = "corMetGlobalMuons"
process.metCorSequence = cms.Sequence(process.metMuonJESCorAK5)

### Cleaning ###################################################################
# flag HB/HE noise
# info: https://twiki.cern.ch/twiki/bin/view/CMS/HcalNoiseInfoLibrary
# uncomment following two lines to filer HCAL noise rather than to flag it!
# process.load('CommonTools/RecoAlgos/HBHENoiseFilter_cfi')
# process.HBHENoiseFilter.maxRBXEMF = cms.double(0.01)
# the next two lines produce the HCAL noise summary flag with an additional cut on RBX
process.load('CommonTools/RecoAlgos/HBHENoiseFilterResultProducer_cfi')
process.HBHENoiseFilterResultProducer.maxRBXEMF = cms.double(0.01)

# ECAL dead cells: this is not a filter. Only a flag is stored.
# Ecal gap boundary energy: specify minimal gap BE here.
### NOTE: THIS PART WAS REDUCED TO ENSURE COMPATIBILITY WITH CMSSW 3_9_7. THIS IS TEMPORARY!
#process.load("PhysicsTools/EcalAnomalousEventFilter/ecalanomalouseventfilter_cfi")
#process.EcalAnomalousEventFilter.FilterAlgo = cms.untracked.string("TuningMode")
#process.EcalAnomalousEventFilter.cutBoundEnergyGapEE = cms.untracked.double(5)
#process.EcalAnomalousEventFilter.cutBoundEnergyGapEB = cms.untracked.double(5)
#process.EcalAnomalousEventFilter.enableGap           = cms.untracked.bool(True)

# See for example DPGAnalysis/Skims/python/MinBiasPDSkim_cfg.py
# require scraping filter
process.scrapingVeto = cms.EDFilter("FilterOutScraping",
     applyfilter = cms.untracked.bool(True),
     debugOn = cms.untracked.bool(False),
     numtrack = cms.untracked.uint32(10),
     thresh = cms.untracked.double(0.25)
     )

# require good primary vertex
process.primaryVertexFilter = cms.EDFilter("GoodVertexFilter",
    vertexCollection = cms.InputTag('offlinePrimaryVertices'),
    minimumNDOF = cms.uint32(4) ,
    maxAbsZ = cms.double(24),
    maxd0 = cms.double(2)
    )
# Cleaning path
# Removed primary vertex cut, can cut offline
# process.cleaning *= process.scrapingVeto*process.primaryVertexFilter


### GenJets ####################################################################
# produce ak5GenJets (collection missing in case of some Spring10 samples)
process.load("RecoJets.Configuration.GenJetParticles_cff")
process.load("RecoJets.Configuration.RecoGenJets_cff")

process.mygenjets = cms.Sequence( process.genParticlesForJets * process.ak5GenJets )

### Analysis configuration #####################################################
process.load("DiLeptonAnalysis.NTupleProducer.ntupleproducer_cfi")
process.analyze.isRealData = cms.untracked.bool(options.runon=='data')
# Synchronise with Jet configuration above (defaults to ak5)
process.analyze.tag_jets  = recoJet_src
process.analyze.jetCorrs  = 'ak5CaloL2L3'
process.analyze.tag_jetID = 'ak5JetID'
# Add some jet collections
process.analyze.jets = (
    # PF jets
    cms.PSet( prefix = cms.untracked.string('PF'),
              tag = cms.untracked.InputTag('ak5PFJets'),
              isPat = cms.untracked.bool(False),
              tag_jetTracks  = cms.untracked.InputTag('ak5JetTracksAssociatorAtVertex'),
              jet_id = cms.untracked.InputTag('ak5JetID'),
              sel_minpt  = cms.double(15.0),
              sel_maxeta = process.analyze.sel_maxjeta,
              corrections = cms.string('ak5PFL2L3'),
              btag_matchdeltaR = cms.double(0.25),
              ),
    # Calo jets (for cross-check)
    cms.PSet( prefix = cms.untracked.string('CA'),
              tag = cms.untracked.InputTag('ak5CaloJets'),
              isPat = cms.untracked.bool(False),
              tag_jetTracks  = cms.untracked.InputTag('ak5JetTracksAssociatorAtVertex'),
              jet_id = cms.untracked.InputTag('ak5JetID'),
              sel_minpt  = process.analyze.sel_mincorjpt,
              sel_maxeta = process.analyze.sel_maxjeta,
              corrections = cms.string('ak5CaloL2L3'),
              btag_matchdeltaR = cms.double(0.25),
              ),
    # PF jets from PF2PAT
    cms.PSet( prefix = cms.untracked.string('PF2PAT'),
              tag = cms.untracked.InputTag('patJetsPF'),
              isPat = cms.untracked.bool(True),
              tag_jetTracks  = cms.untracked.InputTag('ak5JetTracksAssociatorAtVertex'),
              jet_id = cms.untracked.InputTag('ak5JetID'),
              sel_minpt  = cms.double(15.0),
              sel_maxeta = process.analyze.sel_maxjeta,
              corrections = cms.string('ak5PFL2L3'),
              btag_matchdeltaR = cms.double(0.25),
              ),
    )

# Add residual correction for running on data (temporary fix)
if options.runon == 'data':
        process.analyze.jetCorrs = process.analyze.jetCorrs.value() + 'Residual'
        for extJet in process.analyze.jets:
            extJet.corrections = extJet.corrections.value() + 'Residual'

# Dump object information for some HLT trigger filters (from confDB)
process.analyze.hlt_labels = ['hltSingleMu3L3Filtered3',
                              'hltSingleMu5L3Filtered5',
                              'hltSingleMu9L3Filtered9',
                              'hltL1NonIsoHLTNonIsoSingleElectronLWEt10EleIdDphiFilteroHLTNonIsoSingleElectronLWEt10PixelMatchFilter',
                              'hltL1NonIsoHLTNonIsoSingleElectronLWEt10EleIdDphiFilter',
                              'hltL1NonIsoHLTNonIsoSingleElectronLWEt15PixelMatchFilter',
                              'hltL1NonIsoHLTNonIsoSinglePhotonEt10HcalIsolFilter']


#### DEBUG #####################################################################
#process.dump = cms.EDAnalyzer("EventContentAnalyzer")
#process.SimpleMemoryCheck = cms.Service("SimpleMemoryCheck",
#    ignoreTotal = cms.untracked.int32(1) # number of events to ignore at start (default is one)
#)
#process.ProfilerService = cms.Service("ProfilerService",
#                                      firstEvent = cms.untracked.int32(2),
#                                      lastEvent = cms.untracked.int32(51),
#                                      paths = cms.untracked.vstring(['p'])
#                                      )

#### Path ######################################################################
process.p = cms.Path(
    process.scrapingVeto * (
        process.HBHENoiseFilterResultProducer
        # + process.EcalAnomalousEventFilter
        + process.mygenjets
        + process.jecCorSequence
        + process.recoJetIdSequence
        + process.simpleEleIdSequence
        + process.metCorSequence
        + process.mybtag
        + process.patPF2PATSequencePF
        + process.analyze
        )
    )

# remove output path from patTemplate_cfg
# process.p.remove(process.out)


# remove ak5GenJets from the path if it will run on data
if options.runon=='data':
    process.p.remove(process.mygenjets)
if options.runon!='data':
    process.p.remove(process.scrapingVeto)

