import FWCore.ParameterSet.Config as cms
import HLTrigger.HLTfilters.hltHighLevel_cfi as hlt
import FWCore.ParameterSet.VarParsing as VarParsing

process = cms.Process("NTupleProducer")

##############PAT part#######################
# import skeleton process
from PhysicsTools.PatAlgos.patTemplate_cfg import *

from PhysicsTools.PatAlgos.tools.cmsswVersionTools import *
## uncomment this line to run on an 35X input sample
run36xOn35xInput(process)

# load the standard PAT config 
process.load("PhysicsTools.PatAlgos.patSequences_cff")

from PhysicsTools.PatAlgos.tools.coreTools import *
removeMCMatching(process, ['All'])

from PhysicsTools.PatAlgos.tools.jetTools import *
process.load('JetMETCorrections.Configuration.DefaultJEC_cff')

# uncomment the following line to add tcMET to the event content
from PhysicsTools.PatAlgos.tools.metTools import *
addTcMET(process, 'TC')
addPfMET(process, 'PF')

# in PAT (iterativeCone5) to ak5 (anti-kt cone = 0.5)
switchJetCollection35X(process, 
                 cms.InputTag('ak5CaloJets'),   
                 doJTA            = True,            
                 doBTagging       = True,            
                 jetCorrLabel     = ('AK5','Calo'),  
                 doType1MET       = True,
                 genJetCollection = cms.InputTag("ak5GenJets"),
                 doJetID          = True,
                 jetIdLabel       = "ak5"
                        )


# uncomment the following lines to add antikt5JPT jets to your PAT output
addJetCollection35X(process,cms.InputTag('JetPlusTrackZSPCorJetAntiKt5'),
                 'AK5','JPT',
                 doJTA        = True,
                 doBTagging   = True,
                 jetCorrLabel = ('AK5','JPT'),
                 doType1MET   = False,
                 doL1Cleaning = False,
                 genJetCollection = cms.InputTag("ak5GenJets"),
                 doJetID      = True,
		 jetIdLabel       = "ak5"
                 )

# uncomment the following lines to add antikt5Pflow jets to your PAT output
addJetCollection35X(process,cms.InputTag('ak5PFJets'),
                 'AK5','PF',
                 doJTA        = True,
                 doBTagging   = True,
                 jetCorrLabel = ('AK5','PF'),
                 doType1MET   = False,
                 doL1Cleaning = False,
                 genJetCollection = cms.InputTag("ak5GenJets"),
                 doJetID      = True,
		 jetIdLabel       = "ak5"
                 )




### Message Logger #############################################################
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.categories.append('NTP')
process.MessageLogger.cerr.NTP = cms.untracked.PSet(
    limit = cms.untracked.int32(-1),
    reportEvery = cms.untracked.int32(1)
    )
process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(False) )

### Parsing of command line parameters #############################################
### (type of run: data, MC; reconstruction: RECO, PAT, PF) #####################
options = VarParsing.VarParsing ('standard') # set 'standard'  options
options.register ('runon', # register 'runon' option
                  'MC', # the default value
                  VarParsing.VarParsing.multiplicity.singleton, # singleton or list
                  VarParsing.VarParsing.varType.string,         # string, int, or float
                  "Type of sample to run on: data (default), MC")
options.register ('recoType',									# register 'recoType' option
                  'PAT',										# the default value
                  VarParsing.VarParsing.multiplicity.singleton, # singleton or list
                  VarParsing.VarParsing.varType.string,         # string, int, or float
                  "Type of reconstruction to use: RECO (default), PAT, PF")

# get and parse the command line arguments
# set NTupleProducer defaults (override the output, files and maxEvents parameter)

options.files= 'file:/hadoop/Other/35XSpring10/TTbar/GEN-SIM-RECO/START3X_V26_S09-v1/0093/FAA0BC3B-924E-DF11-AFAF-0017A4770C2C.root'
#options.files= 'file:/hadoop/phedex/oviedo/cms/store/data/Commissioning10/MinimumBias/RECO/May6thPDSkim2_SD_Mu-v1/0137/F62CD0CB-D95D-DF11-A830-00261894389E.root'
#options.files='file:/hadoop/Other/TTbarJets-madgraph_20KEventsSincro/reco_7TeV_2_1.root','file:/hadoop/Other/TTbarJets-madgraph_20KEventsSincro/reco_7TeV_1_1.root'
#options.files = 'file:/pool/fanae20/Zmumu_Summer10.root'

options.maxEvents = 10 # If it is different from -1, string "_numEventXX" will be added to the output file name

# Now parse arguments from command line (might overwrite defaults)
options.parseArguments()
options.output='NTupleProducer_36X_'+options.runon+'_'+options.recoType+'.root'

### Running conditions #########################################################
# See https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideFrontierConditions
# to check what cff to use
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
if options.runon=='data':
    # CMSSW_3_6_X:
    process.GlobalTag.globaltag = "GR_R_36X_V12::All"
else:
    # CMSSW_3_6_X:
    process.GlobalTag.globaltag = "START36_V10::All"

### b-tagging ##################################################################
# Simple SV and TrackCounting based algos
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
process.load("RecoBTau.JetTagComputer.jetTagRecord_cfi")
process.load("RecoBTag.SecondaryVertex.simpleSecondaryVertex2TrkES_cfi")
process.load("RecoBTag.SecondaryVertex.simpleSecondaryVertex3TrkES_cfi")
process.load("RecoBTag.SecondaryVertex.simpleSecondaryVertexHighEffBJetTags_cfi")
process.load("RecoBTag.SecondaryVertex.simpleSecondaryVertexHighPurBJetTags_cfi")

### Input/Output ###############################################################
# Input
process.source = cms.Source("PoolSource",
      fileNames = cms.untracked.vstring(options.files)
#Enable if you see duplicate error      duplicateCheckMode = cms.untracked.string("noDuplicateCheck")
)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents) )
# Output
process.TFileService = cms.Service("TFileService",
# keep track of the type of data source and reco type in the ntuple file name
	fileName = cms.string(options.output),  
	closeFileFast = cms.untracked.bool(True)
)


#### Produce JPT jets #########################################################
process.load('RecoJets.Configuration.RecoJPTJets_cff')
					
					
#### Parameterisation for Jet Corrections and JES ME Corrections ###############
recoJet_src = "ak5CaloJets"
genJet_src = "ak5GenJets"

# Jet ID: add the ones you want
#process.load('RecoJets.Configuration.JetIDProducers_cff')
#process.recoJetIdSequence = cms.Sequence( process.ak5JetID )

#### Jet Corrections ###########################################################
# See https://twiki.cern.ch/twiki/bin/view/CMS/WorkBookJetEnergyCorrections
# to check what cff to use
#process.load('JetMETCorrections.Configuration.DefaultJEC_cff')
#process.jecCorSequence = cms.Sequence(
#    process.ak5CaloJetsL2L3*process.ak5PFJetsL2L3*process.ak5JPTJetsL2L3
#    )
### NB: also check the analysis input below.		

### JES MET Corrections ########################################################
from JetMETCorrections.Type1MET.MetType1Corrections_cff import metJESCorAK5CaloJet

process.metMuonJESCorAK5 = metJESCorAK5CaloJet.clone()
process.metMuonJESCorAK5.inputUncorJetsLabel = recoJet_src
#[FR: already default] process.metMuonJESCorAK5.corrector = "ak5CaloL2L3"
process.metMuonJESCorAK5.inputUncorMetLabel = "corMetGlobalMuons"
process.metCorSequence = cms.Sequence(process.metMuonJESCorAK5)

### Cleaning ###################################################################
# flag HB/HE noise
process.load('CommonTools/RecoAlgos/HBHENoiseFilterResultProducer_cfi')
process.cleaning = cms.Sequence(process.HBHENoiseFilterResultProducer)

# See for example DPGAnalysis/Skims/python/MinBiasPDSkim_cfg.py
if options.runon=='data':
    # require physics declared
    process.load('HLTrigger.special.hltPhysicsDeclared_cfi')
    process.hltPhysicsDeclared.L1GtReadoutRecordTag = 'gtDigis'

    # require scraping filter
    process.scrapingVeto = cms.EDFilter("FilterOutScraping",
                                        applyfilter = cms.untracked.bool(True),
                                        debugOn = cms.untracked.bool(False),
                                        numtrack = cms.untracked.uint32(10),
                                        thresh = cms.untracked.double(0.2)
                                        )

    # require good primary vertex
    process.primaryVertexFilter = cms.EDFilter("GoodVertexFilter",
                                               vertexCollection = cms.InputTag('offlinePrimaryVertices'),
                                               minimumNDOF = cms.uint32(4) ,
                                               maxAbsZ = cms.double(15),
                                               maxd0 = cms.double(2)
                                               )
    # Cleaning path
    process.cleaning *= process.scrapingVeto*process.hltPhysicsDeclared*process.primaryVertexFilter
        

### GenJets ####################################################################
# produce ak5GenJets (collection missing in case of some Spring10 samples)
process.load("RecoJets.Configuration.GenJetParticles_cff")
process.load("RecoJets.Configuration.RecoGenJets_cff")	

process.mygenjets = cms.Sequence( process.genParticlesForJets
								* process.ak5GenJets )

### Analysis configuration #####################################################
process.load("DiLeptonAnalysis.NTupleProducer.ntuplepatproducer_cfi")
process.analyze.isRealData = cms.untracked.bool(options.runon=='data')
# Synchronise with Jet configuration above (defaults to ak5)
process.analyze.tag_jets  = recoJet_src
process.analyze.jetCorrs  = 'ak5CaloL2L3'
process.analyze.tag_jetID = 'ak5JetID'
# Add some jet collections
process.analyze.jets = (
    # PF jets
    cms.PSet( prefix = cms.untracked.string('PF'),
              tag = cms.untracked.InputTag('patJetsAK5PF::PAT'),
	      tag_jetTracks  = cms.untracked.InputTag('ak5JetTracksAssociatorAtVertex'),
	      jet_id = cms.untracked.InputTag('ak5JetID'),			  
              sel_minpt  = process.analyze.sel_mincorjpt,
              sel_maxeta = process.analyze.sel_maxjeta,
              corrections = cms.string('ak5PFL2L3'),
              ),
	# JPT
	cms.PSet( prefix = cms.untracked.string('JPT'),
              tag = cms.untracked.InputTag('patJetsAK5JPT::PAT'),
  	      tag_jetTracks  = cms.untracked.InputTag('ak5JetTracksAssociatorAtVertex'),
	      jet_id = cms.untracked.InputTag('ak5JetID'),
              sel_minpt  = process.analyze.sel_mincorjpt,
              sel_maxeta = process.analyze.sel_maxjeta,
              corrections = cms.string('ak5JPTL2L3'),
              ),		  		  
    # Calo jets (for cross-check)
    cms.PSet( prefix = cms.untracked.string('CA'),
              tag = cms.untracked.InputTag('patJets::PAT'),
	      tag_jetTracks  = cms.untracked.InputTag('ak5JetTracksAssociatorAtVertex'),
              jet_id = cms.untracked.InputTag('ak5JetID'),
              sel_minpt  = process.analyze.sel_mincorjpt,
              sel_maxeta = process.analyze.sel_maxjeta,
              corrections = cms.string('ak5CaloL2L3'),
              ),
    )

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
process.mybtag = cms.Sequence(	process.simpleSecondaryVertexHighPurBJetTags
								* process.simpleSecondaryVertexHighEffBJetTags )

process.p = cms.Path(
    process.recoJPTJets  
    + process.mygenjets
    + process.cleaning
    + process.patDefaultSequence
#    + process.jecCorSequence
#    + process.recoJetIdSequence
#    + process.metCorSequence
#    + process.mybtag
	#	+ process.dump
    + process.analyze
    )

# remove ak5GenJets from the path if it will run on data
if options.runon=='data':
	process.p.remove(process.mygenjets)

