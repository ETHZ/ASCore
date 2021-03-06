import FWCore.ParameterSet.Config as cms
import HLTrigger.HLTfilters.hltHighLevel_cfi as hlt
import FWCore.ParameterSet.VarParsing as VarParsing

process = cms.Process("NTupleProducer")
process.load("SimGeneral.HepPDTESSource.pdt_cfi")

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
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(False))
process.MessageLogger.suppressWarning.append('TriggerResultsFilter')
process.MessageLogger.suppressInfo.append('TriggerResultsFilter')

### Parsing of command line parameters #############################################
### (type of run: data, MC; reconstruction: RECO, PAT, PF) #####################
options = VarParsing.VarParsing ('standard') # set 'standard'  options
options.register ('runon', # register 'runon' option
                  'data',  # the default value
                  VarParsing.VarParsing.multiplicity.singleton, # singleton or list
                  VarParsing.VarParsing.varType.string,         # string, int, or float
                  "Type of sample to run on: data (default), MC")
options.register ('ModelScan', # register 'runon' option
                  False,  # the default value
                  VarParsing.VarParsing.multiplicity.singleton, # singleton or list
                  VarParsing.VarParsing.varType.bool,         # string, int, or float
                  "If you are dealing with a model scan, set this to True, otherwise to False (default)")
options.register ('FastSim', # register 'runon' option
                  False,  # the default value
                  VarParsing.VarParsing.multiplicity.singleton, # singleton or list
                  VarParsing.VarParsing.varType.bool,         # string, int, or float
                  "If you are dealing with a FastSim (but not a model scan!), set this to True, otherwise to False (default)")
options.register ('doPhotonStuff',
                  False,
                  VarParsing.VarParsing.multiplicity.singleton, # singleton or list
                  VarParsing.VarParsing.varType.bool,         # string, int, or float
                  "If you want to run the photon-analyses specific configuration, set to True, otherwise False (default)")
options.register ('is2011',
                  False,
                  VarParsing.VarParsing.multiplicity.singleton, # singleton or list
                  VarParsing.VarParsing.varType.bool,         # string, int, or float
                  "If you are running on the 2011 data or MC rerecoes, set to True, otherwise False (default)")
options.register ('GlobalTag',
                  '',
                  VarParsing.VarParsing.multiplicity.singleton, # singleton or list
                  VarParsing.VarParsing.varType.string,         # string, int, or float
                  "Global tag (if unspecified, use config. defaults)")
# get and parse the command line arguments
# set NTupleProducer defaults (override the output, files and maxEvents parameter)
#options.files= 'file:////shome/mdunser/files/isoSynchFile_DoubleMu191700.root'
#options.files= 'file://///shome/pablom/tmp/newCode/CMSSW_5_2_5_patch1/src/DiLeptonAnalysis/NTupleProducer/A8922572-9D84-E111-88B9-003048F024FE.root'
#options.files= 'file:////scratch/mdunser/files/DoubleMu_Run2012D.root'
#options.files= 'file:////shome/casal/files/trackingPOGevents.root'
options.files= 'file:////scratch/mdunser/files/pickevents.root'
#options.files= 'file:////shome/mdunser/files/JetHT_Run2012C_v1.root'
#options.files='file:////scratch/fronga/RelValTTbarLepton_EE4E6727-2C7A-E111-A4E8-002354EF3BCE.root'

options.maxEvents = -1# If it is different from -1, string "_numEventXX" will be added to the output file name 
# Now parse arguments from command line (might overwrite defaults)
options.parseArguments()
if (options.runon!="data"):
   options.runon="MC"
options.output='NTupleProducer_53X_'+options.runon+'.root'

### Running conditions #########################################################
# See https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideFrontierConditions
# to check what cff to use
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

# try this in the future instead of a global tag. still gives some errors at the moment (apr17)
#from Configuration.AlCa.autoCond import autoCond
#process.GlobalTag.globaltag = cms.string( autoCond[ 'startup' ] )
if options.GlobalTag!='':
   process.GlobalTag.globaltag = options.GlobalTag
else:
   if options.runon=='data':
      if (options.is2011):
         process.GlobalTag.globaltag = "FT_R_53_LV3::All" 
      else:
         process.GlobalTag.globaltag = "FT53_V21A_AN6::All"
   else:
      if (options.is2011):
         process.GlobalTag.globaltag = "START53_LV4::All"
      else:
         process.GlobalTag.globaltag = "START53_V7N::All"


### Input/Output ###############################################################
# Input
process.source = cms.Source("PoolSource",
      fileNames = cms.untracked.vstring(options.files)
      # Enable if you see duplicate error:
      # duplicateCheckMode = cms.untracked.string("noDuplicateCheck")
)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents) )
process.ModelScan = cms.untracked.PSet( input = cms.untracked.bool(options.ModelScan) )
process.FastSim = cms.untracked.PSet( input = cms.untracked.bool(options.FastSim) )

# Output
process.out = cms.OutputModule("PoolOutputModule",
                               fileName = cms.untracked.string(options.output),
                               #SelectEvents   = cms.untracked.PSet( SelectEvents = cms.vstring('p') ),
                               splitLevel = cms.untracked.int32(99),  # Turn on split level (smaller files???)
                               dropMetaData = cms.untracked.string('ALL'), # Get rid of metadata related to dropped collections
                               outputCommands = cms.untracked.vstring() # Will be overwritten by PAT: we overwrite at the end
                               )

    
### Jet/MET Corrections ##########################################################
# See https://twiki.cern.ch/twiki/bin/view/CMS/WorkBookJetEnergyCorrections
# Charged hadron subtraction (put that first, it loads many other things)
process.load("PhysicsTools.PatAlgos.patSequences_cff")
from PhysicsTools.PatAlgos.tools.pfTools import *
pfPostfix= 'PFCHS'
# ---  CHS enabled ---
usePF2PAT(process,runPF2PAT=True, jetAlgo='AK5', runOnMC=(options.runon!='data'), postfix=pfPostfix,
          jetCorrections=('AK5PFchs', ['L1FastJet', 'L2Relative', 'L3Absolute']),  # residual set for data later on
          pvCollection=cms.InputTag('goodVertices'))
process.pfPileUpPFCHS.checkClosestZVertex = False
process.pfNoPileUpPFCHS.enable = True ## pfNoPU enabled
# ---  CHS diabled ---
# usePF2PAT(process,runPF2PAT=True, jetAlgo='AK5', runOnMC=(options.runon!='data'), postfix=pfPostfix,
#           jetCorrections=('AK5PF', ['L1FastJet', 'L2Relative', 'L3Absolute']),  # residual set for data later on
#           pvCollection=cms.InputTag('goodVertices'))
# process.pfPileUpPFCHS.checkClosestZVertex = False
# process.pfNoPileUpPFCHS.enable = False ## pfNoPU enabled
# --------- disable top-projections: do old style cleaning --------------------
process.pfNoElectronPFCHS.enable = False 
process.pfNoMuonPFCHS.enable = False
process.pfNoTauPFCHS.enable = False

# L1 Fast-jet corrections (only PF jets)
process.load('JetMETCorrections.Configuration.DefaultJEC_cff')
process.load('RecoJets.Configuration.RecoPFJets_cff')

process.kt6PFJetsForIso = process.kt4PFJets.clone( rParam = 0.6, doRhoFastjet = True )
process.kt6PFJetsForIso.Rho_EtaMax = cms.double(2.5)

# MET corrections
from JetMETCorrections.Type1MET.pfMETCorrections_cff import pfJetMETcorr,pfType1CorrectedMet
process.pfJetMETcorr = pfJetMETcorr.clone()
process.pfType1CorrectedMet = pfType1CorrectedMet.clone()
if options.runon == 'data':
    process.pfJetMETcorr.jetCorrLabel = process.pfJetMETcorr.jetCorrLabel.value()+'Residual'


from JetMETCorrections.Type1MET.MetType1Corrections_cff import metJESCorAK5CaloJet
process.metMuonJESCorAK5 = metJESCorAK5CaloJet.clone()
process.metMuonJESCorAK5.inputUncorJetsLabel = "ak5CaloJets"
process.metMuonJESCorAK5.inputUncorMetLabel = "corMetGlobalMuons"

process.metCorSequence = cms.Sequence(process.pfJetMETcorr + process.pfType1CorrectedMet + process.metMuonJESCorAK5)


######### Get a list of good primary vertices, in 42x, these are DAF vertices ################
process.goodVertices = cms.EDFilter("VertexSelector",
	src = cms.InputTag("offlinePrimaryVertices"),
	cut = cms.string("!isFake && ndof > 4 && abs(z) <= 24 && position.Rho <= 2"),
	filter = cms.bool(True)
)

### Cleaning ###################################################################
### See https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFilters

## The iso-based HBHE noise filter ___________________________________________||
process.load('CommonTools.RecoAlgos.HBHENoiseFilterResultProducer_cfi')

## The HCAL laser filter _____________________________________________________||
process.load("RecoMET.METFilters.hcalLaserEventFilter_cfi")
process.hcalLaserEventFilter.vetoByRunEventNumber=False
process.hcalLaserEventFilter.vetoByHBHEOccupancy=True
process.hcalLaserEventFilter.taggingMode=True

# new recommendation: https://twiki.cern.ch/twiki/bin/view/CMS/PdmVKnowFeatures#HCAL_laser_events_in_prompt_2012
process.load("EventFilter.HcalRawToDigi.hcallasereventfilter2012_cfi")

## The ECAL dead cell trigger primitive filter _______________________________||
process.load('RecoMET.METFilters.EcalDeadCellTriggerPrimitiveFilter_cfi')
## For AOD and RECO recommendation to use recovered rechits
process.EcalDeadCellTriggerPrimitiveFilter.taggingMode = True

## The EE bad SuperCrystal filter ____________________________________________||
process.load('RecoMET.METFilters.eeBadScFilter_cfi')
process.eeBadScFilter.taggingMode = True

## The tracking failure filter _______________________________________________||
process.load('RecoMET.METFilters.trackingFailureFilter_cfi')
process.trackingFailureFilter.taggingMode = True

## Laser correction filter (Run2012 A+B, Jul13 ReReco)
process.load('RecoMET.METFilters.ecalLaserCorrFilter_cfi')
process.ecalLaserCorrFilter.taggingMode = True

## The tracking POG filters __________________________________________________||
process.load('RecoMET.METFilters.trackingPOGFilters_cff')
process.manystripclus53X.taggedMode = cms.untracked.bool(True)
process.manystripclus53X.forcedValue = cms.untracked.bool(False)
process.toomanystripclus53X.taggedMode = cms.untracked.bool(True)
process.toomanystripclus53X.forcedValue = cms.untracked.bool(False)
process.logErrorTooManyClusters.taggedMode = cms.untracked.bool(True)
process.logErrorTooManyClusters.forcedValue = cms.untracked.bool(False)

### GenJets ####################################################################
# produce ak5GenJets (collection missing in case of some Spring10 samples)
process.load("RecoJets.Configuration.GenJetParticles_cff")
process.load("RecoJets.Configuration.RecoGenJets_cff")
process.mygenjets = cms.Sequence( process.genParticlesForJets * process.ak5GenJets )

### Analysis configuration #####################################################
process.load("DiLeptonAnalysis.NTupleProducer.ntupleproducer_cfi")
process.analyze.isRealData      = (options.runon=='data')
process.analyze.isModelScan     = options.ModelScan
process.analyze.isFastSim     = options.FastSim
process.analyze.tag_doPhotonStuff = options.doPhotonStuff
process.analyze.tag_vertex = 'goodVertices'
if options.doPhotonStuff==True:
   process.analyze.tag_vertex = process.analyze.tag_vertex_withbs
if (options.is2011):
   process.analyze.tag_regressionVersion = cms.int32(8)
else:
   process.analyze.tag_regressionVersion = cms.int32(5)

# Add some jet collections
glist_btags = cms.vstring('trackCountingHighEffBJetTags','trackCountingHighPurBJetTags',
		          'simpleSecondaryVertexHighEffBJetTags','simpleSecondaryVertexHighPurBJetTags',
			  'combinedSecondaryVertexBJetTags', 'combinedSecondaryVertexMVABJetTags',
			  'jetBProbabilityBJetTags','jetProbabilityBJetTags')

process.analyze.jets = (
   # Calo jets
     cms.PSet( prefix = cms.string('CA'),
               tag = cms.InputTag('ak5CaloJets'),
               isPat = cms.bool(False),
               tag_jetTracks  = cms.InputTag('ak5JetTracksAssociatorAtVertex'),
               jet_id = cms.InputTag('ak5JetID'),
               sel_minpt  = process.analyze.sel_mincorjpt,
               sel_maxeta = process.analyze.sel_maxjeta,
               corrections = cms.string('ak5CaloL2L3'),
               ), 
     # pf jets with CHS
     cms.PSet( prefix = cms.string('PFCHS'),
               tag = cms.InputTag('patJetsPFCHS'),
               isPat = cms.bool(True),
               sel_minpt  = process.analyze.sel_mincorjpt,
               sel_maxeta = process.analyze.sel_maxjeta,
	       list_btags = glist_btags,
               ), 
     )
# Add PF candidates
process.analyze.pfCandidates = (
     cms.PSet( prefix = cms.string('PFC'),
     tag = cms.InputTag('particleFlow'),
     sel_minpt = cms.double(5.0),
     sel_maxeta = cms.double(5.0), # Not actually used
     ),
)

# Add taus
process.analyze.leptons = (
    cms.PSet( type = cms.string('tau'),
              prefix = cms.string('Tau'),
              tag = cms.InputTag('selectedNewTaus'),
              sel_minpt = process.analyze.sel_minelpt,
              sel_maxeta = process.analyze.sel_maxeleta,
              maxnobjs = cms.uint32(20)
              ),)
    
# # Add residual correction for running on data
# # taken from local sqlite file. see: https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookJetEnergyCorrections#JetEnCor2011V2 
if options.runon == 'data':
    process.analyze.jetCorrs = process.analyze.jetCorrs.value() + 'Residual'
    for extJet in process.analyze.jets:
        if hasattr(extJet,'corrections'): extJet.corrections = extJet.corrections.value() + 'Residual'
    process.patJetCorrFactorsPFCHS.levels.extend( ['L2L3Residual'] )

# Corrected jet collection
if options.runon == 'data':
   process.ak5PFJetsCorrected = process.ak5PFJetsL1FastL2L3Residual.clone()
else:
   process.ak5PFJetsCorrected = process.ak5PFJetsL1FastL2L3.clone()
              

#### Steve Mrenna's Photon - Parton DR match #######################	      
process.printGenParticles = cms.EDAnalyzer("ParticleListDrawer",
	src = cms.InputTag("partonGenJets"),
	maxEventsToPrint = cms.int32(10)
)
#
process.printPhotons = cms.EDAnalyzer("ParticleListDrawer",
     src = cms.InputTag("myphotons"),
     maxEventsToPrint = cms.untracked.int32(10)
)
#
process.printPartons = cms.EDAnalyzer("ParticleListDrawer",
     src = cms.InputTag("myPhotonJetMatch"),
     maxEventsToPrint = cms.untracked.int32(10)
)
#
process.load('DiLeptonAnalysis.NTupleProducer.photonPartonMatch_cfi')

#### DEBUG TOOLS ###################################################
#process.dump = cms.EDAnalyzer("EventContentAnalyzer")
#process.SimpleMemoryCheck = cms.Service("SimpleMemoryCheck",
##                                        ignoreTotal = cms.untracked.int32(1) # number of events to ignore at start (default is one)
#                                        )
#process.ProfilerService = cms.Service("ProfilerService",
#                                      firstEvent = cms.untracked.int32(2),
#                                      lastEvent = cms.untracked.int32(51),
#                                      paths = cms.untracked.vstring(['p'])
#                                      )
#process.Tracer = cms.Service("Tracer")


###############################################################################
### B-tagging general configuration ###########################################
### Need to re-do b-tagging when not using PF2PAT
process.load("RecoBTag.Configuration.RecoBTag_cff")

from RecoJets.JetAssociationProducers.ak5JTA_cff import *
from RecoBTag.Configuration.RecoBTag_cff import *

# create a new jets and tracks association
import RecoJets.JetAssociationProducers.ak5JTA_cff
process.newJetTracksAssociatorAtVertex = RecoJets.JetAssociationProducers.ak5JTA_cff.ak5JetTracksAssociatorAtVertex.clone()
process.newJetTracksAssociatorAtVertex.jets   = 'ak5PFJets'
process.newJetTracksAssociatorAtVertex.tracks = 'generalTracks'

process.newImpactParameterTagInfos = RecoBTag.Configuration.RecoBTag_cff.impactParameterTagInfos.clone()
process.newImpactParameterTagInfos.jetTracks = "newJetTracksAssociatorAtVertex"
process.newTrackCountingHighEffBJetTags          = RecoBTag.Configuration.RecoBTag_cff.trackCountingHighEffBJetTags.clone()
process.newTrackCountingHighEffBJetTags.tagInfos = cms.VInputTag( cms.InputTag("newImpactParameterTagInfos") )
process.newTrackCountingHighPurBJetTags         = RecoBTag.Configuration.RecoBTag_cff.trackCountingHighPurBJetTags.clone()
process.newTrackCountingHighPurBJetTags.tagInfos= cms.VInputTag( cms.InputTag("newImpactParameterTagInfos") )
process.newJetProbabilityBPFJetTags          = RecoBTag.Configuration.RecoBTag_cff.jetProbabilityBJetTags.clone()
process.newJetProbabilityBPFJetTags.tagInfos = cms.VInputTag( cms.InputTag("newImpactParameterTagInfos") )
process.newJetBProbabilityBPFJetTags          = RecoBTag.Configuration.RecoBTag_cff.jetBProbabilityBJetTags.clone()
process.newJetBProbabilityBPFJetTags.tagInfos = cms.VInputTag( cms.InputTag("newImpactParameterTagInfos") )

process.newSecondaryVertexTagInfos                 = RecoBTag.Configuration.RecoBTag_cff.secondaryVertexTagInfos.clone()
process.newSecondaryVertexTagInfos.trackIPTagInfos = "newImpactParameterTagInfos"
process.newSimpleSecondaryVertexHighEffBJetTags          = RecoBTag.Configuration.RecoBTag_cff.simpleSecondaryVertexHighEffBJetTags.clone()
process.newSimpleSecondaryVertexHighEffBJetTags.tagInfos = cms.VInputTag( cms.InputTag("newSecondaryVertexTagInfos") )
process.newSimpleSecondaryVertexHighPurBJetTags         = RecoBTag.Configuration.RecoBTag_cff.simpleSecondaryVertexHighPurBJetTags.clone()
process.newSimpleSecondaryVertexHighPurBJetTags.tagInfos= cms.VInputTag( cms.InputTag("newSecondaryVertexTagInfos") )
process.newCombinedSecondaryVertexBPFJetTags          = RecoBTag.Configuration.RecoBTag_cff.combinedSecondaryVertexBJetTags.clone()
process.newCombinedSecondaryVertexBPFJetTags.tagInfos = cms.VInputTag( cms.InputTag("newImpactParameterTagInfos"), cms.InputTag("newSecondaryVertexTagInfos") )
process.newCombinedSecondaryVertexMVABPFJetTags          = RecoBTag.Configuration.RecoBTag_cff.combinedSecondaryVertexMVABJetTags.clone()
process.newCombinedSecondaryVertexMVABPFJetTags.tagInfos = cms.VInputTag( cms.InputTag("newImpactParameterTagInfos"), cms.InputTag("newSecondaryVertexTagInfos") )

process.newJetTracksAssociator = cms.Sequence( process.newJetTracksAssociatorAtVertex )

process.newJetBtaggingIP = cms.Sequence(
    process.newImpactParameterTagInfos * (
       process.newTrackCountingHighEffBJetTags +
       process.newTrackCountingHighPurBJetTags +
       process.newJetProbabilityBPFJetTags +
       process.newJetBProbabilityBPFJetTags )
    )

process.newJetBtaggingSV = cms.Sequence(
    process.newImpactParameterTagInfos *
    process.newSecondaryVertexTagInfos * (
       process.newSimpleSecondaryVertexHighEffBJetTags +
       process.newSimpleSecondaryVertexHighPurBJetTags +
       process.newCombinedSecondaryVertexBPFJetTags +
       process.newCombinedSecondaryVertexMVABPFJetTags )
    )

process.newJetBtagging = cms.Sequence(
    process.newJetBtaggingIP +
    process.newJetBtaggingSV )

process.newBtaggingSequence = cms.Sequence( process.newJetTracksAssociator * process.newJetBtagging )



#same thing for PF:
import RecoJets.JetAssociationProducers.ak5JTA_cff
process.newPFJetTracksAssociatorAtVertex        = RecoJets.JetAssociationProducers.ak5JTA_cff.ak5JetTracksAssociatorAtVertex.clone()
process.newPFJetTracksAssociatorAtVertex.jets   = 'ak5PFJets'
process.newPFJetTracksAssociatorAtVertex.tracks = 'generalTracks'

process.newPFImpactParameterTagInfos = RecoBTag.Configuration.RecoBTag_cff.impactParameterTagInfos.clone()
process.newPFImpactParameterTagInfos.jetTracks = "newPFJetTracksAssociatorAtVertex"
process.newPFTrackCountingHighEffBJetTags = RecoBTag.Configuration.RecoBTag_cff.trackCountingHighEffBJetTags.clone()
process.newPFTrackCountingHighEffBJetTags.tagInfos = cms.VInputTag( cms.InputTag("newPFImpactParameterTagInfos") )
process.newPFTrackCountingHighPurBJetTags = RecoBTag.Configuration.RecoBTag_cff.trackCountingHighPurBJetTags.clone()
process.newPFTrackCountingHighPurBJetTags.tagInfos = cms.VInputTag( cms.InputTag("newPFImpactParameterTagInfos") )
process.newPFJetProbabilityBPFJetTags = RecoBTag.Configuration.RecoBTag_cff.jetProbabilityBJetTags.clone()
process.newPFJetProbabilityBPFJetTags.tagInfos = cms.VInputTag( cms.InputTag("newPFImpactParameterTagInfos") )
process.newPFJetBProbabilityBPFJetTags = RecoBTag.Configuration.RecoBTag_cff.jetBProbabilityBJetTags.clone()
process.newPFJetBProbabilityBPFJetTags.tagInfos = cms.VInputTag( cms.InputTag("newPFImpactParameterTagInfos") )


process.newPFSecondaryVertexTagInfos = RecoBTag.Configuration.RecoBTag_cff.secondaryVertexTagInfos.clone()
process.newPFSecondaryVertexTagInfos.trackIPTagInfos = "newPFImpactParameterTagInfos"
process.newPFSimpleSecondaryVertexHighEffBJetTags = RecoBTag.Configuration.RecoBTag_cff.simpleSecondaryVertexHighEffBJetTags.clone()
process.newPFSimpleSecondaryVertexHighEffBJetTags.tagInfos = cms.VInputTag( cms.InputTag("newPFSecondaryVertexTagInfos") )
process.newPFSimpleSecondaryVertexHighPurBJetTags = RecoBTag.Configuration.RecoBTag_cff.simpleSecondaryVertexHighPurBJetTags.clone()
process.newPFSimpleSecondaryVertexHighPurBJetTags.tagInfos = cms.VInputTag( cms.InputTag("newPFSecondaryVertexTagInfos") )
process.newPFCombinedSecondaryVertexBPFJetTags = RecoBTag.Configuration.RecoBTag_cff.combinedSecondaryVertexBJetTags.clone()
process.newPFCombinedSecondaryVertexBPFJetTags.tagInfos = cms.VInputTag( cms.InputTag("newPFImpactParameterTagInfos"), cms.InputTag("newPFSecondaryVertexTagInfos") )
process.newPFCombinedSecondaryVertexMVABPFJetTags = RecoBTag.Configuration.RecoBTag_cff.combinedSecondaryVertexMVABJetTags.clone()
process.newPFCombinedSecondaryVertexMVABPFJetTags.tagInfos = cms.VInputTag( cms.InputTag("newPFImpactParameterTagInfos"), cms.InputTag("newPFSecondaryVertexTagInfos") )


process.newPFJetTracksAssociator = cms.Sequence( process.newPFJetTracksAssociatorAtVertex )

process.newPFJetBtaggingIP = cms.Sequence(
    process.newPFImpactParameterTagInfos * (
       process.newPFTrackCountingHighEffBJetTags +
       process.newPFTrackCountingHighPurBJetTags +
       process.newPFJetProbabilityBPFJetTags +
       process.newPFJetBProbabilityBPFJetTags )
    )

process.newPFJetBtaggingSV = cms.Sequence(
    process.newPFImpactParameterTagInfos *
    process.newPFSecondaryVertexTagInfos * (
       process.newPFSimpleSecondaryVertexHighEffBJetTags +
       process.newPFSimpleSecondaryVertexHighPurBJetTags +
       process.newPFCombinedSecondaryVertexBPFJetTags +
       process.newPFCombinedSecondaryVertexMVABPFJetTags )
    )

process.newPFJetBtagging = cms.Sequence(
    process.newPFJetBtaggingIP +
    process.newPFJetBtaggingSV )

process.newPFBtaggingSequence = cms.Sequence(
    process.newPFJetTracksAssociator *
       process.newPFJetBtagging )


# Now retrieve them in the NTupleProducer
process.analyze.tag_btags = ['newPFTrackCountingHighEffBJetTags',
                             'newPFTrackCountingHighPurBJetTags',
                             'newPFSimpleSecondaryVertexHighEffBJetTags',
                             'newPFSimpleSecondaryVertexHighPurBJetTags',
                             'newPFCombinedSecondaryVertexBPFJetTags',
                             'newPFCombinedSecondaryVertexMVABPFJetTags',
                             'newPFJetProbabilityBPFJetTags',
                             'newPFJetBProbabilityBPFJetTags']
			     
			     
# parton flavour
process.load("PhysicsTools.JetMCAlgos.CaloJetsMCFlavour_cfi")

import FWCore.ParameterSet.Config as cms

process.AK5PFbyRef = process.AK5byRef.clone(
  jets = "ak5PFJets"
)

process.AK5PFbyValAlgo = process.AK5byValAlgo.clone(
  srcByReference = "AK5PFbyRef"
)

##########################################################################
### PF isolation settings ################################################
from CommonTools.ParticleFlow.Tools.pfIsolation import setupPFElectronIso, setupPFPhotonIso
process.eleIsoSequence = setupPFElectronIso(process, 'gsfElectrons', newpostfix='Standard')
process.analyze.tag_elepfisosEvent = ['elPFIsoValueCharged03PFIdStandard' , 'elPFIsoValueCharged04PFIdStandard'  ,
                                      'elPFIsoValueNeutral03PFIdStandard' , 'elPFIsoValueNeutral04PFIdStandard'  ,
                                      'elPFIsoValueGamma03PFIdStandard'   , 'elPFIsoValueGamma04PFIdStandard'  ]

process.phoIsoSequence = setupPFPhotonIso(process, 'photons')

process.load("DiLeptonAnalysis.NTupleProducer.leptonPFIsolation_cff")
process.pfPileUp.Vertices = 'goodVertices'
process.analyze.tag_muonpfisosCustom = ['muonPFIsoChHad03'   , 'muonPFIsoChHad04'   ,
                                        'muonPFIsoNHad03'    , 'muonPFIsoNHad04'    ,
                                        'muonPFIsoPhoton03'  , 'muonPFIsoPhoton04'  ,
                                        'muonRadPFIsoChHad03', 'muonRadPFIsoChHad04',
                                        'muonRadPFIsoNHad03' , 'muonRadPFIsoNHad04' ,
                                        'muonRadPFIsoPhoton03', 'muonRadPFIsoPhoton04' ]
process.analyze.tag_elepfisosCustom  = ['electronPFIsoChHad03'    , 'electronPFIsoChHad04'     ,
                                        'electronPFIsoNHad03'     , 'electronPFIsoNHad04'      ,
                                        'electronPFIsoPhoton03'   , 'electronPFIsoPhoton04'    ,
                                        'electronRadPFIsoChHad03' , 'electronRadPFIsoChHad04'  ,
                                        'electronRadPFIsoNHad03'  , 'electronRadPFIsoNHad04'   ,
                                        'electronRadPFIsoPhoton03', 'electronRadPFIsoPhoton04' ]

####################################################################################
### Tau ID
process.load("RecoTauTag.Configuration.RecoPFTauTag_cff")
from PhysicsTools.PatAlgos.patSequences_cff import tauIsoDepositPFCandidates,tauIsoDepositPFChargedHadrons,tauIsoDepositPFNeutralHadrons,tauIsoDepositPFGammas,patTaus
process.tauIsoDepositPFCandidates = tauIsoDepositPFCandidates.clone()
process.tauIsoDepositPFChargedHadrons = tauIsoDepositPFChargedHadrons.clone()
process.tauIsoDepositPFNeutralHadrons = tauIsoDepositPFNeutralHadrons.clone()
process.tauIsoDepositPFGammas = tauIsoDepositPFGammas.clone()
process.patTaus = patTaus.clone()
### Remove MC matching
process.patTaus.addGenJetMatch   = False
process.patTaus.embedGenJetMatch = False
process.patTaus.genJetMatch      = ''
process.patTaus.genParticleMatch = ''
process.patTaus.addGenMatch = False
process.patTaus.embedGenMatch    = False

process.patTaus.tauIDSources = cms.PSet(
    patTaus.tauIDSources,
    byVLooseChargedIsolation = cms.InputTag("hpsPFTauDiscriminationByVLooseChargedIsolation"),
    byLooseChargedIsolation = cms.InputTag("hpsPFTauDiscriminationByLooseChargedIsolation"),
    byMediumChargedIsolation = cms.InputTag("hpsPFTauDiscriminationByMediumChargedIsolation"),
    byTightChargedIsolation = cms.InputTag("hpsPFTauDiscriminationByTightChargedIsolation"),
    byVLooseIsolation = cms.InputTag("hpsPFTauDiscriminationByVLooseIsolation"),
    byLooseIsolation = cms.InputTag("hpsPFTauDiscriminationByLooseIsolation"),
    byMediumIsolation = cms.InputTag("hpsPFTauDiscriminationByMediumIsolation"),
    byTightIsolation = cms.InputTag("hpsPFTauDiscriminationByTightIsolation"),
    byVLooseIsolationDBSumPtCorr = cms.InputTag("hpsPFTauDiscriminationByVLooseIsolationDBSumPtCorr"),
    byLooseIsolationDBSumPtCorr = cms.InputTag("hpsPFTauDiscriminationByLooseIsolationDBSumPtCorr"),
    byMediumIsolationDBSumPtCorr = cms.InputTag("hpsPFTauDiscriminationByMediumIsolationDBSumPtCorr"),
    byTightIsolationDBSumPtCorr = cms.InputTag("hpsPFTauDiscriminationByTightIsolationDBSumPtCorr"),
    byVLooseCombinedIsolationDBSumPtCorr = cms.InputTag("hpsPFTauDiscriminationByVLooseCombinedIsolationDBSumPtCorr"),
    byLooseCombinedIsolationDBSumPtCorr = cms.InputTag("hpsPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr"),
    byMediumCombinedIsolationDBSumPtCorr = cms.InputTag("hpsPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr"),
    byTightCombinedIsolationDBSumPtCorr = cms.InputTag("hpsPFTauDiscriminationByTightCombinedIsolationDBSumPtCorr"),
    #new 2012 ID's
    byLooseCombinedIsolationDBSumPtCorr3Hits = cms.InputTag("hpsPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr3Hits"),
    byMediumCombinedIsolationDBSumPtCorr3Hits = cms.InputTag("hpsPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr3Hits"),
    byTightCombinedIsolationDBSumPtCorr3Hits = cms.InputTag("hpsPFTauDiscriminationByTightCombinedIsolationDBSumPtCorr3Hits"),
    )
#Only an obvious and loose selection
process.selectedNewTaus = cms.EDFilter("PATTauSelector",
                                       src = cms.InputTag("patTaus"),
                                       cut = cms.string("tauID('decayModeFinding')")
                                       )
process.newTaus = cms.Sequence(process.tauIsoDepositPFCandidates+process.tauIsoDepositPFChargedHadrons+process.tauIsoDepositPFNeutralHadrons+process.tauIsoDepositPFGammas+process.patTaus*process.selectedNewTaus)


###############################################################################################################
########################################Special for Tag&Probe##################################################
###############################################################################################################

process.superClusters = cms.EDProducer("SuperClusterMerger",
   src = cms.VInputTag(cms.InputTag("correctedHybridSuperClusters","", "RECO"),
                       cms.InputTag("correctedMulti5x5SuperClustersWithPreshower","", "RECO") )
)

process.superClusterCands = cms.EDProducer("ConcreteEcalCandidateProducer",
   src = cms.InputTag("superClusters"),
   particleType = cms.int32(11),
)

process.goodSuperClusters = cms.EDFilter("CandViewSelector",
      src = cms.InputTag("superClusterCands"),
      cut = cms.string("abs(eta)<2.5 && !(1.4442< abs(eta) <1.566) && et> 10"),
      filter = cms.bool(False)
)

process.JetsToRemoveFromSuperCluster = cms.EDFilter("CaloJetSelector",
    src = cms.InputTag("ak5CaloJets"),
    cut = cms.string('pt>5 && energyFractionHadronic > 0.15')
)

process.goodSuperClustersClean = cms.EDProducer("CandViewCleaner",
    srcObject = cms.InputTag("goodSuperClusters"),
    module_label = cms.string(""),
    srcObjectsToRemove = cms.VInputTag(cms.InputTag("JetsToRemoveFromSuperCluster")),
    deltaRMin = cms.double(0.1)
)

process.sc_sequence = cms.Sequence(
    process.superClusters +
    process.superClusterCands +
    process.goodSuperClusters +
    process.JetsToRemoveFromSuperCluster +
    process.goodSuperClustersClean
)

# Quark/gluon discrimination
process.load('QuarkGluonTagger.EightTeV.QGTagger_RecoJets_cff')
process.QGTagger.srcJets = cms.InputTag("ak5PFJets")
if options.runon == 'data':
   process.QGTagger.jec = cms.untracked.string('ak5PFL1FastL2L3Residual')
else:
   process.QGTagger.jec = cms.untracked.string('ak5PFL1FastL2L3')
process.kt6PFJetsForQGSyst = process.kt4PFJets.clone( rParam = 0.6, doRhoFastjet = True )
process.kt6PFJetsForQGSyst.Rho_EtaMax = cms.double(2.4)

# Beam scraping veto
process.scrapingVeto = cms.EDFilter("FilterOutScraping",
                                    applyfilter = cms.untracked.bool  (True),
                                    debugOn     = cms.untracked.bool  (False),
                                    numtrack    = cms.untracked.uint32(10),
                                    thresh      = cms.untracked.double(0.25)
                                    )


# Diphoton trigger filter
process.diphotonTriggerSelection2012 = cms.EDFilter( "TriggerResultsFilter",
                                                     triggerConditions = cms.vstring(
   'HLT_Photon36_R9Id85_OR_CaloId10_Iso50_Photon22_R9Id85_OR_CaloId10_Iso50_v*'
   ),
                                                     hltResults = cms.InputTag( "TriggerResults", "", "HLT" ),
                                                     l1tResults = cms.InputTag( "gtDigis" ),
                                                     l1tIgnoreMask = cms.bool( False ),
                                                     l1techIgnorePrescales = cms.bool( False ),
                                                     daqPartitions = cms.uint32( 1 ),
                                                     throw = cms.bool( False )
                                                     )

# Diphoton trigger filter
process.diphotonTriggerSelection2011 = cms.EDFilter( "TriggerResultsFilter",
                                                     triggerConditions = cms.vstring(
   'HLT_Photon20_R9Id_Photon18_R9Id OR HLT_Photon20_R9Id_Photon18_R9Id_v*',
   'HLT_Photon26_CaloIdL_IsoVL_Photon18 OR HLT_Photon26_CaloIdL_IsoVL_Photon18_v*',
   'HLT_Photon26_CaloIdL_IsoVL_Photon18_CaloIdL_IsoVL OR HLT_Photon26_CaloIdL_IsoVL_Photon18_CaloIdL_IsoVL_v*',
   'HLT_Photon26_CaloIdL_IsoVL_Photon18_R9Id OR HLT_Photon26_CaloIdL_IsoVL_Photon18_R9Id_v*',
   'HLT_Photon26_CaloIdXL_IsoXL_Photon18_CaloIdXL_IsoXL OR HLT_Photon26_CaloIdXL_IsoXL_Photon18_CaloIdXL_IsoXL_v*',
   'HLT_Photon26_CaloIdXL_IsoXL_Photon18_R9Id OR HLT_Photon26_CaloIdXL_IsoXL_Photon18_R9Id_v*',
   'HLT_Photon26_IsoVL_Photon18 OR HLT_Photon26_IsoVL_Photon18_v*',
   'HLT_Photon26_IsoVL_Photon18_IsoVL OR HLT_Photon26_IsoVL_Photon18_IsoVL_v*',
   'HLT_Photon26_Photon18 OR HLT_Photon26_Photon18_v*',
   'HLT_Photon26_R9Id_Photon18_CaloIdL_IsoVL OR HLT_Photon26_R9Id_Photon18_CaloIdL_IsoVL_v*',
   'HLT_Photon26_R9Id_Photon18_CaloIdXL_IsoXL OR HLT_Photon26_R9Id_Photon18_CaloIdXL_IsoXL_v*',
   'HLT_Photon26_R9Id_Photon18_R9Id OR HLT_Photon26_R9Id_Photon18_R9Id_v*',
   'HLT_Photon32_CaloIdL_Photon26_CaloIdL OR HLT_Photon32_CaloIdL_Photon26_CaloIdL_v*',
   'HLT_Photon36_CaloIdL_IsoVL_Photon22_CaloIdL_IsoVL OR HLT_Photon36_CaloIdL_IsoVL_Photon22_CaloIdL_IsoVL_v*',
   'HLT_Photon36_CaloIdL_IsoVL_Photon22_R9Id OR HLT_Photon36_CaloIdL_IsoVL_Photon22_R9Id_v*',
   'HLT_Photon36_CaloIdL_Photon22_CaloIdL OR HLT_Photon36_CaloIdL_Photon22_CaloIdL_v*',
   'HLT_Photon36_R9Id_Photon22_CaloIdL_IsoVL OR HLT_Photon36_R9Id_Photon22_CaloIdL_IsoVL_v*',
   'HLT_Photon36_R9Id_Photon22_R9Id OR HLT_Photon36_R9Id_Photon22_R9Id_v*' ),
                                                     hltResults = cms.InputTag( "TriggerResults", "", "HLT" ),
                                                     l1tResults = cms.InputTag( "gtDigis" ),
                                                     l1tIgnoreMask = cms.bool( False ),
                                                     l1techIgnorePrescales = cms.bool( False ),
                                                     daqPartitions = cms.uint32( 1 ),
                                                     throw = cms.bool( False )
                                                     )

if (options.is2011):
   process.diphotonTriggerSelection = process.diphotonTriggerSelection2011
else:
   process.diphotonTriggerSelection = process.diphotonTriggerSelection2012


# Preselection: at least one photon passing preselection

process.preselPhotons = cms.EDFilter(
   "PhotonRefSelector",
   src = cms.InputTag( "photons" ),
   cut = cms.string(
   "(hadronicOverEm<0.15 && (abs(superCluster.eta)<2.5) && !(1.4442<abs(superCluster.eta)<1.566) && (superCluster.energy*sin(superCluster.position.theta)>15))"
   " && ((abs(superCluster.eta)>1.5) || (sigmaIetaIeta<0.017)) "
   " && ((abs(superCluster.eta)<1.5) || (sigmaIetaIeta<0.040)) "
   )
   )
process.photonFilter = cms.EDFilter("CandViewCountFilter",
                                    src = cms.InputTag("preselPhotons"),
                                    minNumber = cms.uint32(1)
                                    )
process.photonPreselFilter = cms.Sequence(process.preselPhotons*process.photonFilter)





#### Path ######################################################################

process.p = cms.Path(
        process.myPartons *
        process.scrapingVeto *
        process.AK5PFbyRef *
        process.AK5PFbyValAlgo *
        process.goodVertices * # Filter
        process.hcallasereventfilter2012 * # Filter
#        process.diphotonTriggerSelection *
#        process.photonPreselFilter *
        (
         (process.photonPartonMatch
#	*process.printGenParticles*process.printPhotons*process.printPartons
	 )
	# cleaning 
	+ process.HBHENoiseFilterResultProducer              # tagging mode
	+ process.hcalLaserEventFilter                       # tagging mode
	+ process.EcalDeadCellTriggerPrimitiveFilter         # tagging mode
	+ process.trackingFailureFilter                      # tagging mode
	+ process.eeBadScFilter                              # tagging mode
        + process.ecalLaserCorrFilter                        # tagging mode
	+ process.trkPOGFilters                              # tagging mode
#	+ process.dump  
	# end cleaning 
	+ process.pfIsolationAllSequence
	+ process.kt6PFJetsForIso
	+ process.newBtaggingSequence
	+ process.newPFBtaggingSequence
       	+ process.mygenjets
       	+ process.metCorSequence
        + process.pfParticleSelectionSequence
 	+ process.eleIsoSequence
        + process.phoIsoSequence
	+ process.recoTauClassicHPSSequence
	+ process.PFTau
	+ process.newTaus
        + process.patPF2PATSequencePFCHS
        + process.sc_sequence
        + process.kt6PFJetsForQGSyst
        + process.QuarkGluonTagger
        + process.ak5PFJetsCorrected
        + process.analyze
       	)
   )
   

process.out.outputCommands = ['drop *', 
                              'keep *_analyze_*_'+process.name_(),
			      'keep *_HBHENoiseFilterResultProducer_*_'+process.name_(),
			      'keep *_EcalDeadCellTriggerPrimitiveFilter_*_'+process.name_(),
			      'keep *_hcalLaserEventFilter_*_'+process.name_(),
			      'keep *_trackingFailureFilter_*_'+process.name_(),
			      'keep *_eeBadScFilter_*_'+process.name_(),
			      'keep *_ecalLaserCorrFilter_*_'+process.name_(),
			      'keep *_manystripclus53X_*_'+process.name_(),
			      'keep *_toomanystripclus53X_*_'+process.name_(),
			      'keep *_logErrorTooManyClusters_*_'+process.name_()
			      ] 


process.outpath = cms.EndPath(process.out)


# remove ak5GenJets from the path if it will run on data
if options.runon=='data':
    process.p.remove(process.mygenjets)
    process.p.remove(process.photonPartonMatch)
    process.p.remove(process.myPartons)
    process.p.remove(process.AK5PFbyRef)
    process.p.remove(process.AK5PFbyValAlgo)
else:
   process.p.remove(process.scrapingVeto)
#   process.p.remove(process.diphotonTriggerSelection)
if options.ModelScan==True:
    process.p.remove(process.hcalLaserEventFilter)
if options.FastSim==True:
    process.p.remove(process.hcalLaserEventFilter)
    process.p.remove(process.HBHENoiseFilterResultProducer)
if options.doPhotonStuff==True:
   process.p.remove(process.PFTau)
   process.p.remove(process.newTaus)
   process.p.remove(process.patPF2PATSequencePFCHS)
else:
   process.p.remove(process.kt6PFJetsForQGSyst)
   process.p.remove(process.QuarkGluonTagger)
