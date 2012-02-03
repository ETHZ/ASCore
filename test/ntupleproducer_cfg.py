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
options.register ('ModelScan', # register 'runon' option
                  False,  # the default value
                  VarParsing.VarParsing.multiplicity.singleton, # singleton or list
                  VarParsing.VarParsing.varType.bool,         # string, int, or float
                  "If you are dealing with a model scan, set this to True, otherwise to False (default)")
# get and parse the command line arguments
# set NTupleProducer defaults (override the output, files and maxEvents parameter)
options.files= 'file:/shome/pnef/SUSY/reco-data/mc/GJets_TuneZ2_200_HT_inf_7TeV-madgraph_AODSIM_PU_S4_START42_V11-v1_0000_00F3A238-FFCC-E011-AA04-0026B94D1AEF.root'
options.maxEvents = -1# If it is different from -1, string "_numEventXX" will be added to the output file name
# Now parse arguments from command line (might overwrite defaults)
options.parseArguments()
options.output='NTupleProducer_42X_'+options.runon+'.root'

### Running conditions #########################################################
# See https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideFrontierConditions
# to check what cff to use
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
if options.runon=='data':
    # CMSSW_4_2
    process.GlobalTag.globaltag = "GR_R_42_V19::All"
else:
    # CMSSW_4_2_X:
    process.GlobalTag.globaltag = "START42_V13::All"


### Input/Output ###############################################################
# Input
process.source = cms.Source("PoolSource",
      fileNames = cms.untracked.vstring(options.files)
      # Enable if you see duplicate error:
      # duplicateCheckMode = cms.untracked.string("noDuplicateCheck")
)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents) )
process.ModelScan = cms.untracked.PSet( input = cms.untracked.bool(options.ModelScan) )
# Output
process.TFileService = cms.Service("TFileService",
# Keep track of the type of data source and reco type in the ntuple file name
	fileName = cms.string(options.output),
	closeFileFast = cms.untracked.bool(True)
)

### Electron ID ##############################################################
process.load("DiLeptonAnalysis.NTupleProducer.simpleEleIdSequence_cff")

### Jet Corrections ##########################################################
# See https://twiki.cern.ch/twiki/bin/view/CMS/WorkBookJetEnergyCorrections
# note: this runs the L1Fast-Jet corrections for PF jets. not applied on Calo
process.load('JetMETCorrections.Configuration.DefaultJEC_cff')
process.load('RecoJets.Configuration.RecoPFJets_cff')
# Turn-on the FastJet density calculation -----------------------
process.kt6PFJets.doRhoFastjet = True
# process.kt6PFJets.Rho_EtaMax   = cms.double(4.4) # this is the default value in 4_2
# process.kt6PFJets.Ghost_EtaMax = cms.double(5.0) # this is the default value in 4_2
# Turn-on the FastJet jet area calculation 
process.ak5PFJets.doAreaFastjet = True
process.ak5PFJets.Rho_EtaMax = process.kt6PFJets.Rho_EtaMax

### JES MET Corrections ########################################################
from JetMETCorrections.Type1MET.MetType1Corrections_cff import metJESCorAK5CaloJet

process.metMuonJESCorAK5 = metJESCorAK5CaloJet.clone()
process.metMuonJESCorAK5.inputUncorJetsLabel = "ak5CaloJets"
process.metMuonJESCorAK5.inputUncorMetLabel = "corMetGlobalMuons"
process.metCorSequence = cms.Sequence(process.metMuonJESCorAK5)

### Cleaning ###################################################################
# flag HB/HE noise
# https://twiki.cern.ch/twiki/bin/view/CMS/HcalNoiseInfoLibrary
# https://twiki.cern.ch/twiki/bin/view/CMS/HBHEAnomalousSignals2011
# the next two lines produce the HCAL noise summary flag
from CommonTools.RecoAlgos.HBHENoiseFilterResultProducer_cfi import *
# std config with cut on iso:
process.HBHENoiseFilterResultProducerIso = HBHENoiseFilterResultProducer.clone() 
# HCAL DPG recomended config:
process.HBHENoiseFilterResultProducerStd = HBHENoiseFilterResultProducer.clone()
process.HBHENoiseFilterResultProducerStd.minIsolatedNoiseSumE        = cms.double(999999.)
process.HBHENoiseFilterResultProducerStd.minNumIsolatedNoiseChannels = cms.int32(999999)
process.HBHENoiseFilterResultProducerStd.minIsolatedNoiseSumEt       = cms.double(999999.)



# RA2 RecHitFilter: tagging mode
# https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFilters#RecovRecHitFilter
process.load('SandBox.Skims.recovRecHitFilter_cfi')
process.recovRecHitFilter.TaggingMode           = cms.bool(True)


# ECAL dead cells: this is not a filter. Only a flag is stored.
from JetMETAnalysis.ecalDeadCellTools.EcalDeadCellEventFilter_cfi import *
process.ecalDeadCellTPfilter                           = EcalDeadCellEventFilter.clone()
process.ecalDeadCellTPfilter.tpDigiCollection          = cms.InputTag("ecalTPSkim")
process.ecalDeadCellTPfilter.etValToBeFlagged          = cms.double(63.75)
process.ecalDeadCellTPfilter.ebReducedRecHitCollection = cms.InputTag("reducedEcalRecHitsEB")
process.ecalDeadCellTPfilter.eeReducedRecHitCollection = cms.InputTag("reducedEcalRecHitsEE")
process.ecalDeadCellTPfilter.taggingMode               = cms.bool( True )

# See for example DPGAnalysis/Skims/python/MinBiasPDSkim_cfg.py
# require scraping filter
process.scrapingVeto = cms.EDFilter("FilterOutScraping",
     applyfilter = cms.untracked.bool(True),
     debugOn = cms.untracked.bool(False),
     numtrack = cms.untracked.uint32(10),
     thresh = cms.untracked.double(0.25)
     )

# Get a list of good primary vertices, in 42x, these are DAF vertices
process.goodVertices = cms.EDFilter("VertexSelector",
	src = cms.InputTag("offlinePrimaryVertices"),
	cut = cms.string("!isFake && ndof > 4 && abs(z) <= 24 && position.Rho <= 2"),
	filter = cms.bool(False)
)

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

# Compute the mean pt per unit area (rho) from the
# various PFchs inputs
# -> note: when using PFnoPileUp, rho must be computed from the kt6-pfjets
#          with charged pf-candidates not coming from the PV already subtracted!
#          since the kt6-pfjets come from the pfNoElectron collection (which
#          is different for each pf2pat branch) we run it three times...
from RecoJets.JetProducers.kt4PFJets_cfi import kt4PFJets
process.kt6PFJetsPF2 = kt4PFJets.clone(
		rParam = cms.double(0.6),
		src = cms.InputTag('pfNoElectronPF2'),
		doAreaFastjet = cms.bool(True),
		doRhoFastjet = cms.bool(True)
	)
process.kt6PFJetsPF3 = kt4PFJets.clone(
		rParam = cms.double(0.6),
		src = cms.InputTag('pfNoElectronPF3'),
		doAreaFastjet = cms.bool(True),
		doRhoFastjet = cms.bool(True)
	)
process.kt6PFJetsPFAntiIso = kt4PFJets.clone(
		rParam = cms.double(0.6),
		src = cms.InputTag('pfNoElectronPFAntiIso'),
		doAreaFastjet = cms.bool(True),
		doRhoFastjet = cms.bool(True)
	)

process.load("RecoTauTag/RecoTau/PFRecoTauProducer_cfi")
process.load("RecoTauTag.Configuration.RecoPFTauTag_cff")

### Configuration in common to all collections
pfPostfixes = [ 'PFAntiIso','PF2','PF3' ]
for pf in pfPostfixes:

    usePF2PAT(process,runPF2PAT=True, jetAlgo='AK5', runOnMC=(options.runon != 'data'), postfix=pf) 

    from PhysicsTools.PatAlgos.tools.pfTools import adaptPFTaus
    adaptPFTaus(process,"hpsPFTau",pf)

    # configure PFnoPU
    getattr(process,'pfPileUp'+pf).Enable              = True
    getattr(process,'pfPileUp'+pf).checkClosestZVertex = cms.bool(False)
    getattr(process,'pfPileUp'+pf).Vertices            = cms.InputTag('goodVertices')
    getattr(process,'pfJets'+pf).doAreaFastjet         = True
    getattr(process,'pfJets'+pf).doRhoFastjet          = False
    getattr(process,'pfNoPileUp'+pf).enable            = False  # !!! this is a switch to turn on and off the pile-up removal !!!!
    # Add the KT6 producer to the sequence (see below)
    getattr(process,"patPF2PATSequence"+pf).replace(
	getattr(process,"pfNoElectron"+pf),
	getattr(process,"pfNoElectron"+pf)*getattr(process,'kt6PFJets'+pf)) # pfNoElectron is src of kt6PFJets
    # Jet corrections 
    getattr(process,'patJetCorrFactors'+pf).levels = ['L1FastJet', 'L2Relative', 'L3Absolute']
    getattr(process,'patJetCorrFactors'+pf).rho    = cms.InputTag('kt6PFJets'+pf,'rho')
#   getattr(process,'patJetCorrFactors'+pf).payload= cms.string('AK5PFchs')  # !!! if pfnoPU is turned ON use these JEC!!!
    getattr(process,'patJetCorrFactors'+pf).payload= cms.string('AK5PF')     # !!! if pfnoPU is turned OFF use these JEC!!!

    # set to false to disable jet to be cleaned from Taus
    getattr(process,"pfNoTau"+pf).enable      = False
    
    # turn to false when running on data and MC (for the moment)
    getattr(process, 'patElectrons'+pf).embedGenMatch = False
    getattr(process, 'patMuons'+pf).embedGenMatch = False

    # PatElectronID
    getattr(process, 'patElectrons'+pf).addElectronID = cms.bool(True)
    getattr(process, 'patElectrons'+pf).electronIDSources = cms.PSet(
	simpleEleId95relIso= cms.InputTag("simpleEleId95relIso"),
	simpleEleId90relIso= cms.InputTag("simpleEleId90relIso"),
	simpleEleId85relIso= cms.InputTag("simpleEleId85relIso"),
      	simpleEleId80relIso= cms.InputTag("simpleEleId80relIso"),
   	simpleEleId70relIso= cms.InputTag("simpleEleId70relIso"),
       	simpleEleId60relIso= cms.InputTag("simpleEleId60relIso"),
    	simpleEleId95cIso  = cms.InputTag("simpleEleId95cIso"),
       	simpleEleId90cIso  = cms.InputTag("simpleEleId90cIso"),
    	simpleEleId85cIso  = cms.InputTag("simpleEleId85cIso"),
        simpleEleId80cIso  = cms.InputTag("simpleEleId80cIso"),
    	simpleEleId70cIso  = cms.InputTag("simpleEleId70cIso"),
        simpleEleId60cIso  = cms.InputTag("simpleEleId60cIso"),
    )

      
    # muon vertex 
    getattr(process, 'pfMuonsFromVertex'+pf).vertices=cms.InputTag("goodVertices") # require muon to come from the good vertices as defined above
    getattr(process, 'pfMuonsFromVertex'+pf).d0Cut   =cms.double(0.02) # transverse IP w.r.t. PV
    getattr(process, 'pfMuonsFromVertex'+pf).dzCut   =cms.double(1.)   # longitudinal IP w.r.t. PV
    getattr(process, 'pfSelectedMuons'+pf).cut = cms.string(
 		"abs( eta ) < 2.4 && pt > 5 && muonRef().isNonnull()" 
    )
    
    # electron vertex
    getattr(process, 'pfElectronsFromVertex'+pf).vertices=cms.InputTag("goodVertices") # require eles to come from the good vertices as defined above
    getattr(process, 'pfElectronsFromVertex'+pf).d0Cut   =cms.double(0.04) # transverse IP w.r.t. PV
    getattr(process, 'pfElectronsFromVertex'+pf).dzCut   =cms.double(1.)   # longitudinal IP w.r.t. PV
    getattr(process, 'pfSelectedElectrons'+pf).cut = cms.string(
 		"abs( eta ) < 2.4 && pt > 5 && gsfTrackRef().isNonnull()" 
		+"&& gsfTrackRef().trackerExpectedHitsInner().numberOfHits() <= 1"       # conv rejection
    )


    # tau cut
    getattr(process, 'pfTaus'+pf).cut = cms.string(
 		"abs( eta ) < 2.4 && pt > 10  && abs( charge ) == 1. && leadPFChargedHadrCand().isNonnull()" 
                )


    # ISOLATION
    getattr(process, 'pfIsolatedMuons'+pf)    .combinedIsolationCut = cms.double(0.20)
    getattr(process, 'pfIsolatedElectrons'+pf).combinedIsolationCut = cms.double(0.20)




### Specific to first PF collection: AntiIsolated electrons and muons: Isolation cuts is set to 2
# ID Cuts
process.pfMuonsIDPFAntiIso = cms.EDFilter("MuonIDPFCandidateSelector", 
		src = cms.InputTag("pfMuonsFromVertexPFAntiIso"), 
		cut = cms.string( 
			"muonID('GlobalMuonPromptTight') && " + 
			"isTrackerMuon &&" + 
			"track.numberOfValidHits > 10 && "+ 
			"track.hitPattern.numberOfValidPixelHits > 0" 
			), 
		) 
process.pfSelectedMuonsPFAntiIso.src = 'pfMuonsIDPFAntiIso' 
process.pfMuonSequencePFAntiIso.replace( process.pfSelectedMuonsPFAntiIso, 
		process.pfMuonsIDPFAntiIso * 
		process.pfSelectedMuonsPFAntiIso 
		) 
process.pfIsolatedMuonsPFAntiIso.combinedIsolationCut     = cms.double(2.0) # set min isolation to 2
process.pfIsolatedElectronsPFAntiIso.combinedIsolationCut = cms.double(2.0) # set min isolation to 2


process.selectedPatTausPFAntiIso.cut = cms.string("tauID('byVLooseIsolation')")


### Specific to second PF collection: TIGHT
# ID cuts
process.pfMuonsIDPF2 = cms.EDFilter("MuonIDPFCandidateSelector", 
		 src = cms.InputTag("pfMuonsFromVertexPF2"), 
		 cut = cms.string( 
			"muonID('GlobalMuonPromptTight') && " + 
			"isTrackerMuon && numberOfMatches > 1 && " + 
			"track.numberOfValidHits > 10 && "+ 
			"track.hitPattern.numberOfValidPixelHits > 0 &&" + 
			"globalTrack.ptError/pt < 0.1"
			), 
		 ) 
process.pfSelectedMuonsPF2.src = 'pfMuonsIDPF2' 
process.pfMuonSequencePF2.replace( process.pfSelectedMuonsPF2, 
		process.pfMuonsIDPF2 * 
		process.pfSelectedMuonsPF2 
		) 
process.pfElectronsIDPF2 = cms.EDFilter("ElectronIDPFCandidateSelector", 
		src = cms.InputTag("pfElectronsFromVertexPF2"), 
		recoGsfElectrons = cms.InputTag('gsfElectrons'), 
		electronIdMap = cms.InputTag('simpleEleId90relIso'), 
		bitsToCheck = cms.vstring('id'), 
		) 
process.pfSelectedElectronsPF2.src = 'pfElectronsIDPF2' 
process.pfElectronSequencePF2.replace( process.pfSelectedElectronsPF2, 
		process.pfElectronsIDPF2 * 
		process.pfSelectedElectronsPF2 
		) 

process.selectedPatTausPF2.cut = cms.string("tauID('againstElectronTight') && tauID('againstMuonTight')")


### Specific to second PF collection: LOOSE
# ID cuts
process.pfMuonsIDPF3 = cms.EDFilter("MuonIDPFCandidateSelector", 
		 src = cms.InputTag("pfMuonsFromVertexPF3"), 
		 cut = cms.string( 
			"muonID('GlobalMuonPromptTight') && " + 
			"isTrackerMuon && " + 
			"track.numberOfValidHits > 10 && "+ 
			"track.hitPattern.numberOfValidPixelHits > 0 " 
			), 
		 ) 
process.pfSelectedMuonsPF3.src = 'pfMuonsIDPF3' 
process.pfMuonSequencePF3.replace( process.pfSelectedMuonsPF3, 
		process.pfMuonsIDPF3 * 
		process.pfSelectedMuonsPF3 
		) 


process.selectedPatTausPF3.cut = cms.string("tauID('againstElectronLoose') && tauID('againstMuonLoose')")



### GenJets ####################################################################
# produce ak5GenJets (collection missing in case of some Spring10 samples)
process.load("RecoJets.Configuration.GenJetParticles_cff")
process.load("RecoJets.Configuration.RecoGenJets_cff")
process.mygenjets = cms.Sequence( process.genParticlesForJets * process.ak5GenJets )

### Analysis configuration #####################################################
process.load("DiLeptonAnalysis.NTupleProducer.ntupleproducer_cfi")
process.analyze.isRealData = cms.untracked.bool(options.runon=='data')
process.analyze.isModelScan = cms.untracked.bool(options.ModelScan)

# Add some jet collections
process.analyze.jets = (
   # Calo jets
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
    cms.PSet( prefix = cms.untracked.string('PF2PATAntiIso'),
              tag = cms.untracked.InputTag('patJetsPFAntiIso'),
              isPat = cms.untracked.bool(True),
              tag_jetTracks  = cms.untracked.InputTag('ak5JetTracksAssociatorAtVertex'),
              sel_minpt  = cms.double(15.0),
              sel_maxeta = process.analyze.sel_maxjeta,
              # The corrections are irrelevant for PF2PAT
              corrections = cms.string(''), 
              btag_matchdeltaR = cms.double(0.25),
              ),
    cms.PSet( prefix = cms.untracked.string('PF2PAT2'),
              tag = cms.untracked.InputTag('patJetsPF2'),
              isPat = cms.untracked.bool(True),
              tag_jetTracks  = cms.untracked.InputTag('ak5JetTracksAssociatorAtVertex'),
              sel_minpt  = cms.double(15.0),
              sel_maxeta = process.analyze.sel_maxjeta,
              # The corrections are irrelevant for PF2PAT
              corrections = cms.string(''), 
              btag_matchdeltaR = cms.double(0.25),
              ),
    cms.PSet( prefix = cms.untracked.string('PF2PAT3'),
              tag = cms.untracked.InputTag('patJetsPF3'),
              isPat = cms.untracked.bool(True),
              tag_jetTracks  = cms.untracked.InputTag('ak5JetTracksAssociatorAtVertex'),
              sel_minpt  = cms.double(15.0),
              sel_maxeta = process.analyze.sel_maxjeta,
              # The corrections are irrelevant for PF2PAT
              corrections = cms.string(''), 
              btag_matchdeltaR = cms.double(0.25),
              ),
    )
# # Add residual correction for running on data
# # taken from local sqlite file. see: https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookJetEnergyCorrections#JetEnCor2011V2 
if options.runon == 'data':
	process.analyze.jetCorrs = process.analyze.jetCorrs.value() + 'Residual'
        for extJet in process.analyze.jets:
        	extJet.corrections = extJet.corrections.value() + 'Residual'
	for pf in pfPostfixes:            
		getattr(process,'patJetCorrFactors'+pf).levels.extend( ['L2L3Residual'] )
	
# Add some PF lepton collections
process.analyze.leptons = (
    # PF Electrons
    cms.PSet( type = cms.untracked.string('electron'),
              prefix = cms.untracked.string('PfElAntiIso'),
              tag = cms.untracked.InputTag('patElectronsPFAntiIso'),
              sel_minpt = process.analyze.sel_minelpt,
              sel_maxeta = process.analyze.sel_maxeleta,
              maxnobjs = cms.untracked.uint32(20)
              ),
    # PF Muons
    cms.PSet( type = cms.untracked.string('muon'),
              prefix = cms.untracked.string('PfMuAntiIso'),
              tag = cms.untracked.InputTag('patMuonsPFAntiIso'),
              sel_minpt = process.analyze.sel_minelpt,
              sel_maxeta = process.analyze.sel_maxeleta,
              maxnobjs = cms.untracked.uint32(20)
              ),
    # PF Taus
    cms.PSet( type = cms.untracked.string('tau'),
              prefix = cms.untracked.string('PfTauAntiIso'),
              tag = cms.untracked.InputTag('selectedPatTausPFAntiIso'),
              sel_minpt = process.analyze.sel_minelpt,
              sel_maxeta = process.analyze.sel_maxeleta,
              maxnobjs = cms.untracked.uint32(20)
              ),
    # PF Electrons
    cms.PSet( type = cms.untracked.string('electron'),
              prefix = cms.untracked.string('PfEl2'),
              tag = cms.untracked.InputTag('patElectronsPF2'),
              sel_minpt = process.analyze.sel_minelpt,
              sel_maxeta = process.analyze.sel_maxeleta,
              maxnobjs = cms.untracked.uint32(20)
              ),
    # PF Muons
    cms.PSet( type = cms.untracked.string('muon'),
              prefix = cms.untracked.string('PfMu2'),
              tag = cms.untracked.InputTag('patMuonsPF2'),
              sel_minpt = process.analyze.sel_minelpt,
              sel_maxeta = process.analyze.sel_maxeleta,
              maxnobjs = cms.untracked.uint32(20)
              ),
    # PF Taus
    cms.PSet( type = cms.untracked.string('tau'),
              prefix = cms.untracked.string('PfTau2'),
              tag = cms.untracked.InputTag('selectedPatTausPF2'),
              sel_minpt = process.analyze.sel_minelpt,
              sel_maxeta = process.analyze.sel_maxeleta,
              maxnobjs = cms.untracked.uint32(20)
              ),
    # PF Electrons
    cms.PSet( type = cms.untracked.string('electron'),
              prefix = cms.untracked.string('PfEl3'),
              tag = cms.untracked.InputTag('patElectronsPF3'),
              sel_minpt = process.analyze.sel_minelpt,
              sel_maxeta = process.analyze.sel_maxeleta,
              maxnobjs = cms.untracked.uint32(20)
              ),
    # PF Muons
    cms.PSet( type = cms.untracked.string('muon'),
              prefix = cms.untracked.string('PfMu3'),
              tag = cms.untracked.InputTag('patMuonsPF3'),
              sel_minpt = process.analyze.sel_minelpt,
              sel_maxeta = process.analyze.sel_maxeleta,
              maxnobjs = cms.untracked.uint32(20)
              ),
    # PF Taus
    cms.PSet( type = cms.untracked.string('tau'),
              prefix = cms.untracked.string('PfTau3'),
              tag = cms.untracked.InputTag('selectedPatTausPF3'),
              sel_minpt = process.analyze.sel_minelpt,
              sel_maxeta = process.analyze.sel_maxeleta,
              maxnobjs = cms.untracked.uint32(20)
              ),
    )

## Colins Bernet's Particle Based Noise Rejection Filter
#process.load('SandBox.Skims.jetIDFailureFilter_cfi')
#process.jetIDFailure.taggingMode   = cms.bool(True) # events are not filtered, but tagged
#process.jetIDFailure.JetSource     = cms.InputTag('patJetsPF3')
              
# RA2 TrackingFailureFilter
# https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFilters
process.load('SandBox.Skims.trackingFailureFilter_cfi')
process.trackingFailureFilter.JetSource             = cms.InputTag('patJetsPF3')
process.trackingFailureFilter.TrackSource           = cms.InputTag('generalTracks')
process.trackingFailureFilter.VertexSource          = cms.InputTag('goodVertices')
process.trackingFailureFilter.DzTrVtxMax            = cms.double(1)
process.trackingFailureFilter.DxyTrVtxMax           = cms.double(0.2)
process.trackingFailureFilter.MinSumPtOverHT        = cms.double(0.10)
process.trackingFailureFilter.taggingMode           = cms.bool(True)

#### Steve Mrenna's Photon - Parton DR match #######################	      
process.printGenParticles = cms.EDAnalyzer("ParticleListDrawer",
	src = cms.InputTag("partonGenJets"),
	maxEventsToPrint = cms.untracked.int32(10)
)
#
process.printPhotons = cms.EDAnalyzer("ParticleListDrawer",
     src = cms.InputTag("photons"),
     maxEventsToPrint = cms.untracked.int32(10)
)
#
process.printPartons = cms.EDAnalyzer("ParticleListDrawer",
     src = cms.InputTag("myPhotonJetMatch"),
     maxEventsToPrint = cms.untracked.int32(10)
)
#
process.load('DiLeptonAnalysis.NTupleProducer.photonPartonMatch_cfi')

#### DEBUG #####################################################################
#process.dump = cms.EDAnalyzer("EventContentAnalyzer")
# process.SimpleMemoryCheck = cms.Service("SimpleMemoryCheck",
#    ignoreTotal = cms.untracked.int32(1) # number of events to ignore at start (default is one)
# )
# process.ProfilerService = cms.Service("ProfilerService",
#                                      firstEvent = cms.untracked.int32(2),
#                                      lastEvent = cms.untracked.int32(51),
#                                      paths = cms.untracked.vstring(['p'])
#                                      )
#process.Tracer = cms.Service("Tracer")
#process.options = cms.untracked.PSet(
# 	wantSummary = cms.untracked.bool(True)
#)
#### Path ######################################################################

process.p = cms.Path(
	process.scrapingVeto*(
	process.goodVertices
	+ (process.photonPartonMatch
#	*process.printGenParticles*process.printPhotons*process.printPartons
	)
       	+ process.HBHENoiseFilterResultProducerIso
       	+ process.HBHENoiseFilterResultProducerStd
	+ process.ecalDeadCellTPfilter
	+ process.recovRecHitFilter
	+ process.kt6PFJets
	+ process.ak5PFJets
       	+ process.mygenjets
       	+ process.simpleEleIdSequence
       	+ process.metCorSequence
       	+ process.patPF2PATSequencePFAntiIso
       	+ process.patPF2PATSequencePF2
       	+ process.patPF2PATSequencePF3
	+ process.trackingFailureFilter
#	+ process.jetIDFailure
#	+ process.dump
	+ process.analyze

       	)
   )
# remove ak5GenJets from the path if it will run on data
if options.runon=='data':
    process.p.remove(process.mygenjets)
    process.p.remove(process.photonPartonMatch)
if options.runon!='data':
    process.p.remove(process.scrapingVeto)
if options.ModelScan==True:
    process.p.remove(process.HBHENoiseFilterResultProducerIso)
    process.p.remove(process.HBHENoiseFilterResultProducerStd)
