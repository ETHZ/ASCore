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

### Switch for type of run (data, MC) and reconstruction (RECO, PAT, PF) #######
runon = 'data'
#runon = 'MC31x'
#runon = 'MC34x'
#recoType = 'RECO'
#recoType = 'PAT'
recoType = 'PF'

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
	# process.GlobalTag.globaltag = "MC_31X_V3::All"

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
		# 900 GeV data
#		'/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/SD_InterestingEvents-Dec19thSkim_341_v1/0006/5A5773A3-B9ED-DE11-A99A-0026189437F8.root',
#		'/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/SD_InterestingEvents-Dec19thSkim_341_v1/0006/CA156CCC-B9ED-DE11-B022-00261894390B.root'
		'/store/data/BeamCommissioning09/MinimumBias/RECO/Dec19thReReco_341_v1/0006/DA397007-23EE-DE11-A754-00261894388D.root'
		# 2360 GeV data
#        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/SD_InterestingEvents-Dec19thSkim_341_v1/0005/B641315C-ACED-DE11-82E1-0030486792B6.root',
#        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/SD_InterestingEvents-Dec19thSkim_341_v1/0005/8EB877E1-A5ED-DE11-9680-002618943826.root',
#        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/SD_InterestingEvents-Dec19thSkim_341_v1/0005/52B38DBB-A1ED-DE11-AAB0-00248C0BE014.root',
    )
elif runon=='MC31x':
    source_fileNames = cms.untracked.vstring(
		'/store/mc/Summer09/LM0/GEN-SIM-RECO/MC_31X_V3_7TeV-v1/0006/C81C39F6-B4BF-DE11-95E4-001CC47BA216.root'   
    )
elif runon=='MC34x':
    source_fileNames = cms.untracked.vstring(
		'/store/mc/Summer09/MinBias/GEN-SIM-RECO/V16D_900GeV-v1/0001/FCD70794-F216-DF11-931A-0015170AC494.root'
    )
# Input
process.source = cms.Source("PoolSource",
	fileNames = source_fileNames,
	#eventsToProcess = cms.untracked.VEventRange('123596:267558','123615:7692681','123615:968397','123596:1488220','123596:6732761')
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

### JES MET Corrections ########################################################
from JetMETCorrections.Type1MET.MetType1Corrections_cff import metJESCorAK5CaloJet

process.metMuonJESCorAK5 = metJESCorAK5CaloJet.clone()
process.metMuonJESCorAK5.inputUncorJetsLabel = recoJet_src
process.metMuonJESCorAK5.corrector = "L2L3JetCorrectorAK5Calo"
process.metMuonJESCorAK5.inputUncorMetLabel = "corMetGlobalMuons"
process.metCorSequence = cms.Sequence(process.metMuonJESCorAK5)

### PF+PAT #####################################################################
process.load("PhysicsTools.PatAlgos.patSequences_cff")
process.content = cms.EDAnalyzer("EventContentAnalyzer")
# Configure PAT to use PF2PAT instead of AOD sources
from PhysicsTools.PatAlgos.tools.pfTools import *
# Just to make it work...
process.out = cms.OutputModule("PoolOutputModule",outputCommands=cms.untracked.vstring())
usePF2PAT(process,runPF2PAT=True, jetAlgo='AK5') 
# Configure PAT to use AK5 jet collection	
from PhysicsTools.PatAlgos.tools.jetTools import *
switchJetCollection(process,cms.InputTag(recoJet_src),
                    doJTA        = True,
                    doBTagging   = True,
                    jetCorrLabel = ('AK5','PF'),
                    doType1MET   = False,
                    genJetCollection=cms.InputTag(genJet_src)
                    )
process.ak5JetID.src = cms.InputTag(recoJet_src)                

if runon=='data':
    process.patDefaultSequence.remove(process.genForPF2PATSequence)
    removeMCMatching(process, 'All')
    process.allLayer1Electrons.embedGenMatch = False
    process.allLayer1Muons.embedGenMatch = False

	# get the 900 GeV or 2360 GeV jet corrections
    switchJECSet( process, "900GeV")
	# switchJECSet( process, "2360GeV")

	# remove the tag infos
    process.allLayer1Jets.addTagInfos = False
	# require jet pt > 10 (L2+L3 corrected)
    process.selectedLayer1Jets.cut = cms.string('pt > 10')
	# look for at least n jets
	# process.countLayer1Jets.minNumber = 0

	# add the trigger information to the configuration
    from PhysicsTools.PatAlgos.tools.trigTools import *
    switchOnTrigger( process )
	# remove the dtrigger matches
    process.patTriggerSequence.remove( process.patTriggerMatcher )
    process.patTriggerEvent.patTriggerMatches = []
else:
    process.allLayer1Electrons.embedGenMatch = True
    process.allLayer1Muons.embedGenMatch = True

### Egamma Isolation ###########################################################
# Keep in line with non-PAT config.
process.eleIsoDepositEcalFromHits.ExtractorPSet.barrelEcalHits = "ecalRecHit:EcalRecHitsEB"
process.eleIsoDepositEcalFromHits.ExtractorPSet.endcapEcalHits = "ecalRecHit:EcalRecHitsEE"
process.eleIsoFromDepsEcalFromHitsByCrystal.deposits[0].vetos =	cms.vstring('EcalBarrel:0.045', 
																			'EcalBarrel:RectangularEtaPhiVeto(-0.02,0.02,-0.5,0.5)',
																			'EcalBarrel:AbsThresholdFromTransverse(0.08)',
																			'EcalEndcaps:AbsThreshold(0.100)',
																			'EcalEndcaps:0.070', 
																			'EcalEndcaps:RectangularEtaPhiVeto(-0.02,0.02,-0.5,0.5)')

# Example configuration
process.eleIsoFromDepsEcalFromHitsByCrystal.deposits[0].deltaR = 0.4

### Analysis configuration #####################################################
process.load("DiLeptonAnalysis.NTupleProducer.ntupleproducer_cfi")
process.analyze.tag_muons     = 'allLayer1Muons'
process.analyze.tag_electrons = 'allLayer1Electrons'
process.analyze.tag_jets      = 'allLayer1Jets'
# use 'photons' instead of 'allLayer1Photons' (photons are not supported in PF+PAT)
process.analyze.tag_photons	  = 'photons' 
process.analyze.isPat		= True
process.analyze.isRealData	= cms.untracked.bool(runon=='data')

#### Path ######################################################################
process.mybtag = cms.Sequence(   process.impactParameterTagInfos
                               * process.simpleSecondaryVertexBJetTags )
process.p = cms.Path(
    process.ak5JetID
    + process.patDefaultSequence
    + process.metCorSequence
    + process.mybtag
    + process.analyze
    )

