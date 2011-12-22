import FWCore.ParameterSet.VarParsing as VarParsing


process = cms.Process("ParticleListDrawer")
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")


options = VarParsing.VarParsing ('standard') # set 'standard'  options
options.register ('nEvents', # register 'runon' option
                  1,  # the default value
                  VarParsing.VarParsing.multiplicity.singleton, # singleton or list
                  VarParsing.VarParsing.varType.int,         # string, int, or float
                  "Number of events to run over (standard:10)")
options.register ('SkipOver', # register 'runon' option
                  0,  # the default value
                  VarParsing.VarParsing.multiplicity.singleton, # singleton or list
                  VarParsing.VarParsing.varType.int,         # string, int, or float
                  "Skip the first ... events (standard:0)")
options.parseArguments()


process.printTree = cms.EDAnalyzer("ParticleListDrawer",
	maxEventsToPrint = cms.untracked.int32(options.nEvents),
	printVertex = cms.untracked.bool(False),
	src = cms.InputTag("genParticles")
)
process.source = cms.Source("PoolSource",
#      fileNames = cms.untracked.vstring('file:/shome/buchmann/smsproblem___7AB88B4C-BB1F-E111-BCA6-D485645943AC.root'),
      fileNames = cms.untracked.vstring('dcap://t3se01.psi.ch:22125/pnfs/psi.ch/cms/trivcat/store/user///fronga/0004D5B4-2AF1-E011-8AA6-001D09675288.root'),
      skipEvents = cms.untracked.uint32(options.SkipOver)
)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.nEvents) )
process.p = cms.Path(
       	process.printTree
   )

