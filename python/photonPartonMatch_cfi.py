import FWCore.ParameterSet.Config as cms

from RecoJets.JetProducers.ak5GenJets_cfi import *

myphotons = cms.EDFilter("CandViewShallowCloneProducer",
  src = cms.InputTag("genParticles"),
  cut = cms.string(" pdgId==22 && mother(0).status==3")
)

# filter to request at least 1 photon with status 3 mother
#countFilter = cms.EDFilter("CandCountFilter",
#	src = cms.InputTag("photons"),
#	minNumber = cms.uint32(1)
#)

partons = cms.EDFilter("CandViewShallowCloneProducer",
  src = cms.InputTag("genParticles"),
  cut = cms.string(" (pdgId==21 || abs(pdgId)<6) && status==2 && pt>0 && numberOfMothers>0 " )
)


partonGenJets = ak5GenJets.clone()
partonGenJets.src = cms.InputTag("partons")

myPhotonJetMatch = cms.EDProducer("CandViewCombiner",
		decay       = cms.string("partonGenJets myphotons"),
		cut         = cms.string("deltaR(daughter(0).eta,daughter(0).phi,daughter(1).eta,daughter(1).phi) < 0.5 ")
)


photonPartonMatch = cms.Sequence(
	myphotons
#	*countFilter
	*partons
	*partonGenJets*myPhotonJetMatch
						 )

