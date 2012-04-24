from CommonTools.ParticleFlow.pfNoPileUp_cff import *
pfPileUp.PFCandidates = "particleFlow"
pfNoPileUp.bottomCollection = "particleFlow"

pfPUSequence = cms.Sequence( pfPileUp * pfNoPileUp )

from MyAnalysis.IsolationTools.electronPFIsolations_cff import *
from MyAnalysis.IsolationTools.muonPFIsolations_cff import *
from MyAnalysis.IsolationTools.electronRadialPFIsolations_cff import *
from MyAnalysis.IsolationTools.muonRadialPFIsolations_cff import *

        
myPFIsolation = cms.Sequence( electronPFIsoChHad02 * electronPFIsoNHad02 * electronPFIsoPhoton02
                              * muonPFIsoChHad02 * muonPFIsoNHad02 * muonPFIsoPhoton02 
                              * electronPFIsoChHad03 * electronPFIsoNHad03 * electronPFIsoPhoton03
                              * muonPFIsoChHad03 * muonPFIsoNHad03 * muonPFIsoPhoton03 
                              * electronPFIsoChHad04 * electronPFIsoNHad04 * electronPFIsoPhoton04
                              * muonPFIsoChHad04 * muonPFIsoNHad04 * muonPFIsoPhoton04 )
        
myPFDirIsolation = cms.Sequence( electronRadPFIsoChHad02 * electronRadPFIsoNHad02 * electronRadPFIsoPhoton02
                                 * muonRadPFIsoChHad02 * muonRadPFIsoNHad02 * muonRadPFIsoPhoton02 
                                 * electronRadPFIsoChHad03 * electronRadPFIsoNHad03 * electronRadPFIsoPhoton03
                                 * muonRadPFIsoChHad03 * muonRadPFIsoNHad03 * muonRadPFIsoPhoton03 
                                 * electronRadPFIsoChHad04 * electronRadPFIsoNHad04 * electronRadPFIsoPhoton04
                                 * muonRadPFIsoChHad04 * muonRadPFIsoNHad04 * muonRadPFIsoPhoton04 )
        
pfIsolationAllSequence = cms.Sequence( pfPUSequence * myPFIsolation * myPFDirIsolation )
