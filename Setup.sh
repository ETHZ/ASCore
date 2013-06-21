mygitaddpkg()
{
REMOTE=git://github.com/ETHZ
PACKNAME=$1
VERSION=$2
PACKNAME_GIT=$(echo ${PACKNAME} | sed 's,/,-,g').git

echo Checking out from the ETHZ github area the package ${PACKNAME} version ${VERSION}
git clone ${REMOTE}/${PACKNAME_GIT} ${PACKNAME}
cd ${PACKNAME}
git checkout --detach ${VERSION}
cd ${OLDPWD}
}

# Pat: from SWGuidePATRecipes
# tags from SWGuidePATReleaseNotes V08-09-11-00
mygitaddpkg DataFormats/PatCandidates   DataFormats-PatCandidates-V06-05-06-03
mygitaddpkg PhysicsTools/PatAlgos       PhysicsTools-PatAlgos-V08-09-47
mygitaddpkg CommonTools/ParticleFlow    CommonTools-ParticleFlow-V00-03-16

# revert back to the tag we proposed originally and is in the official recipe:
mygitaddpkg RecoParticleFlow/PFProducer RecoParticleFlow-PFProducer-V15-02-06

# Type1MET
mygitaddpkg JetMETCorrections/Type1MET  JetMETCorrections-Type1MET-V04-06-09
mygitaddpkg PhysicsTools/PatUtils PhysicsTools-PatUtils-V03-09-26
mygitaddpkg CommonTools/RecoUtils CommonTools-RecoUtils-V00-00-13
mygitaddpkg DataFormats/StdDictionaries DataFormats-StdDictionaries-V00-02-15 

# Updated Tau discriminators (not including advanced MVA isolation: breaks CHS postfix)
mygitaddpkg  RecoTauTag/RecoTau         RecoTauTag-RecoTau-V01-04-23
mygitaddpkg  RecoTauTag/Configuration 	RecoTauTag-Configuration-V01-04-10
mygitaddpkg  CondFormats/EgammaObjects	CondFormats-EgammaObjects-V00-04-01

#needed at the moment for AgainstMuonXXX2 --> will be in PhysicsTools/PatAlgos       V08-09-51
cd PhysicsTools/PatAlgos/python/producersLayer1
git checkout PhysicsTools-PatAlgos-V08-09-51 tauProducer_cfi.py
cd ${OLDPWD}
cd PhysicsTools/PatAlgos/python/tools
git checkout PhysicsTools-PatAlgos-V08-09-51 tauTools.py
cd ${OLDPWD}

#  MET Filters (including tracking POG filters)
mygitaddpkg    RecoMET/METAnalyzers                     V00-00-08   
mygitaddpkg    RecoMET/METFilters			RecoMET-METFilters-V00-00-13   
mygitaddpkg    CommonTools/RecoAlgos			CommonTools-RecoAlgos-V00-03-23   
mygitaddpkg    DPGAnalysis/Skims			DPGAnalysis-Skims-V01-00-11-01
mygitaddpkg    DPGAnalysis/SiStripTools			DPGAnalysis-SiStripTools-V00-11-17   
mygitaddpkg    DataFormats/TrackerCommon		DataFormats-TrackerCommon-V00-00-08   
mygitaddpkg    RecoLocalTracker/SubCollectionProducers	RecoLocalTracker-SubCollectionProducers-V01-09-05   
mygitaddpkg    EventFilter/HcalRawToDigi		EventFilter-HcalRawToDigi-V01-02-10   

# Parton Flavour
mygitaddpkg PhysicsTools/JetMCAlgos PhysicsTools-JetMCAlgos-V00-13-10

# ECAL
mygitaddpkg RecoEcal/EgammaCoreTools RecoEcal-EgammaCoreTools-V05-08-26
mygitaddpkg Patches-OldReleases master
cp Patches-OldReleases/RecoEcal/EgammaCoreTools/plugins/EcalClusterLocalContCorrection.h RecoEcal/EgammaCoreTools/plugins/EcalClusterLocalContCorrection.h
cp Patches-OldReleases/RecoEcal/EgammaCoreTools/plugins/EcalClusterLocalContCorrection.cc RecoEcal/EgammaCoreTools/plugins/EcalClusterLocalContCorrection.cc
rm -rf Patches-OldReleases

mygitaddpkg h2gglobe vertex_mva_v4
mv h2gglobe/VertexAnalysis .
rm -rf h2gglobe/*
mv VertexAnalysis h2gglobe/VertexAnalysis

mygitaddpkg emanuele/MyAnalysis/IsolationTools V04-00-01
mv emanuele/MyAnalysis .
rm -rf emanuele
mygitaddpkg sixie/Muon/MuonAnalysisTools V00-00-10
mv sixie/Muon .
rm -rf sixie

mygitaddpkg EGamma/EGammaAnalysisTools V00-00-21
cd EGamma/EGammaAnalysisTools/interface
git checkout V00-00-22 PFIsolationEstimator.h
cd ${OLDPWD}
cd EGamma/EGammaAnalysisTools/src
git checkout V00-00-22 PFIsolationEstimator.cc
cd ${OLDPWD}


mygitaddpkg SCFootprintRemoval V00-02d
mkdir PFIsolation
mv SCFootprintRemoval PFIsolation/SuperClusterFootprintRemoval

mygitaddpkg ASCore EDMdev
mkdir DiLeptonAnalysis
mv ASCore DiLeptonAnalysis/NTupleProducer


