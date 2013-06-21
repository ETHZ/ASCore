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

mygitaddpkg RecoLuminosity/LumiDB RecoLuminosity-LumiDB-V03-00-00
mygitaddpkg MuonAnalysis/MuonAssociators MuonAnalysis-MuonAssociators-V01-13-00 
mygitaddpkg RecoJets/Configuration RecoJets-Configuration-V02-04-17
mygitaddpkg JetMETAnalysis/ecalDeadCellTools Colin_TaggingMode_June30
mygitaddpkg RecoEcal/EgammaCoreTools RecoEcal-EgammaCoreTools-V05-07-00
mygitaddpkg RecoParticleFlow/PFClusterTools RecoParticleFlow-PFClusterTools-V12-01-01
mygitaddpkg RecoTauTag/RecoTau RecoTauTag-RecoTau-V01-02-07-02
mygitaddpkg RecoTauTag/TauTagTools RecoTauTag-TauTagTools-V01-02-00
mygitaddpkg RecoTauTag/Configuration RecoTauTag-Configuration-V01-02-09

mygitaddpkg h2gglobe vertex_mva_v4
mv h2gglobe/VertexAnalysis .
rm -rf h2gglobe/*
mv VertexAnalysis h2gglobe/VertexAnalysis

mygitaddpkg SCFootprintRemoval V00-02d
mkdir PFIsolation
mv SCFootprintRemoval PFIsolation/SuperClusterFootprintRemoval

mygitaddpkg Patches-OldReleases master
cp Patches-OldReleases/RecoEcal/EgammaCoreTools/plugins/EcalClusterLocalContCorrection.h RecoEcal/EgammaCoreTools/plugins/EcalClusterLocalContCorrection.h
cp Patches-OldReleases/RecoEcal/EgammaCoreTools/plugins/EcalClusterLocalContCorrection.cc RecoEcal/EgammaCoreTools/plugins/EcalClusterLocalContCorrection.cc
mkdir -p SandBox/Skims
cp Patches-OldReleases/SandBox/Skims/plugins/BuildFile.xml SandBox/Skims/plugins/BuildFile.xml
cp Patches-OldReleases/SandBox/Skims/plugins/RecovRecHitFilter.cc SandBox/Skims/plugins/RecovRecHitFilter.cc
cp Patches-OldReleases/SandBox/Skims/python/recovRecHitFilter_cfi.py SandBox/Skims/python/recovRecHitFilter_cfi.py
rm -rf Patches-OldReleases

mygitaddpkg ASCore branch_42X_noEDM
mkdir DiLeptonAnalysis
mv ASCore DiLeptonAnalysis/NTupleProducer
