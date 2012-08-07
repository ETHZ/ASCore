#! /bin/bash

export SCRAM_ARCH=slc5_amd64_gcc462
eval `scramv1 ru -sh`

cd $CMSSW_BASE/src/

#####################################################################################################################################
################################################                      ###############################################################
################################################     <recipe>         ###############################################################
################################################                      ###############################################################
#####################################################################################################################################



# Pat: from SWGuidePATRecipes
# tags from SWGuidePATReleaseNotes V08-09-11-00
addpkg DataFormats/PatCandidates   V06-05-01 
addpkg PhysicsTools/PatAlgos       V08-09-21
addpkg CommonTools/ParticleFlow    V00-03-15

# revert back to the tag we proposed originally and is in the official recipe:
addpkg RecoParticleFlow/PFProducer V15-01-11

# Type1MET
addpkg JetMETCorrections/Type1MET  V04-06-09
addpkg PhysicsTools/PatUtils V03-09-23
addpkg CommonTools/RecoUtils V00-00-12

# Updated Tau discriminators (not including advanced MVA isolation: breaks CHS postfix)
cvs co -r V01-04-17 RecoTauTag/RecoTau       
cvs co -r V01-04-03 RecoTauTag/Configuration 
cvs co -r V00-04-01 CondFormats/EgammaObjects

#  MET Filters
cvs co -r V00-00-08      RecoMET/METAnalyzers                             
cvs co -r V00-00-07      RecoMET/METFilters

# ECAL
addpkg RecoEcal/EgammaCoreTools V05-08-22
cvs export -d tmpexportdir -rHEAD UserCode/peruzzi/RecoEcal/EgammaCoreTools/plugins/EcalClusterLocalContCorrection.h
cvs export -d tmpexportdir -rHEAD UserCode/peruzzi/RecoEcal/EgammaCoreTools/plugins/EcalClusterLocalContCorrection.cc
mv -f tmpexportdir/EcalClusterLocalContCorrection.h RecoEcal/EgammaCoreTools/plugins
mv -f tmpexportdir/EcalClusterLocalContCorrection.cc RecoEcal/EgammaCoreTools/plugins
rmdir tmpexportdir
cvs co -r vertex_mva_v4 -d h2gglobe/VertexAnalysis UserCode/HiggsAnalysis/HiggsTo2photons/h2gglobe/VertexAnalysis

# Isolation
cvs co -r V04-00-01 -d MyAnalysis/IsolationTools UserCode/emanuele/MyAnalysis/IsolationTools
# alternate code for photon isolation
# https://twiki.cern.ch/twiki/bin/view/CMS/EgammaPFBasedIsolation#Alternate_code_to_calculate_PF_I
cvs co -r V00-00-21 -d EGamma/EGammaAnalysisTools UserCode/EGamma/EGammaAnalysisTools
cvs up -r 1.13 EGamma/EGammaAnalysisTools/interface/PFIsolationEstimator.h
cvs up -r 1.20 EGamma/EGammaAnalysisTools/src/PFIsolationEstimator.cc

#####################################################################################################################################
################################################                      ###############################################################
################################################     </recipe>        ###############################################################
################################################                      ###############################################################
#####################################################################################################################################

echo "Everything has been set up. You can compile now (scramv1 b -j2) or modify the recipe to your likings"

exit 0
