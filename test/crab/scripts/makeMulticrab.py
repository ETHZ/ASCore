#!/usr/bin/python

## just put the samples in the right category and hit ./makeMulticrab.py
## perhaps check the number of jobs per file later on

## RARES
rareSamples = [
'/TTZJets_8TeV-madgraph_v2/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',
'/TTWJets_8TeV-madgraph/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',
'/TTWWJets_8TeV-madgraph/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',
'/ZZZNoGstarJets_8TeV-madgraph/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',
'/WWWJets_8TeV-madgraph/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',
'/WZJetsTo3LNu_TuneZ2_8TeV-madgraph-tauola/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',
'/WWGJets_8TeV-madgraph/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',
'/WWZNoGstarJets_8TeV-madgraph/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',
'/WZZNoGstarJets_8TeV-madgraph/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',
'/WWJetsTo2L2Nu_TuneZ2star_8TeV-madgraph-tauola/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',
]

## EWK
ewkSamples = [
'/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',
'/DYJetsToLL_M-10To50filter_8TeV-madgraph/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',
'/WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',
]

## TOP
topSamples = [
'/TT_CT10_TuneZ2star_8TeV-powheg-tauola/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',
'/TT_CT10_TuneZ2star_8TeV-powheg-tauola/Summer12_DR53X-PU_S10_START53_V7A-v2/AODSIM',
'/T_s-channel_TuneZ2star_8TeV-powheg-tauola/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',
'/T_t-channel_TuneZ2star_8TeV-powheg-tauola/Summer12_DR53X-PU_S10_START53_V7A-v3/AODSIM',
'/T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',

'/Tbar_s-channel_TuneZ2star_8TeV-powheg-tauola/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',
'/Tbar_t-channel_TuneZ2star_8TeV-powheg-tauola/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',
'/Tbar_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',
]

## QCD
qcdSamples = [
''
]

f = open('multicrab.cfg', 'w')


f. write(
'''
################################################################################
# Section for multicrab: now has just the template crab.cfg, but more
# keys might appear in the future
[MULTICRAB]
cfg=crab_ntupleproducer_mc_T2.cfg
################################################################################

################################################################################
# Section in common for all dataset
# General idea: you define all the parameter in the template (crab.cfg), 
# but you might want to change the template values for all dataset.
# The general syntax is that you first put the crab.cfg [SECTION] and
# the the crab.cfg [key], with a "." in between, exactly as you would do
# to pass to CRAB keys via command line.
# Any parameters can be set or changed
[COMMON]
GRID.se_black_list = T3
GRID.ce_black_list = T3

CMSSW.total_number_of_events=-1
#CMSSW.number_of_jobs = 450
#CMSSW.events_per_job = 25000

# Add a section for each dataset you want to access (or, more precisely,
# any task you want to create).
# The name of the section will be used as USER.ui_working_dir

#######################
##### RARE SAMPLES ####
#######################
'''
)

for sample in rareSamples:
	f.write('['+sample.lstrip('/').rstrip('/AODSIM').replace('/','-')+']\n')
	f.write('CMSSW.datasetpath='+sample+'\n')
	f.write('CMSSW.number_of_jobs = 300\n')
	f.write('\n')
	
f.write('''

#######################
##### EWK  SAMPLES ####
#######################

'''
)

for sample in ewkSamples:
	f.write('['+sample.lstrip('/').rstrip('/AODSIM').replace('/','-')+']\n')
	f.write('CMSSW.datasetpath='+sample+'\n')
	f.write('CMSSW.number_of_jobs = 1200\n')
	f.write('\n')
	

f.write('''

#######################
##### TOP  SAMPLES ####
#######################

'''
)

for sample in topSamples:
	f.write('['+sample.lstrip('/').rstrip('/AODSIM').replace('/','-')+']\n')
	f.write('CMSSW.datasetpath='+sample+'\n')
	f.write('CMSSW.number_of_jobs = 300\n')
	f.write('\n')

f.write('''

#######################
##### QCD  SAMPLES ####
#######################

'''
)

for sample in qcdSamples:
        f.write('['+sample.lstrip('/').rstrip('/AODSIM').replace('/','-')+']\n')
        f.write('CMSSW.datasetpath='+sample+'\n')
        f.write('CMSSW.number_of_jobs = 300\n')
        f.write('\n')


f.close()
