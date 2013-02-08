#!/usr/bin/python

## just put the samples in the right category and hit ./makeMulticrab.py
## perhaps check the number of jobs per file later on

## RARES
rareSamples = [
#ttX
'/TTZJets_8TeV-madgraph_v2/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',
'/TTWJets_8TeV-madgraph/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',
'/TTGJets_8TeV-madgraph/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',
'/TTH_Inclusive_M-125_8TeV_pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',
'/TTWWJets_8TeV-madgraph/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',
'/TBZToLL_4F_TuneZ2star_8TeV-madgraph-tauola/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM',
##triboson
'/WWWJets_8TeV-madgraph/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',
'/WWZNoGstarJets_8TeV-madgraph/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',
'/WZZNoGstarJets_8TeV-madgraph/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',
'/ZZZNoGstarJets_8TeV-madgraph/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',
'/WWGJets_8TeV-madgraph/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',
#diboson
'/WWJetsTo2L2Nu_TuneZ2star_8TeV-madgraph-tauola/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',
'/ZZJetsTo2L2Q_TuneZ2star_8TeV-madgraph-tauola/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',
'/ZZJetsTo4L_TuneZ2star_8TeV-madgraph-tauola/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',
'/ZZJetsTo2L2Nu_TuneZ2star_8TeV-madgraph-tauola/Summer12_DR53X-PU_S10_START53_V7A-v3/AODSIM',
#WZ
'/WZJetsTo2L2Q_TuneZ2star_8TeV-madgraph-tauola/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',
'/WZJetsTo3LNu_TuneZ2_8TeV-madgraph-tauola/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',
#same sign W
'/WpWpqq_8TeV-madgraph/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',
'/WmWmqq_8TeV-madgraph/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',
'/WW_DoubleScattering_8TeV-pythia8/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',
#others
'/GVJets_Incl_8TeV-madgraph/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM',
'/GJets_HT-200To400_8TeV-madgraph_v2/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',
'/GJets_HT-400ToInf_8TeV-madgraph_v3/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM',
]

## EWK
ewkSamples = [
'/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',
'/DYJetsToLL_M-10To50filter_8TeV-madgraph/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',
'/WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',
'/WJetsToLNu_HT-200To250_8TeV-madgraph/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM',
'/WJetsToLNu_HT-250To300_8TeV-madgraph/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',
'/WJetsToLNu_HT-300To400_8TeV-madgraph/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',
'/WJetsToLNu_HT-400ToInf_8TeV-madgraph/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',
'/WJetsToLNu_HT-250To300_8TeV-madgraph_v2/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',
'/WJetsToLNu_HT-300To400_8TeV-madgraph_v2/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',
'/WJetsToLNu_HT-400ToInf_8TeV-madgraph_v2/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',
'/ZJetsToNuNu_100_HT_200_TuneZ2Star_8TeV_madgraph/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM',
'/ZJetsToNuNu_200_HT_400_TuneZ2Star_8TeV_madgraph/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',
'/ZJetsToNuNu_400_HT_inf_TuneZ2Star_8TeV_madgraph/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',
'/ZJetsToNuNu_100_HT_200_TuneZ2Star_8TeV_madgraph_ext/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM',
'/ZJetsToNuNu_200_HT_400_TuneZ2Star_8TeV_madgraph_ext/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',
'/ZJetsToNuNu_400_HT_inf_TuneZ2Star_8TeV_madgraph_ext/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',
]

## TOP
topSamples = [
'/TT_CT10_TuneZ2star_8TeV-powheg-tauola/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',
'/TT_CT10_TuneZ2star_8TeV-powheg-tauola/Summer12_DR53X-PU_S10_START53_V7A-v2/AODSIM',
#single top
'/T_s-channel_TuneZ2star_8TeV-powheg-tauola/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',
'/T_t-channel_TuneZ2star_8TeV-powheg-tauola/Summer12_DR53X-PU_S10_START53_V7A-v3/AODSIM',
'/T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',
#single tbar
'/Tbar_s-channel_TuneZ2star_8TeV-powheg-tauola/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',
'/Tbar_t-channel_TuneZ2star_8TeV-powheg-tauola/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',
'/Tbar_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',
]

## QCD
qcdSamples = [
'/QCD_Pt-50to80_TuneZ2star_8TeV_pythia6/Summer12_DR53X-PU_S10_START53_V7A-v2/AODSIM',
'/QCD_Pt-80to120_TuneZ2star_8TeV_pythia6/Summer12_DR53X-PU_S10_START53_V7A-v3/AODSIM',
'/QCD_Pt-120to170_TuneZ2star_8TeV_pythia6/Summer12_DR53X-PU_S10_START53_V7A-v3/AODSIM',
'/QCD_Pt-170to300_TuneZ2star_8TeV_pythia6/Summer12_DR53X-PU_S10_START53_V7A-v2/AODSIM',
'/QCD_Pt-170to300_TuneZ2star_8TeV_pythia6_v2/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',
'/QCD_Pt-300to470_TuneZ2star_8TeV_pythia6/Summer12_DR53X-PU_S10_START53_V7A-v2/AODSIM',
'/QCD_Pt-300to470_TuneZ2star_8TeV_pythia6_v2/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',
'/QCD_Pt-300to470_TuneZ2star_8TeV_pythia6_v3/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',
'/QCD_Pt-470to600_TuneZ2star_8TeV_pythia6/Summer12_DR53X-PU_S10_START53_V7A-v2/AODSIM',
'/QCD_Pt-600to800_TuneZ2star_8TeV_pythia6/Summer12_DR53X-PU_S10_START53_V7A-v2/AODSIM',
'/QCD_Pt-800to1000_TuneZ2star_8TeV_pythia6/Summer12_DR53X-PU_S10_START53_V7A-v2/AODSIM',
'/QCD_Pt-1000to1400_TuneZ2star_8TeV_pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',
'/QCD_Pt-1400to1800_TuneZ2star_8TeV_pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',
'/QCD_Pt-1800_TuneZ2star_8TeV_pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',
'/QCD_HT-100To250_TuneZ2star_8TeV-madgraph-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',
'/QCD_HT-250To500_TuneZ2star_8TeV-madgraph-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',
'/QCD_HT-500To1000_TuneZ2star_8TeV-madgraph-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',
'/QCD_HT-1000ToInf_TuneZ2star_8TeV-madgraph-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',
# lepton and b- enriched samples
'/QCD_Pt_20_MuEnrichedPt_15_TuneZ2star_8TeV_pythia6/Summer12_DR53X-PU_S10_START53_V7A-v3/AODSIM',
'/QCD_Pt_20_30_EMEnriched_TuneZ2star_8TeV_pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',
'/QCD_Pt_30_80_EMEnriched_TuneZ2star_8TeV_pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',
'/QCD_Pt_80_170_EMEnriched_TuneZ2star_8TeV_pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',
'/QCD_Pt_170_250_EMEnriched_TuneZ2star_8TeV_pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',
'/QCD_Pt_250_350_EMEnriched_TuneZ2star_8TeV_pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',
'/QCD_Pt_350_EMEnriched_TuneZ2star_8TeV_pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',
'/QCD_Pt-50To150_bEnriched_TuneZ2star_8TeV-pythia6-evtgen/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',
'/QCD_Pt-150_bEnriched_TuneZ2star_8TeV-pythia6-evtgen/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',
]


f = open('multicrab.cfg', 'w')

def writeDatasets(slist, name, njobs=1000):
	f.write('''

#######################
##### %s  SAMPLES ####
#######################

''' % (name)
	)

	for sample in slist:
		f.write('# ['+sample.lstrip('/').rstrip('/AODSIM').replace('/','-')+']\n')
		f.write('# CMSSW.datasetpath='+sample+'\n')
		f.write('# CMSSW.number_of_jobs = 1200\n')
		f.write('\n')



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

'''
)

writeDatasets(rareSamples, 'RARE',  300)
writeDatasets(ewkSamples , 'EWK' , 1200)
writeDatasets(topSamples , 'TOP' ,  600)
writeDatasets(qcdSamples , 'QCD' ,  600)


f.close()
