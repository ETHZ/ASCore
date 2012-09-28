import os, ROOT, helper, sys
from math import *

print 'This script produces the datacards for some regions and two models'

xsecFile = ROOT.TFile('C1N2_8TeV_xsecs.root', 'READ', 'xsecFile')
xsecs = xsecFile.Get('errors')


for region in ['HT0MET200lV','HT0MET120NJ2bVlV']:
	infileName  = '../SMSresults_' + region +'.root'
	canvasName  = 'exclusion_HT0_MET200lV_PT20'

	if region == 'HT0MET200lV' :
		chanID = 181
	elif region == 'HT0MET120NJ2bVlV':
		chanID = 182
	else :
		chanID = -1
		
	if chanID == -1 : print 'SOMETHING WENT WRONG!!!, have a look at the regions'
	# thirdVeto = False
	# if 'lV' in infileName:
	# 	thirdVeto = True
	# thirdString = ''
	# if thirdVeto:
	# 	thirdString+='_thirdVeto'
	# outTextName = 'textEfficiencies_HT0_MET200_paper'+thirdString+'.txt'
	
	# if not os.path.isfile(infileName):
	# 	sys.exit('not a file')
	
	# change to nice style for the histo
	mystyle = helper.set1DStyle()
	mystyle.cd()

	statLumi = 0.045
	muEff    = 0.1
	elEff    = 0.1

	cnt=0
	cnt_exp=0
	cnt_exphi=0
	cnt_explo=0

	print 'Running on region '+region
	
	for model in ['TChi','TChiRight']:
		print '... for the model'+model
		for x in ['5','95']:
			print '...with x'+x
			datacardfile = open('SS_'+region+'_'+model+'_datacard_x'+x+'.txt', 'w')
			datacardfile.write('#scan <m1> <m2> <channel id> <yield in pb> <total yield uncertainty> <uncorrelated uncertainty> <muon eff uncertainty>')
			datacardfile.write('<electron eff uncertainty> <tau eff uncertainty> <trigger eff uncertainty> <JES uncertainty> <theory uncertainties> \n')
			datacardfile.write('# Datacard generated for '+model+' and x'+x+'\n')
			datacardfile.write('# ==================================================================================================\n')

			xv = float(x)/100.
			resultFile = ROOT.TFile(infileName, 'READ', 'resultFile')
			resultFile.cd()

			### READING HISTOGRAMS FROM THE FILE TO GET UNCERTAINTIES
			effHistoMet     = resultFile.Get(model+'_eff_x'+x+'_metsmear')
			effHistoMetDown = resultFile.Get(model+'_eff_x'+x+'_METdown' )
			effHistoMetUp   = resultFile.Get(model+'_eff_x'+x+'_METup'   )
			effHistoLepDown = resultFile.Get(model+'_eff_x'+x+'_lepdown' )
			effHistoLepUp   = resultFile.Get(model+'_eff_x'+x+'_lepup'   )
			effHistoJER     = resultFile.Get(model+'_eff_x'+x+'_JER'     )
			effHistoJESUp   = resultFile.Get(model+'_eff_x'+x+'_JESup'   )
			effHistoJESDown = resultFile.Get(model+'_eff_x'+x+'_JESdown' )
			effHistoNorm    = resultFile.Get(model+'_eff_x'+x+'_norm'    )
			yieldHisto      = resultFile.Get(model+'_yield'+x)
			if model == 'TChi'     : countHisto      = resultFile.Get('ModelCount'+x)
			if model == 'TChiRight': countHisto      = resultFile.Get('RightHandedSlepCount'+x)
			nPassHisto      = resultFile.Get(model+'_nPass_x'+x+'_norm')

			# smoothing the limit (unphysical but pretty)
			nContSmooths = 0
			# smoothing the yields (physical but ugly)
			if model == 'TChiSlepSlep': nYieldSmooths = 0
			else : nYieldSmooths = 0

			#yield smoothing
			for n in range(nYieldSmooths):
				yieldHisto.Smooth()

			effErrStat     = ROOT.TH2D (model+x+'effErrStat'      , model+x+'effErrStat'      , 32, 0, 800, 32, 0, 800)
			effErrMetSmear = ROOT.TH2D (model+x+'effErrMetSmear'  , model+x+'effErrMetSmear'  , 32, 0, 800, 32, 0, 800)
			effErrMet      = ROOT.TH2D (model+x+'effErrMet'       , model+x+'effErrMet'       , 32, 0, 800, 32, 0, 800)
			effErrLep      = ROOT.TH2D (model+x+'effErrLep'       , model+x+'effErrLep'       , 32, 0, 800, 32, 0, 800)
			effErrJES      = ROOT.TH2D (model+x+'effErrJES'       , model+x+'effErrJES'       , 32, 0, 800, 32, 0, 800)
			effErrJER      = ROOT.TH2D (model+x+'effErrJER'       , model+x+'effErrJER'       , 32, 0, 800, 32, 0, 800)
			effErrTot      = ROOT.TH2D (model+x+'effErrTot'       , model+x+'effErrTot'       , 32, 0, 800, 32, 0, 800)
			effErrApprox   = ROOT.TH2D (model+x+'effErrApprox'    , model+x+'effErrApprox'    , 32, 0, 800, 32, 0, 800)

			# fill 1 in every bin which is excluded
			for mChi in range(0, 1050, 50):
				for mLSP in range(0, 1000, 50):
					if mLSP > mChi : continue
					bin = yieldHisto.FindBin(mChi, mLSP)
					binYield = yieldHisto.GetBinContent(bin)
					eff        = float(effHistoNorm.GetBinContent(bin))
					# statistical studies
					count      = int(countHisto.GetBinContent(bin))
					nPass      = float(nPassHisto.GetBinContent(bin))
					if not count: continue
					if not eff  : continue

					# systematics studies
					effMetSmear = float(effHistoMet.GetBinContent(bin))
					effMetUp    = float(effHistoMetUp.GetBinContent(bin))
					effMetDown  = float(effHistoMetDown.GetBinContent(bin))
					effLepUp    = float(effHistoLepUp.GetBinContent(bin))
					effLepDown  = float(effHistoLepDown.GetBinContent(bin))
					effJESUp    = float(effHistoJESUp.GetBinContent(bin))
					effJESDown  = float(effHistoJESDown.GetBinContent(bin))
					effJER      = float(effHistoJER.GetBinContent(bin))
					if eff == 0 or effMetSmear == 0 or effMetUp == 0 or effMetDown == 0 or effLepUp == 0 or effLepDown == 0 or effJESUp == 0 or effJESDown == 0 or effJER == 0: continue

					if abs(eff - effLepUp) > abs(eff - effLepDown): lepErr = effLepUp
					else: lepErr = effLepDown
					if abs(eff - effMetUp) > abs(eff - effMetDown): metErr = effMetUp
					else: metErr = effMetDown
					if abs(eff - effJESUp) > abs(eff - effJESDown): jesErr = effJESUp
					else: jesErr = effJESDown

					systErrMetS = abs(1 - effMetSmear / eff)
					systErrMet  = abs(1 - metErr      / eff)
					systErrLep  = abs(1 - lepErr      / eff)
					systErrJES  = abs(1 - jesErr      / eff)
					systErrJER  = abs(1 - effJER      / eff)
					systErrTau  = 0.
					systErrTrig = 0.
					systErrTh   = float(xsecs.GetBinContent(xsecs.FindBin(mChi)))

					systErrTot  = sqrt(systErrMet*systErrMet + systErrLep*systErrLep + systErrJES*systErrJES + systErrJER*systErrJER)
					statErrTot  = yieldHisto.GetBinError(bin)

					effErrMetSmear.Fill(mChi, mLSP, systErrMetS)
					effErrMet     .Fill(mChi, mLSP, systErrMet )
					effErrLep     .Fill(mChi, mLSP, systErrLep )
					effErrJES     .Fill(mChi, mLSP, systErrJES )
					effErrJER     .Fill(mChi, mLSP, systErrJER )
					effErrTot     .Fill(mChi, mLSP, systErrTot )
					effErrStat    .Fill(mChi, mLSP, statErrTot )

					weight = 1./9200.

					# total uncertainty from lepton and MET scaling as well as statistical
					newLine1 = 'scan %3.0f %3.0f %3d %6.5e %6.5e %6.5e '%(mChi, mLSP, chanID, binYield*weight, systErrTot*weight, statErrTot*weight)
					newLine2 = '%6.5e %6.5e %6.5e %6.5e %6.5e %6.5e\n'  %(systErrLep*weight, systErrLep*weight, systErrTau*weight, systErrTrig*weight, systErrJES*weight, systErrTh*weight)
					datacardfile.write(newLine1)
					datacardfile.write(newLine2)
			datacardfile.close()
	print 'DONE!!'
print '\n\n EVERYTHING is done... closing...'	
