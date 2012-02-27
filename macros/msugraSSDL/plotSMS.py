import os, math, ROOT, helper

# this is a python script to create nice exclusion plots for the TChiNuSlep simplified
# model. it is easily expandable to other models or the msugra model, although there
# already esists a c++ script for msugra.

# infileName is the root file with the efficiency and yield histogram named SMSeff and SMSyield
infileName  = '../SMSresults_HT0JV_MET120_PT20_10.root'
# the countour histograms will be saved in the outfile (not really necessary)
outfileName = 'SMScontours_HT0JV_MET120_PT20_10.root'
# the final canvas with the nice histogram will be saved in canvasName.root and canvasName.pdf
canvasName  = 'plots/SMSexclusion_HT0JV_MET120_PT20_10'
#if HT = -1 --> Jet-veto
HT  = -1
MET = 120
# the following limits are from the LandS tool, insert your numbers here
observed =  6.62886
expected =  8.99075
exphi    = 13.11700
explo    =  6.53056

#infileName  = '../SMSresults_HT0_MET200_PT20_10.root'
#outfileName = 'SMScontours_HT0_MET200_PT20_10.root'
#canvasName  = 'plots/SMSexclusion_HT0_MET200_PT20_10'
##if HT = -1 --> Jet-veto
#HT  = 0
#MET = 200
## the following limits are from the LandS tool, insert your numbers here
#observed =  8.41036
#expected =  8.99802
#exphi    = 13.13110
#explo    =  6.51938


# change to nice style for the histo
mystyle = helper.set1DStyle()
mystyle.cd()

# get the file with yield and efficiency
resultFile = ROOT.TFile(infileName, 'READ', 'resultFile')
resultFile.cd()
effHisto   = resultFile.Get('SMSeff')
yieldHisto = resultFile.Get('SMSyield')

# smoothing the limit (unphysical but pretty)
nContSmooths = 1
# smoothing the yields (physical but ugly)
nYieldSmooths = 0

#yield smoothing
for n in range(nYieldSmooths):
	yield_.Smooth()

# create new histogram for observed, expected, +- 1 sigma
obs_contour   = ROOT.TH2D ('obs_exclusion'   , 'obs_exclusion'   , 50, 0, 500, 50, 0, 500)
exp_contour   = ROOT.TH2D ('exp_exclusion'   , 'exp_exclusion'   , 50, 0, 500, 50, 0, 500)
explo_contour = ROOT.TH2D ('explo_exclusion' , 'explo_exclusion' , 50, 0, 500, 50, 0, 500)
exphi_contour = ROOT.TH2D ('exphi_exclusion' , 'exphi_exclusion' , 50, 0, 500, 50, 0, 500)
newEffHisto   = ROOT.TH2D ('niceEff'         , 'niceEff'         , 50, 0, 500, 50, 0, 500)

# fill 1 in every bin that is excluded
for mChi in range(0, 510, 10):
	for mLSP in range(0, 510, 10):
		bin = yieldHisto.FindBin(mChi, mLSP)
		binYield = yieldHisto.GetBinContent(bin)
		newEffHisto.Fill(mChi, mLSP, 100.*effHisto.GetBinContent(bin))
		if (binYield > observed): obs_contour.Fill(mChi, mLSP, 1)
		if (binYield > expected): exp_contour.Fill(mChi, mLSP, 1)
		if (binYield > exphi)   : exphi_contour.Fill(mChi, mLSP, 1)
		if (binYield > explo)   : explo_contour.Fill(mChi, mLSP, 1)

# this is for filling the holes in the plot
for mChi in range(0, 510, 10):
	for mLSP in range(0, 510, 10):
		if mLSP > mChi : continue
		bin = newEffHisto.FindBin(mChi, mLSP)
		binLeft  = newEffHisto.FindBin(mChi-10, mLSP)
		binRight = newEffHisto.FindBin(mChi+10, mLSP)
		binDown  = newEffHisto.FindBin(mChi, mLSP-10)
		binUp    = newEffHisto.FindBin(mChi, mLSP+10)
		binContent      = newEffHisto.GetBinContent(bin)
		binContentLeft  = newEffHisto.GetBinContent(binLeft)
		binContentRight = newEffHisto.GetBinContent(binRight)
		binContentDown  = newEffHisto.GetBinContent(binDown)
		binContentUp    = newEffHisto.GetBinContent(binUp)
		average = (binContentLeft + binContentRight + binContentDown + binContentUp ) /4.
		if binContent == 0:
			if (binContentLeft != 0. and binContentRight != 0. and binContentUp != 0. and binContentDown != 0.):
				print 'Bin at mChi:', mChi, 'and mLSP:', mLSP
				newEffHisto.Fill(mChi, mLSP, average)

# write the histograms into a file, not really necessary if you trust the other part of the code
outfile = ROOT.TFile(outfileName, 'RECREATE', '')
outfile.cd()
obs_contour.Write()
exp_contour.Write()
explo_contour.Write()
exphi_contour.Write()

for i in range(nContSmooths):
	obs_contour  .Smooth()
	exp_contour  .Smooth()
	explo_contour.Smooth()
	exphi_contour.Smooth()

# get the contour graphs from the histograms
grph       = helper.getContourGraph(obs_contour)
grph_exp   = helper.getContourGraph(exp_contour)
grph_explo = helper.getContourGraph(explo_contour)
grph_exphi = helper.getContourGraph(exphi_contour)

# cosmetics
grph.SetLineColor(2)
grph.SetLineWidth(3)

grph_exp.SetLineColor(1)
grph_exp.SetLineWidth(3)

grph_explo.SetLineColor(1)
grph_explo.SetLineWidth(2)
grph_explo.SetLineStyle(2)

grph_exphi.SetLineColor(1)
grph_exphi.SetLineWidth(2)
grph_exphi.SetLineStyle(2)

legend = helper.makeLegend(0.17, 0.65, 0.38, 0.85 )
legend.AddEntry(grph       , 'Observed limit'       , 'l')
legend.AddEntry(grph_exp   , 'Expected limit'       , 'l')
legend.AddEntry(grph_exphi , 'Expected #pm 1 #sigma', 'l')

## this would create a shaded error between the plus and minus one sigma bands
#if grph_exphi.GetN() >= grph_explo.GetN(): n = grph_explo.GetN()
#else: n = grph_exphi.GetN()
#shade = ROOT.TGraph(2*n)
#for i in range(n):
#	shade.SetPoint(i   , grph_explo.GetX()[i]     , grph_explo.GetY()[i])
#	shade.SetPoint(n+i , grph_explo.GetX()[n-i-1] , grph_exphi.GetY()[n-i-1])
#
#shade.SetFillStyle(3003)
#shade.SetFillColor(39)
#shade.Draw('f')

ROOT.gStyle.SetPadRightMargin(0.18)
ROOT.gStyle.SetPadLeftMargin(0.13)
ROOT.gStyle.SetPadBottomMargin(0.13)
helper.useNiceColorPalette()
canv = ROOT.TCanvas('bla', 'bla', 800, 600)
canv.cd()

newEffHisto.GetXaxis().SetTitle('m_{#chi^{#pm}_{2}} (GeV)')
newEffHisto.GetXaxis().SetNdivisions(505)
newEffHisto.GetYaxis().SetTitle('m_{LSP} (GeV)')
newEffHisto.GetYaxis().SetNdivisions(505)
newEffHisto.GetZaxis().SetTitle('Signal eff. #times acc. (%)')
newEffHisto.GetZaxis().SetTitleOffset(0.8)
newEffHisto.GetZaxis().SetRangeUser(0., 16.)
newEffHisto.Draw('colz')
grph.Draw('same')
grph_exp.Draw('same')
grph_explo.Draw('same')
grph_exphi.Draw('same')
legend.Draw('same')

if HT == 0 : helper.drawLatex('#splitline{E_{T}^{miss} > %3.0f GeV}{N_{Jets} #geq 0}' %(MET)        , 0.45 , 0.70 , 0.04, 62)
if HT == -1: helper.drawLatex('#splitline{E_{T}^{miss} > %3.0f GeV}{N_{Jets} = 0 }' %(MET)        , 0.45 , 0.70 , 0.04, 62)
helper.drawLatex('#tilde#chi^{#pm}_{1} - #tilde#chi^{0}_{2} Production'               , 0.17 , 0.48 , 0.03)
helper.drawLatex('m_{#tilde#chi^{#pm}_{1}} = m_{#tilde#chi^{0}_{2}}'                  , 0.17 , 0.44 , 0.03)
helper.drawLatex('m_{#tildel} = m_{#tilde#nu}'                                        , 0.17 , 0.40 , 0.03)
helper.drawLatex('m_{#tildel} = 0.5#upoint(#tilde#chi^{0}_{1} - m_{LSP})'             , 0.17 , 0.36 , 0.03)
helper.drawLatex('N_{UL}^{obs} (N_{UL}^{exp}) = %4.2f (%4.2f)' %(observed , expected) , 0.45 , 0.80 , 0.035)
helper.drawTopLine(lumi = 4680*1.06, False)

# save the canvas as .root and .pdf
canv.SaveAs(canvasName+'.pdf', '')
canv.SaveAs(canvasName+'.root', '')
