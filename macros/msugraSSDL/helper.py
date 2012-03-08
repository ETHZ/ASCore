import ROOT, commands
from array import array

def getLimits():
	limits = {}
	for uncert in [5, 10, 15, 20, 25, 30, 35, 40, 45]:
		combString = './combine card_HT0_MET200_sigUncert'+str(uncert)+'.txt '
		fileString = '-M HybridNew --freq --grid=/scratch/mdunser/Limits/Mar6/limit_HT0_MET200_sigUncert'+str(uncert)+'.root'
		limits[uncert] = {}
		limits[uncert]['obs']   = commands.getoutput(combString+fileString)                           .split('\n')[-2].split()[3]
		limits[uncert]['exp']   = commands.getoutput(combString+fileString+' --expectedFromGrid 0.50').split('\n')[-2].split()[3]
		limits[uncert]['exphi'] = commands.getoutput(combString+fileString+' --expectedFromGrid 0.84').split('\n')[-2].split()[3]
		limits[uncert]['explo'] = commands.getoutput(combString+fileString+' --expectedFromGrid 0.16').split('\n')[-2].split()[3]
	return limits
	

def getContourGraph(histo, n=10):
	histo.Draw('cont list')
	ROOT.gPad.Update()
	object = ROOT.gROOT.GetListOfSpecials().FindObject('contours')
	list   = object.At(n)
	graph  = list.First().Clone()
	return graph


def useNiceColorPalette( NCont = 999):
	stops = [0.00, 0.34, 0.61, 0.84, 1.00]
	red   = [0.00, 0.00, 0.87, 1.00, 0.51]
	green = [0.00, 0.81, 1.00, 0.20, 0.00]
	blue  = [0.51, 1.00, 0.12, 0.00, 0.00]
	
	s = array('d', stops)
	r = array('d', red)
	g = array('d', green)
	b = array('d', blue)
	
	nstops = len(s)
	ROOT.TColor.CreateGradientColorTable(nstops, s, r, g, b, NCont)
	ROOT.gStyle.SetNumberContours( NCont )

def drawLatex(string, x, y, size = 0.04, font=42):
	lat = ROOT.TLatex()
	lat.SetNDC(ROOT.kTRUE)
	lat.SetTextColor(ROOT.kBlack)
	lat.SetTextSize(size)
	lat.SetTextFont(font)
	lat.DrawLatex(x, y, string)

def drawTopLine(lumi, pbFb):
	latex = ROOT.TLatex()
	drawLatex('CMS Preliminary', 0.13, 0.92, 0.05, 62)    
	if pbFb   =='pb': drawLatex('L_{int.} = %4.0f pb^{-1}' %(lumi)      , 0.70, 0.92)
	elif pbFb =='fb': drawLatex('L_{int.} = %2.1f fb^{-1}' %(lumi/1000.), 0.70, 0.92)

def makeLegend(a,b,c,d,e=''):
	legend = ROOT.TLegend(a,b,c,d,e)
	legend.SetFillColor(ROOT.kWhite)
	legend.SetTextFont(42)
	legend.SetBorderSize(1)
	legend.SetMargin(0.35)
	legend.SetTextSize(0.03)
	#legend.SetTextFont(42)
	#stuff.append(legend)
	return legend

bold  = "\033[1m"
reset = "\033[0;0m"
colors = {'QCD'       : 38,
          'TTJets'    : ROOT.kRed-3,
          'WJets'     : ROOT.kYellow,
          'singleTop' : ROOT.kOrange+4,
          'DYJets'    : ROOT.kGreen+3,
          'data'      : ROOT.kBlack
          }

def makeCanvas():
  canvas = ROOT.TCanvas('eff_canvas','test')
  ROOT.gROOT.SetStyle('Plain')
  ROOT.gPad.UseCurrentStyle()
  ROOT.gROOT.ForceStyle()
  canvas.SetFillColor(ROOT.kWhite)
  canvas.SetGrid()
  canvas.GetFrame().SetFillColor(ROOT.kWhite)
  canvas.GetFrame().SetBorderSize(0)
  canvas.SetRightMargin(0.15)
  canvas.SetHighLightColor(ROOT.kWhite) 
  return canvas

def printCanvas( c, plotpath ):
  c.Print( plotpath + c.GetName() + '.root' )
#  c.Print( plotpath + c.GetName() + '.C'    )
  c.Print( plotpath + c.GetName() + '.png'  )
#  c.Print( plotpath + c.GetName() + '.eps'  )
#  c.Print( plotpath + c.GetName() + '.pdf'  )


def set1DStyle():
	#ROOT.gStyle = ROOT.gStyle
	style1d = ROOT.TStyle("mystyle","Style for P-TDR")
	# For the canvas:
	style1d.SetCanvasBorderMode(0)
	style1d.SetCanvasColor(ROOT.kWhite)
	style1d.SetCanvasDefH(600) #Height of canvas
	style1d.SetCanvasDefW(600) #Width of canvas
	style1d.SetCanvasDefX(0)   #POsition on screen
	style1d.SetCanvasDefY(0)
	# For the Pad:
	style1d.SetPadBorderMode(0)
	# style1d.SetPadBorderSize(Width_t size = 1)
	style1d.SetPadColor(ROOT.kWhite)
	style1d.SetPadGridX(False)
	style1d.SetPadGridY(False)
	style1d.SetGridColor(0)
	style1d.SetGridStyle(3)
	style1d.SetGridWidth(1)
	# For the frame:
	style1d.SetFrameBorderMode(0)
	style1d.SetFrameBorderSize(1)
	style1d.SetFrameFillColor(0)
	style1d.SetFrameFillStyle(0)
	style1d.SetFrameLineColor(1)
	style1d.SetFrameLineStyle(1)
	style1d.SetFrameLineWidth(1)
	# For the histo:
	# style1d.SetHistFillColor(1)
	# style1d.SetHistFillStyle(0)
	style1d.SetHistLineColor(1)
	style1d.SetHistLineStyle(0)
	style1d.SetHistLineWidth(1)
	# style1d.SetLegoInnerR(Float_t rad = 0.5)
	# style1d.SetNumberContours(Int_t number = 20)
	style1d.SetEndErrorSize(2)
	# style1d.SetErrorMarker(20)
	style1d.SetErrorX(0.)
	style1d.SetMarkerStyle(20)
	# For the fit/function:
	style1d.SetOptFit(1)
	style1d.SetFitFormat("5.4g")
	style1d.SetFuncColor(2)
	style1d.SetFuncStyle(1)
	style1d.SetFuncWidth(1)
	# For the date:
	style1d.SetOptDate(0)
	# style1d.SetDateX(Float_t x = 0.01)
	# style1d.SetDateY(Float_t y = 0.01)
	# For the statistics box:
	style1d.SetOptFile(0)
	style1d.SetOptStat(0) # To display the mean and RMS:   SetOptStat("mr")
	style1d.SetStatColor(ROOT.kWhite)
	style1d.SetStatFont(42)
	style1d.SetStatFontSize(0.025)
	style1d.SetStatTextColor(1)
	style1d.SetStatFormat("6.4g")
	style1d.SetStatBorderSize(1)
	style1d.SetStatH(0.1)
	style1d.SetStatW(0.15)
	# style1d.SetStatStyle(Style_t style = 1001)
	# style1d.SetStatX(Float_t x = 0)
	# style1d.SetStatY(Float_t y = 0)
	# Margins:
	style1d.SetPadTopMargin(0.1)
	style1d.SetPadBottomMargin(0.1)
	#style1d.SetPadLeftMargin(0.15)
	style1d.SetPadRightMargin(0.15)
	# For the Global title:
	style1d.SetOptTitle(0)
	style1d.SetTitleFont(42)
	style1d.SetTitleColor(1)
	style1d.SetTitleTextColor(1)
	style1d.SetTitleFillColor(10)
	style1d.SetTitleFontSize(0.05)
	# style1d.SetTitleH(0) # Set the height of the title box
	# style1d.SetTitleW(0) # Set the width of the title box
	# style1d.SetTitleX(0) # Set the position of the title box
	# style1d.SetTitleY(0.985) # Set the position of the title box
	# style1d.SetTitleStyle(Style_t style = 1001)
	# style1d.SetTitleBorderSize(2)
	# For the axis titles:
	style1d.SetTitleColor(1, "XYZ")
	style1d.SetTitleFont(42, "XYZ")
	style1d.SetTitleSize(0.06, "XYZ")
	# style1d.SetTitleXSize(Float_t size = 0.02) # Another way to set the size?
	# style1d.SetTitleYSize(Float_t size = 0.02)
	style1d.SetTitleXOffset(0.9)
	style1d.SetTitleYOffset(1.0)
	# style1d.SetTitleOffset(1.1, "Y") # Another way to set the Offset
	# For the axis labels:
	style1d.SetLabelColor(1, "XYZ")
	style1d.SetLabelFont(42, "XYZ")
	style1d.SetLabelOffset(0.007, "XYZ")
	style1d.SetLabelSize(0.05, "XYZ")
	# For the axis:
	style1d.SetAxisColor(1, "XYZ")
	style1d.SetStripDecimals(ROOT.kTRUE)
	style1d.SetTickLength(0.03, "XYZ")
	style1d.SetNdivisions(510, "XYZ")
	style1d.SetPadTickX(1)  # To get tick marks on the opposite side of the frame
	style1d.SetPadTickY(1)
	# Change for log plots:
	style1d.SetOptLogx(0)
	style1d.SetOptLogy(0)
	style1d.SetOptLogz(0)
	# Postscript options:
	style1d.SetPaperSize(20.,20.)
	#style1d.SetLineScalePS(Float_t scale = 3)
	#style1d.SetLineStyleString(Int_t i, const char* text)
	#style1d.SetHeaderPS(const char* header)
	#style1d.SetTitlePS(const char* pstitle)
	#style1d.SetBarOffset(Float_t baroff = 0.5)
	#style1d.SetBarWidth(Float_t barwidth = 0.5)
	#style1d.SetPaintTextFormat(const char* format = "g")
	#style1d.SetPalette(Int_t ncolors = 0, Int_t* colors = 0)
	#style1d.SetTimeOffset(Double_t toffset)
	#style1d.SetHistMinimumZero(ROOT.kTRUE)
	style1d.cd()
	return style1d

def set2DStyle():
  #ROOT.gStyle = ROOT.gStyle
  style2d = ROOT.TStyle("ROOT.gStyle","Style for P-TDR")
# For the canvas:
  style2d.SetCanvasBorderMode(0)
  style2d.SetCanvasColor(ROOT.kWhite)
  style2d.SetCanvasDefH(600) #Height of canvas
  style2d.SetCanvasDefW(600) #Width of canvas
  style2d.SetCanvasDefX(0)   #POsition on screen
  style2d.SetCanvasDefY(0)
# For the Pad:
  style2d.SetPadBorderMode(0)
  # style2d.SetPadBorderSize(Width_t size = 1)
  style2d.SetPadColor(ROOT.kWhite)
  style2d.SetPadGridX(False)
  style2d.SetPadGridY(False)
  style2d.SetGridColor(0)
  style2d.SetGridStyle(3)
  style2d.SetGridWidth(1)
# For the frame:
  style2d.SetFrameBorderMode(0)
  style2d.SetFrameBorderSize(1)
  style2d.SetFrameFillColor(0)
  style2d.SetFrameFillStyle(0)
  style2d.SetFrameLineColor(1)
  style2d.SetFrameLineStyle(1)
  style2d.SetFrameLineWidth(1)
# For the histo:
  # style2d.SetHistFillColor(1)
  # style2d.SetHistFillStyle(0)
  style2d.SetHistLineColor(1)
  style2d.SetHistLineStyle(0)
  style2d.SetHistLineWidth(1)
  # style2d.SetLegoInnerR(Float_t rad = 0.5)
  # style2d.SetNumberContours(Int_t number = 20)
  style2d.SetEndErrorSize(2)
#  style2d.SetErrorMarker(20)
  style2d.SetErrorX(0.)
  style2d.SetMarkerStyle(20)
  #For the fit/function:
  style2d.SetOptFit(1)
  style2d.SetFitFormat("5.4g")
  style2d.SetFuncColor(2)
  style2d.SetFuncStyle(1)
  style2d.SetFuncWidth(1)
#For the date:
  style2d.SetOptDate(0)
  # style2d.SetDateX(Float_t x = 0.01)
  # style2d.SetDateY(Float_t y = 0.01)
# For the statistics box:
  style2d.SetOptFile(0)
  style2d.SetOptStat(0) # To display the mean and RMS:   SetOptStat("mr")
  style2d.SetStatColor(ROOT.kWhite)
  style2d.SetStatFont(42)
  style2d.SetStatFontSize(0.025)
  style2d.SetStatTextColor(1)
  style2d.SetStatFormat("6.4g")
  style2d.SetStatBorderSize(1)
  style2d.SetStatH(0.1)
  style2d.SetStatW(0.15)
  # style2d.SetStatStyle(Style_t style = 1001)
  # style2d.SetStatX(Float_t x = 0)
  # style2d.SetStatY(Float_t y = 0)
# Margins:
  style2d.SetPadTopMargin(0.1)
  style2d.SetPadBottomMargin(0.15)
  style2d.SetPadLeftMargin(0.15)
  style2d.SetPadRightMargin(0.20)
# For the Global title:
  style2d.SetOptTitle(0)
  style2d.SetTitleFont(42)
  style2d.SetTitleColor(1)
  style2d.SetTitleTextColor(1)
  style2d.SetTitleFillColor(10)
  style2d.SetTitleFontSize(0.05)
  # style2d.SetTitleH(0) # Set the height of the title box
  # style2d.SetTitleW(0) # Set the width of the title box
  # style2d.SetTitleX(0) # Set the position of the title box
  # style2d.SetTitleY(0.985) # Set the position of the title box
  # style2d.SetTitleStyle(Style_t style = 1001)
  # style2d.SetTitleBorderSize(2)
  # For the axis titles:
  style2d.SetTitleColor(1, "XYZ")
  style2d.SetTitleFont(42, "XYZ")
  style2d.SetTitleSize(0.06, "XYZ")
  # style2d.SetTitleXSize(Float_t size = 0.02) # Another way to set the size?
  # style2d.SetTitleYSize(Float_t size = 0.02)
  style2d.SetTitleXOffset(0.9)
  style2d.SetTitleYOffset(1.0)
  # style2d.SetTitleOffset(1.1, "Y") # Another way to set the Offset
# For the axis labels:
  style2d.SetLabelColor(1, "XYZ")
  style2d.SetLabelFont(42, "XYZ")
  style2d.SetLabelOffset(0.007, "XYZ")
  style2d.SetLabelSize(0.05, "XYZ")
# For the axis:
  style2d.SetAxisColor(1, "XYZ")
  style2d.SetStripDecimals(ROOT.kTRUE)
  style2d.SetTickLength(0.03, "XYZ")
  style2d.SetNdivisions(510, "XYZ")
  style2d.SetPadTickX(1)  # To get tick marks on the opposite side of the frame
  style2d.SetPadTickY(1)
# Change for log plots:
  style2d.SetOptLogx(0)
  style2d.SetOptLogy(0)
  style2d.SetOptLogz(0)
# Postscript options:
  style2d.SetPaperSize(20.,20.)
  # style2d.SetLineScalePS(Float_t scale = 3)
  # style2d.SetLineStyleString(Int_t i, const char* text)
  # style2d.SetHeaderPS(const char* header)
  # style2d.SetTitlePS(const char* pstitle)
  # style2d.SetBarOffset(Float_t baroff = 0.5)
  # style2d.SetBarWidth(Float_t barwidth = 0.5)
  # style2d.SetPaintTextFormat(const char* format = "g")
  # style2d.SetPalette(Int_t ncolors = 0, Int_t* colors = 0)
  # style2d.SetTimeOffset(Double_t toffset)
  # style2d.SetHistMinimumZero(ROOT.kTRUE)
  style2d.cd()
  return style2d

