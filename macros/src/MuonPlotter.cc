/*****************************************************************************
*   Collection of tools for producing plots for Muon Fake Rate Analysis      *
*                                                                            *
*                                                  (c) 2010 Benjamin Stieger *
*****************************************************************************/
#include "MuonPlotter.hh"
#include "helper/AnaClass.hh"
#include "helper/Utilities.hh"
#include "helper/FPRatios.hh"

#include "TLorentzVector.h"
#include "TGraphAsymmErrors.h"

using namespace std;

int gSWITCH = 0; // 0: high/pt, low HT, hybiso, 1: low/pt, high HT, stand iso

// Binning ///////////////////////////////////////////////////////////////////////
static const int gNPtbins = 5;
static const double gPtbins[gNPtbins+1] = {20., 30., 40., 50., 65., 80.};
static const int gNPt2bins = 6;
static const double gPt2bins[gNPt2bins+1] = {10., 20., 30., 40., 50., 65., 80.};

// static const int gNPtbins = 7;
// static const double gPtbins[gNPtbins+1] = {5., 10., 20., 30., 40., 50., 65., 80.};
// static const int gNPt2bins = 7;
// static const double gPt2bins[gNPt2bins+1] = {5., 10., 20., 30., 40., 50., 65., 80.};

static const int gNEtabins = 1;
static const double gEtabins[gNEtabins+1] = {-2.4, 2.4};

//////////////////////////////////////////////////////////////////////////////////


//____________________________________________________________________________
MuonPlotter::MuonPlotter(){
// Default constructor, no samples are set
}

//____________________________________________________________________________
MuonPlotter::MuonPlotter(TString outputdir){
// Explicit constructor with output directory
	setOutputDir(outputdir);
}

//____________________________________________________________________________
MuonPlotter::MuonPlotter(TString outputdir, TString outputfile){
// Explicit constructor with output directory and output file
	setOutputDir(outputdir);
	setOutputFile(outputfile);
}

//____________________________________________________________________________
MuonPlotter::~MuonPlotter(){}

//____________________________________________________________________________
void MuonPlotter::init(TString filename){
	if(fVerbose > 0) cout << "------------------------------------" << endl;
	if(fVerbose > 0) cout << "Initializing MuonPlotter ... " << endl;
	Util::SetStyle();
	loadSamples(filename);
	readVarNames("anavarnames.dat");
	fOutputFileName = fOutputDir + "Yields.root";

	// fLumiNorm = fSamples[0].lumi; // Normalize everything to this lumi in /pb
	fLumiNorm = 100; // Normalize everything to this lumi in /pb
	fBinWidthScale = 10.; // Normalize Y axis to this binwidth

	fMinPt1 = 20.;
	fMinPt2 = 10.;

	// Prevent root from adding histograms to current file
	TH1::AddDirectory(kFALSE);

	fMCBG.push_back(TTbar);
	fMCBG.push_back(WJets);
	fMCBG.push_back(ZJets);
	fMCBG.push_back(VVJets);
	fMCBG.push_back(QCD15);
	fMCBG.push_back(QCD30);
	fMCBG.push_back(QCD80);
	fMCBG.push_back(QCD170);

	fMCBGSig = fMCBG;
	fMCBGSig.push_back(LM0);

	fMuData.push_back(MuA);
	fMuData.push_back(MuB);
	fEGData.push_back(EGA);
	fEGData.push_back(EGB);
	fJMData.push_back(JMA);
	fJMData.push_back(JMB);

	fAllSamples.push_back(MuA);
	fAllSamples.push_back(MuB);
	fAllSamples.push_back(EGA);
	fAllSamples.push_back(EGB);
	fAllSamples.push_back(JMA);
	fAllSamples.push_back(JMB);
	fAllSamples.push_back(TTbar);
	fAllSamples.push_back(WJets);
	fAllSamples.push_back(ZJets);
	fAllSamples.push_back(VVJets);
	fAllSamples.push_back(QCD15);
	fAllSamples.push_back(QCD30);
	fAllSamples.push_back(QCD80);
	fAllSamples.push_back(QCD170);
	fAllSamples.push_back(LM0);

	bookHistos();
}

//____________________________________________________________________________
void MuonPlotter::loadSamples(const char* filename){
	char buffer[200];
	ifstream IN(filename);

	char ParName[100], StringValue[1000];
	float ParValue;

	if(fVerbose > 0) cout << "------------------------------------" << endl;
	if(fVerbose > 0) cout << "Sample File  " << filename << endl;
	int counter(0);

	while( IN.getline(buffer, 200, '\n') ){
		// ok = false;
		if (buffer[0] == '#') continue; // Skip lines commented with '#'
		if( !strcmp(buffer, "SAMPLE")){
			sample s;

			IN.getline(buffer, 200, '\n');
			sscanf(buffer, "Name\t%s", StringValue);
			s.name = TString(StringValue);

			IN.getline(buffer, 200, '\n');
			sscanf(buffer, "SName\t%s", StringValue);
			s.sname = TString(StringValue);

			IN.getline(buffer, 200, '\n');
			sscanf(buffer, "File\t%s", StringValue);
			TFile *f = TFile::Open(StringValue);
			s.file = f;
			s.tree = (TTree*)f->Get("Analysis");
			if(s.tree == NULL){ cout << " Tree not found in file!" << endl; break; }

			IN.getline(buffer, 200, '\n');
			sscanf(buffer, "Lumi\t%f", &ParValue);
			s.lumi = ParValue;

			IN.getline(buffer, 200, '\n');
			sscanf(buffer, "IsData\t%f", &ParValue);
			s.isdata = (bool)ParValue;

			if(fVerbose > 0){
				cout << " ---- " << endl;
				cout << "  New sample added: " << s.name << endl;
				cout << "   Sample no.  " << counter << endl;
				cout << "   Short name: " << s.sname << endl;
				cout << "   File:       " << (s.file)->GetName() << endl;
				cout << "   Events:     " << s.tree->GetEntries() << endl;
				cout << "   Lumi:       " << s.lumi << endl;
				cout << "   IsData:     " << s.isdata << endl;
			}
			fSampleMap[s.sname] = counter;
			fSamples.push_back(s);
			counter++;
		}
	}
	if(fVerbose > 0) cout << "------------------------------------" << endl;
}

//____________________________________________________________________________
void MuonPlotter::makePlots(){
	if(readHistos(fOutputFileName) != 0) return;
	
	// float scale = fSamples[MuA].lumi + fSamples[MuB].lumi;
	fLumiNorm = 22.;
	// printYields(fLumiNorm);

	makeIntPrediction(Muon);
	makeIntPrediction(Electron);
	makeIntPrediction(EMu);

	// makefRatioPlots();
	// makepRatioPlots();

	// makeDiffPredictionPlots();

	// fLumiNorm = 22.;
	// makeIsoVsPtPlot(fMCBG, 0, &MuonPlotter::isSigSupMuEvent, &MuonPlotter::isLooseMuon, fMuData, 0, &MuonPlotter::isSigSupMuEventTRG, &MuonPlotter::isLooseMuon, "IsoVsPt_SigSuppressed", false);

}

//____________________________________________________________________________
void MuonPlotter::makeDiffPredictionPlots(){
	// Fill the ratios
	fLumiNorm = 1000.;
	// printYields(fLumiNorm);

	cout << "Producing prediction for :" << endl;
	for(size_t i = 0; i < fMCBGSig.size(); ++i){
		int ind = fMCBGSig[i];
		cout << " " << fSamples[ind].sname << flush;
	}
	cout << endl;

	// fillfRatio(fJMData, 0);
	// fillpRatio(fMuData, 0);
	// 
	// makeSSPredictionPlots(fMuData);

	fillfRatio(fMCBGSig, 0);
	fillpRatio(fMCBGSig, 0);
	
	makeSSPredictionPlots(fMCBGSig);
}

//____________________________________________________________________________
void MuonPlotter::makefRatioPlots(){
	// TH1D *h_fdata  = fillRatioPt(fMuData, 0, &MuonPlotter::isSigSupMuEventTRG, &MuonPlotter::isLooseMuon);      // JetMET Dataset (Single Muon Selection)
	// TH1D *h_fttbar = fillRatioPt(TTbar,   0, &MuonPlotter::isGoodMuEvent,                  &MuonPlotter::isFakeTTbarMuon); // TTbarJets MC
	// TH1D *h_fallmc = fillRatioPt(fMCBG,   0, &MuonPlotter::isSigSupMuEvent,    &MuonPlotter::isLooseMuon);      // QCD MC
	fLumiNorm = fSamples[JMA].lumi + fSamples[JMB].lumi;
	TH1D *h_fdata1 = fillRatioPt(fJMData, 1);      // JetMET Dataset (Single Muon Selection)
	TH1D *h_fdata2 = fillRatioPt(fMuData, 1);
	TH1D *h_fallmc = fillRatioPt(fMCBG,   1);      // QCD MC
	TH1D *h_fttbar = fillRatioPt(TTbar,   0, &MuonPlotter::isGoodMuEvent, &MuonPlotter::isFakeTTbarMuon); // TTbarJets MC
	h_fdata1->SetName("fRatioData");
	h_fdata2->SetName("fRatioDataMu");
	h_fttbar->SetName("fRatioTTbar");
	h_fallmc->SetName("fRatioAllMC");

	// setPlottingRange(h_fdata1, h_fttbar, h_fallmc);

	h_fdata1->SetMinimum(0.);
	h_fdata2->SetMinimum(0.);
	h_fttbar->SetMinimum(0.);
	h_fallmc->SetMinimum(0.);

	h_fdata1->SetMarkerColor(kBlack);
	h_fdata2->SetMarkerColor(kBlue);
	h_fttbar->SetMarkerColor(kBlue);
	h_fallmc->SetMarkerColor(kRed);

	h_fdata1->SetMarkerStyle(20);
	h_fdata2->SetMarkerStyle(20);
	h_fttbar->SetMarkerStyle(20);
	h_fallmc->SetMarkerStyle(20);

	// h_fdata1->SetMarkerSize(2);
	// h_fdata2->SetMarkerSize(2);
	// h_fttbar->SetMarkerSize(2);
	// h_fallmc->SetMarkerSize(2);

	h_fdata1->SetLineWidth(2);
	h_fdata2->SetLineWidth(2);
	h_fttbar->SetLineWidth(2);
	h_fallmc->SetLineWidth(2);

	h_fdata1->SetLineColor(kBlack);
	h_fdata2->SetLineColor(kBlue);
	h_fttbar->SetLineColor(kBlue);
	h_fallmc->SetLineColor(kRed);

	h_fdata1->SetFillColor(kBlack);
	h_fdata2->SetFillColor(kBlue);
	h_fttbar->SetFillColor(kBlue);
	h_fallmc->SetFillColor(kRed);

	plotRatioOverlay3H(h_fdata1, "Data (Jet, L = 21.7 pb^{-1})", h_fttbar, "t#bar{t} Fake GenMatch", h_fallmc, "QCD, t#bar{t}+jets, V+jets");
	// setPlottingRange(h_fdata1, h_fdata2);
	plotRatioOverlay2H(h_fdata1, "Data (Jet, L = 21.7 pb^{-1})", h_fdata2, "Data (Muon, L = 21.3 pb^{-1})");
	// plotRatioOverlay3H(h_fdata2, "Data (Mu, L = 21.3 pb^{-1})", h_fttbar, "t#bar{t} Fake GenMatch", h_fallmc, "QCD, t#bar{t}+jets, V+jets");
}

//____________________________________________________________________________
void MuonPlotter::makepRatioPlots(){
	// TH1D *h_pdata  = fillRatioPt(fMuData, 0, &MuonPlotter::isZMuMuEventTRG, &MuonPlotter::isLooseMuon); // Mu Dataset (Di Muon Selection)
	// TH1D *h_pttbar = fillRatioPt(TTbar,   0, &MuonPlotter::isGoodMuEvent, &MuonPlotter::isPromptTTbarMuon); // TTbar
	// TH1D *h_pallmc = fillRatioPt(fMCBG,   0, &MuonPlotter::isZMuMuEvent,    &MuonPlotter::isLooseMuon); // all MC
	fLumiNorm = fSamples[MuA].lumi + fSamples[MuB].lumi;
	TH1D *h_pdata  = fillRatioPt(fMuData, 2);
	TH1D *h_pallmc = fillRatioPt(fMCBG,   2);
	TH1D *h_pttbar = fillRatioPt(TTbar,   0, &MuonPlotter::isGoodMuEvent, &MuonPlotter::isPromptTTbarMuon); // TTbar
	h_pdata ->SetName("pRatioData");
	h_pttbar->SetName("pRatioTTbar");
	h_pallmc->SetName("pRatioAllMC");

	setPlottingRange(h_pdata, h_pttbar, h_pallmc);

	h_pdata ->Draw("goff");
	h_pttbar->Draw("goff");
	h_pallmc->Draw("goff");

	h_pdata ->SetMinimum(0.);
	h_pttbar->SetMinimum(0.);
	h_pallmc->SetMinimum(0.);
	h_pttbar->SetMaximum(1.3);
	h_pallmc->SetMaximum(1.3);

	h_pdata ->SetLineWidth(2);
	h_pttbar->SetLineWidth(2);
	h_pallmc->SetLineWidth(2);

	h_pdata ->SetMarkerColor(kBlack);
	h_pttbar->SetMarkerColor(kBlue);
	h_pallmc->SetMarkerColor(kRed);

	h_pdata ->SetMarkerStyle(20);
	h_pttbar->SetMarkerStyle(20);
	h_pallmc->SetMarkerStyle(20);

	// h_pdata ->SetMarkerSize(2);
	// h_pttbar->SetMarkerSize(2);
	// h_pallmc->SetMarkerSize(2);

	h_pdata ->SetLineColor(kBlack);
	h_pttbar->SetLineColor(kBlue);
	h_pallmc->SetLineColor(kRed);

	h_pdata ->SetFillColor(kBlack);
	h_pttbar->SetFillColor(kBlue);
	h_pallmc->SetFillColor(kRed);

	h_pdata ->SetDrawOption("E1");
	h_pttbar->SetDrawOption("E1");
	h_pallmc->SetDrawOption("E1");

	// plotRatioOverlay3H(h_pdata, "Data (Jet, L = 21.7 pb^{-1})", h_pttbar, "t#bar{t} Prompt GenMatch", h_pallmc, "QCD, t#bar{t}+jets, V+jets");
	plotRatioOverlay3H(h_pdata, "Data (Mu, L = 21.3 pb^{-1})", h_pttbar, "t#bar{t} Prompt GenMatch", h_pallmc, "QCD, t#bar{t}+jets, V+jets");
}

//____________________________________________________________________________
void MuonPlotter::makeIsolationPlots(){
	const int nbins = 20;
	TH1D *h_prompt      = new TH1D("h_prompt",  "Isolation for prompt Muons in WJets", nbins, 0., 1.);
	TH1D *h_fakew       = new TH1D("h_fakew",   "Isolation for fake Muons in WJets events", nbins, 0., 1.);
	TH1D *h_fakeqcd     = new TH1D("h_fakeqcd", "Isolation for fake Muons in QCD events", nbins, 0., 1.);
	TH1D *h_wtau        = new TH1D("h_wtau",   "Isolation for fake Muons from tau in WJets events", nbins, 0., 1.);
	TH1D *h_wnotau      = new TH1D("h_wnotau", "Isolation for fake Muons not from tau in WJets events", nbins, 0., 1.);

	TH1D *h_ttp      = new TH1D("h_ttp", "Isolation for prompt Muons in ttbar", nbins, 0., 1.);
	TH1D *h_ttf      = new TH1D("h_ttf", "Isolation for non prompt Muons in ttbar", nbins, 0., 1.);
	TH1D *h_ttftau   = new TH1D("h_ttftau", "Isolation for Muons from tau in ttbar", nbins, 0., 1.);
	TH1D *h_ttfnotau = new TH1D("h_ttfnotau", "Isolation for Muons not from tau in ttbar", nbins, 0., 1.);

	TH1D *h_qcdb   = new TH1D("h_qcdb",   "Isolation for muons from bottom hadrons in QCD", nbins, 0., 1.);
	TH1D *h_qcdpik = new TH1D("h_qcdpik", "Isolation for muons from pions/kaons in QCD", nbins, 0., 1.);
	TH1D *h_ttb    = new TH1D("h_ttb",    "Isolation for muons from bottom hadrons in ttbar", nbins, 0., 1.);

	TH1D *h_zjets  = new TH1D("h_zjets",  "Isolation for muons from Z boson decays", nbins, 0., 1.);

	h_fakeqcd->SetXTitle(convertVarName("MuIso[0]"));
	h_qcdb->SetXTitle(convertVarName("MuIso[0]"));
	h_prompt->SetXTitle(convertVarName("MuIso[0]"));
	h_ttp->SetXTitle(convertVarName("MuIso[0]"));
	h_zjets->SetXTitle(convertVarName("MuIso[0]"));

	fSamples[1].tree->Project("h_prompt",  "MuIso[0]", "abs(MuGenMoID[0])==24");
	fSamples[1].tree->Project("h_fakew",   "MuIso[1]", "abs(MuGenMoID[1])!=24&&abs(MuGenGMoID[1])!=24");
	fSamples[0].tree->Project("h_fakeqcd", "MuIso[0]", "");
	fSamples[1].tree->Project("h_wtau",    "MuIso[1]", "abs(MuGenMoID[1])==15");
	fSamples[1].tree->Project("h_wnotau",  "MuIso[1]", "abs(MuGenMoID[1])!=24&&abs(MuGenMoID[1])!=15");
	fSamples[2].tree->Project("h_ttp",     "MuIso[0]", "abs(MuGenMoID[0])==24");
	fSamples[2].tree->Project("h_ttf",     "MuIso[0]", "abs(MuGenMoID[0])!=24");
	fSamples[2].tree->Project("h_ttftau",  "MuIso[0]", "abs(MuGenMoID[0])==15");
	fSamples[2].tree->Project("h_ttfnotau","MuIso[0]", "abs(MuGenMoID[0])!=24&&abs(MuGenMoID[0])!=15");
	fSamples[2].tree->Project("h_ttb",     "MuIso[1]", "MuGenMoType[1]==15||MuGenMoType[1]==17||MuGenMoType[1]==21");

	fSamples[4].tree->Project("h_qcdb",    "MuIso[0]", "MuGenMoType[0]==15||MuGenMoType[0]==17||MuGenMoType[0]==21");
	fSamples[4].tree->Project("h_qcdpik",  "MuIso[0]", "MuGenMoType[0]==11||MuGenMoType[0]==12||MuGenMoType[0]==13");

	fSamples[5].tree->Project("h_zjets",   "MuIso[0]", "abs(MuGenMoID[0]==23)");

	h_fakeqcd->SetMinimum(0);
	printHisto(h_fakeqcd, "QCDIso", "Isolation of QCD Muons", "hist");

	plotOverlay3H(h_fakeqcd, "Fake in QCD", h_prompt, "Prompt", h_fakew, "Fake in WJets", false, 0.15);
	plotOverlay3H(h_fakeqcd, "QCD", h_wtau, "W: tau", h_wnotau, "W: No tau", false, 0.15);

	plotOverlay3H(h_prompt, "Prompt W", h_ttp, "Prompt ttbar", h_ttftau, "ttbar tau", false, 0.15);
	plotOverlay3H(h_fakeqcd, "QCD", h_wnotau, "W: no tau", h_ttfnotau, "ttbar: no tau", false, 0.15);

	plotOverlay3H(h_qcdb, "QCD: b", h_qcdpik, "QCD: #pi/K", h_ttb, "ttbar: b", false, 0.15);

	plotOverlay2H(h_fakeqcd, "QCD", h_ttfnotau, "ttbar", false, 0.15);
	plotOverlay2H(h_zjets, "ZJets", h_ttp,      "ttbar", true, 0.15);
}

//____________________________________________________________________________
void MuonPlotter::makePtPlots(){
	const int nbins = 40;
	TH1D *h_prompt   = new TH1D("h_prompt",  "Pt for prompt Muons in WJets",                   nbins, 0., 100.);
	TH1D *h_fakew    = new TH1D("h_fakew",   "Pt for fake Muons in WJets events",              nbins, 0., 100.);
	TH1D *h_fakeqcd  = new TH1D("h_fakeqcd", "Pt for fake Muons in QCD events",                nbins, 0., 100.);
	TH1D *h_wtau     = new TH1D("h_wtau",    "Pt for fake Muons from tau in WJets events",     nbins, 0., 100.);
	TH1D *h_wnotau   = new TH1D("h_wnotau",  "Pt for fake Muons not from tau in WJets events", nbins, 0., 100.);
	TH1D *h_ttp      = new TH1D("h_ttp",     "Pt for prompt Muons in ttbar",                   nbins, 0., 100.);
	TH1D *h_ttp2     = new TH1D("h_ttp2",    "Pt for prompt Muons in ttbar",                   nbins, 0., 200.);
	TH1D *h_ttf      = new TH1D("h_ttf",     "Pt for non prompt Muons in ttbar",               nbins, 0., 100.);
	TH1D *h_ttftau   = new TH1D("h_ttftau",  "Pt for Muons from tau in ttbar",                 nbins, 0., 100.);
	TH1D *h_ttfnotau = new TH1D("h_ttfnotau","Pt for Muons not from tau in ttbar",             nbins, 0., 100.);
	TH1D *h_qcdb     = new TH1D("h_qcdb",   "Pt for muons from bottom hadrons in QCD",         nbins, 0., 100.);
	TH1D *h_qcdpik   = new TH1D("h_qcdpik", "Pt for muons from pions/kaons in QCD",            nbins, 0., 100.);
	TH1D *h_ttb      = new TH1D("h_ttb",    "Pt for muons from bottom hadrons in ttbar",       nbins, 0., 100.);
	TH1D *h_z        = new TH1D("h_z",      "Pt for muons from Z boson decays",                nbins, 0., 100.);

	h_fakeqcd->SetXTitle(convertVarName("MuPt[0]"));
	h_qcdb->SetXTitle(convertVarName("MuPt[0]"));
	h_prompt->SetXTitle(convertVarName("MuPt[0]"));
	h_ttp2->SetXTitle(convertVarName("MuPt[0]"));

	fSamples[6].tree->Project("h_prompt",  "MuPt[0]", "abs(MuGenMoID[0])==24");
	fSamples[5].tree->Project("h_z",       "MuPt[0]", "abs(MuGenMoID[0])==23");
	fSamples[1].tree->Project("h_fakew",   "MuPt[1]", "abs(MuGenMoID[1])!=24&&abs(MuGenGMoID[1])!=24");
	fSamples[0].tree->Project("h_fakeqcd", "MuPt[0]", "");
	fSamples[1].tree->Project("h_wtau",    "MuPt[1]", "abs(MuGenMoID[1])==15");
	fSamples[1].tree->Project("h_wnotau",  "MuPt[1]", "abs(MuGenMoID[1])!=24&&abs(MuGenMoID[1])!=15");
	fSamples[2].tree->Project("h_ttp",     "MuPt[0]", "abs(MuGenMoID[0])==24");
	fSamples[2].tree->Project("h_ttp2",    "MuPt[0]", "abs(MuGenMoID[0])==24");
	fSamples[2].tree->Project("h_ttf",     "MuPt[0]", "abs(MuGenMoID[0])!=24");
	fSamples[2].tree->Project("h_ttftau",  "MuPt[0]", "abs(MuGenMoID[0])==15");
	fSamples[2].tree->Project("h_ttfnotau","MuPt[0]", "abs(MuGenMoID[0])!=24&&abs(MuGenMoID[0])!=15");
	fSamples[2].tree->Project("h_ttb",     "MuPt[1]", "MuGenMoType[1]==15||MuGenMoType[1]==17||MuGenMoType[1]==21");
	fSamples[4].tree->Project("h_qcdb",    "MuPt[0]", "MuGenMoType[0]==15||MuGenMoType[0]==17||MuGenMoType[0]==21");
	fSamples[4].tree->Project("h_qcdpik",  "MuPt[0]", "MuGenMoType[0]==11||MuGenMoType[0]==12||MuGenMoType[0]==13");

	cout << h_fakeqcd->GetEntries() << " " << h_prompt->GetEntries() << " " << h_fakew->GetEntries() << endl;

	// h_ttp2->SetMinimum(0);
	printHisto(h_ttp2, "TTbarPt", "Pt of prompt ttbar Muons", "hist");

	plotOverlay3H(h_fakeqcd, "QCD",         h_prompt, "W-jets",       h_ttp,      "ttbar",         true);
	plotOverlay3H(h_z,       "Z-jets",      h_prompt, "W-jets",       h_ttp,      "ttbar",         true);
	plotOverlay3H(h_fakeqcd, "QCD",         h_wnotau, "W-jets",       h_ttfnotau, "ttbar",         true);

	// plotOverlay3H(h_fakeqcd, "Fake in QCD", h_prompt, "Prompt",       h_fakew,    "Fake in WJets", true);
	// plotOverlay3H(h_fakeqcd, "QCD",         h_wtau,   "W: tau",       h_wnotau,   "W: No tau",     true);
	// plotOverlay3H(h_prompt,  "Prompt W",    h_ttp,    "Prompt ttbar", h_ttftau,   "ttbar tau",     true);
	// plotOverlay3H(h_fakeqcd, "QCD",         h_wnotau, "W: no tau",    h_ttfnotau, "ttbar: no tau", true);
	// plotOverlay3H(h_qcdb,    "QCD: b",      h_qcdpik, "QCD: #pi/K",   h_ttb,      "ttbar: b",      true);

	// plotOverlay2H(h_fakeqcd, "QCD", h_fakew, "WJets", false, 0.15);
}

//____________________________________________________________________________
void MuonPlotter::makeIsoVsPtPlot(int sample1, int muon1, TCut c1, int sample2, int muon2, TCut c2, TString outputname, bool logy){
	const int nbins = 20;
	TH2D *h2_sig = new TH2D("h2_sig", "Isolation vs Pt for fake muons in signal", nbins, 0., 1., gNPtbins, gPtbins);
	TH2D *h2_bg  = new TH2D("h2_bg", "Isolation vs Pt for muons in background",   nbins, 0., 1., gNPtbins, gPtbins);
	h2_sig->SetXTitle(convertVarName("MuIso[0]"));
	h2_sig->SetYTitle(convertVarName("MuPt[0]"));
	h2_bg->SetXTitle(convertVarName("MuIso[0]"));
	h2_bg->SetYTitle(convertVarName("MuPt[0]"));

	fSamples[sample1].tree->Project("h2_bg", Form("MuPt[%d]:MuIso[%d]", muon1, muon1), c1);
	fSamples[sample2].tree->Project("h2_sig", Form("MuPt[%d]:MuIso[%d]", muon2, muon2), c2);

	TLatex *lat = new TLatex();
	lat->SetNDC(kTRUE);
	lat->SetTextColor(kBlack);
	lat->SetTextSize(0.06);

	TCanvas *c_temp = new TCanvas("IsoVsPt", "Isolating in Pt bins", 0, 0, 1200, 800);
	c_temp->Divide(3,2);
	c_temp->cd(6);
	if(logy) gPad->SetLogz(1);
	h2_bg->DrawCopy("colz");
	lat->DrawLatex(0.11,0.92, fSamples[sample1].sname);

	for(size_t i = 1; i <= 5; ++i){
		c_temp->cd(i);
		gStyle->SetOptStat(1111);
		TH1D *h1 = h2_bg->ProjectionX("h2_bgx",i, i);	
		TH1D *h2 = h2_sig->ProjectionX("h2_sigx",i, i);
		h1->SetXTitle(convertVarName("MuIso[0]"));
		h1->SetLineWidth(2);
		h1->SetFillColor(15);
		h1->SetFillStyle(1001);
		h2->SetLineWidth(2);
		h2->SetLineColor(kBlue);
		h2->SetFillColor(kBlue);
		h2->SetFillStyle(3004);

		if(logy) gPad->SetLogy(1);
		gPad->SetFillStyle(0);
		h1->Sumw2();
		h2->Sumw2();

		// Scaling
		if(h1->GetEntries() > 0 ) h1->Scale(1.0/h1->Integral());
		if(h2->GetEntries() > 0 ) h2->Scale(1.0/h2->Integral());

		// Determine plotting range
		double max1 = h1->GetMaximum();
		double max2 = h2->GetMaximum();
		double max  = (max1>max2)?max1:max2;
		if(logy) max = 5*max;
		else max = 1.05*max;
		h1->SetMaximum(max);
		h2->SetMaximum(max);
		if(!logy){
			h1->SetMinimum(0.0);
			h2->SetMinimum(0.0);
		}

		TLegend *leg = new TLegend(0.55,0.75,0.75,0.88);
		if(i == 1){
			leg->AddEntry(h1, fSamples[sample1].sname,"f");
			leg->AddEntry(h2, fSamples[sample2].sname,"f");
			leg->SetFillStyle(0);
			leg->SetTextFont(42);
			leg->SetBorderSize(0);
		}

		TH1D *h1_temp = (TH1D*)h1->DrawCopy("hist");
		TH1D *h2_temp = (TH1D*)h2->DrawCopy("histsames");
		gPad->Update();
		TPaveStats *s1 = (TPaveStats*)h1_temp->GetListOfFunctions()->FindObject("stats");
		TPaveStats *s2 = (TPaveStats*)h2_temp->GetListOfFunctions()->FindObject("stats");
		s2->SetTextColor(kBlue); s2->SetLineColor(kBlue);
		s2->SetY1NDC(s1->GetY1NDC() - (s1->GetY2NDC() - s1->GetY1NDC()));
		s2->SetY2NDC(s1->GetY1NDC());

		if(i==1) leg->Draw();

		double min1 = h1->GetYaxis()->GetXmin();
		double min2 = h2->GetYaxis()->GetXmin();
		double min  = (min1<min2)?min1:min2;

		TLine *l1 = new TLine(0.15, min, 0.15, max);
		l1->SetLineColor(kRed);
		l1->SetLineWidth(2);
		l1->Draw();

		lat->SetTextColor(kBlack);
		lat->SetTextSize(0.05);
		lat->DrawLatex(0.11,0.92, Form("p_{T}(#mu) %3.0f - %3.0f GeV", gPtbins[i-1], gPtbins[i]));

		int bin0 = h1->FindBin(0.00);
		int bin15 = h1->FindBin(0.15);
		int bin1 = h1->FindBin(1.00);
		float f1 = h1->Integral(bin0, bin15) / h1->Integral(bin0, bin1);
		float f2 = h2->Integral(bin0, bin15) / h2->Integral(bin0, bin1);
		lat->SetTextSize(0.04);
		lat->DrawLatex(0.55,0.905, Form("ratio = %4.2f", f1));
		lat->SetTextColor(kBlue);
		lat->DrawLatex(0.55,0.945, Form("ratio = %4.2f", f2));

		gPad->RedrawAxis();
	}
	Util::PrintNoEPS(c_temp, outputname, fOutputDir, fOutputFile);
}

//____________________________________________________________________________
void MuonPlotter::makeIsoVsPtPlot(int sample1, int muon1, bool(MuonPlotter::*eventSelector1)(), bool(MuonPlotter::*muonSelector1)(int), int sample2, int muon2, bool(MuonPlotter::*eventSelector2)(), bool(MuonPlotter::*muonSelector2)(int), TString outputname, bool logy){
	vector<int> samples1; samples1.push_back(sample1);
	vector<int> samples2; samples2.push_back(sample2);
	makeIsoVsPtPlot(samples1, muon1, eventSelector1, muonSelector1, samples2, muon2, eventSelector2, muonSelector2, outputname, logy);
}
void MuonPlotter::makeIsoVsPtPlot(vector<int> samples1, int muon1, bool(MuonPlotter::*eventSelector1)(), bool(MuonPlotter::*muonSelector1)(int), vector<int> samples2, int muon2, bool(MuonPlotter::*eventSelector2)(), bool(MuonPlotter::*muonSelector2)(int), TString outputname, bool logy){
	const int nbins = 30;
	TH2D *h2_sig = new TH2D("h2_sig", "Isolation vs Pt for muons in data",       nbins, 0., 1., gNPt2bins, gPt2bins);
	TH2D *h2_bg  = new TH2D("h2_bg",  "Isolation vs Pt for fake muons in ttbar", nbins, 0., 1., gNPt2bins, gPt2bins);
	h2_sig->SetXTitle(convertVarName("MuIso[0]"));
	h2_bg ->SetXTitle(convertVarName("MuIso[0]"));
	h2_sig->SetYTitle(convertVarName("MuPt[0]"));
	h2_bg ->SetYTitle(convertVarName("MuPt[0]"));
	
	TTree *tree = NULL;
	for(size_t i = 0; i < samples1.size(); ++i){
		int index = samples1[i];
	
		tree = fSamples[index].tree;
		float scale = fLumiNorm / fSamples[index].lumi;
		if(fSamples[index].isdata) scale = 1.;
		tree->ResetBranchAddresses();
		Init(tree);
		if (fChain == 0) return;
		Long64_t nentries = fChain->GetEntriesFast();
		Long64_t nbytes = 0, nb = 0;
		for (Long64_t jentry=0; jentry<nentries;jentry++) {
			Long64_t ientry = LoadTree(jentry);
			if (ientry < 0) break;
			nb = fChain->GetEntry(jentry);   nbytes += nb;

			if((*this.*eventSelector1)() == false) continue;
			if((*this.*muonSelector1)(muon1) == false) continue;
			h2_sig->Fill(MuIso[muon1], MuPt[muon1], scale);
		}
	}

	for(size_t i = 0; i < samples1.size(); ++i){
		int index = samples2[i];

		tree = fSamples[index].tree;
		float scale = fLumiNorm / fSamples[index].lumi;
		if(fSamples[index].isdata) scale = 1.;
		tree->ResetBranchAddresses();
		Init(tree);
		if (fChain == 0) return;
		Long64_t nentries = fChain->GetEntriesFast();
		Long64_t nbytes = 0, nb = 0;
		for (Long64_t jentry=0; jentry<nentries;jentry++) {
			Long64_t ientry = LoadTree(jentry);
			if (ientry < 0) break;
			nb = fChain->GetEntry(jentry);   nbytes += nb;

			if((*this.*eventSelector2)() == false) continue;
			if((*this.*muonSelector2)(muon2) == false) continue;
			h2_bg->Fill(MuIso[muon2], MuPt[muon2], scale);
		}
	}

	// fSamples[sample1].tree->Project("h2_bg", Form("MuPt[%d]:MuIso[%d]", muon1, muon1), c1);
	// fSamples[sample2].tree->Project("h2_sig", Form("MuPt[%d]:MuIso[%d]", muon2, muon2), c2);

	TLatex *lat = new TLatex();
	lat->SetNDC(kTRUE);
	lat->SetTextColor(kBlack);
	lat->SetTextSize(0.06);

	TCanvas *c_temp = new TCanvas("IsoVsPt", "Isolating in Pt bins", 0, 0, 1200, 800);
	c_temp->Divide(3,2);
	// TCanvas *c_temp = new TCanvas("IsoVsPt", "Isolating in Pt bins", 0, 0, 1200, 1200);
	// c_temp->Divide(3,3);
	// c_temp->cd(9);
	// if(logy) gPad->SetLogz(1);
	// h2_sig->DrawCopy("colz");
	// lat->DrawLatex(0.11,0.92, fSamples[TTbar].sname);

	for(size_t i = 1; i <= gNPt2bins; ++i){
		c_temp->cd(i);
		gStyle->SetOptStat(1111);
		TH1D *h1 = h2_sig->ProjectionX(Form("h2_sigx_%d", i), i, i);	
		TH1D *h2 = h2_bg ->ProjectionX(Form("h2_bgx_%d",  i), i, i);
		h1->SetXTitle(convertVarName("MuIso[0]"));
		h1->SetLineWidth(2);
		h1->SetFillColor(15);
		h1->SetFillStyle(1001);
		h2->SetLineWidth(2);
		h2->SetLineColor(kBlue);
		h2->SetFillColor(kBlue);
		h2->SetFillStyle(3004);

		if(logy) gPad->SetLogy(1);
		gPad->SetFillStyle(0);
		h1->Sumw2();
		h2->Sumw2();

		// Scaling
		if(h1->GetEntries() > 0 ) h1->Scale(1.0/h1->Integral());
		if(h2->GetEntries() > 0 ) h2->Scale(1.0/h2->Integral());

		// setPlottingRange(h1, h2);

		// Determine plotting range
		double max1 = h1->GetMaximum();
		double max2 = h2->GetMaximum();
		double max  = (max1>max2)?max1:max2;
		if(logy) max = 5*max;
		else max = 1.05*max;
		h1->SetMaximum(max);
		h2->SetMaximum(max);
		if(!logy){
			h1->SetMinimum(0.0);
			h2->SetMinimum(0.0);
		}

		TLegend *leg = new TLegend(0.45,0.75,0.65,0.88);
		if(i == 1){
			leg->AddEntry(h1, fSamples[TTbar].sname,"f");
			leg->AddEntry(h2, fSamples[MuB].sname,"f");
			leg->SetFillStyle(0);
			leg->SetTextFont(42);
			leg->SetBorderSize(0);
		}

		TH1D *h1_temp = (TH1D*)h1->DrawCopy("hist");
		TH1D *h2_temp = (TH1D*)h2->DrawCopy("histsames");
		h1_temp->SetName(Form("h1_%d",i));
		h2_temp->SetName(Form("h2_%d",i));
		gPad->Update();
		TPaveStats *s1 = (TPaveStats*)h1_temp->GetListOfFunctions()->FindObject("stats");
		TPaveStats *s2 = (TPaveStats*)h2_temp->GetListOfFunctions()->FindObject("stats");
		s2->SetTextColor(kBlue); s2->SetLineColor(kBlue);
		s2->SetY1NDC(s1->GetY1NDC() - (s1->GetY2NDC() - s1->GetY1NDC()));
		s2->SetY2NDC(s1->GetY1NDC());

		if(i==1) leg->Draw();

		// double max1 = h1->GetYaxis()->GetXmax();
		// double max2 = h2->GetYaxis()->GetXmax();
		// double max  = (max1>max2)?max1:max2;
		double min1 = h1->GetYaxis()->GetXmin();
		double min2 = h2->GetYaxis()->GetXmin();
		double min  = (min1<min2)?min1:min2;

		TLine *l1 = new TLine(0.15, min, 0.15, max);
		l1->SetLineColor(kRed);
		l1->SetLineWidth(2);
		l1->Draw();

		lat->SetTextColor(kBlack);
		lat->SetTextSize(0.05);
		lat->DrawLatex(0.11,0.92, Form("p_{T}(#mu) %3.0f - %3.0f GeV", gPt2bins[i-1], gPt2bins[i]));

		int bin0 = h1->FindBin(0.00);
		int bin15 = h1->FindBin(0.15);
		int bin1 = h1->FindBin(1.00);
		float f1 = h1->Integral(bin0, bin15) / h1->Integral(bin0, bin1);
		float f2 = h2->Integral(bin0, bin15) / h2->Integral(bin0, bin1);
		lat->SetTextSize(0.04);
		lat->DrawLatex(0.55,0.905, Form("ratio = %4.2f", f1));
		lat->SetTextColor(kBlue);
		lat->DrawLatex(0.55,0.945, Form("ratio = %4.2f", f2));

		gPad->RedrawAxis();
	}
	Util::PrintNoEPS(c_temp, outputname, fOutputDir, fOutputFile);
}

//____________________________________________________________________________
void MuonPlotter::makeIsoVsNJetsPlot(int sample1, int muon1, TCut c1, int sample2, int muon2, TCut c2, TString outputname, bool logy){
	const int nbins = 20;
	const int Nnjetbins = 4;
	const double njetbins[Nnjetbins+1] = {0.,1.,2.,3.,4.};
	TH2D *h2_sig = new TH2D("h2_sig", "Isolation vs NJets for fake muons in signal", nbins, 0., 1., Nnjetbins, njetbins);
	TH2D *h2_bg  = new TH2D("h2_bg", "Isolation vs NJets for muons in background",   nbins, 0., 1., Nnjetbins, njetbins);
	h2_sig->SetXTitle(convertVarName("MuIso[0]"));
	h2_sig->SetYTitle("NJets");
	h2_bg->SetXTitle(convertVarName("MuIso[0]"));
	h2_bg->SetYTitle("NJets");

	fSamples[sample1].tree->Project("h2_bg",  Form("NJets:MuIso[%d]", muon1), c1);
	fSamples[sample2].tree->Project("h2_sig", Form("NJets:MuIso[%d]", muon2), c2);


	TLatex *lat = new TLatex();
	lat->SetNDC(kTRUE);
	lat->SetTextSize(0.06);
	lat->SetTextColor(kBlack);

	TCanvas *c_temp = new TCanvas("IsoVsNJets", "Isolating in bins of jet multiplicity", 0, 0, 1200, 800);
	c_temp->Divide(3,2);
	c_temp->cd(5);
	if(logy) gPad->SetLogz(1);
	h2_sig->DrawCopy("colz");
	lat->DrawLatex(0.11,0.92, fSamples[sample2].sname);
	c_temp->cd(6);
	if(logy) gPad->SetLogz(1);
	h2_bg->DrawCopy("colz");
	lat->DrawLatex(0.11,0.92, fSamples[sample1].sname);

	for(size_t i = 1; i <= Nnjetbins; ++i){
		c_temp->cd(i);
		gStyle->SetOptStat(1111);
		int lastbin = i;
		if(i == Nnjetbins) lastbin = i+1;
		TH1D *h1 = h2_bg->ProjectionX("h2_bgx",i, lastbin);
		TH1D *h2 = h2_sig->ProjectionX("h2_sigx",i, lastbin);
		h1->SetXTitle(convertVarName("MuIso[0]"));
		h1->SetLineWidth(2);
		h1->SetFillColor(15);
		h1->SetFillStyle(1001);
		h2->SetLineWidth(2);
		h2->SetLineColor(kBlue);
		h2->SetFillColor(kBlue);
		h2->SetFillStyle(3004);

		if(logy) gPad->SetLogy(1);
		gPad->SetFillStyle(0);
		h1->Sumw2();
		h2->Sumw2();
		if(h1->GetEntries() > 0 ) h1->Scale(1.0/h1->Integral());
		if(h2->GetEntries() > 0 ) h2->Scale(1.0/h2->Integral());

		// Determine plotting range
		double max1 = h1->GetMaximum();
		double max2 = h2->GetMaximum();
		double max  = (max1>max2)?max1:max2;
		if(logy) max = 5*max;
		else max = 1.05*max;
		h1->SetMaximum(max);
		h2->SetMaximum(max);
		if(!logy){
			h1->SetMinimum(0.0);
			h2->SetMinimum(0.0);
		}

		TLegend *leg = new TLegend(0.55,0.75,0.75,0.88);
		if(i == 1){
			leg->AddEntry(h1, fSamples[sample1].sname,"f");
			leg->AddEntry(h2, fSamples[sample2].sname,"f");
			leg->SetFillStyle(0);
			leg->SetTextFont(42);
			leg->SetBorderSize(0);
		}

		TH1D *h1_temp = (TH1D*)h1->DrawCopy("hist");
		TH1D *h2_temp = (TH1D*)h2->DrawCopy("histsames");
		gPad->Update();
		TPaveStats *s1 = (TPaveStats*)h1_temp->GetListOfFunctions()->FindObject("stats");
		TPaveStats *s2 = (TPaveStats*)h2_temp->GetListOfFunctions()->FindObject("stats");
		s2->SetTextColor(kBlue); s2->SetLineColor(kBlue);
		s2->SetY1NDC(s1->GetY1NDC() - (s1->GetY2NDC() - s1->GetY1NDC()));
		s2->SetY2NDC(s1->GetY1NDC());

		if(i==1) leg->Draw();

		double min1 = h1->GetYaxis()->GetXmin();
		double min2 = h2->GetYaxis()->GetXmin();
		double min  = (min1<min2)?min1:min2;

		TLine *l1 = new TLine(0.15, min, 0.15, max);
		l1->SetLineColor(kRed);
		l1->SetLineWidth(2);
		l1->Draw();

		lat->SetTextColor(kBlack);
		lat->SetTextSize(0.05);
		if(i < 4) lat->DrawLatex(0.11,0.92, Form("NJets = %1.0f",  njetbins[i-1]));
		else      lat->DrawLatex(0.11,0.92, Form("NJets >= %1.0f", njetbins[i-1]));

		int bin0 = h1->FindBin(0.00);
		int bin15 = h1->FindBin(0.15);
		int bin1 = h1->FindBin(1.00);
		float f1 = h1->Integral(bin0, bin15) / h1->Integral(bin0, bin1);
		float f2 = h2->Integral(bin0, bin15) / h2->Integral(bin0, bin1);
		lat->SetTextSize(0.04);
		lat->DrawLatex(0.55,0.905, Form("ratio = %4.2f", f1));
		lat->SetTextColor(kBlue);
		lat->DrawLatex(0.55,0.945, Form("ratio = %4.2f", f2));

		gPad->RedrawAxis();
	}
	Util::PrintNoEPS(c_temp, outputname, fOutputDir, fOutputFile);
}

//____________________________________________________________________________
void MuonPlotter::produceRatio(int sample, int muon, bool(MuonPlotter::*eventSelector)(), bool(MuonPlotter::*muonSelector)(int), TH2D *&h_2d, TH1D *&h_pt, TH1D *&h_eta, bool output){
	vector<int> samples; samples.push_back(sample);
	produceRatio(samples, muon, eventSelector, muonSelector, h_2d, h_pt, h_eta, output);
}
void MuonPlotter::produceRatio(vector<int> samples, int muon, bool(MuonPlotter::*eventSelector)(), bool(MuonPlotter::*muonSelector)(int), TH2D *&h_2d, TH1D *&h_pt, TH1D *&h_eta, bool output){
// Base function for production of all ratios
/*
	TODO Fix treatment of statistical errors and luminosity scaling here!
*/
	gStyle->SetOptStat(0);
	h_2d->Sumw2();
	h_pt->Sumw2();
	h_eta->Sumw2();

	TH2D *H_ntight = new TH2D("NTight", "NTight Muons", h_2d->GetNbinsX(), h_2d->GetXaxis()->GetXbins()->GetArray(), h_2d->GetNbinsY(),  h_2d->GetYaxis()->GetXbins()->GetArray());
	TH2D *H_nloose = new TH2D("NLoose", "NLoose Muons", h_2d->GetNbinsX(), h_2d->GetXaxis()->GetXbins()->GetArray(), h_2d->GetNbinsY(),  h_2d->GetYaxis()->GetXbins()->GetArray());
	H_ntight->Sumw2();
	H_nloose->Sumw2();

	if(fVerbose>2) cout << "---------------" << endl;
	for(size_t i = 0; i < samples.size(); ++i){
		int sample = samples[i];
		TTree *tree = fSamples[sample].tree;
		if(fVerbose>2) cout << "Producing ratios for " << fSamples[sample].sname << endl;
		tree->ResetBranchAddresses();
		Init(tree);
		if (fChain == 0) return;
		Long64_t nentries = fChain->GetEntriesFast();

		float scale = fLumiNorm / fSamples[sample].lumi;

		Long64_t nbytes = 0, nb = 0;
		for (Long64_t jentry=0; jentry<nentries;jentry++) {
			Long64_t ientry = LoadTree(jentry);
			if (ientry < 0) break;
			nb = fChain->GetEntry(jentry);   nbytes += nb;

			if((*this.*eventSelector)() == false) continue;
			if((*this.*muonSelector)(muon) == false) continue;

			if(isLooseMuon(muon)) H_nloose->Fill(MuPt[muon], MuEta[muon], scale); // Tight or loose
			if(isTightMuon(muon)) H_ntight->Fill(MuPt[muon], MuEta[muon], scale); // Tight
		}

		if(fVerbose>2) cout << " Tight entries so far: " << H_ntight->GetEntries() << " / " << H_ntight->Integral() << endl;
		if(fVerbose>2) cout << " Loose entries so far: " << H_nloose->GetEntries() << " / " << H_nloose->Integral() << endl;
		if(fVerbose>2) cout << "  Ratio: " << (double)H_ntight->GetEntries()/(double)H_nloose->GetEntries() << endl;
	}
	h_2d->Divide(H_ntight, H_nloose);

	TH1D *hmuloosept  = H_nloose->ProjectionX();
	TH1D *hmulooseeta = H_nloose->ProjectionY();
	TH1D *hmutightpt  = H_ntight->ProjectionX();
	TH1D *hmutighteta = H_ntight->ProjectionY();
	h_pt ->Divide(hmutightpt, hmuloosept);
	h_eta->Divide(hmutighteta, hmulooseeta);

	// TGraphAsymmErrors *asymErrors = new TGraphAsymmErrors(h_pt);
	// asymErrors->SetName("asymErrors");
	// asymErrors->BayesDivide(hmutightpt, hmuloosept);

	// TCanvas *c1 = makeCanvas("asymErrors");
	// c1->cd();
	// h_pt->SetLineColor(kBlue);
	// h_pt->Draw();
	// asymErrors->Draw("same");
	// Util::PrintNoEPS(c1, "test", fOutputDir);

	h_pt ->SetXTitle(convertVarName("MuPt[1]"));
	h_eta->SetXTitle(convertVarName("MuEta[1]"));
	h_pt ->SetYTitle("# Tight / # Loose");
	h_eta->SetYTitle("# Tight / # Loose");
	h_2d->SetXTitle(convertVarName("MuPt[1]"));
	h_2d->SetYTitle(convertVarName("MuEta[1]"));
	delete H_ntight, H_nloose, hmuloosept, hmulooseeta, hmutightpt, hmutighteta;
	TString name = "";
	for(size_t i = 0; i < samples.size(); ++i){
		int sample = samples[i];
		name += h_2d->GetName();
		name += "_";
		name += fSamples[sample].sname;
	}
	if(output){
		printHisto(h_2d,  TString("Ratio")    + name, "Fake Ratio vs pt/eta", "colz");
		printHisto(h_pt,  TString("RatioPt")  + name, "Fake Ratio vs pt",     "PE1");
		printHisto(h_eta, TString("RatioEta") + name, "Fake Ratio vs eta",    "PE1");
	}
}

//____________________________________________________________________________
vector<double> MuonPlotter::produceRatio(vector<int> samples, int muon, bool(MuonPlotter::*eventSelector)(), bool(MuonPlotter::*muonSelector)(int)){
// Will return two numbers, the ratio and the error
	vector<double> ratios;
	TH1D *H_ntight = new TH1D("NTight", "NTight Muons",1, 0, 10);
	TH1D *H_nloose = new TH1D("NLoose", "NLoose Muons",1, 0, 10);
	TH1D *H_ratio = new TH1D("TightLooseRatio", "Tight to Loose Ratio",1, 0, 10);
	H_ntight->Sumw2();
	H_nloose->Sumw2();

	if(fVerbose>2) cout << "---------------" << endl;
	for(size_t i = 0; i < samples.size(); ++i){
		int sample = samples[i];
		TTree *tree = fSamples[sample].tree;
		if(fVerbose>2) cout << "Producing ratios for " << fSamples[sample].sname << endl;
		tree->ResetBranchAddresses();
		Init(tree);
		if (fChain == 0) return ratios;
		Long64_t nentries = fChain->GetEntriesFast();
		float scale = fLumiNorm / fSamples[sample].lumi;
		Long64_t nbytes = 0, nb = 0;
		for (Long64_t jentry=0; jentry<nentries;jentry++) {
			Long64_t ientry = LoadTree(jentry);
			if (ientry < 0) break;
			nb = fChain->GetEntry(jentry);   nbytes += nb;

			if((*this.*eventSelector)() == false) continue;
			if((*this.*muonSelector)(muon) == false) continue;

			if(isLooseMuon(muon)) H_nloose->Fill(1, scale); // Tight or loose
			if(isTightMuon(muon)) H_ntight->Fill(1, scale); // Tight
		}
		if(fVerbose>2) cout << " Tight entries so far: " << H_ntight->GetEntries() << endl;
		if(fVerbose>2) cout << " Loose entries so far: " << H_nloose->GetEntries() << endl;
		if(fVerbose>2) cout << "  Ratio so far: " << (double)H_ntight->GetEntries()/(double)H_nloose->GetEntries() << endl;
	}
	H_ratio->Divide(H_ntight, H_nloose);
	ratios.push_back(H_ratio->GetBinContent(1));
	ratios.push_back(H_ratio->GetBinError(1));
	return ratios;
}

//____________________________________________________________________________
void MuonPlotter::fillfRatio(int sample, int muon){ fillfRatio(sample, muon, gNPt2bins, gPt2bins, gNEtabins, gEtabins); }
void MuonPlotter::fillfRatio(vector<int> samples, int muon){ fillfRatio(samples, muon, gNPt2bins, gPt2bins, gNEtabins, gEtabins); }
void MuonPlotter::fillfRatio(int sample, int muon, const int nptbins, const double *ptbins, const int netabins, const double *etabins){
	vector<int> samples; samples.push_back(sample);
	fillfRatio(samples, muon, nptbins, ptbins, netabins, etabins);
}
void MuonPlotter::fillfRatio(vector<int> samples, int muon, const int nptbins, const double *ptbins, const int netabins, const double *etabins){
	gStyle->SetOptStat(0);
	fH2D_fRatio    = new TH2D("fRatio", "Ratio of tight to loose Muons vs Pt vs Eta", nptbins, ptbins, netabins, etabins);
	fH1D_fRatioPt  = new TH1D("fRatioPt",  "Ratio of tight to loose Muons vs Pt", nptbins, ptbins);
	fH1D_fRatioEta = new TH1D("fRatioEta", "Ratio of tight to loose Muons vs Eta", netabins, etabins);

	fH1D_fRatioPt->SetXTitle(convertVarName("MuPt[0]"));
	fH1D_fRatioEta->SetXTitle(convertVarName("MuEta[0]"));
	fH2D_fRatio->SetXTitle(convertVarName("MuPt[0]"));
	fH2D_fRatio->SetYTitle(convertVarName("MuEta[0]"));
	fH1D_fRatioPt ->SetYTitle("# Tight / # Loose");
	fH1D_fRatioEta->SetYTitle("# Tight / # Loose");

	// produceRatio(samples, muon, &MuonPlotter::isSigSupEvent, &MuonPlotter::isLooseMuon, fH2D_fRatio, fH1D_fRatioPt, fH1D_fRatioEta, true);
	calculateRatio(samples, Muon, 1, fH2D_fRatio, fH1D_fRatioPt, fH1D_fRatioEta, false);
}

//____________________________________________________________________________
void MuonPlotter::fillpRatio(int sample, int muon){ fillpRatio(sample, muon, gNPt2bins, gPt2bins, gNEtabins, gEtabins); }
void MuonPlotter::fillpRatio(vector<int> samples, int muon){ fillpRatio(samples, muon, gNPt2bins, gPt2bins, gNEtabins, gEtabins); }
void MuonPlotter::fillpRatio(int sample, int muon, const int nptbins, const double *ptbins, const int netabins, const double *etabins){
	vector<int> samples; samples.push_back(sample);
	fillpRatio(samples, muon, nptbins, ptbins, netabins, etabins);
}
void MuonPlotter::fillpRatio(vector<int> samples, int muon, const int nptbins, const double *ptbins, const int netabins, const double *etabins){
	gStyle->SetOptStat(0);
	// This is supposed to be run on a ZJets selection!
	fH2D_pRatio    = new TH2D("pRatio", "Ratio of tight to loose Muons vs Pt vs Eta", nptbins, ptbins, netabins, etabins);
	fH1D_pRatioPt  = new TH1D("pRatioPt",  "Ratio of tight to loose Muons vs Pt", nptbins, ptbins);
	fH1D_pRatioEta = new TH1D("pRatioEta", "Ratio of tight to loose Muons vs Eta", netabins, etabins);
	fH1D_pRatioPt->SetXTitle(convertVarName("MuPt[0]"));
	fH1D_pRatioEta->SetXTitle(convertVarName("MuEta[0]"));
	fH2D_pRatio->SetXTitle(convertVarName("MuPt[0]"));
	fH2D_pRatio->SetYTitle(convertVarName("MuEta[0]"));
	fH1D_pRatioPt ->SetYTitle("# Tight / # Loose");
	fH1D_pRatioEta->SetYTitle("# Tight / # Loose");

	// produceRatio(samples, muon, &MuonPlotter::isZMuMuEvent, &MuonPlotter::isLooseMuon, fH2D_pRatio, fH1D_pRatioPt, fH1D_pRatioEta, true);
	calculateRatio(samples, Muon, 2, fH2D_pRatio, fH1D_pRatioPt, fH1D_pRatioEta, false);
}

//____________________________________________________________________________
void MuonPlotter::plotRatio(int sample, int muon, bool(MuonPlotter::*eventSelector)(), bool(MuonPlotter::*muonSelector)(int), TString tag){
	vector<int> samples; samples.push_back(sample);
	plotRatio(samples, muon, eventSelector, muonSelector, tag);
}
void MuonPlotter::plotRatio(vector<int> samples, int muon, bool(MuonPlotter::*eventSelector)(), bool(MuonPlotter::*muonSelector)(int), TString tag){
	gStyle->SetOptStat(0);
	TH2D *h_2d  = new TH2D("Ratio",    Form("Ratio of tight to loose Muons vs Pt vs Eta : %s", tag.Data()), gNPtbins, gPtbins, gNEtabins, gEtabins);
	TH1D *h_pt  = new TH1D("RatioPt",  Form("Ratio of tight to loose Muons vs Pt : %s",        tag.Data()), gNPtbins, gPtbins);
	TH1D *h_eta = new TH1D("RatioEta", Form("Ratio of tight to loose Muons vs Eta : %s",       tag.Data()), gNEtabins, gEtabins);
	h_pt->SetXTitle(convertVarName("MuPt[0]"));
	h_eta->SetXTitle(convertVarName("MuEta[0]"));
	h_2d->SetXTitle(convertVarName("MuPt[0]"));
	h_2d->SetYTitle(convertVarName("MuEta[0]"));
	h_pt ->SetYTitle("# Tight / # Loose");
	h_eta->SetYTitle("# Tight / # Loose");
	h_pt->GetYaxis()->SetTitleOffset(1.2);
	h_eta->GetYaxis()->SetTitleOffset(1.2);

	produceRatio(samples, muon, eventSelector, muonSelector, h_2d, h_pt, h_eta, true);
}

//____________________________________________________________________________
TH1D* MuonPlotter::fillRatioPt(int sample, int muon, bool(MuonPlotter::*eventSelector)(), bool(MuonPlotter::*muonSelector)(int), bool output){
	vector<int> samples; samples.push_back(sample);
	return fillRatioPt(samples, muon, eventSelector, muonSelector, output);
}
TH1D* MuonPlotter::fillRatioPt(vector<int> samples, int muon, bool(MuonPlotter::*eventSelector)(), bool(MuonPlotter::*muonSelector)(int), bool output){
	gStyle->SetOptStat(0);
	TH2D *h_2d  = new TH2D("Ratio",    "Ratio of tight to loose Muons vs Pt vs Eta", gNPt2bins, gPt2bins, gNEtabins, gEtabins);
	TH1D *h_pt  = new TH1D("RatioPt",  "Ratio of tight to loose Muons vs Pt",        gNPt2bins, gPt2bins);
	TH1D *h_eta = new TH1D("RatioEta", "Ratio of tight to loose Muons vs Eta",       gNEtabins, gEtabins);

	h_pt->SetXTitle(convertVarName("MuPt[0]"));
	h_pt ->SetYTitle("# Tight / # Loose");
	h_pt->GetYaxis()->SetTitleOffset(1.2);

	produceRatio(samples, muon, eventSelector, muonSelector, h_2d, h_pt, h_eta, output);
	return h_pt;
}
TH1D* MuonPlotter::fillRatioPt(vector<int> samples, int muon, bool(MuonPlotter::*eventSelector)(), bool(MuonPlotter::*muonSelector)(int), const int nptbins, const double* ptbins, const int netabins, const double* etabins, bool output){
	gStyle->SetOptStat(0);
	TH2D *h_2d  = new TH2D("Ratio",    "Ratio of tight to loose Muons vs Pt vs Eta", nptbins, ptbins, netabins, etabins);
	TH1D *h_pt  = new TH1D("RatioPt",  "Ratio of tight to loose Muons vs Pt",        nptbins, ptbins);
	TH1D *h_eta = new TH1D("RatioEta", "Ratio of tight to loose Muons vs Eta",       netabins, etabins);

	h_pt->SetXTitle(convertVarName("MuPt[0]"));
	h_pt ->SetYTitle("# Tight / # Loose");
	h_pt->GetYaxis()->SetTitleOffset(1.2);

	produceRatio(samples, muon, eventSelector, muonSelector, h_2d, h_pt, h_eta, output);
	return h_pt;
}

//____________________________________________________________________________
TH2D* MuonPlotter::fillRatio(int sample, int muon, bool(MuonPlotter::*eventSelector)(), bool(MuonPlotter::*muonSelector)(int), const int nptbins, const double* ptbins, const int netabins, const double* etabins){
	vector<int> samples; samples.push_back(sample);
	return fillRatio(samples, muon, eventSelector, muonSelector, nptbins, ptbins, netabins, etabins);
}
TH2D* MuonPlotter::fillRatio(vector<int> samples, int muon, bool(MuonPlotter::*eventSelector)(), bool(MuonPlotter::*muonSelector)(int), const int nptbins, const double* ptbins, const int netabins, const double* etabins){
	gStyle->SetOptStat(0);
	TH2D *h_2d  = new TH2D("Ratio",    "Ratio of tight to loose Muons vs Pt vs Eta", nptbins, ptbins, netabins, etabins);
	TH1D *h_pt  = new TH1D("RatioPt",  "Ratio of tight to loose Muons vs Pt",        nptbins, ptbins);
	TH1D *h_eta = new TH1D("RatioEta", "Ratio of tight to loose Muons vs Eta",       netabins, etabins);

	h_2d->SetXTitle(convertVarName("MuPt[0]"));
	h_2d->SetYTitle(convertVarName("MuEta[0]"));

	produceRatio(samples, muon, eventSelector, muonSelector, h_2d, h_pt, h_eta, false);
	return h_2d;
}

//____________________________________________________________________________
TH1D* MuonPlotter::fillRatioPt(int sample, int forp, bool output){
	vector<int> samples; samples.push_back(sample);
	return fillRatioPt(samples, forp, output);
}
TH1D* MuonPlotter::fillRatioPt(vector<int> samples, int forp, bool output){
	gStyle->SetOptStat(0);
	TH2D *h_2d  = new TH2D("Ratio",    "Ratio of tight to loose Muons vs Pt vs Eta", gNPt2bins, gPt2bins, gNEtabins, gEtabins);
	TH1D *h_pt  = new TH1D("RatioPt",  "Ratio of tight to loose Muons vs Pt",        gNPt2bins, gPt2bins);
	TH1D *h_eta = new TH1D("RatioEta", "Ratio of tight to loose Muons vs Eta",       gNEtabins, gEtabins);

	h_pt->SetXTitle(convertVarName("MuPt[0]"));
	h_pt->SetYTitle("# Tight / # Loose");
	h_pt->GetYaxis()->SetTitleOffset(1.2);

	calculateRatio(samples, Muon, forp, h_2d, h_pt, h_eta, output);
	return h_pt;
};

//____________________________________________________________________________
void MuonPlotter::calculateRatio(vector<int> samples, gChannel chan, int forp, TH2D*& h_2d, bool output){
	TH1D *h_dummy1 = new TH1D("dummy1", "dummy1", 1, 0.,1.);
	TH1D *h_dummy2 = new TH1D("dummy2", "dummy2", 1, 0.,1.);
	calculateRatio(samples, chan, forp, h_2d, h_dummy1, h_dummy2, output);
	delete h_dummy1, h_dummy2;
}
void MuonPlotter::calculateRatio(vector<int> samples, gChannel chan, int forp, TH2D*& h_2d, TH1D*& h_pt, TH1D*&h_eta, bool output){
/*
TODO Fix treatment of statistical errors and luminosity scaling here!
*/
	gStyle->SetOptStat(0);
	h_2d->Sumw2();
	h_pt->Sumw2();
	h_eta->Sumw2();

	TH2D *H_ntight = new TH2D("NTight", "NTight Muons", h_2d->GetNbinsX(), h_2d->GetXaxis()->GetXbins()->GetArray(), h_2d->GetNbinsY(),  h_2d->GetYaxis()->GetXbins()->GetArray());
	TH2D *H_nloose = new TH2D("NLoose", "NLoose Muons", h_2d->GetNbinsX(), h_2d->GetXaxis()->GetXbins()->GetArray(), h_2d->GetNbinsY(),  h_2d->GetYaxis()->GetXbins()->GetArray());
	H_ntight->Sumw2();
	H_nloose->Sumw2();

	if(fVerbose>2) cout << "---------------" << endl;
	for(size_t i = 0; i < samples.size(); ++i){
		int sample = samples[i];
		TTree *tree = fSamples[sample].tree;
		if(fVerbose>2) cout << "Calculating ratios for " << fSamples[sample].sname << endl;

		float scale = fLumiNorm / fSamples[sample].lumi;

		channel *cha;
		if(chan == Muon)     cha = &fSamples[sample].mumu;
		if(chan == Electron) cha = &fSamples[sample].ee;

		if(forp == 1){
			H_ntight->Add(cha->fhistos.h_ntight, scale);
			H_nloose->Add(cha->fhistos.h_nloose, scale);			
		}
		if(forp == 2){
			H_ntight->Add(cha->phistos.h_ntight, scale);
			H_nloose->Add(cha->phistos.h_nloose, scale);			
		}

		if(fVerbose>2) cout << " Tight entries so far: " << H_ntight->GetEntries() << " / " << H_ntight->Integral() << endl;
		if(fVerbose>2) cout << " Loose entries so far: " << H_nloose->GetEntries() << " / " << H_nloose->Integral() << endl;
		if(fVerbose>2) cout << "  Ratio so far       : " << (double)H_ntight->GetEntries()/(double)H_nloose->GetEntries() << endl;
	}
	h_2d->Divide(H_ntight, H_nloose);

	TH1D *hmuloosept  = H_nloose->ProjectionX();
	TH1D *hmulooseeta = H_nloose->ProjectionY();
	TH1D *hmutightpt  = H_ntight->ProjectionX();
	TH1D *hmutighteta = H_ntight->ProjectionY();

	h_pt ->Divide(hmutightpt, hmuloosept);
	h_eta->Divide(hmutighteta, hmulooseeta);

	// TGraphAsymmErrors *asymErrors = new TGraphAsymmErrors(h_pt);
	// asymErrors->SetName("asymErrors");
	// asymErrors->BayesDivide(hmutightpt, hmuloosept);

	// TCanvas *c1 = makeCanvas("asymErrors");
	// c1->cd();
	// h_pt->SetLineColor(kBlue);
	// h_pt->Draw();
	// asymErrors->Draw("same");
	// Util::PrintNoEPS(c1, "test", fOutputDir);

	h_pt ->SetXTitle(convertVarName("MuPt[1]"));
	h_eta->SetXTitle(convertVarName("MuEta[1]"));
	h_pt ->SetYTitle("# Tight / # Loose");
	h_eta->SetYTitle("# Tight / # Loose");
	h_2d->SetXTitle(convertVarName("MuPt[1]"));
	h_2d->SetYTitle(convertVarName("MuEta[1]"));
	delete H_ntight, H_nloose, hmuloosept, hmulooseeta, hmutightpt, hmutighteta;
	TString name = "";
	for(size_t i = 0; i < samples.size(); ++i){
		int sample = samples[i];
		name += h_2d->GetName();
		name += "_";
		name += fSamples[sample].sname;
	}
	if(output){
		printHisto(h_2d,  TString("Ratio")    + name, "Fake Ratio vs pt/eta", "colz");
		printHisto(h_pt,  TString("RatioPt")  + name, "Fake Ratio vs pt",     "PE1");
		printHisto(h_eta, TString("RatioEta") + name, "Fake Ratio vs eta",    "PE1");
	}
}
void MuonPlotter::calculateRatio(vector<int> samples, gChannel chan, int forp, float &ratio, float &ratioe, bool output){
	double ntight(0.), nloose(0.);
	double ntighte2(0.), nloosee2(0.);
	for(size_t i = 0; i < samples.size(); ++i){
		int index = samples[i];
		float scale = fLumiNorm/fSamples[index].lumi; // Normalize all
		channel *cha;
		if(chan == Muon)     cha = &fSamples[index].mumu;
		if(chan == Electron) cha = &fSamples[index].ee;
		if(forp == 1){
			ntight += scale * cha->numbers.nsst;
			nloose += scale * cha->numbers.nssl;
		}
		if(forp == 2){
			ntight += scale * cha->numbers.nzt;
			nloose += scale * cha->numbers.nzl;
		}
	}
	ratio = ntight/nloose;
	ratioe = TMath::Sqrt( ntight*(1.0-ntight/nloose) ) / nloose;                  // Binomial
	// ratioe = TMath::Sqrt( ntight*ntight*(nloose+ntight)/(nloose*nloose*nloose) ); // Poissonian
}

//____________________________________________________________________________
void MuonPlotter::makeSSPredictionPlots(vector<int> samples){
	// Need filled ratios before calling this function!

	// Prediction: /////////////////////////////////////////////////////////////////////
	TH1D *H_nsigpred = new TH1D("Nsigpred", "Predicted N_sig in Pt1 bins",       gNPtbins,  gPtbins);
	TH1D *H_nfppred  = new TH1D("Nfppred",  "Predicted N_fp in Pt1 bins",        gNPtbins,  gPtbins);
	TH1D *H_nffpred  = new TH1D("Nffpred",  "Predicted N_ff in Pt1 bins",        gNPtbins,  gPtbins);
	TH1D *H_nFpred   = new TH1D("NFpred",   "Total predicted fakes in Pt1 bins", gNPtbins,  gPtbins);
	bool output = false;
	for(size_t i = 0; i < samples.size(); ++i){
		int index = samples[i];
		float scale = fLumiNorm/fSamples[samples[i]].lumi;
		vector<TH1D*> prediction = NsigPredFromFPRatios(index, output);
		H_nsigpred->Add(prediction[0], scale);
		H_nfppred ->Add(prediction[1], scale);
		H_nffpred ->Add(prediction[2], scale);
	}

	// Observation: ////////////////////////////////////////////////////////////////////
	TH1D *H_nsigobs  = new TH1D("Nsigobs", "Observed N_sig in Pt1 bins",  gNPtbins,  gPtbins);
	vector<int> lm0sample; lm0sample.push_back(LM0);
	NObs(H_nsigobs, lm0sample, &MuonPlotter::isGenMatchedSUSYDiLepEvent);

	TH1D *H_nt2obs  = new TH1D("Nt2obs", "Observed N_t2 in Pt1 bins",  gNPtbins,  gPtbins);
	NObs(H_nt2obs, samples);
	// NObs(H_nt2obs, samples, &MuonPlotter::isSSTTEvent);

	TH1D *H_nt2obsttbar = new TH1D("Nt2obsttbar", "Observed N_t2 in Pt1 bins, ttbar only",  gNPtbins,  gPtbins);
	vector<int> ttbarsample; ttbarsample.push_back(TTbar);
	NObs(H_nt2obsttbar, fMCBG);
	// NObs(H_nt2obsttbar, ttbarsample, &MuonPlotter::isSSTTEvent);	

	// Output
	H_nsigobs->SetXTitle(convertVarName("MuPt[0]"));
	H_nsigobs->SetYTitle(Form("Events / %2.0f GeV", fBinWidthScale));

	H_nFpred->Add(H_nfppred);
	H_nFpred->Add(H_nffpred);
	H_nFpred->SetXTitle(convertVarName("MuPt[0]"));
	H_nFpred->SetYTitle(Form("Events / %2.0f GeV", fBinWidthScale));
	H_nt2obs->SetXTitle(convertVarName("MuPt[0]"));
	H_nt2obs->SetYTitle(Form("Events / %2.0f GeV", fBinWidthScale));
	H_nt2obsttbar->SetXTitle(convertVarName("MuPt[0]"));
	H_nt2obsttbar->SetYTitle(Form("Events / %2.0f GeV", fBinWidthScale));

	// Normalize to binwidth
	H_nsigpred    = normHistBW(H_nsigpred,    fBinWidthScale);
	H_nsigobs     = normHistBW(H_nsigobs,     fBinWidthScale);
	H_nfppred     = normHistBW(H_nfppred,     fBinWidthScale);
	H_nffpred     = normHistBW(H_nffpred,     fBinWidthScale);
	H_nFpred      = normHistBW(H_nFpred,      fBinWidthScale);
	H_nt2obs      = normHistBW(H_nt2obs,      fBinWidthScale);
	H_nt2obsttbar = normHistBW(H_nt2obsttbar, fBinWidthScale);

	H_nt2obs->SetFillColor(kBlue);
	H_nt2obs->SetLineColor(kBlue);
	H_nt2obs->SetFillStyle(3004);
	H_nt2obs->SetLineWidth(2);

	H_nsigpred->SetFillColor(8);
	H_nsigpred->SetMarkerColor(8);
	H_nsigpred->SetMarkerStyle(20);
	H_nsigpred->SetLineColor(8);
	H_nsigpred->SetLineWidth(2);

	H_nfppred->SetFillColor(kRed);
	H_nfppred->SetMarkerColor(kRed);
	H_nfppred->SetMarkerStyle(20);
	H_nfppred->SetLineColor(kRed);
	H_nfppred->SetLineWidth(2);

	H_nffpred->SetFillColor(13);
	H_nffpred->SetMarkerColor(13);
	H_nffpred->SetLineColor(13);
	H_nffpred->SetMarkerStyle(20);
	H_nffpred->SetLineWidth(2);

	vector<TH1D*> hists;
	hists.push_back(H_nt2obs);
	hists.push_back(H_nsigpred);
	hists.push_back(H_nfppred);
	hists.push_back(H_nffpred);
	// setPlottingRange(hists);

	plotOverlay4H(H_nt2obs, "N_{ t2}", H_nsigpred, "N_{ pp}" , H_nfppred, "N_{ fp}", H_nffpred, "N_{ ff}");

	H_nFpred->SetMinimum(0.);
	H_nt2obs->SetMinimum(0.);
	H_nsigobs->SetMinimum(0.);
	H_nsigpred->SetMinimum(0.);
	H_nt2obsttbar->SetMinimum(0.);
	// setPlottingRange(H_nsigobs, H_nsigpred);

	H_nsigpred->SetMarkerColor(kRed);
	H_nFpred->SetMarkerColor(kRed);
	H_nFpred->SetMarkerStyle(20);
	H_nsigobs->SetMaximum(14.);
	H_nsigpred->SetMaximum(14.);
	H_nsigpred->SetMinimum(0.);
	plotPredOverlay2HWithRatio(H_nsigobs, "Observed N_{sig}", H_nsigpred, "Predicted N_{sig}");
	plotPredOverlay3HWithRatio(H_nFpred, "Predicted Fakes", H_nt2obs,  "Obs. N_{t2} (QCD, t#bar{t}+jets, V+jets, LM0)", H_nt2obsttbar, "Obs. N_{t2} (QCD, t#bar{t}+jets, V+jets)", false, false);
}

//____________________________________________________________________________
void MuonPlotter::NObs(TH1D *&hist, vector<int> samples, bool(MuonPlotter::*eventSelector)()){
	hist->Sumw2();
	for(size_t i = 0; i < samples.size(); ++i){
		int index = samples[i];
		float scale = fLumiNorm / fSamples[index].lumi;

		TTree *tree = fSamples[index].tree;
		tree->ResetBranchAddresses();
		Init(tree);
		if (fChain == 0) return;
		Long64_t nentries = fChain->GetEntriesFast();
		Long64_t nbytes = 0, nb = 0;
		for (Long64_t jentry=0; jentry<nentries;jentry++) {
			Long64_t ientry = LoadTree(jentry);
			if (ientry < 0) break;
			nb = fChain->GetEntry(jentry);   nbytes += nb;

			if((*this.*eventSelector)() == false) continue;
			hist->Fill(MuPt[0], scale);
		}
	}	
}

//____________________________________________________________________________
void MuonPlotter::NObs(TH1D *&hist, vector<int> samples){
	hist->Sumw2();
	for(size_t i = 0; i < samples.size(); ++i){
		int index = samples[i];
		float scale = fLumiNorm / fSamples[index].lumi;
		channel cha = fSamples[index].mumu;
		TH1D *h_temp = cha.nthistos.h_nt2_pt;
		for(int j = 1; j <= h_temp->GetNbinsX(); ++j){
			double bincenter = h_temp->GetBinCenter(j);
			double content = h_temp->GetBinContent(j);
			int newbin = hist->FindBin(bincenter);

			double oldcontent = hist->GetBinContent(newbin);
			double olderror = hist->GetBinError(newbin);
			double newerror = TMath::Sqrt(olderror*olderror + content*scale*scale);
			double newcontent = oldcontent + content*scale;
			hist->SetBinContent(newbin, newcontent);
			hist->SetBinError(newbin, newerror);
			
		}
	}	
}

//____________________________________________________________________________
void MuonPlotter::makeIntPrediction(gChannel chan){
	if(chan == Muon)     makeIntPredictionMuMu();
	if(chan == Electron) makeIntPredictionEE();
	if(chan == EMu)      makeIntPredictionEMu();
}

//____________________________________________________________________________
void MuonPlotter::makeIntPredictionMuMu(){
	cout << "-----------------------------------" << endl;
	cout << " Producing prediction for MuMu channel" << endl << endl;
	bool data = true; // Use ratios from data or mc?

	// Which samples to use for nt2/nt1/nt0 input?
	// vector<int> inputsamples = fMCBG;
	// vector<int> inputsamples = fMCBGSig;
	vector<int> inputsamples = fMuData;

	// Which luminosity to use?
	// fLumiNorm = fSamples[MuA].lumi + fSamples[MuB].lumi;
	fLumiNorm = 22.;

	// Dummy binning
	const int nptbins = 1;
	const double ptbins[nptbins+1] = {5., 1000.};
	const int netabins = 1;
	const double etabins[netabins+1] = {-2.5, 2.5};

	// Fill the ratios
	TH2D *H_fratio_allmc = new TH2D("fRatioAllMC", "fRatio for all MC", nptbins, ptbins, netabins, etabins);
	TH2D *H_fratio_data  = new TH2D("fRatioData", "fRatio for Mu Data", nptbins, ptbins, netabins, etabins);
	TH2D *H_pratio_allmc = new TH2D("pRatioAllMC", "pRatio for all MC", nptbins, ptbins, netabins, etabins);
	TH2D *H_pratio_data  = new TH2D("pRatioData", "pRatio for Mu Data", nptbins, ptbins, netabins, etabins);

	TH2D *H_fratio_ttbar = fillRatio(TTbar, 0, &MuonPlotter::isGoodMuEvent, &MuonPlotter::isFakeTTbarMuon, nptbins, ptbins, netabins, etabins);
	TH2D *H_pratio_ttbar = fillRatio(TTbar, 0, &MuonPlotter::isGoodMuEvent, &MuonPlotter::isPromptTTbarMuon, nptbins, ptbins, netabins, etabins);
	H_fratio_ttbar->SetName("fRatioTTbar");
	H_pratio_ttbar->SetName("pRatioTTbar");

	float fratio_allmc(0.), fratio_allmc_e(0.);
	float fratio_data(0.), fratio_data_e(0.);
	float pratio_allmc(0.), pratio_allmc_e(0.);
	float pratio_data(0.), pratio_data_e(0.);

	float fratio_ttbar   = H_fratio_ttbar->GetBinContent(1,1);
	float fratio_ttbar_e = H_fratio_ttbar->GetBinError(1,1);

	float pratio_ttbar   = H_pratio_ttbar->GetBinContent(1,1);
	float pratio_ttbar_e = H_pratio_ttbar->GetBinError(1,1);

	calculateRatio(fMCBG,   Muon, 1, fratio_allmc, fratio_allmc_e);
	calculateRatio(fMuData, Muon, 1, fratio_data,  fratio_data_e);
	calculateRatio(fMCBG,   Muon, 2, pratio_allmc, pratio_allmc_e);
	calculateRatio(fMuData, Muon, 2, pratio_data,  pratio_data_e);

	H_fratio_allmc->SetBinContent(1,1,fratio_allmc);
	H_fratio_allmc->SetBinError  (1,1,fratio_allmc_e);
	H_fratio_data ->SetBinContent(1,1,fratio_data);
	H_fratio_data ->SetBinError  (1,1,fratio_data_e);

	H_pratio_allmc->SetBinContent(1,1,pratio_allmc);
	H_pratio_allmc->SetBinError  (1,1,pratio_allmc_e);
	H_pratio_data ->SetBinContent(1,1,pratio_data);
	H_pratio_data ->SetBinError  (1,1,pratio_data_e);

	cout << "  fRatio from all MC         = " << fratio_allmc << " +/- " << fratio_allmc_e << endl;
	cout << "  fRatio from ttbar genmatch = " << fratio_ttbar << " +/- " << fratio_ttbar_e << endl;
	cout << "  fRatio from data           = " << fratio_data  << " +/- " << fratio_data_e << endl;
	cout << " ------------------------------------" << endl;
	cout << "  pRatio from all MC         = " << pratio_allmc << " +/- " << pratio_allmc_e << endl;
	cout << "  pRatio from ttbar genmatch = " << pratio_ttbar << " +/- " << pratio_ttbar_e << endl;
	cout << "  pRatio from data           = " << pratio_data  << " +/- " << pratio_data_e << endl;
	cout << " ------------------------------------" << endl;

	// Add systematics to ratios
	float deviation(0.), olderror(0.), newerror(0.);
	deviation = fabs(fratio_allmc - fratio_ttbar);
	// Add to mc ratios
	olderror = H_fratio_allmc->GetBinError(1,1);
	newerror = sqrt(olderror*olderror + deviation*deviation);
	H_fratio_allmc->SetBinError(1,1,newerror);
	fratio_allmc_e = newerror;
	// Add to data ratios
	olderror = H_fratio_data->GetBinError(1,1);
	newerror = sqrt(olderror*olderror + deviation*deviation);
	H_fratio_data->SetBinError(1,1,newerror);
	fratio_data_e = newerror;
	cout << "  new fRatio (all MC)        = " << H_fratio_allmc->GetBinContent(1,1) << " +/- " << H_fratio_allmc->GetBinError(1,1) << endl;
	cout << "  new fRatio (data)          = " << H_fratio_data ->GetBinContent(1,1) << " +/- " << H_fratio_data ->GetBinError(1,1) << endl;

	deviation  = fabs(pratio_allmc - pratio_ttbar);
	// Add to mc ratios
	olderror = pratio_allmc_e;
	newerror = sqrt(olderror*olderror + deviation*deviation);
	H_pratio_allmc->SetBinError(1,1,newerror);
	pratio_allmc_e = newerror;
	// Add to data ratios
	olderror = pratio_data_e;
	newerror = sqrt(olderror*olderror + deviation*deviation);
	H_pratio_data->SetBinError(1,1,newerror);
	pratio_data_e = newerror;
	cout << "  new pRatio (all MC)        = " << H_pratio_allmc->GetBinContent(1,1) << " +/- " << H_pratio_allmc->GetBinError(1,1) << endl;
	cout << "  new pRatio (data)          = " << H_pratio_data ->GetBinContent(1,1) << " +/- " << H_pratio_data ->GetBinError(1,1) << endl;
	cout << " ------------------------------------" << endl;

	double nt2(0.), nt1(0.), nt0(0.);
	double nt2_e2(0.), nt1_e2(0.), nt0_e2(0.);
	for(size_t i = 0; i < inputsamples.size(); ++i){
		int index = inputsamples[i];
		float scale = fLumiNorm/fSamples[index].lumi; // Normalize all
		channel cha = fSamples[index].mumu;
		if(data) scale = 1.;
		nt2 += scale * cha.numbers.nt2;
		nt1 += scale * cha.numbers.nt10;
		nt0 += scale * cha.numbers.nt0;
		nt2_e2 += scale*scale * cha.numbers.nt2;
		nt1_e2 += scale*scale * cha.numbers.nt10;
		nt0_e2 += scale*scale * cha.numbers.nt0;
	}
	if(data){
		cout << "  Found " << nt2 << " events in Nt2" << endl;
		cout << "  Found " << nt1 << " events in Nt1" << endl;
		cout << "  Found " << nt0 << " events in Nt0" << endl;
	}
	if(!data){
		cout << "  Found " << nt2 << " +/- " << TMath::Sqrt(nt2_e2) << " events in Nt2" << endl;
		cout << "  Found " << nt1 << " +/- " << TMath::Sqrt(nt1_e2) << " events in Nt1" << endl;
		cout << "  Found " << nt0 << " +/- " << TMath::Sqrt(nt0_e2) << " events in Nt0" << endl;		
	}
	cout << " ------------------------------------" << endl;

	// Make prediction
	fFPRatios = new FPRatios();
	fFPRatios->SetVerbose(fVerbose);
	if(data){
		fFPRatios->SetMuFratios(H_fratio_data);
		fFPRatios->SetMuPratios(H_pratio_data);
	}
	else{
		fFPRatios->SetMuFratios(H_fratio_allmc);
		fFPRatios->SetMuPratios(H_pratio_allmc);
	}
	vector<double> nt;
	nt.push_back(nt0); nt.push_back(nt1); nt.push_back(nt2);
	fFPRatios->NevtTopol(0, 2, nt);

	vector<double> vpt, veta;
	vpt.push_back(30.); vpt.push_back(30.); // Fake pts and etas
	veta.push_back(0.); veta.push_back(0.);

	vector<double> nevFP = fFPRatios->NevtPass(vpt, veta);
	vector<double> nevFPEstat = fFPRatios->NevtPassErrStat();
	vector<double> nevFPEsyst = fFPRatios->NevtPassErrSyst();

	cout << "  Prediction for Npp: " << nevFP[0] << " +/- " << nevFPEstat[0] << " (stat) +/- " << nevFPEsyst[0] << " (syst)" << endl;
	cout << "  Prediction for Nfp: " << nevFP[1] << " +/- " << nevFPEstat[1] << " (stat) +/- " << nevFPEsyst[1] << " (syst)" << endl;
	cout << "  Prediction for Nff: " << nevFP[2] << " +/- " << nevFPEstat[2] << " (stat) +/- " << nevFPEsyst[2] << " (syst)" << endl;
	cout << "  Total fakes:        " << nevFP[1]+nevFP[2] << " +/- " << TMath::Sqrt(nevFPEstat[1]*nevFPEstat[1] + nevFPEstat[2]*nevFPEstat[2]) << " (stat) +/- " << TMath::Sqrt(nevFPEsyst[1]*nevFPEsyst[1] + nevFPEsyst[2]*nevFPEsyst[2])<< " (syst)" << endl;
	cout << " ------------------------------------" << endl;

	// Get observation
	// TH1D *H_nsigobs = new TH1D("Nsigobs", "Observed N_sig in Pt1 bins", nptbins, ptbins);
	// NObs(H_nsigobs, inputsamples, &MuonPlotter::isGenMatchedSUSYDiLepEvent);

	// cout << "  Observation from LM0:      " << H_nsigobs->GetBinContent(1) << " +/- " << H_nsigobs->GetBinError(1) << endl;
	cout << "  Nt2 observed (TTbar only):  " << fLumiNorm/fSamples[TTbar].lumi *fSamples[TTbar].mumu.numbers.nt2 << endl;
	cout << "  Nt2 observed (VVjets only): " << fLumiNorm/fSamples[VVJets].lumi*fSamples[VVJets].mumu.numbers.nt2 << endl;
	cout << "  Nt2 observed (Wjets only):  " << fLumiNorm/fSamples[WJets].lumi *fSamples[WJets].mumu.numbers.nt2 << endl;
	cout << "  Nt2 observed (Zjets only):  " << fLumiNorm/fSamples[ZJets].lumi *fSamples[ZJets].mumu.numbers.nt2 << endl;
	float nt2qcd = fLumiNorm/fSamples[QCD15].lumi*fSamples[QCD15].mumu.numbers.nt2;
	nt2qcd += fLumiNorm/fSamples[QCD30].lumi*fSamples[QCD30].mumu.numbers.nt2;
	nt2qcd += fLumiNorm/fSamples[QCD80].lumi*fSamples[QCD80].mumu.numbers.nt2;
	nt2qcd += fLumiNorm/fSamples[QCD170].lumi*fSamples[QCD170].mumu.numbers.nt2;
	cout << "  Nt2 observed (QCD only):    " << nt2qcd << endl;
	cout << " ------------------------------------" << endl;
}

//____________________________________________________________________________
void MuonPlotter::makeIntPredictionEE(){
	cout << "-----------------------------------" << endl;
	cout << " Producing prediction for EE channel" << endl << endl;
	bool data = true; // Use ratios from data or mc?

	// Which samples to use for nt2/nt1/nt0 input?
	// vector<int> inputsamples = fMCBG;
	// vector<int> inputsamples = fMCBGSig;
	vector<int> inputsamples = fEGData;

	// Which luminosity to use?
	// fLumiNorm = fSamples[MuA].lumi + fSamples[MuB].lumi;
	fLumiNorm = 22.;

	// Dummy binning
	const int nptbins = 1;
	const double ptbins[nptbins+1] = {5., 1000.};
	const int netabins = 1;
	const double etabins[netabins+1] = {-2.5, 2.5};

	// Fill the ratios
	TH2D *H_el_fratio_allmc = new TH2D("ElfRatioAllMC", "fRatio for all MC", nptbins, ptbins, netabins, etabins);
	TH2D *H_el_fratio_data  = new TH2D("ElfRatioData",  "fRatio for Mu Data", nptbins, ptbins, netabins, etabins);
	TH2D *H_el_pratio_allmc = new TH2D("ElpRatioAllMC", "pRatio for all MC", nptbins, ptbins, netabins, etabins);
	TH2D *H_el_pratio_data  = new TH2D("ElpRatioData",  "pRatio for Mu Data", nptbins, ptbins, netabins, etabins);

	// TH2D *H_fratio_ttbar = fillRatio(TTbar, 0, &MuonPlotter::isGoodMuEvent, &MuonPlotter::isFakeTTbarMuon, nptbins, ptbins, netabins, etabins);
	// TH2D *H_pratio_ttbar = fillRatio(TTbar, 0, &MuonPlotter::isGoodMuEvent, &MuonPlotter::isPromptTTbarMuon, nptbins, ptbins, netabins, etabins);
	// H_fratio_ttbar->SetName("fRatioTTbar");
	// H_pratio_ttbar->SetName("pRatioTTbar");

	float elfratio_allmc(0.), elfratio_allmc_e(0.);
	float elfratio_data(0.),  elfratio_data_e(0.);
	float elpratio_allmc(0.), elpratio_allmc_e(0.);
	float elpratio_data(0.),  elpratio_data_e(0.);

	// float fratio_ttbar   = H_fratio_ttbar->GetBinContent(1,1);
	// float fratio_ttbar_e = H_fratio_ttbar->GetBinError(1,1);

	// float pratio_ttbar   = H_pratio_ttbar->GetBinContent(1,1);
	// float pratio_ttbar_e = H_pratio_ttbar->GetBinError(1,1);

	calculateRatio(fMCBG,   Electron, 1, elfratio_allmc, elfratio_allmc_e);
	calculateRatio(fMuData, Electron, 1, elfratio_data,  elfratio_data_e);
	calculateRatio(fMCBG,   Electron, 2, elpratio_allmc, elpratio_allmc_e);
	calculateRatio(fMuData, Electron, 2, elpratio_data,  elpratio_data_e);

	H_el_fratio_allmc->SetBinContent(1,1,elfratio_allmc);
	H_el_fratio_allmc->SetBinError  (1,1,elfratio_allmc_e);
	H_el_fratio_data ->SetBinContent(1,1,elfratio_data);
	H_el_fratio_data ->SetBinError  (1,1,elfratio_data_e);

	H_el_pratio_allmc->SetBinContent(1,1,elpratio_allmc);
	H_el_pratio_allmc->SetBinError  (1,1,elpratio_allmc_e);
	H_el_pratio_data ->SetBinContent(1,1,elpratio_data);
	H_el_pratio_data ->SetBinError  (1,1,elpratio_data_e);

	cout << "  Electron fRatio from all MC     = " << elfratio_allmc << " +/- " << elfratio_allmc_e << endl;
	// cout << "  fRatio from ttbar genmatch = " << elfratio_ttbar << " +/- " << elfratio_ttbar_e << endl;
	cout << "  Electron fRatio from data       = " << elfratio_data  << " +/- " << elfratio_data_e << endl;
	cout << " ------------------------------------" << endl;
	cout << "  Electron pRatio from all MC     = " << elpratio_allmc << " +/- " << elpratio_allmc_e << endl;
	// cout << "  pRatio from ttbar genmatch = " << elpratio_ttbar << " +/- " << elspratio_ttbar_e << endl;
	cout << "  Electron pRatio from data       = " << elpratio_data  << " +/- " << elpratio_data_e << endl;
	cout << " ------------------------------------" << endl;

	// // Add systematics to ratios
	// float deviation(0.), olderror(0.), newerror(0.);
	// deviation = fabs(fratio_allmc - fratio_ttbar);
	// // Add to mc ratios
	// olderror = H_fratio_allmc->GetBinError(1,1);
	// newerror = sqrt(olderror*olderror + deviation*deviation);
	// H_fratio_allmc->SetBinError(1,1,newerror);
	// fratio_allmc_e = newerror;
	// // Add to data ratios
	// olderror = H_fratio_data->GetBinError(1,1);
	// newerror = sqrt(olderror*olderror + deviation*deviation);
	// H_fratio_data->SetBinError(1,1,newerror);
	// fratio_data_e = newerror;
	// cout << "  new fRatio (all MC)        = " << H_fratio_allmc->GetBinContent(1,1) << " +/- " << H_fratio_allmc->GetBinError(1,1) << endl;
	// cout << "  new fRatio (data)          = " << H_fratio_data ->GetBinContent(1,1) << " +/- " << H_fratio_data ->GetBinError(1,1) << endl;
	// 
	// deviation  = fabs(pratio_allmc - pratio_ttbar);
	// // Add to mc ratios
	// olderror = pratio_allmc_e;
	// newerror = sqrt(olderror*olderror + deviation*deviation);
	// H_pratio_allmc->SetBinError(1,1,newerror);
	// pratio_allmc_e = newerror;
	// // Add to data ratios
	// olderror = pratio_data_e;
	// newerror = sqrt(olderror*olderror + deviation*deviation);
	// H_pratio_data->SetBinError(1,1,newerror);
	// pratio_data_e = newerror;
	// cout << "  new pRatio (all MC)        = " << H_pratio_allmc->GetBinContent(1,1) << " +/- " << H_pratio_allmc->GetBinError(1,1) << endl;
	// cout << "  new pRatio (data)          = " << H_pratio_data ->GetBinContent(1,1) << " +/- " << H_pratio_data ->GetBinError(1,1) << endl;
	// cout << " ------------------------------------" << endl;

	double nt2(0.),    nt1(0.),    nt0(0.);
	double nt2_e2(0.), nt1_e2(0.), nt0_e2(0.);
	for(size_t i = 0; i < inputsamples.size(); ++i){
		int index = inputsamples[i];
		float scale = fLumiNorm/fSamples[index].lumi; // Normalize all
		channel *cha = &fSamples[index].ee;
		if(data) scale = 1.;
		nt2 += scale * cha->numbers.nt2;
		nt1 += scale * cha->numbers.nt10;
		nt0 += scale * cha->numbers.nt0;
		nt2_e2 += scale*scale * cha->numbers.nt2;
		nt1_e2 += scale*scale * cha->numbers.nt10;
		nt0_e2 += scale*scale * cha->numbers.nt0;
	}
	if(data){
		cout << "  Found " << nt2 << " events in Nt2" << endl;
		cout << "  Found " << nt1 << " events in Nt1" << endl;
		cout << "  Found " << nt0 << " events in Nt0" << endl;
	}
	if(!data){
		cout << "  Found " << nt2 << " +/- " << TMath::Sqrt(nt2_e2) << " events in Nt2" << endl;
		cout << "  Found " << nt1 << " +/- " << TMath::Sqrt(nt1_e2) << " events in Nt1" << endl;
		cout << "  Found " << nt0 << " +/- " << TMath::Sqrt(nt0_e2) << " events in Nt0" << endl;		
	}
	cout << " ------------------------------------" << endl;

	// Make prediction
	fFPRatios = new FPRatios();
	fFPRatios->SetVerbose(fVerbose);
	if(data){
		fFPRatios->SetElFratios(H_el_fratio_data);
		fFPRatios->SetElPratios(H_el_pratio_data);
	}
	else{
		fFPRatios->SetElFratios(H_el_fratio_allmc);
		fFPRatios->SetElPratios(H_el_pratio_allmc);
	}
	vector<double> nt;
	nt.push_back(nt0);
	nt.push_back(nt1);
	nt.push_back(nt2);
	
	fFPRatios->NevtTopol(2, 0, nt);

	vector<double> vpt, veta;
	vpt.push_back(30.); vpt.push_back(30.); // Fake pts and etas (first electron then muon)
	veta.push_back(0.); veta.push_back(0.);

	vector<double> nevFP = fFPRatios->NevtPass(vpt, veta);
	vector<double> nevFPEstat = fFPRatios->NevtPassErrStat();
	vector<double> nevFPEsyst = fFPRatios->NevtPassErrSyst();

	cout << "  Prediction for Npp: " << nevFP[0] << " +/- " << nevFPEstat[0] << " (stat) +/- " << nevFPEsyst[0] << " (syst)" << endl;
	cout << "  Prediction for Nfp: " << nevFP[1] << " +/- " << nevFPEstat[1] << " (stat) +/- " << nevFPEsyst[1] << " (syst)" << endl;
	cout << "  Prediction for Nff: " << nevFP[2] << " +/- " << nevFPEstat[2] << " (stat) +/- " << nevFPEsyst[2] << " (syst)" << endl;
	cout << "  Total fakes:        " << nevFP[1]+nevFP[2] << " +/- " << TMath::Sqrt(nevFPEstat[1]*nevFPEstat[1] + nevFPEstat[2]*nevFPEstat[2]) << " (stat) +/- " << TMath::Sqrt(nevFPEsyst[1]*nevFPEsyst[1] + nevFPEsyst[2]*nevFPEsyst[2])<< " (syst)" << endl;
	cout << " ------------------------------------" << endl;

	// Get observation
	// TH1D *H_nsigobs = new TH1D("Nsigobs", "Observed N_sig in Pt1 bins", nptbins, ptbins);
	// NObs(H_nsigobs, inputsamples, &MuonPlotter::isGenMatchedSUSYDiLepEvent);

	// cout << "  Observation from LM0:      " << H_nsigobs->GetBinContent(1) << " +/- " << H_nsigobs->GetBinError(1) << endl;
	cout << "  Nt2 observed (TTbar only):  " << fLumiNorm/fSamples[TTbar].lumi *fSamples[TTbar].ee.numbers.nt2 << endl;
	cout << "  Nt2 observed (VVjets only): " << fLumiNorm/fSamples[VVJets].lumi*fSamples[VVJets].ee.numbers.nt2 << endl;
	cout << "  Nt2 observed (Wjets only):  " << fLumiNorm/fSamples[WJets].lumi *fSamples[WJets].ee.numbers.nt2 << endl;
	cout << "  Nt2 observed (Zjets only):  " << fLumiNorm/fSamples[ZJets].lumi *fSamples[ZJets].ee.numbers.nt2 << endl;
	float nt2qcd = fLumiNorm/fSamples[QCD15].lumi*fSamples[QCD15].ee.numbers.nt2;
	nt2qcd += fLumiNorm/fSamples[QCD30].lumi*fSamples[QCD30].ee.numbers.nt2;
	nt2qcd += fLumiNorm/fSamples[QCD80].lumi*fSamples[QCD80].ee.numbers.nt2;
	nt2qcd += fLumiNorm/fSamples[QCD170].lumi*fSamples[QCD170].ee.numbers.nt2;
	cout << "  Nt2 observed (QCD only):    " << nt2qcd << endl;
	cout << " ------------------------------------" << endl;
}

//____________________________________________________________________________
void MuonPlotter::makeIntPredictionEMu(){
	cout << "-----------------------------------" << endl;
	cout << " Producing prediction for EMu channel" << endl << endl;
	bool data = true; // Use ratios from data or mc?

	// Which samples to use for nt2/nt1/nt0 input?
	// vector<int> inputsamples = fMCBG;
	// vector<int> inputsamples = fMCBGSig;
	vector<int> inputsamples = fMuData;

	// Which luminosity to use?
	// fLumiNorm = fSamples[MuA].lumi + fSamples[MuB].lumi;
	fLumiNorm = 22.;

	// Dummy binning
	const int nptbins = 1;
	const double ptbins[nptbins+1] = {5., 1000.};
	const int netabins = 1;
	const double etabins[netabins+1] = {-2.5, 2.5};

	// Fill the ratios
	TH2D *H_mu_fratio_allmc = new TH2D("MufRatioAllMC", "fRatio for all MC", nptbins, ptbins, netabins, etabins);
	TH2D *H_mu_fratio_data  = new TH2D("MufRatioData",  "fRatio for Mu Data", nptbins, ptbins, netabins, etabins);
	TH2D *H_mu_pratio_allmc = new TH2D("MupRatioAllMC", "pRatio for all MC", nptbins, ptbins, netabins, etabins);
	TH2D *H_mu_pratio_data  = new TH2D("MupRatioData",  "pRatio for Mu Data", nptbins, ptbins, netabins, etabins);
	TH2D *H_el_fratio_allmc = new TH2D("ElfRatioAllMC", "fRatio for all MC", nptbins, ptbins, netabins, etabins);
	TH2D *H_el_fratio_data  = new TH2D("ElfRatioData",  "fRatio for Mu Data", nptbins, ptbins, netabins, etabins);
	TH2D *H_el_pratio_allmc = new TH2D("ElpRatioAllMC", "pRatio for all MC", nptbins, ptbins, netabins, etabins);
	TH2D *H_el_pratio_data  = new TH2D("ElpRatioData",  "pRatio for Mu Data", nptbins, ptbins, netabins, etabins);

	// TH2D *H_fratio_ttbar = fillRatio(TTbar, 0, &MuonPlotter::isGoodMuEvent, &MuonPlotter::isFakeTTbarMuon, nptbins, ptbins, netabins, etabins);
	// TH2D *H_pratio_ttbar = fillRatio(TTbar, 0, &MuonPlotter::isGoodMuEvent, &MuonPlotter::isPromptTTbarMuon, nptbins, ptbins, netabins, etabins);
	// H_fratio_ttbar->SetName("fRatioTTbar");
	// H_pratio_ttbar->SetName("pRatioTTbar");

	float mufratio_allmc(0.), mufratio_allmc_e(0.);
	float mufratio_data(0.),  mufratio_data_e(0.);
	float mupratio_allmc(0.), mupratio_allmc_e(0.);
	float mupratio_data(0.),  mupratio_data_e(0.);
	float elfratio_allmc(0.), elfratio_allmc_e(0.);
	float elfratio_data(0.),  elfratio_data_e(0.);
	float elpratio_allmc(0.), elpratio_allmc_e(0.);
	float elpratio_data(0.),  elpratio_data_e(0.);

	// float fratio_ttbar   = H_fratio_ttbar->GetBinContent(1,1);
	// float fratio_ttbar_e = H_fratio_ttbar->GetBinError(1,1);

	// float pratio_ttbar   = H_pratio_ttbar->GetBinContent(1,1);
	// float pratio_ttbar_e = H_pratio_ttbar->GetBinError(1,1);

	calculateRatio(fMCBG,   Muon, 1, mufratio_allmc, mufratio_allmc_e);
	calculateRatio(fMuData, Muon, 1, mufratio_data,  mufratio_data_e);
	calculateRatio(fMCBG,   Muon, 2, mupratio_allmc, mupratio_allmc_e);
	calculateRatio(fMuData, Muon, 2, mupratio_data,  mupratio_data_e);

	calculateRatio(fMCBG,   Electron, 1, elfratio_allmc, elfratio_allmc_e);
	calculateRatio(fMuData, Electron, 1, elfratio_data,  elfratio_data_e);
	calculateRatio(fMCBG,   Electron, 2, elpratio_allmc, elpratio_allmc_e);
	calculateRatio(fMuData, Electron, 2, elpratio_data,  elpratio_data_e);

	H_mu_fratio_allmc->SetBinContent(1,1,mufratio_allmc);
	H_mu_fratio_allmc->SetBinError  (1,1,mufratio_allmc_e);
	H_mu_fratio_data ->SetBinContent(1,1,mufratio_data);
	H_mu_fratio_data ->SetBinError  (1,1,mufratio_data_e);
	H_el_fratio_allmc->SetBinContent(1,1,elfratio_allmc);
	H_el_fratio_allmc->SetBinError  (1,1,elfratio_allmc_e);
	H_el_fratio_data ->SetBinContent(1,1,elfratio_data);
	H_el_fratio_data ->SetBinError  (1,1,elfratio_data_e);

	H_mu_pratio_allmc->SetBinContent(1,1,mupratio_allmc);
	H_mu_pratio_allmc->SetBinError  (1,1,mupratio_allmc_e);
	H_mu_pratio_data ->SetBinContent(1,1,mupratio_data);
	H_mu_pratio_data ->SetBinError  (1,1,mupratio_data_e);
	H_el_pratio_allmc->SetBinContent(1,1,elpratio_allmc);
	H_el_pratio_allmc->SetBinError  (1,1,elpratio_allmc_e);
	H_el_pratio_data ->SetBinContent(1,1,elpratio_data);
	H_el_pratio_data ->SetBinError  (1,1,elpratio_data_e);

	cout << "  Muon fRatio from all MC         = " << mufratio_allmc << " +/- " << mufratio_allmc_e << endl;
	// cout << "  fRatio from ttbar genmatch = " << mufratio_ttbar << " +/- " << mufratio_ttbar_e << endl;
	cout << "  Muon fRatio from data           = " << mufratio_data  << " +/- " << mufratio_data_e << endl;
	cout << " ------------------------------------" << endl;
	cout << "  Electron fRatio from all MC     = " << elfratio_allmc << " +/- " << elfratio_allmc_e << endl;
	// cout << "  fRatio from ttbar genmatch = " << elfratio_ttbar << " +/- " << elfratio_ttbar_e << endl;
	cout << "  Electron fRatio from data       = " << elfratio_data  << " +/- " << elfratio_data_e << endl;
	cout << " ------------------------------------" << endl;
	cout << "  Muon pRatio from all MC         = " << mupratio_allmc << " +/- " << mupratio_allmc_e << endl;
	// cout << "  pRatio from ttbar genmatch = " << mupratio_ttbar << " +/- " << mupratio_ttbar_e << endl;
	cout << "  Muon pRatio from data           = " << mupratio_data  << " +/- " << mupratio_data_e << endl;
	cout << " ------------------------------------" << endl;
	cout << "  Electron pRatio from all MC     = " << elpratio_allmc << " +/- " << elpratio_allmc_e << endl;
	// cout << "  pRatio from ttbar genmatch = " << elpratio_ttbar << " +/- " << elspratio_ttbar_e << endl;
	cout << "  Electron pRatio from data       = " << elpratio_data  << " +/- " << elpratio_data_e << endl;
	cout << " ------------------------------------" << endl;

	// // Add systematics to ratios
	// float deviation(0.), olderror(0.), newerror(0.);
	// deviation = fabs(fratio_allmc - fratio_ttbar);
	// // Add to mc ratios
	// olderror = H_fratio_allmc->GetBinError(1,1);
	// newerror = sqrt(olderror*olderror + deviation*deviation);
	// H_fratio_allmc->SetBinError(1,1,newerror);
	// fratio_allmc_e = newerror;
	// // Add to data ratios
	// olderror = H_fratio_data->GetBinError(1,1);
	// newerror = sqrt(olderror*olderror + deviation*deviation);
	// H_fratio_data->SetBinError(1,1,newerror);
	// fratio_data_e = newerror;
	// cout << "  new fRatio (all MC)        = " << H_fratio_allmc->GetBinContent(1,1) << " +/- " << H_fratio_allmc->GetBinError(1,1) << endl;
	// cout << "  new fRatio (data)          = " << H_fratio_data ->GetBinContent(1,1) << " +/- " << H_fratio_data ->GetBinError(1,1) << endl;
	// 
	// deviation  = fabs(pratio_allmc - pratio_ttbar);
	// // Add to mc ratios
	// olderror = pratio_allmc_e;
	// newerror = sqrt(olderror*olderror + deviation*deviation);
	// H_pratio_allmc->SetBinError(1,1,newerror);
	// pratio_allmc_e = newerror;
	// // Add to data ratios
	// olderror = pratio_data_e;
	// newerror = sqrt(olderror*olderror + deviation*deviation);
	// H_pratio_data->SetBinError(1,1,newerror);
	// pratio_data_e = newerror;
	// cout << "  new pRatio (all MC)        = " << H_pratio_allmc->GetBinContent(1,1) << " +/- " << H_pratio_allmc->GetBinError(1,1) << endl;
	// cout << "  new pRatio (data)          = " << H_pratio_data ->GetBinContent(1,1) << " +/- " << H_pratio_data ->GetBinError(1,1) << endl;
	// cout << " ------------------------------------" << endl;

	double nt2(0.),    nt10(0.),    nt01(0.),    nt0(0.);
	double nt2_e2(0.), nt10_e2(0.), nt01_e2(0.), nt0_e2(0.);
	for(size_t i = 0; i < inputsamples.size(); ++i){
		int index = inputsamples[i];
		float scale = fLumiNorm/fSamples[index].lumi; // Normalize all
		channel *cha = &fSamples[index].emu;
		if(data) scale = 1.;
		nt2  += scale * cha->numbers.nt2;
		nt10 += scale * cha->numbers.nt10;
		nt01 += scale * cha->numbers.nt01;
		nt0  += scale * cha->numbers.nt0;
		nt2_e2  += scale*scale * cha->numbers.nt2;
		nt10_e2 += scale*scale * cha->numbers.nt10;
		nt01_e2 += scale*scale * cha->numbers.nt01;
		nt0_e2  += scale*scale * cha->numbers.nt0;
	}
	if(data){
		cout << "  Found " << nt2  << " events in Nt2" << endl;
		cout << "  Found " << nt10 << " events in Nt10" << endl;
		cout << "  Found " << nt01 << " events in Nt01" << endl;
		cout << "  Found " << nt0  << " events in Nt0" << endl;
	}
	if(!data){
		cout << "  Found " << nt2  << " +/- " << TMath::Sqrt(nt2_e2)  << " events in Nt2" << endl;
		cout << "  Found " << nt10 << " +/- " << TMath::Sqrt(nt10_e2) << " events in Nt10" << endl;
		cout << "  Found " << nt01 << " +/- " << TMath::Sqrt(nt01_e2) << " events in Nt01" << endl;
		cout << "  Found " << nt0  << " +/- " << TMath::Sqrt(nt0_e2)  << " events in Nt0" << endl;		
	}
	cout << " ------------------------------------" << endl;

	// Make prediction
	fFPRatios = new FPRatios();
	fFPRatios->SetVerbose(fVerbose);
	if(data){
		fFPRatios->SetMuFratios(H_mu_fratio_data);
		fFPRatios->SetMuPratios(H_mu_pratio_data);
		fFPRatios->SetElFratios(H_el_fratio_data);
		fFPRatios->SetElPratios(H_el_pratio_data);
	}
	else{
		fFPRatios->SetMuFratios(H_mu_fratio_allmc);
		fFPRatios->SetMuPratios(H_mu_pratio_allmc);
		fFPRatios->SetElFratios(H_el_fratio_allmc);
		fFPRatios->SetElPratios(H_el_pratio_allmc);
	}
	vector<double> nt;
	nt.push_back(nt0);
	nt.push_back(nt01); // e passes
	nt.push_back(nt10); // mu passes
	nt.push_back(nt2);
	
	fFPRatios->NevtTopol(1, 1, nt);

	vector<double> vpt, veta;
	vpt.push_back(30.); vpt.push_back(30.); // Fake pts and etas (first electron then muon)
	veta.push_back(0.); veta.push_back(0.);

	vector<double> nevFP = fFPRatios->NevtPass(vpt, veta);
	vector<double> nevFPEstat = fFPRatios->NevtPassErrStat();
	vector<double> nevFPEsyst = fFPRatios->NevtPassErrSyst();

	cout << "  Prediction for Npp: " << nevFP[0] << " +/- " << nevFPEstat[0] << " (stat) +/- " << nevFPEsyst[0] << " (syst)" << endl;
	cout << "  Prediction for Nfp: " << nevFP[1] << " +/- " << nevFPEstat[1] << " (stat) +/- " << nevFPEsyst[1] << " (syst)" << endl;
	cout << "  Prediction for Nff: " << nevFP[2] << " +/- " << nevFPEstat[2] << " (stat) +/- " << nevFPEsyst[2] << " (syst)" << endl;
	cout << "  Total fakes:        " << nevFP[1]+nevFP[2] << " +/- " << TMath::Sqrt(nevFPEstat[1]*nevFPEstat[1] + nevFPEstat[2]*nevFPEstat[2]) << " (stat) +/- " << TMath::Sqrt(nevFPEsyst[1]*nevFPEsyst[1] + nevFPEsyst[2]*nevFPEsyst[2])<< " (syst)" << endl;
	cout << " ------------------------------------" << endl;

	// Get observation
	// TH1D *H_nsigobs = new TH1D("Nsigobs", "Observed N_sig in Pt1 bins", nptbins, ptbins);
	// NObs(H_nsigobs, inputsamples, &MuonPlotter::isGenMatchedSUSYDiLepEvent);

	// cout << "  Observation from LM0:      " << H_nsigobs->GetBinContent(1) << " +/- " << H_nsigobs->GetBinError(1) << endl;
	cout << "  Nt2 observed (TTbar only):  " << fLumiNorm/fSamples[TTbar].lumi *fSamples[TTbar].emu.numbers.nt2 << endl;
	cout << "  Nt2 observed (VVjets only): " << fLumiNorm/fSamples[VVJets].lumi*fSamples[VVJets].emu.numbers.nt2 << endl;
	cout << "  Nt2 observed (Wjets only):  " << fLumiNorm/fSamples[WJets].lumi *fSamples[WJets].emu.numbers.nt2 << endl;
	cout << "  Nt2 observed (Zjets only):  " << fLumiNorm/fSamples[ZJets].lumi *fSamples[ZJets].emu.numbers.nt2 << endl;
	float nt2qcd = fLumiNorm/fSamples[QCD15].lumi*fSamples[QCD15].emu.numbers.nt2;
	nt2qcd += fLumiNorm/fSamples[QCD30].lumi*fSamples[QCD30].emu.numbers.nt2;
	nt2qcd += fLumiNorm/fSamples[QCD80].lumi*fSamples[QCD80].emu.numbers.nt2;
	nt2qcd += fLumiNorm/fSamples[QCD170].lumi*fSamples[QCD170].emu.numbers.nt2;
	cout << "  Nt2 observed (QCD only):    " << nt2qcd << endl;
	cout << " ------------------------------------" << endl;
}

//____________________________________________________________________________
vector<TH1D*> MuonPlotter::NsigPredFromFPRatios(const int sample, bool output){
	///////////////////////////////////////////////////////////////////////////
	// Note: Careful, this is only a workaround at the moment, and is only
	//       really valid for one single eta bin!
	//       In the future we should rewrite the interface to FPRatios and
	//       pass it the ratios directly to have full control
	///////////////////////////////////////////////////////////////////////////
	if(fVerbose > 2) cout << "MuonPlotter::NsigPredFromFPRatios ==> Predicting Nsig from " << fSamples[sample].sname << endl;
	vector<TH1D*> res;

	TH1D *H_dummy    = new TH1D("Dummyhist", "Dummy", gNPt2bins, gPt2bins);
	TH1D *H_nsigpred = new TH1D(Form("Nsigpred_%s", fSamples[sample].sname.Data()), "Predicted N_sig in Pt1 bins", gNPtbins, gPtbins);
	TH1D *H_nfppred  = new TH1D(Form("Nfppred_%s",  fSamples[sample].sname.Data()), "Predicted N_fp in Pt1 bins",  gNPtbins, gPtbins);
	TH1D *H_nffpred  = new TH1D(Form("Nffpred_%s",  fSamples[sample].sname.Data()), "Predicted N_ff in Pt1 bins",  gNPtbins, gPtbins);	
	TH2D *H_nt2mes   = new TH2D(Form("Nt2mes_%s",   fSamples[sample].sname.Data()), "Measured N_t2 in Pt1/2 bins", gNPt2bins, gPt2bins, gNPt2bins, gPt2bins);
	TH2D *H_nt1mes   = new TH2D(Form("Nt1mes_%s",   fSamples[sample].sname.Data()), "Measured N_t1 in Pt1/2 bins", gNPt2bins, gPt2bins, gNPt2bins, gPt2bins);
	TH2D *H_nt0mes   = new TH2D(Form("Nt0mes_%s",   fSamples[sample].sname.Data()), "Measured N_t0 in Pt1/2 bins", gNPt2bins, gPt2bins, gNPt2bins, gPt2bins);

	// // Fill histograms from tree
	// TTree *tree = fSamples[sample].tree;
	// Init(tree);
	// if (fChain == 0) return res;
	// Long64_t nentries = fChain->GetEntriesFast();
	// Long64_t nbytes = 0, nb = 0;
	// for (Long64_t jentry=0; jentry<nentries;jentry++) {
	// 	Long64_t ientry = LoadTree(jentry);
	// 	if (ientry < 0) break;
	// 	nb = fChain->GetEntry(jentry);   nbytes += nb;
	// 
	// 	if(isSSLLEvent() == false) continue;
	// 
	// 	if(  isTightMuon(0) &&  isTightMuon(1) ) H_nt2mes->Fill(MuPt[0], MuPt[1]); // Tight-tight
	// 	if(  isTightMuon(0) && !isTightMuon(1) ) H_nt1mes->Fill(MuPt[0], MuPt[1]); // Tight-loose
	// 	if( !isTightMuon(0) &&  isTightMuon(1) ) H_nt1mes->Fill(MuPt[1], MuPt[0]); // Loose-tight
	// 	if( !isTightMuon(0) && !isTightMuon(1) ) H_nt0mes->Fill(MuPt[0], MuPt[1]); // Loose-loose
	// }

	channel cha = fSamples[sample].mumu;
	H_nt2mes = cha.nthistos.h_nt2;
	H_nt1mes = cha.nthistos.h_nt10;
	H_nt0mes = cha.nthistos.h_nt0;

	if(fVerbose > 2){
		cout << " Found " << H_nt2mes->GetEntries() << " events with tight-tight (Nt2)" << endl;
		cout << " Found " << H_nt1mes->GetEntries() << " events with tight-loose (Nt1)" << endl;
		cout << " Found " << H_nt0mes->GetEntries() << " events with loose-loose (Nt0)" << endl;
	}

	for(size_t i = 1; i <= gNPt2bins; ++i){
		double pt1 = H_dummy->GetBinCenter(i);
		if(fVerbose > 2){
			cout << "=======================================" << endl;
			cout << "Pt1 = " << pt1 << endl;
		}
		double eta1 = 0.0;
		double npppred(0.0), npppredEstat2(0.0), npppredEsyst2(0.0);
		double nfppred(0.0), nfppredEstat2(0.0), nfppredEsyst2(0.0);
		double nffpred(0.0), nffpredEstat2(0.0), nffpredEsyst2(0.0);
		for(size_t j = 1; j <= gNPt2bins; ++j){
			if(fVerbose > 2) cout << " --------" << endl;
			double pt2 = H_dummy->GetBinCenter(j);
			// double pt2 = H_nsigpred->GetBinCenter(j);
			double eta2 = 0.0;
			double nt2 = H_nt2mes->GetBinContent(i,j);
			double nt1 = H_nt1mes->GetBinContent(i,j);
			double nt0 = H_nt0mes->GetBinContent(i,j);

			if(fVerbose > 2) cout << "   Pt2 = " << pt2 << endl;
			if(fVerbose > 2) cout << "   nt2: " << nt2 << "  nt1: " << nt1 << "  nt0: " << nt0 << endl;

			fFPRatios = new FPRatios();
			fFPRatios->SetVerbose(fVerbose);
			fFPRatios->SetMuFratios(fH2D_fRatio);
			fFPRatios->SetMuPratios(fH2D_pRatio);
			vector<double> nt;
			nt.push_back(nt0); nt.push_back(nt1); nt.push_back(nt2);
			fFPRatios->NevtTopol(0, 2, nt);

			vector<double> vpt, veta;
			vpt.push_back(pt1); vpt.push_back(pt2);
			veta.push_back(eta1); veta.push_back(eta2);

			vector<double> nevFP = fFPRatios->NevtPass(vpt, veta);
			vector<double> nevFPEstat = fFPRatios->NevtPassErrStat();
			vector<double> nevFPEsyst = fFPRatios->NevtPassErrSyst();
			npppred += nevFP[0];
			nfppred += nevFP[1];
			nffpred += nevFP[2];
			npppredEstat2 += nevFPEstat[0]*nevFPEstat[0];
			npppredEsyst2 += nevFPEsyst[0]*nevFPEsyst[0];
			nfppredEstat2 += nevFPEstat[1]*nevFPEstat[1];
			nfppredEsyst2 += nevFPEsyst[1]*nevFPEsyst[1];
			nffpredEstat2 += nevFPEstat[2]*nevFPEstat[2];
			nffpredEsyst2 += nevFPEsyst[2]*nevFPEsyst[2];
			delete fFPRatios;
		}
		Int_t bin = H_nsigpred->FindBin(pt1);
		H_nsigpred->SetBinContent(bin, npppred);
		H_nfppred ->SetBinContent(bin, nfppred);
		H_nffpred ->SetBinContent(bin, nffpred);
		H_nsigpred->SetBinError(bin, sqrt(npppredEstat2 + npppredEsyst2));
		H_nfppred ->SetBinError(bin, sqrt(nfppredEstat2 + nfppredEsyst2));
		H_nffpred ->SetBinError(bin, sqrt(nffpredEstat2 + nffpredEsyst2));
	}	

	if(fVerbose > 2) cout << " Predict " << H_nsigpred->Integral() << " signal events (Nsig = p^2*Npp) from this sample" << endl;
	if(fVerbose > 2) cout << " Predict " << H_nfppred->Integral() << " fake-prompt events (f*p*Nfp) from this sample" << endl;
	if(fVerbose > 2) cout << " Predict " << H_nffpred->Integral() << " fake-fake events (f*f*Nff) from this sample" << endl;
	// Output
	H_nsigpred->SetXTitle(convertVarName("MuPt[1]"));
	if(output){
		H_nt2mes->SetXTitle(convertVarName("MuPt[1]"));
		H_nt2mes->SetYTitle(convertVarName("MuPt[1]"));
		H_nt1mes->SetXTitle(convertVarName("MuPt[1]"));
		H_nt1mes->SetYTitle(convertVarName("MuPt[1]"));
		H_nt0mes->SetXTitle(convertVarName("MuPt[1]"));
		H_nt0mes->SetYTitle(convertVarName("MuPt[1]"));
		// if(H_nsigpred->GetMinimum() < 0) H_nsigpred->SetMaximum(0);
		// if(H_nsigpred->GetMinimum() > 0) H_nsigpred->SetMinimum(0);
		printHisto(H_nsigpred, H_nsigpred->GetName(), H_nsigpred->GetTitle());
		printHisto(H_nt2mes, H_nt2mes->GetName(), H_nt2mes->GetTitle(), "colz");
		printHisto(H_nt1mes, H_nt1mes->GetName(), H_nt1mes->GetTitle(), "colz");
		printHisto(H_nt0mes, H_nt0mes->GetName(), H_nt0mes->GetTitle(), "colz");
	}
	res.push_back(H_nsigpred);
	res.push_back(H_nfppred);
	res.push_back(H_nffpred);
	return res;
}

//____________________________________________________________________________
void MuonPlotter::fillYields(){ fillYields(fAllSamples); }
void MuonPlotter::fillYields(int sample){ vector<int> samples; samples.push_back(sample); fillYields(samples); }
void MuonPlotter::fillYields(vector<int> samples){
	for(size_t i = 0; i < samples.size(); ++i){
		int index = samples[i];
		if(fVerbose > 1){
			cout << "-------------------" << endl;
			cout << " Filling yields for " << fSamples[index].name << endl;
		}

		TTree *tree = fSamples[index].tree;
		const bool isdata = fSamples[index].isdata;

		channel mumu = fSamples[index].mumu;
		channel ee   = fSamples[index].ee;
		channel emu  = fSamples[index].emu;
		
		tree->ResetBranchAddresses();
		Init(tree);
		if (fChain == 0) return;
		Long64_t nentries = fChain->GetEntriesFast();
		Long64_t nbytes = 0, nb = 0;
		for (Long64_t jentry=0; jentry<nentries;jentry++) {
			Long64_t ientry = LoadTree(jentry);
			if (ientry < 0) break;
			nb = fChain->GetEntry(jentry);   nbytes += nb;

			// MuMu Channel
			if(isSSLLMuEventTRG()){
				if(  isTightMuon(0) &&  isTightMuon(1) ){ // Tight-tight
					mumu.nthistos.h_nt2    ->Fill(MuPt[0], MuPt[1]);
					mumu.nthistos.h_nt2_pt ->Fill(MuPt[0]);
					mumu.nthistos.h_nt2_eta->Fill(MuEta[0]);
				}
				if(  isTightMuon(0) && !isTightMuon(1) ){ // Tight-loose
					mumu.nthistos.h_nt10    ->Fill(MuPt[0], MuPt[1]);
					mumu.nthistos.h_nt10_pt ->Fill(MuPt[0]);
					mumu.nthistos.h_nt10_eta->Fill(MuEta[0]);
				}
				if( !isTightMuon(0) &&  isTightMuon(1) ){ // Loose-tight
					mumu.nthistos.h_nt10    ->Fill(MuPt[1], MuPt[0]);
					mumu.nthistos.h_nt10_pt ->Fill(MuPt[1]);
					mumu.nthistos.h_nt10_eta->Fill(MuEta[1]);
				}
				if( !isTightMuon(0) && !isTightMuon(1) ){ // Loose-loose
					mumu.nthistos.h_nt0    ->Fill(MuPt[0], MuPt[1]);
					mumu.nthistos.h_nt0_pt ->Fill(MuPt[0]);
					mumu.nthistos.h_nt0_eta->Fill(MuEta[0]);
				}
			}
			if(isSigSupMuEventTRG()){ // f Ratio
				if( isLooseMuon(0) ){
					mumu.fhistos.h_nloose    ->Fill(MuPt[0], MuEta[0]);
					mumu.fhistos.h_nloose_pt ->Fill(MuPt[0]);
					mumu.fhistos.h_nloose_eta->Fill(MuEta[0]);
				}
				if( isTightMuon(0) ){
					mumu.fhistos.h_ntight    ->Fill(MuPt[0], MuEta[0]);
					mumu.fhistos.h_ntight_pt ->Fill(MuPt[0]);
					mumu.fhistos.h_ntight_eta->Fill(MuEta[0]);
				}
			}
			if(isZMuMuEventTRG()){ // p Ratio
				if( isLooseMuon(0) ){
					mumu.phistos.h_nloose    ->Fill(MuPt[0], MuEta[0]);
					mumu.phistos.h_nloose_pt ->Fill(MuPt[0]);
					mumu.phistos.h_nloose_eta->Fill(MuEta[0]);
				}
				if( isTightMuon(0) ){
					mumu.phistos.h_ntight    ->Fill(MuPt[0], MuEta[0]);
					mumu.phistos.h_ntight_pt ->Fill(MuPt[0]);
					mumu.phistos.h_ntight_eta->Fill(MuEta[0]);
				}
			}				
			
			// EE Channel
			if(isSSLLElEventTRG()){
				if(  isTightElectron(0) &&  isTightElectron(1) ){ // Tight-tight
					ee.nthistos.h_nt2    ->Fill(ElPt[0], ElPt[1]);
					ee.nthistos.h_nt2_pt ->Fill(ElPt[0]);
					ee.nthistos.h_nt2_eta->Fill(ElEta[0]);
				}
				if(  isTightElectron(0) && !isTightElectron(1) ){ // Tight-loose
					ee.nthistos.h_nt10    ->Fill(ElPt[0], ElPt[1]);
					ee.nthistos.h_nt10_pt ->Fill(ElPt[0]);
					ee.nthistos.h_nt10_eta->Fill(ElEta[0]);
				}
				if( !isTightElectron(0) &&  isTightElectron(1) ){ // Loose-tight
					ee.nthistos.h_nt10    ->Fill(ElPt[1], ElPt[0]);
					ee.nthistos.h_nt10_pt ->Fill(ElPt[1]);
					ee.nthistos.h_nt10_eta->Fill(ElEta[1]);
				}
				if( !isTightElectron(0) && !isTightElectron(1) ){ // Loose-loose
					ee.nthistos.h_nt0    ->Fill(ElPt[0], ElPt[1]);
					ee.nthistos.h_nt0_pt ->Fill(ElPt[0]);
					ee.nthistos.h_nt0_eta->Fill(ElEta[0]);
				}
			}
			if(isSigSupElEventTRG()){ // f Ratio
				if( isLooseElectron(0) ){
					ee.fhistos.h_nloose    ->Fill(ElPt[0], ElEta[0]);
					ee.fhistos.h_nloose_pt ->Fill(ElPt[0]);
					ee.fhistos.h_nloose_eta->Fill(ElEta[0]);
				}
				if( isTightElectron(0) ){
					ee.fhistos.h_ntight    ->Fill(ElPt[0], ElEta[0]);
					ee.fhistos.h_ntight_pt ->Fill(ElPt[0]);
					ee.fhistos.h_ntight_eta->Fill(ElEta[0]);
				}
			}
			if(isZElElEventTRG()){ // p Ratio
				if( isLooseElectron(0) ){
					ee.phistos.h_nloose    ->Fill(ElPt[0], ElEta[0]);
					ee.phistos.h_nloose_pt ->Fill(ElPt[0]);
					ee.phistos.h_nloose_eta->Fill(ElEta[0]);
				}
				if( isTightElectron(0) ){
					ee.phistos.h_ntight    ->Fill(ElPt[0], ElEta[0]);
					ee.phistos.h_ntight_pt ->Fill(ElPt[0]);
					ee.phistos.h_ntight_eta->Fill(ElEta[0]);
				}
			}
						
			// EMu Channel
			if(isSSLLElMuEventTRG()){
				if(  isTightElectron(0) &&  isTightMuon(0) ){ // Tight-tight
					emu.nthistos.h_nt2    ->Fill(MuPt [0], ElPt[0]);
					emu.nthistos.h_nt2_pt ->Fill(MuPt [0]);
					emu.nthistos.h_nt2_eta->Fill(MuEta[0]);
				}
				if( !isTightElectron(0) &&  isTightMuon(0) ){ // Tight-loose
					emu.nthistos.h_nt10    ->Fill(MuPt [0], ElPt[0]);
					emu.nthistos.h_nt10_pt ->Fill(MuPt [0]);
					emu.nthistos.h_nt10_eta->Fill(MuEta[0]);
				}
				if(  isTightElectron(0) && !isTightMuon(0) ){ // Loose-tight
					emu.nthistos.h_nt01    ->Fill(MuPt [0], ElPt[0]);
					emu.nthistos.h_nt01_pt ->Fill(MuPt [0]);
					emu.nthistos.h_nt01_eta->Fill(MuEta[0]);
				}
				if( !isTightElectron(0) && !isTightMuon(0) ){ // Loose-loose
					emu.nthistos.h_nt0    ->Fill(ElPt [0], MuPt[0]);
					emu.nthistos.h_nt0_pt ->Fill(MuPt [0]);
					emu.nthistos.h_nt0_eta->Fill(MuEta[0]);
				}
			}
		}
		
		// Calculate ratios
		mumu.fhistos.h_ratio    ->Divide(mumu.fhistos.h_ntight    , mumu.fhistos.h_nloose);
		mumu.fhistos.h_ratio_pt ->Divide(mumu.fhistos.h_ntight_pt , mumu.fhistos.h_nloose_pt);
		mumu.fhistos.h_ratio_eta->Divide(mumu.fhistos.h_ntight_eta, mumu.fhistos.h_nloose_eta);

		ee.fhistos.h_ratio    ->Divide(ee.fhistos.h_ntight    , ee.fhistos.h_nloose);
		ee.fhistos.h_ratio_pt ->Divide(ee.fhistos.h_ntight_pt , ee.fhistos.h_nloose_pt);
		ee.fhistos.h_ratio_eta->Divide(ee.fhistos.h_ntight_eta, ee.fhistos.h_nloose_eta);

		numberset numbers;
		numbers.nt2  = mumu.nthistos.h_nt2->GetEntries();
		numbers.nt10 = mumu.nthistos.h_nt10->GetEntries();
		numbers.nt01 = mumu.nthistos.h_nt01->GetEntries();
		numbers.nt0  = mumu.nthistos.h_nt0->GetEntries();
		numbers.nsst = mumu.fhistos.h_ntight->GetEntries();
		numbers.nssl = mumu.fhistos.h_nloose->GetEntries();
		numbers.nzt  = mumu.phistos.h_ntight->GetEntries();
		numbers.nzl  = mumu.phistos.h_nloose->GetEntries();
		fSamples[index].mumu.numbers = numbers;

		numbers.nt2  = emu.nthistos.h_nt2->GetEntries();
		numbers.nt10 = emu.nthistos.h_nt10->GetEntries();
		numbers.nt01 = emu.nthistos.h_nt01->GetEntries();
		numbers.nt0  = emu.nthistos.h_nt0->GetEntries();
		numbers.nsst = emu.fhistos.h_ntight->GetEntries();
		numbers.nssl = emu.fhistos.h_nloose->GetEntries();
		numbers.nzt  = emu.phistos.h_ntight->GetEntries();
		numbers.nzl  = emu.phistos.h_nloose->GetEntries();
		fSamples[index].emu.numbers = numbers;

		numbers.nt2  = ee.nthistos.h_nt2->GetEntries();
		numbers.nt10 = ee.nthistos.h_nt10->GetEntries();
		numbers.nt01 = ee.nthistos.h_nt01->GetEntries();
		numbers.nt0  = ee.nthistos.h_nt0->GetEntries();
		numbers.nsst = ee.fhistos.h_ntight->GetEntries();
		numbers.nssl = ee.fhistos.h_nloose->GetEntries();
		numbers.nzt  = ee.phistos.h_ntight->GetEntries();
		numbers.nzl  = ee.phistos.h_nloose->GetEntries();
		fSamples[index].ee.numbers = numbers;
	}
	writeHistos();
}

//____________________________________________________________________________
void MuonPlotter::printYields(){ printYields(fAllSamples); }
void MuonPlotter::printYields(float scale){ printYields(fAllSamples, scale); }
void MuonPlotter::printYields(int sample){ vector<int> samples; samples.push_back(sample); printYields(samples); }
void MuonPlotter::printYields(vector<int> samples, float lumiscale){
	cout << "-----------------------" << endl;
	if(lumiscale > -1.0) cout << " Numbers scaled to " << lumiscale << " /pb" << endl;
	for(size_t i = 0; i < samples.size(); ++i){
		int index = samples[i];
		channel cha = fSamples[index].mumu;
		numberset numbers = cha.numbers;
		cout << " Sample: " << fSamples[index].sname << endl;
		float scale = lumiscale / fSamples[index].lumi;
		if(scale < 0) scale = 1;
		if(fSamples[index].isdata) scale = 1;
		cout << "   Nt2 = "       << scale*numbers.nt2 <<  "  Nt1 = "       << scale*numbers.nt10 << "  Nt0 = " << scale*numbers.nt0 << endl;
		cout << "   Nss tight = " << scale*numbers.nsst << "  Nss loose = " << scale*numbers.nssl << endl;
		cout << "   NZ tight  = " << scale*numbers.nzt <<  "  NZ  loose = " << scale*numbers.nzl << endl;		
		cout << endl;
	}

	cout << "-----------------------" << endl;
}

void MuonPlotter::bookHistos(){
	for(size_t i = 0; i < fSamples.size(); ++i){
		TString name = fSamples[i].sname;
		fSamples[i].mumu.nthistos.h_nt2       = new TH2D(name + "_MuMu_NT2",        "NT2",        gNPt2bins,  gPt2bins, gNPt2bins, gPt2bins);
		fSamples[i].mumu.nthistos.h_nt2_pt    = new TH1D(name + "_MuMu_NT2_pt",     "NT2",        gNPt2bins,  gPt2bins);
		fSamples[i].mumu.nthistos.h_nt2_eta   = new TH1D(name + "_MuMu_NT2_eta",    "NT2",        gNEtabins, gEtabins);
		fSamples[i].mumu.nthistos.h_nt10      = new TH2D(name + "_MuMu_NT10",       "NT10",       gNPt2bins,  gPt2bins, gNPt2bins, gPt2bins);
		fSamples[i].mumu.nthistos.h_nt10_pt   = new TH1D(name + "_MuMu_NT10_pt",    "NT10 vs pt", gNPt2bins,  gPt2bins);
		fSamples[i].mumu.nthistos.h_nt10_eta  = new TH1D(name + "_MuMu_NT10_eta",   "NT10 vs eta",gNEtabins, gEtabins);
		fSamples[i].mumu.nthistos.h_nt01      = new TH2D(name + "_MuMu_NT01",       "NT01",       gNPt2bins,  gPt2bins, gNPt2bins, gPt2bins);
		fSamples[i].mumu.nthistos.h_nt01_pt   = new TH1D(name + "_MuMu_NT01_pt",    "NT01 vs pt", gNPt2bins,  gPt2bins);
		fSamples[i].mumu.nthistos.h_nt01_eta  = new TH1D(name + "_MuMu_NT01_eta",   "NT01 vs eta",gNEtabins, gEtabins);
		fSamples[i].mumu.nthistos.h_nt0       = new TH2D(name + "_MuMu_NT0",        "NT0",        gNPt2bins,  gPt2bins, gNPt2bins, gPt2bins);
		fSamples[i].mumu.nthistos.h_nt0_pt    = new TH1D(name + "_MuMu_NT0_pt",     "NT0 vs pt",  gNPt2bins,  gPt2bins);
		fSamples[i].mumu.nthistos.h_nt0_eta   = new TH1D(name + "_MuMu_NT0_eta",    "NT0 vs eta", gNEtabins, gEtabins);
		fSamples[i].mumu.fhistos.h_ntight     = new TH2D(name + "_MuMu_fTight",     "NTight Muons for sig. supp. selection",      gNPt2bins, gPt2bins, gNEtabins, gEtabins);
		fSamples[i].mumu.fhistos.h_nloose     = new TH2D(name + "_MuMu_fLoose",     "NLoose Muons for sig. supp. selection",      gNPt2bins, gPt2bins, gNEtabins, gEtabins);
		fSamples[i].mumu.fhistos.h_ntight_pt  = new TH1D(name + "_MuMu_fTight_pt",  "NTight Muons for sig. supp. selection",      gNPt2bins, gPt2bins);
		fSamples[i].mumu.fhistos.h_nloose_pt  = new TH1D(name + "_MuMu_fLoose_pt",  "NLoose Muons for sig. supp. selection",      gNPt2bins, gPt2bins);
		fSamples[i].mumu.fhistos.h_ntight_eta = new TH1D(name + "_MuMu_fTight_eta", "NTight Muons for sig. supp. selection",      gNEtabins, gEtabins);
		fSamples[i].mumu.fhistos.h_nloose_eta = new TH1D(name + "_MuMu_fLoose_eta", "NLoose Muons for sig. supp. selection",      gNEtabins, gEtabins);
		fSamples[i].mumu.fhistos.h_ratio      = new TH2D(name + "_MuMu_fRatio",     "Tight/Loose Ratio for sig. supp. selection", gNPt2bins, gPt2bins, gNEtabins, gEtabins);
		fSamples[i].mumu.fhistos.h_ratio_pt   = new TH1D(name + "_MuMu_fRatio_pt",  "Tight/Loose Ratio for sig. supp. selection", gNPt2bins, gPt2bins);
		fSamples[i].mumu.fhistos.h_ratio_eta  = new TH1D(name + "_MuMu_fRatio_eta", "Tight/Loose Ratio for sig. supp. selection", gNEtabins, gEtabins);
		fSamples[i].mumu.phistos.h_ntight     = new TH2D(name + "_MuMu_pTight",     "NTight Muons for Z decay selection",      gNPt2bins, gPt2bins, gNEtabins, gEtabins);
		fSamples[i].mumu.phistos.h_nloose     = new TH2D(name + "_MuMu_pLoose",     "NLoose Muons for Z decay selection",      gNPt2bins, gPt2bins, gNEtabins, gEtabins);
		fSamples[i].mumu.phistos.h_ntight_pt  = new TH1D(name + "_MuMu_pTight_pt",  "NTight Muons for Z decay selection",      gNPt2bins, gPt2bins);
		fSamples[i].mumu.phistos.h_nloose_pt  = new TH1D(name + "_MuMu_pLoose_pt",  "NLoose Muons for Z decay selection",      gNPt2bins, gPt2bins);
		fSamples[i].mumu.phistos.h_ntight_eta = new TH1D(name + "_MuMu_pTight_eta", "NTight Muons for Z decay selection",      gNEtabins, gEtabins);
		fSamples[i].mumu.phistos.h_nloose_eta = new TH1D(name + "_MuMu_pLoose_eta", "NLoose Muons for Z decay selection",      gNEtabins, gEtabins);		
		fSamples[i].mumu.phistos.h_ratio      = new TH2D(name + "_MuMu_pRatio",     "Tight/Loose Ratio for Z decay selection", gNPt2bins, gPt2bins, gNEtabins, gEtabins);
		fSamples[i].mumu.phistos.h_ratio_pt   = new TH1D(name + "_MuMu_pRatio_pt",  "Tight/Loose Ratio for Z decay selection", gNPt2bins, gPt2bins);
		fSamples[i].mumu.phistos.h_ratio_eta  = new TH1D(name + "_MuMu_pRatio_eta", "Tight/Loose Ratio for Z decay selection", gNEtabins, gEtabins);

		fSamples[i].mumu.nthistos.h_nt2->Sumw2();  fSamples[i].mumu.nthistos.h_nt2_pt->Sumw2();  fSamples[i].mumu.nthistos.h_nt2_eta->Sumw2();
		fSamples[i].mumu.nthistos.h_nt10->Sumw2(); fSamples[i].mumu.nthistos.h_nt10_pt->Sumw2(); fSamples[i].mumu.nthistos.h_nt10_eta->Sumw2();
		fSamples[i].mumu.nthistos.h_nt01->Sumw2(); fSamples[i].mumu.nthistos.h_nt01_pt->Sumw2(); fSamples[i].mumu.nthistos.h_nt01_eta->Sumw2();
		fSamples[i].mumu.nthistos.h_nt0->Sumw2();  fSamples[i].mumu.nthistos.h_nt0_pt->Sumw2();  fSamples[i].mumu.nthistos.h_nt0_eta->Sumw2();

		fSamples[i].mumu.fhistos.h_ntight     ->Sumw2(); fSamples[i].mumu.fhistos.h_nloose     ->Sumw2();
		fSamples[i].mumu.fhistos.h_ntight_pt  ->Sumw2(); fSamples[i].mumu.fhistos.h_nloose_pt  ->Sumw2();
		fSamples[i].mumu.fhistos.h_ntight_eta ->Sumw2(); fSamples[i].mumu.fhistos.h_nloose_eta ->Sumw2();
		fSamples[i].mumu.fhistos.h_ratio_pt   ->Sumw2(); fSamples[i].mumu.fhistos.h_ratio_eta  ->Sumw2();
		fSamples[i].mumu.fhistos.h_ratio      ->Sumw2();

		fSamples[i].mumu.phistos.h_ntight     ->Sumw2(); fSamples[i].mumu.phistos.h_nloose     ->Sumw2();
		fSamples[i].mumu.phistos.h_ntight_pt  ->Sumw2(); fSamples[i].mumu.phistos.h_nloose_pt  ->Sumw2();
		fSamples[i].mumu.phistos.h_ntight_eta ->Sumw2(); fSamples[i].mumu.phistos.h_nloose_eta ->Sumw2();
		fSamples[i].mumu.phistos.h_ratio_pt   ->Sumw2(); fSamples[i].mumu.phistos.h_ratio_eta  ->Sumw2();
		fSamples[i].mumu.phistos.h_ratio      ->Sumw2(); 

		fSamples[i].emu.nthistos.h_nt2       = new TH2D(name + "_EMu_NT2",        "NT2",        gNPt2bins,  gPt2bins, gNPt2bins, gPt2bins);
		fSamples[i].emu.nthistos.h_nt2_pt    = new TH1D(name + "_EMu_NT2_pt",     "NT2",        gNPt2bins,  gPt2bins);
		fSamples[i].emu.nthistos.h_nt2_eta   = new TH1D(name + "_EMu_NT2_eta",    "NT2",        gNEtabins, gEtabins);
		fSamples[i].emu.nthistos.h_nt10      = new TH2D(name + "_EMu_NT10",       "NT10",       gNPt2bins,  gPt2bins, gNPt2bins, gPt2bins);
		fSamples[i].emu.nthistos.h_nt10_pt   = new TH1D(name + "_EMu_NT10_pt",    "NT10 vs pt", gNPt2bins,  gPt2bins);
		fSamples[i].emu.nthistos.h_nt10_eta  = new TH1D(name + "_EMu_NT10_eta",   "NT10 vs eta",gNEtabins, gEtabins);
		fSamples[i].emu.nthistos.h_nt01      = new TH2D(name + "_EMu_NT01",       "NT01",       gNPt2bins,  gPt2bins, gNPt2bins, gPt2bins);
		fSamples[i].emu.nthistos.h_nt01_pt   = new TH1D(name + "_EMu_NT01_pt",    "NT01 vs pt", gNPt2bins,  gPt2bins);
		fSamples[i].emu.nthistos.h_nt01_eta  = new TH1D(name + "_EMu_NT01_eta",   "NT01 vs eta",gNEtabins, gEtabins);
		fSamples[i].emu.nthistos.h_nt0       = new TH2D(name + "_EMu_NT0",        "NT0",        gNPt2bins,  gPt2bins, gNPt2bins, gPt2bins);
		fSamples[i].emu.nthistos.h_nt0_pt    = new TH1D(name + "_EMu_NT0_pt",     "NT0 vs pt",  gNPt2bins,  gPt2bins);
		fSamples[i].emu.nthistos.h_nt0_eta   = new TH1D(name + "_EMu_NT0_eta",    "NT0 vs eta", gNEtabins, gEtabins);
		fSamples[i].emu.fhistos.h_ntight     = new TH2D(name + "_EMu_fTight",     "NTight Muons for sig. supp. selection",      gNPt2bins, gPt2bins, gNEtabins, gEtabins);
		fSamples[i].emu.fhistos.h_nloose     = new TH2D(name + "_EMu_fLoose",     "NLoose Muons for sig. supp. selection",      gNPt2bins, gPt2bins, gNEtabins, gEtabins);
		fSamples[i].emu.fhistos.h_ntight_pt  = new TH1D(name + "_EMu_fTight_pt",  "NTight Muons for sig. supp. selection",      gNPt2bins, gPt2bins);
		fSamples[i].emu.fhistos.h_nloose_pt  = new TH1D(name + "_EMu_fLoose_pt",  "NLoose Muons for sig. supp. selection",      gNPt2bins, gPt2bins);
		fSamples[i].emu.fhistos.h_ntight_eta = new TH1D(name + "_EMu_fTight_eta", "NTight Muons for sig. supp. selection",      gNEtabins, gEtabins);
		fSamples[i].emu.fhistos.h_nloose_eta = new TH1D(name + "_EMu_fLoose_eta", "NLoose Muons for sig. supp. selection",      gNEtabins, gEtabins);
		fSamples[i].emu.fhistos.h_ratio      = new TH2D(name + "_EMu_fRatio",     "Tight/Loose Ratio for sig. supp. selection", gNPt2bins, gPt2bins, gNEtabins, gEtabins);
		fSamples[i].emu.fhistos.h_ratio_pt   = new TH1D(name + "_EMu_fRatio_pt",  "Tight/Loose Ratio for sig. supp. selection", gNPt2bins, gPt2bins);
		fSamples[i].emu.fhistos.h_ratio_eta  = new TH1D(name + "_EMu_fRatio_eta", "Tight/Loose Ratio for sig. supp. selection", gNEtabins, gEtabins);
		fSamples[i].emu.phistos.h_ntight     = new TH2D(name + "_EMu_pTight",     "NTight Muons for Z decay selection",      gNPt2bins, gPt2bins, gNEtabins, gEtabins);
		fSamples[i].emu.phistos.h_nloose     = new TH2D(name + "_EMu_pLoose",     "NLoose Muons for Z decay selection",      gNPt2bins, gPt2bins, gNEtabins, gEtabins);
		fSamples[i].emu.phistos.h_ntight_pt  = new TH1D(name + "_EMu_pTight_pt",  "NTight Muons for Z decay selection",      gNPt2bins, gPt2bins);
		fSamples[i].emu.phistos.h_nloose_pt  = new TH1D(name + "_EMu_pLoose_pt",  "NLoose Muons for Z decay selection",      gNPt2bins, gPt2bins);
		fSamples[i].emu.phistos.h_ntight_eta = new TH1D(name + "_EMu_pTight_eta", "NTight Muons for Z decay selection",      gNEtabins, gEtabins);
		fSamples[i].emu.phistos.h_nloose_eta = new TH1D(name + "_EMu_pLoose_eta", "NLoose Muons for Z decay selection",      gNEtabins, gEtabins);		
		fSamples[i].emu.phistos.h_ratio      = new TH2D(name + "_EMu_pRatio",     "Tight/Loose Ratio for Z decay selection", gNPt2bins, gPt2bins, gNEtabins, gEtabins);
		fSamples[i].emu.phistos.h_ratio_pt   = new TH1D(name + "_EMu_pRatio_pt",  "Tight/Loose Ratio for Z decay selection", gNPt2bins, gPt2bins);
		fSamples[i].emu.phistos.h_ratio_eta  = new TH1D(name + "_EMu_pRatio_eta", "Tight/Loose Ratio for Z decay selection", gNEtabins, gEtabins);

		fSamples[i].emu.nthistos.h_nt2->Sumw2();       fSamples[i].emu.nthistos.h_nt2_pt->Sumw2();  fSamples[i].emu.nthistos.h_nt2_eta->Sumw2();
		fSamples[i].emu.nthistos.h_nt10->Sumw2();      fSamples[i].emu.nthistos.h_nt10_pt->Sumw2(); fSamples[i].emu.nthistos.h_nt10_eta->Sumw2();
		fSamples[i].emu.nthistos.h_nt01->Sumw2();      fSamples[i].emu.nthistos.h_nt01_pt->Sumw2(); fSamples[i].emu.nthistos.h_nt01_eta->Sumw2();
		fSamples[i].emu.nthistos.h_nt0->Sumw2();       fSamples[i].emu.nthistos.h_nt0_pt->Sumw2();  fSamples[i].emu.nthistos.h_nt0_eta->Sumw2();

		fSamples[i].emu.fhistos.h_ntight    ->Sumw2(); fSamples[i].emu.fhistos.h_nloose     ->Sumw2();
		fSamples[i].emu.fhistos.h_ntight_pt ->Sumw2(); fSamples[i].emu.fhistos.h_nloose_pt  ->Sumw2();
		fSamples[i].emu.fhistos.h_ntight_eta->Sumw2(); fSamples[i].emu.fhistos.h_nloose_eta ->Sumw2();
		fSamples[i].emu.fhistos.h_ratio_pt  ->Sumw2(); fSamples[i].emu.fhistos.h_ratio_eta  ->Sumw2();
		fSamples[i].emu.fhistos.h_ratio     ->Sumw2();

		fSamples[i].emu.phistos.h_ntight    ->Sumw2(); fSamples[i].emu.phistos.h_nloose     ->Sumw2();
		fSamples[i].emu.phistos.h_ntight_pt ->Sumw2(); fSamples[i].emu.phistos.h_nloose_pt  ->Sumw2();
		fSamples[i].emu.phistos.h_ntight_eta->Sumw2(); fSamples[i].emu.phistos.h_nloose_eta ->Sumw2();
		fSamples[i].emu.phistos.h_ratio_pt  ->Sumw2(); fSamples[i].emu.phistos.h_ratio_eta  ->Sumw2();
		fSamples[i].emu.phistos.h_ratio     ->Sumw2(); 

		fSamples[i].ee.nthistos.h_nt2       = new TH2D(name + "_EE_NT2",        "NT2",        gNPt2bins,  gPt2bins, gNPt2bins, gPt2bins);
		fSamples[i].ee.nthistos.h_nt2_pt    = new TH1D(name + "_EE_NT2_pt",     "NT2",        gNPt2bins,  gPt2bins);
		fSamples[i].ee.nthistos.h_nt2_eta   = new TH1D(name + "_EE_NT2_eta",    "NT2",        gNEtabins, gEtabins);
		fSamples[i].ee.nthistos.h_nt10      = new TH2D(name + "_EE_NT10",       "NT10",       gNPt2bins,  gPt2bins, gNPt2bins, gPt2bins);
		fSamples[i].ee.nthistos.h_nt10_pt   = new TH1D(name + "_EE_NT10_pt",    "NT10 vs pt", gNPt2bins,  gPt2bins);
		fSamples[i].ee.nthistos.h_nt10_eta  = new TH1D(name + "_EE_NT10_eta",   "NT10 vs eta",gNEtabins, gEtabins);
		fSamples[i].ee.nthistos.h_nt01      = new TH2D(name + "_EE_NT01",       "NT01",       gNPt2bins,  gPt2bins, gNPt2bins, gPt2bins);
		fSamples[i].ee.nthistos.h_nt01_pt   = new TH1D(name + "_EE_NT01_pt",    "NT01 vs pt", gNPt2bins,  gPt2bins);
		fSamples[i].ee.nthistos.h_nt01_eta  = new TH1D(name + "_EE_NT01_eta",   "NT01 vs eta",gNEtabins, gEtabins);
		fSamples[i].ee.nthistos.h_nt0       = new TH2D(name + "_EE_NT0",        "NT0",        gNPt2bins,  gPt2bins, gNPt2bins, gPt2bins);
		fSamples[i].ee.nthistos.h_nt0_pt    = new TH1D(name + "_EE_NT0_pt",     "NT0 vs pt",  gNPt2bins,  gPt2bins);
		fSamples[i].ee.nthistos.h_nt0_eta   = new TH1D(name + "_EE_NT0_eta",    "NT0 vs eta", gNEtabins, gEtabins);
		fSamples[i].ee.fhistos.h_ntight     = new TH2D(name + "_EE_fTight",     "NTight Muons for sig. supp. selection",      gNPt2bins, gPt2bins, gNEtabins, gEtabins);
		fSamples[i].ee.fhistos.h_nloose     = new TH2D(name + "_EE_fLoose",     "NLoose Muons for sig. supp. selection",      gNPt2bins, gPt2bins, gNEtabins, gEtabins);
		fSamples[i].ee.fhistos.h_ntight_pt  = new TH1D(name + "_EE_fTight_pt",  "NTight Muons for sig. supp. selection",      gNPt2bins, gPt2bins);
		fSamples[i].ee.fhistos.h_nloose_pt  = new TH1D(name + "_EE_fLoose_pt",  "NLoose Muons for sig. supp. selection",      gNPt2bins, gPt2bins);
		fSamples[i].ee.fhistos.h_ntight_eta = new TH1D(name + "_EE_fTight_eta", "NTight Muons for sig. supp. selection",      gNEtabins, gEtabins);
		fSamples[i].ee.fhistos.h_nloose_eta = new TH1D(name + "_EE_fLoose_eta", "NLoose Muons for sig. supp. selection",      gNEtabins, gEtabins);
		fSamples[i].ee.fhistos.h_ratio      = new TH2D(name + "_EE_fRatio",     "Tight/Loose Ratio for sig. supp. selection", gNPt2bins, gPt2bins, gNEtabins, gEtabins);
		fSamples[i].ee.fhistos.h_ratio_pt   = new TH1D(name + "_EE_fRatio_pt",  "Tight/Loose Ratio for sig. supp. selection", gNPt2bins, gPt2bins);
		fSamples[i].ee.fhistos.h_ratio_eta  = new TH1D(name + "_EE_fRatio_eta", "Tight/Loose Ratio for sig. supp. selection", gNEtabins, gEtabins);
		fSamples[i].ee.phistos.h_ntight     = new TH2D(name + "_EE_pTight",     "NTight Muons for Z decay selection",      gNPt2bins, gPt2bins, gNEtabins, gEtabins);
		fSamples[i].ee.phistos.h_nloose     = new TH2D(name + "_EE_pLoose",     "NLoose Muons for Z decay selection",      gNPt2bins, gPt2bins, gNEtabins, gEtabins);
		fSamples[i].ee.phistos.h_ntight_pt  = new TH1D(name + "_EE_pTight_pt",  "NTight Muons for Z decay selection",      gNPt2bins, gPt2bins);
		fSamples[i].ee.phistos.h_nloose_pt  = new TH1D(name + "_EE_pLoose_pt",  "NLoose Muons for Z decay selection",      gNPt2bins, gPt2bins);
		fSamples[i].ee.phistos.h_ntight_eta = new TH1D(name + "_EE_pTight_eta", "NTight Muons for Z decay selection",      gNEtabins, gEtabins);
		fSamples[i].ee.phistos.h_nloose_eta = new TH1D(name + "_EE_pLoose_eta", "NLoose Muons for Z decay selection",      gNEtabins, gEtabins);		
		fSamples[i].ee.phistos.h_ratio      = new TH2D(name + "_EE_pRatio",     "Tight/Loose Ratio for Z decay selection", gNPt2bins, gPt2bins, gNEtabins, gEtabins);
		fSamples[i].ee.phistos.h_ratio_pt   = new TH1D(name + "_EE_pRatio_pt",  "Tight/Loose Ratio for Z decay selection", gNPt2bins, gPt2bins);
		fSamples[i].ee.phistos.h_ratio_eta  = new TH1D(name + "_EE_pRatio_eta", "Tight/Loose Ratio for Z decay selection", gNEtabins, gEtabins);

		fSamples[i].ee.nthistos.h_nt2->Sumw2();       fSamples[i].ee.nthistos.h_nt2_pt->Sumw2();  fSamples[i].ee.nthistos.h_nt2_eta->Sumw2();
		fSamples[i].ee.nthistos.h_nt10->Sumw2();      fSamples[i].ee.nthistos.h_nt10_pt->Sumw2(); fSamples[i].ee.nthistos.h_nt10_eta->Sumw2();
		fSamples[i].ee.nthistos.h_nt01->Sumw2();      fSamples[i].ee.nthistos.h_nt01_pt->Sumw2(); fSamples[i].ee.nthistos.h_nt01_eta->Sumw2();
		fSamples[i].ee.nthistos.h_nt0->Sumw2();       fSamples[i].ee.nthistos.h_nt0_pt->Sumw2();  fSamples[i].ee.nthistos.h_nt0_eta->Sumw2();

		fSamples[i].ee.fhistos.h_ntight    ->Sumw2(); fSamples[i].ee.fhistos.h_nloose     ->Sumw2();
		fSamples[i].ee.fhistos.h_ntight_pt ->Sumw2(); fSamples[i].ee.fhistos.h_nloose_pt  ->Sumw2();
		fSamples[i].ee.fhistos.h_ntight_eta->Sumw2(); fSamples[i].ee.fhistos.h_nloose_eta ->Sumw2();
		fSamples[i].ee.fhistos.h_ratio_pt  ->Sumw2(); fSamples[i].ee.fhistos.h_ratio_eta  ->Sumw2();
		fSamples[i].ee.fhistos.h_ratio     ->Sumw2();

		fSamples[i].ee.phistos.h_ntight    ->Sumw2(); fSamples[i].ee.phistos.h_nloose     ->Sumw2();
		fSamples[i].ee.phistos.h_ntight_pt ->Sumw2(); fSamples[i].ee.phistos.h_nloose_pt  ->Sumw2();
		fSamples[i].ee.phistos.h_ntight_eta->Sumw2(); fSamples[i].ee.phistos.h_nloose_eta ->Sumw2();
		fSamples[i].ee.phistos.h_ratio_pt  ->Sumw2(); fSamples[i].ee.phistos.h_ratio_eta  ->Sumw2();
		fSamples[i].ee.phistos.h_ratio     ->Sumw2();
	}
}

void MuonPlotter::writeHistos(){
	TFile *pFile = new TFile(fOutputFileName, "RECREATE");
	pFile->cd();
	for(size_t i = 0; i < fSamples.size(); ++i){
		TDirectory* cdir = Util::FindOrCreate(fSamples[i].sname, pFile);
		cdir->cd();

		for(size_t ch = 0; ch < 3; ++ch){ // Loop over channels, mumu, emu, ee
			channel cha;
			if(ch == Muon)     cha = fSamples[i].mumu;
			if(ch == EMu)      cha = fSamples[i].emu;
			if(ch == Electron) cha = fSamples[i].ee;
			cha.nthistos.h_nt2      ->Write(cha.nthistos.h_nt2      ->GetName(), TObject::kWriteDelete);
			cha.nthistos.h_nt2_pt   ->Write(cha.nthistos.h_nt2_pt   ->GetName(), TObject::kWriteDelete);
			cha.nthistos.h_nt2_eta  ->Write(cha.nthistos.h_nt2_eta  ->GetName(), TObject::kWriteDelete);
			cha.nthistos.h_nt10     ->Write(cha.nthistos.h_nt10     ->GetName(), TObject::kWriteDelete);
			cha.nthistos.h_nt10_pt  ->Write(cha.nthistos.h_nt10_pt  ->GetName(), TObject::kWriteDelete);
			cha.nthistos.h_nt10_eta ->Write(cha.nthistos.h_nt10_eta ->GetName(), TObject::kWriteDelete);
			cha.nthistos.h_nt01     ->Write(cha.nthistos.h_nt01     ->GetName(), TObject::kWriteDelete);
			cha.nthistos.h_nt01_pt  ->Write(cha.nthistos.h_nt01_pt  ->GetName(), TObject::kWriteDelete);
			cha.nthistos.h_nt01_eta ->Write(cha.nthistos.h_nt01_eta ->GetName(), TObject::kWriteDelete);
			cha.nthistos.h_nt0      ->Write(cha.nthistos.h_nt0      ->GetName(), TObject::kWriteDelete);
			cha.nthistos.h_nt0_pt   ->Write(cha.nthistos.h_nt0_pt   ->GetName(), TObject::kWriteDelete);
			cha.nthistos.h_nt0_eta  ->Write(cha.nthistos.h_nt0_eta  ->GetName(), TObject::kWriteDelete);
			cha.fhistos.h_ntight    ->Write(cha.fhistos.h_ntight    ->GetName(), TObject::kWriteDelete);
			cha.fhistos.h_nloose    ->Write(cha.fhistos.h_nloose    ->GetName(), TObject::kWriteDelete);
			cha.fhistos.h_ntight_pt ->Write(cha.fhistos.h_ntight_pt ->GetName(), TObject::kWriteDelete);
			cha.fhistos.h_nloose_pt ->Write(cha.fhistos.h_nloose_pt ->GetName(), TObject::kWriteDelete);
			cha.fhistos.h_ntight_eta->Write(cha.fhistos.h_ntight_eta->GetName(), TObject::kWriteDelete);
			cha.fhistos.h_nloose_eta->Write(cha.fhistos.h_nloose_eta->GetName(), TObject::kWriteDelete);
			cha.fhistos.h_ratio     ->Write(cha.fhistos.h_ratio     ->GetName(), TObject::kWriteDelete);
			cha.fhistos.h_ratio_pt  ->Write(cha.fhistos.h_ratio_pt  ->GetName(), TObject::kWriteDelete);
			cha.fhistos.h_ratio_eta ->Write(cha.fhistos.h_ratio_eta ->GetName(), TObject::kWriteDelete);
			cha.phistos.h_ntight    ->Write(cha.phistos.h_ntight    ->GetName(), TObject::kWriteDelete);
			cha.phistos.h_nloose    ->Write(cha.phistos.h_nloose    ->GetName(), TObject::kWriteDelete);
			cha.phistos.h_ntight_pt ->Write(cha.phistos.h_ntight_pt ->GetName(), TObject::kWriteDelete);
			cha.phistos.h_nloose_pt ->Write(cha.phistos.h_nloose_pt ->GetName(), TObject::kWriteDelete);
			cha.phistos.h_ntight_eta->Write(cha.phistos.h_ntight_eta->GetName(), TObject::kWriteDelete);
			cha.phistos.h_nloose_eta->Write(cha.phistos.h_nloose_eta->GetName(), TObject::kWriteDelete);
			cha.phistos.h_ratio     ->Write(cha.phistos.h_ratio     ->GetName(), TObject::kWriteDelete);
			cha.phistos.h_ratio_pt  ->Write(cha.phistos.h_ratio_pt  ->GetName(), TObject::kWriteDelete);
			cha.phistos.h_ratio_eta ->Write(cha.phistos.h_ratio_eta ->GetName(), TObject::kWriteDelete);
		}
	}
	pFile->Write();
	pFile->Close();
}

int MuonPlotter::readHistos(TString filename){
	TFile *pFile = TFile::Open(filename, "READ");
	if(pFile == NULL){
		cout << "File " << filename << " does not exist!" << endl;
		return 1;
	}
	pFile->cd();
	for(size_t i = 0; i < fSamples.size(); ++i){
		for(size_t ch = 0; ch < 3; ++ch){ // Loop over channels, mumu, emu, ee
			channel *cha;
			if(ch == 0) cha = &fSamples[i].mumu;
			if(ch == 1) cha = &fSamples[i].emu;
			if(ch == 2) cha = &fSamples[i].ee;
		
			cha->nthistos.h_nt2       = (TH2D*)pFile->Get(fSamples[i].sname + "/" + cha->nthistos.h_nt2      ->GetName());
			cha->nthistos.h_nt2_pt    = (TH1D*)pFile->Get(fSamples[i].sname + "/" + cha->nthistos.h_nt2_pt   ->GetName());
			cha->nthistos.h_nt2_eta   = (TH1D*)pFile->Get(fSamples[i].sname + "/" + cha->nthistos.h_nt2_eta  ->GetName());
			cha->nthistos.h_nt10      = (TH2D*)pFile->Get(fSamples[i].sname + "/" + cha->nthistos.h_nt10     ->GetName());
			cha->nthistos.h_nt10_pt   = (TH1D*)pFile->Get(fSamples[i].sname + "/" + cha->nthistos.h_nt10_pt  ->GetName());
			cha->nthistos.h_nt10_eta  = (TH1D*)pFile->Get(fSamples[i].sname + "/" + cha->nthistos.h_nt10_eta ->GetName());
			cha->nthistos.h_nt01      = (TH2D*)pFile->Get(fSamples[i].sname + "/" + cha->nthistos.h_nt01     ->GetName());
			cha->nthistos.h_nt01_pt   = (TH1D*)pFile->Get(fSamples[i].sname + "/" + cha->nthistos.h_nt01_pt  ->GetName());
			cha->nthistos.h_nt01_eta  = (TH1D*)pFile->Get(fSamples[i].sname + "/" + cha->nthistos.h_nt01_eta ->GetName());
			cha->nthistos.h_nt0       = (TH2D*)pFile->Get(fSamples[i].sname + "/" + cha->nthistos.h_nt0      ->GetName());
			cha->nthistos.h_nt0_pt    = (TH1D*)pFile->Get(fSamples[i].sname + "/" + cha->nthistos.h_nt0_pt   ->GetName());
			cha->nthistos.h_nt0_eta   = (TH1D*)pFile->Get(fSamples[i].sname + "/" + cha->nthistos.h_nt0_eta  ->GetName());
			cha->fhistos.h_ntight     = (TH2D*)pFile->Get(fSamples[i].sname + "/" + cha->fhistos.h_ntight    ->GetName());
			cha->fhistos.h_nloose     = (TH2D*)pFile->Get(fSamples[i].sname + "/" + cha->fhistos.h_nloose    ->GetName());
			cha->fhistos.h_ntight_pt  = (TH1D*)pFile->Get(fSamples[i].sname + "/" + cha->fhistos.h_ntight_pt ->GetName());
			cha->fhistos.h_nloose_pt  = (TH1D*)pFile->Get(fSamples[i].sname + "/" + cha->fhistos.h_nloose_pt ->GetName());
			cha->fhistos.h_ntight_eta = (TH1D*)pFile->Get(fSamples[i].sname + "/" + cha->fhistos.h_ntight_eta->GetName());
			cha->fhistos.h_nloose_eta = (TH1D*)pFile->Get(fSamples[i].sname + "/" + cha->fhistos.h_nloose_eta->GetName());
			cha->fhistos.h_ratio      = (TH2D*)pFile->Get(fSamples[i].sname + "/" + cha->fhistos.h_ratio     ->GetName());
			cha->fhistos.h_ratio_pt   = (TH1D*)pFile->Get(fSamples[i].sname + "/" + cha->fhistos.h_ratio_pt  ->GetName());
			cha->fhistos.h_ratio_eta  = (TH1D*)pFile->Get(fSamples[i].sname + "/" + cha->fhistos.h_ratio_eta ->GetName());
			cha->phistos.h_ntight     = (TH2D*)pFile->Get(fSamples[i].sname + "/" + cha->phistos.h_ntight    ->GetName());
			cha->phistos.h_nloose     = (TH2D*)pFile->Get(fSamples[i].sname + "/" + cha->phistos.h_nloose    ->GetName());
			cha->phistos.h_ntight_pt  = (TH1D*)pFile->Get(fSamples[i].sname + "/" + cha->phistos.h_ntight_pt ->GetName());
			cha->phistos.h_nloose_pt  = (TH1D*)pFile->Get(fSamples[i].sname + "/" + cha->phistos.h_nloose_pt ->GetName());
			cha->phistos.h_ntight_eta = (TH1D*)pFile->Get(fSamples[i].sname + "/" + cha->phistos.h_ntight_eta->GetName());
			cha->phistos.h_nloose_eta = (TH1D*)pFile->Get(fSamples[i].sname + "/" + cha->phistos.h_nloose_eta->GetName());
			cha->phistos.h_ratio      = (TH2D*)pFile->Get(fSamples[i].sname + "/" + cha->phistos.h_ratio     ->GetName());
			cha->phistos.h_ratio_pt   = (TH1D*)pFile->Get(fSamples[i].sname + "/" + cha->phistos.h_ratio_pt  ->GetName());
			cha->phistos.h_ratio_eta  = (TH1D*)pFile->Get(fSamples[i].sname + "/" + cha->phistos.h_ratio_eta ->GetName());

			numberset numbers;
			numbers.nt2  = cha->nthistos.h_nt2->GetEntries();
			numbers.nt10 = cha->nthistos.h_nt10->GetEntries();
			numbers.nt01 = cha->nthistos.h_nt01->GetEntries();
			numbers.nt0  = cha->nthistos.h_nt0->GetEntries();
			numbers.nsst = cha->fhistos.h_ntight->GetEntries();
			numbers.nssl = cha->fhistos.h_nloose->GetEntries();
			numbers.nzt  = cha->phistos.h_ntight->GetEntries();
			numbers.nzl  = cha->phistos.h_nloose->GetEntries();

			cha->numbers = numbers;
		}
	}
	return 0;
}

//////////////////////////////////////////////////////////////////////////////
// Event Selections:
//____________________________________________________________________________
bool MuonPlotter::isGoodEvent(){
	// Some global cuts, select events with >1 jets
	if(!passesNJetCut(2)) return false;
	return true;
}

//____________________________________________________________________________
bool MuonPlotter::isGoodMuEvent(){
	// Ask for >0 loose muons, if 2 muons ask for second to be loose too
	if(!isGoodEvent()) return false;
	if(NMus < 1) return false;
	if(isLooseMuon(0) == false) return false;
	if(NMus > 1) if(isLooseMuon(1) == false) return false;
	return true;
}

//____________________________________________________________________________
bool MuonPlotter::isGoodElEvent(){
	// Ask for >0 loose electrons, if 2 electrons ask for second to be loose too
	if(!isGoodEvent()) return false;
	if(NEls < 1) return false;
	if(isLooseElectron(0) == false) return false;
	if(NEls > 1) if(isLooseElectron(1) == false) return false;
	return true;
}

//____________________________________________________________________________
bool MuonPlotter::isGoodElMuEvent(){
	// Ask for >0 loose electrons and muons
	if(!isGoodEvent()) return false;
	if(NEls < 1 || NMus < 1) return false;
	if(isLooseElectron(0) == false) return false;
	if(isLooseMuon(0) == false) return false;
	return true;
}

//____________________________________________________________________________
bool MuonPlotter::passesNJetCut(int cut){
	int njets(0);
	for(size_t i = 0; i < NJets; ++i) if(isGoodJet(i)) njets++;
	return njets>=cut;
}

//____________________________________________________________________________
bool MuonPlotter::passesHTCut(float cut){
	float ht(0.);
	for(size_t i = 0; i < NJets; ++i) if(isGoodJet(i)) ht += JetPt[i];
	return ht >= cut;
}

//____________________________________________________________________________
bool MuonPlotter::passesMETCut(float cut){
	if(pfMET < cut) return false;
	return true;
}

//____________________________________________________________________________
bool MuonPlotter::passesZVeto(float dm){
// Checks if any combination of opposite sign, same flavor tight leptons (e or mu)
// has invariant mass closer than dm to the Z mass, returns true if none found
	const float MMU = 0.1057;
	const float MEL = 0.0005;
	const float MZ  = 91.2;

	if(NMus > 1){
		// First muon
		for(size_t i = 0; i < NMus-1; ++i){
			if(isTightMuon(i)){
				TLorentzVector pmu1, pmu2;
				pmu1.SetPtEtaPhiM(MuPt[i], MuEta[i], MuPhi[i], MMU);

				// Second muon
				for(size_t j = i+1; j < NMus; ++j){ 
					if(isTightMuon(j) && (MuCharge[i] != MuCharge[j]) ){
						pmu2.SetPtEtaPhiM(MuPt[j], MuEta[j], MuPhi[j], MMU);
						if(fabs((pmu1+pmu2).M() - MZ) < dm) return false;
					}
				}
			}
		}
	}
	
	if(NEls > 1){
		// First electron
		for(size_t i = 0; i < NEls-1; ++i){
			if(isTightElectron(i)){
				TLorentzVector pel1, pel2;
				pel1.SetPtEtaPhiM(ElPt[i], ElEta[i], ElPhi[i], MEL);

				// Second electron
				for(size_t j = i+1; j < NEls; ++j){
					if(isTightElectron(j) && (ElCh[i] != ElCh[j]) ){
						pel2.SetPtEtaPhiM(ElPt[j], ElEta[j], ElPhi[j], MEL);
						if(fabs((pel1+pel2).M() - MZ) < dm) return false;
					}
				}
			}
		}		
	}
	return true;
}

//____________________________________________________________________________
bool MuonPlotter::passesZVeto(float lower, float upper){
// Checks if any combination of opposite sign, same flavor tight leptons (e or mu)
// has invariant mass closer than dm to the Z mass, returns true if none found
	const float MMU = 0.1057;
	const float MEL = 0.0005;
	const float MZ  = 91.2;

	if(NMus > 1){
		// First muon
		for(size_t i = 0; i < NMus-1; ++i){
			if(isTightMuon(i)){
				TLorentzVector pmu1, pmu2;
				pmu1.SetPtEtaPhiM(MuPt[i], MuEta[i], MuPhi[i], MMU);

				// Second muon
				for(size_t j = i+1; j < NMus; ++j){ 
					if(isTightMuon(j) && (MuCharge[i] != MuCharge[j]) ){
						pmu2.SetPtEtaPhiM(MuPt[j], MuEta[j], MuPhi[j], MMU);
						float mass = (pmu1+pmu2).M();
						if(mass < lower || mass > upper) return false;
					}
				}
			}
		}
	}

	if(NEls > 1){
		// First electron
		for(size_t i = 0; i < NEls-1; ++i){
			if(isTightElectron(i)){
				TLorentzVector pel1, pel2;
				pel1.SetPtEtaPhiM(ElPt[i], ElEta[i], ElPhi[i], MEL);

				// Second electron
				for(size_t j = i+1; j < NEls; ++j){
					if(isTightElectron(j) && (ElCh[i] != ElCh[j]) ){
						pel2.SetPtEtaPhiM(ElPt[j], ElEta[j], ElPhi[j], MEL);
						float mass = (pel1+pel2).M();
						if(mass < lower || mass > upper) return false;
					}
				}
			}
		}
	}
	return true;
}

//____________________________________________________________________________
bool MuonPlotter::passesMllEventVeto(float cut){
// Checks if any combination of opposite sign, same flavor tight leptons (e or mu)
// has invariant mass smaller than cut, returns true if none found
	const float MMU = 0.1057;
	const float MEL = 0.0005;
	const float MZ  = 91.2;

	if(NMus > 1){
		// First muon
		for(size_t i = 0; i < NMus-1; ++i){
			if(isTightMuon(i)){
				TLorentzVector pmu1, pmu2;
				pmu1.SetPtEtaPhiM(MuPt[i], MuEta[i], MuPhi[i], MMU);

				// Second muon
				for(size_t j = i+1; j < NMus; ++j){ 
					if(isTightMuon(j) && (MuCharge[i] != MuCharge[j]) ){
						pmu2.SetPtEtaPhiM(MuPt[j], MuEta[j], MuPhi[j], MMU);
						if((pmu1+pmu2).M() < cut) return false;
					}
				}
			}
		}		
	}

	if(NEls > 1){
		// First electron
		for(size_t i = 0; i < NEls-1; ++i){
			if(isTightElectron(i)){
				TLorentzVector pel1, pel2;
				pel1.SetPtEtaPhiM(ElPt[i], ElEta[i], ElPhi[i], MEL);

				// Second electron
				for(size_t j = i+1; j < NEls; ++j){
					if(isTightElectron(j) && (ElCh[i] != ElCh[j]) ){
						pel2.SetPtEtaPhiM(ElPt[j], ElEta[j], ElPhi[j], MEL);
						if((pel1+pel2).M() < cut) return false;
					}
				}
			}
		}
	}
	return true;
}

//____________________________________________________________________________
bool MuonPlotter::isMuTriggeredEvent(){
	if(HLT_Mu9 == 0 &&
	   HLT_Mu11 == 0 &&
	   HLT_Mu13_v1 == 0 &&
	   HLT_Mu15 == 0 &&
	   HLT_Mu15_v1 == 0 &&
	   HLT_DoubleMu0 == 0 &&
	   HLT_DoubleMu3 == 0 &&
	   HLT_DoubleMu3_v2 == 0 &&
	   HLT_DoubleMu5_v2 == 0
	   ) return false;
	return true;
}

//____________________________________________________________________________
bool MuonPlotter::isElTriggeredEvent(){
	// Leptonic triggers from UCSD/UCSB/FNAL list
	if (Run==1)                     return  HLT_Ele10_LW_L1R;
	if (Run>1       && Run<138000)  return (HLT_Ele10_LW_L1R            || HLT_Ele10_SW_L1R           || HLT_Ele15_LW_L1R           || HLT_DoubleEle5_SW_L1R);
	if (Run>=138000 && Run<=141900) return (HLT_Ele15_LW_L1R            || HLT_Ele15_SW_L1R           || HLT_Ele10_LW_EleId_L1R     || HLT_DoubleEle5_SW_L1R);
	if (Run>141900)                 return (HLT_Ele10_SW_EleId_L1R      || HLT_Ele15_SW_CaloEleId_L1R || HLT_Ele15_SW_EleId_L1R ||
	                                        HLT_Ele17_SW_LooseEleId_L1R || HLT_Ele17_SW_CaloEleId_L1R || HLT_Ele17_SW_EleId_L1R || 
	                                        HLT_Ele17_SW_TightEleId_L1R || HLT_Ele17_SW_TighterEleId_L1R_v1 || HLT_Ele20_SW_L1R ||
	                                        HLT_Ele22_SW_TighterEleId_L1R_v2        || HLT_Ele22_SW_TighterEleId_L1R_v3 || HLT_Ele27_SW_TightCaloEleIdTrack_L1R_v1 ||
	                                        HLT_Ele32_SW_TightCaloEleIdTrack_L1R_v1 || HLT_Ele32_SW_TighterEleId_L1R_v2 ||
	                                        HLT_DoubleEle10_SW_L1R      || HLT_DoubleEle15_SW_L1R_v1  || HLT_DoubleEle17_SW_L1R_v1);
	return false;
}

//____________________________________________________________________________
bool MuonPlotter::isJetTriggeredEvent(){
	if(HLT_Jet15U == 0  && 
	   HLT_Jet30U == 0  && 
	   HLT_Jet50U == 0  &&
	   HLT_Jet70U == 0  &&
	   HLT_Jet100U == 0 &&
	   HLT_Jet100U_v2 == 0 &&
	   HLT_Jet100U_v3 == 0
	) return false;
	return true;
}

//____________________________________________________________________________
bool MuonPlotter::isHTTriggeredEvent(){
	if(HLT_HT100U == 0  &&
	   HLT_HT120U == 0  &&
	   HLT_HT130U == 0  &&
	   HLT_HT140U == 0  &&
	   HLT_HT150U == 0  &&
	   HLT_HT150U_v3 == 0 &&
	   HLT_HT200U == 0
	) return false;
	return true;
}

//____________________________________________________________________________
bool MuonPlotter::isSigSupMuEvent(){
	if(isGoodMuEvent() == false) return false;
	if(MuMT > 20.) return false;
	if(pfMET > 20.) return false;
	if(NMus > 1) return false;
	return true;
}

//____________________________________________________________________________
bool MuonPlotter::isSigSupMuEventTRG(){
	if(!isJetTriggeredEvent() && !isHTTriggeredEvent()) return false;
	if(!isSigSupMuEvent()) return false;
	return true;
}

//____________________________________________________________________________
bool MuonPlotter::isZMuMuEvent(){
	if(isGoodMuEvent() == false) return false;
	if(NMus != 2) return false;
	if(isLooseMuon(0) == false || isLooseMuon(1) == false) return false;
	if(MuCharge[0] == MuCharge[1]) return false;

	// Z mass window cut
	TLorentzVector p1, p2;
	p1.SetPtEtaPhiM(MuPt[0], MuEta[0], MuPhi[0], 0.1057);
	p2.SetPtEtaPhiM(MuPt[1], MuEta[1], MuPhi[1], 0.1057);
	double m = (p1+p2).M();
	if(fabs(91.2 - m) > 15.) return false;

	if(pfMET > 20.) return false;

	return true;
}

//____________________________________________________________________________
bool MuonPlotter::isZMuMuEventTRG(){
	if(isMuTriggeredEvent() == false) return false;
	if(isZMuMuEvent() == false) return false;
	return true;
}

//____________________________________________________________________________
bool MuonPlotter::isZElElEvent(){
	if(isGoodElEvent() == false) return false;
	if(NEls != 2) return false;
	if(isLooseElectron(0) == false || isLooseElectron(1) == false) return false;
	if(ElCh[0] == ElCh[1]) return false;

	// Z mass window cut
	TLorentzVector p1, p2;
	p1.SetPtEtaPhiM(ElPt[0], ElEta[0], ElPhi[0], 0.0005);
	p2.SetPtEtaPhiM(ElPt[1], ElEta[1], ElPhi[1], 0.0005);
	double m = (p1+p2).M();
	if(fabs(91.2 - m) > 15.) return false;

	if(pfMET > 20.) return false;

	return true;
}

//____________________________________________________________________________
bool MuonPlotter::isZElElEventTRG(){
	if(isElTriggeredEvent() == false) return false;
	if(isZElElEvent() == false) return false;
	return true;
}

//____________________________________________________________________________
bool MuonPlotter::isSigSupElEvent(){
	if(!isGoodElEvent()) return false;
	if(ElMT[0] > 20.) return false;
	if(pfMET > 20.)   return false;
	if(NEls > 1)      return false;
	return true;
}

//____________________________________________________________________________
bool MuonPlotter::isSigSupElEventTRG(){
	if( HLT_Ele10_LW_L1R == 0 &&
	    HLT_Ele10_SW_L1R == 0 &&
	    HLT_Ele15_LW_L1R == 0 &&
		HLT_Ele15_SW_L1R == 0 ) return false;
	if(!isSigSupMuEvent()) return false;
	return true;
}

//____________________________________________________________________________
bool MuonPlotter::isGenMatchedSUSYDiLepEvent(){
	if(isGoodMuEvent() == false) return false;
	// if(isMuTriggeredEvent() == false) return false;
	if(!isSSTTMuEvent()) return false;
	if(isPromptSUSYMuon(0) && isPromptSUSYMuon(1)){
		if(isTightMuon(0) == 1 && isTightMuon(1) == 1) return true;
	}
	return false;
}

//____________________________________________________________________________
bool MuonPlotter::isSSLLMuEvent(){
	// This should include all the cuts for the final selection
	if(!isGoodMuEvent()) return false; // >1 jets, >0 loose muons
	if(NMus < 2) return false;         // >1 muons

	if(gSWITCH == 0) if(!passesZVeto(76., 106.)) return false; // no Zs in event
	if(gSWITCH == 1) if(!passesZVeto(15.))       return false;
	
	if(gSWITCH == 0) if(!passesMllEventVeto(12.)) return false; // no low mass OSSF pairs
	if(gSWITCH == 1) if(!passesMllEventVeto(5.) ) return false;
	
	if(gSWITCH == 0) if(!passesHTCut(60.) )  return false;    // ht cut
	if(gSWITCH == 1) if(!passesHTCut(300.))  return false;
	
	if(gSWITCH == 0) if(!passesMETCut(30.) ) return false;    // met cut
	if(gSWITCH == 1) if(!passesMETCut(30.))  return false;

	if(MuCharge[0] != MuCharge[1]) return false;              // SS

	if(!isGoodPrimMuon(0) || !isGoodSecMuon(1)) return false; // pt cuts
	return true;
}

//____________________________________________________________________________
bool MuonPlotter::isSSLLMuEventTRG(){
	if(gSWITCH == 0) if(!isMuTriggeredEvent()) return false;
	if(gSWITCH == 1) if(!isHTTriggeredEvent()) return false;
	if(!isSSLLMuEvent()) return false;
	return true;
}

//____________________________________________________________________________
bool MuonPlotter::isSSTTMuEvent(){
	if(!isSSLLMuEvent()) return false;
	if(!isTightMuon(0) || !isTightMuon(1)) return false;
	return true;
}

//____________________________________________________________________________
bool MuonPlotter::isSSTTMuEventTRG(){
	if(!isSSLLMuEventTRG()) return false;
	if(!isTightMuon(0) || !isTightMuon(1)) return false;
	return true;
}

//____________________________________________________________________________
bool MuonPlotter::isSSLLElEvent(){
	// This should include all the cuts for the final selection
	if(!isGoodElEvent()) return false; // >1 jets, >0 loose eles
	if(NEls < 2) return false;         // >1 eles

	if(gSWITCH == 0) if(!passesZVeto(76., 106.)) return false; // no Zs in event
	if(gSWITCH == 1) if(!passesZVeto(15.))       return false;

	if(gSWITCH == 0) if(!passesMllEventVeto(12.)) return false; // no low mass OSSF pairs
	if(gSWITCH == 1) if(!passesMllEventVeto(5.) ) return false;

	if(gSWITCH == 0) if(!passesHTCut(60.) )  return false;    // ht cut
	if(gSWITCH == 1) if(!passesHTCut(300.))  return false;

	if(gSWITCH == 0) if(!passesMETCut(30.) ) return false;    // met cut
	if(gSWITCH == 1) if(!passesMETCut(30.))  return false;

	if(ElCh[0] != ElCh[1]) return false;              // SS
	if(!isGoodPrimElectron(0) || !isGoodSecElectron(1)) return false; // pt cuts
	return true;
}

//____________________________________________________________________________
bool MuonPlotter::isSSLLElEventTRG(){
	if(gSWITCH == 0) if(!isElTriggeredEvent()) return false;
	if(gSWITCH == 1) if(!isHTTriggeredEvent()) return false;
	if(!isSSLLElEvent()) return false;
	return true;
}

//____________________________________________________________________________
bool MuonPlotter::isSSTTElEvent(){
	if(!isSSLLElEvent()) return false;
	if(!isTightElectron(0) || !isTightElectron(1)) return false;
	return true;
}

//____________________________________________________________________________
bool MuonPlotter::isSSTTElEventTRG(){
	if(!isSSLLElEventTRG()) return false;
	if(!isTightElectron(0) || !isTightElectron(1)) return false;
	return true;
}

//____________________________________________________________________________
bool MuonPlotter::isSSLLElMuEvent(){
	// This should include all the cuts for the final selection
	if(!isGoodElMuEvent()) return false; // >1 jets, >0 loose eles

	if(gSWITCH == 0) if(!passesZVeto(76., 106.)) return false; // no Zs in event
	if(gSWITCH == 1) if(!passesZVeto(15.))       return false;

	if(gSWITCH == 0) if(!passesMllEventVeto(12.)) return false; // no low mass OSSF pairs
	if(gSWITCH == 1) if(!passesMllEventVeto(5.) ) return false;

	if(gSWITCH == 0) if(!passesHTCut(60.) )  return false;    // ht cut
	if(gSWITCH == 1) if(!passesHTCut(300.))  return false;

	if(gSWITCH == 0) if(!passesMETCut(20.) ) return false;    // met cut
	if(gSWITCH == 1) if(!passesMETCut(30.))  return false;

	if(ElCh[0] != MuCharge[0]) return false;              // SS
	if(!isGoodPrimElectron(0) || !isGoodSecMuon(0) && 
	   !isGoodSecElectron(0)  || !isGoodPrimMuon(0)) return false; // pt cuts
	return true;
}

//____________________________________________________________________________
bool MuonPlotter::isSSLLElMuEventTRG(){
	// For UCSD/SB/FNAL just use OR of all lepton triggers
	if(gSWITCH == 0) if(!isElTriggeredEvent() && !isMuTriggeredEvent()) return false;
	if(gSWITCH == 1) if(!isHTTriggeredEvent()) return false;
	if(!isSSLLElMuEvent()) return false;
	return true;
}

//____________________________________________________________________________
bool MuonPlotter::isSSTTElMuEvent(){
	if(!isSSLLElMuEvent()) return false;
	if(!isTightElectron(0) || !isTightMuon(0)) return false;
	return true;
}

//____________________________________________________________________________
bool MuonPlotter::isSSTTElMuEventTRG(){
	if(!isSSLLElMuEventTRG()) return false;
	if(!isTightElectron(0) || !isTightMuon(0)) return false;
	return true;
}

//////////////////////////////////////////////////////////////////////////////
// Object selections:
//////////////////////////////////////////////////////////////////////////////
// Muons
//____________________________________________________________________________
bool MuonPlotter::isGoodMuon(int muon){
	float ptcut(5.);
	if(gSWITCH == 0) ptcut = 10.;
	if(gSWITCH == 1) ptcut = 5.;
	if(MuPt[muon] < ptcut) return false;
	return true;
}

//____________________________________________________________________________
bool MuonPlotter::isLooseMuon(int muon){
	if(isGoodMuon(muon) == false)  return false;
	if(gSWITCH == 0) if(MuIsoHybrid[muon] > 1.00) return false;
	if(gSWITCH == 1) if(MuIso[muon]       > 1.00) return false;
	return true;
}

//____________________________________________________________________________
bool MuonPlotter::isTightMuon(int muon){
	if(isGoodMuon(muon) == false)  return false;
	if(isLooseMuon(muon) == false) return false;
	if(gSWITCH == 0) if(MuIsoHybrid[muon] > 0.10) return false;
	if(gSWITCH == 1) if(MuIso[muon]       > 0.15) return false;
	return true;
}

//____________________________________________________________________________
bool MuonPlotter::isGoodPrimMuon(int muon){
	if(isLooseMuon(muon) == false) return false;
	if(gSWITCH == 0) if(MuPt[muon] < 20.) return false;
	if(gSWITCH == 1) if(MuPt[muon] < 5. ) return false;
	return true;
}

//____________________________________________________________________________
bool MuonPlotter::isGoodSecMuon(int muon){
	if(isLooseMuon(muon) == false) return false;
	if(gSWITCH == 0) if(MuPt[muon] < 10.) return false;
	if(gSWITCH == 1) if(MuPt[muon] < 5. ) return false;
	return true;
}

//____________________________________________________________________________
bool MuonPlotter::isFakeTTbarMuon(int muon){
	if(isLooseMuon(muon) == false) return false;
	if(abs(MuGenMoID[muon]) == 24 || abs(MuGenMoID[muon]) == 15) return false;
	return true;
}

//____________________________________________________________________________
bool MuonPlotter::isPromptTTbarMuon(int muon){
	if(isLooseMuon(muon) == false) return false;
	if(abs(MuGenMoID[muon] == 24 && abs(MuGenGMoID[muon]) == 6))  return true;
	if(abs(MuGenMoID[muon] == 15 && abs(MuGenGMoID[muon]) == 24)) return true;
	return false;
}

//____________________________________________________________________________
bool MuonPlotter::isPromptSUSYMuon(int muon){
	if(isLooseMuon(muon) == false) return false;
	if( abs(MuGenMoType[muon]) == 9 || abs(MuGenMoType[muon]) == 4  || abs(MuGenMoType[muon]) == 2 ) return true;
	return false;
}

//////////////////////////////////////////////////////////////////////////////
// Electrons
//____________________________________________________________________________
bool MuonPlotter::isGoodElectron(int ele){
	return true;
}

//____________________________________________________________________________
bool MuonPlotter::isLooseElectron(int ele){
	// All electrons are already loose in the high-pt selection (hybiso)
	if(gSWITCH == 1){
		if( fabs(ElEta[ele]) < 1.479 ) if(ElRelIso[ele] > 1.00) return false;
		else                           if(ElRelIso[ele] > 0.60) return false;		
	}
	return true;
}

//____________________________________________________________________________
bool MuonPlotter::isTightElectron(int ele){
	if(!isLooseElectron(ele)) return false;
	if(ElIsGoodElId_WP80[ele]    != 1) return false;
	if(ElIsConvertedEl_WP80[ele] == 1) return false;
	if(ElChIsCons[ele] != 1) return false;
	
	if(gSWITCH == 0) if(ElHybRelIso[ele] > 0.10) return false;
	if(gSWITCH == 1) if(ElRelIso[ele]    > 0.15) return false;
	
	return true;
}

//____________________________________________________________________________
bool MuonPlotter::isGoodPrimElectron(int ele){
	if(isLooseElectron(ele) == false) return false;
	if(gSWITCH == 0) if(ElPt[ele] < 20.) return false;
	if(gSWITCH == 1) if(ElPt[ele] < 10. ) return false;
	return true;
}

//____________________________________________________________________________
bool MuonPlotter::isGoodSecElectron(int ele){
	if(isLooseElectron(ele) == false) return false;
	if(gSWITCH == 0) if(ElPt[ele] < 10.) return false;
	if(gSWITCH == 1) if(ElPt[ele] < 10. ) return false;
	return true;
}

//____________________________________________________________________________
bool MuonPlotter::isGoodJet(int jet){
	float minDR = 0.4;
	for(size_t imu = 0; imu < NMus; ++imu){
		if(!isTightMuon(imu)) continue;
		if(Util::GetDeltaR(MuEta[imu], JetEta[jet], MuPhi[imu], JetPhi[jet]) > minDR ) continue;
		return false;
	}
	for(size_t iel = 0; iel < NEls; ++iel){
		if(!isTightElectron(iel)) continue;
		if(Util::GetDeltaR(ElEta[iel], JetEta[jet], ElPhi[iel], JetPhi[jet]) > minDR ) continue;
		return false;
	}
	return true;
}
