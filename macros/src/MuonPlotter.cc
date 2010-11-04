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

enum gSamples{MuA, MuB, EGA, EGB, JMA, JMB, TTbar, WJets, ZJets, VVJets, QCD15, QCD30, QCD80, QCD170, LM0};

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
MuonPlotter::~MuonPlotter(){
	fOutputFile->Close();
	delete fOutputFile;
}

//____________________________________________________________________________
void MuonPlotter::init(TString filename){
	if(fVerbose > 0) cout << "------------------------------------" << endl;
	if(fVerbose > 0) cout << "Initializing MuonPlotter ... " << endl;
	Util::SetStyle();
	loadSamples(filename);
	readVarNames("anavarnames.dat");

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
			sscanf(buffer, "Color\t%f", &ParValue);
			s.color = ParValue;

			if(fVerbose > 0){
				cout << " ---- " << endl;
				cout << "  New sample added: " << s.name << endl;
				cout << "   Sample no.  " << counter << endl;
				cout << "   Short name: " << s.sname << endl;
				cout << "   File:       " << (s.file)->GetName() << endl;
				cout << "   Events:     " << s.tree->GetEntries() << endl;
				cout << "   Lumi:       " << s.lumi << endl;
				cout << "   Color:      " << s.color << endl;
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
	// cout << " TTbar sample is number: " << fSampleMap["TTbar"] << endl;
	// cout << " Name of sample TTbar is " << fSamples[TTbar].sname << endl;

	makePrediction();

	// printYields(fMCBGSig, false);
	// printYields(fMuData, true);

	// makefRatioPlots();
	// makepRatioPlots();
	// make100pbMCPredictionPlots();
}

//____________________________________________________________________________
void MuonPlotter::make100pbMCPredictionPlots(){
	// Fill the ratios
	fLumiNorm = 1000.;

	cout << "Producing prediction for :" << endl;
	for(size_t i = 0; i < fMCBGSig.size(); ++i){
		int ind = fMCBGSig[i];
		cout << " " << fSamples[ind].sname << flush;
	}
	cout << endl;

	fillfRatio(fMCBGSig, 0);
	fillpRatio(fMCBGSig, 0);

	makeSSPredictionPlots(fMCBGSig);
}

//____________________________________________________________________________
void MuonPlotter::makefRatioPlots(){
	TH1D *h_fdata  = fillRatioPt(fMuData, 0, &MuonPlotter::isSignalSuppressedEventTRG, &MuonPlotter::isLooseMuon);      // JetMET Dataset (Single Muon Selection)
	TH1D *h_fttbar = fillRatioPt(TTbar,   0, &MuonPlotter::isGoodEvent,                &MuonPlotter::isFakeTTbarMuon); // TTbarJets MC
	TH1D *h_fallmc = fillRatioPt(fMCBG,   0, &MuonPlotter::isSignalSuppressedEvent,    &MuonPlotter::isLooseMuon);      // QCD MC
	h_fdata ->SetName("fRatioData");
	h_fttbar->SetName("fRatioTTbar");
	h_fallmc->SetName("fRatioAllMC");
	
	setPlottingRange(h_fdata, h_fttbar, h_fallmc);

	h_fdata ->SetMinimum(0.);
	h_fttbar->SetMinimum(0.);
	h_fallmc->SetMinimum(0.);
	
	h_fdata ->SetMarkerColor(kBlack);
	h_fttbar->SetMarkerColor(kBlue);
	h_fallmc->SetMarkerColor(kRed);

	h_fdata ->SetMarkerStyle(20);
	h_fttbar->SetMarkerStyle(20);
	h_fallmc->SetMarkerStyle(20);

	// h_fdata ->SetMarkerSize(2);
	// h_fttbar->SetMarkerSize(2);
	// h_fallmc->SetMarkerSize(2);
	
	h_fdata ->SetLineWidth(2);
	h_fttbar->SetLineWidth(2);
	h_fallmc->SetLineWidth(2);
	
	h_fdata ->SetLineColor(kBlack);
	h_fttbar->SetLineColor(kBlue);
	h_fallmc->SetLineColor(kRed);
	
	h_fdata ->SetFillColor(kBlack);
	h_fttbar->SetFillColor(kBlue);
	h_fallmc->SetFillColor(kRed);

	plotRatioOverlay3H(h_fdata, "Data (Jet, L = 21.7 pb^{-1})", h_fttbar, "t#bar{t} Fake GenMatch", h_fallmc, "QCD, t#bar{t}+jets, V+jets");
	// plotRatioOverlay3H(h_fdata, "Data (Mu, L = 21.3 pb^{-1})", h_fttbar, "t#bar{t} Fake GenMatch", h_fallmc, "QCD, t#bar{t}+jets, V+jets");
}

//____________________________________________________________________________
void MuonPlotter::makepRatioPlots(){
	TH1D *h_pdata  = fillRatioPt(fMuData, 0, &MuonPlotter::isZEventTRG, &MuonPlotter::isLooseMuon); // Mu Dataset (Di Muon Selection)
	TH1D *h_pttbar = fillRatioPt(TTbar,   0, &MuonPlotter::isGoodEvent, &MuonPlotter::isPromptTTbarMuon); // TTbar
	TH1D *h_pallmc = fillRatioPt(fMCBG,   0, &MuonPlotter::isZEvent,    &MuonPlotter::isLooseMuon); // all MC
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

	plotRatioOverlay2H(h_pttbar, "t#bar{t} Prompt GenMatch", h_pallmc, "QCD, t#bar{t}+jets, V+jets");
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
	const int nbins = 20;
	TH2D *h2_sig = new TH2D("h2_sig", "Isolation vs Pt for muons in data",       nbins, 0., 1., gNPtbins, gPtbins);
	TH2D *h2_bg  = new TH2D("h2_bg",  "Isolation vs Pt for fake muons in ttbar", nbins, 0., 1., gNPtbins, gPtbins);
	h2_sig->SetXTitle(convertVarName("MuIso[0]"));
	h2_bg ->SetXTitle(convertVarName("MuIso[0]"));
	h2_sig->SetYTitle(convertVarName("MuPt[0]"));
	h2_bg ->SetYTitle(convertVarName("MuPt[0]"));
		
	TTree *tree = fSamples[sample1].tree;
	float scale = fLumiNorm / fSamples[sample1].lumi;
	scale = 1.;
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

	tree = fSamples[sample2].tree;
	scale = fLumiNorm / fSamples[sample2].lumi;
	scale = 1.;
	tree->ResetBranchAddresses();
	Init(tree);
	if (fChain == 0) return;
	nentries = fChain->GetEntriesFast();
	nbytes = 0, nb = 0;
	for (Long64_t jentry=0; jentry<nentries;jentry++) {
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		nb = fChain->GetEntry(jentry);   nbytes += nb;
		
		if((*this.*eventSelector2)() == false) continue;
		if((*this.*muonSelector2)(muon2) == false) continue;
		h2_bg->Fill(MuIso[muon2], MuPt[muon2], scale);
	}

	// fSamples[sample1].tree->Project("h2_bg", Form("MuPt[%d]:MuIso[%d]", muon1, muon1), c1);
	// fSamples[sample2].tree->Project("h2_sig", Form("MuPt[%d]:MuIso[%d]", muon2, muon2), c2);

	TLatex *lat = new TLatex();
	lat->SetNDC(kTRUE);
	lat->SetTextColor(kBlack);
	lat->SetTextSize(0.06);
	
	TCanvas *c_temp = new TCanvas("IsoVsPt", "Isolating in Pt bins", 0, 0, 1200, 800);
	c_temp->Divide(3,2);
	c_temp->cd(6);
	if(logy) gPad->SetLogz(1);
	h2_sig->DrawCopy("colz");
	lat->DrawLatex(0.11,0.92, fSamples[sample1].sname);

	for(size_t i = 1; i <= 5; ++i){
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

//____________________________________________________________________________
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

		// TH2D *H_ntight_temp = new TH2D(Form("NTight-%s",fSamples[sample].sname.Data()), "NTight Muons", h_2d->GetNbinsX(), h_2d->GetXaxis()->GetXbins()->GetArray(), h_2d->GetNbinsY(),  h_2d->GetYaxis()->GetXbins()->GetArray());
		// TH2D *H_nloose_temp = new TH2D(Form("NLoose-%s",fSamples[sample].sname.Data()), "NLoose Muons", h_2d->GetNbinsX(), h_2d->GetXaxis()->GetXbins()->GetArray(), h_2d->GetNbinsY(),  h_2d->GetYaxis()->GetXbins()->GetArray());
		// H_ntight_temp->Sumw2();
		// H_nloose_temp->Sumw2();

		float scale = fLumiNorm / fSamples[sample].lumi;

		// nentries = (Long64_t)nentries*scale;

		Long64_t nbytes = 0, nb = 0;
		for (Long64_t jentry=0; jentry<nentries;jentry++) {
			Long64_t ientry = LoadTree(jentry);
			if (ientry < 0) break;
			nb = fChain->GetEntry(jentry);   nbytes += nb;

			if((*this.*eventSelector)() == false) continue;
			if((*this.*muonSelector)(muon) == false) continue;
			
			// H_nloose_temp->Fill(MuPt[muon], MuEta[muon]); // Tight or loose
			// if(MuTight[muon]) H_ntight_temp->Fill(MuPt[muon], MuEta[muon]); // Tight
    		
			if(isLooseMuon(muon)) H_nloose->Fill(MuPt[muon], MuEta[muon], scale); // Tight or loose
			if(isTightMuon(muon)) H_ntight->Fill(MuPt[muon], MuEta[muon], scale); // Tight
		}
	
		// // Rescale errors:
		// H_ntight_temp->SetBinError(1,1, H_ntight_temp->GetBinError(1,1) * 1.0 / sqrt(scale));
		// H_nloose_temp->SetBinError(1,1, H_nloose_temp->GetBinError(1,1) * 1.0 / sqrt(scale));
		// 
		// H_ntight->Add(H_ntight_temp, scale);
		// H_nloose->Add(H_nloose_temp, scale);
	
		if(fVerbose>2) cout << " Tight entries so far: " << H_ntight->GetEntries() << endl;
		if(fVerbose>2) cout << " Loose entries so far: " << H_nloose->GetEntries() << endl;
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

	// h_pt ->SetMinimum(0);
	// h_eta->SetMinimum(0);
	// h_pt ->SetMaximum(0.3);
	// h_eta->SetMaximum(0.3);
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

	produceRatio(samples, muon, &MuonPlotter::isSignalSuppressedEvent, &MuonPlotter::isLooseMuon, fH2D_fRatio, fH1D_fRatioPt, fH1D_fRatioEta, true);
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

	produceRatio(samples, muon, &MuonPlotter::isZEvent, &MuonPlotter::isLooseMuon, fH2D_pRatio, fH1D_pRatioPt, fH1D_pRatioEta, true);
}

//____________________________________________________________________________
void MuonPlotter::plotRatio(int sample, int muon, bool(MuonPlotter::*eventSelector)(), bool(MuonPlotter::*muonSelector)(int), TString tag){
	vector<int> samples; samples.push_back(sample);
	plotRatio(samples, muon, eventSelector, muonSelector, tag);
}

//____________________________________________________________________________
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

//____________________________________________________________________________
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

//____________________________________________________________________________
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
void MuonPlotter::makeSSNsigPredictionPlots(){
	// Fake Ratios: ////////////////////////////////////////////////////////////////////
	fillfRatio(0, 0); // QCD-Pt20
	fillpRatio(7, 0); // Zjets

	// Prediction: /////////////////////////////////////////////////////////////////////
	TH1D *H_nsigpred = new TH1D("Nsigpred", "Predicted N_sig in Pt1 bins", gNPtbins,  gPtbins);
	bool output = true;
	H_nsigpred->Add(NsigPredFromFPRatios(2, output)[0], fLumiNorm/fSamples[2].lumi); // LM0
	H_nsigpred->Add(NsigPredFromFPRatios(4, output)[0], fLumiNorm/fSamples[4].lumi); // ttbar

	// Observation: ////////////////////////////////////////////////////////////////////
	TH1D *H_nsigobs = new TH1D("Nsigobs", "Observed N_sig in Pt1 bins",  gNPtbins,  gPtbins);
	vector<int> samples;
	samples.push_back(2);
	NObs(H_nsigobs, samples, &MuonPlotter::isGenMatchedSUSYDiLepEvent);

	// Output
	H_nsigobs->SetXTitle(convertVarName("MuPt[1]"));
	H_nsigobs->SetYTitle(Form("Events / %2.0f GeV", fBinWidthScale));

	// Normalize to binwidth
	H_nsigobs  = normHistBW(H_nsigobs,  fBinWidthScale);
	H_nsigpred = normHistBW(H_nsigpred, fBinWidthScale);
	
	plotPredOverlay2HWithRatio(H_nsigobs, "Observed", H_nsigpred, "Predicted");
	plotPredOverlay2HWithRatio(H_nsigobs, "Observed", H_nsigpred, "Predicted", true);
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
	NObs(H_nsigobs, samples, &MuonPlotter::isGenMatchedSUSYDiLepEvent);

	TH1D *H_nt2obs  = new TH1D("Nt2obs", "Observed N_t2 in Pt1 bins",  gNPtbins,  gPtbins);
	NObs(H_nt2obs, samples, &MuonPlotter::isSSTTEvent);

	TH1D *H_nt2obsttbar = new TH1D("Nt2obsttbar", "Observed N_t2 in Pt1 bins, ttbar only",  gNPtbins,  gPtbins);
	vector<int> ttbarsample; ttbarsample.push_back(4);
	NObs(H_nt2obsttbar, ttbarsample, &MuonPlotter::isSSTTEvent);	

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
	setPlottingRange(hists);
	
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
	plotPredOverlay3HWithRatio(H_nFpred, "Predicted Fakes", H_nt2obs,  "Obs. N_{t2} (QCD, t#bar{t}+jets, V+jets, LM0)", H_nt2obsttbar, "Obs. N_{t2} (t#bar{t}+jets)", false, false);
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
void MuonPlotter::makePrediction(){
	// Memory:
	// > Further signal suppression for data?
	// > Use InclusiveMu15 instead of QCD?

	bool data = true; // Use ratios from data or mc?
	
	// Which samples to use for nt2/nt1/nt0 input?
	// vector<int> inputsamples = fMCBG;
	// vector<int> inputsamples = fMCBGSig;
	vector<int> inputsamples = fMuData;

	// Which luminosity to use?
	fLumiNorm = fSamples[MuA].lumi + fSamples[MuB].lumi;
	// fLumiNorm = 30.;
	// fLumiNorm = 100.;
	// fLumiNorm = 1000.;
	
	// Dummy binning
	const int nptbins = 1;
	const double ptbins[nptbins+1] = {10., 1000.};
	const int netabins = 1;
	const double etabins[netabins+1] = {-2.5, 2.5};

	// Fill the ratios
	TH2D *fratio_allmc = fillRatio(fMCBG, 0, &MuonPlotter::isSignalSuppressedEvent, &MuonPlotter::isLooseMuon, nptbins, ptbins, netabins, etabins);
	cout << "  fRatio from all MC         = " << fratio_allmc->GetBinContent(1,1) << " +/- " << fratio_allmc->GetBinError(1,1) << endl;

	TH2D *fratio_ttbar = fillRatio(TTbar,    0, &MuonPlotter::isGoodEvent,             &MuonPlotter::isFakeTTbarMuon, nptbins, ptbins, netabins, etabins);
	cout << "  fRatio from ttbar genmatch = " << fratio_ttbar->GetBinContent(1,1) << " +/- " << fratio_ttbar->GetBinError(1,1) << endl;

	TH2D *fratio_data = fillRatio(fMuData,    0, &MuonPlotter::isSignalSuppressedEvent, &MuonPlotter::isLooseMuon, nptbins, ptbins, netabins, etabins);
	cout << "  fRatio from data           = " << fratio_data->GetBinContent(1,1) << " +/- " << fratio_data->GetBinError(1,1) << endl;

	cout << " ------------------------------------" << endl;
	TH2D *pratio_allmc = fillRatio(fMCBG, 0,  &MuonPlotter::isZEvent, &MuonPlotter::isLooseMuon, nptbins, ptbins, netabins, etabins);
	cout << "  pRatio from all MC         = " << pratio_allmc->GetBinContent(1,1) << " +/- " << pratio_allmc->GetBinError(1,1) << endl;

	TH2D *pratio_ttbar = fillRatio(TTbar,    0,  &MuonPlotter::isGoodEvent,            &MuonPlotter::isPromptTTbarMuon, nptbins, ptbins, netabins, etabins);
	cout << "  pRatio from ttbar genmatch = " << pratio_ttbar->GetBinContent(1,1) << " +/- " << pratio_ttbar->GetBinError(1,1) << endl;

	TH2D *pratio_data = fillRatio(fMuData,    0,  &MuonPlotter::isZEvent, &MuonPlotter::isLooseMuon, nptbins, ptbins, netabins, etabins);
	cout << "  pRatio from data           = " << pratio_data->GetBinContent(1,1) << " +/- " << pratio_data->GetBinError(1,1) << endl;
	cout << " ------------------------------------" << endl;

	// Add systematics to ratios
	float ratioallmc(0.), ratiottbar(0.), deviation(0.), olderror(0.), newerror(0.);
	ratioallmc = fratio_allmc->GetBinContent(1,1);
	ratiottbar = fratio_ttbar->GetBinContent(1,1);
	deviation  = fabs(ratioallmc - ratiottbar);
	// Add to mc ratios
	olderror   = fratio_allmc->GetBinError(1,1);
	newerror   = sqrt(olderror*olderror + deviation*deviation);
	fratio_allmc->SetBinError(1,1,newerror);
	// Add to data ratios
	olderror = fratio_data->GetBinError(1,1);
	newerror   = sqrt(olderror*olderror + deviation*deviation);
	fratio_data->SetBinError(1,1,newerror);
	cout << "  new fRatio (all MC)        = " << fratio_allmc->GetBinContent(1,1) << " +/- " << fratio_allmc->GetBinError(1,1) << endl;
	cout << "  new fRatio (data)          = " << fratio_data->GetBinContent(1,1) << " +/- " << fratio_data->GetBinError(1,1) << endl;

	ratioallmc = pratio_allmc->GetBinContent(1,1);
	ratiottbar = pratio_ttbar->GetBinContent(1,1);
	deviation  = fabs(ratioallmc - ratiottbar);
	// Add to mc ratios
	olderror   = pratio_allmc->GetBinError(1,1);
	newerror   = sqrt(olderror*olderror + deviation*deviation);
	pratio_allmc->SetBinError(1,1,newerror);
	// Add to data ratios
	olderror   = pratio_data->GetBinError(1,1);
	newerror   = sqrt(olderror*olderror + deviation*deviation);
	pratio_data->SetBinError(1,1,newerror);
	cout << "  new pRatio (all MC)        = " << pratio_allmc->GetBinContent(1,1) << " +/- " << pratio_allmc->GetBinError(1,1) << endl;
	cout << "  new pRatio (data)          = " << pratio_data->GetBinContent(1,1) << " +/- " << pratio_data->GetBinError(1,1) << endl;
	cout << " ------------------------------------" << endl;

	// Get Nt0/Nt1/Nt2
	TH1D *h_nt0 = new TH1D("nt0", "nt0", 1, 0, 10.);
	TH1D *h_nt1 = new TH1D("nt1", "nt1", 1, 0, 10.);
	TH1D *h_nt2 = new TH1D("nt2", "nt2", 1, 0, 10.);
	for(size_t i = 0; i < inputsamples.size(); ++i){
		int index = inputsamples[i];
		float scale = fLumiNorm/fSamples[index].lumi; // Normalize all
		// Fill histograms from tree
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

			if(isSSLLEvent() == false) continue;

			if(  isTightMuon(0) &&  isTightMuon(1) ) h_nt2->Fill(1); // Tight-tight
			if(  isTightMuon(0) && !isTightMuon(1) ) h_nt1->Fill(1); // Tight-loose
			if( !isTightMuon(0) &&  isTightMuon(1) ) h_nt1->Fill(1); // Loose-tight
			if( !isTightMuon(0) && !isTightMuon(1) ) h_nt0->Fill(1); // Loose-loose
			// if(  isTightMuon(0) &&  isTightMuon(1) ) h_nt2->Fill(1, scale); // Tight-tight
			// if(  isTightMuon(0) && !isTightMuon(1) ) h_nt1->Fill(1, scale); // Tight-loose
			// if( !isTightMuon(0) &&  isTightMuon(1) ) h_nt1->Fill(1, scale); // Loose-tight
			// if( !isTightMuon(0) && !isTightMuon(1) ) h_nt0->Fill(1, scale); // Loose-loose
	   }
	}
	double nt2 = h_nt2->GetBinContent(1);
	double nt1 = h_nt1->GetBinContent(1);
	double nt0 = h_nt0->GetBinContent(1);
	cout << "  Found " << nt2 << " +/- " << h_nt2->GetBinError(1) << " events in Nt2" << endl;
	cout << "  Found " << nt1 << " +/- " << h_nt1->GetBinError(1) << " events in Nt1" << endl;
	cout << "  Found " << nt0 << " +/- " << h_nt0->GetBinError(1) << " events in Nt0" << endl;
	cout << " ------------------------------------" << endl;
	
	// Make prediction
	fFPRatios = new FPRatios();
	fFPRatios->SetVerbose(fVerbose);
	if(data){
		fFPRatios->SetMuFratios(fratio_data);
		fFPRatios->SetMuPratios(pratio_data);		
	}
	else{
		fFPRatios->SetMuFratios(fratio_allmc);
		fFPRatios->SetMuPratios(pratio_allmc);		
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
	cout << " ------------------------------------" << endl;
	
	// Get observation
	TH1D *H_nsigobs = new TH1D("Nsigobs", "Observed N_sig in Pt1 bins", nptbins, ptbins);
	NObs(H_nsigobs, inputsamples, &MuonPlotter::isGenMatchedSUSYDiLepEvent);

	TH1D *H_nt2obs = new TH1D("Nt2obsttbar", "Observed N_t2 in Pt1 bins, ttbar only",  nptbins,  ptbins);
	NObs(H_nt2obs, fMCBG, &MuonPlotter::isSSTTEvent);	

	cout << "  Observation from LM0:      " << H_nsigobs->GetBinContent(1) << " +/- " << H_nsigobs->GetBinError(1) << endl;
	cout << "  Nt2 observed from SM only: " << H_nt2obs->GetBinContent(1) << " +/- " << H_nt2obs->GetBinError(1) << endl;
	
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

	// Fill histograms from tree
	TTree *tree = fSamples[sample].tree;
	Init(tree);
   if (fChain == 0) return res;
   Long64_t nentries = fChain->GetEntriesFast();
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;

		if(isSSLLEvent() == false) continue;
		
		if(  MuTight[0] &&  MuTight[1] ) H_nt2mes->Fill(MuPt[0], MuPt[1]); // Tight-tight
		if(  MuTight[0] && !MuTight[1] ) H_nt1mes->Fill(MuPt[0], MuPt[1]); // Tight-loose
		if( !MuTight[0] &&  MuTight[1] ) H_nt1mes->Fill(MuPt[1], MuPt[0]); // Loose-tight
		if( !MuTight[0] && !MuTight[1] ) H_nt0mes->Fill(MuPt[0], MuPt[1]); // Loose-loose
   }
	
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
void MuonPlotter::printYields(vector<int> samples, bool trigger){
	for(size_t i = 0; i < samples.size(); ++i){
		int index = samples[i];

		// Counters
		int nt2(0), nt1(0), nt0(0);
		int ntsigsupp(0), nlsigsupp(0);
		int nzt(0), nzl(0);

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

			if(!trigger){
				if(isSSLLEvent()){
					if(  isTightMuon(0) &&  isTightMuon(1) ) nt2++; // Tight-tight
					if(  isTightMuon(0) && !isTightMuon(1) ) nt1++; // Tight-loose
					if( !isTightMuon(0) &&  isTightMuon(1) ) nt1++; // Loose-tight
					if( !isTightMuon(0) && !isTightMuon(1) ) nt0++; // Loose-loose				
				}
				if(isSignalSuppressedEvent()){
					if( isLooseMuon(0) ) nlsigsupp++;
					if( isTightMuon(0) ) ntsigsupp++;
				}
				if(isZEvent()){
					if( isLooseMuon(0) ) nzl++;
					if( isTightMuon(0) ) nzt++;
				}
			}
			if(trigger){
				if(isSSLLEventTRG()){
					if(  isTightMuon(0) &&  isTightMuon(1) ) nt2++; // Tight-tight
					if(  isTightMuon(0) && !isTightMuon(1) ) nt1++; // Tight-loose
					if( !isTightMuon(0) &&  isTightMuon(1) ) nt1++; // Loose-tight
					if( !isTightMuon(0) && !isTightMuon(1) ) nt0++; // Loose-loose				
				}
				if(isSignalSuppressedEventTRG()){
					if( isLooseMuon(0) ) nlsigsupp++;
					if( isTightMuon(0) ) ntsigsupp++;
				}
				if(isZEventTRG()){
					if( isLooseMuon(0) ) nzl++;
					if( isTightMuon(0) ) nzt++;
				}				
			}
		}
		
		cout << "-----------------------" << endl;
		cout << " Sample: " << fSamples[index].sname << endl;
		cout << "   Nt2 = "       << nt2 <<       "  Nt1 = "       << nt1 << "  Nt0 = " << nt0 << endl;
		cout << "   Nss tight = " << ntsigsupp << "  Nss loose = " << nlsigsupp << endl;
		cout << "   NZ tight  = " << nzt <<       "  NZ  loose = " << nzl << endl;
		cout << endl;
		float scale = (fSamples[MuA].lumi + fSamples[MuB].lumi) / fSamples[index].lumi;
		cout << "  scaled to " << fSamples[MuA].lumi + fSamples[MuB].lumi << " /pb" << endl;
		cout << "   Nt2 = "       << scale*nt2 <<       "  Nt1 = "       << scale*nt1 << "  Nt0 = " << scale*nt0 << endl;
		cout << "   Nss tight = " << scale*ntsigsupp << "  Nss loose = " << scale*nlsigsupp << endl;
		cout << "   NZ tight  = " << scale*nzt <<       "  NZ  loose = " << scale*nzl << endl;		
		cout << endl;
		scale = 100. / fSamples[index].lumi;
		cout << "  scaled to 100/pb" << " /pb" << endl;
		cout << "   Nt2 = "       << scale*nt2 <<       "  Nt1 = "       << scale*nt1 << "  Nt0 = " << scale*nt0 << endl;
		cout << "   Nss tight = " << scale*ntsigsupp << "  Nss loose = " << scale*nlsigsupp << endl;
		cout << "   NZ tight  = " << scale*nzt <<       "  NZ  loose = " << scale*nzl << endl;		
	}
	cout << "-----------------------" << endl;
}


//////////////////////////////////////////////////////////////////////////////
// Event Selections:
//____________________________________________________________________________
bool MuonPlotter::isGoodEvent(){
	// Some global cuts
	if(NMus < 1) return false;
	// if(NJets < 2) return false;
	if(!passesNJetCut(2)) return false;
	if(isLooseMuon(0) == false) return false;
	if(NMus > 1) if(isLooseMuon(1) == false) return false;
	return true;
}

//____________________________________________________________________________
bool MuonPlotter::passesNJetCut(int cut){
	int njets(0);
	for(size_t i = 0; i < NJets; ++i){
		bool good(true);
		for(size_t imu = 0; imu < NMus; ++imu){
			if(!isTightMuon(imu)) continue;
			if(Util::GetDeltaR(MuEta[imu], JetEta[i], MuPhi[imu], JetPhi[i]) > 0.4 ) continue;
			good = false;
		}
		for(size_t iel = 0; iel < NEls; ++iel){
			if(!isTightElectron(iel)) continue;
			if(Util::GetDeltaR(ElEta[iel], JetEta[i], ElPhi[iel], JetPhi[i]) > 0.4 ) continue;
			good = false;
		}
		if(good) njets++;
	}
	return njets>=cut;
}

//____________________________________________________________________________
bool MuonPlotter::isMuTriggeredEvent(){
	if(HLT_Mu9 == 0 &&
	   HLT_Mu11 == 0 &&
	   HLT_Mu15 == 0 &&
	   HLT_DoubleMu0 == 0 &&
	   HLT_DoubleMu3 == 0
	   ) return false;
	return true;
}

//____________________________________________________________________________
bool MuonPlotter::isJetTriggeredEvent(){
	if(HLT_Jet15U == 0  && 
	   HLT_Jet30U == 0  && 
	   HLT_Jet50U == 0  &&
	   HLT_Jet70U == 0  &&
	   HLT_Jet100U == 0 &&
	   HLT_HT100U == 0  &&
	   HLT_HT120U == 0  &&
	   HLT_HT140U == 0  &&
	   HLT_HT150U == 0
	   ) return false;
	return true;
}

//____________________________________________________________________________
bool MuonPlotter::isSignalSuppressedEvent(){
	if(isGoodEvent() == false) return false;
	if(MuMT > 20.) return false;
	if(pfMET > 20.) return false;
	if(NMus > 1) return false;
	return true;
}

//____________________________________________________________________________
bool MuonPlotter::isSignalSuppressedEventTRG(){
	// if(isMuTriggeredEvent() == false) return false;
	if(isJetTriggeredEvent() == false) return false;
	if(isSignalSuppressedEvent() == false) return false;
	return true;
}

//____________________________________________________________________________
bool MuonPlotter::isZEvent(){
	// if(isGoodEvent() == false) return false;
	if(NJets < 2) return false;
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
bool MuonPlotter::isZEventTRG(){
	if(isMuTriggeredEvent() == false) return false;
	if(isZEvent() == false) return false;
	return true;
}

//____________________________________________________________________________
bool MuonPlotter::isGenMatchedSUSYDiLepEvent(){
	if(isGoodEvent() == false) return false;
	// if(isMuTriggeredEvent() == false) return false;
	if(NMus < 2) return false;
	if(isGoodPrimMuon(0) == false || isGoodSecMuon(1) == false) return false;
	if(isPromptSUSYMuon(0) && isPromptSUSYMuon(1)){
		if(MuTight[0] == 1 && MuTight[1] == 1) return true;
	}
	return false;
}

//____________________________________________________________________________
bool MuonPlotter::isSSLLEvent(){
	if(isGoodEvent() == false) return false;
	if(NMus < 2) return false;
	if(MuCharge[0] != MuCharge[1]) return false;
	if(isGoodPrimMuon(0) == false || isGoodSecMuon(1) == false) return false;
	return true;
}

//____________________________________________________________________________
bool MuonPlotter::isSSLLEventTRG(){
	if(isMuTriggeredEvent() == false) return false;
	if(isSSLLEvent() == false) return false;
	return true;
}

//____________________________________________________________________________
bool MuonPlotter::isSSTTEvent(){
	if(isGoodEvent() == false) return false;
	// if(isMuTriggeredEvent() == false) return false;
	if(NMus < 2) return false;
	if(MuCharge[0] != MuCharge[1]) return false;
	if(isGoodPrimMuon(0) == false || isGoodSecMuon(1) == false) return false;
	if(MuTight[0] == 1 && MuTight[1] == 1){
		if(MuCharge[0] == MuCharge[1]) return true;
	}
	return false;
}

// //____________________________________________________________________________
// bool MuonPlotter::isAddZVetoEvent(){
// 	// Some global cuts
// 	if(NMus > 1) return false;
// 	if(NJets < 2) return false;
// 	if(isLooseMuon(0) == false) return false;
// 	if(NMus > 1) if(isLooseMuon(1) == false) return false;
// 	return true;
// }

//////////////////////////////////////////////////////////////////////////////
// Object selections:
//____________________________________________________________________________
bool MuonPlotter::isGoodMuon(int muon){
	if(MuPt[muon] < 10.) return false;
	if(fabs(MuEta[muon]) > 2.4) return false;
	if(fabs(MuD0BS[muon]) > 0.02) return false;
	return true;
}

//____________________________________________________________________________
bool MuonPlotter::isLooseMuon(int muon){
	if(isGoodMuon(muon) == false)  return false;
	// if(MuIso[muon] > 1.00) return false;
	if(MuIsoHybrid[muon] > 1.00) return false;
	return true;
}

//____________________________________________________________________________
bool MuonPlotter::isTightMuon(int muon){
	if(isGoodMuon(muon) == false)  return false;
	if(isLooseMuon(muon) == false) return false;
	// if(MuIso[muon] > 0.15) return false;
	if(MuIsoHybrid[muon] > 0.1) return false;
	return true;
}

//____________________________________________________________________________
bool MuonPlotter::isLooseNoTightMuon(int muon){
	if(isLooseMuon(muon) == true && isTightMuon(muon) == false) return true;
	return false;
}

//____________________________________________________________________________
bool MuonPlotter::isGoodPrimMuon(int muon){
	if(isLooseMuon(muon) == false) return false;
	if(MuPt[muon] < fMinPt1) return false;
	return true;
}

//____________________________________________________________________________
bool MuonPlotter::isGoodSecMuon(int muon){
	if(isLooseMuon(muon) == false) return false;
	if(MuPt[muon] < fMinPt2) return false;
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

//____________________________________________________________________________
bool MuonPlotter::isTightElectron(int ele){
	if(ElTight[ele] != 1) return false;
	return true;
}

