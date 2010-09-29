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

using namespace std;

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
	// loadFakeRatio(0,0); // QCD Sample, First muon
	readVarNames("anavarnames.dat");

	fNJetsMin = 2;     // Minimal number of jets cut
	fLumiNorm = 1000.; // Normalize everything to this lumi in /pb
	fBinWidthScale = 10.; // Normalize Y axis to this binwidth
	
	// Prevent root from adding histograms to current file
	TH1::AddDirectory(kFALSE);
}

//____________________________________________________________________________
MuonPlotter::ratio_histos MuonPlotter::getRatios(TFile *file){
	ratio_histos histos;
	histos.munt2           = (TH2D*)file->Get("MuNt2");
	histos.munt1           = (TH2D*)file->Get("MuNt1");
	histos.munt0           = (TH2D*)file->Get("MuNt0");
	histos.mu1tight        = (TH2D*)file->Get("Mu1Tight");
	histos.mu1loose        = (TH2D*)file->Get("Mu1Loose");
	histos.mu1loosenotight = (TH2D*)file->Get("Mu1LooseNoTight");
	histos.mu1ratio1       = (TH2D*)file->Get("Mu1Ratio1");
	histos.mu1ratio1pt     = (TH1D*)file->Get("Mu1Ratio1Pt");
	histos.mu1ratio1eta    = (TH1D*)file->Get("Mu1Ratio1Eta");
	histos.mu1ratio2       = (TH2D*)file->Get("Mu1Ratio2");
	histos.mu1ratio2pt     = (TH1D*)file->Get("Mu1Ratio2Pt");
	histos.mu1ratio2eta    = (TH1D*)file->Get("Mu1Ratio2Eta");
	histos.mu2tight        = (TH2D*)file->Get("Mu2Tight");
	histos.mu2loose        = (TH2D*)file->Get("Mu2Loose");
	histos.mu2loosenotight = (TH2D*)file->Get("Mu2LooseNoTight");
	histos.mu2ratio1       = (TH2D*)file->Get("Mu2Ratio1");
	histos.mu2ratio1pt     = (TH1D*)file->Get("Mu2Ratio1Pt");
	histos.mu2ratio1eta    = (TH1D*)file->Get("Mu2Ratio1Eta");
	histos.mu2ratio2       = (TH2D*)file->Get("Mu2Ratio2");
	histos.mu2ratio2pt     = (TH1D*)file->Get("Mu2Ratio2Pt");
	histos.mu2ratio2eta    = (TH1D*)file->Get("Mu2Ratio2Eta");
	return histos;
}

//____________________________________________________________________________
void MuonPlotter::loadSamples(const char* filename){
	char buffer[200];
	ifstream IN(filename);

	char ParName[100], StringValue[1000];
	float ParValue;

	// bool ok(false);

	if(fVerbose > 0) cout << "------------------------------------" << endl;
	if(fVerbose > 0) cout << "Sample File  " << filename << endl;

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
			s.tree = (TTree*)f->Get("MuonAnalysis");
			s.histos = getRatios(f);

			IN.getline(buffer, 200, '\n');
			sscanf(buffer, "Lumi\t%f", &ParValue);
			s.lumi = ParValue;

			IN.getline(buffer, 200, '\n');
			sscanf(buffer, "Color\t%f", &ParValue);
			s.color = ParValue;

			if(fVerbose > 0){
				cout << " ---- " << endl;
				cout << "  New sample added: " << s.name << endl;
				cout << "   Short name: " << s.sname << endl;
				cout << "   File:       " << (s.file)->GetName() << endl;
				cout << "   Events:     " << s.tree->GetEntries() << endl;
				cout << "   Lumi:       " << s.lumi << endl;
				cout << "   Color:      " << s.color << endl;
			}
			fSamples.push_back(s);
		}
	}
	if(fVerbose > 0) cout << "------------------------------------" << endl;
}

//____________________________________________________________________________
void MuonPlotter::makePlots(){
	// makeWJetsNt2PredFromTree();
	// makeWJetsNt2PredFromHistos();

	// // qcd vs ttbarf
	// makeIsoVsPtPlot(4, 0, "", 2, 1, "abs(MuGenMoID[1])!=24&&abs(MuGenMoID[1])!=15", "IsovsPt_QCD");
	// // zjets vs ttbarp
	// makeIsoVsPtPlot(5, 0, "", 2, 0, "(abs(MuGenMoID[0])==24||abs(MuGenMoID[0])==15)&&(abs(MuGenMoID[1])==24||abs(MuGenMoID[1])==15)", "IsovsPt_ZJets", true);
	// 
	// // qcd vs ttbarf
	// makeIsoVsNJetsPlot(4, 0, "", 2, 1, "abs(MuGenMoID[1])!=24&&abs(MuGenMoID[1])!=15", "IsovsNJets_QCD");
	// // zjets vs ttbarp
	// makeIsoVsNJetsPlot(5, 0, "", 2, 1, "abs(MuGenMoID[1])==24||abs(MuGenMoID[1])==15", "IsovsNJets_ZJets", true);

	// makeIsolationPlots();
	// makePtPlots();

	vector<int> samples;
	// samples.push_back(1); // QCD
	samples.push_back(2); // LM0
	// samples.push_back(3); // LM1
	samples.push_back(4); // TTbar
	// samples.push_back(5); // WJets
	// samples.push_back(6); // ZJets
	// samples.push_back(8); // WWJets
	// samples.push_back(9); // ZZJets
	// samples.push_back(10); // WZJets
	makeSSPredictionPlots(samples);

	// makeSSNsigPredictionPlots();
	// makeTTbarNsigPredictionPlots();
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
void MuonPlotter::loadFakeRatio(int sample, int muon){
	if(muon == 0){
		fH2D_fRatio    = fSamples[sample].histos.mu1ratio1;		
		fH1D_fRatioPt  = fSamples[sample].histos.mu1ratio1pt;		
		fH1D_fRatioEta = fSamples[sample].histos.mu1ratio1eta;		
	}
	else if(muon == 1){
		fH2D_fRatio    = fSamples[sample].histos.mu2ratio1;		
		fH1D_fRatioPt  = fSamples[sample].histos.mu2ratio1pt;		
		fH1D_fRatioEta = fSamples[sample].histos.mu2ratio1eta;		
	}
}

//____________________________________________________________________________
void MuonPlotter::producefRatioFromTree(int sample, int muon, const int nptbins, const double *ptbins, const int netabins, const double *etabins){
	gStyle->SetOptStat(0);
	fH2D_fRatio    = new TH2D("fRatio", "Ratio of tight to loose Muons vs Pt vs Eta", nptbins, ptbins, netabins, etabins);
	fH1D_fRatioPt  = new TH1D("fRatioPt",  "Ratio of tight to loose Muons vs Pt", nptbins, ptbins);
	fH1D_fRatioEta = new TH1D("fRatioEta", "Ratio of tight to loose Muons vs Eta", netabins, etabins);
	fH2D_fRatio->Sumw2();
	fH1D_fRatioPt->Sumw2();
	fH1D_fRatioEta->Sumw2();

	TH2D *H_ntight = new TH2D(Form("NTight_%s", fSamples[sample].sname.Data()),  "NTight Muons", nptbins, ptbins, netabins, etabins);
	TH2D *H_nloose = new TH2D(Form("NLoose_%s", fSamples[sample].sname.Data()),  "NLoose Muons", nptbins, ptbins, netabins, etabins);
	H_ntight->Sumw2();
	H_nloose->Sumw2();

	TTree *tree = fSamples[sample].tree;
	double MuPt[2], MuEta[2];
	int MuTight[2];
	int NJets;
	tree->SetBranchAddress("MuPt", &MuPt);
	tree->SetBranchAddress("MuEta", &MuEta);
	tree->SetBranchAddress("MuTight", &MuTight);
	tree->SetBranchAddress("NJets", &NJets);
	Long64_t nentries = tree->GetEntriesFast();
	Long64_t nbytes = 0, nb = 0;
	for( Long64_t jentry = 0; jentry < nentries; jentry++ ){
		tree->GetEntry(jentry);
		if(NJets < fNJetsMin) continue;
		H_nloose->Fill(MuPt[muon], MuEta[muon]); // Tight or loose
		if(MuTight[muon]) H_ntight->Fill(MuPt[muon], MuEta[muon]); // Tight
	}

	fH2D_fRatio->Divide(H_ntight, H_nloose);
	fD_fRatio  = H_ntight->GetEntries()/H_nloose->GetEntries();
	fD_fRatioE = 0.1 * fD_fRatio;

	TH1D *hmuloosept  = H_nloose->ProjectionX();
	TH1D *hmulooseeta = H_nloose->ProjectionY();
	TH1D *hmutightpt  = H_ntight->ProjectionX();
	TH1D *hmutighteta = H_ntight->ProjectionY();

	fH1D_fRatioPt ->Divide(hmutightpt, hmuloosept);
	fH1D_fRatioEta->Divide(hmutighteta, hmulooseeta);
	fH1D_fRatioPt ->SetMinimum(0);
	fH1D_fRatioEta->SetMinimum(0);
	fH1D_fRatioPt ->SetMaximum(0.3);
	fH1D_fRatioEta->SetMaximum(0.3);
	fH1D_fRatioPt ->SetXTitle(convertVarName("MuPt[1]"));
	fH1D_fRatioEta->SetXTitle(convertVarName("MuEta[1]"));
	fH1D_fRatioPt ->SetYTitle("# Tight / # Loose");
	fH1D_fRatioEta->SetYTitle("# Tight / # Loose");
	fH2D_fRatio->SetXTitle(convertVarName("MuPt[1]"));
	fH2D_fRatio->SetYTitle(convertVarName("MuEta[1]"));
	delete H_ntight, H_nloose, hmuloosept, hmulooseeta, hmutightpt, hmutighteta;

	printHisto(fH2D_fRatio,    TString("fRatio") + fSamples[sample].sname,    "Fake Ratio vs pt/eta", "colz");
	printHisto(fH1D_fRatioPt,  TString("fRatioPt") + fSamples[sample].sname,  "Fake Ratio vs pt", "PE1");
	printHisto(fH1D_fRatioEta, TString("fRatioEta") + fSamples[sample].sname, "Fake Ratio vs eta", "PE1");
	if(fVerbose > 0) cout << " " << fSamples[sample].sname <<  ": Average f ratio is " << fD_fRatio << " +/- " << fD_fRatioE << endl;
}

//____________________________________________________________________________
void MuonPlotter::producepRatioFromTree(int sample, int muon, const int nptbins, const double *ptbins, const int netabins, const double *etabins){
	gStyle->SetOptStat(0);
	// This is supposed to be run on a ZJets selection!
	fH2D_pRatio    = new TH2D("pRatio", "Ratio of tight to loose Muons vs Pt vs Eta", nptbins, ptbins, netabins, etabins);
	fH1D_pRatioPt  = new TH1D("pRatioPt",  "Ratio of tight to loose Muons vs Pt", nptbins, ptbins);
	fH1D_pRatioEta = new TH1D("pRatioEta", "Ratio of tight to loose Muons vs Eta", netabins, etabins);
	fH2D_pRatio->Sumw2();
	fH1D_pRatioPt->Sumw2();
	fH1D_pRatioEta->Sumw2();

	TH2D *H_ntight = new TH2D(Form("NTight_%s", fSamples[sample].sname.Data()),  "NTight Muons", nptbins, ptbins, netabins, etabins);
	TH2D *H_nloose = new TH2D(Form("NLoose_%s", fSamples[sample].sname.Data()),  "NLoose Muons", nptbins, ptbins, netabins, etabins);
	H_ntight->Sumw2();
	H_nloose->Sumw2();

	TH1D *H_minv = new TH1D("H_minv", "Invariant Mass of Muons for p Ratio", 100, 50, 120);
	TTree *tree = fSamples[sample].tree;
	double MuPt[2], MuEta[2], MuPhi[2];
	int MuTight[2];
	int NJets;
	tree->SetBranchAddress("MuPt", &MuPt);
	tree->SetBranchAddress("MuEta", &MuEta);
	tree->SetBranchAddress("MuPhi", &MuPhi);
	tree->SetBranchAddress("MuTight", &MuTight);
	tree->SetBranchAddress("NJets", &NJets);
	Long64_t nentries = tree->GetEntriesFast();
	Long64_t nbytes = 0, nb = 0;
	for( Long64_t jentry = 0; jentry < nentries; jentry++ ){
		tree->GetEntry(jentry);
		if(NJets < fNJetsMin) continue;
		TLorentzVector p1, p2;
		p1.SetPtEtaPhiM(MuPt[0], MuEta[0], MuPhi[0], 0.1057);
		p2.SetPtEtaPhiM(MuPt[1], MuEta[1], MuPhi[1], 0.1057);
		double m = (p1+p2).M();
		H_minv->Fill(m);
		if(fabs(91.2 - m) > 10.) continue; // Z mass window cut: 10 GeV

		H_nloose->Fill(MuPt[muon], MuEta[muon]); // Tight or loose
		if(MuTight[muon]) H_ntight->Fill(MuPt[muon], MuEta[muon]); // Tight
	}

	fH2D_pRatio->Divide(H_ntight, H_nloose);
	fD_pRatio  = H_ntight->GetEntries()/H_nloose->GetEntries();
	fD_pRatioE = 0.1 * fD_pRatio;

	TH1D *hmuloosept  = H_nloose->ProjectionX();
	TH1D *hmulooseeta = H_nloose->ProjectionY();
	TH1D *hmutightpt  = H_ntight->ProjectionX();
	TH1D *hmutighteta = H_ntight->ProjectionY();

	fH1D_pRatioPt ->Divide(hmutightpt, hmuloosept);
	fH1D_pRatioEta->Divide(hmutighteta, hmulooseeta);
	fH1D_pRatioPt ->SetMinimum(0.5);
	fH1D_pRatioEta->SetMinimum(0.5);
	fH1D_pRatioPt ->SetXTitle(convertVarName("MuPt[1]"));
	fH1D_pRatioEta->SetXTitle(convertVarName("MuEta[1]"));
	fH1D_pRatioPt ->SetYTitle("# Tight / # Loose");
	fH1D_pRatioEta->SetYTitle("# Tight / # Loose");
	fH2D_pRatio->SetXTitle(convertVarName("MuPt[1]"));
	fH2D_pRatio->SetYTitle(convertVarName("MuEta[1]"));

	H_minv->SetXTitle("m_{#mu#mu} [GeV]");
	printHisto(H_minv, TString("DiMuInvMass_") + fSamples[sample].sname,    "Invariant mass of OS DiMuon Pairs");

	delete H_ntight, H_nloose, hmuloosept, hmulooseeta, hmutightpt, hmutighteta, H_minv;

	printHisto(fH2D_pRatio,    TString("pRatio") + fSamples[sample].sname,    "Fake Ratio vs pt/eta", "colz");
	printHisto(fH1D_pRatioPt,  TString("pRatioPt") + fSamples[sample].sname,  "Fake Ratio vs pt", "PE1");
	printHisto(fH1D_pRatioEta, TString("pRatioEta") + fSamples[sample].sname, "Fake Ratio vs eta", "PE1");
	if(fVerbose > 0) cout << " " << fSamples[sample].sname <<  ": Average p ratio is " << fD_pRatio << " +/- " << fD_pRatioE << endl;
}

//____________________________________________________________________________
void MuonPlotter::producefRatioFromTreeGenMatch(int sample, int muon, const int nptbins, const double *ptbins, const int netabins, const double *etabins){
	gStyle->SetOptStat(0);
	fH2D_fRatio    = new TH2D("fRatio", "Ratio of tight to loose Muons vs Pt vs Eta", nptbins, ptbins, netabins, etabins);
	fH1D_fRatioPt  = new TH1D("fRatioPt",  "Ratio of tight to loose Muons vs Pt", nptbins, ptbins);
	fH1D_fRatioEta = new TH1D("fRatioEta", "Ratio of tight to loose Muons vs Eta", netabins, etabins);
	fH2D_fRatio->Sumw2();
	fH1D_fRatioPt->Sumw2();
	fH1D_fRatioEta->Sumw2();
	fH1D_fRatioPt->SetXTitle(convertVarName("MuPt[1]"));
	fH1D_fRatioEta->SetXTitle(convertVarName("MuEta[1]"));
	fH2D_fRatio->SetXTitle(convertVarName("MuPt[1]"));
	fH2D_fRatio->SetYTitle(convertVarName("MuEta[1]"));
	fH1D_fRatioPt ->SetYTitle("# Tight / # Loose");
	fH1D_fRatioEta->SetYTitle("# Tight / # Loose");

	TH2D *H_ntight = new TH2D(Form("NTight_%s", fSamples[sample].sname.Data()),  "NTight Muons", nptbins, ptbins, netabins, etabins);
	TH2D *H_nloose = new TH2D(Form("NLoose_%s", fSamples[sample].sname.Data()),  "NLoose Muons", nptbins, ptbins, netabins, etabins);
	H_ntight->Sumw2();
	H_nloose->Sumw2();

	TTree *tree = fSamples[sample].tree;
	double MuPt[2], MuEta[2];
	int MuTight[2], MuGenMoID[2];
	tree->SetBranchAddress("MuPt", &MuPt);
	tree->SetBranchAddress("MuEta", &MuEta);
	tree->SetBranchAddress("MuTight", &MuTight);
	tree->SetBranchAddress("MuGenMoID", &MuGenMoID);
	Long64_t nentries = tree->GetEntriesFast();
	Long64_t nbytes = 0, nb = 0;
	for( Long64_t jentry = 0; jentry < nentries; jentry++ ){
		tree->GetEntry(jentry);
		if(abs(MuGenMoID[muon]) == 24 || abs(MuGenMoID[muon]) == 15) continue;
		H_nloose->Fill(MuPt[muon], MuEta[muon]); // Tight or loose
		if(MuTight[muon]) H_ntight->Fill(MuPt[muon], MuEta[muon]); // Tight
	}

	fH2D_fRatio->Divide(H_ntight, H_nloose);
	fD_fRatio  = H_ntight->GetEntries()/H_nloose->GetEntries();
	fD_fRatioE = 0.1 * fD_fRatio;

	TH1D *hmuloosept  = H_nloose->ProjectionX();
	TH1D *hmulooseeta = H_nloose->ProjectionY();
	TH1D *hmutightpt  = H_ntight->ProjectionX();
	TH1D *hmutighteta = H_ntight->ProjectionY();

	fH1D_fRatioPt->Divide(hmutightpt, hmuloosept);
	fH1D_fRatioEta->Divide(hmutighteta, hmulooseeta);
	fH1D_fRatioPt->SetMinimum(0);
	fH1D_fRatioEta->SetMinimum(0);
	delete H_ntight, H_nloose, hmuloosept, hmulooseeta, hmutightpt, hmutighteta;

	printHisto(fH2D_fRatio,    TString("fRatio") + fSamples[sample].sname,    "Fake Ratio vs pt/eta", "colz");
	printHisto(fH1D_fRatioPt,  TString("fRatioPt") + fSamples[sample].sname,  "Fake Ratio vs pt", "PE1");
	printHisto(fH1D_fRatioEta, TString("fRatioEta") + fSamples[sample].sname, "Fake Ratio vs eta", "PE1");
	if(fVerbose > 0) cout << " " << fSamples[sample].sname <<  ": Average f ratio is " << fD_fRatio << " +/- " << fD_fRatioE << endl;
}

//____________________________________________________________________________
void MuonPlotter::producepRatioFromTreeGenMatch(int sample, int muon, const int nptbins, const double *ptbins, const int netabins, const double *etabins){
	gStyle->SetOptStat(0);
	fH2D_pRatio    = new TH2D("pRatio", "Ratio of tight to loose Muons vs Pt vs Eta", nptbins, ptbins, netabins, etabins);
	fH1D_pRatioPt  = new TH1D("pRatioPt",  "Ratio of tight to loose Muons vs Pt", nptbins, ptbins);
	fH1D_pRatioEta = new TH1D("pRatioEta", "Ratio of tight to loose Muons vs Eta", netabins, etabins);
	fH2D_pRatio->Sumw2();
	fH1D_pRatioPt->Sumw2();
	fH1D_pRatioEta->Sumw2();
	fH1D_pRatioPt->SetXTitle(convertVarName("MuPt[1]"));
	fH1D_pRatioEta->SetXTitle(convertVarName("MuEta[1]"));
	fH2D_pRatio->SetXTitle(convertVarName("MuPt[1]"));
	fH2D_pRatio->SetYTitle(convertVarName("MuEta[1]"));
	fH1D_pRatioPt ->SetYTitle("# Tight / # Loose");
	fH1D_pRatioEta->SetYTitle("# Tight / # Loose");

	TH2D *H_ntight = new TH2D(Form("NTight_%s", fSamples[sample].sname.Data()),  "NTight Muons", nptbins, ptbins, netabins, etabins);
	TH2D *H_nloose = new TH2D(Form("NLoose_%s", fSamples[sample].sname.Data()),  "NLoose Muons", nptbins, ptbins, netabins, etabins);
	H_ntight->Sumw2();
	H_nloose->Sumw2();

	TTree *tree = fSamples[sample].tree;
	double MuPt[2], MuEta[2];
	int MuTight[2], MuGenMoID[2];
	tree->SetBranchAddress("MuPt", &MuPt);
	tree->SetBranchAddress("MuEta", &MuEta);
	tree->SetBranchAddress("MuTight", &MuTight);
	tree->SetBranchAddress("MuGenMoID", &MuGenMoID);
	Long64_t nentries = tree->GetEntriesFast();
	Long64_t nbytes = 0, nb = 0;
	for( Long64_t jentry = 0; jentry < nentries; jentry++ ){
		tree->GetEntry(jentry);
		if(abs(MuGenMoID[muon]) != 24 && abs(MuGenMoID[muon]) != 15) continue;
		H_nloose->Fill(MuPt[muon], MuEta[muon]); // Tight or loose
		if(MuTight[muon]) H_ntight->Fill(MuPt[muon], MuEta[muon]); // Tight
	}

	fH2D_pRatio->Divide(H_ntight, H_nloose);
	fD_pRatio  = H_ntight->GetEntries()/H_nloose->GetEntries();
	fD_pRatioE = 0.1 * fD_pRatio;

	TH1D *hmuloosept  = H_nloose->ProjectionX();
	TH1D *hmulooseeta = H_nloose->ProjectionY();
	TH1D *hmutightpt  = H_ntight->ProjectionX();
	TH1D *hmutighteta = H_ntight->ProjectionY();

	fH1D_pRatioPt->Divide(hmutightpt, hmuloosept);
	fH1D_pRatioEta->Divide(hmutighteta, hmulooseeta);
	fH1D_pRatioPt->SetMinimum(0);
	fH1D_pRatioEta->SetMinimum(0);
	delete H_ntight, H_nloose, hmuloosept, hmulooseeta, hmutightpt, hmutighteta;

	printHisto(fH2D_pRatio,    TString("pRatio") + fSamples[sample].sname,    "Fake Ratio vs pt/eta", "colz");
	printHisto(fH1D_pRatioPt,  TString("pRatioPt") + fSamples[sample].sname,  "Fake Ratio vs pt", "PE1");
	printHisto(fH1D_pRatioEta, TString("pRatioEta") + fSamples[sample].sname, "Fake Ratio vs eta", "PE1");
	if(fVerbose > 0) cout << " " << fSamples[sample].sname <<  ": Average p ratio is " << fD_pRatio << " +/- " << fD_pRatioE << endl;
}

//____________________________________________________________________________
double MuonPlotter::getFakeRatio(double pt, double eta){
	if(fH2D_fRatio == NULL || fH2D_fRatio->GetEntries() == 0){
		cout << "MuonPlotter::getFakeRatio ==> fH2D_fRatio not filled!" << endl;
		return -1.0;
	}
	return fH2D_fRatio->GetBinContent(fH2D_fRatio->FindBin(pt, eta));
}

//____________________________________________________________________________
double MuonPlotter::getFakeRatioPt(double pt){
	if(fH1D_fRatioPt == NULL || fH1D_fRatioPt->GetEntries() == 0){
		cout << "MuonPlotter::getFakeRatioPt ==> fH1D_fRatioPt not filled!" << endl;
		return -1.0;
	}
	return fH1D_fRatioPt->GetBinContent(fH1D_fRatioPt->FindBin(pt));
}

//____________________________________________________________________________
double MuonPlotter::getFakeRatioEta(double eta){
	if(fH1D_fRatioEta == NULL || fH1D_fRatioEta->GetEntries() == 0){
		cout << "MuonPlotter::getFakeRatioEta ==> fH1D_fRatioEta not filled!" << endl;
		return -1.0;
	}
	return fH1D_fRatioEta->GetBinContent(fH1D_fRatioEta->FindBin(eta));
}

//____________________________________________________________________________
void MuonPlotter::makeTTbarNsigPredictionPlots(){
	// Fake Ratios: ////////////////////////////////////////////////////////////////////
	producefRatioFromTreeGenMatch(2, 1, gNPtbins, gPtbins, gNEtabins, gEtabins);
	producepRatioFromTreeGenMatch(2, 1, gNPtbins, gPtbins, gNEtabins, gEtabins);

	// producefRatioFromTree(4, 0, gNPtbins, gPtbins, gNEtabins, gEtabins); // QCD-Pt20
	// producepRatioFromTree(5, 0, gNPtbins, gPtbins, gNEtabins, gEtabins); // Zjets

	TH1D *H_nsigpred = new TH1D("Nsigpred", "Predicted N_sig in Pt1 bins", gNPtbins,  gPtbins);

	// Prediction: /////////////////////////////////////////////////////////////////////
	// float wjets_ttbar_scale = fLumiNorm/fSamples[1].lumi;
	// float qcd_ttbar_scale = fLumiNorm/fSamples[3].lumi;   // = 40'223

	H_nsigpred->Add(NsigPredFromFPRatios(2, true)[0]);
	// H_nsigpred->Add(NsigPredFromFPRatios(1, true)[0], wjets_ttbar_scale);
	// H_nsigpred->Add(NsigPredFromFPRatios(3, true)[0], qcd_ttbar_scale);

	// Observation: ////////////////////////////////////////////////////////////////////
	TH1D *H_nsigobs  = new TH1D("Nsigobs",  "Observed N_sig in Pt1 bins",  gNPtbins,  gPtbins);
	NsigObsTTbar(H_nsigobs);

	// cout << "----------------------------------" << endl;
	// cout << " Observed:  " << H_nsigobs->GetBinContent(1)  << " +/- " << H_nsigobs->GetBinError(1) << endl;
	// cout << " Predicted: " << H_nsigpred->GetBinContent(1) << " +/- " << H_nsigpred->GetBinError(1) << endl;
	// cout << "----------------------------------" << endl;

	// Output
	H_nsigobs->SetXTitle(convertVarName("MuPt[1]"));
	H_nsigobs->SetYTitle(Form("Events / %2.0f GeV", fBinWidthScale));
	if(H_nsigobs->GetMaximum() > 1000.) H_nsigobs->GetYaxis()->SetTitleOffset(1.15);
	H_nsigobs = normHistBW(H_nsigobs,   fBinWidthScale);
	H_nsigpred = normHistBW(H_nsigpred, fBinWidthScale);
	plotPredOverlay2HWithRatio(H_nsigobs, "Observed", H_nsigpred, "Predicted");
	plotPredOverlay2HWithRatio(H_nsigobs, "Observed", H_nsigpred, "Predicted", true);
}

//____________________________________________________________________________
void MuonPlotter::makeSSNsigPredictionPlots(){
	// Fake Ratios: ////////////////////////////////////////////////////////////////////
	producefRatioFromTree(0, 0, gNPtbins, gPtbins, gNEtabins, gEtabins); // QCD-Pt20
	producepRatioFromTree(7, 0, gNPtbins, gPtbins, gNEtabins, gEtabins); // Zjets
	// producefRatioFromTreeGenMatch(2, 1, gNPtbins, gPtbins, gNEtabins, gEtabins);
	// producepRatioFromTreeGenMatch(2, 1, gNPtbins, gPtbins, gNEtabins, gEtabins);


	// Prediction: /////////////////////////////////////////////////////////////////////
	TH1D *H_nsigpred = new TH1D("Nsigpred", "Predicted N_sig in Pt1 bins", gNPtbins,  gPtbins);
	bool output = true;
	H_nsigpred->Add(NsigPredFromFPRatios(2, output)[0], fLumiNorm/fSamples[2].lumi); // LM0
	H_nsigpred->Add(NsigPredFromFPRatios(4, output)[0], fLumiNorm/fSamples[4].lumi); // ttbar

	// Observation: ////////////////////////////////////////////////////////////////////
	TH1D *H_nsigobs = new TH1D("Nsigobs", "Observed N_sig in Pt1 bins",  gNPtbins,  gPtbins);
	vector<int> samples;
	samples.push_back(2);
	NsigObsSSLM0(H_nsigobs, samples);

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
	// Fake Ratios: ////////////////////////////////////////////////////////////////////
	producefRatioFromTree(0, 0, gNPtbins, gPtbins, gNEtabins, gEtabins); // QCD-Pt20
	producepRatioFromTree(7, 0, gNPtbins, gPtbins, gNEtabins, gEtabins); // Zjets
	fH1D_fRatioPt->SetMaximum(1.25);
	// plotRatioOverlay2H(fH1D_fRatioPt, "f-Ratio (QCD_Pt20)", fH1D_pRatioPt, "p-Ratio (ZJets)");

	// producefRatioFromTreeGenMatch(4, 1, gNPtbins, gPtbins, gNEtabins, gEtabins);
	// producepRatioFromTreeGenMatch(4, 1, gNPtbins, gPtbins, gNEtabins, gEtabins);
	// fH1D_fRatioPt->SetMaximum(1.1);
	// plotRatioOverlay2H(fH1D_fRatioPt, "f-Ratio (TTbar)", fH1D_pRatioPt, "p-Ratio (TTbar)");

	// Prediction: /////////////////////////////////////////////////////////////////////
	TH1D *H_nsigpred = new TH1D("Nsigpred", "Predicted N_sig in Pt1 bins",       gNPtbins,  gPtbins);
	TH1D *H_nfppred  = new TH1D("Nfppred",  "Predicted N_fp in Pt1 bins",        gNPtbins,  gPtbins);
	TH1D *H_nffpred  = new TH1D("Nffpred",  "Predicted N_ff in Pt1 bins",        gNPtbins,  gPtbins);
	TH1D *H_nFpred   = new TH1D("NFpred",   "Total predicted fakes in Pt1 bins", gNPtbins,  gPtbins);
	bool output = true;
	for(size_t i = 0; i < samples.size(); ++i){
		int index = samples[i];
		float scale = fLumiNorm/fSamples[samples[i]].lumi; // Normalize all to 100/pb
		vector<TH1D*> prediction = NsigPredFromFPRatios(index, output);
		H_nsigpred->Add(prediction[0], scale);
		H_nfppred ->Add(prediction[1], scale);
		H_nffpred ->Add(prediction[2], scale);
	}

	// Observation: ////////////////////////////////////////////////////////////////////
	TH1D *H_nsigobs  = new TH1D("Nsigobs", "Observed N_sig in Pt1 bins",  gNPtbins,  gPtbins);
	NsigObsSSLM0(H_nsigobs, samples);
	TH1D *H_nt2obs  = new TH1D("Nt2obs", "Observed N_t2 in Pt1 bins",  gNPtbins,  gPtbins);
	NT2ObsSS(H_nt2obs, samples);
	TH1D *H_nt2obsttbar = new TH1D("Nt2obsttbar", "Observed N_t2 in Pt1 bins, ttbar only",  gNPtbins,  gPtbins);
	vector<int> ttbarsample; ttbarsample.push_back(4);
	NT2ObsSS(H_nt2obsttbar, ttbarsample);

	// Output
	H_nsigobs->SetXTitle(convertVarName("MuPt[0]"));
	H_nsigobs->SetYTitle(Form("Events / %2.0f GeV", fBinWidthScale));
	H_nsigobs = normHistBW(H_nsigobs, fBinWidthScale);
	H_nsigpred = normHistBW(H_nsigpred, fBinWidthScale);

	H_nFpred->Add(H_nfppred);
	H_nFpred->Add(H_nffpred);
	H_nFpred->SetXTitle(convertVarName("MuPt[0]"));
	H_nFpred->SetYTitle(Form("Events / %2.0f GeV", fBinWidthScale));
	H_nt2obs->SetXTitle(convertVarName("MuPt[0]"));
	H_nt2obs->SetYTitle(Form("Events / %2.0f GeV", fBinWidthScale));
	H_nt2obsttbar->SetXTitle(convertVarName("MuPt[0]"));
	H_nt2obsttbar->SetYTitle(Form("Events / %2.0f GeV", fBinWidthScale));

	H_nFpred = normHistBW(H_nFpred, fBinWidthScale);
	H_nt2obs = normHistBW(H_nt2obs, fBinWidthScale);
	H_nt2obsttbar = normHistBW(H_nt2obsttbar, fBinWidthScale);

	// printHisto(H_nsigobs, H_nsigobs->GetName(), H_nsigobs->GetTitle());
	// printHisto(H_nsigpred, H_nsigpred->GetName(), H_nsigpred->GetTitle());
	// printHisto(H_nt2obs, H_nt2obs->GetName(), H_nt2obs->GetTitle());
	// printHisto(H_nFpred, H_nFpred->GetName(), H_nFpred->GetTitle());

	// plotPredOverlay2HWithRatio(H_nsigobs, "Observed N_{ sig}", H_nsigpred, "Predicted N_{ sig}");
	// plotPredOverlay2HWithRatio(H_nt2obs,  "Observed N_{ t2}",  H_nFpred,   "Predicted Fakes", false, false);
	// plotPredOverlay2HWithRatio(H_nsigobs, "Observed N_{ sig}", H_nsigpred, "Predicted N_{ sig}", true);
	// plotPredOverlay2HWithRatio(H_nt2obs,  "Observed N_{ t2}",  H_nFpred,   "Predicted Fakes", true, false);

	plotPredOverlay3HWithRatio(H_nFpred, "Predicted Fakes", H_nt2obs,  "Obs. N_{ t2} (TTbar + LM0)", H_nt2obsttbar, "Obs. N_{ t2} (TTbar only)", false, false);
	plotPredOverlay3HWithRatio(H_nFpred, "Predicted Fakes", H_nt2obs,  "Obs. N_{ t2} (TTbar + LM0)", H_nt2obsttbar, "Obs. N_{ t2} (TTbar only)", true, false);
}

//____________________________________________________________________________
void MuonPlotter::NsigObsTTbar(TH1D *&hist){
	// Get Nsig observation
	TTree *tree = fSamples[2].tree;
	tree->ResetBranchAddresses();
	double MuPt[2];
	int MuGenMoID[2], MuGenGMoID[2];
	int MuTight[2];
	int NJets;
	tree->SetBranchAddress("MuPt", &MuPt);
	tree->SetBranchAddress("MuGenMoID",  &MuGenMoID);
	tree->SetBranchAddress("MuGenGMoID", &MuGenGMoID);
	tree->SetBranchAddress("MuTight", &MuTight);
	tree->SetBranchAddress("NJets", &NJets);
	Long64_t nentries = tree->GetEntriesFast();
	Long64_t nbytes = 0, nb = 0;
	for( Long64_t jentry = 0; jentry < nentries; jentry++ ){
		Long64_t ientry = tree->LoadTree(jentry);
		if (ientry < 0) break;
		nb = tree->GetEntry(jentry);   nbytes += nb;
		if(NJets < fNJetsMin) continue;
		// // Truth matching to (direct) ttbar, t->W->mu
		// if( abs(MuGenMoID[0]) == 24 && abs(MuGenGMoID[0]) == 6 
		//  && abs(MuGenMoID[1]) == 24 && abs(MuGenGMoID[1]) == 6 ){

		// Truth matching to ttbar, including t->W->tau->mu
		if( ( ( abs(MuGenMoID[0]) == 24 && abs(MuGenGMoID[0]) == 6  )
			||( abs(MuGenMoID[0]) == 15 && abs(MuGenGMoID[0]) == 24 ) )
			&& ( ( abs(MuGenMoID[1]) == 24 && abs(MuGenGMoID[1]) == 6  )
		||( abs(MuGenMoID[1]) == 15 && abs(MuGenGMoID[1]) == 24 ) ) ){

		// // Truth matching to prompt prompt, W->mu
		// if( abs(MuGenMoID[0]) == 24 && abs(MuGenMoID[1]) == 24 ){

			// Apply tight cut!
			if(MuTight[0] != 1 || MuTight[1] != 1) continue;
			hist->Fill(MuPt[0]);
		}
	}
}

//____________________________________________________________________________
void MuonPlotter::NsigObsSSLM0(TH1D *&hist, vector<int> samples){
	hist->Sumw2();
	for(size_t i = 0; i < samples.size(); ++i){
		// Get Nsig observation for SS LM0
		TTree *tree = fSamples[samples[i]].tree;
		tree->ResetBranchAddresses();
		float scale = fLumiNorm / fSamples[samples[i]].lumi;
		double MuPt[2];
		int MuGenMoType[2];
		int MuTight[2];
		int NJets;
		tree->SetBranchAddress("MuPt", &MuPt);
		tree->SetBranchAddress("MuGenMoType",  &MuGenMoType);
		tree->SetBranchAddress("MuTight", &MuTight);
		tree->SetBranchAddress("NJets", &NJets);
		Long64_t nentries = tree->GetEntriesFast();
		Long64_t nbytes = 0, nb = 0;
		for( Long64_t jentry = 0; jentry < nentries; jentry++ ){
			Long64_t ientry = tree->LoadTree(jentry);
			if (ientry < 0) break;
			nb = tree->GetEntry(jentry);   nbytes += nb;
			if(NJets < fNJetsMin) continue;

			// Truth matching, including t->W->tau->mu
			if( (abs(MuGenMoType[0]) == 9 || abs(MuGenMoType[0]) == 4  || abs(MuGenMoType[0]) == 2)
			 && (abs(MuGenMoType[1]) == 9 || abs(MuGenMoType[1]) == 4  || abs(MuGenMoType[1]) == 2) ){
				// Apply tight cut!
				if(MuTight[0] != 1 || MuTight[1] != 1) continue;
				hist->Fill(MuPt[0], scale);
			}
		}
	}
}

//____________________________________________________________________________
void MuonPlotter::NT2ObsSS(TH1D *&hist, vector<int> samples){
	// Get Nsig observation for SS LM0
	hist->Sumw2();
	for(size_t i = 0; i < samples.size(); ++i){
		int index = samples[i];
		TTree *tree = fSamples[index].tree;
		tree->ResetBranchAddresses();
		float scale = fLumiNorm / fSamples[index].lumi;
		double MuPt[2];
		int MuGenMoType[2];
		int MuTight[2];
		int NJets;
		tree->SetBranchAddress("MuPt", &MuPt);
		tree->SetBranchAddress("MuTight", &MuTight);
		tree->SetBranchAddress("NJets", &NJets);
		Long64_t nentries = tree->GetEntriesFast();
		Long64_t nbytes = 0, nb = 0;
		for( Long64_t jentry = 0; jentry < nentries; jentry++ ){
			Long64_t ientry = tree->LoadTree(jentry);
			if (ientry < 0) break;
			nb = tree->GetEntry(jentry);   nbytes += nb;
			if(NJets < fNJetsMin) continue;
			// Apply tight cut!
			// cout << "MuPt[0]=" << MuPt[0] << "  MuTight[0]=" << MuTight[0] << "  MuPt[1]=" << MuPt[1] << "  MuTight[1]=" << MuTight[1] << endl;
			if(MuTight[0] != 1 || MuTight[1] != 1) continue;
			// cout << "  Filling for " << fSamples[index].sname << " ..." << endl;
			hist->Fill(MuPt[0], scale);
		}
	}
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

	// const int NPtbins = 3;
	// const double Ptbins[NPtbins+1] = {10., 15., 20., 35.};
	TH1D *H_nsigpred = new TH1D(Form("Nsigpred_%s", fSamples[sample].sname.Data()), "Predicted N_sig in Pt1 bins", gNPtbins, gPtbins);
	TH1D *H_nfppred = new TH1D(Form("Nfppred_%s", fSamples[sample].sname.Data()), "Predicted N_fp in Pt1 bins", gNPtbins, gPtbins);
	TH1D *H_nffpred = new TH1D(Form("Nffpred_%s", fSamples[sample].sname.Data()), "Predicted N_ff in Pt1 bins", gNPtbins, gPtbins);
	
	TH2D *H_nt2mes =   new TH2D(Form("Nt2mes_%s",   fSamples[sample].sname.Data()), "Measured N_t2 in Pt1/2 bins", gNPtbins, gPtbins, gNPtbins, gPtbins);
	TH2D *H_nt1mes =   new TH2D(Form("Nt1mes_%s",   fSamples[sample].sname.Data()), "Measured N_t1 in Pt1/2 bins", gNPtbins, gPtbins, gNPtbins, gPtbins);
	TH2D *H_nt0mes =   new TH2D(Form("Nt0mes_%s",   fSamples[sample].sname.Data()), "Measured N_t0 in Pt1/2 bins", gNPtbins, gPtbins, gNPtbins, gPtbins);

	// Fill histograms from tree
	TTree *tree = fSamples[sample].tree;
	double MuPt[2];
	int MuTight[2];
	int NJets;
	tree->SetBranchAddress("MuPt", &MuPt);
	tree->SetBranchAddress("MuTight", &MuTight);
	tree->SetBranchAddress("NJets", &NJets);
	Long64_t nentries = tree->GetEntriesFast();
	Long64_t nbytes = 0, nb = 0;
	for( Long64_t jentry = 0; jentry < nentries; jentry++ ){
		tree->GetEntry(jentry);
		if(NJets < fNJetsMin) continue;
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

	for(size_t i = 1; i <= gNPtbins; ++i){
		double pt1 = H_nsigpred->GetBinCenter(i);
		if(fVerbose > 2){
			cout << "=======================================" << endl;
			cout << "Pt1 = " << pt1 << endl;
		}
		double eta1 = 0.0;
		double npppred(0.0), npppredEstat2(0.0), npppredEsyst2(0.0);
		double nfppred(0.0), nfppredEstat2(0.0), nfppredEsyst2(0.0);
		double nffpred(0.0), nffpredEstat2(0.0), nffpredEsyst2(0.0);
		for(size_t j = 1; j <= gNPtbins; ++j){
			if(fVerbose > 2) cout << " --------" << endl;
			double pt2 = H_nsigpred->GetBinCenter(j);
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
		H_nsigpred->SetBinContent(i, npppred);
		H_nfppred ->SetBinContent(i, nfppred);
		H_nffpred ->SetBinContent(i, nffpred);
		H_nsigpred->SetBinError(i, sqrt(npppredEstat2 + npppredEsyst2));
		H_nfppred ->SetBinError(i, sqrt(nfppredEstat2 + nfppredEsyst2));
		H_nffpred ->SetBinError(i, sqrt(nffpredEstat2 + nffpredEsyst2));
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
	vector<TH1D*> res;
	res.push_back(H_nsigpred);
	res.push_back(H_nfppred);
	res.push_back(H_nffpred);
	return res;
}

//____________________________________________________________________________
TH1D* MuonPlotter::NsigPredFromTree(const int sample, const int nbins, const double *ptbins, bool output){
	if(fVerbose > 2) cout << "MuonPlotter::NsigPredFromTree ==> Predicting Nsig from " << fSamples[sample].sname << endl;

	TH1D *H_nsigpred = new TH1D(Form("Nsigpred_%s", fSamples[sample].sname.Data()),  "Predicted N_sig in Pt1 bins",nbins, ptbins);
	TH2D *H_nt2mes = new TH2D(Form("Nt2mes_%s", fSamples[sample].sname.Data()),  "Measured N_t2 in Pt1/2 bins", nbins, ptbins, nbins, ptbins);
	TH2D *H_nt1mes = new TH2D(Form("Nt1mes_%s", fSamples[sample].sname.Data()),  "Measured N_t1 in Pt1/2 bins", nbins, ptbins, nbins, ptbins);
	TH2D *H_nt0mes = new TH2D(Form("Nt0mes_%s", fSamples[sample].sname.Data()),  "Measured N_t0 in Pt1/2 bins", nbins, ptbins, nbins, ptbins);

	// Fill histograms from tree
	TTree *tree = fSamples[sample].tree;
	double MuPt[2];
	int MuTight[2];
	tree->SetBranchAddress("MuPt", &MuPt);
	tree->SetBranchAddress("MuTight", &MuTight);
	Long64_t nentries = tree->GetEntriesFast();
	Long64_t nbytes = 0, nb = 0;
	for( Long64_t jentry = 0; jentry < nentries; jentry++ ){
		tree->GetEntry(jentry);
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
	// Calculate Npp prediction
	// double p = 1.0;
	double p = 0.97;
	for(size_t i = 1; i <= nbins; ++i){
		// double f1 = fH1D_fRatioPt->GetBinContent(i);
		double f1 = getFakeRatioPt(H_nsigpred->GetBinCenter(i));
		double npppred(0.0);		
		for(size_t j = 1; j <= nbins; ++j){
			// double f2 = fH1D_fRatioPt->GetBinContent(j);
			// double f2 = getFakeRatioPt(H_nsigpred->GetBinCenter(j));
			double f2 = 0.108696;
			double nt2 = H_nt2mes->GetBinContent(i,j);
			double nt1 = H_nt1mes->GetBinContent(i,j);
			double nt0 = H_nt0mes->GetBinContent(i,j);
			cout << "   nt2: " << nt2 << "  nt1: " << nt1 << "  nt0: " << nt0 << endl;

			npppred += 1.0/(p-f2) * ( (1.0-f1)*(1.0-f2) * nt2 - f2*(1.0-f1) * nt1 + f1*f2 * nt0 );
			cout << "   f = " << f2 << endl;

		}
		npppred *= 1.0/(p-f1);
		cout << "   p*p*npppred = " << p*p*npppred << endl;
		// cout << i << " : " << npppred << endl;
		H_nsigpred->SetBinContent(i, p*p*npppred); // N_sig = p^2 * N_pp
	}

	if(fVerbose > 2) cout << " Predict " << H_nsigpred->Integral() << " signal events (Nsig) from this sample" << endl;	
	// Output
	H_nsigpred->SetXTitle(convertVarName("MuPt[1]"));
	if(output){
		H_nt2mes->SetXTitle(convertVarName("MuPt[1]"));
		H_nt2mes->SetYTitle(convertVarName("MuPt[1]"));
		H_nt1mes->SetXTitle(convertVarName("MuPt[1]"));
		H_nt1mes->SetYTitle(convertVarName("MuPt[1]"));
		H_nt0mes->SetXTitle(convertVarName("MuPt[1]"));
		H_nt0mes->SetYTitle(convertVarName("MuPt[1]"));
		printHisto(H_nsigpred, H_nsigpred->GetName(), H_nsigpred->GetTitle());
		printHisto(H_nt2mes, H_nt2mes->GetName(), H_nt2mes->GetTitle(), "colz");
		printHisto(H_nt1mes, H_nt1mes->GetName(), H_nt1mes->GetTitle(), "colz");
		printHisto(H_nt0mes, H_nt0mes->GetName(), H_nt0mes->GetTitle(), "colz");
	}
	return H_nsigpred;
}

//____________________________________________________________________________
void MuonPlotter::makeWJetsNt2PredFromTree(){
	// ROOT Memory cleanup...
	if( TH1D *h = (TH1D*)gROOT->FindObject("Nt2obs")) h->Delete();
	if( TH1D *h = (TH1D*)gROOT->FindObject("Nt2pred")) h->Delete();
	if( TH1D *h = (TH1D*)gROOT->FindObject("Nt1mes")) h->Delete();
	if( TH1D *h = (TH1D*)gROOT->FindObject("Nt0mes")) h->Delete();

	// Binning
	const int NPtbins = 3;
	const int NEtabins = 5;
	const double Ptbins[NPtbins+1] = {10., 15., 20., 35.};
	const double Etabins[NEtabins+1] = {-2.4, -1.4, -0.5, 0.5, 1.4, 2.4};

	producefRatioFromTree(0, 0, NPtbins, Ptbins, NEtabins, Etabins);

	TH1D *H_nt2obs  = new TH1D("Nt2obs",  "Observed N_t2 in Pt bins", NPtbins,  Ptbins);
	TH1D *H_nt2pred = new TH1D("Nt2pred",  "Predicted N_t2 in Pt bins", NPtbins,  Ptbins);
	TH2D *H_nt1mes = new TH2D("Nt1mes",  "Measured N_t1 in Pt1/2 bins", NPtbins,  Ptbins, NPtbins, Ptbins);
	TH2D *H_nt0mes = new TH2D("Nt0mes",  "Measured N_t0 in Pt1/2 bins", NPtbins,  Ptbins, NPtbins, Ptbins);

	double MuPt[2];
	int MuTight[2];

	TTree *tree = fSamples[1].tree;
	tree->SetBranchAddress("MuPt", &MuPt);
	tree->SetBranchAddress("MuTight", &MuTight);

	Long64_t nentries = tree->GetEntriesFast();

	Long64_t nbytes = 0, nb = 0;
	for( Long64_t jentry = 0; jentry < nentries; jentry++ ){
		tree->GetEntry(jentry);
		// MuTight is 0 for loosenotight and 1 for tight
		if( MuTight[0] == 1 && MuTight[1] == 1 ) H_nt2obs->Fill(MuPt[0]);
		if( MuTight[0] == 1 && MuTight[1] == 0 ) H_nt1mes->Fill(MuPt[0], MuPt[1]);
		if( MuTight[0] == 0 && MuTight[1] == 1 ) H_nt1mes->Fill(MuPt[1], MuPt[0]);
		if( MuTight[1] == 0 && MuTight[0] == 0 ) H_nt0mes->Fill(MuPt[0], MuPt[1]);
	}

	for(size_t i = 1; i <= NPtbins; ++i){
		double f1 = fH1D_fRatioPt->GetBinContent(i);
		double nt2pred(0.0);		
		for(size_t j = 1; j <= NPtbins; ++j){
			double f2 = fH1D_fRatioPt->GetBinContent(j);
			double nt1 = H_nt1mes->GetBinContent(i,j);
			double nt0 = H_nt0mes->GetBinContent(i,j);
			nt2pred += f2 * nt1 - f1*f2/(1.0-f2) * nt0;
		}
		nt2pred *= 1.0/(1.0 - f1);
		// cout << i << " : " << nt2pred << endl;
		H_nt2pred->SetBinContent(i, nt2pred);
	}
	H_nt2obs->SetXTitle(convertVarName("MuPt[1]"));
	H_nt2pred->SetXTitle(convertVarName("MuPt[1]"));
	// H_nt2obs = normHistBW(H_nt2obs);
	// H_nt2pred = normHistBW(H_nt2pred);
	plotPredOverlay2H(H_nt2obs, "Observed", H_nt2pred, "Predicted");


	printHisto(H_nt1mes, H_nt1mes->GetName(), H_nt1mes->GetTitle(), "colz");
	printHisto(H_nt0mes, H_nt0mes->GetName(), H_nt0mes->GetTitle(), "colz");
	printHisto(fH2D_fRatio,    "fRatio",    "Fake Ratio vs pt/eta", "colz");
	printHisto(fH1D_fRatioPt,  "fRatioPt",  "Fake Ratio vs pt", "PE1");
	printHisto(fH1D_fRatioEta, "fRatioEta", "Fake Ratio vs eta", "PE1");
}

//____________________________________________________________________________
void MuonPlotter::makeWJetsNt2PredFromHistos(){
	// ROOT Memory cleanup...
	if( TH1D *h = (TH1D*)gROOT->FindObject("Nt2obs")) h->Delete();
	if( TH1D *h = (TH1D*)gROOT->FindObject("Nt2pred")) h->Delete();
	if( TH1D *h = (TH1D*)gROOT->FindObject("Nt1mes")) h->Delete();
	if( TH1D *h = (TH1D*)gROOT->FindObject("Nt0mes")) h->Delete();
	const int nbinsx = fH1D_fRatioPt->GetNbinsX();
	// double *xbins = getBinning(fH1D_fRatioPt);
	double xbins[nbinsx+1];
	for(size_t i = 0; i < nbinsx+1; ++i){
		xbins[i] = fH1D_fRatioPt->GetBinLowEdge(i+1);
	}

	TH1D *H_nt2obs  = new TH1D("Nt2obs",  "Observed N_t2 in Pt bins", nbinsx,  xbins);
	TH1D *H_nt2pred = new TH1D("Nt2pred",  "Predicted N_t2 in Pt bins", nbinsx,  xbins);
	TH2D *H_nt1mes = new TH2D("Nt1mes",  "Measured N_t1 in Pt1/2 bins", nbinsx,  xbins, nbinsx, xbins);
	TH2D *H_nt0mes = new TH2D("Nt0mes",  "Measured N_t0 in Pt1/2 bins", nbinsx,  xbins, nbinsx, xbins);

	H_nt2obs = fSamples[1].histos.munt2->ProjectionX("Nt2obs");
	H_nt1mes = fSamples[1].histos.munt1;
	H_nt0mes = fSamples[1].histos.munt0;

	for(size_t i = 1; i <= nbinsx; ++i){
		double f1 = fH1D_fRatioPt->GetBinContent(i);
		double nt2pred(0.0);		
		for(size_t j = 1; j <= nbinsx; ++j){
			double f2 = fH1D_fRatioPt->GetBinContent(j);
			double nt1 = H_nt1mes->GetBinContent(i,j);
			double nt0 = H_nt0mes->GetBinContent(i,j);
			nt2pred += f2 * nt1 - f1*f2/(1.0-f2) * nt0;
		}
		nt2pred *= 1.0/(1.0 - f1);
		// cout << i << " : " << nt2pred << endl;
		H_nt2pred->SetBinContent(i, nt2pred);
	}
	H_nt2obs->SetName("Nt2obs_FromHistos");
	H_nt2pred->SetName("Nt2pred_FromHistos");
	H_nt2obs->SetXTitle(convertVarName("MuPt[1]"));
	H_nt2pred->SetXTitle(convertVarName("MuPt[1]"));
	plotPredOverlay2H(H_nt2obs, "Observed", H_nt2pred, "Predicted");
}

//____________________________________________________________________________
void MuonPlotter::makeSSFRPredictionPlots(){
	// fH1D_SSMuPtObs  = new TH1D("SSMuPtObs", "SS Muons Pt Observed",  gNPtbins, gPtbins);
	// fH1D_SSMuPtPred = new TH1D("SSMuPtPred", "SS Muons Pt Predicted", gNPtbins, gPtbins);
	// fH1D_SSMuEtaObs  = new TH1D("SSMuEtaObs", "SS Muons Eta Observed",  gNEtabins, gEtabins);
	// fH1D_SSMuEtaPred = new TH1D("SSMuEtaPred", "SS Muons Eta Predicted", gNEtabins, gEtabins);
	// 
	// // These two should contain purely FAKE second muons
	// TH1D *hloosenotightpt = new TH1D("loosenotightpt", "loosenotight pt", gNPtbins, gPtbins);
	// TH1D *hloosenotighteta = new TH1D("loosenotighteta", "loosenotight eta", gNEtabins, gEtabins);
	// 
	// TTree *t = fSamples[6].tree;
	// TCut c   = fSamples[6].cut; // c here ensures gen match of first muon to W
	// 
	// TCut muid1 = "abs(MuD0[0])<0.005&&MuCaloComp[0]>0.5&&MuSegmComp[0]>0.5";
	// TCut muid2 = "abs(MuD0[1])<0.005&&MuCaloComp[1]>0.5&&MuSegmComp[1]>0.5";
	// TCut common = muid1 && muid2 && c;
	// TCut loose = "MuIso[1]<1.00" && common;
	// TCut tight = "MuIso[1]<0.15" && common;
	// TCut loosenotight = "MuIso[1]>0.15&&MuIso[1]<1.0" && common;
	// 
	// t->Project("loosenotightpt" , "MuPt[1]",  loosenotight);
	// t->Project("loosenotighteta", "MuEta[1]", loosenotight);
	// 
	// t->Project("SSMuPtObs" , "MuPt[1]",  tight);
	// t->Project("SSMuEtaObs", "MuEta[1]", tight);
	// 
	// // printHisto(fH1D_SSMuPtObs, "Test", "Test of Observation");
	// 
	// fH1D_SSMuPtPred->Multiply(fH1D_MuRatio2Pt, hloosenotightpt);
	// fH1D_SSMuEtaPred->Multiply(fH1D_MuRatio2Eta, hloosenotighteta);
	// 
	// fH1D_SSMuPtPred->SetDrawOption("PE1");
	// fH1D_SSMuPtPred->SetFillColor(kRed);
	// fH1D_SSMuPtPred->SetLineColor(kRed);
	// fH1D_SSMuEtaPred->SetDrawOption("PE1");
	// fH1D_SSMuEtaPred->SetFillColor(kRed);
	// fH1D_SSMuEtaPred->SetLineColor(kRed);
	// fH1D_SSMuPtObs->SetDrawOption("hist");
	// fH1D_SSMuPtObs->SetFillColor(kBlue);
	// fH1D_SSMuPtObs->SetLineColor(kBlue);
	// fH1D_SSMuEtaObs->SetDrawOption("hist");
	// fH1D_SSMuEtaObs->SetFillColor(kBlue);
	// fH1D_SSMuEtaObs->SetLineColor(kBlue);
	// 
	// plotOverlay2H(fH1D_SSMuPtPred, "Predicted", fH1D_SSMuPtObs, "Observed", false);
	// plotOverlay2H(fH1D_SSMuEtaPred, "Predicted", fH1D_SSMuEtaObs, "Observed", false);
}

//____________________________________________________________________________
void MuonPlotter::makeWJetsFRPredictionPlotsOld(){
	ratio_histos histos = fSamples[1].histos;
	ratio_histos qcd_histos = fSamples[0].histos;

	TH1D *WJetsMuPtPred  = new TH1D(*(histos.mu2tight->ProjectionX("WJetsMuPtPred")));
	TH1D *WJetsMuEtaPred = new TH1D(*(histos.mu2tight->ProjectionY("WJetsMuEtaPred")));
	TH1D *WJetsMuPtObs   = histos.mu2tight->ProjectionX("WJetsMuPtObs");
	TH1D *WJetsMuEtaObs  = histos.mu2tight->ProjectionY("WJetsMuEtaObs");

	// These two should contain purely FAKE second muons
	TH1D *hloosenotightpt  = histos.mu2loosenotight->ProjectionX();
	TH1D *hloosenotighteta = histos.mu2loosenotight->ProjectionY();

	WJetsMuPtPred->Multiply( qcd_histos.mu1ratio2pt,  hloosenotightpt);
	WJetsMuEtaPred->Multiply(qcd_histos.mu1ratio2eta, hloosenotighteta);

	WJetsMuPtPred->SetXTitle(convertVarName("MuPt[1]"));
	WJetsMuEtaPred->SetXTitle(convertVarName("MuEta[1]"));
	WJetsMuPtObs->SetXTitle(convertVarName("MuPt[1]"));
	WJetsMuEtaObs->SetXTitle(convertVarName("MuEta[1]"));

	plotOverlay2H(WJetsMuPtObs,  "Observed", WJetsMuPtPred,  "Predicted");
	plotOverlay2H(WJetsMuEtaObs, "Observed", WJetsMuEtaPred, "Predicted");

	cout << "----------------------" << endl;
	for(size_t i = 1; i <= WJetsMuPtObs->GetNbinsX(); ++i){
		cout << "Observed:  " << WJetsMuPtObs->GetBinContent(i) << " +/- " << WJetsMuPtObs->GetBinError(i) << endl;
		cout << "Predicted: " << WJetsMuPtPred->GetBinContent(i) << " +/- " << WJetsMuPtPred->GetBinError(i) << endl;
	}
	cout << "----------------------" << endl;
	for(size_t i = 1; i <= WJetsMuEtaObs->GetNbinsX(); ++i){
		cout << "Observed:  " << WJetsMuEtaObs->GetBinContent(i) << " +/- " << WJetsMuEtaObs->GetBinError(i) << endl;
		cout << "Predicted: " << WJetsMuEtaPred->GetBinContent(i) << " +/- " << WJetsMuEtaPred->GetBinError(i) << endl;
	}
	cout << "----------------------" << endl;
}

//____________________________________________________________________________
void MuonPlotter::plotOrigin(int s1, int s2, int s3){
	// if(fSamples.size() < 3){
	// 	cout << " MuonPlotter::makeOverlay3 ==> Samples not defined!" << endl;
	// 	return;
	// }
	// TString sn1 = fSamples[s1].sname;
	// TCut c1     = fSamples[s1].cut;
	// TTree *t1   = fSamples[s1].tree;
	// 
	// TString sn2 = fSamples[s2].sname;
	// TCut c2     = fSamples[s2].cut;
	// TTree *t2   = fSamples[s2].tree;
	// 
	// TString sn3 = fSamples[s3].sname;
	// TCut c3     = fSamples[s3].cut;
	// TTree *t3   = fSamples[s3].tree;
	// 
	// char name[100];
	// char arg[100];
	// double xbins[23] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23};
	// sprintf(arg, "%s", "MuGenMoType[0]");
	// sprintf(name, "%s_%s", arg, sn1.Data());
	// TH1D *h1 = drawTree1D(arg, c1, name, 22, xbins, t1);
	// sprintf(name, "%s_%s", arg, sn2.Data());
	// TH1D *h2 = drawTree1D(arg, c2, name, 22, xbins, t2);
	// sprintf(name, "%s_%s", arg, sn3.Data());
	// TH1D *h3 = drawTree1D(arg, c3, name, 22, xbins, t3);
	// 
	// h1->GetXaxis()->SetBinLabel(1 ,"Quark");
	// h1->GetXaxis()->SetBinLabel(2 ,"Lepton");
	// h1->GetXaxis()->SetBinLabel(3 ,"Excited Particle");
	// h1->GetXaxis()->SetBinLabel(4 ,"Gauge Boson");
	// h1->GetXaxis()->SetBinLabel(5 ,"Diquark");
	// h1->GetXaxis()->SetBinLabel(6 ,"Technicolor Particle");
	// h1->GetXaxis()->SetBinLabel(7 ,"R-Hadron");
	// h1->GetXaxis()->SetBinLabel(8 ,"Special Particle");
	// h1->GetXaxis()->SetBinLabel(9 ,"Susy Particle");
	// h1->GetXaxis()->SetBinLabel(10,"Kaluza-Klein Excitation");
	// h1->GetXaxis()->SetBinLabel(11,"Light I=1 Meson");
	// h1->GetXaxis()->SetBinLabel(12,"Light I=0 Meson");
	// h1->GetXaxis()->SetBinLabel(13,"s Meson");
	// h1->GetXaxis()->SetBinLabel(14,"c Meson");
	// h1->GetXaxis()->SetBinLabel(15,"b Meson");
	// h1->GetXaxis()->SetBinLabel(16,"c#bar{c} Meson");
	// h1->GetXaxis()->SetBinLabel(17,"b#bar{b} Meson");
	// h1->GetXaxis()->SetBinLabel(18,"Light Baryon");
	// h1->GetXaxis()->SetBinLabel(19,"s Baryon");
	// h1->GetXaxis()->SetBinLabel(20,"c Baryon");
	// h1->GetXaxis()->SetBinLabel(21,"b Baryon");
	// h1->GetXaxis()->SetBinLabel(22,"Pentaquark");
	// 
	// // Remove labels from empty bins
	// for(size_t i = 1; i < 23; ++i){
	// 	if( h1->GetBinContent(i)==0 && h2->GetBinContent(i)==0 && h3->GetBinContent(i)==0) h1->GetXaxis()->SetBinLabel(i,"");
	// }
	// 
	// gStyle->SetOptStat("");
	// 
	// h1->SetLineWidth(2);
	// h1->SetFillColor(fSamples[s1].color);
	// h1->SetFillStyle(1001);
	// h2->SetLineWidth(2);
	// h2->SetLineColor(fSamples[s2].color);
	// h2->SetFillColor(fSamples[s2].color);
	// h2->SetFillStyle(3005);
	// h3->SetLineWidth(2);
	// h3->SetLineColor(fSamples[s3].color);
	// h3->SetFillColor(fSamples[s3].color);
	// h3->SetFillStyle(3004);
	// 
	// char canvtitle[100], canvname[100];
	// sprintf(canvtitle,"%s vs %s vs %s", h1->GetName(), h2->GetName(), h3->GetName());
	// sprintf(canvname,"%s:%s:%s", h1->GetName(), h2->GetName(), h3->GetName());
	// TCanvas *col = new TCanvas(canvname, canvtitle, 0, 0, 900, 700);
	// col->SetFillStyle(0);
	// col->SetFrameFillStyle(0);
	// col->cd();
	// gPad->SetFillStyle(0);
	// h1->Sumw2();
	// h2->Sumw2();
	// h3->Sumw2();
	// 
	// h1->Scale(1.0/h1->Integral());
	// h2->Scale(1.0/h2->Integral());
	// h3->Scale(1.0/h3->Integral());
	// 
	// // Determine plotting range
	// double max1 = h1->GetMaximum();
	// double max2 = h2->GetMaximum();
	// double max3 = h3->GetMaximum();
	// double max12 = (max1>max2)?max1:max2;
	// double max = (max12>max3)?max12:max3;
	// max = 1.05*max;
	// h1->SetMaximum(max);
	// h2->SetMaximum(max);
	// h3->SetMaximum(max);
	// 
	// TLegend *leg = new TLegend(0.65,0.69,0.886,0.88);
	// leg->AddEntry(h1, sn1,"f");
	// leg->AddEntry(h2, sn2,"f");
	// leg->AddEntry(h3, sn3,"f");
	// leg->SetFillStyle(0);
	// leg->SetTextFont(42);
	// leg->SetBorderSize(0);
	// 
	// h1->DrawCopy("hist");
	// h2->DrawCopy("histsame");
	// h3->DrawCopy("histsame");
	// leg->Draw();
	// 
	// gPad->RedrawAxis();
	// TString outputname = TString(h1->GetName()) + "_" + TString(h2->GetName()) + "_" + TString(h3->GetName());
	// Util::PrintNoEPS(col, outputname, fOutputDir, fOutputFile);
}

//____________________________________________________________________________
void MuonPlotter::fillMuControlPlots(TCut loose, TCut tight, int s){
	// gStyle->SetOptStat(111111);
	// TH2D *MuLoose =        new TH2D("hMuLoose", "Number of Loose Muons vs Pt and Eta", gNPtbins, gPtbins, gNEtabins, gEtabins);
	// TH2D *MuTight =        new TH2D("hMuTight", "Number of Tight Muons vs Pt and Eta", gNPtbins, gPtbins, gNEtabins, gEtabins);
	// TH2D *MuLooseNoTight = new TH2D("hMuLooseNoTight", "Number of Loose but not Tight Muons vs Pt and Eta", gNPtbins, gPtbins, gNEtabins, gEtabins);
	// 
	// MuTight->SetXTitle(convertVarName("MuPt[1]"));
	// MuTight->SetYTitle(convertVarName("MuEta[1]"));
	// MuLoose->SetXTitle(convertVarName("MuPt[1]"));
	// MuLoose->SetYTitle(convertVarName("MuEta[1]"));
	// MuLooseNoTight->SetXTitle(convertVarName("MuPt[1]"));
	// MuLooseNoTight->SetYTitle(convertVarName("MuEta[1]"));
	// 
	// TTree *t   = fSamples[s].tree;
	// t->Project("hMuLoose", "MuEta[1]:MuPt[1]", loose);
	// t->Project("hMuTight", "MuEta[1]:MuPt[1]", tight);
	// 
	// MuLoose->Sumw2();
	// MuTight->Sumw2();
	// MuLooseNoTight->Sumw2();
	// 
	// MuLooseNoTight->Add(MuLoose);
	// MuLooseNoTight->Add(MuTight, -1);
	// 
	// MuLoose->SetMinimum(0);
	// MuTight->SetMinimum(0);
	// MuLooseNoTight->SetMinimum(0);
	// 
	// Option_t *drawopt = "text colz";
	// 
	// TString tag = fSamples[s].sname;
	// printHisto(MuTight, tag + "-MuTight", "Tight Muons", drawopt);
	// printHisto(MuLoose, tag + "-MuLoose", "Loose Muons", drawopt);
	// printHisto(MuLooseNoTight, tag + "-MuLooseNoTight", "Loose No Tight Muons", drawopt);
}
