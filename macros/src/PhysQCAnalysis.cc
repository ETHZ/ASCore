#include "PhysQCAnalysis.hh"
#include "TreeReader.hh"
#include "AnaClass.hh"

using namespace std;

PhysQCAnalysis::PhysQCAnalysis(TreeReader *tr) : UserAnalysisBase(tr){
	fTR = tr;
	fTlat = new TLatex();
	fAC = new AnaClass();
}

PhysQCAnalysis::~PhysQCAnalysis(){
}

void PhysQCAnalysis::Begin(const char* filename){
	fAC->readVarNames("varnames.dat");
	fAC->setOutputDir(fOutputDir);
	fAC->setGlobalTag(fTag);
}

void PhysQCAnalysis::Analyze(){
	
}

void PhysQCAnalysis::End(){
}

void PhysQCAnalysis::MakePlots(TString plotlist, TTree *tree){
	fAC->plotPlotList(plotlist, tree, "");
	fAC->plotEID("", tree, "");
}

void PhysQCAnalysis::PlotTriggerStats(){
	const int nentries = fTR->fChain->GetEntries();
	TFile *f = fTR->fChain->GetCurrentFile();
	f->cd("analyze");
	TH1I *hlt_stats = (TH1I*)gDirectory->Get("HLTTriggerStats");
	hlt_stats->GetXaxis()->LabelsOption("v");
	hlt_stats->GetXaxis()->SetLabelSize(0.035);
	hlt_stats->SetMaximum(nentries + 0.05*nentries);
	TH1I *l1p_stats = (TH1I*)gDirectory->Get("L1PhysTriggerStats");
	l1p_stats->GetXaxis()->LabelsOption("v");
	l1p_stats->GetXaxis()->SetLabelSize(0.033);
	l1p_stats->SetMaximum(nentries + 0.05*nentries);
	TH1I *l1t_stats = (TH1I*)gDirectory->Get("L1TechTriggerStats");
	l1t_stats->SetMaximum(nentries + 0.05*nentries);

	fTlat->SetTextColor(kBlack);
	fTlat->SetNDC(kTRUE);
	fTlat->SetTextSize(0.04);
	gStyle->SetOptStat(0);

	TString entries = Form("# Entries = %d", nentries);
	TLine *l1;

	TCanvas *canv;
	TString tempstring;
	tempstring = "HLT Trigger Bits (0-50)";
	canv = new TCanvas("HLTStats1", tempstring , 0, 0, 900, 700);
	canv->SetBottomMargin(0.50);
	gPad->SetGridy();
	hlt_stats->GetXaxis()->SetRange(0,50);
	hlt_stats->DrawCopy();
	l1 = new TLine(0, nentries, 50, nentries);
	l1->SetLineColor(kRed);
	l1->SetLineWidth(2);
	l1->SetLineStyle(2);
	l1->Draw();
	fTlat->DrawLatex(0.11,0.92, tempstring);
	fTlat->DrawLatex(0.60,0.92, entries);
	printPNG(canv, "HLTStats1", fOutputDir);
	printEPS(canv, "HLTStats1", fOutputDir);

	tempstring = "HLT Trigger Bits (51-100)";
	canv = new TCanvas("HLTStats2", tempstring , 0, 0, 900, 700);
	canv->SetBottomMargin(0.50);
	gPad->SetGridy();
	hlt_stats->GetXaxis()->SetRange(51,100);
	hlt_stats->DrawCopy();
	l1 = new TLine(50, nentries, 100, nentries);
	l1->SetLineColor(kRed);
	l1->SetLineWidth(2);
	l1->SetLineStyle(2);
	l1->Draw();
	fTlat->DrawLatex(0.11,0.92, tempstring);
	fTlat->DrawLatex(0.60,0.92, entries);	
	printPNG(canv, "HLTStats2", fOutputDir);
	printEPS(canv, "HLTStats2", fOutputDir);
	
	tempstring = "L1 Phys Bits (0-63)";
	canv = new TCanvas("L1PStats1", tempstring , 0, 0, 900, 700);
	canv->SetBottomMargin(0.42);
	gPad->SetGridy();
	l1p_stats->GetXaxis()->SetRange(0,63);
	l1p_stats->DrawCopy();
	l1 = new TLine(0, nentries, 63, nentries);
	l1->SetLineColor(kRed);
	l1->SetLineWidth(2);
	l1->SetLineStyle(2);
	l1->Draw();
	fTlat->DrawLatex(0.11,0.92, tempstring);
	fTlat->DrawLatex(0.60,0.92, entries);
	printPNG(canv, "L1PStats1", fOutputDir);
	printEPS(canv, "L1PStats1", fOutputDir);
	
	tempstring = "L1 Phys Bits (64-128)";
	canv = new TCanvas("L1PStats2", tempstring , 0, 0, 900, 700);
	canv->SetBottomMargin(0.42);
	gPad->SetGridy();
	l1p_stats->GetXaxis()->SetRange(64,128);
	l1p_stats->DrawCopy();
	l1 = new TLine(63, nentries, 128, nentries);
	l1->SetLineColor(kRed);
	l1->SetLineWidth(2);
	l1->SetLineStyle(2);
	l1->Draw();
	fTlat->DrawLatex(0.11,0.92, tempstring);
	fTlat->DrawLatex(0.60,0.92, entries);
	printPNG(canv, "L1PStats2", fOutputDir);
	printEPS(canv, "L1PStats2", fOutputDir);
	
	tempstring = "L1 Tech Bits";
	canv = new TCanvas("L1TStats", tempstring , 0, 0, 900, 700);
	gPad->SetGridy();
	l1t_stats->GetXaxis()->SetRange(0,64);
	l1t_stats->DrawCopy();
	l1 = new TLine(0, nentries, 64, nentries);
	l1->SetLineColor(kRed);
	l1->SetLineWidth(2);
	l1->SetLineStyle(2);
	l1->Draw();
	fTlat->DrawLatex(0.11,0.92, tempstring);
	fTlat->DrawLatex(0.60,0.92, entries);
	printPNG(canv, "L1TStats", fOutputDir);
	printEPS(canv, "L1TStats", fOutputDir);
}
