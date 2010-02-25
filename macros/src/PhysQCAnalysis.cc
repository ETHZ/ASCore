#include "PhysQCAnalysis.hh"
#include "base/TreeReader.hh"
#include "helper/AnaClass.hh"
#include "base/TreeAnalyzerBase.hh"
#include "TreeCleaner.hh"
#include "helper/Utilities.hh"

using namespace std;

PhysQCAnalysis::PhysQCAnalysis(TreeReader *tr, TreeCleaner *tc) : UserAnalysisBase(tr){
	fTR = tr;
	fTC = tc;
	fTlat = new TLatex();
	fAC = new AnaClass();
	const unsigned int nmubins = 100;
	Util::SetStyle();
	fMuHistos[0] = new TH1D("MuDeltaPOverP", "Mu DeltaP/P",    nmubins, 0., 0.25 );
	fMuHistos[1] = new TH1D("MuNChi2", "Mu Chi2/Ndof",         nmubins, 0., 20. );
	fMuHistos[2] = new TH1D("MuNTkHits", "Mu # Tracker Hits",  nmubins, 0., 25. );
	
	const unsigned int nevbins = 100;
	fMETHistos[0] = new TH1D("METR12", "MET R12", nevbins, 0, 6.283);
	fMETHistos[1] = new TH1D("METR21", "MET R21", nevbins, 0, 6.283);
}

PhysQCAnalysis::~PhysQCAnalysis(){
}

void PhysQCAnalysis::Begin(){
	fAC->readVarNames("varnames.dat");
	fAC->setOutputDir(fOutputDir);
	fAC->setGlobalTag(fTag);
}

void PhysQCAnalysis::Analyze(){
	for(size_t i = 0; i < fTR->NMus; ++i){
		fMuHistos[0]->Fill(fTR->MuPtE[i]/fTR->MuPt[i]);
		fMuHistos[1]->Fill(fTR->MuNChi2[i]);
		fMuHistos[2]->Fill(fTR->MuNTkHits[i]);
	}
	fMETHistos[0]->Fill(fTC->fR12);
	fMETHistos[1]->Fill(fTC->fR21);
}

void PhysQCAnalysis::End(){
	TCanvas *canv;
	TString tempstring1, tempstring2, outputdir;
	outputdir = fOutputDir + "Cleaning/Mu/";
	tempstring1 = "Mu DeltaP/P";
	tempstring2 = "MuDeltaPOverP";
	canv = new TCanvas(tempstring2, tempstring1 , 0, 0, 900, 700);
	fMuHistos[0]->DrawCopy();
	Util::PrintPNG(canv, tempstring2, outputdir);
	Util::PrintEPS(canv, tempstring2, outputdir);
	
	tempstring1 = "Mu NChi2";
	tempstring2 = "MuNChi2";
	canv = new TCanvas(tempstring2, tempstring1 , 0, 0, 900, 700);
	fMuHistos[1]->DrawCopy();
	Util::PrintPNG(canv, tempstring2, outputdir);
	Util::PrintEPS(canv, tempstring2, outputdir);
	
	tempstring1 = "Mu NTkHits";
	tempstring2 = "MuNTkHits";
	canv = new TCanvas(tempstring2, tempstring1 , 0, 0, 900, 700);
	fMuHistos[2]->DrawCopy();
	Util::PrintPNG(canv, tempstring2, outputdir);
	Util::PrintEPS(canv, tempstring2, outputdir);

	outputdir = fOutputDir + "Cleaning/MET/";
	tempstring1 = "MET R12";
	tempstring2 = "METR12";
	canv = new TCanvas(tempstring2, tempstring1 , 0, 0, 900, 700);
	fMETHistos[0]->DrawCopy();
	Util::PrintPNG(canv, tempstring2, outputdir);
	Util::PrintEPS(canv, tempstring2, outputdir);

	tempstring1 = "MET R21";
	tempstring2 = "METR21";
	canv = new TCanvas(tempstring2, tempstring1 , 0, 0, 900, 700);
	fMETHistos[1]->DrawCopy();
	Util::PrintPNG(canv, tempstring2, outputdir);
	Util::PrintEPS(canv, tempstring2, outputdir);
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
	Util::PrintPNG(canv, "HLTStats1", fOutputDir);
	Util::PrintEPS(canv, "HLTStats1", fOutputDir);

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
	Util::PrintPNG(canv, "HLTStats2", fOutputDir);
	Util::PrintEPS(canv, "HLTStats2", fOutputDir);
	
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
	Util::PrintPNG(canv, "L1PStats1", fOutputDir);
	Util::PrintEPS(canv, "L1PStats1", fOutputDir);
	
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
	Util::PrintPNG(canv, "L1PStats2", fOutputDir);
	Util::PrintEPS(canv, "L1PStats2", fOutputDir);
	
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
	Util::PrintPNG(canv, "L1TStats", fOutputDir);
	Util::PrintEPS(canv, "L1TStats", fOutputDir);
}
