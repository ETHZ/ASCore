#include "PhysQCAnalysis.hh"
#include "base/TreeReader.hh"
#include "helper/AnaClass.hh"
#include "base/TreeAnalyzerBase.hh"
#include "TreeCleaner.hh"
#include "helper/Utilities.hh"

using namespace std;

PhysQCAnalysis::PhysQCAnalysis(TreeReader *tr, TreeCleaner *tc) : UserAnalysisBase(tr){
	fTC = tc;
	fTlat = new TLatex();
	fAC = new AnaClass();
	Util::SetStyle();

	const unsigned int nmubins = 100;
	fMuHistos[0] = new TH1D("MuDeltaPOverP", "Mu DeltaP/P",    nmubins, 0., 0.25 );
	fMuHistos[1] = new TH1D("MuNChi2", "Mu Chi2/Ndof",         nmubins, 0., 20. );
	fMuHistos[2] = new TH1D("MuNTkHits", "Mu # Tracker Hits",  nmubins, 0., 25. );

	const unsigned int nelbins = 100;
	fElHistos[0]  = new TH1D("ElHcalOverEcalBar", "El H/E Bar",  nmubins, 0., 0.05);
	fElHistos[1]  = new TH1D("ElSigmaIetaIetaBar", "El sigmaIetaIeta Bar",  nmubins, 0., 0.05);
	fElHistos[2]  = new TH1D("ElDeltaPhiSeedClusterAtCaloBar", "El DeltaPhiSeedClusterAtCalo Bar",  nmubins, -0.1, 0.1);
	fElHistos[3]  = new TH1D("ElDeltaEtaSeedClusterAtCaloBar", "El DeltaEtaSeedClusterAtCalo Bar",  nmubins, -0.05, 0.05);
	fElHistos[4]  = new TH1D("ElDeltaPhiSuperClusterAtVtxBar", "El DeltaPhiSuperClusterAtVtx Bar",  nmubins, -0.1, 0.1);
	fElHistos[5]  = new TH1D("ElDeltaEtaSuperClusterAtVtxBar", "El DeltaEtaSuperClusterAtVtx Bar",  nmubins, -0.05, 0.05);
	fElHistos[6]  = new TH1D("ElESuperClusterOverPBar", "El ESuperClusterOverP Bar",  nmubins, 0., 5.);
	fElHistos[7]  = new TH1D("ElHcalOverEcalEnd", "El H/E End",  nmubins, 0., 0.05);
	fElHistos[8]  = new TH1D("ElSigmaIetaIetaEnd", "El sigmaIetaIeta End",  nmubins, 0., 0.05);
	fElHistos[9]  = new TH1D("ElDeltaPhiSeedClusterAtCaloEnd", "El DeltaPhiSeedClusterAtCalo End",  nmubins, -0.1, 0.1);
	fElHistos[10] = new TH1D("ElDeltaEtaSeedClusterAtCaloEnd", "El DeltaEtaSeedClusterAtCalo End",  nmubins, -0.05, 0.05);
	fElHistos[11] = new TH1D("ElDeltaPhiSuperClusterAtVtxEnd", "El DeltaPhiSuperClusterAtVtx End",  nmubins, -0.1, 0.1);
	fElHistos[12] = new TH1D("ElDeltaEtaSuperClusterAtVtxEnd", "El DeltaEtaSuperClusterAtVtx End",  nmubins, -0.05, 0.05);
	fElHistos[13] = new TH1D("ElESuperClusterOverPEnd", "El ESuperClusterOverP End",  nmubins, 0., 5.);

	const unsigned int nevbins = 100;
	fMETHistos[0] = new TH1D("METR12", "MET R12", nevbins, 0, 6.283);
	fMETHistos[1] = new TH1D("METR21", "MET R21", nevbins, 0, 6.283);
	fMETHistos[3] = new TH1D("EvtEmFrac", "Evt Em Frac", nevbins, 0, 1.05);
	fMETHistos[4] = new TH1D("EvtChFrac", "Evt Ch Frac", nevbins, 0, 1.05);
	fMETDphi12    = new TH2D("METDphi12", "MET Dphi 12", nevbins, 0, 3.142, nevbins, 0, 3.142);
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

	for(size_t i = 0; i < fTR->NEles; ++i){
		if (fTR->ElEta[i] < 1.479) {
			fElHistos[0]->Fill(fTR->ElHcalOverEcal[i]);
			fElHistos[1]->Fill(fTR->ElSigmaIetaIeta[i]);
			fElHistos[2]->Fill(fTR->ElDeltaPhiSeedClusterAtCalo[i]);
			fElHistos[3]->Fill(fTR->ElDeltaEtaSeedClusterAtCalo[i]);
			fElHistos[4]->Fill(fTR->ElDeltaPhiSuperClusterAtVtx[i]);
			fElHistos[5]->Fill(fTR->ElDeltaEtaSuperClusterAtVtx[i]);
			fElHistos[6]->Fill(fTR->ElESuperClusterOverP[i]);
		} else {
			fElHistos[7]->Fill(fTR->ElHcalOverEcal[i]);
			fElHistos[8]->Fill(fTR->ElSigmaIetaIeta[i]);
			fElHistos[9]->Fill(fTR->ElDeltaPhiSeedClusterAtCalo[i]);
			fElHistos[10]->Fill(fTR->ElDeltaEtaSeedClusterAtCalo[i]);
			fElHistos[11]->Fill(fTR->ElDeltaPhiSuperClusterAtVtx[i]);
			fElHistos[12]->Fill(fTR->ElDeltaEtaSuperClusterAtVtx[i]);
			fElHistos[13]->Fill(fTR->ElESuperClusterOverP[i]);
		}
	}

	fMETHistos[0]->Fill(fTC->fR12);
	fMETHistos[1]->Fill(fTC->fR21);
	double METPhi = fTR->MuCorrMETphi;
	double MET = fTR->MuCorrMET;
	double METBadJetmin = 50.;
	if (fTR->NJets > 1 && MET > METBadJetmin) {
		double dPhiMJ1 = Util::DeltaPhi(fTR->JPhi[0], METPhi);
		double dPhiMJ2 = Util::DeltaPhi(fTR->JPhi[1], METPhi);
		fMETDphi12->Fill(dPhiMJ1, dPhiMJ2);
	}
	double fracEm, fracCh;
	GetEvtEmChFrac(fracEm, fracCh);
	fMETHistos[3]->Fill(fracEm);
	fMETHistos[4]->Fill(fracCh);
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

	outputdir = fOutputDir + "Cleaning/El/";
	tempstring1 = "El H/E Bar";
	tempstring2 = "ElHcalOverEcalBar";
	canv = new TCanvas(tempstring2, tempstring1 , 0, 0, 900, 700);
	fElHistos[0]->DrawCopy();
	Util::PrintPNG(canv, tempstring2, outputdir);
	Util::PrintEPS(canv, tempstring2, outputdir);

	tempstring1 = "El SigmaIetaIeta Bar";
	tempstring2 = "ElSigmaIetaIetaBar";
	canv = new TCanvas(tempstring2, tempstring1 , 0, 0, 900, 700);
	fElHistos[1]->DrawCopy();
	Util::PrintPNG(canv, tempstring2, outputdir);
	Util::PrintEPS(canv, tempstring2, outputdir);

	tempstring1 = "El DeltaPhiSeedClusterAtCalo Bar";
	tempstring2 = "ElDeltaPhiSeedClusterAtCaloBar";
	canv = new TCanvas(tempstring2, tempstring1 , 0, 0, 900, 700);
	fElHistos[2]->DrawCopy();
	Util::PrintPNG(canv, tempstring2, outputdir);
	Util::PrintEPS(canv, tempstring2, outputdir);

	tempstring1 = "El DeltaEtaSeedClusterAtCalo Bar";
	tempstring2 = "ElDeltaEtaSeedClusterAtCaloBar";
	canv = new TCanvas(tempstring2, tempstring1 , 0, 0, 900, 700);
	fElHistos[3]->DrawCopy();
	Util::PrintPNG(canv, tempstring2, outputdir);
	Util::PrintEPS(canv, tempstring2, outputdir);

	tempstring1 = "El DeltaPhiSuperClusterAtVtx Bar";
	tempstring2 = "ElDeltaPhiSuperClusterAtVtxBar";
	canv = new TCanvas(tempstring2, tempstring1 , 0, 0, 900, 700);
	fElHistos[4]->DrawCopy();
	Util::PrintPNG(canv, tempstring2, outputdir);
	Util::PrintEPS(canv, tempstring2, outputdir);

	tempstring1 = "El DeltaEtaSuperClusterAtVtx Bar";
	tempstring2 = "ElDeltaEtaSuperClusterAtVtxBar";
	canv = new TCanvas(tempstring2, tempstring1 , 0, 0, 900, 700);
	fElHistos[5]->DrawCopy();
	Util::PrintPNG(canv, tempstring2, outputdir);
	Util::PrintEPS(canv, tempstring2, outputdir);

	tempstring1 = "El ESuperClusterOverP Bar";
	tempstring2 = "ElESuperClusterOverPBar";
	canv = new TCanvas(tempstring2, tempstring1 , 0, 0, 900, 700);
	fElHistos[6]->DrawCopy();
	Util::PrintPNG(canv, tempstring2, outputdir);
	Util::PrintEPS(canv, tempstring2, outputdir);

	tempstring1 = "El H/E End";
	tempstring2 = "ElHcalOverEcalEnd";
	canv = new TCanvas(tempstring2, tempstring1 , 0, 0, 900, 700);
	fElHistos[7]->DrawCopy();
	Util::PrintPNG(canv, tempstring2, outputdir);
	Util::PrintEPS(canv, tempstring2, outputdir);

	tempstring1 = "El SigmaIetaIeta End";
	tempstring2 = "ElSigmaIetaIetaEnd";
	canv = new TCanvas(tempstring2, tempstring1 , 0, 0, 900, 700);
	fElHistos[8]->DrawCopy();
	Util::PrintPNG(canv, tempstring2, outputdir);
	Util::PrintEPS(canv, tempstring2, outputdir);

	tempstring1 = "El DeltaPhiSeedClusterAtCalo End";
	tempstring2 = "ElDeltaPhiSeedClusterAtCaloEnd";
	canv = new TCanvas(tempstring2, tempstring1 , 0, 0, 900, 700);
	fElHistos[9]->DrawCopy();
	Util::PrintPNG(canv, tempstring2, outputdir);
	Util::PrintEPS(canv, tempstring2, outputdir);

	tempstring1 = "El DeltaEtaSeedClusterAtCalo End";
	tempstring2 = "ElDeltaEtaSeedClusterAtCaloEnd";
	canv = new TCanvas(tempstring2, tempstring1 , 0, 0, 900, 700);
	fElHistos[10]->DrawCopy();
	Util::PrintPNG(canv, tempstring2, outputdir);
	Util::PrintEPS(canv, tempstring2, outputdir);

	tempstring1 = "El DeltaPhiSuperClusterAtVtx End";
	tempstring2 = "ElDeltaPhiSuperClusterAtVtxEnd";
	canv = new TCanvas(tempstring2, tempstring1 , 0, 0, 900, 700);
	fElHistos[11]->DrawCopy();
	Util::PrintPNG(canv, tempstring2, outputdir);
	Util::PrintEPS(canv, tempstring2, outputdir);

	tempstring1 = "El DeltaEtaSuperClusterAtVtx End";
	tempstring2 = "ElDeltaEtaSuperClusterAtVtxEnd";
	canv = new TCanvas(tempstring2, tempstring1 , 0, 0, 900, 700);
	fElHistos[12]->DrawCopy();
	Util::PrintPNG(canv, tempstring2, outputdir);
	Util::PrintEPS(canv, tempstring2, outputdir);

	tempstring1 = "El ESuperClusterOverP End";
	tempstring2 = "ElESuperClusterOverPEnd";
	canv = new TCanvas(tempstring2, tempstring1 , 0, 0, 900, 700);
	fElHistos[13]->DrawCopy();
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

	tempstring1 = "MET Dphi 12";
	tempstring2 = "METDphi12";
	canv = new TCanvas(tempstring2, tempstring1 , 0, 0, 900, 700);
	fMETDphi12->DrawCopy();
	Util::PrintPNG(canv, tempstring2, outputdir);
	Util::PrintEPS(canv, tempstring2, outputdir);

	tempstring1 = "Evt Em Frac";
	tempstring2 = "EvtEmFrac";
	canv = new TCanvas(tempstring2, tempstring1 , 0, 0, 900, 700);
	fMETHistos[3]->DrawCopy();
	Util::PrintPNG(canv, tempstring2, outputdir);
	Util::PrintEPS(canv, tempstring2, outputdir);

	tempstring1 = "Evt Ch Frac";
	tempstring2 = "EvtChFrac";
	canv = new TCanvas(tempstring2, tempstring1 , 0, 0, 900, 700);
	fMETHistos[4]->DrawCopy();
	Util::PrintPNG(canv, tempstring2, outputdir);
	Util::PrintEPS(canv, tempstring2, outputdir);

}

void PhysQCAnalysis::PrintInfoStart(int nEntries){
// output the text file with run info

	ofstream file;
	file.open(fOutputDir + "/info.txt");
	time_t tim;
	time(&tim);
//	cout << ctime(&tim) << endl;
	file << ctime(&tim) << endl;
	file << nEntries << endl;

	file.close();
	cout << " Info file made in " << fOutputDir << endl;

}

void PhysQCAnalysis::GetEvtEmChFrac(double & fracEm, double & fracCh){
// Computes the event EM and Charged fractions

	double pt_track = 0.;
	double et_em = 0.;
	double et_had = 0.;
	for( int i = 0; i < fTR->NMus; ++i ){
		if(fTR->MuGood[i] != 0) continue;
		pt_track += fTR->MuPt[i];
	}
	for( int i = 0; i < fTR->NEles; ++i ){
		if(fTR->ElGood[i] != 0) continue;
		pt_track += fTR->ElPt[i];
		et_em += fTR->ElEt[i];
	}
	for( int i = 0; i < fTR->NJets; ++i ){
		if(fTR->JGood[i] != 0) continue;
		pt_track += fTR->JChfrac[i] * fTR->JPt[i];
		et_em    += fTR->JEMfrac[i] * fTR->JEt[i];
		et_had   += (1.-fTR->JEMfrac[i]) * fTR->JEt[i];
	}

	fracCh = 0.;
	fracEm = 0.;
	if( et_em + et_had <= 0. ){
		if( fTR->NMus < 1 ) return; // bad prim vtx will trigger this...
		fracCh = 1.;
		fracEm = 1.;
	} else {
		fracCh = pt_track / (et_em + et_had);
		fracEm = et_em / (et_em + et_had);
	}
	return;

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
