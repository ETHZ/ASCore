#include "SignificanceAnalysis.hh"
#include "base/TreeReader.hh"
#include "helper/Utilities.hh"

using namespace std;

SignificanceAnalysis::SignificanceAnalysis(TreeReader *tr) : UserAnalysisBase(tr){
}

SignificanceAnalysis::~SignificanceAnalysis(){
}

void SignificanceAnalysis::Begin(const char* filename){
	fSignHistsFile = new TFile(fOutputDir + TString(filename), "RECREATE");
	fSignHistsFile->cd();
	fNBinsEta[0] = 10;
	fNBinsEta[1] = 10;
	fNBinsEta[2] = 10;
	fNBinsEta[3] = 20;
	fNBinsEta[4] = 20;
	fNBinsPhi = 20;
	const float etamin[5] = {-3.0, -2.5, -5.0, -10.0, -10.0};
	const float etamax[5] = { 3.0,  2.5,  5.0,  10.0,  10.0};
	const float phimin = -3.1416, phimax = 3.1416;
	TString partname[5];
	partname[0] = "electrons";
	partname[1] = "muons";
	partname[2] = "jets";
	partname[3] = "ecalMET";
	partname[4] = "hcalMET";

	for(int i=0; i<5; i++){
		fH_ptdev[i]    = new TH2D(Form("h_ptdev%d", i),    Form("pT deviation in eta vs phi for %s", partname[i].Data()),  fNBinsEta[i], etamin[i], etamax[i], fNBinsPhi, phimin, phimax);
		fH_ptsum[i]    = new TH2D(Form("h_ptsum%d", i),    Form("pT sum in eta vs phi for %s", partname[i].Data()),        fNBinsEta[i], etamin[i], etamax[i], fNBinsPhi, phimin, phimax);
		fH_pt2sum[i]   = new TH2D(Form("h_pt2sum%d", i),   Form("pT square sum in eta vs phi for %s", partname[i].Data()), fNBinsEta[i], etamin[i], etamax[i], fNBinsPhi, phimin, phimax);
		fH_ptevt[i]    = new TH2I(Form("h_ptevt%d", i),    Form("Number of ev. in eta vs phi for %s", partname[i].Data()), fNBinsEta[i], etamin[i], etamax[i], fNBinsPhi, phimin, phimax);
		fH_ptavg[i]    = new TH2D(Form("h_ptavg%d", i),    Form("pT average in eta vs phi for %s", partname[i].Data()),    fNBinsEta[i], etamin[i], etamax[i], fNBinsPhi, phimin, phimax);
		fH_ptsumeta[i] = new TH1D(Form("h_ptsumeta%d", i), Form("pT sum in eta for %s", partname[i].Data()),               fNBinsEta[i], etamin[i], etamax[i]);
		fH_ptevteta[i] = new TH1I(Form("h_ptevteta%d", i), Form("Number of ev. in eta for %s", partname[i].Data()),        fNBinsEta[i], etamin[i], etamax[i]);
		fH_ptdev[i]->SetXTitle("#eta");
		fH_ptsum[i]->SetXTitle("#eta");
		fH_ptevt[i]->SetXTitle("#eta");
		fH_ptavg[i]->SetXTitle("#eta");
		fH_ptdev[i]->SetYTitle("#phi");
		fH_ptsum[i]->SetYTitle("#phi");
		fH_ptevt[i]->SetYTitle("#phi");
		fH_ptavg[i]->SetYTitle("#phi");

		fH_ptdev[i]->SetStats(false);
		fH_ptavg[i]->SetStats(false);
		fH_ptsum[i]->SetStats(false);
		fH_ptevt[i]->SetStats(false);
	}
}

void SignificanceAnalysis::Analyze(){
	FillSigHistos(0);
	FillSigHistos(1);
	FillSigHistos(2);
	FillSigHistos(3);
	FillSigHistos(4);
}

void SignificanceAnalysis::FillSigHistos(int part){
// Makes 2D plots of pT
// in eta vs phi for
//  part = 0: electrons
//       = 1: muons
//       = 2: jets (default)
//       = 3: ECAL recoil
//       = 4: HCAL recoil
	fSignHistsFile->cd();
	double EcalEta(0.);
	double HcalEta(0.);

	// electrons
	if( part == 0 ){
		for( int ip = 0; ip < fTR->NEles; ++ ip ){
			//////////////////////////////////
			// Your electron selection here //
			//////////////////////////////////			
			fH_ptsum[part]   ->Fill(fTR->ElEta[ip], fTR->ElPhi[ip], fTR->ElPt[ip]);
			fH_pt2sum[part]  ->Fill(fTR->ElEta[ip], fTR->ElPhi[ip], fTR->ElPt[ip]*fTR->ElPt[ip]);
			fH_ptevt[part]   ->Fill(fTR->ElEta[ip], fTR->ElPhi[ip]);
			fH_ptsumeta[part]->Fill(fTR->ElEta[ip], fTR->ElPt[ip]);
			fH_ptevteta[part]->Fill(fTR->ElEta[ip]);
		}
	}
	// muons
	else if( part == 1 ){
		for( int ip = 0; ip < fTR->NMus; ++ ip ){
			//////////////////////////////////
			// Your muon selection here     //
			//////////////////////////////////
			if(fTR->MuIsGlobalMuon[ip] == 0) continue;
			fH_ptsum[part]   ->Fill(fTR->MuEta[ip], fTR->MuPhi[ip], fTR->MuPt[ip]);
			fH_pt2sum[part]  ->Fill(fTR->MuEta[ip], fTR->MuPhi[ip], fTR->MuPt[ip]*fTR->MuPt[ip]);
			fH_ptevt[part]   ->Fill(fTR->MuEta[ip], fTR->MuPhi[ip]);
			fH_ptsumeta[part]->Fill(fTR->MuEta[ip], fTR->MuPt[ip]);
			fH_ptevteta[part]->Fill(fTR->MuEta[ip]);
		}
	}
	// jets
	else if( part == 2 ){
		for( int ip = 0; ip < fTR->NJets; ++ ip ){
			//////////////////////////////////
			// Your jet selection here      //
			//////////////////////////////////
			fH_ptsum[part]   ->Fill(fTR->JEta[ip], fTR->JPhi[ip], fTR->JPt[ip]);
			fH_pt2sum[part]  ->Fill(fTR->JEta[ip], fTR->JPhi[ip], fTR->JPt[ip]*fTR->JPt[ip]);
			fH_ptevt[part]   ->Fill(fTR->JEta[ip], fTR->JPhi[ip]);
			fH_ptsumeta[part]->Fill(fTR->JEta[ip], fTR->JPt[ip]);
			fH_ptevteta[part]->Fill(fTR->JEta[ip]);
		}
	}
	// ECAL MET
	else if( part == 3 ){
		// this is to protect against NAN
		if( fTR->ECALEsumx != fTR->ECALEsumx) return;
		EcalEta = getEta(fTR->ECALEsumx, fTR->ECALEsumy, fTR->ECALEsumz);
		fH_ptsum[part]   ->Fill(EcalEta, fTR->ECALMETPhi, fTR->ECALMET);
		fH_pt2sum[part]  ->Fill(EcalEta, fTR->ECALMETPhi, fTR->ECALMET*fTR->ECALMET);
		fH_ptevt[part]   ->Fill(EcalEta, fTR->ECALMETPhi);
		fH_ptsumeta[part]->Fill(EcalEta, fTR->ECALMET);
		fH_ptevteta[part]->Fill(EcalEta);
	}
	// HCAL MET
	else if( part == 4 ){
	// this is to protect against NAN
		if( fTR->HCALEsumx != fTR->HCALEsumx ) return;
		HcalEta = getEta(fTR->HCALEsumx, fTR->HCALEsumy, fTR->HCALEsumz);
		fH_ptsum[part]   ->Fill(HcalEta, fTR->HCALMETPhi, fTR->HCALMET);
		fH_pt2sum[part]  ->Fill(HcalEta, fTR->HCALMETPhi, fTR->HCALMET*fTR->HCALMET);
		fH_ptevt[part]   ->Fill(HcalEta, fTR->HCALMETPhi);
		fH_ptsumeta[part]->Fill(HcalEta, fTR->HCALMET);
		fH_ptevteta[part]->Fill(HcalEta);
	}
}

void SignificanceAnalysis::End(){
	fSignHistsFile->cd();
	float zminDev = -4., zmaxDev = 4.;

	// Calculate averages and deviations
	for(size_t part = 0; part < 5; ++part){
		for (int i = 0; i < fNBinsEta[part]; ++i) {
			float ptaverEta = 0.;
			if (fH_ptevteta[part]->GetBinContent(i+1) > 0) {
				ptaverEta = fH_ptsumeta[part]->GetBinContent(i+1) / (float)fH_ptevteta[part]->GetBinContent(i+1);
			}
			for (int j = 0; j < fNBinsPhi; ++j) {
				float ptaver = 0.;
				float ptdev = 0.;
				int nptij = (int)fH_ptevt[part]->GetBinContent(i+1, j+1);
				if (nptij > 0) {
					float ptsumij  = fH_ptsum[part]->GetBinContent(i+1, j+1);
					float pt2sumij = fH_pt2sum[part]->GetBinContent(i+1, j+1);
					ptaver = ptsumij / (float)nptij;
					float pterr = sqrt(pt2sumij - (float)nptij*ptaver*ptaver);
					if (pterr <= 0.) {pterr = 0.1;}
					ptdev = (ptsumij - (float)nptij*ptaverEta) / pterr;
				}
				if (ptdev > zmaxDev) {ptdev = zmaxDev;}
				if (ptdev < zminDev) {ptdev = zminDev;}
				fH_ptdev[part]->SetBinContent(i+1, j+1, ptdev);
				fH_ptavg[part]->SetBinContent(i+1, j+1, ptaver);
			}
		}
	}

	TString pnames[5];
	pnames[0] = "el";
	pnames[1] = "mu";
	pnames[2] = "jet";
	pnames[3] = "eMET";
	pnames[4] = "hMET";

	TString pnamel[5];
	pnamel[0] = "Electrons";
	pnamel[1] = "Muons";
	pnamel[2] = "Jets";
	pnamel[3] = "ECALMET";
	pnamel[4] = "HCALMET";

	fTlat->SetTextColor(kBlack);
	fTlat->SetNDC(kTRUE);
	fTlat->SetTextSize(0.04);

	// TStyle *style = gROOT->GetStyle("ETHStyle");
	// const int ncol1 = 11;
	// int colors1[ncol1]= {10, 16, 5, 28, 29, 8, 4, 9, 45, 46, 2};
	// const int ncol2 = 8;
	// // int colors2[ncol2]= {6, 28, 5, 29, 8, 7, 4, 2};
	// int colors2[ncol2]= {4, 28, 5, 29, 8, 7, 47, 2}; // new colors from luc: 21/10/09
	// const int ncont2 = 9;
	// double contours2[ncont2] = { -4.,-3.,-2.,-1., 0., 1., 2., 3., 4.};

	// Plot the histograms
	TString subdir = "SignificancePlots";
	TCanvas *canv;
	for(size_t i = 0; i < 5; ++i){
		// gStyle->SetPalette(ncol2, colors2);
		TString canvname = "PTDev_" + pnames[i];
		TString canvtitle = "pT Deviation for " + pnamel[i];
		canv = new TCanvas(canvname, canvtitle, 0, 0, 900, 700);
		canv->SetRightMargin(0.15);
		// fH_ptdev[i]->SetContour(ncont2, contours2);
		fH_ptdev[i]->SetMinimum(-4);
		fH_ptdev[i]->SetMaximum(4);
		fH_ptdev[i]->DrawCopy("colz");
		fTlat->DrawLatex(0.11,0.92, canvtitle);
		Util::Print(canv, fTag + "_" + canvname, fOutputDir+subdir);

		// gStyle->SetPalette(ncol1, colors1);

		canvname = "PTSum_" + pnames[i];
		canvtitle = "pT sum for " + pnamel[i];
		canv = new TCanvas(canvname, canvtitle, 0, 0, 900, 700);
		canv->SetRightMargin(0.15);
		fH_ptsum[i]->SetMinimum(0);
		fH_ptsum[i]->DrawCopy("lego2 z");
		fTlat->DrawLatex(0.11,0.92, canvtitle);
		Util::Print(canv, fTag + "_" + canvname, fOutputDir+subdir);

		canvname = "PTEvt_" + pnames[i];
		canvtitle = "Event multiplicity for " + pnamel[i];
		canv = new TCanvas(canvname, canvtitle, 0, 0, 900, 700);
		canv->SetRightMargin(0.15);
		fH_ptevt[i]->SetMinimum(0);
		fH_ptevt[i]->DrawCopy("lego2 z");
		fTlat->DrawLatex(0.11,0.92, canvtitle);
		Util::Print(canv, fTag + "_" + canvname, fOutputDir+subdir);

		canvname = "PTAvg_" + pnames[i];
		canvtitle = "pT Average for " + pnamel[i];
		canv = new TCanvas(canvname, canvtitle, 0, 0, 900, 700);
		canv->SetRightMargin(0.15);
		fH_ptavg[i]->SetMinimum(0);
		fH_ptavg[i]->DrawCopy("lego2 z");
		fTlat->DrawLatex(0.11,0.92, canvtitle);
		Util::Print(canv, fTag + "_" + canvname, fOutputDir+subdir);
	}

	// Write the histograms
	for(size_t i = 0; i < 5; ++i){
		fH_ptdev[i]->Write();
		fH_ptsum[i]->Write();
		fH_ptevt[i]->Write();
		fH_ptavg[i]->Write();
	}
	fSignHistsFile->Close();
}

double SignificanceAnalysis::getEta(double x, double y, double z){
	if(fabs(z) <1.0e-5) return 0;
	double theta = atan(sqrt(x*x+y*y)/z);
	if(theta < 0.) theta = theta + 3.141592654;
	return -log(tan(0.5*theta));
}
