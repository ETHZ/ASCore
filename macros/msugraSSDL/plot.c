{
// get LO cross sections out of file
TFile * res_ = new TFile("results_HT450_MET120_PT20_10.root", "READ", "res_");
res_->cd();
TH2D * xsec_nlo_ = (TH2D *) res_->Get("nlo_xsec");
TH2D * yield_ = (TH2D *) res_->Get("msugra_yield");

// smoothing the yields (physical but ugly)
int nYieldSmooths(0);
// smoothing the limit (unphysical but pretty)
int nContSmooths(0);

while (nYieldSmooths > 0) {
	yield_   -> Smooth();
	nYieldSmooths--;
}

int m0_bins(150) , m0_min(0)  , m0_max(3000);
int m12_bins(50) , m12_min(0) , m12_max(1000);

// the following limits are from the LandS tool, insert your numbers here
float observed = 6.07407;
float expected = 6.32644;
float exphi    = 3.74196;
float explo    = 13.5196;


TH2D * obs_contour   = new TH2D ("obs_exclusion"   , "obs_exclusion"   , m0_bins , m0_min, m0_max , m12_bins, m12_min, m12_max);
TH2D * exp_contour   = new TH2D ("exp_exclusion"   , "exp_exclusion"   , m0_bins , m0_min, m0_max , m12_bins, m12_min, m12_max);
TH2D * explo_contour = new TH2D ("explo_exclusion" , "explo_exclusion" , m0_bins , m0_min, m0_max , m12_bins, m12_min, m12_max);
TH2D * exphi_contour = new TH2D ("exphi_exclusion" , "exphi_exclusion" , m0_bins , m0_min, m0_max , m12_bins, m12_min, m12_max);

for (int m0 = 20; m0 <= m0_max; m0+=20){
	for (int m12 = 20; m12 <= m12_max; m12+=20){
		int bin = yield_->FindBin(m0, m12);
		float yield = yield_->GetBinContent(bin);
		if (yield > observed) {
			obs_contour->Fill(m0, m12, 1);
		}
		if (yield > expected) {
			exp_contour->Fill(m0, m12, 1);
		}
		if (yield > exphi) {
			exphi_contour->Fill(m0, m12, 1);
		}
		if (yield > explo) {
			explo_contour->Fill(m0, m12, 1);
		}
	}
}

TFile * outfile = new TFile("plots_crossCheck.root", "RECREATE", "");
obs_contour->Write();
exp_contour->Write();
explo_contour->Write();
exphi_contour->Write();


// ----------- COSMETICS --------------
// filling the holes in the plots
// ------------------------------------
for (int m0 = 20; m0 <= m0_max; m0+=20){
	for (int m12 = 20; m12 <= m12_max; m12+=20){
		int bin = obs_contour->FindBin(m0-10, m12-10);
		if ( (m12 < 200 && m0<1000) || (m0 < 600 && m12 < 300) || (m12 < 220) || (m0 > 2200 && m12 < 220) || (m0 > 2500 && m12 < 260)) {
			obs_contour->SetBinContent  (bin , 1);
			exp_contour->SetBinContent  (bin , 1);
			explo_contour->SetBinContent(bin , 1);
			exphi_contour->SetBinContent(bin , 1);
		}
	}
}


while (nContSmooths > 0) {
	obs_contour   -> Smooth();
	exp_contour   -> Smooth();
	explo_contour -> Smooth();
	exphi_contour -> Smooth();
	nContSmooths--;
}

obs_contour->Draw("cont list");
gPad->Update();
TList * obj   = (TList *) gROOT.GetListOfSpecials().FindObject("contours");
TList * list  = (TList *) obj->At(10);
TGraph * grph = (TGraph*) list->First()->Clone();

exp_contour->Draw("cont list");
gPad->Update();
TList * obj_exp   = (TList *) gROOT.GetListOfSpecials().FindObject("contours");
TList * list_exp  = (TList *) obj_exp->At(10);
TGraph * grph_exp = (TGraph*) list_exp->First()->Clone();

explo_contour->Draw("cont list");
gPad->Update();
TList * obj_explo   = (TList *) gROOT.GetListOfSpecials().FindObject("contours");
TList * list_explo  = (TList *) obj_explo->At(10);
TGraph * grph_explo = (TGraph*) list_explo->First()->Clone();

exphi_contour->Draw("cont list");
gPad->Update();
TList * obj_exphi   = (TList *) gROOT.GetListOfSpecials().FindObject("contours");
TList * list_exphi  = (TList *) obj_exphi->At(10);
TGraph * grph_exphi = (TGraph*) list_exphi->First()->Clone();

TFile * template = new TFile("msugraTanb10Template.root", "READ", "template");
TCanvas * canv = (TCanvas *) template->Get("GridCanvas");
template->Close();
canv->Draw();
grph->SetLineColor(2);
grph->SetLineWidth(3);

grph_exp->SetLineColor(1);
grph_exp->SetLineWidth(3);

grph_explo->SetLineColor(1);
grph_explo->SetLineWidth(2);
grph_explo->SetLineStyle(2);

grph_exphi->SetLineColor(1);
grph_exphi->SetLineWidth(2);
grph_exphi->SetLineStyle(2);

//int n;
//grph_exphi->GetN()>grph_explo->GetN() ? n=grph_explo->GetN() : n=grph_exphi->GetN();
//TGraph * shade = new TGraph(2*n);
//for (i=0;i<n;i++) {
//	shade->SetPoint(i   , grph_explo->GetX()[i]     , grph_explo->GetY()[i]);
//	shade->SetPoint(n+i , grph_explo->GetX()[n-i-1] , grph_exphi->GetY()[n-i-1]);
//}
//
//shade->SetFillStyle(3003);
//shade->SetFillColor(39);
//shade->Draw("f");

grph->Draw("same");
grph_exp->Draw("same");
grph_explo->Draw("same");
grph_exphi->Draw("same");

//file_->cd();
canv->SaveAs("plots/exclusionPlot_crossCheck_withHiggsCombination.pdf", "");

}
