{
//-----------------------------------------------------------------------------------------------------------------------------------
// this shared library which is included below contains the template msugra exculsion plot. 
// it can be obtained by following the instructions here: https://twiki.cern.ch/twiki/bin/view/CMS/SUSYLimitTools2011
// note however that i (marc) adapted the given macro. it didn't work first and also the axes and all that is slightly different
// so i suggest you use the version which is checked into cvs
//-----------------------------------------------------------------------------------------------------------------------------------
gSystem->Load("./ExclusionPlot_C.so");
//-----------------------------------------------------------------------------------------------------------------------------------
// all the files which will be used should be present at runtime. 
// this can be achieved by running the SSDLPlotter with the according functions enabled.
//-----------------------------------------------------------------------------------------------------------------------------------

// get LO cross sections out of file
TFile * res_ = new TFile("res.root", "READ", "res_");
res_->cd();
TH2D * xsec_lo_ = (TH2D *) res_->Get("lo_xsec");

int m0_bins(140), m0_min(220), m0_max(3000);
int m12_bins(46), m12_min(100), m12_max(1000);

// get the calculated limits in the form of TH2Ds
TFile * file_ = new TFile("limit.root", "READ", "file_");
file_->cd();
TH2D * obs_limit_    = (TH2D *) file_->Get("msugra_obslimit");
TH2D * exp_limit_    = (TH2D *) file_->Get("msugra_explimit");
TH2D * explo_limit_  = (TH2D *) file_->Get("msugra_explo");
TH2D * exphi_limit_  = (TH2D *) file_->Get("msugra_exphi");

TH2D * obs_contour   = new TH2D ("obs_exclusion"   , "obs_exclusion"   , m0_bins , m0_min-10 ,m0_max+10 , m12_bins , m12_min-10 , m12_max+10);
TH2D * exp_contour   = new TH2D ("exp_exclusion"   , "exp_exclusion"   , m0_bins , m0_min-10 ,m0_max+10 , m12_bins , m12_min-10 , m12_max+10);
TH2D * explo_contour = new TH2D ("explo_exclusion" , "explo_exclusion" , m0_bins , m0_min-10 ,m0_max+10 , m12_bins , m12_min-10 , m12_max+10);
TH2D * exphi_contour = new TH2D ("exphi_exclusion" , "exphi_exclusion" , m0_bins , m0_min-10 ,m0_max+10 , m12_bins , m12_min-10 , m12_max+10);

for (int m0bin = 1; m0bin<=m0_bins; m0bin++){
	for (int m12bin = 1; m12bin<=m12_bins; m12bin++){
		float lo  = xsec_lo_->GetBinContent(m0bin, m12bin);
		int m0 = m0_min + (m0bin-1)*20;
		int m12 = m12_min + (m12bin-1)*20;

		//if ( m12 > 500 || ( m0 > 1000 && m12 > 400) ) {
		//	if (lo < obs) obs_contour->Fill(m0, m12, -1);
		//	if (lo < exp) exp_contour->Fill(m0, m12, -1);
		//	if (lo < explo) explo_contour->Fill(m0, m12, -1);
		//	if (lo < exphi) exphi_contour->Fill(m0, m12, -1);
		//	continue;
		//}
		//float lo  = 1E9*xsec_lo_->GetBinContent(m0/20., m12/20.);
		//float lo  = xsec_lo_->GetBinContent(m0/20., m12/20.);
// observed limit
		float obs = obs_limit_->GetBinContent(m0bin, m12bin);
		if (lo > obs && obs!=0) obs_contour->Fill(m0, m12, 1);
		if (lo < obs) obs_contour->Fill(m0, m12, -1);
// expected limit
		float exp = exp_limit_->GetBinContent(m0bin, m12bin);
		if (lo > exp && exp!=0) exp_contour->Fill(m0, m12, 1);
		if (lo < exp) exp_contour->Fill(m0, m12, -1);
// expected -1 sigma
		float explo = explo_limit_->GetBinContent(m0bin, m12bin);
		if (lo > explo && explo!=0) explo_contour->Fill(m0, m12, 1);
		if (lo < explo) explo_contour->Fill(m0, m12, -1);
// expected +1 sigma
		float exphi = exphi_limit_->GetBinContent(m0bin, m12bin);
		if (lo > exphi && exphi!=0) exphi_contour->Fill(m0, m12, 1);
		if (lo < exphi) exphi_contour->Fill(m0, m12, -1);
	}
}

TFile * outfile = new TFile("plots.root", "RECREATE", "");
obs_contour->Write();
exp_contour->Write();
explo_contour->Write();
exphi_contour->Write();

// filling the holes in the plots COSMETICS
//for (int m0 = 20; m0<= 2000; m0+=20){
//	for (int m12 = 20; m12<= 760; m12+=20){
//		int m0bin  = m0/20;
//		int m12bin = m12/20;
//		if ( (m12 < 150) || (m12 < 240 && m0 <1000)  || (m12 < 300 && m0 <600)  || (m12 < 200 && m0 <800) || (m12 < 200 && m0 <1300) || (m12 < 250 && m0 > 800 && m0 < 1050) )
//			obs_contour->SetBinContent(m0bin, m12bin, 1);
//		if ( (m12 < 150) || (m12 < 200 && m0 <1050)  || (m12 < 300 && m0 <600) || (m12 < 200 && m0 <1400) )
//			exp_contour->SetBinContent(m0bin, m12bin, 1);
//		if ( (m12 < 150) || (m12 < 250 && m0 <1400)  || (m12 < 320 && m0 <600) )
//			explo_contour->SetBinContent(m0bin, m12bin, 1);
//		if ( (m12 < 150) || (m12 < 200 && m0 <1050)  || (m12 < 300 && m0 <500) )
//			exphi_contour->SetBinContent(m0bin, m12bin, 1);
//	}
//}

//int n_smooths(1);
//
//while (n_smooths > 0) {
//	obs_contour   -> Smooth();
//	exp_contour   -> Smooth();
//	explo_contour -> Smooth();
//	exphi_contour -> Smooth();
//	n_smooths--;
//}
//
//obs_contour->Draw("cont list");
//gPad.Update();
//TList * obj   = (TList *) gROOT.GetListOfSpecials().FindObject("contours");
//TList * list  = (TList *) obj->At(10);
//TGraph * grph = (TGraph*) list->First()->Clone();

//exp_contour->Draw("cont list");
//gPad.Update();
//TList * obj_exp   = (TList *) gROOT.GetListOfSpecials().FindObject("contours");
//TList * list_exp  = (TList *) obj_exp->At(10);
//TGraph * grph_exp = (TGraph*) list_exp->First()->Clone();
//
//explo_contour->Draw("cont list");
//gPad.Update();
//TList * obj_explo   = (TList *) gROOT.GetListOfSpecials().FindObject("contours");
//TList * list_explo  = (TList *) obj_explo->At(10);
//TGraph * grph_explo = (TGraph*) list_explo->First()->Clone();
//
//exphi_contour->Draw("cont list");
//gPad.Update();
//TList * obj_exphi   = (TList *) gROOT.GetListOfSpecials().FindObject("contours");
//TList * list_exphi  = (TList *) obj_exphi->At(10);
//TGraph * grph_exphi = (TGraph*) list_exphi->First()->Clone();

//gStyle->SetOptStat(0);
//TCanvas * canv;
//canv = ExclusionPlot();
//grph->SetLineColor(2);
//grph->SetLineWidth(3);
//
//grph_exp->SetLineColor(1);
//grph_exp->SetLineWidth(3);
//
//grph_explo->SetLineColor(1);
//grph_explo->SetLineWidth(2);
//grph_explo->SetLineStyle(2);
//
//grph_exphi->SetLineColor(1);
//grph_exphi->SetLineWidth(2);
//grph_exphi->SetLineStyle(2);
//
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
//
//grph->Draw("same");
//grph_exp->Draw("same");
//grph_explo->Draw("same");
//grph_exphi->Draw("same");

//canv->Write();

}
