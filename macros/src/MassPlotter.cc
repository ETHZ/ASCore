/*****************************************************************************
*   Small Class to make plots for MassAnalysis                               *
*****************************************************************************/

#include "MassPlotter.hh"

#include "TLorentzVector.h"
#include "TGraphAsymmErrors.h"
#include "helper/Utilities.hh"
#include "helper/MassAnalysisTreeClass.h"

#include "TLatex.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TLine.h"
#include "TArrow.h"
#include "TBox.h"
#include "TMath.h"

#include "TString.h"
#include "TObject.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TH2.h"
#include "TH3.h"
#include "TF1.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TCut.h"
#include "TTreeFormula.h"
#include "TStyle.h"
#include "TRandom.h"
#include "TROOT.h"
#include "TVirtualPad.h"
#include "TClonesArray.h"

#include "TGraph2D.h"
#include "TMatrixD.h"
#include "TMatrixT.h"
#include "TVectorD.h"
#include "TPaveStats.h"
#include <sstream>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <vector>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <map>
#include <time.h> // access to date/time


using namespace std;

//____________________________________________________________________________
MassPlotter::MassPlotter(){
// Default constructor, no samples are set
}

//____________________________________________________________________________
MassPlotter::MassPlotter(TString outputdir){
// Explicit constructor with output directory
	setOutputDir(outputdir);
}

//____________________________________________________________________________
MassPlotter::MassPlotter(TString outputdir, TString outputfile){
// Explicit constructor with output directory and output file
	setOutputDir(outputdir);
	setOutputFile(outputfile);
}

//____________________________________________________________________________
MassPlotter::~MassPlotter(){
	fOutputFile->Close();
	delete fOutputFile;
}

//____________________________________________________________________________
void MassPlotter::init(TString filename){
	if(fVerbose > 0) cout << "------------------------------------" << endl;
	if(fVerbose > 0) cout << "Initializing MassPlotter ... " << endl;
	Util::SetStyle();
	loadSamples(filename);

}

//___________________________________________________________________________
void MassPlotter::makePlots(){

//	double dPhisplit[4]={0., 1.0, 2.0, 3.142};
//	MakeMT2PredictionAndPlots(true,  dPhisplit, 4);  // cleaned, fudgefactor 4
//	MakeMT2PredictionAndPlots(false, dPhisplit, 3);  // not cleaned, fudgefactor 3

	double dPhisplit[4]={0., 1.0, 1.5, 3.142};
//	MakeMT2PredictionAndPlots(true,  dPhisplit, 3);  // cleaned, fudgefactor 2
	MakeMT2PredictionAndPlots(false, dPhisplit, 2.5);  // not cleaned, fudgefactor 2

}
//__________________________________________________________________________

void MassPlotter::MakeMT2PredictionAndPlots(bool cleaned , double dPhisplit[], double fudgefactor){
	TString cleanflag;
	if(cleaned) cleanflag="cleaned";
	else        cleanflag="notcleaned";
	std::ostringstream o0, o1, o2, o3;
	o0 << dPhisplit[0];
	o1 << dPhisplit[1];
	o2 << dPhisplit[2];
	o3 << dPhisplit[3];
	TString dPhirange1 = (TString) o0.str() + "to" +(TString) o1.str();
	TString dPhirange2 = (TString) o2.str() + "to" +(TString) o3.str();
	TString dPhirange  = (TString) o0.str()+(TString) o1.str() + (TString) o2.str() +(TString) o3.str();
	
	std::vector<sample>  QCDSamples;
	std::vector<sample>  Samples_NOTData;
	std::vector<sample>  QCDandDataSamples;
	for(int i=0; i< fSamples.size(); ++i){
		if(fSamples[i].sname=="QCD"){
			QCDSamples.push_back(fSamples[i]);
			QCDandDataSamples.push_back(fSamples[i]);
		}
		if(fSamples[i].type!="data"){
			Samples_NOTData.push_back(fSamples[i]);
		}
		if(fSamples[i].type=="data"){
			QCDandDataSamples.push_back(fSamples[i]);
		}
	}

	double ZerotoPi[4]    ={0., 3.142, 3.142, 3.142};
	double controlsplit[4]={dPhisplit[2], dPhisplit[3], dPhisplit[3], dPhisplit[3]};

//	MakePlot(QCDSamples,        "PseudoJetMT2", gNMT2bins, gMT2bins, cleaned, true, "DPhiMhtMpt", ZerotoPi,     "none", "QCD_"+cleanflag+"_0toPi"             , "none", fudgefactor);
//	MakePlot(QCDSamples,        "PseudoJetMT2", gNMT2bins, gMT2bins, cleaned, true, "DPhiMhtMpt", controlsplit, "none", "QCD_"+cleanflag+"_"+dPhirange2       , "none", fudgefactor);
//	MakePlot(QCDandDataSamples, "PseudoJetMT2", gNMT2bins, gMT2bins, cleaned, true, "DPhiMhtMpt", dPhisplit,    "none", "QCDdata_"+cleanflag+"_"+dPhirange2   , "data", fudgefactor);
	MakePlot(fSamples  ,        "PseudoJetMT2", gNMT2bins, gMT2bins, cleaned, true, "DPhiMhtMpt", ZerotoPi,     "none", "all_"+cleanflag+"_0toPi"             , "none", fudgefactor);
	MakePlot(fSamples  ,        "PseudoJetMT2", gNMT2bins, gMT2bins, cleaned, true, "DPhiMhtMpt", controlsplit, "none", "all_"+cleanflag+"_"+dPhirange2       , "none", fudgefactor);
	MakePlot(fSamples  ,        "PseudoJetMT2", gNMT2bins, gMT2bins, cleaned, true, "DPhiMhtMpt", dPhisplit,    "none", "all_"+cleanflag+"_"+dPhirange1       , "QCD" , fudgefactor);
	MakePlot(fSamples  ,        "PseudoJetMT2", gNMT2bins, gMT2bins, cleaned, true, "DPhiMhtMpt", dPhisplit,    "none", "all_"+cleanflag+"_"+dPhirange1       , "data", fudgefactor);
	MakePlot(fSamples  ,        "PseudoJetMT2", gNMT2bins, gMT2bins, cleaned, true, "DPhiMhtMpt", dPhisplit,    "none", "all_"+cleanflag+"_"+dPhirange1       , "none", fudgefactor);


	ABCD_MT2("DPhiMhtMpt", dPhisplit, "none", gNMT2bins, gMT2bins, "QCD_"+cleanflag+dPhirange ,cleaned,  "QCD",  "mc");



//	MakePlot(Samples_NOTData,        "PseudoJetMT2", gNMT2bins,     gMT2bins,     cleaned, true, "DPhiMhtMpt", ZerotoPi,    "none",   "all_nodata"+cleanflag       , "none", fudgefactor);
//	MakePlot(fSamples  ,        "PseudoJetMT2", gNMT2bins,     gMT2bins,     cleaned, true, "DPhiMhtMpt", ZerotoPi,    "none",   "all_"+cleanflag       , "none", fudgefactor);
//	MakePlot(fSamples  ,        "PseudoJetMT2", gNMT2Leptonbins,     gMT2Leptonbins,     cleaned, true, "DPhiMhtMpt", ZerotoPi,    "none",   "all_"+cleanflag       , "none", fudgefactor);
//	MakePlot(fSamples  ,        "PseudoJetMT2", gNJPtbins,     gJPtbins,     cleaned, true, "DPhiMhtMpt", ZerotoPi,    "none",   "all_"+cleanflag       , "none", fudgefactor);
//	MakePlot(fSamples  ,        "PseudoJetMCT", gNMT2bins,     gMT2bins,     cleaned, true, "DPhiMhtMpt", ZerotoPi,    "none",   "all_"+cleanflag       , "none", fudgefactor);
//	MakePlot(fSamples  ,        "PseudoJetMT2", gNMT2Normbins, gMT2Normbins, cleaned, true, "DPhiMhtMpt", ZerotoPi,    "overHT", "all_overHT"+cleanflag , "none", fudgefactor);



}
//________________________________________________________________________
void MassPlotter::MakePlot(std::vector<sample> Samples, TString branch_name, const int nbins, const double bins[], bool cleaned, bool logflag, 
		           TString cut_branch_name, double ysplit[], TString option, TString version, TString prediction, double factor ){
	
	double lower_cut =  ysplit[0];
	double upper_cut =  ysplit[1];
	THStack h_stack(branch_name, "");
	TH1D*   h_data       = new TH1D("data"  , "", nbins, bins );
	TH1D*   h_mc_sum     = new TH1D("mc_sum", "", nbins, bins );
	TH1D*   h_mc_sum_susy     = new TH1D("mc_sum_susy", "", nbins, bins );
	TH1D*   h_susy  = new TH1D("mc_w_susy", "", nbins, bins);	

	// h_data
	h_data -> Sumw2();
	h_data -> SetMarkerStyle(1);
	h_data -> SetMarkerColor(kRed);
	h_data -> SetLineColor(kRed);
	h_data -> SetStats(false);

	h_mc_sum -> SetFillStyle(3004);
	h_mc_sum -> SetFillColor(kBlack);
	h_mc_sum -> SetStats(0);
	
	h_mc_sum    ->Sumw2();
	h_susy 	    ->Sumw2();
	h_mc_sum_susy -> Sumw2();


	// vector of all histos
	vector<TH1D*> h_samples;


	// legend
	TLegend* Legend1 = new TLegend(.6,.5,.9,.88);
	
	// get prediction
	TH1D* h_prediction    = new TH1D("prediction",    "", nbins, bins);
	if(prediction!="none"){
		if(prediction == "QCD"){	
			h_prediction=GetPrediction(branch_name, gNMT2bins, gMT2bins, cleaned, cut_branch_name, ysplit[2], ysplit[3], false, factor );
		}else if(prediction == "data"){
			h_prediction=GetPrediction(branch_name, gNMT2bins, gMT2bins, cleaned, cut_branch_name, ysplit[2], ysplit[3], true,  factor );
		}
	}
	h_prediction -> SetFillColor(401);
	h_prediction -> SetLineColor(402);
	for(size_t i = 0; i < Samples.size(); ++i){
		h_samples.push_back(new TH1D(branch_name+"_"+Samples[i].name, "", nbins, bins));
		h_samples[i] -> Sumw2();
		h_samples[i] -> SetFillColor(Samples[i].color);
		h_samples[i] -> SetLineColor(Samples[i].color);
		h_samples[i] -> SetMarkerColor(Samples[i].color);
		h_samples[i] -> SetStats(false);
		if(Samples[i].type == "susy" ){
			h_samples[i] -> SetLineColor(kBlack);
			h_samples[i] -> SetLineStyle(kDotted);
		}
		TTree *tree = Samples[i].tree;
		Double_t weight = Samples[i].xsection * Samples[i].kfact * Samples[i].lumi / (Samples[i].nevents);
		if(fVerbose>2) cout << "MakePlot: looping over " << Samples[i].sname << endl;
		if(fVerbose>2) cout << "           sample has weight " << weight << endl; 
		tree->ResetBranchAddresses();	
		Init(tree);
		
		TBranch* branch = fChain->GetBranch(branch_name);
		Double_t yvalue;	
		fChain->SetBranchAddress(branch_name, &yvalue, &branch);

		TBranch* cut_branch = fChain ->GetBranch(cut_branch_name);
		Double_t cut_branch_value;
		fChain ->SetBranchAddress(cut_branch_name, &cut_branch_value, &cut_branch);

		if (fChain == 0) {cout << "fChain ==0" << endl; return;}
		Long64_t nbytes = 0, nb = 0;
		Long64_t nentries = fChain->GetEntriesFast();
		
		for (Long64_t jentry=0; jentry<nentries;jentry++) {
			Long64_t ientry = LoadTree(jentry);
			nb = fChain->GetEntry(jentry);   nbytes += nb;
				
			// cleaning & cuts
			if(cut_branch_value > upper_cut) continue;
			if(cut_branch_value < lower_cut) continue;
			if(LeptConfig !=LeptConfigCut  ) continue; //lepton veto
			if(Vectorsumpt > VectorSumPtCut) continue;
			if( NJets < NJetsCut )           continue;
			if( NJets > MaxNJetsCut )        continue;
			if(cleaned){
				if(R12R21 == 0)      continue;
			}	
			if(branch_name=="PseudoJetMT2"){
				if(yvalue      < 0  ) continue;
			} else if(PseudoJetMT2 < 0  ) continue;      
			
			// fill histo
			if(option == "none")  {
				h_samples[i]                ->Fill(yvalue     , weight);
				if(Samples[i].type == "mc"  &&  Samples[i].sname == "QCD" && yvalue < 50){h_prediction -> Fill(yvalue     , weight);}
			}
			if(option == "overHT")  {
				h_samples[i]                ->Fill(yvalue/HT     , weight);
			}
		}
	}
	if(prediction!="none"){
		h_stack.Add(h_prediction);
		Legend1 ->AddEntry(h_prediction, "data-driven QCD", "f");
		h_mc_sum->Add(h_prediction);
	}
	for(int i=0; i<Samples.size(); ++i){
		if(Samples[i].type=="mc" || Samples[i].type=="susy") h_mc_sum_susy->Add(h_samples[i]);
		if((Samples[i].sname=="QCD" && prediction=="none")){
			h_samples[i] -> SetMinimum(0.01);
			h_stack.Add(h_samples[i]);
			h_mc_sum->Add(h_samples[i]);
			Legend1 ->AddEntry(h_samples[i], Samples[i].sname, "f");
		}
		if(Samples[i].sname!="QCD" && Samples[i].type!="data"){
			h_stack.Add(h_samples[i]);
			if(Samples[i].type!="susy") h_mc_sum->Add(h_samples[i]);
			if(Samples[i].type=="susy") h_susy  ->Add(h_samples[i]);
			if(Samples[i].name != "PhotonJets_Pt40to100-V01-03-01" && Samples[i].name != "PhotonJets_Pt100to200-V01-03-01"){
				Legend1 ->AddEntry(h_samples[i], Samples[i].sname, "f");
			}
		}
		if(Samples[i].type == "data"){
			h_data -> Add(h_samples[i]);
		}

	}
	if(h_data->Integral()>0){
		Legend1     ->AddEntry(h_data, "data", "l");
	}


	TString xtitle = branch_name + " (GeV)";
//	TString xtitle = "leading JPT (GeV)";
	if(option =="overHT") xtitle = branch_name + " / HT (GeV)";
	TString ytitle = "events";

	// print histo without QCD prediction
	printHisto(h_stack, h_data,   Legend1 , branch_name+"_"+version, "hist", logflag, xtitle, ytitle);
	plotRatioStack(h_stack,  h_mc_sum, h_data,  true, false, branch_name+"_"+version+"ratio", Legend1, xtitle, ytitle);


	if(prediction!="none"){
		// print histo with SM prediction
		printHisto(h_stack, h_data, h_mc_sum, Legend1 , branch_name+"_"+version+"_prediction_from_"+prediction, "hist", logflag, xtitle, ytitle);
	}
	

	// _____________________________________________________
	double data_error      =0;
	double mc_sum_error    =0;
	double susy_error      = 0;
	double sm_w_susy_error = 0;
	for(int bin=16; bin <= h_mc_sum->GetNbinsX(); ++bin){
		data_error      += h_data   ->GetBinError(bin)*h_data  ->GetBinError(bin);
		mc_sum_error    += h_mc_sum ->GetBinError(bin)*h_mc_sum->GetBinError(bin);
		susy_error      += h_susy   ->GetBinError(bin) * h_susy ->GetBinError(bin);
		sm_w_susy_error += h_mc_sum_susy   ->GetBinError(bin) * h_mc_sum_susy->GetBinError(bin);
	}
	data_error       = sqrt(data_error);
	mc_sum_error     = sqrt(mc_sum_error);
	susy_error       = sqrt(susy_error);
	sm_w_susy_error  = sqrt(sm_w_susy_error);

	cout << "____________________________________________________" << endl;
	cout << "Prediction for QCD based on "   << prediction   << " cleaned "   << cleaned        << endl;
	cout << " " << cut_branch_name << " in " << ysplit[0]   <<  " and " << ysplit[1] << endl; 
	cout << "  # events with MT2 > 135:                          " << endl;
	cout << "  MC data driven:           " << h_mc_sum->Integral(16, 100)    <<  " plusminus " << mc_sum_error    << endl;
	cout << "  data:                     " << h_data  ->Integral(16, 100)    <<  " plusminus " << data_error      << endl;
	cout << "  SUSY:                     " << h_susy->Integral(16, 100)      <<  " plusminus " << susy_error      << endl; 
	cout << "  SM with SUSY:             " << h_mc_sum_susy->Integral(16, 100)      <<  " plusminus " << sm_w_susy_error      << endl; 
	cout << "  h_mc_sum in bin 16        " << h_mc_sum->GetBinContent(16) <<endl;
	// _________________________________________________________
	for(int i=0; i<h_samples.size(); ++i){
		delete h_samples[i];
	}
	delete h_prediction;
	delete h_mc_sum;
	delete Legend1;
	delete h_data;
	delete h_susy;
	delete h_mc_sum_susy;
}

//_____________________________________________________________________________________
TH1D* MassPlotter::GetPrediction(TString branch_name, const int nbins, const double bins[], bool cleaned, TString cut_branch_name, double lower_cut, double upper_cut, bool data, double factor ){

	TH1D*   h_prediction       = new TH1D("prediction", "", nbins, bins );
	TH1D*   h_fudge            = new TH1D("fudge_factor", "", nbins, bins);
	
	// fill h_fudge
	for(int i=1; i<=nbins+1; ++i ){
		if(h_fudge->GetBinLowEdge(i) <50) continue;
		h_fudge->SetBinContent(i, factor);
		h_fudge->SetBinError(i, 0.5*factor);
	}


	// h_data
	h_prediction -> Sumw2();
	h_prediction -> SetStats(false);

	for(size_t i = 0; i < fSamples.size(); ++i){
		if(data==false && fSamples[i].sname != "QCD") continue;
		if(data==true  && fSamples[i].type  != "data") continue;

		TTree *tree = fSamples[i].tree;
		Double_t weight = fSamples[i].xsection * fSamples[i].kfact * fSamples[i].lumi / (fSamples[i].nevents);
		tree->ResetBranchAddresses();	
		Init(tree);
		
		TBranch* branch = fChain->GetBranch(branch_name);
		Double_t yvalue;	
		fChain->SetBranchAddress(branch_name, &yvalue, &branch);

		TBranch* cut_branch = fChain ->GetBranch(cut_branch_name);
		Double_t cut_branch_value;
		fChain ->SetBranchAddress(cut_branch_name, &cut_branch_value, &cut_branch);

		if (fChain == 0) {cout << "fChain ==0" << endl; return NULL;}
		Long64_t nbytes = 0, nb = 0;
		Long64_t nentries = fChain->GetEntriesFast();
		
		for (Long64_t jentry=0; jentry<nentries;jentry++) {
			Long64_t ientry = LoadTree(jentry);
			nb = fChain->GetEntry(jentry);   nbytes += nb;
				
			// cleaning & cuts
			if(cut_branch_value > upper_cut) continue;
			if(cut_branch_value < lower_cut) continue;
			if(LeptConfig != LeptConfigCut ) continue; //lepton veto
			if(Vectorsumpt > VectorSumPtCut) continue;
			if( NJets < NJetsCut )            continue;
			if( NJets > MaxNJetsCut )        continue;
			if(cleaned){
				if(R12R21 == 0)      continue;
			}	
			if(branch_name=="PseudoJetMT2"){
				if(yvalue      < 0  ) continue;
			} else if(PseudoJetMT2 < 0  ) continue;      
			// fill histo
			if(yvalue > 50) h_prediction     ->Fill(yvalue     , weight);
		}
	}
	h_prediction->Multiply(h_fudge);
	return h_prediction;	
}

//_______________________________________________________________________
void MassPlotter::FillRatioHistHT(TString version, bool cleaned){

	TH1D *h_HT_ratio_b    = new TH1D("HT_ratio_b_"+version    , "", gNHTbins, gHTbins);
	TH1D *h_HT_lower_b    = new TH1D("HT_lower_y_band_b"      , "", gNHTbins, gHTbins);
	TH1D *h_HT_upper_b    = new TH1D("HT_upper_y_band_b"      , "", gNHTbins, gHTbins);
	
	h_HT_ratio_b   ->Sumw2();
	h_HT_lower_b   ->Sumw2();
	h_HT_upper_b   ->Sumw2();

	for(size_t i = 0; i < fSamples.size(); ++i){
		if(fSamples[i].type != "mc") continue;
		if(fSamples[i].sname!= "QCD") continue;
		
		TTree *tree = fSamples[i].tree;
		Double_t weight = fSamples[i].xsection * fSamples[i].kfact * fSamples[i].lumi / (fSamples[i].nevents);
		if(fVerbose>2) cout << "FillRatioHist: looping over " << fSamples[i].sname << endl;
		if(fVerbose>2) cout << "               sample has weight " << weight << endl; 
		
		tree->ResetBranchAddresses();
		Init(tree);
		if (fChain == 0) {cout << "fChain ==0" << endl; return;}
		Long64_t nbytes = 0, nb = 0;
		Long64_t nentries = fChain->GetEntriesFast();
		
		for (Long64_t jentry=0; jentry<nentries;jentry++) {
			Long64_t ientry = LoadTree(jentry);
			nb = fChain->GetEntry(jentry);   nbytes += nb;
		
			// cuts & selection
			if( Vectorsumpt > VectorSumPtCut) continue;
			if( NJets < NJetsCut )            continue;
			if( NJets > MaxNJetsCut )        continue;
			if( LeptConfig !=LeptConfigCut  ) continue;
			if(cleaned){
				if( R12R21 ==0      ) continue;
			}
			// fill histo
			if(PseudoJetMT2/HT < 0.3 ){h_HT_lower_b->Fill(HT);}
			if(PseudoJetMT2/HT > 0.3 ){h_HT_upper_b->Fill(HT);}
		}
	}
	h_HT_ratio_b->Divide( h_HT_upper_b, h_HT_lower_b);	

	TCanvas* c1 = new TCanvas("c1","", 20,100,1000,700);
	
	gPad->SetLogy(1);
	h_HT_ratio_b->Draw("E1");
	Util::PrintNoEPS(c1, "HT_ratio_"+version, fOutputDir, fOutputFile);


	delete h_HT_ratio_b;
	delete h_HT_lower_b;
	delete h_HT_upper_b;
	delete c1;
} 
//_______________________________________________________________________
void MassPlotter::FillRatioHistdPhi(TString cut_branch_name, double lower_cut, double middle_cut, double upper_cut, TString version, bool cleaned){

	TH1D *h_dPhi_ratio_b    = new TH1D("dPhi_ratio_b_"+version    , "", gNDPhibins, gDPhibins);
	TH1D *h_dPhi_lower_b    = new TH1D("dPhi_lower_y_band_b"      , "", gNDPhibins, gDPhibins);
	TH1D *h_dPhi_upper_b    = new TH1D("dPhi_upper_y_band_b"      , "", gNDPhibins, gDPhibins);
	
	TH1D *h_dPhi_ratio_s    = new TH1D("dPhi_ratio_s_"+version    , "", gNDPhibins, gDPhibins);
	TH1D *h_dPhi_lower_s    = new TH1D("dPhi_lower_y_band_s"      , "", gNDPhibins, gDPhibins);
	TH1D *h_dPhi_upper_s    = new TH1D("dPhi_upper_y_band_s"      , "", gNDPhibins, gDPhibins);

	h_dPhi_ratio_b   ->Sumw2();
	h_dPhi_lower_b   ->Sumw2();
	h_dPhi_upper_b   ->Sumw2();

	h_dPhi_ratio_s   ->Sumw2();
	h_dPhi_lower_s   ->Sumw2();
	h_dPhi_upper_s   ->Sumw2();

	for(size_t i = 0; i < fSamples.size(); ++i){
		if(fSamples[i].type != "mc") continue;
		if(fSamples[i].sname!= "QCD") continue;
		
		TTree *tree = fSamples[i].tree;
		Double_t weight = fSamples[i].xsection * fSamples[i].kfact * fSamples[i].lumi / (fSamples[i].nevents);
		if(fVerbose>2) cout << "FillRatioHist: looping over " << fSamples[i].sname << endl;
		if(fVerbose>2) cout << "               sample has weight " << weight << endl; 
		
		tree->ResetBranchAddresses();
		Init(tree);
		if (fChain == 0) {cout << "fChain ==0" << endl; return;}
		Long64_t nbytes = 0, nb = 0;
		Long64_t nentries = fChain->GetEntriesFast();
		
		TBranch* cut_branch = fChain ->GetBranch(cut_branch_name);
		Double_t cut_branch_value;
		fChain ->SetBranchAddress(cut_branch_name, &cut_branch_value, &cut_branch);
		
		for (Long64_t jentry=0; jentry<nentries;jentry++) {
			Long64_t ientry = LoadTree(jentry);
			nb = fChain->GetEntry(jentry);   nbytes += nb;
		

			// cuts & selection
			if(cleaned){
		//		if( R12R21 ==0      ) continue;
				if( Vectorsumpt > VectorSumPtCut) continue;
				if( NJets < NJetsCut )            continue;
				if( NJets > MaxNJetsCut )        continue;
				if( LeptConfig !=LeptConfigCut  ) continue;
				if( MPT < 5         ) continue;
				if( MHT < 10       )  continue;
			}
			// fill histo
			if ( DPhiMhtMpt < 1.0             && cut_branch_value > middle_cut && cut_branch_value < upper_cut ){
				h_dPhi_lower_s -> Fill(DPhiMhtMpt, weight);
			}
			if ( DPhiMhtMpt < 1.0             && cut_branch_value > lower_cut && cut_branch_value < middle_cut){
				h_dPhi_lower_b -> Fill(DPhiMhtMpt, weight);
			}
			if(  DPhiMhtMpt > (TMath::Pi()-1) && cut_branch_value > middle_cut && cut_branch_value < upper_cut){
				h_dPhi_upper_s -> Fill(DPhiMhtMpt-2*(DPhiMhtMpt-(TMath::Pi()/2.)), weight);
			}
			if(  DPhiMhtMpt > (TMath::Pi()-1) && cut_branch_value > lower_cut && cut_branch_value < middle_cut){
				h_dPhi_upper_b -> Fill(DPhiMhtMpt-2*(DPhiMhtMpt-(TMath::Pi()/2.)), weight);
			}
		}
	}
	h_dPhi_ratio_b->Divide(h_dPhi_lower_b, h_dPhi_upper_b);	
	h_dPhi_ratio_s->Divide(h_dPhi_lower_s, h_dPhi_upper_s);

	TCanvas* c1 = new TCanvas("c1","", 20,100,1000,700);
	
	h_dPhi_ratio_s->Draw("E1");
	h_dPhi_ratio_b->Draw("same,E1");
	h_dPhi_ratio_s->SetLineColor(kRed);

	Util::PrintNoEPS(c1, "dPhi_ratio_"+version, fOutputDir, fOutputFile);


	delete h_dPhi_ratio_s;
	delete h_dPhi_lower_s;
	delete h_dPhi_upper_s;
	delete h_dPhi_ratio_b;
	delete h_dPhi_lower_b;
	delete h_dPhi_upper_b;
	delete c1;
} 


//_________________________________________________________________________
void MassPlotter::ABCD_MT2(TString branch_name, double ysplit[], TString option, const int nbins, const double bins[], bool cleaned, TString sname, TString type){
	ABCD_MT2(branch_name, ysplit, option, nbins, bins, "1", cleaned, sname, type);
}

//__________________________________________________________________________
void MassPlotter::ABCD_MT2(TString branch_name, double ysplit[], TString option, const int nbins, const double bins[],  TString version, bool cleaned, TString sname, TString type){
	TH2D *h_ABCD_MT2            = new TH2D("ABCD_MT2_"+branch_name+"_"+version           , "",  nbins, bins, 100, ysplit[0], ysplit[3]);
	TH1D *h_ABCD_lower_y_band   = new TH1D("ABCD_lower_y_band_"+branch_name+"_"+version  , "",  nbins, bins);
	TH1D *h_ABCD_upper_y_band   = new TH1D("ABCD_upper_y_band_"+branch_name+"_"+version  , "",  nbins, bins);

	h_ABCD_MT2          ->Sumw2();
	h_ABCD_lower_y_band ->Sumw2();
	h_ABCD_upper_y_band ->Sumw2();

	for(size_t i = 0; i < fSamples.size(); ++i){
		if(type  !="none"  && fSamples[i].type != type  ) continue;
		if(sname !="none"  && fSamples[i].sname!= sname ) continue;

		TTree *tree = fSamples[i].tree;
		Double_t weight = fSamples[i].xsection * fSamples[i].kfact * fSamples[i].lumi / (fSamples[i].nevents);
		if(fVerbose>2) cout << "ABCD: looping over " << fSamples[i].sname << endl;
		if(fVerbose>2) cout << "      sample has weight " << weight << endl; 
		
		tree->ResetBranchAddresses();
		Init(tree);
		if (fChain == 0) {cout << "fChain ==0" << endl; return;}
		Long64_t nbytes = 0, nb = 0;
		Long64_t nentries = fChain->GetEntriesFast();
		
		TBranch* branch = fChain->GetBranch(branch_name);
		Double_t yvalue;	
		fChain->SetBranchAddress(branch_name, &yvalue, &branch);
		
		for (Long64_t jentry=0; jentry<nentries;jentry++) {
			Long64_t ientry = LoadTree(jentry);
			nb = fChain->GetEntry(jentry);   nbytes += nb;
			
			// cuts & selection
			if( PseudoJetMT2 <0 ) continue;
			if( Vectorsumpt > VectorSumPtCut) continue;
			if( NJets < NJetsCut )            continue;
			if( NJets > MaxNJetsCut )        continue;
			if( LeptConfig !=LeptConfigCut  ) continue;
			if( cleaned){
				if( R12R21 ==0      ) continue;
			}
			// fill histo
			double xvalue=-999.99;
			if      (option == "none")  {xvalue = PseudoJetMT2;}
			else if (option == "overHT"){xvalue = PseudoJetMT2/HT;} 
			if(  yvalue > ysplit[0] && yvalue < ysplit[1] ){
				h_ABCD_lower_y_band -> Fill(xvalue, weight);
				h_ABCD_MT2          -> Fill(xvalue, yvalue, weight);
			}
			if(  yvalue > ysplit[2] && yvalue < ysplit[3] ){
				h_ABCD_upper_y_band -> Fill(xvalue, weight);
				h_ABCD_MT2          -> Fill(xvalue, yvalue, weight);
			}

		}
	}

		
	h_ABCD_upper_y_band->SetLineColor(kBlack);
	h_ABCD_upper_y_band->SetLineStyle(kBlack);	

	h_ABCD_lower_y_band->SetMarkerStyle(1);
	h_ABCD_lower_y_band->SetLineColor(kRed);
	h_ABCD_lower_y_band->SetLineStyle(kDashed);
	h_ABCD_lower_y_band->SetMarkerColor(kRed);
	
	// title and legend for ratio plot
 	TLegend* Legend1 = new TLegend(.55,.6,.9,.85);
	std::ostringstream o0, o1, o2, o3;
	o0 << ysplit[0];
	o1 << ysplit[1];
	o2 << ysplit[2];
	o3 << ysplit[3];
	TString leg1 = "MT2 for DPhi(MPT, MHT) in ("+(TString) o0.str()+ ", " + (TString) o1.str() + ")"; 
	TString leg2 = "MT2 for DPhi(MPT, MHT) in ("+(TString) o2.str()+ ", " + (TString) o3.str() + ")"; 

	Legend1 -> AddEntry(h_ABCD_lower_y_band, leg1 , "l");	
	Legend1 -> AddEntry(h_ABCD_upper_y_band, leg2 , "l");	
	
	TString xtitle = "MT2 (GeV)";
	TString ytitle = "events (normalized)";

	// make plots
	printHisto(h_ABCD_MT2         ,   "MT2_pseudojet_vs_"+branch_name+"_"+version                , "colz", false);
	printHisto(h_ABCD_lower_y_band,   "MT2_pseudojet_vs_"+branch_name+"_"+version+"_lower_y_band", "hist", true);
	printHisto(h_ABCD_upper_y_band,   "MT2_pseudojet_vs_"+branch_name+"_"+version+"_upper_y_band", "hist", true);

	plotRatio(h_ABCD_lower_y_band, h_ABCD_upper_y_band, true, false, "MT2_pseudojet_vs_"+branch_name+"_"+version, Legend1, xtitle, ytitle);

	delete h_ABCD_MT2;
	delete h_ABCD_lower_y_band;
	delete h_ABCD_upper_y_band;

}


//_________________________________________________________________________________
void MassPlotter::plotRatioStack(THStack hstack, TH1* h1_orig, TH1* h2_orig, bool logflag, bool normalize, TString name, TLegend* leg, TString xtitle, TString ytitle){
	// define canvas and pads 
	TH1D *h1 = (TH1D*)h1_orig->Clone("h1_copy");
	TH1D *h2 = (TH1D*)h2_orig->Clone("h2_copy");

	h1->SetStats(0);	
	h2->SetStats(0);	
	
	TCanvas* c1 = new TCanvas("c1","", 20,100,1000,700);
	c1 -> cd();
	
	float border = 0.3;
 	float scale = (1-border)/border;
 
 	TPad *p_plot  = new TPad("plotpad",  "Pad containing the overlay plot", 0.00, border, 1.00, 1.00, 0, 0);
 	p_plot->SetBottomMargin(0);
 	p_plot->Draw();
 	TPad *p_ratio = new TPad("ratiopad", "Pad containing the ratio",        0.00, 0.00, 1.00, border, 0, 0);
 	p_ratio->SetTopMargin(0);
 	p_ratio->SetBottomMargin(0.35);
 	p_ratio->Draw();
 
	// draw overlay plot
 	p_plot ->cd();

	if(logflag) gPad->SetLogy(1);
	gPad->SetFillStyle(0);
		
	// Scaling
	if(normalize){
		h1->Scale(1.0/h1->Integral());
		h2->Scale(1.0/h2->Integral());
	}
	
	// Determine plotting range
	double max1 = h1->GetMaximum();
	double max2 = h2->GetMaximum();
	double max  = (max1>max2)?max1:max2;
	if(logflag) max = 5*max, "3";
	else max = 1.05*max;
	h1->SetMaximum(max);
	h2->SetMaximum(max);
//	h1->SetMinimum(0.000000001);
//	h2->SetMinimum(0.000000001);

	hstack.SetMinimum(0.02);
	hstack.Draw("hist");
	h2      ->Draw("same");

	// title
	TLatex lat;
	lat.SetNDC(1);
	lat.SetTextColor(4);
	lat.SetTextSize(0.05);

	// y axis title 
	lat.SetTextAlign(33); 
	lat.SetTextColor(1);
	lat.SetTextAngle(90);
	lat.DrawLatex(0.035, 0.9, ytitle);
		
 	p_plot ->Draw();
	gPad->RedrawAxis();

	if(leg != NULL ){
		leg -> SetFillColor(0);
		leg -> SetBorderSize(0);
		leg -> Draw();
	} 
	
	// draw the ratio plot
 	p_ratio ->cd();
	gPad->SetLogy();
 
	TH1D *h_ratio = (TH1D*)h2_orig->Clone("h2_copy");	

	h_ratio->GetXaxis()->SetTitleSize(scale * h1->GetXaxis()->GetTitleSize());
	h_ratio->GetYaxis()->SetTitleSize(scale * h1->GetYaxis()->GetTitleSize());
	h_ratio->GetXaxis()->SetLabelSize(scale * h1->GetXaxis()->GetLabelSize());
	h_ratio->GetYaxis()->SetLabelSize(scale * h1->GetYaxis()->GetLabelSize());
	h_ratio->GetXaxis()->SetTickLength(scale * h1->GetXaxis()->GetTickLength());
	h_ratio->GetYaxis()->SetTickLength(h1->GetYaxis()->GetTickLength());
	h_ratio->SetLineWidth(2);
	h_ratio->SetFillColor(kBlue);
	h_ratio->SetLineColor(kBlue);
	h_ratio ->SetStats(0);
	h_ratio ->SetMarkerStyle(20);
	h_ratio ->SetMarkerSize(0.1);

 	h_ratio ->Divide(h2, h1);
 	h_ratio ->SetMinimum(0.1);
 	h_ratio ->SetMaximum(8.0);
	h_ratio ->GetYaxis()->SetTitleOffset(h1->GetYaxis()->GetTitleOffset());

	h_ratio ->DrawCopy("E2");
 
	TLine *l3 = new TLine(h1->GetXaxis()->GetXmin(), 1.00, h1->GetXaxis()->GetXmax(), 1.00);
	l3->SetLineWidth(2);
	l3->SetLineStyle(7);
	l3->Draw();
	// ratio title
	lat.SetTextAlign(33);
	lat.SetTextColor(1);
	lat.SetTextSize(0.1);
	lat.SetNDC(1);
	lat.SetTextAngle(90);
	lat.DrawLatex(0.035, 1.0, "data / MC");
	
	// x axis title
	lat.SetTextAngle(0);
	lat.DrawLatex(0.9, 0.15, xtitle);
	gPad->RedrawAxis();
 	p_ratio ->Draw();
 	c1->Update();

	TString save=name+"_ratio";
	Util::Print(c1, save, fOutputDir, fOutputFile);	

	delete h1;
	delete h2;
	delete h_ratio;
	delete p_plot;
	delete p_ratio;
	delete c1;

}
//_________________________________________________________________________________
void MassPlotter::plotRatio(TH1* h1_orig, TH1* h2_orig, bool logflag, bool normalize, TString name, TLegend* leg, TString xtitle, TString ytitle){
	// define canvas and pads 
	TH1D *h1 = (TH1D*)h1_orig->Clone("h1_copy");
	TH1D *h2 = (TH1D*)h2_orig->Clone("h2_copy");

	h1->SetStats(0);	
	h2->SetStats(0);	
	
	TCanvas* c1 = new TCanvas("c1","", 20,100,1000,700);
	c1 -> cd();
	
	float border = 0.4;
 	float scale = (1-border)/border;
 
 	TPad *p_plot  = new TPad("plotpad",  "Pad containing the overlay plot", 0.00, border, 1.00, 1.00, 0, 0);
 	p_plot->SetBottomMargin(0);
 	p_plot->Draw();
 	TPad *p_ratio = new TPad("ratiopad", "Pad containing the ratio",        0.00, 0.00, 1.00, border, 0, 0);
 	p_ratio->SetTopMargin(0);
 	p_ratio->SetBottomMargin(0.35);
 	p_ratio->Draw();
 
	// draw overlay plot
 	p_plot ->cd();

	if(logflag) gPad->SetLogy(1);
	gPad->SetFillStyle(0);
		
	// Scaling
	if(normalize){
		h1->Scale(1.0/h1->Integral());
		h2->Scale(1.0/h2->Integral());
	}
	
	// Determine plotting range
	double max1 = h1->GetMaximum();
	double max2 = h2->GetMaximum();
	double max  = (max1>max2)?max1:max2;
	if(logflag) max = 5*max, "3";
	else max = 1.05*max;
	h1->SetMaximum(max);
	h2->SetMaximum(max);
//	h1->SetMinimum(0.000000001);
//	h2->SetMinimum(0.000000001);
	h2     ->Draw();
	h1     ->Draw("same");
	
	// title
	TLatex lat;
	lat.SetNDC(1);
	lat.SetTextColor(4);
	lat.SetTextSize(0.05);

	// y axis title 
	lat.SetTextAlign(33); 
	lat.SetTextColor(1);
	lat.SetTextAngle(90);
	lat.DrawLatex(0.035, 0.9, ytitle);
		
 	p_plot ->Draw();
	gPad->RedrawAxis();

	if(leg != NULL ){
		leg -> SetFillColor(0);
		leg -> SetBorderSize(0);
		leg -> Draw();
	} 
	
	// draw the ratio plot
 	p_ratio ->cd();
 
	TH1D *h_ratio = (TH1D*)h2_orig->Clone("h2_copy");	

	h_ratio->GetXaxis()->SetTitleSize(scale * h1->GetXaxis()->GetTitleSize());
	h_ratio->GetYaxis()->SetTitleSize(scale * h1->GetYaxis()->GetTitleSize());
	h_ratio->GetXaxis()->SetLabelSize(scale * h1->GetXaxis()->GetLabelSize());
	h_ratio->GetYaxis()->SetLabelSize(scale * h1->GetYaxis()->GetLabelSize());
	h_ratio->GetXaxis()->SetTickLength(scale * h1->GetXaxis()->GetTickLength());
	h_ratio->GetYaxis()->SetTickLength(h1->GetYaxis()->GetTickLength());
	h_ratio->SetLineWidth(2);
	h_ratio->SetFillColor(kBlue);
	h_ratio->SetLineColor(kBlue);
	h_ratio ->SetStats(0);
	h_ratio ->SetMarkerStyle(20);
	h_ratio ->SetMarkerSize(0.1);

 	h_ratio ->Divide(h1, h2);
 	h_ratio ->SetMinimum(0.2);
 	h_ratio ->SetMaximum(15.0);
	h_ratio ->GetYaxis()->SetTitleOffset(h1->GetYaxis()->GetTitleOffset());

	h_ratio ->DrawCopy("E2");
 
	TLine *l3 = new TLine(h1->GetXaxis()->GetXmin(), 1.00, h1->GetXaxis()->GetXmax(), 1.00);
	l3->SetLineWidth(2);
	l3->SetLineStyle(7);
	l3->Draw();
	// ratio title
	lat.SetTextAlign(33);
	lat.SetTextColor(1);
	lat.SetTextSize(0.07);
	lat.SetNDC(1);
	lat.SetTextAngle(90);
	lat.DrawLatex(0.035, 1.0, "lower / upper");
	
	// x axis title
	lat.SetTextAngle(0);
	lat.DrawLatex(0.9, 0.15, xtitle);
	gPad->SetLogy(1);
	gPad->RedrawAxis();
 	p_ratio ->Draw();
 	c1->Update();

	TString save=name+"_ratio";
	Util::PrintNoEPS(c1, save, fOutputDir, fOutputFile);	

	delete h1;
	delete h2;
	delete h_ratio;
	delete p_plot;
	delete p_ratio;
	delete c1;

}
//___________________________________________________________________________
void MassPlotter::printHisto(TH1* h_orig, TString canvname,  Option_t *drawopt, bool logflag){
	TH1D* h = (TH1D*)h_orig->Clone("h1_copy");
	TCanvas *col = new TCanvas(canvname, "", 0, 0, 900, 700);
	col->SetFillStyle(0);
	col->SetFrameFillStyle(0);
	col->cd();
	gPad->SetFillStyle(0);
	gStyle->SetPalette(1);
	if(logflag) gPad->SetLogy(1);
	h->SetMinimum(0.01);
	gPad->SetRightMargin(0.15);
	gPad->SetLogz();
	h->SetStats(false);
	h->DrawCopy(drawopt);
	gPad->RedrawAxis();
	Util::PrintNoEPS(col, canvname, fOutputDir, fOutputFile);
	h->Write();
	delete col;
	delete h;

}

//____________________________________________________________________________
void MassPlotter::printHisto(THStack h, TString canvname, Option_t *drawopt, bool logflag){
	
	TCanvas *col = new TCanvas(canvname, "", 0, 0, 900, 700);
	col->SetFillStyle(0);
	col->SetFrameFillStyle(0);
	col->cd();
	gPad->SetFillStyle(0);
	if(logflag) gPad->SetLogy(1);
	h.Draw(drawopt);
	gPad->RedrawAxis();
	col ->Update();
	gPad->RedrawAxis();
	Util::PrintNoEPS(col, canvname, fOutputDir, fOutputFile);
	delete col;

}
//____________________________________________________________________________
void MassPlotter::printHisto(THStack h, TH1* h_data, TH1* h_mc_sum, TLegend* leg,  TString canvname, Option_t *drawopt, bool logflag, TString xtitle, TString ytitle){

	TCanvas *col = new TCanvas(canvname, "", 0, 0, 900, 700);
	col->SetFillStyle(0);
	col->SetFrameFillStyle(0);
	col->cd();
	gPad->SetFillStyle(0);
	if(logflag) {
		gPad->SetLogy(1);
		h.SetMinimum(0.01);
		h_mc_sum -> SetMinimum(0.01);
		h_data->SetMinimum(0.01);
	}else{
		h.SetMinimum(0);
	}

	h.Draw(drawopt);
	h_mc_sum -> Draw("same, E2");
	if(h_data->Integral()>0) {
		h_data       ->Draw("same");
	}
	if(leg != NULL ){
		leg -> SetFillColor(0);
		leg -> SetBorderSize(0);
		leg -> Draw();
	} 
	gPad->RedrawAxis();
	col ->Update();
	// title
	TLatex lat;
	lat.SetNDC(1);
	lat.SetTextColor(4);
	lat.SetTextSize(0.04);

	// y axis title 
	lat.SetTextAlign(33); 
	lat.SetTextColor(1);
	lat.SetTextAngle(90);
	lat.DrawLatex(0.03, 0.9, ytitle);
	
	// x axis title
	lat.SetTextAngle(0);
	lat.DrawLatex(0.9, 0.05, xtitle);
	gPad->RedrawAxis();
	Util::PrintNoEPS(col, canvname, fOutputDir, fOutputFile);
	delete col;

}
//____________________________________________________________________________
void MassPlotter::printHisto(THStack h, TH1* h_data, TLegend* leg,  TString canvname, Option_t *drawopt, bool logflag , TString xtitle, TString ytitle){

	TCanvas *col = new TCanvas(canvname, "", 0, 0, 900, 700);
	col->SetFillStyle(0);
	col->SetFrameFillStyle(0);
	col->cd();
	gPad->SetFillStyle(0);
	if(logflag) {
		gPad->SetLogy(1);
		h.SetMinimum(0.01);
		h_data->SetMinimum(0.01);
	}else{
		h.SetMinimum(0);
	}

	h.Draw(drawopt);
	if(h_data->Integral()>0) {
		h_data  ->Draw("same, E1");
	}
	if(leg != NULL ){
		leg -> SetFillColor(0);
		leg -> SetBorderSize(0);
		leg -> Draw();
	} 
	gPad->RedrawAxis();
	col ->Update();

	// title
	TLatex lat;
	lat.SetNDC(1);
	lat.SetTextColor(4);
	lat.SetTextSize(0.04);

	// y axis title 
	lat.SetTextAlign(33); 
	lat.SetTextColor(1);
	lat.SetTextAngle(90);
	lat.DrawLatex(0.03, 0.9, ytitle);
	
	// x axis title
	lat.SetTextAngle(0);
	lat.DrawLatex(0.9, 0.05, xtitle);
	gPad->RedrawAxis();
	col ->Update();
	Util::PrintNoEPS(col, canvname, fOutputDir, fOutputFile);
	delete col;

}
//____________________________________________________________________________
void MassPlotter::loadSamples(const char* filename){
	char buffer[200];
	ifstream IN(filename);

	char ParName[100], StringValue[1000];
	float ParValue;

	if(fVerbose > 0) cout << "------------------------------------" << endl;
	if(fVerbose > 0) cout << "Sample File  " << filename << endl;
	int counter(0);
	
	while( IN.getline(buffer, 200, '\n') ){
		// ok = false;
		if (buffer[0] == '#') {
			continue; // Skip lines commented with '#'
		}
		if( !strcmp(buffer, "GENERAL") ){
			IN.getline(buffer, 200, '\n');
			sscanf(buffer, "Path\t%s", StringValue);
			fPath = StringValue;	
			cout << fPath << endl;
			
			if(fVerbose >0){
				cout << " ----  " << endl;
				cout << "  Path " << fPath << endl;
			}

		}
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
			TString file =fPath+StringValue;
			TFile *f = TFile::Open(file);
			s.file = f;
			s.tree = (TTree*)f->Get("MassAnalysis");
			
			IN.getline(buffer, 200, '\n');
			sscanf(buffer, "Xsection\t%f", &ParValue);
			s.xsection = ParValue;
			
			IN.getline(buffer, 200, '\n');
			sscanf(buffer, "Nevents\t%f", &ParValue);
			s.nevents = ParValue;
			
			IN.getline(buffer, 200, '\n');
			sscanf(buffer, "Kfact\t%f", &ParValue);
			s.kfact = ParValue;
			
			IN.getline(buffer, 200, '\n');
			sscanf(buffer, "Lumi\t%f", &ParValue);
			s.lumi = ParValue;

			IN.getline(buffer, 200, '\n');
			sscanf(buffer, "Type\t%s", &StringValue);
			s.type = StringValue;
			
			IN.getline(buffer, 200, '\n');
			sscanf(buffer, "Color\t%f", &ParValue);
			s.color = ParValue;

			if(fVerbose > 0){
				cout << " ---- " << endl;
				cout << "  New sample added: " << s.name << endl;
				cout << "   Sample no.      " << counter << endl;
				cout << "   Short name:     " << s.sname << endl;
				cout << "   File:           " << (s.file)->GetName() << endl;
				cout << "   Events:         " << s.nevents  << endl;
				cout << "   Events in tree: " << s.tree->GetEntries() << endl; 
				cout << "   Xsection:       " << s.xsection << endl;
				cout << "   Lumi:           " << s.lumi << endl;
				cout << "   kfactor:        " << s.kfact << endl;
				cout << "   type:           " << s.type << endl;
				cout << "   Color:          " << s.color << endl;
			}
			fSamples.push_back(s);
			counter++;
		}
	}
	if(fVerbose > 0) cout << "------------------------------------" << endl;
}

