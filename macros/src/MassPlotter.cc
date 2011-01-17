/*****************************************************************************
*   Small Class to make plots for MassAnalysis                               *
*****************************************************************************/

#include "MassPlotter.hh"

#include "TLorentzVector.h"
#include "TGraphAsymmErrors.h"
#include "helper/Utilities.hh"
#include "helper/Monitor.hh"

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
	//Util::SetStyle();
	loadSamples(filename);

}

//___________________________________________________________________________
void MassPlotter::makePlots(){
	
	double dPhisplit[4]={0., 1.0, 1.5, 3.142};
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

	std::ostringstream cutStream;
	cutStream  
		  << "misc.LeptConfig==9"                           << "&&"
		  << "NJetsIDLoose > 1"                             << "&&"
		  << "jet[0].IsGoodPFJet(100,2.4,1)==1"             << "&&"
		  << "jet[1].IsGoodPFJet(100,2.4,1)==1"             << "&&"
		  << "misc.Vectorsumpt<70"                          << "&&"
	//	  << "GetNjets(50.,2.4,0)-GetNjets(50.,2.4,1)==0"   << "&&"
	//	  << "GetMinR12R21(0,20,5.0,1)>0.5"                 << "&&" 
		  << "MinMetJetDPhi(0,20,5.0,1)>0.3"                << "&&" 
	//	  << "GetMT2Hemi(0,false,1,20,3,1)>90"              << "&&"
	//	  << "GetMT2Hemi(0,false,1,20,3,1)<180"             << "&&"
	//	  << "misc.MT2 > 100"                      << "&&"
	//	  << "pfmet[0].Pt()>100"                            << "&&" 
		  << "misc.EcalDeadCellBEFlag==1"                   << "&&"
	//	  << "GetMT2Hemi(0,false,1,20,3,1)/pfmet[0].Pt()>0.75" << "&&"
	//	  << "pfmet[0].Pt()>100"                            << "&&"
		  << "GetMT2HemiNoISR(false,4,2,2,1)>100"           << "&&"
		  << "misc.HBHENoiseFlag == 1"                   ;
	
	TString cuts = cutStream.str().c_str();
	
	MakePlot(fSamples,"GetMT2HemiNoISR(false,4,2,2,1)" ,cuts, -3, 0, "MT2" ,gNMT2bins,gMT2bins,   
		 false,  true , true,   true,  true,  false);

}

//________________________________________________________________________

void MassPlotter::makePlot(TString var, TString cuts, int njets, int nleps, TString xtitle,
			   const int nbins, const double min, const double max,
			   bool cleaned, bool logflag, bool composited, bool ratio, bool stacked, 
			   bool overlaySUSY, float overlayScale){

  MakePlot(fSamples, var, cuts, njets, nleps, xtitle, nbins, min, max, cleaned, logflag, composited, ratio, stacked, overlaySUSY, overlayScale);
  

}

//________________________________________________________________________

void MassPlotter::MakePlot(TString var, TString cuts, int njets, int nleps, TString xtitle, 
			   const int nbins,  const double *bins,
			   bool cleaned, bool logflag, bool composited, bool ratio, bool stacked, 
			   bool overlaySUSY, float overlayScale){

  MakePlot(fSamples, var, cuts, njets, nleps, xtitle, nbins, bins, cleaned, logflag, composited, ratio, stacked, overlaySUSY, overlayScale);


}

//________________________________________________________________________

void MassPlotter::PrintCutFlow(int njets, int nleps){
  
  Monitor counters[fSamples.size()];
  
  TString  cnames[12] = {"QCD", "W+jets", "Z+jets", "Top","Other", "Total Bkg.", "data", "LM0", "LM1", "LM4", "LM9", "LM13"};
  Monitor  ccount[12], ccount_100[12];

  for(size_t i = 0; i < fSamples.size(); ++i){

    Double_t weight = fSamples[i].xsection * fSamples[i].kfact * fSamples[i].lumi / (fSamples[i].nevents);
    if(fVerbose>2) cout << "PrintCutFlow: looping over " << fSamples[i].name << endl;
    if(fVerbose>2) cout << "              sample has weight " << weight << " and " << fSamples[i].tree->GetEntries() << " entries" << endl; 

    fMT2tree = new MT2tree();
    fSamples[i].tree->SetBranchAddress("MT2tree", &fMT2tree);
    Long64_t nentries =  fSamples[i].tree->GetEntries();
    Long64_t nbytes = 0, nb = 0;
    int nev =0;
    for (Long64_t jentry=0; jentry<nentries;jentry++) {
      nb =  fSamples[i].tree->GetEntry(jentry);   nbytes += nb;
      fSamples[i].tree->SetBranchAddress("MT2tree", &fMT2tree);
      
      if ( fVerbose>2 && jentry % 50000 == 0 )  cout << "+++ Proccessing event " << jentry << endl;


       bool isMT2gt100 = fMT2tree->misc.MT2 > 100.;

      //if( fMT2tree->NJetsIDLoose < 1 || !fMT2tree->jet[0].IsGoodPFJet(100, 2.4, 1) )  continue;
      //if( fMT2tree->NJetsIDLoose < 2 || !fMT2tree->jet[1].IsGoodPFJet(100, 2.4, 1) )  continue;
      if( fMT2tree->NJetsIDLoose < 1 || !fMT2tree->misc.Jet0Pass )  continue;
      if( fMT2tree->NJetsIDLoose < 2 || !fMT2tree->misc.Jet1Pass )  continue;

      if     (njets>1) {
	if(fMT2tree->NJetsIDLoose != njets) continue;
	TString text = TString::Format("All events (jets == %d)",njets);
	counters[i].fill(text.Data(),weight);
	FillMonitor(ccount, fSamples[i].sname, fSamples[i].type, text, weight);
	if(isMT2gt100)     FillMonitor(ccount_100, fSamples[i].sname, fSamples[i].type, text, weight);
      }
      else if(njets <-1) {
	if(fMT2tree->NJetsIDLoose < abs(njets)) continue;
	TString text = TString::Format("All events(jets >= %d)",abs(njets));
	counters[i].fill(text.Data(),weight);
	FillMonitor(ccount, fSamples[i].sname, fSamples[i].type, text, weight);
	if(isMT2gt100)     FillMonitor(ccount_100, fSamples[i].sname, fSamples[i].type, text, weight);
      }

      //if( fMT2tree-> MinMetJetDPhi(0,20) < 0.3 )  continue;
      if( fMT2tree->misc.MinMetJetDPhi < 0.3 )  continue;
      counters[i].fill("Minimum DPhi(MET,jet) > 0.3",weight);
      FillMonitor(ccount, fSamples[i].sname, fSamples[i].type, "Minimum DPhi(MET,jet) > 0.3", weight);
      if(isMT2gt100)     FillMonitor(ccount_100, fSamples[i].sname, fSamples[i].type, "Minimum DPhi(MET,jet) > 0.3", weight);

      if( fMT2tree->misc.HBHENoiseFlag != 1 )  continue;
      counters[i].fill("HBHE noise veto",weight);
      FillMonitor(ccount, fSamples[i].sname, fSamples[i].type, "HBHE noise veto", weight);
      if(isMT2gt100)     FillMonitor(ccount_100, fSamples[i].sname, fSamples[i].type, "HBHE noise veto", weight);

      if( fMT2tree->misc.EcalDeadCellBEFlag != 1 )  continue;
      counters[i].fill("Boundary energy veto (dead ecal)",weight);
      FillMonitor(ccount, fSamples[i].sname, fSamples[i].type, "Boundary energy veto (dead ecal)", weight);
      if(isMT2gt100)     FillMonitor(ccount_100, fSamples[i].sname, fSamples[i].type, "Boundary energy veto (dead ecal)", weight);

      if( fMT2tree->misc.Vectorsumpt > 70. )  continue;
      counters[i].fill("VectorSumPt < 70",weight);
      FillMonitor(ccount, fSamples[i].sname, fSamples[i].type, "VectorSumPt < 70", weight);
      if(isMT2gt100)     FillMonitor(ccount_100, fSamples[i].sname, fSamples[i].type, "VectorSumPt < 70", weight);

      //if( !fMT2tree->PassJetID(50, 6. ,1) )  continue;
      if( !fMT2tree->misc.PassJetID )  continue;
      counters[i].fill("jets > 50GeV failing PFID event veto",weight);
      FillMonitor(ccount, fSamples[i].sname, fSamples[i].type, "jets > 50GeV failing PFID event veto", weight);
      if(isMT2gt100)     FillMonitor(ccount_100, fSamples[i].sname, fSamples[i].type, "jets > 50GeV failing PFID event veto", weight);

      // Only considering (so far): no lepton requirement (<0);  1 lepton (==1), lepton veto (otherwise)
      if(nleps==1){ // exactly 1 lepton
	if( fMT2tree->misc.LeptConfig != 0 && fMT2tree->misc.LeptConfig != 1)  continue;
	counters[i].fill("NLeps == 1",weight);
	FillMonitor(ccount, fSamples[i].sname, fSamples[i].type, "NLeps == 1", weight);
	if(isMT2gt100)     FillMonitor(ccount_100, fSamples[i].sname, fSamples[i].type, "NLeps == 1", weight);

      }
      else if( !(nleps<0) ) { // lepton veto (no lepton requirement if nleps < 0)
	if( fMT2tree->misc.LeptConfig != 9 )  continue;
	counters[i].fill("Lepton veto",weight);
	FillMonitor(ccount, fSamples[i].sname, fSamples[i].type, "Lepton veto", weight);
	if(isMT2gt100)     FillMonitor(ccount_100, fSamples[i].sname, fSamples[i].type, "Lepton veto", weight);
      }
      
      if( fMT2tree->misc.MT2 < 80. )  continue;
      counters[i].fill("MT2 > 80 GeV" ,weight);
      FillMonitor(ccount, fSamples[i].sname, fSamples[i].type, "MT2 > 80 GeV", weight);
      if( fMT2tree->misc.MT2 < 100. )  continue;
      counters[i].fill("MT2 > 100 GeV",weight);
      FillMonitor(ccount, fSamples[i].sname, fSamples[i].type, "MT2 > 100 GeV", weight);
      if( fMT2tree->misc.MT2 < 120. )  continue;
      counters[i].fill("MT2 > 120 GeV",weight);
      FillMonitor(ccount, fSamples[i].sname, fSamples[i].type, "MT2 > 120 GeV", weight);
      if( fMT2tree->misc.MT2 < 135. )  continue;
      counters[i].fill("MT2 > 135 GeV",weight);
      FillMonitor(ccount, fSamples[i].sname, fSamples[i].type, "MT2 > 135 GeV", weight);
      if( fMT2tree->misc.MT2 < 150. )  continue;
      counters[i].fill("MT2 > 150 GeV",weight);
      FillMonitor(ccount, fSamples[i].sname, fSamples[i].type, "MT2 > 150 GeV", weight);
      if( fMT2tree->misc.MT2 < 165. )  continue;
      counters[i].fill("MT2 > 165 GeV",weight);
      FillMonitor(ccount, fSamples[i].sname, fSamples[i].type, "MT2 > 165 GeV", weight);
      if( fMT2tree->misc.MT2 < 180. )  continue;
      counters[i].fill("MT2 > 180 GeV",weight);
      FillMonitor(ccount, fSamples[i].sname, fSamples[i].type, "MT2 > 180 GeV", weight);
      if( fMT2tree->misc.MT2 < 200. )  continue;
      counters[i].fill("MT2 > 200 GeV",weight);
      FillMonitor(ccount, fSamples[i].sname, fSamples[i].type, "MT2 > 200 GeV", weight);
      if( fMT2tree->misc.MT2 < 225. )  continue;
      counters[i].fill("MT2 > 225 GeV",weight);
      FillMonitor(ccount, fSamples[i].sname, fSamples[i].type, "MT2 > 225 GeV", weight);
      if( fMT2tree->misc.MT2 < 250. )  continue;
      counters[i].fill("MT2 > 250 GeV",weight);
      FillMonitor(ccount, fSamples[i].sname, fSamples[i].type, "MT2 > 250 GeV", weight);
      if( fMT2tree->misc.MT2 < 275. )  continue;
      counters[i].fill("MT2 > 275 GeV",weight);
      FillMonitor(ccount, fSamples[i].sname, fSamples[i].type, "MT2 > 275 GeV", weight);
      if( fMT2tree->misc.MT2 < 300. )  continue;
      counters[i].fill("MT2 > 300 GeV",weight);
      FillMonitor(ccount, fSamples[i].sname, fSamples[i].type, "MT2 > 300 GeV", weight);
      
    }
    delete fMT2tree;
  }

  std::cout << setfill('=') << std::setw(70) << "" << std::endl;
  std::cout << "Statistics - by sample" << std::endl;
  std::cout << setfill('_') << std::setw(70) << "" << setfill(' ') << std::endl;
  for ( size_t i = 0; i < fSamples.size(); ++i){
    std::cout << "++++  " << fSamples[i].name << std::endl;  
    std::cout << setfill('-') << std::setw(70) << "" << setfill(' ') << std::endl;    
    counters[i].print();
    std::cout << setfill('_') << std::setw(70) << "" << setfill(' ') << std::endl;    
  }  

  std::cout << setfill('=') << std::setw(70) << "" << std::endl;
  std::cout << "Statistics - by process" << std::endl;
  std::cout << setfill('_') << std::setw(70) << "" << setfill(' ') << std::endl;
  for ( size_t i = 0; i < 12; ++i){
    std::cout << "++++  " << cnames[i] << std::endl;  
    std::cout << setfill('-') << std::setw(70) << "" << setfill(' ') << std::endl;    
    ccount[i].print();
    std::cout << setfill('_') << std::setw(70) << "" << setfill(' ') << std::endl;    
  }  

  std::cout << setfill('=') << std::setw(70) << "" << std::endl;
  std::cout << "Statistics - by process (MT2 > 100GeV)" << std::endl;
  std::cout << setfill('_') << std::setw(70) << "" << setfill(' ') << std::endl;
  for ( size_t i = 0; i < 12; ++i){
    std::cout << "++++  " << cnames[i] << std::endl;  
    std::cout << setfill('-') << std::setw(70) << "" << setfill(' ') << std::endl;    
    ccount_100[i].print();
    std::cout << setfill('_') << std::setw(70) << "" << setfill(' ') << std::endl;    
  }  

}

//________________________________________________________________________

void MassPlotter::FillMonitor(Monitor *count, TString sname, TString type, TString cut, float weight){
  if     (sname == "QCD"    ) 	count[ 0].fill(cut.Data(),weight);
  else if(sname == "Wtolnu" ) 	count[ 1].fill(cut.Data(),weight);
  else if(sname == "DY"     ) 	count[ 2].fill(cut.Data(),weight);
  else if(sname == "Top"    )	count[ 3].fill(cut.Data(),weight);
  else if(sname == "VV" || sname == "VGamma" || sname == "Photons") 	
                                count[ 4].fill(cut.Data(),weight);
  if     (type  == "mc"     )   count[ 5].fill(cut.Data(),weight);
  else if(type  == "data"   )	count[ 6].fill(cut.Data(),weight);
  else if(sname == "LM0"    )	count[ 7].fill(cut.Data(),weight);
  else if(sname == "LM1"    )	count[ 8].fill(cut.Data(),weight);
  else if(sname == "LM4"    )	count[ 9].fill(cut.Data(),weight);
  else if(sname == "LM9"    )	count[10].fill(cut.Data(),weight);
  else if(sname == "LM13"    )	count[11].fill(cut.Data(),weight);
}

//________________________________________________________________________

void MassPlotter::plotSig(TString var, TString cuts, TString xtitle, int nbins, double min, double max, bool cleaned, int type ){

        TString varname = Util::removeFunnyChar(var.Data());

	TH1D*    h_mc_sum    = new TH1D   (varname+"mc_sum", "", nbins, min, max );
	TH1D*    h_susy      = new TH1D   (varname+"susy"  , "", nbins, min, max );	
	h_mc_sum    ->Sumw2();
	h_susy      ->Sumw2();
	// vector of all histos
	vector<TH1D*> h_samples;

	for(size_t i = 0; i < fSamples.size(); ++i){
	        if(!(fSamples[i].type=="susy") && !(fSamples[i].type=="mc")) continue;

		h_samples.push_back(new TH1D(varname+"_"+fSamples[i].name, "", nbins, min, max));
		h_samples[i] -> Sumw2();

		Double_t weight = fSamples[i].xsection * fSamples[i].kfact * fSamples[i].lumi / (fSamples[i].nevents);

		if(fVerbose>2) cout << "MakePlot: looping over " << fSamples[i].sname << endl;
		if(fVerbose>2) cout << "           sample has weight " << weight << " and " << fSamples[i].tree->GetEntries() << " entries" << endl; 
	
		TString variable  = TString::Format("%s>>%s",var.Data(),h_samples[i]->GetName());
		TString selection = TString::Format("(%f) * (%s)",weight,cuts.Data());
		  
		if(fVerbose>2) cout << "+++++ Drawing " << variable  << endl
				    << "\twith cuts: "  << selection << endl;

		int nev = fSamples[i].tree->Draw(variable.Data(),selection.Data(),"goff");
		
		if(fVerbose>2) cout << "\tevents found : "  <<  nev << endl
				    << "\t->Integral() : "  <<  h_samples[i]->Integral() << endl;
		
		if( fSamples[i].type=="mc"        )   h_mc_sum->Add(h_samples[i]);
		else if( fSamples[i].type=="susy" )   h_susy  ->Add(h_samples[i]);
	}
	Float_t  x[nbins], y[nbins];
	for (int i = 1; i <=nbins+1; i++){
	  x[i-1] = h_mc_sum->GetBinLowEdge(i);
	  float s = h_susy  ->Integral(i,nbins+1);
	  float b = h_mc_sum->Integral(i,nbins+1);
	  switch (type){
	  case 0:
	    y[i-1] = s/sqrt(b);
	    break;
	  case 1:
	    y[i-1] = s/sqrt(s+b);
	    break;
	  case 2:
	    y[i-1] = s/b;
	    break;
	  }
	}
	TGraph *sig = new TGraph(nbins+1,x,y);
	sig->SetTitle("");
	sig->GetXaxis()->SetTitle(xtitle);
	sig->SetMarkerStyle(20);
	switch (type){
	case 0:
	  sig->GetYaxis()->SetTitle("S/#sqrt{B}");
	  break;
	case 1:
	  sig->GetYaxis()->SetTitle("S/#sqrt{S+B}");
	  break;
	case 2:
	  sig->GetYaxis()->SetTitle("S/B");
	  break;
	}
	sig->Draw("ACP");
}

//________________________________________________________________________

void MassPlotter::MakePlot(std::vector<sample> Samples, TString var, TString cuts, int njets, int nleps,
			   TString xtitle, const int nbins, const double min, const double max,
			   bool cleaned, bool logflag, bool composited, bool ratio, 
			   bool stacked, bool overlaySUSY, float overlayScale){

  double bins[nbins];
  bins[0] = min;
  for(int i=1; i<=nbins; i++)
    bins[i] = min+i*(max-min)/nbins;
  MakePlot(Samples, var, cuts, njets, nleps, xtitle, nbins, bins, cleaned, logflag, composited, ratio, stacked, overlaySUSY, overlayScale);

}

//________________________________________________________________________

void MassPlotter::MakePlot(std::vector<sample> Samples, TString var, TString cuts, int njets, int nleps,
			   TString xtitle, const int nbins, const double *bins, 
			   bool cleaned, bool logflag, bool composited, bool ratio, 
			   bool stacked, bool overlaySUSY, float overlayScale){

        TString varname = Util::removeFunnyChar(var.Data());

	TString nJets = "NJetsIDLoose";
	nJets += njets < 0 ? ">=" : "==";
	nJets += TString::Format("%d",abs(njets));

	TString  nLeps = nleps < 0 ? "" : (nleps==1 ? " && (misc.LeptConfig == 0 || misc.LeptConfig == 1)" : "&& misc.LeptConfig ==  9");

	THStack* h_stack     = new THStack(varname, "");
  	TH1D*    h_data      = new TH1D   (varname+"data"  , "", nbins, bins );
	TH1D*    h_mc_sum    = new TH1D   (varname+"mc_sum", "", nbins, bins );
	TH1D*    h_susy      = new TH1D   (varname+"susy"  , "", nbins, bins );	

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
	h_susy      ->Sumw2();


	// vector of all histos
	vector<TH1D*> h_samples;

	// arrays of composited stuff
	TH1D    *h_composited[7];
	//TString  cnames[5] = {"QCD", "W/Z/#gamma production", "Top production", "susy", "data"};
	//int      ccolor[5] = {401, 418, 602, 0, 632};
	TString  cnames[7] = {"QCD", "W+jets", "Z+jets", "Top", "Other", "susy", "data"};
	int      ccolor[7] = {401, 417, 419, 600, 603, 0, 632};
	for (int i=0; i<7; i++){
	  h_composited[i] = new TH1D(varname+"_"+cnames[i], "", nbins, bins);
	  h_composited[i] -> Sumw2();
	  h_composited[i] -> SetFillColor  (stacked ? ccolor[i] : 0);
	  h_composited[i] -> SetLineColor  (ccolor[i]);
	  if (!stacked) {
	    h_composited[i] -> SetLineWidth(4);
	  }
	  h_composited[i] -> SetMarkerColor(ccolor[i]);
	  h_composited[i] -> SetStats(false);
	}

	// legend
	TLegend* Legend1 = new TLegend(.71,.54,.91,.92);
	//TLegend* Legend1 = new TLegend(.3,.5,.6,.88);

	for(size_t i = 0; i < Samples.size(); ++i){
		h_samples.push_back(new TH1D(varname+"_"+Samples[i].name, "", nbins, bins));
		h_samples[i] -> Sumw2();
		h_samples[i] -> SetFillColor(stacked ? Samples[i].color : 0);
		h_samples[i] -> SetLineColor(Samples[i].color);
		if (!stacked) {
		  h_samples[i] -> SetLineWidth(4);
		}
		h_samples[i] -> SetMarkerColor(Samples[i].color);
		h_samples[i] -> SetStats(false);
		if(Samples[i].type == "susy" ){
			h_samples[i] -> SetLineColor(kBlack);
			h_samples[i] -> SetLineStyle(kDotted);
			h_composited[5] -> SetLineColor(kBlack);
			h_composited[5] -> SetLineStyle(kDotted);
		}
		Double_t weight = Samples[i].xsection * Samples[i].kfact * Samples[i].lumi / (Samples[i].nevents);
		if(fVerbose>2) cout << "MakePlot: looping over " << Samples[i].sname << endl;
		if(fVerbose>2) cout << "           sample has weight " << weight << " and " << Samples[i].tree->GetEntries() << " entries" << endl; 
	
		TString variable  = TString::Format("%s>>%s",var.Data(),h_samples[i]->GetName());
		TString theCuts = nJets + nLeps + "&&" + cuts;
		TString selection = TString::Format("(%f) * (%s)",weight,theCuts.Data());
		  
		if(fVerbose>2) cout << "+++++ Drawing " << variable  << endl
				    << "\twith cuts: "  << selection << endl;

		int nev = Samples[i].tree->Draw(variable.Data(),selection.Data(),"goff");

		// Add underflow & overflow bins
		// This failed for older ROOT version when the first(last) bin is empty
		// and there are underflow (overflow) events --- must check whether this 
		// is still the case
		h_samples[i]->SetBinContent(1,
				    h_samples[i]->GetBinContent(0) + h_samples[i]->GetBinContent(1));
		h_samples[i]->SetBinError(1,
				  sqrt(h_samples[i]->GetBinError(0)*h_samples[i]->GetBinError(0)+
				       h_samples[i]->GetBinError(1)*h_samples[i]->GetBinError(1) ));
		h_samples[i]->SetBinContent(h_samples[i]->GetNbinsX(),
					    h_samples[i]->GetBinContent(h_samples[i]->GetNbinsX()  )+ 
					    h_samples[i]->GetBinContent(h_samples[i]->GetNbinsX()+1) );
		h_samples[i]->SetBinError(h_samples[i]->GetNbinsX(),
					  sqrt(h_samples[i]->GetBinError(h_samples[i]->GetNbinsX()  )*
					       h_samples[i]->GetBinError(h_samples[i]->GetNbinsX()  )+
					       h_samples[i]->GetBinError(h_samples[i]->GetNbinsX()+1)*
					       h_samples[i]->GetBinError(h_samples[i]->GetNbinsX()+1)  ));
		
		if(fVerbose>2) cout << "\tevents found : "  <<  nev << endl
				    << "\t->Integral() : "  <<  h_samples[i]->Integral() << endl;
		
		if (Samples[i].sname.Contains("QCD")) {
		  h_composited[0]->Add(h_samples[i]);
		}
		else if(Samples[i].sname.Contains("Top") ){
		  h_composited[3]->Add(h_samples[i]);
		}
		else if(Samples[i].sname=="Wtolnu") {
		  h_composited[1]->Add(h_samples[i]);
		}
		else if(Samples[i].sname=="DY") {
		  h_composited[2]->Add(h_samples[i]);
		}
		else if(Samples[i].sname=="VV" || Samples[i].sname=="VGamma" || Samples[i].sname=="Photons") {
		  h_composited[4]->Add(h_samples[i]);
		}
		else if(Samples[i].type.Contains("susy")){
		  h_composited[5]->Add(h_samples[i]);
		  cnames[5] = Samples[i].sname;
		}
		else if(Samples[i].type.Contains("data")){
		  h_composited[6]->Add(h_samples[i]);
		}

	}

	if (composited){
	  for (int i=0; i<6; ++i){
	    if (!stacked)  {
	      h_composited[i]->Scale(1./h_composited[i]->Integral());
	      h_composited[i] -> SetMinimum(0.00005);
	      //if (fVerbose>2) cout << " h_composited[" << i<< "]->Integral = " << h_composited[i]->Integral()  << endl;
	    }
	    if( i!=5 || !overlaySUSY) h_stack  -> Add(h_composited[i]);
	    if( i==5)                 h_susy   -> Add(h_composited[i]);
	    else        	      h_mc_sum -> Add(h_composited[i]);
	    if( i==5 && overlaySUSY)
	      Legend1 ->AddEntry(h_composited[i], cnames[i] + (overlayScale ? 
							      TString::Format(" x %.0f",overlayScale) :
							      " scaled to data")   , "f");  
	    else
	      Legend1 ->AddEntry(h_composited[i], cnames[i], stacked ? "f" : "l");
	  }
	  h_data->Add(h_composited[6]);
	  if(h_data->Integral()>0 && stacked){
	    Legend1     ->AddEntry(h_data, "data", "l");
	  }
	  if(fVerbose > 2 && composited) {
		           cout << "------------------------------------"                << endl
	                        << "QCD Integral:      " << h_composited[0]->Integral()  << endl
			        << "W+jets Integral:   " << h_composited[1]->Integral()  << endl
			        << "Z+jets Integral:   " << h_composited[2]->Integral()  << endl
			        << "Top Integral:      " << h_composited[3]->Integral()  << endl
			        << "Other Integral:    " << h_composited[4]->Integral()  << endl
			        << "SUSY:              " << h_composited[5]->Integral()  << endl
			        << "Data:              " << h_composited[6]->Integral()  << endl;
	  }
	}
	else{
	  for(int i=0; i<Samples.size(); ++i){
	    if (!stacked) {
	      h_samples[i]->Scale(1./h_samples[i]->Integral());
	      h_samples[i] -> SetMinimum(0.0001);
	    }
	    if((Samples[i].sname=="QCD" )){
	      h_samples[i] -> SetMinimum(0.01);
	      h_stack->Add(h_samples[i]);
	      h_mc_sum->Add(h_samples[i]);
	      Legend1 ->AddEntry(h_samples[i], Samples[i].sname, "f");
	    }
	    if(Samples[i].sname!="QCD" && Samples[i].type!="data"){
	      if(Samples[i].type=="susy") h_susy->Add(h_samples[i]);
	      if(Samples[i].type!="susy" || !overlaySUSY) h_stack->Add(h_samples[i]);
	      if(Samples[i].type!="susy") h_mc_sum->Add(h_samples[i]);
	      if(Samples[i].type=="susy" && overlaySUSY){
	         Legend1 ->AddEntry(h_samples[i], Samples[i].sname + (overlayScale ? 
								       TString::Format(" x %.0f",overlayScale) :
								       " scaled to data")   , "f") ;
	      }else{
		 Legend1 ->AddEntry(h_samples[i], Samples[i].sname, stacked ? "f" : "l");
	      }	
	    }
	    if(Samples[i].type == "data"){
	      h_data -> Add(h_samples[i]);
	    }
	    
	  }
	  if(h_data->Integral()>0){
	    Legend1     ->AddEntry(h_data, "data", "l");
	  }
	}

	TString ytitle = "Events";

	TString oname = var+ "_CUT_";
	oname += njets < 0 ? TString::Format("ge%dJets",abs(njets)) : TString::Format("%dJets",abs(njets));
	oname += nleps < 0 ? "_anyLep_" : (nleps == 1 ? "_1Lep_" : "_0Lep_");
	oname += cuts;
	oname.ReplaceAll(">=" ,".ge");
	oname.ReplaceAll("<=" ,".le");
	oname.ReplaceAll(">" ,".gt");
	oname.ReplaceAll("<" ,".lt");
	oname.ReplaceAll("==",".eq");
	oname.ReplaceAll("!=",".ne");
	oname.ReplaceAll("&&","_");
	oname.ReplaceAll("||","_");
	oname.ReplaceAll("misc.","");
	oname.ReplaceAll("LeptConfig","LepCfg");
	oname.ReplaceAll("Vectorsumpt","VSPT");
	oname.ReplaceAll("EcalDeadCellBEFlag","BEFlg");
	oname.ReplaceAll("HBHENoiseFlag","HBHEFlg");
	oname.ReplaceAll("NJetsIDLoose","NJLoose");
	oname.ReplaceAll("isPFIDLoose","isJLoose");
	oname.ReplaceAll("IsGoodPFJet","IsGoodPFJ");
	oname.ReplaceAll("MinMetJetDPhi","MinDPhi");
	oname.ReplaceAll(",","-");
        TString outname = Util::removeFunnyChar(oname.Data());

	if(!stacked) {
	  outname = outname + (cleaned ? "_cleaned" : "") + (logflag ? "_log" : "") + "_shape";
	  printHisto(h_stack, h_data, h_mc_sum, Legend1 , outname, "hist", logflag, xtitle, ytitle, njets, nleps, false);

	}
	else if (!overlaySUSY){
	  outname = outname + (cleaned ? "_cleaned" : "") + (logflag ? "_log" : "") + (composited ? "_comp" : "");
	  printHisto(h_stack, h_data, h_mc_sum, Legend1 , outname, "hist", logflag, xtitle, ytitle, njets, nleps);
	  if (ratio) 
	    plotRatioStack(h_stack,  h_mc_sum, h_data,  logflag, false, outname, Legend1, xtitle, ytitle, njets, nleps);
	}
	else {
	  outname =  outname + (cleaned ? "_cleaned" : "") + (logflag ? "_log" : "") + (composited ? "_comp" : "") + "_overlay";	  
	  printHisto(h_stack, h_data, h_mc_sum, h_susy, Legend1 , outname, "hist", logflag, xtitle, ytitle, njets, nleps, overlayScale);
	  if (ratio) 
	    plotRatioStack(h_stack,  h_mc_sum, h_data, h_susy,  logflag, false, outname, Legend1, xtitle, ytitle, njets, nleps, overlayScale);

	}
	
// 	for(int i=0; i<h_samples.size(); ++i){
// 		delete h_samples[i];
// 	}
// 	delete h_mc_sum;
// 	delete Legend1;
// 	delete h_data;
// 	delete h_susy;
	
}

//________________________________________________________________________
void MassPlotter::MakePlot(std::vector<sample> Samples, TString branch_name, const int nbins, const double bins[], bool cleaned, bool logflag, 
		           TString cut_branch_name, double ysplit[], TString option, TString version, TString prediction, double factor ){
	
	double lower_cut =  ysplit[0];
	double upper_cut =  ysplit[1];
	THStack* h_stack     = new THStack(branch_name, "");
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
	TLegend* Legend1 = new TLegend(.6,.5,.89,.88);
	
	// get prediction
	TH1D* h_prediction    = new TH1D("predict",    "", nbins, bins);
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
		Double_t weight = Samples[i].xsection * Samples[i].kfact * Samples[i].lumi / (Samples[i].nevents);
		if(fVerbose>2) cout << "MakePlot: looping over " << Samples[i].sname << endl;
		if(fVerbose>2) cout << "           sample has weight " << weight << " and " << Samples[i].tree->GetEntries() << " entries" << endl; 
	
		fMT2tree = new MT2tree();
		Samples[i].tree->SetBranchAddress("MT2tree", &fMT2tree);
		Long64_t nentries =  Samples[i].tree->GetEntries();
		Long64_t nbytes = 0, nb = 0;
		int nev =0;
		for (Long64_t jentry=0; jentry<nentries;jentry++) {
			nb =  Samples[i].tree->GetEntry(jentry);   nbytes += nb;
		        Samples[i].tree->SetBranchAddress("MT2tree", &fMT2tree);

			if ( fVerbose>2 && jentry % 5000 == 0 )  cout << "+++ Proccessing event " << jentry << endl;
	
			// cleaning & cuts
			if( fMT2tree->misc.DPhiMhtMpt > upper_cut) continue;
			if( fMT2tree->misc.DPhiMhtMpt < lower_cut) continue;
			if( fMT2tree->misc.LeptConfig !=LeptConfigCut  ) continue; //lepton veto
			if( fMT2tree->misc.Vectorsumpt > VectorSumPtCut) continue;
			if( fMT2tree->NJets < NJetsCut )           continue;
			if( fMT2tree->NJets > MaxNJetsCut )        continue;
			if( fMT2tree->misc.MT2 <0)        continue;			
			// fill histo
			if(option == "none")  {
				h_samples[i]                ->Fill(fMT2tree->misc.MT2     , weight);
				if(Samples[i].type == "mc"  &&  Samples[i].sname == "QCD" && fMT2tree->misc.MT2 < 50){h_prediction -> Fill(fMT2tree->misc.MT2     , weight);}
			}
			nev++;
		}
		if(fVerbose>2) cout << "\tEvents found:  " << nev << endl
				    << "\t->Integral() : "  <<  h_samples[i]->Integral() << endl;

		delete fMT2tree;
	}
	if(prediction!="none"){
		h_stack ->Add(h_prediction);
		Legend1 ->AddEntry(h_prediction, "data-driven QCD", "f");
		h_mc_sum->Add(h_prediction);
	}
	for(int i=0; i<Samples.size(); ++i){
		if(Samples[i].type=="mc" || Samples[i].type=="susy") h_mc_sum_susy->Add(h_samples[i]);
		if((Samples[i].sname=="QCD" && prediction=="none")){
			h_samples[i] -> SetMinimum(0.01);
			h_stack ->Add(h_samples[i]);
			h_mc_sum->Add(h_samples[i]);
			Legend1 ->AddEntry(h_samples[i], Samples[i].sname, "f");
		}
		if(Samples[i].sname!="QCD" && Samples[i].type!="data"){
			h_stack->Add(h_samples[i]);
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

		Double_t weight = fSamples[i].xsection * fSamples[i].kfact * fSamples[i].lumi / (fSamples[i].nevents);

		fMT2tree = new MT2tree();
		fSamples[i].tree->SetBranchAddress("MT2tree", &fMT2tree);
		Long64_t nentries =  fSamples[i].tree->GetEntries();
		Long64_t nbytes = 0, nb = 0;
		
		for (Long64_t jentry=0; jentry<nentries;jentry++) {
			nb =  fSamples[i].tree->GetEntry(jentry);   nbytes += nb;
				
			// cleaning & cuts
			if( fMT2tree->misc.DPhiMhtMpt > upper_cut) continue;
			if( fMT2tree->misc.DPhiMhtMpt < lower_cut) continue;
			if( fMT2tree->misc.LeptConfig !=LeptConfigCut  ) continue; //lepton veto
			if( fMT2tree->misc.Vectorsumpt > VectorSumPtCut) continue;
			if( fMT2tree->NJets < NJetsCut )           continue;
			if( fMT2tree->NJets > MaxNJetsCut )        continue;
			if( fMT2tree->misc.MT2 <0)        continue;			
			// fill histo
			if(fMT2tree->misc.MT2 > 50) h_prediction     ->Fill(fMT2tree->misc.MT2     , weight);
		}
		delete fMT2tree;
	}
	h_prediction->Multiply(h_fudge);
	return h_prediction;	
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

		fMT2tree = new MT2tree();
		fSamples[i].tree->SetBranchAddress("MT2tree", &fMT2tree);
		Long64_t nentries =  fSamples[i].tree->GetEntries();
		Long64_t nbytes = 0, nb = 0;
		
		Double_t weight = fSamples[i].xsection * fSamples[i].kfact * fSamples[i].lumi / (fSamples[i].nevents);
		if(fVerbose>2) cout << "ABCD: looping over " << fSamples[i].sname << endl;
		if(fVerbose>2) cout << "      sample has weight " << weight << endl; 
		
		for (Long64_t jentry=0; jentry<nentries;jentry++) {
			nb =  fSamples[i].tree->GetEntry(jentry);   nbytes += nb;
			
			// cleaning & cuts
			if( fMT2tree->misc.LeptConfig !=LeptConfigCut  ) continue; //lepton veto
			if( fMT2tree->misc.Vectorsumpt > VectorSumPtCut) continue;
			if( fMT2tree->NJets < NJetsCut )           continue;
			if( fMT2tree->NJets > MaxNJetsCut )        continue;
			if( fMT2tree->misc.MT2 <0)        continue;			
			
			// fill histo
			if(  fMT2tree->misc.DPhiMhtMpt > ysplit[0] && fMT2tree->misc.DPhiMhtMpt < ysplit[1] ){
				h_ABCD_lower_y_band -> Fill(fMT2tree->misc.MT2, weight);
				h_ABCD_MT2          -> Fill(fMT2tree->misc.MT2, fMT2tree->misc.DPhiMhtMpt, weight);
			}
			if(  fMT2tree->misc.DPhiMhtMpt > ysplit[2] && fMT2tree->misc.DPhiMhtMpt < ysplit[3] ){
				h_ABCD_upper_y_band -> Fill(fMT2tree->misc.MT2, weight);
				h_ABCD_MT2          -> Fill(fMT2tree->misc.MT2, fMT2tree->misc.DPhiMhtMpt, weight);
			}

		}
		delete fMT2tree;
	}

		
	h_ABCD_upper_y_band->SetLineColor(kBlack);
	h_ABCD_upper_y_band->SetLineStyle(kBlack);	

	h_ABCD_lower_y_band->SetMarkerStyle(1);
	h_ABCD_lower_y_band->SetLineColor(kRed);
	h_ABCD_lower_y_band->SetLineStyle(kDashed);
	h_ABCD_lower_y_band->SetMarkerColor(kRed);
	
	// title and legend for ratio plot
 	TLegend* Legend1 = new TLegend(.55,.6,.89,.88);
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
void MassPlotter::plotRatioStack(THStack* hstack, TH1* h1_orig, TH1* h2_orig, bool logflag, bool normalize, TString name, TLegend* leg, TString xtitle, TString ytitle,int njets, int nleps){
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
	if(logflag) max = 2*max;
	else max = 1.05*max;
	h1->SetMaximum(max);
	h2->SetMaximum(max);
	hstack->SetMaximum(max);
//	h1->SetMinimum(0.000000001);
//	h2->SetMinimum(0.000000001);

	hstack->SetMinimum(0.02);
	hstack ->Draw("hist");
	h2     ->Draw("same");

	// title
	TLatex lat;
	lat.SetNDC(1);
	lat.SetTextColor(4);
	lat.SetTextSize(0.1);

	// y axis title 
	lat.SetTextAlign(33); 
	lat.SetTextColor(1);
	lat.SetTextAngle(90);
	lat.DrawLatex(0.035, 0.9, ytitle);
		
	TLatex TitleBox;
	TitleBox.SetNDC();
	TitleBox.SetTextSize(0.05);
	TString text = njets < 0 ? TString::Format("#geq %d Jets",abs(njets)) : TString::Format("%d Jets",abs(njets));
	text += nleps == 1 ? ", 1 Lepton" : "";
	TitleBox.DrawLatex(0.18,0.943,text.Data());

 	p_plot ->Draw();
	gPad->RedrawAxis();

	if(leg != NULL ){
		leg -> SetFillColor(0);
		leg -> SetBorderSize(0);
		leg -> Draw();
	} 
	
	// draw the ratio plot
 	p_ratio ->cd();
	//gPad->SetLogy();
 
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
 	h_ratio ->SetMinimum(0.4);
 	h_ratio ->SetMaximum(3.0);
	h_ratio ->GetYaxis()->SetTitleOffset(h1->GetYaxis()->GetTitleOffset());

	h_ratio ->DrawCopy("E2");
 
	TLine *l3 = new TLine(h1->GetXaxis()->GetXmin(), 1.00, h1->GetXaxis()->GetXmax(), 1.00);
	l3->SetLineWidth(2);
	l3->SetLineStyle(7);
	l3->Draw();
	// ratio title
	lat.SetTextAlign(33);
	lat.SetTextColor(1);
	lat.SetTextSize(0.15);
	lat.SetNDC(1);
	lat.SetTextAngle(90);
	lat.DrawLatex(0.035, 1.0, "data / MC");
	
	// x axis title
	lat.SetTextSize(0.2);
	lat.SetTextAngle(0);
	float ypos = xtitle.Contains("slash") || xtitle.Contains("_") || xtitle.Contains("^") ? 0.23 : 0.16;
	lat.DrawLatex(0.9, ypos, xtitle);
	gPad->RedrawAxis();
 	p_ratio ->Draw();
 	c1->Update();

	TString save=name+"_ratio";
	Util::Print(c1, save, fOutputDir, fOutputFile);	

// 	delete h1;
// 	delete h2;
// 	delete h_ratio;
// 	delete p_plot;
// 	delete p_ratio;
// 	delete c1;

}
//_________________________________________________________________________________
void MassPlotter::plotRatioStack(THStack* hstack, TH1* h1_orig, TH1* h2_orig, TH1* h3, bool logflag, bool normalize, TString name, TLegend* leg, TString xtitle, TString ytitle,int njets, int nleps, float overlayScale){
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
	if(logflag) max = 2*max;
	else max = 1.05*max;
	h1->SetMaximum(max);
	h2->SetMaximum(max);
	hstack->SetMaximum(max);
//	h1->SetMinimum(0.000000001);
//	h2->SetMinimum(0.000000001);

	hstack->SetMinimum(0.02);
	hstack->Draw("hist");
	h2    ->Draw("same");
	h3->Scale(overlayScale ? overlayScale : h2->Integral() / h3->Integral());
	h3->SetFillColor(0);
	h3->SetLineStyle(kDotted);
	h3->Draw("samehist");

	// title
	TLatex lat;
	lat.SetNDC(1);
	lat.SetTextColor(4);
	lat.SetTextSize(0.1);

	// y axis title 
	lat.SetTextAlign(33); 
	lat.SetTextColor(1);
	lat.SetTextAngle(90);
	lat.DrawLatex(0.035, 0.9, ytitle);
		
	TLatex TitleBox;
	TitleBox.SetNDC();
	TitleBox.SetTextSize(0.05);
	TString text = njets < 0 ? TString::Format("#geq %d Jets",abs(njets)) : TString::Format("%d Jets",abs(njets));
	text += nleps == 1 ? ", 1 Lepton" : "";
	TitleBox.DrawLatex(0.18,0.943,text.Data());

 	p_plot ->Draw();
	gPad->RedrawAxis();

	if(leg != NULL ){
		leg -> SetFillColor(0);
		leg -> SetBorderSize(0);
		leg -> Draw();
	} 
	
	// draw the ratio plot
 	p_ratio ->cd();
	//gPad->SetLogy();
 
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
 	h_ratio ->SetMinimum(0.4);
 	h_ratio ->SetMaximum(3.0);
	h_ratio ->GetYaxis()->SetTitleOffset(h1->GetYaxis()->GetTitleOffset());

	h_ratio ->DrawCopy("E2");
 
	TLine *l3 = new TLine(h1->GetXaxis()->GetXmin(), 1.00, h1->GetXaxis()->GetXmax(), 1.00);
	l3->SetLineWidth(2);
	l3->SetLineStyle(7);
	l3->Draw();
	// ratio title
	lat.SetTextAlign(33);
	lat.SetTextColor(1);
	lat.SetTextSize(0.15);
	lat.SetNDC(1);
	lat.SetTextAngle(90);
	lat.DrawLatex(0.035, 1.0, "data / MC");
	
	// x axis title
	lat.SetTextSize(0.2);
	lat.SetTextAngle(0);
	float ypos = xtitle.Contains("slash") || xtitle.Contains("_") || xtitle.Contains("^") ? 0.23 : 0.16;
	lat.DrawLatex(0.9, ypos, xtitle);
	gPad->RedrawAxis();
 	p_ratio ->Draw();
 	c1->Update();

	TString save=name+"_ratio";
	Util::Print(c1, save, fOutputDir, fOutputFile);	

// 	delete h1;
// 	delete h2;
// 	delete h_ratio;
// 	delete p_plot;
// 	delete p_ratio;
// 	delete c1;

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
	if(logflag) max = 2*max;
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
	lat.SetTextSize(0.07);

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
 	h_ratio ->SetMinimum(0.4);
 	h_ratio ->SetMaximum(3.0);
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
	float ypos = xtitle.Contains("slash") || xtitle.Contains("_") || xtitle.Contains("^") ? 0.23 : 0.16;
	lat.DrawLatex(0.9, ypos, xtitle);
	//gPad->SetLogy(1);
	gPad->RedrawAxis();
 	p_ratio ->Draw();
 	c1->Update();

	TString save=name+"_ratio";
	Util::PrintNoEPS(c1, save, fOutputDir, fOutputFile);	

// 	delete h1;
// 	delete h2;
// 	delete h_ratio;
// 	delete p_plot;
// 	delete p_ratio;
// 	delete c1;

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
void MassPlotter::printHisto(THStack* h, TString canvname, Option_t *drawopt, bool logflag){
	
	TCanvas *col = new TCanvas(canvname, "", 0, 0, 900, 700);
	col->SetFillStyle(0);
	col->SetFrameFillStyle(0);
	col->cd();
	gPad->SetFillStyle(0);
	if(logflag) gPad->SetLogy(1);
	h->Draw(drawopt);
	gPad->RedrawAxis();
	col ->Update();
	gPad->RedrawAxis();
	Util::PrintNoEPS(col, canvname, fOutputDir, fOutputFile);
	delete col;

}
//____________________________________________________________________________
void MassPlotter::printHisto(THStack* h, TH1* h_data, TH1* h_mc_sum, TLegend* leg,  TString canvname, Option_t *drawopt, bool logflag, TString xtitle, TString ytitle,int njets, int nleps, bool stacked){

	TCanvas *col = new TCanvas(canvname, "", 0, 0, 900, 700);
	col->SetFillStyle(0);
	col->SetFrameFillStyle(0);
	col->cd();
	gPad->SetFillStyle(0);
	if(logflag) {
		gPad->SetLogy(1);
		h        -> SetMinimum(stacked ? 0.05 : 0.0001);
		h_mc_sum -> SetMinimum(0.05);
		h_data   -> SetMinimum(0.05);
	}else{
		h->SetMinimum(0);
	}

	// Determine plotting range
	double max1 = h_data->GetMaximum();
	double max2 = h_mc_sum->GetMaximum();
	double max  = (max1>max2)?max1:max2;
	if(logflag) max = 2*max;
	else max = 1.05*max;
	h_data  ->SetMaximum(max);
	h_mc_sum->SetMaximum(max);
	h       ->SetMaximum(max);

	h->Draw(stacked ? drawopt : "histnostack");
	//h_mc_sum -> Draw("same, E2");
	if(h_data->Integral()>0 && stacked) {
		h_data       ->Draw("same");
	}
	if(leg != NULL ){
		leg -> SetFillColor(0);
		leg -> SetBorderSize(0);
		leg -> Draw();
	} 
	gPad->RedrawAxis();
	col ->Update();

	TLatex TitleBox;
	TitleBox.SetNDC();
	TitleBox.SetTextSize(0.05);
	TString text = njets < 0 ? TString::Format("#geq %d Jets",abs(njets)) : TString::Format("%d Jets",abs(njets));
	text += nleps == 1 ? ", 1 Lepton" : "";
	TitleBox.DrawLatex(0.18,0.943,text.Data());

	// title
	TLatex lat;
	lat.SetNDC(1);
	lat.SetTextColor(4);
	lat.SetTextSize(0.07);

	// y axis title 
	lat.SetTextAlign(33); 
	lat.SetTextColor(1);
	lat.SetTextAngle(90);
	lat.DrawLatex(0.03, 0.9, ytitle);
	
	// x axis title
	lat.SetTextAngle(0);
	float ypos = xtitle.Contains("slash") || xtitle.Contains("_") || xtitle.Contains("^") ? 0.1 : 0.06;
	lat.DrawLatex(0.9, ypos, xtitle);
	gPad->RedrawAxis();
	Util::PrintNoEPS(col, canvname, fOutputDir, fOutputFile);
	Util::PrintEPS(col, canvname, fOutputDir);
// 	delete col;

}
//____________________________________________________________________________
void MassPlotter::printHisto(THStack* h, TH1* h_data, TH1* h_mc_sum, TH1* h_susy, TLegend* leg,  TString canvname, Option_t *drawopt, bool logflag, TString xtitle, TString ytitle,int njets, int nleps, float overlayScale){

	TCanvas *col = new TCanvas(canvname, "", 0, 0, 900, 700);
	col->SetFillStyle(0);
	col->SetFrameFillStyle(0);
	col->cd();
	gPad->SetFillStyle(0);
	if(logflag) {
		gPad->SetLogy(1);
		h        -> SetMinimum(0.05);
		h_mc_sum -> SetMinimum(0.05);
		h_data   -> SetMinimum(0.05);
		h_susy   -> SetMinimum(0.05);
	}else{
		h->SetMinimum(0);
	}

	// Determine plotting range
	double max1 = h_data->GetMaximum();
	double max2 = h_mc_sum->GetMaximum();
	double max  = (max1>max2)?max1:max2;
	if(logflag) max = 2*max;
	else max = 1.05*max;
	h_data  ->SetMaximum(max);
	h_mc_sum->SetMaximum(max);
	h       ->SetMaximum(max);

	h->Draw(drawopt);
	//h_mc_sum -> Draw("same, E2");
	if(h_data->Integral()>0) {
		h_data       ->Draw("same");
	}
	h_susy->Scale(overlayScale ? overlayScale : h_data->Integral() / h_susy->Integral());
	h_susy->SetLineStyle(kDotted);
	h_susy->SetFillColor(0);
	h_susy->Draw("samehist");
	if(leg != NULL ){
		leg -> SetFillColor(0);
		leg -> SetBorderSize(0);
		leg -> Draw();
	} 
	gPad->RedrawAxis();
	col ->Update();

	TLatex TitleBox;
	TitleBox.SetNDC();
	TitleBox.SetTextSize(0.05);
	TString text = njets < 0 ? TString::Format("#geq %d Jets",abs(njets)) : TString::Format("%d Jets",abs(njets));
	text += nleps == 1 ? ", 1 Lepton" : "";
	TitleBox.DrawLatex(0.18,0.943,text.Data());

	// title
	TLatex lat;
	lat.SetNDC(1);
	lat.SetTextColor(4);
	lat.SetTextSize(0.07);

	// y axis title 
	lat.SetTextAlign(33); 
	lat.SetTextColor(1);
	lat.SetTextAngle(90);
	lat.DrawLatex(0.03, 0.9, ytitle);
	
	// x axis title
	lat.SetTextAngle(0);
	float ypos = xtitle.Contains("slash") || xtitle.Contains("_") || xtitle.Contains("^") ? 0.1 : 0.06;
	lat.DrawLatex(0.9, ypos, xtitle);
	gPad->RedrawAxis();
	Util::PrintNoEPS(col, canvname, fOutputDir, fOutputFile);
	Util::PrintEPS(col, canvname, fOutputDir);
// 	delete col;

}
//____________________________________________________________________________
void MassPlotter::printHisto(THStack* h, TH1* h_data, TLegend* leg,  TString canvname, Option_t *drawopt, bool logflag , TString xtitle, TString ytitle){

	TCanvas *col = new TCanvas(canvname, "", 0, 0, 900, 700);
	col->SetFillStyle(0);
	col->SetFrameFillStyle(0);
	col->cd();
	gPad->SetFillStyle(0);
	if(logflag) {
		gPad->SetLogy(1);
		h->SetMinimum(0.01);
		h_data->SetMinimum(0.01);
	}else{
		h->SetMinimum(0);
	}

	// Determine plotting range
	double max1 = h_data->GetMaximum();
	double max2 = h     ->GetMaximum();
	double max  = (max1>max2)?max1:max2;
	if(logflag) max = 2*max;
	else max = 1.05*max;
	h_data  ->SetMaximum(max);
	h       ->SetMaximum(max);

	h->Draw(drawopt);
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
	lat.SetTextSize(0.07);

	// y axis title 
	lat.SetTextAlign(33); 
	lat.SetTextColor(1);
	lat.SetTextAngle(90);
	lat.DrawLatex(0.03, 0.9, ytitle);
	
	// x axis title
	lat.SetTextAngle(0);
	float ypos = xtitle.Contains("slash") || xtitle.Contains("_") || xtitle.Contains("^") ? 0.1 : 0.06;
	lat.DrawLatex(0.9, ypos, xtitle);
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
			s.tree = (TTree*)f->Get("MassTree");
			
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
			sscanf(buffer, "Type\t%s", StringValue);
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

