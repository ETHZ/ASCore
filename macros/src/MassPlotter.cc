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
void MassPlotter::makeSmallCopy(int nevents, int sample){
	if(fSamples.size()<sample) return;
	cout << "making small copy of tree " << fSamples[sample].name << endl;
	TFile *newfile = new TFile(fOutputDir+"/"+fSamples[sample].name+"_small.root","recreate");
	TTree *newtree = fSamples[sample].tree->CloneTree(0);
	Int_t nentries = (Int_t)fSamples[sample].tree->GetEntries();

	for (Int_t i=0;i<nentries; i++) {
		if(i > nevents) continue;
		fSamples[sample].tree->GetEntry(i);
		newtree->Fill();
	}


	newfile->Write();

	delete newfile;

}
//___________________________________________________________________________
void MassPlotter::makePlots(){
	double dPhisplit[4]={0., 1.0, 1.5, 3.142};
	MakeMT2PredictionAndPlots(false, dPhisplit, 2.5);  // not cleaned, fudgefactor 2
}

// _________________________________________________________________________
void MassPlotter::makeZnunu(){
	// get samples
	vector<sample> Samples;
	cout << fSamples.size() << endl;
	for(int i=0; i<fSamples.size(); ++i){
	//	if(fSamples[i].sname=="DYToLL"           )   Samples.push_back(fSamples[i]); 
	//	if(fSamples[i].sname=="DYToNuNu"         ) Samples.push_back(fSamples[i]); 
		if(fSamples[i].sname.Contains("QCD")     ) Samples.push_back(fSamples[i]); 
	}
	
	// -------------
	// reco and geom efficiency
	for(int i=0; i<fSamples.size(); ++i){
//		if(fSamples[i].type =="data"    )     PrintZllEfficiency(fSamples[i], true, "ele");
//		if(fSamples[i].sname=="DYToLL")       PrintZllEfficiency(fSamples[i], false, "ele");	
//		if(fSamples[i].sname=="DYToNuNu")     PrintZllEfficiency(fSamples[i], false, "neutrinos");	
	}

	// --------------
	// comparing Z->nunu with Z->ll
	std::ostringstream cutStream;
	cutStream  
		  << "misc.HT       >300"                           << "&&"
	//	  << "misc.MET      >30"                            << "&&"
		  << "misc.Jet0Pass == 1"                           << "&&"
		  << "misc.Jet1Pass == 1"                           << "&&"
		  << "misc.PassJetID == 1"                          << "&&"
	//	  << "misc.Vectorsumpt<70"                          << "&&"
	//	  << "misc.MinMetJetDPhi>0.3"                       << "&&" 
	//	  << "misc.EcalDeadCellBEFlag==1"                   << "&&"
		  << "misc.HBHENoiseFlag == 1"                   ;
	
	TString cuts = cutStream.str().c_str();

	// declare a few maps	
	// to be used for Z->nunu                               to be used for Z->ll
	RemoveLeptMap["GetMT2Hemi(0,false,1,20,3,1)"]                = "GetMT2Hemi(0,false,1,20,3,3)";
//	RemoveLeptMap["misc.MET"]                                    = "GetMETPlusLepts(1)";
//	RemoveLeptMap["pfmet[0].Pt()"]                               = "GetDiLeptonPt()";
	RemoveLeptMap["pfmet[0].Pt()"]                               = "GetMETPlusGenLepts(0,1,1,1113,23,0,100,0,10000)";
	RemoveLeptMap["misc.MET"]                                    = "GetMHT(0,15,5)";
	RemoveLeptMap["GetMHT(0,20,5)"]                              = "GetMHT(0,20,3)";
	RemoveLeptMap["MinMetJetDPhi(0,20,2.4,1)"]                   = "MinMetJetDPhi(0,20,2.4,2)";
	RemoveLeptMap["GetPseudoJetsMETmindPhi(2,3,0,20,2.4,1)"]     = "GetPseudoJetsMETmindPhi(2,3,0,20,2.4,2)";

	TString isZtoll    = "GetDiLeptonInvMass(0,1,0,10,true) > 60 && GetDiLeptonInvMass(0,1,0,10,true) < 120";
	TString removed_ee = "GetMETPlusGenLepts(0,1,1,1113,23,0,100,0,10000)>=0";
	TString isZtoll2   = "GetDiLeptonPt()>-1";
//                            variable for Z->nunu                     cuts,  optcut, replace_cut   xtitle        bins                 add_underflow   logflag  scale_factor normalized
//	CompSamples(Samples, "NJetsIDLoose",                           cuts, isZtoll,       false,      "NJetsIDLoose", 10 , 0 , 10,         false,          true,    1,           true);
//	CompSamples(Samples, "misc.HT"     ,                           cuts, "_",           false,      "HT"          , 50, 50, 500,         false,          true,    1,           true);
//	CompSamples(Samples, "GetGenMET(0,1,14,23,10,2.4,60,120)",     cuts, "_",           true,       "GetMET"      , gNMT2bins, gMT2bins, false,          true,    1,           true);
//	CompSamples(Samples, "GetMHT(1,20,2.4)"                  ,     cuts, isZtoll,       false,      "MHT"         , gNMT2bins, gMT2bins, false,          true,    1,           true);
//	CompSamples(Samples, "pfmet[0].Pt()"                     ,     cuts, removed_ee,    true,       "MET"         , gNMT2bins, gMT2bins, false,          true,    1,           true);
//	CompSamples(Samples, "pfmet[0].Pt()"                     ,     cuts, removed_ee,    true,       "MET"         , 100, 0, 500        , false,          true,    1,           true);
//	CompSamples(Samples, "GetMT2Hemi(0,false,1,20,3,1)"      ,     cuts, "_",           true,       "MT2"         , gNMT2bins, gMT2bins ,false,          true,    1,           true);
	CompSamples(Samples, "GetPseudoJetsMETmindPhi(2,3,0,20,2.4,1)",cuts, "_",           true,       "minJetMETdPhi" ,50,        0, 3.2  ,false,          true,    1,           true);
	
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
		  << "misc.MET > 30"                                << "&&"
		  << "misc.HT > 300"                                << "&&"
		  << "misc.Jet0Pass == 1"                           << "&&"
		  << "misc.Jet1Pass == 1"                           << "&&"
		  << "misc.PassJetID == 1"                          << "&&"
		  << "misc.Vectorsumpt<70"                          << "&&"
		  << "misc.MinMetJetDPhi>0.3"                       << "&&" 
		  << "misc.EcalDeadCellBEFlag==1"                   << "&&"
		  << "misc.HBHENoiseFlag == 1"                   ;

	TString cuts = cutStream.str().c_str();
	
	//                 variable                            cuts njets  nlepts title     bins               cleaned  log  composite    ratio  stacked overlay 
	MakePlot(fSamples,"misc.MT2" ,                         cuts, -2,   0,      "MT2" , gNMT2bins, gMT2bins , false,  true ,  true,      true,  true,  false);

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

// ________________________________________________________________________
void MassPlotter::PrintZllEfficiency(sample Sample , bool data, std::string lept){
	Long64_t nevents=1000000000;

      	std::cout << setfill('=') << std::setw(70) << "" << std::endl;
	cout << "printing kinetic & geom. acceptance for Z->ll for sample: \n"
	     << Sample.name << endl;	
      	std::cout << setfill('-') << std::setw(70) << "" << std::endl;

        enum counters_t { count_begin, all=count_begin, presel,  HCAL_ECAL_noise, VectorSumPt ,PassJetID, MinMetJetDPhi, MET, count_end };
 	Monitor counters[count_end];
	TString lablx[count_end] = {"all events", "presel", "Cal_Noise", "VectorSumPt", "PassJetID", "MinMetJetDPhi", "MET"};

	int pid, flavour; 
	if(lept == "ele" )            {pid = 11; flavour = 1;} 
	else if(lept == "muo" )       {pid = 13; flavour = 2;} 
	else if(lept != "neutrinos")  {cout << "choose appropriate lepton flavour" << endl; return;}

    	fMT2tree = new MT2tree();
    	Sample.tree->SetBranchAddress("MT2tree", &fMT2tree);
    	Long64_t nentries =  Sample.tree->GetEntries();
    	Long64_t nbytes = 0, nb = 0;
    	for (Long64_t jentry=0; jentry<min(nentries, nevents);jentry++) {
      		nb = Sample.tree->GetEntry(jentry);   nbytes += nb;
      		Sample.tree->SetBranchAddress("MT2tree", &fMT2tree);
      
      		if ( fVerbose>2 && jentry % 50000 == 0 )  cout << "+++ Proccessing event " << jentry << endl;

		// check if event has Z->e+e- within acceptance	
		
		bool Zee(false);         
		bool ZeeAccepted(false); 
		bool ZeeReco(false); 
		if(! data) Zee              = fMT2tree->IsGenOSDiLepton(pid,23,0,10,0,10000);
		if(  data) Zee              = (fMT2tree->GetDiLeptonInvMass(0,1,flavour,10,1) > 60 && fMT2tree->GetDiLeptonInvMass(0,1,flavour,10,1) < 120);
		if(! data) ZeeAccepted      = fMT2tree->IsGenOSDiLepton(pid,23,10,2.4,60,120);
		if(! data) ZeeReco          = (fMT2tree->GetDiLeptonInvMass(0,1,flavour,10,1) > 60 && fMT2tree->GetDiLeptonInvMass(0,1,flavour,10,1) < 120);
		if(lept=="neutrinos"){ Zee  = true;}

		// all events
		if(Zee)         counters[all].fill("all events");
		if(ZeeAccepted) counters[all].fill("Zee within acceptance");
		if(ZeeReco)     counters[all].fill("Zee Reco");
		
		// presel
		if(fMT2tree->misc.HT  < 300)             continue;
		if(fMT2tree->misc.Jet0Pass  ==0       )  continue;
		if(fMT2tree->misc.Jet1Pass  ==0       )  continue;
		if(Zee)         counters[presel].fill("presel");
		if(ZeeAccepted) counters[presel].fill("Zee within acceptance");
		if(ZeeReco)     counters[presel].fill("Zee Reco");
		
		// ECAL HCAL Noise
		if(fMT2tree->misc.EcalDeadCellBEFlag  == 0)  continue;
		if(fMT2tree->misc.HBHENoiseFlag       == 0)  continue;
		if(Zee)         counters[HCAL_ECAL_noise].fill("Cal_Noise");
		if(ZeeAccepted) counters[HCAL_ECAL_noise].fill("Zee within acceptance");
		if(ZeeReco)     counters[HCAL_ECAL_noise].fill("Zee Reco");
		
		// VectorSumPt
		if(fMT2tree->misc.Vectorsumpt  > 70       )  continue;
		if(Zee)         counters[VectorSumPt].fill("VectorSumPt");
		if(ZeeAccepted) counters[VectorSumPt].fill("Zee within acceptance");
		if(ZeeReco)     counters[VectorSumPt].fill("Zee Reco");
		
		// PassJetID
		if(fMT2tree->misc.PassJetID  ==0      )  continue;
		if(Zee)         counters[PassJetID].fill("PassJetID");
		if(ZeeAccepted) counters[PassJetID].fill("Zee within acceptance");
		if(ZeeReco)     counters[PassJetID].fill("Zee Reco");
		
		// MinMetJetDPhi
		if(fMT2tree->misc.MinMetJetDPhi  <0.3      )  continue;
		if(Zee)         counters[MinMetJetDPhi].fill("MinMetJetDPhi");
		if(ZeeAccepted) counters[MinMetJetDPhi].fill("Zee within acceptance");
		if(ZeeReco)     counters[MinMetJetDPhi].fill("Zee Reco");
		
		// MET
		if(fMT2tree->misc.MET<30                   )  continue;
		if(Zee)         counters[MET].fill("MET");
		if(ZeeAccepted) counters[MET].fill("Zee within acceptance");
		if(ZeeReco)     counters[MET].fill("Zee Reco");
      	}

	// print stats     
      	std::cout << setfill('=') << std::setw(70) << "" << std::endl;
       	std::cout << "Statistics" << std::endl;
       	std::cout << setfill('-') << std::setw(70) << "" << setfill(' ') << std::endl;
       	for ( counters_t iCount=count_begin; iCount<count_end; iCount = counters_t(iCount+1) ) {
         	counters[iCount].print();
         	std::cout << setfill('-') << std::setw(70) << "" << setfill(' ') << std::endl;  
       	}
  	std::cout << setfill('=') << std::setw(70) << "" << std::endl;

	if(data) return;

	// fill histo
	TH1D* h_Zaccept_num = new TH1D("h_Zaccept_num", "", count_end, 0., (double) count_end );
	TH1D* h_Zaccept_den = new TH1D("h_Zaccept_den", "", count_end, 0., (double) count_end );
	TH1D* h_Zaccept     = new TH1D("h_Zaccept"    , "", count_end, 0., (double) count_end );
	TH1D* h_ZReco_num   = new TH1D("h_ZReco_num"  , "", count_end, 0., (double) count_end );
	TH1D* h_ZReco_den   = new TH1D("h_ZReco_den"  , "", count_end, 0., (double) count_end );
	TH1D* h_ZReco       = new TH1D("h_ZReco"      , "", count_end, 0., (double) count_end );
	for(int i=0; i<count_end; ++i){
		h_Zaccept    ->GetXaxis()->SetBinLabel(i+1, lablx[i]);
		h_Zaccept_num->SetBinContent(i+1,counters[i].counts("Zee within acceptance"));
		h_Zaccept_den->SetBinContent(i+1,counters[i].counts((string) lablx[i]));
		h_ZReco      ->GetXaxis()->SetBinLabel(i+1, lablx[i]);
		h_ZReco_num  ->SetBinContent(i+1,counters[i].counts("Zee Reco"));
		h_ZReco_den  ->SetBinContent(i+1,counters[i].counts("Zee within acceptance"));
	}
	h_Zaccept    ->Sumw2();
	h_Zaccept    ->Divide(h_Zaccept_num,h_Zaccept_den);	
	h_ZReco      ->Sumw2();
	h_ZReco      ->Divide(h_ZReco_num,h_ZReco_den);	
	TCanvas *col  = new TCanvas("GenDYToLLAcceptance", "", 0, 0, 900, 700);
	TCanvas *col2 = new TCanvas("GenDYToLLReco"      , "", 0, 0, 900, 700);
	
	col -> cd();
	h_Zaccept    ->SetLineColor(kBlue);
	h_Zaccept    ->Draw("E");
	string name  ="GenDYToLLAcceptance"+(string) Sample.name; 
	Util::PrintEPS(col, name, fOutputDir);
	
	col2 -> cd();
	h_ZReco    ->SetLineColor(kBlue);
	h_ZReco    ->Draw("E");
	string name2  ="GenDYToLLReco"+(string) Sample.name; 
	Util::PrintEPS(col2, name2, fOutputDir);

	delete h_Zaccept;
	delete h_Zaccept_den;
	delete h_Zaccept_num;
	delete h_ZReco;
	delete h_ZReco_den;
	delete h_ZReco_num;
	delete col;
	delete col2;
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
//       if( fMT2tree->misc.MinMetJetDPhi < 0.3 )  continue;
//       counters[i].fill("Minimum DPhi(MET,jet) > 0.3",weight);
//       FillMonitor(ccount, fSamples[i].sname, fSamples[i].type, "Minimum DPhi(MET,jet) > 0.3", weight);
//       if(isMT2gt100)     FillMonitor(ccount_100, fSamples[i].sname, fSamples[i].type, "Minimum DPhi(MET,jet) > 0.3", weight);

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
//__________________________________________________________________________
void MassPlotter::CompSamples(TString var, TString cuts, TString optcut, bool RemoveLepts, 
		              TString xtitle, const int nbins, const double *bins, bool add_underflow, bool logflag, double scale_factor, bool normalize){

  	CompSamples(fSamples, var, cuts, optcut, RemoveLepts, xtitle, nbins, bins, add_underflow, logflag, scale_factor, normalize);
}
//___________________________________________________________________________________________________________________
void MassPlotter::compSamples(TString var, TString cuts, TString optcut, bool RemoveLepts,
			   TString xtitle, const int nbins, const double min, const double max, bool add_underflow, bool logflag, double scale_factor, bool normalize){

  	CompSamples(fSamples, var, cuts, optcut, RemoveLepts, xtitle, nbins, min, max, add_underflow, logflag, scale_factor, normalize);

}
// _______________________________________________________________________
void MassPlotter::CompSamples(std::vector<sample> Samples, TString var, TString cuts, TString optcut, bool RemoveLepts,
			   TString xtitle, const int nbins, const double min, const double max, bool add_underflow, bool logflag, double scale_factor, bool normalize){

  	double bins[nbins];
  	bins[0] = min;
  	for(int i=1; i<=nbins; i++)
    	bins[i] = min+i*(max-min)/nbins;
  	CompSamples(Samples, var, cuts, optcut, RemoveLepts, xtitle, nbins, bins, add_underflow, logflag, scale_factor, normalize);

}
// __________________________________________________________________________
void MassPlotter::CompSamples(std::vector<sample> Samples, TString var, TString cuts, TString optcut, bool RemoveLepts,
			   TString xtitle, const int nbins, const double *bins, bool add_underflow, bool logflag, double scale_factor, bool normalize){

        TString varname = Util::removeFunnyChar(var.Data());
	TH1D*    histos[2];
	cout << Samples.size() << endl;
	histos[0] = new TH1D(varname+Samples[0].sname  , "", nbins, bins);
	histos[1] = new TH1D(varname+Samples[1].sname  , "", nbins, bins);

	histos[0] -> Sumw2(); 
	histos[1] -> Sumw2();
	histos[0] -> SetLineColor(kBlue); 
	histos[1] -> SetLineColor(kRed); 
	histos[0] -> SetLineWidth(4);
	histos[1] -> SetLineWidth(4);

	// legend
	TLegend* Legend1 = new TLegend(.71,.54,.91,.92);

	for(int i=0; i<2; ++i){
		Double_t weight = scale_factor * Samples[i].xsection * Samples[i].kfact * Samples[i].lumi / (Samples[i].nevents);
		if(fVerbose>2) cout << "GetHisto: looping over " << Samples[i].sname << endl;
		if(fVerbose>2) cout << "           sample has weight " << weight << " and " << Samples[i].tree->GetEntries() << " entries" << endl; 
		// exchange variable to plot with the corresponding one with removed leptons;	
		TString variable; 
		TString selection; 
		if(RemoveLepts && (Samples[i].sname == "DYToLL" || Samples[i].sname=="QCD_2")){
			MapType::iterator iter = RemoveLeptMap.begin();
			iter = RemoveLeptMap.find(var);
			if (iter != RemoveLeptMap.end() ){
				variable = iter -> second;
			} else {cout << "found RemoveLepts==true, but no corresponding map" << endl; variable = var;}
		}else {
			variable = var;
		}
		if(optcut!="_" && Samples[i].sname == "DYToLL"){
			TString newcuts = cuts + "&&" + optcut;  
			selection      = TString::Format("(%f) * (%s)",weight,newcuts.Data());
		} else {
			selection      = TString::Format("(%f) * (%s)",weight,cuts.Data());
		}

		variable          = TString::Format("%s>>%s",variable.Data(),histos[i]->GetName());
		
			  
		if(fVerbose>2) cout << "+++++ Drawing " << variable  << endl
				    << "\twith cuts: "  << selection << endl;

		int nev = Samples[i].tree->Draw(variable.Data(),selection.Data(),"goff");

		// Add underflow & overflow bins
		// This failed for older ROOT version when the first(last) bin is empty
		// and there are underflow (overflow) events --- must check whether this 
		// is still the case
		if(add_underflow) {
			histos[i]->SetBinContent(1,
				    histos[i]->GetBinContent(0) + histos[i]->GetBinContent(1));
			histos[i]->SetBinError(1,
				  sqrt(histos[i]->GetBinError(0)*histos[i]->GetBinError(0)+
				       histos[i]->GetBinError(1)*histos[i]->GetBinError(1) ));
		}
		histos[i]->SetBinContent(histos[i]->GetNbinsX(),
					    histos[i]->GetBinContent(histos[i]->GetNbinsX()  )+ 
					    histos[i]->GetBinContent(histos[i]->GetNbinsX()+1) );
		histos[i]->SetBinError(histos[i]->GetNbinsX(),
					  sqrt(histos[i]->GetBinError(histos[i]->GetNbinsX()  )*
					       histos[i]->GetBinError(histos[i]->GetNbinsX()  )+
					       histos[i]->GetBinError(histos[i]->GetNbinsX()+1)*
					       histos[i]->GetBinError(histos[i]->GetNbinsX()+1)  ));
		
		if(fVerbose>2) cout << "\tevents found : "  <<  nev << endl
				    << "\t->Integral() : "  <<  histos[i]->Integral() << endl;
	     	Legend1 ->AddEntry(histos[i], Samples[i].sname, "l"); 
	}
	
	plotRatio(histos[0], histos[1], logflag, normalize, varname+"_"+Samples[0].sname+"_"+Samples[1].sname, Legend1 , xtitle, "");
	
	
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


//__________________________________________________________________________
void MassPlotter::abcd_MT2(TString var, TString basecut, TString upper_cut, TString lower_cut, const int nbins,const double min, const double max, double fit_min, double fit_max){

  double bins[nbins];
  bins[0] = min;
  for(int i=1; i<=nbins; i++)
    bins[i] = min+i*(max-min)/nbins;
  ABCD_MT2(var, basecut, upper_cut, lower_cut, nbins, bins);
  
}

//__________________________________________________________________________
void MassPlotter::ABCD_MT2(TString var, TString basecut, TString upper_cut, TString lower_cut, const int nbins, const double *bins, double fit_min, double fit_max){
  TH2D *h_ABCD_MT2_qcd           = new TH2D("ABCD_MT2_"+var+"_qcd"           , "",  100, 0, 600, 180, 0, TMath::Pi());  h_ABCD_MT2_qcd          ->Sumw2();
  TH2D *h_ABCD_MT2_susy          = new TH2D("ABCD_MT2_"+var+"_susy"          , "",  100, 0, 600, 180, 0, TMath::Pi());  h_ABCD_MT2_susy         ->Sumw2();
  TH1D *h_ABCD_upper_y_band_susy = new TH1D("ABCD_upper_y_band_"+var+"_susy" , "",  nbins, bins);		        h_ABCD_upper_y_band_susy->Sumw2();
  TH1D *h_ABCD_lower_y_band_data = new TH1D("ABCD_lower_y_band_"+var+"_data" , "",  nbins, bins);		        h_ABCD_lower_y_band_data->Sumw2();
  TH1D *h_ABCD_upper_y_band_data = new TH1D("ABCD_upper_y_band_"+var+"_data" , "",  nbins, bins);		        h_ABCD_upper_y_band_data->Sumw2();
  TH1D *h_ABCD_lower_y_band_qcd  = new TH1D("ABCD_lower_y_band_"+var+"_qcd"  , "",  nbins, bins);		        h_ABCD_lower_y_band_qcd ->Sumw2();
  TH1D *h_ABCD_upper_y_band_qcd  = new TH1D("ABCD_upper_y_band_"+var+"_qcd"  , "",  nbins, bins);		        h_ABCD_upper_y_band_qcd ->Sumw2();
  TH1D *h_ABCD_lower_y_band_mc   = new TH1D("ABCD_lower_y_band_"+var+"_mc"   , "",  nbins, bins);		        h_ABCD_lower_y_band_mc  ->Sumw2();
  TH1D *h_ABCD_upper_y_band_mc   = new TH1D("ABCD_upper_y_band_"+var+"_mc"   , "",  nbins, bins);		        h_ABCD_upper_y_band_mc  ->Sumw2();
  TH1D *h_ABCD_lower_y_band_sub  = new TH1D("ABCD_lower_y_band_"+var+"_sub"  , "",  nbins, bins);		        h_ABCD_lower_y_band_sub ->Sumw2();
  TH1D *h_ABCD_upper_y_band_sub  = new TH1D("ABCD_upper_y_band_"+var+"_sub"  , "",  nbins, bins);		        h_ABCD_upper_y_band_sub ->Sumw2();
  TH1D *ratio_data               = new TH1D("ratio_"+var+"_data"             , "",  nbins, bins);		        ratio_data              ->Sumw2();
  TH1D *ratio_qcd                = new TH1D("ratio_"+var+"_qcd"              , "",  nbins, bins);		        ratio_qcd               ->Sumw2();
  TH1D *ratio_sub                = new TH1D("ratio_"+var+"_sub"              , "",  nbins, bins);		        ratio_sub               ->Sumw2();
  
  for(size_t i = 0; i < fSamples.size(); ++i){
    
    Long64_t nentries =  fSamples[i].tree->GetEntries();
    Long64_t nbytes = 0, nb = 0;
    
    Double_t weight = fSamples[i].xsection * fSamples[i].kfact * fSamples[i].lumi / (fSamples[i].nevents);
    if(fVerbose>2) cout << "ABCD: looping over " << fSamples[i].name << endl;
    if(fVerbose>2) cout << "      sample has weight " << weight << endl; 
    
//     TEntryList *elist = new TEntryList("elist","elist");
//     int nev = fSamples[i].tree->Draw(">>elist",basecut.Data(),"entrylist");
//     fSamples[i].tree->SetEntryList(elist);
    
    TString selection = TString::Format("(%f) * (%s)"      ,weight,basecut.Data());
    TString sel_up    = TString::Format("(%f) * (%s && %s)",weight,basecut.Data(),upper_cut.Data());
    TString sel_lo    = TString::Format("(%f) * (%s && %s)",weight,basecut.Data(),lower_cut.Data());

    TString variable;
    TString var2 = "misc.MT2";
    if (fSamples[i].type == "data"){
      variable  = TString::Format("%s>>+%s",var2.Data(),h_ABCD_lower_y_band_data->GetName());
      int nevL = fSamples[i].tree->Draw(variable,sel_lo,"goff");
      variable  = TString::Format("%s>>+%s",var2.Data(),h_ABCD_upper_y_band_data->GetName());
      int nevU = fSamples[i].tree->Draw(variable,sel_up,"goff");
      if(fVerbose>2) cout << "\tData lower, events found : "  <<  nevL << endl
			  << "\t->Integral() : "  <<  h_ABCD_lower_y_band_data->Integral() + h_ABCD_lower_y_band_data->GetBinContent(nbins+1) << endl;
      if(fVerbose>2) cout << "\tData upper, events found : "  <<  nevU << endl
			  << "\t->Integral() : "  <<  h_ABCD_upper_y_band_data->Integral() + h_ABCD_upper_y_band_data->GetBinContent(nbins+1) << endl;
    }
    else if (fSamples[i].type == "mc" && fSamples[i].name == "QCD_Pt_120to170"){   // WARNING: checking statistical fluctuation in this pt bin!!!
      variable  = TString::Format("%s:misc.MT2>>+%s",var.Data(),h_ABCD_MT2_qcd->GetName());
      TString sel_u    = TString::Format("(%f) * (%s && misc.MT2<100&& %s)",weight,basecut.Data(),upper_cut.Data());
      TString sel_l    = TString::Format("(%f) * (%s && misc.MT2<100&& %s)",weight,basecut.Data(),lower_cut.Data());
      int nev2d= fSamples[i].tree->Draw(variable,selection,"goff");
      variable  = TString::Format("%s>>+%s",var2.Data(),h_ABCD_lower_y_band_qcd->GetName());
      int nevL = fSamples[i].tree->Draw(variable,sel_l,"goff");
      variable  = TString::Format("%s>>+%s",var2.Data(),h_ABCD_upper_y_band_qcd->GetName());
      int nevU = fSamples[i].tree->Draw(variable,sel_u,"goff");
      if(fVerbose>2) cout << "\t" << fSamples[i].name << " lower, events found : "  <<  nevL << endl
			  << "\t->Integral() : "  <<  h_ABCD_lower_y_band_qcd->Integral() + h_ABCD_lower_y_band_qcd->GetBinContent(nbins+1) << endl;
      if(fVerbose>2) cout << "\t" << fSamples[i].name << " upper, events found : "  <<  nevU << endl
			  << "\t->Integral() : "  <<  h_ABCD_upper_y_band_qcd->Integral() + h_ABCD_upper_y_band_qcd->GetBinContent(nbins+1) << endl;
    }
    else if (fSamples[i].type == "mc" && fSamples[i].sname == "QCD"){
      variable  = TString::Format("%s:misc.MT2>>+%s",var.Data(),h_ABCD_MT2_qcd->GetName());
      int nev2d= fSamples[i].tree->Draw(variable,selection,"goff");
      variable  = TString::Format("%s>>+%s",var2.Data(),h_ABCD_lower_y_band_qcd->GetName());
      int nevL = fSamples[i].tree->Draw(variable,sel_lo,"goff");
      variable  = TString::Format("%s>>+%s",var2.Data(),h_ABCD_upper_y_band_qcd->GetName());
      int nevU = fSamples[i].tree->Draw(variable,sel_up,"goff");
      if(fVerbose>2) cout << "\t" << fSamples[i].name << " lower, events found : "  <<  nevL << endl
			  << "\t->Integral() : "  <<  h_ABCD_lower_y_band_qcd->Integral() + h_ABCD_lower_y_band_qcd->GetBinContent(nbins+1) << endl;
      if(fVerbose>2) cout << "\t" << fSamples[i].name << " upper, events found : "  <<  nevU << endl
			  << "\t->Integral() : "  <<  h_ABCD_upper_y_band_qcd->Integral() + h_ABCD_upper_y_band_qcd->GetBinContent(nbins+1) << endl;
    }
    else if (fSamples[i].type == "mc"){
      variable  = TString::Format("%s>>+%s",var2.Data(),h_ABCD_lower_y_band_mc->GetName());
      int nevL = fSamples[i].tree->Draw(variable,sel_lo,"goff");
      variable  = TString::Format("%s>>+%s",var2.Data(),h_ABCD_upper_y_band_mc->GetName());
      int nevU = fSamples[i].tree->Draw(variable,sel_up,"goff");
      if(fVerbose>2) cout << "\t" << fSamples[i].name << " lower, events found : "  <<  nevL << endl
			  << "\t->Integral() : "  <<  h_ABCD_lower_y_band_mc->Integral() + h_ABCD_lower_y_band_mc->GetBinContent(nbins+1) << endl;
      if(fVerbose>2) cout << "\t" << fSamples[i].name << " upper, events found : "  <<  nevU << endl
			  << "\t->Integral() : "  <<  h_ABCD_upper_y_band_mc->Integral() + h_ABCD_upper_y_band_mc->GetBinContent(nbins+1) << endl;
    }
    else if (fSamples[i].type == "susy"){
      variable  = TString::Format("%s:misc.MT2>>+%s",var.Data(),h_ABCD_MT2_susy->GetName());
      int nev2d= fSamples[i].tree->Draw(variable,selection,"goff");
      variable  = TString::Format("%s>>+%s",var2.Data(),h_ABCD_upper_y_band_susy->GetName());
      int nevU = fSamples[i].tree->Draw(variable,sel_up,"goff");
      if(fVerbose>2) cout << "\t" << fSamples[i].name << " upper, events found : "  <<  nevU << endl
			  << "\t->Integral() : "  <<  h_ABCD_upper_y_band_susy->Integral() + h_ABCD_upper_y_band_susy->GetBinContent(nbins+1) << endl;
    }
  }
  
  
  h_ABCD_lower_y_band_sub->Add(h_ABCD_lower_y_band_data,h_ABCD_lower_y_band_mc,1,-1);
  h_ABCD_upper_y_band_sub->Add(h_ABCD_upper_y_band_data,h_ABCD_upper_y_band_mc,1,-1);

  for (int i= 0; i<=nbins; i++){
    if (h_ABCD_lower_y_band_sub->GetBinContent(i)<0)  h_ABCD_lower_y_band_sub->SetBinContent(i, 0.);
    if (h_ABCD_upper_y_band_sub->GetBinContent(i)<0)  h_ABCD_upper_y_band_sub->SetBinContent(i, 0.);
  }

  ratio_qcd ->Divide(h_ABCD_upper_y_band_qcd ,h_ABCD_lower_y_band_qcd );
  ratio_data->Divide(h_ABCD_upper_y_band_data,h_ABCD_lower_y_band_data);
  ratio_sub ->Divide(h_ABCD_upper_y_band_sub ,h_ABCD_lower_y_band_sub );
  
  h_ABCD_lower_y_band_data->SetLineColor(1);  h_ABCD_upper_y_band_data->SetLineColor(2);
  h_ABCD_lower_y_band_qcd ->SetLineColor(1);  h_ABCD_upper_y_band_qcd ->SetLineColor(2);
  h_ABCD_lower_y_band_sub ->SetLineColor(1);  h_ABCD_upper_y_band_sub ->SetLineColor(2);
  h_ABCD_lower_y_band_sub ->SetLineStyle(3);  h_ABCD_upper_y_band_sub ->SetLineStyle(3);   ratio_sub->SetLineStyle(3);
  ratio_qcd ->SetMarkerStyle(24);   ratio_qcd ->SetMarkerSize(0.8);
  ratio_data->SetMarkerStyle(24);   ratio_data->SetMarkerSize(0.8);
  ratio_sub ->SetMarkerStyle(24);   ratio_sub ->SetMarkerSize(0.8);
  ratio_qcd ->SetXTitle("M_{T2}");   ratio_qcd ->SetYTitle("ratio");
  ratio_sub ->SetXTitle("M_{T2}");   ratio_sub ->SetYTitle("ratio");

  TF1 *f_qcd = new TF1("f_qcd","pol0(0)+expo(1)",bins[0], bins[nbins]);   f_qcd->SetLineColor(8);
  TF1 *f_sub = new TF1("f_sub","pol0(0)+expo(1)",bins[0], bins[nbins]);   f_sub->SetLineColor(8);

  gStyle->SetPalette(1);
  gStyle->SetOptFit (0);

  TCanvas *can = new TCanvas("can", "2d qcd", 0, 0, 900, 700);
  can->SetLogz(1);
  h_ABCD_MT2_qcd->SetMinimum(0.0001);
  h_ABCD_MT2_qcd->Draw("colz");
  TCanvas *ccan = new TCanvas("ccan", "2d susy", 0, 0, 900, 700);
  ccan->SetLogz(1);
  h_ABCD_MT2_susy->SetMinimum(0.0001);
  h_ABCD_MT2_susy->Draw("colz");
  TCanvas *can2 = new TCanvas("can2", "qcd", 0, 0, 900, 700);
  can2->SetLogy(1);
  h_ABCD_lower_y_band_qcd ->Draw();  
  h_ABCD_upper_y_band_qcd ->Draw("same");
  TCanvas *can3 = new TCanvas("can3", "data w/ and w/o substraction", 0, 0, 900, 700);
  can3->SetLogy(1);
  h_ABCD_lower_y_band_data->Draw();  
  h_ABCD_upper_y_band_data->Draw("same");
  h_ABCD_lower_y_band_sub ->Draw("same");  
  h_ABCD_upper_y_band_sub ->Draw("same");
  TCanvas *can4 = new TCanvas("can4", "ratio qcd", 0, 0, 900, 700);
  can4->SetLogy(1);
  ratio_qcd->Draw("E1");
  f_qcd->FixParameter(0,0.006);
  ratio_qcd->Fit("f_qcd","0","",fit_min,fit_max);
  f_qcd->Draw("same");
  TLine *l_qcd_min = new TLine(fit_min,f_qcd->GetYaxis()->GetXmin(),fit_min,f_qcd->GetMaximum()); l_qcd_min->SetLineStyle(2); l_qcd_min->Draw("same");
  TLine *l_qcd_max = new TLine(fit_max,f_qcd->GetYaxis()->GetXmin(),fit_max,f_qcd->GetMaximum()); l_qcd_max->SetLineStyle(2); l_qcd_max->Draw("same");

  TCanvas *can5 = new TCanvas("can5", "ratio data", 0, 0, 900, 700);
  can5->SetLogy(1);
  ratio_sub ->Draw("PE1");
  ratio_data->Draw("sameE1");
  //float p0 = f_qcd->GetParameter(0);
  f_sub->FixParameter(0,0.006);
  ratio_sub->Fit("f_sub","0","",fit_min,fit_max);
  f_sub->Draw("same");
  TLine *l_sub_min = new TLine(fit_min,f_sub->GetYaxis()->GetXmin(),fit_min,f_sub->GetMaximum()); l_sub_min->SetLineStyle(2); l_sub_min->Draw("same");
  TLine *l_sub_max = new TLine(fit_max,f_sub->GetYaxis()->GetXmin(),fit_max,f_sub->GetMaximum()); l_sub_max->SetLineStyle(2); l_sub_max->Draw("same");


  PrintABCDPredictions(var, basecut, upper_cut, lower_cut, f_qcd,f_sub);
  
}

//_________________________________________________________________________________
void MassPlotter::PrintABCDPredictions(TString var, TString basecut, TString upper_cut, TString lower_cut, TF1* func_qcd, TF1* func_sub){  
  TH1D *h_pred_lower_y_band_data = new TH1D("pred_lower_y_band_"+var+"_data" , "", 180,100,1000);		        h_pred_lower_y_band_data->Sumw2();
  TH1D *h_pred_upper_y_band_data = new TH1D("pred_upper_y_band_"+var+"_data" , "", 180,100,1000);		        h_pred_upper_y_band_data->Sumw2();
  TH1D *h_pred_lower_y_band_qcd  = new TH1D("pred_lower_y_band_"+var+"_qcd"  , "", 180,100,1000);		        h_pred_lower_y_band_qcd ->Sumw2();
  TH1D *h_pred_upper_y_band_qcd  = new TH1D("pred_upper_y_band_"+var+"_qcd"  , "", 180,100,1000);		        h_pred_upper_y_band_qcd ->Sumw2();
  TH1D *h_pred_lower_y_band_mc   = new TH1D("pred_lower_y_band_"+var+"_mc"   , "", 180,100,1000);		        h_pred_lower_y_band_mc  ->Sumw2();
  TH1D *h_pred_upper_y_band_mc   = new TH1D("pred_upper_y_band_"+var+"_mc"   , "", 180,100,1000);		        h_pred_upper_y_band_mc  ->Sumw2();
  TH1D *h_pred_lower_y_band_sub  = new TH1D("pred_lower_y_band_"+var+"_sub"  , "", 180,100,1000);		        h_pred_lower_y_band_sub ->Sumw2();
  TH1D *h_pred_upper_y_band_sub  = new TH1D("pred_upper_y_band_"+var+"_sub"  , "", 180,100,1000);		        h_pred_upper_y_band_sub ->Sumw2();
  TF1 *f_qcd = new TF1("f_qcd","pol0(0)+expo(1)",100,1000);
  TF1 *f_sub = new TF1("f_sub","pol0(0)+expo(1)",100,1000);
  f_qcd->SetParameters(func_qcd->GetParameter(0),func_qcd->GetParameter(1),func_qcd->GetParameter(2));
  f_sub->SetParameters(func_sub->GetParameter(0),func_sub->GetParameter(1),func_sub->GetParameter(2));

  for(size_t i = 0; i < fSamples.size(); ++i){
   
    Long64_t nentries =  fSamples[i].tree->GetEntries();
    Long64_t nbytes = 0, nb = 0;
    
    Double_t weight = fSamples[i].xsection * fSamples[i].kfact * fSamples[i].lumi / (fSamples[i].nevents);
    if(fVerbose>2) cout << "ABCD: looping over " << fSamples[i].name << endl;
    if(fVerbose>2) cout << "      sample has weight " << weight << endl; 

    TString selection = TString::Format("(%f) * (%s)"      ,weight,basecut.Data());
    TString sel_up    = TString::Format("(%f) * (%s && %s)",weight,basecut.Data(),upper_cut.Data());
    TString sel_lo    = TString::Format("(%f) * (%s && %s)",weight,basecut.Data(),lower_cut.Data());

    TString variable;
    TString var2 = "misc.MT2";
    if (fSamples[i].type == "data"){
      variable  = TString::Format("%s>>+%s",var2.Data(),h_pred_lower_y_band_data->GetName());
      int nevL = fSamples[i].tree->Draw(variable,sel_lo,"goff");
      variable  = TString::Format("%s>>+%s",var2.Data(),h_pred_upper_y_band_data->GetName());
      int nevU = fSamples[i].tree->Draw(variable,sel_up,"goff");
    }
    else if (fSamples[i].type == "mc" && fSamples[i].name == "QCD_Pt_120to170"){   // WARNING: checking statistical fluctuation in this pt bin!!!
      TString sel_u    = TString::Format("(%f) * (%s && misc.MT2<100&& %s)",weight,basecut.Data(),upper_cut.Data());
      TString sel_l    = TString::Format("(%f) * (%s && misc.MT2<100&& %s)",weight,basecut.Data(),lower_cut.Data());
      variable  = TString::Format("%s>>+%s",var2.Data(),h_pred_lower_y_band_qcd->GetName());
      int nevL = fSamples[i].tree->Draw(variable,sel_l,"goff");
      variable  = TString::Format("%s>>+%s",var2.Data(),h_pred_upper_y_band_qcd->GetName());
      int nevU = fSamples[i].tree->Draw(variable,sel_u,"goff");
    }
    else if (fSamples[i].type == "mc" && fSamples[i].sname == "QCD"){
      variable  = TString::Format("%s>>+%s",var2.Data(),h_pred_lower_y_band_qcd->GetName());
      int nevL = fSamples[i].tree->Draw(variable,sel_lo,"goff");
      variable  = TString::Format("%s>>+%s",var2.Data(),h_pred_upper_y_band_qcd->GetName());
      int nevU = fSamples[i].tree->Draw(variable,sel_up,"goff");
    }
    else if (fSamples[i].type == "mc"){
      variable  = TString::Format("%s>>+%s",var2.Data(),h_pred_lower_y_band_mc->GetName());
      int nevL = fSamples[i].tree->Draw(variable,sel_lo,"goff");
      variable  = TString::Format("%s>>+%s",var2.Data(),h_pred_upper_y_band_mc->GetName());
      int nevU = fSamples[i].tree->Draw(variable,sel_up,"goff");
    }
  }
  
  
  h_pred_lower_y_band_sub->Add(h_pred_lower_y_band_data,h_pred_lower_y_band_mc,1,-1);
  h_pred_upper_y_band_sub->Add(h_pred_upper_y_band_data,h_pred_upper_y_band_mc,1,-1);

  for (int i= 0; i<=181; i++){
    if (h_pred_lower_y_band_sub->GetBinContent(i)<0)  h_pred_lower_y_band_sub->SetBinContent(i, 0.);
    if (h_pred_upper_y_band_sub->GetBinContent(i)<0)  h_pred_upper_y_band_sub->SetBinContent(i, 0.);
  }

  TH1D *h_pred_lower_y_band_data_2 = (TH1D*)h_pred_lower_y_band_data->Clone("h_pred_lower_y_band_data_2");       h_pred_lower_y_band_data->Sumw2();
  TH1D *h_pred_lower_y_band_qcd_2  = (TH1D*)h_pred_lower_y_band_qcd ->Clone("h_pred_lower_y_band_qcd_2" );       h_pred_lower_y_band_qcd ->Sumw2();
  TH1D *h_pred_lower_y_band_sub_2  = (TH1D*)h_pred_lower_y_band_sub ->Clone("h_pred_lower_y_band_sub_2" );       h_pred_lower_y_band_sub ->Sumw2();

  h_pred_lower_y_band_data  ->Multiply(f_qcd);
  h_pred_lower_y_band_qcd   ->Multiply(f_qcd);
  h_pred_lower_y_band_sub   ->Multiply(f_qcd);
  h_pred_lower_y_band_data_2->Multiply(f_sub);
  h_pred_lower_y_band_qcd_2 ->Multiply(f_sub);
  h_pred_lower_y_band_sub_2 ->Multiply(f_sub);

  double yield;
  double error;
  cout << "QCD prediction:" << endl;
  yield = Util::IntegralAndError(h_pred_upper_y_band_qcd, 1,181, error);
  cout << "MT2>100 = " << yield  << " +/- " << error;
  yield = Util::IntegralAndError(h_pred_upper_y_band_qcd,11,181, error);
  cout << "\tMT2>150 = " << yield  << " +/- " << error;
  yield = Util::IntegralAndError(h_pred_upper_y_band_qcd,21,181, error);
  cout << "\tMT2>200 = " << yield  << " +/- " << error << endl
       << "------------------------------------" << endl;
  cout << "QCD estimation from QCD (fit to QCD):" << endl;
  yield = Util::IntegralAndError(h_pred_lower_y_band_qcd, 1,181, error);
  cout << "\tMT2>100 = " << yield  << " +/- " << error;
  yield = Util::IntegralAndError(h_pred_lower_y_band_qcd,11,181, error);
  cout << "\tMT2>150 = " << yield  << " +/- " << error;
  yield = Util::IntegralAndError(h_pred_lower_y_band_qcd,21,181, error);
  cout << "\tMT2>200 = " << yield  << " +/- " << error << endl
       << "------------------------------------" << endl;
  cout << "QCD estimation from QCD (fit to data):" << endl;
  yield = Util::IntegralAndError(h_pred_lower_y_band_qcd_2, 1,181, error);
  cout << "\tMT2>100 = " << yield  << " +/- " << error;
  yield = Util::IntegralAndError(h_pred_lower_y_band_qcd_2,11,181, error);
  cout << "\tMT2>150 = " << yield  << " +/- " << error;
  yield = Util::IntegralAndError(h_pred_lower_y_band_qcd_2,21,181, error);
  cout << "\tMT2>200 = " << yield  << " +/- " << error << endl
       << "------------------------------------" << endl;
  cout << "QCD estimation from data (fit to QCD):" << endl;
  yield = Util::IntegralAndError(h_pred_lower_y_band_data, 1,181, error);
  cout << "\tMT2>100 = " << yield  << " +/- " << error;
  yield = Util::IntegralAndError(h_pred_lower_y_band_data,11,181, error);
  cout << "\tMT2>150 = " << yield  << " +/- " << error;
  yield = Util::IntegralAndError(h_pred_lower_y_band_data,21,181, error);
  cout << "\tMT2>200 = " << yield  << " +/- " << error << endl
       << "------------------------------------" << endl;
  cout << "QCD estimation from data (fit to data):" << endl;
  yield = Util::IntegralAndError(h_pred_lower_y_band_data_2, 1,181, error);
  cout << "\tMT2>100 = " << yield  << " +/- " << error;
  yield = Util::IntegralAndError(h_pred_lower_y_band_data_2,11,181, error);
  cout << "\tMT2>150 = " << yield  << " +/- " << error;
  yield = Util::IntegralAndError(h_pred_lower_y_band_data_2,21,181, error);
  cout << "\tMT2>200 = " << yield  << " +/- " << error << endl
       << "------------------------------------" << endl;
  cout << "QCD estimation from EWK subtracted data (fit to QCD):" << endl;
  yield = Util::IntegralAndError(h_pred_lower_y_band_sub, 1,181, error);
  cout << "\tMT2>100 = " << yield  << " +/- " << error;
  yield = Util::IntegralAndError(h_pred_lower_y_band_sub,11,181, error);
  cout << "\tMT2>150 = " << yield  << " +/- " << error;
  yield = Util::IntegralAndError(h_pred_lower_y_band_sub,21,181, error);
  cout << "\tMT2>200 = " << yield  << " +/- " << error << endl
       << "------------------------------------" << endl;
  cout << "QCD estimation from EWK subtracted data (fit to data):" << endl;
  yield = Util::IntegralAndError(h_pred_lower_y_band_sub_2, 1,181, error);
  cout << "\tMT2>100 = " << yield  << " +/- " << error;
  yield = Util::IntegralAndError(h_pred_lower_y_band_sub_2,11,181, error);
  cout << "\tMT2>150 = " << yield  << " +/- " << error;
  yield = Util::IntegralAndError(h_pred_lower_y_band_sub_2,21,181, error);
  cout << "\tMT2>200 = " << yield  << " +/- " << error << endl
       << "------------------------------------" << endl;

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
	
	float border = 0.3;
 	float scale = (1-border)/border;
 
 	TPad *p_plot  = new TPad("plotpad",  "Pad containing the overlay plot", 0.00, border, 1.00, 1.00, 0, 0);
 	p_plot->SetBottomMargin(0);
 	p_plot->Draw();
 	TPad *p_ratio = new TPad("ratiopad", "Pad containing the ratio",        0.00, 0.00, 1.00, border, 0, 0);
 	p_ratio->SetTopMargin(0);
 	p_ratio->SetBottomMargin(0.35);
 	p_ratio->Draw();
 
	gPad->SetFillStyle(0);
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

	gPad->SetFillStyle(0);
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
	lat.DrawLatex(0.035, 1.0, "ratio");
	
	// x axis title
	lat.SetTextAngle(0);
	float ypos = xtitle.Contains("slash") || xtitle.Contains("_") || xtitle.Contains("^") ? 0.23 : 0.16;
	lat.DrawLatex(0.9, ypos, xtitle);
	//gPad->SetLogy(1);
	gPad->RedrawAxis();
 	p_ratio ->Draw();
 	c1->Update();

	TString save=name+"_ratio";
	Util::PrintEPS(c1, save, fOutputDir);	

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
	if(stacked){
	  h_data  ->SetMaximum(max);
	  h_mc_sum->SetMaximum(max);
	  h       ->SetMaximum(max);
	}

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

