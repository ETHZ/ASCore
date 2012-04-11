#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TFile.h"
#include "TH2D.h"
#include "TChain.h"
#include "THStack.h"
#include "TLegend.h"
#include "TPaveText.h"
#include <iostream>
#include <fstream>

#include "DrawBase.h"

#include "../include/SSDLPlotter.hh"


std::pair<TH1F*,TH1F*> getHistoPassingCuts( TTree* tree, std::vector<std::string> names, std::vector<float> cutsMin, std::vector<float> cutsMax );
float computeZBi( float obs, float b_pred, float b_pred_err );


int main( int argc, char* argv[] ) {

  std::string selectionType = "Apr10_Iso005_NoZVeto_jet20";
  if( argc>1 ) {
    std::string selectionType_str(argv[1]);
    selectionType = selectionType_str;
  }


  DrawBase* db = new DrawBase("OPT_ZBi");

  SSDLPlotter* plotter = new SSDLPlotter();
  std::string outputdir = "/shome/pandolf/CMSSW_4_2_8/src/DiLeptonAnalysis/NTupleProducer/macros/" + selectionType;
  plotter->setOutputDir(outputdir);
  plotter->init("/shome/pandolf/CMSSW_4_2_8/src/DiLeptonAnalysis/NTupleProducer/macros/DataCard_SSDL.dat");


  //TFile* yieldsFile = TFile::Open("/shome/pandolf/CMSSW_4_2_8/src/DiLeptonAnalysis/NTupleProducer/macros/Apr5_Iso015_NoZVeto_3jets/SSDLYields.root");
  TFile* yieldsFile = TFile::Open("/shome/pandolf/CMSSW_4_2_8/src/DiLeptonAnalysis/NTupleProducer/macros/OPT_ttW/opt_ttW.root");

  //TTree* sigEvents = (TTree*)yieldsFile->Get("SigEvents");
  TTree* sigEvents = (TTree*)yieldsFile->Get("tree_opt");

  std::string optcutsdir = "optcuts_" + selectionType;
  std::string ZBiFileName = optcutsdir + "/ZBiScan.txt";

  ofstream ofs_ZBi(ZBiFileName.c_str());
  ofs_ZBi << "Expected for 5 fb-1:" << std::endl;
  ofs_ZBi << "Seff   \tS     \tB  \tZBi" << std::endl;

  TGraphErrors* gr_ZBi = new TGraphErrors(0);
  float ZBi_max = 0.;
  float effS_ZBi_max = 0.;
  float effMax = 0.;

  float nTotal_s = 1089608.; //hardwired!!! HORRIBLE!!

  for( unsigned iEff=1; iEff<10; ++iEff ) {

    char infileName[300];
    sprintf( infileName, "%s/cuts_Seff%d.txt", optcutsdir.c_str(), iEff*10);
    ifstream ifs(infileName);
    std::cout << "-> Opening Seff file: " << infileName << std::endl;
  
    std::vector<std::string> varNames;
    std::vector<float> cutsMin;
    std::vector<float> cutsMax;

    while( ifs.good() && !ifs.eof() ) {

      std::string varName;
      float cutMin, cutMax;

      ifs >> varName >> cutMin >> cutMax;

      varNames.push_back( varName );
      cutsMin.push_back( cutMin );
      cutsMax.push_back( cutMax );

    } //while file is good
  
    ifs.close();

    // eliminate last element (last line is read and is empty):
    varNames.pop_back();
    cutsMin.pop_back();
    cutsMax.pop_back();

    std::pair<TH1F*, TH1F*> sig_bg = getHistoPassingCuts( sigEvents, varNames, cutsMin, cutsMax );
    TH1F* h1_signal = sig_bg.first;
    TH1F* h1_bg  = sig_bg.second;

    h1_signal->SetFillColor( 46 );
    h1_bg->SetFillColor( 38 );


    float s = h1_signal->Integral(0, h1_signal->GetNbinsX());
    float b = h1_bg->Integral(0, h1_bg->GetNbinsX());


    int min_NJets = 0;
    int min_NBJets = 0;
    int min_NBJets_med = 0;
    float min_ptLept1 = 0.;
    float min_ptLept2 = 0.;
    float min_met = 0.;

    for( unsigned iVar=0; iVar<varNames.size(); iVar++ ) {

      if( varNames[iVar]=="NJ"     ) min_NJets      = cutsMin[iVar];
      if( varNames[iVar]=="NbJ"    ) min_NBJets     = cutsMin[iVar];
      if( varNames[iVar]=="NbJmed" ) min_NBJets_med = cutsMin[iVar];
      if( varNames[iVar]=="pT1"    ) min_ptLept1    = cutsMin[iVar];
      if( varNames[iVar]=="pT2"    ) min_ptLept2    = cutsMin[iVar];
      if( varNames[iVar]=="MET"    ) min_met        = cutsMin[iVar];

    }


    SSDLPrediction ssdlpred =  plotter->makePredictionSignalEvents(0., 10000., 0., 10000., min_NJets, min_NBJets, min_NBJets_med, min_ptLept1, min_ptLept2);

    float b_pred = ssdlpred.bg_mm + ssdlpred.bg_em + ssdlpred.bg_ee;
    float b_pred_err = sqrt( ssdlpred.bg_mm_err*ssdlpred.bg_mm_err + ssdlpred.bg_em_err*ssdlpred.bg_em_err + ssdlpred.bg_ee_err*ssdlpred.bg_ee_err );
    float obs = b_pred + ssdlpred.s_mm + + ssdlpred.s_em + + ssdlpred.s_ee;

    float ZBi = computeZBi( obs, b_pred, b_pred_err );

    float effS = (float)h1_signal->GetEntries()/nTotal_s;

    if( effS > effMax )
      effMax = effS;

    gr_ZBi->SetPoint( iEff-1, 100.*effS, ZBi );

    if( ZBi > ZBi_max ) {
      ZBi_max = ZBi;
      effS_ZBi_max = effS;
    }

    float yMax = h1_signal->GetMaximum() + h1_bg->GetMaximum();
    yMax*=1.5;

    THStack* stack = new THStack();
    stack->Add( h1_bg );
    stack->Add( h1_signal );

    TH2D* h2_axes = new TH2D("axes", "", 3, -0.5, 2.5, 10, 0., yMax);
    h2_axes->GetXaxis()->SetLabelSize(0.085);
    h2_axes->GetXaxis()->SetBinLabel(1, "#mu#mu");
    h2_axes->GetXaxis()->SetBinLabel(2, "e#mu");
    h2_axes->GetXaxis()->SetBinLabel(3, "ee");
    h2_axes->SetYTitle("Events");


    TLegend* legend = new TLegend(0.6, 0.75, 0.88, 0.88);
    legend->SetFillColor(0);
    legend->SetTextSize(0.035);
    legend->AddEntry( h1_signal, "Signal", "F");
    legend->AddEntry( h1_bg, "Background", "F");

    char canvasName[250];
    sprintf( canvasName, "optcuts/yieldPlot_Seff%d.eps", iEff*10);

    TPaveText* label = new TPaveText( 0.15, 0.65, 0.45, 0.85, "brNDC");
    label->SetFillColor(0);
    label->SetTextSize(0.035);
    label->AddText("L = 5 fb^{-1}");
    char signalLabel[100];
    sprintf( signalLabel, "s = %.2f (%d%%)", s, (int)(((float)h1_signal->GetEntries()/nTotal_s)*100) );
    label->AddText( signalLabel );
    char bgLabel[100];
    sprintf( bgLabel, "b = %.2f", b);
    label->AddText( bgLabel );
    char signifLabel[100];
    sprintf( signifLabel, "ZBi = %.2f", ZBi);
    label->AddText( signifLabel );

    
    TCanvas* c1 = new TCanvas("c1", "c1", 600., 600.);
    c1->cd();
    h2_axes->Draw();
    stack->Draw("histo same");
    legend->Draw("same");
    label->Draw("same");
    //gPad->RedrawAxis();
    c1->SaveAs(canvasName);

    delete c1;
    delete legend;
    delete h2_axes;
    delete stack;
    

    ofs_ZBi << effS << "\t" << s << "\t" << b << "\t" << ZBi << std::endl;

    delete h1_signal;
    delete h1_bg;

std::cout << "### " << iEff << std::endl;
  } // for iEff

  std::cout << "> > >   BEST ZBi: " << ZBi_max << std::endl;
  std::cout << "> > >   signal eff: " << effS_ZBi_max << std::endl;

  ofs_ZBi.close();

  gr_ZBi->SetMarkerSize(2.);
  gr_ZBi->SetMarkerStyle(29);
  gr_ZBi->SetMarkerColor(kRed+3);


  TH2D* h2_axes_gr = new TH2D("axes_gr", "", 10, 0., 1.3*effMax*100., 10, 0., 1.6*ZBi_max ); 
  //TH2D* h2_axes_gr = new TH2D("axes_gr", "", 10, 0., 1., 10, 0., 5.);
  h2_axes_gr->SetYTitle("ZBi (5 fb^{-1})");
  h2_axes_gr->SetXTitle("Signal Efficiency [%]");


  TCanvas* c_gr = new TCanvas("c_gr", "c_gr", 600., 600.);
  c_gr->cd();

  
  h2_axes_gr->Draw();
  gr_ZBi->Draw("P same");

  char ZBi_vs_Seff_name[250];
  sprintf(ZBi_vs_Seff_name, "%s/ZBi_vs_Seff.eps", optcutsdir.c_str() );
  c_gr->SaveAs(ZBi_vs_Seff_name);

}



std::pair<TH1F*,TH1F*> getHistoPassingCuts( TTree* tree, std::vector<std::string> names, std::vector<float> cutsMin, std::vector<float> cutsMax ) {


  std::string base_selection = "eventWeight*(";
  for( unsigned iVar=0; iVar<names.size(); ++iVar) {
    char additionalCut[300];
    sprintf( additionalCut, "%s > %f && %s < %f && ", names[iVar].c_str(), cutsMin[iVar], names[iVar].c_str(), cutsMax[iVar] );
    std::string additionalCut_str(additionalCut);
    base_selection += additionalCut_str;
  }

  std::string signal_selection = base_selection + "  ( SName==\"TTbarW\" || SName==\"TTbarZ\" ) )";
  std::string bg_selection     = base_selection + " !( SName==\"TTbarW\" || SName==\"TTbarZ\" ) )";

  TH1F* h1_signal = new TH1F("signal", "", 3, -0.5, 2.5);
  h1_signal->Sumw2();
  TH1F* h1_bg = new TH1F("bg", "", 3, -0.5, 2.5);
  h1_bg->Sumw2();

  tree->Project( "signal", "Flavor",  signal_selection.c_str() );
  tree->Project( "bg",  "Flavor", bg_selection.c_str() );

  std::pair<TH1F*, TH1F*> returnPair;
  returnPair.first = h1_signal;
  returnPair.second = h1_bg;

  return returnPair;

}




float computeZBi( float obs, float b_pred, float b_pred_err ) {

  float tau = b_pred / ( b_pred_err*b_pred_err );
  float n_off = tau*b_pred;
  float P_Bi = TMath::BetaIncomplete( 1./(1.+tau), obs, n_off+1. );
  float Z_Bi = sqrt(2)*TMath::ErfInverse( 1 - 2.*P_Bi );

  return Z_Bi;

}
