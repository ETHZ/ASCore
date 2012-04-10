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



TH1F* getHistoPassingCuts( std::string histoName, TTree* tree, int nbtags, std::vector<std::string> names, std::vector<float> cutsMin, std::vector<float> cutsMax, float massMin, float massMax );


int main( int argc, char* argv[] ) {


  TFile* yieldsFile = TFile::Open("/shome/pandolf/CMSSW_4_2_8/src/DiLeptonAnalysis/NTupleProducer/macros/Apr5_Iso015_NoZVeto_3jets/SSDLYields.root");

  TTree* sigEvents = (TTree*)yieldsFile->Get("SigEvents");

  std::string ZBiFileName = "optcuts/ZBiScan.txt";

  ofstream ofs_ZBi(ZBiFileName);
  ofs_Zbi << "Expected for 5 fb-1:" << std::endl;
  ofs_ZBi << "Seff   \tS     \tB  \tZBi" << std::endl;

  TGraphErrors* gr_ZBi = new TGraphErrors(0);
  float ZBi_max = 0.;
  float effS_Zbi_max = 0.;


  for( unsigned iEff=1; iEff<10; ++iEff ) {

    char infileName[300];
    sprintf( infileName, "optcuts/cuts_Seff%d.txt", iEff*10);
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
    TH1F* h1_sig = sig_bg.first;
    TH1F* h1_bg  = sig_bg.second;

    h1_signal->SetFillColor( 46 );
    h1_bg->SetFillColor( 38 );


    float s = h1_signal->Integral(0, h1_signal->GetNbinsX());
    float b = h1_bg->Integral(0, h1_bg->GetNbinsX());


    float Zbi = computeZBi();


    float effS = (float)h1_signal->GetEntries()/nTotal_s;


    gr_Zbi->SetPoint( iEff-1, 100.*effS, ZBi );

    if( ZBi > ZBi_max ) {
      ZBi_max = ZBi;
      effS_Zbi_max = effS;
    }

    float yMax = h1_signal->GetMaximum() + h1_bg->GetMaximum();
    yMax*=1.5;

    THStack* stack = new THStack();
    stack->Add( h1_bg );
    stack->Add( h1_signal );

    TH2D* h2_axes = new TH2D("axes", "", 10, massMin, massMax, 10, 0., yMax);
    h2_axes->SetXTitle("ZZ Invariant Mass [GeV/c^{2}]");
    h2_axes->SetYTitle("Events / fb^{-1}");
    h2_axes->GetXaxis()->SetTitleOffset(1.1);
    h2_axes->GetYaxis()->SetTitleOffset(1.5);


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
    c1->SetLeftMargin(0.12);
    c1->SetBottomMargin(0.12);
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
    

    ofs_sign << effS << "\t" << s << "\t" << b << "\t" << ZBi << std::endl;

    delete h1_signal;
    delete h1_bg;

std::cout << "### " << iEff << std::endl;
  } // for iEff

  std::cout << "> > >   BEST ZBi: " << ZBi_max << std::endl;
  std::cout << "> > >   signal eff: " << effS_ZBi_max << std::endl;

  ofs_sign.close();

  graph->SetMarkerSize(2.);
  graph->SetMarkerStyle(29);
  graph->SetMarkerColor(kRed+3);

  graphUL->SetMarkerSize(2.);
  graphUL->SetMarkerStyle(29);
  graphUL->SetMarkerColor(kRed+3);

  graphUL_bg30->SetMarkerSize(2.);
  graphUL_bg30->SetMarkerStyle(20);
  graphUL_bg30->SetMarkerColor(kOrange+1);

  TH2D* h2_axes_gr = new TH2D("axes_gr", "", 10, 0., 1.3*effMax*100., 10, 0., 1.6*signMax ); 
  //TH2D* h2_axes_gr = new TH2D("axes_gr", "", 10, 0., 1., 10, 0., 5.);
  h2_axes_gr->SetYTitle("ZBi (5 fb^{-1})");
  h2_axes_gr->SetXTitle("Signal Efficiency [%]");
  h2_axes_gr->GetXaxis()->SetTitleOffset(1.1);
  h2_axes_gr->GetYaxis()->SetTitleOffset(1.5);


  TCanvas* c_gr = new TCanvas("c_gr", "c_gr", 600., 600.);
  c_gr->SetLeftMargin(0.12);
  c_gr->SetBottomMargin(0.12);
  c_gr->cd();

  
  h2_axes_gr->Draw();
  gr_Zbi->Draw("P same");

  char ZBi_vs_Seff_name[250];
  sprintf(ZBi_vs_Seff_name, "ZBi_vs_Seff.eps" );
  c_gr->SaveAs(ZBi_vs_Seff_name);

}



TH1F* getHistoPassingCuts( TTree* tree, std::vector<std::string> names, std::vector<float> cutsMin, std::vector<float> cutsMax ) {


  TH1F* histo = new TH1F(histoName.c_str(), "", 50, massMin, massMax);
  histo->Sumw2();


  Float_t eventWeight;
  tree->SetBranchAddress( "eventWeight", &eventWeight );

  Float_t absEtaLept1;
  tree->SetBranchAddress( "absEtaLept1", &absEtaLept1 );

  std::vector<float> variables(names.size());
  int index_mZZ=-1;
  int index_mZjj=-1;
  int index_mZll=-1;
  int index_ptLept1=-1;


  for( unsigned i=0; i<names.size(); ++i ) {
    //std::cout << "::getHistoPassingCuts:: Setting Branch Address: '" << names[i] << "'" << std::endl;
    if( names[i]=="QGLikelihoodJet1_T_QGLikelihoodJet2" ) continue;
  //  tree->SetBranchAddress("QGLikelihoodJet1"
  //} else {
      tree->SetBranchAddress( names[i].c_str(), &(variables[i]) );
//std::cout << "set " << names[i] << std::endl;
  //}
    if( names[i]=="mZZ" )
      index_mZZ = i;
    if( names[i]=="mZll" )
      index_mZll = i;
    if( names[i]=="mZjj" )
      index_mZjj = i;
    if( names[i]=="ptLept1" )
      index_ptLept1 = i;
  }

  Int_t leptType;
  tree->SetBranchAddress( "leptType", &leptType );

  Int_t nBTags;
  tree->SetBranchAddress( "nBTags", &nBTags );

  Float_t mZZ;
  if( index_mZZ<0 )
    tree->SetBranchAddress( "mZZ", &mZZ );

  Float_t mZjj;
  if( index_mZjj<0 )
    tree->SetBranchAddress( "mZjj", &mZjj );

  Float_t mZll;
  if( index_mZll<0 )
    tree->SetBranchAddress( "mZll", &mZll );

  Float_t ptLept1;
  if( index_ptLept1<0 )
    tree->SetBranchAddress( "ptLept1", &ptLept1 );


  int nentries = tree->GetEntries();

  //std::cout << "::getHistoPassingCuts:: Begin Loop." << std::endl;

  for( unsigned iEntry=0; iEntry<nentries; ++iEntry) {
 
    tree->GetEntry(iEntry);

    bool pass = true;

    for( unsigned iVar=0; iVar<variables.size() /*&& pass*/; ++iVar) {

 //   // preselection:
 //   if( absEtaLept1>2.1 ) pass=false;

 //   if( names[iVar]=="mZZ" ) {
 //     if( variables[iVar]<190. ) pass=false;
 //   }
 //   else if( names[iVar]=="ptLept1" ) {
 //     if( variables[iVar]<35.0324 ) pass=false;
 //   }
 //   else if( names[iVar]=="deltaRll" ) {
 //     if( variables[iVar]>2.14868 ) pass=false;
 //   }
 //   else if( names[iVar]=="ptJet1" ) {
 //     if( variables[iVar]<34.2011 ) pass=false;
 //   }
 //   else if( names[iVar]=="ptJet2" ) {
 //     if( variables[iVar]<22.5958 ) pass=false;
 //   }
 //   else if( names[iVar]=="mZjj" ) {
 //     if( variables[iVar]<56.9411 || variables[iVar]>112.121 ) pass=false;
 //   }


      // selection:
//std::cout << "requiring " << names[iVar] << ">" << cutsMin[iVar] << " && " << names[iVar] << "<" << cutsMax[iVar] << std::endl;
      if( variables[iVar]<cutsMin[iVar] || variables[iVar]>cutsMax[iVar] )
        pass = false;

    }

    //if( leptType != 1 ) pass = false;

    if( index_ptLept1<0 )
      if( ptLept1<40. ) pass=false;

    if( nBTags != nbtags ) pass = false;


    if( !pass ) continue;

    if( index_mZZ<0 )
      histo->Fill( mZZ, 1000.*eventWeight ); //1 fb-1
    else
      histo->Fill( variables[index_mZZ], 1000.*eventWeight ); //1 fb-1

  } // for entries

  //std::cout << "::getHistoPassingCuts:: End Loop." << std::endl;

  return histo;

} // getHistoPassingCuts
