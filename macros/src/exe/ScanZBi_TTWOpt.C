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

#include "UserCode/pandolf/CommonTools/DrawBase.h"
#include "UserCode/pandolf/CommonTools/StatTools.h"

#include "SSDLPlotter.hh"





float fake_lumi_ = 9100.;





void makeDatacard( TTWZPrediction ttwzpred, float lumiSF, float lumiSF_fake, const std::string& optcutsdir, int iEff );
std::pair<TH1F*,TH1F*> getHistoPassingCuts( TTree* tree, std::vector<std::string> names, std::vector<float> cutsMin, std::vector<float> cutsMax );


int main( int argc, char* argv[] ) {

  if( argc < 2 ) {
    std::cout << "USAGE: ./ScanZBi_TTWOpt [selectionType] [charge (\"plus\" or \"minus\" or \"all\")]" << std::endl;
    exit(999);
  }



  std::string selectionType = "Nov23";
  if( argc>1 ) {
    std::string selectionType_str(argv[1]);
    selectionType = selectionType_str;
  }

  std::string charge = "all";
  if( argc>2 ) {
    std::string charge_str(argv[2]);
    charge = charge_str;
  }

  if(charge!="plus" && charge!="minus" && charge!="all" ) {
    std::cout << "only \"plus\" and \"minus\" and \"all\"  are allowed values for charge." << std::endl;
    exit(777);
  }



  // this is an additional selection
  // that can be set by hand
  // to compare to the opt results
  int min_NJets_h = 3;
  int min_NBJets_h = 1;
  int min_NBJets_med_h = 1;
  float min_ptLept1_h = 40.;
  float min_ptLept2_h = 40.;
  float min_met_h = 0.;
  float min_ht_h = 285.;

  int nEffStep = 10;
  if( min_NJets_h>0 && min_NBJets_h>0 || min_NBJets_med_h>0 ||
      min_ptLept1_h>0 || min_ptLept2_h>0 || 
      min_met_h>0 || min_ht_h>0 ) 
    nEffStep+=1;


  float lumi = 19466.;
  //float lumi = 9100.;


  // this sets the style:
  DrawBase* db = new DrawBase("OPT_ZBi");
  db->set_lumiOnRightSide();
  db->set_lumiNormalization(lumi);

  TPaveText* label_sqrt = db->get_labelSqrt();
    

  //std::string dir = "/shome/mdunser/workspace/CMSSW_5_2_5/src/DiLeptonAnalysis/NTupleProducer/macros/plots/Dec20_muPFIso0p05_elPFIso0p05_jet20_withZveto";
  std::string dir = "/shome/lbaeni/top/CMSSW_5_3_2_patch4/src/DiLeptonAnalysis/NTupleProducer/macros/plots/Feb22/";
  std::string config = dir + "/dumperconfig.cfg";

  SSDLPlotter* plotter = new SSDLPlotter(config);
  //std::string outputdir = "/shome/pandolf/CMSSW_4_2_8/src/DiLeptonAnalysis/NTupleProducer/macros/" + selectionType;
  std::string outputdir = "/shome/pandolf/CMSSW_5_3_7_patch5_OPTttW_2/src/DiLeptonAnalysis/NTupleProducer/macros/";
  plotter->setVerbose(1);
  plotter->fDO_OPT = false;
  plotter->setOutputDir(outputdir);
  plotter->init("/shome/lbaeni/top/CMSSW_5_3_2_patch4/src/DiLeptonAnalysis/NTupleProducer/macros/DataCard_SSDL.dat");
  plotter->readHistos(plotter->fOutputFileName);
  plotter->fillRatios(plotter->fMuData,    plotter->fEGData, 0);
  plotter->fillRatios(plotter->fMCBGMuEnr, plotter->fMCBG, 1);
  plotter->storeWeightedPred(plotter->gRegion[plotter->gBaseRegion]);
  //plotter->doAnalysis();



  std::string optcutsdir = "OPT_ttW/optcuts_" + selectionType + "_" + charge;
  std::string ZBiFileName = optcutsdir + "/ZBiScan.txt";

  ofstream ofs_ZBi(ZBiFileName.c_str());
  ofs_ZBi << "Expected for " << lumi/1000. << " fb-1:" << std::endl;
  ofs_ZBi << "Seff   \tS     \tB +- s(B)\tZBi" << std::endl;

  TGraphErrors* gr_ZBi = new TGraphErrors(0);
  float ZBi_max = 0.;
  float effS_ZBi_max = 0.;
  float effMax = 0.;


  for( unsigned iEff=1; iEff<=nEffStep; ++iEff ) {

    if( iEff==11 ) {
      std::cout << std::endl;
      std::cout << "-> Cross checking this selection: " << std::endl;
      std::cout << "NJets >= " << min_NJets_h << std::endl;
      std::cout << "NBJets >= " << min_NBJets_h << std::endl;
      std::cout << "NBJets_med >= " << min_NBJets_med_h << std::endl;
      std::cout << "pt(Lept1) >= " << min_ptLept1_h << std::endl;
      std::cout << "pt(Lept2) >= " << min_ptLept2_h << std::endl;
      std::cout << "MET >= " << min_met_h << std::endl;
      std::cout << "HT >= " << min_ht_h << std::endl;
    }


    int min_NJets = 0;
    int min_NBJets = 0;
    int min_NBJets_med = 0;
    float min_ptLept1 = 0.;
    float min_ptLept2 = 0.;
    float min_met = 0.;
    float min_ht = 0.;


    if( iEff<=10 ) {

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


      for( unsigned iVar=0; iVar<varNames.size(); iVar++ ) {

        if( varNames[iVar]=="NJ"     ) min_NJets      = cutsMin[iVar];
        if( varNames[iVar]=="NbJ"    ) min_NBJets     = cutsMin[iVar];
        if( varNames[iVar]=="NbJmed" ) min_NBJets_med = cutsMin[iVar];
        if( varNames[iVar]=="pT1"    ) min_ptLept1    = cutsMin[iVar];
        if( varNames[iVar]=="pT2"    ) min_ptLept2    = cutsMin[iVar];
        if( varNames[iVar]=="MET"    ) min_met        = cutsMin[iVar];
        if( varNames[iVar]=="HT"     ) min_ht         = cutsMin[iVar];

      }


    } else { //cross check selection

      std::cout << "-> Checking cross check selection." << std::endl;
      
      min_NJets       = min_NJets_h;
      min_NBJets      = min_NBJets_h;
      min_NBJets_med  = min_NBJets_med_h;
      min_ptLept1     = min_ptLept1_h;
      min_ptLept2     = min_ptLept2_h;
      min_met         = min_met_h;
      min_ht          = min_ht_h;

    }

    

    int charge_int = 0;
    if( charge=="plus")  charge_int = 1;
    if( charge=="minus") charge_int = -1;

    TTWZPrediction ttwzpred =  plotter->makePredictionSignalEvents(min_ht, 10000., min_met, 10000., min_NJets, min_NBJets, min_NBJets_med, min_ptLept1, min_ptLept2, charge_int, true);

    float lumi_SF = lumi/plotter->fLumiNorm;
    float lumi_SF_fake = lumi/fake_lumi_;

    std::cout << "lumi_SF: "  << lumi_SF << std::endl;
    std::cout << "lumi_SF_fake: " <<  lumi_SF_fake << std::endl;

    makeDatacard( ttwzpred, lumi_SF, lumi_SF_fake, optcutsdir, iEff );

    float b_pred_mm = (ttwzpred.rare_mm+ttwzpred.wz_mm+ttwzpred.ttz_mm)*lumi_SF + (ttwzpred.fake_mm                 )*lumi_SF_fake;
    float b_pred_em = (ttwzpred.rare_em+ttwzpred.wz_em+ttwzpred.ttz_em)*lumi_SF + (ttwzpred.fake_em+ttwzpred.cmid_em)*lumi_SF_fake;
    float b_pred_ee = (ttwzpred.rare_ee+ttwzpred.wz_ee+ttwzpred.ttz_ee)*lumi_SF + (ttwzpred.fake_ee+ttwzpred.cmid_ee)*lumi_SF_fake;

    float b_pred_mm_err = (ttwzpred.rare_err_mm+ttwzpred.wz_err_mm)*lumi_SF + (ttwzpred.fake_err_mm                     )*lumi_SF_fake;
    float b_pred_em_err = (ttwzpred.rare_err_em+ttwzpred.wz_err_em)*lumi_SF + (ttwzpred.fake_err_em+ttwzpred.cmid_err_em)*lumi_SF_fake;
    float b_pred_ee_err = (ttwzpred.rare_err_ee+ttwzpred.wz_err_ee)*lumi_SF + (ttwzpred.fake_err_ee+ttwzpred.cmid_err_ee)*lumi_SF_fake;

    float s_mm = ttwzpred.ttw_mm*lumi_SF;
    float s_em = ttwzpred.ttw_em*lumi_SF;
    float s_ee = ttwzpred.ttw_ee*lumi_SF;

    float b_pred = b_pred_mm + b_pred_em + b_pred_ee;
    float b_pred_err = sqrt( b_pred_mm_err*b_pred_mm_err + b_pred_em_err*b_pred_em_err + b_pred_ee_err*b_pred_ee_err );
    float s = s_mm + s_em + s_ee;
    float obs = b_pred + s;

    float ZBi = StatTools::computeZBi( obs, b_pred, b_pred_err );


    TH1F* h1_signal = new TH1F("signal", "", 3, -0.5, 2.5);
    h1_signal->SetBinContent( 1, s_mm+b_pred_mm );
    h1_signal->SetBinContent( 2, s_em+b_pred_em );
    h1_signal->SetBinContent( 3, s_ee+b_pred_ee );

    TH1F* h1_bg = new TH1F("bg", "", 3, -0.5, 2.5);
    h1_bg->SetBinContent( 1, b_pred_mm );
    h1_bg->SetBinContent( 2, b_pred_em );
    h1_bg->SetBinContent( 3, b_pred_ee );
    h1_bg->SetBinError( 1, b_pred_mm_err );
    h1_bg->SetBinError( 2, b_pred_em_err );
    h1_bg->SetBinError( 3, b_pred_ee_err );
    h1_bg->SetLineWidth(1);
    h1_bg->SetFillColor(12);
    h1_bg->SetFillStyle(3005);


    h1_signal->SetFillColor( 46 );
    h1_bg->SetFillColor( 38 );


    float effS = (s_mm+s_em+s_ee)/(0.232*lumi*0.22*0.22*0.67);

    if( effS > effMax )
      effMax = effS;

    if( iEff<=10 )
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
    sprintf( canvasName, "%s/yieldPlot_Seff%d.eps", optcutsdir.c_str(), iEff*10);

    //TPaveText* label = new TPaveText( 0.15, 0.65, 0.45, 0.85, "brNDC");
    //label->SetFillColor(0);
    //label->SetTextSize(0.035);
    //label->AddText("L = 5 fb^{-1}");
    //char signalLabel[100];
    //sprintf( signalLabel, "s = %.2f (%d%%)", s, (int)(((float)h1_signal->GetEntries()/nTotal_s)*100) );
    //label->AddText( signalLabel );
    //char bgLabel[100];
    //sprintf( bgLabel, "b = %.2f", b_pred);
    //label->AddText( bgLabel );
    //char signifLabel[100];
    //sprintf( signifLabel, "ZBi = %.2f", ZBi);
    //label->AddText( signifLabel );

    char b_text[500];
    sprintf( b_text, "BG = %.2f #pm %.2f", b_pred, b_pred_err );
    char obs_text[500];
    sprintf( obs_text, "OBS = %.2f", obs );
    char Zbi_text[100];
    sprintf( Zbi_text, "ZBi = %.3f", ZBi );
    TPaveText* label_ZBi = new TPaveText( 0.23, 0.73, 0.5, 0.88, "brNDC" );
    label_ZBi->SetFillColor(0);
    label_ZBi->SetTextSize(0.035);
    label_ZBi->AddText(b_text);
    label_ZBi->AddText(obs_text);
    label_ZBi->AddText(Zbi_text);

    TCanvas* c1 = new TCanvas("c1", "c1", 600., 600.);
    c1->cd();
    h2_axes->Draw();
    //h1_signal->Draw("same");
    stack->Draw("histo same");
    h1_bg->Draw("0 E2 same");
    legend->Draw("same");
    //label->Draw("same");
    label_ZBi->Draw("same");
    label_sqrt->Draw("same");
    gPad->RedrawAxis();
    c1->SaveAs(canvasName);

    delete c1;
    delete legend;
    delete h2_axes;
    //delete stack;
    

    ofs_ZBi << effS << "\t" << s << "\t" << b_pred << " +- " << b_pred_err << "\t" << ZBi << std::endl;

    delete h1_signal;
    delete h1_bg;

std::cout << "### " << iEff << "   ZBi: " << ZBi << std::endl;
  } // for iEff

  std::cout << "> > >   BEST ZBi: " << ZBi_max << std::endl;
  std::cout << "> > >   signal eff: " << effS_ZBi_max << std::endl;

  ofs_ZBi.close();

  db->resetStyle();

  gr_ZBi->SetMarkerSize(2.);
  gr_ZBi->SetMarkerStyle(21);
  gr_ZBi->SetMarkerColor(kRed+3);


  TH2D* h2_axes_gr = new TH2D("axes_gr", "", 10, 0., 1.1*effMax*100., 10, 0., 1.6*ZBi_max ); 
  //TH2D* h2_axes_gr = new TH2D("axes_gr", "", 10, 0., 1., 10, 0., 5.);
  h2_axes_gr->SetYTitle("Expected ZBi (20 fb^{-1})");
  h2_axes_gr->SetXTitle("Signal Efficiency [%]");


  TCanvas* c_gr = new TCanvas("c_gr", "c_gr", 600., 600.);
  c_gr->cd();

  
  h2_axes_gr->Draw();
  gr_ZBi->Draw("P same");
  label_sqrt->Draw("same");

  char ZBi_vs_Seff_name[250];
  sprintf(ZBi_vs_Seff_name, "%s/ZBi_vs_Seff_%s.eps", optcutsdir.c_str(), charge.c_str() );
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

  std::string signal_selection = base_selection + " ( SName==\"TTbarW\" ) )";
  std::string bg_selection     = base_selection + " ( SName!=\"TTbarW\" ) )";

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



void makeDatacard( TTWZPrediction ttwzpred, float lumiSF, float lumiSF_fake, const std::string& optcutsdir, int iEff ) {
  std::cout << "ttwzpred.fake_mm: " <<  ttwzpred.fake_mm << std::endl;
  std::cout << "ttwzpred.fake_em: " <<  ttwzpred.fake_em << std::endl;
  std::cout << "ttwzpred.fake_ee: " <<  ttwzpred.fake_ee << std::endl;

  // this is the reference datacard:
  std::string refDatacard = "OPT_ttW/datacard_ssdl_3channels_ttW.txt";
  std::ifstream ifs(refDatacard.c_str());
  

  // this is the output datacard:
  char datacardName[500];
  sprintf( datacardName, "%s/datacard_eff%d.txt", optcutsdir.c_str(), iEff*10 );
  std::ofstream datacard(datacardName);


  // read reference and copy to output, except for rate line:
  while( ifs.good() )  {
    char thisLine[1024];
    ifs.getline(thisLine, 1024);
    TString thisLine_tstr(thisLine);
    if( thisLine_tstr.BeginsWith("rate") ) {
      datacard << Form("rate\t\t%5.3f\t\t%5.3f\t\t%5.3f\t\t%5.3f\t\t%5.3f\t\t%5.3f\t\t%5.3f\t\t%5.3f\t\t%5.3f\t\t%5.3f\t\t%5.3f\t\t%5.3f\t\t%5.3f\t\t%5.3f\t\t%5.3f\t\t%5.3f\t\t%5.3f\t\t%5.3f",
                    ttwzpred.ttw_mm*lumiSF, ttwzpred.ttz_mm*lumiSF, ttwzpred.fake_mm*lumiSF_fake, 0.0                         , ttwzpred.wz_mm*lumiSF, ttwzpred.rare_mm*lumiSF,
                    ttwzpred.ttw_em*lumiSF, ttwzpred.ttz_em*lumiSF, ttwzpred.fake_em*lumiSF_fake, ttwzpred.cmid_em*lumiSF_fake, ttwzpred.wz_em*lumiSF, ttwzpred.rare_em*lumiSF,
                    ttwzpred.ttw_ee*lumiSF, ttwzpred.ttz_ee*lumiSF, ttwzpred.fake_ee*lumiSF_fake, ttwzpred.cmid_ee*lumiSF_fake, ttwzpred.wz_ee*lumiSF, ttwzpred.rare_ee*lumiSF) << std::endl;
    } else {
      datacard << thisLine << std::endl;
    }
  } // while ifs good


  datacard.close();

  std::cout << "-> Created datacard: " << datacardName << std::endl;
    

//datacard <<      "#=========================================================================================" << std::endl;
//datacard <<      "# Systematics table for ttW/Z analysis, same-sign channel, subchannels" << std::endl;
//datacard << Form("# Generated on: %s ", asctime(timeinfo)) << std::endl;
//datacard <<      "# Copy between the dashed lines for datacard" << std::endl;
//datacard <<      "#-----------------------------------------------------------------------------------------" << std::endl;
//datacard <<      "imax 3" << std::endl;
//datacard <<      "jmax 5" << std::endl;
//datacard <<      "kmax *" << std::endl;
//datacard << std::endl << std::endl;
//datacard <<      "bin\t\t1\t2\t3" << std::endl;
//datacard <<  "observation\t10\t10\t10" << std::endl; // doesnt really matter
//datacard << std::endl << std::endl;
//datacard <<      "bin\t\t1\t\t1\t\t1\t\t1\t\t1\t\t1\t\t2\t\t2\t\t2\t\t2\t\t2\t\t2\t\t3\t\t3\t\t3\t\t3\t\t3\t\t3" << std::endl;
//datacard <<      "process\t\tttW\t\tttZ\t\tfake\t\tcmid\t\twz\t\trare\t\tttW\t\tttZ\t\tfake\t\tcmid\t\twz\t\trare\t\tttW\t\tttZ\t\tfake\t\tcmid\t\twz\t\trare" << std::endl;
//datacard <<      "process\t\t0\t\t1\t\t2\t\t3\t\t4\t\t5\t\t0\t\t1\t\t2\t\t3\t\t4\t\t5\t\t0\t\t1\t\t2\t\t3\t\t4\t\t5" << std::endl;

//datacard << Form("rate\t\t%5.3f\t\t%5.3f\t\t%5.3f\t\t%5.3f\t\t%5.3f\t\t%5.3f\t\t%5.3f\t\t%5.3f\t\t%5.3f\t\t%5.3f\t\t%5.3f\t\t%5.3f\t\t%5.3f\t\t%5.3f\t\t%5.3f\t\t%5.3f\t\t%5.3f\t\t%5.3f",
//              ttwzpred.ttw_mm*lumiSF, ttwzpred.ttz_mm*lumiSF, ttwzpred.fake_mm*lumiSF, 0.0                    , ttwzpred.wz_mm*lumiSF, ttwzpred.rare_mm*lumiSF,
//              ttwzpred.ttw_em*lumiSF, ttwzpred.ttz_em*lumiSF, ttwzpred.fake_em*lumiSF, ttwzpred.cmid_em*lumiSF, ttwzpred.wz_em*lumiSF, ttwzpred.rare_em*lumiSF,
//              ttwzpred.ttw_ee*lumiSF, ttwzpred.ttz_ee*lumiSF, ttwzpred.fake_ee*lumiSF, ttwzpred.cmid_ee*lumiSF, ttwzpred.wz_ee*lumiSF, ttwzpred.rare_ee*lumiSF) << std::endl;
//datacard << std::endl << std::endl;
//datacard <<      "#syst" << std::endl;
//float lumiError = 1.022;
//datacard <<      "lumi     lnN\t"+lumiError+"\t\t"+lumiError+"\t\t"+lumiError+"\t\t"+lumiError+"\t\t"+lumiError+"\t\t"+lumiError+"\t\t"+lumiError+"\t\t"+lumiError+"\t\t"+lumiError+"\t\t"+lumiError+"\t\t"+lumiError+"\t\t"+lumiError+"\t\t"+lumiError+"\t\t"+lumiError+"\t\t"+lumiError+"\t\t"+lumiError+"\t\t"+lumiError+"\t\t"+lumiError+"" << std::endl;

//datacard << "bgUncfak lnN\t-\t\t-\t\t2.084\t\t-\t\t -    \t\t-    \t\t-\t\t-\t\t 1.564\t\t-    \t\t-    \t\t-    \t\t-\t\t-\t\t 1.681\t\t-    \t\t-    \t\t-    " << std::endl;
//datacard << "bgUnccmi lnN\t-\t\t-\t\t-    \t\t-\t\t -    \t\t-    \t\t-\t\t-\t\t -    \t\t1.060\t\t-    \t\t-    \t\t-\t\t-\t\t -    \t\t1.090\t\t-    \t\t-    " << std::endl;
//datacard << "bgUncwz  lnN\t-\t\t-\t\t-    \t\t-\t\t 1.271\t\t-    \t\t-\t\t-\t\t -    \t\t-    \t\t1.198\t\t-    \t\t-\t\t-\t\t -    \t\t-    \t\t1.278\t\t-    " << std::endl;
//datacard << "bgUncrar lnN\t-\t\t-\t\t-    \t\t-\t\t -    \t\t2.065\t\t-\t\t-\t\t -    \t\t-    \t\t-    \t\t1.673\t\t-\t\t-\t\t -    \t\t-    \t\t-    \t\t1.781" << std::endl;

//datacard << "lept     lnN\t0.938/1.014\t0.938/1.014\t	-    		-		0.974/1.032	0.960/1.072	0.975/1.022	0.975/1.022	-		-		0.959/1.061	0.990/1.019	0.982/1.018	0.982/1.018	-		-		1.000/1.000	0.985/1.013
//datacard << "btag     lnN\t0.967/1.013\t0.967/1.013\t	-    		-		0.918/1.094	0.992/1.082	0.989/1.012	0.989/1.012	-		-		0.975/1.037	0.974/0.997	0.973/1.015	0.973/1.015	-		-		0.995/1.028	0.991/1.007
//datacard << "jes      lnN\t0.968/1.067\t0.968/1.067\t	-    		-		0.862/1.189	0.942/1.058	0.966/1.065	0.966/1.065	-		-		0.939/1.176	0.950/1.066	0.970/1.088	0.970/1.088	-		-		0.886/1.235	0.988/1.018
//datacard << "jer      lnN\t1.012      \t1.012		-    		-		0.996		0.990		0.993		0.993		-		-		1.085		1.024		0.987		0.987		-		-		0.995		1.000
//datacard << "pu       lnN\t1.030      \t1.030		-    		-		1.030		1.030		1.030		1.030		-		-		1.030		1.030		1.030		1.030		-		-		1.030		1.030
//datacard << "matching lnN\t1.015/0.998\t1.015/0.998	-    		-		1.015/0.998	1.015/0.998	1.015/0.998	1.015/0.998	-		-		1.015/0.998	1.015/0.998	1.015/0.998	1.015/0.998	-		-		1.015/0.998	1.015/0.998
//datacard << "scale    lnN\t1.023/0.966\t1.023/0.966	-    		-		1.023/0.966	1.023/0.966	1.023/0.966	1.023/0.966	-		-		1.023/0.966	1.023/0.966	1.023/0.966	1.023/0.966	-		-		1.023/0.966	1.023/0.966


}
