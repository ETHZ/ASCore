#include <sstream>
#include <string>
#include <iomanip>
#include "ExclusionPlot.hh"
 
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TStyle.h"
#include "TLine.h"
#include "TCanvas.h"
#include "TMarker.h"
#include <vector>
#include "TMath.h"

//Produce the limit plot with the command: root -l ExclusionPlot.C+

TCanvas * ExclusionPlot(){
  gStyle->SetPalette(1);

  gStyle->SetPadTickX(1);  // To get tick marks on the opposite side of the frame
  gStyle->SetPadTickY(1);

  gStyle->SetOptTitle(0);
  gStyle->SetTitleFont(42);
  gStyle->SetTitleColor(1);
  gStyle->SetTitleTextColor(1);
  gStyle->SetTitleFillColor(10);
  gStyle->SetTitleColor(1, "XYZ");
  gStyle->SetTitleFont(42, "XYZ");



  Int_t tanBeta = 10;
  Bool_t plotLO = false;
   
  TCanvas * canvas = CommandMSUGRA("SSDL_msugra_tanb10.root",tanBeta, plotLO);
  return canvas;
}


TCanvas* CommandMSUGRA(TString plotName_X,Int_t tanBeta_, Bool_t plotLO_){
  
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1); 
  gStyle->SetTextFont(42);
  gStyle->SetFrameBorderMode(0);

  //convert tanb value to string
  std::stringstream tmp;
  tmp << tanBeta_;
  std::cout << "my tan beta is " << tanBeta_ << std::endl;
  TString tanb( tmp.str() );
  
  
  // Output file
  //std::cout << " create " << plotName_X << std::endl;
  //TFile* output = new TFile( plotName_X, "RECREATE" );
  //if ( !output || output->IsZombie() ) { std::cout << " zombie alarm output is a zombie " << std::endl; }
  

  //set old exclusion Limits
  TGraph* LEP_ch = set_lep_ch(tanBeta_);
  TGraph* LEP_sl = set_lep_sl(tanBeta_);//slepton curve
  TGraph* TEV_sg_cdf = set_tev_sg_cdf(tanBeta_);//squark gluino cdf
  TGraph* TEV_sg_d0 = set_tev_sg_d0(tanBeta_);//squark gluino d0
  //  TGraph* TEV_tlp_cdf = set_tev_tlp_cdf(tanBeta_);//trilepton cdf
  //  TGraph* TEV_tlp_d0 = set_tev_tlp_d0(tanBeta_);//trilepton d0
  TGraph* stau   = set_tev_stau(tanBeta_);//stau 
  TGraph* NoEWSB = set_NoEWSB(tanBeta_); 

  TGraph* TEV_sn_d0_1 = set_sneutrino_d0_1(tanBeta_);
  TGraph* TEV_sn_d0_2 = set_sneutrino_d0_2(tanBeta_);

  //constant ssqquark and gluino lines
  TF1* lnsq[10];
  TF1* lngl[10];
  
  TLatex* sq_text[10];
  TLatex* gl_text[10];

  
  for(int i = 0; i < 6; i++){
    lnsq[i] = constant_squark(tanBeta_,i);
    sq_text[i] = constant_squark_text(i,*lnsq[i],tanBeta_);
    lngl[i] = constant_gluino(tanBeta_,i);
    gl_text[i] = constant_gluino_text(i,*lngl[i]);
  }
  

  //Legends
  TLegend* legst  = makeStauLegend(0.05,tanBeta_);
  TLegend* legNoEWSB  = makeNoEWSBLegend(0.05,tanBeta_);
  TLegend* legexp = makeExpLegend( *TEV_sg_cdf,*TEV_sg_d0,*LEP_ch,*LEP_sl,*TEV_sn_d0_1,0.035,tanBeta_);
  
 
  //make Canvas
  TCanvas* cvsSys = new TCanvas("cvsnm","cvsnm",0,0,800,600);
  gStyle->SetOptTitle(0);
  cvsSys->SetFillColor(0);
  cvsSys->GetPad(0)->SetRightMargin(0.07);
  cvsSys->Range(-120.5298,26.16437,736.0927,750);
  //  cvsSys->Range(-50.5298,26.16437,736.0927,500);
  cvsSys->SetFillColor(0);
  cvsSys->SetBorderMode(0);
  cvsSys->GetPad(0)->SetBorderMode(0);
  cvsSys->GetPad(0)->SetBorderSize(2);
  cvsSys->GetPad(0)->SetLeftMargin(0.1407035);
  cvsSys->GetPad(0)->SetTopMargin(0.08);
  cvsSys->GetPad(0)->SetBottomMargin(0.13);

  cvsSys->SetTitle("tan#beta="+tanb);
 
  //output->cd();
  
//and now  the exclusion limits


  TGraph* SSdilep;
  TGraphErrors* OSdilep;
  TGraphErrors* RA1;

  TGraphErrors* RA1_old;
  TGraphErrors* RA5_old;
  TGraphErrors* RA6_old;


  
  if (tanBeta_ == 10) {
    std::cout << "tanb ==10" << std::endl;
    //SSdilep = SSdilep_NLO();
    //OSdilep = OSdilep_NLO();
    //RA1 = RA1_NLO();

    //RA1_old = getRA1Observed_NLO_tanBeta10();
    //RA5_old = getRA5Observed_NLO_tanBeta10();
    //RA6_old = getRA6Observed_NLO_tanBeta10();

  }



  double m0min = 0;
  if (tanBeta_ == 50) m0min=200;
  //TH2D* hist = new TH2D("h","h",100,m0min,1000,100,120,700);
  TH2D* hist = new TH2D("h","h",100,m0min,2000,100,120,700);
  hist->Draw();  
  hist->GetXaxis()->SetTitle("m_{0} (GeV)");
  hist->GetYaxis()->SetTitle("m_{1/2} (GeV)");
  hist->GetXaxis()->SetTitleOffset(.9);
  hist->GetXaxis()->SetTitleSize(0.06);
  hist->GetYaxis()->SetTitleOffset(1.0);
  hist->GetYaxis()->SetTitleSize(0.06);

  hist->GetXaxis()->SetNdivisions(506);
  //  if (tanBeta_ == 50)  hist->GetXaxis()->SetNdivisions(504);
  hist->GetYaxis()->SetNdivisions(506);

  int col[]={2,3,4};



  //TLegend* myleg;

  //if( plotLO_ ) myleg = new TLegend(0.3,0.65,0.65,0.8,NULL,"brNDC");
  //else          myleg = new TLegend(0.25,0.76,0.44,0.91,NULL,"brNDC");


  //myleg->SetFillColor(0); 
  //myleg->SetShadowColor(0);
  //myleg->SetTextSize(0.04);
  //myleg->SetBorderSize(0);

  //TLegendEntry *entry=myleg->AddEntry("RA1","2011 Limits","l");
  //entry->SetLineColor(1);
  //entry->SetLineStyle(1);
  //entry->SetLineWidth(3);

  //entry=myleg->AddEntry("sRA1","2010 Limits","l");
  //entry->SetLineColor(1);
  //entry->SetLineStyle(2);
  //entry->SetLineWidth(3);

  
  //constant squark and gluino mass contours
  for (int it=0;it<5;it++) {   
    lngl[it]->Draw("same");   
    lnsq[it]->Draw("same");
    sq_text[it]->Draw();
    gl_text[it]->Draw();
  }
  

  //exclusion limits previous experiments
  if(tanBeta_ == 3){
    TEV_sn_d0_1->Draw("fsame");
    TEV_sn_d0_2->Draw("fsame");
  }
  LEP_ch->Draw("fsame");
  if (tanBeta_ != 50) LEP_sl->Draw("fsame");

  //remove CDF/D0 excluded regions
  TEV_sg_cdf->Draw("fsame");
  TEV_sg_d0->Draw("same");  
  TEV_sg_d0->Draw("fsame");


  //other labels
  Double_t xpos = 0;
  Double_t xposi = 0;
  Double_t ypos = 0;
  if(tanBeta_ == 50) xposi = 100;
  if(tanBeta_ == 50) xpos = 200;
  if(tanBeta_ == 50) ypos = -10;
  
  //TLatex* lumilabel = new TLatex(750 +xposi + 100,767.-154,"#sqrt{s} = 7 TeV, #scale[0.65]{#int}Ldt = 0.98 fb^{-1}");
  //TLatex* lumilabel = new TLatex(450,767.-154+100,"#sqrt{s} = 7 TeV,   Ldt = 1 fb^{-1}");
  TLatex* lumilabel = new TLatex(0,720.,"3.2 fb^{-1}, #sqrt{s} = 7 TeV");
  //TLatex* integral_symbol = new TLatex(577 +xposi + 100,767.-145+100,"#int");

  lumilabel->SetTextSize(0.05);
  //integral_symbol->SetTextSize(0.03);
  lumilabel->Draw("same");
  //integral_symbol->Draw("same");

  //TLatex* cmslabel = new TLatex(10.,767.-154+100,"               ");
  //TLatex* cmslabel = new TLatex(10.,767.-154+100,"CMS Preliminary");
  //cmslabel->SetTextSize(0.05);
  //cmslabel->Draw("same");

  TString text_tanBeta;
  text_tanBeta =  "tan#beta = "+tanb+",  A_{0} = 0,  #mu > 0";
  TLatex* cmssmpars = new TLatex(200,650,text_tanBeta);

  cmssmpars->SetTextSize(0.04);
  cmssmpars->Draw("same");

  TLatex* lep_chargino = new TLatex(250,135,"LEP2 #tilde{#chi}_{1}^{#pm}");
  lep_chargino->SetTextSize(0.03);
  lep_chargino->SetTextFont(42);
  //  lep_chargino->Draw("same");

  TLatex* lep_slepton = new TLatex(26,190,"LEP2 #tilde{#font[12]{l}}^{#pm}");
  lep_slepton->SetTextSize(0.03);
  lep_slepton->SetTextAngle(-83);
  lep_slepton->SetTextFont(42);
  //  lep_slepton->Draw("same");



  //LM points
  TMarker* LM0 = new TMarker(200.,160.,20);
  TMarker* LM1 = new TMarker(60.,250.,20);
  TMarker* LM3 = new TMarker(330.,240.,20);
  TMarker* LM6 = new TMarker(80.,400.,20);
    
  LM0->SetMarkerSize(1.2);
  LM1->SetMarkerSize(1.2);
    
  TLatex* tLM0 = new TLatex(205.,160.," LM0");
  tLM0->SetTextSize(0.035);
    
  TLatex* tLM1 = new TLatex(80.,245.,"LM1");
  tLM1->SetTextSize(0.035);
  
  //TLatex* tLM3 = new TLatex(350.,235.,"LM3 (tan#beta=20)");
  TLatex* tLM3 = new TLatex(350.,235.,"LM3");
  tLM3->SetTextSize(0.035);
  
  TLatex* tLM6 = new TLatex(100.,395.,"LM6");
  tLM6->SetTextSize(0.035);
  
  //  if (tanBeta_ != 50){
  //  LM0->Draw("same");   
  //  tLM0->Draw("same");
  //  LM1->Draw("same");   
  //  tLM1->Draw("same");
  // }

  /*
  if (tanBeta_ == 10){ 
    LM1->Draw("same");
    tLM1->Draw("same");
    LM3->Draw("same");
    tLM3->Draw("same");
    LM6->Draw("same");
    tLM6->Draw("same");
  }
  */



  //stau=LSP contour
  stau->Draw("fsame");
  NoEWSB->Draw("fsame");
 
  //legends
  legexp->Draw();
  legst->Draw();
  //legNoEWSB->Draw();
  //myleg->Draw();

  hist->Draw("sameaxis");
  cvsSys->RedrawAxis();
  cvsSys->Update();
  //cvsSys->Write();
  
  // plots are made here
  //if( plotLO_ ){
  //  cvsSys->SaveAs("ExclusionLimit_tanb"+tanb+"_LO.pdf");
  //  cvsSys->SaveAs("ExclusionLimit_tanb"+tanb+"_LO.png");
  //}else{
  //  cvsSys->SaveAs("ExclusionLimit_tanb"+tanb+".eps");
  //  cvsSys->SaveAs("ExclusionLimit_tanb"+tanb+".pdf");
  //  cvsSys->SaveAs("ExclusionLimit_tanb"+tanb+".png");
  //}                 
  
  return cvsSys;
  //output->Write();
  //output->Close();
  //delete output; 
  
}


void setPlottingStyle(TH1F& hsig){
  
  hsig.SetStats(kFALSE);
  
  hsig.SetAxisRange(80,750,"Y");
  hsig.SetAxisRange(0,520,"X");
  hsig.SetAxisRange(200,520,"X");

  hsig.GetXaxis()->SetTitle("m_{0} (GeV)");
  hsig.GetYaxis()->SetTitle("m_{1/2} (GeV)");
  hsig.GetYaxis()->SetTitleOffset(0.8);
  hsig.GetYaxis()->SetTitleSize(0.06);
  hsig.GetYaxis()->SetLabelSize(0.06);
  hsig.GetXaxis()->SetTitleOffset(0.9);
  hsig.GetXaxis()->SetTitleSize(0.06);
  hsig.GetXaxis()->SetLabelSize(0.06);

  hsig.SetLineWidth(1);  
  hsig.SetLineColor(kBlue);  
  
}


TGraph* set_sneutrino_d0_1(Int_t tanBeta){
  double sn_m0[14]= {0,  0, 48, 55, 80, 90,100,105,109,105,100, 72, 55,0};
  double sn_m12[14]={0,140,210,220,237,241,242,241,230,220,210,170,150,0};

  TGraph* sn_d0_gr = new TGraph(14,sn_m0,sn_m12);

  sn_d0_gr->SetFillColor(kGreen+3);
  sn_d0_gr->SetFillStyle(1001);

  return sn_d0_gr;
}

TGraph* set_sneutrino_d0_2(Int_t tanBeta){
  double sn_m0[9]= {0, 45, 75,115,130,150,163,185,0};
  double sn_m12[9]={0,140,170,213,202,183,168,140,0};

  TGraph* sn_d0_gr_2 = new TGraph(9,sn_m0,sn_m12);

  sn_d0_gr_2->SetFillColor(kGreen+3);
  sn_d0_gr_2->SetFillStyle(1001);

  return sn_d0_gr_2;
}

TGraph* set_lep_ch(Int_t tanBeta){
  //if(tanBeta == 10) 
  return set_lep_ch_tanBeta10();
  //if(tanBeta == 40) return set_lep_ch_tanBeta40();
}

TGraph* set_lep_ch_tanBeta10(){

double ch_m0[] ={0  ,100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000,2000,0.};  
double ch_m12[]={162,162, 161, 160, 160, 159, 158, 157, 156, 155, 154, 154, 153, 152, 150, 149, 148,  147 , 146 , 146 , 146 , 147 , 148 , 149 , 151 , 154 ,159, 0., 0.}; 

  TGraph* ch_gr = new TGraph(29,ch_m0,ch_m12);

  ch_gr->SetFillColor(3);
  ch_gr->SetLineColor(3);
  //  ch_gr->SetLineWidth(3);
  ch_gr->SetFillStyle(1001);

  return ch_gr;

}

TGraph* set_lep_ch_tanBeta40(){

  double ch_m0[] = {0,   240, 400, 500.,700.,800.,1000,1200., 1300., 1400., 1450., 1500., 1550., 1600., 1650., 1700., 1750., 1800., 2150., 2150.,0};
  double ch_m12[]= {140, 140, 139, 138, 137, 136, 136,  137,   138,  139,   140,   142,   143,   145,   147,   150,   154,   158,  134,   0,   0};

  TGraph* ch_gr = new TGraph(21,ch_m0,ch_m12);

  ch_gr->SetFillColor(3);
  ch_gr->SetLineColor(3);
  ch_gr->SetFillStyle(1001);

  return ch_gr;

}

TGraph* my_observed(Int_t tanBeta){

// taken from Frederic Ronga's calculation

    double st_m0_tanBeta10[] =  {0,   10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110 , 130,  147, 0, 170, 190, 210, 230, 250, 270, 290, 310, 330, 350, 370, 390, 410, 430, 450, 470, 490, 510, 530, 550, 570, 590, 610, 630, 6500};
    double st_m12_tanBeta10[] = {213,220,240,275,312,351,393,435,476,518, 559, 600., 682., 750.,750, 842., 921., 999., 1076, 1152, 1228, 1304, 1378, 1453, 1527, 1600, 1673, 1746, 1818, 1890, 1962, 2034, 2105, 2175, 2246, 2316, 2386, 2456, 2526, 2595}; 



    TGraph* st_gr_tanBeta10 = new TGraph(15,st_m0_tanBeta10,st_m12_tanBeta10);
    
    st_gr_tanBeta10->SetFillColor(40);
    st_gr_tanBeta10->SetFillStyle(1001);


    //if(tanBeta == 10)
    return st_gr_tanBeta10;
    //if(tanBeta == 40)return st_gr_tanBeta40;
}


TGraph* set_tev_stau(Int_t tanBeta){

// taken from Frederic Ronga's calculation

    double st_m0_tanBeta10[] =  {0,   10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110 , 130,  147, 0, 170, 190, 210, 230, 250, 270, 290, 310, 330, 350, 370, 390, 410, 430, 450, 470, 490, 510, 530, 550, 570, 590, 610, 630, 6500};
    double st_m12_tanBeta10[] = {213,220,240,275,312,351,393,435,476,518, 559, 600., 682., 750.,750, 842., 921., 999., 1076, 1152, 1228, 1304, 1378, 1453, 1527, 1600, 1673, 1746, 1818, 1890, 1962, 2034, 2105, 2175, 2246, 2316, 2386, 2456, 2526, 2595}; 



    TGraph* st_gr_tanBeta10 = new TGraph(15,st_m0_tanBeta10,st_m12_tanBeta10);
    
    st_gr_tanBeta10->SetFillColor(40);
    st_gr_tanBeta10->SetFillStyle(1001);


    //if(tanBeta == 10)
    return st_gr_tanBeta10;
    //if(tanBeta == 40)return st_gr_tanBeta40;
}

TGraph* set_NoEWSB(Int_t tanBeta){

    double st_m0_tanBeta10[]  = { 0, 70, 90,  80,  60, 1000, 1090, 1170, 1240, 1320, 1370, 1440, 1500, 1550, 1610, 1660, 1720, 1780, 1830, 1860, 1920, 1970 , 2000, 2000, 0};
    double st_m12_tanBeta10[] = {10.,10.,20., 30., 40., 50.,  60.,  70.,  80.,  90.,  100., 110., 120., 130., 140., 150., 160., 170., 180., 190., 200., 210., 220,   0,   0}; 

//// Needs to be modified
    double st_m0_tanBeta40[] = { 0, 70, 90,  80,  60, 1000, 1090, 1170, 1240, 1320, 1370, 1440, 1500, 1550, 1610, 1660, 1720, 1780, 1830, 1860, 1920, 1970 , 2000, 2000, 0};
    double st_m12_tanBeta40[]= {10.,10.,20., 30., 40., 50.,  60.,  70.,  80.,  90.,  100., 110., 120., 130., 140., 150., 160., 170., 180., 190., 200., 210., 220,   0,   0};

    TGraph* st_gr_tanBeta10 = new TGraph(25,st_m0_tanBeta10,st_m12_tanBeta10);
    TGraph* st_gr_tanBeta40 = new TGraph(25,st_m0_tanBeta40,st_m12_tanBeta40);

    std::cout << " tan beta in the set_NoEWSB function: " << tanBeta << std::endl;
    
    st_gr_tanBeta40->SetFillColor(40);
    st_gr_tanBeta40->SetFillStyle(1001);
    
    st_gr_tanBeta10->SetFillColor(40);
    st_gr_tanBeta10->SetFillStyle(1001);

    //if(tanBeta == 10) 
    return st_gr_tanBeta10;
    //if(tanBeta == 40) return st_gr_tanBeta40;
}





TGraph* set_lep_sl(Int_t tanBeta){


  //contour from D0 trilepton paper (PLB 680 (2009) 34-43)

  double *sl_m0 = 0;
  double *sl_m12 = 0;
  int n = 0;

  double sl_m0_3[] ={0,  0, 10, 20, 30, 40, 50, 60, 70, 77,88,95};
  double sl_m12_3[]={0,245,242,239,232,222,209,189,165,140,60,0};
  int n_3 = 12;

  double sl_m0_10[]={ 0,  0, 11, 20, 24, 49, 70, 82,88,90};
  double sl_m12_10[]={0,240,237,233,230,200,150,100,50,0};
  int n_10 = 10;

  if (tanBeta==3){
    sl_m0 = sl_m0_3;
    sl_m12 = sl_m12_3;
    n = n_3;
  }
  //CMS PTDR-II
  //* Selectron_R line mass=99, ISASUGRA7.69, A0=0, m_top=175, tan(beta]=10
  if (tanBeta==10 || tanBeta==50){
    sl_m0 = sl_m0_10;
    sl_m12 = sl_m12_10;
    n = n_10;
  }

  TGraph* lep_sl = new TGraph(n,sl_m0,sl_m12);

  lep_sl->SetFillColor(5);
  lep_sl->SetLineColor(5);
  lep_sl->SetFillStyle(1001);
  
  return lep_sl;
}


TGraph* set_tev_sg_cdf(Int_t tanBeta){

  //New CHF from CDF plot in ICHEP2010 talk (E. Halkiadakis)
  double sg_m0[]= {0,  0, 30, 75,150,185,225,310,360,400,430,500,600,600};
  double sg_m12[]={0,162,168,170,160,150,130,120,109,108,100, 96, 95,  0};
  int np=14;

  TGraph* sg_gr = new TGraph(np,sg_m0,sg_m12);

  //  gStyle->SetHatchesLineWidth(3);

  sg_gr->SetFillColor(2);
  sg_gr->SetLineColor(2);
  //  sg_gr->SetLineWidth(3);
  sg_gr->SetFillStyle(1001); 

  return sg_gr;

}

TGraph* set_tev_sg_d0(Int_t tanBeta){

  double sgd_m0[]= {0,  0, 30, 80,150,240,320,400,500,600,600,0};
  double sgd_m12[]={0,167,166,162,156,138,121,109,105,105,  0,0};
  int npd=12;

  TGraph* sgd_gr = new TGraph(npd,sgd_m0,sgd_m12);

  gStyle->SetHatchesLineWidth(3);

  sgd_gr->SetFillColor(kMagenta+3);
  sgd_gr->SetLineColor(kMagenta+3);
  sgd_gr->SetLineWidth(3);
  sgd_gr->SetFillStyle(3335);

  return sgd_gr;

}



//From Sanjay
TF1* constant_squark(int tanBeta,int i){
//---lines of constant gluino/squark.
// Min squark mass from 1st and 2nd generations using fit for tanbeta = 10.

  double coef1[] = {2.67058e+04, 6.39642e+04, 1.16565e+05, 1.95737e+05, 2.86190e+05};
  double coef2[] = {1.98772e-01, 2.11242e-01, 2.17734e-01, 2.39535e-01, 2.39768e-01};
  double coef3[] = {2.67058e+04, 6.39641e+04, 1.16565e+05, 1.95736e+05, 2.86189e+05};
 
  char hname[200];

  sprintf(hname,"lnsq_%i",i);
  TF1* lnsq = new TF1(hname,"sqrt([0]-x*x*[1]+[2])",0,2000);
  lnsq->SetParameter(0,coef1[i-1]);
  lnsq->SetParameter(1,coef2[i-1]);
  lnsq->SetParameter(2,coef3[i-1]);
  lnsq->SetLineWidth(1);
  lnsq->SetLineColor(kGray);

  return lnsq;
}


TF1* constant_gluino(int tanBeta,int i){
//---lines of constant gluino/squark
  char hname[200];
  sprintf(hname,"lngl_%i",i);

  double coef1[] = {201.77, 311.027, 431.582, 553.895, 676.137};
  double coef2[] = {-0.0146608, -0.01677, -0.022244, -0.0271851, -0.0292212};
   
  TF1* lngl = new TF1(hname,"[0]+x*[1]",0,2000);
  lngl->SetParameter(0,coef1[i-1]);
  lngl->SetParameter(1,coef2[i-1]);
  lngl->SetLineWidth(1);
  lngl->SetLineColor(kGray);

  return lngl;
}


TLatex* constant_squark_text(Int_t it,TF1& lnsq,Int_t tanBeta_){
  char legnm[200];
  it--; //For Sanjay's code

  sprintf(legnm,"#font[92]{#tilde{q}(%i)GeV}",500+250*it);
  Double_t place_x = 170;
  if(tanBeta_ == 50)place_x = 290;
  TLatex* t3 = new TLatex(place_x+10*it,lnsq.Eval(-50+place_x+10*it)+5,legnm);
  t3->SetTextSize(0.03);
  t3->SetTextAngle(-15+it*2);
  t3->SetTextColor(kGray+2);

  return t3;
}


TLatex* constant_gluino_text(Int_t it,TF1& lngl){ //, Int_t tanBeta_){
  char legnm[200];
  Int_t tanBeta_ = 10;
  it--; //For Sanjay's code

  sprintf(legnm,"#font[12]{#tilde{g}}#font[92]{(%i)GeV}",500+250*it);

  Double_t place_x = 423;
  Double_t place_y = 18;
  if (tanBeta_ == 10 ) {
    place_x = 825;
    place_y = -15;
  }
  if (tanBeta_ == 50 ) {
    place_x = 543;
    place_y = 13;
  }
  TLatex* t4 = new TLatex(place_x,place_y+lngl.Eval(480),legnm);
  t4->SetTextSize(0.03);
  t4->SetTextAlign(13);
  t4->SetTextColor(kGray+2);

  return t4;
}



TLegend* makeStauLegend(Double_t txtsz,Int_t tanBeta_){
  Double_t ypos_1 = 0.78;
  Double_t ypos_2 = 0.80;
  Double_t xpos_1 = 0.16;
  Double_t xpos_2 = 0.17;
  if(tanBeta_ == 40){
    xpos_1 = 0.17;
    xpos_2 = 0.18;
    ypos_1 = 0.76;
    ypos_2 = 0.78;

  }
    
  TLegend* legst = new TLegend(xpos_1+0.01,ypos_1,xpos_2+0.01,ypos_2);
  //TLegend* legst = new TLegend(xpos_1+0.025,ypos_1,xpos_2+0.025,ypos_2);
  legst->SetHeader("#tilde{#tau} = LSP");
  legst->SetFillStyle(0);
  legst->SetBorderSize(0);
  legst->SetTextSize(0.03);
  legst->SetTextAngle(80);

  return legst;
}

TLegend* makeNoEWSBLegend(Double_t txtsz,Int_t tanBeta_){
  Double_t ypos_1 = 0.10+0.02;
  Double_t ypos_2 = 0.20+0.02;
  Double_t xpos_1 = 0.82;
  Double_t xpos_2 = 0.92;
  if(tanBeta_ == 40){
    xpos_1 = 0.10;
    xpos_2 = 0.20;
    ypos_1 = 0.85;
    ypos_2 = 0.95;

  }

  TLegend* legst = new TLegend(xpos_1,ypos_1,xpos_2,ypos_2);
  legst->SetHeader("NoEWSB");
  legst->SetFillStyle(0);
  legst->SetBorderSize(0);
  legst->SetTextSize(0.03);
  legst->SetTextAngle(20);

  return legst;
}


TLegend* makeExpLegend(TGraph& sg_gr, TGraph& sgd_gr,TGraph& ch_gr,TGraph& sl_gr,TGraph& tev_sn,Double_t txtsz,Int_t tanbeta){

  //TLegend* legexp = new TLegend(0.61,0.65,0.91,0.9,NULL,"brNDC");
  //TLegend* legexp = new TLegend(0.70,0.70,0.91,0.9,NULL,"brNDC");
  TLegend* legexp = new TLegend(0.61,0.65,0.99,0.9,NULL,"brNDC");


  legexp->SetFillColor(0);
  legexp->SetShadowColor(0);
  legexp->SetTextSize(txtsz);
  legexp->SetBorderSize(0);

  sg_gr.SetLineColor(1);
  
  legexp->AddEntry(&sg_gr,"CDF  #tilde{#font[12]{g}}, #tilde{#font[12]{q}}, #scale[0.8]{tan#beta=5, #mu<0}","f"); 

  legexp->AddEntry(&sgd_gr,"D0   #tilde{#font[12]{g}}, #tilde{#font[12]{q}}, #scale[0.8]{tan#beta=3, #mu<0}","f");  

  ch_gr.SetLineColor(1);
  legexp->AddEntry(&ch_gr,"LEP2   #tilde{#chi}_{1}^{#pm}","f");  
  
  sl_gr.SetLineColor(1);
  if(tanbeta != 50) legexp->AddEntry(&sl_gr,"LEP2   #tilde{#font[12]{l}}^{#pm}","f"); 
  if(tanbeta == 3) legexp->AddEntry(&tev_sn,"D0  #chi^{#pm}_{1}, #chi^{0}_{2}","f");  
 

  return legexp;

}

