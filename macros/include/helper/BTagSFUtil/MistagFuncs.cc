
#include "helper/BTagSFUtil/MistagFuncs.h"
#include <cmath>


MistagFuncs::MistagFuncs(){}



TF1* MistagFuncs::GetMistagmean(TString tagger, TString TaggerStrength, float Etamin, float Etamax)
{
  return GetMistag("mean",tagger,TaggerStrength,Etamin,Etamax);
}
TF1* MistagFuncs::GetMistagmin(TString tagger, TString TaggerStrength, float Etamin, float Etamax)
{
  return GetMistag("min",tagger,TaggerStrength,Etamin,Etamax);
}
TF1* MistagFuncs::GetMistagmax(TString tagger, TString TaggerStrength, float Etamin, float Etamax)
{
  return GetMistag("max",tagger,TaggerStrength,Etamin,Etamax);
}

TF1* MistagFuncs::plotmean(TString tagger, TString TaggerStrength, float Etamin, float Etamax, TString opt  , int col , float lineWidth , int lineStyle )
{
  TF1* f = GetMistagmean(tagger,TaggerStrength,Etamin,Etamax);
  if( f != NULL )
  {
    f->SetLineColor(col);
    //f->SetMinimum(0.4);
    //f->SetMaximum(1.6);
    f->SetLineWidth(lineWidth);
    f->SetLineStyle(lineStyle);
    f->Draw(opt);
  }
  //else std::cout << "NULL pointer returned... Function seems not to exist" << std::endl;
  return f;
}

TF1* MistagFuncs::plotmin(TString tagger, TString TaggerStrength, float Etamin, float Etamax, TString opt , int col , float lineWidth , int lineStyle )
{
  TF1* f = GetMistagmin(tagger,TaggerStrength,Etamin,Etamax);
  if( f != NULL )
  {
    f->SetLineColor(col);
    f->SetLineWidth(lineWidth);
    f->SetLineStyle(lineStyle);
    f->Draw(opt);
  }
  //else std::cout << "NULL pointer returned... Function seems not to exist" << std::endl;
  return f;
}
TF1* MistagFuncs::plotmax(TString tagger, TString TaggerStrength, float Etamin, float Etamax, TString opt, int col, float lineWidth , int lineStyle )
{
  TF1* f = GetMistagmax(tagger,TaggerStrength,Etamin,Etamax);
  if( f != NULL )
  {
    f->SetLineColor(col);
    f->SetLineWidth(lineWidth);
    f->SetLineStyle(lineStyle);
    f->Draw(opt);
  }
  //else std::cout << "NULL pointer returned... Function seems not to exist" << std::endl;
  return f;
}

void MistagFuncs::plotmean(TCanvas *yourC, int yourzone, TString tagger, TString TaggerStrength)
{
 TString legTitle = tagger + TaggerStrength;
 //TCanvas *cWork = new TCanvas("cWork", "plots",200,10,700,750);
 yourC->SetFillColor(10);
 yourC->SetFillStyle(4000);
 yourC->SetBorderSize(2);

 yourC->cd(yourzone);
 yourC->cd(yourzone)->SetFillColor(10);
 yourC->cd(yourzone)->SetFillStyle(4000);
 yourC->cd(yourzone)->SetBorderSize(2);
  TF1 *fmean, *fmin, *fmax;
  TF1* f[10];
  TLegend* leg= new TLegend(0.60,0.56,0.89,0.89);
    leg->SetBorderSize(0);
    leg->SetFillColor(kWhite);
    leg->SetTextFont(42);
    leg->SetHeader(legTitle);
  float etamin[10], etamax[10]; 
  int N=1;
  etamin[0] = 0.0; etamax[0] = 2.4;

  if( TaggerStrength == "L" )
  {
    N = 4;
    etamin[1] = 0.0; etamax[1] = 0.5;
    etamin[2] = 0.5; etamax[2] = 1.0;
    etamin[3] = 1.0; etamax[3] = 1.5;
    etamin[4] = 1.5; etamax[4] = 2.4;
  }
  else if( TaggerStrength == "M" )
  {
    N = 3;
    etamin[1] = 0.0; etamax[1] = 0.8;
    etamin[2] = 0.8; etamax[2] = 1.6;
    etamin[3] = 1.6; etamax[3] = 2.4;
  }
  else if( TaggerStrength == "T" )
  {
    N = 1;
    etamin[1] = 0.0; etamax[1] = 2.4;
  }

  //etamin = 0.0; etamax = 2.4;
/*
  fmean = plotmean(tagger,TaggerStrength,etamin[0], etamax[0], "", 1, 2, 1);
  leg->AddEntry(fmean,"Mean(Mistag)","l");
  fmin = plotmin(tagger,TaggerStrength,etamin[0], etamax[0], "same", 1, 2, 2);
  leg->AddEntry(fmin,"Min(Mistag)","l");
  fmax = plotmax(tagger,TaggerStrength,etamin[0], etamax[0], "same", 1, 2, 2);
  leg->AddEntry(fmax,"Max(Mistag)","l");
*/

    f[1] = plotmean(tagger,TaggerStrength,etamin[1], etamax[1], "", 1, 1);
    //TString rangeEta = Form("Mean(Mistag(%1.1f #leq #eta %1.1f))",etamin[1],etamax[1]);
    TString rangeEta = Form("Mistag(%1.1f #leq #eta %1.1f)",etamin[1],etamax[1]);
    leg->AddEntry(f[1],rangeEta,"l");
  for( int i = 2; i <= N; ++i )
  {
    f[i] = plotmean(tagger,TaggerStrength,etamin[i], etamax[i], "same", i, 1);
    //TString rangeEta = Form("Mean(Mistag(%1.1f #leq #eta %1.1f))",etamin[i],etamax[i]);
    TString rangeEta = Form("Mistag(%1.1f #leq #eta %1.1f)",etamin[i],etamax[i]);
    leg->AddEntry(f[i],rangeEta,"l");
  }
  //leg->AddEntry(gg," gluon jets","P");
  leg->Draw();
  //return cWork;
}

TCanvas * MistagFuncs::plotmean(TString tagger, TString TaggerStrength)
{
 TCanvas *cWork = new TCanvas("cWork", "plots",200,10,700,750);
 plotmean(cWork, 0, tagger, TaggerStrength);

 return cWork;
}

TCanvas *MistagFuncs::plotmean(TString selecter)
{
 TCanvas *cWork = NULL; 
 if( selecter == "L" )
 {
   cWork = new TCanvas("cWork", "plots",200,10,700,500);
   cWork->Divide(2,2);
   plotmean(cWork, 1, "JP", selecter);
   plotmean(cWork, 2, "JBP", selecter);
   plotmean(cWork, 3, "CSV", selecter);
   plotmean(cWork, 4, "TCHE", selecter);
 }
 else if( selecter == "M" )
 {
   cWork = new TCanvas("cWork", "plots",200,10,700,750);
   cWork->Divide(2,3);
   plotmean(cWork, 1, "JP", selecter);
   plotmean(cWork, 2, "JBP", selecter);
   plotmean(cWork, 3, "CSV", selecter);
   plotmean(cWork, 4, "TCHE", selecter);
   plotmean(cWork, 5, "TCHP", selecter);
   plotmean(cWork, 6, "SSVHE", selecter);
 }
 else if( selecter == "T" )
 {
   cWork = new TCanvas("cWork", "plots",200,10,700,750);
   cWork->Divide(2,3);
   plotmean(cWork, 1, "JP", selecter);
   plotmean(cWork, 2, "JBP", selecter);
   plotmean(cWork, 3, "CSV", selecter);
   //plotmean(cWork, 4, "TCHE", selecter);
   plotmean(cWork, 5, "TCHP", selecter);
   plotmean(cWork, 6, "SSVHP", selecter);
 }
 else if( selecter == "TCHE" )
 {
   cWork = new TCanvas("cWork", "plots",200,10,700,500);
   cWork->Divide(1,2);
   plotmean(cWork, 1, selecter, "L");
   plotmean(cWork, 2, selecter, "M");
 }
 else if( selecter == "TCHP" )
 {
   cWork = new TCanvas("cWork", "plots",200,10,700,500);
   cWork->Divide(1,2);
   plotmean(cWork, 1, selecter, "M");
   plotmean(cWork, 2, selecter, "T");
 }
 else if( selecter == "SSVHE" )
 {
   cWork = new TCanvas("cWork", "plots",200,10,700,250);
   plotmean(cWork, 0, selecter, "M");
 }
 else if( selecter == "SSVHP" )
 {
   cWork = new TCanvas("cWork", "plots",200,10,700,250);
   plotmean(cWork, 0, selecter, "T");
 }
 else
 {
   cWork = new TCanvas("cWork", "plots",200,10,700,750);
   cWork->Divide(1,3);
   plotmean(cWork, 1, selecter, "L");
   plotmean(cWork, 2, selecter, "M");
   plotmean(cWork, 3, selecter, "T");
 }

 cWork->WaitPrimitive();
 cWork->SaveAs("MistagFunc_"+selecter+".pdf");
 return cWork;
}

TCanvas * MistagFuncs::plotmean()
{
 TCanvas *cWork = new TCanvas("cWork", "plots",200,10,700,750);
 cWork->Divide(3,6, 0.002, 0.002);
 cWork->SetFillColor(10);
 cWork->SetFillStyle(4000);
 cWork->SetBorderSize(1);
 for( int i = 0; i < 3*6; ++i )
 {
   cWork->cd(i+1)->SetFillColor(10);
   cWork->cd(i+1)->SetFillStyle(4000);
   cWork->cd(i+1)->SetBorderSize(1);
 }
   plotmean(cWork, 1, "JP", "L");
   plotmean(cWork, 2, "JP", "M");
   plotmean(cWork, 3, "JP", "T");
   plotmean(cWork, 4, "JBP", "L");
   plotmean(cWork, 5, "JBP", "M");
   plotmean(cWork, 6, "JBP", "T");
   plotmean(cWork, 7, "CSV", "L");
   plotmean(cWork, 8, "CSV", "M");
   plotmean(cWork, 9, "CSV", "T");
   plotmean(cWork, 10, "TCHE", "L");
   plotmean(cWork, 11, "TCHE", "M");
   plotmean(cWork, 14, "TCHP", "M");
   plotmean(cWork, 15, "TCHP", "T");
   plotmean(cWork, 17, "SSVHE", "M");
   plotmean(cWork, 18, "SSVHP", "T");

 return cWork;
}
/*
TCanvas *plotmean(TString tagger, TString TaggerStrength)
{
 TString legTitle = tagger + TaggerStrength;
 TCanvas *cWork = new TCanvas("cWork", "plots",200,10,700,750);
 cWork->SetFillColor(10);
 cWork->SetFillStyle(4000);
 cWork->SetBorderSize(2);
  TF1 *fmean, *fmin, *fmax;
  TF1* f[10];
  TLegend* leg= new TLegend(0.66,0.54,0.96,0.94);
    leg->SetBorderSize(0);
    leg->SetFillColor(kWhite);
    leg->SetTextFont(42);
    leg->SetHeader(legTitle);
  float etamin[10], etamax[10]; 
  int N=1;
  etamin[0] = 0.0; etamax[0] = 2.4;

  if( TaggerStrength == "L" )
  {
    N = 5;
    etamin[1] = 0.0; etamax[1] = 0.5;
    etamin[2] = 0.5; etamax[2] = 1.0;
    etamin[3] = 1.0; etamax[3] = 1.5;
    etamin[4] = 1.5; etamax[4] = 2.4;
  }
  else if( TaggerStrength == "M" )
  {
    N = 4;
    etamin[1] = 0.0; etamax[1] = 0.8;
    etamin[2] = 0.8; etamax[2] = 1.6;
    etamin[3] = 1.6; etamax[3] = 2.4;
  }
  else if( TaggerStrength == "M" )
  {
    N = 1;
  }

  //etamin = 0.0; etamax = 2.4;
  fmean = plotmean(tagger,TaggerStrength,etamin[0], etamax[0], "", 1, 3, 1);
  leg->AddEntry(fmean,"Mean(SF)","l");
  fmin = plotmin(tagger,TaggerStrength,etamin[0], etamax[0], "same", 1, 3, 2);
  leg->AddEntry(fmin,"Min(SF)","l");
  fmax = plotmax(tagger,TaggerStrength,etamin[0], etamax[0], "same", 1, 3, 2);
  leg->AddEntry(fmax,"Max(SF)","l");
  for( int i = 1; i < N; ++i )
  {
    f[i] = plotmean(tagger,TaggerStrength,etamin[i], etamax[i], "same", i+1, 2);
    TString rangeEta = Form("Mean(SF(%1.1f #leq #eta %1.1f))",etamin[i],etamax[i]);
    leg->AddEntry(f[i],rangeEta,"l");
  }
  //leg->AddEntry(gg," gluon jets","P");
  leg->Draw();
  return cWork;
}
*/


TF1* MistagFuncs::GetMistag(TString meanminmax, TString tagger, TString TaggerStrength, Float_t Etamin, Float_t Etamax)
{
  TF1 *tmpMistag = NULL;
  TString Atagger = tagger+TaggerStrength;
  TString sEtamin = Form("%1.1f",Etamin);
  TString sEtamax = Form("%1.1f",Etamax);
  //std::cout << sEtamin << std::endl;
  //std::cout << sEtamax << std::endl;

// Definition of functions from plot33New.C ----------------------

if( Atagger == "CSVL" && sEtamin == "0.0" && sEtamax == "0.5")
{
if( meanminmax == "mean" ) tmpMistag = new TF1("Mistag","242534*(((1+(0.0182863*x))+(4.50105e-05*(x*x)))/(1+(108569*x)))", 20.,670.);
}
if( Atagger == "CSVL" && sEtamin == "0.0" && sEtamax == "2.4")
{
if( meanminmax == "mean" ) tmpMistag = new TF1("Mistag","18168.8*(((1+(0.020356*x))+(2.73475e-05*(x*x)))/(1+(5239.42*x)))", 20.,670.);
}
if( Atagger == "CSVL" && sEtamin == "0.5" && sEtamax == "1.0")
{
if( meanminmax == "mean" ) tmpMistag = new TF1("Mistag","129.938*(((1+(0.0197657*x))+(4.73472e-05*(x*x)))/(1+(55.2415*x)))", 20.,670.);
}
if( Atagger == "CSVL" && sEtamin == "1.0" && sEtamax == "1.5")
{
if( meanminmax == "mean" ) tmpMistag = new TF1("Mistag","592.214*(((1+(0.00671207*x))+(6.46109e-05*(x*x)))/(1+(134.318*x)))", 20.,670.);
}
if( Atagger == "CSVL" && sEtamin == "1.5" && sEtamax == "2.4")
{
if( meanminmax == "mean" ) tmpMistag = new TF1("Mistag","93329*(((1+(0.0219705*x))+(3.76566e-05*(x*x)))/(1+(18245.1*x)))", 20.,670.);
}
if( Atagger == "CSVM" && sEtamin == "0.0" && sEtamax == "0.8")
{
if( meanminmax == "mean" ) tmpMistag = new TF1("Mistag","(0.00967751+(2.54564e-05*x))+(-6.92256e-10*(x*x))", 20.,670.);
}
if( Atagger == "CSVM" && sEtamin == "0.0" && sEtamax == "2.4")
{
if( meanminmax == "mean" ) tmpMistag = new TF1("Mistag","(0.0113428+(5.18983e-05*x))+(-2.59881e-08*(x*x))", 20.,670.);
}
if( Atagger == "CSVM" && sEtamin == "0.8" && sEtamax == "1.6")
{
if( meanminmax == "mean" ) tmpMistag = new TF1("Mistag","(0.00974141+(5.09503e-05*x))+(2.0641e-08*(x*x))", 20.,670.);
}
if( Atagger == "CSVM" && sEtamin == "1.6" && sEtamax == "2.4")
{
if( meanminmax == "mean" ) tmpMistag = new TF1("Mistag","(0.013595+(0.000104538*x))+(-1.36087e-08*(x*x))", 20.,670.);
}
if( Atagger == "CSVT" && sEtamin == "0.0" && sEtamax == "2.4")
{
if( meanminmax == "mean" ) tmpMistag = new TF1("Mistag","0.00315116*(((1+(-0.00769281*x))+(2.58066e-05*(x*x)))+(-2.02149e-08*(x*(x*x))))", 20.,670.);
}
if( Atagger == "JBPL" && sEtamin == "0.0" && sEtamax == "0.5")
{
if( meanminmax == "mean" ) tmpMistag = new TF1("Mistag","(0.0277261+(0.000808207*x))+(-6.44146e-07*(x*x))", 20.,670.);
}
if( Atagger == "JBPL" && sEtamin == "0.0" && sEtamax == "2.4")
{
if( meanminmax == "mean" ) tmpMistag = new TF1("Mistag","(0.0250569+(0.000856242*x))+(-6.98246e-07*(x*x))", 20.,670.);
}
if( Atagger == "JBPL" && sEtamin == "0.5" && sEtamax == "1.0")
{
if( meanminmax == "mean" ) tmpMistag = new TF1("Mistag","(0.0278926+(0.000827697*x))+(-7.01497e-07*(x*x))", 20.,670.);
}
if( Atagger == "JBPL" && sEtamin == "1.0" && sEtamax == "1.5")
{
if( meanminmax == "mean" ) tmpMistag = new TF1("Mistag","(0.0221411+(0.000900444*x))+(-6.52873e-07*(x*x))", 20.,670.);
}
if( Atagger == "JBPL" && sEtamin == "1.5" && sEtamax == "2.4")
{
if( meanminmax == "mean" ) tmpMistag = new TF1("Mistag","(0.0227045+(0.000808122*x))+(-5.67134e-07*(x*x))", 20.,670.);
}
if( Atagger == "JBPM" && sEtamin == "0.0" && sEtamax == "0.8")
{
if( meanminmax == "mean" ) tmpMistag = new TF1("Mistag","(((0.00206106+(0.000105851*x))+(2.691e-08*(x*x)))+(-4.34651e-11*(x*(x*x))))+(-6.73107e-14*(x*(x*(x*x))))", 20.,670.);
}
if( Atagger == "JBPM" && sEtamin == "0.0" && sEtamax == "2.4")
{
if( meanminmax == "mean" ) tmpMistag = new TF1("Mistag","(((0.00205623+(7.46059e-05*x))+(2.14145e-07*(x*x)))+(-3.0738e-10*(x*(x*x))))+(4.42103e-14*(x*(x*(x*x))))", 20.,670.);
}
if( Atagger == "JBPM" && sEtamin == "0.8" && sEtamax == "1.6")
{
if( meanminmax == "mean" ) tmpMistag = new TF1("Mistag","(((0.00318438+(4.40327e-05*x))+(3.46922e-07*(x*x)))+(-3.93396e-10*(x*(x*x))))+(3.94283e-14*(x*(x*(x*x))))", 20.,670.);
}
if( Atagger == "JBPM" && sEtamin == "1.6" && sEtamax == "2.4")
{
if( meanminmax == "mean" ) tmpMistag = new TF1("Mistag","(((0.00209833+(4.27753e-05*x))+(1.96076e-07*(x*x)))+(6.19275e-11*(x*(x*x))))+(-2.63318e-13*(x*(x*(x*x))))", 20.,670.);
}
if( Atagger == "JBPT" && sEtamin == "0.0" && sEtamax == "2.4")
{
if( meanminmax == "mean" ) tmpMistag = new TF1("Mistag","(-3.36681e-05+(1.37292e-05*x))+(1.78479e-08*(x*x))", 20.,670.);
}
if( Atagger == "JPL" && sEtamin == "0.0" && sEtamax == "0.5")
{
if( meanminmax == "mean" ) tmpMistag = new TF1("Mistag","(0.060001+(0.000332202*x))+(-2.36709e-07*(x*x))", 20.,670.);
}
if( Atagger == "JPL" && sEtamin == "0.0" && sEtamax == "2.4")
{
if( meanminmax == "mean" ) tmpMistag = new TF1("Mistag","(0.055398+(0.00044581*x))+(-3.6301e-07*(x*x))", 20.,670.);
}
if( Atagger == "JPL" && sEtamin == "0.5" && sEtamax == "1.0")
{
if( meanminmax == "mean" ) tmpMistag = new TF1("Mistag","(0.0597675+(0.000370979*x))+(-2.94673e-07*(x*x))", 20.,670.);
}
if( Atagger == "JPL" && sEtamin == "1.0" && sEtamax == "1.5")
{
if( meanminmax == "mean" ) tmpMistag = new TF1("Mistag","(0.0483728+(0.000528418*x))+(-3.17825e-07*(x*x))", 20.,670.);
}
if( Atagger == "JPL" && sEtamin == "1.5" && sEtamax == "2.4")
{
if( meanminmax == "mean" ) tmpMistag = new TF1("Mistag","(0.0463159+(0.000546644*x))+(-3.40486e-07*(x*x))", 20.,670.);
}
if( Atagger == "JPM" && sEtamin == "0.0" && sEtamax == "0.8")
{
if( meanminmax == "mean" ) tmpMistag = new TF1("Mistag","(0.00727084+(4.48901e-05*x))+(-4.42894e-09*(x*x))", 20.,670.);
}
if( Atagger == "JPM" && sEtamin == "0.0" && sEtamax == "2.4")
{
if( meanminmax == "mean" ) tmpMistag = new TF1("Mistag","(0.00461602+(6.31828e-05*x))+(-1.22468e-08*(x*x))", 20.,670.);
}
if( Atagger == "JPM" && sEtamin == "0.8" && sEtamax == "1.6")
{
if( meanminmax == "mean" ) tmpMistag = new TF1("Mistag","(0.00389156+(6.35508e-05*x))+(1.54183e-08*(x*x))", 20.,670.);
}
if( Atagger == "JPM" && sEtamin == "1.6" && sEtamax == "2.4")
{
if( meanminmax == "mean" ) tmpMistag = new TF1("Mistag","(0.0032816+(4.18867e-05*x))+(7.44912e-08*(x*x))", 20.,670.);
}
if( Atagger == "JPT" && sEtamin == "0.0" && sEtamax == "2.4")
{
if( meanminmax == "mean" ) tmpMistag = new TF1("Mistag","(0.000379966+(8.30969e-06*x))+(1.10364e-08*(x*x))", 20.,670.);
}
if( Atagger == "SSVHEM" && sEtamin == "0.0" && sEtamax == "0.8")
{
if( meanminmax == "mean" ) tmpMistag = new TF1("Mistag","(((0.000547883+(0.00023023*x))+(-7.31792e-07*(x*x)))+(1.15659e-09*(x*(x*x))))+(-7.00641e-13*(x*(x*(x*x))))", 20.,670.);
}
if( Atagger == "SSVHEM" && sEtamin == "0.0" && sEtamax == "2.4")
{
if( meanminmax == "mean" ) tmpMistag = new TF1("Mistag","(((-0.000420178+(0.00029105*x))+(-8.9398e-07*(x*x)))+(1.35401e-09*(x*(x*x))))+(-7.93826e-13*(x*(x*(x*x))))", 20.,670.);
}
if( Atagger == "SSVHEM" && sEtamin == "0.8" && sEtamax == "1.6")
{
if( meanminmax == "mean" ) tmpMistag = new TF1("Mistag","(((0.000615562+(0.000240254*x))+(-7.00237e-07*(x*x)))+(1.2566e-09*(x*(x*x))))+(-8.59011e-13*(x*(x*(x*x))))", 20.,670.);
}
if( Atagger == "SSVHEM" && sEtamin == "1.6" && sEtamax == "2.4")
{
if( meanminmax == "mean" ) tmpMistag = new TF1("Mistag","(((0.000372388+(0.000309735*x))+(-4.35952e-07*(x*x)))+(3.63763e-10*(x*(x*x))))+(-2.11993e-13*(x*(x*(x*x))))", 20.,670.);
}
if( Atagger == "SSVHPT" && sEtamin == "0.0" && sEtamax == "2.4")
{
if( meanminmax == "mean" ) tmpMistag = new TF1("Mistag","(-2.9605e-05+(2.35624e-05*x))+(-1.77552e-08*(x*x))", 20.,670.);
}
if( Atagger == "TCHEL" && sEtamin == "0.0" && sEtamax == "0.5")
{
if( meanminmax == "mean" ) tmpMistag = new TF1("Mistag","(((-0.0235318+(0.00268868*x))+(-6.47688e-06*(x*x)))+(7.92087e-09*(x*(x*x))))+(-4.06519e-12*(x*(x*(x*x))))", 20.,670.);
}
if( Atagger == "TCHEL" && sEtamin == "0.0" && sEtamax == "2.4")
{
if( meanminmax == "mean" ) tmpMistag = new TF1("Mistag","(((-0.0276197+(0.00291907*x))+(-7.51594e-06*(x*x)))+(9.82128e-09*(x*(x*x))))+(-5.33759e-12*(x*(x*(x*x))))", 20.,670.);
}
if( Atagger == "TCHEL" && sEtamin == "0.5" && sEtamax == "1.0")
{
if( meanminmax == "mean" ) tmpMistag = new TF1("Mistag","(((-0.0257274+(0.00289337*x))+(-7.48879e-06*(x*x)))+(9.84928e-09*(x*(x*x))))+(-5.40844e-12*(x*(x*(x*x))))", 20.,670.);
}
if( Atagger == "TCHEL" && sEtamin == "1.0" && sEtamax == "1.5")
{
if( meanminmax == "mean" ) tmpMistag = new TF1("Mistag","(((-0.0310046+(0.00307803*x))+(-7.94145e-06*(x*x)))+(1.06889e-08*(x*(x*x))))+(-6.08971e-12*(x*(x*(x*x))))", 20.,670.);
}
if( Atagger == "TCHEL" && sEtamin == "1.5" && sEtamax == "2.4")
{
if( meanminmax == "mean" ) tmpMistag = new TF1("Mistag","(((-0.0274561+(0.00301096*x))+(-8.89588e-06*(x*x)))+(1.40142e-08*(x*(x*x))))+(-8.95723e-12*(x*(x*(x*x))))", 20.,670.);
}
if( Atagger == "TCHEM" && sEtamin == "0.0" && sEtamax == "0.8")
{
if( meanminmax == "mean" ) tmpMistag = new TF1("Mistag","(0.000919586+(0.00026266*x))+(-1.75723e-07*(x*x))", 20.,670.);
}
if( Atagger == "TCHEM" && sEtamin == "0.0" && sEtamax == "2.4")
{
if( meanminmax == "mean" ) tmpMistag = new TF1("Mistag","(-0.00256163+(0.000332759*x))+(-2.39887e-07*(x*x))", 20.,670.);
}
if( Atagger == "TCHEM" && sEtamin == "0.8" && sEtamax == "1.6")
{
if( meanminmax == "mean" ) tmpMistag = new TF1("Mistag","(-0.00364137+(0.000350371*x))+(-1.89967e-07*(x*x))", 20.,670.);
}
if( Atagger == "TCHEM" && sEtamin == "1.6" && sEtamax == "2.4")
{
if( meanminmax == "mean" ) tmpMistag = new TF1("Mistag","(-0.00483904+(0.000367751*x))+(-1.36152e-07*(x*x))", 20.,670.);
}
if( Atagger == "TCHPM" && sEtamin == "0.0" && sEtamax == "0.8")
{
if( meanminmax == "mean" ) tmpMistag = new TF1("Mistag","(((-0.00464673+(0.000247485*x))+(9.13236e-07*(x*x)))+(-2.49994e-09*(x*(x*x))))+(1.65678e-12*(x*(x*(x*x))))", 20.,670.);
}
if( Atagger == "TCHPM" && sEtamin == "0.0" && sEtamax == "2.4")
{
if( meanminmax == "mean" ) tmpMistag = new TF1("Mistag","(((-0.00609327+(0.00029725*x))+(9.46617e-07*(x*x)))+(-2.71065e-09*(x*(x*x))))+(1.77968e-12*(x*(x*(x*x))))", 20.,670.);
}
if( Atagger == "TCHPM" && sEtamin == "0.8" && sEtamax == "1.6")
{
if( meanminmax == "mean" ) tmpMistag = new TF1("Mistag","(((-0.0060878+(0.000297422*x))+(1.13369e-06*(x*x)))+(-2.84945e-09*(x*(x*x))))+(1.64721e-12*(x*(x*(x*x))))", 20.,670.);
}
if( Atagger == "TCHPM" && sEtamin == "1.6" && sEtamax == "2.4")
{
if( meanminmax == "mean" ) tmpMistag = new TF1("Mistag","(((-0.00836219+(0.000391889*x))+(2.78156e-07*(x*x)))+(-6.14017e-10*(x*(x*x))))+(-1.30592e-13*(x*(x*(x*x))))", 20.,670.);
}
if( Atagger == "TCHPT" && sEtamin == "0.0" && sEtamax == "2.4")
{
if( meanminmax == "mean" ) tmpMistag = new TF1("Mistag","(-0.00101+(4.70405e-05*x))+(8.3338e-09*(x*x))", 20.,670.);
}

// End of definition of functions from plot33New.C ---------------
  if( tmpMistag == NULL ) std::cout << "NULL pointer returned... Function seems not to exist" << std::endl;

  return tmpMistag;
}





float MistagFuncs::GetMistag_fast(const TString& tagger, const TString& TaggerStrength, Float_t pt, Float_t eta ) {

  float x = pt;
  if( x < 20. ) x = 20.;
  if( x > 670. ) x = 670.;

  float tmpMistag = 0.;

  TString Atagger = tagger+TaggerStrength;

// Definition of functions from plot33New.C ----------------------

if( Atagger == "CSVL" ) {
  if( fabs(eta) >= 0.0 && fabs(eta) < 0.5)
  {
   tmpMistag = 242534*(((1+(0.0182863*x))+(4.50105e-05*(x*x)))/(1+(108569*x)));
  }
  else if( fabs(eta) >= 0.5 && fabs(eta) < 1.0)
  {
   tmpMistag = 129.938*(((1+(0.0197657*x))+(4.73472e-05*(x*x)))/(1+(55.2415*x)));
  }
  else if( fabs(eta) >= 1.0 && fabs(eta) < 1.5)
  {
   tmpMistag = 592.214*(((1+(0.00671207*x))+(6.46109e-05*(x*x)))/(1+(134.318*x)));
  }
  else if( fabs(eta) >= 1.5 && fabs(eta) < 2.4)
  {
   tmpMistag = 93329*(((1+(0.0219705*x))+(3.76566e-05*(x*x)))/(1+(18245.1*x)));
  }
} 


else if( Atagger == "CSVM" ) {
  if( fabs(eta) >= 0.0 && fabs(eta) < 0.8)
  {
   tmpMistag = (0.00967751+(2.54564e-05*x))+(-6.92256e-10*(x*x));
  }
  else if( fabs(eta) >= 0.8 && fabs(eta) < 1.6)
  {
   tmpMistag = (0.00974141+(5.09503e-05*x))+(2.0641e-08*(x*x));
  }
  else if( fabs(eta) >= 1.6 && fabs(eta) < 2.4)
  {
   tmpMistag = (0.013595+(0.000104538*x))+(-1.36087e-08*(x*x));
  }
}


else if( Atagger == "CSVT" ) {
 tmpMistag = 0.00315116*(((1+(-0.00769281*x))+(2.58066e-05*(x*x)))+(-2.02149e-08*(x*(x*x))));
}


else if( Atagger == "TCHEL" ) {
  if( fabs(eta) >= 0.0 && fabs(eta) < 0.5)
  {
   tmpMistag = (((-0.0235318+(0.00268868*x))+(-6.47688e-06*(x*x)))+(7.92087e-09*(x*(x*x))))+(-4.06519e-12*(x*(x*(x*x))));
  }
  else if( fabs(eta) >= 0.5 && fabs(eta) < 1.0)
  {
   tmpMistag = (((-0.0257274+(0.00289337*x))+(-7.48879e-06*(x*x)))+(9.84928e-09*(x*(x*x))))+(-5.40844e-12*(x*(x*(x*x))));
  }
  else if( fabs(eta) >= 1.0 && fabs(eta) < 1.5)
  {
   tmpMistag = (((-0.0310046+(0.00307803*x))+(-7.94145e-06*(x*x)))+(1.06889e-08*(x*(x*x))))+(-6.08971e-12*(x*(x*(x*x))));
  }
  else if( fabs(eta) >= 1.5 && fabs(eta) < 2.4)
  {
   tmpMistag = (((-0.0274561+(0.00301096*x))+(-8.89588e-06*(x*x)))+(1.40142e-08*(x*(x*x))))+(-8.95723e-12*(x*(x*(x*x))));
  }
}


else if( Atagger == "TCHEM" ) {
  if( fabs(eta) >= 0.0 && fabs(eta) < 0.8)
  {
   tmpMistag = (0.000919586+(0.00026266*x))+(-1.75723e-07*(x*x));
  }
  else if( fabs(eta) >= 0.8 && fabs(eta) < 1.6)
  {
   tmpMistag = (-0.00364137+(0.000350371*x))+(-1.89967e-07*(x*x));
  }
  else if( fabs(eta) >= 1.6 && fabs(eta) < 2.4)
  {
   tmpMistag = (-0.00483904+(0.000367751*x))+(-1.36152e-07*(x*x));
  }
}


else if( Atagger == "JBPL" ) {
  if( fabs(eta) >= 0.0 && fabs(eta) < 0.5 ) 
  {
   tmpMistag = (0.0277261+(0.000808207*x))+(-6.44146e-07*(x*x));
  }
  else if( fabs(eta) >= 0.5 && fabs(eta) < 1.0)
  {
   tmpMistag = (0.0278926+(0.000827697*x))+(-7.01497e-07*(x*x));
  }
  else if( fabs(eta) >= 1.0 && fabs(eta) < 1.5)
  {
   tmpMistag = (0.0221411+(0.000900444*x))+(-6.52873e-07*(x*x));
  }
  else if( fabs(eta) >= 1.5 && fabs(eta) < 2.4)
  {
   tmpMistag = (0.0227045+(0.000808122*x))+(-5.67134e-07*(x*x));
  }

}


else if( Atagger == "JBPM" )
{
  if( fabs(eta) >= 0.0 && fabs(eta) < 0.8)
  {
   tmpMistag = (((0.00206106+(0.000105851*x))+(2.691e-08*(x*x)))+(-4.34651e-11*(x*(x*x))))+(-6.73107e-14*(x*(x*(x*x))));
  }
  else if( fabs(eta) >= 0.8 && fabs(eta) < 1.6)
  {
   tmpMistag = (((0.00318438+(4.40327e-05*x))+(3.46922e-07*(x*x)))+(-3.93396e-10*(x*(x*x))))+(3.94283e-14*(x*(x*(x*x))));
  }
  else if( fabs(eta) >= 1.6 && fabs(eta) < 2.4)
  {
   tmpMistag = (((0.00209833+(4.27753e-05*x))+(1.96076e-07*(x*x)))+(6.19275e-11*(x*(x*x))))+(-2.63318e-13*(x*(x*(x*x))));
  }
}


else if( Atagger == "JBPT" ) {
 tmpMistag = (-3.36681e-05+(1.37292e-05*x))+(1.78479e-08*(x*x));
}


else if( Atagger == "JPL" ) {
  if( fabs(eta) >= 0.0 && fabs(eta) < 0.5)
  {
   tmpMistag = (0.060001+(0.000332202*x))+(-2.36709e-07*(x*x));
  }
  else if( fabs(eta) >= 0.5 && fabs(eta) < 1.0)
  {
   tmpMistag = (0.0597675+(0.000370979*x))+(-2.94673e-07*(x*x));
  }
  else if( fabs(eta) >= 1.0 && fabs(eta) < 1.5)
  {
   tmpMistag = (0.0483728+(0.000528418*x))+(-3.17825e-07*(x*x));
  }
  else if( fabs(eta) >= 1.5 && fabs(eta) < 2.4)
  {
   tmpMistag = (0.0463159+(0.000546644*x))+(-3.40486e-07*(x*x));
  }
}


else if( Atagger == "JPM" ) {
  if( fabs(eta) >= 0.0 && fabs(eta) < 0.8)
  {
   tmpMistag = (0.00727084+(4.48901e-05*x))+(-4.42894e-09*(x*x));
  }
  else if( fabs(eta) >= 0.8 && fabs(eta) < 1.6)
  {
   tmpMistag = (0.00389156+(6.35508e-05*x))+(1.54183e-08*(x*x));
  }
  else if( fabs(eta) >= 1.6 && fabs(eta) < 2.4)
  {
   tmpMistag = (0.0032816+(4.18867e-05*x))+(7.44912e-08*(x*x));
  }
}


else if( Atagger == "JPT" ) {
 tmpMistag = (0.000379966+(8.30969e-06*x))+(1.10364e-08*(x*x));
}


else if( Atagger == "SSVHEM" ) {
  if( fabs(eta) >= 0.0 && fabs(eta) < 0.8)
  {
   tmpMistag = (((0.000547883+(0.00023023*x))+(-7.31792e-07*(x*x)))+(1.15659e-09*(x*(x*x))))+(-7.00641e-13*(x*(x*(x*x))));
  }
  else if( fabs(eta) >= 0.8 && fabs(eta) < 1.6)
  {
   tmpMistag = (((0.000615562+(0.000240254*x))+(-7.00237e-07*(x*x)))+(1.2566e-09*(x*(x*x))))+(-8.59011e-13*(x*(x*(x*x))));
  }
  else if( fabs(eta) >= 1.6 && fabs(eta) < 2.4)
  {
   tmpMistag = (((0.000372388+(0.000309735*x))+(-4.35952e-07*(x*x)))+(3.63763e-10*(x*(x*x))))+(-2.11993e-13*(x*(x*(x*x))));
  }
}


else if( Atagger == "SSVHPT" ) {
 tmpMistag = (-2.9605e-05+(2.35624e-05*x))+(-1.77552e-08*(x*x));
}



else if( Atagger == "TCHPM" ) {
  if( fabs(eta) >= 0.0 && fabs(eta) < 0.8)
  {
   tmpMistag = (((-0.00464673+(0.000247485*x))+(9.13236e-07*(x*x)))+(-2.49994e-09*(x*(x*x))))+(1.65678e-12*(x*(x*(x*x))));
  }
  else if( fabs(eta) >= 0.8 && fabs(eta) < 1.6)
  {
   tmpMistag = (((-0.0060878+(0.000297422*x))+(1.13369e-06*(x*x)))+(-2.84945e-09*(x*(x*x))))+(1.64721e-12*(x*(x*(x*x))));
  }
  else if( fabs(eta) >= 1.6 && fabs(eta) < 2.4)
  {
   tmpMistag = (((-0.00836219+(0.000391889*x))+(2.78156e-07*(x*x)))+(-6.14017e-10*(x*(x*x))))+(-1.30592e-13*(x*(x*(x*x))));
  }
}


else if( Atagger == "TCHPT" ) {
 tmpMistag = (-0.00101+(4.70405e-05*x))+(8.3338e-09*(x*x));
}


// End of definition of functions from plot33New.C ---------------

  return tmpMistag;
}

