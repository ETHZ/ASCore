#include <iostream>
#include <fstream>
#include "helper/BTagSFUtil/BTagSFUtil.h"
#include "TMath.h"




BTagSFUtil::BTagSFUtil( int seed ) {

  rand_ = new TRandom3(seed);
  setSFFileName("");

}


BTagSFUtil::BTagSFUtil( const std::string& btagAlgo, int seed ) {

  rand_ = new TRandom3(seed);
  setSFFileName("");

  btagAlgo_ = btagAlgo;

  init_loose_ = false;
  init_medium_ = false;
  init_tight_ = false;

  // this function used if no SF has to be applied
  f1_one_ = new TF1("f1_one", "1.", 0., 10000.);

}


BTagSFUtil::~BTagSFUtil() {

  delete rand_;
  delete f1_one_;

}



void BTagSFUtil::init( const std::string& wp ) {

  std::cout << "[BTagSFUtil::init] Initializing working point: '" << wp << "'." << std::endl;

  InitSFLight(wp);
  InitMistag(wp);
  InitSFb(wp);

  if( wp=="L" ) init_loose_ = true;
  else if( wp=="M" ) init_medium_ = true;
  else if( wp=="T" ) init_tight_ = true;

}



void  BTagSFUtil::InitSFLight( const std::string& wp )
{


  TString Atagger = btagAlgo_+wp;

  float ptMin = 20.;
  float ptMax = 670.;
// Definition of functions from plot33New.C ----------------------

if( Atagger == "CSVL" ) {

  etaBins_SFLight_loose_.push_back(0.);
  etaBins_SFLight_loose_.push_back(0.5);
  etaBins_SFLight_loose_.push_back(1.);
  etaBins_SFLight_loose_.push_back(1.5);
  etaBins_SFLight_loose_.push_back(2.4);

  TF1* sfFunct_0 = new TF1("sfFunct_loose_0", "((1.07536+(0.000175506*x))+(-8.63317e-07*(x*x)))+(3.27516e-10*(x*(x*x)))", ptMin, ptMax);
  TF1* sfFunct_min_0 = new TF1("sfFunct_loose_min_0", "((0.994425+(-8.66392e-05*x))+(-3.03813e-08*(x*x)))+(-3.52151e-10*(x*(x*x)))", ptMin, ptMax);
  TF1* sfFunct_max_0 = new TF1("sfFunct_loose_max_0", "((1.15628+(0.000437668*x))+(-1.69625e-06*(x*x)))+(1.00718e-09*(x*(x*x)))", ptMin, ptMax);

  TF1* sfFunct_1 = new TF1("sfFunct_loose_1", "((1.07846+(0.00032458*x))+(-1.30258e-06*(x*x)))+(8.50608e-10*(x*(x*x)))", ptMin, ptMax);
  TF1* sfFunct_min_1 = new TF1("sfFunct_loose_min_1", "((0.998088+(6.94916e-05*x))+(-4.82731e-07*(x*x)))+(1.63506e-10*(x*(x*x)))", ptMin, ptMax);
  TF1* sfFunct_max_1 = new TF1("sfFunct_loose_max_1", "((1.15882+(0.000579711*x))+(-2.12243e-06*(x*x)))+(1.53771e-09*(x*(x*x)))", ptMin, ptMax);

  TF1* sfFunct_2 = new TF1("sfFunct_loose_2", "((1.08294+(0.000474818*x))+(-1.43857e-06*(x*x)))+(1.13308e-09*(x*(x*x)))", ptMin, ptMax);
  TF1* sfFunct_min_2 = new TF1("sfFunct_loose_min_2", "((1.00294+(0.000289844*x))+(-7.9845e-07*(x*x)))+(5.38525e-10*(x*(x*x)))", ptMin, ptMax);
  TF1* sfFunct_max_2 = new TF1("sfFunct_loose_max_2", "((1.16292+(0.000659848*x))+(-2.07868e-06*(x*x)))+(1.72763e-09*(x*(x*x)))", ptMin, ptMax);

  TF1* sfFunct_3 = new TF1("sfFunct_loose_3", "((1.0617+(0.000173654*x))+(-5.29009e-07*(x*x)))+(5.55931e-10*(x*(x*x)))", ptMin, ptMax);
  TF1* sfFunct_min_3 = new TF1("sfFunct_loose_min_3", "((0.979816+(0.000138797*x))+(-3.14503e-07*(x*x)))+(2.38124e-10*(x*(x*x)))", ptMin, ptMax);
  TF1* sfFunct_max_3 = new TF1("sfFunct_loose_max_3", "((1.14357+(0.00020854*x))+(-7.43519e-07*(x*x)))+(8.73742e-10*(x*(x*x)))", ptMin, ptMax);

  sfLightFunct_loose_.push_back(sfFunct_0);
  sfLightFunct_loose_min_.push_back(sfFunct_min_0);
  sfLightFunct_loose_max_.push_back(sfFunct_max_0);

  sfLightFunct_loose_.push_back(sfFunct_1);
  sfLightFunct_loose_min_.push_back(sfFunct_min_1);
  sfLightFunct_loose_max_.push_back(sfFunct_max_1);

  sfLightFunct_loose_.push_back(sfFunct_2);
  sfLightFunct_loose_min_.push_back(sfFunct_min_2);
  sfLightFunct_loose_max_.push_back(sfFunct_max_2);

  sfLightFunct_loose_.push_back(sfFunct_3);
  sfLightFunct_loose_min_.push_back(sfFunct_min_3);
  sfLightFunct_loose_max_.push_back(sfFunct_max_3);


}


else if( Atagger == "CSVM" ) {

  etaBins_SFLight_medium_.push_back(0.);
  etaBins_SFLight_medium_.push_back(0.8);
  etaBins_SFLight_medium_.push_back(1.6);
  etaBins_SFLight_medium_.push_back(2.4);

  TF1* sfFunct_0 = new TF1("sfFunct_medium", "((1.06182+(0.000617034*x))+(-1.5732e-06*(x*x)))+(3.02909e-10*(x*(x*x)))", ptMin, ptMax);
  TF1* sfFunct_min_0 = new TF1("sfFunct_medium_min_0", "((0.972455+(7.51396e-06*x))+(4.91857e-07*(x*x)))+(-1.47661e-09*(x*(x*x)))", ptMin, ptMax);
  TF1* sfFunct_max_0 = new TF1("sfFunct_medium_max_0", "((1.15116+(0.00122657*x))+(-3.63826e-06*(x*x)))+(2.08242e-09*(x*(x*x)))", ptMin, ptMax);

  TF1* sfFunct_1 = new TF1("sfFunct_medium_1", "((1.111+(-9.64191e-06*x))+(1.80811e-07*(x*x)))+(-5.44868e-10*(x*(x*x)))", ptMin, ptMax);
  TF1* sfFunct_min_1 = new TF1("sfFunct_medium_min_1", "((1.02055+(-0.000378856*x))+(1.49029e-06*(x*x)))+(-1.74966e-09*(x*(x*x)))", ptMin, ptMax);
  TF1* sfFunct_max_1 = new TF1("sfFunct_medium_max_1", "((1.20146+(0.000359543*x))+(-1.12866e-06*(x*x)))+(6.59918e-10*(x*(x*x)))", ptMin, ptMax);

  TF1* sfFunct_2 = new TF1("sfFunct_medium_2", "((1.08498+(-0.000701422*x))+(3.43612e-06*(x*x)))+(-4.11794e-09*(x*(x*x)))", ptMin, ptMax);
  TF1* sfFunct_min_2 = new TF1("sfFunct_medium_min_2", "((0.983476+(-0.000607242*x))+(3.17997e-06*(x*x)))+(-4.01242e-09*(x*(x*x)))", ptMin, ptMax);
  TF1* sfFunct_max_2 = new TF1("sfFunct_medium_max_2", "((1.18654+(-0.000795808*x))+(3.69226e-06*(x*x)))+(-4.22347e-09*(x*(x*x)))", ptMin, ptMax);

  sfLightFunct_medium_.push_back(sfFunct_0);
  sfLightFunct_medium_min_.push_back(sfFunct_min_0);
  sfLightFunct_medium_max_.push_back(sfFunct_max_0);

  sfLightFunct_medium_.push_back(sfFunct_1);
  sfLightFunct_medium_min_.push_back(sfFunct_min_1);
  sfLightFunct_medium_max_.push_back(sfFunct_max_1);

  sfLightFunct_medium_.push_back(sfFunct_2);
  sfLightFunct_medium_min_.push_back(sfFunct_min_2);
  sfLightFunct_medium_max_.push_back(sfFunct_max_2);

}


else if( Atagger == "CSVT" )
{

  etaBins_SFLight_tight_.push_back(0.);
  etaBins_SFLight_tight_.push_back(2.4);

  TF1* sfFunct = new TF1("sfFunct_tight", "((0.948463+(0.00288102*x))+(-7.98091e-06*(x*x)))+(5.50157e-09*(x*(x*x)))", ptMin, ptMax);
  TF1* sfFunct_min = new TF1("sfFunct_tight_min", "((0.899715+(0.00102278*x))+(-2.46335e-06*(x*x)))+(9.71143e-10*(x*(x*x)))", ptMin, ptMax);
  TF1* sfFunct_max = new TF1("sfFunct_tight_max", "((0.997077+(0.00473953*x))+(-1.34985e-05*(x*x)))+(1.0032e-08*(x*(x*x)))", ptMin, ptMax);

  sfLightFunct_tight_.push_back(sfFunct);
  sfLightFunct_tight_min_.push_back(sfFunct_min);
  sfLightFunct_tight_max_.push_back(sfFunct_max);

}



else if( Atagger == "TCHEL" )
{

  etaBins_SFLight_loose_.push_back(0.);
  etaBins_SFLight_loose_.push_back(0.5);
  etaBins_SFLight_loose_.push_back(1.);
  etaBins_SFLight_loose_.push_back(1.5);
  etaBins_SFLight_loose_.push_back(2.4);

  TF1* sfFunct_0 = new TF1("sfFunct_loose_0", "(1.13615*((1+(-0.00119852*x))+(1.17888e-05*(x*x))))+(-9.8581e-08*(x*(x*(x/(1+(0.00689317*x))))))", ptMin, ptMax);
  TF1* sfFunct_min_0 = new TF1("sfFunct_loose_min_0", "(1.0369*((1+(-0.000945578*x))+(7.73273e-06*(x*x))))+(-4.47791e-08*(x*(x*(x/(1+(0.00499343*x))))))", ptMin, ptMax);
  TF1* sfFunct_max_0 = new TF1("sfFunct_loose_max_0", "(1.22179*((1+(-0.000946228*x))+(7.37821e-06*(x*x))))+(-4.8451e-08*(x*(x*(x/(1+(0.0047976*x))))))", ptMin, ptMax);

  TF1* sfFunct_1 = new TF1("sfFunct_loose_1", "(1.13277*((1+(-0.00084146*x))+(3.80313e-06*(x*x))))+(-8.75061e-09*(x*(x*(x/(1+(0.00118695*x))))))", ptMin, ptMax);
  TF1* sfFunct_min_1 = new TF1("sfFunct_loose_min_1", "(0.983748*((1+(7.13613e-05*x))+(-1.08648e-05*(x*x))))+(2.96162e-06*(x*(x*(x/(1+(0.282104*x))))))", ptMin, ptMax);
  TF1* sfFunct_max_1 = new TF1("sfFunct_loose_max_1", "(1.22714*((1+(-0.00085562*x))+(3.74425e-06*(x*x))))+(-8.91028e-09*(x*(x*(x/(1+(0.00109346*x))))))", ptMin, ptMax);

  TF1* sfFunct_2 = new TF1("sfFunct_loose_2", "(1.17163*((1+(-0.000828475*x))+(3.0769e-06*(x*x))))+(-4.692e-09*(x*(x*(x/(1+(0.000337759*x))))))", ptMin, ptMax);
  TF1* sfFunct_min_2 = new TF1("sfFunct_loose_min_2", "(1.0698*((1+(-0.000731877*x))+(2.56922e-06*(x*x))))+(-3.0318e-09*(x*(x*(x/(1+(5.04118e-05*x))))))", ptMin, ptMax);
  TF1* sfFunct_max_2 = new TF1("sfFunct_loose_max_2", "(1.27351*((1+(-0.000911891*x))+(3.5465e-06*(x*x))))+(-6.69625e-09*(x*(x*(x/(1+(0.000590847*x))))))", ptMin, ptMax);

  TF1* sfFunct_3 = new TF1("sfFunct_loose_3", "(1.14554*((1+(-0.000128043*x))+(4.10899e-07*(x*x))))+(-2.07565e-10*(x*(x*(x/(1+(-0.00118618*x))))))", ptMin, ptMax);
  TF1* sfFunct_min_3 = new TF1("sfFunct_loose_min_3", "(1.04766*((1+(-6.87499e-05*x))+(2.2454e-07*(x*x))))+(-1.18395e-10*(x*(x*(x/(1+(-0.00128734*x))))))", ptMin, ptMax);
  TF1* sfFunct_max_3 = new TF1("sfFunct_loose_max_3", "(1.24367*((1+(-0.000182494*x))+(5.92637e-07*(x*x))))+(-3.3745e-10*(x*(x*(x/(1+(-0.00107694*x))))))", ptMin, ptMax);

  sfLightFunct_loose_.push_back(sfFunct_0);
  sfLightFunct_loose_min_.push_back(sfFunct_min_0);
  sfLightFunct_loose_max_.push_back(sfFunct_max_0);

  sfLightFunct_loose_.push_back(sfFunct_1);
  sfLightFunct_loose_min_.push_back(sfFunct_min_1);
  sfLightFunct_loose_max_.push_back(sfFunct_max_1);

  sfLightFunct_loose_.push_back(sfFunct_2);
  sfLightFunct_loose_min_.push_back(sfFunct_min_2);
  sfLightFunct_loose_max_.push_back(sfFunct_max_2);

  sfLightFunct_loose_.push_back(sfFunct_3);
  sfLightFunct_loose_min_.push_back(sfFunct_min_3);
  sfLightFunct_loose_max_.push_back(sfFunct_max_3);

}



else if( Atagger == "TCHEM" ) 
{

  etaBins_SFLight_medium_.push_back(0.);
  etaBins_SFLight_medium_.push_back(0.8);
  etaBins_SFLight_medium_.push_back(1.6);
  etaBins_SFLight_medium_.push_back(2.4);

  TF1* sfFunct_0 = new TF1("sfFunct_medium_0", "(1.2875*((1+(-0.000356371*x))+(1.08081e-07*(x*x))))+(-6.89998e-11*(x*(x*(x/(1+(-0.0012139*x))))))", ptMin, ptMax);
  TF1* sfFunct_min_0 = new TF1("sfFunct_medium_min_0", "(1.11418*((1+(-0.000442274*x))+(1.53463e-06*(x*x))))+(-4.93683e-09*(x*(x*(x/(1+(0.00152436*x))))))", ptMin, ptMax);
  TF1* sfFunct_max_0 = new TF1("sfFunct_medium_max_0", "(1.47515*((1+(-0.000484868*x))+(2.36817e-07*(x*x))))+(-2.05073e-11*(x*(x*(x/(1+(-0.00142819*x))))))", ptMin, ptMax);

  TF1* sfFunct_1 = new TF1("sfFunct_medium_1", "(1.24986*((1+(-0.00039734*x))+(5.37486e-07*(x*x))))+(-1.74023e-10*(x*(x*(x/(1+(-0.00112954*x))))))", ptMin, ptMax);
  TF1* sfFunct_min_1 = new TF1("sfFunct_medium_min_1", "(1.08828*((1+(-0.000208737*x))+(1.50487e-07*(x*x))))+(-2.54249e-11*(x*(x*(x/(1+(-0.00141477*x))))))", ptMin, ptMax);
  TF1* sfFunct_max_1 = new TF1("sfFunct_medium_max_1", "(1.41211*((1+(-0.000559603*x))+(9.50754e-07*(x*x))))+(-5.81148e-10*(x*(x*(x/(1+(-0.000787359*x))))))", ptMin, ptMax);

  TF1* sfFunct_2 = new TF1("sfFunct_medium_2", "(1.10763*((1+(-0.000105805*x))+(7.11718e-07*(x*x))))+(-5.3001e-10*(x*(x*(x/(1+(-0.000821215*x))))))", ptMin, ptMax);
  TF1* sfFunct_min_2 = new TF1("sfFunct_medium_min_2", "(0.958079*((1+(0.000327804*x))+(-4.09511e-07*(x*x))))+(-1.95933e-11*(x*(x*(x/(1+(-0.00143323*x))))))", ptMin, ptMax);
  TF1* sfFunct_max_2 = new TF1("sfFunct_medium_max_2", "(1.26236*((1+(-0.000524055*x))+(2.08863e-06*(x*x))))+(-2.29473e-09*(x*(x*(x/(1+(-0.000276268*x))))))", ptMin, ptMax);

  sfLightFunct_medium_.push_back(sfFunct_0);
  sfLightFunct_medium_min_.push_back(sfFunct_min_0);
  sfLightFunct_medium_max_.push_back(sfFunct_max_0);

  sfLightFunct_medium_.push_back(sfFunct_1);
  sfLightFunct_medium_min_.push_back(sfFunct_min_1);
  sfLightFunct_medium_max_.push_back(sfFunct_max_1);

  sfLightFunct_medium_.push_back(sfFunct_2);
  sfLightFunct_medium_min_.push_back(sfFunct_min_2);
  sfLightFunct_medium_max_.push_back(sfFunct_max_2);

}



else if( Atagger == "JPL" ) 
{

  etaBins_SFLight_loose_.push_back(0.);
  etaBins_SFLight_loose_.push_back(0.5);
  etaBins_SFLight_loose_.push_back(1.);
  etaBins_SFLight_loose_.push_back(1.5);
  etaBins_SFLight_loose_.push_back(2.4);

  TF1* sfFunct_0 = new TF1("sfFunct_loose_0", "((1.02571+(-0.000391686*x))+(1.01948e-06*(x*x)))+(-1.16475e-09*(x*(x*x)))", ptMin, ptMax);
  TF1* sfFunct_min_0 = new TF1("sfFunct_loose_min_0", "((0.931859+(-0.00045457*x))+(1.25431e-06*(x*x)))+(-1.36433e-09*(x*(x*x)))", ptMin, ptMax);
  TF1* sfFunct_max_0 = new TF1("sfFunct_loose_max_0", "((1.11958+(-0.00032886*x))+(7.84649e-07*(x*x)))+(-9.65161e-10*(x*(x*x)))", ptMin, ptMax);

  TF1* sfFunct_1 = new TF1("sfFunct_loose_1", "((1.03375+(-0.00068776*x))+(2.13443e-06*(x*x)))+(-2.24163e-09*(x*(x*x)))", ptMin, ptMax);
  TF1* sfFunct_min_1 = new TF1("sfFunct_loose_min_1", "((0.936905+(-0.000681017*x))+(2.13885e-06*(x*x)))+(-2.22607e-09*(x*(x*x)))", ptMin, ptMax);
  TF1* sfFunct_max_1 = new TF1("sfFunct_loose_max_1", "((1.13063+(-0.000694616*x))+(2.13001e-06*(x*x)))+(-2.25719e-09*(x*(x*x)))", ptMin, ptMax);

  TF1* sfFunct_2 = new TF1("sfFunct_loose_2", "((1.03597+(-0.000778058*x))+(3.02129e-06*(x*x)))+(-3.0478e-09*(x*(x*x)))", ptMin, ptMax);
  TF1* sfFunct_min_2 = new TF1("sfFunct_loose_min_2", "((0.938438+(-0.00074623*x))+(2.89732e-06*(x*x)))+(-2.92483e-09*(x*(x*x)))", ptMin, ptMax);
  TF1* sfFunct_max_2 = new TF1("sfFunct_loose_max_2", "((1.13355+(-0.000810039*x))+(3.14525e-06*(x*x)))+(-3.17077e-09*(x*(x*x)))", ptMin, ptMax);

  TF1* sfFunct_3 = new TF1("sfFunct_loose_3", "((0.95897+(-0.000111286*x))+(1.6091e-06*(x*x)))+(-2.18387e-09*(x*(x*x)))", ptMin, ptMax);
  TF1* sfFunct_min_3 = new TF1("sfFunct_loose_min_3", "((0.867768+(-9.92078e-05*x))+(1.46903e-06*(x*x)))+(-2.02118e-09*(x*(x*x)))", ptMin, ptMax);
  TF1* sfFunct_max_3 = new TF1("sfFunct_loose_max_3", "((1.0502+(-0.000123474*x))+(1.74917e-06*(x*x)))+(-2.34655e-09*(x*(x*x)))", ptMin, ptMax);

  sfLightFunct_loose_.push_back(sfFunct_0);
  sfLightFunct_loose_min_.push_back(sfFunct_min_0);
  sfLightFunct_loose_max_.push_back(sfFunct_max_0);

  sfLightFunct_loose_.push_back(sfFunct_1);
  sfLightFunct_loose_min_.push_back(sfFunct_min_1);
  sfLightFunct_loose_max_.push_back(sfFunct_max_1);

  sfLightFunct_loose_.push_back(sfFunct_2);
  sfLightFunct_loose_min_.push_back(sfFunct_min_2);
  sfLightFunct_loose_max_.push_back(sfFunct_max_2);

  sfLightFunct_loose_.push_back(sfFunct_3);
  sfLightFunct_loose_min_.push_back(sfFunct_min_3);
  sfLightFunct_loose_max_.push_back(sfFunct_max_3);

}


else if( Atagger == "JPM" ) {

  etaBins_SFLight_medium_.push_back(0.);
  etaBins_SFLight_medium_.push_back(0.8);
  etaBins_SFLight_medium_.push_back(1.6);
  etaBins_SFLight_medium_.push_back(2.4);

  TF1* sfFunct_0 = new TF1("sfFunct_medium_0", "((0.970028+(0.00118179*x))+(-4.23119e-06*(x*x)))+(3.61065e-09*(x*(x*x)))", ptMin, ptMax);
  TF1* sfFunct_min_0 = new TF1("sfFunct_medium_min_0", "((0.840326+(0.000626372*x))+(-2.08293e-06*(x*x)))+(1.57604e-09*(x*(x*x)))", ptMin, ptMax);
  TF1* sfFunct_max_0 = new TF1("sfFunct_medium_max_0", "((1.09966+(0.00173739*x))+(-6.37946e-06*(x*x)))+(5.64527e-09*(x*(x*x)))", ptMin, ptMax);

  TF1* sfFunct_1 = new TF1("sfFunct_medium_1", "((0.918387+(0.000898595*x))+(-2.00643e-06*(x*x)))+(1.26486e-09*(x*(x*x)))", ptMin, ptMax);
  TF1* sfFunct_min_1 = new TF1("sfFunct_medium_min_1", "((0.790843+(0.000548016*x))+(-6.70941e-07*(x*x)))+(1.90355e-11*(x*(x*x)))", ptMin, ptMax);
  TF1* sfFunct_max_1 = new TF1("sfFunct_medium_max_1", "((1.0459+(0.00124924*x))+(-3.34192e-06*(x*x)))+(2.51068e-09*(x*(x*x)))", ptMin, ptMax);

  TF1* sfFunct_2 = new TF1("sfFunct_medium_2", "((0.790103+(0.00117865*x))+(-2.07334e-06*(x*x)))+(1.42608e-09*(x*(x*x)))", ptMin, ptMax);
  TF1* sfFunct_min_2 = new TF1("sfFunct_medium_min_2", "((0.667144+(0.00105593*x))+(-1.43608e-06*(x*x)))+(5.24039e-10*(x*(x*x)))", ptMin, ptMax);
  TF1* sfFunct_max_2 = new TF1("sfFunct_medium_max_2", "((0.913027+(0.00130143*x))+(-2.71061e-06*(x*x)))+(2.32812e-09*(x*(x*x)))", ptMin, ptMax);

  sfLightFunct_medium_.push_back(sfFunct_0);
  sfLightFunct_medium_min_.push_back(sfFunct_min_0);
  sfLightFunct_medium_max_.push_back(sfFunct_max_0);

  sfLightFunct_medium_.push_back(sfFunct_1);
  sfLightFunct_medium_min_.push_back(sfFunct_min_1);
  sfLightFunct_medium_max_.push_back(sfFunct_max_1);

  sfLightFunct_medium_.push_back(sfFunct_2);
  sfLightFunct_medium_min_.push_back(sfFunct_min_2);
  sfLightFunct_medium_max_.push_back(sfFunct_max_2);

}



else if( Atagger == "JPT" )
{

  etaBins_SFLight_tight_.push_back(0.);
  etaBins_SFLight_tight_.push_back(2.4);

  TF1* sfFunct = new TF1("sfFunct_tight", "((0.831392+(0.00269525*x))+(-7.33391e-06*(x*x)))+(5.73942e-09*(x*(x*x)))", ptMin, ptMax);
  TF1* sfFunct_min = new TF1("sfFunct_tight_min", "((0.671888+(0.0020106*x))+(-5.03177e-06*(x*x)))+(3.74225e-09*(x*(x*x)))", ptMin, ptMax);
  TF1* sfFunct_max = new TF1("sfFunct_tight_max", "((0.990774+(0.00338018*x))+(-9.63606e-06*(x*x)))+(7.73659e-09*(x*(x*x)))", ptMin, ptMax);

  sfLightFunct_tight_.push_back(sfFunct);
  sfLightFunct_tight_min_.push_back(sfFunct_min);
  sfLightFunct_tight_max_.push_back(sfFunct_max);

}


//else if( Atagger == "JBPL" )
//{
//  if( fabs(eta) >=0.0 && fabs(eta) <0.5)
//  {
//  sfFunct = new TF1("sfFunct", "((0.996303+(-0.00049586*x))+(1.48662e-06*(x*x)))+(-1.60955e-09*(x*(x*x)))", ptMin, ptMax);
//  sfFunct_min = new TF1("sfFunct_min", "((0.909313+(-0.000483037*x))+(1.48507e-06*(x*x)))+(-1.60327e-09*(x*(x*x)))", ptMin, ptMax);
//  sfFunct_max = new TF1("sfFunct_max", "((1.08332+(-0.000508763*x))+(1.48816e-06*(x*x)))+(-1.61583e-09*(x*(x*x)))", ptMin, ptMax);
//  }
//  else if( fabs(eta) >=0.5 && fabs(eta) <1.0)
//  {
//  sfFunct = new TF1("sfFunct", "((1.01607+(-0.000958122*x))+(3.12318e-06*(x*x)))+(-3.13777e-09*(x*(x*x)))", ptMin, ptMax);
//  sfFunct_min = new TF1("sfFunct_min", "((0.925793+(-0.000877501*x))+(2.88538e-06*(x*x)))+(-2.9089e-09*(x*(x*x)))", ptMin, ptMax);
//  sfFunct_max = new TF1("sfFunct_max", "((1.10639+(-0.0010389*x))+(3.36098e-06*(x*x)))+(-3.36665e-09*(x*(x*x)))", ptMin, ptMax);
//  }
//  else if( fabs(eta) >=1.0 && fabs(eta) <1.5)
//  {
//  sfFunct = new TF1("sfFunct", "((1.04234+(-0.00109152*x))+(3.71686e-06*(x*x)))+(-3.57219e-09*(x*(x*x)))", ptMin, ptMax);
//  sfFunct_min = new TF1("sfFunct_min", "((0.947786+(-0.000985917*x))+(3.39659e-06*(x*x)))+(-3.28635e-09*(x*(x*x)))", ptMin, ptMax);
//  sfFunct_max = new TF1("sfFunct_max", "((1.13696+(-0.00119731*x))+(4.03713e-06*(x*x)))+(-3.85803e-09*(x*(x*x)))", ptMin, ptMax);
//  }
//  else if( fabs(eta) >=1.5 && fabs(eta) <2.4)
//  {
//  sfFunct = new TF1("sfFunct", "((0.960685+(-0.000514241*x))+(2.69297e-06*(x*x)))+(-3.12123e-09*(x*(x*x)))", ptMin, ptMax);
//  sfFunct_min = new TF1("sfFunct_min", "((0.875356+(-0.000455763*x))+(2.42337e-06*(x*x)))+(-2.83637e-09*(x*(x*x)))", ptMin, ptMax);
//  sfFunct_max = new TF1("sfFunct_max", "((1.04606+(-0.000572874*x))+(2.96257e-06*(x*x)))+(-3.40609e-09*(x*(x*x)))", ptMin, ptMax);
//  }
//}
//
//
//else if( Atagger == "JBPM" ) 
//{
//  if( fabs(eta) >=0.0 && fabs(eta) <0.8)
//  {
//  sfFunct = new TF1("sfFunct", "((0.932447+(0.000285676*x))+(-1.03771e-06*(x*x)))+(4.52275e-10*(x*(x*x)))", ptMin, ptMax);
//  sfFunct_min = new TF1("sfFunct_min", "((0.822274+(0.000138316*x))+(-4.14616e-07*(x*x)))+(-9.7638e-11*(x*(x*x)))", ptMin, ptMax);
//  sfFunct_max = new TF1("sfFunct_max", "((1.0426+(0.000433059*x))+(-1.6608e-06*(x*x)))+(1.00219e-09*(x*(x*x)))", ptMin, ptMax);
//  }
//  else if( fabs(eta) >=0.8 && fabs(eta) <1.6)
//  {
//  sfFunct = new TF1("sfFunct", "((0.924959+(0.000170347*x))+(-1.56056e-07*(x*x)))+(-2.06751e-10*(x*(x*x)))", ptMin, ptMax);
//  sfFunct_min = new TF1("sfFunct_min", "((0.822394+(-2.61379e-05*x))+(6.08356e-07*(x*x)))+(-9.28476e-10*(x*(x*x)))", ptMin, ptMax);
//  sfFunct_max = new TF1("sfFunct_max", "((1.02752+(0.000366822*x))+(-9.20467e-07*(x*x)))+(5.14974e-10*(x*(x*x)))", ptMin, ptMax);
//  }
//  else if( fabs(eta) >=1.6 && fabs(eta) <2.4)
//  {
//  sfFunct = new TF1("sfFunct", "((0.846053+(0.000224848*x))+(2.87503e-07*(x*x)))+(-5.93182e-10*(x*(x*x)))", ptMin, ptMax);
//  sfFunct_min = new TF1("sfFunct_min", "((0.714511+(0.000568422*x))+(-7.56289e-07*(x*x)))+(2.61634e-10*(x*(x*x)))", ptMin, ptMax);
//  sfFunct_max = new TF1("sfFunct_max", "((0.977599+(-0.000118755*x))+(1.3313e-06*(x*x)))+(-1.448e-09*(x*(x*x)))", ptMin, ptMax);
//  }
//}
//
//
//
//else if( Atagger == "JBPT" ) {
//  sfFunct = new TF1("sfFunct", "((0.771257+(0.00238891*x))+(-6.2112e-06*(x*x)))+(4.33595e-09*(x*(x*x)))", ptMin, ptMax);
//  sfFunct_min = new TF1("sfFunct_min", "((0.666101+(0.00163462*x))+(-3.92728e-06*(x*x)))+(2.48081e-09*(x*(x*x)))", ptMin, ptMax);
//  sfFunct_max = new TF1("sfFunct_max", "((0.87631+(0.00314343*x))+(-8.49513e-06*(x*x)))+(6.19109e-09*(x*(x*x)))", ptMin, ptMax);
//}




//else if( Atagger == "SSVHEM" ) 
//{
//  if( fabs(eta) >=0.0 && fabs(eta) <0.8)
//  {
//  sfFunct = new TF1("sfFunct", "((0.86318+(0.000801639*x))+(-1.64119e-06*(x*x)))+(2.59121e-10*(x*(x*x)))", ptMin, ptMax);
//  sfFunct_min = new TF1("sfFunct_min", "((0.790364+(0.000463086*x))+(-4.35934e-07*(x*x)))+(-9.08296e-10*(x*(x*x)))", ptMin, ptMax);
//  sfFunct_max = new TF1("sfFunct_max", "((0.935969+(0.0011402*x))+(-2.84645e-06*(x*x)))+(1.42654e-09*(x*(x*x)))", ptMin, ptMax);
//  }
//  else if( fabs(eta) >=0.8 && fabs(eta) <1.6)
//  {
//  sfFunct = new TF1("sfFunct", "((0.958973+(-0.000269555*x))+(1.381e-06*(x*x)))+(-1.87744e-09*(x*(x*x)))", ptMin, ptMax);
//  sfFunct_min = new TF1("sfFunct_min", "((0.865771+(-0.000279908*x))+(1.34144e-06*(x*x)))+(-1.75588e-09*(x*(x*x)))", ptMin, ptMax);
//  sfFunct_max = new TF1("sfFunct_max", "((1.0522+(-0.000259296*x))+(1.42056e-06*(x*x)))+(-1.999e-09*(x*(x*x)))", ptMin, ptMax);
//  }
//  else if( fabs(eta) >=1.6 && fabs(eta) <2.4)
//  {
//  sfFunct = new TF1("sfFunct", "((0.923033+(-0.000898227*x))+(4.74565e-06*(x*x)))+(-6.11053e-09*(x*(x*x)))", ptMin, ptMax);
//  sfFunct_min = new TF1("sfFunct_min", "((0.828021+(-0.000731926*x))+(4.19613e-06*(x*x)))+(-5.81379e-09*(x*(x*x)))", ptMin, ptMax);
//  sfFunct_max = new TF1("sfFunct_max", "((1.01812+(-0.00106483*x))+(5.29518e-06*(x*x)))+(-6.40728e-09*(x*(x*x)))", ptMin, ptMax);
//  }
//}
//
//
//
//else if( Atagger == "SSVHPT" )
//{
//  sfFunct = new TF1("sfFunct", "((0.97409+(0.000646241*x))+(-2.86294e-06*(x*x)))+(2.79484e-09*(x*(x*x)))", ptMin, ptMax);
//  sfFunct_min = new TF1("sfFunct_min", "((0.807222+(0.00103676*x))+(-3.6243e-06*(x*x)))+(3.17368e-09*(x*(x*x)))", ptMin, ptMax);
//  sfFunct_max = new TF1("sfFunct_max", "((1.14091+(0.00025586*x))+(-2.10157e-06*(x*x)))+(2.41599e-09*(x*(x*x)))", ptMin, ptMax);
//}
//
//
//
//
//
//else if( Atagger == "TCHPM" )
//{
//  if( fabs(eta) >=0.0 && fabs(eta) <0.8)
//  {
//  sfFunct = new TF1("sfFunct", "((1.27011+(-0.000869141*x))+(2.49796e-06*(x*x)))+(-2.62962e-09*(x*(x*x)))", ptMin, ptMax);
//  sfFunct_min = new TF1("sfFunct_min", "((1.12949+(-0.000678492*x))+(2.02219e-06*(x*x)))+(-2.21675e-09*(x*(x*x)))", ptMin, ptMax);
//  sfFunct_max = new TF1("sfFunct_max", "((1.41077+(-0.00105992*x))+(2.97373e-06*(x*x)))+(-3.0425e-09*(x*(x*x)))", ptMin, ptMax);
//  }
//  else if( fabs(eta) >=0.8 && fabs(eta) <1.6)
//  {
//  sfFunct = new TF1("sfFunct", "((1.36167+(-0.00153237*x))+(4.54567e-06*(x*x)))+(-4.38874e-09*(x*(x*x)))", ptMin, ptMax);
//  sfFunct_min = new TF1("sfFunct_min", "((1.21289+(-0.00126411*x))+(3.81676e-06*(x*x)))+(-3.75847e-09*(x*(x*x)))", ptMin, ptMax);
//  sfFunct_max = new TF1("sfFunct_max", "((1.51053+(-0.00180085*x))+(5.27457e-06*(x*x)))+(-5.01901e-09*(x*(x*x)))", ptMin, ptMax);
//  }
//  else if( fabs(eta) >=1.6 && fabs(eta) <2.4)
//  {
//  sfFunct = new TF1("sfFunct", "((1.22696+(0.000249231*x))+(9.55279e-08*(x*x)))+(-1.04034e-09*(x*(x*x)))", ptMin, ptMax);
//  sfFunct_min = new TF1("sfFunct_min", "((1.07572+(0.00055366*x))+(-9.55796e-07*(x*x)))+(-3.73943e-11*(x*(x*x)))", ptMin, ptMax);
//  sfFunct_max = new TF1("sfFunct_max", "((1.3782+(-5.52498e-05*x))+(1.14685e-06*(x*x)))+(-2.04329e-09*(x*(x*x)))", ptMin, ptMax);
//  }
//}
//
//
//else if( Atagger == "TCHPT" )
//{
//  sfFunct = new TF1("sfFunct", "((1.20711+(0.000681067*x))+(-1.57062e-06*(x*x)))+(2.83138e-10*(x*(x*x)))", ptMin, ptMax);
//  sfFunct_min = new TF1("sfFunct_min", "((1.03418+(0.000428273*x))+(-5.43024e-07*(x*x)))+(-6.18061e-10*(x*(x*x)))", ptMin, ptMax);
//  sfFunct_max = new TF1("sfFunct_max", "((1.38002+(0.000933875*x))+(-2.59821e-06*(x*x)))+(1.18434e-09*(x*(x*x)))", ptMin, ptMax);
//}

}


void BTagSFUtil::InitMistag(const std::string& wp) {

  TString Atagger = btagAlgo_+wp;

  float ptMin = 20.;
  float ptMax = 670.;

// Definition of functions from plot33New.C ----------------------

if( Atagger == "CSVL" ) {

  etaBins_mistag_loose_.push_back(0.);
  etaBins_mistag_loose_.push_back(0.5);
  etaBins_mistag_loose_.push_back(1.);
  etaBins_mistag_loose_.push_back(1.5);
  etaBins_mistag_loose_.push_back(2.4);

  TF1* tmpMistag_0 = new TF1("tmpMistag_loose_0", "242534*(((1+(0.0182863*x))+(4.50105e-05*(x*x)))/(1+(108569*x)))", ptMin, ptMax);
  TF1* tmpMistag_1 = new TF1("tmpMistag_loose_0", "129.938*(((1+(0.0197657*x))+(4.73472e-05*(x*x)))/(1+(55.2415*x)))", ptMin, ptMax);
  TF1* tmpMistag_2 = new TF1("tmpMistag_loose_0", "592.214*(((1+(0.00671207*x))+(6.46109e-05*(x*x)))/(1+(134.318*x)))", ptMin, ptMax);
  TF1* tmpMistag_3 = new TF1("tmpMistag_loose_0", "93329*(((1+(0.0219705*x))+(3.76566e-05*(x*x)))/(1+(18245.1*x)))", ptMin, ptMax);

  mistagFunct_loose_.push_back(tmpMistag_0);
  mistagFunct_loose_.push_back(tmpMistag_1);
  mistagFunct_loose_.push_back(tmpMistag_2);
  mistagFunct_loose_.push_back(tmpMistag_3);

} 


else if( Atagger == "CSVM" ) {

  etaBins_mistag_medium_.push_back(0.);
  etaBins_mistag_medium_.push_back(0.8);
  etaBins_mistag_medium_.push_back(1.6);
  etaBins_mistag_medium_.push_back(2.4);

  TF1* tmpMistag_0 = new TF1("tmpMistag_medium_0", "(0.00967751+(2.54564e-05*x))+(-6.92256e-10*(x*x))", ptMin, ptMax);
  TF1* tmpMistag_1 = new TF1("tmpMistag_medium_1", "(0.00974141+(5.09503e-05*x))+(2.0641e-08*(x*x))", ptMin, ptMax);
  TF1* tmpMistag_2 = new TF1("tmpMistag_medium_2", "(0.013595+(0.000104538*x))+(-1.36087e-08*(x*x))", ptMin, ptMax);
  
  mistagFunct_medium_.push_back(tmpMistag_0);
  mistagFunct_medium_.push_back(tmpMistag_1);
  mistagFunct_medium_.push_back(tmpMistag_2);

}


else if( Atagger == "CSVT" ) {

  etaBins_mistag_tight_.push_back(0.);
  etaBins_mistag_tight_.push_back(2.4);

  TF1* tmpMistag = new TF1("tmpMistag_tight", "0.00315116*(((1+(-0.00769281*x))+(2.58066e-05*(x*x)))+(-2.02149e-08*(x*(x*x))))", ptMin, ptMax);

  mistagFunct_tight_.push_back(tmpMistag);

}


else if( Atagger == "TCHEL" ) {

  etaBins_mistag_loose_.push_back(0.);
  etaBins_mistag_loose_.push_back(0.5);
  etaBins_mistag_loose_.push_back(1.);
  etaBins_mistag_loose_.push_back(1.5);
  etaBins_mistag_loose_.push_back(2.4);

  TF1* tmpMistag_0 = new TF1("tmpMistag_loose_0", "(((-0.0235318+(0.00268868*x))+(-6.47688e-06*(x*x)))+(7.92087e-09*(x*(x*x))))+(-4.06519e-12*(x*(x*(x*x))))", ptMin, ptMax);
  TF1* tmpMistag_1 = new TF1("tmpMistag_loose_1", "(((-0.0257274+(0.00289337*x))+(-7.48879e-06*(x*x)))+(9.84928e-09*(x*(x*x))))+(-5.40844e-12*(x*(x*(x*x))))", ptMin, ptMax);
  TF1* tmpMistag_2 = new TF1("tmpMistag_loose_2", "(((-0.0310046+(0.00307803*x))+(-7.94145e-06*(x*x)))+(1.06889e-08*(x*(x*x))))+(-6.08971e-12*(x*(x*(x*x))))", ptMin, ptMax);
  TF1* tmpMistag_3 = new TF1("tmpMistag_loose_3", "(((-0.0274561+(0.00301096*x))+(-8.89588e-06*(x*x)))+(1.40142e-08*(x*(x*x))))+(-8.95723e-12*(x*(x*(x*x))))", ptMin, ptMax);

  mistagFunct_loose_.push_back(tmpMistag_0);
  mistagFunct_loose_.push_back(tmpMistag_1);
  mistagFunct_loose_.push_back(tmpMistag_2);
  mistagFunct_loose_.push_back(tmpMistag_3);

}


else if( Atagger == "TCHEM" ) {

  etaBins_mistag_medium_.push_back(0.);
  etaBins_mistag_medium_.push_back(0.8);
  etaBins_mistag_medium_.push_back(1.6);
  etaBins_mistag_medium_.push_back(2.4);

  TF1* tmpMistag_0 = new TF1("tmpMistag_medium_0", "(0.000919586+(0.00026266*x))+(-1.75723e-07*(x*x))", ptMin, ptMax);
  TF1* tmpMistag_1 = new TF1("tmpMistag_medium_1", "(-0.00364137+(0.000350371*x))+(-1.89967e-07*(x*x))", ptMin, ptMax);
  TF1* tmpMistag_2 = new TF1("tmpMistag_medium_2", "(-0.00483904+(0.000367751*x))+(-1.36152e-07*(x*x))", ptMin, ptMax);

  mistagFunct_medium_.push_back(tmpMistag_0);
  mistagFunct_medium_.push_back(tmpMistag_1);
  mistagFunct_medium_.push_back(tmpMistag_2);

}



else if( Atagger == "JPL" ) {

  etaBins_mistag_loose_.push_back(0.);
  etaBins_mistag_loose_.push_back(0.5);
  etaBins_mistag_loose_.push_back(1.);
  etaBins_mistag_loose_.push_back(1.5);
  etaBins_mistag_loose_.push_back(2.4);

  TF1* tmpMistag_0 = new TF1("tmpMistag_loose_0", "(0.060001+(0.000332202*x))+(-2.36709e-07*(x*x))", ptMin, ptMax);
  TF1* tmpMistag_1 = new TF1("tmpMistag_loose_1", "(0.0597675+(0.000370979*x))+(-2.94673e-07*(x*x))", ptMin, ptMax);
  TF1* tmpMistag_2 = new TF1("tmpMistag_loose_2", "(0.0483728+(0.000528418*x))+(-3.17825e-07*(x*x))", ptMin, ptMax);
  TF1* tmpMistag_3 = new TF1("tmpMistag_loose_3", "(0.0463159+(0.000546644*x))+(-3.40486e-07*(x*x))", ptMin, ptMax);

  mistagFunct_loose_.push_back(tmpMistag_0);
  mistagFunct_loose_.push_back(tmpMistag_1);
  mistagFunct_loose_.push_back(tmpMistag_2);
  mistagFunct_loose_.push_back(tmpMistag_3);
  
}


else if( Atagger == "JPM" ) {

  etaBins_mistag_medium_.push_back(0.);
  etaBins_mistag_medium_.push_back(0.8);
  etaBins_mistag_medium_.push_back(1.6);
  etaBins_mistag_medium_.push_back(2.4);

  TF1* tmpMistag_0 = new TF1("tmpMistag_medium_0", "(0.00727084+(4.48901e-05*x))+(-4.42894e-09*(x*x))", ptMin, ptMax);
  TF1* tmpMistag_1 = new TF1("tmpMistag_medium_1", "(0.00389156+(6.35508e-05*x))+(1.54183e-08*(x*x))", ptMin, ptMax);
  TF1* tmpMistag_2 = new TF1("tmpMistag_medium_2", "(0.0032816+(4.18867e-05*x))+(7.44912e-08*(x*x))", ptMin, ptMax);

  mistagFunct_medium_.push_back(tmpMistag_0);
  mistagFunct_medium_.push_back(tmpMistag_1);
  mistagFunct_medium_.push_back(tmpMistag_2);

}


else if( Atagger == "JPT" ) {

  etaBins_mistag_tight_.push_back(0.);
  etaBins_mistag_tight_.push_back(2.4);

  TF1* tmpMistag = new TF1("tmpMistag_tight", "(0.000379966+(8.30969e-06*x))+(1.10364e-08*(x*x))", ptMin, ptMax);
 
  mistagFunct_tight_.push_back(tmpMistag);

}



//else if( Atagger == "JBPL" ) {
//  if( fabs(eta) >= 0.0 && fabs(eta) < 0.5 ) 
//  {
//   tmpMistag = (0.0277261+(0.000808207*x))+(-6.44146e-07*(x*x));
//  }
//  else if( fabs(eta) >= 0.5 && fabs(eta) < 1.0)
//  {
//   tmpMistag = (0.0278926+(0.000827697*x))+(-7.01497e-07*(x*x));
//  }
//  else if( fabs(eta) >= 1.0 && fabs(eta) < 1.5)
//  {
//   tmpMistag = (0.0221411+(0.000900444*x))+(-6.52873e-07*(x*x));
//  }
//  else if( fabs(eta) >= 1.5 && fabs(eta) < 2.4)
//  {
//   tmpMistag = (0.0227045+(0.000808122*x))+(-5.67134e-07*(x*x));
//  }
//
//}
//
//
//else if( Atagger == "JBPM" )
//{
//  if( fabs(eta) >= 0.0 && fabs(eta) < 0.8)
//  {
//   tmpMistag = (((0.00206106+(0.000105851*x))+(2.691e-08*(x*x)))+(-4.34651e-11*(x*(x*x))))+(-6.73107e-14*(x*(x*(x*x))));
//  }
//  else if( fabs(eta) >= 0.8 && fabs(eta) < 1.6)
//  {
//   tmpMistag = (((0.00318438+(4.40327e-05*x))+(3.46922e-07*(x*x)))+(-3.93396e-10*(x*(x*x))))+(3.94283e-14*(x*(x*(x*x))));
//  }
//  else if( fabs(eta) >= 1.6 && fabs(eta) < 2.4)
//  {
//   tmpMistag = (((0.00209833+(4.27753e-05*x))+(1.96076e-07*(x*x)))+(6.19275e-11*(x*(x*x))))+(-2.63318e-13*(x*(x*(x*x))));
//  }
//}
//
//
//else if( Atagger == "JBPT" ) {
// tmpMistag = (-3.36681e-05+(1.37292e-05*x))+(1.78479e-08*(x*x));
//}
//
//
//else if( Atagger == "SSVHEM" ) {
//  if( fabs(eta) >= 0.0 && fabs(eta) < 0.8)
//  {
//   tmpMistag = (((0.000547883+(0.00023023*x))+(-7.31792e-07*(x*x)))+(1.15659e-09*(x*(x*x))))+(-7.00641e-13*(x*(x*(x*x))));
//  }
//  else if( fabs(eta) >= 0.8 && fabs(eta) < 1.6)
//  {
//   tmpMistag = (((0.000615562+(0.000240254*x))+(-7.00237e-07*(x*x)))+(1.2566e-09*(x*(x*x))))+(-8.59011e-13*(x*(x*(x*x))));
//  }
//  else if( fabs(eta) >= 1.6 && fabs(eta) < 2.4)
//  {
//   tmpMistag = (((0.000372388+(0.000309735*x))+(-4.35952e-07*(x*x)))+(3.63763e-10*(x*(x*x))))+(-2.11993e-13*(x*(x*(x*x))));
//  }
//}
//
//
//else if( Atagger == "SSVHPT" ) {
// tmpMistag = (-2.9605e-05+(2.35624e-05*x))+(-1.77552e-08*(x*x));
//}
//
//
//
//else if( Atagger == "TCHPM" ) {
//  if( fabs(eta) >= 0.0 && fabs(eta) < 0.8)
//  {
//   tmpMistag = (((-0.00464673+(0.000247485*x))+(9.13236e-07*(x*x)))+(-2.49994e-09*(x*(x*x))))+(1.65678e-12*(x*(x*(x*x))));
//  }
//  else if( fabs(eta) >= 0.8 && fabs(eta) < 1.6)
//  {
//   tmpMistag = (((-0.0060878+(0.000297422*x))+(1.13369e-06*(x*x)))+(-2.84945e-09*(x*(x*x))))+(1.64721e-12*(x*(x*(x*x))));
//  }
//  else if( fabs(eta) >= 1.6 && fabs(eta) < 2.4)
//  {
//   tmpMistag = (((-0.00836219+(0.000391889*x))+(2.78156e-07*(x*x)))+(-6.14017e-10*(x*(x*x))))+(-1.30592e-13*(x*(x*(x*x))));
//  }
//}
//
//
//else if( Atagger == "TCHPT" ) {
// tmpMistag = (-0.00101+(4.70405e-05*x))+(8.3338e-09*(x*x));
//}


// End of definition of functions from plot33New.C ---------------

}



void BTagSFUtil::InitSFb( const std::string& wp ) {

  float ptMin = 30.;
  float ptMax = 670.;

  if( wp=="L" ) {
    if(btagAlgo_ == "CSV")   SFbFunct_loose_ = new TF1("SFb_loose", "1.02658*((1.+(0.0195388*x))/(1.+(0.0209145*x)))", ptMin, ptMax);
    if(btagAlgo_ == "TCHE")  SFbFunct_loose_ = new TF1("SFb_loose", "0.603913*((1.+(0.286361*x))/(1.+(0.170474*x)))", ptMin, ptMax);
    if(btagAlgo_ == "JP")    SFbFunct_loose_ = new TF1("SFb_loose", "0.969851*((1.+(-6.06362e-05*x))/(1.+(-0.000156638*x)))", ptMin, ptMax);
  } else if( wp=="M" ) {
    if(btagAlgo_ == "TCHE")  SFbFunct_medium_ = new TF1("SFb_medium", "0.932251*((1.+(0.00335634*x))/(1.+(0.00305994*x)))", ptMin, ptMax);
    if(btagAlgo_ == "CSV")   SFbFunct_medium_ = new TF1("SFb_medium", "0.6981*((1.+(0.414063*x))/(1.+(0.300155*x)))", ptMin, ptMax);
    if(btagAlgo_ == "JP")    SFbFunct_medium_ = new TF1("SFb_medium", "0.90806*((1.+(0.000236997*x))/(1.+(5.49455e-05*x)))", ptMin, ptMax);
  }

}


TF1* BTagSFUtil::GetFunctionEtaBins( float eta, const std::vector<float>& etaBins, const std::vector<TF1*>& functions ) const {
   
  int iEta=-1;

  for( unsigned i=0; i<etaBins.size()-1 && iEta<0; ++i ) {
    if( fabs(eta)>=etaBins[i] && fabs(eta)<etaBins[i+1] )
      iEta=i;
  }

  TF1* returnFunction;
  if( iEta<0 ) returnFunction = f1_one_;
  else         returnFunction = functions[iEta];

  return returnFunction;

}



float BTagSFUtil::GetSFLight( float pt, float eta, const std::string& wp, const std::string& meanminmax ) {

  checkInit(wp);

  TF1* thisFunct;

  if( wp=="L" ) {
    if( meanminmax=="mean" )
      thisFunct = GetFunctionEtaBins( eta, etaBins_SFLight_loose_, sfLightFunct_loose_ );
    else if( meanminmax=="min" )
      thisFunct = GetFunctionEtaBins( eta, etaBins_SFLight_loose_, sfLightFunct_loose_min_ );
    else if( meanminmax=="max" )
      thisFunct = GetFunctionEtaBins( eta, etaBins_SFLight_loose_, sfLightFunct_loose_max_ );
  } else if( wp=="M" ) {
    if( meanminmax=="mean" )
      thisFunct = GetFunctionEtaBins( eta, etaBins_SFLight_medium_, sfLightFunct_medium_ );
    else if( meanminmax=="min" )
      thisFunct = GetFunctionEtaBins( eta, etaBins_SFLight_medium_, sfLightFunct_medium_min_ );
    else if( meanminmax=="max" )
      thisFunct = GetFunctionEtaBins( eta, etaBins_SFLight_medium_, sfLightFunct_medium_max_ );
  } else if( wp=="T" ) {
    if( meanminmax=="mean" )
      thisFunct = GetFunctionEtaBins( eta, etaBins_SFLight_tight_, sfLightFunct_tight_ );
    else if( meanminmax=="min" )
      thisFunct = GetFunctionEtaBins( eta, etaBins_SFLight_tight_, sfLightFunct_tight_min_ );
    else if( meanminmax=="max" )
      thisFunct = GetFunctionEtaBins( eta, etaBins_SFLight_tight_, sfLightFunct_tight_max_ );
  }
 
  if( pt<thisFunct->GetXmin() ) pt = thisFunct->GetXmin();
  if( pt>thisFunct->GetXmax() ) pt = thisFunct->GetXmax();

  return thisFunct->Eval(pt);

}




float BTagSFUtil::GetMistag( float pt, float eta, const std::string& wp, const std::string& meanminmax ) {

  checkInit(wp);

  // meanminmax not used for now

  TF1* thisFunct;

  if( wp=="L" ) {
    thisFunct = GetFunctionEtaBins( eta, etaBins_mistag_loose_, mistagFunct_loose_ );
  } else if( wp=="M" ) {
    thisFunct = GetFunctionEtaBins( eta, etaBins_mistag_medium_,mistagFunct_medium_ );
  } else if( wp=="T" ) {
    thisFunct = GetFunctionEtaBins( eta, etaBins_mistag_tight_, mistagFunct_tight_ );
  }
 
  
  if( pt<thisFunct->GetXmin() ) pt = thisFunct->GetXmin();
  if( pt>thisFunct->GetXmax() ) pt = thisFunct->GetXmax();


  // no difference between mean/min/max for now
  return thisFunct->Eval(pt);

}




float BTagSFUtil::GetSFb( float pt, float eta, const std::string& wp, const std::string& meanMinMax ) {

  checkInit(wp);

  TF1* thisFunct;

  if( wp=="L" )
    thisFunct = SFbFunct_loose_;
  else if( wp=="M" )
    thisFunct = SFbFunct_medium_;

  // no difference between mean/min/max for now
  // no eta bins
  if( pt<thisFunct->GetMinimumX() ) pt = thisFunct->GetMinimumX();
  if( pt>thisFunct->GetMaximumX() ) pt = thisFunct->GetMaximumX();

  float SFb_err = 0.;
  if( meanMinMax!="mean" ) { 
    SFb_err = getSFb_err( pt, wp );
    if( meanMinMax=="min" ) SFb_err = -SFb_err;
  }

  float SFb = thisFunct->Eval(pt) + SFb_err;

  return SFb;

}


float BTagSFUtil::getSFb_err( float pt, const std::string wp ) {

  // uncertainties taken from https://twiki.cern.ch/twiki/pub/CMS/BtagPOG/SFb-mujet_payload.txt
  
  std::vector<float> ptBins_min;
  ptBins_min.push_back(30.);
  ptBins_min.push_back(40.);
  ptBins_min.push_back(50.);
  ptBins_min.push_back(60.);
  ptBins_min.push_back(70.);
  ptBins_min.push_back(80.);
  ptBins_min.push_back(100.);
  ptBins_min.push_back(120.);
  ptBins_min.push_back(160.);
  ptBins_min.push_back(210.);
  ptBins_min.push_back(260.);
  ptBins_min.push_back(320.);
  ptBins_min.push_back(400.);
  ptBins_min.push_back(500.);
  ptBins_min.push_back(670.);


  std::vector<float> SFb_err;

  TString Atagger = btagAlgo_+wp;


  if( Atagger=="TCHEL" ) {

    SFb_err.push_back(0.0244956);
    SFb_err.push_back(0.0237293);
    SFb_err.push_back(0.0180131);
    SFb_err.push_back(0.0182411);
    SFb_err.push_back(0.0184592);
    SFb_err.push_back(0.0106444);
    SFb_err.push_back(0.011073);
    SFb_err.push_back(0.0106296);
    SFb_err.push_back(0.0175259);
    SFb_err.push_back(0.0161566);
    SFb_err.push_back(0.0158973);
    SFb_err.push_back(0.0186782);
    SFb_err.push_back(0.0371113);
    SFb_err.push_back(0.0289788 );

  } else if( Atagger=="TCHEM" ) {

    SFb_err.push_back(0.0311456);
    SFb_err.push_back(0.0303825);
    SFb_err.push_back(0.0209488);
    SFb_err.push_back(0.0216987);
    SFb_err.push_back(0.0227149);
    SFb_err.push_back(0.0260294);
    SFb_err.push_back(0.0205766);
    SFb_err.push_back(0.0227065);
    SFb_err.push_back(0.0260481);
    SFb_err.push_back(0.0278001);
    SFb_err.push_back(0.0295361);
    SFb_err.push_back(0.0306555);
    SFb_err.push_back(0.0367805);
    SFb_err.push_back(0.0527368 );

  } else if( Atagger=="CSVL" ) {

    SFb_err.push_back(0.0188743);
    SFb_err.push_back(0.0161816);
    SFb_err.push_back(0.0139824);
    SFb_err.push_back(0.0152644);
    SFb_err.push_back(0.0161226);
    SFb_err.push_back(0.0157396);
    SFb_err.push_back(0.0161619);
    SFb_err.push_back(0.0168747);
    SFb_err.push_back(0.0257175);
    SFb_err.push_back(0.026424);
    SFb_err.push_back(0.0264928);
    SFb_err.push_back(0.0315127);
    SFb_err.push_back(0.030734);
    SFb_err.push_back(0.0438259);

  } else if( Atagger=="CSVM" ) {

    SFb_err.push_back(0.0295675);
    SFb_err.push_back(0.0295095);
    SFb_err.push_back(0.0210867);
    SFb_err.push_back(0.0219349);
    SFb_err.push_back(0.0227033);
    SFb_err.push_back(0.0204062);
    SFb_err.push_back(0.0185857);
    SFb_err.push_back(0.0256242);
    SFb_err.push_back(0.0383341);
    SFb_err.push_back(0.0409675);
    SFb_err.push_back(0.0420284);
    SFb_err.push_back(0.0541299);
    SFb_err.push_back(0.0578761);
    SFb_err.push_back(0.0655432 );

  } else if( Atagger=="CSVT" ) {

    SFb_err.push_back(0.0364717);
    SFb_err.push_back(0.0362281);
    SFb_err.push_back(0.0232876);
    SFb_err.push_back(0.0249618);
    SFb_err.push_back(0.0261482);
    SFb_err.push_back(0.0290466);
    SFb_err.push_back(0.0300033);
    SFb_err.push_back(0.0453252);
    SFb_err.push_back(0.0685143);
    SFb_err.push_back(0.0653621);
    SFb_err.push_back(0.0712586);
    SFb_err.push_back(0.094589);
    SFb_err.push_back(0.0777011);
    SFb_err.push_back(0.0866563 );

  } else if( Atagger=="JPL" ) {

    SFb_err.push_back(0.0250319);
    SFb_err.push_back(0.0250197);
    SFb_err.push_back(0.0212994);
    SFb_err.push_back(0.0225867);
    SFb_err.push_back(0.0239025);
    SFb_err.push_back(0.026476);
    SFb_err.push_back(0.0264219);
    SFb_err.push_back(0.0156582);
    SFb_err.push_back(0.0222798);
    SFb_err.push_back(0.0223169);
    SFb_err.push_back(0.0225454);
    SFb_err.push_back(0.0405975);
    SFb_err.push_back(0.0405668);
    SFb_err.push_back(0.0415829 );

  } else if( Atagger=="JPM" ) {

    SFb_err.push_back(0.0352594);
    SFb_err.push_back(0.0353008);
    SFb_err.push_back(0.0299008);
    SFb_err.push_back(0.0276606);
    SFb_err.push_back(0.0292312);
    SFb_err.push_back(0.0336607);
    SFb_err.push_back(0.0284701);
    SFb_err.push_back(0.029544);
    SFb_err.push_back(0.0358872);
    SFb_err.push_back(0.0367869);
    SFb_err.push_back(0.0375048);
    SFb_err.push_back(0.0597367);
    SFb_err.push_back(0.0653152);
    SFb_err.push_back(0.074242 );

  } else if( Atagger=="JPT" ) {

    SFb_err.push_back(0.0475813);
    SFb_err.push_back(0.0472359);
    SFb_err.push_back(0.0378328);
    SFb_err.push_back(0.0334787);
    SFb_err.push_back(0.034681);
    SFb_err.push_back(0.0398312);
    SFb_err.push_back(0.0481646);
    SFb_err.push_back(0.0392262);
    SFb_err.push_back(0.0463086);
    SFb_err.push_back(0.0534565);
    SFb_err.push_back(0.0545823);
    SFb_err.push_back(0.102505);
    SFb_err.push_back(0.113198);
    SFb_err.push_back(0.138116 );

  } // if atagger


  // find pt bin:
  int ptBin = -1;
  for( unsigned ipt=0; ipt<ptBins_min.size()-1; ++ipt ) {
    if( pt>=ptBins_min[ipt] && pt<=ptBins_min[ipt+1] ) { // <= so that 670 is ok
      ptBin = ipt;
      break;
    }
  }

  return SFb_err[ptBin];

}



void BTagSFUtil::checkInit( const std::string& wp ) {

  bool notOK = (wp=="L" && !init_loose_ ) || 
               (wp=="M" && !init_medium_ ) || 
               (wp=="T" && !init_tight_ );

  if( notOK ) init(wp);

}




/***********************************/
/*           NEW FUNCTION          */
/***********************************/


void BTagSFUtil::SF(const std::string& btagAlgo, const std::string& wp, float pt, float eta){
  int NL = 4; 
  int NM = 3;

  //binning in eta is the same for SFl and mistag, but it changes with the wp

  float etamin_L[4] = {0.0, 0.5, 1.0, 1.5} ;
  float etamax_L[4] = {0.5, 1.0, 1.5, 2.4} ;
  float etamin_M[3] = {0.0, 0.8, 1.6} ;
  float etamax_M[3] = {0.8, 1.6, 2.4} ;
  float etamin(0), etamax(0);
  // SF for b's
  //  float b_SF_Medium = h2_Medium_BTAGBEFFCORR_->GetBinContent( iBin_pt, iBin_eta );
  //  float b_SF_Loose = h2_Loose_BTAGBEFFCORR_->GetBinContent( iBin_pt, iBin_eta );
  //  float SFb_;
  //  float SFlight_, Mistag_;
  SFlightFuncs sfl_func;
  MistagFuncs mt_func;   
  float pt_sfb;
  pt_sfb = pt;

  if(pt>670){ 
    //std::cout<<"WARNING: Mistagging rate for pt greater than 670 are not defined, we are going to treat this case as a pt=670 case"<<std::endl;
    pt = 670;}
  if(pt<30) {
    pt_sfb = 30;
  }

  if(wp == "L"){ 
    if(btagAlgo == "CSV")   SFb_ = 1.02658*((1.+(0.0195388*pt_sfb))/(1.+(0.0209145*pt_sfb)));
    if(btagAlgo == "TCHE") SFb_ = 0.603913*((1.+(0.286361*pt_sfb))/(1.+(0.170474*pt_sfb)));
    if(btagAlgo == "JP") SFb_ = 0.969851*((1.+(-6.06362e-05*pt_sfb))/(1.+(-0.000156638*pt_sfb)));
    for(int i=0;i<NL;++i){
      if ((TMath::Abs(eta) >etamin_L[i] || TMath::Abs(eta) == etamin_L[i]) && TMath::Abs(eta) <etamax_L[i]){
	etamin = etamin_L[i];
	etamax = etamax_L[i];
      }
    }
  }// end L 

 if(wp == "M"){ 
   if(btagAlgo == "TCHE") SFb_ = 0.932251*((1.+(0.00335634*pt_sfb))/(1.+(0.00305994*pt_sfb)));
   if(btagAlgo == "CSV")  SFb_ = 0.6981*((1.+(0.414063*pt_sfb))/(1.+(0.300155*pt_sfb)));
   if(btagAlgo == "JP") SFb_ = 0.90806*((1.+(0.000236997*pt_sfb))/(1.+(5.49455e-05*pt_sfb)));
    for(int i=0;i<NM;++i){
      if (( TMath::Abs(eta) >etamin_M[i] || TMath::Abs(eta) == etamin_L[i]) && TMath::Abs(eta) <etamax_M[i]){
	etamin = etamin_M[i];
	etamax = etamax_M[i];
      }
    }
  }// end M 

 // uncomment for debugging
   //std::cout<<"ETA: "<<eta<<std::endl;
   //std::cout<<"ETAMIN AND ETAMAX ARE: "<<etamin<<" ,"<<etamax<<std::endl;
  TF1* SFlight_func = sfl_func.GetSFlmean(TString(btagAlgo),TString(wp),etamin, etamax);
  SFlight_ = SFlight_func->Eval(pt);
  TF1* Mistag_func = mt_func.GetMistagmean(TString(btagAlgo),TString(wp),etamin, etamax);
  Mistag_ =  Mistag_func->Eval(pt);

  delete SFlight_func;
  delete Mistag_func;

}


/***********************************/
/* THIS FUNCTION IS NOT UP TO DATE */
/***********************************/

void BTagSFUtil::set_fileMedium( TFile* file ) {

  fileMedium_ = file;

  if( fileMedium_==0 ) {
    std::cout << "WARNING!!! File '" << file->GetName() << "' does not exist!! Will not set histograms!!" << std::endl;
    return;
  }

  h2_Medium_BTAGBEFFCORR_ = (TH2D*)fileMedium_->Get("BTAGBEFFCORR");
  h2_Medium_BTAGLEFFCORR_ = (TH2D*)fileMedium_->Get("BTAGLEFFCORR");
  h2_Medium_BTAGLEFF_     = (TH2D*)fileMedium_->Get("BTAGLEFF");

  if( h2_Medium_BTAGBEFFCORR_ == 0 ) {
    std::cout << "WARNING!!! Didn't find histogram 'BTAGBEFFCORR' in file '" << fileMedium_->GetName() << "'!!! Program will crash soon!" << std::endl;
  }
  if( h2_Medium_BTAGLEFFCORR_ == 0 ) {
    std::cout << "WARNING!!! Didn't find histogram 'BTAGLEFFCORR' in file '" << fileMedium_->GetName() << "'!!! Program will crash soon!" << std::endl;
  }
  if( h2_Medium_BTAGLEFF_     == 0 ) {
    std::cout << "WARNING!!! Didn't find histogram 'BTAGLEFF' in file '" << fileMedium_->GetName() << "'!!! Program will crash soon!" << std::endl;
  }

}


/***********************************/
/* THIS FUNCTION IS NOT UP TO DATE */
/***********************************/

void BTagSFUtil::set_fileLoose( TFile* file ) {

  fileLoose_ = file;

  if( fileLoose_==0 ) {
    std::cout << "WARNING!!! File '" << file->GetName() << "' does not exist!! Will not set histograms!!" << std::endl;
    return;
  }

  h2_Loose_BTAGBEFFCORR_ = (TH2D*)fileLoose_->Get("BTAGBEFFCORR");
  h2_Loose_BTAGLEFFCORR_ = (TH2D*)fileLoose_->Get("BTAGLEFFCORR");
  h2_Loose_BTAGLEFF_     = (TH2D*)fileLoose_->Get("BTAGLEFF");

  if( h2_Loose_BTAGBEFFCORR_ == 0 ) {
    std::cout << "WARNING!!! Didn't find histogram 'BTAGBEFFCORR' in file '" << fileLoose_->GetName() << "'!!! Program will crash soon!" << std::endl;
  }
  if( h2_Loose_BTAGLEFFCORR_ == 0 ) {
    std::cout << "WARNING!!! Didn't find histogram 'BTAGLEFFCORR' in file '" << fileLoose_->GetName() << "'!!! Program will crash soon!" << std::endl;
  }
  if( h2_Loose_BTAGLEFF_     == 0 ) {
    std::cout << "WARNING!!! Didn't find histogram 'BTAGLEFF' in file '" << fileLoose_->GetName() << "'!!! Program will crash soon!" << std::endl;
  }

}

/***********************************/
/* THIS FUNCTION IS NOT UP TO DATE */
/***********************************/

void BTagSFUtil::modifyBTagsWithSF( bool& isBTagged_loose, bool& isBTagged_medium, int pdgIdPart, float Btageff_SF_l, float Btageff_SF_m, float Btagmistag_SF_l,float Btagmistag_SF_m, float Btagmistag_eff_l, float Btagmistag_eff_m) {
  

  // b quarks and c quarks:
  if( abs( pdgIdPart ) == 5 ||  abs( pdgIdPart ) == 4) { 

    // no need to downgrade if is a b and not tagged
    if( !isBTagged_loose ) return;

    float coin = rand_->Uniform(1.);
    
    if( isBTagged_medium ){ 
      if( coin > Btageff_SF_m ) {isBTagged_medium=false; isBTagged_loose=false;} //turn medium and loose off, 
    }
    else if( isBTagged_loose && !isBTagged_medium ){
      if( coin > Btageff_SF_l ) {isBTagged_loose=false; }//
    }


  // light quarks:
  } else if( abs( pdgIdPart)>0 ) { //in data it is 0 (save computing time)

    // no need to upgrade if is light and medium tagged
    if( isBTagged_medium ) return;

    float  Btagmistag_SF = 1.0;
    float  Btagmistag_eff = 1.0;

    if( isBTagged_loose ) {
      Btagmistag_SF = Btagmistag_SF_m;
      Btagmistag_eff = Btagmistag_eff_m;
    } else  {
      Btagmistag_SF = Btagmistag_SF_l;
      Btagmistag_eff = Btagmistag_eff_l;
    }

    float mistagPercent = ( Btagmistag_SF*Btagmistag_eff - Btagmistag_eff ) / ( 1. - Btagmistag_eff );

    float coin = rand_->Uniform(1.);

    // for light quarks, the jet has to be upgraded:

    if( coin < mistagPercent ) {
      if( !isBTagged_loose ) {isBTagged_loose = true;}
      else if( !isBTagged_medium ) {isBTagged_medium = true; }

    }
    
  } //if light quark

} //modifyBTagsWithSF






void BTagSFUtil::modifyBTagsWithSF_fast( bool& isBTagged_loose, bool& isBTagged_medium, float jetpt, float jeteta, int pdgIdPart, const std::string& meanMinMax ) {
    if( meanMinMax!="mean" && meanMinMax!="min" && meanMinMax!="max" ) {
    std::cout << "[BTagSFUtil::modifyBTagsWithSF_fast] ERROR! meanMinMax can only be equal to 'mean', 'min' or 'max'. Please fix your code." << std::endl;
    exit(777);
  }

  float  b_SF_Medium = GetSFb(jetpt, jeteta, "M", meanMinMax);
  float light_SF_Medium = GetSFLight(jetpt, jeteta, "M", meanMinMax);
  float light_eff_Medium = GetMistag(jetpt, jeteta, "M", meanMinMax);
  float  b_SF_Loose = GetSFb(jetpt, jeteta, "L", meanMinMax);
  float light_SF_Loose = GetSFLight(jetpt, jeteta, "L", meanMinMax);
  float light_eff_Loose = GetMistag(jetpt, jeteta, "L", meanMinMax);


    
  // b quarks and c quarks:
  if( abs( pdgIdPart ) == 5 ||  abs( pdgIdPart ) == 4) { 

    // no need to downgrade if is a b/c and not tagged
    if( !isBTagged_loose ) return;

    float b_eff = 1.;
    if(abs( pdgIdPart ) == 4) b_eff /= 5.; // c's have 20% eff wrt b's

    if( isBTagged_medium ) { 
      bool newBTag = this->applySF( isBTagged_medium, b_SF_Medium, b_eff );
      if( !newBTag ) {
        isBTagged_medium = newBTag;
        isBTagged_loose = newBTag; //need to downgrade also loose
      }
    } else if( isBTagged_loose && !isBTagged_medium ) {
      bool newBTag = this->applySF( isBTagged_loose, b_SF_Loose, b_eff );
      isBTagged_loose = newBTag;
    }
      


//  float coin = rand_->Uniform(1.);

//  if( isBTagged_medium ){ 
//    if( coin > b_SF_Medium ) {isBTagged_medium=false; isBTagged_loose=false;} //turn medium and loose off, 

//  }
//  else if( isBTagged_loose && !isBTagged_medium ){
//    if( coin > b_SF_Loose ) isBTagged_loose=false; //
//  }


  // light quarks:
  } else if( abs( pdgIdPart)>0 ) { //in data it is 0 (save computing time)

    // no need to upgrade if is light and medium tagged
    if( isBTagged_medium ) return;

    // for light quarks, the jet has to be upgraded:
    if( !isBTagged_loose ) {
      bool newBTag = this->applySF( isBTagged_loose, light_SF_Loose, light_eff_Loose );
      isBTagged_loose = newBTag;
    } else if( !isBTagged_medium ) {
      bool newBTag = this->applySF( isBTagged_medium, light_SF_Medium, light_eff_Medium );
      isBTagged_medium = newBTag;
    }

//  float mistagPercent_Loose = ( light_SF_Loose*light_eff_Loose - light_eff_Loose ) / ( 1. - light_eff_Loose );
//  float mistagPercent_Medium = ( light_SF_Medium*light_eff_Medium - light_eff_Medium ) / ( 1. - light_eff_Medium );


//  float coin = rand_->Uniform(1.);

//  // for light quarks, the jet has to be upgraded:
//  if( !isBTagged_loose ) {
//    if( coin < mistagPercent_Loose ) isBTagged_loose = true;
//  } else if( !isBTagged_medium ) {
//    if( coin < mistagPercent_Medium ) isBTagged_medium = true; 
//  }

    
  } //if light quark

} //modifyBTagsWithSF_fast




bool BTagSFUtil::applySF(bool isBTagged, float Btag_SF, float Btag_eff) const {
  
  bool newBTag = isBTagged;

  if (Btag_SF == 1) return newBTag; //no correction needed 

  //throw die
  float coin = rand_->Uniform(1.);    
  
  if(Btag_SF > 1){  // use this if SF>1

    if( !isBTagged ) {

      //fraction of jets that need to be upgraded
      float mistagPercent = (1.0 - Btag_SF) / (1.0 - (Btag_SF/Btag_eff) );

      //upgrade to tagged
      if( coin < mistagPercent ) {newBTag = true;}
    }

  }else{  // use this if SF<1
      
    //downgrade tagged to untagged
    if( isBTagged && coin > Btag_SF ) {newBTag = false;}

  }

  return newBTag;
}





void BTagSFUtil::modifyBTagsWithSF( const std::string& btagAlgo, bool& isBTagged_loose, bool& isBTagged_medium, float jetpt, float jeteta, int pdgIdPart, float sysSF) {
  

  SF(btagAlgo, "M", jetpt, jeteta);
  float  b_SF_Medium = SFb_;
  float light_SF_Medium = SFlight_;
  float light_eff_Medium = Mistag_;
  SF(btagAlgo, "L", jetpt, jeteta);
  float  b_SF_Loose = SFb_;
  float light_SF_Loose = SFlight_;
  float light_eff_Loose = Mistag_;


  // b quarks and c quarks:
  if( abs( pdgIdPart ) == 5 ||  abs( pdgIdPart ) == 4) { 

    // no need to downgrade if is a b and not tagged
    if( !isBTagged_loose ) return;

    //int iBin_pt = (jetpt<240.) ? h2_Medium_BTAGBEFFCORR_->GetXaxis()->FindBin(jetpt) : h2_Medium_BTAGBEFFCORR_->GetXaxis()->FindBin(239.);
    //int iBin_eta = h2_Medium_BTAGBEFFCORR_->GetYaxis()->FindBin(jeteta);

    // SF for b's

    /* old recipe */
    // float b_SF_Medium = h2_Medium_BTAGBEFFCORR_->GetBinContent( iBin_pt, iBin_eta );
    // float b_SF_Loose = h2_Loose_BTAGBEFFCORR_->GetBinContent( iBin_pt, iBin_eta );



    float coin = rand_->Uniform(1.);
    
    if( isBTagged_medium ){ 
      if( coin > b_SF_Medium ) {isBTagged_medium=false; isBTagged_loose=false;} //turn medium and loose off, 
    }
    else if( isBTagged_loose && !isBTagged_medium ){
      if( coin > b_SF_Loose ) isBTagged_loose=false; //
    }


  // light quarks:
  } else if( abs( pdgIdPart)>0 ) { //in data it is 0 (save computing time)

    // no need to upgrade if is light and medium tagged
    if( isBTagged_medium ) return;

    //int iBin_pt = (jetpt<500.) ? h2_Medium_BTAGBEFFCORR_->GetXaxis()->FindBin(jetpt) : h2_Medium_BTAGBEFFCORR_->GetXaxis()->FindBin(499.);
    //int iBin_eta = h2_Medium_BTAGBEFFCORR_->GetYaxis()->FindBin(jeteta);

    // SF for light quarks
    //    float light_SF_Medium = h2_Medium_BTAGLEFFCORR_->GetBinContent( iBin_pt, iBin_eta );
    //    float light_SF_Loose = h2_Loose_BTAGLEFFCORR_->GetBinContent( iBin_pt, iBin_eta );



    //float light_eff_Medium = h2_Medium_BTAGLEFF_->GetBinContent( iBin_pt, iBin_eta );
    //float light_eff_Loose = h2_Loose_BTAGLEFF_->GetBinContent( iBin_pt, iBin_eta );


    float mistagPercent_Loose = ( light_SF_Loose*light_eff_Loose - light_eff_Loose ) / ( 1. - light_eff_Loose );
    float mistagPercent_Medium = ( light_SF_Medium*light_eff_Medium - light_eff_Medium ) / ( 1. - light_eff_Medium );


    float coin = rand_->Uniform(1.);

    // for light quarks, the jet has to be upgraded:
    if( !isBTagged_loose ) {
      if( coin < mistagPercent_Loose ) isBTagged_loose = true;
    } else if( !isBTagged_medium ) {
      if( coin < mistagPercent_Medium ) isBTagged_medium = true; 
    }

    
  } //if light quark

} //modifyBTagsWithSF



/***********************************/
/* THIS FUNCTION IS NOT UP TO DATE */
/***********************************/

//This method of getSF uses the functional form of the mistag SF based on 2011 data. We do not need an eta information becasuse the dependence is flat in eta. 
BTagScaleFactor BTagSFUtil::getSF( const std::string& type, float jetpt ) {

  BTagScaleFactor btsf;

  float SFlight = -1;

  if ( type == "medium"){
  
    SFlight = 1.23344 - 0.00033012*jetpt - 0.000000227661*jetpt*jetpt + 0.00000000154075*jetpt*jetpt*jetpt;
    btsf.SF = SFlight;

  }
  else{
    
    SFlight = 1.14062 - 0.000456396*jetpt + 0.0000000838135*jetpt*jetpt;
    btsf.SF = SFlight;

  }

  btsf.SF_err = 0.;
  btsf.eff = 1.;
  btsf.eff_err = 0.;

  return btsf;


}


/***********************************/
/* THIS FUNCTION IS NOT UP TO DATE */
/***********************************/

BTagScaleFactor BTagSFUtil::getSF( const std::string& fileName, float jetpt, float jeteta ) {

   BTagScaleFactor btsf;

   bool foundSF = false;

   ifstream ifs(fileName.c_str());

   if(!ifs.good()){
     std::cout<<"WARNING ! Didn't find file "<<fileName.c_str()<<std::endl;
   }

   while( ifs.good() && !foundSF ) {

     float etaMin, etaMax, ptMin, ptMax, eff, eff_err, SF, SF_err;
     ifs >> etaMin >> etaMax >> ptMin >> ptMax >> eff >> eff_err >> SF >> SF_err;

     if( fabs(jeteta)>=etaMin && fabs(jeteta)<etaMax && jetpt>=ptMin && (jetpt<ptMax||ptMax==999.) ) {
       btsf.SF = SF;
       btsf.SF_err = SF_err;
       btsf.eff = eff;
       btsf.eff_err = eff_err;
       foundSF = true;
     }

   } // while

   ifs.close();

   if( !foundSF ) {
     std::cout << "WARNING! Didn't find SF in file '" << fileName << "'. Setting it to 1." << std::endl;
     btsf.SF = 1.;
     btsf.SF_err = 0.;
     btsf.eff = 1.;
     btsf.eff_err = 0.;
   }

   return btsf;

}
