#include <iostream>
#include <cstdlib>

#include "RooHistError.h"

#include "DrawBase.h"

#include "../include/SSDLPlotter.hh"
#include "ZBiCalculator.h"


float computeZBi( float obs, float b_pred, float b_pred_err );


int main( int argc, char* argv[] ) {

  std::string selectionType = "Apr10_Iso005_NoZVeto_jet20";
  if( argc>1 ) {
    std::string selectionType_str(argv[1]);
    selectionType = selectionType_str;
  }


  // this sets the style:
  DrawBase* db = new DrawBase("OPT_ZBi");
  db->set_lumiOnRightSide();
  db->set_lumiNormalization(4980.);


  SSDLPlotter* plotter = new SSDLPlotter();
  //std::string outputdir = "/shome/pandolf/CMSSW_4_2_8/src/DiLeptonAnalysis/NTupleProducer/macros/" + selectionType;
  //std::string outputdir = "/shome/pandolf/CMSSW_4_2_8/src/DiLeptonAnalysis/NTupleProducer/macros/Apr12_Iso005_NoZVeto_jet20_new";
  std::string outputdir = "/shome/pandolf/CMSSW_4_2_8/src/DiLeptonAnalysis/NTupleProducer/macros/Apr23_Iso005_NoZVeto_jet20";
  plotter->setVerbose(1);
  plotter->fDO_OPT = false;
  plotter->setOutputDir(outputdir);
  plotter->init("/shome/pandolf/CMSSW_4_2_8/src/DiLeptonAnalysis/NTupleProducer/macros/DataCard_SSDL.dat");
  plotter->readHistos(plotter->fOutputFileName);
  plotter->storeWeightedPred();
  //plotter->doAnalysis();

  int min_NJets = 3;
  int min_NBJets = 1;
  int min_NBJets_med = 1;
  float min_ptLept1 = 55.;
  float min_ptLept2 = 30.;
  float min_met = 0.;
  float min_ht = 100.;


  SSDLPrediction ssdlpred =  plotter->makePredictionSignalEvents(min_ht, 7000., min_met, 7000., min_NJets, min_NBJets, min_NBJets_med, min_ptLept1, min_ptLept2, true);

  float b_pred_mm = ssdlpred.bg_mm;
  float b_pred_em = ssdlpred.bg_em;
  float b_pred_ee = ssdlpred.bg_ee;

  float b_pred_mm_err = ssdlpred.bg_mm_err;
  float b_pred_em_err = ssdlpred.bg_em_err;
  float b_pred_ee_err = ssdlpred.bg_ee_err;

  float b_pred     = ssdlpred.bg;
  float b_pred_err = ssdlpred.bg_err;

  float s_ttw_mm = ssdlpred.s_ttw_mm;
  float s_ttw_em = ssdlpred.s_ttw_em;
  float s_ttw_ee = ssdlpred.s_ttw_ee;

  float s_ttz_mm = ssdlpred.s_ttz_mm;
  float s_ttz_em = ssdlpred.s_ttz_em;
  float s_ttz_ee = ssdlpred.s_ttz_ee;

  //int ns_ttw_mm = ssdlpred.ns_ttw_mm;
  //int ns_ttw_em = ssdlpred.ns_ttw_em;
  //int ns_ttw_ee = ssdlpred.ns_ttw_ee;

  //int ns_ttz_mm = ssdlpred.ns_ttz_mm;
  //int ns_ttz_em = ssdlpred.ns_ttz_em;
  //int ns_ttz_ee = ssdlpred.ns_ttz_ee;

  int obs_mm = ssdlpred.obs_mm;
  int obs_em = ssdlpred.obs_em;
  int obs_ee = ssdlpred.obs_ee;

  float s_ttw = s_ttw_mm + s_ttw_em + s_ttw_ee;
  float s_ttz = s_ttz_mm + s_ttz_em + s_ttz_ee;
  float s = s_ttw + s_ttz;

  //int ns_ttw = ns_ttw_mm + ns_ttw_em + ns_ttw_ee;
  //int ns_ttz = ns_ttz_mm + ns_ttz_em + ns_ttz_ee;
  //int ns = ns_ttw + ns_ttz;

  int obs = obs_mm + obs_em + obs_ee;

  // stat error on observed:
  double obs_plus, obs_minus;
  RooHistError::instance().getPoissonInterval(obs,obs_minus,obs_plus,1.);
  double obs_errPlus = obs_plus-obs;
  double obs_errMinus = obs-obs_minus;


  float ZBi = ZBiCalculator::computeZBi( s+b_pred, b_pred, b_pred_err );
  float ZBi_obs = ZBiCalculator::computeZBi( obs, b_pred, b_pred_err );

  float obs_ttWZ = obs - b_pred;

  float nTotal_ttw = 1089608.; //hardwired!!! HORRIBLE!!
  float nTotal_ttz = 1467136.; //hardwired!!! HORRIBLE!!

  float lumi_pb = 4980.;
  float crossSection_ttw = 0.1633;
  float crossSection_ttz = 0.139;

  float n_passed_ttw = nTotal_ttw*(s_ttw_mm+s_ttw_em+s_ttw_ee)/(crossSection_ttw*lumi_pb);
  float n_passed_ttz = nTotal_ttz*(s_ttz_mm+s_ttz_em+s_ttz_ee)/(crossSection_ttz*lumi_pb);

  float eff_ttwz = (n_passed_ttw+n_passed_ttz)/(nTotal_ttw+nTotal_ttz);

  float crossSection = obs_ttWZ / ( lumi_pb*eff_ttwz );

  float xsecErr_stat_plus = obs_errPlus / ( lumi_pb*eff_ttwz );
  float xsecErr_stat_minus = obs_errMinus / ( lumi_pb*eff_ttwz );

  float lumi_err = 0.025*lumi_pb;
  float xsecErr_syst_lumi = lumi_err * crossSection/lumi_pb;

  float xsecErr_syst = xsecErr_syst_lumi; //just lumi for now

  float xsecErr_tot_plus = sqrt( xsecErr_stat_plus*xsecErr_stat_plus + xsecErr_syst*xsecErr_syst );
  float xsecErr_tot_minus = sqrt( xsecErr_stat_minus*xsecErr_stat_minus + xsecErr_syst*xsecErr_syst );

  float significance = crossSection/xsecErr_tot_minus;

  std::cout << std::endl << std::endl << std::endl;
  std::cout << "=========================================" << std::endl;
  std::cout << std::endl;
  std::cout << " Expected BG: " << b_pred << " +- " << b_pred_err << std::endl;
  std::cout << " Expected Signal: " << s << " (eff=" << eff_ttwz*100. << "%)" << std::endl;
  std::cout << " Expected B+S: " << s+b_pred << std::endl;
  std::cout << " Expected ZBi: " << ZBi << std::endl;
  std::cout << " Observed Events: " << obs << std::endl;
  std::cout << " ZBi (observed): " << ZBi_obs << std::endl;
  std::cout << " Observed Signal (BG subtracted): " << obs_ttWZ << std::endl;
  std::cout << " Expected Cross Section: " << (crossSection_ttw+crossSection_ttz) << " pb" << std::endl;
  std::cout << " Measured Cross Section: " << std::endl;
  std::cout << " " <<  crossSection << "  +" << xsecErr_stat_plus << "/-" << xsecErr_stat_minus << " (stat)   +-" << xsecErr_syst << " (syst)  pb" << std::endl;
  std::cout << " " <<  crossSection << "  +" << xsecErr_tot_plus << "/-" << xsecErr_tot_minus << " pb" << std::endl;
  std::cout << std::endl;
  std::cout << "=========================================" << std::endl;
  std::cout << std::endl;

  return 0;

}


