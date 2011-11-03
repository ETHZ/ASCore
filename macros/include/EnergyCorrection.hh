#ifndef __ENERGYCORRECTION__
#define __ENERGYCORRECTION__

#include <string>
#include "TFile.h"                                                      
#include "TH1F.h"				                       
#include "TMath.h"
#include "TH3F.h"
#include <iostream>
#include "TF1.h"

#include "base/TreeReader.hh"

#define DBG 0

// CORR TUNING: "photons" or "electrons"

// CORR MODES PHOTONS
// 0: nothing
// 5: SC corr FIT + crack + local
// 6: SC corr FIT + crack (no local)

// CORR MODES ELECTRONS
// 0: default el energy from CMSSW
// 20: SC energy (with el corrections from CMSSW)
// 15: SC corr FIT + crack + local
// 16: SC corr FIT + crack (no local)

// ------------------------------------------------------------------------------------
/**
 * \class EnergyCorrection
 *
 * Defines the interface to be implemented by any smearing function on photon informations
 *
 *
 */

class EnergyCorrection
{
public:
  // ! C-TOR 
  EnergyCorrection(TString tuning);
	
  // ! D-TOR
  virtual ~EnergyCorrection();
	
  double get_correctedenergy(TreeReader *fTR, int photon_index, int mode);

  float getEtaCorrectionBarrel(float eta);

  bool isInPhiCracks(double phi, double eta);

private:

  bool forphotons;
  
  double f5x5( double iEta );


  static const Double_t etaCrackMin = 1.44; 
  static const Double_t etaCrackMax = 1.56;
  static const Int_t    nBinsEta              = 14; 
  Double_t       leftEta  [nBinsEta];
  Double_t       rightEta [nBinsEta];



  double getPho_correctedenergy(TreeReader *fTR, int photon_index, int mode);
  double getEl_correctedenergy(TreeReader *fTR, int photon_index, int mode);


  Double_t xcorr[nBinsEta];
  Double_t par0[nBinsEta];
  Double_t par1[nBinsEta];
  Double_t par2[nBinsEta];
  Double_t par3[nBinsEta];
  Double_t par4[nBinsEta];

  Double_t applyScCorrectionsBrEta_photons(Double_t eta, Double_t sigmaPhiSigmaEta);
  Double_t applyScCorrectionsET_EB_photons(Double_t ET);
  Double_t applyScCorrectionsET_EE_photons(Double_t ET);
  Double_t applyScCorrectionsE_EE_photons(Double_t E);

  Double_t applyScCorrectionsBrEta_electrons(Double_t eta, Double_t sigmaPhiSigmaEta);
  Double_t applyScCorrectionsET_EB_electrons(Double_t ET);
  Double_t applyScCorrectionsET_EE_electrons(Double_t ET);
  Double_t applyScCorrectionsE_EE_electrons(Double_t E);


};

#endif

