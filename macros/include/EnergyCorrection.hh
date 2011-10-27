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
// 1: SC corr + crack + local - NOT IMPLEMENTED
// 2: SC corr + crack (no local) - NOT IMPLEMENTED
// 5: SC corr FIT + crack + local
// 6: SC corr FIT + crack (no local)

// CORR MODES ELECTRONS
// 0: default el energy from CMSSW
// 20: SC energy + Ceta + preshower
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

private:

  bool forphotons;
  
	double f5x5( double iEta );
	float getEtaCorrectionBarrel(float eta);
	bool isInPhiCracks(double phi, double eta);
	
	static const Int_t    nBinsEta = 14;                    
	Double_t       leftEta  [nBinsEta];                   
	Double_t       rightEta [nBinsEta];                   
	static const Int_t    nBinsBr = 18;                     
	Double_t       leftBr  [nBinsBr];                      
	Double_t       rightBr [nBinsBr];                      
	Double_t brbins  [2*nBinsBr];    
	TH1F *h_corr[nBinsEta];   
	static const Int_t    nBinsET = 14;             
	Double_t       leftET  [nBinsET];                    
	Double_t       rightET [nBinsET];                    
	Double_t       ETBins  [nBinsET*2];     
	TH1F *h_CBET_EB;
	TH1F *h_CBET_EE;

	double getPho_correctedenergy(TreeReader *fTR, int photon_index, int mode);
	double getEl_correctedenergy(TreeReader *fTR, int photon_index, int mode);

  /*
  // binned corrections
	Double_t applyScCorrectionsBrEta(Double_t eta, Double_t sigmaPhiSigmaEta);
	Double_t applyScCorrectionsET_EB(Double_t ET);
	Double_t applyScCorrectionsET_EE(Double_t ET);
  */

  // fitted corrections
	Double_t xcorr[nBinsEta];
	TF1 *fcorr[nBinsEta]; 
	Double_t applyScCorrectionsBrEta_FIT(Double_t eta, Double_t sigmaPhiSigmaEta);
	Double_t applyScCorrectionsET_EB_FIT(Double_t ET);
	Double_t applyScCorrectionsET_EE_FIT(Double_t ET);

	Double_t applyScCorrectionsET_EB_FIT_electrons(Double_t ET);
	Double_t applyScCorrectionsET_EE_FIT_electrons(Double_t ET);

};

#endif

