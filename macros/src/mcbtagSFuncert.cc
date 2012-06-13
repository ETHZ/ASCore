#include <math.h>
#include <stdlib.h>
#include <set>
#include "TDatabasePDG.h"
#include "Math/VectorUtil.h"
//#include "CMS2.h"
#include "Math/LorentzVector.h"
#include "TSystem.h"
#include "TFile.h"
#include "TAxis.h"
#include <iostream>
#include "mcbtagSFuncert.h"

//------------------------------------------------------------------------
//------------------------------------------------------------------------
// Code by Claudio - Jan 1, 2012
// Added function for >=3 btags - Jan 16,2012
//
// The function btagEventWeight returns the "event weight" for an event with 
// two or more btags due to the "btagging scale factors". 
// This event_weight is meant to multiply the evt_scale1fb variable and
// whatever other scale factors (eg: lepton efficiency scale factors)
// Note: this only applies if the btag scale factors are < 1
//
// The function btagEventUncertainty returns the "event uncertainty" due to 
// btag scale factors for an event with two or more btags.  It is a relative uncertainty.
// It should be used as follows:
// - Call this function for each event which passes the full selection, including
//   the >=2 btag requirement
// - Multiply the "event uncertainty" by the "event weight" defined above as
//   well evt_scale1fb, etc.  Let's call this quantity "scaled event uncertainty"
//   which is an absolute uncertainty.
// - The sum of "scaled event uncertainty" over all events passing the cuts
//   is the systematic uncertainty on the yield in 1 inv fb.
//
// Added Jan 16, 2012:
// btagEventWeight3 and btagEventUncertainty3 are the versions of the 
// above that are appropriate for >=3 btags.
//
// Added by FG on May 27, 2012:
// overload event weight, uncertainty functions to use "official" scale
// factors, uncertainties as recorded in Tools/btagEff_BTV.h
//
// Mods by Claudio June 7, 2012
// Blow up the errors by a factor of 1.5 according to ICHEP BTV reccommendation 
//-----------------------------------------------------------------------
//----------------------------------------------------------------------
// The btagging scale factor and its uncertainty as a function of pt is 
// hardwired in the functions btagScaleFactor and btagScaleFactorError  
// these values are from the BTV based on 2011 mu-jet data
// https://twiki.cern.ch/twiki/pub/CMS/BtagPOG/SFb-mujet_payload.txt
using namespace std;
double btagScaleFactor(double jetpt, std::string algo) {
    if (algo == "CSVM") {
        float pt = max(min(jetpt, 670.),30.);
        return (0.6981*((1.+(0.414063*pt))/(1.+(0.300155*pt))));
    }
    else {
        if (jetpt < 240.) return 0.96;
        return 0.96;
    }
}
double btagScaleFactorError(double jetpt, std::string algo) {
    if (algo == "CSVM") {
        double ptmin[] = {30, 40, 50, 60, 70, 80, 100, 120, 160, 210, 260, 320, 400, 500};
        double ptmax[] = {40, 50, 60, 70, 80,100, 120, 160, 210, 260, 320, 400, 500, 670};
        double SFb_error[] = {
            0.0295675,
            0.0295095,
            0.0210867,
            0.0219349,
            0.0227033,
            0.0204062,
            0.0185857,
            0.0256242,
            0.0383341,
            0.0409675,
            0.0420284,
            0.0541299,
            0.0578761,
            0.0655432 };
        
	double fudgeFactor=1.5;
        const unsigned int nbins = sizeof(ptmin)/sizeof(double);
        if (jetpt < ptmin[0]) return 0.12*fudgeFactor;
        if (jetpt > ptmax[nbins-1]) return 2*SFb_error[nbins-1]*fudgeFactor;
        for (unsigned int idx = 0; idx < nbins; idx++) {
            if (jetpt > ptmin[idx] && jetpt < ptmax[idx])
                return fudgeFactor*SFb_error[idx];
        }
    }
    else return 0.04;
}



//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
// Two simple utility functions: the min pt for btag jets and the eta range
// This are used to determine status=3 taggable jets
double getMinBtagPt()  {return 40.;}
double getMaxBtagEta() {return 2.4;}
//--------------------------------------------------------------------------
// In order to calculate the "event uncertainty" we need the actual values
// of the btagging efficiencies (for data).  These come from some database 
// or some plots or something.  
// Note: these are meant to be the efficiencies for jets in the fiducial region,
// i.e., something like abs(eta)<2.5.

// The btag efficiency does not need to be perfect, since it is only used for
// calculating uncertainties
double btagEff(double jetpt) {
  return 0.643;
}
//------------------------------------------------------------------------
// Here comes btagEventWeight
// Inputs:
// nbjets = number of reconstructed tagged jets (must be between 2 and 4; if
//          there are 5 or more btag jets, set nbjets=4 and only pass the
//          forst four to this function)
// pt1    = pt of the first  btagged jet
// pt2    = pt of the second btagged jet
// pt3    = pt of the third  btagged jet
// pt4    = pt of the fourth btegged jet
// (note: these do not need to be truth matched)
// Returns a negative number if something goes wrong
double btagEventWeight(int nbjets, double pt1, double pt2, double pt3, double pt4){

  // protect against bad input values
  if (nbjets < 2 || nbjets > 4) {
    std::cout << "Illegal nbjets = " << nbjets << " in btagEventWeight" <<endl;
    return -1.;
  }

  // get the scale factors for these jets;
  double sf1 = btagScaleFactor(pt1);
  double sf2 = btagScaleFactor(pt2);
  double sf3 = btagScaleFactor(pt3);
  double sf4 = btagScaleFactor(pt4);

  // calculate the event weight for the two jet case
  if (nbjets == 2) return sf1*sf2;

  // this is the 3 jet case
  if (nbjets == 3) {
    double temp = sf1*sf2+sf1*sf3+sf2*sf3-2*sf1*sf2*sf3;
    if (temp < 0.) {
      std::cout << "Something wrong in btagEventWeight (3 jet case)" << std::endl;
      return -1.;
    }
    return temp;
  }

  // we get here if there are 4 jets
  double temp = sf1*sf2 + sf1*sf3 + sf1*sf4 + sf2*sf3 + sf2*sf4 + sf3*sf4;
  temp = temp + 3*sf1*sf2*sf3*sf4;
  temp = temp - 2*(sf1*sf2*sf3 + sf1*sf2*sf4 + sf1*sf3*sf4 + sf2*sf3*sf4);
  if (temp < 0.) {
    std::cout << "Something wrong in btagEventWeight (4 jet case)" << std::endl;
    return -1.;
  }
  return temp;
}
//------------------------------------------------------------------------
//------------------------------------------------------------------------
// Here comes btagEventUncertainty.
// Note this is quite approximate, but should be good enough 
// as an uncertainty.
// Inputs:
// nbjets = number of b quarks at status = 3 (must be between 2 and 4)
// pt1    = pt of the first  b quark
// pt2    = pt of the second b quark
// pt3    = pt of the third  b quark
// pt4    = pt of the fourth b quark
// eta1   = eta of the first  b quark
// eta2   = eta of the second b quark
// eta3   = eta of the third  b quark
// eta4   = eta of the fourth b quark
double btagEventUncertainty(int nbjets, double pt1, double eta1, double pt2, double eta2, double pt3, double eta3, double pt4, double eta4) {
  // protect against bad input values
  if (nbjets < 2 || nbjets > 4) {
    std::cout << "Illegal nbjets = " << nbjets << " in btagEventUncertainty" <<endl;
    return -1.;
  }

  // Count fiducial jets and load arrays of efficiencies and errors

  double minpt  = getMinBtagPt();
  double etacut = getMaxBtagEta();
  int mynbjet  = 0;
  double eff[4];
  double effErr[4];
  if (pt1 > minpt && fabs(eta1) < etacut) {
    eff[mynbjet]    =  btagEff(pt1);
    effErr[mynbjet] =  eff[mynbjet]*btagScaleFactorError(pt1)/btagScaleFactor(pt1);
    mynbjet++;
  }
  if (pt2 > minpt && fabs(eta2) < etacut) {
    eff[mynbjet]    =  btagEff(pt2);
    effErr[mynbjet] =  eff[mynbjet]*btagScaleFactorError(pt2)/btagScaleFactor(pt2);
    mynbjet++;
  }
  if (nbjets >= 3 && pt3 > minpt && fabs(eta3) < etacut) {
    eff[mynbjet]    =  btagEff(pt3);
    effErr[mynbjet] =  eff[mynbjet]*btagScaleFactorError(pt3)/btagScaleFactor(pt3);
    mynbjet++;
  }
  if (nbjets == 4 && pt4 > minpt && fabs(eta4) < etacut) {
    eff[mynbjet]    =  btagEff(pt4);
    effErr[mynbjet] =  eff[mynbjet]*btagScaleFactorError(pt4)/btagScaleFactor(pt4);
    mynbjet++;
  }

  // If we find < 2 bquarks, we must be dealing with edge effects 
  // We will pretend that we have 2 bjets at threshold
  if (mynbjet < 2) {
    eff[1]    = btagEff(minpt);
    effErr[1] = eff[1]*btagScaleFactorError(minpt)/btagScaleFactor(minpt);
    if (mynbjet == 0) {
      eff[0]    = eff[1];
      effErr[1] = effErr[0];
    }
    mynbjet = 2;
  }


  // Here for two jets
  if (mynbjet == 2) return (eff[0]*effErr[1]+eff[1]*effErr[0])/(eff[0]*eff[1]);

  // Here for three jets
  if (mynbjet == 3) {
    double eps = eff[0]*eff[1] + eff[0]*eff[2] + eff[1]*eff[2] - 2*eff[0]*eff[1]*eff[2];
    if (eps < 0) {
      std::cout << "Negative eps in btagEventUncertainty (3 jet case)" << std::endl;
      return 0.;
    }
    double temp = (eff[1]+eff[2]-2*eff[1]*eff[2])*effErr[0];
    temp       = (eff[0]+eff[2]-2*eff[0]*eff[2])*effErr[1] + temp; 
    temp       = (eff[0]+eff[1]-2*eff[0]*eff[1])*effErr[2] + temp; 
    return fabs(temp)/eps;
  }
  
  // Here for four jets
  if (mynbjet == 4) {
    double eps = 0.;
    for (int i=0; i<4; i++) {
      for (int j=i+1; j<4; j++) {
	eps = eps + eff[i]*eff[j];
      }
    }
    eps = eps + 3*eff[0]*eff[1]*eff[2]*eff[3];
    eps = eps - 2*(eff[0]*eff[1]*eff[2] + eff[0]*eff[1]*eff[3] + 
                   eff[0]*eff[2]*eff[3] + eff[1]*eff[2]*eff[3]);
    if (eps < 0) {
      std::cout << "Negative eps in btagEventUncertainty (4 jet case)" << std::endl;
      return 0.;
    }
    double temp =  effErr[0]*(eff[1] + eff[2] + eff[3] + 3*eff[1]*eff[2]*eff[3]);
    temp = temp + effErr[1]*(eff[0] + eff[2] + eff[3] + 3*eff[0]*eff[2]*eff[3]);
    temp = temp + effErr[2]*(eff[0] + eff[1] + eff[3] + 3*eff[0]*eff[1]*eff[3]);
    temp = temp + effErr[3]*(eff[0] + eff[1] + eff[2] + 3*eff[0]*eff[1]*eff[2]);
    temp = temp - 2 * effErr[0] * (eff[1]*eff[2] + eff[1]*eff[3] + eff[2]*eff[3]);
    temp = temp - 2 * effErr[1] * (eff[0]*eff[2] + eff[0]*eff[3] + eff[2]*eff[3]);
    temp = temp - 2 * effErr[2] * (eff[0]*eff[1] + eff[0]*eff[3] + eff[1]*eff[3]);
    temp = temp - 2 * effErr[3] * (eff[0]*eff[1] + eff[0]*eff[2] + eff[1]*eff[2]);
    return fabs(temp)/eps;
  }

  // We should never get here!!!
  std::cout << "Something wrong in btagEventUncertainty: mynbjet = " << mynbjet << std::endl; 
  return 0;
}
//------------------------------------------------------------------------
// Here comes btagEventWeight3
// Inputs:
// nbjets = number of reconstructed tagged jets (must be between 3 and 4; if
//          there are 5 or more btag jets, set nbjets=4 and only pass the
//          forst four to this function)
// pt1    = pt of the first  btagged jet
// pt2    = pt of the second btagged jet
// pt3    = pt of the third  btagged jet
// pt4    = pt of the fourth btegged jet
// (note: these do not need to be truth matched)
// Returns a negative number if something goes wrong
double btagEventWeight3(int nbjets, double pt1, double pt2, double pt3, double pt4){

  // protect against bad input values
  if (nbjets < 3 || nbjets > 4) {
    std::cout << "Illegal nbjets = " << nbjets << " in btagEventWeight3" <<endl;
    return -1.;
  }

  // get the scale factors for these jets;
  double sf1 = btagScaleFactor(pt1);
  double sf2 = btagScaleFactor(pt2);
  double sf3 = btagScaleFactor(pt3);
  double sf4 = btagScaleFactor(pt4);

  // calculate the event weight for the three jet case
  if (nbjets == 3) return sf1*sf2*sf3;

  // this is the 4 jet case
  if (nbjets == 4) {
    double temp = sf1*sf2*sf3*sf4;
    temp = temp + (1-sf1)*sf2*sf3*sf4;
    temp = temp + (1-sf2)*sf1*sf3*sf4;
    temp = temp + (1-sf3)*sf1*sf2*sf4;
    temp = temp + (1-sf4)*sf1*sf2*sf3;
    if (temp < 0.) {
      std::cout << "Something wrong in btagEventWeight3 (4 jet case)" << std::endl;
      return -1.;
    }
    return temp;
  }
  // should never get here
  return -1.;
}
//------------------------------------------------------------------------
//------------------------------------------------------------------------
// Here comes btagEventUncertainty3
// Inputs:
// nbjets = number of b quarks at status = 3 (must be between 3 and 4)
// pt1    = pt of the first  b quark
// pt2    = pt of the second b quark
// pt3    = pt of the third  b quark
// pt4    = pt of the fourth b quark
// eta1   = eta of the first  b quark
// eta2   = eta of the second b quark
// eta3   = eta of the third  b quark
// eta4   = eta of the fourth b quark
double btagEventUncertainty3(int nbjets, double pt1, double eta1, double pt2, double eta2, double pt3, double eta3, double pt4, double eta4) {
  // protect against bad input values
  if (nbjets < 3 || nbjets > 4) {
    std::cout << "Illegal nbjets = " << nbjets << " in btagEventUncertainty3" <<endl;
    return -1.;
  }

  // Count fiducial jets and load arrays of efficiencies and errors
  double minpt  = getMinBtagPt();
  double etacut = getMaxBtagEta();
  int mynbjet  = 0;
  double eff[4];
  double effErr[4];
  if (pt1 > minpt && fabs(eta1) < etacut) {
    eff[mynbjet]    =  btagEff(pt1);
    effErr[mynbjet] =  eff[mynbjet]*btagScaleFactorError(pt1)/btagScaleFactor(pt1);
    mynbjet++;
  }
  if (pt2 > minpt && fabs(eta2) < etacut) {
    eff[mynbjet]    =  btagEff(pt2);
    effErr[mynbjet] =  eff[mynbjet]*btagScaleFactorError(pt2)/btagScaleFactor(pt2);
    mynbjet++;
  }
  if (nbjets >= 3 && pt3 > minpt && fabs(eta3) < etacut) {
    eff[mynbjet]    =  btagEff(pt3);
    effErr[mynbjet] =  eff[mynbjet]*btagScaleFactorError(pt3)/btagScaleFactor(pt3);
    mynbjet++;
  }
  if (nbjets == 4 && pt4 > minpt && fabs(eta4) < etacut) {
    eff[mynbjet]    =  btagEff(pt4);
    effErr[mynbjet] =  eff[mynbjet]*btagScaleFactorError(pt4)/btagScaleFactor(pt4);
    mynbjet++;
  }

  // If we find < 3 bquarks, we must be dealing with edge effects 
  // We will pretend that we have 3 bjets at threshold
  if (mynbjet < 3) {
    eff[2]    = btagEff(minpt);
    effErr[2] = eff[2]*btagScaleFactorError(minpt)/btagScaleFactor(minpt);
    if (mynbjet < 2) {
      eff[1]    = btagEff(minpt);
      effErr[1] = eff[1]*btagScaleFactorError(minpt)/btagScaleFactor(minpt);
      if (mynbjet == 0) {
          eff[0]    = btagEff(minpt);
          effErr[0] = eff[0]*btagScaleFactorError(minpt)/btagScaleFactor(minpt);
      }
    }
    mynbjet = 3;
  }

  // Here for three jets
  if (mynbjet == 3) return (eff[0]*eff[1]*effErr[2]+eff[0]*eff[2]*effErr[1]+eff[1]*eff[2]*effErr[3])/(eff[0]*eff[1]*eff[2]);

  // Here for four jets  
  if (mynbjet == 4) {
    double eps = eff[0]*eff[1]*eff[2]*eff[3];
    eps = eps + (1-eff[0])*eff[1]*eff[2]*eff[3];
    eps = eps + (1-eff[1])*eff[0]*eff[2]*eff[3];
    eps = eps + (1-eff[2])*eff[0]*eff[1]*eff[3];
    eps = eps + (1-eff[3])*eff[0]*eff[1]*eff[2];
    if (eps < 0) {
      std::cout << "Negative eps in btagEventUncertainty (4 jet case)" << std::endl;
      return 0.;
    }
    double temp =  (eff[1]*eff[2] + eff[1]*eff[3] + eff[2]*eff[3] - 3*eff[1]*eff[2]*eff[3]) * effErr[0];
    temp = temp + (eff[0]*eff[2] + eff[0]*eff[3] + eff[2]*eff[3] - 3*eff[0]*eff[2]*eff[3]) * effErr[1];
    temp = temp + (eff[1]*eff[0] + eff[1]*eff[3] + eff[0]*eff[3] - 3*eff[0]*eff[1]*eff[3]) * effErr[2];
    temp = temp + (eff[1]*eff[0] + eff[1]*eff[2] + eff[0]*eff[2] - 3*eff[0]*eff[1]*eff[2]) * effErr[3];

    return fabs(temp)/eps;
  }

  // We should never get here!!!
  std::cout << "Something wrong in btagEventUncertainty3: mynbjet = " << mynbjet << std::endl; 
  return 0;
}
