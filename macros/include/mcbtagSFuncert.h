#ifndef MCBTAGSFUNCERT_H
#define MCBTAGSFUNCERT_H

enum SMSFastSim {
  SMS_None = 0,
  SMS_T1tttt = 1
};
// if useFastSim is enabled, follow prescription at https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagPOG#2011_FastSim_Correction_Factors
double btagScaleFactor(double jetpt, std::string algo = "CSVM", bool useFastSim = false);
double btagScaleFactorError(double jetpt, std::string algo = "CSVM", bool useFastSim = false, SMSFastSim systType  = SMS_None);

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
// Two simple utility functions: the min pt for btag jets and the eta range
// This are used to determine status=3 taggable jets
double getMinBtagPt();
double getMaxBtagEta();

//--------------------------------------------------------------------------
// In order to calculate the "event uncertainty" we need the actual values
// of the btagging efficiencies (for data).  These come from some database 
// or some plots or something.  
// THE FUNCTION BELOW IS JUST A PLACE HOLDER.  PLEASE FIX IT.
// Note: these are meant to be the efficiencies for jets in the fiducial region,
// i.e., something like abs(eta)<2.5.
// The btag efficiency does not need to be perfect, since it is only used for
// calculating uncertainties
double btagEff(double jetpt);

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
double btagEventWeight(int nbjets, double pt1, double pt2, double pt3=0., double pt4=0., bool useFastSim = false);

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
double btagEventUncertainty(int nbjets, double pt1, double eta1, double pt2, double eta2, double pt3=0., double eta3=0., double pt4=0., double eta4=0., bool useFastSim = false, SMSFastSim systType = SMS_None);

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
double btagEventWeight3(int nbjets, double pt1, double pt2, double pt3, double pt4=0., 
			bool useFastSim = false);

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
double btagEventUncertainty3(int nbjets, double pt1, double eta1, double pt2, double eta2, double pt3, double eta3, double pt4=0., double eta4=0.,
			     bool useFastSim = false, SMSFastSim systType = SMS_None);

#endif
