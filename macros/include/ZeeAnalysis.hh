#ifndef ZeeAnalysis_hh
#define ZeeAnalysis_hh


#include <cmath>
#include <string>
#include <iostream>
#include <fstream>

#include <TH1D.h>

#include "base/TreeReader.hh"
#include "base/UserAnalysisBase.hh"

#include "TLorentzVector.h"

#include "EnergyCorrection.hh"

class ZeeAnalysis : public UserAnalysisBase{
public:
	ZeeAnalysis(TreeReader *tr = NULL);
	virtual ~ZeeAnalysis();

	void Begin();
	void Analyze();
	void End();

private:

  EnergyCorrection *elecorr;
	
	// file for histograms:
	TFile *fHistFile;
	

	TH1D *fHElPt;
	TH1D *fHElPtCorr;
  TH1D *fHInvMass;
  TH1D *fHInvMassCorr;

};
#endif
