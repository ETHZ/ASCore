#ifndef HggAnalysis_hh
#define HggAnalysis_hh


#include <cmath>
#include <string>
#include <iostream>
#include <fstream>

#include <TH1D.h>

#include "base/TreeReader.hh"
#include "base/UserAnalysisBase.hh"

#include "TLorentzVector.h"

#include "EnergyCorrection.hh"

class HggAnalysis : public UserAnalysisBase{
public:
	HggAnalysis(TreeReader *tr = NULL);
	virtual ~HggAnalysis();

	void Begin();
	void Analyze();
	void End();

private:

  EnergyCorrection *phocorr;
	
	// file for histograms:
	TFile *fHistFile;
	

	TH1D *fHPhoPt;
	TH1D *fHPhoPtCorr;
  TH1D *fHInvMass;
  TH1D *fHInvMassCorr;

};
#endif
