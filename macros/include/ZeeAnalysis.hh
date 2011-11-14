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
  ZeeAnalysis(TreeReader *tr = NULL, std::string dataType="data");
	virtual ~ZeeAnalysis();

	void Begin();
	void Analyze();
	void End();

private:

  EnergyCorrection *elecorr;
  TLorentzVector CorrElectron(TreeReader *fTR, int i, int mode);
	
	// file for histograms:
	TFile *fHistFile;
  ofstream *myfile[10][10];

  std::string fDataType_;
  bool isdata;

  TH1D *fHElPt;
  TH1D *fHElPtCorr;

  TH1D *fHInvMass0;
  TH1D *fHInvMass15;
  TH1D *fHInvMass16;
  TH1D *fHInvMass17;
  TH1D *fHInvMass18;
  TH1D *fHInvMass20;
  TH1D *fHInvMassEgen;

  /*
  TH1D *fHErecEGen17cat1;
  TH1D *fHErecEGen20cat1;
  TH1D *fHErecEGen17cat2;
  TH1D *fHErecEGen20cat2;
  TH1D *fHErecEGen17cat3;
  TH1D *fHErecEGen20cat3;
  */

  TH1F *fHNumPU;
  TH1F *fHNumVtx;

};
#endif
