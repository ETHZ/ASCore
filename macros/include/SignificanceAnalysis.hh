#ifndef SignificanceAnalysis_hh
#define SignificanceAnalysis_hh


#include <vector>
#include <cmath>
#include <string>
#include <iostream>
#include <fstream>

#include <TH2D.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TStyle.h>
#include <TLatex.h>
#include <TLeaf.h>
#include <TBranch.h>

#include "base/TreeReader.hh"
#include "base/UserAnalysisBase.hh"

class SignificanceAnalysis : public UserAnalysisBase{
public:
	SignificanceAnalysis(TreeReader *tr = 0);
	virtual ~SignificanceAnalysis();

	void Begin(const char* filename = "SignificancePlots.root");
	void Analyze();
	void FillSigHistos(int part);
	void End();

private:
	double getEta(double x, double y, double z);

	// Significance Plots:
	TFile *fSignHistsFile;
	int fNBinsEta[5];
	int fNBinsPhi;
	TH2D *fH_ptdev[5];
	TH2D *fH_ptsum[5];
	TH2D *fH_pt2sum[5];
	TH2I *fH_ptevt[5];
	TH2D *fH_ptavg[5];
	TH1D *fH_ptsumeta[5];
	TH1I *fH_ptevteta[5];

};
#endif
