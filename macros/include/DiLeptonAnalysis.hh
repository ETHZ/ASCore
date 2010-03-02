#ifndef DiLeptonAnalysis_hh
#define DiLeptonAnalysis_hh


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
#include "helper/Davismt2.h"

class DiLeptonAnalysis : public UserAnalysisBase{
public:
	DiLeptonAnalysis(TreeReader *tr = 0);
	virtual ~DiLeptonAnalysis();

	void Begin(const char* filename = "DiLepTree.root");
	void Analyze();
	void Reset();
	void End();

	Davismt2 *fMT2;

private:
	TFile *fDiLepTreeFile_;
	TTree *fDiLepTree_;

	int fTRunNumber;
	int fTEventNumber;
	int fTLumiSection;

	int fTMu1charge;
	int fTMu2charge;
	int fTNqualmu;
	double fTMu1pt;
	double fTMu2pt;
	double fTMu1eta;
	double fTMu2eta;
	double fTMu1iso;
	double fTMu2iso;
	double fTMu1d0;
	double fTMu2d0;
	double fTMu1ntkhits;
	double fTMu2ntkhits;
	double fTMuminv;
	double fTMumt2_50;
	double fTMumt2_100;

	int fTNqualel;
	int fTEl1charge;
	int fTEl2charge;
	double fTEl1pt;
	double fTEl2pt;
	double fTEl1eta;
	double fTEl2eta;
	double fTEl1iso;
	double fTEl2iso;
	double fTEl1d0;
	double fTEl2d0;
	double fTElminv;
	double fTElmt2_50;
	double fTElmt2_100;

};
#endif
