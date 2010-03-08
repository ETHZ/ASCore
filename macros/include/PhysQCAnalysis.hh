#ifndef PhysQCAnalysis_hh
#define PhysQCAnalysis_hh


#include <vector>
#include <cmath>
#include <string>
#include <iostream>
#include <fstream>
#include <time.h>

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
#include "helper/AnaClass.hh"
#include "base/UserAnalysisBase.hh"
#include "TreeCleaner.hh"

class PhysQCAnalysis : public UserAnalysisBase{
public:
	PhysQCAnalysis(TreeReader *tr = NULL, TreeCleaner *tc = NULL);
	virtual ~PhysQCAnalysis();

	void Begin();
	void Analyze();
	void End();
	void MakePlots(TString, TTree*);
	void PlotTriggerStats();
	void PrintInfoStart(int nEntries);
	double DeltaPhi(double v1, double v2);
	void GetEvtEmChFrac(double & fracEm, double & fracCh);
	
	TreeCleaner *fTC;
	AnaClass *fAC;

private:
	
	TH1D *fMuHistos[3];
	TH1D *fElHistos[14];
	TH1D *fMETHistos[5];
	TH2D *fMETDphi12;

};
#endif
