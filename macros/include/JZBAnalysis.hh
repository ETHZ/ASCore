#ifndef JZBAnalysis_hh
#define JZBAnalysis_hh


#include <cmath>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>

#include <TH1D.h>
#include <TH1F.h>
#include <TH2F.h>
#include "TLorentzVector.h"

#include "base/TreeReader.hh"
#include "base/UserAnalysisBase.hh"
#include "helper/Monitor.hh"




class JZBAnalysis : public UserAnalysisBase{
public:
	JZBAnalysis(TreeReader *tr = NULL);
	virtual ~JZBAnalysis();
	bool IsCustomMu(int);
	bool IsCustomEl(int);
	bool IsCustomJet(int index);

	string outputFileName_; // public name of the output file name

	void Begin();
	void Analyze();
	void End();



private:

  	template<class T> std::string any2string(T i);
	// file for histograms:
	TFile *fHistFile;
	

	TH1F *fHMee[20];

	TH2F *fHElectronPtEta;
	TH2F *fHElectronIDPtEta;
	TH2F *fHElectronIDIsoPtEta;
	TH2F *fHMeeDPhi;
	TH2F *fHMeePt;
	TH2F *fHMDPhiPt;
	TH2F *fHMZPtJ1Pt;
	TH2F *fHMZPtuJ1Pt;
};
#endif
