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


struct lepton {
  TLorentzVector p;
  int charge;
  int type; //0==electron 1==muon 2==tau 
  int index;
};


class JZBAnalysis : public UserAnalysisBase{
public:
	JZBAnalysis(TreeReader *tr = NULL, bool isData = false);
	virtual ~JZBAnalysis();
	const bool IsCustomMu(const int);
	const bool IsCustomEl(const int);
	const bool IsCustomJet(const int index);
        const bool passElTriggers(void);
        const bool passMuTriggers(void);

	string outputFileName_; // public name of the output file name

	void Begin();
	void Analyze();
	void End();



private:

        enum counters_t { count_begin, EV=count_begin, MU, EL, JE, count_end };
        Monitor counters[count_end];

        vector<lepton> sortLeptonsByPt(vector<lepton>&);

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

        bool fIsData_;
};
#endif
