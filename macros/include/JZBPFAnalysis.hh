#ifndef JZBPFAnalysis_hh
#define JZBPFAnalysis_hh


#include <cmath>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>

#include <TH1D.h>
#include <TH1F.h>
#include <TH2F.h>
#include "TLorentzVector.h"
#include <TRandom.h>

#include "base/TreeReader.hh"
#include "base/UserAnalysisBase.hh"
#include "helper/Monitor.hh"


struct PFlepton {
  TLorentzVector p;
  TLorentzVector recop;
  int charge; 
  int recocharge;
  int type; //0==electron 1==muon 2==tau 
  int index;
  int recoindex;
  float genPt;
};


class JZBPFAnalysis : public UserAnalysisBase{
public:
  JZBPFAnalysis(TreeReader *tr=NULL, std::string dataType="mc", bool fullCleaning=false, bool isModelScan=false);
  virtual ~JZBPFAnalysis();
  const bool IsCustomMu(const int);
  const bool IsCustomEl(const int);
  const bool IsCustomPfMu(const int);
  const bool IsCustomPfEl(const int);
  const bool IsCustomJet(const int index);
  const bool passElTriggers(void);
  const bool passEMuTriggers(void);
  const bool passMuTriggers(void);
  int getRecoMuIndex(TLorentzVector);
  int getRecoElIndex(TLorentzVector);
  string outputFileName_; // public name of the output file name

  void Begin(TFile *f);
  void Analyze();
  void End(TFile *f);

  // Fill generator information
  void GeneratorInfo();

private:

  enum counters_t { count_begin, EV=count_begin, TR, MU, PFMU, EL, PFEL, JE, PJ, count_end };
  Monitor counters[count_end];

  vector<PFlepton> sortLeptonsByPt(vector<PFlepton>&);

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

  std::string fDataType_;
  bool fFullCleaning_;
  bool fisModelScan;

  TRandom* rand_;

};
#endif
