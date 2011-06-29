#ifndef RunEfficiency_hh
#define RunEfficiency_hh


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



struct t_lepton {
  TLorentzVector p, pfp;
  float d0PV, dzPV, d0BS, dzBS;
  int charge, pfcharge;
  int type; //0==electron 1==muon 2==tau 
  int index, pfindex;
  float genPt;
  int tag, probe, pprobe;
  int tagReco, probeReco, pprobeReco;
  int tagIso, probeIso, pprobeIso;
  int tagID, probeID, pprobeID;
};



class RunEfficiency : public UserAnalysisBase{
public:
  RunEfficiency(TreeReader *tr=NULL, std::string dataType="mc", bool fullCleaning=false);
  virtual ~RunEfficiency();
  const bool IsCustomMu(const int);
  const bool IsCustomEl(const int);
  const bool IsCustomPfMu(const int);
  const bool IsCustomPfEl(const int);
  const int getPFMuIndex(const int);
  const int getPFElIndex(const int);
  const bool MuPassingTag(const int);
  const bool MuPassingProbe(const int);
  const bool MuPassingPProbe(const int, const int);
  const bool ElPassingTag(const int);
  const bool ElPassingProbe(const int);
  const bool ElPassingPProbe(const int, const int);
  const bool IsCustomJet(const int index);
  const bool passElTriggers(void);
  const bool passEMuTriggers(void);
  const bool passMuTriggers(void);
  const bool IsGoodBasicPFJetPAT3(int, double, double);
  void setInfo();
  string outputFileName_; // public name of the output file name

  void Begin();
  void Analyze();
  void End();

  // Fill generator information
  bool GeneratorInfo();

private:

  enum counters_t { count_begin, EV=count_begin, TR, MU, EL, JE, PJ, count_end };
  Monitor counters[count_end];

  vector<t_lepton> sortLeptonsByPt(vector<t_lepton>&);

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

  TRandom* rand_;

};
#endif
