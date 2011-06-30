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
  int tag[5], probe[5], pprobe[5];
  int tagReco[5], probeReco[5], pprobeReco[5];
  int tagIso[5], probeIso[5], pprobeIso[5];
  int tagID[5], probeID[5], pprobeID[5];
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
  const bool MuPassingTag(int, int);
  const bool MuPassingIsoTag(int, int);
  const bool MuPassingIsoProbe(int, int);
  const bool MuPassingIsoPProbe(int, int, int);
  const bool MuPassingIDTag(int, int);
  const bool MuPassingIDProbe(int, int);
  const bool MuPassingRecoTag(int, int);
  const bool MuPassingRecoProbe(int, int, int);
  const bool MuPassingRecoPProbe(int, int, int);
  const bool MuPassingProbe(int,int);
  const bool MuPassingPProbe(const int, const int);
  const bool ElPassingTag(const int);
  const bool ElPassingTag(const int, const int);
  const bool ElPassingProbe(const int);
  const bool ElPassingPProbe(const int, const int);

  const bool ElPassingRecoTag(int,int);
  const bool ElPassingRecoProbe(int);
  const bool ElPassingRecoProbe(int,int);
  const bool ElPassingRecoPProbe(const int, const int);

  const bool ElPassingIsoTag(int,int);
  const bool ElPassingIsoProbe(int,int);
  const bool ElPassingIsoPProbe(const int, const int);

  const bool ElPassingIDTag(int,int);
  const bool ElPassingIDProbe(int,int);
  const bool ElPassingIDPProbe(int, int, int);

  const bool IsCustomJet(const int index);
  const bool passElTriggers(int);
  const bool passEMuTriggers(int);
  const bool passMuTriggers(int);
  const bool passAnyMT2Trigger();
  const bool MuPassingRecoProbe(int,int);
  const bool MuPassingIDPProbe(int, int,int);
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
