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
#include <TRandom.h>

#include "base/TreeReader.hh"
#include "base/UserAnalysisBase.hh"
#include "helper/Monitor.hh"


struct lepton {
  TLorentzVector p;
  int charge;
  int type; //0==electron 1==muon 2==tau 
  int index;
  float genPt;
  float iso;
  bool ElCInfoIsGsfCtfCons;
  bool ElCInfoIsGsfCtfScPixCons;
  bool ElCInfoIsGsfScPixCons;
};


class JZBAnalysis : public UserAnalysisBase{
public:
  JZBAnalysis(TreeReader *tr=NULL, std::string dataType="mc", bool fullCleaning=false, bool isModelScan=false, bool makeSmall=false, bool doGenInfo=false);
  virtual ~JZBAnalysis();
  const bool IsCustomMu(const int);
  const bool IsCustomEl(const int);
  const bool IsCustomPfMu(const int, const int);
  const bool IsCustomPfEl(const int, const int);
  const bool IsCustomJet(const int);
  const bool IsConvertedPhoton( const int eIndex );
  const bool passTriggers(std::vector<std::string>& triggerPaths);

  string outputFileName_; // public name of the output file name

  void Begin(TFile *f);
  void Analyze();
  void End(TFile *f);

  // Fill generator information
  void GeneratorInfo();

private:

  //enum counters_t { count_begin, EV=count_begin, TR, MU, EL, JE, PJ, count_end };
  enum counters_t { count_begin, EV=count_begin, TR, MU, PFMU, EL, PFEL, JE, PJ, count_end };

  Monitor counters[count_end];

  vector<lepton> sortLeptonsByPt(vector<lepton>&);
  // Add a set of paths to path container
  void addPath(std::vector<std::string>& paths,std::string base, 
               unsigned int start, unsigned int end);

  template<class T> std::string any2string(T i);

  std::string fDataType_;
  bool fFullCleaning_;
  bool fisModelScan;
  bool fdoGenInfo;
  bool fmakeSmall;

  std::vector<std::string> elTriggerPaths, muTriggerPaths, emTriggerPaths;

  TRandom* rand_;

};
#endif
