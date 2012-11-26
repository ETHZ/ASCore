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
  JZBAnalysis(TreeReader *tr=NULL, std::string dataType="mc", bool fullCleaning=false, bool isModelScan=false, bool makeSmall=false, bool doGenInfo=false, vector<string> fileList=vector<string>());
  virtual ~JZBAnalysis();
  const bool IsCustomMu2012(const int);
  const bool IsCustomEl2012(const int);
  const bool IsCustomPhoton2012(const int);
  const float EffArea(float); //Used for the calculation of Electron isolation
  const float IndividualEffArea(float abseta, string type);
  const bool IsCustomJet(const int);
  const bool IsConvertedPhoton( const int eIndex );
  const bool passTriggers(std::vector<std::string>& triggerPaths);
  const bool passFilters(int& bits);
  const float GetLeptonWeight(int id1, float phi1, float eta1, int id2, float phi2, float eta2, float &EffErr);
  const float GetMuonWeight(float eta1, float pt1, float &EffErr);
  const float GetElectronWeight(float eta1, float pt1, float &EffErr);
  int DoFSRRecovery(TLorentzVector &tmpVector,vector<lepton> &photons);
  void StoreAllPhotons(vector<lepton> &photons, lepton &lepton1, lepton &lepton2);
  float ElPFIso(int index);
  float MuPFIso(int index);
  float PhoPFIso(int index);
  bool DecaysToTaus(bool, TreeReader*);
  int DetermineFlavor(bool, TreeReader*);
  bool ShouldPhotonBeMerged(lepton &photon, float dR);
  int ExtractFileNumber(string fileName);
  float GetBWeight(string WP,int JetFlavor, float JetPt, float JetEta, float &Uncert);
  string outputFileName_; // public name of the output file name

  void Begin(TFile *f);
  void Analyze();
  void End(TFile *f);

  // Fill generator information
  void GeneratorInfo();
  
  void DoFSRStudy(bool fdoGenInfo,TreeReader *fTR);

private:

  //enum counters_t { count_begin, EV=count_begin, TR, MU, EL, JE, PJ, count_end };
  enum counters_t { count_begin, EV=count_begin, TR, MU, PFMU, EL, PFEL, JE, PJ, PH, count_end };

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
  int fFile;

  TFile *CSVT_CorrectionFile;
  TFile *CSVM_CorrectionFile;
  TFile *CSVL_CorrectionFile;
  
  TH2F *CSVT_EfficiencyCorrection;
  TH2F *CSVT_MisTagCorrection;
  TH2F *CSVT_EfficiencyCorrectionUncert;
  TH2F *CSVT_MisTagCorrectionUncert;
  
  TH2F *CSVM_EfficiencyCorrection;
  TH2F *CSVM_MisTagCorrection;
  TH2F *CSVM_EfficiencyCorrectionUncert;
  TH2F *CSVM_MisTagCorrectionUncert;
  
  TH2F *CSVL_EfficiencyCorrection;
  TH2F *CSVL_MisTagCorrection;
  TH2F *CSVL_EfficiencyCorrectionUncert;
  TH2F *CSVL_MisTagCorrectionUncert;
  
  std::vector<std::string> elTriggerPaths, muTriggerPaths, emTriggerPaths, meTriggerPaths, metTriggerPaths, htTriggerPaths;

  TRandom* rand_;

};
#endif
