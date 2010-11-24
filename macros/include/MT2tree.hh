#ifndef MT2tree_hh
#define MT2tree_hh

#include "TObject.h"
#include "TLorentzVector.h"
#include "TVector3.h"

enum {m_jetSize = 20, m_eleSize = 10, m_muoSize = 10};

class MT2Misc : public TObject {

public:
  MT2Misc();
  virtual ~MT2Misc();

  void Reset();
  
  Int_t    Run;
  Int_t    Event;
  Int_t    LumiSection;
  Int_t    LeptConfig;
  Int_t    NJetsEta5Pt20;
  Bool_t   IsCleanJetEvent;
  Double_t PseudoJetMT2;
  Double_t PseudoJetMCT;
  Double_t PseudoJet1Pt;
  Double_t PseudoJet2Pt;
  Double_t PseudojetAlphaT;
  Double_t Vectorsumpt;
  Double_t PFMETsign;
  Double_t DPhiMhtMpt;
  Double_t HT;

  ClassDef(MT2Misc, 3)
};

class MT2Jet : public TObject {

public:
  MT2Jet();
  virtual ~MT2Jet();

  void Reset();
  void SetLV(const TLorentzVector v);

  TLorentzVector lv;

  Double_t bTagProbTCHE;
  Double_t bTagProbTCHP;
  Double_t bTagProbSSVHE;
  Double_t bTagProbSSVHP;

  Bool_t isPFIDLoose;
  Bool_t isPFIDMedium;
  Bool_t isPFIDTight;

  Double_t ChHadFrac;
  Double_t NeuHadFrac;
  Double_t ChEmFrac;
  Double_t NeuEmFrac;
  Int_t    ChMult;
  Int_t    NeuMult;
  Int_t    NConstituents;

  Int_t    inHemisphere;

  ClassDef(MT2Jet, 3)
};

class MT2tree : public TObject {

public:
  MT2tree();
  virtual ~MT2tree();

  void Reset();

  void SetNJets (int n);
  void SetNJetsIDLoose  (int n);
  void SetNJetsIDMedium (int n);
  void SetNJetsIDTight  (int n);
  void SetNEles (int n);
  void SetNMuons(int n);
  
  // My functions here
  Double_t PseudoJetDPhi ();
  Double_t PseudoJetAngle();
  Double_t JetsDPhi(int j1=1, int j2=0);
  Double_t MetJetDPhi(int ijet = 0);
  Double_t MinMetJetDPhi();
  Int_t    GetNJets(double minJPt=20, std::string PFJID="loose");
  Int_t    MinMetJetDPhiIndex();
  Double_t GetMT2(double testmass=0 , bool massive=false);
  Double_t GetMT2Leading(double testmass=0, bool massive=true, std::string PFJID="loose");
  Double_t GetMT2Hemi(double testmass=0, bool massive=false, std::string PFJID="loose", double minJPt=20, int hemi_association=3);
  Double_t CalcMT2(double testmass, bool massive, 
		  TLorentzVector visible1, TLorentzVector visible2, TLorentzVector MET );
  Double_t GetMCT(bool massive=false);

  MT2Misc misc;
  Int_t   NJets;
  Int_t   NJetsIDLoose;
  Int_t   NJetsIDMedium;
  Int_t   NJetsIDTight;
  Int_t   NEles;
  Int_t   NMuons;
  MT2Jet         jet[m_jetSize];
  TLorentzVector ele[m_eleSize];
  TLorentzVector muo[m_muoSize];
  TLorentzVector pfmet[2];
  TLorentzVector MPT[2];
  TLorentzVector MHTloose[2];
  TLorentzVector pseudoJets[2];

  
  ClassDef(MT2tree, 3)
};

#endif
