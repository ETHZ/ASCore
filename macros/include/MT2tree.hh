#ifndef MT2tree_hh
#define MT2tree_hh

#include "TObject.h"
#include "TLorentzVector.h"
#include "TVector3.h"

enum {m_jetSize = 40, m_eleSize = 15, m_muoSize = 15};

// MT2Misc ----------------------------------
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
  Bool_t   HBHENoiseFlag;
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

// MT2Jet ----------------------------------
class MT2Jet : public TObject {

public:
  MT2Jet();
  virtual ~MT2Jet();

  void Reset();
  void SetLV(const TLorentzVector v);
  Bool_t IsGoodPFJet(double minJPt=20., double maxJEta=2.4, int PFJID=1); // PFJID: 1 - loose, 2 - medium, 3 - tight

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

// MT2Elec ----------------------------------
class MT2Elec : public TObject {

public:
  MT2Elec();
  virtual ~MT2Elec();

  void Reset();
  void SetLV(const TLorentzVector v);

  TLorentzVector lv;

  Bool_t isTight;

  ClassDef(MT2Elec, 1)
};

// MT2Muon ----------------------------------
class MT2Muon : public TObject {

public:
  MT2Muon();
  virtual ~MT2Muon();

  void Reset();
  void SetLV(const TLorentzVector v);

  TLorentzVector lv;

  Bool_t isTight;

  ClassDef(MT2Muon, 1)
};


// MT2tree ----------------------------------
class MT2tree : public TObject {

public:
  MT2tree();
  virtual ~MT2tree();

  void Reset();

  void SetNJets         (int n);
  void SetNJetsIDLoose  (int n);
  void SetNJetsIDMedium (int n);
  void SetNJetsIDTight  (int n);
  void SetNEles         (int n);
  void SetNElesLoose    (int n);
  void SetNMuons        (int n);
  void SetNMuonsLoose   (int n);
  
  // My functions here
  // NJets
  Int_t    GetNJets(double minJPt=20, double maxJEta=5., int PFJID=1);  // PFJETID includes ptmin (20GeV) and etamax (2.4)
  Int_t    GetNjets(double minJPt=20, double maxJEta=5., int PFJID=0);  // PFJETID not depends on pt and eta
  Int_t    GetJetIndex(int ijet=0, int PFJID=1);
  Double_t JetPt      (int ijet=0, int PFJID=1);
  // dPhi and friends
  Double_t PseudoJetDPhi ();
  Double_t PseudoJetAngle();
  Double_t JetsDPhi(int j1=1, int j2=0, int PFJID=1);
  Double_t MetJetDPhi(int ijet = 0, int PFJID=1, int met=1);
  Double_t MinMetJetDPhi     (int PFJID=1, double minJPt=20, double maxJEta=6., int met=1);
  Int_t    MinMetJetDPhiIndex(int PFJID=1, double minJPt=20, double maxJEta=6., int met=1);
  // MT2 & friends
  Double_t GetMT2(double testmass=0 , bool massive=false, int met=1);
  Double_t GetMT2Leading(double testmass=0, bool massive=true, int PFJID=1, int met=1);
  Double_t GetMT2Hemi(double testmass=0, bool massive=false, int PFJID=1, 
		  double minJPt=20, int hemi_association=3, int met=1);
  Double_t CalcMT2(double testmass, bool massive, 
		  TLorentzVector visible1, TLorentzVector visible2, TLorentzVector MET );
  Double_t GetMCT(bool massive=false, int met=1);

  MT2Misc misc;
  Int_t   NJets;
  Int_t   NJetsIDLoose;
  Int_t   NJetsIDMedium;
  Int_t   NJetsIDTight;
  Int_t   NEles;
  Int_t   NElesLoose;
  Int_t   NMuons;
  Int_t   NMuonsLoose;
  MT2Jet         jet[m_jetSize];
  MT2Elec        ele[m_eleSize];
  MT2Muon        muo[m_muoSize];
  TLorentzVector pfmet[2];
  TLorentzVector MPT[2];
  TLorentzVector MHT[2];
  TLorentzVector MHTloose[2];
  TLorentzVector pseudoJets[2];

  
  ClassDef(MT2tree, 5)
};

#endif
