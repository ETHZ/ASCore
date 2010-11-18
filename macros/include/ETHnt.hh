#ifndef ETHnt_hh
#define ETHnt_hh

#include "TObject.h"
#include "TLorentzVector.h"

enum {m_jetSize = 20, m_eleSize = 10, m_muoSize = 10};

class ETHMisc : public TObject {

public:
  ETHMisc();
  virtual ~ETHMisc();

  void Reset();
  
  Int_t    Run;
  Int_t    Event;
  Int_t    LumiSection;
  Float_t  Weight;
  Int_t    LeptConfig;
  Bool_t   IsCleanMultiJetEvent;
  Bool_t   IsCleanJetEvent;
  Int_t    R12R21;
  Int_t    NJetsPt50Eta25;
  Double_t PseudoJetMT2;
  Double_t PseudoJetMCT;
  Double_t PseudoJet1Pt;
  Double_t PseudoJet2Pt;
  Double_t PseudojetAlphaT;
  Double_t Vectorsumpt;
  Double_t PFMET;
  Double_t PFMETphi;
  Double_t PFMETsign;
  Double_t LeadingJetEta;
  Double_t DPhiJ1Met;
  Double_t DPhiJ2Met;
  Double_t PseudoJetMT2AxisdPhi;
  Double_t R1221min;
  Double_t MPT_sel;
  Double_t MPT;
  Double_t MHT;
  Double_t HT;
  Double_t DPhiMhtMpt;

  ClassDef(ETHMisc, 1)
};

class ETHJet : public TObject {

public:
  ETHJet();
  virtual ~ETHJet();

  void Reset();
  void SetLV(const TLorentzVector v);

  TLorentzVector lv;

  Double_t bTagProbTCHE ;
  Double_t bTagProbTCHP ;
  Double_t bTagProbSSVHE;
  Double_t bTagProbSSVHP;
  Int_t    inHemisphere;

  ClassDef(ETHJet, 1)
};

class ETHnt : public TObject {

public:
  ETHnt();
  virtual ~ETHnt();

  void Reset();

  void SetNJets (int n);
  void SetNEles (int n);
  void SetNMuons(int n);
  
  // My functions here
  Double_t PseudoJetDPhi (); // is this = misc.PseudoJetMT2AxisdPhi?
  Double_t PseudoJetAngle(); // or is this one?

  ETHMisc misc;
  Int_t   NJets;
  Int_t   NEles;
  Int_t   NMuons;
  ETHJet         jet[m_jetSize];
  TLorentzVector ele[m_eleSize];
  TLorentzVector muo[m_muoSize];
  TLorentzVector pfmet[2];
  TLorentzVector pseudoJets[2];

  
  ClassDef(ETHnt, 1)
};

#endif
