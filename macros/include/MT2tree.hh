#ifndef MT2tree_hh
#define MT2tree_hh

#include "TObject.h"
#include "TLorentzVector.h"
#include "TVector3.h"

enum {m_jetSize = 40, m_eleSize = 15, m_muoSize = 15, m_genleptSize=30, m_hemiSize=4};

// MT2Misc ----------------------------------
class MT2Misc : public TObject {

public:
  MT2Misc();
  virtual ~MT2Misc();

  void Reset();
  
  Int_t    Run;
  Int_t    Event;
  Int_t    LumiSection;
  Int_t    NVertices;  
  Int_t    LeptConfig;
  Int_t    NJetsEta5Pt20;
  Int_t    Jet0Pass;
  Int_t    Jet1Pass;
  Int_t    PassJetID;
  Int_t    EcalDeadCellBEFlag;
  Int_t    NECALGapClusters;
  Double_t EcalGapClusterSize[50];
  Double_t EcalGapBE[50];
  Double_t MT2;
  Double_t MT2leading;
  Double_t MT2noISR;
  Double_t MCT;
  Double_t AlphaT;
  Double_t MET;
  Double_t METPhi;
  Double_t Vectorsumpt;
  Double_t PFMETsign;
  Double_t DPhiMhtMpt;
  Double_t MinMetJetDPhi;
  Double_t HT;
  Bool_t   HBHENoiseFlag;

  ClassDef(MT2Misc, 6)
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

  ClassDef(MT2Jet, 5)
};

// MT2Hemi ---------------------------
class MT2Hemi : public TObject {

public:
  MT2Hemi();
  virtual ~MT2Hemi();

  void Reset();
  Int_t          seed_method;
  Int_t          assoc_method;
  Double_t       MT2;
  Double_t       MCT;
  Double_t       AlphaT;
  Double_t       minDHT;
  Double_t       maxDR;
  Double_t       dPhi;

  Int_t          jindices1  [m_jetSize];
  Int_t          jindices2  [m_jetSize];
  Int_t          eleindices1[m_eleSize];
  Int_t          eleindices2[m_eleSize];
  Int_t          muoindices1[m_muoSize];
  Int_t          muoindices2[m_muoSize];
  TLorentzVector lv1;
  TLorentzVector lv2;
  TLorentzVector UTM;

  ClassDef(MT2Hemi, 3)
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


// MT2GenLept ----------------------------------
class MT2GenLept : public TObject {

public:
  MT2GenLept();
  virtual ~MT2GenLept();

  void Reset();

  TLorentzVector lv;
  Int_t          ID;
  Int_t          MID;
  Int_t          MStatus;
  Int_t          GMID;
  Int_t          GMStatus;


  ClassDef(MT2GenLept, 1)
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
  Int_t    GetNjets(double minJPt=20, double maxJEta=5., int PFJID=0);  // PFJETID not depends on pt and eta
  Int_t    GetJetIndex(int ijet=0, int PFJID=1);
  Double_t JetPt      (int ijet=0, int PFJID=1);

  // dPhi and friends
  Bool_t   PassJetID(double minJPt=50, double maxJEta=5.0, int PFJID=1);
  Double_t JetsDPhi(int j1=1, int j2=0, int PFJID=1);
  Double_t JetsInvMass(int j1=0, int j2=1);
  Double_t MetJetDPhi(int ijet = 0, int PFJID=1, int met=1);
  Double_t GetMinR12R21      (int PFJID=1, double minJPt=20, double maxJEta=6., int met=1); // electrons and muons not considered for minDPhi
  Double_t MinMetJetDPhi     (int PFJID=1, double minJPt=20, double maxJEta=6., int met=1); // electrons and muons not considered for minDPhi
  Int_t    MinMetJetDPhiIndex(int PFJID=1, double minJPt=20, double maxJEta=6., int met=1); // electrons and muons not considered for minDPhi
  // MT2 & friends
  Int_t    JetIsInHemi(int jindex=0, int hemi_seed=2, int hemi_association=3, float MaxDR=0);
  Double_t GetMT2Leading(double testmass=0, bool massive=true, int PFJID=1, int met=1);
  Double_t GetMT2Hemi(double testmass=0, bool massive=false, int PFJID=1, 
		  double minJPt=20, int hemi_association=3, int met=1);
  Double_t GetMT2HemiNoISR(bool massive = false, int hemi_seed=4, int hemi_association=2, float MaxDR=0, int met=1);
  Double_t CalcMT2(double testmass, bool massive, 
		  TLorentzVector visible1, TLorentzVector visible2, TLorentzVector MET );
  Double_t GetMCT(bool massive=false, int met=1);

  Int_t   NJets;
  Int_t   NJetsIDLoose;
  Int_t   NJetsIDMedium;
  Int_t   NJetsIDTight;
  Int_t   NEles;
  Int_t   NElesLoose;
  Int_t   NMuons;
  Int_t   NMuonsLoose;
  MT2Misc        misc;
  MT2Jet         jet[m_jetSize];
  MT2Hemi        hemi[m_hemiSize];
  MT2Elec        ele[m_eleSize];
  MT2Muon        muo[m_muoSize];
  MT2GenLept     genlept[m_genleptSize];
  TLorentzVector pfmet[2];
  TLorentzVector MPT[2];
  TLorentzVector MHT[2];
  TLorentzVector MHTloose[2];

  
  ClassDef(MT2tree, 8)
};

#endif
