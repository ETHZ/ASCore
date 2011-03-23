#ifndef MT2tree_hh
#define MT2tree_hh

#include "TObject.h"
#include "TLorentzVector.h"
#include "TVector3.h"

enum {m_jetSize = 40, m_eleSize = 15, m_muoSize = 15, m_genleptSize=30, m_hemiSize=10};

// MT2Misc ----------------------------------
class MT2Misc : public TObject {

public:
  MT2Misc();
  virtual ~MT2Misc();

  void Reset();
  
  Bool_t   HBHENoiseFlag;
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

  ClassDef(MT2Misc, 10)
};


// MT2Znunu --------------------------------
class MT2Znunu : public TObject {

public:
	MT2Znunu();
	virtual ~MT2Znunu();
	void Reset();

  	Int_t    NJetsToRemoveMuo;
  	Int_t    NJetsToRemoveEle;
	Int_t    NJetsIDLoose_matched;
	Double_t GenZee_mll;
	Double_t GenZee_mll_acc; 
	Double_t GenZmumu_mll;
	Double_t GenZmumu_mll_acc;
	Double_t GenZnunu_e_mll;
	Double_t GenZnunu_e_mll_acc;
	Double_t GenZnunu_mu_mll;
	Double_t GenZnunu_mu_mll_acc;
	Double_t GenZnunu_tau_mll;
	Double_t GenZnunu_tau_mll_acc;
	Double_t RecoOSee_mll;
	Double_t RecoOSmumu_mll;
	Double_t HTmatched;
	Double_t METplusLeptsPt;
	Double_t METplusLeptsPtReco;
	Double_t MinMetplusLeptJetDPhi;
	Double_t MinMetplusLeptJetDPhiReco;
	Double_t PassJetID_matched;
	Double_t Jet1Pass_matched;
	Double_t Jet0Pass_matched;
	Double_t Vectorsumpt_matched;

	ClassDef(MT2Znunu, 2);
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

  ClassDef(MT2Hemi, 4)
};


// MT2Elec ----------------------------------
class MT2Elec : public TObject {

public:
  MT2Elec();
  virtual ~MT2Elec();

  void Reset();
  void SetLV(const TLorentzVector v);

  TLorentzVector lv;

  Double_t MT;
  Double_t Iso;
  Int_t    Charge;

  ClassDef(MT2Elec, 5)
};

// MT2Muon ----------------------------------
class MT2Muon : public TObject {

public:
  MT2Muon();
  virtual ~MT2Muon();

  void Reset();
  void SetLV(const TLorentzVector v);

  TLorentzVector lv;

  Double_t MT;
  Double_t Iso;
  Int_t    Charge;

  ClassDef(MT2Muon, 5)
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
  Double_t       MT;


  ClassDef(MT2GenLept, 3)
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
  void SetNMuons        (int n);
  
  // My functions here
  // NJets
  Int_t    GetNjets (double minJPt=20, double maxJEta=5., int PFJID=0);  // PFJETID not depends on pt and eta
  Int_t    GetJetIndex(int ijet=0, int PFJID=1, double minJPt=20, double maxJEta=2.4);
  Int_t    GetJetIndexByEta(int ijet=0, int PFJID=1, double minJPt=20, double maxJEta=2.4);
  Int_t    GetNBtags (int algo=3, double value=2., double minJPt=20, double maxJEta=2.4, int PFJID=1);  // algo - 0:TCHE, 1:TCHP, 2:SSVHE, 3:SSVHP
  Double_t JetPt      (int ijet=0, int PFJID=1, double minJPt=20, double maxJEta=2.4);

  // HT, MHT, ...
  Double_t GetHT         (int PFJID=0, double minJPt=50, double maxJEta=2.4);
  TLorentzVector GetMHTlv(int PFJID=1, double minJPt=20, double maxJEta=2.4);
  Double_t GetMHT        (int PFJID=1, double minJPt=20, double maxJEta=2.4);
  Double_t GetMHTPhi     (int PFJID=1, double minJPt=20, double maxJEta=2.4);
  Double_t GetMHTminusMET(int PFJID=1, double minJPt=20, double maxJEta=2.4);
  // dPhi and friends
  Bool_t   PassJetID(double minJPt=50, double maxJEta=5.0, int PFJID=1);
  Double_t JetsDPhi(int j1=1, int j2=0, int PFJID=1);
  Double_t JetsInvMass(int j1=0, int j2=1);
  Double_t MetJetDPhi(int ijet = 0, int PFJID=1, int met=1);
  Double_t GetMinR12R21      (int PFJID=1, double minJPt=20, double maxJEta=6., int met=1); // electrons and muons not considered for minDPhi
  Double_t MinMetJetDPhi     (int PFJID=1, double minJPt=20, double maxJEta=6., int met=1); // electrons and muons not considered for minDPhi
  Int_t    MinMetJetDPhiIndex(int PFJID=1, double minJPt=20, double maxJEta=6., int met=1); // electrons and muons not considered for minDPhi
  Double_t MaxMetJetDPhi     (int PFJID=1, double minJPt=20, double maxJEta=6., int met=1); // electrons and muons not considered for minDPhi
  Int_t    MaxMetJetDPhiIndex(int PFJID=1, double minJPt=20, double maxJEta=6., int met=1); // electrons and muons not considered for minDPhi
  Double_t GetPseudoJetsdPhi(int hemi_seed=2, int hemi_association=3, 
		                   int PFJID=1, double minJPt=20, double maxJEta=2.4);
  Double_t GetPseudoJetsdPhiMinDHT(int PFJID=0, double minJPt=20, double maxJEta=2.4);
  Double_t GetPseudoJetMetDPhi(int hemi_index=0, int pj=1, int whichmet=1, double met=30);

  // MT2 & friends
  Int_t    JetIsInHemi(int jindex=0, int hemi_seed=2, int hemi_association=3, float MaxDR=0);
  Double_t GetMT2Leading(double testmass=0, bool massive=true, int PFJID=1, int met=1);
  Double_t GetMT2Hemi(double testmass=0, bool massive=false, int PFJID=1, 
		  double minJPt=20, double maxJEta=2.4, int hemi_association=3, int met=1);
  Double_t GetMT2HemiMinDHT(double testmass=0, bool massive=false, int PFJID=0, double minJPt=20, double maxJEta=2.4, int met=1);
  Double_t GetMT2HemiNoISR(bool massive = false, int hemi_seed=4, int hemi_association=2, float MaxDR=0, int met=1);
  Double_t CalcMT2(double testmass, bool massive, 
		  TLorentzVector visible1, TLorentzVector visible2, TLorentzVector MET );
  Double_t SimpleMT2(bool pseudo=true);
  Double_t GetMCT(bool massive=false, int met=1);
  Double_t GetSqrtS(double testmass=0, bool massive=true,int PFJID=1, double minJPt=20, double maxJEta=2.4,int met=1);
  Double_t GetMaxHemiMass(int hemi_index=0);

  // Leptons
  Double_t GenOSDiLeptonInvMass(unsigned int pid=11, unsigned int mother=23, double pt=10, double eta=2.4);
  Double_t GetDiLeptonInvMass(int same_sign=0, int same_flavour=1, int flavour=0, double pt=10, bool exclDiLept=false);
  Bool_t   IsDiLeptonMll(int same_sign=0, int same_flavour=1, int flavour=0, double pt=10, bool exclDiLept=false, double lower_mass=71, double upper_mass=111);
  Bool_t   IsGenOSDiLepton(unsigned int pid=11, unsigned int mother=23, double pt=10, double eta=2.4, double lower_mass=71, double upper_mass=111);
  TLorentzVector GetMETPlusLeptsLV(int OSDiLeptFromZ =1);
  Double_t GetMETPlusLepts(int OSDiLeptFromZ =1);
  Double_t GetMETPlusGenLepts(int met, int RemoveOSSFDiLepts=0, int require_cuts=1, unsigned int pid=11, 
		              unsigned int mother=23, double pt=10, double eta=2.4, double lower_mass=71, double upper_mass=111);
  Double_t GetDiLeptonPt(int same_sign=0, int same_flavour=1, int flavour=0, double pt=10, double lower_mass=71, double upper_mass=111);
  Double_t GetGenLeptPt(int which, int pid, int mother, double pt, double eta);
  Double_t GetGenLeptEta(int which, int pid, int mother, double pt, double eta);
  Int_t    GetGenLeptIndex(int which, int pid, int mother, double pt, double eta);
  Bool_t   GenLeptFromW(int pid, double pt, double eta);

  Int_t     NJets;
  Int_t     NJetsIDLoose;
  Int_t     NJetsIDMedium;
  Int_t     NJetsIDTight;
  Int_t     NEles;
  Int_t     NElesLoose;
  Int_t     NMuons;
  Int_t     NMuonsLoose;
  Int_t     NGenLepts;

  MT2Misc        misc;
  MT2Znunu       Znunu;
  MT2Jet         jet[m_jetSize];
  MT2Hemi        hemi[m_hemiSize];
  MT2Elec        ele[m_eleSize];
  MT2Muon        muo[m_muoSize];
  MT2GenLept     genlept[m_genleptSize];
  TLorentzVector pfmet[2];
  TLorentzVector genmet[2];
  TLorentzVector MPT[2];
  TLorentzVector MHT[2];
  TLorentzVector MHTloose[2];

  
  ClassDef(MT2tree, 15)
};

#endif
