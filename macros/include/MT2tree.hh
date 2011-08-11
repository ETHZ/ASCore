#ifndef MT2tree_hh
#define MT2tree_hh

#include "TObject.h"
#include "TLorentzVector.h"
#include "TVector3.h"

enum {m_jetSize = 40, m_genjetSize = 20, m_eleSize = 5, m_muoSize = 5, m_genleptSize=30, m_hemiSize=4};
// ele and muo are 10 previous to MT2_V00-06-00
// hemi size was 10 previos to MT2_V00-05-00

// MT2Misc ----------------------------------
class MT2Misc : public TObject {

public:
  MT2Misc();
  virtual ~MT2Misc();

  void Reset();
  
  Bool_t   HBHENoiseFlag;
  Bool_t   CrazyHCAL;
  Bool_t   NegativeJEC;
  Bool_t   isData;
  Bool_t   BadEcalTP;
  Bool_t   BadEcalBE;
  Bool_t   CSCTightHaloID;
  Int_t    Run;
  Int_t    Event;
  Int_t    LumiSection;
  Int_t    LeptConfig;
  Int_t    Jet0Pass;
  Int_t    Jet1Pass;
  Int_t    PassJetID;
  Int_t    PassJetID20;
  Double_t MT2;
  Double_t MT2all;
  Double_t MT2leading;
  Double_t MT2noISR;
  Double_t MCT;
  Double_t AlphaT;
  Double_t MET;
  Double_t METPhi;
  Double_t LeadingJPt;
  Double_t SecondJPt;
  Double_t Vectorsumpt;
  Double_t VectorsumptAll;
  Double_t PFMETsign;
  Double_t DPhiMhtMpt;
  Double_t MinMetJetDPhi;
  Double_t HT;
  Double_t caloHT40;    
  Double_t caloHT40_ID;
  Double_t caloHT50;  
  Double_t caloHT50_ID;
  Double_t caloMHT30;  
  Double_t caloMHT30_ID;
  Double_t caloMHT40;  
  Double_t caloMHT40_ID;
  Double_t TrackingFailure;
  Double_t TrackingFailurePVtx;
  
  ClassDef(MT2Misc, 22)
};

// ----------------------------------------
class MT2PileUp : public TObject {

public:
	MT2PileUp();
	virtual ~MT2PileUp();
	void Reset();

	Int_t    PUnumInt;
	Int_t    PUnumIntEarly;
	Int_t    PUnumIntLate;
	Int_t    isS3;
	Double_t PtHat;
	Double_t Weight;
  	Int_t    NVertices;  // good reco vertices
	Double_t Rho;

	ClassDef(MT2PileUp, 3);
};

// --------------------------------
class MT2Trigger : public TObject {

public:
	MT2Trigger();
	virtual ~MT2Trigger();
	void Reset();

	// HT with DPhi
	Bool_t HLT_HT500_JetPt60_DPhi2p94_v1;
	Bool_t HLT_HT550_JetPt60_DPhi2p94_v1;

	// HT
	Bool_t HLT_HT150_v2;
	Bool_t HLT_HT150_v3;
	Bool_t HLT_HT160_v2;
	Bool_t HLT_HT200_v2;
	Bool_t HLT_HT200_v3;
	Bool_t HLT_HT240_v2;
	Bool_t HLT_HT250_v2;
	Bool_t HLT_HT250_v3;
	Bool_t HLT_HT260_v2;
	Bool_t HLT_HT300_v2;
	Bool_t HLT_HT300_v3;
	Bool_t HLT_HT300_v4;
	Bool_t HLT_HT300_v5;
	Bool_t HLT_HT350_v2;
	Bool_t HLT_HT350_v3;
	Bool_t HLT_HT350_v4;
	Bool_t HLT_HT360_v2;
	Bool_t HLT_HT400_v2;
	Bool_t HLT_HT400_v3;
	Bool_t HLT_HT400_v4;
	Bool_t HLT_HT400_v5;
	Bool_t HLT_HT400_v6;
	Bool_t HLT_HT400_v7;
	Bool_t HLT_HT400_v8;
	Bool_t HLT_HT440_v2;
	Bool_t HLT_HT450_v2;
	Bool_t HLT_HT450_v3;
	Bool_t HLT_HT450_v4;
	Bool_t HLT_HT450_v5;
	Bool_t HLT_HT450_v6;
	Bool_t HLT_HT450_v7;
	Bool_t HLT_HT450_v8;
	Bool_t HLT_HT500_v2;
	Bool_t HLT_HT500_v3;
	Bool_t HLT_HT500_v4;
	Bool_t HLT_HT500_v5;
	Bool_t HLT_HT500_v6;
	Bool_t HLT_HT500_v7;
	Bool_t HLT_HT500_v8;
	Bool_t HLT_HT550_v2;
	Bool_t HLT_HT550_v3;
	Bool_t HLT_HT550_v4;
	Bool_t HLT_HT550_v5;
	Bool_t HLT_HT550_v6;
	Bool_t HLT_HT550_v7;
	Bool_t HLT_HT550_v8;
	Bool_t HLT_HT600_v1;
	// HT_MHT
	Bool_t HLT_HT250_MHT60_v2;
	Bool_t HLT_HT250_MHT60_v3;
	Bool_t HLT_HT250_MHT60_v4;
	Bool_t HLT_HT250_MHT60_v5;
	Bool_t HLT_HT250_MHT60_v6;
	Bool_t HLT_HT250_MHT70_v1;
	Bool_t HLT_HT250_MHT70_v2;
	Bool_t HLT_HT250_MHT70_v3;
	Bool_t HLT_HT250_MHT70_v4;
	Bool_t HLT_HT250_MHT90_v1;
	Bool_t HLT_HT250_MHT90_v2;
	Bool_t HLT_HT260_MHT60_v2;
	Bool_t HLT_HT300_MHT75_v4;
	Bool_t HLT_HT300_MHT75_v5;
	Bool_t HLT_HT300_MHT75_v7;
	Bool_t HLT_HT300_MHT75_v8;
	Bool_t HLT_HT300_MHT80_v1;
	Bool_t HLT_HT300_MHT80_v2;
	Bool_t HLT_HT300_MHT90_v1;
	Bool_t HLT_HT300_MHT90_v2;
	Bool_t HLT_HT350_MHT70_v1;
	Bool_t HLT_HT350_MHT70_v2;
	Bool_t HLT_HT350_MHT80_v1;
	Bool_t HLT_HT350_MHT80_v2;
	// Muons
	Bool_t HLT_DoubleMu3_HT160_v2;
	Bool_t HLT_Mu8_Jet40_v2;
	Bool_t HLT_DoubleMu3_v3;
	
	ClassDef(MT2Trigger, 10);
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
	Double_t caloMHT30_matched;
	Double_t caloMHT30ID_matched;
	Double_t caloMHT30_matchedReco;
	Double_t caloMHT30ID_matchedReco;
	Double_t caloHT50_matched;
	Double_t caloHT50ID_matched;
	Double_t caloHT50_matchedReco;
	Double_t caloHT50ID_matchedReco;
	Double_t HTmatched;
	Double_t METplusLeptsPt;
	Double_t METplusLeptsPtReco;
	Double_t MinMetplusLeptJetDPhi;
	Double_t MinMetplusLeptJetDPhiReco;
	Double_t PassJetID_matched;
	Double_t Jet1Pass_matched;
	Double_t Jet0Pass_matched;
	Double_t LeadingJPt_matched;
	Double_t SecondJPt_matched;
	Double_t Vectorsumpt_matched;

	ClassDef(MT2Znunu, 6);
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

  Bool_t isPATPFIDLoose;
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

  Double_t Scale;          // scale factor from JE correction
  Double_t L1FastJetScale; // correction factor from raw to L1FastJetcorrected
  Double_t Area;

  Int_t    Flavour;   // JetFlavour for MC
  
  Bool_t   isTau;      // has to be *ALWAYS FALSE* starting from V02-01-01
  Bool_t   isTauMatch; // tells you if pf-jet is matched to a tau
  Double_t TauDR;
  Double_t TauDPt;
  Int_t    NTauMatch;

  ClassDef(MT2Jet, 10)
};

// MT2GenJet -------------------------
class MT2GenJet : public TObject {

public:
  MT2GenJet();
  virtual ~MT2GenJet();

  void Reset();
  TLorentzVector lv;
  
  Int_t    JetMatchIndex;
  Double_t DeltaR;

  ClassDef(MT2GenJet, 1)
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
  Int_t    ID95;
  Int_t    ID90;
  Int_t    CAJ_n90;
  Int_t    CAJ_n90Hits;

  ClassDef(MT2Elec, 7)
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
  Int_t    NMatches;
  Double_t PtErr;

  ClassDef(MT2Muon, 6)
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
  Double_t       CAJ_n90;
  Double_t       CAJ_n90Hits;


  ClassDef(MT2GenLept, 3)
};

// MT2tree ----------------------------------
class MT2tree : public TObject {

public:
  MT2tree();
  virtual ~MT2tree();

  void Reset();

  void SetNJets         (int n);
  void SetNGenJets      (int n);
  void SetNJetsIDLoose  (int n);
  void SetNJetsIDMedium (int n);
  void SetNJetsIDTight  (int n);
  void SetNJetsAcc      (int n);
  void SetNBJets        (int n);
  void SetNEles         (int n);
  void SetNMuons        (int n);
  void SetNTaus         (int n);
  
  // My functions here
  // NJets
  Int_t    GetNjets   (double minJPt=20, double maxJEta=5., int PFJID=0);  // PFJETID not depends on pt and eta
  Int_t    GetJetIndex(int ijet=0, int PFJID=0, double minJPt=20, double maxJEta=2.4);
  Int_t    GetJetIndexByEta(int ijet=0, int PFJID=0, double minJPt=20, double maxJEta=2.4);
  Int_t    GetNBtags  (int algo=3, double value=2., double minJPt=20, double maxJEta=2.4, int PFJID=0);  // algo - 0:TCHE, 1:TCHP, 2:SSVHE, 3:SSVHP
  Double_t JetPt      (int ijet=0, int PFJID=0, double minJPt=20, double maxJEta=2.4);

  Int_t GetQcdType(bool matchOnPhi=true, int PFID=0);

  // HT, MHT, ...
  Double_t GetHT         (int PFJID=0, double minJPt=50, double maxJEta=2.4);
  TLorentzVector GetMHTlv(int PFJID=0, double minJPt=20, double maxJEta=2.4, bool inclLepts=false);
  Double_t GetMHT        (int PFJID=0, double minJPt=20, double maxJEta=2.4, bool inclLepts=false);
  Double_t GetMHTPhi     (int PFJID=0, double minJPt=20, double maxJEta=2.4, bool inclLepts=false);
  Double_t GetMHTminusMET(int PFJID=0, double minJPt=20, double maxJEta=2.4, bool inclLepts=false);
  // dPhi and friends
  Bool_t   PassJetID(double minJPt=50, double maxJEta=5.0, int PFJID=1);
  Double_t JetsDPhi(int j1=1, int j2=0, int PFJID=0);
  Double_t JetsInvMass(int j1=0, int j2=1);
  Double_t MetJetDPhi(int ijet = 0, int PFJID=0, int met=1);
  Bool_t   PassMinMetJetDPhi03();
  Double_t GetMinR12R21      (int PFJID=0, double minJPt=20, double maxJEta=6., int met=1); // electrons and muons not considered for minDPhi
  Double_t MinMetJetDPhi     (int PFJID=0, double minJPt=20, double maxJEta=6., int met=1); // electrons and muons not considered for minDPhi
  Int_t    MinMetJetDPhiIndex(int PFJID=0, double minJPt=20, double maxJEta=6., int met=1); // electrons and muons not considered for minDPhi
  Double_t MaxMetJetDPhi     (int PFJID=0, double minJPt=20, double maxJEta=6., int met=1); // electrons and muons not considered for minDPhi
  Int_t    MaxMetJetDPhiIndex(int PFJID=0, double minJPt=20, double maxJEta=6., int met=1); // electrons and muons not considered for minDPhi
  Double_t MinMetJetDPhiL2L3 ();
  Double_t GetPseudoJetsdPhi(int hemi_seed=2, int hemi_association=3, int PFJID=0, double minJPt=20, double maxJEta=2.4);
  Double_t GetPseudoJetsdPhiMinDHT(int PFJID=0, double minJPt=20, double maxJEta=2.4);
  Double_t GetPseudoJetMetDPhi(int hemi_index=0, int pj=1, int whichmet=1, double met=30);
  Double_t PseudoJetMetDPhi();

  // MT2 & friends
  Double_t GetMT2Leading(double testmass=0, bool massive=true, int PFJID=0, int met=1);
  Double_t GetMT2Hemi(double testmass=0, bool massive=false, int PFJID=0, double minJPt=20, double maxJEta=2.4, int hemi_association=3, int met=1);
  Double_t GetMT2HemiMinDHT(double testmass=0, bool massive=false, int PFJID=0, double minJPt=20, double maxJEta=2.4, int met=1);
  Double_t GetMT2HemiNoISR(bool massive = false, int hemi_seed=4, int hemi_association=2, float MaxDR=0, int met=1);
  Double_t CalcMT2(double testmass, bool massive, 
		  TLorentzVector visible1, TLorentzVector visible2, TLorentzVector MET );
  Double_t TopDileptonMT2(bool massive=true, bool moreJets=false);
  Double_t WDileptonMT2(bool massive=true);
  Double_t SimpleMT2(bool pseudo=true);
  Double_t GetMCT(bool massive=false, int met=1);
  Double_t GetSqrtS(double testmass=0, bool massive=true,int PFJID=0, double minJPt=20, double maxJEta=2.4,int met=1);
  Double_t GetMaxHemiMass(int hemi_index=0);
  Double_t PseudoJetPtRatio();
  Double_t HemiMassTop();
  Double_t RemoveAndRecalcMT2(int option=0, float prob=0.01, bool hiDphi=true, float dphi=0.0, bool selection=false);

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
  Bool_t   GenLeptFromW(int pid, double pt, double eta, bool includeTaus);
  Double_t GetLeptPt(int index);
  Double_t ElClosestJet();
  Int_t    WDecayMode();
  Int_t    TopDecayMode();
  Bool_t   TopDecayModeResult(Int_t nlepts);
  Bool_t   SLTopAccept(double pt, double eta);
  Double_t SLTopEta(double pt);
  Double_t LeptJetDR(int pid, int index, bool bjet, int ID);

  //Bosons
  Double_t GetGenVPt(int pid);

  // LostLepton estimate
  Bool_t   LostLeptonChanges(); 

  Int_t     NJets;
  Int_t     NGenJets;
  Int_t     NJetsIDLoose;
  Int_t     NJetsIDMedium;
  Int_t     NJetsIDTight;
  Int_t     NJetsAcc;
  Int_t     NBJets;
  Int_t     NEles;
  Int_t     NMuons;
  Int_t     NTaus;
  Int_t     NGenLepts;

  MT2Misc        misc;
  MT2Znunu       Znunu;
  MT2PileUp      pileUp;
  MT2Trigger     trigger;
  MT2Jet         jet[m_jetSize];
  MT2GenJet      genjet[m_genjetSize];
  MT2Hemi        hemi[m_hemiSize];
  MT2Elec        ele[m_eleSize];
  MT2Muon        muo[m_muoSize];
  MT2GenLept     genlept[m_genleptSize];
  TLorentzVector pfmet[2];
  TLorentzVector genmet[2];
  TLorentzVector MPT[2];
  TLorentzVector MHT[2];
  TLorentzVector MHTloose[2];

  
  ClassDef(MT2tree, 18)
};

#endif
