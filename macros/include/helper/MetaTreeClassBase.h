//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Sep 17 20:23:44 2012 by ROOT version 5.32/00
// from TTree Analysis/AnalysisTree
// found on file: tautest.root
//////////////////////////////////////////////////////////

#ifndef MetaTreeClassBase_h
#define MetaTreeClassBase_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class MetaTreeClassBase {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           Run;
   Int_t           Event;
   Int_t           LumiSec;
   Float_t         m0;
   Float_t         m12;
   Int_t           process;
   Float_t         mGlu;
   Float_t         mLSP;
   Int_t           isTChiSlepSnu;
   Int_t           isRightHanded;
   Int_t           HLT_MU8;
   Int_t           HLT_MU8_PS;
   Int_t           HLT_MU17;
   Int_t           HLT_MU17_PS;
   Int_t           HLT_ELE17_TIGHT;
   Int_t           HLT_ELE17_TIGHT_PS;
   Int_t           HLT_ELE17_JET30_TIGHT;
   Int_t           HLT_ELE17_JET30_TIGHT_PS;
   Int_t           HLT_ELE8_TIGHT;
   Int_t           HLT_ELE8_TIGHT_PS;
   Int_t           HLT_ELE8_JET30_TIGHT;
   Int_t           HLT_ELE8_JET30_TIGHT_PS;
   Int_t           HLT_MU17_MU8;
   Int_t           HLT_MU17_MU8_PS;
   Int_t           HLT_MU17_TKMU8;
   Int_t           HLT_MU17_TKMU8_PS;
   Int_t           HLT_ELE17_ELE8_TIGHT;
   Int_t           HLT_ELE17_ELE8_TIGHT_PS;
   Int_t           HLT_MU8_ELE17_TIGHT;
   Int_t           HLT_MU8_ELE17_TIGHT_PS;
   Int_t           HLT_MU17_ELE8_TIGHT;
   Int_t           HLT_MU17_ELE8_TIGHT_PS;
   Float_t         Rho;
   Int_t           NVrtx;
   Float_t         PUWeight;
   Int_t           NMus;
   Int_t           IsSignalMuon[5];   //[NMus]
   Float_t         MuPt[5];   //[NMus]
   Float_t         MuEta[5];   //[NMus]
   Float_t         MuPhi[5];   //[NMus]
   Int_t           MuCharge[5];   //[NMus]
   Float_t         MuPFIso[5];   //[NMus]
   Float_t         MuPFChIso[5];   //[NMus]
   Float_t         MuPFNeIso[5];   //[NMus]
   Float_t         MuRadIso[5];   //[NMus]
   Float_t         MuD0[5];   //[NMus]
   Float_t         MuDz[5];   //[NMus]
   Float_t         MuEMVetoEt[5];   //[NMus]
   Float_t         MuHadVetoEt[5];   //[NMus]
   Int_t           MuPassesTightID[5];   //[NMus]
   Float_t         MuPtE[5];   //[NMus]
   Int_t           MuGenID[5];   //[NMus]
   Int_t           MuGenMID[5];   //[NMus]
   Int_t           MuGenGMID[5];   //[NMus]
   Int_t           MuGenType[5];   //[NMus]
   Int_t           MuGenMType[5];   //[NMus]
   Int_t           MuGenGMType[5];   //[NMus]
   Float_t         MuMT[5];   //[NMus]
   Int_t           NEls;
   Int_t           IsSignalElectron[5];   //[NEls]
   Int_t           ElCharge[5];   //[NEls]
   Int_t           ElChIsCons[5];   //[NEls]
   Float_t         ElPt[5];   //[NEls]
   Float_t         ElEta[5];   //[NEls]
   Float_t         ElPhi[5];   //[NEls]
   Float_t         ElD0[5];   //[NEls]
   Float_t         ElD0Err[5];   //[NEls]
   Float_t         ElDz[5];   //[NEls]
   Float_t         ElDzErr[5];   //[NEls]
   Float_t         ElPFIso[5];   //[NEls]
   Float_t         ElPFChIso[5];   //[NEls]
   Float_t         ElPFNeIso[5];   //[NEls]
   Float_t         ElRadIso[5];   //[NEls]
   Float_t         ElEcalRecHitSumEt[5];   //[NEls]
   Float_t         ElHcalTowerSumEt[5];   //[NEls]
   Float_t         ElTkSumPt[5];   //[NEls]
   Float_t         ElDPhi[5];   //[NEls]
   Float_t         ElDEta[5];   //[NEls]
   Float_t         ElSigmaIetaIeta[5];   //[NEls]
   Float_t         ElHoverE[5];   //[NEls]
   Float_t         ElEPthing[5];   //[NEls]
   Int_t           ElIsGoodElId_LooseWP[5];   //[NEls]
   Int_t           ElIsGoodElId_MediumWP[5];   //[NEls]
   Int_t           ElIsGoodTriggerEl[5];   //[NEls]
   Int_t           ElGenID[5];   //[NEls]
   Int_t           ElGenMID[5];   //[NEls]
   Int_t           ElGenGMID[5];   //[NEls]
   Int_t           ElGenType[5];   //[NEls]
   Int_t           ElGenMType[5];   //[NEls]
   Int_t           ElGenGMType[5];   //[NEls]
   Float_t         ElMT[5];   //[NEls]
   Int_t           NTaus;
   Int_t           TauCharge[5];   //[NTaus]
   Float_t         TauPt[5];   //[NTaus]
   Float_t         TauEta[5];   //[NTaus]
   Float_t         TauPhi[5];   //[NTaus]
   Float_t         TauMVAElRej[5];   //[NTaus]
   Float_t         TauTightMuRej[5];   //[NTaus]
   Float_t         TauLCombIsoDB[5];   //[NTaus]
   Float_t         pfMET;
   Float_t         pfMETPhi;
   Float_t         pfMETType1;
   Float_t         pfMETType1Phi;
   Int_t           NJets;
   Float_t         JetPt[25];   //[NJets]
   Float_t         JetEta[25];   //[NJets]
   Float_t         JetPhi[25];   //[NJets]
   Float_t         JetEnergy[25];   //[NJets]
   Float_t         JetCSVBTag[25];   //[NJets]
   Float_t         JetProbBTag[25];   //[NJets]
   Float_t         JetArea[25];   //[NJets]
   Float_t         JetJEC[25];   //[NJets]
   Int_t           JetPartonID[25];   //[NJets]
   Float_t         JetGenPt[25];   //[NJets]
   Float_t         JetGenEta[25];   //[NJets]
   Float_t         JetGenPhi[25];   //[NJets]

   // List of branches
   TBranch        *b_Run;   //!
   TBranch        *b_Event;   //!
   TBranch        *b_LumiSec;   //!
   TBranch        *b_m0;   //!
   TBranch        *b_m12;   //!
   TBranch        *b_process;   //!
   TBranch        *b_mGlu;   //!
   TBranch        *b_mLSP;   //!
   TBranch        *b_isTChiSlepSnu;   //!
   TBranch        *b_isRightHanded;   //!
   TBranch        *b_HLT_MU8;   //!
   TBranch        *b_HLT_MU8_PS;   //!
   TBranch        *b_HLT_MU17;   //!
   TBranch        *b_HLT_MU17_PS;   //!
   TBranch        *b_HLT_ELE17_TIGHT;   //!
   TBranch        *b_HLT_ELE17_TIGHT_PS;   //!
   TBranch        *b_HLT_ELE17_JET30_TIGHT;   //!
   TBranch        *b_HLT_ELE17_JET30_TIGHT_PS;   //!
   TBranch        *b_HLT_ELE8_TIGHT;   //!
   TBranch        *b_HLT_ELE8_TIGHT_PS;   //!
   TBranch        *b_HLT_ELE8_JET30_TIGHT;   //!
   TBranch        *b_HLT_ELE8_JET30_TIGHT_PS;   //!
   TBranch        *b_HLT_MU17_MU8;   //!
   TBranch        *b_HLT_MU17_MU8_PS;   //!
   TBranch        *b_HLT_MU17_TKMU8;   //!
   TBranch        *b_HLT_MU17_TKMU8_PS;   //!
   TBranch        *b_HLT_ELE17_ELE8_TIGHT;   //!
   TBranch        *b_HLT_ELE17_ELE8_TIGHT_PS;   //!
   TBranch        *b_HLT_MU8_ELE17_TIGHT;   //!
   TBranch        *b_HLT_MU8_ELE17_TIGHT_PS;   //!
   TBranch        *b_HLT_MU17_ELE8_TIGHT;   //!
   TBranch        *b_HLT_MU17_ELE8_TIGHT_PS;   //!
   TBranch        *b_Rho;   //!
   TBranch        *b_NVrtx;   //!
   TBranch        *b_PUWeight;   //!
   TBranch        *b_NMus;   //!
   TBranch        *b_IsSignalMuon;   //!
   TBranch        *b_MuPt;   //!
   TBranch        *b_MuEta;   //!
   TBranch        *b_MuPhi;   //!
   TBranch        *b_MuCharge;   //!
   TBranch        *b_MuPFIso;   //!
   TBranch        *b_MuPFChIso;   //!
   TBranch        *b_MuPFNeIso;   //!
   TBranch        *b_MuRadIso;   //!
   TBranch        *b_MuD0;   //!
   TBranch        *b_MuDz;   //!
   TBranch        *b_MuEMVetoEt;   //!
   TBranch        *b_MuHadVetoEt;   //!
   TBranch        *b_MuPassesTightID;   //!
   TBranch        *b_MuPtE;   //!
   TBranch        *b_MuGenID;   //!
   TBranch        *b_MuGenMID;   //!
   TBranch        *b_MuGenGMID;   //!
   TBranch        *b_MuGenType;   //!
   TBranch        *b_MuGenMType;   //!
   TBranch        *b_MuGenGMType;   //!
   TBranch        *b_MuMT;   //!
   TBranch        *b_NEls;   //!
   TBranch        *b_IsSignalElectron;   //!
   TBranch        *b_ElCharge;   //!
   TBranch        *b_ElChIsCons;   //!
   TBranch        *b_ElPt;   //!
   TBranch        *b_ElEta;   //!
   TBranch        *b_ElPhi;   //!
   TBranch        *b_ElD0;   //!
   TBranch        *b_ElD0Err;   //!
   TBranch        *b_ElDz;   //!
   TBranch        *b_ElDzErr;   //!
   TBranch        *b_ElPFIso;   //!
   TBranch        *b_ElPFChIso;   //!
   TBranch        *b_ElPFNeIso;   //!
   TBranch        *b_ElRadIso;   //!
   TBranch        *b_ElEcalRecHitSumEt;   //!
   TBranch        *b_ElHcalTowerSumEt;   //!
   TBranch        *b_ElTkSumPt;   //!
   TBranch        *b_ElDPhi;   //!
   TBranch        *b_ElDEta;   //!
   TBranch        *b_ElSigmaIetaIeta;   //!
   TBranch        *b_ElHoverE;   //!
   TBranch        *b_ElEPthing;   //!
   TBranch        *b_ElIsGoodElId_LooseWP;   //!
   TBranch        *b_ElIsGoodElId_MediumWP;   //!
   TBranch        *b_ElIsGoodTriggerEl;   //!
   TBranch        *b_ElGenID;   //!
   TBranch        *b_ElGenMID;   //!
   TBranch        *b_ElGenGMID;   //!
   TBranch        *b_ElGenType;   //!
   TBranch        *b_ElGenMType;   //!
   TBranch        *b_ElGenGMType;   //!
   TBranch        *b_ElMT;   //!
   TBranch        *b_NTaus;   //!
   TBranch        *b_TauCharge;   //!
   TBranch        *b_TauPt;   //!
   TBranch        *b_TauEta;   //!
   TBranch        *b_TauPhi;   //!
   TBranch        *b_TauMVAElRej;   //!
   TBranch        *b_TauTightMuRej;   //!
   TBranch        *b_TauLCombIsoDB;   //!
   TBranch        *b_pfMET;   //!
   TBranch        *b_pfMETPhi;   //!
   TBranch        *b_pfMETType1;   //!
   TBranch        *b_pfMETType1Phi;   //!
   TBranch        *b_NJets;   //!
   TBranch        *b_JetPt;   //!
   TBranch        *b_JetEta;   //!
   TBranch        *b_JetPhi;   //!
   TBranch        *b_JetEnergy;   //!
   TBranch        *b_JetCSVBTag;   //!
   TBranch        *b_JetProbBTag;   //!
   TBranch        *b_JetArea;   //!
   TBranch        *b_JetJEC;   //!
   TBranch        *b_JetPartonID;   //!
   TBranch        *b_JetGenPt;   //!
   TBranch        *b_JetGenEta;   //!
   TBranch        *b_JetGenPhi;   //!

   MetaTreeClassBase(TTree *tree=0);
   virtual ~MetaTreeClassBase();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef MetaTreeClassBase_cxx
MetaTreeClassBase::MetaTreeClassBase(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("tautest.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("tautest.root");
      }
      f->GetObject("Analysis",tree);

   }
   Init(tree);
}

MetaTreeClassBase::~MetaTreeClassBase()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t MetaTreeClassBase::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t MetaTreeClassBase::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void MetaTreeClassBase::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("Run", &Run, &b_Run);
   fChain->SetBranchAddress("Event", &Event, &b_Event);
   fChain->SetBranchAddress("LumiSec", &LumiSec, &b_LumiSec);
   fChain->SetBranchAddress("m0", &m0, &b_m0);
   fChain->SetBranchAddress("m12", &m12, &b_m12);
   fChain->SetBranchAddress("process", &process, &b_process);
   fChain->SetBranchAddress("mGlu", &mGlu, &b_mGlu);
   fChain->SetBranchAddress("mLSP", &mLSP, &b_mLSP);
   fChain->SetBranchAddress("isTChiSlepSnu", &isTChiSlepSnu, &b_isTChiSlepSnu);
   fChain->SetBranchAddress("isRightHanded", &isRightHanded, &b_isRightHanded);
   fChain->SetBranchAddress("HLT_MU8", &HLT_MU8, &b_HLT_MU8);
   fChain->SetBranchAddress("HLT_MU8_PS", &HLT_MU8_PS, &b_HLT_MU8_PS);
   fChain->SetBranchAddress("HLT_MU17", &HLT_MU17, &b_HLT_MU17);
   fChain->SetBranchAddress("HLT_MU17_PS", &HLT_MU17_PS, &b_HLT_MU17_PS);
   fChain->SetBranchAddress("HLT_ELE17_TIGHT", &HLT_ELE17_TIGHT, &b_HLT_ELE17_TIGHT);
   fChain->SetBranchAddress("HLT_ELE17_TIGHT_PS", &HLT_ELE17_TIGHT_PS, &b_HLT_ELE17_TIGHT_PS);
   fChain->SetBranchAddress("HLT_ELE17_JET30_TIGHT", &HLT_ELE17_JET30_TIGHT, &b_HLT_ELE17_JET30_TIGHT);
   fChain->SetBranchAddress("HLT_ELE17_JET30_TIGHT_PS", &HLT_ELE17_JET30_TIGHT_PS, &b_HLT_ELE17_JET30_TIGHT_PS);
   fChain->SetBranchAddress("HLT_ELE8_TIGHT", &HLT_ELE8_TIGHT, &b_HLT_ELE8_TIGHT);
   fChain->SetBranchAddress("HLT_ELE8_TIGHT_PS", &HLT_ELE8_TIGHT_PS, &b_HLT_ELE8_TIGHT_PS);
   fChain->SetBranchAddress("HLT_ELE8_JET30_TIGHT", &HLT_ELE8_JET30_TIGHT, &b_HLT_ELE8_JET30_TIGHT);
   fChain->SetBranchAddress("HLT_ELE8_JET30_TIGHT_PS", &HLT_ELE8_JET30_TIGHT_PS, &b_HLT_ELE8_JET30_TIGHT_PS);
   fChain->SetBranchAddress("HLT_MU17_MU8", &HLT_MU17_MU8, &b_HLT_MU17_MU8);
   fChain->SetBranchAddress("HLT_MU17_MU8_PS", &HLT_MU17_MU8_PS, &b_HLT_MU17_MU8_PS);
   fChain->SetBranchAddress("HLT_MU17_TKMU8", &HLT_MU17_TKMU8, &b_HLT_MU17_TKMU8);
   fChain->SetBranchAddress("HLT_MU17_TKMU8_PS", &HLT_MU17_TKMU8_PS, &b_HLT_MU17_TKMU8_PS);
   fChain->SetBranchAddress("HLT_ELE17_ELE8_TIGHT", &HLT_ELE17_ELE8_TIGHT, &b_HLT_ELE17_ELE8_TIGHT);
   fChain->SetBranchAddress("HLT_ELE17_ELE8_TIGHT_PS", &HLT_ELE17_ELE8_TIGHT_PS, &b_HLT_ELE17_ELE8_TIGHT_PS);
   fChain->SetBranchAddress("HLT_MU8_ELE17_TIGHT", &HLT_MU8_ELE17_TIGHT, &b_HLT_MU8_ELE17_TIGHT);
   fChain->SetBranchAddress("HLT_MU8_ELE17_TIGHT_PS", &HLT_MU8_ELE17_TIGHT_PS, &b_HLT_MU8_ELE17_TIGHT_PS);
   fChain->SetBranchAddress("HLT_MU17_ELE8_TIGHT", &HLT_MU17_ELE8_TIGHT, &b_HLT_MU17_ELE8_TIGHT);
   fChain->SetBranchAddress("HLT_MU17_ELE8_TIGHT_PS", &HLT_MU17_ELE8_TIGHT_PS, &b_HLT_MU17_ELE8_TIGHT_PS);
   fChain->SetBranchAddress("Rho", &Rho, &b_Rho);
   fChain->SetBranchAddress("NVrtx", &NVrtx, &b_NVrtx);
   fChain->SetBranchAddress("PUWeight", &PUWeight, &b_PUWeight);
   fChain->SetBranchAddress("NMus", &NMus, &b_NMus);
   fChain->SetBranchAddress("IsSignalMuon", IsSignalMuon, &b_IsSignalMuon);
   fChain->SetBranchAddress("MuPt", MuPt, &b_MuPt);
   fChain->SetBranchAddress("MuEta", MuEta, &b_MuEta);
   fChain->SetBranchAddress("MuPhi", MuPhi, &b_MuPhi);
   fChain->SetBranchAddress("MuCharge", MuCharge, &b_MuCharge);
   fChain->SetBranchAddress("MuPFIso", MuPFIso, &b_MuPFIso);
   fChain->SetBranchAddress("MuPFChIso", MuPFChIso, &b_MuPFChIso);
   fChain->SetBranchAddress("MuPFNeIso", MuPFNeIso, &b_MuPFNeIso);
   fChain->SetBranchAddress("MuRadIso", MuRadIso, &b_MuRadIso);
   fChain->SetBranchAddress("MuD0", MuD0, &b_MuD0);
   fChain->SetBranchAddress("MuDz", MuDz, &b_MuDz);
   fChain->SetBranchAddress("MuEMVetoEt", MuEMVetoEt, &b_MuEMVetoEt);
   fChain->SetBranchAddress("MuHadVetoEt", MuHadVetoEt, &b_MuHadVetoEt);
   fChain->SetBranchAddress("MuPassesTightID", MuPassesTightID, &b_MuPassesTightID);
   fChain->SetBranchAddress("MuPtE", MuPtE, &b_MuPtE);
   fChain->SetBranchAddress("MuGenID", MuGenID, &b_MuGenID);
   fChain->SetBranchAddress("MuGenMID", MuGenMID, &b_MuGenMID);
   fChain->SetBranchAddress("MuGenGMID", MuGenGMID, &b_MuGenGMID);
   fChain->SetBranchAddress("MuGenType", MuGenType, &b_MuGenType);
   fChain->SetBranchAddress("MuGenMType", MuGenMType, &b_MuGenMType);
   fChain->SetBranchAddress("MuGenGMType", MuGenGMType, &b_MuGenGMType);
   fChain->SetBranchAddress("MuMT", MuMT, &b_MuMT);
   fChain->SetBranchAddress("NEls", &NEls, &b_NEls);
   fChain->SetBranchAddress("IsSignalElectron", IsSignalElectron, &b_IsSignalElectron);
   fChain->SetBranchAddress("ElCharge", ElCharge, &b_ElCharge);
   fChain->SetBranchAddress("ElChIsCons", ElChIsCons, &b_ElChIsCons);
   fChain->SetBranchAddress("ElPt", ElPt, &b_ElPt);
   fChain->SetBranchAddress("ElEta", ElEta, &b_ElEta);
   fChain->SetBranchAddress("ElPhi", ElPhi, &b_ElPhi);
   fChain->SetBranchAddress("ElD0", ElD0, &b_ElD0);
   fChain->SetBranchAddress("ElD0Err", ElD0Err, &b_ElD0Err);
   fChain->SetBranchAddress("ElDz", ElDz, &b_ElDz);
   fChain->SetBranchAddress("ElDzErr", ElDzErr, &b_ElDzErr);
   fChain->SetBranchAddress("ElPFIso", ElPFIso, &b_ElPFIso);
   fChain->SetBranchAddress("ElPFChIso", ElPFChIso, &b_ElPFChIso);
   fChain->SetBranchAddress("ElPFNeIso", ElPFNeIso, &b_ElPFNeIso);
   fChain->SetBranchAddress("ElRadIso", ElRadIso, &b_ElRadIso);
   fChain->SetBranchAddress("ElEcalRecHitSumEt", ElEcalRecHitSumEt, &b_ElEcalRecHitSumEt);
   fChain->SetBranchAddress("ElHcalTowerSumEt", ElHcalTowerSumEt, &b_ElHcalTowerSumEt);
   fChain->SetBranchAddress("ElTkSumPt", ElTkSumPt, &b_ElTkSumPt);
   fChain->SetBranchAddress("ElDPhi", ElDPhi, &b_ElDPhi);
   fChain->SetBranchAddress("ElDEta", ElDEta, &b_ElDEta);
   fChain->SetBranchAddress("ElSigmaIetaIeta", ElSigmaIetaIeta, &b_ElSigmaIetaIeta);
   fChain->SetBranchAddress("ElHoverE", ElHoverE, &b_ElHoverE);
   fChain->SetBranchAddress("ElEPthing", ElEPthing, &b_ElEPthing);
   fChain->SetBranchAddress("ElIsGoodElId_LooseWP", ElIsGoodElId_LooseWP, &b_ElIsGoodElId_LooseWP);
   fChain->SetBranchAddress("ElIsGoodElId_MediumWP", ElIsGoodElId_MediumWP, &b_ElIsGoodElId_MediumWP);
   fChain->SetBranchAddress("ElIsGoodTriggerEl", ElIsGoodTriggerEl, &b_ElIsGoodTriggerEl);
   fChain->SetBranchAddress("ElGenID", ElGenID, &b_ElGenID);
   fChain->SetBranchAddress("ElGenMID", ElGenMID, &b_ElGenMID);
   fChain->SetBranchAddress("ElGenGMID", ElGenGMID, &b_ElGenGMID);
   fChain->SetBranchAddress("ElGenType", ElGenType, &b_ElGenType);
   fChain->SetBranchAddress("ElGenMType", ElGenMType, &b_ElGenMType);
   fChain->SetBranchAddress("ElGenGMType", ElGenGMType, &b_ElGenGMType);
   fChain->SetBranchAddress("ElMT", ElMT, &b_ElMT);
   fChain->SetBranchAddress("NTaus", &NTaus, &b_NTaus);
   fChain->SetBranchAddress("TauCharge", TauCharge, &b_TauCharge);
   fChain->SetBranchAddress("TauPt", TauPt, &b_TauPt);
   fChain->SetBranchAddress("TauEta", TauEta, &b_TauEta);
   fChain->SetBranchAddress("TauPhi", TauPhi, &b_TauPhi);
   fChain->SetBranchAddress("TauMVAElRej", TauMVAElRej, &b_TauMVAElRej);
   fChain->SetBranchAddress("TauTightMuRej", TauTightMuRej, &b_TauTightMuRej);
   fChain->SetBranchAddress("TauLCombIsoDB", TauLCombIsoDB, &b_TauLCombIsoDB);
   fChain->SetBranchAddress("pfMET", &pfMET, &b_pfMET);
   fChain->SetBranchAddress("pfMETPhi", &pfMETPhi, &b_pfMETPhi);
   fChain->SetBranchAddress("pfMETType1", &pfMETType1, &b_pfMETType1);
   fChain->SetBranchAddress("pfMETType1Phi", &pfMETType1Phi, &b_pfMETType1Phi);
   fChain->SetBranchAddress("NJets", &NJets, &b_NJets);
   fChain->SetBranchAddress("JetPt", JetPt, &b_JetPt);
   fChain->SetBranchAddress("JetEta", JetEta, &b_JetEta);
   fChain->SetBranchAddress("JetPhi", JetPhi, &b_JetPhi);
   fChain->SetBranchAddress("JetEnergy", JetEnergy, &b_JetEnergy);
   fChain->SetBranchAddress("JetCSVBTag", JetCSVBTag, &b_JetCSVBTag);
   fChain->SetBranchAddress("JetProbBTag", JetProbBTag, &b_JetProbBTag);
   fChain->SetBranchAddress("JetArea", JetArea, &b_JetArea);
   fChain->SetBranchAddress("JetJEC", JetJEC, &b_JetJEC);
   fChain->SetBranchAddress("JetPartonID", JetPartonID, &b_JetPartonID);
   fChain->SetBranchAddress("JetGenPt", JetGenPt, &b_JetGenPt);
   fChain->SetBranchAddress("JetGenEta", JetGenEta, &b_JetGenEta);
   fChain->SetBranchAddress("JetGenPhi", JetGenPhi, &b_JetGenPhi);
   Notify();
}

Bool_t MetaTreeClassBase::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void MetaTreeClassBase::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t MetaTreeClassBase::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef MetaTreeClassBase_cxx
