//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Nov 15 17:07:04 2011 by ROOT version 5.27/06b
// from TTree Analysis/AnalysisTree
// found on file: /scratch/stiegerb/SSDLTrees/2011B/Nov15/MC/LM5_SUSY_sftsht_7TeV-pythia6.root
//////////////////////////////////////////////////////////

#ifndef MetaTreeClassBase_h
#define MetaTreeClassBase_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

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
   Int_t           HLT_MU8_JET40;
   Int_t           HLT_MU8_JET40_PS;
   Int_t           HLT_ELE8_JET40;
   Int_t           HLT_ELE8_JET40_PS;
   Int_t           HLT_DOUBLEMU7;
   Int_t           HLT_DOUBLEMU7_PS;
   Int_t           HLT_MU13_MU8;
   Int_t           HLT_MU13_MU8_PS;
   Int_t           HLT_MU17_MU8;
   Int_t           HLT_MU17_MU8_PS;
   Int_t           HLT_ELE17_ELE8;
   Int_t           HLT_ELE17_ELE8_PS;
   Int_t           HLT_ELE17_ELE8_TIGHT;
   Int_t           HLT_ELE17_ELE8_TIGHT_PS;
   Int_t           HLT_MU17_ELE8;
   Int_t           HLT_MU17_ELE8_PS;
   Int_t           HLT_MU8_ELE17;
   Int_t           HLT_MU8_ELE17_PS;
   Int_t           HLT_MU8_ELE17_TIGHT;
   Int_t           HLT_MU8_ELE17_TIGHT_PS;
   Int_t           HLT_MU17_ELE8_TIGHT;
   Int_t           HLT_MU17_ELE8_TIGHT_PS;
   Int_t           HLT_DOUBLEELE8_HT160;
   Int_t           HLT_DOUBLEELE8_HT160_PS;
   Int_t           HLT_DOUBLEELE8_HT160_TIGHT;
   Int_t           HLT_DOUBLEELE8_HT160_TIGHT_PS;
   Int_t           HLT_DOUBLEMU3_HT160;
   Int_t           HLT_DOUBLEMU3_HT160_PS;
   Int_t           HLT_MU3_ELE8_HT160;
   Int_t           HLT_MU3_ELE8_HT160_PS;
   Int_t           HLT_MU3_ELE8_HT160_TIGHT;
   Int_t           HLT_MU3_ELE8_HT160_TIGHT_PS;
   Int_t           HLT_DOUBLEMU3_MASS4_HT150;
   Int_t           HLT_DOUBLEMU3_MASS4_HT150_PS;
   Int_t           HLT_DOUBLEELE8_CALOIDT_TRKIDVL_MASS4_HT150;
   Int_t           HLT_DOUBLEELE8_CALOIDT_TRKIDVL_MASS4_HT150_PS;
   Int_t           HLT_MU5_ELE8_CALOIDT_TRKIDVL_MASS4_HT150;
   Int_t           HLT_MU5_ELE8_CALOIDT_TRKIDVL_MASS4_HT150_PS;
   Float_t         Rho;
   Int_t           NVrtx;
   Float_t         PUWeight;
   Int_t           NMus;
   Int_t           IsSignalMuon[5];   //[NMus]
   Float_t         MuPt[5];   //[NMus]
   Float_t         MuEta[5];   //[NMus]
   Float_t         MuPhi[5];   //[NMus]
   Int_t           MuCharge[5];   //[NMus]
   Float_t         MuIso[5];   //[NMus]
   Float_t         MuD0[5];   //[NMus]
   Float_t         MuDz[5];   //[NMus]
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
   Float_t         ElRelIso[5];   //[NEls]
   Float_t         ElEcalRecHitSumEt[5];   //[NEls]
   Float_t         ElHcalTowerSumEt[5];   //[NEls]
   Float_t         ElTkSumPt[5];   //[NEls]
   Float_t         ElDPhi[5];   //[NEls]
   Float_t         ElDEta[5];   //[NEls]
   Float_t         ElSigmaIetaIeta[5];   //[NEls]
   Float_t         ElHoverE[5];   //[NEls]
   Int_t           ElIsGoodElId_WP80[5];   //[NEls]
   Int_t           ElIsGoodElId_WP90[5];   //[NEls]
   Int_t           ElGenID[5];   //[NEls]
   Int_t           ElGenMID[5];   //[NEls]
   Int_t           ElGenGMID[5];   //[NEls]
   Int_t           ElGenType[5];   //[NEls]
   Int_t           ElGenMType[5];   //[NEls]
   Int_t           ElGenGMType[5];   //[NEls]
   Float_t         ElMT[5];   //[NEls]
   Float_t         tcMET;
   Float_t         tcMETPhi;
   Float_t         pfMET;
   Float_t         pfMETPhi;
   Int_t           NJets;
   Float_t         JetPt[20];   //[NJets]
   Float_t         JetEta[20];   //[NJets]
   Float_t         JetPhi[20];   //[NJets]
   Float_t         JetSSVHPBTag[20];   //[NJets]
   Float_t         JetArea[20];   //[NJets]

   // List of branches
   TBranch        *b_Run;   //!
   TBranch        *b_Event;   //!
   TBranch        *b_LumiSec;   //!
   TBranch        *b_m0;   //!
   TBranch        *b_m12;   //!
   TBranch        *b_process;   //!
   TBranch        *b_HLT_MU8_JET40;   //!
   TBranch        *b_HLT_MU8_JET40_PS;   //!
   TBranch        *b_HLT_ELE8_JET40;   //!
   TBranch        *b_HLT_ELE8_JET40_PS;   //!
   TBranch        *b_HLT_DOUBLEMU7;   //!
   TBranch        *b_HLT_DOUBLEMU7_PS;   //!
   TBranch        *b_HLT_MU13_MU8;   //!
   TBranch        *b_HLT_MU13_MU8_PS;   //!
   TBranch        *b_HLT_MU17_MU8;   //!
   TBranch        *b_HLT_MU17_MU8_PS;   //!
   TBranch        *b_HLT_ELE17_ELE8;   //!
   TBranch        *b_HLT_ELE17_ELE8_PS;   //!
   TBranch        *b_HLT_ELE17_ELE8_TIGHT;   //!
   TBranch        *b_HLT_ELE17_ELE8_TIGHT_PS;   //!
   TBranch        *b_HLT_MU17_ELE8;   //!
   TBranch        *b_HLT_MU17_ELE8_PS;   //!
   TBranch        *b_HLT_MU8_ELE17;   //!
   TBranch        *b_HLT_MU8_ELE17_PS;   //!
   TBranch        *b_HLT_MU8_ELE17_TIGHT;   //!
   TBranch        *b_HLT_MU8_ELE17_TIGHT_PS;   //!
   TBranch        *b_HLT_MU17_ELE8_TIGHT;   //!
   TBranch        *b_HLT_MU17_ELE8_TIGHT_PS;   //!
   TBranch        *b_HLT_DOUBLEELE8_HT160;   //!
   TBranch        *b_HLT_DOUBLEELE8_HT160_PS;   //!
   TBranch        *b_HLT_DOUBLEELE8_HT160_TIGHT;   //!
   TBranch        *b_HLT_DOUBLEELE8_HT160_TIGHT_PS;   //!
   TBranch        *b_HLT_DOUBLEMU3_HT160;   //!
   TBranch        *b_HLT_DOUBLEMU3_HT160_PS;   //!
   TBranch        *b_HLT_MU3_ELE8_HT160;   //!
   TBranch        *b_HLT_MU3_ELE8_HT160_PS;   //!
   TBranch        *b_HLT_MU3_ELE8_HT160_TIGHT;   //!
   TBranch        *b_HLT_MU3_ELE8_HT160_TIGHT_PS;   //!
   TBranch        *b_HLT_DOUBLEMU3_MASS4_HT150;   //!
   TBranch        *b_HLT_DOUBLEMU3_MASS4_HT150_PS;   //!
   TBranch        *b_HLT_DOUBLEELE8_CALOIDT_TRKIDVL_MASS4_HT150;   //!
   TBranch        *b_HLT_DOUBLEELE8_CALOIDT_TRKIDVL_MASS4_HT150_PS;   //!
   TBranch        *b_HLT_MU5_ELE8_CALOIDT_TRKIDVL_MASS4_HT150;   //!
   TBranch        *b_HLT_MU5_ELE8_CALOIDT_TRKIDVL_MASS4_HT150_PS;   //!
   TBranch        *b_Rho;   //!
   TBranch        *b_NVrtx;   //!
   TBranch        *b_PUWeight;   //!
   TBranch        *b_NMus;   //!
   TBranch        *b_IsSignalMuon;   //!
   TBranch        *b_MuPt;   //!
   TBranch        *b_MuEta;   //!
   TBranch        *b_MuPhi;   //!
   TBranch        *b_MuCharge;   //!
   TBranch        *b_MuIso;   //!
   TBranch        *b_MuD0;   //!
   TBranch        *b_MuDz;   //!
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
   TBranch        *b_ElRelIso;   //!
   TBranch        *b_ElEcalRecHitSumEt;   //!
   TBranch        *b_ElHcalTowerSumEt;   //!
   TBranch        *b_ElTkSumPt;   //!
   TBranch        *b_ElDPhi;   //!
   TBranch        *b_ElDEta;   //!
   TBranch        *b_ElSigmaIetaIeta;   //!
   TBranch        *b_ElHoverE;   //!
   TBranch        *b_ElIsGoodElId_WP80;   //!
   TBranch        *b_ElIsGoodElId_WP90;   //!
   TBranch        *b_ElGenID;   //!
   TBranch        *b_ElGenMID;   //!
   TBranch        *b_ElGenGMID;   //!
   TBranch        *b_ElGenType;   //!
   TBranch        *b_ElGenMType;   //!
   TBranch        *b_ElGenGMType;   //!
   TBranch        *b_ElMT;   //!
   TBranch        *b_tcMET;   //!
   TBranch        *b_tcMETPhi;   //!
   TBranch        *b_pfMET;   //!
   TBranch        *b_pfMETPhi;   //!
   TBranch        *b_NJets;   //!
   TBranch        *b_JetPt;   //!
   TBranch        *b_JetEta;   //!
   TBranch        *b_JetPhi;   //!
   TBranch        *b_JetSSVHPBTag;   //!
   TBranch        *b_JetArea;   //!

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
MetaTreeClassBase::MetaTreeClassBase(TTree *tree)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/scratch/stiegerb/SSDLTrees/2011B/Nov15/MC/LM5_SUSY_sftsht_7TeV-pythia6.root");
      if (!f) {
         f = new TFile("/scratch/stiegerb/SSDLTrees/2011B/Nov15/MC/LM5_SUSY_sftsht_7TeV-pythia6.root");
      }
      tree = (TTree*)gDirectory->Get("Analysis");

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
   if (!fChain->InheritsFrom(TChain::Class()))  return centry;
   TChain *chain = (TChain*)fChain;
   if (chain->GetTreeNumber() != fCurrent) {
      fCurrent = chain->GetTreeNumber();
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
   fChain->SetBranchAddress("HLT_MU8_JET40", &HLT_MU8_JET40, &b_HLT_MU8_JET40);
   fChain->SetBranchAddress("HLT_MU8_JET40_PS", &HLT_MU8_JET40_PS, &b_HLT_MU8_JET40_PS);
   fChain->SetBranchAddress("HLT_ELE8_JET40", &HLT_ELE8_JET40, &b_HLT_ELE8_JET40);
   fChain->SetBranchAddress("HLT_ELE8_JET40_PS", &HLT_ELE8_JET40_PS, &b_HLT_ELE8_JET40_PS);
   fChain->SetBranchAddress("HLT_DOUBLEMU7", &HLT_DOUBLEMU7, &b_HLT_DOUBLEMU7);
   fChain->SetBranchAddress("HLT_DOUBLEMU7_PS", &HLT_DOUBLEMU7_PS, &b_HLT_DOUBLEMU7_PS);
   fChain->SetBranchAddress("HLT_MU13_MU8", &HLT_MU13_MU8, &b_HLT_MU13_MU8);
   fChain->SetBranchAddress("HLT_MU13_MU8_PS", &HLT_MU13_MU8_PS, &b_HLT_MU13_MU8_PS);
   fChain->SetBranchAddress("HLT_MU17_MU8", &HLT_MU17_MU8, &b_HLT_MU17_MU8);
   fChain->SetBranchAddress("HLT_MU17_MU8_PS", &HLT_MU17_MU8_PS, &b_HLT_MU17_MU8_PS);
   fChain->SetBranchAddress("HLT_ELE17_ELE8", &HLT_ELE17_ELE8, &b_HLT_ELE17_ELE8);
   fChain->SetBranchAddress("HLT_ELE17_ELE8_PS", &HLT_ELE17_ELE8_PS, &b_HLT_ELE17_ELE8_PS);
   fChain->SetBranchAddress("HLT_ELE17_ELE8_TIGHT", &HLT_ELE17_ELE8_TIGHT, &b_HLT_ELE17_ELE8_TIGHT);
   fChain->SetBranchAddress("HLT_ELE17_ELE8_TIGHT_PS", &HLT_ELE17_ELE8_TIGHT_PS, &b_HLT_ELE17_ELE8_TIGHT_PS);
   fChain->SetBranchAddress("HLT_MU17_ELE8", &HLT_MU17_ELE8, &b_HLT_MU17_ELE8);
   fChain->SetBranchAddress("HLT_MU17_ELE8_PS", &HLT_MU17_ELE8_PS, &b_HLT_MU17_ELE8_PS);
   fChain->SetBranchAddress("HLT_MU8_ELE17", &HLT_MU8_ELE17, &b_HLT_MU8_ELE17);
   fChain->SetBranchAddress("HLT_MU8_ELE17_PS", &HLT_MU8_ELE17_PS, &b_HLT_MU8_ELE17_PS);
   fChain->SetBranchAddress("HLT_MU8_ELE17_TIGHT", &HLT_MU8_ELE17_TIGHT, &b_HLT_MU8_ELE17_TIGHT);
   fChain->SetBranchAddress("HLT_MU8_ELE17_TIGHT_PS", &HLT_MU8_ELE17_TIGHT_PS, &b_HLT_MU8_ELE17_TIGHT_PS);
   fChain->SetBranchAddress("HLT_MU17_ELE8_TIGHT", &HLT_MU17_ELE8_TIGHT, &b_HLT_MU17_ELE8_TIGHT);
   fChain->SetBranchAddress("HLT_MU17_ELE8_TIGHT_PS", &HLT_MU17_ELE8_TIGHT_PS, &b_HLT_MU17_ELE8_TIGHT_PS);
   fChain->SetBranchAddress("HLT_DOUBLEELE8_HT160", &HLT_DOUBLEELE8_HT160, &b_HLT_DOUBLEELE8_HT160);
   fChain->SetBranchAddress("HLT_DOUBLEELE8_HT160_PS", &HLT_DOUBLEELE8_HT160_PS, &b_HLT_DOUBLEELE8_HT160_PS);
   fChain->SetBranchAddress("HLT_DOUBLEELE8_HT160_TIGHT", &HLT_DOUBLEELE8_HT160_TIGHT, &b_HLT_DOUBLEELE8_HT160_TIGHT);
   fChain->SetBranchAddress("HLT_DOUBLEELE8_HT160_TIGHT_PS", &HLT_DOUBLEELE8_HT160_TIGHT_PS, &b_HLT_DOUBLEELE8_HT160_TIGHT_PS);
   fChain->SetBranchAddress("HLT_DOUBLEMU3_HT160", &HLT_DOUBLEMU3_HT160, &b_HLT_DOUBLEMU3_HT160);
   fChain->SetBranchAddress("HLT_DOUBLEMU3_HT160_PS", &HLT_DOUBLEMU3_HT160_PS, &b_HLT_DOUBLEMU3_HT160_PS);
   fChain->SetBranchAddress("HLT_MU3_ELE8_HT160", &HLT_MU3_ELE8_HT160, &b_HLT_MU3_ELE8_HT160);
   fChain->SetBranchAddress("HLT_MU3_ELE8_HT160_PS", &HLT_MU3_ELE8_HT160_PS, &b_HLT_MU3_ELE8_HT160_PS);
   fChain->SetBranchAddress("HLT_MU3_ELE8_HT160_TIGHT", &HLT_MU3_ELE8_HT160_TIGHT, &b_HLT_MU3_ELE8_HT160_TIGHT);
   fChain->SetBranchAddress("HLT_MU3_ELE8_HT160_TIGHT_PS", &HLT_MU3_ELE8_HT160_TIGHT_PS, &b_HLT_MU3_ELE8_HT160_TIGHT_PS);
   fChain->SetBranchAddress("HLT_DOUBLEMU3_MASS4_HT150", &HLT_DOUBLEMU3_MASS4_HT150, &b_HLT_DOUBLEMU3_MASS4_HT150);
   fChain->SetBranchAddress("HLT_DOUBLEMU3_MASS4_HT150_PS", &HLT_DOUBLEMU3_MASS4_HT150_PS, &b_HLT_DOUBLEMU3_MASS4_HT150_PS);
   fChain->SetBranchAddress("HLT_DOUBLEELE8_CALOIDT_TRKIDVL_MASS4_HT150", &HLT_DOUBLEELE8_CALOIDT_TRKIDVL_MASS4_HT150, &b_HLT_DOUBLEELE8_CALOIDT_TRKIDVL_MASS4_HT150);
   fChain->SetBranchAddress("HLT_DOUBLEELE8_CALOIDT_TRKIDVL_MASS4_HT150_PS", &HLT_DOUBLEELE8_CALOIDT_TRKIDVL_MASS4_HT150_PS, &b_HLT_DOUBLEELE8_CALOIDT_TRKIDVL_MASS4_HT150_PS);
   fChain->SetBranchAddress("HLT_MU5_ELE8_CALOIDT_TRKIDVL_MASS4_HT150", &HLT_MU5_ELE8_CALOIDT_TRKIDVL_MASS4_HT150, &b_HLT_MU5_ELE8_CALOIDT_TRKIDVL_MASS4_HT150);
   fChain->SetBranchAddress("HLT_MU5_ELE8_CALOIDT_TRKIDVL_MASS4_HT150_PS", &HLT_MU5_ELE8_CALOIDT_TRKIDVL_MASS4_HT150_PS, &b_HLT_MU5_ELE8_CALOIDT_TRKIDVL_MASS4_HT150_PS);
   fChain->SetBranchAddress("Rho", &Rho, &b_Rho);
   fChain->SetBranchAddress("NVrtx", &NVrtx, &b_NVrtx);
   fChain->SetBranchAddress("PUWeight", &PUWeight, &b_PUWeight);
   fChain->SetBranchAddress("NMus", &NMus, &b_NMus);
   fChain->SetBranchAddress("MuPt", MuPt, &b_MuPt);
   fChain->SetBranchAddress("MuEta", MuEta, &b_MuEta);
   fChain->SetBranchAddress("MuPhi", MuPhi, &b_MuPhi);
   fChain->SetBranchAddress("MuCharge", MuCharge, &b_MuCharge);
   fChain->SetBranchAddress("MuIso", MuIso, &b_MuIso);
   fChain->SetBranchAddress("MuD0", MuD0, &b_MuD0);
   fChain->SetBranchAddress("MuDz", MuDz, &b_MuDz);
   fChain->SetBranchAddress("MuPtE", MuPtE, &b_MuPtE);
   fChain->SetBranchAddress("MuGenID", MuGenID, &b_MuGenID);
   fChain->SetBranchAddress("MuGenMID", MuGenMID, &b_MuGenMID);
   fChain->SetBranchAddress("MuGenGMID", MuGenGMID, &b_MuGenGMID);
   fChain->SetBranchAddress("MuGenType", MuGenType, &b_MuGenType);
   fChain->SetBranchAddress("MuGenMType", MuGenMType, &b_MuGenMType);
   fChain->SetBranchAddress("MuGenGMType", MuGenGMType, &b_MuGenGMType);
   fChain->SetBranchAddress("MuMT", MuMT, &b_MuMT);
   fChain->SetBranchAddress("NEls", &NEls, &b_NEls);
   fChain->SetBranchAddress("ElCharge", ElCharge, &b_ElCharge);
   fChain->SetBranchAddress("ElChIsCons", ElChIsCons, &b_ElChIsCons);
   fChain->SetBranchAddress("ElPt", ElPt, &b_ElPt);
   fChain->SetBranchAddress("ElEta", ElEta, &b_ElEta);
   fChain->SetBranchAddress("ElPhi", ElPhi, &b_ElPhi);
   fChain->SetBranchAddress("ElD0", ElD0, &b_ElD0);
   fChain->SetBranchAddress("ElD0Err", ElD0Err, &b_ElD0Err);
   fChain->SetBranchAddress("ElDz", ElDz, &b_ElDz);
   fChain->SetBranchAddress("ElDzErr", ElDzErr, &b_ElDzErr);
   fChain->SetBranchAddress("ElRelIso", ElRelIso, &b_ElRelIso);
   fChain->SetBranchAddress("ElEcalRecHitSumEt", ElEcalRecHitSumEt, &b_ElEcalRecHitSumEt);
   fChain->SetBranchAddress("ElHcalTowerSumEt", ElHcalTowerSumEt, &b_ElHcalTowerSumEt);
   fChain->SetBranchAddress("ElTkSumPt", ElTkSumPt, &b_ElTkSumPt);
   fChain->SetBranchAddress("ElDPhi", ElDPhi, &b_ElDPhi);
   fChain->SetBranchAddress("ElDEta", ElDEta, &b_ElDEta);
   fChain->SetBranchAddress("ElSigmaIetaIeta", ElSigmaIetaIeta, &b_ElSigmaIetaIeta);
   fChain->SetBranchAddress("ElHoverE", ElHoverE, &b_ElHoverE);
   fChain->SetBranchAddress("ElIsGoodElId_WP80", ElIsGoodElId_WP80, &b_ElIsGoodElId_WP80);
   fChain->SetBranchAddress("ElIsGoodElId_WP90", ElIsGoodElId_WP90, &b_ElIsGoodElId_WP90);
   fChain->SetBranchAddress("ElGenID", ElGenID, &b_ElGenID);
   fChain->SetBranchAddress("ElGenMID", ElGenMID, &b_ElGenMID);
   fChain->SetBranchAddress("ElGenGMID", ElGenGMID, &b_ElGenGMID);
   fChain->SetBranchAddress("ElGenType", ElGenType, &b_ElGenType);
   fChain->SetBranchAddress("ElGenMType", ElGenMType, &b_ElGenMType);
   fChain->SetBranchAddress("ElGenGMType", ElGenGMType, &b_ElGenGMType);
   fChain->SetBranchAddress("ElMT", ElMT, &b_ElMT);
   fChain->SetBranchAddress("tcMET", &tcMET, &b_tcMET);
   fChain->SetBranchAddress("tcMETPhi", &tcMETPhi, &b_tcMETPhi);
   fChain->SetBranchAddress("pfMET", &pfMET, &b_pfMET);
   fChain->SetBranchAddress("pfMETPhi", &pfMETPhi, &b_pfMETPhi);
   fChain->SetBranchAddress("NJets", &NJets, &b_NJets);
   fChain->SetBranchAddress("JetPt", JetPt, &b_JetPt);
   fChain->SetBranchAddress("JetEta", JetEta, &b_JetEta);
   fChain->SetBranchAddress("JetPhi", JetPhi, &b_JetPhi);
   fChain->SetBranchAddress("JetSSVHPBTag", JetSSVHPBTag, &b_JetSSVHPBTag);
   fChain->SetBranchAddress("JetArea", JetArea, &b_JetArea);
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
