//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Oct 11 16:19:16 2010 by ROOT version 5.25/04
// from TTree MuonAnalysis/MuonFakeAnalysisTree1
// found on file: MetaTrees/V01-04-01/Mu-MTree/SingleMuonSelection.root
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
   Int_t           LumiSection;
   Int_t           HLTMu9;
   Int_t           HLTDoubleMu3;
   Int_t           HLTJet30U;
   Int_t           HLTJet50U;
   Int_t           HLTJet70U;
   Int_t           HLTJet100U;
   Int_t           NJets;
   Double_t        JPt[8];   //[NJets]
   Double_t        JEta[8];   //[NJets]
   Double_t        JPhi[8];   //[NJets]
   Double_t        HT;
   Double_t        MHT;
   Double_t        MET;
   Double_t        MT;
   Double_t        Minv;
   Double_t        SumET;
   Int_t           NMus;
   Double_t        MuPt[2];
   Double_t        MuEta[2];
   Double_t        MuPhi[2];
   Int_t           MuCharge[2];
   Int_t           MuTight[2];
   Double_t        MuIso[2];
   Double_t        MuDRJet[2];
   Double_t        MuDPhiJet[2];
   Double_t        MuDRHardestJet[2];
   Double_t        MuCaloComp[2];
   Double_t        MuSegmComp[2];
   Double_t        MuOuterRad[2];
   Double_t        MuNChi2[2];
   Int_t           MuNTkHits[2];
   Double_t        MuD0[2];
   Double_t        MuDz[2];
   Double_t        MuPtE[2];
   Int_t           MuGenID[2];
   Int_t           MuGenMoID[2];
   Int_t           MuGenGMoID[2];
   Int_t           MuGenType[2];
   Int_t           MuGenMoType[2];
   Int_t           MuGenGMoType[2];

   // List of branches
   TBranch        *b_Run;   //!
   TBranch        *b_Event;   //!
   TBranch        *b_LumiSection;   //!
   TBranch        *b_HLTMu9;   //!
   TBranch        *b_HLTDoubleMu3;   //!
   TBranch        *b_HLTJet30U;   //!
   TBranch        *b_HLTJet50U;   //!
   TBranch        *b_HLTJet70U;   //!
   TBranch        *b_HLTJet100U;   //!
   TBranch        *b_NJets;   //!
   TBranch        *b_JPt;   //!
   TBranch        *b_JEta;   //!
   TBranch        *b_JPhi;   //!
   TBranch        *b_HT;   //!
   TBranch        *b_MHT;   //!
   TBranch        *b_MET;   //!
   TBranch        *b_MT;   //!
   TBranch        *b_Minv;   //!
   TBranch        *b_SumET;   //!
   TBranch        *b_NMus;   //!
   TBranch        *b_MuPt;   //!
   TBranch        *b_MuEta;   //!
   TBranch        *b_MuPhi;   //!
   TBranch        *b_MuCharge;   //!
   TBranch        *b_MuTight;   //!
   TBranch        *b_MuIso;   //!
   TBranch        *b_MuDRJet;   //!
   TBranch        *b_MuDPhiJet;   //!
   TBranch        *b_MuDRHardestJet;   //!
   TBranch        *b_MuCaloComp;   //!
   TBranch        *b_MuSegmComp;   //!
   TBranch        *b_MuOuterRad;   //!
   TBranch        *b_MuNChi2;   //!
   TBranch        *b_MuNTkHits;   //!
   TBranch        *b_MuD0;   //!
   TBranch        *b_MuDz;   //!
   TBranch        *b_MuPtE;   //!
   TBranch        *b_MuGenID;   //!
   TBranch        *b_MuGenMoID;   //!
   TBranch        *b_MuGenGMoID;   //!
   TBranch        *b_MuGenType;   //!
   TBranch        *b_MuGenMoType;   //!
   TBranch        *b_MuGenGMoType;   //!

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
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("MetaTrees/V01-04-01/Mu-MTree/SingleMuonSelection.root");
      if (!f) {
         f = new TFile("MetaTrees/V01-04-01/Mu-MTree/SingleMuonSelection.root");
      }
      tree = (TTree*)gDirectory->Get("MuonAnalysis");

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
   fChain->SetBranchAddress("LumiSection", &LumiSection, &b_LumiSection);
   fChain->SetBranchAddress("HLTMu9", &HLTMu9, &b_HLTMu9);
   fChain->SetBranchAddress("HLTDoubleMu3", &HLTDoubleMu3, &b_HLTDoubleMu3);
   fChain->SetBranchAddress("HLTJet30U", &HLTJet30U, &b_HLTJet30U);
   fChain->SetBranchAddress("HLTJet50U", &HLTJet50U, &b_HLTJet50U);
   fChain->SetBranchAddress("HLTJet70U", &HLTJet70U, &b_HLTJet70U);
   fChain->SetBranchAddress("HLTJet100U", &HLTJet100U, &b_HLTJet100U);
   fChain->SetBranchAddress("NJets", &NJets, &b_NJets);
   fChain->SetBranchAddress("JPt", JPt, &b_JPt);
   fChain->SetBranchAddress("JEta", JEta, &b_JEta);
   fChain->SetBranchAddress("JPhi", JPhi, &b_JPhi);
   fChain->SetBranchAddress("HT", &HT, &b_HT);
   fChain->SetBranchAddress("MHT", &MHT, &b_MHT);
   fChain->SetBranchAddress("MET", &MET, &b_MET);
   fChain->SetBranchAddress("MT", &MT, &b_MT);
   fChain->SetBranchAddress("Minv", &Minv, &b_Minv);
   fChain->SetBranchAddress("SumET", &SumET, &b_SumET);
   fChain->SetBranchAddress("NMus", &NMus, &b_NMus);
   fChain->SetBranchAddress("MuPt", MuPt, &b_MuPt);
   fChain->SetBranchAddress("MuEta", MuEta, &b_MuEta);
   fChain->SetBranchAddress("MuPhi", MuPhi, &b_MuPhi);
   fChain->SetBranchAddress("MuCharge", MuCharge, &b_MuCharge);
   fChain->SetBranchAddress("MuTight", MuTight, &b_MuTight);
   fChain->SetBranchAddress("MuIso", MuIso, &b_MuIso);
   fChain->SetBranchAddress("MuDRJet", MuDRJet, &b_MuDRJet);
   fChain->SetBranchAddress("MuDPhiJet", MuDPhiJet, &b_MuDPhiJet);
   fChain->SetBranchAddress("MuDRHardestJet", MuDRHardestJet, &b_MuDRHardestJet);
   fChain->SetBranchAddress("MuCaloComp", MuCaloComp, &b_MuCaloComp);
   fChain->SetBranchAddress("MuSegmComp", MuSegmComp, &b_MuSegmComp);
   fChain->SetBranchAddress("MuOuterRad", MuOuterRad, &b_MuOuterRad);
   fChain->SetBranchAddress("MuNChi2", MuNChi2, &b_MuNChi2);
   fChain->SetBranchAddress("MuNTkHits", MuNTkHits, &b_MuNTkHits);
   fChain->SetBranchAddress("MuD0", MuD0, &b_MuD0);
   fChain->SetBranchAddress("MuDz", MuDz, &b_MuDz);
   fChain->SetBranchAddress("MuPtE", MuPtE, &b_MuPtE);
   fChain->SetBranchAddress("MuGenID", MuGenID, &b_MuGenID);
   fChain->SetBranchAddress("MuGenMoID", MuGenMoID, &b_MuGenMoID);
   fChain->SetBranchAddress("MuGenGMoID", MuGenGMoID, &b_MuGenGMoID);
   fChain->SetBranchAddress("MuGenType", MuGenType, &b_MuGenType);
   fChain->SetBranchAddress("MuGenMoType", MuGenMoType, &b_MuGenMoType);
   fChain->SetBranchAddress("MuGenGMoType", MuGenGMoType, &b_MuGenGMoType);
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
