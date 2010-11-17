//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Nov 16 10:07:37 2010 by ROOT version 5.22/00d
// from TTree MassAnalysis/MassAnalysisTree
// found on file: /shome/pnef/SUSY/SUSY_macros/analyzed/RunLeptJetMultAnalyzer/20101115_Jet_RunB_V01-08-00/Jet_RunB_PromptReco_V01-08-00/Mass_histos.root
//////////////////////////////////////////////////////////

#ifndef MassAnalysisTreeClass_h
#define MassAnalysisTreeClass_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

class MassAnalysisTreeClass {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           Run;
   Int_t           Event;
   Int_t           LumiSection;
   Int_t           HCALNoiseFlag;
   Float_t         Weight;
   Int_t           NJets;
   Double_t        JPt[11];   //[NJets]
   Double_t        DPhiJetsMet[11];   //[NJets]
   Int_t           NElecs;
   Int_t           NMuons;
   Int_t           LeptConfig;
   Int_t           IsCleanMultiJetEvent;
   Int_t           IsCleanJetEvent;
   Int_t           R12R21;
   Int_t           NJetsPt50Eta25;
   Double_t        PseudoJetMT2;
   Double_t        PseudoJetMT2massive;
   Double_t        PseudoJetMT2simple;
   Double_t        PseudoJetMCT;
   Double_t        PseudoJetMCTMinusMass;
   Double_t        PseudoJetMCTMassless;
   Double_t        PseudoJet1Pt;
   Double_t        PseudoJet2Pt;
   Double_t        PseudojetAlphaT;
   Double_t        Vectorsumpt;
   Double_t        PFMET;
   Double_t        PFMETsign;
   Double_t        LeadingJetEta;
   Double_t        DPhiJ1Met;
   Double_t        DPhiJ2Met;
   Double_t        PseudoJetMT2AxisdPhi;
   Double_t        R1221min;
   Double_t        MPT_sel;
   Double_t        MPT;
   Double_t        MHT;
   Double_t        HT;
   Double_t        DPhiMhtMpt;

   // List of branches
   TBranch        *b_Run;   //!
   TBranch        *b_Event;   //!
   TBranch        *b_LumiSection;   //!
   TBranch        *b_HCALNoiseFlag;   //!
   TBranch        *b_Weight;   //!
   TBranch        *b_NJets;   //!
   TBranch        *b_JPt;   //!
   TBranch        *b_DPhiJetsMet;   //!
   TBranch        *b_NElecs;   //!
   TBranch        *b_NMuons;   //!
   TBranch        *b_LeptConfig;   //!
   TBranch        *b_IsCleanMultiJetEvent;   //!
   TBranch        *b_IsCleanJetEvent;   //!
   TBranch        *b_R12R21;   //!
   TBranch        *b_NJetsPt50Eta25;   //!
   TBranch        *b_PseudoJetMT2;   //!
   TBranch        *b_PseudoJetMT2massive;   //!
   TBranch        *b_PseudoJetMT2simple;   //!
   TBranch        *b_PseudoJetMCT;   //!
   TBranch        *b_PseudoJetMCTMinusMass;   //!
   TBranch        *b_PseudoJetMCTMassless;   //!
   TBranch        *b_PseudoJet1Pt;   //!
   TBranch        *b_PseudoJet2Pt;   //!
   TBranch        *b_PseudojetAlphaT;   //!
   TBranch        *b_Vectorsumpt;   //!
   TBranch        *b_PFMET;   //!
   TBranch        *b_PFMETsign;   //!
   TBranch        *b_LeadingJetEta;   //!
   TBranch        *b_DPhiJ1Met;   //!
   TBranch        *b_DPhiJ2Met;   //!
   TBranch        *b_PseudoJetMT2AxisdPhi;   //!
   TBranch        *b_R1221min;   //!
   TBranch        *b_MPT_sel;   //!
   TBranch        *b_MPT;   //!
   TBranch        *b_MHT;   //!
   TBranch        *b_HT;   //!
   TBranch        *b_DPhiMhtMpt;   //!

   MassAnalysisTreeClass(TTree *tree=0);
   virtual ~MassAnalysisTreeClass();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef MassAnalysisTreeClass_cxx
MassAnalysisTreeClass::MassAnalysisTreeClass(TTree *tree)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/shome/pnef/SUSY/SUSY_macros/analyzed/RunLeptJetMultAnalyzer/20101115_Jet_RunB_V01-08-00/Jet_RunB_PromptReco_V01-08-00/Mass_histos.root");
      if (!f) {
         f = new TFile("/shome/pnef/SUSY/SUSY_macros/analyzed/RunLeptJetMultAnalyzer/20101115_Jet_RunB_V01-08-00/Jet_RunB_PromptReco_V01-08-00/Mass_histos.root");
      }
      tree = (TTree*)gDirectory->Get("MassAnalysis");

   }
   Init(tree);
}

MassAnalysisTreeClass::~MassAnalysisTreeClass()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t MassAnalysisTreeClass::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t MassAnalysisTreeClass::LoadTree(Long64_t entry)
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

void MassAnalysisTreeClass::Init(TTree *tree)
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
   fChain->SetBranchAddress("HCALNoiseFlag", &HCALNoiseFlag, &b_HCALNoiseFlag);
   fChain->SetBranchAddress("Weight", &Weight, &b_Weight);
   fChain->SetBranchAddress("NJets", &NJets, &b_NJets);
   fChain->SetBranchAddress("JPt", JPt, &b_JPt);
   fChain->SetBranchAddress("DPhiJetsMet", DPhiJetsMet, &b_DPhiJetsMet);
   fChain->SetBranchAddress("NElecs", &NElecs, &b_NElecs);
   fChain->SetBranchAddress("NMuons", &NMuons, &b_NMuons);
   fChain->SetBranchAddress("LeptConfig", &LeptConfig, &b_LeptConfig);
   fChain->SetBranchAddress("IsCleanMultiJetEvent", &IsCleanMultiJetEvent, &b_IsCleanMultiJetEvent);
   fChain->SetBranchAddress("IsCleanJetEvent", &IsCleanJetEvent, &b_IsCleanJetEvent);
   fChain->SetBranchAddress("R12R21", &R12R21, &b_R12R21);
   fChain->SetBranchAddress("NJetsPt50Eta25", &NJetsPt50Eta25, &b_NJetsPt50Eta25);
   fChain->SetBranchAddress("PseudoJetMT2", &PseudoJetMT2, &b_PseudoJetMT2);
   fChain->SetBranchAddress("PseudoJetMT2massive", &PseudoJetMT2massive, &b_PseudoJetMT2massive);
   fChain->SetBranchAddress("PseudoJetMT2simple", &PseudoJetMT2simple, &b_PseudoJetMT2simple);
   fChain->SetBranchAddress("PseudoJetMCT", &PseudoJetMCT, &b_PseudoJetMCT);
   fChain->SetBranchAddress("PseudoJetMCTMinusMass", &PseudoJetMCTMinusMass, &b_PseudoJetMCTMinusMass);
   fChain->SetBranchAddress("PseudoJetMCTMassless", &PseudoJetMCTMassless, &b_PseudoJetMCTMassless);
   fChain->SetBranchAddress("PseudoJet1Pt", &PseudoJet1Pt, &b_PseudoJet1Pt);
   fChain->SetBranchAddress("PseudoJet2Pt", &PseudoJet2Pt, &b_PseudoJet2Pt);
   fChain->SetBranchAddress("PseudojetAlphaT", &PseudojetAlphaT, &b_PseudojetAlphaT);
   fChain->SetBranchAddress("Vectorsumpt", &Vectorsumpt, &b_Vectorsumpt);
   fChain->SetBranchAddress("PFMET", &PFMET, &b_PFMET);
   fChain->SetBranchAddress("PFMETsign", &PFMETsign, &b_PFMETsign);
   fChain->SetBranchAddress("LeadingJetEta", &LeadingJetEta, &b_LeadingJetEta);
   fChain->SetBranchAddress("DPhiJ1Met", &DPhiJ1Met, &b_DPhiJ1Met);
   fChain->SetBranchAddress("DPhiJ2Met", &DPhiJ2Met, &b_DPhiJ2Met);
   fChain->SetBranchAddress("PseudoJetMT2AxisdPhi", &PseudoJetMT2AxisdPhi, &b_PseudoJetMT2AxisdPhi);
   fChain->SetBranchAddress("R1221min", &R1221min, &b_R1221min);
   fChain->SetBranchAddress("MPT_sel", &MPT_sel, &b_MPT_sel);
   fChain->SetBranchAddress("MPT", &MPT, &b_MPT);
   fChain->SetBranchAddress("MHT", &MHT, &b_MHT);
   fChain->SetBranchAddress("HT", &HT, &b_HT);
   fChain->SetBranchAddress("DPhiMhtMpt", &DPhiMhtMpt, &b_DPhiMhtMpt);
   Notify();
}

Bool_t MassAnalysisTreeClass::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void MassAnalysisTreeClass::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t MassAnalysisTreeClass::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef MassAnalysisTreeClass_cxx
