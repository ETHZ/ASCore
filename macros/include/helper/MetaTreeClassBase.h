//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Apr  5 12:44:54 2011 by ROOT version 5.27/06b
// from TTree Analysis/AnalysisTree
// found on file: /scratch/stiegerb/SSDLTrees/SingleMu_Run2011A_AOD.root
//////////////////////////////////////////////////////////

#ifndef MetaTreeClassBase_h
#define MetaTreeClassBase_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

using namespace std;

class MetaTreeClassBase {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           Run;
   Int_t           Event;
   Int_t           LumiSec;
   Int_t           NHLTPaths;
   Int_t           HLTResults[214];   //[NHLTPaths]
   Int_t           HLTPrescales[214];   //[NHLTPaths]
   vector<string>  *HLTNames;
   Int_t           NMus;
   Float_t         MuPt[4];   //[NMus]
   Float_t         MuEta[4];   //[NMus]
   Float_t         MuPhi[4];   //[NMus]
   Int_t           MuCharge[4];   //[NMus]
   Int_t           MuTight[4];   //[NMus]
   Float_t         MuIso[4];   //[NMus]
   Float_t         MuIsoHybrid[4];   //[NMus]
   Float_t         MuD0[4];   //[NMus]
   Float_t         MuDz[4];   //[NMus]
   Float_t         MuD0BS[4];   //[NMus]
   Float_t         MuDzBS[4];   //[NMus]
   Float_t         MuPtE[4];   //[NMus]
   Int_t           MuGenID[4];   //[NMus]
   Int_t           MuGenMoID[4];   //[NMus]
   Int_t           MuGenGMoID[4];   //[NMus]
   Int_t           MuGenType[4];   //[NMus]
   Int_t           MuGenMoType[4];   //[NMus]
   Int_t           MuGenGMoType[4];   //[NMus]
   Float_t         MuMT[4];   //[NMus]
   Int_t           NEls;
   Int_t           ElCh[2];   //[NEls]
   Int_t           ElChIsCons[2];   //[NEls]
   Int_t           ElEcalDriven[2];   //[NEls]
   Float_t         ElCaloEnergy[2];   //[NEls]
   Float_t         ElPt[2];   //[NEls]
   Float_t         ElEta[2];   //[NEls]
   Float_t         ElPhi[2];   //[NEls]
   Float_t         ElD0[2];   //[NEls]
   Float_t         ElD0Err[2];   //[NEls]
   Float_t         ElDz[2];   //[NEls]
   Float_t         ElDzErr[2];   //[NEls]
   Float_t         ElEoverP[2];   //[NEls]
   Float_t         ElHoverE[2];   //[NEls]
   Float_t         ElSigmaIetaIeta[2];   //[NEls]
   Float_t         ElDeltaPhiSuperClusterAtVtx[2];   //[NEls]
   Float_t         ElDeltaEtaSuperClusterAtVtx[2];   //[NEls]
   Float_t         ElRelIso[2];   //[NEls]
   Int_t           ElIsGoodElId_WP80[2];   //[NEls]
   Int_t           ElIsGoodElId_WP90[2];   //[NEls]
   Int_t           ElIsGoodElId_WP95[2];   //[NEls]
   Float_t         ElS4OverS1[2];   //[NEls]
   Int_t           ElGenID[2];   //[NEls]
   Int_t           ElGenMID[2];   //[NEls]
   Int_t           ElGenGMID[2];   //[NEls]
   Int_t           ElGenType[2];   //[NEls]
   Int_t           ElGenMType[2];   //[NEls]
   Int_t           ElGenGMType[2];   //[NEls]
   Float_t         ElHybRelIso[2];   //[NEls]
   Float_t         ElMT[2];   //[NEls]
   Float_t         tcMET;
   Float_t         pfMET;
   Int_t           NJets;
   Float_t         JetPt[9];   //[NJets]
   Float_t         JetEta[9];   //[NJets]
   Float_t         JetPhi[9];   //[NJets]

   // List of branches
   TBranch        *b_Run;   //!
   TBranch        *b_Event;   //!
   TBranch        *b_LumiSec;   //!
   TBranch        *b_NHLTPaths;   //!
   TBranch        *b_HLTResults;   //!
   TBranch        *b_HLTPrescales;   //!
   TBranch        *b_HLTNames;   //!
   TBranch        *b_NMus;   //!
   TBranch        *b_MuPt;   //!
   TBranch        *b_MuEta;   //!
   TBranch        *b_MuPhi;   //!
   TBranch        *b_MuCharge;   //!
   TBranch        *b_MuTight;   //!
   TBranch        *b_MuIso;   //!
   TBranch        *b_MuIsoHybrid;   //!
   TBranch        *b_MuD0;   //!
   TBranch        *b_MuDz;   //!
   TBranch        *b_MuD0BS;   //!
   TBranch        *b_MuDzBS;   //!
   TBranch        *b_MuPtE;   //!
   TBranch        *b_MuGenID;   //!
   TBranch        *b_MuGenMoID;   //!
   TBranch        *b_MuGenGMoID;   //!
   TBranch        *b_MuGenType;   //!
   TBranch        *b_MuGenMoType;   //!
   TBranch        *b_MuGenGMoType;   //!
   TBranch        *b_MuMT;   //!
   TBranch        *b_NEls;   //!
   TBranch        *b_ElCh;   //!
   TBranch        *b_ElChIsCons;   //!
   TBranch        *b_ElEcalDriven;   //!
   TBranch        *b_ElCaloEnergy;   //!
   TBranch        *b_ElPt;   //!
   TBranch        *b_ElEta;   //!
   TBranch        *b_ElPhi;   //!
   TBranch        *b_ElD0;   //!
   TBranch        *b_ElD0Err;   //!
   TBranch        *b_ElDz;   //!
   TBranch        *b_ElDzErr;   //!
   TBranch        *b_ElEoverP;   //!
   TBranch        *b_ElHoverE;   //!
   TBranch        *b_ElSigmaIetaIeta;   //!
   TBranch        *b_ElDeltaPhiSuperClusterAtVtx;   //!
   TBranch        *b_ElDeltaEtaSuperClusterAtVtx;   //!
   TBranch        *b_ElRelIso;   //!
   TBranch        *b_ElIsGoodElId_WP80;   //!
   TBranch        *b_ElIsGoodElId_WP90;   //!
   TBranch        *b_ElIsGoodElId_WP95;   //!
   TBranch        *b_ElS4OverS1;   //!
   TBranch        *b_ElGenID;   //!
   TBranch        *b_ElGenMID;   //!
   TBranch        *b_ElGenGMID;   //!
   TBranch        *b_ElGenType;   //!
   TBranch        *b_ElGenMType;   //!
   TBranch        *b_ElGenGMType;   //!
   TBranch        *b_ElHybRelIso;   //!
   TBranch        *b_ElMT;   //!
   TBranch        *b_tcMET;   //!
   TBranch        *b_pfMET;   //!
   TBranch        *b_NJets;   //!
   TBranch        *b_JetPt;   //!
   TBranch        *b_JetEta;   //!
   TBranch        *b_JetPhi;   //!

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
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/scratch/stiegerb/SSDLTrees/SingleMu_Run2011A_AOD.root");
      if (!f) {
         f = new TFile("/scratch/stiegerb/SSDLTrees/SingleMu_Run2011A_AOD.root");
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

   // Set object pointer
   HLTNames = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("Run", &Run, &b_Run);
   fChain->SetBranchAddress("Event", &Event, &b_Event);
   fChain->SetBranchAddress("LumiSec", &LumiSec, &b_LumiSec);
   fChain->SetBranchAddress("NHLTPaths", &NHLTPaths, &b_NHLTPaths);
   fChain->SetBranchAddress("HLTResults", HLTResults, &b_HLTResults);
   fChain->SetBranchAddress("HLTPrescales", HLTPrescales, &b_HLTPrescales);
   fChain->SetBranchAddress("HLTNames", &HLTNames, &b_HLTNames);
   fChain->SetBranchAddress("NMus", &NMus, &b_NMus);
   fChain->SetBranchAddress("MuPt", MuPt, &b_MuPt);
   fChain->SetBranchAddress("MuEta", MuEta, &b_MuEta);
   fChain->SetBranchAddress("MuPhi", MuPhi, &b_MuPhi);
   fChain->SetBranchAddress("MuCharge", MuCharge, &b_MuCharge);
   fChain->SetBranchAddress("MuTight", MuTight, &b_MuTight);
   fChain->SetBranchAddress("MuIso", MuIso, &b_MuIso);
   fChain->SetBranchAddress("MuIsoHybrid", MuIsoHybrid, &b_MuIsoHybrid);
   fChain->SetBranchAddress("MuD0", MuD0, &b_MuD0);
   fChain->SetBranchAddress("MuDz", MuDz, &b_MuDz);
   fChain->SetBranchAddress("MuD0BS", MuD0BS, &b_MuD0BS);
   fChain->SetBranchAddress("MuDzBS", MuDzBS, &b_MuDzBS);
   fChain->SetBranchAddress("MuPtE", MuPtE, &b_MuPtE);
   fChain->SetBranchAddress("MuGenID", MuGenID, &b_MuGenID);
   fChain->SetBranchAddress("MuGenMoID", MuGenMoID, &b_MuGenMoID);
   fChain->SetBranchAddress("MuGenGMoID", MuGenGMoID, &b_MuGenGMoID);
   fChain->SetBranchAddress("MuGenType", MuGenType, &b_MuGenType);
   fChain->SetBranchAddress("MuGenMoType", MuGenMoType, &b_MuGenMoType);
   fChain->SetBranchAddress("MuGenGMoType", MuGenGMoType, &b_MuGenGMoType);
   fChain->SetBranchAddress("MuMT", MuMT, &b_MuMT);
   fChain->SetBranchAddress("NEls", &NEls, &b_NEls);
   fChain->SetBranchAddress("ElCh", ElCh, &b_ElCh);
   fChain->SetBranchAddress("ElChIsCons", ElChIsCons, &b_ElChIsCons);
   fChain->SetBranchAddress("ElEcalDriven", ElEcalDriven, &b_ElEcalDriven);
   fChain->SetBranchAddress("ElCaloEnergy", ElCaloEnergy, &b_ElCaloEnergy);
   fChain->SetBranchAddress("ElPt", ElPt, &b_ElPt);
   fChain->SetBranchAddress("ElEta", ElEta, &b_ElEta);
   fChain->SetBranchAddress("ElPhi", ElPhi, &b_ElPhi);
   fChain->SetBranchAddress("ElD0", ElD0, &b_ElD0);
   fChain->SetBranchAddress("ElD0Err", ElD0Err, &b_ElD0Err);
   fChain->SetBranchAddress("ElDz", ElDz, &b_ElDz);
   fChain->SetBranchAddress("ElDzErr", ElDzErr, &b_ElDzErr);
   fChain->SetBranchAddress("ElEoverP", ElEoverP, &b_ElEoverP);
   fChain->SetBranchAddress("ElHoverE", ElHoverE, &b_ElHoverE);
   fChain->SetBranchAddress("ElSigmaIetaIeta", ElSigmaIetaIeta, &b_ElSigmaIetaIeta);
   fChain->SetBranchAddress("ElDeltaPhiSuperClusterAtVtx", ElDeltaPhiSuperClusterAtVtx, &b_ElDeltaPhiSuperClusterAtVtx);
   fChain->SetBranchAddress("ElDeltaEtaSuperClusterAtVtx", ElDeltaEtaSuperClusterAtVtx, &b_ElDeltaEtaSuperClusterAtVtx);
   fChain->SetBranchAddress("ElRelIso", ElRelIso, &b_ElRelIso);
   fChain->SetBranchAddress("ElIsGoodElId_WP80", ElIsGoodElId_WP80, &b_ElIsGoodElId_WP80);
   fChain->SetBranchAddress("ElIsGoodElId_WP90", ElIsGoodElId_WP90, &b_ElIsGoodElId_WP90);
   fChain->SetBranchAddress("ElIsGoodElId_WP95", ElIsGoodElId_WP95, &b_ElIsGoodElId_WP95);
   fChain->SetBranchAddress("ElS4OverS1", ElS4OverS1, &b_ElS4OverS1);
   fChain->SetBranchAddress("ElGenID", ElGenID, &b_ElGenID);
   fChain->SetBranchAddress("ElGenMID", ElGenMID, &b_ElGenMID);
   fChain->SetBranchAddress("ElGenGMID", ElGenGMID, &b_ElGenGMID);
   fChain->SetBranchAddress("ElGenType", ElGenType, &b_ElGenType);
   fChain->SetBranchAddress("ElGenMType", ElGenMType, &b_ElGenMType);
   fChain->SetBranchAddress("ElGenGMType", ElGenGMType, &b_ElGenGMType);
   fChain->SetBranchAddress("ElHybRelIso", ElHybRelIso, &b_ElHybRelIso);
   fChain->SetBranchAddress("ElMT", ElMT, &b_ElMT);
   fChain->SetBranchAddress("tcMET", &tcMET, &b_tcMET);
   fChain->SetBranchAddress("pfMET", &pfMET, &b_pfMET);
   fChain->SetBranchAddress("NJets", &NJets, &b_NJets);
   fChain->SetBranchAddress("JetPt", JetPt, &b_JetPt);
   fChain->SetBranchAddress("JetEta", JetEta, &b_JetEta);
   fChain->SetBranchAddress("JetPhi", JetPhi, &b_JetPhi);
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
