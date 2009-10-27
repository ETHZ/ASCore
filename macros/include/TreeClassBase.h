//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Oct 20 15:37:53 2009 by ROOT version 5.22/00a
// from TTree Analysis/ETHZAnalysisTree
// found on file: NTuples/SD_Mu9/InclMu15.root
//////////////////////////////////////////////////////////

#ifndef TreeClassBase_h
#define TreeClassBase_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

class TreeClassBase {
	public :
	TTree          *fChain;   //!pointer to the analyzed TTree or TChain
	Int_t           fCurrent; //!current Tree number in a TChain

	// Declaration of leaf types
	Int_t           Run;
	Int_t           Event;
	Int_t           LumiSection;
	Int_t           SigProcID;
	Double_t        ExtXSecLO;
	Double_t        IntXSec;
	Double_t        Weight;
	Int_t           TrigResults[200];
	Double_t        PrimVtxx;
	Double_t        PrimVtxy;
	Double_t        PrimVtxz;
	Double_t        PrimVtxxE;
	Double_t        PrimVtxyE;
	Double_t        PrimVtxzE;
	Double_t        PrimVtxNChi2;
	Double_t        Beamspotx;
	Double_t        Beamspoty;
	Double_t        Beamspotz;
	Int_t           NMus;
	Double_t        MuPx[20];   //[NMus]
	Double_t        MuPy[20];   //[NMus]
	Double_t        MuPz[20];   //[NMus]
	Double_t        MuPt[20];   //[NMus]
	Double_t        MuE[20];   //[NMus]
	Double_t        MuEt[20];   //[NMus]
	Double_t        MuEta[20];   //[NMus]
	Double_t        MuPhi[20];   //[NMus]
	Int_t           MuCharge[20];   //[NMus]
	Double_t        MuPtsum[20];   //[NMus]
	Double_t        MuEtsum[20];   //[NMus]
	Double_t        MuIso[20];   //[NMus]
	Double_t        MuEem[20];   //[NMus]
	Double_t        MuEhad[20];   //[NMus]
	Double_t        MuD0BS[20];   //[NMus]
	Double_t        MuD0PV[20];   //[NMus]
	Double_t        MuD0E[20];   //[NMus]
	Double_t        MuDzBS[20];   //[NMus]
	Double_t        MuDzPV[20];   //[NMus]
	Double_t        MuDzE[20];   //[NMus]
	Double_t        MuNChi2[20];   //[NMus]
	Int_t           MuNGlHits[20];   //[NMus]
	Int_t           MuNMuHits[20];   //[NMus]
	Int_t           MuNTkHits[20];   //[NMus]
	Int_t           MuNMatches[20];   //[NMus]
	Int_t           MuNChambers[20];   //[NMus]
	Double_t        MuCaloComp[20];   //[NMus]
	Double_t        MuSegmComp[20];   //[NMus]
	Int_t           MuTrackerMu[20];   //[NMus]
	Int_t           MuGMPT[20];   //[NMus]
	Int_t           MuID[20];   //[NMus]
	Int_t           MuMID[20];   //[NMus]
	Int_t           NEles;
	Double_t        ElPx[20];   //[NEles]
	Double_t        ElPy[20];   //[NEles]
	Double_t        ElPz[20];   //[NEles]
	Double_t        ElPt[20];   //[NEles]
	Double_t        ElE[20];   //[NEles]
	Double_t        ElEt[20];   //[NEles]
	Double_t        ElEta[20];   //[NEles]
	Double_t        ElPhi[20];   //[NEles]
	Double_t        ElD0BS[20];   //[NEles]
	Double_t        ElD0PV[20];   //[NEles]
	Double_t        ElD0E[20];   //[NEles]
	Double_t        ElDzBS[20];   //[NEles]
	Double_t        ElDzPV[20];   //[NEles]
	Double_t        ElDzE[20];   //[NEles]
	Double_t        ElIso[20];   //[NEles]
	Double_t        ElPtSum[20];   //[NEles]
	Double_t        ElEtSum[20];   //[NEles]
	Double_t        ElNChi2[20];   //[NEles]
	Int_t           ElCharge[20];   //[NEles]
	Int_t           ElID[20][4];   //[NEles]
	Int_t           ElInGap[20];   //[NEles]
	Int_t           ElEcalDriven[20];   //[NEles]
	Int_t           ElTrackerDriven[20];   //[NEles]
	Int_t           ElBasicClustersSize[20];   //[NEles]
	Double_t        Elfbrem[20];   //[NEles]
	Double_t        ElHcalOverEcal[20];   //[NEles]
	Double_t        ElE5x5[20];   //[NEles]
	Double_t        ElE2x5Max[20];   //[NEles]
	Double_t        ElSigmaIetaIeta[20];   //[NEles]
	Double_t        ElDeltaPhiSeedClusterAtCalo[20];   //[NEles]
	Double_t        ElDeltaEtaSeedClusterAtCalo[20];   //[NEles]
	Double_t        ElDeltaPhiSuperClusterAtVtx[20];   //[NEles]
	Double_t        ElDeltaEtaSuperClusterAtVtx[20];   //[NEles]
	Double_t        ElCaloEnergy[20];   //[NEles]
	Double_t        ElTrkMomAtVtx[20];   //[NEles]
	Double_t        ElESuperClusterOverP[20];   //[NEles]
	Int_t           NJets;
	Double_t        JPx[20];   //[NJets]
	Double_t        JPy[20];   //[NJets]
	Double_t        JPz[20];   //[NJets]
	Double_t        JPt[20];   //[NJets]
	Double_t        JE[20];   //[NJets]
	Double_t        JEt[20];   //[NJets]
	Double_t        JEta[20];   //[NJets]
	Double_t        JPhi[20];   //[NJets]
	Double_t        JEMfrac[20];   //[NJets]
	Double_t        JID_HPD[20];   //[NJets]
	Double_t        JID_RBX[20];   //[NJets]
	Double_t        JID_n90Hits[20];   //[NJets]
	Double_t        JID_SubDet1[20];   //[NJets]
	Double_t        JID_SubDet2[20];   //[NJets]
	Double_t        JID_SubDet3[20];   //[NJets]
	Double_t        JID_SubDet4[20];   //[NJets]
	Double_t        JID_resEMF[20];   //[NJets]
	Double_t        JID_HCALTow[20];   //[NJets]
	Double_t        JID_ECALTow[20];   //[NJets]
	Double_t        JbTagProb[20];   //[NJets]
	Double_t        JChfrac[20];   //[NJets]
	Int_t           JNAssoTracks[20];   //[NJets]
	Double_t        JEcorr[20];   //[NJets]
	Double_t        JeMinDR[20];   //[NJets]
	Double_t        JVtxx[20];   //[NJets]
	Double_t        JVtxy[20];   //[NJets]
	Double_t        JVtxz[20];   //[NJets]
	Double_t        JVtxExx[20];   //[NJets]
	Double_t        JVtxEyx[20];   //[NJets]
	Double_t        JVtxEyy[20];   //[NJets]
	Double_t        JVtxEzy[20];   //[NJets]
	Double_t        JVtxEzz[20];   //[NJets]
	Double_t        JVtxEzx[20];   //[NJets]
	Double_t        JVtxNChi2[20];   //[NJets]
	Double_t        TrkPtSumx;
	Double_t        TrkPtSumy;
	Double_t        TrkPtSum;
	Double_t        TrkPtSumPhi;
	Double_t        ECALEsumx;
	Double_t        ECALEsumy;
	Double_t        ECALEsumz;
	Double_t        ECALMET;
	Double_t        ECALMETPhi;
	Double_t        HCALEsumx;
	Double_t        HCALEsumy;
	Double_t        HCALEsumz;
	Double_t        HCALMET;
	Double_t        HCALMETPhi;
	Double_t        RawMET;
	Double_t        RawMETpx;
	Double_t        RawMETpy;
	Double_t        RawMETphi;
	Double_t        MuCorrMET;
	Double_t        MuCorrMETpx;
	Double_t        MuCorrMETpy;
	Double_t        MuCorrMETphi;
	Double_t        TCMET;
	Double_t        TCMETpx;
	Double_t        TCMETpy;
	Double_t        TCMETphi;
	Double_t        MuJESCorrMET;
	Double_t        MuJESCorrMETpx;
	Double_t        MuJESCorrMETpy;
	Double_t        MuJESCorrMETphi;
	Double_t        PFMET;
	Double_t        PFMETpx;
	Double_t        PFMETpy;
	Double_t        PFMETphi;

	// List of branches
	TBranch        *b_Run;   //!
	TBranch        *b_Event;   //!
	TBranch        *b_LumiSection;   //!
	TBranch        *b_SigProcID;   //!
	TBranch        *b_ExtXSecLO;   //!
	TBranch        *b_IntXSec;   //!
	TBranch        *b_Weight;   //!
	TBranch        *b_TrigResults;   //!
	TBranch        *b_PrimVtxx;   //!
	TBranch        *b_PrimVtxy;   //!
	TBranch        *b_PrimVtxz;   //!
	TBranch        *b_PrimVtxxE;   //!
	TBranch        *b_PrimVtxyE;   //!
	TBranch        *b_PrimVtxzE;   //!
	TBranch        *b_PrimVtxNChi2;   //!
	TBranch        *b_Beamspotx;   //!
	TBranch        *b_Beamspoty;   //!
	TBranch        *b_Beamspotz;   //!
	TBranch        *b_NMus;   //!
	TBranch        *b_MuPx;   //!
	TBranch        *b_MuPy;   //!
	TBranch        *b_MuPz;   //!
	TBranch        *b_MuPt;   //!
	TBranch        *b_MuE;   //!
	TBranch        *b_MuEt;   //!
	TBranch        *b_MuEta;   //!
	TBranch        *b_MuPhi;   //!
	TBranch        *b_MuCharge;   //!
	TBranch        *b_MuPtsum;   //!
	TBranch        *b_MuEtsum;   //!
	TBranch        *b_MuIso;   //!
	TBranch        *b_MuEem;   //!
	TBranch        *b_MuEhad;   //!
	TBranch        *b_MuD0BS;   //!
	TBranch        *b_MuD0PV;   //!
	TBranch        *b_MuD0E;   //!
	TBranch        *b_MuDzBS;   //!
	TBranch        *b_MuDzPV;   //!
	TBranch        *b_MuDzE;   //!
	TBranch        *b_MuNChi2;   //!
	TBranch        *b_MuNGlHits;   //!
	TBranch        *b_MuNMuHits;   //!
	TBranch        *b_MuNTkHits;   //!
	TBranch        *b_MuNMatches;   //!
	TBranch        *b_MuNChambers;   //!
	TBranch        *b_MuCaloComp;   //!
	TBranch        *b_MuSegmComp;   //!
	TBranch        *b_MuTrackerMu;   //!
	TBranch        *b_MuGMPT;   //!
	TBranch        *b_MuID;   //!
	TBranch        *b_MuMID;   //!
	TBranch        *b_NEles;   //!
	TBranch        *b_ElPx;   //!
	TBranch        *b_ElPy;   //!
	TBranch        *b_ElPz;   //!
	TBranch        *b_ElPt;   //!
	TBranch        *b_ElE;   //!
	TBranch        *b_ElEt;   //!
	TBranch        *b_ElEta;   //!
	TBranch        *b_ElPhi;   //!
	TBranch        *b_ElD0BS;   //!
	TBranch        *b_ElD0PV;   //!
	TBranch        *b_ElD0E;   //!
	TBranch        *b_ElDzBS;   //!
	TBranch        *b_ElDzPV;   //!
	TBranch        *b_ElDzE;   //!
	TBranch        *b_ElIso;   //!
	TBranch        *b_ElPtSum;   //!
	TBranch        *b_ElEtSum;   //!
	TBranch        *b_ElNChi2;   //!
	TBranch        *b_ElCharge;   //!
	TBranch        *b_ElID;   //!
	TBranch        *b_ElInGap;   //!
	TBranch        *b_ElEcalDriven;   //!
	TBranch        *b_ElTrackerDriven;   //!
	TBranch        *b_ElBasicClustersSize;   //!
	TBranch        *b_Elfbrem;   //!
	TBranch        *b_ElHcalOverEcal;   //!
	TBranch        *b_ElE5x5;   //!
	TBranch        *b_ElE2x5Max;   //!
	TBranch        *b_ElSigmaIetaIeta;   //!
	TBranch        *b_ElDeltaPhiSeedClusterAtCalo;   //!
	TBranch        *b_ElDeltaEtaSeedClusterAtCalo;   //!
	TBranch        *b_ElDeltaPhiSuperClusterAtVtx;   //!
	TBranch        *b_ElDeltaEtaSuperClusterAtVtx;   //!
	TBranch        *b_ElCaloEnergy;   //!
	TBranch        *b_ElTrkMomAtVtx;   //!
	TBranch        *b_ElESuperClusterOverP;   //!
	TBranch        *b_NJets;   //!
	TBranch        *b_JPx;   //!
	TBranch        *b_JPy;   //!
	TBranch        *b_JPz;   //!
	TBranch        *b_JPt;   //!
	TBranch        *b_JE;   //!
	TBranch        *b_JEt;   //!
	TBranch        *b_JEta;   //!
	TBranch        *b_JPhi;   //!
	TBranch        *b_JEMfrac;   //!
	TBranch        *b_JID_HPD;   //!
	TBranch        *b_JID_RBX;   //!
	TBranch        *b_JID_n90Hits;   //!
	TBranch        *b_JID_SubDet1;   //!
	TBranch        *b_JID_SubDet2;   //!
	TBranch        *b_JID_SubDet3;   //!
	TBranch        *b_JID_SubDet4;   //!
	TBranch        *b_JID_resEMF;   //!
	TBranch        *b_JID_HCALTow;   //!
	TBranch        *b_JID_ECALTow;   //!
	TBranch        *b_JbTagProb;   //!
	TBranch        *b_JChfrac;   //!
	TBranch        *b_JNAssoTracks;   //!
	TBranch        *b_JEcorr;   //!
	TBranch        *b_JeMinDR;   //!
	TBranch        *b_JVtxx;   //!
	TBranch        *b_JVtxy;   //!
	TBranch        *b_JVtxz;   //!
	TBranch        *b_JVtxExx;   //!
	TBranch        *b_JVtxEyx;   //!
	TBranch        *b_JVtxEyy;   //!
	TBranch        *b_JVtxEzy;   //!
	TBranch        *b_JVtxEzz;   //!
	TBranch        *b_JVtxEzx;   //!
	TBranch        *b_JVtxNChi2;   //!
	TBranch        *b_TrkPtSumx;   //!
	TBranch        *b_TrkPtSumy;   //!
	TBranch        *b_TrkPtSum;   //!
	TBranch        *b_TrkPtSumPhi;   //!
	TBranch        *b_ECALEsumx;   //!
	TBranch        *b_ECALEsumy;   //!
	TBranch        *b_ECALEsumz;   //!
	TBranch        *b_ECALMET;   //!
	TBranch        *b_ECALMETPhi;   //!
	TBranch        *b_HCALEsumx;   //!
	TBranch        *b_HCALEsumy;   //!
	TBranch        *b_HCALEsumz;   //!
	TBranch        *b_HCALMET;   //!
	TBranch        *b_HCALMETPhi;   //!
	TBranch        *b_RawMET;   //!
	TBranch        *b_RawMETpx;   //!
	TBranch        *b_RawMETpy;   //!
	TBranch        *b_RawMETphi;   //!
	TBranch        *b_MuCorrMET;   //!
	TBranch        *b_MuCorrMETpx;   //!
	TBranch        *b_MuCorrMETpy;   //!
	TBranch        *b_MuCorrMETphi;   //!
	TBranch        *b_TCMET;   //!
	TBranch        *b_TCMETpx;   //!
	TBranch        *b_TCMETpy;   //!
	TBranch        *b_TCMETphi;   //!
	TBranch        *b_MuJESCorrMET;   //!
	TBranch        *b_MuJESCorrMETpx;   //!
	TBranch        *b_MuJESCorrMETpy;   //!
	TBranch        *b_MuJESCorrMETphi;   //!
	TBranch        *b_PFMET;   //!
	TBranch        *b_PFMETpx;   //!
	TBranch        *b_PFMETpy;   //!
	TBranch        *b_PFMETphi;   //!

	TreeClassBase(TTree *tree=0);
	virtual ~TreeClassBase();
	virtual Int_t    Cut(Long64_t entry);
	virtual Int_t    GetEntry(Long64_t entry);
	virtual Long64_t LoadTree(Long64_t entry);
	virtual void     Init(TTree *tree);
	virtual void     Loop();
	virtual Bool_t   Notify();
	virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef TreeClassBase_cxx
TreeClassBase::TreeClassBase(TTree *tree)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
	if (tree == 0) {
		TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("NTuples/SD_Mu9/InclMu15.root");
		if (!f) {
			f = new TFile("NTuples/SD_Mu9/InclMu15.root");
			f->cd("NTuples/SD_Mu9/InclMu15.root:/analyze");
		}
		tree = (TTree*)gDirectory->Get("Analysis");

	}
	Init(tree);
}

TreeClassBase::~TreeClassBase()
{
	if (!fChain) return;
	delete fChain->GetCurrentFile();
}

Int_t TreeClassBase::GetEntry(Long64_t entry)
{
// Read contents of entry.
	if (!fChain) return 0;
	return fChain->GetEntry(entry);
}

Long64_t TreeClassBase::LoadTree(Long64_t entry)
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

void TreeClassBase::Init(TTree *tree)
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
	fChain->SetBranchAddress("SigProcID", &SigProcID, &b_SigProcID);
	fChain->SetBranchAddress("ExtXSecLO", &ExtXSecLO, &b_ExtXSecLO);
	fChain->SetBranchAddress("IntXSec", &IntXSec, &b_IntXSec);
	fChain->SetBranchAddress("Weight", &Weight, &b_Weight);
	fChain->SetBranchAddress("TrigResults", TrigResults, &b_TrigResults);
	fChain->SetBranchAddress("PrimVtxx", &PrimVtxx, &b_PrimVtxx);
	fChain->SetBranchAddress("PrimVtxy", &PrimVtxy, &b_PrimVtxy);
	fChain->SetBranchAddress("PrimVtxz", &PrimVtxz, &b_PrimVtxz);
	fChain->SetBranchAddress("PrimVtxxE", &PrimVtxxE, &b_PrimVtxxE);
	fChain->SetBranchAddress("PrimVtxyE", &PrimVtxyE, &b_PrimVtxyE);
	fChain->SetBranchAddress("PrimVtxzE", &PrimVtxzE, &b_PrimVtxzE);
	fChain->SetBranchAddress("PrimVtxNChi2", &PrimVtxNChi2, &b_PrimVtxNChi2);
	fChain->SetBranchAddress("Beamspotx", &Beamspotx, &b_Beamspotx);
	fChain->SetBranchAddress("Beamspoty", &Beamspoty, &b_Beamspoty);
	fChain->SetBranchAddress("Beamspotz", &Beamspotz, &b_Beamspotz);
	fChain->SetBranchAddress("NMus", &NMus, &b_NMus);
	fChain->SetBranchAddress("MuPx", MuPx, &b_MuPx);
	fChain->SetBranchAddress("MuPy", MuPy, &b_MuPy);
	fChain->SetBranchAddress("MuPz", MuPz, &b_MuPz);
	fChain->SetBranchAddress("MuPt", MuPt, &b_MuPt);
	fChain->SetBranchAddress("MuE", MuE, &b_MuE);
	fChain->SetBranchAddress("MuEt", MuEt, &b_MuEt);
	fChain->SetBranchAddress("MuEta", MuEta, &b_MuEta);
	fChain->SetBranchAddress("MuPhi", MuPhi, &b_MuPhi);
	fChain->SetBranchAddress("MuCharge", MuCharge, &b_MuCharge);
	fChain->SetBranchAddress("MuPtsum", MuPtsum, &b_MuPtsum);
	fChain->SetBranchAddress("MuEtsum", MuEtsum, &b_MuEtsum);
	fChain->SetBranchAddress("MuIso", MuIso, &b_MuIso);
	fChain->SetBranchAddress("MuEem", MuEem, &b_MuEem);
	fChain->SetBranchAddress("MuEhad", MuEhad, &b_MuEhad);
	fChain->SetBranchAddress("MuD0BS", MuD0BS, &b_MuD0BS);
	fChain->SetBranchAddress("MuD0PV", MuD0PV, &b_MuD0PV);
	fChain->SetBranchAddress("MuD0E", MuD0E, &b_MuD0E);
	fChain->SetBranchAddress("MuDzBS", MuDzBS, &b_MuDzBS);
	fChain->SetBranchAddress("MuDzPV", MuDzPV, &b_MuDzPV);
	fChain->SetBranchAddress("MuDzE", MuDzE, &b_MuDzE);
	fChain->SetBranchAddress("MuNChi2", MuNChi2, &b_MuNChi2);
	fChain->SetBranchAddress("MuNGlHits", MuNGlHits, &b_MuNGlHits);
	fChain->SetBranchAddress("MuNMuHits", MuNMuHits, &b_MuNMuHits);
	fChain->SetBranchAddress("MuNTkHits", MuNTkHits, &b_MuNTkHits);
	fChain->SetBranchAddress("MuNMatches", MuNMatches, &b_MuNMatches);
	fChain->SetBranchAddress("MuNChambers", MuNChambers, &b_MuNChambers);
	fChain->SetBranchAddress("MuCaloComp", MuCaloComp, &b_MuCaloComp);
	fChain->SetBranchAddress("MuSegmComp", MuSegmComp, &b_MuSegmComp);
	fChain->SetBranchAddress("MuTrackerMu", MuTrackerMu, &b_MuTrackerMu);
	fChain->SetBranchAddress("MuGMPT", MuGMPT, &b_MuGMPT);
	fChain->SetBranchAddress("MuID", MuID, &b_MuID);
	fChain->SetBranchAddress("MuMID", MuMID, &b_MuMID);
	fChain->SetBranchAddress("NEles", &NEles, &b_NEles);
	fChain->SetBranchAddress("ElPx", ElPx, &b_ElPx);
	fChain->SetBranchAddress("ElPy", ElPy, &b_ElPy);
	fChain->SetBranchAddress("ElPz", ElPz, &b_ElPz);
	fChain->SetBranchAddress("ElPt", ElPt, &b_ElPt);
	fChain->SetBranchAddress("ElE", ElE, &b_ElE);
	fChain->SetBranchAddress("ElEt", ElEt, &b_ElEt);
	fChain->SetBranchAddress("ElEta", ElEta, &b_ElEta);
	fChain->SetBranchAddress("ElPhi", ElPhi, &b_ElPhi);
	fChain->SetBranchAddress("ElD0BS", ElD0BS, &b_ElD0BS);
	fChain->SetBranchAddress("ElD0PV", ElD0PV, &b_ElD0PV);
	fChain->SetBranchAddress("ElD0E", ElD0E, &b_ElD0E);
	fChain->SetBranchAddress("ElDzBS", ElDzBS, &b_ElDzBS);
	fChain->SetBranchAddress("ElDzPV", ElDzPV, &b_ElDzPV);
	fChain->SetBranchAddress("ElDzE", ElDzE, &b_ElDzE);
	fChain->SetBranchAddress("ElIso", ElIso, &b_ElIso);
	fChain->SetBranchAddress("ElPtSum", ElPtSum, &b_ElPtSum);
	fChain->SetBranchAddress("ElEtSum", ElEtSum, &b_ElEtSum);
	fChain->SetBranchAddress("ElNChi2", ElNChi2, &b_ElNChi2);
	fChain->SetBranchAddress("ElCharge", ElCharge, &b_ElCharge);
	fChain->SetBranchAddress("ElID", ElID, &b_ElID);
	fChain->SetBranchAddress("ElInGap", ElInGap, &b_ElInGap);
	fChain->SetBranchAddress("ElEcalDriven", ElEcalDriven, &b_ElEcalDriven);
	fChain->SetBranchAddress("ElTrackerDriven", ElTrackerDriven, &b_ElTrackerDriven);
	fChain->SetBranchAddress("ElBasicClustersSize", ElBasicClustersSize, &b_ElBasicClustersSize);
	fChain->SetBranchAddress("Elfbrem", Elfbrem, &b_Elfbrem);
	fChain->SetBranchAddress("ElHcalOverEcal", ElHcalOverEcal, &b_ElHcalOverEcal);
	fChain->SetBranchAddress("ElE5x5", ElE5x5, &b_ElE5x5);
	fChain->SetBranchAddress("ElE2x5Max", ElE2x5Max, &b_ElE2x5Max);
	fChain->SetBranchAddress("ElSigmaIetaIeta", ElSigmaIetaIeta, &b_ElSigmaIetaIeta);
	fChain->SetBranchAddress("ElDeltaPhiSeedClusterAtCalo", ElDeltaPhiSeedClusterAtCalo, &b_ElDeltaPhiSeedClusterAtCalo);
	fChain->SetBranchAddress("ElDeltaEtaSeedClusterAtCalo", ElDeltaEtaSeedClusterAtCalo, &b_ElDeltaEtaSeedClusterAtCalo);
	fChain->SetBranchAddress("ElDeltaPhiSuperClusterAtVtx", ElDeltaPhiSuperClusterAtVtx, &b_ElDeltaPhiSuperClusterAtVtx);
	fChain->SetBranchAddress("ElDeltaEtaSuperClusterAtVtx", ElDeltaEtaSuperClusterAtVtx, &b_ElDeltaEtaSuperClusterAtVtx);
	fChain->SetBranchAddress("ElCaloEnergy", ElCaloEnergy, &b_ElCaloEnergy);
	fChain->SetBranchAddress("ElTrkMomAtVtx", ElTrkMomAtVtx, &b_ElTrkMomAtVtx);
	fChain->SetBranchAddress("ElESuperClusterOverP", ElESuperClusterOverP, &b_ElESuperClusterOverP);
	fChain->SetBranchAddress("NJets", &NJets, &b_NJets);
	fChain->SetBranchAddress("JPx", JPx, &b_JPx);
	fChain->SetBranchAddress("JPy", JPy, &b_JPy);
	fChain->SetBranchAddress("JPz", JPz, &b_JPz);
	fChain->SetBranchAddress("JPt", JPt, &b_JPt);
	fChain->SetBranchAddress("JE", JE, &b_JE);
	fChain->SetBranchAddress("JEt", JEt, &b_JEt);
	fChain->SetBranchAddress("JEta", JEta, &b_JEta);
	fChain->SetBranchAddress("JPhi", JPhi, &b_JPhi);
	fChain->SetBranchAddress("JEMfrac", JEMfrac, &b_JEMfrac);
	fChain->SetBranchAddress("JID_HPD", JID_HPD, &b_JID_HPD);
	fChain->SetBranchAddress("JID_RBX", JID_RBX, &b_JID_RBX);
	fChain->SetBranchAddress("JID_n90Hits", JID_n90Hits, &b_JID_n90Hits);
	fChain->SetBranchAddress("JID_SubDet1", JID_SubDet1, &b_JID_SubDet1);
	fChain->SetBranchAddress("JID_SubDet2", JID_SubDet2, &b_JID_SubDet2);
	fChain->SetBranchAddress("JID_SubDet3", JID_SubDet3, &b_JID_SubDet3);
	fChain->SetBranchAddress("JID_SubDet4", JID_SubDet4, &b_JID_SubDet4);
	fChain->SetBranchAddress("JID_resEMF", JID_resEMF, &b_JID_resEMF);
	fChain->SetBranchAddress("JID_HCALTow", JID_HCALTow, &b_JID_HCALTow);
	fChain->SetBranchAddress("JID_ECALTow", JID_ECALTow, &b_JID_ECALTow);
	fChain->SetBranchAddress("JbTagProb", JbTagProb, &b_JbTagProb);
	fChain->SetBranchAddress("JChfrac", JChfrac, &b_JChfrac);
	fChain->SetBranchAddress("JNAssoTracks", JNAssoTracks, &b_JNAssoTracks);
	fChain->SetBranchAddress("JEcorr", JEcorr, &b_JEcorr);
	fChain->SetBranchAddress("JeMinDR", JeMinDR, &b_JeMinDR);
	fChain->SetBranchAddress("JVtxx", JVtxx, &b_JVtxx);
	fChain->SetBranchAddress("JVtxy", JVtxy, &b_JVtxy);
	fChain->SetBranchAddress("JVtxz", JVtxz, &b_JVtxz);
	fChain->SetBranchAddress("JVtxExx", JVtxExx, &b_JVtxExx);
	fChain->SetBranchAddress("JVtxEyx", JVtxEyx, &b_JVtxEyx);
	fChain->SetBranchAddress("JVtxEyy", JVtxEyy, &b_JVtxEyy);
	fChain->SetBranchAddress("JVtxEzy", JVtxEzy, &b_JVtxEzy);
	fChain->SetBranchAddress("JVtxEzz", JVtxEzz, &b_JVtxEzz);
	fChain->SetBranchAddress("JVtxEzx", JVtxEzx, &b_JVtxEzx);
	fChain->SetBranchAddress("JVtxNChi2", JVtxNChi2, &b_JVtxNChi2);
	fChain->SetBranchAddress("TrkPtSumx", &TrkPtSumx, &b_TrkPtSumx);
	fChain->SetBranchAddress("TrkPtSumy", &TrkPtSumy, &b_TrkPtSumy);
	fChain->SetBranchAddress("TrkPtSum", &TrkPtSum, &b_TrkPtSum);
	fChain->SetBranchAddress("TrkPtSumPhi", &TrkPtSumPhi, &b_TrkPtSumPhi);
	fChain->SetBranchAddress("ECALEsumx", &ECALEsumx, &b_ECALEsumx);
	fChain->SetBranchAddress("ECALEsumy", &ECALEsumy, &b_ECALEsumy);
	fChain->SetBranchAddress("ECALEsumz", &ECALEsumz, &b_ECALEsumz);
	fChain->SetBranchAddress("ECALMET", &ECALMET, &b_ECALMET);
	fChain->SetBranchAddress("ECALMETPhi", &ECALMETPhi, &b_ECALMETPhi);
	fChain->SetBranchAddress("HCALEsumx", &HCALEsumx, &b_HCALEsumx);
	fChain->SetBranchAddress("HCALEsumy", &HCALEsumy, &b_HCALEsumy);
	fChain->SetBranchAddress("HCALEsumz", &HCALEsumz, &b_HCALEsumz);
	fChain->SetBranchAddress("HCALMET", &HCALMET, &b_HCALMET);
	fChain->SetBranchAddress("HCALMETPhi", &HCALMETPhi, &b_HCALMETPhi);
	fChain->SetBranchAddress("RawMET", &RawMET, &b_RawMET);
	fChain->SetBranchAddress("RawMETpx", &RawMETpx, &b_RawMETpx);
	fChain->SetBranchAddress("RawMETpy", &RawMETpy, &b_RawMETpy);
	fChain->SetBranchAddress("RawMETphi", &RawMETphi, &b_RawMETphi);
	fChain->SetBranchAddress("MuCorrMET", &MuCorrMET, &b_MuCorrMET);
	fChain->SetBranchAddress("MuCorrMETpx", &MuCorrMETpx, &b_MuCorrMETpx);
	fChain->SetBranchAddress("MuCorrMETpy", &MuCorrMETpy, &b_MuCorrMETpy);
	fChain->SetBranchAddress("MuCorrMETphi", &MuCorrMETphi, &b_MuCorrMETphi);
	fChain->SetBranchAddress("TCMET", &TCMET, &b_TCMET);
	fChain->SetBranchAddress("TCMETpx", &TCMETpx, &b_TCMETpx);
	fChain->SetBranchAddress("TCMETpy", &TCMETpy, &b_TCMETpy);
	fChain->SetBranchAddress("TCMETphi", &TCMETphi, &b_TCMETphi);
	fChain->SetBranchAddress("MuJESCorrMET", &MuJESCorrMET, &b_MuJESCorrMET);
	fChain->SetBranchAddress("MuJESCorrMETpx", &MuJESCorrMETpx, &b_MuJESCorrMETpx);
	fChain->SetBranchAddress("MuJESCorrMETpy", &MuJESCorrMETpy, &b_MuJESCorrMETpy);
	fChain->SetBranchAddress("MuJESCorrMETphi", &MuJESCorrMETphi, &b_MuJESCorrMETphi);
	fChain->SetBranchAddress("PFMET", &PFMET, &b_PFMET);
	fChain->SetBranchAddress("PFMETpx", &PFMETpx, &b_PFMETpx);
	fChain->SetBranchAddress("PFMETpy", &PFMETpy, &b_PFMETpy);
	fChain->SetBranchAddress("PFMETphi", &PFMETphi, &b_PFMETphi);
	Notify();
}

Bool_t TreeClassBase::Notify()
{
	// The Notify() function is called when a new file is opened. This
	// can be either for a new TTree in a TChain or when when a new TTree
	// is started when using PROOF. It is normally not necessary to make changes
	// to the generated code, but the routine can be extended by the
	// user if needed. The return value is currently not used.

	return kTRUE;
}

void TreeClassBase::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
	if (!fChain) return;
	fChain->Show(entry);
}
Int_t TreeClassBase::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
	return 1;
}
#endif // #ifdef TreeClassBase_cxx
