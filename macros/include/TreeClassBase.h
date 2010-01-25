//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Jan 25 19:05:48 2010 by ROOT version 5.22/00d
// from TTree Analysis/ETHZAnalysisTree
// found on file: TTBar4Jets_NTP.root
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
   Int_t           HLTResults[200];
   Int_t           L1PhysResults[128];
   Int_t           L1TechResults[64];
   Int_t           PrimVtxGood;
   Double_t        PrimVtxx;
   Double_t        PrimVtxy;
   Double_t        PrimVtxz;
   Double_t        PrimVtxxE;
   Double_t        PrimVtxyE;
   Double_t        PrimVtxzE;
   Double_t        PrimVtxNChi2;
   Int_t           PrimVtxNTracks;
   Double_t        PrimVtxPtSum;
   Double_t        Beamspotx;
   Double_t        Beamspoty;
   Double_t        Beamspotz;
   Int_t           NCaloTowers;
   Int_t           GoodEvent;
   Int_t           MaxMuExceed;
   Int_t           MaxElExceed;
   Int_t           MaxJetExceed;
   Int_t           MaxUncJetExceed;
   Int_t           MaxTrkExceed;
   Int_t           MaxPhotonsExceed;
   Int_t           NMus;
   Int_t           MuGood[5];   //[NMus]
   Int_t           MuIsIso[5];   //[NMus]
   Double_t        MuPx[5];   //[NMus]
   Double_t        MuPy[5];   //[NMus]
   Double_t        MuPz[5];   //[NMus]
   Double_t        MuPt[5];   //[NMus]
   Double_t        MuPtE[5];   //[NMus]
   Double_t        MuE[5];   //[NMus]
   Double_t        MuEt[5];   //[NMus]
   Double_t        MuEta[5];   //[NMus]
   Double_t        MuPhi[5];   //[NMus]
   Int_t           MuCharge[5];   //[NMus]
   Double_t        MuPtsum[5];   //[NMus]
   Double_t        MuEtsum[5];   //[NMus]
   Double_t        MuIso[5];   //[NMus]
   Double_t        MuEem[5];   //[NMus]
   Double_t        MuEhad[5];   //[NMus]
   Double_t        MuD0BS[5];   //[NMus]
   Double_t        MuD0PV[5];   //[NMus]
   Double_t        MuD0E[5];   //[NMus]
   Double_t        MuDzBS[5];   //[NMus]
   Double_t        MuDzPV[5];   //[NMus]
   Double_t        MuDzE[5];   //[NMus]
   Double_t        MuNChi2[5];   //[NMus]
   Int_t           MuNGlHits[5];   //[NMus]
   Int_t           MuNMuHits[5];   //[NMus]
   Int_t           MuNTkHits[5];   //[NMus]
   Int_t           MuNMatches[5];   //[NMus]
   Int_t           MuNChambers[5];   //[NMus]
   Double_t        MuCaloComp[5];   //[NMus]
   Double_t        MuSegmComp[5];   //[NMus]
   Int_t           MuTrackerMu[5];   //[NMus]
   Int_t           MuGMPT[5];   //[NMus]
   Int_t           MuID[5];   //[NMus]
   Int_t           MuMID[5];   //[NMus]
   Int_t           NEles;
   Int_t           ElGood[8];   //[NEles]
   Int_t           ElIsIso[8];   //[NEles]
   Double_t        ElPx[8];   //[NEles]
   Double_t        ElPy[8];   //[NEles]
   Double_t        ElPz[8];   //[NEles]
   Double_t        ElPt[8];   //[NEles]
   Double_t        ElPtE[8];   //[NEles]
   Double_t        ElE[8];   //[NEles]
   Double_t        ElEt[8];   //[NEles]
   Double_t        ElEta[8];   //[NEles]
   Double_t        ElPhi[8];   //[NEles]
   Double_t        ElD0BS[8];   //[NEles]
   Double_t        ElD0PV[8];   //[NEles]
   Double_t        ElD0E[8];   //[NEles]
   Double_t        ElDzBS[8];   //[NEles]
   Double_t        ElDzPV[8];   //[NEles]
   Double_t        ElDzE[8];   //[NEles]
   Double_t        ElIso[8];   //[NEles]
   Double_t        ElPtSum[8];   //[NEles]
   Double_t        ElEtSum[8];   //[NEles]
   Double_t        ElNChi2[8];   //[NEles]
   Int_t           ElCharge[8];   //[NEles]
   Int_t           ElIDTight[8];   //[NEles]
   Int_t           ElIDLoose[8];   //[NEles]
   Int_t           ElIDRobustTight[8];   //[NEles]
   Int_t           ElIDRobustLoose[8];   //[NEles]
   Int_t           ElInGap[8];   //[NEles]
   Int_t           ElEcalDriven[8];   //[NEles]
   Int_t           ElTrackerDriven[8];   //[NEles]
   Int_t           ElBasicClustersSize[8];   //[NEles]
   Double_t        Elfbrem[8];   //[NEles]
   Double_t        ElHcalOverEcal[8];   //[NEles]
   Double_t        ElE5x5[8];   //[NEles]
   Double_t        ElE2x5Max[8];   //[NEles]
   Double_t        ElSigmaIetaIeta[8];   //[NEles]
   Double_t        ElDeltaPhiSeedClusterAtCalo[8];   //[NEles]
   Double_t        ElDeltaEtaSeedClusterAtCalo[8];   //[NEles]
   Double_t        ElDeltaPhiSuperClusterAtVtx[8];   //[NEles]
   Double_t        ElDeltaEtaSuperClusterAtVtx[8];   //[NEles]
   Double_t        ElCaloEnergy[8];   //[NEles]
   Double_t        ElTrkMomAtVtx[8];   //[NEles]
   Double_t        ElESuperClusterOverP[8];   //[NEles]
   Int_t           ElIsInJet[8];   //[NEles]
   Double_t        ElSharedPx[8];   //[NEles]
   Double_t        ElSharedPy[8];   //[NEles]
   Double_t        ElSharedPz[8];   //[NEles]
   Double_t        ElSharedEnergy[8];   //[NEles]
   Int_t           ElDuplicateEl[8];   //[NEles]
   Int_t           ElID[8];   //[NEles]
   Int_t           ElMID[8];   //[NEles]
   Int_t           NPhotons;
   Int_t           PhotonGood[13];   //[NPhotons]
   Int_t           PhotonIsIso[13];   //[NPhotons]
   Double_t        PhotonPt[13];   //[NPhotons]
   Double_t        PhotonPx[13];   //[NPhotons]
   Double_t        PhotonPy[13];   //[NPhotons]
   Double_t        PhotonPz[13];   //[NPhotons]
   Double_t        PhotonEta[13];   //[NPhotons]
   Double_t        PhotonPhi[13];   //[NPhotons]
   Double_t        PhotonEnergy[13];   //[NPhotons]
   Double_t        PhotonCaloPositionX[13];   //[NPhotons]
   Double_t        PhotonCaloPositionY[13];   //[NPhotons]
   Double_t        PhotonCaloPositionZ[13];   //[NPhotons]
   Double_t        PhotonHoverE[13];   //[NPhotons]
   Double_t        PhotonH1overE[13];   //[NPhotons]
   Double_t        PhotonH2overE[13];   //[NPhotons]
   Double_t        PhotonSigmaIetaIeta[13];   //[NPhotons]
   Int_t           PhotonHasPixSeed[13];   //[NPhotons]
   Int_t           PhotonHasConvTrks[13];   //[NPhotons]
   Int_t           NJets;
   Int_t           JGood[15];   //[NJets]
   Double_t        JPx[15];   //[NJets]
   Double_t        JPy[15];   //[NJets]
   Double_t        JPz[15];   //[NJets]
   Double_t        JPt[15];   //[NJets]
   Double_t        JE[15];   //[NJets]
   Double_t        JEt[15];   //[NJets]
   Double_t        JEta[15];   //[NJets]
   Double_t        JPhi[15];   //[NJets]
   Double_t        JEMfrac[15];   //[NJets]
   Int_t           JNConstituents[15];   //[NJets]
   Double_t        JID_HPD[15];   //[NJets]
   Double_t        JID_RBX[15];   //[NJets]
   Double_t        JID_n90Hits[15];   //[NJets]
   Double_t        JID_SubDet1[15];   //[NJets]
   Double_t        JID_SubDet2[15];   //[NJets]
   Double_t        JID_SubDet3[15];   //[NJets]
   Double_t        JID_SubDet4[15];   //[NJets]
   Double_t        JID_resEMF[15];   //[NJets]
   Double_t        JID_HCALTow[15];   //[NJets]
   Double_t        JID_ECALTow[15];   //[NJets]
   Double_t        JEtaEMrms[15];   //[NJets]
   Double_t        JEtaHADrms[15];   //[NJets]
   Double_t        JPhiEMrms[15];   //[NJets]
   Double_t        JPhiHADrms[15];   //[NJets]
   Double_t        JbTagProb[15];   //[NJets]
   Double_t        JChfrac[15];   //[NJets]
   Int_t           JNAssoTracks[15];   //[NJets]
   Double_t        Jtrk1px[15];   //[NJets]
   Double_t        Jtrk1py[15];   //[NJets]
   Double_t        Jtrk1pz[15];   //[NJets]
   Double_t        Jtrk2px[15];   //[NJets]
   Double_t        Jtrk2py[15];   //[NJets]
   Double_t        Jtrk2pz[15];   //[NJets]
   Double_t        Jtrk3px[15];   //[NJets]
   Double_t        Jtrk3py[15];   //[NJets]
   Double_t        Jtrk3pz[15];   //[NJets]
   Double_t        JEcorr[15];   //[NJets]
   Double_t        JeMinDR[15];   //[NJets]
   Double_t        JVtxx[15];   //[NJets]
   Double_t        JVtxy[15];   //[NJets]
   Double_t        JVtxz[15];   //[NJets]
   Double_t        JVtxExx[15];   //[NJets]
   Double_t        JVtxEyx[15];   //[NJets]
   Double_t        JVtxEyy[15];   //[NJets]
   Double_t        JVtxEzy[15];   //[NJets]
   Double_t        JVtxEzz[15];   //[NJets]
   Double_t        JVtxEzx[15];   //[NJets]
   Double_t        JVtxNChi2[15];   //[NJets]
   Int_t           NTracks;
   Int_t           TrkGood[9];   //[NTracks]
   Double_t        TrkPt[9];   //[NTracks]
   Double_t        TrkEta[9];   //[NTracks]
   Double_t        TrkPhi[9];   //[NTracks]
   Double_t        TrkNChi2[9];   //[NTracks]
   Double_t        TrkNHits[9];   //[NTracks]
   Double_t        TrkPtSumx;
   Double_t        TrkPtSumy;
   Double_t        TrkPtSum;
   Double_t        TrkPtSumPhi;
   Double_t        SumEt;
   Double_t        ECALSumEt;
   Double_t        HCALSumEt;
   Double_t        ECALEsumx;
   Double_t        ECALEsumy;
   Double_t        ECALEsumz;
   Double_t        ECALMET;
   Double_t        ECALMETPhi;
   Double_t        ECALMETEta;
   Double_t        HCALEsumx;
   Double_t        HCALEsumy;
   Double_t        HCALEsumz;
   Double_t        HCALMET;
   Double_t        HCALMETPhi;
   Double_t        HCALMETeta;
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
   TBranch        *b_HLTResults;   //!
   TBranch        *b_L1PhysResults;   //!
   TBranch        *b_L1TechResults;   //!
   TBranch        *b_PrimVtxGood;   //!
   TBranch        *b_PrimVtxx;   //!
   TBranch        *b_PrimVtxy;   //!
   TBranch        *b_PrimVtxz;   //!
   TBranch        *b_PrimVtxxE;   //!
   TBranch        *b_PrimVtxyE;   //!
   TBranch        *b_PrimVtxzE;   //!
   TBranch        *b_PrimVtxNChi2;   //!
   TBranch        *b_PrimVtxNTracks;   //!
   TBranch        *b_PrimVtxPtSum;   //!
   TBranch        *b_Beamspotx;   //!
   TBranch        *b_Beamspoty;   //!
   TBranch        *b_Beamspotz;   //!
   TBranch        *b_NCaloTowers;   //!
   TBranch        *b_GoodEvent;   //!
   TBranch        *b_MaxMuExceed;   //!
   TBranch        *b_MaxElExceed;   //!
   TBranch        *b_MaxJetExceed;   //!
   TBranch        *b_MaxUncJetExceed;   //!
   TBranch        *b_MaxTrkExceed;   //!
   TBranch        *b_MaxPhotonsExceed;   //!
   TBranch        *b_NMus;   //!
   TBranch        *b_MuGood;   //!
   TBranch        *b_MuIsIso;   //!
   TBranch        *b_MuPx;   //!
   TBranch        *b_MuPy;   //!
   TBranch        *b_MuPz;   //!
   TBranch        *b_MuPt;   //!
   TBranch        *b_MuPtE;   //!
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
   TBranch        *b_ElGood;   //!
   TBranch        *b_ElIsIso;   //!
   TBranch        *b_ElPx;   //!
   TBranch        *b_ElPy;   //!
   TBranch        *b_ElPz;   //!
   TBranch        *b_ElPt;   //!
   TBranch        *b_ElPtE;   //!
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
   TBranch        *b_ElIDTight;   //!
   TBranch        *b_ElIDLoose;   //!
   TBranch        *b_ElIDRobustTight;   //!
   TBranch        *b_ElIDRobustLoose;   //!
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
   TBranch        *b_ElIsInJet;   //!
   TBranch        *b_ElSharedPx;   //!
   TBranch        *b_ElSharedPy;   //!
   TBranch        *b_ElSharedPz;   //!
   TBranch        *b_ElSharedEnergy;   //!
   TBranch        *b_ElDuplicateEl;   //!
   TBranch        *b_ElID;   //!
   TBranch        *b_ElMID;   //!
   TBranch        *b_NPhotons;   //!
   TBranch        *b_PhotonGood;   //!
   TBranch        *b_PhotonIsIso;   //!
   TBranch        *b_PhotonPt;   //!
   TBranch        *b_PhotonPx;   //!
   TBranch        *b_PhotonPy;   //!
   TBranch        *b_PhotonPz;   //!
   TBranch        *b_PhotonEta;   //!
   TBranch        *b_PhotonPhi;   //!
   TBranch        *b_PhotonEnergy;   //!
   TBranch        *b_PhotonCaloPositionX;   //!
   TBranch        *b_PhotonCaloPositionY;   //!
   TBranch        *b_PhotonCaloPositionZ;   //!
   TBranch        *b_PhotonHoverE;   //!
   TBranch        *b_PhotonH1overE;   //!
   TBranch        *b_PhotonH2overE;   //!
   TBranch        *b_PhotonSigmaIetaIeta;   //!
   TBranch        *b_PhotonHasPixSeed;   //!
   TBranch        *b_PhotonHasConvTrks;   //!
   TBranch        *b_NJets;   //!
   TBranch        *b_JGood;   //!
   TBranch        *b_JPx;   //!
   TBranch        *b_JPy;   //!
   TBranch        *b_JPz;   //!
   TBranch        *b_JPt;   //!
   TBranch        *b_JE;   //!
   TBranch        *b_JEt;   //!
   TBranch        *b_JEta;   //!
   TBranch        *b_JPhi;   //!
   TBranch        *b_JEMfrac;   //!
   TBranch        *b_JNConstituents;   //!
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
   TBranch        *b_JEtaEMrms;   //!
   TBranch        *b_JEtaHADrms;   //!
   TBranch        *b_JPhiEMrms;   //!
   TBranch        *b_JPhiHADrms;   //!
   TBranch        *b_JbTagProb;   //!
   TBranch        *b_JChfrac;   //!
   TBranch        *b_JNAssoTracks;   //!
   TBranch        *b_Jtrk1px;   //!
   TBranch        *b_Jtrk1py;   //!
   TBranch        *b_Jtrk1pz;   //!
   TBranch        *b_Jtrk2px;   //!
   TBranch        *b_Jtrk2py;   //!
   TBranch        *b_Jtrk2pz;   //!
   TBranch        *b_Jtrk3px;   //!
   TBranch        *b_Jtrk3py;   //!
   TBranch        *b_Jtrk3pz;   //!
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
   TBranch        *b_NTracks;   //!
   TBranch        *b_TrkGood;   //!
   TBranch        *b_TrkPt;   //!
   TBranch        *b_TrkEta;   //!
   TBranch        *b_TrkPhi;   //!
   TBranch        *b_TrkNChi2;   //!
   TBranch        *b_TrkNHits;   //!
   TBranch        *b_TrkPtSumx;   //!
   TBranch        *b_TrkPtSumy;   //!
   TBranch        *b_TrkPtSum;   //!
   TBranch        *b_TrkPtSumPhi;   //!
   TBranch        *b_SumEt;   //!
   TBranch        *b_ECALSumEt;   //!
   TBranch        *b_HCALSumEt;   //!
   TBranch        *b_ECALEsumx;   //!
   TBranch        *b_ECALEsumy;   //!
   TBranch        *b_ECALEsumz;   //!
   TBranch        *b_ECALMET;   //!
   TBranch        *b_ECALMETPhi;   //!
   TBranch        *b_ECALMETEta;   //!
   TBranch        *b_HCALEsumx;   //!
   TBranch        *b_HCALEsumy;   //!
   TBranch        *b_HCALEsumz;   //!
   TBranch        *b_HCALMET;   //!
   TBranch        *b_HCALMETPhi;   //!
   TBranch        *b_HCALMETEta;   //!
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
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("TTBar4Jets_NTP.root");
      if (!f) {
         f = new TFile("TTBar4Jets_NTP.root");
         f->cd("TTBar4Jets_NTP.root:/analyze");
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
   fChain->SetBranchAddress("HLTResults", HLTResults, &b_HLTResults);
   fChain->SetBranchAddress("L1PhysResults", L1PhysResults, &b_L1PhysResults);
   fChain->SetBranchAddress("L1TechResults", L1TechResults, &b_L1TechResults);
   fChain->SetBranchAddress("PrimVtxGood", &PrimVtxGood, &b_PrimVtxGood);
   fChain->SetBranchAddress("PrimVtxx", &PrimVtxx, &b_PrimVtxx);
   fChain->SetBranchAddress("PrimVtxy", &PrimVtxy, &b_PrimVtxy);
   fChain->SetBranchAddress("PrimVtxz", &PrimVtxz, &b_PrimVtxz);
   fChain->SetBranchAddress("PrimVtxxE", &PrimVtxxE, &b_PrimVtxxE);
   fChain->SetBranchAddress("PrimVtxyE", &PrimVtxyE, &b_PrimVtxyE);
   fChain->SetBranchAddress("PrimVtxzE", &PrimVtxzE, &b_PrimVtxzE);
   fChain->SetBranchAddress("PrimVtxNChi2", &PrimVtxNChi2, &b_PrimVtxNChi2);
   fChain->SetBranchAddress("PrimVtxNTracks", &PrimVtxNTracks, &b_PrimVtxNTracks);
   fChain->SetBranchAddress("PrimVtxPtSum", &PrimVtxPtSum, &b_PrimVtxPtSum);
   fChain->SetBranchAddress("Beamspotx", &Beamspotx, &b_Beamspotx);
   fChain->SetBranchAddress("Beamspoty", &Beamspoty, &b_Beamspoty);
   fChain->SetBranchAddress("Beamspotz", &Beamspotz, &b_Beamspotz);
   fChain->SetBranchAddress("NCaloTowers", &NCaloTowers, &b_NCaloTowers);
   fChain->SetBranchAddress("GoodEvent", &GoodEvent, &b_GoodEvent);
   fChain->SetBranchAddress("MaxMuExceed", &MaxMuExceed, &b_MaxMuExceed);
   fChain->SetBranchAddress("MaxElExceed", &MaxElExceed, &b_MaxElExceed);
   fChain->SetBranchAddress("MaxJetExceed", &MaxJetExceed, &b_MaxJetExceed);
   fChain->SetBranchAddress("MaxUncJetExceed", &MaxUncJetExceed, &b_MaxUncJetExceed);
   fChain->SetBranchAddress("MaxTrkExceed", &MaxTrkExceed, &b_MaxTrkExceed);
   fChain->SetBranchAddress("MaxPhotonsExceed", &MaxPhotonsExceed, &b_MaxPhotonsExceed);
   fChain->SetBranchAddress("NMus", &NMus, &b_NMus);
   fChain->SetBranchAddress("MuGood", MuGood, &b_MuGood);
   fChain->SetBranchAddress("MuIsIso", MuIsIso, &b_MuIsIso);
   fChain->SetBranchAddress("MuPx", MuPx, &b_MuPx);
   fChain->SetBranchAddress("MuPy", MuPy, &b_MuPy);
   fChain->SetBranchAddress("MuPz", MuPz, &b_MuPz);
   fChain->SetBranchAddress("MuPt", MuPt, &b_MuPt);
   fChain->SetBranchAddress("MuPtE", MuPtE, &b_MuPtE);
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
   fChain->SetBranchAddress("ElGood", ElGood, &b_ElGood);
   fChain->SetBranchAddress("ElIsIso", ElIsIso, &b_ElIsIso);
   fChain->SetBranchAddress("ElPx", ElPx, &b_ElPx);
   fChain->SetBranchAddress("ElPy", ElPy, &b_ElPy);
   fChain->SetBranchAddress("ElPz", ElPz, &b_ElPz);
   fChain->SetBranchAddress("ElPt", ElPt, &b_ElPt);
   fChain->SetBranchAddress("ElPtE", ElPtE, &b_ElPtE);
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
   fChain->SetBranchAddress("ElIDTight", ElIDTight, &b_ElIDTight);
   fChain->SetBranchAddress("ElIDLoose", ElIDLoose, &b_ElIDLoose);
   fChain->SetBranchAddress("ElIDRobustTight", ElIDRobustTight, &b_ElIDRobustTight);
   fChain->SetBranchAddress("ElIDRobustLoose", ElIDRobustLoose, &b_ElIDRobustLoose);
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
   fChain->SetBranchAddress("ElIsInJet", ElIsInJet, &b_ElIsInJet);
   fChain->SetBranchAddress("ElSharedPx", ElSharedPx, &b_ElSharedPx);
   fChain->SetBranchAddress("ElSharedPy", ElSharedPy, &b_ElSharedPy);
   fChain->SetBranchAddress("ElSharedPz", ElSharedPz, &b_ElSharedPz);
   fChain->SetBranchAddress("ElSharedEnergy", ElSharedEnergy, &b_ElSharedEnergy);
   fChain->SetBranchAddress("ElDuplicateEl", ElDuplicateEl, &b_ElDuplicateEl);
   fChain->SetBranchAddress("ElID", ElID, &b_ElID);
   fChain->SetBranchAddress("ElMID", ElMID, &b_ElMID);
   fChain->SetBranchAddress("NPhotons", &NPhotons, &b_NPhotons);
   fChain->SetBranchAddress("PhotonGood", PhotonGood, &b_PhotonGood);
   fChain->SetBranchAddress("PhotonIsIso", PhotonIsIso, &b_PhotonIsIso);
   fChain->SetBranchAddress("PhotonPt", PhotonPt, &b_PhotonPt);
   fChain->SetBranchAddress("PhotonPx", PhotonPx, &b_PhotonPx);
   fChain->SetBranchAddress("PhotonPy", PhotonPy, &b_PhotonPy);
   fChain->SetBranchAddress("PhotonPz", PhotonPz, &b_PhotonPz);
   fChain->SetBranchAddress("PhotonEta", PhotonEta, &b_PhotonEta);
   fChain->SetBranchAddress("PhotonPhi", PhotonPhi, &b_PhotonPhi);
   fChain->SetBranchAddress("PhotonEnergy", PhotonEnergy, &b_PhotonEnergy);
   fChain->SetBranchAddress("PhotonCaloPositionX", PhotonCaloPositionX, &b_PhotonCaloPositionX);
   fChain->SetBranchAddress("PhotonCaloPositionY", PhotonCaloPositionY, &b_PhotonCaloPositionY);
   fChain->SetBranchAddress("PhotonCaloPositionZ", PhotonCaloPositionZ, &b_PhotonCaloPositionZ);
   fChain->SetBranchAddress("PhotonHoverE", PhotonHoverE, &b_PhotonHoverE);
   fChain->SetBranchAddress("PhotonH1overE", PhotonH1overE, &b_PhotonH1overE);
   fChain->SetBranchAddress("PhotonH2overE", PhotonH2overE, &b_PhotonH2overE);
   fChain->SetBranchAddress("PhotonSigmaIetaIeta", PhotonSigmaIetaIeta, &b_PhotonSigmaIetaIeta);
   fChain->SetBranchAddress("PhotonHasPixSeed", PhotonHasPixSeed, &b_PhotonHasPixSeed);
   fChain->SetBranchAddress("PhotonHasConvTrks", PhotonHasConvTrks, &b_PhotonHasConvTrks);
   fChain->SetBranchAddress("NJets", &NJets, &b_NJets);
   fChain->SetBranchAddress("JGood", JGood, &b_JGood);
   fChain->SetBranchAddress("JPx", JPx, &b_JPx);
   fChain->SetBranchAddress("JPy", JPy, &b_JPy);
   fChain->SetBranchAddress("JPz", JPz, &b_JPz);
   fChain->SetBranchAddress("JPt", JPt, &b_JPt);
   fChain->SetBranchAddress("JE", JE, &b_JE);
   fChain->SetBranchAddress("JEt", JEt, &b_JEt);
   fChain->SetBranchAddress("JEta", JEta, &b_JEta);
   fChain->SetBranchAddress("JPhi", JPhi, &b_JPhi);
   fChain->SetBranchAddress("JEMfrac", JEMfrac, &b_JEMfrac);
   fChain->SetBranchAddress("JNConstituents", JNConstituents, &b_JNConstituents);
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
   fChain->SetBranchAddress("JEtaEMrms", JEtaEMrms, &b_JEtaEMrms);
   fChain->SetBranchAddress("JEtaHADrms", JEtaHADrms, &b_JEtaHADrms);
   fChain->SetBranchAddress("JPhiEMrms", JPhiEMrms, &b_JPhiEMrms);
   fChain->SetBranchAddress("JPhiHADrms", JPhiHADrms, &b_JPhiHADrms);
   fChain->SetBranchAddress("JbTagProb", JbTagProb, &b_JbTagProb);
   fChain->SetBranchAddress("JChfrac", JChfrac, &b_JChfrac);
   fChain->SetBranchAddress("JNAssoTracks", JNAssoTracks, &b_JNAssoTracks);
   fChain->SetBranchAddress("Jtrk1px", Jtrk1px, &b_Jtrk1px);
   fChain->SetBranchAddress("Jtrk1py", Jtrk1py, &b_Jtrk1py);
   fChain->SetBranchAddress("Jtrk1pz", Jtrk1pz, &b_Jtrk1pz);
   fChain->SetBranchAddress("Jtrk2px", Jtrk2px, &b_Jtrk2px);
   fChain->SetBranchAddress("Jtrk2py", Jtrk2py, &b_Jtrk2py);
   fChain->SetBranchAddress("Jtrk2pz", Jtrk2pz, &b_Jtrk2pz);
   fChain->SetBranchAddress("Jtrk3px", Jtrk3px, &b_Jtrk3px);
   fChain->SetBranchAddress("Jtrk3py", Jtrk3py, &b_Jtrk3py);
   fChain->SetBranchAddress("Jtrk3pz", Jtrk3pz, &b_Jtrk3pz);
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
   fChain->SetBranchAddress("NTracks", &NTracks, &b_NTracks);
   fChain->SetBranchAddress("TrkGood", TrkGood, &b_TrkGood);
   fChain->SetBranchAddress("TrkPt", TrkPt, &b_TrkPt);
   fChain->SetBranchAddress("TrkEta", TrkEta, &b_TrkEta);
   fChain->SetBranchAddress("TrkPhi", TrkPhi, &b_TrkPhi);
   fChain->SetBranchAddress("TrkNChi2", TrkNChi2, &b_TrkNChi2);
   fChain->SetBranchAddress("TrkNHits", TrkNHits, &b_TrkNHits);
   fChain->SetBranchAddress("TrkPtSumx", &TrkPtSumx, &b_TrkPtSumx);
   fChain->SetBranchAddress("TrkPtSumy", &TrkPtSumy, &b_TrkPtSumy);
   fChain->SetBranchAddress("TrkPtSum", &TrkPtSum, &b_TrkPtSum);
   fChain->SetBranchAddress("TrkPtSumPhi", &TrkPtSumPhi, &b_TrkPtSumPhi);
   fChain->SetBranchAddress("SumEt", &SumEt, &b_SumEt);
   fChain->SetBranchAddress("ECALSumEt", &ECALSumEt, &b_ECALSumEt);
   fChain->SetBranchAddress("HCALSumEt", &HCALSumEt, &b_HCALSumEt);
   fChain->SetBranchAddress("ECALEsumx", &ECALEsumx, &b_ECALEsumx);
   fChain->SetBranchAddress("ECALEsumy", &ECALEsumy, &b_ECALEsumy);
   fChain->SetBranchAddress("ECALEsumz", &ECALEsumz, &b_ECALEsumz);
   fChain->SetBranchAddress("ECALMET", &ECALMET, &b_ECALMET);
   fChain->SetBranchAddress("ECALMETPhi", &ECALMETPhi, &b_ECALMETPhi);
   fChain->SetBranchAddress("ECALMETEta", &ECALMETEta, &b_ECALMETEta);
   fChain->SetBranchAddress("HCALEsumx", &HCALEsumx, &b_HCALEsumx);
   fChain->SetBranchAddress("HCALEsumy", &HCALEsumy, &b_HCALEsumy);
   fChain->SetBranchAddress("HCALEsumz", &HCALEsumz, &b_HCALEsumz);
   fChain->SetBranchAddress("HCALMET", &HCALMET, &b_HCALMET);
   fChain->SetBranchAddress("HCALMETPhi", &HCALMETPhi, &b_HCALMETPhi);
   fChain->SetBranchAddress("HCALMETeta", &HCALMETeta, &b_HCALMETEta);
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
