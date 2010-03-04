//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Mar  4 18:43:49 2010 by ROOT version 5.22/00d
// from TTree Analysis/ETHZAnalysisTree
// found on file: NTupleProducer_34X_MC34x_RECO.root
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
   Int_t           NMusTot;
   Int_t           MuGood[20];   //[NMus]
   Int_t           MuIsIso[20];   //[NMus]
   Double_t        MuPx[20];   //[NMus]
   Double_t        MuPy[20];   //[NMus]
   Double_t        MuPz[20];   //[NMus]
   Double_t        MuPt[20];   //[NMus]
   Double_t        MuPtE[20];   //[NMus]
   Double_t        MuE[20];   //[NMus]
   Double_t        MuEt[20];   //[NMus]
   Double_t        MuEta[20];   //[NMus]
   Double_t        MuPhi[20];   //[NMus]
   Int_t           MuCharge[20];   //[NMus]
   Double_t        MuRelIso03[20];   //[NMus]
   Double_t        MuIso03SumPt[20];   //[NMus]
   Double_t        MuIso03EmEt[20];   //[NMus]
   Double_t        MuIso03HadEt[20];   //[NMus]
   Double_t        MuIso05SumPt[20];   //[NMus]
   Double_t        MuIso05EmEt[20];   //[NMus]
   Double_t        MuIso05HadEt[20];   //[NMus]
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
   Int_t           MuGenID[20];   //[NMus]
   Int_t           MuGenStatus[20];   //[NMus]
   Int_t           MuGenCharge[20];   //[NMus]
   Double_t        MuGenPt[20];   //[NMus]
   Double_t        MuGenEta[20];   //[NMus]
   Double_t        MuGenPhi[20];   //[NMus]
   Double_t        MuGenE[20];   //[NMus]
   Int_t           MuGenMID[20];   //[NMus]
   Int_t           MuGenMStatus[20];   //[NMus]
   Int_t           MuGenMCharge[20];   //[NMus]
   Double_t        MuGenMPt[20];   //[NMus]
   Double_t        MuGenMEta[20];   //[NMus]
   Double_t        MuGenMPhi[20];   //[NMus]
   Double_t        MuGenME[20];   //[NMus]
   Int_t           MuGenGMID[20];   //[NMus]
   Int_t           MuGenGMStatus[20];   //[NMus]
   Int_t           MuGenGMCharge[20];   //[NMus]
   Double_t        MuGenGMPt[20];   //[NMus]
   Double_t        MuGenGMEta[20];   //[NMus]
   Double_t        MuGenGMPhi[20];   //[NMus]
   Double_t        MuGenGME[20];   //[NMus]
   Int_t           NEles;
   Int_t           NElesTot;
   Int_t           ElGood[20];   //[NEles]
   Int_t           ElIsIso[20];   //[NEles]
   Int_t           ElChargeMisIDProb[20];   //[NEles]
   Double_t        ElPx[20];   //[NEles]
   Double_t        ElPy[20];   //[NEles]
   Double_t        ElPz[20];   //[NEles]
   Double_t        ElPt[20];   //[NEles]
   Double_t        ElPtE[20];   //[NEles]
   Double_t        ElE[20];   //[NEles]
   Double_t        ElEt[20];   //[NEles]
   Double_t        ElEta[20];   //[NEles]
   Double_t        ElTheta[20];   //[NEles]
   Double_t        ElPhi[20];   //[NEles]
   Double_t        ElD0BS[20];   //[NEles]
   Double_t        ElD0PV[20];   //[NEles]
   Double_t        ElD0E[20];   //[NEles]
   Double_t        ElDzBS[20];   //[NEles]
   Double_t        ElDzPV[20];   //[NEles]
   Double_t        ElDzE[20];   //[NEles]
   Double_t        ElRelIso04[20];   //[NEles]
   Double_t        ElDR03TkSumPt[20];   //[NEles]
   Double_t        ElDR04TkSumPt[20];   //[NEles]
   Double_t        ElDR03EcalRecHitSumEt[20];   //[NEles]
   Double_t        ElDR04EcalRecHitSumEt[20];   //[NEles]
   Double_t        ElDR03HcalTowerSumEt[20];   //[NEles]
   Double_t        ElDR04HcalTowerSumEt[20];   //[NEles]
   Double_t        ElNChi2[20];   //[NEles]
   Int_t           ElCharge[20];   //[NEles]
   Int_t           ElCInfoIsGsfCtfCons[20];   //[NEles]
   Int_t           ElCInfoIsGsfCtfScPixCons[20];   //[NEles]
   Int_t           ElCInfoIsGsfScPixCons[20];   //[NEles]
   Int_t           ElScPixCharge[20];   //[NEles]
   Double_t        ElClosestCtfTrackPt[20];   //[NEles]
   Double_t        ElClosestCtfTrackEta[20];   //[NEles]
   Double_t        ElClosestCtfTrackPhi[20];   //[NEles]
   Int_t           ElClosestCtfTrackCharge[20];   //[NEles]
   Int_t           ElIDTight[20];   //[NEles]
   Int_t           ElIDLoose[20];   //[NEles]
   Int_t           ElIDRobustTight[20];   //[NEles]
   Int_t           ElIDRobustLoose[20];   //[NEles]
   Int_t           ElInGap[20];   //[NEles]
   Int_t           ElEcalDriven[20];   //[NEles]
   Int_t           ElTrackerDriven[20];   //[NEles]
   Int_t           ElBasicClustersSize[20];   //[NEles]
   Double_t        Elfbrem[20];   //[NEles]
   Double_t        ElHcalOverEcal[20];   //[NEles]
   Double_t        ElE1x5[20];   //[NEles]
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
   Int_t           ElIsInJet[20];   //[NEles]
   Double_t        ElSharedPx[20];   //[NEles]
   Double_t        ElSharedPy[20];   //[NEles]
   Double_t        ElSharedPz[20];   //[NEles]
   Double_t        ElSharedEnergy[20];   //[NEles]
   Int_t           ElDuplicateEl[20];   //[NEles]
   Int_t           ElGenID[20];   //[NEles]
   Int_t           ElGenStatus[20];   //[NEles]
   Int_t           ElGenCharge[20];   //[NEles]
   Double_t        ElGenPt[20];   //[NEles]
   Double_t        ElGenEta[20];   //[NEles]
   Double_t        ElGenPhi[20];   //[NEles]
   Double_t        ElGenE[20];   //[NEles]
   Int_t           ElGenMID[20];   //[NEles]
   Int_t           ElGenMStatus[20];   //[NEles]
   Int_t           ElGenMCharge[20];   //[NEles]
   Double_t        ElGenMPt[20];   //[NEles]
   Double_t        ElGenMEta[20];   //[NEles]
   Double_t        ElGenMPhi[20];   //[NEles]
   Double_t        ElGenME[20];   //[NEles]
   Int_t           ElGenGMID[20];   //[NEles]
   Int_t           ElGenGMStatus[20];   //[NEles]
   Int_t           ElGenGMCharge[20];   //[NEles]
   Double_t        ElGenGMPt[20];   //[NEles]
   Double_t        ElGenGMEta[20];   //[NEles]
   Double_t        ElGenGMPhi[20];   //[NEles]
   Double_t        ElGenGME[20];   //[NEles]
   Int_t           NPhotons;
   Int_t           NPhotonsTot;
   Int_t           PhoGood[20];   //[NPhotons]
   Int_t           PhoIsIso[20];   //[NPhotons]
   Double_t        PhoPt[20];   //[NPhotons]
   Double_t        PhoPx[20];   //[NPhotons]
   Double_t        PhoPy[20];   //[NPhotons]
   Double_t        PhoPz[20];   //[NPhotons]
   Double_t        PhoEta[20];   //[NPhotons]
   Double_t        PhoPhi[20];   //[NPhotons]
   Double_t        PhoEnergy[20];   //[NPhotons]
   Double_t        PhoIso03Ecal[20];   //[NPhotons]
   Double_t        PhoIso03Hcal[20];   //[NPhotons]
   Double_t        PhoIso03TrkSolid[20];   //[NPhotons]
   Double_t        PhoIso03TrkHollow[20];   //[NPhotons]
   Double_t        PhoIso03[20];   //[NPhotons]
   Double_t        PhoCaloPositionX[20];   //[NPhotons]
   Double_t        PhoCaloPositionY[20];   //[NPhotons]
   Double_t        PhoCaloPositionZ[20];   //[NPhotons]
   Double_t        PhoHoverE[20];   //[NPhotons]
   Double_t        PhoH1overE[20];   //[NPhotons]
   Double_t        PhoH2overE[20];   //[NPhotons]
   Double_t        PhoSigmaIetaIeta[20];   //[NPhotons]
   Int_t           PhoHasPixSeed[20];   //[NPhotons]
   Int_t           PhoHasConvTrks[20];   //[NPhotons]
   Int_t           PhoIsInJet[20];   //[NPhotons]
   Double_t        PhoSharedPx[20];   //[NPhotons]
   Double_t        PhoSharedPy[20];   //[NPhotons]
   Double_t        PhoSharedPz[20];   //[NPhotons]
   Double_t        PhoSharedEnergy[20];   //[NPhotons]
   Int_t           NJets;
   Int_t           NJetsTot;
   Int_t           JGood[50];   //[NJets]
   Double_t        JPx[50];   //[NJets]
   Double_t        JPy[50];   //[NJets]
   Double_t        JPz[50];   //[NJets]
   Double_t        JPt[50];   //[NJets]
   Double_t        JE[50];   //[NJets]
   Double_t        JEt[50];   //[NJets]
   Double_t        JEta[50];   //[NJets]
   Double_t        JPhi[50];   //[NJets]
   Double_t        JEMfrac[50];   //[NJets]
   Int_t           JNConstituents[50];   //[NJets]
   Double_t        JID_HPD[50];   //[NJets]
   Double_t        JID_RBX[50];   //[NJets]
   Double_t        JID_n90Hits[50];   //[NJets]
   Double_t        JID_SubDet1[50];   //[NJets]
   Double_t        JID_SubDet2[50];   //[NJets]
   Double_t        JID_SubDet3[50];   //[NJets]
   Double_t        JID_SubDet4[50];   //[NJets]
   Double_t        JID_resEMF[50];   //[NJets]
   Double_t        JID_HCALTow[50];   //[NJets]
   Double_t        JID_ECALTow[50];   //[NJets]
   Double_t        JEtaEMrms[50];   //[NJets]
   Double_t        JEtaHADrms[50];   //[NJets]
   Double_t        JPhiEMrms[50];   //[NJets]
   Double_t        JPhiHADrms[50];   //[NJets]
   Double_t        JbTagProb[50];   //[NJets]
   Double_t        JChfrac[50];   //[NJets]
   Double_t        JMass[50];   //[NJets]
   Int_t           JNAssoTracks[50];   //[NJets]
   Double_t        Jtrk1px[50];   //[NJets]
   Double_t        Jtrk1py[50];   //[NJets]
   Double_t        Jtrk1pz[50];   //[NJets]
   Double_t        Jtrk2px[50];   //[NJets]
   Double_t        Jtrk2py[50];   //[NJets]
   Double_t        Jtrk2pz[50];   //[NJets]
   Double_t        Jtrk3px[50];   //[NJets]
   Double_t        Jtrk3py[50];   //[NJets]
   Double_t        Jtrk3pz[50];   //[NJets]
   Double_t        JEcorr[50];   //[NJets]
   Double_t        JeMinDR[50];   //[NJets]
   Double_t        JVtxx[50];   //[NJets]
   Double_t        JVtxy[50];   //[NJets]
   Double_t        JVtxz[50];   //[NJets]
   Double_t        JVtxExx[50];   //[NJets]
   Double_t        JVtxEyx[50];   //[NJets]
   Double_t        JVtxEyy[50];   //[NJets]
   Double_t        JVtxEzy[50];   //[NJets]
   Double_t        JVtxEzz[50];   //[NJets]
   Double_t        JVtxEzx[50];   //[NJets]
   Double_t        JVtxNChi2[50];   //[NJets]
   Double_t        JGenPt[50];   //[NJets]
   Double_t        JGenEta[50];   //[NJets]
   Double_t        JGenPhi[50];   //[NJets]
   Double_t        JGenE[50];   //[NJets]
   Double_t        JGenEmE[50];   //[NJets]
   Double_t        JGenHadE[50];   //[NJets]
   Double_t        JGenInvE[50];   //[NJets]
   Int_t           NTracks;
   Int_t           NTracksTot;
   Int_t           TrkGood[500];   //[NTracks]
   Double_t        TrkPt[500];   //[NTracks]
   Double_t        TrkEta[500];   //[NTracks]
   Double_t        TrkPhi[500];   //[NTracks]
   Double_t        TrkNChi2[500];   //[NTracks]
   Double_t        TrkNHits[500];   //[NTracks]
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
   Double_t        RawMETemEtFrac;
   Double_t        RawMETemEtInEB;
   Double_t        RawMETemEtInEE;
   Double_t        RawMETemEtInHF;
   Double_t        RawMEThadEtFrac;
   Double_t        RawMEThadEtInHB;
   Double_t        RawMEThadEtInHE;
   Double_t        RawMEThadEtInHF;
   Double_t        RawMETSignificance;
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
   TBranch        *b_NMusTot;   //!
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
   TBranch        *b_MuRelIso03;   //!
   TBranch        *b_MuIso03SumPt;   //!
   TBranch        *b_MuIso03EmEt;   //!
   TBranch        *b_MuIso03HadEt;   //!
   TBranch        *b_MuIso05SumPt;   //!
   TBranch        *b_MuIso05EmEt;   //!
   TBranch        *b_MuIso05HadEt;   //!
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
   TBranch        *b_MuGenID;   //!
   TBranch        *b_MuGenStatus;   //!
   TBranch        *b_MuGenCharge;   //!
   TBranch        *b_MuGenPt;   //!
   TBranch        *b_MuGenEta;   //!
   TBranch        *b_MuGenPhi;   //!
   TBranch        *b_MuGenE;   //!
   TBranch        *b_MuGenMID;   //!
   TBranch        *b_MuGenMStatus;   //!
   TBranch        *b_MuGenMCharge;   //!
   TBranch        *b_MuGenMPt;   //!
   TBranch        *b_MuGenMEta;   //!
   TBranch        *b_MuGenMPhi;   //!
   TBranch        *b_MuGenME;   //!
   TBranch        *b_MuGenGMID;   //!
   TBranch        *b_MuGenGMStatus;   //!
   TBranch        *b_MuGenGMCharge;   //!
   TBranch        *b_MuGenGMPt;   //!
   TBranch        *b_MuGenGMEta;   //!
   TBranch        *b_MuGenGMPhi;   //!
   TBranch        *b_MuGenGME;   //!
   TBranch        *b_NEles;   //!
   TBranch        *b_NElesTot;   //!
   TBranch        *b_ElGood;   //!
   TBranch        *b_ElIsIso;   //!
   TBranch        *b_ElChargeMisIDProb;   //!
   TBranch        *b_ElPx;   //!
   TBranch        *b_ElPy;   //!
   TBranch        *b_ElPz;   //!
   TBranch        *b_ElPt;   //!
   TBranch        *b_ElPtE;   //!
   TBranch        *b_ElE;   //!
   TBranch        *b_ElEt;   //!
   TBranch        *b_ElEta;   //!
   TBranch        *b_ElTheta;   //!
   TBranch        *b_ElPhi;   //!
   TBranch        *b_ElD0BS;   //!
   TBranch        *b_ElD0PV;   //!
   TBranch        *b_ElD0E;   //!
   TBranch        *b_ElDzBS;   //!
   TBranch        *b_ElDzPV;   //!
   TBranch        *b_ElDzE;   //!
   TBranch        *b_ElRelIso04;   //!
   TBranch        *b_ElDR03TkSumPt;   //!
   TBranch        *b_ElDR04TkSumPt;   //!
   TBranch        *b_ElDR03EcalRecHitSumEt;   //!
   TBranch        *b_ElDR04EcalRecHitSumEt;   //!
   TBranch        *b_ElDR03HcalTowerSumEt;   //!
   TBranch        *b_ElDR04HcalTowerSumEt;   //!
   TBranch        *b_ElNChi2;   //!
   TBranch        *b_ElCharge;   //!
   TBranch        *b_ElCInfoIsGsfCtfCons;   //!
   TBranch        *b_ElCInfoIsGsfCtfScPixCons;   //!
   TBranch        *b_ElCInfoIsGsfScPixCons;   //!
   TBranch        *b_ElScPixCharge;   //!
   TBranch        *b_ElClosestCtfTrackPt;   //!
   TBranch        *b_ElClosestCtfTrackEta;   //!
   TBranch        *b_ElClosestCtfTrackPhi;   //!
   TBranch        *b_ElClosestCtfTrackCharge;   //!
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
   TBranch        *b_ElE1x5;   //!
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
   TBranch        *b_ElGenID;   //!
   TBranch        *b_ElGenStatus;   //!
   TBranch        *b_ElGenCharge;   //!
   TBranch        *b_ElGenPt;   //!
   TBranch        *b_ElGenEta;   //!
   TBranch        *b_ElGenPhi;   //!
   TBranch        *b_ElGenE;   //!
   TBranch        *b_ElGenMID;   //!
   TBranch        *b_ElGenMStatus;   //!
   TBranch        *b_ElGenMCharge;   //!
   TBranch        *b_ElGenMPt;   //!
   TBranch        *b_ElGenMEta;   //!
   TBranch        *b_ElGenMPhi;   //!
   TBranch        *b_ElGenME;   //!
   TBranch        *b_ElGenGMID;   //!
   TBranch        *b_ElGenGMStatus;   //!
   TBranch        *b_ElGenGMCharge;   //!
   TBranch        *b_ElGenGMPt;   //!
   TBranch        *b_ElGenGMEta;   //!
   TBranch        *b_ElGenGMPhi;   //!
   TBranch        *b_ElGenGME;   //!
   TBranch        *b_NPhotons;   //!
   TBranch        *b_NPhotonsTot;   //!
   TBranch        *b_PhoGood;   //!
   TBranch        *b_PhoIsIso;   //!
   TBranch        *b_PhoPt;   //!
   TBranch        *b_PhoPx;   //!
   TBranch        *b_PhoPy;   //!
   TBranch        *b_PhoPz;   //!
   TBranch        *b_PhoEta;   //!
   TBranch        *b_PhoPhi;   //!
   TBranch        *b_PhoEnergy;   //!
   TBranch        *b_PhoIso03Ecal;   //!
   TBranch        *b_PhoIso03Hcal;   //!
   TBranch        *b_PhoIso03TrkSolid;   //!
   TBranch        *b_PhoIso03TrkHollow;   //!
   TBranch        *b_PhoIso03;   //!
   TBranch        *b_PhoCaloPositionX;   //!
   TBranch        *b_PhoCaloPositionY;   //!
   TBranch        *b_PhoCaloPositionZ;   //!
   TBranch        *b_PhoHoverE;   //!
   TBranch        *b_PhoH1overE;   //!
   TBranch        *b_PhoH2overE;   //!
   TBranch        *b_PhoSigmaIetaIeta;   //!
   TBranch        *b_PhoHasPixSeed;   //!
   TBranch        *b_PhoHasConvTrks;   //!
   TBranch        *b_PhoIsInJet;   //!
   TBranch        *b_PhoSharedPx;   //!
   TBranch        *b_PhoSharedPy;   //!
   TBranch        *b_PhoSharedPz;   //!
   TBranch        *b_PhoSharedEnergy;   //!
   TBranch        *b_NJets;   //!
   TBranch        *b_NJetsTot;   //!
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
   TBranch        *b_JMass;   //!
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
   TBranch        *b_JGenPt;   //!
   TBranch        *b_JGenEta;   //!
   TBranch        *b_JGenPhi;   //!
   TBranch        *b_JGenE;   //!
   TBranch        *b_JGenEmE;   //!
   TBranch        *b_JGenHadE;   //!
   TBranch        *b_JGenInvE;   //!
   TBranch        *b_NTracks;   //!
   TBranch        *b_NTracksTot;   //!
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
   TBranch        *b_RawMETemEtFrac;   //!
   TBranch        *b_RawMETemEtInEB;   //!
   TBranch        *b_RawMETemEtInEE;   //!
   TBranch        *b_RawMETemEtInHF;   //!
   TBranch        *b_RawMEThadEtFrac;   //!
   TBranch        *b_RawMEThadEtInHB;   //!
   TBranch        *b_RawMEThadEtInHE;   //!
   TBranch        *b_RawMEThadEtInHF;   //!
   TBranch        *b_RawMETSignificance;   //!
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
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("NTupleProducer_34X_MC34x_RECO.root");
      if (!f) {
         f = new TFile("NTupleProducer_34X_MC34x_RECO.root");
         f->cd("NTupleProducer_34X_MC34x_RECO.root:/analyze");
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
   fChain->SetBranchAddress("NMusTot", &NMusTot, &b_NMusTot);
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
   fChain->SetBranchAddress("MuRelIso03", MuRelIso03, &b_MuRelIso03);
   fChain->SetBranchAddress("MuIso03SumPt", MuIso03SumPt, &b_MuIso03SumPt);
   fChain->SetBranchAddress("MuIso03EmEt", MuIso03EmEt, &b_MuIso03EmEt);
   fChain->SetBranchAddress("MuIso03HadEt", MuIso03HadEt, &b_MuIso03HadEt);
   fChain->SetBranchAddress("MuIso05SumPt", MuIso05SumPt, &b_MuIso05SumPt);
   fChain->SetBranchAddress("MuIso05EmEt", MuIso05EmEt, &b_MuIso05EmEt);
   fChain->SetBranchAddress("MuIso05HadEt", MuIso05HadEt, &b_MuIso05HadEt);
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
   fChain->SetBranchAddress("MuGenID", MuGenID, &b_MuGenID);
   fChain->SetBranchAddress("MuGenStatus", MuGenStatus, &b_MuGenStatus);
   fChain->SetBranchAddress("MuGenCharge", MuGenCharge, &b_MuGenCharge);
   fChain->SetBranchAddress("MuGenPt", MuGenPt, &b_MuGenPt);
   fChain->SetBranchAddress("MuGenEta", MuGenEta, &b_MuGenEta);
   fChain->SetBranchAddress("MuGenPhi", MuGenPhi, &b_MuGenPhi);
   fChain->SetBranchAddress("MuGenE", MuGenE, &b_MuGenE);
   fChain->SetBranchAddress("MuGenMID", MuGenMID, &b_MuGenMID);
   fChain->SetBranchAddress("MuGenMStatus", MuGenMStatus, &b_MuGenMStatus);
   fChain->SetBranchAddress("MuGenMCharge", MuGenMCharge, &b_MuGenMCharge);
   fChain->SetBranchAddress("MuGenMPt", MuGenMPt, &b_MuGenMPt);
   fChain->SetBranchAddress("MuGenMEta", MuGenMEta, &b_MuGenMEta);
   fChain->SetBranchAddress("MuGenMPhi", MuGenMPhi, &b_MuGenMPhi);
   fChain->SetBranchAddress("MuGenME", MuGenME, &b_MuGenME);
   fChain->SetBranchAddress("MuGenGMID", MuGenGMID, &b_MuGenGMID);
   fChain->SetBranchAddress("MuGenGMStatus", MuGenGMStatus, &b_MuGenGMStatus);
   fChain->SetBranchAddress("MuGenGMCharge", MuGenGMCharge, &b_MuGenGMCharge);
   fChain->SetBranchAddress("MuGenGMPt", MuGenGMPt, &b_MuGenGMPt);
   fChain->SetBranchAddress("MuGenGMEta", MuGenGMEta, &b_MuGenGMEta);
   fChain->SetBranchAddress("MuGenGMPhi", MuGenGMPhi, &b_MuGenGMPhi);
   fChain->SetBranchAddress("MuGenGME", MuGenGME, &b_MuGenGME);
   fChain->SetBranchAddress("NEles", &NEles, &b_NEles);
   fChain->SetBranchAddress("NElesTot", &NElesTot, &b_NElesTot);
   fChain->SetBranchAddress("ElGood", ElGood, &b_ElGood);
   fChain->SetBranchAddress("ElIsIso", ElIsIso, &b_ElIsIso);
   fChain->SetBranchAddress("ElChargeMisIDProb", ElChargeMisIDProb, &b_ElChargeMisIDProb);
   fChain->SetBranchAddress("ElPx", ElPx, &b_ElPx);
   fChain->SetBranchAddress("ElPy", ElPy, &b_ElPy);
   fChain->SetBranchAddress("ElPz", ElPz, &b_ElPz);
   fChain->SetBranchAddress("ElPt", ElPt, &b_ElPt);
   fChain->SetBranchAddress("ElPtE", ElPtE, &b_ElPtE);
   fChain->SetBranchAddress("ElE", ElE, &b_ElE);
   fChain->SetBranchAddress("ElEt", ElEt, &b_ElEt);
   fChain->SetBranchAddress("ElEta", ElEta, &b_ElEta);
   fChain->SetBranchAddress("ElTheta", ElTheta, &b_ElTheta);
   fChain->SetBranchAddress("ElPhi", ElPhi, &b_ElPhi);
   fChain->SetBranchAddress("ElD0BS", ElD0BS, &b_ElD0BS);
   fChain->SetBranchAddress("ElD0PV", ElD0PV, &b_ElD0PV);
   fChain->SetBranchAddress("ElD0E", ElD0E, &b_ElD0E);
   fChain->SetBranchAddress("ElDzBS", ElDzBS, &b_ElDzBS);
   fChain->SetBranchAddress("ElDzPV", ElDzPV, &b_ElDzPV);
   fChain->SetBranchAddress("ElDzE", ElDzE, &b_ElDzE);
   fChain->SetBranchAddress("ElRelIso04", ElRelIso04, &b_ElRelIso04);
   fChain->SetBranchAddress("ElDR03TkSumPt", ElDR03TkSumPt, &b_ElDR03TkSumPt);
   fChain->SetBranchAddress("ElDR04TkSumPt", ElDR04TkSumPt, &b_ElDR04TkSumPt);
   fChain->SetBranchAddress("ElDR03EcalRecHitSumEt", ElDR03EcalRecHitSumEt, &b_ElDR03EcalRecHitSumEt);
   fChain->SetBranchAddress("ElDR04EcalRecHitSumEt", ElDR04EcalRecHitSumEt, &b_ElDR04EcalRecHitSumEt);
   fChain->SetBranchAddress("ElDR03HcalTowerSumEt", ElDR03HcalTowerSumEt, &b_ElDR03HcalTowerSumEt);
   fChain->SetBranchAddress("ElDR04HcalTowerSumEt", ElDR04HcalTowerSumEt, &b_ElDR04HcalTowerSumEt);
   fChain->SetBranchAddress("ElNChi2", ElNChi2, &b_ElNChi2);
   fChain->SetBranchAddress("ElCharge", ElCharge, &b_ElCharge);
   fChain->SetBranchAddress("ElCInfoIsGsfCtfCons", ElCInfoIsGsfCtfCons, &b_ElCInfoIsGsfCtfCons);
   fChain->SetBranchAddress("ElCInfoIsGsfCtfScPixCons", ElCInfoIsGsfCtfScPixCons, &b_ElCInfoIsGsfCtfScPixCons);
   fChain->SetBranchAddress("ElCInfoIsGsfScPixCons", ElCInfoIsGsfScPixCons, &b_ElCInfoIsGsfScPixCons);
   fChain->SetBranchAddress("ElScPixCharge", ElScPixCharge, &b_ElScPixCharge);
   fChain->SetBranchAddress("ElClosestCtfTrackPt", ElClosestCtfTrackPt, &b_ElClosestCtfTrackPt);
   fChain->SetBranchAddress("ElClosestCtfTrackEta", ElClosestCtfTrackEta, &b_ElClosestCtfTrackEta);
   fChain->SetBranchAddress("ElClosestCtfTrackPhi", ElClosestCtfTrackPhi, &b_ElClosestCtfTrackPhi);
   fChain->SetBranchAddress("ElClosestCtfTrackCharge", ElClosestCtfTrackCharge, &b_ElClosestCtfTrackCharge);
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
   fChain->SetBranchAddress("ElE1x5", ElE1x5, &b_ElE1x5);
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
   fChain->SetBranchAddress("ElGenID", ElGenID, &b_ElGenID);
   fChain->SetBranchAddress("ElGenStatus", ElGenStatus, &b_ElGenStatus);
   fChain->SetBranchAddress("ElGenCharge", ElGenCharge, &b_ElGenCharge);
   fChain->SetBranchAddress("ElGenPt", ElGenPt, &b_ElGenPt);
   fChain->SetBranchAddress("ElGenEta", ElGenEta, &b_ElGenEta);
   fChain->SetBranchAddress("ElGenPhi", ElGenPhi, &b_ElGenPhi);
   fChain->SetBranchAddress("ElGenE", ElGenE, &b_ElGenE);
   fChain->SetBranchAddress("ElGenMID", ElGenMID, &b_ElGenMID);
   fChain->SetBranchAddress("ElGenMStatus", ElGenMStatus, &b_ElGenMStatus);
   fChain->SetBranchAddress("ElGenMCharge", ElGenMCharge, &b_ElGenMCharge);
   fChain->SetBranchAddress("ElGenMPt", ElGenMPt, &b_ElGenMPt);
   fChain->SetBranchAddress("ElGenMEta", ElGenMEta, &b_ElGenMEta);
   fChain->SetBranchAddress("ElGenMPhi", ElGenMPhi, &b_ElGenMPhi);
   fChain->SetBranchAddress("ElGenME", ElGenME, &b_ElGenME);
   fChain->SetBranchAddress("ElGenGMID", ElGenGMID, &b_ElGenGMID);
   fChain->SetBranchAddress("ElGenGMStatus", ElGenGMStatus, &b_ElGenGMStatus);
   fChain->SetBranchAddress("ElGenGMCharge", ElGenGMCharge, &b_ElGenGMCharge);
   fChain->SetBranchAddress("ElGenGMPt", ElGenGMPt, &b_ElGenGMPt);
   fChain->SetBranchAddress("ElGenGMEta", ElGenGMEta, &b_ElGenGMEta);
   fChain->SetBranchAddress("ElGenGMPhi", ElGenGMPhi, &b_ElGenGMPhi);
   fChain->SetBranchAddress("ElGenGME", ElGenGME, &b_ElGenGME);
   fChain->SetBranchAddress("NPhotons", &NPhotons, &b_NPhotons);
   fChain->SetBranchAddress("NPhotonsTot", &NPhotonsTot, &b_NPhotonsTot);
   fChain->SetBranchAddress("PhoGood", PhoGood, &b_PhoGood);
   fChain->SetBranchAddress("PhoIsIso", PhoIsIso, &b_PhoIsIso);
   fChain->SetBranchAddress("PhoPt", PhoPt, &b_PhoPt);
   fChain->SetBranchAddress("PhoPx", PhoPx, &b_PhoPx);
   fChain->SetBranchAddress("PhoPy", PhoPy, &b_PhoPy);
   fChain->SetBranchAddress("PhoPz", PhoPz, &b_PhoPz);
   fChain->SetBranchAddress("PhoEta", PhoEta, &b_PhoEta);
   fChain->SetBranchAddress("PhoPhi", PhoPhi, &b_PhoPhi);
   fChain->SetBranchAddress("PhoEnergy", PhoEnergy, &b_PhoEnergy);
   fChain->SetBranchAddress("PhoIso03Ecal", PhoIso03Ecal, &b_PhoIso03Ecal);
   fChain->SetBranchAddress("PhoIso03Hcal", PhoIso03Hcal, &b_PhoIso03Hcal);
   fChain->SetBranchAddress("PhoIso03TrkSolid", PhoIso03TrkSolid, &b_PhoIso03TrkSolid);
   fChain->SetBranchAddress("PhoIso03TrkHollow", PhoIso03TrkHollow, &b_PhoIso03TrkHollow);
   fChain->SetBranchAddress("PhoIso03", PhoIso03, &b_PhoIso03);
   fChain->SetBranchAddress("PhoCaloPositionX", PhoCaloPositionX, &b_PhoCaloPositionX);
   fChain->SetBranchAddress("PhoCaloPositionY", PhoCaloPositionY, &b_PhoCaloPositionY);
   fChain->SetBranchAddress("PhoCaloPositionZ", PhoCaloPositionZ, &b_PhoCaloPositionZ);
   fChain->SetBranchAddress("PhoHoverE", PhoHoverE, &b_PhoHoverE);
   fChain->SetBranchAddress("PhoH1overE", PhoH1overE, &b_PhoH1overE);
   fChain->SetBranchAddress("PhoH2overE", PhoH2overE, &b_PhoH2overE);
   fChain->SetBranchAddress("PhoSigmaIetaIeta", PhoSigmaIetaIeta, &b_PhoSigmaIetaIeta);
   fChain->SetBranchAddress("PhoHasPixSeed", PhoHasPixSeed, &b_PhoHasPixSeed);
   fChain->SetBranchAddress("PhoHasConvTrks", PhoHasConvTrks, &b_PhoHasConvTrks);
   fChain->SetBranchAddress("PhoIsInJet", PhoIsInJet, &b_PhoIsInJet);
   fChain->SetBranchAddress("PhoSharedPx", PhoSharedPx, &b_PhoSharedPx);
   fChain->SetBranchAddress("PhoSharedPy", PhoSharedPy, &b_PhoSharedPy);
   fChain->SetBranchAddress("PhoSharedPz", PhoSharedPz, &b_PhoSharedPz);
   fChain->SetBranchAddress("PhoSharedEnergy", PhoSharedEnergy, &b_PhoSharedEnergy);
   fChain->SetBranchAddress("NJets", &NJets, &b_NJets);
   fChain->SetBranchAddress("NJetsTot", &NJetsTot, &b_NJetsTot);
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
   fChain->SetBranchAddress("JMass", JMass, &b_JMass);
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
   fChain->SetBranchAddress("JGenPt", JGenPt, &b_JGenPt);
   fChain->SetBranchAddress("JGenEta", JGenEta, &b_JGenEta);
   fChain->SetBranchAddress("JGenPhi", JGenPhi, &b_JGenPhi);
   fChain->SetBranchAddress("JGenE", JGenE, &b_JGenE);
   fChain->SetBranchAddress("JGenEmE", JGenEmE, &b_JGenEmE);
   fChain->SetBranchAddress("JGenHadE", JGenHadE, &b_JGenHadE);
   fChain->SetBranchAddress("JGenInvE", JGenInvE, &b_JGenInvE);
   fChain->SetBranchAddress("NTracks", &NTracks, &b_NTracks);
   fChain->SetBranchAddress("NTracksTot", &NTracksTot, &b_NTracksTot);
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
   fChain->SetBranchAddress("RawMETemEtFrac", &RawMETemEtFrac, &b_RawMETemEtFrac);
   fChain->SetBranchAddress("RawMETemEtInEB", &RawMETemEtInEB, &b_RawMETemEtInEB);
   fChain->SetBranchAddress("RawMETemEtInEE", &RawMETemEtInEE, &b_RawMETemEtInEE);
   fChain->SetBranchAddress("RawMETemEtInHF", &RawMETemEtInHF, &b_RawMETemEtInHF);
   fChain->SetBranchAddress("RawMEThadEtFrac", &RawMEThadEtFrac, &b_RawMEThadEtFrac);
   fChain->SetBranchAddress("RawMEThadEtInHB", &RawMEThadEtInHB, &b_RawMEThadEtInHB);
   fChain->SetBranchAddress("RawMEThadEtInHE", &RawMEThadEtInHE, &b_RawMEThadEtInHE);
   fChain->SetBranchAddress("RawMEThadEtInHF", &RawMEThadEtInHF, &b_RawMEThadEtInHF);
   fChain->SetBranchAddress("RawMETSignificance", &RawMETSignificance, &b_RawMETSignificance);
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
