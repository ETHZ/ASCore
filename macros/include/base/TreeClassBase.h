//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sat Apr  2 18:29:06 2011 by ROOT version 5.27/06b
// from TTree Analysis/ETHZAnalysisTree
// found on file: /scratch/fronga/NTupleProducer_38X_data_74_1_015.root
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
   Float_t         PtHat;
   Int_t           SigProcID;
   Float_t         PDFScalePDF;
   Int_t           PDFID1;
   Int_t           PDFID2;
   Float_t         PDFx1;
   Float_t         PDFx2;
   Float_t         PDFxPDF1;
   Float_t         PDFxPDF2;
   Float_t         ExtXSecLO;
   Float_t         IntXSec;
   Int_t           PUnumInteractions;
   Float_t         PUzPositions[1];   //[PUnumInteractions]
   Float_t         PUsumPtLowPt[1];   //[PUnumInteractions]
   Float_t         PUsumPtHighPt[1];   //[PUnumInteractions]
   Float_t         PUnTrksLowPt[1];   //[PUnumInteractions]
   Float_t         PUnTrksHighPt[1];   //[PUnumInteractions]
   Float_t         Weight;
   Int_t           HLTResults[200];
   Int_t           L1PhysResults[128];
   Int_t           L1TechResults[64];
   Int_t           HLTPrescale[200];
   Int_t           NHLTObjs;
   Int_t           HLTObjectID[7][10];
   Float_t         HLTObjectPt[7][10];
   Float_t         HLTObjectEta[7][10];
   Float_t         HLTObjectPhi[7][10];
   Int_t           PrimVtxGood;
   Float_t         PrimVtxx;
   Float_t         PrimVtxy;
   Float_t         PrimVtxz;
   Float_t         PrimVtxRho;
   Float_t         PrimVtxxE;
   Float_t         PrimVtxyE;
   Float_t         PrimVtxzE;
   Float_t         PrimVtxNChi2;
   Float_t         PrimVtxNdof;
   Int_t           PrimVtxIsFake;
   Float_t         PrimVtxPtSum;
   Float_t         Beamspotx;
   Float_t         Beamspoty;
   Float_t         Beamspotz;
   Int_t           NCaloTowers;
   Int_t           GoodEvent;
   Int_t           MaxMuExceed;
   Int_t           MaxElExceed;
   Int_t           MaxJetExceed;
   Int_t           MaxUncJetExceed;
   Int_t           MaxTrkExceed;
   Int_t           MaxPhotonsExceed;
   Int_t           MaxGenLepExceed;
   Int_t           MaxGenJetExceed;
   Int_t           MaxVerticesExceed;
   Int_t           HBHENoiseFlag;
   Int_t           NGenLeptons;
   Int_t           GenLeptonID[100];   //[NGenLeptons]
   Float_t         GenLeptonPt[100];   //[NGenLeptons]
   Float_t         GenLeptonEta[100];   //[NGenLeptons]
   Float_t         GenLeptonPhi[100];   //[NGenLeptons]
   Int_t           GenLeptonMID[100];   //[NGenLeptons]
   Int_t           GenLeptonMStatus[100];   //[NGenLeptons]
   Float_t         GenLeptonMPt[100];   //[NGenLeptons]
   Float_t         GenLeptonMEta[100];   //[NGenLeptons]
   Float_t         GenLeptonMPhi[100];   //[NGenLeptons]
   Int_t           GenLeptonGMID[100];   //[NGenLeptons]
   Int_t           GenLeptonGMStatus[100];   //[NGenLeptons]
   Float_t         GenLeptonGMPt[100];   //[NGenLeptons]
   Float_t         GenLeptonGMEta[100];   //[NGenLeptons]
   Float_t         GenLeptonGMPhi[100];   //[NGenLeptons]
   Int_t           NGenJets;
   Float_t         GenJetPt[1];   //[NGenJets]
   Float_t         GenJetEta[1];   //[NGenJets]
   Float_t         GenJetPhi[1];   //[NGenJets]
   Float_t         GenJetE[1];   //[NGenJets]
   Float_t         GenJetEmE[1];   //[NGenJets]
   Float_t         GenJetHadE[1];   //[NGenJets]
   Float_t         GenJetInvE[1];   //[NGenJets]
   Int_t           NVrtx;
   Float_t         VrtxX[25];   //[NVrtx]
   Float_t         VrtxY[25];   //[NVrtx]
   Float_t         VrtxZ[25];   //[NVrtx]
   Float_t         VrtxXE[25];   //[NVrtx]
   Float_t         VrtxYE[25];   //[NVrtx]
   Float_t         VrtxZE[25];   //[NVrtx]
   Float_t         VrtxNdof[25];   //[NVrtx]
   Float_t         VrtxChi2[25];   //[NVrtx]
   Float_t         VrtxNtrks[25];   //[NVrtx]
   Float_t         VrtxSumPt[25];   //[NVrtx]
   Int_t           VrtxIsFake[25];   //[NVrtx]
   Int_t           NMus;
   Int_t           NMusTot;
   Int_t           NGMus;
   Int_t           NTMus;
   Int_t           MuGood[30];   //[NMus]
   Int_t           MuIsIso[30];   //[NMus]
   Int_t           MuIsGlobalMuon[30];   //[NMus]
   Int_t           MuIsTrackerMuon[30];   //[NMus]
   Float_t         MuPx[30];   //[NMus]
   Float_t         MuPy[30];   //[NMus]
   Float_t         MuPz[30];   //[NMus]
   Float_t         MuPt[30];   //[NMus]
   Float_t         MuInnerTkPt[30];   //[NMus]
   Float_t         MuPtE[30];   //[NMus]
   Float_t         MuE[30];   //[NMus]
   Float_t         MuEt[30];   //[NMus]
   Float_t         MuEta[30];   //[NMus]
   Float_t         MuPhi[30];   //[NMus]
   Int_t           MuCharge[30];   //[NMus]
   Float_t         MuRelIso03[30];   //[NMus]
   Float_t         MuIso03SumPt[30];   //[NMus]
   Float_t         MuIso03EmEt[30];   //[NMus]
   Float_t         MuIso03HadEt[30];   //[NMus]
   Float_t         MuIso03EMVetoEt[30];   //[NMus]
   Float_t         MuIso03HadVetoEt[30];   //[NMus]
   Float_t         MuIso05SumPt[30];   //[NMus]
   Float_t         MuIso05EmEt[30];   //[NMus]
   Float_t         MuIso05HadEt[30];   //[NMus]
   Float_t         MuEem[30];   //[NMus]
   Float_t         MuEhad[30];   //[NMus]
   Float_t         MuD0BS[30];   //[NMus]
   Float_t         MuD0PV[30];   //[NMus]
   Float_t         MuD0E[30];   //[NMus]
   Float_t         MuDzBS[30];   //[NMus]
   Float_t         MuDzPV[30];   //[NMus]
   Float_t         MuDzE[30];   //[NMus]
   Float_t         MuNChi2[30];   //[NMus]
   Int_t           MuNGlHits[30];   //[NMus]
   Int_t           MuNMuHits[30];   //[NMus]
   Int_t           MuNTkHits[30];   //[NMus]
   Int_t           MuNPxHits[30];   //[NMus]
   Float_t         MuInnerTkNChi2[30];   //[NMus]
   Int_t           MuNMatches[30];   //[NMus]
   Int_t           MuNChambers[30];   //[NMus]
   Float_t         MuCaloComp[30];   //[NMus]
   Float_t         MuSegmComp[30];   //[NMus]
   Int_t           MuIsGMPT[30];   //[NMus]
   Int_t           MuIsGMTkChiComp[30];   //[NMus]
   Int_t           MuIsGMStaChiComp[30];   //[NMus]
   Int_t           MuIsGMTkKinkTight[30];   //[NMus]
   Int_t           MuIsAllStaMuons[30];   //[NMus]
   Int_t           MuIsAllTrkMuons[30];   //[NMus]
   Int_t           MuIsTrkMuonArbitrated[30];   //[NMus]
   Int_t           MuIsAllArbitrated[30];   //[NMus]
   Int_t           MuIsTMLSLoose[30];   //[NMus]
   Int_t           MuIsTMLSTight[30];   //[NMus]
   Int_t           MuIsTM2DCompLoose[30];   //[NMus]
   Int_t           MuIsTM2DCompTight[30];   //[NMus]
   Int_t           MuIsTMOneStationLoose[30];   //[NMus]
   Int_t           MuIsTMOneStationTight[30];   //[NMus]
   Int_t           MuIsTMLSOptLowPtLoose[30];   //[NMus]
   Int_t           MuIsTMLSAngLoose[30];   //[NMus]
   Int_t           MuIsTMLastStationAngTight[30];   //[NMus]
   Int_t           MuIsTMOneStationAngTight[30];   //[NMus]
   Int_t           MuIsTMOneStationAngLoose[30];   //[NMus]
   Int_t           MuGenID[30];   //[NMus]
   Int_t           MuGenStatus[30];   //[NMus]
   Float_t         MuGenPt[30];   //[NMus]
   Float_t         MuGenEta[30];   //[NMus]
   Float_t         MuGenPhi[30];   //[NMus]
   Float_t         MuGenE[30];   //[NMus]
   Int_t           MuGenMID[30];   //[NMus]
   Int_t           MuGenMStatus[30];   //[NMus]
   Float_t         MuGenMPt[30];   //[NMus]
   Float_t         MuGenMEta[30];   //[NMus]
   Float_t         MuGenMPhi[30];   //[NMus]
   Float_t         MuGenME[30];   //[NMus]
   Int_t           MuGenGMID[30];   //[NMus]
   Int_t           MuGenGMStatus[30];   //[NMus]
   Float_t         MuGenGMPt[30];   //[NMus]
   Float_t         MuGenGMEta[30];   //[NMus]
   Float_t         MuGenGMPhi[30];   //[NMus]
   Float_t         MuGenGME[30];   //[NMus]
   Int_t           NPfMus;
   Int_t           NPfMusTot;
   Float_t         PfMuPx[1];   //[NPfMus]
   Float_t         PfMuPy[1];   //[NPfMus]
   Float_t         PfMuPz[1];   //[NPfMus]
   Float_t         PfMuPt[1];   //[NPfMus]
   Float_t         PfMuPtE[1];   //[NPfMus]
   Float_t         PfMuE[1];   //[NPfMus]
   Float_t         PfMuEt[1];   //[NPfMus]
   Float_t         PfMuEta[1];   //[NPfMus]
   Float_t         PfMuPhi[1];   //[NPfMus]
   Int_t           PfMuCharge[1];   //[NPfMus]
   Float_t         PfMuParticleIso[1];   //[NPfMus]
   Float_t         PfMuChargedHadronIso[1];   //[NPfMus]
   Float_t         PfMuNeutralHadronIso[1];   //[NPfMus]
   Float_t         PfMuPhotonIso[1];   //[NPfMus]
   Int_t           NEles;
   Int_t           NElesTot;
   Int_t           ElGood[20];   //[NEles]
   Int_t           ElIsIso[20];   //[NEles]
   Int_t           ElChargeMisIDProb[20];   //[NEles]
   Float_t         ElPx[20];   //[NEles]
   Float_t         ElPy[20];   //[NEles]
   Float_t         ElPz[20];   //[NEles]
   Float_t         ElPt[20];   //[NEles]
   Float_t         ElPtE[20];   //[NEles]
   Float_t         ElE[20];   //[NEles]
   Float_t         ElEt[20];   //[NEles]
   Float_t         ElEta[20];   //[NEles]
   Float_t         ElTheta[20];   //[NEles]
   Float_t         ElSCEta[20];   //[NEles]
   Float_t         ElPhi[20];   //[NEles]
   Float_t         ElGsfTkPt[20];   //[NEles]
   Float_t         ElGsfTkEta[20];   //[NEles]
   Float_t         ElGsfTkPhi[20];   //[NEles]
   Float_t         ElTrkMomentumError[20];   //[NEles]
   Float_t         ElEcalEnergyError[20];   //[NEles]
   Float_t         ElEleMomentumError[20];   //[NEles]
   Int_t           ElNBrems[20];   //[NEles]
   Float_t         ElD0BS[20];   //[NEles]
   Float_t         ElD0PV[20];   //[NEles]
   Float_t         ElD0E[20];   //[NEles]
   Float_t         ElDzBS[20];   //[NEles]
   Float_t         ElDzPV[20];   //[NEles]
   Float_t         ElDzE[20];   //[NEles]
   Float_t         ElRelIso03[20];   //[NEles]
   Float_t         ElRelIso04[20];   //[NEles]
   Float_t         ElDR03TkSumPt[20];   //[NEles]
   Float_t         ElDR04TkSumPt[20];   //[NEles]
   Float_t         ElDR03EcalRecHitSumEt[20];   //[NEles]
   Float_t         ElDR04EcalRecHitSumEt[20];   //[NEles]
   Float_t         ElDR03HcalTowerSumEt[20];   //[NEles]
   Float_t         ElDR04HcalTowerSumEt[20];   //[NEles]
   Float_t         ElNChi2[20];   //[NEles]
   Int_t           ElCharge[20];   //[NEles]
   Int_t           ElCInfoIsGsfCtfCons[20];   //[NEles]
   Int_t           ElCInfoIsGsfCtfScPixCons[20];   //[NEles]
   Int_t           ElCInfoIsGsfScPixCons[20];   //[NEles]
   Int_t           ElScPixCharge[20];   //[NEles]
   Float_t         ElClosestCtfTrackPt[20];   //[NEles]
   Float_t         ElClosestCtfTrackEta[20];   //[NEles]
   Float_t         ElClosestCtfTrackPhi[20];   //[NEles]
   Int_t           ElClosestCtfTrackCharge[20];   //[NEles]
   Float_t         ElIDMva[20];   //[NEles]
   Int_t           ElIDTight[20];   //[NEles]
   Int_t           ElIDLoose[20];   //[NEles]
   Int_t           ElIDRobustTight[20];   //[NEles]
   Int_t           ElIDRobustLoose[20];   //[NEles]
   Int_t           ElIDsimpleWPrelIso[20];   //[NEles]
   Int_t           ElIDsimpleWP80relIso[20];   //[NEles]
   Int_t           ElIDsimpleWP85relIso[20];   //[NEles]
   Int_t           ElIDsimpleWP90relIso[20];   //[NEles]
   Int_t           ElIDsimpleWP95relIso[20];   //[NEles]
   Int_t           ElInGap[20];   //[NEles]
   Int_t           ElEcalDriven[20];   //[NEles]
   Int_t           ElTrackerDriven[20];   //[NEles]
   Int_t           ElBasicClustersSize[20];   //[NEles]
   Float_t         Elfbrem[20];   //[NEles]
   Float_t         ElHcalOverEcal[20];   //[NEles]
   Float_t         ElE1x5[20];   //[NEles]
   Float_t         ElE5x5[20];   //[NEles]
   Float_t         ElE2x5Max[20];   //[NEles]
   Float_t         ElSigmaIetaIeta[20];   //[NEles]
   Float_t         ElDeltaPhiSeedClusterAtCalo[20];   //[NEles]
   Float_t         ElDeltaEtaSeedClusterAtCalo[20];   //[NEles]
   Float_t         ElDeltaPhiSuperClusterAtVtx[20];   //[NEles]
   Float_t         ElDeltaEtaSuperClusterAtVtx[20];   //[NEles]
   Float_t         ElCaloEnergy[20];   //[NEles]
   Float_t         ElTrkMomAtVtx[20];   //[NEles]
   Float_t         ElESuperClusterOverP[20];   //[NEles]
   Int_t           ElNumberOfMissingInnerHits[20];   //[NEles]
   Float_t         ElConvPartnerTrkDist[20];   //[NEles]
   Float_t         ElConvPartnerTrkDCot[20];   //[NEles]
   Float_t         ElConvPartnerTrkPt[20];   //[NEles]
   Float_t         ElConvPartnerTrkEta[20];   //[NEles]
   Float_t         ElConvPartnerTrkPhi[20];   //[NEles]
   Float_t         ElConvPartnerTrkCharge[20];   //[NEles]
   Int_t           ElScSeedSeverity[20];   //[NEles]
   Float_t         ElE1OverE9[20];   //[NEles]
   Float_t         ElS4OverS1[20];   //[NEles]
   Int_t           ElGenID[20];   //[NEles]
   Int_t           ElGenStatus[20];   //[NEles]
   Float_t         ElGenPt[20];   //[NEles]
   Float_t         ElGenEta[20];   //[NEles]
   Float_t         ElGenPhi[20];   //[NEles]
   Float_t         ElGenE[20];   //[NEles]
   Int_t           ElGenMID[20];   //[NEles]
   Int_t           ElGenMStatus[20];   //[NEles]
   Float_t         ElGenMPt[20];   //[NEles]
   Float_t         ElGenMEta[20];   //[NEles]
   Float_t         ElGenMPhi[20];   //[NEles]
   Float_t         ElGenME[20];   //[NEles]
   Int_t           ElGenGMID[20];   //[NEles]
   Int_t           ElGenGMStatus[20];   //[NEles]
   Float_t         ElGenGMPt[20];   //[NEles]
   Float_t         ElGenGMEta[20];   //[NEles]
   Float_t         ElGenGMPhi[20];   //[NEles]
   Float_t         ElGenGME[20];   //[NEles]
   Int_t           NPfEls;
   Int_t           NPfElsTot;
   Float_t         PfElPx[1];   //[NPfEls]
   Float_t         PfElPy[1];   //[NPfEls]
   Float_t         PfElPz[1];   //[NPfEls]
   Float_t         PfElPt[1];   //[NPfEls]
   Float_t         PfElPtE[1];   //[NPfEls]
   Float_t         PfElE[1];   //[NPfEls]
   Float_t         PfElEt[1];   //[NPfEls]
   Float_t         PfElEta[1];   //[NPfEls]
   Float_t         PfElPhi[1];   //[NPfEls]
   Int_t           PfElCharge[1];   //[NPfEls]
   Float_t         PfElParticleIso[1];   //[NPfEls]
   Float_t         PfElChargedHadronIso[1];   //[NPfEls]
   Float_t         PfElNeutralHadronIso[1];   //[NPfEls]
   Float_t         PfElPhotonIso[1];   //[NPfEls]
   Int_t           NPfTaus;
   Int_t           NPfTausTot;
   Float_t         PfTauPx[1];   //[NPfTaus]
   Float_t         PfTauPy[1];   //[NPfTaus]
   Float_t         PfTauPz[1];   //[NPfTaus]
   Float_t         PfTauPt[1];   //[NPfTaus]
   Float_t         PfTauPtE[1];   //[NPfTaus]
   Float_t         PfTauE[1];   //[NPfTaus]
   Float_t         PfTauEt[1];   //[NPfTaus]
   Float_t         PfTauEta[1];   //[NPfTaus]
   Float_t         PfTauPhi[1];   //[NPfTaus]
   Int_t           PfTauCharge[1];   //[NPfTaus]
   Float_t         PfTauParticleIso[1];   //[NPfTaus]
   Float_t         PfTauChargedHadronIso[1];   //[NPfTaus]
   Float_t         PfTauNeutralHadronIso[1];   //[NPfTaus]
   Float_t         PfTauPhotonIso[1];   //[NPfTaus]
   Int_t           NPhotons;
   Int_t           NPhotonsTot;
   Int_t           PhoGood[50];   //[NPhotons]
   Int_t           PhoIsIso[50];   //[NPhotons]
   Float_t         PhoPt[50];   //[NPhotons]
   Float_t         PhoPx[50];   //[NPhotons]
   Float_t         PhoPy[50];   //[NPhotons]
   Float_t         PhoPz[50];   //[NPhotons]
   Float_t         PhoEta[50];   //[NPhotons]
   Float_t         PhoPhi[50];   //[NPhotons]
   Float_t         PhoEnergy[50];   //[NPhotons]
   Float_t         PhoIso03Ecal[50];   //[NPhotons]
   Float_t         PhoIso03Hcal[50];   //[NPhotons]
   Float_t         PhoIso03TrkSolid[50];   //[NPhotons]
   Float_t         PhoIso03TrkHollow[50];   //[NPhotons]
   Float_t         PhoIso03[50];   //[NPhotons]
   Float_t         PhoIso04Ecal[50];   //[NPhotons]
   Float_t         PhoIso04Hcal[50];   //[NPhotons]
   Float_t         PhoIso04TrkSolid[50];   //[NPhotons]
   Float_t         PhoIso04TrkHollow[50];   //[NPhotons]
   Float_t         PhoIso04[50];   //[NPhotons]
   Float_t         PhoR9[50];   //[NPhotons]
   Float_t         PhoCaloPositionX[50];   //[NPhotons]
   Float_t         PhoCaloPositionY[50];   //[NPhotons]
   Float_t         PhoCaloPositionZ[50];   //[NPhotons]
   Float_t         PhoHoverE[50];   //[NPhotons]
   Float_t         PhoH1overE[50];   //[NPhotons]
   Float_t         PhoH2overE[50];   //[NPhotons]
   Float_t         PhoSigmaIetaIeta[50];   //[NPhotons]
   Float_t         PhoSCRawEnergy[50];   //[NPhotons]
   Float_t         PhoSCEtaWidth[50];   //[NPhotons]
   Float_t         PhoSCSigmaPhiPhi[50];   //[NPhotons]
   Int_t           PhoHasPixSeed[50];   //[NPhotons]
   Int_t           PhoHasConvTrks[50];   //[NPhotons]
   Int_t           PhoScSeedSeverity[50];   //[NPhotons]
   Float_t         PhoE1OverE9[50];   //[NPhotons]
   Float_t         PhoS4OverS1[50];   //[NPhotons]
   Int_t           NJets;
   Int_t           NJetsTot;
   Int_t           JGood[100];   //[NJets]
   Float_t         JPx[100];   //[NJets]
   Float_t         JPy[100];   //[NJets]
   Float_t         JPz[100];   //[NJets]
   Float_t         JPt[100];   //[NJets]
   Float_t         JE[100];   //[NJets]
   Float_t         JEt[100];   //[NJets]
   Float_t         JEta[100];   //[NJets]
   Float_t         JPhi[100];   //[NJets]
   Float_t         JEcorr[100];   //[NJets]
   Float_t         JEtaRms[100];   //[NJets]
   Float_t         JPhiRms[100];   //[NJets]
   Int_t           JNConstituents[100];   //[NJets]
   Int_t           JNAssoTracks[100];   //[NJets]
   Int_t           JNNeutrals[100];   //[NJets]
   Float_t         JChargedEmFrac[100];   //[NJets]
   Float_t         JNeutralEmFrac[100];   //[NJets]
   Float_t         JChargedHadFrac[100];   //[NJets]
   Float_t         JNeutralHadFrac[100];   //[NJets]
   Float_t         JChargedMuEnergyFrac[100];   //[NJets]
   Float_t         JeMinDR[100];   //[NJets]
   Float_t         JbTagProbTkCntHighEff[100];   //[NJets]
   Float_t         JbTagProbTkCntHighPur[100];   //[NJets]
   Float_t         JbTagProbSimpSVHighEff[100];   //[NJets]
   Float_t         JbTagProbSimpSVHighPur[100];   //[NJets]
   Float_t         JMass[100];   //[NJets]
   Float_t         Jtrk1px[100];   //[NJets]
   Float_t         Jtrk1py[100];   //[NJets]
   Float_t         Jtrk1pz[100];   //[NJets]
   Float_t         Jtrk2px[100];   //[NJets]
   Float_t         Jtrk2py[100];   //[NJets]
   Float_t         Jtrk2pz[100];   //[NJets]
   Float_t         Jtrk3px[100];   //[NJets]
   Float_t         Jtrk3py[100];   //[NJets]
   Float_t         Jtrk3pz[100];   //[NJets]
   Float_t         JVtxx[100];   //[NJets]
   Float_t         JVtxy[100];   //[NJets]
   Float_t         JVtxz[100];   //[NJets]
   Float_t         JVtxExx[100];   //[NJets]
   Float_t         JVtxEyx[100];   //[NJets]
   Float_t         JVtxEyy[100];   //[NJets]
   Float_t         JVtxEzy[100];   //[NJets]
   Float_t         JVtxEzz[100];   //[NJets]
   Float_t         JVtxEzx[100];   //[NJets]
   Float_t         JVtxNChi2[100];   //[NJets]
   Float_t         JGenPt[100];   //[NJets]
   Float_t         JGenEta[100];   //[NJets]
   Float_t         JGenPhi[100];   //[NJets]
   Float_t         JGenE[100];   //[NJets]
   Float_t         JGenEmE[100];   //[NJets]
   Float_t         JGenHadE[100];   //[NJets]
   Float_t         JGenInvE[100];   //[NJets]
   Int_t           CANJets;
   Double_t        CAJPx[100];   //[CANJets]
   Double_t        CAJPy[100];   //[CANJets]
   Double_t        CAJPz[100];   //[CANJets]
   Double_t        CAJPt[100];   //[CANJets]
   Double_t        CAJE[100];   //[CANJets]
   Double_t        CAJEt[100];   //[CANJets]
   Double_t        CAJEta[100];   //[CANJets]
   Double_t        CAJPhi[100];   //[CANJets]
   Double_t        CAJScale[100];   //[CANJets]
   Double_t        CAJbTagProbTkCntHighEff[100];   //[CANJets]
   Double_t        CAJbTagProbTkCntHighPur[100];   //[CANJets]
   Double_t        CAJbTagProbSimpSVHighEff[100];   //[CANJets]
   Double_t        CAJbTagProbSimpSVHighPur[100];   //[CANJets]
   Double_t        CAJID_HPD[100];   //[CANJets]
   Double_t        CAJID_RBX[100];   //[CANJets]
   Double_t        CAJID_n90Hits[100];   //[CANJets]
   Double_t        CAJID_resEMF[100];   //[CANJets]
   Double_t        CAJEMfrac[100];   //[CANJets]
   Int_t           CAJNAssoTracks[100];   //[CANJets]
   Double_t        CAJChfrac[100];   //[CANJets]
   Int_t           CAJNConstituents[100];   //[CANJets]
   Int_t           PF2PATNJets;
   Double_t        PF2PATJPx[7];   //[PF2PATNJets]
   Double_t        PF2PATJPy[7];   //[PF2PATNJets]
   Double_t        PF2PATJPz[7];   //[PF2PATNJets]
   Double_t        PF2PATJPt[7];   //[PF2PATNJets]
   Double_t        PF2PATJE[7];   //[PF2PATNJets]
   Double_t        PF2PATJEt[7];   //[PF2PATNJets]
   Double_t        PF2PATJEta[7];   //[PF2PATNJets]
   Double_t        PF2PATJPhi[7];   //[PF2PATNJets]
   Double_t        PF2PATJScale[7];   //[PF2PATNJets]
   Double_t        PF2PATJbTagProbTkCntHighEff[7];   //[PF2PATNJets]
   Double_t        PF2PATJbTagProbTkCntHighPur[7];   //[PF2PATNJets]
   Double_t        PF2PATJbTagProbSimpSVHighEff[7];   //[PF2PATNJets]
   Double_t        PF2PATJbTagProbSimpSVHighPur[7];   //[PF2PATNJets]
   Int_t           PF2PATJChMult[7];   //[PF2PATNJets]
   Int_t           PF2PATJNeuMult[7];   //[PF2PATNJets]
   Double_t        PF2PATJChHadfrac[7];   //[PF2PATNJets]
   Double_t        PF2PATJNeuHadfrac[7];   //[PF2PATNJets]
   Double_t        PF2PATJChEmfrac[7];   //[PF2PATNJets]
   Double_t        PF2PATJNeuEmfrac[7];   //[PF2PATNJets]
   Int_t           PF2PATJNConstituents[7];   //[PF2PATNJets]
   Int_t           NTracks;
   Int_t           NTracksTot;
   Int_t           TrkGood[500];   //[NTracks]
   Float_t         TrkPt[500];   //[NTracks]
   Float_t         TrkEta[500];   //[NTracks]
   Float_t         TrkPhi[500];   //[NTracks]
   Float_t         TrkNChi2[500];   //[NTracks]
   Float_t         TrkNHits[500];   //[NTracks]
   Float_t         TrkPtSumx;
   Float_t         TrkPtSumy;
   Float_t         TrkPtSum;
   Float_t         TrkPtSumPhi;
   Float_t         SumEt;
   Float_t         ECALSumEt;
   Float_t         HCALSumEt;
   Float_t         ECALEsumx;
   Float_t         ECALEsumy;
   Float_t         ECALEsumz;
   Float_t         ECALMET;
   Float_t         ECALMETPhi;
   Float_t         ECALMETEta;
   Float_t         HCALEsumx;
   Float_t         HCALEsumy;
   Float_t         HCALEsumz;
   Float_t         HCALMET;
   Float_t         HCALMETPhi;
   Float_t         HCALMETeta;
   Float_t         RawMET;
   Float_t         RawMETpx;
   Float_t         RawMETpy;
   Float_t         RawMETphi;
   Float_t         RawMETemEtFrac;
   Float_t         RawMETemEtInEB;
   Float_t         RawMETemEtInEE;
   Float_t         RawMETemEtInHF;
   Float_t         RawMEThadEtFrac;
   Float_t         RawMEThadEtInHB;
   Float_t         RawMEThadEtInHE;
   Float_t         RawMEThadEtInHF;
   Float_t         RawMETSignificance;
   Float_t         GenMET;
   Float_t         GenMETpx;
   Float_t         GenMETpy;
   Float_t         GenMETphi;
   Float_t         TCMET;
   Float_t         TCMETpx;
   Float_t         TCMETpy;
   Float_t         TCMETphi;
   Float_t         TCMETSignificance;
   Float_t         MuJESCorrMET;
   Float_t         MuJESCorrMETpx;
   Float_t         MuJESCorrMETpy;
   Float_t         MuJESCorrMETphi;
   Float_t         PFMET;
   Float_t         PFMETpx;
   Float_t         PFMETpy;
   Float_t         PFMETphi;
   Float_t         PFMETSignificance;
   Float_t         METR12;
   Float_t         METR21;

   // List of branches
   TBranch        *b_Run;   //!
   TBranch        *b_Event;   //!
   TBranch        *b_LumiSection;   //!
   TBranch        *b_PtHat;   //!
   TBranch        *b_SigProcID;   //!
   TBranch        *b_PDFScalePDF;   //!
   TBranch        *b_PDFID1;   //!
   TBranch        *b_PDFID2;   //!
   TBranch        *b_PDFx1;   //!
   TBranch        *b_PDFx2;   //!
   TBranch        *b_PDFxPDF1;   //!
   TBranch        *b_PDFxPDF2;   //!
   TBranch        *b_ExtXSecLO;   //!
   TBranch        *b_IntXSec;   //!
   TBranch        *b_PUnumInteractions;   //!
   TBranch        *b_PUzPositions;   //!
   TBranch        *b_PUsumPtLowPt;   //!
   TBranch        *b_PUsumPtHighPt;   //!
   TBranch        *b_PUnTrksLowPt;   //!
   TBranch        *b_PUnTrksHighPt;   //!
   TBranch        *b_Weight;   //!
   TBranch        *b_HLTResults;   //!
   TBranch        *b_L1PhysResults;   //!
   TBranch        *b_L1TechResults;   //!
   TBranch        *b_HLTPrescale;   //!
   TBranch        *b_NHLTObjs;   //!
   TBranch        *b_HLTObjectID;   //!
   TBranch        *b_HLTObjectPt;   //!
   TBranch        *b_HLTObjectEta;   //!
   TBranch        *b_HLTObjectPhi;   //!
   TBranch        *b_PrimVtxGood;   //!
   TBranch        *b_PrimVtxx;   //!
   TBranch        *b_PrimVtxy;   //!
   TBranch        *b_PrimVtxz;   //!
   TBranch        *b_PrimVtxRho;   //!
   TBranch        *b_PrimVtxxE;   //!
   TBranch        *b_PrimVtxyE;   //!
   TBranch        *b_PrimVtxzE;   //!
   TBranch        *b_PrimVtxNChi2;   //!
   TBranch        *b_PrimVtxNdof;   //!
   TBranch        *b_PrimVtxIsFake;   //!
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
   TBranch        *b_MaxGenLepExceed;   //!
   TBranch        *b_MaxGenJetExceed;   //!
   TBranch        *b_MaxVerticesExceed;   //!
   TBranch        *b_HBHENoiseFlag;   //!
   TBranch        *b_NGenLeptons;   //!
   TBranch        *b_GenLeptonID;   //!
   TBranch        *b_GenLeptonPt;   //!
   TBranch        *b_GenLeptonEta;   //!
   TBranch        *b_GenLeptonPhi;   //!
   TBranch        *b_GenLeptonMID;   //!
   TBranch        *b_GenLeptonMStatus;   //!
   TBranch        *b_GenLeptonMPt;   //!
   TBranch        *b_GenLeptonMEta;   //!
   TBranch        *b_GenLeptonMPhi;   //!
   TBranch        *b_GenLeptonGMID;   //!
   TBranch        *b_GenLeptonGMStatus;   //!
   TBranch        *b_GenLeptonGMPt;   //!
   TBranch        *b_GenLeptonGMEta;   //!
   TBranch        *b_GenLeptonGMPhi;   //!
   TBranch        *b_NGenJets;   //!
   TBranch        *b_GenJetPt;   //!
   TBranch        *b_GenJetEta;   //!
   TBranch        *b_GenJetPhi;   //!
   TBranch        *b_GenJetE;   //!
   TBranch        *b_GenJetEmE;   //!
   TBranch        *b_GenJetHadE;   //!
   TBranch        *b_GenJetInvE;   //!
   TBranch        *b_NVrtx;   //!
   TBranch        *b_VrtxX;   //!
   TBranch        *b_VrtxY;   //!
   TBranch        *b_VrtxZ;   //!
   TBranch        *b_VrtxXE;   //!
   TBranch        *b_VrtxYE;   //!
   TBranch        *b_VrtxZE;   //!
   TBranch        *b_VrtxNdof;   //!
   TBranch        *b_VrtxChi2;   //!
   TBranch        *b_VrtxNtrks;   //!
   TBranch        *b_VrtxSumPt;   //!
   TBranch        *b_VrtxIsFake;   //!
   TBranch        *b_NMus;   //!
   TBranch        *b_NMusTot;   //!
   TBranch        *b_NGMus;   //!
   TBranch        *b_NTMus;   //!
   TBranch        *b_MuGood;   //!
   TBranch        *b_MuIsIso;   //!
   TBranch        *b_MuIsGlobalMuon;   //!
   TBranch        *b_MuIsTrackerMuon;   //!
   TBranch        *b_MuPx;   //!
   TBranch        *b_MuPy;   //!
   TBranch        *b_MuPz;   //!
   TBranch        *b_MuPt;   //!
   TBranch        *b_MuInnerTkPt;   //!
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
   TBranch        *b_MuIso03EMVetoEt;   //!
   TBranch        *b_MuIso03HadVetoEt;   //!
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
   TBranch        *b_MuNPxHits;   //!
   TBranch        *b_MuInnerTkNChi2;   //!
   TBranch        *b_MuNMatches;   //!
   TBranch        *b_MuNChambers;   //!
   TBranch        *b_MuCaloComp;   //!
   TBranch        *b_MuSegmComp;   //!
   TBranch        *b_MuIsGMPT;   //!
   TBranch        *b_MuIsGMTkChiComp;   //!
   TBranch        *b_MuIsGMStaChiComp;   //!
   TBranch        *b_MuIsGMTkKinkTight;   //!
   TBranch        *b_MuIsAllStaMuons;   //!
   TBranch        *b_MuIsAllTrkMuons;   //!
   TBranch        *b_MuIsTrkMuonArbitrated;   //!
   TBranch        *b_MuIsAllArbitrated;   //!
   TBranch        *b_MuIsTMLSLoose;   //!
   TBranch        *b_MuIsTMLSTight;   //!
   TBranch        *b_MuIsTM2DCompLoose;   //!
   TBranch        *b_MuIsTM2DCompTight;   //!
   TBranch        *b_MuIsTMOneStationLoose;   //!
   TBranch        *b_MuIsTMOneStationTight;   //!
   TBranch        *b_MuIsTMLSOptLowPtLoose;   //!
   TBranch        *b_MuIsTMLSAngLoose;   //!
   TBranch        *b_MuIsTMLastStationAngTight;   //!
   TBranch        *b_MuIsTMOneStationAngTight;   //!
   TBranch        *b_MuIsTMOneStationAngLoose;   //!
   TBranch        *b_MuGenID;   //!
   TBranch        *b_MuGenStatus;   //!
   TBranch        *b_MuGenPt;   //!
   TBranch        *b_MuGenEta;   //!
   TBranch        *b_MuGenPhi;   //!
   TBranch        *b_MuGenE;   //!
   TBranch        *b_MuGenMID;   //!
   TBranch        *b_MuGenMStatus;   //!
   TBranch        *b_MuGenMPt;   //!
   TBranch        *b_MuGenMEta;   //!
   TBranch        *b_MuGenMPhi;   //!
   TBranch        *b_MuGenME;   //!
   TBranch        *b_MuGenGMID;   //!
   TBranch        *b_MuGenGMStatus;   //!
   TBranch        *b_MuGenGMPt;   //!
   TBranch        *b_MuGenGMEta;   //!
   TBranch        *b_MuGenGMPhi;   //!
   TBranch        *b_MuGenGME;   //!
   TBranch        *b_NPfMus;   //!
   TBranch        *b_NPfMusTot;   //!
   TBranch        *b_PfMuPx;   //!
   TBranch        *b_PfMuPy;   //!
   TBranch        *b_PfMuPz;   //!
   TBranch        *b_PfMuPt;   //!
   TBranch        *b_PfMuPtE;   //!
   TBranch        *b_PfMuE;   //!
   TBranch        *b_PfMuEt;   //!
   TBranch        *b_PfMuEta;   //!
   TBranch        *b_PfMuPhi;   //!
   TBranch        *b_PfMuCharge;   //!
   TBranch        *b_PfMuParticleIso;   //!
   TBranch        *b_PfMuChargedHadronIso;   //!
   TBranch        *b_PfMuNeutralHadronIso;   //!
   TBranch        *b_PfMuPhotonIso;   //!
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
   TBranch        *b_ElSCEta;   //!
   TBranch        *b_ElPhi;   //!
   TBranch        *b_ElGsfTkPt;   //!
   TBranch        *b_ElGsfTkEta;   //!
   TBranch        *b_ElGsfTkPhi;   //!
   TBranch        *b_ElTrkMomentumError;   //!
   TBranch        *b_ElEcalEnergyError;   //!
   TBranch        *b_ElEleMomentumError;   //!
   TBranch        *b_ElNBrems;   //!
   TBranch        *b_ElD0BS;   //!
   TBranch        *b_ElD0PV;   //!
   TBranch        *b_ElD0E;   //!
   TBranch        *b_ElDzBS;   //!
   TBranch        *b_ElDzPV;   //!
   TBranch        *b_ElDzE;   //!
   TBranch        *b_ElRelIso03;   //!
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
   TBranch        *b_ElIDMva;   //!
   TBranch        *b_ElIDTight;   //!
   TBranch        *b_ElIDLoose;   //!
   TBranch        *b_ElIDRobustTight;   //!
   TBranch        *b_ElIDRobustLoose;   //!
   TBranch        *b_ElIDsimpleWPrelIso;   //!
   TBranch        *b_ElIDsimpleWP80relIso;   //!
   TBranch        *b_ElIDsimpleWP85relIso;   //!
   TBranch        *b_ElIDsimpleWP90relIso;   //!
   TBranch        *b_ElIDsimpleWP95relIso;   //!
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
   TBranch        *b_ElNumberOfMissingInnerHits;   //!
   TBranch        *b_ElConvPartnerTrkDist;   //!
   TBranch        *b_ElConvPartnerTrkDCot;   //!
   TBranch        *b_ElConvPartnerTrkPt;   //!
   TBranch        *b_ElConvPartnerTrkEta;   //!
   TBranch        *b_ElConvPartnerTrkPhi;   //!
   TBranch        *b_ElConvPartnerTrkCharge;   //!
   TBranch        *b_ElScSeedSeverity;   //!
   TBranch        *b_ElE1OverE9;   //!
   TBranch        *b_ElS4OverS1;   //!
   TBranch        *b_ElGenID;   //!
   TBranch        *b_ElGenStatus;   //!
   TBranch        *b_ElGenPt;   //!
   TBranch        *b_ElGenEta;   //!
   TBranch        *b_ElGenPhi;   //!
   TBranch        *b_ElGenE;   //!
   TBranch        *b_ElGenMID;   //!
   TBranch        *b_ElGenMStatus;   //!
   TBranch        *b_ElGenMPt;   //!
   TBranch        *b_ElGenMEta;   //!
   TBranch        *b_ElGenMPhi;   //!
   TBranch        *b_ElGenME;   //!
   TBranch        *b_ElGenGMID;   //!
   TBranch        *b_ElGenGMStatus;   //!
   TBranch        *b_ElGenGMPt;   //!
   TBranch        *b_ElGenGMEta;   //!
   TBranch        *b_ElGenGMPhi;   //!
   TBranch        *b_ElGenGME;   //!
   TBranch        *b_NPfEls;   //!
   TBranch        *b_NPfElsTot;   //!
   TBranch        *b_PfElPx;   //!
   TBranch        *b_PfElPy;   //!
   TBranch        *b_PfElPz;   //!
   TBranch        *b_PfElPt;   //!
   TBranch        *b_PfElPtE;   //!
   TBranch        *b_PfElE;   //!
   TBranch        *b_PfElEt;   //!
   TBranch        *b_PfElEta;   //!
   TBranch        *b_PfElPhi;   //!
   TBranch        *b_PfElCharge;   //!
   TBranch        *b_PfElParticleIso;   //!
   TBranch        *b_PfElChargedHadronIso;   //!
   TBranch        *b_PfElNeutralHadronIso;   //!
   TBranch        *b_PfElPhotonIso;   //!
   TBranch        *b_NPfTaus;   //!
   TBranch        *b_NPfTausTot;   //!
   TBranch        *b_PfTauPx;   //!
   TBranch        *b_PfTauPy;   //!
   TBranch        *b_PfTauPz;   //!
   TBranch        *b_PfTauPt;   //!
   TBranch        *b_PfTauPtE;   //!
   TBranch        *b_PfTauE;   //!
   TBranch        *b_PfTauEt;   //!
   TBranch        *b_PfTauEta;   //!
   TBranch        *b_PfTauPhi;   //!
   TBranch        *b_PfTauCharge;   //!
   TBranch        *b_PfTauParticleIso;   //!
   TBranch        *b_PfTauChargedHadronIso;   //!
   TBranch        *b_PfTauNeutralHadronIso;   //!
   TBranch        *b_PfTauPhotonIso;   //!
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
   TBranch        *b_PhoIso04Ecal;   //!
   TBranch        *b_PhoIso04Hcal;   //!
   TBranch        *b_PhoIso04TrkSolid;   //!
   TBranch        *b_PhoIso04TrkHollow;   //!
   TBranch        *b_PhoIso04;   //!
   TBranch        *b_PhoR9;   //!
   TBranch        *b_PhoCaloPositionX;   //!
   TBranch        *b_PhoCaloPositionY;   //!
   TBranch        *b_PhoCaloPositionZ;   //!
   TBranch        *b_PhoHoverE;   //!
   TBranch        *b_PhoH1overE;   //!
   TBranch        *b_PhoH2overE;   //!
   TBranch        *b_PhoSigmaIetaIeta;   //!
   TBranch        *b_PhoSCRawEnergy;   //!
   TBranch        *b_PhoSCEtaWidth;   //!
   TBranch        *b_PhoSCSigmaPhiPhi;   //!
   TBranch        *b_PhoHasPixSeed;   //!
   TBranch        *b_PhoHasConvTrks;   //!
   TBranch        *b_PhoScSeedSeverity;   //!
   TBranch        *b_PhoE1OverE9;   //!
   TBranch        *b_PhoS4OverS1;   //!
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
   TBranch        *b_JEcorr;   //!
   TBranch        *b_JEtaRms;   //!
   TBranch        *b_JPhiRms;   //!
   TBranch        *b_JNConstituents;   //!
   TBranch        *b_JNAssoTracks;   //!
   TBranch        *b_JNNeutrals;   //!
   TBranch        *b_JChargedEmFrac;   //!
   TBranch        *b_JNeutralEmFrac;   //!
   TBranch        *b_JChargedHadFrac;   //!
   TBranch        *b_JNeutralHadFrac;   //!
   TBranch        *b_JChargedMuEnergyFrac;   //!
   TBranch        *b_JeMinDR;   //!
   TBranch        *b_JbTagProbTkCntHighEff;   //!
   TBranch        *b_JbTagProbTkCntHighPur;   //!
   TBranch        *b_JbTagProbSimpSVHighEff;   //!
   TBranch        *b_JbTagProbSimpSVHighPur;   //!
   TBranch        *b_JMass;   //!
   TBranch        *b_Jtrk1px;   //!
   TBranch        *b_Jtrk1py;   //!
   TBranch        *b_Jtrk1pz;   //!
   TBranch        *b_Jtrk2px;   //!
   TBranch        *b_Jtrk2py;   //!
   TBranch        *b_Jtrk2pz;   //!
   TBranch        *b_Jtrk3px;   //!
   TBranch        *b_Jtrk3py;   //!
   TBranch        *b_Jtrk3pz;   //!
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
   TBranch        *b_CANJets;   //!
   TBranch        *b_CAJPx;   //!
   TBranch        *b_CAJPy;   //!
   TBranch        *b_CAJPz;   //!
   TBranch        *b_CAJPt;   //!
   TBranch        *b_CAJE;   //!
   TBranch        *b_CAJEt;   //!
   TBranch        *b_CAJEta;   //!
   TBranch        *b_CAJPhi;   //!
   TBranch        *b_CAJScale;   //!
   TBranch        *b_CAJbTagProbTkCntHighEff;   //!
   TBranch        *b_CAJbTagProbTkCntHighPur;   //!
   TBranch        *b_CAJbTagProbSimpSVHighEff;   //!
   TBranch        *b_CAJbTagProbSimpSVHighPur;   //!
   TBranch        *b_CAJID_HPD;   //!
   TBranch        *b_CAJID_RBX;   //!
   TBranch        *b_CAJID_n90Hits;   //!
   TBranch        *b_CAJID_resEMF;   //!
   TBranch        *b_CAJEMfrac;   //!
   TBranch        *b_CAJNAssoTracks;   //!
   TBranch        *b_CAJChfrac;   //!
   TBranch        *b_CAJNConstituents;   //!
   TBranch        *b_PF2PATNJets;   //!
   TBranch        *b_PF2PATJPx;   //!
   TBranch        *b_PF2PATJPy;   //!
   TBranch        *b_PF2PATJPz;   //!
   TBranch        *b_PF2PATJPt;   //!
   TBranch        *b_PF2PATJE;   //!
   TBranch        *b_PF2PATJEt;   //!
   TBranch        *b_PF2PATJEta;   //!
   TBranch        *b_PF2PATJPhi;   //!
   TBranch        *b_PF2PATJScale;   //!
   TBranch        *b_PF2PATJbTagProbTkCntHighEff;   //!
   TBranch        *b_PF2PATJbTagProbTkCntHighPur;   //!
   TBranch        *b_PF2PATJbTagProbSimpSVHighEff;   //!
   TBranch        *b_PF2PATJbTagProbSimpSVHighPur;   //!
   TBranch        *b_PF2PATJChMult;   //!
   TBranch        *b_PF2PATJNeuMult;   //!
   TBranch        *b_PF2PATJChHadfrac;   //!
   TBranch        *b_PF2PATJNeuHadfrac;   //!
   TBranch        *b_PF2PATJChEmfrac;   //!
   TBranch        *b_PF2PATJNeuEmfrac;   //!
   TBranch        *b_PF2PATJNConstituents;   //!
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
   TBranch        *b_GenMET;   //!
   TBranch        *b_GenMETpx;   //!
   TBranch        *b_GenMETpy;   //!
   TBranch        *b_GenMETphi;   //!
   TBranch        *b_TCMET;   //!
   TBranch        *b_TCMETpx;   //!
   TBranch        *b_TCMETpy;   //!
   TBranch        *b_TCMETphi;   //!
   TBranch        *b_TCMETSignificance;   //!
   TBranch        *b_MuJESCorrMET;   //!
   TBranch        *b_MuJESCorrMETpx;   //!
   TBranch        *b_MuJESCorrMETpy;   //!
   TBranch        *b_MuJESCorrMETphi;   //!
   TBranch        *b_PFMET;   //!
   TBranch        *b_PFMETpx;   //!
   TBranch        *b_PFMETpy;   //!
   TBranch        *b_PFMETphi;   //!
   TBranch        *b_PFMETSignificance;   //!
   TBranch        *b_METR12;   //!
   TBranch        *b_METR21;   //!

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
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/scratch/fronga/NTupleProducer_38X_data_74_1_015.root");
      if (!f) {
         f = new TFile("/scratch/fronga/NTupleProducer_38X_data_74_1_015.root");
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
   fChain->SetBranchAddress("PtHat", &PtHat, &b_PtHat);
   fChain->SetBranchAddress("SigProcID", &SigProcID, &b_SigProcID);
   fChain->SetBranchAddress("PDFScalePDF", &PDFScalePDF, &b_PDFScalePDF);
   fChain->SetBranchAddress("PDFID1", &PDFID1, &b_PDFID1);
   fChain->SetBranchAddress("PDFID2", &PDFID2, &b_PDFID2);
   fChain->SetBranchAddress("PDFx1", &PDFx1, &b_PDFx1);
   fChain->SetBranchAddress("PDFx2", &PDFx2, &b_PDFx2);
   fChain->SetBranchAddress("PDFxPDF1", &PDFxPDF1, &b_PDFxPDF1);
   fChain->SetBranchAddress("PDFxPDF2", &PDFxPDF2, &b_PDFxPDF2);
   fChain->SetBranchAddress("ExtXSecLO", &ExtXSecLO, &b_ExtXSecLO);
   fChain->SetBranchAddress("IntXSec", &IntXSec, &b_IntXSec);
   fChain->SetBranchAddress("PUnumInteractions", &PUnumInteractions, &b_PUnumInteractions);
   fChain->SetBranchAddress("PUzPositions", &PUzPositions, &b_PUzPositions);
   fChain->SetBranchAddress("PUsumPtLowPt", &PUsumPtLowPt, &b_PUsumPtLowPt);
   fChain->SetBranchAddress("PUsumPtHighPt", &PUsumPtHighPt, &b_PUsumPtHighPt);
   fChain->SetBranchAddress("PUnTrksLowPt", &PUnTrksLowPt, &b_PUnTrksLowPt);
   fChain->SetBranchAddress("PUnTrksHighPt", &PUnTrksHighPt, &b_PUnTrksHighPt);
   fChain->SetBranchAddress("Weight", &Weight, &b_Weight);
   fChain->SetBranchAddress("HLTResults", HLTResults, &b_HLTResults);
   fChain->SetBranchAddress("L1PhysResults", L1PhysResults, &b_L1PhysResults);
   fChain->SetBranchAddress("L1TechResults", L1TechResults, &b_L1TechResults);
   fChain->SetBranchAddress("HLTPrescale", HLTPrescale, &b_HLTPrescale);
   fChain->SetBranchAddress("NHLTObjs", &NHLTObjs, &b_NHLTObjs);
   fChain->SetBranchAddress("HLTObjectID", HLTObjectID, &b_HLTObjectID);
   fChain->SetBranchAddress("HLTObjectPt", HLTObjectPt, &b_HLTObjectPt);
   fChain->SetBranchAddress("HLTObjectEta", HLTObjectEta, &b_HLTObjectEta);
   fChain->SetBranchAddress("HLTObjectPhi", HLTObjectPhi, &b_HLTObjectPhi);
   fChain->SetBranchAddress("PrimVtxGood", &PrimVtxGood, &b_PrimVtxGood);
   fChain->SetBranchAddress("PrimVtxx", &PrimVtxx, &b_PrimVtxx);
   fChain->SetBranchAddress("PrimVtxy", &PrimVtxy, &b_PrimVtxy);
   fChain->SetBranchAddress("PrimVtxz", &PrimVtxz, &b_PrimVtxz);
   fChain->SetBranchAddress("PrimVtxRho", &PrimVtxRho, &b_PrimVtxRho);
   fChain->SetBranchAddress("PrimVtxxE", &PrimVtxxE, &b_PrimVtxxE);
   fChain->SetBranchAddress("PrimVtxyE", &PrimVtxyE, &b_PrimVtxyE);
   fChain->SetBranchAddress("PrimVtxzE", &PrimVtxzE, &b_PrimVtxzE);
   fChain->SetBranchAddress("PrimVtxNChi2", &PrimVtxNChi2, &b_PrimVtxNChi2);
   fChain->SetBranchAddress("PrimVtxNdof", &PrimVtxNdof, &b_PrimVtxNdof);
   fChain->SetBranchAddress("PrimVtxIsFake", &PrimVtxIsFake, &b_PrimVtxIsFake);
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
   fChain->SetBranchAddress("MaxGenLepExceed", &MaxGenLepExceed, &b_MaxGenLepExceed);
   fChain->SetBranchAddress("MaxGenJetExceed", &MaxGenJetExceed, &b_MaxGenJetExceed);
   fChain->SetBranchAddress("MaxVerticesExceed", &MaxVerticesExceed, &b_MaxVerticesExceed);
   fChain->SetBranchAddress("HBHENoiseFlag", &HBHENoiseFlag, &b_HBHENoiseFlag);
   fChain->SetBranchAddress("NGenLeptons", &NGenLeptons, &b_NGenLeptons);
   fChain->SetBranchAddress("GenLeptonID", &GenLeptonID, &b_GenLeptonID);
   fChain->SetBranchAddress("GenLeptonPt", &GenLeptonPt, &b_GenLeptonPt);
   fChain->SetBranchAddress("GenLeptonEta", &GenLeptonEta, &b_GenLeptonEta);
   fChain->SetBranchAddress("GenLeptonPhi", &GenLeptonPhi, &b_GenLeptonPhi);
   fChain->SetBranchAddress("GenLeptonMID", &GenLeptonMID, &b_GenLeptonMID);
   fChain->SetBranchAddress("GenLeptonMStatus", &GenLeptonMStatus, &b_GenLeptonMStatus);
   fChain->SetBranchAddress("GenLeptonMPt", &GenLeptonMPt, &b_GenLeptonMPt);
   fChain->SetBranchAddress("GenLeptonMEta", &GenLeptonMEta, &b_GenLeptonMEta);
   fChain->SetBranchAddress("GenLeptonMPhi", &GenLeptonMPhi, &b_GenLeptonMPhi);
   fChain->SetBranchAddress("GenLeptonGMID", &GenLeptonGMID, &b_GenLeptonGMID);
   fChain->SetBranchAddress("GenLeptonGMStatus", &GenLeptonGMStatus, &b_GenLeptonGMStatus);
   fChain->SetBranchAddress("GenLeptonGMPt", &GenLeptonGMPt, &b_GenLeptonGMPt);
   fChain->SetBranchAddress("GenLeptonGMEta", &GenLeptonGMEta, &b_GenLeptonGMEta);
   fChain->SetBranchAddress("GenLeptonGMPhi", &GenLeptonGMPhi, &b_GenLeptonGMPhi);
   fChain->SetBranchAddress("NGenJets", &NGenJets, &b_NGenJets);
   fChain->SetBranchAddress("GenJetPt", &GenJetPt, &b_GenJetPt);
   fChain->SetBranchAddress("GenJetEta", &GenJetEta, &b_GenJetEta);
   fChain->SetBranchAddress("GenJetPhi", &GenJetPhi, &b_GenJetPhi);
   fChain->SetBranchAddress("GenJetE", &GenJetE, &b_GenJetE);
   fChain->SetBranchAddress("GenJetEmE", &GenJetEmE, &b_GenJetEmE);
   fChain->SetBranchAddress("GenJetHadE", &GenJetHadE, &b_GenJetHadE);
   fChain->SetBranchAddress("GenJetInvE", &GenJetInvE, &b_GenJetInvE);
   fChain->SetBranchAddress("NVrtx", &NVrtx, &b_NVrtx);
   fChain->SetBranchAddress("VrtxX", VrtxX, &b_VrtxX);
   fChain->SetBranchAddress("VrtxY", VrtxY, &b_VrtxY);
   fChain->SetBranchAddress("VrtxZ", VrtxZ, &b_VrtxZ);
   fChain->SetBranchAddress("VrtxXE", VrtxXE, &b_VrtxXE);
   fChain->SetBranchAddress("VrtxYE", VrtxYE, &b_VrtxYE);
   fChain->SetBranchAddress("VrtxZE", VrtxZE, &b_VrtxZE);
   fChain->SetBranchAddress("VrtxNdof", VrtxNdof, &b_VrtxNdof);
   fChain->SetBranchAddress("VrtxChi2", VrtxChi2, &b_VrtxChi2);
   fChain->SetBranchAddress("VrtxNtrks", VrtxNtrks, &b_VrtxNtrks);
   fChain->SetBranchAddress("VrtxSumPt", VrtxSumPt, &b_VrtxSumPt);
   fChain->SetBranchAddress("VrtxIsFake", VrtxIsFake, &b_VrtxIsFake);
   fChain->SetBranchAddress("NMus", &NMus, &b_NMus);
   fChain->SetBranchAddress("NMusTot", &NMusTot, &b_NMusTot);
   fChain->SetBranchAddress("NGMus", &NGMus, &b_NGMus);
   fChain->SetBranchAddress("NTMus", &NTMus, &b_NTMus);
   fChain->SetBranchAddress("MuGood", MuGood, &b_MuGood);
   fChain->SetBranchAddress("MuIsIso", MuIsIso, &b_MuIsIso);
   fChain->SetBranchAddress("MuIsGlobalMuon", MuIsGlobalMuon, &b_MuIsGlobalMuon);
   fChain->SetBranchAddress("MuIsTrackerMuon", MuIsTrackerMuon, &b_MuIsTrackerMuon);
   fChain->SetBranchAddress("MuPx", MuPx, &b_MuPx);
   fChain->SetBranchAddress("MuPy", MuPy, &b_MuPy);
   fChain->SetBranchAddress("MuPz", MuPz, &b_MuPz);
   fChain->SetBranchAddress("MuPt", MuPt, &b_MuPt);
   fChain->SetBranchAddress("MuInnerTkPt", MuInnerTkPt, &b_MuInnerTkPt);
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
   fChain->SetBranchAddress("MuIso03EMVetoEt", MuIso03EMVetoEt, &b_MuIso03EMVetoEt);
   fChain->SetBranchAddress("MuIso03HadVetoEt", MuIso03HadVetoEt, &b_MuIso03HadVetoEt);
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
   fChain->SetBranchAddress("MuNPxHits", MuNPxHits, &b_MuNPxHits);
   fChain->SetBranchAddress("MuInnerTkNChi2", MuInnerTkNChi2, &b_MuInnerTkNChi2);
   fChain->SetBranchAddress("MuNMatches", MuNMatches, &b_MuNMatches);
   fChain->SetBranchAddress("MuNChambers", MuNChambers, &b_MuNChambers);
   fChain->SetBranchAddress("MuCaloComp", MuCaloComp, &b_MuCaloComp);
   fChain->SetBranchAddress("MuSegmComp", MuSegmComp, &b_MuSegmComp);
   fChain->SetBranchAddress("MuIsGMPT", MuIsGMPT, &b_MuIsGMPT);
   fChain->SetBranchAddress("MuIsGMTkChiComp", MuIsGMTkChiComp, &b_MuIsGMTkChiComp);
   fChain->SetBranchAddress("MuIsGMStaChiComp", MuIsGMStaChiComp, &b_MuIsGMStaChiComp);
   fChain->SetBranchAddress("MuIsGMTkKinkTight", MuIsGMTkKinkTight, &b_MuIsGMTkKinkTight);
   fChain->SetBranchAddress("MuIsAllStaMuons", MuIsAllStaMuons, &b_MuIsAllStaMuons);
   fChain->SetBranchAddress("MuIsAllTrkMuons", MuIsAllTrkMuons, &b_MuIsAllTrkMuons);
   fChain->SetBranchAddress("MuIsTrkMuonArbitrated", MuIsTrkMuonArbitrated, &b_MuIsTrkMuonArbitrated);
   fChain->SetBranchAddress("MuIsAllArbitrated", MuIsAllArbitrated, &b_MuIsAllArbitrated);
   fChain->SetBranchAddress("MuIsTMLSLoose", MuIsTMLSLoose, &b_MuIsTMLSLoose);
   fChain->SetBranchAddress("MuIsTMLSTight", MuIsTMLSTight, &b_MuIsTMLSTight);
   fChain->SetBranchAddress("MuIsTM2DCompLoose", MuIsTM2DCompLoose, &b_MuIsTM2DCompLoose);
   fChain->SetBranchAddress("MuIsTM2DCompTight", MuIsTM2DCompTight, &b_MuIsTM2DCompTight);
   fChain->SetBranchAddress("MuIsTMOneStationLoose", MuIsTMOneStationLoose, &b_MuIsTMOneStationLoose);
   fChain->SetBranchAddress("MuIsTMOneStationTight", MuIsTMOneStationTight, &b_MuIsTMOneStationTight);
   fChain->SetBranchAddress("MuIsTMLSOptLowPtLoose", MuIsTMLSOptLowPtLoose, &b_MuIsTMLSOptLowPtLoose);
   fChain->SetBranchAddress("MuIsTMLSAngLoose", MuIsTMLSAngLoose, &b_MuIsTMLSAngLoose);
   fChain->SetBranchAddress("MuIsTMLastStationAngTight", MuIsTMLastStationAngTight, &b_MuIsTMLastStationAngTight);
   fChain->SetBranchAddress("MuIsTMOneStationAngTight", MuIsTMOneStationAngTight, &b_MuIsTMOneStationAngTight);
   fChain->SetBranchAddress("MuIsTMOneStationAngLoose", MuIsTMOneStationAngLoose, &b_MuIsTMOneStationAngLoose);
   fChain->SetBranchAddress("MuGenID", MuGenID, &b_MuGenID);
   fChain->SetBranchAddress("MuGenStatus", MuGenStatus, &b_MuGenStatus);
   fChain->SetBranchAddress("MuGenPt", MuGenPt, &b_MuGenPt);
   fChain->SetBranchAddress("MuGenEta", MuGenEta, &b_MuGenEta);
   fChain->SetBranchAddress("MuGenPhi", MuGenPhi, &b_MuGenPhi);
   fChain->SetBranchAddress("MuGenE", MuGenE, &b_MuGenE);
   fChain->SetBranchAddress("MuGenMID", MuGenMID, &b_MuGenMID);
   fChain->SetBranchAddress("MuGenMStatus", MuGenMStatus, &b_MuGenMStatus);
   fChain->SetBranchAddress("MuGenMPt", MuGenMPt, &b_MuGenMPt);
   fChain->SetBranchAddress("MuGenMEta", MuGenMEta, &b_MuGenMEta);
   fChain->SetBranchAddress("MuGenMPhi", MuGenMPhi, &b_MuGenMPhi);
   fChain->SetBranchAddress("MuGenME", MuGenME, &b_MuGenME);
   fChain->SetBranchAddress("MuGenGMID", MuGenGMID, &b_MuGenGMID);
   fChain->SetBranchAddress("MuGenGMStatus", MuGenGMStatus, &b_MuGenGMStatus);
   fChain->SetBranchAddress("MuGenGMPt", MuGenGMPt, &b_MuGenGMPt);
   fChain->SetBranchAddress("MuGenGMEta", MuGenGMEta, &b_MuGenGMEta);
   fChain->SetBranchAddress("MuGenGMPhi", MuGenGMPhi, &b_MuGenGMPhi);
   fChain->SetBranchAddress("MuGenGME", MuGenGME, &b_MuGenGME);
   fChain->SetBranchAddress("NPfMus", &NPfMus, &b_NPfMus);
   fChain->SetBranchAddress("NPfMusTot", &NPfMusTot, &b_NPfMusTot);
   fChain->SetBranchAddress("PfMuPx", &PfMuPx, &b_PfMuPx);
   fChain->SetBranchAddress("PfMuPy", &PfMuPy, &b_PfMuPy);
   fChain->SetBranchAddress("PfMuPz", &PfMuPz, &b_PfMuPz);
   fChain->SetBranchAddress("PfMuPt", &PfMuPt, &b_PfMuPt);
   fChain->SetBranchAddress("PfMuPtE", &PfMuPtE, &b_PfMuPtE);
   fChain->SetBranchAddress("PfMuE", &PfMuE, &b_PfMuE);
   fChain->SetBranchAddress("PfMuEt", &PfMuEt, &b_PfMuEt);
   fChain->SetBranchAddress("PfMuEta", &PfMuEta, &b_PfMuEta);
   fChain->SetBranchAddress("PfMuPhi", &PfMuPhi, &b_PfMuPhi);
   fChain->SetBranchAddress("PfMuCharge", &PfMuCharge, &b_PfMuCharge);
   fChain->SetBranchAddress("PfMuParticleIso", &PfMuParticleIso, &b_PfMuParticleIso);
   fChain->SetBranchAddress("PfMuChargedHadronIso", &PfMuChargedHadronIso, &b_PfMuChargedHadronIso);
   fChain->SetBranchAddress("PfMuNeutralHadronIso", &PfMuNeutralHadronIso, &b_PfMuNeutralHadronIso);
   fChain->SetBranchAddress("PfMuPhotonIso", &PfMuPhotonIso, &b_PfMuPhotonIso);
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
   fChain->SetBranchAddress("ElSCEta", ElSCEta, &b_ElSCEta);
   fChain->SetBranchAddress("ElPhi", ElPhi, &b_ElPhi);
   fChain->SetBranchAddress("ElGsfTkPt", ElGsfTkPt, &b_ElGsfTkPt);
   fChain->SetBranchAddress("ElGsfTkEta", ElGsfTkEta, &b_ElGsfTkEta);
   fChain->SetBranchAddress("ElGsfTkPhi", ElGsfTkPhi, &b_ElGsfTkPhi);
   fChain->SetBranchAddress("ElTrkMomentumError", ElTrkMomentumError, &b_ElTrkMomentumError);
   fChain->SetBranchAddress("ElEcalEnergyError", ElEcalEnergyError, &b_ElEcalEnergyError);
   fChain->SetBranchAddress("ElEleMomentumError", ElEleMomentumError, &b_ElEleMomentumError);
   fChain->SetBranchAddress("ElNBrems", ElNBrems, &b_ElNBrems);
   fChain->SetBranchAddress("ElD0BS", ElD0BS, &b_ElD0BS);
   fChain->SetBranchAddress("ElD0PV", ElD0PV, &b_ElD0PV);
   fChain->SetBranchAddress("ElD0E", ElD0E, &b_ElD0E);
   fChain->SetBranchAddress("ElDzBS", ElDzBS, &b_ElDzBS);
   fChain->SetBranchAddress("ElDzPV", ElDzPV, &b_ElDzPV);
   fChain->SetBranchAddress("ElDzE", ElDzE, &b_ElDzE);
   fChain->SetBranchAddress("ElRelIso03", ElRelIso03, &b_ElRelIso03);
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
   fChain->SetBranchAddress("ElIDMva", ElIDMva, &b_ElIDMva);
   fChain->SetBranchAddress("ElIDTight", ElIDTight, &b_ElIDTight);
   fChain->SetBranchAddress("ElIDLoose", ElIDLoose, &b_ElIDLoose);
   fChain->SetBranchAddress("ElIDRobustTight", ElIDRobustTight, &b_ElIDRobustTight);
   fChain->SetBranchAddress("ElIDRobustLoose", ElIDRobustLoose, &b_ElIDRobustLoose);
   fChain->SetBranchAddress("ElIDsimpleWPrelIso", ElIDsimpleWPrelIso, &b_ElIDsimpleWPrelIso);
   fChain->SetBranchAddress("ElIDsimpleWP80relIso", ElIDsimpleWP80relIso, &b_ElIDsimpleWP80relIso);
   fChain->SetBranchAddress("ElIDsimpleWP85relIso", ElIDsimpleWP85relIso, &b_ElIDsimpleWP85relIso);
   fChain->SetBranchAddress("ElIDsimpleWP90relIso", ElIDsimpleWP90relIso, &b_ElIDsimpleWP90relIso);
   fChain->SetBranchAddress("ElIDsimpleWP95relIso", ElIDsimpleWP95relIso, &b_ElIDsimpleWP95relIso);
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
   fChain->SetBranchAddress("ElNumberOfMissingInnerHits", ElNumberOfMissingInnerHits, &b_ElNumberOfMissingInnerHits);
   fChain->SetBranchAddress("ElConvPartnerTrkDist", ElConvPartnerTrkDist, &b_ElConvPartnerTrkDist);
   fChain->SetBranchAddress("ElConvPartnerTrkDCot", ElConvPartnerTrkDCot, &b_ElConvPartnerTrkDCot);
   fChain->SetBranchAddress("ElConvPartnerTrkPt", ElConvPartnerTrkPt, &b_ElConvPartnerTrkPt);
   fChain->SetBranchAddress("ElConvPartnerTrkEta", ElConvPartnerTrkEta, &b_ElConvPartnerTrkEta);
   fChain->SetBranchAddress("ElConvPartnerTrkPhi", ElConvPartnerTrkPhi, &b_ElConvPartnerTrkPhi);
   fChain->SetBranchAddress("ElConvPartnerTrkCharge", ElConvPartnerTrkCharge, &b_ElConvPartnerTrkCharge);
   fChain->SetBranchAddress("ElScSeedSeverity", ElScSeedSeverity, &b_ElScSeedSeverity);
   fChain->SetBranchAddress("ElE1OverE9", ElE1OverE9, &b_ElE1OverE9);
   fChain->SetBranchAddress("ElS4OverS1", ElS4OverS1, &b_ElS4OverS1);
   fChain->SetBranchAddress("ElGenID", ElGenID, &b_ElGenID);
   fChain->SetBranchAddress("ElGenStatus", ElGenStatus, &b_ElGenStatus);
   fChain->SetBranchAddress("ElGenPt", ElGenPt, &b_ElGenPt);
   fChain->SetBranchAddress("ElGenEta", ElGenEta, &b_ElGenEta);
   fChain->SetBranchAddress("ElGenPhi", ElGenPhi, &b_ElGenPhi);
   fChain->SetBranchAddress("ElGenE", ElGenE, &b_ElGenE);
   fChain->SetBranchAddress("ElGenMID", ElGenMID, &b_ElGenMID);
   fChain->SetBranchAddress("ElGenMStatus", ElGenMStatus, &b_ElGenMStatus);
   fChain->SetBranchAddress("ElGenMPt", ElGenMPt, &b_ElGenMPt);
   fChain->SetBranchAddress("ElGenMEta", ElGenMEta, &b_ElGenMEta);
   fChain->SetBranchAddress("ElGenMPhi", ElGenMPhi, &b_ElGenMPhi);
   fChain->SetBranchAddress("ElGenME", ElGenME, &b_ElGenME);
   fChain->SetBranchAddress("ElGenGMID", ElGenGMID, &b_ElGenGMID);
   fChain->SetBranchAddress("ElGenGMStatus", ElGenGMStatus, &b_ElGenGMStatus);
   fChain->SetBranchAddress("ElGenGMPt", ElGenGMPt, &b_ElGenGMPt);
   fChain->SetBranchAddress("ElGenGMEta", ElGenGMEta, &b_ElGenGMEta);
   fChain->SetBranchAddress("ElGenGMPhi", ElGenGMPhi, &b_ElGenGMPhi);
   fChain->SetBranchAddress("ElGenGME", ElGenGME, &b_ElGenGME);
   fChain->SetBranchAddress("NPfEls", &NPfEls, &b_NPfEls);
   fChain->SetBranchAddress("NPfElsTot", &NPfElsTot, &b_NPfElsTot);
   fChain->SetBranchAddress("PfElPx", &PfElPx, &b_PfElPx);
   fChain->SetBranchAddress("PfElPy", &PfElPy, &b_PfElPy);
   fChain->SetBranchAddress("PfElPz", &PfElPz, &b_PfElPz);
   fChain->SetBranchAddress("PfElPt", &PfElPt, &b_PfElPt);
   fChain->SetBranchAddress("PfElPtE", &PfElPtE, &b_PfElPtE);
   fChain->SetBranchAddress("PfElE", &PfElE, &b_PfElE);
   fChain->SetBranchAddress("PfElEt", &PfElEt, &b_PfElEt);
   fChain->SetBranchAddress("PfElEta", &PfElEta, &b_PfElEta);
   fChain->SetBranchAddress("PfElPhi", &PfElPhi, &b_PfElPhi);
   fChain->SetBranchAddress("PfElCharge", &PfElCharge, &b_PfElCharge);
   fChain->SetBranchAddress("PfElParticleIso", &PfElParticleIso, &b_PfElParticleIso);
   fChain->SetBranchAddress("PfElChargedHadronIso", &PfElChargedHadronIso, &b_PfElChargedHadronIso);
   fChain->SetBranchAddress("PfElNeutralHadronIso", &PfElNeutralHadronIso, &b_PfElNeutralHadronIso);
   fChain->SetBranchAddress("PfElPhotonIso", &PfElPhotonIso, &b_PfElPhotonIso);
   fChain->SetBranchAddress("NPfTaus", &NPfTaus, &b_NPfTaus);
   fChain->SetBranchAddress("NPfTausTot", &NPfTausTot, &b_NPfTausTot);
   fChain->SetBranchAddress("PfTauPx", PfTauPx, &b_PfTauPx);
   fChain->SetBranchAddress("PfTauPy", PfTauPy, &b_PfTauPy);
   fChain->SetBranchAddress("PfTauPz", PfTauPz, &b_PfTauPz);
   fChain->SetBranchAddress("PfTauPt", PfTauPt, &b_PfTauPt);
   fChain->SetBranchAddress("PfTauPtE", PfTauPtE, &b_PfTauPtE);
   fChain->SetBranchAddress("PfTauE", PfTauE, &b_PfTauE);
   fChain->SetBranchAddress("PfTauEt", PfTauEt, &b_PfTauEt);
   fChain->SetBranchAddress("PfTauEta", PfTauEta, &b_PfTauEta);
   fChain->SetBranchAddress("PfTauPhi", PfTauPhi, &b_PfTauPhi);
   fChain->SetBranchAddress("PfTauCharge", PfTauCharge, &b_PfTauCharge);
   fChain->SetBranchAddress("PfTauParticleIso", PfTauParticleIso, &b_PfTauParticleIso);
   fChain->SetBranchAddress("PfTauChargedHadronIso", PfTauChargedHadronIso, &b_PfTauChargedHadronIso);
   fChain->SetBranchAddress("PfTauNeutralHadronIso", PfTauNeutralHadronIso, &b_PfTauNeutralHadronIso);
   fChain->SetBranchAddress("PfTauPhotonIso", PfTauPhotonIso, &b_PfTauPhotonIso);
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
   fChain->SetBranchAddress("PhoIso04Ecal", PhoIso04Ecal, &b_PhoIso04Ecal);
   fChain->SetBranchAddress("PhoIso04Hcal", PhoIso04Hcal, &b_PhoIso04Hcal);
   fChain->SetBranchAddress("PhoIso04TrkSolid", PhoIso04TrkSolid, &b_PhoIso04TrkSolid);
   fChain->SetBranchAddress("PhoIso04TrkHollow", PhoIso04TrkHollow, &b_PhoIso04TrkHollow);
   fChain->SetBranchAddress("PhoIso04", PhoIso04, &b_PhoIso04);
   fChain->SetBranchAddress("PhoR9", PhoR9, &b_PhoR9);
   fChain->SetBranchAddress("PhoCaloPositionX", PhoCaloPositionX, &b_PhoCaloPositionX);
   fChain->SetBranchAddress("PhoCaloPositionY", PhoCaloPositionY, &b_PhoCaloPositionY);
   fChain->SetBranchAddress("PhoCaloPositionZ", PhoCaloPositionZ, &b_PhoCaloPositionZ);
   fChain->SetBranchAddress("PhoHoverE", PhoHoverE, &b_PhoHoverE);
   fChain->SetBranchAddress("PhoH1overE", PhoH1overE, &b_PhoH1overE);
   fChain->SetBranchAddress("PhoH2overE", PhoH2overE, &b_PhoH2overE);
   fChain->SetBranchAddress("PhoSigmaIetaIeta", PhoSigmaIetaIeta, &b_PhoSigmaIetaIeta);
   fChain->SetBranchAddress("PhoSCRawEnergy", PhoSCRawEnergy, &b_PhoSCRawEnergy);
   fChain->SetBranchAddress("PhoSCEtaWidth", PhoSCEtaWidth, &b_PhoSCEtaWidth);
   fChain->SetBranchAddress("PhoSCSigmaPhiPhi", PhoSCSigmaPhiPhi, &b_PhoSCSigmaPhiPhi);
   fChain->SetBranchAddress("PhoHasPixSeed", PhoHasPixSeed, &b_PhoHasPixSeed);
   fChain->SetBranchAddress("PhoHasConvTrks", PhoHasConvTrks, &b_PhoHasConvTrks);
   fChain->SetBranchAddress("PhoScSeedSeverity", PhoScSeedSeverity, &b_PhoScSeedSeverity);
   fChain->SetBranchAddress("PhoE1OverE9", PhoE1OverE9, &b_PhoE1OverE9);
   fChain->SetBranchAddress("PhoS4OverS1", PhoS4OverS1, &b_PhoS4OverS1);
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
   fChain->SetBranchAddress("JEcorr", JEcorr, &b_JEcorr);
   fChain->SetBranchAddress("JEtaRms", JEtaRms, &b_JEtaRms);
   fChain->SetBranchAddress("JPhiRms", JPhiRms, &b_JPhiRms);
   fChain->SetBranchAddress("JNConstituents", JNConstituents, &b_JNConstituents);
   fChain->SetBranchAddress("JNAssoTracks", JNAssoTracks, &b_JNAssoTracks);
   fChain->SetBranchAddress("JNNeutrals", JNNeutrals, &b_JNNeutrals);
   fChain->SetBranchAddress("JChargedEmFrac", JChargedEmFrac, &b_JChargedEmFrac);
   fChain->SetBranchAddress("JNeutralEmFrac", JNeutralEmFrac, &b_JNeutralEmFrac);
   fChain->SetBranchAddress("JChargedHadFrac", JChargedHadFrac, &b_JChargedHadFrac);
   fChain->SetBranchAddress("JNeutralHadFrac", JNeutralHadFrac, &b_JNeutralHadFrac);
   fChain->SetBranchAddress("JChargedMuEnergyFrac", JChargedMuEnergyFrac, &b_JChargedMuEnergyFrac);
   fChain->SetBranchAddress("JeMinDR", JeMinDR, &b_JeMinDR);
   fChain->SetBranchAddress("JbTagProbTkCntHighEff", JbTagProbTkCntHighEff, &b_JbTagProbTkCntHighEff);
   fChain->SetBranchAddress("JbTagProbTkCntHighPur", JbTagProbTkCntHighPur, &b_JbTagProbTkCntHighPur);
   fChain->SetBranchAddress("JbTagProbSimpSVHighEff", JbTagProbSimpSVHighEff, &b_JbTagProbSimpSVHighEff);
   fChain->SetBranchAddress("JbTagProbSimpSVHighPur", JbTagProbSimpSVHighPur, &b_JbTagProbSimpSVHighPur);
   fChain->SetBranchAddress("JMass", JMass, &b_JMass);
   fChain->SetBranchAddress("Jtrk1px", Jtrk1px, &b_Jtrk1px);
   fChain->SetBranchAddress("Jtrk1py", Jtrk1py, &b_Jtrk1py);
   fChain->SetBranchAddress("Jtrk1pz", Jtrk1pz, &b_Jtrk1pz);
   fChain->SetBranchAddress("Jtrk2px", Jtrk2px, &b_Jtrk2px);
   fChain->SetBranchAddress("Jtrk2py", Jtrk2py, &b_Jtrk2py);
   fChain->SetBranchAddress("Jtrk2pz", Jtrk2pz, &b_Jtrk2pz);
   fChain->SetBranchAddress("Jtrk3px", Jtrk3px, &b_Jtrk3px);
   fChain->SetBranchAddress("Jtrk3py", Jtrk3py, &b_Jtrk3py);
   fChain->SetBranchAddress("Jtrk3pz", Jtrk3pz, &b_Jtrk3pz);
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
   fChain->SetBranchAddress("CANJets", &CANJets, &b_CANJets);
   fChain->SetBranchAddress("CAJPx", CAJPx, &b_CAJPx);
   fChain->SetBranchAddress("CAJPy", CAJPy, &b_CAJPy);
   fChain->SetBranchAddress("CAJPz", CAJPz, &b_CAJPz);
   fChain->SetBranchAddress("CAJPt", CAJPt, &b_CAJPt);
   fChain->SetBranchAddress("CAJE", CAJE, &b_CAJE);
   fChain->SetBranchAddress("CAJEt", CAJEt, &b_CAJEt);
   fChain->SetBranchAddress("CAJEta", CAJEta, &b_CAJEta);
   fChain->SetBranchAddress("CAJPhi", CAJPhi, &b_CAJPhi);
   fChain->SetBranchAddress("CAJScale", CAJScale, &b_CAJScale);
   fChain->SetBranchAddress("CAJbTagProbTkCntHighEff", CAJbTagProbTkCntHighEff, &b_CAJbTagProbTkCntHighEff);
   fChain->SetBranchAddress("CAJbTagProbTkCntHighPur", CAJbTagProbTkCntHighPur, &b_CAJbTagProbTkCntHighPur);
   fChain->SetBranchAddress("CAJbTagProbSimpSVHighEff", CAJbTagProbSimpSVHighEff, &b_CAJbTagProbSimpSVHighEff);
   fChain->SetBranchAddress("CAJbTagProbSimpSVHighPur", CAJbTagProbSimpSVHighPur, &b_CAJbTagProbSimpSVHighPur);
   fChain->SetBranchAddress("CAJID_HPD", CAJID_HPD, &b_CAJID_HPD);
   fChain->SetBranchAddress("CAJID_RBX", CAJID_RBX, &b_CAJID_RBX);
   fChain->SetBranchAddress("CAJID_n90Hits", CAJID_n90Hits, &b_CAJID_n90Hits);
   fChain->SetBranchAddress("CAJID_resEMF", CAJID_resEMF, &b_CAJID_resEMF);
   fChain->SetBranchAddress("CAJEMfrac", CAJEMfrac, &b_CAJEMfrac);
   fChain->SetBranchAddress("CAJNAssoTracks", CAJNAssoTracks, &b_CAJNAssoTracks);
   fChain->SetBranchAddress("CAJChfrac", CAJChfrac, &b_CAJChfrac);
   fChain->SetBranchAddress("CAJNConstituents", CAJNConstituents, &b_CAJNConstituents);
   fChain->SetBranchAddress("PF2PATNJets", &PF2PATNJets, &b_PF2PATNJets);
   fChain->SetBranchAddress("PF2PATJPx", PF2PATJPx, &b_PF2PATJPx);
   fChain->SetBranchAddress("PF2PATJPy", PF2PATJPy, &b_PF2PATJPy);
   fChain->SetBranchAddress("PF2PATJPz", PF2PATJPz, &b_PF2PATJPz);
   fChain->SetBranchAddress("PF2PATJPt", PF2PATJPt, &b_PF2PATJPt);
   fChain->SetBranchAddress("PF2PATJE", PF2PATJE, &b_PF2PATJE);
   fChain->SetBranchAddress("PF2PATJEt", PF2PATJEt, &b_PF2PATJEt);
   fChain->SetBranchAddress("PF2PATJEta", PF2PATJEta, &b_PF2PATJEta);
   fChain->SetBranchAddress("PF2PATJPhi", PF2PATJPhi, &b_PF2PATJPhi);
   fChain->SetBranchAddress("PF2PATJScale", PF2PATJScale, &b_PF2PATJScale);
   fChain->SetBranchAddress("PF2PATJbTagProbTkCntHighEff", PF2PATJbTagProbTkCntHighEff, &b_PF2PATJbTagProbTkCntHighEff);
   fChain->SetBranchAddress("PF2PATJbTagProbTkCntHighPur", PF2PATJbTagProbTkCntHighPur, &b_PF2PATJbTagProbTkCntHighPur);
   fChain->SetBranchAddress("PF2PATJbTagProbSimpSVHighEff", PF2PATJbTagProbSimpSVHighEff, &b_PF2PATJbTagProbSimpSVHighEff);
   fChain->SetBranchAddress("PF2PATJbTagProbSimpSVHighPur", PF2PATJbTagProbSimpSVHighPur, &b_PF2PATJbTagProbSimpSVHighPur);
   fChain->SetBranchAddress("PF2PATJChMult", PF2PATJChMult, &b_PF2PATJChMult);
   fChain->SetBranchAddress("PF2PATJNeuMult", PF2PATJNeuMult, &b_PF2PATJNeuMult);
   fChain->SetBranchAddress("PF2PATJChHadfrac", PF2PATJChHadfrac, &b_PF2PATJChHadfrac);
   fChain->SetBranchAddress("PF2PATJNeuHadfrac", PF2PATJNeuHadfrac, &b_PF2PATJNeuHadfrac);
   fChain->SetBranchAddress("PF2PATJChEmfrac", PF2PATJChEmfrac, &b_PF2PATJChEmfrac);
   fChain->SetBranchAddress("PF2PATJNeuEmfrac", PF2PATJNeuEmfrac, &b_PF2PATJNeuEmfrac);
   fChain->SetBranchAddress("PF2PATJNConstituents", PF2PATJNConstituents, &b_PF2PATJNConstituents);
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
   fChain->SetBranchAddress("GenMET", &GenMET, &b_GenMET);
   fChain->SetBranchAddress("GenMETpx", &GenMETpx, &b_GenMETpx);
   fChain->SetBranchAddress("GenMETpy", &GenMETpy, &b_GenMETpy);
   fChain->SetBranchAddress("GenMETphi", &GenMETphi, &b_GenMETphi);
   fChain->SetBranchAddress("TCMET", &TCMET, &b_TCMET);
   fChain->SetBranchAddress("TCMETpx", &TCMETpx, &b_TCMETpx);
   fChain->SetBranchAddress("TCMETpy", &TCMETpy, &b_TCMETpy);
   fChain->SetBranchAddress("TCMETphi", &TCMETphi, &b_TCMETphi);
   fChain->SetBranchAddress("TCMETSignificance", &TCMETSignificance, &b_TCMETSignificance);
   fChain->SetBranchAddress("MuJESCorrMET", &MuJESCorrMET, &b_MuJESCorrMET);
   fChain->SetBranchAddress("MuJESCorrMETpx", &MuJESCorrMETpx, &b_MuJESCorrMETpx);
   fChain->SetBranchAddress("MuJESCorrMETpy", &MuJESCorrMETpy, &b_MuJESCorrMETpy);
   fChain->SetBranchAddress("MuJESCorrMETphi", &MuJESCorrMETphi, &b_MuJESCorrMETphi);
   fChain->SetBranchAddress("PFMET", &PFMET, &b_PFMET);
   fChain->SetBranchAddress("PFMETpx", &PFMETpx, &b_PFMETpx);
   fChain->SetBranchAddress("PFMETpy", &PFMETpy, &b_PFMETpy);
   fChain->SetBranchAddress("PFMETphi", &PFMETphi, &b_PFMETphi);
   fChain->SetBranchAddress("PFMETSignificance", &PFMETSignificance, &b_PFMETSignificance);
   fChain->SetBranchAddress("METR12", &METR12, &b_METR12);
   fChain->SetBranchAddress("METR21", &METR21, &b_METR21);
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
