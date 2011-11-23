//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Nov 23 19:23:26 2011 by ROOT version 5.27/06b
// from TTree Analysis/ETHZAnalysisTree
// found on file: dcap://t3se01.psi.ch:22125/pnfs/psi.ch/cms/trivcat/store/user/susy/ntuples/mc/V02-05-02/G_Pt-300to470_Tune23_7TeV_herwigpp_Summer11_PU_S3/NTupleProducer_42X_MC_259_1_Le7.root
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
   Float_t         QCDPartonicHT;
   Int_t           SigProcID;
   Float_t         PDFScalePDF;
   Int_t           PDFID1;
   Int_t           PDFID2;
   Float_t         PDFx1;
   Float_t         PDFx2;
   Float_t         PDFxPDF1;
   Float_t         PDFxPDF2;
   Float_t         GenWeight;
   Float_t         ExtXSecLO;
   Float_t         IntXSec;
   Float_t         pdfW[100];
   Float_t         pdfWsum;
   Int_t           NPdfs;
   Int_t           PUnumInteractions;
   Int_t           PUnumFilled;
   Int_t           PUOOTnumInteractionsEarly;
   Int_t           PUOOTnumInteractionsLate;
   Float_t         PUzPositions[50];   //[PUnumFilled]
   Float_t         PUsumPtLowPt[50];   //[PUnumFilled]
   Float_t         PUsumPtHighPt[50];   //[PUnumFilled]
   Float_t         PUnTrksLowPt[50];   //[PUnumFilled]
   Float_t         PUnTrksHighPt[50];   //[PUnumFilled]
   Float_t         Rho;
   Float_t         RhoPFnoPU;
   Float_t         Weight;
   Int_t           HLTResults[500];
   Int_t           HLTPrescale[500];
   Int_t           L1PhysResults[128];
   Int_t           L1TechResults[64];
   Int_t           NHLTObjs;
   Int_t           HLTObjectID[7][10];
   Float_t         HLTObjectPt[7][10];
   Float_t         HLTObjectEta[7][10];
   Float_t         HLTObjectPhi[7][10];
   Float_t         PUWeightTotal;
   Float_t         PUWeightInTime;
   Float_t         MassGlu;
   Float_t         MassChi;
   Float_t         MassLSP;
   Float_t         M0;
   Float_t         M12;
   Float_t         signMu;
   Float_t         A0;
   Int_t           process;
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
   Int_t           MaxGenPhoExceed;
   Int_t           MaxGenJetExceed;
   Int_t           MaxVerticesExceed;
   Int_t           HBHENoiseFlag;
   Int_t           HBHENoiseFlagIso;
   Int_t           CSCTightHaloID;
   Int_t           EcalDeadTPFilterFlag;
   Int_t           RecovRecHitFilterFlag;
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
   Int_t           NGenPhotons;
   Float_t         GenPhotonPt[100];   //[NGenPhotons]
   Float_t         GenPhotonEta[100];   //[NGenPhotons]
   Float_t         GenPhotonPhi[100];   //[NGenPhotons]
   Int_t           GenPhotonMotherID[100];   //[NGenPhotons]
   Int_t           GenPhotonMotherStatus[100];   //[NGenPhotons]
   Int_t           NGenJets;
   Float_t         GenJetPt[100];   //[NGenJets]
   Float_t         GenJetEta[100];   //[NGenJets]
   Float_t         GenJetPhi[100];   //[NGenJets]
   Float_t         GenJetE[100];   //[NGenJets]
   Float_t         GenJetEmE[100];   //[NGenJets]
   Float_t         GenJetHadE[100];   //[NGenJets]
   Float_t         GenJetInvE[100];   //[NGenJets]
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
   Int_t           NEBhits;
   Float_t         EBrechitE[20];   //[NEBhits]
   Float_t         EBrechitPt[20];   //[NEBhits]
   Float_t         EBrechitEta[20];   //[NEBhits]
   Float_t         EBrechitPhi[20];   //[NEBhits]
   Float_t         EBrechitChi2[20];   //[NEBhits]
   Float_t         EBrechitTime[20];   //[NEBhits]
   Float_t         EBrechitE4oE1[20];   //[NEBhits]
   Float_t         EBrechitE2oE9[20];   //[NEBhits]
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
   Int_t           ElSCindex[20];   //[NEles]
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
   Float_t         PhoSigmaEtaEta[50];   //[NPhotons]
   Float_t         PhoE1x5[50];   //[NPhotons]
   Float_t         PhoE2x5[50];   //[NPhotons]
   Float_t         PhoE3x3[50];   //[NPhotons]
   Float_t         PhoE5x5[50];   //[NPhotons]
   Float_t         PhomaxEnergyXtal[50];   //[NPhotons]
   Float_t         PhoIso03HcalDepth1[50];   //[NPhotons]
   Float_t         PhoIso03HcalDepth2[50];   //[NPhotons]
   Float_t         PhoIso04HcalDepth1[50];   //[NPhotons]
   Float_t         PhoIso04HcalDepth2[50];   //[NPhotons]
   Int_t           PhoIso03nTrksSolid[50];   //[NPhotons]
   Int_t           PhoIso03nTrksHollow[50];   //[NPhotons]
   Int_t           PhoIso04nTrksSolid[50];   //[NPhotons]
   Int_t           PhoIso04nTrksHollow[50];   //[NPhotons]
   Int_t           PhoisEB[50];   //[NPhotons]
   Int_t           PhoisEE[50];   //[NPhotons]
   Int_t           PhoisEBEtaGap[50];   //[NPhotons]
   Int_t           PhoisEBPhiGap[50];   //[NPhotons]
   Int_t           PhoisEERingGap[50];   //[NPhotons]
   Int_t           PhoisEEDeeGap[50];   //[NPhotons]
   Int_t           PhoisEBEEGap[50];   //[NPhotons]
   Int_t           PhoisPFlowPhoton[50];   //[NPhotons]
   Int_t           PhoisStandardPhoton[50];   //[NPhotons]
   Int_t           PhoMCmatchindex[50];   //[NPhotons]
   Int_t           PhoMCmatchexitcode[50];   //[NPhotons]
   Float_t         Pho_ChargedHadronIso[50];   //[NPhotons]
   Float_t         Pho_NeutralHadronIso[50];   //[NPhotons]
   Float_t         Pho_PhotonIso[50];   //[NPhotons]
   Int_t           Pho_isPFPhoton[50];   //[NPhotons]
   Int_t           Pho_isPFElectron[50];   //[NPhotons]
   Int_t           PhotSCindex[50];   //[NPhotons]
   Float_t         Pho_Cone04PhotonIso_dR0_dEta0_pt0[50];   //[NPhotons]
   Float_t         Pho_Cone04PhotonIso_dR0_dEta0_pt5[50];   //[NPhotons]
   Float_t         Pho_Cone04PhotonIso_dR8_dEta0_pt0[50];   //[NPhotons]
   Float_t         Pho_Cone04PhotonIso_dR8_dEta0_pt5[50];   //[NPhotons]
   Float_t         Pho_Cone01PhotonIso_dR045EB070EE_dEta015_pt08EB1EE_mvVtx[50];   //[NPhotons]
   Float_t         Pho_Cone02PhotonIso_dR045EB070EE_dEta015_pt08EB1EE_mvVtx[50];   //[NPhotons]
   Float_t         Pho_Cone03PhotonIso_dR045EB070EE_dEta015_pt08EB1EE_mvVtx[50];   //[NPhotons]
   Float_t         Pho_Cone04PhotonIso_dR045EB070EE_dEta015_pt08EB1EE_mvVtx[50];   //[NPhotons]
   Float_t         Pho_Cone04NeutralHadronIso_dR0_dEta0_pt0[50];   //[NPhotons]
   Float_t         Pho_Cone04NeutralHadronIso_dR0_dEta0_pt5[50];   //[NPhotons]
   Float_t         Pho_Cone04NeutralHadronIso_dR0_dEta0_pt0_nocracks[50];   //[NPhotons]
   Float_t         Pho_Cone04NeutralHadronIso_dR0_dEta0_pt5_nocracks[50];   //[NPhotons]
   Float_t         Pho_Cone04NeutralHadronIso_dR7_dEta0_pt0[50];   //[NPhotons]
   Float_t         Pho_Cone04NeutralHadronIso_dR7_dEta0_pt5[50];   //[NPhotons]
   Float_t         Pho_Cone01NeutralHadronIso_dR0_dEta0_pt0_mvVtx[50];   //[NPhotons]
   Float_t         Pho_Cone02NeutralHadronIso_dR0_dEta0_pt0_mvVtx[50];   //[NPhotons]
   Float_t         Pho_Cone03NeutralHadronIso_dR0_dEta0_pt0_mvVtx[50];   //[NPhotons]
   Float_t         Pho_Cone04NeutralHadronIso_dR0_dEta0_pt0_mvVtx[50];   //[NPhotons]
   Float_t         Pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_dz0_old[50];   //[NPhotons]
   Float_t         Pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_PFnoPU_old[50];   //[NPhotons]
   Float_t         Pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_dz0_old[50];   //[NPhotons]
   Float_t         Pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_PFnoPU_old[50];   //[NPhotons]
   Float_t         Pho_Cone01ChargedHadronIso_dR0_dEta0_pt0_dz0[50];   //[NPhotons]
   Float_t         Pho_Cone01ChargedHadronIso_dR0_dEta0_pt0_dz1_dxy01[50];   //[NPhotons]
   Float_t         Pho_Cone01ChargedHadronIso_dR0_dEta0_pt0_PFnoPU[50];   //[NPhotons]
   Float_t         Pho_Cone01ChargedHadronIso_dR015_dEta0_pt0_dz0[50];   //[NPhotons]
   Float_t         Pho_Cone01ChargedHadronIso_dR015_dEta0_pt0_dz1_dxy01[50];   //[NPhotons]
   Float_t         Pho_Cone01ChargedHadronIso_dR015_dEta0_pt0_PFnoPU[50];   //[NPhotons]
   Float_t         Pho_Cone02ChargedHadronIso_dR0_dEta0_pt0_dz0[50];   //[NPhotons]
   Float_t         Pho_Cone02ChargedHadronIso_dR0_dEta0_pt0_dz1_dxy01[50];   //[NPhotons]
   Float_t         Pho_Cone02ChargedHadronIso_dR0_dEta0_pt0_PFnoPU[50];   //[NPhotons]
   Float_t         Pho_Cone02ChargedHadronIso_dR015_dEta0_pt0_dz0[50];   //[NPhotons]
   Float_t         Pho_Cone02ChargedHadronIso_dR015_dEta0_pt0_dz1_dxy01[50];   //[NPhotons]
   Float_t         Pho_Cone02ChargedHadronIso_dR015_dEta0_pt0_PFnoPU[50];   //[NPhotons]
   Float_t         Pho_Cone03ChargedHadronIso_dR0_dEta0_pt0_dz0[50];   //[NPhotons]
   Float_t         Pho_Cone03ChargedHadronIso_dR0_dEta0_pt0_dz1_dxy01[50];   //[NPhotons]
   Float_t         Pho_Cone03ChargedHadronIso_dR0_dEta0_pt0_PFnoPU[50];   //[NPhotons]
   Float_t         Pho_Cone03ChargedHadronIso_dR015_dEta0_pt0_dz0[50];   //[NPhotons]
   Float_t         Pho_Cone03ChargedHadronIso_dR015_dEta0_pt0_dz1_dxy01[50];   //[NPhotons]
   Float_t         Pho_Cone03ChargedHadronIso_dR015_dEta0_pt0_PFnoPU[50];   //[NPhotons]
   Float_t         Pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_dz0[50];   //[NPhotons]
   Float_t         Pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_dz1_dxy01[50];   //[NPhotons]
   Float_t         Pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_PFnoPU[50];   //[NPhotons]
   Float_t         Pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_dz0[50];   //[NPhotons]
   Float_t         Pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_dz1_dxy01[50];   //[NPhotons]
   Float_t         Pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_PFnoPU[50];   //[NPhotons]
   Int_t           NSuperClusters;
   Float_t         SCRaw[100];   //[NSuperClusters]
   Float_t         SCPre[100];   //[NSuperClusters]
   Float_t         SCEnergy[100];   //[NSuperClusters]
   Float_t         SCEta[100];   //[NSuperClusters]
   Float_t         SCPhi[100];   //[NSuperClusters]
   Float_t         SCPhiWidth[100];   //[NSuperClusters]
   Float_t         SCEtaWidth[100];   //[NSuperClusters]
   Float_t         SCBrem[100];   //[NSuperClusters]
   Float_t         SCR9[100];   //[NSuperClusters]
   Float_t         SCcrackcorrseed[100];   //[NSuperClusters]
   Float_t         SCcrackcorr[100];   //[NSuperClusters]
   Float_t         SClocalcorrseed[100];   //[NSuperClusters]
   Float_t         SClocalcorr[100];   //[NSuperClusters]
   Float_t         SCcrackcorrseedfactor[100];   //[NSuperClusters]
   Float_t         SClocalcorrseedfactor[100];   //[NSuperClusters]
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
   Float_t         JArea[100];   //[NJets]
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
   Int_t           JGenJetIndex[100];   //[NJets]
   Int_t           CANJets;
   Double_t        CAJPx[100];   //[CANJets]
   Double_t        CAJPy[100];   //[CANJets]
   Double_t        CAJPz[100];   //[CANJets]
   Double_t        CAJPt[100];   //[CANJets]
   Double_t        CAJE[100];   //[CANJets]
   Double_t        CAJEt[100];   //[CANJets]
   Double_t        CAJEta[100];   //[CANJets]
   Double_t        CAJPhi[100];   //[CANJets]
   Int_t           CAJFlavour[100];   //[CANJets]
   Double_t        CAJScale[100];   //[CANJets]
   Double_t        CAJL1FastJetScale[100];   //[CANJets]
   Double_t        CAJArea[100];   //[CANJets]
   Double_t        CAJbTagProbTkCntHighEff[100];   //[CANJets]
   Double_t        CAJbTagProbTkCntHighPur[100];   //[CANJets]
   Double_t        CAJbTagProbSimpSVHighEff[100];   //[CANJets]
   Double_t        CAJbTagProbSimpSVHighPur[100];   //[CANJets]
   Int_t           CAJIDLoose[100];   //[CANJets]
   Double_t        CAJID_HPD[100];   //[CANJets]
   Double_t        CAJID_RBX[100];   //[CANJets]
   Double_t        CAJID_n90Hits[100];   //[CANJets]
   Double_t        CAJID_resEMF[100];   //[CANJets]
   Double_t        CAJEMfrac[100];   //[CANJets]
   Int_t           CAJNAssoTracks[100];   //[CANJets]
   Double_t        CAJChfrac[100];   //[CANJets]
   Int_t           CAJNConstituents[100];   //[CANJets]
   Int_t           CAJn90[100];   //[CANJets]
   Int_t           PF2PATAntiIsoNJets;
   Double_t        PF2PATAntiIsoJPx[100];   //[PF2PATAntiIsoNJets]
   Double_t        PF2PATAntiIsoJPy[100];   //[PF2PATAntiIsoNJets]
   Double_t        PF2PATAntiIsoJPz[100];   //[PF2PATAntiIsoNJets]
   Double_t        PF2PATAntiIsoJPt[100];   //[PF2PATAntiIsoNJets]
   Double_t        PF2PATAntiIsoJE[100];   //[PF2PATAntiIsoNJets]
   Double_t        PF2PATAntiIsoJEt[100];   //[PF2PATAntiIsoNJets]
   Double_t        PF2PATAntiIsoJEta[100];   //[PF2PATAntiIsoNJets]
   Double_t        PF2PATAntiIsoJPhi[100];   //[PF2PATAntiIsoNJets]
   Int_t           PF2PATAntiIsoJFlavour[100];   //[PF2PATAntiIsoNJets]
   Double_t        PF2PATAntiIsoJScale[100];   //[PF2PATAntiIsoNJets]
   Double_t        PF2PATAntiIsoJL1FastJetScale[100];   //[PF2PATAntiIsoNJets]
   Double_t        PF2PATAntiIsoJArea[100];   //[PF2PATAntiIsoNJets]
   Double_t        PF2PATAntiIsoJbTagProbTkCntHighEff[100];   //[PF2PATAntiIsoNJets]
   Double_t        PF2PATAntiIsoJbTagProbTkCntHighPur[100];   //[PF2PATAntiIsoNJets]
   Double_t        PF2PATAntiIsoJbTagProbSimpSVHighEff[100];   //[PF2PATAntiIsoNJets]
   Double_t        PF2PATAntiIsoJbTagProbSimpSVHighPur[100];   //[PF2PATAntiIsoNJets]
   Int_t           PF2PATAntiIsoJIDLoose[100];   //[PF2PATAntiIsoNJets]
   Int_t           PF2PATAntiIsoJChMult[100];   //[PF2PATAntiIsoNJets]
   Int_t           PF2PATAntiIsoJNeuMult[100];   //[PF2PATAntiIsoNJets]
   Double_t        PF2PATAntiIsoJChHadfrac[100];   //[PF2PATAntiIsoNJets]
   Double_t        PF2PATAntiIsoJNeuHadfrac[100];   //[PF2PATAntiIsoNJets]
   Double_t        PF2PATAntiIsoJChEmfrac[100];   //[PF2PATAntiIsoNJets]
   Double_t        PF2PATAntiIsoJNeuEmfrac[100];   //[PF2PATAntiIsoNJets]
   Double_t        PF2PATAntiIsoJChMufrac[100];   //[PF2PATAntiIsoNJets]
   Int_t           PF2PATAntiIsoJNConstituents[100];   //[PF2PATAntiIsoNJets]
   Int_t           PF2PAT2NJets;
   Double_t        PF2PAT2JPx[100];   //[PF2PAT2NJets]
   Double_t        PF2PAT2JPy[100];   //[PF2PAT2NJets]
   Double_t        PF2PAT2JPz[100];   //[PF2PAT2NJets]
   Double_t        PF2PAT2JPt[100];   //[PF2PAT2NJets]
   Double_t        PF2PAT2JE[100];   //[PF2PAT2NJets]
   Double_t        PF2PAT2JEt[100];   //[PF2PAT2NJets]
   Double_t        PF2PAT2JEta[100];   //[PF2PAT2NJets]
   Double_t        PF2PAT2JPhi[100];   //[PF2PAT2NJets]
   Int_t           PF2PAT2JFlavour[100];   //[PF2PAT2NJets]
   Double_t        PF2PAT2JScale[100];   //[PF2PAT2NJets]
   Double_t        PF2PAT2JL1FastJetScale[100];   //[PF2PAT2NJets]
   Double_t        PF2PAT2JArea[100];   //[PF2PAT2NJets]
   Double_t        PF2PAT2JbTagProbTkCntHighEff[100];   //[PF2PAT2NJets]
   Double_t        PF2PAT2JbTagProbTkCntHighPur[100];   //[PF2PAT2NJets]
   Double_t        PF2PAT2JbTagProbSimpSVHighEff[100];   //[PF2PAT2NJets]
   Double_t        PF2PAT2JbTagProbSimpSVHighPur[100];   //[PF2PAT2NJets]
   Int_t           PF2PAT2JIDLoose[100];   //[PF2PAT2NJets]
   Int_t           PF2PAT2JChMult[100];   //[PF2PAT2NJets]
   Int_t           PF2PAT2JNeuMult[100];   //[PF2PAT2NJets]
   Double_t        PF2PAT2JChHadfrac[100];   //[PF2PAT2NJets]
   Double_t        PF2PAT2JNeuHadfrac[100];   //[PF2PAT2NJets]
   Double_t        PF2PAT2JChEmfrac[100];   //[PF2PAT2NJets]
   Double_t        PF2PAT2JNeuEmfrac[100];   //[PF2PAT2NJets]
   Double_t        PF2PAT2JChMufrac[100];   //[PF2PAT2NJets]
   Int_t           PF2PAT2JNConstituents[100];   //[PF2PAT2NJets]
   Int_t           PF2PAT3NJets;
   Double_t        PF2PAT3JPx[100];   //[PF2PAT3NJets]
   Double_t        PF2PAT3JPy[100];   //[PF2PAT3NJets]
   Double_t        PF2PAT3JPz[100];   //[PF2PAT3NJets]
   Double_t        PF2PAT3JPt[100];   //[PF2PAT3NJets]
   Double_t        PF2PAT3JE[100];   //[PF2PAT3NJets]
   Double_t        PF2PAT3JEt[100];   //[PF2PAT3NJets]
   Double_t        PF2PAT3JEta[100];   //[PF2PAT3NJets]
   Double_t        PF2PAT3JPhi[100];   //[PF2PAT3NJets]
   Int_t           PF2PAT3JFlavour[100];   //[PF2PAT3NJets]
   Double_t        PF2PAT3JScale[100];   //[PF2PAT3NJets]
   Double_t        PF2PAT3JL1FastJetScale[100];   //[PF2PAT3NJets]
   Double_t        PF2PAT3JArea[100];   //[PF2PAT3NJets]
   Double_t        PF2PAT3JbTagProbTkCntHighEff[100];   //[PF2PAT3NJets]
   Double_t        PF2PAT3JbTagProbTkCntHighPur[100];   //[PF2PAT3NJets]
   Double_t        PF2PAT3JbTagProbSimpSVHighEff[100];   //[PF2PAT3NJets]
   Double_t        PF2PAT3JbTagProbSimpSVHighPur[100];   //[PF2PAT3NJets]
   Int_t           PF2PAT3JIDLoose[100];   //[PF2PAT3NJets]
   Int_t           PF2PAT3JChMult[100];   //[PF2PAT3NJets]
   Int_t           PF2PAT3JNeuMult[100];   //[PF2PAT3NJets]
   Double_t        PF2PAT3JChHadfrac[100];   //[PF2PAT3NJets]
   Double_t        PF2PAT3JNeuHadfrac[100];   //[PF2PAT3NJets]
   Double_t        PF2PAT3JChEmfrac[100];   //[PF2PAT3NJets]
   Double_t        PF2PAT3JNeuEmfrac[100];   //[PF2PAT3NJets]
   Double_t        PF2PAT3JChMufrac[100];   //[PF2PAT3NJets]
   Int_t           PF2PAT3JNConstituents[100];   //[PF2PAT3NJets]
   Int_t           PfMuAntiIsoMaxLepExc;
   Int_t           PfMuAntiIsoNObjsTot;
   Int_t           PfMuAntiIsoNObjs;
   Double_t        PfMuAntiIsoPx[20];   //[PfMuAntiIsoNObjs]
   Double_t        PfMuAntiIsoPy[20];   //[PfMuAntiIsoNObjs]
   Double_t        PfMuAntiIsoPz[20];   //[PfMuAntiIsoNObjs]
   Double_t        PfMuAntiIsoPt[20];   //[PfMuAntiIsoNObjs]
   Double_t        PfMuAntiIsoE[20];   //[PfMuAntiIsoNObjs]
   Double_t        PfMuAntiIsoEt[20];   //[PfMuAntiIsoNObjs]
   Double_t        PfMuAntiIsoEta[20];   //[PfMuAntiIsoNObjs]
   Double_t        PfMuAntiIsoPhi[20];   //[PfMuAntiIsoNObjs]
   Int_t           PfMuAntiIsoCharge[20];   //[PfMuAntiIsoNObjs]
   Double_t        PfMuAntiIsoParticleIso[20];   //[PfMuAntiIsoNObjs]
   Double_t        PfMuAntiIsoChargedHadronIso[20];   //[PfMuAntiIsoNObjs]
   Double_t        PfMuAntiIsoNeutralHadronIso[20];   //[PfMuAntiIsoNObjs]
   Double_t        PfMuAntiIsoPhotonIso[20];   //[PfMuAntiIsoNObjs]
   Double_t        PfMuAntiIsoPtErr[20];   //[PfMuAntiIsoNObjs]
   Int_t           PfMuAntiIsoNMatches[20];   //[PfMuAntiIsoNObjs]
   Int_t           PfMu2MaxLepExc;
   Int_t           PfMu2NObjsTot;
   Int_t           PfMu2NObjs;
   Double_t        PfMu2Px[20];   //[PfMu2NObjs]
   Double_t        PfMu2Py[20];   //[PfMu2NObjs]
   Double_t        PfMu2Pz[20];   //[PfMu2NObjs]
   Double_t        PfMu2Pt[20];   //[PfMu2NObjs]
   Double_t        PfMu2E[20];   //[PfMu2NObjs]
   Double_t        PfMu2Et[20];   //[PfMu2NObjs]
   Double_t        PfMu2Eta[20];   //[PfMu2NObjs]
   Double_t        PfMu2Phi[20];   //[PfMu2NObjs]
   Int_t           PfMu2Charge[20];   //[PfMu2NObjs]
   Double_t        PfMu2ParticleIso[20];   //[PfMu2NObjs]
   Double_t        PfMu2ChargedHadronIso[20];   //[PfMu2NObjs]
   Double_t        PfMu2NeutralHadronIso[20];   //[PfMu2NObjs]
   Double_t        PfMu2PhotonIso[20];   //[PfMu2NObjs]
   Double_t        PfMu2PtErr[20];   //[PfMu2NObjs]
   Int_t           PfMu2NMatches[20];   //[PfMu2NObjs]
   Int_t           PfMu3MaxLepExc;
   Int_t           PfMu3NObjsTot;
   Int_t           PfMu3NObjs;
   Double_t        PfMu3Px[20];   //[PfMu3NObjs]
   Double_t        PfMu3Py[20];   //[PfMu3NObjs]
   Double_t        PfMu3Pz[20];   //[PfMu3NObjs]
   Double_t        PfMu3Pt[20];   //[PfMu3NObjs]
   Double_t        PfMu3E[20];   //[PfMu3NObjs]
   Double_t        PfMu3Et[20];   //[PfMu3NObjs]
   Double_t        PfMu3Eta[20];   //[PfMu3NObjs]
   Double_t        PfMu3Phi[20];   //[PfMu3NObjs]
   Int_t           PfMu3Charge[20];   //[PfMu3NObjs]
   Double_t        PfMu3ParticleIso[20];   //[PfMu3NObjs]
   Double_t        PfMu3ChargedHadronIso[20];   //[PfMu3NObjs]
   Double_t        PfMu3NeutralHadronIso[20];   //[PfMu3NObjs]
   Double_t        PfMu3PhotonIso[20];   //[PfMu3NObjs]
   Double_t        PfMu3PtErr[20];   //[PfMu3NObjs]
   Int_t           PfMu3NMatches[20];   //[PfMu3NObjs]
   Int_t           PfElAntiIsoMaxLepExc;
   Int_t           PfElAntiIsoNObjsTot;
   Int_t           PfElAntiIsoNObjs;
   Double_t        PfElAntiIsoPx[20];   //[PfElAntiIsoNObjs]
   Double_t        PfElAntiIsoPy[20];   //[PfElAntiIsoNObjs]
   Double_t        PfElAntiIsoPz[20];   //[PfElAntiIsoNObjs]
   Double_t        PfElAntiIsoPt[20];   //[PfElAntiIsoNObjs]
   Double_t        PfElAntiIsoE[20];   //[PfElAntiIsoNObjs]
   Double_t        PfElAntiIsoEt[20];   //[PfElAntiIsoNObjs]
   Double_t        PfElAntiIsoEta[20];   //[PfElAntiIsoNObjs]
   Double_t        PfElAntiIsoPhi[20];   //[PfElAntiIsoNObjs]
   Int_t           PfElAntiIsoCharge[20];   //[PfElAntiIsoNObjs]
   Double_t        PfElAntiIsoParticleIso[20];   //[PfElAntiIsoNObjs]
   Double_t        PfElAntiIsoChargedHadronIso[20];   //[PfElAntiIsoNObjs]
   Double_t        PfElAntiIsoNeutralHadronIso[20];   //[PfElAntiIsoNObjs]
   Double_t        PfElAntiIsoPhotonIso[20];   //[PfElAntiIsoNObjs]
   Int_t           PfElAntiIsoID95[20];   //[PfElAntiIsoNObjs]
   Int_t           PfElAntiIsoID90[20];   //[PfElAntiIsoNObjs]
   Int_t           PfElAntiIsoID85[20];   //[PfElAntiIsoNObjs]
   Int_t           PfElAntiIsoID80[20];   //[PfElAntiIsoNObjs]
   Int_t           PfEl2MaxLepExc;
   Int_t           PfEl2NObjsTot;
   Int_t           PfEl2NObjs;
   Double_t        PfEl2Px[20];   //[PfEl2NObjs]
   Double_t        PfEl2Py[20];   //[PfEl2NObjs]
   Double_t        PfEl2Pz[20];   //[PfEl2NObjs]
   Double_t        PfEl2Pt[20];   //[PfEl2NObjs]
   Double_t        PfEl2E[20];   //[PfEl2NObjs]
   Double_t        PfEl2Et[20];   //[PfEl2NObjs]
   Double_t        PfEl2Eta[20];   //[PfEl2NObjs]
   Double_t        PfEl2Phi[20];   //[PfEl2NObjs]
   Int_t           PfEl2Charge[20];   //[PfEl2NObjs]
   Double_t        PfEl2ParticleIso[20];   //[PfEl2NObjs]
   Double_t        PfEl2ChargedHadronIso[20];   //[PfEl2NObjs]
   Double_t        PfEl2NeutralHadronIso[20];   //[PfEl2NObjs]
   Double_t        PfEl2PhotonIso[20];   //[PfEl2NObjs]
   Int_t           PfEl2ID95[20];   //[PfEl2NObjs]
   Int_t           PfEl2ID90[20];   //[PfEl2NObjs]
   Int_t           PfEl2ID85[20];   //[PfEl2NObjs]
   Int_t           PfEl2ID80[20];   //[PfEl2NObjs]
   Int_t           PfEl3MaxLepExc;
   Int_t           PfEl3NObjsTot;
   Int_t           PfEl3NObjs;
   Double_t        PfEl3Px[20];   //[PfEl3NObjs]
   Double_t        PfEl3Py[20];   //[PfEl3NObjs]
   Double_t        PfEl3Pz[20];   //[PfEl3NObjs]
   Double_t        PfEl3Pt[20];   //[PfEl3NObjs]
   Double_t        PfEl3E[20];   //[PfEl3NObjs]
   Double_t        PfEl3Et[20];   //[PfEl3NObjs]
   Double_t        PfEl3Eta[20];   //[PfEl3NObjs]
   Double_t        PfEl3Phi[20];   //[PfEl3NObjs]
   Int_t           PfEl3Charge[20];   //[PfEl3NObjs]
   Double_t        PfEl3ParticleIso[20];   //[PfEl3NObjs]
   Double_t        PfEl3ChargedHadronIso[20];   //[PfEl3NObjs]
   Double_t        PfEl3NeutralHadronIso[20];   //[PfEl3NObjs]
   Double_t        PfEl3PhotonIso[20];   //[PfEl3NObjs]
   Int_t           PfEl3ID95[20];   //[PfEl3NObjs]
   Int_t           PfEl3ID90[20];   //[PfEl3NObjs]
   Int_t           PfEl3ID85[20];   //[PfEl3NObjs]
   Int_t           PfEl3ID80[20];   //[PfEl3NObjs]
   Int_t           PfTauAntiIsoMaxLepExc;
   Int_t           PfTauAntiIsoNObjsTot;
   Int_t           PfTauAntiIsoNObjs;
   Double_t        PfTauAntiIsoPx[20];   //[PfTauAntiIsoNObjs]
   Double_t        PfTauAntiIsoPy[20];   //[PfTauAntiIsoNObjs]
   Double_t        PfTauAntiIsoPz[20];   //[PfTauAntiIsoNObjs]
   Double_t        PfTauAntiIsoPt[20];   //[PfTauAntiIsoNObjs]
   Double_t        PfTauAntiIsoE[20];   //[PfTauAntiIsoNObjs]
   Double_t        PfTauAntiIsoEt[20];   //[PfTauAntiIsoNObjs]
   Double_t        PfTauAntiIsoEta[20];   //[PfTauAntiIsoNObjs]
   Double_t        PfTauAntiIsoPhi[20];   //[PfTauAntiIsoNObjs]
   Int_t           PfTauAntiIsoCharge[20];   //[PfTauAntiIsoNObjs]
   Double_t        PfTauAntiIsoParticleIso[20];   //[PfTauAntiIsoNObjs]
   Double_t        PfTauAntiIsoChargedHadronIso[20];   //[PfTauAntiIsoNObjs]
   Double_t        PfTauAntiIsoNeutralHadronIso[20];   //[PfTauAntiIsoNObjs]
   Double_t        PfTauAntiIsoPhotonIso[20];   //[PfTauAntiIsoNObjs]
   Int_t           PfTauAntiIsoDecayMode[20];   //[PfTauAntiIsoNObjs]
   Int_t           PfTau2MaxLepExc;
   Int_t           PfTau2NObjsTot;
   Int_t           PfTau2NObjs;
   Double_t        PfTau2Px[20];   //[PfTau2NObjs]
   Double_t        PfTau2Py[20];   //[PfTau2NObjs]
   Double_t        PfTau2Pz[20];   //[PfTau2NObjs]
   Double_t        PfTau2Pt[20];   //[PfTau2NObjs]
   Double_t        PfTau2E[20];   //[PfTau2NObjs]
   Double_t        PfTau2Et[20];   //[PfTau2NObjs]
   Double_t        PfTau2Eta[20];   //[PfTau2NObjs]
   Double_t        PfTau2Phi[20];   //[PfTau2NObjs]
   Int_t           PfTau2Charge[20];   //[PfTau2NObjs]
   Double_t        PfTau2ParticleIso[20];   //[PfTau2NObjs]
   Double_t        PfTau2ChargedHadronIso[20];   //[PfTau2NObjs]
   Double_t        PfTau2NeutralHadronIso[20];   //[PfTau2NObjs]
   Double_t        PfTau2PhotonIso[20];   //[PfTau2NObjs]
   Int_t           PfTau2DecayMode[20];   //[PfTau2NObjs]
   Int_t           PfTau3MaxLepExc;
   Int_t           PfTau3NObjsTot;
   Int_t           PfTau3NObjs;
   Double_t        PfTau3Px[20];   //[PfTau3NObjs]
   Double_t        PfTau3Py[20];   //[PfTau3NObjs]
   Double_t        PfTau3Pz[20];   //[PfTau3NObjs]
   Double_t        PfTau3Pt[20];   //[PfTau3NObjs]
   Double_t        PfTau3E[20];   //[PfTau3NObjs]
   Double_t        PfTau3Et[20];   //[PfTau3NObjs]
   Double_t        PfTau3Eta[20];   //[PfTau3NObjs]
   Double_t        PfTau3Phi[20];   //[PfTau3NObjs]
   Int_t           PfTau3Charge[20];   //[PfTau3NObjs]
   Double_t        PfTau3ParticleIso[20];   //[PfTau3NObjs]
   Double_t        PfTau3ChargedHadronIso[20];   //[PfTau3NObjs]
   Double_t        PfTau3NeutralHadronIso[20];   //[PfTau3NObjs]
   Double_t        PfTau3PhotonIso[20];   //[PfTau3NObjs]
   Int_t           PfTau3DecayMode[20];   //[PfTau3NObjs]
   Int_t           NTracks;
   Int_t           NTracksTot;
   Int_t           TrkGood[500];   //[NTracks]
   Float_t         TrkPt[500];   //[NTracks]
   Float_t         TrkEta[500];   //[NTracks]
   Float_t         TrkPhi[500];   //[NTracks]
   Float_t         TrkNChi2[500];   //[NTracks]
   Float_t         TrkNHits[500];   //[NTracks]
   Float_t         TrkVtxDz[500];   //[NTracks]
   Float_t         TrkVtxDxy[500];   //[NTracks]
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
   Float_t         PFSumEt;
   Float_t         PFMETPAT;
   Float_t         PFMETPATpx;
   Float_t         PFMETPATpy;
   Float_t         PFMETPATphi;
   Float_t         PFMETPATSignificance;
   Float_t         METR12;
   Float_t         METR21;

   // List of branches
   TBranch        *b_Run;   //!
   TBranch        *b_Event;   //!
   TBranch        *b_LumiSection;   //!
   TBranch        *b_PtHat;   //!
   TBranch        *b_QCDPartonicHT;   //!
   TBranch        *b_SigProcID;   //!
   TBranch        *b_PDFScalePDF;   //!
   TBranch        *b_PDFID1;   //!
   TBranch        *b_PDFID2;   //!
   TBranch        *b_PDFx1;   //!
   TBranch        *b_PDFx2;   //!
   TBranch        *b_PDFxPDF1;   //!
   TBranch        *b_PDFxPDF2;   //!
   TBranch        *b_GenWeight;   //!
   TBranch        *b_ExtXSecLO;   //!
   TBranch        *b_IntXSec;   //!
   TBranch        *b_pdfW;   //!
   TBranch        *b_pdfWsum;   //!
   TBranch        *b_NPdfs;   //!
   TBranch        *b_PUnumInteractions;   //!
   TBranch        *b_PUnumFilled;   //!
   TBranch        *b_PUOOTnumInteractionsEarly;   //!
   TBranch        *b_PUOOTnumInteractionsLate;   //!
   TBranch        *b_PUzPositions;   //!
   TBranch        *b_PUsumPtLowPt;   //!
   TBranch        *b_PUsumPtHighPt;   //!
   TBranch        *b_PUnTrksLowPt;   //!
   TBranch        *b_PUnTrksHighPt;   //!
   TBranch        *b_Rho;   //!
   TBranch        *b_RhoPFnoPU;   //!
   TBranch        *b_Weight;   //!
   TBranch        *b_HLTResults;   //!
   TBranch        *b_HLTPrescale;   //!
   TBranch        *b_L1PhysResults;   //!
   TBranch        *b_L1TechResults;   //!
   TBranch        *b_NHLTObjs;   //!
   TBranch        *b_HLTObjectID;   //!
   TBranch        *b_HLTObjectPt;   //!
   TBranch        *b_HLTObjectEta;   //!
   TBranch        *b_HLTObjectPhi;   //!
   TBranch        *b_PUWeightTotal;   //!
   TBranch        *b_PUWeightInTime;   //!
   TBranch        *b_MassGlu;   //!
   TBranch        *b_MassChi;   //!
   TBranch        *b_MassLSP;   //!
   TBranch        *b_M0;   //!
   TBranch        *b_M12;   //!
   TBranch        *b_signMu;   //!
   TBranch        *b_A0;   //!
   TBranch        *b_process;   //!
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
   TBranch        *b_MaxGenPhoExceed;   //!
   TBranch        *b_MaxGenJetExceed;   //!
   TBranch        *b_MaxVerticesExceed;   //!
   TBranch        *b_HBHENoiseFlag;   //!
   TBranch        *b_HBHENoiseFlagIso;   //!
   TBranch        *b_CSCTightHaloID;   //!
   TBranch        *b_EcalDeadTPFilterFlag;   //!
   TBranch        *b_RecovRecHitFilterFlag;   //!
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
   TBranch        *b_NGenPhotons;   //!
   TBranch        *b_GenPhotonPt;   //!
   TBranch        *b_GenPhotonEta;   //!
   TBranch        *b_GenPhotonPhi;   //!
   TBranch        *b_GenPhotonMotherID;   //!
   TBranch        *b_GenPhotonMotherStatus;   //!
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
   TBranch        *b_NEBhits;   //!
   TBranch        *b_EBrechitE;   //!
   TBranch        *b_EBrechitPt;   //!
   TBranch        *b_EBrechitEta;   //!
   TBranch        *b_EBrechitPhi;   //!
   TBranch        *b_EBrechitChi2;   //!
   TBranch        *b_EBrechitTime;   //!
   TBranch        *b_EBrechitE4oE1;   //!
   TBranch        *b_EBrechitE2oE9;   //!
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
   TBranch        *b_ElSCindex;   //!
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
   TBranch        *b_PhoSigmaEtaEta;   //!
   TBranch        *b_PhoE1x5;   //!
   TBranch        *b_PhoE2x5;   //!
   TBranch        *b_PhoE3x3;   //!
   TBranch        *b_PhoE5x5;   //!
   TBranch        *b_PhomaxEnergyXtal;   //!
   TBranch        *b_PhoIso03HcalDepth1;   //!
   TBranch        *b_PhoIso03HcalDepth2;   //!
   TBranch        *b_PhoIso04HcalDepth1;   //!
   TBranch        *b_PhoIso04HcalDepth2;   //!
   TBranch        *b_PhoIso03nTrksSolid;   //!
   TBranch        *b_PhoIso03nTrksHollow;   //!
   TBranch        *b_PhoIso04nTrksSolid;   //!
   TBranch        *b_PhoIso04nTrksHollow;   //!
   TBranch        *b_PhoisEB;   //!
   TBranch        *b_PhoisEE;   //!
   TBranch        *b_PhoisEBEtaGap;   //!
   TBranch        *b_PhoisEBPhiGap;   //!
   TBranch        *b_PhoisEERingGap;   //!
   TBranch        *b_PhoisEEDeeGap;   //!
   TBranch        *b_PhoisEBEEGap;   //!
   TBranch        *b_PhoisPFlowPhoton;   //!
   TBranch        *b_PhoisStandardPhoton;   //!
   TBranch        *b_PhoMCmatchindex;   //!
   TBranch        *b_PhoMCmatchexitcode;   //!
   TBranch        *b_Pho_ChargedHadronIso;   //!
   TBranch        *b_Pho_NeutralHadronIso;   //!
   TBranch        *b_Pho_PhotonIso;   //!
   TBranch        *b_Pho_isPFPhoton;   //!
   TBranch        *b_Pho_isPFElectron;   //!
   TBranch        *b_PhotSCindex;   //!
   TBranch        *b_Pho_Cone04PhotonIso_dR0_dEta0_pt0;   //!
   TBranch        *b_Pho_Cone04PhotonIso_dR0_dEta0_pt5;   //!
   TBranch        *b_Pho_Cone04PhotonIso_dR8_dEta0_pt0;   //!
   TBranch        *b_Pho_Cone04PhotonIso_dR8_dEta0_pt5;   //!
   TBranch        *b_Pho_Cone01PhotonIso_dR045EB070EE_dEta015_pt08EB1EE_mvVtx;   //!
   TBranch        *b_Pho_Cone02PhotonIso_dR045EB070EE_dEta015_pt08EB1EE_mvVtx;   //!
   TBranch        *b_Pho_Cone03PhotonIso_dR045EB070EE_dEta015_pt08EB1EE_mvVtx;   //!
   TBranch        *b_Pho_Cone04PhotonIso_dR045EB070EE_dEta015_pt08EB1EE_mvVtx;   //!
   TBranch        *b_Pho_Cone04NeutralHadronIso_dR0_dEta0_pt0;   //!
   TBranch        *b_Pho_Cone04NeutralHadronIso_dR0_dEta0_pt5;   //!
   TBranch        *b_Pho_Cone04NeutralHadronIso_dR0_dEta0_pt0_nocracks;   //!
   TBranch        *b_Pho_Cone04NeutralHadronIso_dR0_dEta0_pt5_nocracks;   //!
   TBranch        *b_Pho_Cone04NeutralHadronIso_dR7_dEta0_pt0;   //!
   TBranch        *b_Pho_Cone04NeutralHadronIso_dR7_dEta0_pt5;   //!
   TBranch        *b_Pho_Cone01NeutralHadronIso_dR0_dEta0_pt0_mvVtx;   //!
   TBranch        *b_Pho_Cone02NeutralHadronIso_dR0_dEta0_pt0_mvVtx;   //!
   TBranch        *b_Pho_Cone03NeutralHadronIso_dR0_dEta0_pt0_mvVtx;   //!
   TBranch        *b_Pho_Cone04NeutralHadronIso_dR0_dEta0_pt0_mvVtx;   //!
   TBranch        *b_Pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_dz0_old;   //!
   TBranch        *b_Pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_PFnoPU_old;   //!
   TBranch        *b_Pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_dz0_old;   //!
   TBranch        *b_Pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_PFnoPU_old;   //!
   TBranch        *b_Pho_Cone01ChargedHadronIso_dR0_dEta0_pt0_dz0;   //!
   TBranch        *b_Pho_Cone01ChargedHadronIso_dR0_dEta0_pt0_dz1_dxy01;   //!
   TBranch        *b_Pho_Cone01ChargedHadronIso_dR0_dEta0_pt0_PFnoPU;   //!
   TBranch        *b_Pho_Cone01ChargedHadronIso_dR015_dEta0_pt0_dz0;   //!
   TBranch        *b_Pho_Cone01ChargedHadronIso_dR015_dEta0_pt0_dz1_dxy01;   //!
   TBranch        *b_Pho_Cone01ChargedHadronIso_dR015_dEta0_pt0_PFnoPU;   //!
   TBranch        *b_Pho_Cone02ChargedHadronIso_dR0_dEta0_pt0_dz0;   //!
   TBranch        *b_Pho_Cone02ChargedHadronIso_dR0_dEta0_pt0_dz1_dxy01;   //!
   TBranch        *b_Pho_Cone02ChargedHadronIso_dR0_dEta0_pt0_PFnoPU;   //!
   TBranch        *b_Pho_Cone02ChargedHadronIso_dR015_dEta0_pt0_dz0;   //!
   TBranch        *b_Pho_Cone02ChargedHadronIso_dR015_dEta0_pt0_dz1_dxy01;   //!
   TBranch        *b_Pho_Cone02ChargedHadronIso_dR015_dEta0_pt0_PFnoPU;   //!
   TBranch        *b_Pho_Cone03ChargedHadronIso_dR0_dEta0_pt0_dz0;   //!
   TBranch        *b_Pho_Cone03ChargedHadronIso_dR0_dEta0_pt0_dz1_dxy01;   //!
   TBranch        *b_Pho_Cone03ChargedHadronIso_dR0_dEta0_pt0_PFnoPU;   //!
   TBranch        *b_Pho_Cone03ChargedHadronIso_dR015_dEta0_pt0_dz0;   //!
   TBranch        *b_Pho_Cone03ChargedHadronIso_dR015_dEta0_pt0_dz1_dxy01;   //!
   TBranch        *b_Pho_Cone03ChargedHadronIso_dR015_dEta0_pt0_PFnoPU;   //!
   TBranch        *b_Pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_dz0;   //!
   TBranch        *b_Pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_dz1_dxy01;   //!
   TBranch        *b_Pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_PFnoPU;   //!
   TBranch        *b_Pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_dz0;   //!
   TBranch        *b_Pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_dz1_dxy01;   //!
   TBranch        *b_Pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_PFnoPU;   //!
   TBranch        *b_NSuperClusters;   //!
   TBranch        *b_SCRaw;   //!
   TBranch        *b_SCPre;   //!
   TBranch        *b_SCEnergy;   //!
   TBranch        *b_SCEta;   //!
   TBranch        *b_SCPhi;   //!
   TBranch        *b_SCPhiWidth;   //!
   TBranch        *b_SCEtaWidth;   //!
   TBranch        *b_SCBrem;   //!
   TBranch        *b_SCR9;   //!
   TBranch        *b_SCcrackcorrseed;   //!
   TBranch        *b_SCcrackcorr;   //!
   TBranch        *b_SClocalcorrseed;   //!
   TBranch        *b_SClocalcorr;   //!
   TBranch        *b_SCcrackcorrseedfactor;   //!
   TBranch        *b_SClocalcorrseedfactor;   //!
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
   TBranch        *b_JArea;   //!
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
   TBranch        *b_JGenJetIndex;   //!
   TBranch        *b_CANJets;   //!
   TBranch        *b_CAJPx;   //!
   TBranch        *b_CAJPy;   //!
   TBranch        *b_CAJPz;   //!
   TBranch        *b_CAJPt;   //!
   TBranch        *b_CAJE;   //!
   TBranch        *b_CAJEt;   //!
   TBranch        *b_CAJEta;   //!
   TBranch        *b_CAJPhi;   //!
   TBranch        *b_CAJFlavour;   //!
   TBranch        *b_CAJScale;   //!
   TBranch        *b_CAJL1FastJetScale;   //!
   TBranch        *b_CAJArea;   //!
   TBranch        *b_CAJbTagProbTkCntHighEff;   //!
   TBranch        *b_CAJbTagProbTkCntHighPur;   //!
   TBranch        *b_CAJbTagProbSimpSVHighEff;   //!
   TBranch        *b_CAJbTagProbSimpSVHighPur;   //!
   TBranch        *b_CAJIDLoose;   //!
   TBranch        *b_CAJID_HPD;   //!
   TBranch        *b_CAJID_RBX;   //!
   TBranch        *b_CAJID_n90Hits;   //!
   TBranch        *b_CAJID_resEMF;   //!
   TBranch        *b_CAJEMfrac;   //!
   TBranch        *b_CAJNAssoTracks;   //!
   TBranch        *b_CAJChfrac;   //!
   TBranch        *b_CAJNConstituents;   //!
   TBranch        *b_CAJn90;   //!
   TBranch        *b_PF2PATAntiIsoNJets;   //!
   TBranch        *b_PF2PATAntiIsoJPx;   //!
   TBranch        *b_PF2PATAntiIsoJPy;   //!
   TBranch        *b_PF2PATAntiIsoJPz;   //!
   TBranch        *b_PF2PATAntiIsoJPt;   //!
   TBranch        *b_PF2PATAntiIsoJE;   //!
   TBranch        *b_PF2PATAntiIsoJEt;   //!
   TBranch        *b_PF2PATAntiIsoJEta;   //!
   TBranch        *b_PF2PATAntiIsoJPhi;   //!
   TBranch        *b_PF2PATAntiIsoJFlavour;   //!
   TBranch        *b_PF2PATAntiIsoJScale;   //!
   TBranch        *b_PF2PATAntiIsoJL1FastJetScale;   //!
   TBranch        *b_PF2PATAntiIsoJArea;   //!
   TBranch        *b_PF2PATAntiIsoJbTagProbTkCntHighEff;   //!
   TBranch        *b_PF2PATAntiIsoJbTagProbTkCntHighPur;   //!
   TBranch        *b_PF2PATAntiIsoJbTagProbSimpSVHighEff;   //!
   TBranch        *b_PF2PATAntiIsoJbTagProbSimpSVHighPur;   //!
   TBranch        *b_PF2PATAntiIsoJIDLoose;   //!
   TBranch        *b_PF2PATAntiIsoJChMult;   //!
   TBranch        *b_PF2PATAntiIsoJNeuMult;   //!
   TBranch        *b_PF2PATAntiIsoJChHadfrac;   //!
   TBranch        *b_PF2PATAntiIsoJNeuHadfrac;   //!
   TBranch        *b_PF2PATAntiIsoJChEmfrac;   //!
   TBranch        *b_PF2PATAntiIsoJNeuEmfrac;   //!
   TBranch        *b_PF2PATAntiIsoJChMufrac;   //!
   TBranch        *b_PF2PATAntiIsoJNConstituents;   //!
   TBranch        *b_PF2PAT2NJets;   //!
   TBranch        *b_PF2PAT2JPx;   //!
   TBranch        *b_PF2PAT2JPy;   //!
   TBranch        *b_PF2PAT2JPz;   //!
   TBranch        *b_PF2PAT2JPt;   //!
   TBranch        *b_PF2PAT2JE;   //!
   TBranch        *b_PF2PAT2JEt;   //!
   TBranch        *b_PF2PAT2JEta;   //!
   TBranch        *b_PF2PAT2JPhi;   //!
   TBranch        *b_PF2PAT2JFlavour;   //!
   TBranch        *b_PF2PAT2JScale;   //!
   TBranch        *b_PF2PAT2JL1FastJetScale;   //!
   TBranch        *b_PF2PAT2JArea;   //!
   TBranch        *b_PF2PAT2JbTagProbTkCntHighEff;   //!
   TBranch        *b_PF2PAT2JbTagProbTkCntHighPur;   //!
   TBranch        *b_PF2PAT2JbTagProbSimpSVHighEff;   //!
   TBranch        *b_PF2PAT2JbTagProbSimpSVHighPur;   //!
   TBranch        *b_PF2PAT2JIDLoose;   //!
   TBranch        *b_PF2PAT2JChMult;   //!
   TBranch        *b_PF2PAT2JNeuMult;   //!
   TBranch        *b_PF2PAT2JChHadfrac;   //!
   TBranch        *b_PF2PAT2JNeuHadfrac;   //!
   TBranch        *b_PF2PAT2JChEmfrac;   //!
   TBranch        *b_PF2PAT2JNeuEmfrac;   //!
   TBranch        *b_PF2PAT2JChMufrac;   //!
   TBranch        *b_PF2PAT2JNConstituents;   //!
   TBranch        *b_PF2PAT3NJets;   //!
   TBranch        *b_PF2PAT3JPx;   //!
   TBranch        *b_PF2PAT3JPy;   //!
   TBranch        *b_PF2PAT3JPz;   //!
   TBranch        *b_PF2PAT3JPt;   //!
   TBranch        *b_PF2PAT3JE;   //!
   TBranch        *b_PF2PAT3JEt;   //!
   TBranch        *b_PF2PAT3JEta;   //!
   TBranch        *b_PF2PAT3JPhi;   //!
   TBranch        *b_PF2PAT3JFlavour;   //!
   TBranch        *b_PF2PAT3JScale;   //!
   TBranch        *b_PF2PAT3JL1FastJetScale;   //!
   TBranch        *b_PF2PAT3JArea;   //!
   TBranch        *b_PF2PAT3JbTagProbTkCntHighEff;   //!
   TBranch        *b_PF2PAT3JbTagProbTkCntHighPur;   //!
   TBranch        *b_PF2PAT3JbTagProbSimpSVHighEff;   //!
   TBranch        *b_PF2PAT3JbTagProbSimpSVHighPur;   //!
   TBranch        *b_PF2PAT3JIDLoose;   //!
   TBranch        *b_PF2PAT3JChMult;   //!
   TBranch        *b_PF2PAT3JNeuMult;   //!
   TBranch        *b_PF2PAT3JChHadfrac;   //!
   TBranch        *b_PF2PAT3JNeuHadfrac;   //!
   TBranch        *b_PF2PAT3JChEmfrac;   //!
   TBranch        *b_PF2PAT3JNeuEmfrac;   //!
   TBranch        *b_PF2PAT3JChMufrac;   //!
   TBranch        *b_PF2PAT3JNConstituents;   //!
   TBranch        *b_PfMuAntiIsoMaxLepExc;   //!
   TBranch        *b_PfMuAntiIsoNObjsTot;   //!
   TBranch        *b_PfMuAntiIsoNObjs;   //!
   TBranch        *b_PfMuAntiIsoPx;   //!
   TBranch        *b_PfMuAntiIsoPy;   //!
   TBranch        *b_PfMuAntiIsoPz;   //!
   TBranch        *b_PfMuAntiIsoPt;   //!
   TBranch        *b_PfMuAntiIsoE;   //!
   TBranch        *b_PfMuAntiIsoEt;   //!
   TBranch        *b_PfMuAntiIsoEta;   //!
   TBranch        *b_PfMuAntiIsoPhi;   //!
   TBranch        *b_PfMuAntiIsoCharge;   //!
   TBranch        *b_PfMuAntiIsoParticleIso;   //!
   TBranch        *b_PfMuAntiIsoChargedHadronIso;   //!
   TBranch        *b_PfMuAntiIsoNeutralHadronIso;   //!
   TBranch        *b_PfMuAntiIsoPhotonIso;   //!
   TBranch        *b_PfMuAntiIsoPtErr;   //!
   TBranch        *b_PfMuAntiIsoNMatches;   //!
   TBranch        *b_PfMu2MaxLepExc;   //!
   TBranch        *b_PfMu2NObjsTot;   //!
   TBranch        *b_PfMu2NObjs;   //!
   TBranch        *b_PfMu2Px;   //!
   TBranch        *b_PfMu2Py;   //!
   TBranch        *b_PfMu2Pz;   //!
   TBranch        *b_PfMu2Pt;   //!
   TBranch        *b_PfMu2E;   //!
   TBranch        *b_PfMu2Et;   //!
   TBranch        *b_PfMu2Eta;   //!
   TBranch        *b_PfMu2Phi;   //!
   TBranch        *b_PfMu2Charge;   //!
   TBranch        *b_PfMu2ParticleIso;   //!
   TBranch        *b_PfMu2ChargedHadronIso;   //!
   TBranch        *b_PfMu2NeutralHadronIso;   //!
   TBranch        *b_PfMu2PhotonIso;   //!
   TBranch        *b_PfMu2PtErr;   //!
   TBranch        *b_PfMu2NMatches;   //!
   TBranch        *b_PfMu3MaxLepExc;   //!
   TBranch        *b_PfMu3NObjsTot;   //!
   TBranch        *b_PfMu3NObjs;   //!
   TBranch        *b_PfMu3Px;   //!
   TBranch        *b_PfMu3Py;   //!
   TBranch        *b_PfMu3Pz;   //!
   TBranch        *b_PfMu3Pt;   //!
   TBranch        *b_PfMu3E;   //!
   TBranch        *b_PfMu3Et;   //!
   TBranch        *b_PfMu3Eta;   //!
   TBranch        *b_PfMu3Phi;   //!
   TBranch        *b_PfMu3Charge;   //!
   TBranch        *b_PfMu3ParticleIso;   //!
   TBranch        *b_PfMu3ChargedHadronIso;   //!
   TBranch        *b_PfMu3NeutralHadronIso;   //!
   TBranch        *b_PfMu3PhotonIso;   //!
   TBranch        *b_PfMu3PtErr;   //!
   TBranch        *b_PfMu3NMatches;   //!
   TBranch        *b_PfElAntiIsoMaxLepExc;   //!
   TBranch        *b_PfElAntiIsoNObjsTot;   //!
   TBranch        *b_PfElAntiIsoNObjs;   //!
   TBranch        *b_PfElAntiIsoPx;   //!
   TBranch        *b_PfElAntiIsoPy;   //!
   TBranch        *b_PfElAntiIsoPz;   //!
   TBranch        *b_PfElAntiIsoPt;   //!
   TBranch        *b_PfElAntiIsoE;   //!
   TBranch        *b_PfElAntiIsoEt;   //!
   TBranch        *b_PfElAntiIsoEta;   //!
   TBranch        *b_PfElAntiIsoPhi;   //!
   TBranch        *b_PfElAntiIsoCharge;   //!
   TBranch        *b_PfElAntiIsoParticleIso;   //!
   TBranch        *b_PfElAntiIsoChargedHadronIso;   //!
   TBranch        *b_PfElAntiIsoNeutralHadronIso;   //!
   TBranch        *b_PfElAntiIsoPhotonIso;   //!
   TBranch        *b_PfElAntiIsoID95;   //!
   TBranch        *b_PfElAntiIsoID90;   //!
   TBranch        *b_PfElAntiIsoID85;   //!
   TBranch        *b_PfElAntiIsoID80;   //!
   TBranch        *b_PfEl2MaxLepExc;   //!
   TBranch        *b_PfEl2NObjsTot;   //!
   TBranch        *b_PfEl2NObjs;   //!
   TBranch        *b_PfEl2Px;   //!
   TBranch        *b_PfEl2Py;   //!
   TBranch        *b_PfEl2Pz;   //!
   TBranch        *b_PfEl2Pt;   //!
   TBranch        *b_PfEl2E;   //!
   TBranch        *b_PfEl2Et;   //!
   TBranch        *b_PfEl2Eta;   //!
   TBranch        *b_PfEl2Phi;   //!
   TBranch        *b_PfEl2Charge;   //!
   TBranch        *b_PfEl2ParticleIso;   //!
   TBranch        *b_PfEl2ChargedHadronIso;   //!
   TBranch        *b_PfEl2NeutralHadronIso;   //!
   TBranch        *b_PfEl2PhotonIso;   //!
   TBranch        *b_PfEl2ID95;   //!
   TBranch        *b_PfEl2ID90;   //!
   TBranch        *b_PfEl2ID85;   //!
   TBranch        *b_PfEl2ID80;   //!
   TBranch        *b_PfEl3MaxLepExc;   //!
   TBranch        *b_PfEl3NObjsTot;   //!
   TBranch        *b_PfEl3NObjs;   //!
   TBranch        *b_PfEl3Px;   //!
   TBranch        *b_PfEl3Py;   //!
   TBranch        *b_PfEl3Pz;   //!
   TBranch        *b_PfEl3Pt;   //!
   TBranch        *b_PfEl3E;   //!
   TBranch        *b_PfEl3Et;   //!
   TBranch        *b_PfEl3Eta;   //!
   TBranch        *b_PfEl3Phi;   //!
   TBranch        *b_PfEl3Charge;   //!
   TBranch        *b_PfEl3ParticleIso;   //!
   TBranch        *b_PfEl3ChargedHadronIso;   //!
   TBranch        *b_PfEl3NeutralHadronIso;   //!
   TBranch        *b_PfEl3PhotonIso;   //!
   TBranch        *b_PfEl3ID95;   //!
   TBranch        *b_PfEl3ID90;   //!
   TBranch        *b_PfEl3ID85;   //!
   TBranch        *b_PfEl3ID80;   //!
   TBranch        *b_PfTauAntiIsoMaxLepExc;   //!
   TBranch        *b_PfTauAntiIsoNObjsTot;   //!
   TBranch        *b_PfTauAntiIsoNObjs;   //!
   TBranch        *b_PfTauAntiIsoPx;   //!
   TBranch        *b_PfTauAntiIsoPy;   //!
   TBranch        *b_PfTauAntiIsoPz;   //!
   TBranch        *b_PfTauAntiIsoPt;   //!
   TBranch        *b_PfTauAntiIsoE;   //!
   TBranch        *b_PfTauAntiIsoEt;   //!
   TBranch        *b_PfTauAntiIsoEta;   //!
   TBranch        *b_PfTauAntiIsoPhi;   //!
   TBranch        *b_PfTauAntiIsoCharge;   //!
   TBranch        *b_PfTauAntiIsoParticleIso;   //!
   TBranch        *b_PfTauAntiIsoChargedHadronIso;   //!
   TBranch        *b_PfTauAntiIsoNeutralHadronIso;   //!
   TBranch        *b_PfTauAntiIsoPhotonIso;   //!
   TBranch        *b_PfTauAntiIsoDecayMode;   //!
   TBranch        *b_PfTau2MaxLepExc;   //!
   TBranch        *b_PfTau2NObjsTot;   //!
   TBranch        *b_PfTau2NObjs;   //!
   TBranch        *b_PfTau2Px;   //!
   TBranch        *b_PfTau2Py;   //!
   TBranch        *b_PfTau2Pz;   //!
   TBranch        *b_PfTau2Pt;   //!
   TBranch        *b_PfTau2E;   //!
   TBranch        *b_PfTau2Et;   //!
   TBranch        *b_PfTau2Eta;   //!
   TBranch        *b_PfTau2Phi;   //!
   TBranch        *b_PfTau2Charge;   //!
   TBranch        *b_PfTau2ParticleIso;   //!
   TBranch        *b_PfTau2ChargedHadronIso;   //!
   TBranch        *b_PfTau2NeutralHadronIso;   //!
   TBranch        *b_PfTau2PhotonIso;   //!
   TBranch        *b_PfTau2DecayMode;   //!
   TBranch        *b_PfTau3MaxLepExc;   //!
   TBranch        *b_PfTau3NObjsTot;   //!
   TBranch        *b_PfTau3NObjs;   //!
   TBranch        *b_PfTau3Px;   //!
   TBranch        *b_PfTau3Py;   //!
   TBranch        *b_PfTau3Pz;   //!
   TBranch        *b_PfTau3Pt;   //!
   TBranch        *b_PfTau3E;   //!
   TBranch        *b_PfTau3Et;   //!
   TBranch        *b_PfTau3Eta;   //!
   TBranch        *b_PfTau3Phi;   //!
   TBranch        *b_PfTau3Charge;   //!
   TBranch        *b_PfTau3ParticleIso;   //!
   TBranch        *b_PfTau3ChargedHadronIso;   //!
   TBranch        *b_PfTau3NeutralHadronIso;   //!
   TBranch        *b_PfTau3PhotonIso;   //!
   TBranch        *b_PfTau3DecayMode;   //!
   TBranch        *b_NTracks;   //!
   TBranch        *b_NTracksTot;   //!
   TBranch        *b_TrkGood;   //!
   TBranch        *b_TrkPt;   //!
   TBranch        *b_TrkEta;   //!
   TBranch        *b_TrkPhi;   //!
   TBranch        *b_TrkNChi2;   //!
   TBranch        *b_TrkNHits;   //!
   TBranch        *b_TrkVtxDz;   //!
   TBranch        *b_TrkVtxDxy;   //!
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
   TBranch        *b_PFSumEt;   //!
   TBranch        *b_PFMETPAT;   //!
   TBranch        *b_PFMETPATpx;   //!
   TBranch        *b_PFMETPATpy;   //!
   TBranch        *b_PFMETPATphi;   //!
   TBranch        *b_PFMETPATSignificance;   //!
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
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("dcap://t3se01.psi.ch:22125/pnfs/psi.ch/cms/trivcat/store/user/susy/ntuples/mc/V02-05-02/G_Pt-300to470_Tune23_7TeV_herwigpp_Summer11_PU_S3/NTupleProducer_42X_MC_259_1_Le7.root");
      if (!f) {
         f = new TFile("dcap://t3se01.psi.ch:22125/pnfs/psi.ch/cms/trivcat/store/user/susy/ntuples/mc/V02-05-02/G_Pt-300to470_Tune23_7TeV_herwigpp_Summer11_PU_S3/NTupleProducer_42X_MC_259_1_Le7.root");
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
   fChain->SetBranchAddress("QCDPartonicHT", &QCDPartonicHT, &b_QCDPartonicHT);
   fChain->SetBranchAddress("SigProcID", &SigProcID, &b_SigProcID);
   fChain->SetBranchAddress("PDFScalePDF", &PDFScalePDF, &b_PDFScalePDF);
   fChain->SetBranchAddress("PDFID1", &PDFID1, &b_PDFID1);
   fChain->SetBranchAddress("PDFID2", &PDFID2, &b_PDFID2);
   fChain->SetBranchAddress("PDFx1", &PDFx1, &b_PDFx1);
   fChain->SetBranchAddress("PDFx2", &PDFx2, &b_PDFx2);
   fChain->SetBranchAddress("PDFxPDF1", &PDFxPDF1, &b_PDFxPDF1);
   fChain->SetBranchAddress("PDFxPDF2", &PDFxPDF2, &b_PDFxPDF2);
   fChain->SetBranchAddress("GenWeight", &GenWeight, &b_GenWeight);
   fChain->SetBranchAddress("ExtXSecLO", &ExtXSecLO, &b_ExtXSecLO);
   fChain->SetBranchAddress("IntXSec", &IntXSec, &b_IntXSec);
   fChain->SetBranchAddress("pdfW", pdfW, &b_pdfW);
   fChain->SetBranchAddress("pdfWsum", &pdfWsum, &b_pdfWsum);
   fChain->SetBranchAddress("NPdfs", &NPdfs, &b_NPdfs);
   fChain->SetBranchAddress("PUnumInteractions", &PUnumInteractions, &b_PUnumInteractions);
   fChain->SetBranchAddress("PUnumFilled", &PUnumFilled, &b_PUnumFilled);
   fChain->SetBranchAddress("PUOOTnumInteractionsEarly", &PUOOTnumInteractionsEarly, &b_PUOOTnumInteractionsEarly);
   fChain->SetBranchAddress("PUOOTnumInteractionsLate", &PUOOTnumInteractionsLate, &b_PUOOTnumInteractionsLate);
   fChain->SetBranchAddress("PUzPositions", PUzPositions, &b_PUzPositions);
   fChain->SetBranchAddress("PUsumPtLowPt", PUsumPtLowPt, &b_PUsumPtLowPt);
   fChain->SetBranchAddress("PUsumPtHighPt", PUsumPtHighPt, &b_PUsumPtHighPt);
   fChain->SetBranchAddress("PUnTrksLowPt", PUnTrksLowPt, &b_PUnTrksLowPt);
   fChain->SetBranchAddress("PUnTrksHighPt", PUnTrksHighPt, &b_PUnTrksHighPt);
   fChain->SetBranchAddress("Rho", &Rho, &b_Rho);
   fChain->SetBranchAddress("RhoPFnoPU", &RhoPFnoPU, &b_RhoPFnoPU);
   fChain->SetBranchAddress("Weight", &Weight, &b_Weight);
   fChain->SetBranchAddress("HLTResults", HLTResults, &b_HLTResults);
   fChain->SetBranchAddress("HLTPrescale", HLTPrescale, &b_HLTPrescale);
   fChain->SetBranchAddress("L1PhysResults", L1PhysResults, &b_L1PhysResults);
   fChain->SetBranchAddress("L1TechResults", L1TechResults, &b_L1TechResults);
   fChain->SetBranchAddress("NHLTObjs", &NHLTObjs, &b_NHLTObjs);
   fChain->SetBranchAddress("HLTObjectID", HLTObjectID, &b_HLTObjectID);
   fChain->SetBranchAddress("HLTObjectPt", HLTObjectPt, &b_HLTObjectPt);
   fChain->SetBranchAddress("HLTObjectEta", HLTObjectEta, &b_HLTObjectEta);
   fChain->SetBranchAddress("HLTObjectPhi", HLTObjectPhi, &b_HLTObjectPhi);
   fChain->SetBranchAddress("PUWeightTotal", &PUWeightTotal, &b_PUWeightTotal);
   fChain->SetBranchAddress("PUWeightInTime", &PUWeightInTime, &b_PUWeightInTime);
   fChain->SetBranchAddress("MassGlu", &MassGlu, &b_MassGlu);
   fChain->SetBranchAddress("MassChi", &MassChi, &b_MassChi);
   fChain->SetBranchAddress("MassLSP", &MassLSP, &b_MassLSP);
   fChain->SetBranchAddress("M0", &M0, &b_M0);
   fChain->SetBranchAddress("M12", &M12, &b_M12);
   fChain->SetBranchAddress("signMu", &signMu, &b_signMu);
   fChain->SetBranchAddress("A0", &A0, &b_A0);
   fChain->SetBranchAddress("process", &process, &b_process);
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
   fChain->SetBranchAddress("MaxGenPhoExceed", &MaxGenPhoExceed, &b_MaxGenPhoExceed);
   fChain->SetBranchAddress("MaxGenJetExceed", &MaxGenJetExceed, &b_MaxGenJetExceed);
   fChain->SetBranchAddress("MaxVerticesExceed", &MaxVerticesExceed, &b_MaxVerticesExceed);
   fChain->SetBranchAddress("HBHENoiseFlag", &HBHENoiseFlag, &b_HBHENoiseFlag);
   fChain->SetBranchAddress("HBHENoiseFlagIso", &HBHENoiseFlagIso, &b_HBHENoiseFlagIso);
   fChain->SetBranchAddress("CSCTightHaloID", &CSCTightHaloID, &b_CSCTightHaloID);
   fChain->SetBranchAddress("EcalDeadTPFilterFlag", &EcalDeadTPFilterFlag, &b_EcalDeadTPFilterFlag);
   fChain->SetBranchAddress("RecovRecHitFilterFlag", &RecovRecHitFilterFlag, &b_RecovRecHitFilterFlag);
   fChain->SetBranchAddress("NGenLeptons", &NGenLeptons, &b_NGenLeptons);
   fChain->SetBranchAddress("GenLeptonID", GenLeptonID, &b_GenLeptonID);
   fChain->SetBranchAddress("GenLeptonPt", GenLeptonPt, &b_GenLeptonPt);
   fChain->SetBranchAddress("GenLeptonEta", GenLeptonEta, &b_GenLeptonEta);
   fChain->SetBranchAddress("GenLeptonPhi", GenLeptonPhi, &b_GenLeptonPhi);
   fChain->SetBranchAddress("GenLeptonMID", GenLeptonMID, &b_GenLeptonMID);
   fChain->SetBranchAddress("GenLeptonMStatus", GenLeptonMStatus, &b_GenLeptonMStatus);
   fChain->SetBranchAddress("GenLeptonMPt", GenLeptonMPt, &b_GenLeptonMPt);
   fChain->SetBranchAddress("GenLeptonMEta", GenLeptonMEta, &b_GenLeptonMEta);
   fChain->SetBranchAddress("GenLeptonMPhi", GenLeptonMPhi, &b_GenLeptonMPhi);
   fChain->SetBranchAddress("GenLeptonGMID", GenLeptonGMID, &b_GenLeptonGMID);
   fChain->SetBranchAddress("GenLeptonGMStatus", GenLeptonGMStatus, &b_GenLeptonGMStatus);
   fChain->SetBranchAddress("GenLeptonGMPt", GenLeptonGMPt, &b_GenLeptonGMPt);
   fChain->SetBranchAddress("GenLeptonGMEta", GenLeptonGMEta, &b_GenLeptonGMEta);
   fChain->SetBranchAddress("GenLeptonGMPhi", GenLeptonGMPhi, &b_GenLeptonGMPhi);
   fChain->SetBranchAddress("NGenPhotons", &NGenPhotons, &b_NGenPhotons);
   fChain->SetBranchAddress("GenPhotonPt", GenPhotonPt, &b_GenPhotonPt);
   fChain->SetBranchAddress("GenPhotonEta", GenPhotonEta, &b_GenPhotonEta);
   fChain->SetBranchAddress("GenPhotonPhi", GenPhotonPhi, &b_GenPhotonPhi);
   fChain->SetBranchAddress("GenPhotonMotherID", GenPhotonMotherID, &b_GenPhotonMotherID);
   fChain->SetBranchAddress("GenPhotonMotherStatus", GenPhotonMotherStatus, &b_GenPhotonMotherStatus);
   fChain->SetBranchAddress("NGenJets", &NGenJets, &b_NGenJets);
   fChain->SetBranchAddress("GenJetPt", GenJetPt, &b_GenJetPt);
   fChain->SetBranchAddress("GenJetEta", GenJetEta, &b_GenJetEta);
   fChain->SetBranchAddress("GenJetPhi", GenJetPhi, &b_GenJetPhi);
   fChain->SetBranchAddress("GenJetE", GenJetE, &b_GenJetE);
   fChain->SetBranchAddress("GenJetEmE", GenJetEmE, &b_GenJetEmE);
   fChain->SetBranchAddress("GenJetHadE", GenJetHadE, &b_GenJetHadE);
   fChain->SetBranchAddress("GenJetInvE", GenJetInvE, &b_GenJetInvE);
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
   fChain->SetBranchAddress("NEBhits", &NEBhits, &b_NEBhits);
   fChain->SetBranchAddress("EBrechitE", EBrechitE, &b_EBrechitE);
   fChain->SetBranchAddress("EBrechitPt", EBrechitPt, &b_EBrechitPt);
   fChain->SetBranchAddress("EBrechitEta", EBrechitEta, &b_EBrechitEta);
   fChain->SetBranchAddress("EBrechitPhi", EBrechitPhi, &b_EBrechitPhi);
   fChain->SetBranchAddress("EBrechitChi2", EBrechitChi2, &b_EBrechitChi2);
   fChain->SetBranchAddress("EBrechitTime", EBrechitTime, &b_EBrechitTime);
   fChain->SetBranchAddress("EBrechitE4oE1", EBrechitE4oE1, &b_EBrechitE4oE1);
   fChain->SetBranchAddress("EBrechitE2oE9", EBrechitE2oE9, &b_EBrechitE2oE9);
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
   fChain->SetBranchAddress("ElSCindex", ElSCindex, &b_ElSCindex);
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
   fChain->SetBranchAddress("PhoSigmaEtaEta", PhoSigmaEtaEta, &b_PhoSigmaEtaEta);
   fChain->SetBranchAddress("PhoE1x5", PhoE1x5, &b_PhoE1x5);
   fChain->SetBranchAddress("PhoE2x5", PhoE2x5, &b_PhoE2x5);
   fChain->SetBranchAddress("PhoE3x3", PhoE3x3, &b_PhoE3x3);
   fChain->SetBranchAddress("PhoE5x5", PhoE5x5, &b_PhoE5x5);
   fChain->SetBranchAddress("PhomaxEnergyXtal", PhomaxEnergyXtal, &b_PhomaxEnergyXtal);
   fChain->SetBranchAddress("PhoIso03HcalDepth1", PhoIso03HcalDepth1, &b_PhoIso03HcalDepth1);
   fChain->SetBranchAddress("PhoIso03HcalDepth2", PhoIso03HcalDepth2, &b_PhoIso03HcalDepth2);
   fChain->SetBranchAddress("PhoIso04HcalDepth1", PhoIso04HcalDepth1, &b_PhoIso04HcalDepth1);
   fChain->SetBranchAddress("PhoIso04HcalDepth2", PhoIso04HcalDepth2, &b_PhoIso04HcalDepth2);
   fChain->SetBranchAddress("PhoIso03nTrksSolid", PhoIso03nTrksSolid, &b_PhoIso03nTrksSolid);
   fChain->SetBranchAddress("PhoIso03nTrksHollow", PhoIso03nTrksHollow, &b_PhoIso03nTrksHollow);
   fChain->SetBranchAddress("PhoIso04nTrksSolid", PhoIso04nTrksSolid, &b_PhoIso04nTrksSolid);
   fChain->SetBranchAddress("PhoIso04nTrksHollow", PhoIso04nTrksHollow, &b_PhoIso04nTrksHollow);
   fChain->SetBranchAddress("PhoisEB", PhoisEB, &b_PhoisEB);
   fChain->SetBranchAddress("PhoisEE", PhoisEE, &b_PhoisEE);
   fChain->SetBranchAddress("PhoisEBEtaGap", PhoisEBEtaGap, &b_PhoisEBEtaGap);
   fChain->SetBranchAddress("PhoisEBPhiGap", PhoisEBPhiGap, &b_PhoisEBPhiGap);
   fChain->SetBranchAddress("PhoisEERingGap", PhoisEERingGap, &b_PhoisEERingGap);
   fChain->SetBranchAddress("PhoisEEDeeGap", PhoisEEDeeGap, &b_PhoisEEDeeGap);
   fChain->SetBranchAddress("PhoisEBEEGap", PhoisEBEEGap, &b_PhoisEBEEGap);
   fChain->SetBranchAddress("PhoisPFlowPhoton", PhoisPFlowPhoton, &b_PhoisPFlowPhoton);
   fChain->SetBranchAddress("PhoisStandardPhoton", PhoisStandardPhoton, &b_PhoisStandardPhoton);
   fChain->SetBranchAddress("PhoMCmatchindex", PhoMCmatchindex, &b_PhoMCmatchindex);
   fChain->SetBranchAddress("PhoMCmatchexitcode", PhoMCmatchexitcode, &b_PhoMCmatchexitcode);
   fChain->SetBranchAddress("Pho_ChargedHadronIso", Pho_ChargedHadronIso, &b_Pho_ChargedHadronIso);
   fChain->SetBranchAddress("Pho_NeutralHadronIso", Pho_NeutralHadronIso, &b_Pho_NeutralHadronIso);
   fChain->SetBranchAddress("Pho_PhotonIso", Pho_PhotonIso, &b_Pho_PhotonIso);
   fChain->SetBranchAddress("Pho_isPFPhoton", Pho_isPFPhoton, &b_Pho_isPFPhoton);
   fChain->SetBranchAddress("Pho_isPFElectron", Pho_isPFElectron, &b_Pho_isPFElectron);
   fChain->SetBranchAddress("PhotSCindex", PhotSCindex, &b_PhotSCindex);
   fChain->SetBranchAddress("Pho_Cone04PhotonIso_dR0_dEta0_pt0", Pho_Cone04PhotonIso_dR0_dEta0_pt0, &b_Pho_Cone04PhotonIso_dR0_dEta0_pt0);
   fChain->SetBranchAddress("Pho_Cone04PhotonIso_dR0_dEta0_pt5", Pho_Cone04PhotonIso_dR0_dEta0_pt5, &b_Pho_Cone04PhotonIso_dR0_dEta0_pt5);
   fChain->SetBranchAddress("Pho_Cone04PhotonIso_dR8_dEta0_pt0", Pho_Cone04PhotonIso_dR8_dEta0_pt0, &b_Pho_Cone04PhotonIso_dR8_dEta0_pt0);
   fChain->SetBranchAddress("Pho_Cone04PhotonIso_dR8_dEta0_pt5", Pho_Cone04PhotonIso_dR8_dEta0_pt5, &b_Pho_Cone04PhotonIso_dR8_dEta0_pt5);
   fChain->SetBranchAddress("Pho_Cone01PhotonIso_dR045EB070EE_dEta015_pt08EB1EE_mvVtx", Pho_Cone01PhotonIso_dR045EB070EE_dEta015_pt08EB1EE_mvVtx, &b_Pho_Cone01PhotonIso_dR045EB070EE_dEta015_pt08EB1EE_mvVtx);
   fChain->SetBranchAddress("Pho_Cone02PhotonIso_dR045EB070EE_dEta015_pt08EB1EE_mvVtx", Pho_Cone02PhotonIso_dR045EB070EE_dEta015_pt08EB1EE_mvVtx, &b_Pho_Cone02PhotonIso_dR045EB070EE_dEta015_pt08EB1EE_mvVtx);
   fChain->SetBranchAddress("Pho_Cone03PhotonIso_dR045EB070EE_dEta015_pt08EB1EE_mvVtx", Pho_Cone03PhotonIso_dR045EB070EE_dEta015_pt08EB1EE_mvVtx, &b_Pho_Cone03PhotonIso_dR045EB070EE_dEta015_pt08EB1EE_mvVtx);
   fChain->SetBranchAddress("Pho_Cone04PhotonIso_dR045EB070EE_dEta015_pt08EB1EE_mvVtx", Pho_Cone04PhotonIso_dR045EB070EE_dEta015_pt08EB1EE_mvVtx, &b_Pho_Cone04PhotonIso_dR045EB070EE_dEta015_pt08EB1EE_mvVtx);
   fChain->SetBranchAddress("Pho_Cone04NeutralHadronIso_dR0_dEta0_pt0", Pho_Cone04NeutralHadronIso_dR0_dEta0_pt0, &b_Pho_Cone04NeutralHadronIso_dR0_dEta0_pt0);
   fChain->SetBranchAddress("Pho_Cone04NeutralHadronIso_dR0_dEta0_pt5", Pho_Cone04NeutralHadronIso_dR0_dEta0_pt5, &b_Pho_Cone04NeutralHadronIso_dR0_dEta0_pt5);
   fChain->SetBranchAddress("Pho_Cone04NeutralHadronIso_dR0_dEta0_pt0_nocracks", Pho_Cone04NeutralHadronIso_dR0_dEta0_pt0_nocracks, &b_Pho_Cone04NeutralHadronIso_dR0_dEta0_pt0_nocracks);
   fChain->SetBranchAddress("Pho_Cone04NeutralHadronIso_dR0_dEta0_pt5_nocracks", Pho_Cone04NeutralHadronIso_dR0_dEta0_pt5_nocracks, &b_Pho_Cone04NeutralHadronIso_dR0_dEta0_pt5_nocracks);
   fChain->SetBranchAddress("Pho_Cone04NeutralHadronIso_dR7_dEta0_pt0", Pho_Cone04NeutralHadronIso_dR7_dEta0_pt0, &b_Pho_Cone04NeutralHadronIso_dR7_dEta0_pt0);
   fChain->SetBranchAddress("Pho_Cone04NeutralHadronIso_dR7_dEta0_pt5", Pho_Cone04NeutralHadronIso_dR7_dEta0_pt5, &b_Pho_Cone04NeutralHadronIso_dR7_dEta0_pt5);
   fChain->SetBranchAddress("Pho_Cone01NeutralHadronIso_dR0_dEta0_pt0_mvVtx", Pho_Cone01NeutralHadronIso_dR0_dEta0_pt0_mvVtx, &b_Pho_Cone01NeutralHadronIso_dR0_dEta0_pt0_mvVtx);
   fChain->SetBranchAddress("Pho_Cone02NeutralHadronIso_dR0_dEta0_pt0_mvVtx", Pho_Cone02NeutralHadronIso_dR0_dEta0_pt0_mvVtx, &b_Pho_Cone02NeutralHadronIso_dR0_dEta0_pt0_mvVtx);
   fChain->SetBranchAddress("Pho_Cone03NeutralHadronIso_dR0_dEta0_pt0_mvVtx", Pho_Cone03NeutralHadronIso_dR0_dEta0_pt0_mvVtx, &b_Pho_Cone03NeutralHadronIso_dR0_dEta0_pt0_mvVtx);
   fChain->SetBranchAddress("Pho_Cone04NeutralHadronIso_dR0_dEta0_pt0_mvVtx", Pho_Cone04NeutralHadronIso_dR0_dEta0_pt0_mvVtx, &b_Pho_Cone04NeutralHadronIso_dR0_dEta0_pt0_mvVtx);
   fChain->SetBranchAddress("Pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_dz0_old", Pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_dz0_old, &b_Pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_dz0_old);
   fChain->SetBranchAddress("Pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_PFnoPU_old", Pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_PFnoPU_old, &b_Pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_PFnoPU_old);
   fChain->SetBranchAddress("Pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_dz0_old", Pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_dz0_old, &b_Pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_dz0_old);
   fChain->SetBranchAddress("Pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_PFnoPU_old", Pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_PFnoPU_old, &b_Pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_PFnoPU_old);
   fChain->SetBranchAddress("Pho_Cone01ChargedHadronIso_dR0_dEta0_pt0_dz0", Pho_Cone01ChargedHadronIso_dR0_dEta0_pt0_dz0, &b_Pho_Cone01ChargedHadronIso_dR0_dEta0_pt0_dz0);
   fChain->SetBranchAddress("Pho_Cone01ChargedHadronIso_dR0_dEta0_pt0_dz1_dxy01", Pho_Cone01ChargedHadronIso_dR0_dEta0_pt0_dz1_dxy01, &b_Pho_Cone01ChargedHadronIso_dR0_dEta0_pt0_dz1_dxy01);
   fChain->SetBranchAddress("Pho_Cone01ChargedHadronIso_dR0_dEta0_pt0_PFnoPU", Pho_Cone01ChargedHadronIso_dR0_dEta0_pt0_PFnoPU, &b_Pho_Cone01ChargedHadronIso_dR0_dEta0_pt0_PFnoPU);
   fChain->SetBranchAddress("Pho_Cone01ChargedHadronIso_dR015_dEta0_pt0_dz0", Pho_Cone01ChargedHadronIso_dR015_dEta0_pt0_dz0, &b_Pho_Cone01ChargedHadronIso_dR015_dEta0_pt0_dz0);
   fChain->SetBranchAddress("Pho_Cone01ChargedHadronIso_dR015_dEta0_pt0_dz1_dxy01", Pho_Cone01ChargedHadronIso_dR015_dEta0_pt0_dz1_dxy01, &b_Pho_Cone01ChargedHadronIso_dR015_dEta0_pt0_dz1_dxy01);
   fChain->SetBranchAddress("Pho_Cone01ChargedHadronIso_dR015_dEta0_pt0_PFnoPU", Pho_Cone01ChargedHadronIso_dR015_dEta0_pt0_PFnoPU, &b_Pho_Cone01ChargedHadronIso_dR015_dEta0_pt0_PFnoPU);
   fChain->SetBranchAddress("Pho_Cone02ChargedHadronIso_dR0_dEta0_pt0_dz0", Pho_Cone02ChargedHadronIso_dR0_dEta0_pt0_dz0, &b_Pho_Cone02ChargedHadronIso_dR0_dEta0_pt0_dz0);
   fChain->SetBranchAddress("Pho_Cone02ChargedHadronIso_dR0_dEta0_pt0_dz1_dxy01", Pho_Cone02ChargedHadronIso_dR0_dEta0_pt0_dz1_dxy01, &b_Pho_Cone02ChargedHadronIso_dR0_dEta0_pt0_dz1_dxy01);
   fChain->SetBranchAddress("Pho_Cone02ChargedHadronIso_dR0_dEta0_pt0_PFnoPU", Pho_Cone02ChargedHadronIso_dR0_dEta0_pt0_PFnoPU, &b_Pho_Cone02ChargedHadronIso_dR0_dEta0_pt0_PFnoPU);
   fChain->SetBranchAddress("Pho_Cone02ChargedHadronIso_dR015_dEta0_pt0_dz0", Pho_Cone02ChargedHadronIso_dR015_dEta0_pt0_dz0, &b_Pho_Cone02ChargedHadronIso_dR015_dEta0_pt0_dz0);
   fChain->SetBranchAddress("Pho_Cone02ChargedHadronIso_dR015_dEta0_pt0_dz1_dxy01", Pho_Cone02ChargedHadronIso_dR015_dEta0_pt0_dz1_dxy01, &b_Pho_Cone02ChargedHadronIso_dR015_dEta0_pt0_dz1_dxy01);
   fChain->SetBranchAddress("Pho_Cone02ChargedHadronIso_dR015_dEta0_pt0_PFnoPU", Pho_Cone02ChargedHadronIso_dR015_dEta0_pt0_PFnoPU, &b_Pho_Cone02ChargedHadronIso_dR015_dEta0_pt0_PFnoPU);
   fChain->SetBranchAddress("Pho_Cone03ChargedHadronIso_dR0_dEta0_pt0_dz0", Pho_Cone03ChargedHadronIso_dR0_dEta0_pt0_dz0, &b_Pho_Cone03ChargedHadronIso_dR0_dEta0_pt0_dz0);
   fChain->SetBranchAddress("Pho_Cone03ChargedHadronIso_dR0_dEta0_pt0_dz1_dxy01", Pho_Cone03ChargedHadronIso_dR0_dEta0_pt0_dz1_dxy01, &b_Pho_Cone03ChargedHadronIso_dR0_dEta0_pt0_dz1_dxy01);
   fChain->SetBranchAddress("Pho_Cone03ChargedHadronIso_dR0_dEta0_pt0_PFnoPU", Pho_Cone03ChargedHadronIso_dR0_dEta0_pt0_PFnoPU, &b_Pho_Cone03ChargedHadronIso_dR0_dEta0_pt0_PFnoPU);
   fChain->SetBranchAddress("Pho_Cone03ChargedHadronIso_dR015_dEta0_pt0_dz0", Pho_Cone03ChargedHadronIso_dR015_dEta0_pt0_dz0, &b_Pho_Cone03ChargedHadronIso_dR015_dEta0_pt0_dz0);
   fChain->SetBranchAddress("Pho_Cone03ChargedHadronIso_dR015_dEta0_pt0_dz1_dxy01", Pho_Cone03ChargedHadronIso_dR015_dEta0_pt0_dz1_dxy01, &b_Pho_Cone03ChargedHadronIso_dR015_dEta0_pt0_dz1_dxy01);
   fChain->SetBranchAddress("Pho_Cone03ChargedHadronIso_dR015_dEta0_pt0_PFnoPU", Pho_Cone03ChargedHadronIso_dR015_dEta0_pt0_PFnoPU, &b_Pho_Cone03ChargedHadronIso_dR015_dEta0_pt0_PFnoPU);
   fChain->SetBranchAddress("Pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_dz0", Pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_dz0, &b_Pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_dz0);
   fChain->SetBranchAddress("Pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_dz1_dxy01", Pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_dz1_dxy01, &b_Pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_dz1_dxy01);
   fChain->SetBranchAddress("Pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_PFnoPU", Pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_PFnoPU, &b_Pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_PFnoPU);
   fChain->SetBranchAddress("Pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_dz0", Pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_dz0, &b_Pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_dz0);
   fChain->SetBranchAddress("Pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_dz1_dxy01", Pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_dz1_dxy01, &b_Pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_dz1_dxy01);
   fChain->SetBranchAddress("Pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_PFnoPU", Pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_PFnoPU, &b_Pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_PFnoPU);
   fChain->SetBranchAddress("NSuperClusters", &NSuperClusters, &b_NSuperClusters);
   fChain->SetBranchAddress("SCRaw", SCRaw, &b_SCRaw);
   fChain->SetBranchAddress("SCPre", SCPre, &b_SCPre);
   fChain->SetBranchAddress("SCEnergy", SCEnergy, &b_SCEnergy);
   fChain->SetBranchAddress("SCEta", SCEta, &b_SCEta);
   fChain->SetBranchAddress("SCPhi", SCPhi, &b_SCPhi);
   fChain->SetBranchAddress("SCPhiWidth", SCPhiWidth, &b_SCPhiWidth);
   fChain->SetBranchAddress("SCEtaWidth", SCEtaWidth, &b_SCEtaWidth);
   fChain->SetBranchAddress("SCBrem", SCBrem, &b_SCBrem);
   fChain->SetBranchAddress("SCR9", SCR9, &b_SCR9);
   fChain->SetBranchAddress("SCcrackcorrseed", SCcrackcorrseed, &b_SCcrackcorrseed);
   fChain->SetBranchAddress("SCcrackcorr", SCcrackcorr, &b_SCcrackcorr);
   fChain->SetBranchAddress("SClocalcorrseed", SClocalcorrseed, &b_SClocalcorrseed);
   fChain->SetBranchAddress("SClocalcorr", SClocalcorr, &b_SClocalcorr);
   fChain->SetBranchAddress("SCcrackcorrseedfactor", SCcrackcorrseedfactor, &b_SCcrackcorrseedfactor);
   fChain->SetBranchAddress("SClocalcorrseedfactor", SClocalcorrseedfactor, &b_SClocalcorrseedfactor);
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
   fChain->SetBranchAddress("JArea", JArea, &b_JArea);
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
   fChain->SetBranchAddress("JGenJetIndex", JGenJetIndex, &b_JGenJetIndex);
   fChain->SetBranchAddress("CANJets", &CANJets, &b_CANJets);
   fChain->SetBranchAddress("CAJPx", CAJPx, &b_CAJPx);
   fChain->SetBranchAddress("CAJPy", CAJPy, &b_CAJPy);
   fChain->SetBranchAddress("CAJPz", CAJPz, &b_CAJPz);
   fChain->SetBranchAddress("CAJPt", CAJPt, &b_CAJPt);
   fChain->SetBranchAddress("CAJE", CAJE, &b_CAJE);
   fChain->SetBranchAddress("CAJEt", CAJEt, &b_CAJEt);
   fChain->SetBranchAddress("CAJEta", CAJEta, &b_CAJEta);
   fChain->SetBranchAddress("CAJPhi", CAJPhi, &b_CAJPhi);
   fChain->SetBranchAddress("CAJFlavour", CAJFlavour, &b_CAJFlavour);
   fChain->SetBranchAddress("CAJScale", CAJScale, &b_CAJScale);
   fChain->SetBranchAddress("CAJL1FastJetScale", CAJL1FastJetScale, &b_CAJL1FastJetScale);
   fChain->SetBranchAddress("CAJArea", CAJArea, &b_CAJArea);
   fChain->SetBranchAddress("CAJbTagProbTkCntHighEff", CAJbTagProbTkCntHighEff, &b_CAJbTagProbTkCntHighEff);
   fChain->SetBranchAddress("CAJbTagProbTkCntHighPur", CAJbTagProbTkCntHighPur, &b_CAJbTagProbTkCntHighPur);
   fChain->SetBranchAddress("CAJbTagProbSimpSVHighEff", CAJbTagProbSimpSVHighEff, &b_CAJbTagProbSimpSVHighEff);
   fChain->SetBranchAddress("CAJbTagProbSimpSVHighPur", CAJbTagProbSimpSVHighPur, &b_CAJbTagProbSimpSVHighPur);
   fChain->SetBranchAddress("CAJIDLoose", CAJIDLoose, &b_CAJIDLoose);
   fChain->SetBranchAddress("CAJID_HPD", CAJID_HPD, &b_CAJID_HPD);
   fChain->SetBranchAddress("CAJID_RBX", CAJID_RBX, &b_CAJID_RBX);
   fChain->SetBranchAddress("CAJID_n90Hits", CAJID_n90Hits, &b_CAJID_n90Hits);
   fChain->SetBranchAddress("CAJID_resEMF", CAJID_resEMF, &b_CAJID_resEMF);
   fChain->SetBranchAddress("CAJEMfrac", CAJEMfrac, &b_CAJEMfrac);
   fChain->SetBranchAddress("CAJNAssoTracks", CAJNAssoTracks, &b_CAJNAssoTracks);
   fChain->SetBranchAddress("CAJChfrac", CAJChfrac, &b_CAJChfrac);
   fChain->SetBranchAddress("CAJNConstituents", CAJNConstituents, &b_CAJNConstituents);
   fChain->SetBranchAddress("CAJn90", CAJn90, &b_CAJn90);
   fChain->SetBranchAddress("PF2PATAntiIsoNJets", &PF2PATAntiIsoNJets, &b_PF2PATAntiIsoNJets);
   fChain->SetBranchAddress("PF2PATAntiIsoJPx", PF2PATAntiIsoJPx, &b_PF2PATAntiIsoJPx);
   fChain->SetBranchAddress("PF2PATAntiIsoJPy", PF2PATAntiIsoJPy, &b_PF2PATAntiIsoJPy);
   fChain->SetBranchAddress("PF2PATAntiIsoJPz", PF2PATAntiIsoJPz, &b_PF2PATAntiIsoJPz);
   fChain->SetBranchAddress("PF2PATAntiIsoJPt", PF2PATAntiIsoJPt, &b_PF2PATAntiIsoJPt);
   fChain->SetBranchAddress("PF2PATAntiIsoJE", PF2PATAntiIsoJE, &b_PF2PATAntiIsoJE);
   fChain->SetBranchAddress("PF2PATAntiIsoJEt", PF2PATAntiIsoJEt, &b_PF2PATAntiIsoJEt);
   fChain->SetBranchAddress("PF2PATAntiIsoJEta", PF2PATAntiIsoJEta, &b_PF2PATAntiIsoJEta);
   fChain->SetBranchAddress("PF2PATAntiIsoJPhi", PF2PATAntiIsoJPhi, &b_PF2PATAntiIsoJPhi);
   fChain->SetBranchAddress("PF2PATAntiIsoJFlavour", PF2PATAntiIsoJFlavour, &b_PF2PATAntiIsoJFlavour);
   fChain->SetBranchAddress("PF2PATAntiIsoJScale", PF2PATAntiIsoJScale, &b_PF2PATAntiIsoJScale);
   fChain->SetBranchAddress("PF2PATAntiIsoJL1FastJetScale", PF2PATAntiIsoJL1FastJetScale, &b_PF2PATAntiIsoJL1FastJetScale);
   fChain->SetBranchAddress("PF2PATAntiIsoJArea", PF2PATAntiIsoJArea, &b_PF2PATAntiIsoJArea);
   fChain->SetBranchAddress("PF2PATAntiIsoJbTagProbTkCntHighEff", PF2PATAntiIsoJbTagProbTkCntHighEff, &b_PF2PATAntiIsoJbTagProbTkCntHighEff);
   fChain->SetBranchAddress("PF2PATAntiIsoJbTagProbTkCntHighPur", PF2PATAntiIsoJbTagProbTkCntHighPur, &b_PF2PATAntiIsoJbTagProbTkCntHighPur);
   fChain->SetBranchAddress("PF2PATAntiIsoJbTagProbSimpSVHighEff", PF2PATAntiIsoJbTagProbSimpSVHighEff, &b_PF2PATAntiIsoJbTagProbSimpSVHighEff);
   fChain->SetBranchAddress("PF2PATAntiIsoJbTagProbSimpSVHighPur", PF2PATAntiIsoJbTagProbSimpSVHighPur, &b_PF2PATAntiIsoJbTagProbSimpSVHighPur);
   fChain->SetBranchAddress("PF2PATAntiIsoJIDLoose", PF2PATAntiIsoJIDLoose, &b_PF2PATAntiIsoJIDLoose);
   fChain->SetBranchAddress("PF2PATAntiIsoJChMult", PF2PATAntiIsoJChMult, &b_PF2PATAntiIsoJChMult);
   fChain->SetBranchAddress("PF2PATAntiIsoJNeuMult", PF2PATAntiIsoJNeuMult, &b_PF2PATAntiIsoJNeuMult);
   fChain->SetBranchAddress("PF2PATAntiIsoJChHadfrac", PF2PATAntiIsoJChHadfrac, &b_PF2PATAntiIsoJChHadfrac);
   fChain->SetBranchAddress("PF2PATAntiIsoJNeuHadfrac", PF2PATAntiIsoJNeuHadfrac, &b_PF2PATAntiIsoJNeuHadfrac);
   fChain->SetBranchAddress("PF2PATAntiIsoJChEmfrac", PF2PATAntiIsoJChEmfrac, &b_PF2PATAntiIsoJChEmfrac);
   fChain->SetBranchAddress("PF2PATAntiIsoJNeuEmfrac", PF2PATAntiIsoJNeuEmfrac, &b_PF2PATAntiIsoJNeuEmfrac);
   fChain->SetBranchAddress("PF2PATAntiIsoJChMufrac", PF2PATAntiIsoJChMufrac, &b_PF2PATAntiIsoJChMufrac);
   fChain->SetBranchAddress("PF2PATAntiIsoJNConstituents", PF2PATAntiIsoJNConstituents, &b_PF2PATAntiIsoJNConstituents);
   fChain->SetBranchAddress("PF2PAT2NJets", &PF2PAT2NJets, &b_PF2PAT2NJets);
   fChain->SetBranchAddress("PF2PAT2JPx", PF2PAT2JPx, &b_PF2PAT2JPx);
   fChain->SetBranchAddress("PF2PAT2JPy", PF2PAT2JPy, &b_PF2PAT2JPy);
   fChain->SetBranchAddress("PF2PAT2JPz", PF2PAT2JPz, &b_PF2PAT2JPz);
   fChain->SetBranchAddress("PF2PAT2JPt", PF2PAT2JPt, &b_PF2PAT2JPt);
   fChain->SetBranchAddress("PF2PAT2JE", PF2PAT2JE, &b_PF2PAT2JE);
   fChain->SetBranchAddress("PF2PAT2JEt", PF2PAT2JEt, &b_PF2PAT2JEt);
   fChain->SetBranchAddress("PF2PAT2JEta", PF2PAT2JEta, &b_PF2PAT2JEta);
   fChain->SetBranchAddress("PF2PAT2JPhi", PF2PAT2JPhi, &b_PF2PAT2JPhi);
   fChain->SetBranchAddress("PF2PAT2JFlavour", PF2PAT2JFlavour, &b_PF2PAT2JFlavour);
   fChain->SetBranchAddress("PF2PAT2JScale", PF2PAT2JScale, &b_PF2PAT2JScale);
   fChain->SetBranchAddress("PF2PAT2JL1FastJetScale", PF2PAT2JL1FastJetScale, &b_PF2PAT2JL1FastJetScale);
   fChain->SetBranchAddress("PF2PAT2JArea", PF2PAT2JArea, &b_PF2PAT2JArea);
   fChain->SetBranchAddress("PF2PAT2JbTagProbTkCntHighEff", PF2PAT2JbTagProbTkCntHighEff, &b_PF2PAT2JbTagProbTkCntHighEff);
   fChain->SetBranchAddress("PF2PAT2JbTagProbTkCntHighPur", PF2PAT2JbTagProbTkCntHighPur, &b_PF2PAT2JbTagProbTkCntHighPur);
   fChain->SetBranchAddress("PF2PAT2JbTagProbSimpSVHighEff", PF2PAT2JbTagProbSimpSVHighEff, &b_PF2PAT2JbTagProbSimpSVHighEff);
   fChain->SetBranchAddress("PF2PAT2JbTagProbSimpSVHighPur", PF2PAT2JbTagProbSimpSVHighPur, &b_PF2PAT2JbTagProbSimpSVHighPur);
   fChain->SetBranchAddress("PF2PAT2JIDLoose", PF2PAT2JIDLoose, &b_PF2PAT2JIDLoose);
   fChain->SetBranchAddress("PF2PAT2JChMult", PF2PAT2JChMult, &b_PF2PAT2JChMult);
   fChain->SetBranchAddress("PF2PAT2JNeuMult", PF2PAT2JNeuMult, &b_PF2PAT2JNeuMult);
   fChain->SetBranchAddress("PF2PAT2JChHadfrac", PF2PAT2JChHadfrac, &b_PF2PAT2JChHadfrac);
   fChain->SetBranchAddress("PF2PAT2JNeuHadfrac", PF2PAT2JNeuHadfrac, &b_PF2PAT2JNeuHadfrac);
   fChain->SetBranchAddress("PF2PAT2JChEmfrac", PF2PAT2JChEmfrac, &b_PF2PAT2JChEmfrac);
   fChain->SetBranchAddress("PF2PAT2JNeuEmfrac", PF2PAT2JNeuEmfrac, &b_PF2PAT2JNeuEmfrac);
   fChain->SetBranchAddress("PF2PAT2JChMufrac", PF2PAT2JChMufrac, &b_PF2PAT2JChMufrac);
   fChain->SetBranchAddress("PF2PAT2JNConstituents", PF2PAT2JNConstituents, &b_PF2PAT2JNConstituents);
   fChain->SetBranchAddress("PF2PAT3NJets", &PF2PAT3NJets, &b_PF2PAT3NJets);
   fChain->SetBranchAddress("PF2PAT3JPx", PF2PAT3JPx, &b_PF2PAT3JPx);
   fChain->SetBranchAddress("PF2PAT3JPy", PF2PAT3JPy, &b_PF2PAT3JPy);
   fChain->SetBranchAddress("PF2PAT3JPz", PF2PAT3JPz, &b_PF2PAT3JPz);
   fChain->SetBranchAddress("PF2PAT3JPt", PF2PAT3JPt, &b_PF2PAT3JPt);
   fChain->SetBranchAddress("PF2PAT3JE", PF2PAT3JE, &b_PF2PAT3JE);
   fChain->SetBranchAddress("PF2PAT3JEt", PF2PAT3JEt, &b_PF2PAT3JEt);
   fChain->SetBranchAddress("PF2PAT3JEta", PF2PAT3JEta, &b_PF2PAT3JEta);
   fChain->SetBranchAddress("PF2PAT3JPhi", PF2PAT3JPhi, &b_PF2PAT3JPhi);
   fChain->SetBranchAddress("PF2PAT3JFlavour", PF2PAT3JFlavour, &b_PF2PAT3JFlavour);
   fChain->SetBranchAddress("PF2PAT3JScale", PF2PAT3JScale, &b_PF2PAT3JScale);
   fChain->SetBranchAddress("PF2PAT3JL1FastJetScale", PF2PAT3JL1FastJetScale, &b_PF2PAT3JL1FastJetScale);
   fChain->SetBranchAddress("PF2PAT3JArea", PF2PAT3JArea, &b_PF2PAT3JArea);
   fChain->SetBranchAddress("PF2PAT3JbTagProbTkCntHighEff", PF2PAT3JbTagProbTkCntHighEff, &b_PF2PAT3JbTagProbTkCntHighEff);
   fChain->SetBranchAddress("PF2PAT3JbTagProbTkCntHighPur", PF2PAT3JbTagProbTkCntHighPur, &b_PF2PAT3JbTagProbTkCntHighPur);
   fChain->SetBranchAddress("PF2PAT3JbTagProbSimpSVHighEff", PF2PAT3JbTagProbSimpSVHighEff, &b_PF2PAT3JbTagProbSimpSVHighEff);
   fChain->SetBranchAddress("PF2PAT3JbTagProbSimpSVHighPur", PF2PAT3JbTagProbSimpSVHighPur, &b_PF2PAT3JbTagProbSimpSVHighPur);
   fChain->SetBranchAddress("PF2PAT3JIDLoose", PF2PAT3JIDLoose, &b_PF2PAT3JIDLoose);
   fChain->SetBranchAddress("PF2PAT3JChMult", PF2PAT3JChMult, &b_PF2PAT3JChMult);
   fChain->SetBranchAddress("PF2PAT3JNeuMult", PF2PAT3JNeuMult, &b_PF2PAT3JNeuMult);
   fChain->SetBranchAddress("PF2PAT3JChHadfrac", PF2PAT3JChHadfrac, &b_PF2PAT3JChHadfrac);
   fChain->SetBranchAddress("PF2PAT3JNeuHadfrac", PF2PAT3JNeuHadfrac, &b_PF2PAT3JNeuHadfrac);
   fChain->SetBranchAddress("PF2PAT3JChEmfrac", PF2PAT3JChEmfrac, &b_PF2PAT3JChEmfrac);
   fChain->SetBranchAddress("PF2PAT3JNeuEmfrac", PF2PAT3JNeuEmfrac, &b_PF2PAT3JNeuEmfrac);
   fChain->SetBranchAddress("PF2PAT3JChMufrac", PF2PAT3JChMufrac, &b_PF2PAT3JChMufrac);
   fChain->SetBranchAddress("PF2PAT3JNConstituents", PF2PAT3JNConstituents, &b_PF2PAT3JNConstituents);
   fChain->SetBranchAddress("PfMuAntiIsoMaxLepExc", &PfMuAntiIsoMaxLepExc, &b_PfMuAntiIsoMaxLepExc);
   fChain->SetBranchAddress("PfMuAntiIsoNObjsTot", &PfMuAntiIsoNObjsTot, &b_PfMuAntiIsoNObjsTot);
   fChain->SetBranchAddress("PfMuAntiIsoNObjs", &PfMuAntiIsoNObjs, &b_PfMuAntiIsoNObjs);
   fChain->SetBranchAddress("PfMuAntiIsoPx", PfMuAntiIsoPx, &b_PfMuAntiIsoPx);
   fChain->SetBranchAddress("PfMuAntiIsoPy", PfMuAntiIsoPy, &b_PfMuAntiIsoPy);
   fChain->SetBranchAddress("PfMuAntiIsoPz", PfMuAntiIsoPz, &b_PfMuAntiIsoPz);
   fChain->SetBranchAddress("PfMuAntiIsoPt", PfMuAntiIsoPt, &b_PfMuAntiIsoPt);
   fChain->SetBranchAddress("PfMuAntiIsoE", PfMuAntiIsoE, &b_PfMuAntiIsoE);
   fChain->SetBranchAddress("PfMuAntiIsoEt", PfMuAntiIsoEt, &b_PfMuAntiIsoEt);
   fChain->SetBranchAddress("PfMuAntiIsoEta", PfMuAntiIsoEta, &b_PfMuAntiIsoEta);
   fChain->SetBranchAddress("PfMuAntiIsoPhi", PfMuAntiIsoPhi, &b_PfMuAntiIsoPhi);
   fChain->SetBranchAddress("PfMuAntiIsoCharge", PfMuAntiIsoCharge, &b_PfMuAntiIsoCharge);
   fChain->SetBranchAddress("PfMuAntiIsoParticleIso", PfMuAntiIsoParticleIso, &b_PfMuAntiIsoParticleIso);
   fChain->SetBranchAddress("PfMuAntiIsoChargedHadronIso", PfMuAntiIsoChargedHadronIso, &b_PfMuAntiIsoChargedHadronIso);
   fChain->SetBranchAddress("PfMuAntiIsoNeutralHadronIso", PfMuAntiIsoNeutralHadronIso, &b_PfMuAntiIsoNeutralHadronIso);
   fChain->SetBranchAddress("PfMuAntiIsoPhotonIso", PfMuAntiIsoPhotonIso, &b_PfMuAntiIsoPhotonIso);
   fChain->SetBranchAddress("PfMuAntiIsoPtErr", PfMuAntiIsoPtErr, &b_PfMuAntiIsoPtErr);
   fChain->SetBranchAddress("PfMuAntiIsoNMatches", PfMuAntiIsoNMatches, &b_PfMuAntiIsoNMatches);
   fChain->SetBranchAddress("PfMu2MaxLepExc", &PfMu2MaxLepExc, &b_PfMu2MaxLepExc);
   fChain->SetBranchAddress("PfMu2NObjsTot", &PfMu2NObjsTot, &b_PfMu2NObjsTot);
   fChain->SetBranchAddress("PfMu2NObjs", &PfMu2NObjs, &b_PfMu2NObjs);
   fChain->SetBranchAddress("PfMu2Px", &PfMu2Px, &b_PfMu2Px);
   fChain->SetBranchAddress("PfMu2Py", &PfMu2Py, &b_PfMu2Py);
   fChain->SetBranchAddress("PfMu2Pz", &PfMu2Pz, &b_PfMu2Pz);
   fChain->SetBranchAddress("PfMu2Pt", &PfMu2Pt, &b_PfMu2Pt);
   fChain->SetBranchAddress("PfMu2E", &PfMu2E, &b_PfMu2E);
   fChain->SetBranchAddress("PfMu2Et", &PfMu2Et, &b_PfMu2Et);
   fChain->SetBranchAddress("PfMu2Eta", &PfMu2Eta, &b_PfMu2Eta);
   fChain->SetBranchAddress("PfMu2Phi", &PfMu2Phi, &b_PfMu2Phi);
   fChain->SetBranchAddress("PfMu2Charge", &PfMu2Charge, &b_PfMu2Charge);
   fChain->SetBranchAddress("PfMu2ParticleIso", &PfMu2ParticleIso, &b_PfMu2ParticleIso);
   fChain->SetBranchAddress("PfMu2ChargedHadronIso", &PfMu2ChargedHadronIso, &b_PfMu2ChargedHadronIso);
   fChain->SetBranchAddress("PfMu2NeutralHadronIso", &PfMu2NeutralHadronIso, &b_PfMu2NeutralHadronIso);
   fChain->SetBranchAddress("PfMu2PhotonIso", &PfMu2PhotonIso, &b_PfMu2PhotonIso);
   fChain->SetBranchAddress("PfMu2PtErr", &PfMu2PtErr, &b_PfMu2PtErr);
   fChain->SetBranchAddress("PfMu2NMatches", &PfMu2NMatches, &b_PfMu2NMatches);
   fChain->SetBranchAddress("PfMu3MaxLepExc", &PfMu3MaxLepExc, &b_PfMu3MaxLepExc);
   fChain->SetBranchAddress("PfMu3NObjsTot", &PfMu3NObjsTot, &b_PfMu3NObjsTot);
   fChain->SetBranchAddress("PfMu3NObjs", &PfMu3NObjs, &b_PfMu3NObjs);
   fChain->SetBranchAddress("PfMu3Px", &PfMu3Px, &b_PfMu3Px);
   fChain->SetBranchAddress("PfMu3Py", &PfMu3Py, &b_PfMu3Py);
   fChain->SetBranchAddress("PfMu3Pz", &PfMu3Pz, &b_PfMu3Pz);
   fChain->SetBranchAddress("PfMu3Pt", &PfMu3Pt, &b_PfMu3Pt);
   fChain->SetBranchAddress("PfMu3E", &PfMu3E, &b_PfMu3E);
   fChain->SetBranchAddress("PfMu3Et", &PfMu3Et, &b_PfMu3Et);
   fChain->SetBranchAddress("PfMu3Eta", &PfMu3Eta, &b_PfMu3Eta);
   fChain->SetBranchAddress("PfMu3Phi", &PfMu3Phi, &b_PfMu3Phi);
   fChain->SetBranchAddress("PfMu3Charge", &PfMu3Charge, &b_PfMu3Charge);
   fChain->SetBranchAddress("PfMu3ParticleIso", &PfMu3ParticleIso, &b_PfMu3ParticleIso);
   fChain->SetBranchAddress("PfMu3ChargedHadronIso", &PfMu3ChargedHadronIso, &b_PfMu3ChargedHadronIso);
   fChain->SetBranchAddress("PfMu3NeutralHadronIso", &PfMu3NeutralHadronIso, &b_PfMu3NeutralHadronIso);
   fChain->SetBranchAddress("PfMu3PhotonIso", &PfMu3PhotonIso, &b_PfMu3PhotonIso);
   fChain->SetBranchAddress("PfMu3PtErr", &PfMu3PtErr, &b_PfMu3PtErr);
   fChain->SetBranchAddress("PfMu3NMatches", &PfMu3NMatches, &b_PfMu3NMatches);
   fChain->SetBranchAddress("PfElAntiIsoMaxLepExc", &PfElAntiIsoMaxLepExc, &b_PfElAntiIsoMaxLepExc);
   fChain->SetBranchAddress("PfElAntiIsoNObjsTot", &PfElAntiIsoNObjsTot, &b_PfElAntiIsoNObjsTot);
   fChain->SetBranchAddress("PfElAntiIsoNObjs", &PfElAntiIsoNObjs, &b_PfElAntiIsoNObjs);
   fChain->SetBranchAddress("PfElAntiIsoPx", PfElAntiIsoPx, &b_PfElAntiIsoPx);
   fChain->SetBranchAddress("PfElAntiIsoPy", PfElAntiIsoPy, &b_PfElAntiIsoPy);
   fChain->SetBranchAddress("PfElAntiIsoPz", PfElAntiIsoPz, &b_PfElAntiIsoPz);
   fChain->SetBranchAddress("PfElAntiIsoPt", PfElAntiIsoPt, &b_PfElAntiIsoPt);
   fChain->SetBranchAddress("PfElAntiIsoE", PfElAntiIsoE, &b_PfElAntiIsoE);
   fChain->SetBranchAddress("PfElAntiIsoEt", PfElAntiIsoEt, &b_PfElAntiIsoEt);
   fChain->SetBranchAddress("PfElAntiIsoEta", PfElAntiIsoEta, &b_PfElAntiIsoEta);
   fChain->SetBranchAddress("PfElAntiIsoPhi", PfElAntiIsoPhi, &b_PfElAntiIsoPhi);
   fChain->SetBranchAddress("PfElAntiIsoCharge", PfElAntiIsoCharge, &b_PfElAntiIsoCharge);
   fChain->SetBranchAddress("PfElAntiIsoParticleIso", PfElAntiIsoParticleIso, &b_PfElAntiIsoParticleIso);
   fChain->SetBranchAddress("PfElAntiIsoChargedHadronIso", PfElAntiIsoChargedHadronIso, &b_PfElAntiIsoChargedHadronIso);
   fChain->SetBranchAddress("PfElAntiIsoNeutralHadronIso", PfElAntiIsoNeutralHadronIso, &b_PfElAntiIsoNeutralHadronIso);
   fChain->SetBranchAddress("PfElAntiIsoPhotonIso", PfElAntiIsoPhotonIso, &b_PfElAntiIsoPhotonIso);
   fChain->SetBranchAddress("PfElAntiIsoID95", PfElAntiIsoID95, &b_PfElAntiIsoID95);
   fChain->SetBranchAddress("PfElAntiIsoID90", PfElAntiIsoID90, &b_PfElAntiIsoID90);
   fChain->SetBranchAddress("PfElAntiIsoID85", PfElAntiIsoID85, &b_PfElAntiIsoID85);
   fChain->SetBranchAddress("PfElAntiIsoID80", PfElAntiIsoID80, &b_PfElAntiIsoID80);
   fChain->SetBranchAddress("PfEl2MaxLepExc", &PfEl2MaxLepExc, &b_PfEl2MaxLepExc);
   fChain->SetBranchAddress("PfEl2NObjsTot", &PfEl2NObjsTot, &b_PfEl2NObjsTot);
   fChain->SetBranchAddress("PfEl2NObjs", &PfEl2NObjs, &b_PfEl2NObjs);
   fChain->SetBranchAddress("PfEl2Px", PfEl2Px, &b_PfEl2Px);
   fChain->SetBranchAddress("PfEl2Py", PfEl2Py, &b_PfEl2Py);
   fChain->SetBranchAddress("PfEl2Pz", PfEl2Pz, &b_PfEl2Pz);
   fChain->SetBranchAddress("PfEl2Pt", PfEl2Pt, &b_PfEl2Pt);
   fChain->SetBranchAddress("PfEl2E", PfEl2E, &b_PfEl2E);
   fChain->SetBranchAddress("PfEl2Et", PfEl2Et, &b_PfEl2Et);
   fChain->SetBranchAddress("PfEl2Eta", PfEl2Eta, &b_PfEl2Eta);
   fChain->SetBranchAddress("PfEl2Phi", PfEl2Phi, &b_PfEl2Phi);
   fChain->SetBranchAddress("PfEl2Charge", PfEl2Charge, &b_PfEl2Charge);
   fChain->SetBranchAddress("PfEl2ParticleIso", PfEl2ParticleIso, &b_PfEl2ParticleIso);
   fChain->SetBranchAddress("PfEl2ChargedHadronIso", PfEl2ChargedHadronIso, &b_PfEl2ChargedHadronIso);
   fChain->SetBranchAddress("PfEl2NeutralHadronIso", PfEl2NeutralHadronIso, &b_PfEl2NeutralHadronIso);
   fChain->SetBranchAddress("PfEl2PhotonIso", PfEl2PhotonIso, &b_PfEl2PhotonIso);
   fChain->SetBranchAddress("PfEl2ID95", PfEl2ID95, &b_PfEl2ID95);
   fChain->SetBranchAddress("PfEl2ID90", PfEl2ID90, &b_PfEl2ID90);
   fChain->SetBranchAddress("PfEl2ID85", PfEl2ID85, &b_PfEl2ID85);
   fChain->SetBranchAddress("PfEl2ID80", PfEl2ID80, &b_PfEl2ID80);
   fChain->SetBranchAddress("PfEl3MaxLepExc", &PfEl3MaxLepExc, &b_PfEl3MaxLepExc);
   fChain->SetBranchAddress("PfEl3NObjsTot", &PfEl3NObjsTot, &b_PfEl3NObjsTot);
   fChain->SetBranchAddress("PfEl3NObjs", &PfEl3NObjs, &b_PfEl3NObjs);
   fChain->SetBranchAddress("PfEl3Px", PfEl3Px, &b_PfEl3Px);
   fChain->SetBranchAddress("PfEl3Py", PfEl3Py, &b_PfEl3Py);
   fChain->SetBranchAddress("PfEl3Pz", PfEl3Pz, &b_PfEl3Pz);
   fChain->SetBranchAddress("PfEl3Pt", PfEl3Pt, &b_PfEl3Pt);
   fChain->SetBranchAddress("PfEl3E", PfEl3E, &b_PfEl3E);
   fChain->SetBranchAddress("PfEl3Et", PfEl3Et, &b_PfEl3Et);
   fChain->SetBranchAddress("PfEl3Eta", PfEl3Eta, &b_PfEl3Eta);
   fChain->SetBranchAddress("PfEl3Phi", PfEl3Phi, &b_PfEl3Phi);
   fChain->SetBranchAddress("PfEl3Charge", PfEl3Charge, &b_PfEl3Charge);
   fChain->SetBranchAddress("PfEl3ParticleIso", PfEl3ParticleIso, &b_PfEl3ParticleIso);
   fChain->SetBranchAddress("PfEl3ChargedHadronIso", PfEl3ChargedHadronIso, &b_PfEl3ChargedHadronIso);
   fChain->SetBranchAddress("PfEl3NeutralHadronIso", PfEl3NeutralHadronIso, &b_PfEl3NeutralHadronIso);
   fChain->SetBranchAddress("PfEl3PhotonIso", PfEl3PhotonIso, &b_PfEl3PhotonIso);
   fChain->SetBranchAddress("PfEl3ID95", PfEl3ID95, &b_PfEl3ID95);
   fChain->SetBranchAddress("PfEl3ID90", PfEl3ID90, &b_PfEl3ID90);
   fChain->SetBranchAddress("PfEl3ID85", PfEl3ID85, &b_PfEl3ID85);
   fChain->SetBranchAddress("PfEl3ID80", PfEl3ID80, &b_PfEl3ID80);
   fChain->SetBranchAddress("PfTauAntiIsoMaxLepExc", &PfTauAntiIsoMaxLepExc, &b_PfTauAntiIsoMaxLepExc);
   fChain->SetBranchAddress("PfTauAntiIsoNObjsTot", &PfTauAntiIsoNObjsTot, &b_PfTauAntiIsoNObjsTot);
   fChain->SetBranchAddress("PfTauAntiIsoNObjs", &PfTauAntiIsoNObjs, &b_PfTauAntiIsoNObjs);
   fChain->SetBranchAddress("PfTauAntiIsoPx", PfTauAntiIsoPx, &b_PfTauAntiIsoPx);
   fChain->SetBranchAddress("PfTauAntiIsoPy", PfTauAntiIsoPy, &b_PfTauAntiIsoPy);
   fChain->SetBranchAddress("PfTauAntiIsoPz", PfTauAntiIsoPz, &b_PfTauAntiIsoPz);
   fChain->SetBranchAddress("PfTauAntiIsoPt", PfTauAntiIsoPt, &b_PfTauAntiIsoPt);
   fChain->SetBranchAddress("PfTauAntiIsoE", PfTauAntiIsoE, &b_PfTauAntiIsoE);
   fChain->SetBranchAddress("PfTauAntiIsoEt", PfTauAntiIsoEt, &b_PfTauAntiIsoEt);
   fChain->SetBranchAddress("PfTauAntiIsoEta", PfTauAntiIsoEta, &b_PfTauAntiIsoEta);
   fChain->SetBranchAddress("PfTauAntiIsoPhi", PfTauAntiIsoPhi, &b_PfTauAntiIsoPhi);
   fChain->SetBranchAddress("PfTauAntiIsoCharge", PfTauAntiIsoCharge, &b_PfTauAntiIsoCharge);
   fChain->SetBranchAddress("PfTauAntiIsoParticleIso", PfTauAntiIsoParticleIso, &b_PfTauAntiIsoParticleIso);
   fChain->SetBranchAddress("PfTauAntiIsoChargedHadronIso", PfTauAntiIsoChargedHadronIso, &b_PfTauAntiIsoChargedHadronIso);
   fChain->SetBranchAddress("PfTauAntiIsoNeutralHadronIso", PfTauAntiIsoNeutralHadronIso, &b_PfTauAntiIsoNeutralHadronIso);
   fChain->SetBranchAddress("PfTauAntiIsoPhotonIso", PfTauAntiIsoPhotonIso, &b_PfTauAntiIsoPhotonIso);
   fChain->SetBranchAddress("PfTauAntiIsoDecayMode", PfTauAntiIsoDecayMode, &b_PfTauAntiIsoDecayMode);
   fChain->SetBranchAddress("PfTau2MaxLepExc", &PfTau2MaxLepExc, &b_PfTau2MaxLepExc);
   fChain->SetBranchAddress("PfTau2NObjsTot", &PfTau2NObjsTot, &b_PfTau2NObjsTot);
   fChain->SetBranchAddress("PfTau2NObjs", &PfTau2NObjs, &b_PfTau2NObjs);
   fChain->SetBranchAddress("PfTau2Px", PfTau2Px, &b_PfTau2Px);
   fChain->SetBranchAddress("PfTau2Py", PfTau2Py, &b_PfTau2Py);
   fChain->SetBranchAddress("PfTau2Pz", PfTau2Pz, &b_PfTau2Pz);
   fChain->SetBranchAddress("PfTau2Pt", PfTau2Pt, &b_PfTau2Pt);
   fChain->SetBranchAddress("PfTau2E", PfTau2E, &b_PfTau2E);
   fChain->SetBranchAddress("PfTau2Et", PfTau2Et, &b_PfTau2Et);
   fChain->SetBranchAddress("PfTau2Eta", PfTau2Eta, &b_PfTau2Eta);
   fChain->SetBranchAddress("PfTau2Phi", PfTau2Phi, &b_PfTau2Phi);
   fChain->SetBranchAddress("PfTau2Charge", PfTau2Charge, &b_PfTau2Charge);
   fChain->SetBranchAddress("PfTau2ParticleIso", PfTau2ParticleIso, &b_PfTau2ParticleIso);
   fChain->SetBranchAddress("PfTau2ChargedHadronIso", PfTau2ChargedHadronIso, &b_PfTau2ChargedHadronIso);
   fChain->SetBranchAddress("PfTau2NeutralHadronIso", PfTau2NeutralHadronIso, &b_PfTau2NeutralHadronIso);
   fChain->SetBranchAddress("PfTau2PhotonIso", PfTau2PhotonIso, &b_PfTau2PhotonIso);
   fChain->SetBranchAddress("PfTau2DecayMode", PfTau2DecayMode, &b_PfTau2DecayMode);
   fChain->SetBranchAddress("PfTau3MaxLepExc", &PfTau3MaxLepExc, &b_PfTau3MaxLepExc);
   fChain->SetBranchAddress("PfTau3NObjsTot", &PfTau3NObjsTot, &b_PfTau3NObjsTot);
   fChain->SetBranchAddress("PfTau3NObjs", &PfTau3NObjs, &b_PfTau3NObjs);
   fChain->SetBranchAddress("PfTau3Px", PfTau3Px, &b_PfTau3Px);
   fChain->SetBranchAddress("PfTau3Py", PfTau3Py, &b_PfTau3Py);
   fChain->SetBranchAddress("PfTau3Pz", PfTau3Pz, &b_PfTau3Pz);
   fChain->SetBranchAddress("PfTau3Pt", PfTau3Pt, &b_PfTau3Pt);
   fChain->SetBranchAddress("PfTau3E", PfTau3E, &b_PfTau3E);
   fChain->SetBranchAddress("PfTau3Et", PfTau3Et, &b_PfTau3Et);
   fChain->SetBranchAddress("PfTau3Eta", PfTau3Eta, &b_PfTau3Eta);
   fChain->SetBranchAddress("PfTau3Phi", PfTau3Phi, &b_PfTau3Phi);
   fChain->SetBranchAddress("PfTau3Charge", PfTau3Charge, &b_PfTau3Charge);
   fChain->SetBranchAddress("PfTau3ParticleIso", PfTau3ParticleIso, &b_PfTau3ParticleIso);
   fChain->SetBranchAddress("PfTau3ChargedHadronIso", PfTau3ChargedHadronIso, &b_PfTau3ChargedHadronIso);
   fChain->SetBranchAddress("PfTau3NeutralHadronIso", PfTau3NeutralHadronIso, &b_PfTau3NeutralHadronIso);
   fChain->SetBranchAddress("PfTau3PhotonIso", PfTau3PhotonIso, &b_PfTau3PhotonIso);
   fChain->SetBranchAddress("PfTau3DecayMode", PfTau3DecayMode, &b_PfTau3DecayMode);
   fChain->SetBranchAddress("NTracks", &NTracks, &b_NTracks);
   fChain->SetBranchAddress("NTracksTot", &NTracksTot, &b_NTracksTot);
   fChain->SetBranchAddress("TrkGood", TrkGood, &b_TrkGood);
   fChain->SetBranchAddress("TrkPt", TrkPt, &b_TrkPt);
   fChain->SetBranchAddress("TrkEta", TrkEta, &b_TrkEta);
   fChain->SetBranchAddress("TrkPhi", TrkPhi, &b_TrkPhi);
   fChain->SetBranchAddress("TrkNChi2", TrkNChi2, &b_TrkNChi2);
   fChain->SetBranchAddress("TrkNHits", TrkNHits, &b_TrkNHits);
   fChain->SetBranchAddress("TrkVtxDz", TrkVtxDz, &b_TrkVtxDz);
   fChain->SetBranchAddress("TrkVtxDxy", TrkVtxDxy, &b_TrkVtxDxy);
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
   fChain->SetBranchAddress("PFSumEt", &PFSumEt, &b_PFSumEt);
   fChain->SetBranchAddress("PFMETPAT", &PFMETPAT, &b_PFMETPAT);
   fChain->SetBranchAddress("PFMETPATpx", &PFMETPATpx, &b_PFMETPATpx);
   fChain->SetBranchAddress("PFMETPATpy", &PFMETPATpy, &b_PFMETPATpy);
   fChain->SetBranchAddress("PFMETPATphi", &PFMETPATphi, &b_PFMETPATphi);
   fChain->SetBranchAddress("PFMETPATSignificance", &PFMETPATSignificance, &b_PFMETPATSignificance);
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
