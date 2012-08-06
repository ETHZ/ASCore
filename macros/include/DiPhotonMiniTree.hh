#ifndef DiPhotonMiniTree_hh
#define DiPhotonMiniTree_hh


#include <cmath>
#include <string>
#include <iostream>
#include <fstream>

#include <TH1.h>

#include "base/TreeReader.hh"
#include "base/UserAnalysisBase.hh"

#include "TLorentzVector.h"

#include "EnergyCorrection.hh"

#include <vector>

#include "TRandom3.h"
#include "helper/Utilities.hh"

#include "Math/Vector3D.h"
#include "TVector3.h"

class DiPhotonMiniTree : public UserAnalysisBase{
public:
  DiPhotonMiniTree(TreeReader *tr = NULL, std::string dataType="data", Float_t aw=-999, Float_t* _kfac=NULL);
  virtual ~DiPhotonMiniTree();

  void Begin();
  void Analyze();
  void End();

private:

  float eff_area_EB;
  float eff_area_EE;

  Float_t* kfactors;

  EnergyCorrection *phocorr;
  TLorentzVector CorrPhoton(TreeReader *fTR, int i, int mode);
	
  std::vector<int> ApplyPixelVeto(TreeReader *fTR, vector<int> passing, bool forelectron=0);
  std::vector<int> PhotonPreSelection(TreeReader *fTR, vector<int> passing);
  std::vector<int> GenLevelIsolationCut(TreeReader *fTR, vector<int> passing);
  std::vector<int> PhotonSelection(TreeReader *fTR, vector<int> passing);
  std::vector<int> SignalSelection(TreeReader *fTR, vector<int> passing);
  std::vector<int> BackgroundSelection(TreeReader *fTR, vector<int> passing);
  std::vector<int> ImpingingTrackSelection(TreeReader *fTR, std::vector<int> passing, bool invert=false);
  bool FindImpingingTrack(TreeReader *fTR, TVector3 caloposition, TVector3 refvertex, int &reference_index_found, bool invert, std::vector<int> removals);
  std::vector<int> GetPFCandRemovals(TreeReader *fTR, int phoqi);
  bool SinglePhotonEventSelection(TreeReader *fTR, std::vector<int> passing);
  bool StandardEventSelection(TreeReader *fTR, std::vector<int> passing);
  bool TriggerSelection();
  int Count_part_isrfsr_gamma(TreeReader *fTR, vector<int> passing);
  void ResetVars();
  void Fillhist_PFPhotonDepositAroundImpingingTrack(int phoqi, int trkindex);  

  void FillLead(int index);
  void FillTrail(int index);

  //  double etaTransformation(float EtaParticle, float Zvertex);
  double phiNorm(float phi);
  bool FindCloseJetsAndPhotons(TreeReader *fTR, float eta, float phi, int phoqi);
  float RandomConePhotonIsolation(TreeReader *fTR, int phoqi);
  float CombinedPFIsolation(float eta, float phi, int phoqi);

  int FindPFCandType(int id);

  TRandom3 *randomgen;
  
  std::string fDataType_;
  bool isdata;

  TH1F *histo_PFPhotonDepositAroundImpingingTrack;

  Float_t AddWeight;

  int impinging_track_pfcand[100];

  TTree* OutputTree[7];

  Float_t event_luminormfactor;
  Float_t event_Kfactor;
  Float_t event_weight;

  Float_t event_rho;
  Int_t event_nPU;
  Int_t event_nRecVtx;
  
 
  Float_t dipho_mgg_photon;
  Float_t dipho_mgg_newCorrNoCrack;
  Float_t dipho_mgg_newCorr;
  Float_t dipho_mgg_newCorrLocal;
  
  Float_t pholead_eta, photrail_eta;
  Float_t pholead_px, photrail_px;
  Float_t pholead_py, photrail_py;
  Float_t pholead_pt, photrail_pt;
  Float_t pholead_pz, photrail_pz;
  Float_t pholead_energy, photrail_energy;

  Float_t pholead_SCeta, photrail_SCeta;
  Float_t pholead_SCphi, photrail_SCphi;
 
  Int_t pholead_PhoHasPixSeed, photrail_PhoHasPixSeed;
  Int_t pholead_PhoHasConvTrks, photrail_PhoHasConvTrks;
  Int_t pholead_PhoScSeedSeverity, photrail_PhoScSeedSeverity;
 
  Float_t pholead_energySCdefault, photrail_energySCdefault;
  Float_t pholead_energyNewCorr, photrail_energyNewCorr;
  Float_t pholead_energyNewCorrLocal, photrail_energyNewCorrLocal;

  Float_t pholead_r9, photrail_r9;
  Float_t pholead_sieie, photrail_sieie;
  Float_t pholead_hoe, photrail_hoe;
  Float_t pholead_brem, photrail_brem;
  Float_t pholead_sigmaPhi, photrail_sigmaPhi;
  Float_t pholead_sigmaEta, photrail_sigmaEta;
  

  Float_t pholead_pho_Cone01PhotonIso_dEta015EB_dR070EE_mvVtx;
  Float_t photrail_pho_Cone01PhotonIso_dEta015EB_dR070EE_mvVtx;
  Float_t pholead_pho_Cone02PhotonIso_dEta015EB_dR070EE_mvVtx;
  Float_t photrail_pho_Cone02PhotonIso_dEta015EB_dR070EE_mvVtx;
  Float_t pholead_pho_Cone03PhotonIso_dEta015EB_dR070EE_mvVtx;
  Float_t photrail_pho_Cone03PhotonIso_dEta015EB_dR070EE_mvVtx;
  Float_t pholead_pho_Cone04PhotonIso_dEta015EB_dR070EE_mvVtx;
  Float_t photrail_pho_Cone04PhotonIso_dEta015EB_dR070EE_mvVtx;
  Float_t pholead_pho_Cone01NeutralHadronIso_mvVtx;
  Float_t photrail_pho_Cone01NeutralHadronIso_mvVtx;
  Float_t pholead_pho_Cone02NeutralHadronIso_mvVtx;
  Float_t photrail_pho_Cone02NeutralHadronIso_mvVtx;
  Float_t pholead_pho_Cone03NeutralHadronIso_mvVtx;
  Float_t photrail_pho_Cone03NeutralHadronIso_mvVtx;
  Float_t pholead_pho_Cone04NeutralHadronIso_mvVtx;
  Float_t photrail_pho_Cone04NeutralHadronIso_mvVtx;
  Float_t pholead_pho_Cone01ChargedHadronIso_dR02_dz02_dxy01;
  Float_t photrail_pho_Cone01ChargedHadronIso_dR02_dz02_dxy01;
  Float_t pholead_pho_Cone02ChargedHadronIso_dR02_dz02_dxy01;
  Float_t photrail_pho_Cone02ChargedHadronIso_dR02_dz02_dxy01;
  Float_t pholead_pho_Cone03ChargedHadronIso_dR02_dz02_dxy01;
  Float_t photrail_pho_Cone03ChargedHadronIso_dR02_dz02_dxy01;
  Float_t pholead_pho_Cone04ChargedHadronIso_dR02_dz02_dxy01;
  Float_t photrail_pho_Cone04ChargedHadronIso_dR02_dz02_dxy01;
  Float_t pholead_pho_Cone03PFCombinedIso;
  Float_t photrail_pho_Cone03PFCombinedIso;
  Float_t pholead_pho_Cone04PFCombinedIso;
  Float_t photrail_pho_Cone04PFCombinedIso;
  Int_t pholead_PhoPassConvSafeElectronVeto;
  Int_t photrail_PhoPassConvSafeElectronVeto;
  Float_t pholead_GenPhotonIsoDR04;
  Float_t photrail_GenPhotonIsoDR04;
 
  Float_t  pholead_PhoIso03Ecal, photrail_PhoIso03Ecal;
  Float_t  pholead_PhoIso03Hcal, photrail_PhoIso03Hcal;
  Float_t  pholead_PhoIso03TrkSolid, photrail_PhoIso03TrkSolid;
  Float_t  pholead_PhoIso03TrkHollow, photrail_PhoIso03TrkHollow;
  Float_t  pholead_PhoIso03, photrail_PhoIso03;
  Float_t  pholead_PhoIso04Ecal, photrail_PhoIso04Ecal;
  Float_t  pholead_PhoIso04Hcal, photrail_PhoIso04Hcal;
  Float_t  pholead_PhoIso04TrkSolid, photrail_PhoIso04TrkSolid;
  Float_t  pholead_PhoIso04TrkHollow, photrail_PhoIso04TrkHollow;
  Float_t  pholead_PhoIso04, photrail_PhoIso04;


  Float_t pholead_PhoE1OverE9, photrail_PhoE1OverE9;
  Float_t pholead_PhoS4OverS1, photrail_PhoS4OverS1;
  Float_t pholead_PhoSigmaEtaEta, photrail_PhoSigmaEtaEta;
  Float_t pholead_PhoE1x5, photrail_PhoE1x5;
  Float_t pholead_PhoE2x5, photrail_PhoE2x5;
  Float_t pholead_PhoE3x3, photrail_PhoE3x3;
  Float_t pholead_PhoE5x5, photrail_PhoE5x5;
  Float_t pholead_PhomaxEnergyXtal, photrail_PhomaxEnergyXtal;
  Float_t pholead_PhoIso03HcalDepth1, photrail_PhoIso03HcalDepth1;
  Float_t pholead_PhoIso03HcalDepth2, photrail_PhoIso03HcalDepth2;
  Float_t pholead_PhoIso04HcalDepth1, photrail_PhoIso04HcalDepth1;
  Float_t pholead_PhoIso04HcalDepth2, photrail_PhoIso04HcalDepth2;
  Int_t pholead_PhoIso03nTrksSolid, photrail_PhoIso03nTrksSolid;
  Int_t pholead_PhoIso03nTrksHollow, photrail_PhoIso03nTrksHollow;
  Int_t pholead_PhoIso04nTrksSolid, photrail_PhoIso04nTrksSolid;
  Int_t pholead_PhoIso04nTrksHollow, photrail_PhoIso04nTrksHollow;
  Float_t pholead_Pho_ChargedHadronIso, photrail_Pho_ChargedHadronIso;
  Float_t pholead_Pho_NeutralHadronIso, photrail_Pho_NeutralHadronIso;
  Float_t pholead_Pho_PhotonIso, photrail_Pho_PhotonIso;
  Int_t pholead_Pho_isPFPhoton, photrail_Pho_isPFPhoton;
  Int_t pholead_Pho_isPFElectron, photrail_Pho_isPFElectron;

  Int_t pholead_hasimpingingtrack, photrail_hasimpingingtrack;

  Int_t pholead_PhoMCmatchindex, photrail_PhoMCmatchindex;
  Int_t pholead_PhoMCmatchexitcode, photrail_PhoMCmatchexitcode;




  TH1F *fHNumPU;
  TH1F *fHNumVtx;

};
#endif
