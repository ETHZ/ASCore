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

#include "TGeoPara.h"
#include "TGeoTube.h"

const int global_size_pfcandarrays = 30;

typedef struct {
  float photon;
  float charged;
  float neutral;
  int nphotoncand;
  int nchargedcand;
  int nneutralcand;
  std::vector<float> photoncandenergies;
  std::vector<float> chargedcandenergies;
  std::vector<float> neutralcandenergies;
  std::vector<float> photoncandets;
  std::vector<float> chargedcandets;
  std::vector<float> neutralcandets;
  std::vector<float> photoncanddetas;
  std::vector<float> chargedcanddetas;
  std::vector<float> neutralcanddetas;
  std::vector<float> photoncanddphis;
  std::vector<float> chargedcanddphis;
  std::vector<float> neutralcanddphis;
} isolations_struct;

typedef struct {
  float dR;
  float dEta;
  float dPhi;
} angular_distances_struct;

class DiPhotonMiniTree : public UserAnalysisBase{
public:
  DiPhotonMiniTree(TreeReader *tr = NULL, std::string dataType="data", Float_t aw=-999, Float_t* _kfac=NULL, Float_t _minthrpfphotoncandEB=0, Float_t _minthrpfphotoncandEE=0);
  virtual ~DiPhotonMiniTree();

  void Begin();
  void Analyze();
  void End();

private:

  float global_linkbyrechit_enlargement;
  float global_minthrpfphotoncandEB;
  float global_minthrpfphotoncandEE;

  TGeoPara eegeom;

  Float_t* kfactors;

  EnergyCorrection *phocorr;
  TLorentzVector CorrPhoton(TreeReader *fTR, int i, int mode);
	
  std::vector<int> ApplyPixelVeto(TreeReader *fTR, vector<int> passing, bool forelectron=0);
  std::vector<int> PhotonPreSelection(TreeReader *fTR, vector<int> passing);
  std::vector<int> GenLevelIsolationCut(TreeReader *fTR, vector<int> passing);
  std::vector<int> PhotonSelection(TreeReader *fTR, vector<int> passing, TString mode="");
  std::vector<int> SignalSelection(TreeReader *fTR, vector<int> passing);
  std::vector<int> BackgroundSelection(TreeReader *fTR, vector<int> passing);
  std::vector<int> ImpingingTrackSelection(TreeReader *fTR, std::vector<int> passing, bool invert=false);
  std::vector<int> GetPFCandRemovals(TreeReader *fTR, int phoqi);
  bool SinglePhotonEventSelection(TreeReader *fTR, std::vector<int> &passing);
  bool StandardEventSelection(TreeReader *fTR, std::vector<int> &passing);
  bool TriggerSelection();
  int Count_part_isrfsr_gamma(TreeReader *fTR, vector<int> passing);
  void ResetVars();
  void Fillhist_PFPhotonDepositAroundImpingingTrack(int phoqi, int trkindex);  
  std::vector<int> GetPFCandInsideFootprint(TreeReader *fTR, int phoqi, float rotation_phi, TString component);
  std::vector<int> GetPFCandWithFootprintRemoval(TreeReader *fTR, int phoqi, float rotation_phi, bool outoffootprint, TString component);
  TVector3 PropagatePFCandToEcal(int pfcandindex, float position, bool isbarrel);
  bool FindImpingingTrack(TreeReader *fTR, int phoqi, int &reference_index_found, bool dofootprintremoval = false, std::vector<int> removals = std::vector<int>());
  float PFIsolation(int phoqi, float rotation_phi, TString component, int *counter = NULL, std::vector<float> *energies = NULL, std::vector<float> *ets = NULL,  std::vector<float> *detas = NULL, std::vector<float> *dphis = NULL, std::vector<int> removals = std::vector<int>());
  angular_distances_struct GetPFCandDeltaRFromSC(TreeReader *fTR, int phoqi, int pfindex, float rotation_phi = 0);
  bool FindCloseJetsAndPhotons(TreeReader *fTR, float rotation_phi, int phoqi, TString mod="");
  std::vector<int> GetPFCandIDedRemovals(TreeReader *fTR, int phoqi);

  std::vector<int> DiPhotonInvariantMassCutSelection(TreeReader *fTR, std::vector<int> passing);

  float SieieRescale(float sieie, bool isbarrel);
  float CalculateSCArea(TreeReader *fTR, int scindex);
  float GetPUEnergy(TreeReader *fTR, TString mode, bool isbarrel);

  void FillLead(int index);
  void FillTrail(int index);

  //  double etaTransformation(float EtaParticle, float Zvertex);
  double phiNorm(float phi);
  isolations_struct RandomConeIsolation(TreeReader *fTR, int phoqi, TString mod="");
  isolations_struct PFConeIsolation(TreeReader *fTR, int phoqi);
  int FindPFCandType(int id);

  int CountChargedHadronsInCone(TreeReader *fTR, int phoqi, std::vector<int> removals = std::vector<int>(), bool skipvetocones=false);
  std::vector<int> NChargedHadronsInConeSelection(TreeReader *fTR, std::vector<int> passing, int minimum=0, int maximum=9999);

  std::vector<int> MuonSelection(TreeReader *fTR, std::vector<int> passing);
  bool DiMuonFromZSelection(TreeReader *fTR, std::vector<int> &passing);
  float PFPhotonIsolationAroundMuon(int muqi, int *counter, std::vector<float> *energies = NULL, std::vector<float> *ets = NULL,  std::vector<float> *detas = NULL, std::vector<float> *dphis = NULL);
  void FillMuonInfo(int index);

  float DeltaPhiSigned(float phi1, float phi2);

  TRandom3 *randomgen;

  float scarea[100];
  float scareaSF[100];
  
  std::string fDataType_;
  bool isdata;

  TH1F *histo_PFPhotonDepositAroundImpingingTrack;

  Float_t AddWeight;

  int impinging_track_pfcand[100];

  TTree* OutputTree[12];

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

  Int_t pholead_Nchargedhadronsincone, photrail_Nchargedhadronsincone;

  Float_t pholead_scarea, photrail_scarea;
  Float_t pholead_scareaSF, photrail_scareaSF;

  Int_t pholead_Npfcandphotonincone;
  Int_t pholead_Npfcandchargedincone;
  Int_t pholead_Npfcandneutralincone;
  Int_t photrail_Npfcandphotonincone;
  Int_t photrail_Npfcandchargedincone;
  Int_t photrail_Npfcandneutralincone;

  Float_t pholead_photonpfcandenergies[global_size_pfcandarrays];
  Float_t pholead_chargedpfcandenergies[global_size_pfcandarrays];
  Float_t pholead_neutralpfcandenergies[global_size_pfcandarrays];
  Float_t photrail_photonpfcandenergies[global_size_pfcandarrays];
  Float_t photrail_chargedpfcandenergies[global_size_pfcandarrays];
  Float_t photrail_neutralpfcandenergies[global_size_pfcandarrays];

  Float_t pholead_photonpfcandets[global_size_pfcandarrays];
  Float_t pholead_chargedpfcandets[global_size_pfcandarrays];
  Float_t pholead_neutralpfcandets[global_size_pfcandarrays];
  Float_t photrail_photonpfcandets[global_size_pfcandarrays];
  Float_t photrail_chargedpfcandets[global_size_pfcandarrays];
  Float_t photrail_neutralpfcandets[global_size_pfcandarrays];

  Float_t pholead_photonpfcanddetas[global_size_pfcandarrays];
  Float_t pholead_chargedpfcanddetas[global_size_pfcandarrays];
  Float_t pholead_neutralpfcanddetas[global_size_pfcandarrays];
  Float_t photrail_photonpfcanddetas[global_size_pfcandarrays];
  Float_t photrail_chargedpfcanddetas[global_size_pfcandarrays];
  Float_t photrail_neutralpfcanddetas[global_size_pfcandarrays];

  Float_t pholead_photonpfcanddphis[global_size_pfcandarrays];
  Float_t pholead_chargedpfcanddphis[global_size_pfcandarrays];
  Float_t pholead_neutralpfcanddphis[global_size_pfcandarrays];
  Float_t photrail_photonpfcanddphis[global_size_pfcandarrays];
  Float_t photrail_chargedpfcanddphis[global_size_pfcandarrays];
  Float_t photrail_neutralpfcanddphis[global_size_pfcandarrays];

  TH1F *fHNumPU;
  TH1F *fHNumPU_noweight;
  TH1F *fHNumVtx;

};
#endif
