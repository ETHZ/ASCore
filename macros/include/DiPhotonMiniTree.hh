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

class DiPhotonMiniTree : public UserAnalysisBase{
public:
  DiPhotonMiniTree(TreeReader *tr = NULL, std::string dataType="data", double aw=-999, double* _kfac=NULL);
  virtual ~DiPhotonMiniTree();

  void Begin();
  void Analyze();
  void End();

private:

  double* kfactors;

  std::vector<int> PhotonSelection(TreeReader *fTR);
  bool TriggerSelection();

  EnergyCorrection *phocorr;
  TLorentzVector CorrPhoton(TreeReader *fTR, int i, int mode);
	
  std::string fDataType_;
  bool isdata;

  double AddWeight;

  TTree* OutputTree;

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
  


  Float_t photrail_Pho_Cone04PhotonIso_dR0_dEta0_pt0,    pholead_Pho_Cone04PhotonIso_dR0_dEta0_pt0;
  Float_t photrail_Pho_Cone04PhotonIso_dR0_dEta0_pt5,    pholead_Pho_Cone04PhotonIso_dR0_dEta0_pt5;
  Float_t photrail_Pho_Cone04PhotonIso_dR8_dEta0_pt0,    pholead_Pho_Cone04PhotonIso_dR8_dEta0_pt0;
  Float_t photrail_Pho_Cone04PhotonIso_dR8_dEta0_pt5,    pholead_Pho_Cone04PhotonIso_dR8_dEta0_pt5;
  Float_t photrail_Pho_Cone01PhotonIso_dR045EB070EE_dEta015_pt08EB1EE_mvVtx,    pholead_Pho_Cone01PhotonIso_dR045EB070EE_dEta015_pt08EB1EE_mvVtx;
  Float_t photrail_Pho_Cone02PhotonIso_dR045EB070EE_dEta015_pt08EB1EE_mvVtx,    pholead_Pho_Cone02PhotonIso_dR045EB070EE_dEta015_pt08EB1EE_mvVtx;
  Float_t photrail_Pho_Cone03PhotonIso_dR045EB070EE_dEta015_pt08EB1EE_mvVtx,    pholead_Pho_Cone03PhotonIso_dR045EB070EE_dEta015_pt08EB1EE_mvVtx;
  Float_t photrail_Pho_Cone04PhotonIso_dR045EB070EE_dEta015_pt08EB1EE_mvVtx,    pholead_Pho_Cone04PhotonIso_dR045EB070EE_dEta015_pt08EB1EE_mvVtx;
  Float_t photrail_Pho_Cone04NeutralHadronIso_dR0_dEta0_pt0,    pholead_Pho_Cone04NeutralHadronIso_dR0_dEta0_pt0;
  Float_t photrail_Pho_Cone04NeutralHadronIso_dR0_dEta0_pt5,    pholead_Pho_Cone04NeutralHadronIso_dR0_dEta0_pt5;
  Float_t photrail_Pho_Cone04NeutralHadronIso_dR0_dEta0_pt0_nocracks,    pholead_Pho_Cone04NeutralHadronIso_dR0_dEta0_pt0_nocracks;
  Float_t photrail_Pho_Cone04NeutralHadronIso_dR0_dEta0_pt5_nocracks,    pholead_Pho_Cone04NeutralHadronIso_dR0_dEta0_pt5_nocracks;
  Float_t photrail_Pho_Cone04NeutralHadronIso_dR7_dEta0_pt0,    pholead_Pho_Cone04NeutralHadronIso_dR7_dEta0_pt0;
  Float_t photrail_Pho_Cone04NeutralHadronIso_dR7_dEta0_pt5,    pholead_Pho_Cone04NeutralHadronIso_dR7_dEta0_pt5;
  Float_t photrail_Pho_Cone01NeutralHadronIso_dR0_dEta0_pt0_mvVtx,    pholead_Pho_Cone01NeutralHadronIso_dR0_dEta0_pt0_mvVtx;
  Float_t photrail_Pho_Cone02NeutralHadronIso_dR0_dEta0_pt0_mvVtx,    pholead_Pho_Cone02NeutralHadronIso_dR0_dEta0_pt0_mvVtx;
  Float_t photrail_Pho_Cone03NeutralHadronIso_dR0_dEta0_pt0_mvVtx,    pholead_Pho_Cone03NeutralHadronIso_dR0_dEta0_pt0_mvVtx;
  Float_t photrail_Pho_Cone04NeutralHadronIso_dR0_dEta0_pt0_mvVtx,    pholead_Pho_Cone04NeutralHadronIso_dR0_dEta0_pt0_mvVtx;
  Float_t photrail_Pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_dz0_old,    pholead_Pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_dz0_old;
  Float_t photrail_Pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_PFnoPU_old,    pholead_Pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_PFnoPU_old;
  Float_t photrail_Pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_dz0_old,    pholead_Pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_dz0_old;
  Float_t photrail_Pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_PFnoPU_old,    pholead_Pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_PFnoPU_old;
  Float_t photrail_Pho_Cone01ChargedHadronIso_dR0_dEta0_pt0_dz0,    pholead_Pho_Cone01ChargedHadronIso_dR0_dEta0_pt0_dz0;
  Float_t photrail_Pho_Cone01ChargedHadronIso_dR0_dEta0_pt0_dz1_dxy01,    pholead_Pho_Cone01ChargedHadronIso_dR0_dEta0_pt0_dz1_dxy01;
  Float_t photrail_Pho_Cone01ChargedHadronIso_dR0_dEta0_pt0_PFnoPU,    pholead_Pho_Cone01ChargedHadronIso_dR0_dEta0_pt0_PFnoPU;
  Float_t photrail_Pho_Cone01ChargedHadronIso_dR015_dEta0_pt0_dz0,    pholead_Pho_Cone01ChargedHadronIso_dR015_dEta0_pt0_dz0;
  Float_t photrail_Pho_Cone01ChargedHadronIso_dR015_dEta0_pt0_dz1_dxy01,    pholead_Pho_Cone01ChargedHadronIso_dR015_dEta0_pt0_dz1_dxy01;
  Float_t photrail_Pho_Cone01ChargedHadronIso_dR015_dEta0_pt0_PFnoPU,    pholead_Pho_Cone01ChargedHadronIso_dR015_dEta0_pt0_PFnoPU;
  Float_t photrail_Pho_Cone02ChargedHadronIso_dR0_dEta0_pt0_dz0,    pholead_Pho_Cone02ChargedHadronIso_dR0_dEta0_pt0_dz0;
  Float_t photrail_Pho_Cone02ChargedHadronIso_dR0_dEta0_pt0_dz1_dxy01,    pholead_Pho_Cone02ChargedHadronIso_dR0_dEta0_pt0_dz1_dxy01;
  Float_t photrail_Pho_Cone02ChargedHadronIso_dR0_dEta0_pt0_PFnoPU,    pholead_Pho_Cone02ChargedHadronIso_dR0_dEta0_pt0_PFnoPU;
  Float_t photrail_Pho_Cone02ChargedHadronIso_dR015_dEta0_pt0_dz0,    pholead_Pho_Cone02ChargedHadronIso_dR015_dEta0_pt0_dz0;
  Float_t photrail_Pho_Cone02ChargedHadronIso_dR015_dEta0_pt0_dz1_dxy01,    pholead_Pho_Cone02ChargedHadronIso_dR015_dEta0_pt0_dz1_dxy01;
  Float_t photrail_Pho_Cone02ChargedHadronIso_dR015_dEta0_pt0_PFnoPU,    pholead_Pho_Cone02ChargedHadronIso_dR015_dEta0_pt0_PFnoPU;
  Float_t photrail_Pho_Cone03ChargedHadronIso_dR0_dEta0_pt0_dz0,    pholead_Pho_Cone03ChargedHadronIso_dR0_dEta0_pt0_dz0;
  Float_t photrail_Pho_Cone03ChargedHadronIso_dR0_dEta0_pt0_dz1_dxy01,    pholead_Pho_Cone03ChargedHadronIso_dR0_dEta0_pt0_dz1_dxy01;
  Float_t photrail_Pho_Cone03ChargedHadronIso_dR0_dEta0_pt0_PFnoPU,    pholead_Pho_Cone03ChargedHadronIso_dR0_dEta0_pt0_PFnoPU;
  Float_t photrail_Pho_Cone03ChargedHadronIso_dR015_dEta0_pt0_dz0,    pholead_Pho_Cone03ChargedHadronIso_dR015_dEta0_pt0_dz0;
  Float_t photrail_Pho_Cone03ChargedHadronIso_dR015_dEta0_pt0_dz1_dxy01,    pholead_Pho_Cone03ChargedHadronIso_dR015_dEta0_pt0_dz1_dxy01;
  Float_t photrail_Pho_Cone03ChargedHadronIso_dR015_dEta0_pt0_PFnoPU,    pholead_Pho_Cone03ChargedHadronIso_dR015_dEta0_pt0_PFnoPU;
  Float_t photrail_Pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_dz0,    pholead_Pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_dz0;
  Float_t photrail_Pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_dz1_dxy01,    pholead_Pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_dz1_dxy01;
  Float_t photrail_Pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_PFnoPU,    pholead_Pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_PFnoPU;
  Float_t photrail_Pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_dz0,    pholead_Pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_dz0;
  Float_t photrail_Pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_dz1_dxy01,    pholead_Pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_dz1_dxy01;
  Float_t photrail_Pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_PFnoPU,    pholead_Pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_PFnoPU;

 

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


  Int_t pholead_PhoMCmatchindex, photrail_PhoMCmatchindex;
  Int_t pholead_PhoMCmatchexitcode, photrail_PhoMCmatchexitcode;




  TH1F *fHNumPU;
  TH1F *fHNumVtx;

};
#endif
