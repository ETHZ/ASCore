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
  DiPhotonMiniTree(TreeReader *tr = NULL, std::string dataType="data");
	virtual ~DiPhotonMiniTree();

	void Begin();
	void Analyze();
	void End();

private:

  std::vector<int> PhotonSelection(TreeReader *fTR);
  bool PhotonID_EGM_10_006_Loose(TreeReader *fTR, int i);
  bool PhotonID_EGM_10_006_Loose_SigmaIetaIeta_Relaxed(TreeReader *fTR, int i);

  bool TriggerSelection();

  EnergyCorrection *phocorr;
  TLorentzVector CorrPhoton(TreeReader *fTR, int i, int mode);
	
  std::string fDataType_;
  bool isdata;

  TFile* fMiniTree;
  TTree* OutputTree;

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
  
  Float_t pholead_energySCdefault, photrail_energySCdefault;
  Float_t pholead_energyNewCorr, photrail_energyNewCorr;
  Float_t pholead_energyNewCorrLocal, photrail_energyNewCorrLocal;

  Float_t pholead_r9, photrail_r9;
  Float_t pholead_sieie, photrail_sieie;
  Float_t pholead_hoe, photrail_hoe;
  Float_t pholead_brem, photrail_brem;
  Float_t pholead_sigmaPhi, photrail_sigmaPhi;
  Float_t pholead_sigmaEta, photrail_sigmaEta;
  

  TH1F *fHNumPU;
  TH1F *fHNumVtx;

};
#endif
