#ifndef ZeeMiniTree_hh
#define ZeeMiniTree_hh


#include <cmath>
#include <string>
#include <iostream>
#include <fstream>

#include <TH1D.h>

#include "base/TreeReader.hh"
#include "base/UserAnalysisBase.hh"

#include "TLorentzVector.h"

#include "EnergyCorrection.hh"

class ZeeMiniTree : public UserAnalysisBase{
public:
  ZeeMiniTree(TreeReader *tr = NULL, std::string dataType="data");
	virtual ~ZeeMiniTree();

	void Begin();
	void Analyze();
	void End();

private:

  EnergyCorrection *elecorr;
  TLorentzVector CorrElectron(TreeReader *fTR, int i, int mode);
	
  std::string fDataType_;
  bool isdata;

  TFile* fMiniTree;
  TTree* OutputTree;

  Float_t event_weight;
  Float_t event_rho;
  Int_t event_nPU;
  Int_t event_nRecVtx;
  
  Int_t diel_cat;
  Float_t diel_mee_electron;
  Float_t diel_mee_SCdefault;
  Float_t diel_mee_newCorrNoCrack;
  Float_t diel_mee_newCorr;
  Float_t diel_mee_newCorrLocal;

  Float_t ellead_eta, eltrail_eta;
  Float_t ellead_px, eltrail_px;
  Float_t ellead_py, eltrail_py;
  Float_t ellead_pt, eltrail_pt;
  Float_t ellead_pz, eltrail_pz;
  Float_t ellead_energy, eltrail_energy;
  Float_t ellead_SCeta, eltrail_SCeta;
  Float_t ellead_fbrem, eltrail_fbrem;
  Float_t ellead_Pin, eltrail_Pin;

  Float_t ellead_SCenergyCetaCorr, eltrail_SCenergyCetaCorr;
  Float_t ellead_energySCdefault, eltrail_energySCdefault;
  Float_t ellead_energyNewCorrNoCrack, eltrail_energyNewCorrNoCrack;
  Float_t ellead_energyNewCorr, eltrail_energyNewCorr;
  Float_t ellead_energyNewCorrLocal, eltrail_energyNewCorrLocal;

  Float_t ellead_r9, eltrail_r9;
  Float_t ellead_sieie, eltrail_sieie;
  Float_t ellead_hoe, eltrail_hoe;
  Float_t ellead_brem, eltrail_brem;
  Float_t ellead_sigmaPhi, eltrail_sigmaPhi;
  Float_t ellead_sigmaEta, eltrail_sigmaEta;
  

  TH1F *fHNumPU;
  TH1F *fHNumVtx;

};
#endif
