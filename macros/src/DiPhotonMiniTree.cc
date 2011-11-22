#include "helper/Utilities.hh"
#include "DiPhotonMiniTree.hh"

#include "DiPhotonPurity.hh"


#include <iostream>
#include <fstream>
#include <assert.h>

using namespace std;

DiPhotonMiniTree::DiPhotonMiniTree(TreeReader *tr, std::string dataType, double aw) : UserAnalysisBase(tr), fDataType_(dataType), AddWeight(aw){
  Util::SetStyle();	
  if (fDataType_ == "mc") isdata=false;
  else if (fDataType_ == "data") isdata=true; 
  else {
    std::cout << "wrong data type" << std::endl;
    assert(1==0);
  }
  phocorr = new EnergyCorrection("photons");
}

DiPhotonMiniTree::~DiPhotonMiniTree(){
  delete phocorr;
}

void DiPhotonMiniTree::Begin(){

  cout << "Begin" << endl;

  // Define the output file of histograms
  const char* filename = "MiniTree_Diphoton.root";
  fMiniTree = new TFile(fOutputDir + TString(filename), "RECREATE");

  fMiniTree->cd();
  OutputTree = new TTree("Tree","Tree");

 OutputTree->Branch("event_luminormfactor",&event_luminormfactor,"event_luminormfactor/F");

  OutputTree->Branch("event_weight",&event_weight,"event_weight/F");
  OutputTree->Branch("event_rho",&event_rho,"event_rho/F");
  OutputTree->Branch("event_nPU",&event_nPU,"event_nPU/I");
  OutputTree->Branch("event_nRecVtx",&event_nRecVtx,"event_nRecVtx/I");

  OutputTree->Branch("dipho_mgg_photon",&dipho_mgg_photon,"dipho_mgg_photon/F");
  OutputTree->Branch("dipho_mgg_newCorr",&dipho_mgg_newCorr,"dipho_mgg_newCorr/F");
  OutputTree->Branch("dipho_mgg_newCorrLocal",&dipho_mgg_newCorrLocal,"dipho_mgg_newCorrLocal/F");

  OutputTree->Branch("pholead_eta",&pholead_eta,"pholead_eta/F");
  OutputTree->Branch("photrail_eta",&photrail_eta,"photrail_eta/F");
  OutputTree->Branch("pholead_px",&pholead_px,"pholead_px/F");
  OutputTree->Branch("photrail_px",&photrail_px,"photrail_px/F");
  OutputTree->Branch("pholead_py",&pholead_py,"pholead_py/F");
  OutputTree->Branch("photrail_py",&photrail_py,"photrail_py/F");
  OutputTree->Branch("pholead_pt",&pholead_pt,"pholead_pt/F");
  OutputTree->Branch("photrail_pt",&photrail_pt,"photrail_pt/F");
  OutputTree->Branch("pholead_pz",&pholead_pz,"pholead_pz/F");
  OutputTree->Branch("photrail_pz",&photrail_pz,"photrail_pz/F");
  OutputTree->Branch("pholead_energy",&pholead_energy,"pholead_energy/F");
  OutputTree->Branch("photrail_energy",&photrail_energy,"photrail_energy/F");
  OutputTree->Branch("pholead_energySCdefault",&pholead_energySCdefault,"pholead_energySCdefault/F");
  OutputTree->Branch("photrail_energySCdefault",&photrail_energySCdefault,"photrail_energySCdefault/F");
  OutputTree->Branch("pholead_energyNewCorr",&pholead_energyNewCorr,"pholead_energyNewCorr/F");
  OutputTree->Branch("photrail_energyNewCorr",&photrail_energyNewCorr,"photrail_energyNewCorr/F");
  OutputTree->Branch("pholead_energyNewCorrLocal",&pholead_energyNewCorrLocal,"pholead_energyNewCorrLocal/F");
  OutputTree->Branch("photrail_energyNewCorrLocal",&photrail_energyNewCorrLocal,"photrail_energyNewCorrLocal/F");
  OutputTree->Branch("pholead_SCeta",&pholead_SCeta,"pholead_SCeta/F");
  OutputTree->Branch("photrail_SCeta",&photrail_SCeta,"photrail_SCeta/F");

  OutputTree->Branch("pholead_r9",&pholead_r9,"pholead_r9/F");
  OutputTree->Branch("photrail_r9",&photrail_r9,"photrail_r9/F");
  OutputTree->Branch("pholead_sieie",&pholead_sieie,"pholead_sieie/F");
  OutputTree->Branch("photrail_sieie",&photrail_sieie,"photrail_sieie/F");
  OutputTree->Branch("pholead_hoe",&pholead_hoe,"pholead_hoe/F");
  OutputTree->Branch("photrail_hoe",&photrail_hoe,"photrail_hoe/F");
  OutputTree->Branch("pholead_brem",&pholead_brem,"pholead_brem/F");
  OutputTree->Branch("photrail_brem",&photrail_brem,"photrail_brem/F");
  OutputTree->Branch("pholead_sigmaPhi",&pholead_sigmaPhi,"pholead_sigmaPhi/F");
  OutputTree->Branch("photrail_sigmaPhi",&photrail_sigmaPhi,"photrail_sigmaPhi/F");
  OutputTree->Branch("pholead_sigmaEta",&pholead_sigmaEta,"pholead_sigmaEta/F");
  OutputTree->Branch("photrail_sigmaEta",&photrail_sigmaEta,"photrail_sigmaEta/F");


  OutputTree->Branch("pholead_Pho_Cone04PhotonIso_dR0_dEta0_pt5",&pholead_Pho_Cone04PhotonIso_dR0_dEta0_pt5,"pholead_Pho_Cone04PhotonIso_dR0_dEta0_pt5/F");
  OutputTree->Branch("pholead_Pho_Cone04PhotonIso_dR8_dEta0_pt0",&pholead_Pho_Cone04PhotonIso_dR8_dEta0_pt0,"pholead_Pho_Cone04PhotonIso_dR8_dEta0_pt0/F");
  OutputTree->Branch("pholead_Pho_Cone04PhotonIso_dR8_dEta0_pt5",&pholead_Pho_Cone04PhotonIso_dR8_dEta0_pt5,"pholead_Pho_Cone04PhotonIso_dR8_dEta0_pt5/F");
  OutputTree->Branch("pholead_Pho_Cone01PhotonIso_dR045EB070EE_dEta015_pt08EB1EE_mvVtx",&pholead_Pho_Cone01PhotonIso_dR045EB070EE_dEta015_pt08EB1EE_mvVtx,"pholead_Pho_Cone01PhotonIso_dR045EB070EE_dEta015_pt08EB1EE_mvVtx/F");
  OutputTree->Branch("pholead_Pho_Cone02PhotonIso_dR045EB070EE_dEta015_pt08EB1EE_mvVtx",&pholead_Pho_Cone02PhotonIso_dR045EB070EE_dEta015_pt08EB1EE_mvVtx,"pholead_Pho_Cone02PhotonIso_dR045EB070EE_dEta015_pt08EB1EE_mvVtx/F");
  OutputTree->Branch("pholead_Pho_Cone03PhotonIso_dR045EB070EE_dEta015_pt08EB1EE_mvVtx",&pholead_Pho_Cone03PhotonIso_dR045EB070EE_dEta015_pt08EB1EE_mvVtx,"pholead_Pho_Cone03PhotonIso_dR045EB070EE_dEta015_pt08EB1EE_mvVtx/F");
  OutputTree->Branch("pholead_Pho_Cone04PhotonIso_dR045EB070EE_dEta015_pt08EB1EE_mvVtx",&pholead_Pho_Cone04PhotonIso_dR045EB070EE_dEta015_pt08EB1EE_mvVtx,"pholead_Pho_Cone04PhotonIso_dR045EB070EE_dEta015_pt08EB1EE_mvVtx/F");
  OutputTree->Branch("pholead_Pho_Cone04NeutralHadronIso_dR0_dEta0_pt0",&pholead_Pho_Cone04NeutralHadronIso_dR0_dEta0_pt0,"pholead_Pho_Cone04NeutralHadronIso_dR0_dEta0_pt0/F");
  OutputTree->Branch("pholead_Pho_Cone04NeutralHadronIso_dR0_dEta0_pt5",&pholead_Pho_Cone04NeutralHadronIso_dR0_dEta0_pt5,"pholead_Pho_Cone04NeutralHadronIso_dR0_dEta0_pt5/F");
  OutputTree->Branch("pholead_Pho_Cone04NeutralHadronIso_dR0_dEta0_pt0_nocracks",&pholead_Pho_Cone04NeutralHadronIso_dR0_dEta0_pt0_nocracks,"pholead_Pho_Cone04NeutralHadronIso_dR0_dEta0_pt0_nocracks/F");
  OutputTree->Branch("pholead_Pho_Cone04NeutralHadronIso_dR0_dEta0_pt5_nocracks",&pholead_Pho_Cone04NeutralHadronIso_dR0_dEta0_pt5_nocracks,"pholead_Pho_Cone04NeutralHadronIso_dR0_dEta0_pt5_nocracks/F");
  OutputTree->Branch("pholead_Pho_Cone04NeutralHadronIso_dR7_dEta0_pt0",&pholead_Pho_Cone04NeutralHadronIso_dR7_dEta0_pt0,"pholead_Pho_Cone04NeutralHadronIso_dR7_dEta0_pt0/F");
  OutputTree->Branch("pholead_Pho_Cone04NeutralHadronIso_dR7_dEta0_pt5",&pholead_Pho_Cone04NeutralHadronIso_dR7_dEta0_pt5,"pholead_Pho_Cone04NeutralHadronIso_dR7_dEta0_pt5/F");
  OutputTree->Branch("pholead_Pho_Cone01NeutralHadronIso_dR0_dEta0_pt0_mvVtx",&pholead_Pho_Cone01NeutralHadronIso_dR0_dEta0_pt0_mvVtx,"pholead_Pho_Cone01NeutralHadronIso_dR0_dEta0_pt0_mvVtx/F");
  OutputTree->Branch("pholead_Pho_Cone02NeutralHadronIso_dR0_dEta0_pt0_mvVtx",&pholead_Pho_Cone02NeutralHadronIso_dR0_dEta0_pt0_mvVtx,"pholead_Pho_Cone02NeutralHadronIso_dR0_dEta0_pt0_mvVtx/F");
  OutputTree->Branch("pholead_Pho_Cone03NeutralHadronIso_dR0_dEta0_pt0_mvVtx",&pholead_Pho_Cone03NeutralHadronIso_dR0_dEta0_pt0_mvVtx,"pholead_Pho_Cone03NeutralHadronIso_dR0_dEta0_pt0_mvVtx/F");
  OutputTree->Branch("pholead_Pho_Cone04NeutralHadronIso_dR0_dEta0_pt0_mvVtx",&pholead_Pho_Cone04NeutralHadronIso_dR0_dEta0_pt0_mvVtx,"pholead_Pho_Cone04NeutralHadronIso_dR0_dEta0_pt0_mvVtx/F");
  OutputTree->Branch("pholead_Pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_dz0_old",&pholead_Pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_dz0_old,"pholead_Pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_dz0_old/F");
  OutputTree->Branch("pholead_Pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_PFnoPU_old",&pholead_Pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_PFnoPU_old,"pholead_Pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_PFnoPU_old/F");
  OutputTree->Branch("pholead_Pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_dz0_old",&pholead_Pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_dz0_old,"pholead_Pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_dz0_old/F");
  OutputTree->Branch("pholead_Pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_PFnoPU_old",&pholead_Pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_PFnoPU_old,"pholead_Pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_PFnoPU_old/F");
  OutputTree->Branch("pholead_Pho_Cone01ChargedHadronIso_dR0_dEta0_pt0_dz0",&pholead_Pho_Cone01ChargedHadronIso_dR0_dEta0_pt0_dz0,"pholead_Pho_Cone01ChargedHadronIso_dR0_dEta0_pt0_dz0/F");
  OutputTree->Branch("pholead_Pho_Cone01ChargedHadronIso_dR0_dEta0_pt0_dz1_dxy01",&pholead_Pho_Cone01ChargedHadronIso_dR0_dEta0_pt0_dz1_dxy01,"pholead_Pho_Cone01ChargedHadronIso_dR0_dEta0_pt0_dz1_dxy01/F");
  OutputTree->Branch("pholead_Pho_Cone01ChargedHadronIso_dR0_dEta0_pt0_PFnoPU",&pholead_Pho_Cone01ChargedHadronIso_dR0_dEta0_pt0_PFnoPU,"pholead_Pho_Cone01ChargedHadronIso_dR0_dEta0_pt0_PFnoPU/F");
  OutputTree->Branch("pholead_Pho_Cone01ChargedHadronIso_dR015_dEta0_pt0_dz0",&pholead_Pho_Cone01ChargedHadronIso_dR015_dEta0_pt0_dz0,"pholead_Pho_Cone01ChargedHadronIso_dR015_dEta0_pt0_dz0/F");
  OutputTree->Branch("pholead_Pho_Cone01ChargedHadronIso_dR015_dEta0_pt0_dz1_dxy01",&pholead_Pho_Cone01ChargedHadronIso_dR015_dEta0_pt0_dz1_dxy01,"pholead_Pho_Cone01ChargedHadronIso_dR015_dEta0_pt0_dz1_dxy01/F");
  OutputTree->Branch("pholead_Pho_Cone01ChargedHadronIso_dR015_dEta0_pt0_PFnoPU",&pholead_Pho_Cone01ChargedHadronIso_dR015_dEta0_pt0_PFnoPU,"pholead_Pho_Cone01ChargedHadronIso_dR015_dEta0_pt0_PFnoPU/F");
  OutputTree->Branch("pholead_Pho_Cone02ChargedHadronIso_dR0_dEta0_pt0_dz0",&pholead_Pho_Cone02ChargedHadronIso_dR0_dEta0_pt0_dz0,"pholead_Pho_Cone02ChargedHadronIso_dR0_dEta0_pt0_dz0/F");
  OutputTree->Branch("pholead_Pho_Cone02ChargedHadronIso_dR0_dEta0_pt0_dz1_dxy01",&pholead_Pho_Cone02ChargedHadronIso_dR0_dEta0_pt0_dz1_dxy01,"pholead_Pho_Cone02ChargedHadronIso_dR0_dEta0_pt0_dz1_dxy01/F");
  OutputTree->Branch("pholead_Pho_Cone02ChargedHadronIso_dR0_dEta0_pt0_PFnoPU",&pholead_Pho_Cone02ChargedHadronIso_dR0_dEta0_pt0_PFnoPU,"pholead_Pho_Cone02ChargedHadronIso_dR0_dEta0_pt0_PFnoPU/F");
  OutputTree->Branch("pholead_Pho_Cone02ChargedHadronIso_dR015_dEta0_pt0_dz0",&pholead_Pho_Cone02ChargedHadronIso_dR015_dEta0_pt0_dz0,"pholead_Pho_Cone02ChargedHadronIso_dR015_dEta0_pt0_dz0/F");
  OutputTree->Branch("pholead_Pho_Cone02ChargedHadronIso_dR015_dEta0_pt0_dz1_dxy01",&pholead_Pho_Cone02ChargedHadronIso_dR015_dEta0_pt0_dz1_dxy01,"pholead_Pho_Cone02ChargedHadronIso_dR015_dEta0_pt0_dz1_dxy01/F");
  OutputTree->Branch("pholead_Pho_Cone02ChargedHadronIso_dR015_dEta0_pt0_PFnoPU",&pholead_Pho_Cone02ChargedHadronIso_dR015_dEta0_pt0_PFnoPU,"pholead_Pho_Cone02ChargedHadronIso_dR015_dEta0_pt0_PFnoPU/F");
  OutputTree->Branch("pholead_Pho_Cone03ChargedHadronIso_dR0_dEta0_pt0_dz0",&pholead_Pho_Cone03ChargedHadronIso_dR0_dEta0_pt0_dz0,"pholead_Pho_Cone03ChargedHadronIso_dR0_dEta0_pt0_dz0/F");
  OutputTree->Branch("pholead_Pho_Cone03ChargedHadronIso_dR0_dEta0_pt0_dz1_dxy01",&pholead_Pho_Cone03ChargedHadronIso_dR0_dEta0_pt0_dz1_dxy01,"pholead_Pho_Cone03ChargedHadronIso_dR0_dEta0_pt0_dz1_dxy01/F");
  OutputTree->Branch("pholead_Pho_Cone03ChargedHadronIso_dR0_dEta0_pt0_PFnoPU",&pholead_Pho_Cone03ChargedHadronIso_dR0_dEta0_pt0_PFnoPU,"pholead_Pho_Cone03ChargedHadronIso_dR0_dEta0_pt0_PFnoPU/F");
  OutputTree->Branch("pholead_Pho_Cone03ChargedHadronIso_dR015_dEta0_pt0_dz0",&pholead_Pho_Cone03ChargedHadronIso_dR015_dEta0_pt0_dz0,"pholead_Pho_Cone03ChargedHadronIso_dR015_dEta0_pt0_dz0/F");
  OutputTree->Branch("pholead_Pho_Cone03ChargedHadronIso_dR015_dEta0_pt0_dz1_dxy01",&pholead_Pho_Cone03ChargedHadronIso_dR015_dEta0_pt0_dz1_dxy01,"pholead_Pho_Cone03ChargedHadronIso_dR015_dEta0_pt0_dz1_dxy01/F");
  OutputTree->Branch("pholead_Pho_Cone03ChargedHadronIso_dR015_dEta0_pt0_PFnoPU",&pholead_Pho_Cone03ChargedHadronIso_dR015_dEta0_pt0_PFnoPU,"pholead_Pho_Cone03ChargedHadronIso_dR015_dEta0_pt0_PFnoPU/F");
  OutputTree->Branch("pholead_Pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_dz0",&pholead_Pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_dz0,"pholead_Pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_dz0/F");
  OutputTree->Branch("pholead_Pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_dz1_dxy01",&pholead_Pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_dz1_dxy01,"pholead_Pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_dz1_dxy01/F");
  OutputTree->Branch("pholead_Pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_PFnoPU",&pholead_Pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_PFnoPU,"pholead_Pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_PFnoPU/F");
  OutputTree->Branch("pholead_Pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_dz0",&pholead_Pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_dz0,"pholead_Pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_dz0/F");
  OutputTree->Branch("pholead_Pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_dz1_dxy01",&pholead_Pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_dz1_dxy01,"pholead_Pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_dz1_dxy01/F");
  OutputTree->Branch("pholead_Pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_PFnoPU",&pholead_Pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_PFnoPU,"pholead_Pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_PFnoPU/F");


  OutputTree->Branch("photrail_Pho_Cone04PhotonIso_dR0_dEta0_pt5",&photrail_Pho_Cone04PhotonIso_dR0_dEta0_pt5,"photrail_Pho_Cone04PhotonIso_dR0_dEta0_pt5/F");
  OutputTree->Branch("photrail_Pho_Cone04PhotonIso_dR8_dEta0_pt0",&photrail_Pho_Cone04PhotonIso_dR8_dEta0_pt0,"photrail_Pho_Cone04PhotonIso_dR8_dEta0_pt0/F");
  OutputTree->Branch("photrail_Pho_Cone04PhotonIso_dR8_dEta0_pt5",&photrail_Pho_Cone04PhotonIso_dR8_dEta0_pt5,"photrail_Pho_Cone04PhotonIso_dR8_dEta0_pt5/F");
  OutputTree->Branch("photrail_Pho_Cone01PhotonIso_dR045EB070EE_dEta015_pt08EB1EE_mvVtx",&photrail_Pho_Cone01PhotonIso_dR045EB070EE_dEta015_pt08EB1EE_mvVtx,"photrail_Pho_Cone01PhotonIso_dR045EB070EE_dEta015_pt08EB1EE_mvVtx/F");
  OutputTree->Branch("photrail_Pho_Cone02PhotonIso_dR045EB070EE_dEta015_pt08EB1EE_mvVtx",&photrail_Pho_Cone02PhotonIso_dR045EB070EE_dEta015_pt08EB1EE_mvVtx,"photrail_Pho_Cone02PhotonIso_dR045EB070EE_dEta015_pt08EB1EE_mvVtx/F");
  OutputTree->Branch("photrail_Pho_Cone03PhotonIso_dR045EB070EE_dEta015_pt08EB1EE_mvVtx",&photrail_Pho_Cone03PhotonIso_dR045EB070EE_dEta015_pt08EB1EE_mvVtx,"photrail_Pho_Cone03PhotonIso_dR045EB070EE_dEta015_pt08EB1EE_mvVtx/F");
  OutputTree->Branch("photrail_Pho_Cone04PhotonIso_dR045EB070EE_dEta015_pt08EB1EE_mvVtx",&photrail_Pho_Cone04PhotonIso_dR045EB070EE_dEta015_pt08EB1EE_mvVtx,"photrail_Pho_Cone04PhotonIso_dR045EB070EE_dEta015_pt08EB1EE_mvVtx/F");
  OutputTree->Branch("photrail_Pho_Cone04NeutralHadronIso_dR0_dEta0_pt0",&photrail_Pho_Cone04NeutralHadronIso_dR0_dEta0_pt0,"photrail_Pho_Cone04NeutralHadronIso_dR0_dEta0_pt0/F");
  OutputTree->Branch("photrail_Pho_Cone04NeutralHadronIso_dR0_dEta0_pt5",&photrail_Pho_Cone04NeutralHadronIso_dR0_dEta0_pt5,"photrail_Pho_Cone04NeutralHadronIso_dR0_dEta0_pt5/F");
  OutputTree->Branch("photrail_Pho_Cone04NeutralHadronIso_dR0_dEta0_pt0_nocracks",&photrail_Pho_Cone04NeutralHadronIso_dR0_dEta0_pt0_nocracks,"photrail_Pho_Cone04NeutralHadronIso_dR0_dEta0_pt0_nocracks/F");
  OutputTree->Branch("photrail_Pho_Cone04NeutralHadronIso_dR0_dEta0_pt5_nocracks",&photrail_Pho_Cone04NeutralHadronIso_dR0_dEta0_pt5_nocracks,"photrail_Pho_Cone04NeutralHadronIso_dR0_dEta0_pt5_nocracks/F");
  OutputTree->Branch("photrail_Pho_Cone04NeutralHadronIso_dR7_dEta0_pt0",&photrail_Pho_Cone04NeutralHadronIso_dR7_dEta0_pt0,"photrail_Pho_Cone04NeutralHadronIso_dR7_dEta0_pt0/F");
  OutputTree->Branch("photrail_Pho_Cone04NeutralHadronIso_dR7_dEta0_pt5",&photrail_Pho_Cone04NeutralHadronIso_dR7_dEta0_pt5,"photrail_Pho_Cone04NeutralHadronIso_dR7_dEta0_pt5/F");
  OutputTree->Branch("photrail_Pho_Cone01NeutralHadronIso_dR0_dEta0_pt0_mvVtx",&photrail_Pho_Cone01NeutralHadronIso_dR0_dEta0_pt0_mvVtx,"photrail_Pho_Cone01NeutralHadronIso_dR0_dEta0_pt0_mvVtx/F");
  OutputTree->Branch("photrail_Pho_Cone02NeutralHadronIso_dR0_dEta0_pt0_mvVtx",&photrail_Pho_Cone02NeutralHadronIso_dR0_dEta0_pt0_mvVtx,"photrail_Pho_Cone02NeutralHadronIso_dR0_dEta0_pt0_mvVtx/F");
  OutputTree->Branch("photrail_Pho_Cone03NeutralHadronIso_dR0_dEta0_pt0_mvVtx",&photrail_Pho_Cone03NeutralHadronIso_dR0_dEta0_pt0_mvVtx,"photrail_Pho_Cone03NeutralHadronIso_dR0_dEta0_pt0_mvVtx/F");
  OutputTree->Branch("photrail_Pho_Cone04NeutralHadronIso_dR0_dEta0_pt0_mvVtx",&photrail_Pho_Cone04NeutralHadronIso_dR0_dEta0_pt0_mvVtx,"photrail_Pho_Cone04NeutralHadronIso_dR0_dEta0_pt0_mvVtx/F");
  OutputTree->Branch("photrail_Pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_dz0_old",&photrail_Pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_dz0_old,"photrail_Pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_dz0_old/F");
  OutputTree->Branch("photrail_Pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_PFnoPU_old",&photrail_Pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_PFnoPU_old,"photrail_Pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_PFnoPU_old/F");
  OutputTree->Branch("photrail_Pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_dz0_old",&photrail_Pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_dz0_old,"photrail_Pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_dz0_old/F");
  OutputTree->Branch("photrail_Pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_PFnoPU_old",&photrail_Pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_PFnoPU_old,"photrail_Pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_PFnoPU_old/F");
  OutputTree->Branch("photrail_Pho_Cone01ChargedHadronIso_dR0_dEta0_pt0_dz0",&photrail_Pho_Cone01ChargedHadronIso_dR0_dEta0_pt0_dz0,"photrail_Pho_Cone01ChargedHadronIso_dR0_dEta0_pt0_dz0/F");
  OutputTree->Branch("photrail_Pho_Cone01ChargedHadronIso_dR0_dEta0_pt0_dz1_dxy01",&photrail_Pho_Cone01ChargedHadronIso_dR0_dEta0_pt0_dz1_dxy01,"photrail_Pho_Cone01ChargedHadronIso_dR0_dEta0_pt0_dz1_dxy01/F");
  OutputTree->Branch("photrail_Pho_Cone01ChargedHadronIso_dR0_dEta0_pt0_PFnoPU",&photrail_Pho_Cone01ChargedHadronIso_dR0_dEta0_pt0_PFnoPU,"photrail_Pho_Cone01ChargedHadronIso_dR0_dEta0_pt0_PFnoPU/F");
  OutputTree->Branch("photrail_Pho_Cone01ChargedHadronIso_dR015_dEta0_pt0_dz0",&photrail_Pho_Cone01ChargedHadronIso_dR015_dEta0_pt0_dz0,"photrail_Pho_Cone01ChargedHadronIso_dR015_dEta0_pt0_dz0/F");
  OutputTree->Branch("photrail_Pho_Cone01ChargedHadronIso_dR015_dEta0_pt0_dz1_dxy01",&photrail_Pho_Cone01ChargedHadronIso_dR015_dEta0_pt0_dz1_dxy01,"photrail_Pho_Cone01ChargedHadronIso_dR015_dEta0_pt0_dz1_dxy01/F");
  OutputTree->Branch("photrail_Pho_Cone01ChargedHadronIso_dR015_dEta0_pt0_PFnoPU",&photrail_Pho_Cone01ChargedHadronIso_dR015_dEta0_pt0_PFnoPU,"photrail_Pho_Cone01ChargedHadronIso_dR015_dEta0_pt0_PFnoPU/F");
  OutputTree->Branch("photrail_Pho_Cone02ChargedHadronIso_dR0_dEta0_pt0_dz0",&photrail_Pho_Cone02ChargedHadronIso_dR0_dEta0_pt0_dz0,"photrail_Pho_Cone02ChargedHadronIso_dR0_dEta0_pt0_dz0/F");
  OutputTree->Branch("photrail_Pho_Cone02ChargedHadronIso_dR0_dEta0_pt0_dz1_dxy01",&photrail_Pho_Cone02ChargedHadronIso_dR0_dEta0_pt0_dz1_dxy01,"photrail_Pho_Cone02ChargedHadronIso_dR0_dEta0_pt0_dz1_dxy01/F");
  OutputTree->Branch("photrail_Pho_Cone02ChargedHadronIso_dR0_dEta0_pt0_PFnoPU",&photrail_Pho_Cone02ChargedHadronIso_dR0_dEta0_pt0_PFnoPU,"photrail_Pho_Cone02ChargedHadronIso_dR0_dEta0_pt0_PFnoPU/F");
  OutputTree->Branch("photrail_Pho_Cone02ChargedHadronIso_dR015_dEta0_pt0_dz0",&photrail_Pho_Cone02ChargedHadronIso_dR015_dEta0_pt0_dz0,"photrail_Pho_Cone02ChargedHadronIso_dR015_dEta0_pt0_dz0/F");
  OutputTree->Branch("photrail_Pho_Cone02ChargedHadronIso_dR015_dEta0_pt0_dz1_dxy01",&photrail_Pho_Cone02ChargedHadronIso_dR015_dEta0_pt0_dz1_dxy01,"photrail_Pho_Cone02ChargedHadronIso_dR015_dEta0_pt0_dz1_dxy01/F");
  OutputTree->Branch("photrail_Pho_Cone02ChargedHadronIso_dR015_dEta0_pt0_PFnoPU",&photrail_Pho_Cone02ChargedHadronIso_dR015_dEta0_pt0_PFnoPU,"photrail_Pho_Cone02ChargedHadronIso_dR015_dEta0_pt0_PFnoPU/F");
  OutputTree->Branch("photrail_Pho_Cone03ChargedHadronIso_dR0_dEta0_pt0_dz0",&photrail_Pho_Cone03ChargedHadronIso_dR0_dEta0_pt0_dz0,"photrail_Pho_Cone03ChargedHadronIso_dR0_dEta0_pt0_dz0/F");
  OutputTree->Branch("photrail_Pho_Cone03ChargedHadronIso_dR0_dEta0_pt0_dz1_dxy01",&photrail_Pho_Cone03ChargedHadronIso_dR0_dEta0_pt0_dz1_dxy01,"photrail_Pho_Cone03ChargedHadronIso_dR0_dEta0_pt0_dz1_dxy01/F");
  OutputTree->Branch("photrail_Pho_Cone03ChargedHadronIso_dR0_dEta0_pt0_PFnoPU",&photrail_Pho_Cone03ChargedHadronIso_dR0_dEta0_pt0_PFnoPU,"photrail_Pho_Cone03ChargedHadronIso_dR0_dEta0_pt0_PFnoPU/F");
  OutputTree->Branch("photrail_Pho_Cone03ChargedHadronIso_dR015_dEta0_pt0_dz0",&photrail_Pho_Cone03ChargedHadronIso_dR015_dEta0_pt0_dz0,"photrail_Pho_Cone03ChargedHadronIso_dR015_dEta0_pt0_dz0/F");
  OutputTree->Branch("photrail_Pho_Cone03ChargedHadronIso_dR015_dEta0_pt0_dz1_dxy01",&photrail_Pho_Cone03ChargedHadronIso_dR015_dEta0_pt0_dz1_dxy01,"photrail_Pho_Cone03ChargedHadronIso_dR015_dEta0_pt0_dz1_dxy01/F");
  OutputTree->Branch("photrail_Pho_Cone03ChargedHadronIso_dR015_dEta0_pt0_PFnoPU",&photrail_Pho_Cone03ChargedHadronIso_dR015_dEta0_pt0_PFnoPU,"photrail_Pho_Cone03ChargedHadronIso_dR015_dEta0_pt0_PFnoPU/F");
  OutputTree->Branch("photrail_Pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_dz0",&photrail_Pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_dz0,"photrail_Pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_dz0/F");
  OutputTree->Branch("photrail_Pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_dz1_dxy01",&photrail_Pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_dz1_dxy01,"photrail_Pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_dz1_dxy01/F");
  OutputTree->Branch("photrail_Pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_PFnoPU",&photrail_Pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_PFnoPU,"photrail_Pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_PFnoPU/F");
  OutputTree->Branch("photrail_Pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_dz0",&photrail_Pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_dz0,"photrail_Pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_dz0/F");
  OutputTree->Branch("photrail_Pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_dz1_dxy01",&photrail_Pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_dz1_dxy01,"photrail_Pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_dz1_dxy01/F");
  OutputTree->Branch("photrail_Pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_PFnoPU",&photrail_Pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_PFnoPU,"photrail_Pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_PFnoPU/F");

  OutputTree->Branch("pholead_PhoIso03Ecal",&pholead_PhoIso03Ecal,"pholead_PhoIso03Ecal/F");
  OutputTree->Branch("pholead_PhoIso03Hcal",&pholead_PhoIso03Hcal,"pholead_PhoIso03Hcal/F");
  OutputTree->Branch("pholead_PhoIso03TrkSolid",&pholead_PhoIso03TrkSolid,"pholead_PhoIso03TrkSolid/F");
  OutputTree->Branch("pholead_PhoIso03TrkHollow",&pholead_PhoIso03TrkHollow,"pholead_PhoIso03TrkHollow/F");
  OutputTree->Branch("pholead_PhoIso03",&pholead_PhoIso03,"pholead_PhoIso03/F");
  OutputTree->Branch("pholead_PhoIso04Ecal",&pholead_PhoIso04Ecal,"pholead_PhoIso04Ecal/F");
  OutputTree->Branch("pholead_PhoIso04Hcal",&pholead_PhoIso04Hcal,"pholead_PhoIso04Hcal/F");
  OutputTree->Branch("pholead_PhoIso04TrkSolid",&pholead_PhoIso04TrkSolid,"pholead_PhoIso04TrkSolid/F");
  OutputTree->Branch("pholead_PhoIso04TrkHollow",&pholead_PhoIso04TrkHollow,"pholead_PhoIso04TrkHollow/F");
  OutputTree->Branch("pholead_PhoIso04",&pholead_PhoIso04,"pholead_PhoIso04/F");


  OutputTree->Branch("photrail_PhoIso03Ecal",&photrail_PhoIso03Ecal,"photrail_PhoIso03Ecal/F");
  OutputTree->Branch("photrail_PhoIso03Hcal",&photrail_PhoIso03Hcal,"photrail_PhoIso03Hcal/F");
  OutputTree->Branch("photrail_PhoIso03TrkSolid",&photrail_PhoIso03TrkSolid,"photrail_PhoIso03TrkSolid/F");
  OutputTree->Branch("photrail_PhoIso03TrkHollow",&photrail_PhoIso03TrkHollow,"photrail_PhoIso03TrkHollow/F");
  OutputTree->Branch("photrail_PhoIso03",&photrail_PhoIso03,"photrail_PhoIso03/F");
  OutputTree->Branch("photrail_PhoIso04Ecal",&photrail_PhoIso04Ecal,"photrail_PhoIso04Ecal/F");
  OutputTree->Branch("photrail_PhoIso04Hcal",&photrail_PhoIso04Hcal,"photrail_PhoIso04Hcal/F");
  OutputTree->Branch("photrail_PhoIso04TrkSolid",&photrail_PhoIso04TrkSolid,"photrail_PhoIso04TrkSolid/F");
  OutputTree->Branch("photrail_PhoIso04TrkHollow",&photrail_PhoIso04TrkHollow,"photrail_PhoIso04TrkHollow/F");
  OutputTree->Branch("photrail_PhoIso04",&photrail_PhoIso04,"photrail_PhoIso04/F");


  OutputTree->Branch("pholead_PhoS4OverS1",&pholead_PhoS4OverS1,"pholead_PhoS4OverS1/F");
  OutputTree->Branch("pholead_PhoSigmaEtaEta",&pholead_PhoSigmaEtaEta,"pholead_PhoSigmaEtaEta/F");
  OutputTree->Branch("pholead_PhoE1x5",&pholead_PhoE1x5,"pholead_PhoE1x5/F");
  OutputTree->Branch("pholead_PhoE2x5",&pholead_PhoE2x5,"pholead_PhoE2x5/F");
  OutputTree->Branch("pholead_PhoE3x3",&pholead_PhoE3x3,"pholead_PhoE3x3/F");
  OutputTree->Branch("pholead_PhoE5x5",&pholead_PhoE5x5,"pholead_PhoE5x5/F");
  OutputTree->Branch("pholead_PhomaxEnergyXtal",&pholead_PhomaxEnergyXtal,"pholead_PhomaxEnergyXtal/F");
  OutputTree->Branch("pholead_PhoIso03HcalDepth1",&pholead_PhoIso03HcalDepth1,"pholead_PhoIso03HcalDepth1/F");
  OutputTree->Branch("pholead_PhoIso03HcalDepth2",&pholead_PhoIso03HcalDepth2,"pholead_PhoIso03HcalDepth2/F");
  OutputTree->Branch("pholead_PhoIso04HcalDepth1",&pholead_PhoIso04HcalDepth1,"pholead_PhoIso04HcalDepth1/F");
  OutputTree->Branch("pholead_PhoIso04HcalDepth2",&pholead_PhoIso04HcalDepth2,"pholead_PhoIso04HcalDepth2/F");
  OutputTree->Branch("pholead_PhoIso03nTrksSolid",&pholead_PhoIso03nTrksSolid,"pholead_PhoIso03nTrksSolid/I");
  OutputTree->Branch("pholead_PhoIso03nTrksHollow",&pholead_PhoIso03nTrksHollow,"pholead_PhoIso03nTrksHollow/I");
  OutputTree->Branch("pholead_PhoIso04nTrksSolid",&pholead_PhoIso04nTrksSolid,"pholead_PhoIso04nTrksSolid/I");
  OutputTree->Branch("pholead_PhoIso04nTrksHollow",&pholead_PhoIso04nTrksHollow,"pholead_PhoIso04nTrksHollow/I");
  OutputTree->Branch("pholead_Pho_ChargedHadronIso",&pholead_Pho_ChargedHadronIso,"pholead_Pho_ChargedHadronIso/F");
  OutputTree->Branch("pholead_Pho_NeutralHadronIso",&pholead_Pho_NeutralHadronIso,"pholead_Pho_NeutralHadronIso/F");
  OutputTree->Branch("pholead_Pho_PhotonIso",&pholead_Pho_PhotonIso,"pholead_Pho_PhotonIso/F");
  OutputTree->Branch("pholead_Pho_isPFPhoton",&pholead_Pho_isPFPhoton,"pholead_Pho_isPFPhoton/I");
  OutputTree->Branch("pholead_Pho_isPFElectron",&pholead_Pho_isPFElectron,"pholead_Pho_isPFElectron/I");

  OutputTree->Branch("photrail_PhoS4OverS1",&photrail_PhoS4OverS1,"photrail_PhoS4OverS1/F");
  OutputTree->Branch("photrail_PhoSigmaEtaEta",&photrail_PhoSigmaEtaEta,"photrail_PhoSigmaEtaEta/F");
  OutputTree->Branch("photrail_PhoE1x5",&photrail_PhoE1x5,"photrail_PhoE1x5/F");
  OutputTree->Branch("photrail_PhoE2x5",&photrail_PhoE2x5,"photrail_PhoE2x5/F");
  OutputTree->Branch("photrail_PhoE3x3",&photrail_PhoE3x3,"photrail_PhoE3x3/F");
  OutputTree->Branch("photrail_PhoE5x5",&photrail_PhoE5x5,"photrail_PhoE5x5/F");
  OutputTree->Branch("photrail_PhomaxEnergyXtal",&photrail_PhomaxEnergyXtal,"photrail_PhomaxEnergyXtal/F");
  OutputTree->Branch("photrail_PhoIso03HcalDepth1",&photrail_PhoIso03HcalDepth1,"photrail_PhoIso03HcalDepth1/F");
  OutputTree->Branch("photrail_PhoIso03HcalDepth2",&photrail_PhoIso03HcalDepth2,"photrail_PhoIso03HcalDepth2/F");
  OutputTree->Branch("photrail_PhoIso04HcalDepth1",&photrail_PhoIso04HcalDepth1,"photrail_PhoIso04HcalDepth1/F");
  OutputTree->Branch("photrail_PhoIso04HcalDepth2",&photrail_PhoIso04HcalDepth2,"photrail_PhoIso04HcalDepth2/F");
  OutputTree->Branch("photrail_PhoIso03nTrksSolid",&photrail_PhoIso03nTrksSolid,"photrail_PhoIso03nTrksSolid/I");
  OutputTree->Branch("photrail_PhoIso03nTrksHollow",&photrail_PhoIso03nTrksHollow,"photrail_PhoIso03nTrksHollow/I");
  OutputTree->Branch("photrail_PhoIso04nTrksSolid",&photrail_PhoIso04nTrksSolid,"photrail_PhoIso04nTrksSolid/I");
  OutputTree->Branch("photrail_PhoIso04nTrksHollow",&photrail_PhoIso04nTrksHollow,"photrail_PhoIso04nTrksHollow/I");
  OutputTree->Branch("photrail_Pho_ChargedHadronIso",&photrail_Pho_ChargedHadronIso,"photrail_Pho_ChargedHadronIso/F");
  OutputTree->Branch("photrail_Pho_NeutralHadronIso",&photrail_Pho_NeutralHadronIso,"photrail_Pho_NeutralHadronIso/F");
  OutputTree->Branch("photrail_Pho_PhotonIso",&photrail_Pho_PhotonIso,"photrail_Pho_PhotonIso/F");
  OutputTree->Branch("photrail_Pho_isPFPhoton",&photrail_Pho_isPFPhoton,"photrail_Pho_isPFPhoton/I");
  OutputTree->Branch("photrail_Pho_isPFElectron",&photrail_Pho_isPFElectron,"photrail_Pho_isPFElectron/I");


  OutputTree->Branch("pholead_PhoMCmatchindex",&pholead_PhoMCmatchindex,"pholead_PhoMCmatchindex/I");
  OutputTree->Branch("pholead_PhoMCmatchexitcode",&pholead_PhoMCmatchexitcode,"pholead_PhoMCmatchexitcode/I");

  OutputTree->Branch("photrail_PhoMCmatchindex",&photrail_PhoMCmatchindex,"photrail_PhoMCmatchindex/I");
  OutputTree->Branch("photrail_PhoMCmatchexitcode",&photrail_PhoMCmatchexitcode,"photrail_PhoMCmatchexitcode/I");


  fHNumPU = new TH1F("NumPU","NumPU",40,0,40);
  fHNumVtx = new TH1F("NumVtx","NumVtx",40,0,40);
	

	
  cout << "Tree and histos created" << endl;

  // for (int i=0;i<10;i++)  cout << "No trigger selection!!!" << endl;	

}

void DiPhotonMiniTree::Analyze(){

  //cout << "Analyze this event" << endl;


  //cout << "A" << endl;

  float weight;
  if (!isdata) weight = GetPUWeight(fTR->PUnumInteractions);
  else weight=1;
  
  event_luminormfactor=AddWeight;

  if (!isdata) fHNumPU->Fill(fTR->PUnumInteractions,weight);
  fHNumVtx->Fill(fTR->NVrtx,weight);

    if (!TriggerSelection()) return;


  std::vector<int> passing = PhotonSelection(fTR);

  if (passing.size()<2) return;

  float invmass0 = (CorrPhoton(fTR,passing.at(0),0)+CorrPhoton(fTR,passing.at(1),0)).M();
  float invmass5 = (CorrPhoton(fTR,passing.at(0),5)+CorrPhoton(fTR,passing.at(1),5)).M();
  float invmass6 = (CorrPhoton(fTR,passing.at(0),6)+CorrPhoton(fTR,passing.at(1),6)).M();

  event_weight = weight;
  event_rho = fTR->Rho;
  if (!isdata) event_nPU = fTR->PUnumInteractions;
  event_nRecVtx = fTR->NVrtx;
  

  dipho_mgg_photon = invmass0;
  dipho_mgg_newCorr = invmass5;
  dipho_mgg_newCorrLocal = invmass6;

  pholead_eta = fTR->PhoEta[passing.at(0)];
  photrail_eta = fTR->PhoEta[passing.at(1)];
  
  pholead_px = fTR->PhoPx[passing.at(0)];
  photrail_px = fTR->PhoPx[passing.at(1)];
  pholead_py = fTR->PhoPy[passing.at(0)];
  photrail_py = fTR->PhoPy[passing.at(1)];
  pholead_pt = fTR->PhoPt[passing.at(0)];
  photrail_pt = fTR->PhoPt[passing.at(1)];
  pholead_pz = fTR->PhoPz[passing.at(0)];
  photrail_pz = fTR->PhoPz[passing.at(1)];
  pholead_energy = fTR->PhoEnergy[passing.at(0)];
  photrail_energy = fTR->PhoEnergy[passing.at(1)];

  pholead_SCeta = fTR->SCEta[fTR->PhotSCindex[passing.at(0)]];
  photrail_SCeta = fTR->SCEta[fTR->PhotSCindex[passing.at(1)]];
 
  pholead_energySCdefault = CorrPhoton(fTR,passing.at(0),0).E();
  photrail_energySCdefault = CorrPhoton(fTR,passing.at(1),0).E();
  pholead_energyNewCorr = CorrPhoton(fTR,passing.at(0),5).E();
  photrail_energyNewCorr = CorrPhoton(fTR,passing.at(1),5).E();
  pholead_energyNewCorrLocal = CorrPhoton(fTR,passing.at(0),6).E();
  photrail_energyNewCorrLocal = CorrPhoton(fTR,passing.at(1),6).E();

  pholead_r9 = fTR->SCR9[fTR->PhotSCindex[passing.at(0)]];
  photrail_r9 = fTR->SCR9[fTR->PhotSCindex[passing.at(1)]];
  pholead_sieie = fTR->PhoSigmaIetaIeta[passing.at(0)];
  photrail_sieie = fTR->PhoSigmaIetaIeta[passing.at(1)];
  pholead_hoe = fTR->PhoHoverE[passing.at(0)];
  photrail_hoe = fTR->PhoHoverE[passing.at(1)];
  pholead_brem =  fTR->SCBrem[fTR->PhotSCindex[passing.at(0)]];
  photrail_brem = fTR->SCBrem[fTR->PhotSCindex[passing.at(1)]];
  pholead_sigmaPhi = fTR->SCPhiWidth[fTR->PhotSCindex[passing.at(0)]];
  photrail_sigmaPhi = fTR->SCPhiWidth[fTR->PhotSCindex[passing.at(1)]];
  pholead_sigmaEta = fTR->SCEtaWidth[fTR->PhotSCindex[passing.at(0)]];
  photrail_sigmaEta = fTR->SCEtaWidth[fTR->PhotSCindex[passing.at(1)]];

  pholead_Pho_Cone04PhotonIso_dR0_dEta0_pt0=fTR->Pho_Cone04PhotonIso_dR0_dEta0_pt0[passing.at(0)];
  pholead_Pho_Cone04PhotonIso_dR0_dEta0_pt5=fTR->Pho_Cone04PhotonIso_dR0_dEta0_pt5[passing.at(0)];
  pholead_Pho_Cone04PhotonIso_dR8_dEta0_pt0=fTR->Pho_Cone04PhotonIso_dR8_dEta0_pt0[passing.at(0)];
  pholead_Pho_Cone04PhotonIso_dR8_dEta0_pt5=fTR->Pho_Cone04PhotonIso_dR8_dEta0_pt5[passing.at(0)];
  pholead_Pho_Cone01PhotonIso_dR045EB070EE_dEta015_pt08EB1EE_mvVtx=fTR->Pho_Cone01PhotonIso_dR045EB070EE_dEta015_pt08EB1EE_mvVtx[passing.at(0)];
  pholead_Pho_Cone02PhotonIso_dR045EB070EE_dEta015_pt08EB1EE_mvVtx=fTR->Pho_Cone02PhotonIso_dR045EB070EE_dEta015_pt08EB1EE_mvVtx[passing.at(0)];
  pholead_Pho_Cone03PhotonIso_dR045EB070EE_dEta015_pt08EB1EE_mvVtx=fTR->Pho_Cone03PhotonIso_dR045EB070EE_dEta015_pt08EB1EE_mvVtx[passing.at(0)];
  pholead_Pho_Cone04PhotonIso_dR045EB070EE_dEta015_pt08EB1EE_mvVtx=fTR->Pho_Cone04PhotonIso_dR045EB070EE_dEta015_pt08EB1EE_mvVtx[passing.at(0)];
  pholead_Pho_Cone04NeutralHadronIso_dR0_dEta0_pt0=fTR->Pho_Cone04NeutralHadronIso_dR0_dEta0_pt0[passing.at(0)];
  pholead_Pho_Cone04NeutralHadronIso_dR0_dEta0_pt5=fTR->Pho_Cone04NeutralHadronIso_dR0_dEta0_pt5[passing.at(0)];
  pholead_Pho_Cone04NeutralHadronIso_dR0_dEta0_pt0_nocracks=fTR->Pho_Cone04NeutralHadronIso_dR0_dEta0_pt0_nocracks[passing.at(0)];
  pholead_Pho_Cone04NeutralHadronIso_dR0_dEta0_pt5_nocracks=fTR->Pho_Cone04NeutralHadronIso_dR0_dEta0_pt5_nocracks[passing.at(0)];
  pholead_Pho_Cone04NeutralHadronIso_dR7_dEta0_pt0=fTR->Pho_Cone04NeutralHadronIso_dR7_dEta0_pt0[passing.at(0)];
  pholead_Pho_Cone04NeutralHadronIso_dR7_dEta0_pt5=fTR->Pho_Cone04NeutralHadronIso_dR7_dEta0_pt5[passing.at(0)];
  pholead_Pho_Cone01NeutralHadronIso_dR0_dEta0_pt0_mvVtx=fTR->Pho_Cone01NeutralHadronIso_dR0_dEta0_pt0_mvVtx[passing.at(0)];
  pholead_Pho_Cone02NeutralHadronIso_dR0_dEta0_pt0_mvVtx=fTR->Pho_Cone02NeutralHadronIso_dR0_dEta0_pt0_mvVtx[passing.at(0)];
  pholead_Pho_Cone03NeutralHadronIso_dR0_dEta0_pt0_mvVtx=fTR->Pho_Cone03NeutralHadronIso_dR0_dEta0_pt0_mvVtx[passing.at(0)];
  pholead_Pho_Cone04NeutralHadronIso_dR0_dEta0_pt0_mvVtx=fTR->Pho_Cone04NeutralHadronIso_dR0_dEta0_pt0_mvVtx[passing.at(0)];
  pholead_Pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_dz0_old=fTR->Pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_dz0_old[passing.at(0)];
  pholead_Pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_PFnoPU_old=fTR->Pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_PFnoPU_old[passing.at(0)];
  pholead_Pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_dz0_old=fTR->Pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_dz0_old[passing.at(0)];
  pholead_Pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_PFnoPU_old=fTR->Pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_PFnoPU_old[passing.at(0)];
  pholead_Pho_Cone01ChargedHadronIso_dR0_dEta0_pt0_dz0=fTR->Pho_Cone01ChargedHadronIso_dR0_dEta0_pt0_dz0[passing.at(0)];
  pholead_Pho_Cone01ChargedHadronIso_dR0_dEta0_pt0_dz1_dxy01=fTR->Pho_Cone01ChargedHadronIso_dR0_dEta0_pt0_dz1_dxy01[passing.at(0)];
  pholead_Pho_Cone01ChargedHadronIso_dR0_dEta0_pt0_PFnoPU=fTR->Pho_Cone01ChargedHadronIso_dR0_dEta0_pt0_PFnoPU[passing.at(0)];
  pholead_Pho_Cone01ChargedHadronIso_dR015_dEta0_pt0_dz0=fTR->Pho_Cone01ChargedHadronIso_dR015_dEta0_pt0_dz0[passing.at(0)];
  pholead_Pho_Cone01ChargedHadronIso_dR015_dEta0_pt0_dz1_dxy01=fTR->Pho_Cone01ChargedHadronIso_dR015_dEta0_pt0_dz1_dxy01[passing.at(0)];
  pholead_Pho_Cone01ChargedHadronIso_dR015_dEta0_pt0_PFnoPU=fTR->Pho_Cone01ChargedHadronIso_dR015_dEta0_pt0_PFnoPU[passing.at(0)];
  pholead_Pho_Cone02ChargedHadronIso_dR0_dEta0_pt0_dz0=fTR->Pho_Cone02ChargedHadronIso_dR0_dEta0_pt0_dz0[passing.at(0)];
  pholead_Pho_Cone02ChargedHadronIso_dR0_dEta0_pt0_dz1_dxy01=fTR->Pho_Cone02ChargedHadronIso_dR0_dEta0_pt0_dz1_dxy01[passing.at(0)];
  pholead_Pho_Cone02ChargedHadronIso_dR0_dEta0_pt0_PFnoPU=fTR->Pho_Cone02ChargedHadronIso_dR0_dEta0_pt0_PFnoPU[passing.at(0)];
  pholead_Pho_Cone02ChargedHadronIso_dR015_dEta0_pt0_dz0=fTR->Pho_Cone02ChargedHadronIso_dR015_dEta0_pt0_dz0[passing.at(0)];
  pholead_Pho_Cone02ChargedHadronIso_dR015_dEta0_pt0_dz1_dxy01=fTR->Pho_Cone02ChargedHadronIso_dR015_dEta0_pt0_dz1_dxy01[passing.at(0)];
  pholead_Pho_Cone02ChargedHadronIso_dR015_dEta0_pt0_PFnoPU=fTR->Pho_Cone02ChargedHadronIso_dR015_dEta0_pt0_PFnoPU[passing.at(0)];
  pholead_Pho_Cone03ChargedHadronIso_dR0_dEta0_pt0_dz0=fTR->Pho_Cone03ChargedHadronIso_dR0_dEta0_pt0_dz0[passing.at(0)];
  pholead_Pho_Cone03ChargedHadronIso_dR0_dEta0_pt0_dz1_dxy01=fTR->Pho_Cone03ChargedHadronIso_dR0_dEta0_pt0_dz1_dxy01[passing.at(0)];
  pholead_Pho_Cone03ChargedHadronIso_dR0_dEta0_pt0_PFnoPU=fTR->Pho_Cone03ChargedHadronIso_dR0_dEta0_pt0_PFnoPU[passing.at(0)];
  pholead_Pho_Cone03ChargedHadronIso_dR015_dEta0_pt0_dz0=fTR->Pho_Cone03ChargedHadronIso_dR015_dEta0_pt0_dz0[passing.at(0)];
  pholead_Pho_Cone03ChargedHadronIso_dR015_dEta0_pt0_dz1_dxy01=fTR->Pho_Cone03ChargedHadronIso_dR015_dEta0_pt0_dz1_dxy01[passing.at(0)];
  pholead_Pho_Cone03ChargedHadronIso_dR015_dEta0_pt0_PFnoPU=fTR->Pho_Cone03ChargedHadronIso_dR015_dEta0_pt0_PFnoPU[passing.at(0)];
  pholead_Pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_dz0=fTR->Pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_dz0[passing.at(0)];
  pholead_Pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_dz1_dxy01=fTR->Pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_dz1_dxy01[passing.at(0)];
  pholead_Pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_PFnoPU=fTR->Pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_PFnoPU[passing.at(0)];
  pholead_Pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_dz0=fTR->Pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_dz0[passing.at(0)];
  pholead_Pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_dz1_dxy01=fTR->Pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_dz1_dxy01[passing.at(0)];
  pholead_Pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_PFnoPU=fTR->Pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_PFnoPU[passing.at(0)];

  photrail_Pho_Cone04PhotonIso_dR0_dEta0_pt0=fTR->Pho_Cone04PhotonIso_dR0_dEta0_pt0[passing.at(1)];
  photrail_Pho_Cone04PhotonIso_dR0_dEta0_pt5=fTR->Pho_Cone04PhotonIso_dR0_dEta0_pt5[passing.at(1)];
  photrail_Pho_Cone04PhotonIso_dR8_dEta0_pt0=fTR->Pho_Cone04PhotonIso_dR8_dEta0_pt0[passing.at(1)];
  photrail_Pho_Cone04PhotonIso_dR8_dEta0_pt5=fTR->Pho_Cone04PhotonIso_dR8_dEta0_pt5[passing.at(1)];
  photrail_Pho_Cone01PhotonIso_dR045EB070EE_dEta015_pt08EB1EE_mvVtx=fTR->Pho_Cone01PhotonIso_dR045EB070EE_dEta015_pt08EB1EE_mvVtx[passing.at(1)];
  photrail_Pho_Cone02PhotonIso_dR045EB070EE_dEta015_pt08EB1EE_mvVtx=fTR->Pho_Cone02PhotonIso_dR045EB070EE_dEta015_pt08EB1EE_mvVtx[passing.at(1)];
  photrail_Pho_Cone03PhotonIso_dR045EB070EE_dEta015_pt08EB1EE_mvVtx=fTR->Pho_Cone03PhotonIso_dR045EB070EE_dEta015_pt08EB1EE_mvVtx[passing.at(1)];
  photrail_Pho_Cone04PhotonIso_dR045EB070EE_dEta015_pt08EB1EE_mvVtx=fTR->Pho_Cone04PhotonIso_dR045EB070EE_dEta015_pt08EB1EE_mvVtx[passing.at(1)];
  photrail_Pho_Cone04NeutralHadronIso_dR0_dEta0_pt0=fTR->Pho_Cone04NeutralHadronIso_dR0_dEta0_pt0[passing.at(1)];
  photrail_Pho_Cone04NeutralHadronIso_dR0_dEta0_pt5=fTR->Pho_Cone04NeutralHadronIso_dR0_dEta0_pt5[passing.at(1)];
  photrail_Pho_Cone04NeutralHadronIso_dR0_dEta0_pt0_nocracks=fTR->Pho_Cone04NeutralHadronIso_dR0_dEta0_pt0_nocracks[passing.at(1)];
  photrail_Pho_Cone04NeutralHadronIso_dR0_dEta0_pt5_nocracks=fTR->Pho_Cone04NeutralHadronIso_dR0_dEta0_pt5_nocracks[passing.at(1)];
  photrail_Pho_Cone04NeutralHadronIso_dR7_dEta0_pt0=fTR->Pho_Cone04NeutralHadronIso_dR7_dEta0_pt0[passing.at(1)];
  photrail_Pho_Cone04NeutralHadronIso_dR7_dEta0_pt5=fTR->Pho_Cone04NeutralHadronIso_dR7_dEta0_pt5[passing.at(1)];
  photrail_Pho_Cone01NeutralHadronIso_dR0_dEta0_pt0_mvVtx=fTR->Pho_Cone01NeutralHadronIso_dR0_dEta0_pt0_mvVtx[passing.at(1)];
  photrail_Pho_Cone02NeutralHadronIso_dR0_dEta0_pt0_mvVtx=fTR->Pho_Cone02NeutralHadronIso_dR0_dEta0_pt0_mvVtx[passing.at(1)];
  photrail_Pho_Cone03NeutralHadronIso_dR0_dEta0_pt0_mvVtx=fTR->Pho_Cone03NeutralHadronIso_dR0_dEta0_pt0_mvVtx[passing.at(1)];
  photrail_Pho_Cone04NeutralHadronIso_dR0_dEta0_pt0_mvVtx=fTR->Pho_Cone04NeutralHadronIso_dR0_dEta0_pt0_mvVtx[passing.at(1)];
  photrail_Pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_dz0_old=fTR->Pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_dz0_old[passing.at(1)];
  photrail_Pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_PFnoPU_old=fTR->Pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_PFnoPU_old[passing.at(1)];
  photrail_Pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_dz0_old=fTR->Pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_dz0_old[passing.at(1)];
  photrail_Pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_PFnoPU_old=fTR->Pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_PFnoPU_old[passing.at(1)];
  photrail_Pho_Cone01ChargedHadronIso_dR0_dEta0_pt0_dz0=fTR->Pho_Cone01ChargedHadronIso_dR0_dEta0_pt0_dz0[passing.at(1)];
  photrail_Pho_Cone01ChargedHadronIso_dR0_dEta0_pt0_dz1_dxy01=fTR->Pho_Cone01ChargedHadronIso_dR0_dEta0_pt0_dz1_dxy01[passing.at(1)];
  photrail_Pho_Cone01ChargedHadronIso_dR0_dEta0_pt0_PFnoPU=fTR->Pho_Cone01ChargedHadronIso_dR0_dEta0_pt0_PFnoPU[passing.at(1)];
  photrail_Pho_Cone01ChargedHadronIso_dR015_dEta0_pt0_dz0=fTR->Pho_Cone01ChargedHadronIso_dR015_dEta0_pt0_dz0[passing.at(1)];
  photrail_Pho_Cone01ChargedHadronIso_dR015_dEta0_pt0_dz1_dxy01=fTR->Pho_Cone01ChargedHadronIso_dR015_dEta0_pt0_dz1_dxy01[passing.at(1)];
  photrail_Pho_Cone01ChargedHadronIso_dR015_dEta0_pt0_PFnoPU=fTR->Pho_Cone01ChargedHadronIso_dR015_dEta0_pt0_PFnoPU[passing.at(1)];
  photrail_Pho_Cone02ChargedHadronIso_dR0_dEta0_pt0_dz0=fTR->Pho_Cone02ChargedHadronIso_dR0_dEta0_pt0_dz0[passing.at(1)];
  photrail_Pho_Cone02ChargedHadronIso_dR0_dEta0_pt0_dz1_dxy01=fTR->Pho_Cone02ChargedHadronIso_dR0_dEta0_pt0_dz1_dxy01[passing.at(1)];
  photrail_Pho_Cone02ChargedHadronIso_dR0_dEta0_pt0_PFnoPU=fTR->Pho_Cone02ChargedHadronIso_dR0_dEta0_pt0_PFnoPU[passing.at(1)];
  photrail_Pho_Cone02ChargedHadronIso_dR015_dEta0_pt0_dz0=fTR->Pho_Cone02ChargedHadronIso_dR015_dEta0_pt0_dz0[passing.at(1)];
  photrail_Pho_Cone02ChargedHadronIso_dR015_dEta0_pt0_dz1_dxy01=fTR->Pho_Cone02ChargedHadronIso_dR015_dEta0_pt0_dz1_dxy01[passing.at(1)];
  photrail_Pho_Cone02ChargedHadronIso_dR015_dEta0_pt0_PFnoPU=fTR->Pho_Cone02ChargedHadronIso_dR015_dEta0_pt0_PFnoPU[passing.at(1)];
  photrail_Pho_Cone03ChargedHadronIso_dR0_dEta0_pt0_dz0=fTR->Pho_Cone03ChargedHadronIso_dR0_dEta0_pt0_dz0[passing.at(1)];
  photrail_Pho_Cone03ChargedHadronIso_dR0_dEta0_pt0_dz1_dxy01=fTR->Pho_Cone03ChargedHadronIso_dR0_dEta0_pt0_dz1_dxy01[passing.at(1)];
  photrail_Pho_Cone03ChargedHadronIso_dR0_dEta0_pt0_PFnoPU=fTR->Pho_Cone03ChargedHadronIso_dR0_dEta0_pt0_PFnoPU[passing.at(1)];
  photrail_Pho_Cone03ChargedHadronIso_dR015_dEta0_pt0_dz0=fTR->Pho_Cone03ChargedHadronIso_dR015_dEta0_pt0_dz0[passing.at(1)];
  photrail_Pho_Cone03ChargedHadronIso_dR015_dEta0_pt0_dz1_dxy01=fTR->Pho_Cone03ChargedHadronIso_dR015_dEta0_pt0_dz1_dxy01[passing.at(1)];
  photrail_Pho_Cone03ChargedHadronIso_dR015_dEta0_pt0_PFnoPU=fTR->Pho_Cone03ChargedHadronIso_dR015_dEta0_pt0_PFnoPU[passing.at(1)];
  photrail_Pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_dz0=fTR->Pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_dz0[passing.at(1)];
  photrail_Pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_dz1_dxy01=fTR->Pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_dz1_dxy01[passing.at(1)];
  photrail_Pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_PFnoPU=fTR->Pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_PFnoPU[passing.at(1)];
  photrail_Pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_dz0=fTR->Pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_dz0[passing.at(1)];
  photrail_Pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_dz1_dxy01=fTR->Pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_dz1_dxy01[passing.at(1)];
  photrail_Pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_PFnoPU=fTR->Pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_PFnoPU[passing.at(1)];

  pholead_PhoIso03Ecal=fTR->PhoIso03Ecal[passing.at(0)];
  pholead_PhoIso03Hcal=fTR->PhoIso03Hcal[passing.at(0)];
  pholead_PhoIso03TrkSolid=fTR->PhoIso03TrkSolid[passing.at(0)];
  pholead_PhoIso03TrkHollow=fTR->PhoIso03TrkHollow[passing.at(0)];
  pholead_PhoIso03=fTR->PhoIso03[passing.at(0)];
  pholead_PhoIso04Ecal=fTR->PhoIso04Ecal[passing.at(0)];
  pholead_PhoIso04Hcal=fTR->PhoIso04Hcal[passing.at(0)];
  pholead_PhoIso04TrkSolid=fTR->PhoIso04TrkSolid[passing.at(0)];
  pholead_PhoIso04TrkHollow=fTR->PhoIso04TrkHollow[passing.at(0)];
  pholead_PhoIso04=fTR->PhoIso04[passing.at(0)];

  photrail_PhoIso03Ecal=fTR->PhoIso03Ecal[passing.at(1)];
  photrail_PhoIso03Hcal=fTR->PhoIso03Hcal[passing.at(1)];
  photrail_PhoIso03TrkSolid=fTR->PhoIso03TrkSolid[passing.at(1)];
  photrail_PhoIso03TrkHollow=fTR->PhoIso03TrkHollow[passing.at(1)];
  photrail_PhoIso03=fTR->PhoIso03[passing.at(1)];
  photrail_PhoIso04Ecal=fTR->PhoIso04Ecal[passing.at(1)];
  photrail_PhoIso04Hcal=fTR->PhoIso04Hcal[passing.at(1)];
  photrail_PhoIso04TrkSolid=fTR->PhoIso04TrkSolid[passing.at(1)];
  photrail_PhoIso04TrkHollow=fTR->PhoIso04TrkHollow[passing.at(1)];
  photrail_PhoIso04=fTR->PhoIso04[passing.at(1)];



  pholead_PhoE1OverE9=fTR->PhoE1OverE9[passing.at(0)];
  pholead_PhoS4OverS1=fTR->PhoS4OverS1[passing.at(0)];
  pholead_PhoSigmaEtaEta=fTR->PhoSigmaEtaEta[passing.at(0)];
  pholead_PhoE1x5=fTR->PhoE1x5[passing.at(0)];
  pholead_PhoE2x5=fTR->PhoE2x5[passing.at(0)];
  pholead_PhoE3x3=fTR->PhoE3x3[passing.at(0)];
  pholead_PhoE5x5=fTR->PhoE5x5[passing.at(0)];
  pholead_PhomaxEnergyXtal=fTR->PhomaxEnergyXtal[passing.at(0)];
  pholead_PhoIso03HcalDepth1=fTR->PhoIso03HcalDepth1[passing.at(0)];
  pholead_PhoIso03HcalDepth2=fTR->PhoIso03HcalDepth2[passing.at(0)];
  pholead_PhoIso04HcalDepth1=fTR->PhoIso04HcalDepth1[passing.at(0)];
  pholead_PhoIso04HcalDepth2=fTR->PhoIso04HcalDepth2[passing.at(0)];
  pholead_PhoIso03nTrksSolid=fTR->PhoIso03nTrksSolid[passing.at(0)];
  pholead_PhoIso03nTrksHollow=fTR->PhoIso03nTrksHollow[passing.at(0)];
  pholead_PhoIso04nTrksSolid=fTR->PhoIso04nTrksSolid[passing.at(0)];
  pholead_PhoIso04nTrksHollow=fTR->PhoIso04nTrksHollow[passing.at(0)];
  pholead_Pho_ChargedHadronIso=fTR->Pho_ChargedHadronIso[passing.at(0)];
  pholead_Pho_NeutralHadronIso=fTR->Pho_NeutralHadronIso[passing.at(0)];
  pholead_Pho_PhotonIso=fTR->Pho_PhotonIso[passing.at(0)];
  pholead_Pho_isPFPhoton=fTR->Pho_isPFPhoton[passing.at(0)];
  pholead_Pho_isPFElectron=fTR->Pho_isPFElectron[passing.at(0)];

  photrail_PhoE1OverE9=fTR->PhoE1OverE9[passing.at(1)];
  photrail_PhoS4OverS1=fTR->PhoS4OverS1[passing.at(1)];
  photrail_PhoSigmaEtaEta=fTR->PhoSigmaEtaEta[passing.at(1)];
  photrail_PhoE1x5=fTR->PhoE1x5[passing.at(1)];
  photrail_PhoE2x5=fTR->PhoE2x5[passing.at(1)];
  photrail_PhoE3x3=fTR->PhoE3x3[passing.at(1)];
  photrail_PhoE5x5=fTR->PhoE5x5[passing.at(1)];
  photrail_PhomaxEnergyXtal=fTR->PhomaxEnergyXtal[passing.at(1)];
  photrail_PhoIso03HcalDepth1=fTR->PhoIso03HcalDepth1[passing.at(1)];
  photrail_PhoIso03HcalDepth2=fTR->PhoIso03HcalDepth2[passing.at(1)];
  photrail_PhoIso04HcalDepth1=fTR->PhoIso04HcalDepth1[passing.at(1)];
  photrail_PhoIso04HcalDepth2=fTR->PhoIso04HcalDepth2[passing.at(1)];
  photrail_PhoIso03nTrksSolid=fTR->PhoIso03nTrksSolid[passing.at(1)];
  photrail_PhoIso03nTrksHollow=fTR->PhoIso03nTrksHollow[passing.at(1)];
  photrail_PhoIso04nTrksSolid=fTR->PhoIso04nTrksSolid[passing.at(1)];
  photrail_PhoIso04nTrksHollow=fTR->PhoIso04nTrksHollow[passing.at(1)];
  photrail_Pho_ChargedHadronIso=fTR->Pho_ChargedHadronIso[passing.at(1)];
  photrail_Pho_NeutralHadronIso=fTR->Pho_NeutralHadronIso[passing.at(1)];
  photrail_Pho_PhotonIso=fTR->Pho_PhotonIso[passing.at(1)];
  photrail_Pho_isPFPhoton=fTR->Pho_isPFPhoton[passing.at(1)];
  photrail_Pho_isPFElectron=fTR->Pho_isPFElectron[passing.at(1)];



  pholead_PhoMCmatchindex=fTR->PhoMCmatchindex[passing.at(0)];
  pholead_PhoMCmatchexitcode=fTR->PhoMCmatchexitcode[passing.at(0)];

  photrail_PhoMCmatchindex=fTR->PhoMCmatchindex[passing.at(1)];
  photrail_PhoMCmatchexitcode=fTR->PhoMCmatchexitcode[passing.at(1)];







  OutputTree->Fill();

 
 
 

}

void DiPhotonMiniTree::End(){
  fMiniTree->cd();
  OutputTree->Write();	
  fHNumPU->Write();
  fHNumVtx->Write();
	
  fMiniTree->Close();

}

TLorentzVector DiPhotonMiniTree::CorrPhoton(TreeReader *fTR, int i, int mode){
  TLorentzVector corr;
  corr.SetPtEtaPhiE(fTR->PhoPt[i],fTR->PhoEta[i],fTR->PhoPhi[i],fTR->PhoEnergy[i]);
  float corrE=phocorr->get_correctedenergy(fTR,i,mode);
  corr.SetE(corrE);
  corr.SetRho(corrE);
  return corr;
};

std::vector<int> DiPhotonMiniTree::PhotonSelection(TreeReader *fTR){

  vector<int> passing;
  for (int i=0; i<fTR->NPhotons; i++){
    passing.push_back(i);
  }

  for (vector<int>::iterator it = passing.begin(); it != passing.end(); ){
    if (fTR->PhotSCindex[*it]==-1) it=passing.erase(it); else it++;
  }

  for (vector<int>::iterator it = passing.begin(); it != passing.end(); ){
    float eta=fTR->SCEta[fTR->PhotSCindex[*it]];
    if (fabs(eta)>1.4442 && fabs(eta)<1.56) it=passing.erase(it); else it++;
  }

  for (vector<int>::iterator it = passing.begin(); it != passing.end(); ){
    float energy=fTR->SCRaw[fTR->PhotSCindex[*it]];
    float eta=fTR->SCEta[fTR->PhotSCindex[*it]];
    if (fabs(eta)<1.4442) energy*=phocorr->getEtaCorrectionBarrel(eta);
    if (fabs(eta)>1.56) energy+=fTR->SCPre[fTR->PhotSCindex[*it]];
    if (energy/cosh(eta)<20 || fabs(eta)>2.5) it=passing.erase(it); else it++;
  }

  return passing;

};

bool DiPhotonMiniTree::TriggerSelection(){
#include "DiPhotonTriggerSelection.cc"
};




