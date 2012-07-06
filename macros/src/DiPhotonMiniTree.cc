#include "helper/Utilities.hh"
#include "DiPhotonMiniTree.hh"

#include "DiPhotonPurity.hh"


#include <iostream>
#include <fstream>
#include <assert.h>

using namespace std;

DiPhotonMiniTree::DiPhotonMiniTree(TreeReader *tr, std::string dataType, Float_t aw, Float_t* _kfac) : UserAnalysisBase(tr), fDataType_(dataType), AddWeight(aw), kfactors(_kfac){
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

  fOutputFile->cd();

  OutputTree[0] = new TTree("Tree_standard_sel","Tree_standard_sel");
  OutputTree[1] = new TTree("Tree_signal_template","Tree_signal_template");
  OutputTree[2] = new TTree("Tree_background_template","Tree_background_template");
  OutputTree[3] = new TTree("Tree_DY_sel","Tree_DY_sel");

  for (int i=0; i<4; i++){

  OutputTree[i]->Branch("event_luminormfactor",&event_luminormfactor,"event_luminormfactor/F");
  OutputTree[i]->Branch("event_Kfactor",&event_Kfactor,"event_Kfactor/F");

  OutputTree[i]->Branch("event_weight",&event_weight,"event_weight/F");
  OutputTree[i]->Branch("event_rho",&event_rho,"event_rho/F");
  OutputTree[i]->Branch("event_nPU",&event_nPU,"event_nPU/I");
  OutputTree[i]->Branch("event_nRecVtx",&event_nRecVtx,"event_nRecVtx/I");

  OutputTree[i]->Branch("dipho_mgg_photon",&dipho_mgg_photon,"dipho_mgg_photon/F");
  OutputTree[i]->Branch("dipho_mgg_newCorr",&dipho_mgg_newCorr,"dipho_mgg_newCorr/F");
  OutputTree[i]->Branch("dipho_mgg_newCorrLocal",&dipho_mgg_newCorrLocal,"dipho_mgg_newCorrLocal/F");

  OutputTree[i]->Branch("pholead_eta",&pholead_eta,"pholead_eta/F");
  OutputTree[i]->Branch("photrail_eta",&photrail_eta,"photrail_eta/F");
  OutputTree[i]->Branch("pholead_px",&pholead_px,"pholead_px/F");
  OutputTree[i]->Branch("photrail_px",&photrail_px,"photrail_px/F");
  OutputTree[i]->Branch("pholead_py",&pholead_py,"pholead_py/F");
  OutputTree[i]->Branch("photrail_py",&photrail_py,"photrail_py/F");
  OutputTree[i]->Branch("pholead_pt",&pholead_pt,"pholead_pt/F");
  OutputTree[i]->Branch("photrail_pt",&photrail_pt,"photrail_pt/F");
  OutputTree[i]->Branch("pholead_pz",&pholead_pz,"pholead_pz/F");
  OutputTree[i]->Branch("photrail_pz",&photrail_pz,"photrail_pz/F");
  OutputTree[i]->Branch("pholead_energy",&pholead_energy,"pholead_energy/F");
  OutputTree[i]->Branch("photrail_energy",&photrail_energy,"photrail_energy/F");
  OutputTree[i]->Branch("pholead_energySCdefault",&pholead_energySCdefault,"pholead_energySCdefault/F");
  OutputTree[i]->Branch("photrail_energySCdefault",&photrail_energySCdefault,"photrail_energySCdefault/F");
  OutputTree[i]->Branch("pholead_energyNewCorr",&pholead_energyNewCorr,"pholead_energyNewCorr/F");
  OutputTree[i]->Branch("photrail_energyNewCorr",&photrail_energyNewCorr,"photrail_energyNewCorr/F");
  OutputTree[i]->Branch("pholead_energyNewCorrLocal",&pholead_energyNewCorrLocal,"pholead_energyNewCorrLocal/F");
  OutputTree[i]->Branch("photrail_energyNewCorrLocal",&photrail_energyNewCorrLocal,"photrail_energyNewCorrLocal/F");
 
  OutputTree[i]->Branch("pholead_SCeta",&pholead_SCeta,"pholead_SCeta/F");
  OutputTree[i]->Branch("photrail_SCeta",&photrail_SCeta,"photrail_SCeta/F");
  OutputTree[i]->Branch("pholead_SCphi",&pholead_SCphi,"pholead_SCphi/F");
  OutputTree[i]->Branch("photrail_SCphi",&photrail_SCphi,"photrail_SCphi/F");
  
  OutputTree[i]->Branch("pholead_PhoHasPixSeed",&pholead_PhoHasPixSeed,"pholead_PhoHasPixSeed/I");
  OutputTree[i]->Branch("pholead_PhoHasConvTrks",&pholead_PhoHasConvTrks,"pholead_PhoHasConvTrks/I");
  OutputTree[i]->Branch("pholead_PhoScSeedSeverity",&pholead_PhoScSeedSeverity,"pholead_PhoScSeedSeverity/I");
  
  OutputTree[i]->Branch("photrail_PhoHasPixSeed",&photrail_PhoHasPixSeed,"photrail_PhoHasPixSeed/I");
  OutputTree[i]->Branch("photrail_PhoHasConvTrks",&photrail_PhoHasConvTrks,"photrail_PhoHasConvTrks/I");
  OutputTree[i]->Branch("photrail_PhoScSeedSeverity",&photrail_PhoScSeedSeverity,"photrail_PhoScSeedSeverity/I");

  OutputTree[i]->Branch("pholead_r9",&pholead_r9,"pholead_r9/F");
  OutputTree[i]->Branch("photrail_r9",&photrail_r9,"photrail_r9/F");
  OutputTree[i]->Branch("pholead_sieie",&pholead_sieie,"pholead_sieie/F");
  OutputTree[i]->Branch("photrail_sieie",&photrail_sieie,"photrail_sieie/F");
  OutputTree[i]->Branch("pholead_hoe",&pholead_hoe,"pholead_hoe/F");
  OutputTree[i]->Branch("photrail_hoe",&photrail_hoe,"photrail_hoe/F");
  OutputTree[i]->Branch("pholead_brem",&pholead_brem,"pholead_brem/F");
  OutputTree[i]->Branch("photrail_brem",&photrail_brem,"photrail_brem/F");
  OutputTree[i]->Branch("pholead_sigmaPhi",&pholead_sigmaPhi,"pholead_sigmaPhi/F");
  OutputTree[i]->Branch("photrail_sigmaPhi",&photrail_sigmaPhi,"photrail_sigmaPhi/F");
  OutputTree[i]->Branch("pholead_sigmaEta",&pholead_sigmaEta,"pholead_sigmaEta/F");
  OutputTree[i]->Branch("photrail_sigmaEta",&photrail_sigmaEta,"photrail_sigmaEta/F");

  OutputTree[i]->Branch("pholead_pho_Cone01PhotonIso_dEta015EB_dR070EE_mvVtx",&pholead_pho_Cone01PhotonIso_dEta015EB_dR070EE_mvVtx,"pholead_pho_Cone01PhotonIso_dEta015EB_dR070EE_mvVtx/F");
  OutputTree[i]->Branch("photrail_pho_Cone01PhotonIso_dEta015EB_dR070EE_mvVtx",&photrail_pho_Cone01PhotonIso_dEta015EB_dR070EE_mvVtx,"photrail_pho_Cone01PhotonIso_dEta015EB_dR070EE_mvVtx/F");
  OutputTree[i]->Branch("pholead_pho_Cone02PhotonIso_dEta015EB_dR070EE_mvVtx",&pholead_pho_Cone02PhotonIso_dEta015EB_dR070EE_mvVtx,"pholead_pho_Cone02PhotonIso_dEta015EB_dR070EE_mvVtx/F");
  OutputTree[i]->Branch("photrail_pho_Cone02PhotonIso_dEta015EB_dR070EE_mvVtx",&photrail_pho_Cone02PhotonIso_dEta015EB_dR070EE_mvVtx,"photrail_pho_Cone02PhotonIso_dEta015EB_dR070EE_mvVtx/F");
  OutputTree[i]->Branch("pholead_pho_Cone03PhotonIso_dEta015EB_dR070EE_mvVtx",&pholead_pho_Cone03PhotonIso_dEta015EB_dR070EE_mvVtx,"pholead_pho_Cone03PhotonIso_dEta015EB_dR070EE_mvVtx/F");
  OutputTree[i]->Branch("photrail_pho_Cone03PhotonIso_dEta015EB_dR070EE_mvVtx",&photrail_pho_Cone03PhotonIso_dEta015EB_dR070EE_mvVtx,"photrail_pho_Cone03PhotonIso_dEta015EB_dR070EE_mvVtx/F");
  OutputTree[i]->Branch("pholead_pho_Cone04PhotonIso_dEta015EB_dR070EE_mvVtx",&pholead_pho_Cone04PhotonIso_dEta015EB_dR070EE_mvVtx,"pholead_pho_Cone04PhotonIso_dEta015EB_dR070EE_mvVtx/F");
  OutputTree[i]->Branch("photrail_pho_Cone04PhotonIso_dEta015EB_dR070EE_mvVtx",&photrail_pho_Cone04PhotonIso_dEta015EB_dR070EE_mvVtx,"photrail_pho_Cone04PhotonIso_dEta015EB_dR070EE_mvVtx/F");
  OutputTree[i]->Branch("pholead_pho_Cone01NeutralHadronIso_mvVtx",&pholead_pho_Cone01NeutralHadronIso_mvVtx,"pholead_pho_Cone01NeutralHadronIso_mvVtx/F");
  OutputTree[i]->Branch("photrail_pho_Cone01NeutralHadronIso_mvVtx",&photrail_pho_Cone01NeutralHadronIso_mvVtx,"photrail_pho_Cone01NeutralHadronIso_mvVtx/F");
  OutputTree[i]->Branch("pholead_pho_Cone02NeutralHadronIso_mvVtx",&pholead_pho_Cone02NeutralHadronIso_mvVtx,"pholead_pho_Cone02NeutralHadronIso_mvVtx/F");
  OutputTree[i]->Branch("photrail_pho_Cone02NeutralHadronIso_mvVtx",&photrail_pho_Cone02NeutralHadronIso_mvVtx,"photrail_pho_Cone02NeutralHadronIso_mvVtx/F");
  OutputTree[i]->Branch("pholead_pho_Cone03NeutralHadronIso_mvVtx",&pholead_pho_Cone03NeutralHadronIso_mvVtx,"pholead_pho_Cone03NeutralHadronIso_mvVtx/F");
  OutputTree[i]->Branch("photrail_pho_Cone03NeutralHadronIso_mvVtx",&photrail_pho_Cone03NeutralHadronIso_mvVtx,"photrail_pho_Cone03NeutralHadronIso_mvVtx/F");
  OutputTree[i]->Branch("pholead_pho_Cone04NeutralHadronIso_mvVtx",&pholead_pho_Cone04NeutralHadronIso_mvVtx,"pholead_pho_Cone04NeutralHadronIso_mvVtx/F");
  OutputTree[i]->Branch("photrail_pho_Cone04NeutralHadronIso_mvVtx",&photrail_pho_Cone04NeutralHadronIso_mvVtx,"photrail_pho_Cone04NeutralHadronIso_mvVtx/F");
  OutputTree[i]->Branch("pholead_pho_Cone01ChargedHadronIso_dR02_dz02_dxy01",&pholead_pho_Cone01ChargedHadronIso_dR02_dz02_dxy01,"pholead_pho_Cone01ChargedHadronIso_dR02_dz02_dxy01/F");
  OutputTree[i]->Branch("photrail_pho_Cone01ChargedHadronIso_dR02_dz02_dxy01",&photrail_pho_Cone01ChargedHadronIso_dR02_dz02_dxy01,"photrail_pho_Cone01ChargedHadronIso_dR02_dz02_dxy01/F");
  OutputTree[i]->Branch("pholead_pho_Cone02ChargedHadronIso_dR02_dz02_dxy01",&pholead_pho_Cone02ChargedHadronIso_dR02_dz02_dxy01,"pholead_pho_Cone02ChargedHadronIso_dR02_dz02_dxy01/F");
  OutputTree[i]->Branch("photrail_pho_Cone02ChargedHadronIso_dR02_dz02_dxy01",&photrail_pho_Cone02ChargedHadronIso_dR02_dz02_dxy01,"photrail_pho_Cone02ChargedHadronIso_dR02_dz02_dxy01/F");
  OutputTree[i]->Branch("pholead_pho_Cone03ChargedHadronIso_dR02_dz02_dxy01",&pholead_pho_Cone03ChargedHadronIso_dR02_dz02_dxy01,"pholead_pho_Cone03ChargedHadronIso_dR02_dz02_dxy01/F");
  OutputTree[i]->Branch("photrail_pho_Cone03ChargedHadronIso_dR02_dz02_dxy01",&photrail_pho_Cone03ChargedHadronIso_dR02_dz02_dxy01,"photrail_pho_Cone03ChargedHadronIso_dR02_dz02_dxy01/F");
  OutputTree[i]->Branch("pholead_pho_Cone04ChargedHadronIso_dR02_dz02_dxy01",&pholead_pho_Cone04ChargedHadronIso_dR02_dz02_dxy01,"pholead_pho_Cone04ChargedHadronIso_dR02_dz02_dxy01/F");
  OutputTree[i]->Branch("photrail_pho_Cone04ChargedHadronIso_dR02_dz02_dxy01",&photrail_pho_Cone04ChargedHadronIso_dR02_dz02_dxy01,"photrail_pho_Cone04ChargedHadronIso_dR02_dz02_dxy01/F");
  OutputTree[i]->Branch("pholead_pho_Cone03PFCombinedIso",&pholead_pho_Cone03PFCombinedIso,"pholead_pho_Cone03PFCombinedIso/F");
  OutputTree[i]->Branch("photrail_pho_Cone03PFCombinedIso",&photrail_pho_Cone03PFCombinedIso,"photrail_pho_Cone03PFCombinedIso/F");
  OutputTree[i]->Branch("pholead_pho_Cone04PFCombinedIso",&pholead_pho_Cone04PFCombinedIso,"pholead_pho_Cone04PFCombinedIso/F");
  OutputTree[i]->Branch("photrail_pho_Cone04PFCombinedIso",&photrail_pho_Cone04PFCombinedIso,"photrail_pho_Cone04PFCombinedIso/F");

  OutputTree[i]->Branch("pholead_PhoPassConvSafeElectronVeto",&pholead_PhoPassConvSafeElectronVeto,"pholead_PhoPassConvSafeElectronVeto/I");
  OutputTree[i]->Branch("photrail_PhoPassConvSafeElectronVeto",&photrail_PhoPassConvSafeElectronVeto,"photrail_PhoPassConvSafeElectronVeto/I");

  OutputTree[i]->Branch("pholead_GenPhotonIsoDR04",&pholead_GenPhotonIsoDR04,"pholead_GenPhotonIsoDR04/F");
  OutputTree[i]->Branch("photrail_GenPhotonIsoDR04",&photrail_GenPhotonIsoDR04,"photrail_GenPhotonIsoDR04/F");

  OutputTree[i]->Branch("pholead_PhoIso03Ecal",&pholead_PhoIso03Ecal,"pholead_PhoIso03Ecal/F");
  OutputTree[i]->Branch("pholead_PhoIso03Hcal",&pholead_PhoIso03Hcal,"pholead_PhoIso03Hcal/F");
  OutputTree[i]->Branch("pholead_PhoIso03TrkSolid",&pholead_PhoIso03TrkSolid,"pholead_PhoIso03TrkSolid/F");
  OutputTree[i]->Branch("pholead_PhoIso03TrkHollow",&pholead_PhoIso03TrkHollow,"pholead_PhoIso03TrkHollow/F");
  OutputTree[i]->Branch("pholead_PhoIso03",&pholead_PhoIso03,"pholead_PhoIso03/F");
  OutputTree[i]->Branch("pholead_PhoIso04Ecal",&pholead_PhoIso04Ecal,"pholead_PhoIso04Ecal/F");
  OutputTree[i]->Branch("pholead_PhoIso04Hcal",&pholead_PhoIso04Hcal,"pholead_PhoIso04Hcal/F");
  OutputTree[i]->Branch("pholead_PhoIso04TrkSolid",&pholead_PhoIso04TrkSolid,"pholead_PhoIso04TrkSolid/F");
  OutputTree[i]->Branch("pholead_PhoIso04TrkHollow",&pholead_PhoIso04TrkHollow,"pholead_PhoIso04TrkHollow/F");
  OutputTree[i]->Branch("pholead_PhoIso04",&pholead_PhoIso04,"pholead_PhoIso04/F");


  OutputTree[i]->Branch("photrail_PhoIso03Ecal",&photrail_PhoIso03Ecal,"photrail_PhoIso03Ecal/F");
  OutputTree[i]->Branch("photrail_PhoIso03Hcal",&photrail_PhoIso03Hcal,"photrail_PhoIso03Hcal/F");
  OutputTree[i]->Branch("photrail_PhoIso03TrkSolid",&photrail_PhoIso03TrkSolid,"photrail_PhoIso03TrkSolid/F");
  OutputTree[i]->Branch("photrail_PhoIso03TrkHollow",&photrail_PhoIso03TrkHollow,"photrail_PhoIso03TrkHollow/F");
  OutputTree[i]->Branch("photrail_PhoIso03",&photrail_PhoIso03,"photrail_PhoIso03/F");
  OutputTree[i]->Branch("photrail_PhoIso04Ecal",&photrail_PhoIso04Ecal,"photrail_PhoIso04Ecal/F");
  OutputTree[i]->Branch("photrail_PhoIso04Hcal",&photrail_PhoIso04Hcal,"photrail_PhoIso04Hcal/F");
  OutputTree[i]->Branch("photrail_PhoIso04TrkSolid",&photrail_PhoIso04TrkSolid,"photrail_PhoIso04TrkSolid/F");
  OutputTree[i]->Branch("photrail_PhoIso04TrkHollow",&photrail_PhoIso04TrkHollow,"photrail_PhoIso04TrkHollow/F");
  OutputTree[i]->Branch("photrail_PhoIso04",&photrail_PhoIso04,"photrail_PhoIso04/F");


  OutputTree[i]->Branch("pholead_PhoS4OverS1",&pholead_PhoS4OverS1,"pholead_PhoS4OverS1/F");
  OutputTree[i]->Branch("pholead_PhoSigmaEtaEta",&pholead_PhoSigmaEtaEta,"pholead_PhoSigmaEtaEta/F");
  OutputTree[i]->Branch("pholead_PhoE1x5",&pholead_PhoE1x5,"pholead_PhoE1x5/F");
  OutputTree[i]->Branch("pholead_PhoE2x5",&pholead_PhoE2x5,"pholead_PhoE2x5/F");
  OutputTree[i]->Branch("pholead_PhoE3x3",&pholead_PhoE3x3,"pholead_PhoE3x3/F");
  OutputTree[i]->Branch("pholead_PhoE5x5",&pholead_PhoE5x5,"pholead_PhoE5x5/F");
  OutputTree[i]->Branch("pholead_PhomaxEnergyXtal",&pholead_PhomaxEnergyXtal,"pholead_PhomaxEnergyXtal/F");
  OutputTree[i]->Branch("pholead_PhoIso03HcalDepth1",&pholead_PhoIso03HcalDepth1,"pholead_PhoIso03HcalDepth1/F");
  OutputTree[i]->Branch("pholead_PhoIso03HcalDepth2",&pholead_PhoIso03HcalDepth2,"pholead_PhoIso03HcalDepth2/F");
  OutputTree[i]->Branch("pholead_PhoIso04HcalDepth1",&pholead_PhoIso04HcalDepth1,"pholead_PhoIso04HcalDepth1/F");
  OutputTree[i]->Branch("pholead_PhoIso04HcalDepth2",&pholead_PhoIso04HcalDepth2,"pholead_PhoIso04HcalDepth2/F");
  OutputTree[i]->Branch("pholead_PhoIso03nTrksSolid",&pholead_PhoIso03nTrksSolid,"pholead_PhoIso03nTrksSolid/I");
  OutputTree[i]->Branch("pholead_PhoIso03nTrksHollow",&pholead_PhoIso03nTrksHollow,"pholead_PhoIso03nTrksHollow/I");
  OutputTree[i]->Branch("pholead_PhoIso04nTrksSolid",&pholead_PhoIso04nTrksSolid,"pholead_PhoIso04nTrksSolid/I");
  OutputTree[i]->Branch("pholead_PhoIso04nTrksHollow",&pholead_PhoIso04nTrksHollow,"pholead_PhoIso04nTrksHollow/I");
  OutputTree[i]->Branch("pholead_Pho_ChargedHadronIso",&pholead_Pho_ChargedHadronIso,"pholead_Pho_ChargedHadronIso/F");
  OutputTree[i]->Branch("pholead_Pho_NeutralHadronIso",&pholead_Pho_NeutralHadronIso,"pholead_Pho_NeutralHadronIso/F");
  OutputTree[i]->Branch("pholead_Pho_PhotonIso",&pholead_Pho_PhotonIso,"pholead_Pho_PhotonIso/F");
  OutputTree[i]->Branch("pholead_Pho_isPFPhoton",&pholead_Pho_isPFPhoton,"pholead_Pho_isPFPhoton/I");
  OutputTree[i]->Branch("pholead_Pho_isPFElectron",&pholead_Pho_isPFElectron,"pholead_Pho_isPFElectron/I");

  OutputTree[i]->Branch("photrail_PhoS4OverS1",&photrail_PhoS4OverS1,"photrail_PhoS4OverS1/F");
  OutputTree[i]->Branch("photrail_PhoSigmaEtaEta",&photrail_PhoSigmaEtaEta,"photrail_PhoSigmaEtaEta/F");
  OutputTree[i]->Branch("photrail_PhoE1x5",&photrail_PhoE1x5,"photrail_PhoE1x5/F");
  OutputTree[i]->Branch("photrail_PhoE2x5",&photrail_PhoE2x5,"photrail_PhoE2x5/F");
  OutputTree[i]->Branch("photrail_PhoE3x3",&photrail_PhoE3x3,"photrail_PhoE3x3/F");
  OutputTree[i]->Branch("photrail_PhoE5x5",&photrail_PhoE5x5,"photrail_PhoE5x5/F");
  OutputTree[i]->Branch("photrail_PhomaxEnergyXtal",&photrail_PhomaxEnergyXtal,"photrail_PhomaxEnergyXtal/F");
  OutputTree[i]->Branch("photrail_PhoIso03HcalDepth1",&photrail_PhoIso03HcalDepth1,"photrail_PhoIso03HcalDepth1/F");
  OutputTree[i]->Branch("photrail_PhoIso03HcalDepth2",&photrail_PhoIso03HcalDepth2,"photrail_PhoIso03HcalDepth2/F");
  OutputTree[i]->Branch("photrail_PhoIso04HcalDepth1",&photrail_PhoIso04HcalDepth1,"photrail_PhoIso04HcalDepth1/F");
  OutputTree[i]->Branch("photrail_PhoIso04HcalDepth2",&photrail_PhoIso04HcalDepth2,"photrail_PhoIso04HcalDepth2/F");
  OutputTree[i]->Branch("photrail_PhoIso03nTrksSolid",&photrail_PhoIso03nTrksSolid,"photrail_PhoIso03nTrksSolid/I");
  OutputTree[i]->Branch("photrail_PhoIso03nTrksHollow",&photrail_PhoIso03nTrksHollow,"photrail_PhoIso03nTrksHollow/I");
  OutputTree[i]->Branch("photrail_PhoIso04nTrksSolid",&photrail_PhoIso04nTrksSolid,"photrail_PhoIso04nTrksSolid/I");
  OutputTree[i]->Branch("photrail_PhoIso04nTrksHollow",&photrail_PhoIso04nTrksHollow,"photrail_PhoIso04nTrksHollow/I");
  OutputTree[i]->Branch("photrail_Pho_ChargedHadronIso",&photrail_Pho_ChargedHadronIso,"photrail_Pho_ChargedHadronIso/F");
  OutputTree[i]->Branch("photrail_Pho_NeutralHadronIso",&photrail_Pho_NeutralHadronIso,"photrail_Pho_NeutralHadronIso/F");
  OutputTree[i]->Branch("photrail_Pho_PhotonIso",&photrail_Pho_PhotonIso,"photrail_Pho_PhotonIso/F");
  OutputTree[i]->Branch("photrail_Pho_isPFPhoton",&photrail_Pho_isPFPhoton,"photrail_Pho_isPFPhoton/I");
  OutputTree[i]->Branch("photrail_Pho_isPFElectron",&photrail_Pho_isPFElectron,"photrail_Pho_isPFElectron/I");

  OutputTree[i]->Branch("pholead_PhoMCmatchindex",&pholead_PhoMCmatchindex,"pholead_PhoMCmatchindex/I");
  OutputTree[i]->Branch("pholead_PhoMCmatchexitcode",&pholead_PhoMCmatchexitcode,"pholead_PhoMCmatchexitcode/I");

  OutputTree[i]->Branch("photrail_PhoMCmatchindex",&photrail_PhoMCmatchindex,"photrail_PhoMCmatchindex/I");
  OutputTree[i]->Branch("photrail_PhoMCmatchexitcode",&photrail_PhoMCmatchexitcode,"photrail_PhoMCmatchexitcode/I");

  }


  fHNumPU = new TH1F("NumPU","NumPU",40,0,40);
  fHNumVtx = new TH1F("NumVtx","NumVtx",40,0,40);
	

	
  cout << "Trees and histos created" << endl;



}

void DiPhotonMiniTree::Analyze(){

  //    cout << "Analyze this event" << endl;



    //    cout << "A" << endl;

  float weight;
  if (!isdata) weight = GetPUWeight(fTR->PUnumInteractions);
  else weight=1;
  
  event_luminormfactor=AddWeight;

  if (!isdata) fHNumPU->Fill(fTR->PUnumInteractions,weight);
  fHNumVtx->Fill(fTR->NVrtx,weight);

  event_weight = weight;
  event_rho = fTR->Rho;
  if (!isdata) event_nPU = fTR->PUnumInteractions;
  event_nRecVtx = fTR->NVrtx;

  //  cout << "B" << endl;

  int nmatched_part_isrfsr_gamma=0;
  {
    std::vector<int> temp;
    for (int i=0; i<fTR->NPhotons; i++){
      temp.push_back(i);
    }
    //    cout << "B1" << endl;
    temp = PhotonPreSelection(fTR,temp);
    //    cout << "B2" << endl;
    temp = ApplyPixelVeto(fTR,temp,0);
    //    cout << "B3" << endl;
    nmatched_part_isrfsr_gamma = Count_part_isrfsr_gamma(fTR,temp);
    //    cout << "B4" << endl;
  }
  if (!isdata) event_Kfactor = kfactors[2-nmatched_part_isrfsr_gamma];
  else event_Kfactor=1;

  //  cout << "C" << endl;

  // 0 = standard selection for data and MC
  // 1 = signal template generation from MC
  // 2 = background template generation from MC
  // 3 = DY selection (standard with reversed pixel veto)

  if (!TriggerSelection()) return;

  std::vector<int> passing_selection[4];

  bool pass[4];

  for (int sel_cat=0; sel_cat<4; sel_cat++){

    pass[sel_cat]=false;

    std::vector<int> passing;

    for (int i=0; i<fTR->NPhotons; i++){
      passing.push_back(i);
    }

    passing = PhotonPreSelection(fTR,passing);

    if (sel_cat!=3) passing = ApplyPixelVeto(fTR,passing,0);
    if (sel_cat==3) passing = ApplyPixelVeto(fTR,passing,1);

    passing = PhotonSelection(fTR,passing);

    if (sel_cat==0 || sel_cat==3){
      pass[sel_cat] = StandardEventSelection(fTR,passing);
    }
    if (sel_cat==1 || sel_cat==2){
      if (sel_cat==1) passing = SignalSelection(fTR,passing);
      if (sel_cat==2) passing = BackgroundSelection(fTR,passing);
      pass[sel_cat] = SinglePhotonEventSelection(fTR,passing);
    }

    passing_selection[sel_cat] = passing;

  }

  //  cout << "D" << endl;

  for (int sel_cat=0; sel_cat<4; sel_cat++){

    if (!pass[sel_cat]) continue;

    std::vector<int> passing = passing_selection[sel_cat];

    int minsize=999;
    if (sel_cat==0 || sel_cat==3) minsize=2;
    if (sel_cat==1 || sel_cat==2) minsize=1;

    if (!(passing_selection[sel_cat].size()>=minsize)){
      std::cout << "Error!!!" << std::endl;
      continue;
    }

    ResetVars();  
  
  pholead_eta = fTR->PhoEta[passing.at(0)];
  pholead_px = fTR->PhoPx[passing.at(0)];
  pholead_py = fTR->PhoPy[passing.at(0)];
  pholead_pt = fTR->PhoPt[passing.at(0)];
  pholead_pz = fTR->PhoPz[passing.at(0)];
  pholead_energy = fTR->PhoEnergy[passing.at(0)];
  pholead_SCeta = fTR->SCEta[fTR->PhotSCindex[passing.at(0)]];
  pholead_SCphi = fTR->SCPhi[fTR->PhotSCindex[passing.at(0)]];
  pholead_PhoHasPixSeed=fTR->PhoHasPixSeed[passing.at(0)];
  pholead_PhoHasConvTrks=fTR->PhoHasConvTrks[passing.at(0)];
  pholead_PhoScSeedSeverity=fTR->PhoScSeedSeverity[passing.at(0)];
  pholead_energySCdefault = CorrPhoton(fTR,passing.at(0),0).E();
  pholead_energyNewCorr = CorrPhoton(fTR,passing.at(0),5).E();
  pholead_energyNewCorrLocal = CorrPhoton(fTR,passing.at(0),6).E();
  pholead_r9 = fTR->SCR9[fTR->PhotSCindex[passing.at(0)]];
  pholead_sieie = fTR->PhoSigmaIetaIeta[passing.at(0)];
  pholead_hoe = fTR->PhoHoverE[passing.at(0)];
  pholead_brem = fTR->SCBrem[fTR->PhotSCindex[passing.at(0)]];
  pholead_sigmaPhi = fTR->SCPhiWidth[fTR->PhotSCindex[passing.at(0)]];
  pholead_sigmaEta = fTR->SCEtaWidth[fTR->PhotSCindex[passing.at(0)]];
  pholead_pho_Cone01PhotonIso_dEta015EB_dR070EE_mvVtx=fTR->pho_Cone01PhotonIso_dEta015EB_dR070EE_mvVtx[passing.at(0)];
  pholead_pho_Cone02PhotonIso_dEta015EB_dR070EE_mvVtx=fTR->pho_Cone02PhotonIso_dEta015EB_dR070EE_mvVtx[passing.at(0)];
  pholead_pho_Cone03PhotonIso_dEta015EB_dR070EE_mvVtx=fTR->pho_Cone03PhotonIso_dEta015EB_dR070EE_mvVtx[passing.at(0)];
  pholead_pho_Cone04PhotonIso_dEta015EB_dR070EE_mvVtx=fTR->pho_Cone04PhotonIso_dEta015EB_dR070EE_mvVtx[passing.at(0)];
  pholead_pho_Cone01NeutralHadronIso_mvVtx=fTR->pho_Cone01NeutralHadronIso_mvVtx[passing.at(0)];
  pholead_pho_Cone02NeutralHadronIso_mvVtx=fTR->pho_Cone02NeutralHadronIso_mvVtx[passing.at(0)];
  pholead_pho_Cone03NeutralHadronIso_mvVtx=fTR->pho_Cone03NeutralHadronIso_mvVtx[passing.at(0)];
  pholead_pho_Cone04NeutralHadronIso_mvVtx=fTR->pho_Cone04NeutralHadronIso_mvVtx[passing.at(0)];
  pholead_pho_Cone01ChargedHadronIso_dR02_dz02_dxy01=fTR->pho_Cone01ChargedHadronIso_dR02_dz02_dxy01[passing.at(0)];
  pholead_pho_Cone02ChargedHadronIso_dR02_dz02_dxy01=fTR->pho_Cone02ChargedHadronIso_dR02_dz02_dxy01[passing.at(0)];
  pholead_pho_Cone03ChargedHadronIso_dR02_dz02_dxy01=fTR->pho_Cone03ChargedHadronIso_dR02_dz02_dxy01[passing.at(0)];
  pholead_pho_Cone04ChargedHadronIso_dR02_dz02_dxy01=fTR->pho_Cone04ChargedHadronIso_dR02_dz02_dxy01[passing.at(0)];
  pholead_pho_Cone03PFCombinedIso=fTR->pho_Cone03PFCombinedIso[passing.at(0)];
  pholead_pho_Cone04PFCombinedIso=fTR->pho_Cone04PFCombinedIso[passing.at(0)];
  pholead_PhoPassConvSafeElectronVeto=fTR->PhoPassConvSafeElectronVeto[passing.at(0)];
  if (fTR->PhoMCmatchindex[passing.at(0)]>=0) pholead_GenPhotonIsoDR04=fTR->GenPhotonIsoDR04[fTR->PhoMCmatchindex[passing.at(0)]];
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
  pholead_PhoMCmatchindex=fTR->PhoMCmatchindex[passing.at(0)];
  pholead_PhoMCmatchexitcode=fTR->PhoMCmatchexitcode[passing.at(0)];

  if (sel_cat==0 || sel_cat==3){
  photrail_eta = fTR->PhoEta[passing.at(1)];
  photrail_px = fTR->PhoPx[passing.at(1)];
  photrail_py = fTR->PhoPy[passing.at(1)];
  photrail_pt = fTR->PhoPt[passing.at(1)];
  photrail_pz = fTR->PhoPz[passing.at(1)];
  photrail_energy = fTR->PhoEnergy[passing.at(1)];
  photrail_SCeta = fTR->SCEta[fTR->PhotSCindex[passing.at(1)]];
  photrail_SCphi = fTR->SCPhi[fTR->PhotSCindex[passing.at(1)]];
  photrail_PhoHasPixSeed=fTR->PhoHasPixSeed[passing.at(1)];
  photrail_PhoHasConvTrks=fTR->PhoHasConvTrks[passing.at(1)];
  photrail_PhoScSeedSeverity=fTR->PhoScSeedSeverity[passing.at(1)];
  photrail_energySCdefault = CorrPhoton(fTR,passing.at(1),0).E();
  photrail_energyNewCorr = CorrPhoton(fTR,passing.at(1),5).E();
  photrail_energyNewCorrLocal = CorrPhoton(fTR,passing.at(1),6).E();
  photrail_r9 = fTR->SCR9[fTR->PhotSCindex[passing.at(1)]];
  photrail_sieie = fTR->PhoSigmaIetaIeta[passing.at(1)];
  photrail_hoe = fTR->PhoHoverE[passing.at(1)];
  photrail_brem = fTR->SCBrem[fTR->PhotSCindex[passing.at(1)]];
  photrail_sigmaPhi = fTR->SCPhiWidth[fTR->PhotSCindex[passing.at(1)]];
  photrail_sigmaEta = fTR->SCEtaWidth[fTR->PhotSCindex[passing.at(1)]];
  photrail_pho_Cone01PhotonIso_dEta015EB_dR070EE_mvVtx=fTR->pho_Cone01PhotonIso_dEta015EB_dR070EE_mvVtx[passing.at(1)];
  photrail_pho_Cone02PhotonIso_dEta015EB_dR070EE_mvVtx=fTR->pho_Cone02PhotonIso_dEta015EB_dR070EE_mvVtx[passing.at(1)];
  photrail_pho_Cone03PhotonIso_dEta015EB_dR070EE_mvVtx=fTR->pho_Cone03PhotonIso_dEta015EB_dR070EE_mvVtx[passing.at(1)];
  photrail_pho_Cone04PhotonIso_dEta015EB_dR070EE_mvVtx=fTR->pho_Cone04PhotonIso_dEta015EB_dR070EE_mvVtx[passing.at(1)];
  photrail_pho_Cone01NeutralHadronIso_mvVtx=fTR->pho_Cone01NeutralHadronIso_mvVtx[passing.at(1)];
  photrail_pho_Cone02NeutralHadronIso_mvVtx=fTR->pho_Cone02NeutralHadronIso_mvVtx[passing.at(1)];
  photrail_pho_Cone03NeutralHadronIso_mvVtx=fTR->pho_Cone03NeutralHadronIso_mvVtx[passing.at(1)];
  photrail_pho_Cone04NeutralHadronIso_mvVtx=fTR->pho_Cone04NeutralHadronIso_mvVtx[passing.at(1)];
  photrail_pho_Cone01ChargedHadronIso_dR02_dz02_dxy01=fTR->pho_Cone01ChargedHadronIso_dR02_dz02_dxy01[passing.at(1)];
  photrail_pho_Cone02ChargedHadronIso_dR02_dz02_dxy01=fTR->pho_Cone02ChargedHadronIso_dR02_dz02_dxy01[passing.at(1)];
  photrail_pho_Cone03ChargedHadronIso_dR02_dz02_dxy01=fTR->pho_Cone03ChargedHadronIso_dR02_dz02_dxy01[passing.at(1)];
  photrail_pho_Cone04ChargedHadronIso_dR02_dz02_dxy01=fTR->pho_Cone04ChargedHadronIso_dR02_dz02_dxy01[passing.at(1)];
  photrail_pho_Cone03PFCombinedIso=fTR->pho_Cone03PFCombinedIso[passing.at(1)];
  photrail_pho_Cone04PFCombinedIso=fTR->pho_Cone04PFCombinedIso[passing.at(1)];
  photrail_PhoPassConvSafeElectronVeto=fTR->PhoPassConvSafeElectronVeto[passing.at(1)];
  if (fTR->PhoMCmatchindex[passing.at(1)]>=0) photrail_GenPhotonIsoDR04=fTR->GenPhotonIsoDR04[fTR->PhoMCmatchindex[passing.at(1)]];
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
  photrail_PhoMCmatchindex=fTR->PhoMCmatchindex[passing.at(1)];
  photrail_PhoMCmatchexitcode=fTR->PhoMCmatchexitcode[passing.at(1)];
  }

  if (sel_cat==0 || sel_cat==3) {
    float invmass0 = (CorrPhoton(fTR,passing.at(0),0)+CorrPhoton(fTR,passing.at(1),0)).M();
    float invmass5 = (CorrPhoton(fTR,passing.at(0),5)+CorrPhoton(fTR,passing.at(1),5)).M();
    float invmass6 = (CorrPhoton(fTR,passing.at(0),6)+CorrPhoton(fTR,passing.at(1),6)).M();
    dipho_mgg_photon = invmass0;
    dipho_mgg_newCorr = invmass5;
    dipho_mgg_newCorrLocal = invmass6;
  }
  


  OutputTree[sel_cat]->Fill();
 
  }
 
  //  cout << "E" << endl;

}

void DiPhotonMiniTree::End(){
  fOutputFile->cd();
  for (int i=0; i<4; i++) OutputTree[i]->Write();	
  fHNumPU->Write();
  fHNumVtx->Write();
	
  fOutputFile->Close();

}

TLorentzVector DiPhotonMiniTree::CorrPhoton(TreeReader *fTR, int i, int mode){
  TLorentzVector corr;
  corr.SetPtEtaPhiE(fTR->PhoPt[i],fTR->PhoEta[i],fTR->PhoPhi[i],fTR->PhoEnergy[i]);
  float corrE=phocorr->get_correctedenergy(fTR,i,mode);
  corr.SetE(corrE);
  corr.SetRho(corrE);
  return corr;
};

std::vector<int> DiPhotonMiniTree::ApplyPixelVeto(TreeReader *fTR, vector<int> passing, bool forelectron){

  for (vector<int>::iterator it = passing.begin(); it != passing.end(); ){ // Pixel veto (conversion safe)
    int wantpixelseed;
    if (forelectron) wantpixelseed=1; else wantpixelseed=0; 
    if (fTR->PhoPassConvSafeElectronVeto[*it]==wantpixelseed) it=passing.erase(it); else it++;
  }

  return passing;

};

std::vector<int> DiPhotonMiniTree::PhotonPreSelection(TreeReader *fTR, std::vector<int> passing){

  for (vector<int>::iterator it = passing.begin(); it != passing.end(); ){ // SC matching
    if (fTR->PhotSCindex[*it]==-1) it=passing.erase(it); else it++;
  }

  for (vector<int>::iterator it = passing.begin(); it != passing.end(); ){ // fiducial region
    float eta=fTR->SCEta[fTR->PhotSCindex[*it]];
    if ((fabs(eta)>1.4442 && fabs(eta)<1.56) || (fabs(eta)>2.5)) it=passing.erase(it); else it++;
  }

  for (vector<int>::iterator it = passing.begin(); it != passing.end(); ){ // Pt cut on RawEnCetaCorr
    float energy=fTR->SCRaw[fTR->PhotSCindex[*it]];
    float eta=fTR->SCEta[fTR->PhotSCindex[*it]];
    if (fabs(eta)<1.4442) energy*=phocorr->getEtaCorrectionBarrel(eta);
    if (fabs(eta)>1.56) energy+=fTR->SCPre[fTR->PhotSCindex[*it]];
    if (energy/cosh(eta)<20) it=passing.erase(it); else it++;
  }

//  for (vector<int>::iterator it = passing.begin(); it != passing.end(); ){ // Remove eta/phi cracks
//    float eta=fTR->SCEta[fTR->PhotSCindex[*it]];
//    float phi=fTR->SCPhi[fTR->PhotSCindex[*it]];
//    if ((phocorr->isInPhiCracks(phi,eta)) || (phocorr->isInEBEtaCracks(eta))) it=passing.erase(it); else it++;
//  }

  // MVA presel from Hgg

  for (vector<int>::iterator it = passing.begin(); it != passing.end(); ){ // HoverE cut
    float r9=fTR->SCR9[fTR->PhotSCindex[*it]];
    float eta=fTR->SCEta[fTR->PhotSCindex[*it]];
    float hoe=fTR->PhoHoverE[*it];
    bool pass=0;
    if (fabs(eta)<1.4442 && r9<0.9 && hoe<0.075) pass=1;
    if (fabs(eta)<1.4442 && r9>0.9 && hoe<0.082) pass=1;
    if (fabs(eta)>1.56 && r9<0.9 && hoe<0.075) pass=1;
    if (fabs(eta)>1.56 && r9>0.9 && hoe<0.075) pass=1;
    if (!pass) it=passing.erase(it); else it++;
  }
	
  for (vector<int>::iterator it = passing.begin(); it != passing.end(); ){ // sieie cut
    float eta=fTR->SCEta[fTR->PhotSCindex[*it]];
    float sieie=fTR->PhoSigmaIetaIeta[*it];
    bool pass=0;
    if (fabs(eta)<1.4442 && sieie<0.014 && sieie>0.001) pass=1; // to add sigmaiphiphi>0.001 in the future
    if (fabs(eta)>1.56 && sieie<0.034) pass=1;
    if (!pass) it=passing.erase(it); else it++;
  }
	
  for (vector<int>::iterator it = passing.begin(); it != passing.end(); ){ // isolation cuts (trigger)
    float r9=fTR->SCR9[fTR->PhotSCindex[*it]];
    bool pass=0;
    float etcorrecaliso=fTR->PhoIso03Ecal[*it]-0.012*fTR->PhoPt[*it];
    float etcorrhcaliso=fTR->PhoIso03Hcal[*it]-0.005*fTR->PhoPt[*it];
    float etcorrtrkiso=fTR->PhoIso03TrkHollow[*it]-0.002*fTR->PhoPt[*it];
    if (r9<0.9 && etcorrecaliso<4 && etcorrhcaliso<4 && etcorrtrkiso<4) pass=1;
    if (r9>0.9 && etcorrecaliso<50 && etcorrhcaliso<50 && etcorrtrkiso<50) pass=1;
    if (!pass) it=passing.erase(it); else it++;
  }

  for (vector<int>::iterator it = passing.begin(); it != passing.end(); ){ // isolation cuts (filter)
    float r9=fTR->SCR9[fTR->PhotSCindex[*it]];
    bool pass=0;
    if(fTR->pho_Cone02ChargedHadronIso_dR02_dz02_dxy01[*it]<4) pass=1;
    if (!pass) it=passing.erase(it); else it++;
  }

  return passing;

};

std::vector<int> DiPhotonMiniTree::GenLevelIsolationCut(TreeReader *fTR, std::vector<int> passing){

  // Genlevel isolation cut for MC
  // selects only photons matched to gen level, with gen-level iso < 5 GeV

  if (isdata) {
    std::cout << "Calling gen level isolation cut on data!!!" << std::endl;
    return passing;
  }
  
  for (std::vector<int>::iterator it = passing.begin(); it != passing.end(); ){
      bool pass=0;
      if (fTR->PhoMCmatchexitcode[*it]==1 || fTR->PhoMCmatchexitcode[*it]==2)
	if(fTR->GenPhotonIsoDR04[fTR->PhoMCmatchindex[*it]]<5)
	  pass=1;
      if (!pass) it=passing.erase(it); else it++;
    }

  return passing;

};

std::vector<int> DiPhotonMiniTree::PhotonSelection(TreeReader *fTR, std::vector<int> passing){

  for (vector<int>::iterator it = passing.begin(); it != passing.end(); ){ // HoverE cut
    float eta=fTR->SCEta[fTR->PhotSCindex[*it]];
    float hoe=fTR->PhoHoverE[*it];
    bool pass=0;
    if (fabs(eta)<1.4442 && hoe<0.05) pass=1;
    if (fabs(eta)>1.56 && hoe<0.05) pass=1;
    if (!pass) it=passing.erase(it); else it++;
  }
	
  for (vector<int>::iterator it = passing.begin(); it != passing.end(); ){ // sieie cut
    float eta=fTR->SCEta[fTR->PhotSCindex[*it]];
    float sieie=fTR->PhoSigmaIetaIeta[*it];
    bool pass=0;
    if (fabs(eta)<1.4442 && sieie<0.011) pass=1;
    if (fabs(eta)>1.56 && sieie<0.030) pass=1;
    if (!pass) it=passing.erase(it); else it++;
  }

  for (vector<int>::iterator it = passing.begin(); it != passing.end(); ){
    bool pass=0;
    if (fTR->pho_Cone04PFCombinedIso[*it]*fTR->PhoPt[*it]<5) pass=1;
    if (!pass) it=passing.erase(it); else it++;
  }

  return passing;

};

std::vector<int> DiPhotonMiniTree::SignalSelection(TreeReader *fTR, std::vector<int> passing){

  for (vector<int>::iterator it = passing.begin(); it != passing.end(); ){
    bool pass=0;
    if (fTR->PhoMCmatchexitcode[*it]==1 || fTR->PhoMCmatchexitcode[*it]==2) pass=1;
    if (!pass) it=passing.erase(it); else it++;
  }

  std::vector<int> passing2 = GenLevelIsolationCut(fTR,passing);

  return passing2;

};

std::vector<int> DiPhotonMiniTree::BackgroundSelection(TreeReader *fTR, std::vector<int> passing){

  for (vector<int>::iterator it = passing.begin(); it != passing.end(); ){
    bool pass=0;
    if (!(fTR->PhoMCmatchexitcode[*it]==1 || fTR->PhoMCmatchexitcode[*it]==2)) pass=1;
    if (!pass) it=passing.erase(it); else it++;
  }

  return passing;

};


bool DiPhotonMiniTree::SinglePhotonEventSelection(TreeReader *fTR, std::vector<int> passing){

  if (passing.size()<1) return false;

  if (fTR->PhoPt[passing.at(0)]<40) return false;

  return true;

};

bool DiPhotonMiniTree::StandardEventSelection(TreeReader *fTR, std::vector<int> passing){

  if (passing.size()<2) return false;

  float invmass0 = (CorrPhoton(fTR,passing.at(0),0)+CorrPhoton(fTR,passing.at(1),0)).M();
  float deta=fTR->PhoEta[passing.at(0)]-fTR->PhoEta[passing.at(1)];
  float dphi=fTR->PhoPhi[passing.at(0)]-fTR->PhoPhi[passing.at(1)];
  if (dphi>TMath::Pi()) dphi=2*TMath::Pi()-dphi;
  if (dphi<-TMath::Pi()) dphi=-2*TMath::Pi()-dphi;
  double dR=sqrt(dphi*dphi+deta*deta);

  if (fTR->PhoPt[passing.at(0)]<40) return false;
  if (fTR->PhoPt[passing.at(1)]<30) return false;
  if (invmass0<80) return false;
  if (dR<0.4) return false;

  return true;
};

bool DiPhotonMiniTree::TriggerSelection(){
#include "DiPhotonTriggerSelection.cc"
};

int DiPhotonMiniTree::Count_part_isrfsr_gamma(TreeReader *fTR, std::vector<int> passing){

  if (isdata){
    std::cout << "calling count of gammas for k-factor on data!!!" << std::endl;
    return 0;
  }

  int res=0;
  for (int i=0; i<2 && i<passing.size(); i++){
    int code = fTR->PhoMCmatchexitcode[passing.at(i)];
    if (code==1 || code==2) res++;
  }
  return res;

};

void DiPhotonMiniTree::ResetVars(){
  dipho_mgg_photon = -999;
  dipho_mgg_newCorr = -999;
  dipho_mgg_newCorrLocal = -999;
  pholead_eta = -999;
  photrail_eta = -999;
  pholead_px = -999;
  photrail_px = -999;
  pholead_py = -999;
  photrail_py = -999;
  pholead_pt = -999;
  photrail_pt = -999;
  pholead_pz = -999;
  photrail_pz = -999;
  pholead_energy = -999;
  photrail_energy = -999;
  pholead_energySCdefault = -999;
  photrail_energySCdefault = -999;
  pholead_energyNewCorr = -999;
  photrail_energyNewCorr = -999;
  pholead_energyNewCorrLocal = -999;
  photrail_energyNewCorrLocal = -999;
  pholead_SCeta = -999;
  photrail_SCeta = -999;
  pholead_SCphi = -999;
  photrail_SCphi = -999;
  pholead_PhoHasPixSeed = -999;
  pholead_PhoHasConvTrks = -999;
  pholead_PhoScSeedSeverity = -999;
  photrail_PhoHasPixSeed = -999;
  photrail_PhoHasConvTrks = -999;
  photrail_PhoScSeedSeverity = -999;
  pholead_r9 = -999;
  photrail_r9 = -999;
  pholead_sieie = -999;
  photrail_sieie = -999;
  pholead_hoe = -999;
  photrail_hoe = -999;
  pholead_brem = -999;
  photrail_brem = -999;
  pholead_sigmaPhi = -999;
  photrail_sigmaPhi = -999;
  pholead_sigmaEta = -999;
  photrail_sigmaEta = -999;
  pholead_pho_Cone01PhotonIso_dEta015EB_dR070EE_mvVtx = -999;
  photrail_pho_Cone01PhotonIso_dEta015EB_dR070EE_mvVtx = -999;
  pholead_pho_Cone02PhotonIso_dEta015EB_dR070EE_mvVtx = -999;
  photrail_pho_Cone02PhotonIso_dEta015EB_dR070EE_mvVtx = -999;
  pholead_pho_Cone03PhotonIso_dEta015EB_dR070EE_mvVtx = -999;
  photrail_pho_Cone03PhotonIso_dEta015EB_dR070EE_mvVtx = -999;
  pholead_pho_Cone04PhotonIso_dEta015EB_dR070EE_mvVtx = -999;
  photrail_pho_Cone04PhotonIso_dEta015EB_dR070EE_mvVtx = -999;
  pholead_pho_Cone01NeutralHadronIso_mvVtx = -999;
  photrail_pho_Cone01NeutralHadronIso_mvVtx = -999;
  pholead_pho_Cone02NeutralHadronIso_mvVtx = -999;
  photrail_pho_Cone02NeutralHadronIso_mvVtx = -999;
  pholead_pho_Cone03NeutralHadronIso_mvVtx = -999;
  photrail_pho_Cone03NeutralHadronIso_mvVtx = -999;
  pholead_pho_Cone04NeutralHadronIso_mvVtx = -999;
  photrail_pho_Cone04NeutralHadronIso_mvVtx = -999;
  pholead_pho_Cone01ChargedHadronIso_dR02_dz02_dxy01 = -999;
  photrail_pho_Cone01ChargedHadronIso_dR02_dz02_dxy01 = -999;
  pholead_pho_Cone02ChargedHadronIso_dR02_dz02_dxy01 = -999;
  photrail_pho_Cone02ChargedHadronIso_dR02_dz02_dxy01 = -999;
  pholead_pho_Cone03ChargedHadronIso_dR02_dz02_dxy01 = -999;
  photrail_pho_Cone03ChargedHadronIso_dR02_dz02_dxy01 = -999;
  pholead_pho_Cone04ChargedHadronIso_dR02_dz02_dxy01 = -999;
  photrail_pho_Cone04ChargedHadronIso_dR02_dz02_dxy01 = -999;
  pholead_pho_Cone03PFCombinedIso = -999;
  photrail_pho_Cone03PFCombinedIso = -999;
  pholead_pho_Cone04PFCombinedIso = -999;
  photrail_pho_Cone04PFCombinedIso = -999;
  pholead_PhoPassConvSafeElectronVeto = -999;
  photrail_PhoPassConvSafeElectronVeto = -999;
  pholead_GenPhotonIsoDR04 = -999;
  photrail_GenPhotonIsoDR04 = -999;
  pholead_PhoIso03Ecal = -999;
  pholead_PhoIso03Hcal = -999;
  pholead_PhoIso03TrkSolid = -999;
  pholead_PhoIso03TrkHollow = -999;
  pholead_PhoIso03 = -999;
  pholead_PhoIso04Ecal = -999;
  pholead_PhoIso04Hcal = -999;
  pholead_PhoIso04TrkSolid = -999;
  pholead_PhoIso04TrkHollow = -999;
  pholead_PhoIso04 = -999;
  photrail_PhoIso03Ecal = -999;
  photrail_PhoIso03Hcal = -999;
  photrail_PhoIso03TrkSolid = -999;
  photrail_PhoIso03TrkHollow = -999;
  photrail_PhoIso03 = -999;
  photrail_PhoIso04Ecal = -999;
  photrail_PhoIso04Hcal = -999;
  photrail_PhoIso04TrkSolid = -999;
  photrail_PhoIso04TrkHollow = -999;
  photrail_PhoIso04 = -999;
  pholead_PhoS4OverS1 = -999;
  pholead_PhoSigmaEtaEta = -999;
  pholead_PhoE1x5 = -999;
  pholead_PhoE2x5 = -999;
  pholead_PhoE3x3 = -999;
  pholead_PhoE5x5 = -999;
  pholead_PhomaxEnergyXtal = -999;
  pholead_PhoIso03HcalDepth1 = -999;
  pholead_PhoIso03HcalDepth2 = -999;
  pholead_PhoIso04HcalDepth1 = -999;
  pholead_PhoIso04HcalDepth2 = -999;
  pholead_PhoIso03nTrksSolid = -999;
  pholead_PhoIso03nTrksHollow = -999;
  pholead_PhoIso04nTrksSolid = -999;
  pholead_PhoIso04nTrksHollow = -999;
  pholead_Pho_ChargedHadronIso = -999;
  pholead_Pho_NeutralHadronIso = -999;
  pholead_Pho_PhotonIso = -999;
  pholead_Pho_isPFPhoton = -999;
  pholead_Pho_isPFElectron = -999;
  photrail_PhoS4OverS1 = -999;
  photrail_PhoSigmaEtaEta = -999;
  photrail_PhoE1x5 = -999;
  photrail_PhoE2x5 = -999;
  photrail_PhoE3x3 = -999;
  photrail_PhoE5x5 = -999;
  photrail_PhomaxEnergyXtal = -999;
  photrail_PhoIso03HcalDepth1 = -999;
  photrail_PhoIso03HcalDepth2 = -999;
  photrail_PhoIso04HcalDepth1 = -999;
  photrail_PhoIso04HcalDepth2 = -999;
  photrail_PhoIso03nTrksSolid = -999;
  photrail_PhoIso03nTrksHollow = -999;
  photrail_PhoIso04nTrksSolid = -999;
  photrail_PhoIso04nTrksHollow = -999;
  photrail_Pho_ChargedHadronIso = -999;
  photrail_Pho_NeutralHadronIso = -999;
  photrail_Pho_PhotonIso = -999;
  photrail_Pho_isPFPhoton = -999;
  photrail_Pho_isPFElectron = -999;
  pholead_PhoMCmatchindex = -999;
  pholead_PhoMCmatchexitcode = -999;
  photrail_PhoMCmatchindex = -999;
  photrail_PhoMCmatchexitcode = -999;
};
