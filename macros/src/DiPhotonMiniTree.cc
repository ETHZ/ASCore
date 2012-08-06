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
  randomgen = new TRandom3(0);

  eff_area_EB = 0.406;
  eff_area_EE = 0.528;
}

DiPhotonMiniTree::~DiPhotonMiniTree(){
  delete phocorr;
  delete randomgen;
}

void DiPhotonMiniTree::Begin(){

  cout << "Begin" << endl;

  fOutputFile->cd();

  OutputTree[0] = new TTree("Tree_standard_sel","Tree_standard_sel");
  OutputTree[1] = new TTree("Tree_signal_template","Tree_signal_template");
  OutputTree[2] = new TTree("Tree_background_template","Tree_background_template");
  OutputTree[3] = new TTree("Tree_DY_sel","Tree_DY_sel");
  OutputTree[4] = new TTree("Tree_randomcone_signal_template","Tree_randomcone_signal_template");
  OutputTree[5] = new TTree("Tree_impinging_track_template", "Tree_impinging_track_template");
  OutputTree[6] = new TTree("Tree_singlephoton_inclusive_sel","Tree_singlephoton_inclusive_sel");

  //  histo_PFPhotonDepositAroundImpingingTrack = new TH1F("PFPhotonDepositAroundImpingingTrack","PFPhotonDepositAroundImpingingTrack",50,0,0.2);

  for (int i=0; i<7; i++){

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

  OutputTree[i]->Branch("pholead_hasimpingingtrack",&pholead_hasimpingingtrack,"pholead_hasimpingingtrack/I");
  OutputTree[i]->Branch("photrail_hasimpingingtrack",&photrail_hasimpingingtrack,"photrail_hasimpingingtrack/I");

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
  if (!isdata) {
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
    event_Kfactor = kfactors[2-nmatched_part_isrfsr_gamma];
  }
  else event_Kfactor=1;

  //  cout << "C" << endl;

  // 0 = standard selection for data and MC
  // 1 = signal template generation from MC
  // 2 = background template generation from MC
  // 3 = DY selection (standard with reversed pixel veto)

  if (!TriggerSelection()) return;

  std::vector<int> passing_selection[7];

  bool pass[7];

  for (int sel_cat=0; sel_cat<7; sel_cat++){

    if (isdata){ // do not run these cats on data
      if (sel_cat==1) continue;
      if (sel_cat==2) continue;
    }

    pass[sel_cat]=false;

    std::vector<int> passing;

    for (int i=0; i<fTR->NPhotons; i++){
      passing.push_back(i);
    }

    passing = PhotonPreSelection(fTR,passing);

    // comment this block for the noselection running
    {
      if (sel_cat!=3) passing = ApplyPixelVeto(fTR,passing,0);
      if (sel_cat==3) passing = ApplyPixelVeto(fTR,passing,1);
      passing = PhotonSelection(fTR,passing);
      if (sel_cat==5) passing = ImpingingTrackSelection(fTR,passing,false); // select impinging tracks
      else if (sel_cat==6) passing=passing;
      else passing = ImpingingTrackSelection(fTR,passing,true); // (inverted selection)
    }

    if (sel_cat==0 || sel_cat==3){
      pass[sel_cat] = StandardEventSelection(fTR,passing);
    }
    if (sel_cat==1 || sel_cat==2){
      if (sel_cat==1) passing = SignalSelection(fTR,passing);
      if (sel_cat==2) passing = BackgroundSelection(fTR,passing);
      pass[sel_cat] = SinglePhotonEventSelection(fTR,passing);
    }
    if (sel_cat==4){ // random cone
      //      passing = SignalSelection(fTR,passing); // uncomment to make random cone only from true photons (only in MC!)
      pass[sel_cat] = SinglePhotonEventSelection(fTR,passing);
    }
    if (sel_cat==5){ // impinging track
      pass[sel_cat] = SinglePhotonEventSelection(fTR,passing);
    }
    if (sel_cat==6){
      pass[sel_cat] = SinglePhotonEventSelection(fTR,passing);
    }

    passing_selection[sel_cat] = passing;

  }

  //  cout << "D" << endl;

  for (int sel_cat=0; sel_cat<7; sel_cat++){

    if (!pass[sel_cat]) continue;

    std::vector<int> passing = passing_selection[sel_cat];

    int minsize=999;
    if (sel_cat==0 || sel_cat==3) minsize=2;
    if (sel_cat==1 || sel_cat==2 || sel_cat==4 || sel_cat==5 || sel_cat==6) minsize=1;

    if (!(passing_selection[sel_cat].size()>=minsize)){
      std::cout << "Error!!!" << std::endl;
      continue;
    }

    if (sel_cat==0 || sel_cat==3){
      ResetVars();
      FillLead(passing.at(0));
      FillTrail(passing.at(1));
      float invmass0 = (CorrPhoton(fTR,passing.at(0),0)+CorrPhoton(fTR,passing.at(1),0)).M();
      float invmass5 = (CorrPhoton(fTR,passing.at(0),5)+CorrPhoton(fTR,passing.at(1),5)).M();
      float invmass6 = (CorrPhoton(fTR,passing.at(0),6)+CorrPhoton(fTR,passing.at(1),6)).M();
      dipho_mgg_photon = invmass0;
      dipho_mgg_newCorr = invmass5;
      dipho_mgg_newCorrLocal = invmass6;
      OutputTree[sel_cat]->Fill();
    }

    if (sel_cat==1 || sel_cat==2 || sel_cat==4 || sel_cat==5 || sel_cat==6){
      for (int i=0; i<passing.size(); i++){
      ResetVars();
      FillLead(passing.at(i));
      if (sel_cat==4) pholead_pho_Cone04PhotonIso_dEta015EB_dR070EE_mvVtx = RandomConePhotonIsolation(fTR,passing.at(i));
      OutputTree[sel_cat]->Fill();
      }
    }

  }
 
  //  cout << "E" << endl;

}



void DiPhotonMiniTree::FillLead(int index){

  pholead_eta = fTR->PhoEta[index];
  pholead_px = fTR->PhoPx[index];
  pholead_py = fTR->PhoPy[index];
  pholead_pt = fTR->PhoPt[index];
  pholead_pz = fTR->PhoPz[index];
  pholead_energy = fTR->PhoEnergy[index];
  pholead_SCeta = fTR->SCEta[fTR->PhotSCindex[index]];
  pholead_SCphi = fTR->SCPhi[fTR->PhotSCindex[index]];
  pholead_PhoHasPixSeed=fTR->PhoHasPixSeed[index];
  pholead_PhoHasConvTrks=fTR->PhoHasConvTrks[index];
  pholead_PhoScSeedSeverity=fTR->PhoScSeedSeverity[index];
  pholead_energySCdefault = CorrPhoton(fTR,index,0).E();
  pholead_energyNewCorr = CorrPhoton(fTR,index,5).E();
  pholead_energyNewCorrLocal = CorrPhoton(fTR,index,6).E();
  pholead_r9 = fTR->SCR9[fTR->PhotSCindex[index]];
  pholead_sieie = fTR->PhoSigmaIetaIeta[index];
  pholead_hoe = fTR->PhoHoverE[index];
  pholead_brem = fTR->SCBrem[fTR->PhotSCindex[index]];
  pholead_sigmaPhi = fTR->SCPhiWidth[fTR->PhotSCindex[index]];
  pholead_sigmaEta = fTR->SCEtaWidth[fTR->PhotSCindex[index]];
  pholead_pho_Cone01PhotonIso_dEta015EB_dR070EE_mvVtx=fTR->pho_Cone01PhotonIso_dEta015EB_dR070EE_mvVtx[index];
  pholead_pho_Cone02PhotonIso_dEta015EB_dR070EE_mvVtx=fTR->pho_Cone02PhotonIso_dEta015EB_dR070EE_mvVtx[index];
  pholead_pho_Cone03PhotonIso_dEta015EB_dR070EE_mvVtx=fTR->pho_Cone03PhotonIso_dEta015EB_dR070EE_mvVtx[index];
  pholead_pho_Cone04PhotonIso_dEta015EB_dR070EE_mvVtx=fTR->pho_Cone04PhotonIso_dEta015EB_dR070EE_mvVtx[index];
  pholead_pho_Cone01NeutralHadronIso_mvVtx=fTR->pho_Cone01NeutralHadronIso_mvVtx[index];
  pholead_pho_Cone02NeutralHadronIso_mvVtx=fTR->pho_Cone02NeutralHadronIso_mvVtx[index];
  pholead_pho_Cone03NeutralHadronIso_mvVtx=fTR->pho_Cone03NeutralHadronIso_mvVtx[index];
  pholead_pho_Cone04NeutralHadronIso_mvVtx=fTR->pho_Cone04NeutralHadronIso_mvVtx[index];
  pholead_pho_Cone01ChargedHadronIso_dR02_dz02_dxy01=fTR->pho_Cone01ChargedHadronIso_dR02_dz02_dxy01[index];
  pholead_pho_Cone02ChargedHadronIso_dR02_dz02_dxy01=fTR->pho_Cone02ChargedHadronIso_dR02_dz02_dxy01[index];
  pholead_pho_Cone03ChargedHadronIso_dR02_dz02_dxy01=fTR->pho_Cone03ChargedHadronIso_dR02_dz02_dxy01[index];
  pholead_pho_Cone04ChargedHadronIso_dR02_dz02_dxy01=fTR->pho_Cone04ChargedHadronIso_dR02_dz02_dxy01[index];
  pholead_pho_Cone03PFCombinedIso=fTR->pho_Cone03PFCombinedIso[index];
  pholead_pho_Cone04PFCombinedIso=fTR->pho_Cone04PFCombinedIso[index];
  pholead_PhoPassConvSafeElectronVeto=fTR->PhoPassConvSafeElectronVeto[index];
  if (fTR->PhoMCmatchindex[index]>=0) pholead_GenPhotonIsoDR04=fTR->GenPhotonIsoDR04[fTR->PhoMCmatchindex[index]];
  pholead_PhoIso03Ecal=fTR->PhoIso03Ecal[index];
  pholead_PhoIso03Hcal=fTR->PhoIso03Hcal[index];
  pholead_PhoIso03TrkSolid=fTR->PhoIso03TrkSolid[index];
  pholead_PhoIso03TrkHollow=fTR->PhoIso03TrkHollow[index];
  pholead_PhoIso03=fTR->PhoIso03[index];
  pholead_PhoIso04Ecal=fTR->PhoIso04Ecal[index];
  pholead_PhoIso04Hcal=fTR->PhoIso04Hcal[index];
  pholead_PhoIso04TrkSolid=fTR->PhoIso04TrkSolid[index];
  pholead_PhoIso04TrkHollow=fTR->PhoIso04TrkHollow[index];
  pholead_PhoIso04=fTR->PhoIso04[index];
  pholead_PhoE1OverE9=fTR->PhoE1OverE9[index];
  pholead_PhoS4OverS1=fTR->PhoS4OverS1[index];
  pholead_PhoSigmaEtaEta=fTR->PhoSigmaEtaEta[index];
  pholead_PhoE1x5=fTR->PhoE1x5[index];
  pholead_PhoE2x5=fTR->PhoE2x5[index];
  pholead_PhoE3x3=fTR->PhoE3x3[index];
  pholead_PhoE5x5=fTR->PhoE5x5[index];
  pholead_PhomaxEnergyXtal=fTR->PhomaxEnergyXtal[index];
  pholead_PhoIso03HcalDepth1=fTR->PhoIso03HcalDepth1[index];
  pholead_PhoIso03HcalDepth2=fTR->PhoIso03HcalDepth2[index];
  pholead_PhoIso04HcalDepth1=fTR->PhoIso04HcalDepth1[index];
  pholead_PhoIso04HcalDepth2=fTR->PhoIso04HcalDepth2[index];
  pholead_PhoIso03nTrksSolid=fTR->PhoIso03nTrksSolid[index];
  pholead_PhoIso03nTrksHollow=fTR->PhoIso03nTrksHollow[index];
  pholead_PhoIso04nTrksSolid=fTR->PhoIso04nTrksSolid[index];
  pholead_PhoIso04nTrksHollow=fTR->PhoIso04nTrksHollow[index];
  pholead_Pho_ChargedHadronIso=fTR->Pho_ChargedHadronIso[index];
  pholead_Pho_NeutralHadronIso=fTR->Pho_NeutralHadronIso[index];
  pholead_Pho_PhotonIso=fTR->Pho_PhotonIso[index];
  pholead_Pho_isPFPhoton=fTR->Pho_isPFPhoton[index];
  pholead_Pho_isPFElectron=fTR->Pho_isPFElectron[index];
  pholead_PhoMCmatchindex=fTR->PhoMCmatchindex[index];
  pholead_PhoMCmatchexitcode=fTR->PhoMCmatchexitcode[index];
  TVector3 photon_position = TVector3(fTR->SCx[fTR->PhotSCindex[index]],fTR->SCy[fTR->PhotSCindex[index]],fTR->SCz[fTR->PhotSCindex[index]]);
  TVector3 phovtx(fTR->PhoVx[index],fTR->PhoVy[index],fTR->PhoVz[index]);
  int a;
  pholead_hasimpingingtrack = FindImpingingTrack(fTR,photon_position,phovtx,a,false,GetPFCandRemovals(fTR,index));

};

void DiPhotonMiniTree::FillTrail(int index){

  photrail_eta = fTR->PhoEta[index];
  photrail_px = fTR->PhoPx[index];
  photrail_py = fTR->PhoPy[index];
  photrail_pt = fTR->PhoPt[index];
  photrail_pz = fTR->PhoPz[index];
  photrail_energy = fTR->PhoEnergy[index];
  photrail_SCeta = fTR->SCEta[fTR->PhotSCindex[index]];
  photrail_SCphi = fTR->SCPhi[fTR->PhotSCindex[index]];
  photrail_PhoHasPixSeed=fTR->PhoHasPixSeed[index];
  photrail_PhoHasConvTrks=fTR->PhoHasConvTrks[index];
  photrail_PhoScSeedSeverity=fTR->PhoScSeedSeverity[index];
  photrail_energySCdefault = CorrPhoton(fTR,index,0).E();
  photrail_energyNewCorr = CorrPhoton(fTR,index,5).E();
  photrail_energyNewCorrLocal = CorrPhoton(fTR,index,6).E();
  photrail_r9 = fTR->SCR9[fTR->PhotSCindex[index]];
  photrail_sieie = fTR->PhoSigmaIetaIeta[index];
  photrail_hoe = fTR->PhoHoverE[index];
  photrail_brem = fTR->SCBrem[fTR->PhotSCindex[index]];
  photrail_sigmaPhi = fTR->SCPhiWidth[fTR->PhotSCindex[index]];
  photrail_sigmaEta = fTR->SCEtaWidth[fTR->PhotSCindex[index]];
  photrail_pho_Cone01PhotonIso_dEta015EB_dR070EE_mvVtx=fTR->pho_Cone01PhotonIso_dEta015EB_dR070EE_mvVtx[index];
  photrail_pho_Cone02PhotonIso_dEta015EB_dR070EE_mvVtx=fTR->pho_Cone02PhotonIso_dEta015EB_dR070EE_mvVtx[index];
  photrail_pho_Cone03PhotonIso_dEta015EB_dR070EE_mvVtx=fTR->pho_Cone03PhotonIso_dEta015EB_dR070EE_mvVtx[index];
  photrail_pho_Cone04PhotonIso_dEta015EB_dR070EE_mvVtx=fTR->pho_Cone04PhotonIso_dEta015EB_dR070EE_mvVtx[index];
  photrail_pho_Cone01NeutralHadronIso_mvVtx=fTR->pho_Cone01NeutralHadronIso_mvVtx[index];
  photrail_pho_Cone02NeutralHadronIso_mvVtx=fTR->pho_Cone02NeutralHadronIso_mvVtx[index];
  photrail_pho_Cone03NeutralHadronIso_mvVtx=fTR->pho_Cone03NeutralHadronIso_mvVtx[index];
  photrail_pho_Cone04NeutralHadronIso_mvVtx=fTR->pho_Cone04NeutralHadronIso_mvVtx[index];
  photrail_pho_Cone01ChargedHadronIso_dR02_dz02_dxy01=fTR->pho_Cone01ChargedHadronIso_dR02_dz02_dxy01[index];
  photrail_pho_Cone02ChargedHadronIso_dR02_dz02_dxy01=fTR->pho_Cone02ChargedHadronIso_dR02_dz02_dxy01[index];
  photrail_pho_Cone03ChargedHadronIso_dR02_dz02_dxy01=fTR->pho_Cone03ChargedHadronIso_dR02_dz02_dxy01[index];
  photrail_pho_Cone04ChargedHadronIso_dR02_dz02_dxy01=fTR->pho_Cone04ChargedHadronIso_dR02_dz02_dxy01[index];
  photrail_pho_Cone03PFCombinedIso=fTR->pho_Cone03PFCombinedIso[index];
  photrail_pho_Cone04PFCombinedIso=fTR->pho_Cone04PFCombinedIso[index];
  photrail_PhoPassConvSafeElectronVeto=fTR->PhoPassConvSafeElectronVeto[index];
  if (fTR->PhoMCmatchindex[index]>=0) photrail_GenPhotonIsoDR04=fTR->GenPhotonIsoDR04[fTR->PhoMCmatchindex[index]];
  photrail_PhoIso03Ecal=fTR->PhoIso03Ecal[index];
  photrail_PhoIso03Hcal=fTR->PhoIso03Hcal[index];
  photrail_PhoIso03TrkSolid=fTR->PhoIso03TrkSolid[index];
  photrail_PhoIso03TrkHollow=fTR->PhoIso03TrkHollow[index];
  photrail_PhoIso03=fTR->PhoIso03[index];
  photrail_PhoIso04Ecal=fTR->PhoIso04Ecal[index];
  photrail_PhoIso04Hcal=fTR->PhoIso04Hcal[index];
  photrail_PhoIso04TrkSolid=fTR->PhoIso04TrkSolid[index];
  photrail_PhoIso04TrkHollow=fTR->PhoIso04TrkHollow[index];
  photrail_PhoIso04=fTR->PhoIso04[index];
  photrail_PhoE1OverE9=fTR->PhoE1OverE9[index];
  photrail_PhoS4OverS1=fTR->PhoS4OverS1[index];
  photrail_PhoSigmaEtaEta=fTR->PhoSigmaEtaEta[index];
  photrail_PhoE1x5=fTR->PhoE1x5[index];
  photrail_PhoE2x5=fTR->PhoE2x5[index];
  photrail_PhoE3x3=fTR->PhoE3x3[index];
  photrail_PhoE5x5=fTR->PhoE5x5[index];
  photrail_PhomaxEnergyXtal=fTR->PhomaxEnergyXtal[index];
  photrail_PhoIso03HcalDepth1=fTR->PhoIso03HcalDepth1[index];
  photrail_PhoIso03HcalDepth2=fTR->PhoIso03HcalDepth2[index];
  photrail_PhoIso04HcalDepth1=fTR->PhoIso04HcalDepth1[index];
  photrail_PhoIso04HcalDepth2=fTR->PhoIso04HcalDepth2[index];
  photrail_PhoIso03nTrksSolid=fTR->PhoIso03nTrksSolid[index];
  photrail_PhoIso03nTrksHollow=fTR->PhoIso03nTrksHollow[index];
  photrail_PhoIso04nTrksSolid=fTR->PhoIso04nTrksSolid[index];
  photrail_PhoIso04nTrksHollow=fTR->PhoIso04nTrksHollow[index];
  photrail_Pho_ChargedHadronIso=fTR->Pho_ChargedHadronIso[index];
  photrail_Pho_NeutralHadronIso=fTR->Pho_NeutralHadronIso[index];
  photrail_Pho_PhotonIso=fTR->Pho_PhotonIso[index];
  photrail_Pho_isPFPhoton=fTR->Pho_isPFPhoton[index];
  photrail_Pho_isPFElectron=fTR->Pho_isPFElectron[index];
  photrail_PhoMCmatchindex=fTR->PhoMCmatchindex[index];
  photrail_PhoMCmatchexitcode=fTR->PhoMCmatchexitcode[index];
  TVector3 photon_position = TVector3(fTR->SCx[fTR->PhotSCindex[index]],fTR->SCy[fTR->PhotSCindex[index]],fTR->SCz[fTR->PhotSCindex[index]]);
  TVector3 phovtx(fTR->PhoVx[index],fTR->PhoVy[index],fTR->PhoVz[index]);
  int a;
  pholead_hasimpingingtrack = FindImpingingTrack(fTR,photon_position,phovtx,a,false,GetPFCandRemovals(fTR,index));

};


void DiPhotonMiniTree::End(){
  fOutputFile->cd();
  for (int i=0; i<7; i++) OutputTree[i]->Write();	
  fHNumPU->Write();
  fHNumVtx->Write();
  //  histo_PFPhotonDepositAroundImpingingTrack->Write();
	
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
    if (fTR->PhotSCindex[*it]<0) it=passing.erase(it); else it++;
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
    float eta=fTR->SCEta[fTR->PhotSCindex[*it]];
    const float eff_area = (fabs(eta)<1.4442) ? eff_area_EB : eff_area_EE;
    const float dR=0.4;
    float puenergy =3.14*dR*dR*eff_area*fTR->Rho;
    if (fTR->pho_Cone04PFCombinedIso[*it]*fTR->PhoPt[*it]-puenergy<5) pass=1;
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


bool DiPhotonMiniTree::SinglePhotonEventSelection(TreeReader *fTR, std::vector<int> &passing){

  for (vector<int>::iterator it = passing.begin(); it != passing.end(); ){
    bool pass=0;
    if (fTR->PhoPt[*it]>40) pass=1;
    if (!pass) it=passing.erase(it); else it++;
  }

  if (passing.size()<1) return false;

  return true;

};

bool DiPhotonMiniTree::StandardEventSelection(TreeReader *fTR, std::vector<int> &passing){

  if (passing.size()<2) return false;

  passing.resize(2); // keep only the first two


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

//double DiPhotonMiniTree::etaTransformation(  float EtaParticle , float Zvertex)  {
//
//  //---Definitions
//  const float pi = 3.1415927;
//
//  //---Definitions for ECAL
//  const float R_ECAL           = 136.5;
//  const float Z_Endcap         = 328.0;
//  const float etaBarrelEndcap  = 1.479; 
//   
//  //---ETA correction
//
//  float Theta = 0.0  ; 
//  float ZEcal = R_ECAL*sinh(EtaParticle)+Zvertex;
//
//  if(ZEcal != 0.0) Theta = atan(R_ECAL/ZEcal);
//  if(Theta<0.0) Theta = Theta+pi ;
//  double ETA = - log(tan(0.5*Theta));
//         
//  if( fabs(ETA) > etaBarrelEndcap )
//    {
//      float Zend = Z_Endcap ;
//      if(EtaParticle<0.0 )  Zend = -Zend ;
//      float Zlen = Zend - Zvertex ;
//      float RR = Zlen/sinh(EtaParticle); 
//      Theta = atan(RR/Zend);
//      if(Theta<0.0) Theta = Theta+pi ;
//      ETA = - log(tan(0.5*Theta));		      
//    } 
//  //---Return the result
//  return ETA;
//  //---end
//}

double DiPhotonMiniTree::phiNorm(float phi) {

  const float pi = 3.1415927;
  const float twopi = 2.0*pi;

  if(phi >  pi) {phi = phi - twopi;}
  if(phi < -pi) {phi = phi + twopi;}

  return phi;
}

bool DiPhotonMiniTree::FindCloseJetsAndPhotons(TreeReader *fTR, float eta, float phi, int phoqi){

  const bool debug=false;
  if (debug) std::cout << "calling FindCloseJetsAndPhotons eta=" << eta << " phi=" << phi << std::endl;

  const float mindR = 0.4;
  bool found=false;

  for (int i=0; i<fTR->NJets; i++){
    if (fTR->JPt[i]<20) continue;
    float dR = Util::GetDeltaR(eta,fTR->JEta[i],phi,fTR->JPhi[i]);
    if (dR<mindR) found=true;
    if (debug) if (dR<mindR) std::cout << "Found jet eta=" << fTR->JEta[i] << " phi=" << fTR->JPhi[i] << std::endl;
  }

  for (int i=0; i<fTR->NPhotons; i++){
    if (fTR->PhoPt[i]<10) continue;
    float dR = Util::GetDeltaR(eta,fTR->PhoEta[i],phi,fTR->PhoPhi[i]);
    if (dR<mindR) found=true;
    if (debug) if (dR<mindR) std::cout << "Found phot eta=" << fTR->PhoEta[i] << " phi=" << fTR->PhoPhi[i] << std::endl;
  }


  const float eff_area = (fabs(eta)<1.4442) ? eff_area_EB : eff_area_EE;
  float puenergy =3.14*0.4*0.4*eff_area*fTR->Rho;
  if (CombinedPFIsolation(eta,phi,phoqi)-puenergy>5) found=true;

  TVector3 photon_position = TVector3(fTR->SCx[fTR->PhotSCindex[phoqi]],fTR->SCy[fTR->PhotSCindex[phoqi]],fTR->SCz[fTR->PhotSCindex[phoqi]]);
  photon_position.SetPtEtaPhi(photon_position.Perp(),eta,phi);
  TVector3 phovtx(fTR->PhoVx[phoqi],fTR->PhoVy[phoqi],fTR->PhoVz[phoqi]);
    
  int a;
  bool found_impinging = FindImpingingTrack(fTR,photon_position,phovtx,a,false,GetPFCandRemovals(fTR,phoqi));
  if (found_impinging) found=true;

  if (debug) std::cout << "returning " << found << std::endl;
  return found;

};

std::vector<int> DiPhotonMiniTree::GetPFCandRemovals(TreeReader *fTR, int phoqi){
  std::vector<int> out;
  if (fTR->Pho_isPFPhoton[phoqi]) out.push_back(fTR->pho_matchedPFPhotonCand[phoqi]);
  if (fTR->Pho_isPFElectron[phoqi]) out.push_back(fTR->pho_matchedPFElectronCand[phoqi]);
  return out;
}

float DiPhotonMiniTree::RandomConePhotonIsolation(TreeReader *fTR, int phoqi){

  float result=0;
  
  //  double sceta = fTR->SCEta[fTR->PhotSCindex[phoqi]];
  //  double scphi = fTR->SCPhi[fTR->PhotSCindex[phoqi]];
  
  TVector3 photon_true_position = TVector3(fTR->SCx[fTR->PhotSCindex[phoqi]],fTR->SCy[fTR->PhotSCindex[phoqi]],fTR->SCz[fTR->PhotSCindex[phoqi]]);
  
  double true_sceta = photon_true_position.Eta();
  double true_scphi = photon_true_position.Phi();

  double rotated_scphi = phiNorm(true_scphi+0.5*TMath::Pi());

  bool isok = !(FindCloseJetsAndPhotons(fTR,true_sceta,rotated_scphi,phoqi));
  if (!isok) {
    rotated_scphi=phiNorm(true_scphi-0.5*TMath::Pi());
    isok=!(FindCloseJetsAndPhotons(fTR,true_sceta,rotated_scphi,phoqi));
  }
  
    int count=0;
    while (!isok && count<10) {
      rotated_scphi = phiNorm(true_scphi+randomgen->Uniform(0.8,6.28-0.8));
      isok=!(FindCloseJetsAndPhotons(fTR,true_sceta,rotated_scphi,phoqi));
      count++;
    }

    if (count==10){
      std::cout << "Error in random cone generation!!!" << std::endl;
      return -999;
    };

    TVector3 photon_rotated_position(photon_true_position.x(),photon_true_position.y(),photon_true_position.z());
    photon_rotated_position.SetPhi(rotated_scphi);

    for (int i=0; i<fTR->NPfCand; i++){
      
      if (fTR->Pho_isPFPhoton[phoqi] && fTR->pho_matchedPFPhotonCand[phoqi]==i) continue;
      if (fTR->Pho_isPFElectron[phoqi] && fTR->pho_matchedPFElectronCand[phoqi]==i) continue;

    if (fTR->PfCandPdgId[i]!=22) continue;

    TVector3 pfvertex(fTR->PfCandVx[i],fTR->PfCandVy[i],fTR->PfCandVz[i]);
    TVector3 rotated_photon_direction = photon_rotated_position-pfvertex;

    float sceta = rotated_photon_direction.Eta();
    float scphi = rotated_photon_direction.Phi();

    double dEta = fTR->PfCandEta[i] - sceta;
    double dPhi = Util::DeltaPhi(fTR->PfCandPhi[i],scphi);

    double pt = fTR->PfCandPt[i];
    double dR = sqrt(dEta*dEta+dPhi*dPhi);

    if (dR>0.4) continue;

    if (fTR->PhoisEB[phoqi]){
      if (fabs(dEta)<0.015) continue;
    }
    else if (fTR->PhoisEE[phoqi]){
      float limit_dR = 0.00864*fabs(sinh(sceta))*4;
      if (dR<limit_dR) continue;
    }
    else {
      std::cout << "Something wrong " << fTR->PhoEta[phoqi] << " " << photon_true_position.Eta() << " " << fTR->PhoisEB[phoqi] << fTR->PhoisEE[phoqi] << std::endl;
      return -999;
    }

    result+=pt;

  } // end pf cand loop

  return result;

};

std::vector<int> DiPhotonMiniTree::ImpingingTrackSelection(TreeReader *fTR, std::vector<int> passing, bool invert){

  for (int i=0; i<100; i++) impinging_track_pfcand[i]=-999;

  for (vector<int>::iterator it = passing.begin(); it != passing.end(); ){
    
    bool found=0;
    int phoqi=*it;
    
    TVector3 photon_position = TVector3(fTR->SCx[fTR->PhotSCindex[phoqi]],fTR->SCy[fTR->PhotSCindex[phoqi]],fTR->SCz[fTR->PhotSCindex[phoqi]]);
    TVector3 phovtx(fTR->PhoVx[phoqi],fTR->PhoVy[phoqi],fTR->PhoVz[phoqi]);
    
    found = FindImpingingTrack(fTR,photon_position,phovtx,impinging_track_pfcand[phoqi],invert,GetPFCandRemovals(fTR,phoqi));
    
    if (!found) it=passing.erase(it); else it++;
    
  }

  return passing;
  
};


bool DiPhotonMiniTree::FindImpingingTrack(TreeReader *fTR, TVector3 caloposition, TVector3 refvertex, int &reference_index_found, bool invert, std::vector<int> removals){

  bool found = false;

  TVector3 photon_position = caloposition;
  TVector3 phovtx = refvertex;

  for (int i=0; i<fTR->NPfCand; i++){
    
    bool removed = false;
    for (int j=0; j<removals.size(); j++) {
      if (i==removals.at(j)) removed=true;
    }
    if (removed) continue;
    
    float id = fTR->PfCandPdgId[i];
    if (!(fabs(id)==211 || fabs(id)==321 || id==999211 || fabs(id)==2212)) continue;
    
    TVector3 pfvertex(fTR->PfCandVx[i],fTR->PfCandVy[i],fTR->PfCandVz[i]);
    TVector3 photon_direction = photon_position-pfvertex;
      
    double sceta = photon_direction.Eta();
    double scphi = photon_direction.Phi();

    double dEta = fTR->PfCandEta[i] - sceta;
    double dPhi = Util::DeltaPhi(fTR->PfCandPhi[i],scphi);
    double dR = sqrt(dEta*dEta+dPhi*dPhi);

    double pt = fTR->PfCandPt[i];

    TVector3 vtxmom(fTR->PfCandTrackRefPx[i],fTR->PfCandTrackRefPy[i],fTR->PfCandTrackRefPz[i]);

    if (vtxmom.x()==-999 || vtxmom.y()==-999 || vtxmom.z()==-999) {
      std::cout << "Something wrong with vtxmom from trackref, fallback" << std::endl;
      vtxmom = TVector3(fTR->PfCandPx[i],fTR->PfCandPy[i],fTR->PfCandPz[i]);
    }
      
    double dxy = ( -(pfvertex.x()-phovtx.x())*vtxmom.y() +(pfvertex.y()-phovtx.y())*vtxmom.x() ) / vtxmom.Perp();
    double dz = (pfvertex.z()-phovtx.z()) - ( (pfvertex.x()-phovtx.x())*vtxmom.x() + (pfvertex.y()-phovtx.y())*vtxmom.y() ) / vtxmom.Perp() * vtxmom.z() / vtxmom.Perp();
    dxy=fabs(dxy);
    dz=fabs(dz);

    if (pt<3) continue;
    if (dR>0.4) continue;

    if (dz>0.2) continue;
    if (dxy>0.1) continue;
    if (dR<0.02) continue;

    if (!(fTR->PfCandHasHitInFirstPixelLayer[i])) continue;

    found=1;
    reference_index_found = i;
    break;
      
  } // end pf cand loop

  if (invert) return !found;
  return found;

};


void DiPhotonMiniTree::Fillhist_PFPhotonDepositAroundImpingingTrack(int phoqi, int trkindex){

  return;

//
//    TVector3 photon_position = TVector3(fTR->SCx[fTR->PhotSCindex[phoqi]],fTR->SCy[fTR->PhotSCindex[phoqi]],fTR->SCz[fTR->PhotSCindex[phoqi]]);
//    TVector3 phovtx(fTR->PhoVx[phoqi],fTR->PhoVy[phoqi],fTR->PhoVz[phoqi]);
//
//    TVector3 trkvertex(fTR->PfCandVx[trkindex],fTR->PfCandVy[trkindex],fTR->PfCandVz[trkindex]);
//    TVector3 photon_direction = photon_position-pfvertex;
//      
//      double sceta = photon_direction.Eta();
//      double scphi = photon_direction.Phi();
//
//    double trketa = pvm_trk.Eta();
//    double trkphi = pvm_trk.Phi();
//
//    double dR_trkpho = Util::GetDeltaR(sceta,trketa,scphi,trkphi);
//    double dEta_trkpho = trketa-sceta;
//    
//    if (dR_trkpho>0.4) return;
//    if (fabs(sceta)<1.479 && fabs(dEta_trkpho)<0.215) return;
//    if (fabs(sceta)>1.479 && dR_trkpho<0.27) return;
//
//    
//    for (int i=0; i<fTR->NPfCand; i++){
//
//      if (fTR->Pho_isPFPhoton[phoqi] && fTR->pho_matchedPFPhotonCand[phoqi]==i) continue;
//      if (fTR->Pho_isPFElectron[phoqi] && fTR->pho_matchedPFElectronCand[phoqi]==i) continue;
//
//      if (fTR->PfCandPdgId[i]!=22) continue;
//
//      TVector3 pfvtx(fTR->PfCandVx[i],fTR->PfCandVy[i],fTR->PfCandVz[i]);
//      TVector3 momentum(fTR->PfCandMomX[i],fTR->PfCandMomY[i],fTR->PfCandMomZ[i]);
//      TVector3 pvm(momentum*r/momentum.R()+pfvtx);
//
//      double dEta_trk = pvm.Eta() - trketa;
//      double dPhi_trk = Util::DeltaPhi(pvm.Phi(),trkphi);
//      double dR_trk = sqrt(dEta_trk*dEta_trk+dPhi_trk*dPhi_trk);
//
//      double dEta_pho = pvm.Eta() - sceta;
//      double dPhi_pho = Util::DeltaPhi(pvm.Phi(),scphi);
//      double dR_pho = sqrt(dEta_pho*dEta_pho+dPhi_pho*dPhi_pho);
//
//      double pt = fTR->PfCandPt[i];
//
//
//      if (dR_trk>0.2) continue;
//
//      bool error=0;
//      if (fabs(sceta)<1.479 && fabs(dEta_pho)<0.015) error=1;
//      if (fabs(sceta)>1.479 && dR_pho<0.07) error=1;
//
//      if (error) {
//	std::cout << "Error in impinging track histo filler!!!" << std::endl;
//	assert(1==0);
//      }
//
//      histo_PFPhotonDepositAroundImpingingTrack->Fill(dR_trk,pt);
//      
//      
//    } // end pf cand loop
//

};

float DiPhotonMiniTree::CombinedPFIsolation(float eta, float phi, int phoqi){

  float result=0;

  TVector3 photon_position = TVector3(fTR->SCx[fTR->PhotSCindex[phoqi]],fTR->SCy[fTR->PhotSCindex[phoqi]],fTR->SCz[fTR->PhotSCindex[phoqi]]);
  TVector3 phovtx(fTR->PhoVx[phoqi],fTR->PhoVy[phoqi],fTR->PhoVz[phoqi]);
 
  for (int i=0; i<fTR->NPfCand; i++){

    if (fTR->Pho_isPFPhoton[phoqi] && fTR->pho_matchedPFPhotonCand[phoqi]==i) continue;
    if (fTR->Pho_isPFElectron[phoqi] && fTR->pho_matchedPFElectronCand[phoqi]==i) continue;

    int type = FindPFCandType(fTR->PfCandPdgId[i]);

    if (!(type==0 || type==1 || type==2)) continue;

    TVector3 pfvertex(fTR->PfCandVx[i],fTR->PfCandVy[i],fTR->PfCandVz[i]);
    TVector3 photon_direction = photon_position-pfvertex;

    float sceta = photon_direction.Eta();
    float scphi = photon_direction.Phi();

    double dEta = fTR->PfCandEta[i] - sceta;
    double dPhi = Util::DeltaPhi(fTR->PfCandPhi[i],scphi);

    double pt = fTR->PfCandPt[i];
    double dR = sqrt(dEta*dEta+dPhi*dPhi);

    double dxy;
    double dz;

    if (type==1){

      TVector3 vtxmom(fTR->PfCandTrackRefPx[i],fTR->PfCandTrackRefPy[i],fTR->PfCandTrackRefPz[i]);

      if (vtxmom.x()==-999 || vtxmom.y()==-999 || vtxmom.z()==-999) {
	std::cout << "Something wrong with vtxmom from trackref, fallback" << std::endl;
	vtxmom = TVector3(fTR->PfCandPx[i],fTR->PfCandPy[i],fTR->PfCandPz[i]);
      }
      
      dxy = ( -(pfvertex.x()-phovtx.x())*vtxmom.y() +(pfvertex.y()-phovtx.y())*vtxmom.x() ) / vtxmom.Perp();
      dz = (pfvertex.z()-phovtx.z()) - ( (pfvertex.x()-phovtx.x())*vtxmom.x() + (pfvertex.y()-phovtx.y())*vtxmom.y() ) / vtxmom.Perp() * vtxmom.z() / vtxmom.Perp();
      dxy=fabs(dxy);
      dz=fabs(dz);

    }

    if (dR>0.4) continue;

    if (type==2){ 

      if (fTR->PhoisEB[phoqi]){
	if (fabs(dEta)<0.015) continue;
      }
      else if (fTR->PhoisEE[phoqi]){
	float limit_dR = 0.00864*fabs(sinh(sceta))*4;
	if (dR<limit_dR) continue;
      }
      else {
	std::cout << "Something wrong" << std::endl;
	return -999;
      }

    }

    if (type==1){

      if (dz>0.2) continue;
      if (dxy>0.1) continue;
      if (dR<0.02) continue;

    }


    result+=pt;

  } // end pf cand loop

  return result;

};


int DiPhotonMiniTree::FindPFCandType(int id){

  int type = -1;

  if (id==111 || id==130 || id==310 || id==2112) type=0; //neutral hadrons
  if (fabs(id)==211 || fabs(id)==321 || id==999211 || fabs(id)==2212) type=1; //charged hadrons
  if (id==22) type=2; //photons
  if (fabs(id)==11) type=3; //electrons
  if (fabs(id)==13) type=4; //muons

  return type;
}
