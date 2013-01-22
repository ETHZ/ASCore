#include "helper/Utilities.hh"
#include "DiPhotonMiniTree.hh"

#include "DiPhotonPurity.hh"


#include <iostream>
#include <fstream>
#include <assert.h>

using namespace std;

DiPhotonMiniTree::DiPhotonMiniTree(TreeReader *tr, std::string dataType, Float_t aw, Float_t* _kfac, Float_t _minthrpfphotoncandEB, Float_t _minthrpfphotoncandEE) : UserAnalysisBase(tr), fDataType_(dataType), AddWeight(aw), kfactors(_kfac), global_minthrpfphotoncandEB(_minthrpfphotoncandEB), global_minthrpfphotoncandEE(_minthrpfphotoncandEE){
  Util::SetStyle();	
  if (fDataType_ == "mc") isdata=false;
  else if (fDataType_ == "data") isdata=true; 
  else {
    std::cout << "wrong data type" << std::endl;
    assert(1==0);
  }
  phocorr = new EnergyCorrection("photons");
  randomgen = new TRandom3(0);

  global_linkbyrechit_enlargement = 0.25; // xtal_size_eff = (1+global_linkbyrechit_enlargement)*xtal_size

  eegeom = TGeoPara(1,1,1,0,0,0);

  // 020616 no cleaning numbers (superseded: should be updated with new ones after removing pf ch iso from presel) XXX
  // EA SET TO ZERO HARDCODED BELOW, these values are ignored;
  float _binsdef_single_gamma_EB_eta[n_templates_EB+1] = {0,0.2,0.4,0.6,0.8,1,1.2,1.4442}; 
  float _binsdef_single_gamma_EE_eta[n_templates_EE+1] = {1.56,1.653,1.8,2,2.2,2.5};
  float _eff_areas_EB_photon_data[n_templates_EB] = {2.601118e-01,2.584915e-01,2.640072e-01,2.656851e-01,2.564615e-01,2.396511e-01,1.645776e-01};
  float _eff_areas_EE_photon_data[n_templates_EE] = {5.783452e-02,8.321881e-02,1.177009e-01,1.422445e-01,1.139434e-01};
  float _eff_areas_EB_charged_data[n_templates_EB] = {2.526058e-02,2.619860e-02,2.118318e-02,2.078713e-02,2.066629e-02,1.908721e-02,1.732135e-02};
  float _eff_areas_EE_charged_data[n_templates_EE] = {1.538706e-02,1.737151e-02,1.410831e-02,1.558452e-02,1.056310e-02};
  float _eff_areas_EB_neutral_data[n_templates_EB] = {9.534780e-02,8.826486e-02,8.199201e-02,8.729101e-02,1.018310e-01,1.244805e-01,1.542784e-01};
  float _eff_areas_EE_neutral_data[n_templates_EE] = {1.547500e-01,1.480421e-01,1.391524e-01,1.546148e-01,2.487922e-01};
  float _eff_areas_EB_photon_MC[n_templates_EB] = {2.728305e-01,2.747126e-01,2.742736e-01,2.689392e-01,2.642000e-01,2.430941e-01,1.701134e-01};
  float _eff_areas_EE_photon_MC[n_templates_EE] = {5.066866e-02,8.511977e-02,1.342678e-01,1.693861e-01,1.448627e-01};
  float _eff_areas_EB_charged_MC[n_templates_EB] = {2.333960e-02,2.682862e-02,2.747902e-02,2.494792e-02,2.206384e-02,1.864648e-02,2.470047e-02};
  float _eff_areas_EE_charged_MC[n_templates_EE] = {2.264766e-02,1.934945e-02,1.834243e-02,2.223225e-02,8.191396e-03};
  float _eff_areas_EB_neutral_MC[n_templates_EB] = {7.527103e-02,7.490141e-02,7.510347e-02,7.806471e-02,9.715760e-02,1.193150e-01,1.449191e-01};
  float _eff_areas_EE_neutral_MC[n_templates_EE] = {1.507419e-01,1.383201e-01,1.489503e-01,1.460166e-01,2.548731e-01};

//  old 020615 w/cleaning numbers
//  float _eff_areas_EB_photon_data[n_templates_EB] = {2.004668e-01,2.000222e-01,2.083325e-01,2.129163e-01,2.082317e-01,1.982015e-01,1.383834e-01}; 
//  float _eff_areas_EE_photon_data[n_templates_EE] = {3.727486e-02,5.494237e-02,7.876623e-02,1.006998e-01,8.432818e-02};
//  float _eff_areas_EB_charged_data[n_templates_EB] = {2.707052e-02,2.715322e-02,2.181513e-02,2.302210e-02,2.296409e-02,2.010706e-02,2.016756e-02};
//  float _eff_areas_EE_charged_data[n_templates_EE] = {1.667700e-02,1.863876e-02,1.442144e-02,1.656359e-02,1.121806e-02}; 
//  float _eff_areas_EB_neutral_data[n_templates_EB] = {1.008816e-01,9.405549e-02,8.752698e-02,9.280508e-02,1.083173e-01,1.327634e-01,1.648043e-01}; 
//  float _eff_areas_EE_neutral_data[n_templates_EE] = {1.665812e-01,1.587098e-01,1.496094e-01,1.669686e-01,2.647152e-01}; 
//  float _eff_areas_EB_photon_MC[n_templates_EB] = {2.049085e-01,2.037967e-01,2.129352e-01,2.154778e-01,2.155464e-01,2.037776e-01,1.407869e-01};
//  float _eff_areas_EE_photon_MC[n_templates_EE] = {4.410657e-02,6.010713e-02,9.360530e-02,1.421708e-01,1.269196e-01};
//  float _eff_areas_EB_charged_MC[n_templates_EB] = {3.131846e-02,3.181332e-02,3.066015e-02,2.876652e-02,2.709936e-02,2.977523e-02,2.038828e-02};
//  float _eff_areas_EE_charged_MC[n_templates_EE] = {2.047554e-02,1.836507e-02,1.890404e-02,1.743312e-02,9.424824e-03}; 
//  float _eff_areas_EB_neutral_MC[n_templates_EB] = {8.039549e-02,8.165860e-02,8.187589e-02,8.591983e-02,9.982728e-02,1.248690e-01,1.519101e-01}; 
//  float _eff_areas_EE_neutral_MC[n_templates_EE] = {1.622648e-01,1.533238e-01,1.456079e-01,1.670986e-01,2.775698e-01}; 


  binsdef_single_gamma_EB_eta.assign(_binsdef_single_gamma_EB_eta,_binsdef_single_gamma_EB_eta+n_templates_EB+1);
  binsdef_single_gamma_EE_eta.assign(_binsdef_single_gamma_EE_eta,_binsdef_single_gamma_EE_eta+n_templates_EE+1);
  eff_areas_EB_photon_data.assign(_eff_areas_EB_photon_data ,_eff_areas_EB_photon_data +n_templates_EB);
  eff_areas_EE_photon_data.assign(_eff_areas_EE_photon_data ,_eff_areas_EE_photon_data +n_templates_EE);
  eff_areas_EB_charged_data.assign(_eff_areas_EB_charged_data ,_eff_areas_EB_charged_data +n_templates_EB);
  eff_areas_EE_charged_data.assign(_eff_areas_EE_charged_data ,_eff_areas_EE_charged_data +n_templates_EE);
  eff_areas_EB_neutral_data.assign(_eff_areas_EB_neutral_data ,_eff_areas_EB_neutral_data +n_templates_EB);
  eff_areas_EE_neutral_data.assign(_eff_areas_EE_neutral_data ,_eff_areas_EE_neutral_data +n_templates_EE);
  eff_areas_EB_photon_MC.assign(_eff_areas_EB_photon_MC ,_eff_areas_EB_photon_MC +n_templates_EB);
  eff_areas_EE_photon_MC.assign(_eff_areas_EE_photon_MC ,_eff_areas_EE_photon_MC +n_templates_EE);
  eff_areas_EB_charged_MC.assign(_eff_areas_EB_charged_MC ,_eff_areas_EB_charged_MC +n_templates_EB);
  eff_areas_EE_charged_MC.assign(_eff_areas_EE_charged_MC ,_eff_areas_EE_charged_MC +n_templates_EE);
  eff_areas_EB_neutral_MC.assign(_eff_areas_EB_neutral_MC ,_eff_areas_EB_neutral_MC +n_templates_EB);
  eff_areas_EE_neutral_MC.assign(_eff_areas_EE_neutral_MC ,_eff_areas_EE_neutral_MC +n_templates_EE);


}

DiPhotonMiniTree::~DiPhotonMiniTree(){
  delete phocorr;
  delete randomgen;
}

void DiPhotonMiniTree::Begin(){

  cout << "Begin" << endl;

  fOutputFile->cd();

  treename[0] = TString("Tree_2Dstandard_selection");
  treename[1] = TString("Tree_1Drandomcone_template");
  treename[2] = TString("Tree_1Dsideband_template");
  treename[3] = TString("Tree_2DZee_pixelvetoreversed_selection");
  treename[4] = TString("Tree_1Dpreselection");
  treename[5] = TString("Tree_1Dselection");
  treename[6] = TString("Tree_2Drandomcone_template");
  treename[7] = TString("Tree_2Drandomconesideband_template");
  treename[8] = TString("Tree_2Dsideband_template");
  treename[9] = TString("Tree_2Dstandard_preselection");
  treename[10] = TString("Tree_2DZmumu_selection");
  treename[11] = TString("Tree_1Dsignal_template");
  treename[12] = TString("Tree_1Dbackground_template");
  treename[13] = TString("Tree_2Dtruesigsig_template");
  treename[14] = TString("Tree_2Dtruesigbkg_template");
  treename[15] = TString("Tree_2Dtruebkgbkg_template");
  treename[16] = TString("Tree_2Drconeplusgenfake_template");
  treename[17] = TString("Tree_2Dgenpromptplussideband_template");

  for (int i=0; i<18; i++){
    OutputTree[i] = new TTree(treename[i].Data(),treename[i].Data());
    is2d[i] = (treename[i].First("2")>-1);
  }

  for (int i=0; i<18; i++){

  OutputTree[i]->Branch("event_luminormfactor",&event_luminormfactor,"event_luminormfactor/F");
  OutputTree[i]->Branch("event_Kfactor",&event_Kfactor,"event_Kfactor/F");

  OutputTree[i]->Branch("event_weight",&event_weight,"event_weight/F");
  OutputTree[i]->Branch("event_weight3D",&event_weight3D,"event_weight3D/F");
  OutputTree[i]->Branch("event_rho",&event_rho,"event_rho/F");
  OutputTree[i]->Branch("event_sigma",&event_sigma,"event_sigma/F");
  OutputTree[i]->Branch("event_nPU",&event_nPU,"event_nPU/I");
  OutputTree[i]->Branch("event_PUOOTnumInteractionsEarly",&event_PUOOTnumInteractionsEarly,"event_PUOOTnumInteractionsEarly/I");
  OutputTree[i]->Branch("event_PUOOTnumInteractionsLate",&event_PUOOTnumInteractionsLate,"event_PUOOTnumInteractionsLate/I");
  OutputTree[i]->Branch("event_nRecVtx",&event_nRecVtx,"event_nRecVtx/I");
  OutputTree[i]->Branch("event_pass12whoissiglike",&event_pass12whoissiglike,"event_pass12whoissiglike/I");

  OutputTree[i]->Branch("event_CSCTightHaloID",&event_CSCTightHaloID,"event_CSCTightHaloID/I");
  OutputTree[i]->Branch("event_NMuons",&event_NMuons,"event_NMuons/I");
  OutputTree[i]->Branch("event_NMuonsTot",&event_NMuonsTot,"event_NMuonsTot/I");

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
  OutputTree[i]->Branch("pholead_Nchargedhadronsincone",&pholead_Nchargedhadronsincone,"pholead_Nchargedhadronsincone/I");
  OutputTree[i]->Branch("photrail_Nchargedhadronsincone",&photrail_Nchargedhadronsincone,"photrail_Nchargedhadronsincone/I");

  OutputTree[i]->Branch("pholead_scarea",&pholead_scarea,"pholead_scarea/F");
  OutputTree[i]->Branch("pholead_scareaSF",&pholead_scareaSF,"pholead_scareaSF/F");
  OutputTree[i]->Branch("photrail_scarea",&photrail_scarea,"photrail_scarea/F");
  OutputTree[i]->Branch("photrail_scareaSF",&photrail_scareaSF,"photrail_scareaSF/F");

  OutputTree[i]->Branch("pholead_Npfcandphotonincone",&pholead_Npfcandphotonincone,"pholead_Npfcandphotonincone/I");
  OutputTree[i]->Branch("pholead_Npfcandchargedincone",&pholead_Npfcandchargedincone,"pholead_Npfcandchargedincone/I");
  OutputTree[i]->Branch("pholead_Npfcandneutralincone",&pholead_Npfcandneutralincone,"pholead_Npfcandneutralincone/I");

  OutputTree[i]->Branch("photrail_Npfcandphotonincone",&photrail_Npfcandphotonincone,"photrail_Npfcandphotonincone/I");
  OutputTree[i]->Branch("photrail_Npfcandchargedincone",&photrail_Npfcandchargedincone,"photrail_Npfcandchargedincone/I");
  OutputTree[i]->Branch("photrail_Npfcandneutralincone",&photrail_Npfcandneutralincone,"photrail_Npfcandneutralincone/I");

  OutputTree[i]->Branch("pholead_photonpfcandenergies",&pholead_photonpfcandenergies,"pholead_photonpfcandenergies[pholead_Npfcandphotonincone]/F");
  OutputTree[i]->Branch("pholead_chargedpfcandenergies",&pholead_chargedpfcandenergies,"pholead_chargedpfcandenergies[pholead_Npfcandchargedincone]/F");
  OutputTree[i]->Branch("pholead_neutralpfcandenergies",&pholead_neutralpfcandenergies,"pholead_neutralpfcandenergies[pholead_Npfcandneutralincone]/F");
  OutputTree[i]->Branch("pholead_photonpfcandets",&pholead_photonpfcandets,"pholead_photonpfcandets[pholead_Npfcandphotonincone]/F");
  OutputTree[i]->Branch("pholead_chargedpfcandets",&pholead_chargedpfcandets,"pholead_chargedpfcandets[pholead_Npfcandchargedincone]/F");
  OutputTree[i]->Branch("pholead_neutralpfcandets",&pholead_neutralpfcandets,"pholead_neutralpfcandets[pholead_Npfcandneutralincone]/F");
  OutputTree[i]->Branch("pholead_photonpfcanddetas",&pholead_photonpfcanddetas,"pholead_photonpfcanddetas[pholead_Npfcandphotonincone]/F");
  OutputTree[i]->Branch("pholead_chargedpfcanddetas",&pholead_chargedpfcanddetas,"pholead_chargedpfcanddetas[pholead_Npfcandchargedincone]/F");
  OutputTree[i]->Branch("pholead_neutralpfcanddetas",&pholead_neutralpfcanddetas,"pholead_neutralpfcanddetas[pholead_Npfcandneutralincone]/F");
  OutputTree[i]->Branch("pholead_photonpfcanddphis",&pholead_photonpfcanddphis,"pholead_photonpfcanddphis[pholead_Npfcandphotonincone]/F");
  OutputTree[i]->Branch("pholead_chargedpfcanddphis",&pholead_chargedpfcanddphis,"pholead_chargedpfcanddphis[pholead_Npfcandchargedincone]/F");
  OutputTree[i]->Branch("pholead_neutralpfcanddphis",&pholead_neutralpfcanddphis,"pholead_neutralpfcanddphis[pholead_Npfcandneutralincone]/F");
  OutputTree[i]->Branch("photrail_photonpfcandenergies",&photrail_photonpfcandenergies,"photrail_photonpfcandenergies[photrail_Npfcandphotonincone]/F");
  OutputTree[i]->Branch("photrail_chargedpfcandenergies",&photrail_chargedpfcandenergies,"photrail_chargedpfcandenergies[photrail_Npfcandchargedincone]/F");
  OutputTree[i]->Branch("photrail_neutralpfcandenergies",&photrail_neutralpfcandenergies,"photrail_neutralpfcandenergies[photrail_Npfcandneutralincone]/F");
  OutputTree[i]->Branch("photrail_photonpfcandets",&photrail_photonpfcandets,"photrail_photonpfcandets[photrail_Npfcandphotonincone]/F");
  OutputTree[i]->Branch("photrail_chargedpfcandets",&photrail_chargedpfcandets,"photrail_chargedpfcandets[photrail_Npfcandchargedincone]/F");
  OutputTree[i]->Branch("photrail_neutralpfcandets",&photrail_neutralpfcandets,"photrail_neutralpfcandets[photrail_Npfcandneutralincone]/F");
  OutputTree[i]->Branch("photrail_photonpfcanddetas",&photrail_photonpfcanddetas,"photrail_photonpfcanddetas[photrail_Npfcandphotonincone]/F");
  OutputTree[i]->Branch("photrail_chargedpfcanddetas",&photrail_chargedpfcanddetas,"photrail_chargedpfcanddetas[photrail_Npfcandchargedincone]/F");
  OutputTree[i]->Branch("photrail_neutralpfcanddetas",&photrail_neutralpfcanddetas,"photrail_neutralpfcanddetas[photrail_Npfcandneutralincone]/F");
  OutputTree[i]->Branch("photrail_photonpfcanddphis",&photrail_photonpfcanddphis,"photrail_photonpfcanddphis[photrail_Npfcandphotonincone]/F");
  OutputTree[i]->Branch("photrail_chargedpfcanddphis",&photrail_chargedpfcanddphis,"photrail_chargedpfcanddphis[photrail_Npfcandchargedincone]/F");
  OutputTree[i]->Branch("photrail_neutralpfcanddphis",&photrail_neutralpfcanddphis,"photrail_neutralpfcanddphis[photrail_Npfcandneutralincone]/F");

  }


  fHNumPU = new TH1F("NumPU_rew","NumPU_rew",50,0,50);
  fHNumPU_noweight = new TH1F("NumPU_noweight","NumPU_noweight",50,0,50);
  fHNumPUTrue = new TH1F("NumPUTrue_rew","NumPUTrue_rew",50,0,50);
  fHNumPUTrue_noweight = new TH1F("NumPUTrue_noweight","NumPUTrue_noweight",50,0,50);
  fHNumVtx = new TH1F("NumVtx_rew","NumVtx_rew",50,0,50);
	

	
  cout << "Trees and histos created" << endl;



}

void DiPhotonMiniTree::Analyze(){

  //  cout << "Analyze this event" << endl;

  if (fTR->NSuperClusters==101) return;

  //  cout << "A" << endl;

  float weight;
  if (!isdata) weight = GetPUWeight(fTR->PUnumInteractions);
  else weight=1;
  float weight3D;
  if (!isdata) weight3D = GetPUWeight3D(fTR->PUOOTnumInteractionsEarly,fTR->PUnumInteractions,fTR->PUOOTnumInteractionsLate);
  else weight3D=1;
  
  event_luminormfactor=AddWeight;

  if (!isdata) {
    fHNumPU->Fill(fTR->PUnumInteractions,weight);
    fHNumPU_noweight->Fill(fTR->PUnumInteractions);
    fHNumPUTrue->Fill(fTR->PUnumTrueInteractions,weight);
    fHNumPUTrue_noweight->Fill(fTR->PUnumTrueInteractions);
  }
  fHNumVtx->Fill(fTR->NVrtx,weight);

  //  return; // RUNNING FOR PU FILE


  // FILTERS
  // scraping veto done at ntuplizer level
  if (!PassPrimaryVertexFilter()) return;
  if (!fTR->HBHENoiseFlagIso) return;
  //  if (fTR->CSCTightHaloID) return;
  //  if (!fTR->EcalDeadTPFilterFlag) return;
  //  if (!fTR->RA2TrackingFailureFilterFlag) return;




  event_weight = weight;
  event_weight3D = weight3D;
  event_rho = fTR->Rho;
  event_sigma = fTR->Sigma;
  if (!isdata) {
    event_nPU = fTR->PUnumInteractions;
    event_PUOOTnumInteractionsEarly = fTR->PUOOTnumInteractionsEarly;
    event_PUOOTnumInteractionsLate = fTR->PUOOTnumInteractionsLate;
  }
  event_nRecVtx = fTR->NVrtx;
  event_CSCTightHaloID = fTR->CSCTightHaloID;
  event_NMuons = fTR->NMus;
  event_NMuonsTot = fTR->NMusTot;

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
  

  // sc area and scale factor calculation for isolation
  {
    float const conearea = TMath::Pi()*0.4*0.4;
    //    std::cout << "nsc " << fTR->NSuperClusters << std::endl;
    for (int i=0; i<100; i++){
      if (i<fTR->NSuperClusters){
	scarea[i] = CalculateSCArea(fTR,i); 
	scareaSF[i] = conearea/(conearea-scarea[i]);
	if (scareaSF[i]<0) std::cout << "SC area larger than isolation cone!!!" << std::endl;
      }
      else {
	scarea[i] = -999;
	scareaSF[i] = -999;
      }
    }
  }

  //  cout << "C" << endl;


  bool passtrigger = TriggerSelection();

  std::vector<int> passing_selection[18];

  bool pass[18];
  int pass12_whoissiglike;

  for (int sel_cat=0; sel_cat<18; sel_cat++){

    if (sel_cat!=10 && !passtrigger) continue; // no trigger for Zmumu selection

    if (isdata){
      if (sel_cat>=11) continue;
    }

    pass[sel_cat]=false;

    std::vector<int> passing;

    if (sel_cat==10){
      for (int i=0; i<fTR->NMus; i++){
        passing.push_back(i);
      }
      passing = MuonSelection(fTR,passing);
    }
    else {
      for (int i=0; i<fTR->NPhotons; i++){
	passing.push_back(i);
      }
      passing = PhotonPreSelection(fTR,passing);
    }


    if (sel_cat==0){
      passing = PhotonSelection(fTR,passing);
      pass[sel_cat] = StandardEventSelection(fTR,passing);
    }
    else if (sel_cat==1){
      pass[sel_cat] = SinglePhotonEventSelection(fTR,passing);
    }
    else if (sel_cat==2){
      passing = PhotonSelection(fTR,passing,"invert_sieie_cut");
      pass[sel_cat] = SinglePhotonEventSelection(fTR,passing);
    }
    else if (sel_cat==3){
      passing = PhotonSelection(fTR,passing,"revert_pixel_veto"); // revert pixel veto already done
      pass[sel_cat] = StandardEventSelection(fTR,passing);
    }
    else if (sel_cat==4){
      pass[sel_cat] = SinglePhotonEventSelection(fTR,passing);
    }
    else if (sel_cat==5){
      passing = PhotonSelection(fTR,passing);
      pass[sel_cat] = SinglePhotonEventSelection(fTR,passing);
    }
    else if (sel_cat==6){
      pass[sel_cat] = StandardEventSelection(fTR,passing);
    }
    else if (sel_cat==7){
      std::vector<int> passing_bkg = PhotonSelection(fTR,passing,"invert_sieie_cut");
      std::vector<int> newpassing;
      if (passing_bkg.size()>=1 && passing.size()>=2){
	int fondo = passing_bkg[0];
	int forcone = (passing[0]!=fondo) ? forcone = passing[0] : forcone = passing[1];
	if (fTR->PhoPt[fondo] > fTR->PhoPt[forcone]) {newpassing.push_back(fondo); newpassing.push_back(forcone); pass12_whoissiglike=1;}
	else {newpassing.push_back(forcone); newpassing.push_back(fondo); pass12_whoissiglike=0;}
      }
      passing=newpassing;
      pass[sel_cat] = StandardEventSelection(fTR,passing);
    }
    else if (sel_cat==8){
      passing = PhotonSelection(fTR,passing,"invert_sieie_cut");
      pass[sel_cat] = StandardEventSelection(fTR,passing);
    }
    else if (sel_cat==9){
      pass[sel_cat] = StandardEventSelection(fTR,passing);
    }
    else if (sel_cat==10){
      pass[sel_cat] = DiMuonFromZSelection(fTR,passing);
    }
    else if (sel_cat==11){
      passing = SignalSelection(fTR,passing);
      passing = PhotonSelection(fTR,passing);
      pass[sel_cat] = SinglePhotonEventSelection(fTR,passing);
    }
    else if (sel_cat==12){
      passing = BackgroundSelection(fTR,passing);
      passing = PhotonSelection(fTR,passing);
      pass[sel_cat] = SinglePhotonEventSelection(fTR,passing);
    }
    else if (sel_cat==13){
      passing = SignalSelection(fTR,passing);
      passing = PhotonSelection(fTR,passing);
      pass[sel_cat] = StandardEventSelection(fTR,passing);
    }
    else if (sel_cat==14){
      passing = PhotonSelection(fTR,passing);
      std::vector<int> passing_sig = SignalSelection(fTR,passing);
      std::vector<int> passing_bkg = BackgroundSelection(fTR,passing);
      std::vector<int> newpassing;
      if (passing_bkg.size()>=1 && passing_sig.size()>=1){
	int fondo = passing_bkg[0];
	int prompt = passing_sig[0];
	if (fTR->PhoPt[fondo] > fTR->PhoPt[prompt]) {newpassing.push_back(fondo); newpassing.push_back(prompt); pass12_whoissiglike=1;}
	else {newpassing.push_back(prompt); newpassing.push_back(fondo); pass12_whoissiglike=0;}
      }
      passing=newpassing;
      pass[sel_cat] = StandardEventSelection(fTR,passing);
    }
    else if (sel_cat==15){
      passing = BackgroundSelection(fTR,passing);
      passing = PhotonSelection(fTR,passing);
      pass[sel_cat] = StandardEventSelection(fTR,passing);
    }
    else if (sel_cat==16){
      std::vector<int> passing_bkg = BackgroundSelection(fTR,passing);
      passing_bkg = PhotonSelection(fTR,passing_bkg);
      std::vector<int> newpassing;
      if (passing_bkg.size()>=1 && passing.size()>=2){
	int fondo = passing_bkg[0];
	int forcone = (passing[0]!=fondo) ? forcone = passing[0] : forcone = passing[1];
	if (fTR->PhoPt[fondo] > fTR->PhoPt[forcone]) {newpassing.push_back(fondo); newpassing.push_back(forcone); pass12_whoissiglike=1;}
	else {newpassing.push_back(forcone); newpassing.push_back(fondo); pass12_whoissiglike=0;}
      }
      passing=newpassing;
      pass[sel_cat] = StandardEventSelection(fTR,passing);
    }
    else if (sel_cat==17){
      std::vector<int> passing_sig = SignalSelection(fTR,passing);
      passing_sig = PhotonSelection(fTR,passing_sig);
      std::vector<int> passing_bkg = PhotonSelection(fTR,passing,"invert_sieie_cut");
      std::vector<int> newpassing;
      if (passing_bkg.size()>=1 && passing_sig.size()>=1){
	int fondo = passing_bkg[0];
	int prompt = passing_sig[0];
	if (fTR->PhoPt[fondo] > fTR->PhoPt[prompt]) {newpassing.push_back(fondo); newpassing.push_back(prompt); pass12_whoissiglike=1;}
	else {newpassing.push_back(prompt); newpassing.push_back(fondo); pass12_whoissiglike=0;}
      }
      passing=newpassing;
      pass[sel_cat] = StandardEventSelection(fTR,passing);      
    }

    passing_selection[sel_cat] = passing;

  }



  //  cout << "D" << endl;

  for (int sel_cat=0; sel_cat<18; sel_cat++){

    if (sel_cat!=10 && !passtrigger) continue; // no trigger for Zmumu selection

    if (isdata){
      if (sel_cat>=11) continue;
    }

    if (!pass[sel_cat]) continue;

    std::vector<int> passing = passing_selection[sel_cat];
    int minsize = (is2d[sel_cat]) ? 2 : 1;

    if (passing_selection[sel_cat].size()<minsize){
      std::cout << "Error!!!" << std::endl;
      continue;
    }

    if (sel_cat==10){
      for (int i=0; i<passing.size(); i++){
	ResetVars();
	FillMuonInfo(passing.at(i));
	OutputTree[sel_cat]->Fill();
      }
    }
    else if (is2d[sel_cat]){
      ResetVars();
      FillLead(passing.at(0));
      FillTrail(passing.at(1));
      bool dofill=true;

      if (sel_cat==7 || sel_cat==14 || sel_cat==16 || sel_cat==17) event_pass12whoissiglike=pass12_whoissiglike;

      if (sel_cat==6 || ((sel_cat==7 || sel_cat==16) && pass12_whoissiglike==0)) {
	  isolations_struct rcone_isos;
	  rcone_isos = RandomConeIsolation(fTR,passing.at(0),"");
	  pholead_pho_Cone04PhotonIso_dEta015EB_dR070EE_mvVtx = rcone_isos.photon;
	  pholead_pho_Cone04NeutralHadronIso_mvVtx = rcone_isos.neutral;
	  pholead_pho_Cone04ChargedHadronIso_dR02_dz02_dxy01 = rcone_isos.charged;
	  pholead_pho_Cone04PFCombinedIso = rcone_isos.photon+rcone_isos.neutral+rcone_isos.charged;
	  pholead_Npfcandphotonincone = rcone_isos.nphotoncand;
	  pholead_Npfcandchargedincone = rcone_isos.nchargedcand;
	  pholead_Npfcandneutralincone = rcone_isos.nneutralcand;
	  for (int i=0; i<pholead_Npfcandphotonincone && i<global_size_pfcandarrays; i++) pholead_photonpfcandenergies[i] = rcone_isos.photoncandenergies.at(i);
	  for (int i=0; i<pholead_Npfcandchargedincone && i<global_size_pfcandarrays; i++) pholead_chargedpfcandenergies[i] = rcone_isos.chargedcandenergies.at(i);
	  for (int i=0; i<pholead_Npfcandneutralincone && i<global_size_pfcandarrays; i++) pholead_neutralpfcandenergies[i] = rcone_isos.neutralcandenergies.at(i);
	  for (int i=0; i<pholead_Npfcandphotonincone && i<global_size_pfcandarrays; i++) pholead_photonpfcandets[i] = rcone_isos.photoncandets.at(i);
	  for (int i=0; i<pholead_Npfcandchargedincone && i<global_size_pfcandarrays; i++) pholead_chargedpfcandets[i] = rcone_isos.chargedcandets.at(i);
	  for (int i=0; i<pholead_Npfcandneutralincone && i<global_size_pfcandarrays; i++) pholead_neutralpfcandets[i] = rcone_isos.neutralcandets.at(i);
	  for (int i=0; i<pholead_Npfcandphotonincone && i<global_size_pfcandarrays; i++) pholead_photonpfcanddetas[i] = rcone_isos.photoncanddetas.at(i);
	  for (int i=0; i<pholead_Npfcandchargedincone && i<global_size_pfcandarrays; i++) pholead_chargedpfcanddetas[i] = rcone_isos.chargedcanddetas.at(i);
	  for (int i=0; i<pholead_Npfcandneutralincone && i<global_size_pfcandarrays; i++) pholead_neutralpfcanddetas[i] = rcone_isos.neutralcanddetas.at(i);
	  for (int i=0; i<pholead_Npfcandphotonincone && i<global_size_pfcandarrays; i++) pholead_photonpfcanddphis[i] = rcone_isos.photoncanddphis.at(i);
	  for (int i=0; i<pholead_Npfcandchargedincone && i<global_size_pfcandarrays; i++) pholead_chargedpfcanddphis[i] = rcone_isos.chargedcanddphis.at(i);
	  for (int i=0; i<pholead_Npfcandneutralincone && i<global_size_pfcandarrays; i++) pholead_neutralpfcanddphis[i] = rcone_isos.neutralcanddphis.at(i);
	  if (rcone_isos.photon==-999 || rcone_isos.neutral==-999 || rcone_isos.charged==-999) dofill=false;
	}
      if (sel_cat==6 || ((sel_cat==7 || sel_cat==16) && pass12_whoissiglike==1)) {
	  isolations_struct rcone_isos;
	  rcone_isos = RandomConeIsolation(fTR,passing.at(1),"");
	  photrail_pho_Cone04PhotonIso_dEta015EB_dR070EE_mvVtx = rcone_isos.photon;
	  photrail_pho_Cone04NeutralHadronIso_mvVtx = rcone_isos.neutral;
	  photrail_pho_Cone04ChargedHadronIso_dR02_dz02_dxy01 = rcone_isos.charged;
	  photrail_pho_Cone04PFCombinedIso = rcone_isos.photon+rcone_isos.neutral+rcone_isos.charged;
	  photrail_Npfcandphotonincone = rcone_isos.nphotoncand;
	  photrail_Npfcandchargedincone = rcone_isos.nchargedcand;
	  photrail_Npfcandneutralincone = rcone_isos.nneutralcand;
	  for (int i=0; i<photrail_Npfcandphotonincone && i<global_size_pfcandarrays; i++) photrail_photonpfcandenergies[i] = rcone_isos.photoncandenergies.at(i);
	  for (int i=0; i<photrail_Npfcandchargedincone && i<global_size_pfcandarrays; i++) photrail_chargedpfcandenergies[i] = rcone_isos.chargedcandenergies.at(i);
	  for (int i=0; i<photrail_Npfcandneutralincone && i<global_size_pfcandarrays; i++) photrail_neutralpfcandenergies[i] = rcone_isos.neutralcandenergies.at(i);
	  for (int i=0; i<photrail_Npfcandphotonincone && i<global_size_pfcandarrays; i++) photrail_photonpfcandets[i] = rcone_isos.photoncandets.at(i);
	  for (int i=0; i<photrail_Npfcandchargedincone && i<global_size_pfcandarrays; i++) photrail_chargedpfcandets[i] = rcone_isos.chargedcandets.at(i);
	  for (int i=0; i<photrail_Npfcandneutralincone && i<global_size_pfcandarrays; i++) photrail_neutralpfcandets[i] = rcone_isos.neutralcandets.at(i);
	  for (int i=0; i<photrail_Npfcandphotonincone && i<global_size_pfcandarrays; i++) photrail_photonpfcanddetas[i] = rcone_isos.photoncanddetas.at(i);
	  for (int i=0; i<photrail_Npfcandchargedincone && i<global_size_pfcandarrays; i++) photrail_chargedpfcanddetas[i] = rcone_isos.chargedcanddetas.at(i);
	  for (int i=0; i<photrail_Npfcandneutralincone && i<global_size_pfcandarrays; i++) photrail_neutralpfcanddetas[i] = rcone_isos.neutralcanddetas.at(i);
	  for (int i=0; i<photrail_Npfcandphotonincone && i<global_size_pfcandarrays; i++) photrail_photonpfcanddphis[i] = rcone_isos.photoncanddphis.at(i);
	  for (int i=0; i<photrail_Npfcandchargedincone && i<global_size_pfcandarrays; i++) photrail_chargedpfcanddphis[i] = rcone_isos.chargedcanddphis.at(i);
	  for (int i=0; i<photrail_Npfcandneutralincone && i<global_size_pfcandarrays; i++) photrail_neutralpfcanddphis[i] = rcone_isos.neutralcanddphis.at(i);
	  if (rcone_isos.photon==-999 || rcone_isos.neutral==-999 || rcone_isos.charged==-999) dofill=false;
	}

      float invmass0 = (CorrPhoton(fTR,passing.at(0),0)+CorrPhoton(fTR,passing.at(1),0)).M();
      float invmass5 = (CorrPhoton(fTR,passing.at(0),5)+CorrPhoton(fTR,passing.at(1),5)).M();
      float invmass6 = (CorrPhoton(fTR,passing.at(0),6)+CorrPhoton(fTR,passing.at(1),6)).M();
      dipho_mgg_photon = invmass0;
      dipho_mgg_newCorr = invmass5;
      dipho_mgg_newCorrLocal = invmass6;
      if (dofill) OutputTree[sel_cat]->Fill();
    }

    else if (!is2d[sel_cat]){

      for (int i=0; i<passing.size(); i++){
      ResetVars();
      FillLead(passing.at(i));
      bool dofill = true;

      if (sel_cat==1) {
	isolations_struct rcone_isos;
	rcone_isos = RandomConeIsolation(fTR,passing.at(i),"");
	pholead_pho_Cone04PhotonIso_dEta015EB_dR070EE_mvVtx = rcone_isos.photon;
	pholead_pho_Cone04NeutralHadronIso_mvVtx = rcone_isos.neutral;
	pholead_pho_Cone04ChargedHadronIso_dR02_dz02_dxy01 = rcone_isos.charged;
	pholead_pho_Cone04PFCombinedIso = rcone_isos.photon+rcone_isos.neutral+rcone_isos.charged;
	pholead_Npfcandphotonincone = rcone_isos.nphotoncand;
	pholead_Npfcandchargedincone = rcone_isos.nchargedcand;
	pholead_Npfcandneutralincone = rcone_isos.nneutralcand;
	for (int i=0; i<pholead_Npfcandphotonincone && i<global_size_pfcandarrays; i++) pholead_photonpfcandenergies[i] = rcone_isos.photoncandenergies.at(i);
	for (int i=0; i<pholead_Npfcandchargedincone && i<global_size_pfcandarrays; i++) pholead_chargedpfcandenergies[i] = rcone_isos.chargedcandenergies.at(i);
	for (int i=0; i<pholead_Npfcandneutralincone && i<global_size_pfcandarrays; i++) pholead_neutralpfcandenergies[i] = rcone_isos.neutralcandenergies.at(i);
	for (int i=0; i<pholead_Npfcandphotonincone && i<global_size_pfcandarrays; i++) pholead_photonpfcandets[i] = rcone_isos.photoncandets.at(i);
	for (int i=0; i<pholead_Npfcandchargedincone && i<global_size_pfcandarrays; i++) pholead_chargedpfcandets[i] = rcone_isos.chargedcandets.at(i);
	for (int i=0; i<pholead_Npfcandneutralincone && i<global_size_pfcandarrays; i++) pholead_neutralpfcandets[i] = rcone_isos.neutralcandets.at(i);
	for (int i=0; i<pholead_Npfcandphotonincone && i<global_size_pfcandarrays; i++) pholead_photonpfcanddetas[i] = rcone_isos.photoncanddetas.at(i);
	for (int i=0; i<pholead_Npfcandchargedincone && i<global_size_pfcandarrays; i++) pholead_chargedpfcanddetas[i] = rcone_isos.chargedcanddetas.at(i);
	for (int i=0; i<pholead_Npfcandneutralincone && i<global_size_pfcandarrays; i++) pholead_neutralpfcanddetas[i] = rcone_isos.neutralcanddetas.at(i);
	for (int i=0; i<pholead_Npfcandphotonincone && i<global_size_pfcandarrays; i++) pholead_photonpfcanddphis[i] = rcone_isos.photoncanddphis.at(i);
	for (int i=0; i<pholead_Npfcandchargedincone && i<global_size_pfcandarrays; i++) pholead_chargedpfcanddphis[i] = rcone_isos.chargedcanddphis.at(i);
	for (int i=0; i<pholead_Npfcandneutralincone && i<global_size_pfcandarrays; i++) pholead_neutralpfcanddphis[i] = rcone_isos.neutralcanddphis.at(i);
	if (rcone_isos.photon==-999 || rcone_isos.neutral==-999 || rcone_isos.charged==-999) dofill=false;
      }

      if (dofill) OutputTree[sel_cat]->Fill();
      }

    }

  }
 

};

void DiPhotonMiniTree::End(){
  fOutputFile->cd();
  for (int i=0; i<18; i++) OutputTree[i]->Write();	
  fHNumPU->Write();
  fHNumPU_noweight->Write();
  fHNumPUTrue->Write();
  fHNumPUTrue_noweight->Write();
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
    if (fTR->PhoHasPixSeed[*it]!=wantpixelseed) it=passing.erase(it); else it++;
    //    if (fTR->PhoPassConvSafeElectronVeto[*it]==wantpixelseed) it=passing.erase(it); else it++;
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
    float sieie=SieieRescale(fTR->PhoSigmaIetaIeta[*it],(bool)(fabs(eta)<1.4442));
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

//  XXX EM ENRICHMENT PRESELECTION NOT APPLIED
//  for (vector<int>::iterator it = passing.begin(); it != passing.end(); ){ // isolation cuts (filter)
//    float r9=fTR->SCR9[fTR->PhotSCindex[*it]];
//    bool pass=0;
//    if(fTR->pho_Cone02ChargedHadronIso_dR02_dz02_dxy01[*it]<4) pass=1;
//    if (!pass) it=passing.erase(it); else it++;
//  }

  return passing;

};

std::vector<int> DiPhotonMiniTree::MuonSelection(TreeReader *fTR, std::vector<int> passing){

  for (vector<int>::iterator it = passing.begin(); it != passing.end(); ){
    float eta=fTR->MuEta[*it];
    if ((fabs(eta)>1.4442 && fabs(eta)<1.56) || (fabs(eta)>2.5)) it=passing.erase(it); else it++;
  }

  for (vector<int>::iterator it = passing.begin(); it != passing.end(); ){
    if (fTR->MuPt[*it]<10) it=passing.erase(it); else it++;
  }

  for (vector<int>::iterator it = passing.begin(); it != passing.end(); ){
    if (!(fTR->MuIsGlobalMuon[*it] || fTR->MuIsTrackerMuon[*it])) it=passing.erase(it); else it++;
  }

//  for (vector<int>::iterator it = passing.begin(); it != passing.end(); ){
//    if (fTR->MuRelIso03[*it]>0.15) it=passing.erase(it); else it++;
//  }

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

std::vector<int> DiPhotonMiniTree::PhotonSelection(TreeReader *fTR, std::vector<int> passing, TString mode){

  if (mode!="" && mode!="invert_sieie_cut" && mode!="revert_pixel_veto"){
    std::cout << "Error" << std::endl;
    return std::vector<int>();
  }

  if (mode=="revert_pixel_veto") passing = ApplyPixelVeto(fTR,passing,1); // for DY pixel veto reverse
  else passing = ApplyPixelVeto(fTR,passing,0);


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
    float sieie=SieieRescale(fTR->PhoSigmaIetaIeta[*it],(bool)(fabs(eta)<1.4442));
    bool pass=0;
    if (mode=="invert_sieie_cut"){ // sieie sideband
      if (fabs(eta)<1.4442 && sieie>0.011 && sieie<0.014) pass=1;
      if (fabs(eta)>1.56 && sieie>0.030 && sieie<0.034) pass=1;
    }
    else {
    if (fabs(eta)<1.4442 && sieie<0.011) pass=1;
    if (fabs(eta)>1.56 && sieie<0.030) pass=1;
    }
    if (!pass) it=passing.erase(it); else it++;
  }

//  for (vector<int>::iterator it = passing.begin(); it != passing.end(); ){
//    bool pass=1;
//    float eta=fTR->SCEta[fTR->PhotSCindex[*it]];
//    float phi=fTR->SCPhi[fTR->PhotSCindex[*it]];
//    for (int i=0; i<fTR->NMus; i++){
//      float mueta = fTR->MuEta[i];
//      float muphi = fTR->MuPhi[i];
//      if (Util::GetDeltaR(mueta,eta,muphi,phi)<0.4) pass=0;
//    }
//    if (!pass) it=passing.erase(it); else it++;
//  }

  for (vector<int>::iterator it = passing.begin(); it != passing.end(); ){
    bool pass=0;
    float eta=fTR->SCEta[fTR->PhotSCindex[*it]];
    float combiso = PFIsolation(*it,0,"combined");
    if (mode=="no_combiso_cut") pass=1; // pass in any case
    else if (mode=="cut_combiso_sideband"){ // selection for sideband
      //      if (combiso> && combiso<) pass=1;
    }
    else if (combiso<999) pass=1;
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
    if (!(fTR->PhoMCmatchexitcode[*it]==1 || fTR->PhoMCmatchexitcode[*it]==2)) {
      pass=1;
    }
    else {
      if (fTR->GenPhotonIsoDR04[fTR->PhoMCmatchindex[*it]]>5) pass=1;
    }
    if (!pass) it=passing.erase(it); else it++;
  }

  return passing;

};

std::vector<int> DiPhotonMiniTree::DiPhotonInvariantMassCutSelection(TreeReader *fTR, std::vector<int> passing){

  for (vector<int>::iterator it = passing.begin(); it != passing.end(); ){
    bool pass=0;
    if (fTR->PhoPt[*it]>25) pass=1;
    if (!pass) it=passing.erase(it); else it++;
  }

  if (passing.size()<2) return std::vector<int>();

  passing.resize(2);

  float invmass0 = (CorrPhoton(fTR,passing.at(0),0)+CorrPhoton(fTR,passing.at(1),0)).M();
  if (fabs(invmass0-91.2)>10) return std::vector<int>();

  return passing;

};

bool DiPhotonMiniTree::SinglePhotonEventSelection(TreeReader *fTR, std::vector<int> &passing){

  for (vector<int>::iterator it = passing.begin(); it != passing.end(); ){
    bool pass=0;
    if (fTR->PhoPt[*it]>25) pass=1;
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
  float dphi=Util::DeltaPhi(fTR->PhoPhi[passing.at(0)],fTR->PhoPhi[passing.at(1)]);
  double dR=sqrt(dphi*dphi+deta*deta);

  if (fTR->PhoPt[passing.at(0)]<40) return false;
  if (fTR->PhoPt[passing.at(1)]<25) return false;
  //  if (invmass0<80) return false;
  if (dR<0.45) return false;

  return true;
};

bool DiPhotonMiniTree::DiMuonFromZSelection(TreeReader *fTR, std::vector<int> &passing){

  if (passing.size()<2) return false;

  passing.resize(2); // keep only the first two

  TLorentzVector mu1(fTR->MuPx[passing.at(0)],fTR->MuPy[passing.at(0)],fTR->MuPz[passing.at(0)],fTR->MuE[passing.at(0)]);
  TLorentzVector mu2(fTR->MuPx[passing.at(1)],fTR->MuPy[passing.at(1)],fTR->MuPz[passing.at(1)],fTR->MuE[passing.at(1)]);

  float invmass0 = (mu1+mu2).M();

  if (fabs(invmass0-91.2)>10) return false;

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


double DiPhotonMiniTree::phiNorm(float phi) {

  const float pi = 3.1415927;
  const float twopi = 2.0*pi;

  if(phi >  pi) {phi = phi - twopi;}
  if(phi < -pi) {phi = phi + twopi;}

  return phi;
}

bool DiPhotonMiniTree::FindCloseJetsAndPhotons(TreeReader *fTR, float rotation_phi, int phoqi, TString mod){

  if (mod!="" && mod!="nocombisocut") {std::cout << "error" << std::endl; return true;}

  TVector3 photon_position = TVector3(fTR->SCx[fTR->PhotSCindex[phoqi]],fTR->SCy[fTR->PhotSCindex[phoqi]],fTR->SCz[fTR->PhotSCindex[phoqi]]);
  if (rotation_phi!=0) {
    TRotation r; r.RotateZ(rotation_phi);
    photon_position *= r;
  }
  double eta = photon_position.Eta();
  double phi = photon_position.Phi();
  
  const bool debug=false;
  if (debug) std::cout << "calling FindCloseJetsAndPhotons eta=" << eta << " phi=" << phi << std::endl;

  const float mindR = 0.8;
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

  if (mod!="nocombisocut") { if (PFIsolation(phoqi,rotation_phi,"combined")>999) found=true; }

  for (int i=0; i<fTR->NMus; i++){
    float mueta = fTR->MuEta[i];
    float muphi = fTR->MuPhi[i];
    if (Util::GetDeltaR(mueta,eta,muphi,phi)<0.4) found=true;
  }

  if (debug) std::cout << "returning " << found << std::endl;
  return found;

};

std::vector<int> DiPhotonMiniTree::GetPFCandIDedRemovals(TreeReader *fTR, int phoqi){
  std::vector<int> out;
  if (fTR->Pho_isPFPhoton[phoqi]) out.push_back(fTR->pho_matchedPFPhotonCand[phoqi]);
  return out;
};

std::vector<int> DiPhotonMiniTree::GetPFCandInsideFootprint(TreeReader *fTR, int phoqi, float rotation_phi, TString component){
  return GetPFCandWithFootprintRemoval(fTR,phoqi,rotation_phi,false,component);
};

std::vector<int> DiPhotonMiniTree::GetPFCandWithFootprintRemoval(TreeReader *fTR, int phoqi, float rotation_phi, bool outoffootprint, TString component){

  if (component!="photon"){
    std::cout << "propagation not implemented for non photon objects!!!" << std::endl;
    return std::vector<int>();
  }

  if (component!="neutral" && component!="charged" && component!="photon" && component!="combined") {
    std::cout << "Wrong choice for component" << std::endl;
    return std::vector<int>();
  }

  int scindex = fTR->PhotSCindex[phoqi];
  
  if (scindex<0) {
    std::cout << "Error in GetPFCandOverlappingSC" << std::endl;
    std::cout << scindex << " " << phoqi << std::endl;
    return std::vector<int>();
  }

  bool isbarrel = fTR->PhoisEB[phoqi];
  int nxtals = fTR->SCNXtals[scindex];

  std::vector<int> result;

  for (int i=0; i<fTR->NPfCand; i++){

    int type = FindPFCandType(fTR->PfCandPdgId[i]);
    if (!(type==0 || type==1 || type==2)) continue;

    if (component=="neutral" && type!=0) continue;
    if (component=="charged" && type!=1) continue;
    if (component=="photon" && type!=2) continue;

    TVector3 sc_position = TVector3(fTR->SCx[scindex],fTR->SCy[scindex],fTR->SCz[scindex]);
    //    TVector3 ecalpfhit = PropagatePFCandToEcal(i,isbarrel ? sc_position.Perp() : sc_position.z(),isbarrel); // approximate, but much faster! Note that inputs are not changed by rotation

    bool inside=false;

    if (fTR->Pho_isPFPhoton[phoqi] && fTR->pho_matchedPFPhotonCand[phoqi]==i) continue;

    for (int j=0; j<nxtals; j++){
      
      TVector3 xtal_position = TVector3(fTR->SCxtalX[scindex][j],fTR->SCxtalY[scindex][j],fTR->SCxtalZ[scindex][j]);

      if (rotation_phi!=0) {
	TRotation r; r.RotateZ(rotation_phi);
	xtal_position *= r;
      }

      TVector3 ecalpfhit = PropagatePFCandToEcal(i,isbarrel ? xtal_position.Perp() : xtal_position.z(), isbarrel); // this would be the most correct

      if (isbarrel){
	float xtalEtaWidth = fTR->SCxtalEtaWidth[scindex][j]*(1+global_linkbyrechit_enlargement);
	float xtalPhiWidth = fTR->SCxtalPhiWidth[scindex][j]*(1+global_linkbyrechit_enlargement);
	if (fabs(ecalpfhit.Eta()-xtal_position.Eta())<xtalEtaWidth/2 && Util::DeltaPhi(ecalpfhit.Phi(),xtal_position.Phi())<xtalPhiWidth/2) inside=true;
      }
      else { // EE
	if (ecalpfhit.z()*xtal_position.z()>0){
	  TVector3 xtal_corners[4];
	  for (int k=0; k<4; k++) xtal_corners[k] = TVector3(fTR->SCxtalfrontX[scindex][j][k],fTR->SCxtalfrontY[scindex][j][k],fTR->SCxtalfrontZ[scindex][j][k]);
	  if (rotation_phi!=0) {
	    TRotation r; r.RotateZ(rotation_phi);
	    for (int k=0; k<4; k++) xtal_corners[k] *= r;
	  }
	  float hitx = ecalpfhit.x();
	  float hity = ecalpfhit.y();
	  float polx[5];
	  float poly[5];
	  for (int k=0; k<4; k++) polx[k] = xtal_corners[k].x();
	  for (int k=0; k<4; k++) poly[k] = xtal_corners[k].y();
	  polx[4]=polx[0]; poly[4]=poly[0]; // closed polygon
	  float centerx = (polx[0]+polx[1]+polx[2]+polx[3])/4;
	  float centery = (poly[0]+poly[1]+poly[2]+poly[3])/4;
	  hitx = centerx + (hitx-centerx)/(1.0+global_linkbyrechit_enlargement);
	  hity = centery + (hity-centery)/(1.0+global_linkbyrechit_enlargement);
	  if (TMath::IsInside(hitx,hity,5,polx,poly)) inside=true;
	}
      }

    }

    if (outoffootprint) inside=!inside;
    if (inside) result.push_back(i);

  }

  return result;

};

TVector3 DiPhotonMiniTree::PropagatePFCandToEcal(int pfcandindex, float position, bool isbarrel){
  // WARNING: this propagates until EE+ or EE- at the given TMath::Abs(position.z()) for isbarrel=0, depending on where the candidate is pointing.

  int i = pfcandindex;

  if (FindPFCandType(fTR->PfCandPdgId[i])!=2) {
    std::cout << "propagation not implemented for non photon objects!!!" << std::endl;
    return TVector3(0,0,0);
  }

  TVector3 pfvertex(fTR->PfCandVx[i],fTR->PfCandVy[i],fTR->PfCandVz[i]);
  TVector3 pfmomentum(fTR->PfCandPx[i],fTR->PfCandPy[i],fTR->PfCandPz[i]);
  pfmomentum = pfmomentum.Unit();
  TVector3 ecalpfhit(0,0,0);
  if (isbarrel){
    double p[3] = {pfvertex.x(),pfvertex.y(),pfvertex.z()};
    double d[3] = {pfmomentum.x(),pfmomentum.y(),pfmomentum.z()};
    double dist = TGeoTube::DistFromInsideS(p,d,0,position,1e+10);
    ecalpfhit = pfvertex + dist*pfmomentum;
  }
  else { // EE
    double dim[6]={1e+10,1e+10,fabs(position),0,0,0};
    double p[3] = {pfvertex.x(),pfvertex.y(),pfvertex.z()};
    double d[3] = {pfmomentum.x(),pfmomentum.y(),pfmomentum.z()};
    eegeom.SetDimensions(dim);
    double dist = eegeom.DistFromInside(p,d);
    ecalpfhit = pfvertex + dist*pfmomentum;
  }

  return ecalpfhit;

};

isolations_struct DiPhotonMiniTree::RandomConeIsolation(TreeReader *fTR, int phoqi, TString mod){

  float result=0;
  const double pi = TMath::Pi();

  double rotation_phi = pi/2;

  bool isok = !(FindCloseJetsAndPhotons(fTR,rotation_phi,phoqi,mod));
  if (!isok) {
    rotation_phi = -pi/2;
    isok=!(FindCloseJetsAndPhotons(fTR,rotation_phi,phoqi,mod));
  }
  
  int count=0;
  while (!isok && count<20) {
    rotation_phi = randomgen->Uniform(0.8,2*pi-0.8);
    isok=!(FindCloseJetsAndPhotons(fTR,rotation_phi,phoqi,mod));
    count++;
  }

  isolations_struct out;
  out.nphotoncand=0; out.nchargedcand=0; out.nneutralcand=0;

  if (count==20){
    std::cout << "Error in random cone generation!!!"  << std::endl;
    out.photon = -999;
    out.charged = -999;
    out.neutral = -999;
    return out;
  };


  out.photon = PFIsolation(phoqi,rotation_phi,"photon",&(out.nphotoncand),&(out.photoncandenergies),&(out.photoncandets),&(out.photoncanddetas),&(out.photoncanddphis));
  out.charged = PFIsolation(phoqi,rotation_phi,"charged",&(out.nchargedcand),&(out.chargedcandenergies),&(out.chargedcandets),&(out.chargedcanddetas),&(out.chargedcanddphis));
  out.neutral = PFIsolation(phoqi,rotation_phi,"neutral",&(out.nneutralcand),&(out.neutralcandenergies),&(out.neutralcandets),&(out.neutralcanddetas),&(out.neutralcanddphis));
  return out;

};

isolations_struct DiPhotonMiniTree::PFConeIsolation(TreeReader *fTR, int phoqi){
  isolations_struct out;
  out.nphotoncand=0; out.nchargedcand=0; out.nneutralcand=0;
  out.photon = PFIsolation(phoqi,0,"photon",&(out.nphotoncand),&(out.photoncandenergies),&(out.photoncandets),&(out.photoncanddetas),&(out.photoncanddphis));
  out.charged = PFIsolation(phoqi,0,"charged",&(out.nchargedcand),&(out.chargedcandenergies),&(out.chargedcandets),&(out.chargedcanddetas),&(out.chargedcanddphis));
  out.neutral = PFIsolation(phoqi,0,"neutral",&(out.nneutralcand),&(out.neutralcandenergies),&(out.neutralcandets),&(out.neutralcanddetas),&(out.neutralcanddphis));
  return out;
};

float DiPhotonMiniTree::PFIsolation(int phoqi, float rotation_phi, TString component, int *counter, std::vector<float> *energies, std::vector<float> *ets, std::vector<float> *detas, std::vector<float> *dphis, std::vector<int> removals){

  if (component=="combined") return PFIsolation(phoqi,rotation_phi,"photon",counter,energies,ets,detas,dphis,removals) \
			       + PFIsolation(phoqi,rotation_phi,"charged",counter,energies,ets,detas,dphis,removals) \
			       + PFIsolation(phoqi,rotation_phi,"neutral",counter,energies,ets,detas,dphis,removals);

  if (component!="neutral" && component!="charged" && component!="photon") {std::cout << "wrong" << std::endl; return -999;} 

  float minimal_pfphotoncand_threshold_EB = (component=="photon") ? global_minthrpfphotoncandEB : 0.0;
  float minimal_pfphotoncand_threshold_EE = (component=="photon") ? global_minthrpfphotoncandEE : 0.0;

  // footprint removal only for photons, charged use Poter's veto cones, nothing for neutrals

  float result=0;
  float scaleresult=1;
  
  TVector3 photon_position = TVector3(fTR->SCx[fTR->PhotSCindex[phoqi]],fTR->SCy[fTR->PhotSCindex[phoqi]],fTR->SCz[fTR->PhotSCindex[phoqi]]);

  if (rotation_phi!=0) {
    TRotation r; r.RotateZ(rotation_phi);
    photon_position *= r;
  } 

  if (component=="photon"){
    std::vector<int> footprint = GetPFCandInsideFootprint(fTR,phoqi,rotation_phi,"photon");
    for (int i=0; i<footprint.size(); i++) removals.push_back(footprint.at(i));
    scaleresult = scareaSF[fTR->PhotSCindex[phoqi]];
  }

  for (int i=0; i<fTR->NPfCand; i++){

    if (fTR->Pho_isPFPhoton[phoqi] && fTR->pho_matchedPFPhotonCand[phoqi]==i) continue;
 
    int type = FindPFCandType(fTR->PfCandPdgId[i]);

    if (!(type==0 || type==1 || type==2)) continue;

    if (component=="neutral" && type!=0) continue;
    if (component=="charged" && type!=1) continue;
    if (component=="photon" && type!=2) continue;

    bool removed = false;
    for (int j=0; j<removals.size(); j++) {
      if (i==removals.at(j)) removed=true;
    }
    if (removed) continue;

    double pt = fTR->PfCandPt[i];
    double energy = fTR->PfCandEnergy[i];

    double dR = 999;
    double dEta = 999;
    double dPhi = 999;

    if (type==2){
      angular_distances_struct angles = GetPFCandDeltaRFromSC(fTR,phoqi,i,rotation_phi);
      dR = angles.dR;
      dEta = angles.dEta;
      dPhi = angles.dPhi;
      if (dR>0.4) continue;
    }
    else if (type==0 || type==1) {

      TVector3 phovtx(fTR->PhoVx[phoqi],fTR->PhoVy[phoqi],fTR->PhoVz[phoqi]);
      TVector3 pfvertex(fTR->PfCandVx[i],fTR->PfCandVy[i],fTR->PfCandVz[i]);
      TVector3 photon_direction = photon_position-pfvertex;
      float sceta = photon_direction.Eta();
      float scphi = photon_direction.Phi();
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


      dEta = fTR->PfCandEta[i] - sceta;
      dPhi = DeltaPhiSigned(fTR->PfCandPhi[i],scphi);
      dR = sqrt(dEta*dEta+dPhi*dPhi);
      
      if (dR>0.4) continue;



//      if (type==2){ 
//	if (fTR->PhoisEB[phoqi]){
//	  if (fabs(dEta)<0.015) continue;
//	}
//	else if (fTR->PhoisEE[phoqi]){
//	  float limit_dR = 0.00864*fabs(sinh(sceta))*4;
//	  if (dR<limit_dR) continue;
//	}
//	else {
//	  std::cout << "Something wrong" << std::endl;
//	  return -999;
//	}
//      }

      if (type==1){
	if (dz>0.2) continue;
	if (dxy>0.1) continue;
	if (dR<0.02) continue;
	scaleresult = 1+2.5e-3;
      }

    }
    else {
      std::cout << "something wrong" << std::endl;
      continue;
    }

    // pfcandidate threshold
    if (type==2){ 
      if (fTR->PhoisEB[phoqi]){
	if (pt<minimal_pfphotoncand_threshold_EB) continue;
      }
      else if (fTR->PhoisEE[phoqi]){
	if (pt<minimal_pfphotoncand_threshold_EE) continue;
      }
      else {
	std::cout << "Something wrong" << std::endl;
	return -999;
      }
    }

    if (counter) (*counter)++;
    if (energies) energies->push_back(energy);
    if (ets) ets->push_back(pt);
    if (detas) detas->push_back(dEta);
    if (dphis) dphis->push_back(dPhi);
    result+=pt;

  } // end pf cand loop


  return result*scaleresult-GetPUEnergy(fTR,component,fTR->PhoisEB[phoqi]);

};

float DiPhotonMiniTree::PFPhotonIsolationAroundMuon(int muqi, int *counter, std::vector<float> *energies, std::vector<float> *ets, std::vector<float> *detas, std::vector<float> *dphis){
  
  float minimal_pfphotoncand_threshold_EB = global_minthrpfphotoncandEB;
  float minimal_pfphotoncand_threshold_EE = global_minthrpfphotoncandEE;

  float result=0;
  
  TLorentzVector mu(fTR->MuPx[muqi],fTR->MuPy[muqi],fTR->MuPz[muqi],fTR->MuE[muqi]);
  bool isbarrel = (fabs(mu.Eta())<1.4442);

  for (int i=0; i<fTR->NPfCand; i++){
    
    int type = FindPFCandType(fTR->PfCandPdgId[i]);
    if (type!=2) continue;
    
    double pt = fTR->PfCandPt[i];
    double energy = fTR->PfCandEnergy[i];
    double dEta = fTR->PfCandEta[i] - mu.Eta();
    double dPhi = DeltaPhiSigned(fTR->PfCandPhi[i],mu.Phi());
    double dR = sqrt(dEta*dEta+dPhi*dPhi);
  
    if (dR>0.4) continue;
    if (dR<0.1) continue;
  
  
    // pfcandidate threshold
    if (isbarrel){
      if (pt<minimal_pfphotoncand_threshold_EB) continue;
    }
    else {
      if (pt<minimal_pfphotoncand_threshold_EE) continue;
    }
    
    if (counter) (*counter)++;
    if (energies) energies->push_back(energy);
    if (ets) ets->push_back(pt);
    if (detas) detas->push_back(dEta);
    if (dphis) dphis->push_back(dPhi);

    result+=pt;
    
  } // end pf cand loop

  const float scaleresult = 1.066;      
  return result*scaleresult;
  
};


angular_distances_struct DiPhotonMiniTree::GetPFCandDeltaRFromSC(TreeReader *fTR, int phoqi, int pfindex, float rotation_phi){

  int i = pfindex;
  int scindex = fTR->PhotSCindex[phoqi];

  if (scindex<0) {
    std::cout << "Error in GetPFCandDeltaRFromSC" << std::endl;
    angular_distances_struct out;
    out.dR = 999;
    out.dEta = 999;
    out.dPhi = 999;
    return out;
  }

  bool isbarrel = fTR->PhoisEB[phoqi];

  TVector3 sc_position = TVector3(fTR->SCx[fTR->PhotSCindex[phoqi]],fTR->SCy[fTR->PhotSCindex[phoqi]],fTR->SCz[fTR->PhotSCindex[phoqi]]);

  if (rotation_phi!=0) {
    TRotation r; r.RotateZ(rotation_phi);
    sc_position *= r;
  }

  TVector3 pfvertex(fTR->PfCandVx[i],fTR->PfCandVy[i],fTR->PfCandVz[i]);
  TVector3 pfmomentum(fTR->PfCandPx[i],fTR->PfCandPy[i],fTR->PfCandPz[i]);
  pfmomentum = pfmomentum.Unit();

  TVector3 ecalpfhit = PropagatePFCandToEcal(i,isbarrel ? sc_position.Perp() : sc_position.z(),isbarrel);

  angular_distances_struct out;
  out.dR = Util::GetDeltaR(sc_position.Eta(),ecalpfhit.Eta(),sc_position.Phi(),ecalpfhit.Phi());
  out.dEta = ecalpfhit.Eta()-sc_position.Eta();
  out.dPhi = DeltaPhiSigned(ecalpfhit.Phi(),sc_position.Phi());

  return out;

}

int DiPhotonMiniTree::FindPFCandType(int id){

  int type = -1;

  if (id==111 || id==130 || id==310 || id==2112) type=0; //neutral hadrons
  if (fabs(id)==211 || fabs(id)==321 || id==999211 || fabs(id)==2212) type=1; //charged hadrons
  if (id==22) type=2; //photons
  if (fabs(id)==11) type=3; //electrons
  if (fabs(id)==13) type=4; //muons

  return type;
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
  pholead_sieie = SieieRescale(fTR->PhoSigmaIetaIeta[index],(bool)(fabs(fTR->SCEta[fTR->PhotSCindex[index]])<1.4442));
  pholead_hoe = fTR->PhoHoverE[index];
  pholead_brem = fTR->SCBrem[fTR->PhotSCindex[index]];
  pholead_sigmaPhi = fTR->SCPhiWidth[fTR->PhotSCindex[index]];
  pholead_sigmaEta = fTR->SCEtaWidth[fTR->PhotSCindex[index]];
  pholead_pho_Cone01PhotonIso_dEta015EB_dR070EE_mvVtx=fTR->pho_Cone01PhotonIso_dEta015EB_dR070EE_mvVtx[index];
  pholead_pho_Cone02PhotonIso_dEta015EB_dR070EE_mvVtx=fTR->pho_Cone02PhotonIso_dEta015EB_dR070EE_mvVtx[index];
  pholead_pho_Cone03PhotonIso_dEta015EB_dR070EE_mvVtx=fTR->pho_Cone03PhotonIso_dEta015EB_dR070EE_mvVtx[index];
  //pholead_pho_Cone04PhotonIso_dEta015EB_dR070EE_mvVtx=fTR->pho_Cone04PhotonIso_dEta015EB_dR070EE_mvVtx[index];
  isolations_struct isos = PFConeIsolation(fTR,index);
  pholead_pho_Cone04PhotonIso_dEta015EB_dR070EE_mvVtx=isos.photon;
  //  std::cout << "debug PFPhotonIso " << PFIsolation(index,-999,"photon") << " " << pholead_pho_Cone04PhotonIso_dEta015EB_dR070EE_mvVtx << std::endl;
  pholead_pho_Cone01NeutralHadronIso_mvVtx=fTR->pho_Cone01NeutralHadronIso_mvVtx[index];
  pholead_pho_Cone02NeutralHadronIso_mvVtx=fTR->pho_Cone02NeutralHadronIso_mvVtx[index];
  pholead_pho_Cone03NeutralHadronIso_mvVtx=fTR->pho_Cone03NeutralHadronIso_mvVtx[index];
  //  pholead_pho_Cone04NeutralHadronIso_mvVtx=fTR->pho_Cone04NeutralHadronIso_mvVtx[index];
  pholead_pho_Cone04NeutralHadronIso_mvVtx=isos.neutral;
  pholead_pho_Cone01ChargedHadronIso_dR02_dz02_dxy01=fTR->pho_Cone01ChargedHadronIso_dR02_dz02_dxy01[index];
  pholead_pho_Cone02ChargedHadronIso_dR02_dz02_dxy01=fTR->pho_Cone02ChargedHadronIso_dR02_dz02_dxy01[index];
  pholead_pho_Cone03ChargedHadronIso_dR02_dz02_dxy01=fTR->pho_Cone03ChargedHadronIso_dR02_dz02_dxy01[index];
  //  pholead_pho_Cone04ChargedHadronIso_dR02_dz02_dxy01=fTR->pho_Cone04ChargedHadronIso_dR02_dz02_dxy01[index];
  pholead_pho_Cone04ChargedHadronIso_dR02_dz02_dxy01=isos.charged;
  pholead_pho_Cone03PFCombinedIso=fTR->pho_Cone03PFCombinedIso[index];
  //  pholead_pho_Cone04PFCombinedIso=fTR->pho_Cone04PFCombinedIso[index];
  pholead_pho_Cone04PFCombinedIso=pholead_pho_Cone04PhotonIso_dEta015EB_dR070EE_mvVtx+pholead_pho_Cone04NeutralHadronIso_mvVtx+pholead_pho_Cone04ChargedHadronIso_dR02_dz02_dxy01;
  pholead_Npfcandphotonincone = isos.nphotoncand;
  pholead_Npfcandchargedincone = isos.nchargedcand;
  pholead_Npfcandneutralincone = isos.nneutralcand;
  for (int i=0; i<pholead_Npfcandphotonincone && i<global_size_pfcandarrays; i++) pholead_photonpfcandenergies[i] = isos.photoncandenergies.at(i);
  for (int i=0; i<pholead_Npfcandchargedincone && i<global_size_pfcandarrays; i++) pholead_chargedpfcandenergies[i] = isos.chargedcandenergies.at(i);
  for (int i=0; i<pholead_Npfcandneutralincone && i<global_size_pfcandarrays; i++) pholead_neutralpfcandenergies[i] = isos.neutralcandenergies.at(i);
  for (int i=0; i<pholead_Npfcandphotonincone && i<global_size_pfcandarrays; i++) pholead_photonpfcandets[i] = isos.photoncandets.at(i);
  for (int i=0; i<pholead_Npfcandchargedincone && i<global_size_pfcandarrays; i++) pholead_chargedpfcandets[i] = isos.chargedcandets.at(i);
  for (int i=0; i<pholead_Npfcandneutralincone && i<global_size_pfcandarrays; i++) pholead_neutralpfcandets[i] = isos.neutralcandets.at(i);
  for (int i=0; i<pholead_Npfcandphotonincone && i<global_size_pfcandarrays; i++) pholead_photonpfcanddetas[i] = isos.photoncanddetas.at(i);
  for (int i=0; i<pholead_Npfcandchargedincone && i<global_size_pfcandarrays; i++) pholead_chargedpfcanddetas[i] = isos.chargedcanddetas.at(i);
  for (int i=0; i<pholead_Npfcandneutralincone && i<global_size_pfcandarrays; i++) pholead_neutralpfcanddetas[i] = isos.neutralcanddetas.at(i);
  for (int i=0; i<pholead_Npfcandphotonincone && i<global_size_pfcandarrays; i++) pholead_photonpfcanddphis[i] = isos.photoncanddphis.at(i);
  for (int i=0; i<pholead_Npfcandchargedincone && i<global_size_pfcandarrays; i++) pholead_chargedpfcanddphis[i] = isos.chargedcanddphis.at(i);
  for (int i=0; i<pholead_Npfcandneutralincone && i<global_size_pfcandarrays; i++) pholead_neutralpfcanddphis[i] = isos.neutralcanddphis.at(i);
  //  std::cout << "debug PFCombinedIso " << PFIsolation(index) << " " << pholead_pho_Cone04PFCombinedIso*pholead_pt << std::endl;
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
  int a;
  pholead_hasimpingingtrack = FindImpingingTrack(fTR,index,a);
  //  pholead_Nchargedhadronsincone = CountChargedHadronsInCone(fTR,index,remove,global_dofootprintremoval);
  pholead_scarea = scarea[fTR->PhotSCindex[index]];
  pholead_scareaSF = scareaSF[fTR->PhotSCindex[index]];
};

float DiPhotonMiniTree::SieieRescale(float sieie, bool isbarrel){
  if (isdata) return sieie; // rescale sieie only in MC
  return isbarrel ? 0.87*sieie+0.0011 : 0.99*sieie;
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
  photrail_sieie = SieieRescale(fTR->PhoSigmaIetaIeta[index],(bool)(fabs(fTR->SCEta[fTR->PhotSCindex[index]])<1.4442));
  photrail_hoe = fTR->PhoHoverE[index];
  photrail_brem = fTR->SCBrem[fTR->PhotSCindex[index]];
  photrail_sigmaPhi = fTR->SCPhiWidth[fTR->PhotSCindex[index]];
  photrail_sigmaEta = fTR->SCEtaWidth[fTR->PhotSCindex[index]];
  photrail_pho_Cone01PhotonIso_dEta015EB_dR070EE_mvVtx=fTR->pho_Cone01PhotonIso_dEta015EB_dR070EE_mvVtx[index];
  photrail_pho_Cone02PhotonIso_dEta015EB_dR070EE_mvVtx=fTR->pho_Cone02PhotonIso_dEta015EB_dR070EE_mvVtx[index];
  photrail_pho_Cone03PhotonIso_dEta015EB_dR070EE_mvVtx=fTR->pho_Cone03PhotonIso_dEta015EB_dR070EE_mvVtx[index];
  //  photrail_pho_Cone04PhotonIso_dEta015EB_dR070EE_mvVtx=fTR->pho_Cone04PhotonIso_dEta015EB_dR070EE_mvVtx[index];
  isolations_struct isos = PFConeIsolation(fTR,index);
  photrail_pho_Cone04PhotonIso_dEta015EB_dR070EE_mvVtx=isos.photon;
  photrail_pho_Cone01NeutralHadronIso_mvVtx=fTR->pho_Cone01NeutralHadronIso_mvVtx[index];
  photrail_pho_Cone02NeutralHadronIso_mvVtx=fTR->pho_Cone02NeutralHadronIso_mvVtx[index];
  photrail_pho_Cone03NeutralHadronIso_mvVtx=fTR->pho_Cone03NeutralHadronIso_mvVtx[index];
  //  photrail_pho_Cone04NeutralHadronIso_mvVtx=fTR->pho_Cone04NeutralHadronIso_mvVtx[index];
  photrail_pho_Cone04NeutralHadronIso_mvVtx=isos.neutral;
  photrail_pho_Cone01ChargedHadronIso_dR02_dz02_dxy01=fTR->pho_Cone01ChargedHadronIso_dR02_dz02_dxy01[index];
  photrail_pho_Cone02ChargedHadronIso_dR02_dz02_dxy01=fTR->pho_Cone02ChargedHadronIso_dR02_dz02_dxy01[index];
  photrail_pho_Cone03ChargedHadronIso_dR02_dz02_dxy01=fTR->pho_Cone03ChargedHadronIso_dR02_dz02_dxy01[index];
  //  photrail_pho_Cone04ChargedHadronIso_dR02_dz02_dxy01=fTR->pho_Cone04ChargedHadronIso_dR02_dz02_dxy01[index];
  photrail_pho_Cone04ChargedHadronIso_dR02_dz02_dxy01=isos.charged;
  photrail_pho_Cone03PFCombinedIso=fTR->pho_Cone03PFCombinedIso[index];
  //  photrail_pho_Cone04PFCombinedIso=fTR->pho_Cone04PFCombinedIso[index];
  photrail_pho_Cone04PFCombinedIso=photrail_pho_Cone04PhotonIso_dEta015EB_dR070EE_mvVtx+photrail_pho_Cone04NeutralHadronIso_mvVtx+photrail_pho_Cone04ChargedHadronIso_dR02_dz02_dxy01;
  photrail_Npfcandphotonincone = isos.nphotoncand;
  photrail_Npfcandchargedincone = isos.nchargedcand;
  photrail_Npfcandneutralincone = isos.nneutralcand;
  for (int i=0; i<photrail_Npfcandphotonincone && i<global_size_pfcandarrays; i++) photrail_photonpfcandenergies[i] = isos.photoncandenergies.at(i);
  for (int i=0; i<photrail_Npfcandchargedincone && i<global_size_pfcandarrays; i++) photrail_chargedpfcandenergies[i] = isos.chargedcandenergies.at(i);
  for (int i=0; i<photrail_Npfcandneutralincone && i<global_size_pfcandarrays; i++) photrail_neutralpfcandenergies[i] = isos.neutralcandenergies.at(i);
  for (int i=0; i<photrail_Npfcandphotonincone && i<global_size_pfcandarrays; i++) photrail_photonpfcandets[i] = isos.photoncandets.at(i);
  for (int i=0; i<photrail_Npfcandchargedincone && i<global_size_pfcandarrays; i++) photrail_chargedpfcandets[i] = isos.chargedcandets.at(i);
  for (int i=0; i<photrail_Npfcandneutralincone && i<global_size_pfcandarrays; i++) photrail_neutralpfcandets[i] = isos.neutralcandets.at(i);
  for (int i=0; i<photrail_Npfcandphotonincone && i<global_size_pfcandarrays; i++) photrail_photonpfcanddetas[i] = isos.photoncanddetas.at(i);
  for (int i=0; i<photrail_Npfcandchargedincone && i<global_size_pfcandarrays; i++) photrail_chargedpfcanddetas[i] = isos.chargedcanddetas.at(i);
  for (int i=0; i<photrail_Npfcandneutralincone && i<global_size_pfcandarrays; i++) photrail_neutralpfcanddetas[i] = isos.neutralcanddetas.at(i);
  for (int i=0; i<photrail_Npfcandphotonincone && i<global_size_pfcandarrays; i++) photrail_photonpfcanddphis[i] = isos.photoncanddphis.at(i);
  for (int i=0; i<photrail_Npfcandchargedincone && i<global_size_pfcandarrays; i++) photrail_chargedpfcanddphis[i] = isos.chargedcanddphis.at(i);
  for (int i=0; i<photrail_Npfcandneutralincone && i<global_size_pfcandarrays; i++) photrail_neutralpfcanddphis[i] = isos.neutralcanddphis.at(i);
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
  int a;
  photrail_hasimpingingtrack = FindImpingingTrack(fTR,index,a);
  //  photrail_Nchargedhadronsincone = CountChargedHadronsInCone(fTR,index,remove,global_dofootprintremoval);
  photrail_scarea = scarea[fTR->PhotSCindex[index]];
  photrail_scareaSF = scareaSF[fTR->PhotSCindex[index]];

};

void DiPhotonMiniTree::FillMuonInfo(int index){

  pholead_eta = fTR->MuEta[index];
  pholead_px = fTR->MuPx[index];
  pholead_py = fTR->MuPy[index];
  pholead_pt = fTR->MuPt[index];
  pholead_pz = fTR->MuPz[index];
  pholead_energy = fTR->MuE[index];
  pholead_SCeta = fTR->MuEta[index];
  pholead_SCphi = fTR->MuPhi[index];
  pholead_pho_Cone03PFCombinedIso = fTR->MuRelIso03[index];
  pholead_Npfcandphotonincone = 0;
  std::vector<float> energies;
  std::vector<float> ets;
  std::vector<float> detas;
  std::vector<float> dphis;
  pholead_pho_Cone04PhotonIso_dEta015EB_dR070EE_mvVtx = PFPhotonIsolationAroundMuon(index,&pholead_Npfcandphotonincone,&energies,&ets,&detas,&dphis);
  for (int i=0; i<pholead_Npfcandphotonincone && i<global_size_pfcandarrays; i++) pholead_photonpfcandenergies[i] = energies.at(i);
  for (int i=0; i<pholead_Npfcandphotonincone && i<global_size_pfcandarrays; i++) pholead_photonpfcandets[i] = ets.at(i);
  for (int i=0; i<pholead_Npfcandphotonincone && i<global_size_pfcandarrays; i++) pholead_photonpfcanddetas[i] = detas.at(i);
  for (int i=0; i<pholead_Npfcandphotonincone && i<global_size_pfcandarrays; i++) pholead_photonpfcanddphis[i] = dphis.at(i);
};

void DiPhotonMiniTree::ResetVars(){

  event_pass12whoissiglike = -999;
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
  pholead_Nchargedhadronsincone = -999;
  photrail_Nchargedhadronsincone = -999;
  pholead_scarea = -999;
  photrail_scarea = -999;
  pholead_scareaSF = -999;
  photrail_scareaSF = -999;
  pholead_Npfcandphotonincone = -999;
  pholead_Npfcandchargedincone = -999;
  pholead_Npfcandneutralincone = -999;
  photrail_Npfcandphotonincone = -999;
  photrail_Npfcandchargedincone = -999;
  photrail_Npfcandneutralincone = -999;
  for (int i=0; i<global_size_pfcandarrays; i++) pholead_photonpfcandenergies[i]=-999;
  for (int i=0; i<global_size_pfcandarrays; i++) pholead_chargedpfcandenergies[i]=-999;
  for (int i=0; i<global_size_pfcandarrays; i++) pholead_neutralpfcandenergies[i]=-999;
  for (int i=0; i<global_size_pfcandarrays; i++) pholead_photonpfcandets[i]=-999;
  for (int i=0; i<global_size_pfcandarrays; i++) pholead_chargedpfcandets[i]=-999;
  for (int i=0; i<global_size_pfcandarrays; i++) pholead_neutralpfcandets[i]=-999;
  for (int i=0; i<global_size_pfcandarrays; i++) pholead_photonpfcanddetas[i]=-999;
  for (int i=0; i<global_size_pfcandarrays; i++) pholead_chargedpfcanddetas[i]=-999;
  for (int i=0; i<global_size_pfcandarrays; i++) pholead_neutralpfcanddetas[i]=-999;
  for (int i=0; i<global_size_pfcandarrays; i++) pholead_photonpfcanddphis[i]=-999;
  for (int i=0; i<global_size_pfcandarrays; i++) pholead_chargedpfcanddphis[i]=-999;
  for (int i=0; i<global_size_pfcandarrays; i++) pholead_neutralpfcanddphis[i]=-999;
  for (int i=0; i<global_size_pfcandarrays; i++) photrail_photonpfcandenergies[i]=-999;
  for (int i=0; i<global_size_pfcandarrays; i++) photrail_chargedpfcandenergies[i]=-999;
  for (int i=0; i<global_size_pfcandarrays; i++) photrail_neutralpfcandenergies[i]=-999;
  for (int i=0; i<global_size_pfcandarrays; i++) photrail_photonpfcandets[i]=-999;
  for (int i=0; i<global_size_pfcandarrays; i++) photrail_chargedpfcandets[i]=-999;
  for (int i=0; i<global_size_pfcandarrays; i++) photrail_neutralpfcandets[i]=-999;
  for (int i=0; i<global_size_pfcandarrays; i++) photrail_photonpfcanddetas[i]=-999;
  for (int i=0; i<global_size_pfcandarrays; i++) photrail_chargedpfcanddetas[i]=-999;
  for (int i=0; i<global_size_pfcandarrays; i++) photrail_neutralpfcanddetas[i]=-999;
  for (int i=0; i<global_size_pfcandarrays; i++) photrail_photonpfcanddphis[i]=-999;
  for (int i=0; i<global_size_pfcandarrays; i++) photrail_chargedpfcanddphis[i]=-999;
  for (int i=0; i<global_size_pfcandarrays; i++) photrail_neutralpfcanddphis[i]=-999;
};

float DiPhotonMiniTree::CalculateSCArea(TreeReader *fTR, int scindex){
  //  std::cout << "call scarea " << scindex << " " << fTR->SCNXtals[scindex] << std::endl; // DEBUUUUG
  if (scindex>=fTR->NSuperClusters) return -999;
  float area=0;
  for (int i=0; i<fTR->SCNXtals[scindex]; i++) area+=fTR->SCxtalEtaWidth[scindex][i]*fTR->SCxtalPhiWidth[scindex][i];
  return area;
};

float DiPhotonMiniTree::GetPUEnergy(TreeReader *fTR, TString mode, float eta){
  return 0; /// XXX

  eta=fabs(eta);
  bool isbarrel = (eta<1.4442);

  int bin = Choose_bin_eta(fabs(eta),isbarrel ? 0 : 1);

  float eff_area = 0;

  std::vector<float> *eff_areas_EB_photon = isdata ?  &eff_areas_EB_photon_data :  &eff_areas_EB_photon_MC;
  std::vector<float> *eff_areas_EB_charged = isdata ? &eff_areas_EB_charged_data : &eff_areas_EB_charged_MC;
  std::vector<float> *eff_areas_EB_neutral = isdata ? &eff_areas_EB_neutral_data : &eff_areas_EB_neutral_MC;
  std::vector<float> *eff_areas_EE_photon = isdata ?  &eff_areas_EB_photon_data :  &eff_areas_EB_photon_MC;
  std::vector<float> *eff_areas_EE_charged = isdata ? &eff_areas_EE_charged_data : &eff_areas_EE_charged_MC;
  std::vector<float> *eff_areas_EE_neutral = isdata ? &eff_areas_EE_neutral_data : &eff_areas_EE_neutral_MC;

  if (mode=="photon") eff_area = isbarrel ?  (*eff_areas_EB_photon)[bin] :  (*eff_areas_EE_photon)[bin];
  if (mode=="charged") eff_area = isbarrel ? (*eff_areas_EB_charged)[bin] : (*eff_areas_EE_charged)[bin];
  if (mode=="neutral") eff_area = isbarrel ? (*eff_areas_EB_neutral)[bin] : (*eff_areas_EE_neutral)[bin];

  if (eff_area==0) std::cout << "Warning: problem in eta-dependent EA" << std::endl;

  return TMath::Pi()*0.4*0.4*eff_area*fTR->Rho;

};

Int_t DiPhotonMiniTree::Choose_bin_eta(float eta, int region){

  eta=fabs(eta);

  int index;

  std::vector<float> *cuts;

  if (region==0) {cuts=&binsdef_single_gamma_EB_eta; index=n_templates_EB;}
  if (region==1) {cuts=&binsdef_single_gamma_EE_eta; index=n_templates_EE;}

  assert (index!=0);

  (*cuts)[index]=9999;

  if (eta<(*cuts)[0]){
    std::cout << "WARNING: called bin choice for out-of-range value " << eta << " cuts[0]=" << (*cuts)[0] << std::endl;
    return -999;
  }

  for (int i=0; i<index; i++) if ((eta>=(*cuts)[i]) && (eta<(*cuts)[i+1])) return i;

  std::cout << "WARNING: called bin choice for out-of-range value " << eta << std::endl;
  return -999;

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

std::vector<int> DiPhotonMiniTree::ImpingingTrackSelection(TreeReader *fTR, std::vector<int> passing, bool invert){

  return std::vector<int>();

//
//  for (int i=0; i<100; i++) impinging_track_pfcand[i]=-999;
//
//  for (vector<int>::iterator it = passing.begin(); it != passing.end(); ){
//    
//    bool found=0;
//    int phoqi=*it;
//    
//    found = FindImpingingTrack(fTR,phoqi,impinging_track_pfcand[phoqi]);
//
//    if (invert) { // selection = 0 impinging tracks
//      if (found) it=passing.erase(it); else it++;
//    }
//
//    else { // selection = 1 impinging track
//
//      if (found){    
//	std::vector<int> remove;
//	remove.push_back(impinging_track_pfcand[phoqi]);
//	float eta=fTR->SCEta[fTR->PhotSCindex[*it]];
//	const float eff_area = (fabs(eta)<1.4442) ? eff_area_EB : eff_area_EE;
//	const float dR=0.4;
//	float puenergy =3.14*dR*dR*eff_area*fTR->Rho;
//	if (PFIsolation(phoqi,0,"combined",remove)-puenergy>5) found=0;
//      }
//      
//      if (!found) it=passing.erase(it); else it++;
//    
//    }
//
//  }
//
//  return passing;
//  
};


bool DiPhotonMiniTree::FindImpingingTrack(TreeReader *fTR, int phoqi, int &reference_index_found, bool dofootprintremoval, std::vector<int> removals){

  return false;

  // WARNING THE FOLLOWING IS NOT UP TO DATE!!!
//  bool found = false;
//
//  TVector3 photon_position = TVector3(fTR->SCx[fTR->PhotSCindex[phoqi]],fTR->SCy[fTR->PhotSCindex[phoqi]],fTR->SCz[fTR->PhotSCindex[phoqi]]);
//
//  if (dofootprintremoval){
//    std::vector<int> footprint = GetPFCandInsideFootprint(fTR,phoqi,0,"charged");
//    for (int i=0; i<footprint.size(); i++) removals.push_back(footprint.at(i));
//  }
//
//  for (int i=0; i<fTR->NPfCand; i++){
//
//    if (fTR->Pho_isPFPhoton[phoqi] && fTR->pho_matchedPFPhotonCand[phoqi]==i) continue;
//
//    int type = FindPFCandType(fTR->PfCandPdgId[i]);
//    if (type!=1) continue;
//
//    bool removed = false;
//    for (int j=0; j<removals.size(); j++) {
//      if (i==removals.at(j)) removed=true;
//    }
//    if (removed) continue;
//
//    double pt = fTR->PfCandPt[i];
//    if (pt<1.5) continue;
//
//    if (!(fTR->PfCandHasHitInFirstPixelLayer[i])) continue;
//
//    if (dofootprintremoval){
//      
//      std::cout << "Wrong! No propagator for charged stuff" << std::endl;
//      continue;
//
//      //      double dR = GetPFCandDeltaRFromSC(fTR,phoqi,i,0);
//      //      if (dR>0.4) continue;
//    }
//
//    if (!dofootprintremoval){ // veto cones
//
//      TVector3 phovtx(fTR->PhoVx[phoqi],fTR->PhoVy[phoqi],fTR->PhoVz[phoqi]);
//      TVector3 pfvertex(fTR->PfCandVx[i],fTR->PfCandVy[i],fTR->PfCandVz[i]);
//      TVector3 photon_direction = photon_position-pfvertex;
//      float sceta = photon_direction.Eta();
//      float scphi = photon_direction.Phi();
//      double dxy;
//      double dz;
//	
//      TVector3 vtxmom(fTR->PfCandTrackRefPx[i],fTR->PfCandTrackRefPy[i],fTR->PfCandTrackRefPz[i]);
//
//      if (vtxmom.x()==-999 || vtxmom.y()==-999 || vtxmom.z()==-999) {
//	std::cout << "Something wrong with vtxmom from trackref, fallback" << std::endl;
//	vtxmom = TVector3(fTR->PfCandPx[i],fTR->PfCandPy[i],fTR->PfCandPz[i]);
//      }
//      
//      dxy = ( -(pfvertex.x()-phovtx.x())*vtxmom.y() +(pfvertex.y()-phovtx.y())*vtxmom.x() ) / vtxmom.Perp();
//      dz = (pfvertex.z()-phovtx.z()) - ( (pfvertex.x()-phovtx.x())*vtxmom.x() + (pfvertex.y()-phovtx.y())*vtxmom.y() ) / vtxmom.Perp() * vtxmom.z() / vtxmom.Perp();
//      dxy=fabs(dxy);
//      dz=fabs(dz);
//
//
//      double dEta = fTR->PfCandEta[i] - sceta;
//      double dPhi = Util::DeltaPhi(fTR->PfCandPhi[i],scphi);
//      double dR = sqrt(dEta*dEta+dPhi*dPhi);
//
//      if (dR>0.4) continue;
//
//      // additional veto sieie-orthogonal
//      if (fTR->PhoisEB[phoqi]){
//	if (fabs(dEta)<0.05) continue;
//	if (fabs(dPhi)<0.05) continue;
//      }
//
//      if (dz>0.2) continue;
//      if (dxy>0.1) continue;
//      if (dR<0.02) continue;
//
//    }
//
//    found=1;
//    reference_index_found = i;
//    break;
//
//  } // end pf cand loop
//
//  return found;
//
};

int DiPhotonMiniTree::CountChargedHadronsInCone(TreeReader *fTR, int phoqi, std::vector<int> removals, bool skipvetocones){

  return 0;

//  int found = 0;
//  
//  TVector3 photon_position = TVector3(fTR->SCx[fTR->PhotSCindex[phoqi]],fTR->SCy[fTR->PhotSCindex[phoqi]],fTR->SCz[fTR->PhotSCindex[phoqi]]);
//  TVector3 phovtx(fTR->PhoVx[phoqi],fTR->PhoVy[phoqi],fTR->PhoVz[phoqi]);
//    
//  for (int i=0; i<fTR->NPfCand; i++){
//    
//    bool removed = false;
//    for (int j=0; j<removals.size(); j++) {
//      if (i==removals.at(j)) removed=true;
//    }
//    if (removed) continue;
//    
//    float id = fTR->PfCandPdgId[i];
//    if (!(fabs(id)==211 || fabs(id)==321 || id==999211 || fabs(id)==2212)) continue;
//    
//    TVector3 pfvertex(fTR->PfCandVx[i],fTR->PfCandVy[i],fTR->PfCandVz[i]);
//    TVector3 photon_direction = photon_position-pfvertex;
//      
//    double sceta = photon_direction.Eta();
//    double scphi = photon_direction.Phi();
//
//    double dEta = fTR->PfCandEta[i] - sceta;
//    double dPhi = Util::DeltaPhi(fTR->PfCandPhi[i],scphi);
//    double dR = sqrt(dEta*dEta+dPhi*dPhi);
//
//    double pt = fTR->PfCandPt[i];
//
//    TVector3 vtxmom(fTR->PfCandTrackRefPx[i],fTR->PfCandTrackRefPy[i],fTR->PfCandTrackRefPz[i]);
//
//    if (vtxmom.x()==-999 || vtxmom.y()==-999 || vtxmom.z()==-999) {
//      std::cout << "Something wrong with vtxmom from trackref, fallback" << std::endl;
//      vtxmom = TVector3(fTR->PfCandPx[i],fTR->PfCandPy[i],fTR->PfCandPz[i]);
//    }
//      
//    double dxy = ( -(pfvertex.x()-phovtx.x())*vtxmom.y() +(pfvertex.y()-phovtx.y())*vtxmom.x() ) / vtxmom.Perp();
//    double dz = (pfvertex.z()-phovtx.z()) - ( (pfvertex.x()-phovtx.x())*vtxmom.x() + (pfvertex.y()-phovtx.y())*vtxmom.y() ) / vtxmom.Perp() * vtxmom.z() / vtxmom.Perp();
//    dxy=fabs(dxy);
//    dz=fabs(dz);
//
//    if (dR>0.4) continue;
//
//    if (!skipvetocones){
//      // additional veto sieie-orthogonal
//      if (fabs(sceta)<1.4442){
//	if (fabs(dEta)<0.05) continue;
//	if (fabs(dPhi)<0.05) continue;
//      }
//      
//      if (dz>0.2) continue;
//      if (dxy>0.1) continue;
//      if (dR<0.02) continue;
//    }
//
//    if (!(fTR->PfCandHasHitInFirstPixelLayer[i])) continue;
//
//    found++;
//      
//  } // end pf cand loop
//
//  return found;
//
};

std::vector<int> DiPhotonMiniTree::NChargedHadronsInConeSelection(TreeReader *fTR, std::vector<int> passing, int minimum, int maximum){

  return std::vector<int>();

//  if (maximum<minimum) {
//    std::cout << "Wrong values for NChargedHadronsInCone selection!!!" << std::endl;
//    return std::vector<int>();
//  }
//
//  for (vector<int>::iterator it = passing.begin(); it != passing.end(); ){
//    bool pass=0;
//    std::vector<int> remove = GetPFCandRemovals(fTR,*it);
//    std::vector<int> footprint = GetPFCandInFootprint(fTR,*it,-999);
//    for (int i=0; i<footprint.size(); i++) remove.push_back(footprint.at(i));
//    int n = CountChargedHadronsInCone(fTR,*it,remove,trueblabla);
//    if (n>=minimum && n<=maximum) pass=1; else pass=0;
//    if (!pass) it=passing.erase(it); else it++;
//  }
//
//  return passing;
//
};

void DiPhotonMiniTree::Fillhist_PFPhotonDepositAroundImpingingTrack(int phoqi, int trkindex){
  return;
};

float DiPhotonMiniTree::DeltaPhiSigned(float phi1, float phi2){
  //copy-paste from reco::deltaPhi
  double result = phi1 - phi2;
  const float pi = TMath::Pi();
  while (result > pi) result -= 2*pi;
  while (result <= -pi) result += 2*pi;
  return result;
};

bool DiPhotonMiniTree::PassPrimaryVertexFilter(){
  for (int i=0; i<fTR->NVrtx; i++){
    if (fTR->VrtxNdof[i]>4 && fabs(fTR->VrtxZ[i])<=24 && TVector3(fTR->VrtxX[i],fTR->VrtxY[i],fTR->VrtxZ[i]).Perp()<=2) return true;
  }
  return false;
};
