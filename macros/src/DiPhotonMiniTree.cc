#include "helper/Utilities.hh"
#include "DiPhotonMiniTree.hh"

#include "DiPhotonPurity.hh"

#include <iostream>
#include <fstream>
#include <assert.h>

using namespace std;

DiPhotonMiniTree::DiPhotonMiniTree(TreeReader *tr, std::string dataType) : UserAnalysisBase(tr), fDataType_(dataType){
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

	fHNumPU = new TH1F("NumPU","NumPU",40,0,40);
	fHNumVtx = new TH1F("NumVtx","NumVtx",40,0,40);
	

	
	  cout << "Tree and histos created" << endl;
	

}

void DiPhotonMiniTree::Analyze(){

  //cout << "Analyze this event" << endl;


  //cout << "A" << endl;

  float weight;
  if (!isdata) weight = GetPUWeight(fTR->PUnumInteractions);
  else weight=1;

  if (!isdata) fHNumPU->Fill(fTR->PUnumInteractions,weight);
  fHNumVtx->Fill(fTR->NVrtx,weight);

  if (isdata && !TriggerSelection()) return;

  std::vector<int> passing = PhotonSelection(fTR);

   bool evtisok=true;

   for (vector<int>::iterator it = passing.begin(); it != passing.end(); it++){
     float eta=fTR->SCEta[fTR->PhotSCindex[*it]];
     float phi=fTR->SCPhi[fTR->PhotSCindex[*it]];
     if ( (fabs(eta)>1.4442 && fabs(eta)<1.56) || (fabs(eta)>2.5) || (phocorr->isInPhiCracks(phi,eta)) || (phocorr->isInEBEtaCracks(eta)) ) evtisok=false;
   }

   if (!evtisok) return;

   if (passing.size()<2) return;

   // TLorentzVector pho[2];
   // for (int i=0;i<2;i++){
   //   pho[i].SetPtEtaPhiE(fTR->PhoPt[passing.at(i)],fTR->PhoEta[passing.at(i)],fTR->PhoPhi[passing.at(i)],fTR->PhoEnergy[passing.at(i)]);
   // }

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
     float energy=fTR->SCRaw[fTR->PhotSCindex[*it]];
     float eta=fTR->SCEta[fTR->PhotSCindex[*it]];
     if (fabs(eta)<1.4442) energy*=phocorr->getEtaCorrectionBarrel(eta);
     if (fabs(eta)>1.56) energy+=fTR->SCPre[fTR->PhotSCindex[*it]];
     if (energy/cosh(eta)<10) it=passing.erase(it); else it++;
   }

     for (vector<int>::iterator it = passing.begin(); it != passing.end(); ){
       if (!PhotonID_EGM_10_006_Loose(fTR,*it)) it=passing.erase(it); else it++;
       //    if (!PhotonID_EGM_10_006_Loose_SigmaIetaIeta_Relaxed(fTR,*it)) it=passing.erase(it); else it++;
     }

     return passing;

};




bool DiPhotonMiniTree::PhotonID_EGM_10_006_Loose(TreeReader *fTR, int i){

  if (fTR->PhoIso04Ecal[i]>4.2) return false;
  if (fTR->PhoIso04Hcal[i]>2.2) return false;
  if (fTR->PhoIso04TrkHollow[i]>2.0) return false;
  if (fTR->PhoHoverE[i]>0.05) return false;


  if (fabs(fTR->PhoEta[i])<1.4442) {
    if (fTR->PhoSigmaIetaIeta[i]>0.01) return false;
  }
  else {
    if (fTR->PhoSigmaIetaIeta[i]>0.03) return false;
  }


  return true;

};

bool DiPhotonMiniTree::PhotonID_EGM_10_006_Loose_SigmaIetaIeta_Relaxed(TreeReader *fTR, int i){

  if (fTR->PhoIso04Ecal[i]>4.2) return false;
  if (fTR->PhoIso04Hcal[i]>2.2) return false;
  if (fTR->PhoIso04TrkHollow[i]>2.0) return false;
  if (fTR->PhoHoverE[i]>0.05) return false;

  /*
  if (fabs(fTR->PhoEta[i])<1.4442) {
    if (fTR->PhoSigmaIetaIeta[i]>0.01) return false;
  }
  else {
    if (fTR->PhoSigmaIetaIeta[i]>0.03) return false;
  }
  */

  return true;

};

bool DiPhotonMiniTree::TriggerSelection(){

  vector<string> triggers;

  triggers.push_back("HLT_Photon26_IsoVL_Photon18_v2");
  triggers.push_back("HLT_Photon20_R9Id_Photon18_R9Id_v2");
  triggers.push_back("HLT_Photon26_Photon18_v2");
  triggers.push_back("HLT_Photon26_IsoVL_Photon18_v2");
  triggers.push_back("HLT_Photon26_IsoVL_Photon18_IsoVL_v2");
  triggers.push_back("HLT_Photon26_CaloIdL_IsoVL_Photon18_v2");
  triggers.push_back("HLT_Photon26_CaloIdL_IsoVL_Photon18_R9Id_v1");
  triggers.push_back("HLT_Photon26_CaloIdL_IsoVL_Photon18_CaloIdL_IsoVL_v2");
  triggers.push_back("HLT_Photon26_R9Id_Photon18_CaloIdL_IsoVL_v1");
  triggers.push_back("HLT_Photon32_CaloIdL_Photon26_CaloIdL_v2");
  triggers.push_back("HLT_Photon36_CaloIdL_Photon22_CaloIdL_v1");

  for (vector<string>::const_iterator it=triggers.begin(); it!=triggers.end(); it++){
    if ( GetHLTPrescale(*it)!=0) {
      cout << "warning: using prescaled trigger!!! " << *it << " " << GetHLTPrescale(*it) << endl;
    }
    if ( GetHLTResult(*it) )        return true;
  }

 return false;

};
