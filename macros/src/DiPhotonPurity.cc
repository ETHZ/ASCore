#include "helper/Utilities.hh"
#include "DiPhotonPurity.hh"


#include <iostream>
#include <fstream>
#include <assert.h>
#include <vector>

using namespace std;

DiPhotonPurity::DiPhotonPurity(TreeReader *tr, std::string dataType, double aw) : UserAnalysisBase(tr), fDataType_(dataType), AddWeight(aw){
	Util::SetStyle();
	
	if (fDataType_ == "mc") isdata=false;
	else if (fDataType_ == "data") isdata=true; 
	else {
	  std::cout << "wrong data type" << std::endl;
	  assert(1==0);
	}
       phocorr = new EnergyCorrection("photons");

}

DiPhotonPurity::~DiPhotonPurity(){
  delete phocorr;
}

void DiPhotonPurity::Begin(){
	// Define the output file of histograms
	const char* filename = "diphoton_purity.root";
	fHistFile = new TFile(fOutputDir + TString(filename), "RECREATE");

	/*
	fHsieie_all = new TH1F("sieie_all","sieie_all",100,0,0.05);
	fHsieie_signal = new TH1F("sieie_signal","sieie_signal",100,0,0.05);
	fHsieie_background = new TH1F("sieie_background","sieie_background",100,0,0.05);
	*/
	fHgginvmass = new TH1F("gginvmass","gginvmass",110,80,300);

	fHNumPU = new TH1F("NumPU","NumPU",40,0,40);
	fHNumVtx = new TH1F("NumVtx","NumVtx",40,0,40);

}

void DiPhotonPurity::Analyze(){

  

  float weight;
  if (!isdata) weight = GetPUWeight(fTR->PUnumInteractions);
  else weight=1;
  weight*=AddWeight;

  if (!isdata && weight<0) weight=1; // TO ALLOW NOT USING PU REW IN MC

  if (!isdata) fHNumPU->Fill(fTR->PUnumInteractions,weight);
  fHNumVtx->Fill(fTR->NVrtx,weight);

   if (!TriggerSelection()) return;

  std::vector<int> passing = PhotonSelection(fTR);

   bool evtisok=true;

   for (vector<int>::iterator it = passing.begin(); it != passing.end(); it++){
     float eta=fTR->SCEta[fTR->PhotSCindex[*it]];
     float phi=fTR->SCPhi[fTR->PhotSCindex[*it]];
     if ( (fabs(eta)>1.4442 && fabs(eta)<1.56) || (fabs(eta)>2.5) || (phocorr->isInPhiCracks(phi,eta)) || (phocorr->isInEBEtaCracks(eta)) ) evtisok=false;
   }

   if (!evtisok) return;



   /*   if (!isdata){
     for (vector<int>::iterator it = passing.begin(); it != passing.end(); ){
       if (fTR->PhoMCmatchexitcode[*it]<0)  it=passing.erase(it); else it++;
     }
   }


   for (vector<int>::iterator it = passing.begin(); it != passing.end(); it++){
     fHsieie_all->Fill(fTR->PhoSigmaIetaIeta[*it],weight);
     if (!isdata){
       if (fTR->PhoMCmatchexitcode[*it]>0) fHsieie_signal->Fill(fTR->PhoSigmaIetaIeta[*it],weight);
       else fHsieie_background->Fill(fTR->PhoSigmaIetaIeta[*it],weight);
     }
   }

   */

   if (passing.size()<2) return;

   TLorentzVector pho[2];
   for (int i=0;i<2;i++){
     pho[i].SetPtEtaPhiE(fTR->PhoPt[passing.at(i)],fTR->PhoEta[passing.at(i)],fTR->PhoPhi[passing.at(i)],fTR->PhoEnergy[passing.at(i)]);
   }

   fHgginvmass->Fill((pho[0]+pho[1]).M(),weight);

}

void DiPhotonPurity::End(){



	fHistFile->cd();	
	/*
	fHsieie_all->Write();
	fHsieie_signal->Write();
	fHsieie_background->Write();
	*/
	fHNumPU->Write();
	fHNumVtx->Write();

	fHgginvmass->Write();

	
	fHistFile->Close();



}


TLorentzVector DiPhotonPurity::CorrPhoton(TreeReader *fTR, int i, int mode){
  TLorentzVector corr;
  corr.SetPtEtaPhiE(fTR->PhoPt[i],fTR->PhoEta[i],fTR->PhoPhi[i],fTR->PhoEnergy[i]);
  float corrE=phocorr->get_correctedenergy(fTR,i,mode);
  corr.SetE(corrE);
  corr.SetRho(corrE);
  return corr;
};


std::vector<int> DiPhotonPurity::PhotonSelection(TreeReader *fTR){

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




bool DiPhotonPurity::PhotonID_EGM_10_006_Loose(TreeReader *fTR, int i){

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

bool DiPhotonPurity::PhotonID_EGM_10_006_Loose_SigmaIetaIeta_Relaxed(TreeReader *fTR, int i){

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

bool DiPhotonPurity::TriggerSelection(){
#include "DiPhotonTriggerSelection.cc"
};
