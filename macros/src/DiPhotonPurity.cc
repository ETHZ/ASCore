#include "helper/Utilities.hh"
#include "DiPhotonPurity.hh"

#include <iostream>
#include <fstream>
#include <assert.h>
#include <vector>

using namespace std;

DiPhotonPurity::DiPhotonPurity(TreeReader *tr, std::string dataType) : UserAnalysisBase(tr), fDataType_(dataType){
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

	fHsieie = new TH1F("sieie","sieie",100,0,0.05);

	fHNumPU = new TH1F("NumPU","NumPU",40,0,40);
	fHNumVtx = new TH1F("NumVtx","NumVtx",40,0,40);

}

void DiPhotonPurity::Analyze(){

  

  float weight;
  if (!isdata) weight = GetPUWeight(fTR->PUnumInteractions);
  else weight=1;

  if (!isdata) fHNumPU->Fill(fTR->PUnumInteractions,weight);
  fHNumVtx->Fill(fTR->NVrtx,weight);



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
     if (energy/cosh(eta)<10 || energy/cosh(eta)>200) it=passing.erase(it); else it++;
   }



     for (vector<int>::iterator it = passing.begin(); it != passing.end(); ){
       if (!PhotonID_EGM_10_006_Loose(fTR,*it)) it=passing.erase(it); else it++;
     }



   bool evtisok=true;

   for (vector<int>::iterator it = passing.begin(); it != passing.end(); it++){
     float eta=fTR->SCEta[fTR->PhotSCindex[*it]];
     float phi=fTR->SCPhi[fTR->PhotSCindex[*it]];
     if ( (fabs(eta)>1.4442 && fabs(eta)<1.56) || (fabs(eta)>2.5) || (phocorr->isInPhiCracks(phi,eta)) || (phocorr->isInEBEtaCracks(eta)) ) evtisok=false;
   }

   if (!evtisok) return;



   for (vector<int>::iterator it = passing.begin(); it != passing.end(); it++){
     fHsieie->Fill(fTR->PhoSigmaIetaIeta[*it],weight);
   }
 
 

}

void DiPhotonPurity::End(){



	fHistFile->cd();	

	fHsieie->Write();

	fHNumPU->Write();
	fHNumVtx->Write();

	
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
