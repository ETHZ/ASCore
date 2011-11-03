#include "helper/Utilities.hh"
#include "ZeeAnalysis.hh"

#include <iostream>
#include <fstream>
#include <assert.h>

using namespace std;

ZeeAnalysis::ZeeAnalysis(TreeReader *tr, std::string dataType) : UserAnalysisBase(tr), fDataType_(dataType){
	Util::SetStyle();	
	if (fDataType_ == "mc") isdata=false;
	else if (fDataType_ == "data") isdata=true; 
	else {
	  std::cout << "wrong data type" << std::endl;
	  assert(1==0);
	}
        elecorr = new EnergyCorrection("electrons");
	for (int i=0; i<6; i++) for (int j=0; j<4; j++) myfile[i][j] = new ofstream();
}

ZeeAnalysis::~ZeeAnalysis(){
  delete elecorr;
  for (int i=0; i<6; i++) for (int j=0; j<4; j++) delete myfile[i][j];
}

void ZeeAnalysis::Begin(){
	// Define the output file of histograms
	const char* filename = "histos_ee.root";
	fHistFile = new TFile(fOutputDir + TString(filename), "RECREATE");

	
	const char* masspointsfile = "masspoints_ee_";
	int codes[6]={0,15,16,17,20,999};
	for (int i=0; i<6; i++)	for (int j=0; j<4; j++){
	TString path=fOutputDir;
	path.Append(masspointsfile);
	path+=codes[i];
	path.Append("_cat");
	path+=j;
	path.Append(".txt");
      	myfile[i][j]->open(path.Data());
	}
	
	// Define the histograms

	fHElPt   = new TH1D("ElPt", "Pt of electrons", 100, 0., 500.);
	fHElPtCorr   = new TH1D("ElPtCorr", "Corrected Pt of electrons", 100, 0., 500.);
	fHInvMass0 = new TH1D("eeInvMass0","Invariant ee mass", 100, 70,110);
	fHInvMass15 = new TH1D("eeInvMass15","Invariant ee mass corr15", 100, 70,110);
	fHInvMass16 = new TH1D("eeInvMass16","Invariant ee mass corr16", 100, 70,110);
	fHInvMass17 = new TH1D("eeInvMass17","Invariant ee mass corr17", 100, 70,110);
	fHInvMass20 = new TH1D("eeInvMass20","Invariant ee mass corr20", 100, 70,110);
	fHInvMassEgen = new TH1D("eeInvMassEgen","Invariant ee mass corrEgen", 100, 70,110);

	fHErecEGen17  = new TH1D("ErecEgen17","ErecEgen17",100,0.5,1.5);
	fHErecEGen20  = new TH1D("ErecEgen20","ErecEgen20",100,0.5,1.5);


}

void ZeeAnalysis::Analyze(){

  /*
  for (int i=0; i<fTR->NEles; i++){
    fHElPt -> Fill(fTR->ElPt[i]);
    float correnergy = elecorr->get_correctedenergy(fTR,i,15);
    fHElPtCorr -> Fill(fTR->ElPt[i]*correnergy/fTR->ElE[i]);
   }
  */

  float weight;
  if (!isdata) weight = GetPUWeight(fTR->PUnumInteractions);
  else weight=1;


  vector<int> passing;
   for (int i=0; i<fTR->NEles; i++){
      passing.push_back(i);
  }

   for (vector<int>::iterator it = passing.begin(); it != passing.end(); ){
     if (IsGoodElId_WP80(*it)) it=passing.erase(it); else it++;
   }

   for (vector<int>::iterator it = passing.begin(); it != passing.end(); ){
     if (fTR->ElSCindex[*it]==-1) it=passing.erase(it); else it++;
   }

   for (vector<int>::iterator it = passing.begin(); it != passing.end(); ){
     float energy=fTR->SCRaw[fTR->ElSCindex[*it]];
     float eta=fTR->SCEta[fTR->ElSCindex[*it]];
     if (fabs(eta)<1.4442) energy*=elecorr->getEtaCorrectionBarrel(eta);
     if (fabs(eta)>1.56) energy+=fTR->SCPre[fTR->ElSCindex[*it]];
     if (energy/cosh(eta)<30) it=passing.erase(it); else it++;
   }

   for (vector<int>::iterator it = passing.begin(); it != passing.end(); ){
     float eta=fTR->SCEta[fTR->ElSCindex[*it]];
     float phi=fTR->SCPhi[fTR->ElSCindex[*it]];
     if ( (fabs(eta)>2.5) ) it=passing.erase(it); else it++;
     //     if ( (fabs(eta)>1.4442 && fabs(eta)<1.56) || (fabs(eta)>2.5) ) it=passing.erase(it); else it++;
     //     if ( (fabs(eta)>1.4442 && fabs(eta)<1.56) || (fabs(eta)>2.5) || (elecorr->isInPhiCracks(phi,eta))) it=passing.erase(it); else it++;
   }

   if (passing.size()<2) return;
   
   if (fTR->ElSCindex[passing.at(0)]==-1 || fTR->ElSCindex[passing.at(1)]==-1) return;

   TLorentzVector elec[2];
   for (int i=0; i<2; i++){
     elec[i].SetPtEtaPhiE(fTR->ElPt[passing.at(i)],fTR->ElEta[passing.at(i)],fTR->ElPhi[passing.at(i)],fTR->ElE[passing.at(i)]);
   }


  bool masswindow=false;
  if (fabs((elec[0]+elec[1]).M()-91.2)<30) masswindow=true;
  if (!masswindow) return;



  float invmass0=(elec[0]+elec[1]).M();
  float invmass15=(CorrElectron(fTR,passing.at(0),15)+CorrElectron(fTR,passing.at(1),15)).M();
  float invmass16=(CorrElectron(fTR,passing.at(0),16)+CorrElectron(fTR,passing.at(1),16)).M();
  float invmass17=(CorrElectron(fTR,passing.at(0),17)+CorrElectron(fTR,passing.at(1),17)).M();
  float invmass20=(CorrElectron(fTR,passing.at(0),20)+CorrElectron(fTR,passing.at(1),20)).M();
  
  TLorentzVector genelec[2];
  for (int i=0; i<2; i++) genelec[i].SetPtEtaPhiE(fTR->ElGenPt[passing.at(i)],fTR->ElGenEta[passing.at(i)],fTR->ElGenPhi[passing.at(i)],fTR->ElGenE[passing.at(i)]);
  
  float invmassEgen=(genelec[0]+genelec[1]).M();

  fHInvMass0->Fill(invmass0,weight);
  fHInvMass15 ->Fill(invmass15,weight);
  fHInvMass16 ->Fill(invmass16,weight);
  fHInvMass17 ->Fill(invmass17,weight);
  fHInvMass20 ->Fill(invmass20,weight);
  fHInvMassEgen->Fill(invmassEgen,weight);



   int cat;
   
   float abseta0 = fTR->SCEta[fTR->ElSCindex[passing.at(0)]];
   float abseta1 = fTR->SCEta[fTR->ElSCindex[passing.at(1)]];
   if (abseta0<1.4442 && abseta1<1.4442) cat=1;
   else if (abseta0>1.56 && abseta1>1.56) cat=3;
   else if (abseta0<1.4442 && abseta1>1.56) cat=2;
   else if (abseta0>1.56 && abseta1<1.4442) cat=2;
   else cat=-1;
   
   
   if (cat!=-1){
     *(myfile[0][0]) << invmass0 << " " << weight << std::endl;
     *(myfile[1][0]) << invmass15 << " " << weight << std::endl;
     *(myfile[2][0]) << invmass16 << " " << weight << std::endl;
     *(myfile[3][0]) << invmass17 << " " << weight << std::endl;
     *(myfile[4][0]) << invmass20 << " " << weight << std::endl;
     *(myfile[5][0]) << invmassEgen << " " << weight << std::endl;
     
     *(myfile[0][cat]) << invmass0 << " " << weight << std::endl;
     *(myfile[1][cat]) << invmass15 << " " << weight << std::endl;
     *(myfile[2][cat]) << invmass16 << " " << weight << std::endl;
     *(myfile[3][cat]) << invmass17 << " " << weight << std::endl;
     *(myfile[4][cat]) << invmass20 << " " << weight << std::endl;
     *(myfile[5][cat]) << invmassEgen << " " << weight << std::endl;
   }
 

   if (cat==1) for (int i=0; i<2; i++) {
     if (fTR->ElGenE[passing.at(i)]<0) continue;
     TLorentzVector correl17 = CorrElectron(fTR,passing.at(i),17);
     TLorentzVector correl20 = CorrElectron(fTR,passing.at(i),20);
     fHErecEGen17->Fill(correl17.E()/fTR->ElGenE[passing.at(i)]);
     fHErecEGen20->Fill(correl20.E()/fTR->ElGenE[passing.at(i)]);
   }
   
 
 
 

}

void ZeeAnalysis::End(){
	fHistFile->cd();	

	fHElPt   ->Write();
	fHElPtCorr -> Write();

	fHInvMass0 -> Write();
	fHInvMass15 -> Write();
	fHInvMass16 -> Write();
	fHInvMass17 -> Write();
	fHInvMass20 -> Write();
	fHInvMassEgen -> Write();

	fHErecEGen17->Write();
	fHErecEGen20->Write();


	fHistFile->Close();
	for (int i=0; i<5; i++) for (int j=0; j<4; j++)  myfile[i][j]->close();
}

TLorentzVector ZeeAnalysis::CorrElectron(TreeReader *fTR, int i, int mode){
  TLorentzVector corr;
  corr.SetPtEtaPhiE(fTR->ElPt[i],fTR->ElEta[i],fTR->ElPhi[i],fTR->ElE[i]);
  corr.SetE(elecorr->get_correctedenergy(fTR,i,mode));
  corr.SetRho(elecorr->get_correctedenergy(fTR,i,mode));
  return corr;
}
