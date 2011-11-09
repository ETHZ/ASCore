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

	fHErecEGen17cat1  = new TH1D("ErecEgen17cat1","ErecEgen17cat1",100,0.5,1.5);
	fHErecEGen20cat1  = new TH1D("ErecEgen20cat1","ErecEgen20cat1",100,0.5,1.5);
	fHErecEGen17cat2  = new TH1D("ErecEgen17cat2","ErecEgen17cat2",100,0.5,1.5);
	fHErecEGen20cat2  = new TH1D("ErecEgen20cat2","ErecEgen20cat2",100,0.5,1.5);
	fHErecEGen17cat3  = new TH1D("ErecEgen17cat3","ErecEgen17cat3",100,0.5,1.5);
	fHErecEGen20cat3  = new TH1D("ErecEgen20cat3","ErecEgen20cat3",100,0.5,1.5);

	fHNumPU = new TH1I("NumPU","NumPU",40,0,40);

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

  fHNumPU->Fill(fTR->NVrtx,weight);


  bool evtisok = true;

  vector<int> passing;
   for (int i=0; i<fTR->NEles; i++){
      passing.push_back(i);
  }

   for (vector<int>::iterator it = passing.begin(); it != passing.end(); ){
     if (!IsGoodElId_WP80(*it)) it=passing.erase(it); else it++;
   }

   for (vector<int>::iterator it = passing.begin(); it != passing.end(); ){
     float energy=fTR->SCRaw[fTR->ElSCindex[*it]];
     float eta=fTR->SCEta[fTR->ElSCindex[*it]];
     if (fabs(eta)<1.4442) energy*=elecorr->getEtaCorrectionBarrel(eta);
     if (fabs(eta)>1.56) energy+=fTR->SCPre[fTR->ElSCindex[*it]];
     if (energy/cosh(eta)<30) it=passing.erase(it); else it++;
   }

   for (vector<int>::iterator it = passing.begin(); it != passing.end(); it++){
     float eta=fTR->SCEta[fTR->ElSCindex[*it]];
     float phi=fTR->SCPhi[fTR->ElSCindex[*it]];
     if ( (fabs(eta)>1.4442 && fabs(eta)<1.56) || (fabs(eta)>2.5) || (elecorr->isInPhiCracks(phi,eta))) evtisok=false;
   }

   for (vector<int>::iterator it = passing.begin(); it != passing.end(); it++){
     if (fTR->ElSCindex[*it]==-1) evtisok=false;
   }

   if (!evtisok) return;

   if (passing.size()<2) return;

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

  //  std::cout << "bla" << std::endl;
  
  TLorentzVector genelec[2];
  float invmassEgen=0;
  if (!isdata){
    for (int i=0; i<2; i++) genelec[i].SetPtEtaPhiE(fTR->ElGenPt[passing.at(i)],fTR->ElGenEta[passing.at(i)],fTR->ElGenPhi[passing.at(i)],fTR->ElGenE[passing.at(i)]);
    invmassEgen=(genelec[0]+genelec[1]).M();
    fHInvMassEgen->Fill(invmassEgen,weight);
  }

  fHInvMass0->Fill(invmass0,weight);
  fHInvMass15 ->Fill(invmass15,weight);
  fHInvMass16 ->Fill(invmass16,weight);
  fHInvMass17 ->Fill(invmass17,weight);
  fHInvMass20 ->Fill(invmass20,weight);




   int cat;
   
   float abseta0 = fabs(fTR->SCEta[fTR->ElSCindex[passing.at(0)]]);
   float abseta1 = fabs(fTR->SCEta[fTR->ElSCindex[passing.at(1)]]);
   if (abseta0<1.4442 && abseta1<1.4442) cat=1;
   else if (abseta0>1.56 && abseta1>1.56) cat=3;
   else if (abseta0<1.4442 && abseta1>1.56) cat=2;
   else if (abseta0>1.56 && abseta1<1.4442) cat=2;
   else cat=-1;
   
   /*
   std::cout << "cat " << cat << std::endl;
   std::cout << "mass17 " << invmass17 << " mass20 " << invmass20 << std::endl;
   //   std::cout << fTR->SCEta[fTR->ElSCindex[passing.at(0)]] << " " << fTR->SCEta[fTR->ElSCindex[passing.at(1)]] << " " << cat << std::endl;
   */

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
 
    for (int i=0; i<2; i++) {
     if (fTR->ElGenE[passing.at(i)]<0) continue;
     TLorentzVector correl17 = CorrElectron(fTR,passing.at(i),17);
     TLorentzVector correl20 = CorrElectron(fTR,passing.at(i),20);
     if (cat==1) fHErecEGen17cat1->Fill(correl17.E()/fTR->ElGenE[passing.at(i)]);
     if (cat==1) fHErecEGen20cat1->Fill(correl20.E()/fTR->ElGenE[passing.at(i)]);
     if (cat==2) fHErecEGen17cat2->Fill(correl17.E()/fTR->ElGenE[passing.at(i)]);
     if (cat==2) fHErecEGen20cat2->Fill(correl20.E()/fTR->ElGenE[passing.at(i)]);
     if (cat==3) fHErecEGen17cat3->Fill(correl17.E()/fTR->ElGenE[passing.at(i)]);
     if (cat==3) fHErecEGen20cat3->Fill(correl20.E()/fTR->ElGenE[passing.at(i)]);
    }

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

	fHErecEGen17cat1->Write();
	fHErecEGen20cat1->Write();
	fHErecEGen17cat2->Write();
	fHErecEGen20cat2->Write();
	fHErecEGen17cat3->Write();
	fHErecEGen20cat3->Write();

	fHNumPU->Write();
	
	fHistFile->Close();
	for (int i=0; i<5; i++) for (int j=0; j<4; j++)  myfile[i][j]->close();
}

TLorentzVector ZeeAnalysis::CorrElectron(TreeReader *fTR, int i, int mode){
  TLorentzVector corr;
  corr.SetPtEtaPhiE(fTR->ElPt[i],fTR->ElEta[i],fTR->ElPhi[i],fTR->ElE[i]);
  float corrE=elecorr->get_correctedenergy(fTR,i,mode);
  corr.SetE(corrE);
  corr.SetRho(corrE);
  return corr;
}
