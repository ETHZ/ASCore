#include "helper/Utilities.hh"
#include "ZeeMiniTree.hh"

#include <iostream>
#include <fstream>
#include <assert.h>

using namespace std;

ZeeMiniTree::ZeeMiniTree(TreeReader *tr, std::string dataType) : UserAnalysisBase(tr), fDataType_(dataType){
	Util::SetStyle();	
	if (fDataType_ == "mc") isdata=false;
	else if (fDataType_ == "data") isdata=true; 
	else {
	  std::cout << "wrong data type" << std::endl;
	  assert(1==0);
	}
        elecorr = new EnergyCorrection("electrons");
}

ZeeMiniTree::~ZeeMiniTree(){
  delete elecorr;
}

void ZeeMiniTree::Begin(){

  cout << "Begin" << endl;

	// Define the output file of histograms
	const char* filename = "MiniTree_Zee.root";
	fMiniTree = new TFile(fOutputDir + TString(filename), "RECREATE");

	fMiniTree->cd();
	OutputTree = new TTree("Tree","Tree");

	OutputTree->Branch("event_weight",&event_weight,"event_weight/F");
	OutputTree->Branch("event_rho",&event_rho,"event_rho/F");
	OutputTree->Branch("event_nPU",&event_nPU,"event_nPU/I");
	OutputTree->Branch("event_nRecVtx",&event_nRecVtx,"event_nRecVtx/I");

	OutputTree->Branch("diel_cat",&diel_cat,"diel_cat/I");
	OutputTree->Branch("diel_mee_electron",&diel_mee_electron,"diel_mee_electron/F");
	OutputTree->Branch("diel_mee_SCdefault",&diel_mee_SCdefault,"diel_mee_SCdefault/F");
	OutputTree->Branch("diel_mee_newCorrNoCrack",&diel_mee_newCorrNoCrack,"diel_mee_newCorrNoCrack/F");
	OutputTree->Branch("diel_mee_newCorr",&diel_mee_newCorr,"diel_mee_newCorr/F");
	OutputTree->Branch("diel_mee_newCorrLocal",&diel_mee_newCorrLocal,"diel_mee_newCorrLocal/F");

	OutputTree->Branch("ellead_eta",&ellead_eta,"ellead_eta/F");
	OutputTree->Branch("eltrail_eta",&eltrail_eta,"eltrail_eta/F");
	OutputTree->Branch("ellead_px",&ellead_px,"ellead_px/F");
	OutputTree->Branch("eltrail_px",&eltrail_px,"eltrail_px/F");
	OutputTree->Branch("ellead_py",&ellead_py,"ellead_py/F");
	OutputTree->Branch("eltrail_py",&eltrail_py,"eltrail_py/F");
	OutputTree->Branch("ellead_pt",&ellead_pt,"ellead_pt/F");
	OutputTree->Branch("eltrail_pt",&eltrail_pt,"eltrail_pt/F");
	OutputTree->Branch("ellead_pz",&ellead_pz,"ellead_pz/F");
	OutputTree->Branch("eltrail_pz",&eltrail_pz,"eltrail_pz/F");
	OutputTree->Branch("ellead_energy",&ellead_energy,"ellead_energy/F");
	OutputTree->Branch("eltrail_energy",&eltrail_energy,"eltrail_energy/F");
	OutputTree->Branch("ellead_energySCdefault",&ellead_energySCdefault,"ellead_energySCdefault/F");
	OutputTree->Branch("eltrail_energySCdefault",&eltrail_energySCdefault,"eltrail_energySCdefault/F");
	OutputTree->Branch("ellead_energyNewCorrNoCrack",&ellead_energyNewCorrNoCrack,"ellead_energyNewCorrNoCrack/F");
	OutputTree->Branch("eltrail_energyNewCorrNoCrack",&eltrail_energyNewCorrNoCrack,"eltrail_energyNewCorrNoCrack/F");
	OutputTree->Branch("ellead_energyNewCorr",&ellead_energyNewCorr,"ellead_energyNewCorr/F");
	OutputTree->Branch("eltrail_energyNewCorr",&eltrail_energyNewCorr,"eltrail_energyNewCorr/F");
	OutputTree->Branch("ellead_energyNewCorrLocal",&ellead_energyNewCorrLocal,"ellead_energyNewCorrLocal/F");
	OutputTree->Branch("eltrail_energyNewCorrLocal",&eltrail_energyNewCorrLocal,"eltrail_energyNewCorrLocal/F");
	OutputTree->Branch("ellead_SCeta",&ellead_SCeta,"ellead_SCeta/F");
	OutputTree->Branch("eltrail_SCeta",&eltrail_SCeta,"eltrail_SCeta/F");
	OutputTree->Branch("ellead_fbrem",&ellead_fbrem,"ellead_fbrem/F");
	OutputTree->Branch("eltrail_fbrem",&eltrail_fbrem,"eltrail_fbrem/F");

	OutputTree->Branch("ellead_r9",&ellead_r9,"ellead_r9/F");
	OutputTree->Branch("eltrail_r9",&eltrail_r9,"eltrail_r9/F");
	OutputTree->Branch("ellead_sieie",&ellead_sieie,"ellead_sieie/F");
	OutputTree->Branch("eltrail_sieie",&eltrail_sieie,"eltrail_sieie/F");
	OutputTree->Branch("ellead_hoe",&ellead_hoe,"ellead_hoe/F");
	OutputTree->Branch("eltrail_hoe",&eltrail_hoe,"eltrail_hoe/F");
	OutputTree->Branch("ellead_brem",&ellead_brem,"ellead_brem/F");
	OutputTree->Branch("eltrail_brem",&eltrail_brem,"eltrail_brem/F");
	OutputTree->Branch("ellead_sigmaPhi",&ellead_sigmaPhi,"ellead_sigmaPhi/F");
	OutputTree->Branch("eltrail_sigmaPhi",&eltrail_sigmaPhi,"eltrail_sigmaPhi/F");
	OutputTree->Branch("ellead_sigmaEta",&ellead_sigmaEta,"ellead_sigmaEta/F");
	OutputTree->Branch("eltrail_sigmaEta",&eltrail_sigmaEta,"eltrail_sigmaEta/F");

	fHNumPU = new TH1F("NumPU","NumPU",40,0,40);
	fHNumVtx = new TH1F("NumVtx","NumVtx",40,0,40);
	

	
	  cout << "Tree and histos created" << endl;
	

}

void ZeeMiniTree::Analyze(){

  //cout << "Analyze this event" << endl;


  //cout << "A" << endl;

  float weight;
  if (!isdata) weight = GetPUWeight(fTR->PUnumInteractions);
  else weight=1;

  if (!isdata) fHNumPU->Fill(fTR->PUnumInteractions,weight);
  fHNumVtx->Fill(fTR->NVrtx,weight);

  //cout << "B" << endl;

  bool evtisok = true;

  vector<int> passing;
   for (int i=0; i<fTR->NEles; i++){
      passing.push_back(i);
  }

   //cout << "C" << endl;

   for (vector<int>::iterator it = passing.begin(); it != passing.end(); ){
     if (!IsGoodElId_WP80(*it)) it=passing.erase(it); else it++;
   }

   //cout << "D" << endl;

   for (vector<int>::iterator it = passing.begin(); it != passing.end(); ){
     float energy=fTR->SCRaw[fTR->ElSCindex[*it]];
     float eta=fTR->SCEta[fTR->ElSCindex[*it]];
     if (fabs(eta)<1.4442) energy*=elecorr->getEtaCorrectionBarrel(eta);
     if (fabs(eta)>1.56) energy+=fTR->SCPre[fTR->ElSCindex[*it]];
     if (energy/cosh(eta)<30) it=passing.erase(it); else it++;
   }

   //cout << "E" << endl;

   for (vector<int>::iterator it = passing.begin(); it != passing.end(); it++){
     float eta=fTR->SCEta[fTR->ElSCindex[*it]];
     float phi=fTR->SCPhi[fTR->ElSCindex[*it]];
     if ( (fabs(eta)>1.4442 && fabs(eta)<1.56) || (fabs(eta)>2.5) || (elecorr->isInPhiCracks(phi,eta))) evtisok=false;
   }

   //cout << "F" << endl;

   for (vector<int>::iterator it = passing.begin(); it != passing.end(); it++){
     if (fTR->ElSCindex[*it]==-1) evtisok=false;
   }

   if (!evtisok) return;

   if (passing.size()<2) return;

   //cout << "G" << endl;

   TLorentzVector elec[2];
   for (int i=0; i<2; i++){
     elec[i].SetPtEtaPhiE(fTR->ElPt[passing.at(i)],fTR->ElEta[passing.at(i)],fTR->ElPhi[passing.at(i)],fTR->ElE[passing.at(i)]);
   }

   //cout << "H" << endl;

  bool masswindow=false;
  if (fabs((elec[0]+elec[1]).M()-91.2)<30) masswindow=true;
  if (!masswindow) return;

  //cout << "electrons selected" << endl;

  
  float invmass0=(elec[0]+elec[1]).M();
  float invmass15=(CorrElectron(fTR,passing.at(0),15)+CorrElectron(fTR,passing.at(1),15)).M();
  float invmass16=(CorrElectron(fTR,passing.at(0),16)+CorrElectron(fTR,passing.at(1),16)).M();
  float invmass17=(CorrElectron(fTR,passing.at(0),17)+CorrElectron(fTR,passing.at(1),17)).M(); //our correction (material+crack)
  float invmass18=(CorrElectron(fTR,passing.at(0),18)+CorrElectron(fTR,passing.at(1),18)).M(); //our correction material only
  float invmass20=(CorrElectron(fTR,passing.at(0),20)+CorrElectron(fTR,passing.at(1),20)).M(); //default SC energy




  //  std::cout << "bla" << std::endl;
  
  TLorentzVector genelec[2];
  float invmassEgen=0;
  if (!isdata){
    for (int i=0; i<2; i++) genelec[i].SetPtEtaPhiE(fTR->ElGenPt[passing.at(i)],fTR->ElGenEta[passing.at(i)],fTR->ElGenPhi[passing.at(i)],fTR->ElGenE[passing.at(i)]);
    invmassEgen=(genelec[0]+genelec[1]).M();
   }

  

   int cat;
   
   float abseta0 = fabs(fTR->SCEta[fTR->ElSCindex[passing.at(0)]]);
   float abseta1 = fabs(fTR->SCEta[fTR->ElSCindex[passing.at(1)]]);
   if (abseta0<1.4442 && abseta1<1.4442) cat=1;
   else if (abseta0>1.56 && abseta1>1.56) cat=3;
   else if (abseta0<1.4442 && abseta1>1.56) cat=2;
   else if (abseta0>1.56 && abseta1<1.4442) cat=2;
   else cat=-1;
 


  //HERE GOES MY STUFF  
  // cout << "here goes my stuff" << endl;

  event_weight = weight;
  event_rho = fTR->Rho;
  if (!isdata) event_nPU = fTR->PUnumInteractions;
  event_nRecVtx = fTR->NVrtx;
  
  diel_cat = cat;
  diel_mee_electron = invmass0;
  diel_mee_SCdefault = invmass20;
  diel_mee_newCorrNoCrack = invmass18;
  diel_mee_newCorr = invmass15;
  diel_mee_newCorrLocal = invmass16;

  ellead_eta = fTR->ElEta[passing.at(0)];
  eltrail_eta = fTR->ElEta[passing.at(1)];
  
  ellead_px = fTR->ElPx[passing.at(0)];
  eltrail_px = fTR->ElPx[passing.at(1)];
  ellead_py = fTR->ElPy[passing.at(0)];
  eltrail_py = fTR->ElPy[passing.at(1)];
  ellead_pt = fTR->ElPt[passing.at(0)];
  eltrail_pt = fTR->ElPt[passing.at(1)];
  ellead_pz = fTR->ElPz[passing.at(0)];
  eltrail_pz = fTR->ElPz[passing.at(1)];
  ellead_energy = fTR->ElE[passing.at(0)];
  eltrail_energy = fTR->ElE[passing.at(1)];

  ellead_SCeta = fTR->SCEta[fTR->ElSCindex[passing.at(0)]];
  eltrail_SCeta = fTR->SCEta[fTR->ElSCindex[passing.at(1)]];
  ellead_fbrem = fTR->Elfbrem[passing.at(0)];
  eltrail_fbrem = fTR->Elfbrem[passing.at(1)];


  ellead_energySCdefault = CorrElectron(fTR,passing.at(0),20).E();
  eltrail_energySCdefault = CorrElectron(fTR,passing.at(1),20).E();
  ellead_energyNewCorrNoCrack = CorrElectron(fTR,passing.at(0),18).E();
  eltrail_energyNewCorrNoCrack = CorrElectron(fTR,passing.at(1),18).E();
  ellead_energyNewCorr = CorrElectron(fTR,passing.at(0),15).E();
  eltrail_energyNewCorr = CorrElectron(fTR,passing.at(1),15).E();
  ellead_energyNewCorrLocal = CorrElectron(fTR,passing.at(0),16).E();
  eltrail_energyNewCorrLocal = CorrElectron(fTR,passing.at(1),16).E();

  ellead_r9 = fTR->SCR9[fTR->ElSCindex[passing.at(0)]];
  eltrail_r9 = fTR->SCR9[fTR->ElSCindex[passing.at(1)]];
  ellead_sieie = fTR->ElSigmaIetaIeta[passing.at(0)];
  eltrail_sieie = fTR->ElSigmaIetaIeta[passing.at(1)];
  ellead_hoe = fTR->ElHcalOverEcal[passing.at(0)];
  eltrail_hoe = fTR->ElHcalOverEcal[passing.at(1)];
  ellead_brem =  fTR->SCBrem[fTR->ElSCindex[passing.at(0)]];
  eltrail_brem = fTR->SCBrem[fTR->ElSCindex[passing.at(1)]];
  ellead_sigmaPhi = fTR->SCPhiWidth[fTR->ElSCindex[passing.at(0)]];
  eltrail_sigmaPhi = fTR->SCPhiWidth[fTR->ElSCindex[passing.at(1)]];
  ellead_sigmaEta = fTR->SCEtaWidth[fTR->ElSCindex[passing.at(0)]];
  eltrail_sigmaEta = fTR->SCEtaWidth[fTR->ElSCindex[passing.at(1)]];

  OutputTree->Fill();

 
 
 

}

void ZeeMiniTree::End(){
	fMiniTree->cd();
	OutputTree->Write();	
	fHNumPU->Write();
	fHNumVtx->Write();
	
	fMiniTree->Close();

}

TLorentzVector ZeeMiniTree::CorrElectron(TreeReader *fTR, int i, int mode){
  TLorentzVector corr;
  corr.SetPtEtaPhiE(fTR->ElPt[i],fTR->ElEta[i],fTR->ElPhi[i],fTR->ElE[i]);
  float corrE=elecorr->get_correctedenergy(fTR,i,mode);
  corr.SetE(corrE);
  corr.SetRho(corrE);
  return corr;
}
