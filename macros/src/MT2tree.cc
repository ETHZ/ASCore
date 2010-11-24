#include "MT2tree.hh"

#include "helper/Davismt2.h"
#include "helper/TMctLib.h"

#include <vector>
#include "helper/Hemisphere.hh"
#include "TLorentzVector.h"
#include "TMath.h"

using std::vector;

MT2Misc::MT2Misc(){
  Reset();
}

MT2Misc::~MT2Misc(){
}

void MT2Misc::Reset() {
  Run                     = -1;	  
  Event		  	  = -1;	  
  LumiSection		  = -1;	  
  LeptConfig		  = -1;	  
  NJetsEta5Pt20           = -1;
  IsCleanJetEvent	  =  0;	  
  PseudoJetMT2		  = -99999.99;
  PseudoJetMCT		  = -99999.99;
  PseudoJet1Pt		  = -99999.99;
  PseudoJet2Pt		  = -99999.99;
  PseudojetAlphaT	  = -99999.99;
  Vectorsumpt		  = -99999.99;
  PFMETsign		  = -99999.99;
  HT			  = -99999.99;
  DPhiMhtMpt              = -99999.99;
}

MT2Jet::MT2Jet(){
  Reset();
}

MT2Jet::~MT2Jet(){
}

void MT2Jet::Reset() {
  lv.SetPxPyPzE(0, 0, 0, 0);
  bTagProbTCHE  = -99999.99;
  bTagProbTCHP  = -99999.99;
  bTagProbSSVHE = -99999.99;
  bTagProbSSVHP = -99999.99;
  inHemisphere  = -1;

  isPFIDLoose   = 0;
  isPFIDMedium  = 0;
  isPFIDTight   = 0;
  ChHadFrac     = -99999.99; 
  NeuHadFrac    = -99999.99; 
  ChEmFrac      = -99999.99;
  NeuEmFrac     = -99999.99; 
  ChMult        = -1; 
  NeuMult       = -1; 
  NConstituents = -1; 

}

void MT2Jet::SetLV(const TLorentzVector v) {
  lv = v;
}

MT2tree::MT2tree(){
  Reset();
}

MT2tree::~MT2tree(){
}

void MT2tree::Reset() {
  NJets         = 0;
  NEles         = 0;
  NMuons        = 0;
  NJetsIDLoose  = 0;
  NJetsIDMedium = 0;
  NJetsIDTight  = 0;

  misc.Reset();
  for (int i = 0; i < m_jetSize; ++i) {
    jet[i].Reset();
  }
  for (int i = 0; i < m_eleSize; ++i) {
    ele[i].SetPxPyPzE(0., 0., 0., 0.);
  }
  for (int i = 0; i < m_muoSize; ++i) {
    muo[i].SetPxPyPzE(0., 0., 0., 0.);
  }
  pfmet     [0].SetPxPyPzE(0., 0., 0., 0.);
  MPT       [0].SetPxPyPzE(0., 0., 0., 0.);
  pseudoJets[0].SetPxPyPzE(0., 0., 0., 0.);
  pseudoJets[1].SetPxPyPzE(0., 0., 0., 0.);
  MHTloose  [0].SetPxPyPzE(0., 0., 0., 0.);
}

void MT2tree::SetNJets(int n) {
  NJets = n;
}

void MT2tree::SetNJetsIDLoose(int n) {
  NJetsIDLoose = n;
}

void MT2tree::SetNJetsIDMedium(int n) {
  NJetsIDMedium = n;
}

void MT2tree::SetNJetsIDTight(int n) {
  NJetsIDTight = n;
}

void MT2tree::SetNEles(int n) {
  NEles = n;
}

void MT2tree::SetNMuons(int n) {
  NMuons = n;
}

Double_t MT2tree::PseudoJetDPhi() {
  if (NJets<2)
    return -999;
  return TMath::Abs(pseudoJets[0].DeltaPhi(pseudoJets[1]));
}

Double_t MT2tree::PseudoJetAngle() {
  if (NJets<2)
    return -999;
  return pseudoJets[0].Angle(pseudoJets[1].Vect());
}

Double_t MT2tree::JetsDPhi(int j1, int j2){
  if (j1>=NJets || j2>=NJets)
    return -999;
  return TMath::Abs(jet[j1].lv.DeltaPhi(jet[j2].lv));
}

Double_t MT2tree::MetJetDPhi(int ijet) {
  if (ijet>=NJets)
    return -999;
  return TMath::Abs(jet[ijet].lv.DeltaPhi(pfmet[0]));
}

Double_t MT2tree::MinMetJetDPhi() {
  if (NJets<1)
    return -999;
  Double_t minDPhi=10;
  for (int i=0; i<NJets; i++){
    Double_t dphi = TMath::Abs(jet[i].lv.DeltaPhi(pfmet[0]));
    if(dphi<minDPhi)  minDPhi=dphi;
  }
  return minDPhi;
}

Int_t MT2tree::MinMetJetDPhiIndex() {
  if (NJets<1)
    return -999;
  Double_t minDPhi=10;
  Int_t    imin;
  for (int i=0; i<NJets; i++){
    Double_t dphi = TMath::Abs(jet[i].lv.DeltaPhi(pfmet[0]));
    if(dphi<minDPhi)  {
      minDPhi=dphi;
      imin = i;
    }
  }
  return imin;
}

Int_t MT2tree::GetNJets(double minJPt, std::string PFJID){
  int njets=0;
  for(int i=0; i<NJets; ++i){
	if(PFJID=="loose"  && jet[i].isPFIDLoose ==false) continue;
	if(PFJID=="medium" && jet[i].isPFIDMedium==false) continue;
	if(PFJID=="tight"  && jet[i].isPFIDTight ==false) continue;
	if(jet[i].lv.Pt() < minJPt )                      continue;
	njets++;
  }
  return njets;
}

Double_t MT2tree::GetMT2(double testmass, bool massive) {
  if (NJetsIDLoose<2)
    return -999;

  return CalcMT2(testmass, massive, pseudoJets[0], pseudoJets[1], pfmet[0]);
}

Double_t MT2tree::GetMT2Leading(double testmass, bool massive, std::string PFJID){
  
  TLorentzVector leadingJets[2]; 
  int index=0; 
  for(int i=0; i<NJets; ++i){
	if(PFJID=="loose"  && jet[i].isPFIDLoose ==false) continue;
	if(PFJID=="medium" && jet[i].isPFIDMedium==false) continue;
	if(PFJID=="tight"  && jet[i].isPFIDTight ==false) continue;
	if(index >1                                     ) continue;
	leadingJets[index] = jet[i].lv;
	index++;
  }
  if(index!=2) return -999;
  return CalcMT2(testmass, massive, leadingJets[0], leadingJets[1], pfmet[0]);

}

Double_t MT2tree::GetMT2Hemi(double testmass, bool massive, std::string PFJID, double minJPt, int hemi_association) {
  vector<float> px, py, pz, E;
  for(int i=0; i<NJets; ++i){
	if(PFJID=="loose"  && jet[i].isPFIDLoose ==false) continue;
	if(PFJID=="medium" && jet[i].isPFIDMedium==false) continue;
	if(PFJID=="tight"  && jet[i].isPFIDTight ==false) continue;
	if(jet[i].lv.Pt() < minJPt )                      continue;
  	px.push_back(jet[i].lv.Px());
	py.push_back(jet[i].lv.Py());
	pz.push_back(jet[i].lv.Pz());
	 E.push_back(jet[i].lv.E());
  }
		
  if (px.size()<2) return -999;
  

  // get hemispheres (seed 2: max inv mass, association method: default 3 = minimal lund distance)
  Hemisphere* hemi = new Hemisphere(px, py, pz, E, 2, hemi_association);
  vector<int> grouping = hemi->getGrouping();

  TLorentzVector pseudojet1(0.,0.,0.,0.);
  TLorentzVector pseudojet2(0.,0.,0.,0.);
	
  for(int i=0; i<px.size(); ++i){
	if(grouping[i]==1){
		pseudojet1.SetPx(pseudojet1.Px() + px[i]);
		pseudojet1.SetPy(pseudojet1.Py() + py[i]);
		pseudojet1.SetPz(pseudojet1.Pz() + pz[i]);
		pseudojet1.SetE( pseudojet1.E()  + E[i]);	
	}else if(grouping[i] == 2){
		pseudojet2.SetPx(pseudojet2.Px() + px[i]);
		pseudojet2.SetPy(pseudojet2.Py() + py[i]);
		pseudojet2.SetPz(pseudojet2.Pz() + pz[i]);
		pseudojet2.SetE( pseudojet2.E()  + E[i]);
	}
  }
  delete hemi;
 
  return CalcMT2(testmass, massive, pseudojet1, pseudojet2, pfmet[0]); 
}

Double_t MT2tree::CalcMT2(double testmass, bool massive, TLorentzVector visible1, TLorentzVector visible2, TLorentzVector MET ){
  
  double pa[3];
  double pb[3];
  double pmiss[3];
  
  pmiss[0] = 0;
  pmiss[1] = MET.Px();
  pmiss[2] = MET.Py();
  
  pa[0] = massive ? visible1.M() : 0;
  pa[1] = visible1.Px();
  pa[2] = visible1.Py();
  
  pb[0] = massive ? visible2.M() : 0;
  pb[1] = visible2.Px();
  pb[2] = visible2.Py();
  
  Davismt2 *mt2 = new Davismt2();
  mt2->set_momenta(pa, pb, pmiss);
  mt2->set_mn(testmass);
  Double_t MT2=mt2->get_mt2();
  delete mt2;
  return MT2;

}

Double_t MT2tree::GetMCT(bool massive) {
  if (NJetsIDLoose<2)
    return -999;

  TLorentzVector DTM(0,0,0,0);
  TLorentzVector p1, p2;
  p1.SetXYZM(pseudoJets[0].Px(), pseudoJets[0].Py(), pseudoJets[0].Pz(), massive ? pseudoJets[0].M() : 0);
  p2.SetXYZM(pseudoJets[1].Px(), pseudoJets[1].Py(), pseudoJets[1].Pz(), massive ? pseudoJets[1].M() : 0);
  TVector2 pmiss;
  pmiss.Set(pfmet[0].Px(), pfmet[0].Py());

  TMctLib *mct = new TMctLib();
  Double_t MCT =  mct -> mctcorr(p1, p2, DTM, pmiss, 7000, 0.);
  delete mct;
  return MCT;
}

ClassImp(MT2Misc)
ClassImp(MT2Jet)
ClassImp(MT2tree)
