#include "MT2tree.hh"

#include "helper/Davismt2.h"
#include "helper/TMctLib.h"

#include <vector>

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
  NJets  = 0;
  NEles  = 0;
  NMuons = 0;
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

}

void MT2tree::SetNJets(int n) {
  NJets = n;
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

Double_t MT2tree::GetMT2(double testmass, bool massive) {
  if (NJets<2)
    return -999;

  double pa[3];
  double pb[3];
  double pmiss[3];
  
  pmiss[0] = 0;
  pmiss[1] = pfmet[0].Px();
  pmiss[2] = pfmet[0].Py();
  
  pa[0] = massive ? pseudoJets[0].M() : 0;
  pa[1] = pseudoJets[0].Px();
  pa[2] = pseudoJets[0].Py();
  
  pb[0] = massive ? pseudoJets[1].M() : 0;
  pb[1] = pseudoJets[1].Px();
  pb[2] = pseudoJets[1].Py();
  
  Davismt2 *mt2 = new Davismt2();
  mt2->set_momenta(pa, pb, pmiss);
  mt2->set_mn(testmass);
  return mt2->get_mt2();
}

Double_t MT2tree::GetMCT(bool massive) {
  if (NJets<2)
    return -999;

  TLorentzVector DTM(0,0,0,0);
  TLorentzVector p1, p2;
  p1.SetXYZM(pseudoJets[0].Px(), pseudoJets[0].Py(), pseudoJets[0].Pz(), massive ? pseudoJets[0].M() : 0);
  p2.SetXYZM(pseudoJets[1].Px(), pseudoJets[1].Py(), pseudoJets[1].Pz(), massive ? pseudoJets[1].M() : 0);
  TVector2 pmiss;
  pmiss.Set(pfmet[0].Px(), pfmet[0].Py());

  TMctLib *mct = new TMctLib();
  return mct -> mctcorr(p1, p2, DTM, pmiss, 7000, 0.);
}

ClassImp(MT2Misc)
ClassImp(MT2Jet)
ClassImp(MT2tree)
