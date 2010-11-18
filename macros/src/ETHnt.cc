#include "ETHnt.hh"

#include <vector>

#include "TLorentzVector.h"

using std::vector;

ETHMisc::ETHMisc(){
  Reset();
}

ETHMisc::~ETHMisc(){
}

void ETHMisc::Reset() {
  Run                     = -1;	  
  Event		  	  = -1;	  
  LumiSection		  = -1;	  
  Weight		  = -99999.99;
  LeptConfig		  = -1;	  
  IsCleanMultiJetEvent	  =  0;	  
  IsCleanJetEvent	  =  0;	  
  R12R21		  = -1;	  
  NJetsPt50Eta25	  = -1;	  
  PseudoJetMT2		  = -99999.99;
  PseudoJetMCT		  = -99999.99;
  PseudoJet1Pt		  = -99999.99;
  PseudoJet2Pt		  = -99999.99;
  PseudojetAlphaT	  = -99999.99;
  Vectorsumpt		  = -99999.99;
  PFMET		  	  = -99999.99;
  PFMETphi		  = -99999.99;
  PFMETsign		  = -99999.99;
  LeadingJetEta	  	  = -99999.99;
  DPhiJ1Met		  = -99999.99;
  DPhiJ2Met		  = -99999.99;
  PseudoJetMT2AxisdPhi	  = -99999.99;
  R1221min		  = -99999.99;
  MPT_sel		  = -99999.99;
  MPT			  = -99999.99;
  MHT			  = -99999.99;
  HT			  = -99999.99;
  DPhiMhtMpt		  = -99999.99;
}

ETHJet::ETHJet(){
  Reset();
}

ETHJet::~ETHJet(){
}

void ETHJet::Reset() {
  lv.SetPxPyPzE(0, 0, 0, 0);
  bTagProbTCHE  = -99999.99;
  bTagProbTCHP  = -99999.99;
  bTagProbSSVHE = -99999.99;
  bTagProbSSVHP = -99999.99;
  inHemisphere  = -1;
}

void ETHJet::SetLV(const TLorentzVector v) {
  lv = v;
}

ETHnt::ETHnt(){
  Reset();
}

ETHnt::~ETHnt(){
}

void ETHnt::Reset() {
  NJets  = 0;
  NEles  = 0;
  NMuons = 0;
  misc.Reset();
  for (int i = 0; i < m_jetSize; ++i) {
    jet[i].Reset();
  }
  for (int i = 0; i < m_eleSize; ++i) {
    ele[i].SetPxPyPzE(0, 0, 0, 0);
  }
  for (int i = 0; i < m_muoSize; ++i) {
    muo[i].SetPxPyPzE(0, 0, 0, 0);
  }
  pfmet     [0].SetPxPyPzE(0, 0, 0, 0);
  pseudoJets[0].SetPxPyPzE(0, 0, 0, 0);
  pseudoJets[1].SetPxPyPzE(0, 0, 0, 0);

}

void ETHnt::SetNJets(int n) {
  NJets = n;
}

void ETHnt::SetNEles(int n) {
  NEles = n;
}

void ETHnt::SetNMuons(int n) {
  NMuons = n;
}

Double_t ETHnt::PseudoJetDPhi() {
  if (NJets<2)
    return -999;
  return pseudoJets[0].DeltaPhi(pseudoJets[1]);
}

Double_t ETHnt::PseudoJetAngle() {
  if (NJets<2)
    return -999;
  return abs(pseudoJets[0].Angle(pseudoJets[1].Vect()));
}

ClassImp(ETHMisc)
ClassImp(ETHJet)
ClassImp(ETHnt)
