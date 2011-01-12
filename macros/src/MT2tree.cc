#include "MT2tree.hh"

#include "helper/Davismt2.h"
#include "helper/TMctLib.h"

#include <vector>
#include "helper/Hemisphere.hh"
#include "helper/Utilities.hh"
#include "TLorentzVector.h"
#include "TMath.h"

using std::vector;

// MT2Misc -----------------------------------
MT2Misc::MT2Misc(){
  Reset();
}

MT2Misc::~MT2Misc(){
}

void MT2Misc::Reset() {
  Run                     = -1;	  
  Event		  	  = -1;	  
  LumiSection		  = -1;	  
  NVertices               = -1;
  LeptConfig		  = -1;	  
  NJetsEta5Pt20           = -1;
  PassJetID               = -1;
  Jet0Pass                = -1;
  Jet1Pass                = -1;
  HBHENoiseFlag           =  0;
  MT2                     = -99999.99;
  MT2leading              = -99999.99;
  MT2noISR                = -99999.99;
  MCT                     = -99999.99;
  AlphaT                  = -99999.99;
  MET                     = -99999.99;
  METPhi                  = -99999.99;
  Vectorsumpt		  = -99999.99;
  PFMETsign		  = -99999.99;
  HT			  = -99999.99;
  DPhiMhtMpt              = -99999.99;
  MinMetJetDPhi           = -99999.99;
  EcalDeadCellBEFlag      = -1;
  NECALGapClusters        = -1;
  for(int i=0; i<50; ++i){
    EcalGapClusterSize[i] = -1;
    EcalGapBE[i]          = -99999.99;    
  }
}

// MT2Jet -----------------------------------
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

Bool_t MT2Jet::IsGoodPFJet(double minJPt, double maxJEta, int PFJID) {

  double pt = lv.Pt();
  double eta = lv.Eta();
  if ( pt < minJPt || fabs(eta) > maxJEta )     return false;
  
  switch (PFJID) {
  case 3:               // TIGHT
    if ( !(NeuEmFrac     < 0.90) ) return false;
    if ( !(NeuHadFrac    < 0.90) ) return false;
    // break;   // No break: medium contains tight -> check common cuts
  case 2:               // MEDIUM
    if ( !(NeuEmFrac     < 0.95) ) return false;
    if ( !(NeuHadFrac    < 0.95) ) return false;
    // break;   // No break: loose contains medium -> check common cuts
  case 1:               // LOOSE
    if ( !(NConstituents > 1) )    return false;
    if ( !(NeuEmFrac     < 0.99) ) return false;
    if ( !(NeuHadFrac    < 0.99) ) return false;
    if ( fabs(eta) < 2.4 ) { // Cuts for |eta|<2.4
      if ( !(ChEmFrac  < 0.99) )  return false;
      if ( !(ChHadFrac > 0.00) )  return false;
      if ( !(ChMult    > 0   ) )  return false;
    }
    break;
  default:
    // None of the above. Do we want any default cut?
    break;
  }

  return true;
}


// MT2Hemisphere -----------------------------------
MT2Hemi::MT2Hemi(){
  Reset();
}

MT2Hemi::~MT2Hemi(){
}

void MT2Hemi::Reset(){
  seed_method    = -1;
  assoc_method   = -1;
  MT2            = -99999.99;
  MCT            = -99999.99;
  AlphaT         = -99999.99;
  minDHT         = -99999.99;
  maxDR          = -99999.99;
  dPhi           = -99999.99;
  for(int i=0; i<m_jetSize; ++i){
  	jindices1[i]  =-1;
  	jindices2[i]  =-1;
  }
  for(int i=0; i<m_eleSize; ++i){
  	eleindices1[i]=-1;
  	eleindices2[i]=-1;
  }
  for(int i=0; i<m_muoSize; ++i){
  	muoindices1[i]=-1;
  	muoindices2[i]=-1;
  }
  
  lv1. SetPxPyPzE(0, 0, 0, 0);
  lv2. SetPxPyPzE(0, 0, 0, 0);
  UTM. SetPxPyPzE(0, 0, 0, 0); 
}


// MT2GenLept -----------------------------------
MT2GenLept::MT2GenLept(){
  Reset();
}

MT2GenLept::~MT2GenLept(){
}

void MT2GenLept::Reset(){
  lv.  SetPxPyPzE(0, 0, 0, 0);

  ID       = -999;
  MID      = -999;
  MStatus  = -999;
  GMID     = -999;
  GMStatus = -999;
}

// MT2Muon -----------------------------------
MT2Muon::MT2Muon(){
  Reset();
}

MT2Muon::~MT2Muon(){
}

void MT2Muon::Reset() {
  lv.SetPxPyPzE(0, 0, 0, 0);
  isTight       = 0;
}

void MT2Muon::SetLV(const TLorentzVector v) {
  lv = v;
}

// MT2Elec -----------------------------------
MT2Elec::MT2Elec(){
  Reset();
}

MT2Elec::~MT2Elec(){
}

void MT2Elec::Reset() {
  lv.SetPxPyPzE(0, 0, 0, 0);
  isTight       = 0;
}

void MT2Elec::SetLV(const TLorentzVector v) {
  lv = v;
}

// MT2tree ----------------------------------
MT2tree::MT2tree(){
  Reset();
}

MT2tree::~MT2tree(){
}

void MT2tree::Reset() {
  NJets         = 0;
  NJetsIDLoose  = 0;
  NJetsIDMedium = 0;
  NJetsIDTight  = 0;
  NEles         = 0;
  NElesLoose    = 0;
  NMuons        = 0;
  NMuonsLoose   = 0;

  misc.Reset();
  for (int i = 0; i < m_jetSize; ++i) {
    jet[i].Reset();
  }
  for (int i = 0; i < m_eleSize; ++i) {
    ele[i].Reset();
  }
  for (int i = 0; i < m_muoSize; ++i) {
    muo[i].Reset();
  }
  for (int i = 0; i < m_genleptSize; ++i) {
    genlept[i].Reset();
  }
  for (int i = 0; i < m_hemiSize; ++i) {
    hemi[i].Reset();
  }
  pfmet     [0].SetPxPyPzE(0., 0., 0., 0.);
  MPT       [0].SetPxPyPzE(0., 0., 0., 0.);
  MHTloose  [0].SetPxPyPzE(0., 0., 0., 0.);
  MHT       [0].SetPxPyPzE(0., 0., 0., 0.);
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

void MT2tree::SetNElesLoose(int n) {
  NElesLoose = n;
}

void MT2tree::SetNMuons(int n) {
  NMuons = n;
}

void MT2tree::SetNMuonsLoose(int n) {
  NMuonsLoose = n;
}

Double_t MT2tree::GetMinR12R21(int PFJID, double minJPt, double maxJEta, int met){
	TLorentzVector MET(0., 0., 0., 0.);
	if(met==1)      MET=pfmet[0];
        else if(met==2) MET=MHT[0];
	else            return -900;

	vector<int> indices;
	for(int i=0; i<NJets; ++i){
		if( jet[i].IsGoodPFJet(minJPt, maxJEta, PFJID)==false) continue;
		indices.push_back(i);
	}
	if(indices.size()<2) return -999;
	double d_phi1 = fabs(Util::DeltaPhi(jet[indices[0]].lv.Phi(), MET.Phi()));	
	double d_phi2 = fabs(Util::DeltaPhi(jet[indices[1]].lv.Phi(), MET.Phi()));	
	double R12    = sqrt( d_phi1*d_phi1 + (TMath::Pi()-d_phi2)*(TMath::Pi()-d_phi2) );
	double R21    = sqrt( d_phi2*d_phi2 + (TMath::Pi()-d_phi1)*(TMath::Pi()-d_phi1) );

	if(R12<R21) return R12;
	else        return R21;
}

Bool_t MT2tree::PassJetID(double minJPt, double maxJEta, int PFJID) {
	int njets=0;
		for(int i=0; i<NJets; ++i){
			if(jet[i].lv.Pt() >= minJPt && fabs(jet[i].lv.Eta()) <= maxJEta &&
				jet[i].IsGoodPFJet(minJPt,maxJEta,PFJID)==false)   return false;
		}
	return true;
}

Double_t MT2tree::JetsInvMass(int j1, int j2){
	if(NJets < 2) return -999.99;
	TLorentzVector sum = jet[j1].lv + jet[j2].lv;
	return sum.M();
}

Double_t MT2tree::JetsDPhi(int j1, int j2, int PFJID){
  std::vector<int> indices;
  for(int i = 0; i<NJets; ++i){
  	if(jet[i].isPFIDLoose ==false && PFJID==1  )continue;
  	if(jet[i].isPFIDMedium==false && PFJID==2  )continue;
  	if(jet[i].isPFIDTight ==false && PFJID==3  )continue;
	indices.push_back(i);
  }
  if (j1>=indices.size() || j2>=indices.size())  return -999;
  return TMath::Abs(jet[indices[j1]].lv.DeltaPhi(jet[indices[j2]].lv));
}

Double_t MT2tree::MetJetDPhi(int ijet, int PFJID, int met) {
  TLorentzVector MET(0., 0., 0., 0.);
  if(met==1)      MET = pfmet[0];
  else if(met==2) MET = MHTloose[0];
  else            return -999;

  std::vector<int> indices;
  for(int i = 0; i<NJets; ++i){
  	if(jet[i].isPFIDLoose ==false && PFJID==1  )continue;
  	if(jet[i].isPFIDMedium==false && PFJID==2  )continue;
  	if(jet[i].isPFIDTight ==false && PFJID==3  )continue;
	indices.push_back(i);
  }
  if (ijet>=indices.size())  return -999;
  return TMath::Abs(jet[indices[ijet]].lv.DeltaPhi(MET));
}

Double_t MT2tree::MinMetJetDPhi(int PFJID, double minJPt, double maxJEta, int met) {
// Attention: electrons and muons are not considered for minDPhi
  TLorentzVector MET(0., 0., 0., 0.);
  if(met==1)      MET = pfmet[0];
  else if(met==2) MET = MHTloose[0];
  else            return -999;

  int index = MinMetJetDPhiIndex(PFJID, minJPt, maxJEta, met);
  return TMath::Abs(jet[index].lv.DeltaPhi(MET));
}

Int_t MT2tree::MinMetJetDPhiIndex(int PFJID, double minJPt, double maxJEta, int met) {
// Attention: electrons and muons are not considered for minDPhi
  TLorentzVector MET(0., 0., 0., 0.);
  if(met==1)      MET = pfmet[0];
  else if(met==2) MET = MHTloose[0];
  else            return -999;

  std::vector<int> indices;
  for(int i = 0; i<NJets; ++i){
	if(jet[i].IsGoodPFJet(minJPt,maxJEta,PFJID)==false) continue;
	indices.push_back(i);
  }
  if(indices.size()<1)  return -999;
  Double_t minDPhi=10;
  Int_t    imin;
  for (int i=0; i<indices.size(); i++){
    Double_t dphi = TMath::Abs(jet[indices[i]].lv.DeltaPhi(MET));
    if(dphi<minDPhi)  {
      minDPhi=dphi;
      imin = indices[i];
    }
  }
  return imin;
}

Int_t MT2tree::GetNjets(double minJPt, double maxJEta, int PFJID){
  int njets=0;
  for(int i=0; i<NJets; ++i){
	if(jet[i].IsGoodPFJet(minJPt,maxJEta,PFJID)==false) continue;
	njets++;
  }
  return njets;
}

Double_t MT2tree::JetPt(int ijet, int PFJID) {
  int index = GetJetIndex(ijet, PFJID);
  if ( index < 0 )       return -999.;
  return jet[index].lv.Pt();
}

Int_t MT2tree::GetJetIndex(int ijet, int PFJID) {
  std::vector<int> indices;
  for(int i = 0; i<NJets; ++i){
  	if     (jet[i].isPFIDLoose ==false && PFJID==1  )continue;
  	else if(jet[i].isPFIDMedium==false && PFJID==2  )continue;
  	else if(jet[i].isPFIDTight ==false && PFJID==3  )continue;
	indices.push_back(i);
  }
  if (ijet>=indices.size())  return -9;
  return indices[ijet];
}



Int_t MT2tree::JetIsInHemi(int jindex, int hemi_seed, int hemi_association, float MaxDR){
  // returns 1 or 2 depending if jet is associated to hemi 1 or 2
  // returns 0 if jet is not associated to either of the two hemispheres

  // hemi_seed   = 2: max invariant mass, 4: two leading jet
  // hemi_assoc  = 2: min invariant mass, 3: minial lund distance, 
  // MaxDR: if >0 RejectISRDRmax(MaxDR) is called for hemi            

  
  vector<float> px, py, pz, E;
  vector<int>   jsel; // contains indices of all jets fed into hemisphere algo
  for(int i=0; i<NJets; ++i){
	if(jet[i].isPFIDLoose ==false) continue;
  	px.push_back(jet[i].lv.Px());
	py.push_back(jet[i].lv.Py());
	pz.push_back(jet[i].lv.Pz());
	 E.push_back(jet[i].lv.E());
	jsel.push_back(i);
  }
		
  if (px.size()<2) return -10;

  // get hemispheres (seed 2: max inv mass, association method: default 3 = minimal lund distance)
  Hemisphere* hemi = new Hemisphere(px, py, pz, E, hemi_seed, hemi_association);
  hemi->SetDebug(0);
  if(MaxDR > 0) hemi->RejectISRDRmax(MaxDR);
  vector<int> grouping = hemi->getGrouping();
  delete hemi;
  
  vector<int>    jused; // contains indices of jets used for pseudojets
  for(int i=0; i<px.size(); ++i){
	if(grouping[i]==1){
		if(jsel[i]==jindex) return 1;
	}else if(grouping[i] == 2){
		if(jsel[i]==jindex) return 2;
	}
  }
  return -1;

}

Double_t MT2tree::GetMT2Leading(double testmass, bool massive, int PFJID, int met){
  TLorentzVector MET(0., 0., 0., 0.);
  if(met==1)      MET = pfmet[0];
  else if(met==2) MET = MHTloose[0];
  else            return -999;

  
  TLorentzVector leadingJets[2]; 
  int index=0; 
  for(int i=0; i<NJets; ++i){
	if(PFJID==1  && jet[i].isPFIDLoose ==false) continue;
	if(PFJID==2  && jet[i].isPFIDMedium==false) continue;
	if(PFJID==3  && jet[i].isPFIDTight ==false) continue;
	if(index >1                               ) continue;
	leadingJets[index] = jet[i].lv;
	index++;
  }
  if(index!=2) return -999;
  cout << "leadingjets " <<leadingJets[0].Pt() << " " << leadingJets[1].Pt() << endl;
  return CalcMT2(testmass, massive, leadingJets[0], leadingJets[1], MET);

}

Double_t MT2tree::GetMT2HemiNoISR(bool massive, int hemi_seed, int hemi_association, float MaxDR, int met) {
  // hemi_seed   = 2: max invariant mass, 4: two leading jet
  // hemi_assoc  = 2: min invariant mass, 3: minial lund distance, 
  // MaxDR: if >0 RejectISRDRmax(MaxDR) is called for hemi            
  // met ==3 option adds jets > jet_for_met_threshold GeV rejected by the chosen hemispheres to the MET
  float jet_for_met_threshold=20;

  TLorentzVector MET(0., 0., 0., 0.);
  if(met==1 || met==3)  MET = pfmet[0];
  else if(met==2)       MET = MHTloose[0];
  else return -999;
 

  vector<float> px, py, pz, E;
  vector<int>   jsel; // contains indices of all jets fed into hemisphere algo
  for(int i=0; i<NJets; ++i){
	if(jet[i].isPFIDLoose ==false) continue;
  	px.push_back(jet[i].lv.Px());
	py.push_back(jet[i].lv.Py());
	pz.push_back(jet[i].lv.Pz());
	 E.push_back(jet[i].lv.E());
	jsel.push_back(i);
  }
		
  if (px.size()<2) return -999;

  // get hemispheres (seed 2: max inv mass, association method: default 3 = minimal lund distance)
  Hemisphere* hemi = new Hemisphere(px, py, pz, E, hemi_seed, hemi_association);
  hemi->SetDebug(0);
  if(MaxDR > 0) hemi->RejectISRDRmax(MaxDR);
  vector<int> grouping = hemi->getGrouping();

  TLorentzVector pseudojet1(0.,0.,0.,0.);
  TLorentzVector pseudojet2(0.,0.,0.,0.);
  vector<int>    jused; // contains indices of jets used for pseudojets
  for(int i=0; i<px.size(); ++i){
	if(grouping[i]==1){
		jused     .push_back(jsel[i]);
		pseudojet1.SetPx(pseudojet1.Px() + px[i]);
		pseudojet1.SetPy(pseudojet1.Py() + py[i]);
		pseudojet1.SetPz(pseudojet1.Pz() + pz[i]);
		pseudojet1.SetE( pseudojet1.E()  + E[i]);	
	}else if(grouping[i] == 2){
		jused     .push_back(jsel[i]);
		pseudojet2.SetPx(pseudojet2.Px() + px[i]);
		pseudojet2.SetPy(pseudojet2.Py() + py[i]);
		pseudojet2.SetPz(pseudojet2.Pz() + pz[i]);
		pseudojet2.SetE( pseudojet2.E()  + E[i]);
	}
  }
  delete hemi;
 
  for(int i=0; i<NJets; ++i){
	bool used(false);  
  	for(int j=0; j<jused.size(); ++j){
		if(jused[j] == i) used=true; 
	}
	TLorentzVector unused_jet(0, 0, 0, 0);
	if(! used && jet[0].lv.Pt()>jet_for_met_threshold && met==3) unused_jet.SetXYZM(jet[i].lv.Px(), jet[i].lv.Py(), 0, 0);
        MET += unused_jet; // adding all jets that were rejected by hemisphere algo to MET	
  }
  return CalcMT2(0, massive, pseudojet1, pseudojet2, MET); 
}

Double_t MT2tree::GetMT2Hemi(double testmass, bool massive, int PFJID, double minJPt, int hemi_association, int met) {
  TLorentzVector MET(0., 0., 0., 0.);
  if(met==1)      MET = pfmet[0];
  else if(met==2) MET = MHTloose[0];
  else            return -999;

  vector<float> px, py, pz, E;
  for(int i=0; i<NJets; ++i){
	if(PFJID==1  && jet[i].isPFIDLoose ==false) continue;
	if(PFJID==2  && jet[i].isPFIDMedium==false) continue;
	if(PFJID==3  && jet[i].isPFIDTight ==false) continue;
	if(jet[i].lv.Pt() < minJPt )                continue;
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
 
  return CalcMT2(testmass, massive, pseudojet1, pseudojet2, MET); 
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

Double_t MT2tree::GetMCT(bool massive, int met) { // FIXME!!
  TLorentzVector MET(0., 0., 0., 0.);
  if(met==1)      MET = pfmet[0];
  else if(met==2) MET = MHTloose[0];
  else            return -999;

  if (NJetsIDLoose<2)
    return -999;

  TLorentzVector DTM(0,0,0,0);
  TLorentzVector p1, p2;
//  p1.SetXYZM(pseudoJets[0].Px(), pseudoJets[0].Py(), pseudoJets[0].Pz(), massive ? pseudoJets[0].M() : 0);
//  p2.SetXYZM(pseudoJets[1].Px(), pseudoJets[1].Py(), pseudoJets[1].Pz(), massive ? pseudoJets[1].M() : 0);
  TVector2 pmiss;
  pmiss.Set(MET.Px(), MET.Py());

  TMctLib *mct = new TMctLib();
  Double_t MCT =  mct -> mctcorr(p1, p2, DTM, pmiss, 7000, 0.);
  delete mct;
  return -999.99;
}

ClassImp(MT2Misc)
ClassImp(MT2Jet)
ClassImp(MT2Elec)
ClassImp(MT2Muon)
ClassImp(MT2Hemi)
ClassImp(MT2GenLept)
ClassImp(MT2tree)
