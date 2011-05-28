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
  HBHENoiseFlag           =  0;
  CrazyHCAL               =  0;
  isData                  =  0;
  Run                     = -1;	  
  Event		  	  = -1;	  
  LumiSection		  = -1;	  
  LeptConfig		  = -1;	  
  PassJetID               = -1;
  PassJetID20             = -1;
  Jet0Pass                = -1;
  Jet1Pass                = -1;
  MT2                     = -99999.99;
  MT2all                  = -99999.99;
  MT2leading              = -99999.99;
  MT2noISR                = -99999.99;
  MCT                     = -99999.99;
  AlphaT                  = -99999.99;
  MET                     = -99999.99;
  METPhi                  = -99999.99;
  LeadingJPt              = -99999.99;
  SecondJPt               = -99999.99;
  Vectorsumpt		  = -99999.99;
  VectorsumptAll 	  = -99999.99;
  PFMETsign		  = -99999.99;
  HT			  = -99999.99;
  caloHT30       	  = -99999.99;
  caloHT40       	  = -99999.99;
  caloHT50       	  = -99999.99;
  caloMHT20       	  = -99999.99;
  caloMHT30       	  = -99999.99;
  caloMHT40       	  = -99999.99;
  DPhiMhtMpt              = -99999.99;
  MinMetJetDPhi           = -99999.99;
  EcalDeadCellBEFlag      = -1;
  NECALGapClusters        = -1;
  for(int i=0; i<50; ++i){
    EcalGapClusterSize[i] = -1;
    EcalGapBE[i]          = -99999.99;    
  }

}

// ------------------------------------------------------------
// MT2PileUp
MT2PileUp::MT2PileUp(){
	Reset();
}

MT2PileUp::~MT2PileUp(){
}

void MT2PileUp::Reset(){
	PUnumInt  = -999;
	PtHat     = -999.99;;
	Weight    = -999.99;
  	NVertices = -1;
}

// ------------------------------------------------------
// MT2trigger
MT2Trigger::MT2Trigger(){
	Reset();
}
MT2Trigger::~MT2Trigger(){
}

void MT2Trigger::Reset(){
	
	// HT
	HLT_HT150_v2            = false;
	HLT_HT150_v3            = false;
	HLT_HT160_v2            = false;
	HLT_HT200_v2            = false;
	HLT_HT200_v3            = false;
	HLT_HT240_v2            = false;
	HLT_HT250_v2            = false;
	HLT_HT250_v3            = false;
	HLT_HT260_v2            = false;
	HLT_HT300_v2            = false;
	HLT_HT300_v3            = false;
	HLT_HT300_v4            = false;
	HLT_HT350_v2            = false;
	HLT_HT350_v3            = false;
	HLT_HT360_v2            = false;
	HLT_HT400_v2            = false;
	HLT_HT400_v3            = false;
	HLT_HT440_v2            = false;
	HLT_HT450_v2            = false;
	HLT_HT450_v3            = false;
	HLT_HT500_v2            = false;
	HLT_HT500_v3            = false;
	HLT_HT550_v2            = false;
	HLT_HT550_v3            = false;
	// HT_MHT
	HLT_HT250_MHT60_v2      = false;
	HLT_HT250_MHT60_v3      = false;
	HLT_HT260_MHT60_v2      = false;
	HLT_HT300_MHT75_v4      = false;
	// QuadJet
	HLT_QuadJet50_BTagIP_v1 = false;
	HLT_QuadJet50_Jet40_v1  = false;
	// Muons
	HLT_DoubleMu3_HT160_v2  = false;
	HLT_Mu8_Jet40_v2        = false;
	HLT_DoubleMu3_v3        = false;
}


// MT2Znunu ------------------------------------
MT2Znunu::MT2Znunu(){
	Reset();
}

MT2Znunu::~MT2Znunu(){
}

void MT2Znunu::Reset(){
	NJetsToRemoveEle          = -999;
	NJetsToRemoveMuo          = -999;
	NJetsIDLoose_matched      = -999;
	PassJetID_matched         = -999;
	Jet1Pass_matched          = -999;
	Jet0Pass_matched          = -999;
	LeadingJPt_matched        = -99999.99;
	SecondJPt_matched         = -99999.99;
	HTmatched                 = -99999.99;
	caloMHT30_matched         = -99999.99;
	caloMHT30_matchedReco     = -99999.99;
	caloHT50_matched          = -99999.99;
	caloHT50_matchedReco      = -99999.99;
	GenZmumu_mll              = -99999.99;
	GenZmumu_mll_acc          = -99999.99;
	GenZee_mll                = -99999.99;
	GenZee_mll_acc            = -99999.99;
	GenZnunu_e_mll            = -99999.99;
	GenZnunu_e_mll_acc        = -99999.99;
	GenZnunu_mu_mll           = -99999.99;
	GenZnunu_mu_mll_acc       = -99999.99;
	GenZnunu_tau_mll          = -99999.99;
	GenZnunu_tau_mll_acc      = -99999.99;
	RecoOSee_mll              = -99999.99;
	RecoOSmumu_mll            = -99999.99;
	METplusLeptsPt            = -99999.99;
	METplusLeptsPtReco        = -99999.99;
	MinMetplusLeptJetDPhi     = -99999.99;
	MinMetplusLeptJetDPhiReco = -99999.99;
	Vectorsumpt_matched       = -99999.99;
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

  isPATPFIDLoose= 0; // PAT PFJetIDLoose regardless of eta and pt
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
  
  isTau         = 0;  // starting from ntuple V02-01-01: this has to be 0! 
  isTauMatch    = 0;  // tell you if the jet is matched to a tau.
  TauDR         = -99999.99;
  TauDPt        = -99999.99;
  NTauMatch     = 0;
}

void MT2Jet::SetLV(const TLorentzVector v) {
  lv = v;
}

Bool_t MT2Jet::IsGoodPFJet(double minJPt, double maxJEta, int PFJID) {

  double pt = lv.Pt();
  double eta = lv.Eta();
  if ( pt < minJPt || fabs(eta) > maxJEta )     return false;
  if ( isTau )                                  return true;  // for now every tau passes the "ID".
  
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
  MT       = -9999.99;
}

// MT2Muon -----------------------------------
MT2Muon::MT2Muon(){
  Reset();
}

MT2Muon::~MT2Muon(){
}

void MT2Muon::Reset() {
  lv.SetPxPyPzE(0, 0, 0, 0);
  MT            = -9999.99;
  Iso           = -9999.99;
  Charge        = -999;
  NMatches      = -999;
  PtErr         = -999.99;
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
  MT            = -9999.99;
  Iso           = -9999.99;
  Charge        = -999;
  ID95          = -999;
  ID90          = -999;
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
  NJets            = 0;
  NTaus            = 0;
  NJetsIDLoose     = 0;
  NJetsIDMedium    = 0;
  NJetsIDTight     = 0;
  NJetsAcc         = 0;
  NEles            = 0;
  NMuons           = 0;
  NGenLepts        = 0;

  misc.Reset();
  Znunu.Reset();
  pileUp.Reset();
  trigger.Reset();

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
  genmet    [0].SetPxPyPzE(0., 0., 0., 0.);
  MPT       [0].SetPxPyPzE(0., 0., 0., 0.);
  MHTloose  [0].SetPxPyPzE(0., 0., 0., 0.);
  MHT       [0].SetPxPyPzE(0., 0., 0., 0.);
}

void MT2tree::SetNJets(int n) {
  NJets = n;
}

void MT2tree::SetNJetsAcc(int n) {
  NJetsAcc = n;
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

void MT2tree::SetNTaus(int n) {
  NTaus = n;
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
  else if(met==2) MET = GetMHTlv(PFJID, minJPt, maxJEta);
  else if(met==3) {
	  double mass = GetDiLeptonInvMass(0,1,0,5.0,1);
	  if(mass > 71 && mass <111) MET = GetMETPlusLeptsLV(1) + pfmet[0] ;
	  else return -888.;
  }
  else         return -999;


  int index = MinMetJetDPhiIndex(PFJID, minJPt, maxJEta, met);
  if(index >=0) {
	  return TMath::Abs(jet[index].lv.DeltaPhi(MET));
  } else  return -999.99;
}

Int_t MT2tree::MinMetJetDPhiIndex(int PFJID, double minJPt, double maxJEta, int met) {
// Attention: electrons and muons are not considered for minDPhi
  TLorentzVector MET(0., 0., 0., 0.);
  if(met==1)      MET = pfmet[0];
  else if(met==2) MET = GetMHTlv(PFJID, minJPt, maxJEta);
  else if(met==3) {
	  double mass = GetDiLeptonInvMass(0,1,0,5.0,1);
	  if(mass > 71 && mass <111) MET = GetMETPlusLeptsLV(1) + pfmet[0] ;
	  else return -888.;
  }
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

Double_t MT2tree::MaxMetJetDPhi(int PFJID, double minJPt, double maxJEta, int met) {
// Attention: electrons and muons are not considered for minDPhi
  TLorentzVector MET(0., 0., 0., 0.);
  if(met==1)      MET = pfmet[0];
  else if(met==2) MET = MHTloose[0];
  else            return -999;

  int index = MaxMetJetDPhiIndex(PFJID, minJPt, maxJEta, met);
  return TMath::Abs(jet[index].lv.DeltaPhi(MET));
}

Int_t MT2tree::MaxMetJetDPhiIndex(int PFJID, double minJPt, double maxJEta, int met) {
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
  Double_t maxDPhi=0.;
  Int_t    imax;
  for (int i=0; i<indices.size(); i++){
    Double_t dphi = TMath::Abs(jet[indices[i]].lv.DeltaPhi(MET));
    if(dphi>maxDPhi)  {
      maxDPhi=dphi;
      imax = indices[i];
    }
  }
  return imax;
}

Bool_t MT2tree::PassMinMetJetDPhi03(){
	if( NJetsIDLoose < 3)                              return true;
	if( NJetsIDLoose >=3  && misc.MinMetJetDPhi > 0.3) return true;
	return false;
}

Int_t MT2tree::GetNjets(double minJPt, double maxJEta, int PFJID){
  int njets=0;
  for(int i=0; i<NJets; ++i){
	if(jet[i].IsGoodPFJet(minJPt,maxJEta,PFJID)==false) continue;
	njets++;
  }
  return njets;
}

Int_t MT2tree::GetNBtags (int algo, double value, double minJPt, double maxJEta, int PFJID){  // algo - 0:TCHE, 1:TCHP, 2:SSVHE, 3:SSVHP
  int nbjets=0;
  for(int i=0; i<NJets; ++i){
    if(jet[i].IsGoodPFJet(minJPt,maxJEta,PFJID)==false || jet[i].isTau ) continue; // FIXME: taus don't have JID and BTAG info
    switch(algo){
    case 0: 
      if( jet[i].bTagProbTCHE < value ) continue;
      break;
    case 1: 
      if( jet[i].bTagProbTCHP < value ) continue;
      break;
    case 2: 
      if( jet[i].bTagProbSSVHE < value ) continue;
      break;
    case 3: 
      if( jet[i].bTagProbSSVHP < value ) continue;
      break;
    default:
      continue;
      break;
    }
    nbjets++; 
  }
  return nbjets;
}

Double_t MT2tree::JetPt(int ijet, int PFJID, double minJPt, double maxJEta) {
  int index = GetJetIndex(ijet, PFJID,minJPt,maxJEta);
  if ( index < 0 )       return -999.;
  return jet[index].lv.Pt();
}

Int_t MT2tree::GetJetIndex(int ijet, int PFJID, double minJPt, double maxJEta) {
  std::vector<int> indices;
  for(int i = 0; i<NJets; ++i){
    if(jet[i].IsGoodPFJet(minJPt,maxJEta,PFJID)==false) continue;
    indices.push_back(i);
  }
  if (ijet>=indices.size())  return -9;
  return indices[ijet];
}

Int_t MT2tree::GetJetIndexByEta(int ijet, int PFJID, double minJPt, double maxJEta) {
  std::vector<int> indices, indx;
  std::vector<double> etas;
  for(int i = 0; i<NJets; ++i){
    if(jet[i].IsGoodPFJet(minJPt,maxJEta,PFJID)==false) continue;
    indices.push_back(i);
    etas.push_back(fabs(jet[i].lv.Eta()));
  }
  if (ijet>=indices.size())  return -9;
  indx = Util::VSort(indices,etas,true);
  return indx[ijet];
}

Double_t MT2tree::GetHT(int PFJID, double minJPt, double maxJEta){
  Double_t ht=0;
  for(int i = 0; i<NJets; ++i){
    if(jet[i].IsGoodPFJet(minJPt,maxJEta,PFJID)==false) continue;
    ht+=jet[i].lv.Pt();
  }
  return ht;
}

TLorentzVector MT2tree::GetMHTlv(int PFJID, double minJPt, double maxJEta, bool inclLepts){
  TLorentzVector mht(0,0,0,0);
  TLorentzVector j(0,0,0,0);
  for(int i = 0; i<NJets; ++i){
    if(jet[i].IsGoodPFJet(minJPt,maxJEta,PFJID)==false) continue;
    j.SetPtEtaPhiM(jet[i].lv.Pt(),0,jet[i].lv.Phi(),0);
    mht-=j;
  }
  if(!inclLepts) return mht;
  // add leptons
  for(int i=0; i<NEles; ++i){
    j.SetPtEtaPhiM(ele[i].lv.Pt(),0, ele[i].lv.Phi(), 0);
    mht-=j;
  }
  for(int i=0; i<NMuons; ++i){
    j.SetPtEtaPhiM(muo[i].lv.Pt(),0, muo[i].lv.Phi(), 0);
    mht-=j;
  }
  return mht;
}

Double_t MT2tree::GetMHT(int PFJID, double minJPt, double maxJEta, bool inclLepts){
  TLorentzVector mht = GetMHTlv(PFJID,minJPt,maxJEta, inclLepts);
  return mht.Pt();
}

Double_t MT2tree::GetMHTPhi(int PFJID, double minJPt, double maxJEta, bool inclLepts){
  TLorentzVector mht = GetMHTlv(PFJID,minJPt,maxJEta, inclLepts);
  return mht.Phi();
}

Double_t MT2tree::GetMHTminusMET(int PFJID, double minJPt, double maxJEta, bool inclLepts){
  TLorentzVector mht = GetMHTlv(PFJID,minJPt,maxJEta, inclLepts);
  return (mht-pfmet[0]).Pt();
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
  //if(index!=2) return -999;
  //cout << "leadingjets " <<leadingJets[0].Pt() << " " << leadingJets[1].Pt() << endl;
  return CalcMT2(testmass, massive, leadingJets[0], leadingJets[1], MET);

}

Double_t MT2tree::GetMT2HemiNoISR(bool massive, int hemi_seed, int hemi_association, float MaxDR, int met) {
  // hemi_seed   = 2: max invariant mass, 4: two leading jet
  // hemi_assoc  = 2: min invariant mass, 3: minial lund distance, 
  // MaxDR: if >0 RejectISRDRmax(MaxDR) is called for hemi            
  // met ==3 option adds jets > jet_for_met_threshold GeV rejected by the chosen hemispheres to the MET
  // met ==4 option takes sets \vec{MET} = -\vec{pj1} - \vec{pj2}, i.e. UTM is zero. 
  float jet_for_met_threshold=20;

  TLorentzVector MET(0., 0., 0., 0.);
  if(met==1 || met==3)  MET = pfmet[0];
  else if(met==2)       MET = MHTloose[0];
  else if(met==4)       MET.SetPxPyPzE(0,0,0,0);
  else return -999;
 

  vector<float> px, py, pz, E;
  vector<int>   jsel; // contains indices of all jets fed into hemisphere algo
  for(int i=0; i<NJets; ++i){
	if(jet[i].IsGoodPFJet(20, 2.4, 1)) continue;
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
  if(met==4) MET = -pseudojet1 - pseudojet2;
  return CalcMT2(0, massive, pseudojet1, pseudojet2, MET); 
}

Double_t MT2tree::GetMT2Hemi(double testmass, bool massive, int PFJID, double minJPt, double maxJEta, int hemi_association, int met) {
  TLorentzVector MET(0., 0., 0., 0.);
  if(met==1)      MET = pfmet[0];
  else if(met==2) MET = MHTloose[0];
  else if(met==3) MET = pfmet[0]; //plus OS dileptons
  else if(met==4) MET.SetPxPyPzE(0,0,0,0);
  else            return -999;

  if( met ==3 ){
  // adding OS dilepton LV to MET
	double dilept_invmass= GetDiLeptonInvMass(0,1,0,5.0,true);
	if(dilept_invmass < 111 && dilept_invmass > 71){
		for(int i=0; i<NEles; ++i){
			MET = MET + ele[i].lv;
		}
		for(int i=0; i<NMuons; ++i){
			MET = MET + muo[i].lv;
		}
	} else {return -111;}
  }

  vector<float> px, py, pz, E;
  for(int i=0; i<NJets; ++i){
	if(jet[i].IsGoodPFJet(minJPt, maxJEta, PFJID) ==false) continue;
  	px.push_back(jet[i].lv.Px());
	py.push_back(jet[i].lv.Py());
	pz.push_back(jet[i].lv.Pz());
	E .push_back(jet[i].lv.E ());
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

  if(met==4) {MET = -pseudojet1 - pseudojet2;} 
  if(MET.Pt()<30) {return -222;}
  return CalcMT2(testmass, massive, pseudojet1, pseudojet2, MET); 
}

Double_t MT2tree::GetMT2HemiMinDHT(double testmass, bool massive, int PFJID, double minJPt, double maxJEta, int met) {
  TLorentzVector MET(0., 0., 0., 0.);
  if(met==1)      MET = pfmet[0];
  else if(met==2) MET = MHTloose[0];
  else if(met==3) MET = pfmet[0]; //plus OS dileptons
  else if(met==4) MET.SetPxPyPzE(0,0,0,0);
  else            return -999;

  if( met ==3 ){
  // adding OS dilepton LV to MET
	double dilept_invmass= GetDiLeptonInvMass(0,1,0,5.0,true);
	if(dilept_invmass < 111 && dilept_invmass > 71){
		for(int i=0; i<NEles; ++i){
			MET = MET + ele[i].lv;
		}
		for(int i=0; i<NMuons; ++i){
			MET = MET + muo[i].lv;
		}
	} else {return -111;}
  }

  
  vector<TLorentzVector> p4s;
  for(int i=0; i<NJets; ++i){
	if(jet[i].IsGoodPFJet(minJPt, maxJEta, PFJID) ==false) continue;
	p4s.push_back(jet[i].lv);
  }

  if(p4s.size() < 2) return -200;

  std::vector<std::vector<double> > ht( 1<<(p4s.size()-1) , std::vector<double>( 2, 0.) ); // initializiert einen vector ht der size 1<<(p4s.size()-1), wobei jeder eintrag ein vector (0,0) ist

  std::vector<std::vector<double> > px( 1<<(p4s.size()-1) , std::vector<double>( 2, 0.) ); 
  std::vector<std::vector<double> > py( 1<<(p4s.size()-1) , std::vector<double>( 2, 0.) ); 
  std::vector<std::vector<double> > pz( 1<<(p4s.size()-1) , std::vector<double>( 2, 0.) ); 
  std::vector<std::vector<double> > E ( 1<<(p4s.size()-1) , std::vector<double>( 2, 0.) ); 

  for(unsigned i=0; i < ht.size(); i++) {                                             
  for(unsigned j=0; j < p4s.size(); j++) {
		ht [i] [(i/(1<<j))%2] += p4s[j].Pt();
		px [i] [(i/(1<<j))%2] += p4s[j].Px();
		py [i] [(i/(1<<j))%2] += p4s[j].Py();
		pz [i] [(i/(1<<j))%2] += p4s[j].Pz();
		E  [i] [(i/(1<<j))%2] += p4s[j].E();
  }  
  }
  std::vector<double> deltaHT; for(unsigned i=0; i<ht.size(); i++) deltaHT.push_back(fabs(ht[i][0]-ht[i][1]));
  const double mDHT = *(std::min_element( deltaHT.begin(), deltaHT.end() ));
  int pos=distance(deltaHT.begin(), min_element(deltaHT.begin(), deltaHT.end()));
  TLorentzVector pj1, pj2;
  pj1.SetPxPyPzE(px[pos][0], py[pos][0], pz[pos][0], E[pos][0]);
  pj2.SetPxPyPzE(px[pos][1], py[pos][1], pz[pos][1], E[pos][1]);

  if(met==4) {MET = -pj1 - pj2;} 
  if(MET.Pt()<30) {return -222;}
  return CalcMT2(testmass, massive, pj1, pj2, MET); 
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

Double_t MT2tree::SimpleMT2(bool pseudo){
  if (!pseudo){
    if(NJets < 2) return -999.99;
    return sqrt(2*jet[0].lv.Pt()*jet[1].lv.Pt()*(1+TMath::Cos(jet[0].lv.DeltaPhi(jet[1].lv))));
  }
  else{
    if(NJetsAcc < 2) return -999.99;
    return sqrt(2*hemi[0].lv1.Pt()*hemi[0].lv2.Pt()*(1+TMath::Cos(hemi[0].lv1.DeltaPhi(hemi[0].lv2))));
  }
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

Double_t MT2tree::GetSqrtS(double testmass, bool massive, int PFJID, double minJPt, double maxJEta,int met){
  TLorentzVector MET(0., 0., 0., 0.);
  if(met==1)      MET = pfmet[0];
  else if(met==2) MET = MHTloose[0];
  else            return -999;

  TLorentzVector obj(0,0,0,0);
  TLorentzVector vis(0,0,0,0);
  for(int i = 0; i<NJets; ++i){
    if(jet[i].IsGoodPFJet(minJPt,maxJEta,PFJID)==false) continue;
    obj.SetPtEtaPhiM(jet[i].lv.Pt(),jet[i].lv.Eta(),jet[i].lv.Phi(),massive ? jet[i].lv.M() : 0);
    vis+=obj;
  }
  for(int i = 0; i<NEles; ++i){
    obj.SetPtEtaPhiM(ele[i].lv.Pt(),ele[i].lv.Eta(),ele[i].lv.Phi(),ele[i].lv.M());
    vis+=obj;
  }
  for(int i = 0; i<NMuons; ++i){
    obj.SetPtEtaPhiM(muo[i].lv.Pt(),muo[i].lv.Eta(),muo[i].lv.Phi(),muo[i].lv.M());
    vis+=obj;
  }
  return sqrt(TMath::Power(vis.M(),2)+TMath::Power(vis.Pt(),2)) + sqrt(TMath::Power(2*testmass,2)+TMath::Power(MET.Pt(),2));
}

Double_t MT2tree::GetMaxHemiMass(int hemi_index){
	return TMath::Max(hemi[hemi_index].lv1.M(),hemi[hemi_index].lv2.M());
}

Double_t MT2tree::GetPseudoJetsdPhi(int hemi_seed, int hemi_association, int PFJID, double minJPt, double maxJEta){

  vector<float> px, py, pz, E;
  for(int i=0; i<NJets; ++i){
	if(jet[i].IsGoodPFJet(minJPt, maxJEta, PFJID) ==false) continue;
  	px.push_back(jet[i].lv.Px());
	py.push_back(jet[i].lv.Py());
	pz.push_back(jet[i].lv.Pz());
	 E.push_back(jet[i].lv.E());
  }
		
  if (px.size()<2) return -999;

  // get hemispheres (seed 2: max inv mass, association method: default 3 = minimal lund distance)
  Hemisphere* hemi = new Hemisphere(px, py, pz, E, hemi_seed, hemi_association);
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
  return Util::DeltaPhi(pseudojet1.Phi(), pseudojet2.Phi());
}

Double_t MT2tree::GetPseudoJetsdPhiMinDHT(int PFJID, double minJPt, double maxJEta){
	// hemispheres minimizing deltaHT
	vector<TLorentzVector> p4s;
	for(int i=0; i<NJets; ++i){
		if(jet[i].IsGoodPFJet(minJPt, maxJEta, PFJID) ==false) continue;
		p4s.push_back(jet[i].lv);
	}

	if(p4s.size() < 2) return -200;

  	std::vector<std::vector<double> > ht( 1<<(p4s.size()-1) , std::vector<double>( 2, 0.) ); // initializiert einen vector ht der size 1<<(p4s.size()-1), wobei jeder eintrag ein vector (0,0) ist
       	
  	std::vector<std::vector<double> > px( 1<<(p4s.size()-1) , std::vector<double>( 2, 0.) ); 
  	std::vector<std::vector<double> > py( 1<<(p4s.size()-1) , std::vector<double>( 2, 0.) ); 
  	std::vector<std::vector<double> > pz( 1<<(p4s.size()-1) , std::vector<double>( 2, 0.) ); 
  	std::vector<std::vector<double> > E ( 1<<(p4s.size()-1) , std::vector<double>( 2, 0.) ); 

  	for(unsigned i=0; i < ht.size(); i++) {                                             
    	for(unsigned j=0; j < p4s.size(); j++) {
			ht [i] [(i/(1<<j))%2] += p4s[j].Pt();
			px [i] [(i/(1<<j))%2] += p4s[j].Px();
			py [i] [(i/(1<<j))%2] += p4s[j].Py();
			pz [i] [(i/(1<<j))%2] += p4s[j].Pz();
			E  [i] [(i/(1<<j))%2] += p4s[j].E();
    	}
  	}
  	std::vector<double> deltaHT; for(unsigned i=0; i<ht.size(); i++) deltaHT.push_back(fabs(ht[i][0]-ht[i][1]));
	const double mDHT = *(std::min_element( deltaHT.begin(), deltaHT.end() ));
	int pos=distance(deltaHT.begin(), min_element(deltaHT.begin(), deltaHT.end()));
	TLorentzVector pj1, pj2;
	pj1.SetPxPyPzE(px[pos][0], py[pos][0], pz[pos][0], E[pos][0]);
	pj2.SetPxPyPzE(px[pos][1], py[pos][1], pz[pos][1], E[pos][1]);

	return Util::DeltaPhi(pj1.Phi(), pj2.Phi());
}

Double_t MT2tree::GetPseudoJetMetDPhi(int hemi_index, int pj, int whichmet, double met){
	if(whichmet==1 && pfmet[0].Pt() < met)                                    return -777;
	if(whichmet==2 && (hemi[hemi_index].lv1+hemi[hemi_index].lv2).Pt() < met) return -777;
	if(pj ==1){
		if(whichmet==1)        return TMath::Abs(hemi[hemi_index].lv1.DeltaPhi(pfmet[0]));
		else if (whichmet ==2) return TMath::Abs(hemi[hemi_index].lv1.DeltaPhi(-hemi[hemi_index].lv1 -hemi[hemi_index].lv2));
		else                   return -111;
	}else if (pj ==2){
		if(whichmet==1)        return TMath::Abs(hemi[hemi_index].lv2.DeltaPhi(pfmet[0]));
		else if (whichmet ==2) return TMath::Abs(hemi[hemi_index].lv2.DeltaPhi(-hemi[hemi_index].lv1 -hemi[hemi_index].lv2));
		else                   return -111;       
	}else {return -999.99;}
}

Double_t MT2tree::GenOSDiLeptonInvMass(unsigned int pid, unsigned int mother, double pt, double eta){
	 // returns true if there are 
	 //  - two particles with abs(ID)=pid
	 //  - they have opposite charge
	 //  - they come from "mother"
	 //  - they have Pt() > pt and |Eta()| < eta
	 //  - their inv mass is in (lower_Minv, upper_Minv)
	vector<int> indices;	
	for(int i=0; i<NGenLepts; ++i){
		if(abs(genlept[i].ID)       !=pid   ) continue;
		if(abs(genlept[i].MID)      !=mother) continue;
		if(    genlept[i].lv.Pt()   < pt    ) continue;
		if(fabs(genlept[i].lv.Eta())> eta   ) continue;
		indices.push_back(i);
	}
	if(indices.size()!=2)                                        return -100;
	if( (genlept[indices[0]].ID) * (genlept[indices[1]].ID) >=0 && pid !=12 && pid != 14 && pid !=16) return -150;
	
	double Mass= (genlept[indices[0]].lv + genlept[indices[1]].lv).M();
	if(Mass<0) Mass = -Mass;
	return Mass;

}


Bool_t   MT2tree::IsGenOSDiLepton(unsigned int pid, unsigned int mother, double pt, double eta, double lower_mass, double upper_mass){
	 double Mass = GenOSDiLeptonInvMass(pid, mother, pt, eta);
	 if(Mass < 0 ) return false;
	 if(Mass < lower_mass || Mass > upper_mass) return false;
	 return true;
}


Double_t MT2tree::GetDiLeptonInvMass(int same_sign, int same_flavour, int flavour, double pt, bool exclDiLept){
	// flavour == 0 : don't care if el or mu
	// flavour == 1 : only electrons
	// flavour == 2 : only muons
	
	struct lepton {TLorentzVector lv; int charge; string flavour;} lepts[NEles + NMuons];	
	int counter = 0;
	for(int i=0; i<NEles; ++i){
		lepts[counter].lv=ele[i].lv; lepts[counter].charge =ele[i].Charge; lepts[counter].flavour ="ele";
		counter++;
	}
	for(int i=0; i<NMuons; ++i){
		lepts[counter].lv=muo[i].lv; lepts[counter].charge =muo[i].Charge; lepts[counter].flavour ="muo"; 
		counter++;
	}
	
	int index1=-1, index2=-1;
	if( (NEles + NMuons) <= 1)                                               return -10;
	if( (NEles + NMuons) >  2 && exclDiLept)                                 return -100;
	else if( (NEles + NMuons) >  2 && !exclDiLept){
		double pt1=0, pt2=0;
		for(int i=0; i<(NEles + NMuons); ++i){
			if(lepts[i].lv.Pt()>pt1){
				pt2 =pt1;
				index2=index1;
				pt1 =lepts[i].lv.Pt();
				index1=i;
			} else if(lepts[i].lv.Pt()>pt2){
				index2=i;
				pt2   =lepts[i].lv.Pt();
			}
		}	
	
	}else if ((NEles + NMuons)==2) {
		index1=0; index2=1; 
	}

	if(  lepts[index1].charge * lepts[index2].charge ==  1 && same_sign ==0)             return -200;
	if(  lepts[index1].charge * lepts[index2].charge == -1 && same_sign ==1)             return -300;
	if(  lepts[index1].lv.Pt() < pt || lepts[index2].lv.Pt() < pt )                      return -400;
	if(  lepts[index1].flavour != lepts[index2].flavour && same_flavour == 1 )           return -500;
	if(  lepts[index1].flavour == lepts[index2].flavour && same_flavour == 0 )           return -600;
	if( (lepts[index1].flavour =="muo" || lepts[index2].flavour=="muo") && flavour == 1) return -910;
	if( (lepts[index1].flavour =="ele" || lepts[index2].flavour=="ele") && flavour == 2) return -920;

	TLorentzVector sum = lepts[index1].lv + lepts[index2].lv;
	double mass = sum.M();
	if(mass < 0) return -mass;
	else return mass;
}

Bool_t MT2tree::IsDiLeptonMll(int same_sign, int same_flavour, int flavour, double pt, bool exclDiLept, double lower_mass, double upper_mass){
	double mass = GetDiLeptonInvMass(same_sign, same_flavour, flavour, pt, exclDiLept);
	if(mass < 0) return false;
	if(mass < lower_mass ||  mass > upper_mass) return false;
	return true;
}


TLorentzVector MT2tree::GetMETPlusLeptsLV(int OSDiLeptFromZ){
	TLorentzVector lv = pfmet[0];
	if(OSDiLeptFromZ ==1) {
		double mass = GetDiLeptonInvMass(0,1,0,5.0,1);
		if(mass < 0 )               return lv;
		if(mass < 71 || mass > 111) return lv; 
	}
	for(int i=0; i<NEles; ++i){
		lv+=ele[i].lv;
	}
	for(int i=0; i<NMuons; ++i){
		lv+=muo[i].lv;
	}
	return lv;
}

Double_t MT2tree::GetMETPlusLepts(int OSDiLeptFromZ){
	if(OSDiLeptFromZ ==1) {
		double mass = GetDiLeptonInvMass(0,1,0,5.0,1);
		if(mass < 0 )               return -2000;
		if(mass < 71 || mass > 111) return -1000; 
	}
	TLorentzVector lv = pfmet[0];
	for(int i=0; i<NEles; ++i){
		lv+=ele[i].lv;
	}
	for(int i=0; i<NMuons; ++i){
		lv+=muo[i].lv;
	}
	return lv.Pt();
}

Double_t MT2tree::GetDiLeptonPt(int same_sign, int same_flavour, int flavour, double pt, double lower_mass, double upper_mass){
	double mass = GetDiLeptonInvMass(same_sign, same_flavour, flavour, pt, 1);
	if(mass < 0 )                              return -2000;
	if(mass < lower_mass || mass > upper_mass) return -1000; 
	TLorentzVector lv(0., 0., 0., 0.);
	for(int i=0; i<NEles; ++i){
		lv+=ele[i].lv;
	}
	for(int i=0; i<NMuons; ++i){
		lv+=muo[i].lv;
	}
	return lv.Pt();
}

Double_t MT2tree::GetMETPlusGenLepts(int met, int RemoveOSSFDiLepts, int require_cuts,  unsigned int pid, unsigned int mother, double pt, double eta, double lower_mass, double upper_mass ){

	TLorentzVector lv;
        if(met==1) lv= genmet[0];
        if(met==0) lv= pfmet[0];


	vector<int> indices;	
	if(RemoveOSSFDiLepts ==1 || require_cuts ==1 ) {
		for(int i=0; i<NGenLepts; ++i){
			if(pid==1113) {if(abs(genlept[i].ID) !=11 && abs(genlept[i].ID)!=13  ) continue;}
			if(pid!=1113) {if(abs(genlept[i].ID) !=pid   ) continue;}
			if(abs(genlept[i].MID)               !=mother) continue;
			if(    genlept[i].lv.Pt()            < pt    ) continue;
			if(fabs(genlept[i].lv.Eta())         > eta   ) continue;
			indices.push_back(i);
		}
		if(indices.size()!=2)                                         return -100;
		if( (genlept[indices[0]].ID) * (genlept[indices[1]].ID) >=0 ) return -150;
		
		double Mass= (genlept[indices[0]].lv + genlept[indices[1]].lv).M();
		if(Mass<0) Mass = -Mass;
		if(Mass > upper_mass || Mass < lower_mass)                    return -200;
		if(RemoveOSSFDiLepts ==1) {
			lv+=genlept[indices[0]].lv;
			lv+=genlept[indices[1]].lv;
		}
	}
	
	return lv.Pt();
}

Int_t   MT2tree::GetGenLeptIndex(int which, int pid, int mother, double pt, double eta){
	vector<int>    indices;
	vector<double> pts, etas;
	for(int i=0; i<NGenLepts; ++i){
		if(pid==1113)        {if(abs(genlept[i].ID) !=11 && abs(genlept[i].ID)!=13  ) continue;}
		else if(pid==121416) {if(abs(genlept[i].ID) !=12 && abs(genlept[i].ID)!=14 && abs(genlept[i].ID)!=16) continue;}
		else if(abs(genlept[i].ID) !=pid   ) continue;
		if(abs(genlept[i].MID)               !=mother) continue;
		if(    genlept[i].lv.Pt()            < pt    ) continue;
		if(fabs(genlept[i].lv.Eta())         > eta   ) continue;
		indices.push_back(i);
		pts.push_back(genlept[i].lv.Pt());
	}
	if(indices.size() < 1 || which >= indices.size()) return -1;
	else indices = Util::VSort(indices, pts);

	return (Int_t) indices[which]; 	
}

Double_t MT2tree::GetGenLeptEta(int which, int pid, int mother, double pt, double eta){
	Int_t index = GetGenLeptIndex(which, pid, mother, pt, eta);	
	if(index ==-1) return -999.99;
	else           return genlept[index].lv.Eta();	
}

Double_t MT2tree::GetGenLeptPt(int which, int pid, int mother, double pt, double eta){
	Int_t index = GetGenLeptIndex(which, pid, mother, pt, eta);	
	if(index ==-1) return -999.99;
	else           return genlept[index].lv.Pt();	
}

Bool_t MT2tree::GenLeptFromW(int pid, double pt, double eta, bool includeTau){
	bool good(false);
	for(int i=0; i<NGenLepts; ++i){
		if(abs(genlept[i].ID) !=pid                         ) continue;
		if( (!includeTau) && abs(genlept[i].MID) !=24       ) continue;
		if(   includeTau  && !((abs(genlept[i].MID)==15 && abs(genlept[i].GMID)==24 ) || abs(genlept[i].MID)==24)) continue;
		if(genlept[i].lv.Pt()       < pt                    ) continue;
		if(fabs(genlept[i].lv.Eta())>eta                    ) continue;
		good=true;
	}
	return good;	
}

Double_t MT2tree::GetLeptPt(int index){
	if(NEles + NMuons ==0) return -999;
	double pt_0 =0;
	double pt_1 =0;
	for(int i=0; i<NEles; ++i){
		if(ele[i].lv.Pt() > pt_0)                          {pt_1 = pt_0; pt_0 = ele[i].lv.Pt();}
		if(ele[i].lv.Pt() < pt_0 && ele[i].lv.Pt() > pt_1) {pt_1 = ele[i].lv.Pt();}
	}
	for(int i=0; i<NMuons; ++i){
		if(muo[i].lv.Pt() > pt_0)                          {pt_1 = pt_0; pt_0 = muo[i].lv.Pt();}
		if(muo[i].lv.Pt() < pt_0 && muo[i].lv.Pt() > pt_1) {pt_1 = muo[i].lv.Pt();}
	}

	if     (index==0) return pt_0;
	else if(index==1) return pt_1;
	else return -1;
}

Double_t MT2tree::ElClosestJet(){
	double dR=1000;
	for(int i=0; i<NEles; ++i){
		for(int j=0; j<NJets; ++j){
			if(ele[i].lv.DeltaR(jet[j].lv) < dR) {dR=ele[i].lv.DeltaR(jet[j].lv);}
		}	
	}
	return dR;
}

Int_t MT2tree::TopDecayMode(){
	// bit map: 
	// 1 = electron 1
	// 2 = electron 2
	// 4 = muon 1
	// 8 = muon 2
	// 16 = tau 1
	// 32 = tau 2
	// 64 = leptonic tau1
	// 128= leptonic tau2
	Int_t bit=0;
	Bool_t acceptance(true);
	for(int i=0; i<NGenLepts; ++i){
		if( abs(genlept[i].ID)==11 && abs(genlept[i].MID)==24 && abs(genlept[i].GMID)==6 ) {
			if( (bit & 1)==0) bit = bit | 1;
			else              bit = bit | 2;
		} 
		if( abs(genlept[i].ID)==13 && abs(genlept[i].MID)==24 && abs(genlept[i].GMID)==6 ) {
			if( (bit & 4)==0) bit = bit | 4;
			else              bit = bit | 8;
		} 
		if( abs(genlept[i].ID)==11 && abs(genlept[i].MID)==15 && abs(genlept[i].GMID)==24 ) {
			if     ( (bit & 16)==0) bit = bit | 16;
			else                    bit = bit | 32;
			if     ( (bit & 1 )==0) bit = bit |  1;
			else                    bit = bit |  2;
			if     ( (bit & 64)==0) bit = bit | 64;
			else                    bit = bit |128;
		} 
		if( abs(genlept[i].ID)==13 && abs(genlept[i].MID)==15 && abs(genlept[i].GMID)==24 ) {
			if     ( (bit & 16)==0) bit = bit | 16;
			else                    bit = bit | 32;
			if     ( (bit & 4 )==0) bit = bit |  4;
			else                    bit = bit |  8;
			if     ( (bit & 64)==0) bit = bit | 64;
			else                    bit = bit |128;
		}
		if( abs(genlept[i].ID)==16 && abs(genlept[i].MID)==24 && abs(genlept[i].GMID)==6 ){
			if     ( (bit & 16)==0) bit = bit | 16;
			else                    bit = bit | 32;
		}
	}
	return bit;
}


Bool_t MT2tree::TopDecayModeResult(Int_t nlepts){
	Int_t bit =TopDecayMode();
	if(nlepts == 1){ // semileptonic without leptonic tau
		if     ( (bit & 2 )==2 || (bit & 8)==8) return false; // more than one e/mu
		if     ( (bit & 64)==64)                return false; // at least one leptonic tau
		if     ( (bit & 1 )==1 && (bit & 4)==0) return true;
		else if( (bit & 1 )==0 && (bit & 4)==4) return true;
		else                                    return false;
	}else if(nlepts == 115){ // semileptonic with leptonic tau
		if     ( (bit & 2 )==2 || (bit & 8)==8) return false; // more than one e/mu
		if     ( (bit & 1 )==1 && (bit & 4)==0) return true;
		else if( (bit & 1 )==0 && (bit & 4)==4) return true;
		else                                    return false;
	}else if(nlepts == 215){ // fully leptonic with leptonic tau
		if     ( (bit & 2 )==2 || (bit & 8)==8) return true; // two eles or two muons
		if     ( (bit & 1 )==1 && (bit & 4)==4) return true; // one ele and one muo
		else                                    return false;
	}else if(nlepts == 2){ // fully leptonic without leptonic tau
		if     ( (bit & 64)==64               ) return false; // leptonic tau
		if     ( (bit & 2 )==2 || (bit & 8)==8) return true; // two eles or two muons
		if     ( (bit & 1 )==1 && (bit & 4)==4) return true; // one ele and one muo
		else                                    return false;
	}else if(nlepts == 0){ // fully hadronic without hadronic tau
		if     ( (bit & 1 )==1 || (bit & 4)==4) return false; // ele or muo
		if     ( (bit & 16)==16               ) return false; // tau
		else                                    return true;
	}else if(nlepts == 15){ // fully hadronic with hadronic tau
		if     ( (bit & 1 )==1 || (bit & 4)==4) return false; // ele or muo
		else                                    return true;
	}else if(nlepts ==11){
		if     ( (bit & 4 )==4 ) return false; // muon
		if     ( (bit & 2 )==2 ) return false; // two electron
		if     ( (bit & 1 )==1 ) return true;  // electron
		else                     return false;
	}else if(nlepts ==13){
		if     ( (bit & 1 )==1 ) return false; // ele
		if     ( (bit & 8 )==8 ) return false; // two muons
		if     ( (bit & 4 )==4 ) return true;  // muon
		else                     return false;
	}
	else                                           return false;
}

Bool_t MT2tree::SLTopAccept(double pt, double eta){
	for(int i=0; i<NGenLepts; ++i){
		if(    abs(genlept[i].ID)  !=11 && abs(genlept[i].ID) !=13                               ) continue;
		if( ! (abs(genlept[i].MID) ==15 && abs(genlept[i].GMID)==24 || abs(genlept[i].MID) ==24 )) continue;
		if(genlept[i].lv.Pt()>pt && fabs(genlept[i].lv.Eta()) < eta)                              return true;	
	}
	return false;
}

Double_t MT2tree::SLTopEta(double pt){
	for(int i=0; i<NGenLepts; ++i){
		if(    abs(genlept[i].ID)  !=11 && abs(genlept[i].ID) !=13                               ) continue;
		if( ! (abs(genlept[i].MID) ==15 && abs(genlept[i].GMID)==24 || abs(genlept[i].MID) ==24 )) continue;
		if(genlept[i].lv.Pt()>pt)      return genlept[i].lv.Eta();	
	}
	return -999.99;
}


Int_t MT2tree::WDecayMode(){
	// bit map:
	// 0 not recognized
	// 1= ele
	// 2= muo
	// 4= tau
	// 8= tau stable (i.e. problem in MC sample)
	Int_t result =0;
	for(int i=0; i<NGenLepts; ++i){
		if( abs(genlept[i].ID)==11 && abs(genlept[i].MID)==24 )                               {result = result | 1;} 
		if( abs(genlept[i].ID)==13 && abs(genlept[i].MID)==24 )                               {result = result | 2;} 
		if( abs(genlept[i].ID)==16 && abs(genlept[i].MID)==24 )                               {result = result | 4;}  // tau neutrino
		if( abs(genlept[i].ID)==11 && abs(genlept[i].MID)==15 && abs(genlept[i].GMID)==24 )   {result = result | 5;} 
		if( abs(genlept[i].ID)==13 && abs(genlept[i].MID)==15 && abs(genlept[i].GMID)==24 )   {result = result | 6;} 
		if( abs(genlept[i].ID)==15 && abs(genlept[i].MID)==24 )                               {result = result | 12;} // stable tau
	}
	return result;
}

ClassImp(MT2Misc)
ClassImp(MT2Znunu)
ClassImp(MT2PileUp)
ClassImp(MT2Trigger)
ClassImp(MT2Jet)
ClassImp(MT2Elec)
ClassImp(MT2Muon)
ClassImp(MT2Hemi)
ClassImp(MT2GenLept)
ClassImp(MT2tree)
