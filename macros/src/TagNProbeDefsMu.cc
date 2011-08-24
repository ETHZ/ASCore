#include "helper/Utilities.hh"
#include "RunEfficiency.hh"
#include "TF1.h"
#include <time.h>
#include <TRandom.h>
#include "TF1.h"


using namespace std;


#define jMax 30  // do not touch this
#define metMax 30
#define rMax 30
#define nAnalysis 3


//____________________________________________________________
const bool RunEfficiency::MuAcceptance(int index) {

  if(!(fTR->MuPt[index]>5) )return false;
  if(!(fabs(fTR->MuEta[index]) < 2.4) ) return false;

  return true;

}


//____________________________________________________________
const bool RunEfficiency::MuPassingTag(int n, int index){

  if(!(fTR->MuPt[index]>20) )return false;
  if(!(fabs(fTR->MuEta[index]) < 2.4) ) return false;

  if(n==0) {
    if ( !fTR->MuIsGMPT[index]) return false;
    if ( !fTR->MuIsGlobalMuon[index]) return false;
    if ( !fTR->MuIsTrackerMuon[index]) return false;
  } else if(n==1) {
    if ( !fTR->MuIsGMPT[index]) return false;
    if ( !fTR->MuIsGlobalMuon[index]) return false;
    if ( !fTR->MuIsTrackerMuon[index]) return false;
  } else if(n==1) {
    if ( !fTR->MuIsGMPT[index]) return false;
    if ( !fTR->MuIsGlobalMuon[index]) return false;
    if ( !fTR->MuIsTrackerMuon[index]) return false;
  }
  return true;
}


//__________________________________________________________________________
const bool RunEfficiency::MuPassingRecoProbe(const int n, const int index){

  // Acceptance cuts
  if(!MuAcceptance(index)) return false;

  if(n == 0) {
    if ( !fTR->MuIsTrackerMuon[index] ) return false;
  } else if(n == 1) {
    if ( !fTR->MuIsTrackerMuon[index] ) return false;
  } else if(n == 2) {
    if ( !fTR->MuIsTrackerMuon[index] ) return false;
  }
  return true;
}


//__________________________________________________________________________
const bool RunEfficiency::MuPassingIDProbe(const int n, const int index){
  return MuPassingRecoPProbe(n, index); 
}


//__________________________________________________________________________
const bool RunEfficiency::MuPassingIsoProbe(const int n, const int index){
  return MuPassingIDPProbe(n, index); 
}


//_________________________________________________________________________________
const bool RunEfficiency::MuPassingRecoPProbe(const int n, const int index){

  if(!MuAcceptance(index)) return false;

  if(n == 0) {
    if ( !fTR->MuIsTrackerMuon[index] ) return false;
    if ( !(fTR->MuNTkHits[index] >= 11) )     return false;
    if ( !(fTR->MuNPxHits[index] > 0) )       return false;
    if ( !(fTR->MuNMatches[index] > 1) )      return false;
    if ( !(fabs(fTR->MuD0PV[index]) < 0.02) ) return false;
    if ( !(fabs(fTR->MuDzPV[index]) < 1.0 ) ) return false;
    if ( !(fTR->MuPtE[index]/fTR->MuPt[index] < 0.1) ) return false;
  } else if(n == 1) {
    if ( !fTR->MuIsTrackerMuon[index] ) return false;
    if ( !(fTR->MuNTkHits[index] >= 11) )     return false;
    if ( !(fTR->MuNPxHits[index] > 0) )       return false;
    if ( !(fTR->MuNMatches[index] > 1) )      return false;
    if ( !(fabs(fTR->MuD0PV[index]) < 0.02) ) return false;
    if ( !(fabs(fTR->MuDzPV[index]) < 1.0 ) ) return false;
    if ( !(fTR->MuPtE[index]/fTR->MuPt[index] < 0.1) ) return false;
  } else if(n == 2) {
    if ( !fTR->MuIsTrackerMuon[index] ) return false;
    if ( !(fTR->MuNTkHits[index] >= 11) )     return false;
    if ( !(fTR->MuNPxHits[index] > 0) )       return false;
    if ( !(fTR->MuNMatches[index] > 1) )      return false;
    if ( !(fabs(fTR->MuD0PV[index]) < 0.02) ) return false;
    if ( !(fabs(fTR->MuDzPV[index]) < 1.0 ) ) return false;
    if ( !(fTR->MuPtE[index]/fTR->MuPt[index] < 0.1) ) return false;
  }
  return true;
  
}


//___________________________________________________________________________
const bool RunEfficiency::MuPassingIDPProbe(int n, int index){
  
  if(!MuAcceptance(index)) return false;
  
  if(n==0) {
    if ( !fTR->MuIsGMPT[index] )        return false;
    if ( !fTR->MuIsGlobalMuon[index] )  return false;
    if ( !fTR->MuIsTrackerMuon[index] ) return false;
  } else if(n==1) {
    if(fTR->MuIsGlobalMuon[index] == 0)  return false;
    if(fTR->MuIsTrackerMuon[index] == 0) return false;
  } else if(n==2) {
    return true; 
  }
  return true;
}


//______________________________________________________________________________
const bool RunEfficiency::MuPassingIsoPProbe(int n, int index, int pfindex){
  
  if(!MuAcceptance(index)) return false;
  
  if(n==0) {
    double hybridIso = fTR->MuRelIso03[index]*fTR->MuPt[index]/std::max((float)20.,fTR->MuPt[index]);
    if ( !(hybridIso < 0.15) ) return false;
    if ( !(fTR->MuPtE[index]/fTR->MuPt[index] < 0.1) ) return false;
  } else if(n==1) {
    if(fTR->MuRelIso03[index] > 0.15) return false;
  } else if(n==2) {
    if(pfindex == -1) return false;
  }
  return true;
} 


//_______________________________________________________________________
const int RunEfficiency::getPFMuIndex(const int recoIndex) {

  int muIndex = 0;
  float px= fTR->MuPx[recoIndex];
  float py= fTR->MuPy[recoIndex];
  float pz= fTR->MuPz[recoIndex];
  float energy =  fTR->MuE[recoIndex];
  TLorentzVector tmpVector(px,py,pz,energy);

  for(muIndex=0;muIndex<fTR->PfMu3NObjs;muIndex++) {
    float pfpx= fTR->PfMu3Px[muIndex];
    float pfpy= fTR->PfMu3Py[muIndex];
    float pfpz= fTR->PfMu3Pz[muIndex];
    float pfenergy =  fTR->PfMu3E[muIndex];
    //if(pfpx == 0 && pfpy == 0 && pfpz == 0) break;
    TLorentzVector tmpPfVector(pfpx,pfpy,pfpz,pfenergy);
    if(tmpPfVector.DeltaR(tmpVector)<0.05) return muIndex; //maybe to tight
  }
  return -1;

}


//__________________________________________________________________________
const bool RunEfficiency::IsCustomPfMu(const int index){
  
  if ( !(fTR->PfMu3Pt[index] > 5) )       return false;
  if ( !(fabs(fTR->PfMu3Eta[index])<2.4) ) return false;

  return true;
}



