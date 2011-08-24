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
const bool RunEfficiency::ElAcceptance(int index) {
  
  if(!(fTR->ElPt[index]>5) )return false;
  if(!(fabs(fTR->ElEta[index]) < 2.4) ) return false;
  return true;

}

//_____________________________________________________________
const bool RunEfficiency::ElPassingTag(int n, int index) {

  //The tag has to have a more constrained acceptance
  if(n == 0) {
    if(!(fTR->ElPt[index]>20) )return false;
    if(!(fabs(fTR->ElEta[index]) < 2.4) ) return false;
    int elIDWP95 = fTR->ElIDsimpleWP95relIso[index];
    if (elIDWP95!=7) return false;
  } else if(n == 1) {
    if(!(fTR->ElPt[index]>20) )return false;
    if(!(fabs(fTR->ElEta[index]) < 2.4) ) return false;
    int elIDWP95 = fTR->ElIDsimpleWP95relIso[index];
    if (elIDWP95!=7) return false;
  } else if(n == 2) {  
    if(!(fTR->ElPt[index]>20) )return false;
    if(!(fabs(fTR->ElEta[index]) < 2.4) ) return false;
    int elIDWP95 = fTR->ElIDsimpleWP95relIso[index];
    if (elIDWP95!=7) return false;
  }

  return true;

}

//_____________________________________________________________
const bool RunEfficiency::ElPassingRecoProbe(int n, int index) {

  if(!ElAcceptance(index)) return false;
  if(n == 0) {
    if(!fTR->ElEcalDriven[index]) return false;
    if(fTR->ElCaloEnergy[index] < 10.) return false;
  } else if(n == 1) {
    if(!fTR->ElEcalDriven[index]) return false;
    if(fTR->ElCaloEnergy[index] < 10.) return false;
  } else if(n == 2) {
    if(!fTR->ElEcalDriven[index]) return false;
    if(fTR->ElCaloEnergy[index] < 10.) return false;
  }

  return true;

}


//_____________________________________________________________
const bool RunEfficiency::ElPassingIDProbe(int n, int index) {
  //The probe for the ID  has to be the passing probe of the reco
  return ElPassingRecoPProbe(n, index);   
}


//_____________________________________________________________
const bool RunEfficiency::ElPassingIsoProbe(int n, int index) {
  //The probe for the Iso  has to be the passing probe of the ID
  return ElPassingIDPProbe(n, index);   
}


//_____________________________________________________________
const bool RunEfficiency::ElPassingRecoPProbe(int n, int index) {
  if(!ElAcceptance(index)) return false;

  if(n == 0) {
    if ( !(fTR->ElNumberOfMissingInnerHits[index] <= 1 ) ) return false;
    if ( !(fabs(fTR->ElD0PV[index]) < 0.04) ) return false;
    if ( !(fabs(fTR->ElDzPV[index]) < 1.0 ) ) return false;
  } else if(n == 1) {
    if ( !(fTR->ElNumberOfMissingInnerHits[index] <= 1 ) ) return false;
    if ( !(fabs(fTR->ElD0PV[index]) < 0.04) ) return false;
    if ( !(fabs(fTR->ElDzPV[index]) < 1.0 ) ) return false;
  } else if(n == 2) {
    //Attention: MT2 is using pf leptons Reco, ID and Iso does not make sense here.
    return true;
  }
  return true;
}
 

//_____________________________________________________________
const bool RunEfficiency::ElPassingIDPProbe(int n, int index) {
  if(!ElAcceptance(index)) return false;

  if(n == 0) {
    int elIDWP95 = fTR->ElIDsimpleWP95relIso[index];
    if(elIDWP95!=7) return false;
  } else if(n == 1) {
    int elIDWP95 = fTR->ElIDsimpleWP95relIso[index];
    if(elIDWP95!=7) return false;
  } else if(n == 2) {
    return true;
  }
  return true;
}


//_____________________________________________________________
const bool RunEfficiency::ElPassingIsoPProbe(int n, int index,int pfindex) {
  if(!ElAcceptance(index)) return false;

  if(n == 0) {
    double hybridIso = fTR->ElRelIso03[index]*fTR->ElPt[index]/std::max((float)20.,fTR->ElPt[index]);
    if ( !(hybridIso < 0.15) ) return false;
  } else if(n == 1) {
    double hybridIso = fTR->ElRelIso03[index]*fTR->ElPt[index]/std::max((float)20.,fTR->ElPt[index]);
    if ( !(hybridIso < 0.15) ) return false;
  } else if(n == 2) {
    if(pfindex == -1) return false;
    if( ! IsCustomPfEl(pfindex) ) return false;
  }
  return true;
}


//_______________________________________________________________________
const int RunEfficiency::getPFElIndex(const int recoIndex) {

  int elIndex = 0;
  float px= fTR->ElPx[recoIndex];
  float py= fTR->ElPy[recoIndex];
  float pz= fTR->ElPz[recoIndex];
  float energy =  fTR->ElE[recoIndex];
  TLorentzVector tmpVector(px,py,pz,energy);

  for(elIndex=0;elIndex<fTR->PfEl3NObjs;elIndex++) {
    float pfpx= fTR->PfEl3Px[elIndex];
    float pfpy= fTR->PfEl3Py[elIndex];
    float pfpz= fTR->PfEl3Pz[elIndex];
    float pfenergy =  fTR->PfEl3E[elIndex];
    //if(pfpx == 0 && pfpy == 0 && pfpz == 0) break;
    TLorentzVector tmpPfVector(pfpx,pfpy,pfpz,pfenergy);
    if(tmpPfVector.DeltaR(tmpVector)<0.05) return elIndex; //maybe to tight
  }
  return -1;

}


//__________________________________________________________________________
const bool RunEfficiency::IsCustomPfEl(const int index){

  if(!(fTR->PfEl3Pt[index]>5) )return false;
  if(!(fabs(fTR->PfEl3Eta[index]) < 2.4) ) return false;
  if(!(fTR->PfElID95[index])) return false;
  
  return true;

}

