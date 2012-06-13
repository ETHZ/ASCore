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
  if(!(fabs(fTR->ElEta[index]) < 2.5) ) return false;
  return true;

}


//_____________________________________________________________
const bool RunEfficiency::ElPassingTag(int n, int index) {

  //The tag has to have a more constrained acceptance
  if(n == 0) {
    if(!(fTR->ElPt[index]>20) )return false;
    if(!(fabs(fTR->ElEta[index]) < 2.5) ) return false;
    int elIDWP90 = fTR->ElIDsimpleWP90relIso[index];
    if(!(elIDWP90&1)) return false;
    if ( fabs(fTR->ElEta[index]) < 1.479 ) { // Barrel
      if ( !(fTR->ElHcalOverEcal[index]<0.1) ) return false;    // 0.12
      if ( !(fTR->ElSigmaIetaIeta[index]<0.011) ) return false; // 0.01
      if ( !(fabs(fTR->ElDeltaPhiSuperClusterAtVtx[index])<0.15) ) return false; // 0.8
      if ( !(fabs(fTR->ElDeltaEtaSuperClusterAtVtx[index])<0.01) ) return false; // 0.007
    } else { // Endcap
      if ( !(fTR->ElHcalOverEcal[index]<0.075) ) return false;  // 0.05
      if ( !(fTR->ElSigmaIetaIeta[index]<0.031) ) return false; // 0.03
      if ( !(fabs(fTR->ElDeltaPhiSuperClusterAtVtx[index])<0.1) ) return false;  // 0.7
      if ( !(fabs(fTR->ElDeltaEtaSuperClusterAtVtx[index])<0.01) ) return false; // 0.009
    }
  } else if(n == 1) {
    if(!(fTR->ElPt[index]>20) )return false;
    if(!(fabs(fTR->ElEta[index]) < 2.5) ) return false;
    int elIDWP90 = fTR->ElIDsimpleWP90relIso[index];
    if(!(elIDWP90&1)) return false;
    if ( fabs(fTR->ElEta[index]) < 1.479 ) { // Barrel
      if ( !(fTR->ElHcalOverEcal[index]<0.1) ) return false;    // 0.12
      if ( !(fTR->ElSigmaIetaIeta[index]<0.011) ) return false; // 0.01
      if ( !(fabs(fTR->ElDeltaPhiSuperClusterAtVtx[index])<0.15) ) return false; // 0.8
      if ( !(fabs(fTR->ElDeltaEtaSuperClusterAtVtx[index])<0.01) ) return false; // 0.007
    } else { // Endcap
      if ( !(fTR->ElHcalOverEcal[index]<0.075) ) return false;  // 0.05
      if ( !(fTR->ElSigmaIetaIeta[index]<0.031) ) return false; // 0.03
      if ( !(fabs(fTR->ElDeltaPhiSuperClusterAtVtx[index])<0.1) ) return false;  // 0.7
      if ( !(fabs(fTR->ElDeltaEtaSuperClusterAtVtx[index])<0.01) ) return false; // 0.009
    }
  } else if(n == 2) {  
    if(!(fTR->ElPt[index]>20) )return false;
    if(!(fabs(fTR->ElEta[index]) < 2.5) ) return false;
    int elIDWP90 = fTR->ElIDsimpleWP90relIso[index];
    if(!(elIDWP90&1)) return false;
    if ( fabs(fTR->ElEta[index]) < 1.479 ) { // Barrel
      if ( !(fTR->ElHcalOverEcal[index]<0.1) ) return false;    // 0.12
      if ( !(fTR->ElSigmaIetaIeta[index]<0.011) ) return false; // 0.01
      if ( !(fabs(fTR->ElDeltaPhiSuperClusterAtVtx[index])<0.15) ) return false; // 0.8
      if ( !(fabs(fTR->ElDeltaEtaSuperClusterAtVtx[index])<0.01) ) return false; // 0.007
    } else { // Endcap
      if ( !(fTR->ElHcalOverEcal[index]<0.075) ) return false;  // 0.05
      if ( !(fTR->ElSigmaIetaIeta[index]<0.031) ) return false; // 0.03
      if ( !(fabs(fTR->ElDeltaPhiSuperClusterAtVtx[index])<0.1) ) return false;  // 0.7
      if ( !(fabs(fTR->ElDeltaEtaSuperClusterAtVtx[index])<0.01) ) return false; // 0.009
    }
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
    if ( !(fTR->ElESuperClusterOverP[index]*fTR->ElTrkMomAtVtx[index]>10) ) return false;
    if ( !(fabs(fTR->ElD0PV[index]) < 0.04) ) return false;
    if ( !(fabs(fTR->ElDzPV[index]) < 1.0 ) ) return false;
  } else if(n == 1) {
    if ( !(fTR->ElNumberOfMissingInnerHits[index] <= 1 ) ) return false;
    if ( !(fTR->ElESuperClusterOverP[index]*fTR->ElTrkMomAtVtx[index]>10) ) return false;
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
    int elIDWP90 = fTR->ElIDsimpleWP90relIso[index];
    if(!(elIDWP90&1)) return false;
    if ( fabs(fTR->ElEta[index]) < 1.479 ) { // Barrel
      if ( !(fTR->ElHcalOverEcal[index]<0.1) ) return false;    // 0.12
      if ( !(fTR->ElSigmaIetaIeta[index]<0.011) ) return false; // 0.01
      if ( !(fabs(fTR->ElDeltaPhiSuperClusterAtVtx[index])<0.15) ) return false; // 0.8
      if ( !(fabs(fTR->ElDeltaEtaSuperClusterAtVtx[index])<0.01) ) return false; // 0.007
    } else { // Endcap
      if ( !(fTR->ElHcalOverEcal[index]<0.075) ) return false;  // 0.05
      if ( !(fTR->ElSigmaIetaIeta[index]<0.031) ) return false; // 0.03
      if ( !(fabs(fTR->ElDeltaPhiSuperClusterAtVtx[index])<0.1) ) return false;  // 0.7
      if ( !(fabs(fTR->ElDeltaEtaSuperClusterAtVtx[index])<0.01) ) return false; // 0.009
    }
  } else if(n == 1) {
    int elIDWP90 = fTR->ElIDsimpleWP90relIso[index];
    if(!(elIDWP90&1)) return false;
    if ( fabs(fTR->ElEta[index]) < 1.479 ) { // Barrel
      if ( !(fTR->ElHcalOverEcal[index]<0.1) ) return false;    // 0.12
      if ( !(fTR->ElSigmaIetaIeta[index]<0.011) ) return false; // 0.01
      if ( !(fabs(fTR->ElDeltaPhiSuperClusterAtVtx[index])<0.15) ) return false; // 0.8
      if ( !(fabs(fTR->ElDeltaEtaSuperClusterAtVtx[index])<0.01) ) return false; // 0.007
    } else { // Endcap
      if ( !(fTR->ElHcalOverEcal[index]<0.075) ) return false;  // 0.05
      if ( !(fTR->ElSigmaIetaIeta[index]<0.031) ) return false; // 0.03
      if ( !(fabs(fTR->ElDeltaPhiSuperClusterAtVtx[index])<0.1) ) return false;  // 0.7
      if ( !(fabs(fTR->ElDeltaEtaSuperClusterAtVtx[index])<0.01) ) return false; // 0.009
    }
  } else if(n == 2) {
    return true;
  }
  return true;
}


//_____________________________________________________________
const bool RunEfficiency::ElPassingIsoPProbe(int n, int index,int pfindex) {
  if(!ElAcceptance(index)) return false;

  if(n == 0) {
    double pedestal = 0.;
    if ( fabs(fTR->ElEta[index]) < 1.479 ) pedestal = 1.0;
    double iso = fTR->ElDR03TkSumPt[index]+std::max(fTR->ElDR03EcalRecHitSumEt[index]-pedestal,0.)+fTR->ElDR03HcalTowerSumEt[index];
    double hybridIso = iso/fTR->ElPt[index]; // Ditched the flat iso below 20GeV (irrelevant anyway)
    if ( !(hybridIso < 0.15) ) return false;
  } else if(n == 1) {
    double pedestal = 0.;
    if ( fabs(fTR->ElEta[index]) < 1.479 ) pedestal = 1.0;
    double iso = fTR->ElDR03TkSumPt[index]+std::max(fTR->ElDR03EcalRecHitSumEt[index]-pedestal,0.)+fTR->ElDR03HcalTowerSumEt[index];
    double hybridIso = iso/fTR->ElPt[index]; // Ditched the flat iso below 20GeV (irrelevant anyway)
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
  if(!(fabs(fTR->PfEl3Eta[index]) < 2.5) ) return false;
  if(!(fTR->PfElID95[index])) return false;
  
  return true;

}

