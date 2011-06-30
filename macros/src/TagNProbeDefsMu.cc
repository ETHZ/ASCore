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

//************************************************************************************

/*

Steps: 		Tag	Probe	Passing Probe
(basic)		1	
Reco		-	2	3
ID		-	3	4
Isolation	-	4	5


1) Tag
2) Reco Probe
3) Passing Reco Probe
- 3) ID Probe
4) Passing ID Probe
- 4) Iso Probe
5) Iso Passing Probe

*/
 
///
/// 1
///

const bool RunEfficiency::MuPassingTag(int n, int index){


  if(n==0) // JZB
  {
	if ( !fTR->MuIsGMPT[index] )        return false;
	if ( !fTR->MuIsGlobalMuon[index] )  return false;
	if ( !fTR->MuIsTrackerMuon[index] ) return false;
  } else if(n==1) // same sign
  {
	if(fTR->MuIsGlobalMuon[index] == 0)  return false;
	if(fTR->MuIsTrackerMuon[index] == 0) return false;
  } else if(n==1) // mt2
  {
	if(getPFMuIndex(index)==-1) return false;
  }

  /*
  if(n == 0) {
    // Quality cuts
    if ( !fTR->MuIsGMPT[index] )        return false;
    if ( !fTR->MuIsGlobalMuon[index] )  return false;
    if ( !fTR->MuIsTrackerMuon[index] ) return false;
  } else if(n == 1) {
    if ( !fTR->MuIsGMPT[index] )        return false;
    if ( !fTR->MuIsGlobalMuon[index] )  return false;
    if ( !fTR->MuIsTrackerMuon[index] ) return false;
  } else if(n == 2) {
    if ( !fTR->MuIsGMPT[index] )        return false;
    if ( !fTR->MuIsGlobalMuon[index] )  return false;
    if ( !fTR->MuIsTrackerMuon[index] ) return false;
  }
*/
  return true;
}

const bool RunEfficiency::MuPassingRecoTag(int n, int index){
  return RunEfficiency::MuPassingTag(n,index);
}

const bool RunEfficiency::MuPassingIDTag(int n, int index) {
  return RunEfficiency::MuPassingTag(n,index);
}

const bool RunEfficiency::MuPassingIsoTag(int n, int index) {
  return RunEfficiency::MuPassingTag(n,index);
}

///--------------------------------------------------------------------------------
///
/// 2
///

const bool RunEfficiency::MuPassingRecoProbe(int n, int index){
  if(n!=2) {
	if ( !(fTR->MuPt[index] > 5) )       return false;
	if ( !(fabs(fTR->MuEta[index])<2.4) ) return false;
  } else {
	if(getPFMuIndex(index)==-1) return false;
  } 
  

  return true;
}

///--------------------------------------------------------------------------------
///
/// 3
///

const bool RunEfficiency::MuPassingRecoPProbe(const int n, const int index, const int pfindex){
  
  if(n!=2) {//jzb&same sign
	if ( !(fTR->MuPt[index] > 5) )       return false;
	if ( !(fabs(fTR->MuEta[index])<2.4) ) return false;
  } 

  if(n == 0) {//JZB
	if ( !(fTR->MuNTkHits[index] >= 11) )     return false;
	if ( !(fTR->MuNPxHits[index] > 0) )       return false;
	if ( !(fTR->MuNMatches[index] > 1) )      return false;
	if ( !(fabs(fTR->MuD0PV[index]) < 0.02) ) return false;
	if ( !(fabs(fTR->MuDzPV[index]) < 1.0 ) ) return false;
  } else if(n == 1) {//SameSign
	if(fTR->MuNChi2[index] > 10)   return false;
	if(fTR->MuPtE[index]/fTR->MuPt[index] > 0.1) return false;
	if(fTR->MuNTkHits[index] < 11) return false;
	if(fTR->MuNMuHits[index] < 1)  return false;
	if(fabs(fTR->MuD0PV[index]) > 0.02)    return false;
	if(fabs(fTR->MuDzPV[index]) > 1.00)    return false;
	if(fTR->MuIso03EMVetoEt[index] > 4.0)  return false;
	if(fTR->MuIso03HadVetoEt[index] > 6.0) return false;
  } else if(n == 2) {//MT2
//	if ( !fTR->MuIsTrackerMuon[index] ) return false;
	if(pfindex == -1) return false;
	cout << "don't know exactly what to do for mt2!" << endl;
  }
  return true;
}

/// 
/// 3a
///

const bool RunEfficiency::MuPassingIDProbe(int n, int index){
  return RunEfficiency::MuPassingRecoPProbe(n,index,getPFMuIndex(index));
} 

///--------------------------------------------------------------------------------
///
/// 4
///


const bool RunEfficiency::MuPassingIDPProbe(int n, int index, int pfindex){
  if(!RunEfficiency::MuPassingIDProbe(n,index)) return false;
  if(n==0) // JZB
  {
	if ( !fTR->MuIsGMPT[index] )        return false;
	if ( !fTR->MuIsGlobalMuon[index] )  return false;
	if ( !fTR->MuIsTrackerMuon[index] ) return false;
  } else if(n==1) // same sign
  {
	if(fTR->MuIsGlobalMuon[index] == 0)  return false;
	if(fTR->MuIsTrackerMuon[index] == 0) return false;
  } else if(n==2) // mt2
  {
	if(pfindex == -1) return false;
	cout << "don't know exactly what to do for mt2!" << endl;
  }
  return true;
}

///
/// 4a
///

const bool RunEfficiency::MuPassingIsoProbe(int n, int index){
  return RunEfficiency::MuPassingIDPProbe(n, index, getPFMuIndex(index));
} 

///--------------------------------------------------------------------------------
///
/// 5
///

const bool RunEfficiency::MuPassingIsoPProbe(int n, int index, int pfindex){
  RunEfficiency::MuPassingIDPProbe(n,index,pfindex);
  if(n==0) // JZB
  {
	double hybridIso = fTR->MuRelIso03[index]*fTR->MuPt[index]/std::max((float)20.,fTR->MuPt[index]);
	if ( !(hybridIso < 0.15) ) return false;
	if ( !(fTR->MuPtE[index]/fTR->MuPt[index] < 0.1) ) return false;
  } else if(n==1) // same sign
  {
	if(fTR->MuRelIso03[index] > 0.15) return false;
  } else if(n==2) // mt2
  {
	if(pfindex == -1) return false;
	cout << "don't know exactly what to do for mt2!" << endl;
  }
  return true;
} 

/*

Steps: 		Tag	Probe	Passing Probe
(basic)		1	
Reco		-	2	3
ID		-	3	4
Isolation	-	4	5


1) Tag
2) Reco Probe
3) Passing Reco Probe
- 3) ID Probe
4) Passing ID Probe
- 4) Iso Probe
5) Iso Passing Probe

*/