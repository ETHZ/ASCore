#include "base/TreeReader.hh"
#include <stdlib.h>

#include "helper/pdgparticle.hh"
#include "base/UserAnalysisBase.hh"
#include "TH1I.h"

using namespace std;

UserAnalysisBase::UserAnalysisBase(TreeReader *tr){
	fTR = tr;
	fTlat = new TLatex();
	fVerbose = false;
}

UserAnalysisBase::~UserAnalysisBase(){
}

void UserAnalysisBase::BeginRun(Int_t& run) {
  // Called when run changes
  // Need to re-create HLT map in case it changed
	
  GetHLTNames(run);
}

void UserAnalysisBase::ReadPDGTable(const char* filename){
// Fills the fPDGMap map from a textfile to associate pdgids with names
	int pdgid(0), type(0);
	string Name, Texname, Typename;
	ifstream IN(filename);
	char buff1[200], buff2[200], buff3[200];
	char readbuff[200];
	// Loop over lines of datafile
	while( IN.getline(readbuff, 200, '\n') ){
		if (readbuff[0] == '#') {continue;} // Skip lines commented with '#'
		sscanf(readbuff, "%d %s %d %s %s", &type, buff1, &pdgid, buff2, buff3);
		// Convert chararrays to strings
		Typename = string(buff1); Name = string(buff2); Texname = string(buff3);
		pdgparticle *p = new pdgparticle(pdgid, Name, Texname, type, Typename);
		// Fill map
		fPDGMap[pdgid] = *p;
	}
}

int UserAnalysisBase::GetPDGParticle(pdgparticle &part, int id){
	if( fPDGMap.empty() ){
		if(fVerbose > 0) cout << "UserAnalysisBase::GetPDGParticle ==> PDGMap not filled!" << endl;
		return -1;
	}
	else{
		map<int, pdgparticle>::iterator it = fPDGMap.find(id);
		if(it == fPDGMap.end()){
			if(fVerbose > 1) cout << "UserAnalysisBase::GetPDGParticle ==> PDGParticle with ID " << id << " not found!" << endl;
			return -1;
		}
		else{
			part = it->second;
			return 1;
		}
	}
	return 0;
}

void UserAnalysisBase::GetHLTNames(Int_t& run){
	
	cout << "GetHLTNames... ";
	TFile *f = fTR->fChain->GetCurrentFile();
	TTree* runTree = (TTree*)f->Get("analyze/RunInfo");
	std::vector<std::string>* HLTNames;
	if ( !runTree ) {
		std::cerr << "!!! UserAnalysisBase::GetHLTNames "
		          << "Couldn't get analyze/RunInfo tree" << std::endl;
		return;
	}
	//TH1I *hlt_stats = (TH1I*)f->Get("analyze/HLTTriggerStats");

	if ( fVerbose>0 ) std::cout << "Retrieving HLTNames for run " << run << std::endl;
	runTree->SetBranchAddress("HLTNames",&HLTNames);
	runTree->GetEntryWithIndex(run);

	for( int i=0; i < HLTNames->size(); i++ ) 
		fHLTLabelMap[(*HLTNames)[i]] = i; 
}

int UserAnalysisBase::GetHLTBit(string theHltName){
	if( fHLTLabelMap.empty() ) return -1;
	else{
		map<string,int>::iterator it = fHLTLabelMap.find(theHltName);
		if(it == fHLTLabelMap.end()){
			if(fVerbose > 0) cout << "UserAnalysisBase::GetHLTBit ==> Bit with name " << theHltName << " not found!" << endl;
			return -1;
		}
		else{
			return it->second;
		}
	}
}

bool UserAnalysisBase::GetHLTResult(string theHltName){
	if( fHLTLabelMap.empty() ) return false;
	else{
		int bit = GetHLTBit(theHltName);
		if(bit == -1){
			if(fVerbose > 0) cout << "UserAnalysisBase::GetHLTResult ==> Bit with name " << theHltName << " not found!" << endl;
			return false;
		}
		else return (bool)fTR->HLTResults[bit];
	}
}

// Object selections:
bool UserAnalysisBase::IsGoodJ_TDL(int index) {
	if( fTR->JPt[index] < 30. ) 			return false;
	if( fabs(fTR->JEta[index]) > 2.5 ) 		return false;
	if (fTR->JEMfrac[index] <= 0.01)  		return false;
	if (fTR->JID_n90Hits[index] <= 1) 		return false;
	if (fTR->JID_HPD[index] >= 0.98)  		return false;

	return true;
}

bool UserAnalysisBase::IsGoodbJ_TDL(int index) {
	if(! IsGoodJ_TDL(index)) return false;
	if(fTR->JbTagProbSimpSVHighPur[index] < 2.0) return false;
	if(fabs(fTR->JEta[index]) > 2.0) return false;
	return true;
}

bool UserAnalysisBase::IsGoodBasicPFJet(int index) {
	// Basic PF jet cleaning and ID cuts
	if(fTR->PFJPt[index] < 30) return false;
	if(fabs(fTR->PFJEta[index]) > 3.0) return false;
	// Loose PF jet ID (WARNING: HF not included in our ntuple)
	// See PhysicsTools/SelectorUtils/interface/PFJetIDSelectionFunctor.h
	if ( !(fTR->PFJNConstituents[index] > 1) )    return false;
	if ( !(fTR->PFJNeuEmfrac[index]     < 0.99) ) return false;
	if ( !(fTR->PFJNeuHadfrac[index]    < 0.99) ) return false;
	if (fabs(fTR->PFJEta[index]) < 2.4 ) { // Cuts for |eta|<2.4
		if ( !(fTR->PFJChEmfrac[index]  < 0.99) )  return false;
		if ( !(fTR->PFJChHadfrac[index] > 0.00) )  return false;
		if ( !(fTR->PFJChMult[index]    > 0   ) )  return false;
	}
	return true;
}

bool UserAnalysisBase::IsGoodBasicJet(int index){
	// Basic Jet cleaning and ID cuts
	if(fTR->JPt[index] < 30) return false;
	if(fabs(fTR->JEta[index]) > 3.0) return false;

	if(fTR->JEt[index] - fTR->JPt[index] < -0.0001 ) return false;
	if(fTR->JID_n90Hits[index] < 2 ) return false;
	if(fTR->JID_HPD[index] > 0.98 ) return false;
	if(fTR->JID_RBX[index] > 0.95 ) return false;
	if(fTR->JEMfrac[index] > 1. ) return false;
	if(fTR->JEMfrac[index] < 0.01 ) return false;
	// Have a linearly decreasing cut value in the transition region where
	// the jet cones intersects with the tracker acceptance, i.e. between
	// eta 1.9 and 2.9
	const double chmin = 0.05;
	double temp = chmin;
	if(fabs(fTR->JEta[index]) > 1.9) temp = chmin * (1. - fabs(fTR->JEta[index]) + 1.9);
	if(fabs(fTR->JEta[index]) > 2.9) temp = 0.;
	if( fTR->JChfrac[index] < temp && fabs(fTR->JEta[index]) < 2.9) return false;
	return true;
}

bool UserAnalysisBase::IsGoodBasicMu(int index){
	// Basic muon cleaning and ID
	if(fTR->MuIsGlobalMuon[index] == 0)  return false;
	if(fTR->MuIsTrackerMuon[index] == 0) return false;

	if(fTR->MuPt[index] < 5)          return false;
	if(fabs(fTR->MuEta[index]) > 2.4) return false;

	if(fTR->MuNChi2[index] > 10)   return false;
	if(fTR->MuNTkHits[index] < 11) return false;
	if(fTR->MuNMuHits[index] < 1)  return false;

	if(fabs(fTR->MuD0PV[index]) > 0.02)    return false;
	if(fTR->MuIso03EMVetoEt[index] > 4.0)  return false;
	if(fTR->MuIso03HadVetoEt[index] > 6.0) return false;

	double iso = fTR->MuRelIso03[index];
	double pt = fTR->MuPt[index];
	double hybiso = iso*pt / std::max(20.,pt);
	if(hybiso > 1.0) return false;

	return true;
}

bool UserAnalysisBase::IsGoodMu_TDL(int index){
	if(fabs(fTR->MuEta[index])>2.5) return false;
	if(fTR->MuIsGMPT[index] == 0) return false;
	if(fTR->MuIsGlobalMuon[index] == 0) return false;
	if(fTR->MuNChi2[index] > 10) return false;
	if(fabs(fTR->MuD0BS[index]) > 0.02) return false;
	if(fTR->MuNTkHits[index] < 11) return false;
	if(fTR->MuPt[index] < 10) return false;
	if(fTR->MuRelIso03[index] > 0.15) return false;
	return true;
}

bool UserAnalysisBase::IsTightMu(int index){
	if(!IsGoodBasicMu(index)) return false;
	if(fTR->MuRelIso03[index] > 0.15) return false;
	return true;
}

bool UserAnalysisBase::IsLooseMu(int index){
	if(!IsGoodBasicMu(index)) return false;
	if(fTR->MuRelIso03[index] > 1.0) return false;
	return true;
}

bool UserAnalysisBase::IsLooseNoTightMu(int index){
	if(IsLooseMu(index) && !IsTightMu(index)) return true;
	else return false;
}

bool UserAnalysisBase::IsGoodEl_TDL(int index){
	// ---- electron selection from Top-Dilepton group  ----
	// ---- El id WP90
	double etaEgapUpperEdge         = 1.5660;
	double etaEgapLowerEdge         = 1.4442;
	
	if( fTR->ElPt[index] < 10. ) 	return false;
	if( fabs(fTR->ElEta[index]) > 2.5 ) return false;
	if( fabs(fTR->ElEta[index]) > etaEgapLowerEdge && fabs(fTR->ElEta[index]) < etaEgapUpperEdge) return false;
	if( fTR->ElD0BS[index] >= 0.04 ) return false;
	
	// conversion rejection && ID WP90
	if(fTR->ElIDsimpleWPrelIso[index]!=5 && fTR->ElIDsimpleWPrelIso[index]!=7) return false;	
	
	// isolation
	double elIsoEcal;
	if( fabs(fTR->ElEta[index]) < etaEgapLowerEdge ){
		elIsoEcal = fTR->ElDR03EcalRecHitSumEt[index] - 1.;
		elIsoEcal = ((0. > elIsoEcal) ? 0. : elIsoEcal);
	} else if(fabs(fTR->ElEta[index]) > etaEgapUpperEdge ){
		elIsoEcal = fTR->ElDR03EcalRecHitSumEt[index];
	}
	double pt = ((20. > fTR->ElEt[index]) ? 20. : fTR->ElEt[index]);
	double elIso = (fTR->ElDR03TkSumPt[index] + elIsoEcal + fTR->ElDR03HcalTowerSumEt[index]) / pt;
	if( elIso > 0.15 ) return false;
	
	// rejection if matched to muon	
	for( int im = 0; im < fTR->NMus; ++im ){
		if( (fTR->MuIsGlobalMuon[im] == 0 || fTR->MuIsTrackerMuon[im] == 0) && ( fTR->MuNTkHits[im] > 10) ) {
			double deltaR = Util::GetDeltaR(fTR->ElEta[index], fTR->MuEta[im], fTR->ElPhi[index], fTR->MuPhi[im]);
			if(deltaR <= 0.1) return false;
		}
	}
	
	return true;
}

bool UserAnalysisBase::IsGoodBasicEl(int index){
	// Basic electron cleaning and ID cuts
	// (corresponding to the 95% eff. WP without isolation cuts
	// and with coversion rejection cuts)
	
	double pt  = fTR->ElPt[index];
	double eta = fTR->ElEta[index];

	// barrel cut values
	double ElecHoverEBarmax         = 0.15;
	double ElecSigmaEtaEtaBarmin    = 0.002;
	double ElecSigmaEtaEtaBarmax    = 0.01;
	double ElecEoverPInBarmin       = 0.;
	double ElecDeltaEtaInBarmax     = 0.007;
	double ElecDeltaPhiInBarmax     = 0.8;
	double ElecDeltaPhiOutBarmax    = 999.0;
	// endcaps cut values
	double ElecHoverEEndmax         = 0.07;
	double ElecSigmaEtaEtaEndmax    = 0.03;
	double ElecEoverPInEndmin       = 0.;
	double ElecDeltaEtaInEndmax     = 999.0; // not used due to misalignement in the EndCaps
	double ElecDeltaPhiInEndmax     = 0.7;
	double ElecDeltaPhiOutEndmax    = 999.0;
	// kinematics + ECAL gap cut values
	double etamax                   = 3.0;
	double etaEgapUpperEdge         = 1.560;
	double etaEgapLowerEdge         = 1.442;
	double pTmin                    = 10.;
	
	// electron ID
	if( fabs(eta) < 1.479 ){ // Barrel
		if( fTR->ElHcalOverEcal[index] > ElecHoverEBarmax ) return false;
		if( fTR->ElESuperClusterOverP[index] < ElecEoverPInBarmin ) return false;
		if( fTR->ElSigmaIetaIeta[index] > ElecSigmaEtaEtaBarmax ) return false;
		if( fTR->ElSigmaIetaIeta[index] < ElecSigmaEtaEtaBarmin ) return false;
		if( fabs(fTR->ElDeltaEtaSuperClusterAtVtx[index]) > ElecDeltaEtaInBarmax ) return false;
		if( fabs(fTR->ElDeltaPhiSuperClusterAtVtx[index]) > ElecDeltaPhiInBarmax ) return false;
		if( fabs(fTR->ElDeltaPhiSeedClusterAtCalo[index]) > ElecDeltaPhiOutBarmax ) return false;
	}
	else{ // EndCap
		if( fTR->ElHcalOverEcal[index] > ElecHoverEEndmax ) return false;
		if( fTR->ElESuperClusterOverP[index] < ElecEoverPInEndmin ) return false;
		if( fTR->ElSigmaIetaIeta[index] > ElecSigmaEtaEtaEndmax ) return false;
		if( fabs(fTR->ElDeltaEtaSuperClusterAtVtx[index]) > ElecDeltaEtaInEndmax ) return false;
		if( fabs(fTR->ElDeltaPhiSuperClusterAtVtx[index]) > ElecDeltaPhiInEndmax ) return false;
		if( fabs(fTR->ElDeltaPhiSeedClusterAtCalo[index]) > ElecDeltaPhiOutEndmax ) return false;
	}

	// electron kinematic cuts + the ECAL gap veto
	if ( (fabs(eta) > etamax) ||
	   ( (fabs(eta) > etaEgapLowerEdge) && (fabs(eta) < etaEgapUpperEdge) ) ) return false;	
	if( pt < pTmin) return false;

	return true;
}

bool UserAnalysisBase::IsGoodElId_WP90(int index){
	// Converted electrons with WP90 cleaning and ID cuts
	
	double pt	= fTR->ElPt[index];
	double eta	= fTR->ElEta[index];
	
	// barrel cut values
	double ElecHoverEBarmax			= 0.12;
	double ElecSigmaEtaEtaBarmin	= 0.002;
	double ElecSigmaEtaEtaBarmax	= 0.01;
	double ElecEoverPInBarmin		= 0.;
	double ElecDeltaEtaInBarmax		= 0.007;
	double ElecDeltaPhiInBarmax		= 0.8;
	double ElecDeltaPhiOutBarmax	= 999.0;
	// endcaps cut values
	double ElecHoverEEndmax			= 0.05;
	double ElecSigmaEtaEtaEndmax	= 0.03;
	double ElecEoverPInEndmin		= 0.;
	double ElecDeltaEtaInEndmax		= 0.009; // not used with 3_6_X data due to misalignement in the EndCaps
	double ElecDeltaPhiInEndmax		= 0.7;
	double ElecDeltaPhiOutEndmax	= 999.0;
	// kinematics + ECAL gap cut values
	double etamax					= 3.0;
	double etaEgapUpperEdge			= 1.560;
	double etaEgapLowerEdge			= 1.442;
	double pTmin					= 10.;
	
	// electron ID
	if( fabs(eta) < 1.479 ){ // Barrel
		if(		 fTR->ElHcalOverEcal[index]					> ElecHoverEBarmax      ) return false;
		if(		 fTR->ElESuperClusterOverP[index]			< ElecEoverPInBarmin    ) return false;
		if(		 fTR->ElSigmaIetaIeta[index]				> ElecSigmaEtaEtaBarmax ) return false;
		if(		 fTR->ElSigmaIetaIeta[index]				< ElecSigmaEtaEtaBarmin ) return false;
		if( fabs(fTR->ElDeltaEtaSuperClusterAtVtx[index])	> ElecDeltaEtaInBarmax  ) return false;
		if( fabs(fTR->ElDeltaPhiSuperClusterAtVtx[index])	> ElecDeltaPhiInBarmax  ) return false;
		if( fabs(fTR->ElDeltaPhiSeedClusterAtCalo[index])	> ElecDeltaPhiOutBarmax ) return false;
	}
	else{ // EndCap
		if( 	 fTR->ElHcalOverEcal[index]					> ElecHoverEEndmax      ) return false;
		if(		 fTR->ElESuperClusterOverP[index]			< ElecEoverPInEndmin    ) return false;
		if(		 fTR->ElSigmaIetaIeta[index]				> ElecSigmaEtaEtaEndmax ) return false;
		if( fabs(fTR->ElDeltaEtaSuperClusterAtVtx[index])	> ElecDeltaEtaInEndmax  ) return false;
		if( fabs(fTR->ElDeltaPhiSuperClusterAtVtx[index])	> ElecDeltaPhiInEndmax  ) return false;
		if( fabs(fTR->ElDeltaPhiSeedClusterAtCalo[index])	> ElecDeltaPhiOutEndmax ) return false;
	}
	
	// electron kinematic cuts + the ECAL gap veto
	if ( (fabs(eta) > etamax) ||
		( (fabs(eta) > etaEgapLowerEdge) && (fabs(eta) < etaEgapUpperEdge) ) )	return false;	
	if( pt < pTmin) return false;
	
	return true;
}

bool UserAnalysisBase::IsGoodElId_WP80(int index){
	// Converted electrons with WP80 cleaning and ID cuts
	
	double pt	= fTR->ElPt[index];
	double eta	= fTR->ElEta[index];
	
	// barrel cut values
	double ElecHoverEBarmax			= 0.04;
	double ElecSigmaEtaEtaBarmin	= 0.002;
	double ElecSigmaEtaEtaBarmax	= 0.01;
	double ElecEoverPInBarmin		= 0.;
	double ElecDeltaEtaInBarmax		= 0.004;
	double ElecDeltaPhiInBarmax		= 0.6;
	double ElecDeltaPhiOutBarmax	= 999.0;
	// endcaps cut values
	double ElecHoverEEndmax			= 0.025;
	double ElecSigmaEtaEtaEndmax	= 0.03;
	double ElecEoverPInEndmin		= 0.;
	double ElecDeltaEtaInEndmax		= 0.007; // not used with 3_6_X data due to misalignement in the EndCaps
	double ElecDeltaPhiInEndmax		= 0.03;
	double ElecDeltaPhiOutEndmax	= 999.0;
	// kinematics + ECAL gap cut values
	double etamax					= 3.0;
	double etaEgapUpperEdge			= 1.560;
	double etaEgapLowerEdge			= 1.442;
	double pTmin					= 5.;
	
	// electron ID
	if( fabs(eta) < 1.479 ){ // Barrel
		if(		 fTR->ElHcalOverEcal[index]					> ElecHoverEBarmax      ) return false;
		if(		 fTR->ElESuperClusterOverP[index]			< ElecEoverPInBarmin    ) return false;
		if(		 fTR->ElSigmaIetaIeta[index]				> ElecSigmaEtaEtaBarmax ) return false;
		if(		 fTR->ElSigmaIetaIeta[index]				< ElecSigmaEtaEtaBarmin ) return false;
		if( fabs(fTR->ElDeltaEtaSuperClusterAtVtx[index])	> ElecDeltaEtaInBarmax  ) return false;
		if( fabs(fTR->ElDeltaPhiSuperClusterAtVtx[index])	> ElecDeltaPhiInBarmax  ) return false;
		if( fabs(fTR->ElDeltaPhiSeedClusterAtCalo[index])	> ElecDeltaPhiOutBarmax ) return false;
	}
	else{ // EndCap
		if( 	 fTR->ElHcalOverEcal[index]					> ElecHoverEEndmax      ) return false;
		if(		 fTR->ElESuperClusterOverP[index]			< ElecEoverPInEndmin    ) return false;
		if(		 fTR->ElSigmaIetaIeta[index]				> ElecSigmaEtaEtaEndmax ) return false;
		if( fabs(fTR->ElDeltaEtaSuperClusterAtVtx[index])	> ElecDeltaEtaInEndmax  ) return false;
		if( fabs(fTR->ElDeltaPhiSuperClusterAtVtx[index])	> ElecDeltaPhiInEndmax  ) return false;
		if( fabs(fTR->ElDeltaPhiSeedClusterAtCalo[index])	> ElecDeltaPhiOutEndmax ) return false;
	}
	
	// electron kinematic cuts + the ECAL gap veto
	if ( (fabs(eta) > etamax) ||
		( (fabs(eta) > etaEgapLowerEdge) && (fabs(eta) < etaEgapUpperEdge) ) )	return false;	
	if( pt < pTmin) return false;
	
	return true;
}

bool UserAnalysisBase::IsConvertedEl_WP90(int index){
	// conversions cut values for WP90
	double ElecConvPartTrackDistmin	= 0.02;
	double ElecConvPartTrackDCotmin	= 0.02;
	double ElecNMissHitsmax			= 1.;
	
	// return true if electron is conversion-like
	if ( fTR->ElNumberOfMissingInnerHits[index]	> ElecNMissHitsmax )		 return true;
	if ( fabs(fTR->ElConvPartnerTrkDist[index])	< ElecConvPartTrackDistmin &&
		fabs(fTR->ElConvPartnerTrkDCot[index])	< ElecConvPartTrackDCotmin ) return true;
	
	return false;
}

bool UserAnalysisBase::IsConvertedEl_WP80(int index){
	// conversions cut values
	double ElecConvPartTrackDistmin	= 0.02;
	double ElecConvPartTrackDCotmin	= 0.02;
	double ElecNMissHitsmax			= 0.;
	
	// return true if electron is conversion-like
	if (  fTR->ElNumberOfMissingInnerHits[index]	> ElecNMissHitsmax )		 return true;
	if (  fabs(fTR->ElConvPartnerTrkDist[index])	< ElecConvPartTrackDistmin &&
		fabs(fTR->ElConvPartnerTrkDCot[index])	< ElecConvPartTrackDCotmin ) return true;
	
	return false;
}

bool UserAnalysisBase::IsIsolatedEl_WP90(int index){
	// combined relative isolation cut values for WP90
	double ElecCombRelIsoEBarmax = 0.10;
	double ElecCombRelIsoEEndmax = 0.07;	
	return IsIsolatedEl(index, ElecCombRelIsoEBarmax, ElecCombRelIsoEEndmax);
}

bool UserAnalysisBase::IsIsolatedEl_WP80(int index){
	// combined relative isolation cut values for WP80
	double ElecCombRelIsoEBarmax = 0.07;
	double ElecCombRelIsoEEndmax = 0.06;
	return IsIsolatedEl(index, ElecCombRelIsoEBarmax, ElecCombRelIsoEEndmax);
}

bool UserAnalysisBase::IsIsolatedEl_RA5(int index){
	// combined relative isolation cut values for WP90
	double ElecCombRelIsoEBarmax = 0.1;
	double ElecCombRelIsoEEndmax = 0.1;	
	return IsIsolatedEl(index, ElecCombRelIsoEBarmax, ElecCombRelIsoEEndmax);
}

bool UserAnalysisBase::IsIsolatedEl_LooseIso(int index){
	// combined relative isolation cut values for loose electrons
	double ElecCombRelIsoEBarmax = 1.;
	double ElecCombRelIsoEEndmax = 0.6;
	return IsIsolatedEl(index, ElecCombRelIsoEBarmax, ElecCombRelIsoEEndmax);
}

bool UserAnalysisBase::IsIsolatedEl(int index, double ElecCombRelIsoEBarmax, double ElecCombRelIsoEEndmax){
	// combined relative electron isolation for given cut values
	if( fabs(fTR->ElEta[index]) < 1.479 ){ // Barrel
		double combRelIso = ( fTR->ElDR03TkSumPt[index] + std::max(0., fTR->ElDR03EcalRecHitSumEt[index] - 1.) + fTR->ElDR03HcalTowerSumEt[index] ) / (std::max(20.,fTR->ElEt[index]));
		if( combRelIso > ElecCombRelIsoEBarmax ) return false;
	}
	else{ // EndCap
		double combRelIso = ( fTR->ElDR03TkSumPt[index] + fTR->ElDR03EcalRecHitSumEt[index] + fTR->ElDR03HcalTowerSumEt[index] ) / std::max(20.,fTR->ElEt[index]);
		if( combRelIso > ElecCombRelIsoEEndmax ) return false;
	}
	return true;
}

bool UserAnalysisBase::IsTightEl(int index){
	// Definition of "Tight electron" (El.Id cuts: WP80%; El.Convers.Reject.: WP80%; El.RelIso: WP90%)
	if(!IsLooseEl(index)) return false;			// corresponding to the defintion of the "Loose Electron"
	
	if(!IsGoodElId_WP80(index)) return false;	// corresponding to the 80% eff. WP without isolation cuts or conversion rejection
	if(IsConvertedEl_WP80(index)) return false;	
	if(!IsIsolatedEl_WP90(index)) return false;
	return true;
}

bool UserAnalysisBase::IsLooseEl(int index){
	// Definition of "Loose electron" (El.Id cuts: WP90%; El.Convers.Reject.: WP90%; El.RelIso: Loose)
	if(!IsGoodElId_WP90(index)) return false;		// corresponding to the 90% eff. WP without isolation cuts or conversion rejection
	if(IsConvertedEl_WP90(index)) return false;	
	if(!IsIsolatedEl_LooseIso(index)) return false;
	return true;
}

bool UserAnalysisBase::IsLooseNoTightEl(int index){
	if(IsLooseEl(index) && !IsTightEl(index)) return true;
	else return false;
}

bool UserAnalysisBase::IsElFromPrimaryVx(int ind){
	// Returns true if the electron is compatible with the primary vertex (false otherwise)
	
	double distVxmax = 10.;
	double d0Vxmax = 0.02;
	
	// take error from vertex and from track extrapolation into account
	double drVxsq = fTR->PrimVtxxE*fTR->PrimVtxxE + fTR->PrimVtxyE*fTR->PrimVtxyE;
	double d0, dd0, dz, ddz;
	if( ind < 0 || !isGoodBasicPrimaryVertex() ) return false;
	
	d0  = fTR->ElD0PV[ind];
	dd0 = sqrt(fTR->ElD0E[ind]*fTR->ElD0E[ind] + drVxsq);
	dz  = fTR->ElDzPV[ind];
	ddz = sqrt(fTR->ElDzE[ind]*fTR->ElDzE[ind] + fTR->PrimVtxzE*fTR->PrimVtxzE);
	
	// test that the distance is not too large
	if( fabs(d0) > distVxmax * dd0 || fabs(dz) > distVxmax * ddz ) return false;
	
	// test that the distance is within d0Vxmax cut
	if( fabs(d0) > d0Vxmax ) return false;
	
	return true;
}

bool UserAnalysisBase::IsGoodBasicPho(int index){
	// Basic photon cleaning / ID cuts (to be completed)

	double pt	= fTR->PhoPt[index];
	double eta	= fTR->PhoEta[index];

	// barrel cut values
	double PhoHoverEBarmax		= 0.04;
	double PhoSigmaEtaEtaBarmin	= 0.002;
	double PhoSigmaEtaEtaBarmax	= 0.01;
	// endcaps cut values
	double PhoHoverEEndmax		= 0.025;
	double PhoSigmaEtaEtaEndmax	= 0.03;
	// kinematics + ECAL gap cut values
	double etamax				= 3.0;
	double etaEgapUpperEdge		= 1.560;
	double etaEgapLowerEdge		= 1.442;
	double pTmin				= 10.;

	// photon ID
	if( fabs(eta) < 1.479 ){ // Barrel
		if(	fTR->ElHcalOverEcal [index]	> PhoHoverEBarmax      ) return false;
		if(	fTR->ElSigmaIetaIeta[index]	> PhoSigmaEtaEtaBarmax ) return false;
		if(	fTR->ElSigmaIetaIeta[index]	< PhoSigmaEtaEtaBarmin ) return false;
	}
	else{ // EndCap
		if(	fTR->ElHcalOverEcal [index]	> PhoHoverEEndmax      ) return false;
		if(	fTR->ElSigmaIetaIeta[index]	> PhoSigmaEtaEtaEndmax ) return false;
	}

	// photon kinematic cuts + the ECAL gap veto
	if ( (fabs(eta) > etamax) ||
		( (fabs(eta) > etaEgapLowerEdge) && (fabs(eta) < etaEgapUpperEdge) ) )	return false;	
	if( pt < pTmin) return false;

	return true;
}

// Event selections:
bool UserAnalysisBase::isGoodBasicPrimaryVertex(void){
	// Basic primary vertex check
	// (corresponding to TLeptonJEts recommendation)
	
	double	dzVxmax			= 24;
	double	dRVxmax			= 2.;
	int		PrimVtxNdofmin	= 5;
	double	chisqVxmax		= 5.0;
	double	sumPtTkfromVxmin= 30.0;
	
	// Check if the reconstructed primary vertex is fake
	if ( fTR->PrimVtxIsFake )	return false;	
	// Check that there are tracks at the Primary Vertex
	if ( fTR->PrimVtxNdof < PrimVtxNdofmin )	return false;
	// Check the chi2/ndof
	if( fTR->PrimVtxNChi2 > chisqVxmax || fTR->PrimVtxNChi2 < 0. ) return false;
	
	// Check compatibility of vertex with beam spot
	double xVx = fTR->PrimVtxx - fTR->Beamspotx;
	double yVx = fTR->PrimVtxy - fTR->Beamspoty;
	double zVx = fTR->PrimVtxz - fTR->Beamspotz;
	double rVx = sqrt(xVx*xVx + yVx*yVx);	
	if (rVx > dRVxmax || fabs(zVx) > dzVxmax) return false;
	
	// Check that there is sufficient Et in the tracks
	if( fTR->PrimVtxPtSum < sumPtTkfromVxmin ) return false;
	
	return true;
}

bool UserAnalysisBase::IsGoodEvent(){
	// Some cuts on the primary vertex
	double pvx = fTR->PrimVtxx;
	double pvy = fTR->PrimVtxy;
	double pvz = fTR->PrimVtxz;
	if(fTR->PrimVtxIsFake) return false;
	if(fabs(pvz) > 24.) return false;
	if(sqrt(pvx*pvx + pvy*pvy) > 2.0) return false; // Wrt 0,0
	if(fTR->PrimVtxNdof < 5) return false;
	return true;
}

bool UserAnalysisBase::IsGoodMuEvent(){
	if(!IsGoodEvent()) return false;
	bool HLT_Mu9       = GetHLTResult("HLT_Mu9");
	bool HLT_Mu11      = GetHLTResult("HLT_Mu11");
	bool HLT_Mu15      = GetHLTResult("HLT_Mu15");
	bool HLT_DoubleMu0 = GetHLTResult("HLT_DoubleMu0");
	bool HLT_DoubleMu3 = GetHLTResult("HLT_DoubleMu3");

	bool HLT_Jet30U    = GetHLTResult("HLT_Jet30U");
	bool HLT_Jet50U    = GetHLTResult("HLT_Jet50U");
	bool HLT_Jet70U    = GetHLTResult("HLT_Jet70U");
	bool HLT_Jet100U   = GetHLTResult("HLT_Jet100U");
	return (HLT_Mu9 ||
	        HLT_Mu11 ||
	        HLT_Mu15 ||
	        HLT_DoubleMu0 ||
	        HLT_DoubleMu3 ||
	        HLT_Jet30U ||
	        HLT_Jet50U ||
	        HLT_Jet70U ||
	        HLT_Jet100U
	);
}

bool UserAnalysisBase::IsGoodElEvent(){
	int run = fTR->Run;
	// signle-e triggers without ElID or Iso cuts - (should be used for FP ratio measurements)
	bool HLT_Ele10_LW_L1R =				GetHLTResult("HLT_Ele10_LW_L1R");
	//bool HLT_Ele10_SW_L1R =				GetHLTResult("HLT_Ele10_SW_L1R"); // unused at the moment
	bool HLT_Ele15_LW_L1R =				GetHLTResult("HLT_Ele15_LW_L1R");
	//bool HLT_Ele15_SW_L1R =				GetHLTResult("HLT_Ele15_SW_L1R"); // unused at the moment
	// double-e triggers without ElID or Iso cuts - (should be used for FP ratio measurements)
	bool HLT_DoubleEle5_SW_L1R =		GetHLTResult("HLT_DoubleEle5_SW_L1R");
	bool HLT_DoubleEle10_SW_L1R =		GetHLTResult("HLT_DoubleEle10_SW_L1R");
	bool HLT_DoubleEle15_SW_L1R_v1 =	GetHLTResult("HLT_DoubleEle15_SW_L1R_v1");
	// e triggers with ElID or Iso cuts - (should not be used for FP ratio measurements?)
	bool HLT_Ele10_LW_EleId_L1R =		GetHLTResult("HLT_Ele10_LW_EleId_L1R");
	bool HLT_Ele10_SW_EleId_L1R =		GetHLTResult("HLT_Ele10_SW_EleId_L1R");
	bool HLT_Ele15_SW_CaloEleId_L1R =	GetHLTResult("HLT_Ele15_SW_CaloEleId_L1R");
	bool HLT_Ele15_SW_EleId_L1R =		GetHLTResult("HLT_Ele15_SW_EleId_L1R");
	bool HLT_Ele17_SW_LooseEleId_L1R =	GetHLTResult("HLT_Ele17_SW_LooseEleId_L1R");
	bool HLT_Ele17_SW_CaloEleId_L1R =	GetHLTResult("HLT_Ele17_SW_CaloEleId_L1R");
	bool HLT_Ele17_SW_EleId_L1R =		GetHLTResult("HLT_Ele17_SW_EleId_L1R");
	bool HLT_Ele17_SW_TightCaloEleId_SC8HE_L1R_v1 = GetHLTResult("HLT_Ele17_SW_TightCaloEleId_SC8HE_L1R_v1");
	//bool HLT_Ele17_SW_TightEleIdIsol_L1R_v1 = GetHLTResult("HLT_Ele17_SW_TightEleIdIsol_L1R_v1"); // unused at the moment

	if (run==1)						return (HLT_Ele10_LW_L1R);
	if (run>1		&& run<138000)	return (HLT_Ele10_LW_L1R						 || HLT_Ele15_LW_L1R					|| HLT_DoubleEle5_SW_L1R);
	if (run>=138000 && run<=141900)	return (HLT_Ele15_LW_L1R						 || HLT_Ele10_LW_EleId_L1R				|| HLT_DoubleEle5_SW_L1R);
	if (run>141900)					return (HLT_Ele10_SW_EleId_L1R					 || HLT_Ele15_SW_CaloEleId_L1R			|| HLT_Ele15_SW_EleId_L1R ||
											HLT_Ele17_SW_LooseEleId_L1R				 || HLT_Ele17_SW_CaloEleId_L1R			|| HLT_Ele17_SW_EleId_L1R || 
											HLT_Ele17_SW_TightCaloEleId_SC8HE_L1R_v1 || /*HLT_Ele17_SW_TightEleIdIsol_L1R_v1  ||*/
											HLT_DoubleEle10_SW_L1R					 || HLT_DoubleEle15_SW_L1R_v1);
	return false;
}

bool UserAnalysisBase::IsGoodHadronicEvent(){
	// HLT and Jet triggers without ElID or Iso cuts - (should be used for FP ratio measurements)
	bool HLT_Jet30U	= GetHLTResult("HLT_Jet30U");
	bool HLT_Jet50U	= GetHLTResult("HLT_Jet50U");
	bool HLT_Jet70U	= GetHLTResult("fTHLT_Jet70U");
	bool HLT_Jet100U= GetHLTResult("fTHLT_Jet100U");
	bool HLT_HT100U	= GetHLTResult("fTHLT_HT100U");
	bool HLT_HT120U	= GetHLTResult("fTHLT_HT120U");
	bool HLT_HT140U	= GetHLTResult("fTHLT_HT140U");
	bool HLT_HT150U	= GetHLTResult("fTHLT_HT150U");
	
	return (HLT_Jet30U || HLT_Jet50U || HLT_Jet70U || HLT_Jet100U || HLT_HT100U || HLT_HT120U || HLT_HT140U || HLT_HT150U);
}

vector<int> UserAnalysisBase::ElectronSelection(){
	// Returns the vector of indices of
	// good electrons sorted by Pt
	vector<int>    selectedObjInd;
	vector<double> selectedObjPt;
	// form the vector of indices
	for(int ind = 0; ind < fTR->NEles; ++ind){
		// selection
		if(!IsLooseEl(ind)) continue;
		// additional kinematic cuts
		if(fabs(fTR->ElEta[ind]) > 2.4) continue;
		if(fTR->ElPt[ind] < 5.) continue;
		// additional cuts
		if(fTR->ElDuplicateEl[ind] >= 0) continue;
		if(!IsElFromPrimaryVx(ind)) continue;
		
		// rejection if matched to any muon	
		for( int im = 0; im < fTR->NMus; ++im ){
			if( (fTR->MuIsGlobalMuon[im] == 0 || fTR->MuIsTrackerMuon[im] == 0) && ( fTR->MuNTkHits[im] > 10) ) {
				double deltaR = Util::GetDeltaR(fTR->ElEta[ind], fTR->MuEta[im], fTR->ElPhi[ind], fTR->MuPhi[im]);
				if(deltaR <= 0.1) continue;
			}
		}
		
		selectedObjInd.push_back(ind);
		selectedObjPt.push_back(fTR->ElPt[ind]);
	}
	return Util::VSort(selectedObjInd, selectedObjPt);
}

vector<int> UserAnalysisBase::JetSelection(){
	// Returns the vector of indices of
	// good jets sorted by Pt
	vector<int>    selectedObjInd;
	vector<double> selectedObjPt;
	// form the vector of indices
	for(int ind = 0; ind < fTR->NJets; ++ind){
		// selection
		if(!IsGoodBasicJet(ind)) continue;

		selectedObjInd.push_back(ind);
		selectedObjPt.push_back(fTR->JPt[ind]);
	}
	return Util::VSort(selectedObjInd, selectedObjPt);
}

vector<int> UserAnalysisBase::PFJetSelection(){
	// Returns the vector of indices of
	// good jets sorted by Pt
	vector<int>    selectedObjInd;
	vector<double> selectedObjPt;
	// form the vector of indices
	for(int ind = 0; ind < fTR->NJets; ++ind){
		// selection
		if(!IsGoodBasicPFJet(ind)) continue;

		selectedObjInd.push_back(ind);
		selectedObjPt.push_back(fTR->JPt[ind]);
	}
	return Util::VSort(selectedObjInd, selectedObjPt);
}

vector<int> UserAnalysisBase::PhotonSelection(){
	// Returns the vector of indices of
	// good photons sorted by Pt
	vector<int>    selectedObjInd;
	vector<double> selectedObjPt;
	// form the vector of indices
	for(int ind = 0; ind < fTR->NPhotons; ++ind){
		// selection
		if(!IsGoodBasicPho(ind)) continue;
		if(fTR->PhoIsElDupl[ind] >= 0) continue;
		selectedObjInd.push_back(ind);
		selectedObjPt.push_back(fTR->PhoPt[ind]);
	}
	return Util::VSort(selectedObjInd, selectedObjPt);
}

vector<int> UserAnalysisBase::MuonSelection(){
	// Returns the vector of indices of
	// good muons sorted by Pt
	vector<int>	selectedObjInd;
	vector<double>	selectedObjPt;
	// form the vector of indices
	for(int ind = 0; ind < fTR->NMus; ++ind){
		// selection
		if(!IsGoodBasicMu(ind)) continue;
		selectedObjInd.push_back(ind);
		selectedObjPt.push_back(fTR->MuPt[ind]);
	}	
	return Util::VSort(selectedObjInd, selectedObjPt);
}

bool UserAnalysisBase::SingleElectronSelection(int &index){
	// Selects events with (at least) one good electron and gives the index
	// of the hardest one in the argument
	if( fTR->NEles < 1 ) return false;
	vector<int> ElInd = ElectronSelection();
	if(ElInd.size() < 1) return false;
	index = ElInd[0];
	return true;
}

bool UserAnalysisBase::DiElectronSelection(int &ind1, int &ind2){
	// Selects events with (at least) two good electrons and gives the indices
	// of the hardest two in the argument
	if( fTR->NEles < 2 ) return false;
	vector<int> ElInd = ElectronSelection();
	if(ElInd.size() < 2) return false;
	ind1 = ElInd[0];
	ind2 = ElInd[1];
	return true;
}

bool UserAnalysisBase::SSDiElectronSelection(int &prim, int &sec){
	// Select events with 2 SS electrons
	if( fTR->NEles < 2 ) return false;
	vector<int> ElInd = ElectronSelection();
	
	vector<int> priElInd;
	vector<int> secElInd;
	for(unsigned iel = 0; iel < ElInd.size(); ++iel){
		if(fTR->ElPt[ElInd[iel]] > 10.)	priElInd.push_back(ElInd[iel]);
		if(fTR->ElPt[ElInd[iel]] > 10.)	secElInd.push_back(ElInd[iel]);
	}
	
	unsigned npriel = priElInd.size();
	unsigned nsecel = secElInd.size();
	
	// Select events with at least one primary electron and at least one other
	if(npriel < 1 || nsecel < 2) return false;
	
	int ind1 = priElInd[0];
	int ind2 = secElInd[0];
	if(ind2 == ind1) ind2 = secElInd[1];
	// Select same sign electrons
	if(fTR->ElCharge[ind1] != fTR->ElCharge[ind2]) return false;
	prim = ind1;
	sec = ind2;
	return true;
}

bool UserAnalysisBase::SingleMuonSelection(int &index){
	// Selects events with (at least) one good muon and gives the index
	// of the hardest one in the argument
	// Assumes the muons are sorted by pt in the ntuple
	if( fTR->NMus < 1 ) return false;
	vector<int> MuInd;
	for(size_t imu = 0; imu < fTR->NMus; ++imu){
		if(!IsGoodBasicMu(imu)) continue;
		MuInd.push_back(imu);
	}
	if(MuInd.size() < 1) return false;
	
	int maxind = MuInd[0];
	for(size_t i = 0; i < MuInd.size(); ++i){
		int ind = MuInd[i];
		if(fTR->MuPt[ind] > fTR->MuPt[maxind]) maxind = ind;
	}
	index = maxind;
	return true;
}

bool UserAnalysisBase::DiMuonSelection(int &ind1, int &ind2, int charge){
	// Selects events with (at least) two good muons and gives the indices
	// of the hardest two in the argument
	// charge is the relative charge, 0 = no cut on charge, 1 = SS, -1 = OS
	if( fTR->NMus < 2 ) return false;
	vector<int> MuInd;
	for(size_t imu = 0; imu < fTR->NMus; ++imu){
		if(fTR->MuPt[imu] < 10.) continue;
		if(!IsGoodBasicMu(imu)) continue;
		MuInd.push_back(imu);
	}
	if(MuInd.size() < 2) return false;

	// Sort by pt
	int maxind(MuInd[0]), maxind2(MuInd[1]);
	for(size_t i = 1; i < MuInd.size(); ++i){
		int index = MuInd[i];
		if(fTR->MuPt[index] > fTR->MuPt[maxind]){ maxind2 = maxind; maxind = index; }
		else if(fTR->MuPt[index] > fTR->MuPt[maxind2]) maxind2 = index;
	}
	
	ind1 = maxind;
	ind2 = maxind2;
	// Charge selection
	if(charge != 0) if(fTR->MuCharge[ind1] * fTR->MuCharge[ind2] != charge) return false;
	return true;
}

bool UserAnalysisBase::SSDiMuonSelection(int &prim, int &sec){
	// Selects events with (at least) two good muons and gives the indices
	// of the hardest two in the argument
	// charge is the relative charge, 0 = no cut on charge, 1 = SS, -1 = OS
	return DiMuonSelection(prim, sec, 1);
}

///////////////////////////////////////////////////////////////////////////////////////
// Cut Stuff //////////////////////////////////////////////////////////////////////////
void UserAnalysisBase::ReadObjCuts(const char* filename){
// Fills the vectors containing object quality cuts for muons, electrons and jets
	ifstream IN(filename);
	float inf(0.), sup(0.);
	char readbuff[200], branch[200];
	cout << "Reading object selection from " << filename << endl;
	// Loop over lines of datafile
	while( IN.getline(readbuff, 200, '\n') ){
		if (readbuff[0] == '#') {continue;} // Skip lines commented with '#'
		sscanf(readbuff, "%s %f %f", branch, &inf, &sup);
		TBranch *tempbranch = NULL;
		if((tempbranch = fTR->fChain->FindBranch(branch)) == NULL){
			cout << "UserAnalysisBase::ReadObjCuts ==> Branch \"" << branch << "\" not found, continuing... " << endl;
			continue;
		}
		Cut cut;
		if (branch[0] == 'M'){
			cut.branch = tempbranch;
			cut.upperbound = sup;
			cut.lowerbound = inf;
			cout << "  Cutting  " << (cut.branch)->GetName() << "\t\tbetween " << cut.lowerbound << " and " << cut.upperbound << endl;
			fMuCuts.push_back(cut);
		}
		else if (branch[0] == 'E'){
			cut.branch = tempbranch;
			cut.upperbound = sup;
			cut.lowerbound = inf;
			cout << "  Cutting  " << (cut.branch)->GetName() << "\t\tbetween " << cut.lowerbound << " and " << cut.upperbound << endl;
			fElCuts.push_back(cut);
		}
		else if (branch[0] == 'J'){
			cut.branch = tempbranch;
			cut.upperbound = sup;
			cut.lowerbound = inf;
			cout << "  Cutting  " << (cut.branch)->GetName() << "\t\tbetween " << cut.lowerbound << " and " << cut.upperbound << endl;
			fJetCuts.push_back(cut);
		}
	}
	cout << " ----------- " << endl;
}

void UserAnalysisBase::ReadEvtSel(const char* filename){
// Fills the vectors containing event selection cuts
	ifstream IN(filename);
	float inf(0.), sup(0.);
	char readbuff[200], branch[200];
	cout << "Reading event selection from " << filename << endl;
	// Loop over lines of datafile
	while( IN.getline(readbuff, 200, '\n') ){
		if (readbuff[0] == '#') {continue;} // Skip lines commented with '#'
		sscanf(readbuff, "%s %f %f", branch, &inf, &sup);
		TBranch *tempbranch = NULL;
		if((tempbranch = fTR->fChain->FindBranch(branch)) != NULL){
			Cut cut;
			cut.branch = tempbranch;
			cut.upperbound = sup;
			cut.lowerbound = inf;
			cout << "  Cutting  " << (cut.branch)->GetName() << "\t\tbetween " << cut.lowerbound << " and " << cut.upperbound << endl;
			fEvtSelCuts.push_back(cut);
		}
		// else if();
		// Maybe implement here additional, more complicated cuts, to be read from the same file
	}
	cout << " ----------- " << endl;
}

bool UserAnalysisBase::IsGoodObj(int i, vector<Cut> *cutVec){
// This loops on the vector of cuts as filled in ReadObjCuts() and
// determines if an object is good
	TBranch *b;
	TLeaf *l;
	double *fval;
	int *ival;
	for(vector<Cut>::const_iterator it = cutVec->begin(); it != cutVec->end(); it++){
		b = it->branch;
		l = b->GetLeaf(b->GetName());
		if(!strcmp(l->GetTypeName(), "Double_t")){
			fval = (Double_t*)l->GetValuePointer();
			// cout << " Object " << i << " has " << it->branch->GetName() << " " << fval[i] << endl;
			if(fval[i] < it->lowerbound) return false;
			if(fval[i] > it->upperbound) return false;
		}
		if(!strcmp(l->GetTypeName(), "Int_t")){
			ival = (Int_t*)l->GetValuePointer();
			// cout << " Object " << i << " has " << it->branch->GetName() << " " << ival[i] << endl;
			if(ival[i] < it->lowerbound) return false;
			if(ival[i] > it->upperbound) return false;
		}
	}
	return true;
}

bool UserAnalysisBase::IsGoodEvt(vector<Cut> *cutVec){
// This loops on the vector of cuts as filled in ReadEvtSel()
// and determines if an event is selected
	TBranch *b;
	TLeaf *l;
	double *fval;
	int *ival;
		// double *val;
	for(vector<Cut>::const_iterator it = cutVec->begin(); it != cutVec->end(); it++){
		b = it->branch;
		l = b->GetLeaf(b->GetName());
		if(!strcmp(l->GetTypeName(), "Double_t")){
			fval = (double*)l->GetValuePointer();
			// cout << " Object " << i << " has " << it->branch->GetName() << " " << fval[i] << endl;
			if(*fval < it->lowerbound) return false;
			if(*fval > it->upperbound) return false;
		}
		if(!strcmp(l->GetTypeName(), "Int_t")){
			ival = (int*)l->GetValuePointer();
			// cout << " Object " << i << " has " << it->branch->GetName() << " " << ival[i] << endl;
			if(*ival < it->lowerbound) return false;
			if(*ival > it->upperbound) return false;
		}
	}
	return true;
}
