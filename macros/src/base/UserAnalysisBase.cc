#include "base/TreeReader.hh"
#include <stdlib.h>

#include "helper/pdgparticle.hh"
#include "helper/Monitor.hh"
#include "base/UserAnalysisBase.hh"
#include "TH1I.h"
#include "TLorentzVector.h"

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

	fHLTLabels.clear();
	fHLTLabels.reserve(HLTNames->size());
	for( int i=0; i < HLTNames->size(); i++ ){
		fHLTLabelMap[(*HLTNames)[i]] = i; 
		fHLTLabels.push_back((*HLTNames)[i]);
	}
}

int UserAnalysisBase::GetHLTBit(string theHltName){
	if( fHLTLabelMap.empty() ) return -1;
	else{
		map<string,int>::iterator it = fHLTLabelMap.find(theHltName);
		if(it == fHLTLabelMap.end()){
			if(fVerbose > 1) cout << "UserAnalysisBase::GetHLTBit ==> Bit with name " << theHltName << " not found!" << endl;
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
			if(fVerbose > 1) cout << "UserAnalysisBase::GetHLTResult ==> Bit with name " << theHltName << " not found!" << endl;
			return false;
		}
		else return (bool)fTR->HLTResults[bit];
	}
}

void UserAnalysisBase::GetEvtEmChFrac(double & fracEm, double & fracCh){
// Computes the event EM and Charged fractions

	int nMuGood = 0;
	double pt_mu = 0.;
	double pt_track = 0.;
	double et_em = 0.;
	double et_had = 0.;
	for( int i = 0; i < fTR->NMus; ++i ){
		if(fTR->MuGood[i] != 0) continue;
		pt_mu += fTR->MuPt[i];
		pt_track += fTR->MuPt[i];
		nMuGood++;
	}
	for( int i = 0; i < fTR->NEles; ++i ){
		if(fTR->ElGood[i] != 0) continue;
		pt_track += fTR->ElPt[i];
		double em = fTR->ElEt[i];
		if (em < 0.) em = 0.;
		et_em += em;
	}
	for( int i = 0; i < fTR->NPhotons; ++i ){
		if(fTR->PhoGood[i] != 0) continue;
		double em = fTR->PhoPt[i];
		if (em < 0.) em = 0.;
		et_em += em;
	}
	for( int i = 0; i < fTR->NJets; ++i ){
		if(fTR->JGood[i] != 0) continue;
		double pt = fTR->JChfrac[i] * fTR->JPt[i];
		if (pt < 0.) pt = 0.;
		double em = fTR->JEMfrac[i] * fTR->JEt[i];
		if (em < 0.) em = 0.;
		double had = fTR->JEt[i] - em;
		if (had < 0.) had = 0.;
		pt_track += pt;
		et_em    += em;
		et_had   += had;
	}

	fracCh = 0.;
	fracEm = 0.;
	if( et_em + et_had <= 0. ){
		if( nMuGood < 1 ) return;
		fracCh = 1.;
		fracEm = 1.;
	} else {
		fracCh = pt_track / (et_em + et_had + pt_mu);
		fracEm = et_em / (et_em + et_had + pt_mu);
	}
	return;

}

///////////////////////////////////////////////////////////////
// Object selections:
// JETS

bool UserAnalysisBase::IsGoodBasicPFJet(int index, double ptcut, double absetacut){
	// Basic PF jet cleaning and ID cuts
	// cut at pt of ptcut (default = 30 GeV)
	// cut at abs(eta) of absetacut (default = 2.5)
	if(fTR->PFJPt[index] < ptcut           ) return false;
	if(fabs(fTR->PFJEta[index]) > absetacut) return false;
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

bool UserAnalysisBase::IsGoodPFJetMedium(int index, double ptcut, double absetacut) {
	// Medium PF JID
	if ( ! IsGoodBasicPFJet(index, ptcut, absetacut) ) return false;
	if ( !(fTR->PFJNeuHadfrac[index] < 0.95)         ) return false;
	if ( !(fTR->PFJNeuEmfrac[index]  < 0.95)         ) return false;
	return true;
}

bool UserAnalysisBase::IsGoodPFJetTight(int index, double ptcut, double absetacut) {
	// Tight PF JID
	if ( ! IsGoodBasicPFJet(index, ptcut, absetacut) ) return false;
	if ( !(fTR->PFJNeuHadfrac[index] < 0.90)         ) return false;
	if ( !(fTR->PFJNeuEmfrac[index]  < 0.90)         ) return false;
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

// MUONS
bool UserAnalysisBase::IsGoodBasicMu(int index){
	// Basic muon cleaning and ID
	if(fTR->MuIsGlobalMuon[index] == 0)  return false;
	if(fTR->MuIsTrackerMuon[index] == 0) return false;

	if(fTR->MuPt[index] < 5)          return false;
	if(fabs(fTR->MuEta[index]) > 2.4) return false;

	if(fTR->MuNChi2[index] > 10)   return false;
	if(fTR->MuNTkHits[index] < 11) return false;
	if(fTR->MuNMuHits[index] < 1)  return false;

	if(fabs(fTR->MuD0BS[index]) > 0.02)    return false;
	if(fabs(fTR->MuDzPV[index]) > 1.00)    return false;

	if(fTR->MuIso03EMVetoEt[index] > 4.0)  return false;
	if(fTR->MuIso03HadVetoEt[index] > 6.0) return false;

	double iso = fTR->MuRelIso03[index];
	double pt = fTR->MuPt[index];
	double hybiso = iso*pt / std::max(20.,pt);
	if(hybiso > 1.0) return false;

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

// ELECTRONS

bool UserAnalysisBase::IsGoodBasicEl(int index){
	// Electrons with WP95 ID and conv. rej.
	// No isolation cut
	// ECAL gap veto
	if(fTR->ElIDsimpleWP95relIso[index] != 5 && fTR->ElIDsimpleWP95relIso[index] != 7) return false; // Ele ID WP 80
	if(IsElInGap(index)) return false;
	if(!IsElFromPrimaryVx(index)) return false;
	
	return true;
}

bool UserAnalysisBase::IsGoodElId_WP90(int index){
	// Electrons with WP90 ID and conv. rej.
	// No isolation cut
	// ECAL gap veto
	if(!IsGoodBasicEl(index)) return false;
	if(fTR->ElIDsimpleWPrelIso[index] != 5 && fTR->ElIDsimpleWPrelIso[index] != 7) return false;
	return true;
}

bool UserAnalysisBase::IsGoodElId_WP80(int index){
	// Electrons with WP80 ID and conv. rej. cuts
	// No isolation cut
	// ECAL gap veto	
	if(!IsGoodBasicEl(index)) return false;
	if(fTR->ElIDsimpleWP80relIso[index] != 5 && fTR->ElIDsimpleWP80relIso[index] != 7) return false; // Ele ID WP 80
	return true;
}

bool UserAnalysisBase::IsIsolatedEl(int index, double ElecCombRelIsoEBarmax, double ElecCombRelIsoEEndmax){
	// hybrid relative electron isolation for given cut values
	if( fabs(fTR->ElEta[index]) < 1.479 ){ // Barrel
		if( hybRelElIso(index) > ElecCombRelIsoEBarmax ) return false;
	}
	else{ // EndCap
		if( hybRelElIso(index) > ElecCombRelIsoEEndmax ) return false;
	}
	return true;
}

double UserAnalysisBase::hybRelElIso(int index){
	// the value of hybrid relative electron isolation
	if( fabs(fTR->ElEta[index]) < 1.479 ) // Barrel
		return ( fTR->ElDR03TkSumPt[index] + std::max(0., fTR->ElDR03EcalRecHitSumEt[index] - 1.) + fTR->ElDR03HcalTowerSumEt[index] ) / (std::max(20.,fTR->ElEt[index]));	
	else // EndCap
		return ( fTR->ElDR03TkSumPt[index] + fTR->ElDR03EcalRecHitSumEt[index] + fTR->ElDR03HcalTowerSumEt[index] ) / std::max(20.,fTR->ElEt[index]);
}

bool UserAnalysisBase::IsTightEl(int index){
	// Definition of "Tight electron" (El.Id cuts: WP80%; El.Convers.Reject.: WP80%; El.RelIso: WP95%)
	if(!IsLooseEl(index)) return false;			// corresponding to the defintion of the "Loose Electron"

	if(!IsGoodElId_WP80(index)) return false;
	if(!IsIsolatedEl(index, 0.15, 0.15)) return false; // corresponds to WP 95 isolation

	return true;
}

bool UserAnalysisBase::IsLooseEl(int index){
	// Definition of "Loose electron" (reco cuts, El.Id, El.Convers.Reject., El.RelIso)
	if(fTR->ElPt[index] < 10.) return false;
	if(fabs(fTR->ElEta[index]) > 2.4) return false;

	if(!fTR->ElEcalDriven[index]) return false;
	if(fTR->ElCaloEnergy[index] < 10.) return false;
	if(fTR->ElDuplicateEl[index] >= 0) return false;
	
    // (El.Id cuts: WP90%; El.Convers.Reject.: WP80%; El.RelIso: Loose)	
	if(!IsGoodElId_WP90(index)) return false;
	if(fTR->ElNumberOfMissingInnerHits[index] > 0) return false; // Ask for tighter WP80 conversion rejection
	if(!IsIsolatedEl(index, 1.0, 0.6)) return false;
	// rejection if matched to any muon which passes basic quility cuts
	for( int im = 0; im < fTR->NMus; ++im ){
		if( (fTR->MuIsGlobalMuon[im] == 1 && fTR->MuIsTrackerMuon[im] == 1) &&
		    (fTR->MuNTkHits[im] > 10) && (fTR->MuNChi2[im] < 11) && (fTR->MuNMuHits[im] > 0)) {
			double deltaR = Util::GetDeltaR(fTR->ElEta[index], fTR->MuEta[im], fTR->ElPhi[index], fTR->MuPhi[im]);
			if(deltaR <= 0.1) return false;
		}
	}
	return true;
}

bool UserAnalysisBase::IsLooseNoTightEl(int index){
	if(IsLooseEl(index) && !IsTightEl(index)) return true;
	else return false;
}

bool UserAnalysisBase::IsElFromPrimaryVx(int ind){
	// Returns true if the electron is compatible with the primary vertex (false otherwise)
	if(fabs(fTR->ElD0BS[ind]) > 0.02) return false;
	if(fabs(fTR->ElDzPV[ind]) > 1.0) return false;
	return true;
}

bool UserAnalysisBase::IsElInGap(int ind){
	// ECAL gap veto
	if ( fabs(fTR->ElEta[ind]) > 1.442 && fabs(fTR->ElEta[ind]) < 1.560 )  return true;
	return false;
}

// PHOTONS
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

///////////////////////////////////////////////////////////////
// Event selections:
bool UserAnalysisBase::IsGoodEvent(){
	// Some cuts on the primary vertex
	double pvx = fTR->PrimVtxx;
	double pvy = fTR->PrimVtxy;
	double pvz = fTR->PrimVtxz;
	if (fTR->PrimVtxIsFake) return false;
	if (fabs(pvz) > 24.) return false;
	if (sqrt(pvx*pvx + pvy*pvy) > 2.0) return false; // Wrt 0,0
	if (fTR->PrimVtxNdof < 5) return false;
	return true;
}

bool UserAnalysisBase::IsGoodMuEvent(){
	if(!IsGoodEvent()) return false;
	bool HLT_Mu9       = GetHLTResult("HLT_Mu9");
	bool HLT_Mu11      = GetHLTResult("HLT_Mu11");
	bool HLT_Mu13_v1   = GetHLTResult("HLT_Mu13_v1");
	bool HLT_Mu15      = GetHLTResult("HLT_Mu15");
	bool HLT_Mu15_v1   = GetHLTResult("HLT_Mu15_v1");
	bool HLT_DoubleMu0 = GetHLTResult("HLT_DoubleMu0");
	bool HLT_DoubleMu3 = GetHLTResult("HLT_DoubleMu3");
	bool HLT_DoubleMu3_v2 = GetHLTResult("HLT_DoubleMu3_v2");
	bool HLT_DoubleMu5_v2 = GetHLTResult("HLT_DoubleMu5_v2");
	
	return (HLT_Mu9     ||
	        HLT_Mu11    ||
	        HLT_Mu13_v1 ||
	        HLT_Mu15    ||
	        HLT_Mu15_v1 ||
	        HLT_DoubleMu0    ||
	        HLT_DoubleMu3    ||
	        HLT_DoubleMu3_v2 ||
	        HLT_DoubleMu5_v2);
}

bool UserAnalysisBase::IsGoodElEvent(){
	// signle-e triggers without ElID or Iso cuts
	bool HLT_Ele10_LW_L1R =				GetHLTResult("HLT_Ele10_LW_L1R");
	bool HLT_Ele10_SW_L1R =				GetHLTResult("HLT_Ele10_SW_L1R");
	bool HLT_Ele15_LW_L1R =				GetHLTResult("HLT_Ele15_LW_L1R");
	bool HLT_Ele15_SW_L1R =				GetHLTResult("HLT_Ele15_SW_L1R");
	bool HLT_Ele20_SW_L1R =				GetHLTResult("HLT_Ele20_SW_L1R");
	// double-e triggers without ElID or Iso cuts
	bool HLT_DoubleEle5_SW_L1R =		GetHLTResult("HLT_DoubleEle5_SW_L1R");
	bool HLT_DoubleEle10_SW_L1R =		GetHLTResult("HLT_DoubleEle10_SW_L1R");
	bool HLT_DoubleEle15_SW_L1R_v1 =	GetHLTResult("HLT_DoubleEle15_SW_L1R_v1");
	bool HLT_DoubleEle17_SW_L1R_v1 =    GetHLTResult("HLT_DoubleEle17_SW_L1R_v1");
	// e triggers with ElID or Iso cuts
	bool HLT_Ele10_LW_EleId_L1R =		GetHLTResult("HLT_Ele10_LW_EleId_L1R");
	bool HLT_Ele10_SW_EleId_L1R =		GetHLTResult("HLT_Ele10_SW_EleId_L1R");
	bool HLT_Ele15_SW_CaloEleId_L1R =	GetHLTResult("HLT_Ele15_SW_CaloEleId_L1R");
	bool HLT_Ele15_SW_EleId_L1R =		GetHLTResult("HLT_Ele15_SW_EleId_L1R");
	bool HLT_Ele17_SW_LooseEleId_L1R =	GetHLTResult("HLT_Ele17_SW_LooseEleId_L1R");
	bool HLT_Ele17_SW_CaloEleId_L1R =	GetHLTResult("HLT_Ele17_SW_CaloEleId_L1R");
	bool HLT_Ele17_SW_EleId_L1R =		GetHLTResult("HLT_Ele17_SW_EleId_L1R");
	bool HLT_Ele17_SW_TightEleId_L1R =	GetHLTResult("HLT_Ele17_SW_TightEleId_L1R");
	bool HLT_Ele17_SW_TightCaloEleId_SC8HE_L1R_v1 = GetHLTResult("HLT_Ele17_SW_TightCaloEleId_SC8HE_L1R_v1");
	bool HLT_Ele17_SW_TightEleIdIsol_L1R_v1 =       GetHLTResult("HLT_Ele17_SW_TightEleIdIsol_L1R_v1");
	bool HLT_Ele17_SW_TighterEleId_L1R_v1 =         GetHLTResult("HLT_Ele17_SW_TighterEleId_L1R_v1");
	bool HLT_Ele22_SW_TighterEleId_L1R_v2 =         GetHLTResult("HLT_Ele22_SW_TighterEleId_L1R_v2");
	bool HLT_Ele22_SW_TighterEleId_L1R_v3 =         GetHLTResult("HLT_Ele22_SW_TighterEleId_L1R_v3");
	bool HLT_Ele27_SW_TightCaloEleIdTrack_L1R_v1 =  GetHLTResult("HLT_Ele27_SW_TightCaloEleIdTrack_L1R_v1");
	bool HLT_Ele32_SW_TightCaloEleIdTrack_L1R_v1 =  GetHLTResult("HLT_Ele32_SW_TightCaloEleIdTrack_L1R_v1");
	bool HLT_Ele32_SW_TighterEleId_L1R_v2 =         GetHLTResult("HLT_Ele32_SW_TighterEleId_L1R_v2");

	return (HLT_Ele10_LW_L1R				|| HLT_Ele15_LW_L1R            || HLT_DoubleEle5_SW_L1R				    ||
			HLT_Ele10_SW_L1R				|| HLT_Ele15_SW_L1R            || HLT_Ele20_SW_L1R                      ||
			HLT_Ele10_LW_EleId_L1R			|| HLT_Ele10_SW_EleId_L1R      || HLT_Ele15_SW_CaloEleId_L1R			|| HLT_Ele15_SW_EleId_L1R ||
            HLT_Ele17_SW_LooseEleId_L1R     || HLT_Ele17_SW_CaloEleId_L1R  || HLT_Ele17_SW_TightCaloEleId_SC8HE_L1R_v1 ||
			HLT_Ele17_SW_EleId_L1R          || HLT_Ele17_SW_TightEleId_L1R || HLT_Ele17_SW_TightEleIdIsol_L1R_v1    || HLT_Ele17_SW_TighterEleId_L1R_v1     || 
			HLT_Ele20_SW_L1R                ||
            HLT_Ele22_SW_TighterEleId_L1R_v2        || HLT_Ele22_SW_TighterEleId_L1R_v3 || HLT_Ele27_SW_TightCaloEleIdTrack_L1R_v1 ||
			HLT_Ele32_SW_TightCaloEleIdTrack_L1R_v1 || HLT_Ele32_SW_TighterEleId_L1R_v2 ||
			HLT_DoubleEle10_SW_L1R                  || HLT_DoubleEle15_SW_L1R_v1        || HLT_DoubleEle17_SW_L1R_v1);
}

bool UserAnalysisBase::IsGoodElEvent_RA5(){
	int run = fTR->Run;
	// signle-e triggers without ElID or Iso cuts
	bool HLT_Ele10_LW_L1R =          GetHLTResult("HLT_Ele10_LW_L1R");
	bool HLT_Ele10_SW_L1R =          GetHLTResult("HLT_Ele10_SW_L1R");
	bool HLT_Ele15_LW_L1R =          GetHLTResult("HLT_Ele15_LW_L1R");
	bool HLT_Ele15_SW_L1R =          GetHLTResult("HLT_Ele15_SW_L1R");
	bool HLT_Ele20_SW_L1R =          GetHLTResult("HLT_Ele20_SW_L1R");
	// double-e triggers without ElID or Iso cuts
	bool HLT_DoubleEle5_SW_L1R =     GetHLTResult("HLT_DoubleEle5_SW_L1R");
	bool HLT_DoubleEle10_SW_L1R =    GetHLTResult("HLT_DoubleEle10_SW_L1R");
	bool HLT_DoubleEle15_SW_L1R_v1 = GetHLTResult("HLT_DoubleEle15_SW_L1R_v1");
	bool HLT_DoubleEle17_SW_L1R_v1 = GetHLTResult("HLT_DoubleEle17_SW_L1R_v1");
	// e triggers with ElID or Iso cuts
	bool HLT_Ele10_LW_EleId_L1R =		GetHLTResult("HLT_Ele10_LW_EleId_L1R");
	bool HLT_Ele10_SW_EleId_L1R =		GetHLTResult("HLT_Ele10_SW_EleId_L1R");
	bool HLT_Ele15_SW_CaloEleId_L1R =	GetHLTResult("HLT_Ele15_SW_CaloEleId_L1R");
	bool HLT_Ele15_SW_EleId_L1R =		GetHLTResult("HLT_Ele15_SW_EleId_L1R");
	bool HLT_Ele17_SW_LooseEleId_L1R =	GetHLTResult("HLT_Ele17_SW_LooseEleId_L1R");
	bool HLT_Ele17_SW_CaloEleId_L1R =	GetHLTResult("HLT_Ele17_SW_CaloEleId_L1R");
	bool HLT_Ele17_SW_EleId_L1R =		GetHLTResult("HLT_Ele17_SW_EleId_L1R");
	bool HLT_Ele17_SW_TightEleId_L1R =  GetHLTResult("HLT_Ele17_SW_TightEleId_L1R");
	bool HLT_Ele17_SW_TighterEleId_L1R_v1 =         GetHLTResult("HLT_Ele17_SW_TighterEleId_L1R_v1");
	bool HLT_Ele22_SW_TighterEleId_L1R_v2 =         GetHLTResult("HLT_Ele22_SW_TighterEleId_L1R_v2");
	bool HLT_Ele22_SW_TighterEleId_L1R_v3 =         GetHLTResult("HLT_Ele22_SW_TighterEleId_L1R_v3");
	bool HLT_Ele27_SW_TightCaloEleIdTrack_L1R_v1 =  GetHLTResult("HLT_Ele27_SW_TightCaloEleIdTrack_L1R_v1");
	bool HLT_Ele32_SW_TightCaloEleIdTrack_L1R_v1 =  GetHLTResult("HLT_Ele32_SW_TightCaloEleIdTrack_L1R_v1");
	bool HLT_Ele32_SW_TighterEleId_L1R_v2 =         GetHLTResult("HLT_Ele32_SW_TighterEleId_L1R_v2");
	
	if (run==1)						return (HLT_Ele10_LW_L1R);
	if (run>1		&& run<138000)	return (HLT_Ele10_LW_L1R            || HLT_Ele10_SW_L1R           || HLT_Ele15_LW_L1R           || HLT_DoubleEle5_SW_L1R);
	if (run>=138000 && run<=141900)	return (HLT_Ele15_LW_L1R            || HLT_Ele15_SW_L1R           || HLT_Ele10_LW_EleId_L1R     || HLT_DoubleEle5_SW_L1R);
	if (run>141900)					return (HLT_Ele10_SW_EleId_L1R      || HLT_Ele15_SW_CaloEleId_L1R || HLT_Ele15_SW_EleId_L1R ||
											HLT_Ele17_SW_LooseEleId_L1R || HLT_Ele17_SW_CaloEleId_L1R || HLT_Ele17_SW_EleId_L1R || 
											HLT_Ele17_SW_TightEleId_L1R || HLT_Ele17_SW_TighterEleId_L1R_v1 || HLT_Ele20_SW_L1R ||
											HLT_Ele22_SW_TighterEleId_L1R_v2        || HLT_Ele22_SW_TighterEleId_L1R_v3 || HLT_Ele27_SW_TightCaloEleIdTrack_L1R_v1 ||
											HLT_Ele32_SW_TightCaloEleIdTrack_L1R_v1 || HLT_Ele32_SW_TighterEleId_L1R_v2 ||
											HLT_DoubleEle10_SW_L1R      || HLT_DoubleEle15_SW_L1R_v1  || HLT_DoubleEle17_SW_L1R_v1);
	return false;
}

bool UserAnalysisBase::IsGoodElEvent_TDL(){
	int run = fTR->Run;
	// signle-e triggers without ElID or Iso cuts
	bool HLT_Ele10_LW_L1R =				GetHLTResult("HLT_Ele10_LW_L1R");
	bool HLT_Ele15_LW_L1R =				GetHLTResult("HLT_Ele15_LW_L1R");
	bool HLT_Ele15_SW_L1R =				GetHLTResult("HLT_Ele15_SW_L1R");
	bool HLT_Ele20_SW_L1R =				GetHLTResult("HLT_Ele20_SW_L1R");
	// double-e triggers without ElID or Iso cuts
	bool HLT_DoubleEle10_SW_L1R =		GetHLTResult("HLT_DoubleEle10_SW_L1R");
	bool HLT_DoubleEle15_SW_L1R_v1 =	GetHLTResult("HLT_DoubleEle15_SW_L1R_v1");
	// e triggers with ElID or Iso cuts
	bool HLT_Ele15_SW_CaloEleId_L1R =	GetHLTResult("HLT_Ele15_SW_CaloEleId_L1R");
	bool HLT_Ele17_SW_CaloEleId_L1R =	GetHLTResult("HLT_Ele17_SW_CaloEleId_L1R");
	bool HLT_Ele17_SW_TightEleId_L1R =	GetHLTResult("HLT_Ele17_SW_TightEleId_L1R");
	bool HLT_Ele17_SW_TightCaloEleId_SC8HE_L1R_v1 = GetHLTResult("HLT_Ele17_SW_TightCaloEleId_SC8HE_L1R_v1");

	if (run==1)						return (HLT_Ele10_LW_L1R);
	if (run>1		&& run<138000)	return (HLT_Ele10_LW_L1R);
	if (run>=138000 && run<141900)	return (HLT_Ele15_LW_L1R);
	if (run>=141900 && run<144000)	return (HLT_Ele15_SW_L1R);
	if (run>=144000 && run<144114)	return (HLT_Ele15_SW_CaloEleId_L1R	|| HLT_Ele20_SW_L1R							|| HLT_DoubleEle10_SW_L1R );
	if (run>=146000 && run<147120)	return (HLT_DoubleEle10_SW_L1R		|| HLT_Ele17_SW_CaloEleId_L1R);
	if (run>=147120)				return (HLT_DoubleEle15_SW_L1R_v1	|| HLT_Ele17_SW_TightCaloEleId_SC8HE_L1R_v1	|| HLT_Ele17_SW_TightEleId_L1R);
	return false;
}

bool UserAnalysisBase::IsGoodElFakesEvent(){
	// signle-e triggers without ElID or Iso cuts to be used for FP ratio measurements
	bool HLT_Ele10_LW_L1R =	GetHLTResult("HLT_Ele10_LW_L1R");
	bool HLT_Ele10_SW_L1R =	GetHLTResult("HLT_Ele10_SW_L1R");
	bool HLT_Ele15_LW_L1R =	GetHLTResult("HLT_Ele15_LW_L1R");
	bool HLT_Ele15_SW_L1R =	GetHLTResult("HLT_Ele15_SW_L1R");
	
	return (HLT_Ele10_LW_L1R || HLT_Ele10_SW_L1R || HLT_Ele15_LW_L1R || HLT_Ele15_SW_L1R);
}

bool UserAnalysisBase::IsGoodHadronicEvent(){
	// HLT and Jet triggers without ElID or Iso cuts - (should be used for FP ratio measurements)
	bool HLT_Jet30U	= GetHLTResult("HLT_Jet30U");
	bool HLT_Jet50U	= GetHLTResult("HLT_Jet50U");
	bool HLT_Jet70U	= GetHLTResult("HLT_Jet70U");
	bool HLT_Jet100U= GetHLTResult("HLT_Jet100U");
	bool HLT_Jet100U_v2= GetHLTResult("HLT_Jet100U_v2");
	bool HLT_Jet100U_v3= GetHLTResult("HLT_Jet100U_v3");
	bool HLT_HT100U	= GetHLTResult("HLT_HT100U");
	bool HLT_HT120U	= GetHLTResult("HLT_HT120U");
	bool HLT_HT130U	= GetHLTResult("HLT_HT130U");
	bool HLT_HT140U	= GetHLTResult("HLT_HT140U");
	bool HLT_HT150U	= GetHLTResult("HLT_HT150U");
	bool HLT_HT150U_v3	        = GetHLTResult("HLT_HT150U_v3");
	bool HLT_HT200U	            = GetHLTResult("HLT_HT200U");
	// HLT cross-triggers
	bool HLT_Mu5_HT70U_v1	    = GetHLTResult("HLT_Mu5_HT70U_v1");
	bool HLT_Mu5_HT100U_v1	    = GetHLTResult("HLT_Mu5_HT100U_v1");
	bool HLT_DoubleMu3_HT50U    = GetHLTResult("HLT_DoubleMu3_HT50U");
	bool HLT_Ele10_HT70U	    = GetHLTResult("HLT_Ele10_HT70U");
	bool HLT_Ele10_HT100U	    = GetHLTResult("HLT_Ele10_HT100U");
	bool HLT_Ele10_EleId_HT70U	= GetHLTResult("HLT_Ele10_EleId_HT70U");
	bool HLT_DoubleEle8_SW_HT70U_LR1_v1	= GetHLTResult("HLT_DoubleEle8_SW_HT70U_LR1_v1");
	bool HLT_Mu5_Ele5	        = GetHLTResult("HLT_Mu5_Ele5");
	bool HLT_Mu3_Ele8_HT70U_v1	= GetHLTResult("HLT_Mu3_Ele8_HT70U_v1");

	return (HLT_Jet30U       || HLT_Jet50U        || HLT_Jet70U            || HLT_Jet100U || HLT_Jet100U_v2 || HLT_Jet100U_v3 ||
			HLT_HT100U       || HLT_HT120U        || HLT_HT130U            || HLT_HT140U  || HLT_HT150U || HLT_HT200U || HLT_HT150U_v3 ||
			HLT_Mu5_HT70U_v1 || HLT_Mu5_HT100U_v1 || HLT_DoubleMu3_HT50U   ||
			HLT_Ele10_HT70U  || HLT_Ele10_HT100U  || HLT_Ele10_EleId_HT70U || HLT_DoubleEle8_SW_HT70U_LR1_v1 ||
			HLT_Mu5_Ele5     || HLT_Mu3_Ele8_HT70U_v1);
}

vector<int> UserAnalysisBase::ElectronSelection(bool(UserAnalysisBase::*eleSelector)(int)){
	// Returns the vector of indices of
	// good electrons sorted by Pt
	if(eleSelector == NULL) eleSelector = &UserAnalysisBase::IsGoodBasicEl;
	vector<int>    selectedObjInd;
	vector<double> selectedObjPt;
	// form the vector of indices
	for(int ind = 0; ind < fTR->NEles; ++ind){
		// selection
		if((*this.*eleSelector)(ind) == false) continue;

		selectedObjInd.push_back(ind);
		selectedObjPt.push_back(fTR->ElPt[ind]);
	}
	return Util::VSort(selectedObjInd, selectedObjPt);
}

vector<int> UserAnalysisBase::JetSelection(bool(UserAnalysisBase::*jetSelector)(int)){
	// Returns the vector of indices of
	// good jets sorted by Pt
	vector<int>    selectedObjInd;
	vector<double> selectedObjPt;
	if(jetSelector == NULL) jetSelector = &UserAnalysisBase::IsGoodBasicJet;
	// form the vector of indices
	for(int ind = 0; ind < fTR->NJets; ++ind){
		// selection
		if((*this.*jetSelector)(ind) == false) continue;

		// additional kinematic cuts
		if(fabs(fTR->JEta[ind]) > 2.5) continue;
		
		selectedObjInd.push_back(ind);
		selectedObjPt.push_back(fTR->JPt[ind]);
	}
	return Util::VSort(selectedObjInd, selectedObjPt);
}

vector<int> UserAnalysisBase::PFJetSelection(double ptcut, double absetacut, bool(UserAnalysisBase::*pfjetSelector)(int, double, double)){
	// Returns the vector of indices of
	// good jets sorted by Pt
	// cut at pt of ptcut (default = 30.)
	// cut at abs(eta) of absetacut (default = 2.5)
	if(pfjetSelector == NULL) pfjetSelector = &UserAnalysisBase::IsGoodBasicPFJet;
	vector<int>    selectedObjInd;
	vector<double> selectedObjPt;
	// form the vector of indices
	for(int ind = 0; ind < fTR->PFNJets; ++ind){
		// selection
		if((*this.*pfjetSelector)(ind, ptcut, absetacut) == false) continue;

		selectedObjInd.push_back(ind);
		selectedObjPt.push_back(fTR->PFJPt[ind]);
	}
	return Util::VSort(selectedObjInd, selectedObjPt);
}

vector<int> UserAnalysisBase::PhotonSelection(bool(UserAnalysisBase::*phoSelector)(int)){
	// Returns the vector of indices of
	// good photons sorted by Pt
	if(phoSelector == NULL) phoSelector = &UserAnalysisBase::IsGoodBasicPho;
	vector<int>    selectedObjInd;
	vector<double> selectedObjPt;
	// form the vector of indices
	for(int ind = 0; ind < fTR->NPhotons; ++ind){
		// selection
		if((*this.*phoSelector)(ind) == false) continue;
		if(fTR->PhoIsElDupl[ind] >= 0) continue;
		selectedObjInd.push_back(ind);
		selectedObjPt.push_back(fTR->PhoPt[ind]);
	}
	return Util::VSort(selectedObjInd, selectedObjPt);
}

vector<int> UserAnalysisBase::MuonSelection(bool(UserAnalysisBase::*muonSelector)(int)){
	// Returns the vector of indices of
	// good muons sorted by Pt
	if(muonSelector == NULL) muonSelector = &UserAnalysisBase::IsGoodBasicMu;
	vector<int>	selectedObjInd;
	vector<double>	selectedObjPt;
	// form the vector of indices
	for(int ind = 0; ind < fTR->NMus; ++ind){
		// selection
		if((*this.*muonSelector)(ind) == false) continue;
		selectedObjInd.push_back(ind);
		selectedObjPt.push_back(fTR->MuPt[ind]);
	}	
	return Util::VSort(selectedObjInd, selectedObjPt);
}

bool UserAnalysisBase::SingleElectronSelection(int &index, bool(UserAnalysisBase::*eleSelector)(int)){
	// Selects events with (at least) one good electron and gives the index
	// of the hardest one in the argument
	if(eleSelector == NULL) eleSelector = &UserAnalysisBase::IsGoodBasicEl;
	if( fTR->NEles < 1 ) return false;
	vector<int> ElInd = ElectronSelection(eleSelector);
	if(ElInd.size() < 1) return false;
	index = ElInd[0];
	return true;
}

bool UserAnalysisBase::DiElectronSelection(int &ind1, int &ind2, int charge, bool(UserAnalysisBase::*eleSelector)(int)){
	// Selects events with (at least) two good electrons and gives the indices
	// of the hardest two in the argument (if selected)
	// charge is the relative charge, 0 = no cut on charge (default), 1 = SS, -1 = OS
	if(eleSelector == NULL) eleSelector = &UserAnalysisBase::IsGoodBasicEl;
	if( fTR->NEles < 2 ) return false;
	vector<int> ElInd = ElectronSelection(eleSelector);
	if(ElInd.size() < 2) return false;
	// Charge selection
	if(charge != 0) if(fTR->ElCharge[ElInd[0]] * fTR->ElCharge[ElInd[1]] != charge) return false;
	ind1 = ElInd[0];
	ind2 = ElInd[1];
	return true;
}

bool UserAnalysisBase::SSDiElectronSelection(int &prim, int &sec, bool(UserAnalysisBase::*eleSelector)(int)){
	if(eleSelector == NULL) eleSelector = &UserAnalysisBase::IsGoodBasicEl;
	return DiElectronSelection(prim, sec, 1, eleSelector);
}

bool UserAnalysisBase::OSDiElectronSelection(int &prim, int &sec, bool(UserAnalysisBase::*eleSelector)(int)){
	if(eleSelector == NULL) eleSelector = &UserAnalysisBase::IsGoodBasicEl;
	return DiElectronSelection(prim, sec, -1, eleSelector);
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

///////////////////////////////////////////////////////////////
// Cut Stuff
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

void UserAnalysisBase::EventPrint(){ 

  char ccharge[] = {'-', '0', '+'};
  TLorentzVector p1(0.,0.,0.,0.), p2(0.,0.,0.,0.), psum(0.,0.,0.,0.);
  double minv = -999.99;

  double etaBar = 2.6;
  double muIsomax = 0.15;
  double elIsomax = 0.15;
  double phoIsomax = 0.15;
  //  double METPhi = fTR->MuJESCorrMETphi;
  //  double MET = fTR->MuJESCorrMET;
  double METPhi = fTR->PFMETphi;
  double MET = fTR->PFMET;

  double METsign = fTR->RawMETSignificance;

  // print the event and primary vertex info
  TString tag = "";
  if (fTR->NJets == 1 && fabs(fTR->JEta[0]) < etaBar) tag = "***";
  if (fTR->NJets >= 2 && (fabs(fTR->JEta[0]) < etaBar || fabs(fTR->JEta[1] < etaBar))) tag = "***";
  cout << endl;
  cout << " Run = " << fTR->Run << ", Lumi sect = " << fTR->LumiSection
       << ", Event = " << fTR->Event << "  " << tag << endl;
  cout << " Number of vertices = " << fTR->NVrtx << endl;
  //  if (fTR->GoodEvent > 0) cout << " -> bad event, tag = " << fTR->GoodEvent << endl;
  if (fTR->HBHENoiseFlag == 0) cout << " -> Noisy HCAL " << endl;
  if (fTR->PrimVtxGood > 0) cout << " -> bad PrimVtx, tag = " << fTR->PrimVtxGood << endl;
  cout << " PrimVtx     x = " << fTR->PrimVtxx << ", y = " << fTR->PrimVtxy << ", z = " << fTR->PrimVtxz << endl;
  cout << " PrimVtx, Ndof = " << fTR->PrimVtxNdof << ", chisq = " << fTR->PrimVtxNChi2 
       << ", TkPtSum = " << fTR->PrimVtxPtSum << endl;
  double fracEm, fracCh;
  GetEvtEmChFrac(fracEm, fracCh);
  cout << " Event fEM = " << fracEm << ", fCh = " << fracCh << endl;
  cout << " caloMET  = " << MET    << ", METPhi = " << METPhi << ", METsignif = " << METsign << endl;
  cout << " TCMET    = " << fTR->TCMET    << ", METPhi = " << fTR->TCMETphi << endl;
  cout << " PFMET    = " << fTR->PFMET    << ", METPhi = " << fTR->PFMETphi << endl;
		
  // print the jets info
  cout << " Number of jets in the ntuple = " << fTR->PFNJets
    //       << ", total number of jets = " << fTR->NJetsTot
       << endl;
  int nJetsCand = 0;
  for (int i = 0; i < fTR->PFNJets; ++i) {
    double dPhiMJ = Util::DeltaPhi(fTR->PFJPhi[i], METPhi);
    cout << " Jet" << i << " Pt = " << fTR->PFJPt[i] << ", Phi = " << fTR->PFJPhi[i]
	 << ", Eta = " << fTR->PFJEta[i] << ", E = " << fTR->PFJE[i] << endl;
    double jmass = sqrt(fTR->PFJE[i]*fTR->PFJE[i]-fTR->PFJPx[i]*fTR->PFJPx[i]-fTR->PFJPy[i]*fTR->PFJPy[i]-fTR->PFJPz[i]*fTR->PFJPz[i]);
    cout << "      " << " Px = " << fTR->PFJPx[i] << " Py = " << fTR->PFJPy[i] << " Pz = " << fTR->PFJPz[i]
	 << " Jet mass = " << jmass << endl;
    //		cout << "      " << " dPhiJM = " << dPhiMJ << ", JetNtrk = " << fTR->JNAssoTracks[i] << endl;
    //		cout << "      " << " JetfEM = " << fTR->JEMfrac[i] << ", JetfCh = " << fTR->JChfrac[i] << endl;
    cout << "      " << " dPhiJM = " << dPhiMJ << ", JetNtrk = " << fTR->PFJChMult[i] << endl;
    cout << "      " << " JetfEMch = " << fTR->PFJChEmfrac[i] << ", JetfEMneu = " << fTR->PFJNeuEmfrac[i]
	 << " JetfHach = " << fTR->PFJChHadfrac[i] << ", JetfHaneu = " << fTR->PFJNeuHadfrac[i] << endl;
    //    cout << "      " << " Jetn90 = " << fTR->JID_n90Hits[i] << ", JetHPD = " << fTR->JID_HPD[i] 
    //	 << ", JetRBX = " << fTR->JID_RBX[i] << endl;
    //    cout << "      " << " Jet Vtx Chisq/ndof = " << fTR->JVtxNChi2[i] << endl;
    cout << "      " << " b-tag TCHE = " << fTR->PFJbTagProbTkCntHighEff[i]<< ", TCHP = " << fTR->PFJbTagProbTkCntHighPur[i] << endl;
    cout << "      " << " b-tag SSVHE = " << fTR->PFJbTagProbSimpSVHighEff[i]<< ", SSVHP = " << fTR->PFJbTagProbSimpSVHighPur[i] << endl;
    //    if (fTR->JGood[i] > 0) cout << "      -> bad jet, tag = " << fTR->JGood[i] << endl;
    //    if (fTR->JGood[i] == 0)
    nJetsCand++;
    if (fTR->PFJPt[i] > 20. && fabs(dPhiMJ-3.141592654) < 0.05) {
      cout << "      -> jet back-to-back with MET" << endl;
    }
  }

  // for multi-jet events
  if (fTR->PFNJets >= 2) {
    double dPhiMJ1 = Util::DeltaPhi(fTR->PFJPhi[0], METPhi);
    double dPhiMJ2 = Util::DeltaPhi(fTR->PFJPhi[1], METPhi);
    double R12 = sqrt(dPhiMJ1*dPhiMJ1 + (TMath::Pi()-dPhiMJ2)*(TMath::Pi()-dPhiMJ2) );
    double R21 = sqrt(dPhiMJ2*dPhiMJ2 + (TMath::Pi()-dPhiMJ1)*(TMath::Pi()-dPhiMJ1) );
    double dPhij12 = Util::DeltaPhi(fTR->PFJPhi[0], fTR->PFJPhi[1]);
    double dRj12 = Util::GetDeltaR(fTR->PFJEta[0], fTR->PFJEta[1], fTR->PFJPhi[0], fTR->PFJPhi[1]);
    cout << " R12    = " << R12 << ", R21    = " << R21 
	 << ", dRj12  = " << dRj12 << ", dPhij12 = " << dPhij12 << endl;
  }
		
  // print muon info
  for (int i = 0; i < fTR->NMus; ++i) {
    cout << " Muon" << i  << " " << ccharge[fTR->MuCharge[i]+1] << " Pt = " << fTR->MuPt[i] << ", Phi = " << fTR->MuPhi[i]
	 << ", Eta = " << fTR->MuEta[i] << ", E = " << fTR->MuE[i] << endl;
    cout << "      " << " Px = " << fTR->MuPx[i] << " Py = " << fTR->MuPy[i] << " Pz = " << fTR->MuPz[i] << endl;
    TString muQual = "";
    if (fTR->MuIsTrackerMuon[i]) muQual = muQual + "TrackerMuon ";
    if (fTR->MuIsGlobalMuon[i]) muQual = muQual + "GlobalMuon ";
    if (muQual != "") cout << "       " << muQual << endl;
    cout << "      " << " doPvx = " << fTR->MuD0PV[i] << " (signif. = " << fTR->MuD0PV[i]/fTR->MuD0E[i] << ")"
	 << " dzPvx = " << fTR->MuDzPV[i] << " (signif. = " << fTR->MuDzPV[i]/fTR->MuDzE[i] << ")" << endl;
    cout << "      " << " chisq = " << fTR->MuNChi2[i] << " #Trk hits = " << fTR->MuNTkHits[i] << endl;
    cout << "      " << " IsoVal03 = " << fTR->MuRelIso03[i] << endl;
    cout << "      " << " IsoVal Trk = " << fTR->MuIso03SumPt[i] << " Em = " << fTR->MuIso03EmEt[i] << " Had = " << fTR->MuIso03HadEt[i] << endl;
    cout << "      " << " Energy assoc to the muon Em = " << fTR->MuEem[i] << " Had = " << fTR->MuEhad[i] << endl;
    //		if (fTR->MuGood[i] > 0) cout << "      -> bad muon, tag = " << fTR->MuGood[i] << endl;
    //		if (fTR->MuIsIso[i] == 1) cout << "      -> muon is isolated " << endl;
    if (fTR->MuRelIso03[i] <= muIsomax) cout << "      -> muon is isolated " << endl;
    else cout << "      -> muon NOT isolated " << endl;
  }
  
  // print electron info
  for (int i = 0; i < fTR->NEles; ++i) {
    cout << " Elec" << i  << " " << ccharge[fTR->ElCharge[i]+1] << " Pt = " << fTR->ElPt[i] << ", Phi = " << fTR->ElPhi[i]
	 << ", Eta = " << fTR->ElEta[i] << ", E = " << fTR->ElE[i] << endl;
    cout << "      " << " Px = " << fTR->ElPx[i] << " Py = " << fTR->ElPy[i] << " Pz = " << fTR->ElPz[i] << endl;
    TString elQual = "";
    if (fTR->ElTrackerDriven[i]) elQual = elQual + "TrackerDriven ";
    if (fTR->ElEcalDriven[i]) elQual = elQual + "EcalDriven ";
    if (elQual != "") cout << "       " << elQual << endl;
    cout << "      " << " doPvx = " << fTR->ElD0PV[i] << " (signif. = " << fTR->ElD0PV[i]/fTR->ElD0E[i] << ")"
	 << " dzPvx = " << fTR->ElDzPV[i] << " (signif. = " << fTR->ElDzPV[i]/fTR->ElDzE[i] << ")" << endl;
    cout << "      " << " HoverE = " << fTR->ElHcalOverEcal[i] << " sigmaIetaIeta = " << fTR->ElSigmaIetaIeta[i] << endl;
    cout << "      " << " DeltaEtain = " << fTR->ElDeltaEtaSuperClusterAtVtx[i]
	 << " DeltaPhiin = " << fTR->ElDeltaPhiSuperClusterAtVtx[i] << endl;
    cout << "      " << " Brem fraction = " << fTR->Elfbrem[i] << endl;
    cout << "      " << " Conversion: missed inner hits = " << fTR->ElNumberOfMissingInnerHits[i]
	 << " dist = " << fTR->ElConvPartnerTrkDist[i] << " cotTheta = " << fTR->ElConvPartnerTrkDCot[i]
	 << " partner charge = " << ccharge[(int)fTR->ElConvPartnerTrkCharge[i]+1] << endl;
    cout << "      " << " GSF-CTF charge consistency = " << fTR->ElCInfoIsGsfCtfCons[i] << " Sc charge = " << fTR->ElScPixCharge[i] << endl;
    //    if (fTR->ElIsInJet[i] >= 0) cout  << "      " << " Is in jet nber = " << fTR->ElIsInJet[i]
    //				      << " sharedE = " << fTR->ElSharedEnergy[i] << endl;
    cout << "      " << " IsoVal03 = " << fTR->ElRelIso03[i] << endl;
    cout << "      " << " IsoVal Trk = " << fTR->ElDR03TkSumPt[i] << " Em = " << fTR->ElDR03EcalRecHitSumEt[i] 
	 << " Had = " << fTR->ElDR03HcalTowerSumEt[i] << endl;
    //		if (fTR->ElGood[i] > 0) cout << "      -> bad electron, tag = " << fTR->ElGood[i] << endl;
    //		if (fTR->ElIsIso[i] == 1) cout << "      -> electron is isolated " << endl;
    if (fTR->ElRelIso03[i] <= elIsomax) cout << "      -> electron is isolated " << endl;
    else cout << "      -> electron NOT isolated " << endl;
  }
		
  // print photon info
  for (int i = 0; i < fTR->NPhotons; ++i) {
    cout << " Phot" << i << " Pt = " << fTR->PhoPt[i] << ", Phi = " << fTR->PhoPhi[i]
		                << ", Eta = " << fTR->PhoEta[i] << ", E = " << fTR->PhoEnergy[i] << endl;
    cout << "      " << " Px = " << fTR->PhoPx[i] << " Py = " << fTR->PhoPy[i] << " Pz = " << fTR->PhoPz[i] << endl;
    cout << "      " << " HoverE = " << fTR->PhoHoverE[i] << " sigmaEtaEta = " << fTR->PhoSigmaIetaIeta[i] << endl;
    //    if (fTR->PhoIsInJet[i] >= 0) cout << "      "  << " Is in jet nber = " << fTR->PhoIsInJet[i]
    //				      << " sharedE = " << fTR->PhoSharedEnergy[i] << endl;
    cout << "      " << " IsoVal = " << fTR->PhoIso03[i] << endl;
    cout << "      " << " IsoVal Trk = " << fTR->PhoIso03TrkSolid[i] << " Em = " << fTR->PhoIso03Ecal[i] 
	 << " Had = " << fTR->PhoIso03Hcal[i] << endl;
    cout << "      " << " Has conversion tracks = " << fTR->PhoHasConvTrks[i] << endl;
    if (fTR->PhoIsElDupl[i] >= 0) cout << "      " << " Photon duplic with Elec" << fTR->PhoIsElDupl[i] << endl;
    //		if (fTR->PhoGood[i] > 0) cout << "      -> bad photon, tag = " << fTR->PhoGood[i] << endl;
    //		if (fTR->PhoIsIso[i] == 1) cout << "      -> photon is isolated " << endl;
    if (fTR->PhoIso03[i] <= phoIsomax) cout << "      -> photon is isolated " << endl;
    else cout << "      -> photon NOT isolated " << endl;
  }

  // jets/ lepton/ photon invariant masses
  for (int i = 0; i < fTR->NJets; ++i) {
    if (fTR->JGood[i] != 0) continue;
    for (int j = i+1; j < fTR->NJets; ++j) {
      if (fTR->JGood[j] != 0) continue;
      p1.SetPxPyPzE(fTR->JPx[i],fTR->JPy[i],fTR->JPz[i],fTR->JE[i]);
      p2.SetPxPyPzE(fTR->JPx[j],fTR->JPy[j],fTR->JPz[j],fTR->JE[j]);
      psum = p1 + p2;
      minv = psum.M();
      cout << " Inv. mass (jet" << i << ", jet" << j << ") = " << minv << endl;
    }
  }
  for (int i = 0; i < fTR->NMus; ++i) {
    if (fTR->MuGood[i] != 0 || fTR->MuIsIso[i] == 0) continue;
    for (int j = i+1; j < fTR->NMus; ++j) {
      if (fTR->MuGood[j] != 0 || fTR->MuIsIso[j] == 0) continue;
      p1.SetPtEtaPhiM(fTR->MuPt[i],fTR->MuEta[i],fTR->MuPhi[i],0.10566);
      p2.SetPtEtaPhiM(fTR->MuPt[j],fTR->MuEta[j],fTR->MuPhi[j],0.10566);
      psum = p1 + p2;
      minv = psum.M();
      cout << " Inv. mass (muon" << i << ", muon" << j << ") = " << minv << endl;
    }
  }
  for (int i = 0; i < fTR->NEles; ++i) {
    if (fTR->ElGood[i] != 0 || fTR->ElIsIso[i] == 0) continue;
    for (int j = i+1; j < fTR->NEles; ++j) {
      if (fTR->ElGood[j] != 0 || fTR->ElIsIso[j] == 0) continue;
      p1.SetPtEtaPhiM(fTR->ElPt[i],fTR->ElEta[i],fTR->ElPhi[i],0.00051);
      p2.SetPtEtaPhiM(fTR->ElPt[j],fTR->ElEta[j],fTR->ElPhi[j],0.00051);
      psum = p1 + p2;
      minv = psum.M();
      cout << " Inv. mass (elec" << i << ", elec" << j << ") = " << minv << endl;
    }
  }
  for (int i = 0; i < fTR->NMus; ++i) {
    if (fTR->MuGood[i] != 0 || fTR->MuIsIso[i] == 0) continue;
    for (int j = 0; j < fTR->NEles; ++j) {
      if (fTR->ElGood[j] != 0 || fTR->ElIsIso[j] == 0) continue;
      p1.SetPtEtaPhiM(fTR->MuPt[i],fTR->MuEta[i],fTR->MuPhi[i],0.10566);
      p2.SetPtEtaPhiM(fTR->ElPt[j],fTR->ElEta[j],fTR->ElPhi[j],0.00051);
      psum = p1 + p2;
      minv = psum.M();
      cout << " Inv. mass (muon" << i << ", elec" << j << ") = " << minv << endl;
    }
  }
  for (int i = 0; i < fTR->NPhotons; ++i) {
    if (fTR->PhoGood[i] != 0 || fTR->PhoIsIso[i] == 0) continue;
    for (int j = i+1; j < fTR->NPhotons; ++j) {
      if (fTR->PhoGood[j] != 0 || fTR->PhoIsIso[j] == 0) continue;
      p1.SetPtEtaPhiM(fTR->PhoPt[i],fTR->PhoEta[i],fTR->PhoPhi[i],0.);
      p2.SetPtEtaPhiM(fTR->PhoPt[j],fTR->PhoEta[j],fTR->PhoPhi[j],0.);
      psum = p1 + p2;
      minv = psum.M();
      cout << " Inv. mass (phot" << i << ", phot" << j << ") = " << minv << endl;
    }
  }

}
