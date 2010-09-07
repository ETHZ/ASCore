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
	TFile *f = fTR->fChain->GetCurrentFile();
        TTree* runTree = (TTree*)f->Get("analyze/RunInfo");
        std::vector<std::string>* HLTNames;
        if ( !runTree ) {
          std::cerr << "!!! UserAnalysisBase::GetHLTNames "
                    << "Coudln't get analyze/RunInfo tree" << std::endl;
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

	if( fTR->JPt[index] < 40. ) 			return false;
	if( fabs(fTR->JEta[index]) > 2.5 ) 		return false;
	if (fTR->JEMfrac[index] <= 0.01)  		return false;
	if (fTR->JID_n90Hits[index] <= 1) 		return false;
	if (fTR->JID_HPD[index] >= 0.98)  		return false;

	return true;
}

bool UserAnalysisBase::IsGoodbJ_TDL(int index) {
	if(! IsGoodJ_TDL(index)) return false;
	if(fTR->JbTagProbSimpSVHighEff[index] < 1.74) return false;

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
	// Basic muon cleaning and ID cuts
	if(fTR->MuIsGlobalMuon[index] == 0) return false;
	if(fTR->MuIsTrackerMuon[index] == 0) return false;
	if(fTR->MuPt[index] < 10) return false;
	if(fTR->MuNChi2[index] > 10) return false;
	if(fTR->MuNTkHits[index] < 11) return false;
	if(fabs(fTR->MuD0BS[index]) > 0.02) return false;
	if(fTR->MuIsGMPT[index] == 0) return false;
	if(fTR->MuRelIso03[index] > 1.0) return false;
	// if(fTR->MuPtE[index]/fTR->MuPt[index] > 0.5) return false; // No effect after other cuts (on ~ 120/nb)
	return true;
}

bool UserAnalysisBase::IsGoodMu_TDL(int index){
	if(!IsGoodBasicMu(index)) return false;
	if(fTR->MuPt[index] < 20) return false;
	if( fabs(fTR->MuEta[index]) > 2.5 ) return false;
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
	
	if( fTR->ElPt[index] < 20. ) 	return false;
	if( fabs(fTR->ElEta[index]) > 2.5 ) return false;
	if( fabs(fTR->ElEta[index]) > etaEgapLowerEdge && fabs(fTR->ElEta[index]) < etaEgapUpperEdge) return false;
	if( fTR->ElD0BS[index] >= 0.04 ) return false;
	
	// conversion rejection
	if( fTR->ElNumberOfMissingInnerHits[index] > 1) return false;
	if(fabs(fTR->ElConvPartnerTrkDist[index]) < 0.02 && fabs(fTR->ElConvPartnerTrkDCot[index]) < 0.02) return false;
	
	// el id
	if(fabs(fTR->ElEta[index]) < etaEgapLowerEdge){// for barrel
		if(fTR->ElSigmaIetaIeta[index] > 0.01) return false;
		if(fabs(fTR->ElDeltaPhiSuperClusterAtVtx[index]) > 0.8) return false;
		if(fabs(fTR->ElDeltaEtaSuperClusterAtVtx[index]) > 0.007) return false;
		if( fTR->ElHcalOverEcal[index] > 0.12 ) return false;
	}else if(fabs(fTR->ElEta[index]) > etaEgapUpperEdge){ // for endcap
		if(fTR->ElSigmaIetaIeta[index] > 0.03) return false;
		if(fabs(fTR->ElDeltaPhiSuperClusterAtVtx[index]) > 0.7) return false;
		if(fabs(fTR->ElDeltaEtaSuperClusterAtVtx[index]) > 0.009) return false;
		if( fTR->ElHcalOverEcal[index] > 0.05 ) return false;		
	}
	
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
	// conversions cut values
	double ElecConvPartTrackDistmin = 0.02;
	double ElecConvPartTrackDCotmin = 0.02;
	double ElecNMissHitsmax         = 1.;
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

	// electrons from conversions
	if ( fTR->ElNumberOfMissingInnerHits[index] > ElecNMissHitsmax ) return false;
	if ( fTR->ElConvPartnerTrkDist[index] < ElecConvPartTrackDistmin &&
	     fTR->ElConvPartnerTrkDCot[index] < ElecConvPartTrackDCotmin ) return false;
	
	// electron kinematic cuts + the ECAL gap veto
	if ( (fabs(eta) > etamax) ||
	   ( (fabs(eta) > etaEgapLowerEdge) && (fabs(eta) < etaEgapUpperEdge) ) ) return false;	
	if( pt < pTmin) return false;

	return true;
}

bool UserAnalysisBase::IsTightEl(int index){
	if(!IsGoodBasicEl(index)) return false;	
	if(!fTR->ElIDsimpleWPrelIso[index]) return false; // corresponding to the 90% eff. WP with isolation cuts
	
	// combined relative isolation cut values (for tight electrons)
	double ElecCombRelIsoEBarmax = 0.15;
	double ElecCombRelIsoEEndmax = 0.1;
	// electron isolation
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

bool UserAnalysisBase::IsLooseEl(int index){
	if(!IsGoodBasicEl(index)) return false;

	// combined relative isolation cut values (relaxed for loose electrons)
	double ElecCombRelIsoEBarmax	= 1.;
	double ElecCombRelIsoEEndmax	= 0.6;	
	// electron isolation
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

bool UserAnalysisBase::IsLooseNoTightEl(int index){
	if(IsLooseEl(index) && !IsTightEl(index)) return true;
	else return false;
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
bool UserAnalysisBase::IsGoodMuEvent(){
	string MuonHLTName = "HLT_Mu9"; // needs also matching HLT object pt>15 GeV/c; this includes MC samples where run=1
	return GetHLTResult(MuonHLTName);	
}

bool UserAnalysisBase::IsGoodElEvent(){
	string ElectronHLTName;
	int run = fTR->Run;	
	if (run<=138000)				ElectronHLTName = "HLT_Ele10_LW_L1R"; // needs also matching HLT object pt>15 GeV/c; this includes MC samples where run=1
	if (run>138000 && run<=141900)	ElectronHLTName = "HLT_Ele15_LW_L1R";
	if (run>141900)					ElectronHLTName = "HLT_Ele15_SW_L1R";
	return GetHLTResult(ElectronHLTName);	
}

vector<int> UserAnalysisBase::ElectronSelection(){
	// Returns the vector of indices of
	// good electrons sorted by Pt
	vector<int>    selectedObjInd;
	vector<double> selectedObjPt;
	// form the vector of indices
	for(int ind = 0; ind < fTR->NEles; ++ind){
		// selection
		if(!IsGoodBasicEl(ind)) continue;
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

vector<int> UserAnalysisBase::PhotonSelection(){
	// Returns the vector of indices of
	// good photons sorted by Pt
	vector<int>    selectedObjInd;
	vector<double> selectedObjPt;
	// form the vector of indices
	for(int ind = 0; ind < fTR->NPhotons; ++ind){
		// selection
		if(!IsGoodBasicPho(ind)) continue;
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
		if(fTR->ElPt[ElInd[iel]] > 20.)	priElInd.push_back(ElInd[iel]);
		if(fTR->ElPt[ElInd[iel]] > 15.)	secElInd.push_back(ElInd[iel]);
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

bool UserAnalysisBase::DiMuonSelection(int &ind1, int &ind2){
	// Selects events with (at least) two good muons and gives the indices
	// of the hardest two in the argument
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
	return true;
}

bool UserAnalysisBase::SSDiMuonSelection(int &prim, int &sec){
	// Select events with either 2 muons or 2 electrons
	if( fTR->NMus < 2 ) return false;

	vector<int> priMuInd;
	vector<int> secMuInd;
	for(size_t imu = 0; imu < fTR->NMus; ++imu){
		// Muon selection
		if(fTR->MuPt[imu] < 10.) continue;
		if(!IsGoodBasicMu(imu)) continue;

		if(fTR->MuPt[imu] > 20. && fTR->MuRelIso03[imu] < 0.1) priMuInd.push_back(imu);
		if(fTR->MuRelIso03[imu] < 2.0)                         secMuInd.push_back(imu);
	}

	unsigned nprimu = priMuInd.size();
	unsigned nsecmu = secMuInd.size();

	// Select events with at least one primary muon and at least one other
	if(nprimu < 1 || nsecmu < 2) return false;

	int ind1 = priMuInd[0];
	int ind2 = secMuInd[0];
	if(ind2 == ind1) ind2 = secMuInd[1];
	// Select same sign muons
	if(fTR->MuCharge[ind1] != fTR->MuCharge[ind2]) return false;
	prim = ind1;
	sec = ind2;
	return true;
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
