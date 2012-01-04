#include "base/TreeReader.hh"
#include <stdlib.h>

#include "helper/pdgparticle.hh"
#include "helper/Monitor.hh"
#include "base/UserAnalysisBase.hh"
#include "TH1I.h"
#include "TLorentzVector.h"
#include "TSystem.h"

#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"

using namespace std;

UserAnalysisBase::UserAnalysisBase(TreeReader *tr){
	fTR = tr;
	fTlat = new TLatex();
	fVerbose = false;
	fDoPileUpReweight = false;
	fDoPileUpReweight3D = false;

        // Put all JES-related stuff between pre-compiler flags
#ifdef DOJES	
	//----------- Correction Object ------------------------------
	vector<JetCorrectorParameters> JetCorPar;
	JetCorrectorParameters *ResJetPar = new JetCorrectorParameters("/shome/pnef/MT2Analysis/Code/JetEnergyCorrection/GR_R_42_V19_AK5PF/GR_R_42_V19_AK5PF_L2L3Residual.txt");
	JetCorrectorParameters *L3JetPar  = new JetCorrectorParameters("/shome/pnef/MT2Analysis/Code/JetEnergyCorrection/GR_R_42_V19_AK5PF/GR_R_42_V19_AK5PF_L3Absolute.txt");
	JetCorrectorParameters *L2JetPar  = new JetCorrectorParameters("/shome/pnef/MT2Analysis/Code/JetEnergyCorrection/GR_R_42_V19_AK5PF/GR_R_42_V19_AK5PF_L2Relative.txt");
	JetCorrectorParameters *L1JetPar  = new JetCorrectorParameters("/shome/pnef/MT2Analysis/Code/JetEnergyCorrection/GR_R_42_V19_AK5PF/GR_R_42_V19_AK5PF_L1FastJet.txt");
	JetCorPar.push_back(*L1JetPar);
	JetCorPar.push_back(*L2JetPar);
	JetCorPar.push_back(*L3JetPar);
	JetCorPar.push_back(*ResJetPar);

	fJetCorrector = new FactorizedJetCorrector(JetCorPar);
        jecUnc = new JetCorrectionUncertainty("/shome/pnef/MT2Analysis/Code/JetEnergyCorrection/GR_R_42_V19_AK5PF/GR_R_42_V19_AK5PF_Uncertainty.txt");
	delete L1JetPar, L2JetPar, L3JetPar, ResJetPar, jecUnc;
#endif

}

UserAnalysisBase::~UserAnalysisBase(){
	if(fDoPileUpReweight) delete fPUWeight;
	if(fDoPileUpReweight3D) delete fPUWeight3D;


#ifdef DOJES       
	delete fJetCorrector;
        delete jecUnc;
#endif

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

	fHLTLabelMap.clear();
	fHLTLabels.clear();
	fHLTLabels.reserve(HLTNames->size());
	for( int i=0; i < HLTNames->size(); i++ ){
		fHLTLabelMap[(*HLTNames)[i]] = i; 
		fHLTLabels.push_back((*HLTNames)[i]);
		if (fVerbose>3) cout << " " << i << " " << (*HLTNames)[i] << endl;
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
		if (bit < 0 ) {
			if(fVerbose > 1) cout << "UserAnalysisBase::GetHLTResult ==> Bit with name " << theHltName << " not found!" << endl;
			return false;
		}
		else return (bool)fTR->HLTResults[bit];
	}
}

int UserAnalysisBase::GetHLTPrescale(string theHltName){
	if( fHLTLabelMap.empty() ) return 0;
	else{
		int bit = GetHLTBit(theHltName);
		if (bit < 0 ) {
			if(fVerbose > 1) cout << "UserAnalysisBase::GetHLTPrescale ==> Bit with name " << theHltName << " not found!" << endl;
			return 0;
		}
		else return fTR->HLTPrescale[bit];
	}
}

void UserAnalysisBase::GetEvtEmChFrac(double & fracEm, double & fracCh){
// Computes the event EM and Charged fractions
       std::cerr << "NEED TO REVISE" << std::endl;
        exit(-1);

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
		double pt = ( fTR->JChargedHadFrac[i] + fTR->JChargedEmFrac[i] ) * fTR->JPt[i];
		if (pt < 0.) pt = 0.;
		double em = ( fTR->JNeutralEmFrac[i] + fTR->JNeutralEmFrac[i] ) * fTR->JEt[i];
		if (em < 0.) em = 0.;
		double had = ( fTR->JChargedHadFrac[i] + fTR->JNeutralHadFrac[i] ) * fTR->JEt[i];
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
	if(fTR->JPt[index] < ptcut           ) return false;
	if(fabs(fTR->JEta[index]) > absetacut) return false;
	// Loose PF jet ID (WARNING: HF not included in our ntuple)
	// See PhysicsTools/SelectorUtils/interface/PFJetIDSelectionFunctor.h
	if ( !(fTR->JNConstituents[index] > 1) )    return false;
	if ( !(fTR->JNeutralEmFrac[index]     < 0.99) ) return false;
	if ( !(fTR->JNeutralHadFrac[index]    < 0.99) ) return false;
	if (fabs(fTR->JEta[index]) < 2.4 ) { // Cuts for |eta|<2.4
		if ( !(fTR->JChargedEmFrac[index]  < 0.99) )  return false;
		if ( !(fTR->JChargedHadFrac[index] > 0.00) )  return false;
		if ( !(fTR->JNAssoTracks[index]    > 0   ) )  return false;
	}
	return true;
}

bool UserAnalysisBase::IsGoodBasicPFJetPAT(int index, double ptcut, double absetacut){
	// Basic PF jet cleaning and ID cuts
	// cut at pt of ptcut (default = 30 GeV)
	// cut at abs(eta) of absetacut (default = 2.5)
//	if(fTR->PF2PATJPt[index] < ptcut           ) return false;
//	if(fabs(fTR->PF2PATJEta[index]) > absetacut) return false;
//	if(fTR->PF2PATJIDLoose[index]    ==0       ) return false;
	return false;
}

bool UserAnalysisBase::IsGoodPFJetMedium(int index, double ptcut, double absetacut) {
	// Medium PF JID
	if ( ! IsGoodBasicPFJet(index, ptcut, absetacut) ) return false;
	if ( !(fTR->JNeutralHadFrac[index] < 0.95)         ) return false;
	if ( !(fTR->JNeutralEmFrac[index]  < 0.95)         ) return false;
	return true;
}

bool UserAnalysisBase::IsGoodPFJetMediumPAT(int index, double ptcut, double absetacut) {
	// Medium PF JID
	if ( ! IsGoodBasicPFJetPAT(index, ptcut, absetacut) ) return false;
//	if ( !(fTR->PF2PATJNeuHadfrac[index] < 0.95)         ) return false;
//	if ( !(fTR->PF2PATJNeuEmfrac[index]  < 0.95)         ) return false;
	return true;
}

bool UserAnalysisBase::IsGoodPFJetTight(int index, double ptcut, double absetacut) {
	// Tight PF JID
	if ( ! IsGoodBasicPFJet(index, ptcut, absetacut) ) return false;
	if ( !(fTR->JNeutralHadFrac[index] < 0.90)         ) return false;
	if ( !(fTR->JNeutralEmFrac[index]  < 0.90)         ) return false;
	return true;
}

bool UserAnalysisBase::IsGoodPFJetTightPAT(int index, double ptcut, double absetacut) {
	// Tight PF JID
	if ( ! IsGoodBasicPFJetPAT(index, ptcut, absetacut) ) return false;
//	if ( !(fTR->PF2PATJNeuHadfrac[index] < 0.90)         ) return false;
//	if ( !(fTR->PF2PATJNeuEmfrac[index]  < 0.90)         ) return false;
	return true;
}

bool UserAnalysisBase::IsGoodBasicJet(int index){
	// Basic Jet cleaning and ID cuts
	if(fTR->JPt[index] < 30) return false;
	if(fabs(fTR->JEta[index]) > 3.0) return false;

	if(fTR->CAJEt[index] - fTR->JPt[index] < -0.0001 ) return false;
	if(fTR->CAJID_n90Hits[index] < 2 ) return false;
	if(fTR->CAJID_HPD[index] > 0.98 ) return false;
	if(fTR->CAJID_RBX[index] > 0.95 ) return false;
	if(fTR->CAJEMfrac[index] > 1. ) return false;
	if(fTR->CAJEMfrac[index] < 0.01 ) return false;
	// Have a linearly decreasing cut value in the transition region where
	// the jet cones intersects with the tracker acceptance, i.e. between
	// eta 1.9 and 2.9
	const double chmin = 0.05;
	double temp = chmin;
	if(fabs(fTR->CAJEta[index]) > 1.9) temp = chmin * (1. - fabs(fTR->CAJEta[index]) + 1.9);
	if(fabs(fTR->CAJEta[index]) > 2.9) temp = 0.;
	if( fTR->CAJChfrac[index] < temp && fabs(fTR->CAJEta[index]) < 2.9) return false;
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

	if(fabs(fTR->MuD0PV[index]) > 0.02)    return false;
	if(fabs(fTR->MuDzPV[index]) > 1.00)    return false;

	if(fTR->MuIso03EMVetoEt[index] > 4.0)  return false;
	if(fTR->MuIso03HadVetoEt[index] > 6.0) return false;

	if(fTR->MuPtE[index]/fTR->MuPt[index] > 0.1) return false;

	if(fTR->MuRelIso03[index] > 1.0) return false;
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
	// Electrons with WP95 ID
	// if(fTR->ElIDsimpleWP95relIso[index] != 5 && fTR->ElIDsimpleWP95relIso[index] != 7) return false;		
	if( fabs(fTR->ElEta[index]) < 1.479 ){ // Barrel
		if(fTR->ElSigmaIetaIeta            [index] > 0.01 ) return false;
		if(fabs(fTR->ElDeltaPhiSuperClusterAtVtx[index]) > 0.80 ) return false;
		if(fabs(fTR->ElDeltaEtaSuperClusterAtVtx[index]) > 0.007) return false;
		if(fTR->ElHcalOverEcal             [index] > 0.15 ) return false;	
	}
	if( fabs(fTR->ElEta[index]) > 1.479 ){ // Endcap
		if(fTR->ElSigmaIetaIeta            [index] > 0.03 ) return false;
		if(fabs(fTR->ElDeltaPhiSuperClusterAtVtx[index]) > 0.70 ) return false;
		if(fabs(fTR->ElDeltaEtaSuperClusterAtVtx[index]) > 0.01) return false;
		// if(fTR->ElHcalOverEcal             [index] > 0.07 ) return false;	
	}
	
	// ECAL gap veto
	if ( fabs(fTR->ElSCEta[index]) > 1.4442 && fabs(fTR->ElSCEta[index]) < 1.566 )  return false;

	if(fabs(fTR->ElD0PV[index]) > 0.02) return false;
	if(fabs(fTR->ElDzPV[index]) > 1.00) return false;

	return true;
}

bool UserAnalysisBase::ElPassesWP80_ConvRej(int index){
	// if(fTR->ElNumberOfMissingInnerHits[index] > 0   ) return false;
	// if(fabs(fTR->ElConvPartnerTrkDist[index]) < 0.02 && fabs(fTR->ElConvPartnerTrkDCot[index]) < 0.02) return false;
	if(fTR->ElIDsimpleWP80relIso[index] < 4) return false;
	return true;
}

bool UserAnalysisBase::IsGoodElId_WP90(int index){
	// Electrons with WP90 ID and WP80 conv. rej.
	if( fabs(fTR->ElEta[index]) < 1.479 ){ // Barrel
		if(fTR->ElSigmaIetaIeta            [index] > 0.01 ) return false;
		if(fabs(fTR->ElDeltaPhiSuperClusterAtVtx[index]) > 0.80 ) return false;
		if(fabs(fTR->ElDeltaEtaSuperClusterAtVtx[index]) > 0.007) return false;
		if(fTR->ElHcalOverEcal             [index] > 0.12 ) return false;	
	}
	if( fabs(fTR->ElEta[index]) > 1.479 ){ // Endcap
		if(fTR->ElSigmaIetaIeta            [index] > 0.03 ) return false;
		if(fabs(fTR->ElDeltaPhiSuperClusterAtVtx[index]) > 0.70 ) return false;
		if(fabs(fTR->ElDeltaEtaSuperClusterAtVtx[index]) > 0.009) return false;
		// if(fTR->ElHcalOverEcal             [index] > 0.15 ) return false;	
	}

	// WP80 conv. rejection
	if(ElPassesWP80_ConvRej(index) == false) return false;
	return true;
}

bool UserAnalysisBase::IsGoodElId_WP80(int index){
	// Electrons with WP80 ID and conv. rej. cuts
	if( fabs(fTR->ElEta[index]) < 1.479 ){ // Barrel
		if(fTR->ElSigmaIetaIeta            [index] > 0.01 ) return false;
		if(fabs(fTR->ElDeltaPhiSuperClusterAtVtx[index]) > 0.06 ) return false;
		if(fabs(fTR->ElDeltaEtaSuperClusterAtVtx[index]) > 0.004) return false;
		if(fTR->ElHcalOverEcal             [index] > 0.04 ) return false;	
	}
	if( fabs(fTR->ElEta[index]) > 1.479 ){ // Endcap
		if(fTR->ElSigmaIetaIeta            [index] > 0.03 ) return false;
		if(fabs(fTR->ElDeltaPhiSuperClusterAtVtx[index]) > 0.03 ) return false;
		if(fabs(fTR->ElDeltaEtaSuperClusterAtVtx[index]) > 0.007) return false;
		// if(fTR->ElHcalOverEcal             [index] > 0.15 ) return false;	
	}
	if(fTR->ElPt[index] < 20.){
		if(fTR->Elfbrem[index] > 0.15) return true;
		if(fabs(fTR->ElSCEta[index]) < 1.0 && fTR->ElESuperClusterOverP[index] > 0.95 ) return true;
		return false;
	}

	// WP80 conv. rejection
	if(ElPassesWP80_ConvRej(index) == false) return false;
	return true;
}

float UserAnalysisBase::relElIso(int index){
	// Apply 1 GeV subtraction in ECAL Barrel isolation
	if( fabs(fTR->ElEta[index]) < 1.479 ){ // Barrel
		double iso = ( fTR->ElDR03TkSumPt[index] + TMath::Max(0., fTR->ElDR03EcalRecHitSumEt[index] - 1.) + fTR->ElDR03HcalTowerSumEt[index] ) / fTR->ElPt[index];
		return iso;
	}
	return fTR->ElRelIso03[index]; // EndCap
}

bool UserAnalysisBase::IsLooseEl(int index){
	if(!IsGoodBasicEl(index))         return false;
	if(fTR->ElPt[index] < 10.)        return false;
	if(fabs(fTR->ElEta[index]) > 2.4) return false;

	if(!fTR->ElEcalDriven[index]) return false;
	
	// Loose identification criteria: WP90
	if(!IsGoodElId_WP90(index)) return false;

	// Loose isolation criteria
	if( fabs(fTR->ElEta[index]) <= 1.479 && relElIso(index) > 1.0 ) return false;
	if( fabs(fTR->ElEta[index]) >  1.479 && relElIso(index) > 0.6 ) return false;

	return true;
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
	if (fTR->PrimVtxNdof < 4.0) return false; // this is a float, cut value at 4.0
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
	for(int ind = 0; ind < fTR->NJets; ++ind){
		// selection
		if((*this.*pfjetSelector)(ind, ptcut, absetacut) == false) continue;

		selectedObjInd.push_back(ind);
		selectedObjPt.push_back(fTR->JPt[ind]);
	}
	return Util::VSort(selectedObjInd, selectedObjPt);
}

vector<int> UserAnalysisBase::PhotonSelection(bool(UserAnalysisBase::*phoSelector)(int)){
        std::cerr << "NEED TO REVISE" << std::endl;
        exit(-1);
	// Returns the vector of indices of
	// good photons sorted by Pt
	if(phoSelector == NULL) phoSelector = &UserAnalysisBase::IsGoodBasicPho;
	vector<int>    selectedObjInd;
	vector<double> selectedObjPt;
	// form the vector of indices
	for(int ind = 0; ind < fTR->NPhotons; ++ind){
		// selection
		if((*this.*phoSelector)(ind) == false) continue;
		//if(fTR->PhoIsElDupl[ind] >= 0) continue;
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
// JEC
#ifdef DOJES
float UserAnalysisBase::GetJetPtNoResidual(int jetindex){
	if(!fIsData) return fTR->JPt[jetindex]; // do nothing for MC (no res corr here)

	float rawpt = fTR->JPt[jetindex]/fTR->JEcorr[jetindex];
	fJetCorrector->setJetEta(fTR->JEta[jetindex]);
	fJetCorrector->setRho(fTR->Rho);
	fJetCorrector->setJetA(fTR->JArea[jetindex]);
	fJetCorrector->setJetPt(rawpt); // IMPORTANT: the correction is a function of the RAW pt

	// The getSubCorrections member function returns the vector of the subcorrections UP to
	// the given level. For example in the example above, factors[0] is the L1 correction
	// and factors[3] is the L1+L2+L3+Residual correction.
	vector<float> factors = fJetCorrector->getSubCorrections();

	// Sanity check: JEcorr should be the full set of corrections applied
	if(fabs(factors[3] - fTR->JEcorr[jetindex]) > 0.000001 && fVerbose > 2) cout << "UserAnalysisBase::GetJetPtNoResidual ==> WARNING: Your JECs don't seem to be consistent!" << endl;

	double l1l2l3scale = factors[2];
	return rawpt*l1l2l3scale;
}
#endif

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
	if (fTR->PF2PAT3NJets == 1 && fabs(fTR->PF2PAT3JEta[0]) < etaBar) tag = "***";
	if (fTR->PF2PAT3NJets >= 2 && (fabs(fTR->PF2PAT3JEta[0]) < etaBar || fabs(fTR->PF2PAT3JEta[1] < etaBar))) tag = "***";
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
//  GetEvtEmChFrac(fracEm, fracCh);
//  cout << " Event fEM = " << fracEm << ", fCh = " << fracCh << endl;
	cout << " TCMET          = " << fTR->TCMET    << ", METPhi = " << fTR->TCMETphi << endl;
	cout << " PFMET          = " << fTR->PFMET    << ", METPhi = " << fTR->PFMETphi << endl;

// print the jets info
	cout << " Number of jets in the ntuple = " << fTR->PF2PAT3NJets << endl;
	int nJetsCand = 0;
	for (int i = 0; i < fTR->PF2PAT3NJets; ++i) {
		double dPhiMJ = Util::DeltaPhi(fTR->PF2PAT3JPhi[i], METPhi);
		cout << " Jet" << i << " Pt = " << fTR->PF2PAT3JPt[i] << ", Phi = " << fTR->PF2PAT3JPhi[i]
			<< ", Eta = " << fTR->PF2PAT3JEta[i] << ", E = " << fTR->PF2PAT3JE[i] << endl;
		double jmass = sqrt(fTR->PF2PAT3JE[i]*fTR->PF2PAT3JE[i]-fTR->PF2PAT3JPx[i]*fTR->PF2PAT3JPx[i]-fTR->PF2PAT3JPy[i]*fTR->PF2PAT3JPy[i]-fTR->PF2PAT3JPz[i]*fTR->PF2PAT3JPz[i]);
		cout << "      " << " Px = " << fTR->PF2PAT3JPx[i] << " Py = " << fTR->PF2PAT3JPy[i] << " Pz = " << fTR->PF2PAT3JPz[i]
			<< " Jet mass = " << jmass << endl;
		cout << "      " << " dPhiJM = " << dPhiMJ <<  endl;

		cout << "      " << " pf-JID    (1=true, 0=false)= " << fTR->PF2PAT3JIDLoose[i] << endl;  
		cout << "      " <<" JetfEMch  = " << fTR->PF2PAT3JChEmfrac[i]  << ", JetfEMneu = "  << fTR->PF2PAT3JNeuEmfrac[i] << endl;
		cout << "      " <<" JetfHach  = " << fTR->PF2PAT3JChHadfrac[i] << ", JetfHaneu = "  << fTR->PF2PAT3JNeuHadfrac[i] <<  endl;
		cout << "      " <<" JetCHMult = " << fTR->PF2PAT3JChMult[i]    << ", JetNeuMult = " << fTR->PF2PAT3JNeuMult[i]    << endl; 
		cout << "      " << " b-tag TCHE = " << fTR->PF2PAT3JbTagProbTkCntHighEff[i]<< ", TCHP = " << fTR->PF2PAT3JbTagProbTkCntHighPur[i] << endl;
		cout << "      " << " b-tag SSVHE = " << fTR->PF2PAT3JbTagProbSimpSVHighEff[i]<< ", SSVHP = " << fTR->PF2PAT3JbTagProbSimpSVHighPur[i] << endl;
		nJetsCand++;
		if (fTR->PF2PAT3JPt[i] > 20. && fabs(dPhiMJ-3.141592654) < 0.05) {
			cout << "      -> jet back-to-back with MET" << endl;
		}
	}

// for multi-jet events
	if (fTR->PF2PAT3NJets >= 2) {
		double dPhiMJ1 = Util::DeltaPhi(fTR->PF2PAT3JPhi[0], METPhi);
		double dPhiMJ2 = Util::DeltaPhi(fTR->PF2PAT3JPhi[1], METPhi);
		double R12 = sqrt(dPhiMJ1*dPhiMJ1 + (TMath::Pi()-dPhiMJ2)*(TMath::Pi()-dPhiMJ2) );
		double R21 = sqrt(dPhiMJ2*dPhiMJ2 + (TMath::Pi()-dPhiMJ1)*(TMath::Pi()-dPhiMJ1) );
		double dPhij12 = Util::DeltaPhi(fTR->PF2PAT3JPhi[0], fTR->PF2PAT3JPhi[1]);
		double dRj12 = Util::GetDeltaR(fTR->PF2PAT3JEta[0], fTR->PF2PAT3JEta[1], fTR->PF2PAT3JPhi[0], fTR->PF2PAT3JPhi[1]);
		cout << " R12    = " << R12 << ", R21    = " << R21 
			<< ", dRj12  = " << dRj12 << ", dPhij12 = " << dPhij12 << endl;
	}

	if(fTR->PfTau3NObjs>0) cout << "PF2PAT Taus: these hadronic taus are contained in above jet collection! " << endl;
	for(int i=0; i<fTR->PfTau3NObjs; ++i){
		TLorentzVector tau(fTR->PfTau3Px[i],fTR->PfTau3Py[i],fTR->PfTau3Pz[i],fTR->PfTau3E[i]);
		cout << " Tau " << i <<  " has mass " << tau.M() << " pt " << tau.Pt() << " Eta " << tau.Eta() << " Phi " << tau.Phi()<< endl;
		cout << "       px " << tau.Px() << ", py " << tau.Py() << ", pz " << tau.Pz() << ", E " << tau.E() << endl; 
		cout << "       reconstruced tau decay mode (5*(number of charged hadrons-1) + number of pi0s.) " << fTR->PfTau3DecayMode[i] << endl;
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
		cout << "      " << " chisq = " << fTR->MuNChi2[i] << " #Trk hits = " << fTR->MuNTkHits[i] << " NMatches " << fTR->MuNMatches[i] << endl;
		cout << "      " << " NPixel hits " << fTR->MuNPxHits[i] << " global track Nhits " << fTR->MuNGlHits[i] << endl;
		cout << "      " << " IsoVal03 = " << fTR->MuRelIso03[i] << endl;
		cout << "      " << " IsoVal Trk = " << fTR->MuIso03SumPt[i] << " Em = " << fTR->MuIso03EmEt[i] << " Had = " << fTR->MuIso03HadEt[i] << endl;
		cout << "      " << " Energy assoc to the muon Em = " << fTR->MuEem[i] << " Had = " << fTR->MuEhad[i] << endl;
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
		cout << "      " << " IsoVal03 = " << fTR->ElRelIso03[i] << endl;
		cout << "      " << " IsoVal Trk = " << fTR->ElDR03TkSumPt[i] << " Em = " << fTR->ElDR03EcalRecHitSumEt[i] 
			<< " Had = " << fTR->ElDR03HcalTowerSumEt[i] << endl;
		if (fTR->ElRelIso03[i] <= elIsomax) cout << "      -> electron is isolated " << endl;
		else cout << "      -> electron NOT isolated " << endl;
	}

// print photon info
	for (int i = 0; i < fTR->NPhotons; ++i) {
		cout << " Phot" << i << " Pt = " << fTR->PhoPt[i] << ", Phi = " << fTR->PhoPhi[i]
			<< ", Eta = " << fTR->PhoEta[i] << ", E = " << fTR->PhoEnergy[i] << endl;
		cout << "      " << " Px = " << fTR->PhoPx[i] << " Py = " << fTR->PhoPy[i] << " Pz = " << fTR->PhoPz[i] << endl;
		cout << "      " << " HoverE = " << fTR->PhoHoverE[i] << " sigmaEtaEta = " << fTR->PhoSigmaIetaIeta[i] << endl;
		cout << "      " << " IsoVal Trk = " << fTR->PhoIso03TrkSolid[i] << " Em = " << fTR->PhoIso03Ecal[i] 
			<< " Had = " << fTR->PhoIso03Hcal[i] << endl;
		cout << "      " << " Has conversion tracks = " << fTR->PhoHasConvTrks[i] << endl;
		if (fTR->PhoIso03[i] <= phoIsomax) cout << "      -> photon is isolated " << endl;
		else cout << "      -> photon NOT isolated " << endl;
	}

// jets/ lepton/ photon invariant masses
	for (int i = 0; i < fTR->PF2PAT3NJets; ++i) {
		for (int j = i+1; j < fTR->PF2PAT3NJets; ++j) {
			p1.SetPxPyPzE(fTR->PF2PAT3JPx[i],fTR->PF2PAT3JPy[i],fTR->PF2PAT3JPz[i],fTR->PF2PAT3JE[i]);
			p2.SetPxPyPzE(fTR->PF2PAT3JPx[j],fTR->PF2PAT3JPy[j],fTR->PF2PAT3JPz[j],fTR->PF2PAT3JE[j]);
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


// ---------------------------------------------
// Pile Up Reweighting
void UserAnalysisBase::SetPileUpSrc(string data_PileUp, string mc_PileUp){
	if(data_PileUp.size() == 0 ) return;
	if(fDoPileUpReweight  == 1 ) {cout << "ERROR in SetPileUpSrc: fPUWeight already initialized" << endl; return; }
	if(mc_PileUp.size() ==0){
		fPUWeight = new PUWeight(data_PileUp.c_str());
	}else {
		fPUWeight = new PUWeight(data_PileUp.c_str(), mc_PileUp.c_str());
	} 
	fDoPileUpReweight = true;
}

// Pile Up Reweighting - 3D
void UserAnalysisBase::SetPileUp3DSrc(string data_PileUp, string mc_PileUp){
	if(data_PileUp.size() == 0 ) return;
	if(fDoPileUpReweight3D == 1 ) {cout << "ERROR in SetPileUpSrc3D: fPUWeight already initialized" << endl; return; }
	//if(mc_PileUp.size() ==0){
	//	fPUWeight = new Lumi3DReWeighting(data_PileUp.c_str());
	//}else {
	fPUWeight3D = new Lumi3DReWeighting( mc_PileUp, data_PileUp,"pileup","pileup");
	fPUWeight3D->weight3D_init(1);
		//} 
	fDoPileUpReweight3D = true;
}

float UserAnalysisBase::GetPUWeight(int nPUinteractions){
	if(! fDoPileUpReweight) return -999.99;
	else return fPUWeight->GetWeight(nPUinteractions); 
}

float UserAnalysisBase::GetPUWeight(int nPUinteractions, int nPUinteractionsLate){
	if(! fDoPileUpReweight) return -999.99;
	else return fPUWeight->weightOOT(nPUinteractions, nPUinteractionsLate);
}

float UserAnalysisBase::GetPUWeight3D( int nPUinteractionsEarly, int nPUinteractions, int nPUinteractionsLate){
  //  return -1;
  if(! fDoPileUpReweight3D) return -999.99;
  else return fPUWeight3D->weight3D(nPUinteractionsEarly ,nPUinteractions , nPUinteractionsLate);
}

