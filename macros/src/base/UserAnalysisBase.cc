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

void UserAnalysisBase::Begin(){}

void UserAnalysisBase::Analyze(){}

void UserAnalysisBase::End(){}

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
		if(fVerbose > 0) cout << "UserAnalysisBase::GetPDGParticle ==> PDGMap not fille!" << endl;
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

void UserAnalysisBase::GetHLTNames(){
	TFile *f = fTR->fChain->GetCurrentFile();
	TH1I *hlt_stats = (TH1I*)f->Get("analyze/HLTTriggerStats");

	for( int i=0; i < hlt_stats->GetNbinsX(); i++ ){
		fHLTLabelMap[hlt_stats->GetXaxis()->GetBinLabel(i+1)] = i;
	}
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
	if(fTR->MuPtE[index]/fTR->MuPt[index] > 0.5) return false;
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

// Event selections:
bool UserAnalysisBase::IsGoodMuEvent(){
	// if(GetHLTResult("")) return false;
	return true;
}

bool UserAnalysisBase::IsGoodElEvent(){
	// if(GetHLTResult("")) return false;
	return true;
}

bool UserAnalysisBase::SingleMuonSelection(int &index){
	// Selects events with (at least) one good muon and gives the index
	// of the hardest one in the argument
	// Assumes the muons are sorted by pt in the ntuple
	if( fTR->NMus < 1 ) return false;
	vector<int> MuInd;
	for(size_t imu = 0; imu < fTR->NMus; ++imu){
		if(fTR->MuPt[imu] < 10.) continue;
		if(!IsGoodBasicMu(imu)) continue;
		MuInd.push_back(imu);
	}
	if(MuInd.size() < 1) return false;
	// Maybe add a sorting here, in case muons are not sorted prior
	index = MuInd[0];
	return true;
}

bool UserAnalysisBase::DiMuonSelection(int &ind1, int &ind2){
	// Selects events with (at least) two good muons and gives the indices
	// of the hardest two in the argument
	// Assumes the muons are sorted by pt in the ntuple
	if( fTR->NMus < 2 ) return false;
	vector<int> MuInd;
	for(size_t imu = 0; imu < fTR->NMus; ++imu){
		if(fTR->MuPt[imu] < 10.) continue;
		if(!IsGoodBasicMu(imu)) continue;
		MuInd.push_back(imu);
	}
	if(MuInd.size() < 2) return false;
	ind1 = MuInd[0];
	ind2 = MuInd[1];
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
