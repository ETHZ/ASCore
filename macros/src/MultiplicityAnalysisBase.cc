#include "helper/Utilities.hh"
#include "MultiplicityAnalysisBase.hh"


using namespace std;

MultiplicityAnalysisBase::MultiplicityAnalysisBase(TreeReader *tr) : UserAnalysisBase(tr){
	Util::SetStyle();	
	fCut_PFMET_min                      = 0;
	fCut_HT_min                         = 0;
	fCut_JPt_hardest_min                = 0;
	fCut_DiLeptInvMass_min              = 0;
	fCut_DiLeptInvMass_max              = 9999.99;
	fCut_DiLeptOSSFInvMass_lowercut     = 9999.99;
	fCut_DiLeptOSSFInvMass_uppercut     =-9999.99;
	fCut_PtHat_max                      = 999999.;
	fCut_Run_min                        = 0;
	fCut_Run_max                        = 9999999;
	
	fRequiredHLT.clear();
	fVetoedHLT.clear();
}

MultiplicityAnalysisBase::~MultiplicityAnalysisBase(){
}



void MultiplicityAnalysisBase::InitializeEvent(){
	GetLeptonJetIndices();
	FindLeptonConfig();
}

void MultiplicityAnalysisBase::GetLeptonJetIndices(){
	fElecs.clear();
	fMuons.clear();
	fJets.clear();
	fBJets.clear();

	for(int i=0; i< fTR->NMus; ++i){
		if(! IsGoodMu_TDL(i) ) continue;
		fMuons.push_back(i);
	}
	
	for(int i=0; i< fTR->NEles; ++i){
		if(! IsGoodEl_TDL(i) ) continue;
		fElecs.push_back(i);
	}
	
	for(int ij=0; ij < fTR->NJets; ++ij){
		if(! IsGoodJ_TDL(ij) ) continue;
		bool JGood(true);
		for(int i=0; i<fMuons.size(); ++i){
			double deltaR = Util::GetDeltaR(fTR->JEta[ij], fTR->MuEta[fMuons[i]], fTR->JPhi[ij], fTR->MuPhi[fMuons[i]]);
			if(deltaR < 0.4)   JGood=false;
		}
		for(int i=0; i<fElecs.size(); ++i){
			double deltaR = Util::GetDeltaR(fTR->JEta[ij], fTR->ElEta[fElecs[i]], fTR->JPhi[ij], fTR->ElPhi[fElecs[i]]);
			if(deltaR < 0.4)   JGood=false;
		}
		if(JGood=false) continue;
		fJets.push_back(ij);
		// bjets
		if(! IsGoodbJ_TDL(ij) ) continue;
		fBJets.push_back(ij);
	}	
}

void MultiplicityAnalysisBase::FindLeptonConfig(){
	fLeptConfig = null;
	if(fElecs.size()==1 && fMuons.size()==0){
		fLeptConfig=e;
	}else if(fElecs.size()==0 && fMuons.size()==1) {
		fLeptConfig=mu;
	}else if(fElecs.size() + fMuons.size() ==2){
		int charge1, charge2;
		if(fElecs.size()==2){
			charge1=fTR->ElCharge[fElecs[0]];
			charge2=fTR->ElCharge[fElecs[1]];
			if(charge1*charge2==1){fLeptConfig=SS_ee;}
			else{fLeptConfig=OS_ee;}
		} else if(fMuons.size()==2){
			charge1=fTR->MuCharge[fMuons[0]];
			charge2=fTR->MuCharge[fMuons[1]];
			if(charge1*charge2==1){fLeptConfig=SS_mumu;}
			else{fLeptConfig=OS_mumu;}
		} else{
			charge1=fTR->ElCharge[fElecs[0]];
			charge2=fTR->MuCharge[fMuons[0]];			
			if(charge1*charge2==1){fLeptConfig=SS_emu;}
			else{fLeptConfig=OS_emu;}
		}
	}
}


bool MultiplicityAnalysisBase::IsGoodEvent(){
	
	// Run
	if(fTR->Run < fCut_Run_min ) {return false;}
	if(fTR->Run > fCut_Run_max ) {return false;}
	
	//PtHat
	if(fTR->PtHat > fCut_PtHat_max ){return false;}
	
	// MET
	if(fTR->PFMET < fCut_PFMET_min){return false;}
	
	// HT
	double HT=0;
	for(int j=0; j<fJets.size(); ++j){
		HT += fTR->JPt[fJets[j]];
	}
	if(HT<fCut_HT_min){return false;}
	
	// hardest jet
	double hardest_jet=0;
	for(int j=0; j<fJets.size(); ++j){
		if(fTR->JPt[fJets[j]] > hardest_jet){
			hardest_jet = fTR->JPt[fJets[j]];
		}
	}
	if(hardest_jet<fCut_JPt_hardest_min){return false;}
	
	// DiLeptonInvMass_min DiLeptonInvMass_max
	if(fLeptConfig==SS_ee || fLeptConfig==OS_ee || fLeptConfig==SS_mumu || fLeptConfig==OS_mumu ||
	   fLeptConfig==SS_emu || fLeptConfig==OS_emu ){
		double invmass = GetDiLeptInvMass();
		
		if(invmass < fCut_DiLeptInvMass_min) {return false;}
		if(invmass > fCut_DiLeptInvMass_max) {return false;}		 
				 
	}
	
	// Zmass cut
	// DiLeptonInvMass_min
	if( fLeptConfig==OS_ee ||fLeptConfig==OS_mumu ){
		double invmass = GetDiLeptInvMass();
		
		if(invmass > fCut_DiLeptOSSFInvMass_lowercut) {return false;}	
		if(invmass < fCut_DiLeptOSSFInvMass_uppercut) {return false;}	
	}
	
	// HLT triggers
	if(fRequiredHLT.size() !=0 ){
		bool HLT_fire(false);
		for(int i=0; i<fRequiredHLT.size(); ++i){
			if( GetHLTResult(fRequiredHLT[i]) ){  // require HLT bit
				HLT_fire = true;
			} 
		}
		if(! HLT_fire) return false;
	}
	if(fVetoedHLT.size() !=0){
		for(int i=0; i<fVetoedHLT.size(); ++i){
			if( GetHLTResult(fVetoedHLT[i]) ){   // veto HLT bit 
				return false;
			} 
		}
	}
	
	return true;	
}


double MultiplicityAnalysisBase::GetDiLeptInvMass(){
	TLorentzVector p1, p2;
	if(fLeptConfig==SS_ee || fLeptConfig==OS_ee){
		p1.SetPtEtaPhiM(fTR->ElPt[fElecs[0]],fTR->ElEta[fElecs[0]],fTR->ElPhi[fElecs[0]], 0. );
		p2.SetPtEtaPhiM(fTR->ElPt[fElecs[1]],fTR->ElEta[fElecs[1]],fTR->ElPhi[fElecs[1]], 0. );
	} else if(fLeptConfig==SS_mumu || fLeptConfig==OS_mumu ){
		p1.SetPtEtaPhiM(fTR->MuPt[fMuons[0]],fTR->MuEta[fMuons[0]],fTR->MuPhi[fMuons[0]], 0.105 );
		p2.SetPtEtaPhiM(fTR->MuPt[fMuons[1]],fTR->MuEta[fMuons[1]],fTR->MuPhi[fMuons[1]], 0.105 );
	} else if(fLeptConfig==SS_emu || fLeptConfig==OS_emu ){
		p1.SetPtEtaPhiM(fTR->MuPt[fMuons[0]],fTR->MuEta[fMuons[0]],fTR->MuPhi[fMuons[0]], 0.105 );
		p2.SetPtEtaPhiM(fTR->ElPt[fElecs[0]],fTR->ElEta[fElecs[0]],fTR->ElPhi[fElecs[0]], 0. );
	} else {
		cout << "ERROR in MultiplicityAnalysisBase::GetDiLeptInvMass" << endl; 
		return -999.99;
	}
	double invmass = (p1+p2).M();
	return invmass;
}

void MultiplicityAnalysisBase::ReadCuts(const char* SetofCuts="multiplicity_cuts/default.dat"){
	
	ifstream IN(SetofCuts);
	char buffer[200];
	char ParName[100];
	char StringValue[100];
	float ParValue;
	int   FlagValue;
	int   IntValue;

	bool verbose(true);
	bool ok(false);
	
	
	while( IN.getline(buffer, 200, '\n') ){
		ok = false;
		if (buffer[0] == '#') {continue;} // Skip lines commented with '#'

		// strings
		sscanf(buffer, "%s %s", ParName, StringValue);
		if( !strcmp(ParName, "SetName") ){
			fSetName = TString(StringValue); ok = true;
			if(verbose){cout << "Reading cut parameters for set: " << fSetName << endl; }
		} else if( !strcmp(ParName, "HLT_required") ){
			fRequiredHLT.push_back(StringValue); ok = true;
		} else if( !strcmp(ParName, "HLT_vetoed") ){
			fVetoedHLT.push_back(StringValue); ok = true;
		}	

		// ints
		sscanf(buffer, "%s %i", ParName, &IntValue);
		if( !strcmp(ParName, "Run_min") ){
			fCut_Run_min = int(IntValue); ok = true;
		} else if( !strcmp(ParName, "Run_max") ){
			fCut_Run_max = int(IntValue); ok = true;
		}

		// floats 
		sscanf(buffer, "%s %f", ParName, &ParValue);
		if( !strcmp(ParName, "PFMET_min") ){
			fCut_PFMET_min            = float(ParValue); ok = true;
		} else if( !strcmp(ParName, "HT_min") ){
			fCut_HT_min               = float(ParValue); ok = true;
		} else if( !strcmp(ParName, "JPt_hardest_min") ){
			fCut_JPt_hardest_min      = float(ParValue); ok = true;			
		} else if( !strcmp(ParName, "DiLeptInvMass_min") ){
			fCut_DiLeptInvMass_min    = float(ParValue); ok = true;
		} 	else if( !strcmp(ParName, "DiLeptInvMass_max") ){
			fCut_DiLeptInvMass_max    = float(ParValue); ok = true;
		} else if( !strcmp(ParName, "DiLeptOSSFInvMass_lowercut") ){
			fCut_DiLeptOSSFInvMass_lowercut    = float(ParValue); ok = true;
		} else if( !strcmp(ParName, "DiLeptOSSFInvMass_uppercut") ){
			fCut_DiLeptOSSFInvMass_uppercut    = float(ParValue); ok = true;
		} else if( !strcmp(ParName, "PtHat_max")){
			fCut_PtHat_max                     = float(ParValue); ok = true;
		}		

		if(!ok) cout << "%% MultiplicityAnalysis::ReadCuts ==> ERROR: Unknown variable " << ParName << endl;
	}	
	if(verbose){
		cout << "setting cuts to: " << endl;
		cout << "  PFMET_min                   " << fCut_PFMET_min                  <<endl;
		cout << "  HT_min                      " << fCut_HT_min                     <<endl;
		cout << "  JPt_hardest_min             " << fCut_JPt_hardest_min            <<endl;
		cout << "  DiLeptInvMass_min           " << fCut_DiLeptInvMass_min          <<endl;
		cout << "  DiLeptInvMass_max           " << fCut_DiLeptInvMass_max          <<endl;		
		cout << "  DiLeptOSSFInvMass_lowercut  " << fCut_DiLeptOSSFInvMass_lowercut <<endl;
		cout << "  DiLeptOSSFInvMass_uppercut  " << fCut_DiLeptOSSFInvMass_uppercut <<endl;
		cout << "  PtHat_max                   " << fCut_PtHat_max                  <<endl;
		cout << "  Run_min                     " << fCut_Run_min                    <<endl;
		cout << "  Run_max                     " << fCut_Run_max                    <<endl;

		for(int i=0; i<fRequiredHLT.size(); ++i){
			cout << "  HLTRequired (logic OR)      " << fRequiredHLT[i]                  <<endl;
		}
		for(int i=0; i<fVetoedHLT.size(); ++i){
			cout << "  HLTVetoed                   " << fVetoedHLT[i]                    <<endl;
		}
		cout << "--------------"    << endl;	
	}			
}
