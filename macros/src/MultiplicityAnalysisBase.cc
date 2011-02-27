#include "helper/Utilities.hh"
#include "MultiplicityAnalysisBase.hh"


using namespace std;

MultiplicityAnalysisBase::MultiplicityAnalysisBase(TreeReader *tr) : UserAnalysisBase(tr){
	Util::SetStyle();	
	fCut_PFMET_min                      = 0;
	fCut_MHT_min                        = 0;
	fCut_HT_min                         = 0;
	fCut_JPt_hardest_min                = 0;
	fCut_JPt_second_min                 = 0;
	fCut_VSPT                           = 9999.99;
	fCut_DiLeptInvMass_min              = 0;
	fCut_DiLeptInvMass_max              = 9999.99;
	fCut_Zselector                      = 0;
	fCut_Zveto                          = 0;
	fCut_DiLeptOSSFInvMass_lowercut     = 9999.99;
	fCut_DiLeptOSSFInvMass_uppercut     =-9999.99;
	fCut_PtHat_max                      = 999999.;
	fCut_Run_min                        = 0;
	fCut_Run_max                        = 9999999;
	fCut_ElectronTrigger                = 0;
	fCut_MuonTrigger                    = 0;

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
	fJetsLoose.clear();
	fJetsMedium.clear();
	fJetsTight.clear();
	fBJets.clear();

	vector<double> mutight;
	for(int i=0; i< fTR->NMus; ++i){
		if(! IsGoodBasicMu(i)           ) continue;
		if(! (fTR->MuPt[i] > 10)        ) continue;
		if(! (fabs(fTR->MuEta[i]) < 2.4)) continue;
		double iso    = fTR->MuRelIso03[i];
		double pt     = fTR->MuPt[i];
		double hybiso = iso*pt / std::max(20.,pt);
		if(! (hybiso < 0.15)            ) continue;
		fMuons.push_back(i);
		mutight.push_back(fTR->MuPt[i]);
	}
	fMuons      = Util::VSort(fMuons     , mutight);
	
	vector<double> eltight;
	for(int i=0; i< fTR->NEles; ++i){
		if(! IsLooseEl(i)               ) continue;
		if(! (fTR->ElPt[i] > 10)        ) continue;
		if(! (fabs(fTR->ElEta[i]) < 2.4)) continue;
		if(! IsIsolatedEl(i, 0.15, 0.15)) continue; //hybiso 
		fElecs.push_back(i);
		eltight.push_back(fTR->ElPt[i]);
	}
	fElecs      = Util::VSort(fElecs     , eltight);

	vector<double> pt1;
	vector<double> pt2;
	vector<double> pfloose;
	vector<double> pfmedium;
	vector<double> pftight;

	fNJets_toremove_ele=0;
	fNJets_toremove_muo=0;
	for(int i=0; i<fElecs.size(); ++i){
		if(fTR->ElPt[fElecs[i]] > 15) fNJets_toremove_ele ++;
	}	
	for(int i=0; i<fMuons.size(); ++i){
		if(fTR->MuPt[fMuons[i]] > 15) fNJets_toremove_muo ++;
	}	

	bool doSel(true);
	for(int ij=0; ij < fTR->PFNJets; ++ij){
		if(fTR->PFJPt[ij] < 15) continue;  // note: ETH ntuple only stores PFJets > 15 GeV (defualt config)

		bool JGood(true);
		for(int i=0; i<fMuons.size(); ++i){
			double deltaR = Util::GetDeltaR(fTR->PFJEta[ij], fTR->MuEta[fMuons[i]], 
					fTR->PFJPhi[ij], fTR->MuPhi[fMuons[i]]);
			if(deltaR < 0.4)   {JGood=false; fNJets_toremove_muo--;}
		}
		for(int i=0; i<fElecs.size(); ++i){
			double deltaR = Util::GetDeltaR(fTR->PFJEta[ij], fTR->ElEta[fElecs[i]], 
					fTR->PFJPhi[ij], fTR->ElPhi[fElecs[i]]);
			if(deltaR < 0.4)   {JGood=false; fNJets_toremove_ele--;}
		}
		if(JGood==false) continue;
		
		fJets.push_back(ij);                 // fJets has all jets except for duplicates with selected leptons
		pt1.push_back(fTR->PFJPt[ij]);
		if(! IsGoodBasicPFJet(ij,  20., 2.4) ) continue;
		fJetsLoose.push_back(ij);
		pfloose.push_back(fTR->PFJPt[ij]);
		if( fTR->PFJbTagProbTkCntHighEff[ij] > 3.3){
			fBJets.push_back(ij);
			pt2.push_back(fTR->JPt[ij]);	
		}
		if(! IsGoodPFJetMedium(ij, 20., 2.4) ) continue;
		fJetsMedium.push_back(ij);
		pfmedium.push_back(fTR->PFJPt[ij]);
		if(! IsGoodPFJetTight(ij,  20., 2.4) ) continue;
		fJetsTight.push_back(ij);
		pftight.push_back(fTR->PFJPt[ij]);
	}
	fJets        = Util::VSort(fJets,       pt1);
	fBJets       = Util::VSort(fBJets,      pt2);
	fJetsLoose   = Util::VSort(fJetsLoose,  pfloose);
	fJetsMedium  = Util::VSort(fJetsMedium, pfmedium);
	fJetsTight   = Util::VSort(fJetsTight,  pftight);
	pfloose.clear();
	pfmedium.clear();
	pftight.clear();
	pt1.clear();
	pt2.clear();
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
	}else if(fElecs.size() + fMuons.size() >2){
		fLeptConfig=multilept;
	}

}


bool MultiplicityAnalysisBase::IsSelectedEvent(){
	// goodevt from UserAnalysisBase
	if(!IsGoodEvent()) {return false;}

	// Run
	if(fTR->Run < fCut_Run_min ) {return false;}
	if(fTR->Run > fCut_Run_max ) {return false;}
	
	//PtHat
	if(fTR->PtHat > fCut_PtHat_max ){return false;}
	
	// MET
	if(fTR->PFMET < fCut_PFMET_min){return false;}
	
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
//	if(fCut_ElectronTrigger == 1){
//		if(! IsGoodElEvent_RA5()) return false;
//	}
//	if(fCut_MuonTrigger ==1 ){
//		if(! IsGoodMuEvent()) return false;
//	}

	// HT
	double HT=0;
	for(int j=0; j<fJetsLoose.size(); ++j){
		if(fTR->PFJPt[fJetsLoose[j]] > 50 && fabs(fTR->PFJEta[fJetsLoose[j]])<2.4){
			HT += fTR->PFJPt[fJetsLoose[j]];
		}
	}
	fHT = HT;
	if(HT<fCut_HT_min){return false;}


	// leading jets
	bool leadingjets(true);
	if(fCut_JPt_hardest_min > 0){
		if(fJets.size() <1) leadingjets=false;
		else if(IsGoodBasicPFJet(fJets[0], fCut_JPt_hardest_min, 2.4) == false ){leadingjets=false;}
	}
	if(fCut_JPt_second_min > 0){
		if(fJets.size() <2) leadingjets=false;
		else if(IsGoodBasicPFJet(fJets[1], fCut_JPt_second_min,  2.4) == false ){leadingjets=false;}
	}
	if(leadingjets == false) return false;
	
	
	// MHT
	TVector3 MHTall(0., 0., 0.);
        for(int i=0; i<fTR->PFNJets; ++i) {
		TVector3 jet;
		jet.SetPtEtaPhi(fTR->PFJPt[i], fTR->PFJEta[i], fTR->PFJPhi[i]);
		MHTall += jet;
	}	
	fMHTall = MHTall.Pt();
	if(fMHTall < fCut_MHT_min) {return false;}

	
	// DiLeptonInvMass_min DiLeptonInvMass_max
	if(fLeptConfig==SS_ee || fLeptConfig==OS_ee || fLeptConfig==SS_mumu || fLeptConfig==OS_mumu ||
	   fLeptConfig==SS_emu || fLeptConfig==OS_emu ){
		double invmass = GetDiLeptInvMass();
		
		if(invmass < fCut_DiLeptInvMass_min) {return false;}
		if(invmass > fCut_DiLeptInvMass_max) {return false;}		 
				 
	}
	
	// Z-bosoon selector or veto
	if( fLeptConfig==OS_ee ||fLeptConfig==OS_mumu ){
		double invmass = GetDiLeptInvMass();
		if( fCut_Zselector==1 && (invmass > fCut_DiLeptOSSFInvMass_lowercut) && (invmass < fCut_DiLeptOSSFInvMass_uppercut) ) {return true; }	
		if( fCut_Zveto    ==1 && (invmass > fCut_DiLeptOSSFInvMass_lowercut) && (invmass < fCut_DiLeptOSSFInvMass_uppercut) ) {return false;}	
	} else if(fCut_Zselector ==1 ) {return false;}

	// VectorSumPt of selected & identified objects(el, mu, jet) and MET
	double px=0;
	double py=0;
	for(int i=0; i<fJetsLoose.size(); ++i){
		px+=fTR->PFJPx[fJetsLoose[i]];
		py+=fTR->PFJPy[fJetsLoose[i]];
	}
	fMHT = sqrt(px*px + py*py); // MHT from all selected jets
	TVector3 MHTvec(0.,0.,0.);
	MHTvec.SetXYZ(px, py, 0.);
	fMHTphi=MHTvec.Phi();

	for(int i=0; i<fMuons.size(); ++i){
		px+=fTR->MuPx[fMuons[i]];
		py+=fTR->MuPy[fMuons[i]];
	}
	for(int i=0; i<fElecs.size(); ++i){
		px+=fTR->ElPx[fElecs[i]];
		py+=fTR->ElPy[fElecs[i]];
	}
	px+=fTR->PFMETpx;
	py+=fTR->PFMETpy;

	fVectorSumPt = sqrt(px*px + py*py);

	if(fVectorSumPt > fCut_VSPT)    return false;

	// ------------------------------------------------------------------------------------------	
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
		} else if( !strcmp(ParName, "Zselector") ){
			fCut_Zselector = int(IntValue); ok = true;
		} else if( !strcmp(ParName, "Zveto") ){
			fCut_Zveto   = int(IntValue); ok = true;
		} else if( !strcmp(ParName, "ElectronTrigger") ){
			fCut_ElectronTrigger   = int(IntValue); ok = true;
		} else if( !strcmp(ParName, "MuonTrigger") ){
			fCut_MuonTrigger   = int(IntValue); ok = true;
		}

		// floats 
		sscanf(buffer, "%s %f", ParName, &ParValue);
		if( !strcmp(ParName, "PFMET_min") ){
			fCut_PFMET_min            = float(ParValue); ok = true;
		} else if( !strcmp(ParName, "MHT_min") ){
			fCut_MHT_min              = float(ParValue); ok = true;
		} else if( !strcmp(ParName, "HT_min") ){
			fCut_HT_min               = float(ParValue); ok = true;
		} else if( !strcmp(ParName, "JPt_hardest_min") ){
			fCut_JPt_hardest_min      = float(ParValue); ok = true;			
		} else if( !strcmp(ParName, "JPt_second_min") ){
			fCut_JPt_second_min      = float(ParValue); ok = true;			
		} else if( !strcmp(ParName, "DiLeptInvMass_min") ){
			fCut_DiLeptInvMass_min    = float(ParValue); ok = true;
		} else if( !strcmp(ParName, "DiLeptInvMass_max") ){
			fCut_DiLeptInvMass_max    = float(ParValue); ok = true;
		} else if( !strcmp(ParName, "DiLeptOSSFInvMass_lowercut") ){
			fCut_DiLeptOSSFInvMass_lowercut    = float(ParValue); ok = true;
		} else if( !strcmp(ParName, "DiLeptOSSFInvMass_uppercut") ){
			fCut_DiLeptOSSFInvMass_uppercut    = float(ParValue); ok = true;
		} else if( !strcmp(ParName, "PtHat_max")){
			fCut_PtHat_max                     = float(ParValue); ok = true;
		} else if( !strcmp(ParName, "VSPT_max")){
			fCut_VSPT                          = float(ParValue); ok = true;
		}  		

		if(!ok) cout << "%% MultiplicityAnalysis::ReadCuts ==> ERROR: Unknown variable " << ParName << endl;
	}	
	if(verbose){
		cout << "setting cuts to: " << endl;
		cout << "  PFMET_min                   " << fCut_PFMET_min                  <<endl;
		cout << "  MHT_min                     " << fCut_MHT_min                    <<endl;
		cout << "  HT_min                      " << fCut_HT_min                     <<endl;
		cout << "  JPt_hardest_min             " << fCut_JPt_hardest_min            <<endl;
		cout << "  JPt_second_min              " << fCut_JPt_second_min             <<endl;
		cout << "  VSPT_max                    " << fCut_VSPT                       <<endl;
		cout << "  DiLeptInvMass_min           " << fCut_DiLeptInvMass_min          <<endl;
		cout << "  DiLeptInvMass_max           " << fCut_DiLeptInvMass_max          <<endl;		
		cout << "  Zveto                       " << fCut_Zveto                      <<endl;
		cout << "  Zselector                   " << fCut_Zselector                  <<endl;
		cout << "  DiLeptOSSFInvMass_lowercut  " << fCut_DiLeptOSSFInvMass_lowercut <<endl;
		cout << "  DiLeptOSSFInvMass_uppercut  " << fCut_DiLeptOSSFInvMass_uppercut <<endl;
		cout << "  PtHat_max                   " << fCut_PtHat_max                  <<endl;
		cout << "  Run_min                     " << fCut_Run_min                    <<endl;
		cout << "  Run_max                     " << fCut_Run_max                    <<endl;
		cout << "  ElectronTrigger             " << fCut_ElectronTrigger            <<endl;
		cout << "  MuonTrigger                 " << fCut_MuonTrigger                <<endl;

		for(int i=0; i<fRequiredHLT.size(); ++i){
			cout << "  HLTRequired (logic OR)      " << fRequiredHLT[i]                  <<endl;
		}
		for(int i=0; i<fVetoedHLT.size(); ++i){
			cout << "  HLTVetoed                   " << fVetoedHLT[i]                    <<endl;
		}
		cout << "--------------"    << endl;	
	}			
}
