#include "helper/Utilities.hh"
#include "MultiplicityAnalysisBase.hh"
#include "TSystem.h"
#include "JetCorrectionUncertainty.h"

using namespace std;

MultiplicityAnalysisBase::MultiplicityAnalysisBase(TreeReader *tr) : UserAnalysisBase(tr){
	Util::SetStyle();	
	fCut_PFMET_min                      = 0;
	fCut_HT_min                         = 0;
	fCut_caloHT50_min                   = 0;
	fCut_caloHT50ID_min                 = 0;
	fCut_caloMHT30_min                  = 0;
	fCut_caloMHT30ID_min                = 0;
	fCut_JPt_hardest_min                = 0;
	fCut_JPt_second_min                 = 0;
	fCut_DiLeptInvMass_min              = 0;
	fCut_DiLeptInvMass_max              = 9999.99;
	fCut_PtHat_max                      = 999999.;
	fCut_Run_min                        = 0;
	fCut_Run_max                        = 9999999;
	fDoJESUncertainty                   = false;
	fJESUpDown                          = 0;

	fRequiredHLT.clear();
	fVetoedHLT.clear();
}

MultiplicityAnalysisBase::~MultiplicityAnalysisBase(){
	delete fJecUncPF;
	delete fJecUncCalo;
}



void MultiplicityAnalysisBase::InitializeEvent(){
	GetLeptonJetIndices();
	FindLeptonConfig();
}

void MultiplicityAnalysisBase::GetLeptonJetIndices(){
	fIsNANObj = false; 

	fElecs.clear();
	fMuons.clear();
	fTaus.clear();
	fJets.clear();
	fJetTaus.reset();

	vector<double> mutight;
	for(int i=0; i< fTR->PfMu3NObjs; ++i){
		if(std::isnan(fTR->PfMu3Pt[i]))  {fIsNANObj = true; continue;} //protection against objects with NAN-Pt
		fMuons.push_back(i);
		mutight.push_back(fTR->PfMu3Pt[i]);
	}
	fMuons      = Util::VSort(fMuons     , mutight);
	
	vector<double> eltight;
	for(int i=0; i< fTR->PfEl3NObjs; ++i){
		if(std::isnan(fTR->PfEl3Pt[i]))  {fIsNANObj = true; continue;} //protection against objects with NAN-Pt
		fElecs.push_back(i);
		eltight.push_back(fTR->PfEl3Pt[i]);
	}
	fElecs      = Util::VSort(fElecs     , eltight);
	
	vector<double> pt1; 
	for(int ij=0; ij < fTR->PF2PAT3NJets; ++ij){
		if(std::isnan(fTR->PF2PAT3JPt[ij])){ fIsNANObj = true; continue; } //protection against objects with NAN-Pt
		if(Jet(ij).Pt() < 20)       continue;  // note: ETH ntuple only stores PF2PAT3Jets > 15 GeV (defualt config)
		fJets.push_back(ij);                   // fJets has all jets except for duplicates with selected leptons
		pt1.push_back(Jet(ij).Pt());
		fJetTaus.index.push_back(ij); fJetTaus.pt.push_back(Jet(ij).Pt()); fJetTaus.isTau.push_back(0); fJetTaus.NObjs++;
	}
	fJets        = Util::VSort(fJets,       pt1);
	pt1.clear();

	// !!!!!!!!!
	// Don't add taus to jets starting from ntuple V02-01-01: jets include taus	
	// !!!!!!!!!
	vector<double> taus;
	for(int i=0; i< fTR->PfTau3NObjs; ++i){
		if(std::isnan(fTR->PfTau3Pt[i])) { fIsNANObj = true; continue;} //protection against objects with NAN-Pt
		if(fTR->PfTau3Pt[i]   < 20    ) continue; // note: taus go up to 2.5 in Eta
		fTaus.push_back(i);
		taus.push_back(fTR->PfTau3Pt[i]);
	}
	fTaus          = Util::VSort(fTaus     , taus);
	
	//sort fJetTaus accorting to Pt
	fJetTaus.index = Util::VSort(fJetTaus.index, fJetTaus.pt);
	fJetTaus.isTau = Util::VSort(fJetTaus.isTau, fJetTaus.pt);
	sort(fJetTaus.pt.begin(), fJetTaus.pt.end());
	reverse(fJetTaus.pt.begin(), fJetTaus.pt.end());

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
			charge1=fTR->PfEl3Charge[fElecs[0]];
			charge2=fTR->PfEl3Charge[fElecs[1]];
			if(charge1*charge2==1){fLeptConfig=SS_ee;}
			else{fLeptConfig=OS_ee;}
		} else if(fMuons.size()==2){
			charge1=fTR->PfMu3Charge[fMuons[0]];
			charge2=fTR->PfMu3Charge[fMuons[1]];
			if(charge1*charge2==1){fLeptConfig=SS_mumu;}
			else{fLeptConfig=OS_mumu;}
		} else{
			charge1=fTR->PfEl3Charge[fElecs[0]];
			charge2=fTR->PfMu3Charge[fMuons[0]];			
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

	// Protection against events with NAN-Pt objects
	if(fIsNANObj) {
		cout << "WARNING: Event " << fTR->Event << " Lumi " << fTR->LumiSection << " Run " << fTR->Run << " has NAN-Pt Objects!! Skipping Event" << endl;
		return false;
	}
	
	//PtHat
	if(fTR->PtHat > fCut_PtHat_max ){return false;}
	
	// MET
	if(MET().Pt() < fCut_PFMET_min){return false;}
	
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

	// HT from jets + taus
	double HT=0;
	for(int j=0; j<fJetTaus.NObjs; ++j){
		if(!fJetTaus.isTau[j]){
			if(Jet(fJetTaus.index[j]).Pt() > 50 && fabs(Jet(fJetTaus.index[j]).Eta())<3.0){
				HT += Jet(fJetTaus.index[j]).Pt();
			}
		}else {
			if(fTR->PfTau3Pt[fJetTaus.index[j]] > 50 && fabs(fTR->PfTau3Eta[fJetTaus.index[j]]) < 3.0){
				HT +=fTR->PfTau3Pt[fJetTaus.index[j]];
			}
		}
	}
	fHT = HT;
	if(HT<fCut_HT_min){return false;}
	
	//caloHT and calo MHT
	TLorentzVector mht30(0,0,0,0), mht30_ID(0,0,0,0);
	fCaloHT50 =0.0, fCaloHT50_ID =0.0;
	for(int j=0; j<fTR->CANJets; ++j){
		if( CAJet(j).Pt()<30 || fabs(CAJet(j).Eta())>3.0 ) continue;
	  	// MHT
		mht30 -= CAJet(j);
		// MHT_ID
		if(fTR->CAJn90[j]>=2 && fTR->CAJEMfrac[j]>=0.000001) mht30_ID -= CAJet(j);	
		// HT
		if( CAJet(j).Pt()<50 ) continue;
		fCaloHT50  += CAJet(j).Pt();
		// HT_ID
		if(fTR->CAJn90[j]>=2 && fTR->CAJEMfrac[j]>=0.000001) fCaloHT50_ID += CAJet(j).Pt();
	}
	fCaloMHT30   =mht30.Pt();
	fCaloMHT30_ID=mht30_ID.Pt();
	if(fCaloHT50      < fCut_caloHT50_min   ) return false;
	if(fCaloHT50_ID   < fCut_caloHT50ID_min ) return false;
	if(fCaloMHT30     < fCut_caloMHT30_min  ) return false;
	if(fCaloMHT30_ID  < fCut_caloMHT30ID_min) return false;

	// leading jets including JID for jets
	bool leadingjets(true);
	if(fCut_JPt_hardest_min > 0){
		if(fJetTaus.NObjs <1) leadingjets=false;
		else if(!fJetTaus.isTau[0]){
			if(! IsGoodBasicPFJetPAT3(fJetTaus.index[0], fCut_JPt_hardest_min, 2.4)) {leadingjets=false;}
		}else {
			if(fTR->PfTau3Pt[fJetTaus.index[0]]   < fCut_JPt_hardest_min || fabs(fTR->PfTau3Eta[fJetTaus.index[0]]) >2.4  ) {leadingjets=false;}
		}
	}
	if(fCut_JPt_second_min > 0){
		if(fJetTaus.NObjs <2) leadingjets=false;
		else if(!fJetTaus.isTau[1]){
			if(! IsGoodBasicPFJetPAT3(fJetTaus.index[1], fCut_JPt_second_min, 2.4)) {leadingjets=false;}
		}else {
			if(fTR->PfTau3Pt[fJetTaus.index[1]]   < fCut_JPt_second_min || fabs(fTR->PfTau3Eta[fJetTaus.index[1]]) >2.4  ) {leadingjets=false;}
		}
	}
	if(leadingjets == false) return false;
	
	// DiLeptonInvMass_min DiLeptonInvMass_max
	if(fLeptConfig==SS_ee || fLeptConfig==OS_ee || fLeptConfig==SS_mumu || fLeptConfig==OS_mumu ||
	   fLeptConfig==SS_emu || fLeptConfig==OS_emu ){
		double invmass = GetDiLeptInvMass();
		
		if(invmass < fCut_DiLeptInvMass_min) {return false;}
		if(invmass > fCut_DiLeptInvMass_max) {return false;}		 
				 
	}

	// flag crazy events (dedicated to filter HCAL "laser" events)
	if(fTR->HCALSumEt > 10000 && fTR->PF2PAT3NJets > 25 ){
		fCrazyHCAL =1;
		cout << "WARNING: crazy HCAL event: Run " << fTR->Run << " LS " << fTR->LumiSection  << " Event " << fTR->Event 
		     << " has HCALSumET " << fTR->HCALSumEt << " njets " << fTR->PF2PAT3NJets << endl;
	}
	else    fCrazyHCAL =0;	

	// flag events with Jets with correction factor from JE correction 
	fNegativeJEC =0; 
	for( int i=0; i<fTR->PF2PAT3NJets; ++i){
		if(fTR->PF2PAT3JScale[i]>=0) continue;
		fNegativeJEC =1;
	}
	if(fNegativeJEC) {
		cout << "WARNING: Event with Jets with negative JEC: Run " << fTR->Run << " LS " << fTR->LumiSection << " Event " << fTR->Event
		     << " N PFJets " << fTR->PF2PAT3NJets << " NCaloJets " << fTR->CANJets << endl; 
	}


	// ------------------------------------------------------------------------------------------	
	return true;	
}

double MultiplicityAnalysisBase::GetDiLeptInvMass(){
	TLorentzVector p1, p2;
	if(fLeptConfig==SS_ee || fLeptConfig==OS_ee){
		p1.SetPtEtaPhiM(fTR->PfEl3Pt[fElecs[0]],fTR->PfEl3Eta[fElecs[0]],fTR->PfEl3Phi[fElecs[0]], 0. );
		p2.SetPtEtaPhiM(fTR->PfEl3Pt[fElecs[1]],fTR->PfEl3Eta[fElecs[1]],fTR->PfEl3Phi[fElecs[1]], 0. );
	} else if(fLeptConfig==SS_mumu || fLeptConfig==OS_mumu ){
		p1.SetPtEtaPhiM(fTR->PfMu3Pt[fMuons[0]],fTR->PfMu3Eta[fMuons[0]],fTR->PfMu3Phi[fMuons[0]], 0.105 );
		p2.SetPtEtaPhiM(fTR->PfMu3Pt[fMuons[1]],fTR->PfMu3Eta[fMuons[1]],fTR->PfMu3Phi[fMuons[1]], 0.105 );
	} else if(fLeptConfig==SS_emu || fLeptConfig==OS_emu ){
		p1.SetPtEtaPhiM(fTR->PfMu3Pt[fMuons[0]],fTR->PfMu3Eta[fMuons[0]],fTR->PfMu3Phi[fMuons[0]], 0.105 );
		p2.SetPtEtaPhiM(fTR->PfEl3Pt[fElecs[0]],fTR->PfEl3Eta[fElecs[0]],fTR->PfEl3Phi[fElecs[0]], 0. );
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
		} else if( !strcmp(ParName, "JESUpDown")){
			fJESUpDown   = int(IntValue); ok = true;
		}

		// floats 
		sscanf(buffer, "%s %f", ParName, &ParValue);
		if( !strcmp(ParName, "PFMET_min") ){
			fCut_PFMET_min            = float(ParValue); ok = true;
		} else if( !strcmp(ParName, "HT_min") ){
			fCut_HT_min               = float(ParValue); ok = true;
		} else if( !strcmp(ParName, "caloHT50_min") ){
			fCut_caloHT50_min         = float(ParValue); ok = true;
		} else if( !strcmp(ParName, "caloHT50ID_min") ){
			fCut_caloHT50ID_min       = float(ParValue); ok = true;
		} else if( !strcmp(ParName, "caloMHT30_min") ){
			fCut_caloMHT30_min        = float(ParValue); ok = true;
		} else if( !strcmp(ParName, "caloMHT30ID_min") ){
			fCut_caloMHT30ID_min      = float(ParValue); ok = true;
		} else if( !strcmp(ParName, "JPt_hardest_min") ){
			fCut_JPt_hardest_min      = float(ParValue); ok = true;			
		} else if( !strcmp(ParName, "JPt_second_min") ){
			fCut_JPt_second_min      = float(ParValue); ok = true;			
		} else if( !strcmp(ParName, "DiLeptInvMass_min") ){
			fCut_DiLeptInvMass_min    = float(ParValue); ok = true;
		} else if( !strcmp(ParName, "DiLeptInvMass_max") ){
			fCut_DiLeptInvMass_max    = float(ParValue); ok = true;
		} else if( !strcmp(ParName, "PtHat_max")){
			fCut_PtHat_max            = float(ParValue); ok = true;
		}  		

		if(!ok) cout << "%% MultiplicityAnalysis::ReadCuts ==> ERROR: Unknown variable " << ParName << endl;
	}	
	if(verbose){
		cout << "setting cuts to: " << endl;
		cout << "  PFMET_min                   " << fCut_PFMET_min                  <<endl;
		cout << "  HT_min                      " << fCut_HT_min                     <<endl;
		cout << "  caloHT50_min                " << fCut_caloHT50_min               <<endl;
		cout << "  caloHT50ID_min              " << fCut_caloHT50ID_min             <<endl;
		cout << "  caloMHT30_min               " << fCut_caloMHT30_min              <<endl;
		cout << "  caloMHT30ID_min             " << fCut_caloMHT30ID_min            <<endl;
		cout << "  JPt_hardest_min             " << fCut_JPt_hardest_min            <<endl;
		cout << "  JPt_second_min              " << fCut_JPt_second_min             <<endl;
		cout << "  DiLeptInvMass_min           " << fCut_DiLeptInvMass_min          <<endl;
		cout << "  DiLeptInvMass_max           " << fCut_DiLeptInvMass_max          <<endl;		
		cout << "  PtHat_max                   " << fCut_PtHat_max                  <<endl;
		cout << "  Run_min                     " << fCut_Run_min                    <<endl;
		cout << "  Run_max                     " << fCut_Run_max                    <<endl;
		cout << "  JESUpDown                   " << fJESUpDown                      <<endl;
		if(fJESUpDown ==-1 || fJESUpDown==1){
		if(fJESUpDown ==1) cout << "  Do up  -scaling of Jets and MET "             <<endl;
		else               cout << "  Do down-scaling of Jets and MET "             <<endl;
		fDoJESUncertainty = true;
		}

		for(int i=0; i<fRequiredHLT.size(); ++i){
			cout << "  HLTRequired (logic OR)      " << fRequiredHLT[i]                  <<endl;
		}
		for(int i=0; i<fVetoedHLT.size(); ++i){
			cout << "  HLTVetoed                   " << fVetoedHLT[i]                    <<endl;
		}
		cout << "--------------"    << endl;	
	}			
	if(fDoJESUncertainty){
		cout << "--------------------- " << endl;
		cout << " -> initialize JetCorrectionUncertainty " << endl;
		Initialize_JetCorrectionUncertainty();
		cout << "--------------------- " << endl;
	}else{
		fJecUncPF   =NULL;
		fJecUncCalo =NULL;
	}
}


// ****************************************************************************************************

bool MultiplicityAnalysisBase::IsGoodBasicPFJetPAT3(int index, double ptcut, double absetacut){
	// Basic PF jet cleaning and ID cuts
	// cut at pt of ptcut (default = 30 GeV)
	// cut at abs(eta) of absetacut (default = 2.5)
	if(Jet(index).Pt() < ptcut                  ) return false;
	if(fabs(fTR->PF2PAT3JEta[index]) > absetacut) return false;
	if(fTR->PF2PAT3JIDLoose[index]    ==0       ) return false;
	if(fTR->PF2PAT3JScale[index]     < 0        ) return false;
	return true;
}

bool MultiplicityAnalysisBase::IsGoodPFJetMediumPAT3(int index, double ptcut, double absetacut) {
	// Medium PF JID
	if ( ! IsGoodBasicPFJetPAT3(index, ptcut, absetacut)  ) return false;
	if ( !(fTR->PF2PAT3JNeuHadfrac[index] < 0.95)         ) return false;
	if ( !(fTR->PF2PAT3JNeuEmfrac[index]  < 0.95)         ) return false;
	return true;
}

bool MultiplicityAnalysisBase::IsGoodPFJetTightPAT3(int index, double ptcut, double absetacut) {
	// Tight PF JID
	if ( ! IsGoodBasicPFJetPAT3(index, ptcut, absetacut)  ) return false;
	if ( !(fTR->PF2PAT3JNeuHadfrac[index] < 0.90)         ) return false;
	if ( !(fTR->PF2PAT3JNeuEmfrac[index]  < 0.90)         ) return false;
	return true;
}

// Jets and JES uncertainty
void MultiplicityAnalysisBase::Initialize_JetCorrectionUncertainty(){
	fJecUncCalo = new JetCorrectionUncertainty("Fall10_AK5Calo_Uncertainty.txt");
	fJecUncPF   = new JetCorrectionUncertainty("Fall10_AK5PF_Uncertainty.txt");
}
// pf
TLorentzVector MultiplicityAnalysisBase::Jet(int index){
	TLorentzVector j(0,0,0,0);
	j.SetPtEtaPhiE(fTR->PF2PAT3JPt[index], fTR->PF2PAT3JEta[index], fTR->PF2PAT3JPhi[index], fTR->PF2PAT3JE[index]);
	if  (! fDoJESUncertainty) return j;
	else                      return MultiplicityAnalysisBase::PFJetJESScaled(j);
}

TLorentzVector MultiplicityAnalysisBase::PFJetJESScaled(TLorentzVector j){
	TLorentzVector j_scaled(0.,0.,0.,0);
	if(fJESUpDown== 1) j_scaled.SetPtEtaPhiM(j.Perp()*(1+GetJECUncertPF(j.Pt(),j.Eta())), j.Eta(), j.Phi(), j.M());
	if(fJESUpDown==-1) j_scaled.SetPtEtaPhiM(j.Perp()*(1-GetJECUncertPF(j.Pt(),j.Eta())), j.Eta(), j.Phi(), j.M());
	return j_scaled;
}

double MultiplicityAnalysisBase::GetJECUncertPF(double pt, double eta){
	// pt must be corrected! not raw  	
	fJecUncPF->setJetPt(pt);
	fJecUncPF->setJetEta(eta);
	double uncert= fJecUncPF->getUncertainty(true);
	if     (fabs(eta) <1.5)  uncert = sqrt(uncert*uncert + 0.02*0.02);
	else if(fabs(eta) <3.0)  uncert = sqrt(uncert*uncert + 0.06*0.06);
	else if(fabs(eta) <5.0)  uncert = sqrt(uncert*uncert + 0.20*0.20);
	return uncert; 
}

// Calo
TLorentzVector MultiplicityAnalysisBase::CAJet(int index){
	TLorentzVector j(0,0,0,0);
	j.SetPtEtaPhiE(fTR->CAJPt[index], fTR->CAJEta[index], fTR->CAJPhi[index], fTR->CAJE[index]);
	if  (! fDoJESUncertainty) return j;
	else                      return MultiplicityAnalysisBase::CAJetJESScaled(j);
}


TLorentzVector MultiplicityAnalysisBase::CAJetJESScaled(TLorentzVector j){
	TLorentzVector j_scaled(0.,0.,0.,0);
	if(fJESUpDown== 1) j_scaled.SetPtEtaPhiM(j.Perp()*(1+GetJECUncertCalo(j.Pt(),j.Eta())), j.Eta(), j.Phi(), j.M());
	if(fJESUpDown==-1) j_scaled.SetPtEtaPhiM(j.Perp()*(1-GetJECUncertCalo(j.Pt(),j.Eta())), j.Eta(), j.Phi(), j.M());
	return j_scaled;
}

double MultiplicityAnalysisBase::GetJECUncertCalo(double pt, double eta){
	// pt must be corrected! not raw  	
	fJecUncCalo->setJetPt(pt);
	fJecUncCalo->setJetEta(eta);
	double uncert= fJecUncCalo->getUncertainty(true);
	if     (fabs(eta) <1.5)  uncert = sqrt(uncert*uncert + 0.02*0.02);
	else if(fabs(eta) <3.0)  uncert = sqrt(uncert*uncert + 0.06*0.06);
	else if(fabs(eta) <5.0)  uncert = sqrt(uncert*uncert + 0.20*0.20);
	return uncert; 
}

// pfMET 
TLorentzVector MultiplicityAnalysisBase::MET(){
	TLorentzVector MET(0,0,0,0);
	if(!fDoJESUncertainty){
		MET.SetPtEtaPhiM(fTR->PFMETPAT, 0., fTR->PFMETPATphi, 0);
		return MET;
	} else{
		if      (fJESUpDown==1){
			MET.SetPtEtaPhiM(fTR->PFMETPAT*1.05, 0., fTR->PFMETPATphi, 0);
		}else if(fJESUpDown==-1){
			MET.SetPtEtaPhiM(fTR->PFMETPAT*0.95, 0., fTR->PFMETPATphi, 0);
		}else{
			cout << " something wrong in met scaling" << endl; exit(1);
		}		
		return MET;
	}
}

