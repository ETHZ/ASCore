#include "helper/Utilities.hh"
#include "MassAnalysis.hh"
#include "TLorentzVector.h"
#include <sstream>

using namespace std;

MassAnalysis::MassAnalysis(TreeReader *tr) : MultiplicityAnalysisBase(tr){
	Util::SetStyle();	
}

MassAnalysis::~MassAnalysis(){
}

void MassAnalysis::Begin(const char* filename){
	// Define the output file of histograms
	fHistFile = new TFile(fOutputDir + TString(filename), "RECREATE");

	// book tree
	fMT2tree = new MT2tree();
	BookTree();

	// define which triggers to fill
	if(fisData){
		// HT
		fTriggerMap["HLT_HT150_v2"]            = &fMT2tree->trigger.HLT_HT150_v2;
		fTriggerMap["HLT_HT150_v3"]            = &fMT2tree->trigger.HLT_HT150_v3;
		fTriggerMap["HLT_HT160_v2"]            = &fMT2tree->trigger.HLT_HT160_v2;
		fTriggerMap["HLT_HT200_v2"]            = &fMT2tree->trigger.HLT_HT200_v2;
		fTriggerMap["HLT_HT200_v3"]            = &fMT2tree->trigger.HLT_HT200_v3;
		fTriggerMap["HLT_HT240_v2"]            = &fMT2tree->trigger.HLT_HT240_v2;
		fTriggerMap["HLT_HT250_v2"]            = &fMT2tree->trigger.HLT_HT250_v2;
		fTriggerMap["HLT_HT250_v3"]            = &fMT2tree->trigger.HLT_HT250_v3;
		fTriggerMap["HLT_HT260_v2"]            = &fMT2tree->trigger.HLT_HT260_v2;
		fTriggerMap["HLT_HT300_v2"]            = &fMT2tree->trigger.HLT_HT300_v2;
		fTriggerMap["HLT_HT300_v3"]            = &fMT2tree->trigger.HLT_HT300_v3;
		fTriggerMap["HLT_HT300_v4"]            = &fMT2tree->trigger.HLT_HT300_v4;
		fTriggerMap["HLT_HT300_v5"]            = &fMT2tree->trigger.HLT_HT300_v5;
		fTriggerMap["HLT_HT350_v2"]            = &fMT2tree->trigger.HLT_HT350_v2;
		fTriggerMap["HLT_HT350_v3"]            = &fMT2tree->trigger.HLT_HT350_v3;
		fTriggerMap["HLT_HT350_v4"]            = &fMT2tree->trigger.HLT_HT350_v4;
		fTriggerMap["HLT_HT360_v2"]            = &fMT2tree->trigger.HLT_HT360_v2;
		fTriggerMap["HLT_HT400_v2"]            = &fMT2tree->trigger.HLT_HT400_v2;
		fTriggerMap["HLT_HT400_v3"]            = &fMT2tree->trigger.HLT_HT400_v3;
		fTriggerMap["HLT_HT400_v4"]            = &fMT2tree->trigger.HLT_HT400_v4;
		fTriggerMap["HLT_HT400_v5"]            = &fMT2tree->trigger.HLT_HT400_v5;
		fTriggerMap["HLT_HT400_v6"]            = &fMT2tree->trigger.HLT_HT400_v6;
		fTriggerMap["HLT_HT400_v7"]            = &fMT2tree->trigger.HLT_HT400_v7;
		fTriggerMap["HLT_HT440_v2"]            = &fMT2tree->trigger.HLT_HT440_v2;
		fTriggerMap["HLT_HT450_v2"]            = &fMT2tree->trigger.HLT_HT450_v2;
		fTriggerMap["HLT_HT450_v3"]            = &fMT2tree->trigger.HLT_HT450_v3;
		fTriggerMap["HLT_HT450_v4"]            = &fMT2tree->trigger.HLT_HT450_v4;
		fTriggerMap["HLT_HT450_v5"]            = &fMT2tree->trigger.HLT_HT450_v5;
		fTriggerMap["HLT_HT450_v6"]            = &fMT2tree->trigger.HLT_HT450_v6;
		fTriggerMap["HLT_HT450_v7"]            = &fMT2tree->trigger.HLT_HT450_v7;
		fTriggerMap["HLT_HT500_v2"]            = &fMT2tree->trigger.HLT_HT500_v2;
		fTriggerMap["HLT_HT500_v3"]            = &fMT2tree->trigger.HLT_HT500_v3;
		fTriggerMap["HLT_HT500_v4"]            = &fMT2tree->trigger.HLT_HT500_v4;
		fTriggerMap["HLT_HT500_v5"]            = &fMT2tree->trigger.HLT_HT500_v5;
		fTriggerMap["HLT_HT500_v6"]            = &fMT2tree->trigger.HLT_HT500_v6;
		fTriggerMap["HLT_HT500_v7"]            = &fMT2tree->trigger.HLT_HT500_v7;
		fTriggerMap["HLT_HT550_v2"]            = &fMT2tree->trigger.HLT_HT550_v2;
		fTriggerMap["HLT_HT550_v3"]            = &fMT2tree->trigger.HLT_HT550_v3;
		fTriggerMap["HLT_HT550_v4"]            = &fMT2tree->trigger.HLT_HT550_v4;
		fTriggerMap["HLT_HT550_v5"]            = &fMT2tree->trigger.HLT_HT550_v5;
		fTriggerMap["HLT_HT550_v6"]            = &fMT2tree->trigger.HLT_HT550_v6;
		fTriggerMap["HLT_HT550_v7"]            = &fMT2tree->trigger.HLT_HT550_v7;
		// MHT_HT
		fTriggerMap["HLT_HT250_MHT60_v2"]      = &fMT2tree->trigger.HLT_HT250_MHT60_v2;
		fTriggerMap["HLT_HT250_MHT60_v3"]      = &fMT2tree->trigger.HLT_HT250_MHT60_v3;
		fTriggerMap["HLT_HT250_MHT60_v4"]      = &fMT2tree->trigger.HLT_HT250_MHT60_v4;
		fTriggerMap["HLT_HT250_MHT70_v1"]      = &fMT2tree->trigger.HLT_HT250_MHT70_v1;
		fTriggerMap["HLT_HT260_MHT60_v2"]      = &fMT2tree->trigger.HLT_HT260_MHT60_v2;
		fTriggerMap["HLT_HT300_MHT75_v4"]      = &fMT2tree->trigger.HLT_HT300_MHT75_v4;
		fTriggerMap["HLT_HT300_MHT75_v5"]      = &fMT2tree->trigger.HLT_HT300_MHT75_v5;
		// QuadJet
		fTriggerMap["HLT_QuadJet50_BTagIP_v1"] = &fMT2tree->trigger.HLT_QuadJet50_BTagIP_v1;
		fTriggerMap["HLT_QuadJet50_Jet40_v1"]  = &fMT2tree->trigger.HLT_QuadJet50_Jet40_v1;
		// Muons
		fTriggerMap["HLT_DoubleMu3_HT160_v2"]  = &fMT2tree->trigger.HLT_DoubleMu3_HT160_v2;
		fTriggerMap["HLT_DoubleMu3_v3"]        = &fMT2tree->trigger.HLT_DoubleMu3_v3;
		fTriggerMap["HLT_Mu8_Jet40_v2"]        = &fMT2tree->trigger.HLT_Mu8_Jet40_v2;
	}


	// initialize fDeadCellFilterBE and fDeadCellFilterTP
	fDeadCellFilterBE.Reset();
	fDeadCellFilterTP.Reset();
	for(int i=0; i<fTPfiles.size(); ++i){
		MassAnalysis::DeadCellParser(fDeadCellFilterTP, fTPfiles[i]);
	}
	for(int i=0; i<fBEfiles.size(); ++i){
		MassAnalysis::DeadCellParser(fDeadCellFilterBE, fBEfiles[i]);
	}


}

void MassAnalysis::Analyze(){	

	// ---------------------------------------------------
	// Initialize fElecs, fJetsLoose, fBJets, fMuons, fLeptConfig 
	InitializeEvent();

	// check if event passes selection cuts and trigger requirements
	if(! IsSelectedEvent()){return;}
	// --------------------------------------------------------------------

	// reset tree variables
	ResetTree();

	// -----------------------------------------------------------
	// fill fHemiObjects
	for(int i=0; i<gNHemispheres; ++i){ fHemiObjects[i].Reset();}
	   // ************************************************* //
	   // changes here affect misc.MT2 type variables !!!!  //
	   // ************************************************  //

	// seed: max inv mass, assoc: minimal lund distance, no dR cut (i.e. third parameter is 0)
	GetMT2Variables(2, 3, 0.0, 20, 2.4, 0,  fHemiObjects[0]); 
	
	// seed: max inv mass, assoc: minimal lund distance, no dR cut (i.e. third parameter is 0), JETID enforced
	GetMT2Variables(2, 3, 0.0, 20, 2.4, 1,  fHemiObjects[1]); 
	
	// minimizing deltaHT
	bool minimizeDHT(true);
	GetMT2Variables(minimizeDHT, 20, 2.4,   fHemiObjects[2]); 
	
	// seed: max inv mass, assoc: minimmal inv mass, JET ID required
	GetMT2Variables(2, 2, 0.0, 20, 2.4, 1,  fHemiObjects[3]); 
	
	// -----------------------------------------------------------


	// Fill Tree !! has to be called at the very end !!
	MassAnalysis::FillTree();
	
	// print interesting events
//	InterestingEvents();
}

// ***********************************************************************************************
void MassAnalysis::BookTree(){
	          
	fATree = new TTree("MassTree", "MassTree");
	fATree->Branch("MT2tree" , "MT2tree" , &fMT2tree);

}

void MassAnalysis::ResetTree(){
        fMT2tree->Reset();
}

void MassAnalysis::FillTree(){
	// check size of jets electrons muons and genleptons
	if(fJetTaus.NObjs > 40) {cout << "ERROR: fJetTaus.NObjs > 40: " << "run " << fTR->Run << " Event " << fTR->Event << " skip event" << endl; return;}
	if(fElecs.size()  > 5 ) {cout << "ERROR: fElecs.size()  >  5: " << "run " << fTR->Run << " Event " << fTR->Event << " skip event" << endl; return;}
	if(fMuons.size()  > 5 ) {cout << "ERROR: fMuons.size()  >  5: " << "run " << fTR->Run << " Event " << fTR->Event << " skip event" << endl; return;}

	// ---------------------------------------------------------------
	// Fill jets 4-momenta & ID's 
	for(int i=0; i<fJetTaus.NObjs; ++i) {
		if(! fJetTaus.isTau[i]){
			fMT2tree->jet[i].lv.SetPtEtaPhiE( Jet(fJetTaus.index[i]).Pt(),  Jet(fJetTaus.index[i]).Eta(), 
						          Jet(fJetTaus.index[i]).Phi(), Jet(fJetTaus.index[i]).E()  ); //SetLV(GetJet4Momenta(fJets[i]));
			// b-tag info now should be available
			fMT2tree->jet[i].bTagProbTCHE  =  fTR->PF2PAT3JbTagProbTkCntHighEff [fJetTaus.index[i]];
			fMT2tree->jet[i].bTagProbTCHP  =  fTR->PF2PAT3JbTagProbTkCntHighPur [fJetTaus.index[i]];
			fMT2tree->jet[i].bTagProbSSVHE =  fTR->PF2PAT3JbTagProbSimpSVHighEff[fJetTaus.index[i]];
			fMT2tree->jet[i].bTagProbSSVHP =  fTR->PF2PAT3JbTagProbSimpSVHighPur[fJetTaus.index[i]];
		  
			// Jet id variables
			if(IsGoodBasicPFJetPAT3 (fJetTaus.index[i], 20., 2.4))  fMT2tree->jet[i].isPFIDLoose   =true;
			if(IsGoodPFJetMediumPAT3(fJetTaus.index[i], 20., 2.4))  fMT2tree->jet[i].isPFIDMedium  =true;
			if(IsGoodPFJetTightPAT3 (fJetTaus.index[i], 20., 2.4))  fMT2tree->jet[i].isPFIDTight   =true;
			if(fTR->PF2PAT3JIDLoose [fJetTaus.index[i]]==1       )  fMT2tree->jet[i].isPATPFIDLoose=true;;  // PF jetID regardless of eta and pt
			fMT2tree->jet[i].ChHadFrac      = fTR->PF2PAT3JChHadfrac     [fJetTaus.index[i]];	
			fMT2tree->jet[i].NeuHadFrac     = fTR->PF2PAT3JNeuHadfrac    [fJetTaus.index[i]];
			fMT2tree->jet[i].ChEmFrac       = fTR->PF2PAT3JChEmfrac      [fJetTaus.index[i]];
			fMT2tree->jet[i].NeuEmFrac      = fTR->PF2PAT3JNeuEmfrac     [fJetTaus.index[i]];
			fMT2tree->jet[i].ChMult         = fTR->PF2PAT3JChMult        [fJetTaus.index[i]];
			fMT2tree->jet[i].NeuMult        = fTR->PF2PAT3JNeuMult       [fJetTaus.index[i]];
			fMT2tree->jet[i].NConstituents  = fTR->PF2PAT3JNConstituents [fJetTaus.index[i]];
			fMT2tree->jet[i].Scale          = fTR->PF2PAT3JScale         [fJetTaus.index[i]];
			fMT2tree->jet[i].L1FastJetScale = fTR->PF2PAT3JL1FastJetScale[fJetTaus.index[i]];
			fMT2tree->jet[i].Area           = fTR->PF2PAT3JArea          [fJetTaus.index[i]];
			fMT2tree->jet[i].isTau          = false;
			if(!fisData){
			fMT2tree->jet[i].Flavour        = fTR->PF2PAT3JFlavour       [fJetTaus.index[i]];
			}
		}else { // this is obsolete starting from ntuple V02-01-01 as taus are included in jets
			fMT2tree->jet[i].lv.SetPxPyPzE( fTR->PfTau3Px[fJetTaus.index[i]],fTR->PfTau3Py[fJetTaus.index[i]],fTR->PfTau3Pz[fJetTaus.index[i]],fTR->PfTau3E[fJetTaus.index[i]]);
			fMT2tree->jet[i].isTau = true;
			cout << "ERROR: there is something wrong: jet should not be a tau: double counting?" << endl;
			exit(1);
		}
	}

	// --------------------------------------------------------------
	// match taus to jets
	for(int i=0; i<fTaus.size(); ++i){
		TLorentzVector tau;
		tau.SetPxPyPzE(fTR->PfTau3Px[fTaus[i]],fTR->PfTau3Py[fTaus[i]],fTR->PfTau3Pz[fTaus[i]],fTR->PfTau3E[fTaus[i]]);
		double mindR =10000; int jindex =-1;
		for(int j=0; j<fJets.size(); ++j){
			double dR = tau.DeltaR(fMT2tree->jet[j].lv);
			if(dR < mindR && dR < 0.5) {mindR=dR; jindex = j;} 	
		}
		if(jindex ==-1) continue;
		fMT2tree->jet[jindex].NTauMatch++;
		if(fMT2tree->jet[jindex].isTauMatch && fMT2tree->jet[jindex].TauDR < mindR) continue;
		fMT2tree->jet[jindex].isTauMatch = 1;
		fMT2tree->jet[jindex].TauDR      = mindR;
		fMT2tree->jet[jindex].TauDPt     = fMT2tree->jet[jindex].lv.Pt()-tau.Pt();
	}

	
	// ---------------------------------------------------------------
	// Set NJets, NElecs, NMuons
	fMT2tree->SetNJets         ((Int_t)fJetTaus.NObjs);
	fMT2tree->SetNGenJets      (fTR->NGenJets > gNGenJets ? gNGenJets: fTR->NGenJets);
	fMT2tree->SetNJetsIDLoose  (fMT2tree->GetNjets(20, 2.4, 1));
	fMT2tree->SetNJetsAcc      (fMT2tree->GetNjets(20, 2.4, 0));
	fMT2tree->SetNBJets        (fMT2tree->GetNBtags(3,2.0,20,2.4,1));
	fMT2tree->SetNEles         ((Int_t)fElecs.size());
	fMT2tree->SetNMuons        ((Int_t)fMuons.size());
	fMT2tree->SetNTaus         ((Int_t)fTaus.size());
	
	// --------------------------------------------------------------
	// Fill GenJets
	for(int i=0; i<fTR->NGenJets; ++i){
		if(i >= gNGenJets) {
			cout << "WARNING: Event " << fTR->Event << " Run " << fTR->Run << " has more than " << gNGenJets << " GenJets " << endl;
			continue;
	       	}
		fMT2tree->genjet[i].lv.SetPtEtaPhiE(fTR->GenJetPt[i], fTR->GenJetEta[i], fTR->GenJetPhi[i], fTR->GenJetE[i]);
		double mindR=999.99;
		int    index=-1;
		for(int j=0; j<fMT2tree->NJets; ++j){
			double dR=fMT2tree->jet[j].lv.DeltaR(fMT2tree->genjet[i].lv);
			if(dR < mindR) {
				mindR = dR;
				index = j;
			}
		}
		fMT2tree->genjet[i].DeltaR        = mindR;
		fMT2tree->genjet[i].JetMatchIndex = index;
	}
	
	
	// -----------------------------------------------------------------
	// Fill leptons 4-momenta & tight_flag & charge
	TLorentzVector METlv;
	METlv.SetPtEtaPhiM(MET().Pt(), 0., MET().Phi(), 0.);

	for(int i=0; i<fElecs.size(); ++i) {
	  	fMT2tree->ele[i].lv.SetPtEtaPhiE(fTR->PfEl3Pt [fElecs[i]], fTR->PfEl3Eta[fElecs[i]], 
				          fTR->PfEl3Phi[fElecs[i]], fTR->PfEl3E  [fElecs[i]]); // = GetEle4Momenta(fElecs[i]);
		fMT2tree->ele[i].MT     = GetMT(fMT2tree->ele[i].lv, 0., METlv, 0.); 
		fMT2tree->ele[i].Charge = fTR->PfEl3Charge[fElecs[i]];
		fMT2tree->ele[i].Iso    = (fTR->PfEl3ChargedHadronIso[fElecs[i]] + fTR->PfEl3NeutralHadronIso[fElecs[i]] + fTR->PfEl3PhotonIso[fElecs[i]])/fTR->PfEl3Pt[fElecs[i]];
		fMT2tree->ele[i].ID90   = fTR->PfEl3ID90[fElecs[i]];
		fMT2tree->ele[i].ID95   = fTR->PfEl3ID95[fElecs[i]];
		
		// match to caloJets
		for(int j=0; j<fTR->CANJets; ++j){
			if(fMT2tree->ele[i].lv.DeltaR(CAJet(j)) <0.4) {
				fMT2tree->ele[i].CAJ_n90     = fTR->CAJn90[j];
				fMT2tree->ele[i].CAJ_n90Hits = fTR->CAJID_n90Hits[j];
			}
		}
	}
	for(int i=0; i<fMuons.size(); ++i) {
	  	fMT2tree->muo[i].lv.SetPtEtaPhiM(fTR->PfMu3Pt [fMuons[i]], fTR->PfMu3Eta[fMuons[i]], 
				          fTR->PfMu3Phi[fMuons[i]], 0.106);                     // = GetMuo4Momenta(fMuons[i]);
		fMT2tree->muo[i].MT       = GetMT(fMT2tree->muo[i].lv, fMT2tree->muo[i].lv.M(), METlv, 0.); 
		fMT2tree->muo[i].Charge   = fTR->PfMu3Charge[fMuons[i]];	
		fMT2tree->muo[i].Iso      = (fTR->PfMu3ChargedHadronIso[fMuons[i]] + fTR->PfMu3NeutralHadronIso[fMuons[i]] + fTR->PfMu3PhotonIso[fMuons[i]])/fTR->PfMu3Pt[fMuons[i]];
		fMT2tree->muo[i].NMatches = fTR->PfMu3NMatches[fMuons[i]];
		fMT2tree->muo[i].PtErr    = fTR->PfMu3PtErr[fMuons[i]];
	}
	
	// ----------------------------------------------------------------
	// fill MT2Hemi
	for(int h=0; h<gNHemispheres; ++h){
		fMT2tree->hemi[h].assoc_method  = fHemiObjects[h].assoc;
		fMT2tree->hemi[h].seed_method   = fHemiObjects[h].seed;
		fMT2tree->hemi[h].maxDR         = fHemiObjects[h].maxDR;
		fMT2tree->hemi[h].MT2           = fHemiObjects[h].MT2;
		fMT2tree->hemi[h].MCT           = fHemiObjects[h].MCT;
		fMT2tree->hemi[h].AlphaT        = fHemiObjects[h].alphaT;
		fMT2tree->hemi[h].minDHT        = fHemiObjects[h].minDHT;
		fMT2tree->hemi[h].dPhi          = fHemiObjects[h].dPhi;
		fMT2tree->hemi[h].lv1           = fHemiObjects[h].pjet1;
		fMT2tree->hemi[h].lv2           = fHemiObjects[h].pjet2;
		fMT2tree->hemi[h].UTM           = fHemiObjects[h].UTM;

		int jcount1=0, jcount2=0, elecount1=0, elecount2=0, muocount1=0, muocount2=0; 
		for(int i=0; i<fHemiObjects[h].objects.size(); ++i){
			HemiObject obj=fHemiObjects[h].objects[i];
			if(obj.type=="jet" && obj.hemi==1){
				fMT2tree->hemi[h].jindices1[jcount1]=obj.index;
				jcount1++;
			} else if(obj.type=="jet" && obj.hemi==2){
				fMT2tree->hemi[h].jindices2[jcount2]=obj.index;
				jcount2++;
			} else if(obj.type=="ele" && obj.hemi==1){
				fMT2tree->hemi[h].eleindices1[elecount1]=obj.index;
				elecount1++;
			} else if(obj.type=="ele" && obj.hemi==2){
				fMT2tree->hemi[h].eleindices2[elecount2]=obj.index;
				elecount2++;
			} else if(obj.type=="muo" && obj.hemi==1){
				fMT2tree->hemi[h].muoindices1[muocount1]=obj.index;
				muocount1++;
			} else if(obj.type=="muo" && obj.hemi==2){
				fMT2tree->hemi[h].muoindices2[muocount2]=obj.index;
				muocount2++;
			}
		}
	}
	
	

	// ---------------------------------------------------------------
	// GenMET	
	fMT2tree->genmet[0].SetPtEtaPhiM(fTR->GenMET, 0., fTR->GenMETphi, 0.);
	// -------------------------------------------------------------------
	// Genleptons
	int NGenLepts=0;
	for(int i=0; i<fTR->NGenLeptons; ++i){
		double mass=0;
		if     (abs(fTR->GenLeptonID[i]) == 15)   mass=1.776; // tau
		else if(abs(fTR->GenLeptonID[i]) == 13)   mass=0.106; // mu
		else if(abs(fTR->GenLeptonID[i]) == 11)   mass=0.;    // el 
		else if(abs(fTR->GenLeptonID[i]) == 12 || 
			abs(fTR->GenLeptonID[i]) == 14 || 
			abs(fTR->GenLeptonID[i]) == 16)   mass=0.;    // nu 
		else if(abs(fTR->GenLeptonID[i]) == 5 )   mass=4.2;   // bottom-quark
		else   continue;
		NGenLepts++;
		if(i >= 30 ) {cout << "ERROR: NGenLepts >=30: skipping remaining genlepts for event " << fTR->Event << endl; continue;}
		fMT2tree->genlept[i].lv.SetPtEtaPhiM(fTR->GenLeptonPt[i], fTR->GenLeptonEta[i], fTR->GenLeptonPhi[i], mass);
		fMT2tree->genlept[i].ID       = fTR->GenLeptonID[i];
		fMT2tree->genlept[i].MID      = fTR->GenLeptonMID[i];
		fMT2tree->genlept[i].MStatus  = fTR->GenLeptonMStatus[i];
		fMT2tree->genlept[i].GMID     = fTR->GenLeptonGMID[i];
		fMT2tree->genlept[i].GMStatus = fTR->GenLeptonGMStatus[i];
		if(abs(fMT2tree->genlept[i].ID) == 11 || abs(fMT2tree->genlept[i].ID) == 13 || abs(fMT2tree->genlept[i].ID) == 15   ){
			fMT2tree->genlept[i].MT = GetMT(fMT2tree->genlept[i].lv, fMT2tree->genlept[i].lv.M(), fMT2tree->genmet[0], 0.);
		}
		if(abs(fTR->GenLeptonID[i]) == 11){// match to caloJets
			for(int j=0; j<fTR->CANJets; ++j){
				if(fMT2tree->genlept[i].lv.DeltaR(CAJet(j)) <0.4) {
					fMT2tree->genlept[i].CAJ_n90     = fTR->CAJn90[j];
					fMT2tree->genlept[i].CAJ_n90Hits = fTR->CAJID_n90Hits[j];
				}
			}
		}
	}
	fMT2tree->NGenLepts = NGenLepts;

	// --------------------------------------------------------------------
	// MET and MPT
	fMT2tree->pfmet[0]=MET();

	// Fill vector sum of tracks
	TVector3 tracks(0.,0.,0.);
	for(int i=0; i< fTR->NTracks; ++i){
		TVector3 track;
		track.SetPtEtaPhi(fabs(fTR->TrkPt[i]), fTR->TrkEta[i], fTR->TrkPhi[i]);
		tracks += track;
	}	
	fMT2tree->MPT[0].SetXYZM(-tracks.Px(), -tracks.Py(), 0, 0);

	// Pile UP info and reco vertices
	if(!fisData){
		fMT2tree->pileUp.PUnumInt          = fTR->PUnumInteractions;
		fMT2tree->pileUp.PUnumIntLate      = fTR->PUOOTnumInteractionsLate;
		fMT2tree->pileUp.PUnumIntEarly     = fTR->PUOOTnumInteractionsEarly;
		fMT2tree->pileUp.PtHat             = fTR->PtHat;
		fMT2tree->pileUp.isS3              = (int) isS3;
		//////////// S3 vs S4
		if(isS3)
		  fMT2tree->pileUp.Weight            = GetPUWeight(fTR->PUnumInteractions, fTR->PUOOTnumInteractionsLate);
		else
		  fMT2tree->pileUp.Weight            = GetPUWeight(fTR->PUnumInteractions);

		fMT2tree->pileUp.Rho               = fTR->Rho;
		
		if(fVerbose > 3) {
			cout << "fTR->PUnumInteractions " <<  fTR->PUnumInteractions << " weight "  
		     	     << " GetPUWeight() " << GetPUWeight(fTR->PUnumInteractions) << endl; 
		}
	}
	int nvertex=0;
	for(int i=0; i<fTR->NVrtx; ++i){
		if(fabs(fTR->VrtxZ[i]) > 24) continue;
		if(sqrt( (fTR->VrtxX[i])*(fTR->VrtxX[i]) + (fTR->VrtxY[i])*(fTR->VrtxY[i])) > 2) continue;
		if(fTR->VrtxNdof[i]<=4) continue;
		nvertex++;
	}
	fMT2tree->pileUp.NVertices=nvertex;
	// _________

	// _________
	// HLT triggers
	for (StringBoolMap::iterator iter = fTriggerMap.begin(); iter != fTriggerMap.end(); ++iter){
		if(GetHLTResult(iter->first)) *iter->second =1;
	}
	// _________
	
	//__________
	// ECAL Dead Cell
	for(int i=0;i<fDeadCellFilterBE.event.size(); ++i){
		if(fTR->Run         !=fDeadCellFilterBE.run[i] )  continue;
		if(fTR->LumiSection !=fDeadCellFilterBE.lumi[i])  continue;
		if(fTR->Event       !=fDeadCellFilterBE.event[i]) continue;
		fMT2tree->misc.BadEcalBE = 1;
	}
	for(int i=0;i<fDeadCellFilterTP.event.size(); ++i){
		if(fTR->Run         !=fDeadCellFilterTP.run[i] )  continue;
		if(fTR->LumiSection !=fDeadCellFilterTP.lumi[i])  continue;
		if(fTR->Event       !=fDeadCellFilterTP.event[i]) continue;
		fMT2tree->misc.BadEcalTP = 1;
	}
	if(fTR->EcalDeadTPFilterFlag==0) fMT2tree->misc.BadEcalTP=1;
	
	// ------------------------------------------------------------------
	// fill misc 
	fMT2tree->misc.isData              = fisData;
	fMT2tree->misc.Run                 = fTR->Run;
	fMT2tree->misc.Event		   = fTR->Event;
	fMT2tree->misc.LumiSection	   = fTR->LumiSection;
	fMT2tree->misc.LeptConfig          = (int) fLeptConfig;
	fMT2tree->misc.HBHENoiseFlag	   = fTR->HBHENoiseFlag;
	fMT2tree->misc.CrazyHCAL           = fCrazyHCAL;
	fMT2tree->misc.NegativeJEC         = fNegativeJEC;
	fMT2tree->misc.HT                  = fHT;
	fMT2tree->misc.PFMETsign	   = (MET().Pt())/sqrt(fTR->SumEt);
	
	fMT2tree->misc.MT2                 = fHemiObjects[1].MT2;    // note: this is a bit dangerous, 
	fMT2tree->misc.MT2all              = fHemiObjects[0].MT2;    
	fMT2tree->misc.MCT                 = fHemiObjects[1].MCT;
	fMT2tree->misc.AlphaT              = fHemiObjects[2].alphaT; // minimizing deltaHT 

	fMT2tree->misc.MET                 = MET().Pt();
	fMT2tree->misc.METPhi              = MET().Phi();

	fMT2tree->misc.LeadingJPt          = (fMT2tree->NJets > 0) ? fMT2tree->jet[0].lv.Pt() : 0;
	fMT2tree->misc.SecondJPt           = (fMT2tree->NJets > 1) ? fMT2tree->jet[1].lv.Pt() : 0;


	// ----------------------------------------------------------------------------------
	// warning: these MT2tree methods are only to be called once all needed variables are filled!
	// warning: hardcoded values!
	fMT2tree->misc.Vectorsumpt	   = fMT2tree->GetMHTminusMET(1, 20, 2.4, true); // including leptons, ID jets only
	fMT2tree->misc.VectorsumptAll	   = fMT2tree->GetMHTminusMET(0, 20, 2.4, true); // including leptons
	fMT2tree->misc.MinMetJetDPhi       = fMT2tree->MinMetJetDPhi(0,20,5.0,1);
	fMT2tree->misc.PassJetID           = fMT2tree->PassJetID(50,2.4,1);
	fMT2tree->misc.PassJetID20         = fMT2tree->PassJetID(20,2.4,1);
	if(fMT2tree->NJets > 0) {
		fMT2tree->misc.Jet0Pass      = (Int_t) fMT2tree->jet[0].IsGoodPFJet(100,2.4,1);
	} else  fMT2tree->misc.Jet0Pass      = 0; 
	if(fMT2tree->NJets > 1) {
		fMT2tree->misc.Jet1Pass      = (Int_t) fMT2tree->jet[1].IsGoodPFJet( 60,2.4,1);
	} else  fMT2tree->misc.Jet1Pass      = 0;
	
	// MHT from jets and taus
	fMT2tree->MHTloose[0]=fMT2tree->GetMHTlv(1, 20, 2.4, true); // only jets satisfying the loose PF-ID and leptons
	fMT2tree->MHT     [0]=fMT2tree->GetMHTlv(0, 20, 2.4, true); // jets and leptons
	
	// DPhiMhtMpt: should be replaced by method to calculate it on the fly. 
	fMT2tree->misc.DPhiMhtMpt=Util::DeltaPhi(fMT2tree->MHT[0].Phi(), fMT2tree->MPT[0].Phi());	
	

	// ---------------------------	
	// calo HT and MHT
	float caHT40=0, caHT40_ID=0;
	TLorentzVector mht40(0,0,0,0), mht40_ID(0,0,0,0);
	for(int j=0; j<fTR->CANJets; ++j){
	  	if( CAJet(j).Pt()<40 || fabs(CAJet(j).Eta())>3.0 ) continue;

	    	caHT40   += CAJet(j).Pt();
	    	mht40    -= CAJet(j);
	  	
	  	if(! (fTR->CAJn90[j]>=2 && fTR->CAJEMfrac[j]>=0.000001)) continue;
		caHT40_ID += CAJet(j).Pt();
		mht40_ID  -= CAJet(j);
	}
	fMT2tree->misc.caloHT40     = caHT40;
	fMT2tree->misc.caloHT40_ID  = caHT40_ID;
	fMT2tree->misc.caloMHT40    = mht40.Pt();
	fMT2tree->misc.caloMHT40_ID = mht40_ID.Pt();
	fMT2tree->misc.caloHT50     = fCaloHT50;
	fMT2tree->misc.caloHT50_ID  = fCaloHT50_ID;
	fMT2tree->misc.caloMHT30    = fCaloMHT30;
	fMT2tree->misc.caloMHT30_ID = fCaloMHT30_ID;

	// RA2 tracking failure
	fMT2tree->misc.TrackingFailure     = fTR->TrkPtSum/fHT;
	fMT2tree->misc.TrackingFailurePVtx = fTR->PrimVtxPtSum/fHT;

	// _________
	// stuff for Z->nunu (close your eyes...)	
	vector<int>    jindi;
	vector<double> jpt;
	bool PassJetID_matched(true);
	double HTmatched=0, vectorsumpt_matched_px=0, vectorsumpt_matched_py=0;
	int    NJetsIDLoose_matched=0;
	for(int i=0; i<fTR->PF2PAT3NJets; ++i){
		if(Jet(i).Pt() < 20) continue;  
		bool jet(true);
		for(int gen=0; gen<fMT2tree->NGenLepts; ++gen){
			if( ((abs(fMT2tree->genlept[gen].ID) == 11 || abs(fMT2tree->genlept[gen].ID)==13) && fMT2tree->genlept[gen].MID ==23) ) {
				double deltaR = Util::GetDeltaR(Jet(i).Eta(), fMT2tree->genlept[gen].lv.Eta(), Jet(i).Phi(), fMT2tree->genlept[gen].lv.Phi());
				if(deltaR < 0.4) jet=false;
			}
		}
		if(  jet == false) continue;	
		if(Jet(i).Pt() > 50 && fabs(Jet(i).Eta())<2.4 && IsGoodBasicPFJetPAT3(i,  50., 2.4)==false ){
			PassJetID_matched  = false;
		}
		jindi.push_back(i); jpt.push_back(Jet(i).Pt());
		if(! IsGoodBasicPFJetPAT3(i,  20., 2.4) ) continue;
		vectorsumpt_matched_px+=Jet(i).Px();
		vectorsumpt_matched_py+=Jet(i).Py();
		NJetsIDLoose_matched++;
		if(Jet(i).Pt()<50 || fabs(Jet(i).Eta())>3 ) continue;
		HTmatched += Jet(i).Pt();
	}
	TLorentzVector met     = fMT2tree->pfmet[0];
	for(int gen=0; gen<fMT2tree->NGenLepts; ++gen){
		if( ((abs(fMT2tree->genlept[gen].ID) == 11 || abs(fMT2tree->genlept[gen].ID)==13) && fMT2tree->genlept[gen].MID ==23) ) {
			vectorsumpt_matched_px+=fMT2tree->genlept[gen].lv.Px();
			vectorsumpt_matched_py+=fMT2tree->genlept[gen].lv.Py();
			met +=fMT2tree->genlept[gen].lv;
		}
	}
	
	float           caHT50_matched      =0.0; 
	float           caHT50ID_matched    =0.0; 
	float           caHT50_matchedReco  =0.0;
	float           caHT50ID_matchedReco=0.0;
	TLorentzVector  mht30_matched       (0,0,0,0); 
	TLorentzVector  mht30ID_matched     (0,0,0,0); 
	TLorentzVector  mht30_matchedReco   (0,0,0,0);
	TLorentzVector  mht30ID_matchedReco (0,0,0,0);
	for(int j=0; j<fTR->CANJets; ++j){
	  	if( CAJet(j).Pt()<30 || fabs(CAJet(j).Eta())>3.0 ) continue;
		bool jet(true);
		for(int gen=0; gen<fMT2tree->NGenLepts; ++gen){
			if( ((abs(fMT2tree->genlept[gen].ID) == 11 || abs(fMT2tree->genlept[gen].ID)==13) && fMT2tree->genlept[gen].MID ==23) ) {
				double deltaR = Util::GetDeltaR(CAJet(j).Eta(), fMT2tree->genlept[gen].lv.Eta(), CAJet(j).Phi(), fMT2tree->genlept[gen].lv.Phi());
				if(deltaR < 0.4) jet=false;
			}
		}
		if(  jet == false) continue;	
	  	mht30_matched    -= CAJet(j);  // MHT30
		if(fTR->CAJn90[j]>=2 && fTR->CAJEMfrac[j]>=0.000001){
		mht30ID_matched  -= CAJet(j); // MHT30_ID
		} 
	  	if(CAJet(j).Pt()>50) {
			caHT50_matched    += CAJet(j).Pt(); //HT50
			if(fTR->CAJn90[j]>=2 && fTR->CAJEMfrac[j]>=0.000001) {
			caHT50ID_matched  += CAJet(j).Pt(); //HT50
			}
		}
	}
	for(int j=0; j<fTR->CANJets; ++j){
	  	if( CAJet(j).Pt()<30 || fabs(CAJet(j).Eta())>3.0 ) continue;
		bool jet(true);
		// remove overlap from reco electrons and muons
		for(int e=0; e<fMT2tree->NEles; ++e){
			double dR=Util::GetDeltaR(CAJet(j).Eta(), fMT2tree->ele[e].lv.Eta(), CAJet(j).Phi(), fMT2tree->ele[e].lv.Phi());
			if(dR < 0.4) jet=false;
		}
		for(int m=0; m<fMT2tree->NMuons; ++m){
			double dR=Util::GetDeltaR(CAJet(j).Eta(), fMT2tree->muo[m].lv.Eta(), CAJet(j).Phi(), fMT2tree->muo[m].lv.Phi());
			if(dR < 0.4) jet=false;
		}
		if(jet==false) continue;
	  	mht30_matchedReco   -= CAJet(j);  // MHT30
		if(fTR->CAJn90[j]>=2 && fTR->CAJEMfrac[j]>=0.000001){
		mht30ID_matchedReco -= CAJet(j);  // MHT30_ID
		}
	  	if(CAJet(j).Pt()>50) {
			caHT50_matchedReco    += CAJet(j).Pt(); //HT50
			if(fTR->CAJn90[j]>=2 && fTR->CAJEMfrac[j]>=0.000001){
			caHT50ID_matchedReco  += CAJet(j).Pt(); //HT50ID
			}
		}
	}


	fMT2tree->Znunu.caloMHT30_matched      =mht30_matched.Pt();
	fMT2tree->Znunu.caloMHT30ID_matched    =mht30ID_matched.Pt();
	fMT2tree->Znunu.caloMHT30_matchedReco  =mht30_matchedReco.Pt();
	fMT2tree->Znunu.caloMHT30ID_matchedReco=mht30ID_matchedReco.Pt();
	fMT2tree->Znunu.caloHT50_matched       =caHT50_matched;
	fMT2tree->Znunu.caloHT50ID_matched     =caHT50ID_matched;
	fMT2tree->Znunu.caloHT50_matchedReco   =caHT50_matchedReco;
	fMT2tree->Znunu.caloHT50ID_matchedReco =caHT50ID_matchedReco;

	double mindPhi=10;
	if(jindi.size()<1){mindPhi = -999.99;}
	else{
		for(int i=0; i<jindi.size(); ++i){
			if(Jet(jindi[i]).Pt()       < 20 ) continue;
			if(fabs(Jet(jindi[i]).Eta())> 5.0) continue;
			double dphi = TMath::Abs(Jet(jindi[i]).DeltaPhi(met));
			if(dphi < mindPhi){ 
				mindPhi = dphi;
			}
		}
	}
	if(mindPhi==10.){
		fMT2tree->Znunu.MinMetplusLeptJetDPhi      = -999.99;
	} else  fMT2tree->Znunu.MinMetplusLeptJetDPhi      =  mindPhi;
	fMT2tree->Znunu.MinMetplusLeptJetDPhiReco          =  fMT2tree->MinMetJetDPhi(0, 20, 5.0 ,3);

	jindi   = Util::VSort(jindi, jpt);
	if(jindi.size() >0){
		fMT2tree->Znunu.Jet0Pass_matched   = (Int_t) IsGoodBasicPFJetPAT3(jindi[0],  100., 2.4);
		fMT2tree->Znunu.LeadingJPt_matched = jpt[0];
	} else  fMT2tree->Znunu.Jet0Pass_matched   =0; 
	if(jindi.size() >1){
		fMT2tree->Znunu.Jet1Pass_matched   = (Int_t) IsGoodBasicPFJetPAT3(jindi[1],   60., 2.4);
		fMT2tree->Znunu.SecondJPt_matched  = jpt[1];
	} else  fMT2tree->Znunu.Jet1Pass_matched   =0;
	fMT2tree->Znunu.PassJetID_matched          = (Int_t) PassJetID_matched;
	fMT2tree->Znunu.Vectorsumpt_matched        = sqrt( pow(vectorsumpt_matched_px+MET().Px(),2) + pow(vectorsumpt_matched_py+MET().Py(),2));
	
	fMT2tree->Znunu.HTmatched                  = HTmatched;
	fMT2tree->Znunu.NJetsIDLoose_matched       = NJetsIDLoose_matched;

	fMT2tree->Znunu.NJetsToRemoveEle           = fNJets_toremove_ele;
	fMT2tree->Znunu.NJetsToRemoveMuo           = fNJets_toremove_muo;

	fMT2tree->Znunu.RecoOSee_mll               = fMT2tree->GetDiLeptonInvMass(0, 1, 1, 5.0, 1); 
	fMT2tree->Znunu.RecoOSmumu_mll             = fMT2tree->GetDiLeptonInvMass(0, 1, 2, 5.0, 1); 

	fMT2tree->Znunu.GenZee_mll                 = fMT2tree->GenOSDiLeptonInvMass(11,23,0,100);
	fMT2tree->Znunu.GenZee_mll_acc             = fMT2tree->GenOSDiLeptonInvMass(11,23,5.0,2.4);
	fMT2tree->Znunu.GenZmumu_mll               = fMT2tree->GenOSDiLeptonInvMass(13,23,0,100);
	fMT2tree->Znunu.GenZmumu_mll_acc           = fMT2tree->GenOSDiLeptonInvMass(13,23,5.0,2.4);

	fMT2tree->Znunu.GenZnunu_e_mll             = fMT2tree->GenOSDiLeptonInvMass(12,23,0,100);
	fMT2tree->Znunu.GenZnunu_e_mll_acc         = fMT2tree->GenOSDiLeptonInvMass(12,23,5.0,2.4);
	fMT2tree->Znunu.GenZnunu_mu_mll            = fMT2tree->GenOSDiLeptonInvMass(14,23,0,100);
	fMT2tree->Znunu.GenZnunu_mu_mll_acc        = fMT2tree->GenOSDiLeptonInvMass(14,23,5.0,2.4);
	fMT2tree->Znunu.GenZnunu_tau_mll           = fMT2tree->GenOSDiLeptonInvMass(16,23,0,100);
	fMT2tree->Znunu.GenZnunu_tau_mll_acc       = fMT2tree->GenOSDiLeptonInvMass(16,23,5.0,2.4);

	fMT2tree->Znunu.METplusLeptsPt             = fMT2tree->GetMETPlusGenLepts(0, 1, 1,  1113, 23, 0, 100, 0, 10000);
	fMT2tree->Znunu.METplusLeptsPtReco         = fMT2tree->GetMETPlusLepts(1);
	
	
	// ----------------------------------------------------------------------
	// fill tree
	fATree            ->Fill();

}

void MassAnalysis::InterestingEvents(){

 	if(fMT2tree->misc.MT2 < 300             ) return; 	
	if(fMT2tree->NJetsIDLoose < 2           ) return; 	
	if(fMT2tree->misc.MET < 30              ) return; 	
	if(fMT2tree->misc.HT < 300              ) return; 	
	if(fMT2tree->misc.Jet0Pass == 0         ) return; 	
	if(fMT2tree->misc.Jet1Pass == 0         ) return; 	
	if(fMT2tree->misc.PassJetID == 0        ) return; 	
	if(fMT2tree->misc.Vectorsumpt>70        ) return; 	
	if(fMT2tree->misc.MinMetJetDPhi<0.3     ) return; 	
	if(fMT2tree->misc.HBHENoiseFlag == 0    ) return; 	
	if(fMT2tree->NEles + fMT2tree->NMuons!=0) return;

	cout << "-------------------------------------------------------------------------------" << endl;	
	EventPrint();
	cout << "++++ MT2 Prtinouts: ++++ " << endl;
	cout << "NEles " << fMT2tree->NEles << endl;
	for(int i=0; i<fMT2tree->NEles; ++i){
		cout << " pt " << fMT2tree->ele[i].lv.Pt()  << " eta " << fMT2tree->ele[i].lv.Eta() << " phi " << fMT2tree->ele[i].lv.Phi() << endl;
	}
	cout << "  MT2 "                  << fMT2tree->hemi[0].MT2     << endl;
	cout << "   association method: " << fMT2tree->hemi[0].assoc_method << " seed method " << fMT2tree->hemi[0].seed_method << endl;
	cout << "   pseudojet 1: pt "     << fMT2tree->hemi[0].lv1.Pt()     << " eta "         << fMT2tree->hemi[0].lv1.Eta()   << " phi " << fMT2tree->hemi[0].lv1.Phi() << endl;
	cout << "   pseudojet 2: pt "     << fMT2tree->hemi[0].lv2.Pt()     << " eta "         << fMT2tree->hemi[0].lv2.Eta()   << " phi " << fMT2tree->hemi[0].lv2.Phi() << endl;
	cout << "   pseudojets dPhi "     << fMT2tree->hemi[0].dPhi         << endl;
	cout << "   sqrt(2*pj1*pj2*(1+cosTheta)) " << sqrt(2*fMT2tree->hemi[0].lv1.Pt()*fMT2tree->hemi[0].lv2.Pt() *(1+cos(fMT2tree->hemi[0].dPhi))) << endl;
	cout << "   |vec{MET} - vec{MHT}| " << fMT2tree->misc.Vectorsumpt   << endl;
	cout << "   minMetJetDPhi "          << fMT2tree->misc.MinMetJetDPhi << endl;
	cout << "   jets in hemi 1: "     << endl;
	for(int i=0; i<10; ++i){
		if(fMT2tree->hemi[0].jindices1[i]!=-1){
			cout << "     index " << fMT2tree->hemi[0].jindices1[i] << endl;
			cout << "     jet " << fMT2tree->hemi[0].jindices1[i] << " pt " << fMT2tree->jet[fMT2tree->hemi[0].jindices1[i]].lv.Pt() 
			     <<     " eta " << fMT2tree->jet[fMT2tree->hemi[0].jindices1[i]].lv.Eta() 
		             <<     " phi " << fMT2tree->jet[fMT2tree->hemi[0].jindices1[i]].lv.Phi()
		             <<     " Mass "<< fMT2tree->jet[fMT2tree->hemi[0].jindices1[i]].lv.M() ;  	     
			if(fMT2tree->jet[fMT2tree->hemi[0].jindices1[i]].isTauMatch) {cout << " jet is matched to Tau ";} 
			cout << endl;	     
	       	}
	}
	cout << "   jets in hemi 2: "     << endl;
	for(int i=0; i<10; ++i){
		if(fMT2tree->hemi[0].jindices2[i]!=-1){
			cout << "     index " << fMT2tree->hemi[0].jindices2[i] << endl;
			cout << "     jet " << fMT2tree->hemi[0].jindices2[i] << " pt " << fMT2tree->jet[fMT2tree->hemi[0].jindices2[i]].lv.Pt() 
			     <<     " eta " << fMT2tree->jet[fMT2tree->hemi[0].jindices2[i]].lv.Eta() 
		             <<     " phi " << fMT2tree->jet[fMT2tree->hemi[0].jindices2[i]].lv.Phi()
		             <<     " Mass "<< fMT2tree->jet[fMT2tree->hemi[0].jindices2[i]].lv.M() ; 
			if(fMT2tree->jet[fMT2tree->hemi[0].jindices2[i]].isTauMatch) {cout << " jet is matched to Tau ";}
		        cout << endl;	     
	       	}
	}
	cout << "   electrons in hemi 1: " << endl;
	for(int i=0; i<5; ++i){
		if(fMT2tree->hemi[0].eleindices1[i]==-1) continue;
		cout << "     ele " << fMT2tree->hemi[0].eleindices1[i] << " pt " << fMT2tree->ele[fMT2tree->hemi[0].eleindices1[i]].lv.Pt() 
		     <<     " eta " << fMT2tree->ele[fMT2tree->hemi[0].eleindices1[i]].lv.Eta() 
		     <<     " phi " << fMT2tree->ele[fMT2tree->hemi[0].eleindices1[i]].lv.Phi() << endl;
	}
	cout << "   electrons in hemi 2: " << endl;
	for(int i=0; i<5; ++i){
		if(fMT2tree->hemi[0].eleindices2[i]==-1) continue;
		cout << "     ele " << fMT2tree->hemi[0].eleindices2[i] << " pt " << fMT2tree->ele[fMT2tree->hemi[0].eleindices2[i]].lv.Pt() 
		     <<     " eta " << fMT2tree->ele[fMT2tree->hemi[0].eleindices2[i]].lv.Eta() 
		     <<     " phi " << fMT2tree->ele[fMT2tree->hemi[0].eleindices2[i]].lv.Phi() << endl;
	}
	cout << "   muons in hemi 1: " << endl;
	for(int i=0; i<5; ++i){
		if(fMT2tree->hemi[0].muoindices1[i]==-1) continue;
		cout << "     muo " << fMT2tree->hemi[0].muoindices1[i] << " pt " << fMT2tree->muo[fMT2tree->hemi[0].muoindices1[i]].lv.Pt() 
		     <<     " eta " << fMT2tree->muo[fMT2tree->hemi[0].muoindices1[i]].lv.Eta() 
		     <<     " phi " << fMT2tree->muo[fMT2tree->hemi[0].muoindices1[i]].lv.Phi() << endl;
	}
	cout << "   muons in hemi 2: " << endl;
	for(int i=0; i<5; ++i){
		if(fMT2tree->hemi[0].muoindices2[i]==-1) continue;
		cout << "     muo " << fMT2tree->hemi[0].muoindices2[i] << " pt " << fMT2tree->muo[fMT2tree->hemi[0].muoindices2[i]].lv.Pt() 
		     <<     " eta " << fMT2tree->muo[fMT2tree->hemi[0].muoindices2[i]].lv.Eta() 
		     <<     " phi " << fMT2tree->muo[fMT2tree->hemi[0].muoindices2[i]].lv.Phi() << endl;
	}
	cout << "  genlepton info " << endl;
	for(int i=0; i<fMT2tree->NGenLepts; ++i){
		if(abs(fMT2tree->genlept[i].MID) == 23 || abs(fMT2tree->genlept[i].MID) == 24 || abs(fMT2tree->genlept[i].MID) == 15 ){
			cout << "    genlepton with PID " << fMT2tree->genlept[i].ID  << " with Mother " << fMT2tree->genlept[i].MID << " and GrandMother " << fMT2tree->genlept[i].GMID 
			     <<              " with Pt  " << fMT2tree->genlept[i].lv.Pt() << " eta " << fMT2tree->genlept[i].lv.Eta() << " phi "  << fMT2tree->genlept[i].lv.Phi() << endl;
		}
	}
	cout << "-------------------------------------------------------------------------------" << endl;	
}

// ************************************************************************************************************************************* 
void MassAnalysis::GetMT2Variables(int hemi_seed, int hemi_assoc, double maxDR, double minJPt, double maxJEta, int PFJID, HemiObjects& hemiobject){

	// clear hemiobject
	hemiobject.objects .clear();
	hemiobject.pjet1   .SetPxPyPzE(0.,0.,0.,0.);
	hemiobject.pjet2   .SetPxPyPzE(0.,0.,0.,0.);
	hemiobject.UTM     .SetPxPyPzE(0.,0.,0.,0.);
	hemiobject.MT2     =-999.99;
	hemiobject.MCT     =-999.99;

	// fill stuff
	hemiobject.assoc   = hemi_assoc;
	hemiobject.seed    = hemi_seed;
	hemiobject.maxDR   = maxDR;
	hemiobject.minDHT  =-999.99;
	hemiobject.alphaT  =-999.99;
	hemiobject.dPhi    =-999.99;

	// protection against crazy HCAL events 
	if(fCrazyHCAL ) return;

	// make pseudojets with hemispheres
	vector<float> px, py, pz, E;
	for(int i=0; i<fJetTaus.NObjs; ++i){
		if(!fJetTaus.isTau[i]){
			if(PFJID==1 && !IsGoodBasicPFJetPAT3(fJetTaus.index[i], 20, 2.4)                             )  continue;
			else if(fabs(Jet(fJetTaus.index[i]).Eta())>2.4 || Jet(fJetTaus.index[i]).Pt()< 20            )  continue;
			if(fTR->PF2PAT3JScale[fJetTaus.index[i]] <0){
				cout << "WARNING: MT2 calculation failed: selected Jet with negative JEC: Run " << fTR->Run<< " Lumi: " << fTR->LumiSection << " Event " << fTR->Event 
			 	     <<": accoc: " << hemi_assoc << " seed " << hemi_seed << " PFJID " << PFJID << endl;
				return;
			}
			px.push_back(Jet(fJetTaus.index[i]).Px());
			py.push_back(Jet(fJetTaus.index[i]).Py());
			pz.push_back(Jet(fJetTaus.index[i]).Pz());
			 E.push_back(Jet(fJetTaus.index[i]).E() );
			 fHemiObject.type="jet"; fHemiObject.index=i; fHemiObject.hemi=0;
			 hemiobject.objects.push_back(fHemiObject);
		}else {
			if(fabs(fTR->PfTau3Eta[fJetTaus.index[i]])  >2.4 || fTR->PfTau3Pt[fJetTaus.index[i]]< 20) continue;
			px.push_back(fTR->PfTau3Px[fJetTaus.index[i]]);
			py.push_back(fTR->PfTau3Py[fJetTaus.index[i]]);
			pz.push_back(fTR->PfTau3Pz[fJetTaus.index[i]]);
			 E.push_back(fTR->PfTau3E [fJetTaus.index[i]]);
			 fHemiObject.type="jet"; fHemiObject.index=i; fHemiObject.hemi=0;
			 hemiobject.objects.push_back(fHemiObject);
		}
	}
	for(int i=0; i< fElecs.size(); ++i){
		px.push_back(fTR->PfEl3Px[fElecs[i]]);
		py.push_back(fTR->PfEl3Py[fElecs[i]]);
		pz.push_back(fTR->PfEl3Pz[fElecs[i]]);
		 E.push_back(fTR->PfEl3E [fElecs[i]]);
		 fHemiObject.type="ele"; fHemiObject.index=i; fHemiObject.hemi=0;
		 hemiobject.objects.push_back(fHemiObject);
	}
	for(int i=0; i< fMuons.size(); ++i){
		px.push_back(fTR->PfMu3Px[fMuons[i]]);
		py.push_back(fTR->PfMu3Py[fMuons[i]]);
		pz.push_back(fTR->PfMu3Pz[fMuons[i]]);
		 E.push_back(fTR->PfMu3E [fMuons[i]]);
		 fHemiObject.type="muo"; fHemiObject.index=i; fHemiObject.hemi=0;
		 hemiobject.objects.push_back(fHemiObject);
	}
	if(px.size()<2) {
		return; // protection against events with only one jet
	}	

	// get hemispheres 
	Hemisphere* hemi      = new Hemisphere(px, py, pz, E, hemi_seed, hemi_assoc);
	if(maxDR>0) hemi      -> RejectISRDRmax(maxDR);  // consider only jets withing dR < maxDR from hemi-axix for pseudojets
	hemi                  -> SetDebug(0);

	vector<int> grouping  = hemi ->getGrouping();
	
	TLorentzVector pseudojet1(0.,0.,0.,0.);
	TLorentzVector pseudojet2(0.,0.,0.,0.);

	double dHT = 0;	
	int counter=0;
	for(int i=0; i<px.size(); ++i){
		if(grouping[i]==1){
			pseudojet1.SetPx(pseudojet1.Px() + px[i]);
			pseudojet1.SetPy(pseudojet1.Py() + py[i]);
			pseudojet1.SetPz(pseudojet1.Pz() + pz[i]);
			pseudojet1.SetE( pseudojet1.E()  + E[i]);
			hemiobject.objects[i].hemi=1;
			dHT += sqrt(px[i]*px[i] + py[i]*py[i]);	
			counter ++;

		}else if(grouping[i] == 2){
			pseudojet2.SetPx(pseudojet2.Px() + px[i]);
			pseudojet2.SetPy(pseudojet2.Py() + py[i]);
			pseudojet2.SetPz(pseudojet2.Pz() + pz[i]);
			pseudojet2.SetE( pseudojet2.E()  + E[i]);
			hemiobject.objects[i].hemi=2;
			dHT -= sqrt(px[i]*px[i] + py[i]*py[i]);	
			counter ++;
		}
	}
	delete hemi;
	if(pseudojet1.Pt()==0 || pseudojet2.Pt()==0 ){
		cout << "++ WARNING : hemi seed " <<  hemi_seed  << " hemi assoc " << hemi_assoc << " n input " << grouping.size() << " n objets in hemi " << counter 
		     << " pseudojet1 pt " << pseudojet1.Pt() << " pseudojet2 pt " << pseudojet2.Pt() << " Event " << fTR->Event  << " Run " << fTR->Run << endl;
		return;
	}

	TLorentzVector pmiss;
	pmiss.SetPx(MET().Px());
	pmiss.SetPy(MET().Py());

	// fill mindHT (minimized for assoc method 9)
	hemiobject.minDHT = dHT;

	// fill UTM
	hemiobject.UTM          = - pmiss - pseudojet1 - pseudojet2;

	// fill MT2
	hemiobject.MT2          = GetMT2(pseudojet1, 0., pseudojet2, 0., pmiss, 0.);
	
	// fill MCT
	TVector2 pmiss_vector2;
	pmiss_vector2.Set(MET().Px(), MET().Py());
	TLorentzVector downstream(0.,0.,0.,0.); // no particles are downstream, i.e. not selected jets are upstream. 
	hemiobject.MCT   =GetMCTcorr(pseudojet1, pseudojet2, downstream, pmiss_vector2);
	
	// fill pseudojets
	hemiobject.pjet1 =pseudojet1;
	hemiobject.pjet2 =pseudojet2;
	hemiobject.dPhi  =Util::DeltaPhi(pseudojet1.Phi(), pseudojet2.Phi());

	// --------------------------------
	// testing: comparing to ToveyMT2	
	fMCT = new TMctLib;
	double MT2Tovey  =fMCT->mt2(pseudojet1, pseudojet2, downstream, pmiss_vector2, 14000, 0);
	delete fMCT;
	// --------------------------------
}

void MassAnalysis::GetMT2Variables(bool minimizeDHT,  double minJPt, double maxJEta, HemiObjects& hemiobject){

	// FIXME: save info about which jet is in which hemisphere e.g. fHemiObject.hemi 
	if(minimizeDHT !=1){ cout << "MassAnalysis::GetMT2Variables: ERROR: got minimizeDHT!=1" << endl; return; }

	// clear hemiobject
	hemiobject.objects .clear();
	hemiobject.pjet1   .SetPxPyPzE(0.,0.,0.,0.);
	hemiobject.pjet2   .SetPxPyPzE(0.,0.,0.,0.);
	hemiobject.UTM     .SetPxPyPzE(0.,0.,0.,0.);
	hemiobject.MT2     =-999.99;
	hemiobject.MCT     =-999.99;
	hemiobject.maxDR   =-999.99;
	hemiobject.dPhi    =-999.99;

	// fill stuff
	hemiobject.assoc   = -1;
	hemiobject.seed    = -1;

	// protection against crazy HCAL events 
	if(fCrazyHCAL ) return;

	// make pseudojets with hemispheres
	vector<TLorentzVector> p4s;
	for(int i=0; i<fJetTaus.NObjs; ++i){
		if(!fJetTaus.isTau[i]){
			if(fabs(Jet(fJetTaus.index[i]).Eta())>2.4 || Jet(fJetTaus.index[i]).Pt() < minJPt) continue;
			if(fTR->PF2PAT3JScale[fJetTaus.index[i]] <0){
				cout << "WARNING: MT2 minimizeDHT calculation failed: selected Jet with negative JEC: Run " << fTR->Run<< " Lumi: " << fTR->LumiSection << " Event " << fTR->Event <<endl;
				return;
			}
			p4s.push_back(Jet(fJetTaus.index[i]));
			fHemiObject.type="jet"; fHemiObject.index=i; fHemiObject.hemi=-1;
			hemiobject.objects.push_back(fHemiObject);
		}else {
			if(fabs(fTR->PfTau3Eta[fJetTaus.index[i]])  >2.4 || fTR->PfTau3Pt[fJetTaus.index[i]] < minJPt) continue;
			TLorentzVector v; 
			v.SetPxPyPzE(fTR->PfTau3Px[fJetTaus.index[i]], fTR->PfTau3Py[fJetTaus.index[i]], fTR->PfTau3Pz[fJetTaus.index[i]],fTR->PfTau3E[fJetTaus.index[i]]);
			fHemiObject.type="jet"; fHemiObject.index=i; fHemiObject.hemi=-1;
			hemiobject.objects.push_back(fHemiObject);
		}
	}
	for(int i=0; i< fElecs.size(); ++i){
		TLorentzVector v;
		v.SetPxPyPzE(fTR->PfEl3Px[fElecs[i]], fTR->PfEl3Py[fElecs[i]], fTR->PfEl3Pz[fElecs[i]],fTR->PfEl3E[fElecs[i]]);
		p4s.push_back(v);
		fHemiObject.type="ele"; fHemiObject.index=i; fHemiObject.hemi=-1;
		hemiobject.objects.push_back(fHemiObject);
	}
	for(int i=0; i< fMuons.size(); ++i){
		TLorentzVector v;
		v.SetPxPyPzE(fTR->PfMu3Px[fMuons[i]], fTR->PfMu3Py[fMuons[i]], fTR->PfMu3Pz[fMuons[i]], fTR->PfMu3E[fMuons[i]]);
		p4s.push_back(v);
		fHemiObject.type="muo"; fHemiObject.index=i; fHemiObject.hemi=-1;
		hemiobject.objects.push_back(fHemiObject);
	}
	
	if(p4s.size()<2) {
		return; // protection against events with only one jet
	}

	// get pseudojets and minDHT
	TLorentzVector pj1, pj2;
	hemiobject.minDHT  =MinDeltaHt_pseudojets(p4s, pj1, pj2);
    	
	std::vector<double> pTs; for(unsigned i=0; i<p4s.size(); i++) pTs.push_back(p4s[i].Pt());
    	const double sumPT = accumulate( pTs.begin(), pTs.end(), double(0) );
	const TLorentzVector sumP4 = accumulate( p4s.begin(), p4s.end(), TLorentzVector() );

    	hemiobject.alphaT  = 0.5 * ( sumPT - hemiobject.minDHT ) / sqrt( sumPT*sumPT - sumP4.Perp2() );
		
	TLorentzVector pmiss;
	pmiss.SetPx(MET().Px());
	pmiss.SetPy(MET().Py());
	
	// fill UTM
	hemiobject.UTM          = - pmiss - pj1 - pj2;

	// fill MT2
	hemiobject.MT2          = GetMT2(pj1, 0., pj2, 0., pmiss, 0.);
	
	// fill MCT
	TVector2 pmiss_vector2;
	pmiss_vector2.Set(MET().Px(), MET().Py());
	TLorentzVector downstream(0.,0.,0.,0.); // no particles are downstream, i.e. not selected jets are upstream. 
	hemiobject.MCT   =GetMCTcorr(pj1, pj2, downstream, pmiss_vector2);
	
	// fill pseudojets
	hemiobject.pjet1 =pj1;
	hemiobject.pjet2 =pj2;
	hemiobject.dPhi  =Util::DeltaPhi(pj1.Phi(), pj2.Phi());

}

// ****************************************************************************************************
double MassAnalysis::GetMT(TLorentzVector lv1, double m1, TLorentzVector lv2, double m2){
	// returns Mt. Not not rely on Pz measurement -> also suitable if one lv == MET  
	double ET_1 = sqrt(m1*m1 + lv1.Perp2());
	double ET_2 = sqrt(m2*m2 + lv2.Perp2());
	double MTsquared = m1*m1 + m2*m2 + 2*ET_1*ET_2 - 2*lv1.Px()*lv2.Px() - 2*lv1.Py()*lv2.Py();
	if(MTsquared < 0) return -sqrt(MTsquared);
	else              return  sqrt(MTsquared);
}

double MassAnalysis::GetMT(TLorentzVector lv1, TLorentzVector lv2){
	return GetMT(lv1, 0., lv2, 0.);
}

// ****************************************************************************************************
double MassAnalysis::GetMCTperp(TLorentzVector p1, TLorentzVector p2, TLorentzVector P_UTM){
	double m1=p1.M();
	double m2=p2.M();
	
	TVector3 p1_T_perp  =  GetMomPerp(p1, P_UTM); 
	TVector3 p2_T_perp  =  GetMomPerp(p2, P_UTM);
	
	double E1_T_perp = sqrt( pow(m1,2) + p1_T_perp.Mag2());
	double E2_T_perp = sqrt( pow(m2,2) + p2_T_perp.Mag2());	
	
//	double MCT_perp= sqrt( pow(m1,2) + pow(m2,2) + 2 *(E1_T_perp*E2_T_perp + p1_T_perp.Dot(p2_T_perp)) );
	double MCT_perp= sqrt( pow(m1,2) + pow(m2,2) + 2 *(E1_T_perp*E2_T_perp + p1_T_perp.Mag()*p2_T_perp.Mag()*cos(p1_T_perp.Angle(p2_T_perp))) );


	// compare to Tovey Impelemntation
	TLorentzVector DTM(0.,0.,0.,0.);
	TVector2       pmiss;
	pmiss.Set(-P_UTM.Px()-p1.Px()-p2.Px(), -P_UTM.Py()-p1.Py()-p2.Py());
	double ToveyMCTperp=GetToveyMCTperp(p1, p2, DTM, pmiss);
	// if(fabs(MCT_perp - ToveyMCTperp)/ToveyMCTperp > 0.01) cout << "MCTperp" << MCT_perp << " Tovey " << ToveyMCTperp << endl;


	return MCT_perp;
}

// *****************************************************************************************************
double MassAnalysis::GetToveyMCTperp(TLorentzVector p1, TLorentzVector p2, TLorentzVector DTM, TVector2 pmiss){
	fMCT = new TMctLib;
	double MCTperp = fMCT -> mcy(p1, p2, DTM, pmiss);
	delete fMCT;
	return MCTperp;
}

// *****************************************************************************************************
double MassAnalysis::GetMCTcorr(TLorentzVector p1, TLorentzVector p2, TLorentzVector DTM, TVector2 pmiss){
	fMCT = new TMctLib;
	double MCTcorr = fMCT -> mctcorr(p1, p2, DTM, pmiss, 7000., 0.);
	delete fMCT;
	return MCTcorr;
}

// ******************************************************************************************************
double MassAnalysis::GetMCTcorr(TLorentzVector p1, TLorentzVector p2, TLorentzVector UTM){
	TVector2 pmiss;
	pmiss.Set(-UTM.Px()-p1.Px()-p2.Px(), -UTM.Py()-p1.Py()-p2.Py());
	TLorentzVector DTM(0.,0.,0.,0.);

	fMCT = new TMctLib;
	double MCTcorr = fMCT -> mctcorr(p1, p2, DTM, pmiss, 7000., 0.);
	delete fMCT;
	return MCTcorr;

}
// ******************************************************************************************************
double MassAnalysis::GetMT2perp(TLorentzVector p1, TLorentzVector p2, TLorentzVector P_UTM, double m_inv){
	TVector3 p1_T_perp =  GetMomPerp(p1, P_UTM);
	TVector3 p2_T_perp =  GetMomPerp(p2, P_UTM);
	
	// double A_T_perp   = 1/2.*(p1_T_perp.Mag()*p2_T_perp.Mag()+p1_T_perp.Dot(p2_T_perp) );
	// above implementation is problematic: A_T_perp is sometimes <0 (presicion problem) 
	// workaround: implementation below
	double A_T_perp   = 1/2.*( p1_T_perp.Mag()*p2_T_perp.Mag() * (1+cos(p1_T_perp.Angle(p2_T_perp)) ));
	double MT2_perp   = sqrt(A_T_perp) + sqrt(A_T_perp + pow(m_inv,2));
	
	return MT2_perp;
}

// *******************************************************************************************************
TVector3 MassAnalysis::GetMomPerp(TLorentzVector p, TLorentzVector P_UTM){
	TVector3 p_T(p.Px(), p.Py(), 0.);
	TVector3 pUTM_T(P_UTM.Px(), P_UTM.Py(), 0.);
	double pUTM_T_norm2=pUTM_T.Mag2();
	
	TVector3 p_T_perp = p_T - 1/pUTM_T_norm2*(p_T.Dot(pUTM_T))*pUTM_T ;	
	return p_T_perp;
}

// ****************************************************************************************************
double MassAnalysis::GetMCT(TLorentzVector p1, TLorentzVector p2){
	double m1=p1.M();
	double m2=p2.M();
	double ET1=sqrt( pow(p1.Pt(),2) +pow(m1,2) );
	double ET2=sqrt( pow(p2.Pt(),2) +pow(m2,2) );
	
	double MCTsquared=pow(m1,2) + pow(m2,2) +2*(ET1*ET2+ (p1.Px() * p2.Px()) + (p1.Py() * p2.Py()));
	
	// Tovey implementation
	fMCT = new TMctLib();
	double Tovey_MCT = fMCT -> mct(p1, p2);
	delete fMCT;
//	if(fabs(Tovey_MCT -sqrt(MCTsquared))/Tovey_MCT > 0.01) cout << "Tovey_MCT "<< Tovey_MCT << " MCT " << sqrt(MCTsquared) << endl;

	return sqrt(MCTsquared);
}

// ****************************************************************************************************
double MassAnalysis::GetAlphaT(std::vector<TLorentzVector>& p4s) {
    if(p4s.size()<2) return 0;
    
    std::vector<double> pTs; for(unsigned i=0; i<p4s.size(); i++) pTs.push_back(p4s[i].Pt());
    const std::vector<double> DHTper( DeltaSumPt_permutations(p4s) );
    
    const double mDHT = *(std::min_element( DHTper.begin(), DHTper.end() ));
    const double sumPT = accumulate( pTs.begin(), pTs.end(), double(0) );
    const TLorentzVector sumP4 = accumulate( p4s.begin(), p4s.end(), TLorentzVector() );

    return 0.5 * ( sumPT - mDHT ) / sqrt( sumPT*sumPT - sumP4.Perp2() );
}

std::vector<double> MassAnalysis::DeltaSumPt_permutations(std::vector<TLorentzVector>& p4s) {
  	std::vector<std::vector<double> > ht( 1<<(p4s.size()-1) , std::vector<double>( 2, 0.) ); // initializiert einen vector ht der size 1<<(p4s.size()-1), wobei jeder eintrag ein vector (0,0) ist
  	for(unsigned i=0; i < ht.size(); i++) {                                                  // 1<<(p4s.size()-1) = 2^(p4s.size() -1)
    	for(unsigned j=0; j < p4s.size(); j++) {
			ht [i] [(i/(1<<j))%2] += p4s[j].Pt();
    	}
  	}
  	std::vector<double> deltaHT; for(unsigned i=0; i<ht.size(); i++) deltaHT.push_back(fabs(ht[i][0]-ht[i][1]));
  	return deltaHT;
}

// ***************************************************************************************************
double MassAnalysis::MinDeltaHt_pseudojets(std::vector<TLorentzVector>& p4s, TLorentzVector& pj1, TLorentzVector& pj2){

  	std::vector<std::vector<double> > ht( 1<<(p4s.size()-1) , std::vector<double>( 2, 0.) ); // initializiert einen vector ht der size 1<<(p4s.size()-1), wobei jeder eintrag ein vector (0,0) ist
												 // 1<<(p4s.size()-1) = 2^(p4s.size() -1)
       	
  	std::vector<std::vector<double> > px( 1<<(p4s.size()-1) , std::vector<double>( 2, 0.) ); 
  	std::vector<std::vector<double> > py( 1<<(p4s.size()-1) , std::vector<double>( 2, 0.) ); 
  	std::vector<std::vector<double> > pz( 1<<(p4s.size()-1) , std::vector<double>( 2, 0.) ); 
  	std::vector<std::vector<double> > E ( 1<<(p4s.size()-1) , std::vector<double>( 2, 0.) ); 

  	for(unsigned i=0; i < ht.size(); i++) {                                             
    	for(unsigned j=0; j < p4s.size(); j++) {
			ht [i] [(i/(1<<j))%2] += p4s[j].Pt();
			px [i] [(i/(1<<j))%2] += p4s[j].Px();
			py [i] [(i/(1<<j))%2] += p4s[j].Py();
			pz [i] [(i/(1<<j))%2] += p4s[j].Pz();
			E  [i] [(i/(1<<j))%2] += p4s[j].E();
    	}
  	}
  	std::vector<double> deltaHT; for(unsigned i=0; i<ht.size(); i++) deltaHT.push_back(fabs(ht[i][0]-ht[i][1]));
	const double mDHT = *(std::min_element( deltaHT.begin(), deltaHT.end() ));
	int pos=distance(deltaHT.begin(), min_element(deltaHT.begin(), deltaHT.end()));
	pj1.SetPxPyPzE(px[pos][0], py[pos][0], pz[pos][0], E[pos][0]);
	pj2.SetPxPyPzE(px[pos][1], py[pos][1], pz[pos][1], E[pos][1]);
	return mDHT;
}


// ****************************************************************************************************
double MassAnalysis::GetMT2(TLorentzVector v1, double mv1, TLorentzVector v2, double mv2, TLorentzVector p_unobs, int m_invisible){
	double pa[3];
	double pb[3];
	double pmiss[3];

	pmiss[1]=p_unobs.Px();
	pmiss[2]=p_unobs.Py();

	pa[0]=mv1;
	pa[1]=v1.Px();
	pa[2]=v1.Py();

	pb[0]=mv2;
	pb[1]=v2.Px();
	pb[2]=v2.Py();
	
	fMT2 = new Davismt2();
	fMT2->set_momenta(pa, pb, pmiss);
	fMT2->set_mn(m_invisible);
	double MT2=fMT2->get_mt2();
	delete fMT2;
	return MT2;
}


// ****************************************************************************************************
vector<TLorentzVector> MassAnalysis::GetLepton4Momenta(){
	vector<TLorentzVector> momenta;
	for(int i=0; i<fElecs.size(); ++i){
		TLorentzVector p;
		p.SetPtEtaPhiE(fTR->PfEl3Pt[fElecs[i]],fTR->PfEl3Eta[fElecs[i]], fTR->PfEl3Phi[fElecs[i]], fTR->PfEl3E[fElecs[i]]);
		momenta.push_back(p);
	}
	for(int i=0; i<fMuons.size(); ++i){
		TLorentzVector p;
		p.SetPtEtaPhiM(fTR->PfMu3Pt[fMuons[i]],fTR->PfMu3Eta[fMuons[i]], fTR->PfMu3Phi[fMuons[i]],0.105);
		momenta.push_back(p);
	}
	return momenta;
}

// ********************************************************************************
// BE file parser
void MassAnalysis::DeadCellParser(DeadCellFilter &DeadCellFilter_, string file_){
	string line;	
	string path="/shome/pnef/SUSY/SUSY_macros/LeptJetMult/ECALDeadCell/";
	string file=path+file_;
	ifstream IN(file.c_str());
	if (!IN.is_open()) {cout << "ERROR: cannot open dead cell file " << file << endl; exit(1);}
	else{
		if(fVerbose>0) cout << "--------------------------"          << endl;
		if(fVerbose>0) cout << "DeadCellParser: read file " << file  << endl;
		while ( ! IN.eof() ){
			getline (IN, line);
			TString Line=  line.c_str();
			if(Line.EndsWith("\",") && Line.BeginsWith("\"")){
				Line.ReplaceAll("\"", "");
				Line.ReplaceAll(",", "");
				TObjArray *p= (TObjArray*) Line.Tokenize(":");
				TString run  =((TObjString*)p->At(0))->GetString();
				TString lumi =((TObjString*)p->At(1))->GetString();
				TString event=((TObjString*)p->At(2))->GetString();
				DeadCellFilter_.run.  push_back(run.Atoi());
				DeadCellFilter_.lumi. push_back(lumi.Atoi());
				DeadCellFilter_.event.push_back(event.Atoi());
			}
		}
		if(fVerbose >0) cout << "DeadCellParser: read " <<  DeadCellFilter_.run.size() << " events to be filtered. " << endl; 
		if(fVerbose >0) cout << "--------------------------"          << endl;
	}
}
// *************************************

// ****************************************************************************************************
void MassAnalysis::End(){
	cout << " *************************************************************** " << endl;
	cout << " MassAnalysis::End()                                             " << endl;
	cout << " *************************************************************** " << endl;
	
	fHistFile->cd();	

	// write tree
	fATree->Write();
	fHistFile                ->Close();

	cout << " MassAnalysis::RealEnd()                                             " << endl;

}

