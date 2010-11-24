#include "RatioAnalysis.hh"
#include "base/TreeReader.hh"
#include "helper/Utilities.hh"
#include <iostream>
#include <sstream>

using namespace std;

RatioAnalysis::RatioAnalysis(TreeReader *tr) : MultiplicityAnalysisBase(tr){
	Util::SetStyle();	
}

RatioAnalysis::~RatioAnalysis(){
}

void RatioAnalysis::Begin(){
	const char* filename = "Ratio_histos.root";
	fMPHistFile = new TFile(fOutputDir + TString(filename), "RECREATE");

	fRatHist  = new TH1D("RatioStats", "Histo of Ratio statistics", 63, 0, 63);	
	fRatHist->SetStats(false);

	fDRJ12  = new TH1D("DRJ12", "Delta R between 2 leading jets; DR(j1,j2)", 100, 0., 5.);
	fPtJets = new TH1D("PtJets", "pT of all leading jets", 100, 0., 250.);
	fPtJ1  = new TH1D("PtJ1", "Pt of leading jet 1 for DR<1", 100, 0., 250.);
	fPtJ2  = new TH1D("PtJ2", "Pt of leading jet 2 for DR<1", 100, 0., 250.);
	fEtaJets = new TH1D("EtaJets", "Eta of 2 leading jets", 50, -3., 3.);
	fRptDPhi2j = new TH2D("RptDPhi2j", "dPhi(j1,j2) versus Ptj2/Ptj1 leading jets; Ptj2/Ptj1; dPhi(j1,j2)", 
			      50, 0., 1., 50, 0., 2.);
	fAlpT2j  = new TH1D("AlpT2j", "AlphaT in 2 jet events; AlphaT(j1,j2)", 100, 0., 3.);
	fBProbJets = new TH1D("BProbJets", "b-tagging prob for all jets", 50, 0., 5.);

	fEta2B = new TH1D("Eta2B", "Eta of bs in 2b events", 50, -3., 3.);
	fdPhi2B = new TH1D("dPhi2B", "Delta Phi of 2 b in 2b events", 50, 0., 3.5);
	fDR2B = new TH1D("DR2B", "dR(b1,b2) ; dR(b1,b2)", 50, 0., 6.);
	//	fEta2Bcand = new TH1D("Eta2Bcand", "Eta of lost b in 1b events", 50, -5., 5.);
	//	fdPhi2Bcand = new TH1D("dPhi2Bcand", "Delta Phi(b,bcand) in 1b events", 50, 0., 3.5);
	fBProb2B = new TH1D("BProb2B", "b-tagging prob for 2b events", 50, 0., 5.);
	//	fBProb2Bcand = new TH1D("BProb2Bcand", "b-tagging prob for b-jet candidate in 1b events", 50, 0., 5.);
	fNjets2B = new TH1D("Njets2B", "Number of jets in 2b events", 10, 0., 10.);
	fMass2B = new TH1D("Mass2B", "Mass of b-jets in 2b events", 50, 0., 50.);
	fNjets2Bnear = new TH1D("Njets2Bnear", "Number of jets in 2b events with dPhi(2b) < 1.", 10, 0., 10.);
	fMET2Bnear = new TH1D("MET2Bnear", "MET in 2b events with dPhi(2b) < 1.", 100, 0., 100.);
	fdPhiJBnear = new TH2D("dPhiJBnear", "dPhi(b, nearest jet) in 2b events with dPhi(2b) < 1.;dPhi(b1,b2);dPhi(b,j)",
			       50, 0., 3.5, 50, 0., 3.5);
	fEta1B = new TH1D("Eta1B", "Eta of b in 1b events", 50, -3., 3.);
	fBProb1B = new TH1D("BProb1B", "b-tagging prob for b in 1b events", 50, 0., 5.);
	fNjets1B = new TH1D("Njets1B", "Number of jets in 1b events", 10, 0., 10.);
	//	fNjets2Bcand = new TH1D("Njets2Bcand", "Number of jets in 1b events with b-candidate", 10, 0., 10.);
	fMass1B = new TH1D("Mass1B", "Mass of b-jet in 1b events", 50, 0., 50.);
	//	fMass2Bcand = new TH1D("Mass2Bcand", "Mass of b-jet candidate in 1b events", 50, 0., 250.);
	fMETdPhi1B = new TH2D("METdPhi1B", "MET vs dPhi(MET,b) in 1b events", 100, 0., 100., 50, 0., 3.5);
	fdPhiMET1B = new TH1D("dPhiMET1B", "dPhi(MET,b) in 1b events for MET > 5.",50, 0., 3.5);

	fJMult = new TH1D("JMult", "Jet multiplicity in events with >= 2 jets", 10, 0., 10.);
	fJMult1B = new TH1D("JMult1B", "Jet multiplicity in 1b events", 10, 0., 10.);
	fJMult2B = new TH1D("JMult2B", "Jet multiplicity in 2b events", 10, 0., 10.);
	fJMult3B = new TH1D("JMult3B", "Jet multiplicity in 3b events", 10, 0., 10.);
	fJMult4B = new TH1D("JMult4B", "Jet multiplicity in 4b events", 10, 0., 10.);


	// set the efficiencies and fake rates
	fEffe = 1.;
	fEffm = 1.;
	fEffb = 0.5;
	fdEffb = 0.05;
	fFake = 0.;
	fFakm = 0.;
	fFakb = 0.005;
	
	// initialize the counters
	counter=0;
	fNEvtLeptons = 0;
	// jet ratios
	fNEvtJets = 0;
	fNEvt1Jets = 0;
	fNEvt2Jets = 0;
	fNEvt3Jets = 0;
	fNEvt4Jets = 0;
	fNEvt5Jets = 0;
	fNEvt6Jets = 0;
	fNEvt7Jets = 0;
	fNEvt8Jets = 0;
	fNEvtBJets = 0;
	fNEvt1BJets = 0;
	fNEvt2BJets = 0;
	fNEvt3BJets = 0;
	fNEvt4BJets = 0;
	fNEvtj12near1 = 0;
	fNEvtj12near2 = 0;
	fNEvtj12near4 = 0;
	// jet multiplicities in 1l events
	fNEvt1J1l = 0;
	fNEvt2J1l = 0;
	fNEvt3J1l = 0;
	fNEvt4J1l = 0;
	fNEvt5J1l = 0;
	fNEvt6J1l = 0;
	fNEvt7J1l = 0;
	fNEvt8J1l = 0;
	// b-jet multiplicities in 1l events
	fNEvt1B1l = 0;
	fNEvt2B1l = 0;
	fNEvt3B1l = 0;
	fNEvt4B1l = 0;

	// lepton ratios
	fNEvtepep = 0;
	fNEvtepen = 0;
	fNEvtenen = 0;	
	fNEvtmpmp = 0;
	fNEvtmpmn = 0;
	fNEvtmnmn = 0;	
	fNEvtepmp = 0;
	fNEvtepmn = 0;
	fNEvtenmp = 0;
	fNEvtenmn = 0;	
	// lepton ratios, events with b-jets
	fNEvtepepB = 0;
	fNEvtepenB = 0;
	fNEvtenenB = 0;	
	fNEvtmpmpB = 0;
	fNEvtmpmnB = 0;
	fNEvtmnmnB = 0;	
	fNEvtepmpB = 0;
	fNEvtepmnB = 0;
	fNEvtenmpB = 0;
	fNEvtenmnB = 0;	
	// jet multiplicities in 2l events
	fNEvt1BOS2l = 0;
	fNEvt2BOS2l = 0;
	fNEvt3BOS2l = 0;
	fNEvt4BOS2l = 0;
	fNEvt1BSS2l = 0;
	fNEvt2BSS2l = 0;
	fNEvt3BSS2l = 0;
	fNEvt4BSS2l = 0;
	fNEvt1BOSem = 0;
	fNEvt2BOSem = 0;
	fNEvt3BOSem = 0;
	fNEvt4BOSem = 0;

	fNevtAnom = 0;

	fPtSum1B = 0.;
	fPtSumsq1B = 0.;
	fPtSum2B = 0.;
	fPtSumsq2B = 0.;

	
	
}

void RatioAnalysis::Analyze(){
	//define an array for leptons (es+mus) with
	//information on pt, eta and category
	//category 1: e+
	//category 2: e-
	//category 3: mu+
	//category 4: mu-
	
	fLeptCat.clear();
  
	// ---------------------------------------------------
	// Initialize fElecs, fJetsLoose, fBJets, fMuons, fLeptConfig 
	InitializeEvent();
	// ----------------------------------------------------
		
	// --------------------------------------------------------------------
	// check if event passes selection cuts and trigger requirements
	if(! IsSelectedEvent()){return;}
	
	
	// -----------------------------------------------------------
	// do counting 
		
	const int NLepts = 40;
	int LeptCat[NLepts];
	double LeptPt[NLepts];
	double LeptEta[NLepts];

	for(int tmp=0; tmp<20; ++tmp){
		LeptCat[tmp]=0;
		LeptPt[tmp]=-999.;
		LeptEta[tmp]=-999.;
	}	
	
	//loop over es and mus and fill the leptons array
	int leptcounter=-1;
	for(int i=0; i<fElecs.size(); ++i){
		leptcounter++;
		LeptPt[leptcounter] =fTR->ElPt[fElecs[i]];
		LeptEta[leptcounter]=fTR->ElEta[fElecs[i]];
		if(fTR->ElCharge[fElecs[i]]>0) {LeptCat[leptcounter]=1;}
		else {LeptCat[leptcounter]=2;}
	}
	for(int i=0; i<fMuons.size(); ++i){
		leptcounter++;
		LeptPt[leptcounter] =fTR->MuPt[fMuons[i]];
		LeptEta[leptcounter]=fTR->MuEta[fMuons[i]];
		if(fTR->MuCharge[fMuons[i]]>0) {LeptCat[leptcounter]=3;}
		else {LeptCat[leptcounter]=4;}
	}

	//total number of leptons in the event
	int NLeptons= fElecs.size()+fMuons.size();
	
	// reorder here ...
	bool changed = false;
	do {
		changed = false;
		for (int i = 0; i < NLeptons; i++){
			for (int j = i+1; j < NLeptons; j++){
				if (LeptPt[i] < LeptPt[j] ){
					int cat=LeptCat[i];
					double pt=LeptPt[i];
					double eta=LeptEta[i];
					LeptCat[i]=LeptCat[j];
					LeptPt[i]=LeptPt[j];
					LeptEta[i]=LeptEta[j];
					LeptCat[j]=cat;
					LeptPt[j]=pt;
					LeptEta[j]=eta;
					changed = true;
				}
			}
		}
	} while (changed);

	// save the lepton categories
	for (int i = 0; i < NLeptons; ++i) {
	  fLeptCat.push_back(LeptCat[i]);
	}




	// count number of jets
	
	fNQJets = fJetsLoose.size();
	fNBJets = fBJets.size();
	
	
	double ptJet1 = 0., ptJet2 = 0.;
	int ij1 = -1, ij2 = -1;
	int ib1 = -1;
	int ib2 = -1;
	
	for(int ij=0; ij<fJetsLoose.size(); ++ij){
		if(fTR->PFJPt[fJetsLoose[ij]] > ptJet1){
			ptJet2 = ptJet1;
		  	ij2 = ij1;
			ptJet1 = fTR->PFJPt[fJetsLoose[ij]];
			ij1 = fJetsLoose[ij];
		} else if( fTR->PFJPt[fJetsLoose[ij]] > ptJet2 ){
			ptJet2 = fTR->PFJPt[fJetsLoose[ij]];
			ij2 = fJetsLoose[ij];
		}
	}
	for(int ij=0; ij<fBJets.size(); ++ij){
		if (ib1 < 0) ib1 = fBJets[ij];
		else if (ib2 < 0) ib2 = fBJets[ij];
	}
	

	double deltaRj12 = 999.;
	if (ij2 >= 0) {
	  deltaRj12 = Util::GetDeltaR(fTR->PFJEta[ij1], fTR->PFJEta[ij2], fTR->JPhi[ij1], fTR->JPhi[ij2]);
	  double dPhi2j = Util::DeltaPhi(fTR->JPhi[ij1], fTR->JPhi[ij2]);
	  fDRJ12->Fill(deltaRj12);
	  fPtJets ->Fill(fTR->PFJPt[ij1]);
	  fPtJets ->Fill(fTR->PFJPt[ij2]);
	  fRptDPhi2j->Fill(fTR->PFJPt[ij2]/fTR->PFJPt[ij1], dPhi2j);

	  if( fJetsLoose.size()==2){
//	  if (fTR->NJets == 2) { 
	    double mT = sqrt((fTR->PFJPt[ij1]+fTR->PFJPt[ij2])*(fTR->PFJPt[ij1]+fTR->PFJPt[ij2])
	      - (fTR->PFJPx[ij1]+fTR->PFJPx[ij2])*(fTR->PFJPx[ij1]+fTR->PFJPx[ij2])
	      - (fTR->PFJPy[ij1]+fTR->PFJPy[ij2])*(fTR->PFJPy[ij1]+fTR->PFJPy[ij2]) );
	    double alpT = fTR->PFJPt[ij2] / mT;
	    fAlpT2j->Fill(alpT);
	  }
	}
	if (deltaRj12 < 1.) {
	  fPtJ1->Fill(fTR->PFJPt[ij1]);
	  fPtJ2->Fill(fTR->PFJPt[ij2]);
	  fEtaJets ->Fill(fTR->PFJEta[ij1]);
	  fEtaJets ->Fill(fTR->PFJEta[ij2]);
	}

	for(int ij=0; ij < fJetsLoose.size(); ++ij){
	  fBProbJets->Fill(fTR->JbTagProbSimpSVHighEff[fJetsLoose[ij]]);
	}

	if (ib1 >= 0 && ib2 >= 0) {
	  fPtSum2B += fTR->PFJPt[ib1];
	  fPtSum2B += fTR->PFJPt[ib2];
	  fPtSumsq2B += fTR->PFJPt[ib1] *  fTR->PFJPt[ib1];
	  fPtSumsq2B += fTR->PFJPt[ib2] *  fTR->PFJPt[ib2];
	  fEta2B->Fill(fTR->PFJEta[ib1]);
	  fEta2B->Fill(fTR->PFJEta[ib2]);
	  double dPhi2b = Util::DeltaPhi(fTR->JPhi[ib1], fTR->JPhi[ib2]);
	  fdPhi2B->Fill(dPhi2b);
	  double dR2B = Util::GetDeltaR(fTR->PFJEta[ib1], fTR->PFJEta[ib2], fTR->JPhi[ib1], fTR->JPhi[ib2]);
	  fDR2B->Fill(dR2B);
	  fBProb2B->Fill(fTR->JbTagProbSimpSVHighEff[ib1]);
	  fBProb2B->Fill(fTR->JbTagProbSimpSVHighEff[ib2]);
	  fNjets2B->Fill(fTR->NJets);
	  double mass1 = sqrt(fTR->PFJE[ib1]*fTR->PFJE[ib1] - fTR->PFJPx[ib1]*fTR->PFJPx[ib1]
			     - fTR->PFJPy[ib1]*fTR->PFJPy[ib1] - fTR->PFJPz[ib1]*fTR->PFJPz[ib1]);
	  fMass2B->Fill(mass1);
	  double mass2 = sqrt(fTR->PFJE[ib2]*fTR->PFJE[ib2] - fTR->PFJPx[ib2]*fTR->PFJPx[ib2]
			     - fTR->PFJPy[ib2]*fTR->PFJPy[ib2] - fTR->PFJPz[ib2]*fTR->PFJPz[ib2]);
	  fMass2B->Fill(mass2);
	  //	  fMass2B->Fill(fTR->JMass[ib1]);
	  //	  fMass2B->Fill(fTR->JMass[ib2]);
	  if (dPhi2b < 1.) {
	    fNjets2Bnear->Fill(fTR->NJets);
	    fMET2Bnear->Fill(fTR->PFMET);
	    double dPhibjmin = 999.;
	    for(int ij=0; ij < fTR->NJets; ++ij){
	      if (ij == ib1 || ij == ib2) continue;
	      double dPhib1J = Util::DeltaPhi(fTR->JPhi[ib1], fTR->JPhi[ij]);
	      double dPhib2J = Util::DeltaPhi(fTR->JPhi[ib2], fTR->JPhi[ij]);
	      if (dPhib1J < dPhibjmin) dPhibjmin = dPhib1J;
	      if (dPhib2J < dPhibjmin) dPhibjmin = dPhib2J;
	    }
	    fdPhiJBnear->Fill(dPhi2b, dPhibjmin);
	  }
	}
	if (ib1 >= 0 && ib2 < 0) {
	  fPtSum1B += fTR->PFJPt[ib1];
	  fPtSumsq1B += fTR->PFJPt[ib1] *  fTR->PFJPt[ib1];
	  fEta1B->Fill(fTR->PFJEta[ib1]);
	  fBProb1B->Fill(fTR->JbTagProbSimpSVHighEff[ib1]);
	  fNjets1B->Fill(fTR->NJets);
	  double mass1 = sqrt(fTR->PFJE[ib1]*fTR->PFJE[ib1] - fTR->PFJPx[ib1]*fTR->PFJPx[ib1]
			     - fTR->PFJPy[ib1]*fTR->PFJPy[ib1] - fTR->PFJPz[ib1]*fTR->PFJPz[ib1]);
	  fMass1B->Fill(mass1);
	  //	  fMass1B->Fill(fTR->JMass[ib1]);
	  double dPhiMET1B = Util::DeltaPhi(fTR->JPhi[ib1], fTR->PFMETphi);
	  fMETdPhi1B->Fill(fTR->PFMET, dPhiMET1B);
	  if (fTR->PFMET > 5.) {
	    fdPhiMET1B->Fill(dPhiMET1B);
	  }
	}
	fJMult->Fill(fTR->NJets);
	if (fNBJets == 1) fJMult1B->Fill(fTR->NJets);
	else if (fNBJets == 2) fJMult2B->Fill(fTR->NJets);
	else if (fNBJets == 3) fJMult3B->Fill(fTR->NJets);
	else if (fNBJets == 4) fJMult4B->Fill(fTR->NJets);



	// save the interesting numbers of events
	// for jet ratios
	if (fNQJets > 0) {
	  fNEvtJets++;
	  if (NLeptons == 1) fNEvtJ1l++;
	} 
	if (fNQJets == 1) {
	  fNEvt1Jets++;
	  if (NLeptons == 1) fNEvt1J1l++;
	} else if (fNQJets == 2) {
	  fNEvt2Jets++;
	  if (NLeptons == 1) fNEvt2J1l++;
	} else if (fNQJets == 3) {
	  fNEvt3Jets++;
	  if (NLeptons == 1) fNEvt3J1l++;
	} else if (fNQJets == 4) {
	  fNEvt4Jets++;
	  if (NLeptons == 1) fNEvt4J1l++;
	} else if (fNQJets == 5) {
	  fNEvt5Jets++;
	  if (NLeptons == 1) fNEvt5J1l++;
	} else if (fNQJets == 6) {
	  fNEvt6Jets++;
	  if (NLeptons == 1) fNEvt6J1l++;
	} else if (fNQJets == 7) {
	  fNEvt7Jets++;
	  if (NLeptons == 1) fNEvt7J1l++;
	} else if (fNQJets >= 8) {
	  fNEvt8Jets++;
	  if (NLeptons == 1) fNEvt8J1l++;
        }
	if (fNBJets > 0) {
	  fNEvtBJets++;
	  if (NLeptons == 1) fNEvtB1l++;
	}
	if (fNBJets == 1) {
	  fNEvt1BJets++;
	  if (NLeptons == 1) fNEvt1B1l++;
        }
	else if (fNBJets == 2) {
	  fNEvt2BJets++;
	  if (NLeptons == 1) fNEvt2B1l++;
	}
	else if (fNBJets == 3) {
	  fNEvt3BJets++;
	  if (NLeptons == 1) fNEvt3B1l++;
	}
	else if (fNBJets == 4) {
	  fNEvt4BJets++;
	  if (NLeptons == 1) fNEvt4B1l++;
	}

	int isAnomJ = 0;
	if (deltaRj12 < 1.) {
	  fNEvtj12near1++;
	  isAnomJ = 1;
	}
	if (deltaRj12 < 2.) fNEvtj12near2++;
	if (deltaRj12 < 4.) fNEvtj12near4++;
	if (isAnomJ != 0) {
	  fNevtAnom++;
	  fIRun.push_back(fTR->Run);
	  fILumi.push_back(fTR->LumiSection);
	  fINber.push_back(fTR->Event); 
	  fINlep.push_back(-1);
	}
	if (fNBJets > 2) {
	  fNevtAnom++;
	  fIRun.push_back(fTR->Run);
	  fILumi.push_back(fTR->LumiSection);
	  fINber.push_back(fTR->Event); 
	  if (fNBJets == 3) fINlep.push_back(-2);
	  if (fNBJets == 4) fINlep.push_back(-3);
	  // cout << " Event with " << fNBJets << " b-jets, flag = " << fINlep[fINlep.size()-1] << endl;
	  // cout << " Run/Lum/Evt = " << fTR->Run << " " << fTR->LumiSection << " "  << fTR->Event << endl;
	}

	// for charge and lepton ratios
	SaveLeptConfs();
	

}

void RatioAnalysis::SaveLeptConfs() {
  // saves the lepton configurations

  	if (fLeptCat.size() < 2) return;
  	
  	// for all dilepton events
  	int isSS2l = 0;
  	int isOS2l = 0;
  	int isOSem = 0;
  	if (fLeptCat.size() == 2) {
  	  	if (fLeptCat[0] == 1 && fLeptCat[1] == 1) fNEvtepep++;
  	  	else if (fLeptCat[0] == 2 && fLeptCat[1] == 2) fNEvtenen++;
  	  	else if (fLeptCat[0] == 1 && fLeptCat[1] == 2) fNEvtepen++;
  	  	else if (fLeptCat[0] == 2 && fLeptCat[1] == 1) fNEvtepen++;
  	  	
  	  	else if (fLeptCat[0] == 3 && fLeptCat[1] == 3) fNEvtmpmp++;
  	  	else if (fLeptCat[0] == 4 && fLeptCat[1] == 4) fNEvtmnmn++;
  	  	else if (fLeptCat[0] == 3 && fLeptCat[1] == 4) fNEvtmpmn++;
  	  	else if (fLeptCat[0] == 4 && fLeptCat[1] == 3) fNEvtmpmn++;
  	  	
  	  	else if (fLeptCat[0] == 1 && fLeptCat[1] == 3) fNEvtepmp++;
  	  	else if (fLeptCat[0] == 3 && fLeptCat[1] == 1) fNEvtepmp++;
  	  	else if (fLeptCat[0] == 2 && fLeptCat[1] == 4) fNEvtenmn++;
  	  	else if (fLeptCat[0] == 4 && fLeptCat[1] == 2) fNEvtenmn++;
  	  	else if (fLeptCat[0] == 1 && fLeptCat[1] == 4) fNEvtepmn++;
  	  	else if (fLeptCat[0] == 4 && fLeptCat[1] == 1) fNEvtepmn++;
  	  	else if (fLeptCat[0] == 2 && fLeptCat[1] == 3) fNEvtenmp++;
  	  	else if (fLeptCat[0] == 3 && fLeptCat[1] == 2) fNEvtenmp++;
  	  	
  	  	// for events with b-jets
  	  	if (fNBJets >= 1.) {
  	  	  	if (fLeptCat[0] == 1 && fLeptCat[1] == 1) {
  	  			fNEvtepepB++;
  	  			isSS2l++;
  	  	  	} else if (fLeptCat[0] == 2 && fLeptCat[1] == 2) {
  	  			fNEvtenenB++;
  	  			isSS2l++;
  	  	  	}else if (fLeptCat[0] == 1 && fLeptCat[1] == 2) {
  	  			fNEvtepenB++;
  	  			isOS2l++;
  	  	  	}else if (fLeptCat[0] == 2 && fLeptCat[1] == 1) {
  	  			fNEvtepenB++;
  	  			isOS2l++;
  	  	  	}
  	  	  
  	  	  	else if (fLeptCat[0] == 3 && fLeptCat[1] == 3) {
  	  			fNEvtmpmpB++;
  	  			isSS2l++;
  	  	  	} else if (fLeptCat[0] == 4 && fLeptCat[1] == 4) {
  	  			fNEvtmnmnB++;
  	  			isSS2l++;
  	  	  	} else if (fLeptCat[0] == 3 && fLeptCat[1] == 4) {
  	  			fNEvtmpmnB++;
  	  			isOS2l++;
  	  	  	} else if (fLeptCat[0] == 4 && fLeptCat[1] == 3) {
    			fNEvtmpmnB++;
  	  			isOS2l++;
  	  	  	}
  	  	
  	  	  	else if (fLeptCat[0] == 1 && fLeptCat[1] == 3) {
  	  			fNEvtepmpB++;
  	  			isSS2l++;
  	  	  	} else if (fLeptCat[0] == 3 && fLeptCat[1] == 1) {
  	  			fNEvtepmpB++;
  	  			isSS2l++;
  	  	  	} else if (fLeptCat[0] == 2 && fLeptCat[1] == 4) {
  	  			fNEvtenmnB++;
  	  			isSS2l++;
  	  	  	} else if (fLeptCat[0] == 4 && fLeptCat[1] == 2) {
  	  			fNEvtenmnB++;
  	  			isSS2l++;
  	  	  	} else if (fLeptCat[0] == 1 && fLeptCat[1] == 4) {
  	  			fNEvtepmnB++;
  	  			isOSem++;
  	  	  	} else if (fLeptCat[0] == 4 && fLeptCat[1] == 1) {
  	  			fNEvtepmnB++;
  	  			isOSem++;
  	  	  	} else if (fLeptCat[0] == 2 && fLeptCat[1] == 3) {
  	  			fNEvtenmpB++;
  	  			isOSem++;
  	  	  	}else if (fLeptCat[0] == 3 && fLeptCat[1] == 2) { 
  	  			fNEvtenmpB++;
  	  			isOSem++;
  	  	  	}
  	  	
  	  	  	if (fNBJets == 1) {
  	  			if (isSS2l != 0) fNEvt1BSS2l++;
  	  			if (isOS2l != 0) fNEvt1BOS2l++;
  	  			if (isOSem != 0) fNEvt1BOSem++;
  	  	  	} else if (fNBJets == 2) {
  	  			if (isSS2l != 0) fNEvt2BSS2l++;
  	  			if (isOS2l != 0) fNEvt2BOS2l++;
				if (isOSem != 0) fNEvt2BOSem++;
  	  	  	} else if (fNBJets == 3) {
  	  			if (isSS2l != 0) fNEvt3BSS2l++;
  	  			if (isOS2l != 0) fNEvt3BOS2l++;
  	  			if (isOSem != 0) fNEvt3BOSem++;
  	 	  	} else if (fNBJets == 4) {
  	  			if (isSS2l != 0) fNEvt4BSS2l++;
  	  			if (isOS2l != 0) fNEvt4BOS2l++;
				if (isOSem != 0) fNEvt4BOSem++;
  	  	  	}
  	  	}
  	}
  	
  	// for events with >2 leptons
  	if (isSS2l != 0 || fLeptCat.size() > 2) {
   		fNevtAnom++;
		fIRun.push_back(fTR->Run);
    		fILumi.push_back(fTR->LumiSection);
		fINber.push_back(fTR->Event); 
  	 	fINlep.push_back(fLeptCat.size());
  	}	
}


void RatioAnalysis::End(){

  cout << " Total triggers processed = " << counter << endl;
  
  PrintJetRatios();
  PrintChaRatios();
  PrintDilRatios();
  PrintAnomEvts();

   SaveRatioHist();

}

void RatioAnalysis::PrintJetRatios() {

cout << endl;
  cout << " Jet ratios " << endl;
  cout << " ========== " << endl;
  cout << endl;
  cout << " Total number of events with jets = " << fNEvtJets << endl;
  if (fNEvtJets <= 0) return;

  // statistics for any jets
  float rat1jbyj = (float)fNEvt1Jets/ (float)fNEvtJets;
  cout << " Number of events with 1 jets = " << fNEvt1Jets << " Ratio to jets = " << rat1jbyj << endl;
  float rat2jbyj = (float)fNEvt2Jets/ (float)fNEvtJets;
  cout << " Number of events with 2 jets = " << fNEvt2Jets << " Ratio to jets = " << rat2jbyj << endl;
  float rat3jbyj = (float)fNEvt3Jets/ (float)fNEvtJets;
  cout << " Number of events with 3 jets = " << fNEvt3Jets << " Ratio to jets = " << rat3jbyj << endl;
  float rat4jbyj = (float)fNEvt4Jets/ (float)fNEvtJets;
  cout << " Number of events with 4 jets = " << fNEvt4Jets << " Ratio to jets = " << rat4jbyj << endl;
  float rat5jbyj = (float)fNEvt5Jets/ (float)fNEvtJets;
  cout << " Number of events with 5 jets = " << fNEvt5Jets << " Ratio to jets = " << rat5jbyj << endl;
  float rat6jbyj = (float)fNEvt6Jets/ (float)fNEvtJets;
  cout << " Number of events with 6 jets = " << fNEvt6Jets << " Ratio to jets = " << rat6jbyj << endl;
  float rat7jbyj = (float)fNEvt7Jets/ (float)fNEvtJets;
  cout << " Number of events with 7 jets = " << fNEvt7Jets << " Ratio to jets = " << rat7jbyj << endl;
  float rat8jbyj = (float)fNEvt8Jets / (float)fNEvtJets;
  cout << " Number of events with >=8 jets = " << fNEvt8Jets << " Ratio to jets = " << rat8jbyj << endl;

  // statistics on acoplanar leading jets
  double ratj12near1 = (float)fNEvtj12near1 / (float)fNEvtJets;
  double ratj12near2 = (float)fNEvtj12near2 / (float)fNEvtJets;
  double ratj12near4 = (float)fNEvtj12near4 / (float)fNEvtJets;
  cout << endl;
  cout << " Number of events with acoplanar leading jets " << endl;
  cout << "  Number of events with DeltaR < 1 = " << fNEvtj12near1 << " Ratio to all jets = " << ratj12near1 << endl;
  cout << "  Number of events with DeltaR < 2 = " << fNEvtj12near2 << " Ratio to all jets = " << ratj12near2 << endl;
  cout << "  Number of events with DeltaR < 4 = " << fNEvtj12near4 << " Ratio to all jets = " << ratj12near4 << endl;

  // statistics for b-jets
  double ratBbyall = 0.;
  if (fNEvtJets > 0) ratBbyall = (double)fNEvtBJets / (double)(fNEvtJets-fNEvt1Jets);
  cout << endl;
  cout << " Number of b-jet events (bare) = " << fNEvtBJets << " Ratio to >=2 jets = " << ratBbyall << endl;
  if (fNEvtBJets <= 0) return;

  // compute the btag efficiencies
  double pTav1B = fPtSum1B / (double)fNEvt1BJets;
  double dpTav1B =  sqrt(fPtSumsq1B / (double)fNEvt1BJets - pTav1B*pTav1B);
  double pTav2B = 0.5 * fPtSum2B / (double)fNEvt2BJets;
  double dpTav2B =  sqrt(0.5 * fPtSumsq2B / (double)fNEvt2BJets - pTav2B*pTav2B);
  cout << "  Average pT in 1b events = " << pTav1B << ", r.m.s. = " << dpTav1B << endl;
  cout << "  Average pT in 2b events = " << pTav2B << ", r.m.s. = " << dpTav2B << endl;
  // b-tagging efficiency
  double effc, effuds;
  GetEfficB (pTav1B, fEffb, fdEffb, effc, effuds);
  double p = fEffb;
  double dp = fdEffb;
  double alp = 0.02 + 0.03 * effc / fEffb;
  double dalp = 0.007;
  fFakb = effuds;
  double palp = p * alp;
  double dpalp = sqrt(dp*dp/(p*p) + dalp*dalp/(alp*alp)) * p*alp;
  cout << "  b-tag effic, p = " << p << " +- " << dp << ", b-tag effic for c-jets = " << effc << endl;
  cout << "  GS prob, alpha = " << alp << " +- " << dalp;
  cout << ", p.alpha = " << p*alp << " +- " << dpalp << endl;
  cout << "  fake b probability (uds) = " << fFakb << endl;

  // compute numbers of b-jet events corrected for fakes
  double nEvt1BJets = fNEvt1BJets
    - (fNEvt1Jets+2.*fNEvt2Jets+3.*fNEvt3Jets+4.*fNEvt4Jets+5.*fNEvt5Jets+6.*fNEvt6Jets+7.*fNEvt7Jets+8.*fNEvt8Jets)*fFakb;
  double nEvt2BJets = fNEvt2BJets
    - (fNEvt2Jets+3.*fNEvt3Jets+6.*fNEvt4Jets+10.*fNEvt5Jets+15.*fNEvt6Jets+20.*fNEvt7Jets+27.*fNEvt8Jets)*fFakb*fFakb;
  double nEvt3BJets = fNEvt3BJets
    - (fNEvt3Jets+4.*fNEvt4Jets+10.*fNEvt5Jets+20.*fNEvt6Jets+35.*fNEvt7Jets+56.*fNEvt8Jets)*fFakb*fFakb*fFakb;
  double nEvt4BJets = fNEvt4BJets
    - (fNEvt4Jets+5.*fNEvt5Jets+15.*fNEvt6Jets+35.*fNEvt7Jets+70.*fNEvt8Jets)*fFakb*fFakb*fFakb*fFakb;
  double nEvtBJets = nEvt1BJets + nEvt2BJets + nEvt3BJets + nEvt4BJets;
  cout << endl;
  cout << " Number of events with 1-4 b-jets (corr fakes = " << fFakb << ") = " << nEvtBJets << endl;
  double rat1bjbyj = 0.;
  if (nEvtBJets > 0.) rat1bjbyj = nEvt1BJets/ nEvtBJets;
  cout << "  Number of events with 1 b-jets uncorr = " << fNEvt1BJets << ", corr = "
       << nEvt1BJets << " Ratio to b-jets = " << rat1bjbyj << endl;
  double rat2bjbyj = 0.;
  if (nEvtBJets > 0.) rat2bjbyj = nEvt2BJets/ nEvtBJets;
  cout << "  Number of events with 2 b-jets uncorr = " << fNEvt2BJets << ", corr = "
       << nEvt2BJets << " Ratio to b-jets = " << rat2bjbyj << endl;
  double rat3bjbyj = 0.;
  if (nEvtBJets > 0.) rat3bjbyj = nEvt3BJets/ nEvtBJets;
  cout << "  Number of events with 3 b-jets uncorr = " << fNEvt3BJets << ", corr = "
       << nEvt3BJets << " Ratio to b-jets = " << rat3bjbyj << endl;
  double rat4bjbyj = 0.;
  if (nEvtBJets > 0.) rat4bjbyj = nEvt4BJets/ nEvtBJets;
  cout << "  Number of events with 4 b-jets uncorr = " << fNEvt4BJets << ", corr = "
       << nEvt4BJets << " Ratio to b-jets = " << rat4bjbyj << endl;

  double ngs[4], dngs[4], ncorr[4], dncorr[4], npair[4],dnpair[4];

  double nbmul[4];
  nbmul[0] = nEvt1BJets;
  nbmul[1] = nEvt2BJets;
  nbmul[2] = nEvt3BJets;
  nbmul[3] = nEvt4BJets;
  double n2j = fNEvtJets  - fNEvt1Jets;

  // compute event numbers using given alpha
  double gfact = 2. * palp;
  ngs[0] = gfact * n2j;
  ngs[1] = gfact * nbmul[0];
  ngs[2] = gfact * nbmul[1];
  ngs[3] = gfact * nbmul[2];
  dngs[0] = 2. * sqrt(dpalp*dpalp*n2j*n2j + palp*palp*n2j);
  dngs[1] = 2. * sqrt(dpalp*dpalp*nbmul[0]*nbmul[0] + palp*palp*nbmul[0]);
  dngs[2] = 2. * sqrt(dpalp*dpalp*nbmul[1]*nbmul[1] + palp*palp*nbmul[1]);
  dngs[3] = 2. * sqrt(dpalp*dpalp*nbmul[2]*nbmul[2] + palp*palp*nbmul[2]);

  ncorr[0] = nbmul[0] - ngs[0];
  ncorr[1] = nbmul[1] - ngs[1];
  ncorr[2] = nbmul[2] - ngs[2];
  ncorr[3] = nbmul[3] - ngs[3];
  dncorr[0] = sqrt(nbmul[0] + dngs[0]*dngs[0]);
  dncorr[1] = sqrt(nbmul[1] + dngs[1]*dngs[1]);
  dncorr[2] = sqrt(nbmul[2] + dngs[2]*dngs[2]);
  dncorr[3] = sqrt(nbmul[3] + dngs[3]*dngs[3]);

  BSolve (p, dp, ncorr, dncorr, npair, dnpair);

  cout << " Nber of GS events 1b = " << ngs[0] << " +- " << dngs[0];
  cout << " Corr nber of 1b events = " << ncorr[0] << " +- " << dncorr[0] << endl;
  cout << " Nber of GS events 2b = " << ngs[1] << " +- " << dngs[1];
  cout << " Corr nber of 2b events = " << ncorr[1] << " +- " << dncorr[1] << endl;
  cout << " Nber of GS events 3b = " << ngs[2] << " +- " << dngs[2];
  cout << " Corr nber of 3b events = " << ncorr[2] << " +- " << dncorr[2] << endl;
  cout << " Nber of GS events 4b = " << ngs[3] << " +- " << dngs[3];
  cout << " Corr nber of 4b events = " << ncorr[3] << " +- " << dncorr[3] << endl;
  cout << " True nber of 1b events = " << npair[0] << " +- " << dnpair[0] << endl;
  cout << " True nber of 2b events = " << npair[1] << " +- " << dnpair[1] << endl;
  cout << " True nber of 3b events = " << npair[2] << " +- " << dnpair[2] << endl;
  cout << " True nber of 4b events = " << npair[3] << " +- " << dnpair[3] << endl;
  // compute alpha from the data
  double alp1, dalp1, dpalp1; 
  GetGSFE (n2j, nbmul, p, dp, alp1, dalp1, ngs, dngs, ncorr, dncorr);

  cout << endl;
  cout << " Alpha computed from data = " << alp1 << " +- " << dalp1 << endl;

  BSolve (p, dp, ncorr, dncorr, npair, dnpair);

  ncorr[0] = nbmul[0] - ngs[0];
  ncorr[1] = nbmul[1] - ngs[1];
  ncorr[2] = nbmul[2] - ngs[2];
  ncorr[3] = nbmul[3] - ngs[3];
  dncorr[0] = sqrt(nbmul[0] + dngs[0]*dngs[0]);
  dncorr[1] = sqrt(nbmul[1] + dngs[1]*dngs[1]);
  dncorr[2] = sqrt(nbmul[2] + dngs[2]*dngs[2]);
  dncorr[3] = sqrt(nbmul[3] + dngs[3]*dngs[3]);
  cout << " Nber of GS events 1b = " << ngs[0] << " +- " << dngs[0];
  cout << " Corr nber of 1b events = " << ncorr[0] << " +- " << dncorr[0] << endl;
  cout << " Nber of GS events 2b = " << ngs[1] << " +- " << dngs[1];
  cout << " Corr nber of 2b events = " << ncorr[1] << " +- " << dncorr[1] << endl;
  cout << " Nber of GS events 3b = " << ngs[2] << " +- " << dngs[2];
  cout << " Corr nber of 3b events = " << ncorr[2] << " +- " << dncorr[2] << endl;
  cout << " Nber of GS events 4b = " << ngs[3] << " +- " << dngs[3];
  cout << " Corr nber of 4b events = " << ncorr[3] << " +- " << dncorr[3] << endl;
  cout << " True nber of 1b events = " << npair[0] << " +- " << dnpair[0] << endl;
  cout << " True nber of 2b events = " << npair[1] << " +- " << dnpair[1] << endl;
  cout << " True nber of 3b events = " << npair[2] << " +- " << dnpair[2] << endl;
  cout << " True nber of 4b events = " << npair[3] << " +- " << dnpair[3] << endl;

  // Jet multiplicities for 1l events
  double nEvt1J1l = fNEvt1J1l;
  double nEvt2J1l = fNEvt2J1l;
  double nEvt3J1l = fNEvt3J1l;
  double nEvt4J1l = fNEvt4J1l;
  double nEvt5J1l = fNEvt5J1l;
  double nEvt6J1l = fNEvt6J1l;
  double nEvt7J1l = fNEvt7J1l;
  double nEvt8J1l = fNEvt8J1l;
  double nEvtJ1l = nEvt1J1l + nEvt2J1l + nEvt3J1l + nEvt4J1l + nEvt5J1l + nEvt6J1l + nEvt7J1l + nEvt8J1l;
  double rat1lbyall = 0.;
  if (fNEvtJets > 0) rat1lbyall = nEvtJ1l / (double)fNEvtJets;
  cout << endl;
  cout << " Number of 1l events with jets = " << nEvtJ1l << " Ratio to all jets = " << rat1lbyall << endl;
  if (nEvtJ1l > 0.) {
    float rat1jby1l = nEvt1J1l / nEvtJ1l;
    cout << "  Number of 1l events with 1 jets = " << nEvt1J1l << " Ratio to 1l+jets = " << rat1jby1l << endl;
    float rat2jby1l = nEvt2J1l / nEvtJ1l;
    cout << "  Number of 1l events with 2 jets = " << nEvt2J1l << " Ratio to 1l+jets = " << rat2jby1l << endl;
    float rat3jby1l = nEvt3J1l / nEvtJ1l;
    cout << "  Number of 1l events with 3 jets = " << nEvt3J1l << " Ratio to 1l+jets = " << rat3jby1l << endl;
    float rat4jby1l = nEvt4J1l / nEvtJ1l;
    cout << "  Number of 1l events with 4 jets = " << nEvt4J1l << " Ratio to 1l+jets = " << rat4jby1l << endl;
    float rat5jby1l = nEvt5J1l / nEvtJ1l;
    cout << "  Number of 1l events with 5 jets = " << nEvt5J1l << " Ratio to 1l+jets = " << rat5jby1l << endl;
    float rat6jby1l = nEvt6J1l / nEvtJ1l;
    cout << "  Number of 1l events with 6 jets = " << nEvt6J1l << " Ratio to 1l+jets = " << rat6jby1l << endl;
    float rat7jby1l = nEvt7J1l / nEvtJ1l;
    cout << "  Number of 1l events with 7 jets = " << nEvt7J1l << " Ratio to 1l+jets = " << rat7jby1l << endl;
    float rat8jby1l = nEvt8J1l / nEvtJ1l;
    cout << "  Number of 1l events with >=8 jets = " << nEvt8J1l << " Ratio to 1l+jets = " << rat8jby1l << endl;

    // Compute number of 1l events with b-jets corrected for fakes
    double nEvt1B1l = fNEvt1B1l
      - (nEvt1J1l+2.*nEvt2J1l+3.*nEvt3J1l+4.*nEvt4J1l+5.*nEvt5J1l+6.*nEvt6J1l+7.*nEvt7J1l+8.*nEvt8J1l)*fFakb;
    double nEvt2B1l = fNEvt2B1l
      - (nEvt2J1l+3.*nEvt3J1l+6.*nEvt4J1l+10.*nEvt5J1l+15.*nEvt6J1l+20.*nEvt7J1l+27.*nEvt8J1l)*fFakb*fFakb;
    double nEvt3B1l = fNEvt3B1l
      - (nEvt3J1l+4.*nEvt4J1l+10.*nEvt5J1l+20.*nEvt6J1l+35.*nEvt7J1l+56.*nEvt8J1l)*fFakb*fFakb*fFakb;
    double nEvt4B1l = fNEvt4B1l
      - (nEvt4J1l+5.*nEvt5J1l+15.*nEvt6J1l+35.*nEvt7J1l+70.*nEvt8J1l)*fFakb*fFakb*fFakb*fFakb;
    double nEvtB1l = nEvt1B1l + nEvt2B1l + nEvt3B1l + nEvt4B1l;
    cout << endl;
    cout << " Number of 1l events with 1-4 b-jets (corr fakes = " << fFakb << ") = " << nEvtB1l << endl;
    double rat1bjby1l = 0.;
    if (nEvtB1l > 0.) rat1bjby1l = nEvt1B1l/ nEvtB1l;
    cout << "  Number of 1l events with 1 b-jets uncorr = " << fNEvt1B1l << ", corr = "
	 << nEvt1B1l << " Ratio to b-jets + 1l = " << rat1bjby1l << endl;
    double rat2bjby1l = 0.;
    if (nEvtB1l > 0.) rat2bjby1l = nEvt2B1l/ nEvtB1l;
    cout << "  Number of 1l events with 2 b-jets uncorr = " << fNEvt2B1l << ", corr = "
	 << nEvt2B1l << " Ratio to b-jets + 1l = " << rat2bjby1l << endl;
    double rat3bjby1l = 0.;
    if (nEvtB1l > 0.) rat3bjby1l = nEvt3B1l/ nEvtB1l;
    cout << "  Number of 1l events with 3 b-jets uncorr = " << fNEvt3B1l << ", corr = "
	 << nEvt3B1l << " Ratio to b-jets + 1l = " << rat3bjby1l << endl;
    double rat4bjby1l = 0.;
    if (nEvtB1l > 0.) rat4bjby1l = nEvt4B1l/ nEvtB1l;
    cout << "  Number of 1l events with 4 b-jets uncorr = " << fNEvt4B1l << ", corr = "
	 << nEvt4B1l << " Ratio to b-jets + 1l = " << rat4bjby1l << endl;
  }

}

void RatioAnalysis::GetEfficB (double pt, double &effb, double &deffb, double &effc, double &effuds) {

  double tabPt[]    = {    30.,     50.,     70.,   100.,   170.,   200.};
  double tabEffb[]  = {   0.25,    0.42,    0.48,   0.55,   0.55,   0.48};
  double tabdEffb[] = {   0.02,    0.03,    0.03,   0.03,   0.03,   0.03};
  double tabEffc[]  = {   0.03,    0.07,    0.08,   0.09,   0.09,   0.09};
  double tabEffu[]  = {0.00018, 0.00048, 0.00076, 0.0011, 0.0018, 0.0020};

  if (pt >= tabPt[0]) {
    int ifnd = -1;
    for (int i = 0; i < 5; ++i) {
      if (pt < tabPt[i+1]) {
	ifnd = i;
	break;
      }
    }
    if (ifnd < 0) ifnd = 4;
    double rat = (pt - tabPt[ifnd]) / (tabPt[ifnd+1] - tabPt[ifnd]);
    effb = tabEffb[ifnd] + rat*(tabEffb[ifnd+1] - tabEffb[ifnd]);
    deffb = tabdEffb[ifnd] + rat*(tabdEffb[ifnd+1] - tabdEffb[ifnd]);
    effc = tabEffc[ifnd] + rat*(tabEffc[ifnd+1] - tabEffc[ifnd]);
    effuds = tabEffu[ifnd] + rat*(tabEffu[ifnd+1] - tabEffu[ifnd]);
  }
  else {
  }

}


void RatioAnalysis::BSolve (double p, double dp, double ncorr[], double dncorr[],
	     double npair[], double dnpair[]) {
  // p     = the b-jet identification efficiency
  // dp    = uncertainty on the b-jet identification efficiency
  // ncorr[]  = corrected number of events
  // dncorr[] = uncertainty on the corrected number of events
  // Returns:
  // npair[]  = true numbers of events (from inversion of the equations)
  // dnpair[] = uncertainty on the true numbers of events

  // compute the numbers of events from PC
  double nEvt1BJets = ncorr[0];
  double nEvt2BJets = ncorr[1];
  double nEvt3BJets = ncorr[2];
  double nEvt4BJets = ncorr[3];
  if (nEvt1BJets < 0.) nEvt1BJets = 0.;
  if (nEvt2BJets < 0.) nEvt2BJets = 0.;
  if (nEvt3BJets < 0.) nEvt3BJets = 0.;
  if (nEvt4BJets < 0.) nEvt4BJets = 0.;
  //  cout << " Nber of corr events 1b = " << nEvt1BJets << ", 2b = " << nEvt2BJets << ", 3b = " << nEvt3BJets
  //       << ", 4b = " << nEvt4BJets << endl;

  // invert the equations to get the true numbers of events
  double prat = (1.-p) / p;
  npair[0] = ( nEvt1BJets - 2.*nEvt2BJets*prat + 3.*nEvt3BJets*(1-p)*(1.-p)/(p*p*p)*prat*prat
		     - 4.*nEvt4BJets*prat*prat*prat ) / p;
  npair[1] = ( nEvt2BJets - 3.*nEvt3BJets*prat + 6.*nEvt4BJets*prat*prat ) / (p*p);
  npair[2] = ( nEvt3BJets - 4.*nEvt4BJets*prat ) / (p*p*p) ;
  npair[3] = nEvt4BJets/(p*p*p*p);

  // compute the uncertainties on the true numbers of events
  double dnEvt1Bsq = dncorr[0]*dncorr[0];
  double dnEvt2Bsq = dncorr[1]*dncorr[1];
  double dnEvt3Bsq = dncorr[2]*dncorr[2];
  double dnEvt4Bsq = dncorr[3]*dncorr[3];
  dnpair[0] = sqrt(dnEvt1Bsq + 4.*dnEvt2Bsq*prat*prat + 9.*dnEvt3Bsq*prat*prat*prat*prat
			 + 16.*dnEvt4Bsq*prat*prat*prat*prat*prat*prat ) / p;
  dnpair[1] = sqrt(dnEvt2Bsq + 9.*dnEvt3Bsq*prat*prat + 36.*dnEvt4Bsq*prat*prat*prat*prat ) / (p*p);
  dnpair[2] = sqrt(dnEvt3Bsq + 16.*dnEvt4Bsq*prat*prat ) / (p*p*p);
  dnpair[3] = sqrt(dnEvt4Bsq ) / (p*p*p*p);

  
  return;

}

void RatioAnalysis::GetGSFE (double n2j, double nbmul[], double p, double dp, 
	      double &alp, double &dalp, double ngs[], double dngs[], double ncorr[], double dncorr[]) {
  // n2j   = the number of events with >= 2 jets
  // nbmul = the number of events in b-jet multiplicities (1, 2, 3, 4)
  // p     = the b-jet identification efficiency
  // dp    = uncertainty on the b-jet identification efficiency
  // Returns:
  // alp   = the gluon splitting probability alpha
  // dalp  = uncertainty on the gluon splitting probability alpha
  // ngs[]    = number of GS+FE events
  // dngs[]   = uncertainty on the number of GS+FE events

  // compute the value of alpha
  double prat = (1.-p) / p;
  double a = nbmul[0] - 2.*prat*nbmul[1] + 3.*prat*prat*nbmul[2] - 4.*prat*prat*prat*nbmul[3];
  double b = n2j - 2.*prat*nbmul[0] + 3.*prat*prat*nbmul[1] - 4.*prat*prat*prat*nbmul[2];
  //  cout << " a = " << a << ", b = " << b << endl;
  alp = a / (2. * p * b);

  // compute the uncertainty on alpha
  double fact = 1. / (2.*p*b);
  double dalpd1b = fact * (1. + 4.*alp*(1.-p));
  double dalpd2b = -fact * prat * (2. + 6.*alp*(1.-p));
  double dalpd3b = fact * prat * prat * (3. + 8.*alp*(1.-p));
  double dalpd4b = -4.*fact * prat * prat * prat;
  double dalpdn2j = -fact * a / b;
  double dadp = -2.*nbmul[1]+6.*prat*nbmul[2]-12.*prat*prat*nbmul[3];
  double dbdp = -2.*nbmul[0]+6.*prat*nbmul[1]-12.*prat*prat*nbmul[2];
  //  cout << " dadp = " << dadp << ", dbdp = " << dbdp << endl;
  double dalpdp = -1./(2*p*p) * (a/b + dadp/(p*b) - a*dbdp/(p*b*b) );
  //  cout << " dalpdp = " << dalpdp << endl;
  dalp = sqrt(dalpd1b*dalpd1b*nbmul[0] + dalpd2b*dalpd2b*nbmul[1] + dalpd3b*dalpd3b*nbmul[2]
	      + dalpd4b*dalpd4b*nbmul[3] + dalpdn2j*dalpdn2j*n2j + dalpdp*dalpdp*dp*dp); 
  double dpalp = sqrt(p*p*dalp*dalp + (alp*alp+2.*alp*p*dalpdp)*dp*dp);
  //  cout << " p.alpha = " << p*alp << " +- " << dpalp << endl;

  // compute the number of events from GS+FE
  double palp = p * alp;
  double gfact = 2. * palp;
  ngs[0] = gfact * n2j;
  ngs[1] = gfact * nbmul[0];
  ngs[2] = gfact * nbmul[1];
  ngs[3] = gfact * nbmul[2];

  // compute the uncertainties on the number of events from GS+FE
  cout << " dpalp = " << dpalp << endl;
  //  dngs[0] = 2. * sqrt(dpalp*dpalp*n2j*n2j + palp*palp*n2j);
  dngs[0] = 2. * sqrt(dpalp*dpalp*n2j*n2j + palp*p*n2j*(alp+2.*dalpdn2j*n2j) );
  dngs[1] = 2. * sqrt(dpalp*dpalp*nbmul[0]*nbmul[0] + palp*p*nbmul[0]*(alp+2.*dalpd1b*nbmul[0]) );
  dngs[2] = 2. * sqrt(dpalp*dpalp*nbmul[1]*nbmul[1] + palp*p*nbmul[1]*(alp+2.*dalpd2b*nbmul[1]) );
  dngs[3] = 2. * sqrt(dpalp*dpalp*nbmul[2]*nbmul[2] + palp*p*nbmul[2]*(alp+2.*dalpd3b*nbmul[2]) );

  // compute the corrected number of events
  ncorr[0] = nbmul[0] - ngs[0];
  ncorr[1] = nbmul[1] - ngs[1];
  ncorr[2] = nbmul[2] - ngs[2];
  ncorr[3] = nbmul[3] - ngs[3];
  dncorr[0] = sqrt(nbmul[0] + dngs[0]*dngs[0]);
  dncorr[1] = sqrt(nbmul[1] + dngs[1]*dngs[1]);
  dncorr[2] = sqrt(nbmul[2] + dngs[2]*dngs[2]);
  dncorr[3] = sqrt(nbmul[3] + dngs[3]*dngs[3]);

  return;
}


void RatioAnalysis::PrintChaRatios() {

  // Charge Ratios for all dilepton events
  cout << endl;
  cout << " Charge ratios for all events of chargino type " << endl;
  cout << " ============================================= " << endl;
  cout << "  Effic_e = " << fEffe << ", Effic_m = " << fEffm << endl;
  cout << "  FakeR_e = " << fFake << ", FakeR_m = " << fFakm << endl;
  cout << endl;

  // compute numbers of charge configs corrected for fakes
  double effee = fEffe*fEffe;
  double effmm = fEffm*fEffm;
  double effem = fEffe*fEffm;
  double fakee = (1.-fFake)*(1.-fFake);
  double fakmm = (1.-fFakm)*(1.-fFakm);
  double fakem = (1.-fFake)*(1.-fFakm);
  double nEvtepep = fNEvtepep*fakee/effee;
  double nEvtmpmp = fNEvtmpmp*fakmm/effmm;
  double nEvtepmp = fNEvtepmp*fakem/effem;
  double nEvtenen = fNEvtenen*fakee/effee;
  double nEvtmnmn = fNEvtmnmn*fakmm/effmm;
  double nEvtenmn = fNEvtenmn*fakem/effem;
  double nEvtepmn = fNEvtepmn*fakem/effem;
  double nEvtenmp = fNEvtenmp*fakem/effem;
  double dnEvtepep = sqrt(fNEvtepep)*fakee/effee;
  double dnEvtmpmp = sqrt(fNEvtmpmp)*fakmm/effmm;
  double dnEvtepmp = sqrt(fNEvtepmp)*fakem/effem;
  double dnEvtenen = sqrt(fNEvtenen)*fakee/effee;
  double dnEvtmnmn = sqrt(fNEvtmnmn)*fakmm/effmm;
  double dnEvtenmn = sqrt(fNEvtenmn)*fakem/effem;
  double dnEvtepmn = sqrt(fNEvtepmn)*fakem/effem;
  double dnEvtenmp = sqrt(fNEvtenmp)*fakem/effem;

  int nChappunc = fNEvtepep + fNEvtmpmp + fNEvtepmp;
  int nChannunc = fNEvtenen + fNEvtmnmn + fNEvtenmn;
  int nChapnunc = 2. * (fNEvtepmn + fNEvtenmp);
  double nChapp = nEvtepep + nEvtmpmp + nEvtepmp;
  double nChann = nEvtenen + nEvtmnmn + nEvtenmn;
  double nChapn = 2. * (nEvtepmn + nEvtenmp);
  double dnChapp = sqrt(dnEvtepep*dnEvtepep + dnEvtmpmp*dnEvtmpmp + dnEvtepmp*dnEvtepmp);
  double dnChann = sqrt(dnEvtenen*dnEvtenen + dnEvtmnmn*dnEvtmnmn + dnEvtenmn*dnEvtenmn);
  double dnChapn = 2. * sqrt(dnEvtepmn*dnEvtepmn + dnEvtenmp*dnEvtenmp);
  cout << "  Events with ++ uncorr = " << nChappunc << ", corr = " << nChapp << " +- " << dnChapp << endl;
  cout << "  Events with -- uncorr = " << nChannunc << ", corr = " << nChann << " +- " << dnChann << endl;
  cout << "  Events with +- uncorr = " << nChapnunc << ", corr = " << nChapn << " +- " << dnChapn << endl;
  double ratSS = 0.;
  if (nChann != 0.) ratSS = nChapp / nChann;
  double dratSS = 0.;
  if (ratSS != 0.) dratSS = sqrt(dnChapp*dnChapp/(nChapp*nChapp) + dnChann*dnChann/(nChann*nChann) ) * ratSS;
  double sumppnn = nChapp + nChann;
  double dsumppnn = sqrt(dnChapp*dnChapp + dnChann*dnChann);
  double ratSSOS = 0.;
  if (nChapn != 0.) ratSSOS = sumppnn / nChapn;
  double dratSSOS = 0.;
  if (ratSSOS != 0.) dratSSOS = sqrt(dsumppnn*dsumppnn/(sumppnn*sumppnn) + dnChapn*dnChapn/(nChapn*nChapn) ) * ratSSOS;
  cout << "  Ratio ++/--        = " << ratSS   << " +- " << dratSS << endl;
  cout << "  Ratio (++ + --)/+- = " << ratSSOS << " +- " << dratSSOS << endl;

  cout << endl;
  cout << "  Events with e+e+ uncorr = " << fNEvtepep << ", corr = " << nEvtepep << " +- " << dnEvtepep << endl;
  cout << "  Events with e-e- uncorr = " << fNEvtenen << ", corr = " << nEvtenen << " +- " << dnEvtenen << endl;
  cout << "  Events with m+m+ uncorr = " << fNEvtmpmp << ", corr = " << nEvtmpmp << " +- " << dnEvtmpmp << endl;
  cout << "  Events with m-m- uncorr = " << fNEvtmnmn << ", corr = " << nEvtmnmn << " +- " << dnEvtmnmn << endl;


  // Ratios for dilepton events with b-jets
  cout << endl;
  cout << " Charge ratios for events with b-jets " << endl;

  // compute numbers of charge configs corrected for fakes
  double nEvtepepB = fNEvtepepB*fakee/effee;
  double nEvtmpmpB = fNEvtmpmpB*fakmm/effmm;
  double nEvtepmpB = fNEvtepmpB*fakem/effem;
  double nEvtenenB = fNEvtenenB*fakee/effee;
  double nEvtmnmnB = fNEvtmnmnB*fakmm/effmm;
  double nEvtenmnB = fNEvtenmnB*fakem/effem;
  double nEvtepmnB = fNEvtepmnB*fakem/effem;
  double nEvtenmpB = fNEvtenmpB*fakem/effem;
  double dnEvtepepB = sqrt(fNEvtepepB)*fakee/effee;
  double dnEvtmpmpB = sqrt(fNEvtmpmpB)*fakmm/effmm;
  double dnEvtepmpB = sqrt(fNEvtepmpB)*fakem/effem;
  double dnEvtenenB = sqrt(fNEvtenenB)*fakee/effee;
  double dnEvtmnmnB = sqrt(fNEvtmnmnB)*fakmm/effmm;
  double dnEvtenmnB = sqrt(fNEvtenmnB)*fakem/effem;
  double dnEvtepmnB = sqrt(fNEvtepmnB)*fakem/effem;
  double dnEvtenmpB = sqrt(fNEvtenmpB)*fakem/effem;

  int nChappuncB = fNEvtepepB + fNEvtmpmpB + fNEvtepmpB;
  int nChannuncB = fNEvtenenB + fNEvtmnmnB + fNEvtenmnB;
  int nChapnuncB = 2. * (fNEvtepmnB + fNEvtenmpB);
  double nChappB = nEvtepepB + nEvtmpmpB + nEvtepmpB;
  double nChannB = nEvtenenB + nEvtmnmnB + nEvtenmnB;
  double nChapnB = 2. * (nEvtepmnB + nEvtenmpB);
  double dnChappB = sqrt(dnEvtepepB*dnEvtepepB + dnEvtmpmpB*dnEvtmpmpB + dnEvtepmpB*dnEvtepmpB);
  double dnChannB = sqrt(dnEvtenenB*dnEvtenenB + dnEvtmnmnB*dnEvtmnmnB + dnEvtenmnB*dnEvtenmnB);
  double dnChapnB = 2. * sqrt(dnEvtepmnB*dnEvtepmnB + dnEvtenmpB*dnEvtenmpB);
  cout << "  Events with ++ uncorr = " << nChappuncB << ", corr = " << nChappB << " +- " << dnChappB << endl;
  cout << "  Events with -- uncorr = " << nChannuncB << ", corr = " << nChannB << " +- " << dnChannB << endl;
  cout << "  Events with +- uncorr = " << nChapnuncB << ", corr = " << nChapnB << " +- " << dnChapnB << endl;
  double ratSSB = 0.;
  if (nChannB != 0.) ratSSB = nChappB / nChannB;
  double dratSSB = 0.;
  if (ratSSB != 0.)
    dratSSB = sqrt(dnChappB*dnChappB/(nChappB*nChappB) + dnChannB*dnChannB/(nChannB*nChannB) ) * ratSSB;
  double sumppnnB = nChappB + nChannB;
  double dsumppnnB = sqrt(dnChappB*dnChappB + dnChannB*dnChannB);
  double ratSSOSB = 0.;
  if (nChapnB != 0.) ratSSOSB = sumppnnB / nChapnB;
  double dratSSOSB = 0.;
  if (ratSSOSB != 0.) 
    dratSSOSB = sqrt(dsumppnnB*dsumppnnB/(sumppnnB*sumppnnB) + dnChapnB*dnChapnB/(nChapnB*nChapnB) ) * ratSSOSB;
  cout << "  Ratio ++/--        = " << ratSSB   << " +- " << dratSSB << endl;
  cout << "  Ratio (++ + --)/+- = " << ratSSOSB << " +- " << dratSSOSB << endl;

}

void RatioAnalysis::PrintDilRatios() {

  // Charginos/neutralinos in all dilepton events
  cout << endl;
  cout << " Charginos/neutralinos in dilepton events " << endl;
  cout << " ======================================== " << endl;
  cout << "  Effic_e = " << fEffe << ", Effic_m = " << fEffm << endl;
  cout << "  FakeR_e = " << fFake << ", FakeR_m = " << fFakm << endl;
  cout << endl;

  // compute numbers of charge configs corrected for fakes
  double effee = fEffe*fEffe;
  double effmm = fEffm*fEffm;
  double effem = fEffe*fEffm;
  double fakee = (1.-fFake)*(1.-fFake);
  double fakmm = (1.-fFakm)*(1.-fFakm);
  double fakem = (1.-fFake)*(1.-fFakm);
  double nEvtepep = fNEvtepep*fakee/effee;
  double nEvtmpmp = fNEvtmpmp*fakmm/effmm;
  double nEvtepmp = fNEvtepmp*fakem/effem;
  double nEvtenen = fNEvtenen*fakee/effee;
  double nEvtmnmn = fNEvtmnmn*fakmm/effmm;
  double nEvtenmn = fNEvtenmn*fakem/effem;
  double nEvtepmn = fNEvtepmn*fakem/effem;
  double nEvtenmp = fNEvtenmp*fakem/effem;
  double nEvtepen = fNEvtepen*fakee/effee;
  double nEvtmpmn = fNEvtmpmn*fakmm/effmm;
  double dnEvtepep = sqrt(fNEvtepep)*fakee/effee;
  double dnEvtmpmp = sqrt(fNEvtmpmp)*fakmm/effmm;
  double dnEvtepmp = sqrt(fNEvtepmp)*fakem/effem;
  double dnEvtenen = sqrt(fNEvtenen)*fakee/effee;
  double dnEvtmnmn = sqrt(fNEvtmnmn)*fakmm/effmm;
  double dnEvtenmn = sqrt(fNEvtenmn)*fakem/effem;
  double dnEvtepmn = sqrt(fNEvtepmn)*fakem/effem;
  double dnEvtenmp = sqrt(fNEvtenmp)*fakem/effem;
  double dnEvtepen = sqrt(fNEvtepen)*fakee/effee;
  double dnEvtmpmn = sqrt(fNEvtmpmn)*fakmm/effmm;
  double nEvtJ2l = fNEvtepep + fNEvtepen + fNEvtenen + fNEvtmpmp + fNEvtmpmn + fNEvtmnmn +
                   fNEvtepmp + fNEvtepmn + fNEvtenmp + fNEvtenmn;

  double rat2lbyall = 0.;
  if (fNEvtJets > 0) rat2lbyall = nEvtJ2l / (double)fNEvtJets;
  cout << " Number of 2l events with jets = " << nEvtJ2l << " Ratio to all jets = " << rat2lbyall << endl;

  int nDilppunc = fNEvtepep + fNEvtmpmp + fNEvtepmp;
  int nDilnnunc = fNEvtenen + fNEvtmnmn + fNEvtenmn;
  int nDilChinounc = nDilppunc + nDilnnunc;
  double nDilpp = nEvtepep + nEvtmpmp + nEvtepmp;
  double nDilnn = nEvtenen + nEvtmnmn + nEvtenmn;
  double nDilChino = nDilpp + nDilnn;
  double dnDilpp = sqrt(dnEvtepep*dnEvtepep + dnEvtmpmp*dnEvtmpmp + dnEvtepmp*dnEvtepmp);
  double dnDilnn = sqrt(dnEvtenen*dnEvtenen + dnEvtmnmn*dnEvtmnmn + dnEvtenmn*dnEvtenmn);
  double dnDilChino = sqrt(dnDilpp*dnDilpp + dnDilnn*dnDilnn);
  cout << "  Events with ++ uncorr = " << nDilppunc << ", corr = " << nDilpp << " +- " << dnDilpp << endl;
  cout << "  Events with -- uncorr = " << nDilnnunc << ", corr = " << nDilnn << " +- " << dnDilnn << endl;
  cout << "  SS Chargino candidates uncorr = " << nDilChinounc << ", corr = " << nDilChino << " +- " << dnDilChino << endl;

  int nDilemunc = fNEvtepmn + fNEvtenmp;
  int nDilNinounc = fNEvtepen + fNEvtmpmn - nDilemunc;
  double nDilepen = nEvtepen;
  double nDilmpmn = nEvtmpmn;
  double nDilpn = nDilepen + nDilmpmn;
  double nDilem = nEvtepmn + nEvtenmp;
  double nDilNino = nDilpn - nDilem;
  double dnDilepen = dnEvtepen;
  double dnDilmpmn = dnEvtmpmn;
  double dnDilpn = sqrt(dnDilepen*dnDilepen + nDilmpmn*dnDilmpmn);
  double dnDilem = sqrt(dnEvtepmn*dnEvtepmn + dnEvtenmp*dnEvtenmp);
  double dnDilNino = sqrt(dnDilpn*dnDilpn + nDilem*dnDilem);
  cout << "  Events with e+e- uncorr  = " << fNEvtepen << ", corr = " << nDilepen << " +- " << dnDilepen << endl;
  cout << "  Events with m+m- uncorr  = " << fNEvtmpmn << ", corr = " << nDilmpmn << " +- " << dnDilmpmn << endl;
  cout << "  Events with OS em uncorr = " << fNEvtepmn + fNEvtenmp << ", corr = " << nDilem << " +- " << dnDilem << endl;
  cout << "  Neutralino candidates uncorr = " << nDilNino << ", corr = " << nDilNino << " +- " << dnDilNino << endl;
  double ratOSembyll = 0.;
  if (nDilpn > 0.) ratOSembyll = nDilem / nDilpn;
  double dratOSembyll = 0.;
  if (ratOSembyll != 0.)
    dratOSembyll = sqrt(dnDilpn*dnDilpn/(nDilpn*nDilpn) + dnDilem*dnDilem/(nDilem*nDilem)) * ratOSembyll;
  cout << "  Ratio OS(em)/(ee+mm) corr = " << ratOSembyll << " +- " << dratOSembyll << endl;


  // Charginos/neutralinos in dilepton events with b-jets
  cout << endl;
  cout << " Charginos/neutralinos in dilepton events with b-jets " << endl;

  // compute numbers of charge configs corrected for fakes
  double nEvtepepB = fNEvtepepB*fakee/effee;
  double nEvtmpmpB = fNEvtmpmpB*fakmm/effmm;
  double nEvtepmpB = fNEvtepmpB*fakem/effem;
  double nEvtenenB = fNEvtenenB*fakee/effee;
  double nEvtmnmnB = fNEvtmnmnB*fakmm/effmm;
  double nEvtenmnB = fNEvtenmnB*fakem/effem;
  double nEvtepmnB = fNEvtepmnB*fakem/effem;
  double nEvtenmpB = fNEvtenmpB*fakem/effem;
  double nEvtepenB = fNEvtepenB*fakee/effee;
  double nEvtmpmnB = fNEvtmpmnB*fakmm/effmm;
  double dnEvtepepB = sqrt(fNEvtepepB)*fakee/effee;
  double dnEvtmpmpB = sqrt(fNEvtmpmpB)*fakmm/effmm;
  double dnEvtepmpB = sqrt(fNEvtepmpB)*fakem/effem;
  double dnEvtenenB = sqrt(fNEvtenenB)*fakee/effee;
  double dnEvtmnmnB = sqrt(fNEvtmnmnB)*fakmm/effmm;
  double dnEvtenmnB = sqrt(fNEvtenmnB)*fakem/effem;
  double dnEvtepmnB = sqrt(fNEvtepmnB)*fakem/effem;
  double dnEvtenmpB = sqrt(fNEvtenmpB)*fakem/effem;
  double dnEvtepenB = sqrt(fNEvtepenB)*fakee/effee;
  double dnEvtmpmnB = sqrt(fNEvtmpmnB)*fakmm/effmm;
  double nEvtJ2lB = fNEvtepepB + fNEvtepenB + fNEvtenenB + fNEvtmpmpB + fNEvtmpmnB + fNEvtmnmnB +
                    fNEvtepmpB + fNEvtepmnB + fNEvtenmpB + fNEvtenmnB;

  double rat2lbyallB = 0.;
  if (fNEvtBJets > 0) rat2lbyallB = nEvtJ2lB / (double)fNEvtBJets;
  cout << " Number of 2l events with b-jets = " << nEvtJ2lB << " Ratio to all b-jets = " << rat2lbyallB << endl;

  int nDilppuncB = fNEvtepepB + fNEvtmpmpB + fNEvtepmpB;
  int nDilnnuncB = fNEvtenenB + fNEvtmnmnB + fNEvtenmnB;
  int nDilChinouncB = nDilppuncB + nDilnnuncB;
  double nDilppB = nEvtepepB + nEvtmpmpB + nEvtepmpB;
  double nDilnnB = nEvtenenB + nEvtmnmnB + nEvtenmnB;
  double nDilChinoB = nDilppB + nDilnnB;
  double dnDilppB = sqrt(dnEvtepepB*dnEvtepepB + dnEvtmpmpB*dnEvtmpmpB + dnEvtepmpB*dnEvtepmpB);
  double dnDilnnB = sqrt(dnEvtenenB*dnEvtenenB + dnEvtmnmnB*dnEvtmnmnB + dnEvtenmnB*dnEvtenmnB);
  double dnDilChinoB = sqrt(dnDilppB*dnDilppB + dnDilnnB*dnDilnnB);
  cout << "  Events with ++ uncorr = " << nDilppuncB << ", corr = " << nDilppB << " +- " << dnDilppB << endl;
  cout << "  Events with -- ucorr = " << nDilnnuncB << ", corr = " << nDilnnB << " +- " << dnDilnnB << endl;
  cout << "  SS Chargino candidates uncorr = " << nDilChinouncB << ", corr = " << nDilChinoB << " +- " << dnDilChinoB << endl;

  int nDilemuncB = fNEvtepmnB + fNEvtenmpB;
  int nDilNinouncB = fNEvtepenB + fNEvtmpmnB - nDilemuncB;
  double nDilepenB = nEvtepenB;
  double nDilmpmnB = nEvtmpmnB;
  double nDilpnB = nDilepenB + nDilmpmnB;
  double nDilemB = nEvtepmnB + nEvtenmpB;
  double nDilNinoB = nDilpnB - nDilemB;
  double dnDilepenB = dnEvtepenB;
  double dnDilmpmnB = dnEvtmpmnB;
  double dnDilpnB = sqrt(dnDilepenB*dnDilepenB + nDilmpmnB*dnDilmpmnB);
  double dnDilemB = sqrt(dnEvtepmnB*dnEvtepmnB + dnEvtenmpB*dnEvtenmpB);
  double dnDilNinoB = sqrt(dnDilpnB*dnDilpnB + nDilemB*dnDilemB);
  cout << "  Events with e+e-  uncorr  = " << fNEvtepenB << ", corr = " << nDilepenB << " +- " << dnDilepenB << endl;
  cout << "  Events with m+m-  uncorr  = " << fNEvtmpmnB << ", corr = " << nDilmpmnB << " +- " << dnDilmpmnB << endl;
  cout << "  Events with OS em uncorr  = " << nDilemuncB << ", corr = " << nDilemB << " +- " << dnDilemB << endl;
  cout << "  Neutralino candidates uncorr = " << nDilNinouncB << ", corr = " << nDilNinoB << " +- " << dnDilNinoB << endl;
  double ratOSembyllB = 0.;
  if (nDilpnB > 0.) ratOSembyllB = nDilemB / nDilpnB;
  double dratOSembyllB = 0.;
  if (ratOSembyllB != 0.) 
    dratOSembyllB = sqrt(dnDilpnB*dnDilpnB/(nDilpnB*nDilpnB) + dnDilemB*dnDilemB/(nDilemB*nDilemB)) * ratOSembyll;
  cout << "  Ratio OS(em)/(ee+mm) = " << ratOSembyllB << " +- " << dratOSembyllB << endl;

  // b-jet multiplicities in 2l SS, OS, em
  double nEvt1BSS2l = fNEvt1BSS2l;
  double nEvt2BSS2l = fNEvt2BSS2l;
  double nEvt3BSS2l = fNEvt3BSS2l;
  double nEvt4BSS2l = fNEvt4BSS2l;
  double nEvt1BOS2l = fNEvt1BOS2l;
  double nEvt2BOS2l = fNEvt2BOS2l;
  double nEvt3BOS2l = fNEvt3BOS2l;
  double nEvt4BOS2l = fNEvt4BOS2l;
  double nEvt1BOSem = fNEvt1BOSem;
  double nEvt2BOSem = fNEvt2BOSem;
  double nEvt3BOSem = fNEvt3BOSem;
  double nEvt4BOSem = fNEvt4BOSem;
  double dnEvt1BSS2l = sqrt(fNEvt1BSS2l);
  double dnEvt2BSS2l = sqrt(fNEvt2BSS2l);
  double dnEvt3BSS2l = sqrt(fNEvt3BSS2l);
  double dnEvt4BSS2l = sqrt(fNEvt4BSS2l);
  double dnEvt1BOS2l = sqrt(fNEvt1BOS2l);
  double dnEvt2BOS2l = sqrt(fNEvt2BOS2l);
  double dnEvt3BOS2l = sqrt(fNEvt3BOS2l);
  double dnEvt4BOS2l = sqrt(fNEvt4BOS2l);
  double dnEvt1BOSem = sqrt(fNEvt1BOSem);
  double dnEvt2BOSem = sqrt(fNEvt2BOSem);
  double dnEvt3BOSem = sqrt(fNEvt3BOSem);
  double dnEvt4BOSem = sqrt(fNEvt4BOSem);
  cout << endl;
  cout << " b-jet multiplicities in dilepton events " << endl;
  cout << "  Events with SS 2l, 1b = " << nEvt1BSS2l << " +- " << dnEvt1BSS2l << endl;
  cout << "  Events with SS 2l, 2b = " << nEvt2BSS2l << " +- " << dnEvt2BSS2l << endl;
  cout << "  Events with SS 2l, 3b = " << nEvt3BSS2l << " +- " << dnEvt3BSS2l << endl;
  cout << "  Events with SS 2l, 4b = " << nEvt4BSS2l << " +- " << dnEvt4BSS2l << endl;
  cout << endl;
  cout << "  Events with OS 2l, 1b = " << nEvt1BOS2l << " +- " << dnEvt1BOS2l << endl;
  cout << "  Events with OS 2l, 2b = " << nEvt2BOS2l << " +- " << dnEvt2BOS2l << endl;
  cout << "  Events with OS 2l, 3b = " << nEvt3BOS2l << " +- " << dnEvt3BOS2l << endl;
  cout << "  Events with OS 2l, 4b = " << nEvt4BOS2l << " +- " << dnEvt4BOS2l << endl;
  cout << endl;
  cout << "  Events with OS em, 1b = " << nEvt1BOSem << " +- " << dnEvt1BOSem << endl;
  cout << "  Events with OS em, 2b = " << nEvt2BOSem << " +- " << dnEvt2BOSem << endl;
  cout << "  Events with OS em, 3b = " << nEvt3BOSem << " +- " << dnEvt3BOSem << endl;
  cout << "  Events with OS em, 4b = " << nEvt4BOSem << " +- " << dnEvt4BOSem << endl;


}



void RatioAnalysis::PrintAnomEvts() {

  cout << endl;
  cout << " Anomalous events list: " << fIRun.size() << " events " << endl;
  cout << endl;

  if (fIRun.size() <= 0) return;

  int nEvtPrnt = fIRun.size();
  if (nEvtPrnt > 100) nEvtPrnt = 100;
  for (int i = 0; i < nEvtPrnt; ++i) {
    cout << " Run/Lumi/Event " << fIRun[i] << " " << fILumi[i] << " " << fINber[i];
    if (fINlep[i] > 0) cout << ", Nlept = " << fINlep[i] << endl;
    if (fINlep[i] == -1) cout << ", Acop leading jets " << endl;
    if (fINlep[i] == -2) cout << ", 3 b-jets " << endl;
    if (fINlep[i] == -3) cout << ", 4 b-jets " << endl;
 }

}

void RatioAnalysis::SaveRatioHist() {

  cout << " saving the plots" << endl;
  const int nRatBins = 63;
  const char* lablx[nRatBins] = {"evj   ", "ev1j  ", "ev2j  ", "ev3j  ", "ev4j  ", "ev5j  ", "ev6j  ", "ev7j  ", "ev8j  ",
				 "acop1 ", "acop2 ", "acop4",
 				 "evb   ", "ev1b  ", "ev2b  ", "ev3b  ", "ev4b  ",
				 "ev1l  ", "ev11l ", "ev21l ", "ev31l ", "ev41l ", "ev51l ", "ev61l ", "ev71l ", "ev81l ",
				 "evb1l ", "ev1b1l", "ev2b1l", "ev3b1l", "ev4b1l",
				 "mpmp  ", "mnmn  ", "epmp  ", "enmn  ", "epep  ", "enen  ", "mpmn  ", "epmn  ", "enmp  ", "epen  ",
				 "mpmpB ", "mnmnB ", "epmpB ", "enmnB ", "epepB ", "enenB ", "mpmnB ", "epmnB ", "enmpB ", "epenB ",
				 "SS2l1b", "SS2l2b", "SS2l3b", "SS2l4b", 
				 "OS2l1b", "OS2l2b", "OS2l3b", "OS2l4b", 
				 "OSem1b", "OSem2b", "OSem3b", "OSem4b"};
 
  int ratBin[nRatBins];
  ratBin[0] = fNEvtJets;
  ratBin[1] = fNEvt1Jets;
  ratBin[2] = fNEvt2Jets;
  ratBin[3] = fNEvt3Jets;
  ratBin[4] = fNEvt4Jets;
  ratBin[5] = fNEvt5Jets;
  ratBin[6] = fNEvt6Jets;
  ratBin[7] = fNEvt7Jets;
  ratBin[8] = fNEvt8Jets;
  ratBin[9] = fNEvtj12near1;
  ratBin[10] = fNEvtj12near2;
  ratBin[11] = fNEvtj12near4;
  ratBin[12] = fNEvtBJets;
  ratBin[13] = fNEvt1BJets;
  ratBin[14] = fNEvt2BJets;
  ratBin[15] = fNEvt3BJets;
  ratBin[16] = fNEvt4BJets;
  ratBin[17] = fNEvtJ1l;
  ratBin[18] = fNEvt1J1l;
  ratBin[19] = fNEvt2J1l;
  ratBin[20] = fNEvt3J1l;
  ratBin[21] = fNEvt4J1l;
  ratBin[22] = fNEvt5J1l;
  ratBin[23] = fNEvt6J1l;
  ratBin[24] = fNEvt7J1l;
  ratBin[25] = fNEvt8J1l;
  ratBin[26] = fNEvtB1l;
  ratBin[27] = fNEvt1B1l;
  ratBin[28] = fNEvt2B1l;
  ratBin[29] = fNEvt3B1l;
  ratBin[30] = fNEvt4B1l;
  ratBin[31] = fNEvtmpmp;
  ratBin[32] = fNEvtmnmn;
  ratBin[33] = fNEvtepmp;
  ratBin[34] = fNEvtenmn;
  ratBin[35] = fNEvtepep; 
  ratBin[36] = fNEvtenen;
  ratBin[37] = fNEvtmpmn;
  ratBin[38] = fNEvtepmn;
  ratBin[39] = fNEvtenmp;
  ratBin[40] = fNEvtepen;
  ratBin[41] = fNEvtmpmpB;
  ratBin[42] = fNEvtmnmnB;
  ratBin[43] = fNEvtepmpB;
  ratBin[44] = fNEvtenmnB;
  ratBin[45] = fNEvtepepB;
  ratBin[46] = fNEvtenenB;
  ratBin[47] = fNEvtmpmnB;
  ratBin[48] = fNEvtepmnB;
  ratBin[49] = fNEvtenmpB;
  ratBin[50] = fNEvtepenB;
  ratBin[51] = fNEvt1BSS2l;
  ratBin[52] = fNEvt2BSS2l;
  ratBin[53] = fNEvt3BSS2l;
  ratBin[54] = fNEvt4BSS2l;
  ratBin[55] = fNEvt1BOS2l;
  ratBin[56] = fNEvt2BOS2l;
  ratBin[57] = fNEvt3BOS2l;
  ratBin[58] = fNEvt4BOS2l;
  ratBin[59] = fNEvt1BOSem;
  ratBin[60] = fNEvt2BOSem;
  ratBin[61] = fNEvt3BOSem;
  ratBin[62] = fNEvt4BOSem;
	
  for (int i = 0; i < nRatBins; ++i) {
    fRatHist->GetXaxis()->SetBinLabel(i+1, lablx[i]);
    fRatHist->SetBinContent(i+1, ratBin[i]);
  }
  cout << " ratio plot filled" << endl;
  /*
  fTlat->SetTextColor(kBlack);
  fTlat->SetNDC(kTRUE);
  fTlat->SetTextSize(0.04);

  TString subdir = "MultiplicityPlots";
  TString canvtitle;
  TCanvas *canv;

  canvtitle = "RatioStats";
  canv = new TCanvas("RatioStats", canvtitle , 0, 0, 900, 700);
  canv->SetRightMargin(0.15);
  gPad->SetLogy();	
  fRatHist->DrawCopy();
  fTlat->DrawLatex(0.11,0.92, canvtitle);
  Util::Print(canv, fTag + "_ratiostats", fOutputDir+subdir);
  cout << " wrote ratio plots eps/png" << endl;
  */
  fMPHistFile->cd();
  cout << " directory to MPHist set" << endl;
  fRatHist->Write();
  fDRJ12->Write();
  fPtJ1->Write();
  fPtJ2->Write();
  fPtJets->Write();
  fEtaJets->Write();
  fRptDPhi2j->Write();
  fAlpT2j->Write();
  fBProbJets->Write();
  fEta2B->Write();
  fdPhi2B->Write();
  fEta2B->Write();
  fDR2B->Write();
  fBProb2B->Write();
  fNjets2B->Write();
  fMass2B->Write();
  fNjets2Bnear->Write();
  fMET2Bnear->Write();
  fdPhiJBnear->Write();
  fEta1B->Write();
  fBProb1B->Write();
  fNjets1B->Write();
  fMass1B->Write();
  fMETdPhi1B->Write();
  fdPhiMET1B->Write();
  fJMult->Write();
  fJMult1B->Write();
  fJMult2B->Write();
  fJMult3B->Write();
  fJMult4B->Write();

  cout << " wrote all plots" << endl;
  fMPHistFile->Close();
	
  cout << "number of events fired:      " << counter << endl;
  cout << "number of events in dataset: " << fTR->GetEntries() << endl;
}


