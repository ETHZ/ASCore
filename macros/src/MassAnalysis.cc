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

void MassAnalysis::Begin(){
	// Define the output file of histograms
	const char* filename = "Mass_histos.root";
	fHistFile = new TFile(fOutputDir + TString(filename), "RECREATE");
	
	fMT2_histos_step = 50;
	fMT2_histos_number = 6;
	
	// Define the histograms
	fHMT2_vs_M_OSDiLept = new TH2D("MT2_vs_M_OSDiLept"   , "MT2_vs_M_OSDiLept"       , 200, 0, 200, 100, 0, 500);
	fHAlpahT_DiJet      = new TH1D("AlpahT_DiJet"        , "alpha_T for dijets"      , 100, 0., 3.);
	fHAlpahT_PseudoJet  = new TH1D("AlpahT_PseudoJet"    , "alpha_T for pseudojets"  , 100, 0., 3.);
	
	fHMCT_DiJet         = new TH1D("MCT_DiJet"          , "di-Jets MCT"              , 100, 0., 500.);
	fHMCT_DiBJet        = new TH1D("MCT_DiBJet"         , "di-B-Jets MCT"            , 100, 0., 500.);
	fHMCT_OSee          = new TH1D("MCT_OSee"           , "OS ee MCT"                ,  40, 0., 400.);
	fHMCT_OSmumu        = new TH1D("MCT_OSmumu"         , "OS mumu MCT"              ,  40, 0., 400.);
	fHMCT_OSemu         = new TH1D("MCT_OSemu"          , "OS emu MCT"               ,  40, 0., 400.);
	                                                                                 
	fHMCTperp_OSee      = new TH1D("MCTperp_OSee"        , "OS ee MCT perp"          ,  40, 0., 400.);	
	fHMCTperp_OSmumu    = new TH1D("MCTperp_OSmumu"      , "OS mumu MCT perp"        ,  40, 0., 400.);	
	fHMCTperp_OSemu     = new TH1D("MCTperp_OSemu"       , "OS emu MCT perp"         ,  40, 0., 400.);	    
	                                             
	fHMCT_TTbar         = new TH1D("MCT_TTbar"           , "MCT for TTbar on GenLevel"      , 50, 0., 200.);
	fHMCTperp_TTbar     = new TH1D("MCTperp_TTbar"       , "TTbar MCT perp"                 , 50, 0., 200.);
	fHMT2_TTbar         = new TH1D("MT2_TTbar"           , "MT2 for TTbar with Gen Info"    , 50, 0., 200.);
	fHMT2perp_TTbar     = new TH1D("MT2perp_TTbar"       , "MT2perp for TTbar with Gen Info", 50, 0., 200.);
	
	fHInvMassOSee      = new TH1D("InvMass_OSee"  ,  " Inv Mass of OS di-leptons ee"   , 30, 0, 300);
	fHInvMassOSmumu    = new TH1D("InvMass_OSmumu",  " Inv Mass of OS di-leptons mumu" , 30, 0, 300);
	fHInvMassOSemu     = new TH1D("InvMass_OSemu" ,  " Inv Mass of OS di-leptons emu"  , 30, 0, 300);
	fHInvMassOSll      = new TH1D("InvMass_OSll"  ,  " Inv Mass of OS di-leptons ll"   , 30, 0, 300);
	fHInvMassSSee      = new TH1D("InvMass_SSee"  ,  " Inv Mass of SS di-leptons ee"   , 30, 0, 300);
	fHInvMassSSmumu    = new TH1D("InvMass_SSmumu",  " Inv Mass of SS di-leptons mumu" , 30, 0, 300);
	fHInvMassSSemu     = new TH1D("InvMass_SSemu" ,  " Inv Mass of SS di-leptons emu"  , 30, 0, 300);
	fHInvMassSSll      = new TH1D("InvMass_SSll"  ,  " Inv Mass of SS di-leptons ll"   , 30, 0, 300);
	
		
	for(int i=0; i<fMT2_histos_number; ++i){
		std::stringstream out;
		int mass = i*fMT2_histos_step;
		out << mass ;
		TString Mass = (TString) out.str();
			
		fHMT2_diBjet[i]                = new TH1D("MT2_diBjet_M"+Mass,            " ",100, 0., 500.);
		fHMT2_dijet[i]                 = new TH1D("MT2_dijet_M"+Mass,             " ",100, 0., 500.);
		fHMT2_pseudojet[i]             = new TH1D("MT2_pseudojet_M"+Mass,         " ",100, 0., 500.);
		fHMT2_diBjetdiLept[i]          = new TH1D("MT2_diBjetdiLept_M"+Mass,      " ",100, 0., 500.);
		fHMT2_PseudoJetsWithB[i]       = new TH1D("MT2_PseudoJetsWithB_M"+Mass,   " ",100, 0., 500.);     
		fHMT2_PseudoJetsWithLeptons[i] = new TH1D("MT2_PseudoJetsLeptons_M"+Mass, " ",100, 0., 500.);
		fHMT2_PseudoJetClean[i]        = new TH1D("MT2_PseudoJetsCleaned_M"+Mass, " ",100, 0., 500.);


		fHMT2_OSll[i]                 = new TH1D("MT2_OSll_M"+Mass,               " ", 40, 0., 400.);
		fHMT2_OSee[i]                 = new TH1D("MT2_OSee_M"+Mass,               " ", 40, 0., 400.);
		fHMT2_OSmumu[i]               = new TH1D("MT2_OSmumu_M"+Mass,             " ", 40, 0., 400.);
		fHMT2_OSemu[i]                = new TH1D("MT2_OSemu_M"+Mass,              " ", 40, 0., 400.);
		fHMT2_OSllminusemu[i]         = new TH1D("MT2_OSllminusemu_M"+Mass,       " ", 40, 0., 400.);
		
		fHMT2_SSll[i]                 = new TH1D("MT2_SSll_M"+Mass,               " ", 20, 0., 400.);
		
		fHMT2perp_OSee[i]             = new TH1D("MT2perp_OSee_M"+Mass,           " ", 40, 0., 400.);
		fHMT2perp_OSmumu[i]           = new TH1D("MT2perp_OSmumu_M"+Mass,         " ", 40, 0., 400.);
		fHMT2perp_OSemu[i]            = new TH1D("MT2perp_OSemu_M"+Mass,          " ", 40, 0., 400.);
	}
}

void MassAnalysis::Analyze(){	
	
	// ---------------------------------------------------
	// Initialize fElecs, fJets, fBJets, fMuons, fLeptConfig 
	InitializeEvent();
	// ----------------------------------------------------
	
	// --------------------------------------------------------------------
	// check if event passes selection cuts and trigger requirements
	if(! IsSelectedEvent()){return;}
	// --------------------------------------------------------------------

	// Mass distributions for OS dileptons
	DiLeptonMasses();

	// Masses for di-jets
	DiJetMasses();

	// Masses for di-bjets
	DiBJetMasses();

	// Masses for pseudojets
	PseudoJetMasses();

	// Masses for TTbar on GenLevel
	MassesForTTbar();

}

// *************************************************************************************************
void MassAnalysis::PseudoJetMasses(){
	if(fJets.size() <3 ) return;

	// ---------------------------------------------
	// MT2
	//    if there are selected leptons, add them to the visible system and fill fHMT2_PseudoJetsWithLeptons
	//    if there are no selected leptons and no jets above 30 that fail the selection cuts: fill fHMT2_PseudoJetClean
	//    if there is at least one b-jet: fill fHMT2_PseudoJetsWithB    
	//    in any case, fill fHMT2_pseudojet
	
	// make pseudojets with hemispheres
	vector<float> px, py, pz, E;
	for(int i=0; i<fJets.size(); ++i){
		px.push_back(fTR->JPx[fJets[i]]);
		py.push_back(fTR->JPy[fJets[i]]);
		pz.push_back(fTR->JPz[fJets[i]]);
		 E.push_back(fTR->JE[fJets[i]]);
	}
		
	vector<TLorentzVector> leptonmomenta = GetLepton4Momenta();
	for(int i=0; i< leptonmomenta.size(); ++i){
		px.push_back(leptonmomenta[i].Px());
		py.push_back(leptonmomenta[i].Py());
		pz.push_back(leptonmomenta[i].Pz());
		 E.push_back(leptonmomenta[i].E());
	}

	// get hemispheres (seed: max inv mass, association method: minimal lund distance)
	fHemisphere = new Hemisphere(px, py, pz, E, 2, 3);
	vector<int> grouping = fHemisphere->getGrouping();

	TLorentzVector pseudojet1=(0.,0.,0.,0.);
	TLorentzVector pseudojet2=(0.,0.,0.,0.);
	
	TLorentzVector pmiss;
	pmiss.SetPx(fTR->PFMETpx);
	pmiss.SetPy(fTR->PFMETpy);

	for(int i=0; i<px.size(); ++i){
		if(grouping[i]==1){
			pseudojet1.SetPx(pseudojet1.Px() + px[i]);
			pseudojet1.SetPy(pseudojet1.Py() + py[i]);
		}else if(grouping[i] == 2){
			pseudojet2.SetPx(pseudojet2.Px() + px[i]);
			pseudojet2.SetPy(pseudojet2.Py() + py[i]);
		}
	}
	delete fHemisphere;
		
	// count the number of jets > 30 GeV regardless of other acceptence and ID cuts
	int jcounter =0;
	for(int i=0; i<fTR->NJets; ++i){
		if(fTR->JPt[i]>30) jcounter++; 
	}
	bool jet_acceptance(true);
	if(jcounter != fJets.size()) jet_acceptance = false; 

	for(int i=0; i<fMT2_histos_number; ++i){
		double MT2=GetMT2(pseudojet1, 0., pseudojet2, 0., pmiss, i*fMT2_histos_step);
		fHMT2_pseudojet[i] -> Fill(MT2);
		if(fLeptConfig == null && jet_acceptance){fHMT2_PseudoJetClean[i]->Fill(MT2);}
		if(fLeptConfig!= null) {fHMT2_PseudoJetsWithLeptons[i] -> Fill(MT2);}
		if(fBJets.size()>0) {fHMT2_PseudoJetsWithB[i]->Fill(MT2);} 
	}



	// ----------------------------------------------
	// AlphaT does not make sense....
	std::vector<TLorentzVector> jets;
	jcounter =0;
	for(int i=0; i<fTR->NJets; ++i){
		if(fTR->JPt[i]>30) jcounter++;
	}
	// requiring no leptons, and no jets above 30 that failed the selection
	if(fLeptConfig == null && fJets.size()!=0 && fTR->JPt[fJets[0]] > 150 && jcounter == fJets.size() ){
		for(int i=0; i<fJets.size(); ++i){
			TLorentzVector v;
	 		v.SetPxPyPzE(fTR->JPx[fJets[i]],fTR->JPy[fJets[i]], fTR->JPz[fJets[i]], fTR->JE[fJets[i]]);
	 		jets.push_back(v);
	 	}
	 	double alphaT= GetAlphaT(jets);
	 	if(fJets.size()>2){  // should be jets.size()
	 		fHAlpahT_PseudoJet->Fill(alphaT);
	 	}
	}

}

// **************************************************************************************************
void MassAnalysis::DiBJetMasses(){
	if(fBJets.size()!=2) return;
		
	TLorentzVector p1, p2;
	p1.SetPxPyPzE(fTR->JPx[fBJets[0]], fTR->JPy[fBJets[0]], fTR->JPz[fBJets[0]], fTR->JE[fBJets[0]]);
	p2.SetPxPyPzE(fTR->JPx[fBJets[1]], fTR->JPy[fBJets[1]], fTR->JPz[fBJets[1]], fTR->JE[fBJets[1]]);

	// -------------------------------------------
	// MCT
	double MCT=GetMCT(p1, p2);
	fHMCT_DiBJet -> Fill(MCT);

	// ------------------------------------------
	// MT2
	TLorentzVector pmiss(0., 0., 0., 0.);
	pmiss.SetPx(fTR->PFMETpx);
	pmiss.SetPy(fTR->PFMETpy);

	for(int i=0; i<fMT2_histos_number; ++i){
		double MT2=GetMT2(p1, 0., p2, 0., pmiss, i*fMT2_histos_step);
		fHMT2_diBjet[i] ->Fill(MT2);
	}
	
	// MT2 in case of additional leptons
	// di-b-jets with leptons summed up to MET
  	// for TTbar events, the parent mass would then be the top!
	if(fLeptConfig == OS_ee || fLeptConfig == OS_emu ||fLeptConfig == OS_mumu ){
		vector<TLorentzVector> leptmomenta = GetLepton4Momenta();
		TLorentzVector p_invisible=pmiss + leptmomenta[0] + leptmomenta[1];
		for(int i=0; i<fMT2_histos_number; ++i){	
			double MT2=GetMT2(p1, 0., p2, 0., p_invisible, i*fMT2_histos_step);
			fHMT2_diBjetdiLept[i] -> Fill(MT2);
		}
	}
}





// ***************************************************************************************************
void MassAnalysis::DiJetMasses(){
	if( fJets.size()!=2 ) return;
	TLorentzVector p1;
	TLorentzVector p2;		
	p1.SetPxPyPzE(fTR->JPx[fJets[0]], fTR->JPy[fJets[0]], fTR->JPz[fJets[0]], fTR->JE[fJets[0]]);
	p2.SetPxPyPzE(fTR->JPx[fJets[1]], fTR->JPy[fJets[1]], fTR->JPz[fJets[1]], fTR->JE[fJets[1]]);
		
	// ------------------------------------------------
	// MCT
	double MCT=GetMCT(p1, p2);
	fHMCT_DiJet -> Fill(MCT);
	
	// -------------------------------------------------
	// MT2
	TLorentzVector pmiss(0., 0., 0., 0.);
	pmiss.SetPx(fTR->PFMETpx);
	pmiss.SetPy(fTR->PFMETpy);
	for(int i=0; i<fMT2_histos_number; ++i){
		double MT2=GetMT2(p1, 0., p2, 0., pmiss, i*fMT2_histos_step);
		fHMT2_dijet[i] ->Fill(MT2);
	}

	// --------------------------------------------------
	// alphaT : this implementation does not make sense...
	std::vector<TLorentzVector> jets;
	int jcounter =0;
	for(int i=0; i<fTR->NJets; ++i){
		        if(fTR->JPt[i]>30) jcounter++;
	}
	// requiring no leptons, and no jets above 30 that failed the selection
	if(fLeptConfig == null && fJets.size()!=0 && fTR->JPt[fJets[0]] > 150 && jcounter == fJets.size() ){
		for(int i=0; i<fJets.size(); ++i){
			TLorentzVector v;
	 		v.SetPxPyPzE(fTR->JPx[fJets[i]],fTR->JPy[fJets[i]], fTR->JPz[fJets[i]], fTR->JE[fJets[i]]);
	 		jets.push_back(v);
	 	}
	 	double alphaT= GetAlphaT(jets);
	 	if(fJets.size()==2){  // should be jets.size()
	 		fHAlpahT_DiJet->Fill(alphaT);
	 	}
	}
}

// ****************************************************************************************************
void MassAnalysis::DiLeptonMasses(){
	if(! (fLeptConfig==OS_ee || fLeptConfig==OS_mumu || fLeptConfig==OS_emu || fLeptConfig==SS_ee || fLeptConfig==SS_mumu || fLeptConfig==SS_emu) ) return;

	// TLorentzVector of di-leptons
	vector<TLorentzVector> dilept4mom = GetLepton4Momenta();
	TLorentzVector p1=dilept4mom[0];
	TLorentzVector p2=dilept4mom[1];

	// define Upstream Transverse Momentum
	TLorentzVector P_UTM(0.,0.,0.,0.);
	for(int i=0; i<fJets.size(); ++i){
		TLorentzVector p;
		p.SetPxPyPzE(fTR->JPx[fJets[i]],fTR->JPy[fJets[i]],fTR->JPz[fJets[i]],fTR->JE[fJets[i]]);
		P_UTM +=p;
	}

	// calculate inv mass, MCT, MCTperp, MT2, MT2perp	
	double MCTperp=-1.;
	double MT2perp[fMT2_histos_number];
	double MT2[fMT2_histos_number];

	double MCT    =GetMCT(p1, p2);
	double mass   =(p1+p2).M();
	
	if(P_UTM != (0.,0.,0.,0.)) {
		MCTperp=GetMCTperp(p1, p2, P_UTM);
	}
	TLorentzVector pmiss(0., 0., 0., 0.);
	pmiss.SetPx(fTR->PFMETpx);
	pmiss.SetPy(fTR->PFMETpy);
	for(int i=0; i<fMT2_histos_number; ++i){
		int m_invisible = i*fMT2_histos_step;
		MT2[i]    =GetMT2(p1, p1.M(), p2, p2.M(), pmiss, m_invisible);
		if(P_UTM != (0.,0.,0.,0.)){
			MT2perp[i]=GetMT2perp(p1, p2, P_UTM, m_invisible);
		}else   MT2perp[i]=-1.;
	}
		
		
	// fill histos
	if(fLeptConfig==OS_ee){
		fHInvMassOSee   -> Fill(mass);
		fHInvMassOSll   -> Fill(mass);	
		if(P_UTM != (0.,0.,0.,0.)) fHMCTperp_OSee -> Fill(MCTperp);
		fHMCT_OSee      -> Fill(MCT);
		for(int i=0; i<fMT2_histos_number; ++i){
			fHMT2_OSll[i]     -> Fill(MT2[i]);
			fHMT2_OSee[i]     -> Fill(MT2[i]);
			if(P_UTM != (0.,0.,0.,0.)) {
				fHMT2perp_OSee[i] -> Fill(MT2perp[i]);
			}
		}
	} else if(fLeptConfig==OS_mumu){
		fHInvMassOSmumu -> Fill(mass);
		fHInvMassOSll   -> Fill(mass);	
		if(P_UTM != (0.,0.,0.,0.)) fHMCTperp_OSmumu -> Fill(MCTperp);
		fHMCT_OSmumu      -> Fill(MCT);
		for(int i=0; i<fMT2_histos_number; ++i){
			fHMT2_OSll[i]       -> Fill(MT2[i]);
			fHMT2_OSmumu[i]     -> Fill(MT2[i]);
			if(P_UTM != (0.,0.,0.,0.)) {
				fHMT2perp_OSmumu[i] -> Fill(MT2perp[i]);
			}
		}
	} else if(fLeptConfig==OS_emu){
		fHInvMassOSemu  -> Fill(mass);
		fHInvMassOSll   -> Fill(mass);	
		if(P_UTM != (0.,0.,0.,0.)) fHMCTperp_OSemu -> Fill(MCTperp);
		fHMCT_OSemu     -> Fill(MCT);
		for(int i=0; i<fMT2_histos_number; ++i){
			fHMT2_OSll[i]       -> Fill(MT2[i]);
			fHMT2_OSemu[i]      -> Fill(MT2[i]);
			if(P_UTM != (0.,0.,0.,0.)) {
				fHMT2perp_OSemu[i]  -> Fill(MT2perp[i]);
			}
		}
	} else if(fLeptConfig==SS_emu || fLeptConfig==SS_ee || fLeptConfig==SS_mumu ){
		for(int i=0; i<fMT2_histos_number; ++i){
			fHMT2_SSll[i]    -> Fill(MT2[i]);
		}
		fHInvMassSSll -> Fill(mass);
		if(fLeptConfig==SS_emu ) fHInvMassSSemu  -> Fill(mass);
		if(fLeptConfig==SS_ee  ) fHInvMassSSee   -> Fill(mass);
		if(fLeptConfig==SS_mumu) fHInvMassSSmumu -> Fill(mass);
	}
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

	return MCT_perp;
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
void MassAnalysis::MassesForTTbar(){
	// MT2 and MCTperp for TTbar events using Gen Info
	vector<TLorentzVector> leptons;
	for(int i=0; i<fTR->NGenLeptons; ++i){
		if(! (abs(fTR->GenLeptonID[i]) == 11 || abs(fTR->GenLeptonID[i]) == 13) ) continue;
		if(! (abs(fTR->GenLeptonMID[i]) == 24) ) continue;
		if(! (abs(fTR->GenLeptonGMID[i]) == 6) ) continue;
		TLorentzVector p;
		p.SetPtEtaPhiM(fTR->GenLeptonPt[i],fTR->GenLeptonEta[i],fTR->GenLeptonPhi[i],0.);
		leptons.push_back(p);
	}
	vector<TLorentzVector> neutrinos;
	for(int i=0; i<fTR->NGenLeptons; ++i){
		if(! (abs(fTR->GenLeptonID[i]) == 12 || abs(fTR->GenLeptonID[i]) == 14) ) continue;
		if(! (abs(fTR->GenLeptonMID[i]) == 24) ) continue;
		if(! (abs(fTR->GenLeptonGMID[i]) == 6) ) continue;
		TLorentzVector p;
		p.SetPtEtaPhiM(fTR->GenLeptonPt[i],fTR->GenLeptonEta[i],fTR->GenLeptonPhi[i],0.);
		neutrinos.push_back(p);
	}
	if(leptons.size()==2 && neutrinos.size()==2){
		// define P_UTM: negative vector sum of leptons and neutrinos
		//               this should be equal to the two b-jets (no quarks in genparticles)
		double Px_UTM = -(leptons[0].Px()+leptons[1].Px()+neutrinos[0].Px()+neutrinos[1].Px() );
		double Py_UTM = -(leptons[0].Py()+leptons[1].Py()+neutrinos[0].Py()+neutrinos[1].Py() );
		TLorentzVector P_UTM(0., 0., 0., 0.);
		P_UTM.SetPx(Px_UTM);
		P_UTM.SetPy(Py_UTM);
		
		// MCTperp, MT2perp, MCT, MT2
		TLorentzVector pmiss = neutrinos[0]+neutrinos[1];
		double MT2    =GetMT2(leptons[0], leptons[0].M(), leptons[1], leptons[1].M(), pmiss, 0.);
		double MCT    =GetMCT(leptons[0], leptons[1]);
		double MCTperp=GetMCTperp(leptons[0], leptons[1], P_UTM);
		double MT2perp=GetMT2perp(leptons[0], leptons[1], P_UTM, 0.);
		fHMCT_TTbar     -> Fill(MCT);
		fHMCTperp_TTbar -> Fill(MCTperp);
		fHMT2perp_TTbar -> Fill(MT2perp);
		fHMT2_TTbar     -> Fill(MT2);
	}
}

// ****************************************************************************************************
double MassAnalysis::GetMCT(TLorentzVector p1, TLorentzVector p2){
	double m1=p1.M();
	double m2=p2.M();
	double ET1=sqrt( pow(p1.Pt(),2) +pow(m1,2) );
	double ET2=sqrt( pow(p2.Pt(),2) +pow(m2,2) );
	
	double MCTsquared=pow(m1,2) + pow(m2,2) +2*(ET1*ET2+ (p1.Px() * p2.Px()) + (p1.Py() * p2.Py()));
	
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
  	for(unsigned i=0; i < ht.size(); i++) {
    	for(unsigned j=0; j < p4s.size(); j++) {
			ht [i] [(i/(1<<j))%2] += p4s[j].Pt();
    	}
  	}
  	std::vector<double> deltaHT; for(unsigned i=0; i<ht.size(); i++) deltaHT.push_back(fabs(ht[i][0]-ht[i][1]));
  	return deltaHT;
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
		p.SetPtEtaPhiE(fTR->ElPt[fElecs[i]],fTR->ElEta[fElecs[i]], fTR->ElPhi[fElecs[i]], fTR->ElE[fElecs[i]]);
		momenta.push_back(p);
	}
	for(int i=0; i<fMuons.size(); ++i){
		TLorentzVector p;
		p.SetPtEtaPhiM(fTR->MuPt[fMuons[i]],fTR->MuEta[fMuons[i]], fTR->MuPhi[fMuons[i]],0.105);
		momenta.push_back(p);
	}
	return momenta;
}

// ****************************************************************************************************
void MassAnalysis::End(){
	fHistFile->cd();	
	
	// fill OS ll minus emu MT2 histo
	for(int i=0; i<fMT2_histos_number; ++i){
		for(int bin=1; bin<= fHMT2_OSee[0]->GetNbinsX(); ++bin){
			double bincont = fHMT2_OSee[i]->GetBinContent(bin)+fHMT2_OSmumu[i]->GetBinContent(bin);
			fHMT2_OSllminusemu[i] ->SetBinContent(bin, bincont);
		}
	}
	for(int i=0; i<fMT2_histos_number; ++i){
		fHMT2_OSllminusemu[i]->Write();
		fHMT2_OSll[i]        ->Write();
		fHMT2_OSee[i]        ->Write();
		fHMT2_OSmumu[i]      ->Write();
		fHMT2_OSemu[i]       ->Write();
		
		fHMT2_SSll[i]        ->Write();
	
		fHMT2_pseudojet[i]             ->Write();
		fHMT2_PseudoJetsWithB[i]       ->Write();
		fHMT2_PseudoJetsWithLeptons[i] ->Write();
		fHMT2_PseudoJetClean[i]        ->Write();
	        fHMT2_diBjetdiLept[i]->Write();
		fHMT2_dijet[i]       ->Write();
		fHMT2_diBjet[i]      ->Write();
		
		fHMT2perp_OSee[i]    ->Write();
		fHMT2perp_OSmumu[i]  ->Write();
		fHMT2perp_OSemu[i]   ->Write();
	}
	fHMT2_vs_M_OSDiLept      ->Write();
	
	fHAlpahT_DiJet           ->Write();
	fHAlpahT_PseudoJet       ->Write();
 
	
	fHMCT_DiJet              ->Write();
	fHMCT_DiBJet             ->Write();
	fHMCT_OSee               ->Write();
	fHMCT_OSmumu             ->Write();
	fHMCT_OSemu              ->Write();
	 
	fHMCTperp_OSee           ->Write();
	fHMCTperp_OSmumu         ->Write();
	fHMCTperp_OSemu          ->Write();
	
	fHMCTperp_TTbar          ->Write();
	fHMT2_TTbar              ->Write();
	fHMT2perp_TTbar          ->Write();
	fHMCT_TTbar              ->Write();
	
	fHInvMassOSee            ->Write();
	fHInvMassOSmumu          ->Write();
	fHInvMassOSemu           ->Write();
	fHInvMassOSll            ->Write();
	fHInvMassSSee            ->Write();
	fHInvMassSSmumu          ->Write();
	fHInvMassSSemu           ->Write();
	fHInvMassSSll            ->Write();
	
	fHistFile                ->Close();
}
