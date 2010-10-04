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
		
		TString hist_name_dijet           = "MT2_dijet_M" + (TString) out.str();	
		TString hist_name_diBjet          = "MT2_diBjet_M" + (TString) out.str();
		TString hist_name_pseudojet       = "MT2_pseudojet_M" + (TString) out.str();		
		TString hist_name_diBjetdiLept    = "MT2_diBjetdiLept_M" + (TString) out.str();
		TString hist_name_PseudoJetsWithB = "MT2_PseudoJetsWithB_M" + (TString) out.str();
		TString hist_name_PseudoJetsWithLeptons = "MT2_PseudoJetsLeptons_M" + (TString) out.str();
		TString hist_name_PseudoJetsCleaned = "MT2_PseudoJetsCleaned_M" + (TString) out.str();
	
		TString hist_name_SSll            = "MT2_SSll_M" + (TString) out.str();		
		TString hist_name_OSll            = "MT2_OSll_M" + (TString) out.str();
		TString hist_name_OSllminusemu    = "MT2_OSllminusemu_M" + (TString) out.str();
		TString hist_name_OSee            = "MT2_OSee_M" + (TString) out.str();
		TString hist_name_OSmumu          = "MT2_OSmumu_M" + (TString) out.str();
		TString hist_name_OSemu           = "MT2_OSemu_M" + (TString) out.str();	
		
		TString hist_name_MT2perp_OSee   = "MT2perp_OSee_M" + (TString) out.str();
		TString hist_name_MT2perp_OSmumu = "MT2perp_OSmumu_M" + (TString) out.str();
		TString hist_name_MT2perp_OSemu  = "MT2perp_OSemu_M" + (TString) out.str();

		TString title_dijet           = "Di-Jet MT2, mass = " + (TString) out.str();
		TString title_diBjet          = "Di-bJet MT2, mass = " + (TString) out.str();
		TString title_diBjetdiLept    = "Di-bJet with di-Lepton MT2, mass = " + (TString) out.str();
		TString title_pseudojet       = "Pseudo-Jet MT2, mass = " + (TString) out.str();
		TString title_PseudoJetsWithB = "Pseudo-Jet MT2 with at least one b jet, mass = " +(TString) out.str();
		TString title_PseudoJetsWithLeptons = "Pseudo-Jet MT2 with at least one lepton, mass = " +(TString) out.str();
		TString title_PseudoJetsClean = "Pseudo-Jet MT2 cleaned, mass = " +(TString) out.str();		
		
		TString title_SSll           = "SS dilept MT2, mass = " + (TString) out.str();
		TString title_OSll           = "OS ll MT2, mass = " + (TString) out.str();
		TString title_OSee           = "OS ee MT2, mass = " + (TString) out.str();
		TString title_OSmumu         = "OS mumu MT2, mass = " + (TString) out.str();
		TString title_OSemu          = "OS emu MT2, mass = " + (TString) out.str();
		TString title_OSllminusemu   = "OS ll minus emu MT2, mass = " + (TString) out.str();
		
		TString title_MT2perp_OSee     ="MT2perp OS ee, mass = " + (TString) out.str();
		TString title_MT2perp_OSemu    ="MT2perp OS emu, mass = " + (TString) out.str();
		TString title_MT2perp_OSmumu   ="MT2perp OS mumu, mass = " + (TString) out.str();
		

		fHMT2_diBjet[i]          = new TH1D(hist_name_diBjet        , title_diBjet        , 100, 0., 500.);
		fHMT2_dijet[i]           = new TH1D(hist_name_dijet         , title_dijet         , 100, 0., 500.);
		fHMT2_pseudojet[i]       = new TH1D(hist_name_pseudojet     , title_pseudojet     , 100, 0., 500.);
		fHMT2_diBjetdiLept[i]    = new TH1D(hist_name_diBjetdiLept  , title_diBjetdiLept  , 100, 0., 500.);
		fHMT2_PseudoJetsWithB[i]       = new TH1D(hist_name_PseudoJetsWithB        , title_PseudoJetsWithB        , 100, 0., 500.);     
		fHMT2_PseudoJetsWithLeptons[i] = new TH1D(hist_name_PseudoJetsWithLeptons  , title_PseudoJetsWithLeptons  , 100, 0., 500.);
		fHMT2_PseudoJetClean[i]        = new TH1D(hist_name_PseudoJetsCleaned      , title_PseudoJetsClean        , 100, 0., 500.);


		fHMT2_OSll[i]            = new TH1D(hist_name_OSll          , title_OSll          , 40, 0., 400.);
		fHMT2_OSee[i]            = new TH1D(hist_name_OSee          , title_OSee          , 40, 0., 400.);
		fHMT2_OSmumu[i]          = new TH1D(hist_name_OSmumu        , title_OSmumu        , 40, 0., 400.);
		fHMT2_OSemu[i]           = new TH1D(hist_name_OSemu         , title_OSemu         , 40, 0., 400.);
		fHMT2_OSllminusemu[i]    = new TH1D(hist_name_OSllminusemu  , title_OSllminusemu  , 40, 0., 400.);
		
		fHMT2_SSll[i]            = new TH1D(hist_name_SSll          , title_SSll          ,  20, 0., 400.);
		
		fHMT2perp_OSee[i]        = new TH1D(hist_name_MT2perp_OSee    , title_MT2perp_OSee    , 40, 0., 400.);
		fHMT2perp_OSmumu[i]      = new TH1D(hist_name_MT2perp_OSmumu  , title_MT2perp_OSmumu  , 40, 0., 400.);
		fHMT2perp_OSemu[i]       = new TH1D(hist_name_MT2perp_OSemu   , title_MT2perp_OSemu   , 40, 0., 400.);
		
	}
	
}

void MassAnalysis::Analyze(){	
	
	// ---------------------------------------------------
	// Initialize fElecs, fJets, fBJets, fMuons, fLeptConfig 
	InitializeEvent();
	// ----------------------------------------------------
		
	// --------------------------------------------------------------------
	// check if event passes selection cuts and trigger requirements
	if(! IsGoodEvent()){return;}
	// --------------------------------------------------------------------

	// ----------------------------------------------------------
	// MT2 plots
	DiLeptonMT2();
	JetMT2();
	
	
	// --------------------------------------------------------
	// AlphaT !!veto leptonic events!!
	std::vector<TLorentzVector> jets;
	TLorentzVector v;
	int jcounter =0;
	for(int i=0; i<fTR->NJets; ++i){
		if(fTR->JPt[i]>30) jcounter++; 
	}
	// requiring no leptons, and no jets above 30 that failed the selection
	if(fLeptConfig == null && fJets.size()!=0 && fTR->JPt[fJets[0]] > 150 && jcounter == fJets.size() ){
		for(int i=0; i<fJets.size(); ++i){
			v.SetPxPyPzE(fTR->JPx[fJets[i]],fTR->JPy[fJets[i]], fTR->JPz[fJets[i]], fTR->JE[fJets[i]]);	
			jets.push_back(v);
		}		

		double alphaT= GetAlphaT(jets);
		if(fJets.size()==2){
			fHAlpahT_DiJet->Fill(alphaT);
		} else if (fJets.size() >2){
			fHAlpahT_PseudoJet->Fill(alphaT);
		}
	}
	
	// -------------------------------------------------------------
	// MCT, MCTperp, MT2perp ---------------------------------------
	TLorentzVector p1;
	TLorentzVector p2;
	if(fJets.size()==2){
		p1.SetPxPyPzE(fTR->JPx[fJets[0]], fTR->JPy[fJets[0]], fTR->JPz[fJets[0]], fTR->JE[fJets[0]]);
		p2.SetPxPyPzE(fTR->JPx[fJets[1]], fTR->JPy[fJets[1]], fTR->JPz[fJets[1]], fTR->JE[fJets[1]]);
		double MCT=GetMCT(p1, p2);
		fHMCT_DiJet -> Fill(MCT);
	}
	if(fBJets.size()==2){
		p1.SetPxPyPzE(fTR->JPx[fBJets[0]], fTR->JPy[fBJets[0]], fTR->JPz[fBJets[0]], fTR->JE[fBJets[0]]);
		p2.SetPxPyPzE(fTR->JPx[fBJets[1]], fTR->JPy[fBJets[1]], fTR->JPz[fBJets[1]], fTR->JE[fBJets[1]]);
		double MCT=GetMCT(p1, p2);
		fHMCT_DiBJet -> Fill(MCT);
	}
	if(fLeptConfig==OS_ee || fLeptConfig==OS_mumu || fLeptConfig==OS_emu){
		
		// define Upstream Transverse Momentum
		TLorentzVector P_UTM(0.,0.,0.,0.);
		for(int i=0; i<fJets.size(); ++i){
			TLorentzVector p;
			p.SetPxPyPzE(fTR->JPx[fJets[i]],fTR->JPy[fJets[i]],fTR->JPz[fJets[i]],fTR->JE[fJets[i]]);
			P_UTM +=p;
		}
		// fill 4momenta
		if(fLeptConfig==OS_ee){
			p1.SetPtEtaPhiE(fTR->ElPt[fElecs[0]], fTR->ElEta[fElecs[0]], fTR->ElPhi[fElecs[0]], fTR->ElE[fElecs[0]]);
			p2.SetPtEtaPhiE(fTR->ElPt[fElecs[1]], fTR->ElEta[fElecs[1]], fTR->ElPhi[fElecs[1]], fTR->ElE[fElecs[1]]);
		} else if(fLeptConfig==OS_mumu){
			p1.SetPtEtaPhiM(fTR->MuPt[fMuons[0]], fTR->MuEta[fMuons[0]], fTR->MuPhi[fMuons[0]], 0.105);
			p2.SetPtEtaPhiM(fTR->MuPt[fMuons[1]], fTR->MuEta[fMuons[1]], fTR->MuPhi[fMuons[1]], 0.105);
		} else if(fLeptConfig==OS_emu){
			p1.SetPtEtaPhiE(fTR->ElPt[fElecs[0]], fTR->ElEta[fElecs[0]], fTR->ElPhi[fElecs[0]], fTR->ElE[fElecs[0]]);
			p2.SetPtEtaPhiM(fTR->MuPt[fMuons[0]], fTR->MuEta[fMuons[0]], fTR->MuPhi[fMuons[0]], 0.105);
		}
		// fill histos
		if(fLeptConfig==OS_ee){
			double MCT=GetMCT(p1, p2);
			fHMCT_OSee     -> Fill(MCT);	
			if(P_UTM != (0.,0.,0.,0.)){
				double MCTperp=GetMCTperp(p1, p2, P_UTM);
				fHMCTperp_OSee -> Fill(MCTperp);
			
				for(int i=0; i<fMT2_histos_number; ++i){
					double m_inv =i*fMT2_histos_step;
					double MT2perp=GetMT2perp(p1, p2, P_UTM, m_inv);
					fHMT2perp_OSee[i] -> Fill(MT2perp);
				}
			}
		} else if(fLeptConfig==OS_mumu){
			double MCT=GetMCT(p1, p2);
			fHMCT_OSmumu     -> Fill(MCT);	
			if(P_UTM != (0.,0.,0.,0.)){
				double MCTperp=GetMCTperp(p1, p2, P_UTM);
				fHMCTperp_OSmumu -> Fill(MCTperp);

				for(int i=0; i<fMT2_histos_number; ++i){
					double m_inv =i*fMT2_histos_step;
					double MT2perp=GetMT2perp(p1, p2, P_UTM, m_inv);
					fHMT2perp_OSmumu[i] -> Fill(MT2perp);
				}
			}
		} else if(fLeptConfig==OS_emu){
			double MCT=GetMCT(p1, p2);
			fHMCT_OSemu     -> Fill(MCT);	
			if(P_UTM != (0.,0.,0.,0.)){
				double MCTperp=GetMCTperp(p1, p2, P_UTM);
				fHMCTperp_OSemu -> Fill(MCTperp);
			
				for(int i=0; i<fMT2_histos_number; ++i){
					double m_inv =i*fMT2_histos_step;
					double MT2perp=GetMT2perp(p1, p2, P_UTM, m_inv);
					fHMT2perp_OSemu[i] -> Fill(MT2perp);
				}
			}	
		}
		
	}
	
	// ----------------------------------------------------------------------------------
	// MCTperp and MT2 for TTbar events
	MassesForTTbar();
	
	// ----------------------------------------------------------------------------------
	// Di-Lepton Inv Masses
	TLorentzVector l1, l2;
	if(fLeptConfig==OS_ee || fLeptConfig==SS_ee){
		l1.SetPtEtaPhiE(fTR->ElPt[fElecs[0]],fTR->ElEta[fElecs[0]],fTR->ElPhi[fElecs[0]],fTR->ElE[fElecs[0]]);
		l2.SetPtEtaPhiE(fTR->ElPt[fElecs[1]],fTR->ElEta[fElecs[1]],fTR->ElPhi[fElecs[1]],fTR->ElE[fElecs[1]]);
		double mass=(l1+l2).M();
		if(fLeptConfig==OS_ee) {
			fHInvMassOSee -> Fill(mass);
			fHInvMassOSll -> Fill(mass);
		} else  {
			fHInvMassSSee -> Fill(mass);
			fHInvMassSSll  -> Fill(mass);
		}
	} else if (fLeptConfig==OS_mumu || fLeptConfig==SS_mumu){
		l1.SetPtEtaPhiM(fTR->MuPt[fMuons[0]],fTR->MuEta[fMuons[0]],fTR->MuPhi[fMuons[0]],0.105);
		l2.SetPtEtaPhiM(fTR->MuPt[fMuons[1]],fTR->MuEta[fMuons[1]],fTR->MuPhi[fMuons[1]],0.105);
		double mass=(l1+l2).M();
		if(fLeptConfig==OS_mumu) {
			fHInvMassOSmumu -> Fill(mass);
			fHInvMassOSll   -> Fill(mass);
		} else {
			fHInvMassSSmumu -> Fill(mass);
			fHInvMassSSll  -> Fill(mass);
		}
	} else if (fLeptConfig==OS_emu || fLeptConfig==SS_emu){
		l1.SetPtEtaPhiM(fTR->MuPt[fMuons[0]],fTR->MuEta[fMuons[0]],fTR->MuPhi[fMuons[0]],0.105);
		l2.SetPtEtaPhiE(fTR->ElPt[fElecs[0]],fTR->ElEta[fElecs[0]],fTR->ElPhi[fElecs[0]],fTR->ElE[fElecs[0]]);
		double mass=(l1+l2).M();
		if(fLeptConfig==OS_emu)  {
			fHInvMassOSemu -> Fill(mass);
			fHInvMassOSll  -> Fill(mass);
		} else {
			fHInvMassSSemu -> Fill(mass);
			fHInvMassSSll  -> Fill(mass);
		}
	}	
}


 // ****************************************************************************************************
double MassAnalysis::GetMCTperp(TLorentzVector p1, TLorentzVector p2, TLorentzVector P_UTM){
	double m1=p1.M();
	double m2=p2.M();
	
	TVector3 p1_T_perp  =  GetMomPerp(p1, P_UTM); 
	TVector3 p2_T_perp  =  GetMomPerp(p2, P_UTM);
	
	double E1_T_perp = pow( pow(m1,2) + p1_T_perp.Mag2() ,0.5);
	double E2_T_perp = pow( pow(m2,2) + p2_T_perp.Mag2() ,0.5);	
	
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
		
		// MCTperp, MT2perp, MCT
		double MCT    =GetMCT(leptons[0], leptons[1]);
		double MCTperp=GetMCTperp(leptons[0], leptons[1], P_UTM);
		double MT2perp=GetMT2perp(leptons[0], leptons[1], P_UTM, 0);
		fHMCT_TTbar     -> Fill(MCT);
		fHMCTperp_TTbar -> Fill(MCTperp);
		fHMT2perp_TTbar -> Fill(MT2perp);
		
		// MT2
		double pa[3], pb[3], pmiss[3];
		pa[0]=0;  // corresponds to mass of particle a
		pa[1]=leptons[0].Px();  // px of particle a
		pa[2]=leptons[0].Py();  // py of particle a
		pb[0]=0;  // corresponds to mass of particle b
		pb[1]=leptons[1].Px();  // px of particle b
		pb[2]=leptons[1].Py();  // py of particle b
		pmiss[0]=0; // this value is ignored in Davis code
		pmiss[1]=neutrinos[0].Px() + neutrinos[1].Px();
		pmiss[2]=neutrinos[0].Py() + neutrinos[1].Py();
		double MT2=GetMT2(pa, pb, pmiss, 0);
		fHMT2_TTbar->Fill(MT2);
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
double MassAnalysis::GetMT2(double pa[], double pb[], double pmiss[], int m_invisible){
	fMT2 = new Davismt2();
	fMT2->set_momenta(pa, pb, pmiss);
	fMT2->set_mn(m_invisible);
	double MT2=fMT2->get_mt2();
	delete fMT2;
	return MT2;
}

// ****************************************************************************************************
void MassAnalysis::JetMT2(){
	double pa[3], pb[3], pmiss[3];
	pa[0]=0;  // corresponds to mass of particle a
	pa[1]=0;  // px of particle a
	pa[2]=0;  // py of particle a
	pb[0]=0;  // corresponds to mass of particle b
	pb[1]=0;  // px of particle b
	pb[2]=0;  // py of particle b
	pmiss[0]=0; // this value is ignored in Davis code
	pmiss[1]=fTR->PFMETpx;
	pmiss[2]=fTR->PFMETpy;
	
	// di-jet
	if(fJets.size()==2){
		pa[1]=fTR->JPx[fJets[0]];
		pa[2]=fTR->JPy[fJets[0]];
		pb[1]=fTR->JPx[fJets[1]];
		pb[2]=fTR->JPy[fJets[1]];
		
		for(int i=0; i<fMT2_histos_number; ++i){
			double MT2=GetMT2(pa, pb, pmiss, i*fMT2_histos_step);
			fHMT2_dijet[i] ->Fill(MT2);
		}
	}
	
	// di B-jets	
	if (fBJets.size()==2 ){
		// inclusive di B-jets without caring about leptons
		pa[1]=fTR->JPx[fBJets[0]];
		pa[2]=fTR->JPy[fBJets[0]];
		pb[1]=fTR->JPx[fBJets[1]];
		pb[2]=fTR->JPy[fBJets[1]];
		
		for(int i=0; i<fMT2_histos_number; ++i){
			double MT2=GetMT2(pa, pb, pmiss, i*fMT2_histos_step);
			fHMT2_diBjet[i] ->Fill(MT2);
		}
		
		// di-b-jets with leptons summed up to MET
  		// for TTbar events, the parent mass would then be the top!
		if(fLeptConfig == OS_ee || fLeptConfig == OS_emu ||fLeptConfig == OS_mumu ){
			vector<TLorentzVector> momenta = GetLepton4Momenta();
			double METplusLept[3];
			METplusLept[0] = 0;
			METplusLept[1] = pmiss[1]+ momenta[0].Px() + momenta[1].Px();
			METplusLept[2] = pmiss[2]+ momenta[0].Py() + momenta[1].Py();
			for(int i=0; i<fMT2_histos_number; ++i){	
				double MT2=GetMT2(pa, pb, METplusLept, i*fMT2_histos_step);
				fHMT2_diBjetdiLept[i] -> Fill(MT2);
			}
		}			
	}
	
	// pseudojets
	//    if there are selected leptons, add them to the visible system and fill fHMT2_PseudoJetsWithLeptons
	//    if there are no selected leptons and no jets above 30 that fail the selection cuts: fill fHMT2_PseudoJetClean
	//    if there is at least one b-jet: fill fHMT2_PseudoJetsWithB    
	//    in any case, fill fHMT2_pseudojet
 	if (fJets.size()>2){
		pa[1]=0;  // px of particle a
		pa[2]=0;  // py of particle a
		pb[1]=0;  // px of particle b
		pb[2]=0;  // py of particle b
	
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

		for(int i=0; i<px.size(); ++i){
			if(grouping[i]==1){
				pa[1] +=px[i];
				pa[2] +=py[i];
			}else if(grouping[i] == 2){
				pb[1] +=px[i];
				pb[2] +=py[i];				
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
			double MT2=GetMT2(pa, pb, pmiss, i*fMT2_histos_step);
			fHMT2_pseudojet[i] -> Fill(MT2);
			if(fLeptConfig == null && jet_acceptance){fHMT2_PseudoJetClean[i]->Fill(MT2);}
			if(fLeptConfig!= null) {fHMT2_PseudoJetsWithLeptons[i] -> Fill(MT2);}
			if(fBJets.size()>0) {fHMT2_PseudoJetsWithB[i]->Fill(MT2);} 
		}
	} 
}

// ****************************************************************************************************
void MassAnalysis::DiLeptonMT2(){
	double pa[3], pb[3], pmiss[3];
	pa[0]=0;  // corresponds to mass of particle a
	pb[0]=0;  // corresponds to mass of particle b
	pmiss[0]=0; // this value is ignored in Davis code
	pmiss[1]=fTR->PFMETpx;
	pmiss[2]=fTR->PFMETpy;
	
	if(fLeptConfig == null || fLeptConfig ==  e || fLeptConfig == mu){return;}
	else if(fLeptConfig==OS_ee || fLeptConfig == SS_ee){
		pa[1]=fTR->ElPx[fElecs[0]];
		pa[2]=fTR->ElPy[fElecs[0]];
		pb[1]=fTR->ElPx[fElecs[1]];
		pb[2]=fTR->ElPy[fElecs[1]];
	}else if(fLeptConfig==OS_mumu || fLeptConfig == SS_mumu){
		pa[1]=fTR->MuPx[fMuons[0]];
		pa[2]=fTR->MuPy[fMuons[0]];
		pb[1]=fTR->MuPx[fMuons[1]];
		pb[2]=fTR->MuPy[fMuons[1]];
	}else if(fLeptConfig==OS_emu || fLeptConfig == SS_emu){
		pa[1]=fTR->ElPx[fElecs[0]];
		pa[2]=fTR->ElPy[fElecs[0]];
		pb[1]=fTR->MuPx[fMuons[0]];
		pb[2]=fTR->MuPy[fMuons[0]];
	}

	fMT2 = new Davismt2();
	fMT2->set_momenta(pa, pb, pmiss);
	
	for(int i=0; i<fMT2_histos_number; ++i){
		fMT2->set_mn(i*fMT2_histos_step);
		double MT2 = fMT2->get_mt2();
		if(fLeptConfig==SS_ee  || fLeptConfig==SS_mumu  || fLeptConfig==SS_emu){
			fHMT2_SSll[i] -> Fill(MT2);
		}else if (fLeptConfig==OS_ee){
			fHMT2_OSll[i] -> Fill(MT2);
			fHMT2_OSee[i] -> Fill(MT2);
		}else if (fLeptConfig==OS_mumu){
			fHMT2_OSll[i]   -> Fill(MT2);
			fHMT2_OSmumu[i] -> Fill(MT2);
		}else if (fLeptConfig==OS_emu){
			fHMT2_OSll[i]   -> Fill(MT2);
			fHMT2_OSemu[i]  -> Fill(MT2);
		}
		if (fLeptConfig==OS_ee  || fLeptConfig==OS_mumu  || fLeptConfig==OS_emu){
			for(int mtest=0; mtest<200; mtest++){
				double mass = (double) mtest;
				fMT2->set_mn(mass);
				double MT2 = fMT2->get_mt2();
				fHMT2_vs_M_OSDiLept->Fill(mass, MT2);
			}
		}
	}
	
	delete fMT2;
	
}

// ****************************************************************************************************
vector<TLorentzVector> MassAnalysis::GetLepton4Momenta(){
	vector<TLorentzVector> momenta;
	for(int i=0; i<fElecs.size(); ++i){
		TLorentzVector p;
		p.SetPtEtaPhiM(fTR->ElPt[fElecs[i]],fTR->ElEta[fElecs[i]], fTR->ElPhi[fElecs[i]], fTR->ElE[fElecs[i]]);
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
