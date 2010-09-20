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
	fMT2_histos_number = 7;
	
	// Define the histograms
	fHMT2_vs_M_OSDiLept = new TH2D("MT2_vs_M_OSDiLept", "MT2_vs_M_OSDiLept", 200, 0, 200, 100, 0, 500);
	for(int i=0; i<fMT2_histos_number; ++i){
		std::stringstream out;
		int mass = i*fMT2_histos_step;
		out << mass ;
		
		TString hist_name_dijet           = "MT2_dijet_M" + (TString) out.str();	
		TString hist_name_diBjet          = "MT2_diBjet_M" + (TString) out.str();
		TString hist_name_pseudojet       = "MT2_pseudojet_M" + (TString) out.str();		
			
		TString hist_name_SSll            = "MT2_SSll_M" + (TString) out.str();		
		TString hist_name_OSll            = "MT2_OSll_M" + (TString) out.str();
		TString hist_name_OSllminusemu    = "MT2_OSllminusemu_M" + (TString) out.str();
		TString hist_name_OSee            = "MT2_OSee_M" + (TString) out.str();
		TString hist_name_OSmumu          = "MT2_OSmumu_M" + (TString) out.str();
		TString hist_name_OSemu           = "MT2_OSemu_M" + (TString) out.str();	

		TString title_dijet           = "Di-Jet MT2, mass = " + (TString) out.str();
		TString title_diBjet          = "Di-bJet MT2, mass = " + (TString) out.str();
		TString title_pseudojet       = "Pseudo-Jet MT2, mass = " + (TString) out.str();
		
		TString title_SSll           = "SS dilept MT2, mass = " + (TString) out.str();
		TString title_OSll           = "OS ll MT2, mass = " + (TString) out.str();
		TString title_OSee           = "OS ee MT2, mass = " + (TString) out.str();
		TString title_OSmumu         = "OS mumu MT2, mass = " + (TString) out.str();
		TString title_OSemu          = "OS emu MT2, mass = " + (TString) out.str();
		TString title_OSllminusemu   = "OS ll minus emu MT2, mass = " + (TString) out.str();

		fHMT2_diBjet[i]          = new TH1D(hist_name_diBjet        , title_diBjet        , 100, 0., 500.);
		fHMT2_dijet[i]           = new TH1D(hist_name_dijet         , title_dijet         , 100, 0., 500.);
		fHMT2_pseudojet[i]       = new TH1D(hist_name_pseudojet     , title_pseudojet     , 100, 0., 500.);
		fHMT2_OSll[i]            = new TH1D(hist_name_OSll          , title_OSll          , 100, 0., 500.);
		fHMT2_OSee[i]            = new TH1D(hist_name_OSee          , title_OSee          , 100, 0., 500.);
		fHMT2_OSmumu[i]          = new TH1D(hist_name_OSmumu        , title_OSmumu        , 100, 0., 500.);
		fHMT2_OSemu[i]           = new TH1D(hist_name_OSemu         , title_OSemu         , 100, 0., 500.);
		fHMT2_OSllminusemu[i]    = new TH1D(hist_name_OSllminusemu  , title_OSllminusemu  , 100, 0., 500.);
		
		fHMT2_SSll[i]            = new TH1D(hist_name_SSll          , title_SSll            , 100, 0., 500.);
		
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
	// Fill MT2 plots
	DiLeptonMT2();
	JetMT2();
	
	
}


double MassAnalysis::GetMT2(double pa[], double pb[], double pmiss[], int m_invisible){
	fMT2 = new Davismt2();
	fMT2->set_momenta(pa, pb, pmiss);
	fMT2->set_mn(m_invisible);
	double MT2=fMT2->get_mt2();
	delete fMT2;
	return MT2;
}

void MassAnalysis::JetMT2(){
	double pa[3], pb[3], pmiss[3];
	pa[0]=0;  // corresponds to mass uf particle a
	pa[1]=0;  // px of particle a
	pa[2]=0;  // py of particle a
	pb[0]=0;  // corresponds to mass uf particle b
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
		pa[1]=fTR->JPx[fBJets[0]];
		pa[2]=fTR->JPy[fBJets[0]];
		pb[1]=fTR->JPx[fBJets[1]];
		pb[2]=fTR->JPy[fBJets[1]];
		
		for(int i=0; i<fMT2_histos_number; ++i){
			double MT2=GetMT2(pa, pb, pmiss, i*fMT2_histos_step);
			fHMT2_diBjet[i] ->Fill(MT2);
		}	
	}
	
	// pseudojets
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
		
		fHemisphere = new Hemisphere(px, py, pz, E, 2, 3);
		vector<int> grouping = fHemisphere->getGrouping();

		for(int i=0; i<fJets.size(); ++i){
			if(grouping[i]==1){
				pa[1] +=px[i];
				pa[2] +=py[i];
			}else if(grouping[i] == 2){
				pb[1] +=px[i];
				pb[2] +=py[i];				
			}
		}
		delete fHemisphere;
		for(int i=0; i<fMT2_histos_number; ++i){
			double MT2=GetMT2(pa, pb, pmiss, i*fMT2_histos_step);
			fHMT2_pseudojet[i] ->Fill(MT2);
		}
	} 
}

void MassAnalysis::DiLeptonMT2(){
	double pa[3], pb[3], pmiss[3];
	pa[0]=0;  // corresponds to mass uf particle a
	pb[0]=0;  // corresponds to mass uf particle b
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
		
		fHMT2_dijet[i]       ->Write();
		fHMT2_diBjet[i]      ->Write();
		fHMT2_pseudojet[i]   ->Write();
	}
	fHMT2_vs_M_OSDiLept      ->Write();
	
	fHistFile->Close();
}
