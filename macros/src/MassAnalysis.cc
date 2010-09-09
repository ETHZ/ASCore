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
	
	// Define the histograms
	fHMT2_vs_M_OSDiLept = new TH2D("MT2_vs_M_OSDiLept", "MT2_vs_M_OSDiLept", 200, 0, 200, 100, 0, 500);
	fMT2_histos_step = 25;
	for(int i=0; i<10; ++i){
		std::stringstream out;
		int mass = i*fMT2_histos_step;
		out << mass ;
		TString hist_name_OS= "MT2_OS_M" + (TString) out.str();
		TString hist_name_SS= "MT2_SS_M" + (TString) out.str();
		TString title_SS = "SS dilept MT2, mass = " + (TString) out.str();
		TString title_OS = "OS dilept MT2, mass = " + (TString) out.str();
		fHMT2_OS[i]     = new TH1D(hist_name_OS , title_OS, 100, 0., 500.);
		fHMT2_SS[i]     = new TH1D(hist_name_SS , title_SS, 100, 0., 500.);
		
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
	MT2();
	
}





void MassAnalysis::MT2(){
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
	
	if(fLeptConfig==SS_ee  || fLeptConfig==SS_mumu  || fLeptConfig==SS_emu) {
		for(int i=0; i<10; ++i){
			fMT2->set_mn(i*fMT2_histos_step);
			double MT2 = fMT2->get_mt2();
			fHMT2_SS[i] -> Fill(MT2);
		}
	} else if(fLeptConfig==OS_ee  || fLeptConfig==OS_mumu || fLeptConfig==OS_emu ) {
		for(int i=0; i<10; ++i){
			fMT2->set_mn(i*fMT2_histos_step);
			double MT2 = fMT2->get_mt2();
			fHMT2_OS[i] -> Fill(MT2);
		}
		for(int mtest=0; mtest<200; mtest++){
			double mass = (double) mtest;
			fMT2->set_mn(mass);
			double MT2 = fMT2->get_mt2();
			fHMT2_vs_M_OSDiLept->Fill(mass, MT2);
		}
	}
	
	delete fMT2;
	
}

void MassAnalysis::End(){
	fHistFile->cd();	
	for(int i=0; i<10; ++i){
		fHMT2_OS[i]->Write();
		fHMT2_SS[i]->Write();
	}
	fHMT2_vs_M_OSDiLept ->Write();
	
	fHistFile->Close();
}
