#ifndef MassAnalysis_hh
#define MassAnalysis_hh


#include <cmath>
#include <string>
#include <iostream>
#include <fstream>

#include <TH1D.h>
#include <TH2D.h>

#include <numeric>


#include "base/TreeReader.hh"
#include "MultiplicityAnalysisBase.hh"
#include "helper/Davismt2.h"
#include "helper/Hemisphere.hh"


class MassAnalysis : public MultiplicityAnalysisBase{
public:
	MassAnalysis(TreeReader *tr = NULL);
	virtual ~MassAnalysis();

	void Begin();
	void Analyze();
	void End();


private:
	void DiLeptonMT2();
	void JetMT2();
	double GetMT2(double pa[], double pb[], double pmiss[], int m_invisible);
	double GetAlphaT(std::vector<TLorentzVector>& p4s);
	std::vector<double> DeltaSumPt_permutations(std::vector<TLorentzVector>& p4s);
	double GetMCT(TLorentzVector p1, TLorentzVector p2);
	double GetMCTperp(TLorentzVector p1, TLorentzVector p2, TLorentzVector P_UTM);
	double GetMT2perp(TLorentzVector p1, TLorentzVector p2, TLorentzVector P_UTM, double m_inv);
	TVector3 GetMomPerp(TLorentzVector p, TLorentzVector P_UTM);
	void MassesForTTbar();
	vector<TLorentzVector> GetLepton4Momenta();



	Davismt2 *fMT2;
	Hemisphere *fHemisphere;
	
	// data members
	int fMT2_histos_step;
  	int fMT2_histos_number;
	
	// file for histograms:
	TFile* fHistFile;

	// histos
	TH1D* fHMT2_SSll[10];
	      
	TH1D* fHMT2_dijet[10];
	TH1D* fHMT2_diBjet[10];
	TH1D* fHMT2_pseudojet[10];
	TH1D* fHMT2_diBjetdiLept[10];
	TH1D* fHMT2_PseudoJetsWithB[10];      
	TH1D* fHMT2_PseudoJetsWithLeptons[10];
	TH1D* fHMT2_PseudoJetClean[10];       

	TH1D* fHMT2_OSll[10];
	TH1D* fHMT2_OSee[10];
	TH1D* fHMT2_OSmumu[10];
	TH1D* fHMT2_OSemu[10];
	TH1D* fHMT2_OSllminusemu[10];
	
	TH1D* fHMT2perp_OSee[10];  
	TH1D* fHMT2perp_OSmumu[10];
	TH1D* fHMT2perp_OSemu[10];
	      
	TH1D* fHAlpahT_DiJet;      
	TH1D* fHAlpahT_PseudoJet; 
	    
	TH1D* fHMCT_DiJet;   
	TH1D* fHMCT_DiBJet;  
	TH1D* fHMCT_OSee;    
	TH1D* fHMCT_OSmumu;  
	TH1D* fHMCT_OSemu;   
	
	TH1D* fHMCTperp_OSee;  
	TH1D* fHMCTperp_OSmumu;
	TH1D* fHMCTperp_OSemu; 

	TH1D* fHMCT_TTbar;	
	TH1D* fHMCTperp_TTbar;	
	TH1D* fHMT2_TTbar;
	TH1D* fHMT2perp_TTbar;
	
	TH1D* fHInvMassOSee;
	TH1D* fHInvMassOSmumu;
	TH1D* fHInvMassOSemu;
	TH1D* fHInvMassOSll;
	
	TH1D* fHInvMassSSee;
	TH1D* fHInvMassSSmumu;
	TH1D* fHInvMassSSemu; 	
	TH1D* fHInvMassSSll;  	
	
	
	TH2D *fHMT2_vs_M_OSDiLept;
	
	
	
};
#endif
