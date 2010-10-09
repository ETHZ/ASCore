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
	void SingleLeptonMasses();
	void DiLeptonMasses();
	void DiJetMasses();
	void DiBJetMasses();
	void MultiBMasses();
	void PseudoJetMasses();
	void OSDiLeptonMasses();
	void MassesForTTbar();
	void ControlPlots();
	void VectorSumPt();

	double GetMT2(TLorentzVector v1, double mv1, TLorentzVector v2, double mv2, TLorentzVector p_unobs, int m_invisible); 
	double GetAlphaT(std::vector<TLorentzVector>& p4s);
	std::vector<double> DeltaSumPt_permutations(std::vector<TLorentzVector>& p4s);
	double GetMCT(TLorentzVector p1, TLorentzVector p2);
	double GetMCTperp(TLorentzVector p1, TLorentzVector p2, TLorentzVector P_UTM);
	double GetMT2perp(TLorentzVector p1, TLorentzVector p2, TLorentzVector P_UTM, double m_inv);
	TVector3 GetMomPerp(TLorentzVector p, TLorentzVector P_UTM);
	vector<TLorentzVector> GetLepton4Momenta();
	bool IsCleanJetEvent();

	Davismt2 *fMT2;
	Hemisphere *fHemisphere;
	
	// data members
	int fMT2_histos_step;
  	int fMT2_histos_number;
	float fVectorSumPt;
	vector<int> interesting_Run;
	vector<int> interesting_Lumi;
	vector<int> interesting_Event;
	vector<string> interesting_Type;

	// file for histograms:
	TFile* fHistFile;

	// histos
	TH1D* fHMT2_SSll[10];
	      
	TH1D* fHMT2_dijet[10];
	TH1D* fHMT2_diBjet[10];
	TH1D* fHMT2_pseudojet[10];
	TH1D* fHMT2_diBjetdiLept[10];
	TH1D* fHMT2_PseudoJetWithB[10];      
	TH1D* fHMT2_PseudoJetWithLeptons[10];
	TH1D* fHMT2_PseudoJetNoLeptons[10];

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
	    
	TH1D* fHMCT_PseudoJet;
	TH1D* fHMCT_PseudoJetWithLeptons;
	TH1D* fHMCT_PseudoJetNoLeptons;
	TH1D* fHMCT_PseudoJetWithB;
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
	TH1D* fHInvMassDijet;	
	TH1D* fHInvMassdiBHemi;
	TH1D* fHInvMassDiBjet;

	TH1D* fHJpt;   
	TH1D* fHJEta;  
	TH1D* fHBJpt; 
	TH1D* fHBJEta;
       	TH1D* fHVectorSumPt;	
	TH2D* fHVectorSumPtvsDiJetMT2;
	TH2D* fHVectorSumPtvsPseudoJetMT2;

	TH1D* fHMT_single_e;
	TH1D* fHMT_single_mu;
	TH1D* fHMT_single_e_nojets; 
	TH1D* fHMT_single_mu_nojets;

	TH2D *fHMT2_vs_M_OSDiLept;
	
	
	
};
#endif
