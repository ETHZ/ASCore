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
#include "helper/TMctLib.h"
#include "helper/Hemisphere.hh"


class MassAnalysis : public MultiplicityAnalysisBase{
public:
	MassAnalysis(TreeReader *tr = NULL);
	virtual ~MassAnalysis();

	void Begin();
	void Analyze();
	void End();
	void SetWeight(float weight);


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
	void PrintEvent();

	void BookTree();
	void FillTree();
	void ResetTree(); 

	double GetMT2(TLorentzVector v1, double mv1, TLorentzVector v2, double mv2, TLorentzVector p_unobs, int m_invisible); 
	double GetAlphaT(std::vector<TLorentzVector>& p4s);
	std::vector<double> DeltaSumPt_permutations(std::vector<TLorentzVector>& p4s);
	double GetMCT(TLorentzVector p1, TLorentzVector p2);
	double GetMCTperp(TLorentzVector p1, TLorentzVector p2, TLorentzVector P_UTM);
	double GetToveyMCTperp(TLorentzVector p1, TLorentzVector p2, TLorentzVector DTM, TVector2 pmiss);
	double GetMCTcorr(TLorentzVector p1, TLorentzVector p2, TLorentzVector DTM, TVector2 pmiss);
	double GetMCTcorr(TLorentzVector p1, TLorentzVector p2, TLorentzVector UTM);
	double GetMT2perp(TLorentzVector p1, TLorentzVector p2, TLorentzVector P_UTM, double m_inv);
	TVector3 GetMomPerp(TLorentzVector p, TLorentzVector P_UTM);
	vector<TLorentzVector> GetLepton4Momenta();

	Davismt2 *fMT2;
	TMctLib  *fMCT;
	Hemisphere *fHemisphere;
	
	// data members
	int fMT2_histos_step;
  	int fMT2_histos_number;
	float fWeight;

	vector<int> interesting_Run;
	vector<int> interesting_Lumi;
	vector<int> interesting_Event;
	vector<double> interesting_value;
	vector<string> interesting_Type;

	// file for histograms:
	TFile* fHistFile;

	// Tree
	TTree* fMassTree;
	static const int gMaxnjets = 30;

	int    fTrun;
	int    fTevent;
	int    fTlumisec;
	float  fTweight;
		
	int    fTnjets;
	int    fTnelecs;
	int    fTnmuons;
	int    fTleptconfig;
	
	int    fTisCleanMultiJetEvent;
	int    fTisCleanJetEvent;
	int    fTnJetsPt50Eta25;
	int    fTr12r21;
	
	double fTpseudoJetMT2;
	double fTpseudoJetMCT;
	double fTpseudoJet1Pt;
	double fTpseudoJet2Pt;
	double fTpseudojetAlphaT;
	double fTleadingJetEta;
	double fTvectorsumpt;
	double fTpfmet;
	double fTpfmetsign;
	double fTmHT;
	double fThT;
	double fTmpt_sel;
	double fTmpt;
	double fTdPhiJ1MET;
	double fTdPhiJ2MET;
	double fTPseudoJetMT2AxisdPhi;
	double fTr1221min;
	double fTdPhiMhtMpt;
	double fTjpt[gMaxnjets];

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
	TH1D* fHMCTcorr_PseudoJet;
	TH1D* fHMCTcorr_PseudoJetWithLeptons;
	TH1D* fHMCTcorr_PseudoJetNoLeptons;
	TH1D* fHMCTcorr_PseudoJetWithB;
	TH1D* fHMCTcorr_DiJet;   
	TH1D* fHMCTcorr_DiBJet;  
	TH1D* fHMCTcorr_OSee;    
	TH1D* fHMCTcorr_OSmumu;  
	TH1D* fHMCTcorr_OSemu;   

	TH1D* fHMCTperp_OSee;  
	TH1D* fHMCTperp_OSmumu;
	TH1D* fHMCTperp_OSemu; 

	TH1D* fHMCT_TTbar;	
	TH1D* fHMCTcorr_TTbar;
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

	TH2D* fHDPhiJ1vsDPhiJ2;
	TH2D* fHDPhiJ1vsDPhiJ2_clean;
	TH1D* fHMHT;
	TH1D* fHMHT_clean;
	TH1D* fHPFMET;
	TH1D* fHPFMET_clean;
	TH1D* fHJEta;  
	TH1D* fHJEta_clean;  
	TH1D* fHJpt;   
	TH1D* fHJpt_clean;  	
	TH1D* fHBJpt; 
	TH1D* fHBJEta;
	TH1D* fHNJets;
	TH1D* fHNJets_clean;
       	TH1D* fHVectorSumPt;	
	TH1D* fHVectorSumPt_clean;
	TH1D* fHLeptConfig;
	TH1D* fHLeptConfig_clean;		
	TH1D* fHPseudoJetMT2AxisdPhi;
	TH1D* fHMPT;
	TH1D* fHMPT_selected;

	TH2D* fHVectorSumPtvsDiJetMT2;

	TH2D* fHPseudoJetMT2vsVectorSumPt;
	TH2D* fHPseudoJetMT2vsMETsign;
	TH2D* fHPseudoJetMT2vsMET;
	TH2D* fHPseudoJetMT2vsAlphaT;
	TH2D* fHPseudoJetMT2vsLeadingJEta;
	TH2D* fHPseudoJetMT2vsMHT;

	TH1D* fHMT_single_e;
	TH1D* fHMT_single_mu;
	TH1D* fHMT_single_e_nojets; 
	TH1D* fHMT_single_mu_nojets;

	TH2D *fHMT2_vs_M_OSDiLept;
	
	TH1D *fHPFandCalo_deltaR;
	TH1D *fHPFJ1J2DeltaR;	
	
};
#endif
