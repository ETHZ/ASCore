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
#include "MT2tree.hh"


class MassAnalysis : public MultiplicityAnalysisBase{
public:
	MassAnalysis(TreeReader *tr = NULL);
	virtual ~MassAnalysis();

	void Begin(const char* filename = "Mass_histos.root");
	void Analyze();
	void End();


private:
	void PseudoJetMasses();

	void BookTree();
	void FillTree();
	void ResetTree(); 

	struct HemiObject {string type; int index; int hemi;} fHemiObject;
	struct HemiObjects{
		vector<HemiObject> objects;
		TLorentzVector pjet1;
		TLorentzVector pjet2;
		TLorentzVector UTM ;
		double MT2;
		double MCT;
		int seed;
		int assoc;	
	} fHemiObjects[3];
	
	void   GetMT2Variables(int hemi_seed, int hemi_assoc, double maxDR, double minJPt, double maxJEta, HemiObjects& hemiobject);
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

	vector<int> interesting_Run;
	vector<int> interesting_Lumi;
	vector<int> interesting_Event;
	vector<double> interesting_value;
	vector<string> interesting_Type;

	// file for histograms:
	TFile* fHistFile;

        // MT2 mini-tree
        MT2tree* fMT2tree;
        TTree* fATree;
	
};
#endif
