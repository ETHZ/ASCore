#ifndef MultiplicityAnalysisBase_hh
#define MultiplicityAnalysisBase_hh


#include <cmath>
#include <string>
#include <iostream>
#include <fstream>
#include <TLorentzVector.h>

#include "base/TreeReader.hh"
#include "base/UserAnalysisBase.hh"

class MultiplicityAnalysisBase : public UserAnalysisBase{
public:
	MultiplicityAnalysisBase(TreeReader *tr = NULL);
	virtual ~MultiplicityAnalysisBase();

	void ReadCuts(const char* SetofCuts);

	bool IsSelectedEvent();
	void InitializeEvent();
	double GetDiLeptInvMass();
	
	vector<int> fElecs;
	vector<int> fTaus;
	vector<int> fMuons;
	vector<int> fJets;
	
	
	struct JetTau {
		int NObjs;
		vector<int> index; 
		vector<double> pt; 
		vector<int> isTau;
		void reset(){
			NObjs=0;
			index.clear();
			pt.clear();
			isTau.clear();
		}
	} fJetTaus;

	enum LeptConfig {
	 	e, mu, OS_emu, OS_ee, OS_mumu, SS_emu, SS_ee, SS_mumu, multilept, null
  	};
	LeptConfig fLeptConfig;
	
	
	//  ---- set of cuts ---
	TString fSetName;
	float fCut_PFMET_min;
	float fCut_HT_min;
	float fCut_JPt_hardest_min;
	float fCut_JPt_second_min;
	float fCut_DiLeptInvMass_min;
	float fCut_DiLeptInvMass_max;
	float fCut_PtHat_max;
	int   fCut_Run_min;
	int   fCut_Run_max;

	// members
	float fHT;
	int   fNJets_toremove_ele;
	int   fNJets_toremove_muo;

private:
	void FindLeptonConfig();
	void GetLeptonJetIndices();

	
	// ---- required and vetoed triggers ----
	std::vector<std::string> fRequiredHLT; 
	std::vector<std::string> fVetoedHLT;

};
#endif
