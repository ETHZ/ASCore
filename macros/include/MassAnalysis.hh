#ifndef MassAnalysis_hh
#define MassAnalysis_hh


#include <cmath>
#include <string>
#include <iostream>
#include <fstream>

#include <TH1D.h>
#include <TH2D.h>


#include "base/TreeReader.hh"
#include "base/UserAnalysisBase.hh"
#include "helper/Davismt2.h"


class MassAnalysis : public UserAnalysisBase{
public:
	MassAnalysis(TreeReader *tr = NULL);
	virtual ~MassAnalysis();

	void Begin();
	void Analyze();
	void End();
	void Reset();


private:
	
	void FindLeptonConfig();
	void MT2();
	
	std::vector<int> jets;
	std::vector<int> bjets;
	std::vector<int> elecs;
	std::vector<int> muons;
		
	enum LeptConfig {
	 	e, mu, OS_emu, OS_ee, OS_mumu, SS_emu, SS_ee, SS_mumu, null
  	};
	LeptConfig fLeptConfig;
	
	Davismt2 *fMT2;
	int fMT2_histos_step;
	
	// file for histograms:
	TFile *fHistFile;
	
	TH1D *fHMT2_OS[10];
	TH1D *fHMT2_SS[10];
	TH2D *fHMT2_vs_M_OSDiLept;
	
	
	
};
#endif
