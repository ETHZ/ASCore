#ifndef MassAnalysis_hh
#define MassAnalysis_hh


#include <cmath>
#include <string>
#include <iostream>
#include <fstream>

#include <TH1D.h>
#include <TH2D.h>


#include "base/TreeReader.hh"
#include "MultiplicityAnalysisBase.hh"
#include "helper/Davismt2.h"


class MassAnalysis : public MultiplicityAnalysisBase{
public:
	MassAnalysis(TreeReader *tr = NULL);
	virtual ~MassAnalysis();

	void Begin();
	void Analyze();
	void End();


private:
	void MT2();

	
	Davismt2 *fMT2;
	int fMT2_histos_step;
	
	// file for histograms:
	TFile *fHistFile;
	
	TH1D *fHMT2_OS[10];
	TH1D *fHMT2_SS[10];
	TH2D *fHMT2_vs_M_OSDiLept;
	
	
	
};
#endif
