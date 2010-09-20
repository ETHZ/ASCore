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
	

	Davismt2 *fMT2;
	Hemisphere *fHemisphere;
	
	// data members
	int fMT2_histos_step;
  	int fMT2_histos_number;
	
	// file for histograms:
	TFile *fHistFile;

	// histos
	TH1D *fHMT2_SSll[10];
	
	TH1D *fHMT2_dijet[10];
	TH1D *fHMT2_diBjet[10];
	TH1D *fHMT2_pseudojet[10];
	
	TH1D *fHMT2_OSll[10];
	TH1D *fHMT2_OSee[10];
	TH1D *fHMT2_OSmumu[10];
	TH1D *fHMT2_OSemu[10];
	TH1D *fHMT2_OSllminusemu[10];
	TH2D *fHMT2_vs_M_OSDiLept;
	
	
	
};
#endif
