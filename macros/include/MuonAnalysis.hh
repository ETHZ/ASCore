#ifndef MuonAnalysis_hh
#define MuonAnalysis_hh

#include <vector>
#include <cmath>
#include <string>
#include <iostream>
#include <fstream>

#include <TH1D.h>
#include <TH2D.h>

#include "base/TreeReader.hh"
#include "base/UserAnalysisBase.hh"

class MuonAnalysis : public UserAnalysisBase{
public:
	MuonAnalysis(TreeReader *tr = 0);
	virtual ~MuonAnalysis();

	virtual void Begin();
	virtual void Analyze();
	virtual void End();

	virtual void FillMuonTree(std::vector<int>);

	virtual void BookTree();
	virtual void ResetTree();
	
private:	
	static const int gMaxnjets = 30;
	static const int gMaxnmus = 5;
  
	TTree *fMuTree;
	int fTrun;
	int fTevent;
	int fTlumisec;
	int fT_HLTMu9;
	int fT_HLTDoubleMu3;
	int fT_HLTJet30U;
	int fT_HLTJet50U;
	int fT_HLTJet70U;
	int fT_HLTJet100U;
	int fTnjets;
	float fTjpt[gMaxnjets];
	float fTjeta[gMaxnjets];
	float fTjphi[gMaxnjets];
	float fTHT;
	float fTMHT;
	float fTMET;
	float fTMT;
	float fTMinv;
	float fTSumET; // Directly from NTupleProducer, i.e. scalar sum of all calotower.et()
	int fTnmus;
	float fTmupt[gMaxnmus];
	float fTmueta[gMaxnmus];
	float fTmuphi[gMaxnmus];
	float fTmuiso[gMaxnmus];
	int fTmucharge[gMaxnmus];
	int fTmutight[gMaxnmus]; // 0 for loose (but not tight), 1 for tight
	float fTDRjet[gMaxnmus];
	float fTDRhardestjet[gMaxnmus];
	float fTDPhijet[gMaxnmus];
	float fTmucalocomp[gMaxnmus];
	float fTmusegmcomp[gMaxnmus];
	float fTmuouterrad[gMaxnmus];
	float fTmunchi2[gMaxnmus];
	int fTmuntkhits[gMaxnmus];
	int fTmunmuhits[gMaxnmus];
	float fTmuemvetoet[gMaxnmus];
	float fTmuhadvetoet[gMaxnmus];
	float fTmud0[gMaxnmus];
	float fTmudz[gMaxnmus];
	float fTmuptE[gMaxnmus];
	int fTmuid[gMaxnmus];
	int fTmumoid[gMaxnmus];
	int fTmugmoid[gMaxnmus];
	int fTmutype[gMaxnmus];
	int fTmumotype[gMaxnmus];
	int fTmugmotype[gMaxnmus];
	
};
#endif
