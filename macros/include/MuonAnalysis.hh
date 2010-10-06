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
	virtual void AnalyzeDi();
	virtual void AnalyzeSS();
	virtual void End();

	virtual void FillMuonTree(int, int = -1);

	virtual void BookTree();
	virtual void ResetTree();
	
private:	
	static const int gMaxnjets = 30;
  
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
	double fTjpt[gMaxnjets];
	double fTjeta[gMaxnjets];
	double fTjphi[gMaxnjets];
	double fTHT;
	double fTMHT;
	double fTMET;
	double fTMT;
	double fTMinv;
	double fTSumET; // Directly from NTupleProducer, i.e. scalar sum of all calotower.et()
	int fTnmus;
	double fTmupt[2];
	double fTmueta[2];
	double fTmuphi[2];
	double fTmuiso[2];
	int fTmucharge[2];
	int fTmutight[2]; // 0 for loose (but not tight), 1 for tight
	double fTDRjet[2];
	double fTDPhijet[2];
	double fTmucalocomp[2];
	double fTmusegmcomp[2];
	double fTmuouterrad[2];
	double fTmunchi2[2];
	int fTmuntkhits[2];
	double fTmud0[2];
	double fTmudz[2];
	double fTmuptE[2];
	int fTmuid[2];
	int fTmumoid[2];
	int fTmugmoid[2];
	int fTmutype[2];
	int fTmumotype[2];
	int fTmugmotype[2];
	
};
#endif
