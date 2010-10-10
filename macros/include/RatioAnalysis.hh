#ifndef RatioAnalysis_hh
#define RatioAnalysis_hh


#include <vector>
#include <cmath>
#include <string>
#include <iostream>
#include <fstream>

#include <TH2D.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TStyle.h>
#include <TLatex.h>
#include <TLeaf.h>
#include <TBranch.h>

#include "base/TreeReader.hh"
#include "MultiplicityAnalysisBase.hh"

class RatioAnalysis : public MultiplicityAnalysisBase{
public:
	RatioAnalysis(TreeReader *tr = 0);
	virtual ~RatioAnalysis();
	
	void Begin();
	void Analyze();
	void End();

	

private:
 	void PrintJetRatios();
  	void GetEfficB (double pt, double &effb, double &deffb, double &effc, double &effuds);
  	void BSolve (double p, double dp, double ncorr[], double dncorr[],
	       double npair[], double dnpair[]);
  	void GetGSFE (double n2j, double nbmul[], double p, double dp, 
		double &alp, double &dalp, double ngs[], double dngs[], double ncorr[], double dncorr[]);
  	void PrintChaRatios();
  	void PrintDilRatios();
	void SaveRatioHist();
	void PrintAnomEvts();
	void SaveLeptConfs();


	std::vector<int> fLeptCat;
	unsigned int fNQJets, fNBJets;
	
	 //  ---- counters ---
  	int counter;
  	int fNEvtJets, fNEvt1Jets, fNEvt2Jets, fNEvt3Jets, fNEvt4Jets, fNEvt5Jets, fNEvt6Jets, fNEvt7Jets, fNEvt8Jets;
  	int fNEvtj12near1, fNEvtj12near2,  fNEvtj12near4;
  	int fNEvtJ1l, fNEvt1J1l, fNEvt2J1l, fNEvt3J1l, fNEvt4J1l, fNEvt5J1l, fNEvt6J1l, fNEvt7J1l, fNEvt8J1l;
  	int fNEvtBJets, fNEvt1BJets, fNEvt2BJets, fNEvt3BJets, fNEvt4BJets;
  	int fNEvtB1l, fNEvt1B1l, fNEvt2B1l, fNEvt3B1l, fNEvt4B1l;
  	int fNEvtLeptons;
  	int fNEvtepep, fNEvtepen, fNEvtenen;
  	int fNEvtmpmp, fNEvtmpmn, fNEvtmnmn;
  	int fNEvtepmp, fNEvtepmn, fNEvtenmp, fNEvtenmn;
	int fNEvtepepB, fNEvtepenB, fNEvtenenB;
  	int fNEvtmpmpB, fNEvtmpmnB, fNEvtmnmnB;
	int fNEvtepmpB, fNEvtepmnB, fNEvtenmpB, fNEvtenmnB;
	int fNEvt1BOS2l, fNEvt2BOS2l, fNEvt3BOS2l, fNEvt4BOS2l;
	int fNEvt1BSS2l, fNEvt2BSS2l, fNEvt3BSS2l, fNEvt4BSS2l;
  	int fNEvt1BOSem, fNEvt2BOSem, fNEvt3BOSem, fNEvt4BOSem;

  	int fNevtAnom;

  	double fPtSum1B;
  	double fPtSumsq1B;
  	double fPtSum2B;
  	double fPtSumsq2B;
  	std::vector<int> fIRun, fILumi, fINber, fINlep;
	
  	// ---- Multiplicity Plots Variables ---
  	TFile *fMPHistFile;
  	TH1D *fRatHist;

  	TH1D *fDRJ12;
  	TH1D *fPtJets;
  	TH1D *fPtJ1;
  	TH1D *fPtJ2;
  	TH1D *fEtaJets;
  	TH2D *fRptDPhi2j;
  	TH1D *fAlpT2j;
  	TH1D *fBProbJets;

  	TH1D *fEta1B;
  	TH1D *fEta2B;
  	TH1D *fdPhi2B;
  	TH1D *fDR2B;
  	TH1D *fBProb2B;
  	TH1D *fBProb1B;
  	TH1D *fNjets2B;
  	TH1D *fNjets1B;
  	TH1D *fMass2B;
  	TH1D *fMass1B;
  	TH2D *fMETdPhi1B;
  	TH1D *fdPhiMET1B;
  	TH1D *fNjets2Bnear;
  	TH1D *fMET2Bnear;
  	TH2D *fdPhiJBnear;
  	TH1D *fJMult;
  	TH1D *fJMult1B;
  	TH1D *fJMult2B;
  	TH1D *fJMult3B;
  	TH1D *fJMult4B;
	
  	// ---- efficiencies and fake rates
  	double fEffe, fEffm, fEffb, fdEffb;
  	double fFake, fFakm, fFakb;
	  
};
#endif
