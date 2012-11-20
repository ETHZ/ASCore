#ifndef SSDLAnalysis_hh
#define SSDLAnalysis_hh

#include <vector>
#include <cmath>
#include <string>
#include <iostream>
#include <fstream>
#include <time.h>

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
#include "helper/Utilities.hh"
#include "base/UserAnalysisBase.hh"
#include "helper/Monitor.hh"
#include "helper/PUWeight.h"

using namespace std;

class SSDLAnalysis : public UserAnalysisBase{
public:
	SSDLAnalysis(TreeReader *tr = 0);
	virtual ~SSDLAnalysis();
	
	void Begin(const char* filename = "SSDLTree.root");
	void Analyze();
	void End();
	
	inline void DoFillEffTree(bool fill){fDoFillEffTree = fill;};
	
	void ReadTriggers(const char* = "HLTPaths_SSDL.dat");
	void AddTriggerBranches();
	bool FillTriggers(); // Returns OR of list of triggers
	void FillAnalysisTree();
	void BookTree();
	void ResetTree();
	
	void FillEffTree();
	void BookEffTree();
	void ResetEffTree();

	bool IsSignalMuon(int);
        bool IsSignalMuon(int, int&, int&, int&);
	bool IsSignalElectron(int);
        bool IsSignalElectron(int, int&, int&, int&);
	bool IsTightMuon(int, int);
	bool IsTightEle(int, int);
	
	int JetPartonMatch(int);
	int GenJetMatch(int);
	
	double corrMuIso(int);
	double corrElIso(int);
	
	const bool AddBranch(const char* name, const char* type, void* address, const char* size = 0);

	struct HLTPathSet{
		TString name;
		vector<string> paths;
	};

private:
	TH1F *fHEvCount;
	
	bool fDoFillEffTree;
	
	Monitor fCounter;
	string fCutnames[4];

	static const int fMaxNjets = 40;
	static const int fMaxNmus  = 5;
	static const int fMaxNeles = 5;
	static const int fMaxNtaus = 5;
	
	static const int nx = 3;
	static const float x_values[nx]; // = {0.05, 0.5, 0.95};

	static TString gBaseDir;
	
	JetCorrectionUncertainty *fJetCorrUnc;
	
	TTree* fAnalysisTree;
	TH2D* fMsugraCount;
	TH2D* fProcessCount[10];
	TH2D* fRightHandedSlepCountAll;
	TH2D* fRightHandedCountAll;
	TH2D* fTChiSlepSlepCountAll;
	TH2D* fTChiSlepSnuCountAll;
	TH2D* fModelCountAll;
	TH2D* fRightHandedSlepCount[nx];
	TH2D* fRightHandedCount[nx];
	TH2D* fTChiSlepSlepCount[nx];
	TH2D* fTChiSlepSnuCount[nx];
	TH2D* fModelCount[nx];
	
	/////////////////////////////////////
	// Tree branches
	int   fTRunNumber;
	int   fTEventNumber;
	int   fTLumiSection;

	float  fTm0;
	float  fTm12;
	int    fTprocess;
	float  fTmGlu;
	float  fTmLSP;
	int    fTisTChiSlepSnu;
	int    fTisRightHanded;

	// PileUP info
	float fTrho;
	int   fTnvrtx;
	float fTpuweight;
	
	// Triggers
	vector<HLTPathSet> fHLTPathSets;
	vector<string> fHLTPaths;
	vector<int>    fHLTResults;
	vector<int>    fHLTPrescales;	

	// JetMET properties
	int   fTnqjets;
	float fTJetpt [fMaxNjets];
	float fTJeteta[fMaxNjets];
	float fTJetphi[fMaxNjets];
	float fTJetenergy[fMaxNjets];
	float fTJetbtag1[fMaxNjets]; // CSV
	float fTJetbtag2[fMaxNjets]; // JetProb
	// MARC float fTJetbtag3[fMaxNjets]; // TCHP
	// MARC float fTJetbtag4[fMaxNjets]; // TCHE
	float fTJetArea[fMaxNjets];
	float fTJetJEC[fMaxNjets];
	int   fTJetPartonID[fMaxNjets];
	float fTJetGenpt [fMaxNjets];
	float fTJetGeneta[fMaxNjets];
	float fTJetGenphi[fMaxNjets];

	float fTpfMET;
	float fTpfMETphi;
	float fTpfMETType1;
	float fTpfMETType1phi;
	
	// Muon properties
	int   fTnqmus;
	int   fTIsSignalMuon  [fMaxNmus];
	float fTmupt          [fMaxNmus];
	float fTmueta         [fMaxNmus];
	float fTmuphi         [fMaxNmus];
	float fTmudetiso      [fMaxNmus];
	float fTmupfiso       [fMaxNmus];
	float fTmupfchiso     [fMaxNmus];
	float fTmupfneiso     [fMaxNmus];
	float fTmuradiso      [fMaxNmus];
	int   fTmucharge      [fMaxNmus];
	float fTmud0          [fMaxNmus];
	float fTmudz          [fMaxNmus];
	float fTmuptE         [fMaxNmus];
	float fTmuEMVetoEt    [fMaxNmus];
	float fTmuHadVetoEt   [fMaxNmus];
	int   fTmuPassesTightID[fMaxNmus];
	int   fTmuid          [fMaxNmus];
	int   fTmumoid        [fMaxNmus];
	int   fTmugmoid       [fMaxNmus];
	int   fTmutype        [fMaxNmus];
	int   fTmumotype      [fMaxNmus];
	int   fTmugmotype     [fMaxNmus];
	float fTmuMT          [fMaxNmus];
	
	// Electron properties
	int   fTnqels;
	int   fTIsSignalElectron [fMaxNeles];
	int   fTElcharge         [fMaxNeles];
	int   fTElChargeIsCons   [fMaxNeles];
	float fTElpt             [fMaxNeles];
	float fTEleta            [fMaxNeles];
	float fTElphi            [fMaxNeles];
	float fTEld0             [fMaxNeles];
	float fTElD0Err          [fMaxNeles];
	float fTEldz             [fMaxNeles];
	float fTElDzErr          [fMaxNeles];
	float fTElDetIso         [fMaxNeles];
	float fTElPFIso          [fMaxNeles];
	float fTElPFchiso        [fMaxNeles];
	float fTElPFneiso        [fMaxNeles];
	float fTElRadIso         [fMaxNeles];
	float fTElMVAIDnoTrig    [fMaxNeles];
	float fTElMVAIDTrig      [fMaxNeles];
	float fTElEcalRecHitSumEt[fMaxNeles];
	float fTElHcalTowerSumEt [fMaxNeles];
	float fTElTkSumPt        [fMaxNeles];
	float fTElDPhi           [fMaxNeles];
	float fTElDEta           [fMaxNeles];
	float fTElSigmaIetaIeta  [fMaxNeles];
	float fTElHoverE         [fMaxNeles];
	float fTElEPthing        [fMaxNeles];
	// MARC int   fTElIsGoodElId_WP80[fMaxNeles];
	// MARC int   fTElIsGoodElId_WP90[fMaxNeles];
	int   fTElIsGoodElId_LooseWP[fMaxNeles];
	int   fTElIsGoodElId_MediumWP[fMaxNeles];
	int   fTElIsGoodTriggerEl[fMaxNeles];
	float fTElMT             [fMaxNeles];
	int   fTElGenID          [fMaxNeles];
	int   fTElGenMID         [fMaxNeles];
	int   fTElGenGMID        [fMaxNeles];
	int   fTElGenType        [fMaxNeles];
	int   fTElGenMType       [fMaxNeles];
	int   fTElGenGMType      [fMaxNeles];

	// Tau properties
	int fTnqtaus;
	int   fTTaucharge      [fMaxNtaus];
	float fTTaupt          [fMaxNtaus];
	float fTTaueta         [fMaxNtaus];
	float fTTauphi         [fMaxNtaus];
	float fTTauMVAElRej    [fMaxNtaus];
	float fTTauTightMuRej  [fMaxNtaus];
	float fTTauLCombIsoDB  [fMaxNtaus];

	TTree *fLepEffTree; // lepton efficiency tree, filled once per lepton
	int   fLETrun;
	int   fLETevent;
	int   fLETlumi;
	float fLETrho;
	int   fLETnvrtx;
	float fLETpuweight;
	int   fLETtype; // mu(0), el(1)
	float fLETpt;
	float fLETeta;
	float fLETphi;
	float fLETiso;
	int   fLETpassed1;
	int   fLETpassed2;
	int   fLETpassed3;
	int   fLETpassed4;
	
};
#endif
