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
	
	void ReadTriggers(const char* = "HLTPaths_SSDL.dat");
	void BookTriggers(vector<string>);
	bool FillTriggers(vector<string>); // returns OR of list of triggers
	void BookTree();
	void ResetTree();
	
	const bool AddBranch(const char* name, const char* type, void* address, const char* size = 0);

private:
	Monitor fCounter;
	string fCutnames[4];

	static const int fMaxNjets = 30;
	static const int fMaxNmus  = 5;
	static const int fMaxNeles = 5;
	
	static const int gMaxhltbits = 300;
	
	static TString gBaseDir;
	
	TTree* fAnalysisTree;
	
	// run/sample properties
	int   fTRunNumber;
	int   fTEventNumber;
	int   fTLumiSection;

	// triggers
	int fNHLTPaths;
	vector<string> fHLTPaths;
	vector<int>    fHLTResults;
	vector<int>    fHLTPrescales;	

	// jet-MET properties
	int   fTnqjets;
	float fTJetpt [fMaxNjets];
	float fTJeteta[fMaxNjets];
	float fTJetphi[fMaxNjets];
	float fTJetbtag[fMaxNjets]; // tight WP: > 2.
	float fTtcMET;
	float fTpfMET;
	
	// event properties
	float fTrho;
	int fTnvrtx;
	float fTpuweight;
	
	//muon properties
	int   fTnqmus;
	float fTmupt          [fMaxNmus];
	float fTmueta         [fMaxNmus];
	float fTmuphi         [fMaxNmus];
	float fTmuiso         [fMaxNmus];
	float fTmuisohyb      [fMaxNmus];
	int   fTmucharge      [fMaxNmus];
	int   fTmutight       [fMaxNmus]; // 0 for loose (but not tight), 1 for tight
	float fTmud0          [fMaxNmus];
	float fTmudz          [fMaxNmus];
	float fTmud0bs        [fMaxNmus];
	float fTmudzbs        [fMaxNmus];
	float fTmuptE         [fMaxNmus];
	int   fTmuid          [fMaxNmus];
	int   fTmumoid        [fMaxNmus];
	int   fTmugmoid       [fMaxNmus];
	int   fTmutype        [fMaxNmus];
	int   fTmumotype      [fMaxNmus];
	int   fTmugmotype     [fMaxNmus];
	float fTmuMT          [fMaxNmus];
	
	// electron properties
	int     fTnqels;
	int   fTElcharge                   [fMaxNeles];
	int   fTElChargeIsCons             [fMaxNeles];
	// int   fTElChargeIsGenCons          [fMaxNeles];
	int   fTElEcalDriven               [fMaxNeles];
	float fTElCaloEnergy               [fMaxNeles];
	float fTElpt                       [fMaxNeles];
	float fTEleta                      [fMaxNeles];
	float fTElphi                      [fMaxNeles];
	float fTEld0                       [fMaxNeles];
	float fTElD0Err                    [fMaxNeles];
	float fTEldz                       [fMaxNeles];
	float fTElDzErr                    [fMaxNeles];
	float fTElEoverP                   [fMaxNeles];
	float fTElHoverE                   [fMaxNeles];
	float fTElSigmaIetaIeta            [fMaxNeles];
	float fTElDeltaPhiSuperClusterAtVtx[fMaxNeles];
	float fTElDeltaEtaSuperClusterAtVtx[fMaxNeles];
	float fTElRelIso                   [fMaxNeles];
	int   fTElIsGoodElId_WP80          [fMaxNeles];
	int   fTElIsGoodElId_WP90          [fMaxNeles];
	int   fTElIsGoodElId_WP95          [fMaxNeles];
	float fTElS4OverS1                 [fMaxNeles];
	float fTElConvPartnerTrkDist       [fMaxNeles];
	float fTElConvPartnerTrkDCot       [fMaxNeles];
	float fTElMT                       [fMaxNeles];
	int   fTElGenID                    [fMaxNeles];
	int   fTElGenMID                   [fMaxNeles];
	int   fTElGenGMID                  [fMaxNeles];
	int   fTElGenType                  [fMaxNeles];
	int   fTElGenMType                 [fMaxNeles];
	int   fTElGenGMType                [fMaxNeles];
	float fTElHybRelIso                [fMaxNeles];
};
#endif
