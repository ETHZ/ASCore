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

using namespace std;

class SSDLAnalysis : public UserAnalysisBase{
public:
	SSDLAnalysis(TreeReader *tr = 0);
	virtual ~SSDLAnalysis();
	
	void           Begin(const char* filename = "SSDLTree.root");
	void           Analyze();
	void           End();
	
	void           BookTree();	
	void           ResetTree();

private:
	static const int fMaxNjets = 30;
	static const int fMaxNmus  = 5;
	static const int fMaxNeles = 5;
	
	static const int gMaxhltbits = 300;
	
	TTree* fAnalysisTree;
	
	// run/sample properties
	int   fTRunNumber;
	int   fTEventNumber;
	int   fTLumiSection;
	// trigger properties
	int fTHLTNPaths;
	int fTHLTres     [gMaxhltbits];
	int fTHLTprescale[gMaxhltbits];
	std::vector<std::string> fTHLTnames;

	// jet-MET properties
	int   fTnqjets;
	float fTJetpt [fMaxNjets];
	float fTJeteta[fMaxNjets];
	float fTJetphi[fMaxNjets];
	float fTtcMET;
	float fTpfMET;
	
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
