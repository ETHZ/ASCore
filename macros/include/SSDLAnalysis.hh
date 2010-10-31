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
#include "helper/Davismt2.h"
#include "helper/Utilities.hh"
#include "base/UserAnalysisBase.hh"

using namespace std;

class SSDLAnalysis : public UserAnalysisBase{
public:
	SSDLAnalysis(TreeReader *tr = 0);
	virtual ~SSDLAnalysis();
	
	Davismt2 *fMT2;
	
	void Begin(const char* filename = "SSDLTree.root");
	void Analyze();
	void End();

	void BookTree();
	void BookTriggerVariables(TTree* tree);
	void BookMuonVariables(TTree* tree);
	void BookElectronVariables(TTree* tree);
	void BookAuxVariables(TTree* tree);
	void BookFPRVariables(TTree* tree);
	void ResetTree();
	void ResetTriggerVariables();
	void ResetMuonVariables();
	void ResetElectronVariables();
	void ResetAuxVariables();
	void ResetFPRVariables();
	
	TLorentzVector jetTotalP(vector<int>& qualJetInd);
	TLorentzVector elTotalP(vector<int>& qualElInd);
	TLorentzVector muTotalP(vector<int>& qualMuInd);
	TLorentzVector phoTotalP(vector<int>& qualPhoInd);
	float jetHT(vector<int>& qualJetInd);	
	float minDRtoJet(float lepEta, float lepPhi);
	void transverseMasses(TLorentzVector p1, TLorentzVector p2, TVector3 jtotPT, float &lepminv, float &lepmtinv, float &lepmct, float &lepmctort, float &lepmctparl, float &lepmt2_0, float &lepmt2_50, float &lepmt2_100, float &lepmt2orth_0, float &lepmt2orth_50, float &lepmt2orth_100);
	void transverseAlphas(vector<int> qualElInd, vector<int> qualMuInd, vector<int> qualPhoInd, vector<int> qualJetInd, float &alphaT_h, float &alphaCT_h, float &alphaT, float &alphaCT, float &alphaT_new, float &alphaCT_new);

	void DumpRunAndTiggerProperties();
	void DumpJetMETProperties(vector<int>& selectedJetInd);
	void DumpPhotonProperties(vector<int>& selectedPhoInd);
	void DumpMuonProperties(vector<int>& selectedMuInd);
	void DumpElectronProperties(vector<int>& selectedElInd, TVector3 jtotPT);
	void DumpFPRatioProperties();
	void DumpElectronLooseAndTighPtAndEta(int elindex, float &elLoosePt, float &elTightPt, float &elLooseEta, float &elTightEta);
	void DumpTwoElectronPtAndEta(int el1index, int el2index, float &el1Pt, float &el2Pt, float &el1Eta, float &el2Eta);

private:
	static const int gMaxNjets = 30;
	static const int gMaxNmus  = 5;
	static const int gMaxNeles = 5;
	static const int gMaxNphos = 5;

	TTree* fAnalysisTree;
	
	// run/sample properties
	int fTRunNumber;
	int fTEventNumber;
	int fTLumiSection;
	float fTextxslo;
	float fTintxs;

	// trigger properties
	int fT_HLTMu9;
	int fT_HLTMu11;
	int fT_HLTMu15;
	int fT_HLTDoubleMu3;
	int fT_HLTDoubleMu0;
	int fTHLT_Jet15U;
	int fTHLT_Jet30U;
	int fTHLT_Jet50U;
	int fTHLT_Jet70U;
	int fTHLT_Jet100U;
	int fTHLT_HT100U;
	int fTHLT_HT120U;
	int fTHLT_HT140U;
	int fTHLT_HT150U;
	int fTHLT_Ele10_LW_L1R;
	int fTHLT_Ele10_SW_L1R;
	int fTHLT_Ele15_LW_L1R;
	int fTHLT_Ele15_SW_L1R;
	int fTHLT_Ele15_SW_CaloEleId_L1R;
	int fTHLT_Ele20_SW_L1R;
	int fTHLT_DoubleEle5_SW_L1R;
	int fTHLT_DoubleEle10_SW_L1R;
	int fTHLT_DoubleEle15_SW_L1R_v1;
	int fTHTL_GoodElEvent;
	int fTHTL_GoodElFakesEvent;
	int fTHTL_GoodHadronicEvent;
	int fTHTL_GoodMuEvent;

	// jet-MET properties
	int fTnqjets;
	float fTJetpt [gMaxNjets];
	float fTJeteta[gMaxNjets];
	float fTJetphi[gMaxNjets];
	float fTHT;
	float fTSumEt;
	float fTtcMET;
	float fTpfMET;
	float fTMuCorrMET;
	float fTdPhiMJ1;
	float fTdPhiMJ2;
	float fTR12;
	float fTR21;
	float fTR12plusR21;

	// photon properties
	int fTnqphos;
	float fTPhopt          [gMaxNphos];
	float fTPhoeta         [gMaxNphos];
	float fTPhophi         [gMaxNphos];
	float fTPhoRelIso      [gMaxNphos];
	float fTPhoDRjet       [gMaxNphos];
	float fTPhoDRhardestjet[gMaxNphos];
	
	//muon properties
	int fTnqmus;
	float fTmupt          [gMaxNmus];
	float fTmueta         [gMaxNmus];
	float fTmuphi         [gMaxNmus];
	float fTmuiso         [gMaxNmus];
	float fTmuisohyb      [gMaxNmus];
	int fTmucharge        [gMaxNmus];
	int fTmutight         [gMaxNmus]; // 0 for loose (but not tight), 1 for tight
	float fTmuDRjet       [gMaxNmus];
	float fTmuDRhardestjet[gMaxNmus];
	float fTmud0          [gMaxNmus];
	float fTmudz          [gMaxNmus];
	float fTmud0bs        [gMaxNmus];
	float fTmudzbs        [gMaxNmus];
	float fTmuptE         [gMaxNmus];
	int fTmuid            [gMaxNmus];
	int fTmumoid          [gMaxNmus];
	int fTmugmoid         [gMaxNmus];
	int fTmutype          [gMaxNmus];
	int fTmumotype        [gMaxNmus];
	int fTmugmotype       [gMaxNmus];
	float fTmuMT;
	float fTmuMinv;
	
	// electron properties
	int fTnqels;
	int fTElcharge                     [gMaxNeles];
	int fTElChargeIsCons               [gMaxNeles];
	int fTElChargeIsGenCons            [gMaxNeles];
	float fTElpt                       [gMaxNeles];
	float fTEleta                      [gMaxNeles];
	float fTElphi                      [gMaxNeles];
	float fTEld0                       [gMaxNeles];
	float fTElD0Err                    [gMaxNeles];
	float fTElEoverP                   [gMaxNeles];
	float fTElHoverE                   [gMaxNeles];
	float fTElSigmaIetaIeta            [gMaxNeles];
	float fTElDeltaPhiSuperClusterAtVtx[gMaxNeles];
	float fTElDeltaEtaSuperClusterAtVtx[gMaxNeles];
	float fTElIDsimpleWP80relIso       [gMaxNeles];
	float fTElIDsimpleWPrelIso         [gMaxNeles];
	float fTElIDsimpleWP95relIso       [gMaxNeles];	
	float fTElRelIso                   [gMaxNeles];
	float fTElDR04TkSumPt              [gMaxNeles];
	float fTElDR04EcalRecHitSumEt      [gMaxNeles];
	float fTElDR04HcalTowerSumEt       [gMaxNeles];
	float fTElS4OverS1                 [gMaxNeles];
	float fTElConvPartnerTrkDist       [gMaxNeles];
	float fTElConvPartnerTrkDCot       [gMaxNeles];
	float fTElChargeMisIDProb          [gMaxNeles];
	float fTElMT                       [gMaxNeles];
	float fTElDRjet                    [gMaxNeles];
	float fTElDRhardestjet             [gMaxNeles];
	int fTElGenID                      [gMaxNeles];
	int fTElGenStatus                  [gMaxNeles];
	int fTElGenMID                     [gMaxNeles];
	int fTElGenMStatus                 [gMaxNeles];
	int fTElGenGMID                    [gMaxNeles];
	int fTElGenGMStatus                [gMaxNeles];
	int fTElTight                      [gMaxNeles];
	float fTElHybRelIso                [gMaxNeles];
	float fTElminv;
	float fTElmtinv;

	// other properties
	float fTElmt2_0;
	float fTElmt2_50;
	float fTElmt2_100;
	float fTElmCT;
	float fTElmCTorth;
	float fTElmCTparl;
	float fTElmT2orth_0;
	float fTElmT2orth_50;
	float fTElmT2orth_100;	
	float fTMumt2_0;
	float fTMumt2_50;
	float fTMumt2_100;
	float fTMumCT;
	float fTMumCTorth;
	float fTMumCTparl;
	float fTMumT2orth_0;
	float fTMumT2orth_50;
	float fTMumT2orth_100;
	float fTalphaT_h;
	float fTalphaCT_h;	
	float fTalphaT;
	float fTalphaCT;	
	float fTalphaT_new;
	float fTalphaCT_new;

	// fake ratio propeties
	int		fTisSE_QCDLike;
	int		fTSE_QCDLike_FakeElGenID;
	float	fTSE_QCDLike_ElLoosePt;
	float	fTSE_QCDLike_ElTightPt;
	float	fTSE_QCDLike_ElLooseEta;
	float	fTSE_QCDLike_ElTightEta;
	
	int		fTisSE_AntiQCDLike;
	int		fTSE_AntiQCDLike_FakeElGenID;
	float	fTSE_AntiQCDLike_ElLoosePt;
	float	fTSE_AntiQCDLike_ElTightPt;
	float	fTSE_AntiQCDLike_ElLooseEta;
	float	fTSE_AntiQCDLike_ElTightEta;
	
	int		fTisDE_ZJetsLike;
	float	fTDE_ZJetsLike_ElLoosePt;
	float	fTDE_ZJetsLike_ElTightPt;
	float	fTDE_ZJetsLike_ElLooseEta;
	float	fTDE_ZJetsLike_ElTightEta;
	float	fTDE_ZJetsLike_PromptElGenLoosePt;
	float	fTDE_ZJetsLike_PromptElGenTightPt;
	float	fTDE_ZJetsLike_PromptElGenLooseEta;
	float	fTDE_ZJetsLike_PromptElGenTightEta;

	int		fTisDE_WJetsLike;
	int		fTDE_WJetsLike_FakeElGenID;
	float	fTDE_WJetsLike_ElLoosePt;
	float	fTDE_WJetsLike_ElTightPt;
	float	fTDE_WJetsLike_ElLooseEta;
	float	fTDE_WJetsLike_ElTightEta;
	float	fTDE_WJetsLike_ElTightMT;
	float	fTDE_WJetsLike_FakeElGenLoosePt;
	float	fTDE_WJetsLike_FakeElGenTightPt;
	float	fTDE_WJetsLike_FakeElGenLooseEta;
	float	fTDE_WJetsLike_FakeElGenTightEta;
	float	fTDE_WJetsLike_PromptElGenLoosePt;
	float	fTDE_WJetsLike_PromptElGenTightPt;
	float	fTDE_WJetsLike_PromptElGenLooseEta;
	float	fTDE_WJetsLike_PromptElGenTightEta;
	float	fTDE_WJetsLike_PromptElGenMT;
	
	int		fTisDE_AntiWJetsLike;
	int		fTDE_AntiWJetsLike_FakeElGenID;
	float	fTDE_AntiWJetsLike_ElLoosePt;
	float	fTDE_AntiWJetsLike_ElTightPt;
	float	fTDE_AntiWJetsLike_ElLooseEta;
	float	fTDE_AntiWJetsLike_ElTightEta;
	float	fTDE_AntiWJetsLike_ElTightMT;
	float	fTDE_AntiWJetsLike_FakeElGenLoosePt;
	float	fTDE_AntiWJetsLike_FakeElGenTightPt;
	float	fTDE_AntiWJetsLike_FakeElGenLooseEta;
	float	fTDE_AntiWJetsLike_FakeElGenTightEta;
	float	fTDE_AntiWJetsLike_PromptElGenLoosePt;
	float	fTDE_AntiWJetsLike_PromptElGenTightPt;
	float	fTDE_AntiWJetsLike_PromptElGenLooseEta;
	float	fTDE_AntiWJetsLike_PromptElGenTightEta;
	float	fTDE_AntiWJetsLike_PromptElGenMT;
	
	int		fTisDE_SignalLike;
	float	fTDE_Ntt_El1Pt;
	float	fTDE_Ntt_El2Pt;
	float	fTDE_Ntl_El1Pt;
	float	fTDE_Ntl_El2Pt;
	float	fTDE_Nlt_El1Pt;
	float	fTDE_Nlt_El2Pt;
	float	fTDE_Nll_El1Pt;
	float	fTDE_Nll_El2Pt;
	float	fTDE_Ntt_El1Eta;
	float	fTDE_Ntt_El2Eta;
	float	fTDE_Ntl_El1Eta;
	float	fTDE_Ntl_El2Eta;
	float	fTDE_Nlt_El1Eta;
	float	fTDE_Nlt_El2Eta;
	float	fTDE_Nll_El1Eta;
	float	fTDE_Nll_El2Eta;
};
#endif
