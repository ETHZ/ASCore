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
	
	void           Begin(const char* filename = "SSDLTree.root");
	void           Analyze();
	void           End();
	
	void           BookTree();
	void           BookRunAndTriggerVariables(TTree* tree);
	void           BookMuonVariables(TTree* tree);
	void           BookElectronVariables(TTree* tree);
	void           BookJetMETVariables(TTree* tree);
	void           BookFPRVariables(TTree* tree);
	
	void           ResetTree();
	void           ResetRunAndTriggerVariables();
	void           ResetMuonVariables();
	void           ResetElectronVariables();
	void           ResetJetMETVariables();
	void           ResetFPRVariables();
	
	TLorentzVector jetTotalP(vector<int>& qualJetInd);
	TLorentzVector elTotalP (vector<int>& qualElInd);
	TLorentzVector muTotalP (vector<int>& qualMuInd);
	TLorentzVector phoTotalP(vector<int>& qualPhoInd);
	float          jetHT(vector<int>& qualJetInd);	
	float          minDRtoJet(float lepEta, float lepPhi);
	void           transverseMasses(TLorentzVector p1, TLorentzVector p2, TVector3 jtotPT, float &lepminv, float &lepmtinv, float &lepmct, float &lepmctort, float &lepmctparl, float &lepmt2_0, float &lepmt2_50, float &lepmt2_100, float &lepmt2orth_0, float &lepmt2orth_50, float &lepmt2orth_100);
	void           transverseAlphas(vector<int> qualElInd, vector<int> qualMuInd, vector<int> qualPhoInd, vector<int> qualJetInd, float &alphaT_h, float &alphaCT_h, float &alphaT, float &alphaCT, float &alphaT_new, float &alphaCT_new);
	
	void           DumpRunAndTiggerProperties();
	void           DumpJetMETProperties  (vector<int>& selectedJetInd);
	void           DumpPhotonProperties  (vector<int>& selectedPhoInd);
	void           DumpMuonProperties    (vector<int>& selectedMuInd);
	void           DumpElectronProperties(vector<int>& selectedElInd, TVector3 jtotPT);
	void           DumpFPRatioProperties();
	void           DumpElectronLooseAndTighPtAndEta(int elindex, float &elLoosePt, float &elTightPt, float &elLooseEta, float &elTightEta);
	void           DumpTwoElectronPtAndEta(int el1index, int el2index, float &el1Pt, float &el2Pt, float &el1Eta, float &el2Eta);
	
private:
	static const int fMaxNjets = 30;
	static const int fMaxNmus  = 5;
	static const int fMaxNeles = 5;
	static const int fMaxNphos = 5;
	
	TTree* fAnalysisTree;
	
	// run/sample properties
	int   fTRunNumber;
	int   fTEventNumber;
	int   fTLumiSection;
	float fTextxslo;
	float fTintxs;
	// trigger properties
	int   fTHLT_Mu9;
	int   fTHLT_Mu11;
	int   fTHLT_Mu13_v1;
	int   fTHLT_Mu15;
	int   fTHLT_Mu15_v1;
	int   fTHLT_DoubleMu0;
	int   fTHLT_DoubleMu3;
	int   fTHLT_DoubleMu3_v2;
	int   fTHLT_DoubleMu5_v2;
	int   fTHLT_Jet15U;
	int   fTHLT_Jet30U;
	int   fTHLT_Jet50U;
	int   fTHLT_Jet70U;
	int   fTHLT_Jet100U;
	int   fTHLT_Jet100U_v2;
	int   fTHLT_Jet100U_v3;
	int   fTHLT_HT100U;
	int   fTHLT_HT120U;
	int   fTHLT_HT140U;
	int   fTHLT_HT150U;
	int   fTHLT_HT150U_v3;
	int   fTHLT_HT130U;
	int   fTHLT_HT200U;
	int   fTHLT_Mu5_HT70U_v1;
	int   fTHLT_Mu5_HT100U_v1;
	int   fTHLT_DoubleMu3_HT50U;
	int   fTHLT_Ele10_HT70U;
	int   fTHLT_Ele10_HT100U;
	int   fTHLT_Ele10_EleId_HT70U;
	int   fTHLT_DoubleEle8_SW_HT70U_LR1_v1;
	int   fTHLT_Mu5_Ele5;
	int   fTHLT_Mu3_Ele8_HT70U_v1;
	int   fTHLT_Ele10_LW_L1R;
	int   fTHLT_Ele10_SW_L1R;
	int   fTHLT_Ele15_LW_L1R;
	int   fTHLT_Ele15_SW_L1R;
	int   fTHLT_Ele20_SW_L1R;
	int   fTHLT_DoubleEle5_SW_L1R;
	int   fTHLT_DoubleEle10_SW_L1R;
	int   fTHLT_DoubleEle15_SW_L1R_v1;
	int   fTHLT_DoubleEle17_SW_L1R_v1;
	int   fTHLT_Ele10_LW_EleId_L1R;
	int   fTHLT_Ele10_SW_EleId_L1R;
	int   fTHLT_Ele15_SW_CaloEleId_L1R;
	int   fTHLT_Ele15_SW_EleId_L1R;
	int   fTHLT_Ele17_SW_LooseEleId_L1R;
	int   fTHLT_Ele17_SW_TightEleId_L1R;
	int   fTHLT_Ele17_SW_CaloEleId_L1R;
	int   fTHLT_Ele17_SW_EleId_L1R;
	int   fTHLT_Ele17_SW_TightCaloEleId_SC8HE_L1R_v1;
	int   fTHLT_Ele17_SW_TightEleIdIsol_L1R_v1;
	int   fTHLT_Ele17_SW_TighterEleId_L1R_v1;
	int   fTHLT_Ele22_SW_TighterEleId_L1R_v2;
	int   fTHLT_Ele22_SW_TighterEleId_L1R_v3;
	int   fTHLT_Ele27_SW_TightCaloEleIdTrack_L1R_v1;
	int   fTHLT_Ele32_SW_TightCaloEleIdTrack_L1R_v1;
	int   fTHLT_Ele32_SW_TighterEleId_L1R_v2;
	int   fTHLT_GoodElEvent;
	int   fTHLT_GoodElEvent_RA5;
	int   fTHLT_GoodElEvent_TDL;
	int   fTHLT_GoodElFakesEvent;
	int   fTHLT_GoodHadronicEvent;
	int   fTHLT_GoodMuEvent;
	
	// jet-MET properties
	int   fTnqjets;
	float fTJetpt [fMaxNjets];
	float fTJeteta[fMaxNjets];
	float fTJetphi[fMaxNjets];
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
	int   fTnqphos;
	float fTPhopt          [fMaxNphos];
	float fTPhoeta         [fMaxNphos];
	float fTPhophi         [fMaxNphos];
	float fTPhoRelIso      [fMaxNphos];
	float fTPhoDRjet       [fMaxNphos];
	float fTPhoDRhardestjet[fMaxNphos];
	
	//muon properties
	int   fTnqmus;
	float fTmupt          [fMaxNmus];
	float fTmueta         [fMaxNmus];
	float fTmuphi         [fMaxNmus];
	float fTmuiso         [fMaxNmus];
	float fTmuisohyb      [fMaxNmus];
	int   fTmucharge      [fMaxNmus];
	int   fTmutight       [fMaxNmus]; // 0 for loose (but not tight), 1 for tight
	float fTmuDRjet       [fMaxNmus];
	float fTmuDRhardestjet[fMaxNmus];
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
	float fTmuMT;
	float fTmuMinv;
	
	// electron properties
	int     fTnqels;
	int   fTElcharge                   [fMaxNeles];
	int   fTElChargeIsCons             [fMaxNeles];
	int   fTElChargeIsGenCons          [fMaxNeles];
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
	float fTElDRjet                    [fMaxNeles];
	float fTElDRhardestjet             [fMaxNeles];
	int   fTElGenID                    [fMaxNeles];
	int   fTElGenMID                   [fMaxNeles];
	int   fTElGenGMID                  [fMaxNeles];
	int   fTElGenType                  [fMaxNeles];
	int   fTElGenMType                 [fMaxNeles];
	int   fTElGenGMType                [fMaxNeles];
	int   fTElTight                    [fMaxNeles];
	float fTElHybRelIso                [fMaxNeles];
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
