#include "SSDLAnalysis.hh"

using namespace std;

const int SSDLAnalysis::gMaxNjets;
const int SSDLAnalysis::gMaxNmus;
const int SSDLAnalysis::gMaxNeles;
const int SSDLAnalysis::gMaxNphos;

SSDLAnalysis::SSDLAnalysis(TreeReader *tr): UserAnalysisBase(tr){
	//SetStyle();	
}

SSDLAnalysis::~SSDLAnalysis(){
}

void SSDLAnalysis::Begin(const char* filename){
	ReadPDGTable();
	BookTree();
}

void SSDLAnalysis::End(){
	fOutputFile->cd();
	fAnalysisTree->Write();
	fOutputFile->Close();
}


void SSDLAnalysis::BookTree(){
	fOutputFile->cd();
	fAnalysisTree = new TTree("Analysis", "AnalysisTree");
	BookTriggerVariables  (fAnalysisTree);
	BookMuonVariables     (fAnalysisTree);
	BookElectronVariables (fAnalysisTree);
	BookAuxVariables      (fAnalysisTree);
	BookFPRVariables      (fAnalysisTree);
}

void SSDLAnalysis::BookTriggerVariables(TTree* tree){
	tree->Branch("HLT_Mu9",                    &fT_HLTMu9,                    "HLT_Mu9/I");
	tree->Branch("HLT_Mu11",                   &fT_HLTMu11,                   "HLT_Mu11/I");
	tree->Branch("HLT_Mu15",                   &fT_HLTMu15,                   "HLT_Mu15/I");
	tree->Branch("HLT_DoubleMu0",              &fT_HLTDoubleMu0,              "HLT_DoubleMu0/I");
	tree->Branch("HLT_DoubleMu3",              &fT_HLTDoubleMu3,              "HLT_DoubleMu3/I");
	tree->Branch("HLT_Jet15U",                 &fTHLT_Jet15U,                 "HLT_Jet15U/I");
	tree->Branch("HLT_Jet30U",                 &fTHLT_Jet30U,                 "HLT_Jet30U/I");
	tree->Branch("HLT_Jet50U",                 &fTHLT_Jet50U,                 "HLT_Jet50U/I");
	tree->Branch("HLT_Jet70U",                 &fTHLT_Jet70U,                 "HLT_Jet70U/I");
	tree->Branch("HLT_Jet100U",                &fTHLT_Jet100U,                "HLT_Jet100U/I");
	tree->Branch("HLT_HT100U",                 &fTHLT_HT100U,                 "HLT_HT100U/I");
	tree->Branch("HLT_HT120U",                 &fTHLT_HT120U,                 "HLT_HT120U/I");
	tree->Branch("HLT_HT140U",                 &fTHLT_HT140U,                 "HLT_HT140U/I");
	tree->Branch("HLT_HT150U",                 &fTHLT_HT150U,                 "HLT_HT150U/I");
	tree->Branch("HLT_Ele10_LW_L1R",           &fTHLT_Ele10_LW_L1R,           "HLT_Ele10_LW_L1R/I");
	tree->Branch("HLT_Ele10_SW_L1R",           &fTHLT_Ele10_SW_L1R,           "HLT_Ele10_SW_L1R/I");
	tree->Branch("HLT_Ele15_LW_L1R",           &fTHLT_Ele15_LW_L1R,           "HLT_Ele15_LW_L1R/I");
	tree->Branch("HLT_Ele15_SW_L1R",           &fTHLT_Ele15_SW_L1R,           "HLT_Ele15_SW_L1R/I");
	tree->Branch("HLT_Ele15_SW_CaloEleId_L1R", &fTHLT_Ele15_SW_CaloEleId_L1R, "HLT_Ele15_SW_CaloEleId_L1R/I");
	tree->Branch("HLT_Ele20_SW_L1R",           &fTHLT_Ele20_SW_L1R,           "HLT_Ele20_SW_L1R/I");
	tree->Branch("HLT_DoubleEle5_SW_L1R",      &fTHLT_DoubleEle5_SW_L1R,      "HLT_DoubleEle5_SW_L1R/I");
	tree->Branch("HLT_DoubleEle10_SW_L1R",     &fTHLT_DoubleEle10_SW_L1R,     "HLT_DoubleEle10_SW_L1R/I");
	tree->Branch("HLT_DoubleEle15_SW_L1R_v1",  &fTHLT_DoubleEle15_SW_L1R_v1,  "HLT_DoubleEle15_SW_L1R_v1/I");
	tree->Branch("HTL_GoodElEvent",            &fTHTL_GoodElEvent,            "HTL_GoodElEvent/I");
	tree->Branch("HTL_GoodElFakesEvent",       &fTHTL_GoodElFakesEvent,       "HTL_GoodElFakesEvent/I");
	tree->Branch("HTL_GoodMuEvent",            &fTHTL_GoodMuEvent,            "HTL_GoodMuEvent/I");
	tree->Branch("HTL_GoodHadronicEvent",      &fTHTL_GoodHadronicEvent,      "HTL_GoodHadronicEvent/I");
}

void SSDLAnalysis::BookMuonVariables(TTree* tree){
	// single-muon properties
	tree->Branch("NMus"          ,&fTnqmus,          "NMus/I");
	tree->Branch("MuPt"          ,&fTmupt,           "MuPt[NMus]/F");
	tree->Branch("MuEta"         ,&fTmueta,          "MuEta[NMus]/F");
	tree->Branch("MuPhi"         ,&fTmuphi,          "MuPhi[NMus]/F");
	tree->Branch("MuCharge"      ,&fTmucharge,       "MuCharge[NMus]/I");
	tree->Branch("MuTight"       ,&fTmutight,        "MuTight[NMus]/I");
	tree->Branch("MuIso"         ,&fTmuiso,          "MuIso[NMus]/F");
	tree->Branch("MuIsoHybrid"   ,&fTmuisohyb,       "MuIsoHybrid[NMus]/F");
	tree->Branch("MuDRJet"       ,&fTmuDRjet,        "MuDRJet[NMus]/F");
	tree->Branch("MuDRHardestJet",&fTmuDRhardestjet, "MuDRHardestJet[NMus]/F");
	tree->Branch("MuD0"          ,&fTmud0,           "MuD0[NMus]/F");
	tree->Branch("MuDz"          ,&fTmudz,           "MuDz[NMus]/F");
	tree->Branch("MuD0BS"        ,&fTmud0bs,         "MuD0BS[NMus]/F");
	tree->Branch("MuDzBS"        ,&fTmudzbs,         "MuDzBS[NMus]/F");
	tree->Branch("MuPtE"         ,&fTmuptE,          "MuPtE[NMus]/F");
	tree->Branch("MuGenID"       ,&fTmuid,           "MuGenID[NMus]/I");
	tree->Branch("MuGenMoID"     ,&fTmumoid,         "MuGenMoID[NMus]/I");
	tree->Branch("MuGenGMoID"    ,&fTmugmoid,        "MuGenGMoID[NMus]/I");
	tree->Branch("MuGenType"     ,&fTmutype,         "MuGenType[NMus]/I");
	tree->Branch("MuGenMoType"   ,&fTmumotype,       "MuGenMoType[NMus]/I");
	tree->Branch("MuGenGMoType"  ,&fTmugmotype,      "MuGenGMoType[NMus]/I");
	tree->Branch("MuMT"          ,&fTmuMT,           "MuMT/F");
	tree->Branch("MuMinv"        ,&fTmuMinv,         "MuMinv/F");
}

void SSDLAnalysis::BookElectronVariables(TTree* tree){
	// single-electron properties
	tree->Branch("NEls",                   &fTnqels,               "NEls/I");
	tree->Branch("ElCh",                   &fTElcharge,            "ElCh[NEls]/I");
	tree->Branch("ElChIsCons",             &fTElChargeIsCons,      "ElChIsCons[NEls]/I");
	tree->Branch("ElChIsGenCons",          &fTElChargeIsGenCons,   "ElChIsGenCons[NEls]/I");
	tree->Branch("ElPt",                   &fTElpt,                "ElPt[NEls]/F");
	tree->Branch("ElEta",                  &fTEleta,               "ElEta[NEls]/F");
	tree->Branch("ElPhi",                  &fTElphi,               "ElPhi[NEls]/F");
	tree->Branch("ElD0",                   &fTEld0,                "ElD0[NEls]/F");
	tree->Branch("ElD0Err",                &fTElD0Err,             "ElD0Err[NEls]/F");
	tree->Branch("ElEoverP",               &fTElEoverP,            "ElEoverP[NEls]/F");
	tree->Branch("ElHoverE",               &fTElHoverE,            "ElHoverE[NEls]/F");
	tree->Branch("ElSigmaIetaIeta",        &fTElSigmaIetaIeta,     "ElSigmaIetaIeta[NEls]/F");
	tree->Branch("ElDeltaPhiSuperClusterAtVtx", &fTElDeltaPhiSuperClusterAtVtx, "ElDeltaPhiSuperClusterAtVtx[NEls]/F");
	tree->Branch("ElDeltaEtaSuperClusterAtVtx", &fTElDeltaEtaSuperClusterAtVtx, "ElDeltaEtaSuperClusterAtVtx[NEls]/F");
	tree->Branch("ElRelIso",               &fTElRelIso,            "ElRelIso[NEls]/F");
	tree->Branch("ElDR04TkSumPt",          &fTElDR04TkSumPt,       "ElDR04TkSumPt[NEls]/F");
	tree->Branch("ElDR04EcalRecHitSumEt",  &fTElDR04EcalRecHitSumEt, "ElDR04EcalRecHitSumEt[NEls]/F");
	tree->Branch("ElDR04HcalTowerSumEt",   &fTElDR04HcalTowerSumEt,"ElDR04HcalTowerSumEt[NEls]/F");
	tree->Branch("ElS4OverS1",             &fTElS4OverS1,          "ElS4OverS1[NEls]/F");
	tree->Branch("ElConvPartnerTrkDist",   &fTElConvPartnerTrkDist,"ElConvPartnerTrkDist[NEls]/F");
	tree->Branch("ElConvPartnerTrkDCot",   &fTElConvPartnerTrkDCot,"ElConvPartnerTrkDCot[NEls]/F");
	tree->Branch("ElChargeMisIDProb",      &fTElChargeMisIDProb,   "ElChargeMisIDProb[NEls]/F");
	tree->Branch("ElDRjet",                &fTElDRjet,             "ElDRjet[NEls]/F");
	tree->Branch("ElDRhardestjet",         &fTElDRhardestjet,      "ElDRhardestjet[NEls]/F");
	tree->Branch("ElGenID",                &fTElGenID,             "ElGenID[NEls]/I");
	tree->Branch("ElGenStatus",            &fTElGenStatus,         "ElGenStatus[NEls]/I");
	tree->Branch("ElGenMID",               &fTElGenMID,            "ElGenMID[NEls]/I");
	tree->Branch("ElGenMStatus",           &fTElGenMStatus,        "ElGenMStatus[NEls]/I");
	tree->Branch("ElGenGMID",              &fTElGenGMID,           "ElGenGMID[NEls]/I");
	tree->Branch("ElGenGMStatus",          &fTElGenGMStatus,       "ElGenGMStatus[NEls]/I");
	tree->Branch("ElTight",                &fTElTight,             "ElTight[NEls]/I");
	tree->Branch("ElHybRelIso",            &fTElHybRelIso,         "ElHybRelIso[NEls]/F");
	tree->Branch("ElMT",                   &fTElMT,                "ElMT[NEls]/F");
	tree->Branch("ElMInv",                 &fTElminv,              "ElMInv/F");
	tree->Branch("ElMTInv",                &fTElmtinv,             "ElMTInv/F");
}

void SSDLAnalysis::BookAuxVariables(TTree* tree){
	// run/sample properties
	tree->Branch("Run",             &fTRunNumber,   "Run/I");
	tree->Branch("Event",           &fTEventNumber, "Event/I");
	tree->Branch("LumiSec",         &fTLumiSection, "LumiSec/I");
	tree->Branch("ExtXSecLO",       &fTextxslo,     "ExtXSecLO/F");
	tree->Branch("IntXSec",         &fTintxs,       "IntXSec/F");
	// event properties
	tree->Branch("HT",              &fTHT,          "HT/F");
	tree->Branch("SumEt",           &fTSumEt,       "SumEt/F");
	tree->Branch("tcMET",           &fTtcMET,       "tcMET/F");
	tree->Branch("pfMET",           &fTpfMET,       "pfMET/F");
	tree->Branch("MuCorrMET",       &fTMuCorrMET,   "MuCorrMET/F");
	// jet-MET properties
	tree->Branch("NJets",           &fTnqjets,      "NJets/I");
	tree->Branch("JetPt",           &fTJetpt,       "JetPt[NJets]/F");
	tree->Branch("JetEta",          &fTJeteta,      "JetEta[NJets]/F");
	tree->Branch("JetPhi",          &fTJetphi,      "JetPhi[NJets]/F");
	tree->Branch("dPhiMJ1",         &fTdPhiMJ1,     "dPhiMJ1/F");
	tree->Branch("dPhiMJ2",         &fTdPhiMJ2,     "dPhiMJ2/F");
	tree->Branch("R12",             &fTR12,         "R12/F");
	tree->Branch("R21",             &fTR21,         "R21/F");
	tree->Branch("R12plusR21",      &fTR12plusR21,  "R12plusR21/F");
	// photon properties
	tree->Branch("NPhos",           &fTnqphos,      "NPhos/I");
	tree->Branch("PhoPt",           &fTPhopt,       "PhoPt[NPhos]/F");
	tree->Branch("PhoEta",          &fTPhoeta,      "PhoEta[NPhos]/F");
	tree->Branch("PhoPhi",          &fTPhophi,      "PhoPhi[NPhos]/F");
	tree->Branch("PhoRelIso",       &fTPhoRelIso,   "PhoRelIso[NPhos]/F");
	tree->Branch("PhoDRjet",        &fTPhoDRjet,    "PhoDRjet[NPhos]/F");
	tree->Branch("PhoDRhardestjet", &fTPhoDRhardestjet,"PhoDRhardestjet[NPhos]/F");
}

void SSDLAnalysis::BookFPRVariables(TTree* tree){
	// other properties
	tree->Branch("AlphaT_h",		&fTalphaT_h,		"AlphaT_h/F");
	tree->Branch("AlphaCT_h",		&fTalphaCT_h,		"AlphaCT_h/F");
	tree->Branch("AlphaT",			&fTalphaT,			"AlphaT/F");
	tree->Branch("AlphaCT",			&fTalphaCT,			"AlphaCT/F");
	tree->Branch("AlphaT_new",		&fTalphaT_new,		"AlphaT_new/F");
	tree->Branch("AlphaCT_new",		&fTalphaCT_new,		"AlphaCT_new/F");
	tree->Branch("ElMT2_0",			&fTElmt2_0,			"ElMT2_0/F");
	tree->Branch("ElMT2_50",		&fTElmt2_50,		"ElMT2_50/F");
	tree->Branch("ElMT2_100",		&fTElmt2_100,		"ElMT2_100/F");
	tree->Branch("ElMCT",			&fTElmCT,			"ElMCT/F");
	tree->Branch("ElMCTparl",		&fTElmCTparl,		"ElMCTparl/F");
	tree->Branch("ElMCTorth",		&fTElmCTorth,		"ElMCTorth/F");
	tree->Branch("ElMT2orth_0",		&fTElmT2orth_0,		"ElMT2orth_0/F");
	tree->Branch("ElMT2orth_50",	&fTElmT2orth_50,	"ElMT2orth_50/F");
	tree->Branch("ElMT2orth_100",	&fTElmT2orth_100,	"ElMT2orth_100/F");
	tree->Branch("MuMT2_0",			&fTMumt2_0,			"MuMT2_0/F");
	tree->Branch("MuMT2_50",		&fTMumt2_50,		"MuMT2_50/F");
	tree->Branch("MuMT2_100",		&fTMumt2_100,		"MuMT2_100/F");
	tree->Branch("MuMCT",			&fTMumCT,			"MuMCT/F");	
	tree->Branch("MuMCTparl",		&fTMumCTparl,		"MuMCTparl/F");
	tree->Branch("MuMCTorth",		&fTMumCTorth,		"MuMCTorth/F");
	tree->Branch("MuMT2orth_0",		&fTMumT2orth_0,		"MuMT2orth_0/F");
	tree->Branch("MuMT2orth_50",	&fTMumT2orth_50,	"MuMT2orth_50/F");
	tree->Branch("MuMT2orth_100",	&fTMumT2orth_100,	"MuMT2orth_100/F");
	
	// fake ratio propeties
	tree->Branch("isSE_QCDLike",				&fTisSE_QCDLike,				"isSE_QCDLike/I");
	tree->Branch("SE_QCDLike_FakeElGenID",		&fTSE_QCDLike_FakeElGenID,		"SE_QCDLike_FakeElGenID/I");
	tree->Branch("SE_QCDLike_ElLoosePt",		&fTSE_QCDLike_ElLoosePt,		"SE_QCDLike_ElLoosePt/F");
	tree->Branch("SE_QCDLike_ElTightPt",		&fTSE_QCDLike_ElTightPt,		"SE_QCDLike_ElTightPt/F");
	tree->Branch("SE_QCDLike_ElLooseEta",		&fTSE_QCDLike_ElLooseEta,		"SE_QCDLike_ElLooseEta/F");
	tree->Branch("SE_QCDLike_ElTightEta",		&fTSE_QCDLike_ElTightEta,		"SE_QCDLike_ElTightEta/F");
	
	tree->Branch("isSE_AntiQCDLike",			&fTisSE_AntiQCDLike,				"isSE_AntiQCDLike/I");
	tree->Branch("SE_AntiQCDLike_FakeElGenID",	&fTSE_AntiQCDLike_FakeElGenID,		"SE_AntiQCDLike_FakeElGenID/I");
	tree->Branch("SE_AntiQCDLike_ElLoosePt",	&fTSE_AntiQCDLike_ElLoosePt,		"SE_AntiQCDLike_ElLoosePt/F");
	tree->Branch("SE_AntiQCDLike_ElTightPt",	&fTSE_AntiQCDLike_ElTightPt,		"SE_AntiQCDLike_ElTightPt/F");
	tree->Branch("SE_AntiQCDLike_ElLooseEta",	&fTSE_AntiQCDLike_ElLooseEta,		"SE_AntiQCDLike_ElLooseEta/F");
	tree->Branch("SE_AntiQCDLike_ElTightEta",	&fTSE_AntiQCDLike_ElTightEta,		"SE_AntiQCDLike_ElTightEta/F");
	
	tree->Branch("isDE_ZJetsLike",				&fTisDE_ZJetsLike,			"isDE_ZJetsLike/I");
	tree->Branch("DE_ZJetsLike_ElLoosePt",		&fTDE_ZJetsLike_ElLoosePt,	"DE_ZJetsLike_ElLoosePt/F");
	tree->Branch("DE_ZJetsLike_ElTightPt",		&fTDE_ZJetsLike_ElTightPt,	"DE_ZJetsLike_ElTightPt/F");
	tree->Branch("DE_ZJetsLike_ElLooseEta",		&fTDE_ZJetsLike_ElLooseEta,	"DE_ZJetsLike_ElLooseEta/F");
	tree->Branch("DE_ZJetsLike_ElTightEta",		&fTDE_ZJetsLike_ElTightEta,	"DE_ZJetsLike_ElTightEta/F");
	tree->Branch("DE_ZJetsLike_PromptElGenLoosePt",	&fTDE_ZJetsLike_PromptElGenLoosePt,	"DE_ZJetsLike_PromptElGenLoosePt/F");
	tree->Branch("DE_ZJetsLike_PromptElGenTightPt",	&fTDE_ZJetsLike_PromptElGenTightPt,	"DE_ZJetsLike_PromptElGenTightPt/F");
	tree->Branch("DE_ZJetsLike_PromptElGenLooseEta",	&fTDE_ZJetsLike_PromptElGenLooseEta,	"DE_ZJetsLike_PromptElGenLooseEta/F");
	tree->Branch("DE_ZJetsLike_PromptElGenTightEta",	&fTDE_ZJetsLike_PromptElGenTightEta,	"DE_ZJetsLike_PromptElGenTightEta/F");
	
	tree->Branch("isDE_WJetsLike",					&fTisDE_WJetsLike,			"isDE_WJetsLike/I");
	tree->Branch("DE_WJetsLike_FakeElGenID",		&fTDE_WJetsLike_FakeElGenID,"DE_WJetsLike_FakeElGenID/I");
	tree->Branch("DE_WJetsLike_ElLoosePt",			&fTDE_WJetsLike_ElLoosePt,	"DE_WJetsLike_ElLoosePt/F");
	tree->Branch("DE_WJetsLike_ElTightPt",			&fTDE_WJetsLike_ElTightPt,	"DE_WJetsLike_ElTightPt/F");
	tree->Branch("DE_WJetsLike_ElLooseEta",			&fTDE_WJetsLike_ElLooseEta,	"DE_WJetsLike_ElLooseEta/F");
	tree->Branch("DE_WJetsLike_ElTightEta",			&fTDE_WJetsLike_ElTightEta,	"DE_WJetsLike_ElTightEta/F");
	tree->Branch("DE_WJetsLike_ElTightMT",			&fTDE_WJetsLike_ElTightMT,	"DE_WJetsLike_ElTightMT/F");
	tree->Branch("DE_WJetsLike_FakeElGenLoosePt",	&fTDE_WJetsLike_FakeElGenLoosePt,	"DE_WJetsLike_FakeElGenLoosePt/F");
	tree->Branch("DE_WJetsLike_FakeElGenTightPt",	&fTDE_WJetsLike_FakeElGenTightPt,	"DE_WJetsLike_FakeElGenTightPt/F");
	tree->Branch("DE_WJetsLike_FakeElGenLooseEta",	&fTDE_WJetsLike_FakeElGenLooseEta,	"DE_WJetsLike_FakeElGenLooseEta/F");
	tree->Branch("DE_WJetsLike_FakeElGenTightEta",	&fTDE_WJetsLike_FakeElGenTightEta,	"DE_WJetsLike_FakeElGenTightEta/F");
	tree->Branch("DE_WJetsLike_PromptElGenLoosePt",	&fTDE_WJetsLike_PromptElGenLoosePt,	"DE_WJetsLike_PromptElGenLoosePt/F");
	tree->Branch("DE_WJetsLike_PromptElGenTightPt",	&fTDE_WJetsLike_PromptElGenTightPt,	"DE_WJetsLike_PromptElGenTightPt/F");
	tree->Branch("DE_WJetsLike_PromptElGenLooseEta",	&fTDE_WJetsLike_PromptElGenLooseEta,	"DE_WJetsLike_PromptElGenLooseEta/F");
	tree->Branch("DE_WJetsLike_PromptElGenTightEta",	&fTDE_WJetsLike_PromptElGenTightEta,	"DE_WJetsLike_PromptElGenTightEta/F");
	tree->Branch("DE_WJetsLike_PromptElGenMT",			&fTDE_WJetsLike_PromptElGenMT,	"DE_WJetsLike_PromptElGenMT/F");
	
	tree->Branch("isDE_AntiWJetsLike",				&fTisDE_AntiWJetsLike,			"isDE_AntiWJetsLike/I");
	tree->Branch("DE_AntiWJetsLike_FakeElGenID",	&fTDE_AntiWJetsLike_FakeElGenID,"DE_AntiWJetsLike_FakeElGenID/I");
	tree->Branch("DE_AntiWJetsLike_ElLoosePt",		&fTDE_AntiWJetsLike_ElLoosePt,	"DE_AntiWJetsLike_ElLoosePt/F");
	tree->Branch("DE_AntiWJetsLike_ElTightPt",		&fTDE_AntiWJetsLike_ElTightPt,	"DE_AntiWJetsLike_ElTightPt/F");
	tree->Branch("DE_AntiWJetsLike_ElLooseEta",		&fTDE_AntiWJetsLike_ElLooseEta,	"DE_AntiWJetsLike_ElLooseEta/F");
	tree->Branch("DE_AntiWJetsLike_ElTightEta",		&fTDE_AntiWJetsLike_ElTightEta,	"DE_AntiWJetsLike_ElTightEta/F");
	tree->Branch("DE_AntiWJetsLike_ElTightMT",		&fTDE_AntiWJetsLike_ElTightMT,	"DE_AntiWJetsLike_ElTightMT/F");
	tree->Branch("DE_AntiWJetsLike_FakeElGenLoosePt",	&fTDE_AntiWJetsLike_FakeElGenLoosePt,	"DE_AntiWJetsLike_FakeElGenLoosePt/F");
	tree->Branch("DE_AntiWJetsLike_FakeElGenTightPt",	&fTDE_AntiWJetsLike_FakeElGenTightPt,	"DE_AntiWJetsLike_FakeElGenTightPt/F");
	tree->Branch("DE_AntiWJetsLike_FakeElGenLooseEta",	&fTDE_AntiWJetsLike_FakeElGenLooseEta,	"DE_AntiWJetsLike_FakeElGenLooseEta/F");
	tree->Branch("DE_AntiWJetsLike_FakeElGenTightEta",	&fTDE_AntiWJetsLike_FakeElGenTightEta,	"DE_AntiWJetsLike_FakeElGenTightEta/F");
	tree->Branch("DE_AntiWJetsLike_PromptElGenLoosePt",	&fTDE_AntiWJetsLike_PromptElGenLoosePt,	"DE_AntiWJetsLike_PromptElGenLoosePt/F");
	tree->Branch("DE_AntiWJetsLike_PromptElGenTightPt",	&fTDE_AntiWJetsLike_PromptElGenTightPt,	"DE_AntiWJetsLike_PromptElGenTightPt/F");
	tree->Branch("DE_AntiWJetsLike_PromptElGenLooseEta",	&fTDE_AntiWJetsLike_PromptElGenLooseEta,"DE_AntiWJetsLike_PromptElGenLooseEta/F");
	tree->Branch("DE_AntiWJetsLike_PromptElGenTightEta",	&fTDE_AntiWJetsLike_PromptElGenTightEta,"DE_AntiWJetsLike_PromptElGenTightEta/F");
	tree->Branch("DE_AntiWJetsLike_PromptElGenMT",			&fTDE_AntiWJetsLike_PromptElGenMT,		"DE_AntiWJetsLike_PromptElGenMT/F");
	
	tree->Branch("isDE_SignalLike",&fTisDE_SignalLike,		"isDE_SignalLike/I");
	tree->Branch("DE_Ntt_El1Pt",	&fTDE_Ntt_El1Pt,		"DE_Ntt_El1Pt/F");
	tree->Branch("DE_Ntt_El2Pt",	&fTDE_Ntt_El2Pt,		"DE_Ntt_El2Pt/F");
	tree->Branch("DE_Ntt_El1Eta",	&fTDE_Ntt_El1Eta,		"DE_Ntt_El1Eta/F");
	tree->Branch("DE_Ntt_El2Eta",	&fTDE_Ntt_El2Eta,		"DE_Ntt_El2Eta/F");
	tree->Branch("DE_Ntl_El1Pt",	&fTDE_Ntl_El1Pt,		"DE_Ntl_El1Pt/F");
	tree->Branch("DE_Ntl_El2Pt",	&fTDE_Ntl_El2Pt,		"DE_Ntl_El2Pt/F");
	tree->Branch("DE_Ntl_El1Eta",	&fTDE_Ntl_El1Eta,		"DE_Ntl_El1Eta/F");
	tree->Branch("DE_Ntl_El2Eta",	&fTDE_Ntl_El2Eta,		"DE_Ntl_El2Eta/F");
	tree->Branch("DE_Nlt_El1Pt",	&fTDE_Nlt_El1Pt,		"DE_Nlt_El1Pt/F");
	tree->Branch("DE_Nlt_El2Pt",	&fTDE_Nlt_El2Pt,		"DE_Nlt_El2Pt/F");
	tree->Branch("DE_Nlt_El1Eta",	&fTDE_Nlt_El1Eta,		"DE_Nlt_El1Eta/F");
	tree->Branch("DE_Nlt_El2Eta",	&fTDE_Nlt_El2Eta,		"DE_Nlt_El2Eta/F");
	tree->Branch("DE_Nll_El1Pt",	&fTDE_Nll_El1Pt,		"DE_Nll_El1Pt/F");
	tree->Branch("DE_Nll_El2Pt",	&fTDE_Nll_El2Pt,		"DE_Nll_El2Pt/F");
	tree->Branch("DE_Nll_El1Eta",	&fTDE_Nll_El1Eta,		"DE_Nll_El1Eta/F");
	tree->Branch("DE_Nll_El2Eta",	&fTDE_Nll_El2Eta,		"DE_Nll_El2Eta/F");	
}

void SSDLAnalysis::Analyze(){
	//--------------------------------------------------------------------------
	// initial event selection: good event trigger, good primary vertex...
	//--------------------------------------------------------------------------
	if( !IsGoodMuEvent() && !IsGoodElEvent() && !IsGoodElFakesEvent() && !IsGoodHadronicEvent() ) return;
	if( !IsGoodEvent() ) return;
	ResetTree();

	//--------------------------------------------------------------------------
	// Form array of indices of good muons, electrons, jets and photons
	//--------------------------------------------------------------------------
	vector<int>	selectedMuInd  = MuonSelection();
	vector<int>	selectedElInd  = ElectronSelection(); 
	vector<int>	selectedJetInd = PFJetSelection();
	vector<int>	selectedPhoInd = PhotonSelection();
	fTnqmus   = selectedMuInd.size();
	fTnqels   = selectedElInd.size();
	fTnqjets  = selectedJetInd.size();
	fTnqphos  = selectedPhoInd.size();
	
	//--------------------------------------------------------------------------
	// Dump basic run and trigger information
	//--------------------------------------------------------------------------
	DumpRunAndTiggerProperties();

	//--------------------------------------------------------------------------
	// Jets and jet-related properties, up to min(fTnqjets, gMaxNjets)
	//--------------------------------------------------------------------------
	DumpJetMETProperties(selectedJetInd);
	// get the total 4-momenta,  total transverse 3-momenta 
	TLorentzVector jtotP = jetTotalP(selectedJetInd);
	TVector3       jtotPT(jtotP.Px(), jtotP.Py(), 0.);	
	
	//--------------------------------------------------------------------------
	// Photons up to min(fTnqphos, gMaxNphos)
	//--------------------------------------------------------------------------
	DumpPhotonProperties(selectedPhoInd);
	
	//--------------------------------------------------------------------------
	// Muons up to min(fTnqmus, gMaxNmus)
	//--------------------------------------------------------------------------
	DumpMuonProperties(selectedMuInd);

	//--------------------------------------------------------------------------
	// Electrons up to min(fTnqels, gMaxNeles)
	//--------------------------------------------------------------------------
	DumpElectronProperties(selectedElInd, jtotPT);

	//--------------------------------------------------------------------------
	// Calculate different con/transverse (leptonic) alphas for the event with
	// (at least two) electron(s), muon(s), photon(s) and jet(s)
	//--------------------------------------------------------------------------
	if (fTnqmus + fTnqels + fTnqjets >= 2)
		transverseAlphas(selectedElInd, selectedMuInd, selectedPhoInd, selectedJetInd, fTalphaT_h, fTalphaCT_h, fTalphaT, fTalphaCT, fTalphaT_new, fTalphaCT_new);

	//--------------------------------------------------------------------------
	// Event flagging and dumping of pt and eta for fake ratio analysis
	//--------------------------------------------------------------------------
	DumpFPRatioProperties();
	
	fAnalysisTree->Fill();
}

void SSDLAnalysis::ResetTree(){
	ResetTriggerVariables	();
	ResetMuonVariables		();
	ResetElectronVariables	();
	ResetAuxVariables		();
	ResetFPRVariables		();
}

void SSDLAnalysis::ResetTriggerVariables(){
	// muon trigger
	fT_HLTMu9				= 0;
	fT_HLTMu11				= 0;
	fT_HLTMu15				= 0;
	fT_HLTDoubleMu3			= 0;
	// hadronic trigger
	fTHLT_Jet15U			= 0;
	fTHLT_Jet30U			= 0;
	fTHLT_Jet50U			= 0;
	fTHLT_Jet70U			= 0;
	fTHLT_Jet100U			= 0;
	fTHLT_HT100U			= 0;
	fTHLT_HT120U			= 0;
	fTHLT_HT140U			= 0;
	fTHLT_HT150U			= 0;
	// electron trigger
	fTHLT_Ele10_LW_L1R		= 0;
	fTHLT_Ele10_SW_L1R		= 0;
	fTHLT_Ele15_LW_L1R		= 0;
	fTHLT_Ele15_SW_L1R		= 0;
	fTHLT_Ele15_SW_CaloEleId_L1R= 0;
	fTHLT_Ele20_SW_L1R		= 0;
	fTHLT_DoubleEle5_SW_L1R	= 0;	
	fTHLT_DoubleEle10_SW_L1R= 0;	
	fTHLT_DoubleEle15_SW_L1R_v1	= 0;	
	// trigger summary
	fTHTL_GoodMuEvent		= 0;
	fTHTL_GoodElEvent		= 0;
	fTHTL_GoodElFakesEvent	= 0;
	fTHTL_GoodHadronicEvent	= 0;
}

void SSDLAnalysis::ResetMuonVariables(){
	fTmuMT      = -999.99;
	fTmuMinv    = -999.99;
	fTnqmus = 0;
	for(int i = 0; i < gMaxNmus; i++){
		fTmupt          [i] = -999.99;
		fTmueta         [i] = -999.99;
		fTmuphi         [i] = -999.99;
		fTmucharge      [i] = -999;
		fTmutight       [i] = -999;
		fTmuiso         [i] = -999.99;
		fTmuisohyb      [i] = -999.99;
		fTmuDRjet       [i] = -999.99;
		fTmuDRhardestjet[i] = -999.99;
		fTmud0          [i] = -999.99;
		fTmudz          [i] = -999.99;
		fTmud0bs        [i] = -999.99;
		fTmudzbs        [i] = -999.99;
		fTmuptE         [i] = -999.99;
		fTmuid          [i] = -999;
		fTmumoid        [i] = -999;
		fTmugmoid       [i] = -999;
		fTmutype        [i] = -999;
		fTmumotype      [i] = -999;
		fTmugmotype     [i] = -999;
	}	
}

void SSDLAnalysis::ResetElectronVariables(){
	// electron properties
	fTnqels = 0;
	for(int i = 0; i < gMaxNeles; i++){
		fTElcharge			[i]	= -999;
		fTElChargeIsCons	[i]	= -999;
		fTElChargeIsGenCons	[i]	= -999;
		fTElpt				[i]	= -999.99;
		fTEleta				[i]	= -999.99;
		fTElphi				[i]	= -999.99;
		fTEld0				[i]	= -999.99;
		fTElD0Err			[i]	= -999.99;
		fTElEoverP			[i]	= -999.99;
		fTElHoverE			[i]	= -999.99;
		fTElSigmaIetaIeta	[i]	= -999.99;
		fTElDeltaPhiSuperClusterAtVtx	[i]	= -999.99;
		fTElDeltaPhiSuperClusterAtVtx	[i]	= -999.99;
		fTElIDsimpleWP80relIso			[i]	= -999.99;
		fTElIDsimpleWPrelIso			[i]	= -999.99;
		fTElIDsimpleWP95relIso			[i]	= -999.99;
		fTElRelIso						[i]	= -999.99;
		fTElDR04TkSumPt					[i]	= -999.99;
		fTElDR04EcalRecHitSumEt			[i]	= -999.99;
		fTElDR04HcalTowerSumEt			[i]	= -999.99;
		fTElS4OverS1					[i]	= -999.99;
		fTElConvPartnerTrkDist			[i]	= -999.99;
		fTElConvPartnerTrkDCot			[i]	= -999.99;
		fTElChargeMisIDProb				[i]	= -999.99;
		fTElMT							[i]	= -999.99;
		fTElDRjet						[i]	= -999.99;
		fTElDRhardestjet				[i]	= -999.99;
		fTElGenID						[i]	= -999;
		fTElGenStatus					[i]	= -999;
		fTElGenMID						[i]	= -999;
		fTElGenMStatus					[i]	= -999;
		fTElGenGMID						[i]	= -999;
		fTElGenGMStatus					[i]	= -999;	
		fTElTight						[i]	= -999;	
		fTElHybRelIso					[i]	= -999.99;	
	}
		
	// di-electron properties	
	fTElminv				= -999.99;
	fTElmtinv				= -999.99;
}

void SSDLAnalysis::ResetAuxVariables(){
	// sample/run
	fTRunNumber		= 0;
	fTEventNumber	= 0;
	fTLumiSection	= 0;
	fTextxslo		= -999.99;
	fTintxs			= -999.99;
	// jet-MET properties
	fTnqjets		= 0;
	for(int i = 0; i < gMaxNjets; i++){
		fTJetpt				[i]	= -999.99;
		fTJeteta			[i]	= -999.99;
		fTJetphi			[i]	= -999.99;
	}
	fTHT			= -999.99;
	fTSumEt			= -999.99;
	fTtcMET			= -999.99;
	fTpfMET			= -999.99;
	fTMuCorrMET		= -999.99;
	fTdPhiMJ1		= -999.99;
	fTdPhiMJ2		= -999.99;
	fTR12			= -999.99;
	fTR21			= -999.99;
	fTR12plusR21	= -999.99;
	// photon properties
	fTnqphos		= 0;
	for(int i = 0; i < gMaxNphos; i++){
		fTPhopt				[i]	= -999.99;
		fTPhoeta			[i]	= -999.99;
		fTPhophi			[i]	= -999.99;
		fTPhoRelIso			[i]	= -999.99;
		fTPhoDRjet			[i]	= -999.99;
		fTPhoDRhardestjet	[i]	= -999.99;
	}
}

void SSDLAnalysis::ResetFPRVariables(){
	// other properties
	fTalphaT_h		= -999.99;
	fTalphaCT_h		= -999.99;
	fTalphaT		= -999.99;
	fTalphaCT		= -999.99;
	fTalphaT_new	= -999.99;
	fTalphaCT_new	= -999.99;
	fTElmt2_0				= -999.99;
	fTElmt2_50				= -999.99;
	fTElmt2_100				= -999.99;
	fTElmCT					= -999.99;
	fTElmCTorth				= -999.99;
	fTElmCTparl				= -999.99;
	fTElmT2orth_0			= -999.99;
	fTElmT2orth_50			= -999.99;
	fTElmT2orth_100			= -999.99;	
	fTMumt2_0				= -999.99;
	fTMumt2_50				= -999.99;
	fTMumt2_100				= -999.99;
	fTMumCT					= -999.99;
	fTMumCTorth				= -999.99;
	fTMumCTparl				= -999.99;
	fTMumT2orth_0			= -999.99;
	fTMumT2orth_50			= -999.99;
	fTMumT2orth_100			= -999.99;
	
	// fake rate propeties
	fTisSE_QCDLike			= -999;
	fTSE_QCDLike_FakeElGenID= -999;
	fTSE_QCDLike_ElLoosePt	= -999.99;
	fTSE_QCDLike_ElTightPt	= -999.99;
	fTSE_QCDLike_ElLooseEta	= -999.99;
	fTSE_QCDLike_ElTightEta	= -999.99;
	
	fTisSE_AntiQCDLike			= -999;
	fTSE_AntiQCDLike_FakeElGenID= -999;
	fTSE_AntiQCDLike_ElLoosePt	= -999.99;
	fTSE_AntiQCDLike_ElTightPt	= -999.99;
	fTSE_AntiQCDLike_ElLooseEta	= -999.99;
	fTSE_AntiQCDLike_ElTightEta	= -999.99;
	
	fTisDE_ZJetsLike			= -999;
	fTDE_ZJetsLike_ElLoosePt	= -999.99;
	fTDE_ZJetsLike_ElTightPt	= -999.99;
	fTDE_ZJetsLike_ElLooseEta	= -999.99;
	fTDE_ZJetsLike_ElTightEta	= -999.99;
	fTDE_ZJetsLike_PromptElGenLoosePt	= -999.99;
	fTDE_ZJetsLike_PromptElGenTightPt	= -999.99;
	fTDE_ZJetsLike_PromptElGenLooseEta	= -999.99;
	fTDE_ZJetsLike_PromptElGenTightEta	= -999.99;
	
	fTisDE_WJetsLike			= -999;
	fTDE_WJetsLike_FakeElGenID	= -999;
	fTDE_WJetsLike_ElLoosePt	= -999.99;
	fTDE_WJetsLike_ElTightPt	= -999.99;
	fTDE_WJetsLike_ElLooseEta	= -999.99;
	fTDE_WJetsLike_ElTightEta	= -999.99;
	fTDE_WJetsLike_FakeElGenLoosePt	= -999.99;
	fTDE_WJetsLike_FakeElGenTightPt	= -999.99;
	fTDE_WJetsLike_FakeElGenLooseEta= -999.99;
	fTDE_WJetsLike_FakeElGenTightEta= -999.99;
	fTDE_WJetsLike_ElTightMT			= -999.99;
	fTDE_WJetsLike_PromptElGenLoosePt	= -999.99;
	fTDE_WJetsLike_PromptElGenTightPt	= -999.99;
	fTDE_WJetsLike_PromptElGenLooseEta	= -999.99;
	fTDE_WJetsLike_PromptElGenTightEta	= -999.99;
	fTDE_WJetsLike_PromptElGenMT		= -999.99;
	
	fTisDE_AntiWJetsLike			= -999;
	fTDE_AntiWJetsLike_FakeElGenID	= -999;
	fTDE_AntiWJetsLike_ElLoosePt	= -999.99;
	fTDE_AntiWJetsLike_ElTightPt	= -999.99;
	fTDE_AntiWJetsLike_ElLooseEta	= -999.99;
	fTDE_AntiWJetsLike_ElTightEta	= -999.99;
	fTDE_AntiWJetsLike_FakeElGenLoosePt	= -999.99;
	fTDE_AntiWJetsLike_FakeElGenTightPt	= -999.99;
	fTDE_AntiWJetsLike_FakeElGenLooseEta= -999.99;
	fTDE_AntiWJetsLike_FakeElGenTightEta= -999.99;
	fTDE_AntiWJetsLike_ElTightMT		= -999.99;
	fTDE_AntiWJetsLike_PromptElGenLoosePt	= -999.99;
	fTDE_AntiWJetsLike_PromptElGenTightPt	= -999.99;
	fTDE_AntiWJetsLike_PromptElGenLooseEta	= -999.99;
	fTDE_AntiWJetsLike_PromptElGenTightEta	= -999.99;
	fTDE_AntiWJetsLike_PromptElGenMT		= -999.99;
	
	fTisDE_SignalLike	= -999;
	fTDE_Ntt_El1Pt = -999.99;
	fTDE_Ntt_El2Pt = -999.99;
	fTDE_Ntl_El1Pt = -999.99;
	fTDE_Ntl_El2Pt = -999.99;
	fTDE_Nlt_El1Pt = -999.99;
	fTDE_Nlt_El2Pt = -999.99;
	fTDE_Nll_El1Pt = -999.99;
	fTDE_Nll_El2Pt = -999.99;		
	fTDE_Ntt_El1Eta = -999.99;
	fTDE_Ntt_El2Eta = -999.99;
	fTDE_Ntl_El1Eta = -999.99;
	fTDE_Ntl_El2Eta = -999.99;
	fTDE_Nlt_El1Eta = -999.99;
	fTDE_Nlt_El2Eta = -999.99;
	fTDE_Nll_El1Eta = -999.99;
	fTDE_Nll_El2Eta = -999.99;
}

TLorentzVector SSDLAnalysis::jetTotalP(vector<int>& selectedJetInd){
	TLorentzVector totP(0., 0., 0., 0.);
	unsigned nqjets = selectedJetInd.size();
	for(size_t i = 0; i < nqjets; ++i) {
		int index_i = selectedJetInd[i];
		TLorentzVector tempP(fTR->PFJPx[index_i], fTR->PFJPy[index_i], fTR->PFJPz[index_i], fTR->PFJE[index_i]);
		totP += tempP;
	}
	return(totP);
}

TLorentzVector SSDLAnalysis::elTotalP(vector<int>& selectedElInd){
	TLorentzVector totP(0., 0., 0., 0.);
	unsigned nqels = selectedElInd.size();
	for(size_t i = 0; i < nqels; ++i) {
		int index_i = selectedElInd[i];
		TLorentzVector tempP(fTR->ElPx[index_i], fTR->ElPy[index_i], fTR->ElPz[index_i], fTR->ElE[index_i]);
		totP += tempP;
	}
	return(totP);
}

TLorentzVector SSDLAnalysis::muTotalP(vector<int>& selectedMuInd){
	TLorentzVector totP(0., 0., 0., 0.);
	unsigned nqmus = selectedMuInd.size();
	for(size_t i = 0; i < nqmus; ++i) {
		int index_i = selectedMuInd[i];
		TLorentzVector tempP(fTR->MuPx[index_i], fTR->MuPy[index_i], fTR->MuPz[index_i], fTR->MuE[index_i]);
		totP += tempP;
	}
	return(totP);
}

TLorentzVector SSDLAnalysis::phoTotalP(vector<int>& selectedPhoInd){
	TLorentzVector totP(0., 0., 0., 0.);
	unsigned nqphos = selectedPhoInd.size();
	for(size_t i = 0; i < nqphos; ++i) {
		int index_i = selectedPhoInd[i];
		TLorentzVector tempP(fTR->PhoPx[index_i], fTR->PhoPy[index_i], fTR->PhoPz[index_i], fTR->PhoEnergy[index_i]);
		totP += tempP;
	}
	return(totP);
}

float SSDLAnalysis::jetHT(vector<int>& selectedJetInd){
	float hT(0.);
	unsigned nqjets = selectedJetInd.size();
	for(size_t i = 0; i < nqjets; ++i) {
		int index_i = selectedJetInd[i];
		hT += fTR->PFJPt[index_i];
	}
	return(hT);
}

float SSDLAnalysis::minDRtoJet(float lepEta, float lepPhi){
	// Note, this can only be called after Jet properties have been dumped
	float mindr(100.);
	for(int j = 0; j < fTnqjets; ++j){
		float dr = Util::GetDeltaR(fTJeteta[j], lepEta, fTJetphi[j], lepPhi); 
		mindr = mindr>dr? dr:mindr;
	}
	return mindr;
}

void SSDLAnalysis::transverseMasses(TLorentzVector p1, TLorentzVector p2, TVector3 jtotPT, float &fTLepminv, float &fTLepmtinv, float &fTLepmCT, float &fTLepmCTorth, float &fTLepmCTparl, float &fTLepmt2_0, float &fTLepmt2_50, float &fTLepmt2_100, float &fTLepmT2orth_0, float &fTLepmT2orth_50, float &fTLepmT2orth_100 ){
	// Calculate invariant mass of the pair, standard MT2 for three test masses (0, 50 and 100 GeV),
	// con-transverse invariant mass for the given pair of 4-momenta (standard, plus orthogonal and paralel to the jtotPT vector),
	// and orthogonal MT2 for three test masses (0, 50 and 100 GeV)
	// Needs testing/comparison with other codes...
	
	// initialize test masses
	float testMass1 = 1.e-3;
	float testMass2 = 50.;
	float testMass3 = 100.;
	
	// invariant mass of the pair
	fTLepminv = (p1+p2).Mag();
	float mtsquare = (p1+p2).Et()*(p1+p2).Et() - (p1+p2).Pt()*(p1+p2).Pt();
	fTLepmtinv = mtsquare < 0.0 ? -TMath::Sqrt(-mtsquare) : TMath::Sqrt(mtsquare);
	
	// MT2 for three test masses
	fMT2 = new Davismt2();
	double pa[3], pb[3], pmiss[3];
	pa[0] = 0.; pb[0] = 0.; pmiss[0] = 0.;
	pa[1] = p1.Px();
	pa[2] = p1.Py();
	pb[1] = p2.Px();
	pb[2] = p2.Py();
	pmiss[1] = fTR->TCMETpx;
	pmiss[2] = fTR->TCMETpy;
	fMT2->set_momenta(pa, pb, pmiss);
	fMT2->set_mn(testMass1);
	fTLepmt2_0 = fMT2->get_mt2();
	fMT2->set_mn(testMass2);
	fTLepmt2_50 = fMT2->get_mt2();
	fMT2->set_mn(testMass3);
	fTLepmt2_100 = fMT2->get_mt2();
	delete fMT2;
	
	// con-transverse invariant masses (standard, orthogonal, paralel)
	TVector3 p1T(p1.Px(), p1.Py(), 0.); 
	TVector3 p2T(p2.Px(), p2.Py(), 0.);
	TVector3 p1Tparl = (1/p1T.Mag2())*( (p1T.Dot(jtotPT))*p1T  ); 
	TVector3 p2Tparl = (1/p2T.Mag2())*( (p2T.Dot(jtotPT))*p2T  ); 
	TVector3 p1Torth = p1T - p1Tparl;
	TVector3 p2Torth = p2T - p2Tparl;
	float e1T = sqrt(p1.M2()+p1T.Mag2());
	float e2T = sqrt(p2.M2()+p2T.Mag2());
	float e1Torth = sqrt(p1.M2()+p1Torth.Mag2());
	float e2Torth = sqrt(p2.M2()+p2Torth.Mag2());
	float e1Tparl = sqrt(p1.M2()+p1Tparl.Mag2());
	float e2Tparl = sqrt(p2.M2()+p2Tparl.Mag2());
	float ATorth	= 0.5*( p1Torth.Mag()*p2Torth.Mag() + p1Torth.Dot(p2Torth) );		
	fTLepmCT			= sqrt( p1.M2()+p2.M2() + 2*( e1T*e2T+p1T*p2T ) );		
	fTLepmCTorth		= sqrt( p1.M2()+p2.M2() + 2*( e1Torth*e2Torth+p1Torth*p2Torth ) );
	fTLepmCTparl		= sqrt( p1.M2()+p2.M2() + 2*( e1Tparl*e2Tparl+p1Tparl*p2Tparl ) );

	// orthogonal MT2 for three test masses (0, 50, 100 GeV)
	fTLepmT2orth_0		= sqrt(ATorth) + sqrt(ATorth+testMass1*testMass1);		
	fTLepmT2orth_50		= sqrt(ATorth) + sqrt(ATorth+testMass2*testMass2);		
	fTLepmT2orth_100	= sqrt(ATorth) + sqrt(ATorth+testMass3*testMass3);		
}

void SSDLAnalysis::transverseAlphas(vector<int> selectedElInd, vector<int> selectedMuInd, vector<int> selectedPhoInd, vector<int> selectedJetInd, float &alphaT_h, float &alphaCT_h, float &alphaT, float &alphaCT, float &alphaT_new, float &alphaCT_new){
	// calculate alphaT and alphaCT variables for given arrays of good leptons and good jets
	// Needs testing/comparison with other codes...

	// check the number of particles (minimum 2)
	unsigned nEls	= selectedElInd.size();
	unsigned nMus	= selectedMuInd.size();
	unsigned nPhos	= selectedPhoInd.size();
	unsigned nJets	= selectedJetInd.size();
	unsigned nTot	= nEls + nMus + nPhos + nJets;
	if (nTot < 2) return;
	
	// set the maximum number of particles for the calculation of alpha variables to be 15	
	unsigned nMax = (nTot<15)?(nTot):(15);
	
	TLorentzVector eltotP	= elTotalP(selectedElInd);
	TLorentzVector mutotP	= muTotalP(selectedMuInd);
	TLorentzVector phototP	= phoTotalP(selectedPhoInd);
	TLorentzVector jtotP	= jetTotalP(selectedJetInd);
	
	TLorentzVector totP = eltotP + mutotP + phototP + jtotP;
	// get Mt and Et for all the particles together
	//	float tot_mT = totP.Mt();
	//	float tot_Et = totP.Et();
	float tot_EtSum = 0.;
	for(unsigned k=0; k < nMax; k++) {
		TLorentzVector tempP( 0. , 0. , 0. , 0. );
		if (k<nEls) {
			int index_i = selectedElInd[k];
			tempP = TLorentzVector(fTR->ElPx[index_i], fTR->ElPy[index_i], fTR->ElPz[index_i], fTR->ElE[index_i]);
		} else if (k<nEls+nMus) {
			int index_i = selectedMuInd[k-nEls];
			tempP = TLorentzVector(fTR->MuPx[index_i], fTR->MuPy[index_i], fTR->MuPz[index_i], fTR->MuE[index_i]);
		} else if (k<nEls+nMus+nPhos) {
			int index_i = selectedPhoInd[k-nEls-nMus];
			tempP = TLorentzVector(fTR->PhoPx[index_i], fTR->PhoPy[index_i], fTR->PhoPz[index_i], fTR->PhoEnergy[index_i]);
		} else {
			int index_i = selectedJetInd[k-nEls-nMus-nPhos];
			tempP = TLorentzVector(fTR->PFJPx[index_i], fTR->PFJPy[index_i], fTR->PFJPz[index_i], fTR->PFJE[index_i]);
		}
		tot_EtSum += ((  ((int)(pow(2.,(int)nMax)) - 1)  & 1<<k) > 0) * tempP.Et();
	}
	
	//	float minDiff_groupP1_EtSum = tot_EtSum;
	float minDiff_groupP1_Et_h = 0.;
	float minDiff_groupP2_Et_h = 0.;
	float minDiff_groupP1P2_EtSum = tot_EtSum;
	float minDiff_mInvCT_h = 0.;
	float minDiff_mInvT_h = 0.;
	float minDiff_groupP1_Et = 0.;
	float minDiff_groupP2_Et = 0.;
	float minDiff_groupP1P2 = tot_EtSum;
	float minDiff_mInvCT = 0.;
	float minDiff_mInvT = 0.;
	float minDiff_mInvCT_new = 0.;
	float minDiff_mInvT_new = 0.;
	
	// find two groups of particles in the system jets+electrons+muons with minimum difference in the scalar sum of their Et
	for(unsigned n=1; n < pow(2.,(int)nMax) - 1; n++) {
		TLorentzVector groupP1( 0. , 0. , 0. , 0. );
		float groupP1_EtSum = 0.;
		float groupP2_EtSum = 0.;
		for(unsigned k=0; k < nMax; k++) {
			TLorentzVector tempP( 0. , 0. , 0. , 0. );
			if (k<nEls) {
				int index_i = selectedElInd[k];
				tempP = TLorentzVector(fTR->ElPx[index_i], fTR->ElPy[index_i], fTR->ElPz[index_i], fTR->ElE[index_i]);
			} else if (k<nEls+nMus) {
				int index_i = selectedMuInd[k-nEls];
				tempP = TLorentzVector(fTR->MuPx[index_i], fTR->MuPy[index_i], fTR->MuPz[index_i], fTR->MuE[index_i]);
			} else if (k<nEls+nMus+nPhos) {
				int index_i = selectedPhoInd[k-nEls-nMus];
				tempP = TLorentzVector(fTR->PhoPx[index_i], fTR->PhoPy[index_i], fTR->PhoPz[index_i], fTR->PhoEnergy[index_i]);
			} else {
				int index_i = selectedJetInd[k-nEls-nMus-nPhos];
				tempP = TLorentzVector(fTR->PFJPx[index_i], fTR->PFJPy[index_i], fTR->PFJPz[index_i], fTR->PFJE[index_i]);
			}
			groupP1 += ((n & 1<<k) > 0) * tempP;
			groupP1_EtSum += ((n & 1<<k) > 0) * tempP.Et();
		}
		TLorentzVector groupP2 = totP-groupP1;		
		groupP2_EtSum = tot_EtSum-groupP1_EtSum;
		float mtsquare;
		
		if (fabs(groupP1_EtSum -  groupP2_EtSum) < fabs(minDiff_groupP1P2_EtSum)) {
			TVector3 p1T (groupP1.Px(), groupP1.Py(), 0.);
			TVector3 p2T (groupP2.Px(), groupP2.Py(), 0.);			
			minDiff_groupP1_Et_h = groupP1_EtSum;
			minDiff_groupP2_Et_h = groupP2_EtSum;
			minDiff_groupP1P2_EtSum = fabs(minDiff_groupP1_Et_h - minDiff_groupP2_Et_h);
			
			mtsquare = pow(tot_EtSum,(int)2) - (p1T - p2T).Mag2();
			minDiff_mInvCT_h =	mtsquare < 0.0 ? -TMath::Sqrt(-mtsquare) : TMath::Sqrt(mtsquare);
			
			mtsquare = pow(tot_EtSum,(int)2) - (p1T + p2T).Mag2();
			minDiff_mInvT_h  =	mtsquare < 0.0 ? -TMath::Sqrt(-mtsquare) : TMath::Sqrt(mtsquare);
		}
		if (fabs(groupP1.Et() -  groupP2.Et()) < fabs(minDiff_groupP1P2)) {
			TVector3 p1T (groupP1.Px(), groupP1.Py(), 0.);
			TVector3 p2T (groupP2.Px(), groupP2.Py(), 0.);			
			float e1T = groupP1.Et();
			float e2T = groupP2.Et();
			minDiff_groupP1_Et = groupP1.Et();
			minDiff_groupP2_Et = groupP2.Et();
			minDiff_groupP1P2 = fabs(minDiff_groupP1_Et - minDiff_groupP2_Et);
			
			mtsquare = pow(e1T+e2T,(int)2) - (p1T - p2T).Mag2();
			minDiff_mInvCT  =	mtsquare < 0.0 ? -TMath::Sqrt(-mtsquare) : TMath::Sqrt(mtsquare);
			
			mtsquare = pow(e1T+e2T,(int)2) - (p1T + p2T).Mag2();
			minDiff_mInvT  =	mtsquare < 0.0 ? -TMath::Sqrt(-mtsquare) : TMath::Sqrt(mtsquare);
			
			mtsquare = groupP1.M2()+groupP2.M2() + 2*( e1T*e2T-p1T.Dot(p2T) );
			minDiff_mInvCT_new  =	mtsquare < 0.0 ? -TMath::Sqrt(-mtsquare) : TMath::Sqrt(mtsquare);
			
			mtsquare = groupP1.M2()+groupP2.M2() + 2*( e1T*e2T+p1T.Dot(p2T) );
			minDiff_mInvT_new  =	mtsquare < 0.0 ? -TMath::Sqrt(-mtsquare) : TMath::Sqrt(mtsquare);
		}
	}
	
	// calculate the alphaT and alphaCT variables for pseudo-dijet system built out of nJets + nEls + nMus particles
	alphaT_h	= ((minDiff_groupP1_Et_h<(minDiff_groupP2_Et_h))?(minDiff_groupP1_Et_h)	:(minDiff_groupP2_Et_h))/minDiff_mInvT_h ;
	alphaCT_h	= ((minDiff_groupP1_Et_h<(minDiff_groupP2_Et_h))?(minDiff_groupP1_Et_h)	:(minDiff_groupP2_Et_h))/minDiff_mInvCT_h ;
	alphaT		= ((minDiff_groupP1_Et<(minDiff_groupP2_Et))	?(minDiff_groupP1_Et)	:(minDiff_groupP2_Et))/minDiff_mInvT ;
	alphaCT		= ((minDiff_groupP1_Et<(minDiff_groupP2_Et))	?(minDiff_groupP1_Et)	:(minDiff_groupP2_Et))/minDiff_mInvCT ;
	alphaT_new	= ((minDiff_groupP1_Et<(minDiff_groupP2_Et))	?(minDiff_groupP1_Et)	:(minDiff_groupP2_Et))/minDiff_mInvT_new ;
	alphaCT_new = ((minDiff_groupP1_Et<(minDiff_groupP2_Et))	?(minDiff_groupP1_Et)	:(minDiff_groupP2_Et))/minDiff_mInvCT_new ;
	
}

void SSDLAnalysis::DumpRunAndTiggerProperties() {
	// event and run info
	fTRunNumber                  = fTR->Run;
	fTEventNumber                = fTR->Event;
	fTLumiSection                = fTR->LumiSection;
	fTextxslo                    = fTR->ExtXSecLO;
	fTintxs                      = fTR->IntXSec;
	// muon triggers
	fT_HLTMu9                    = GetHLTResult("HLT_Mu9");
	fT_HLTMu11                   = GetHLTResult("HLT_Mu11");
	fT_HLTMu15                   = GetHLTResult("HLT_Mu15");
	fT_HLTDoubleMu0              = GetHLTResult("HLT_DoubleMu0");
	fT_HLTDoubleMu3              = GetHLTResult("HLT_DoubleMu3");
	// e triggers without ElID or Iso cuts - (should be used for FP ratio measurements)	
	fTHLT_Ele10_LW_L1R           = GetHLTResult("HLT_Ele10_LW_L1R");
	fTHLT_Ele10_SW_L1R           = GetHLTResult("HLT_Ele10_SW_L1R");
	fTHLT_Ele15_LW_L1R           = GetHLTResult("HLT_Ele15_LW_L1R");
	fTHLT_Ele15_SW_L1R           = GetHLTResult("HLT_Ele15_SW_L1R");
	fTHLT_Ele15_SW_CaloEleId_L1R = GetHLTResult("HLT_Ele15_SW_CaloEleId_L1R");
	fTHLT_Ele20_SW_L1R           = GetHLTResult("HLT_Ele20_SW_L1R");
	fTHLT_DoubleEle5_SW_L1R      = GetHLTResult("HLT_DoubleEle5_SW_L1R");
	fTHLT_DoubleEle10_SW_L1R     = GetHLTResult("HLT_DoubleEle10_SW_L1R");
	fTHLT_DoubleEle15_SW_L1R_v1  = GetHLTResult("HLT_DoubleEle15_SW_L1R_v1");
	// hadronic triggers
	fTHLT_Jet30U                 = GetHLTResult("HLT_Jet30U"); 
	fTHLT_Jet50U                 = GetHLTResult("HLT_Jet50U"); 
	fTHLT_Jet70U                 = GetHLTResult("HLT_Jet70U"); 
	fTHLT_Jet100U                = GetHLTResult("HLT_Jet100U"); 
	fTHLT_HT100U                 = GetHLTResult("HLT_HT100U"); 
	fTHLT_HT120U                 = GetHLTResult("HLT_HT120U"); 
	fTHLT_HT140U                 = GetHLTResult("HLT_HT140U"); 
	fTHLT_HT150U                 = GetHLTResult("HLT_HT150U"); 
	// summary of muon/electron/hadronic trigger requirements
	fTHTL_GoodMuEvent            = IsGoodMuEvent();
	fTHTL_GoodElEvent            = IsGoodElEvent();
	fTHTL_GoodElFakesEvent       = IsGoodElFakesEvent();
	fTHTL_GoodHadronicEvent      = IsGoodHadronicEvent();
}

void SSDLAnalysis::DumpJetMETProperties(vector<int>& selectedJetInd){
	// Dump basic jet and MET properties
	int jetindex(-1);
	int nqjets = selectedJetInd.size();
	for(int ind=0; ind<std::min(nqjets,gMaxNjets); ind++){
		jetindex = selectedJetInd[ind];
		// dump properties
		fTJetpt  [ind] = fTR->PFJPt  [jetindex];
		fTJeteta [ind] = fTR->PFJEta [jetindex];
		fTJetphi [ind] = fTR->JPhi   [jetindex];
	}
	// get SumEt and METs and and Ht of all "good" jets
	fTSumEt     = fTR->SumEt;
	fTtcMET     = fTR->TCMET;
	fTpfMET     = fTR->PFMET;
	fTMuCorrMET = fTR->MuCorrMET;	
	fTHT        = jetHT(selectedJetInd);
	// get dPhi and R12, R21 variables for J1, J2 and MET
	if( nqjets >= 2 ){
		float METBadJetmin = 20.;
		float METPhi = fTR->PFMETphi;
		float MET = fTR->PFMET;
		if( MET > METBadJetmin ){
			fTdPhiMJ1    = Util::DeltaPhi(fTR->PFJPhi[selectedJetInd[0]], METPhi);
			fTdPhiMJ2    = Util::DeltaPhi(fTR->PFJPhi[selectedJetInd[1]], METPhi);
			const float pi     = TMath::Pi();
			fTR12        = sqrt(fTdPhiMJ1*fTdPhiMJ1 + (pi-fTdPhiMJ2)*(pi-fTdPhiMJ2) );
			fTR21        = sqrt(fTdPhiMJ2*fTdPhiMJ2 + (pi-fTdPhiMJ1)*(pi-fTdPhiMJ1) );
			fTR12plusR21 = fTR12 + fTR21;
		}
	}
}

void SSDLAnalysis::DumpPhotonProperties(vector<int>& selectedPhoInd){
	// Dump basic photon properties
	int phoindex(-1);
	int nqphos	= selectedPhoInd.size();
	for(int ind=0; ind<std::min(nqphos,gMaxNphos); ind++){
		phoindex = selectedPhoInd[ind];
		// dump properties
		fTPhopt		[ind] = fTR->PhoPt		[phoindex];
		fTPhoeta	[ind] = fTR->PhoEta		[phoindex];
		fTPhophi    [ind] = fTR->PhoPhi		[phoindex];
		fTPhoRelIso	[ind] = fTR->PhoIso03	[phoindex];
		fTPhoDRjet[ind]			= minDRtoJet(fTPhoeta[ind], fTPhophi[ind]);
		fTPhoDRhardestjet[ind]	= Util::GetDeltaR(fTJeteta[0], fTPhoeta[ind], fTJetphi[0], fTPhophi[ind]);		
	}		
}

void SSDLAnalysis::DumpMuonProperties(vector<int>& selectedMuInd){
	int nqmus = selectedMuInd.size();
	if( nqmus < 1 ) return;
	for(int i = 0; i < std::min(nqmus, gMaxNmus); ++i){
		int index = selectedMuInd[i];
		fTmupt       [i] = fTR->MuPt[index];
		fTmueta      [i] = fTR->MuEta[index];
		fTmuphi      [i] = fTR->MuPhi[index];
		fTmucharge   [i] = fTR->MuCharge[index];
		if(IsTightMu(index))        fTmutight[i] = 1;
		if(IsLooseNoTightMu(index)) fTmutight[i] = 0;
		fTmuiso      [i] = fTR->MuRelIso03[index];
		if(fTmupt[i] > 20.) fTmuisohyb[i] = fTmuiso[i];
		else fTmuisohyb[i] = fTmuiso[i]*fTmupt[i] / 20.;
		fTmud0       [i] = fTR->MuD0PV[index];
		fTmudz       [i] = fTR->MuDzPV[index];
		fTmud0bs     [i] = fTR->MuD0BS[index];
		fTmudzbs     [i] = fTR->MuDzBS[index];
		fTmuptE      [i] = fTR->MuPtE[index];
		
		fTmuDRjet		[i]	= minDRtoJet(fTmueta[i], fTmuphi[i]);
		fTmuDRhardestjet[i] = Util::GetDeltaR(fTJeteta[0], fTmueta[i], fTJetphi[0], fTmuphi[i]);

		fTmuid     [i] = fTR->MuGenID  [index];
		fTmumoid   [i] = fTR->MuGenMID [index];
		fTmugmoid  [i] = fTR->MuGenGMID[index];
		pdgparticle mu, mo, gmo;
		GetPDGParticle(mu,  abs(fTR->MuGenID  [index]));
		GetPDGParticle(mo,  abs(fTR->MuGenMID [index]));
		GetPDGParticle(gmo, abs(fTR->MuGenGMID[index]));
		fTmutype   [i] = mu.get_type();
		fTmumotype [i] = mo.get_type();
		fTmugmotype[i] = gmo.get_type();
	}
	
	// Calculate invariant mass, if more than 1 muon
	TLorentzVector pmu1;
	pmu1.SetXYZM(fTR->MuPx[selectedMuInd[0]], fTR->MuPy[selectedMuInd[0]], fTR->MuPz[selectedMuInd[0]], 0.105);
	if(selectedMuInd.size() > 1){
		TLorentzVector pmu2;
		pmu2.SetXYZM(fTR->MuPx[selectedMuInd[1]], fTR->MuPy[selectedMuInd[1]], fTR->MuPz[selectedMuInd[1]], 0.105);
		TLorentzVector pdimu = pmu1 + pmu2;
		fTmuMinv = pdimu.Mag();
	}

	// Calculate mT:
	double ETlept = sqrt(pmu1.M2() + pmu1.Perp2());
	double METpx  = fTR->PFMETpx;
	double METpy  = fTR->PFMETpy;
	fTmuMT        = sqrt( 2*fTR->PFMET*ETlept - pmu1.Px()*METpx - pmu1.Py()*METpy );
}

void SSDLAnalysis::DumpElectronProperties(vector<int>& selectedElInd, TVector3 jtotPT){
	// Dump basic electron properties
	int elindex(-1);
	int nqels	= selectedElInd.size();
	TLorentzVector p[nqels];
	/////// BUG ?? //////////
	TLorentzVector p_MET(fTR->PFMETpx, fTR->PFMETpy, 0, fTR->PFMET);
	for(int ind=0; ind<std::min(nqels,gMaxNeles); ind++){
		elindex = selectedElInd[ind];
		fTElpt							[ind] = fTR->ElPt					[elindex];
		fTElcharge						[ind] = fTR->ElCharge				[elindex];
		fTElChargeIsCons				[ind] = fTR->ElCInfoIsGsfCtfScPixCons[elindex];
		fTElChargeIsGenCons				[ind] = (fTR->ElCharge[elindex])==(fTR->ElGenCharge[elindex]);
		fTEleta							[ind] = fTR->ElEta					[elindex];
		fTElphi							[ind] = fTR->ElPhi					[elindex];
		fTEld0							[ind] = fTR->ElD0PV					[elindex];
		fTElD0Err						[ind] = fTR->ElD0E					[elindex];		
		fTElEoverP						[ind] = fTR->ElESuperClusterOverP	[elindex];
		fTElHoverE						[ind] = fTR->ElHcalOverEcal			[elindex];
		fTElSigmaIetaIeta				[ind] = fTR->ElSigmaIetaIeta		[elindex];
		fTElDeltaPhiSuperClusterAtVtx	[ind] = fTR->ElDeltaPhiSuperClusterAtVtx[elindex];
		fTElDeltaEtaSuperClusterAtVtx	[ind] = fTR->ElDeltaEtaSuperClusterAtVtx[elindex];
		fTElIDsimpleWP80relIso			[ind] = fTR->ElIDsimpleWP80relIso	[elindex];
		fTElIDsimpleWPrelIso			[ind] = fTR->ElIDsimpleWPrelIso		[elindex];
		fTElIDsimpleWP95relIso			[ind] = fTR->ElIDsimpleWP95relIso	[elindex];		
		fTElRelIso						[ind] = fTR->ElRelIso03				[elindex];
		fTElDR04TkSumPt					[ind] = fTR->ElDR03TkSumPt			[elindex];
		fTElDR04EcalRecHitSumEt			[ind] = fTR->ElDR03EcalRecHitSumEt	[elindex];
		fTElDR04HcalTowerSumEt			[ind] = fTR->ElDR03HcalTowerSumEt	[elindex];		
		fTElS4OverS1					[ind] = fTR->ElS4OverS1				[elindex];
		fTElConvPartnerTrkDist			[ind] = fTR->ElConvPartnerTrkDist	[elindex];
		fTElConvPartnerTrkDCot			[ind] = fTR->ElConvPartnerTrkDCot	[elindex];
		fTElChargeMisIDProb				[ind] = fTR->ElChargeMisIDProb		[elindex];		
		fTElGenID						[ind] = fTR->ElGenID				[elindex];
		fTElGenStatus					[ind] = fTR->ElGenStatus			[elindex];
		fTElGenMID						[ind] = fTR->ElGenMID				[elindex];
		fTElGenMStatus					[ind] = fTR->ElGenMStatus			[elindex];
		fTElGenGMID						[ind] = fTR->ElGenGMID				[elindex];
		fTElGenGMStatus					[ind] = fTR->ElGenGMStatus			[elindex];
		
		p[ind] = TLorentzVector(fTR->ElPx[elindex], fTR->ElPy[elindex], fTR->ElPz[elindex], fTR->ElE[elindex]);
		/////// BUG ?? //////////
		float mtsquare = (p[ind]+p_MET).Et()*(p[ind]+p_MET).Et() - (p[ind]+p_MET).Pt()*(p[ind]+p_MET).Pt();
		fTElMT			[ind]	= mtsquare < 0.0 ? -TMath::Sqrt(-mtsquare) : TMath::Sqrt(mtsquare);		
		fTElDRjet		[ind]	= minDRtoJet(fTEleta[ind], fTElphi[ind]);
		fTElDRhardestjet[ind]	= Util::GetDeltaR(fTJeteta[0], fTEleta[ind], fTJetphi[0], fTElphi[ind]);	
		fTElTight		[ind]	= IsTightEl(elindex);
		fTElHybRelIso	[ind]	= hybRelElIso(elindex);
	}		

	// calculate different con/transverse/invariant masses for the pair of two hardest electrons
	if (nqels>=2) {
		transverseMasses(p[0], p[1], jtotPT, fTElminv, fTElmtinv, fTElmCT, fTElmCTorth, fTElmCTparl, fTElmt2_0, fTElmt2_50, fTElmt2_100, fTElmT2orth_0, fTElmT2orth_50, fTElmT2orth_100 );	
	}
}

void SSDLAnalysis::DumpFPRatioProperties(){
	// Event flagging and dumping of pt and eta for fake ratio analysis
	int el1index(-1), el2index(-1);
	bool  singleElectronSelection	= SingleElectronSelection(el1index);
	bool  diElectronSelection		= DiElectronSelection(el1index, el2index);
	
	// QCD-like event (MET < 20.) with only one loose/tight electron
	if (singleElectronSelection && !diElectronSelection &&
		!((fTElmtinv>76.)&&(fTElmtinv< 106.)) &&
		(fTpfMET<20.)) {
		fTisSE_QCDLike = true;
		fTSE_QCDLike_FakeElGenID		= fTR->ElGenID[el1index];
		// store pt and eta if event is QCD-like
		DumpElectronLooseAndTighPtAndEta(el1index,	fTSE_QCDLike_ElLoosePt, fTSE_QCDLike_ElTightPt, fTSE_QCDLike_ElLooseEta, fTSE_QCDLike_ElTightEta);
	}
	// Anti-QCD-like event (MET > 30.) with only one loose/tight electron
	if (singleElectronSelection && !diElectronSelection &&
		!((fTElmtinv>76.)&&(fTElmtinv< 106.)) &&
		(fTpfMET>30.)) {
		fTisSE_AntiQCDLike = true;
		fTSE_AntiQCDLike_FakeElGenID	= fTR->ElGenID[el1index];
		// store pt and eta if event is Anti-QCD-like
		DumpElectronLooseAndTighPtAndEta(el1index,	fTSE_AntiQCDLike_ElLoosePt, fTSE_AntiQCDLike_ElTightPt, fTSE_AntiQCDLike_ElLooseEta, fTSE_AntiQCDLike_ElTightEta);
	}
	
	// ZJets-like event (76. < mT(El1,El2) < 106.) with one tight and one loose/tight electron
	if (diElectronSelection && fTnqels<3 &&
		(IsTightEl(el1index) || IsTightEl(el2index)) &&
		(fTElmtinv>76.)&&(fTElmtinv< 106.)) {
		fTisDE_ZJetsLike = true;
		
		// Loose and Tight (tag and probe)
		// decide which electron is the looser
		int elLooserIndex(-1);
		if (IsLooseEl(el1index) && IsTightEl(el2index)) elLooserIndex = el1index;
		if (IsTightEl(el1index) && IsLooseEl(el2index)) elLooserIndex = el2index;
		// store pt and eta of the "looser" electron if event is ZJets-like
		DumpElectronLooseAndTighPtAndEta(elLooserIndex,	fTDE_ZJetsLike_ElLoosePt, fTDE_ZJetsLike_ElTightPt, fTDE_ZJetsLike_ElLooseEta, fTDE_ZJetsLike_ElTightEta);
		
		// Prompt info (Tag and Probe)
		int elPromptTagIndex(-1);
		// check if electron is from Z
		bool el1IsFromZ = abs(fTR->ElGenID[el1index])==11 && abs(fTR->ElGenMID[el1index])==23; 
		bool el2IsFromZ = abs(fTR->ElGenID[el2index])==11 && abs(fTR->ElGenMID[el2index])==23;
		if (el1IsFromZ && el2IsFromZ) {
			// decide which electron is the tag one (from Z)
			elPromptTagIndex = el2index;
			// store pt and eta if event is ZJets-like and gen metched as coming from the Z
			DumpElectronLooseAndTighPtAndEta(elPromptTagIndex,	fTDE_ZJetsLike_PromptElGenLoosePt, fTDE_ZJetsLike_PromptElGenTightPt, fTDE_ZJetsLike_PromptElGenLooseEta, fTDE_ZJetsLike_PromptElGenTightEta);
		}
	}
	
	bool  ssdiElectronSelection = SSDiElectronSelection(el1index, el2index);
	// TTbar&WJets-like event (30. < MET < 80.) with one tight and one loose/tight electron
	if (ssdiElectronSelection && fTnqels<3 &&
		(IsTightEl(el1index) || IsTightEl(el2index)) &&
		!((fTElmtinv>76.)&&(fTElmtinv< 106.)) &&
		(fTpfMET>30.)&&(fTpfMET< 80.)) {
		fTisDE_WJetsLike = true;
		
		// Loose and Tight
		// decide which electron is the looser one and store the MT for the tighter one
		int elLooserIndex(-1);
		if (IsLooseEl(el1index) && IsTightEl(el2index)) {
			elLooserIndex = el1index;
			fTDE_WJetsLike_ElTightMT = fTElMT[2];
		}
		if (IsTightEl(el1index) && IsLooseEl(el2index)) {
			elLooserIndex = el2index;
			fTDE_WJetsLike_ElTightMT = fTElMT[1];
		}
		// store pt and eta of the "looser" electron if event is TTbar&WJets-like
		DumpElectronLooseAndTighPtAndEta(elLooserIndex,	fTDE_WJetsLike_ElLoosePt, fTDE_WJetsLike_ElTightPt, fTDE_WJetsLike_ElLooseEta, fTDE_WJetsLike_ElTightEta);
		
		// Fake (non-prompt) and Prompt
		int elFakeIndex(-1);
		int elPromptIndex(-1);
		// check if electron is from W
		bool el1IsFromW = abs(fTR->ElGenID[el1index])==11 && abs(fTR->ElGenMID[el1index])==24; 
		bool el2IsFromW = abs(fTR->ElGenID[el2index])==11 && abs(fTR->ElGenMID[el2index])==24;
		if (el1IsFromW ^ el2IsFromW) {
			// decide which electron is the prompt one (from W) and store the MT for that one
			if ( !el1IsFromW && el2IsFromW ) {
				elFakeIndex = el1index;
				elPromptIndex = el2index;
				fTDE_WJetsLike_PromptElGenMT = fTElMT[2];				
			}
			if ( el1IsFromW && !el2IsFromW ) {
				elFakeIndex = el2index;
				elPromptIndex = el1index;
				fTDE_WJetsLike_PromptElGenMT = fTElMT[1];
			}
			fTDE_WJetsLike_FakeElGenID		= fTR->ElGenID[elFakeIndex];
			// store pt and eta if event is TTbar&WJets-like and gen metched as coming from the W
			DumpElectronLooseAndTighPtAndEta(elFakeIndex,	fTDE_WJetsLike_FakeElGenLoosePt,	fTDE_WJetsLike_FakeElGenTightPt,	fTDE_WJetsLike_FakeElGenLooseEta,	fTDE_WJetsLike_FakeElGenTightEta);
			DumpElectronLooseAndTighPtAndEta(elPromptIndex,	fTDE_WJetsLike_PromptElGenLoosePt,	fTDE_WJetsLike_PromptElGenTightPt,	fTDE_WJetsLike_PromptElGenLooseEta, fTDE_WJetsLike_PromptElGenTightEta);
		}		
	}
	// Anti-TTbar&WJets-like event (MET > 80.) with one tight and one loose/tight electron
	if (ssdiElectronSelection &&
		(IsTightEl(el1index) || IsTightEl(el2index)) &&
		!((fTElmtinv>76.)&&(fTElmtinv< 106.)) &&
		(fTpfMET> 80.)) {
		fTisDE_AntiWJetsLike = true;
		
		// Loose and Tight
		// decide which electron is the looser one and store the MT for the tighter one
		int elLooserIndex(-1);
		if (IsLooseEl(el1index) && IsTightEl(el2index)) {
			elLooserIndex = el1index;
			fTDE_AntiWJetsLike_ElTightMT = fTElMT[2];
		}
		if (IsTightEl(el1index) && IsLooseEl(el2index)) {
			elLooserIndex = el2index;
			fTDE_AntiWJetsLike_ElTightMT = fTElMT[1];
		}
		// store pt and eta if event is TTbar&WJets-like
		DumpElectronLooseAndTighPtAndEta(elLooserIndex,	fTDE_AntiWJetsLike_ElLoosePt, fTDE_AntiWJetsLike_ElTightPt, fTDE_AntiWJetsLike_ElLooseEta, fTDE_AntiWJetsLike_ElTightEta);
		
		// Fake (non-prompt) and Prompt		
		int elFakeIndex(-1);
		int elPromptIndex(-1);				
		// check if electron is from W
		bool el1IsFromW = abs(fTR->ElGenID[el1index])==11 && abs(fTR->ElGenMID[el1index])==24; 
		bool el2IsFromW = abs(fTR->ElGenID[el2index])==11 && abs(fTR->ElGenMID[el2index])==24; 
		
		if (el1IsFromW ^ el2IsFromW) {
			// decide which electron is the prompt one (from W) and store the MT for that one
			if ( !el1IsFromW && el2IsFromW ) {
				elFakeIndex = el1index;
				elPromptIndex = el2index;
				fTDE_AntiWJetsLike_PromptElGenMT = fTElMT[2];				
			}
			if ( el1IsFromW && !el2IsFromW ) {
				elFakeIndex = el2index;
				elPromptIndex = el1index;
				fTDE_AntiWJetsLike_PromptElGenMT = fTElMT[1];				
			}
			fTDE_AntiWJetsLike_FakeElGenID		= fTR->ElGenID[elFakeIndex];
			// store pt and eta if event is TTbar&WJets-like and gen metched as coming from the W
			DumpElectronLooseAndTighPtAndEta(elFakeIndex,	fTDE_AntiWJetsLike_FakeElGenLoosePt,	fTDE_AntiWJetsLike_FakeElGenTightPt,	fTDE_AntiWJetsLike_FakeElGenLooseEta,	fTDE_AntiWJetsLike_FakeElGenTightEta);
			DumpElectronLooseAndTighPtAndEta(elPromptIndex,	fTDE_AntiWJetsLike_PromptElGenLoosePt,	fTDE_AntiWJetsLike_PromptElGenTightPt,	fTDE_AntiWJetsLike_PromptElGenLooseEta, fTDE_AntiWJetsLike_PromptElGenTightEta);
		}
	}

	// Signal-like event (MET>80.) with two loose/tight electrons
	if (ssdiElectronSelection &&
		(fTpfMET> 80.)) {
		fTisDE_SignalLike = true;
		if(fVerbose > 0) cout <<	"Signal-like event [" << "RunNumber: " << fTRunNumber << ", EventNumber: " << fTEventNumber << ", LumiSection: " << fTLumiSection << "]" << std::endl <<
		"                  [nqmus: " << fTnqmus << ", nqels: " << fTnqels << ", nqjets: " << fTnqjets << "]" << endl;
		if (IsTightEl(el1index) && IsTightEl(el2index))				
			DumpTwoElectronPtAndEta(el1index, el2index, fTDE_Ntt_El1Pt, fTDE_Ntt_El2Pt, fTDE_Ntt_El1Eta, fTDE_Ntt_El2Eta);
		if (IsTightEl(el1index) && IsLooseNoTightEl(el2index))		
			DumpTwoElectronPtAndEta(el1index, el2index, fTDE_Ntl_El1Pt, fTDE_Ntl_El2Pt, fTDE_Ntl_El1Eta, fTDE_Ntl_El2Eta);
		if (IsLooseNoTightEl(el1index) && IsTightEl(el2index))		
			DumpTwoElectronPtAndEta(el1index, el2index, fTDE_Nlt_El1Pt, fTDE_Nlt_El2Pt, fTDE_Nlt_El1Eta, fTDE_Nlt_El2Eta);
		if (IsLooseNoTightEl(el1index) && IsLooseNoTightEl(el2index))
			DumpTwoElectronPtAndEta(el1index, el2index, fTDE_Nll_El1Pt, fTDE_Nll_El2Pt, fTDE_Nll_El1Eta, fTDE_Nll_El2Eta);
	}	
}	

void SSDLAnalysis::DumpElectronLooseAndTighPtAndEta(int elindex, float &elLoosePt, float &elTightPt, float &elLooseEta, float &elTightEta) {
	if (IsLooseEl(elindex))	elLoosePt	= fTR->ElPt [elindex];
	if (IsTightEl(elindex))	elTightPt	= fTR->ElPt [elindex];
	if (IsLooseEl(elindex))	elLooseEta	= fTR->ElEta[elindex];
	if (IsTightEl(elindex))	elTightEta	= fTR->ElEta[elindex];
}

void SSDLAnalysis::DumpTwoElectronPtAndEta(int el1index, int el2index, float &el1Pt, float &el2Pt, float &el1Eta, float &el2Eta) {
	el1Pt		= fTR->ElPt [el1index];
	el2Pt		= fTR->ElPt [el2index];
	el1Eta		= fTR->ElEta [el1index];
	el2Eta		= fTR->ElEta [el2index];
}
