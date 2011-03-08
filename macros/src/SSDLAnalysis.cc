#include "SSDLAnalysis.hh"

using namespace std;

const int SSDLAnalysis::fMaxNjets;
const int SSDLAnalysis::fMaxNmus;
const int SSDLAnalysis::fMaxNeles;
const int SSDLAnalysis::gMaxhltbits;

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


    // run/sample properties
	fAnalysisTree->Branch("Run",                &fTRunNumber,         "Run/I");
	fAnalysisTree->Branch("Event",              &fTEventNumber,       "Event/I");
	fAnalysisTree->Branch("LumiSec",            &fTLumiSection,       "LumiSec/I");
	// HLT triggers
	fAnalysisTree->Branch("NHLTPaths",          &fTHLTNPaths         ,"NHLTPaths/I");
	fAnalysisTree->Branch("HLTResults",         &fTHLTres            ,"HLTResults[NHLTPaths]/I");
	fAnalysisTree->Branch("HLTPrescales",       &fTHLTprescale       ,"HLTPrescale[NHLTPaths]/I");
	fAnalysisTree->Branch("HLTNames",           &fTHLTnames);

	// single-muon properties
	fAnalysisTree->Branch("NMus"          ,&fTnqmus,          "NMus/I");
	fAnalysisTree->Branch("MuPt"          ,&fTmupt,           "MuPt[NMus]/F");
	fAnalysisTree->Branch("MuEta"         ,&fTmueta,          "MuEta[NMus]/F");
	fAnalysisTree->Branch("MuPhi"         ,&fTmuphi,          "MuPhi[NMus]/F");
	fAnalysisTree->Branch("MuCharge"      ,&fTmucharge,       "MuCharge[NMus]/I");
	fAnalysisTree->Branch("MuTight"       ,&fTmutight,        "MuTight[NMus]/I");
	fAnalysisTree->Branch("MuIso"         ,&fTmuiso,          "MuIso[NMus]/F");
	fAnalysisTree->Branch("MuIsoHybrid"   ,&fTmuisohyb,       "MuIsoHybrid[NMus]/F");
	fAnalysisTree->Branch("MuD0"          ,&fTmud0,           "MuD0[NMus]/F");
	fAnalysisTree->Branch("MuDz"          ,&fTmudz,           "MuDz[NMus]/F");
	fAnalysisTree->Branch("MuD0BS"        ,&fTmud0bs,         "MuD0BS[NMus]/F");
	fAnalysisTree->Branch("MuDzBS"        ,&fTmudzbs,         "MuDzBS[NMus]/F");
	fAnalysisTree->Branch("MuPtE"         ,&fTmuptE,          "MuPtE[NMus]/F");
	fAnalysisTree->Branch("MuGenID"       ,&fTmuid,           "MuGenID[NMus]/I");
	fAnalysisTree->Branch("MuGenMoID"     ,&fTmumoid,         "MuGenMoID[NMus]/I");
	fAnalysisTree->Branch("MuGenGMoID"    ,&fTmugmoid,        "MuGenGMoID[NMus]/I");
	fAnalysisTree->Branch("MuGenType"     ,&fTmutype,         "MuGenType[NMus]/I");
	fAnalysisTree->Branch("MuGenMoType"   ,&fTmumotype,       "MuGenMoType[NMus]/I");
	fAnalysisTree->Branch("MuGenGMoType"  ,&fTmugmotype,      "MuGenGMoType[NMus]/I");
	fAnalysisTree->Branch("MuMT"          ,&fTmuMT,           "MuMT[NMus]/F");

	// single-electron properties
	fAnalysisTree->Branch("NEls",                   &fTnqels,               "NEls/I");
	fAnalysisTree->Branch("ElCh",                   &fTElcharge,            "ElCh[NEls]/I");
	fAnalysisTree->Branch("ElChIsCons",             &fTElChargeIsCons,      "ElChIsCons[NEls]/I");
	fAnalysisTree->Branch("ElChIsGenCons",          &fTElChargeIsGenCons,   "ElChIsGenCons[NEls]/I");
	fAnalysisTree->Branch("ElEcalDriven",           &fTElEcalDriven,        "ElEcalDriven[NEls]/I");
	fAnalysisTree->Branch("ElCaloEnergy",           &fTElCaloEnergy,        "ElCaloEnergy[NEls]/F");
	fAnalysisTree->Branch("ElPt",                   &fTElpt,                "ElPt[NEls]/F");
	fAnalysisTree->Branch("ElEta",                  &fTEleta,               "ElEta[NEls]/F");
	fAnalysisTree->Branch("ElPhi",                  &fTElphi,               "ElPhi[NEls]/F");
	fAnalysisTree->Branch("ElD0",                   &fTEld0,                "ElD0[NEls]/F");
	fAnalysisTree->Branch("ElD0Err",                &fTElD0Err,             "ElD0Err[NEls]/F");
	fAnalysisTree->Branch("ElDz",                   &fTEldz,                "ElDz[NEls]/F");
	fAnalysisTree->Branch("ElDzErr",                &fTElDzErr,             "ElDzErr[NEls]/F");
	fAnalysisTree->Branch("ElEoverP",               &fTElEoverP,            "ElEoverP[NEls]/F");
	fAnalysisTree->Branch("ElHoverE",               &fTElHoverE,            "ElHoverE[NEls]/F");
	fAnalysisTree->Branch("ElSigmaIetaIeta",        &fTElSigmaIetaIeta,     "ElSigmaIetaIeta[NEls]/F");
	fAnalysisTree->Branch("ElDeltaPhiSuperClusterAtVtx", &fTElDeltaPhiSuperClusterAtVtx, "ElDeltaPhiSuperClusterAtVtx[NEls]/F");
	fAnalysisTree->Branch("ElDeltaEtaSuperClusterAtVtx", &fTElDeltaEtaSuperClusterAtVtx, "ElDeltaEtaSuperClusterAtVtx[NEls]/F");
	fAnalysisTree->Branch("ElRelIso",               &fTElRelIso,            "ElRelIso[NEls]/F");
	fAnalysisTree->Branch("ElIsGoodElId_WP80",      &fTElIsGoodElId_WP80,   "ElIsGoodElId_WP80[NEls]/I");
	fAnalysisTree->Branch("ElIsGoodElId_WP90",      &fTElIsGoodElId_WP90,   "ElIsGoodElId_WP90[NEls]/I");
	fAnalysisTree->Branch("ElIsGoodElId_WP95",      &fTElIsGoodElId_WP95,   "ElIsGoodElId_WP95[NEls]/I");
	fAnalysisTree->Branch("ElS4OverS1",             &fTElS4OverS1,          "ElS4OverS1[NEls]/F");
	fAnalysisTree->Branch("ElGenID",                &fTElGenID,             "ElGenID[NEls]/I");
	fAnalysisTree->Branch("ElGenMID",               &fTElGenMID,            "ElGenMID[NEls]/I");
	fAnalysisTree->Branch("ElGenGMID",              &fTElGenGMID,           "ElGenGMID[NEls]/I");
	fAnalysisTree->Branch("ElGenType",              &fTElGenType,           "ElGenType[NEls]/I");
	fAnalysisTree->Branch("ElGenMType",             &fTElGenMType,          "ElGenMType[NEls]/I");
	fAnalysisTree->Branch("ElGenGMType",            &fTElGenGMType,         "ElGenGMType[NEls]/I");
	fAnalysisTree->Branch("ElHybRelIso",            &fTElHybRelIso,         "ElHybRelIso[NEls]/F");
	fAnalysisTree->Branch("ElMT",                   &fTElMT,                "ElMT[NEls]/F");

	// event properties
	fAnalysisTree->Branch("tcMET",           &fTtcMET,       "tcMET/F");
	fAnalysisTree->Branch("pfMET",           &fTpfMET,       "pfMET/F");
	// jet-MET properties
	fAnalysisTree->Branch("NJets",           &fTnqjets,      "NJets/I");
	fAnalysisTree->Branch("JetPt",           &fTJetpt,       "JetPt[NJets]/F");
	fAnalysisTree->Branch("JetEta",          &fTJeteta,      "JetEta[NJets]/F");
	fAnalysisTree->Branch("JetPhi",          &fTJetphi,      "JetPhi[NJets]/F");
}

void SSDLAnalysis::Analyze(){
	// initial event selection: good event trigger, good primary vertex...
	if( !IsGoodMuEvent() && !IsGoodElEvent() && !IsGoodElFakesEvent() && !IsGoodHadronicEvent()) return;
	if( !IsGoodEvent() ) return;
	ResetTree();
	
	// Do object selections
	vector<int> selectedMuInd  = MuonSelection(&UserAnalysisBase::IsGoodBasicMu);
	vector<int> selectedElInd  = ElectronSelection(&UserAnalysisBase::IsLooseEl);
	vector<int> selectedJetInd = PFJetSelection(30., 2.5, &UserAnalysisBase::IsGoodBasicPFJet);
	fTnqmus  = std::min((int) selectedMuInd .size(), fMaxNmus);
	fTnqels  = std::min((int) selectedElInd .size(), fMaxNeles);
	fTnqjets = std::min((int) selectedJetInd.size(), fMaxNjets);
	
	// Require at least one loose lepton
	if( (fTnqmus + fTnqels) < 1 ) return;
	
	// event and run info
	fTRunNumber                  = fTR->Run;
	fTEventNumber                = fTR->Event;
	fTLumiSection                = fTR->LumiSection;
		
	// Get trigger results
	fTHLTnames = fHLTLabels;
	fTHLTNPaths = fHLTLabels.size();
	for(unsigned int i = 0; i < fHLTLabels.size(); i++ ){
		fTHLTres[i]      = fTR->HLTResults[i];
		fTHLTprescale[i] = fTR->HLTPrescale[i];
	}
		
	// Dump basic jet and MET properties
	int jetindex(-1);
	int nqjets = selectedJetInd.size();
	for(int ind=0; ind<std::min(nqjets,fMaxNjets); ind++){
		jetindex = selectedJetInd[ind];
		// dump properties
		fTJetpt  [ind] = fTR->PFJPt [jetindex];
		fTJeteta [ind] = fTR->PFJEta[jetindex];
		fTJetphi [ind] = fTR->PFJPhi[jetindex];
	}
	// get METs
	fTtcMET     = fTR->TCMET;
	fTpfMET     = fTR->PFMET;
	// get dPhi and R12, R21 variables for J1, J2 and MET
	if( nqjets >= 2 ){
		float METBadJetmin = 20.;
		float METPhi = fTR->PFMETphi;
		float MET = fTR->PFMET;
	}
	
	int nqmus = selectedMuInd.size();
	if( nqmus < 1 ) return;
	for(int i = 0; i < std::min(nqmus, fMaxNmus); ++i){
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

		// Calculate mT:
		TLorentzVector pmu;
		pmu.SetXYZM(fTR->MuPx[index], fTR->MuPy[index], fTR->MuPz[index], 0.105);	
		double ETlept = sqrt(pmu.M2() + pmu.Perp2());
		double METpx  = fTR->PFMETpx;
		double METpy  = fTR->PFMETpy;
		fTmuMT[i]     = sqrt( 2*fTR->PFMET*ETlept - pmu.Px()*METpx - pmu.Py()*METpy );
	}

	// Dump basic electron properties
	int elindex(-1);
	int nqels = selectedElInd.size();
	TLorentzVector p[nqels];
	TLorentzVector p_MET(fTR->PFMETpx, fTR->PFMETpy, 0, fTR->PFMET);
	for(int ind=0; ind<std::min(nqels,fMaxNeles); ind++){
		elindex = selectedElInd[ind];
		fTElcharge                   [ind] = fTR->ElCharge                [elindex];
		fTElChargeIsCons             [ind] = fTR->ElCInfoIsGsfCtfScPixCons[elindex];
		fTElChargeIsGenCons          [ind] = (fTR->ElCharge[elindex])==(fTR->ElGenCharge[elindex]);
		fTElEcalDriven               [ind] = fTR->ElEcalDriven           [elindex];
		fTElCaloEnergy               [ind] = fTR->ElCaloEnergy           [elindex];
		fTElpt                       [ind] = fTR->ElPt                   [elindex];
		fTEleta                      [ind] = fTR->ElEta                  [elindex];
		fTElphi                      [ind] = fTR->ElPhi                  [elindex];
		fTEld0                       [ind] = fTR->ElD0PV                 [elindex];
		fTElD0Err                    [ind] = fTR->ElD0E                  [elindex];
		fTEldz                       [ind] = fTR->ElDzPV                 [elindex];
		fTElDzErr                    [ind] = fTR->ElDzE                  [elindex];
		fTElEoverP                   [ind] = fTR->ElESuperClusterOverP   [elindex];
		fTElHoverE                   [ind] = fTR->ElHcalOverEcal         [elindex];
		fTElSigmaIetaIeta            [ind] = fTR->ElSigmaIetaIeta        [elindex];
		fTElDeltaPhiSuperClusterAtVtx[ind] = fTR->ElDeltaPhiSuperClusterAtVtx[elindex];
		fTElDeltaEtaSuperClusterAtVtx[ind] = fTR->ElDeltaEtaSuperClusterAtVtx[elindex];
		fTElRelIso                   [ind] = fTR->ElRelIso03                 [elindex];
		fTElS4OverS1                 [ind] = fTR->ElS4OverS1                 [elindex];
		fTElGenID                    [ind] = fTR->ElGenID                    [elindex];
		fTElGenMID                   [ind] = fTR->ElGenMID                   [elindex];
		fTElGenGMID                  [ind] = fTR->ElGenGMID                  [elindex];
		
		pdgparticle el, emo, egmo;
		GetPDGParticle(el,   abs(fTR->ElGenID  [elindex]));
		GetPDGParticle(emo,  abs(fTR->ElGenMID [elindex]));
		GetPDGParticle(egmo, abs(fTR->ElGenGMID[elindex]));
		fTElGenType  [ind] = el.get_type();
		fTElGenMType [ind] = emo.get_type();
		fTElGenGMType[ind] = egmo.get_type();
		
		p[ind] = TLorentzVector(fTR->ElPx[elindex], fTR->ElPy[elindex], fTR->ElPz[elindex], fTR->ElE[elindex]);
		float mtsquare = (p[ind]+p_MET).Et()*(p[ind]+p_MET).Et() - (p[ind]+p_MET).Pt()*(p[ind]+p_MET).Pt();
		fTElMT             [ind] = mtsquare < 0.0 ? -TMath::Sqrt(-mtsquare) : TMath::Sqrt(mtsquare);
		fTElHybRelIso      [ind] = hybRelElIso(elindex);
		fTElIsGoodElId_WP80[ind] = IsGoodElId_WP80(elindex);
		fTElIsGoodElId_WP90[ind] = IsGoodElId_WP90(elindex);
		fTElIsGoodElId_WP95[ind] = IsGoodBasicEl(elindex);
	}

	fAnalysisTree->Fill();
}

void SSDLAnalysis::ResetTree(){
	// sample/run
	fTRunNumber                  = 0;
	fTEventNumber                = 0;
	fTLumiSection                = 0;

	fTHLTNPaths	= 0;
	for(size_t i = 0; i < gMaxhltbits; ++i){
		fTHLTres     [i] = -1;
		fTHLTprescale[i] = -1;
	}
	fTHLTnames.clear();
	fTHLTnames.resize(fHLTLabelMap.size());

	// muon properties
	fTnqmus = 0;
	for(int i = 0; i < fMaxNmus; i++){
		fTmupt          [i] = -999.99;
		fTmueta         [i] = -999.99;
		fTmuphi         [i] = -999.99;
		fTmucharge      [i] = -999;
		fTmutight       [i] = -999;
		fTmuiso         [i] = -999.99;
		fTmuisohyb      [i] = -999.99;
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
		fTmuMT          [i] = -999.99;
	}

	// electron properties
	fTnqels = 0;
	for(int i = 0; i < fMaxNeles; i++){
		fTElcharge          [i] = -999;
		fTElChargeIsCons    [i] = -999;
		fTElEcalDriven      [i] = -999;
		fTElChargeIsGenCons [i] = -999;
		fTElCaloEnergy      [i] = -999.99;
		fTElpt              [i] = -999.99;
		fTEleta             [i] = -999.99;
		fTElphi             [i] = -999.99;
		fTEld0              [i] = -999.99;
		fTElD0Err           [i] = -999.99;
		fTEldz              [i] = -999.99;
		fTElDzErr           [i] = -999.99;
		fTElEoverP          [i] = -999.99;
		fTElHoverE          [i] = -999.99;
		fTElSigmaIetaIeta   [i] = -999.99;
		fTElDeltaPhiSuperClusterAtVtx[i] = -999.99;
		fTElDeltaPhiSuperClusterAtVtx[i] = -999.99;
		fTElRelIso                   [i] = -999.99;
		fTElIsGoodElId_WP80          [i] = -999;
		fTElIsGoodElId_WP90          [i] = -999;
		fTElIsGoodElId_WP95          [i] = -999;
		fTElS4OverS1                 [i] = -999.99;
		fTElMT                       [i] = -999.99;
		fTElGenID                    [i] = -999;
		fTElGenMID                   [i] = -999;
		fTElGenGMID                  [i] = -999;
		fTElGenType                  [i] = -999;
		fTElGenMType                 [i] = -999;
		fTElGenGMType                [i] = -999;
		fTElHybRelIso                [i] = -999.99;
	}

	// jet-MET properties
	fTnqjets = 0;
	for(int i = 0; i < fMaxNjets; i++){
		fTJetpt [i] = -999.99;
		fTJeteta[i] = -999.99;
		fTJetphi[i] = -999.99;
	}
	fTtcMET      = -999.99;
	fTpfMET      = -999.99;
}
