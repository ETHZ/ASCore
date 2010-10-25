#include "MuonAnalysis.hh"
#include "base/TreeReader.hh"
#include "helper/Utilities.hh"

#include <TH1D.h>
#include <TH2D.h>
#include <TVector3.h>
#include <TLorentzVector.h>

using namespace std;

MuonAnalysis::MuonAnalysis(TreeReader *tr) : UserAnalysisBase(tr){
}

MuonAnalysis::~MuonAnalysis(){
}

void MuonAnalysis::Begin(){
	ReadPDGTable();
	BookTree();
}

void MuonAnalysis::Analyze(){
	if(!IsGoodMuEvent()) return; // Trigger, etc.
	FillMuonTree();
}

void MuonAnalysis::End(){
	fOutputFile->cd();
	fMuTree->Write();
	fOutputFile->Close();
}

void MuonAnalysis::FillMuonTree(){
	vector<int> goodMuons = MuonSelection();
	if( goodMuons.size() < 1 ) return;
	ResetTree();
	fTrun     = fTR->Run;
	fTevent   = fTR->Event;
	fTlumisec = fTR->LumiSection;
	fThbhenoiseflag = fTR->HBHENoiseFlag;

	// Trigger info
	fT_HLTMu9       = GetHLTResult("HLT_Mu9")       ? 1:0;
	fT_HLTMu11      = GetHLTResult("HLT_Mu11")      ? 1:0;
	fT_HLTMu15      = GetHLTResult("HLT_Mu15")   ? 1:0;
	fT_HLTDoubleMu0 = GetHLTResult("HLT_DoubleMu0") ? 1:0;
	fT_HLTDoubleMu3 = GetHLTResult("HLT_DoubleMu3") ? 1:0;
	fT_HLTJet30U    = GetHLTResult("HLT_Jet30U")    ? 1:0;
	fT_HLTJet50U    = GetHLTResult("HLT_Jet50U")    ? 1:0;
	fT_HLTJet70U    = GetHLTResult("HLT_Jet70U")    ? 1:0;
	fT_HLTJet100U   = GetHLTResult("HLT_Jet100U")   ? 1:0;

	// Jet Variables
	vector<int> goodJets = PFJetSelection();
	fTnjets = goodJets.size();
	double ht = 0;
	TVector3 mht(0.,0.,0.);
	for(size_t i = 0; i < fTnjets; ++i){
		int index = goodJets[i];
		fTjpt [i] = fTR->PFJPt[index];
		fTjeta[i] = fTR->PFJEta[index];
		fTjphi[i] = fTR->PFJPhi[index];
		ht += fTR->PFJPt[index];
		mht(0) += fTR->PFJPx[index];
		mht(1) += fTR->PFJPy[index];
		mht(2) += fTR->PFJPz[index];
	}
	fTHT    = ht;
	fTMHT   = mht.Mag();
	fTPFMET   = fTR->PFMET;
	fTTCMET   = fTR->TCMET;
	fTSumET = fTR->SumEt;

	// Muon Variables
	for(size_t i = 0; i < goodMuons.size(); ++i){
		int index = goodMuons[i];
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

		double mindr(100.);
		double mindphi(100.);
		for(size_t j = 0; j < fTnjets; ++j){
			double dr   = Util::GetDeltaR(fTjeta[j], fTmueta[i], fTjphi[j], fTmuphi[i]); 
			double dphi = Util::DeltaPhi( fTjphi[j], fTmuphi[i]); 
			mindr = mindr>dr? dr:mindr;
			mindphi = mindphi>dphi? dphi:mindphi;
		}
		fTDRjet[i] = mindr;
		fTDPhijet[i] = mindphi;

		fTDRhardestjet[i] = Util::GetDeltaR(fTjeta[0], fTmueta[i], fTjphi[0], fTmuphi[i]);

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
	fTnmus = goodMuons.size();
	
	// Calculate invariant mass, if more than 1 muon
	TLorentzVector pmu1;
	pmu1.SetXYZM(fTR->MuPx[goodMuons[0]], fTR->MuPy[goodMuons[0]], fTR->MuPz[goodMuons[0]], 0.105);
	if(goodMuons.size() > 1){	
		TLorentzVector pmu2;
		pmu2.SetXYZM(fTR->MuPx[goodMuons[1]], fTR->MuPy[goodMuons[1]], fTR->MuPz[goodMuons[1]], 0.105);
		TLorentzVector pdimu = pmu1 + pmu2;
		fTMinv = pdimu.Mag();
	}

	// Calculate mT:
	double ETlept = sqrt(pmu1.M2() + pmu1.Perp2());
	double METpx  = fTR->PFMETpx;
	double METpy  = fTR->PFMETpy;
	fTMT          = sqrt( 2*fTPFMET*ETlept - pmu1.Px()*METpx - pmu1.Py()*METpy );
	
	fMuTree->Fill();
}

void MuonAnalysis::BookTree(){
	fOutputFile->cd();
	fMuTree = new TTree("MuonAnalysis", "MuonAnalysisTree");
	fMuTree->Branch("Run"           ,&fTrun,          "Run/I");
	fMuTree->Branch("Event"         ,&fTevent,        "Event/I");
	fMuTree->Branch("LumiSection"   ,&fTlumisec,      "LumiSection/I");
	fMuTree->Branch("HBHENoiseFlag" ,&fThbhenoiseflag,"HBHENoiseFlag/I");
	fMuTree->Branch("HLTMu9"        ,&fT_HLTMu9,      "HLTMu9/I");
	fMuTree->Branch("HLTMu11"       ,&fT_HLTMu11,     "HLTMu11/I");
	fMuTree->Branch("HLTMu15"       ,&fT_HLTMu15,     "HLTMu15/I");
	fMuTree->Branch("HLTDoubleMu0"  ,&fT_HLTDoubleMu0,"HLTDoubleMu0/I");
	fMuTree->Branch("HLTDoubleMu3"  ,&fT_HLTDoubleMu3,"HLTDoubleMu3/I");
	fMuTree->Branch("HLTJet30U"     ,&fT_HLTJet30U,   "HLTJet30U/I");
	fMuTree->Branch("HLTJet50U"     ,&fT_HLTJet50U,   "HLTJet50U/I");
	fMuTree->Branch("HLTJet70U"     ,&fT_HLTJet70U,   "HLTJet70U/I");
	fMuTree->Branch("HLTJet100U"    ,&fT_HLTJet100U,  "HLTJet100U/I");

	fMuTree->Branch("HT"            ,&fTHT,           "HT/F");
	fMuTree->Branch("MHT"           ,&fTMHT,          "MHT/F");
	fMuTree->Branch("PFMET"         ,&fTPFMET,        "PFMET/F");
	fMuTree->Branch("TCMET"         ,&fTTCMET,        "TCMET/F");
	fMuTree->Branch("SumET"         ,&fTSumET,        "SumET/F");
	fMuTree->Branch("MT"            ,&fTMT,           "MT/F");
	fMuTree->Branch("Minv"          ,&fTMinv,         "Minv/F");

	fMuTree->Branch("NJets"         ,&fTnjets,        "NJets/I");
	fMuTree->Branch("JPt"           ,&fTjpt,          "JPt[NJets]/F");
	fMuTree->Branch("JEta"          ,&fTjeta,         "JEta[NJets]/F");
	fMuTree->Branch("JPhi"          ,&fTjphi,         "JPhi[NJets]/F");

	fMuTree->Branch("NMus"          ,&fTnmus,         "NMus/I");
	fMuTree->Branch("MuPt"          ,&fTmupt,         "MuPt[NMus]/F");
	fMuTree->Branch("MuEta"         ,&fTmueta,        "MuEta[NMus]/F");
	fMuTree->Branch("MuPhi"         ,&fTmuphi,        "MuPhi[NMus]/F");
	fMuTree->Branch("MuCharge"      ,&fTmucharge,     "MuCharge[NMus]/I");
	fMuTree->Branch("MuTight"       ,&fTmutight,      "MuTight[NMus]/I");
	fMuTree->Branch("MuIso"         ,&fTmuiso,        "MuIso[NMus]/F");
	fMuTree->Branch("MuIsoHybrid"   ,&fTmuisohyb,     "MuIsoHybrid[NMus]/F");
	fMuTree->Branch("MuDRJet"       ,&fTDRjet,        "MuDRJet[NMus]/F");
	fMuTree->Branch("MuDPhiJet"     ,&fTDPhijet,      "MuDPhiJet[NMus]/F");
	fMuTree->Branch("MuDRHardestJet",&fTDRhardestjet, "MuDRHardestJet[NMus]/F");
	fMuTree->Branch("MuD0"          ,&fTmud0,         "MuD0[NMus]/F");
	fMuTree->Branch("MuDz"          ,&fTmudz,         "MuDz[NMus]/F");
	fMuTree->Branch("MuD0BS"        ,&fTmud0bs,       "MuD0BS[NMus]/F");
	fMuTree->Branch("MuDzBS"        ,&fTmudzbs,       "MuDzBS[NMus]/F");
	fMuTree->Branch("MuPtE"         ,&fTmuptE,        "MuPtE[NMus]/F");
	fMuTree->Branch("MuGenID"       ,&fTmuid,         "MuGenID[NMus]/I");
	fMuTree->Branch("MuGenMoID"     ,&fTmumoid,       "MuGenMoID[NMus]/I");
	fMuTree->Branch("MuGenGMoID"    ,&fTmugmoid,      "MuGenGMoID[NMus]/I");
	fMuTree->Branch("MuGenType"     ,&fTmutype,       "MuGenType[NMus]/I");
	fMuTree->Branch("MuGenMoType"   ,&fTmumotype,     "MuGenMoType[NMus]/I");
	fMuTree->Branch("MuGenGMoType"  ,&fTmugmotype,    "MuGenGMoType[NMus]/I");
}

void MuonAnalysis::ResetTree(){
	fTrun     = 0;
	fTevent   = 0;
	fTlumisec = 0;
	fThbhenoiseflag = 0;
	fTnjets   = 0;
	fT_HLTMu9       = 0;
	fT_HLTMu11      = 0;
	fT_HLTMu15      = 0;
	fT_HLTDoubleMu0 = 0;
	fT_HLTDoubleMu3 = 0;
	fT_HLTJet30U    = 0;
	fT_HLTJet50U    = 0;
	fT_HLTJet70U    = 0;
	fT_HLTJet100U   = 0;	
	fTHT      = -999.99;
	fTMHT     = -999.99;
	fTPFMET   = -999.99;
	fTTCMET   = -999.99;
	fTSumET   = -999.99;
	fTMT      = -999.99;
	fTMinv    = -999.99;
	for(int i = 0; i < gMaxnjets; i++){
		fTjpt[i]  = -999.99;
		fTjeta[i] = -999.99;		
		fTjphi[i] = -999.99;
	}
	fTnmus = 0;
	for(int i = 0; i < gMaxnmus; i++){
		fTmupt        [i] = -999.99;
		fTmueta       [i] = -999.99;
		fTmuphi       [i] = -999.99;
		fTmucharge    [i] = -999;
		fTmutight     [i] = -999;
		fTmuiso       [i] = -999.99;
		fTmuisohyb    [i] = -999.99;
		fTDRjet       [i] = -999.99;
		fTDPhijet     [i] = -999.99;
		fTDRhardestjet[i] = -999.99;
		fTmud0        [i] = -999.99;
		fTmudz        [i] = -999.99;
		fTmud0bs      [i] = -999.99;
		fTmudzbs      [i] = -999.99;
		fTmuptE       [i] = -999.99;
		fTmuid        [i] = -999;
		fTmumoid      [i] = -999;
		fTmugmoid     [i] = -999;
		fTmutype      [i] = -999;
		fTmumotype    [i] = -999;
		fTmugmotype   [i] = -999;
	}
}
