#include "MuonAnalysis.hh"
#include "base/TreeReader.hh"
#include "helper/Utilities.hh"

#include <TH1D.h>
#include <TH2D.h>
#include <TVector3.h>

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
	int muind(-1);
	if(!SingleMuonSelection(muind)) return;
	FillMuonTree(muind);
}

void MuonAnalysis::AnalyzeDi(){
	int muind1(-1), muind2(-1);
	if(!DiMuonSelection(muind1, muind2, -1)) return;
	FillMuonTree(muind1, muind2);
}

void MuonAnalysis::AnalyzeSS(){
	int muind1(-1), muind2(-1);
	if(!SSDiMuonSelection(muind1, muind2)) return;
	FillMuonTree(muind1, muind2);
}

void MuonAnalysis::End(){
	fOutputFile->cd();
	fMuTree->Write();
	fOutputFile->Close();
}

void MuonAnalysis::FillMuonTree(int ind1, int ind2){
	ResetTree();

	// Jet Variables
	vector<int> goodJets = JetSelection();
	fTnjets = goodJets.size();
	double ht = 0;
	TVector3 mht(0.,0.,0.);
	for(size_t i = 0; i < fTnjets; ++i){
		int index = goodJets[i];
		fTjpt [i] = fTR->JPt[index];
		fTjeta[i] = fTR->JEta[index];
		fTjphi[i] = fTR->JPhi[index];
		ht += fTR->JPt[index];
		mht(0) += fTR->JPx[index];
		mht(1) += fTR->JPy[index];
		mht(2) += fTR->JPz[index];
	}
	fTHT    = ht;
	fTMHT   = mht.Mag();
	fTMET   = fTR->PFMET;
	fTSumET = fTR->SumEt;
	
	// Muon Variables
	int nmus(0);
	for(size_t imu = 0; imu < fTR->NMus; ++imu) if(IsGoodBasicMu(imu)) nmus++;
	int indices[2] = {ind1, ind2};
	for(size_t i = 0; i < 2; ++i){
		int index = indices[i];
		if(index < 0) continue;
		fTmupt      [i] = fTR->MuPt[index];
		fTmueta     [i] = fTR->MuEta[index];
		fTmuphi     [i] = fTR->MuPhi[index];
		fTmucharge  [i] = fTR->MuCharge[index];
		if(IsTightMu(index))        fTmutight[i] = 1;
		if(IsLooseNoTightMu(index)) fTmutight[i] = 0;
		fTmuiso     [i] = fTR->MuRelIso03[index];
		fTmucalocomp[i] = fTR->MuCaloComp[index];
		fTmusegmcomp[i] = fTR->MuSegmComp[index];
		fTmuouterrad[i] = fTR->MuOutPosRadius[index];
		fTmunchi2   [i] = fTR->MuNChi2[index];
		fTmuntkhits [i] = fTR->MuNTkHits[index];
		fTmud0      [i] = fTR->MuD0PV[index];
		fTmudz      [i] = fTR->MuDzPV[index];
		fTmuptE     [i] = fTR->MuPtE[index];

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
	fTnmus = nmus;
	
	fMuTree->Fill();
}

void MuonAnalysis::BookTree(){
	fOutputFile->cd();
	fMuTree = new TTree("MuonAnalysis", "MuonFakeAnalysisTree1");
	fMuTree->Branch("Run"           ,&fTrun,       "Run/I");
	fMuTree->Branch("Event"         ,&fTevent,     "Event/I");
	fMuTree->Branch("LumiSection"   ,&fTlumisec,   "LumiSection/I");
	fMuTree->Branch("NJets"         ,&fTnjets,     "NJets/I");
	fMuTree->Branch("JPt"           ,&fTjpt,       "JPt[NJets]/D");
	fMuTree->Branch("JEta"          ,&fTjeta,      "JEta[NJets]/D");
	fMuTree->Branch("JPhi"          ,&fTjphi,      "JPhi[NJets]/D");
	fMuTree->Branch("HT"            ,&fTHT,        "HT/D");
	fMuTree->Branch("MHT"           ,&fTMHT,       "MHT/D");
	fMuTree->Branch("MET"           ,&fTMET,       "MET/D");
	fMuTree->Branch("SumET"         ,&fTSumET,     "SumET/D");
	fMuTree->Branch("NMus"          ,&fTnmus,      "NMus/I");
	fMuTree->Branch("MuPt"          ,&fTmupt,      "MuPt[2]/D");
	fMuTree->Branch("MuEta"         ,&fTmueta,     "MuEta[2]/D");
	fMuTree->Branch("MuPhi"         ,&fTmuphi,     "MuPhi[2]/D");
	fMuTree->Branch("MuCharge"      ,&fTmucharge,  "MuCharge[2]/I");
	fMuTree->Branch("MuTight"       ,&fTmutight,   "MuTight[2]/I");
	fMuTree->Branch("MuIso"         ,&fTmuiso,     "MuIso[2]/D");
	fMuTree->Branch("MuDRJet"       ,&fTDRjet,     "MuDRJet[2]/D");
	fMuTree->Branch("MuDPhiJet"     ,&fTDPhijet,   "MuDPhiJet[2]/D");
	fMuTree->Branch("MuCaloComp"    ,&fTmucalocomp,"MuCaloComp[2]/D");
	fMuTree->Branch("MuSegmComp"    ,&fTmusegmcomp,"MuSegmComp[2]/D");
	fMuTree->Branch("MuOuterRad"    ,&fTmuouterrad,"MuOuterRad[2]/D");
	fMuTree->Branch("MuNChi2"       ,&fTmunchi2,   "MuNChi2[2]/D");
	fMuTree->Branch("MuNTkHits"     ,&fTmuntkhits, "MuNTkHits[2]/I");
	fMuTree->Branch("MuD0"          ,&fTmud0,      "MuD0[2]/D");
	fMuTree->Branch("MuDz"          ,&fTmudz,      "MuDz[2]/D");
	fMuTree->Branch("MuPtE"         ,&fTmuptE,     "MuPtE[2]/D");
	fMuTree->Branch("MuGenID"       ,&fTmuid,      "MuGenID[2]/I");
	fMuTree->Branch("MuGenMoID"     ,&fTmumoid,    "MuGenMoID[2]/I");
	fMuTree->Branch("MuGenGMoID"    ,&fTmugmoid,   "MuGenGMoID[2]/I");
	fMuTree->Branch("MuGenType"     ,&fTmutype,    "MuGenType[2]/I");
	fMuTree->Branch("MuGenMoType"   ,&fTmumotype,  "MuGenMoType[2]/I");
	fMuTree->Branch("MuGenGMoType"  ,&fTmugmotype, "MuGenGMoType[2]/I");
}

void MuonAnalysis::ResetTree(){
	fTrun     = 0;
	fTevent   = 0;
	fTlumisec = 0;
	fTnjets   = 0;
	fTHT      = -999.99;
	fTMHT     = -999.99;
	fTMET     = -999.99;
	fTSumET   = -999.99;
	for(int i = 0; i<gMaxnjets; i++){
		fTjpt[i]  = -999.99;
		fTjeta[i] = -999.99;		
		fTjphi[i] = -999.99;
	}
	fTnmus = 0;
	for(int i = 0; i<2; i++){
		fTmupt      [i] = -999.99;
		fTmueta     [i] = -999.99;
		fTmuphi     [i] = -999.99;
		fTmucharge  [i] = -999;
		fTmutight   [i] = -999;
		fTmuiso     [i] = -999.99;
		fTDRjet     [i] = -999.99;
		fTDPhijet   [i] = -999.99;
		fTmucalocomp[i] = -999.99;
		fTmusegmcomp[i] = -999.99;
		fTmuouterrad[i] = -999.99;
		fTmunchi2   [i] = -999.99;
		fTmuntkhits [i] = -999;
		fTmud0      [i] = -999.99;
		fTmudz      [i] = -999.99;
		fTmuptE     [i] = -999.99;
		fTmuid      [i] = -999;
		fTmumoid    [i] = -999;
		fTmugmoid   [i] = -999;
		fTmutype    [i] = -999;
		fTmumotype  [i] = -999;
		fTmugmotype [i] = -999;
	}
}
