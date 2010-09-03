#include "DiLeptonAnalysis.hh"
#include "base/TreeReader.hh"

using namespace std;

DiLeptonAnalysis::DiLeptonAnalysis(TreeReader *tr) : UserAnalysisBase(tr){
}

DiLeptonAnalysis::~DiLeptonAnalysis(){
}

void DiLeptonAnalysis::Begin(const char* filename){
	fDiLepTreeFile_ = new TFile(fOutputDir + TString(filename), "RECREATE");
	fDiLepTree_ = new TTree("Analysis", "NanoAnalysisTree");
	fDiLepTree_->Branch("Run",           &fTRunNumber,      "Run/I");
	fDiLepTree_->Branch("Event",         &fTEventNumber,    "Event/I");
	fDiLepTree_->Branch("LumiSec",       &fTLumiSection,    "LumiSec/I");
	fDiLepTree_->Branch("NQMus",         &fTNqualmu,        "NQMus/I");
	fDiLepTree_->Branch("Mu1Ch",         &fTMu1charge,      "Mu1Ch/I");
	fDiLepTree_->Branch("Mu2Ch",         &fTMu2charge,      "Mu2Ch/I");
	fDiLepTree_->Branch("Mu1Pt",         &fTMu1pt,          "Mu1Pt/D");
	fDiLepTree_->Branch("Mu2Pt",         &fTMu2pt,          "Mu2Pt/D");
	fDiLepTree_->Branch("Mu1Eta",        &fTMu1eta,         "Mu1Eta/D");
	fDiLepTree_->Branch("Mu2Eta",        &fTMu2eta,         "Mu2Eta/D");
	fDiLepTree_->Branch("Mu1Iso",        &fTMu1iso,         "Mu1Iso/D");
	fDiLepTree_->Branch("Mu2Iso",        &fTMu2iso,         "Mu2Iso/D");
	fDiLepTree_->Branch("Mu1D0",         &fTMu1d0,          "Mu1D0/D");
	fDiLepTree_->Branch("Mu2D0",         &fTMu2d0,          "Mu2D0/D");
	fDiLepTree_->Branch("Mu1NTkHits",    &fTMu1ntkhits,     "Mu1NTkHits/D");
	fDiLepTree_->Branch("Mu2NTkHits",    &fTMu2ntkhits,     "Mu2NTkHits/D");
	fDiLepTree_->Branch("MuMInv",        &fTMuminv,         "MuMInv/D");
	fDiLepTree_->Branch("MuMT2_50",      &fTMumt2_50,       "MuMT2_50/D");
	fDiLepTree_->Branch("MuMT2_100",     &fTMumt2_100,      "MuMT2_100/D");
	fDiLepTree_->Branch("NQEls",         &fTNqualel,        "NQEls/I");
	fDiLepTree_->Branch("El1Ch",         &fTEl1charge,      "El1Ch/I");
	fDiLepTree_->Branch("El2Ch",         &fTEl2charge,      "El2Ch/I");
	fDiLepTree_->Branch("El1Pt",         &fTEl1pt,          "El1Pt/D");
	fDiLepTree_->Branch("El2Pt",         &fTEl2pt,          "El2Pt/D");
	fDiLepTree_->Branch("El1Eta",        &fTEl1eta,         "El1Eta/D");
	fDiLepTree_->Branch("El2Eta",        &fTEl2eta,         "El2Eta/D");
	fDiLepTree_->Branch("El1Iso",        &fTEl1iso,         "El1Iso/D");
	fDiLepTree_->Branch("El2Iso",        &fTEl2iso,         "El2Iso/D");
	fDiLepTree_->Branch("El1D0",         &fTEl1d0,          "El1D0/D");
	fDiLepTree_->Branch("El2D0",         &fTEl2d0,          "El2D0/D");
	fDiLepTree_->Branch("ElMInv",        &fTElminv,         "ElMInv/D");
	fDiLepTree_->Branch("ElMT2_50",      &fTElmt2_50,       "ElMT2_50/D");
	fDiLepTree_->Branch("ElMT2_100",     &fTElmt2_100,      "ElMT2_100/D");
}


void DiLeptonAnalysis::Analyze(){
	Reset();
	// Select events with either 2 muons or 2 electrons
	if( fTR->NMus < 2 && fTR->NEles < 2 ) return;

	vector<int> qualMuInd;
	for(size_t imu = 0; imu < fTR->NMus; ++imu){
		// Muon selection
		if(fTR->MuIsGlobalMuon[imu] == 0) continue;
		if(fTR->MuGood[imu] != 0) continue;
		if(fTR->MuPt[imu] < 10) continue;
		if(fabs(fTR->MuEta[imu]) > 2.4) continue;
		qualMuInd.push_back(imu);
	}
	vector<int> qualElInd;
	for(size_t iel = 0; iel < fTR->NEles; ++iel){
		// Electron selection
		if(fTR->ElGood[iel] != 0) continue;
		if(fTR->ElPt[iel] < 10) continue;
		if(fabs(fTR->ElEta[iel]) > 2.4) continue;		
		qualElInd.push_back(iel);
	}

	int nqmus = qualMuInd.size();
	int nqels = qualElInd.size();

	// Select events with either 2 qualified muons or electrons
	if(nqmus < 2 && nqels < 2) return;

	fTNqualmu = nqmus;
	fTNqualel = nqels;

	fMT2 = new Davismt2();
	double pa[3], pb[3], pmiss[3];

	int lep1index(-1), lep2index(-1);
	// DiMuons:
	if(nqmus > 1){
		// Find the two hardest muons
		double maxmupt = 0.;
		for(size_t i = 0; i < nqmus; ++i){
			int index = qualMuInd[i];
			if(fTR->MuPt[index] < maxmupt) continue;
			maxmupt = fTR->MuPt[index];
			lep1index = index;
		}
		maxmupt = 0.;
		for(size_t i = 0; i < nqmus; ++i){
			int index = qualMuInd[i];
			if(index == lep1index) continue;
			if(fTR->MuPt[index] < maxmupt) continue;
			maxmupt = fTR->MuPt[index];
			lep2index = index;
		}

		fTMu1pt       = fTR->MuPt[lep1index];
		fTMu2pt       = fTR->MuPt[lep2index];
		fTMu1charge   = fTR->MuCharge[lep1index];
		fTMu2charge   = fTR->MuCharge[lep2index];
		fTMu1eta      = fTR->MuEta[lep1index];
		fTMu2eta      = fTR->MuEta[lep2index];
		fTMu1iso      = fTR->MuRelIso03[lep1index];
		fTMu2iso      = fTR->MuRelIso03[lep2index];
		fTMu1d0       = fTR->MuD0BS[lep1index];
		fTMu2d0       = fTR->MuD0BS[lep2index];
		fTMu1ntkhits  = fTR->MuNTkHits[lep1index];
		fTMu2ntkhits  = fTR->MuNTkHits[lep2index];

		// Calculate invariant mass of the pair
		TLorentzVector p1(fTR->MuPx[lep1index], fTR->MuPy[lep1index], fTR->MuPz[lep1index], fTR->MuE[lep1index]);
		TLorentzVector p2(fTR->MuPx[lep2index], fTR->MuPy[lep2index], fTR->MuPz[lep2index], fTR->MuE[lep2index]);
		fTMuminv = (p1+p2).Mag();

		pa[0] = 0.; pb[0] = 0.; pmiss[0] = 0.;
		pa[1] = fTR->MuPx[lep1index];
		pa[2] = fTR->MuPz[lep1index];
		pb[1] = fTR->MuPx[lep2index];
		pb[2] = fTR->MuPy[lep2index];
		pmiss[1] = fTR->TCMETpx;
		pmiss[2] = fTR->TCMETpy;
		fMT2->set_momenta(pa, pb, pmiss);
		fMT2->set_mn(50.);
		fTMumt2_50 = fMT2->get_mt2();
		fMT2->set_mn(100.);
		fTMumt2_100 = fMT2->get_mt2();
	}else if(nqmus > 0){ // Meaning there is only one muon
		fTMu1pt       = fTR->MuPt[qualMuInd[0]];
		fTMu1charge   = fTR->MuCharge[qualMuInd[0]];
		fTMu1eta      = fTR->MuEta[qualMuInd[0]];
		fTMu1iso      = fTR->MuRelIso03[qualMuInd[0]];
		fTMu1d0       = fTR->MuD0BS[qualMuInd[0]];
		fTMu1ntkhits  = fTR->MuNTkHits[qualMuInd[0]];
	}

	lep1index = -1;
	lep2index = -1;

	// DiElectrons
	if(nqels > 1){
		// Find the two hardest muons
		double maxelpt = 0.;
		for(size_t i = 0; i < nqels; ++i){
			int index = qualElInd[i];
			if(fTR->ElPt[index] < maxelpt) continue;
			maxelpt = fTR->ElPt[index];
			lep1index = index;
		}
		maxelpt = 0.;
		for(size_t i = 0; i < nqels; ++i){
			int index = qualElInd[i];
			if(index == lep1index) continue;
			if(fTR->ElPt[index] < maxelpt) continue;
			maxelpt = fTR->ElPt[index];
			lep2index = index;
		}

		fTEl1charge   = fTR->ElCharge[lep1index];
		fTEl2charge   = fTR->ElCharge[lep2index];
		fTEl1pt       = fTR->ElPt[lep1index];
		fTEl2pt       = fTR->ElPt[lep2index];
		fTEl1eta      = fTR->ElEta[lep1index];
		fTEl2eta      = fTR->ElEta[lep2index];
		fTEl1iso      = fTR->ElRelIso04[lep1index];
		fTEl2iso      = fTR->ElRelIso04[lep2index];
		fTEl1d0       = fTR->ElD0BS[lep1index];
		fTEl2d0       = fTR->ElD0BS[lep2index];

		// Calculate invariant mass of the pair
		TLorentzVector p1(fTR->ElPx[lep1index], fTR->ElPy[lep1index], fTR->ElPz[lep1index], fTR->ElE[lep1index]);
		TLorentzVector p2(fTR->ElPx[lep2index], fTR->ElPy[lep2index], fTR->ElPz[lep2index], fTR->ElE[lep2index]);
		fTElminv = (p1+p2).Mag();

		pa[0] = 0.; pb[0] = 0.; pmiss[0] = 0.;
		pa[1] = fTR->ElPx[lep1index];
		pa[2] = fTR->ElPy[lep1index];
		pb[1] = fTR->ElPx[lep2index];
		pb[2] = fTR->ElPy[lep2index];
		pmiss[1] = fTR->TCMETpx;
		pmiss[2] = fTR->TCMETpy;
		fMT2->set_momenta(pa, pb, pmiss);
		fMT2->set_mn(50.);
		fTElmt2_50 = fMT2->get_mt2();
		fMT2->set_mn(100.);
		fTElmt2_100 = fMT2->get_mt2();
	}else if (nqels>0){ // Meaning there is only one electron
		fTEl1pt       = fTR->ElPt[qualElInd[0]];
		fTEl1charge   = fTR->ElCharge[qualElInd[0]];
		fTEl1eta      = fTR->ElEta[qualElInd[0]];
		fTEl1iso      = fTR->ElRelIso04[qualElInd[0]];
		fTEl1d0       = fTR->ElD0BS[qualElInd[0]];
	}
	fDiLepTree_->Fill();
	delete fMT2;
}

void DiLeptonAnalysis::Reset(){
	fTRunNumber   = -999;
	fTEventNumber = -999;
	fTLumiSection = -999;
	fTNqualmu     = -999;
	fTMu1charge   = -999;
	fTMu2charge   = -999;
	fTMu1pt       = -999.99;
	fTMu2pt       = -999.99;
	fTMu1eta      = -999.99;
	fTMu2eta      = -999.99;
	fTMu1iso      = -999.99;
	fTMu2iso      = -999.99;
	fTMu1d0       = -999.99;
	fTMu2d0       = -999.99;
	fTMu1ntkhits  = -999.99;
	fTMu2ntkhits  = -999.99;
	fTMuminv      = -999.99;
	fTMumt2_50    = -999.99;
	fTMumt2_100   = -999.99;
	fTNqualel     = -999;
	fTEl1charge   = -999;
	fTEl2charge   = -999;
	fTEl1pt       = -999.99;
	fTEl2pt       = -999.99;
	fTEl1eta      = -999.99;
	fTEl2eta      = -999.99;
	fTEl1iso      = -999.99;
	fTEl2iso      = -999.99;
	fTEl1d0       = -999.99;
	fTEl2d0       = -999.99;
	fTElminv      = -999.99;
	fTElmt2_50    = -999.99;
	fTElmt2_100   = -999.99;
}

void DiLeptonAnalysis::End(){
	fDiLepTreeFile_->cd();
	fDiLepTree_->Write();
	fDiLepTreeFile_->Close();
}

