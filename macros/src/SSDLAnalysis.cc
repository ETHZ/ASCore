#include "SSDLAnalysis.hh"
#include "helper/Monitor.hh"
#include "helper/PUWeight.h"


using namespace std;

const int SSDLAnalysis::fMaxNjets;
const int SSDLAnalysis::fMaxNmus;
const int SSDLAnalysis::fMaxNeles;

TString SSDLAnalysis::gBaseDir = "/shome/stiegerb/Workspace/cmssw/CMSSW_4_1_3/src/DiLeptonAnalysis/NTupleProducer/macros/";

SSDLAnalysis::SSDLAnalysis(TreeReader *tr): UserAnalysisBase(tr){
	//SetStyle();
	fHLTPaths.clear();

	fCutnames[0] = "All events";
	fCutnames[1] = " ... passes primary Vertex cuts";
	fCutnames[2] = " ... passes triggers (data only)";
	fCutnames[3] = " ... has at least one loose lepton (mc only)";


	fCounter.fill(fCutnames[0], 0.);
	fCounter.fill(fCutnames[1], 0.);
	fCounter.fill(fCutnames[2], 0.);
	fCounter.fill(fCutnames[3], 0.);
}

SSDLAnalysis::~SSDLAnalysis(){
}

void SSDLAnalysis::Begin(const char* filename){
	string pileupsrc = string(gBaseDir + "data_pileup.root");
	SetPileUpSrc(pileupsrc);
	ReadTriggers(gBaseDir + "HLTPaths_SSDL.dat");
	ReadPDGTable(gBaseDir + "pdgtable.txt");
	BookTree();
}

void SSDLAnalysis::End(){
	fOutputFile->cd();
	fAnalysisTree->Write();
	fOutputFile->Close();
	fCounter.print();
}

void SSDLAnalysis::ReadTriggers(const char* triggerfile){
	// Read in a bunch of HLT paths to be stored in the mini tree
	ifstream IN(triggerfile);
	char buffer[200];
	char StringValue[100];
	bool ok(false);
	
	
	while( IN.getline(buffer, 200, '\n') ){
		// ok = false;
		if (buffer[0] == '#') continue; // Skip lines commented with '#'
		if( !strcmp(buffer, "PATHSET")){
			HLTPathSet ps;

			IN.getline(buffer, 200, '\n');
			sscanf(buffer, "Name\t%s", StringValue);
			ps.name = TString(StringValue);
			while( IN.getline(buffer, 200, '\n') && IN.gcount() > 1 && buffer[0] != '#' ){
				// End declaration of paths of a set with either an empty line or a comment
				sscanf(buffer, "Path\t%s", StringValue);
				ps.paths.push_back(string(StringValue));
			}
			fHLTPathSets.push_back(ps);
			ok = true;
		}

		if(!ok) cout << "%% SSDLAnalysis::ReadTriggers ==> ERROR: Reading failed!" << endl;
	}
	if(fVerbose > 0){
		cout << "Adding Trigger path sets:" << endl;
		for(size_t i = 0; i < fHLTPathSets.size(); ++i){
			HLTPathSet ps = fHLTPathSets[i];
			cout << " " << ps.name << endl;
			for(size_t j = 0; j < ps.paths.size(); ++j){
				cout << "   > " << ps.paths[j] << endl;
			}
		}
		cout << "--------------" << endl;
	}

	fHLTResults.resize(fHLTPathSets.size());
	fHLTPrescales.resize(fHLTPathSets.size());
}

void SSDLAnalysis::BookTriggers(){
	for(int i = 0; i < fHLTPathSets.size(); i++){
		HLTPathSet ps = fHLTPathSets[i];
		TString prescalename = ps.name + "_PS";
		if(AddBranch(ps.name.Data(),      "I", &fHLTResults[i])   == false ) exit(-1);
		if(AddBranch(prescalename.Data(), "I", &fHLTPrescales[i]) == false ) exit(-1);
	}
}

bool SSDLAnalysis::FillTriggers(){
	// Returns OR of trigger results, i.e. true if ANY of them passed
	bool accept = false;
	if(fHLTPathSets.size() == 0 ) return false;

	for(int i = 0; i < fHLTPathSets.size(); i++){ // loop over path sets
		HLTPathSet ps = fHLTPathSets[i];

		for(size_t j = 0; j < ps.paths.size(); ++j){ // loop over paths
			fHLTResults[i]   = -1;
			fHLTPrescales[i] = -1;
			string path = ps.paths[j];

			if(GetHLTBit(path) == -1) continue; // Bit not found
			bool triggered   = GetHLTResult(path);
			fHLTResults[i]   = triggered ? 1:0;
			fHLTPrescales[i] = GetHLTPrescale(path);
			accept = accept || triggered;
			break; // break loop after the bit is found
		}
	}
	return accept;
}

const bool SSDLAnalysis::AddBranch( const char* name, const char* type, void* address, const char* size ){
	// This is copied from FillerBase.cc
  
  // Form input
  std::string fullname(name);
  
  std::string branchType(fullname);
  if ( size ) // Size needs to be pre-fixed
    branchType += "[" + std::string(size) + "]";
  branchType += "/"+std::string(type);

  // Declare branch
  TBranch* b = fAnalysisTree->Branch(fullname.c_str(),address,branchType.c_str());

  return !(b==0); // return 1 if branch was successfully created, 0 otherwise

}

void SSDLAnalysis::BookTree(){
	fOutputFile->cd();
	fAnalysisTree = new TTree("Analysis", "AnalysisTree");

    // run/sample properties
	fAnalysisTree->Branch("Run",              &fTRunNumber,         "Run/I");
	fAnalysisTree->Branch("Event",            &fTEventNumber,       "Event/I");
	fAnalysisTree->Branch("LumiSec",          &fTLumiSection,       "LumiSec/I");

	// HLT triggers
	BookTriggers();
	
	// event properties
	fAnalysisTree->Branch("Rho",           &fTrho,       "Rho/F");
	fAnalysisTree->Branch("NVrtx",         &fTnvrtx,     "NVrtx/I");
	fAnalysisTree->Branch("PUWeight",      &fTpuweight,  "PUWeight/F");

	
	// single-muon properties
	fAnalysisTree->Branch("NMus"          ,&fTnqmus,          "NMus/I");
	fAnalysisTree->Branch("MuPt"          ,&fTmupt,           "MuPt[NMus]/F");
	fAnalysisTree->Branch("MuEta"         ,&fTmueta,          "MuEta[NMus]/F");
	fAnalysisTree->Branch("MuPhi"         ,&fTmuphi,          "MuPhi[NMus]/F");
	fAnalysisTree->Branch("MuCharge"      ,&fTmucharge,       "MuCharge[NMus]/I");
	fAnalysisTree->Branch("MuIso"         ,&fTmuiso,          "MuIso[NMus]/F");
	fAnalysisTree->Branch("MuD0"          ,&fTmud0,           "MuD0[NMus]/F");
	fAnalysisTree->Branch("MuDz"          ,&fTmudz,           "MuDz[NMus]/F");
	fAnalysisTree->Branch("MuPtE"         ,&fTmuptE,          "MuPtE[NMus]/F");
	fAnalysisTree->Branch("MuGenID"       ,&fTmuid,           "MuGenID[NMus]/I");
	fAnalysisTree->Branch("MuGenMID"      ,&fTmumoid,         "MuGenMID[NMus]/I");
	fAnalysisTree->Branch("MuGenGMID"     ,&fTmugmoid,        "MuGenGMID[NMus]/I");
	fAnalysisTree->Branch("MuGenType"     ,&fTmutype,         "MuGenType[NMus]/I");
	fAnalysisTree->Branch("MuGenMType"    ,&fTmumotype,       "MuGenMType[NMus]/I");
	fAnalysisTree->Branch("MuGenGMType"   ,&fTmugmotype,      "MuGenGMType[NMus]/I");
	fAnalysisTree->Branch("MuMT"          ,&fTmuMT,           "MuMT[NMus]/F");

	// single-electron properties
	fAnalysisTree->Branch("NEls",                   &fTnqels,               "NEls/I");
	fAnalysisTree->Branch("ElCharge",               &fTElcharge,            "ElCh[NEls]/I");
	fAnalysisTree->Branch("ElChIsCons",             &fTElChargeIsCons,      "ElChIsCons[NEls]/I");
	fAnalysisTree->Branch("ElPt",                   &fTElpt,                "ElPt[NEls]/F");
	fAnalysisTree->Branch("ElEta",                  &fTEleta,               "ElEta[NEls]/F");
	fAnalysisTree->Branch("ElPhi",                  &fTElphi,               "ElPhi[NEls]/F");
	fAnalysisTree->Branch("ElD0",                   &fTEld0,                "ElD0[NEls]/F");
	fAnalysisTree->Branch("ElD0Err",                &fTElD0Err,             "ElD0Err[NEls]/F");
	fAnalysisTree->Branch("ElDz",                   &fTEldz,                "ElDz[NEls]/F");
	fAnalysisTree->Branch("ElDzErr",                &fTElDzErr,             "ElDzErr[NEls]/F");
	fAnalysisTree->Branch("ElRelIso",               &fTElRelIso,            "ElRelIso[NEls]/F");
	fAnalysisTree->Branch("ElEcalRecHitSumEt",      &fTElEcalRecHitSumEt,   "ElEcalRecHitSumEt[NEls]/F");
	fAnalysisTree->Branch("ElIsGoodElId_WP80",      &fTElIsGoodElId_WP80,   "ElIsGoodElId_WP80[NEls]/I");
	fAnalysisTree->Branch("ElIsGoodElId_WP90",      &fTElIsGoodElId_WP90,   "ElIsGoodElId_WP90[NEls]/I");
	fAnalysisTree->Branch("ElGenID",                &fTElGenID,             "ElGenID[NEls]/I");
	fAnalysisTree->Branch("ElGenMID",               &fTElGenMID,            "ElGenMID[NEls]/I");
	fAnalysisTree->Branch("ElGenGMID",              &fTElGenGMID,           "ElGenGMID[NEls]/I");
	fAnalysisTree->Branch("ElGenType",              &fTElGenType,           "ElGenType[NEls]/I");
	fAnalysisTree->Branch("ElGenMType",             &fTElGenMType,          "ElGenMType[NEls]/I");
	fAnalysisTree->Branch("ElGenGMType",            &fTElGenGMType,         "ElGenGMType[NEls]/I");
	fAnalysisTree->Branch("ElMT",                   &fTElMT,                "ElMT[NEls]/F");

	// jet-MET properties
	fAnalysisTree->Branch("tcMET",         &fTtcMET,    "tcMET/F");
	fAnalysisTree->Branch("tcMETPhi",      &fTtcMETphi, "tcMETPhi/F");
	fAnalysisTree->Branch("pfMET",         &fTpfMET,    "pfMET/F");
	fAnalysisTree->Branch("pfMETPhi",      &fTpfMETphi, "pfMETPhi/F");
	fAnalysisTree->Branch("NJets",         &fTnqjets,   "NJets/I");
	fAnalysisTree->Branch("JetPt",         &fTJetpt,    "JetPt[NJets]/F");
	fAnalysisTree->Branch("JetEta",        &fTJeteta,   "JetEta[NJets]/F");
	fAnalysisTree->Branch("JetPhi",        &fTJetphi,   "JetPhi[NJets]/F");
	fAnalysisTree->Branch("JetSSVHPBTag",  &fTJetbtag,  "JetSSVHPBTag[NJets]/F"); // tight WP: > 2.
	fAnalysisTree->Branch("JetArea",       &fTJetArea,  "JetArea[NJets]/F");
}

void SSDLAnalysis::Analyze(){
	fCounter.fill(fCutnames[0]);
	// initial event selection: good event trigger, good primary vertex...
	if( !IsGoodEvent() ) return;
	fCounter.fill(fCutnames[1]);
	ResetTree();
	
	// Trigger selection
	// if(fIsData && FillTriggers(fHLTPaths) == false) return;
	FillTriggers();
	fCounter.fill(fCutnames[2]);

	// Do object selections
	vector<int> selectedMuInd  = MuonSelection(&UserAnalysisBase::IsGoodBasicMu);
	vector<int> selectedElInd  = ElectronSelection(&UserAnalysisBase::IsLooseEl);
	vector<int> selectedJetInd = PFJetSelection(40., 2.5, &UserAnalysisBase::IsGoodBasicPFJet);
	fTnqmus  = std::min((int) selectedMuInd .size(), fMaxNmus);
	fTnqels  = std::min((int) selectedElInd .size(), fMaxNeles);
	fTnqjets = std::min((int) selectedJetInd.size(), fMaxNjets);
	
	// Require at least one loose lepton
	// if(!fIsData && (fTnqmus + fTnqels) < 1 ) return;
	if( (fTnqmus + fTnqels) < 1 ) return;
	fCounter.fill(fCutnames[3]);

	// event and run info
	fTRunNumber   = fTR->Run;
	fTEventNumber = fTR->Event;
	fTLumiSection = fTR->LumiSection;

	// Dump basic jet and MET properties
	int jetindex(-1);
	int nqjets = selectedJetInd.size();
	for(int ind=0; ind<std::min(nqjets,fMaxNjets); ind++){
		jetindex = selectedJetInd[ind];
		// dump properties
		fTJetpt  [ind] = fTR->JPt [jetindex];
		fTJeteta [ind] = fTR->JEta[jetindex];
		fTJetphi [ind] = fTR->JPhi[jetindex];
		fTJetbtag[ind] = fTR->JbTagProbTkCntHighPur[jetindex];
		fTJetArea[ind] = fTR->JArea[jetindex];
	}
	// get METs
	fTtcMET     = fTR->TCMET;
	fTtcMETphi  = fTR->TCMETphi;
	fTpfMET     = fTR->PFMET;
	fTpfMETphi  = fTR->PFMETphi;

	// PU correction
	fTrho   = fTR->Rho;
	fTnvrtx = fTR->NVrtx;
	if(!fIsData) fTpuweight = GetPUWeight(fTR->PUnumInteractions);
	else fTpuweight = 1.;

	int nqmus = selectedMuInd.size();
	for(int i = 0; i < std::min(nqmus, fMaxNmus); ++i){
		int index = selectedMuInd[i];
		fTmupt       [i] = fTR->MuPt[index];
		fTmueta      [i] = fTR->MuEta[index];
		fTmuphi      [i] = fTR->MuPhi[index];
		fTmucharge   [i] = fTR->MuCharge[index];
		fTmuiso      [i] = fTR->MuRelIso03[index];
		fTmud0       [i] = fTR->MuD0PV[index];
		fTmudz       [i] = fTR->MuDzPV[index];
		fTmuptE      [i] = fTR->MuPtE[index];
		
		if(fIsData == false){ // mc truth information
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
		} else{
			fTmuid     [i] = -888;
			fTmumoid   [i] = -888;
			fTmugmoid  [i] = -888;
			fTmutype   [i] = -888;
			fTmumotype [i] = -888;
			fTmugmotype[i] = -888;
		}
		
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
		fTElcharge           [ind] = fTR->ElCharge                [elindex];
		fTElChargeIsCons     [ind] = fTR->ElCInfoIsGsfCtfScPixCons[elindex];
		fTElpt               [ind] = fTR->ElPt                   [elindex];
		fTEleta              [ind] = fTR->ElEta                  [elindex];
		fTElphi              [ind] = fTR->ElPhi                  [elindex];
		fTEld0               [ind] = fTR->ElD0PV                 [elindex];
		fTElD0Err            [ind] = fTR->ElD0E                  [elindex];
		fTEldz               [ind] = fTR->ElDzPV                 [elindex];
		fTElDzErr            [ind] = fTR->ElDzE                  [elindex];
		fTElRelIso           [ind] = relElIso(elindex); // correct by 1 GeV in ecal for barrel
		fTElEcalRecHitSumEt  [ind] = fTR->ElDR03EcalRecHitSumEt  [elindex];
		
		if(fIsData == false){ // mc truth information		
			fTElGenID  [ind] = fTR->ElGenID  [elindex];
			fTElGenMID [ind] = fTR->ElGenMID [elindex];
			fTElGenGMID[ind] = fTR->ElGenGMID[elindex];
			pdgparticle el, emo, egmo;
			GetPDGParticle(el,   abs(fTR->ElGenID  [elindex]));
			GetPDGParticle(emo,  abs(fTR->ElGenMID [elindex]));
			GetPDGParticle(egmo, abs(fTR->ElGenGMID[elindex]));
			fTElGenType  [ind] = el.get_type();
			fTElGenMType [ind] = emo.get_type();
			fTElGenGMType[ind] = egmo.get_type();
		} else{
			fTElGenID    [ind] = -888;
			fTElGenMID   [ind] = -888;
			fTElGenGMID  [ind] = -888;
			fTElGenType  [ind] = -888;
			fTElGenMType [ind] = -888;
			fTElGenGMType[ind] = -888;   
		}
		
		// Calculate mT:
		TLorentzVector pel;
		pel.SetXYZM(fTR->ElPx[elindex], fTR->ElPy[elindex], fTR->ElPz[elindex], 0.0005);
		double ETlept = sqrt(pel.M2() + pel.Perp2());
		double METpx  = fTR->PFMETpx;
		double METpy  = fTR->PFMETpy;
		fTElMT[ind]   = sqrt( 2*fTR->PFMET*ETlept - pel.Px()*METpx - pel.Py()*METpy );
		
		fTElIsGoodElId_WP80[ind] = IsGoodElId_WP80(elindex);
		fTElIsGoodElId_WP90[ind] = IsGoodElId_WP90(elindex);
	}

	fAnalysisTree->Fill();
}

void SSDLAnalysis::ResetTree(){
	// sample/run
	fTRunNumber                  = 0;
	fTEventNumber                = 0;
	fTLumiSection                = 0;

	for(size_t i = 0; i < fHLTPathSets.size(); ++i){
		fHLTResults[i]   -2;
		fHLTPrescales[i] -2;
	}
	
	fTrho      = -999.99;
	fTnvrtx    = -999;
	fTpuweight = -999.99;
	
	// muon properties
	fTnqmus = 0;
	for(int i = 0; i < fMaxNmus; i++){
		fTmupt          [i] = -999.99;
		fTmueta         [i] = -999.99;
		fTmuphi         [i] = -999.99;
		fTmucharge      [i] = -999;
		fTmuiso         [i] = -999.99;
		fTmud0          [i] = -999.99;
		fTmudz          [i] = -999.99;
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
		fTElpt              [i] = -999.99;
		fTEleta             [i] = -999.99;
		fTElphi             [i] = -999.99;
		fTEld0              [i] = -999.99;
		fTElD0Err           [i] = -999.99;
		fTEldz              [i] = -999.99;
		fTElDzErr           [i] = -999.99;
		fTElRelIso          [i] = -999.99;
		fTElEcalRecHitSumEt [i] = -999.99;
		fTElIsGoodElId_WP80 [i] = -999;
		fTElIsGoodElId_WP90 [i] = -999;
		fTElMT              [i] = -999.99;
		fTElGenID           [i] = -999;
		fTElGenMID          [i] = -999;
		fTElGenGMID         [i] = -999;
		fTElGenType         [i] = -999;
		fTElGenMType        [i] = -999;
		fTElGenGMType       [i] = -999;
	}

	// jet-MET properties
	fTnqjets = 0;
	for(int i = 0; i < fMaxNjets; i++){
		fTJetpt [i]  = -999.99;
		fTJeteta[i]  = -999.99;
		fTJetphi[i]  = -999.99;
		fTJetbtag[i] = -999.99;
		fTJetArea[i] = -999.99;
	}
	fTtcMET      = -999.99;
	fTtcMETphi   = -999.99;
	fTpfMET      = -999.99;
	fTpfMETphi   = -999.99;
}
