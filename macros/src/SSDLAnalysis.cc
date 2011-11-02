#include "SSDLAnalysis.hh"
#include "helper/Monitor.hh"
#include "helper/PUWeight.h"


using namespace std;

const int SSDLAnalysis::fMaxNjets;
const int SSDLAnalysis::fMaxNmus;
const int SSDLAnalysis::fMaxNeles;

TString SSDLAnalysis::gBaseDir = "/shome/stiegerb/Workspace/cmssw/CMSSW_4_1_3/src/DiLeptonAnalysis/NTupleProducer/macros/";
// TString SSDLAnalysis::gBaseDir = "/home/stiegerb/Workspace/cmssw/CMSSW_4_1_3/src/DiLeptonAnalysis/NTupleProducer/macros/";

//____________________________________________________________________________
SSDLAnalysis::SSDLAnalysis(TreeReader *tr): UserAnalysisBase(tr){
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

//____________________________________________________________________________
void SSDLAnalysis::Begin(const char* filename){
	// string pileupsrc = string(gBaseDir + "data_pileup.root");
	// SetPileUpSrc(pileupsrc);
	// ReadTriggers(gBaseDir + "HLTPaths_SSDL.dat");
	// ReadPDGTable(gBaseDir + "pdgtable.txt");
	SetPileUpSrc("/shome/stiegerb/Workspace/cmssw/CMSSW_4_1_3/src/DiLeptonAnalysis/NTupleProducer/macros/data_pileup.root");
	ReadTriggers("/shome/stiegerb/Workspace/cmssw/CMSSW_4_1_3/src/DiLeptonAnalysis/NTupleProducer/macros/HLTPaths_SSDL.dat");
	ReadPDGTable("/shome/stiegerb/Workspace/cmssw/CMSSW_4_1_3/src/DiLeptonAnalysis/NTupleProducer/macros/pdgtable.txt");
	BookTree();
	if(!fIsData && fDoFillEffTree) BookEffTree();
}

//____________________________________________________________________________
void SSDLAnalysis::End(){
	fOutputFile->cd();
	fAnalysisTree->Write();
	if(fDoFillEffTree && !fIsData) fLepEffTree->Write();
	fOutputFile->Close();
	fCounter.print();
}

//____________________________________________________________________________
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
void SSDLAnalysis::AddTriggerBranches(){
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

//____________________________________________________________________________
void SSDLAnalysis::BookTree(){
	fOutputFile->cd();
	fAnalysisTree = new TTree("Analysis", "AnalysisTree");

    // run/sample properties
	fAnalysisTree->Branch("Run",              &fTRunNumber,         "Run/I");
	fAnalysisTree->Branch("Event",            &fTEventNumber,       "Event/I");
	fAnalysisTree->Branch("LumiSec",          &fTLumiSection,       "LumiSec/I");

	// HLT triggers
	AddTriggerBranches();
	
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
	fAnalysisTree->Branch("MuCorrIso"     ,&fTmuciso,         "MuCorrIso[NMus]/F");
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
	fAnalysisTree->Branch("ElCorrIso",              &fTElCorrIso,           "ElCorrIso[NEls]/F");
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

void SSDLAnalysis::BookEffTree(){
	fOutputFile->cd();
	fLepEffTree = new TTree("LeptonEfficiency", "LeptonEfficiencyTree");

    // run/sample properties
	fLepEffTree->Branch("Run",              &fLETrun      ,       "Run/I");
	fLepEffTree->Branch("Event",            &fLETevent    ,       "Event/I");
	fLepEffTree->Branch("LumiSec",          &fLETlumi     ,       "LumiSec/I");
	fLepEffTree->Branch("Rho",              &fLETrho      ,       "Rho/F");
	fLepEffTree->Branch("NVrtx",            &fLETnvrtx    ,       "NVrtx/I");
	fLepEffTree->Branch("PUWeight",         &fLETpuweight ,       "PUWeight/F");
	fLepEffTree->Branch("Type",             &fLETtype     ,       "Type/I");
	fLepEffTree->Branch("Pt",               &fLETpt       ,       "Pt/F");
	fLepEffTree->Branch("Eta",              &fLETeta      ,       "Eta/F");
	fLepEffTree->Branch("Phi",              &fLETphi      ,       "Phi/F");
	fLepEffTree->Branch("Iso",              &fLETiso      ,       "Iso/F");
	fLepEffTree->Branch("Pass1",            &fLETpassed1  ,       "Pass1/I");
	fLepEffTree->Branch("Pass2",            &fLETpassed2  ,       "Pass2/I");
	fLepEffTree->Branch("Pass3",            &fLETpassed3  ,       "Pass3/I");
	fLepEffTree->Branch("Pass4",            &fLETpassed4  ,       "Pass4/I");
}

//____________________________________________________________________________
void SSDLAnalysis::Analyze(){
	FillAnalysisTree();
	if(fDoFillEffTree && !fIsData) FillEffTree();
}
void SSDLAnalysis::FillAnalysisTree(){
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
	vector<int> selectedMuInd  = MuonSelection(           &UserAnalysisBase::IsVeryLooseMu);
	vector<int> selectedElInd  = ElectronSelection(       &UserAnalysisBase::IsVeryLooseEl);
	vector<int> selectedJetInd = PFJetSelection(40., 2.5, &UserAnalysisBase::IsGoodBasicPFJet);
	fTnqmus  = std::min( (int)selectedMuInd .size(), fMaxNmus );
	fTnqels  = std::min( (int)selectedElInd .size(), fMaxNeles);
	fTnqjets = std::min( (int)selectedJetInd.size(), fMaxNjets);

	// Require at least one loose lepton
	if( (fTnqmus + fTnqels) < 1 ) return;
	fCounter.fill(fCutnames[3]);

	// Event and run info
	fTRunNumber   = fTR->Run;
	fTEventNumber = fTR->Event;
	fTLumiSection = fTR->LumiSection;

	// Dump basic jet and MET properties
	for(int ind = 0; ind < fTnqjets; ind++){
		int jetindex = selectedJetInd[ind];
		fTJetpt  [ind] = fTR->JPt [jetindex];
		fTJeteta [ind] = fTR->JEta[jetindex];
		fTJetphi [ind] = fTR->JPhi[jetindex];
		fTJetbtag[ind] = fTR->JbTagProbTkCntHighPur[jetindex];
		fTJetArea[ind] = fTR->JArea[jetindex];
	}

	// Get METs
	fTtcMET     = fTR->TCMET;
	fTtcMETphi  = fTR->TCMETphi;
	fTpfMET     = fTR->PFMET;
	fTpfMETphi  = fTR->PFMETphi;

	// PU correction
	fTrho   = fTR->Rho;
	fTnvrtx = fTR->NVrtx;
	if(!fIsData) fTpuweight = GetPUWeight(fTR->PUnumInteractions);
	else fTpuweight = 1.;

	// Dump muon properties
	for(int i = 0; i < fTnqmus; ++i){
		int index = selectedMuInd[i];
		fTmupt    [i] = fTR->MuPt      [index];
		fTmueta   [i] = fTR->MuEta     [index];
		fTmuphi   [i] = fTR->MuPhi     [index];
		fTmucharge[i] = fTR->MuCharge  [index];
		fTmuiso   [i] = fTR->MuRelIso03[index];
		fTmuciso  [i] = corrMuIso(index);
		fTmud0    [i] = fTR->MuD0PV    [index];
		fTmudz    [i] = fTR->MuDzPV    [index];
		fTmuptE   [i] = fTR->MuPtE     [index];
		
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

	// Dump electron properties
	for(int ind = 0; ind < fTnqels; ind++){
		int elindex = selectedElInd[ind];
		fTElcharge           [ind] = fTR->ElCharge                [elindex];
		fTElChargeIsCons     [ind] = fTR->ElCInfoIsGsfCtfScPixCons[elindex];
		fTElpt               [ind] = fTR->ElPt                    [elindex];
		fTEleta              [ind] = fTR->ElEta                   [elindex];
		fTElphi              [ind] = fTR->ElPhi                   [elindex];
		fTEld0               [ind] = fTR->ElD0PV                  [elindex];
		fTElD0Err            [ind] = fTR->ElD0E                   [elindex];
		fTEldz               [ind] = fTR->ElDzPV                  [elindex];
		fTElDzErr            [ind] = fTR->ElDzE                   [elindex];
		fTElRelIso           [ind] = relElIso(elindex); // correct by 1 GeV in ecal for barrel
		fTElCorrIso          [ind] = corrElIso(elindex);
		fTElEcalRecHitSumEt  [ind] = fTR->ElDR03EcalRecHitSumEt   [elindex];
		
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
		}
		else{
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
void SSDLAnalysis::FillEffTree(){
	if(fIsData){
		if(fVerbose > 0) cout << "Trying to fill MC lepton efficiencies on data, Stupid..." << endl;
		return;
	}
	ResetEffTree();

	fLETevent    = fTR->Run;
	fLETrun      = fTR->Event;
	fLETlumi     = fTR->LumiSection;
	fLETrho      = fTR->Rho;
	fLETnvrtx    = fTR->NVrtx;
	fLETpuweight = GetPUWeight(fTR->PUnumInteractions);

	// Muon loop
	for(int i = 0; i < fTR->NMus; ++i){
		if(IsSignalMuon(i) == false) continue; // match to signal mu
		if(fTR->MuPt[i] < 5.) continue; // pt cut
		fLETtype   = 0; // mu
		fLETpt     = fTR->MuPt      [i];
		fLETeta    = fTR->MuEta     [i];
		fLETphi    = fTR->MuPhi     [i];
		fLETiso    = fTR->MuRelIso03[i];
		
		fLETpassed1 = IsTightMuon(1,i)?1:0;
		fLETpassed2 = IsTightMuon(2,i)?1:0;
		fLETpassed3 = IsTightMuon(3,i)?1:0;
		fLETpassed4 = IsTightMuon(4,i)?1:0;
		fLepEffTree->Fill();
	}

	// Electron loop
	for(int i = 0; i < fTR->NEles; ++i){
		if(IsSignalElectron(i) == false) continue; // match to signal el
		if(fTR->ElPt[i] < 5.) continue; // pt cut
		fLETtype   = 1; // mu
		fLETpt     = fTR->ElPt      [i];
		fLETeta    = fTR->ElEta     [i];
		fLETphi    = fTR->ElPhi     [i];
		fLETiso    = relElIso(i);
		
		fLETpassed1 = IsTightEle(1,i)?1:0;
		fLETpassed2 = IsTightEle(2,i)?1:0;
		fLETpassed3 = IsTightEle(3,i)?1:0;
		fLETpassed4 = IsTightEle(4,i)?1:0;
		fLepEffTree->Fill();
	}
}

//____________________________________________________________________________
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
		fTmuciso        [i] = -999.99;
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
		fTElCorrIso         [i] = -999.99;
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
void SSDLAnalysis::ResetEffTree(){
	fLETevent      = -999;
	fLETrun        = -999;
	fLETlumi       = -999;
	fLETrho        = -999.99;
	fLETnvrtx      = -999;
	fLETpuweight   = -999.99;
	fLETtype       = -999;
	fLETpt         = -999.99;
	fLETeta        = -999.99;
	fLETphi        = -999.99;
	fLETiso        = -999.99;
	fLETpassed1    = -999;
	fLETpassed2    = -999;
	fLETpassed3    = -999;
	fLETpassed4    = -999;
}

//____________________________________________________________________________
// Some same-sign specific object selections
bool SSDLAnalysis::IsSignalMuon(int index){
	int matched = -1;
	// Match to a gen muon
	float mindr(100.);
	for(size_t i = 0; i < fTR->NGenLeptons; ++i){
		if(abs(fTR->GenLeptonID[i]) != 13) continue; // muons
		float eta = fTR->GenLeptonEta[i];
		float phi = fTR->GenLeptonPhi[i];
		float DR = Util::GetDeltaR(eta, fTR->MuEta[index], phi, fTR->MuPhi[index]);

		if(DR > mindr) continue; // minimize DR
		mindr = DR;
		matched = i;
	}

	if(mindr > 0.2) return false; // max distance
	if(matched < 0) return false; // match unsuccessful

	pdgparticle mo, gmo;
	GetPDGParticle(mo,  abs(fTR->GenLeptonMID [matched]));
	GetPDGParticle(gmo, abs(fTR->GenLeptonGMID[matched]));
	
	if(mo.get_type() > 10  || gmo.get_type() > 10)  return false; // mother is SM hadron
	if(mo.get_type() == 9  || gmo.get_type() == 9)  return true; // matched to susy
	if(mo.get_type() == 4 && abs(mo.get_pdgid()) != 21 && abs(mo.get_pdgid()) != 22) return true; // mother is gauge or higgs, but not gamma or gluon
	return false; // default
}
bool SSDLAnalysis::IsSignalElectron(int index){
	int matched = -1;
	// Match to a gen electron
	float mindr(100.);
	for(size_t i = 0; i < fTR->NGenLeptons; ++i){
		if(abs(fTR->GenLeptonID[i]) != 11) continue; // electrons
		float eta = fTR->GenLeptonEta[i];
		float phi = fTR->GenLeptonPhi[i];
		float DR = Util::GetDeltaR(eta, fTR->ElEta[index], phi, fTR->ElPhi[index]);

		if(DR > mindr) continue; // minimize DR
		mindr = DR;
		matched = i;
	}

	if(mindr > 0.2) return false; // max distance
	if(matched < 0) return false; // match unsuccessful

	pdgparticle mo, gmo;
	GetPDGParticle(mo,  abs(fTR->GenLeptonMID [matched]));
	GetPDGParticle(gmo, abs(fTR->GenLeptonGMID[matched]));
	
	if(mo.get_type() > 10  || gmo.get_type() > 10)  return false; // mother is SM hadron
	if(mo.get_type() == 9  || gmo.get_type() == 9)  return true; // matched to susy
	if(mo.get_type() == 4 && abs(mo.get_pdgid()) != 21 && abs(mo.get_pdgid()) != 22) return true; // mother is gauge or higgs, but not gamma or gluon
	return false; // default
}
bool SSDLAnalysis::IsTightMuon(int toggle, int index){
	if(toggle == 1){
		if(IsGoodBasicMu(index) == false) return false;
		if(fTR->MuPtE[index]/fTR->MuPt[index] > 0.1) return false;
		if(fTR->MuRelIso03[index] > 0.15) return false;
		return true;
	}
	if(toggle == 2){
		if(IsGoodBasicMu(index) == false) return false;
		if(fTR->MuPtE[index]/fTR->MuPt[index] > 0.1) return false;
		if(corrMuIso(index) > 0.15) return false;
		return true;
	}
	if(toggle == 3){
		if(IsGoodBasicMu(index) == false) return false;
		if(fTR->MuPtE[index]/fTR->MuPt[index] > 0.1) return false;
		if(fTR->MuRelIso03[index] > 0.1) return false;
		return true;
	}
	if(toggle == 4){
		if(IsGoodBasicMu(index) == false) return false;
		if(fTR->MuPtE[index]/fTR->MuPt[index] > 0.1) return false;
		if(corrMuIso(index) > 0.1) return false;
		return true;
	}
	cout << "Choose your toggle!" << endl;
	return false;
}
bool SSDLAnalysis::IsTightEle(int toggle, int index){
	if(toggle == 1){
		if(IsLooseEl(index) == false) return false;
		if(fTR->ElDR03EcalRecHitSumEt[index]/fTR->ElPt[index] > 0.2) return false;
		if(fTR->ElCInfoIsGsfCtfScPixCons[index] != 1) return false;
		if(IsGoodElId_WP80(index) == false) return false;
		if(relElIso(index) > 0.15) return false;
		return true;		
	}
	if(toggle == 2){
		if(IsLooseEl(index) == false) return false;
		if(fTR->ElDR03EcalRecHitSumEt[index]/fTR->ElPt[index] > 0.2) return false;
		if(fTR->ElCInfoIsGsfCtfScPixCons[index] != 1) return false;
		if(IsGoodElId_WP80(index) == false) return false;
		if(corrElIso(index) > 0.15) return false;
		return true;
	}
	if(toggle == 3){
		if(IsLooseEl(index) == false) return false;
		if(fTR->ElDR03EcalRecHitSumEt[index]/fTR->ElPt[index] > 0.2) return false;
		if(fTR->ElCInfoIsGsfCtfScPixCons[index] != 1) return false;
		if(IsGoodElId_WP80(index) == false) return false;
		if(relElIso(index) > 0.1) return false;
		return true;		
	}
	if(toggle == 4){
		if(IsLooseEl(index) == false) return false;
		if(fTR->ElDR03EcalRecHitSumEt[index]/fTR->ElPt[index] > 0.2) return false;
		if(fTR->ElCInfoIsGsfCtfScPixCons[index] != 1) return false;
		if(IsGoodElId_WP80(index) == false) return false;
		if(corrElIso(index) > 0.1) return false;
		return true;
	}
	cout << "Choose your toggle!" << endl;
	return false;
}

//____________________________________________________________________________
// Same-sign specific isolation
double SSDLAnalysis::corrMuIso(int index){
	double newiso = fTR->MuRelIso03[index] - TMath::Log(fTR->MuPt[index]) * (fTR->NVrtx - 1) / (30. * fTR->MuPt[index]);
	if(newiso < 0.) newiso = 0.; // cut off at 0
	return newiso;
}
double SSDLAnalysis::corrElIso(int index){
	double newiso = relElIso(index) - TMath::Log(fTR->ElPt[index]) * (fTR->NVrtx - 1) / (30. * fTR->ElPt[index]);
	if(newiso < 0.) newiso = 0.; // cut off at 0
	return newiso;
}
