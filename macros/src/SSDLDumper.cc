/*****************************************************************************
*   Collection of tools for producing plots for same-sign dilepton analysis  *
*                                                                            *
*                                                  (c) 2010 Benjamin Stieger *
******************************************************************************
* This runs on SSDLTrees, applies all object and event selections for signal *
* and control regions of the SSDL analysis, and stores the relevant numbers  *
* in a set of histograms in a ROOT file.                                     *
*****************************************************************************/
#include "SSDLDumper.hh"

#include "helper/AnaClass.hh"
#include "helper/Utilities.hh"
#include "helper/FPRatios.hh"
#include "helper/Davismt2.h"
#include "helper/FakeRatios.hh"
#include "helper/Monitor.hh"

#include "TLorentzVector.h"
#include "TGraphAsymmErrors.h"
#include "TEfficiency.h"

#include "TWbox.h"
#include "TMultiGraph.h"
#include "TGaxis.h"

#include <iostream>
#include <iomanip>
#include <time.h> // access to date/time

using namespace std;

//////////////////////////////////////////////////////////////////////////////////
// Global parameters:

static const float gMMU = 0.1057;
static const float gMEL = 0.0005;
static const float gMZ  = 91.2;

static const double gStatBetaAlpha = 1.;
static const double gStatBetaBeta  = 1.;

// Regions ///////////////////////////////////////////////////////////////////////
float SSDLDumper::Region::minMu1pt[2] = {20.,   5.};
float SSDLDumper::Region::minMu2pt[2] = {10.,   5.};
float SSDLDumper::Region::minEl1pt[2] = {20.,  10.};
float SSDLDumper::Region::minEl2pt[2] = {10.,  10.};

TString SSDLDumper::Region::sname [SSDLDumper::gNREGIONS] = {"HT80MET30", "HT80MET20to50", "HT80MET100", "HT200MET30", "HT200MET120", "HT400MET50", "HT400MET120", "HT400MET0",  "HT80METInv30", "HT0MET120"};
float SSDLDumper::Region::minHT   [SSDLDumper::gNREGIONS] = {        80.,             80.,          80.,         200.,          200.,         400.,          400.,         400.,            80.,          0.};
float SSDLDumper::Region::maxHT   [SSDLDumper::gNREGIONS] = {      7000.,           7000.,        7000.,        7000.,         7000.,        7000.,         7000.,        7000.,          7000.,       7000.};
float SSDLDumper::Region::minMet  [SSDLDumper::gNREGIONS] = {        30.,             20.,         100.,          30.,          120.,          50.,          120.,           0.,             0.,        120.};
float SSDLDumper::Region::maxMet  [SSDLDumper::gNREGIONS] = {      7000.,             50.,        7000.,        7000.,         7000.,        7000.,         7000.,        7000.,            30.,       7000.};
int   SSDLDumper::Region::minNjets[SSDLDumper::gNREGIONS] = {         2 ,              2 ,           2 ,           2 ,            2 ,           2 ,            2 ,           2 ,             2 ,          0 };

// Muon Binning //////////////////////////////////////////////////////////////////
double SSDLDumper::gMuPtbins [gNMuPtbins+1]  = {5., 10., 15., 20., 25., 35., 45, 60.};
double SSDLDumper::gMuPt2bins[gNMuPt2bins+1] = {5., 10., 15., 20., 25., 35., 45, 60.};
double SSDLDumper::gMuEtabins[gNMuEtabins+1] = {0., 1.0, 1.479, 2.0, 2.5};

// Electron Binning //////////////////////////////////////////////////////////////
double SSDLDumper::gElPtbins [gNElPtbins+1]  = {10., 15., 20., 25., 35., 55.};
double SSDLDumper::gElPt2bins[gNElPt2bins+1] = {10., 15., 20., 25., 35., 55.};
double SSDLDumper::gElEtabins[gNElEtabins+1] = {0., 1.0, 1.479, 2.0, 2.5};
//////////////////////////////////////////////////////////////////////////////////

// NVrtx Binning //////////////////////////////////////////////////////////////
double SSDLDumper::gNVrtxBins[gNNVrtxBins+1]  = {0, 2, 4, 6, 8, 10, 12, 14, 16, 18};
//////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////
TString SSDLDumper::gKinSelNames[gNKinSels] = {"LL", "TT", "Sig"};
TString SSDLDumper::KinPlots::var_name[SSDLDumper::gNKinVars] = {"HT", "MET", "NJets", "Pt1", "Pt2", "InvMassSF", "InvMassMM", "InvMassEE", "InvMassEM", "MT2"};
int     SSDLDumper::KinPlots::nbins[SSDLDumper::gNKinVars]    = { 20 ,   20 ,      8 ,   20 ,   20 ,        30  ,        30  ,        30  ,        30  ,   20 };
float   SSDLDumper::KinPlots::xmin[SSDLDumper::gNKinVars]     = {100.,    0.,      0.,   10.,   10.,        20. ,        20. ,        20. ,        20. ,    0.};
float   SSDLDumper::KinPlots::xmax[SSDLDumper::gNKinVars]     = {800.,  210.,      8.,  200.,  100.,       300. ,       300. ,       300. ,       300. ,  100.};
TString SSDLDumper::KinPlots::axis_label[SSDLDumper::gNKinVars] = {"H_{T} (GeV)",
                                                                     "E_{T}^{miss} (GeV)",
                                                                     "N_{Jets}",
                                                                     "P_{T} (l_{1}) (GeV)",
                                                                     "P_{T} (l_{2}) (GeV)",
                                                                     "m_{ll} (SF) (GeV)",
                                                                     "m_{#mu#mu} (GeV)",
                                                                     "m_{ee} (GeV)",
                                                                     "m_{ll} (OF) (GeV)",
                                                                     "M_{T2} (GeV)"};

//////////////////////////////////////////////////////////////////////////////////
double SSDLDumper::gDiffHTBins[gNDiffHTBins+1]   = { 0., 100., 200.,  300., 400., 900.};
double SSDLDumper::gDiffMETBins[gNDiffMETBins+1] = {30.,  60.,  90., 200.};
double SSDLDumper::gDiffNJBins[gNDiffNJBins+1]   = { 0.,   1.,   2.,    3.,   4.,   8.}; // fill NJets + 0.5 to hit the right bin
double SSDLDumper::gDiffMT2Bins[gNDiffMT2Bins+1] = { 0.,  25.,  50., 150.};
double SSDLDumper::gDiffPT1Bins[gNDiffPT1Bins+1] = { 20., 40., 60., 80., 100., 200.};
double SSDLDumper::gDiffPT2Bins[gNDiffPT2Bins+1] = { 10., 20., 30., 40.,  50., 100.};

//////////////////////////////////////////////////////////////////////////////////
TString SSDLDumper::DiffPredYields::var_name[SSDLDumper::gNDiffVars] = {"HT", "MET", "NJets", "MT2", "PT1", "PT2"};
int     SSDLDumper::DiffPredYields::nbins[SSDLDumper::gNDiffVars]    = {gNDiffHTBins, gNDiffMETBins, gNDiffNJBins, gNDiffMT2Bins, gNDiffPT1Bins, gNDiffPT2Bins};
double* SSDLDumper::DiffPredYields::bins[SSDLDumper::gNDiffVars]     = {gDiffHTBins,  gDiffMETBins,  gDiffNJBins,  gDiffMT2Bins,  gDiffPT1Bins,  gDiffPT2Bins};
TString SSDLDumper::DiffPredYields::axis_label[SSDLDumper::gNDiffVars] = {"H_{T} (GeV)",
                                                                           "E_{T}^{miss} (GeV)",
                                                                           "N_{Jets}",
                                                                           "M_{T2} (GeV)",
                                                                           "P_{T} (l_{1}) (GeV)",
                                                                           "P_{T} (l_{2}) (GeV)"};

//////////////////////////////////////////////////////////////////////////////////
TString SSDLDumper::FRatioPlots::var_name[SSDLDumper::gNRatioVars] = {"NJets",  "HT", "MaxJPt", "NVertices", "ClosJetPt", "AwayJetPt", "NBJets", "MET",  "MT"};
int     SSDLDumper::FRatioPlots::nbins[SSDLDumper::gNRatioVars]    = {     7 ,   20 ,      20 ,         9  ,        20  ,        20  ,       3 ,   10 ,   10 };
float   SSDLDumper::FRatioPlots::xmin[SSDLDumper::gNRatioVars]     = {     1.,   50.,      30.,         0. ,        30. ,        30. ,       0.,    0.,    0.};
float   SSDLDumper::FRatioPlots::xmax[SSDLDumper::gNRatioVars]     = {     8.,  500.,     300.,        18. ,       150. ,       300. ,       3.,   50.,  100.};

//////////////////////////////////////////////////////////////////////////////////
TString SSDLDumper::IsoPlots::sel_name[SSDLDumper::gNSels] = {"Base", "SigSup"};
int     SSDLDumper::IsoPlots::nbins[SSDLDumper::gNSels]    = {20, 20};


TString SSDLDumper::gEMULabel[2] = {"mu", "el"};
TString SSDLDumper::gChanLabel[3] = {"MM", "EM", "EE"}; // make SURE this is the same order as gChannel enum!
TString SSDLDumper::gHiLoLabel[3] = {"HighPt", "LowPt", "TauChan"};

//____________________________________________________________________________
SSDLDumper::SSDLDumper(){}

//____________________________________________________________________________
SSDLDumper::~SSDLDumper(){
	if(fOutputFile != NULL && fOutputFile->IsOpen()) fOutputFile->Close();
	fChain = 0;
}

//____________________________________________________________________________
void SSDLDumper::init(TString inputfile, TString sname, int datamc, int chan){
	fSample = new Sample(inputfile, sname, datamc, chan);
	fSamples.push_back(fSample);

	if(fVerbose > 0) cout << "------------------------------------" << endl;
	if(fVerbose > 0) cout << " Initializing SSDLDumper ... " << endl;
	if(fVerbose > 0) cout << "   Running on:      " << fSample->location << endl;
	if(fVerbose > 0) cout << "   Naming it:       " << fSample->sname << endl;
	if(fVerbose > 0) cout << "   Is data/mc:      " << fSample->datamc << endl;
	if(fVerbose > 0) cout << "------------------------------------" << endl;
	init();
}
void SSDLDumper::init(TString datacard){
	readDatacard(datacard);
	if(fVerbose > 0) cout << "------------------------------------" << endl;
	if(fVerbose > 0) cout << " Initializing SSDLDumper ... " << endl;
	if(fVerbose > 0) cout << "   Running on datacard " << datacard << endl;
	if(fVerbose > 0) cout << "------------------------------------" << endl;
	init();
}

void SSDLDumper::init(){
	Util::SetStyle();
	fDoCounting = false;  // Disable counters by default

	resetHypLeptons();
	initCutNames();
	
	// Cuts:
	fC_minMu1pt = 20.;
	fC_minMu2pt = 10.;
	fC_minEl1pt = 20.;
	fC_minEl2pt = 10.;
	fC_minHT    = 80.;
	fC_minMet   = 30.;
	fC_maxHT    = 7000.;
	fC_maxMet   = 7000.;
	fC_minNjets = 2;
	
	fC_maxMet_Control = 20.;
	fC_maxMt_Control  = 20.;

	// Prevent root from adding histograms to current file
	TH1::AddDirectory(kFALSE);
}
void SSDLDumper::InitMC(TTree *tree){
// Copied from MetaTreeClassBase, remove a few branches that are not in older version of minitrees
	
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("Run", &Run, &b_Run);
   fChain->SetBranchAddress("Event", &Event, &b_Event);
   fChain->SetBranchAddress("LumiSec", &LumiSec, &b_LumiSec);
   fChain->SetBranchAddress("Rho", &Rho, &b_Rho);
   fChain->SetBranchAddress("NVrtx", &NVrtx, &b_NVrtx);
   fChain->SetBranchAddress("PUWeight", &PUWeight, &b_PUWeight);
   fChain->SetBranchAddress("NMus", &NMus, &b_NMus);
   fChain->SetBranchAddress("MuPt", MuPt, &b_MuPt);
   fChain->SetBranchAddress("MuEta", MuEta, &b_MuEta);
   fChain->SetBranchAddress("MuPhi", MuPhi, &b_MuPhi);
   fChain->SetBranchAddress("MuCharge", MuCharge, &b_MuCharge);
   fChain->SetBranchAddress("MuIso", MuIso, &b_MuIso);
   fChain->SetBranchAddress("MuD0", MuD0, &b_MuD0);
   fChain->SetBranchAddress("MuDz", MuDz, &b_MuDz);
   fChain->SetBranchAddress("MuPtE", MuPtE, &b_MuPtE);
   fChain->SetBranchAddress("MuGenID", MuGenID, &b_MuGenID);
   fChain->SetBranchAddress("MuGenMID", MuGenMID, &b_MuGenMID);
   fChain->SetBranchAddress("MuGenGMID", MuGenGMID, &b_MuGenGMID);
   fChain->SetBranchAddress("MuGenType", MuGenType, &b_MuGenType);
   fChain->SetBranchAddress("MuGenMType", MuGenMType, &b_MuGenMType);
   fChain->SetBranchAddress("MuGenGMType", MuGenGMType, &b_MuGenGMType);
   fChain->SetBranchAddress("MuMT", MuMT, &b_MuMT);
   fChain->SetBranchAddress("NEls", &NEls, &b_NEls);
   fChain->SetBranchAddress("ElCharge", ElCharge, &b_ElCharge);
   fChain->SetBranchAddress("ElChIsCons", ElChIsCons, &b_ElChIsCons);
   fChain->SetBranchAddress("ElPt", ElPt, &b_ElPt);
   fChain->SetBranchAddress("ElEta", ElEta, &b_ElEta);
   fChain->SetBranchAddress("ElPhi", ElPhi, &b_ElPhi);
   fChain->SetBranchAddress("ElD0", ElD0, &b_ElD0);
   fChain->SetBranchAddress("ElD0Err", ElD0Err, &b_ElD0Err);
   fChain->SetBranchAddress("ElDz", ElDz, &b_ElDz);
   fChain->SetBranchAddress("ElDzErr", ElDzErr, &b_ElDzErr);
   fChain->SetBranchAddress("ElRelIso", ElRelIso, &b_ElRelIso);
   fChain->SetBranchAddress("ElEcalRecHitSumEt", ElEcalRecHitSumEt, &b_ElEcalRecHitSumEt);
   fChain->SetBranchAddress("ElIsGoodElId_WP80", ElIsGoodElId_WP80, &b_ElIsGoodElId_WP80);
   fChain->SetBranchAddress("ElIsGoodElId_WP90", ElIsGoodElId_WP90, &b_ElIsGoodElId_WP90);
   fChain->SetBranchAddress("ElGenID", ElGenID, &b_ElGenID);
   fChain->SetBranchAddress("ElGenMID", ElGenMID, &b_ElGenMID);
   fChain->SetBranchAddress("ElGenGMID", ElGenGMID, &b_ElGenGMID);
   fChain->SetBranchAddress("ElGenType", ElGenType, &b_ElGenType);
   fChain->SetBranchAddress("ElGenMType", ElGenMType, &b_ElGenMType);
   fChain->SetBranchAddress("ElGenGMType", ElGenGMType, &b_ElGenGMType);
   fChain->SetBranchAddress("ElMT", ElMT, &b_ElMT);
   fChain->SetBranchAddress("tcMET", &tcMET, &b_tcMET);
   fChain->SetBranchAddress("tcMETPhi", &tcMETPhi, &b_tcMETPhi);
   fChain->SetBranchAddress("pfMET", &pfMET, &b_pfMET);
   fChain->SetBranchAddress("pfMETPhi", &pfMETPhi, &b_pfMETPhi);
   fChain->SetBranchAddress("NJets", &NJets, &b_NJets);
   fChain->SetBranchAddress("JetPt", JetPt, &b_JetPt);
   fChain->SetBranchAddress("JetEta", JetEta, &b_JetEta);
   fChain->SetBranchAddress("JetPhi", JetPhi, &b_JetPhi);
   fChain->SetBranchAddress("JetSSVHPBTag", JetSSVHPBTag, &b_JetSSVHPBTag);
   fChain->SetBranchAddress("JetArea", JetArea, &b_JetArea);
   Notify();
}

void SSDLDumper::readDatacard(TString cardfile){
	char buffer[1000];
	ifstream IN(cardfile);

	char inputfile[100], sname[100];
	int datamc, chan, color;
	float lumi;

	if(fVerbose > 2) cout << "------------------------------------" << endl;
	if(fVerbose > 2) cout << "Reading datacard  " << cardfile << endl;
	int counter(0);

	while( IN.getline(buffer, 200, '\n') ){
		lumi = 1.0;
		color = 1;
		if (buffer[0] == '#') continue; // Skip lines commented with '#'
		if( sscanf(buffer, "%s\t%s\t%d\t%d\t%f\t%d", sname, inputfile, &datamc, &chan, &lumi, &color) > 3){
			Sample *S = new Sample(inputfile, sname, datamc, chan, lumi, color);
			fSamples.push_back(S);
			if(fVerbose > 2){
				cout << " ---- " << endl;
				cout << "  New sample added: " << S->sname << endl;
				cout << "   Input:      " << S->location << endl;
				cout << "   DataMC:     " << S->datamc << endl;
				cout << "   Channel:    " << S->chansel << endl;
				cout << "   Lumi:       " << S->lumi << endl;
				cout << "   Color:      " << S->color << endl;
			}
		}
		else{
			cout << " SSDLDumper::readDatacard ==> Wrong dataformat in datacard " << cardfile << " detected, aborting..." << endl;
			exit(1);
		}
	}
	if(fVerbose > 2) cout << "------------------------------------" << endl;
	IN.close();
}

//____________________________________________________________________________
const int     SSDLDumper::getNPtBins (gChannel chan){
	if(chan == Muon || chan == ElMu) return gNMuPtbins;
	if(chan == Elec) return gNElPtbins;
}
const double *SSDLDumper::getPtBins  (gChannel chan){
	if(chan == Muon || chan == ElMu) return gMuPtbins;
	if(chan == Elec) return gElPtbins;
}
const int     SSDLDumper::getNPt2Bins(gChannel chan){
	if(chan == Muon || chan == ElMu) return gNMuPt2bins;
	if(chan == Elec) return gNElPt2bins;
}
const double *SSDLDumper::getPt2Bins (gChannel chan){
	if(chan == Muon || chan == ElMu) return gMuPt2bins;
	if(chan == Elec) return gElPt2bins;
}
const int     SSDLDumper::getNEtaBins(gChannel chan){
	if(chan == Muon || chan == ElMu) return gNMuEtabins;
	if(chan == Elec)            return gNElEtabins;
}
const double *SSDLDumper::getEtaBins (gChannel chan){
	if(chan == Muon || chan == ElMu) return gMuEtabins;
	if(chan == Elec)            return gElEtabins;
}

//____________________________________________________________________________
void SSDLDumper::loop(){
	for(size_t i = 0; i < fSamples.size(); ++i){
		fSample = fSamples[i]; // TODO: Clean this up, call the triggers with an argument
		fOutputFileName = fOutputDir + fSample->sname + "_Yields.root";
		loopEvents(fSample);
	}
}
void SSDLDumper::loopEvents(Sample *S){
	fDoCounting = true;
	if(S->datamc == 0){
		TString eventfilename  = fOutputDir + S->sname + "_SignalEvents.txt";
		fOUTSTREAM.open(eventfilename.Data(), ios::trunc);		
	}

	TFile *pFile = new TFile(fOutputFileName, "RECREATE");
		
	bookHistos(S);
	
	TTree *tree = S->getTree();

	// Stuff to execute for each sample BEFORE looping on the events
	initCounters();

	// Event loop
	tree->ResetBranchAddresses();
	// Init(tree);
	if(S->datamc == 0) Init(tree);
	else InitMC(tree);

	if (fChain == 0) return;
	Long64_t nentries = fChain->GetEntriesFast();
	Long64_t nbytes = 0, nb = 0;
	for (Long64_t jentry=0; jentry<nentries;jentry++) {
		printProgress(jentry, nentries, S->sname);

		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		nb = fChain->GetEntry(jentry);   nbytes += nb;

		fCounter[Muon].fill(fMMCutNames[0]);
		fCounter[ElMu].fill(fEMCutNames[0]);
		fCounter[Elec].fill(fEECutNames[0]);

		// Select mutually exclusive runs for Jet and MultiJet datasets
		if(!isGoodRun(S)) continue;

		fCounter[Muon].fill(fMMCutNames[1]);
		fCounter[ElMu].fill(fEMCutNames[1]);
		fCounter[Elec].fill(fEECutNames[1]);

		fillKinPlots(S, HighPt);
		fillKinPlots(S, LowPt);
		for(gRegion r = region_begin; r < gNREGIONS; r=gRegion(r+1)) fillYields(S, r, HighPt);
		for(gRegion r = region_begin; r < gNREGIONS; r=gRegion(r+1)) fillYields(S, r, LowPt);

		fillDiffYields(S);
		fillRatioPlots(S);
		fillMuIsoPlots(S);
		fillElIsoPlots(S);

	}

	// Stuff to execute for each sample AFTER looping on the events
	fillCutFlowHistos(S);

	writeHistos(S, pFile);
	writeSigGraphs(S, Muon, HighPt, pFile);
	writeSigGraphs(S, Muon, LowPt,  pFile);
	writeSigGraphs(S, Elec, HighPt, pFile);
	writeSigGraphs(S, Elec, LowPt,  pFile);
	writeSigGraphs(S, ElMu, HighPt, pFile);
	writeSigGraphs(S, ElMu, LowPt,  pFile);

	deleteHistos(S);
	S->cleanUp();

	pFile->Write();
	pFile->Close();

	if(S->datamc == 0) fOUTSTREAM.close();
	fDoCounting = false;
}

//____________________________________________________________________________
void SSDLDumper::fillYields(Sample *S, gRegion reg, gHiLoSwitch hilo){
	///////////////////////////////////////////////////
	// Set custom event selections here:
	fC_minMu1pt = Region::minMu1pt[hilo];
	fC_minMu2pt = Region::minMu2pt[hilo];
	fC_minEl1pt = Region::minEl1pt[hilo];
	fC_minEl2pt = Region::minEl2pt[hilo];
	fC_minHT    = Region::minHT   [reg];
	fC_maxHT    = Region::maxHT   [reg];
	fC_minMet   = Region::minMet  [reg];
	fC_maxMet   = Region::maxMet  [reg];
	fC_minNjets = Region::minNjets[reg];

	///////////////////////////////////////////////////
	// SS YIELDS
	// MuMu Channel
	resetHypLeptons();
	fDoCounting = false;
	if(reg == Baseline && hilo == HighPt) fDoCounting = true;
	fCurrentChannel = Muon;
	float puweight = PUWeight;
	if(S->datamc == 4) puweight = 1; // fix for samples with no pileup
	int mu1(-1), mu2(-1);
	if(mumuSignalTrigger()){ // Trigger selection
		if(fDoCounting) fCounter[Muon].fill(fMMCutNames[2]);
		if(isSSLLMuEvent(mu1, mu2)){ // Same-sign loose-loose di muon event
			if(  isTightMuon(mu1) &&  isTightMuon(mu2) ){ // Tight-tight
				if(fDoCounting) fCounter[Muon].fill(fMMCutNames[14]); // ... first muon passes tight cut
				if(fDoCounting) fCounter[Muon].fill(fMMCutNames[15]); // ... second muon passes tight cut
				if(fDoCounting) fCounter[Muon].fill(fMMCutNames[16]); // ... both muons pass tight cut
				S->region[reg][hilo].mm.nt20_pt ->Fill(MuPt [mu1], MuPt [mu2], puweight);
				S->region[reg][hilo].mm.nt20_eta->Fill(fabs(MuEta[mu1]), fabs(MuEta[mu2]), puweight);
				if(S->datamc == 0 && reg == Baseline && hilo == HighPt){
					fOUTSTREAM << Form("%12s: MuMu - run %6.0d / ls %5.0d / ev %11.0d - HT(#J/#bJ) %6.2f(%1d/%1d) MET %6.2f MT2 %6.2f Pt1 %6.2f Pt2 %6.2f Charge %2d", S->sname.Data(), Run, LumiSec, Event, getHT(), getNJets(), getNBTags(), pfMET, getMT2(mu1,mu2,1), MuPt[mu1], MuPt[mu2], MuCharge[mu1]) << endl ;
				}
				if(S->datamc > 0 ){
					S->region[reg][hilo].mm.nt11_origin->Fill(muIndexToBin(mu1)-0.5, muIndexToBin(mu2)-0.5, puweight);
					if(isPromptMuon(mu1) && isPromptMuon(mu2)) S->region[reg][hilo].mm.nt2pp_pt->Fill(MuPt[mu1], MuPt[mu2], puweight);
					if(isPromptMuon(mu1) &&   isFakeMuon(mu2)) S->region[reg][hilo].mm.nt2pf_pt->Fill(MuPt[mu1], MuPt[mu2], puweight);
					if(  isFakeMuon(mu1) && isPromptMuon(mu2)) S->region[reg][hilo].mm.nt2fp_pt->Fill(MuPt[mu1], MuPt[mu2], puweight);
					if(  isFakeMuon(mu1) &&   isFakeMuon(mu2)) S->region[reg][hilo].mm.nt2ff_pt->Fill(MuPt[mu1], MuPt[mu2], puweight);
				}
			}
			if(  isTightMuon(mu1) && !isTightMuon(mu2) ){ // Tight-loose
				if(fDoCounting) fCounter[Muon].fill(fMMCutNames[14]); // ... first muon passes tight cut
				S->region[reg][hilo].mm.nt10_pt ->Fill(MuPt [mu1], MuPt [mu2], puweight);
				S->region[reg][hilo].mm.nt10_eta->Fill(fabs(MuEta[mu1]), fabs(MuEta[mu2]), puweight);
				if(S->datamc > 0) S->region[reg][hilo].mm.nt10_origin->Fill(muIndexToBin(mu1)-0.5, muIndexToBin(mu2)-0.5, puweight);
			}
			if( !isTightMuon(mu1) &&  isTightMuon(mu2) ){ // Loose-tight
				if(fDoCounting) fCounter[Muon].fill(fMMCutNames[15]); // ... second muon passes tight cut
				S->region[reg][hilo].mm.nt10_pt ->Fill(MuPt [mu2], MuPt [mu1], puweight); // tight one always in x axis; fill same again
				S->region[reg][hilo].mm.nt10_eta->Fill(fabs(MuEta[mu2]), fabs(MuEta[mu1]), puweight);
				if(S->datamc > 0) S->region[reg][hilo].mm.nt10_origin->Fill(muIndexToBin(mu2)-0.5, muIndexToBin(mu1)-0.5, puweight);
			}
			if( !isTightMuon(mu1) && !isTightMuon(mu2) ){ // Loose-loose
				S->region[reg][hilo].mm.nt00_pt ->Fill(MuPt [mu1], MuPt [mu2], puweight);
				S->region[reg][hilo].mm.nt00_eta->Fill(fabs(MuEta[mu1]), fabs(MuEta[mu2]), puweight);
				if(S->datamc > 0) S->region[reg][hilo].mm.nt00_origin->Fill(muIndexToBin(mu1)-0.5, muIndexToBin(mu2)-0.5, puweight);
			}
			if(S->datamc > 0){
				if(isPromptMuon(mu1) && isPromptMuon(mu2)) S->region[reg][hilo].mm.npp_pt->Fill(MuPt[mu1], MuPt[mu2], puweight);
				if(isPromptMuon(mu1) &&   isFakeMuon(mu2)) S->region[reg][hilo].mm.npf_pt->Fill(MuPt[mu1], MuPt[mu2], puweight);
				if(  isFakeMuon(mu1) && isPromptMuon(mu2)) S->region[reg][hilo].mm.nfp_pt->Fill(MuPt[mu1], MuPt[mu2], puweight);
				if(  isFakeMuon(mu1) &&   isFakeMuon(mu2)) S->region[reg][hilo].mm.nff_pt->Fill(MuPt[mu1], MuPt[mu2], puweight);			
			}
		}
		resetHypLeptons();
	}
	if(singleMuTrigger() && isSigSupMuEvent()){
		if( isTightMuon(0) ){
			S->region[reg][hilo].mm.fntight->Fill(MuPt[0], fabs(MuEta[0]), singleMuPrescale() * puweight);
			if(S->datamc > 0) S->region[reg][hilo].mm.sst_origin->Fill(muIndexToBin(0)-0.5, puweight);
		}
		if( isLooseMuon(0) ){
			S->region[reg][hilo].mm.fnloose->Fill(MuPt[0], fabs(MuEta[0]), singleMuPrescale() * puweight);
			if(S->datamc > 0) S->region[reg][hilo].mm.ssl_origin->Fill(muIndexToBin(0)-0.5, puweight);
		}
	}
	if(doubleMuTrigger() && isZMuMuEvent()){
		if( isTightMuon(0) ){
			S->region[reg][hilo].mm.pntight->Fill(MuPt[0], fabs(MuEta[0]), puweight);
			if(S->datamc > 0) S->region[reg][hilo].mm.zt_origin->Fill(muIndexToBin(0)-0.5, puweight);
		}
		if( isLooseMuon(0) ){
			S->region[reg][hilo].mm.pnloose->Fill(MuPt[0], fabs(MuEta[0]), puweight);
			if(S->datamc > 0) S->region[reg][hilo].mm.zl_origin->Fill(muIndexToBin(0)-0.5, puweight);
		}
	}				

	// EE Channel
	fCurrentChannel = Elec;
	int el1(-1), el2(-1);
	if(elelSignalTrigger()){
		if(fDoCounting) fCounter[Elec].fill(fEECutNames[2]);
		if( isSSLLElEvent(el1, el2) ){
			if(  isTightElectron(el1) &&  isTightElectron(el2) ){ // Tight-tight
				if(fDoCounting) fCounter[Elec].fill(fEECutNames[14]); // " ... first electron passes tight cut
				if(fDoCounting) fCounter[Elec].fill(fEECutNames[15]); // " ... second electron passes tight cut
				if(fDoCounting) fCounter[Elec].fill(fEECutNames[16]); // " ... both electrons pass tight cut
				S->region[reg][hilo].ee.nt20_pt ->Fill(ElPt [el1], ElPt [el2], puweight);
				S->region[reg][hilo].ee.nt20_eta->Fill(fabs(ElEta[el1]), fabs(ElEta[el2]), puweight);
				if(S->datamc == 0 && reg == Baseline && hilo == HighPt){
					fOUTSTREAM << Form("%12s: ElEl - run %6.0d / ls %5.0d / ev %11.0d - HT(#J/#bJ) %6.2f(%1d/%1d) MET %6.2f MT2 %6.2f Pt1 %6.2f Pt2 %6.2f Charge %2d", S->sname.Data(), Run, LumiSec, Event, getHT(), getNJets(), getNBTags(), pfMET, getMT2(el1,el2,2), ElPt[el1], ElPt[el2], ElCharge[el1]) << endl ;
				}
				if(S->datamc > 0 ){
					S->region[reg][hilo].ee.nt11_origin->Fill(elIndexToBin(el1)-0.5, elIndexToBin(el2)-0.5, puweight);
					if(isPromptElectron(el1) && isPromptElectron(el2)){
						S->region[reg][hilo].ee.nt2pp_pt->Fill(ElPt[el1], ElPt[el2], puweight);
						if(!isChargeMatchedElectron(el1) || !isChargeMatchedElectron(el2)){
							S->region[reg][hilo].ee.nt2pp_cm_pt->Fill(ElPt[el1], ElPt[el2], puweight);							
						}
					}
					if(isPromptElectron(el1) &&   isFakeElectron(el2)) S->region[reg][hilo].ee.nt2pf_pt->Fill(ElPt[el1], ElPt[el2], puweight);
					if(  isFakeElectron(el1) && isPromptElectron(el2)) S->region[reg][hilo].ee.nt2fp_pt->Fill(ElPt[el1], ElPt[el2], puweight);
					if(  isFakeElectron(el1) &&   isFakeElectron(el2)) S->region[reg][hilo].ee.nt2ff_pt->Fill(ElPt[el1], ElPt[el2], puweight);
				}
			}
			if(  isTightElectron(el1) && !isTightElectron(el2) ){ // Tight-loose
				if(fDoCounting) fCounter[Elec].fill(fEECutNames[14]);
				S->region[reg][hilo].ee.nt10_pt ->Fill(ElPt [el1], ElPt [el2], puweight);
				S->region[reg][hilo].ee.nt10_eta->Fill(fabs(ElEta[el1]), fabs(ElEta[el2]), puweight);
				if(S->datamc > 0) S->region[reg][hilo].ee.nt10_origin->Fill(elIndexToBin(el1)-0.5, elIndexToBin(el2)-0.5, puweight);
			}
			if( !isTightElectron(el1) &&  isTightElectron(el2) ){ // Loose-tight
				if(fDoCounting) fCounter[Elec].fill(fEECutNames[15]);
				S->region[reg][hilo].ee.nt10_pt ->Fill(ElPt [el2], ElPt [el1], puweight); // tight one always in x axis; fill same again
				S->region[reg][hilo].ee.nt10_eta->Fill(fabs(ElEta[el2]), fabs(ElEta[el2]), puweight);
				if(S->datamc > 0) S->region[reg][hilo].ee.nt10_origin->Fill(elIndexToBin(el2)-0.5, elIndexToBin(el1)-0.5, puweight);
			}
			if( !isTightElectron(el1) && !isTightElectron(el2) ){ // Loose-loose
				S->region[reg][hilo].ee.nt00_pt ->Fill(ElPt [el1], ElPt [el2], puweight);
				S->region[reg][hilo].ee.nt00_eta->Fill(fabs(ElEta[el1]), fabs(ElEta[el2]), puweight);
				if(S->datamc > 0) S->region[reg][hilo].ee.nt00_origin->Fill(elIndexToBin(el1)-0.5, elIndexToBin(el2)-0.5, puweight);
			}
			if(S->datamc > 0 ){
				if(isPromptElectron(el1) && isPromptElectron(el2)){
					S->region[reg][hilo].ee.npp_pt->Fill(ElPt[el1], ElPt[el2], puweight);
					if(!isChargeMatchedElectron(el1) || !isChargeMatchedElectron(el2)){
						S->region[reg][hilo].ee.npp_cm_pt->Fill(ElPt[el1], ElPt[el2], puweight);							
					}
				}
				if(isPromptElectron(el1) &&   isFakeElectron(el2)) S->region[reg][hilo].ee.npf_pt->Fill(ElPt[el1], ElPt[el2], puweight);
				if(  isFakeElectron(el1) && isPromptElectron(el2)) S->region[reg][hilo].ee.nfp_pt->Fill(ElPt[el1], ElPt[el2], puweight);
				if(  isFakeElectron(el1) &&   isFakeElectron(el2)) S->region[reg][hilo].ee.nff_pt->Fill(ElPt[el1], ElPt[el2], puweight);
			}
		}
		resetHypLeptons();
	}
	if(singleElTrigger() && isSigSupElEvent()){
		if( isTightElectron(0) ){
			S->region[reg][hilo].ee.fntight->Fill(ElPt[0], fabs(ElEta[0]), singleElPrescale() * puweight);
			if(S->datamc > 0) S->region[reg][hilo].ee.sst_origin->Fill(elIndexToBin(0)-0.5, puweight);
		}
		if( isLooseElectron(0) ){
			S->region[reg][hilo].ee.fnloose->Fill(ElPt[0], fabs(ElEta[0]), singleElPrescale() * puweight);
			if(S->datamc > 0) S->region[reg][hilo].ee.ssl_origin->Fill(elIndexToBin(0)-0.5, puweight);
		}
	}
	int elind;
	if(doubleElTrigger() && isZElElEvent(elind)){
		if( isTightElectron(elind) ){
			S->region[reg][hilo].ee.pntight->Fill(ElPt[elind], fabs(ElEta[elind]), puweight);
			if(S->datamc > 0) S->region[reg][hilo].ee.zt_origin->Fill(elIndexToBin(elind)-0.5, puweight);
		}
		if( isLooseElectron(elind) ){
			S->region[reg][hilo].ee.pnloose->Fill(ElPt[elind], fabs(ElEta[elind]), puweight);
			if(S->datamc > 0) S->region[reg][hilo].ee.zl_origin->Fill(elIndexToBin(elind)-0.5, puweight);
		}
	}

	// EMu Channel
	fCurrentChannel = ElMu;
	int mu(-1), el(-1);
	if(elmuSignalTrigger()){
		if(fDoCounting) fCounter[ElMu].fill(fEMCutNames[2]);
		if( isSSLLElMuEvent(mu, el) ){
			if(  isTightElectron(el) &&  isTightMuon(mu) ){ // Tight-tight
				if(fDoCounting) fCounter[ElMu].fill(fEMCutNames[14]);
				if(fDoCounting) fCounter[ElMu].fill(fEMCutNames[15]);
				if(fDoCounting) fCounter[ElMu].fill(fEMCutNames[16]);
				S->region[reg][hilo].em.nt20_pt ->Fill(MuPt [mu], ElPt [el], puweight);
				S->region[reg][hilo].em.nt20_eta->Fill(fabs(MuEta[mu]), fabs(ElEta[el]), puweight);
				if(S->datamc == 0 && reg == Baseline && hilo == HighPt){
					fOUTSTREAM << Form("%12s: ElMu - run %6.0d / ls %5.0d / ev %11.0d - HT(#J/#bJ) %6.2f(%1d/%1d) MET %6.2f MT2 %6.2f Pt1 %6.2f Pt2 %6.2f Charge %2d", S->sname.Data(), Run, LumiSec, Event, getHT(), getNJets(), getNBTags(), pfMET, getMT2(mu,el,3), MuPt[mu], ElPt[el], ElCharge[el]) << endl;
				}
				
				if(S->datamc > 0){
					S->region[reg][hilo].em.nt11_origin->Fill(muIndexToBin(mu)-0.5, elIndexToBin(el)-0.5, puweight);
					if(isPromptMuon(mu) && isPromptElectron(el)) S->region[reg][hilo].em.nt2pp_pt->Fill(MuPt[mu], ElPt[el], puweight);
					if(isPromptMuon(mu) &&   isFakeElectron(el)) S->region[reg][hilo].em.nt2pf_pt->Fill(MuPt[mu], ElPt[el], puweight);
					if(  isFakeMuon(mu) && isPromptElectron(el)) S->region[reg][hilo].em.nt2fp_pt->Fill(MuPt[mu], ElPt[el], puweight);
					if(  isFakeMuon(mu) &&   isFakeElectron(el)) S->region[reg][hilo].em.nt2ff_pt->Fill(MuPt[mu], ElPt[el], puweight);
				}
			}
			if( !isTightElectron(el) &&  isTightMuon(mu) ){ // Tight-loose
				if(fDoCounting) fCounter[ElMu].fill(fEMCutNames[14]);
				S->region[reg][hilo].em.nt10_pt ->Fill(MuPt [mu], ElPt [el], puweight);
				S->region[reg][hilo].em.nt10_eta->Fill(fabs(MuEta[mu]), fabs(ElEta[el]), puweight);
				if(S->datamc > 0) S->region[reg][hilo].em.nt10_origin->Fill(muIndexToBin(mu)-0.5, elIndexToBin(el)-0.5, puweight);
			}
			if(  isTightElectron(el) && !isTightMuon(mu) ){ // Loose-tight
				if(fDoCounting) fCounter[ElMu].fill(fEMCutNames[15]);
				S->region[reg][hilo].em.nt01_pt ->Fill(MuPt [mu], ElPt [el], puweight); // muon always in x axis for e/mu
				S->region[reg][hilo].em.nt01_eta->Fill(fabs(MuEta[mu]), fabs(ElEta[el]), puweight);
				if(S->datamc > 0) S->region[reg][hilo].em.nt01_origin->Fill(muIndexToBin(mu)-0.5, elIndexToBin(el)-0.5, puweight);
			}
			if( !isTightElectron(0) && !isTightMuon(0) ){ // Loose-loose
				S->region[reg][hilo].em.nt00_pt ->Fill(MuPt [mu], ElPt [el], puweight);
				S->region[reg][hilo].em.nt00_eta->Fill(fabs(MuEta[mu]), fabs(ElEta[el]), puweight);
				if(S->datamc > 0) S->region[reg][hilo].em.nt00_origin->Fill(muIndexToBin(mu)-0.5, elIndexToBin(el)-0.5, puweight);
			}
			if(S->datamc > 0){
				if(isPromptMuon(mu) && isPromptElectron(el)) S->region[reg][hilo].em.npp_pt->Fill(MuPt[mu], ElPt[el], puweight);
				if(isPromptMuon(mu) &&   isFakeElectron(el)) S->region[reg][hilo].em.npf_pt->Fill(MuPt[mu], ElPt[el], puweight);
				if(  isFakeMuon(mu) && isPromptElectron(el)) S->region[reg][hilo].em.nfp_pt->Fill(MuPt[mu], ElPt[el], puweight);
				if(  isFakeMuon(mu) &&   isFakeElectron(el)) S->region[reg][hilo].em.nff_pt->Fill(MuPt[mu], ElPt[el], puweight);
			}
		}
		resetHypLeptons();
	}

	///////////////////////////////////////////////////
	// OS YIELDS
	fDoCounting = false;
	fChargeSwitch = 1;

	// EE Channel
	fCurrentChannel = Elec;
	if(elelSignalTrigger()){
		if( isSSLLElEvent(el1, el2) ){ // this selects now OS events with the exact same cuts
			if(  isTightElectron(el1) &&  isTightElectron(el2) ){ // Tight-tight
				if( isBarrelElectron(el1) &&  isBarrelElectron(el2)) S->region[reg][hilo].ee.nt20_OS_BB_pt->Fill(ElPt[el1], ElPt[el2], puweight);
				if(!isBarrelElectron(el1) && !isBarrelElectron(el2)) S->region[reg][hilo].ee.nt20_OS_EE_pt->Fill(ElPt[el1], ElPt[el2], puweight);
				if( isBarrelElectron(el1) && !isBarrelElectron(el2)) S->region[reg][hilo].ee.nt20_OS_EB_pt->Fill(ElPt[el1], ElPt[el2], puweight);
				if(!isBarrelElectron(el1) &&  isBarrelElectron(el2)) S->region[reg][hilo].ee.nt20_OS_EB_pt->Fill(ElPt[el1], ElPt[el2], puweight);
			}
		}
		resetHypLeptons();
	}

	// EMu Channel
	fCurrentChannel = ElMu;
	if(elmuSignalTrigger()){
		if( isSSLLElMuEvent(mu, el) ){ // this selects now OS events with the exact same cuts
			if(  isTightElectron(el) &&  isTightMuon(mu) ){ // Tight-tight
				if( isBarrelElectron(el)) S->region[reg][hilo].em.nt20_OS_BB_pt->Fill(MuPt[mu], ElPt[el], puweight);
				if(!isBarrelElectron(el)) S->region[reg][hilo].em.nt20_OS_EE_pt->Fill(MuPt[mu], ElPt[el], puweight);
			}
		}
		resetHypLeptons();
	}
	fChargeSwitch = 0;
	fDoCounting = false;
	resetHypLeptons();
}
void SSDLDumper::fillDiffYields(Sample *S){
	///////////////////////////////////////////////////
	// Set custom event selections here:
	fC_minMu1pt = Region::minMu1pt[HighPt];
	fC_minMu2pt = Region::minMu2pt[HighPt];
	fC_minEl1pt = Region::minEl1pt[HighPt];
	fC_minEl2pt = Region::minEl2pt[HighPt];
	fC_minMet = 30.;
	fC_maxMet = 7000.;
	fC_maxHT  = 7000.;
	///////////////////////////////////////////////////
	// SS YIELDS
	// MuMu Channel
	resetHypLeptons();
	fDoCounting = false;
	fCurrentChannel = Muon;
	float puweight = PUWeight;
	if(S->datamc == 4) puweight = 1; // fix for samples with no pileup
	int mu1(-1), mu2(-1);
	if(mumuSignalTrigger()){ // Trigger selection
		fC_minHT  = 0.;
		fC_minNjets = 0;
		if(isSSLLMuEvent(mu1, mu2)){ // Same-sign loose-loose di muon event
			float HT  = getHT();
			float NJ  = getNJets() + 0.5;
			if(  isTightMuon(mu1) &&  isTightMuon(mu2) ){ // Tight-tight
				S->diffyields[Muon].hnt11[0]->Fill(HT,    puweight);
				S->diffyields[Muon].hnt11[2]->Fill(NJ,    puweight);
			}
			if(  isTightMuon(mu1) && !isTightMuon(mu2) ){ // Tight-loose
				S->diffyields[Muon].hnt10[0]->Fill(HT,    puweight);
				S->diffyields[Muon].hnt10[2]->Fill(NJ,    puweight);
			}
			if( !isTightMuon(mu1) &&  isTightMuon(mu2) ){ // Loose-tight
				S->diffyields[Muon].hnt01[0]->Fill(HT,    puweight);
				S->diffyields[Muon].hnt01[2]->Fill(NJ,    puweight);
			}
			if( !isTightMuon(mu1) && !isTightMuon(mu2) ){ // Loose-loose
				S->diffyields[Muon].hnt00[0]->Fill(HT,    puweight);
				S->diffyields[Muon].hnt00[2]->Fill(NJ,    puweight);
			}
		}
		resetHypLeptons();
		// Ask for 2 jets / 80 GeV HT when obtaining yields for met and mt2 plots
		fC_minHT =  80.;
		fC_minNjets = 2;
		if(isSSLLMuEvent(mu1, mu2)){ // Same-sign loose-loose di muon event
			float MT2 = getMT2(mu1, mu2, 1);
			if(  isTightMuon(mu1) &&  isTightMuon(mu2) ){ // Tight-tight
				S->diffyields[Muon].hnt11[1]->Fill(pfMET, puweight);
				S->diffyields[Muon].hnt11[3]->Fill(MT2,   puweight);
				S->diffyields[Muon].hnt11[4]->Fill(MuPt[mu1], puweight);
				S->diffyields[Muon].hnt11[5]->Fill(MuPt[mu2], puweight);
			}
			if(  isTightMuon(mu1) && !isTightMuon(mu2) ){ // Tight-loose
				S->diffyields[Muon].hnt10[1]->Fill(pfMET, puweight);
				S->diffyields[Muon].hnt10[3]->Fill(MT2,   puweight);
				S->diffyields[Muon].hnt10[4]->Fill(MuPt[mu1], puweight);
				S->diffyields[Muon].hnt10[5]->Fill(MuPt[mu2], puweight);
			}
			if( !isTightMuon(mu1) &&  isTightMuon(mu2) ){ // Loose-tight
				S->diffyields[Muon].hnt01[1]->Fill(pfMET, puweight);
				S->diffyields[Muon].hnt01[3]->Fill(MT2,   puweight);
				S->diffyields[Muon].hnt01[4]->Fill(MuPt[mu1], puweight);
				S->diffyields[Muon].hnt01[5]->Fill(MuPt[mu2], puweight);
			}
			if( !isTightMuon(mu1) && !isTightMuon(mu2) ){ // Loose-loose
				S->diffyields[Muon].hnt00[1]->Fill(pfMET, puweight);
				S->diffyields[Muon].hnt00[3]->Fill(MT2,   puweight);
				S->diffyields[Muon].hnt00[4]->Fill(MuPt[mu1], puweight);
				S->diffyields[Muon].hnt00[5]->Fill(MuPt[mu2], puweight);
			}
		}
		resetHypLeptons();
	}

	// EE Channel
	fCurrentChannel = Elec;
	int el1(-1), el2(-1);
	if(elelSignalTrigger()){
		fC_minHT  = 0.;
		fC_minNjets = 0;
		if( isSSLLElEvent(el1, el2) ){
			float HT  = getHT();
			float NJ  = getNJets() + 0.5;
			if(  isTightElectron(el1) &&  isTightElectron(el2) ){ // Tight-tight
				S->diffyields[Elec].hnt11[0]->Fill(HT,    puweight);
				S->diffyields[Elec].hnt11[2]->Fill(NJ,    puweight);
			}
			if(  isTightElectron(el1) && !isTightElectron(el2) ){ // Tight-loose
				S->diffyields[Elec].hnt10[0]->Fill(HT,    puweight);
				S->diffyields[Elec].hnt10[2]->Fill(NJ,    puweight);
			}
			if( !isTightElectron(el1) &&  isTightElectron(el2) ){ // Loose-tight
				S->diffyields[Elec].hnt01[0]->Fill(HT,    puweight);
				S->diffyields[Elec].hnt01[2]->Fill(NJ,    puweight);
			}
			if( !isTightElectron(el1) && !isTightElectron(el2) ){ // Loose-loose
				S->diffyields[Elec].hnt00[0]->Fill(HT,    puweight);
				S->diffyields[Elec].hnt00[2]->Fill(NJ,    puweight);
			}
		}
		resetHypLeptons();
		// Ask for 2 jets / 80 GeV HT when obtaining yields for met and mt2 plots
		fC_minHT =  80.;
		fC_minNjets = 2;
		if( isSSLLElEvent(el1, el2) ){
			float MT2 = getMT2(el1, el2, 2);
			if(  isTightElectron(el1) &&  isTightElectron(el2) ){ // Tight-tight
				S->diffyields[Elec].hnt11[1]->Fill(pfMET, puweight);
				S->diffyields[Elec].hnt11[3]->Fill(MT2,   puweight);
				S->diffyields[Elec].hnt11[4]->Fill(ElPt[el1], puweight);
				S->diffyields[Elec].hnt11[5]->Fill(ElPt[el2], puweight);
			}
			if(  isTightElectron(el1) && !isTightElectron(el2) ){ // Tight-loose
				S->diffyields[Elec].hnt10[1]->Fill(pfMET, puweight);
				S->diffyields[Elec].hnt10[3]->Fill(MT2,   puweight);
				S->diffyields[Elec].hnt10[4]->Fill(ElPt[el1], puweight);
				S->diffyields[Elec].hnt10[5]->Fill(ElPt[el2], puweight);
			}
			if( !isTightElectron(el1) &&  isTightElectron(el2) ){ // Loose-tight
				S->diffyields[Elec].hnt01[1]->Fill(pfMET, puweight);
				S->diffyields[Elec].hnt01[3]->Fill(MT2,   puweight);
				S->diffyields[Elec].hnt01[4]->Fill(ElPt[el1], puweight);
				S->diffyields[Elec].hnt01[5]->Fill(ElPt[el2], puweight);
			}
			if( !isTightElectron(el1) && !isTightElectron(el2) ){ // Loose-loose
				S->diffyields[Elec].hnt00[1]->Fill(pfMET, puweight);
				S->diffyields[Elec].hnt00[3]->Fill(MT2,   puweight);
				S->diffyields[Elec].hnt00[4]->Fill(ElPt[el1], puweight);
				S->diffyields[Elec].hnt00[5]->Fill(ElPt[el2], puweight);
			}
		}
		resetHypLeptons();
	}

	// EMu Channel
	fCurrentChannel = ElMu;
	int mu(-1), el(-1);
	if(elmuSignalTrigger()){
		fC_minHT  = 0.;
		fC_minNjets = 0;
		if( isSSLLElMuEvent(mu, el) ){
			float HT  = getHT();
			float NJ  = getNJets() + 0.5;
			if(  isTightElectron(el) &&  isTightMuon(mu) ){ // Tight-tight
				S->diffyields[ElMu].hnt11[0]->Fill(HT,    puweight);
				S->diffyields[ElMu].hnt11[2]->Fill(NJ,    puweight);
			}
			if( !isTightElectron(el) &&  isTightMuon(mu) ){ // Tight-loose
				S->diffyields[ElMu].hnt10[0]->Fill(HT,    puweight);
				S->diffyields[ElMu].hnt10[2]->Fill(NJ,    puweight);
			}
			if(  isTightElectron(el) && !isTightMuon(mu) ){ // Loose-tight
				S->diffyields[ElMu].hnt01[0]->Fill(HT,    puweight);
				S->diffyields[ElMu].hnt01[2]->Fill(NJ,    puweight);
			}
			if( !isTightElectron(0) && !isTightMuon(0) ){ // Loose-loose
				S->diffyields[ElMu].hnt00[0]->Fill(HT,    puweight);
				S->diffyields[ElMu].hnt00[2]->Fill(NJ,    puweight);
			}
		}
		resetHypLeptons();
		fC_minHT =  80.;
		fC_minNjets = 2;		
		if( isSSLLElMuEvent(mu, el) ){
			float MT2 = getMT2(mu, el, 3);
			float ptmax = MuPt[mu];
			float ptmin = ElPt[el];
			if(ptmin > ptmax){
				ptmin = MuPt[mu];
				ptmax = ElPt[el];
			}
			if(  isTightElectron(el) &&  isTightMuon(mu) ){ // Tight-tight
				S->diffyields[ElMu].hnt11[1]->Fill(pfMET, puweight);
				S->diffyields[ElMu].hnt11[3]->Fill(MT2,   puweight);
				S->diffyields[ElMu].hnt11[4]->Fill(ptmax, puweight);
				S->diffyields[ElMu].hnt11[5]->Fill(ptmin, puweight);
			}
			if( !isTightElectron(el) &&  isTightMuon(mu) ){ // Tight-loose
				S->diffyields[ElMu].hnt10[1]->Fill(pfMET, puweight);
				S->diffyields[ElMu].hnt10[3]->Fill(MT2,   puweight);
				S->diffyields[ElMu].hnt10[4]->Fill(ptmax, puweight);
				S->diffyields[ElMu].hnt10[5]->Fill(ptmin, puweight);
			}
			if(  isTightElectron(el) && !isTightMuon(mu) ){ // Loose-tight
				S->diffyields[ElMu].hnt01[1]->Fill(pfMET, puweight);
				S->diffyields[ElMu].hnt01[3]->Fill(MT2,   puweight);
				S->diffyields[ElMu].hnt01[4]->Fill(ptmax, puweight);
				S->diffyields[ElMu].hnt01[5]->Fill(ptmin, puweight);
			}
			if( !isTightElectron(0) && !isTightMuon(0) ){ // Loose-loose
				S->diffyields[ElMu].hnt00[1]->Fill(pfMET, puweight);
				S->diffyields[ElMu].hnt00[3]->Fill(MT2,   puweight);
				S->diffyields[ElMu].hnt00[4]->Fill(ptmax, puweight);
				S->diffyields[ElMu].hnt00[5]->Fill(ptmin, puweight);
			}
		}
		resetHypLeptons();
	}

	///////////////////////////////////////////////////
	// OS YIELDS
	fChargeSwitch = 1;

	// EE Channel
	fCurrentChannel = Elec;
	if(elelSignalTrigger()){
		fC_minHT    = 0.;
		fC_minNjets = 0;
		if( isSSLLElEvent(el1, el2) ){ // this selects now OS events with the exact same cuts
			if(  isTightElectron(el1) &&  isTightElectron(el2) ){ // Tight-tight
				float HT  = getHT();
				float NJ  = getNJets() + 0.5;
				if( isBarrelElectron(el1) &&  isBarrelElectron(el2)){
					S->diffyields[Elec].hnt2_os_BB[0]->Fill(HT, puweight);
					S->diffyields[Elec].hnt2_os_BB[2]->Fill(NJ, puweight);
				}
				if(!isBarrelElectron(el1) && !isBarrelElectron(el2)){
					S->diffyields[Elec].hnt2_os_EE[0]->Fill(HT, puweight);
					S->diffyields[Elec].hnt2_os_EE[2]->Fill(NJ, puweight);
				}
				if( isBarrelElectron(el1) && !isBarrelElectron(el2)){
					S->diffyields[Elec].hnt2_os_EB[0]->Fill(HT, puweight);
					S->diffyields[Elec].hnt2_os_EB[2]->Fill(NJ, puweight);
				}
				if(!isBarrelElectron(el1) &&  isBarrelElectron(el2)){
					S->diffyields[Elec].hnt2_os_EB[0]->Fill(HT, puweight);
					S->diffyields[Elec].hnt2_os_EB[2]->Fill(NJ, puweight);
				}
			}
		}
		resetHypLeptons();
		fC_minHT    = 80.;
		fC_minNjets = 2;
		if( isSSLLElEvent(el1, el2) ){ // this selects now OS events with the exact same cuts
			if(  isTightElectron(el1) &&  isTightElectron(el2) ){ // Tight-tight
				float MT2 = getMT2(el1, el2, 2);
				if( isBarrelElectron(el1) &&  isBarrelElectron(el2)){
					S->diffyields[Elec].hnt2_os_BB[1]->Fill(pfMET, puweight);
					S->diffyields[Elec].hnt2_os_BB[3]->Fill(MT2  , puweight);
					S->diffyields[Elec].hnt2_os_BB[4]->Fill(ElPt[el1], puweight);
					S->diffyields[Elec].hnt2_os_BB[5]->Fill(ElPt[el2], puweight);
				}
				if(!isBarrelElectron(el1) && !isBarrelElectron(el2)){
					S->diffyields[Elec].hnt2_os_EE[1]->Fill(pfMET, puweight);
					S->diffyields[Elec].hnt2_os_EE[3]->Fill(MT2  , puweight);
					S->diffyields[Elec].hnt2_os_EE[4]->Fill(ElPt[el1], puweight);
					S->diffyields[Elec].hnt2_os_EE[5]->Fill(ElPt[el2], puweight);
				}
				if( isBarrelElectron(el1) && !isBarrelElectron(el2)){
					S->diffyields[Elec].hnt2_os_EB[1]->Fill(pfMET, puweight);
					S->diffyields[Elec].hnt2_os_EB[3]->Fill(MT2  , puweight);
					S->diffyields[Elec].hnt2_os_EB[4]->Fill(ElPt[el1], puweight);
					S->diffyields[Elec].hnt2_os_EB[5]->Fill(ElPt[el2], puweight);
				}
				if(!isBarrelElectron(el1) &&  isBarrelElectron(el2)){
					S->diffyields[Elec].hnt2_os_EB[1]->Fill(pfMET, puweight);
					S->diffyields[Elec].hnt2_os_EB[3]->Fill(MT2  , puweight);
					S->diffyields[Elec].hnt2_os_EB[4]->Fill(ElPt[el1], puweight);
					S->diffyields[Elec].hnt2_os_EB[5]->Fill(ElPt[el2], puweight);
				}
			}
		}
		resetHypLeptons();
	}

	// EMu Channel
	fCurrentChannel = ElMu;
	if(elmuSignalTrigger()){
		fC_minHT    = 0.;
		fC_minNjets = 0;
		if( isSSLLElMuEvent(mu, el) ){ // this selects now OS events with the exact same cuts
			if(  isTightElectron(el) &&  isTightMuon(mu) ){ // Tight-tight
				float HT  = getHT();
				float NJ  = getNJets() + 0.5;
				if( isBarrelElectron(el)){
					S->diffyields[ElMu].hnt2_os_BB[0]->Fill(HT, puweight);
					S->diffyields[ElMu].hnt2_os_BB[2]->Fill(NJ, puweight);
				}
				if(!isBarrelElectron(el)){
					S->diffyields[ElMu].hnt2_os_EE[0]->Fill(HT, puweight);
					S->diffyields[ElMu].hnt2_os_EE[2]->Fill(NJ, puweight);
				}
			}
		}
		resetHypLeptons();
		fC_minHT    = 80.;
		fC_minNjets = 2;
		if( isSSLLElMuEvent(mu, el) ){ // this selects now OS events with the exact same cuts
			if(  isTightElectron(el) &&  isTightMuon(mu) ){ // Tight-tight
				float MT2 = getMT2(mu, el, 3);
				float ptmax = MuPt[mu];
				float ptmin = ElPt[el];
				if(ptmin > ptmax){
					ptmin = MuPt[mu];
					ptmax = ElPt[el];
				}
				if( isBarrelElectron(el)){
					S->diffyields[ElMu].hnt2_os_BB[1]->Fill(pfMET, puweight);
					S->diffyields[ElMu].hnt2_os_BB[3]->Fill(MT2  , puweight);
					S->diffyields[ElMu].hnt2_os_BB[4]->Fill(ptmax, puweight);
					S->diffyields[ElMu].hnt2_os_BB[5]->Fill(ptmin, puweight);
				}
				if(!isBarrelElectron(el)){
					S->diffyields[ElMu].hnt2_os_EE[1]->Fill(pfMET, puweight);
					S->diffyields[ElMu].hnt2_os_EE[3]->Fill(MT2  , puweight);
					S->diffyields[ElMu].hnt2_os_EE[4]->Fill(ptmax, puweight);
					S->diffyields[ElMu].hnt2_os_EE[5]->Fill(ptmin, puweight);
				}
			}
		}
		resetHypLeptons();
	}
	fChargeSwitch = 0;
	resetHypLeptons();
}
void SSDLDumper::fillRatioPlots(Sample *S){
	resetHypLeptons();
	fDoCounting = false;
	fChargeSwitch = 0;
	float puweight = PUWeight;
	if(S->datamc == 4) puweight = 1; // fix for samples with no pileup

	// Reset event selections to baseline:
	fC_minMu1pt = Region::minMu1pt[HighPt];
	fC_minMu2pt = Region::minMu2pt[HighPt];
	fC_minEl1pt = Region::minEl1pt[HighPt];
	fC_minEl2pt = Region::minEl2pt[HighPt];

	fC_minHT    = Region::minHT   [Baseline];
	fC_minMet   = Region::minMet  [Baseline];
	fC_maxHT    = Region::maxHT   [Baseline];
	fC_maxMet   = Region::maxMet  [Baseline];
	fC_minNjets = Region::minNjets[Baseline];

	FRatioPlots *RP0 = &S->ratioplots[0];
	FRatioPlots *RP1 = &S->ratioplots[1];

	fCurrentChannel = Muon;
	if(singleMuTrigger()){
		if(isSigSupMuEvent()){
			if( isTightMuon(0) ){
				RP0->ntight[0]->Fill(getNJets(),               puweight);
				RP0->ntight[1]->Fill(getHT(),                  puweight);
				RP0->ntight[2]->Fill(getMaxJPt(),              puweight);
				RP0->ntight[3]->Fill(NVrtx,                    puweight);
				RP0->ntight[4]->Fill(getClosestJetPt(0, Muon), puweight);
				RP0->ntight[5]->Fill(getAwayJetPt(0, Muon),    puweight);
				RP0->ntight[6]->Fill(getNBTags(),              puweight);
			}
			if( isLooseMuon(0) ){
				RP0->nloose[0]->Fill(getNJets(),               puweight);
				RP0->nloose[1]->Fill(getHT(),                  puweight);
				RP0->nloose[2]->Fill(getMaxJPt(),              puweight);
				RP0->nloose[3]->Fill(NVrtx,                    puweight);
				RP0->nloose[4]->Fill(getClosestJetPt(0, Muon), puweight);
				RP0->nloose[5]->Fill(getAwayJetPt(0, Muon),    puweight);
				RP0->nloose[6]->Fill(getNBTags(),              puweight);
			}
		}
		fC_maxMet_Control = 1000.;
		if(isSigSupMuEvent()){
			if( isTightMuon(0) ){
				RP0->ntight[7]->Fill(pfMET,                    puweight);
			}
			if( isLooseMuon(0) ){
				RP0->nloose[7]->Fill(pfMET,                    puweight);
			}
		}		
		fC_maxMet_Control = 20.;
		fC_maxMt_Control = 1000.;
		if(isSigSupMuEvent()){
			if( isTightMuon(0) ){
				RP0->ntight[8]->Fill(MuMT[0],                  puweight);
			}
			if( isLooseMuon(0) ){
				RP0->nloose[8]->Fill(MuMT[0],                  puweight);
			}
		}		
		fC_maxMt_Control = 20.;
	}
	resetHypLeptons();
	if(singleElTrigger()){
		if(isSigSupElEvent()){
			if( isTightElectron(0) ){
				RP1->ntight[0]->Fill(getNJets(),                   puweight);
				RP1->ntight[1]->Fill(getHT(),                      puweight);
				RP1->ntight[2]->Fill(getMaxJPt(),                  puweight);
				RP1->ntight[3]->Fill(NVrtx,                        puweight);
				RP1->ntight[4]->Fill(getClosestJetPt(0, Elec), puweight);
				RP1->ntight[5]->Fill(getAwayJetPt(0, Elec),    puweight);
				RP1->ntight[6]->Fill(getNBTags(),                  puweight);
			}
			if( isLooseElectron(0) ){
				RP1->nloose[0]->Fill(getNJets(),                   puweight);
				RP1->nloose[1]->Fill(getHT(),                      puweight);
				RP1->nloose[2]->Fill(getMaxJPt(),                  puweight);
				RP1->nloose[3]->Fill(NVrtx,                        puweight);
				RP1->nloose[4]->Fill(getClosestJetPt(0, Elec), puweight);
				RP1->nloose[5]->Fill(getAwayJetPt(0, Elec),    puweight);
				RP1->nloose[6]->Fill(getNBTags(),                  puweight);
			}
		}
		fC_maxMet_Control = 1000.;
		if(isSigSupElEvent()){
			if( isTightElectron(0) ){
				RP1->ntight[7]->Fill(pfMET,                    puweight);
			}
			if( isLooseElectron(0) ){
				RP1->nloose[7]->Fill(pfMET,                    puweight);
			}
		}		
		fC_maxMet_Control = 20.;
		fC_maxMt_Control = 1000.;
		if(isSigSupElEvent()){
			if( isTightElectron(0) ){
				RP1->ntight[8]->Fill(ElMT[0],                  puweight);
			}
			if( isLooseElectron(0) ){
				RP1->nloose[8]->Fill(ElMT[0],                  puweight);
			}
		}		
		fC_maxMt_Control = 20.;
	}
	resetHypLeptons();
}
void SSDLDumper::fillKinPlots(Sample *S, gHiLoSwitch hilo){
	resetHypLeptons();
	fDoCounting = false;
	KinPlots *KP0 = &S->kinplots[0][hilo];
	KinPlots *KP1 = &S->kinplots[1][hilo];
	KinPlots *KP2 = &S->kinplots[2][hilo];
	int ind1(-1), ind2(-1);
	int mu(-1), el(-1); // for e/mu channel, easier readability

	float puweight = PUWeight;
	if(S->datamc == 4) puweight = 1; // fix for samples with no pileup

	///////////////////////////////////////////////////
	// Set custom event selections here:
	fC_minMu1pt = Region::minMu1pt[hilo];
	fC_minMu2pt = Region::minMu2pt[hilo];
	fC_minEl1pt = Region::minEl1pt[hilo];
	fC_minEl2pt = Region::minEl2pt[hilo];
	fC_minHT    = Region::minHT   [Baseline];
	fC_maxHT    = Region::maxHT   [Baseline];
	fC_minMet   = Region::minMet  [Baseline];
	fC_maxMet   = Region::maxMet  [Baseline];
	fC_minNjets = Region::minNjets[Baseline];
	if(hilo == LowPt) fC_minHT = 200.;

	////////////////////////////////////////////////////////////////////////////////////////////////////////
	// MUMU CHANNEL:  //////////////////////////////////////////////////////////////////////////////////////
	if(mumuSignalTrigger() && abs(isSSLLEvent(ind1, ind2)) == 1){ // trigger && select loose mu/mu pair
		if(MuPt[ind1] > fC_minMu1pt && MuPt[ind2] > fC_minMu2pt){ // pt cuts

			// Fill histos
			KP0->hmetvsht->Fill(getHT(), pfMET, puweight);
			
			KP0->hvar[0]->Fill(getHT(),    puweight);
			KP0->hvar[1]->Fill(pfMET,      puweight);
			KP0->hvar[2]->Fill(getNJets(), puweight);
			KP0->hvar[3]->Fill(MuPt[ind1], puweight);
			KP0->hvar[4]->Fill(MuPt[ind2], puweight);
			TLorentzVector p1, p2;
			p1.SetPtEtaPhiM(MuPt[ind1], MuEta[ind1], MuPhi[ind1], gMMU);
			p2.SetPtEtaPhiM(MuPt[ind2], MuEta[ind2], MuPhi[ind2], gMMU);
			float mass = (p1+p2).M();
			KP0->hvar[5]->Fill(mass, puweight); // SF
			KP0->hvar[6]->Fill(mass, puweight); // MM
			KP0->hvar[9]->Fill(getMT2(ind1, ind2, 1), puweight);

			if(isTightMuon(ind1) && isTightMuon(ind2)){ // tight-tight
				KP1->hmetvsht->Fill(getHT(), pfMET, puweight);
				KP1->hvar[0]->Fill(getHT(),               puweight);
				KP1->hvar[1]->Fill(pfMET,                 puweight);
				KP1->hvar[2]->Fill(getNJets(),            puweight);
				KP1->hvar[3]->Fill(MuPt[ind1],            puweight);
				KP1->hvar[4]->Fill(MuPt[ind2],            puweight);
				KP1->hvar[5]->Fill(mass,                  puweight); // SF
				KP1->hvar[6]->Fill(mass,                  puweight); // MM
				KP1->hvar[9]->Fill(getMT2(ind1, ind2, 1), puweight);
				if(isSSLLMuEvent(ind1, ind2)){ // signal region
					KP2->hmetvsht->Fill(getHT(), pfMET, puweight);
					if(hilo == HighPt){ // Store signal events
						fSigEv_HI_MM_HT .push_back(getHT());
						fSigEv_HI_MM_MET.push_back(pfMET);
					} else{
						fSigEv_LO_MM_HT .push_back(getHT());
						fSigEv_LO_MM_MET.push_back(pfMET);						
					}					
					
					KP2->hvar[0]->Fill(getHT(),    puweight);
					KP2->hvar[1]->Fill(pfMET,      puweight);
					KP2->hvar[2]->Fill(getNJets(), puweight);
					KP2->hvar[3]->Fill(MuPt[ind1], puweight);
					KP2->hvar[4]->Fill(MuPt[ind2], puweight);
					KP2->hvar[5]->Fill(mass,       puweight); // SF
					KP2->hvar[6]->Fill(mass,       puweight); // MM					
					KP2->hvar[9]->Fill(getMT2(ind1, ind2, 1), puweight);
				}
			}
		}
		resetHypLeptons();
	}

	////////////////////////////////////////////////////////////////////////////////////////////////////////
	// EE CHANNEL:  ////////////////////////////////////////////////////////////////////////////////////////
	else if(elelSignalTrigger() && abs(isSSLLEvent(ind1, ind2)) == 2){ // trigger && select loose e/e pair
		if(ElPt[ind1] > fC_minEl1pt && ElPt[ind2] > fC_minEl2pt){ // pt cuts
			// Fill histos
			KP0->hmetvsht->Fill(getHT(), pfMET, puweight);
			KP0->hvar[0]->Fill(getHT(),    puweight);
			KP0->hvar[1]->Fill(pfMET,      puweight);
			KP0->hvar[2]->Fill(getNJets(), puweight);
			KP0->hvar[3]->Fill(ElPt[ind1], puweight);
			KP0->hvar[4]->Fill(ElPt[ind2], puweight);
			TLorentzVector p1, p2;
			p1.SetPtEtaPhiM(ElPt[ind1], ElEta[ind1], ElPhi[ind1], gMEL);
			p2.SetPtEtaPhiM(ElPt[ind2], ElEta[ind2], ElPhi[ind2], gMEL);
			float mass = (p1+p2).M();
			KP0->hvar[5]->Fill(mass, puweight); // SF
			KP0->hvar[7]->Fill(mass, puweight); // MM
			KP0->hvar[9]->Fill(getMT2(ind1, ind2, 2), puweight);
			if(isTightElectron(ind1) && isTightElectron(ind2)){ // tight-tight
				KP1->hmetvsht->Fill(getHT(), pfMET, puweight);
				KP1->hvar[0]->Fill(getHT(),               puweight);
				KP1->hvar[1]->Fill(pfMET,                 puweight);
				KP1->hvar[2]->Fill(getNJets(),            puweight);
				KP1->hvar[3]->Fill(ElPt[ind1],            puweight);
				KP1->hvar[4]->Fill(ElPt[ind2],            puweight);
				KP1->hvar[5]->Fill(mass,                  puweight); // SF
				KP1->hvar[7]->Fill(mass,                  puweight); // MM
				KP1->hvar[9]->Fill(getMT2(ind1, ind2, 2), puweight);
				if(isSSLLElEvent(ind1, ind2)){
					KP2->hmetvsht->Fill(getHT(), pfMET, puweight);
					if(hilo == HighPt){ // Store signal events
						fSigEv_HI_EE_HT .push_back(getHT());
						fSigEv_HI_EE_MET.push_back(pfMET);
					} else{
						fSigEv_LO_EE_HT .push_back(getHT());
						fSigEv_LO_EE_MET.push_back(pfMET);						
					}
					KP2->hvar[0]->Fill(getHT(),               puweight);
					KP2->hvar[1]->Fill(pfMET,                 puweight);
					KP2->hvar[2]->Fill(getNJets(),            puweight);
					KP2->hvar[3]->Fill(ElPt[ind1],            puweight);
					KP2->hvar[4]->Fill(ElPt[ind2],            puweight);
					KP2->hvar[5]->Fill(mass,                  puweight); // SF
					KP2->hvar[7]->Fill(mass,                  puweight); // MM							
					KP2->hvar[9]->Fill(getMT2(ind1, ind2, 2), puweight);
				}
			}
		}
		resetHypLeptons();
	}

	////////////////////////////////////////////////////////////////////////////////////////////////////////
	// EMU CHANNEL:  ///////////////////////////////////////////////////////////////////////////////////////
	else if(elmuSignalTrigger() && abs(isSSLLEvent(mu, el)) == 3){ // trigger && select loose e/mu pair
		if( (MuPt[mu] > fC_minMu1pt && ElPt[el] > fC_minEl2pt) || (MuPt[mu] > fC_minMu2pt && ElPt[el] > fC_minEl1pt) ){
			// Fill histos
			KP0->hmetvsht->Fill(getHT(), pfMET, puweight);
			KP0->hvar[0]->Fill(getHT(),    puweight);
			KP0->hvar[1]->Fill(pfMET,      puweight);
			KP0->hvar[2]->Fill(getNJets(), puweight);
			float ptmax = MuPt[mu];
			float ptmin = ElPt[el];
			if(ptmin > ptmax){
				ptmin = MuPt[mu];
				ptmax = ElPt[el];
			}
			KP0->hvar[3]->Fill(ptmax, puweight);
			KP0->hvar[4]->Fill(ptmin, puweight);
			TLorentzVector p1, p2;
			p1.SetPtEtaPhiM(MuPt[mu], MuEta[mu], MuPhi[mu], gMMU);
			p2.SetPtEtaPhiM(ElPt[el], ElEta[el], ElPhi[el], gMEL);
			float mass = (p1+p2).M();
			KP0->hvar[8]->Fill(mass,                  puweight); // EM
			KP0->hvar[9]->Fill(getMT2(mu, el, 3), puweight);
			if(isTightMuon(mu) && isTightElectron(el)){ // tight-tight
				KP1->hmetvsht->Fill(getHT(), pfMET, puweight);
				KP1->hvar[0]->Fill(getHT(),               puweight);
				KP1->hvar[1]->Fill(pfMET,                 puweight);
				KP1->hvar[2]->Fill(getNJets(),            puweight);
				KP1->hvar[3]->Fill(ptmax,                 puweight);
				KP1->hvar[4]->Fill(ptmin,                 puweight);
				KP1->hvar[8]->Fill(mass,                  puweight); // EM
				KP1->hvar[9]->Fill(getMT2(mu, el, 3), puweight);
				if( isSSLLElMuEvent(mu, el) ){
					KP2->hmetvsht->Fill(getHT(), pfMET, puweight);
					if(hilo == HighPt){ // Store signal events
						fSigEv_HI_EM_HT .push_back(getHT());
						fSigEv_HI_EM_MET.push_back(pfMET);
					} else{
						fSigEv_LO_EM_HT .push_back(getHT());
						fSigEv_LO_EM_MET.push_back(pfMET);						
					}
					KP2->hvar[0]->Fill(getHT(),               puweight);
					KP2->hvar[1]->Fill(pfMET,                 puweight);
					KP2->hvar[2]->Fill(getNJets(),            puweight);
					KP2->hvar[3]->Fill(ptmax,                 puweight);
					KP2->hvar[4]->Fill(ptmin,                 puweight);
					KP2->hvar[8]->Fill(mass,                  puweight); // EM						
					KP2->hvar[9]->Fill(getMT2(mu, el, 3), puweight);
				}
			}
		}
	}
	resetHypLeptons();
	return;
}
void SSDLDumper::fillMuIsoPlots(Sample *S){
	resetHypLeptons();
	fDoCounting = false;
	IsoPlots *IP0 = &S->isoplots[0]; // mu

	// Reset event selections to baseline:
	fC_minMu1pt = Region::minMu1pt[HighPt];
	fC_minMu2pt = Region::minMu2pt[HighPt];
	fC_minEl1pt = Region::minEl1pt[HighPt];
	fC_minEl2pt = Region::minEl2pt[HighPt];

	fC_minHT    = Region::minHT   [Baseline];
	fC_maxHT    = Region::maxHT   [Baseline];
	fC_minMet   = Region::minMet  [Baseline];
	fC_maxMet   = Region::maxMet  [Baseline];
	fC_minNjets = Region::minNjets[Baseline];

	int muind1(-1), muind2(-1);
	if(hasLooseMuons(muind1, muind2) > 0){
		setHypLepton1(muind1, Muon);
		// Common trigger selection
		if(!singleMuTrigger()) return;
		float prescale = singleMuPrescale();
		float puweight = PUWeight;
		if(S->datamc == 4) puweight = 1; // fix for samples with no pileup
		float scale = prescale * puweight;

		// Common event selections
		if(!passesJet50Cut()) return; // make trigger 100% efficient

		// Common object selections
		if(!isLooseMuon(muind1)) return;
		if(MuPt[muind1] < fC_minMu2pt) return;
		if(MuPt[muind1] > gMuPt2bins[gNMuPt2bins]) return;

		////////////////////////////////////////////////////
		// MOST LOOSE SELECTION
		IP0->hiso[0]->Fill(MuIso[muind1], scale);
		for(size_t k = 0; k < gNMuPt2bins; ++k){
			if(MuPt[muind1] < gMuPt2bins[k]) continue;
			if(MuPt[muind1] > gMuPt2bins[k+1]) continue;
			IP0->hiso_pt[0][k]->Fill(MuIso[muind1], scale);
		}
		for(size_t k = 0; k < gNNVrtxBins; ++k){
			if(NVrtx <  gNVrtxBins[k]) continue;
			if(NVrtx >= gNVrtxBins[k+1]) continue;
			IP0->hiso_nv[0][k]->Fill(MuIso[muind1], scale);
		}

		////////////////////////////////////////////////////
		// SIGNAL SUPPRESSED SELECTION
		if(isSigSupMuEvent()){
			IP0->hiso[1]->Fill(MuIso[muind1], scale);
			for(size_t k = 0; k < gNMuPt2bins; ++k){
				if(MuPt[muind1] < gMuPt2bins[k]) continue;
				if(MuPt[muind1] > gMuPt2bins[k+1]) continue;
				IP0->hiso_pt[1][k]->Fill(MuIso[muind1], scale);
			}
			for(size_t k = 0; k < gNNVrtxBins; ++k){
				if(NVrtx <  gNVrtxBins[k]) continue;
				if(NVrtx >= gNVrtxBins[k+1]) continue;
				IP0->hiso_nv[1][k]->Fill(MuIso[muind1], scale);
			}
		}
		resetHypLeptons();
		////////////////////////////////////////////////////
	}
	return;
}
void SSDLDumper::fillElIsoPlots(Sample *S){
	resetHypLeptons();
	fDoCounting = false;
	IsoPlots *IP = &S->isoplots[1]; // mu

	// Reset event selections to baseline:
	fC_minMu1pt = Region::minMu1pt[HighPt];
	fC_minMu2pt = Region::minMu2pt[HighPt];
	fC_minEl1pt = Region::minEl1pt[HighPt];
	fC_minEl2pt = Region::minEl2pt[HighPt];

	fC_minHT    = Region::minHT   [Baseline];
	fC_maxHT    = Region::maxHT   [Baseline];
	fC_minMet   = Region::minMet  [Baseline];
	fC_maxMet   = Region::maxMet  [Baseline];
	fC_minNjets = Region::minNjets[Baseline];

	int elind1(-1), elind2(-1);
	if(hasLooseElectrons(elind1, elind2) > 0){
		setHypLepton1(elind1, Elec);
		// Common trigger selection
		if(!singleElTrigger()) return;
		float prescale = singleElPrescale();
		float puweight = PUWeight;
		if(S->datamc == 4) puweight = 1; // fix for samples with no pileup
		float scale = prescale * puweight;

		// Common event selections
		if(!passesJet50Cut()) return; // make trigger 100% efficient

		// Common object selections
		if(!isLooseElectron(elind1)) return;
		// if(ElIsGoodElId_WP80[elind1] != 1) return false; // apply tight ID for the iso plots?

		if(ElPt[elind1] < fC_minEl2pt) return;
		if(ElPt[elind1] > gElPt2bins[gNElPt2bins]) return;

		////////////////////////////////////////////////////
		// MOST LOOSE SELECTION
		IP->hiso[0]->Fill(ElRelIso[elind1], scale);
		for(size_t k = 0; k < gNElPt2bins; ++k){
			if(ElPt[elind1] < gElPt2bins[k]) continue;
			if(ElPt[elind1] > gElPt2bins[k+1]) continue;
			IP->hiso_pt[0][k]->Fill(ElRelIso[elind1], scale);
		}
		for(size_t k = 0; k < gNNVrtxBins; ++k){
			if(NVrtx <  gNVrtxBins[k]) continue;
			if(NVrtx >= gNVrtxBins[k+1]) continue;
			IP->hiso_nv[0][k]->Fill(ElRelIso[elind1], scale);
		}

		////////////////////////////////////////////////////
		// SIGNAL SUPPRESSED SELECTION
		if(isSigSupElEvent()){
			IP->hiso[1]->Fill(ElRelIso[elind1], scale);
			for(size_t k = 0; k < gNElPt2bins; ++k){
				if(ElPt[elind1] < gElPt2bins[k]) continue;
				if(ElPt[elind1] > gElPt2bins[k+1]) continue;
				IP->hiso_pt[1][k]->Fill(ElRelIso[elind1], scale);
			}
			for(size_t k = 0; k < gNNVrtxBins; ++k){
				if(NVrtx <  gNVrtxBins[k]) continue;
				if(NVrtx >= gNVrtxBins[k+1]) continue;
				IP->hiso_nv[1][k]->Fill(ElRelIso[elind1], scale);
			}
		}
		////////////////////////////////////////////////////
		resetHypLeptons();
	}
	return;
}

//____________________________________________________________________________
void SSDLDumper::storeNumbers(Sample *S, gChannel chan, gRegion reg){
	Channel *C;
	if(chan == Muon) C = &S->region[reg][HighPt].mm;
	if(chan == Elec) C = &S->region[reg][HighPt].ee;
	if(chan == ElMu) C = &S->region[reg][HighPt].em;
	S->numbers[reg][chan].nt2  = C->nt20_pt->GetEntries();
	S->numbers[reg][chan].nt10 = C->nt10_pt->GetEntries();
	S->numbers[reg][chan].nt01 = C->nt01_pt->GetEntries();
	S->numbers[reg][chan].nt0  = C->nt00_pt->GetEntries();
	if(chan != ElMu){
		S->numbers[reg][chan].nsst = C->fntight->GetEntries();
		S->numbers[reg][chan].nssl = C->fnloose->GetEntries();
		S->numbers[reg][chan].nzt  = C->pntight->GetEntries();
		S->numbers[reg][chan].nzl  = C->pnloose->GetEntries();
	}
}

//____________________________________________________________________________
void SSDLDumper::initCutNames(){
	// Muon channel
	fMMCutNames.push_back("All events"); //                            = 0
	fMMCutNames.push_back(" ... is good run"); //                      = 1
	fMMCutNames.push_back(" ... passes triggers"); //                  = 2
	fMMCutNames.push_back(" ... has 1 loose muon"); //                 = 3
	fMMCutNames.push_back(" ... has 2 loose muons"); //                = 4
	fMMCutNames.push_back(" ... has same-sign muons"); //              = 5
	fMMCutNames.push_back(" ... passes Z veto"); //                    = 6
	fMMCutNames.push_back(" ... passes Minv veto"); //                 = 7
	fMMCutNames.push_back(" ... has one jet > 50 GeV"); //             = 8
	fMMCutNames.push_back(" ... passes NJets cut"); //                 = 9
	fMMCutNames.push_back(" ... passes HT cut"); //                    = 10
	fMMCutNames.push_back(" ... passes MET cut"); //                   = 11
	fMMCutNames.push_back(" ... second muon passes pt cut"); //        = 12
	fMMCutNames.push_back(" ... first muon passes pt cut"); //         = 13
	fMMCutNames.push_back(" ... first muon passes tight cut"); //      = 14
	fMMCutNames.push_back(" ... second muon passes tight cut"); //     = 15
	fMMCutNames.push_back(" ... both muons pass tight cut"); //        = 16

	// Electron channel
	fEECutNames.push_back("All events"); //                            = 0
	fEECutNames.push_back(" ... is good run"); //                      = 1
	fEECutNames.push_back(" ... passes triggers"); //                  = 2
	fEECutNames.push_back(" ... has 1 loose electron"); //             = 3
	fEECutNames.push_back(" ... has 2 loose electrons"); //            = 4
	fEECutNames.push_back(" ... has same-sign electrons"); //          = 5
	fEECutNames.push_back(" ... passes Z veto"); //                    = 6
	fEECutNames.push_back(" ... passes Minv veto"); //                 = 7
	fEECutNames.push_back(" ... has one jet > 50 GeV"); //             = 8
	fEECutNames.push_back(" ... passes NJets cut"); //                 = 9
	fEECutNames.push_back(" ... passes HT cut"); //                    = 10
	fEECutNames.push_back(" ... passes MET cut"); //                   = 11
	fEECutNames.push_back(" ... second electron passes pt cut"); //    = 12
	fEECutNames.push_back(" ... first electron passes pt cut"); //     = 13
	fEECutNames.push_back(" ... first electron passes tight cut"); //  = 14
	fEECutNames.push_back(" ... second electron passes tight cut"); // = 15
	fEECutNames.push_back(" ... both electrons pass tight cut"); //    = 16
	
	// E-Mu channel
	fEMCutNames.push_back("All events"); //                            = 0
	fEMCutNames.push_back(" ... is good run"); //                      = 1
	fEMCutNames.push_back(" ... passes triggers"); //                  = 2
	fEMCutNames.push_back(" ... has a loose muon"); //                 = 3
	fEMCutNames.push_back(" ... has a loose electron"); //             = 4
	fEMCutNames.push_back(" ... has both"); //                         = 5
	fEMCutNames.push_back(" ... has same-sign electron muon pair"); // = 6
	fEMCutNames.push_back(" ... passes Z veto"); //                    = 7
	fEMCutNames.push_back(" ... has one jet > 50 GeV"); //             = 8
	fEMCutNames.push_back(" ... passes NJets cut"); //                 = 9
	fEMCutNames.push_back(" ... passes HT cut"); //                    = 10
	fEMCutNames.push_back(" ... passes MET cut"); //                   = 11
	fEMCutNames.push_back(" ... muon passes pt cut"); //               = 12
	fEMCutNames.push_back(" ... electron passes pt cut"); //           = 13
	fEMCutNames.push_back(" ... muon passes tight cut"); //            = 14
	fEMCutNames.push_back(" ... electron passes tight cut"); //        = 15
	fEMCutNames.push_back(" ... both e and mu pass tight cuts"); //    = 16
}
void SSDLDumper::initCounters(){
	fCounter[Muon].fill(fMMCutNames[0],  0.);
	fCounter[Muon].fill(fMMCutNames[1],  0.);
	fCounter[Muon].fill(fMMCutNames[2],  0.);
	fCounter[Muon].fill(fMMCutNames[3],  0.);
	fCounter[Muon].fill(fMMCutNames[4],  0.);
	fCounter[Muon].fill(fMMCutNames[5],  0.);
	fCounter[Muon].fill(fMMCutNames[6],  0.);
	fCounter[Muon].fill(fMMCutNames[7],  0.);
	fCounter[Muon].fill(fMMCutNames[8],  0.);
	fCounter[Muon].fill(fMMCutNames[9],  0.);
	fCounter[Muon].fill(fMMCutNames[10], 0.);
	fCounter[Muon].fill(fMMCutNames[11], 0.);
	fCounter[Muon].fill(fMMCutNames[12], 0.);
	fCounter[Muon].fill(fMMCutNames[13], 0.);
	fCounter[Muon].fill(fMMCutNames[14], 0.);
	fCounter[Muon].fill(fMMCutNames[15], 0.);
	fCounter[Muon].fill(fMMCutNames[16], 0.);

	fCounter[Elec].fill(fEECutNames[0],  0.);
	fCounter[Elec].fill(fEECutNames[1],  0.);
	fCounter[Elec].fill(fEECutNames[2],  0.);
	fCounter[Elec].fill(fEECutNames[3],  0.);
	fCounter[Elec].fill(fEECutNames[4],  0.);
	fCounter[Elec].fill(fEECutNames[5],  0.);
	fCounter[Elec].fill(fEECutNames[6],  0.);
	fCounter[Elec].fill(fEECutNames[7],  0.);
	fCounter[Elec].fill(fEECutNames[8],  0.);
	fCounter[Elec].fill(fEECutNames[9],  0.);
	fCounter[Elec].fill(fEECutNames[10], 0.);
	fCounter[Elec].fill(fEECutNames[11], 0.);
	fCounter[Elec].fill(fEECutNames[12], 0.);
	fCounter[Elec].fill(fEECutNames[13], 0.);
	fCounter[Elec].fill(fEECutNames[14], 0.);
	fCounter[Elec].fill(fEECutNames[15], 0.);
	fCounter[Elec].fill(fEECutNames[16], 0.);

	fCounter[ElMu].fill(fEMCutNames[0],  0.);
	fCounter[ElMu].fill(fEMCutNames[1],  0.);
	fCounter[ElMu].fill(fEMCutNames[2],  0.);
	fCounter[ElMu].fill(fEMCutNames[3],  0.);
	fCounter[ElMu].fill(fEMCutNames[4],  0.);
	fCounter[ElMu].fill(fEMCutNames[5],  0.);
	fCounter[ElMu].fill(fEMCutNames[6],  0.);
	fCounter[ElMu].fill(fEMCutNames[7],  0.);
	fCounter[ElMu].fill(fEMCutNames[8],  0.);
	fCounter[ElMu].fill(fEMCutNames[9],  0.);
	fCounter[ElMu].fill(fEMCutNames[10], 0.);
	fCounter[ElMu].fill(fEMCutNames[11], 0.);
	fCounter[ElMu].fill(fEMCutNames[12], 0.);
	fCounter[ElMu].fill(fEMCutNames[13], 0.);
	fCounter[ElMu].fill(fEMCutNames[14], 0.);
	fCounter[ElMu].fill(fEMCutNames[15], 0.);
	fCounter[ElMu].fill(fEMCutNames[16], 0.);
}
void SSDLDumper::fillCutFlowHistos(Sample *S){
	for(int i=0; i<fMMCutNames.size(); i++) S->cutFlowHisto[Muon]->SetBinContent(i+1, fCounter[Muon].counts(fMMCutNames[i]));
	for(int i=0; i<fEECutNames.size(); i++) S->cutFlowHisto[Elec]->SetBinContent(i+1, fCounter[Elec].counts(fEECutNames[i]));
	for(int i=0; i<fEMCutNames.size(); i++) S->cutFlowHisto[ElMu]->SetBinContent(i+1, fCounter[ElMu].counts(fEMCutNames[i]));	
}
void SSDLDumper::printCutFlow(gChannel chan, gSample indmin, gSample indmax){
	vector<string> names;
	if(chan == Muon) names = fMMCutNames;
	if(chan == Elec) names = fEECutNames;
	if(chan == ElMu) names = fEMCutNames;

	fOUTSTREAM << "------------------------------------------";
	for(int i = 0; i <= indmax-indmin; i++) fOUTSTREAM << "--------------";
	fOUTSTREAM << endl;

	fOUTSTREAM << " Cutname                                 | ";
	for(gSample i = indmin; i <= indmax; i=gSample(i+1)) fOUTSTREAM << setw(11) << fSamples[i]->sname << " | ";
	fOUTSTREAM << endl;

	fOUTSTREAM << "------------------------------------------";
	for(int i = 0; i <= indmax-indmin; i++) fOUTSTREAM << "--------------";
	fOUTSTREAM << endl;

	for( int c = 0; c < fMMCutNames.size(); c++ ){
		fOUTSTREAM << setw(40) << names[c] << " | ";

		for(gSample i = indmin; i <= indmax; i=gSample(i+1)){
			fOUTSTREAM << setw(11) << setprecision(11) << fSamples[i]->cutFlowHisto[chan]->GetBinContent(c+1) << " | ";
		}
		fOUTSTREAM << endl;
	}
	fOUTSTREAM << "------------------------------------------";
	for(int i = 0; i <= indmax-indmin; i++) fOUTSTREAM << "--------------";
	fOUTSTREAM << endl;	
}
void SSDLDumper::printCutFlows(TString filename){
	fOUTSTREAM.open(filename.Data(), ios::trunc);
	fOUTSTREAM << " Printing Cutflow for Mu/Mu channel..." << endl;
	// printCutFlow(Muon, DoubleMu1, EleHad2);
	printCutFlow(Muon, DoubleMu1, DoubleMu4);
	printCutFlow(Muon, TTJets, GJets200);
	printCutFlow(Muon, GVJets, ttbarW);
	printCutFlow(Muon, LM0, LM6);
	printCutFlow(Muon, LM7, LM13);
	printCutFlow(Muon, QCDMuEnr10, QCD470);
	printCutFlow(Muon, QCD600, QCD1000MG);
	fOUTSTREAM << endl << endl;
	
	fOUTSTREAM << " Printing Cutflow for E/Mu channel..." << endl;
	printCutFlow(ElMu, MuEG1, MuEG4);
	printCutFlow(ElMu, TTJets, GJets200);
	printCutFlow(ElMu, GVJets, ttbarW);
	printCutFlow(ElMu, LM0, LM6);
	printCutFlow(ElMu, LM7, LM13);
	printCutFlow(ElMu, QCDMuEnr10, QCD470);
	printCutFlow(ElMu, QCD600, QCD1000MG);
	fOUTSTREAM << endl << endl;
	
	fOUTSTREAM << " Printing Cutflow for E/E channel..." << endl;
	printCutFlow(Elec, DoubleEle1, DoubleEle4);
	printCutFlow(Elec, TTJets, GJets200);
	printCutFlow(Elec, GVJets, ttbarW);
	printCutFlow(Elec, LM0, LM6);
	printCutFlow(Elec, LM7, LM13);
	printCutFlow(Elec, QCDMuEnr10, QCD470);
	printCutFlow(Elec, QCD600, QCD1000MG);
	fOUTSTREAM << endl << endl;

	fOUTSTREAM.close();
}

//____________________________________________________________________________
void SSDLDumper::bookHistos(Sample *S){
	// Cut flow histos
	S->cutFlowHisto[Muon] = new TH1D("MMCutFlow", "MMCutFlow", fMMCutNames.size(), 0, fMMCutNames.size());
	S->cutFlowHisto[Elec] = new TH1D("EECutFlow", "EECutFlow", fEECutNames.size(), 0, fEECutNames.size());
	S->cutFlowHisto[ElMu] = new TH1D("EMCutFlow", "EMCutFlow", fEMCutNames.size(), 0, fEMCutNames.size());
	for(int i=0; i<fMMCutNames.size(); i++) S->cutFlowHisto[Muon]->GetXaxis()->SetBinLabel(i+1, TString(fMMCutNames[i]));
	for(int i=0; i<fEECutNames.size(); i++) S->cutFlowHisto[Elec]->GetXaxis()->SetBinLabel(i+1, TString(fEECutNames[i]));
	for(int i=0; i<fEMCutNames.size(); i++) S->cutFlowHisto[ElMu]->GetXaxis()->SetBinLabel(i+1, TString(fEMCutNames[i]));	
	
	// Histos for differential yields
	for(size_t k = 0; k < gNCHANNELS; ++k){
		TString name;
		for(size_t j = 0; j < gNDiffVars; ++j){
			name = Form("%s_%s_NT11_%s", S->sname.Data(), DiffPredYields::var_name[j].Data(), gChanLabel[k].Data());
			S->diffyields[k].hnt11[j] = new TH1D(name, DiffPredYields::var_name[j].Data(), DiffPredYields::nbins[j], DiffPredYields::bins[j]);
			S->diffyields[k].hnt11[j]->SetFillColor(S->color);
			S->diffyields[k].hnt11[j]->SetXTitle(DiffPredYields::axis_label[j]);
			S->diffyields[k].hnt11[j]->Sumw2();
			name = Form("%s_%s_NT10_%s", S->sname.Data(), DiffPredYields::var_name[j].Data(), gChanLabel[k].Data());
			S->diffyields[k].hnt10[j] = new TH1D(name, DiffPredYields::var_name[j].Data(), DiffPredYields::nbins[j], DiffPredYields::bins[j]);
			S->diffyields[k].hnt10[j]->SetFillColor(S->color);
			S->diffyields[k].hnt10[j]->SetXTitle(DiffPredYields::axis_label[j]);
			S->diffyields[k].hnt10[j]->Sumw2();
			name = Form("%s_%s_NT01_%s", S->sname.Data(), DiffPredYields::var_name[j].Data(), gChanLabel[k].Data());
			S->diffyields[k].hnt01[j] = new TH1D(name, DiffPredYields::var_name[j].Data(), DiffPredYields::nbins[j], DiffPredYields::bins[j]);
			S->diffyields[k].hnt01[j]->SetFillColor(S->color);
			S->diffyields[k].hnt01[j]->SetXTitle(DiffPredYields::axis_label[j]);
			S->diffyields[k].hnt01[j]->Sumw2();
			name = Form("%s_%s_NT00_%s", S->sname.Data(), DiffPredYields::var_name[j].Data(), gChanLabel[k].Data());
			S->diffyields[k].hnt00[j] = new TH1D(name, DiffPredYields::var_name[j].Data(), DiffPredYields::nbins[j], DiffPredYields::bins[j]);
			S->diffyields[k].hnt00[j]->SetFillColor(S->color);
			S->diffyields[k].hnt00[j]->SetXTitle(DiffPredYields::axis_label[j]);
			S->diffyields[k].hnt00[j]->Sumw2();

			if(k == Muon) continue;
			name = Form("%s_%s_NT11_OS_BB_%s", S->sname.Data(), DiffPredYields::var_name[j].Data(), gChanLabel[k].Data());
			S->diffyields[k].hnt2_os_BB[j] = new TH1D(name, DiffPredYields::var_name[j].Data(), DiffPredYields::nbins[j], DiffPredYields::bins[j]);
			S->diffyields[k].hnt2_os_BB[j]->SetFillColor(S->color);
			S->diffyields[k].hnt2_os_BB[j]->SetXTitle(DiffPredYields::axis_label[j]);
			S->diffyields[k].hnt2_os_BB[j]->Sumw2();
			name = Form("%s_%s_NT11_OS_EE_%s", S->sname.Data(), DiffPredYields::var_name[j].Data(), gChanLabel[k].Data());
			S->diffyields[k].hnt2_os_EE[j] = new TH1D(name, DiffPredYields::var_name[j].Data(), DiffPredYields::nbins[j], DiffPredYields::bins[j]);
			S->diffyields[k].hnt2_os_EE[j]->SetFillColor(S->color);
			S->diffyields[k].hnt2_os_EE[j]->SetXTitle(DiffPredYields::axis_label[j]);
			S->diffyields[k].hnt2_os_EE[j]->Sumw2();
			if(k == ElMu) continue;
			name = Form("%s_%s_NT11_OS_EB_%s", S->sname.Data(), DiffPredYields::var_name[j].Data(), gChanLabel[k].Data());
			S->diffyields[k].hnt2_os_EB[j] = new TH1D(name, DiffPredYields::var_name[j].Data(), DiffPredYields::nbins[j], DiffPredYields::bins[j]);
			S->diffyields[k].hnt2_os_EB[j]->SetFillColor(S->color);
			S->diffyields[k].hnt2_os_EB[j]->SetXTitle(DiffPredYields::axis_label[j]);
			S->diffyields[k].hnt2_os_EB[j]->Sumw2();
		}
	}	

	// Kinematical histos
	for(size_t hilo = 0; hilo < 2; ++hilo){
		for(size_t k = 0; k < gNKinSels; ++k){
			TString name = Form("%s_%s_%s_HTvsMET", S->sname.Data(), gKinSelNames[k].Data(), gHiLoLabel[hilo].Data());
			S->kinplots[k][hilo].hmetvsht = new TH2D(name, "HTvsMET", KinPlots::nHTBins, KinPlots::HTmin, KinPlots::HTmax, KinPlots::nMETBins, KinPlots::METmin, KinPlots::METmax);
			for(size_t j = 0; j < gNKinVars; ++j){
				name = Form("%s_%s_%s_%s", S->sname.Data(), gKinSelNames[k].Data(), gHiLoLabel[hilo].Data(), KinPlots::var_name[j].Data());
				S->kinplots[k][hilo].hvar[j] = new TH1D(name, KinPlots::var_name[j].Data(), KinPlots::nbins[j], KinPlots::xmin[j], KinPlots::xmax[j]);
				S->kinplots[k][hilo].hvar[j]->SetFillColor(S->color);
				S->kinplots[k][hilo].hvar[j]->SetXTitle(KinPlots::axis_label[j]);
				S->kinplots[k][hilo].hvar[j]->Sumw2();
			}
		}
	}

	for(size_t l = 0; l < 2; ++l){
		// Isolation histos
		for(size_t j = 0; j < gNSels; ++j){
			TString name = Form("%s_%s_%siso", S->sname.Data(), IsoPlots::sel_name[j].Data(), gEMULabel[l].Data());
			S->isoplots[l].hiso[j] = new TH1D(name, Form("%siso", gEMULabel[l].Data()), IsoPlots::nbins[j], 0., 1.);
			S->isoplots[l].hiso[j]->SetFillColor(S->color);
			S->isoplots[l].hiso[j]->SetXTitle("RelIso");
			S->isoplots[l].hiso[j]->Sumw2();
			for(int k = 0; k < gNMuPt2bins; ++k){
				name = Form("%s_%s_%siso_pt%d", S->sname.Data(), IsoPlots::sel_name[j].Data(), gEMULabel[l].Data(), k);
				S->isoplots[l].hiso_pt[j][k] = new TH1D(name, Form("%siso_pt%d", gEMULabel[l].Data(), k), IsoPlots::nbins[j], 0., 1.);
				S->isoplots[l].hiso_pt[j][k]->SetFillColor(S->color);
				S->isoplots[l].hiso_pt[j][k]->SetXTitle("RelIso");
				S->isoplots[l].hiso_pt[j][k]->Sumw2();
			}
			for(int k = 0; k < gNNVrtxBins; ++k){
				name = Form("%s_%s_%siso_nv%d", S->sname.Data(), IsoPlots::sel_name[j].Data(), gEMULabel[l].Data(), k);
				S->isoplots[l].hiso_nv[j][k] = new TH1D(name, Form("%siso_nv%d", gEMULabel[l].Data(), k), IsoPlots::nbins[j], 0., 1.);
				S->isoplots[l].hiso_nv[j][k]->SetFillColor(S->color);
				S->isoplots[l].hiso_nv[j][k]->SetXTitle("RelIso");
				S->isoplots[l].hiso_nv[j][k]->Sumw2();
			}
		}

		// Ratio histos
		for(size_t j = 0; j < gNRatioVars; ++j){
			TString name = Form("%s_%s_ntight_%s", S->sname.Data(), gEMULabel[l].Data(), FRatioPlots::var_name[j].Data());
			S->ratioplots[l].ntight[j] = new TH1D(name, "ntight", FRatioPlots::nbins[j], FRatioPlots::xmin[j], FRatioPlots::xmax[j]);
			S->ratioplots[l].ntight[j]->SetFillColor(S->color);
			S->ratioplots[l].ntight[j]->SetXTitle(FRatioPlots::var_name[j]);
			S->ratioplots[l].ntight[j]->Sumw2();
			name = Form("%s_%s_nloose_%s", S->sname.Data(), gEMULabel[l].Data(), FRatioPlots::var_name[j].Data());
			S->ratioplots[l].nloose[j] = new TH1D(name, "nloose", FRatioPlots::nbins[j], FRatioPlots::xmin[j], FRatioPlots::xmax[j]);
			S->ratioplots[l].nloose[j]->SetFillColor(S->color);
			S->ratioplots[l].nloose[j]->SetXTitle(FRatioPlots::var_name[j]);
			S->ratioplots[l].nloose[j]->Sumw2();
		}
	}

	for(size_t hilo = 0; hilo < 2; ++hilo){
		for(gRegion r = region_begin; r < gNREGIONS; r=gRegion(r+1)){
			Region *R = &S->region[r][hilo];
			for(gChannel c = channels_begin; c < gNCHANNELS; c=gChannel(c+1)){
				Channel *C;
				if(c == Muon) C = &R->mm;
				if(c == Elec) C = &R->ee;
				if(c == ElMu) C = &R->em;
				TString rootname = S->sname + "_" + Region::sname[r] + "_" + gChanLabel[c] + "_" + gHiLoLabel[hilo];
				// Yields common for all channels and data-mc:
				C->nt20_pt  = new TH2D(rootname + "_NT20_pt",  "NT20_pt",  getNPt2Bins(c), getPt2Bins(c), getNPt2Bins(c), getPt2Bins(c)); C->nt20_pt ->Sumw2();
				C->nt10_pt  = new TH2D(rootname + "_NT10_pt",  "NT10_pt",  getNPt2Bins(c), getPt2Bins(c), getNPt2Bins(c), getPt2Bins(c)); C->nt10_pt ->Sumw2();
				C->nt01_pt  = new TH2D(rootname + "_NT01_pt",  "NT01_pt",  getNPt2Bins(c), getPt2Bins(c), getNPt2Bins(c), getPt2Bins(c)); C->nt01_pt ->Sumw2();
				C->nt00_pt  = new TH2D(rootname + "_NT00_pt",  "NT00_pt",  getNPt2Bins(c), getPt2Bins(c), getNPt2Bins(c), getPt2Bins(c)); C->nt00_pt ->Sumw2();
				C->nt20_eta = new TH2D(rootname + "_NT20_eta", "NT20_eta", getNEtaBins(c), getEtaBins(c), getNEtaBins(c), getEtaBins(c)); C->nt20_eta->Sumw2();
				C->nt10_eta = new TH2D(rootname + "_NT10_eta", "NT10_eta", getNEtaBins(c), getEtaBins(c), getNEtaBins(c), getEtaBins(c)); C->nt10_eta->Sumw2();
				C->nt01_eta = new TH2D(rootname + "_NT01_eta", "NT01_eta", getNEtaBins(c), getEtaBins(c), getNEtaBins(c), getEtaBins(c)); C->nt01_eta->Sumw2();
				C->nt00_eta = new TH2D(rootname + "_NT00_eta", "NT00_eta", getNEtaBins(c), getEtaBins(c), getNEtaBins(c), getEtaBins(c)); C->nt00_eta->Sumw2();

				// MC truth info
				if(S->datamc > 0){
					C->npp_pt   = new TH2D(rootname + "_NPP_pt",   "NPP_pt",   getNPt2Bins(c), getPt2Bins(c), getNPt2Bins(c), getPt2Bins(c)); C->npp_pt->Sumw2();
					C->nfp_pt   = new TH2D(rootname + "_NFP_pt",   "NFP_pt",   getNPt2Bins(c), getPt2Bins(c), getNPt2Bins(c), getPt2Bins(c)); C->nfp_pt->Sumw2();
					C->npf_pt   = new TH2D(rootname + "_NPF_pt",   "NPF_pt",   getNPt2Bins(c), getPt2Bins(c), getNPt2Bins(c), getPt2Bins(c)); C->npf_pt->Sumw2();
					C->nff_pt   = new TH2D(rootname + "_NFF_pt",   "NFF_pt",   getNPt2Bins(c), getPt2Bins(c), getNPt2Bins(c), getPt2Bins(c)); C->nff_pt->Sumw2();
					C->nt2pp_pt = new TH2D(rootname + "_NT2PP_pt", "NT2PP_pt", getNPt2Bins(c), getPt2Bins(c), getNPt2Bins(c), getPt2Bins(c)); C->nt2pp_pt->Sumw2();
					C->nt2fp_pt = new TH2D(rootname + "_NT2FP_pt", "NT2FP_pt", getNPt2Bins(c), getPt2Bins(c), getNPt2Bins(c), getPt2Bins(c)); C->nt2fp_pt->Sumw2();
					C->nt2pf_pt = new TH2D(rootname + "_NT2PF_pt", "NT2PF_pt", getNPt2Bins(c), getPt2Bins(c), getNPt2Bins(c), getPt2Bins(c)); C->nt2pf_pt->Sumw2();
					C->nt2ff_pt = new TH2D(rootname + "_NT2FF_pt", "NT2FF_pt", getNPt2Bins(c), getPt2Bins(c), getNPt2Bins(c), getPt2Bins(c)); C->nt2ff_pt->Sumw2();

					C->nt11_origin = new TH2D(rootname + "_NT20_Origin",  "NT2Origin",  15, 0, 15, 15, 0, 15);
					C->nt10_origin = new TH2D(rootname + "_NT10_Origin",  "NT1Origin",  15, 0, 15, 15, 0, 15);
					C->nt01_origin = new TH2D(rootname + "_NT01_Origin",  "NT01Origin", 15, 0, 15, 15, 0, 15);
					C->nt00_origin = new TH2D(rootname + "_NT00_Origin",  "NT0Origin",  15, 0, 15, 15, 0, 15);
					label2OriginAxes(C->nt11_origin->GetXaxis(), C->nt11_origin->GetYaxis(), c);
					label2OriginAxes(C->nt10_origin->GetXaxis(), C->nt10_origin->GetYaxis(), c);
					label2OriginAxes(C->nt01_origin->GetXaxis(), C->nt01_origin->GetYaxis(), c);
					label2OriginAxes(C->nt00_origin->GetXaxis(), C->nt00_origin->GetYaxis(), c);					
				}

				// Charge misid truth
				if(c != Muon){
					C->npp_cm_pt   = new TH2D(rootname + "_NPP_CM_pt",   "NPP_CM_pt",   getNPt2Bins(c), getPt2Bins(c), getNPt2Bins(c), getPt2Bins(c)); C->npp_cm_pt->Sumw2();
					C->nt2pp_cm_pt = new TH2D(rootname + "_NT2PP_CM_pt", "NT2PP_CM_pt", getNPt2Bins(c), getPt2Bins(c), getNPt2Bins(c), getPt2Bins(c)); C->nt2pp_cm_pt->Sumw2();					
				}

					// OS Yields
				if(c == Elec){
					C->nt20_OS_BB_pt = new TH2D(rootname + "_NT20_OS_BB_pt",  "NT20_OS_BB_pt",  getNPt2Bins(c), getPt2Bins(c), getNPt2Bins(c), getPt2Bins(c)); C->nt20_OS_BB_pt ->Sumw2();
	                   C->nt20_OS_EE_pt = new TH2D(rootname + "_NT20_OS_EE_pt",  "NT20_OS_EE_pt",  getNPt2Bins(c), getPt2Bins(c), getNPt2Bins(c), getPt2Bins(c)); C->nt20_OS_EE_pt ->Sumw2();
	                   C->nt20_OS_EB_pt = new TH2D(rootname + "_NT20_OS_EB_pt",  "NT20_OS_EB_pt",  getNPt2Bins(c), getPt2Bins(c), getNPt2Bins(c), getPt2Bins(c)); C->nt20_OS_EB_pt ->Sumw2();
				}
				if(c == ElMu){
					C->nt20_OS_BB_pt = new TH2D(rootname + "_NT20_OS_BB_pt",  "NT20_OS_BB_pt",  getNPt2Bins(Muon), getPt2Bins(Muon), getNPt2Bins(Elec), getPt2Bins(Elec)); C->nt20_OS_BB_pt ->Sumw2();
	                   C->nt20_OS_EE_pt = new TH2D(rootname + "_NT20_OS_EE_pt",  "NT20_OS_EE_pt",  getNPt2Bins(Muon), getPt2Bins(Muon), getNPt2Bins(Elec), getPt2Bins(Elec)); C->nt20_OS_EE_pt ->Sumw2();
				}

				// Ratios
				if(c != ElMu){
					C->fntight  = new TH2D(rootname + "_fNTight",  "fNTight",  getNPt2Bins(c), getPt2Bins(c), getNEtaBins(c), getEtaBins(c)); C->fntight ->Sumw2();
					C->fnloose  = new TH2D(rootname + "_fNLoose",  "fNLoose",  getNPt2Bins(c), getPt2Bins(c), getNEtaBins(c), getEtaBins(c)); C->fnloose ->Sumw2();
					C->pntight  = new TH2D(rootname + "_pNTight",  "pNTight",  getNPt2Bins(c), getPt2Bins(c), getNEtaBins(c), getEtaBins(c)); C->pntight ->Sumw2();
					C->pnloose  = new TH2D(rootname + "_pNLoose",  "pNLoose",  getNPt2Bins(c), getPt2Bins(c), getNEtaBins(c), getEtaBins(c)); C->pnloose ->Sumw2();

					if(S->datamc > 0){
						C->sst_origin = new TH1D(rootname + "_fTOrigin", "fTOrigin", 15, 0, 15); C->sst_origin->Sumw2();
						C->ssl_origin = new TH1D(rootname + "_fLOrigin", "fLOrigin", 15, 0, 15); C->ssl_origin->Sumw2();
						C->zt_origin  = new TH1D(rootname + "_pTOrigin", "pTOrigin", 15, 0, 15); C->zt_origin ->Sumw2();
						C->zl_origin  = new TH1D(rootname + "_pLOrigin", "pLOrigin", 15, 0, 15); C->zl_origin ->Sumw2();
						labelOriginAxis(C->sst_origin->GetXaxis() , c);
						labelOriginAxis(C->ssl_origin->GetXaxis() , c);
						labelOriginAxis(C->zt_origin->GetXaxis()  , c);
						labelOriginAxis(C->zl_origin->GetXaxis()  , c);
					}
				}
			}
		}
	}
}
void SSDLDumper::deleteHistos(Sample *S){
	delete S->cutFlowHisto[Muon];
	delete S->cutFlowHisto[Elec];
	delete S->cutFlowHisto[ElMu];
	
	// Kinematical histos
	for(size_t hilo = 0; hilo < 2; ++hilo){
		for(size_t k = 0; k < gNKinSels; ++k){
			delete S->kinplots[k][hilo].hmetvsht;
			for(size_t j = 0; j < gNKinVars; ++j) delete S->kinplots[k][hilo].hvar[j];
		}
	}

	// Histos for differential yields
	for(size_t k = 0; k < gNCHANNELS; ++k){
		for(size_t j = 0; j < gNDiffVars; ++j){
			delete S->diffyields[k].hnt11[j];
			delete S->diffyields[k].hnt10[j];
			delete S->diffyields[k].hnt01[j];
			delete S->diffyields[k].hnt00[j];
			if(k == Muon) continue;
			delete S->diffyields[k].hnt2_os_BB[j];
			delete S->diffyields[k].hnt2_os_EE[j];
			if(k == ElMu) continue;
			delete S->diffyields[k].hnt2_os_EB[j];
		}
	}

	for(size_t l = 0; l < 2; ++l){
		// Isolation histos
		for(size_t j = 0; j < gNSels; ++j){
			delete S->isoplots[l].hiso[j];
			for(int k = 0; k < gNMuPt2bins; ++k){
				delete S->isoplots[l].hiso_pt[j][k];
			}
			for(int k = 0; k < gNNVrtxBins; ++k){
				delete S->isoplots[l].hiso_nv[j][k];
			}
		}

		// Ratio histos
		for(size_t j = 0; j < gNRatioVars; ++j){
			delete S->ratioplots[l].ntight[j];
			delete S->ratioplots[l].nloose[j];
		}
	}

	for(size_t hilo = 0; hilo < 2; ++hilo){
		for(gRegion r = region_begin; r < gNREGIONS; r=gRegion(r+1)){
			Region *R = &S->region[r][hilo];
			for(gChannel c = channels_begin; c < gNCHANNELS; c=gChannel(c+1)){
				Channel *C;
				if(c == Muon) C = &R->mm;
				if(c == Elec) C = &R->ee;
				if(c == ElMu) C = &R->em;

				delete C->nt20_pt;
				delete C->nt10_pt;
				delete C->nt01_pt;
				delete C->nt00_pt;
				delete C->nt20_eta;
				delete C->nt10_eta;
				delete C->nt01_eta;
				delete C->nt00_eta;

				// MC truth info
				if(S->datamc > 0){
					delete C->npp_pt;
					delete C->nfp_pt;
					delete C->npf_pt;
					delete C->nff_pt;
					delete C->nt2pp_pt;
					delete C->nt2fp_pt;
					delete C->nt2pf_pt;
					delete C->nt2ff_pt;

					delete C->nt11_origin;
					delete C->nt10_origin;
					delete C->nt01_origin;
					delete C->nt00_origin;
				}

				// Charge misid truth
				if(c != Muon){
					delete C->npp_cm_pt;
					delete C->nt2pp_cm_pt;
				}

					// OS Yields
				if(c == Elec){
					delete C->nt20_OS_BB_pt;
					delete C->nt20_OS_EE_pt;
					delete C->nt20_OS_EB_pt;
				}
				if(c == ElMu){
					delete C->nt20_OS_BB_pt;
					delete C->nt20_OS_EE_pt;
				}

				// Ratios
				if(c != ElMu){
					delete C->fntight;
					delete C->fnloose;
					delete C->pntight;
					delete C->pnloose;

					if(S->datamc > 0){
						delete C->sst_origin;
						delete C->ssl_origin;
						delete C->zt_origin;
						delete C->zl_origin;
					}
				}
			}
		}
	}
}
void SSDLDumper::writeHistos(Sample *S, TFile *pFile){
	pFile->cd();
	TDirectory* sdir = Util::FindOrCreate(S->sname, pFile);
	sdir->cd();

	TString temp;
	TDirectory *rdir;

	// Cut Flows
	S->cutFlowHisto[Muon]->Write(S->cutFlowHisto[Muon]->GetName(), TObject::kWriteDelete);
	S->cutFlowHisto[Elec]->Write(S->cutFlowHisto[Elec]->GetName(), TObject::kWriteDelete);
	S->cutFlowHisto[ElMu]->Write(S->cutFlowHisto[ElMu]->GetName(), TObject::kWriteDelete);

	// Histos for differential yields
	for(size_t k = 0; k < gNCHANNELS; ++k){
		temp = S->sname + "/DiffYields/";
		rdir = Util::FindOrCreate(temp, pFile);
		rdir->cd();
		for(size_t j = 0; j < gNDiffVars; ++j){
			S->diffyields[k].hnt11[j]->Write(S->diffyields[k].hnt11[j]->GetName(), TObject::kWriteDelete);
			S->diffyields[k].hnt10[j]->Write(S->diffyields[k].hnt10[j]->GetName(), TObject::kWriteDelete);
			S->diffyields[k].hnt01[j]->Write(S->diffyields[k].hnt01[j]->GetName(), TObject::kWriteDelete);
			S->diffyields[k].hnt00[j]->Write(S->diffyields[k].hnt00[j]->GetName(), TObject::kWriteDelete);
			if(k == Muon) continue;
			S->diffyields[k].hnt2_os_BB[j]->Write(S->diffyields[k].hnt2_os_BB[j]->GetName(), TObject::kWriteDelete);
			S->diffyields[k].hnt2_os_EE[j]->Write(S->diffyields[k].hnt2_os_EE[j]->GetName(), TObject::kWriteDelete);
			if(k == ElMu) continue;
			S->diffyields[k].hnt2_os_EB[j]->Write(S->diffyields[k].hnt2_os_EB[j]->GetName(), TObject::kWriteDelete);
		}
	}	

	// Kinematic histos
	for(size_t hilo = 0; hilo < 2; ++hilo){
		temp = S->sname + "/KinPlots/" + gHiLoLabel[hilo] + "/";
		rdir = Util::FindOrCreate(temp, pFile);
		rdir->cd();
		for(size_t k = 0; k < gNKinSels; ++k){
			KinPlots *kp = &S->kinplots[k][hilo];
			kp->hmetvsht->Write(kp->hmetvsht->GetName(), TObject::kWriteDelete);
			for(size_t j = 0; j < gNKinVars; ++j) kp->hvar[j]->Write(kp->hvar[j]->GetName(), TObject::kWriteDelete);
		}
	}

	// Isolation histos
	temp = S->sname + "/IsoPlots/";
	rdir = Util::FindOrCreate(temp, pFile);
	rdir->cd();
	for(size_t l = 0; l < 2; ++l){
		IsoPlots *ip = &S->isoplots[l];
		for(size_t j = 0; j < gNSels; ++j){
			ip->hiso[j]->Write(ip->hiso[j]->GetName(), TObject::kWriteDelete);
			for(int k = 0; k < gNMuPt2bins; ++k) ip->hiso_pt[j][k]->Write(ip->hiso_pt[j][k]->GetName(), TObject::kWriteDelete);
			for(int k = 0; k < gNNVrtxBins; ++k) ip->hiso_nv[j][k]->Write(ip->hiso_nv[j][k]->GetName(), TObject::kWriteDelete);
		}
	}

	// Ratio histos
	temp = S->sname + "/FRatioPlots/";
	rdir = Util::FindOrCreate(temp, pFile);
	rdir->cd();
	for(size_t l = 0; l < 2; ++l){
		FRatioPlots *rp = &S->ratioplots[l];
		for(size_t j = 0; j < gNRatioVars; ++j){
			rp->ntight[j]->Write(rp->ntight[j]->GetName(), TObject::kWriteDelete);
			rp->nloose[j]->Write(rp->nloose[j]->GetName(), TObject::kWriteDelete);
		}
	}

	// Yields
	for(size_t hilo = 0; hilo < 2; ++hilo){
		for(gRegion r = region_begin; r < gNREGIONS; r=gRegion(r+1)){
			Region *R = &S->region[r][hilo];
			TString temp = S->sname + "/" + Region::sname[r] + "/" + gHiLoLabel[hilo] + "/";
			TDirectory* rdir = Util::FindOrCreate(temp, pFile);
			rdir->cd();

			for(gChannel ch = channels_begin; ch < gNCHANNELS; ch=gChannel(ch+1)){ // Loop over channels, mumu, emu, ee
				Channel *C;
				if(ch == Muon)     C = &R->mm;
				if(ch == Elec) C = &R->ee;
				if(ch == ElMu)      C = &R->em;
				C->nt20_pt    ->Write(C->nt20_pt    ->GetName(), TObject::kWriteDelete);
				C->nt10_pt    ->Write(C->nt10_pt    ->GetName(), TObject::kWriteDelete);
				C->nt01_pt    ->Write(C->nt01_pt    ->GetName(), TObject::kWriteDelete);
				C->nt00_pt    ->Write(C->nt00_pt    ->GetName(), TObject::kWriteDelete);
				C->nt20_eta   ->Write(C->nt20_eta   ->GetName(), TObject::kWriteDelete);
				C->nt10_eta   ->Write(C->nt10_eta   ->GetName(), TObject::kWriteDelete);
				C->nt01_eta   ->Write(C->nt01_eta   ->GetName(), TObject::kWriteDelete);
				C->nt00_eta   ->Write(C->nt00_eta   ->GetName(), TObject::kWriteDelete);
			
				if(ch == Elec || ch == ElMu){
					C->nt20_OS_BB_pt->Write(C->nt20_OS_BB_pt->GetName(), TObject::kWriteDelete);
					C->nt20_OS_EE_pt->Write(C->nt20_OS_EE_pt->GetName(), TObject::kWriteDelete);
					if(ch == Elec) C->nt20_OS_EB_pt->Write(C->nt20_OS_EB_pt->GetName(), TObject::kWriteDelete);					
				}
			
				if(S->datamc > 0){
					C->npp_pt     ->Write(C->npp_pt     ->GetName(), TObject::kWriteDelete);
					C->nfp_pt     ->Write(C->nfp_pt     ->GetName(), TObject::kWriteDelete);
					C->npf_pt     ->Write(C->npf_pt     ->GetName(), TObject::kWriteDelete);
					C->nff_pt     ->Write(C->nff_pt     ->GetName(), TObject::kWriteDelete);
					C->nt2pp_pt   ->Write(C->nt2pp_pt   ->GetName(), TObject::kWriteDelete);
					C->nt2fp_pt   ->Write(C->nt2fp_pt   ->GetName(), TObject::kWriteDelete);
					C->nt2pf_pt   ->Write(C->nt2pf_pt   ->GetName(), TObject::kWriteDelete);
					C->nt2ff_pt   ->Write(C->nt2ff_pt   ->GetName(), TObject::kWriteDelete);
					C->nt11_origin->Write(C->nt11_origin->GetName(), TObject::kWriteDelete);
					C->nt10_origin->Write(C->nt10_origin->GetName(), TObject::kWriteDelete);
					C->nt01_origin->Write(C->nt01_origin->GetName(), TObject::kWriteDelete);
					C->nt00_origin->Write(C->nt00_origin->GetName(), TObject::kWriteDelete);
					if(ch != Muon){
						C->npp_cm_pt  ->Write(C->npp_cm_pt  ->GetName(), TObject::kWriteDelete);
						C->nt2pp_cm_pt->Write(C->nt2pp_cm_pt->GetName(), TObject::kWriteDelete);						
					}
				}
				if(ch != ElMu){
					C->fntight    ->Write(C->fntight    ->GetName(), TObject::kWriteDelete);
					C->fnloose    ->Write(C->fnloose    ->GetName(), TObject::kWriteDelete);
					C->pntight    ->Write(C->pntight    ->GetName(), TObject::kWriteDelete);
					C->pnloose    ->Write(C->pnloose    ->GetName(), TObject::kWriteDelete);
					if(S->datamc > 0){
						C->sst_origin ->Write(C->sst_origin ->GetName(), TObject::kWriteDelete);
						C->ssl_origin ->Write(C->ssl_origin ->GetName(), TObject::kWriteDelete);
						C->zt_origin  ->Write(C->zt_origin  ->GetName(), TObject::kWriteDelete);
						C->zl_origin  ->Write(C->zl_origin  ->GetName(), TObject::kWriteDelete);						
					}
				}
			}
		}
	}
}
void SSDLDumper::writeSigGraphs(Sample *S, gChannel chan, gHiLoSwitch hilo, TFile *pFile){
	TString channame = "MM";
	if(chan == Elec) channame = "EE";
	if(chan == ElMu)      channame = "EM";
	vector<float> ht;
	vector<float> met;
	if(chan == Muon){
		if(hilo == HighPt){
			ht  = fSigEv_HI_MM_HT;
			met = fSigEv_HI_MM_MET;
		}
		if(hilo == LowPt){
			ht  = fSigEv_LO_MM_HT;
			met = fSigEv_LO_MM_MET;
		}
	}
	if(chan == Elec){
		if(hilo == HighPt){
			ht  = fSigEv_HI_EE_HT;
			met = fSigEv_HI_EE_MET;
		}
		if(hilo == LowPt){
			ht  = fSigEv_LO_EE_HT;
			met = fSigEv_LO_EE_MET;
		}
	}
	if(chan == ElMu){
		if(hilo == HighPt){
			ht  = fSigEv_HI_EM_HT;
			met = fSigEv_HI_EM_MET;
		}
		if(hilo == LowPt){
			ht  = fSigEv_LO_EM_HT;
			met = fSigEv_LO_EM_MET;
		}
	}
	
	const int nsig = ht.size();
	float ht_a [nsig];
	float met_a[nsig];
	for(size_t i = 0; i < ht.size(); ++i){
		ht_a[i] = ht[i];
		met_a[i] = met[i];
	}
	
	TGraph *sigevents = new TGraph(nsig, ht_a, met_a);
	sigevents->SetName(Form("%s_%s_%s_SigEvents", S->sname.Data(), channame.Data(), gHiLoLabel[hilo].Data()));

	pFile->cd();
	TString dirname = S->sname + "/SigGraphs/";
	TDirectory* dir = Util::FindOrCreate(dirname, pFile);
	dir->cd();
	sigevents->Write(sigevents->GetName(), TObject::kWriteDelete);
	delete sigevents;
}
int  SSDLDumper::readHistos(TString filename){
	TFile *pFile = TFile::Open(filename, "READ");
	if(pFile == NULL){
		cout << "File " << filename << " does not exist!" << endl;
		exit(1);
	}

	pFile->cd();
	if(gNSAMPLES != fSamples.size()){
		cout << "Mismatch in number of samples! Help!" << endl;
		exit(1);
	}

	TString getname;
	
	for(gSample i = sample_begin; i < gNSAMPLES; i=gSample(i+1)){
		Sample *S = fSamples[i];

		// Cut flow histos
		S->cutFlowHisto[Muon] = (TH1D*)pFile->Get(S->sname + "/MMCutFlow");
		S->cutFlowHisto[Elec] = (TH1D*)pFile->Get(S->sname + "/EECutFlow");
		S->cutFlowHisto[ElMu] = (TH1D*)pFile->Get(S->sname + "/EMCutFlow");

		// Histos for differential yields
		for(size_t k = 0; k < gNCHANNELS; ++k){
			TString name;
			for(size_t j = 0; j < gNDiffVars; ++j){
				getname = Form("%s_%s_NT11_%s", S->sname.Data(), DiffPredYields::var_name[j].Data(), gChanLabel[k].Data());
				S->diffyields[k].hnt11[j] = (TH1D*)pFile->Get(S->sname + "/DiffYields/" + getname);
				getname = Form("%s_%s_NT10_%s", S->sname.Data(), DiffPredYields::var_name[j].Data(), gChanLabel[k].Data());
				S->diffyields[k].hnt10[j] = (TH1D*)pFile->Get(S->sname + "/DiffYields/" + getname);
				getname = Form("%s_%s_NT01_%s", S->sname.Data(), DiffPredYields::var_name[j].Data(), gChanLabel[k].Data());
				S->diffyields[k].hnt01[j] = (TH1D*)pFile->Get(S->sname + "/DiffYields/" + getname);
				getname = Form("%s_%s_NT00_%s", S->sname.Data(), DiffPredYields::var_name[j].Data(), gChanLabel[k].Data());
				S->diffyields[k].hnt00[j] = (TH1D*)pFile->Get(S->sname + "/DiffYields/" + getname);
				if(k == Muon) continue;
				getname = Form("%s_%s_NT11_OS_BB_%s", S->sname.Data(), DiffPredYields::var_name[j].Data(), gChanLabel[k].Data());
				S->diffyields[k].hnt2_os_BB[j] = (TH1D*)pFile->Get(S->sname + "/DiffYields/" + getname);
				getname = Form("%s_%s_NT11_OS_EE_%s", S->sname.Data(), DiffPredYields::var_name[j].Data(), gChanLabel[k].Data());
				S->diffyields[k].hnt2_os_EE[j] = (TH1D*)pFile->Get(S->sname + "/DiffYields/" + getname);
				if(k == ElMu) continue;
				getname = Form("%s_%s_NT11_OS_EB_%s", S->sname.Data(), DiffPredYields::var_name[j].Data(), gChanLabel[k].Data());
				S->diffyields[k].hnt2_os_EB[j] = (TH1D*)pFile->Get(S->sname + "/DiffYields/" + getname);
			}
		}	


		// Kinematic histos
		for(size_t hilo = 0; hilo < 2; ++hilo){
			for(size_t k = 0; k < gNKinSels; ++k){
				KinPlots *kp = &S->kinplots[k][hilo];
				getname = Form("%s_%s_%s_HTvsMET", S->sname.Data(), gKinSelNames[k].Data(), gHiLoLabel[hilo].Data());
				kp->hmetvsht = (TH2D*)pFile->Get(S->sname + "/KinPlots/" + gHiLoLabel[hilo] + "/" + getname);
				for(size_t j = 0; j < gNKinVars; ++j){
					getname = Form("%s_%s_%s_%s", S->sname.Data(), gKinSelNames[k].Data(), gHiLoLabel[hilo].Data(), KinPlots::var_name[j].Data());
					kp->hvar[j] = (TH1D*)pFile->Get(S->sname + "/KinPlots/" + gHiLoLabel[hilo] + "/" + getname);
					kp->hvar[j]->SetFillColor(S->color);
				}
			}
		}

		for(size_t lep = 0; lep < 2; ++lep){ // e-mu loop
			// Isolation histos
			IsoPlots *ip = &S->isoplots[lep];
			for(size_t j = 0; j < gNSels; ++j){
				getname = Form("%s_%s_%siso", S->sname.Data(), IsoPlots::sel_name[j].Data(), gEMULabel[lep].Data());
				ip->hiso[j] = (TH1D*)pFile->Get(S->sname + "/IsoPlots/" + getname);
				ip->hiso[j]->SetFillColor(S->color);
				for(int k = 0; k < gNMuPt2bins; ++k){
					getname = Form("%s_%s_%siso_pt%d", S->sname.Data(), IsoPlots::sel_name[j].Data(), gEMULabel[lep].Data(), k);
					ip->hiso_pt[j][k] = (TH1D*)pFile->Get(S->sname + "/IsoPlots/" + getname);
					ip->hiso_pt[j][k]->SetFillColor(S->color);
				}
				for(int k = 0; k < gNNVrtxBins; ++k){
					getname = Form("%s_%s_%siso_nv%d", S->sname.Data(), IsoPlots::sel_name[j].Data(), gEMULabel[lep].Data(), k);
					ip->hiso_nv[j][k] = (TH1D*)pFile->Get(S->sname + "/IsoPlots/" + getname);
					ip->hiso_nv[j][k]->SetFillColor(S->color);
				}
			}

			// Ratio histos
			FRatioPlots *rp = &S->ratioplots[lep];
			for(size_t j = 0; j < gNRatioVars; ++j){
				getname = Form("%s_%s_ntight_%s", S->sname.Data(), gEMULabel[lep].Data(), FRatioPlots::var_name[j].Data());
				rp->ntight[j] = (TH1D*)pFile->Get(S->sname + "/FRatioPlots/" + getname);
				getname = Form("%s_%s_nloose_%s", S->sname.Data(), gEMULabel[lep].Data(), FRatioPlots::var_name[j].Data());
				rp->nloose[j] = (TH1D*)pFile->Get(S->sname + "/FRatioPlots/" + getname);
			}
		}

		// Yields
		for(size_t hilo = 0; hilo < 2; ++hilo){
			for(gRegion r = region_begin; r < gNREGIONS; r=gRegion(r+1)){ // Loop over regions
				Region *R = &S->region[r][hilo];
				for(gChannel ch = channels_begin; ch < gNCHANNELS; ch=gChannel(ch+1)){ // Loop over channels, mumu, emu, ee
					Channel *C;
					if(ch == Muon)     C = &R->mm;
					if(ch == Elec) C = &R->ee;
					if(ch == ElMu)      C = &R->em;
					TString root = S->sname +"/"+ Region::sname[r] +"/"+ gHiLoLabel[hilo] +"/"+ S->sname +"_"+ Region::sname[r] +"_"+ gChanLabel[ch] +"_"+ gHiLoLabel[hilo];
					C->nt20_pt  = (TH2D*)pFile->Get(root + "_NT20_pt");
					C->nt10_pt  = (TH2D*)pFile->Get(root + "_NT10_pt");
					C->nt01_pt  = (TH2D*)pFile->Get(root + "_NT01_pt");
					C->nt00_pt  = (TH2D*)pFile->Get(root + "_NT00_pt");
					C->nt20_eta = (TH2D*)pFile->Get(root + "_NT20_eta");
					C->nt10_eta = (TH2D*)pFile->Get(root + "_NT10_eta");
					C->nt01_eta = (TH2D*)pFile->Get(root + "_NT01_eta");
					C->nt00_eta = (TH2D*)pFile->Get(root + "_NT00_eta");

					if(S->datamc > 0){
						C->npp_pt      = (TH2D*)pFile->Get(root + "_NPP_pt");
						C->nfp_pt      = (TH2D*)pFile->Get(root + "_NFP_pt");
						C->npf_pt      = (TH2D*)pFile->Get(root + "_NPF_pt");
						C->nff_pt      = (TH2D*)pFile->Get(root + "_NFF_pt");
						C->nt2pp_pt    = (TH2D*)pFile->Get(root + "_NT2PP_pt");
						C->nt2fp_pt    = (TH2D*)pFile->Get(root + "_NT2FP_pt");
						C->nt2pf_pt    = (TH2D*)pFile->Get(root + "_NT2PF_pt");
						C->nt2ff_pt    = (TH2D*)pFile->Get(root + "_NT2FF_pt");
						C->nt11_origin = (TH2D*)pFile->Get(root + "_NT20_Origin");
						C->nt10_origin = (TH2D*)pFile->Get(root + "_NT10_Origin");
						C->nt01_origin = (TH2D*)pFile->Get(root + "_NT01_Origin");
						C->nt00_origin = (TH2D*)pFile->Get(root + "_NT00_Origin");

						if(ch != Muon){
							C->npp_cm_pt   = (TH2D*)pFile->Get(root + "_NPP_CM_pt");
							C->nt2pp_cm_pt = (TH2D*)pFile->Get(root + "_NT2PP_CM_pt");						
						}
					}
					if(ch == Elec || ch == ElMu){
						C->nt20_OS_BB_pt = (TH2D*)pFile->Get(root + "_NT20_OS_BB_pt");
						C->nt20_OS_EE_pt = (TH2D*)pFile->Get(root + "_NT20_OS_EE_pt");
						if(ch == Elec) C->nt20_OS_EB_pt = (TH2D*)pFile->Get(root + "_NT20_OS_EB_pt");
					}

					if(ch != ElMu){
						C->fntight     = (TH2D*)pFile->Get(root + "_fNTight");
						C->fnloose     = (TH2D*)pFile->Get(root + "_fNLoose");
						C->pntight     = (TH2D*)pFile->Get(root + "_pNTight");
						C->pnloose     = (TH2D*)pFile->Get(root + "_pNLoose");
						if(S->datamc > 0){
							C->sst_origin  = (TH1D*)pFile->Get(root + "_fTOrigin");
							C->ssl_origin  = (TH1D*)pFile->Get(root + "_fLOrigin");
							C->zt_origin   = (TH1D*)pFile->Get(root + "_pTOrigin");
							C->zl_origin   = (TH1D*)pFile->Get(root + "_pLOrigin");
						}
					}
				}
				if(hilo == HighPt){
					storeNumbers(S, Muon, r);
					storeNumbers(S, Elec, r);
					storeNumbers(S, ElMu, r);					
				}
			}
		}
	}
	return 0;
}
int  SSDLDumper::readSigGraphs(TString filename){
	TFile *pFile = TFile::Open(filename, "READ");
	if(pFile == NULL){
		cout << "File " << filename << " does not exist!" << endl;
		return 1;
	}
	
	TString channame, getname;
	Color_t color[3] = {kBlack, kBlue, kRed};
	// Size_t  size [3] = {1.6, 1.8, 1.5};
	Size_t size = 1.5;
	Style_t style[3] = {8, 23, 21};
	
	for(gSample i = sample_begin; i < gNSAMPLES; i=gSample(i+1)){
		Sample *S = fSamples[i];
		for(size_t hilo = 0; hilo < 2; ++hilo){
			for(gChannel ch = channels_begin; ch < gNCHANNELS; ch=gChannel(ch+1)){ // Loop over channels, mumu, emu, ee
				if(ch == Muon)     channame = "MM";
				if(ch == Elec) channame = "EE";
				if(ch == ElMu)      channame = "EM";
				getname = Form("%s/SigGraphs/%s_%s_%s_SigEvents", S->sname.Data(), S->sname.Data(), channame.Data(), gHiLoLabel[hilo].Data());
				S->sigevents[ch][hilo] = (TGraph*)pFile->Get(getname);
				S->sigevents[ch][hilo]->SetMarkerColor(color[ch]);
				S->sigevents[ch][hilo]->SetMarkerStyle(style[ch]);
				S->sigevents[ch][hilo]->SetMarkerSize(size);
			}
		}
	}
	return 0;
}

//////////////////////////////////////////////////////////////////////////////
// Geninfo stuff
//____________________________________________________________________________
int SSDLDumper::muIndexToBin(int ind){
	// For the origin histograms
	// return the bin to fill for each id/type
	int id    = abs(MuGenID[ind]);
	int mid   = abs(MuGenMID[ind]);
	int mtype = abs(MuGenMType[ind]);
	if(id  != 13)                                 return 1; // mis id
	if(mid == 24)                                 return 2; // W
	if(mid == 23)                                 return 3; // Z
	if(mtype == 2)                                return 4; // tau
	if(mtype == 11 || mtype == 12 || mtype == 18) return 5; // light hadrons
	if(mtype == 13 || mtype == 19)                return 6; // strange hadrons
	if(mtype == 14 || mtype == 16 || mtype == 20) return 7; // charmed hadrons
	if(mtype == 15 || mtype == 17 || mtype == 21) return 8; // bottom hadrons
	if(mtype == 91 || mtype == 92)                return 9; // pythia strings
	return 15;                                              // uid
}
int SSDLDumper::elIndexToBin(int ind){
	// For the origin histograms
	// return the bin to fill for each id/type
	int id    = abs(ElGenID[ind]);
	int type  = abs(ElGenType[ind]);
	int mid   = abs(ElGenMID[ind]);
	int mtype = abs(ElGenMType[ind]);
	if(id  != 11){                                 // mis id
		if(type == 0 || type == 2)                 return 1;  // mis-match
		if(id == 22)                               return 2;  // gamma
		if(type == 11 || type == 12 || type == 13 ||
		   type == 18 || type == 19)               return 3;  // Hadr. fake
		return 15;                                            // uid
	}
	if(mid == 24)                                  return 4;  // W
	if(mid == 23)                                  return 5;  // Z
	if(mtype == 2)                                 return 6;  // tau

	if(mtype == 11 || mtype == 12 || mtype == 18)  return 7;  // light hadrons
	if(mtype == 13 || mtype == 19)                 return 8;  // strange hadrons
	if(mtype == 14 || mtype == 16 || mtype == 20)  return 9;  // charmed hadrons
	if(mtype == 15 || mtype == 17 || mtype == 21)  return 10; // bottom hadrons
	if(mtype == 91 || mtype == 92)                 return 11; // pythia strings
	return 15;                                                // uid
}
TString SSDLDumper::muBinToLabel(int bin){
	// For the origin histograms
	// return the bin label for each bin
	switch( bin ){
		case 1:  return "Fake";
		case 2:  return "W";
		case 3:  return "Z";
		case 4:  return "#tau";
		case 5:  return "Light had.";
		case 6:  return "Strange had.";
		case 7:  return "Charmed had.";
		case 8:  return "Bottom had.";
		case 9:  return "QCD String";
		case 10: return "";
		case 11: return "";
		case 12: return "";
		case 13: return "";
		case 14: return "";
		case 15: return "Unidentified";
		default: return "?";
	}
}
TString SSDLDumper::elBinToLabel(int bin){
	// For the origin histograms
	// return the bin label for each bin
	switch( bin ){
		case 1:  return "Mismatch (#mu, #nu, etc.)";
		case 2:  return "Gamma fake / Conversion";
		case 3:  return "Hadronic fake";
		case 4:  return "W";
		case 5:  return "Z";
		case 6:  return "#tau";
		case 7:  return "Light had.";
		case 8:  return "Strange had.";
		case 9:  return "Charmed had.";
		case 10: return "Bottom had.";
		case 11: return "QCD string";
		case 12: return "";
		case 13: return "";
		case 14: return "";
		case 15: return "Unidentified";
		default: return "?";
	}
}
void SSDLDumper::labelOriginAxis(TAxis *axis, gChannel chan){
	if(chan == Muon){
		axis->SetTitle("#mu Origin");
		for(size_t i = 1; i <= 15; ++i){
			axis->SetBinLabel(i, muBinToLabel(i));
		}		
	}
	if(chan == Elec){
		axis->SetTitle("e Origin");
		for(size_t i = 1; i <= 15; ++i){
			axis->SetBinLabel(i, elBinToLabel(i));
		}
	}
	return;
}
void SSDLDumper::label2OriginAxes(TAxis *axis1, TAxis *axis2, gChannel chan){
	if(chan == Muon || chan == Elec){
		labelOriginAxis(axis1, chan);
		labelOriginAxis(axis2, chan);		
	}
	if(chan == ElMu){
		labelOriginAxis(axis1, Muon);
		labelOriginAxis(axis2, Elec);		
	}
	return;
}

//////////////////////////////////////////////////////////////////////////////
// Trigger stuff:
//____________________________________________________________________________
bool SSDLDumper::mumuSignalTrigger(){
	return ( doubleMuTrigger() || doubleMuHTTrigger() );
}
bool SSDLDumper::elelSignalTrigger(){
	return ( doubleElTrigger() || doubleElHTTrigger() );
}
bool SSDLDumper::elmuSignalTrigger(){
	return ( eMuTrigger() || eMuHTTrigger() );
}

//____________________________________________________________________________
// The following triggers will always return true for MC samples, and they
// will always return false if they are called in the 'wrong' dataset
bool SSDLDumper::singleMuTrigger(){
	// Pretend MC samples always fire trigger
	if(fSample->datamc > 0) return true;
	return (HLT_MU8_JET40 > 0);
}
float SSDLDumper::singleMuPrescale(){
	// Pretend MC samples have prescale 1.
	if(fSample->datamc > 0) return 1.;
	// Get the prescale factor for whichever of these triggers fired
	// Only correct if they are mutually exclusive!
	if(HLT_MU8_JET40_PS > 0) return HLT_MU8_JET40_PS;
	return 1;
}
bool SSDLDumper::singleElTrigger(){
	// Pretend MC samples always fire trigger
	if(fSample->datamc > 0) return true;
	return (HLT_ELE8_JET40 > 0);
}
float SSDLDumper::singleElPrescale(){
	// Pretend MC samples have prescale 1.
	if(fSample->datamc > 0) return 1.;
	// Get the prescale factor for whichever of these triggers fired
	// Only correct if they are mutually exclusive!
	if( HLT_ELE8_JET40_PS > 0 ) return HLT_ELE8_JET40_PS;
	return 1.;
}

bool SSDLDumper::doubleMuTrigger(){
	// Pretend MC samples always fire trigger
	if(fSample->datamc > 0) return true;
	// Only look for mumu events
	if(fSample->chansel != -1 && fSample->chansel != 0) return false;
	return ( (HLT_DOUBLEMU7 > 0) || 
	         (HLT_MU13_MU8  > 0) );
}
bool SSDLDumper::doubleElTrigger(){
	// Pretend MC samples always fire trigger
	if(fSample->datamc > 0) return true;
	// Only look for elel events
	if(fSample->chansel != -1 && fSample->chansel != 1) return false;

	return ( (HLT_ELE17_ELE8       > 0) ||
	         (HLT_ELE17_ELE8_TIGHT > 0) );
}
bool SSDLDumper::doubleMuHTTrigger(){
	// Pretend MC samples always fire trigger
	if(fSample->datamc > 0) return true;
	// Only look for mumu events
	if(fSample->chansel != -1 && fSample->chansel != 0) return false;

	return (HLT_DOUBLEMU3_HT160 > 0);
}
bool SSDLDumper::doubleElHTTrigger(){
	// Pretend MC samples always fire trigger
	if(fSample->datamc > 0) return true;
	// Only look for elel events
	if(fSample->chansel != -1 && fSample->chansel != 1) return false;

	return ( (HLT_DOUBLEELE8_HT160 > 0) ||
	         (HLT_DOUBLEELE8_HT160_TIGHT > 0) );
}
bool SSDLDumper::eMuTrigger(){
	// Pretend MC samples always fire trigger
	if(fSample->datamc > 0) return true;
	// Only look for emu events
	if(fSample->chansel != -1 && fSample->chansel != 2) return false;
	return ( (HLT_MU17_ELE8       > 0) || 
	         (HLT_MU8_ELE17       > 0) ||
	         (HLT_MU17_ELE8_TIGHT > 0) ||
	         (HLT_MU8_ELE17_TIGHT > 0) );
}
bool SSDLDumper::eMuHTTrigger(){
	// Pretend MC samples always fire trigger
	if(fSample->datamc > 0) return true;
	// Only look for emu events
	if(fSample->chansel != -1 && fSample->chansel != 2) return false;

	return ( (HLT_MU3_ELE8_HT160 > 0) ||
	         (HLT_MU3_ELE8_HT160_TIGHT > 0) );
}

//////////////////////////////////////////////////////////////////////////////
// Helper functions:
//____________________________________________________________________________
int SSDLDumper::getNJets(){
	int njets(0);
	for(size_t i = 0; i < NJets; ++i) if(isGoodJet(i)) njets++;
	return njets;
}
int SSDLDumper::getNBTags(){
	int ntags(0);
	for(size_t i = 0; i < NJets; ++i) if(isGoodJet(i) && JetSSVHPBTag[i] > 2.0) ntags++;
	return ntags;
}
float SSDLDumper::getHT(){
	float ht(0.);
	for(size_t i = 0; i < NJets; ++i) if(isGoodJet(i)) ht += JetPt[i];
	return ht;
}
float SSDLDumper::getMT2(int ind1, int ind2, int toggle){
	// Calculate MT2 variable for two leptons and missing energy,
	// assuming zero testmass
	// Toggle switches between mumu (1), ee(2), emu(3)
	double pa[3];
	double pb[3];
	double pmiss[3];

	TLorentzVector pmet, pl1, pl2;

	if(toggle == 1){ // mumu
		pl1.SetPtEtaPhiM(MuPt[ind1], MuEta[ind1], MuPhi[ind1], gMMU);
		pl2.SetPtEtaPhiM(MuPt[ind2], MuEta[ind2], MuPhi[ind2], gMMU);			
	}
	if(toggle == 2){ // ee
		pl1.SetPtEtaPhiM(ElPt[ind1], ElEta[ind1], ElPhi[ind1], gMEL);
		pl2.SetPtEtaPhiM(ElPt[ind2], ElEta[ind2], ElPhi[ind2], gMEL);			
	}
	if(toggle == 3){ // emu
		pl1.SetPtEtaPhiM(MuPt[ind1], MuEta[ind1], MuPhi[ind1], gMMU);
		pl2.SetPtEtaPhiM(ElPt[ind2], ElEta[ind2], ElPhi[ind2], gMEL);			
	}
	
	pmet.SetPtEtaPhiM(pfMET, 0., pfMETPhi, 0.);
	pmiss[0] = 0.; // irrelevant
	pmiss[1] = pmet.Px();
	pmiss[2] = pmet.Py();

	pa[0] = 0.;
	pa[1] = pl1.Px();
	pa[2] = pl1.Py();

	pb[0] = 0.;
	pb[1] = pl2.Px();
	pb[2] = pl2.Py();
	
	Davismt2 *DavisMT2 = new Davismt2();
	DavisMT2->set_verbose(0);
	DavisMT2->set_momenta(pa, pb, pmiss);
	DavisMT2->set_mn(0.); // testmass
	double MT2 = DavisMT2->get_mt2();
	delete DavisMT2;
	return MT2;
}
int   SSDLDumper::getClosestJet(int ind, gChannel chan){
// Get index of the closest jet
	float lepeta = (chan==Muon)?MuEta[ind]:ElEta[ind];
	float lepphi = (chan==Muon)?MuPhi[ind]:ElPhi[ind];
	
	float mindr = 999.;
	int cljetindex = -1;
	for(size_t i = 0; i < NJets; ++i){
		if(isGoodJet(i) == false) continue;
		float dr = Util::GetDeltaR(lepeta, JetEta[i], lepphi, JetPhi[i]);
		if(dr > mindr) continue;
		mindr = dr;
		cljetindex = i;
	}
	return cljetindex;
}
float SSDLDumper::getClosestJetPt(int ind, gChannel chan){
// Get the pt of the closest jet
	int jind = getClosestJet(ind, chan);
	if(jind > -1) return JetPt[jind];
	return -1;
}
float SSDLDumper::getClosestJetDR(int ind, gChannel chan){
// Get delta R to the closest jet
	float lepeta = (chan==Muon)?MuEta[ind]:ElEta[ind];
	float lepphi = (chan==Muon)?MuPhi[ind]:ElPhi[ind];
	int jind = getClosestJet(ind, chan);
	if(jind > -1) return Util::GetDeltaR(lepeta, JetEta[jind], lepphi, JetPhi[jind]);
	return -1;
}
float SSDLDumper::getSecondClosestJetDR(int ind, gChannel chan){
// Get the pt of the closest jet
	float lepeta = (chan==Muon)?MuEta[ind]:ElEta[ind];
	float lepphi = (chan==Muon)?MuPhi[ind]:ElPhi[ind];
	
	if(NJets < 2) return 10.;
	
	float mindr  = 888.;
	float mindr2 = 999.;
	for(size_t i = 0; i < NJets; ++i){
		float dr = Util::GetDeltaR(lepeta, JetEta[i], lepphi, JetPhi[i]);
		if(dr < mindr){
			mindr2 = mindr;
			mindr = dr;			
		}
		else if(dr < mindr2){
			mindr2 = dr;
		}
	}
	return mindr2;
}
float SSDLDumper::getAwayJetPt(int ind, gChannel chan){
// Get the pt of away jet
// DR > 0.1, choose hardest
	float lepeta = (chan==Muon)?MuEta[ind]:ElEta[ind];
	float lepphi = (chan==Muon)?MuPhi[ind]:ElPhi[ind];
	
	if(NJets < 1) return 0.;
	
	float mindr  = 888.;
	for(size_t i = 0; i < NJets; ++i){
		if(Util::GetDeltaR(lepeta, JetEta[i], lepphi, JetPhi[i]) < 1.0) continue;
		if(!isGoodJet(i)) continue;
		return JetPt[i]; // assume sorted by pt, so this will return the hardest one
	}
	return 0.;
}
float SSDLDumper::getMaxJPt(){
	float maxpt(0.);
	for(size_t i = 0; i < NJets; ++i){
		if(!isGoodJet(i)) continue;
		if(JetPt[i] < maxpt) continue;
		maxpt = JetPt[i];
	}
	return maxpt;
}

void SSDLDumper::resetHypLeptons(){
	TLorentzVector vec(0., 0., 0., 0.);
	fHypLepton1 = lepton(vec, 0, -1, -1);
	fHypLepton2 = lepton(vec, 0, -1, -1);
}
void SSDLDumper::setHypLepton1(int index, gChannel chan){
	TLorentzVector vec(0., 0., 0., 0.);
	if(chan == Muon){
		vec.SetPtEtaPhiM(MuPt[index], MuEta[index], MuPhi[index], gMMU);
		fHypLepton1 = lepton(vec, MuCharge[index], 0, index);
	}
	else if(chan == Elec){
		vec.SetPtEtaPhiM(ElPt[index], ElEta[index], ElPhi[index], gMEL);
		fHypLepton1 = lepton(vec, ElCharge[index], 1, index);
	}
	else exit(-1);
}
void SSDLDumper::setHypLepton2(int index, gChannel chan){
	TLorentzVector vec(0., 0., 0., 0.);
	if(chan == Muon){
		vec.SetPtEtaPhiM(MuPt[index], MuEta[index], MuPhi[index], gMMU);
		fHypLepton2 = lepton(vec, MuCharge[index], 0, index);
	}
	else if(chan == Elec){
		vec.SetPtEtaPhiM(ElPt[index], ElEta[index], ElPhi[index], gMEL);
		fHypLepton2 = lepton(vec, ElCharge[index], 1, index);
	}
	else exit(-1);
}

//////////////////////////////////////////////////////////////////////////////
// Event Selections:
//____________________________________________________________________________
bool SSDLDumper::isGoodEvent(){
	// Global cuts?
	return true;
}
bool SSDLDumper::isGoodMuEvent(){
	// Ask for >0 loose muons, if 2 muons ask for second to be loose too
	// return number of loose muons
	if(!isGoodEvent()) return false;
	if(NMus < 1) return false;
	if(isLooseMuon(0) == false) return false;
	if(NMus > 1) if(isLooseMuon(1) == false) return false;
	return true;
}

int SSDLDumper::hasLooseMuons(int &mu1, int &mu2){
	// Returns the number of loose muons and fills their indices in mu1 and mu2
	// Assumes the muons are sorted by pt in the minitree
	vector<int> loosemus;
	mu1 = -1;
	mu2 = -1;
	for(size_t i = 0; i < NMus; ++i) if(isLooseMuon(i)) loosemus.push_back(i);
	if(loosemus.size() > 0) mu1 = loosemus[0];
	if(loosemus.size() > 1) mu2 = loosemus[1];
	return loosemus.size();
}
int SSDLDumper::hasLooseMuons(){
	int ind1(-1), ind2(-1);
	return hasLooseMuons(ind1, ind2);
}
int SSDLDumper::hasLooseElectrons(int &el1, int &el2){
	// Returns the number of loose electrons and fills their indices in el1 and el2
	// Assumes the electrons are sorted by pt in the minitree
	vector<int> looseels;
	el1 = -1;
	el2 = -1;
	for(size_t i = 0; i < NEls; ++i) if(isLooseElectron(i)) looseels.push_back(i);
	if(looseels.size() > 0) el1 = looseels[0];
	if(looseels.size() > 1) el2 = looseels[1];
	return looseels.size();
}
int SSDLDumper::hasLooseElectrons(){
	int ind1(-1), ind2(-1);
	return hasLooseElectrons(ind1, ind2);
}

//____________________________________________________________________________
int SSDLDumper::isSSLLEvent(int &ind1, int &ind2){
	// Check first if there is a tight-tight pair in the event
	int res = isSSEvent(ind1, &SSDLDumper::isTightMuon, ind2, &SSDLDumper::isTightElectron);
	if(res > 0) return res;
	return isSSEvent(ind1, &SSDLDumper::isLooseMuon, ind2, &SSDLDumper::isLooseElectron);
}
int SSDLDumper::isOSLLEvent(int &ind1, int &ind2){
	// Check first if there is a tight-tight pair in the event
	int res = isOSEvent(ind1, &SSDLDumper::isTightMuon, ind2, &SSDLDumper::isTightElectron);
	if(res > 0) return res;
	else return isOSEvent(ind1, &SSDLDumper::isLooseMuon, ind2, &SSDLDumper::isLooseElectron);
}
int SSDLDumper::isSSEvent(int &ind1, bool(SSDLDumper::*muonSelector)(int), int &ind2, bool(SSDLDumper::*eleSelector)(int)){
	// Looks for a SS pair of leptons with given object selectors
	// Return the channel: 0 = none found
	//                     1 / -1 = mu+mu+ / mu-mu- pair
	//                     2 / -2 = e+e+   / e-e-   pair
	//                     3 / -3 = mu+e+  / mu-e-  pair
	// The indices in the argument given are sorted by pt unless
	// it's a e/mu event when they are sorted such that the muon
	// is ind1

	// The pair selected is the one with hardest pt1 + pt2
	vector<lepton> tmp_Leptons_p;
	vector<lepton> tmp_Leptons_m;

	// First store all loose leptons in two vectors according to their charges
	TLorentzVector plep;
	for(size_t i = 0; i < NMus; ++i){
		if((*this.*muonSelector)(i) == false) continue;
		if(MuCharge[i] == 1 ){
			plep.SetPtEtaPhiM(MuPt[i], MuEta[i], MuPhi[i], gMMU);
			lepton tmpLepton(plep, 1, 0, i);
			tmp_Leptons_p.push_back(tmpLepton);
		}
		if(MuCharge[i] == -1){
			plep.SetPtEtaPhiM(MuPt[i], MuEta[i], MuPhi[i], gMMU);
			lepton tmpLepton(plep, -1, 0, i);
			tmp_Leptons_m.push_back(tmpLepton);
		}
	}
	for(size_t i = 0; i < NEls; ++i){
		if((*this.*eleSelector)(i) == false) continue;
		if(ElCharge[i] == 1 ){
			plep.SetPtEtaPhiM(ElPt[i], ElEta[i], ElPhi[i], gMEL);
			lepton tmpLepton(plep, 1, 1, i);
			tmp_Leptons_p.push_back(tmpLepton);
		}
		if(ElCharge[i] == -1){
			plep.SetPtEtaPhiM(ElPt[i], ElEta[i], ElPhi[i], gMEL);
			lepton tmpLepton(plep, -1, 1, i);
			tmp_Leptons_m.push_back(tmpLepton);
		}
	}

	// Check for at least one pair
	if(tmp_Leptons_m.size() < 2 && tmp_Leptons_p.size() < 2) return 0;

	/////////////////////////////////////////////////////////////////////////
	// Sort these vectors by type and pt
	vector<lepton> Leptons_p = sortLeptonsByTypeAndPt(tmp_Leptons_p);
	vector<lepton> Leptons_m = sortLeptonsByTypeAndPt(tmp_Leptons_m);

	// Proceed to select one ss pair
	double ptsum1(0.), ptsum2(0.);
	int typesum1(9), typesum2(9); // 0 for mm, 1 for em, 2 for ee
	int select(0); // switch between the two possible pairs
	if(Leptons_p.size() > 1){
		ptsum1   = Leptons_p[0].p.Pt() + Leptons_p[1].p.Pt();
		typesum1 = Leptons_p[0].type   + Leptons_p[1].type;
	}
	if(Leptons_m.size() > 1){
		ptsum2   = Leptons_m[0].p.Pt() + Leptons_m[1].p.Pt();
		typesum2 = Leptons_m[0].type   + Leptons_m[1].type;
	}
	// Selection logic:
	if(typesum1 < typesum2)                     select = 1;  // first pair had more muons
	if(typesum1 > typesum2)                     select = -1; // second pair has more muons
	if(typesum1 == typesum2 && ptsum1 > ptsum2) select = 1;  // both have same #muons, select by ptsum
	if(typesum1 == typesum2 && ptsum1 < ptsum2) select = -1;

	vector<lepton> selectedPair;
	if(select == 1){ // positive
		selectedPair.push_back(Leptons_p[0]);
		selectedPair.push_back(Leptons_p[1]);
	}
	if(select == -1){ // negative
		selectedPair.push_back(Leptons_m[0]);
		selectedPair.push_back(Leptons_m[1]);
	}
	/////////////////////////////////////////////////////////////////////////

	int result = 0;
	if(selectedPair[0].type == 0 && selectedPair[1].type == 0) result = 1; // mu/mu
	if(selectedPair[0].type == 1 && selectedPair[1].type == 1) result = 2; // el/el
	if(selectedPair[0].type == 0 && selectedPair[1].type == 1) result = 3; // mu/el
	result *= select; // Add charge to result

	// Return values
	ind1 = selectedPair[0].index;
	ind2 = selectedPair[1].index;
	return result;
}
int SSDLDumper::isOSEvent(int &ind1, bool(SSDLDumper::*muonSelector)(int), int &ind2, bool(SSDLDumper::*eleSelector)(int)){
	// Looks for a OS pair of leptons with given object selectors
	// Return the channel: 0 = none found
	//                     1 / -1 = mu+mu- / mu-mu+ pair
	//                     2 / -2 = e+e-   / e-e+   pair
	//                     3 / -3 = mu+e-  / mu-e+  pair
	// The indices in the argument given are sorted by pt unless
	// it's a e/mu event when they are sorted such that the muon
	// is ind1
	// Return value has sign of harder of the two, or the mu if
	// it's a mu/e pair
	// The event selected is the one with hardest pt1 + pt2
	vector<lepton> tmp_Leptons;

	// First store all loose leptons in a vector
	for(size_t i = 0; i < NMus; ++i){
		if((*this.*muonSelector)(i) == false) continue;
		TLorentzVector pmu;
		pmu.SetPtEtaPhiM(MuPt[i], MuEta[i], MuPhi[i], gMMU);
		lepton tmpLepton(pmu, MuCharge[i], 0, i);
		tmp_Leptons.push_back(tmpLepton);
	}
	for(size_t i = 0; i < NEls; ++i){
		if((*this.*eleSelector)(i) == false) continue;
		TLorentzVector p;
		p.SetPtEtaPhiM(ElPt[i], ElEta[i], ElPhi[i], gMEL);
		lepton tmpLepton(p, ElCharge[i], 1, i);
		tmp_Leptons.push_back(tmpLepton);
	}

	// Sort these vector by their flavor and their pt
	vector<lepton> v_Leptons;
	v_Leptons = sortLeptonsByTypeAndPt(tmp_Leptons);

	// Proceed to select one os pair
	if(v_Leptons.size() < 2) return 0;

	vector<lepton> selectedPair;
	selectedPair.push_back(v_Leptons[0]);
	for(size_t i = 1; i < v_Leptons.size(); ++i){ // look for the next lep with opp. charge
		if(selectedPair[0].charge == v_Leptons[i].charge) continue;
		selectedPair.push_back(v_Leptons[i]);
		break;
	}
	if(selectedPair.size() < 2) return 0;

	int result = 0;
	if(selectedPair[0].type == 0 && selectedPair[1].type == 0) result = 1; // mu/mu
	if(selectedPair[0].type == 1 && selectedPair[1].type == 1) result = 2; // el/el
	if(selectedPair[0].type == 0 && selectedPair[1].type == 1) result = 3; // mu/el
	result *= selectedPair[0].charge; // Add charge to result

	// Result
	ind1 = selectedPair[0].index;
	ind2 = selectedPair[1].index;
	return result;
}

//____________________________________________________________________________
bool momentumComparator(SSDLDumper::lepton i, SSDLDumper::lepton j){ return (i.p.Pt()>j.p.Pt()); }
bool momentumAndTypeComparator(SSDLDumper::lepton i, SSDLDumper::lepton j){
	// If e/mu, return the muon irrespective of type
	if(i.type != j.type) return i.type < j.type;
	return momentumComparator(i,j);
}
vector<SSDLDumper::lepton> SSDLDumper::sortLeptonsByPt(vector<lepton>& leptons){
	vector<lepton> theLep = leptons;
	sort (theLep.begin(), theLep.end(), momentumComparator);
	return theLep;
}
vector<SSDLDumper::lepton> SSDLDumper::sortLeptonsByTypeAndPt(vector<lepton>& leptons){
	vector<lepton> theLep = leptons;
	sort (theLep.begin(), theLep.end(), momentumAndTypeComparator);
	return theLep;
}

//____________________________________________________________________________
bool SSDLDumper::passesNJetCut(int cut){
	return getNJets() >= cut;
}
bool SSDLDumper::passesJet50Cut(){
	// Return true if event contains one good jet with pt > 50
	for(size_t i = 0; i < NJets; ++i) if(isGoodJet(i, 50)) return true;
	return false;
}

//____________________________________________________________________________
bool SSDLDumper::passesHTCut(float min, float max){
	return (getHT() >= min && getHT() < max);
}
bool SSDLDumper::passesMETCut(float min, float max){
	return (pfMET >= min && pfMET < max);
}
bool SSDLDumper::passesZVeto(bool(SSDLDumper::*muonSelector)(int), bool(SSDLDumper::*eleSelector)(int), float dm){
// Checks if any combination of opposite sign, same flavor leptons (e or mu)
// has invariant mass closer than dm to the Z mass, returns true if none found
// Default for dm is 15 GeV
	if(NMus > 1){
		for(size_t i = 0; i < NMus-1; ++i){
			if((*this.*muonSelector)(i)){
				TLorentzVector pmu1, pmu2;
				pmu1.SetPtEtaPhiM(MuPt[i], MuEta[i], MuPhi[i], gMMU);

				// Second muon
				for(size_t j = i+1; j < NMus; ++j){ 
					if((*this.*muonSelector)(j) && (MuCharge[i] != MuCharge[j]) ){
						pmu2.SetPtEtaPhiM(MuPt[j], MuEta[j], MuPhi[j], gMMU);
						if(fabs((pmu1+pmu2).M() - gMZ) < dm) return false;
					}
				}
			}
		}
	}

	if(NEls > 1){
		for(size_t i = 0; i < NEls-1; ++i){
			if((*this.*eleSelector)(i)){
				TLorentzVector pel1, pel2;
				pel1.SetPtEtaPhiM(ElPt[i], ElEta[i], ElPhi[i], gMEL);

				// Second electron
				for(size_t j = i+1; j < NEls; ++j){
					if((*this.*eleSelector)(j) && (ElCharge[i] != ElCharge[j]) ){
						pel2.SetPtEtaPhiM(ElPt[j], ElEta[j], ElPhi[j], gMEL);
						if(fabs((pel1+pel2).M() - gMZ) < dm) return false;
					}
				}
			}
		}		
	}
	return true;
}
bool SSDLDumper::passesZVeto(float dm){
	return passesZVeto(&SSDLDumper::isTightMuon, &SSDLDumper::isTightElectron, dm);
	// return passesZVeto(&SSDLDumper::isLooseMuon, &SSDLDumper::isLooseElectron, dm);
}
bool SSDLDumper::passesMllEventVeto(float cut){
// Checks if any combination of opposite sign, same flavor tight leptons (e or mu)
// has invariant mass smaller than cut, returns true if none found
	if(NMus > 1){
		// First muon
		for(size_t i = 0; i < NMus-1; ++i){
			if(isTightMuon(i)){
				TLorentzVector pmu1, pmu2;
				pmu1.SetPtEtaPhiM(MuPt[i], MuEta[i], MuPhi[i], gMMU);

				// Second muon
				for(size_t j = i+1; j < NMus; ++j){ 
					if(isTightMuon(j)){
						pmu2.SetPtEtaPhiM(MuPt[j], MuEta[j], MuPhi[j], gMMU);
						if((pmu1+pmu2).M() < cut) return false;
					}
				}
			}
		}
	}

	if(NEls > 1){ 
	// First electron
		for(size_t i = 0; i < NEls-1; ++i){
			if(isTightElectron(i)){
				TLorentzVector pel1, pel2;
				pel1.SetPtEtaPhiM(ElPt[i], ElEta[i], ElPhi[i], gMEL);

				// Second electron
				for(size_t j = i+1; j < NEls; ++j){
					if(isTightElectron(j)){
						pel2.SetPtEtaPhiM(ElPt[j], ElEta[j], ElPhi[j], gMEL);
						if((pel1+pel2).M() < cut) return false;
					}
				}
			}
		}
	}
	return true;
}
bool SSDLDumper::passesMllEventVeto(int ind1, int ind2, int toggle, float cut){
	// Calculate invariant mass of pair, return false if it's smaller than cut
	TLorentzVector pmet, pl1, pl2;

	if(toggle == 1){ // mumu
		pl1.SetPtEtaPhiM(MuPt[ind1], MuEta[ind1], MuPhi[ind1], gMMU);
		pl2.SetPtEtaPhiM(MuPt[ind2], MuEta[ind2], MuPhi[ind2], gMMU);			
	}
	if(toggle == 2){ // ee
		pl1.SetPtEtaPhiM(ElPt[ind1], ElEta[ind1], ElPhi[ind1], gMEL);
		pl2.SetPtEtaPhiM(ElPt[ind2], ElEta[ind2], ElPhi[ind2], gMEL);			
	}
	if(toggle == 3){ // emu
		pl1.SetPtEtaPhiM(MuPt[ind1], MuEta[ind1], MuPhi[ind1], gMMU);
		pl2.SetPtEtaPhiM(ElPt[ind2], ElEta[ind2], ElPhi[ind2], gMEL);			
	}
	if( (pl1+pl2).M() < cut ) return false;
	return true;
}

//____________________________________________________________________________
bool SSDLDumper::isGoodRun(Sample *S){
	// Select runs such that JetB and MultiJet datasets are mutually exclusive
	// if(gSample(sample) == JMB)      if(Run > 147195) return false;
	// if(gSample(sample) == MultiJet) if(Run < 147196) return false;
	// if(sample == JMB)      if(Run > 148058) return false;
	// if(sample == MultiJet) if(Run < 148822) return false;
	return true;
}

//____________________________________________________________________________
bool SSDLDumper::isSigSupMuEvent(){
	int mu1(-1), mu2(-1);
	if(hasLooseMuons(mu1, mu2) < 1) return false;
	setHypLepton1(mu1, Muon);
	if(!passesJet50Cut())           return false;
	if(!passesNJetCut(1))           return false;
	// if(!passesNJetCut(fC_minNjets)) return false;
	if(MuMT[0] > fC_maxMt_Control)  return false;
	if(pfMET > fC_maxMet_Control)   return false;
	if(NMus > 1)                    return false;
	return true;
}
bool SSDLDumper::isZMuMuEvent(){
	int mu1(-1), mu2(-1);
	if(hasLooseMuons(mu1, mu2) < 2)  return false;
	// if(!isTightMuon(0) && !isTightMuon(1)) return false; // at least one tight

	if(MuCharge[mu1] == MuCharge[mu2]) return false; // os

	// Z mass window cut
	TLorentzVector p1, p2;
	p1.SetPtEtaPhiM(MuPt[mu1], MuEta[mu1], MuPhi[mu1], 0.1057);
	p2.SetPtEtaPhiM(MuPt[mu2], MuEta[mu2], MuPhi[mu2], 0.1057);
	double m = (p1+p2).M();
	if(fabs(gMZ - m) > 15.) return false;

	setHypLepton1(mu1, Muon);
	setHypLepton2(mu2, Muon);

	if(pfMET > 20.) return false;

	if(!passesNJetCut(2)) return false;
	return true;
}
bool SSDLDumper::isSigSupElEvent(){
	int el1(-1), el2(-1);
	if(hasLooseElectrons(el1, el2) < 1) return false;
	setHypLepton1(el1, Elec);
	if(!passesJet50Cut())               return false;
	if(!passesNJetCut(1))               return false;
	// if(!passesNJetCut(fC_minNjets))     return false;
	if(ElMT[0] > fC_maxMt_Control)      return false;
	if(pfMET > fC_maxMet_Control)       return false;
	if(NEls > 1)                        return false;
	return true;
}
bool SSDLDumper::isZElElEvent(int &elind){
	int el1(-1), el2(-1);
	if(hasLooseElectrons(el1, el2) < 2)  return false;
	if(!isTightElectron(el1) && !isTightElectron(el2)) return false; // at least one tight

	if(ElCharge[el1] == ElCharge[el2]) return false; // os

	// Z mass window cut
	TLorentzVector p1, p2;
	p1.SetPtEtaPhiM(ElPt[el1], ElEta[el1], ElPhi[el1], gMEL);
	p2.SetPtEtaPhiM(ElPt[el2], ElEta[el2], ElPhi[el2], gMEL);
	double m = (p1+p2).M();
	if(fabs(gMZ - m) > 15.) return false;

	setHypLepton1(el1, Elec);
	setHypLepton2(el2, Elec);


	// If only the harder one tight or both tight, return the softer one
	// If only the softer one tight, return the harder one
	elind = el2;
	if(isTightElectron(el2) && !isTightElectron(el1)) elind = el1;

	if(pfMET > 20.) return false;
	if(!passesNJetCut(2)) return false;
	return true;
}

//____________________________________________________________________________
bool SSDLDumper::isGenMatchedSUSYDiLepEvent(){
	int ind1(-1), ind2(-1);
	return isGenMatchedSUSYDiLepEvent(ind1, ind2);
}
bool SSDLDumper::isGenMatchedSUSYDiLepEvent(int &mu1, int &mu2){
	if(!isSSLLMuEvent(mu1, mu2)) return false;
	if(!isTightMuon(mu1) || !isTightMuon(mu2)) return false;
	if(isPromptMuon(mu1) && isPromptMuon(mu2)) return true;
	return false;
}
bool SSDLDumper::isGenMatchedSUSYEEEvent(){
	int ind1(-1), ind2(-1);
	if(!isSSLLElEvent(ind1, ind2)) return false;
	if(!isTightElectron(ind1) || !isTightElectron(ind2)) return false;
	if(isPromptElectron(ind1) && isPromptElectron(ind2)) return true;
	return false;
}
bool SSDLDumper::isGenMatchedSUSYEMuEvent(){
	int muind(-1), elind(-1);
	if(!isSSLLElMuEvent(muind, elind)) return false;
	if(!isTightElectron(elind) || !isTightMuon(muind)) return false;
	if(isPromptMuon(muind) && isPromptElectron(elind)) return true;
	return false;
}

//____________________________________________________________________________
bool SSDLDumper::isSSLLMuEvent(int& mu1, int& mu2){
	// This should include all the cuts for the final selection
	int nmus = hasLooseMuons(mu1, mu2);
	if(nmus < 1) return false; // >0 loose muons
	if(fDoCounting) fCounter[Muon].fill(fMMCutNames[3]);
	if(nmus < 2) return false; // >1 loose muons
	if(fDoCounting) fCounter[Muon].fill(fMMCutNames[4]);

	if(fChargeSwitch == 0){
		if(abs(isSSLLEvent(mu1, mu2)) != 1) return false;
		if(fDoCounting) fCounter[Muon].fill(fMMCutNames[5]);
	}
	if(fChargeSwitch == 1){
		if(abs(isOSLLEvent(mu1, mu2)) != 1) return false;
	}

	if(fChargeSwitch == 0){
		// Need to remove the Z veto in OS case, since these OS events are exactly what 
		// contributes to the SS yield. They would NOT fire the Z veto, since they are
		// misidd as SS events
		if(!passesZVeto()) return false; // no Zs in event
		if(fDoCounting) fCounter[Muon].fill(fMMCutNames[6]);
	}

	if(!passesMllEventVeto(mu1, mu2, 1, 5.)) return false; // no low mass OSSF pairs
	if(fDoCounting) fCounter[Muon].fill(fMMCutNames[7]);

	// Define hypothesis leptons
	setHypLepton1(mu1, Muon);
	setHypLepton1(mu2, Muon);

	// if(!passesJet50Cut()) return false; // one jet with pt > 50
	// if(fDoCounting) fCounter[Muon].fill(fMMCutNames[8]);

	if(!passesNJetCut(fC_minNjets) ) return false;    // njets cut
	if(fDoCounting) fCounter[Muon].fill(fMMCutNames[9]);

	if(!passesHTCut(fC_minHT, fC_maxHT) )  return false;    // ht cut
	if(fDoCounting) fCounter[Muon].fill(fMMCutNames[10]);

	if(!passesMETCut(fC_minMet, fC_maxMet) ) return false;    // met cut
	if(fDoCounting) fCounter[Muon].fill(fMMCutNames[11]);

	if(!isGoodSecMuon(mu2)) return false; // pt cuts
	if(fDoCounting) fCounter[Muon].fill(fMMCutNames[12]);

	if(!isGoodPrimMuon(mu1)) return false;
	if(fDoCounting) fCounter[Muon].fill(fMMCutNames[13]);

	return true;
}
bool SSDLDumper::isSSLLElEvent(int& el1, int& el2){
	// This should include all the cuts for the final selection
	int nels = hasLooseElectrons(el1, el2);
	if(nels < 1) return false; // >0 eles;
	if(fDoCounting) fCounter[Elec].fill(fEECutNames[3]);
	if(nels < 2) return false; // >1 eles
	if(fDoCounting) fCounter[Elec].fill(fEECutNames[4]);

	if(fChargeSwitch == 0){
		if(abs(isSSLLEvent(el1, el2)) != 2) return false;
		if(fDoCounting) fCounter[Elec].fill(fEECutNames[5]);
	}
	if(fChargeSwitch == 1){
		if(abs(isOSLLEvent(el1, el2)) != 2) return false;
	}

	if(fChargeSwitch == 0){
		// Need to remove the Z veto in OS case, since these OS events are exactly what 
		// contributes to the SS yield. They would NOT fire the Z veto, since they are
		// misidd as SS events
		if(!passesZVeto()) return false; // no Zs in event
		if(fDoCounting) fCounter[Elec].fill(fEECutNames[6]);
	}

	if(!passesMllEventVeto(el1, el2, 2, 5.)) return false; // no low mass OSSF pairs
	if(fDoCounting) fCounter[Elec].fill(fEECutNames[7]);

	// Define hypothesis leptons
	setHypLepton1(el1, Elec);
	setHypLepton1(el2, Elec);

	// if(!passesJet50Cut()) return false;
	// if(fDoCounting) fCounter[Elec].fill(fEECutNames[8]);

	if(!passesNJetCut(fC_minNjets) ) return false;    // njets cut
	if(fDoCounting) fCounter[Elec].fill(fEECutNames[9]);

	if(!passesHTCut(fC_minHT, fC_maxHT) )  return false;    // ht cut
	if(fDoCounting) fCounter[Elec].fill(fEECutNames[10]);

	if(!passesMETCut(fC_minMet, fC_maxMet) ) return false;    // met cut
	if(fDoCounting) fCounter[Elec].fill(fEECutNames[11]);

	if(!isGoodSecElectron(el2)) return false; // pt cuts
	if(fDoCounting) fCounter[Elec].fill(fEECutNames[12]);

	if(!isGoodPrimElectron(el1)) return false;
	if(fDoCounting) fCounter[Elec].fill(fEECutNames[13]);

	return true;
}
bool SSDLDumper::isSSLLElMuEvent(int& mu, int& el){
	// This should include all the cuts for the final selection
	int nmus = hasLooseMuons(mu, el);
	if(nmus > 0 && fDoCounting) fCounter[ElMu].fill(fEMCutNames[3]);
	int nels = hasLooseElectrons(el, mu);
	if(nels > 0 && fDoCounting) fCounter[ElMu].fill(fEMCutNames[4]);
	if(nels < 1 || nmus < 1) return false;
	if(nels > 0 && fDoCounting) fCounter[ElMu].fill(fEMCutNames[5]);

	if(fChargeSwitch == 0){
		if(abs(isSSLLEvent(mu, el)) != 3) return false;
		if(fDoCounting) fCounter[ElMu].fill(fEMCutNames[6]);
	}
	if(fChargeSwitch == 1){
		if(abs(isOSLLEvent(mu, el)) != 3) return false;
		if(fDoCounting) fCounter[ElMu].fill(fEMCutNames[6]);
	}

	if(fChargeSwitch == 0){
		// Need to remove the Z veto in OS case, since these OS events are exactly what 
		// contributes to the SS yield. They would NOT fire the Z veto, since they are
		// misidd as SS events
		if(!passesZVeto()) return false;
		if(fDoCounting) fCounter[ElMu].fill(fEMCutNames[7]);
	}

	// Define hypothesis leptons
	setHypLepton1(mu, Muon);
	setHypLepton1(el, Elec);

	// if(!passesJet50Cut()) return false;
	// if(fDoCounting) fCounter[ElMu].fill(fEMCutNames[8]);

	if(!passesNJetCut(fC_minNjets) ) return false;    // njets cut
	if(fDoCounting) fCounter[ElMu].fill(fEMCutNames[9]);

	if(!passesHTCut(fC_minHT, fC_maxHT) )  return false;    // ht cut
	if(fDoCounting) fCounter[ElMu].fill(fEMCutNames[10]);

	if(!passesMETCut(fC_minMet, fC_maxMet) ) return false;    // met cut
	if(fDoCounting) fCounter[ElMu].fill(fEMCutNames[11]);

	if(MuPt[mu] > ElPt[el]){
		if(!isGoodPrimMuon(mu))    return false;
		if(fDoCounting) fCounter[ElMu].fill(fEMCutNames[12]);
		if(!isGoodSecElectron(el)) return false;
		if(fDoCounting) fCounter[ElMu].fill(fEMCutNames[13]);
	}
	else if(MuPt[mu] < ElPt[el]){
		if(!isGoodPrimElectron(el)) return false;
		if(fDoCounting) fCounter[ElMu].fill(fEMCutNames[13]);
		if(!isGoodSecMuon(mu))      return false;
		if(fDoCounting) fCounter[ElMu].fill(fEMCutNames[12]);
	}
	return true;
}

//////////////////////////////////////////////////////////////////////////////
// Object selections:
//////////////////////////////////////////////////////////////////////////////
// Muons
//____________________________________________________________________________
bool SSDLDumper::isGoodMuon(int muon, float ptcut){
	if(muon >= NMus) return false; // Sanity check
	if(ptcut < 0.) ptcut = fC_minMu2pt;
	if(MuPt[muon] < ptcut) return false;
	if(MuPtE[muon]/MuPt[muon] > 0.1) return false;
	return true;
}
bool SSDLDumper::isLooseMuon(int muon){
	if(isGoodMuon(muon) == false)  return false;
	if(MuIso[muon] > 1.00) return false;
	return true;
}
bool SSDLDumper::isTightMuon(int muon){
	if(isGoodMuon(muon) == false)  return false;
	if(isLooseMuon(muon) == false) return false;
	if(MuIso[muon] > 0.15) return false;
	// if(MuIso[muon] > 0.1) return false;
	return true;
}
bool SSDLDumper::isGoodPrimMuon(int muon, float ptcut){
	if(ptcut < 0.) ptcut = fC_minMu1pt;
	if(isLooseMuon(muon) == false) return false;
	if(MuPt[muon] < ptcut) return false;
	return true;
}
bool SSDLDumper::isGoodSecMuon(int muon, float ptcut){
	if(ptcut < 0.) ptcut = fC_minMu2pt;
	if(isLooseMuon(muon) == false) return false;
	if(MuPt[muon] < ptcut) return false;
	return true;
}

//____________________________________________________________________________
bool SSDLDumper::isFakeMuon(int muon){
	if(isLooseMuon(muon) == false) return false;
	if(MuGenMType[muon] == 2)  return false;
	if(MuGenMType[muon] == 4)  return false;
	return true;
}
bool SSDLDumper::isPromptMuon(int muon){
	if(isLooseMuon(muon) == false) return false;
	if(abs(MuGenID[muon]) != 13) return false; // not matched to mu
	if(MuGenMType[muon] == 4) return true; // W/Z -> mu
	if(abs(MuGenMID[muon]) == 15 && MuGenGMType[muon] == 4) return true; // W/Z -> tau -> mu
	if(MuGenMType[muon] == 9) return true; // susy particle
	return false;
}
bool SSDLDumper::isChargeMatchedMuon(int mu){
	if(MuGenID[mu] ==  13) return MuCharge[mu] < 0; // muon (-)
	if(MuGenID[mu] == -13) return MuCharge[mu] > 0; // anti-muon (+)
	return true;
}

//////////////////////////////////////////////////////////////////////////////
// Electrons
//____________________________________________________________________________
bool SSDLDumper::isGoodElectron(int ele, float ptcut){
// Some common electron ID cuts
	if(ele >= NEls) return false; // Sanity check
	if(ptcut < 0.) ptcut = fC_minEl2pt;
	if(ElPt[ele] < ptcut) return false;

	if(ElEcalRecHitSumEt[ele]/ElPt[ele] > 0.2) return false; // when using "CaloIsoVL" triggers

	// Reject electrons closer than 0.1 in DR to tight muons
	for(size_t i = 0; i < NMus; ++i){
		if(!isTightMuon(i)) continue;
		if(Util::GetDeltaR(MuEta[i], ElEta[ele], MuPhi[i], ElPhi[ele]) > 0.1 ) continue;
		return false;
	}
	
	
	return true;
}
bool SSDLDumper::isLooseElectron(int ele){
	// All electrons are already loose in the high-pt selection (hybiso)
	if(isGoodElectron(ele) == false) return false;
	if( fabs(ElEta[ele]) < 1.479 ) if(ElRelIso[ele] > 1.00) return false;
	else                           if(ElRelIso[ele] > 0.60) return false;		
	if(ElChIsCons[ele] != 1) return false;
	return true;
}
bool SSDLDumper::isTightElectron(int ele){
	if(!isLooseElectron(ele))       return false;
	if(ElIsGoodElId_WP80[ele] != 1) return false;

	if(ElRelIso[ele]    > 0.15) return false;
	// if(ElRelIso[ele]    > 0.1) return false;

	return true;
}
bool SSDLDumper::isGoodPrimElectron(int ele, float ptcut){
	if(ptcut < 0.) ptcut = fC_minEl1pt;
	if(isLooseElectron(ele) == false) return false;
	if(ElPt[ele] < ptcut) return false;
	return true;
}
bool SSDLDumper::isGoodSecElectron(int ele, float ptcut){
	if(ptcut < 0.) ptcut = fC_minEl2pt;
	if(isLooseElectron(ele) == false) return false;
	if(ElPt[ele] < ptcut) return false;
	return true;
}

//____________________________________________________________________________
bool SSDLDumper::isFakeElectron(int ele){
	if(isLooseElectron(ele) == false) return false;
	if(ElGenMType[ele] == 2) return false;
	if(ElGenMType[ele] == 4) return false;
	return true;
}
bool SSDLDumper::isPromptElectron(int ele){
	if(isLooseElectron(ele) == false) return false;
	if(abs(ElGenID[ele]) != 11) return false; // not matched to e
	if(ElGenMType[ele] == 4) return true; // W/Z -> el
	if(abs(ElGenMID[ele]) == 15 && ElGenGMType[ele] == 4) return true; // W/Z -> tau -> el
	if(ElGenMType[ele] == 9) return true; // susy particle
	return false;
}
bool SSDLDumper::isChargeMatchedElectron(int ele){
	if(ElGenID[ele] ==  11) return ElCharge[ele] < 0; // electron
	if(ElGenID[ele] == -11) return ElCharge[ele] > 0; // positron
	return true;
}

bool SSDLDumper::isBarrelElectron(int ele){
	// true if in barrel, false if in endcap
	if(fabs(ElEta[ele]) < 1.479 ) return true;
	return false;
}

//////////////////////////////////////////////////////////////////////////////
// Jets
//____________________________________________________________________________
bool SSDLDumper::isGoodJet(int jet, float pt){
	if(jet >= NJets) return false; // Sanity check
	float minDR = 0.4;

	// Remove jets close to hypothesis leptons
	if(fHypLepton1.index > -1) if(Util::GetDeltaR(fHypLepton1.p.Eta(), JetEta[jet], fHypLepton1.p.Phi(), JetPhi[jet]) < minDR) return false;
	if(fHypLepton2.index > -1) if(Util::GetDeltaR(fHypLepton2.p.Eta(), JetEta[jet], fHypLepton2.p.Phi(), JetPhi[jet]) < minDR) return false;

	// Remove jets close to all tight leptons
	for(size_t imu = 0; imu < NMus; ++imu){
		if(!isTightMuon(imu)) continue;
		if(!isGoodSecMuon(imu)) continue; // pt > 10
		if(Util::GetDeltaR(MuEta[imu], JetEta[jet], MuPhi[imu], JetPhi[jet]) > minDR ) continue;
		return false;
	}
	for(size_t iel = 0; iel < NEls; ++iel){
		if(!isTightElectron(iel)) continue;
		if(!isGoodSecElectron(iel)) continue;
		if(Util::GetDeltaR(ElEta[iel], JetEta[jet], ElPhi[iel], JetPhi[jet]) > minDR ) continue;
		return false;
	}
	if(JetPt[jet] < pt) return false;
	return true;
}

