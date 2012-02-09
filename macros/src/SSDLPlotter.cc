/*****************************************************************************
*   Collection of tools for producing plots for same-sign dilepton analysis  *
*                                                                            *
*                                                  (c) 2010 Benjamin Stieger *
******************************************************************************
* This class runs on ROOT files from the SSDLDumper class, reads in all the
* histograms and produces the final plots.
*****************************************************************************/
#include "SSDLPlotter.hh"

#include "helper/AnaClass.hh"
#include "helper/Utilities.hh"
#include "helper/FPRatios.hh"
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

TString SSDLPlotter::gHiLoLabel[3] = {"HighPt", "LowPt", "TauChan"};

// Charge misid probability (from Hamed)
double SSDLPlotter::gEChMisIDB   = 0.0002;
double SSDLPlotter::gEChMisIDB_E = 0.0001;
double SSDLPlotter::gEChMisIDE   = 0.0028;
double SSDLPlotter::gEChMisIDE_E = 0.0004;

//____________________________________________________________________________
SSDLPlotter::SSDLPlotter(){
// Default constructor, no samples are set
}

//____________________________________________________________________________
SSDLPlotter::SSDLPlotter(TString outputdir){
// Explicit constructor with output directory
	setOutputDir(outputdir);
}

//____________________________________________________________________________
SSDLPlotter::SSDLPlotter(TString outputdir, TString outputfile){
// Explicit constructor with output directory and output file
	setOutputDir(outputdir);
	setOutputFile(outputfile);
}

//____________________________________________________________________________
SSDLPlotter::~SSDLPlotter(){
	if(fOutputFile != NULL && fOutputFile->IsOpen()) fOutputFile->Close();
	fChain = 0;
}

//____________________________________________________________________________
void SSDLPlotter::init(TString filename){
	if(fVerbose > 0) cout << "------------------------------------" << endl;
	if(fVerbose > 0) cout << "Initializing SSDLPlotter ... " << endl;
	Util::SetStyle();
	gStyle->SetOptStat(0);

	readDatacard(filename);

	readVarNames("anavarnames.dat");
	fOutputFileName = fOutputDir + "SSDLYields.root";
	fLatex = new TLatex();
	fLatex->SetNDC(kTRUE);
	fLatex->SetTextColor(kBlack);
	fLatex->SetTextSize(0.04);

	fBinWidthScale = 10.; // Normalize Y axis to this binwidth
	fDoCounting = false; // Disable counters by default

	resetHypLeptons();
	initCutNames();
	
	// Luminosity
	// fLumiNorm = 2096.; // Pre 2011B
	// fLumiNorm = 3200.; // Including 2011B (1.014 /fb)
	fLumiNorm = 4680.; // Full 2011B
	// fLumiNorm = 1000.; // Rare Yields Sync
	// fLumiNorm = 1014.; // Only 2011B

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

	fMCBG.push_back(TTJets);
	fMCBG.push_back(TJets_t);
	fMCBG.push_back(TJets_tW);
	fMCBG.push_back(TJets_s);
	fMCBG.push_back(WJets);
	fMCBG.push_back(DYJets);
	fMCBG.push_back(GJets40);
	fMCBG.push_back(GJets100);
	fMCBG.push_back(GJets200);
	fMCBG.push_back(WW);
	fMCBG.push_back(WZ);
	fMCBG.push_back(ZZ);
	fMCBG.push_back(GVJets);

	fMCBG.push_back(DPSWW);
	fMCBG.push_back(TTbarW);
	fMCBG.push_back(TTbarZ);
	fMCBG.push_back(TTbarG);
	fMCBG.push_back(WpWp);
	fMCBG.push_back(WmWm);
	fMCBG.push_back(WWZ);
	fMCBG.push_back(WZZ);
	fMCBG.push_back(WWG);
	fMCBG.push_back(WWW);
	fMCBG.push_back(ZZZ);

	fMCBG.push_back(QCD15);
	fMCBG.push_back(QCD30);
	fMCBG.push_back(QCD50);
	fMCBG.push_back(QCD80);
	fMCBG.push_back(QCD120);
	fMCBG.push_back(QCD170);
	fMCBG.push_back(QCD300);
	fMCBG.push_back(QCD470);
	fMCBG.push_back(QCD600);
	fMCBG.push_back(QCD800);
	fMCBG.push_back(QCD1000);
	fMCBG.push_back(QCD1400);
	fMCBG.push_back(QCD1800);

	fMCBGSig = fMCBG;
	fMCBGSig.push_back(LM6);

	fMCBGMuEnr.push_back(TTJets);
	fMCBGMuEnr.push_back(TJets_t);
	fMCBGMuEnr.push_back(TJets_tW);
	fMCBGMuEnr.push_back(TJets_s);
	fMCBGMuEnr.push_back(WJets);
	fMCBGMuEnr.push_back(DYJets);
	fMCBGMuEnr.push_back(GJets40);
	fMCBGMuEnr.push_back(GJets100);
	fMCBGMuEnr.push_back(GJets200);
	fMCBGMuEnr.push_back(WW);
	fMCBGMuEnr.push_back(WZ);
	fMCBGMuEnr.push_back(ZZ);
	fMCBGMuEnr.push_back(GVJets);
	fMCBGMuEnr.push_back(DPSWW);
	fMCBGMuEnr.push_back(TTbarW);
	fMCBGMuEnr.push_back(TTbarZ);
	fMCBGMuEnr.push_back(TTbarG);
	fMCBGMuEnr.push_back(WpWp);
	fMCBGMuEnr.push_back(WmWm);
	fMCBGMuEnr.push_back(WWZ);
	fMCBGMuEnr.push_back(WZZ);
	fMCBGMuEnr.push_back(WWG);
	fMCBGMuEnr.push_back(WWW);
	fMCBGMuEnr.push_back(ZZZ);
	fMCBGMuEnr.push_back(QCDMuEnr10);

	fMCBGMuEnrSig = fMCBGMuEnr;
	fMCBGMuEnrSig.push_back(LM6);

	fMCRareSM.push_back(WZ);
	fMCRareSM.push_back(ZZ);
	fMCRareSM.push_back(GVJets);
	fMCRareSM.push_back(DPSWW);
	fMCRareSM.push_back(TTbarW);
	fMCRareSM.push_back(TTbarZ);
	fMCRareSM.push_back(TTbarG);
	fMCRareSM.push_back(WpWp);
	fMCRareSM.push_back(WmWm);
	fMCRareSM.push_back(WWZ);
	fMCRareSM.push_back(WZZ);
	fMCRareSM.push_back(WWG);
	fMCRareSM.push_back(WWW);
	fMCRareSM.push_back(ZZZ);

	fMuData    .push_back(DoubleMu1);
	fMuData    .push_back(DoubleMu2);
	fMuData    .push_back(DoubleMu3);
	fMuData    .push_back(DoubleMu4);
	fMuData    .push_back(DoubleMu5);
	fMuHadData .push_back(MuHad1);
	fMuHadData .push_back(MuHad2);
	fEGData    .push_back(DoubleEle1);
	fEGData    .push_back(DoubleEle2);
	fEGData    .push_back(DoubleEle3);
	fEGData    .push_back(DoubleEle4);
	fEGData    .push_back(DoubleEle5);
	fEleHadData.push_back(EleHad1);
	fEleHadData.push_back(EleHad2);
	fMuEGData  .push_back(MuEG1);
	fMuEGData  .push_back(MuEG2);
	fMuEGData  .push_back(MuEG3);
	fMuEGData  .push_back(MuEG4);
	fMuEGData  .push_back(MuEG5);

	fHighPtData.push_back(DoubleMu1);
	fHighPtData.push_back(DoubleMu2);
	fHighPtData.push_back(DoubleMu3);
	fHighPtData.push_back(DoubleMu4);
	fHighPtData.push_back(DoubleMu5);
	fHighPtData.push_back(DoubleEle1);
	fHighPtData.push_back(DoubleEle2);
	fHighPtData.push_back(DoubleEle3);
	fHighPtData.push_back(DoubleEle4);
	fHighPtData.push_back(DoubleEle5);
	fHighPtData.push_back(MuEG1);
	fHighPtData.push_back(MuEG2);
	fHighPtData.push_back(MuEG3);
	fHighPtData.push_back(MuEG4);
	fHighPtData.push_back(MuEG5);

	fLowPtData.push_back(MuHad1);
	fLowPtData.push_back(MuHad2);
	fLowPtData.push_back(EleHad1);
	fLowPtData.push_back(EleHad2);
	fLowPtData.push_back(MuEG1);
	fLowPtData.push_back(MuEG2);
	fLowPtData.push_back(MuEG3);
	fLowPtData.push_back(MuEG4);
}

//____________________________________________________________________________
const int     SSDLPlotter::getNPtBins (gChannel chan){
	if(chan == Muon || chan == ElMu) return gNMuPtbins;
	if(chan == Elec) return gNElPtbins;
}
const double *SSDLPlotter::getPtBins  (gChannel chan){
	if(chan == Muon || chan == ElMu) return gMuPtbins;
	if(chan == Elec) return gElPtbins;
}
const int     SSDLPlotter::getNPt2Bins(gChannel chan){
	if(chan == Muon || chan == ElMu) return gNMuPt2bins;
	if(chan == Elec) return gNElPt2bins;
}
const double *SSDLPlotter::getPt2Bins (gChannel chan){
	if(chan == Muon || chan == ElMu) return gMuPt2bins;
	if(chan == Elec) return gElPt2bins;
}
const int     SSDLPlotter::getNEtaBins(gChannel chan){
	if(chan == Muon || chan == ElMu) return gNMuEtabins;
	if(chan == Elec)            return gNElEtabins;
}
const double *SSDLPlotter::getEtaBins (gChannel chan){
	if(chan == Muon || chan == ElMu) return gMuEtabins;
	if(chan == Elec)            return gElEtabins;
}

//____________________________________________________________________________
void SSDLPlotter::doAnalysis(){	
	// sandBox();
	// return;
	
	if(readHistos(fOutputFileName) != 0) return;
	storeWeightedPred();

	// makePileUpPlots(true); // loops on all data!
	
	// printCutFlows(fOutputDir + "CutFlow.txt");
	// makeOriginPlots(Baseline);
	// makeOriginPlots(HT80MET302b);
	// makeOriginPlots(HT200MET30);
	// makeOriginPlots(HT200MET302b);
	// printOrigins(Baseline);
	// printOrigins(HT200MET30);
	// printOrigins(HT200MET302b);

	// makeMuIsolationPlots(false); // if true, loops on TTbar sample
	// makeElIsolationPlots(false); // if true, loops on TTbar sample
	// makeElIdPlots();
	// makeNT2KinPlots();
	// makeMETvsHTPlot(fMuData, fEGData, fMuEGData, HighPt);
	// makeMETvsHTPlot(fMuHadData, fEleHadData, fMuEGData, LowPt);
	// makeMETvsHTPlotPRL();
	// makeMETvsHTPlotTau();
	// makePRLPlot1();

	// makeRatioPlots(Muon);
	// makeRatioPlots(Elec);
	// makeNTightLoosePlots(Muon);
	// makeNTightLoosePlots(Elec);

	// makeFRvsPtPlots(Muon, SigSup);
	// makeFRvsPtPlots(Elec, SigSup);
	// makeFRvsPtPlots(Muon, ZDecay);
	// makeFRvsPtPlots(Elec, ZDecay);
	// makeFRvsEtaPlots(Muon);
	// makeFRvsEtaPlots(Elec);
	
	makeAllClosureTests();
	makeAllIntPredictions();

	// makeDiffPrediction();
	// makeRelIsoTTSigPlots();
	// load_msugraInfo("/scratch/mdunser/120131/msugraScan_2.root");
}

//____________________________________________________________________________
void SSDLPlotter::sandBox(){
	fOutputSubDir = "sandbox/";
	vector<int> samples;
	// samples.push_back(DoubleMu1);
	samples.push_back(DoubleMu2);
	// samples.push_back(QCDMuEnr10);

	TH1D *hdphi1_da = new TH1D("hdphi1_da", "hdphi1_da", 20, 0., 3.1416);
	TH1D *hdphi2_da = new TH1D("hdphi2_da", "hdphi2_da", 20, 0., 3.1416);
	TH1D *hdphi1_mc = new TH1D("hdphi1_mc", "hdphi1_mc", 20, 0., 3.1416);
	TH1D *hdphi2_mc = new TH1D("hdphi2_mc", "hdphi2_mc", 20, 0., 3.1416);
	
	hdphi1_da->Sumw2();
	hdphi2_da->Sumw2();
	hdphi1_mc->Sumw2();
	hdphi2_mc->Sumw2();
	
	for(size_t i = 0; i < samples.size(); ++i){
		Sample *S = fSamples[samples[i]];
		
		TTree *tree = S->getTree();
		tree->ResetBranchAddresses();
		if(S->datamc < 1) Init(tree);
		if(S->datamc > 0) InitMC(tree);
		for (Long64_t jentry=0; jentry<tree->GetEntriesFast();jentry++) {
			tree->GetEntry(jentry);
			printProgress(jentry, tree->GetEntriesFast(), S->name);

			if(singleMuTrigger() && isSigSupMuEvent()){
				setHypLepton1(0, Muon);
				int jind = getClosestJet(0, Muon);
				float dphi = Util::GetDeltaR(MuEta[0], JetEta[jind], MuPhi[0], JetPhi[jind]);
				if(getNJets() == 1){
					if(S->datamc == 0) hdphi1_da->Fill(dphi, singleMuPrescale());
					if(S->datamc > 0 ) hdphi1_mc->Fill(dphi);
				}
				if(getNJets() > 1){
					if(S->datamc == 0) hdphi2_da->Fill(dphi, singleMuPrescale());
					if(S->datamc > 0 ) hdphi2_mc->Fill(dphi);					
				}
			}
		}
		S->cleanUp();
	}
	
	printObject(hdphi1_da, "DR_1Jet_Data", "PEX");
	printObject(hdphi2_da, "DR_2Jet_Data", "PEX");
	// printObject(hdphi1_mc, "DPhi_1Jet_MC", "PEX");
	// printObject(hdphi2_mc, "DPhi_2Jet_MC", "PEX");
}

//____________________________________________________________________________
void SSDLPlotter::makePileUpPlots(bool write){
	fOutputSubDir = "PileUp/";

	TH1D *smu_nvertices, *dmu_nvertices, *sel_nvertices, *del_nvertices, *mue_nvertices;
	TH1D *mu_ntight, *mu_nloose, *mu_ratio, *el_ntight, *el_nloose, *el_ratio;

	if(!write){
		TFile *file = TFile::Open(fOutputDir + fOutputSubDir + "histos.root");
		smu_nvertices = (TH1D*)file->Get("smu_nvertices");
		dmu_nvertices = (TH1D*)file->Get("dmu_nvertices");
		sel_nvertices = (TH1D*)file->Get("sel_nvertices");
		del_nvertices = (TH1D*)file->Get("del_nvertices");
		mue_nvertices = (TH1D*)file->Get("mue_nvertices");
		mu_ratio      = (TH1D*)file->Get("mu_ratio");
		el_ratio      = (TH1D*)file->Get("el_ratio");
	}
	else{
		vector<gSample> samples;
		// samples.push_back(DoubleMu1);
		// samples.push_back(DoubleMu2);
		// samples.push_back(DoubleMu3);
		// samples.push_back(DoubleMu4);
		samples.push_back(DoubleMu5);
		// samples.push_back(MuEG1);
		// samples.push_back(MuEG2);
		// samples.push_back(MuEG3);
		// samples.push_back(MuEG4);
		samples.push_back(MuEG5);
		// samples.push_back(DoubleEle1);
		// samples.push_back(DoubleEle2);
		// samples.push_back(DoubleEle3);
		// samples.push_back(DoubleEle4);
		samples.push_back(DoubleEle5);

		smu_nvertices = new TH1D("smu_nvertices", "smu_nvertices", FRatioPlots::nbins[3], FRatioPlots::xmin[3], FRatioPlots::xmax[3]);
		dmu_nvertices = new TH1D("dmu_nvertices", "dmu_nvertices", FRatioPlots::nbins[3], FRatioPlots::xmin[3], FRatioPlots::xmax[3]);
		sel_nvertices = new TH1D("sel_nvertices", "sel_nvertices", FRatioPlots::nbins[3], FRatioPlots::xmin[3], FRatioPlots::xmax[3]);
		del_nvertices = new TH1D("del_nvertices", "del_nvertices", FRatioPlots::nbins[3], FRatioPlots::xmin[3], FRatioPlots::xmax[3]);
		mue_nvertices = new TH1D("mue_nvertices", "mue_nvertices", FRatioPlots::nbins[3], FRatioPlots::xmin[3], FRatioPlots::xmax[3]);

		mu_ntight = new TH1D("mu_ntight", "ntight", FRatioPlots::nbins[3], FRatioPlots::xmin[3], FRatioPlots::xmax[3]);
		mu_nloose = new TH1D("mu_nloose", "nloose", FRatioPlots::nbins[3], FRatioPlots::xmin[3], FRatioPlots::xmax[3]);
		mu_ratio  = new TH1D("mu_ratio",  "ratio",  FRatioPlots::nbins[3], FRatioPlots::xmin[3], FRatioPlots::xmax[3]);
		el_ntight = new TH1D("el_ntight", "ntight", FRatioPlots::nbins[3], FRatioPlots::xmin[3], FRatioPlots::xmax[3]);
		el_nloose = new TH1D("el_nloose", "nloose", FRatioPlots::nbins[3], FRatioPlots::xmin[3], FRatioPlots::xmax[3]);
		el_ratio  = new TH1D("el_ratio",  "ratio",  FRatioPlots::nbins[3], FRatioPlots::xmin[3], FRatioPlots::xmax[3]);

		mu_ntight    ->Sumw2();
		mu_nloose    ->Sumw2();
		mu_ratio     ->Sumw2();
		el_ntight    ->Sumw2();
		el_nloose    ->Sumw2();
		el_ratio     ->Sumw2();
		smu_nvertices->Sumw2();
		dmu_nvertices->Sumw2();
		sel_nvertices->Sumw2();
		del_nvertices->Sumw2();
		mue_nvertices->Sumw2();
		mu_ntight    ->SetXTitle("N_{Vertices}");
		mu_nloose    ->SetXTitle("N_{Vertices}");
		mu_ratio     ->SetXTitle("N_{Vertices}");
		el_ntight    ->SetXTitle("N_{Vertices}");
		el_nloose    ->SetXTitle("N_{Vertices}");
		el_ratio     ->SetXTitle("N_{Vertices}");
		smu_nvertices->SetXTitle("N_{Vertices}");
		dmu_nvertices->SetXTitle("N_{Vertices}");
		sel_nvertices->SetXTitle("N_{Vertices}");
		del_nvertices->SetXTitle("N_{Vertices}");
		mue_nvertices->SetXTitle("N_{Vertices}");
		mu_ntight    ->SetYTitle("N_{Events}");
		mu_nloose    ->SetYTitle("N_{Events}");
		mu_ratio     ->SetYTitle("N_{Events}");
		el_ntight    ->SetYTitle("N_{Events}");
		el_nloose    ->SetYTitle("N_{Events}");
		el_ratio     ->SetYTitle("N_{Events}");
		smu_nvertices->SetYTitle("N_{Events}");
		dmu_nvertices->SetYTitle("N_{Events}");
		sel_nvertices->SetYTitle("N_{Events}");
		del_nvertices->SetYTitle("N_{Events}");
		mue_nvertices->SetYTitle("N_{Events}");
		mu_ntight    ->GetYaxis()->SetTitleOffset(1.2);
		mu_nloose    ->GetYaxis()->SetTitleOffset(1.2);
		mu_ratio     ->GetYaxis()->SetTitleOffset(1.2);
		el_ntight    ->GetYaxis()->SetTitleOffset(1.2);
		el_nloose    ->GetYaxis()->SetTitleOffset(1.2);
		el_ratio     ->GetYaxis()->SetTitleOffset(1.2);
		smu_nvertices->GetYaxis()->SetTitleOffset(1.2);
		dmu_nvertices->GetYaxis()->SetTitleOffset(1.2);
		sel_nvertices->GetYaxis()->SetTitleOffset(1.2);
		del_nvertices->GetYaxis()->SetTitleOffset(1.2);
		mue_nvertices->GetYaxis()->SetTitleOffset(1.2);

		for(size_t i = 0; i < samples.size(); ++i){
			Sample *S = fSamples[samples[i]];
			fSample = S; // necessary for triggers to work properly
			fCurrentSample = samples[i];
			resetHypLeptons();
			fDoCounting = false;
		
			TTree *tree = S->getTree();
			tree->ResetBranchAddresses();
			if(S->datamc < 1) Init(tree);
			if(S->datamc > 0) InitMC(tree);
			for (Long64_t jentry=0; jentry<tree->GetEntriesFast();jentry++) {
				tree->GetEntry(jentry);
				printProgress(jentry, tree->GetEntriesFast(), S->sname);

				if(mumuSignalTrigger()) dmu_nvertices->Fill(NVrtx);
				if(elelSignalTrigger()) del_nvertices->Fill(NVrtx);
				if(elmuSignalTrigger()) mue_nvertices->Fill(NVrtx);
				if(singleMuTrigger()){
					smu_nvertices->Fill(NVrtx);
					if(isSigSupMuEvent()){
						if(isTightMuon(0)) mu_ntight->Fill(NVrtx);
						if(isLooseMuon(0)) mu_nloose->Fill(NVrtx);					
					}
				}
				resetHypLeptons();
				if(singleElTrigger()){
					sel_nvertices->Fill(NVrtx);
				 	if(isSigSupElEvent()){
						if(isTightElectron(0)) el_ntight->Fill(NVrtx);
						if(isLooseElectron(0)) el_nloose->Fill(NVrtx);
					}
				}
			}
			S->cleanUp();
		}
	
		mu_ratio->Divide(mu_ntight, mu_nloose, 1., 1., "B");
		mu_ratio->GetYaxis()->SetRangeUser(0., 0.4);
		el_ratio->Divide(el_ntight, el_nloose, 1., 1., "B");
		el_ratio->GetYaxis()->SetRangeUser(0., 0.4);

		dmu_nvertices->Scale(1./dmu_nvertices->Integral());
		smu_nvertices->Scale(1./smu_nvertices->Integral());
		del_nvertices->Scale(1./del_nvertices->Integral());
		sel_nvertices->Scale(1./sel_nvertices->Integral());
		mue_nvertices->Scale(1./mue_nvertices->Integral());

	}

	mu_ratio     ->SetYTitle("N_{Events} (Normalized)");
	el_ratio     ->SetYTitle("N_{Events} (Normalized)");
	smu_nvertices->SetYTitle("N_{Events} (Normalized)");
	dmu_nvertices->SetYTitle("N_{Events} (Normalized)");
	sel_nvertices->SetYTitle("N_{Events} (Normalized)");
	del_nvertices->SetYTitle("N_{Events} (Normalized)");
	mue_nvertices->SetYTitle("N_{Events} (Normalized)");

	mu_ratio     ->GetYaxis()->SetTitleOffset(1.3);
	el_ratio     ->GetYaxis()->SetTitleOffset(1.3);
	smu_nvertices->GetYaxis()->SetTitleOffset(1.3);
	dmu_nvertices->GetYaxis()->SetTitleOffset(1.3);
	sel_nvertices->GetYaxis()->SetTitleOffset(1.3);
	del_nvertices->GetYaxis()->SetTitleOffset(1.3);
	mue_nvertices->GetYaxis()->SetTitleOffset(1.3);

	// Color_t colors[6] = {1, 1, 1, 1, 1, 1};
	Color_t colors[5] = {31, 41, 51, 61, 81};
	// Color_t colors[6] = {1, 12, 39, 38, 32, 30};//, 29};
	Style_t styles[6] = {24, 25, 26, 27, 32, 30};//, 29};

	// dmu_nvertices->SetMarkerStyle(styles[0]);
	// smu_nvertices->SetMarkerStyle(styles[1]);
	// del_nvertices->SetMarkerStyle(styles[2]);
	// sel_nvertices->SetMarkerStyle(styles[3]);
	// mue_nvertices->SetMarkerStyle(styles[4]);
	// dmu_nvertices->SetMarkerColor(colors[0]);
	// smu_nvertices->SetMarkerColor(colors[1]);
	// del_nvertices->SetMarkerColor(colors[2]);
	// sel_nvertices->SetMarkerColor(colors[3]);
	// mue_nvertices->SetMarkerColor(colors[4]);
	dmu_nvertices->SetLineColor(colors[0]);
	smu_nvertices->SetLineColor(colors[1]);
	del_nvertices->SetLineColor(colors[2]);
	sel_nvertices->SetLineColor(colors[3]);
	mue_nvertices->SetLineColor(colors[4]);
	// dmu_nvertices->SetMarkerSize(1.5);
	// smu_nvertices->SetMarkerSize(1.5);
	// del_nvertices->SetMarkerSize(1.5);
	// sel_nvertices->SetMarkerSize(1.5);
	// mue_nvertices->SetMarkerSize(1.5);
	dmu_nvertices->SetLineWidth(2);
	smu_nvertices->SetLineWidth(2);
	del_nvertices->SetLineWidth(2);
	sel_nvertices->SetLineWidth(2);
	mue_nvertices->SetLineWidth(2);
	dmu_nvertices->SetFillStyle(0);
	smu_nvertices->SetFillStyle(0);
	del_nvertices->SetFillStyle(0);
	sel_nvertices->SetFillStyle(0);
	mue_nvertices->SetFillStyle(0);

	dmu_nvertices->GetYaxis()->SetRangeUser(0., 0.4);
	smu_nvertices->GetYaxis()->SetRangeUser(0., 0.4);
	del_nvertices->GetYaxis()->SetRangeUser(0., 0.4);
	sel_nvertices->GetYaxis()->SetRangeUser(0., 0.4);
	mue_nvertices->GetYaxis()->SetRangeUser(0., 0.4);

	mu_ratio->SetMarkerStyle(20);
	mu_ratio->SetMarkerColor(kBlue);
	mu_ratio->SetLineColor(kBlue);
	mu_ratio->SetMarkerSize(1.8);
	mu_ratio->SetLineWidth(2);
	el_ratio->SetMarkerStyle(21);
	el_ratio->SetMarkerColor(kRed);
	el_ratio->SetLineColor(kRed);
	el_ratio->SetMarkerSize(1.8);
	el_ratio->SetLineWidth(2);

	// printObject(mu_ratio, "MuRatio", "PEX");
	// printObject(el_ratio, "ElRatio", "PEX");
	// printObject(dmu_nvertices, "NVertices_DMu", "L");
	// printObject(smu_nvertices, "NVertices_SMu", "L");
	// printObject(del_nvertices, "NVertices_DEl", "L");
	// printObject(sel_nvertices, "NVertices_SEl", "L");
	// printObject(mue_nvertices, "NVertices_MuE", "L");

	// TLegend *leg = new TLegend(0.35,0.15,0.70,0.40);
	TLegend *leg = new TLegend(0.50,0.60,0.88,0.88);
	leg->AddEntry(dmu_nvertices, Form("DoubleMu Trig., Mean = %4.2f", dmu_nvertices->GetMean()),  "l");
	leg->AddEntry(smu_nvertices, Form("SingleMu Trig., Mean = %4.2f", smu_nvertices->GetMean()),  "l");
	leg->AddEntry(del_nvertices, Form("DoubleEle Trig., Mean = %4.2f", del_nvertices->GetMean()), "l");
	leg->AddEntry(sel_nvertices, Form("SingleEle Trig., Mean = %4.2f", sel_nvertices->GetMean()), "l");
	leg->AddEntry(mue_nvertices, Form("MuEle Trig., Mean = %4.2f", mue_nvertices->GetMean()),   "l");
	leg->AddEntry(mu_ratio, Form("TL Ratio (Muons)"), "p");
	leg->AddEntry(el_ratio, Form("TL Ratio (Electrons)"),  "p");
	leg->SetFillStyle(0);
	leg->SetTextFont(42);
	// leg->SetTextSize(0.03);
	leg->SetBorderSize(0);

	TCanvas *c_temp = new TCanvas("C_temp", "HT vs MET in Data vs MC", 0, 0, 800, 600);
	c_temp->cd();
	dmu_nvertices->Draw("axis");
	dmu_nvertices->DrawCopy("hist L same");
	smu_nvertices->DrawCopy("hist L same");
	del_nvertices->DrawCopy("hist L same");
	sel_nvertices->DrawCopy("hist L same");
	mue_nvertices->DrawCopy("hist L same");
	mu_ratio->DrawCopy("PE X0 same");
	el_ratio->DrawCopy("PE X0 same");
	TGaxis *axis = new TGaxis(18, 0, 18, 0.4, 0, 0.4, 510, "+L");
	axis->SetLabelFont(42);
	axis->SetTitleFont(42);
	axis->SetTitleOffset(1.2);
	axis->SetTitle("TL Ratio");
	axis->Draw();
	leg->Draw();
	drawTopLine();
	// fLatex->DrawLatex(0.10,0.92, fSamples[sample]->name);
	// Util::PrintNoEPS( c_temp, "PileUp", fOutputDir + fOutputSubDir, NULL);
	Util::PrintPDF(   c_temp, "PileUp", fOutputDir + fOutputSubDir);

	if(write){
		TFile *file = new TFile(fOutputDir + fOutputSubDir + "histos.root", "RECREATE");
		mu_ratio->Write(mu_ratio->GetName(), TObject::kWriteDelete);
		el_ratio->Write(el_ratio->GetName(), TObject::kWriteDelete);
		dmu_nvertices->Write(dmu_nvertices->GetName(), TObject::kWriteDelete);
		smu_nvertices->Write(smu_nvertices->GetName(), TObject::kWriteDelete);
		del_nvertices->Write(del_nvertices->GetName(), TObject::kWriteDelete);
		sel_nvertices->Write(sel_nvertices->GetName(), TObject::kWriteDelete);
		mue_nvertices->Write(mue_nvertices->GetName(), TObject::kWriteDelete);
		file->Close();
	}

	// Cleanup
	delete leg, c_temp;
	if(write) delete mu_ntight, mu_nloose, el_ntight, el_nloose, mu_ratio, el_ratio;
	else delete mu_ratio, el_ratio;
	delete dmu_nvertices, smu_nvertices, del_nvertices, sel_nvertices, mue_nvertices;
	fOutputSubDir = "";
}

//____________________________________________________________________________
void SSDLPlotter::makeNT012Plots(vector<int> mcsamples, gChannel chan, gRegion reg, gHiLoSwitch hilo){
	TString name;
	if(chan == Muon)     name = "MuMu";
	if(chan == Elec) name = "ElEl";
	if(chan == ElMu)      name = "ElMu";

	fOutputSubDir = name + "Predictions";

	THStack *hnt2_stack = new THStack(Form("%s_nt2_stack", name.Data()), "Observed Nt2");
	THStack *hnt1_stack = new THStack(Form("%s_nt1_stack", name.Data()), "Observed Nt1");
	THStack *hnt0_stack = new THStack(Form("%s_nt0_stack", name.Data()), "Observed Nt0");
	const unsigned int nmcsamples = mcsamples.size();
	TH1D* hnt2[nmcsamples];
	TH1D* hnt1[nmcsamples];
	TH1D* hnt0[nmcsamples];

	for(size_t i = 0; i < mcsamples.size(); ++i){
		Sample *S = fSamples[mcsamples[i]];
		float scale = fLumiNorm / S->lumi;
		Channel *cha;
		if(chan == Muon)     cha = &S->region[reg][hilo].mm;
		if(chan == Elec) cha = &S->region[reg][hilo].ee;
		if(chan == ElMu)      cha = &S->region[reg][hilo].em;
		hnt2[i] = (TH1D*)(cha->nt20_pt->ProjectionX())->Clone();
		hnt1[i] = (TH1D*)(cha->nt10_pt->ProjectionX())->Clone();
		hnt0[i] = (TH1D*)(cha->nt00_pt->ProjectionX())->Clone();
		hnt2[i]->Scale(scale);
		hnt1[i]->Scale(scale);
		hnt0[i]->Scale(scale);

		hnt2[i]->SetFillColor(S->color);
		hnt1[i]->SetFillColor(S->color);
		hnt0[i]->SetFillColor(S->color);

		hnt2[i] = normHistBW(hnt2[i], fBinWidthScale);
		hnt1[i] = normHistBW(hnt1[i], fBinWidthScale);
		hnt0[i] = normHistBW(hnt0[i], fBinWidthScale);

		hnt2_stack->Add(hnt2[i]);
		hnt1_stack->Add(hnt1[i]);
		hnt0_stack->Add(hnt0[i]);
	}

	hnt2_stack->Draw();
	hnt2_stack->GetXaxis()->SetTitle(convertVarName("MuPt[0]"));
	hnt1_stack->Draw();
	hnt1_stack->GetXaxis()->SetTitle(convertVarName("MuPt[0]"));
	hnt0_stack->Draw();
	hnt0_stack->GetXaxis()->SetTitle(convertVarName("MuPt[0]"));

	TCanvas *c2 = new TCanvas(Form("%s_Nt2", name.Data()), "Observed Nt2", 0, 0, 800, 600);
	TCanvas *c1 = new TCanvas(Form("%s_Nt1", name.Data()), "Observed Nt1", 0, 0, 800, 600);
	TCanvas *c0 = new TCanvas(Form("%s_Nt0", name.Data()), "Observed Nt0", 0, 0, 800, 600);

	// TLegend *leg = new TLegend(0.13,0.60,0.38,0.88);
	TLegend *leg = new TLegend(0.75,0.60,0.89,0.88);
	for(size_t i = 0; i < nmcsamples; ++i){
		int index = mcsamples[i];
		leg->AddEntry(hnt2[i], (fSamples[index]->sname).Data(), "f");
	}
	leg->SetFillStyle(0);
	leg->SetTextFont(42);
	leg->SetBorderSize(0);

	c2->cd();
	hnt2_stack->Draw("hist");
	leg->Draw();

	c1->cd();
	hnt1_stack->Draw("hist");
	leg->Draw();

	c0->cd();
	hnt0_stack->Draw("hist");
	leg->Draw();


	Util::PrintNoEPS(c2, name + "ObservedNt2", fOutputDir + fOutputSubDir, fOutputFile);
	Util::PrintNoEPS(c1, name + "ObservedNt1", fOutputDir + fOutputSubDir, fOutputFile);
	Util::PrintNoEPS(c0, name + "ObservedNt0", fOutputDir + fOutputSubDir, fOutputFile);
	delete c2, c1, c0;
	// for(size_t i = 0; i < nmcsamples; ++i){
	// 	delete hnt2[i];
	// 	delete hnt1[i];
	// 	delete hnt0[i];
	// }
	fOutputSubDir = "";
}
void SSDLPlotter::makeNT012Plots(gChannel chan, vector<int> mcsamples, bool(SSDLPlotter::*eventSelector)(int&, int&), TString tag){
	const bool read = false;

	THStack *hnt20_stack = new THStack("nt20_stack", "Observed Nt20");
	THStack *hnt10_stack = new THStack("nt10_stack", "Observed Nt10");
	THStack *hnt01_stack = new THStack("nt01_stack", "Observed Nt01");
	THStack *hnt00_stack = new THStack("nt00_stack", "Observed Nt00");
	const unsigned int nmcsamples = mcsamples.size();
	TH1D *hnt20[nmcsamples];
	TH1D *hnt10[nmcsamples];
	TH1D *hnt01[nmcsamples];
	TH1D *hnt00[nmcsamples];

	if(read){
		hnt20_stack = (THStack*)fOutputFile->Get(fOutputDir + "/nt20_stack");
		hnt10_stack = (THStack*)fOutputFile->Get(fOutputDir + "/nt10_stack");
		hnt01_stack = (THStack*)fOutputFile->Get(fOutputDir + "/nt01_stack");
		hnt00_stack = (THStack*)fOutputFile->Get(fOutputDir + "/nt00_stack");
		TList *list20 = hnt20_stack->GetHists();
		TList *list10 = hnt10_stack->GetHists();
		TList *list01 = hnt01_stack->GetHists();
		TList *list00 = hnt00_stack->GetHists();
		if(list20->GetSize() != nmcsamples || list10->GetSize() != nmcsamples || list00->GetSize() != nmcsamples) return;
		for(size_t i = 0; i < nmcsamples; ++i){
			int index = mcsamples[i];
			hnt20[i] = (TH1D*)list20->At(i);
			hnt20[i]->SetFillColor(fSamples[index]->color);
			hnt10[i] = (TH1D*)list10->At(i);
			hnt10[i]->SetFillColor(fSamples[index]->color);
			hnt01[i] = (TH1D*)list01->At(i);
			hnt01[i]->SetFillColor(fSamples[index]->color);
			hnt00[i] = (TH1D*)list00->At(i);
			hnt00[i]->SetFillColor(fSamples[index]->color);
		}
	}

	if(!read){
		TTree *tree = NULL;
		for(size_t i = 0; i < mcsamples.size(); ++i){
			int index = mcsamples[i];
			tree = fSamples[index]->getTree();
			hnt20[i] = new TH1D(Form("nt20_%s", fSamples[index]->sname.Data()), "Observed Nt20", getNPt2Bins(Muon), getPt2Bins(Muon));
			hnt10[i] = new TH1D(Form("nt10_%s", fSamples[index]->sname.Data()), "Observed Nt10", getNPt2Bins(Muon), getPt2Bins(Muon));
			hnt01[i] = new TH1D(Form("nt01_%s", fSamples[index]->sname.Data()), "Observed Nt01", getNPt2Bins(Muon), getPt2Bins(Muon));
			hnt00[i] = new TH1D(Form("nt00_%s", fSamples[index]->sname.Data()), "Observed Nt00", getNPt2Bins(Muon), getPt2Bins(Muon));
			hnt20[i]->SetFillColor(fSamples[index]->color);
			hnt10[i]->SetFillColor(fSamples[index]->color);
			hnt01[i]->SetFillColor(fSamples[index]->color);
			hnt00[i]->SetFillColor(fSamples[index]->color);
			hnt20[i]->Sumw2();
			hnt10[i]->Sumw2();
			hnt01[i]->Sumw2();
			hnt00[i]->Sumw2();
			float scale = fLumiNorm / fSamples[index]->lumi;
			if(fSamples[index]->datamc == 0) scale = 1;
			tree->ResetBranchAddresses();
			Init(tree);
			if (fChain == 0) return;
			Long64_t nentries = fChain->GetEntriesFast();
			Long64_t nbytes = 0, nb = 0;
			for (Long64_t jentry=0; jentry<nentries;jentry++) {
				Long64_t ientry = LoadTree(jentry);
				if (ientry < 0) break;
				if(fVerbose > 1) printProgress(jentry, nentries, fSamples[index]->name);
				nb = fChain->GetEntry(jentry);   nbytes += nb;

				int ind1(-1), ind2(-1);
				if((*this.*eventSelector)(ind1, ind2) == false) continue;

				if(chan == Muon){
					if( isTightMuon(ind1) &&  isTightMuon(ind2)) hnt20[i]->Fill(MuPt[ind1], scale);
					if( isTightMuon(ind1) && !isTightMuon(ind2)) hnt10[i]->Fill(MuPt[ind1], scale);
					if(!isTightMuon(ind1) &&  isTightMuon(ind2)) hnt10[i]->Fill(MuPt[ind2], scale);
					if(!isTightMuon(ind1) && !isTightMuon(ind2)) hnt00[i]->Fill(MuPt[ind1], scale);
				}
				if(chan == Elec){
					if( isTightElectron(ind1) &&  isTightElectron(ind2)) hnt20[i]->Fill(ElPt[ind1], scale);
					if( isTightElectron(ind1) && !isTightElectron(ind2)) hnt10[i]->Fill(ElPt[ind1], scale);
					if(!isTightElectron(ind1) &&  isTightElectron(ind2)) hnt10[i]->Fill(ElPt[ind2], scale);
					if(!isTightElectron(ind1) && !isTightElectron(ind2)) hnt00[i]->Fill(ElPt[ind1], scale);
				}
				if(chan == ElMu){
					if( isTightMuon(ind1) &&  isTightElectron(ind2)) hnt20[i]->Fill(MuPt[ind1], scale);
					if( isTightMuon(ind1) && !isTightElectron(ind2)) hnt10[i]->Fill(MuPt[ind1], scale);
					if(!isTightMuon(ind1) &&  isTightElectron(ind2)) hnt01[i]->Fill(MuPt[ind1], scale);
					if(!isTightMuon(ind1) && !isTightElectron(ind2)) hnt00[i]->Fill(MuPt[ind1], scale);
				}
			}
			hnt20_stack->Add(hnt20[i]);
			hnt10_stack->Add(hnt10[i]);
			hnt01_stack->Add(hnt01[i]);
			hnt00_stack->Add(hnt00[i]);
			if(fVerbose > 1) cout << endl;
			fSamples[index]->cleanUp();
		}
	}

	hnt20_stack->Draw();
	hnt20_stack->GetXaxis()->SetTitle(convertVarName("MuPt[0]"));
	hnt10_stack->Draw();
	hnt10_stack->GetXaxis()->SetTitle(convertVarName("MuPt[0]"));
	hnt01_stack->Draw();
	hnt01_stack->GetXaxis()->SetTitle(convertVarName("MuPt[0]"));
	hnt00_stack->Draw();
	hnt00_stack->GetXaxis()->SetTitle(convertVarName("MuPt[0]"));

	TCanvas *c20 = new TCanvas("Nt20", "Observed Nt20", 0, 0, 800, 600);
	TCanvas *c10 = new TCanvas("Nt10", "Observed Nt10", 0, 0, 800, 600);
	TCanvas *c01 = new TCanvas("Nt01", "Observed Nt01", 0, 0, 800, 600);
	TCanvas *c00 = new TCanvas("Nt00", "Observed Nt00", 0, 0, 800, 600);

	// TLegend *leg = new TLegend(0.13,0.60,0.38,0.88);
	TLegend *leg = new TLegend(0.75,0.60,0.89,0.88);
	for(size_t i = 0; i < nmcsamples; ++i){
		int index = mcsamples[i];
		leg->AddEntry(hnt20[i], fSamples[index]->sname.Data(), "f");
	}
	leg->SetFillStyle(0);
	leg->SetTextFont(42);
	leg->SetBorderSize(0);

	c20->cd();
	hnt20_stack->Draw("hist");
	leg->Draw();

	c10->cd();
	hnt10_stack->Draw("hist");
	leg->Draw();

	c01->cd();
	hnt01_stack->Draw("hist");
	leg->Draw();

	c00->cd();
	hnt00_stack->Draw("hist");
	leg->Draw();

	Util::PrintNoEPS(c20, tag + "ObservedNt20", fOutputDir + fOutputSubDir, fOutputFile);
	Util::PrintNoEPS(c10, tag + "ObservedNt10", fOutputDir + fOutputSubDir, fOutputFile);
	Util::PrintNoEPS(c01, tag + "ObservedNt01", fOutputDir + fOutputSubDir, fOutputFile);
	Util::PrintNoEPS(c00, tag + "ObservedNt00", fOutputDir + fOutputSubDir, fOutputFile);
}
//____________________________________________________________________________
void SSDLPlotter::makeRelIsoTTSigPlots(){
	char cmd[100];
    sprintf(cmd,"mkdir -p %s%s", fOutputDir.Data(), fOutputSubDir.Data());
    system(cmd);

	TLatex *lat = new TLatex();
	lat->SetNDC(kTRUE);
	lat->SetTextColor(kBlack);
	lat->SetTextSize(0.04);

	TLatex *el_lat = new TLatex();
	el_lat->SetNDC(kTRUE);
	el_lat->SetTextColor(kBlack);
	el_lat->SetTextSize(0.04);

	// Create histograms
	TH1D * hiso_ttbar  = new TH1D("MuIsoTTbar_"   , "Muon Isolation in TTbar for " , IsoPlots::nbins[0], 0., 1.);
	TH1D * hiso_signal = new TH1D("MuIsoSignal_"  , "Muon Isolation in Signal for ", IsoPlots::nbins[0], 0., 1.);
	hiso_ttbar ->Sumw2();
	hiso_signal->Sumw2();
	TH1D * el_hiso_ttbar  = new TH1D("ElIsoTTbar_"   , "Electron Isolation in TTbar for " , IsoPlots::nbins[0], 0., 0.6);
	TH1D * el_hiso_signal = new TH1D("ElIsoSignal_"  , "Electron Isolation in Signal for ", IsoPlots::nbins[0], 0., 0.6);
	el_hiso_ttbar ->Sumw2();
	el_hiso_signal->Sumw2();

	////////////////////////////////////////////////////
	// Fill ttbar histos
	TTree *ttbar_tree = fSamples[TTJets]->getTree();
	ttbar_tree->ResetBranchAddresses();
	InitMC(ttbar_tree);

	// Event loop muons and electrons
	if (fChain == 0) return;
	Long64_t tt_nentries = fChain->GetEntriesFast();
	for (Long64_t jentry=0; jentry<tt_nentries;jentry++) {
		printProgress(jentry, tt_nentries, fSamples[TTJets]->name);

		ttbar_tree->GetEntry(jentry);
		//Long64_t ientry = LoadTree(jentry);
		//if (ientry < 0) break;
		
		int muind1(-1), muind2(-1);
		if(hasLooseMuons(muind1, muind2) > 0) {
			// Common event selections
			if(passesJet50Cut()) { // make trigger 100% efficient
				// Common object selections
				if( isLooseMuon(muind1) ) { //&& isLooseMuon(muind2)){
					if(MuPt[muind1] > fC_minMu2pt) {
					// match muons to susy or vector boson
						if( (MuGenMType[muind1] == 9  || MuGenGMType[muind1] == 9 ) || (MuGenMType[muind1]  == 4 && abs(MuGenMID[muind1]) != 21 && abs(MuGenMID[muind1]) != 22) ) {
							hiso_ttbar->Fill(MuIso[muind1]);
						}
						//if( (MuGenMType[muind2] == 9  || MuGenGMType[muind2] == 9 ) || (MuGenMType[muind2]  == 4 && abs(MuGenMID[muind2]) != 21 && abs(MuGenMID[muind2]) != 22) ) {
						//	hiso_ttbar->Fill(MuIso[muind2]);
						//}
					}
				}
			}
		}
		int elind1(-1), elind2(-1);
		if(hasLooseElectrons(elind1, elind2) > 0) {
			// Common event selections
			if(passesJet50Cut()) { // make trigger 100% efficient
				// Common object selections
				if(isLooseElectron(elind1) ) { //&& isLooseElectron(elind2)) {
					if(ElPt[elind1] > fC_minEl2pt) {
						// match electrons to susy particle or vector boson
						if( (ElGenMType[elind1] == 9  || ElGenGMType[elind1] == 9 ) || (ElGenMType[elind1]  == 4 && abs(ElGenMID[elind1]) != 21 && abs(ElGenMID[elind1]) != 22) ) {
							el_hiso_ttbar->Fill(ElRelIso[elind1]);
						}
						//if( (ElGenMType[elind2] == 9  || ElGenGMType[elind2] == 9 ) || (ElGenMType[elind2]  == 4 && abs(ElGenMID[elind2]) != 21 && abs(ElGenMID[elind2]) != 22) ) {
						//	el_hiso_ttbar->Fill(ElRelIso[elind2]);
						//}
					}
				}
			}
		}
	} // end loop over all events
	fSamples[TTJets]->cleanUp();
	cout << endl;
	////////////////////////////////////////////////////

	////////////////////////////////////////////////////
	// Fill signal histos
	TTree *signal_tree = fSamples[LM6]->getTree();
	signal_tree->ResetBranchAddresses();
	InitMC(signal_tree);
	// Event loop

	// Event loop muons and electrons
	if (fChain == 0) return;
	Long64_t sig_nentries = fChain->GetEntriesFast();
	for (Long64_t jentry=0; jentry<sig_nentries;jentry++) {
		printProgress(jentry, sig_nentries, fSamples[LM6]->name);

		signal_tree->GetEntry(jentry);
		//Long64_t ientry = LoadTree(jentry);
		//if (ientry < 0) break;
		
		int muind1(-1), muind2(-1);
		if(hasLooseMuons(muind1, muind2) > 0) {
			// Common event selections
			if(passesJet50Cut()) { // make trigger 100% efficient
				// Common object selections
				if( isLooseMuon(muind1) ) { //&& isLooseMuon(muind2)){
					if(MuPt[muind1] > fC_minMu2pt) {
					// match muons to susy or vector boson
						if( (MuGenMType[muind1] == 9  || MuGenGMType[muind1] == 9 ) || (MuGenMType[muind1]  == 4 && abs(MuGenMID[muind1]) != 21 && abs(MuGenMID[muind1]) != 22) ) {
							hiso_signal->Fill(MuIso[muind1]);
						}
						//if( (MuGenMType[muind2] == 9  || MuGenGMType[muind2] == 9 ) || (MuGenMType[muind2]  == 4 && abs(MuGenMID[muind2]) != 21 && abs(MuGenMID[muind2]) != 22) ) {
						//	hiso_signal->Fill(MuIso[muind2]);
						//}
					}
				}
			}
		}
		int elind1(-1), elind2(-1);
		if(hasLooseElectrons(elind1, elind2) > 0) {
			// Common event selections
			if(passesJet50Cut()) { // make trigger 100% efficient
				// Common object selections
				if(isLooseElectron(elind1) ){ //&& isLooseElectron(elind2)) {
					if(ElPt[elind1] > fC_minEl2pt) {
						// match electrons to susy particle or vector boson
						if( (ElGenMType[elind1] == 9  || ElGenGMType[elind1] == 9 ) || (ElGenMType[elind1]  == 4 && abs(ElGenMID[elind1]) != 21 && abs(ElGenMID[elind1]) != 22) ) {
							el_hiso_signal->Fill(ElRelIso[elind1]);
						}
						//if( (ElGenMType[elind2] == 9  || ElGenGMType[elind2] == 9 ) || (ElGenMType[elind2]  == 4 && abs(ElGenMID[elind2]) != 21 && abs(ElGenMID[elind2]) != 22) ) {
						//	el_hiso_signal->Fill(ElRelIso[elind2]);
						//}
					}
				}
			}
		}
	} // end loop over all events
	fSamples[LM6]->cleanUp();
	cout << endl;
	////////////////////////////////////////////////////
	
	//--------------------------------------------------------------------------------------
	// Format and calculate muon histos
	hiso_ttbar->SetXTitle("rel Iso muons");
	hiso_ttbar->SetLineWidth(2);
	//hiso_ttbar->SetFillColor(kWhite);
	//hiso_ttbar->SetLineColor(kRed);
	hiso_ttbar->SetMarkerStyle(20);
	hiso_ttbar->SetMarkerColor(kBlue-3);
	hiso_ttbar->SetMarkerSize(1.15);
		
	//hiso_signal[i]->SetXTitle(convertVarName("MuIso[0]"));
	hiso_signal->SetLineWidth(2);
	//hiso_signal->SetLineColor(kBlue);
	hiso_signal->SetMarkerStyle(21);
	hiso_signal->SetMarkerColor(kRed+1);
	hiso_signal->SetMarkerSize(1.15);
		
	//double max1 = hiso_mc_s[i]->GetMaximum();
	//double max2 = hiso_data[i]->GetMaximum();
	//double max = max1>max2?max1:max2;
	//hiso_mc_s[i]->SetMaximum(1.5*max);
	//hiso_data[i]->SetMaximum(1.5*max);

	int bin0   = hiso_ttbar->FindBin(0.0);
	int bin015 = hiso_ttbar->FindBin(0.15) - 1; // bins start at lower edge...
	int bin1   = hiso_ttbar->FindBin(1.0)  - 1;
	float ratio_signal = hiso_signal ->Integral(bin0, bin015) / hiso_signal ->Integral(bin0, bin1);
	float ratio_ttbar  = hiso_ttbar  ->Integral(bin0, bin015) / hiso_ttbar  ->Integral(bin0, bin1);

	
	hiso_signal ->Scale(1/hiso_signal->Integral());
	hiso_ttbar  ->Scale(1/hiso_ttbar ->Integral());
	

	TCanvas *c_temp = new TCanvas("MuIso" , "Muon Isolation in TTJets and Signal", 0, 0, 800, 600);
	c_temp->cd();

	TLegend *leg = new TLegend(0.30,0.65,0.55,0.88);
	// TLegend *leg = new TLegend(0.75,0.60,0.89,0.88);
	leg->AddEntry(hiso_signal, "Signal","p");
	leg->AddEntry(hiso_ttbar , "TTbar","p");
	leg->SetFillStyle(0);
	leg->SetTextFont(42);
	leg->SetBorderSize(0);

	gPad->SetLogy();
	hiso_ttbar->Draw("pe");
	hiso_signal->Draw("same");
	leg->Draw();
	lat->DrawLatex(0.75,0.92, Form("L_{int.} = %2.1f fb^{-1}", fLumiNorm/1000.));
	lat->SetTextColor(kRed+1);
	lat->DrawLatex(0.70,0.70, Form("R^{T/L}_{Signal} = %4.2f", ratio_signal));
	lat->SetTextColor(kBlue-3);
	lat->DrawLatex(0.70,0.60, Form("R^{T/L}_{TTbar} = %4.2f", ratio_ttbar));
	

	Util::PrintPDF(c_temp, "MuIso" , fOutputDir + fOutputSubDir);

	//--------------------------------------------------------------------------------------
	// Format and calculate electron histos
	el_hiso_ttbar->SetXTitle("rel Iso electrons");
	el_hiso_ttbar->SetLineWidth(2);
	//hiso_ttbar->SetFillColor(kWhite);
	//hiso_ttbar->SetLineColor(kRed);
	el_hiso_ttbar->SetMarkerStyle(20);
	el_hiso_ttbar->SetMarkerColor(kBlue-3);
	el_hiso_ttbar->SetMarkerSize(1.15);

	//hiso_signal[i]->SetXTitle(convertVarName("MuIso[0]"));
	el_hiso_signal->SetLineWidth(2);
	//hiso_signal->SetLineColor(kBlue);
	el_hiso_signal->SetMarkerStyle(21);
	el_hiso_signal->SetMarkerColor(kRed+1);
	el_hiso_signal->SetMarkerSize(1.15);
		
	//double max1 = el_hiso_mc_s[i]->GetMaximum();
	//double max2 = el_hiso_data[i]->GetMaximum();
	//double max = max1>max2?max1:max2;
	//el_hiso_mc_s[i]->SetMaximum(1.5*max);
	//el_hiso_data[i]->SetMaximum(1.5*max);

	int el_bin0   = el_hiso_ttbar->FindBin(0.0);
	int el_bin015 = el_hiso_ttbar->FindBin(0.15) - 1; // bins start at lower edge...
	//int el_bin1   = el_hiso_ttbar->FindBin(1.0)  - 1;
	int el_bin1   = el_hiso_ttbar->FindBin(0.6)  - 1;
	float el_ratio_signal = el_hiso_signal ->Integral(el_bin0, el_bin015) / el_hiso_signal ->Integral(el_bin0, el_bin1);
	float el_ratio_ttbar  = el_hiso_ttbar  ->Integral(el_bin0, el_bin015) / el_hiso_ttbar  ->Integral(el_bin0, el_bin1);

	
	el_hiso_signal ->Scale(1/el_hiso_signal->Integral());
	el_hiso_ttbar  ->Scale(1/el_hiso_ttbar ->Integral());
	

	TCanvas *el_temp = new TCanvas("ElIso" , "Electron Isolation in TTJets and Signal", 0, 0, 800, 600);
	el_temp->cd();

	TLegend *el_leg = new TLegend(0.30,0.65,0.55,0.88);
	el_leg->AddEntry(el_hiso_signal, "Signal","p");
	el_leg->AddEntry(el_hiso_ttbar , "TTbar","p");
	el_leg->SetFillStyle(0);
	el_leg->SetTextFont(42);
	el_leg->SetBorderSize(0);

	gPad->SetLogy();
	el_hiso_ttbar->Draw("");
	el_hiso_signal->Draw("same");
	el_leg->Draw();
	el_lat->DrawLatex(0.75,0.92, Form("L_{int.} = %2.1f fb^{-1}", fLumiNorm/1000.));
	el_lat->SetTextColor(kRed+1);
	el_lat->DrawLatex(0.70,0.70, Form("R^{T/L}_{Signal} = %4.2f", el_ratio_signal));
	el_lat->SetTextColor(kBlue-3);
	el_lat->DrawLatex(0.70,0.60, Form("R^{T/L}_{TTbar} = %4.2f", el_ratio_ttbar));
	

	// Util::PrintNoEPS(c_temp, "MuIso" + IsoPlots::sel_name[i], fOutputDir + fOutputSubDir, NULL);
	Util::PrintPDF(el_temp, "ElIso" , fOutputDir + fOutputSubDir);

}

//____________________________________________________________________________
void SSDLPlotter::makeMuIsolationPlots(bool dottbar){
	char cmd[100];
    sprintf(cmd,"mkdir -p %s%s", fOutputDir.Data(), fOutputSubDir.Data());
    system(cmd);

	TH1D    *hiso_data [gNSels];
	TH1D    *hiso_mc   [gNSels];
	TH1D	*hiso_ttbar[gNSels];
	TH1D    *hiso_qcd  [gNSels];
	TH1D    *hiso_ttj  [gNSels];
	TH1D    *hiso_ewk  [gNSels];
	TH1D    *hiso_rare [gNSels];
	TH1D    *hiso_db   [gNSels];
	THStack *hiso_mc_s [gNSels];

	TH1D    *hiso_data_pt [gNSels][gNMuPt2bins];
	TH1D    *hiso_mc_pt   [gNSels][gNMuPt2bins];
	TH1D    *hiso_ttbar_pt[gNSels][gNMuPt2bins];
	TH1D    *hiso_qcd_pt  [gNSels][gNMuPt2bins];
	TH1D    *hiso_ttj_pt  [gNSels][gNMuPt2bins];
	TH1D    *hiso_ewk_pt  [gNSels][gNMuPt2bins];
	TH1D    *hiso_rare_pt [gNSels][gNMuPt2bins];
	TH1D    *hiso_db_pt   [gNSels][gNMuPt2bins];
	THStack *hiso_mc_pt_s [gNSels][gNMuPt2bins];

	TH1D    *hiso_data_nv [gNSels][gNNVrtxBins];
	TH1D    *hiso_mc_nv   [gNSels][gNNVrtxBins];
	TH1D    *hiso_ttbar_nv[gNSels][gNNVrtxBins];
	TH1D    *hiso_qcd_nv  [gNSels][gNNVrtxBins];
	TH1D    *hiso_ttj_nv  [gNSels][gNNVrtxBins];
	TH1D    *hiso_ewk_nv  [gNSels][gNNVrtxBins];
	TH1D    *hiso_rare_nv [gNSels][gNNVrtxBins];
	TH1D    *hiso_db_nv   [gNSels][gNNVrtxBins];
	THStack *hiso_mc_nv_s [gNSels][gNNVrtxBins];

	TLatex *lat = new TLatex();
	lat->SetNDC(kTRUE);
	lat->SetTextColor(kBlack);
	lat->SetTextSize(0.04);

	// Create histograms
	for(size_t i = 0; i < gNSels; ++i){
		hiso_data [i] = new TH1D("MuIsoData_"          + IsoPlots::sel_name[i], "Muon Isolation in Data for "  + IsoPlots::sel_name[i], IsoPlots::nbins[i], 0., 1.0);
		hiso_mc   [i] = new TH1D("MuIsoMC_"            + IsoPlots::sel_name[i], "Muon Isolation in MC for "    + IsoPlots::sel_name[i], IsoPlots::nbins[i], 0., 1.0);
		hiso_ttbar[i] = new TH1D("MuIsoTTbar_"         + IsoPlots::sel_name[i], "Muon Isolation in TTbar for " + IsoPlots::sel_name[i], IsoPlots::nbins[i], 0., 1.0);
		hiso_qcd  [i] = new TH1D("MuIsoQCD_"           + IsoPlots::sel_name[i], "Muon Isolation in QCD for "   + IsoPlots::sel_name[i], IsoPlots::nbins[i], 0., 1.0);
		hiso_ttj  [i] = new TH1D("MuIsoTop_"           + IsoPlots::sel_name[i], "Muon Isolation in Top for "   + IsoPlots::sel_name[i], IsoPlots::nbins[i], 0., 1.0);
		hiso_ewk  [i] = new TH1D("MuIsoEWK_"           + IsoPlots::sel_name[i], "Muon Isolation in EWK for "   + IsoPlots::sel_name[i], IsoPlots::nbins[i], 0., 1.0);
		hiso_rare [i] = new TH1D("MuIsorare_"          + IsoPlots::sel_name[i], "Muon Isolation in rare for "  + IsoPlots::sel_name[i], IsoPlots::nbins[i], 0., 1.0);
		hiso_db   [i] = new TH1D("MuIsoDB_"            + IsoPlots::sel_name[i], "Muon Isolation in DB for "    + IsoPlots::sel_name[i], IsoPlots::nbins[i], 0., 1.0);
		hiso_mc_s[i]  = new THStack("MuIsoMC_stacked_" + IsoPlots::sel_name[i], "Muon Isolation in MC for "    + IsoPlots::sel_name[i]);
		hiso_data [i] ->Sumw2();
		hiso_mc   [i] ->Sumw2();
		hiso_ttbar[i] ->Sumw2();
		hiso_qcd  [i] ->Sumw2();
		hiso_ttj  [i] ->Sumw2();
		hiso_ewk  [i] ->Sumw2();
		hiso_rare [i] ->Sumw2();
		hiso_db   [i] ->Sumw2();

		for(int k = 0; k < gNMuPt2bins; ++k){
			hiso_data_pt [i][k] = new TH1D(Form("MuIsoData_%s_Pt%d"         , IsoPlots::sel_name[i].Data(), k), "Muon Isolation in Data for "  + IsoPlots::sel_name[i], IsoPlots::nbins[i], 0., 1.0);
			hiso_mc_pt   [i][k] = new TH1D(Form("MuIsoMC_%s_Pt%d"           , IsoPlots::sel_name[i].Data(), k), "Muon Isolation in MC for "    + IsoPlots::sel_name[i], IsoPlots::nbins[i], 0., 1.0);
			hiso_ttbar_pt[i][k] = new TH1D(Form("MuIsoTTbar_%s_Pt%d"        , IsoPlots::sel_name[i].Data(), k), "Muon Isolation in TTbar for " + IsoPlots::sel_name[i], IsoPlots::nbins[i], 0., 1.0);
			hiso_qcd_pt  [i][k] = new TH1D(Form("MuIsoQCD_%s_Pt%d"          , IsoPlots::sel_name[i].Data(), k), "Muon Isolation in QCD for "   + IsoPlots::sel_name[i], IsoPlots::nbins[i], 0., 1.0);
			hiso_ttj_pt  [i][k] = new TH1D(Form("MuIsoTop_%s_Pt%d"          , IsoPlots::sel_name[i].Data(), k), "Muon Isolation in Top for "   + IsoPlots::sel_name[i], IsoPlots::nbins[i], 0., 1.0);
			hiso_ewk_pt  [i][k] = new TH1D(Form("MuIsoEWK_%s_Pt%d"          , IsoPlots::sel_name[i].Data(), k), "Muon Isolation in EWK for "   + IsoPlots::sel_name[i], IsoPlots::nbins[i], 0., 1.0);
			hiso_rare_pt [i][k] = new TH1D(Form("MuIsorare_%s_Pt%d"         , IsoPlots::sel_name[i].Data(), k), "Muon Isolation in rare for "  + IsoPlots::sel_name[i], IsoPlots::nbins[i], 0., 1.0);
			hiso_db_pt   [i][k] = new TH1D(Form("MuIsoDB_%s_Pt%d"           , IsoPlots::sel_name[i].Data(), k), "Muon Isolation in DB for "    + IsoPlots::sel_name[i], IsoPlots::nbins[i], 0., 1.0);
			hiso_mc_pt_s [i][k] = new THStack(Form("MuIsoMC_stacked_%s_Pt%d", IsoPlots::sel_name[i].Data(), k), "Muon Isolation in MC for "    + IsoPlots::sel_name[i]);
			hiso_data_pt [i][k]->Sumw2();
			hiso_mc_pt   [i][k]->Sumw2();
			hiso_ttbar_pt[i][k]->Sumw2();
			hiso_qcd_pt  [i][k]->Sumw2();
			hiso_ttj_pt  [i][k]->Sumw2();
			hiso_ewk_pt  [i][k]->Sumw2();
			hiso_rare_pt [i][k]->Sumw2();
			hiso_db_pt   [i][k]->Sumw2();
		}
		for(int k = 0; k < gNNVrtxBins; ++k){
			hiso_data_nv [i][k] = new TH1D(Form("MuIsoData_%s_NVtx%d"         , IsoPlots::sel_name[i].Data(), k), "Muon Isolation in Data for "  + IsoPlots::sel_name[i], IsoPlots::nbins[i], 0., 1.0);
			hiso_mc_nv   [i][k] = new TH1D(Form("MuIsoMC_%s_NVtx%d"           , IsoPlots::sel_name[i].Data(), k), "Muon Isolation in MC for "    + IsoPlots::sel_name[i], IsoPlots::nbins[i], 0., 1.0);
			hiso_ttbar_nv[i][k] = new TH1D(Form("MuIsoTTbar_%s_NVtx%d"        , IsoPlots::sel_name[i].Data(), k), "Muon Isolation in TTbar for " + IsoPlots::sel_name[i], IsoPlots::nbins[i], 0., 1.0);
			hiso_qcd_nv  [i][k] = new TH1D(Form("MuIsoQCD_%s_NVtx%d"          , IsoPlots::sel_name[i].Data(), k), "Muon Isolation in QCD for "   + IsoPlots::sel_name[i], IsoPlots::nbins[i], 0., 1.0);
			hiso_ttj_nv  [i][k] = new TH1D(Form("MuIsoTop_%s_NVtx%d"          , IsoPlots::sel_name[i].Data(), k), "Muon Isolation in Top for "   + IsoPlots::sel_name[i], IsoPlots::nbins[i], 0., 1.0);
			hiso_ewk_nv  [i][k] = new TH1D(Form("MuIsoEWK_%s_NVtx%d"          , IsoPlots::sel_name[i].Data(), k), "Muon Isolation in EWK for "   + IsoPlots::sel_name[i], IsoPlots::nbins[i], 0., 1.0);
			hiso_rare_nv [i][k] = new TH1D(Form("MuIsorare_%s_NVtx%d"         , IsoPlots::sel_name[i].Data(), k), "Muon Isolation in rare for "  + IsoPlots::sel_name[i], IsoPlots::nbins[i], 0., 1.0);
			hiso_db_nv   [i][k] = new TH1D(Form("MuIsoDB_%s_NVtx%d"           , IsoPlots::sel_name[i].Data(), k), "Muon Isolation in DB for "    + IsoPlots::sel_name[i], IsoPlots::nbins[i], 0., 1.0);
			hiso_mc_nv_s [i][k] = new THStack(Form("MuIsoMC_stacked_%s_NVtx%d", IsoPlots::sel_name[i].Data(), k), "Muon Isolation in MC for "    + IsoPlots::sel_name[i]);
			hiso_data_nv [i][k]->Sumw2();
			hiso_mc_nv   [i][k]->Sumw2();
			hiso_ttbar_nv[i][k]->Sumw2();
			hiso_qcd_nv  [i][k]->Sumw2();
			hiso_ttj_nv  [i][k]->Sumw2();
			hiso_ewk_nv  [i][k]->Sumw2();
			hiso_rare_nv [i][k]->Sumw2();
			hiso_db_nv   [i][k]->Sumw2();
		}
	}

	for(size_t i = 0; i < gNSels; ++i){
		hiso_qcd [i]->SetFillColor(kYellow-4);
		hiso_db  [i]->SetFillColor(kSpring-9);
		hiso_ewk [i]->SetFillColor(kGreen +1);
		hiso_ttj [i]->SetFillColor(kAzure-5);
		hiso_rare[i]->SetFillColor(kAzure+8);

		for(int k = 0; k < gNMuPt2bins; ++k){
			hiso_qcd_pt [i][k]->SetFillColor(kYellow-4);
			hiso_db_pt  [i][k]->SetFillColor(kSpring-9);
			hiso_ewk_pt [i][k]->SetFillColor(kGreen +1);
			hiso_ttj_pt [i][k]->SetFillColor(kAzure-5);
			hiso_rare_pt[i][k]->SetFillColor(kAzure+8);
		}

		for(int k = 0; k < gNNVrtxBins; ++k){
			hiso_qcd_nv [i][k]->SetFillColor(kYellow-4);
			hiso_db_nv  [i][k]->SetFillColor(kSpring-9);
			hiso_ewk_nv [i][k]->SetFillColor(kGreen +1);
			hiso_ttj_nv [i][k]->SetFillColor(kAzure-5);
			hiso_rare_nv[i][k]->SetFillColor(kAzure+8);
		}
	}

	if(dottbar){
		
		////////////////////////////////////////////////////
		// Fill ttbar histos
		// Sample loop
		TTree *tree = fSamples[TTJets]->getTree();

		// Event loop
		tree->ResetBranchAddresses();
		InitMC(tree);

		if (fChain == 0) return;
		Long64_t nentries = fChain->GetEntriesFast();
		Long64_t nbytes = 0, nb = 0;
		for (Long64_t jentry=0; jentry<nentries;jentry++) {
			printProgress(jentry, nentries, fSamples[TTJets]->name);

			Long64_t ientry = LoadTree(jentry);
			if (ientry < 0) break;
			nb = fChain->GetEntry(jentry);   nbytes += nb;
		
			int muind1(-1), muind2(-1);
			if(hasLooseMuons(muind1, muind2) < 1) continue;

			// Common event selections
			if(!passesJet50Cut()) continue; // make trigger 100% efficient

			// Common object selections
			if(!isLooseMuon(muind1)) continue;
			if(MuPt[muind1] < fC_minMu2pt) continue;
			if(MuPt[muind1] > gMuPt2bins[gNMuPt2bins]) continue;

			// Select genmatched fake muons
			if(isPromptMuon(muind1)) continue;

			////////////////////////////////////////////////////
			// MOST LOOSE SELECTION
			hiso_ttbar[0]->Fill(MuIso[muind1]);
			for(size_t k = 0; k < gNMuPt2bins; ++k){
				if(MuPt[muind1] < gMuPt2bins[k]) continue;
				if(MuPt[muind1] > gMuPt2bins[k+1]) continue;
				hiso_ttbar_pt[0][k]->Fill(MuIso[muind1]);
			}
			for(size_t k = 0; k < gNNVrtxBins; ++k){
				if(NVrtx < gNVrtxBins[k]) continue;
				if(NVrtx > gNVrtxBins[k+1]) continue;
				hiso_ttbar_nv[0][k]->Fill(MuIso[muind1]);
			}

			////////////////////////////////////////////////////
			// SIGNAL SUPPRESSED SELECTION
			if(isSigSupMuEvent()){
				hiso_ttbar[1]->Fill(MuIso[muind1]);
				for(size_t k = 0; k < gNMuPt2bins; ++k){
					if(MuPt[muind1] < gMuPt2bins[k]) continue;
					if(MuPt[muind1] > gMuPt2bins[k+1]) continue;
					hiso_ttbar_pt[1][k]->Fill(MuIso[muind1]);
				}
				for(size_t k = 0; k < gNNVrtxBins; ++k){
					if(NVrtx < gNVrtxBins[k]) continue;
					if(NVrtx > gNVrtxBins[k+1]) continue;
					hiso_ttbar_nv[1][k]->Fill(MuIso[muind1]);
				}
			}
			// ////////////////////////////////////////////////////
			// // SIGNAL SELECTION
			// if(isSSLLMuEvent(muind1, muind2)){
			// 	int fakemu = muind1;
			// 	if(isPromptMuon(muind1) &&  isPromptMuon(muind2)) continue;
			// 	if(isPromptMuon(muind1) && !isPromptMuon(muind2)) fakemu = muind2;
			// 	
			// 	hiso_ttbar[1]->Fill(MuIso[fakemu]);
			// 	for(size_t k = 0; k < gNMuPt2bins; ++k){
			// 		if(MuPt[fakemu] < gMuPt2bins[k]) continue;
			// 		if(MuPt[fakemu] > gMuPt2bins[k+1]) continue;
			// 		hiso_ttbar_pt[1][k]->Fill(MuIso[fakemu]);
			// 	}
			// 	for(size_t k = 0; k < gNNVrtxBins; ++k){
			// 		if(NVrtx < gNVrtxBins[k]) continue;
			// 		if(NVrtx > gNVrtxBins[k+1]) continue;
			// 		hiso_ttbar_nv[1][k]->Fill(MuIso[fakemu]);
			// 	}
			// }
			////////////////////////////////////////////////////
		}
		fSamples[TTJets]->cleanUp();
		cout << endl;
	}
	////////////////////////////////////////////////////
	
	// Create plots
	vector<int> mcsamples = fMCBGMuEnr;
	// vector<int> datasamples = fJMData;
	vector<int> datasamples = fMuData;

	for(size_t i = 0; i < gNSels; ++i){
		fOutputSubDir = "Isolation/Muons/";
		hiso_data[i]->SetXTitle(convertVarName("MuIso[0]"));
		hiso_data[i]->SetLineWidth(3);
		hiso_data[i]->SetLineColor(kBlack);
		hiso_data[i]->SetMarkerStyle(8);
		hiso_data[i]->SetMarkerColor(kBlack);
		hiso_data[i]->SetMarkerSize(1.2);
		
		hiso_ttbar[i]->SetXTitle(convertVarName("MuIso[0]"));
		hiso_ttbar[i]->SetLineWidth(3);
		hiso_ttbar[i]->SetLineColor(kRed);
		hiso_ttbar[i]->SetMarkerStyle(23);
		hiso_ttbar[i]->SetMarkerColor(kRed);
		hiso_ttbar[i]->SetMarkerSize(1.3);
		
		for(int k = 0; k < gNMuPt2bins; ++k){
			hiso_data_pt[i][k]->SetXTitle(convertVarName("MuIso[0]"));
			hiso_data_pt[i][k]->SetLineWidth(3);
			hiso_data_pt[i][k]->SetLineColor(kBlack);
			hiso_data_pt[i][k]->SetMarkerStyle(8);
			hiso_data_pt[i][k]->SetMarkerColor(kBlack);
			hiso_data_pt[i][k]->SetMarkerSize(1.2);

			hiso_ttbar_pt[i][k]->SetXTitle(convertVarName("MuIso[0]"));
			hiso_ttbar_pt[i][k]->SetLineWidth(3);
			hiso_ttbar_pt[i][k]->SetLineColor(kRed);
			hiso_ttbar_pt[i][k]->SetMarkerStyle(23);
			hiso_ttbar_pt[i][k]->SetMarkerColor(kRed);
			hiso_ttbar_pt[i][k]->SetMarkerSize(1.3);
		}
		for(int k = 0; k < gNNVrtxBins; ++k){
			hiso_data_nv[i][k]->SetXTitle(convertVarName("MuIso[0]"));
			hiso_data_nv[i][k]->SetLineWidth(3);
			hiso_data_nv[i][k]->SetLineColor(kBlack);
			hiso_data_nv[i][k]->SetMarkerStyle(8);
			hiso_data_nv[i][k]->SetMarkerColor(kBlack);
			hiso_data_nv[i][k]->SetMarkerSize(1.2);

			hiso_ttbar_nv[i][k]->SetXTitle(convertVarName("MuIso[0]"));
			hiso_ttbar_nv[i][k]->SetLineWidth(3);
			hiso_ttbar_nv[i][k]->SetLineColor(kRed);
			hiso_ttbar_nv[i][k]->SetMarkerStyle(23);
			hiso_ttbar_nv[i][k]->SetMarkerColor(kRed);
			hiso_ttbar_nv[i][k]->SetMarkerSize(1.3);
		}

		// Apply weights to MC histos
		for(size_t j = 0; j < gNSAMPLES; ++j){
			Sample *S = fSamples[j];
			float lumiscale = fLumiNorm / S->lumi;
			if(S->datamc == 0) continue;
			S->isoplots[0].hiso[i]->Scale(lumiscale);
			for(size_t k = 0; k < gNMuPt2bins; ++k){
				S->isoplots[0].hiso_pt[i][k]->Scale(lumiscale);
			}
			for(size_t k = 0; k < gNNVrtxBins; ++k){
				S->isoplots[0].hiso_nv[i][k]->Scale(lumiscale);
			}
		}

		// Fill data histo
		for(size_t j = 0; j < datasamples.size(); ++j){
			Sample *S = fSamples[datasamples[j]];
			hiso_data[i]->Add(S->isoplots[0].hiso[i]);
			hiso_data[i]->SetXTitle(convertVarName("MuIso[0]"));
			for(int k = 0; k < gNMuPt2bins; ++k){
				hiso_data_pt[i][k]->Add(S->isoplots[0].hiso_pt[i][k]);
				hiso_data_pt[i][k]->SetXTitle(convertVarName("MuIso[0]"));
			}
			for(int k = 0; k < gNNVrtxBins; ++k){
				hiso_data_nv[i][k]->Add(S->isoplots[0].hiso_nv[i][k]);
				hiso_data_nv[i][k]->SetXTitle(convertVarName("MuIso[0]"));
			}
		}

		// Scale to get equal integrals
		float intscale(0.);
		float intscale_pt[gNMuPt2bins];
		float intscale_nv[gNNVrtxBins];
		for(size_t j = 0; j < mcsamples.size(); ++j){
			Sample *S = fSamples[mcsamples[j]];
			intscale += S->isoplots[0].hiso[i]->Integral();
			for(int k = 0; k < gNMuPt2bins; ++k){
				intscale_pt[k] += S->isoplots[0].hiso_pt[i][k]->Integral();
			}
			for(int k = 0; k < gNNVrtxBins; ++k){
				intscale_nv[k] += S->isoplots[0].hiso_nv[i][k]->Integral();
			}
		}
		intscale = hiso_data[i]->Integral() / intscale;
		for(size_t j = 0; j < gNMuPt2bins; ++j) intscale_pt[j] = hiso_data_pt[i][j]->Integral() / intscale_pt[j];
		for(size_t j = 0; j < gNNVrtxBins; ++j) intscale_nv[j] = hiso_data_nv[i][j]->Integral() / intscale_nv[j];
		
		for(size_t j = 0; j < mcsamples.size(); ++j){
			Sample *S = fSamples[mcsamples[j]];
			S->isoplots[0].hiso[i]->Scale(intscale);
			for(int k = 0; k < gNMuPt2bins; ++k) S->isoplots[0].hiso_pt[i][k]->Scale(intscale_pt[k]);
			for(int k = 0; k < gNNVrtxBins; ++k) S->isoplots[0].hiso_nv[i][k]->Scale(intscale_nv[k]);
		}
		if(dottbar){
			hiso_ttbar[i]->Scale(hiso_data[i]->Integral() / hiso_ttbar[i]->Integral());
			for(int k = 0; k < gNMuPt2bins; ++k) hiso_ttbar_pt[i][k]->Scale(hiso_data_pt[i][k]->Integral() / hiso_ttbar_pt[i][k]->Integral());
			for(int k = 0; k < gNNVrtxBins; ++k) hiso_ttbar_nv[i][k]->Scale(hiso_data_nv[i][k]->Integral() / hiso_ttbar_nv[i][k]->Integral());			
		}
		

		// Fill MC stacks
		for(size_t j = 0; j < mcsamples.size();   ++j){
			Sample *S = fSamples[mcsamples[j]];
			TString s_name = S->sname;
			hiso_mc  [i]->Add(S->isoplots[0].hiso[i]);
			// sampleType Function: QCD = 1 , Top = 2, EWK = 3 , Rare = 4 , DB = 5
			if ( sampleType(s_name) == 1) hiso_qcd[i] ->Add( S->isoplots[0].hiso[i] );
			if ( sampleType(s_name) == 2) hiso_ttj[i] ->Add( S->isoplots[0].hiso[i] );
			if ( sampleType(s_name) == 3) hiso_ewk[i] ->Add( S->isoplots[0].hiso[i] );
			if ( sampleType(s_name) == 4) hiso_rare[i]->Add( S->isoplots[0].hiso[i] );
			if ( sampleType(s_name) == 5) hiso_db[i]  ->Add( S->isoplots[0].hiso[i] );

			for(int k = 0; k < gNMuPt2bins; ++k){
				if ( sampleType(s_name) == 1) hiso_qcd_pt [i][k]->Add( S->isoplots[0].hiso_pt[i][k] );
				if ( sampleType(s_name) == 2) hiso_ttj_pt [i][k]->Add( S->isoplots[0].hiso_pt[i][k] );
				if ( sampleType(s_name) == 3) hiso_ewk_pt [i][k]->Add( S->isoplots[0].hiso_pt[i][k] );
				if ( sampleType(s_name) == 4) hiso_rare_pt[i][k]->Add( S->isoplots[0].hiso_pt[i][k] );
				if ( sampleType(s_name) == 5) hiso_db_pt  [i][k]->Add( S->isoplots[0].hiso_pt[i][k] );
			}
			for(int k = 0; k < gNNVrtxBins; ++k){
				if ( sampleType(s_name) == 1) hiso_qcd_nv [i][k]->Add( S->isoplots[0].hiso_nv[i][k] );
				if ( sampleType(s_name) == 2) hiso_ttj_nv [i][k]->Add( S->isoplots[0].hiso_nv[i][k] );
				if ( sampleType(s_name) == 3) hiso_ewk_nv [i][k]->Add( S->isoplots[0].hiso_nv[i][k] );
				if ( sampleType(s_name) == 4) hiso_rare_nv[i][k]->Add( S->isoplots[0].hiso_nv[i][k] );
				if ( sampleType(s_name) == 5) hiso_db_nv  [i][k]->Add( S->isoplots[0].hiso_nv[i][k] );
			}

		}
		hiso_mc_s[i]->Add(hiso_qcd[i]);
		hiso_mc_s[i]->Add(hiso_db[i]);
		hiso_mc_s[i]->Add(hiso_ewk[i]);
		hiso_mc_s[i]->Add(hiso_rare[i]);
		hiso_mc_s[i]->Add(hiso_ttj[i]);
		hiso_mc_s[i]->Draw("goff");
		hiso_mc_s[i]->GetXaxis()->SetTitle(convertVarName("MuIso[0]"));

		for(int k = 0; k < gNMuPt2bins; ++k){
			hiso_mc_pt_s[i][k]->Add(hiso_qcd_pt[i][k]);
			hiso_mc_pt_s[i][k]->Add(hiso_db_pt[i][k]);
			hiso_mc_pt_s[i][k]->Add(hiso_ewk_pt[i][k]);
			hiso_mc_pt_s[i][k]->Add(hiso_rare_pt[i][k]);
			hiso_mc_pt_s[i][k]->Add(hiso_ttj_pt[i][k]);
			hiso_mc_pt_s[i][k]->Draw("goff");
			hiso_mc_pt_s[i][k]->GetXaxis()->SetTitle(convertVarName("MuIso[0]"));
		}

		for(int k = 0; k < gNNVrtxBins; ++k){
			hiso_mc_nv_s[i][k]->Add(hiso_qcd_nv[i][k]);
			hiso_mc_nv_s[i][k]->Add(hiso_db_nv[i][k]);
			hiso_mc_nv_s[i][k]->Add(hiso_ewk_nv[i][k]);
			hiso_mc_nv_s[i][k]->Add(hiso_rare_nv[i][k]);
			hiso_mc_nv_s[i][k]->Add(hiso_ttj_nv[i][k]);
			hiso_mc_nv_s[i][k]->Draw("goff");
			hiso_mc_nv_s[i][k]->GetXaxis()->SetTitle(convertVarName("MuIso[0]"));
		}
		// Fill MC stacks
		//for(size_t j = 0; j < mcsamples.size();   ++j){
		//	Sample *S = fSamples[mcsamples[j]];
		//	hiso_mc  [i]->Add(S->isoplots[0].hiso[i]);
		//	hiso_mc_s[i]->Add(S->isoplots[0].hiso[i]);
		//	hiso_mc_s[i]->Draw("goff");
		//	hiso_mc_s[i]->GetXaxis()->SetTitle(convertVarName("MuIso[0]"));
		//	for(int k = 0; k < gNMuPt2bins; ++k){
		//		hiso_mc_pt  [i][k]->Add(S->isoplots[0].hiso_pt[i][k]);
		//		hiso_mc_pt_s[i][k]->Add(S->isoplots[0].hiso_pt[i][k]);
		//		hiso_mc_pt_s[i][k]->Draw("goff");
		//		hiso_mc_pt_s[i][k]->GetXaxis()->SetTitle(convertVarName("MuIso[0]"));
		//	}
		//	for(int k = 0; k < gNNVrtxBins; ++k){
		//		hiso_mc_nv  [i][k]->Add(S->isoplots[0].hiso_nv[i][k]);
		//		hiso_mc_nv_s[i][k]->Add(S->isoplots[0].hiso_nv[i][k]);
		//		hiso_mc_nv_s[i][k]->Draw("goff");
		//		hiso_mc_nv_s[i][k]->GetXaxis()->SetTitle(convertVarName("MuIso[0]"));
		//	}
		//}

		double max1 = hiso_mc_s[i]->GetMaximum();
		double max2 = hiso_data[i]->GetMaximum();
		double max = max1>max2?max1:max2;
		hiso_mc_s[i]->SetMaximum(1.5*max);
		hiso_data[i]->SetMaximum(1.5*max);

		int bin0   = hiso_data[i]->FindBin(0.0);
		int bin015 = hiso_data[i]->FindBin(0.15) - 1; // bins start at lower edge...
		int bin1   = hiso_data[i]->FindBin(1.0)  - 1;
		printf("bin 0: %3.0i bin015: %3.0i bin1: %3.0i \n", bin0, bin015, bin1);
		float ratio_data  = hiso_data[i] ->Integral(bin0, bin015) / hiso_data[i] ->Integral(bin0, bin1);
		float ratio_mc    = hiso_mc[i]   ->Integral(bin0, bin015) / hiso_mc[i]   ->Integral(bin0, bin1);
		float ratio_ttbar(0.);
		if(dottbar) ratio_ttbar = hiso_ttbar[i]->Integral(bin0, bin015) / hiso_ttbar[i]->Integral(bin0, bin1);

		TCanvas *c_temp = new TCanvas("MuIso" + IsoPlots::sel_name[i], "Muon Isolation in Data vs MC", 0, 0, 800, 600);
		c_temp->cd();

		TLegend *leg = new TLegend(0.15,0.65,0.40,0.88);
		// TLegend *leg = new TLegend(0.75,0.60,0.89,0.88);
		leg->AddEntry(hiso_data[i], "Data","p");
		if(dottbar) leg->AddEntry(hiso_ttbar[i], "TTbar fake","p");
		leg->AddEntry(hiso_ttj[i],  "Top","f");
		leg->AddEntry(hiso_rare[i], "Rare SM","f");
		leg->AddEntry(hiso_ewk[i],  "Single Boson","f");
		leg->AddEntry(hiso_db[i],   "Di-Boson","f");
		leg->AddEntry(hiso_qcd[i],  "QCD","f");
		leg->SetFillStyle(0);
		leg->SetTextFont(42);
		leg->SetBorderSize(0);

		// gPad->SetLogy();
		hiso_mc_s[i]->Draw("hist");
		if(dottbar) hiso_ttbar[i]->DrawCopy("PE X0 same");
		hiso_data[i]->DrawCopy("PE X0 same");
		leg->Draw();
		lat->DrawLatex(0.70,0.92, Form("L_{int.} = %2.1f fb^{-1}", fLumiNorm/1000.));
		lat->DrawLatex(0.75,0.85, Form("R^{T/L}_{Data}  = %4.2f", ratio_data));
		lat->DrawLatex(0.75,0.80, Form("R^{T/L}_{MC}   = %4.2f", ratio_mc));
		lat->SetTextColor(kRed);
		if(dottbar) lat->DrawLatex(0.75,0.75, Form("R^{T/L}_{TTbar} = %4.2f", ratio_ttbar));
		lat->SetTextColor(kBlack);
		

		// Util::PrintNoEPS(c_temp, "MuIso" + IsoPlots::sel_name[i], fOutputDir + fOutputSubDir, NULL);
		Util::PrintPDF(c_temp, "MuIso" + IsoPlots::sel_name[i], fOutputDir + fOutputSubDir);

		for(int k = 0; k < gNMuPt2bins; ++k){
			fOutputSubDir = "Isolation/Muons/PtBinned/";
			ratio_data  = hiso_data_pt[i][k] ->Integral(bin0, bin015) / hiso_data_pt[i][k] ->Integral(bin0, bin1);
			ratio_mc    = hiso_mc_pt[i][k]   ->Integral(bin0, bin015) / hiso_mc_pt[i][k]   ->Integral(bin0, bin1);
			ratio_ttbar = 0.;
			if(dottbar) hiso_ttbar_pt[i][k]->Integral(bin0, bin015) / hiso_ttbar_pt[i][k]->Integral(bin0, bin1);

			TCanvas *c_temp = new TCanvas(Form("MuIso%s_pt_%d", IsoPlots::sel_name[i].Data(), k), "Muon Isolation in Data vs MC", 0, 0, 800, 600);
			c_temp->cd();

			max1 = hiso_mc_pt_s[i][k]->GetMaximum();
			max2 = hiso_data_pt[i][k]->GetMaximum();
			max = max1>max2?max1:max2;
			hiso_mc_pt_s[i][k]->SetMaximum(1.5*max);
			hiso_data_pt[i][k]->SetMaximum(1.5*max);

			TLegend *leg_pt = new TLegend(0.15,0.65,0.40,0.88);
			// TLegend *leg = new TLegend(0.75,0.60,0.89,0.88);
			leg_pt->AddEntry(hiso_data_pt[i][k], "Data","p");
			if(dottbar) leg_pt->AddEntry(hiso_ttbar_pt[i][k], "TTbar fake","p");
			leg_pt->AddEntry(hiso_ttj_pt  [i][k], "Top","f");
			leg_pt->AddEntry(hiso_rare_pt [i][k], "Rare SM","f");
			leg_pt->AddEntry(hiso_ewk_pt  [i][k], "Single Boson","f");
			leg_pt->AddEntry(hiso_db_pt   [i][k], "Di-Boson","f");
			leg_pt->AddEntry(hiso_qcd_pt  [i][k], "QCD","f");
			leg_pt->SetFillStyle(0);
			leg_pt->SetTextFont(42);
			leg_pt->SetBorderSize(0);

			// gPad->SetLogy();
			hiso_mc_pt_s[i][k]->Draw("hist");
			if(dottbar) hiso_ttbar_pt[i][k]->DrawCopy("PE X0 same");
			hiso_data_pt[i][k]->DrawCopy("PE X0 same");
			leg_pt->Draw();
			lat->DrawLatex(0.20,0.92, Form("p_{T}(#mu) %3.0f - %3.0f GeV", getPt2Bins(Muon)[k], getPt2Bins(Muon)[k+1]));
			lat->DrawLatex(0.70,0.92, Form("L_{int.} = %2.1f fb^{-1}", fLumiNorm/1000.));
			lat->DrawLatex(0.75,0.85, Form("R^{T/L}_{Data}  = %4.2f", ratio_data));
			lat->DrawLatex(0.75,0.80, Form("R^{T/L}_{MC}   = %4.2f", ratio_mc));
			lat->SetTextColor(kRed);
			if(dottbar) lat->DrawLatex(0.75,0.75, Form("R^{T/L}_{TTbar} = %4.2f", ratio_ttbar));
			lat->SetTextColor(kBlack);

			// Util::PrintNoEPS(c_temp, Form("MuIso%s_pt_%d", IsoPlots::sel_name[i].Data(), k), fOutputDir + fOutputSubDir, NULL);
			Util::PrintPDF(c_temp, Form("MuIso%s_pt_%d", IsoPlots::sel_name[i].Data(), k), fOutputDir + fOutputSubDir);
		}
		for(int k = 0; k < gNNVrtxBins; ++k){
			fOutputSubDir = "Isolation/Muons/NVrtxBinned/";
			ratio_data  = hiso_data_nv[i][k] ->Integral(bin0, bin015) / hiso_data_nv[i][k] ->Integral(bin0, bin1);
			ratio_mc    = hiso_mc_nv[i][k]   ->Integral(bin0, bin015) / hiso_mc_nv[i][k]   ->Integral(bin0, bin1);
			ratio_ttbar = 0.;
			if(dottbar) hiso_ttbar_nv[i][k]->Integral(bin0, bin015) / hiso_ttbar_nv[i][k]->Integral(bin0, bin1);

			TCanvas *c_temp = new TCanvas(Form("MuIso%s_nv_%d", IsoPlots::sel_name[i].Data(), k), "Muon Isolation in Data vs MC", 0, 0, 800, 600);
			c_temp->cd();

			max1 = hiso_mc_nv_s[i][k]->GetMaximum();
			max2 = hiso_data_nv[i][k]->GetMaximum();
			max = max1>max2?max1:max2;
			hiso_mc_nv_s[i][k]->SetMaximum(1.5*max);
			hiso_data_nv[i][k]->SetMaximum(1.5*max);

			TLegend *leg_nv = new TLegend(0.15,0.65,0.40,0.88);
			// TLegend *leg = new TLegend(0.75,0.60,0.89,0.88);
			leg_nv->AddEntry(hiso_data_nv[i][k], "Data","p");
			if(dottbar) leg_nv->AddEntry(hiso_ttbar_nv[i][k], "TTbar fake","p");
			leg_nv->AddEntry(hiso_ttj_nv  [i][k],  "Top","f");
			leg_nv->AddEntry(hiso_rare_nv [i][k],  "Rare SM","f");
			leg_nv->AddEntry(hiso_ewk_nv  [i][k],  "Single Boson","f");
			leg_nv->AddEntry(hiso_db_nv   [i][k],  "Di-Boson","f");
			leg_nv->AddEntry(hiso_qcd_nv  [i][k],  "QCD","f");
			leg_nv->SetFillStyle(0);
			leg_nv->SetTextFont(42);
			leg_nv->SetBorderSize(0);

			// gPad->SetLogy();
			hiso_mc_nv_s[i][k]->Draw("hist");
			if(dottbar) hiso_ttbar_nv[i][k]->DrawCopy("PE X0 same");
			hiso_data_nv[i][k]->DrawCopy("PE X0 same");
			leg_nv->Draw();
			lat->DrawLatex(0.20,0.92, Form("N_{Vrtx.} %2.0f - %2.0f", gNVrtxBins[k], gNVrtxBins[k+1]));
			lat->DrawLatex(0.70,0.92, Form("L_{int.} = %2.1f fb^{-1}", fLumiNorm/1000.));
			lat->DrawLatex(0.75,0.85, Form("R^{T/L}_{Data}  = %4.2f", ratio_data));
			lat->DrawLatex(0.75,0.80, Form("R^{T/L}_{MC}   = %4.2f", ratio_mc));
			lat->SetTextColor(kRed);
			if(dottbar) lat->DrawLatex(0.75,0.75, Form("R^{T/L}_{TTbar} = %4.2f", ratio_ttbar));
			lat->SetTextColor(kBlack);

			// Util::PrintNoEPS(c_temp, Form("MuIso%s_nv_%d", IsoPlots::sel_name[i].Data(), k), fOutputDir + fOutputSubDir, NULL);
			Util::PrintPDF(c_temp, Form("MuIso%s_nv_%d", IsoPlots::sel_name[i].Data(), k), fOutputDir + fOutputSubDir);
		}
	}
}
void SSDLPlotter::makeElIsolationPlots(bool dottbar){
	char cmd[100];
    sprintf(cmd,"mkdir -p %s%s", fOutputDir.Data(), fOutputSubDir.Data());
    system(cmd);

	TH1D    *hiso_data [gNSels];
	TH1D    *hiso_mc   [gNSels];
	TH1D	*hiso_ttbar[gNSels];
	TH1D    *hiso_qcd  [gNSels];
	TH1D    *hiso_ttj  [gNSels];
	TH1D    *hiso_ewk  [gNSels];
	TH1D    *hiso_rare [gNSels];
	TH1D    *hiso_db   [gNSels];
	THStack *hiso_mc_s [gNSels];

	TH1D    *hiso_data_pt [gNSels][gNElPt2bins];
	TH1D    *hiso_mc_pt   [gNSels][gNElPt2bins];
	TH1D    *hiso_ttbar_pt[gNSels][gNElPt2bins];
	TH1D    *hiso_qcd_pt  [gNSels][gNElPt2bins];
	TH1D    *hiso_ttj_pt  [gNSels][gNElPt2bins];
	TH1D    *hiso_ewk_pt  [gNSels][gNElPt2bins];
	TH1D    *hiso_rare_pt [gNSels][gNElPt2bins];
	TH1D    *hiso_db_pt   [gNSels][gNElPt2bins];
	THStack *hiso_mc_pt_s [gNSels][gNElPt2bins];

	TH1D    *hiso_data_nv [gNSels][gNNVrtxBins];
	TH1D    *hiso_mc_nv   [gNSels][gNNVrtxBins];
	TH1D    *hiso_ttbar_nv[gNSels][gNNVrtxBins];
	TH1D    *hiso_qcd_nv  [gNSels][gNNVrtxBins];
	TH1D    *hiso_ttj_nv  [gNSels][gNNVrtxBins];
	TH1D    *hiso_ewk_nv  [gNSels][gNNVrtxBins];
	TH1D    *hiso_rare_nv [gNSels][gNNVrtxBins];
	TH1D    *hiso_db_nv   [gNSels][gNNVrtxBins];
	THStack *hiso_mc_nv_s [gNSels][gNNVrtxBins];

	TLatex *lat = new TLatex();
	lat->SetNDC(kTRUE);
	lat->SetTextColor(kBlack);
	lat->SetTextSize(0.04);

	// Create histograms
	for(size_t i = 0; i < gNSels; ++i){
		hiso_data [i] = new TH1D("ElIsoData_"          + IsoPlots::sel_name[i], "Electron Isolation in Data for "  + IsoPlots::sel_name[i], IsoPlots::nbins[i], 0., 0.6);
		hiso_mc   [i] = new TH1D("ElIsoMC_"            + IsoPlots::sel_name[i], "Electron Isolation in MC for "    + IsoPlots::sel_name[i], IsoPlots::nbins[i], 0., 0.6);
		hiso_ttbar[i] = new TH1D("ElIsoTTbar_"         + IsoPlots::sel_name[i], "Electron Isolation in TTbar for " + IsoPlots::sel_name[i], IsoPlots::nbins[i], 0., 0.6);
		hiso_qcd  [i] = new TH1D("ElIsoQCD_"           + IsoPlots::sel_name[i], "Electron Isolation in QCD   for " + IsoPlots::sel_name[i], IsoPlots::nbins[i], 0., 0.6);
		hiso_ttj  [i] = new TH1D("ElIsoTop_"           + IsoPlots::sel_name[i], "Electron Isolation in Top   for " + IsoPlots::sel_name[i], IsoPlots::nbins[i], 0., 0.6);
		hiso_ewk  [i] = new TH1D("ElIsoEWK_"           + IsoPlots::sel_name[i], "Electron Isolation in EWK   for " + IsoPlots::sel_name[i], IsoPlots::nbins[i], 0., 0.6);
		hiso_rare [i] = new TH1D("ElIsoRARE_"          + IsoPlots::sel_name[i], "Electron Isolation in RARE  for " + IsoPlots::sel_name[i], IsoPlots::nbins[i], 0., 0.6);
		hiso_db   [i] = new TH1D("ElIsoDB_"            + IsoPlots::sel_name[i], "Electron Isolation in DB    for " + IsoPlots::sel_name[i], IsoPlots::nbins[i], 0., 0.6);
		hiso_mc_s [i] = new THStack("ElIsoMC_stacked_" + IsoPlots::sel_name[i], "Electron Isolation in MC for "    + IsoPlots::sel_name[i]);
		hiso_data [i] ->Sumw2();
		hiso_mc   [i] ->Sumw2();
		hiso_ttbar[i] ->Sumw2();
		hiso_qcd  [i] ->Sumw2();
		hiso_ttj  [i] ->Sumw2();
		hiso_ewk  [i] ->Sumw2();
		hiso_rare [i] ->Sumw2();
		hiso_db   [i] ->Sumw2();

		for(int k = 0; k < gNElPt2bins; ++k){
			hiso_data_pt[i][k]  = new TH1D(Form("ElIsoData_%s_Pt%d",          IsoPlots::sel_name[i].Data(), k), "Electron Isolation in Data for "  + IsoPlots::sel_name[i], IsoPlots::nbins[i], 0., 0.6);
			hiso_mc_pt  [i][k]  = new TH1D(Form("ElIsoMC_%s_Pt%d",            IsoPlots::sel_name[i].Data(), k), "Electron Isolation in MC for "    + IsoPlots::sel_name[i], IsoPlots::nbins[i], 0., 0.6);
			hiso_ttbar_pt[i][k] = new TH1D(Form("ElIsoTTbar_%s_Pt%d",         IsoPlots::sel_name[i].Data(), k), "Electron Isolation in TTbar for " + IsoPlots::sel_name[i], IsoPlots::nbins[i], 0., 0.6);
			hiso_qcd_pt [i][k]  = new TH1D(Form("ElIsoQCD_%s_Pt%d",           IsoPlots::sel_name[i].Data(), k), "Electron Isolation in QCD   for " + IsoPlots::sel_name[i], IsoPlots::nbins[i], 0., 0.6);
			hiso_ttj_pt [i][k]  = new TH1D(Form("ElIsoTop_%s_Pt%d",           IsoPlots::sel_name[i].Data(), k), "Electron Isolation in Top   for " + IsoPlots::sel_name[i], IsoPlots::nbins[i], 0., 0.6);
			hiso_ewk_pt [i][k]  = new TH1D(Form("ElIsoEWK_%s_Pt%d",           IsoPlots::sel_name[i].Data(), k), "Electron Isolation in EWK   for " + IsoPlots::sel_name[i], IsoPlots::nbins[i], 0., 0.6);
			hiso_rare_pt[i][k]  = new TH1D(Form("ElIsoRARE_%s_Pt%d",          IsoPlots::sel_name[i].Data(), k), "Electron Isolation in RARE  for " + IsoPlots::sel_name[i], IsoPlots::nbins[i], 0., 0.6);
			hiso_db_pt  [i][k]  = new TH1D(Form("ElIsoDB_%s_Pt%d",            IsoPlots::sel_name[i].Data(), k), "Electron Isolation in DB    for " + IsoPlots::sel_name[i], IsoPlots::nbins[i], 0., 0.6);
			hiso_mc_pt_s[i][k]  = new THStack(Form("ElIsoMC_stacked_%s_Pt%d", IsoPlots::sel_name[i].Data(), k), "Electron Isolation in MC for "    + IsoPlots::sel_name[i]);
			hiso_data_pt [i][k]->Sumw2();
			hiso_mc_pt   [i][k]->Sumw2();
			hiso_ttbar_pt[i][k]->Sumw2();
			hiso_qcd_pt  [i][k]->Sumw2();
			hiso_ttj_pt  [i][k]->Sumw2();
			hiso_ewk_pt  [i][k]->Sumw2();
			hiso_rare_pt [i][k]->Sumw2();
			hiso_db_pt   [i][k]->Sumw2();
		}
		for(int k = 0; k < gNNVrtxBins; ++k){
			hiso_data_nv[i][k]  = new TH1D(Form("ElIsoData_%s_NVtx%d",      IsoPlots::sel_name[i].Data(), k), "Electron Isolation in Data for "  + IsoPlots::sel_name[i], IsoPlots::nbins[i], 0., 0.6);
			hiso_mc_nv  [i][k]  = new TH1D(Form("ElIsoMC_%s_NVtx%d",        IsoPlots::sel_name[i].Data(), k), "Electron Isolation in MC for "    + IsoPlots::sel_name[i], IsoPlots::nbins[i], 0., 0.6);
			hiso_ttbar_nv[i][k] = new TH1D(Form("ElIsoTTbar_%s_NVtx%d",     IsoPlots::sel_name[i].Data(), k), "Electron Isolation in TTbar for " + IsoPlots::sel_name[i], IsoPlots::nbins[i], 0., 0.6);
			hiso_qcd_nv [i][k]  = new TH1D(Form("ElIsoQCD_%s_NVtx%d",       IsoPlots::sel_name[i].Data(), k), "Electron Isolation in QCD   for " + IsoPlots::sel_name[i], IsoPlots::nbins[i], 0., 0.6);
			hiso_ttj_nv [i][k]  = new TH1D(Form("ElIsoTop_%s_NVtx%d",       IsoPlots::sel_name[i].Data(), k), "Electron Isolation in Top   for " + IsoPlots::sel_name[i], IsoPlots::nbins[i], 0., 0.6);
			hiso_ewk_nv [i][k]  = new TH1D(Form("ElIsoEWK_%s_NVtx%d",       IsoPlots::sel_name[i].Data(), k), "Electron Isolation in EWK   for " + IsoPlots::sel_name[i], IsoPlots::nbins[i], 0., 0.6);
			hiso_rare_nv[i][k]  = new TH1D(Form("ElIsoRARE_%s_NVtx%d",      IsoPlots::sel_name[i].Data(), k), "Electron Isolation in RARE  for " + IsoPlots::sel_name[i], IsoPlots::nbins[i], 0., 0.6);
			hiso_db_nv  [i][k]  = new TH1D(Form("ElIsoDB_%s_NVtx%d",        IsoPlots::sel_name[i].Data(), k), "Electron Isolation in DB    for " + IsoPlots::sel_name[i], IsoPlots::nbins[i], 0., 0.6);
			hiso_mc_nv_s[i][k]  = new THStack(Form("ElIsoMC_stacked_%s_NVtx%d", IsoPlots::sel_name[i].Data(), k), "Electron Isolation in MC for "    + IsoPlots::sel_name[i]);
			hiso_data_nv [i][k]->Sumw2();
			hiso_mc_nv   [i][k]->Sumw2();
			hiso_ttbar_nv[i][k]->Sumw2();
			hiso_qcd_nv  [i][k]->Sumw2();
			hiso_ttj_nv  [i][k]->Sumw2();
			hiso_ewk_nv  [i][k]->Sumw2();
			hiso_rare_nv [i][k]->Sumw2();
			hiso_db_nv   [i][k]->Sumw2();
		}
	}

	for(size_t i = 0; i < gNSels; ++i){
		hiso_qcd [i]->SetFillColor(kYellow-4);
		hiso_db  [i]->SetFillColor(kSpring-9);
		hiso_ewk [i]->SetFillColor(kGreen +1);
		hiso_ttj [i]->SetFillColor(kAzure-5);
		hiso_rare[i]->SetFillColor(kAzure+8);

		for(int k = 0; k < gNElPt2bins; ++k){
			hiso_qcd_pt [i][k]->SetFillColor(kYellow-4);
			hiso_db_pt  [i][k]->SetFillColor(kSpring-9);
			hiso_ewk_pt [i][k]->SetFillColor(kGreen +1);
			hiso_ttj_pt [i][k]->SetFillColor(kAzure-5);
			hiso_rare_pt[i][k]->SetFillColor(kAzure+8);
		}

		for(int k = 0; k < gNNVrtxBins; ++k){
			hiso_qcd_nv [i][k]->SetFillColor(kYellow-4);
			hiso_db_nv  [i][k]->SetFillColor(kSpring-9);
			hiso_ewk_nv [i][k]->SetFillColor(kGreen +1);
			hiso_ttj_nv [i][k]->SetFillColor(kAzure-5);
			hiso_rare_nv[i][k]->SetFillColor(kAzure+8);
		}
	}

	////////////////////////////////////////////////////
	if(dottbar){
		// Fill ttbar histos
		// Sample loop
		TTree *tree = fSamples[TTJets]->getTree();

		// Event loop
		tree->ResetBranchAddresses();
		InitMC(tree);

		if (fChain == 0) return;
		Long64_t nentries = fChain->GetEntriesFast();
		Long64_t nbytes = 0, nb = 0;
		for (Long64_t jentry=0; jentry<nentries;jentry++) {
			printProgress(jentry, nentries, fSamples[TTJets]->name);

			Long64_t ientry = LoadTree(jentry);
			if (ientry < 0) break;
			nb = fChain->GetEntry(jentry);   nbytes += nb;

			int elind1(-1), elind2(-1);
			if(hasLooseElectrons(elind1, elind2) < 1) continue;

			// Common event selections
			if(!passesJet50Cut()) continue; // make trigger 100% efficient

			// Common object selections
			if(!isLooseElectron(elind1)) continue;
			// if(ElIsGoodElId_WP80[elind1] != 1) return false; // apply tight ID for the iso plots?

			// Select genmatched fake muons
			if(ElGenMType[elind1] == 2 || ElGenMType[elind1] == 4) continue;
			// Exclude also conversions here?

			if(ElPt[elind1] < fC_minEl2pt) continue;
			if(ElPt[elind1] > gElPt2bins[gNElPt2bins]) continue;


			////////////////////////////////////////////////////
			// MOST LOOSE SELECTION
			hiso_ttbar[0]->Fill(ElRelIso[elind1]);
			for(size_t k = 0; k < gNElPt2bins; ++k){
				if(ElPt[elind1] < gElPt2bins[k]) continue;
				if(ElPt[elind1] > gElPt2bins[k+1]) continue;
				hiso_ttbar_pt[0][k]->Fill(ElRelIso[elind1]);
			}
			for(size_t k = 0; k < gNNVrtxBins; ++k){
				if(NVrtx < gNVrtxBins[k]) continue;
				if(NVrtx > gNVrtxBins[k+1]) continue;
				hiso_ttbar_nv[0][k]->Fill(ElRelIso[elind1]);
			}

			////////////////////////////////////////////////////
			// SIGNAL SUPPRESSED SELECTION
			if(isSigSupElEvent()){
				hiso_ttbar[1]->Fill(ElRelIso[elind1]);
				for(size_t k = 0; k < gNElPt2bins; ++k){
					if(ElPt[elind1] < gElPt2bins[k]) continue;
					if(ElPt[elind1] > gElPt2bins[k+1]) continue;
					hiso_ttbar_pt[1][k]->Fill(ElRelIso[elind1]);
				}
				for(size_t k = 0; k < gNNVrtxBins; ++k){
					if(NVrtx < gNVrtxBins[k]) continue;
					if(NVrtx > gNVrtxBins[k+1]) continue;
					hiso_ttbar_nv[1][k]->Fill(ElRelIso[elind1]);
				}
			}
			////////////////////////////////////////////////////
		}
		fSamples[TTJets]->cleanUp();
		cout << endl;
	}
	////////////////////////////////////////////////////

	// Create plots
	vector<int> mcsamples = fMCBG;
	vector<int> datasamples = fEGData;

	for(size_t i = 0; i < gNSels; ++i){
		fOutputSubDir = "Isolation/Electrons/";
		hiso_data[i]->SetXTitle(convertVarName("ElRelIso[0]"));
		hiso_data[i]->SetLineWidth(3);
		hiso_data[i]->SetLineColor(kBlack);
		hiso_data[i]->SetMarkerStyle(8);
		hiso_data[i]->SetMarkerColor(kBlack);
		hiso_data[i]->SetMarkerSize(1.2);

		hiso_ttbar[i]->SetXTitle(convertVarName("ElRelIso[0]"));
		hiso_ttbar[i]->SetLineWidth(3);
		hiso_ttbar[i]->SetLineColor(kRed);
		hiso_ttbar[i]->SetMarkerStyle(23);
		hiso_ttbar[i]->SetMarkerColor(kRed);
		hiso_ttbar[i]->SetMarkerSize(1.3);
		
		for(int k = 0; k < gNElPt2bins; ++k){
			hiso_data_pt[i][k]->SetXTitle(convertVarName("ElRelIso[0]"));
			hiso_data_pt[i][k]->SetLineWidth(3);
			hiso_data_pt[i][k]->SetLineColor(kBlack);
			hiso_data_pt[i][k]->SetMarkerStyle(8);
			hiso_data_pt[i][k]->SetMarkerColor(kBlack);
			hiso_data_pt[i][k]->SetMarkerSize(1.2);

			hiso_ttbar_pt[i][k]->SetXTitle(convertVarName("ElRelIso[0]"));
			hiso_ttbar_pt[i][k]->SetLineWidth(3);
			hiso_ttbar_pt[i][k]->SetLineColor(kRed);
			hiso_ttbar_pt[i][k]->SetMarkerStyle(23);
			hiso_ttbar_pt[i][k]->SetMarkerColor(kRed);
			hiso_ttbar_pt[i][k]->SetMarkerSize(1.3);
		}
		for(int k = 0; k < gNNVrtxBins; ++k){
			hiso_data_nv[i][k]->SetXTitle(convertVarName("ElRelIso[0]"));
			hiso_data_nv[i][k]->SetLineWidth(3);
			hiso_data_nv[i][k]->SetLineColor(kBlack);
			hiso_data_nv[i][k]->SetMarkerStyle(8);
			hiso_data_nv[i][k]->SetMarkerColor(kBlack);
			hiso_data_nv[i][k]->SetMarkerSize(1.2);

			hiso_ttbar_nv[i][k]->SetXTitle(convertVarName("ElRelIso[0]"));
			hiso_ttbar_nv[i][k]->SetLineWidth(3);
			hiso_ttbar_nv[i][k]->SetLineColor(kRed);
			hiso_ttbar_nv[i][k]->SetMarkerStyle(23);
			hiso_ttbar_nv[i][k]->SetMarkerColor(kRed);
			hiso_ttbar_nv[i][k]->SetMarkerSize(1.3);
		}

		// Apply weights to MC histos
		for(size_t j = 0; j < gNSAMPLES; ++j){
			Sample *S = fSamples[j];
			float lumiscale = fLumiNorm / S->lumi;
			if(S->datamc == 0) continue;
			S->isoplots[1].hiso[i]->Scale(lumiscale);
			for(size_t k = 0; k < gNElPt2bins; ++k){
				S->isoplots[1].hiso_pt[i][k]->Scale(lumiscale);
			}
			for(size_t k = 0; k < gNNVrtxBins; ++k){
				S->isoplots[1].hiso_nv[i][k]->Scale(lumiscale);
			}
		}

		// Fill data histo
		for(size_t j = 0; j < datasamples.size(); ++j){
			Sample *S = fSamples[datasamples[j]];
			hiso_data[i]->Add(S->isoplots[1].hiso[i]);
			hiso_data[i]->SetXTitle(convertVarName("ElRelIso[0]"));
			for(int k = 0; k < gNElPt2bins; ++k){
				hiso_data_pt[i][k]->Add(S->isoplots[1].hiso_pt[i][k]);
				hiso_data_pt[i][k]->SetXTitle(convertVarName("ElRelIso[0]"));
			}
			for(int k = 0; k < gNNVrtxBins; ++k){
				hiso_data_nv[i][k]->Add(S->isoplots[1].hiso_nv[i][k]);
				hiso_data_nv[i][k]->SetXTitle(convertVarName("ElRelIso[0]"));
			}
		}

		// Scale to get equal integrals
		float intscale(0.);
		float intscale_pt[gNElPt2bins];
		float intscale_nv[gNNVrtxBins];
		for(size_t j = 0; j < mcsamples.size();   ++j){
			Sample *S = fSamples[mcsamples[j]];
			intscale += S->isoplots[1].hiso[i]->Integral();
			for(int k = 0; k < gNElPt2bins; ++k){
				intscale_pt[k] += S->isoplots[1].hiso_pt[i][k]->Integral();
			}
			for(int k = 0; k < gNNVrtxBins; ++k){
				intscale_nv[k] += S->isoplots[1].hiso_nv[i][k]->Integral();
			}
		}
		intscale = hiso_data[i]->Integral() / intscale;
		for(size_t j = 0; j < gNElPt2bins; ++j) intscale_pt[j] = hiso_data_pt[i][j]->Integral() / intscale_pt[j];
		for(size_t j = 0; j < gNNVrtxBins; ++j) intscale_nv[j] = hiso_data_nv[i][j]->Integral() / intscale_nv[j];
		
		for(size_t j = 0; j < mcsamples.size();   ++j){
			Sample *S = fSamples[mcsamples[j]];			
			S->isoplots[1].hiso[i]->Scale(intscale);
			for(int k = 0; k < gNElPt2bins; ++k){
				S->isoplots[1].hiso_pt[i][k]->Scale(intscale_pt[k]);
			}
			for(int k = 0; k < gNNVrtxBins; ++k){
				S->isoplots[1].hiso_nv[i][k]->Scale(intscale_nv[k]);
			}
		}
		if(dottbar){
			hiso_ttbar[i]->Scale(hiso_data[i]->Integral() / hiso_ttbar[i]->Integral());
			for(int k = 0; k < gNElPt2bins; ++k) hiso_ttbar_pt[i][k]->Scale(hiso_data_pt[i][k]->Integral() / hiso_ttbar_pt[i][k]->Integral());
			for(int k = 0; k < gNNVrtxBins; ++k) hiso_ttbar_nv[i][k]->Scale(hiso_data_nv[i][k]->Integral() / hiso_ttbar_nv[i][k]->Integral());			
		}

		// Fill MC stacks
		for(size_t j = 0; j < mcsamples.size();   ++j){
			Sample *S = fSamples[mcsamples[j]];
			TString s_name = S->sname;
			hiso_mc  [i]->Add(S->isoplots[1].hiso[i]);
			// sampleType Function: QCD = 1 , Top = 2, EWK = 3 , Rare = 4 , DB = 5
			if ( sampleType(s_name) == 1) hiso_qcd[i] ->Add( S->isoplots[1].hiso[i] );
			if ( sampleType(s_name) == 2) hiso_ttj[i] ->Add( S->isoplots[1].hiso[i] );
			if ( sampleType(s_name) == 3) hiso_ewk[i] ->Add( S->isoplots[1].hiso[i] );
			if ( sampleType(s_name) == 4) hiso_rare[i]->Add( S->isoplots[1].hiso[i] );
			if ( sampleType(s_name) == 5) hiso_db[i]  ->Add( S->isoplots[1].hiso[i] );

			for(int k = 0; k < gNMuPt2bins; ++k){
				if ( sampleType(s_name) == 1) hiso_qcd_pt [i][k]->Add( S->isoplots[1].hiso_pt[i][k] );
				if ( sampleType(s_name) == 2) hiso_ttj_pt [i][k]->Add( S->isoplots[1].hiso_pt[i][k] );
				if ( sampleType(s_name) == 3) hiso_ewk_pt [i][k]->Add( S->isoplots[1].hiso_pt[i][k] );
				if ( sampleType(s_name) == 4) hiso_rare_pt[i][k]->Add( S->isoplots[1].hiso_pt[i][k] );
				if ( sampleType(s_name) == 5) hiso_db_pt  [i][k]->Add( S->isoplots[1].hiso_pt[i][k] );
			}
			for(int k = 0; k < gNNVrtxBins; ++k){
				if ( sampleType(s_name) == 1) hiso_qcd_nv [i][k]->Add( S->isoplots[1].hiso_nv[i][k] );
				if ( sampleType(s_name) == 2) hiso_ttj_nv [i][k]->Add( S->isoplots[1].hiso_nv[i][k] );
				if ( sampleType(s_name) == 3) hiso_ewk_nv [i][k]->Add( S->isoplots[1].hiso_nv[i][k] );
				if ( sampleType(s_name) == 4) hiso_rare_nv[i][k]->Add( S->isoplots[1].hiso_nv[i][k] );
				if ( sampleType(s_name) == 5) hiso_db_nv  [i][k]->Add( S->isoplots[1].hiso_nv[i][k] );
			}

		}
		hiso_mc_s[i]->Add(hiso_qcd[i]);
		hiso_mc_s[i]->Add(hiso_db[i]);
		hiso_mc_s[i]->Add(hiso_ewk[i]);
		hiso_mc_s[i]->Add(hiso_rare[i]);
		hiso_mc_s[i]->Add(hiso_ttj[i]);
		hiso_mc_s[i]->Draw("goff");
		hiso_mc_s[i]->GetXaxis()->SetTitle(convertVarName("ElRelIso[0]"));

		for(int k = 0; k < gNElPt2bins; ++k){
			hiso_mc_pt_s[i][k]->Add(hiso_qcd_pt[i][k]);
			hiso_mc_pt_s[i][k]->Add(hiso_db_pt[i][k]);
			hiso_mc_pt_s[i][k]->Add(hiso_ewk_pt[i][k]);
			hiso_mc_pt_s[i][k]->Add(hiso_rare_pt[i][k]);
			hiso_mc_pt_s[i][k]->Add(hiso_ttj_pt[i][k]);
			hiso_mc_pt_s[i][k]->Draw("goff");
			hiso_mc_pt_s[i][k]->GetXaxis()->SetTitle(convertVarName("ElRelIso[0]"));
		}

		for(int k = 0; k < gNNVrtxBins; ++k){
			hiso_mc_nv_s[i][k]->Add(hiso_qcd_nv[i][k]);
			hiso_mc_nv_s[i][k]->Add(hiso_db_nv[i][k]);
			hiso_mc_nv_s[i][k]->Add(hiso_ewk_nv[i][k]);
			hiso_mc_nv_s[i][k]->Add(hiso_rare_nv[i][k]);
			hiso_mc_nv_s[i][k]->Add(hiso_ttj_nv[i][k]);
			hiso_mc_nv_s[i][k]->Draw("goff");
			hiso_mc_nv_s[i][k]->GetXaxis()->SetTitle(convertVarName("ElRelIso[0]"));
		}
		//// Fill MC stacks
		//for(size_t j = 0; j < mcsamples.size();   ++j){
		//	Sample *S = fSamples[mcsamples[j]];			
		//	hiso_mc  [i]->Add(S->isoplots[1].hiso[i]);
		//	hiso_mc_s[i]->Add(S->isoplots[1].hiso[i]);
		//	hiso_mc_s[i]->Draw("goff");
		//	hiso_mc_s[i]->GetXaxis()->SetTitle(convertVarName("ElRelIso[0]"));
		//	for(int k = 0; k < gNElPt2bins; ++k){
		//		hiso_mc_pt  [i][k]->Add(S->isoplots[1].hiso_pt[i][k]);
		//		hiso_mc_pt_s[i][k]->Add(S->isoplots[1].hiso_pt[i][k]);
		//		hiso_mc_pt_s[i][k]->Draw("goff");
		//		hiso_mc_pt_s[i][k]->GetXaxis()->SetTitle(convertVarName("ElRelIso[0]"));
		//	}
		//	for(int k = 0; k < gNNVrtxBins; ++k){
		//		hiso_mc_nv  [i][k]->Add(S->isoplots[1].hiso_nv[i][k]);
		//		hiso_mc_nv_s[i][k]->Add(S->isoplots[1].hiso_nv[i][k]);
		//		hiso_mc_nv_s[i][k]->Draw("goff");
		//		hiso_mc_nv_s[i][k]->GetXaxis()->SetTitle(convertVarName("ElRelIso[0]"));
		//	}
		//}


		double max1 = hiso_mc_s[i]->GetMaximum();
		double max2 = hiso_data[i]->GetMaximum();
		double max = max1>max2?max1:max2;

		hiso_mc_s[i]->SetMaximum(1.2*max);
		hiso_data[i]->SetMaximum(1.2*max);

		int bin0   = hiso_data[i]->FindBin(0.0);
		int bin015 = hiso_data[i]->FindBin(0.15) - 1; // bins start at lower edge...
		int bin06  = hiso_data[i]->FindBin(0.6)  - 1;
		float ratio_data  = hiso_data[i] ->Integral(bin0, bin015) / hiso_data[i] ->Integral(bin0, bin06);
		float ratio_mc    = hiso_mc[i]   ->Integral(bin0, bin015) / hiso_mc[i]   ->Integral(bin0, bin06);
		float ratio_ttbar = 0.;
		if(dottbar) hiso_ttbar[i]->Integral(bin0, bin015) / hiso_ttbar[i]->Integral(bin0, bin06);

		TCanvas *c_temp = new TCanvas("ElIso" + IsoPlots::sel_name[i], "Electron Isolation in Data vs MC", 0, 0, 800, 600);
		c_temp->cd();

		TLegend *leg = new TLegend(0.70,0.30,0.95,0.53);
		// TLegend *leg = new TLegend(0.75,0.60,0.89,0.88);
		leg->AddEntry(hiso_data[i], "Data","p");
		if(dottbar) leg->AddEntry(hiso_ttbar[i], "TTbar fake","p");
		//for(size_t j = 0; j < mcsamples.size(); ++j) leg->AddEntry(fSamples[mcsamples[j]]->isoplots[1].hiso[i], fSamples[mcsamples[j]]->sname.Data(), "f");
		leg->AddEntry(hiso_ttj[i],  "Top","f");
		leg->AddEntry(hiso_rare[i], "Rare SM","f");
		leg->AddEntry(hiso_ewk[i],  "Single Boson","f");
		leg->AddEntry(hiso_db[i],   "Di-Boson","f");
		leg->AddEntry(hiso_qcd[i],  "QCD","f");
		leg->SetFillStyle(0);
		leg->SetTextFont(42);
		leg->SetBorderSize(0);

		// gPad->SetLogy();
		hiso_mc_s[i]->Draw("hist");
		if(dottbar) hiso_ttbar[i]->DrawCopy("PE X0 same");
		hiso_data[i]->DrawCopy("PE X0 same");
		leg->Draw();
		lat->DrawLatex(0.70,0.92, Form("L_{int.} = %2.1f fb^{-1}", fLumiNorm/1000.));
		lat->DrawLatex(0.75,0.85, Form("R^{T/L}_{Data} = %4.2f", ratio_data));
		lat->DrawLatex(0.75,0.80, Form("R^{T/L}_{MC}  = %4.2f", ratio_mc));
		lat->SetTextColor(kRed);
		if(dottbar) lat->DrawLatex(0.75,0.75, Form("R^{T/L}_{TTbar} = %4.2f", ratio_ttbar));
		lat->SetTextColor(kBlack);

		// Util::PrintNoEPS(c_temp, "Iso" + IsoPlots::sel_name[i], fOutputDir + fOutputSubDir, NULL);
		Util::PrintPDF(c_temp, "ElIso" + IsoPlots::sel_name[i], fOutputDir + fOutputSubDir);

		for(int k = 0; k < gNElPt2bins; ++k){
			fOutputSubDir = "Isolation/Electrons/PtBinned/";
			double max1 = hiso_mc_pt_s[i][k]->GetMaximum();
			double max2 = hiso_data_pt[i][k]->GetMaximum();
			double max = max1>max2?max1:max2;

			hiso_mc_pt_s[i][k]->SetMaximum(1.2*max);
			hiso_data_pt[i][k]->SetMaximum(1.2*max);
						
			ratio_data = hiso_data_pt[i][k]  ->Integral(bin0, bin015) / hiso_data_pt[i][k] ->Integral(bin0, bin06);
			ratio_mc   = hiso_mc_pt[i][k]    ->Integral(bin0, bin015) / hiso_mc_pt[i][k]   ->Integral(bin0, bin06);
			ratio_ttbar = 0.;
			if(dottbar) ratio_ttbar = hiso_ttbar_pt[i][k]->Integral(bin0, bin015) / hiso_ttbar_pt[i][k]->Integral(bin0, bin06);

			TCanvas *c_temp = new TCanvas(Form("ElIso%s_pt_%d", IsoPlots::sel_name[i].Data(), k), "Electron Isolation in Data vs MC", 0, 0, 800, 600);
			c_temp->cd();

			TLegend *leg_pt = new TLegend(0.70,0.30,0.95,0.53);
			// TLegend *leg_pt = new TLegend(0.75,0.60,0.89,0.88);
			leg_pt->AddEntry(hiso_data_pt[i][k], "Data","p");
			if(dottbar) leg_pt->AddEntry(hiso_ttbar_pt[i][k], "TTbar fake","p");
			//for(size_t j = 0; j < mcsamples.size(); ++j) leg_pt->AddEntry(fSamples[mcsamples[j]]->isoplots[1].hiso_pt[i][k], fSamples[mcsamples[j]]->sname.Data(), "f");
			leg_pt->AddEntry(hiso_ttj_pt  [i][k], "Top","f");
			leg_pt->AddEntry(hiso_rare_pt [i][k], "Rare SM","f");
			leg_pt->AddEntry(hiso_ewk_pt  [i][k], "Single Boson","f");
			leg_pt->AddEntry(hiso_db_pt   [i][k], "Di-Boson","f");
			leg_pt->AddEntry(hiso_qcd_pt  [i][k], "QCD","f");
			leg_pt->SetFillStyle(0);
			leg_pt->SetTextFont(42);
			leg_pt->SetBorderSize(0);

			// gPad->SetLogy();
			hiso_mc_pt_s[i][k]->Draw("hist");
			if(dottbar) hiso_ttbar_pt[i][k]->DrawCopy("PE X0 same");
			hiso_data_pt[i][k]->DrawCopy("PE X0 same");
			leg_pt->Draw();
			lat->DrawLatex(0.20,0.92, Form("p_{T}(e) %3.0f - %3.0f GeV", getPt2Bins(Elec)[k], getPt2Bins(Elec)[k+1]));
			lat->DrawLatex(0.70,0.92, Form("L_{int.} = %2.1f fb^{-1}", fLumiNorm/1000.));
			lat->DrawLatex(0.75,0.85, Form("R^{T/L}_{Data} = %4.2f", ratio_data));
			lat->DrawLatex(0.75,0.80, Form("R^{T/L}_{MC}  = %4.2f", ratio_mc));
			lat->SetTextColor(kRed);
			if(dottbar) lat->DrawLatex(0.75,0.75, Form("R^{T/L}_{TTbar} = %4.2f", ratio_ttbar));
			lat->SetTextColor(kBlack);

			// Util::PrintNoEPS(c_temp, Form("ElIso%s_pt_%d", IsoPlots::sel_name[i].Data(), k), fOutputDir + fOutputSubDir, NULL);
			Util::PrintPDF(c_temp, Form("ElIso%s_pt_%d", IsoPlots::sel_name[i].Data(), k), fOutputDir + fOutputSubDir);
		}
		for(int k = 0; k < gNNVrtxBins; ++k){
			fOutputSubDir = "Isolation/Electrons/NVrtxBinned/";
			double max1 = hiso_mc_nv_s[i][k]->GetMaximum();
			double max2 = hiso_data_nv[i][k]->GetMaximum();
			double max = max1>max2?max1:max2;

			hiso_mc_nv_s[i][k]->SetMaximum(1.2*max);
			hiso_data_nv[i][k]->SetMaximum(1.2*max);
						
			ratio_data = hiso_data_nv[i][k]  ->Integral(bin0, bin015) / hiso_data_nv[i][k] ->Integral(bin0, bin06);
			ratio_mc   = hiso_mc_nv[i][k]    ->Integral(bin0, bin015) / hiso_mc_nv[i][k]   ->Integral(bin0, bin06);
			ratio_ttbar = 0.;
			if(dottbar) ratio_ttbar = hiso_ttbar_nv[i][k]->Integral(bin0, bin015) / hiso_ttbar_nv[i][k]->Integral(bin0, bin06);

			TCanvas *c_temp = new TCanvas(Form("ElIso%s_nv_%d", IsoPlots::sel_name[i].Data(), k), "Electron Isolation in Data vs MC", 0, 0, 800, 600);
			c_temp->cd();

			TLegend *leg_nv = new TLegend(0.70,0.30,0.95,0.53);
			// TLegend *leg_nv = new TLegend(0.75,0.60,0.89,0.88);
			leg_nv->AddEntry(hiso_data_nv[i][k], "Data","p");
			if(dottbar) leg_nv->AddEntry(hiso_ttbar_nv[i][k], "TTbar fake","p");
			//for(size_t j = 0; j < mcsamples.size(); ++j) leg_nv->AddEntry(fSamples[mcsamples[j]]->isoplots[1].hiso_nv[i][k], fSamples[mcsamples[j]]->sname.Data(), "f");
			leg_nv->AddEntry(hiso_ttj_nv  [i][k],  "Top","f");
			leg_nv->AddEntry(hiso_rare_nv [i][k],  "Rare SM","f");
			leg_nv->AddEntry(hiso_ewk_nv  [i][k],  "Single Boson","f");
			leg_nv->AddEntry(hiso_db_nv   [i][k],  "Di-Boson","f");
			leg_nv->AddEntry(hiso_qcd_nv  [i][k],  "QCD","f");
			leg_nv->SetFillStyle(0);
			leg_nv->SetTextFont(42);
			leg_nv->SetBorderSize(0);

			// gPad->SetLogy();
			hiso_mc_nv_s[i][k]->Draw("hist");
			if(dottbar) hiso_ttbar_nv[i][k]->DrawCopy("PE X0 same");
			hiso_data_nv[i][k]->DrawCopy("PE X0 same");
			leg_nv->Draw();
			lat->DrawLatex(0.20,0.92, Form("N_{Vrtx.} %2.0f - %2.0f", gNVrtxBins[k], gNVrtxBins[k+1]));
			lat->DrawLatex(0.70,0.92, Form("L_{int.} = %2.1f fb^{-1}", fLumiNorm/1000.));
			lat->DrawLatex(0.75,0.85, Form("R^{T/L}_{Data} = %4.2f", ratio_data));
			lat->DrawLatex(0.75,0.80, Form("R^{T/L}_{MC}  = %4.2f", ratio_mc));
			lat->SetTextColor(kRed);
			if(dottbar) lat->DrawLatex(0.75,0.75, Form("R^{T/L}_{TTbar} = %4.2f", ratio_ttbar));
			lat->SetTextColor(kBlack);

			// Util::PrintNoEPS(c_temp, Form("ElIso%s_nv_%d", IsoPlots::sel_name[i].Data(), k), fOutputDir + fOutputSubDir, NULL);
			Util::PrintPDF(c_temp, Form("ElIso%s_nv_%d", IsoPlots::sel_name[i].Data(), k), fOutputDir + fOutputSubDir);
		}
	}
}

void SSDLPlotter::makeElIdPlots(){
	char cmd[100];
    sprintf(cmd,"mkdir -p %s%s", fOutputDir.Data(), fOutputSubDir.Data());
    system(cmd);

	vector<int> mcsamples = fMCBG;
	vector<int> datasamples = fEGData;
	
	TH1D    *data [4] [gNSels];
	TH1D    *qcd  [4] [gNSels];
	TH1D    *ttj  [4] [gNSels];
	TH1D    *ewk  [4] [gNSels];
	TH1D    *rare [4] [gNSels];
	TH1D    *db   [4] [gNSels];
	TH1D	*ttbar[4] [gNSels];
	THStack *mc_s [4] [gNSels];

	float intscale[4];
	double datamax[4], mcmax[4], max[4];
	double axismin[4], axismax[4];

	TLatex *lat = new TLatex();
	lat->SetNDC(kTRUE);
	lat->SetTextColor(kBlack);
	lat->SetTextSize(0.04);

	// loop k is over the 4 ID variables: HoE, sigmaIetaIeta, dEta, dPhi
	for (int k =0 ; k<4 ; k++){
		TString varName;
		switch (k) {
			case 0: 
				varName = "ElHoverE";
				axismin[k] = 0.;
				axismax[k] = 0.15;
				break ;
			case 1: 
				varName = "ElSigmaIetaIeta";
				axismin[k] = 0.;
				axismax[k] = 0.035;
				break ;
			case 2: 
				varName = "ElDEta";
				axismin[k] = -0.01;
				axismax[k] =  0.01;
				break ;
			case 3: 
				varName = "ElDPhi";
				axismin[k] = -0.15;
				axismax[k] =  0.15;
				break ;
		}
		for(size_t i = 0; i < gNSels; ++i){
			data [k][i]  = new TH1D("El"+varName+"Data_"          + IdPlots::sel_name[i], "Electron "+varName+" in Data for "  + IdPlots::sel_name[i], IdPlots::nbins[i], axismin[k], axismax[k]);
			qcd  [k][i]  = new TH1D("El"+varName+"MC_"            + IdPlots::sel_name[i], "Electron "+varName+" in MC for "    + IdPlots::sel_name[i], IdPlots::nbins[i], axismin[k], axismax[k]);
			ttj  [k][i]  = new TH1D("El"+varName+"MC_"            + IdPlots::sel_name[i], "Electron "+varName+" in MC for "    + IdPlots::sel_name[i], IdPlots::nbins[i], axismin[k], axismax[k]);
			ewk  [k][i]  = new TH1D("El"+varName+"MC_"            + IdPlots::sel_name[i], "Electron "+varName+" in MC for "    + IdPlots::sel_name[i], IdPlots::nbins[i], axismin[k], axismax[k]);
			rare [k][i]  = new TH1D("El"+varName+"MC_"            + IdPlots::sel_name[i], "Electron "+varName+" in MC for "    + IdPlots::sel_name[i], IdPlots::nbins[i], axismin[k], axismax[k]);
			db   [k][i]  = new TH1D("El"+varName+"MC_"            + IdPlots::sel_name[i], "Electron "+varName+" in MC for "    + IdPlots::sel_name[i], IdPlots::nbins[i], axismin[k], axismax[k]);
			mc_s [k][i]  = new THStack("El"+varName+"MC_stacked_" + IdPlots::sel_name[i], "Electron "+varName+" in MC for "    + IdPlots::sel_name[i]);

			data [k][i]->Sumw2();
			qcd  [k][i]->Sumw2();
			ttj  [k][i]->Sumw2();
			ewk  [k][i]->Sumw2();
			rare [k][i]->Sumw2();
			db   [k][i]->Sumw2();

			qcd  [k][i]->SetFillColor(kYellow-4);
			ttj  [k][i]->SetFillColor(kAzure-5);
			ewk  [k][i]->SetFillColor(kGreen +1);
			rare [k][i]->SetFillColor(kAzure+8);
			db   [k][i]->SetFillColor(kSpring-9);

			fOutputSubDir = "ElectronIDPlots/";
			data[k][i]->SetXTitle(convertVarName(varName));
			data[k][i]->SetLineWidth(3);
			data[k][i]->SetLineColor(kBlack);
			data[k][i]->SetMarkerStyle(8);
			data[k][i]->SetMarkerColor(kBlack);
			data[k][i]->SetMarkerSize(1.2);

			// Apply weights to MC histos
			for(size_t j = 0; j < gNSAMPLES; ++j){
				Sample *S = fSamples[j];
				float lumiscale = fLumiNorm / S->lumi;
				if(S->datamc == 0) continue;
				switch (k) {
					case 0: S->idplots.hhoe[i]->Scale(lumiscale); break;
					case 1: S->idplots.hsiesie[i]->Scale(lumiscale); break;
					case 2: S->idplots.hdeta[i]->Scale(lumiscale); break;
					case 3: S->idplots.hdphi[i]->Scale(lumiscale); break;
				}
			}

			// Fill data histo
			for(size_t j = 0; j < datasamples.size(); ++j){
				Sample *S = fSamples[datasamples[j]];
				switch (k) {
					case 0: data[k][i]->Add(S->idplots.hhoe[i]); break;
					case 1: data[k][i]->Add(S->idplots.hsiesie[i]); break;
					case 2: data[k][i]->Add(S->idplots.hdeta[i]); break;
					case 3: data[k][i]->Add(S->idplots.hdphi[i]); break;
				}
			}

			// Scale to get equal integrals
			intscale[k] = 0.;
			for(size_t j = 0; j < mcsamples.size();   ++j){
				Sample *S = fSamples[mcsamples[j]];
				switch (k) {
					case 0: intscale[k] += S->idplots.hhoe[i]->Integral(); break;
					case 1: intscale[k] += S->idplots.hsiesie[i]->Integral(); break;
					case 2: intscale[k] += S->idplots.hdeta[i]->Integral(); break;
					case 3: intscale[k] += S->idplots.hdphi[i]->Integral(); break;
				}
			}
			intscale[k] = data[k][i]->Integral() / intscale[k];
			
			for(size_t j = 0; j < mcsamples.size();   ++j){
				Sample *S = fSamples[mcsamples[j]];			
				switch (k) {
					case 0: S -> idplots.hhoe[i]    -> Scale(intscale[k]); break;
					case 1: S -> idplots.hsiesie[i] -> Scale(intscale[k]); break;
					case 2: S -> idplots.hdeta[i]   -> Scale(intscale[k]); break;
					case 3: S -> idplots.hdphi[i]   -> Scale(intscale[k]); break;
				}
			}

			// Fill MC stacks
			for(size_t j = 0; j < mcsamples.size();   ++j){
				Sample *S = fSamples[mcsamples[j]];			
				TString s_name = S->sname;
				TH1D * histo;
				switch (k) {
					case 0: histo = S->idplots.hhoe[i]   ; break;
					case 1: histo = S->idplots.hsiesie[i]; break;
					case 2: histo = S->idplots.hdeta[i]  ; break;
					case 3: histo = S->idplots.hdphi[i]  ; break;
				}
				if ( sampleType(s_name) == 1) qcd [k][i]->Add( histo );
				if ( sampleType(s_name) == 2) ttj [k][i]->Add( histo );
				if ( sampleType(s_name) == 3) ewk [k][i]->Add( histo );
				if ( sampleType(s_name) == 4) rare[k][i]->Add( histo );
				if ( sampleType(s_name) == 5) db  [k][i]->Add( histo );
				delete histo;
			} // end loop over MC samples

			mc_s[k][i]->Add(qcd [k][i]);
			mc_s[k][i]->Add(db  [k][i]);
			mc_s[k][i]->Add(ewk [k][i]);
			mc_s[k][i]->Add(rare[k][i]);
			mc_s[k][i]->Add(ttj [k][i]);
			mc_s[k][i]->Draw("goff");
			mc_s[k][i]->GetXaxis()->SetTitle(convertVarName(varName));

			datamax[k] = data[k][i]->GetMaximum();
			datamax[k] = mc_s[k][i]->GetMaximum();
			max[k] = datamax[k] > mcmax[k] ? datamax[k] : mcmax[k];

			mc_s[k][i]->SetMaximum(1.2*max[k]);
			data[k][i]->SetMaximum(1.2*max[k]);

			// H over E
			TCanvas *temp = new TCanvas(varName + IdPlots::sel_name[i], "Electron "+varName+" in Data vs MC", 0, 0, 800, 600);
			temp->cd();

			TLegend *leg = new TLegend(0.70,0.65,0.89,0.88);
			leg->AddEntry(data[k][i], "Data","p");
			leg->AddEntry(ttj[k][i],  "Top","f");
			leg->AddEntry(rare[k][i], "Rare SM","f");
			leg->AddEntry(ewk[k][i],  "Single Boson","f");
			leg->AddEntry(db[k][i],   "Di-Boson","f");
			leg->AddEntry(qcd[k][i],  "QCD","f");

			leg->SetFillStyle(0);
			leg->SetTextFont(42);
			leg->SetBorderSize(0);

			gPad->SetLogy();
			mc_s[k][i]->Draw("hist");
			data[k][i]->DrawCopy("PE X0 same");
			leg->Draw();
			lat->DrawLatex(0.70,0.92, Form("L_{int.} = %2.1f fb^{-1}", fLumiNorm/1000.));
			lat->SetTextColor(kRed);
			//lat->DrawLatex(0.75,0.75, Form("R^{T/L}_{TTbar} = %4.2f", ratio_ttbar));
			lat->SetTextColor(kBlack);

			Util::PrintPDF(temp, varName + IdPlots::sel_name[i], fOutputDir + fOutputSubDir);
			delete temp, leg;
		} // end loop over region selections
	} //end loop over ID variables
} //end function

inline int SSDLPlotter::sampleType(TString s_name){
	// returns an integer corresponding to the MC sample type:
	// QCD = 1, Top = 2, EWK = 3, RARE = 4 and DiBoson = 5
	// if none matches, returns 0
	if ( (s_name.Contains("QCD")) ||
	     (s_name) == "MuEnr10" )
	return 1;
	if ( (s_name.Contains("SingleT")) ||
	     (s_name) == "TTJets" )
	return 2;
	if ( (s_name.Contains("DYJets")) ||
	     (s_name.Contains("GJets"))  ||
	     (s_name) == "WJets" )
	return 3;
	if ( (s_name) == "TTbarW"    ||
	     (s_name) == "TTbarZ"    ||
	     (s_name) == "TTbarG"    ||
	     (s_name) == "DPSWW"     ||
	     (s_name) == "WWZ"       ||
	     (s_name) == "WZZ"       ||
	     (s_name) == "WZZ"       ||
	     (s_name) == "WWG"       ||
	     (s_name) == "ZZZ"       ||
	     (s_name) == "WWW"       ||
	     (s_name) == "W+W+"      ||
	     (s_name) == "W-W-")
	return 4;
	if ( (s_name.Contains("GVJets"))    ||
	     (s_name.Contains("WWTo2L2Nu")) ||
	     (s_name.Contains("WZTo3LNu"))  ||
	     (s_name.Contains("ZZTo4L")) )
	return 5;
	else {
		cout << "ERROR: This should not happen! Check the SSDLPlotter::sampleType() function to have "<< s_name << "in one of the lists!" << endl;
		return 0;
	}
}

void SSDLPlotter::makeNT2KinPlots(gHiLoSwitch hilo){
	TString selname[3] = {"LooseLoose", "TightTight", "Signal"};

	for(size_t s = 0; s < 3; ++s){ // loop on selections
		fOutputSubDir = "KinematicPlots/" + gHiLoLabel[hilo] + "/" + selname[s];
		char cmd[100];
	    sprintf(cmd,"mkdir -p %s%s", fOutputDir.Data(), fOutputSubDir.Data());
	    system(cmd);
		
		TH1D    *hvar_data[gNKinVars];

		TH1D    *hvar_qcd  [gNKinVars];
		TH1D    *hvar_ttj  [gNKinVars];
		TH1D    *hvar_ewk  [gNKinVars];
		TH1D    *hvar_rare [gNKinVars];
		TH1D    *hvar_db   [gNKinVars];

		THStack *hvar_mc_s[gNKinVars];

		TLatex *lat = new TLatex();
		lat->SetNDC(kTRUE);
		lat->SetTextColor(kBlack);
		lat->SetTextSize(0.04);

		// Create histograms
		for(size_t i = 0; i < gNKinVars; ++i){
			hvar_data[i] = new TH1D("Data_"          + KinPlots::var_name[i], KinPlots::var_name[i] + " in Data", KinPlots::nbins[i], KinPlots::xmin[i], KinPlots::xmax[i]);

			hvar_qcd [i] = new TH1D("QCD_"           + KinPlots::var_name[i], KinPlots::var_name[i] + " in MC", KinPlots::nbins[i], KinPlots::xmin[i], KinPlots::xmax[i]);
			hvar_ttj [i] = new TH1D("TTjets_"        + KinPlots::var_name[i], KinPlots::var_name[i] + " in MC", KinPlots::nbins[i], KinPlots::xmin[i], KinPlots::xmax[i]);
			hvar_ewk [i] = new TH1D("EWK_"           + KinPlots::var_name[i], KinPlots::var_name[i] + " in MC", KinPlots::nbins[i], KinPlots::xmin[i], KinPlots::xmax[i]);
			hvar_rare[i] = new TH1D("Rare_"          + KinPlots::var_name[i], KinPlots::var_name[i] + " in MC", KinPlots::nbins[i], KinPlots::xmin[i], KinPlots::xmax[i]);
			hvar_db[i]   = new TH1D("DB_"            + KinPlots::var_name[i], KinPlots::var_name[i] + " in MC", KinPlots::nbins[i], KinPlots::xmin[i], KinPlots::xmax[i]);

			hvar_mc_s[i] = new THStack("MC_stacked_" + KinPlots::var_name[i], KinPlots::var_name[i] + " in MC");
		}

		// Adjust overflow bins:
		for(size_t i = 0; i < gNKinVars; ++i){
			for(gSample j = sample_begin; j < gNSAMPLES; j=gSample(j+1)){
				Int_t nbins     = fSamples[j]->kinplots[s][hilo].hvar[i]->GetNbinsX();
				Double_t binc   = fSamples[j]->kinplots[s][hilo].hvar[i]->GetBinContent(nbins);
				Double_t overfl = fSamples[j]->kinplots[s][hilo].hvar[i]->GetBinContent(nbins+1);
				fSamples[j]->kinplots[s][hilo].hvar[i]->SetBinContent(nbins, binc + overfl);
			}
		}

		vector<int> mcsamples   = fMCBGMuEnr;
		vector<int> datasamples = fHighPtData;
		//////////////////////////////////////////////////////////
		// Make kin plots
		for(size_t i = 0; i < gNKinVars; ++i){
			// Create plots
			if(i != 6) mcsamples = fMCBG;

			hvar_data[i]->SetXTitle(KinPlots::axis_label[i]);
			hvar_data[i]->SetLineWidth(3);
			hvar_data[i]->SetLineColor(kBlack);
			hvar_data[i]->SetMarkerStyle(8);
			hvar_data[i]->SetMarkerColor(kBlack);
			hvar_data[i]->SetMarkerSize(1.2);

			// Scale by luminosity
			for(size_t j = 0; j < gNSAMPLES; ++j){
				float lumiscale = fLumiNorm / fSamples[j]->lumi;
				if(fSamples[j]->datamc == 0) continue;
				fSamples[j]->kinplots[s][hilo].hvar[i]->Scale(lumiscale);
			}

			// Fill data histo
			for(size_t j = 0; j < datasamples.size(); ++j){
				Sample *S = fSamples[datasamples[j]];
				hvar_data[i]->Add(S->kinplots[s][hilo].hvar[i]);
				hvar_data[i]->SetXTitle(KinPlots::axis_label[i]);
			}

			// Scale to get equal integrals
			float intscale(0.);
			for(size_t j = 0; j < mcsamples.size();   ++j){
				Sample *S = fSamples[mcsamples[j]];
				intscale += S->kinplots[s][hilo].hvar[i]->Integral();
			}
			intscale = hvar_data[i]->Integral() / intscale;
			
			for(size_t j = 0; j < mcsamples.size();   ++j){
				Sample *S = fSamples[mcsamples[j]];
				S->kinplots[s][hilo].hvar[i]->Scale(intscale);
			}

			hvar_qcd [i]->SetFillColor(kYellow-4);
			hvar_db  [i]->SetFillColor(kSpring-9);
			hvar_ewk [i]->SetFillColor(kGreen +1);
			hvar_ttj [i]->SetFillColor(kAzure-5);
			hvar_rare[i]->SetFillColor(kAzure+8);
			// hvar_ttj [i]->SetFillColor(kAzure +1);
			// hvar_rare[i]->SetFillColor(kViolet+5);


			// Fill MC stacks
			for(size_t j = 0; j < mcsamples.size();   ++j){
				Sample *S = fSamples[mcsamples[j]];
				TString s_name = S->sname;
				// sampleType Function: QCD = 1 , Top = 2, EWK = 3 , Rare = 4 , DB = 5
				if ( sampleType(s_name) == 1) hvar_qcd [i]->Add( S->kinplots[s][hilo].hvar[i] );
				if ( sampleType(s_name) == 2) hvar_ttj [i]->Add( S->kinplots[s][hilo].hvar[i] );
				if ( sampleType(s_name) == 3) hvar_ewk [i]->Add( S->kinplots[s][hilo].hvar[i] );
				if ( sampleType(s_name) == 4) hvar_rare[i]->Add( S->kinplots[s][hilo].hvar[i] );
				if ( sampleType(s_name) == 5) hvar_db  [i]->Add( S->kinplots[s][hilo].hvar[i] );
			}
			hvar_mc_s[i]->Add(hvar_qcd[i]);
			hvar_mc_s[i]->Add(hvar_db[i]);
			hvar_mc_s[i]->Add(hvar_ewk[i]);
			hvar_mc_s[i]->Add(hvar_rare[i]);
			hvar_mc_s[i]->Add(hvar_ttj[i]);
			hvar_mc_s[i]->Draw("goff");
			hvar_mc_s[i]->GetXaxis()->SetTitle(KinPlots::axis_label[i]);

			double max1 = hvar_mc_s[i]->GetMaximum();
			double max2 = hvar_data[i]->GetMaximum();
			double max = max1>max2?max1:max2;
			// hvar_mc_s[i]->SetMaximum(1.5*max);
			// hvar_data[i]->SetMaximum(1.5*max);
			// hvar_mc_s[i]->SetMinimum(0.);
			// hvar_data[i]->SetMinimum(0.);
			hvar_mc_s[i]->SetMaximum(5.*max);
			hvar_data[i]->SetMaximum(5.*max);
			hvar_mc_s[i]->SetMinimum(0.5);
			hvar_data[i]->SetMinimum(0.5);

			TCanvas *c_temp = new TCanvas("C_" + KinPlots::var_name[i], KinPlots::var_name[i] + " in Data vs MC", 0, 0, 800, 600);
			c_temp->cd();

			// TLegend *leg = new TLegend(0.15,0.50,0.40,0.88);
			// TLegend *leg = new TLegend(0.70,0.30,0.90,0.68);
			TLegend *leg = new TLegend(0.70,0.65,0.89,0.88);
			leg->AddEntry(hvar_data[i], "Data","p");
			leg->AddEntry(hvar_ttj[i],  "Top","f");
			leg->AddEntry(hvar_rare[i], "Rare SM","f");
			leg->AddEntry(hvar_ewk[i],  "Single Boson","f");
			leg->AddEntry(hvar_db[i],   "Di-Boson","f");
			leg->AddEntry(hvar_qcd[i],  "QCD","f");
	
			// for(size_t j = 0; j < mcsamples.size(); ++j){ 
			// 	leg->AddEntry(fSamples[mcsamples[j]]->kinplots[s][hilo].hvar[i], fSamples[mcsamples[j]]->sname.Data(), "f");
			// }
			leg->SetFillStyle(0);
			leg->SetTextFont(42);
			leg->SetBorderSize(0);

			gPad->SetLogy();
			hvar_mc_s[i]->Draw("hist");
			hvar_data[i]->DrawCopy("PE X0 same");
			leg->Draw();
			
			drawTopLine();
			// lat->DrawLatex(0.70,0.92, Form("L_{int.} = %2.1f fb^{-1}", fLumiNorm/1000.));
			// lat->DrawLatex(0.11,0.92, selname[s]);

			if(i < 5)  lat->DrawLatex(0.14,0.85, "ee/e#mu/#mu#mu");
			if(i == 5) lat->DrawLatex(0.14,0.85, "ee/#mu#mu");
			if(i == 6) lat->DrawLatex(0.14,0.85, "#mu#mu");
			if(i == 7) lat->DrawLatex(0.14,0.85, "ee");
			if(i == 8) lat->DrawLatex(0.14,0.85, "e#mu");
			if(i > 8)  lat->DrawLatex(0.14,0.85, "ee/e#mu/#mu#mu");

			// Util::PrintNoEPS(c_temp, KinPlots::var_name[i], fOutputDir + fOutputSubDir, NULL);
			Util::PrintPDF(c_temp, KinPlots::var_name[i], fOutputDir + fOutputSubDir);
			delete c_temp;
			delete leg;
		}
	}	
}
void SSDLPlotter::makeMETvsHTPlot(vector<int> mmsamples, vector<int> eesamples, vector<int> emsamples, gHiLoSwitch hilo){
	if(readSigGraphs(fOutputFileName) != 0) return;
	TString hiloname[2] = {"p_{T}(l_{1}/l_{2}) > 20/10 GeV", "p_{T}(#mu/e) > 5/10 GeV"};

	fOutputSubDir = "KinematicPlots/HTvsMET/";
	char cmd[100];
    sprintf(cmd,"mkdir -p %s%s", fOutputDir.Data(), fOutputSubDir.Data());
    system(cmd);
	
	const float htmax = 1300.;
	const float metmax = 250.;

	// Create histograms
	TH2D *hmetvsht_da_mm = new TH2D("Data_HTvsMET_mm", "Data_HTvsMET_mm", 100, 0., htmax, 100, 0., metmax);
	TH2D *hmetvsht_da_ee = new TH2D("Data_HTvsMET_ee", "Data_HTvsMET_ee", 100, 0., htmax, 100, 0., metmax);
	TH2D *hmetvsht_da_em = new TH2D("Data_HTvsMET_em", "Data_HTvsMET_em", 100, 0., htmax, 100, 0., metmax);
	// TH2D *hmetvsht_da_mm = new TH2D("Data_HTvsMET_mm", "Data_HTvsMET_mm", KinPlots::nHTBins, KinPlots::HTmin, KinPlots::HTmax, KinPlots::nMETBins, KinPlots::METmin, KinPlots::METmax);
	// TH2D *hmetvsht_da_ee = new TH2D("Data_HTvsMET_ee", "Data_HTvsMET_ee", KinPlots::nHTBins, KinPlots::HTmin, KinPlots::HTmax, KinPlots::nMETBins, KinPlots::METmin, KinPlots::METmax);
	// TH2D *hmetvsht_da_em = new TH2D("Data_HTvsMET_em", "Data_HTvsMET_em", KinPlots::nHTBins, KinPlots::HTmin, KinPlots::HTmax, KinPlots::nMETBins, KinPlots::METmin, KinPlots::METmax);

	// vector<int> mcsamples   = fMCBGMuEnr;
	// const gSample sig = LM1;

	//////////////////////////////////////////////////////////
	// Make MET vs HT plot:
	hmetvsht_da_mm->SetMarkerStyle(8);
	hmetvsht_da_mm->SetMarkerColor(kBlack);
	hmetvsht_da_mm->SetMarkerSize(1.5);
	hmetvsht_da_mm->GetYaxis()->SetTitleOffset(1.4);

	hmetvsht_da_ee->SetMarkerStyle(21);
	hmetvsht_da_ee->SetMarkerColor(kRed);
	hmetvsht_da_ee->SetMarkerSize(1.4);
	hmetvsht_da_ee->GetYaxis()->SetTitleOffset(1.4);

	hmetvsht_da_em->SetMarkerStyle(23);
	hmetvsht_da_em->SetMarkerColor(kBlue);
	hmetvsht_da_em->SetMarkerSize(1.7  );
	hmetvsht_da_em->GetYaxis()->SetTitleOffset(1.4);

	// TGraphs:
	TMultiGraph *gmetvsht_da_mm = new TMultiGraph("HTvsMET_mm", "HTvsMET_mm");
	TMultiGraph *gmetvsht_da_ee = new TMultiGraph("HTvsMET_ee", "HTvsMET_ee");
	TMultiGraph *gmetvsht_da_em = new TMultiGraph("HTvsMET_em", "HTvsMET_em");

	// Fill data histo
	// for(size_t i = 0; i < mmsamples.size(); ++i) hmetvsht_da_mm->Add(fSamples[mmsamples[i]]->kinplots[s][hilo].hmetvsht);
	// for(size_t i = 0; i < eesamples.size(); ++i) hmetvsht_da_ee->Add(fSamples[eesamples[i]]->kinplots[s][hilo].hmetvsht);
	// for(size_t i = 0; i < emsamples.size(); ++i) hmetvsht_da_em->Add(fSamples[emsamples[i]]->kinplots[s][hilo].hmetvsht);

	for(size_t i = 0; i < mmsamples.size(); ++i) gmetvsht_da_mm->Add(fSamples[mmsamples[i]]->sigevents[Muon][hilo]);
	for(size_t i = 0; i < eesamples.size(); ++i) gmetvsht_da_ee->Add(fSamples[eesamples[i]]->sigevents[Elec][hilo]);
	for(size_t i = 0; i < emsamples.size(); ++i) gmetvsht_da_em->Add(fSamples[emsamples[i]]->sigevents[ElMu][hilo]);

	hmetvsht_da_mm->SetXTitle(KinPlots::axis_label[0]);
	hmetvsht_da_mm->SetYTitle(KinPlots::axis_label[1]);
	hmetvsht_da_ee->SetXTitle(KinPlots::axis_label[0]);
	hmetvsht_da_ee->SetYTitle(KinPlots::axis_label[1]);
	hmetvsht_da_em->SetXTitle(KinPlots::axis_label[0]);
	hmetvsht_da_em->SetYTitle(KinPlots::axis_label[1]);

	TLegend *leg = new TLegend(0.80,0.70,0.95,0.88);
	leg->AddEntry(hmetvsht_da_mm, "#mu#mu","p");
	leg->AddEntry(hmetvsht_da_ee, "ee","p");
	leg->AddEntry(hmetvsht_da_em, "e#mu","p");
	leg->SetFillStyle(0);
	leg->SetTextFont(42);
	leg->SetTextSize(0.05);
	leg->SetBorderSize(0);

	// Special effects:
	const float lowerht = hilo==HighPt? 80.:200.;
	// const float lowerht = hilo==HighPt? 0.:200.;

	TWbox *lowhtbox  = new TWbox(0., 0., lowerht, metmax, kBlack, 0, 0);
	TWbox *lowmetbox = new TWbox(lowerht, 0., htmax,    30., kBlack, 0, 0);
	lowhtbox ->SetFillColor(12);
	lowmetbox->SetFillColor(12);
	lowhtbox ->SetFillStyle(3005);
	lowmetbox->SetFillStyle(3005);
	TLine *boxborder1 = new TLine(lowerht,30.,lowerht,metmax);
	TLine *boxborder2 = new TLine(lowerht,30.,htmax,30.);
	boxborder1->SetLineWidth(1);
	boxborder2->SetLineWidth(1);
	boxborder1->SetLineColor(14);
	boxborder2->SetLineColor(14);

	TLine *sig1x = new TLine(lowerht, 120., lowerht, metmax); // met 120 ht 80
	TLine *sig1y = new TLine(lowerht, 120., htmax,   120.); 
	TLine *sig2x = new TLine(200., 120., 200.,  metmax);      // met 120 ht 200
	TLine *sig2y = new TLine(200., 120., htmax, 120.);
	TLine *sig3x = new TLine(450.,  50., 450.,  metmax);      // met 50  ht 450
	TLine *sig3y = new TLine(450.,  50., htmax,  50.);
	TLine *sig4x = new TLine(450., 120., 450.,  metmax);      // met 120 ht 450
	TLine *sig4y = new TLine(450., 120., htmax, 120.);

	sig1x->SetLineWidth(2);
	sig1y->SetLineWidth(2);
	sig2x->SetLineWidth(2);
	sig2y->SetLineWidth(2);
	sig3x->SetLineWidth(1);
	sig3y->SetLineWidth(1);
	sig4x->SetLineWidth(3);
	sig4y->SetLineWidth(3);

	sig1x->SetLineStyle(3);
	sig1y->SetLineStyle(3);
	sig2x->SetLineStyle(2);
	sig2y->SetLineStyle(2);
	sig3x->SetLineStyle(1);
	sig3y->SetLineStyle(1);
	// sig4x->SetLineStyle(1);
	// sig4y->SetLineStyle(1);

	// float legymax = hilo==HighPt?0.54:0.50;
	// TLegend *regleg = new TLegend(0.70,0.37,0.88,0.54);
	// TLegend *regleg = new TLegend(0.70,0.45,0.88,0.62);
	TLegend *regleg = new TLegend(0.70,0.29,0.88,0.46);
	regleg->AddEntry(sig4x, "Search Region 1","l");
	regleg->AddEntry(sig2x, "Search Region 2","l");
	regleg->AddEntry(sig3x, "Search Region 3","l");
	if(hilo != LowPt) regleg->AddEntry(sig1x, "Search Region 4","l");
	regleg->SetFillStyle(0);
	regleg->SetTextFont(42);
	regleg->SetTextSize(0.03);
	regleg->SetBorderSize(0);
	

	TCanvas *c_temp = new TCanvas("C_HTvsMET", "HT vs MET in Data vs MC", 0, 0, 600, 600);
	c_temp->cd();
	c_temp->SetRightMargin(0.03);
	c_temp->SetLeftMargin(0.13);

	hmetvsht_da_mm->DrawCopy("axis");

	lowhtbox ->Draw();
	lowmetbox->Draw();
	boxborder1->Draw();
	boxborder2->Draw();

	if(hilo != LowPt) sig1x->Draw();
	if(hilo != LowPt) sig1y->Draw();
	sig2x->Draw();
	sig2y->Draw();
	sig3x->Draw();
	sig3y->Draw();
	sig4x->Draw();
	sig4y->Draw();

	// Graphs
	gmetvsht_da_ee->Draw("P");
	gmetvsht_da_em->Draw("P");
	gmetvsht_da_mm->Draw("P");
	
	leg->Draw();
	regleg->Draw();
	drawTopLine();
	TPaveText *pave = new TPaveText(0.16, 0.83, 0.53, 0.88, "NDC");
	pave->SetFillColor(0);
	pave->SetFillStyle(1001);
	pave->SetBorderSize(0);
	pave->SetMargin(0.05);
	pave->SetTextFont(42);
	pave->SetTextSize(0.04);
	pave->SetTextAlign(12);
	pave->AddText(hiloname[hilo]);
	pave->Draw();
	gPad->RedrawAxis();

	// Util::PrintNoEPS(c_temp, "HTvsMET_" + gHiLoLabel[hilo], fOutputDir + fOutputSubDir, NULL);
	Util::PrintPDF(c_temp, "HTvsMET_" + gHiLoLabel[hilo], fOutputDir + fOutputSubDir);
	// Util::SaveAsMacro(c_temp, "HTvsMET_" + gHiLoLabel[hilo], fOutputDir + fOutputSubDir);
	delete c_temp;
	delete leg, regleg;
	delete hmetvsht_da_mm, hmetvsht_da_ee, hmetvsht_da_em;//, hmetvsht_mc;
	delete gmetvsht_da_mm, gmetvsht_da_ee, gmetvsht_da_em;//, hmetvsht_mc;
}
void SSDLPlotter::makeMETvsHTPlotPRL(){
	if(readSigGraphs(fOutputFileName) != 0) return;
	gHiLoSwitch hilo = HighPt;
	TString hiloname[2] = {"p_{T}(l_{1}/l_{2}) > 20/10 GeV", "p_{T}(#mu/e) > 5/10 GeV"};
	fOutputSubDir = "PRL";
	char cmd[100];
    sprintf(cmd,"mkdir -p %s%s", fOutputDir.Data(), fOutputSubDir.Data());
    system(cmd);
	
	const float htmax = 1500.;
	const float metmax = 250.;
	
	Color_t col_mm = kBlack;
	Color_t col_ee = kRed;
	Color_t col_em = kBlue;
	// Color_t col_mt = kMagenta;
	// Color_t col_et = kGreen;
	// Color_t col_tt = kCyan;
	// Color_t col_mm = kAzure   +2;
	// Color_t col_ee = kPink    +2;
	// Color_t col_em = kOrange  +2;
	Color_t col_mt = kGreen;
	Color_t col_et = kMagenta;
	Color_t col_tt = kCyan;

	// Create histograms
	TH2D *hmetvsht_da_mm = new TH2D("Data_HTvsMET_mm", "Data_HTvsMET_mm", 100, 0., htmax, 100, 0., metmax);
	TH2D *hmetvsht_da_ee = new TH2D("Data_HTvsMET_ee", "Data_HTvsMET_ee", 100, 0., htmax, 100, 0., metmax);
	TH2D *hmetvsht_da_em = new TH2D("Data_HTvsMET_em", "Data_HTvsMET_em", 100, 0., htmax, 100, 0., metmax);

	//////////////////////////////////////////////////////////
	// Make MET vs HT plot:
	hmetvsht_da_mm->SetMarkerStyle(8);
	hmetvsht_da_mm->SetMarkerColor(col_mm);
	hmetvsht_da_mm->SetMarkerSize(1.5);
	hmetvsht_da_mm->GetYaxis()->SetTitleOffset(1.4);

	hmetvsht_da_ee->SetMarkerStyle(21);
	hmetvsht_da_ee->SetMarkerColor(col_ee);
	hmetvsht_da_ee->SetMarkerSize(1.4);
	hmetvsht_da_ee->GetYaxis()->SetTitleOffset(1.4);

	hmetvsht_da_em->SetMarkerStyle(23);
	hmetvsht_da_em->SetMarkerColor(col_em);
	hmetvsht_da_em->SetMarkerSize(1.7  );
	hmetvsht_da_em->GetYaxis()->SetTitleOffset(1.4);

///////////// HARDCODED

	// Updated numbers from Ronny with 4.7/fb from Dec 14
	const int nmmev = 83;
	const int nemev = 89;
	const int neeev = 29;
	float ht_mm [nmmev] = {435.9, 436.4, 204.86, 209.26, 218.65, 219.79, 221.71, 222.29, 259.24, 270.02, 489.07, 1070.04, 200.766, 202.494, 204.392, 204.452, 206.212, 208.623, 208.657, 210.096, 211.624, 212.403, 213.249, 213.868, 214.341, 217.729, 220.158, 222.854, 224.574, 226.234, 226.762, 228.338, 228.398, 229.865, 235.782, 245.168, 249.944, 252.454, 253.169, 253.518, 254.717, 256.236, 257.237, 257.834, 258.778, 262.353, 264.248, 267.277, 276.444, 277.564, 279.782, 280.798, 284.485, 284.787, 284.968, 291.472, 293.045, 294.471, 299.879, 314.186, 315.412, 316.189, 318.984, 319.601, 320.678, 331.128, 334.147, 334.864, 354.292, 358.104, 369.434, 372.221, 375.212, 382.892, 393.964, 400.209, 432.051, 455.478, 539.545, 585.687, 587.099, 592.905, 646.496 };
	float met_mm[nmmev] = {62.904, 203, 78.257, 58.0265, 39.8026, 63.0605, 48.5812, 66.7972, 40.6448, 70.0999, 57.5302, 62.992, 74.0824, 95.5969, 70.9175, 44.4644, 83.0277, 30.0451, 51.7086, 46.5073, 53.1209, 30.4686, 112.292, 37.6446, 39.8278, 45.5598, 78.2578, 126.779, 76.8364, 47.6917, 30.5832, 114.587, 48.3581, 30.1828, 37.4236, 58.5595, 205.792, 34.5121, 122.463, 44.8432, 72.3296, 71.6319, 34.209, 70.4086, 53.2947, 43.4222, 120.754, 107.896, 35.0554, 30.6752, 35.0923, 73.8093, 38.8426, 49.9091, 31.0167, 40.1132, 143.384, 37.0943, 62.2534, 59.5283, 64.3967, 35.416, 31.1259, 53.6233, 41.3591, 33.9659, 30.4681, 157.391, 42.3864, 75.1158, 43.2589, 51.0459, 30.7741, 106.485, 186.437, 73.8847, 69.9304, 62.8467, 150.782, 35.0852, 44.8593, 58.09, 143.069};

	float ht_em [nemev] = {304.5, 246.19, 258.43, 261.28, 324.42, 360.08, 417.45, 426.82, 1166.17 , 200.856, 201.015, 203.856, 215.111, 217.713, 224.065, 224.366, 230.359, 230.402, 230.636, 231.013, 231.351, 231.549, 232.082, 235.257, 235.387, 236.783, 239.412, 243.326, 247.511, 247.792, 253.246, 253.919, 255.472, 257.601, 259.813, 260.531, 264.287, 271.378, 273.446, 274.043, 275.092, 279.396, 280.382, 280.612, 280.877, 284.684, 284.782, 286.638, 290.427, 290.471, 291.739, 296.215, 299.392, 300.756, 303.904, 304.033, 316.425, 318.618, 320.556, 320.681, 321.842, 323.792, 325.135, 331.484, 346.067, 347.603, 356.005, 358.059, 359.211, 374.363, 377.818, 377.982, 379.805, 392.301, 397.064, 405.976, 411.237, 420.223, 421.674, 428.888, 433.076, 458.401, 482.122, 493.013, 520.802, 573.502, 616.057, 641.947, 687.421 };
	float met_em[nemev] = {36.6628, 113.665, 53.0968, 72.1098, 56.9286, 57.1148, 86.4964, 52.7702, 148.799, 49.1676, 96.2492, 43.4862, 116.699, 83.6971, 138.175, 45.47, 46.4088, 100.221, 40.2768, 45.8226, 73.4375, 188.366, 86.1323, 52.378, 69.2557, 96.2093, 64.9008, 37.4178, 64.1492, 53.2864, 70.9379, 62.5218, 109.433, 57.4553, 109.035, 105.571, 60.1917, 42.4493, 34.1574, 36.9665, 65.4461, 78.0987, 84.8157, 73.4878, 100.28, 93.9971, 38.5, 135.902, 53.6195, 64.6249, 53.5031, 150.137, 83.4427, 36.4459, 44.6219, 177.228, 95.2148, 61.8755, 38.5762, 115.654, 80.1506, 67.6595, 124.552, 42.7002, 33.9257, 77.7492, 40.6771, 171.337, 62.6187, 54.3684, 121.968, 35.6684, 77.5705, 101.652, 134.948, 133.954, 50.0136, 123.933, 33.6699, 55.3003, 38.1614, 68.5016, 120.909, 80.8768, 32.2697, 172.382, 57.0807, 57.2559, 92.6551};

	float ht_ee [neeev] = {264.,  444., 795.,  204.7, 215.7,  227.6, 229.6, 236.6, 242.1,  246.6, 250.1, 253.8, 255.4, 256.6, 257.6, 268.4, 279.1, 279.6, 294.2, 294.9, 298.6, 356.9, 362.4, 415.1, 481.2,  498.4, 506.1, 561.3, 563.7};
	float met_ee[neeev] = {59.5, 111.9, 42.25, 39.28, 39.23, 146.6, 203.3,  37.62, 66.36, 108.3,  60.66, 46.39, 41.77, 76.5,  64.3, 196.7,  30.48, 36.88, 56.58, 34.11, 33.43, 93.17, 34.58, 37.06, 61.27, 152.5,  53.51, 30.08, 83.83};

	// TGraphs:
	TGraph *gmetvsht_da_mm_lowpt = new TGraph(nmmev, ht_mm, met_mm);
	TGraph *gmetvsht_da_em_lowpt = new TGraph(nemev, ht_em, met_em);
	TGraph *gmetvsht_da_ee_lowpt = new TGraph(neeev, ht_ee, met_ee);
	gmetvsht_da_mm_lowpt->SetName("HTvsMET_mm_lowpt");
	gmetvsht_da_em_lowpt->SetName("HTvsMET_em_lowpt");
	gmetvsht_da_ee_lowpt->SetName("HTvsMET_ee_lowpt");
	gmetvsht_da_mm_lowpt->SetMarkerColor(col_mm);
	gmetvsht_da_mm_lowpt->SetMarkerStyle(8);
	gmetvsht_da_mm_lowpt->SetMarkerSize(1.5);
	gmetvsht_da_em_lowpt->SetMarkerColor(col_em);
	gmetvsht_da_em_lowpt->SetMarkerStyle(23);
	gmetvsht_da_em_lowpt->SetMarkerSize(1.5);
	gmetvsht_da_ee_lowpt->SetMarkerColor(col_ee);
	gmetvsht_da_ee_lowpt->SetMarkerStyle(21);
	gmetvsht_da_ee_lowpt->SetMarkerSize(1.5);

///////////// FROM FILE

	// TGraphs:
	TMultiGraph *gmetvsht_da_mm = new TMultiGraph("HTvsMET_mm", "HTvsMET_mm");
	TMultiGraph *gmetvsht_da_ee = new TMultiGraph("HTvsMET_ee", "HTvsMET_ee");
	TMultiGraph *gmetvsht_da_em = new TMultiGraph("HTvsMET_em", "HTvsMET_em");

	vector<int> mmsamples = fMuData;
	vector<int> eesamples = fEGData;
	vector<int> emsamples = fMuEGData;
	for(size_t i = 0; i < mmsamples.size(); ++i) gmetvsht_da_mm->Add(fSamples[mmsamples[i]]->sigevents[Muon][HighPt]);
	for(size_t i = 0; i < eesamples.size(); ++i) gmetvsht_da_ee->Add(fSamples[eesamples[i]]->sigevents[Elec][HighPt]);
	for(size_t i = 0; i < emsamples.size(); ++i) gmetvsht_da_em->Add(fSamples[emsamples[i]]->sigevents[ElMu][HighPt]);

///////////////////////

	hmetvsht_da_mm->SetXTitle(KinPlots::axis_label[0]);
	hmetvsht_da_mm->SetYTitle(KinPlots::axis_label[1]);
	hmetvsht_da_ee->SetXTitle(KinPlots::axis_label[0]);
	hmetvsht_da_ee->SetYTitle(KinPlots::axis_label[1]);
	hmetvsht_da_em->SetXTitle(KinPlots::axis_label[0]);
	hmetvsht_da_em->SetYTitle(KinPlots::axis_label[1]);

	// TLegend *leg = new TLegend(0.80,0.70,0.95,0.88);
	TLegend *leg = new TLegend(0.67,0.70,0.82,0.88);
	leg->AddEntry(hmetvsht_da_mm, "#mu#mu","p");
	leg->AddEntry(hmetvsht_da_ee, "ee","p");
	leg->AddEntry(hmetvsht_da_em, "e#mu","p");
	leg->SetFillStyle(0);
	leg->SetTextFont(42);
	leg->SetTextSize(0.05);
	leg->SetBorderSize(0);

	//////////////////////////////////////////////////////////
	// TAUS //////////////////////////////////////////////////
	// Create histograms
	TH2D *hmetvsht_da_mt = new TH2D("Data_HTvsMET_mt", "Data_HTvsMET_mt", 100, 0., htmax, 100, 0., metmax);
	TH2D *hmetvsht_da_et = new TH2D("Data_HTvsMET_et", "Data_HTvsMET_et", 100, 0., htmax, 100, 0., metmax);
	TH2D *hmetvsht_da_tt = new TH2D("Data_HTvsMET_tt", "Data_HTvsMET_tt", 100, 0., htmax, 100, 0., metmax);

	//////////////////////////////////////////////////////////
	// Make MET vs HT plot:
	hmetvsht_da_mt->SetMarkerStyle(22);
	hmetvsht_da_mt->SetMarkerSize(1.8);
	hmetvsht_da_mt->SetMarkerColor(col_mt);
	hmetvsht_da_mt->GetYaxis()->SetTitleOffset(1.4);

	hmetvsht_da_et->SetMarkerStyle(33);
	hmetvsht_da_et->SetMarkerSize(2.3);
	hmetvsht_da_et->SetMarkerColor(col_et);
	// hmetvsht_da_et->SetMarkerColor(46);

	hmetvsht_da_tt->SetMarkerStyle(34);
	hmetvsht_da_tt->SetMarkerSize(1.8);
	hmetvsht_da_tt->SetMarkerColor(col_tt);
	// hmetvsht_da_tt->SetMarkerColor(38);


	// Updated numbers with 3.2/fb from Nov 9
	const int nmtev = 22;
	float a_mt_ht [nmtev] = {454.149, 380.437, 363.305, 549.702, 361.937, 355.314, 363.292, 440.719, 640.22, 487.826, 481.437, 357.928, 464.144, 532.001, 1022.7, 378.033, 382.573, 368.67,  428.258, 864.134, 430.531, 398.942};
	float a_mt_met[nmtev] = {144.284,  91.627, 217.503, 143.927, 81.9502, 109.763,  89.345, 111.706, 82.201, 161.936, 120.827, 113.490,  85.155, 107.131,  82.53, 118.257, 107.807, 133.973, 108.852, 134.839, 107.308, 142.304};
	const int netev = 12;
	float a_et_ht [netev] = {421.656, 421.634, 537.589, 444.368, 388.334, 393.389, 418.707, 399.363, 393.997, 555.166, 612.39, 735.811};
	float a_et_met[netev] = {103.668, 93.5886, 83.3471, 194.835, 80.4797, 88.2911, 117.529, 165.110, 83.4841,  88.946,  88.53, 183.931};

	TGraph *gmetvsht_da_mt = new TGraph(nmtev, a_mt_ht, a_mt_met);
	gmetvsht_da_mt->SetName("Data_HTvsMET_mt_graph");
	gmetvsht_da_mt->SetMarkerStyle(hmetvsht_da_mt->GetMarkerStyle());
	gmetvsht_da_mt->SetMarkerSize( hmetvsht_da_mt->GetMarkerSize());
	gmetvsht_da_mt->SetMarkerColor(hmetvsht_da_mt->GetMarkerColor());

	TGraph *gmetvsht_da_et = new TGraph(netev, a_et_ht, a_et_met);
	gmetvsht_da_et->SetName("Data_HTvsMET_et_graph");
	gmetvsht_da_et->SetMarkerStyle(hmetvsht_da_et->GetMarkerStyle());
	gmetvsht_da_et->SetMarkerSize( hmetvsht_da_et->GetMarkerSize());
	gmetvsht_da_et->SetMarkerColor(hmetvsht_da_et->GetMarkerColor());

	hmetvsht_da_mt->SetXTitle(KinPlots::axis_label[0]);
	hmetvsht_da_mt->SetYTitle(KinPlots::axis_label[1]);

	// TLegend *leg = new TLegend(0.80,0.82,0.95,0.88);
	TLegend *leg2 = new TLegend(0.80,0.70,0.95,0.88);
	leg2->AddEntry(hmetvsht_da_mt, "#mu#tau","p");
	leg2->AddEntry(hmetvsht_da_et, "e#tau","p");
	leg2->AddEntry(hmetvsht_da_tt, "#tau#tau","p");
	leg2->SetFillStyle(0);
	leg2->SetTextFont(42);
	leg2->SetTextSize(0.05);
	leg2->SetBorderSize(0);

	//////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////

	// Special effects:
	const float lowerht = hilo==HighPt? 80.:200.;

	TWbox *lowhtbox  = new TWbox(0., 0., lowerht, metmax, kBlack, 0, 0);
	TWbox *lowmetbox = new TWbox(lowerht, 0., htmax,    30., kBlack, 0, 0);
	lowhtbox ->SetFillColor(12);
	lowmetbox->SetFillColor(12);
	lowhtbox ->SetFillStyle(3005);
	lowmetbox->SetFillStyle(3005);
	TLine *boxborder1 = new TLine(lowerht,30.,lowerht, metmax);
	TLine *boxborder2 = new TLine(lowerht,30.,htmax,30.);
	boxborder1->SetLineWidth(1);
	boxborder2->SetLineWidth(1);
	boxborder1->SetLineColor(14);
	boxborder2->SetLineColor(14);

	TLine *sig1x = new TLine(lowerht, 120., lowerht, metmax); // met 120 ht 80
	TLine *sig1y = new TLine(lowerht, 120., htmax,   120.); 
	TLine *sig2x = new TLine(200., 120., 200.,  metmax);      // met 120 ht 200
	TLine *sig2y = new TLine(200., 120., htmax, 120.);
	TLine *sig3x = new TLine(450.,  50., 450.,  metmax);      // met 50  ht 450
	TLine *sig3y = new TLine(450.,  50., htmax,  50.);
	TLine *sig4x = new TLine(450., 120., 450.,  metmax);      // met 120 ht 450
	TLine *sig4y = new TLine(450., 120., htmax, 120.);

	sig1x->SetLineWidth(2);
	sig1y->SetLineWidth(2);
	sig2x->SetLineWidth(2);
	sig2y->SetLineWidth(2);
	sig3x->SetLineWidth(1);
	sig3y->SetLineWidth(1);
	sig4x->SetLineWidth(3);
	sig4y->SetLineWidth(3);

	sig1x->SetLineStyle(3);
	sig1y->SetLineStyle(3);
	sig2x->SetLineStyle(2);
	sig2y->SetLineStyle(2);
	sig3x->SetLineStyle(1);
	sig3y->SetLineStyle(1);
	// sig4x->SetLineStyle(1);
	// sig4y->SetLineStyle(1);

	// TLegend *regleg = new TLegend(0.70,0.47,0.88,0.6);
	// TLegend *regleg = new TLegend(0.67,0.51,0.87,0.68);
	TLegend *regleg = new TLegend(0.76,0.29,0.93,0.46);
	regleg->AddEntry(sig4x, "Region 1","l");
	regleg->AddEntry(sig2x, "Region 2","l");
	regleg->AddEntry(sig3x, "Region 3","l");
	regleg->AddEntry(sig1x, "Region 4","l");
	regleg->SetFillStyle(0);
	regleg->SetTextFont(42);
	regleg->SetTextSize(0.03);
	regleg->SetBorderSize(0);
	

	TCanvas *c_temp = new TCanvas("C_HTvsMET", "HT vs MET in Data vs MC", 0, 0, 600, 600);
	c_temp->cd();
	c_temp->SetRightMargin(0.05);
	c_temp->SetLeftMargin(0.13);

	hmetvsht_da_mm->DrawCopy("axis");

	lowhtbox ->Draw();
	lowmetbox->Draw();
	boxborder1->Draw();
	boxborder2->Draw();

	if(hilo != LowPt) sig1x->Draw();
	if(hilo != LowPt) sig1y->Draw();
	sig2x->Draw();
	sig2y->Draw();
	sig3x->Draw();
	sig3y->Draw();
	sig4x->Draw();
	sig4y->Draw();

	// Graphs
	gmetvsht_da_ee->Draw("P");
	gmetvsht_da_em->Draw("P");
	gmetvsht_da_mm->Draw("P");
	
	gmetvsht_da_ee_lowpt->Draw("P");
	gmetvsht_da_em_lowpt->Draw("P");
	gmetvsht_da_mm_lowpt->Draw("P");
	
	gmetvsht_da_mt->Draw("P");
	gmetvsht_da_et->Draw("P");
	
	leg->Draw();
	leg2->Draw();
	regleg->Draw();

	drawTopLine();
	// TPaveText *pave = new TPaveText(0.16, 0.83, 0.55, 0.88, "NDC");
	// pave->SetFillColor(0);
	// pave->SetFillStyle(1001);
	// pave->SetBorderSize(0);
	// pave->SetMargin(0.05);
	// pave->SetTextFont(42);
	// pave->SetTextSize(0.04);
	// pave->SetTextAlign(12);
	// pave->AddText(hiloname[hilo]);
	// pave->Draw();
	gPad->RedrawAxis();

	// Util::PrintNoEPS(c_temp, "HTvsMET_" + gHiLoLabel[hilo], fOutputDir + fOutputSubDir, NULL);
	Util::PrintPDF(c_temp, "HTvsMET_PRL", fOutputDir + fOutputSubDir);
	// Util::SaveAsMacro(c_temp, "HTvsMET_" + gHiLoLabel[hilo], fOutputDir + fOutputSubDir);
	delete c_temp;
	delete leg, regleg;
	delete hmetvsht_da_mm, hmetvsht_da_ee, hmetvsht_da_em;//, hmetvsht_mc;
	delete gmetvsht_da_mm, gmetvsht_da_ee, gmetvsht_da_em;//, hmetvsht_mc;
}
void SSDLPlotter::makeMETvsHTPlotTau(){
	fOutputSubDir = "KinematicPlots/HTvsMET/";
	char cmd[100];
    sprintf(cmd,"mkdir -p %s%s", fOutputDir.Data(), fOutputSubDir.Data());
    system(cmd);
	
	const float htmax = 1100.;
	const float metmax = 300.;

	// Create histograms
	TH2D *hmetvsht_da_mt = new TH2D("Data_HTvsMET_mt", "Data_HTvsMET_mt", 100, 0., htmax, 100, 0., metmax);
	TH2D *hmetvsht_da_et = new TH2D("Data_HTvsMET_et", "Data_HTvsMET_et", 100, 0., htmax, 100, 0., metmax);
	TH2D *hmetvsht_da_tt = new TH2D("Data_HTvsMET_tt", "Data_HTvsMET_tt", 100, 0., htmax, 100, 0., metmax);

	//////////////////////////////////////////////////////////
	// Make MET vs HT plot:
	hmetvsht_da_mt->SetMarkerStyle(22);
	hmetvsht_da_mt->SetMarkerSize(1.8);
	hmetvsht_da_mt->SetMarkerColor(51);
	hmetvsht_da_mt->GetYaxis()->SetTitleOffset(1.4);

	hmetvsht_da_et->SetMarkerStyle(33);
	hmetvsht_da_et->SetMarkerSize(2.3);
	hmetvsht_da_et->SetMarkerColor(46);

	hmetvsht_da_tt->SetMarkerStyle(34);
	hmetvsht_da_tt->SetMarkerSize(1.8);
	hmetvsht_da_tt->SetMarkerColor(38);


	// Updated numbers with 3.2/fb from Nov 9
	const int nmtev = 16;
	float a_mt_ht [nmtev] = {454.149, 380.437, 363.305, 549.702, 361.937, 355.314, 363.292, 440.719, 640.22, 487.826, 481.437, 357.928, 464.144, 532.001, 1022.7, 378.033};
	float a_mt_met[nmtev] = {144.284,  91.627, 217.503, 143.927, 81.9502, 109.763,  89.345, 111.706, 82.201, 161.936, 120.827, 113.490,  85.155, 107.131,  82.53, 118.257};
	const int netev = 10;
	float a_et_ht [netev] = {421.656, 421.634, 537.589, 444.368, 388.334, 393.389, 418.707, 399.363, 393.997, 555.166};
	float a_et_met[netev] = {103.668, 93.5886, 83.3471, 194.835, 80.4797, 88.2911, 117.529, 165.110, 83.4841,  88.946};

	TGraph *gmetvsht_da_mt = new TGraph(nmtev, a_mt_ht, a_mt_met);
	gmetvsht_da_mt->SetName("Data_HTvsMET_mt_graph");
	gmetvsht_da_mt->SetMarkerStyle(hmetvsht_da_mt->GetMarkerStyle());
	gmetvsht_da_mt->SetMarkerSize( hmetvsht_da_mt->GetMarkerSize());
	gmetvsht_da_mt->SetMarkerColor(hmetvsht_da_mt->GetMarkerColor());

	TGraph *gmetvsht_da_et = new TGraph(netev, a_et_ht, a_et_met);
	gmetvsht_da_et->SetName("Data_HTvsMET_et_graph");
	gmetvsht_da_et->SetMarkerStyle(hmetvsht_da_et->GetMarkerStyle());
	gmetvsht_da_et->SetMarkerSize( hmetvsht_da_et->GetMarkerSize());
	gmetvsht_da_et->SetMarkerColor(hmetvsht_da_et->GetMarkerColor());

	hmetvsht_da_mt->SetXTitle(KinPlots::axis_label[0]);
	hmetvsht_da_mt->SetYTitle(KinPlots::axis_label[1]);

	// TLegend *leg = new TLegend(0.80,0.82,0.95,0.88);
	TLegend *leg = new TLegend(0.80,0.70,0.95,0.88);
	leg->AddEntry(hmetvsht_da_mt, "#mu#tau","p");
	leg->AddEntry(hmetvsht_da_et, "e#tau","p");
	leg->AddEntry(hmetvsht_da_tt, "#tau#tau","p");
	leg->SetFillStyle(0);
	leg->SetTextFont(42);
	leg->SetTextSize(0.05);
	leg->SetBorderSize(0);

	// Special effects:
	const float lowerht = 350.;
	// HT350met80 and HT400met120).
	TWbox *lowhtbox  = new TWbox(0., 0., lowerht, metmax, kBlack, 0, 0);
	TWbox *lowmetbox = new TWbox(lowerht, 0., htmax,    80., kBlack, 0, 0);
	lowhtbox ->SetFillColor(12);
	lowmetbox->SetFillColor(12);
	lowhtbox ->SetFillStyle(3005);
	lowmetbox->SetFillStyle(3005);
	TLine *boxborder1 = new TLine(lowerht, 80., lowerht, metmax);
	TLine *boxborder2 = new TLine(lowerht, 80., htmax,   80.);
	boxborder1->SetLineWidth(1);
	boxborder2->SetLineWidth(1);
	boxborder1->SetLineColor(14);
	boxborder2->SetLineColor(14);

	// TLine *sig1x = new TLine(lowerht, 100., lowerht, 400.);
	// TLine *sig1y = new TLine(lowerht, 100., htmax,   100.);
	// TLine *sig2x = new TLine(lowerht, 120., lowerht, 400.);
	// TLine *sig2y = new TLine(lowerht, 120., htmax,   120.);
	// TLine *sig3x = new TLine(400.,  50., 400.,  400.);
	// TLine *sig3y = new TLine(400.,  50., htmax,  50.);
	TLine *sig4x = new TLine(450., 120., 450.,  metmax);
	TLine *sig4y = new TLine(450., 120., htmax, 120.);

	// sig1x->SetLineWidth(2);
	// sig1y->SetLineWidth(2);
	// sig2x->SetLineWidth(2);
	// sig2y->SetLineWidth(2);
	// sig3x->SetLineWidth(2);
	// sig3y->SetLineWidth(2);
	sig4x->SetLineWidth(3);
	sig4y->SetLineWidth(3);

	// sig1x->SetLineStyle(2);
	// sig1y->SetLineStyle(2);
	// sig2x->SetLineStyle(2);
	// sig2y->SetLineStyle(2);
	// sig3x->SetLineStyle(2);
	// sig3y->SetLineStyle(2);
	// sig4x->SetLineStyle(2);
	// sig4y->SetLineStyle(2);

	TLegend *regleg = new TLegend(0.70,0.60,0.88,0.65);
	regleg->AddEntry(sig4x, "Search Region 1","l");
	regleg->SetFillStyle(0);
	regleg->SetTextFont(42);
	regleg->SetTextSize(0.03);
	regleg->SetBorderSize(0);

	TCanvas *c_temp = new TCanvas("C_HTvsMET", "HT vs MET in Data vs MC", 0, 0, 600, 600);
	c_temp->cd();
	c_temp->SetRightMargin(0.03);
	c_temp->SetLeftMargin(0.13);

	hmetvsht_da_mt->Draw("axis");

	lowhtbox ->Draw();
	lowmetbox->Draw();
	boxborder1->Draw();
	boxborder2->Draw();

	// if(hilo != LowPt) sig1x->Draw();
	// if(hilo != LowPt) sig1y->Draw();
	// sig2x->Draw();
	// sig2y->Draw();
	// sig3x->Draw();
	// sig3y->Draw();
	sig4x->Draw();
	sig4y->Draw();

	// Graphs
	gmetvsht_da_mt->Draw("P");
	gmetvsht_da_et->Draw("P");
	
	leg->Draw();
	regleg->Draw();
	drawTopLine();
	TPaveText *pave = new TPaveText(0.16, 0.83, 0.55, 0.88, "NDC");
	pave->SetFillColor(0);
	pave->SetFillStyle(1001);
	pave->SetBorderSize(0);
	pave->SetMargin(0.05);
	pave->SetTextFont(42);
	pave->SetTextSize(0.04);
	pave->SetTextAlign(12);
	pave->AddText("p_{T}(#mu/e/#tau) > 5/10/15 GeV");
	pave->Draw();
	gPad->RedrawAxis();

	// Util::PrintNoEPS( c_temp, "HTvsMET_TauChan", fOutputDir + fOutputSubDir, NULL);
	Util::PrintPDF(   c_temp, "HTvsMET_TauChan", fOutputDir + fOutputSubDir);
	// Util::SaveAsMacro(c_temp, "HTvsMET_TauChan", fOutputDir + fOutputSubDir);
	delete c_temp;
	delete leg, regleg;
	delete hmetvsht_da_mt;
	delete gmetvsht_da_mt;
	fOutputSubDir = "";
}
void SSDLPlotter::makeFRvsPtPlots(gChannel chan, gFPSwitch fp){
	fOutputSubDir = "Ratios/";
	char cmd[100];
    sprintf(cmd,"mkdir -p %s%s", fOutputDir.Data(), fOutputSubDir.Data());
    system(cmd);

	TString pfname = "Non-prompt ";
	if(fp == ZDecay) pfname = "Prompt ";
	
	TString name;
	if(chan == Muon) name = "Muons";
	if(chan == Elec) name = "Electrons";

	TH1D *h_dummy1 = new TH1D("dummy1", "dummy1", getNEtaBins(chan), getEtaBins(chan));
	TH2D *h_dummy2 = new TH2D("dummy2", "dummy2", getNPt2Bins(chan), getPt2Bins(chan), getNEtaBins(chan), getEtaBins(chan));

	TH1D *h_ptratio_data = new TH1D("Ratio_data", "Tight/Loose Ratio in data", getNPt2Bins(chan), getPt2Bins(chan));
	TH1D *h_ptratio_mc   = new TH1D("Ratio_mc",   "Tight/Loose Ratio in MC",   getNPt2Bins(chan), getPt2Bins(chan));

	vector<int> datasamples;
	vector<int> mcsamples;

	if(chan == Muon){
		datasamples = fMuData;
		mcsamples   = fMCBGMuEnr;
	}
	if(chan == Elec){
		datasamples = fEGData;
		mcsamples   = fMCBG;
	}

	calculateRatio(datasamples, chan, fp, h_dummy2, h_ptratio_data, h_dummy1);
	calculateRatio(mcsamples,   chan, fp, h_dummy2, h_ptratio_mc,   h_dummy1);

	float maximum = 0.8;
	if(fp == ZDecay) maximum = 1.1;
	h_ptratio_data->SetMaximum(maximum);
	h_ptratio_mc  ->SetMaximum(maximum);
	h_ptratio_data->SetMinimum(0.0);
	h_ptratio_mc  ->SetMinimum(0.0);

	if(chan == Muon) h_ptratio_mc->SetXTitle(convertVarName("MuPt[0]"));
	if(chan == Elec) h_ptratio_mc->SetXTitle(convertVarName("ElPt[0]"));
	h_ptratio_mc->GetYaxis()->SetTitleOffset(1.2);
	h_ptratio_mc->SetYTitle("N_{Tight}/N_{Loose}");

	h_ptratio_data->SetMarkerColor(kBlack);
	h_ptratio_data->SetMarkerStyle(20);
	h_ptratio_data->SetMarkerSize(1.5);
	h_ptratio_data->SetLineWidth(2);
	h_ptratio_data->SetLineColor(kBlack);
	h_ptratio_data->SetFillColor(kBlack);
	
	h_ptratio_mc  ->SetMarkerColor(kRed);
	h_ptratio_mc  ->SetMarkerStyle(23);
	h_ptratio_mc  ->SetMarkerSize(1.5);
	h_ptratio_mc  ->SetLineWidth(2);
	h_ptratio_mc  ->SetLineColor(kRed);
	h_ptratio_mc  ->SetFillColor(kRed);

	TLatex *lat = new TLatex();
	lat->SetNDC(kTRUE);
	lat->SetTextColor(kBlack);
	lat->SetTextSize(0.04);
	
	TLegend *leg;
	if(fp == SigSup) leg = new TLegend(0.15,0.75,0.35,0.88);
	if(fp == ZDecay) leg = new TLegend(0.15,0.15,0.35,0.28);
	leg->AddEntry(h_ptratio_data, "Data",       "p");
	leg->AddEntry(h_ptratio_mc,   "Simulation", "p");
	leg->SetTextSize(0.04);
	leg->SetFillStyle(0);
	leg->SetTextFont(42);
	leg->SetBorderSize(0);

	TCanvas *c_temp = new TCanvas("C_PtRatioPlot", "fRatio vs Pt in Data vs MC", 0, 0, 800, 600);
	c_temp->cd();

	h_ptratio_mc->DrawCopy("PE X0");
	h_ptratio_data->Draw("PE X0 same");
	leg->Draw();
	// lat->DrawLatex(0.70,0.92, Form("L_{int.} = %2.1f fb^{-1}", fLumiNorm/1000.));
	lat->SetTextSize(0.04);
	if(fp == SigSup) lat->DrawLatex(0.62,0.85, pfname + name);
	if(fp == ZDecay) lat->DrawLatex(0.67,0.15, pfname + name);
	double ymean(0.), yrms(0.);
	getWeightedYMeanRMS(h_ptratio_data, ymean, yrms);
	drawTopLine();
	lat->SetTextSize(0.03);
	// lat->DrawLatex(0.25,0.92, Form("Mean ratio: %4.2f #pm %4.2f", ymean, yrms));

	TString fpname = "F";
	if(fp == ZDecay) fpname = "P";
	
	// Util::PrintNoEPS( c_temp, fpname + "Ratio_" + name + "_Pt", fOutputDir + fOutputSubDir, NULL);
	Util::PrintPDF(c_temp, fpname + "Ratio_" + name + "_Pt", fOutputDir + fOutputSubDir);
	delete h_ptratio_mc, h_ptratio_data;
	delete c_temp, lat, leg;
	fOutputSubDir = "";
}
void SSDLPlotter::makeFRvsPtPlotsForPAS(gChannel chan){
	Util::SetTDRStyle();
	fOutputSubDir = "Ratios/forPAS/";
	char cmd[100];
    sprintf(cmd,"mkdir -p %s%s", fOutputDir.Data(), fOutputSubDir.Data());
    system(cmd);

	TString name;
	if(chan == Muon)     name = "Muons";
	if(chan == Elec) name = "Electrons";

	TH1D *h_dummy1 = new TH1D("dummy1", "dummy1", getNEtaBins(chan), getEtaBins(chan));
	TH2D *h_dummy2 = new TH2D("dummy2", "dummy2", getNPt2Bins(chan), getPt2Bins(chan), getNEtaBins(chan), getEtaBins(chan));

	TH1D *h_ratio_A2 = new TH1D("Ratio_data2", "Tight/Loose Ratio in data for A2", getNPt2Bins(chan), getPt2Bins(chan));
	TH1D *h_ratio_A1 = new TH1D("Ratio_data1", "Tight/Loose Ratio in data for A1", getNPt2Bins(chan), getPt2Bins(chan));

	vector<int> datasamples;

	if(chan == Muon)     datasamples = fMuData;
	if(chan == Elec) datasamples = fEGData;

	float A1MBins [5] = {0.31, 0.25, 0.21, 0.18, 0.20};
	float A1MBinsE[5] = {0.002, 0.002, 0.002, 0.003, 0.001};

	float A1EBins [5] = {0.22, 0.21, 0.14, 0.22, 0.35};
	float A1EBinsE[5] = {0.02, 0.02, 0.02, 0.02, 0.04};

	for(size_t i = 0; i < 5; ++i){
		if(chan == Muon){
			h_ratio_A1->SetBinContent(i+1, A1MBins[i]);
			h_ratio_A1->SetBinError(i+1,   A1MBinsE[i]);
		}
		if(chan == Elec){
			h_ratio_A1->SetBinContent(i+1, A1EBins[i]);
			h_ratio_A1->SetBinError(i+1,   A1EBinsE[i]);
		}
	}

	calculateRatio(datasamples, chan, SigSup, h_dummy2, h_ratio_A2, h_dummy1);

	h_ratio_A2->GetXaxis()->SetTitle("p_{T} (GeV)");
	h_ratio_A2->GetYaxis()->SetTitle("TL Ratio");
	h_ratio_A2->SetMarkerStyle(20);
	h_ratio_A2->SetMarkerSize(1.6);
	h_ratio_A2->SetMarkerColor(kBlue);
	h_ratio_A2->SetLineColor(kBlue);

	h_ratio_A1->GetXaxis()->SetTitle("p_{T} (GeV)");
	h_ratio_A1->GetYaxis()->SetTitle("TL Ratio");
	h_ratio_A1->SetMarkerStyle(23);
	h_ratio_A1->SetMarkerSize(1.8);
	h_ratio_A1->SetMarkerColor(kBlack);
	h_ratio_A1->SetLineColor(kBlack);

	h_ratio_A2->GetYaxis()->SetRangeUser(0., 0.7);
	h_ratio_A1->GetYaxis()->SetRangeUser(0., 0.7);

	TLegend *leg = new TLegend(0.21,0.58,0.47,0.78);
	leg->AddEntry(h_ratio_A1, "Method A1","p");
	leg->AddEntry(h_ratio_A2, "Method A2","p");
	leg->SetTextSize(0.05);
	// leg->SetTextFont(42);
	leg->SetFillStyle(0);
	leg->SetBorderSize(0);	

	TCanvas *c_temp = new TCanvas();
	c_temp->cd();
	c_temp->SetRightMargin(0.05);

	h_ratio_A2->Draw("PE");
	h_ratio_A1->Draw("PE same");
	leg->Draw();
	TLatex lat;
	lat.SetNDC(kTRUE);
	lat.SetTextSize(0.05);
	lat.DrawLatex(0.23, 0.88, "CMS Preliminary");
	lat.DrawLatex(0.70, 0.88, name);
	lat.DrawLatex(0.23, 0.81, "L_{int.} = 0.98 fb^{-1},   #sqrt{s} = 7 TeV");

	Util::PrintNoEPS( c_temp, "FRatio_" + name + "_Pt_A1vsA2", fOutputDir + fOutputSubDir, NULL);
	Util::PrintPDF(   c_temp, "FRatio_" + name + "_Pt_A1vsA2", fOutputDir + fOutputSubDir);
	Util::SaveAsMacro(c_temp, "FRatio_" + name + "_Pt_A1vsA2", fOutputDir + fOutputSubDir);
	delete h_ratio_A1, h_ratio_A2;
	delete c_temp;
	fOutputSubDir = "";
}
void SSDLPlotter::makeFRvsEtaPlots(gChannel chan){
	fOutputSubDir = "Ratios/";
	char cmd[100];
    sprintf(cmd,"mkdir -p %s%s", fOutputDir.Data(), fOutputSubDir.Data());
    system(cmd);

	TString name;
	if(chan == Muon)     name = "Muons";
	if(chan == Elec) name = "Electrons";

	TString pfname = "Non-prompt ";

	TH1D *h_dummy1 = new TH1D("dummy1", "dummy1", getNPt2Bins(chan), getPt2Bins(chan));
	TH2D *h_dummy2 = new TH2D("dummy2", "dummy2", getNPt2Bins(chan), getPt2Bins(chan), getNEtaBins(chan), getEtaBins(chan));

	TH1D *h_etaratio_data = new TH1D("Ratio_data", "Tight/Loose Ratio in data", getNEtaBins(chan), getEtaBins(chan));
	TH1D *h_etaratio_mc   = new TH1D("Ratio_mc",   "Tight/Loose Ratio in MC",   getNEtaBins(chan), getEtaBins(chan));
	// h_etaratio_data->Sumw2();
	// h_etaratio_mc  ->Sumw2();

	vector<int> datasamples;
	vector<int> mcsamples;

	if(chan == Muon){
		datasamples = fMuData;
		mcsamples   = fMCBGMuEnr;
	}
	if(chan == Elec){
		datasamples = fEGData;
		mcsamples   = fMCBG;
	}

	calculateRatio(datasamples, chan, SigSup, h_dummy2, h_dummy1, h_etaratio_data);
	calculateRatio(mcsamples,   chan, SigSup, h_dummy2, h_dummy1, h_etaratio_mc);

	float max = 0.4;
	if(chan==Elec) max = 0.8;
	h_etaratio_data->SetMaximum(max);
	h_etaratio_mc  ->SetMaximum(max);
	h_etaratio_data->SetMinimum(0.0);
	h_etaratio_mc  ->SetMinimum(0.0);

	if(chan == Muon) h_etaratio_mc->SetXTitle(convertVarName("MuEta[0]"));
	if(chan == Elec) h_etaratio_mc->SetXTitle(convertVarName("ElEta[0]"));
	h_etaratio_mc->GetYaxis()->SetTitleOffset(1.2);
	h_etaratio_mc->SetYTitle("N_{Tight}/N_{Loose}");
	
	h_etaratio_data->SetMarkerColor(kBlack);
	h_etaratio_data->SetMarkerStyle(20);
	h_etaratio_data->SetMarkerSize(1.5);
	h_etaratio_data->SetLineWidth(2);
	h_etaratio_data->SetLineColor(kBlack);
	h_etaratio_data->SetFillColor(kBlack);
	
	h_etaratio_mc  ->SetMarkerColor(kRed);
	h_etaratio_mc  ->SetMarkerStyle(23);
	h_etaratio_mc  ->SetMarkerSize(1.5);
	h_etaratio_mc  ->SetLineWidth(2);
	h_etaratio_mc  ->SetLineColor(kRed);
	h_etaratio_mc  ->SetFillColor(kRed);

	// // h_etaratio_data->GetXaxis()->SetTitle("p_{T} (GeV)");
	// h_etaratio_data->GetXaxis()->SetTitle("#left|#eta#right|");
	// h_etaratio_data->GetYaxis()->SetTitle("TL Ratio");
	// h_etaratio_data->SetMarkerStyle(22);
	// h_etaratio_data->SetMarkerSize(1.4);
	// h_etaratio_data->SetMarkerColor(kBlack);
	// h_etaratio_data->SetLineColor(kBlack);
	// h_etaratio_data->GetYaxis()->SetRangeUser(0., 0.5);

	TLatex *lat = new TLatex();
	lat->SetNDC(kTRUE);
	lat->SetTextColor(kBlack);
	lat->SetTextSize(0.04);
	
	TLegend *leg = new TLegend(0.15,0.75,0.35,0.88);
	leg->AddEntry(h_etaratio_data, "Data",       "p");
	leg->AddEntry(h_etaratio_mc,   "Simulation", "p");
	leg->SetFillStyle(0);
	leg->SetTextFont(42);
	leg->SetTextSize(0.04);
	leg->SetBorderSize(0);

	TCanvas *c_temp = new TCanvas("C_EtaRatioPlot", "fRatio vs Eta in Data vs MC", 0, 0, 800, 600);
	// TCanvas *c_temp = new TCanvas();
	c_temp->cd();

	// gPad->SetLogy();
	// h_etaratio_data->Draw("PE");
	// TLatex lat;
	// lat.SetNDC(kTRUE);
	//     lat.SetTextSize(0.028);
	//     lat.DrawLatex(0.23, 0.88, "CMS Preliminary");
	//     lat.DrawLatex(0.23, 0.79, "#int L dt = XXX pb^{-1},   #sqrt{s} = 7 TeV");
	//     lat.DrawLatex(0.83, 0.88, name);
	
	h_etaratio_mc->DrawCopy("PE X0");
	h_etaratio_data->Draw("PE X0 same");
	leg->Draw();
	drawTopLine();
	lat->SetTextSize(0.04);
	lat->DrawLatex(0.62,0.85, pfname + name);
	double ymean(0.), yrms(0.);
	getWeightedYMeanRMS(h_etaratio_data, ymean, yrms);
	lat->SetTextSize(0.03);
	// lat->DrawLatex(0.25,0.92, Form("Mean ratio: %4.2f #pm %4.2f", ymean, yrms));
	
	// Util::PrintNoEPS( c_temp, "FRatio_" + name + "_Eta", fOutputDir + fOutputSubDir, NULL);
	Util::PrintPDF(   c_temp, "FRatio_" + name + "_Eta", fOutputDir + fOutputSubDir);
	delete h_etaratio_mc, h_etaratio_data;
	// delete c_temp;
	delete c_temp, lat, leg;
	fOutputSubDir = "";
}
void SSDLPlotter::makeFRvsEtaPlotsForPAS(gChannel chan){
	Util::SetTDRStyle();
	fOutputSubDir = "Ratios/forPAS";
	char cmd[100];
    sprintf(cmd,"mkdir -p %s%s", fOutputDir.Data(), fOutputSubDir.Data());
    system(cmd);

	TString name;
	if(chan == Muon)     name = "Muons";
	if(chan == Elec) name = "Electrons";

	TH1D *h_dummy1 = new TH1D("dummy1", "dummy1", getNPt2Bins(chan), getPt2Bins(chan));
	TH2D *h_dummy2 = new TH2D("dummy2", "dummy2", getNPt2Bins(chan), getPt2Bins(chan), getNEtaBins(chan), getEtaBins(chan));

	TH1D *h_ratio_A2 = new TH1D("Ratio_data2", "Tight/Loose Ratio in data for A2", getNEtaBins(chan), getEtaBins(chan));
	TH1D *h_ratio_A1 = new TH1D("Ratio_data1", "Tight/Loose Ratio in data for A1", getNEtaBins(chan), getEtaBins(chan));

	vector<int> datasamples;

	if(chan == Muon)     datasamples = fMuData;
	if(chan == Elec) datasamples = fEGData;

	calculateRatio(datasamples, chan, SigSup, h_dummy2, h_dummy1, h_ratio_A2);

	float A1MBins [4] = {0.2, 0.23, 0.25, 0.26};
	float A1MBinsE[4] = {0.001, 0.002, 0.002, 0.002};

	float A1EBins [4] = {0.2, 0.22, 0.23, 0.28};
	float A1EBinsE[4] = {0.01, 0.02, 0.02, 0.03};

	for(size_t i = 0; i < 4; ++i){
		if(chan == Muon){
			h_ratio_A1->SetBinContent(i+1, A1MBins[i]);
			h_ratio_A1->SetBinError(i+1,   A1MBinsE[i]);
		}
		if(chan == Elec){
			h_ratio_A1->SetBinContent(i+1, A1EBins[i]);
			h_ratio_A1->SetBinError(i+1,   A1EBinsE[i]);
		}
	}

	h_ratio_A2->SetMaximum(0.4);
	h_ratio_A2->SetMinimum(0.0);
	h_ratio_A1->SetMaximum(0.4);
	h_ratio_A1->SetMinimum(0.0);

	// h_ratio_data->GetXaxis()->SetTitle("p_{T} (GeV)");
	h_ratio_A2->GetXaxis()->SetTitle("#left|#eta#right|");
	h_ratio_A2->GetYaxis()->SetTitle("TL Ratio");
	h_ratio_A2->SetMarkerStyle(20);
	h_ratio_A2->SetMarkerSize(1.6);
	h_ratio_A2->SetMarkerColor(kBlue);
	h_ratio_A2->SetLineColor(kBlue);

	h_ratio_A1->GetXaxis()->SetTitle("#left|#eta#right|");
	h_ratio_A1->GetYaxis()->SetTitle("TL Ratio");
	h_ratio_A1->SetMarkerStyle(23);
	h_ratio_A1->SetMarkerSize(1.8);
	h_ratio_A1->SetMarkerColor(kBlack);
	h_ratio_A1->SetLineColor(kBlack);

	h_ratio_A2->GetYaxis()->SetRangeUser(0., 0.5);
	h_ratio_A1->GetYaxis()->SetRangeUser(0., 0.5);

	TLegend *leg = new TLegend(0.21,0.58,0.47,0.78);
	leg->AddEntry(h_ratio_A1, "Method A1","p");
	leg->AddEntry(h_ratio_A2, "Method A2","p");
	leg->SetTextSize(0.05);
	// leg->SetTextFont(42);
	leg->SetFillStyle(0);
	leg->SetBorderSize(0);	

	TCanvas *c_temp = new TCanvas();
	c_temp->cd();
	c_temp->SetRightMargin(0.05);

	h_ratio_A2->Draw("PE");
	h_ratio_A1->Draw("PE same");
	leg->Draw();
	TLatex lat;
	lat.SetNDC(kTRUE);
	lat.SetTextSize(0.05);
	lat.DrawLatex(0.23, 0.88, "CMS Preliminary");
	lat.DrawLatex(0.70, 0.88, name);
	lat.DrawLatex(0.23, 0.81, "L_{int.} = 0.98 fb^{-1},   #sqrt{s} = 7 TeV");
	
	Util::PrintNoEPS( c_temp, "FRatio_" + name + "_Eta_A1vsA2", fOutputDir + fOutputSubDir, NULL);
	Util::PrintPDF(   c_temp, "FRatio_" + name + "_Eta_A1vsA2", fOutputDir + fOutputSubDir);
	Util::SaveAsMacro(c_temp, "FRatio_" + name + "_Eta_A1vsA2", fOutputDir + fOutputSubDir);
	delete h_ratio_A1, h_ratio_A2;
	delete c_temp;
	fOutputSubDir = "";
}
void SSDLPlotter::makeRatioPlots(gChannel chan){
	TString name;
	if(chan == Muon)     name = "Muons";
	if(chan == Elec) name = "Electrons";

	fOutputSubDir = "Ratios/" + name + "/";
	char cmd[100];
    sprintf(cmd,"mkdir -p %s%s", fOutputDir.Data(), fOutputSubDir.Data());
    system(cmd);

	vector<int> datasamples;
	vector<int> mcsamples;

	if(chan == Muon){
		datasamples = fMuData;
		mcsamples   = fMCBGMuEnr;
	}
	if(chan == Elec){
		datasamples = fEGData;
		mcsamples   = fMCBG;
	}

	// Customization
	TString axis_name[gNRatioVars] = {"N_{Jets}",  "H_{T} (GeV)", "P_{T}(Hardest Jet) (GeV)", "N_{Vertices}", "p_{T}(Closest Jet) (GeV)", "p_{T}(Away Jet) (GeV)", "N_{BJets}", "E_{T}^{miss} (GeV)", "m_{T} (GeV)"};

	for(size_t i = 0; i < gNRatioVars; ++i){
		TH1D *h_ratio_data = getFRatio(datasamples, chan, i);
		TH1D *h_ratio_mc   = getFRatio(mcsamples,   chan, i);
		h_ratio_data->SetName(Form("FRatio_%s_data", FRatioPlots::var_name[i].Data()));
		h_ratio_mc  ->SetName(Form("FRatio_%s_mc",   FRatioPlots::var_name[i].Data()));

		float max = 0.4;
		if(i==8) max = 1.0;
		if(chan==Elec) max = 0.8;
		h_ratio_data->SetMaximum(max);
		h_ratio_mc  ->SetMaximum(max);
		h_ratio_data->SetMinimum(0.0);
		h_ratio_mc  ->SetMinimum(0.0);

		h_ratio_data->SetXTitle(axis_name[i]);
	    h_ratio_mc  ->SetXTitle(axis_name[i]);
		h_ratio_data->GetYaxis()->SetTitleOffset(1.2);
	    h_ratio_mc  ->GetYaxis()->SetTitleOffset(1.2);
		h_ratio_data->SetYTitle("N_{Tight}/N_{Loose}");
	    h_ratio_mc  ->SetYTitle("N_{Tight}/N_{Loose}");

		h_ratio_data->SetMarkerColor(kBlack);
		h_ratio_data->SetMarkerStyle(20);
		h_ratio_data->SetMarkerSize(1.5);
		h_ratio_data->SetLineWidth(2);
		h_ratio_data->SetLineColor(kBlack);
		h_ratio_data->SetFillColor(kBlack);

		h_ratio_mc  ->SetMarkerColor(kRed);
		h_ratio_mc  ->SetMarkerStyle(23);
		h_ratio_mc  ->SetMarkerSize(1.5);
		h_ratio_mc  ->SetLineWidth(2);
		h_ratio_mc  ->SetLineColor(kRed);
		h_ratio_mc  ->SetFillColor(kRed);

		TLatex *lat = new TLatex();
		lat->SetNDC(kTRUE);
		lat->SetTextColor(kBlack);
		lat->SetTextSize(0.04);
	
		TLegend *leg = new TLegend(0.15,0.75,0.35,0.88);
		leg->AddEntry(h_ratio_data, "Data",         "p");
		leg->AddEntry(h_ratio_mc,   "Simulation",   "p");
		leg->SetFillStyle(0);
		leg->SetTextFont(42);
		leg->SetBorderSize(0);

		TCanvas *c_temp = new TCanvas("C_Temp", "fRatio", 0, 0, 800, 600);
		c_temp->cd();

		// gPad->SetLogy();
		h_ratio_mc  ->DrawCopy("PE X0");
		h_ratio_data->DrawCopy("PE X0 same");
		leg->Draw();
		lat->DrawLatex(0.70,0.92, Form("L_{int.} = %2.1f fb^{-1}", fLumiNorm/1000.));
		lat->DrawLatex(0.11,0.92, name);
		double ymean(0.), yrms(0.);
		getWeightedYMeanRMS(h_ratio_data, ymean, yrms);
		lat->SetTextSize(0.03);
		lat->DrawLatex(0.25,0.92, Form("Mean ratio: %4.2f #pm %4.2f", ymean, yrms));
	
		// Util::PrintNoEPS(c_temp, "FRatio_" + name + "_" + FRatioPlots::var_name[i], fOutputDir + fOutputSubDir, NULL);
		Util::PrintPDF(  c_temp, "FRatio_" + name + "_" + FRatioPlots::var_name[i], fOutputDir + fOutputSubDir);
		delete c_temp, leg, lat;
		delete h_ratio_data, h_ratio_mc;
	}
	fOutputSubDir = "";
}
void SSDLPlotter::makeNTightLoosePlots(gChannel chan){
	TString name;
	if(chan == Muon)     name = "Muons";
	if(chan == Elec) name = "Electrons";

	fOutputSubDir = "Ratios/" + name + "/NTightLoose/";
	char cmd[100];
    sprintf(cmd,"mkdir -p %s%s", fOutputDir.Data(), fOutputSubDir.Data());
    system(cmd);

	vector<int> datasamples;
	vector<int> mcsamples;

	if(chan == Muon){
		datasamples = fMuData;
		mcsamples   = fMCBGMuEnr;
	}
	if(chan == Elec){
		datasamples = fEGData;
		mcsamples   = fMCBG;
	}

	// Customization
	TString axis_name[gNRatioVars] = {"N_{Jets}",  "H_{T} (GeV)", "P_{T}(Hardest Jet) (GeV)", "N_{Vertices}", "p_{T}(Closest Jet) (GeV)", "p_{T}(Away Jet) (GeV)", "N_{BJets}", "E_{T}^{miss} (GeV)", "m_{T} (GeV)"};

	for(size_t i = 0; i < gNRatioVars; ++i){
		THStack *hsntight = new THStack(Form("NTight_%s", FRatioPlots::var_name[i].Data()), "Stack of tight");
		THStack *hsnloose = new THStack(Form("NLoose_%s", FRatioPlots::var_name[i].Data()), "Stack of loose");
		const unsigned int nmcsamples = mcsamples.size();

		// TLegend *leg = new TLegend(0.13,0.60,0.38,0.88);
		TLegend *leg = new TLegend(0.75,0.60,0.89,0.88);
		for(size_t j = 0; j < mcsamples.size(); ++j){
			Sample *S = fSamples[mcsamples[j]];
			FRatioPlots *rat;
			if(chan == Muon) rat = &S->ratioplots[0];
			if(chan == Elec) rat = &S->ratioplots[1];
			rat->ntight[i]->SetFillColor(S->color);
			rat->nloose[i]->SetFillColor(S->color);
			float scale = fLumiNorm / S->lumi;
			rat->ntight[i]->Scale(scale);
			rat->nloose[i]->Scale(scale);
			hsntight->Add(rat->ntight[i]);
			hsnloose->Add(rat->nloose[i]);

			if(rat->nloose[i]->Integral() < 1 ) continue;
			leg->AddEntry(rat->ntight[i], S->sname.Data(), "f");
		}
		hsntight->Draw();
		hsntight->GetXaxis()->SetTitle(axis_name[i]);
		hsnloose->Draw();
		hsnloose->GetXaxis()->SetTitle(axis_name[i]);
		leg->SetFillStyle(0);
		leg->SetTextFont(42);
		leg->SetBorderSize(0);

		TCanvas *c_tight = new TCanvas(Form("NTight_%s", FRatioPlots::var_name[i].Data()), "Tight Stack", 0, 0, 800, 600);
		TCanvas *c_loose = new TCanvas(Form("NLoose_%s", FRatioPlots::var_name[i].Data()), "Loose Stack", 0, 0, 800, 600);


		c_tight->cd();
		gPad->SetLogy();
		hsntight->Draw("hist");
		leg->Draw();

		c_loose->cd();
		gPad->SetLogy();
		hsnloose->Draw("hist");
		leg->Draw();

		Util::PrintNoEPS(c_tight, Form("NTight_%s", FRatioPlots::var_name[i].Data()), fOutputDir + fOutputSubDir, fOutputFile);
		Util::PrintNoEPS(c_loose, Form("NLoose_%s", FRatioPlots::var_name[i].Data()), fOutputDir + fOutputSubDir, fOutputFile);	
		delete hsntight, hsnloose, c_tight, c_loose, leg;
	}
	fOutputSubDir = "";
}

void SSDLPlotter::makePRLPlot1(){
	FakeRatios *FR = new FakeRatios();
	const int nchans = 8;
	TString axis_labels_1[nchans] = {
		"H_{T} >  80",
		"H_{T} > 200",
		"",
		"H_{T} > 450",
		"",
		"H_{T} > 450",
		"",
		""
	};
	TString axis_labels_2[nchans] = {
		"E_{T}^{miss} > 120",
		"E_{T}^{miss} > 120",
		"",
		"E_{T}^{miss} > 50",
		"",
		"E_{T}^{miss} > 120",
		"",
		""
	};
	TH1D    *h_obs      = new TH1D("h_observed",   "Observed number of events", nchans, 0., (double)nchans);
	TH1D    *h_pred_sf  = new TH1D("h_pred_sfake", "Predicted single fakes",    nchans, 0., (double)nchans);
	TH1D    *h_pred_df  = new TH1D("h_pred_dfake", "Predicted double fakes",    nchans, 0., (double)nchans);
	TH1D    *h_pred_cm  = new TH1D("h_pred_chmid", "Predicted charge mis id",   nchans, 0., (double)nchans);
	TH1D    *h_pred_mc  = new TH1D("h_pred_mc",    "Predicted WW/WZ/ZZ",        nchans, 0., (double)nchans);
	TH1D    *h_pred_tot = new TH1D("h_pred_tot",   "Total Prediction",          nchans, 0., (double)nchans);
	THStack *hs_pred    = new THStack("hs_predicted", "Predicted number of events");
	
	const double p_df [nchans] = { 0.1, 0.5, 0.1, 0.8, 0.2, 0.1, 0.0, 0.0};
	const double p_sf [nchans] = {14.4, 15., 8.4, 7.7, 4.4, 2.2, 1.1, 4.0};
	const double p_cm [nchans] = { 0.4, 0.3, 0.3, 0.2, 0.2, 0.1, 0.0, 0.3};
	const double p_mc [nchans] = { 5.8, 5.1, 4.8, 3.0, 2.9, 1.5, 1.4, 4.0};
	const double p_E  [nchans] = { 9.0, 10., 6.4, 6.0, 4.6, 2.7, 2.3, 0.6};
	const double n_obs[nchans] = {19. , 20., 17., 11., 8. , 3. , 2. , 5. };
	
	h_obs->SetMarkerColor(kBlack);
	h_obs->SetMarkerStyle(20);
	h_obs->SetMarkerSize(2.5);
	h_obs->SetLineWidth(2);
	h_obs->SetLineColor(kBlack);
	h_obs->SetFillColor(kBlack);
	
	h_pred_sf->SetLineWidth(1);
	h_pred_df->SetLineWidth(1);
	h_pred_cm->SetLineWidth(1);
	h_pred_mc->SetLineWidth(1);

	Color_t col_sf = 50;
	Color_t col_df = 38;
	Color_t col_cm = 42;
	Color_t col_mc = 31;
	h_pred_sf->SetLineColor(col_sf);
	h_pred_sf->SetFillColor(col_sf);
	h_pred_df->SetLineColor(col_df);
	h_pred_df->SetFillColor(col_df);
	h_pred_cm->SetLineColor(col_cm);
	h_pred_cm->SetFillColor(col_cm);
	h_pred_mc->SetLineColor(col_mc);
	h_pred_mc->SetFillColor(col_mc);

	h_pred_tot  ->SetLineWidth(1);
	h_pred_tot  ->SetFillColor(12);
	h_pred_tot  ->SetFillStyle(3005);
	
	// Add numbers:
	for(size_t i = 0; i < nchans; ++i){
		h_obs     ->SetBinContent(i+1, n_obs[i]);
		h_pred_sf ->SetBinContent(i+1, p_sf [i]);
		h_pred_df ->SetBinContent(i+1, p_df [i]);
		h_pred_cm ->SetBinContent(i+1, p_cm [i]);
		h_pred_mc ->SetBinContent(i+1, p_mc [i]);
		h_pred_sf ->GetXaxis()->SetBinLabel(i+1, "");
		h_pred_tot->SetBinError(i+1, p_E[i]);
		h_obs     ->SetBinError(i+1, FR->getEStat(n_obs[i]));
	}
	
	h_pred_tot->Add(h_pred_sf);
	h_pred_tot->Add(h_pred_df);
	h_pred_tot->Add(h_pred_cm);
	h_pred_tot->Add(h_pred_mc);

	hs_pred->Add(h_pred_sf);
	hs_pred->Add(h_pred_df);
	hs_pred->Add(h_pred_cm);
	hs_pred->Add(h_pred_mc);

	// double max = h_obs->Integral()*0.5;
	double max = 43.;
	h_obs     ->SetMaximum(max);
	h_pred_sf ->SetMaximum(max);
	h_pred_df ->SetMaximum(max);
	h_pred_cm ->SetMaximum(max);
	h_pred_mc ->SetMaximum(max);
	h_pred_tot->SetMaximum(max);
	hs_pred   ->SetMaximum(max);
	
	hs_pred->Draw("goff");
	hs_pred->GetYaxis()->SetTitle("Events");
	hs_pred->GetYaxis()->SetTitleOffset(1.1);
	hs_pred->GetYaxis()->SetTitleSize(0.05);
	hs_pred->GetYaxis()->SetTickLength(0.02);
	// hs_pred->GetXaxis()->SetLabelOffset(0.01);
	// hs_pred->GetXaxis()->SetLabelFont(42);
	// hs_pred->GetXaxis()->SetLabelSize(0.05);
	// hs_pred->GetXaxis()->LabelsOption("v");
	for(size_t i = 0; i < nchans; ++i) hs_pred->GetXaxis()->SetBinLabel(i+1, "");
	
	TLegend *leg = new TLegend(0.58,0.66,0.93,0.89);
	// TLegend *leg = new TLegend(0.27,0.65,0.62,0.88);
	// TLegend *leg = new TLegend(0.15,0.65,0.50,0.88);
	leg->AddEntry(h_obs,      "Observed",         "p");
	leg->AddEntry(h_pred_mc,  "Irreducible (MC)", "f");
	leg->AddEntry(h_pred_cm,  "Charge MisID",     "f");
	leg->AddEntry(h_pred_df,  "Double Fakes",     "f");
	leg->AddEntry(h_pred_sf,  "Single Fakes",     "f");
	leg->AddEntry(h_pred_tot, "Total Uncertainty","f");
	leg->SetFillStyle(0);
	leg->SetTextFont(42);
	// leg->SetFillColor(kWhite);
	// leg->SetTextSize(0.05);
	leg->SetBorderSize(0);
	
	TCanvas *c_temp = new TCanvas("C_ObsPred", "Observed vs Predicted", 0, 0, 600, 600);
	c_temp->cd();
	c_temp->SetBottomMargin(0.2);
	c_temp->SetLeftMargin(0.12);
	c_temp->SetRightMargin(0.05);
	hs_pred->Draw("hist");
	h_pred_tot->DrawCopy("0 E2 same");
	h_obs->DrawCopy("PE X0 same");
	leg->Draw();

	for(size_t i = 0; i < nchans; ++i){
		float bincenter = hs_pred->GetXaxis()->GetBinCenter(i+1);
		// float shift = 0.01;
		// float split = 0.03;
		// float x1 = bincenter - split + shift;
		// float x2 = bincenter + split + shift;
		float shift = 0.08;
		float split = 0.22;
		float x1 = bincenter - split + shift;
		float x2 = bincenter + split + shift;
		TLatex *label = new TLatex();
		// label->SetNDC(0);
		label->SetTextColor(kBlack);
		label->SetTextSize(0.035);
		// label->SetTextAlign();
		label->SetTextAngle(90);
		// label->DrawLatex(x1, -9., axis_labels_1[i]);
		// label->DrawLatex(x2, -9., axis_labels_2[i]);
		label->DrawLatex(x1, -11., axis_labels_1[i]);
		label->DrawLatex(x2, -11., axis_labels_2[i]);
	}
	TLine *l1 = new TLine(1.0, max*1.05, 1.0, -11.);
	TLine *l2 = new TLine(3.0, max*1.05, 3.0, -11.);
	TLine *l3 = new TLine(5.0, max*0.65, 5.0, -11.);
	l1->Draw();
	l2->Draw();
	l3->Draw();

	TLatex *lat = new TLatex();
	// lat->SetNDC(kTRUE);
	lat->SetTextColor(kBlack);
	lat->SetTextSize(0.03);
	
	float shift = 0.38;
	float split = 0.17;
	lat->SetTextAngle(90);
	lat->DrawLatex(1.0 - shift - split, 32.0, "High p_{T}");
	lat->DrawLatex(1.0 - shift + split, 32.0, "(ee/#mu#mu/e#mu)");
	lat->DrawLatex(2.0 - shift - split, 33.0, "Low p_{T}" );
	lat->DrawLatex(2.0 - shift + split, 33.0, "(ee/#mu#mu/e#mu)" );
	lat->DrawLatex(3.0 - shift - split, 22.0, "High p_{T}");
	lat->DrawLatex(3.0 - shift + split, 22.0, "(ee/#mu#mu/e#mu)");
	lat->DrawLatex(4.0 - shift - split, 19.5, "Low p_{T}" );
	lat->DrawLatex(4.0 - shift + split, 19.5, "(ee/#mu#mu/e#mu)" );
	lat->DrawLatex(5.0 - shift - split, 14.0, "High p_{T}");
	lat->DrawLatex(5.0 - shift + split, 14.0, "(ee/#mu#mu/e#mu)");
	lat->DrawLatex(6.0 - shift - split,  8.0, "Low p_{T}" );
	lat->DrawLatex(6.0 - shift + split,  8.0, "(ee/#mu#mu/e#mu)" );
	lat->DrawLatex(7.0 - shift - split,  6.0, "High p_{T}");
	lat->DrawLatex(7.0 - shift + split,  6.0, "(ee/#mu#mu/e#mu)");
	lat->DrawLatex(8.0 - shift - split, 12.0, "Tau channels");
	lat->DrawLatex(8.0 - shift + split, 12.0, "(e#tau/#mu#tau/#tau#tau)");
	// float shift = 0.4;
	// lat->DrawLatex(1.0 - shift, 32.0, "High p_{T} (e/#mu)");
	// lat->DrawLatex(2.0 - shift, 33.0, "Low p_{T} (e/#mu)" );
	// lat->DrawLatex(3.0 - shift, 22.0, "High p_{T} (e/#mu)");
	// lat->DrawLatex(4.0 - shift, 19.5, "Low p_{T} (e/#mu)" );
	// lat->DrawLatex(5.0 - shift, 14.0, "High p_{T} (e/#mu)");
	// lat->DrawLatex(6.0 - shift,  8.0, "Low p_{T} (e/#mu)" );
	// lat->DrawLatex(7.0 - shift,  6.0, "High p_{T} (e/#mu)");
	// lat->DrawLatex(8.0 - shift, 12.0, "e#tau/#mu#tau/#tau#tau");
	drawTopLine();
	
	gPad->RedrawAxis();
	fOutputSubDir = "PRL";
	Util::PrintPDF(c_temp, "ObsPred_MultiChan", fOutputDir + fOutputSubDir);
	delete c_temp;	
	delete h_obs, h_pred_sf, h_pred_df, h_pred_cm, h_pred_mc, h_pred_tot, hs_pred;
	delete FR;	
}
void SSDLPlotter::makeIsoVsMETPlot(gSample sample){
	fOutputSubDir = "IsoVsMETPlots/";
	fCurrentSample = sample;
	fCurrentChannel = Muon;
	const int nisobins = 20;
	const int nmetbins = 7;
	float metbins[nmetbins+1] = {0., 5., 10., 15., 20., 30., 40., 1000.};
	float ratio[nmetbins];
	// Color_t colors[nmetbins] = {1, 1, 1, 1, 1};
	Color_t colors[nmetbins] = {1, 12, 39, 38, 32, 30, 29};
	// Color_t colors[nmetbins] = {52, 63, 74, 85, 96};
	Style_t styles[nmetbins] = {20, 21, 22, 23, 34, 33, 29};
	TH1F* hiso[nmetbins];
	
	TTree *tree = fSamples[sample]->getTree();
	for(int i = 0; i < nmetbins; ++i){
		TString name = Form("h%d",i+1);
		hiso[i] = new TH1F(name, Form("MET%.0fto%.0f", metbins[i], metbins[i+1]), nisobins, 0., 1.);
		hiso[i]->Sumw2();
		hiso[i]->SetLineWidth(1);
		hiso[i]->SetLineColor(colors[i]);
		hiso[i]->SetMarkerColor(colors[i]);
		hiso[i]->SetMarkerStyle(styles[i]);
		hiso[i]->SetMarkerSize(1.3);
		hiso[i]->SetFillStyle(0);
		hiso[i]->SetXTitle("RelIso(#mu)");
	}
		
	tree->ResetBranchAddresses();
	if(fSamples[sample]->datamc < 1) Init(tree);
	if(fSamples[sample]->datamc > 0) InitMC(tree);
	fSample = fSamples[sample];
	for (Long64_t jentry=0; jentry<tree->GetEntriesFast();jentry++) {
		tree->GetEntry(jentry);
		printProgress(jentry, tree->GetEntriesFast(), fSamples[sample]->name);
		if(singleMuTrigger() == false) continue;
		int mu1(-1), mu2(-1);
		if(hasLooseMuons(mu1, mu2) < 1) continue;
		setHypLepton1(mu1, Muon);
		if(!passesJet50Cut())           continue;
		if(!passesNJetCut(1))           continue;
		if(MuMT[mu1] > 20.)             continue;
		if(NMus > 1)                    continue;
		for(size_t j = 0; j < nmetbins; ++j) if(pfMET > metbins[j] && pfMET < metbins[j+1]) hiso[j]->Fill(MuIso[mu1], singleMuPrescale());
	}

	for(int i = 0; i < nmetbins; ++i){
		ratio[i] = hiso[i]->Integral(1, hiso[i]->FindBin(0.1499))/hiso[i]->Integral(1, nisobins);
		hiso[i]->Scale(1./hiso[i]->Integral());
	}

	TLegend *leg = new TLegend(0.35,0.15,0.70,0.40);
	// TLegend *leg = new TLegend(0.60,0.65,0.88,0.88);
	for(int i = 0; i < nmetbins; ++i) leg->AddEntry(hiso[i], Form("T/L = %5.2f: %.0f < E_{T}^{miss} < %.0f", ratio[i], metbins[i], metbins[i+1]), "p");
	leg->SetFillStyle(0);
	leg->SetTextFont(42);
	// leg->SetTextSize(0.03);
	leg->SetBorderSize(0);

	TLine *line = new TLine(0.15, 0.00, 0.15, 0.08);
	line->SetLineStyle(3);

	TCanvas *c_temp = new TCanvas("C_HTvsMET", "HT vs MET in Data vs MC", 0, 0, 800, 600);
	c_temp->cd();
	hiso[0]->Draw("axis");
	hiso[0]->GetYaxis()->SetRangeUser(0., 0.08);
	for(int i = 0; i < nmetbins; ++i) hiso[i]->DrawCopy("PE X0 same");
	leg->Draw();
	line->Draw();
	fLatex->DrawLatex(0.10,0.92, fSamples[sample]->name);
	Util::PrintNoEPS( c_temp, "IsoVsMET_" + fSamples[sample]->sname, fOutputDir + fOutputSubDir, NULL);
	Util::PrintPDF(   c_temp, "IsoVsMET_" + fSamples[sample]->sname, fOutputDir + fOutputSubDir);

	// Cleanup
	delete leg, c_temp;
	for(int i = 0; i < nmetbins; ++i) hiso[i]->Delete();
	fOutputSubDir = "";
	fSamples[sample]->cleanUp();
}

//____________________________________________________________________________
void SSDLPlotter::fillRatios(vector<int> musamples, vector<int> elsamples, int datamc){
	if(datamc == 0){
		fH1D_MufRatio = fillRatioPt(Muon, musamples, SigSup, false);
		fH1D_MupRatio = fillRatioPt(Muon, musamples, ZDecay, false);
		fH1D_ElfRatio = fillRatioPt(Elec, elsamples, SigSup, false);
		fH1D_ElpRatio = fillRatioPt(Elec, elsamples, ZDecay, false);
	}
	if(datamc == 1){
		fH1D_MufRatio_MC = fillRatioPt(Muon, musamples, SigSup, false);
		fH1D_MupRatio_MC = fillRatioPt(Muon, musamples, ZDecay, false);
		fH1D_ElfRatio_MC = fillRatioPt(Elec, elsamples, SigSup, false);
		fH1D_ElpRatio_MC = fillRatioPt(Elec, elsamples, ZDecay, false);
	}
}
TH1D* SSDLPlotter::fillRatioPt(gChannel chan, int sample, gFPSwitch fp, bool output){
	vector<int> samples; samples.push_back(sample);
	return fillRatioPt(chan, samples, fp);
}
TH1D* SSDLPlotter::fillRatioPt(gChannel chan, vector<int> samples, gFPSwitch fp, bool output){
	gStyle->SetOptStat(0);
	TString shortname[2] = {"Mu", "El"};
	TString longname[2] = {"Muons", "Electrons"};
	int muelswitch = 0;
	if(chan == Elec) muelswitch = 1;
	TH2D *h_2d  = new TH2D(Form("%sRatio",   shortname[muelswitch].Data()), Form("Ratio of tight to loose %s vs Pt vs Eta", longname[muelswitch].Data()), getNPt2Bins(chan), getPt2Bins(chan), getNEtaBins(chan), getEtaBins(chan));
	TH1D *h_pt  = new TH1D(Form("%sRatioPt" ,shortname[muelswitch].Data()), Form("Ratio of tight to loose %s vs Pt"       , longname[muelswitch].Data()), getNPt2Bins(chan), getPt2Bins(chan));
	TH1D *h_eta = new TH1D(Form("%sRatioEta",shortname[muelswitch].Data()), Form("Ratio of tight to loose %s vs Eta"      , longname[muelswitch].Data()), getNEtaBins(chan), getEtaBins(chan));

	h_pt->SetXTitle("p_{#perp} (GeV)");
	h_pt->SetYTitle("# Tight / # Loose");
	h_pt->GetYaxis()->SetTitleOffset(1.2);

	calculateRatio(samples, chan, fp, h_2d, h_pt, h_eta, output);
	delete h_2d, h_eta;
	return h_pt;
}

//____________________________________________________________________________
void SSDLPlotter::calculateRatio(vector<int> samples, gChannel chan, gFPSwitch fp, TH2D*& h_2d, bool output){
	TH1D *h_dummy1 = new TH1D("dummy1", "dummy1", 1, 0.,1.);
	TH1D *h_dummy2 = new TH1D("dummy2", "dummy2", 1, 0.,1.);
	calculateRatio(samples, chan, fp, h_2d, h_dummy1, h_dummy2, output);
	delete h_dummy1, h_dummy2;
}
void SSDLPlotter::calculateRatio(vector<int> samples, gChannel chan, gFPSwitch fp, TH2D*& h_2d, TH1D*& h_pt, TH1D*&h_eta, bool output){
/*
TODO Fix treatment of statistical errors and luminosity scaling here!
*/
	gStyle->SetOptStat(0);

	h_2d->Sumw2();
	h_pt->Sumw2();
	h_eta->Sumw2();

	TH2D *H_ntight = new TH2D("NTight", "NTight Muons", h_2d->GetNbinsX(), h_2d->GetXaxis()->GetXbins()->GetArray(), h_2d->GetNbinsY(),  h_2d->GetYaxis()->GetXbins()->GetArray());
	TH2D *H_nloose = new TH2D("NLoose", "NLoose Muons", h_2d->GetNbinsX(), h_2d->GetXaxis()->GetXbins()->GetArray(), h_2d->GetNbinsY(),  h_2d->GetYaxis()->GetXbins()->GetArray());
	H_ntight->Sumw2();
	H_nloose->Sumw2();

	getPassedTotal(samples, chan, fp, H_ntight, H_nloose, output);
	h_2d->Divide(H_ntight, H_nloose, 1., 1., "B");

	TH1D *hmuloosept  = H_nloose->ProjectionX();
	TH1D *hmulooseeta = H_nloose->ProjectionY();
	TH1D *hmutightpt  = H_ntight->ProjectionX();
	TH1D *hmutighteta = H_ntight->ProjectionY();

	h_pt ->Divide(hmutightpt,  hmuloosept,  1., 1., "B"); // binomial
	h_eta->Divide(hmutighteta, hmulooseeta, 1., 1., "B"); // weights are ignored
	delete H_ntight, H_nloose, hmuloosept, hmulooseeta, hmutightpt, hmutighteta;
	TString name = "";
	for(size_t i = 0; i < samples.size(); ++i){
		int sample = samples[i];
		name += h_2d->GetName();
		name += "_";
		name += fSamples[sample]->sname;
	}
	if(output){
		printObject(h_2d,  TString("Ratio")    + name, "colz");
		printObject(h_pt,  TString("RatioPt")  + name, "PE1");
		printObject(h_eta, TString("RatioEta") + name, "PE1");
	}
}
void SSDLPlotter::calculateRatio(vector<int> samples, gChannel chan, gFPSwitch fp, float &ratio, float &ratioe){
	double ntight(0.), nloose(0.);
	double ntighte2(0.), nloosee2(0.);
	vector<int> v_ntight, v_nloose;
	vector<float> v_scale;
	vector<TString> v_name;
	for(size_t i = 0; i < samples.size(); ++i){
		int index = samples[i];
		Sample *S = fSamples[index];

		int ntight_sam(0), nloose_sam(0);
		v_name.push_back(S->sname);

		float scale = fLumiNorm/S->lumi; // Normalize all
		if(S->datamc == 0) scale = 1;
		if(fp == SigSup){
			ntight += scale * S->numbers[Baseline][chan].nsst;
			nloose += scale * S->numbers[Baseline][chan].nssl;

			ntight_sam += S->numbers[Baseline][chan].nsst;
			nloose_sam += S->numbers[Baseline][chan].nssl;
		}
		if(fp == ZDecay){
			ntight += scale * S->numbers[Baseline][chan].nzt;
			nloose += scale * S->numbers[Baseline][chan].nzl;

			ntight_sam += S->numbers[Baseline][chan].nzt;
			nloose_sam += S->numbers[Baseline][chan].nzl;
		}
		v_ntight.push_back(ntight_sam);
		v_nloose.push_back(nloose_sam);
		v_scale.push_back(scale);
	}
	ratioWithBinomErrors(ntight, nloose, ratio, ratioe);
	// ratioWithPoissErrors(ntight, nloose, ratio, ratioe);
	if(fVerbose > 2){
		cout << "--------------------------------------------------------" << endl;
		TString s_forp;
		if(fp == SigSup) s_forp = "f Ratio";
		if(fp == ZDecay) s_forp = "p Ratio";
		TString s_channel;
		if(chan == Muon)     s_channel = "Muon";
		if(chan == Elec) s_channel = "Electron";
		cout << " Calling calculateRatio for " << s_forp << " in " << s_channel << " channel..." << endl;
		for(size_t i = 0; i < v_ntight.size(); ++i){
			cout << setw(9) << v_name[i] << ": ";
			cout << " Nt:    " << setw(8) << setprecision(2) << v_ntight[i];
			cout << " Nl:    " << setw(8) << setprecision(2) << v_nloose[i];
			cout << " Scale: " << v_scale[i] << endl;
		}
		cout << "--------------------------------------------------------" << endl;
		cout << " Total: ";
		cout << " Nt:    " << setw(8) << ntight;
		cout << " Nl:    " << setw(8) << nloose;
		cout << " Ratio: " << ratio << endl;
		cout << "--------------------------------------------------------" << endl;
	}
}
void SSDLPlotter::calculateRatio(vector<int> samples, gChannel chan, gFPSwitch fp, float &ratio, float &ratioeup, float &ratioelow){
	// Careful, this method only takes integer numbers for passed/total events, therefore
	// only makes sense for application on data right now.
	int ntight(0), nloose(0);
	float ntighte2(0.), nloosee2(0.);
	vector<int> v_ntight, v_nloose;
	vector<TString> v_name;
	for(size_t i = 0; i < samples.size(); ++i){
		Sample *S = fSamples[samples[i]];

		int ntight_sam(0), nloose_sam(0);
		v_name.push_back(S->sname);

		float scale = fLumiNorm/S->lumi; // Normalize all
		if(S->datamc == 0) scale = 1;
		if(fp == SigSup){
			ntight += scale * S->numbers[Baseline][chan].nsst;
			nloose += scale * S->numbers[Baseline][chan].nssl;

			ntight_sam += S->numbers[Baseline][chan].nsst;
			nloose_sam += S->numbers[Baseline][chan].nssl;
		}
		if(fp == ZDecay){
			ntight += S->numbers[Baseline][chan].nzt;
			nloose += S->numbers[Baseline][chan].nzl;

			ntight_sam += S->numbers[Baseline][chan].nzt;
			nloose_sam += S->numbers[Baseline][chan].nzl;
		}
		v_ntight.push_back(ntight_sam);
		v_nloose.push_back(nloose_sam);
	}
	ratioWithAsymmCPErrors(ntight, nloose, ratio, ratioeup, ratioelow);
	if(fVerbose > 2){
		cout << "--------------------------------------------------------" << endl;
		TString s_forp;
		if(fp == SigSup) s_forp = "f Ratio";
		if(fp == ZDecay) s_forp = "p Ratio";
		TString s_channel;
		if(chan == Muon)     s_channel = "Muon";
		if(chan == Elec) s_channel = "Electron";
		cout << " Calling calculateRatio for " << s_forp << " in " << s_channel << " channel..." << endl;
		for(size_t i = 0; i < v_ntight.size(); ++i){
			cout << setw(9) << v_name[i] << ": ";
			cout << " Nt:    " << setw(8) << setprecision(2) << v_ntight[i];
			cout << " Nl:    " << setw(8) << setprecision(2) << v_nloose[i];
			cout << endl;
		}
		cout << "--------------------------------------------------------" << endl;
		cout << " Total: ";
		cout << " Nt:    " << setw(8) << ntight;
		cout << " Nl:    " << setw(8) << nloose;
		cout << " Ratio: " << ratio << " + " << ratioeup << " - " << ratioelow << endl;
		cout << "--------------------------------------------------------" << endl;
	}
}

//____________________________________________________________________________
void SSDLPlotter::getPassedTotal(vector<int> samples, gChannel chan, gFPSwitch fp, TH2D*& h_passed, TH2D*& h_total, bool output, gHiLoSwitch hilo){
	// toggle: choose binning: 0: pt, default, 1: nvrtx, 2: closest jet pt
	if(fVerbose>2) cout << "---------------" << endl;
	for(size_t i = 0; i < samples.size(); ++i){
		Sample *S = fSamples[samples[i]];

		float scale = fLumiNorm / S->lumi;
		if(S->datamc == 0) scale = 1;

		Channel *C;
		if(chan == Muon) C = &S->region[Baseline][hilo].mm;
		if(chan == Elec) C = &S->region[Baseline][hilo].ee;
		TH2D *ntight, *nloose;
		if(fp == SigSup){
			ntight = C->fntight;
			nloose = C->fnloose;
		} else if(fp == ZDecay){
			ntight = C->pntight;
			nloose = C->pnloose;
		}
		h_passed->Add(ntight, scale);
		h_total ->Add(nloose, scale);
	}
	TString name = "";
	for(size_t i = 0; i < samples.size(); ++i){
		int sample = samples[i];
		if(i > 0) name += "_";
		name += fSamples[sample]->sname;
	}
	if(output){
		printObject(h_passed, TString("Passed") + name, "colz");
		printObject(h_total,  TString("Total")  + name, "colz");
	}	
}
TH1D* SSDLPlotter::getFRatio(vector<int> samples, gChannel chan, int ratiovar, bool output){
	gStyle->SetOptStat(0);

	TH1D *ntight = new TH1D("h_NTight",  "NTight",  FRatioPlots::nbins[ratiovar], FRatioPlots::xmin[ratiovar], FRatioPlots::xmax[ratiovar]);
	TH1D *nloose = new TH1D("h_NLoose",  "NLoose",  FRatioPlots::nbins[ratiovar], FRatioPlots::xmin[ratiovar], FRatioPlots::xmax[ratiovar]);
	TH1D *ratio  = new TH1D("h_TLRatio", "TLRatio", FRatioPlots::nbins[ratiovar], FRatioPlots::xmin[ratiovar], FRatioPlots::xmax[ratiovar]);
	ntight->Sumw2(); nloose->Sumw2(); ratio->Sumw2();

	for(size_t i = 0; i < samples.size(); ++i){
		Sample *S = fSamples[samples[i]];

		float scale = fLumiNorm / S->lumi;
		if(S->datamc == 0) scale = 1;

		FRatioPlots *RP;
		if(chan == Muon)     RP = &S->ratioplots[0];
		if(chan == Elec) RP = &S->ratioplots[1];
		ntight->Add(RP->ntight[ratiovar], scale);
		nloose->Add(RP->nloose[ratiovar], scale);
	}
	
	ratio->Divide(ntight, nloose, 1., 1., "B");

	delete ntight, nloose;
	return ratio;
}

//____________________________________________________________________________
void SSDLPlotter::ratioWithBinomErrors(float ntight, float nloose, float &ratio, float &error){
	ratio = ntight/nloose;
	error = TMath::Sqrt( ntight*(1.0-ntight/nloose) ) / nloose;                  // Binomial
}
void SSDLPlotter::ratioWithPoissErrors(float ntight, float nloose, float &ratio, float &error){
	ratio = ntight/nloose;
	error = TMath::Sqrt( ntight*ntight*(nloose+ntight)/(nloose*nloose*nloose) ); // Poissonian	
}
void SSDLPlotter::ratioWithAsymmCPErrors(int passed, int total, float &ratio, float &upper, float &lower){
	TEfficiency *eff = new TEfficiency("TempEfficiency", "TempEfficiency", 1, 0., 1.);
	eff->SetStatisticOption(TEfficiency::kFCP); // Frequentist Clopper Pearson = default
	eff->SetConfidenceLevel(0.683); // 1-sigma = default
	if( eff->SetTotalEvents(1, total) && eff->SetPassedEvents(1, passed) ){
		ratio = eff->GetEfficiency(1);
		upper = eff->GetEfficiencyErrorUp(1);
		lower = eff->GetEfficiencyErrorLow(1);
	}
	else{
		ratio = 1;
		upper = 1;
		lower = 0;
	};
	delete eff;
}

//____________________________________________________________________________
void SSDLPlotter::makeAllIntPredictions(){
	TString outputdir = Util::MakeOutputDir(fOutputDir + "IntPredictions");
	fOutputSubDir = "IntPredictions/";
	time_t rawtime;
	struct tm* timeinfo;
	time(&rawtime);
	timeinfo = localtime(&rawtime);
	// access a chararray containing the date with asctime(timeinfo)
	
	TString tablefilename = outputdir + "Table2.tex";
	// TString didarfilename = fOutputDir + fOutputSubDir + "forDidar.txt";
	TString notetable     = outputdir + "NoteTable.tex";
	fOUTSTREAM.open(tablefilename.Data(), ios::trunc);
	fOUTSTREAM << "==========================================================================================================" << endl;
	fOUTSTREAM << " Table 2 inputs from ETH Analysis" << endl;
	fOUTSTREAM << Form(" Generated on: %s ", asctime(timeinfo)) << endl;
	fOUTSTREAM << endl;
	
	// fOUTSTREAM2.open(didarfilename.Data(), ios::trunc);
	// fOUTSTREAM2 << "////////////////////////////////////////////////////////////////////" << endl;
	// fOUTSTREAM2 << "// Plot inputs from ETH Analysis" << endl;
	// fOUTSTREAM2 << Form("// Generated on: %s ", asctime(timeinfo)) << endl;
	// fOUTSTREAM2 << "// Format is {ee, mm, em, total}" << endl;
	// fOUTSTREAM2 << "// Errors are on sum of backgrounds" << endl;
	// fOUTSTREAM2 << endl;

	fOUTSTREAM3.open(notetable.Data(), ios::trunc);
	fOUTSTREAM3 << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
	fOUTSTREAM3 << Form("%%%% Generated on: %s ", asctime(timeinfo)) << endl;
	fOUTSTREAM3 << "%% Format is tot, (ee, mm, em)" << endl;
	fOUTSTREAM3 << endl;

	for(size_t i = 0; i < gNREGIONS; ++i){
		TString outputname = outputdir + "DataPred_" + Region::sname[i] + ".txt";
		makeIntPrediction(outputname, gRegion(i));
	}

	fOUTSTREAM.close();
	fOUTSTREAM2.close();
	fOutputSubDir = "";
}
void SSDLPlotter::makeIntPrediction(TString filename, gRegion reg, gHiLoSwitch hilo){
	ofstream OUT(filename.Data(), ios::trunc);

	TLatex *lat = new TLatex();
	lat->SetNDC(kTRUE);
	lat->SetTextColor(kBlack);
	lat->SetTextSize(0.04);

	vector<int> musamples;
	vector<int> elsamples;
	vector<int> emusamples;
	
	const float RareESyst  = 0.5;
	const float RareESyst2 = RareESyst*RareESyst;
	
	const float FakeESyst  = 0.5;
	const float FakeESyst2 = FakeESyst*FakeESyst;

	musamples = fMuData;
	elsamples = fEGData;
	emusamples = fMuEGData;

	bool diffratio = true; // use flat or differential ratios

	OUT << "/////////////////////////////////////////////////////////////////////////////" << endl;
	OUT << " Producing integrated predictions for region " << Region::sname[reg] << endl;
	OUT << "  scaling MC to " << fLumiNorm << " /pb" << endl << endl;

	///////////////////////////////////////////////////////////////////////////////////
	// RATIOS /////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////
	float mufratio_data(0.),  mufratio_data_e(0.);
	float mupratio_data(0.),  mupratio_data_e(0.);
	float elfratio_data(0.),  elfratio_data_e(0.);
	float elpratio_data(0.),  elpratio_data_e(0.);
	float mufratio_allmc(0.), mufratio_allmc_e(0.);
	float mupratio_allmc(0.), mupratio_allmc_e(0.);
	float elfratio_allmc(0.), elfratio_allmc_e(0.);
	float elpratio_allmc(0.), elpratio_allmc_e(0.);

	calculateRatio(fMuData, Muon, SigSup, mufratio_data, mufratio_data_e);
	calculateRatio(fMuData, Muon, ZDecay, mupratio_data, mupratio_data_e);

	calculateRatio(fEGData, Elec, SigSup, elfratio_data, elfratio_data_e);
	calculateRatio(fEGData, Elec, ZDecay, elpratio_data, elpratio_data_e);

	calculateRatio(fMCBGMuEnr, Muon, SigSup, mufratio_allmc, mufratio_allmc_e);
	calculateRatio(fMCBGMuEnr, Muon, ZDecay, mupratio_allmc, mupratio_allmc_e);
	calculateRatio(fMCBG,      Elec, SigSup, elfratio_allmc, elfratio_allmc_e);
	calculateRatio(fMCBG,      Elec, ZDecay, elpratio_allmc, elpratio_allmc_e);

	///////////////////////////////////////////////////////////////////////////////////
	// OBSERVATIONS ///////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////
	float nt2_mm(0.), nt10_mm(0.), nt0_mm(0.);
	float nt2_em(0.), nt10_em(0.), nt01_em(0.), nt0_em(0.);
	float nt2_ee(0.), nt10_ee(0.), nt0_ee(0.);

	// FR Predictions from event-by-event weights (pre stored)
	float npp_mm(0.), npf_mm(0.), nff_mm(0.);
	float npp_em(0.), npf_em(0.), nfp_em(0.), nff_em(0.);
	float npp_ee(0.), npf_ee(0.), nff_ee(0.);

	// OS yields
	float nt2_ee_BB_os(0.), nt2_ee_EE_os(0.), nt2_ee_EB_os(0.);
	float nt2_em_BB_os(0.), nt2_em_EE_os(0.);

	for(size_t i = 0; i < musamples.size(); ++i){
		Sample *S = fSamples[musamples[i]];
		nt2_mm  += S->numbers[reg][Muon].nt2;
		nt10_mm += S->numbers[reg][Muon].nt10;
		nt0_mm  += S->numbers[reg][Muon].nt0;
		
		if(diffratio){
			npp_mm += S->numbers[reg][Muon].npp;
			npf_mm += S->numbers[reg][Muon].npf + S->numbers[reg][Muon].nfp;
			nff_mm += S->numbers[reg][Muon].nff;			
		}
	}
	for(size_t i = 0; i < emusamples.size(); ++i){
		Sample *S = fSamples[emusamples[i]];
		nt2_em  += S->numbers[reg][ElMu].nt2;
		nt10_em += S->numbers[reg][ElMu].nt10;
		nt01_em += S->numbers[reg][ElMu].nt01;
		nt0_em  += S->numbers[reg][ElMu].nt0;

		if(diffratio){
			npp_em += S->numbers[reg][ElMu].npp;
			npf_em += S->numbers[reg][ElMu].npf;
			nfp_em += S->numbers[reg][ElMu].nfp;
			nff_em += S->numbers[reg][ElMu].nff;
		}

		nt2_em_BB_os += S->region[reg][hilo].em.nt20_OS_BB_pt->GetEntries(); // ele in barrel
		nt2_em_EE_os += S->region[reg][hilo].em.nt20_OS_EE_pt->GetEntries(); // ele in endcal
	}
	for(size_t i = 0; i < elsamples.size(); ++i){
		Sample *S = fSamples[elsamples[i]];
		nt2_ee  += S->numbers[reg][Elec].nt2;
		nt10_ee += S->numbers[reg][Elec].nt10;
		nt0_ee  += S->numbers[reg][Elec].nt0;

		if(diffratio){
			npp_ee += S->numbers[reg][Elec].npp;
			npf_ee += S->numbers[reg][Elec].npf + S->numbers[reg][Muon].nfp;
			nff_ee += S->numbers[reg][Elec].nff;
		}
		
		nt2_ee_BB_os += S->region[reg][hilo].ee.nt20_OS_BB_pt->GetEntries(); // both in barrel
		nt2_ee_EE_os += S->region[reg][hilo].ee.nt20_OS_EE_pt->GetEntries(); // both in endcal
		nt2_ee_EB_os += S->region[reg][hilo].ee.nt20_OS_EB_pt->GetEntries(); // one barrel, one endcap
	}

	///////////////////////////////////////////////////////////////////////////////////
	// PRINTOUT ///////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////
	OUT << "---------------------------------------------------------------------------------------------------------" << endl;
	OUT << "         RATIOS  ||     Mu-fRatio      |     Mu-pRatio      ||     El-fRatio      |     El-pRatio      ||" << endl;
	OUT << "---------------------------------------------------------------------------------------------------------" << endl;
	OUT << setw(16) << "        allMC    ||";
	OUT << setw(7)  << setprecision(3) << mufratio_allmc << " +/- " << setw(7) << setprecision(3) << mufratio_allmc_e << " |";
	OUT << setw(7)  << setprecision(3) << mupratio_allmc << " +/- " << setw(7) << setprecision(3) << mupratio_allmc_e << " ||";
	OUT << setw(7)  << setprecision(3) << elfratio_allmc << " +/- " << setw(7) << setprecision(3) << elfratio_allmc_e << " |";
	OUT << setw(7)  << setprecision(3) << elpratio_allmc << " +/- " << setw(7) << setprecision(3) << elpratio_allmc_e << " ||";
	OUT << endl;
	OUT << setw(16) << "  data stat only ||";
	OUT << setw(7)  << setprecision(3) << mufratio_data  << " +/- " << setw(7) << setprecision(3) << mufratio_data_e  << " |";
	OUT << setw(7)  << setprecision(3) << mupratio_data  << " +/- " << setw(7) << setprecision(3) << mupratio_data_e  << " ||";
	OUT << setw(7)  << setprecision(3) << elfratio_data  << " +/- " << setw(7) << setprecision(3) << elfratio_data_e  << " |";
	OUT << setw(7)  << setprecision(3) << elpratio_data  << " +/- " << setw(7) << setprecision(3) << elpratio_data_e  << " ||";
	OUT << endl;
	OUT << "---------------------------------------------------------------------------------------------------------" << endl << endl;
	OUT << "-------------------------------------------------------------------------------------------------------------" << endl;
	OUT << "                 |           Mu/Mu          |                E/Mu               |           E/E            ||" << endl;
	OUT << "         YIELDS  |   Ntt  |   Nt1  |   Nll  |   Ntt  |   Ntl  |   Nlt  |   Nll  |   Ntt  |   Nt1  |   Nll  ||" << endl;
	OUT << "-------------------------------------------------------------------------------------------------------------" << endl;
	float nt2sum_mm(0.), nt10sum_mm(0.), nt0sum_mm(0.);
	float nt2sum_em(0.), nt10sum_em(0.), nt01sum_em(0.), nt0sum_em(0.);
	float nt2sum_ee(0.), nt10sum_ee(0.), nt0sum_ee(0.);

	const float mmTrigScale = 0.92;
	const float emTrigScale = 0.95;
	const float eeTrigScale = 1.00;

	for(size_t i = 0; i < fMCBG.size(); ++i){
		Sample *S = fSamples[fMCBG[i]];
		float scale = fLumiNorm / S->lumi;

		float temp_nt2_mm  = mmTrigScale*scale*S->region[reg][hilo].mm.nt20_pt->Integral(0, getNPt2Bins(Muon)+1); nt2sum_mm  += temp_nt2_mm ;
		float temp_nt10_mm = mmTrigScale*scale*S->region[reg][hilo].mm.nt10_pt->Integral(0, getNPt2Bins(Muon)+1); nt10sum_mm += temp_nt10_mm;
		float temp_nt0_mm  = mmTrigScale*scale*S->region[reg][hilo].mm.nt01_pt->Integral(0, getNPt2Bins(Muon)+1); nt0sum_mm  += temp_nt0_mm ;
		float temp_nt2_em  = emTrigScale*scale*S->region[reg][hilo].em.nt20_pt->Integral(0, getNPt2Bins(ElMu)+1); nt2sum_em  += temp_nt2_em ;
		float temp_nt10_em = emTrigScale*scale*S->region[reg][hilo].em.nt10_pt->Integral(0, getNPt2Bins(ElMu)+1); nt10sum_em += temp_nt10_em;
		float temp_nt01_em = emTrigScale*scale*S->region[reg][hilo].em.nt01_pt->Integral(0, getNPt2Bins(ElMu)+1); nt01sum_em += temp_nt01_em;
		float temp_nt0_em  = emTrigScale*scale*S->region[reg][hilo].em.nt00_pt->Integral(0, getNPt2Bins(ElMu)+1); nt0sum_em  += temp_nt0_em ;
		float temp_nt2_ee  = eeTrigScale*scale*S->region[reg][hilo].ee.nt20_pt->Integral(0, getNPt2Bins(Elec)+1); nt2sum_ee  += temp_nt2_ee ;
		float temp_nt10_ee = eeTrigScale*scale*S->region[reg][hilo].ee.nt10_pt->Integral(0, getNPt2Bins(Elec)+1); nt10sum_ee += temp_nt10_ee;
		float temp_nt0_ee  = eeTrigScale*scale*S->region[reg][hilo].ee.nt01_pt->Integral(0, getNPt2Bins(Elec)+1); nt0sum_ee  += temp_nt0_ee ;

		TString tempname = S->sname;
		OUT << Form("%16s & %6.2f & %6.2f & %6.2f & %6.2f & %6.2f & %6.2f & %6.2f & %6.2f & %6.2f & %6.2f \\\\\n", (tempname.ReplaceAll("_","\\_")).Data(),
		temp_nt2_mm , temp_nt10_mm, temp_nt0_mm,
		temp_nt2_em , temp_nt10_em, temp_nt01_em, temp_nt0_em,
		temp_nt2_ee , temp_nt10_ee, temp_nt0_ee);
	}	
	OUT << "\\hline" << endl;
	OUT << Form("%16s & %6.2f & %6.2f & %6.2f & %6.2f & %6.2f & %6.2f & %6.2f & %6.2f & %6.2f & %6.2f \\\\\n", "MC sum",
	nt2sum_mm ,	nt10sum_mm,	nt0sum_mm ,
	nt2sum_em ,	nt10sum_em,	nt01sum_em,	nt0sum_em ,
	nt2sum_ee ,	nt10sum_ee,	nt0sum_ee);
	OUT << "\\hline" << endl;

	for(gSample i = sample_begin; i < gNSAMPLES; i=gSample(i+1)){
		Sample *S = fSamples[i];
		if(S->datamc != 2) continue;
		float scale = fLumiNorm / S->lumi;

		float temp_nt2_mm  = mmTrigScale*scale*S->region[reg][hilo].mm.nt20_pt->Integral(0, getNPt2Bins(Muon)+1);
		float temp_nt10_mm = mmTrigScale*scale*S->region[reg][hilo].mm.nt10_pt->Integral(0, getNPt2Bins(Muon)+1);
		float temp_nt0_mm  = mmTrigScale*scale*S->region[reg][hilo].mm.nt01_pt->Integral(0, getNPt2Bins(Muon)+1);
		float temp_nt2_em  = emTrigScale*scale*S->region[reg][hilo].em.nt20_pt->Integral(0, getNPt2Bins(ElMu)+1);
		float temp_nt10_em = emTrigScale*scale*S->region[reg][hilo].em.nt10_pt->Integral(0, getNPt2Bins(ElMu)+1);
		float temp_nt01_em = emTrigScale*scale*S->region[reg][hilo].em.nt01_pt->Integral(0, getNPt2Bins(ElMu)+1);
		float temp_nt0_em  = emTrigScale*scale*S->region[reg][hilo].em.nt00_pt->Integral(0, getNPt2Bins(ElMu)+1);
		float temp_nt2_ee  = eeTrigScale*scale*S->region[reg][hilo].ee.nt20_pt->Integral(0, getNPt2Bins(Elec)+1);
		float temp_nt10_ee = eeTrigScale*scale*S->region[reg][hilo].ee.nt10_pt->Integral(0, getNPt2Bins(Elec)+1);
		float temp_nt0_ee  = eeTrigScale*scale*S->region[reg][hilo].ee.nt01_pt->Integral(0, getNPt2Bins(Elec)+1);

		TString tempname = S->sname;
		OUT << Form("%16s & %6.2f & %6.2f & %6.2f & %6.2f & %6.2f & %6.2f & %6.2f & %6.2f & %6.2f & %6.2f \\\\\n", (tempname.ReplaceAll("_","\\_")).Data(),
		temp_nt2_mm , temp_nt10_mm, temp_nt0_mm,
		temp_nt2_em , temp_nt10_em, temp_nt01_em, temp_nt0_em,
		temp_nt2_ee , temp_nt10_ee, temp_nt0_ee);
	}	
	OUT << "\\hline" << endl;
	OUT << Form("%16s & %6.0f & %6.0f & %6.0f & %6.0f & %6.0f & %6.0f & %6.0f & %6.0f & %6.0f & %6.0f \\\\\n", "Data",
	nt2_mm, nt10_mm, nt0_mm, nt2_em, nt10_em, nt01_em, nt0_em, nt2_ee, nt10_ee, nt0_ee);
	OUT << "\\hline" << endl;	
	OUT << endl;

	///////////////////////////////////////////////////////////////////////////////////
	// PREDICTIONS ////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////
	FakeRatios *FR = new FakeRatios();
	FR->setNToyMCs(100);
	FR->setAddESyst(0.5);
	// FR->setAddESyst(0.0);

	// FR->setMFRatio(mufratio_data, 0.10);
	// FR->setEFRatio(elfratio_data, 0.10);
	// FR->setMPRatio(mupratio_data, 0.05);
	// FR->setEPRatio(elpratio_data, 0.05);
	FR->setMFRatio(mufratio_data, mufratio_data_e); // set error to pure statistical of ratio
	FR->setEFRatio(elfratio_data, elfratio_data_e);
	FR->setMPRatio(mupratio_data, mupratio_data_e);
	FR->setEPRatio(elpratio_data, elpratio_data_e);

	FR->setMMNtl(nt2_mm, nt10_mm, nt0_mm);
	FR->setEENtl(nt2_ee, nt10_ee, nt0_ee);
	FR->setEMNtl(nt2_em, nt10_em, nt01_em, nt0_em);

	float nF_mm = 0;
	float nF_em = 0;
	float nF_ee = 0;
	float nSF   = 0;
	float nDF   = 0;
	float nF    = 0;

	// Event-by-event differential ratios:
	if(diffratio){
		nF_mm = npf_mm + nff_mm;
		nF_em = npf_em+nfp_em+nff_em;
		nF_ee = npf_ee+nff_ee;
		nSF   = npf_mm + npf_em + nfp_em + npf_ee;
		nDF   = nff_mm + nff_em + nff_ee;
		nF    = nF_mm + nF_em + nF_ee;
	}
	// Flat ratios:
	else{
		npp_mm = FR->getMMNpp();
		npf_mm = FR->getMMNpf();
		nff_mm = FR->getMMNff();
		npp_ee = FR->getEENpp();
		npf_ee = FR->getEENpf();
		nff_ee = FR->getEENff();
		npp_em = FR->getEMNpp();
		npf_em = FR->getEMNpf();
		nfp_em = FR->getEMNfp();
		nff_em = FR->getEMNff();
		nF_mm  = FR->getMMTotFakes();
		nF_em  = FR->getEMTotFakes();
		nF_ee  = FR->getEETotFakes();
		nSF    = FR->getTotSingleFakes();
		nDF    = FR->getTotDoubleFakes();
		nF     = FR->getTotFakes();		
	}

	OUT << "  Fake Predictions:" << endl;
	OUT << "------------------------------------------------------------------------------------------" << endl;
	OUT << "                 |          Mu/Mu        |         El/El         |          El/Mu        |" << endl;
	OUT << "------------------------------------------------------------------------------------------" << endl;
	OUT << " Npp             |" << Form(" %5.1f ± %5.1f ± %5.1f | %5.1f ± %5.1f ± %5.1f | %5.1f ± %5.1f ± %5.1f |",
	npp_mm, FR->getMMNppEStat(), FakeESyst*npp_mm,
	npp_ee, FR->getEENppEStat(), FakeESyst*npp_ee, 
	npp_em, FR->getEMNppEStat(), FakeESyst*npp_em) << endl;
	OUT << " Npf             |" << Form(" %5.1f ± %5.1f ± %5.1f | %5.1f ± %5.1f ± %5.1f | %5.1f ± %5.1f ± %5.1f |",
	npf_mm, FR->getMMNpfEStat(), FakeESyst*npf_mm,
	npf_ee, FR->getEENpfEStat(), FakeESyst*npf_ee, 
	npf_em, FR->getEMNpfEStat(), FakeESyst*npf_em) << endl;
	OUT << " Nfp             |" << Form("    -                  |    -                  | %5.1f ± %5.1f ± %5.1f |",
	nfp_em, FR->getEMNfpEStat(), FakeESyst*nfp_em) << endl;
	OUT << " Nff             |" << Form(" %5.1f ± %5.1f ± %5.1f | %5.1f ± %5.1f ± %5.1f | %5.1f ± %5.1f ± %5.1f |",
	nff_mm, FR->getMMNffEStat(), FakeESyst*nff_mm,
	nff_ee, FR->getEENffEStat(), FakeESyst*nff_ee, 
	nff_em, FR->getEMNffEStat(), FakeESyst*nff_em) << endl;
	OUT << "------------------------------------------------------------------------------------------" << endl;
	OUT << " Total Fakes     |" << Form(" %5.1f ± %5.1f ± %5.1f | %5.1f ± %5.1f ± %5.1f | %5.1f ± %5.1f ± %5.1f |",
	nF_mm, FR->getMMTotEStat(), FakeESyst*nF_mm,
	nF_ee, FR->getEETotEStat(), FakeESyst*nF_ee, 
	nF_em, FR->getEMTotEStat(), FakeESyst*nF_em) << endl;
	OUT << "------------------------------------------------------------------------------------------" << endl;
	OUT << " (Value ± E_stat ± E_syst) " << endl;
	OUT << "//////////////////////////////////////////////////////////////////////////////////////////" << endl;
	OUT << endl;

	///////////////////////////////////////////////////////////////////////////////////
	// E-CHARGE MISID /////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////
	float nt2_ee_chmid(0.), nt2_ee_chmid_e1(0.), nt2_ee_chmid_e2(0.);
	float nt2_em_chmid(0.), nt2_em_chmid_e1(0.), nt2_em_chmid_e2(0.);
	
	// Abbreviations
	float fb  = gEChMisIDB;
	float fbE = gEChMisIDB_E;
	float fe  = gEChMisIDE;
	float feE = gEChMisIDE_E;

	// Simple error propagation assuming error on number of events is sqrt(N)
	nt2_ee_chmid    = 2*fb*nt2_ee_BB_os + 2*fe*nt2_ee_EE_os + (fb+fe)*nt2_ee_EB_os;
	nt2_ee_chmid_e1 = sqrt( (4*fb*fb*FR->getEStat2(nt2_ee_BB_os)) + (4*fe*fe*FR->getEStat2(nt2_ee_EE_os)) + (fb+fe)*(fb+fe)*FR->getEStat2(nt2_ee_EB_os) ); // stat only
	nt2_ee_chmid_e2 = sqrt( (4*nt2_ee_BB_os*nt2_ee_BB_os*fbE*fbE) + (4*nt2_ee_EE_os*nt2_ee_EE_os*feE*feE) + (fbE*fbE+feE*feE)*nt2_ee_EB_os*nt2_ee_EB_os ); // syst only

	nt2_em_chmid    = fb*nt2_em_BB_os + fe*nt2_em_EE_os;
	nt2_em_chmid_e1 = sqrt( fb*fb*FR->getEStat2(nt2_em_BB_os) + fe*fe*FR->getEStat2(nt2_em_EE_os) );
	nt2_em_chmid_e2 = sqrt( nt2_em_BB_os*nt2_em_BB_os * fbE*fbE + nt2_em_EE_os*nt2_em_EE_os * feE*feE );

	///////////////////////////////////////////////////////////////////////////////////
	// PRINTOUT ///////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////
	OUT << "--------------------------------------------------------------" << endl;
	OUT << "       E-ChMisID  ||       Barrel       |       Endcap      ||" << endl;
	OUT << "--------------------------------------------------------------" << endl;
	OUT << "                  ||";
	OUT << setw(7)  << setprecision(2) << fb  << " +/- " << setw(7) << setprecision(3) << fbE  << " |";
	OUT << setw(7)  << setprecision(2) << fe  << " +/- " << setw(7) << setprecision(3) << feE  << " ||";
	OUT << endl;
	OUT << "--------------------------------------------------------------" << endl << endl;

	OUT << "-----------------------------------------------------------------------" << endl;
	OUT << "                 ||       E/Mu        ||             E/E             ||" << endl;
	OUT << "      OS-YIELDS  ||   N_B   |   N_E   ||   N_BB  |   N_EB  |   N_EE  ||" << endl;
	OUT << "-----------------------------------------------------------------------" << endl;

	float mc_os_em_bb_sum(0.), mc_os_em_ee_sum(0.);
	float mc_os_ee_bb_sum(0.), mc_os_ee_eb_sum(0.), mc_os_ee_ee_sum(0.);

	for(size_t i = 0; i < fMCBG.size(); ++i){
		Sample *S = fSamples[fMCBG[i]];
		float scale = fLumiNorm / S->lumi;
		
		mc_os_em_bb_sum += emTrigScale*scale*S->region[reg][hilo].em.nt20_OS_BB_pt->GetEntries();
		mc_os_em_ee_sum += emTrigScale*scale*S->region[reg][hilo].em.nt20_OS_EE_pt->GetEntries();
		mc_os_ee_bb_sum += eeTrigScale*scale*S->region[reg][hilo].ee.nt20_OS_BB_pt->GetEntries();
		mc_os_ee_eb_sum += eeTrigScale*scale*S->region[reg][hilo].ee.nt20_OS_EB_pt->GetEntries();
		mc_os_ee_ee_sum += eeTrigScale*scale*S->region[reg][hilo].ee.nt20_OS_EE_pt->GetEntries();

		OUT << setw(16) << S->sname << " || ";
		OUT << setw(7)  << setprecision(2) << emTrigScale*scale*S->region[reg][hilo].em.nt20_OS_BB_pt->GetEntries() << " | ";
		OUT << setw(7)  << setprecision(2) << emTrigScale*scale*S->region[reg][hilo].em.nt20_OS_EE_pt->GetEntries() << " || ";
		OUT << setw(7)  << setprecision(2) << eeTrigScale*scale*S->region[reg][hilo].ee.nt20_OS_BB_pt->GetEntries() << " | ";
		OUT << setw(7)  << setprecision(2) << eeTrigScale*scale*S->region[reg][hilo].ee.nt20_OS_EB_pt->GetEntries() << " | ";
		OUT << setw(7)  << setprecision(2) << eeTrigScale*scale*S->region[reg][hilo].ee.nt20_OS_EE_pt->GetEntries() << " || ";
		OUT << endl;
	}	
	OUT << "-----------------------------------------------------------------------" << endl;
	OUT << setw(16) << "MC sum" << " || ";
	OUT << setw(7) << Form("%5.1f", mc_os_em_bb_sum ) << " | ";
	OUT << setw(7) << Form("%5.1f", mc_os_em_ee_sum ) << " || ";
	OUT << setw(7) << Form("%5.1f", mc_os_ee_bb_sum ) << " | ";
	OUT << setw(7) << Form("%5.1f", mc_os_ee_eb_sum ) << " | ";
	OUT << setw(7) << Form("%5.1f", mc_os_ee_ee_sum ) << " || ";
	OUT << endl;
	OUT << "-----------------------------------------------------------------------" << endl;
	OUT << setw(16) << "data"  << " || ";
	OUT << setw(7) << Form("%5.0f", nt2_em_BB_os ) << " | ";
	OUT << setw(7) << Form("%5.0f", nt2_em_EE_os ) << " || ";
	OUT << setw(7) << Form("%5.0f", nt2_ee_BB_os ) << " | ";
	OUT << setw(7) << Form("%5.0f", nt2_ee_EB_os ) << " | ";
	OUT << setw(7) << Form("%5.0f", nt2_ee_EE_os ) << " || ";
	OUT << endl;
	OUT << "-----------------------------------------------------------------------" << endl;
	OUT << setw(16) << "pred. SS contr."  << " || ";
	OUT << setw(7) << Form("%6.4f",   fb   * nt2_em_BB_os ) << " | ";
	OUT << setw(7) << Form("%6.4f",   fe   * nt2_em_EE_os ) << " || ";
	OUT << setw(7) << Form("%6.4f", 2*fb   * nt2_ee_BB_os ) << " | ";
	OUT << setw(7) << Form("%6.4f", (fb+fe)* nt2_ee_EB_os ) << " | ";
	OUT << setw(7) << Form("%6.4f", 2*fe   * nt2_ee_EE_os ) << " || ";
	OUT << endl;
	OUT << "-----------------------------------------------------------------------" << endl << endl;
	OUT << "/////////////////////////////////////////////////////////////////////////////" << endl;
	OUT << "----------------------------------------------------------------------------------------------" << endl;
	OUT << "       SUMMARY   ||         Mu/Mu         ||         E/Mu          ||          E/E          ||" << endl;
	OUT << "==============================================================================================" << endl;
	OUT << Form("%16s || %5.2f ± %5.2f ± %5.2f || %5.2f ± %5.2f ± %5.2f || %5.2f ± %5.2f ± %5.2f ||\n", "pred. fakes",
	nF_mm, FR->getMMTotEStat(), FakeESyst*nF_mm,
	nF_em, FR->getEMTotEStat(), FakeESyst*nF_em,
	nF_ee, FR->getEETotEStat(), FakeESyst*nF_ee);
	OUT << Form("%16s ||                       || %5.2f ± %5.2f ± %5.2f || %5.2f ± %5.2f ± %5.2f ||\n", "pred. chmisid",
	nt2_em_chmid, nt2_em_chmid_e1, nt2_em_chmid_e2, nt2_ee_chmid, nt2_ee_chmid_e1, nt2_ee_chmid_e2);

	OUT << "----------------------------------------------------------------------------------------------" << endl;
	float nt2_rare_mc_mm(0.),    nt2_rare_mc_em(0.),    nt2_rare_mc_ee(0.);
	float nt2_rare_mc_mm_e1(0.), nt2_rare_mc_em_e1(0.), nt2_rare_mc_ee_e1(0.);
	for(size_t i = 0; i < fMCRareSM.size(); ++i){
		Sample *S = fSamples[fMCRareSM[i]];
		float scale = fLumiNorm/S->lumi;

		float temp_nt2_mm = mmTrigScale*scale*S->region[reg][hilo].mm.nt20_pt->Integral(0, getNPt2Bins(Muon)+1);
		float temp_nt2_em = emTrigScale*scale*S->region[reg][hilo].em.nt20_pt->Integral(0, getNPt2Bins(ElMu)+1);
		float temp_nt2_ee = eeTrigScale*scale*S->region[reg][hilo].ee.nt20_pt->Integral(0, getNPt2Bins(Elec)+1);

		nt2_rare_mc_mm += temp_nt2_mm;
		nt2_rare_mc_em += temp_nt2_em;
		nt2_rare_mc_ee += temp_nt2_ee;

		nt2_rare_mc_mm_e1 += mmTrigScale*mmTrigScale*scale*scale * FR->getEStat2(S->region[reg][hilo].mm.nt20_pt->GetEntries()); // for stat error take actual entries, not pileup weighted integral...
		nt2_rare_mc_em_e1 += emTrigScale*emTrigScale*scale*scale * FR->getEStat2(S->region[reg][hilo].em.nt20_pt->GetEntries());
		nt2_rare_mc_ee_e1 += eeTrigScale*eeTrigScale*scale*scale * FR->getEStat2(S->region[reg][hilo].ee.nt20_pt->GetEntries());

		OUT << Form("%16s || %5.2f ± %5.2f         || %5.2f ± %5.2f         || %5.2f ± %5.2f         ||\n", S->sname.Data(),
		temp_nt2_mm, mmTrigScale*scale*FR->getEStat(S->region[reg][hilo].mm.nt20_pt->GetEntries()),
		temp_nt2_em, emTrigScale*scale*FR->getEStat(S->region[reg][hilo].em.nt20_pt->GetEntries()),
		temp_nt2_ee, eeTrigScale*scale*FR->getEStat(S->region[reg][hilo].ee.nt20_pt->GetEntries()));
	}
	OUT << "----------------------------------------------------------------------------------------------" << endl;
	OUT << Form("%16s || %5.2f ± %5.2f ± %5.2f || %5.2f ± %5.2f ± %5.2f || %5.2f ± %5.2f ± %5.2f ||\n", "Rare SM (Sum)",
	nt2_rare_mc_mm, sqrt(nt2_rare_mc_mm_e1), RareESyst*nt2_rare_mc_mm,
	nt2_rare_mc_em, sqrt(nt2_rare_mc_em_e1), RareESyst*nt2_rare_mc_em,
	nt2_rare_mc_ee, sqrt(nt2_rare_mc_ee_e1), RareESyst*nt2_rare_mc_ee);
	OUT << "----------------------------------------------------------------------------------------------" << endl;
	// Just add different errors in quadrature (they are independent)
	float mm_tot_sqerr1 = FR->getMMTotEStat()*FR->getMMTotEStat() + nt2_rare_mc_mm_e1;
	float em_tot_sqerr1 = FR->getEMTotEStat()*FR->getEMTotEStat() + nt2_em_chmid_e1*nt2_em_chmid_e1 + nt2_rare_mc_em_e1;
	float ee_tot_sqerr1 = FR->getEETotEStat()*FR->getEETotEStat() + nt2_ee_chmid_e1*nt2_ee_chmid_e1 + nt2_rare_mc_ee_e1;
	float mm_tot_sqerr2 = nF_mm*nF_mm*FakeESyst2 + RareESyst2*nt2_rare_mc_mm*nt2_rare_mc_mm;
	float em_tot_sqerr2 = nF_em*nF_em*FakeESyst2 + nt2_em_chmid_e2*nt2_em_chmid_e2 + RareESyst2*nt2_rare_mc_em*nt2_rare_mc_em;
	float ee_tot_sqerr2 = nF_ee*nF_ee*FakeESyst2 + nt2_ee_chmid_e2*nt2_ee_chmid_e2 + RareESyst2*nt2_rare_mc_ee*nt2_rare_mc_ee;
	OUT << Form("%16s || %5.2f ± %5.2f ± %5.2f || %5.2f ± %5.2f ± %5.2f || %5.2f ± %5.2f ± %5.2f ||\n", "tot. backgr.",
	nF_mm + nt2_rare_mc_mm,                sqrt(mm_tot_sqerr1), sqrt(mm_tot_sqerr2),
	nF_em + nt2_em_chmid + nt2_rare_mc_em, sqrt(em_tot_sqerr1), sqrt(em_tot_sqerr2),
	nF_ee + nt2_ee_chmid + nt2_rare_mc_ee, sqrt(ee_tot_sqerr1), sqrt(ee_tot_sqerr2));
	OUT << Form("%16s || %5.2f ± %5.2f         || %5.2f ± %5.2f         || %5.2f ± %5.2f         ||\n", "",
	nF_mm + nt2_rare_mc_mm,                sqrt(mm_tot_sqerr1 + mm_tot_sqerr2),
	nF_em + nt2_em_chmid + nt2_rare_mc_em, sqrt(em_tot_sqerr1 + em_tot_sqerr2),
	nF_ee + nt2_ee_chmid + nt2_rare_mc_ee, sqrt(ee_tot_sqerr1 + ee_tot_sqerr2));
	OUT << "==============================================================================================" << endl;
	OUT << Form("%16s || %5.0f                 || %5.0f                 || %5.0f                 ||\n", "observed", nt2_mm, nt2_em, nt2_ee);
	OUT << "----------------------------------------------------------------------------------------------" << endl;
	OUT << Form("%16s || %5.2f                 || %5.2f                 || %5.2f                 ||\n", "MC sum", nt2sum_mm, nt2sum_em, nt2sum_ee);
	OUT << "==============================================================================================" << endl;
	OUT << "combined observed: ";
	OUT << setw(5) << left << Form("%2.0f", nt2_mm+nt2_em+nt2_ee ) << endl;
	OUT << "        predicted: ";
	float tot_pred        = nF + nt2_rare_mc_mm + nt2_em_chmid + nt2_rare_mc_em + nt2_ee_chmid + nt2_rare_mc_ee;
	float comb_tot_sqerr1 = FR->getTotEStat()*FR->getTotEStat() + nt2_rare_mc_mm_e1 + nt2_em_chmid_e1*nt2_em_chmid_e1 + nt2_rare_mc_em_e1 + nt2_ee_chmid_e1*nt2_ee_chmid_e1 + nt2_rare_mc_ee_e1;
	float comb_tot_sqerr2 = nF*nF*FakeESyst2 + RareESyst2*(nt2_rare_mc_mm + nt2_rare_mc_em + nt2_rare_mc_ee)*(nt2_rare_mc_mm + nt2_rare_mc_em + nt2_rare_mc_ee) + nt2_em_chmid_e2*nt2_em_chmid_e2 + nt2_ee_chmid_e2*nt2_ee_chmid_e2;
	OUT << setw(5) << left << Form("%5.2f", tot_pred ) << " ± ";
	OUT << setw(5) << Form("%5.2f", sqrt(comb_tot_sqerr1)) << " ± ";
	OUT << setw(5) << Form("%5.2f", sqrt(comb_tot_sqerr2)) << endl;
	OUT << "==============================================================================================" << endl;
	OUT.close();
	
	///////////////////////////////////////////////////////////////////////////////////
	//  OUTPUT FOR PAS TABLE  /////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////
	fOUTSTREAM << Region::sname[reg] << endl;
	fOUTSTREAM << Form("{\\bf predicted BG b} & {\\boldmath $%5.2f\\pm %5.2f$} & {\\boldmath $%5.2f \\pm %5.2f$} & {\\boldmath $%5.2f\\pm %5.2f$} & {\\boldmath $%5.2f\\pm %5.2f$} & \\\\ \n",
	nF_ee + nt2_ee_chmid + nt2_rare_mc_ee, sqrt(ee_tot_sqerr1 + ee_tot_sqerr2),
	nF_mm + nt2_rare_mc_mm,                sqrt(mm_tot_sqerr1 + mm_tot_sqerr2),
	nF_em + nt2_em_chmid + nt2_rare_mc_em, sqrt(em_tot_sqerr1 + em_tot_sqerr2),
	tot_pred, sqrt(comb_tot_sqerr1 + comb_tot_sqerr2));
	fOUTSTREAM << Form("{\\bf observed} & {\\bf %2.0f} & {\\bf %2.0f} & {\\bf %2.0f} & {\\bf %2.0f}  & {\\bf XX} \\\\ \\hline \n", nt2_ee, nt2_mm, nt2_em, nt2_ee+nt2_mm+nt2_em);
	
	///////////////////////////////////////////////////////////////////////////////////
	//  OUTPUT FOR DIDARS PLOT  ///////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////
	// fOUTSTREAM2 << "// " + Region::sname[reg] << endl;
	// fOUTSTREAM2 << Form("float %s_SS_ETH[4]    = {%6.3f, %6.3f, %6.3f, %6.3f }; \n", Region::sname[reg].Data(), nt2_rare_mc_ee, nt2_rare_mc_mm, nt2_rare_mc_em, nt2_rare_mc_ee + nt2_rare_mc_mm + nt2_rare_mc_em);
	// fOUTSTREAM2 << Form("float %s_OS_ETH[4]    = {%6.3f, %6.3f, %6.3f, %6.3f }; \n", Region::sname[reg].Data(), nt2_ee_chmid, 0.00, nt2_em_chmid, nt2_ee_chmid + nt2_em_chmid);
	// fOUTSTREAM2 << Form("float %s_PF_ETH[4]    = {%6.3f, %6.3f, %6.3f, %6.3f }; \n", Region::sname[reg].Data(), FR->getEENpf(), FR->getMMNpf(), FR->getEMNpf()+FR->getEMNfp(), FR->getEENpf() + FR->getMMNpf() + FR->getEMNpf() + FR->getEMNfp());
	// fOUTSTREAM2 << Form("float %s_FF_ETH[4]    = {%6.3f, %6.3f, %6.3f, %6.3f }; \n", Region::sname[reg].Data(), FR->getEENff(), FR->getMMNff(), FR->getEMNff(), FR->getEENff()+FR->getMMNff()+FR->getEMNff());
	// fOUTSTREAM2 << Form("float %s_Error_ETH[4] = {%6.3f, %6.3f, %6.3f, %6.3f }; \n", Region::sname[reg].Data(), sqrt(ee_tot_sqerr1 + ee_tot_sqerr2), sqrt(mm_tot_sqerr1 + mm_tot_sqerr2), sqrt(em_tot_sqerr1 + em_tot_sqerr2), sqrt(comb_tot_sqerr1 + comb_tot_sqerr2));
	// fOUTSTREAM2 << endl;
	
	///////////////////////////////////////////////////////////////////////////////////
	//  OUTPUT FOR ANALYSIS NOTE  /////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////
	fOUTSTREAM3 << "%% " + Region::sname[reg] << endl;
	fOUTSTREAM3 << "-----------------------------------------------------------------" << endl;
	fOUTSTREAM3 << Form("DF:  %6.1f ± %6.1f  ( %5.1f±%5.1f | %5.1f±%5.1f | %5.1f±%5.1f )\n",
	nff_em + nff_mm + nff_ee, sqrt(FR->getTotDoubleEStat()*FR->getTotDoubleEStat() + nDF*nDF*FakeESyst2),
	nff_ee, sqrt(FR->getEENffEStat()*FR->getEENffEStat()+nff_ee*nff_ee*FakeESyst2),
	nff_mm, sqrt(FR->getMMNffEStat()*FR->getMMNffEStat()+nff_mm*nff_mm*FakeESyst2),
	nff_em, sqrt(FR->getEMNffEStat()*FR->getEMNffEStat()+nff_em*nff_em*FakeESyst2));
	fOUTSTREAM3 << Form("SF:  %6.1f ± %6.1f  ( %5.1f±%5.1f | %5.1f±%5.1f | %5.1f±%5.1f )\n",
	npf_em + nfp_em + npf_mm + npf_ee, sqrt(FR->getTotSingleEStat()*FR->getTotSingleEStat() + nSF*nSF*FakeESyst2),
	npf_ee,          sqrt(FR->getEENpfEStat()   *FR->getEENpfEStat()    +  npf_ee*npf_ee*FakeESyst2),
	npf_mm,          sqrt(FR->getMMNpfEStat()   *FR->getMMNpfEStat()    +  npf_mm*npf_mm*FakeESyst2),
	npf_em + nfp_em, sqrt(FR->getEMSingleEStat()*FR->getEMSingleEStat() + (npf_em+nfp_em)*(npf_em+nfp_em)*FakeESyst2));
	fOUTSTREAM3 << Form("CM:  %6.1f ± %6.1f  ( %5.1f±%5.1f |   -         | %5.1f±%5.1f )\n",
	nt2_ee_chmid + nt2_em_chmid, sqrt(nt2_ee_chmid_e1*nt2_ee_chmid_e1 + nt2_ee_chmid_e2*nt2_ee_chmid_e2 + nt2_em_chmid_e1*nt2_em_chmid_e1 + nt2_em_chmid_e2*nt2_em_chmid_e2),
	nt2_ee_chmid, sqrt(nt2_ee_chmid_e1*nt2_ee_chmid_e1 + nt2_ee_chmid_e2*nt2_ee_chmid_e2),
	nt2_em_chmid, sqrt(nt2_em_chmid_e1*nt2_em_chmid_e1 + nt2_em_chmid_e2*nt2_em_chmid_e2));
	fOUTSTREAM3 << Form("MC:  %6.1f ± %6.1f  ( %5.1f±%5.1f | %5.1f±%5.1f | %5.1f±%5.1f )\n",
	nt2_rare_mc_ee + nt2_rare_mc_mm + nt2_rare_mc_em, sqrt(nt2_rare_mc_ee_e1 + nt2_rare_mc_mm_e1 + nt2_rare_mc_em_e1 + RareESyst2*(nt2_rare_mc_ee + nt2_rare_mc_mm + nt2_rare_mc_em)*(nt2_rare_mc_ee + nt2_rare_mc_mm + nt2_rare_mc_em)),
	nt2_rare_mc_ee, sqrt(nt2_rare_mc_ee_e1 + RareESyst2*nt2_rare_mc_ee*nt2_rare_mc_ee),
	nt2_rare_mc_mm, sqrt(nt2_rare_mc_mm_e1 + RareESyst2*nt2_rare_mc_mm*nt2_rare_mc_mm),
	nt2_rare_mc_em, sqrt(nt2_rare_mc_em_e1 + RareESyst2*nt2_rare_mc_em*nt2_rare_mc_em));
	fOUTSTREAM3 << Form("Tot: %6.1f ± %6.1f  ( %5.1f±%5.1f | %5.1f±%5.1f | %5.1f±%5.1f )\n",
	tot_pred, sqrt(comb_tot_sqerr1 + comb_tot_sqerr2),
	nF_ee + nt2_rare_mc_ee + nt2_ee_chmid, sqrt(ee_tot_sqerr1 + ee_tot_sqerr2),
	nF_mm + nt2_rare_mc_mm, sqrt(mm_tot_sqerr1 + mm_tot_sqerr2),
	nF_em + nt2_rare_mc_em + nt2_em_chmid, sqrt(em_tot_sqerr1 + em_tot_sqerr2));
	fOUTSTREAM3 << Form("Obs: %4.0f             ( %3.0f         | %3.0f         | %3.0f         )\n", nt2_mm+nt2_em+nt2_ee, nt2_ee, nt2_mm, nt2_em);
	fOUTSTREAM3 << "-----------------------------------------------------------------" << endl;
	fOUTSTREAM3 << endl;
	
	///////////////////////////////////////////////////////////////////////////////////
	//  OUTPUT AS PLOT  ///////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////
	TH1D    *h_obs        = new TH1D("h_observed",   "Observed number of events",  3, 0., 3.);
	TH1D    *h_pred_sfake = new TH1D("h_pred_sfake", "Predicted single fakes", 3, 0., 3.);
	TH1D    *h_pred_dfake = new TH1D("h_pred_dfake", "Predicted double fakes", 3, 0., 3.);
	TH1D    *h_pred_chmid = new TH1D("h_pred_chmid", "Predicted charge mis id", 3, 0., 3.);
	TH1D    *h_pred_mc    = new TH1D("h_pred_mc",    "Predicted WW/WZ/ZZ", 3, 0., 3.);
	TH1D    *h_pred_tot   = new TH1D("h_pred_tot",   "Total Prediction", 3, 0., 3.);
	THStack *hs_pred      = new THStack("hs_predicted", "Predicted number of events");
	
	h_obs->SetMarkerColor(kBlack);
	h_obs->SetMarkerStyle(20);
	h_obs->SetMarkerSize(2.5);
	h_obs->SetLineWidth(2);
	h_obs->SetLineColor(kBlack);
	h_obs->SetFillColor(kBlack);
	
	h_pred_sfake->SetLineWidth(1);
	h_pred_dfake->SetLineWidth(1);
	h_pred_chmid->SetLineWidth(1);
	h_pred_mc   ->SetLineWidth(1);
	h_pred_sfake->SetLineColor(50);
	h_pred_sfake->SetFillColor(50);
	h_pred_dfake->SetLineColor(38);
	h_pred_dfake->SetFillColor(38);
	h_pred_chmid->SetLineColor(42);
	h_pred_chmid->SetFillColor(42);
	h_pred_mc   ->SetLineColor(31);
	h_pred_mc   ->SetFillColor(31);
	h_pred_tot  ->SetLineWidth(1);
	// h_pred_tot  ->SetFillColor(kBlack);
	// h_pred_tot  ->SetFillStyle(3013);
	h_pred_tot  ->SetFillColor(12);
	h_pred_tot  ->SetFillStyle(3005);
	
	// Add numbers:
	h_obs->SetBinContent(1, nt2_ee);
	h_obs->SetBinContent(2, nt2_mm);
	h_obs->SetBinContent(3, nt2_em);
	h_obs->SetBinError(1, FR->getEStat(nt2_ee)); // FIXME
	h_obs->SetBinError(2, FR->getEStat(nt2_mm)); // FIXME
	h_obs->SetBinError(3, FR->getEStat(nt2_em)); // FIXME
	
	h_pred_sfake->SetBinContent(1, npf_ee);
	h_pred_sfake->SetBinContent(2, npf_mm);
	h_pred_sfake->SetBinContent(3, npf_em+nfp_em);
	h_pred_sfake->GetXaxis()->SetBinLabel(1, "ee");
	h_pred_sfake->GetXaxis()->SetBinLabel(2, "#mu#mu");
	h_pred_sfake->GetXaxis()->SetBinLabel(3, "e#mu");
	
	h_pred_dfake->SetBinContent(1, nff_ee);
	h_pred_dfake->SetBinContent(2, nff_mm);
	h_pred_dfake->SetBinContent(3, nff_em);
	
	h_pred_chmid->SetBinContent(1, nt2_ee_chmid);
	h_pred_chmid->SetBinContent(2, 0.);
	h_pred_chmid->SetBinContent(3, nt2_em_chmid);
	
	h_pred_mc->SetBinContent(1, nt2_rare_mc_ee);
	h_pred_mc->SetBinContent(2, nt2_rare_mc_mm);
	h_pred_mc->SetBinContent(3, nt2_rare_mc_em);
	
	h_pred_tot->Add(h_pred_sfake);
	h_pred_tot->Add(h_pred_dfake);
	h_pred_tot->Add(h_pred_chmid);
	h_pred_tot->Add(h_pred_mc);
	h_pred_tot->SetBinError(1, sqrt(ee_tot_sqerr1 + ee_tot_sqerr2));
	h_pred_tot->SetBinError(2, sqrt(mm_tot_sqerr1 + mm_tot_sqerr2));
	h_pred_tot->SetBinError(3, sqrt(em_tot_sqerr1 + em_tot_sqerr2));
	
	hs_pred->Add(h_pred_sfake);
	hs_pred->Add(h_pred_dfake);
	hs_pred->Add(h_pred_chmid);
	hs_pred->Add(h_pred_mc);
	
	double max = h_obs->Integral();
	h_obs       ->SetMaximum(max>1?max+1:1.);
	h_pred_sfake->SetMaximum(max>1?max+1:1.);
	h_pred_dfake->SetMaximum(max>1?max+1:1.);
	h_pred_chmid->SetMaximum(max>1?max+1:1.);
	h_pred_mc   ->SetMaximum(max>1?max+1:1.);
	h_pred_tot  ->SetMaximum(max>1?max+1:1.);
	hs_pred     ->SetMaximum(max>1?max+1:1.);
	
	hs_pred->Draw("goff");
	hs_pred->GetXaxis()->SetBinLabel(1, "ee");
	hs_pred->GetXaxis()->SetBinLabel(2, "#mu#mu");
	hs_pred->GetXaxis()->SetBinLabel(3, "e#mu");
	hs_pred->GetXaxis()->SetLabelOffset(0.01);
	hs_pred->GetXaxis()->SetLabelFont(42);
	hs_pred->GetXaxis()->SetLabelSize(0.1);
	
	TLegend *leg = new TLegend(0.15,0.65,0.50,0.88);
	leg->AddEntry(h_obs,        "Observed","p");
	leg->AddEntry(h_pred_sfake, "Single Fakes","f");
	leg->AddEntry(h_pred_dfake, "Double Fakes","f");
	leg->AddEntry(h_pred_chmid, "Charge MisID","f");
	leg->AddEntry(h_pred_mc,    "Irreducible (MC)","f");
	leg->AddEntry(h_pred_tot,    "Total Uncertainty","f");
	leg->SetFillStyle(0);
	leg->SetTextFont(42);
	// leg->SetTextSize(0.05);
	leg->SetBorderSize(0);
	
	TCanvas *c_temp = new TCanvas("C_ObsPred", "Observed vs Predicted", 0, 0, 600, 600);
	c_temp->cd();
	
	hs_pred->Draw("hist");
	h_pred_tot->DrawCopy("0 E2 same");
	h_obs->DrawCopy("PE X0 same");
	leg->Draw();
	
	lat->SetTextSize(0.03);
	lat->DrawLatex(0.16,0.60, Form("H_{T} > %.0f GeV, N_{Jets} #geq %1d", Region::minHT[reg], Region::minNjets[reg]));
	if(reg != Control) lat->DrawLatex(0.16,0.55, Form("E_{T}^{miss} > %.0f GeV", Region::minMet[reg]));
	if(reg == Control) lat->DrawLatex(0.16,0.55, Form("E_{T}^{miss} > %.0f GeV, < %.0f GeV", Region::minMet[reg], Region::maxMet[reg]));
	drawTopLine();
	
	gPad->RedrawAxis();
	// Util::PrintNoEPS(c_temp, "ObsPred_" + Region::sname[reg], fOutputDir + fOutputSubDir, NULL);
	Util::PrintPDF(c_temp,   "ObsPred_" + Region::sname[reg], fOutputDir + fOutputSubDir);
	delete c_temp;	
	delete h_obs, h_pred_sfake, h_pred_dfake, h_pred_chmid, h_pred_mc, h_pred_tot, hs_pred;
	delete FR;
}

void SSDLPlotter::makeDiffPrediction(){
	fOutputSubDir = "DiffPredictionPlots/";
	TLatex *lat = new TLatex();
	lat->SetNDC(kTRUE);
	lat->SetTextColor(kBlack);
	lat->SetTextSize(0.04);

	vector<int> musamples;
	vector<int> elsamples;
	vector<int> emusamples;

	musamples = fMuData;
	elsamples = fEGData;
	emusamples = fMuEGData;

	///////////////////////////////////////////////////////////////////////////////////
	// RATIOS /////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////
	float mufratio_data(0.),  mufratio_data_e(0.);
	float mupratio_data(0.),  mupratio_data_e(0.);
	float elfratio_data(0.),  elfratio_data_e(0.);
	float elpratio_data(0.),  elpratio_data_e(0.);

	calculateRatio(fMuData, Muon, SigSup, mufratio_data, mufratio_data_e);
	calculateRatio(fMuData, Muon, ZDecay, mupratio_data, mupratio_data_e);

	calculateRatio(fEGData, Elec, SigSup, elfratio_data, elfratio_data_e);
	calculateRatio(fEGData, Elec, ZDecay, elpratio_data, elpratio_data_e);

	//{"HT1", "HT2", "MET1", "MET2", "NJets", "MT2", "PT1", "PT2", "NBJets"};
	float binwidthscale[gNDiffVars] = {100., 100., 30., 30., 1., 25., 20., 10., 1.};

	// Loop on the different variables
	for(size_t j = 0; j < gNDiffVars; ++j){
		TString varname    = DiffPredYields::var_name[j];
		const int nbins    = DiffPredYields::nbins[j];
		const double *bins = DiffPredYields::bins[j];

		///////////////////////////////////////////////////////////////////////////////////
		// OBSERVATIONS ///////////////////////////////////////////////////////////////////
		///////////////////////////////////////////////////////////////////////////////////
		TH1D *nt11 = new TH1D(Form("NT11_%s", varname.Data()), varname, nbins, bins); nt11->Sumw2();

		TH1D *nt11_mm = new TH1D(Form("NT11_MM_%s", varname.Data()), varname, nbins, bins); nt11_mm->Sumw2();
		TH1D *nt10_mm = new TH1D(Form("NT10_MM_%s", varname.Data()), varname, nbins, bins); nt10_mm->Sumw2();
		TH1D *nt01_mm = new TH1D(Form("NT01_MM_%s", varname.Data()), varname, nbins, bins); nt01_mm->Sumw2();
		TH1D *nt00_mm = new TH1D(Form("NT00_MM_%s", varname.Data()), varname, nbins, bins); nt00_mm->Sumw2();
		TH1D *nt11_ee = new TH1D(Form("NT11_EE_%s", varname.Data()), varname, nbins, bins); nt11_ee->Sumw2();
		TH1D *nt10_ee = new TH1D(Form("NT10_EE_%s", varname.Data()), varname, nbins, bins); nt10_ee->Sumw2();
		TH1D *nt01_ee = new TH1D(Form("NT01_EE_%s", varname.Data()), varname, nbins, bins); nt01_ee->Sumw2();
		TH1D *nt00_ee = new TH1D(Form("NT00_EE_%s", varname.Data()), varname, nbins, bins); nt00_ee->Sumw2();
		TH1D *nt11_em = new TH1D(Form("NT11_EM_%s", varname.Data()), varname, nbins, bins); nt11_em->Sumw2();
		TH1D *nt10_em = new TH1D(Form("NT10_EM_%s", varname.Data()), varname, nbins, bins); nt10_em->Sumw2();
		TH1D *nt01_em = new TH1D(Form("NT01_EM_%s", varname.Data()), varname, nbins, bins); nt01_em->Sumw2();
		TH1D *nt00_em = new TH1D(Form("NT00_EM_%s", varname.Data()), varname, nbins, bins); nt00_em->Sumw2();

		// OS yields
		TH1D *nt2_os_ee_bb = new TH1D(Form("NT2_OS_EE_BB_%s", varname.Data()), varname, nbins, bins); nt2_os_ee_bb->Sumw2();
		TH1D *nt2_os_ee_eb = new TH1D(Form("NT2_OS_EE_EB_%s", varname.Data()), varname, nbins, bins); nt2_os_ee_eb->Sumw2();
		TH1D *nt2_os_ee_ee = new TH1D(Form("NT2_OS_EE_EE_%s", varname.Data()), varname, nbins, bins); nt2_os_ee_ee->Sumw2();
		TH1D *nt2_os_em_bb = new TH1D(Form("NT2_OS_EM_BB_%s", varname.Data()), varname, nbins, bins); nt2_os_em_bb->Sumw2();
		TH1D *nt2_os_em_ee = new TH1D(Form("NT2_OS_EM_EE_%s", varname.Data()), varname, nbins, bins); nt2_os_em_ee->Sumw2();

		for(size_t i = 0; i < musamples.size(); ++i){
			Sample *S = fSamples[musamples[i]];
			nt11_mm->Add(S->diffyields[Muon].hnt11[j]);
			nt10_mm->Add(S->diffyields[Muon].hnt10[j]);
			nt01_mm->Add(S->diffyields[Muon].hnt01[j]);
			nt00_mm->Add(S->diffyields[Muon].hnt00[j]);
		}
		for(size_t i = 0; i < elsamples.size(); ++i){
			Sample *S = fSamples[elsamples[i]];
			nt11_ee->Add(S->diffyields[Elec].hnt11[j]);
			nt10_ee->Add(S->diffyields[Elec].hnt10[j]);
			nt01_ee->Add(S->diffyields[Elec].hnt01[j]);
			nt00_ee->Add(S->diffyields[Elec].hnt00[j]);

			nt2_os_ee_bb->Add(S->diffyields[Elec].hnt2_os_BB[j]);
			nt2_os_ee_eb->Add(S->diffyields[Elec].hnt2_os_EB[j]);
			nt2_os_ee_ee->Add(S->diffyields[Elec].hnt2_os_EE[j]);
		}
		for(size_t i = 0; i < emusamples.size(); ++i){
			Sample *S = fSamples[emusamples[i]];
			nt11_em->Add(S->diffyields[ElMu].hnt11[j]);
			nt10_em->Add(S->diffyields[ElMu].hnt10[j]);
			nt01_em->Add(S->diffyields[ElMu].hnt01[j]);
			nt00_em->Add(S->diffyields[ElMu].hnt00[j]);

			nt2_os_em_bb->Add(S->diffyields[ElMu].hnt2_os_BB[j]);
			nt2_os_em_ee->Add(S->diffyields[ElMu].hnt2_os_EE[j]);
		}
		
		nt11->Add(nt11_mm);
		nt11->Add(nt11_ee);
		nt11->Add(nt11_em);

		
		///////////////////////////////////////////////////////////////////////////////////
		// Errors /////////////////////////////////////////////////////////////////////////
		///////////////////////////////////////////////////////////////////////////////////
		// Save the squared sum of total error in the bin errors of these histos,
		// then apply square root before plotting
		TH1D *totbg    = new TH1D(Form("TotBG_%s",    varname.Data()), varname, nbins, bins); totbg   ->Sumw2();
		TH1D *totbg_mm = new TH1D(Form("TotBG_mm_%s", varname.Data()), varname, nbins, bins); totbg_mm->Sumw2();
		TH1D *totbg_em = new TH1D(Form("TotBG_em_%s", varname.Data()), varname, nbins, bins); totbg_em->Sumw2();
		TH1D *totbg_ee = new TH1D(Form("TotBG_ee_%s", varname.Data()), varname, nbins, bins); totbg_ee->Sumw2();

		///////////////////////////////////////////////////////////////////////////////////
		// MC Predictions /////////////////////////////////////////////////////////////////
		///////////////////////////////////////////////////////////////////////////////////
		// TH1D *nt11_mc = new TH1D(Form("NT11_MC_%s", varname.Data()), varname, nbins, bins); nt11_mc->Sumw2();
		// 
		// TH1D *nt11_mm_mc = new TH1D(Form("NT11_MM_MC_%s", varname.Data()), varname, nbins, bins); nt11_mm_mc->Sumw2();
		// TH1D *nt11_ee_mc = new TH1D(Form("NT11_EE_MC_%s", varname.Data()), varname, nbins, bins); nt11_ee_mc->Sumw2();
		// TH1D *nt11_em_mc = new TH1D(Form("NT11_EM_MC_%s", varname.Data()), varname, nbins, bins); nt11_em_mc->Sumw2();
		// 
		// for(size_t i = 0; i < fMCBG.size(); ++i){
		// 	Sample *S = fSamples[fMCBG[i]];
		// 	float scale = fLumiNorm / S->lumi;
		// 	nt11_mm_mc->Add(S->diffyields[Muon].hnt11[j], scale);
		// 	nt11_ee_mc->Add(S->diffyields[Elec].hnt11[j], scale);
		// 	nt11_em_mc->Add(S->diffyields[ElMu].hnt11[j], scale);
		// }
		// 
		// nt11_mc->Add(nt11_mm_mc);
		// nt11_mc->Add(nt11_ee_mc);
		// nt11_mc->Add(nt11_em_mc);

		///////////////////////////////////////////////////////////////////////////////////
		// RARE SM MC /////////////////////////////////////////////////////////////////////
		///////////////////////////////////////////////////////////////////////////////////
		TH1D *nt11_ss = new TH1D(Form("NT11_SS_%s", varname.Data()), varname, nbins, bins); nt11_ss->Sumw2();

		TH1D *nt11_mm_ss = new TH1D(Form("NT11_MM_SS_%s", varname.Data()), varname, nbins, bins); nt11_mm_ss->Sumw2();
		TH1D *nt11_ee_ss = new TH1D(Form("NT11_EE_SS_%s", varname.Data()), varname, nbins, bins); nt11_ee_ss->Sumw2();
		TH1D *nt11_em_ss = new TH1D(Form("NT11_EM_SS_%s", varname.Data()), varname, nbins, bins); nt11_em_ss->Sumw2();

		for(size_t i = 0; i < fMCRareSM.size(); ++i){
			FakeRatios *FR = new FakeRatios();
			Sample *S = fSamples[fMCRareSM[i]];
			float scale = fLumiNorm / S->lumi;
			nt11_mm_ss->Add(S->diffyields[Muon].hnt11[j], scale);
			nt11_ee_ss->Add(S->diffyields[Elec].hnt11[j], scale);
			nt11_em_ss->Add(S->diffyields[ElMu].hnt11[j], scale);

			// Errors
			for(size_t b = 0; b < nbins; ++b){
				float ss_mm = S->diffyields[Muon].hnt11[j]->GetBinContent(b+1);
				float ss_ee = S->diffyields[Elec].hnt11[j]->GetBinContent(b+1);
				float ss_em = S->diffyields[ElMu].hnt11[j]->GetBinContent(b+1);

				float esyst2_mm = 0.25 * ss_mm*ss_mm*scale*scale;
				float esyst2_ee = 0.25 * ss_ee*ss_ee*scale*scale;
				float esyst2_em = 0.25 * ss_em*ss_em*scale*scale;

				float estat2_mm = scale*scale*FR->getEStat2(ss_mm);
				float estat2_ee = scale*scale*FR->getEStat2(ss_ee);
				float estat2_em = scale*scale*FR->getEStat2(ss_em);

				float prev    = totbg   ->GetBinError(b+1);
				float prev_mm = totbg_mm->GetBinError(b+1);
				float prev_em = totbg_em->GetBinError(b+1);
				float prev_ee = totbg_ee->GetBinError(b+1);

				totbg   ->SetBinError(b+1, prev    + esyst2_mm + esyst2_ee + esyst2_em + estat2_mm + estat2_ee + estat2_em);
				totbg_mm->SetBinError(b+1, prev_mm + esyst2_mm + estat2_mm);
				totbg_em->SetBinError(b+1, prev_em + esyst2_em + estat2_em);
				totbg_ee->SetBinError(b+1, prev_ee + esyst2_ee + estat2_ee);
			}
			delete FR;
		}

		nt11_ss->Add(nt11_mm_ss);
		nt11_ss->Add(nt11_ee_ss);
		nt11_ss->Add(nt11_em_ss);

		totbg   ->Add(nt11_ss);
		totbg_mm->Add(nt11_mm_ss);
		totbg_em->Add(nt11_em_ss);
		totbg_ee->Add(nt11_ee_ss);

		///////////////////////////////////////////////////////////////////////////////////
		// FAKE PREDICTIONS ///////////////////////////////////////////////////////////////
		///////////////////////////////////////////////////////////////////////////////////
		TH1D *nt11_sf = new TH1D(Form("NT11_SF_%s", varname.Data()), varname, nbins, bins); nt11_sf->Sumw2();
		TH1D *nt11_df = new TH1D(Form("NT11_DF_%s", varname.Data()), varname, nbins, bins); nt11_df->Sumw2();

		TH1D *nt11_mm_sf = new TH1D(Form("NT11_MM_SF_%s", varname.Data()), varname, nbins, bins); nt11_mm_sf->Sumw2();
		TH1D *nt11_ee_sf = new TH1D(Form("NT11_EE_SF_%s", varname.Data()), varname, nbins, bins); nt11_ee_sf->Sumw2();
		TH1D *nt11_em_sf = new TH1D(Form("NT11_EM_SF_%s", varname.Data()), varname, nbins, bins); nt11_em_sf->Sumw2();
		TH1D *nt11_mm_df = new TH1D(Form("NT11_MM_DF_%s", varname.Data()), varname, nbins, bins); nt11_mm_df->Sumw2();
		TH1D *nt11_ee_df = new TH1D(Form("NT11_EE_DF_%s", varname.Data()), varname, nbins, bins); nt11_ee_df->Sumw2();
		TH1D *nt11_em_df = new TH1D(Form("NT11_EM_DF_%s", varname.Data()), varname, nbins, bins); nt11_em_df->Sumw2();

		// Differential ratios
		for(size_t i = 0; i < musamples.size(); ++i){
			Sample *S = fSamples[musamples[i]];
			nt11_mm_sf->Add(S->diffyields[Muon].hnpf[j]);
			nt11_mm_sf->Add(S->diffyields[Muon].hnfp[j]);
			nt11_mm_df->Add(S->diffyields[Muon].hnff[j]);
		}
		for(size_t i = 0; i < elsamples.size(); ++i){
			Sample *S = fSamples[elsamples[i]];
			nt11_ee_sf->Add(S->diffyields[Elec].hnpf[j]);
			nt11_ee_sf->Add(S->diffyields[Elec].hnfp[j]);
			nt11_ee_df->Add(S->diffyields[Elec].hnff[j]);
		}
		for(size_t i = 0; i < emusamples.size(); ++i){
			Sample *S = fSamples[emusamples[i]];
			nt11_em_sf->Add(S->diffyields[ElMu].hnpf[j]);
			nt11_em_sf->Add(S->diffyields[ElMu].hnfp[j]);
			nt11_em_df->Add(S->diffyields[ElMu].hnff[j]);
		}

		for(size_t i = 0; i < nbins; ++i){
			const float FakeESyst2 = 0.25;
			FakeRatios *FR = new FakeRatios();
			FR->setNToyMCs(100); // speedup
			FR->setAddESyst(0.5); // additional systematics

			FR->setMFRatio(mufratio_data, mufratio_data_e); // set error to pure statistical of ratio
			FR->setEFRatio(elfratio_data, elfratio_data_e);
			FR->setMPRatio(mupratio_data, mupratio_data_e);
			FR->setEPRatio(elpratio_data, elpratio_data_e);

			FR->setMMNtl(nt11_mm->GetBinContent(i+1), nt10_mm->GetBinContent(i+1) + nt01_mm->GetBinContent(i+1), nt00_mm->GetBinContent(i+1));
			FR->setEENtl(nt11_ee->GetBinContent(i+1), nt10_ee->GetBinContent(i+1) + nt01_ee->GetBinContent(i+1), nt00_ee->GetBinContent(i+1));
			FR->setEMNtl(nt11_em->GetBinContent(i+1), nt10_em->GetBinContent(i+1),  nt01_em->GetBinContent(i+1), nt00_em->GetBinContent(i+1));
			
			// nt11_mm_sf->SetBinContent(i+1, FR->getMMNpf());
			// nt11_ee_sf->SetBinContent(i+1, FR->getEENpf());
			// nt11_em_sf->SetBinContent(i+1, FR->getEMNpf() + FR->getEMNfp());
			// nt11_mm_df->SetBinContent(i+1, FR->getMMNff());
			// nt11_ee_df->SetBinContent(i+1, FR->getEENff());
			// nt11_em_df->SetBinContent(i+1, FR->getEMNff());
			
			float mm_tot_fakes = nt11_mm_sf->GetBinContent(i+1) + nt11_mm_df->GetBinContent(i+1);
			float ee_tot_fakes = nt11_ee_sf->GetBinContent(i+1) + nt11_ee_df->GetBinContent(i+1);
			float em_tot_fakes = nt11_em_sf->GetBinContent(i+1) + nt11_em_df->GetBinContent(i+1);
			float tot_fakes = mm_tot_fakes + ee_tot_fakes + em_tot_fakes;
			
			// Errors (add total errors of fakes)
			float esyst2_mm  = FakeESyst2*mm_tot_fakes*mm_tot_fakes;
			float esyst2_ee  = FakeESyst2*ee_tot_fakes*ee_tot_fakes;
			float esyst2_em  = FakeESyst2*em_tot_fakes*em_tot_fakes;
			float esyst2_tot = FakeESyst2*tot_fakes*tot_fakes;
			float estat2_mm  = FR->getMMTotEStat()*FR->getMMTotEStat();
			float estat2_ee  = FR->getEETotEStat()*FR->getEETotEStat();
			float estat2_em  = FR->getEMTotEStat()*FR->getEMTotEStat();
			float estat2_tot = FR->getTotEStat()  *FR->getTotEStat();

			float prev    = totbg   ->GetBinError(i+1);
			float prev_mm = totbg_mm->GetBinError(i+1);
			float prev_em = totbg_em->GetBinError(i+1);
			float prev_ee = totbg_ee->GetBinError(i+1);

			totbg   ->SetBinError(i+1, prev    + esyst2_tot + estat2_tot);
			totbg_mm->SetBinError(i+1, prev_mm + esyst2_mm + estat2_mm);
			totbg_em->SetBinError(i+1, prev_em + esyst2_em + estat2_em);
			totbg_ee->SetBinError(i+1, prev_ee + esyst2_ee + estat2_ee);
			
			delete FR;
		}
		
		nt11_sf->Add(nt11_mm_sf);
		nt11_sf->Add(nt11_ee_sf);
		nt11_sf->Add(nt11_em_sf);

		nt11_df->Add(nt11_mm_df);
		nt11_df->Add(nt11_ee_df);
		nt11_df->Add(nt11_em_df);

		totbg   ->Add(nt11_sf);
		totbg   ->Add(nt11_df);
		totbg_mm->Add(nt11_mm_sf);
		totbg_mm->Add(nt11_mm_df);
		totbg_em->Add(nt11_em_sf);
		totbg_em->Add(nt11_em_df);
		totbg_ee->Add(nt11_ee_sf);
		totbg_ee->Add(nt11_ee_df);

		///////////////////////////////////////////////////////////////////////////////////
		// E-CHARGE MISID /////////////////////////////////////////////////////////////////
		///////////////////////////////////////////////////////////////////////////////////
		TH1D *nt11_cm = new TH1D(Form("NT11_CM_%s", varname.Data()), varname, nbins, bins); nt11_cm->Sumw2();

		TH1D *nt11_ee_cm = new TH1D(Form("NT11_EE_CM_%s", varname.Data()), varname, nbins, bins); nt11_ee_cm->Sumw2();
		TH1D *nt11_em_cm = new TH1D(Form("NT11_EM_CM_%s", varname.Data()), varname, nbins, bins); nt11_em_cm->Sumw2();

		// Abbreviations
		float fb  = gEChMisIDB;
		float fbE = gEChMisIDB_E;
		float fe  = gEChMisIDE;
		float feE = gEChMisIDE_E;

		for(size_t i = 0; i < nbins; ++i){
			float nt2_ee_BB_os = nt2_os_ee_bb->GetBinContent(i+1);
			float nt2_ee_EB_os = nt2_os_ee_eb->GetBinContent(i+1);
			float nt2_ee_EE_os = nt2_os_ee_ee->GetBinContent(i+1);
			float nt2_em_BB_os = nt2_os_em_bb->GetBinContent(i+1);
			float nt2_em_EE_os = nt2_os_em_ee->GetBinContent(i+1);
			
			// Errors
			FakeRatios *FR = new FakeRatios();

			// Simple error propagation assuming error on number of events is FR->getEStat2()
			nt11_ee_cm->SetBinContent(i+1, 2*fb*nt2_ee_BB_os + 2*fe*nt2_ee_EE_os + (fb+fe)*nt2_ee_EB_os);
			float nt11_ee_cm_e1 = sqrt( (4*fb*fb*FR->getEStat2(nt2_ee_BB_os)) + (4*fe*fe*FR->getEStat2(nt2_ee_EE_os)) + (fb+fe)*(fb+fe)*FR->getEStat2(nt2_ee_EB_os) ); // stat only
			float nt11_ee_cm_e2 = sqrt( (4*nt2_ee_BB_os*nt2_ee_BB_os*fbE*fbE) + (4*nt2_ee_EE_os*nt2_ee_EE_os*feE*feE) + (fbE*fbE+feE*feE)*nt2_ee_EB_os*nt2_ee_EB_os ); // syst only

			nt11_em_cm->SetBinContent(i+i, fb*nt2_em_BB_os + fe*nt2_em_EE_os);
			float nt11_em_cm_e1 = sqrt( fb*fb*FR->getEStat2(nt2_em_BB_os) + fe*fe*FR->getEStat2(nt2_em_EE_os) );
			float nt11_em_cm_e2 = sqrt( nt2_em_BB_os*nt2_em_BB_os * fbE*fbE + nt2_em_EE_os*nt2_em_EE_os * feE*feE );
			
			float esyst2_ee  = nt11_ee_cm_e2*nt11_ee_cm_e2;
			float esyst2_em  = nt11_em_cm_e2*nt11_em_cm_e2;
			float esyst2_tot = nt11_ee_cm_e2*nt11_ee_cm_e2 + nt11_em_cm_e2*nt11_em_cm_e2;
			float estat2_ee  = nt11_ee_cm_e1*nt11_ee_cm_e1;
			float estat2_em  = nt11_em_cm_e1*nt11_em_cm_e1;
			float estat2_tot = nt11_ee_cm_e1*nt11_ee_cm_e1 + nt11_em_cm_e1*nt11_em_cm_e1;

			float prev    = totbg   ->GetBinError(i+1);
			float prev_em = totbg_em->GetBinError(i+1);
			float prev_ee = totbg_ee->GetBinError(i+1);

			totbg   ->SetBinError(i+1, prev    + esyst2_tot + estat2_tot);
			totbg_em->SetBinError(i+1, prev_em + esyst2_em  + estat2_em);
			totbg_ee->SetBinError(i+1, prev_ee + esyst2_ee  + estat2_ee);
			delete FR;
		}

		nt11_cm->Add(nt11_ee_cm);
		nt11_cm->Add(nt11_em_cm);

		totbg   ->Add(nt11_cm);
		totbg_em->Add(nt11_em_cm);
		totbg_ee->Add(nt11_ee_cm);

		///////////////////////////////////////////////////////////////////////////////////
		// SIGNAL /////////////////////////////////////////////////////////////////////////
		///////////////////////////////////////////////////////////////////////////////////
		gSample sigsam = LM11;
		TH1D *nt11_sig    = new TH1D(Form("NT11_Sig%s",    varname.Data()), varname, nbins, bins); nt11_sig   ->Sumw2();
		TH1D *nt11_mm_sig = new TH1D(Form("NT11_mm_Sig%s", varname.Data()), varname, nbins, bins); nt11_mm_sig->Sumw2();
		TH1D *nt11_em_sig = new TH1D(Form("NT11_em_Sig%s", varname.Data()), varname, nbins, bins); nt11_em_sig->Sumw2();
		TH1D *nt11_ee_sig = new TH1D(Form("NT11_ee_Sig%s", varname.Data()), varname, nbins, bins); nt11_ee_sig->Sumw2();
		nt11_mm_sig->Add(fSamples[sigsam]->diffyields[Muon].hnt11[j], fLumiNorm / fSamples[LM4]->lumi);
		nt11_ee_sig->Add(fSamples[sigsam]->diffyields[Elec].hnt11[j], fLumiNorm / fSamples[LM4]->lumi);
		nt11_em_sig->Add(fSamples[sigsam]->diffyields[ElMu].hnt11[j], fLumiNorm / fSamples[LM4]->lumi);
		nt11_sig   ->Add(fSamples[sigsam]->diffyields[Muon].hnt11[j], fLumiNorm / fSamples[LM4]->lumi);
		nt11_sig   ->Add(fSamples[sigsam]->diffyields[Elec].hnt11[j], fLumiNorm / fSamples[LM4]->lumi);
		nt11_sig   ->Add(fSamples[sigsam]->diffyields[ElMu].hnt11[j], fLumiNorm / fSamples[LM4]->lumi);

		///////////////////////////////////////////////////////////////////////////////////
		// OUTPUT /////////////////////////////////////////////////////////////////////////
		///////////////////////////////////////////////////////////////////////////////////
		nt11->SetMarkerColor(kBlack);
		nt11->SetMarkerStyle(20);
		nt11->SetMarkerSize(2.0);
		nt11->SetLineWidth(2);
		nt11->SetLineColor(kBlack);
		nt11->SetFillColor(kBlack);
		nt11_mm->SetMarkerColor(kBlack);
		nt11_mm->SetMarkerStyle(20);
		nt11_mm->SetMarkerSize(2.0);
		nt11_mm->SetLineWidth(2);
		nt11_mm->SetLineColor(kBlack);
		nt11_mm->SetFillColor(kBlack);
		nt11_ee->SetMarkerColor(kBlack);
		nt11_ee->SetMarkerStyle(20);
		nt11_ee->SetMarkerSize(2.0);
		nt11_ee->SetLineWidth(2);
		nt11_ee->SetLineColor(kBlack);
		nt11_ee->SetFillColor(kBlack);
		nt11_em->SetMarkerColor(kBlack);
		nt11_em->SetMarkerStyle(20);
		nt11_em->SetMarkerSize(2.0);
		nt11_em->SetLineWidth(2);
		nt11_em->SetLineColor(kBlack);
		nt11_em->SetFillColor(kBlack);

		nt11_sf->SetLineWidth(1);
		nt11_df->SetLineWidth(1);
		nt11_ss->SetLineWidth(1);
		nt11_sf->SetLineColor(50);
		nt11_sf->SetFillColor(50);
		nt11_df->SetLineColor(38);
		nt11_df->SetFillColor(38);
		nt11_cm->SetLineColor(42);
		nt11_cm->SetFillColor(42);
		nt11_ss->SetLineColor(31);
		nt11_ss->SetFillColor(31);

		nt11_mm_sf->SetLineWidth(1);
		nt11_mm_df->SetLineWidth(1);
		nt11_mm_ss->SetLineWidth(1);
		nt11_mm_sf->SetLineColor(50);
		nt11_mm_sf->SetFillColor(50);
		nt11_mm_df->SetLineColor(38);
		nt11_mm_df->SetFillColor(38);
		nt11_mm_ss->SetLineColor(31);
		nt11_mm_ss->SetFillColor(31);

		nt11_ee_sf->SetLineWidth(1);
		nt11_ee_df->SetLineWidth(1);
		nt11_ee_ss->SetLineWidth(1);
		nt11_ee_sf->SetLineColor(50);
		nt11_ee_sf->SetFillColor(50);
		nt11_ee_df->SetLineColor(38);
		nt11_ee_df->SetFillColor(38);
		nt11_ee_cm->SetLineColor(42);
		nt11_ee_cm->SetFillColor(42);
		nt11_ee_ss->SetLineColor(31);
		nt11_ee_ss->SetFillColor(31);

		nt11_em_sf->SetLineWidth(1);
		nt11_em_df->SetLineWidth(1);
		nt11_em_ss->SetLineWidth(1);
		nt11_em_sf->SetLineColor(50);
		nt11_em_sf->SetFillColor(50);
		nt11_em_df->SetLineColor(38);
		nt11_em_df->SetFillColor(38);
		nt11_em_cm->SetLineColor(42);
		nt11_em_cm->SetFillColor(42);
		nt11_em_ss->SetLineColor(31);
		nt11_em_ss->SetFillColor(31);
		
		nt11_sig   ->SetLineWidth(2);
		nt11_mm_sig->SetLineWidth(2);
		nt11_em_sig->SetLineWidth(2);
		nt11_ee_sig->SetLineWidth(2);

		nt11_sig   ->SetLineColor(kBlue);
		nt11_mm_sig->SetLineColor(kBlue);
		nt11_em_sig->SetLineColor(kBlue);
		nt11_ee_sig->SetLineColor(kBlue);
		
		nt11_sig   ->SetFillStyle(0);
		nt11_mm_sig->SetFillStyle(0);
		nt11_em_sig->SetFillStyle(0);
		nt11_ee_sig->SetFillStyle(0);
		
		totbg   ->SetLineWidth(1);
		totbg_mm->SetLineWidth(1);
		totbg_em->SetLineWidth(1);
		totbg_ee->SetLineWidth(1);


		totbg   ->SetFillColor(12);
		totbg_mm->SetFillColor(12);
		totbg_em->SetFillColor(12);
		totbg_ee->SetFillColor(12);
		totbg   ->SetFillStyle(3005);
		totbg_mm->SetFillStyle(3005);
		totbg_em->SetFillStyle(3005);
		totbg_ee->SetFillStyle(3005);

		// Take square root of sum of squared errors:
		// (stored the SQUARED errors before)
		for(size_t i = 0; i < nbins; ++i){
			float prev    = totbg   ->GetBinError(i+1);
			float prev_mm = totbg_mm->GetBinError(i+1);
			float prev_em = totbg_em->GetBinError(i+1);
			float prev_ee = totbg_ee->GetBinError(i+1);

			totbg   ->SetBinError(i+1, sqrt(prev)   );
			totbg_mm->SetBinError(i+1, sqrt(prev_mm));
			totbg_em->SetBinError(i+1, sqrt(prev_em));
			totbg_ee->SetBinError(i+1, sqrt(prev_ee));
		}

		// Normalize everything to binwidth
		nt11_sf    = normHistBW(nt11_sf, binwidthscale[j]);
		nt11_df    = normHistBW(nt11_df, binwidthscale[j]);
		nt11_ss    = normHistBW(nt11_ss, binwidthscale[j]);
		nt11_cm    = normHistBW(nt11_cm, binwidthscale[j]);

		nt11_mm_sf = normHistBW(nt11_mm_sf, binwidthscale[j]);
		nt11_mm_df = normHistBW(nt11_mm_df, binwidthscale[j]);
		nt11_mm_ss = normHistBW(nt11_mm_ss, binwidthscale[j]);

		nt11_ee_sf = normHistBW(nt11_ee_sf, binwidthscale[j]);
		nt11_ee_df = normHistBW(nt11_ee_df, binwidthscale[j]);
		nt11_ee_ss = normHistBW(nt11_ee_ss, binwidthscale[j]);
		nt11_ee_cm = normHistBW(nt11_ee_cm, binwidthscale[j]);
		
		nt11_em_sf = normHistBW(nt11_em_sf, binwidthscale[j]);
		nt11_em_df = normHistBW(nt11_em_df, binwidthscale[j]);
		nt11_em_ss = normHistBW(nt11_em_ss, binwidthscale[j]);
		nt11_em_cm = normHistBW(nt11_em_cm, binwidthscale[j]);
		
		totbg      = normHistBW(totbg,    binwidthscale[j]);
		totbg_mm   = normHistBW(totbg_mm, binwidthscale[j]);
		totbg_em   = normHistBW(totbg_em, binwidthscale[j]);
		totbg_ee   = normHistBW(totbg_ee, binwidthscale[j]);

		nt11       = normHistBW(nt11,    binwidthscale[j]);
		nt11_mm    = normHistBW(nt11_mm, binwidthscale[j]);
		nt11_em    = normHistBW(nt11_em, binwidthscale[j]);
		nt11_ee    = normHistBW(nt11_ee, binwidthscale[j]);

		// Fill stacks
		THStack *nt11_tot    = new THStack("NT11_TotalBG", "NT11_TotalBG");
		THStack *nt11_mm_tot = new THStack("NT11_MM_TotalBG", "NT11_MM_TotalBG");
		THStack *nt11_ee_tot = new THStack("NT11_EE_TotalBG", "NT11_EE_TotalBG");
		THStack *nt11_em_tot = new THStack("NT11_EM_TotalBG", "NT11_EM_TotalBG");
		
		nt11_tot->Add(nt11_sf);
		nt11_tot->Add(nt11_df);
		nt11_tot->Add(nt11_ss);
		nt11_tot->Add(nt11_cm);

		nt11_mm_tot->Add(nt11_mm_sf);
		nt11_mm_tot->Add(nt11_mm_df);
		nt11_mm_tot->Add(nt11_mm_ss);

		nt11_ee_tot->Add(nt11_ee_sf);
		nt11_ee_tot->Add(nt11_ee_df);
		nt11_ee_tot->Add(nt11_ee_ss);
		nt11_ee_tot->Add(nt11_ee_cm);

		nt11_em_tot->Add(nt11_em_sf);
		nt11_em_tot->Add(nt11_em_df);
		nt11_em_tot->Add(nt11_em_ss);
		nt11_em_tot->Add(nt11_em_cm);

		// Signal
		// nt11_tot->Add(nt11_sig);
		// nt11_mm_tot->Add(nt11_mm_sig);
		// nt11_ee_tot->Add(nt11_ee_sig);
		// nt11_em_tot->Add(nt11_em_sig);

		TString ytitle = Form("Events / %3.0f GeV", binwidthscale[j]);
		if(j==4 || j==8) ytitle = "Events";
		nt11_tot->Draw("goff");
		nt11_tot->GetXaxis()->SetTitle(DiffPredYields::axis_label[j]);
		nt11_tot->GetYaxis()->SetTitle(ytitle);
		nt11_mm_tot->Draw("goff");
		nt11_mm_tot->GetXaxis()->SetTitle(DiffPredYields::axis_label[j]);
		nt11_mm_tot->GetYaxis()->SetTitle(ytitle);
		nt11_ee_tot->Draw("goff");
		nt11_ee_tot->GetXaxis()->SetTitle(DiffPredYields::axis_label[j]);
		nt11_ee_tot->GetYaxis()->SetTitle(ytitle);
		nt11_em_tot->Draw("goff");
		nt11_em_tot->GetXaxis()->SetTitle(DiffPredYields::axis_label[j]);
		nt11_em_tot->GetYaxis()->SetTitle(ytitle);

		nt11_tot   ->SetMinimum(0.5*nt11   ->GetMinimum());
		nt11_mm_tot->SetMinimum(0.5*nt11_mm->GetMinimum());
		nt11_ee_tot->SetMinimum(0.5*nt11_ee->GetMinimum());
		nt11_em_tot->SetMinimum(0.5*nt11_em->GetMinimum());

		double max = nt11->Integral();
		nt11    ->SetMaximum(max>1?max+1:1.);
		nt11_sf ->SetMaximum(max>1?max+1:1.);
		nt11_df ->SetMaximum(max>1?max+1:1.);
		nt11_cm ->SetMaximum(max>1?max+1:1.);
		nt11_ss ->SetMaximum(max>1?max+1:1.);
		nt11_tot->SetMaximum(max>1?max+1:1.);

		double max_mm = nt11_mm->Integral();
		nt11_mm    ->SetMaximum(max_mm>1?max_mm+1:1.);
		nt11_mm_sf ->SetMaximum(max_mm>1?max_mm+1:1.);
		nt11_mm_df ->SetMaximum(max_mm>1?max_mm+1:1.);
		nt11_mm_ss ->SetMaximum(max_mm>1?max_mm+1:1.);
		nt11_mm_tot->SetMaximum(max_mm>1?max_mm+1:1.);

		double max_ee = nt11_ee->Integral();
		nt11_ee    ->SetMaximum(max_ee>1?max_ee+1:1.);
		nt11_ee_sf ->SetMaximum(max_ee>1?max_ee+1:1.);
		nt11_ee_df ->SetMaximum(max_ee>1?max_ee+1:1.);
		nt11_ee_cm ->SetMaximum(max_ee>1?max_ee+1:1.);
		nt11_ee_ss ->SetMaximum(max_ee>1?max_ee+1:1.);
		nt11_ee_tot->SetMaximum(max_ee>1?max_ee+1:1.);

		double max_em = nt11_em->Integral();
		nt11_em    ->SetMaximum(max_em>1?max_em+1:1.);
		nt11_em_sf ->SetMaximum(max_em>1?max_em+1:1.);
		nt11_em_df ->SetMaximum(max_em>1?max_em+1:1.);
		nt11_em_cm ->SetMaximum(max_em>1?max_em+1:1.);
		nt11_em_ss ->SetMaximum(max_em>1?max_em+1:1.);
		nt11_em_tot->SetMaximum(max_em>1?max_em+1:1.);
		
		fOutputSubDir = "DiffPredictionPlots/";
		/////////////////////////////////////////////////////////////////
		TLegend *leg = new TLegend(0.60,0.67,0.90,0.88);
		leg->AddEntry(nt11,    "Observed","p");
		leg->AddEntry(nt11_sf, "Single Fakes","f");
		leg->AddEntry(nt11_df, "Double Fakes","f");
		leg->AddEntry(nt11_ss, "Irreducible (MC)","f");
		leg->AddEntry(nt11_cm, "Charge MisID","f");
		leg->AddEntry(totbg,   "Total Uncertainty","f");
		// leg->AddEntry(nt11_sig,fSamples[sigsam]->sname,"l");
		leg->SetFillStyle(0);
		leg->SetTextFont(42);
		leg->SetBorderSize(0);
		
		TCanvas *c_temp = new TCanvas("C_ObsPred_" + varname, "Observed vs Predicted", 0, 0, 800, 600);
		c_temp->cd();
		gPad->SetLogy();
		
		nt11_tot->Draw("hist");
		// nt11_error->DrawCopy("X0 E1 same");
		nt11->DrawCopy("PE X0 same");
		totbg->DrawCopy("0 E2 same");
		// nt11_sig->DrawCopy("hist same");
		leg->Draw();
		lat->SetTextSize(0.04);
		lat->DrawLatex(0.55,0.92, "#mu#mu/ee/e#mu");
		drawDiffCuts(j);
		drawTopLine();
		
		gPad->RedrawAxis();
		Util::PrintPDF(c_temp, "ObsPred_" + varname, fOutputDir + fOutputSubDir);
		gPad->SetLogy(0);
		float minopt, maxopt;
		vector<TH1D*> histvec;
		histvec.push_back(nt11);
		histvec.push_back(totbg);
		getPlottingRange(minopt, maxopt, histvec, 0.1);
		nt11_tot->SetMinimum(0);
		nt11_tot->SetMaximum(maxopt);
		Util::PrintPDF(c_temp, "ObsPred_" + varname + "_lin", fOutputDir + fOutputSubDir + "lin/");

		fOutputSubDir = "DiffPredictionPlots/IndividualChannels/";
		/////////////////////////////////////////////////////////////////
		TLegend *leg_mm = new TLegend(0.60,0.67,0.90,0.88);
		leg_mm->AddEntry(nt11_mm,    "Observed","p");
		leg_mm->AddEntry(nt11_mm_sf, "Single Fakes","f");
		leg_mm->AddEntry(nt11_mm_df, "Double Fakes","f");
		leg_mm->AddEntry(nt11_mm_ss, "Irreducible (MC)","f");
		leg_mm->AddEntry(totbg_mm,   "Total Uncertainty","f");
		// leg_mm->AddEntry(nt11_mm_sig,fSamples[sigsam]->sname,"l");
		leg_mm->SetFillStyle(0);
		leg_mm->SetTextFont(42);
		leg_mm->SetBorderSize(0);
		
		c_temp = new TCanvas("C_ObsPred_MM_" + varname, "Observed vs Predicted", 0, 0, 800, 600);
		c_temp->cd();
		gPad->SetLogy();
		
		nt11_mm_tot->Draw("hist");
		// nt11_error->DrawCopy("X0 E1 same");
		nt11_mm->DrawCopy("PE X0 same");
		totbg_mm->DrawCopy("0 E2 same");
		// nt11_mm_sig->DrawCopy("hist same");
		leg_mm->Draw();
		lat->SetTextSize(0.04);
		lat->DrawLatex(0.65,0.92, "#mu#mu");
		drawDiffCuts(j);
		drawTopLine();
		
		gPad->RedrawAxis();
		Util::PrintPDF(c_temp, varname + "_MM_ObsPred", fOutputDir + fOutputSubDir);

		histvec.clear();
		histvec.push_back(nt11_mm);
		histvec.push_back(totbg_mm);
		getPlottingRange(minopt, maxopt, histvec, 0.1);
		nt11_mm_tot->SetMinimum(0);
		nt11_mm_tot->SetMaximum(maxopt);
		gPad->SetLogy(0);
		Util::PrintPDF(c_temp, varname + "_MM_ObsPred_lin", fOutputDir + fOutputSubDir + "lin/");

		/////////////////////////////////////////////////////////////////
		TLegend *leg_ee = new TLegend(0.60,0.67,0.90,0.88);
		leg_ee->AddEntry(nt11_ee,    "Observed","p");
		leg_ee->AddEntry(nt11_ee_sf, "Single Fakes","f");
		leg_ee->AddEntry(nt11_ee_df, "Double Fakes","f");
		leg_ee->AddEntry(nt11_ee_ss, "Irreducible (MC)","f");
		leg_ee->AddEntry(nt11_ee_cm, "Charge MisID","f");
		leg_ee->AddEntry(totbg_ee,   "Total Uncertainty","f");
		// leg_mm->AddEntry(nt11_ee_sig,fSamples[sigsam]->sname,"l");
		leg_ee->SetFillStyle(0);
		leg_ee->SetTextFont(42);
		leg_ee->SetBorderSize(0);
		
		c_temp = new TCanvas("C_ObsPred_EE_" + varname, "Observed vs Predicted", 0, 0, 800, 600);
		c_temp->cd();
		gPad->SetLogy();
		
		nt11_ee_tot->Draw("hist");
		// nt11_error->DrawCopy("X0 E1 same");
		nt11_ee->DrawCopy("PE X0 same");
		totbg_ee->DrawCopy("0 E2 same");
		// nt11_ee_sig->DrawCopy("hist same");
		leg_ee->Draw();
		lat->SetTextSize(0.04);
		lat->DrawLatex(0.65,0.92, "ee");
		drawDiffCuts(j);
		drawTopLine();
		
		gPad->RedrawAxis();
		Util::PrintPDF(c_temp, varname + "_EE_ObsPred", fOutputDir + fOutputSubDir);

		histvec.clear();
		histvec.push_back(nt11_ee);
		histvec.push_back(totbg_ee);
		getPlottingRange(minopt, maxopt, histvec, 0.1);
		nt11_ee_tot->SetMinimum(0);
		nt11_ee_tot->SetMaximum(maxopt);
		gPad->SetLogy(0);
		Util::PrintPDF(c_temp, varname + "_EE_ObsPred_lin", fOutputDir + fOutputSubDir + "lin/");

		/////////////////////////////////////////////////////////////////
		TLegend *leg_em = new TLegend(0.60,0.67,0.90,0.88);
		leg_em->AddEntry(nt11_em,    "Observed","p");
		leg_em->AddEntry(nt11_em_sf, "Single Fakes","f");
		leg_em->AddEntry(nt11_em_df, "Double Fakes","f");
		leg_em->AddEntry(nt11_em_ss, "Irreducible (MC)","f");
		leg_em->AddEntry(nt11_em_cm, "Charge MisID","f");
		leg_em->AddEntry(totbg_em,   "Total Uncertainty","f");
		// leg_mm->AddEntry(nt11_em_sig,fSamples[sigsam]->sname,"l");
		leg_em->SetFillStyle(0);
		leg_em->SetTextFont(42);
		leg_em->SetBorderSize(0);
		
		c_temp = new TCanvas("C_ObsPred_EM_" + varname, "Observed vs Predicted", 0, 0, 800, 600);
		c_temp->cd();
		gPad->SetLogy();
		
		nt11_em_tot->Draw("hist");
		totbg_em->DrawCopy("0 E2 same");
		nt11_em->DrawCopy("PE X0 same");
		// nt11_em_sig->DrawCopy("hist same");
		leg_em->Draw();
		lat->SetTextSize(0.04);
		lat->DrawLatex(0.65,0.92, "e#mu");
		drawDiffCuts(j);
		drawTopLine();
		
		gPad->RedrawAxis();
		Util::PrintPDF(c_temp, varname + "_EM_ObsPred", fOutputDir + fOutputSubDir);

		histvec.clear();
		histvec.push_back(nt11_em);
		histvec.push_back(totbg_em);
		getPlottingRange(minopt, maxopt, histvec, 0.1);
		nt11_em_tot->SetMinimum(0);
		nt11_em_tot->SetMaximum(maxopt);
		gPad->SetLogy(0);
		Util::PrintPDF(c_temp, varname + "_EM_ObsPred_lin", fOutputDir + fOutputSubDir + "lin/");

		// Cleanup
		delete c_temp, leg, leg_mm, leg_em, leg_ee;
		delete nt11, nt11_mm, nt11_ee, nt11_em;
		delete nt11_sig, nt11_mm_sig, nt11_ee_sig, nt11_em_sig;
		delete nt10_mm, nt10_em, nt10_ee, nt01_mm, nt01_em, nt01_ee, nt00_mm, nt00_em, nt00_ee;
		delete nt2_os_ee_bb, nt2_os_ee_eb, nt2_os_ee_ee, nt2_os_em_bb, nt2_os_em_ee;
		delete nt11_ss, nt11_mm_ss, nt11_em_ss, nt11_ee_ss;
		delete nt11_sf, nt11_mm_sf, nt11_em_sf, nt11_ee_sf;
		delete nt11_df, nt11_mm_df, nt11_em_df, nt11_ee_df;
		delete nt11_cm, nt11_em_cm, nt11_ee_cm;
		// delete nt11_mc, nt11_mm_mc, nt11_ee_mc, nt11_em_mc;
		delete nt11_tot, nt11_mm_tot, nt11_ee_tot, nt11_em_tot;
		delete totbg, totbg_mm, totbg_em, totbg_ee;
	}
	fOutputSubDir = "";
}

//____________________________________________________________________________
void SSDLPlotter::makeAllClosureTests(){
	TString outputdir = Util::MakeOutputDir(fOutputDir + "MCClosure");
	for(size_t i = 0; i < gNREGIONS; ++i){
		TString outputname = outputdir + "MCClosure_" + Region::sname[i] + ".txt";
		makeIntMCClosure(fMCBG, outputname, gRegion(i));
	}
	for(size_t i = 0; i < gNREGIONS; ++i){
		TString outputname = outputdir + "MCClosure_Sig_" + Region::sname[i] + ".txt";
		makeIntMCClosure(fMCBGSig, outputname, gRegion(i));
	}
	fOutputSubDir = "";
}
void SSDLPlotter::makeIntMCClosure(vector<int> samples, TString filename, gRegion reg, gHiLoSwitch hilo){
	ofstream OUT(filename.Data(), ios::trunc);
	const int nsamples = samples.size();

	OUT << "/////////////////////////////////////////////////////////////////////////////" << endl;
	OUT << " Producing integrated predictions" << endl;
	OUT << "  scaling MC to " << fLumiNorm << " /pb" << endl << endl;

	///////////////////////////////////////////////////////////////////////////////////
	// RATIOS /////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////
	float mf(0.), mf_e(0.), mp(0.), mp_e(0.), ef(0.), ef_e(0.), ep(0.), ep_e(0.);

	calculateRatio(fMCBG, Muon, SigSup, mf, mf_e);
	calculateRatio(fMCBG, Muon, ZDecay, mp, mp_e);
	calculateRatio(fMCBG, Elec, SigSup, ef, ef_e);
	calculateRatio(fMCBG, Elec, ZDecay, ep, ep_e);

	///////////////////////////////////////////////////////////////////////////////////
	// OBSERVATIONS ///////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////
	vector<float> ntt_mm, ntl_mm, nll_mm;
	vector<float> ntt_ee, ntl_ee, nll_ee;
	vector<float> ntt_em, ntl_em, nlt_em, nll_em;

	vector<float> npp_mm, npf_mm, nfp_mm, nff_mm;
	vector<float> npp_ee, npf_ee, nfp_ee, nff_ee;
	vector<float> npp_em, npf_em, nfp_em, nff_em;

	vector<float> npp_tt_mm, npf_tt_mm, nfp_tt_mm, nff_tt_mm;
	vector<float> npp_tt_ee, npf_tt_ee, nfp_tt_ee, nff_tt_ee;
	vector<float> npp_tt_em, npf_tt_em, nfp_tt_em, nff_tt_em;

	// OS yields
	vector<float> ntt_os_BB_ee, ntt_os_EE_ee, ntt_os_EB_ee;
	vector<float> ntt_os_BB_em, ntt_os_EE_em;
	// Charge misid
	vector<float> npp_tt_cm_ee, npp_cm_ee;
	vector<float> npp_tt_cm_em, npp_cm_em;

	vector<float> scales;
	vector<TString> names;
	for(size_t i = 0; i < samples.size(); ++i){
		Sample *S = fSamples[samples[i]];
		float scale = fLumiNorm / S->lumi;
		names.push_back(S->sname);
		scales.push_back(scale);
		ntt_mm.push_back(S->numbers[reg][Muon].nt2);
		ntl_mm.push_back(S->numbers[reg][Muon].nt10);
		nll_mm.push_back(S->numbers[reg][Muon].nt0);

		ntt_em.push_back(S->numbers[reg][ElMu].nt2);
		ntl_em.push_back(S->numbers[reg][ElMu].nt10);
		nlt_em.push_back(S->numbers[reg][ElMu].nt01);
		nll_em.push_back(S->numbers[reg][ElMu].nt0);

		ntt_ee.push_back(S->numbers[reg][Elec].nt2);
		ntl_ee.push_back(S->numbers[reg][Elec].nt10);
		nll_ee.push_back(S->numbers[reg][Elec].nt0);

		npp_mm.push_back(S->region[reg][hilo].mm.npp_pt->GetEntries());
		npf_mm.push_back(S->region[reg][hilo].mm.npf_pt->GetEntries());
		nfp_mm.push_back(S->region[reg][hilo].mm.nfp_pt->GetEntries());
		nff_mm.push_back(S->region[reg][hilo].mm.nff_pt->GetEntries());

		npp_em.push_back(S->region[reg][hilo].em.npp_pt->GetEntries());
		npf_em.push_back(S->region[reg][hilo].em.npf_pt->GetEntries());
		nfp_em.push_back(S->region[reg][hilo].em.nfp_pt->GetEntries());
		nff_em.push_back(S->region[reg][hilo].em.nff_pt->GetEntries());

		npp_ee.push_back(S->region[reg][hilo].ee.npp_pt->GetEntries());
		npf_ee.push_back(S->region[reg][hilo].ee.npf_pt->GetEntries());
		nfp_ee.push_back(S->region[reg][hilo].ee.nfp_pt->GetEntries());
		nff_ee.push_back(S->region[reg][hilo].ee.nff_pt->GetEntries());

		npp_tt_mm.push_back(S->region[reg][hilo].mm.nt2pp_pt->GetEntries());
		npf_tt_mm.push_back(S->region[reg][hilo].mm.nt2pf_pt->GetEntries());
		nfp_tt_mm.push_back(S->region[reg][hilo].mm.nt2fp_pt->GetEntries());
		nff_tt_mm.push_back(S->region[reg][hilo].mm.nt2ff_pt->GetEntries());

		npp_tt_em.push_back(S->region[reg][hilo].em.nt2pp_pt->GetEntries());
		npf_tt_em.push_back(S->region[reg][hilo].em.nt2pf_pt->GetEntries());
		nfp_tt_em.push_back(S->region[reg][hilo].em.nt2fp_pt->GetEntries());
		nff_tt_em.push_back(S->region[reg][hilo].em.nt2ff_pt->GetEntries());

		npp_tt_ee.push_back(S->region[reg][hilo].ee.nt2pp_pt->GetEntries());
		npf_tt_ee.push_back(S->region[reg][hilo].ee.nt2pf_pt->GetEntries());
		nfp_tt_ee.push_back(S->region[reg][hilo].ee.nt2fp_pt->GetEntries());
		nff_tt_ee.push_back(S->region[reg][hilo].ee.nt2ff_pt->GetEntries());
		
		ntt_os_BB_em.push_back(S->region[reg][hilo].em.nt20_OS_BB_pt->GetEntries()); // ele in barrel
		ntt_os_EE_em.push_back(S->region[reg][hilo].em.nt20_OS_EE_pt->GetEntries()); // ele in endcal
		ntt_os_BB_ee.push_back(S->region[reg][hilo].ee.nt20_OS_BB_pt->GetEntries()); // both in barrel
		ntt_os_EE_ee.push_back(S->region[reg][hilo].ee.nt20_OS_EE_pt->GetEntries()); // both in endcal
		ntt_os_EB_ee.push_back(S->region[reg][hilo].ee.nt20_OS_EB_pt->GetEntries()); // one barrel, one endcap
		
		npp_tt_cm_ee.push_back(scale*S->region[reg][hilo].ee.nt2pp_cm_pt->GetEntries());
		npp_cm_ee   .push_back(scale*S->region[reg][hilo].ee.npp_cm_pt->GetEntries());
		npp_tt_cm_em.push_back(scale*S->region[reg][hilo].em.nt2pp_cm_pt->GetEntries());
		npp_cm_em   .push_back(scale*S->region[reg][hilo].em.npp_cm_pt->GetEntries());
	}

	///////////////////////////////////////////////////////////////////////////////////
	// PREDICTIONS ////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////
	FakeRatios *FR = new FakeRatios();
	FR->setNToyMCs(100);
	FR->setAddESyst(0.0);

	///////////////////////////////////////////////////////////////////////////////////
	// Fiddle with ratios by hand...
	// mf = 0.067;
	// ef = 0.25;
	// mp = 0.90;
	// ep = 0.90;
	// Cheating ratios:
	// mf = 0.04;
	// ef = 0.17;
	// mp = 0.90;
	// ep = 0.90;
	///////////////////////////////////////////////////////////////////////////////////

	FR->setMFRatio(mf, mf_e); // set error to pure statistical of ratio
	FR->setEFRatio(ef, ef_e);
	FR->setMPRatio(mp, mp_e);
	FR->setEPRatio(ep, ep_e);

	vector<float> npp_pred_mm,    npf_pred_mm,    nff_pred_mm;
	vector<float> npp_pred_mm_e1, npf_pred_mm_e1, nff_pred_mm_e1;
	vector<float> npp_pred_ee,    npf_pred_ee,    nff_pred_ee;
	vector<float> npp_pred_ee_e1, npf_pred_ee_e1, nff_pred_ee_e1;
	vector<float> npp_pred_em,    npf_pred_em,    nfp_pred_em,    nff_pred_em;
	vector<float> npp_pred_em_e1, npf_pred_em_e1, nfp_pred_em_e1, nff_pred_em_e1;

	vector<float> nF_pred_mm_e1, nF_pred_em_e1, nF_pred_ee_e1; // combined stat errors on fakes

	for(size_t i = 0; i < nsamples; ++i){
		FR->setMMNtl(ntt_mm[i], ntl_mm[i], nll_mm[i]);
		FR->setEENtl(ntt_ee[i], ntl_ee[i], nll_ee[i]);
		FR->setEMNtl(ntt_em[i], ntl_em[i], nlt_em[i], nll_em[i]);

		Sample *S = fSamples[samples[i]];
		npp_pred_mm.push_back(S->numbers[reg][Muon].npp);
		npf_pred_mm.push_back(S->numbers[reg][Muon].npf + S->numbers[reg][Muon].nfp);
		nff_pred_mm.push_back(S->numbers[reg][Muon].nff);
		npp_pred_ee.push_back(S->numbers[reg][Elec].npp);
		npf_pred_ee.push_back(S->numbers[reg][Elec].npf + S->numbers[reg][Muon].nfp);
		nff_pred_ee.push_back(S->numbers[reg][Elec].nff);
		npp_pred_em.push_back(S->numbers[reg][ElMu].npp);
		npf_pred_em.push_back(S->numbers[reg][ElMu].npf);
		nfp_pred_em.push_back(S->numbers[reg][ElMu].nfp);
		nff_pred_em.push_back(S->numbers[reg][ElMu].nff);


		npp_pred_mm_e1.push_back(FR->getMMNppEStat());
		npf_pred_mm_e1.push_back(FR->getMMNpfEStat());
		nff_pred_mm_e1.push_back(FR->getMMNffEStat());
		nF_pred_mm_e1 .push_back(FR->getMMTotEStat());

		npp_pred_ee_e1.push_back(FR->getEENppEStat());
		npf_pred_ee_e1.push_back(FR->getEENpfEStat());
		nff_pred_ee_e1.push_back(FR->getEENffEStat());
		nF_pred_ee_e1 .push_back(FR->getEETotEStat());

		npp_pred_em_e1.push_back(FR->getEMNppEStat());
		npf_pred_em_e1.push_back(FR->getEMNpfEStat());
		nfp_pred_em_e1.push_back(FR->getEMNfpEStat());
		nff_pred_em_e1.push_back(FR->getEMNffEStat());
		nF_pred_em_e1 .push_back(FR->getEMTotEStat());

		// npp_pred_mm   .push_back(FR->getMMNpp());
		// npp_pred_mm_e1.push_back(FR->getMMNppEStat());
		// npf_pred_mm   .push_back(FR->getMMNpf());
		// npf_pred_mm_e1.push_back(FR->getMMNpfEStat());
		// nff_pred_mm   .push_back(FR->getMMNff());
		// nff_pred_mm_e1.push_back(FR->getMMNffEStat());
		// nF_pred_mm_e1 .push_back(FR->getMMTotEStat());
		// 
		// npp_pred_ee   .push_back(FR->getEENpp());
		// npp_pred_ee_e1.push_back(FR->getEENppEStat());
		// npf_pred_ee   .push_back(FR->getEENpf());
		// npf_pred_ee_e1.push_back(FR->getEENpfEStat());
		// nff_pred_ee   .push_back(FR->getEENff());
		// nff_pred_ee_e1.push_back(FR->getEENffEStat());
		// nF_pred_ee_e1 .push_back(FR->getEETotEStat());
		// 
		// npp_pred_em   .push_back(FR->getEMNpp());
		// npp_pred_em_e1.push_back(FR->getEMNppEStat());
		// npf_pred_em   .push_back(FR->getEMNpf());
		// npf_pred_em_e1.push_back(FR->getEMNpfEStat());
		// nfp_pred_em   .push_back(FR->getEMNfp());
		// nfp_pred_em_e1.push_back(FR->getEMNfpEStat());
		// nff_pred_em   .push_back(FR->getEMNff());
		// nff_pred_em_e1.push_back(FR->getEMNffEStat());
		// nF_pred_em_e1 .push_back(FR->getEMTotEStat());
	}
	
	// Charge MisID Predictions
	// Abbreviations
	float fb  = gEChMisIDB;
	float fbE = gEChMisIDB_E;
	float fe  = gEChMisIDE;
	float feE = gEChMisIDE_E;
	
	vector<float> ntt_cm_ee, ntt_cm_em;
	for(size_t i = 0; i < nsamples; ++i){
		ntt_cm_ee.push_back(2*fb*ntt_os_BB_ee[i] + 2*fe*ntt_os_EE_ee[i] + (fb+fe)*ntt_os_EB_ee[i]);
		ntt_cm_em.push_back(  fb*ntt_os_BB_em[i] +   fe*ntt_os_EE_em[i]);
	}

	// Sums
	float ntt_sum_mm(0.), ntl_sum_mm(0.), nll_sum_mm(0.);
	float ntt_sum_em(0.), ntl_sum_em(0.), nlt_sum_em(0.), nll_sum_em(0.);
	float ntt_sum_ee(0.), ntl_sum_ee(0.), nll_sum_ee(0.);

	float npp_sum_mm(0.), npf_sum_mm(0.), nff_sum_mm(0.);
	float npp_sum_em(0.), npf_sum_em(0.), nfp_sum_em(0.), nff_sum_em(0.);
	float npp_sum_ee(0.), npf_sum_ee(0.), nff_sum_ee(0.);

	float npp_pred_sum_mm(0.), npf_pred_sum_mm(0.), nff_pred_sum_mm(0.);
	float npp_pred_sum_em(0.), npf_pred_sum_em(0.), nfp_pred_sum_em(0.), nff_pred_sum_em(0.);
	float npp_pred_sum_ee(0.), npf_pred_sum_ee(0.), nff_pred_sum_ee(0.);

	float ntt_cm_sum_ee(0.), ntt_cm_sum_em(0.);

	float npp_cm_sum_ee(0.),    npp_cm_sum_em(0.);
	float npp_tt_cm_sum_ee(0.), npp_tt_cm_sum_em(0.);

	float ntt_rare_mm(0.), ntt_rare_em(0.), ntt_rare_ee(0.);
	for(size_t i = 0; i < fMCRareSM.size(); ++i){
		Sample *S = fSamples[fMCRareSM[i]];
		float scale = fLumiNorm/S->lumi;
		ntt_rare_mm += scale*S->numbers[reg][Muon].nt2;
		ntt_rare_em += scale*S->numbers[reg][ElMu].nt2;
		ntt_rare_ee += scale*S->numbers[reg][Elec].nt2;
	}



	// Squared errors
	float npp_pred_sum_mm_e1(0.), npf_pred_sum_mm_e1(0.), nff_pred_sum_mm_e1(0.);
	float npp_pred_sum_em_e1(0.), npf_pred_sum_em_e1(0.), nfp_pred_sum_em_e1(0.), nff_pred_sum_em_e1(0.);
	float npp_pred_sum_ee_e1(0.), npf_pred_sum_ee_e1(0.), nff_pred_sum_ee_e1(0.);

	// Combined stat. errors
	float nF_pred_sum_mm_e1(0.), nF_pred_sum_em_e1(0.), nF_pred_sum_ee_e1(0.);

	for(size_t i = 0; i < nsamples; ++i){
		ntt_sum_mm         += scales[i] * ntt_mm[i];
		ntl_sum_mm         += scales[i] * ntl_mm[i];
		nll_sum_mm         += scales[i] * nll_mm[i];
		npp_sum_mm         += scales[i] * npp_mm[i];
		npf_sum_mm         += scales[i] * (npf_mm[i]+nfp_mm[i]);
		nff_sum_mm         += scales[i] * nff_mm[i];
		npp_pred_sum_mm    += scales[i] * npp_pred_mm[i];
		npf_pred_sum_mm    += scales[i] * npf_pred_mm[i];
		nff_pred_sum_mm    += scales[i] * nff_pred_mm[i];
		npp_pred_sum_mm_e1 += scales[i]*scales[i] * npp_pred_mm_e1[i]*npp_pred_mm_e1[i];
		npf_pred_sum_mm_e1 += scales[i]*scales[i] * npf_pred_mm_e1[i]*npf_pred_mm_e1[i];
		nff_pred_sum_mm_e1 += scales[i]*scales[i] * nff_pred_mm_e1[i]*nff_pred_mm_e1[i];
		nF_pred_sum_mm_e1  += scales[i]*scales[i] * nF_pred_mm_e1[i]*nF_pred_mm_e1[i];

		ntt_sum_ee         += scales[i] * ntt_ee[i];
		ntl_sum_ee         += scales[i] * ntl_ee[i];
		nll_sum_ee         += scales[i] * nll_ee[i];
		npp_sum_ee         += scales[i] * npp_ee[i];
		npf_sum_ee         += scales[i] * (npf_ee[i]+nfp_ee[i]);
		nff_sum_ee         += scales[i] * nff_ee[i];
		npp_pred_sum_ee    += scales[i] * npp_pred_ee[i];
		npf_pred_sum_ee    += scales[i] * npf_pred_ee[i];
		nff_pred_sum_ee    += scales[i] * nff_pred_ee[i];
		npp_pred_sum_ee_e1 += scales[i]*scales[i] * npp_pred_ee_e1[i]*npp_pred_ee_e1[i];
		npf_pred_sum_ee_e1 += scales[i]*scales[i] * npf_pred_ee_e1[i]*npf_pred_ee_e1[i];
		nff_pred_sum_ee_e1 += scales[i]*scales[i] * nff_pred_ee_e1[i]*nff_pred_ee_e1[i];
		nF_pred_sum_ee_e1  += scales[i]*scales[i] * nF_pred_ee_e1[i]*nF_pred_ee_e1[i];

		ntt_sum_em         += scales[i] * ntt_em[i];
		ntl_sum_em         += scales[i] * ntl_em[i];
		nlt_sum_em         += scales[i] * nlt_em[i];
		nll_sum_em         += scales[i] * nll_em[i];
		npp_sum_em         += scales[i] * npp_em[i];
		npf_sum_em         += scales[i] * npf_em[i];
		nfp_sum_em         += scales[i] * nfp_em[i];
		nff_sum_em         += scales[i] * nff_em[i];
		npp_pred_sum_em    += scales[i] * npp_pred_em[i];
		npf_pred_sum_em    += scales[i] * npf_pred_em[i];
		nfp_pred_sum_em    += scales[i] * nfp_pred_em[i];
		nff_pred_sum_em    += scales[i] * nff_pred_em[i];
		npp_pred_sum_em_e1 += scales[i]*scales[i] * npp_pred_em_e1[i]*npp_pred_em_e1[i];
		npf_pred_sum_em_e1 += scales[i]*scales[i] * npf_pred_em_e1[i]*npf_pred_em_e1[i];
		nfp_pred_sum_em_e1 += scales[i]*scales[i] * nfp_pred_em_e1[i]*nfp_pred_em_e1[i];
		nff_pred_sum_em_e1 += scales[i]*scales[i] * nff_pred_em_e1[i]*nff_pred_em_e1[i];
		nF_pred_sum_em_e1  += scales[i]*scales[i] * nF_pred_em_e1[i]*nF_pred_em_e1[i];
		
		ntt_cm_sum_ee   += scales[i] * ntt_cm_ee[i];
		ntt_cm_sum_em   += scales[i] * ntt_cm_em[i];

		npp_cm_sum_ee    += scales[i] * npp_cm_ee[i];
		npp_cm_sum_em    += scales[i] * npp_cm_em[i];
		npp_tt_cm_sum_ee += scales[i] * npp_tt_cm_ee[i];
		npp_tt_cm_sum_em += scales[i] * npp_tt_cm_em[i];
	}


	///////////////////////////////////////////////////////////////////////////////////
	// PRINTOUT ///////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////
	OUT << "-------------------------------------------------------------------------------------------------" << endl;
	OUT << "         RATIOS  ||    Mu-fRatio     |    Mu-pRatio     ||    El-fRatio     |    El-pRatio     ||" << endl;
	OUT << "-------------------------------------------------------------------------------------------------" << endl;
	OUT << setw(16) << "           MC    ||";
	OUT << setw(7)  << setprecision(2) << mf << " ± " << setw(7) << setprecision(2) << mf_e << " |";
	OUT << setw(7)  << setprecision(2) << mp << " ± " << setw(7) << setprecision(2) << mp_e << " ||";
	OUT << setw(7)  << setprecision(2) << ef << " ± " << setw(7) << setprecision(2) << ef_e << " |";
	OUT << setw(7)  << setprecision(2) << ep << " ± " << setw(7) << setprecision(2) << ep_e << " ||";
	OUT << endl;
	OUT << "-------------------------------------------------------------------------------------------------" << endl << endl;

	OUT << "==========================================================================================================================" << endl;
	OUT << "                 ||            Mu/Mu            ||                   E/Mu                ||             E/E             ||" << endl;
	OUT << "          YIELDS ||   Ntt   |   Ntl   |   Nll   ||   Ntt   |   Ntl   |   Nlt   |   Nll   ||   Ntt   |   Ntl   |   Nll   ||" << endl;
	OUT << "--------------------------------------------------------------------------------------------------------------------------" << endl;
	for(size_t i = 0; i < nsamples; ++i){
		OUT << setw(16) << names[i] << " || ";
		OUT << setw(7)  << Form("%6.3f", scales[i]*ntt_mm[i]) << " | ";
		OUT << setw(7)  << Form("%6.3f", scales[i]*ntl_mm[i]) << " | ";
		OUT << setw(7)  << Form("%6.3f", scales[i]*nll_mm[i]) << " || ";
		OUT << setw(7)  << Form("%6.3f", scales[i]*ntt_em[i]) << " | ";
		OUT << setw(7)  << Form("%6.3f", scales[i]*ntl_em[i]) << " | ";
		OUT << setw(7)  << Form("%6.3f", scales[i]*nlt_em[i]) << " | ";
		OUT << setw(7)  << Form("%6.3f", scales[i]*nll_em[i]) << " || ";
		OUT << setw(7)  << Form("%6.3f", scales[i]*ntt_ee[i]) << " | ";
		OUT << setw(7)  << Form("%6.3f", scales[i]*ntl_ee[i]) << " | ";
		OUT << setw(7)  << Form("%6.3f", scales[i]*nll_ee[i]) << " || ";
		OUT << endl;
	}	
	OUT << "--------------------------------------------------------------------------------------------------------------------------" << endl;
	OUT << setw(16) << "Sum"  << " || ";
	OUT << setw(7) << Form("%6.3f", ntt_sum_mm) << " | ";
	OUT << setw(7) << Form("%6.3f", ntl_sum_mm) << " | ";
	OUT << setw(7) << Form("%6.3f", nll_sum_mm) << " || ";
	OUT << setw(7) << Form("%6.3f", ntt_sum_em) << " | ";
	OUT << setw(7) << Form("%6.3f", ntl_sum_em) << " | ";
	OUT << setw(7) << Form("%6.3f", nlt_sum_em) << " | ";
	OUT << setw(7) << Form("%6.3f", nll_sum_em) << " || ";
	OUT << setw(7) << Form("%6.3f", ntt_sum_ee) << " | ";
	OUT << setw(7) << Form("%6.3f", ntl_sum_ee) << " | ";
	OUT << setw(7) << Form("%6.3f", nll_sum_ee) << " || ";
	OUT << endl;
	OUT << setw(16) << "Channels sum"  << " || ";
	OUT << Form("                     %6.3f ||                               %6.3f ||                      %6.3f || ",
	ntt_sum_mm+ntl_sum_mm+nll_sum_mm, ntt_sum_em+ntl_sum_em+nlt_sum_em+nll_sum_em, ntt_sum_ee+ntl_sum_ee+nll_sum_ee) << endl;
	OUT << "==========================================================================================================================" << endl;
	OUT << endl;
	OUT << "==========================================================================================================================" << endl;
	OUT << "                 ||            Mu/Mu            ||                   E/Mu                ||             E/E             ||" << endl;
	OUT << "           TRUTH ||   Npp   |   Nfp   |   Nff   ||   Npp   |   Npf   |   Nfp   |   Nff   ||   Npp   |   Npf   |   Nff   ||" << endl;
	OUT << "--------------------------------------------------------------------------------------------------------------------------" << endl;
	for(size_t i = 0; i < nsamples; ++i){
		OUT << setw(16) << names[i] << " || ";
		OUT << setw(7)  << Form("%6.3f", scales[i]*npp_mm[i]) << " | ";
		OUT << setw(7)  << Form("%6.3f", scales[i]*(npf_mm[i]+nfp_mm[i])) << " | ";
		OUT << setw(7)  << Form("%6.3f", scales[i]*nff_mm[i]) << " || ";
		OUT << setw(7)  << Form("%6.3f", scales[i]*npp_em[i]) << " | ";
		OUT << setw(7)  << Form("%6.3f", scales[i]*npf_em[i]) << " | ";
		OUT << setw(7)  << Form("%6.3f", scales[i]*nfp_em[i]) << " | ";
		OUT << setw(7)  << Form("%6.3f", scales[i]*nff_em[i]) << " || ";
		OUT << setw(7)  << Form("%6.3f", scales[i]*npp_ee[i]) << " | ";
		OUT << setw(7)  << Form("%6.3f", scales[i]*npf_ee[i]+nfp_ee[i]) << " | ";
		OUT << setw(7)  << Form("%6.3f", scales[i]*nff_ee[i]) << " || ";
		OUT << endl;
	}
	OUT << "--------------------------------------------------------------------------------------------------------------------------" << endl;
	OUT << setw(16) << "Sum"  << " || ";
	OUT << setw(7) << Form("%6.3f", npp_sum_mm) << " | ";
	OUT << setw(7) << Form("%6.3f", npf_sum_mm) << " | ";
	OUT << setw(7) << Form("%6.3f", nff_sum_mm) << " || ";
	OUT << setw(7) << Form("%6.3f", npp_sum_em) << " | ";
	OUT << setw(7) << Form("%6.3f", npf_sum_em) << " | ";
	OUT << setw(7) << Form("%6.3f", nfp_sum_em) << " | ";
	OUT << setw(7) << Form("%6.3f", nff_sum_em) << " || ";
	OUT << setw(7) << Form("%6.3f", npp_sum_ee) << " | ";
	OUT << setw(7) << Form("%6.3f", npf_sum_ee) << " | ";
	OUT << setw(7) << Form("%6.3f", nff_sum_ee) << " || ";
	OUT << endl;
	OUT << setw(16) << "Channels sum"  << " || ";
	OUT << Form("                     %6.3f ||                               %6.3f ||                      %6.3f || ",
	npp_sum_mm+npf_sum_mm+nff_sum_mm, npp_sum_em+npf_sum_em+nfp_sum_em+nff_sum_em, npp_sum_ee+npf_sum_ee+nff_sum_ee) << endl;
	OUT << "==========================================================================================================================" << endl;
	OUT << endl;

	OUT << "==========================================================================================================================" << endl;
	OUT << "                 ||            Mu/Mu            ||                   E/Mu                ||             E/E             ||" << endl;
	OUT << "     TRUTH IN TT ||   Npp   |   Nfp   |   Nff   ||   Npp   |   Npf   |   Nfp   |   Nff   ||   Npp   |   Npf   |   Nff   ||" << endl;
	OUT << "--------------------------------------------------------------------------------------------------------------------------" << endl;
	for(size_t i = 0; i < nsamples; ++i){
		OUT << setw(16) << names[i] << " || ";
		OUT << setw(7)  << Form("%6.3f", mp*mp*scales[i]*npp_mm[i]) << " | ";
		OUT << setw(7)  << Form("%6.3f", mp*mf*scales[i]*(npf_mm[i]+nfp_mm[i])) << " | ";
		OUT << setw(7)  << Form("%6.3f", mf*mf*scales[i]*nff_mm[i]) << " || ";
		OUT << setw(7)  << Form("%6.3f", mp*ep*scales[i]*npp_em[i]) << " | ";
		OUT << setw(7)  << Form("%6.3f", mp*ef*scales[i]*npf_em[i]) << " | ";
		OUT << setw(7)  << Form("%6.3f", mf*ep*scales[i]*nfp_em[i]) << " | ";
		OUT << setw(7)  << Form("%6.3f", mf*ef*scales[i]*nff_em[i]) << " || ";
		OUT << setw(7)  << Form("%6.3f", ep*ep*scales[i]*npp_ee[i]) << " | ";
		OUT << setw(7)  << Form("%6.3f", ep*ef*scales[i]*(npf_ee[i]+nfp_ee[i])) << " | ";
		OUT << setw(7)  << Form("%6.3f", ef*ef*scales[i]*nff_ee[i]) << " || ";
		OUT << endl;
	}
	OUT << "--------------------------------------------------------------------------------------------------------------------------" << endl;
	OUT << setw(16) << "Npf Sum"  << " || ";
	OUT << setw(7) << Form("%6.3f", mp*mp*npp_sum_mm) << " | ";
	OUT << setw(7) << Form("%6.3f", mp*mf*npf_sum_mm) << " | ";
	OUT << setw(7) << Form("%6.3f", mf*mf*nff_sum_mm) << " || ";
	OUT << setw(7) << Form("%6.3f", mp*ep*npp_sum_em) << " | ";
	OUT << setw(7) << Form("%6.3f", mp*ef*npf_sum_em) << " | ";
	OUT << setw(7) << Form("%6.3f", mf*ep*nfp_sum_em) << " | ";
	OUT << setw(7) << Form("%6.3f", mf*ef*nff_sum_em) << " || ";
	OUT << setw(7) << Form("%6.3f", ep*ep*npp_sum_ee) << " | ";
	OUT << setw(7) << Form("%6.3f", ep*ef*npf_sum_ee) << " | ";
	OUT << setw(7) << Form("%6.3f", ef*ef*nff_sum_ee) << " || ";
	OUT << endl;
	OUT << "==========================================================================================================================" << endl;
	OUT << endl;

	OUT << "==========================================================================================================================" << endl;
	OUT << "                 ||            Mu/Mu            ||                   E/Mu                ||             E/E             ||" << endl;
	OUT << "     PRED. IN TT ||   Npp   |   Nfp   |   Nff   ||   Npp   |   Npf   |   Nfp   |   Nff   ||   Npp   |   Npf   |   Nff   ||" << endl;
	OUT << "--------------------------------------------------------------------------------------------------------------------------" << endl;
	for(size_t i = 0; i < nsamples; ++i){
		OUT << setw(16) << names[i] << " || ";
		OUT << setw(7)  << Form("%6.3f", scales[i]*npp_pred_mm[i]) << " | ";
		OUT << setw(7)  << Form("%6.3f", scales[i]*npf_pred_mm[i]) << " | ";
		OUT << setw(7)  << Form("%6.3f", scales[i]*nff_pred_mm[i]) << " || ";
		OUT << setw(7)  << Form("%6.3f", scales[i]*npp_pred_em[i]) << " | ";
		OUT << setw(7)  << Form("%6.3f", scales[i]*npf_pred_em[i]) << " | ";
		OUT << setw(7)  << Form("%6.3f", scales[i]*nfp_pred_em[i]) << " | ";
		OUT << setw(7)  << Form("%6.3f", scales[i]*nff_pred_em[i]) << " || ";
		OUT << setw(7)  << Form("%6.3f", scales[i]*npp_pred_ee[i]) << " | ";
		OUT << setw(7)  << Form("%6.3f", scales[i]*npf_pred_ee[i]) << " | ";
		OUT << setw(7)  << Form("%6.3f", scales[i]*nff_pred_ee[i]) << " || ";
		OUT << endl;
	}
	OUT << "--------------------------------------------------------------------------------------------------------------------------" << endl;
	OUT << setw(16) << "Pred Sum"  << " || ";
	OUT << setw(7) << Form("%6.3f", npp_pred_sum_mm) << " | ";
	OUT << setw(7) << Form("%6.3f", npf_pred_sum_mm) << " | ";
	OUT << setw(7) << Form("%6.3f", nff_pred_sum_mm) << " || ";
	OUT << setw(7) << Form("%6.3f", npp_pred_sum_em) << " | ";
	OUT << setw(7) << Form("%6.3f", npf_pred_sum_em) << " | ";
	OUT << setw(7) << Form("%6.3f", nfp_pred_sum_em) << " | ";
	OUT << setw(7) << Form("%6.3f", nff_pred_sum_em) << " || ";
	OUT << setw(7) << Form("%6.3f", npp_pred_sum_ee) << " | ";
	OUT << setw(7) << Form("%6.3f", npf_pred_sum_ee) << " | ";
	OUT << setw(7) << Form("%6.3f", nff_pred_sum_ee) << " || ";
	OUT << endl;
	OUT << "==========================================================================================================================" << endl;
	OUT << endl;

	OUT << "=========================================================" << endl;
	OUT << "                 ||       E/Mu      ||       E/E       ||" << endl;
	OUT << "    CHARGE MISID ||  Pred  |  Truth ||  Pred  |  Truth ||" << endl;
	OUT << "---------------------------------------------------------" << endl;
	for(size_t i = 0; i < nsamples; ++i){
		OUT << setw(16) << names[i] << " || ";
		OUT << setw(7)  << Form("%6.3f | %6.3f || ", scales[i]*ntt_cm_em[i], scales[i]*npp_tt_cm_em[i]);
		OUT << setw(7)  << Form("%6.3f | %6.3f || ", scales[i]*ntt_cm_ee[i], scales[i]*npp_tt_cm_ee[i]);
		OUT << endl;
	}	
	OUT << "---------------------------------------------------------" << endl;
	OUT << setw(16) << "Sum"  << " || ";
	OUT << setw(7)  << Form("%6.3f | %6.3f || ", ntt_cm_sum_em, npp_tt_cm_sum_em);
	OUT << setw(7)  << Form("%6.3f | %6.3f || ", ntt_cm_sum_ee, npp_tt_cm_sum_ee) << endl;
	OUT << "=========================================================" << endl;
	OUT << endl;

	OUT << "==========================================================================================================================" << endl;
	OUT << "                 ||            Mu/Mu            ||                   E/Mu                ||             E/E             ||" << endl;
	OUT << "     PRED. IN TT ||   Npp   |   Nfp   |   Nff   ||   Npp   |   Npf   |   Nfp   |   Nff   ||   Npp   |   Npf   |   Nff   ||" << endl;
	OUT << "--------------------------------------------------------------------------------------------------------------------------" << endl;
	OUT << setw(16) << "Npf Truth"  << " || ";
	OUT << setw(7) << Form("%6.3f", mp*mp*npp_sum_mm) << " | ";
	OUT << setw(7) << Form("%6.3f", mp*mf*npf_sum_mm) << " | ";
	OUT << setw(7) << Form("%6.3f", mf*mf*nff_sum_mm) << " || ";
	OUT << setw(7) << Form("%6.3f", mp*ep*npp_sum_em) << " | ";
	OUT << setw(7) << Form("%6.3f", mp*ef*npf_sum_em) << " | ";
	OUT << setw(7) << Form("%6.3f", mf*ep*nfp_sum_em) << " | ";
	OUT << setw(7) << Form("%6.3f", mf*ef*nff_sum_em) << " || ";
	OUT << setw(7) << Form("%6.3f", ep*ep*npp_sum_ee) << " | ";
	OUT << setw(7) << Form("%6.3f", ep*ef*npf_sum_ee) << " | ";
	OUT << setw(7) << Form("%6.3f", ef*ef*nff_sum_ee) << " || ";
	OUT << endl;
	OUT << "--------------------------------------------------------------------------------------------------------------------------" << endl;
	OUT << setw(16) << "FR Prediction"  << " || ";
	OUT << setw(7) << Form("%6.3f", npp_pred_sum_mm) << " | ";
	OUT << setw(7) << Form("%6.3f", npf_pred_sum_mm) << " | ";
	OUT << setw(7) << Form("%6.3f", nff_pred_sum_mm) << " || ";
	OUT << setw(7) << Form("%6.3f", npp_pred_sum_em) << " | ";
	OUT << setw(7) << Form("%6.3f", npf_pred_sum_em) << " | ";
	OUT << setw(7) << Form("%6.3f", nfp_pred_sum_em) << " | ";
	OUT << setw(7) << Form("%6.3f", nff_pred_sum_em) << " || ";
	OUT << setw(7) << Form("%6.3f", npp_pred_sum_ee) << " | ";
	OUT << setw(7) << Form("%6.3f", npf_pred_sum_ee) << " | ";
	OUT << setw(7) << Form("%6.3f", nff_pred_sum_ee) << " || ";
	OUT << endl;
	OUT << "==========================================================================================================================" << endl;
	OUT << setw(16) << "Pred. Fakes"  << " || ";
	OUT << setw(7) << Form("%6.3f", npf_pred_sum_mm+nff_pred_sum_mm) << " |                   || ";
	OUT << setw(7) << Form("%6.3f", npf_pred_sum_em+nfp_pred_sum_em+nff_pred_sum_em) << " |                             || ";
	OUT << setw(7) << Form("%6.3f", npf_pred_sum_ee+nff_pred_sum_ee) << " |                   || ";
	OUT << endl;
	OUT << setw(16) << "Pred. Ch-MID"  << " ||         |                   || ";
	OUT << setw(7) << Form("%6.3f", ntt_cm_sum_em) << " |                             || ";
	OUT << setw(7) << Form("%6.3f", ntt_cm_sum_ee) << " |                   || ";
	OUT << endl;
	OUT << setw(16) << "Pred. Rare SM"  << " || ";
	OUT << setw(7) << Form("%6.3f", ntt_rare_mm) << " |                   || ";
	OUT << setw(7) << Form("%6.3f", ntt_rare_em) << " |                             || ";
	OUT << setw(7) << Form("%6.3f", ntt_rare_ee) << " |                   || ";
	OUT << endl;
	OUT << "--------------------------------------------------------------------------------------------------------------------------" << endl;
	OUT << setw(16) << "Total BG Pred."  << " || ";
	float tot_bg_mm = npf_pred_sum_mm                +nff_pred_sum_mm              +ntt_rare_mm;
	float tot_bg_em = npf_pred_sum_em+nfp_pred_sum_em+nff_pred_sum_em+ntt_cm_sum_em+ntt_rare_em;
	float tot_bg_ee = npf_pred_sum_ee                +nff_pred_sum_ee+ntt_cm_sum_ee+ntt_rare_ee;
	float tot_bg = tot_bg_mm + tot_bg_em + tot_bg_ee;
	float ntt_sum = ntt_sum_mm+ntt_sum_em+ntt_sum_ee;
	OUT << setw(7) << Form("%6.3f", tot_bg_mm) << " | P/O  | P-O/P |    || ";
	OUT << setw(7) << Form("%6.3f", tot_bg_em) << " | P/O  | P-O/P |              || ";
	OUT << setw(7) << Form("%6.3f", tot_bg_ee) << " | P/O  | P-O/P |    || ";
	OUT << endl;
	OUT << "--------------------------------------------------------------------------------------------------------------------------" << endl;
	OUT << setw(16) << "Observed"  << " || ";
	OUT << setw(7) << Form(" %6.3f | %4.2f | %5.2f |", ntt_sum_mm, tot_bg_mm/ntt_sum_mm, (tot_bg_mm-ntt_sum_mm)/tot_bg_mm) << "    || ";
	OUT << setw(7) << Form(" %6.3f | %4.2f | %5.2f |", ntt_sum_em, tot_bg_em/ntt_sum_em, (tot_bg_em-ntt_sum_em)/tot_bg_em) << "              || ";
	OUT << setw(7) << Form(" %6.3f | %4.2f | %5.2f |", ntt_sum_ee, tot_bg_ee/ntt_sum_ee, (tot_bg_ee-ntt_sum_ee)/tot_bg_ee) << "    || ";
	OUT << endl;
	OUT << "==========================================================================================================================" << endl;
	OUT << Form("     Comb. Fakes || %6.3f |              ||",
	npf_pred_sum_mm+nff_pred_sum_mm + 
	npf_pred_sum_em+nfp_pred_sum_em+nff_pred_sum_em +
	npf_pred_sum_ee+nff_pred_sum_ee) << endl;
	OUT << Form("    Comb. Ch-MID || %6.3f |              ||", ntt_cm_sum_em+ntt_cm_sum_ee) << endl;
	OUT << Form("   Comb. Rare SM || %6.3f |              ||", ntt_rare_mm+ntt_rare_em+ntt_rare_ee) << endl;
	OUT << "---------------------------------------------" << endl;
	OUT << Form("  Total BG Comb. || %6.3f | P/O  | P-O/P ||", tot_bg) << endl;
	OUT << "---------------------------------------------" << endl;
	OUT << Form("        Observed || %6.3f | %4.2f | %5.2f ||", ntt_sum, tot_bg/ntt_sum, (tot_bg-ntt_sum)/tot_bg) << endl;
	OUT << "=============================================" << endl;
	OUT << endl;

	

	// OUT << "===========================================================================================================================================" << endl;
	// OUT << "  All predictions (Npp / Npf / (Nfp) / Nff):                                                                                              |" << endl;
	// for(size_t i = 0; i < nsamples; ++i){
	// 	OUT << setw(16) << left << names[i];
	// 	OUT << Form("  MM || %9.3f ± %7.3f (stat) | %9.3f ± %7.3f (stat) | %9.3f ± %7.3f (stat) |                            |",
	// 	npp_pred_mm[i], npp_pred_mm_e1[i], npf_pred_mm[i], npf_pred_mm_e1[i], nff_pred_mm[i], nff_pred_mm_e1[i]) << endl;
	// 	OUT << " scale = " << setw(7) << setprecision(2) << scales[i];
	// 	OUT << Form("  EM || %9.3f ± %7.3f (stat) | %9.3f ± %7.3f (stat) | %9.3f ± %7.3f (stat) | %9.3f ± %7.3f (stat) |",
	// 	npp_pred_em[i], npp_pred_em_e1[i], npf_pred_em[i], npf_pred_em_e1[i], nfp_pred_em[i], nfp_pred_em_e1[i], nff_pred_mm[i], nff_pred_em_e1[i]) << endl;
	// 	OUT << Form("                  EE || %9.3f ± %7.3f (stat) | %9.3f ± %7.3f (stat) | %9.3f ± %7.3f (stat) |                            |",
	// 	npp_pred_ee[i], npp_pred_ee_e1[i], npf_pred_ee[i], npf_pred_ee_e1[i], nff_pred_ee[i], nff_pred_ee_e1[i]) << endl;
	// }
	// OUT << "===========================================================================================================================================" << endl;
	// OUT << endl;

	// OUT << "==========================================================================================================================" << endl;
	// OUT << "  PREDICTIONS (in tt window)" << endl;
	// OUT << "--------------------------------------------------------------" << endl;
	// OUT << " Mu/Mu Channel:" << endl;
	// OUT << "  Npp*pp:        " <<  Form("%6.3f", npp_pred_sum_mm) << " ± " << Form("%6.3f", sqrt(npp_pred_sum_mm_e1)) << " (stat)" << endl;
	// OUT << "  Nfp*fp:        " <<  Form("%6.3f", npf_pred_sum_mm) << " ± " << Form("%6.3f", sqrt(npf_pred_sum_mm_e1)) << " (stat)" << endl;
	// OUT << "  Nff*ff:        " <<  Form("%6.3f", nff_pred_sum_mm) << " ± " << Form("%6.3f", sqrt(nff_pred_sum_mm_e1)) << " (stat)" << endl;
	// OUT << "  Total fakes:   " <<  Form("%6.3f", npf_pred_sum_mm+nff_pred_sum_mm) << " ± " << Form("%6.3f", sqrt(nF_pred_sum_mm_e1)) << " (stat)" << endl;
	// OUT << "--------------------------------------------------------------" << endl;
	// OUT << " E/Mu Channel:" << endl;
	// OUT << "  Npp*pp:        " <<  Form("%6.3f", npp_pred_sum_em) << " ± " << Form("%6.3f", sqrt(npp_pred_sum_em_e1)) << " (stat)" << endl;
	// OUT << "  Nfp*fp:        " <<  Form("%6.3f", npf_pred_sum_em) << " ± " << Form("%6.3f", sqrt(npf_pred_sum_em_e1)) << " (stat)" << endl;
	// OUT << "  Npf*pf:        " <<  Form("%6.3f", nfp_pred_sum_em) << " ± " << Form("%6.3f", sqrt(nfp_pred_sum_em_e1)) << " (stat)" << endl;
	// OUT << "  Nff*ff:        " <<  Form("%6.3f", nff_pred_sum_em) << " ± " << Form("%6.3f", sqrt(nff_pred_sum_em_e1)) << " (stat)" << endl;
	// OUT << "  Total fakes:   " <<  Form("%6.3f", npf_pred_sum_em+nfp_pred_sum_em+nff_pred_sum_em) << " ± " << Form("%6.3f", sqrt(nF_pred_sum_em_e1)) << " (stat)" << endl;
	// OUT << "--------------------------------------------------------------" << endl;
	// OUT << " E/E Channel:" << endl;
	// OUT << "  Npp*pp:        " <<  Form("%6.3f", npp_pred_sum_ee) << " ± " << Form("%6.3f", sqrt(npp_pred_sum_ee_e1)) << " (stat)" << endl;
	// OUT << "  Nfp*fp:        " <<  Form("%6.3f", npf_pred_sum_ee) << " ± " << Form("%6.3f", sqrt(npf_pred_sum_ee_e1)) << " (stat)" << endl;
	// OUT << "  Nff*ff:        " <<  Form("%6.3f", nff_pred_sum_ee) << " ± " << Form("%6.3f", sqrt(nff_pred_sum_ee_e1)) << " (stat)" << endl;
	// OUT << "  Total fakes:   " <<  Form("%6.3f", npf_pred_sum_ee+nff_pred_sum_ee) << " ± " << Form("%6.3f", sqrt(nF_pred_sum_ee_e1)) << " (stat)" << endl;
	// OUT << "==========================================================================================================================" << endl;
	// OUT << "/////////////////////////////////////////////////////////////////////////////" << endl;
	
	OUT.close();
	delete FR;
}
void SSDLPlotter::makeTTbarClosure(){
	TString filename = "TTbarClosure.txt";
	ofstream OUT(fOutputDir + filename.Data(), ios::trunc);

	///////////////////////////////////////////////////////////////////////////////////
	// RATIOS /////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////
	float mufratio_allmc(0.), mufratio_allmc_e(0.);
	float mupratio_allmc(0.), mupratio_allmc_e(0.);
	float elfratio_allmc(0.), elfratio_allmc_e(0.);
	float elpratio_allmc(0.), elpratio_allmc_e(0.);

	vector<int> ttjets; ttjets.push_back(TTJets);
	// calculateRatio(ttjets, Muon,     SigSup, mufratio_allmc, mufratio_allmc_e);
	// calculateRatio(ttjets, Muon,     ZDecay, mupratio_allmc, mupratio_allmc_e);
	// calculateRatio(ttjets, Elec, SigSup, elfratio_allmc, elfratio_allmc_e);
	// calculateRatio(ttjets, Elec, ZDecay, elpratio_allmc, elpratio_allmc_e);
	calculateRatio(fMCBGMuEnr, Muon, SigSup, mufratio_allmc, mufratio_allmc_e);
	calculateRatio(fMCBGMuEnr, Muon, ZDecay, mupratio_allmc, mupratio_allmc_e);
	calculateRatio(fMCBG,      Elec, SigSup, elfratio_allmc, elfratio_allmc_e);
	calculateRatio(fMCBG,      Elec, ZDecay, elpratio_allmc, elpratio_allmc_e);

	///////////////////////////////////////////////////////////////////////////////////
	// OBSERVATIONS ///////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////
	float nt_mm(0.), nl_mm(0.),  np_mm(0.), nf_mm(0.);
	float nt_em(0.), nl_em(0.),  np_em(0.), nf_em(0.);
	float nt_me(0.), nl_me(0.),  np_me(0.), nf_me(0.);
	float nt_ee(0.), nl_ee(0.),  np_ee(0.), nf_ee(0.);

	Sample *S = fSamples[TTJets];
	TTree *tree = S->getTree();
	fCurrentSample = TTJets;

	// Event loop
	tree->ResetBranchAddresses();
	// Init(tree);
	if(S->datamc == 0) Init(tree);
	else InitMC(tree);

	if (fChain == 0) return;
	Long64_t nentries = fChain->GetEntriesFast();
	Long64_t nbytes = 0, nb = 0;
	for (Long64_t jentry=0; jentry<nentries;jentry++) {
		printProgress(jentry, nentries, S->name);

		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		nb = fChain->GetEntry(jentry);   nbytes += nb;

		int ind1(-1), ind2(-1);
		
		// MuMu
		fCurrentChannel = Muon;
		if(isSSLLMuEvent(ind1, ind2)){
			if(isPromptMuon(ind1)){
				if( isTightMuon(ind2)){
					nt_mm++;	
					if( isPromptMuon(ind2)) np_mm++;
					if(!isPromptMuon(ind2)) nf_mm++;
				}
				if(!isTightMuon(ind2))  nl_mm++;
			}
			else if(isPromptMuon(ind2)){
				if( isTightMuon(ind1)){
					nt_mm++;
					if( isPromptMuon(ind1)) np_mm++;
					if(!isPromptMuon(ind1)) nf_mm++;
				}
				if(!isTightMuon(ind1))  nl_mm++;				
			}			
		}
		
		
		// EE
		fCurrentChannel = Elec;
		if(isSSLLElEvent(ind1, ind2)){
			if(isPromptElectron(ind1)){
				if( isTightElectron(ind2)){
					nt_ee++;
					if( isPromptElectron(ind2)) np_ee++;
					if(!isPromptElectron(ind2)) nf_ee++;
				}
				if(!isTightElectron(ind2))  nl_ee++;
			}
			else if(isPromptElectron(ind2)){
				if( isTightElectron(ind1)){
					nt_ee++;
					if( isPromptElectron(ind1)) np_ee++;
					if(!isPromptElectron(ind1)) nf_ee++;
				}
				if(!isTightElectron(ind1))  nl_ee++;				
			}			
		}
		
		// EMu
		fCurrentChannel = ElMu;
		if(isSSLLElMuEvent(ind1, ind2)){
			if(isPromptElectron(ind2)){
				if( isTightMuon(ind1)){
					nt_em++;
					if( isPromptMuon(ind1)) np_em++;
					if(!isPromptMuon(ind1)) nf_em++;
				}
				if(!isTightMuon(ind1))  nl_em++;				
			}			
			else if(isPromptMuon(ind1)){
				if( isTightElectron(ind2)){
					nt_me++;
					if( isPromptElectron(ind2)) np_me++;
					if(!isPromptElectron(ind2)) nf_me++;
				}
				if(!isTightElectron(ind2))  nl_me++;
			}
		}
	}
	cout << endl;
	S->cleanUp();
	// Scale by luminosity:
	float scale = fLumiNorm / fSamples[TTJets]->lumi;
	nt_mm*=scale; nl_mm*=scale; np_mm*=scale; nf_mm*=scale;
	nt_em*=scale; nl_em*=scale; np_em*=scale; nf_em*=scale;
	nt_me*=scale; nl_me*=scale; np_me*=scale; nf_me*=scale;
	nt_ee*=scale; nl_ee*=scale; np_ee*=scale; nf_ee*=scale;	
	
	///////////////////////////////////////////////////////////////////////////////////
	// PREDICTIONS ////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////
	FPRatios *fMMFPRatios = new FPRatios();
	FPRatios *fMEFPRatios = new FPRatios();
	FPRatios *fEMFPRatios = new FPRatios();
	FPRatios *fEEFPRatios = new FPRatios();

	fMMFPRatios->SetVerbose(fVerbose);
	fMEFPRatios->SetVerbose(fVerbose);
	fEMFPRatios->SetVerbose(fVerbose);
	fEEFPRatios->SetVerbose(fVerbose);

	fMMFPRatios->SetMuFratios(mufratio_allmc, mufratio_allmc_e);
	fMMFPRatios->SetMuPratios(mupratio_allmc, mupratio_allmc_e);

	fEEFPRatios->SetElFratios(elfratio_allmc, elfratio_allmc_e);
	fEEFPRatios->SetElPratios(elpratio_allmc, elpratio_allmc_e);

	fEMFPRatios->SetMuFratios(mufratio_allmc, mufratio_allmc_e);
	fEMFPRatios->SetMuPratios(mupratio_allmc, mupratio_allmc_e);
	fMEFPRatios->SetElFratios(elfratio_allmc, elfratio_allmc_e);
	fMEFPRatios->SetElPratios(elpratio_allmc, elpratio_allmc_e);

	vector<double> vpt, veta;
	vpt.push_back(30.); // Fake pts and etas (first electron then muon)
	veta.push_back(0.);


	// MuMu Channel
	vector<double> v_nt_mm;
	v_nt_mm.push_back(nl_mm);
	v_nt_mm.push_back(nt_mm);

	fMMFPRatios->NevtTopol(0, 1, v_nt_mm);

	vector<double> MM_Nev   = fMMFPRatios->NevtPass(vpt, veta);
	vector<double> MM_Estat = fMMFPRatios->NevtPassErrStat();
	vector<double> MM_Esyst = fMMFPRatios->NevtPassErrSyst();

	// EM Channel
	vector<double> v_nt_em;
	v_nt_em.push_back(nl_em);
	v_nt_em.push_back(nt_em);
	fEMFPRatios->NevtTopol(0, 1, v_nt_em);
	vector<double> EM_Nev   = fEMFPRatios->NevtPass(vpt, veta);
	vector<double> EM_Estat = fEMFPRatios->NevtPassErrStat();
	vector<double> EM_Esyst = fEMFPRatios->NevtPassErrSyst();

	// ME Channel
	vector<double> v_nt_me;
	v_nt_me.push_back(nl_me);
	v_nt_me.push_back(nt_me);
	fMEFPRatios->NevtTopol(1, 0, v_nt_me);
	vector<double> ME_Nev   = fMEFPRatios->NevtPass(vpt, veta);
	vector<double> ME_Estat = fMEFPRatios->NevtPassErrStat();
	vector<double> ME_Esyst = fMEFPRatios->NevtPassErrSyst();


	// EE Channel
	vector<double> v_nt_ee;
	v_nt_ee.push_back(nl_ee);
	v_nt_ee.push_back(nt_ee);
	fEEFPRatios->NevtTopol(1, 0, v_nt_ee);
	vector<double> EE_Nev   = fEEFPRatios->NevtPass(vpt, veta);
	vector<double> EE_Estat = fEEFPRatios->NevtPassErrStat();
	vector<double> EE_Esyst = fEEFPRatios->NevtPassErrSyst();
	
	///////////////////////////////////////////////////////////////////////////////////
	// PRINTOUT ///////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////
	OUT << "-----------------------------------------------------------------------------------------------------" << endl;
	OUT << "     RATIOS  ||     Mu-fRatio      |     Mu-pRatio      ||     El-fRatio      |     El-pRatio      ||" << endl;
	OUT << "-----------------------------------------------------------------------------------------------------" << endl;
	OUT << setw(12) << "    allMC    ||";
	OUT << setw(7)  << setprecision(2) << mufratio_allmc << " +/- " << setw(7) << setprecision(2) << mufratio_allmc_e << " |";
	OUT << setw(7)  << setprecision(2) << mupratio_allmc << " +/- " << setw(7) << setprecision(2) << mupratio_allmc_e << " ||";
	OUT << setw(7)  << setprecision(2) << elfratio_allmc << " +/- " << setw(7) << setprecision(2) << elfratio_allmc_e << " |";
	OUT << setw(7)  << setprecision(2) << elpratio_allmc << " +/- " << setw(7) << setprecision(2) << elpratio_allmc_e << " ||";
	OUT << endl;
	OUT << "-----------------------------------------------------------------------------------------------------" << endl << endl;
	OUT << "---------------------------------------------------------------------------------------------------" << endl;
	OUT << "             ||       Mu/Mu       ||        E/Mu       ||       Mu/E        ||        E/E        ||" << endl;
	OUT << "             ||    Nt   |    Nl   ||    Nt   |   Nl    ||    Nt   |   Nl    ||    Nt   |   Nl    ||" << endl;
	OUT << "---------------------------------------------------------------------------------------------------" << endl;
	OUT << setw(12) << "Observed" << " || ";
	OUT << setw(7)  << numbForm(nt_mm) << " | ";
	OUT << setw(7)  << numbForm(nl_mm) << " || ";
	OUT << setw(7)  << numbForm(nt_em) << " | ";
	OUT << setw(7)  << numbForm(nl_em) << " || ";
	OUT << setw(7)  << numbForm(nt_me) << " | ";
	OUT << setw(7)  << numbForm(nl_me) << " || ";
	OUT << setw(7)  << numbForm(nt_ee) << " | ";
	OUT << setw(7)  << numbForm(nl_ee) << " || ";
	OUT << endl;
	OUT << "---------------------------------------------------------------------------------------------------" << endl;
	OUT << setw(12) << "Truth ratio" << " || ";
	OUT << setw(17)  << numbForm(nt_mm/nl_mm) << " || ";
	OUT << setw(17)  << numbForm(nt_em/nl_em) << " || ";
	OUT << setw(17)  << numbForm(nt_me/nl_me) << " || ";
	OUT << setw(17)  << numbForm(nt_ee/nl_ee) << " || ";
	OUT << endl;
	OUT << "---------------------------------------------------------------------------------------------------" << endl;
	OUT << endl;
	OUT << "---------------------------------------------------------------------------------------------------" << endl;
	OUT << "             ||    Np   |    Nf   ||    Np   |    Nf   ||    Np   |    Nf   ||    Np   |    Nf   ||" << endl;
	OUT << "---------------------------------------------------------------------------------------------------" << endl;
	OUT << setw(12) << "MC Truth" << " || ";
	OUT << setw(7)  << numbForm(np_mm) << " | ";
	OUT << setw(7)  << numbForm(nf_mm) << " || ";
	OUT << setw(7)  << numbForm(np_em) << " | ";
	OUT << setw(7)  << numbForm(nf_em) << " || ";
	OUT << setw(7)  << numbForm(np_me) << " | ";
	OUT << setw(7)  << numbForm(nf_me) << " || ";
	OUT << setw(7)  << numbForm(np_ee) << " | ";
	OUT << setw(7)  << numbForm(nf_ee) << " || ";
	OUT << endl;
	OUT << setw(12) << "Predicted" << " || ";
	OUT << setw(7)  << numbForm(MM_Nev[0]) << " | ";
	OUT << setw(7)  << numbForm(MM_Nev[1]) << " || ";
	OUT << setw(7)  << numbForm(EM_Nev[0]) << " | ";
	OUT << setw(7)  << numbForm(EM_Nev[1]) << " || ";
	OUT << setw(7)  << numbForm(ME_Nev[0]) << " | ";
	OUT << setw(7)  << numbForm(ME_Nev[1]) << " || ";
	OUT << setw(7)  << numbForm(EE_Nev[0]) << " | ";
	OUT << setw(7)  << numbForm(EE_Nev[1]) << " || ";
	OUT << endl;
	OUT << "---------------------------------------------------------------------------------------------------" << endl;
	
	OUT.close();
	delete fMMFPRatios;
	delete fEMFPRatios;
	delete fMEFPRatios;
	delete fEEFPRatios;
}

void SSDLPlotter::storeWeightedPred(){
	fillRatios(fMuData, fEGData, 0);
	fillRatios(fMCBGMuEnr, fMCBG, 1);
	TFile *pFile = TFile::Open(fOutputFileName);
	TTree *sigtree; getObjectSafe(pFile, "SigEvents", sigtree);

	string *sname = 0;
	int   stype, flav, cat, njets, nbjets;
	float puweight, pT1, pT2, HT, MET, MT2;

	sigtree->SetBranchAddress("SName",    &sname);
	sigtree->SetBranchAddress("SType",    &stype);
	sigtree->SetBranchAddress("PUWeight", &puweight);
	sigtree->SetBranchAddress("Flavor",   &flav);
	sigtree->SetBranchAddress("pT1",      &pT1);
	sigtree->SetBranchAddress("pT2",      &pT2);
	sigtree->SetBranchAddress("TLCat",    &cat);
	sigtree->SetBranchAddress("HT",       &HT);
	sigtree->SetBranchAddress("MET",      &MET);
	sigtree->SetBranchAddress("MT2",      &MT2);
	sigtree->SetBranchAddress("NJ",       &njets);
	sigtree->SetBranchAddress("NbJ",      &nbjets);
	
	FakeRatios *FR = new FakeRatios();

	float npp(0.), npf(0.), nfp(0.), nff(0.);
	float f1(0.), f2(0.), p1(0.), p2(0.);
	
	for( int i = 0; i < sigtree->GetEntries(); i++ ){
		sigtree->GetEntry(i);
		Sample *S = fSampleMap[TString(*sname)];
		int datamc = S->datamc;
		
		gChannel chan = gChannel(flav);
		f1 = getFRatio(chan, pT1, S->datamc);
		f2 = getFRatio(chan, pT2, S->datamc);
		p1 = getPRatio(chan, pT1, S->datamc);
		p2 = getPRatio(chan, pT2, S->datamc);
		if(chan == ElMu){
			f1 = getFRatio(Muon, pT1, S->datamc);
			f2 = getFRatio(Elec, pT2, S->datamc);
			p1 = getPRatio(Muon, pT1, S->datamc);
			p2 = getPRatio(Elec, pT2, S->datamc);			
		}

		// Get the weights (don't depend on event selection)
		npp = FR->getWpp(FakeRatios::gTLCat(cat), f1, f2, p1, p2);
		npf = FR->getWpf(FakeRatios::gTLCat(cat), f1, f2, p1, p2);
		nfp = FR->getWfp(FakeRatios::gTLCat(cat), f1, f2, p1, p2);
		nff = FR->getWff(FakeRatios::gTLCat(cat), f1, f2, p1, p2);			

		// Store them in the right places for the different purposes
		// Integrated predictions
		for(gRegion r = region_begin; r < gNREGIONS; r = gRegion(r+1)){
			// Event Selection:
			if(HT     < Region::minHT    [r] || HT  > Region::maxHT [r]) continue;
			if(MET    < Region::minMet   [r] || MET > Region::maxMet[r]) continue;
			if(njets  < Region::minNjets [r]) continue;
			if(nbjets < Region::minNbjets[r]) continue;

			S->numbers[r][chan].npp += puweight * npp;
			S->numbers[r][chan].npf += puweight * npf;
			S->numbers[r][chan].nfp += puweight * nfp;
			S->numbers[r][chan].nff += puweight * nff;
		}
		
		// Differential predictions
		if(MET > 50.){
			fillWithoutOF(S->diffyields[chan].hnpp[0], HT, puweight * npp);
			fillWithoutOF(S->diffyields[chan].hnpf[0], HT, puweight * npf);
			fillWithoutOF(S->diffyields[chan].hnfp[0], HT, puweight * nfp);
			fillWithoutOF(S->diffyields[chan].hnff[0], HT, puweight * nff);

			fillWithoutOF(S->diffyields[chan].hnpp[4], njets + 0.5, puweight * npp);
			fillWithoutOF(S->diffyields[chan].hnpf[4], njets + 0.5, puweight * npf);
			fillWithoutOF(S->diffyields[chan].hnfp[4], njets + 0.5, puweight * nfp);
			fillWithoutOF(S->diffyields[chan].hnff[4], njets + 0.5, puweight * nff);			
		}
		if(MET > 120.){
			fillWithoutOF(S->diffyields[chan].hnpp[1], HT, puweight * npp);
			fillWithoutOF(S->diffyields[chan].hnpf[1], HT, puweight * npf);
			fillWithoutOF(S->diffyields[chan].hnfp[1], HT, puweight * nfp);
			fillWithoutOF(S->diffyields[chan].hnff[1], HT, puweight * nff);			
		}
		if(HT > 80. && MET > 30. && njets > 1){
			fillWithoutOF(S->diffyields[chan].hnpp[5], MT2, puweight * npp);
			fillWithoutOF(S->diffyields[chan].hnpf[5], MT2, puweight * npf);
			fillWithoutOF(S->diffyields[chan].hnfp[5], MT2, puweight * nfp);
			fillWithoutOF(S->diffyields[chan].hnff[5], MT2, puweight * nff);

			fillWithoutOF(S->diffyields[chan].hnpp[6], pT1, puweight * npp);
			fillWithoutOF(S->diffyields[chan].hnpf[6], pT1, puweight * npf);
			fillWithoutOF(S->diffyields[chan].hnfp[6], pT1, puweight * nfp);
			fillWithoutOF(S->diffyields[chan].hnff[6], pT1, puweight * nff);

			fillWithoutOF(S->diffyields[chan].hnpp[7], pT2, puweight * npp);
			fillWithoutOF(S->diffyields[chan].hnpf[7], pT2, puweight * npf);
			fillWithoutOF(S->diffyields[chan].hnfp[7], pT2, puweight * nfp);
			fillWithoutOF(S->diffyields[chan].hnff[7], pT2, puweight * nff);
		}
		if(HT > 200. && MET > 30. && njets > 1){
			fillWithoutOF(S->diffyields[chan].hnpp[2], MET, puweight * npp);
			fillWithoutOF(S->diffyields[chan].hnpf[2], MET, puweight * npf);
			fillWithoutOF(S->diffyields[chan].hnfp[2], MET, puweight * nfp);
			fillWithoutOF(S->diffyields[chan].hnff[2], MET, puweight * nff);

			fillWithoutOF(S->diffyields[chan].hnpp[8], nbjets + 0.5, puweight * npp);
			fillWithoutOF(S->diffyields[chan].hnpf[8], nbjets + 0.5, puweight * npf);
			fillWithoutOF(S->diffyields[chan].hnfp[8], nbjets + 0.5, puweight * nfp);
			fillWithoutOF(S->diffyields[chan].hnff[8], nbjets + 0.5, puweight * nff);
		}
		if(HT > 450. && MET > 30. && njets > 1){
			fillWithoutOF(S->diffyields[chan].hnpp[3], MET, puweight * npp);
			fillWithoutOF(S->diffyields[chan].hnpf[3], MET, puweight * npf);
			fillWithoutOF(S->diffyields[chan].hnfp[3], MET, puweight * nfp);
			fillWithoutOF(S->diffyields[chan].hnff[3], MET, puweight * nff);
		}
	}
	delete FR;
}
float SSDLPlotter::getFRatio(gChannel chan, float pt, int datamc){
	// if(chan == Muon) return 0.0672; // flat ratios
	// if(chan == Elec) return 0.224;

	const float mu_flatout = 35.;
	const float el_flatout = 40.;
	
	if(chan == Muon){
		TH1D *histo          = fH1D_MufRatio;
		if(datamc > 0) histo = fH1D_MufRatio_MC;
		if(!histo){
			cerr << "SSDLPlotter::getPRatio ==> Warning: ratio histo not filled, exiting" << endl;
			exit(-1);
		}
		if(pt >= mu_flatout){
			Int_t asym_val_bin = histo->FindBin(mu_flatout);
			return histo->GetBinContent(asym_val_bin - 1); // findbin(35) returns bin from 35-45, want the previous one
		}
		Int_t binnumb = histo->FindBin(pt);
		return histo->GetBinContent(binnumb);
	}
	if(chan == Elec){
		TH1D *histo          = fH1D_ElfRatio;
		if(datamc > 0) histo = fH1D_ElfRatio_MC;
		if(!histo){
			cerr << "SSDLPlotter::getPRatio ==> Warning: ratio histo not filled, exiting" << endl;
			exit(-1);
		}
		if(pt >= el_flatout){
			Int_t asym_val_bin = histo->FindBin(el_flatout);
			return histo->GetBinContent(asym_val_bin - 1);
		}
		Int_t binnumb = histo->FindBin(pt);
		return histo->GetBinContent(binnumb);
	}
}
float SSDLPlotter::getPRatio(gChannel chan, float pt, int datamc){
	// if(chan == Muon) return 0.941; // flat ratios
	// if(chan == Elec) return 0.922;
	
	if(chan == Muon){
		TH1D *histo          = fH1D_MupRatio;
		if(datamc > 0) histo = fH1D_MupRatio_MC;
		if(!histo){
			cerr << "SSDLPlotter::getPRatio ==> Warning: ratio histo not filled, exiting" << endl;
			exit(-1);
		}
		Int_t binnumb = histo->FindBin(pt);
		return histo->GetBinContent(binnumb);
	}
	if(chan == Elec){
		TH1D *histo          = fH1D_ElpRatio;
		if(datamc > 0) histo = fH1D_ElpRatio_MC;
		if(!histo){
			cerr << "SSDLPlotter::getPRatio ==> Warning: ratio histo not filled, exiting" << endl;
			exit(-1);
		}
		Int_t binnumb = histo->FindBin(pt);
		return histo->GetBinContent(binnumb);
	}
}

//____________________________________________________________________________
void SSDLPlotter::printYields(gChannel chan, float lumiscale){
	cout << setfill('-') << setw(97) << "-" << endl;
	TString name = "Mu/Mu";
	if(chan == Elec) name = "E/E";
	if(chan == ElMu) name = "E/Mu";
	cout << " Printing yields for " << name << " channel..." << endl;
	if(lumiscale > -1.0) cout << " Numbers scaled to " << lumiscale << " /pb" << endl;
	cout << "           Name |   Nt2   |   Nt10  |   Nt01  |   Nt0   |   Nsst  |   Nssl  |   NZ t  |   NZ l  |" << endl;
	cout << setfill('-') << setw(97) << "-" << endl;
	cout << setfill(' ');
	for(gSample i = sample_begin; i < gNSAMPLES; i=gSample(i+1)){
		Sample *S = fSamples[i];
		NumberSet numbers = S->numbers[Baseline][chan];
		float scale = lumiscale / S->lumi;
		if(S->datamc == 0 || scale < 0) scale = 1;
		cout << setw(15) << S->sname << " |";
		cout << setw(8)  << setprecision(3) << scale*numbers.nt2  << " |";
		cout << setw(8)  << setprecision(3) << scale*numbers.nt10 << " |";
		cout << setw(8)  << setprecision(3) << scale*numbers.nt01 << " |";
		cout << setw(8)  << setprecision(3) << scale*numbers.nt0  << " |";
		cout << setw(8)  << setprecision(3) << scale*numbers.nsst << " |"; 
		cout << setw(8)  << setprecision(3) << scale*numbers.nssl << " |";
		cout << setw(8)  << setprecision(3) << scale*numbers.nzt  << " |";
		cout << setw(8)  << setprecision(3) << scale*numbers.nzl  << " |";
		cout << endl;
	}

	cout << setfill('-') << setw(97) << "-" << endl;
}
void SSDLPlotter::printYieldsShort(float luminorm){
	vector<int> musamples;
	vector<int> elsamples;
	vector<int> emusamples;

	cout << "---------------------------------------------------" << endl;
	cout << "Printing yields" << endl;
	musamples = fMuData;
	elsamples = fEGData;
	emusamples = fMuEGData;

	float nt20[gNCHANNELS],    nt10[gNCHANNELS],    nt01[gNCHANNELS],    nt00[gNCHANNELS];

	// float nt2_mumu(0.),    nt10_mumu(0.),    nt0_mumu(0.);
	// float nt2_emu(0.),    nt10_emu(0.),    nt01_emu(0.),    nt0_emu(0.);
	// float nt2_ee(0.),    nt10_ee(0.),    nt0_ee(0.);

	for(size_t i = 0; i < musamples.size(); ++i){
		Sample *S = fSamples[musamples[i]];
		nt20[Muon]    += S->numbers[Baseline][Muon].nt2;
		nt10[Muon]    += S->numbers[Baseline][Muon].nt10;
		nt00[Muon]    += S->numbers[Baseline][Muon].nt0;
	}

	for(size_t i = 0; i < emusamples.size(); ++i){
		Sample *S = fSamples[emusamples[i]];
		nt20[ElMu]    += S->numbers[Baseline][ElMu].nt2;
		nt10[ElMu]    += S->numbers[Baseline][ElMu].nt10;
		nt01[ElMu]    += S->numbers[Baseline][ElMu].nt01;
		nt00[ElMu]    += S->numbers[Baseline][ElMu].nt0;
	}		

	for(size_t i = 0; i < elsamples.size(); ++i){
		Sample *S = fSamples[elsamples[i]];
		nt20[Elec]    += S->numbers[Baseline][Elec].nt2;
		nt10[Elec]    += S->numbers[Baseline][Elec].nt10;
		nt00[Elec]    += S->numbers[Baseline][Elec].nt0;
	}		

	// for(size_t i = 0; i < musamples.size(); ++i){
	// 	int index = musamples[i];
	// 	Channel *mumu = fSamples[index]->mm;
	// 	nt2_mumu  += mumu->numbers.nt2;
	// 	nt10_mumu += mumu->numbers.nt10;
	// 	nt0_mumu  += mumu->numbers.nt0;
	// 	nt2_mumu_e2  += mumu->numbers.nt2;
	// 	nt10_mumu_e2 += mumu->numbers.nt10;
	// 	nt0_mumu_e2  += mumu->numbers.nt0;
	// }		
	// for(size_t i = 0; i < emusamples.size(); ++i){
	// 	int index = emusamples[i];
	// 	Channel *emu = fSamples[index]->em;
	// 	nt2_emu  += emu->numbers.nt2;
	// 	nt10_emu += emu->numbers.nt10;
	// 	nt01_emu += emu->numbers.nt01;
	// 	nt0_emu  += emu->numbers.nt0;
	// 	nt2_emu_e2  += emu->numbers.nt2;
	// 	nt10_emu_e2 += emu->numbers.nt10;
	// 	nt01_emu_e2 += emu->numbers.nt01;
	// 	nt0_emu_e2  += emu->numbers.nt0;
	// }		
	// for(size_t i = 0; i < elsamples.size(); ++i){
	// 	int index = elsamples[i];
	// 	Channel *ee = fSamples[index]->ee;
	// 	nt2_ee  += ee->numbers.nt2;
	// 	nt10_ee += ee->numbers.nt10;
	// 	nt0_ee  += ee->numbers.nt0;
	// 	nt2_ee_e2  += ee->numbers.nt2;
	// 	nt10_ee_e2 += ee->numbers.nt10;
	// 	nt0_ee_e2  += ee->numbers.nt0;
	// }
	cout << "-------------------------------------------------------------------------------------------------------------------" << endl;
	cout << "          ||            Mu/Mu            ||                   E/Mu                ||             E/E             ||" << endl;
	cout << "  YIELDS  ||   Nt2   |   Nt1   |   Nt0   ||   Nt2   |   Nt10  |   Nt01  |   Nt0   ||   Nt2   |   Nt1   |   Nt0   ||" << endl;
	cout << "-------------------------------------------------------------------------------------------------------------------" << endl;

	float nt2sum_mumu(0.), nt10sum_mumu(0.), nt0sum_mumu(0.);
	float nt2sum_emu(0.), nt10sum_emu(0.), nt01sum_emu(0.), nt0sum_emu(0.);
	float nt2sum_ee(0.), nt10sum_ee(0.), nt0sum_ee(0.);
	for(size_t i = 0; i < fMCBG.size(); ++i){
		int index = fMCBG[i];
		Sample *S = fSamples[index];
		float scale = luminorm / S->lumi;
		if(luminorm < 0) scale = 1;
		nt2sum_mumu  += scale*S->numbers[Baseline][Muon]    .nt2;
		nt10sum_mumu += scale*S->numbers[Baseline][Muon]    .nt10;
		nt0sum_mumu  += scale*S->numbers[Baseline][Muon]    .nt0;
		nt2sum_emu   += scale*S->numbers[Baseline][ElMu]     .nt2;
		nt10sum_emu  += scale*S->numbers[Baseline][ElMu]     .nt10;
		nt01sum_emu  += scale*S->numbers[Baseline][ElMu]     .nt01;
		nt0sum_emu   += scale*S->numbers[Baseline][ElMu]     .nt0;
		nt2sum_ee    += scale*S->numbers[Baseline][Elec].nt2;
		nt10sum_ee   += scale*S->numbers[Baseline][Elec].nt10;
		nt0sum_ee    += scale*S->numbers[Baseline][Elec].nt0;

		cout << setw(9) << S->sname << " || ";
		cout << setw(7) << scale*S->numbers[Baseline][Muon]    .nt2  << " | ";
		cout << setw(7) << scale*S->numbers[Baseline][Muon]    .nt10 << " | ";
		cout << setw(7) << scale*S->numbers[Baseline][Muon]    .nt0  << " || ";
		cout << setw(7) << scale*S->numbers[Baseline][ElMu]     .nt2  << " | ";
		cout << setw(7) << scale*S->numbers[Baseline][ElMu]     .nt10 << " | ";
		cout << setw(7) << scale*S->numbers[Baseline][ElMu]     .nt01 << " | ";
		cout << setw(7) << scale*S->numbers[Baseline][ElMu]     .nt0  << " || ";
		cout << setw(7) << scale*S->numbers[Baseline][Elec].nt2  << " | ";
		cout << setw(7) << scale*S->numbers[Baseline][Elec].nt10 << " | ";
		cout << setw(7) << scale*S->numbers[Baseline][Elec].nt0  << " || ";
		cout << endl;
	}
	if(luminorm > 0){
		cout << "-------------------------------------------------------------------------------------------------------------------" << endl;
		cout << setw(9) << "MC sum" << " || ";
		cout << setw(7) << nt2sum_mumu  << " | ";
		cout << setw(7) << nt10sum_mumu << " | ";
		cout << setw(7) << nt0sum_mumu  << " || ";
		cout << setw(7) << nt2sum_emu   << " | ";
		cout << setw(7) << nt10sum_emu  << " | ";
		cout << setw(7) << nt01sum_emu  << " | ";
		cout << setw(7) << nt0sum_emu   << " || ";
		cout << setw(7) << nt2sum_ee    << " | ";
		cout << setw(7) << nt10sum_ee   << " | ";
		cout << setw(7) << nt0sum_ee    << " || ";
		cout << endl;
	}
	cout << "-------------------------------------------------------------------------------------------------------------------" << endl;
	cout << setw(9) << fSamples[LM0]->sname << " || ";
	float scale = luminorm / fSamples[LM0]->lumi;
	if(luminorm < 0) scale = 1;
	cout << setw(7) << scale * fSamples[LM0]->numbers[Baseline][Muon]    .nt2  << " | ";
	cout << setw(7) << scale * fSamples[LM0]->numbers[Baseline][Muon]    .nt10 << " | ";
	cout << setw(7) << scale * fSamples[LM0]->numbers[Baseline][Muon]    .nt0  << " || ";
	cout << setw(7) << scale * fSamples[LM0]->numbers[Baseline][ElMu]     .nt2  << " | ";
	cout << setw(7) << scale * fSamples[LM0]->numbers[Baseline][ElMu]     .nt10 << " | ";
	cout << setw(7) << scale * fSamples[LM0]->numbers[Baseline][ElMu]     .nt01 << " | ";
	cout << setw(7) << scale * fSamples[LM0]->numbers[Baseline][ElMu]     .nt0  << " || ";
	cout << setw(7) << scale * fSamples[LM0]->numbers[Baseline][Elec].nt2  << " | ";
	cout << setw(7) << scale * fSamples[LM0]->numbers[Baseline][Elec].nt10 << " | ";
	cout << setw(7) << scale * fSamples[LM0]->numbers[Baseline][Elec].nt0  << " || ";
	cout << endl;
	cout << "-------------------------------------------------------------------------------------------------------------------" << endl;
	cout << setw(9) << "data"  << " || ";
	cout << setw(7) << nt20[Muon]     << " | ";
	cout << setw(7) << nt10[Muon]     << " | ";
	cout << setw(7) << nt00[Muon]     << " || ";
	cout << setw(7) << nt20[ElMu]      << " | ";
	cout << setw(7) << nt10[ElMu]      << " | ";
	cout << setw(7) << nt01[ElMu]      << " | ";
	cout << setw(7) << nt00[ElMu]      << " || ";
	cout << setw(7) << nt20[Elec] << " | ";
	cout << setw(7) << nt10[Elec] << " | ";
	cout << setw(7) << nt00[Elec] << " || ";
	cout << endl;
	cout << "-------------------------------------------------------------------------------------------------------------------" << endl;

}

//////////////////////////////////////////////////////////////////////////////
// Geninfo stuff
//____________________________________________________________________________
void SSDLPlotter::makeOriginPlots(gRegion reg){
	gStyle->SetPaintTextFormat("5.2f");
	// useNiceColorPalette();
	bool hasBjets = (Region::minNbjets[reg] > 0);

	fOutputSubDir = "Origins/";
	// make the histograms first. one for ttjets and one for all mc (without signal). this in each channel
	TH2D    *horigin_tt [gNCHANNELS];
	TH2D    *horigin_mc [gNCHANNELS];
	int nbins(15);

	std::vector<int> mcsamples;
	std::vector<int>::const_iterator sampleInd;

	TLatex *lat = new TLatex();
	lat->SetNDC(kTRUE);
	lat->SetTextColor(kBlack);
	lat->SetTextSize(0.04);
	lat->SetTextAlign(12);

	for (int i=0; i < gNCHANNELS; i++){
		if(i == 0) mcsamples = fMCBGMuEnr;
		if(i == 1) mcsamples = fMCBG;
		if(i == 2) mcsamples = fMCBG;
		horigin_tt [i] = new TH2D("OriginHistoTTJ_" + SSDLDumper::gChanLabel[i], "Origin Histogram for TTJets "    + SSDLDumper::gChanLabel[i], nbins, 0, nbins, nbins, 0, nbins);
		horigin_mc [i] = new TH2D("OriginHistoMC_"  + SSDLDumper::gChanLabel[i], "Origin Histogram for total MC "  + SSDLDumper::gChanLabel[i], nbins, 0, nbins, nbins, 0, nbins);
		horigin_tt [i]->Sumw2();
		horigin_mc [i]->Sumw2();
		for (sampleInd=mcsamples.begin(); sampleInd != mcsamples.end(); sampleInd ++){ // sample loop
			Sample *sample = fSamples[*sampleInd];
			TString s_name = sample->sname;
			float scale = fLumiNorm / fSamples[*sampleInd]->lumi;
			if(i == 0) horigin_mc[i]->Add(sample->region[reg][HighPt].mm.nt11_origin, scale);
			if(i == 1) horigin_mc[i]->Add(sample->region[reg][HighPt].em.nt11_origin, scale);
			if(i == 2) horigin_mc[i]->Add(sample->region[reg][HighPt].ee.nt11_origin, scale);
			if (s_name == "TTJets") {
				if(i == 0) horigin_tt[i]->Add(sample->region[reg][HighPt].mm.nt11_origin, scale);
				if(i == 1) horigin_tt[i]->Add(sample->region[reg][HighPt].em.nt11_origin, scale);
				if(i == 2) horigin_tt[i]->Add(sample->region[reg][HighPt].ee.nt11_origin, scale);
			}
		} // end sample loop
		TAxis* xAxis_tt = horigin_tt[i]->GetXaxis();
		TAxis* yAxis_tt = horigin_tt[i]->GetYaxis();
		TAxis* xAxis_mc = horigin_mc[i]->GetXaxis();
		TAxis* yAxis_mc = horigin_mc[i]->GetYaxis();

		horigin_mc[i]->GetZaxis()->SetRangeUser(0., 0.25);
		horigin_tt[i]->GetZaxis()->SetRangeUser(0., 0.25);

		if (i == 0){
			for (int bin = 1; bin<=nbins; bin++) {
				xAxis_tt->SetBinLabel(bin, SSDLDumper::muBinToLabel(bin));
				yAxis_tt->SetBinLabel(bin, SSDLDumper::muBinToLabel(bin));
				xAxis_mc->SetBinLabel(bin, SSDLDumper::muBinToLabel(bin));
				yAxis_mc->SetBinLabel(bin, SSDLDumper::muBinToLabel(bin));
			}
			horigin_mc[i] = mirrorHisto(horigin_mc[i]);
			horigin_tt[i] = mirrorHisto(horigin_tt[i]);
		}
		if (i == 1){
			for (int bin = 1; bin<=nbins; bin++) {
				xAxis_tt->SetBinLabel(bin, SSDLDumper::muBinToLabel(bin));
				yAxis_tt->SetBinLabel(bin, SSDLDumper::elBinToLabel(bin));
				xAxis_mc->SetBinLabel(bin, SSDLDumper::muBinToLabel(bin));
				yAxis_mc->SetBinLabel(bin, SSDLDumper::elBinToLabel(bin));
			}
		}
		if (i == 2){
			for (int bin = 1; bin<=nbins; bin++) {
				xAxis_tt->SetBinLabel(bin, SSDLDumper::elBinToLabel(bin));
				yAxis_tt->SetBinLabel(bin, SSDLDumper::elBinToLabel(bin));
				xAxis_mc->SetBinLabel(bin, SSDLDumper::elBinToLabel(bin));
				yAxis_mc->SetBinLabel(bin, SSDLDumper::elBinToLabel(bin));
			}
			horigin_mc[i] = mirrorHisto(horigin_mc[i]);
			horigin_tt[i] = mirrorHisto(horigin_tt[i]);
		}
		
		horigin_tt[i]->GetXaxis()->LabelsOption("v");
		horigin_tt[i]->GetYaxis()->LabelsOption("d");
		horigin_mc[i]->GetXaxis()->LabelsOption("v");
		horigin_mc[i]->GetYaxis()->LabelsOption("d");
		
		TCanvas *c_temp = new TCanvas("Origin" + SSDLDumper::gChanLabel[i], "Origin plot in region " + SSDLDumper::Region::sname[reg], 0, 0, 600, 600);
		c_temp->cd();
		gPad->SetRightMargin(0.05);
		gPad->SetLeftMargin(0.3);
		gPad->SetBottomMargin(0.3);
		gPad->SetGrid(1,1);

		float latX = 0.035;
		float latY = 0.18;
		// total MC histos
		float total_mc = horigin_mc[i]->Integral();
		horigin_mc[i]->Scale(100./total_mc);
		horigin_mc[i]->Draw("col text");
		lat->SetTextSize(0.04);
		lat->DrawLatex(0.05, 0.93, "Origins in "+SSDLDumper::gChanLabel[i]+" channel");
		lat->DrawLatex(latX, latY, Form("#splitline{H_{T} > %.0f GeV}{E_{T}^{miss} > %.0f GeV}", Region::minHT[reg], Region::minMet[reg]));
		if (hasBjets) lat->DrawLatex(latX, latY-0.09, Form("N_{b-jets} #geq %1d", Region::minNbjets[reg]));
		lat->SetTextSize(0.03);
		lat->DrawLatex(0.75, 0.93, Form("Exp. # ev.: %5.2f", total_mc));
		Util::PrintPDF(c_temp, "Origin_" + SSDLDumper::gChanLabel[i] + "_" + SSDLDumper::Region::sname[reg], fOutputDir + fOutputSubDir);

		// TTJets only histos
		float total_tt = horigin_tt[i]->Integral();
		horigin_tt[i]->Scale(100./total_tt);
		horigin_tt[i]->Draw("col text");
		lat->SetTextSize(0.04);
		lat->DrawLatex(0.05, 0.93, "Origins in "+SSDLDumper::gChanLabel[i]+" channel (ttbar only)");
		lat->DrawLatex(latX, latY, Form("#splitline{H_{T} > %.0f GeV}{E_{T}^{miss} > %.0f GeV}", Region::minHT[reg], Region::minMet[reg]));
		if (hasBjets) lat->DrawLatex(latX, latY-0.09, Form("N_{b-jets} #geq %1d", Region::minNbjets[reg]));
		lat->SetTextSize(0.03);
		lat->DrawLatex(0.75, 0.93, Form("Exp. # ev.: %5.2f", total_tt));
		Util::PrintPDF(c_temp, "Origin_TTJets_" + SSDLDumper::gChanLabel[i] + "_" + SSDLDumper::Region::sname[reg], fOutputDir + fOutputSubDir);

		delete c_temp;
	} // end channel loop
	fOutputSubDir = "";
}
void SSDLPlotter::printOrigins(gRegion reg){
	TString filename = fOutputDir + "Origins_"+SSDLDumper::Region::sname[reg]+".txt";
	fOUTSTREAM.open(filename.Data(), ios::trunc);
	printMuOriginTable(reg);
	fOUTSTREAM << endl << endl;
	printElOriginTable(reg);
	fOUTSTREAM << endl << endl;
	printEMuOriginTable(reg);
	fOUTSTREAM << endl;
	fOUTSTREAM.close();
}
void SSDLPlotter::printMuOriginTable(gRegion reg){
	fOUTSTREAM << "-------------------------------------------" << endl;
	fOUTSTREAM << " Printing origins for Mu/Mu channel..." << endl;

	printMuOriginHeader("NT20");
	for(gSample i = sample_begin; i < gNSAMPLES; i=gSample(i+1)){
		if(i > QCDMuEnr10) continue;
		print2MuOriginsFromSample(fSamples[i], 2, reg);
	}
	fOUTSTREAM << "-------------------------------------------------------------------------------------------------------------------------" << endl;
	printOriginSummary2L(fMCBGMuEnr, 2, Muon, reg);
	fOUTSTREAM << "=========================================================================================================================" << endl << endl;
	printMuOriginHeader("NT10");
	for(gSample i = sample_begin; i < gNSAMPLES; i=gSample(i+1)){
		if(i > QCDMuEnr10) continue;
		print2MuOriginsFromSample(fSamples[i], 1, reg);
	}
	fOUTSTREAM << "-------------------------------------------------------------------------------------------------------------------------" << endl;
	printOriginSummary2L(fMCBGMuEnr, 1, Muon, reg);
	fOUTSTREAM << "=========================================================================================================================" << endl << endl;
	printMuOriginHeader("NT00");
	for(gSample i = sample_begin; i < gNSAMPLES; i=gSample(i+1)){
		if(i > QCDMuEnr10) continue;
		print2MuOriginsFromSample(fSamples[i], 0, reg);
	}
	fOUTSTREAM << "-------------------------------------------------------------------------------------------------------------------------" << endl;
	printOriginSummary2L(fMCBGMuEnr, 0, Muon, reg);
	fOUTSTREAM << "=========================================================================================================================" << endl << endl;

	printMuOriginHeader("SSTi");
	for(gSample i = sample_begin; i < gNSAMPLES; i=gSample(i+1)){
		if(i > QCDMuEnr10) continue;
		printMuOriginFromSample(fSamples[i], 1, reg);
	}
	fOUTSTREAM << "-------------------------------------------------------------------------------------------------------------------------" << endl;
	printOriginSummary(fMCBGMuEnr, 1, Muon, reg);
	fOUTSTREAM << "=========================================================================================================================" << endl << endl;
	printMuOriginHeader("SSLo");
	for(gSample i = sample_begin; i < gNSAMPLES; i=gSample(i+1)){
		if(i > QCDMuEnr10) continue;
		printMuOriginFromSample(fSamples[i], 2, reg);
	}
	fOUTSTREAM << "-------------------------------------------------------------------------------------------------------------------------" << endl;
	printOriginSummary(fMCBGMuEnr, 2, Muon, reg);
	fOUTSTREAM << "=========================================================================================================================" << endl << endl;
	printMuOriginHeader("Z Ti");
	for(gSample i = sample_begin; i < gNSAMPLES; i=gSample(i+1)){
		if(i > QCDMuEnr10) continue;
		printMuOriginFromSample(fSamples[i], 3, reg);
	}
	fOUTSTREAM << "-------------------------------------------------------------------------------------------------------------------------" << endl;
	printOriginSummary(fMCBGMuEnr, 3, Muon, reg);
	fOUTSTREAM << "=========================================================================================================================" << endl << endl;
	printMuOriginHeader("Z Lo");
	for(gSample i = sample_begin; i < gNSAMPLES; i=gSample(i+1)){
		if(i > QCDMuEnr10) continue;
		printMuOriginFromSample(fSamples[i], 4, reg);
	}
	fOUTSTREAM << "-------------------------------------------------------------------------------------------------------------------------" << endl;
	printOriginSummary(fMCBGMuEnr, 4, Muon, reg);
	fOUTSTREAM << "=========================================================================================================================" << endl << endl;
}
void SSDLPlotter::printElOriginTable(gRegion reg){
	fOUTSTREAM << "-------------------------------------------" << endl;
	fOUTSTREAM << " Printing origins for E/E channel..." << endl;

	printElOriginHeader("NT20");
	for(gSample i = sample_begin; i < gNSAMPLES; i=gSample(i+1)){
		print2ElOriginsFromSample(fSamples[i], 2, reg);
	}
	fOUTSTREAM << "-------------------------------------------------------------------------------------------------------------------------" << endl;
	printOriginSummary2L(fMCBG, 2, Elec, reg);
	fOUTSTREAM << "=============================================================================================================================================" << endl << endl;
	printElOriginHeader("NT10");
	for(gSample i = sample_begin; i < gNSAMPLES; i=gSample(i+1)){
		print2ElOriginsFromSample(fSamples[i], 1, reg);
	}
	fOUTSTREAM << "-------------------------------------------------------------------------------------------------------------------------" << endl;
	printOriginSummary2L(fMCBG, 1, Elec, reg);
	fOUTSTREAM << "=============================================================================================================================================" << endl << endl;
	printElOriginHeader("NT00");
	for(gSample i = sample_begin; i < gNSAMPLES; i=gSample(i+1)){
		print2ElOriginsFromSample(fSamples[i], 0, reg);
	}
	fOUTSTREAM << "-------------------------------------------------------------------------------------------------------------------------" << endl;
	printOriginSummary2L(fMCBG, 0, Elec, reg);
	fOUTSTREAM << "=============================================================================================================================================" << endl << endl;

	printElOriginHeader("SSTi");
	for(gSample i = sample_begin; i < gNSAMPLES; i=gSample(i+1)){
		printElOriginFromSample(fSamples[i], 1, reg);
	}
	fOUTSTREAM << "-------------------------------------------------------------------------------------------------------------------------" << endl;
	printOriginSummary(fMCBG, 1, Elec, reg);
	fOUTSTREAM << "=============================================================================================================================================" << endl << endl;
	printElOriginHeader("SSLo");
	for(gSample i = sample_begin; i < gNSAMPLES; i=gSample(i+1)){
		printElOriginFromSample(fSamples[i], 2, reg);
	}
	fOUTSTREAM << "-------------------------------------------------------------------------------------------------------------------------" << endl;
	printOriginSummary(fMCBG, 2, Elec, reg);
	fOUTSTREAM << "=============================================================================================================================================" << endl << endl;
	printElOriginHeader("Z Ti");
	for(gSample i = sample_begin; i < gNSAMPLES; i=gSample(i+1)){
		printElOriginFromSample(fSamples[i], 3, reg);
	}
	fOUTSTREAM << "-------------------------------------------------------------------------------------------------------------------------" << endl;
	printOriginSummary(fMCBG, 3, Elec, reg);
	fOUTSTREAM << "=============================================================================================================================================" << endl << endl;
	printElOriginHeader("Z Lo");
	for(gSample i = sample_begin; i < gNSAMPLES; i=gSample(i+1)){
		printElOriginFromSample(fSamples[i], 4, reg);
	}
	fOUTSTREAM << "-------------------------------------------------------------------------------------------------------------------------" << endl;
	printOriginSummary(fMCBG, 4, Elec, reg);
	fOUTSTREAM << "=============================================================================================================================================" << endl;
}
void SSDLPlotter::printEMuOriginTable(gRegion reg){
	fOUTSTREAM << "-------------------------------------------" << endl;
	fOUTSTREAM << " Printing origins for E/Mu channel..." << endl;

	printEMuOriginHeader("NT20");
	for(gSample i = sample_begin; i < gNSAMPLES; i=gSample(i+1)){
		printEMuOriginsFromSample(fSamples[i], 2, reg);
	}
	fOUTSTREAM << "---------------------------------------------------------------------------------------------------------------------------------------------" << endl;
	printOriginSummary2L(fMCBG, 2, ElMu, reg);
	fOUTSTREAM << "=============================================================================================================================================" << endl << endl;
	printEMuOriginHeader("NT10");
	for(gSample i = sample_begin; i < gNSAMPLES; i=gSample(i+1)){
		printEMuOriginsFromSample(fSamples[i], 1, reg);
	}
	fOUTSTREAM << "---------------------------------------------------------------------------------------------------------------------------------------------" << endl;
	printOriginSummary2L(fMCBG, 1, ElMu, reg);
	fOUTSTREAM << "=============================================================================================================================================" << endl << endl;
	printEMuOriginHeader("NT01");
	for(gSample i = sample_begin; i < gNSAMPLES; i=gSample(i+1)){
		printEMuOriginsFromSample(fSamples[i], 10, reg);
	}
	fOUTSTREAM << "---------------------------------------------------------------------------------------------------------------------------------------------" << endl;
	printOriginSummary2L(fMCBG, 10, ElMu, reg);
	fOUTSTREAM << "=============================================================================================================================================" << endl << endl;
	printEMuOriginHeader("NT00");
	for(gSample i = sample_begin; i < gNSAMPLES; i=gSample(i+1)){
		printEMuOriginsFromSample(fSamples[i], 0, reg);
	}
	fOUTSTREAM << "---------------------------------------------------------------------------------------------------------------------------------------------" << endl;
	printOriginSummary2L(fMCBG, 0, ElMu, reg);
	fOUTSTREAM << "=============================================================================================================================================" << endl << endl;
}

void SSDLPlotter::printMuOriginHeader(TString name){
	fOUTSTREAM << "-------------------------------------------------------------------------------------------------------------------------" << endl;
	fOUTSTREAM << "    " << name;
	fOUTSTREAM << "            | Fakes   | W       | Z       | tau     | ud had. | s hadr. | c hadr. | b hadr. | strings | unid    |" << endl;
	fOUTSTREAM << "=========================================================================================================================" << endl;	
}
void SSDLPlotter::printElOriginHeader(TString name){
	fOUTSTREAM << "---------------------------------------------------------------------------------------------------------------------------------------------" << endl;
	fOUTSTREAM << "    " << name;
	fOUTSTREAM << "            |           Fakes             |                                  real electrons from...                                 |" << endl;
	fOUTSTREAM << "                    | Mismatc | Gam/Con | Hadron. | W       | Z       | tau     | ud had. | s hadr. | c hadr. | b hadr. | strings | unid    |" << endl;
	fOUTSTREAM << "=============================================================================================================================================" << endl;
}
void SSDLPlotter::printEMuOriginHeader(TString name){
	TString descr = "";
	if(name == "NT20") descr = "both tight";
	if(name == "NT10") descr = "muon tight";
	if(name == "NT01") descr = "elec tight";
	if(name == "NT00") descr = "none tight";
	fOUTSTREAM << "---------------------------------------------------------------------------------------------------------------------------------------------" << endl;
	fOUTSTREAM << "    " << name;
	fOUTSTREAM << "            |           Fakes             |                             real electrons/muons from...                                |" << endl;
	fOUTSTREAM << "      " << descr;
	fOUTSTREAM << "    | Mismatc | Gam/Con | Hadron. | W       | Z       | tau     | ud had. | s hadr. | c hadr. | b hadr. | strings | unid    |" << endl;
	fOUTSTREAM << "=============================================================================================================================================" << endl;
}

void SSDLPlotter::printMuOriginFromSample(Sample *S, int toggle, gRegion reg, gHiLoSwitch hilo){
	if(S->datamc == 0 || S->datamc == 2 ) return;
	if(toggle != 1 && toggle != 2 && toggle != 3 && toggle != 4) return;
	Channel *C = &S->region[reg][hilo].mm;

	TH1D *histo;
	if(toggle == 1) histo = (TH1D*)C->sst_origin->Clone();
	if(toggle == 2) histo = (TH1D*)C->ssl_origin->Clone();
	if(toggle == 3) histo = (TH1D*)C->zt_origin->Clone();
	if(toggle == 4) histo = (TH1D*)C->zl_origin->Clone();

	if(histo->Integral() == 0.) return;
	histo->Scale(100./histo->Integral());

	fOUTSTREAM << setw(15) << S->sname;
	fOUTSTREAM << "     |";
	for(size_t i = 1; i <= 9; ++i)	fOUTSTREAM << setw(6)  << setprecision(3) << Form(" %6.2f%%", histo->GetBinContent(i)) << " |";
	fOUTSTREAM << setw(6)  << setprecision(3) << Form(" %6.2f%%", histo->GetBinContent(15)) << " |" << endl;
}
void SSDLPlotter::printElOriginFromSample(Sample *S, int toggle, gRegion reg, gHiLoSwitch hilo){
	if(S->datamc == 0 || S->datamc == 2 ) return;
	if(toggle != 1 && toggle != 2 && toggle != 3 && toggle != 4) return;
	Channel *C = &S->region[reg][hilo].ee;

	TH1D *histo;
	if(toggle == 1) histo = (TH1D*)C->sst_origin->Clone();
	if(toggle == 2) histo = (TH1D*)C->ssl_origin->Clone();
	if(toggle == 3) histo = (TH1D*)C->zt_origin->Clone();
	if(toggle == 4) histo = (TH1D*)C->zl_origin->Clone();

	if(histo->Integral() == 0.) return;
	histo->Scale(100./histo->Integral());

	fOUTSTREAM << setw(15) << S->sname;
	fOUTSTREAM << "     |";
	for(size_t i = 1; i <= 11; ++i)	fOUTSTREAM << setw(6)  << setprecision(3) << Form(" %6.2f%%", histo->GetBinContent(i)) << " |";
	fOUTSTREAM << setw(6)  << setprecision(3) << Form(" %6.2f%%", histo->GetBinContent(15)) << " |" << endl;
}

void SSDLPlotter::print2MuOriginsFromSample(Sample *S, int toggle, gRegion reg, gHiLoSwitch hilo){
	if(S->datamc == 0) return;
	if(toggle != 0 && toggle != 1 && toggle != 2) return;
	Channel *C = &S->region[reg][hilo].mm;

	TH2D *histo2d;
	if(toggle == 0) histo2d = C->nt00_origin;
	if(toggle == 1) histo2d = C->nt10_origin;
	if(toggle == 2) histo2d = C->nt11_origin;

	TH1D *histo = histo2d->ProjectionX();
	if(histo->Integral() == 0.) return;
	histo->Scale(100./histo->Integral());

	fOUTSTREAM << setw(15) << S->sname;

	// first muon
	fOUTSTREAM << " mu1 |";
	for(size_t i = 1; i <= 9; ++i)	fOUTSTREAM << setw(6)  << setprecision(3) << Form(" %6.2f%%", histo->GetBinContent(i)) << " |";
	fOUTSTREAM << setw(6)  << setprecision(3) << Form(" %6.2f%%", histo->GetBinContent(15)) << " |" << endl;

	// second muon
	fOUTSTREAM << "                mu2 |";
	histo = histo2d->ProjectionY();
	histo->Scale(100./histo->Integral());
	for(size_t i = 1; i <= 9; ++i)	fOUTSTREAM << setw(6)  << setprecision(3) << Form(" %6.2f%%", histo->GetBinContent(i)) << " |";
	fOUTSTREAM << setw(6)  << setprecision(3) << Form(" %6.2f%%", histo->GetBinContent(15)) << " |" << endl;			
}
void SSDLPlotter::print2ElOriginsFromSample(Sample *S, int toggle, gRegion reg, gHiLoSwitch hilo){
	if(S->datamc == 0) return;
	if(toggle != 0 && toggle != 1 && toggle != 2) return;
	Channel *C = &S->region[reg][hilo].ee;

	TH2D *histo2d;
	if(toggle == 0) histo2d = C->nt00_origin;
	if(toggle == 1) histo2d = C->nt10_origin;
	if(toggle == 2) histo2d = C->nt11_origin;

	TH1D *histo = histo2d->ProjectionX();
	if(histo->Integral() == 0.) return;
	histo->Scale(100./histo->Integral());

	fOUTSTREAM << setw(15) << S->sname;

	// first electron
	fOUTSTREAM << " el1 |";
	for(size_t i = 1; i <= 11; ++i)	fOUTSTREAM << setw(6)  << setprecision(3) << Form(" %6.2f%%", histo->GetBinContent(i)) << " |";
	fOUTSTREAM << setw(6)  << setprecision(3) << Form(" %6.2f%%", histo->GetBinContent(15)) << " |" << endl;

	// second electron
	fOUTSTREAM << "                el2 |";
	histo = histo2d->ProjectionY();
	histo->Scale(100./histo->Integral());
	for(size_t i = 1; i <= 11; ++i)	fOUTSTREAM << setw(6)  << setprecision(3) << Form(" %6.2f%%", histo->GetBinContent(i)) << " |";
	fOUTSTREAM << setw(6)  << setprecision(3) << Form(" %6.2f%%", histo->GetBinContent(15)) << " |" << endl;			
}
void SSDLPlotter::printEMuOriginsFromSample(Sample *S, int toggle, gRegion reg, gHiLoSwitch hilo){
	if(S->datamc == 0) return;
	if(toggle != 0 && toggle != 1 && toggle != 2 && toggle != 10) return;
	Channel *C = &S->region[reg][hilo].em;

	TH2D *histo2d;
	if(toggle == 0)  histo2d = C->nt00_origin;
	if(toggle == 1)  histo2d = C->nt10_origin;
	if(toggle == 10) histo2d = C->nt01_origin;
	if(toggle == 2)  histo2d = C->nt11_origin;

	TH1D *histo = histo2d->ProjectionX();
	if(histo->Integral() == 0.) return;
	histo->Scale(100./histo->Integral());

	// muon
	fOUTSTREAM << setw(15) << S->sname;
	fOUTSTREAM << " mu  |         |         |";
	for(size_t i = 1; i <= 9; ++i)	fOUTSTREAM << setw(6)  << setprecision(3) << Form(" %6.2f%%", histo->GetBinContent(i)) << " |";
	fOUTSTREAM << setw(6)  << setprecision(3) << Form(" %6.2f%%", histo->GetBinContent(15)) << " |" << endl;

	// electron
	fOUTSTREAM << "                el  |";
	histo = histo2d->ProjectionY();
	histo->Scale(100./histo->Integral());
	for(size_t i = 1; i <= 11; ++i)	fOUTSTREAM << setw(6)  << setprecision(3) << Form(" %6.2f%%", histo->GetBinContent(i)) << " |";
	fOUTSTREAM << setw(6)  << setprecision(3) << Form(" %6.2f%%", histo->GetBinContent(15)) << " |" << endl;			
}

void SSDLPlotter::printOriginSummary(vector<int> samples, int toggle, gChannel chan, gRegion reg, gHiLoSwitch hilo){
	if(toggle != 1 && toggle != 2 && toggle != 3 && toggle != 4) return;
	TH1D *histosum = new TH1D("SST_Origin_Sum", "SSTOrigin",  15, 0, 15);
	histosum->Sumw2();
	for(size_t i = 0; i < samples.size(); ++i){
		Sample *S = fSamples[samples[i]];
		Channel *C;
		if(chan == Muon) C = &S->region[reg][hilo].mm;
		if(chan == Elec) C = &S->region[reg][hilo].ee;

		TH1D *histo;
		if(toggle == 1) histo = (TH1D*)C->sst_origin->Clone();
		if(toggle == 2) histo = (TH1D*)C->ssl_origin->Clone();
		if(toggle == 3) histo = (TH1D*)C->zt_origin->Clone();
		if(toggle == 4) histo = (TH1D*)C->zl_origin->Clone();

		float scale = fLumiNorm / S->lumi;
		//histosum->Add(histo, scale);
		if (!S->datamc == 2) { cout << S->name<< endl; histosum->Add(histo, scale);}
	}
	histosum->Scale(100./histosum->Integral());
	fOUTSTREAM << " Weighted Sum       |";
	int ncols = 9;
	if(chan == Elec) ncols = 11;
	for(size_t i = 1; i <= ncols; ++i)	fOUTSTREAM << setw(6)  << setprecision(3) << Form(" %6.2f%%", histosum->GetBinContent(i)) << " |";
	fOUTSTREAM << setw(6)  << setprecision(3) << Form(" %6.2f%%", histosum->GetBinContent(15)) << " |" << endl;
	delete histosum;
}
void SSDLPlotter::printOriginSummary2L(vector<int> samples, int toggle, gChannel chan, gRegion reg, gHiLoSwitch hilo){
	if(toggle != 0 && toggle != 1 && toggle != 10 && toggle != 2) return;
	TH1D *histosum1 = new TH1D("SST_Origin_Sum1", "SSTOrigin",  15, 0, 15);
	TH1D *histosum2 = new TH1D("SST_Origin_Sum2", "SSTOrigin",  15, 0, 15);
	histosum1->Sumw2();
	histosum2->Sumw2();
	for(size_t i = 0; i < samples.size(); ++i){
		Sample *S = fSamples[samples[i]];
		Channel *C;
		if(chan == Muon)     C = &S->region[reg][hilo].mm;
		if(chan == Elec) C = &S->region[reg][hilo].ee;
		if(chan == ElMu)      C = &S->region[reg][hilo].em;

		TH2D *histo2d;
		if(toggle == 0)  histo2d = C->nt00_origin;
		if(toggle == 1)  histo2d = C->nt10_origin;
		if(toggle == 10) histo2d = C->nt01_origin;
		if(toggle == 2)  histo2d = C->nt11_origin;

		float scale = fLumiNorm / S->lumi;
		if (S->datamc != 2) {
			histosum1->Add(histo2d->ProjectionX(), scale);
			histosum2->Add(histo2d->ProjectionY(), scale);
		}
	}
	histosum1->Scale(100./histosum1->Integral());
	histosum2->Scale(100./histosum2->Integral());
	if(chan == Muon)     fOUTSTREAM << " Weighted Sum: mu1  |";
	if(chan == ElMu)      fOUTSTREAM << " Weighted Sum: muon |         |         |";
	if(chan == Elec) fOUTSTREAM << " Weighted Sum: el1  |";
	int ncols = 9;
	if(chan == Elec) ncols = 11;
	for(size_t i = 1; i <= ncols; ++i)	fOUTSTREAM << setw(6)  << setprecision(3) << Form(" %6.2f%%", histosum1->GetBinContent(i)) << " |";
	fOUTSTREAM << setw(6)  << setprecision(3) << Form(" %6.2f%%", histosum1->GetBinContent(15)) << " |" << endl;

	if(chan == Muon)     fOUTSTREAM << "               mu2  |";
	if(chan == ElMu)      fOUTSTREAM << "               ele  |";
	if(chan == Elec) fOUTSTREAM << "               el2  |";
	if(chan == ElMu) ncols = 11;	
	for(size_t i = 1; i <= ncols; ++i)	fOUTSTREAM << setw(6)  << setprecision(3) << Form(" %6.2f%%", histosum2->GetBinContent(i)) << " |";
	fOUTSTREAM << setw(6)  << setprecision(3) << Form(" %6.2f%%", histosum2->GetBinContent(15)) << " |" << endl;
	delete histosum1;
	delete histosum2;
}

//____________________________________________________________________________
void SSDLPlotter::drawTopLine(){
	fLatex->SetTextFont(62);
	fLatex->SetTextSize(0.05);
	fLatex->DrawLatex(0.13,0.92, "CMS Preliminary");	
	fLatex->SetTextFont(42);
	fLatex->SetTextSize(0.04);
	fLatex->DrawLatex(0.70,0.92, Form("L_{int.} = %2.1f fb^{-1}", fLumiNorm/1000.));
	// fLatex->DrawLatex(0.70,0.92, Form("L_{int.} = %4.0f pb^{-1}", fLumiNorm));
	return;
}
void SSDLPlotter::drawDiffCuts(int j){
	fLatex->SetTextFont(42);
	fLatex->SetTextSize(0.03);
	if(j==8)       fLatex->DrawLatex(0.13,0.85, "E_{T}^{miss} > 30 GeV, H_{T} > 200 GeV, N_{Jets} #geq 2");
	if(j>4&&j<8)   fLatex->DrawLatex(0.13,0.85, "E_{T}^{miss} > 30 GeV, H_{T} > 80 GeV, N_{Jets} #geq 2");
	if(j==0||j==4) fLatex->DrawLatex(0.13,0.85, "E_{T}^{miss} > 50 GeV");
	if(j==1)       fLatex->DrawLatex(0.13,0.85, "E_{T}^{miss} > 120 GeV");
	if(j==2)       fLatex->DrawLatex(0.13,0.85, "H_{T} > 200 GeV, N_{Jets} #geq 2");
	if(j==3)       fLatex->DrawLatex(0.13,0.85, "H_{T} > 450 GeV, N_{Jets} #geq 2");
	return;
}

void SSDLPlotter::load_kfacs(TFile * results){
    ifstream IN("msugraSSDL/nlo_kfactors.txt");

    TH2D *kfac_[10];
    char buffer[1000];
    for (int i = 0 ; i< 10; i++){
      kfac_[i]  = new TH2D(Form("kfac_%i", i), Form("kfac_%i", i), gM0bins, gM0min+10, gM0max+10, gM12bins, gM12min+10, gM12max+10);
    }

    while( IN.getline(buffer, 1000, '\n') ){
      if (buffer[0] == '#') continue; // Skip lines commented with '#'
      float p[10];
      float m0_(-1), m12_(-1);
      sscanf(buffer, "%f | %f | %f | %f | %f | %f | %f | %f | %f | %f | %f | %f", &m0_, &m12_, &p[0], &p[1], &p[2], &p[3], &p[4], &p[5], &p[6], &p[7], &p[8], &p[9]);
      if(fVerbose > 1) cout << Form("m0: %4.0f m12: %4.0f (1) %1.5f (2) %1.5f (3) %1.5f (4) %1.5f (5) %1.5f (6) %1.5f (7) %1.5f (8) %1.5f (9) %1.5f (10) %1.5f", m0_, m12_, p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7], p[8], p[9]) << endl;
      for (int i = 0 ; i< 10; i++){
          kfac_[i]->Fill(m0_, m12_, p[i]);
      }
    }
	results->cd();
    for (int i   = 0 ; i< 10; i++){
      kfac_[i]->Write();
    }
}
void SSDLPlotter::load_loxsecs(TFile * results){
	ifstream IN("msugraSSDL/xsec_lo.txt");
	results->cd();
	TH2D * lo_xsec   = new TH2D("lo_xsec", "lo_xsec", gM0bins, gM0min+10, gM0max+10, gM12bins, gM12min+10, gM12max+10);
	
	char buffer[1000];
	while( IN.getline(buffer, 1000, '\n') ){
	  if (buffer[0] == '#') continue; // Skip lines commented with '#'
	  float xsec;
	  float m0_(-1), m12_(-1);
	  sscanf(buffer, "%f | %f | %f", &m0_, &m12_, &xsec);
	  if(fVerbose > 1) cout << Form("m0: %4.0f m12: %4.0f xsec: %10.5f", m0_, m12_, xsec) << endl;
	  lo_xsec->Fill(m0_, m12_, xsec);
	}
	results->cd();
	lo_xsec->Write();
}

void SSDLPlotter::load_nloxsecs(TFile * results) {
	results->cd();
	TH2D * nlo_xsec   = new TH2D("nlo_xsec", "nlo_xsec", gM0bins, gM0min+10, gM0max+10, gM12bins, gM12min+10, gM12max+10);
	
	ifstream IN("msugraSSDL/xsec_nlo.txt");
	char buffer[1000];
	while( IN.getline(buffer, 1000, '\n') ){
		if (buffer[0] == '#') continue; // Skip lines commented with '#'
    	float p[10];
		float m0_(-1), m12_(-1), nloXsec(0);
    	sscanf(buffer, "%f | %f | %f | %f | %f | %f | %f | %f | %f | %f | %f | %f", &m0_, &m12_, &p[0], &p[1], &p[2], &p[3], &p[4], &p[5], &p[6], &p[7], &p[8], &p[9]);
		for (int i=0; i<10; i++) {
			nloXsec += p[i];
		}
		nlo_xsec->Fill(m0_, m12_, nloXsec);
	}
	results->cd();
	nlo_xsec->Write();
}

void SSDLPlotter::load_msugraInfo( const char * filestring){
	// ATTENTION: it is important that m0, m12 and process are known by the InitMC function!! This function is in SSDLDumper now, so no problem

	fC_minMu1pt = 20.;
	fC_minMu2pt = 10.;
	fC_minEl1pt = 20.;
	fC_minEl2pt = 10.;
	fC_minHT    = 450.;
	fC_minMet   = 120.;
	fC_maxHT    = 7000.;
	fC_maxMet   = 7000.;
	fC_minNjets = 2;
	
	fOUTSTREAM.open("msugraSSDL/output.txt", ios::trunc);
	TFile * res_ = new TFile(Form("msugraSSDL/results_HT%3.0f_MET%3.0f_PT%2.0f_%2.0f.root",fC_minHT, fC_minMet, fC_minMu1pt, fC_minMu2pt), "RECREATE", "res_");

	// comment the next two lines in case you already have the necessary files at hand.
	// you can also just let them in, it doesn't really affect the performance much
	SSDLPlotter::load_kfacs(res_);
	SSDLPlotter::load_loxsecs(res_);
	SSDLPlotter::load_nloxsecs(res_);

	TH2D  * nPass_  = new TH2D("msugra_nPass"  , "msugra_nPass"  , gM0bins , gM0min+10 , gM0max+10 , gM12bins , gM12min+10 , gM12max+10);
	TH2D  * wPass_  = new TH2D("msugra_wPass"  , "msugra_wPass"  , gM0bins , gM0min+10 , gM0max+10 , gM12bins , gM12min+10 , gM12max+10);
	TH2D  * wCount_ = new TH2D("msugra_wcount" , "msugra_wcount" , gM0bins , gM0min+10 , gM0max+10 , gM12bins , gM12min+10 , gM12max+10);
	TH2D  * yield_  = new TH2D("msugra_yield"  , "msugra_yield"  , gM0bins , gM0min+10 , gM0max+10 , gM12bins , gM12min+10 , gM12max+10);

	TH2D * kfacs[10];
	for (int i = 0; i < 10; i++) {
		kfacs[i] = (TH2D *) res_->Get(Form("kfac_%i", i));
	}

	TFile * file_ = new TFile(filestring, "READ", "file_"); // example file: "/scratch/mdunser/111111_msugra/msugra_tan10.root"
	if ( file_->IsOpen() ) cout << "File is open " << endl;
	TH2D * count_ = (TH2D *) file_->Get("msugra_count");

	TH2D * kCounts[10];
	for (int i = 0; i < 10; i++) {
		kCounts[i] = (TH2D *) file_->Get(Form("msugra_count_process%i", i+1));
	}
	
	// this is no longer necessary, if we use the Higgs tool, i.e. if we don't need efficiencies but yields
	//for (int m0_ = gM0min+20; m0_<=gM0max; m0_+=20){
	//	for (int m12_ = gM12min+20; m12_<=gM12max; m12_+=20){
	//		int bin = wCount_->FindBin(m0_, m12_);
	//		float binCount = 0.;
	//		float binFill =0.;
	//		for (int i = 0; i < 10; i++) {
	//			float kWeight = kfacs[i]->GetBinContent(bin);
	//			float kCount  = kCounts[i]->GetBinContent(bin);
	//			binCount += kCount;
	//			binFill += kWeight*kCount;
	//		}
	//		wCount_->Fill(m0_, m12_, binFill);
	//	}
	//}
//////--------------------------------------------------------
	ifstream IN("msugraSSDL/xsec_nlo.txt");
	TH2D *nlo_[10];
	char buffer[1000];
	for (int i = 0 ; i< 10; i++){
		nlo_[i]  = new TH2D(Form("nlo_%i", i), Form("nlo_%i", i), gM0bins, gM0min+10, gM0max+10, gM12bins, gM12min+10, gM12max+10);
	}
	
	while( IN.getline(buffer, 1000, '\n') ){
	  if (buffer[0] == '#') continue; // Skip lines commented with '#'
	  float p[10];
	  float m0_(-1), m12_(-1);
	  sscanf(buffer, "%f | %f | %f | %f | %f | %f | %f | %f | %f | %f | %f | %f", &m0_, &m12_, &p[0], &p[1], &p[2], &p[3], &p[4], &p[5], &p[6], &p[7], &p[8], &p[9]);
	  if(fVerbose > 1) cout << Form("m0: %4.0f m12: %4.0f (1) %1.5f (2) %1.5f (3) %1.5f (4) %1.5f (5) %1.5f (6) %1.5f (7) %1.5f (8) %1.5f (9) %1.5f (10) %1.5f", m0_, m12_, p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7], p[8], p[9]) << endl;
	  for (int i = 0 ; i< 10; i++){
	      nlo_[i]->Fill(m0_, m12_, p[i]);
	  }
	}
//////--------------------------------------------------------

	TTree * tree_= (TTree *) file_->Get("Analysis");
	tree_->ResetBranchAddresses();

	InitMC(tree_);
	double tot_events = tree_->GetEntriesFast();
	cout << "Total Number of entries: " << tot_events << endl;
	int n_tot = 0;
	float tightTot(0);
	float signalTot(0);
	float nEE(0), nEM(0), nMM(0);

	for (Long64_t jentry=0; jentry<tree_->GetEntriesFast();jentry++) {
		tree_->GetEntry(jentry);
		printProgress(jentry, tot_events, "CMSSM Scan");

		int mu1(-1), mu2(-1);
		if( isSSLLMuEvent(mu1, mu2) ){ // Same-sign loose-loose di muon event
			if(isTightMuon(mu1) &&  isTightMuon(mu2) ){ // Tight-tight
				tightTot++;
				//if ( !(tree_->GetLeaf("IsSignalMuon")->GetValue(mu1) == 1 && tree_->GetLeaf("IsSignalMuon")->GetValue(mu2) == 1) ) continue;
				if ( IsSignalMuon[mu1] != 1 || IsSignalMuon[mu2] != 1 ) continue;
				signalTot++;
				fOUTSTREAM << Form("MuMu - ev %11.0d - m0 %4.0f - m12 %4.0f - process %2.0i - HT(#J/#bJ) %6.2f(%1d/%1d) MET %6.2f Pt1 %6.2f Pt2 %6.2f Charge %2d", Event, m0, m12, process, getHT(), getNJets(), getNBTags(), pfMET, MuPt[mu1], MuPt[mu2], MuCharge[mu1]) << endl ;
				n_tot++;
				int bin = nlo_[process-1]->FindBin(m0, m12);
				float weight = fLumiNorm * nlo_[process-1]->GetBinContent(bin)/ kCounts[process-1]->GetBinContent(bin);
				fOUTSTREAM << " m0: " << m0 << " m12: " << m12 << " process " << process << " weight: " << weight << endl;
				nPass_->Fill(m0, m12);
				yield_->Fill(m0, m12, weight);
				nMM++;
				continue;
			}
			resetHypLeptons();
		}
		int mu(-1), el(-1);
		if( isSSLLElMuEvent(mu, el) ){
			if(  isTightElectron(el) &&  isTightMuon(mu) ){ // Tight-tight
				tightTot++;
				//if ( !(tree_->GetLeaf("IsSignalElectron")->GetValue(el) == 1 && tree_->GetLeaf("IsSignalMuon")->GetValue(mu) == 1) ) continue;
				if ( IsSignalMuon[mu] != 1 || IsSignalElectron[el] != 1 ) continue;
				signalTot++;
				fOUTSTREAM << Form("ElMu - ev %11.0d - m0 %4.0f - m12 %4.0f - process %2.0i - HT(#J/#bJ) %6.2f(%1d/%1d) MET %6.2f Pt1 %6.2f Pt2 %6.2f Charge %2d", Event, m0, m12, process, getHT(), getNJets(), getNBTags(), pfMET, MuPt[mu], ElPt[el], MuCharge[mu]) << endl ;
				n_tot++;
				int bin = nlo_[process-1]->FindBin(m0, m12);
				float weight = fLumiNorm * nlo_[process-1]->GetBinContent(bin)/ kCounts[process-1]->GetBinContent(bin);
				fOUTSTREAM << " m0 " << m0 << " m12: " << m12 << " process " << process << " weight: " << weight << endl;
				nPass_->Fill(m0, m12);
				yield_->Fill(m0, m12, weight);
				nEM++;
				continue;
			}
			resetHypLeptons();
		}
		int el1(-1), el2(-1);
		if( isSSLLElEvent(el1, el2) ){
			if(  isTightElectron(el1) &&  isTightElectron(el2) ){ // Tight-tight
				tightTot++;
				//if ( !(tree_->GetLeaf("IsSignalElectron")->GetValue(el1) == 1 && tree_->GetLeaf("IsSignalElectron")->GetValue(el2) == 1) ) continue;
				if ( IsSignalElectron[el1] != 1 || IsSignalElectron[el2] != 1 ) continue;
				signalTot++;
				fOUTSTREAM << Form("ElEl - ev %11.0d - m0 %4.0f - m12 %4.0f - process %2.0i - HT(#J/#bJ) %6.2f(%1d/%1d) MET %6.2f Pt1 %6.2f Pt2 %6.2f Charge %2d", Event, m0, m12, process, getHT(), getNJets(), getNBTags(), pfMET, ElPt[el1], ElPt[el2], ElCharge[el1]) << endl ;
				n_tot++;
				int bin = nlo_[process-1]->FindBin(m0, m12);
				float weight = fLumiNorm * nlo_[process-1]->GetBinContent(bin)/ kCounts[process-1]->GetBinContent(bin);
				fOUTSTREAM << " m0 " << m0 << " m12: " << m12 << " process " << process << " weight: " << weight << endl;
				nPass_->Fill(m0, m12);
				yield_->Fill(m0, m12, weight);
				nEE++;
			}
		}
	}
	cout << "Total Number of SS events (with MET > 120 GeV and HT > 450 GeV): " << n_tot << endl;
	cout << "nEE: " << nEE << " nEM: " << nEM << " nMM: " << nMM << endl;
	fOUTSTREAM << "Total Number of SS events (with MET > 120 GeV and HT > 450 GeV): " << n_tot << endl;
	fOUTSTREAM << "Total number of tight pairs: " << tightTot << " total number of signal pairs: " << signalTot << " resulting efficiency: " << signalTot/tightTot << endl;
	fOUTSTREAM << "nEE: " << nEE << " nEM: " << nEM << " nMM: " << nMM << endl;

	TH2D *eff_   = new TH2D("msugra_eff"   , "msugra_eff"  , gM0bins, gM0min+10, gM0max+10, gM12bins, gM12min+10, gM12max+10);
	eff_->Divide(nPass_, count_, 1., 1.);

	res_->cd();

	for (int i = 0 ; i< 10; i++){
	    nlo_[i]->Write();
	}
	yield_->Write();
	count_->Write();
	//wCount_->Write();
	nPass_->Write();
	eff_->Write();

}
