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
static const bool gEWKino = true;
static const bool gRunSMSscan = false;
static const bool gRatiosFromTTbar = false;
static const bool gFullDataBlind = false;

static const float gMMU = 0.1057;
static const float gMEL = 0.0005;
static const float gMZ  = 91.2;

static const double gStatBetaAlpha = 1.;
static const double gStatBetaBeta  = 1.;

TString SSDLPlotter::gHiLoLabel[3] = {"HighPt", "LowPt", "TauChan"};

double SSDLPlotter::gEChMisIDB   = 0.;
double SSDLPlotter::gEChMisIDB_E = 0.;
double SSDLPlotter::gEChMisIDE   = 0.;
double SSDLPlotter::gEChMisIDE_E = 0.;

float SSDLPlotter::gMMTrigScale = 0.872;
float SSDLPlotter::gEMTrigScale = 0.925;
float SSDLPlotter::gEETrigScale = 0.954;

//____________________________________________________________________________
SSDLPlotter::SSDLPlotter(TString configfile):SSDLDumper(configfile){ 
// Default constructor, no samples are set
	fDO_OPT=false;
}
SSDLPlotter::SSDLPlotter(TString outputdir, TString outputfile):SSDLDumper(){
// Explicit constructor with output directory and output file
	fDO_OPT=false;
	setOutputDir(outputdir);
	setOutputFile(outputfile);
}
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

	readVarNames("varnames.dat");
	// MARC readVarNames("anavarnames.dat");
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
	//	fLumiNorm = 3950.; // HPA ICHEP dataset
	//	fLumiNorm = 9185.6; //HCP analysis
	if (gFullDataBlind) fLumiNorm = 19466.;   // Moriond 2013 dataset
	else fLumiNorm = 9100.;	// ttV 2012 analysis
//	fLumiNorm = 9100.;	// ttV 2012 analysis
//	fLumiNorm = 19466.;   // Moriond 2013 dataset

	// Cuts:
	fC_minMu1pt = 20.;
	fC_minMu2pt = 20.;
	fC_minEl1pt = 20.;
	fC_minEl2pt = 20.;
	fC_minHT    = 80.;
	fC_minMet   = 0.;
	fC_maxHT    = 8000.;
	fC_maxMet   = 8000.;
	fC_minNjets = 2;
	
	fC_maxMet_Control = 20.;
	fC_maxMt_Control  = 20.;

	// Prevent root from adding histograms to current file
	TH1::AddDirectory(kFALSE);

	fMCBG.push_back(TTJets);
	fMCBG.push_back(SingleT_t);
	fMCBG.push_back(SingleTbar_t);
	fMCBG.push_back(SingleT_tW);
	fMCBG.push_back(SingleTbar_tW);
	fMCBG.push_back(SingleT_s);
	fMCBG.push_back(SingleTbar_s);
	fMCBG.push_back(WJets);
	fMCBG.push_back(DYJets);
	// MARC fMCBG.push_back(GJets40);
	fMCBG.push_back(GJets200);
	fMCBG.push_back(GJets400);
	fMCBG.push_back(WW);
	fMCBG.push_back(WZ);
	fMCBG.push_back(ZZ);
	// MARC fMCBG.push_back(GVJets);
	fMCBG.push_back(DPSWW);
	fMCBG.push_back(TTbarH);
	fMCBG.push_back(TTbarW);
	fMCBG.push_back(TTbarZ);
	fMCBG.push_back(TTbarG);
	fMCBG.push_back(TbZ);
	fMCBG.push_back(WpWp);
	fMCBG.push_back(WmWm);
	fMCBG.push_back(WWZ);
	fMCBG.push_back(WZZ);
	fMCBG.push_back(WWG);
	fMCBG.push_back(WWW);
	fMCBG.push_back(ZZZ);
	// QCD:
	fMCBG.push_back(QCDMuEnr15);
	fMCBG.push_back(QCDEM20);
	fMCBG.push_back(QCDEM30);
	fMCBG.push_back(QCDEM80);
	fMCBG.push_back(QCDEM170);
	fMCBG.push_back(QCDEM250);
	fMCBG.push_back(QCD80);
	fMCBG.push_back(QCD120);
	fMCBG.push_back(QCD170);
	fMCBG.push_back(QCD300);
	fMCBG.push_back(QCD470);
	fMCBG.push_back(QCD600);
	fMCBG.push_back(QCD800);

	fMCBGNoQCDNoGJets.push_back(TTJets);
	fMCBGNoQCDNoGJets.push_back(SingleT_t);
	fMCBGNoQCDNoGJets.push_back(SingleTbar_t);
	fMCBGNoQCDNoGJets.push_back(SingleT_tW);
	fMCBGNoQCDNoGJets.push_back(SingleTbar_tW);
	fMCBGNoQCDNoGJets.push_back(SingleT_s);
	fMCBGNoQCDNoGJets.push_back(SingleTbar_s);
	fMCBGNoQCDNoGJets.push_back(WJets);
	fMCBGNoQCDNoGJets.push_back(DYJets);
	fMCBGNoQCDNoGJets.push_back(WW);
	fMCBGNoQCDNoGJets.push_back(WZ);
	fMCBGNoQCDNoGJets.push_back(ZZ);
	fMCBGNoQCDNoGJets.push_back(DPSWW);
	fMCBGNoQCDNoGJets.push_back(TTbarH);
	fMCBGNoQCDNoGJets.push_back(TTbarW);
	fMCBGNoQCDNoGJets.push_back(TTbarZ);
	fMCBGNoQCDNoGJets.push_back(TTbarG);
	fMCBGNoQCDNoGJets.push_back(TbZ);
	fMCBGNoQCDNoGJets.push_back(WpWp);
	fMCBGNoQCDNoGJets.push_back(WmWm);
	fMCBGNoQCDNoGJets.push_back(WWZ);
	fMCBGNoQCDNoGJets.push_back(WZZ);
	fMCBGNoQCDNoGJets.push_back(WWG);
	fMCBGNoQCDNoGJets.push_back(WWW);
	fMCBGNoQCDNoGJets.push_back(ZZZ);

	// input samples for closure tests
	fClosureSamples.push_back(TTJets);
	fClosureSamples.push_back(SingleT_t);
	fClosureSamples.push_back(SingleTbar_t);
	fClosureSamples.push_back(SingleT_tW);
	fClosureSamples.push_back(SingleTbar_tW);
	fClosureSamples.push_back(SingleT_s);
	fClosureSamples.push_back(SingleTbar_s);
	fClosureSamples.push_back(DYJets);
	fClosureSamples.push_back(WW);
	fClosureSamples.push_back(WJets);
	
	fTTJets.push_back(TTJets);
	
	fWJets.push_back(WJets);
	
	fDYJets.push_back(DYJets);
	
	fMCBGMuEnr.push_back(TTJets);
	fMCBGMuEnr.push_back(SingleT_t);
	fMCBGMuEnr.push_back(SingleTbar_t);
	fMCBGMuEnr.push_back(SingleT_tW);
	fMCBGMuEnr.push_back(SingleTbar_tW);
	fMCBGMuEnr.push_back(SingleT_s);
	fMCBGMuEnr.push_back(SingleTbar_s);
	fMCBGMuEnr.push_back(WJets);
	fMCBGMuEnr.push_back(DYJets);
	// MARC fMCBGMuEnr.push_back(GJets40);
	// MARC fMCBGMuEnr.push_back(GJets100);
	fMCBGMuEnr.push_back(GJets200);
	fMCBGMuEnr.push_back(GJets400);
	fMCBGMuEnr.push_back(WW);
	fMCBGMuEnr.push_back(WZ);
	fMCBGMuEnr.push_back(ZZ);
	// MARC fMCBGMuEnr.push_back(GVJets);
	fMCBGMuEnr.push_back(DPSWW);
	fMCBGMuEnr.push_back(TTbarH);
	fMCBGMuEnr.push_back(TTbarW);
	fMCBGMuEnr.push_back(TTbarZ);
	fMCBGMuEnr.push_back(TTbarG);
	fMCBGMuEnr.push_back(TbZ);
	fMCBGMuEnr.push_back(WpWp);
	fMCBGMuEnr.push_back(WmWm);
	fMCBGMuEnr.push_back(WWZ);
	fMCBGMuEnr.push_back(WZZ);
	fMCBGMuEnr.push_back(WWG);
	fMCBGMuEnr.push_back(WWW);
	fMCBGMuEnr.push_back(ZZZ);
	fMCBGMuEnr.push_back(QCD80);
	fMCBGMuEnr.push_back(QCD120);
	fMCBGMuEnr.push_back(QCD170);
	fMCBGMuEnr.push_back(QCD300);
	fMCBGMuEnr.push_back(QCD470);
	fMCBGMuEnr.push_back(QCD600);
	fMCBGMuEnr.push_back(QCD800);
	fMCBGMuEnr.push_back(QCDMuEnr15);
	
	fMCBGEMEnr.push_back(TTJets);
	fMCBGEMEnr.push_back(SingleT_t);
	fMCBGEMEnr.push_back(SingleTbar_t);
	fMCBGEMEnr.push_back(SingleT_tW);
	fMCBGEMEnr.push_back(SingleTbar_tW);
	fMCBGEMEnr.push_back(SingleT_s);
	fMCBGEMEnr.push_back(SingleTbar_s);
	fMCBGEMEnr.push_back(WJets);
	fMCBGEMEnr.push_back(DYJets);
	// MARC fMCBGEMEnr.push_back(GJets40);
	// MARC fMCBGEMEnr.push_back(GJets100);
	fMCBGEMEnr.push_back(GJets200);
	fMCBGEMEnr.push_back(GJets400);
	fMCBGEMEnr.push_back(WW);
	fMCBGEMEnr.push_back(WZ);
	fMCBGEMEnr.push_back(ZZ);
	// MARC fMCBGEMEnr.push_back(GVJets);
	fMCBGEMEnr.push_back(DPSWW);
	fMCBGEMEnr.push_back(TTbarH);
	fMCBGEMEnr.push_back(TTbarW);
	fMCBGEMEnr.push_back(TTbarZ);
	fMCBGEMEnr.push_back(TTbarG);
	fMCBGEMEnr.push_back(TbZ);
	fMCBGEMEnr.push_back(WpWp);
	fMCBGEMEnr.push_back(WmWm);
	fMCBGEMEnr.push_back(WWZ);
	fMCBGEMEnr.push_back(WZZ);
	fMCBGEMEnr.push_back(WWG);
	fMCBGEMEnr.push_back(WWW);
	fMCBGEMEnr.push_back(ZZZ);

	fMCBGEMEnr.push_back(QCDEM20);
	fMCBGEMEnr.push_back(QCDEM30);
	fMCBGEMEnr.push_back(QCDEM80);
	fMCBGEMEnr.push_back(QCDEM170);
	fMCBGEMEnr.push_back(QCDEM250);
	fMCBGEMEnr.push_back(QCD80);
	fMCBGEMEnr.push_back(QCD120);
	fMCBGEMEnr.push_back(QCD170);
	fMCBGEMEnr.push_back(QCD300);
	fMCBGEMEnr.push_back(QCD470);
	fMCBGEMEnr.push_back(QCD600);
	fMCBGEMEnr.push_back(QCD800);
	fMuEnr.push_back(QCDMuEnr15);
	
	fEMEnr.push_back(QCDEM20);
	fEMEnr.push_back(QCDEM30);
	fEMEnr.push_back(QCDEM80);
	fEMEnr.push_back(QCDEM170);
	fEMEnr.push_back(QCDEM250);

	fMCBGMuEnrSig = fMCBGMuEnr;


	fMCRareSM.push_back(ZZ);
	// MARC fMCRareSM.push_back(GVJets);
	fMCRareSM.push_back(DPSWW);
	if(gEWKino)  fMCRareSM.push_back(TTbarH);
	if(gEWKino)  fMCRareSM.push_back(TTbarW);
	if(gEWKino)  fMCRareSM.push_back(TTbarZ);
	fMCRareSM.push_back(TTbarG);
	fMCRareSM.push_back(TbZ);
	fMCRareSM.push_back(WpWp);
	fMCRareSM.push_back(WmWm);
	fMCRareSM.push_back(WWZ);
	fMCRareSM.push_back(WZZ);
	fMCRareSM.push_back(WWG);
	fMCRareSM.push_back(WWW);
	fMCRareSM.push_back(ZZZ);

	fMuData    .push_back(DoubleMu1);
	fMuData    .push_back(DoubleMu1a);
	fMuData    .push_back(DoubleMu2);
	fMuData    .push_back(DoubleMu3);
	fMuData    .push_back(DoubleMu4);
	fMuData    .push_back(DoubleMu5);
	fMuData    .push_back(DoubleMu5a);
	// fMuHadData .push_back(MuHad1);
	// fMuHadData .push_back(MuHad2);
	fEGData    .push_back(DoubleEle1);
	fEGData    .push_back(DoubleEle1a);
	fEGData    .push_back(DoubleEle2);
	fEGData    .push_back(DoubleEle3);
	fEGData    .push_back(DoubleEle4);
	fEGData    .push_back(DoubleEle5);
	fEGData    .push_back(DoubleEle5a);
	// fEleHadData.push_back(EleHad1);
	// fEleHadData.push_back(EleHad2);
	fMuEGData  .push_back(MuEG1);
	fMuEGData  .push_back(MuEG1a);
	fMuEGData  .push_back(MuEG2);
	fMuEGData  .push_back(MuEG3);
	fMuEGData  .push_back(MuEG4);
	fMuEGData  .push_back(MuEG5);
	fMuEGData  .push_back(MuEG5a);

	fHighPtData.push_back(DoubleMu1);
	fHighPtData.push_back(DoubleMu1a);
	fHighPtData.push_back(DoubleMu2);
	fHighPtData.push_back(DoubleMu3);
	fHighPtData.push_back(DoubleMu4);
	fHighPtData.push_back(DoubleMu5);
	fHighPtData.push_back(DoubleMu5a);
	fHighPtData.push_back(DoubleEle1);
	fHighPtData.push_back(DoubleEle1a);
	fHighPtData.push_back(DoubleEle2);
	fHighPtData.push_back(DoubleEle3);
	fHighPtData.push_back(DoubleEle4);
	fHighPtData.push_back(DoubleEle5);
	fHighPtData.push_back(DoubleEle5a);
	fHighPtData.push_back(MuEG1);
	fHighPtData.push_back(MuEG1a);
	fHighPtData.push_back(MuEG2);
	fHighPtData.push_back(MuEG3);
	fHighPtData.push_back(MuEG4);
	fHighPtData.push_back(MuEG5);
	fHighPtData.push_back(MuEG5a);

	// fLowPtData.push_back(MuHad1);
	// fLowPtData.push_back(MuHad2);
	// fLowPtData.push_back(EleHad1);
	// fLowPtData.push_back(EleHad2);
	// fLowPtData.push_back(MuEG1);
	// fLowPtData.push_back(MuEG2);
	// fLowPtData.push_back(MuEG3);
	// fLowPtData.push_back(MuEG4);
}
void SSDLPlotter::doSMSscans(TString region, TString file){
        // This macro runs over the the SMS scans 
	// TString pathtofile = ""; //"dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/";
	// TString filename   = "/mdunser/SSDLTrees/2012/Oct15/SMS-TChiSlepSnu_Mchargino-100to1000_mLSP-0to975_8TeV-Pythia6Z_4.root";  
	// TString filename   = "/shome//mdunser/xsecfiles/output_38.root";  
	
	cout << "Running Over region: " << region << endl;
	cout << "On file: " << file << endl;
	// scanSMS(pathtofile + filename, gRegion[region]);
	scanModelGeneric(file , gRegion[region]);
	
}
void SSDLPlotter::doAnalysis(){
	// sandBox();
	// pythiaMadgraph(true);
	// pythiaMadgraph(false);
        // return;
       
  //	if (gRunSMSscan) return; //DO NOT RUN THE ANALYSIS IF RUNNING THE SCAN
        if(readHistos(fOutputFileName) != 0) return;
	fillRatios(fMuData, fEGData, 0);
	fillRatios(fMCBGMuEnr, fMCBGEMEnr, 1);
	storeWeightedPred(gRegion[gBaseRegion]);


//SANTI	makeOSSFEfficiency();
//SANTI	// ------------------- the following is for testing and closure things
	// fill fake ratios to run closure tests
//	fillRatios(fMuData, fEGData, 0);
//	fillRatios(fMuEnr, fEMEnr, fDYJets, fDYJets, 1);
//	fillRatios(fMuData, fEGData, 0);
//	fillRatios(fMuEnr, fEMEnr, fMCBGMuEnr, fMCBGEMEnr, 0);
//	fillRatios(fMuEnr, fEMEnr, fMCBGMuEnr, fMCBGEMEnr, 1);
//	fillRatios(fMCBGMuEnr, fMCBGEMEnr, 0);
//	fillRatios(fMCBGMuEnr, fMCBGEMEnr, 1);
//	fillRatios(fMuData, fEGData, 0);
//	fillRatios(fMuData, fEGData, 1);
//	fillRatios(fTTJets, fTTJets, 0);
//	fillRatios(fTTJets, fTTJets, 1);
//	storeWeightedPred(gRegion[gBaseRegion]);
	

	// makeROCCurve(); return;
	// makePileUpPlots(); // loops on all data!
	
	// printCutFlows(fOutputDir + "CutFlow.txt");

	// makeOriginPlots(gRegion["Baseline"]);
	// printOrigins(gRegion["Baseline"]);
	// makeOriginPlots(HT0MET120);
	// printOrigins(HT0MET120);
	// makeOriginPlots(HT0MET120lV);
	// printOrigins(HT0MET120lV);

	// makeMuIsolationPlots(false); // if true, loops on TTbar sample
	// makeElIsolationPlots(false); // if true, loops on TTbar sample
	// makeElIdPlots();
	// makeNT2KinPlots(false);
	// makeNT2KinPlots(true);
	//makeMETvsHTPlot(fMuData, fEGData, fMuEGData, HighPt);

	// makeMETvsHTPlotPRL();
	//	makeMETvsHTPlot0HT();
	// makeMETvsHTPlotTau();

//	makeRatioPlots(Muon);
//	makeRatioPlots(Elec);
//	make2DRatioPlots(Muon);
//	make2DRatioPlots(Elec);
//	// // makeNTightLoosePlots(Muon);
//	// // makeNTightLoosePlots(Elec);
//	// 
//	makeFRvsPtPlots(Muon, SigSup);
//	makeFRvsPtPlots(Elec, SigSup);
//	makeFRvsPtPlots(Muon, ZDecay);
//	makeFRvsPtPlots(Elec, ZDecay);
//	makeFRvsNVPlots(Muon, ZDecay); 
//	makeFRvsNVPlots(Elec, ZDecay);
//	makeFRvsEtaPlots(Muon);
//	makeFRvsEtaPlots(Elec);
//	makeChMidvsPtPlots();
//
//	makeAllClosureTestsTTW();
//	makeAllIntPredictions();
//	makeAllClosureTests();
//
	// makeDiffPrediction();
	// makeTTWDiffPredictions();
	// makeTTWIntPredictions();
	// printAllYieldTables();
	

// 	makeTTWIntPredictionsSigEvent(285., 8000., 0., 8000., 3, 1, 1, 40., 40., 0, true);
 	makeTTWIntPredictionsSigEvent(268.235, 8000., 0., 8000., 3, 1, 1, 35.5956, 35.5956, 0, true);
// 	makeTTWIntPredictionsSigEvent(285., 8000., 0., 8000., 3, 1, 1, 40., 40.,-1, true);
// 	makeTTWIntPredictionsSigEvent(285., 8000., 0., 8000., 3, 1, 1, 40., 40.,+1, true);
//	makePredictionSignalEvents(285., 8000., 0., 8000., 3, 1, 1, 40., 40.,-1, true);
//	makePredictionSignalEvents(285., 8000., 0., 8000., 3, 1, 1, 40., 40.,+1, true);
//	makePredictionSignalEvents(285., 8000., 0., 8000., 3, 1, 1, 40., 40., 0, true);
//	makePredictionSignalEvents(100., 8000., 0., 8000., 3, 1, 1, 55., 35., 0, true);
//	makePredictionSignalEvents(175., 8000., 0., 7000., 3, 1, 1, 35., 35., 0, true);
	//makePredictionSignalEvents(280., 8000., 0., 7000., 3, 1, 1, 40., 40., 0, true);
	// makeRelIsoTTSigPlots();
	
//	makeFakeGenIDTables();
//	makeMIDIsolationPlots(Muon, SigSup);
//	makeMIDIsolationPlots(Muon, Sig);
}

//____________________________________________________________________________
void SSDLPlotter::sandBox(){
	FakeRatios *FR = new FakeRatios();
	FR->setIsMC(true);
	FR->setNToyMCs(100);
	FR->setAddESyst(0.0);
	
	FR->setNGen(59574249);
	
	FR->setMFRatio(0.085, 5.8e-05);
	FR->setMPRatio(0.930, 0.0016);
	FR->setEFRatio(0.190, 0.00011);
	FR->setEPRatio(0.900, 0.0022);

	// TTJets || 196 | 4600 | 768 || 367 | 905 | 3972 | 676 || 179 | 741 | 109 || 
	FR->setMMNtl(196, 4600, 768);
	FR->setEMNtl(367, 905, 3972, 676);
	FR->setEENtl(179, 741, 109);

	
	cout << Form("  MM || %9.3f Â± %7.3f (stat) | %9.3f Â± %7.3f (stat) | %9.3f Â± %7.3f (stat) |                            |",
	FR->getMMNpp(), FR->getMMNppEStat(), FR->getMMNpf(), FR->getMMNpfEStat(), FR->getMMNff(), FR->getMMNffEStat()) << endl;
	cout << Form("  EM || %9.3f Â± %7.3f (stat) | %9.3f Â± %7.3f (stat) | %9.3f Â± %7.3f (stat) | %9.3f Â± %7.3f (stat) |",
	FR->getEMNpp(), FR->getEMNppEStat(), FR->getEMNpf(), FR->getEMNpfEStat(), FR->getEMNfp(), FR->getEMNfpEStat(), FR->getEMNff(), FR->getEMNffEStat()) << endl;
	cout << Form("  EE || %9.3f Â± %7.3f (stat) | %9.3f Â± %7.3f (stat) | %9.3f Â± %7.3f (stat) |                            |",
	FR->getEENpp(), FR->getEENppEStat(), FR->getEENpf(), FR->getEENpfEStat(), FR->getEENff(), FR->getEENffEStat()) << endl;
	
	delete FR;
}
void SSDLPlotter::pythiaMadgraph(bool pythia){

	fOutputSubDir = "IntPredictions/";
	TString jvString = "";

	const float RareESyst  = 0.5;
	const float RareESyst2 = RareESyst*RareESyst;
	
	const float FakeESyst  = 0.5;
	const float FakeESyst2 = FakeESyst*FakeESyst;

	const float WZESyst  = 0.2;
	const float WZESyst2 = WZESyst*WZESyst;

	TFile * pFile;
	if (pythia) pFile = TFile::Open("/shome/mdunser/workspace/CMSSW_5_2_5/src/DiLeptonAnalysis/NTupleProducer/macros/plots/Jul25_WZ__WZ_invZveto_jet70/WZTo3LNuPythia_Yields.root");
	else pFile = TFile::Open("/shome/mdunser/workspace/CMSSW_5_2_5/src/DiLeptonAnalysis/NTupleProducer/macros/plots/Sep10_WZ_invZveto_jet70/WZTo3LNu_dataAndMadgraph_Yields.root");

	// ///////////////////////////////////////////////////////////////////////////////////
	// // RATIOS /////////////////////////////////////////////////////////////////////////
	// ///////////////////////////////////////////////////////////////////////////////////
	// float mufratio_data(0.),  mufratio_data_e(0.);
	// float mupratio_data(0.),  mupratio_data_e(0.);
	// float elfratio_data(0.),  elfratio_data_e(0.);
	// float elpratio_data(0.),  elpratio_data_e(0.);

	// calculateRatio(fMuData, Muon, SigSup, mufratio_data, mufratio_data_e);
	// calculateRatio(fMuData, Muon, ZDecay, mupratio_data, mupratio_data_e);

	// calculateRatio(fEGData, Elec, SigSup, elfratio_data, elfratio_data_e);
	// calculateRatio(fEGData, Elec, ZDecay, elpratio_data, elpratio_data_e);

	// ///////////////////////////////////////////////////////////////////////////////////
	// // OBSERVATIONS ///////////////////////////////////////////////////////////////////
	// ///////////////////////////////////////////////////////////////////////////////////
	// float nt2_mm(0.), nt10_mm(0.), nt0_mm(0.);
	// float nt2_em(0.), nt10_em(0.), nt01_em(0.), nt0_em(0.);
	// float nt2_ee(0.), nt10_ee(0.), nt0_ee(0.);

	// // FR Predictions from event-by-event weights (pre stored)
	// float npp_mm(0.), npf_mm(0.), nff_mm(0.);
	// float npp_em(0.), npf_em(0.), nfp_em(0.), nff_em(0.);
	// float npp_ee(0.), npf_ee(0.), nff_ee(0.);

	// // OS yields
	// float nt2_ee_BB_os(0.), nt2_ee_EE_os(0.), nt2_ee_EB_os(0.);
	// float nt2_em_BB_os(0.), nt2_em_EE_os(0.);

	// float nt2_rare_mc_mm(0.),    nt2_rare_mc_em(0.),    nt2_rare_mc_ee(0.);
	// float nt2_rare_mc_mm_e2(0.), nt2_rare_mc_em_e2(0.), nt2_rare_mc_ee_e2(0.);

	// float nt2_wz_mc_mm(0.),    nt2_wz_mc_em(0.),    nt2_wz_mc_ee(0.);
	// float nt2_wz_mc_mm_e2(0.), nt2_wz_mc_em_e2(0.), nt2_wz_mc_ee_e2(0.);

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	TTree *sigtree; getObjectSafe(pFile, "SigEvents", sigtree);
	
	string *sname = 0;
	int flag;
	int   SType, Flavor, TLCat, NJ, NbJ, NbJmed;
	float puweight, pT1, pT2, HT, MET, MT2, SLumi;
	float eta1, eta2, mll;
	int   event, run;

	sigtree->SetBranchAddress("SystFlag", &flag);
	sigtree->SetBranchAddress("Event",    &event);
	sigtree->SetBranchAddress("Run",      &run);
	sigtree->SetBranchAddress("SName",    &sname);
	sigtree->SetBranchAddress("SType",    &SType);
	sigtree->SetBranchAddress("PUWeight", &puweight);
	sigtree->SetBranchAddress("SLumi",    &SLumi);
	sigtree->SetBranchAddress("Flavor",   &Flavor);
	sigtree->SetBranchAddress("pT1",      &pT1);
	sigtree->SetBranchAddress("pT2",      &pT2);
	sigtree->SetBranchAddress("eta1",     &eta1);
	sigtree->SetBranchAddress("eta2",     &eta2);
	sigtree->SetBranchAddress("TLCat",    &TLCat);
	sigtree->SetBranchAddress("HT",       &HT);
	sigtree->SetBranchAddress("MET",      &MET);
	sigtree->SetBranchAddress("MT2",      &MT2);
	sigtree->SetBranchAddress("NJ",       &NJ);
	sigtree->SetBranchAddress("NbJ",      &NbJ);
	sigtree->SetBranchAddress("NbJmed",   &NbJmed);
	sigtree->SetBranchAddress("Mll",      &mll);

	// FakeRatios *FR = new FakeRatios();

	TH1D  * WZmadgraph_allMet_   = new TH1D("WZmadgraph_allMet_"  , "WZmadgraph_allMet_" , 50, 0. , 350.);
	TH1D  * WZmadgraph_allHT_    = new TH1D("WZmadgraph_allHT_"   , "WZmadgraph_allHT_"  , 50, 0. , 350.);
	TH1D  * WZmadgraph_allnJets_ = new TH1D("WZmadgraph_allnJets_", "WZmadgraph_allnJets", 10, 0 , 10);
	TH1D  * WZmadgraph_bnJets_ = new TH1D("WZmadgraph_bnJets_", "WZmadgraph_bnJets", 10, 0 , 10);
	TH1D  * WZmadgraph_bnJetsMet50_ = new TH1D("WZmadgraph_bnJetsMet50_", "WZmadgraph_bnJetsMet50", 10, 0 , 10);
	TH1D  * WZmadgraph_allnJetsMet50_ = new TH1D("WZmadgraph_allnJetsMet50_", "WZmadgraph_allnJetsMet50", 10, 0 , 10);
	TH1D  * WZmadgraph_allnJetsMet100_ = new TH1D("WZmadgraph_allnJetsMet100_", "WZmadgraph_allnJetsMet100", 10, 0 , 10);
	TH1D  * WZmadgraph_allnJetsMet120_ = new TH1D("WZmadgraph_allnJetsMet120_", "WZmadgraph_allnJetsMet120", 10, 0 , 10);
	TH1D  * WZmadgraph_allnJetsMet150_ = new TH1D("WZmadgraph_allnJetsMet150_", "WZmadgraph_allnJetsMet150", 10, 0 , 10);
	TH1D  * WZmadgraph_allnJetsMet200_ = new TH1D("WZmadgraph_allnJetsMet200_", "WZmadgraph_allnJetsMet200", 10, 0 , 10);

	TH1D  * WZpythia_allMet_   = new TH1D("WZpythia_allMet_"  , "WZpythia_allMet_" , 50, 0. , 350.);
	TH1D  * WZpythia_allHT_    = new TH1D("WZpythia_allHT_"   , "WZpythia_allHT_"  , 50, 0. , 350.);
	TH1D  * WZpythia_allnJets_ = new TH1D("WZpythia_allnJets_", "WZpythia_allnJets", 10, 0 , 10);
	TH1D  * WZpythia_bnJets_ = new TH1D("WZpythia_bnJets_", "WZpythia_bnJets", 10, 0 , 10);
	TH1D  * WZpythia_bnJetsMet50_ = new TH1D("WZpythia_bnJetsMet50_", "WZpythia_bnJetsMet50", 10, 0 , 10);
	TH1D  * WZpythia_allnJetsMet50_ = new TH1D("WZpythia_allnJetsMet50_", "WZpythia_allnJetsMet50", 10, 0 , 10);
	TH1D  * WZpythia_allnJetsMet100_ = new TH1D("WZpythia_allnJetsMet100_", "WZpythia_allnJetsMet100", 10, 0 , 10);
	TH1D  * WZpythia_allnJetsMet120_ = new TH1D("WZpythia_allnJetsMet120_", "WZpythia_allnJetsMet120", 10, 0 , 10);
	TH1D  * WZpythia_allnJetsMet150_ = new TH1D("WZpythia_allnJetsMet150_", "WZpythia_allnJetsMet150", 10, 0 , 10);
	TH1D  * WZpythia_allnJetsMet200_ = new TH1D("WZpythia_allnJetsMet200_", "WZpythia_allnJetsMet200", 10, 0 , 10);

	TH1D  * WZdata_allMet_   = new TH1D("WZdata_allMet_"  , "WZdata_allMet_" , 50, 0. , 350.);
	TH1D  * WZdata_allHT_    = new TH1D("WZdata_allHT_"   , "WZdata_allHT_"  , 50, 0. , 350.);
	TH1D  * WZdata_allnJets_ = new TH1D("WZdata_allnJets_", "WZdata_allnJets", 10, 0 , 10);
	TH1D  * WZdata_bnJets_ = new TH1D("WZdata_bnJets_", "WZdata_bnJets", 10, 0 , 10);
	TH1D  * WZdata_bnJetsMet50_ = new TH1D("WZdata_bnJetsMet50_", "WZdata_bnJetsMet50", 10, 0 , 10);
	TH1D  * WZdata_allnJetsMet50_ = new TH1D("WZdata_allnJetsMet50_", "WZdata_allnJetsMet50", 10, 0 , 10);
	TH1D  * WZdata_allnJetsMet100_ = new TH1D("WZdata_allnJetsMet100_", "WZdata_allnJetsMet100", 10, 0 , 10);
	TH1D  * WZdata_allnJetsMet120_ = new TH1D("WZdata_allnJetsMet120_", "WZdata_allnJetsMet120", 10, 0 , 10);
	TH1D  * WZdata_allnJetsMet150_ = new TH1D("WZdata_allnJetsMet150_", "WZdata_allnJetsMet150", 10, 0 , 10);
	TH1D  * WZdata_allnJetsMet200_ = new TH1D("WZdata_allnJetsMet200_", "WZdata_allnJetsMet200", 10, 0 , 10);

	float minPt1 = 20.;
	float minPt2 = 10.;
	float trigScale[3] = {gMMTrigScale, gEMTrigScale, gEETrigScale};
	for( int i = 0; i < sigtree->GetEntries(); i++ ){
		sigtree->GetEntry(i);
		
		// if( flag != systflag ) continue;
		
		if ( mll < 8.) continue;

		gChannel chan = gChannel(Flavor);
		if(chan == ElMu || Flavor == 3){
			if(pT1 > pT2){
				if(pT1 < minPt1) continue;
				if(pT2 < minPt2) continue;
			}
			if(pT1 < pT2){
				if(pT1 < minPt2) continue;
				if(pT2 < minPt1) continue;
			}
		}
		else{
			if(pT1 < minPt1) continue;
			if(pT2 < minPt2) continue;
		}

		// GET ALL DATA EVENTS
		if(SType < 3) {             // 0,1,2 are DoubleMu, DoubleEle, MuEG
			if (Flavor < 3) {
				Sample *S = fSampleMap[TString(*sname)];

				// float npp(0.) , npf(0.) , nfp(0.) , nff(0.);
				// float f1(0.)  , f2(0.)  , p1(0.)  , p2(0.);
				// f1 = getFRatio(chan, pT1, eta1, S->datamc);
				// f2 = getFRatio(chan, pT2, eta2, S->datamc);
				// p1 = getPRatio(chan, pT1, S->datamc);
				// p2 = getPRatio(chan, pT2, S->datamc);
				// if(chan == ElMu){
				// 	f1 = getFRatio(Muon, pT1, eta1, S->datamc);
				// 	f2 = getFRatio(Elec, pT2, eta2, S->datamc);
				// 	p1 = getPRatio(Muon, pT1, S->datamc);
				// 	p2 = getPRatio(Elec, pT2, S->datamc);
				// }
				// Get the weights (don't depend on event selection)
				// npp = FR->getWpp(FakeRatios::gTLCat(TLCat), f1, f2, p1, p2);
				// npf = FR->getWpf(FakeRatios::gTLCat(TLCat), f1, f2, p1, p2);
				// nfp = FR->getWfp(FakeRatios::gTLCat(TLCat), f1, f2, p1, p2);
				// nff = FR->getWff(FakeRatios::gTLCat(TLCat), f1, f2, p1, p2);			

				if (Flavor == 0) {      // MUMU
				//	npp_mm += npp;
				//	npf_mm += (npf+nfp);
				//	nff_mm += nff;
				//	if (TLCat == 0)               nt2_mm++;
				// 	if (TLCat == 1 || TLCat == 2) nt10_mm++;
				//	if (TLCat == 3) {
						WZdata_allMet_   ->Fill(MET);
						WZdata_allHT_    ->Fill(HT);
						WZdata_allnJets_ ->Fill(NJ);
						WZdata_bnJets_   ->Fill(NbJ);
						if (MET >  50. ) { WZdata_allnJetsMet50_ ->Fill(NJ); WZdata_bnJetsMet50_ ->Fill(NbJ); }
						if (MET > 100. ) WZdata_allnJetsMet100_ ->Fill(NJ);
						if (MET > 120. ) WZdata_allnJetsMet120_ ->Fill(NJ);
						if (MET > 200. ) WZdata_allnJetsMet200_ ->Fill(NJ);
						//nt0_mm++;
				//	}
				}
				if (Flavor == 1) {       // E-MU
				//	npp_em += npp;
				//	npf_em += npf;
				//	nfp_em += nfp;
				//	nff_em += nff;
				//	if (TLCat == 0) nt2_em++;
				// 	if (TLCat == 1) nt10_em++;
				// 	if (TLCat == 2) nt01_em++;
				// 	if (TLCat == 3) {
						WZdata_allMet_   ->Fill(MET);
						WZdata_allHT_    ->Fill(HT);
						WZdata_allnJets_ ->Fill(NJ);
						WZdata_bnJets_   ->Fill(NbJ);
						if (MET >  50. ) { WZdata_allnJetsMet50_ ->Fill(NJ); WZdata_bnJetsMet50_ ->Fill(NbJ); }
						if (MET > 100. ) WZdata_allnJetsMet100_ ->Fill(NJ);
						if (MET > 120. ) WZdata_allnJetsMet120_ ->Fill(NJ);
						if (MET > 200. ) WZdata_allnJetsMet200_ ->Fill(NJ);
						//nt0_em++;
				//	}
				}
				if (Flavor == 2) {       // E-E
				//	npp_ee += npp;
				//	npf_ee += (nfp+npf);
				//	nff_ee += nff;
				//	if (TLCat == 0)               nt2_ee++;
				// 	if (TLCat == 1 || TLCat == 2) nt10_ee++;
				// 	if (TLCat == 3) {
						WZdata_allMet_   ->Fill(MET);
						WZdata_allHT_    ->Fill(HT);
						WZdata_allnJets_ ->Fill(NJ);
						WZdata_bnJets_   ->Fill(NbJ);
						if (MET >  50. ) { WZdata_allnJetsMet50_ ->Fill(NJ); WZdata_bnJetsMet50_ ->Fill(NbJ); }
						if (MET > 100. ) WZdata_allnJetsMet100_ ->Fill(NJ);
						if (MET > 120. ) WZdata_allnJetsMet120_ ->Fill(NJ);
						if (MET > 200. ) WZdata_allnJetsMet200_ ->Fill(NJ);
						//nt0_ee++;
				//	}
				}
			}
		}
		// GET RARE MC EVENTS
		if (SType == 15 && TLCat == 0) { // tight-tight rare MC events
			if (!pythia) {
				if (*sname != "WZTo3LNuMadgraph") continue;
			}
			else {
				if (*sname != "WZTo3LNuPythia") continue;
			}
			if (Flavor > 2) continue;
			float scale = 1; //fLumiNorm / SLumi;
			float tmp_nt2_rare_mc    = puweight*trigScale[Flavor]*scale;
			// float tmp_nt2_rare_mc_e2 = tmp_nt2_rare_mc*tmp_nt2_rare_mc;
				if (!pythia) {
					WZmadgraph_allMet_   ->Fill(MET, tmp_nt2_rare_mc);
					WZmadgraph_allHT_    ->Fill(HT, tmp_nt2_rare_mc);
					WZmadgraph_allnJets_ ->Fill(NJ, tmp_nt2_rare_mc);
					WZmadgraph_bnJets_   ->Fill(NbJ, tmp_nt2_rare_mc);
					if (MET >  50. ) { WZmadgraph_allnJetsMet50_ ->Fill(NJ, tmp_nt2_rare_mc); WZmadgraph_bnJetsMet50_ ->Fill(NbJ, tmp_nt2_rare_mc); }
					if (MET > 100. ) WZmadgraph_allnJetsMet100_ ->Fill(NJ, tmp_nt2_rare_mc);
					if (MET > 120. ) WZmadgraph_allnJetsMet120_ ->Fill(NJ, tmp_nt2_rare_mc);
					if (MET > 200. ) WZmadgraph_allnJetsMet200_ ->Fill(NJ, tmp_nt2_rare_mc);
				}
				else {
					WZpythia_allMet_   ->Fill(MET, tmp_nt2_rare_mc);
					WZpythia_allHT_    ->Fill(HT, tmp_nt2_rare_mc);
					WZpythia_allnJets_ ->Fill(NJ, tmp_nt2_rare_mc);
					WZpythia_bnJets_   ->Fill(NbJ, tmp_nt2_rare_mc);
					if (MET >  50. ) { WZpythia_allnJetsMet50_ ->Fill(NJ, tmp_nt2_rare_mc); WZpythia_bnJetsMet50_ ->Fill(NbJ, tmp_nt2_rare_mc); }
					if (MET > 100. ) WZpythia_allnJetsMet100_ ->Fill(NJ, tmp_nt2_rare_mc);
					if (MET > 120. ) WZpythia_allnJetsMet120_ ->Fill(NJ, tmp_nt2_rare_mc);
					if (MET > 200. ) WZpythia_allnJetsMet200_ ->Fill(NJ, tmp_nt2_rare_mc);
				}
		} // end rare mc events
		
	}

	if (!pythia) {
		TFile * res_ = new TFile("WZmadgraph.root", "RECREATE", "res_");
		res_   -> cd();
		WZmadgraph_allnJets_ ->Write();
		WZmadgraph_bnJets_ ->Write();
		WZmadgraph_allMet_   ->Write();
		WZmadgraph_allHT_    ->Write();
		WZmadgraph_bnJetsMet50_  ->Write();
		WZmadgraph_allnJetsMet50_  ->Write();
		WZmadgraph_allnJetsMet100_ ->Write();
		WZmadgraph_allnJetsMet120_ ->Write();
		WZmadgraph_allnJetsMet200_ ->Write();
		WZdata_allnJets_ ->Write();
		WZdata_bnJets_ ->Write();
		WZdata_allMet_   ->Write();
		WZdata_allHT_    ->Write();
		WZdata_bnJetsMet50_  ->Write();
		WZdata_allnJetsMet50_  ->Write();
		WZdata_allnJetsMet100_ ->Write();
		WZdata_allnJetsMet120_ ->Write();
		WZdata_allnJetsMet200_ ->Write();
	}
	else {
		TFile * resPyth_ = new TFile("WZpythia.root", "RECREATE", "resPyth_");
		resPyth_   -> cd();
		WZpythia_allnJets_ ->Write();
		WZpythia_bnJets_ ->Write();
		WZpythia_allMet_   ->Write();
		WZpythia_allHT_    ->Write();
		WZpythia_bnJetsMet50_  ->Write();
		WZpythia_allnJetsMet50_  ->Write();
		WZpythia_allnJetsMet100_ ->Write();
		WZpythia_allnJetsMet120_ ->Write();
		WZpythia_allnJetsMet200_ ->Write();
	}
}
void SSDLPlotter::makeROCCurve(){

	cout << "making the ROC curves. might take a while..." << endl;
	fOutputSubDir = "ElectronIDPlots/";
	ofstream OUT(fOutputDir+fOutputSubDir+"ROCCurveInfo.txt", ios::trunc);

	TLatex *lat = new TLatex();
	lat->SetNDC(kTRUE);
	lat->SetTextColor(kBlack);
	lat->SetTextSize(0.04);

	TFile *pFile = TFile::Open(fOutputFileName);
	TTree *sigtree; getObjectSafe(pFile, "SigEvents", sigtree);
	string *sname = 0;
	int flag;
	int   SType, Flavor, TLCat, NJ, NbJ, NbJmed;
	float puweight, pT1, pT2, HT, MET, MT2, SLumi;
	float eta1, eta2, mll, mvaid1, mvaid2, medwp1, medwp2;
	float pfiso1, pfiso2;

	sigtree->SetBranchAddress("SystFlag", &flag);
	sigtree->SetBranchAddress("SName",    &sname);
	sigtree->SetBranchAddress("SType",    &SType);
	sigtree->SetBranchAddress("PUWeight", &puweight);
	sigtree->SetBranchAddress("SLumi",    &SLumi);
	sigtree->SetBranchAddress("Flavor",   &Flavor);
	sigtree->SetBranchAddress("pT1",      &pT1);
	sigtree->SetBranchAddress("pT2",      &pT2);
	sigtree->SetBranchAddress("eta1",     &eta1);
	sigtree->SetBranchAddress("eta2",     &eta2);
	sigtree->SetBranchAddress("TLCat",    &TLCat);
	sigtree->SetBranchAddress("HT",       &HT);
	sigtree->SetBranchAddress("MET",      &MET);
	sigtree->SetBranchAddress("MT2",      &MT2);
	sigtree->SetBranchAddress("NJ",       &NJ);
	sigtree->SetBranchAddress("NbJ",      &NbJ);
	sigtree->SetBranchAddress("NbJmed",   &NbJmed);
	sigtree->SetBranchAddress("Mll",      &mll);
	sigtree->SetBranchAddress("PFIso1",   &pfiso1);
	sigtree->SetBranchAddress("PFIso2",   &pfiso2);
	sigtree->SetBranchAddress("MVAID1",   &mvaid1);
	sigtree->SetBranchAddress("MVAID2",   &mvaid2);
	sigtree->SetBranchAddress("medWP1",   &medwp1);
	sigtree->SetBranchAddress("medWP2",   &medwp2);


	std::vector< std::pair < float , float> > etas;
	std::vector< std::pair < float , float> >::const_iterator eit;
	etas.push_back( make_pair(0.    , 0.8  ) );
	etas.push_back( make_pair(0.8   , 1.479) );
	etas.push_back( make_pair(1.479 , 2.5  ) );
	
	for (eit = etas.begin(); eit != etas.end(); eit++) {
		std::vector < std::pair < float, float >  > points_mvaid; // fill ttw eff in the first float, (1-ttj eff) in second float
	
		float ttw_medwp_eff, ttj_medwp_eff;
		for (float cut = 0.5; cut <= 1.; cut+=0.002) {

			ttw_medwp_eff = 0.; ttj_medwp_eff = 0.;
			int n_ttw_tot(0), n_ttw_pass(0);
			int n_ttw_medwp_pass(0);
			int n_ttj_tot(0), n_ttj_pass(0);
			int n_ttj_medwp_pass(0);

			for( int i = 0; i < sigtree->GetEntries(); i++ ){
				sigtree->GetEntry(i);
				
				if (*sname != "TTJets" && *sname != "TTbarW" ) continue; // just look at relevant samples
				if (Flavor != 2) continue; // only interested in events with two electrons
				// preselection
				if (pfiso1 > 0.05 || pfiso2 > 0.05 ) continue;
				if (NJ  < 3)      continue;

				// pt cuts
				gChannel chan = gChannel(Flavor);
				if(chan == ElMu){
					if(pT1 > pT2){
						if(pT1 < 20.) continue;
						if(pT2 < 20.) continue;
					}
					if(pT1 < pT2){
						if(pT1 < 20.) continue;
						if(pT2 < 20.) continue;
					}
				}
				else{
					if(pT1 < 20.) continue;
					if(pT2 < 20.) continue;
				}

				// get TTbarW events
				if (*sname == "TTbarW") {
					// first lepton
					if (fabs(eta1) > eit->first  && fabs(eta1) < eit->second) {
						n_ttw_tot ++;
						if (mvaid1 > cut) n_ttw_pass++;
						if (medwp1 == 1 ) n_ttw_medwp_pass++;
					}
					// second lepton
					if (fabs(eta2) > eit->first  && fabs(eta2) < eit->second) {
						n_ttw_tot ++;
						if (mvaid2 > cut) n_ttw_pass++;
						if (medwp2 == 1 ) n_ttw_medwp_pass++;
					}
				}
				// get TTJets events
				if (*sname == "TTJets" ) {
					// first lepton
					if (fabs(eta1) > eit->first  && fabs(eta1) < eit->second) {
						n_ttj_tot ++;
						if (mvaid1 > cut) n_ttj_pass++;
						if (medwp1 == 1 ) n_ttj_medwp_pass++;
					}
					// second lepton
					if (fabs(eta2) > eit->first  && fabs(eta2) < eit->second) {
						n_ttj_tot ++;
						if (mvaid2 > cut) n_ttj_pass++;
						if (medwp2 == 1 ) n_ttj_medwp_pass++;
					}
				}
			} // end loop on events in the tree

			// calculate the efficiencies
			ttw_medwp_eff = (float) n_ttw_medwp_pass/n_ttw_tot;
			ttj_medwp_eff = (float) n_ttj_medwp_pass/n_ttj_tot;
			float ttw_mvaid_eff = (float) n_ttw_pass/n_ttw_tot;
			float ttj_mvaid_eff = (float) n_ttj_pass/n_ttj_tot;

			points_mvaid.push_back( make_pair(ttw_mvaid_eff, 1. - ttj_mvaid_eff) );
			OUT << Form("Values: %4.3f \t %4.3f \t %4.3f --- sig eff: %4.3f || (1 - fake eff): %4.3f", cut, cut, cut, ttw_mvaid_eff, ttj_mvaid_eff) << endl;
			OUT << Form("         medium WP              --- sig eff: %4.3f || (1 - fake eff): %4.3f", ttw_medwp_eff, ttj_medwp_eff) << endl;
			// cout<< Form("Values: %4.3f \t %4.3f \t %4.3f --- sig eff: %4.3f || (1 - fake eff): %4.3f", cut, cut, cut, ttw_mvaid_eff, ttj_mvaid_eff) << endl;
			// cout<< Form("         medium WP              --- sig eff: %4.3f || (1 - fake eff): %4.3f", ttw_medwp_eff, ttj_medwp_eff) << endl;

		}

	// CREATING THE HISTGRAMS AND COSMETICS
	// ------------------------------------
	TH2D * h_medwp = new TH2D("medium WP ROC", "medium WP ROC", 50, 0., 1., 50, 0., 1.);
	TH2D * h_mvaid = new TH2D("MVA ID ROC"   , "MVA ID ROC"   , 50, 0., 1., 50, 0., 1.);
	
	// some cosmetics for the drawing
	h_medwp->SetMarkerStyle(21);
	h_medwp->SetMarkerColor(8);
	h_medwp->SetMarkerSize(1.1);

	h_mvaid->SetMarkerStyle(20);
	h_mvaid->SetMarkerColor(9);
	h_mvaid->SetMarkerSize(1.1);
	
	h_mvaid->GetXaxis()->SetTitle("signal eff.");
	h_mvaid->GetYaxis()->SetTitle("1 - ttj eff.");

	// FILLING THE HISTGRAMS
	// ---------------------
	// fill the medium WP histogram, it's just one point...
	h_medwp->Fill(ttw_medwp_eff, 1. - ttj_medwp_eff);
	// loop over all the mva ID points and fill them into the histogram
	std::vector< std::pair < float, float > >::const_iterator it;
	for (it = points_mvaid.begin(); it != points_mvaid.end(); it++) {
		h_mvaid->Fill(it->first, it->second);
	}

	TLegend *leg = new TLegend(0.15,0.62,0.50,0.88);
	leg->AddEntry(h_medwp,   "medium WP","p");
	leg->AddEntry(h_mvaid,   "MVA ID"   ,"p");
	leg->SetFillStyle(0);
	leg->SetTextFont(42);
	// leg->SetTextSize(0.05);
	leg->SetBorderSize(0);

	TCanvas *c_temp = new TCanvas("C_ObsPred", "Observed vs Predicted", 0, 0, 600, 600);
	c_temp->cd();
	
	h_mvaid->Draw("");
	h_medwp->Draw("same");
	leg->Draw();
	lat->DrawLatex(0.25,0.40, Form("%3.2f < |#eta| < %3.2f", eit->first, eit->second));
	
	drawTopLine(0.56);
	
	gPad->RedrawAxis();
	Util::PrintPDF(c_temp, Form("ROCcurve_eta%d", (int) eit->second), fOutputDir + fOutputSubDir);
	delete c_temp;	
	delete h_medwp, h_mvaid;

	}
	OUT.close();

}
void SSDLPlotter::plotWeightedHT(){
	fOutputSubDir = "sandbox/";
	//vector<int> samples;
	//samples.push_back(TTJets);

	TFile * file_ = new TFile("/scratch/mdunser/SSDLTrees/ttjets.root", "READ", "file_");
	TTree * tree = (TTree *) file_->Get("Analysis");
	tree->ResetBranchAddresses();

	TH1D *normalHT   = new TH1D("normalHT"   , "normalHT"   , 25 , 0. , 200.);
	TH1D *weightedHT = new TH1D("weightedHT" , "weightedHT" , 25 , 0. , 200.);
	
	normalHT->Sumw2();
	weightedHT->Sumw2();
	
//	for(size_t i = 0; i < samples.size(); ++i){
//		Sample *S = fSamples[samples[i]];
		
		//TTree *tree = S->getTree();
		//tree->ResetBranchAddresses();
		//Init(tree);
		Init(tree);
		for (Long64_t jentry=0; jentry<tree->GetEntriesFast();jentry++) {
			tree->GetEntry(jentry);
			printProgress(jentry, tree->GetEntriesFast(), "TTJets");
			normalHT->Fill(getHT());
			weightedHT->Fill(getWeightedHT());

		}
		//S->cleanUp();
	//}
	normalHT->GetYaxis()->SetRangeUser(0., 1.5*normalHT->GetMaximum());
	weightedHT->GetYaxis()->SetRangeUser(0., 1.5*normalHT->GetMaximum());
	TFile * res_ = new TFile(Form(fOutputDir+fOutputSubDir+"htCompare.root"), "RECREATE", "res_");
	res_   -> cd();
	normalHT->Write();
	weightedHT->Write();
	
	//printObject(normalHT  , "normalHT"  , "LEX");
	//printObject(weightedHT, "weightedHT", "LEX");
}

//____________________________________________________________________________
//SANTI  void SSDLPlotter::makePileUpPlots(bool write){
//SANTI  	fOutputSubDir = "PileUp/";
//SANTI  
//SANTI  	TH1D *smu_nvertices, *dmu_nvertices, *sel_nvertices, *del_nvertices, *mue_nvertices;
//SANTI  	TH1D *mu_ntight, *mu_nloose, *mu_ratio, *el_ntight, *el_nloose, *el_ratio;
//SANTI  
//SANTI  	if(!write){
//SANTI  		TFile *file = TFile::Open(fOutputDir + fOutputSubDir + "histos.root");
//SANTI  		smu_nvertices = (TH1D*)file->Get("smu_nvertices");
//SANTI  		dmu_nvertices = (TH1D*)file->Get("dmu_nvertices");
//SANTI  		sel_nvertices = (TH1D*)file->Get("sel_nvertices");
//SANTI  		del_nvertices = (TH1D*)file->Get("del_nvertices");
//SANTI  		mue_nvertices = (TH1D*)file->Get("mue_nvertices");
//SANTI  		mu_ratio      = (TH1D*)file->Get("mu_ratio");
//SANTI  		el_ratio      = (TH1D*)file->Get("el_ratio");
//SANTI  	}
//SANTI  	else{
//SANTI  		vector<gSample> samples;
//SANTI  		samples.push_back(DoubleMu1);
//SANTI  		// samples.push_back(DoubleMu2);
//SANTI  		// samples.push_back(DoubleMu3);
//SANTI  		// samples.push_back(DoubleMu4);
//SANTI  		// samples.push_back(DoubleMu5);
//SANTI  		samples.push_back(MuEG1);
//SANTI  		// samples.push_back(MuEG2);
//SANTI  		// samples.push_back(MuEG3);
//SANTI  		// samples.push_back(MuEG4);
//SANTI  		// samples.push_back(MuEG5);
//SANTI  		samples.push_back(DoubleEle1);
//SANTI  		// samples.push_back(DoubleEle2);
//SANTI  		// samples.push_back(DoubleEle3);
//SANTI  		// samples.push_back(DoubleEle4);
//SANTI  		// samples.push_back(DoubleEle5);
//SANTI  
//SANTI  		smu_nvertices = new TH1D("smu_nvertices", "smu_nvertices", FRatioPlots::nbins[3], FRatioPlots::xmin[3], FRatioPlots::xmax[3]);
//SANTI  		dmu_nvertices = new TH1D("dmu_nvertices", "dmu_nvertices", FRatioPlots::nbins[3], FRatioPlots::xmin[3], FRatioPlots::xmax[3]);
//SANTI  		sel_nvertices = new TH1D("sel_nvertices", "sel_nvertices", FRatioPlots::nbins[3], FRatioPlots::xmin[3], FRatioPlots::xmax[3]);
//SANTI  		del_nvertices = new TH1D("del_nvertices", "del_nvertices", FRatioPlots::nbins[3], FRatioPlots::xmin[3], FRatioPlots::xmax[3]);
//SANTI  		mue_nvertices = new TH1D("mue_nvertices", "mue_nvertices", FRatioPlots::nbins[3], FRatioPlots::xmin[3], FRatioPlots::xmax[3]);
//SANTI  
//SANTI  		mu_ntight = new TH1D("mu_ntight", "ntight", FRatioPlots::nbins[3], FRatioPlots::xmin[3], FRatioPlots::xmax[3]);
//SANTI  		mu_nloose = new TH1D("mu_nloose", "nloose", FRatioPlots::nbins[3], FRatioPlots::xmin[3], FRatioPlots::xmax[3]);
//SANTI  		mu_ratio  = new TH1D("mu_ratio",  "ratio",  FRatioPlots::nbins[3], FRatioPlots::xmin[3], FRatioPlots::xmax[3]);
//SANTI  		el_ntight = new TH1D("el_ntight", "ntight", FRatioPlots::nbins[3], FRatioPlots::xmin[3], FRatioPlots::xmax[3]);
//SANTI  		el_nloose = new TH1D("el_nloose", "nloose", FRatioPlots::nbins[3], FRatioPlots::xmin[3], FRatioPlots::xmax[3]);
//SANTI  		el_ratio  = new TH1D("el_ratio",  "ratio",  FRatioPlots::nbins[3], FRatioPlots::xmin[3], FRatioPlots::xmax[3]);
//SANTI  
//SANTI  		mu_ntight    ->Sumw2();
//SANTI  		mu_nloose    ->Sumw2();
//SANTI  		mu_ratio     ->Sumw2();
//SANTI  		el_ntight    ->Sumw2();
//SANTI  		el_nloose    ->Sumw2();
//SANTI  		el_ratio     ->Sumw2();
//SANTI  		smu_nvertices->Sumw2();
//SANTI  		dmu_nvertices->Sumw2();
//SANTI  		sel_nvertices->Sumw2();
//SANTI  		del_nvertices->Sumw2();
//SANTI  		mue_nvertices->Sumw2();
//SANTI  		mu_ntight    ->SetXTitle("N_{Vertices}");
//SANTI  		mu_nloose    ->SetXTitle("N_{Vertices}");
//SANTI  		mu_ratio     ->SetXTitle("N_{Vertices}");
//SANTI  		el_ntight    ->SetXTitle("N_{Vertices}");
//SANTI  		el_nloose    ->SetXTitle("N_{Vertices}");
//SANTI  		el_ratio     ->SetXTitle("N_{Vertices}");
//SANTI  		smu_nvertices->SetXTitle("N_{Vertices}");
//SANTI  		dmu_nvertices->SetXTitle("N_{Vertices}");
//SANTI  		sel_nvertices->SetXTitle("N_{Vertices}");
//SANTI  		del_nvertices->SetXTitle("N_{Vertices}");
//SANTI  		mue_nvertices->SetXTitle("N_{Vertices}");
//SANTI  		mu_ntight    ->SetYTitle("N_{Events}");
//SANTI  		mu_nloose    ->SetYTitle("N_{Events}");
//SANTI  		mu_ratio     ->SetYTitle("N_{Events}");
//SANTI  		el_ntight    ->SetYTitle("N_{Events}");
//SANTI  		el_nloose    ->SetYTitle("N_{Events}");
//SANTI  		el_ratio     ->SetYTitle("N_{Events}");
//SANTI  		smu_nvertices->SetYTitle("N_{Events}");
//SANTI  		dmu_nvertices->SetYTitle("N_{Events}");
//SANTI  		sel_nvertices->SetYTitle("N_{Events}");
//SANTI  		del_nvertices->SetYTitle("N_{Events}");
//SANTI  		mue_nvertices->SetYTitle("N_{Events}");
//SANTI  		mu_ntight    ->GetYaxis()->SetTitleOffset(1.2);
//SANTI  		mu_nloose    ->GetYaxis()->SetTitleOffset(1.2);
//SANTI  		mu_ratio     ->GetYaxis()->SetTitleOffset(1.2);
//SANTI  		el_ntight    ->GetYaxis()->SetTitleOffset(1.2);
//SANTI  		el_nloose    ->GetYaxis()->SetTitleOffset(1.2);
//SANTI  		el_ratio     ->GetYaxis()->SetTitleOffset(1.2);
//SANTI  		smu_nvertices->GetYaxis()->SetTitleOffset(1.2);
//SANTI  		dmu_nvertices->GetYaxis()->SetTitleOffset(1.2);
//SANTI  		sel_nvertices->GetYaxis()->SetTitleOffset(1.2);
//SANTI  		del_nvertices->GetYaxis()->SetTitleOffset(1.2);
//SANTI  		mue_nvertices->GetYaxis()->SetTitleOffset(1.2);
//SANTI  
//SANTI  		for(size_t i = 0; i < samples.size(); ++i){
//SANTI  			Sample *S = fSamples[samples[i]];
//SANTI  			fSample = S; // necessary for triggers to work properly
//SANTI  			fCurrentSample = samples[i];
//SANTI  			resetHypLeptons();
//SANTI  			fDoCounting = false;
//SANTI  		
//SANTI  			TTree *tree = S->getTree();
//SANTI  			tree->ResetBranchAddresses();
//SANTI  			Init(tree);
//SANTI  			for (Long64_t jentry=0; jentry<tree->GetEntriesFast();jentry++) {
//SANTI  				tree->GetEntry(jentry);
//SANTI  				printProgress(jentry, tree->GetEntriesFast(), S->sname);
//SANTI  
//SANTI  				if(mumuSignalTrigger()) dmu_nvertices->Fill(NVrtx);
//SANTI  				if(elelSignalTrigger()) del_nvertices->Fill(NVrtx);
//SANTI  				if(elmuSignalTrigger()) mue_nvertices->Fill(NVrtx);
//SANTI  				if(singleMuTrigger()){
//SANTI  					smu_nvertices->Fill(NVrtx);
//SANTI  					if(isSigSupMuEvent()){
//SANTI  						if(isTightMuon(0)) mu_ntight->Fill(NVrtx);
//SANTI  						if(isLooseMuon(0)) mu_nloose->Fill(NVrtx);					
//SANTI  					}
//SANTI  				}
//SANTI  				resetHypLeptons();
//SANTI  				if(singleElTrigger()){
//SANTI  					sel_nvertices->Fill(NVrtx);
//SANTI  				 	if(isSigSupElEvent()){
//SANTI  						if(isTightElectron(0)) el_ntight->Fill(NVrtx);
//SANTI  						if(isLooseElectron(0)) el_nloose->Fill(NVrtx);
//SANTI  					}
//SANTI  				}
//SANTI  			}
//SANTI  			S->cleanUp();
//SANTI  		}
//SANTI  	
//SANTI  		mu_ratio->Divide(mu_ntight, mu_nloose, 1., 1., "B");
//SANTI  		mu_ratio->GetYaxis()->SetRangeUser(0., 0.4);
//SANTI  		el_ratio->Divide(el_ntight, el_nloose, 1., 1., "B");
//SANTI  		el_ratio->GetYaxis()->SetRangeUser(0., 0.4);
//SANTI  
//SANTI  		dmu_nvertices->Scale(1./dmu_nvertices->Integral());
//SANTI  		smu_nvertices->Scale(1./smu_nvertices->Integral());
//SANTI  		del_nvertices->Scale(1./del_nvertices->Integral());
//SANTI  		sel_nvertices->Scale(1./sel_nvertices->Integral());
//SANTI  		mue_nvertices->Scale(1./mue_nvertices->Integral());
//SANTI  
//SANTI  	}
//SANTI  
//SANTI  	mu_ratio     ->SetYTitle("N_{Events} (Normalized)");
//SANTI  	el_ratio     ->SetYTitle("N_{Events} (Normalized)");
//SANTI  	smu_nvertices->SetYTitle("N_{Events} (Normalized)");
//SANTI  	dmu_nvertices->SetYTitle("N_{Events} (Normalized)");
//SANTI  	sel_nvertices->SetYTitle("N_{Events} (Normalized)");
//SANTI  	del_nvertices->SetYTitle("N_{Events} (Normalized)");
//SANTI  	mue_nvertices->SetYTitle("N_{Events} (Normalized)");
//SANTI  
//SANTI  	mu_ratio     ->GetYaxis()->SetTitleOffset(1.3);
//SANTI  	el_ratio     ->GetYaxis()->SetTitleOffset(1.3);
//SANTI  	smu_nvertices->GetYaxis()->SetTitleOffset(1.3);
//SANTI  	dmu_nvertices->GetYaxis()->SetTitleOffset(1.3);
//SANTI  	sel_nvertices->GetYaxis()->SetTitleOffset(1.3);
//SANTI  	del_nvertices->GetYaxis()->SetTitleOffset(1.3);
//SANTI  	mue_nvertices->GetYaxis()->SetTitleOffset(1.3);
//SANTI  
//SANTI  	// Color_t colors[6] = {1, 1, 1, 1, 1, 1};
//SANTI  	Color_t colors[5] = {31, 41, 51, 61, 81};
//SANTI  	// Color_t colors[6] = {1, 12, 39, 38, 32, 30};//, 29};
//SANTI  	Style_t styles[6] = {24, 25, 26, 27, 32, 30};//, 29};
//SANTI  
//SANTI  	// dmu_nvertices->SetMarkerStyle(styles[0]);
//SANTI  	// smu_nvertices->SetMarkerStyle(styles[1]);
//SANTI  	// del_nvertices->SetMarkerStyle(styles[2]);
//SANTI  	// sel_nvertices->SetMarkerStyle(styles[3]);
//SANTI  	// mue_nvertices->SetMarkerStyle(styles[4]);
//SANTI  	// dmu_nvertices->SetMarkerColor(colors[0]);
//SANTI  	// smu_nvertices->SetMarkerColor(colors[1]);
//SANTI  	// del_nvertices->SetMarkerColor(colors[2]);
//SANTI  	// sel_nvertices->SetMarkerColor(colors[3]);
//SANTI  	// mue_nvertices->SetMarkerColor(colors[4]);
//SANTI  	dmu_nvertices->SetLineColor(colors[0]);
//SANTI  	smu_nvertices->SetLineColor(colors[1]);
//SANTI  	del_nvertices->SetLineColor(colors[2]);
//SANTI  	sel_nvertices->SetLineColor(colors[3]);
//SANTI  	mue_nvertices->SetLineColor(colors[4]);
//SANTI  	// dmu_nvertices->SetMarkerSize(1.5);
//SANTI  	// smu_nvertices->SetMarkerSize(1.5);
//SANTI  	// del_nvertices->SetMarkerSize(1.5);
//SANTI  	// sel_nvertices->SetMarkerSize(1.5);
//SANTI  	// mue_nvertices->SetMarkerSize(1.5);
//SANTI  	dmu_nvertices->SetLineWidth(2);
//SANTI  	smu_nvertices->SetLineWidth(2);
//SANTI  	del_nvertices->SetLineWidth(2);
//SANTI  	sel_nvertices->SetLineWidth(2);
//SANTI  	mue_nvertices->SetLineWidth(2);
//SANTI  	dmu_nvertices->SetFillStyle(0);
//SANTI  	smu_nvertices->SetFillStyle(0);
//SANTI  	del_nvertices->SetFillStyle(0);
//SANTI  	sel_nvertices->SetFillStyle(0);
//SANTI  	mue_nvertices->SetFillStyle(0);
//SANTI  
//SANTI  	dmu_nvertices->GetYaxis()->SetRangeUser(0., 0.4);
//SANTI  	smu_nvertices->GetYaxis()->SetRangeUser(0., 0.4);
//SANTI  	del_nvertices->GetYaxis()->SetRangeUser(0., 0.4);
//SANTI  	sel_nvertices->GetYaxis()->SetRangeUser(0., 0.4);
//SANTI  	mue_nvertices->GetYaxis()->SetRangeUser(0., 0.4);
//SANTI  
//SANTI  	mu_ratio->SetMarkerStyle(20);
//SANTI  	mu_ratio->SetMarkerColor(kBlue);
//SANTI  	mu_ratio->SetLineColor(kBlue);
//SANTI  	mu_ratio->SetMarkerSize(1.8);
//SANTI  	mu_ratio->SetLineWidth(2);
//SANTI  	el_ratio->SetMarkerStyle(21);
//SANTI  	el_ratio->SetMarkerColor(kRed);
//SANTI  	el_ratio->SetLineColor(kRed);
//SANTI  	el_ratio->SetMarkerSize(1.8);
//SANTI  	el_ratio->SetLineWidth(2);
//SANTI  
//SANTI  	// printObject(mu_ratio, "MuRatio", "PEX");
//SANTI  	// printObject(el_ratio, "ElRatio", "PEX");
//SANTI  	// printObject(dmu_nvertices, "NVertices_DMu", "L");
//SANTI  	// printObject(smu_nvertices, "NVertices_SMu", "L");
//SANTI  	// printObject(del_nvertices, "NVertices_DEl", "L");
//SANTI  	// printObject(sel_nvertices, "NVertices_SEl", "L");
//SANTI  	// printObject(mue_nvertices, "NVertices_MuE", "L");
//SANTI  
//SANTI  	// TLegend *leg = new TLegend(0.35,0.15,0.70,0.40);
//SANTI  	TLegend *leg = new TLegend(0.50,0.60,0.88,0.88);
//SANTI  	leg->AddEntry(dmu_nvertices, Form("DoubleMu Trig., Mean = %4.2f", dmu_nvertices->GetMean()),  "l");
//SANTI  	leg->AddEntry(smu_nvertices, Form("SingleMu Trig., Mean = %4.2f", smu_nvertices->GetMean()),  "l");
//SANTI  	leg->AddEntry(del_nvertices, Form("DoubleEle Trig., Mean = %4.2f", del_nvertices->GetMean()), "l");
//SANTI  	leg->AddEntry(sel_nvertices, Form("SingleEle Trig., Mean = %4.2f", sel_nvertices->GetMean()), "l");
//SANTI  	leg->AddEntry(mue_nvertices, Form("MuEle Trig., Mean = %4.2f", mue_nvertices->GetMean()),   "l");
//SANTI  	leg->AddEntry(mu_ratio, Form("TL Ratio (Muons)"), "p");
//SANTI  	leg->AddEntry(el_ratio, Form("TL Ratio (Electrons)"),  "p");
//SANTI  	leg->SetFillStyle(0);
//SANTI  	leg->SetTextFont(42);
//SANTI  	// leg->SetTextSize(0.03);
//SANTI  	leg->SetBorderSize(0);
//SANTI  
//SANTI  	TCanvas *c_temp = new TCanvas("C_temp", "HT vs MET in Data vs MC", 0, 0, 800, 600);
//SANTI  	c_temp->cd();
//SANTI  	dmu_nvertices->Draw("axis");
//SANTI  	dmu_nvertices->DrawCopy("hist L same");
//SANTI  	smu_nvertices->DrawCopy("hist L same");
//SANTI  	del_nvertices->DrawCopy("hist L same");
//SANTI  	sel_nvertices->DrawCopy("hist L same");
//SANTI  	mue_nvertices->DrawCopy("hist L same");
//SANTI  	mu_ratio->DrawCopy("PE X0 same");
//SANTI  	el_ratio->DrawCopy("PE X0 same");
//SANTI  	TGaxis *axis = new TGaxis(18, 0, 18, 0.4, 0, 0.4, 510, "+L");
//SANTI  	axis->SetLabelFont(42);
//SANTI  	axis->SetTitleFont(42);
//SANTI  	axis->SetTitleOffset(1.2);
//SANTI  	axis->SetTitle("TL Ratio");
//SANTI  	axis->Draw();
//SANTI  	leg->Draw();
//SANTI  	drawTopLine();
//SANTI  	// fLatex->DrawLatex(0.10,0.92, fSamples[sample]->name);
//SANTI  	// Util::PrintNoEPS( c_temp, "PileUp", fOutputDir + fOutputSubDir, NULL);
//SANTI  	Util::PrintPDF(   c_temp, "PileUp", fOutputDir + fOutputSubDir);
//SANTI  
//SANTI  	if(write){
//SANTI  		TFile *file = new TFile(fOutputDir + fOutputSubDir + "histos.root", "RECREATE");
//SANTI  		mu_ratio->Write(mu_ratio->GetName(), TObject::kWriteDelete);
//SANTI  		el_ratio->Write(el_ratio->GetName(), TObject::kWriteDelete);
//SANTI  		dmu_nvertices->Write(dmu_nvertices->GetName(), TObject::kWriteDelete);
//SANTI  		smu_nvertices->Write(smu_nvertices->GetName(), TObject::kWriteDelete);
//SANTI  		del_nvertices->Write(del_nvertices->GetName(), TObject::kWriteDelete);
//SANTI  		sel_nvertices->Write(sel_nvertices->GetName(), TObject::kWriteDelete);
//SANTI  		mue_nvertices->Write(mue_nvertices->GetName(), TObject::kWriteDelete);
//SANTI  		file->Close();
//SANTI  	}
//SANTI  
//SANTI  	// Cleanup
//SANTI  	delete leg, c_temp;
//SANTI  	if(write) delete mu_ntight, mu_nloose, el_ntight, el_nloose, mu_ratio, el_ratio;
//SANTI  	else delete mu_ratio, el_ratio;
//SANTI  	delete dmu_nvertices, smu_nvertices, del_nvertices, sel_nvertices, mue_nvertices;
//SANTI  	fOutputSubDir = "";
//SANTI  }
void SSDLPlotter::makePileUpPlots(bool RunPeriod){
	fOutputSubDir = "PileUp/";
	
	if (RunPeriod) return;
	
	TH1D *hel_dtrig, *hel_strig, *hel_ssdl, *hel_ntight, *hel_nloose; 
	TH1D *hmu_dtrig, *hmu_strig, *hmu_ssdl, *hmu_ntight, *hmu_nloose; 
	TH1D *h_ratio;
	
	hel_dtrig  = new TH1D("hel_dtrig" , "DoubleTrigger_Elec", 50, 0.5, 49.5);
	hel_strig  = new TH1D("hel_strig" , "SingleTrigger_Elec", 50, 0.5, 49.5);
	hel_ssdl   = new TH1D("hel_ssdl"  , "SSLLEvent_Elec"    , 50, 0.5, 49.5);
	hel_ntight = new TH1D("hel_ntight", "SigSupTight_Elec"  , 50, 0.5, 49.5);
	hel_nloose = new TH1D("hel_nloose", "SigSupLoose_Elec"  , 50, 0.5, 49.5);

	hmu_dtrig  = new TH1D("hmu_dtrig" , "DoubleTrigger_Muon", 50, 0.5, 49.5);
	hmu_strig  = new TH1D("hmu_strig" , "SingleTrigger_Muon", 50, 0.5, 49.5);
	hmu_ssdl   = new TH1D("hmu_ssdl"  , "SSLLEvent_Muon"    , 50, 0.5, 49.5);
	hmu_ntight = new TH1D("hmu_ntight", "SigSupTight_Muon"  , 50, 0.5, 49.5);
	hmu_nloose = new TH1D("hmu_nloose", "SigSupLoose_Muon"  , 50, 0.5, 49.5);
	
	h_ratio    = new TH1D("h_ratio"   , "Ratio"             , 50, 0.5, 49.5);
	
	hel_dtrig ->Sumw2();    
	hel_strig ->Sumw2();	
	hel_ssdl  ->Sumw2();	
	hel_ntight->Sumw2();	
	hel_nloose->Sumw2();	
				
	hmu_dtrig ->Sumw2();	
	hmu_strig ->Sumw2();	
	hmu_ssdl  ->Sumw2();	
	hmu_ntight->Sumw2();	
	hmu_nloose->Sumw2();	
	
	h_ratio   ->Sumw2();
	
	hel_dtrig ->SetYTitle("N_{Events}");    hel_dtrig ->GetXaxis()->SetTitleOffset(1.2);
	hel_strig ->SetYTitle("N_{Events}");	hel_strig ->GetXaxis()->SetTitleOffset(1.2);
	hel_ssdl  ->SetYTitle("N_{Events}");	hel_ssdl  ->GetXaxis()->SetTitleOffset(1.2);
	hel_ntight->SetYTitle("N_{Events}");	hel_ntight->GetXaxis()->SetTitleOffset(1.2);
	hel_nloose->SetYTitle("N_{Events}");	hel_nloose->GetXaxis()->SetTitleOffset(1.2);

	hmu_dtrig ->SetYTitle("N_{Events}");	hmu_dtrig ->GetXaxis()->SetTitleOffset(1.2);
	hmu_strig ->SetYTitle("N_{Events}");	hmu_strig ->GetXaxis()->SetTitleOffset(1.2);
	hmu_ssdl  ->SetYTitle("N_{Events}");	hmu_ssdl  ->GetXaxis()->SetTitleOffset(1.2);
	hmu_ntight->SetYTitle("N_{Events}");	hmu_ntight->GetXaxis()->SetTitleOffset(1.2);
	hmu_nloose->SetYTitle("N_{Events}");	hmu_nloose->GetXaxis()->SetTitleOffset(1.2);
	
	for(size_t j = 1; j <= hmu_dtrig ->GetNbinsX(); ++j)   hmu_dtrig ->GetXaxis()->SetBinLabel(j, "");
	for(size_t j = 1; j <= hmu_strig ->GetNbinsX(); ++j)   hmu_strig ->GetXaxis()->SetBinLabel(j, "");
	for(size_t j = 1; j <= hmu_ssdl  ->GetNbinsX(); ++j)   hmu_ssdl  ->GetXaxis()->SetBinLabel(j, "");
	for(size_t j = 1; j <= hmu_ntight->GetNbinsX(); ++j)   hmu_ntight->GetXaxis()->SetBinLabel(j, "");
	for(size_t j = 1; j <= hmu_nloose->GetNbinsX(); ++j)   hmu_nloose->GetXaxis()->SetBinLabel(j, "");

	for(size_t j = 1; j <= hel_dtrig ->GetNbinsX(); ++j)   hel_dtrig ->GetXaxis()->SetBinLabel(j, "");
	for(size_t j = 1; j <= hel_strig ->GetNbinsX(); ++j)   hel_strig ->GetXaxis()->SetBinLabel(j, "");
	for(size_t j = 1; j <= hel_ssdl  ->GetNbinsX(); ++j)   hel_ssdl  ->GetXaxis()->SetBinLabel(j, "");
	for(size_t j = 1; j <= hel_ntight->GetNbinsX(); ++j)   hel_ntight->GetXaxis()->SetBinLabel(j, "");
	for(size_t j = 1; j <= hel_nloose->GetNbinsX(); ++j)   hel_nloose->GetXaxis()->SetBinLabel(j, "");

	// Fill Plots
	vector<int> musamples = fMuData;
	vector<int> elsamples = fEGData;
	for (size_t j=0; j<musamples.size(); j++){
		Sample *S = fSamples[j];
		hmu_dtrig ->Add(S->puplots[0].hdtrig );
		hmu_strig ->Add(S->puplots[0].hstrig );
		hmu_ssdl  ->Add(S->puplots[0].hssdl  );
		hmu_ntight->Add(S->puplots[0].hntight);
		hmu_nloose->Add(S->puplots[0].hnloose);
	}

	for (size_t j=0; j<elsamples.size(); j++){
		Sample *S = fSamples[j];
		hel_dtrig ->Add(S->puplots[1].hdtrig );
		hel_strig ->Add(S->puplots[1].hstrig );
		hel_ssdl  ->Add(S->puplots[1].hssdl  );
		hel_ntight->Add(S->puplots[1].hntight);
		hel_nloose->Add(S->puplots[1].hnloose);
	}
	
	hmu_dtrig ->Scale(1./std::max(1.0,(double)hmu_dtrig ->Integral()));
	hmu_strig ->Scale(1./std::max(1.0,(double)hmu_strig ->Integral()));
	hmu_ssdl  ->Scale(1./std::max(1.0,(double)hmu_ssdl  ->Integral()));
	hmu_ntight->Scale(1./std::max(1.0,(double)hmu_ntight->Integral()));
	hmu_nloose->Scale(1./std::max(1.0,(double)hmu_nloose->Integral()));

	hel_dtrig ->Scale(1./std::max(1.0,(double)hel_dtrig ->Integral()));
	hel_strig ->Scale(1./std::max(1.0,(double)hel_strig ->Integral()));
	hel_ssdl  ->Scale(1./std::max(1.0,(double)hel_ssdl  ->Integral()));
	hel_ntight->Scale(1./std::max(1.0,(double)hel_ntight->Integral()));
	hel_nloose->Scale(1./std::max(1.0,(double)hel_nloose->Integral()));

	// Now make the plots prettier
	hel_dtrig ->SetMarkerColor(kRed+1 );  hel_dtrig ->SetLineColor(kRed+1 );
	hel_strig ->SetMarkerColor(kBlue-4);  hel_strig ->SetLineColor(kBlue-4);
	hel_ssdl  ->SetMarkerColor(kRed+1 );  hel_ssdl  ->SetLineColor(kRed+1 );
	hel_nloose->SetMarkerColor(kBlue-4);  hel_nloose->SetLineColor(kBlue-4);

	hmu_dtrig ->SetMarkerColor(kRed+1 );  hmu_dtrig ->SetLineColor(kRed+1 );
	hmu_strig ->SetMarkerColor(kBlue-4);  hmu_strig ->SetLineColor(kBlue-4);
	hmu_ssdl  ->SetMarkerColor(kRed+1 );  hmu_ssdl  ->SetLineColor(kRed+1 );
	hmu_nloose->SetMarkerColor(kBlue-4);  hmu_nloose->SetLineColor(kBlue-4);

	// Plot everything
	float border = 0.3;
	float scale = (1-border)/border;
	TString title = "N_{Vertices}";
	h_ratio->SetXTitle(title);
	h_ratio->SetYTitle("");
	h_ratio->GetXaxis()->SetTitleSize(scale * 0.04);
	h_ratio->GetXaxis()->SetLabelSize(scale * hmu_dtrig->GetXaxis()->GetLabelSize());
	h_ratio->GetYaxis()->SetLabelSize(scale * hmu_dtrig->GetYaxis()->GetLabelSize());
	h_ratio->GetXaxis()->SetTickLength(scale * hmu_dtrig->GetXaxis()->GetTickLength());
	h_ratio->GetYaxis()->SetTickLength(hmu_dtrig->GetYaxis()->GetTickLength());
	
	h_ratio->SetFillStyle(1001);
	h_ratio->SetLineWidth(1);
	h_ratio->SetFillColor(  kGray+1);
	h_ratio->SetLineColor(  kGray+1);
	h_ratio->SetMarkerColor(kGray+1);

	TCanvas *c_temp = new TCanvas("C_temp", "Pile-up", 0, 0, 800, 600);
	c_temp->cd();
	drawTopLine();

	TPad *p_plot  = new TPad("plotpad",  "Pad containing the plot", 0.00, border, 1.00, 1.00, 0, 0);
	p_plot->SetBottomMargin(0.015);
	p_plot->Draw();
	TPad *p_ratio = new TPad("ratiopad", "Pad containing the ratio", 0.00, 0.00, 1.00, border, 0, 0);
	p_ratio->SetTopMargin(0.025);
	p_ratio->SetBottomMargin(0.35);
	p_ratio->Draw();
	
	p_ratio->cd();
	h_ratio->Divide(hel_dtrig,hel_strig);
	h_ratio->SetMaximum(1.99);
	h_ratio->SetMinimum(0.0);
	h_ratio->DrawCopy("P E2");
	TLine *l3 = new TLine(0.5, 1.00, 49.5, 1.00);
	l3->SetLineWidth(2);
	l3->SetLineStyle(7);
	l3->Draw();
	gPad->RedrawAxis();

	// electrons
	TLegend *leg1 = new TLegend(0.50,0.60,0.88,0.88);
	leg1->AddEntry(hel_dtrig, Form("DoubleEl Trig., Mean = %4.2f", hel_dtrig->GetMean()), "l");
	leg1->AddEntry(hel_strig, Form("SingleEl Trig., Mean = %4.2f", hel_strig->GetMean()), "l");
	leg1->SetFillStyle(0);
	leg1->SetTextFont(42);
	leg1->SetBorderSize(0);
	
	p_plot->cd();
	hel_dtrig->DrawCopy("P");
	hel_strig->DrawCopy("P SAME");
	leg1->Draw();
	gPad->RedrawAxis();
	Util::PrintPDF(   c_temp, "NVrtx_Trigger_Elec", fOutputDir + fOutputSubDir);
	
	//  muons
	TLegend *leg2 = new TLegend(0.50,0.60,0.88,0.88);
	leg2->AddEntry(hmu_dtrig, Form("DoubleMu Trig., Mean = %4.2f", hmu_dtrig->GetMean()), "l");
	leg2->AddEntry(hmu_strig, Form("SingleMu Trig., Mean = %4.2f", hmu_strig->GetMean()), "l");
	leg2->SetFillStyle(0);
	leg2->SetTextFont(42);
	leg2->SetBorderSize(0);
	
	p_plot->cd();
	hmu_dtrig->DrawCopy("P");
	hmu_strig->DrawCopy("P SAME");
	leg2->Draw();
	
	p_ratio->cd();
	h_ratio->Divide(hmu_dtrig,hmu_strig);
	h_ratio->SetMaximum(1.99);
	h_ratio->SetMinimum(0.0);
	h_ratio->DrawCopy("E2 ");
	l3->Draw();
	gPad->RedrawAxis();
	Util::PrintPDF(   c_temp, "NVrtx_Trigger_Muon", fOutputDir + fOutputSubDir);


	// electrons
	TLegend *leg3 = new TLegend(0.50,0.60,0.88,0.88);
	leg3->AddEntry(hel_ssdl  , Form("SSDL El Event.  , Mean = %4.2f", hel_ssdl->GetMean()), "l");
	leg3->AddEntry(hel_nloose, Form("SigSup El Event., Mean = %4.2f", hel_nloose->GetMean()), "l");
	leg3->SetFillStyle(0);
	leg3->SetTextFont(42);
	leg3->SetBorderSize(0);
	
	p_plot->cd();
	hel_ssdl->DrawCopy("P");
	hel_nloose->DrawCopy("P SAME");
	leg3->Draw();
	gPad->RedrawAxis();

	p_ratio->cd();
	h_ratio->Divide(hel_ssdl,hel_nloose);
	h_ratio->SetMaximum(1.99);
	h_ratio->SetMinimum(0.0);
	h_ratio->DrawCopy("E2 ");
	l3->Draw();
	gPad->RedrawAxis();
	Util::PrintPDF(   c_temp, "NVrtx_SSDLvsSigSup_Elec", fOutputDir + fOutputSubDir);
	
	//  muons
	TLegend *leg4 = new TLegend(0.50,0.60,0.88,0.88);
	leg4->AddEntry(hmu_ssdl  , Form("SSDL Mu Event.  , Mean = %4.2f", hmu_ssdl->GetMean()), "l");
	leg4->AddEntry(hmu_nloose, Form("SigSup Mu Event., Mean = %4.2f", hmu_nloose->GetMean()), "l");
	leg4->SetFillStyle(0);
	leg4->SetTextFont(42);
	leg4->SetBorderSize(0);
	
	p_plot->cd();
	hmu_ssdl->DrawCopy("P");
	hmu_nloose->DrawCopy("P SAME");
	leg4->Draw();
	
	p_ratio->cd();
	h_ratio->Divide(hmu_ssdl,hmu_nloose);
	h_ratio->SetMaximum(1.99);
	h_ratio->SetMinimum(0.0);
	h_ratio->DrawCopy("E2 ");
	l3->Draw();
	gPad->RedrawAxis();
	Util::PrintPDF(   c_temp, "NVrtx_SSDLvsSigSup_Muon", fOutputDir + fOutputSubDir);
     


	// Cleanup
	delete leg1,leg2;
	delete p_plot, p_ratio;
	delete c_temp;
}

//____________________________________________________________________________
void SSDLPlotter::makeNT012Plots(vector<int> mcsamples, gChannel chan, int reg, gHiLoSwitch hilo){
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
		float scale = fLumiNorm / S->getLumi();
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
			hnt20[i] = new TH1D(Form("nt20_%s", fSamples[index]->sname.Data()), "Observed Nt20", getNFPtBins(Muon), getFPtBins(Muon));
			hnt10[i] = new TH1D(Form("nt10_%s", fSamples[index]->sname.Data()), "Observed Nt10", getNFPtBins(Muon), getFPtBins(Muon));
			hnt01[i] = new TH1D(Form("nt01_%s", fSamples[index]->sname.Data()), "Observed Nt01", getNFPtBins(Muon), getFPtBins(Muon));
			hnt00[i] = new TH1D(Form("nt00_%s", fSamples[index]->sname.Data()), "Observed Nt00", getNFPtBins(Muon), getFPtBins(Muon));
			hnt20[i]->SetFillColor(fSamples[index]->color);
			hnt10[i]->SetFillColor(fSamples[index]->color);
			hnt01[i]->SetFillColor(fSamples[index]->color);
			hnt00[i]->SetFillColor(fSamples[index]->color);
			hnt20[i]->Sumw2();
			hnt10[i]->Sumw2();
			hnt01[i]->Sumw2();
			hnt00[i]->Sumw2();
			float scale = fLumiNorm / fSamples[index]->getLumi();
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
// void SSDLPlotter::makeRelIsoTTSigPlots(){
// 	char cmd[100];
//     sprintf(cmd,"mkdir -p %s%s", fOutputDir.Data(), fOutputSubDir.Data());
//     system(cmd);

// 	TLatex *lat = new TLatex();
// 	lat->SetNDC(kTRUE);
// 	lat->SetTextColor(kBlack);
// 	lat->SetTextSize(0.04);

// 	TLatex *el_lat = new TLatex();
// 	el_lat->SetNDC(kTRUE);
// 	el_lat->SetTextColor(kBlack);
// 	el_lat->SetTextSize(0.04);

// 	// Create histograms
// 	TH1D * hiso_ttbar  = new TH1D("MuIsoTTbar_"   , "Muon Isolation in TTbar for " , IsoPlots::nbins[0], 0., 1.);
// 	TH1D * hiso_signal = new TH1D("MuIsoSignal_"  , "Muon Isolation in Signal for ", IsoPlots::nbins[0], 0., 1.);
// 	hiso_ttbar ->Sumw2();
// 	hiso_signal->Sumw2();
// 	TH1D * el_hiso_ttbar  = new TH1D("ElIsoTTbar_"   , "Electron Isolation in TTbar for " , IsoPlots::nbins[0], 0., 0.6);
// 	TH1D * el_hiso_signal = new TH1D("ElIsoSignal_"  , "Electron Isolation in Signal for ", IsoPlots::nbins[0], 0., 0.6);
// 	el_hiso_ttbar ->Sumw2();
//  	el_hiso_signal->Sumw2();

// 	////////////////////////////////////////////////////
// 	// Fill ttbar histos
// 	TTree *ttbar_tree = fSamples[TTJets]->getTree();
// 	ttbar_tree->ResetBranchAddresses();
// 	Init(ttbar_tree);

// 	// Event loop muons and electrons
// 	if (fChain == 0) return;
// 	Long64_t tt_nentries = fChain->GetEntriesFast();
// 	for (Long64_t jentry=0; jentry<tt_nentries;jentry++) {
// 		printProgress(jentry, tt_nentries, fSamples[TTJets]->name);

// 		ttbar_tree->GetEntry(jentry);
// 		//Long64_t ientry = LoadTree(jentry);
// 		//if (ientry < 0) break;
		
// 		int muind1(-1), muind2(-1);
// 		if(hasLooseMuons(muind1, muind2) > 0) {
// 			// Common event selections
// 			if(passesJet50Cut()) { // make trigger 100% efficient
// 				// Common object selections
// 				if( isLooseMuon(muind1) ) { //&& isLooseMuon(muind2)){
// 					if(MuPt[muind1] > fC_minMu2pt) {
// 					// match muons to susy or vector boson
// 						if( (MuGenMType[muind1] == 9  || MuGenGMType[muind1] == 9 ) || (MuGenMType[muind1]  == 4 && abs(MuGenMID[muind1]) != 21 && abs(MuGenMID[muind1]) != 22) ) {
// 							hiso_ttbar->Fill(MuPFIso[muind1]); // MARC
// 						}
// 						//if( (MuGenMType[muind2] == 9  || MuGenGMType[muind2] == 9 ) || (MuGenMType[muind2]  == 4 && abs(MuGenMID[muind2]) != 21 && abs(MuGenMID[muind2]) != 22) ) {
// 						//	hiso_ttbar->Fill(MuIso[muind2]);
// 						//}
// 					}
// 				}
// 			}
// 		}
// 		int elind1(-1), elind2(-1);
// 		if(hasLooseElectrons(elind1, elind2) > 0) {
// 			// Common event selections
// 			if(passesJet50Cut()) { // make trigger 100% efficient
// 				// Common object selections
// 				if(isLooseElectron(elind1) ) { //&& isLooseElectron(elind2)) {
// 					if(ElPt[elind1] > fC_minEl2pt) {
// 						// match electrons to susy particle or vector boson
// 						if( (ElGenMType[elind1] == 9  || ElGenGMType[elind1] == 9 ) || (ElGenMType[elind1]  == 4 && abs(ElGenMID[elind1]) != 21 && abs(ElGenMID[elind1]) != 22) ) {
// 							el_hiso_ttbar->Fill(ElPFIso[elind1]); // MARC
// 						}
// 						//if( (ElGenMType[elind2] == 9  || ElGenGMType[elind2] == 9 ) || (ElGenMType[elind2]  == 4 && abs(ElGenMID[elind2]) != 21 && abs(ElGenMID[elind2]) != 22) ) {
// 						//	el_hiso_ttbar->Fill(ElPFIso[elind2]);
// 						//}
// 					}
// 				}
// 			}
// 		}
// 	} // end loop over all events
// 	fSamples[TTJets]->cleanUp();
// 	cout << endl;
// 	////////////////////////////////////////////////////

// 	////////////////////////////////////////////////////
// 	// Fill signal histos
// 	TTree *signal_tree = fSamples[LM6]->getTree();
// 	signal_tree->ResetBranchAddresses();
// 	Init(signal_tree);
// 	// Event loop

// 	// Event loop muons and electrons
// 	if (fChain == 0) return;
// 	Long64_t sig_nentries = fChain->GetEntriesFast();
// 	for (Long64_t jentry=0; jentry<sig_nentries;jentry++) {
// 		printProgress(jentry, sig_nentries, fSamples[LM6]->name);

// 		signal_tree->GetEntry(jentry);
// 		//Long64_t ientry = LoadTree(jentry);
// 		//if (ientry < 0) break;
		
// 		int muind1(-1), muind2(-1);
// 		if(hasLooseMuons(muind1, muind2) > 0) {
// 			// Common event selections
// 			if(passesJet50Cut()) { // make trigger 100% efficient
// 				// Common object selections
// 				if( isLooseMuon(muind1) ) { //&& isLooseMuon(muind2)){
// 					if(MuPt[muind1] > fC_minMu2pt) {
// 					// match muons to susy or vector boson
// 						if( (MuGenMType[muind1] == 9  || MuGenGMType[muind1] == 9 ) || (MuGenMType[muind1]  == 4 && abs(MuGenMID[muind1]) != 21 && abs(MuGenMID[muind1]) != 22) ) {
// 							hiso_signal->Fill(MuPFIso[muind1]); // MARC
// 						}
// 						//if( (MuGenMType[muind2] == 9  || MuGenGMType[muind2] == 9 ) || (MuGenMType[muind2]  == 4 && abs(MuGenMID[muind2]) != 21 && abs(MuGenMID[muind2]) != 22) ) {
// 						//	hiso_signal->Fill(MuIso[muind2]);
// 						//}
// 					}
// 				}
// 			}
// 		}
// 		int elind1(-1), elind2(-1);
// 		if(hasLooseElectrons(elind1, elind2) > 0) {
// 			// Common event selections
// 			if(passesJet50Cut()) { // make trigger 100% efficient
// 				// Common object selections
// 				if(isLooseElectron(elind1) ){ //&& isLooseElectron(elind2)) {
// 					if(ElPt[elind1] > fC_minEl2pt) {
// 						// match electrons to susy particle or vector boson
// 						if( (ElGenMType[elind1] == 9  || ElGenGMType[elind1] == 9 ) || (ElGenMType[elind1]  == 4 && abs(ElGenMID[elind1]) != 21 && abs(ElGenMID[elind1]) != 22) ) {
// 							el_hiso_signal->Fill(ElPFIso[elind1]); // MARC
// 						}
// 						//if( (ElGenMType[elind2] == 9  || ElGenGMType[elind2] == 9 ) || (ElGenMType[elind2]  == 4 && abs(ElGenMID[elind2]) != 21 && abs(ElGenMID[elind2]) != 22) ) {
// 						//	el_hiso_signal->Fill(ElPFIso[elind2]);
// 						//}
// 					}
// 				}
// 			}
// 		}
// 	} // end loop over all events
// 	fSamples[LM6]->cleanUp();
// 	cout << endl;
// 	////////////////////////////////////////////////////
	
// 	//--------------------------------------------------------------------------------------
// 	// Format and calculate muon histos
// 	hiso_ttbar->SetXTitle("rel Iso muons");
// 	hiso_ttbar->SetLineWidth(2);
// 	//hiso_ttbar->SetFillColor(kWhite);
// 	//hiso_ttbar->SetLineColor(kRed);
// 	hiso_ttbar->SetMarkerStyle(20);
// 	hiso_ttbar->SetMarkerColor(kBlue-3);
// 	hiso_ttbar->SetMarkerSize(1.15);
		
// 	//hiso_signal[i]->SetXTitle(convertVarName("MuIso[0]"));
// 	hiso_signal->SetLineWidth(2);
// 	//hiso_signal->SetLineColor(kBlue);
// 	hiso_signal->SetMarkerStyle(21);
// 	hiso_signal->SetMarkerColor(kRed+1);
// 	hiso_signal->SetMarkerSize(1.15);
		
// 	//double max1 = hiso_mc_s[i]->GetMaximum();
// 	//double max2 = hiso_data[i]->GetMaximum();
// 	//double max = max1>max2?max1:max2;
// 	//hiso_mc_s[i]->SetMaximum(1.5*max);
// 	//hiso_data[i]->SetMaximum(1.5*max);

// 	int bin0   = hiso_ttbar->FindBin(0.0);
// 	int bin015 = hiso_ttbar->FindBin(0.15) - 1; // bins start at lower edge...
// 	int bin1   = hiso_ttbar->FindBin(1.0)  - 1;
// 	float ratio_signal = hiso_signal ->Integral(bin0, bin015) / hiso_signal ->Integral(bin0, bin1);
// 	float ratio_ttbar  = hiso_ttbar  ->Integral(bin0, bin015) / hiso_ttbar  ->Integral(bin0, bin1);

	
// 	hiso_signal ->Scale(1/hiso_signal->Integral());
// 	hiso_ttbar  ->Scale(1/hiso_ttbar ->Integral());
	

// 	TCanvas *c_temp = new TCanvas("MuIso" , "Muon Isolation in TTJets and Signal", 0, 0, 800, 600);
// 	c_temp->cd();

// 	TLegend *leg = new TLegend(0.30,0.65,0.55,0.88);
// 	// TLegend *leg = new TLegend(0.75,0.60,0.89,0.88);
// 	leg->AddEntry(hiso_signal, "Signal","p");
// 	leg->AddEntry(hiso_ttbar , "TTbar","p");
// 	leg->SetFillStyle(0);
// 	leg->SetTextFont(42);
// 	leg->SetBorderSize(0);

// 	gPad->SetLogy();
// 	hiso_ttbar->Draw("pe");
// 	hiso_signal->Draw("same");
// 	leg->Draw();
// 	lat->DrawLatex(0.75,0.92, Form("L_{int.} = %2.1f fb^{-1}", fLumiNorm/1000.));
// 	lat->SetTextColor(kRed+1);
// 	lat->DrawLatex(0.70,0.70, Form("R^{T/L}_{Signal} = %4.2f", ratio_signal));
// 	lat->SetTextColor(kBlue-3);
// 	lat->DrawLatex(0.70,0.60, Form("R^{T/L}_{TTbar} = %4.2f", ratio_ttbar));
	

// 	Util::PrintPDF(c_temp, "MuIso" , fOutputDir + fOutputSubDir);

// 	//--------------------------------------------------------------------------------------
// 	// Format and calculate electron histos
// 	el_hiso_ttbar->SetXTitle("rel Iso electrons");
// 	el_hiso_ttbar->SetLineWidth(2);
// 	//hiso_ttbar->SetFillColor(kWhite);
// 	//hiso_ttbar->SetLineColor(kRed);
// 	el_hiso_ttbar->SetMarkerStyle(20);
// 	el_hiso_ttbar->SetMarkerColor(kBlue-3);
// 	el_hiso_ttbar->SetMarkerSize(1.15);

// 	//hiso_signal[i]->SetXTitle(convertVarName("MuIso[0]"));
// 	el_hiso_signal->SetLineWidth(2);
// 	//hiso_signal->SetLineColor(kBlue);
// 	el_hiso_signal->SetMarkerStyle(21);
// 	el_hiso_signal->SetMarkerColor(kRed+1);
// 	el_hiso_signal->SetMarkerSize(1.15);
		
// 	//double max1 = el_hiso_mc_s[i]->GetMaximum();
// 	//double max2 = el_hiso_data[i]->GetMaximum();
// 	//double max = max1>max2?max1:max2;
// 	//el_hiso_mc_s[i]->SetMaximum(1.5*max);
// 	//el_hiso_data[i]->SetMaximum(1.5*max);

// 	int el_bin0   = el_hiso_ttbar->FindBin(0.0);
// 	int el_bin015 = el_hiso_ttbar->FindBin(0.15) - 1; // bins start at lower edge...
// 	//int el_bin1   = el_hiso_ttbar->FindBin(1.0)  - 1;
// 	int el_bin1   = el_hiso_ttbar->FindBin(0.6)  - 1;
// 	float el_ratio_signal = el_hiso_signal ->Integral(el_bin0, el_bin015) / el_hiso_signal ->Integral(el_bin0, el_bin1);
// 	float el_ratio_ttbar  = el_hiso_ttbar  ->Integral(el_bin0, el_bin015) / el_hiso_ttbar  ->Integral(el_bin0, el_bin1);

	
// 	el_hiso_signal ->Scale(1/el_hiso_signal->Integral());
// 	el_hiso_ttbar  ->Scale(1/el_hiso_ttbar ->Integral());
	

// 	TCanvas *el_temp = new TCanvas("ElIso" , "Electron Isolation in TTJets and Signal", 0, 0, 800, 600);
// 	el_temp->cd();

// 	TLegend *el_leg = new TLegend(0.30,0.65,0.55,0.88);
// 	el_leg->AddEntry(el_hiso_signal, "Signal","p");
// 	el_leg->AddEntry(el_hiso_ttbar , "TTbar","p");
// 	el_leg->SetFillStyle(0);
// 	el_leg->SetTextFont(42);
// 	el_leg->SetBorderSize(0);

// 	gPad->SetLogy();
// 	el_hiso_ttbar->Draw("");
// 	el_hiso_signal->Draw("same");
// 	el_leg->Draw();
// 	el_lat->DrawLatex(0.75,0.92, Form("L_{int.} = %2.1f fb^{-1}", fLumiNorm/1000.));
// 	el_lat->SetTextColor(kRed+1);
// 	el_lat->DrawLatex(0.70,0.70, Form("R^{T/L}_{Signal} = %4.2f", el_ratio_signal));
// 	el_lat->SetTextColor(kBlue-3);
// 	el_lat->DrawLatex(0.70,0.60, Form("R^{T/L}_{TTbar} = %4.2f", el_ratio_ttbar));
	

// 	// Util::PrintNoEPS(c_temp, "MuIso" + IsoPlots::sel_name[i], fOutputDir + fOutputSubDir, NULL);
// 	Util::PrintPDF(el_temp, "ElIso" , fOutputDir + fOutputSubDir);
// }

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

	TH1D    *hiso_data_pt [gNSels][gNMuFPtBins];
	TH1D    *hiso_mc_pt   [gNSels][gNMuFPtBins];
	TH1D    *hiso_ttbar_pt[gNSels][gNMuFPtBins];
	TH1D    *hiso_qcd_pt  [gNSels][gNMuFPtBins];
	TH1D    *hiso_ttj_pt  [gNSels][gNMuFPtBins];
	TH1D    *hiso_ewk_pt  [gNSels][gNMuFPtBins];
	TH1D    *hiso_rare_pt [gNSels][gNMuFPtBins];
	TH1D    *hiso_db_pt   [gNSels][gNMuFPtBins];
	THStack *hiso_mc_pt_s [gNSels][gNMuFPtBins];

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

		for(int k = 0; k < gNMuFPtBins; ++k){
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
		
		for(int k = 0; k < gNMuFPtBins; ++k){
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
// LUKAS		TTree *tree = fSamples[TTJets]->getTree();
		
		
// LUKAS		TFile * file_ = new TFile("dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/mdunser/SSDLTrees/2012/Oct26/TT_CT10_TuneZ2star_8TeV-powheg-tauola-Summer12_DR53X-PU_S10_START53_V7A-v1/output_0.root", "READ", "file_");
		TFile * file_ = new TFile("/shome/lbaeni/top/CMSSW_5_3_2_patch4/src/DiLeptonAnalysis/NTupleProducer/macros/TT_CT10_TuneZ2star_8TeV-powheg-tauola-Summer12_DR53X-PU_S10_START53_V7A-v1_output_0.root", "READ", "file_");
		TTree * tree = (TTree *) file_->Get("Analysis");

		// Event loop
		tree->ResetBranchAddresses();
		Init(tree);

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
			if(MuPt[muind1] > gMuFPtBins[gNMuFPtBins]) continue;

			// Select genmatched fake muons
			if(isPromptMuon(muind1)) continue;

			////////////////////////////////////////////////////
			// MOST LOOSE SELECTION
			hiso_ttbar[0]->Fill(MuPFIso[muind1]);
			for(size_t k = 0; k < gNMuFPtBins; ++k){
				if(MuPt[muind1] < gMuFPtBins[k]) continue;
				if(MuPt[muind1] > gMuFPtBins[k+1]) continue;
				hiso_ttbar_pt[0][k]->Fill(MuPFIso[muind1]);
			}
			for(size_t k = 0; k < gNNVrtxBins; ++k){
				if(NVrtx < gNVrtxBins[k]) continue;
				if(NVrtx > gNVrtxBins[k+1]) continue;
				hiso_ttbar_nv[1][k]->Fill(MuPFIso[muind1]); // MARC
			}

			////////////////////////////////////////////////////
			// SIGNAL SUPPRESSED SELECTION
			if(isSigSupMuEvent()){
				hiso_ttbar[1]->Fill(MuPFIso[muind1]);
				for(size_t k = 0; k < gNMuFPtBins; ++k){
					if(MuPt[muind1] < gMuFPtBins[k]) continue;
					if(MuPt[muind1] > gMuFPtBins[k+1]) continue;
					hiso_ttbar_pt[1][k]->Fill(MuPFIso[muind1]); // MARC
				}
				for(size_t k = 0; k < gNNVrtxBins; ++k){
					if(NVrtx < gNVrtxBins[k]) continue;
					if(NVrtx > gNVrtxBins[k+1]) continue;
					hiso_ttbar_nv[1][k]->Fill(MuPFIso[muind1]); // MARC
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
			// 	for(size_t k = 0; k < gNMuFPtBins; ++k){
			// 		if(MuPt[fakemu] < gMuFPtBins[k]) continue;
			// 		if(MuPt[fakemu] > gMuFPtBins[k+1]) continue;
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
	} // end dottbar
	////////////////////////////////////////////////////
	
	// Create plots
//	vector<int> mcsamples = fMCBG;
	vector<int> mcsamples = fMCBGMuEnr;
	// vector<int> datasamples = fJMData;
	vector<int> datasamples = fMuData;

	for(size_t i = 0; i < gNSels; ++i){
		fOutputSubDir = "Isolation/Muons/";
		hiso_data[i]->SetXTitle(convertVarName("MuPFIso[0]"));
		hiso_data[i]->SetLineWidth(3);
		hiso_data[i]->SetLineColor(kBlack);
		hiso_data[i]->SetMarkerStyle(8);
		hiso_data[i]->SetMarkerColor(kBlack);
		hiso_data[i]->SetMarkerSize(1.2);
		
		hiso_ttbar[i]->SetXTitle(convertVarName("MuPFIso[0]"));
		hiso_ttbar[i]->SetLineWidth(3);
		hiso_ttbar[i]->SetLineColor(kRed);
		hiso_ttbar[i]->SetMarkerStyle(23);
		hiso_ttbar[i]->SetMarkerColor(kRed);
		hiso_ttbar[i]->SetMarkerSize(1.3);
		
		for(int k = 0; k < gNMuFPtBins; ++k){
			hiso_data_pt[i][k]->SetXTitle(convertVarName("MuPFIso[0]"));
			hiso_data_pt[i][k]->SetLineWidth(3);
			hiso_data_pt[i][k]->SetLineColor(kBlack);
			hiso_data_pt[i][k]->SetMarkerStyle(8);
			hiso_data_pt[i][k]->SetMarkerColor(kBlack);
			hiso_data_pt[i][k]->SetMarkerSize(1.2);

			hiso_ttbar_pt[i][k]->SetXTitle(convertVarName("MuPFIso[0]"));
			hiso_ttbar_pt[i][k]->SetLineWidth(3);
			hiso_ttbar_pt[i][k]->SetLineColor(kRed);
			hiso_ttbar_pt[i][k]->SetMarkerStyle(23);
			hiso_ttbar_pt[i][k]->SetMarkerColor(kRed);
			hiso_ttbar_pt[i][k]->SetMarkerSize(1.3);
		}
		for(int k = 0; k < gNNVrtxBins; ++k){
			hiso_data_nv[i][k]->SetXTitle(convertVarName("MuPFIso[0]"));
			hiso_data_nv[i][k]->SetLineWidth(3);
			hiso_data_nv[i][k]->SetLineColor(kBlack);
			hiso_data_nv[i][k]->SetMarkerStyle(8);
			hiso_data_nv[i][k]->SetMarkerColor(kBlack);
			hiso_data_nv[i][k]->SetMarkerSize(1.2);

			hiso_ttbar_nv[i][k]->SetXTitle(convertVarName("MuPFIso[0]"));
			hiso_ttbar_nv[i][k]->SetLineWidth(3);
			hiso_ttbar_nv[i][k]->SetLineColor(kRed);
			hiso_ttbar_nv[i][k]->SetMarkerStyle(23);
			hiso_ttbar_nv[i][k]->SetMarkerColor(kRed);
			hiso_ttbar_nv[i][k]->SetMarkerSize(1.3);
		}

		// Apply weights to MC histos
//		cout << "Apply weights to MC histos" << endl;
		for(size_t j = 0; j < gNSAMPLES; ++j){
			Sample *S = fSamples[j];
//			cout << "at sample: " << S->sname << "\txsec: " << S->xsec << "\tget lumi of this sample: " << S->getLumi() << "\tintegral of the iso histo: " << S->isoplots[0].hiso[i]->Integral() << endl;

			float lumiscale = fLumiNorm / S->getLumi();
			if(S->datamc == 0) continue;
			S->isoplots[0].hiso[i]->Scale(lumiscale);
			for(size_t k = 0; k < gNMuFPtBins; ++k){
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
			hiso_data[i]->SetXTitle(convertVarName("MuPFIso[0]"));
			for(int k = 0; k < gNMuFPtBins; ++k){
				hiso_data_pt[i][k]->Add(S->isoplots[0].hiso_pt[i][k]);
				hiso_data_pt[i][k]->SetXTitle(convertVarName("MuPFIso[0]"));
			}
			for(int k = 0; k < gNNVrtxBins; ++k){
				hiso_data_nv[i][k]->Add(S->isoplots[0].hiso_nv[i][k]);
				hiso_data_nv[i][k]->SetXTitle(convertVarName("MuPFIso[0]"));
			}
		}

		// Scale to get equal integrals
//		cout << "Scale to get equal integrals" << endl;
		float intscale(0.);
		float intscale_pt[gNMuFPtBins];
		float intscale_nv[gNNVrtxBins];
		for(size_t j = 0; j < mcsamples.size(); ++j){
			Sample *S = fSamples[mcsamples[j]];
			intscale += S->isoplots[0].hiso[i]->Integral();
//			cout << "at sample: " << S->sname << "\txsec: " << S->xsec << "\tget lumi of this sample: " << S->getLumi() << "\tintegral of the iso histo: " << S->isoplots[0].hiso[i]->Integral() << endl;
			for(int k = 0; k < gNMuFPtBins; ++k){
				intscale_pt[k] += S->isoplots[0].hiso_pt[i][k]->Integral();
			}
			for(int k = 0; k < gNNVrtxBins; ++k){
				intscale_nv[k] += S->isoplots[0].hiso_nv[i][k]->Integral();
			}
		}
		intscale = hiso_data[i]->Integral() / intscale;
		for(size_t j = 0; j < gNMuFPtBins; ++j) intscale_pt[j] = hiso_data_pt[i][j]->Integral() / intscale_pt[j];
		for(size_t j = 0; j < gNNVrtxBins; ++j) intscale_nv[j] = hiso_data_nv[i][j]->Integral() / intscale_nv[j];
		
		for(size_t j = 0; j < mcsamples.size(); ++j){
			Sample *S = fSamples[mcsamples[j]];
			S->isoplots[0].hiso[i]->Scale(intscale);
			for(int k = 0; k < gNMuFPtBins; ++k) S->isoplots[0].hiso_pt[i][k]->Scale(intscale_pt[k]);
			for(int k = 0; k < gNNVrtxBins; ++k) S->isoplots[0].hiso_nv[i][k]->Scale(intscale_nv[k]);
		}
		if(dottbar){
			hiso_ttbar[i]->Scale(hiso_data[i]->Integral() / hiso_ttbar[i]->Integral());
			for(int k = 0; k < gNMuFPtBins; ++k) hiso_ttbar_pt[i][k]->Scale(hiso_data_pt[i][k]->Integral() / hiso_ttbar_pt[i][k]->Integral());
			for(int k = 0; k < gNNVrtxBins; ++k) hiso_ttbar_nv[i][k]->Scale(hiso_data_nv[i][k]->Integral() / hiso_ttbar_nv[i][k]->Integral());			
		}
		

		// Fill MC stacks
//		cout << "Fill MC stacks" << endl;
		for(size_t j = 0; j < mcsamples.size();   ++j){
			Sample *S = fSamples[mcsamples[j]];
			TString s_name = S->sname;
//			cout << "at sample: " << S->sname << "\txsec: " << S->xsec << "\tget lumi of this sample: " << S->getLumi() << "\tintegral of the iso histo: " << S->isoplots[0].hiso[i]->Integral() << endl;
			hiso_mc  [i]->Add(S->isoplots[0].hiso[i]);
			// sample type: QCD = 1 , Top = 2, EWK = 3 , Rare = 4 , DB = 5
			if ( S->getType() == 1) hiso_qcd[i] ->Add( S->isoplots[0].hiso[i] );
			if ( S->getType() == 2) hiso_ttj[i] ->Add( S->isoplots[0].hiso[i] );
			if ( S->getType() == 3) hiso_ewk[i] ->Add( S->isoplots[0].hiso[i] );
			if ( S->getType() == 4) hiso_rare[i]->Add( S->isoplots[0].hiso[i] );
			if ( S->getType() == 5) hiso_db[i]  ->Add( S->isoplots[0].hiso[i] );

			for(int k = 0; k < gNMuFPtBins; ++k){
				if ( S->getType() == 1) hiso_qcd_pt [i][k]->Add( S->isoplots[0].hiso_pt[i][k] );
				if ( S->getType() == 2) hiso_ttj_pt [i][k]->Add( S->isoplots[0].hiso_pt[i][k] );
				if ( S->getType() == 3) hiso_ewk_pt [i][k]->Add( S->isoplots[0].hiso_pt[i][k] );
				if ( S->getType() == 4) hiso_rare_pt[i][k]->Add( S->isoplots[0].hiso_pt[i][k] );
				if ( S->getType() == 5) hiso_db_pt  [i][k]->Add( S->isoplots[0].hiso_pt[i][k] );
			}
			for(int k = 0; k < gNNVrtxBins; ++k){
				if ( S->getType() == 1) hiso_qcd_nv [i][k]->Add( S->isoplots[0].hiso_nv[i][k] );
				if ( S->getType() == 2) hiso_ttj_nv [i][k]->Add( S->isoplots[0].hiso_nv[i][k] );
				if ( S->getType() == 3) hiso_ewk_nv [i][k]->Add( S->isoplots[0].hiso_nv[i][k] );
				if ( S->getType() == 4) hiso_rare_nv[i][k]->Add( S->isoplots[0].hiso_nv[i][k] );
				if ( S->getType() == 5) hiso_db_nv  [i][k]->Add( S->isoplots[0].hiso_nv[i][k] );
			}

		}
		hiso_mc_s[i]->Add(hiso_qcd[i]);
		hiso_mc_s[i]->Add(hiso_db[i]);
		hiso_mc_s[i]->Add(hiso_ewk[i]);
		hiso_mc_s[i]->Add(hiso_rare[i]);
		hiso_mc_s[i]->Add(hiso_ttj[i]);
		hiso_mc_s[i]->Draw("goff");
		hiso_mc_s[i]->GetXaxis()->SetTitle(convertVarName("MuPFIso[0]"));

		for(int k = 0; k < gNMuFPtBins; ++k){
			hiso_mc_pt_s[i][k]->Add(hiso_qcd_pt[i][k]);
			hiso_mc_pt_s[i][k]->Add(hiso_db_pt[i][k]);
			hiso_mc_pt_s[i][k]->Add(hiso_ewk_pt[i][k]);
			hiso_mc_pt_s[i][k]->Add(hiso_rare_pt[i][k]);
			hiso_mc_pt_s[i][k]->Add(hiso_ttj_pt[i][k]);
			hiso_mc_pt_s[i][k]->Draw("goff");
			hiso_mc_pt_s[i][k]->GetXaxis()->SetTitle(convertVarName("MuPFIso[0]"));
		}

		for(int k = 0; k < gNNVrtxBins; ++k){
			hiso_mc_nv_s[i][k]->Add(hiso_qcd_nv[i][k]);
			hiso_mc_nv_s[i][k]->Add(hiso_db_nv[i][k]);
			hiso_mc_nv_s[i][k]->Add(hiso_ewk_nv[i][k]);
			hiso_mc_nv_s[i][k]->Add(hiso_rare_nv[i][k]);
			hiso_mc_nv_s[i][k]->Add(hiso_ttj_nv[i][k]);
			hiso_mc_nv_s[i][k]->Draw("goff");
			hiso_mc_nv_s[i][k]->GetXaxis()->SetTitle(convertVarName("MuPFIso[0]"));
		}
		// Fill MC stacks
		//for(size_t j = 0; j < mcsamples.size();   ++j){
		//	Sample *S = fSamples[mcsamples[j]];
		//	hiso_mc  [i]->Add(S->isoplots[0].hiso[i]);
		//	hiso_mc_s[i]->Add(S->isoplots[0].hiso[i]);
		//	hiso_mc_s[i]->Draw("goff");
		//	hiso_mc_s[i]->GetXaxis()->SetTitle(convertVarName("MuPFIso[0]"));
		//	for(int k = 0; k < gNMuFPtBins; ++k){
		//		hiso_mc_pt  [i][k]->Add(S->isoplots[0].hiso_pt[i][k]);
		//		hiso_mc_pt_s[i][k]->Add(S->isoplots[0].hiso_pt[i][k]);
		//		hiso_mc_pt_s[i][k]->Draw("goff");
		//		hiso_mc_pt_s[i][k]->GetXaxis()->SetTitle(convertVarName("MuPFIso[0]"));
		//	}
		//	for(int k = 0; k < gNNVrtxBins; ++k){
		//		hiso_mc_nv  [i][k]->Add(S->isoplots[0].hiso_nv[i][k]);
		//		hiso_mc_nv_s[i][k]->Add(S->isoplots[0].hiso_nv[i][k]);
		//		hiso_mc_nv_s[i][k]->Draw("goff");
		//		hiso_mc_nv_s[i][k]->GetXaxis()->SetTitle(convertVarName("MuPFIso[0]"));
		//	}
		//}

		double max1 = hiso_mc_s[i]->GetMaximum();
		double max2 = hiso_data[i]->GetMaximum();
		double max = max1>max2?max1:max2;
		hiso_mc_s[i]->SetMaximum(1.5*max);
		hiso_data[i]->SetMaximum(1.5*max);

		int bin0   = hiso_data[i]->FindBin(0.0);
		int bin015 = hiso_data[i]->FindBin(0.1) - 1; // bins start at lower edge...
		int bin1   = hiso_data[i]->FindBin(1.0)  - 1;
		// printf("bin 0: %3.0i bin015: %3.0i bin1: %3.0i \n", bin0, bin015, bin1);
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
		lat->DrawLatex(0.70,0.92, Form("L_{int.} = %5.2f fb^{-1}", fLumiNorm/1000.));
		lat->DrawLatex(0.75,0.85, Form("R^{T/L}_{Data}  = %4.2f", ratio_data));
		lat->DrawLatex(0.75,0.80, Form("R^{T/L}_{MC}   = %4.2f", ratio_mc));
		lat->SetTextColor(kRed);
		if(dottbar) lat->DrawLatex(0.75,0.75, Form("R^{T/L}_{TTbar} = %4.2f", ratio_ttbar));
		lat->SetTextColor(kBlack);
		

		// Util::PrintNoEPS(c_temp, "MuIso" + IsoPlots::sel_name[i], fOutputDir + fOutputSubDir, NULL);
		Util::PrintPDF(c_temp, "MuIso" + IsoPlots::sel_name[i], fOutputDir + fOutputSubDir);

		for(int k = 0; k < gNMuFPtBins; ++k){
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
			lat->DrawLatex(0.20,0.92, Form("p_{T}(#mu) %3.0f - %3.0f GeV", getFPtBins(Muon)[k], getFPtBins(Muon)[k+1]));
			lat->DrawLatex(0.70,0.92, Form("L_{int.} = %5.2f fb^{-1}", fLumiNorm/1000.));
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
			lat->DrawLatex(0.70,0.92, Form("L_{int.} = %5.2f fb^{-1}", fLumiNorm/1000.));
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

	TH1D    *hiso_data_pt [gNSels][gNElFPtBins];
	TH1D    *hiso_mc_pt   [gNSels][gNElFPtBins];
	TH1D    *hiso_ttbar_pt[gNSels][gNElFPtBins];
	TH1D    *hiso_qcd_pt  [gNSels][gNElFPtBins];
	TH1D    *hiso_ttj_pt  [gNSels][gNElFPtBins];
	TH1D    *hiso_ewk_pt  [gNSels][gNElFPtBins];
	TH1D    *hiso_rare_pt [gNSels][gNElFPtBins];
	TH1D    *hiso_db_pt   [gNSels][gNElFPtBins];
	THStack *hiso_mc_pt_s [gNSels][gNElFPtBins];

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

		for(int k = 0; k < gNElFPtBins; ++k){
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

		for(int k = 0; k < gNElFPtBins; ++k){
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
		Init(tree);

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
			if(ElPt[elind1] > gElFPtBins[gNElFPtBins]) continue;


			////////////////////////////////////////////////////
			// MOST LOOSE SELECTION
			hiso_ttbar[0]->Fill(ElPFIso[elind1]);
			for(size_t k = 0; k < gNElFPtBins; ++k){
				if(ElPt[elind1] < gElFPtBins[k]) continue;
				if(ElPt[elind1] > gElFPtBins[k+1]) continue;
				hiso_ttbar_pt[0][k]->Fill(ElPFIso[elind1]); // MARC
			}
			for(size_t k = 0; k < gNNVrtxBins; ++k){
				if(NVrtx < gNVrtxBins[k]) continue;
				if(NVrtx > gNVrtxBins[k+1]) continue;
				hiso_ttbar_nv[1][k]->Fill(ElPFIso[elind1]); // MARC
			}

			////////////////////////////////////////////////////
			// SIGNAL SUPPRESSED SELECTION
			if(isSigSupElEvent()){
				hiso_ttbar[1]->Fill(ElPFIso[elind1]);
				for(size_t k = 0; k < gNElFPtBins; ++k){
					if(ElPt[elind1] < gElFPtBins[k]) continue;
					if(ElPt[elind1] > gElFPtBins[k+1]) continue;
					hiso_ttbar_pt[1][k]->Fill(ElPFIso[elind1]); // MARC
				}
				for(size_t k = 0; k < gNNVrtxBins; ++k){
					if(NVrtx < gNVrtxBins[k]) continue;
					if(NVrtx > gNVrtxBins[k+1]) continue;
					hiso_ttbar_nv[1][k]->Fill(ElPFIso[elind1]); // MARC
				}
			}
			////////////////////////////////////////////////////
		}
	// attention, this only loops on one of the ttbar samples
		fSamples[TTJets]->cleanUp();
		cout << endl;
	}
	////////////////////////////////////////////////////

	// Create plots
//	vector<int> mcsamples = fMCBG;
	vector<int> mcsamples = fMCBGEMEnr;
	vector<int> datasamples = fEGData;

	for(size_t i = 0; i < gNSels; ++i){
		fOutputSubDir = "Isolation/Electrons/";
		hiso_data[i]->SetXTitle(convertVarName("ElPFIso[0]"));
		hiso_data[i]->SetLineWidth(3);
		hiso_data[i]->SetLineColor(kBlack);
		hiso_data[i]->SetMarkerStyle(8);
		hiso_data[i]->SetMarkerColor(kBlack);
		hiso_data[i]->SetMarkerSize(1.2);

		hiso_ttbar[i]->SetXTitle(convertVarName("ElPFIso[0]"));
		hiso_ttbar[i]->SetLineWidth(3);
		hiso_ttbar[i]->SetLineColor(kRed);
		hiso_ttbar[i]->SetMarkerStyle(23);
		hiso_ttbar[i]->SetMarkerColor(kRed);
		hiso_ttbar[i]->SetMarkerSize(1.3);
		
		for(int k = 0; k < gNElFPtBins; ++k){
			hiso_data_pt[i][k]->SetXTitle(convertVarName("ElPFIso[0]"));
			hiso_data_pt[i][k]->SetLineWidth(3);
			hiso_data_pt[i][k]->SetLineColor(kBlack);
			hiso_data_pt[i][k]->SetMarkerStyle(8);
			hiso_data_pt[i][k]->SetMarkerColor(kBlack);
			hiso_data_pt[i][k]->SetMarkerSize(1.2);

			hiso_ttbar_pt[i][k]->SetXTitle(convertVarName("ElPFIso[0]"));
			hiso_ttbar_pt[i][k]->SetLineWidth(3);
			hiso_ttbar_pt[i][k]->SetLineColor(kRed);
			hiso_ttbar_pt[i][k]->SetMarkerStyle(23);
			hiso_ttbar_pt[i][k]->SetMarkerColor(kRed);
			hiso_ttbar_pt[i][k]->SetMarkerSize(1.3);
		}
		for(int k = 0; k < gNNVrtxBins; ++k){
			hiso_data_nv[i][k]->SetXTitle(convertVarName("ElPFIso[0]"));
			hiso_data_nv[i][k]->SetLineWidth(3);
			hiso_data_nv[i][k]->SetLineColor(kBlack);
			hiso_data_nv[i][k]->SetMarkerStyle(8);
			hiso_data_nv[i][k]->SetMarkerColor(kBlack);
			hiso_data_nv[i][k]->SetMarkerSize(1.2);

			hiso_ttbar_nv[i][k]->SetXTitle(convertVarName("ElPFIso[0]"));
			hiso_ttbar_nv[i][k]->SetLineWidth(3);
			hiso_ttbar_nv[i][k]->SetLineColor(kRed);
			hiso_ttbar_nv[i][k]->SetMarkerStyle(23);
			hiso_ttbar_nv[i][k]->SetMarkerColor(kRed);
			hiso_ttbar_nv[i][k]->SetMarkerSize(1.3);
		}

		// Apply weights to MC histos
		for(size_t j = 0; j < gNSAMPLES; ++j){
			Sample *S = fSamples[j];
			float lumiscale = fLumiNorm / S->getLumi();
			if(S->datamc == 0) continue;
			S->isoplots[1].hiso[i]->Scale(lumiscale);
			for(size_t k = 0; k < gNElFPtBins; ++k){
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
			hiso_data[i]->SetXTitle(convertVarName("ElPFIso[0]"));
			for(int k = 0; k < gNElFPtBins; ++k){
				hiso_data_pt[i][k]->Add(S->isoplots[1].hiso_pt[i][k]);
				hiso_data_pt[i][k]->SetXTitle(convertVarName("ElPFIso[0]"));
			}
			for(int k = 0; k < gNNVrtxBins; ++k){
				hiso_data_nv[i][k]->Add(S->isoplots[1].hiso_nv[i][k]);
				hiso_data_nv[i][k]->SetXTitle(convertVarName("ElPFIso[0]"));
			}
		}

		// Scale to get equal integrals
		float intscale(0.);
		float intscale_pt[gNElFPtBins];
		float intscale_nv[gNNVrtxBins];
		for(size_t j = 0; j < mcsamples.size();   ++j){
			Sample *S = fSamples[mcsamples[j]];
			intscale += S->isoplots[1].hiso[i]->Integral();
			for(int k = 0; k < gNElFPtBins; ++k){
				intscale_pt[k] += S->isoplots[1].hiso_pt[i][k]->Integral();
			}
			for(int k = 0; k < gNNVrtxBins; ++k){
				intscale_nv[k] += S->isoplots[1].hiso_nv[i][k]->Integral();
			}
		}
		intscale = hiso_data[i]->Integral() / intscale;
		for(size_t j = 0; j < gNElFPtBins; ++j) intscale_pt[j] = hiso_data_pt[i][j]->Integral() / intscale_pt[j];
		for(size_t j = 0; j < gNNVrtxBins; ++j) intscale_nv[j] = hiso_data_nv[i][j]->Integral() / intscale_nv[j];
		
		for(size_t j = 0; j < mcsamples.size();   ++j){
			Sample *S = fSamples[mcsamples[j]];			
			S->isoplots[1].hiso[i]->Scale(intscale);
			for(int k = 0; k < gNElFPtBins; ++k){
				S->isoplots[1].hiso_pt[i][k]->Scale(intscale_pt[k]);
			}
			for(int k = 0; k < gNNVrtxBins; ++k){
				S->isoplots[1].hiso_nv[i][k]->Scale(intscale_nv[k]);
			}
		}
		if(dottbar){
			hiso_ttbar[i]->Scale(hiso_data[i]->Integral() / hiso_ttbar[i]->Integral());
			for(int k = 0; k < gNElFPtBins; ++k) hiso_ttbar_pt[i][k]->Scale(hiso_data_pt[i][k]->Integral() / hiso_ttbar_pt[i][k]->Integral());
			for(int k = 0; k < gNNVrtxBins; ++k) hiso_ttbar_nv[i][k]->Scale(hiso_data_nv[i][k]->Integral() / hiso_ttbar_nv[i][k]->Integral());			
		}

		// Fill MC stacks
		for(size_t j = 0; j < mcsamples.size();   ++j){
			Sample *S = fSamples[mcsamples[j]];
			TString s_name = S->sname;
			hiso_mc  [i]->Add(S->isoplots[1].hiso[i]);
			// sample type: QCD = 1 , Top = 2, EWK = 3 , Rare = 4 , DB = 5
			if ( S->getType() == 1) hiso_qcd[i] ->Add( S->isoplots[1].hiso[i] );
			if ( S->getType() == 2) hiso_ttj[i] ->Add( S->isoplots[1].hiso[i] );
			if ( S->getType() == 3) hiso_ewk[i] ->Add( S->isoplots[1].hiso[i] );
			if ( S->getType() == 4) hiso_rare[i]->Add( S->isoplots[1].hiso[i] );
			if ( S->getType() == 5) hiso_db[i]  ->Add( S->isoplots[1].hiso[i] );

			for(int k = 0; k < gNMuFPtBins; ++k){
				if ( S->getType() == 1) hiso_qcd_pt [i][k]->Add( S->isoplots[1].hiso_pt[i][k] );
				if ( S->getType() == 2) hiso_ttj_pt [i][k]->Add( S->isoplots[1].hiso_pt[i][k] );
				if ( S->getType() == 3) hiso_ewk_pt [i][k]->Add( S->isoplots[1].hiso_pt[i][k] );
				if ( S->getType() == 4) hiso_rare_pt[i][k]->Add( S->isoplots[1].hiso_pt[i][k] );
				if ( S->getType() == 5) hiso_db_pt  [i][k]->Add( S->isoplots[1].hiso_pt[i][k] );
			}
			for(int k = 0; k < gNNVrtxBins; ++k){
				if ( S->getType() == 1) hiso_qcd_nv [i][k]->Add( S->isoplots[1].hiso_nv[i][k] );
				if ( S->getType() == 2) hiso_ttj_nv [i][k]->Add( S->isoplots[1].hiso_nv[i][k] );
				if ( S->getType() == 3) hiso_ewk_nv [i][k]->Add( S->isoplots[1].hiso_nv[i][k] );
				if ( S->getType() == 4) hiso_rare_nv[i][k]->Add( S->isoplots[1].hiso_nv[i][k] );
				if ( S->getType() == 5) hiso_db_nv  [i][k]->Add( S->isoplots[1].hiso_nv[i][k] );
			}

		}
		hiso_mc_s[i]->Add(hiso_qcd[i]);
		hiso_mc_s[i]->Add(hiso_db[i]);
		hiso_mc_s[i]->Add(hiso_ewk[i]);
		hiso_mc_s[i]->Add(hiso_rare[i]);
		hiso_mc_s[i]->Add(hiso_ttj[i]);
		hiso_mc_s[i]->Draw("goff");
		hiso_mc_s[i]->GetXaxis()->SetTitle(convertVarName("ElPFIso[0]"));

		for(int k = 0; k < gNElFPtBins; ++k){
			hiso_mc_pt_s[i][k]->Add(hiso_qcd_pt[i][k]);
			hiso_mc_pt_s[i][k]->Add(hiso_db_pt[i][k]);
			hiso_mc_pt_s[i][k]->Add(hiso_ewk_pt[i][k]);
			hiso_mc_pt_s[i][k]->Add(hiso_rare_pt[i][k]);
			hiso_mc_pt_s[i][k]->Add(hiso_ttj_pt[i][k]);
			hiso_mc_pt_s[i][k]->Draw("goff");
			hiso_mc_pt_s[i][k]->GetXaxis()->SetTitle(convertVarName("ElPFIso[0]"));
		}

		for(int k = 0; k < gNNVrtxBins; ++k){
			hiso_mc_nv_s[i][k]->Add(hiso_qcd_nv[i][k]);
			hiso_mc_nv_s[i][k]->Add(hiso_db_nv[i][k]);
			hiso_mc_nv_s[i][k]->Add(hiso_ewk_nv[i][k]);
			hiso_mc_nv_s[i][k]->Add(hiso_rare_nv[i][k]);
			hiso_mc_nv_s[i][k]->Add(hiso_ttj_nv[i][k]);
			hiso_mc_nv_s[i][k]->Draw("goff");
			hiso_mc_nv_s[i][k]->GetXaxis()->SetTitle(convertVarName("ElPFIso[0]"));
		}
		//// Fill MC stacks
		//for(size_t j = 0; j < mcsamples.size();   ++j){
		//	Sample *S = fSamples[mcsamples[j]];			
		//	hiso_mc  [i]->Add(S->isoplots[1].hiso[i]);
		//	hiso_mc_s[i]->Add(S->isoplots[1].hiso[i]);
		//	hiso_mc_s[i]->Draw("goff");
		//	hiso_mc_s[i]->GetXaxis()->SetTitle(convertVarName("ElPFIso[0]"));
		//	for(int k = 0; k < gNElFPtBins; ++k){
		//		hiso_mc_pt  [i][k]->Add(S->isoplots[1].hiso_pt[i][k]);
		//		hiso_mc_pt_s[i][k]->Add(S->isoplots[1].hiso_pt[i][k]);
		//		hiso_mc_pt_s[i][k]->Draw("goff");
		//		hiso_mc_pt_s[i][k]->GetXaxis()->SetTitle(convertVarName("ElPFIso[0]"));
		//	}
		//	for(int k = 0; k < gNNVrtxBins; ++k){
		//		hiso_mc_nv  [i][k]->Add(S->isoplots[1].hiso_nv[i][k]);
		//		hiso_mc_nv_s[i][k]->Add(S->isoplots[1].hiso_nv[i][k]);
		//		hiso_mc_nv_s[i][k]->Draw("goff");
		//		hiso_mc_nv_s[i][k]->GetXaxis()->SetTitle(convertVarName("ElPFIso[0]"));
		//	}
		//}


		double max1 = hiso_mc_s[i]->GetMaximum();
		double max2 = hiso_data[i]->GetMaximum();
		double max = max1>max2?max1:max2;

		hiso_mc_s[i]->SetMaximum(1.2*max);
		hiso_data[i]->SetMaximum(1.2*max);

		int bin0   = hiso_data[i]->FindBin(0.0);
		int bin015 = hiso_data[i]->FindBin(0.09) - 1; // bins start at lower edge...
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
		lat->DrawLatex(0.70,0.92, Form("L_{int.} = %5.2f fb^{-1}", fLumiNorm/1000.));
		lat->DrawLatex(0.75,0.85, Form("R^{T/L}_{Data} = %4.2f", ratio_data));
		lat->DrawLatex(0.75,0.80, Form("R^{T/L}_{MC}  = %4.2f", ratio_mc));
		lat->SetTextColor(kRed);
		if(dottbar) lat->DrawLatex(0.75,0.75, Form("R^{T/L}_{TTbar} = %4.2f", ratio_ttbar));
		lat->SetTextColor(kBlack);

		// Util::PrintNoEPS(c_temp, "Iso" + IsoPlots::sel_name[i], fOutputDir + fOutputSubDir, NULL);
		Util::PrintPDF(c_temp, "ElIso" + IsoPlots::sel_name[i], fOutputDir + fOutputSubDir);

		for(int k = 0; k < gNElFPtBins; ++k){
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
			lat->DrawLatex(0.20,0.92, Form("p_{T}(e) %3.0f - %3.0f GeV", getFPtBins(Elec)[k], getFPtBins(Elec)[k+1]));
			lat->DrawLatex(0.70,0.92, Form("L_{int.} = %5.2f fb^{-1}", fLumiNorm/1000.));
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
			lat->DrawLatex(0.70,0.92, Form("L_{int.} = %5.2f fb^{-1}", fLumiNorm/1000.));
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

	vector<int> mcsamples = fMCBGEMEnr;
	vector<int> datasamples = fEGData;
	
	int nids = 6;
	TH1D    *data [nids] [gNSels];
	TH1D    *qcd  [nids] [gNSels];
	TH1D    *ttj  [nids] [gNSels];
	TH1D    *ewk  [nids] [gNSels];
	TH1D    *rare [nids] [gNSels];
	TH1D    *db   [nids] [gNSels];
	TH1D	*ttbar[nids] [gNSels];
	THStack *mc_s [nids] [gNSels];

	float intscale[nids];
	double datamax[nids], mcmax[nids], max[nids];
	double axismin[nids], axismax[nids];

	TLatex *lat = new TLatex();
	lat->SetNDC(kTRUE);
	lat->SetTextColor(kBlack);
	lat->SetTextSize(0.04);

	// loop k is over the electron ID variables: HoE, sigmaIetaIeta, dEta, dPhi, mvaID, medium WP pass bool
	for (int k =0 ; k<nids ; k++){
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
			case 4: 
				varName = "ElMVAID";
				axismin[k] = -1.;
				axismax[k] =  1.;
				break ;
			case 5: 
				varName = "ElmedWP";
				axismin[k] =  0.;
				axismax[k] =  2.;
				break ;
		}
		for(size_t i = 0; i < gNSels; ++i){
			int nbins;
			if (k != 5) nbins = IdPlots::nbins[i] ;
			else nbins = 3;
			data [k][i]  = new TH1D("El"+varName+"Data_"          + IdPlots::sel_name[i], "Electron "+varName+" in Data for "  + IdPlots::sel_name[i], nbins, axismin[k], axismax[k]);
			qcd  [k][i]  = new TH1D("El"+varName+"MC_"            + IdPlots::sel_name[i], "Electron "+varName+" in MC for "    + IdPlots::sel_name[i], nbins, axismin[k], axismax[k]);
			ttj  [k][i]  = new TH1D("El"+varName+"MC_"            + IdPlots::sel_name[i], "Electron "+varName+" in MC for "    + IdPlots::sel_name[i], nbins, axismin[k], axismax[k]);
			ewk  [k][i]  = new TH1D("El"+varName+"MC_"            + IdPlots::sel_name[i], "Electron "+varName+" in MC for "    + IdPlots::sel_name[i], nbins, axismin[k], axismax[k]);
			rare [k][i]  = new TH1D("El"+varName+"MC_"            + IdPlots::sel_name[i], "Electron "+varName+" in MC for "    + IdPlots::sel_name[i], nbins, axismin[k], axismax[k]);
			db   [k][i]  = new TH1D("El"+varName+"MC_"            + IdPlots::sel_name[i], "Electron "+varName+" in MC for "    + IdPlots::sel_name[i], nbins, axismin[k], axismax[k]);
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
				float lumiscale = fLumiNorm / S->getLumi();
				if(S->datamc == 0) continue;
				switch (k) {
					case 0: S->idplots.hhoe[i]->Scale(lumiscale); break;
					case 1: S->idplots.hsiesie[i]->Scale(lumiscale); break;
					case 2: S->idplots.hdeta[i]->Scale(lumiscale); break;
					case 3: S->idplots.hdphi[i]->Scale(lumiscale); break;
					case 4: S->idplots.hmvaid[i]->Scale(lumiscale); break;
					case 5: S->idplots.hmedwp[i]->Scale(lumiscale); break;
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
					case 4: data[k][i]->Add(S->idplots.hmvaid[i]); break;
					case 5: data[k][i]->Add(S->idplots.hmedwp[i]); break;
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
					case 4: intscale[k] += S->idplots.hmvaid[i]->Integral(); break;
					case 5: intscale[k] += S->idplots.hmedwp[i]->Integral(); break;
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
					case 4: S -> idplots.hmvaid[i]  -> Scale(intscale[k]); break;
					case 5: S -> idplots.hmedwp[i]  -> Scale(intscale[k]); break;
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
					case 4: histo = S->idplots.hmvaid[i] ; break;
					case 5: histo = S->idplots.hmedwp[i] ; break;
				}
				if ( S->getType() == 1) qcd [k][i]->Add( histo );
				if ( S->getType() == 2) ttj [k][i]->Add( histo );
				if ( S->getType() == 3) ewk [k][i]->Add( histo );
				if ( S->getType() == 4) rare[k][i]->Add( histo );
				if ( S->getType() == 5) db  [k][i]->Add( histo );
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

			if (k != 5) gPad->SetLogy();
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

void SSDLPlotter::makeNT2KinPlots(bool loglin){
	TString selname[3] = {"LooseLoose", "TightTight", "Signal"};
	for(size_t s = 0; s < 3; ++s){ // loop on selections
		fOutputSubDir = "KinematicPlots/" + selname[s];
		if(loglin) fOutputSubDir += "/log/";
		char cmd[100];
		sprintf(cmd,"mkdir -p %s%s", fOutputDir.Data(), fOutputSubDir.Data());
		system(cmd);
		
		TH1D    *hvar_data[gNKinVars];

		TH1D    *hvar_qcd  [gNKinVars];
		TH1D    *hvar_ttj  [gNKinVars];
		TH1D    *hvar_ewk  [gNKinVars];
		TH1D    *hvar_rare [gNKinVars];
		TH1D    *hvar_db   [gNKinVars];
		TH1D    *hvar_tot  [gNKinVars];

		TH1D    *hvar_rat  [gNKinVars]; // ratio
		THStack *hvar_mc_s[gNKinVars];

		TLatex *lat = new TLatex();
		lat->SetNDC(kTRUE);
		lat->SetTextColor(kBlack);
		lat->SetTextSize(0.04);

		// Create histograms
		for(size_t i = 0; i < gNKinVars; ++i){
			hvar_data[i] = new TH1D("Data_"    + KinPlots::var_name[i], KinPlots::var_name[i] + " in Data", KinPlots::nbins[i], KinPlots::xmin[i], KinPlots::xmax[i]);
			hvar_qcd [i] = new TH1D("QCD_"     + KinPlots::var_name[i], KinPlots::var_name[i] + " in MC",   KinPlots::nbins[i], KinPlots::xmin[i], KinPlots::xmax[i]);
			hvar_ttj [i] = new TH1D("TTjets_"  + KinPlots::var_name[i], KinPlots::var_name[i] + " in MC",   KinPlots::nbins[i], KinPlots::xmin[i], KinPlots::xmax[i]);
			hvar_ewk [i] = new TH1D("EWK_"     + KinPlots::var_name[i], KinPlots::var_name[i] + " in MC",   KinPlots::nbins[i], KinPlots::xmin[i], KinPlots::xmax[i]);
			hvar_rare[i] = new TH1D("Rare_"    + KinPlots::var_name[i], KinPlots::var_name[i] + " in MC",   KinPlots::nbins[i], KinPlots::xmin[i], KinPlots::xmax[i]);
			hvar_db[i]   = new TH1D("DB_"      + KinPlots::var_name[i], KinPlots::var_name[i] + " in MC",   KinPlots::nbins[i], KinPlots::xmin[i], KinPlots::xmax[i]);
			hvar_tot[i]  = new TH1D("Tot_"     + KinPlots::var_name[i], KinPlots::var_name[i] + " in MC",   KinPlots::nbins[i], KinPlots::xmin[i], KinPlots::xmax[i]);
			hvar_rat[i]  = new TH1D("Ratio_"   + KinPlots::var_name[i], KinPlots::var_name[i] + " Data/MC", KinPlots::nbins[i], KinPlots::xmin[i], KinPlots::xmax[i]);

			hvar_mc_s[i] = new THStack("MC_stacked_" + KinPlots::var_name[i], KinPlots::var_name[i] + " in MC");
		}

		// // Adjust overflow bins:
		// for(size_t i = 0; i < gNKinVars; ++i){
		// 	for(gSample j = sample_begin; j < gNSAMPLES; j=gSample(j+1)){
		// 		Int_t nbins     = fSamples[j]->kinplots[s][HighPt].hvar[i]->GetNbinsX();
		// 		Double_t binc   = fSamples[j]->kinplots[s][HighPt].hvar[i]->GetBinContent(nbins);
		// 		Double_t overfl = fSamples[j]->kinplots[s][HighPt].hvar[i]->GetBinContent(nbins+1);
		// 		fSamples[j]->kinplots[s][HighPt].hvar[i]->SetBinContent(nbins, binc + overfl);
		// 	}
		// }

		vector<int> mcsamples   = fMCBG; //NoQCDNoGJets;
		vector<int> datasamples = fHighPtData;
		//////////////////////////////////////////////////////////
		// Make kin plots
		for(size_t i = 0; i < gNKinVars; ++i){
			// Create plots
			bool intlabel = false;
			if(i == 2 || i == 10 || i == 11) intlabel = true;

			hvar_data[i]->SetXTitle(KinPlots::axis_label[i]);
			hvar_data[i]->SetLineWidth(3);
			hvar_data[i]->SetLineColor(kBlack);
			hvar_data[i]->SetMarkerStyle(8);
			hvar_data[i]->SetMarkerColor(kBlack);
			hvar_data[i]->SetMarkerSize(1.3);

			// Scale by luminosity
			for(size_t j = 0; j < gNSAMPLES; ++j){
				float lumiscale = fLumiNorm / fSamples[j]->getLumi();
				if(fSamples[j]->datamc == 0) continue;
				fSamples[j]->kinplots[s][HighPt].hvar[i]->Scale(lumiscale);
			}

			// Fill data histo
			for(size_t j = 0; j < datasamples.size(); ++j){
				Sample *S = fSamples[datasamples[j]];
				hvar_data[i]->Add(S->kinplots[s][HighPt].hvar[i]);
				hvar_data[i]->SetXTitle(KinPlots::axis_label[i]);
			}

			// // Scale to get equal integrals
			// float intscale(0.);
			// for(size_t j = 0; j < mcsamples.size();   ++j){
			// 	Sample *S = fSamples[mcsamples[j]];
			// 	intscale += S->kinplots[s][HighPt].hvar[i]->Integral();
			// }
			// intscale = hvar_data[i]->Integral() / intscale;
			// 
			// for(size_t j = 0; j < mcsamples.size();   ++j){
			// 	Sample *S = fSamples[mcsamples[j]];
			// 	S->kinplots[s][HighPt].hvar[i]->Scale(intscale);
			// }

			hvar_qcd [i]->SetFillColor(kYellow-4);
			hvar_db  [i]->SetFillColor(kSpring-9);
			hvar_ewk [i]->SetFillColor(kAzure+8);
			hvar_ttj [i]->SetFillColor(kAzure-5);
			hvar_rare[i]->SetFillColor(kGreen+1);
			// hvar_ttj [i]->SetFillColor(kAzure +1);
			// hvar_rare[i]->SetFillColor(kViolet+ 5);


			// Fill MC stacks
			for(size_t j = 0; j < mcsamples.size();   ++j){
			  Sample *S = fSamples[mcsamples[j]];
			  TString s_name = S->sname;
			  //			  if (s_name == "TTJets2") continue;
			  // sample type: QCD = 1 , Top = 2, EWK = 3 , Rare = 4 , DB = 5
 			  if ( S->getType() == 1)                       hvar_qcd [i]->Add( S->kinplots[s][HighPt].hvar[i] ); // ttZ
 			  if ( S->getType() == 5 && S->getProc() == 7 ) hvar_db  [i]->Add( S->kinplots[s][HighPt].hvar[i] ); // ttW
 			  if ( S->getType() == 2)                       hvar_ttj [i]->Add( S->kinplots[s][HighPt].hvar[i] ); // top + qcd
 			  if ( S->getType() == 3 || 
			       (S->getType() == 5 && S->getProc() != 7)) hvar_ewk [i]->Add( S->kinplots[s][HighPt].hvar[i] ); // single boson
			  if ( S->getType() == 4)                       hvar_rare[i]->Add( S->kinplots[s][HighPt].hvar[i] ); // di boson
 			}
 			hvar_mc_s[i]->Add(hvar_ttj[i]);
 			hvar_mc_s[i]->Add(hvar_ewk[i]);
 			hvar_mc_s[i]->Add(hvar_rare[i]);
 			hvar_mc_s[i]->Add(hvar_db[i]);
 			hvar_mc_s[i]->Add(hvar_qcd[i]);
 			hvar_mc_s[i]->Draw("goff");
			hvar_mc_s[i]->GetXaxis()->SetTitle(KinPlots::axis_label[i]);
			//			if(intlabel) for(size_t j = 1; j <= hvar_data[i]->GetNbinsX(); ++j)            hvar_data[i]->GetXaxis()->SetBinLabel(j, Form("%d", j-1));
			//			if(intlabel) for(size_t j = 1; j <= hvar_mc_s[i]->GetXaxis()->GetNbins(); ++j) hvar_mc_s[i]->GetXaxis()->SetBinLabel(j, Form("%d", j-1));
 			for(size_t j = 1; j <= hvar_data[i]->GetNbinsX(); ++j)            hvar_data[i]->GetXaxis()->SetBinLabel(j, "");
 			for(size_t j = 1; j <= hvar_mc_s[i]->GetXaxis()->GetNbins(); ++j) hvar_mc_s[i]->GetXaxis()->SetBinLabel(j, "");

			hvar_tot[i]->Add(hvar_qcd[i]);
			hvar_tot[i]->Add(hvar_db[i]);
			hvar_tot[i]->Add(hvar_ewk[i]);
			hvar_tot[i]->Add(hvar_rare[i]);
			hvar_tot[i]->Add(hvar_ttj[i]);

			float binwidth = hvar_data[i]->GetBinWidth(1);
			TString ytitle = Form("Events / %3.0f GeV", binwidth);
			if(intlabel) ytitle = "Events"; // Njets, Nbjets
			hvar_mc_s[i]->GetYaxis()->SetTitle(ytitle);
			if(i==2) hvar_mc_s[i]->GetYaxis()->SetTitleOffset(1.25);
			hvar_mc_s[i]->GetYaxis()->SetTitleSize(0.04);

			double max1 = hvar_mc_s[i]->GetMaximum();
			double max2 = hvar_data[i]->GetMaximum();
			double max = max1>max2?max1:max2;
			if(loglin){ // logarithmic
				hvar_mc_s[i]->SetMaximum(5.*max);
				hvar_data[i]->SetMaximum(5.*max);
				hvar_mc_s[i]->SetMinimum(0.5);
				hvar_data[i]->SetMinimum(0.5);
			}
			else{ // linear
				hvar_mc_s[i]->SetMaximum(1.5*max);
				hvar_data[i]->SetMaximum(1.5*max);
				hvar_mc_s[i]->SetMinimum(0.);
				hvar_data[i]->SetMinimum(0.);				
			}
			
			// Ratio plot:
			float border = 0.3;
			float scale = (1-border)/border;
						
			hvar_rat[i]->SetXTitle(hvar_data[i]->GetXaxis()->GetTitle());
			hvar_rat[i]->SetYTitle("");
			hvar_rat[i]->GetXaxis()->SetTitleSize(scale * 0.04);
			hvar_rat[i]->GetXaxis()->SetLabelSize(scale * hvar_data[i]->GetXaxis()->GetLabelSize());
			hvar_rat[i]->GetYaxis()->SetLabelSize(scale * hvar_data[i]->GetYaxis()->GetLabelSize());
			hvar_rat[i]->GetXaxis()->SetTickLength(scale * hvar_data[i]->GetXaxis()->GetTickLength());
			hvar_rat[i]->GetYaxis()->SetTickLength(hvar_data[i]->GetYaxis()->GetTickLength());
			if(intlabel){
				hvar_rat[i]->GetXaxis()->SetLabelSize(scale*0.06);
				hvar_rat[i]->GetXaxis()->SetLabelOffset(0.02);
				hvar_rat[i]->GetXaxis()->SetTitleOffset(1.20);
				for(size_t j = 1; j <= hvar_rat[i]->GetXaxis()->GetNbins(); ++j) hvar_rat[i]->GetXaxis()->SetBinLabel(j, Form("%d", j-1));
			}

			hvar_rat[i]->SetFillStyle(1001);
			hvar_rat[i]->SetLineWidth(1);
			hvar_rat[i]->SetFillColor(  kGray+1);
			hvar_rat[i]->SetLineColor(  kGray+1);
			hvar_rat[i]->SetMarkerColor(kGray+1);

			hvar_rat[i]->Divide(hvar_data[i], hvar_tot[i]);
			
			// Canvas
			TCanvas *c_temp = new TCanvas("C_" + KinPlots::var_name[i], KinPlots::var_name[i] + " in Data vs MC", 0, 0, 600, 600);
			c_temp->cd();

			TPad *p_plot  = new TPad("plotpad",  "Pad containing the plot", 0.00, border, 1.00, 1.00, 0, 0);
			p_plot->SetBottomMargin(0.015);
			p_plot->Draw();
			TPad *p_ratio = new TPad("ratiopad", "Pad containing the ratio", 0.00, 0.00, 1.00, border, 0, 0);
			p_ratio->SetTopMargin(0.025);
			p_ratio->SetBottomMargin(0.35);
			p_ratio->Draw();

			p_ratio->cd();
			hvar_rat[i]->GetYaxis()->SetNdivisions(505);
			// setPlottingRange(hvar_rat[i], 0.3);
			hvar_rat[i]->SetMaximum(1.99);
			hvar_rat[i]->SetMinimum(0.0);
			hvar_rat[i]->DrawCopy("E2 ");
			TLine *l3 = new TLine(hvar_data[i]->GetXaxis()->GetXmin(), 1.00, hvar_data[i]->GetXaxis()->GetXmax(), 1.00);
			l3->SetLineWidth(2);
			l3->SetLineStyle(7);
			l3->Draw();
			gPad->RedrawAxis();
			p_ratio->Draw();

			FakeRatios *FR = new FakeRatios();
			TGraphAsymmErrors* gr_obs = FR->getGraphPoissonErrors(hvar_data[i]);
			gr_obs->SetMarkerColor(kBlack);
			gr_obs->SetMarkerStyle(8);
			gr_obs->SetMarkerSize(1.2);
			gr_obs->SetLineWidth(2);
			gr_obs->SetLineColor(kBlack);
			gr_obs->SetFillColor(kBlack);
			
			p_plot->cd();
			if(loglin) p_plot->SetLogy(1);

			// TLegend *leg = new TLegend(0.15,0.50,0.40,0.88);
			// TLegend *leg = new TLegend(0.70,0.30,0.90,0.68);
			TLegend *leg = new TLegend(0.70,0.62,0.89,0.88);
// 			leg->AddEntry(hvar_data[i], "Data",      "p");
// 			leg->AddEntry(hvar_qcd[i],  "ttZ",       "f");
// 			leg->AddEntry(hvar_db[i],   "ttW",       "f");
// 			leg->AddEntry(hvar_rare[i], "Diboson",   "f");
// 			leg->AddEntry(hvar_ewk[i],  "Single boson", "f");
// 			leg->AddEntry(hvar_ttj[i],  "Top",       "f");
// 			leg->AddEntry(hvar_rat[i],  "Ratio",     "f");

			leg->AddEntry(hvar_data[i], "Data",      "p");
			leg->AddEntry(hvar_qcd[i],  "QCD",       "f");
			leg->AddEntry(hvar_db[i],   "WZ",       "f");
			leg->AddEntry(hvar_rare[i], "Rare SM",   "f");
			leg->AddEntry(hvar_ewk[i],  "W/Z/WW/ZZ", "f");
			leg->AddEntry(hvar_ttj[i],  "Top",       "f");
			leg->AddEntry(hvar_rat[i],  "Ratio",     "f");
	
			leg->SetFillStyle(0);
			leg->SetTextFont(42);
			leg->SetBorderSize(0);

			// if(loglin) gPad->SetLogy();
			hvar_mc_s[i]->Draw("hist");
			// hvar_data[i]->DrawCopy("PE X0 same");
			gr_obs->Draw("P same");
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

			p_plot->Draw();
			gPad->RedrawAxis();

			c_temp->Update();


			// Util::PrintNoEPS(c_temp, KinPlots::var_name[i], fOutputDir + fOutputSubDir, NULL);
			Util::PrintPDF(c_temp, KinPlots::var_name[i], fOutputDir + fOutputSubDir);
			Util::PrintPNG(c_temp, KinPlots::var_name[i], fOutputDir + fOutputSubDir);
						
			// Undo scaling
			for(size_t j = 0; j < gNSAMPLES; ++j){
				float lumiscale = fLumiNorm / fSamples[j]->getLumi();
				if(fSamples[j]->datamc == 0) continue;
				fSamples[j]->kinplots[s][HighPt].hvar[i]->Scale(1./lumiscale);
			}
			
			delete FR, gr_obs;
			delete c_temp;
			delete leg;
		}
	}	
}
void SSDLPlotter::makeMETvsHTPlot(vector<int> mmsamples, vector<int> eesamples, vector<int> emsamples, gHiLoSwitch hilo){
	if(readSigGraphs(fOutputFileName) != 0) return;
	TString hiloname[2] = {"p_{T}(l_{1}/l_{2}) > 20/20 GeV", "p_{T}(#mu/e) > 5/10 GeV"};

	fOutputSubDir = "KinematicPlots/HTvsMET/";
	char cmd[100];
	sprintf(cmd,"mkdir -p %s%s", fOutputDir.Data(), fOutputSubDir.Data());
	system(cmd);
	
	const float htmax = 1300.;
	const float metmax = 250.;

	Color_t col_mm = kBlack;
	Color_t col_ee = kRed;
	Color_t col_em = kBlue;

	// Create histograms
	TH2D *hmetvsht_da_mm = new TH2D("Data_HTvsMET_mm", "Data_HTvsMET_mm", 100, 0., htmax, 100, 0., metmax);
	TH2D *hmetvsht_da_ee = new TH2D("Data_HTvsMET_ee", "Data_HTvsMET_ee", 100, 0., htmax, 100, 0., metmax);
	TH2D *hmetvsht_da_em = new TH2D("Data_HTvsMET_em", "Data_HTvsMET_em", 100, 0., htmax, 100, 0., metmax);
	//TH2D *hmetvsht_da_mm = new TH2D("Data_HTvsMET_mm", "Data_HTvsMET_mm", KinPlots::nHTBins, KinPlots::HTmin, KinPlots::HTmax, KinPlots::nMETBins, KinPlots::METmin, KinPlots::METmax);
	//TH2D *hmetvsht_da_ee = new TH2D("Data_HTvsMET_ee", "Data_HTvsMET_ee", KinPlots::nHTBins, KinPlots::HTmin, KinPlots::HTmax, KinPlots::nMETBins, KinPlots::METmin, KinPlots::METmax);
	//TH2D *hmetvsht_da_em = new TH2D("Data_HTvsMET_em", "Data_HTvsMET_em", KinPlots::nHTBins, KinPlots::HTmin, KinPlots::HTmax, KinPlots::nMETBins, KinPlots::METmin, KinPlots::METmax);

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

// 	TGraph *gmetvsht_da_mm_ex1 = getSigEventGraph(Muon, 0., 450., 0., 30.);
// 	TGraph *gmetvsht_da_ee_ex1 = getSigEventGraph(Elec, 0., 450., 0., 30.);
// 	TGraph *gmetvsht_da_em_ex1 = getSigEventGraph(ElMu, 0., 450., 0., 30.);
	TGraph *gmetvsht_da_mm_ex2 = getSigEventGraph(Muon, 0., 80.,  0., 7000.);
	TGraph *gmetvsht_da_ee_ex2 = getSigEventGraph(Elec, 0., 80.,  0., 7000.);
	TGraph *gmetvsht_da_em_ex2 = getSigEventGraph(ElMu, 0., 80.,  0., 7000.);

// 	gmetvsht_da_mm_ex1->SetMarkerColor(16);
// 	gmetvsht_da_ee_ex1->SetMarkerColor(16);
// 	gmetvsht_da_em_ex1->SetMarkerColor(16);
// 	gmetvsht_da_mm_ex1->SetMarkerSize(1.5);
// 	gmetvsht_da_ee_ex1->SetMarkerSize(1.5);
// 	gmetvsht_da_em_ex1->SetMarkerSize(1.5);
	gmetvsht_da_mm_ex2->SetMarkerColor(16);
	gmetvsht_da_ee_ex2->SetMarkerColor(16);
	gmetvsht_da_em_ex2->SetMarkerColor(16);
	gmetvsht_da_mm_ex2->SetMarkerSize(1.7);
	gmetvsht_da_ee_ex2->SetMarkerSize(1.7);
	gmetvsht_da_em_ex2->SetMarkerSize(1.7);
	
	//////////////////////////////////////////////////////////
	// Custom added for 0 < MET < 50:
	// 	TGraph *gmetvsht_da_mm_lowmet = getSigEventGraph(Muon, 450., 7000., 0., 50.);
	// 	TGraph *gmetvsht_da_ee_lowmet = getSigEventGraph(Elec, 450., 7000., 0., 50.);
	// 	TGraph *gmetvsht_da_em_lowmet = getSigEventGraph(ElMu, 450., 7000., 0., 50.);
	
	// 	gmetvsht_da_mm_lowmet->SetMarkerColor(col_mm);
	// 	gmetvsht_da_mm_lowmet->SetMarkerStyle(8);
	// 	gmetvsht_da_mm_lowmet->SetMarkerSize(1.5);
	// 	gmetvsht_da_em_lowmet->SetMarkerColor(col_em);
	// 	gmetvsht_da_em_lowmet->SetMarkerStyle(23);
	// 	gmetvsht_da_em_lowmet->SetMarkerSize(1.5);
	// 	gmetvsht_da_ee_lowmet->SetMarkerColor(col_ee);
	// 	gmetvsht_da_ee_lowmet->SetMarkerStyle(21);
	// 	gmetvsht_da_ee_lowmet->SetMarkerSize(1.5);
	//////////////////////////////////////////////////////////

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

	TWbox *lowhtbox  = new TWbox(0., 0., lowerht, metmax, kBlack, 0, 0);
	TWbox *lowmetbox = new TWbox(lowerht, 0., htmax,  0., kBlack, 0, 0);
	lowhtbox ->SetFillColor(12);
	lowmetbox->SetFillColor(12);
	lowhtbox ->SetFillStyle(3005);
	lowmetbox->SetFillStyle(3005);
	TLine *boxborder1 = new TLine(lowerht,0.,lowerht,metmax);
	TLine *boxborder2 = new TLine(lowerht,0.,htmax,0.);
	TLine *boxborder3 = new TLine(htmax, 0.,htmax,0.);
	boxborder1->SetLineWidth(1);
	boxborder2->SetLineWidth(1);
	boxborder3->SetLineWidth(1);
	boxborder1->SetLineColor(14);
	boxborder2->SetLineColor(14);
	boxborder3->SetLineColor(14);

	//HT80   MET0
	//HT200  MET50 
	//HT320  MET0
	//HT200  MET50
	//HT200  MET120
	//HT320  MET50
	//HT320  MET120

	TLine *sig1x = new TLine(lowerht, 30., lowerht,  metmax);   //HT80  MET30
	TLine *sig1y = new TLine(lowerht, 30.,   htmax,     30.); 
	TLine *sig2x = new TLine(200.,    50.,    200.,  metmax);   //HT200 MET50
	TLine *sig2y = new TLine(200.,    50.,   htmax,     50.);           
	TLine *sig3x = new TLine(320.,     0.,    320.,  metmax);   //HT320 MET0
	TLine *sig4x = new TLine(200.,   120.,    200.,  metmax);   //HT200 MET120
	TLine *sig4y = new TLine(200.,   120.,   htmax,    120.);
	TLine *sig5x = new TLine(320.,    50.,    320.,  metmax);   //HT320 MET50
	TLine *sig5y = new TLine(320.,    50.,   htmax,     50.);
	TLine *sig6x = new TLine(320.,   120.,    320.,  metmax);   //HT320 MET120
	TLine *sig6y = new TLine(320.,   120.,   htmax,    120.);

	sig1x->SetLineWidth(2);
	sig1y->SetLineWidth(2);
	sig2x->SetLineWidth(2);
	sig2y->SetLineWidth(2);
	sig3x->SetLineWidth(1);
	sig4x->SetLineWidth(3);
	sig4y->SetLineWidth(3);
	sig5x->SetLineWidth(3);
	sig5y->SetLineWidth(3);
	sig6x->SetLineWidth(3);
	sig6y->SetLineWidth(3);

	sig1x->SetLineStyle(3);
	sig1y->SetLineStyle(3);
	sig2x->SetLineStyle(2);
	sig2y->SetLineStyle(2);
	sig3x->SetLineStyle(1);
	// sig4x->SetLineStyle(1);
	// sig4y->SetLineStyle(1);
	sig5x->SetLineStyle(4);
	sig5y->SetLineStyle(4);

	sig6x->SetLineStyle(5);
	sig6y->SetLineStyle(5);


	// float legymax = hilo==HighPt?0.54:0.50;
	// TLegend *regleg = new TLegend(0.70,0.37,0.88,0.54);
	// TLegend *regleg = new TLegend(0.70,0.45,0.88,0.62);
	// TLegend *regleg = new TLegend(0.70,0.29,0.88,0.46);
	TLegend *regleg = new TLegend(0.70,0.27,0.88,0.47);
	regleg->AddEntry(sig1x, "Search Region 1","l");
	regleg->AddEntry(sig2x, "Search Region 2","l");
	regleg->AddEntry(sig3x, "Search Region 3","l");
	regleg->AddEntry(sig4x, "Search Region 4","l");
	regleg->AddEntry(sig5x, "Search Region 5","l");
	regleg->AddEntry(sig6x, "Search Region 6","l");
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
	boxborder3->Draw();

	if(hilo != LowPt) sig1x->Draw();
	if(hilo != LowPt) sig1y->Draw();
	sig2x->Draw();
	sig2y->Draw();
	sig3x->Draw();
	sig4x->Draw();
	sig4y->Draw();
	sig5x->Draw();
	sig5y->Draw();
	sig6x->Draw();
	sig6y->Draw();

	// Graphs
	gmetvsht_da_ee->Draw("P");
	gmetvsht_da_em->Draw("P");
	gmetvsht_da_mm->Draw("P");

	gmetvsht_da_mm_ex2->Draw("P");
	gmetvsht_da_ee_ex2->Draw("P");
	gmetvsht_da_em_ex2->Draw("P");
	
// 	gmetvsht_da_mm_lowmet->Draw("P");
// 	gmetvsht_da_ee_lowmet->Draw("P");
// 	gmetvsht_da_em_lowmet->Draw("P");

	leg->Draw();
	regleg->Draw();
	drawTopLine();
	// TPaveText *pave = new TPaveText(0.16, 0.83, 0.53, 0.88, "NDC");
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

	//////////////////////////////////////////////////////////
	// Custom added for 0 < MET < 50:
	TGraph *gmetvsht_da_mm_lowmet = getSigEventGraph(Muon, 450., 7000., 0., 50.);
	TGraph *gmetvsht_da_ee_lowmet = getSigEventGraph(Elec, 450., 7000., 0., 50.);
	TGraph *gmetvsht_da_em_lowmet = getSigEventGraph(ElMu, 450., 7000., 0., 50.);

	gmetvsht_da_mm_lowmet->SetMarkerColor(col_mm);
	gmetvsht_da_mm_lowmet->SetMarkerStyle(8);
	gmetvsht_da_mm_lowmet->SetMarkerSize(1.5);
	gmetvsht_da_em_lowmet->SetMarkerColor(col_em);
	gmetvsht_da_em_lowmet->SetMarkerStyle(23);
	gmetvsht_da_em_lowmet->SetMarkerSize(1.5);
	gmetvsht_da_ee_lowmet->SetMarkerColor(col_ee);
	gmetvsht_da_ee_lowmet->SetMarkerStyle(21);
	gmetvsht_da_ee_lowmet->SetMarkerSize(1.5);
	//////////////////////////////////////////////////////////

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
	// const float lowerht = 350.; // tau only
	const float minmet = 30.;

	TWbox *lowhtbox  = new TWbox(0.,      0., lowerht, metmax, kBlack, 0, 0);
	TWbox *lowmetbox = new TWbox(lowerht, 0., 450.,   minmet, kBlack, 0, 0);
	lowhtbox ->SetFillColor(12);
	lowmetbox->SetFillColor(12);
	lowhtbox ->SetFillStyle(3005);
	lowmetbox->SetFillStyle(3005);
	TLine *boxborder1 = new TLine(lowerht, minmet, lowerht, metmax);
	TLine *boxborder2 = new TLine(lowerht, minmet, 450.,   minmet);
	TLine *boxborder3 = new TLine(450., 0., 450., 30.);
	boxborder1->SetLineWidth(1);
	boxborder2->SetLineWidth(1);
	boxborder3->SetLineWidth(1);
	boxborder1->SetLineColor(14);
	boxborder2->SetLineColor(14);
	boxborder3->SetLineColor(14);

	TLine *sig1x = new TLine(lowerht, 120., lowerht, metmax); // met 120 ht 80
	TLine *sig1y = new TLine(lowerht, 120., htmax,   120.); 
	TLine *sig2x = new TLine(200., 120., 200.,  metmax);      // met 120 ht 200
	TLine *sig2y = new TLine(200., 120., htmax, 120.);
	TLine *sig3x = new TLine(450.,  50., 450.,  metmax);      // met 50  ht 450
	TLine *sig3y = new TLine(450.,  50., htmax,  50.);
	TLine *sig4x = new TLine(450., 120., 450.,  metmax);      // met 120 ht 450
	TLine *sig4y = new TLine(450., 120., htmax, 120.);
	TLine *sig5x = new TLine(450.,   0., 450.,  metmax);      // met 0 ht 450

	sig1x->SetLineWidth(2);
	sig1y->SetLineWidth(2);
	sig2x->SetLineWidth(2);
	sig2y->SetLineWidth(2);
	sig3x->SetLineWidth(1);
	sig3y->SetLineWidth(1);
	sig4x->SetLineWidth(3);
	sig4y->SetLineWidth(3);
	sig5x->SetLineWidth(2);

	sig1x->SetLineStyle(3);
	sig1y->SetLineStyle(3);
	sig2x->SetLineStyle(2);
	sig2y->SetLineStyle(2);
	sig3x->SetLineStyle(1);
	sig3y->SetLineStyle(1);
	// sig4x->SetLineStyle(1);
	// sig4y->SetLineStyle(1);
	sig5x->SetLineStyle(4);

	// TLegend *regleg = new TLegend(0.70,0.47,0.88,0.6);
	// TLegend *regleg = new TLegend(0.67,0.51,0.87,0.68);
	// TLegend *regleg = new TLegend(0.76,0.29,0.93,0.46);
	TLegend *regleg = new TLegend(0.76,0.26,0.93,0.48);
	regleg->AddEntry(sig1x, "Region 1","l");
	regleg->AddEntry(sig2x, "Region 2","l");
	regleg->AddEntry(sig3x, "Region 3","l");
	regleg->AddEntry(sig4x, "Region 4","l");
	regleg->AddEntry(sig5x, "Region 5","l");
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
	boxborder3->Draw();

	if(hilo != LowPt) sig1x->Draw();
	if(hilo != LowPt) sig1y->Draw();
	sig2x->Draw();
	sig2y->Draw();
	sig3x->Draw();
	sig3y->Draw();
	sig4x->Draw();
	sig4y->Draw();
	sig5x->Draw();

	// Graphs
	gmetvsht_da_ee->Draw("P");
	gmetvsht_da_em->Draw("P");
	gmetvsht_da_mm->Draw("P");
	
	gmetvsht_da_ee_lowpt->Draw("P");
	gmetvsht_da_em_lowpt->Draw("P");
	gmetvsht_da_mm_lowpt->Draw("P");
	
	gmetvsht_da_mt->Draw("P");
	gmetvsht_da_et->Draw("P");
	
	gmetvsht_da_ee_lowmet->Draw("P");
	gmetvsht_da_em_lowmet->Draw("P");
	gmetvsht_da_mm_lowmet->Draw("P");
	
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
void SSDLPlotter::makeMETvsHTPlot0HT(){
	fOutputSubDir = "NoHT";
	char cmd[100];
	sprintf(cmd,"mkdir -p %s%s", fOutputDir.Data(), fOutputSubDir.Data());
	system(cmd);
	
	const float htmax = 1500.;
	const float metmax = 550.;
	const float metmin = 50.;
	const float njetsmax = 5.25;
	
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
	TH2D *hmetvsht_da_mm = new TH2D("Data_HTvsMET_mm", "Data_HTvsMET_mm", 100, 0., htmax, 100, metmin, metmax);
	TH2D *hmetvsht_da_ee = new TH2D("Data_HTvsMET_ee", "Data_HTvsMET_ee", 100, 0., htmax, 100, metmin, metmax);
	TH2D *hmetvsht_da_em = new TH2D("Data_HTvsMET_em", "Data_HTvsMET_em", 100, 0., htmax, 100, metmin, metmax);

	TH2D *hmetvsnjets_da_mm = new TH2D("Data_NjetsvsMET_mm", "Data_NjetsvsMET_mm", njetsmax, 0., njetsmax, 100, metmin, metmax);
	TH2D *hmetvsnjets_da_ee = new TH2D("Data_NjetsvsMET_ee", "Data_NjetsvsMET_ee", njetsmax, 0., njetsmax, 100, metmin, metmax);
	TH2D *hmetvsnjets_da_em = new TH2D("Data_NjetsvsMET_em", "Data_NjetsvsMET_em", njetsmax, 0., njetsmax, 100, metmin, metmax);

	TH2D *hmetvsnjetswithb_da_mm = new TH2D("Data_NjetsvsMETwithb_mm", "Data_NjetsvsMETwithb_mm", njetsmax, 0., njetsmax, 100, metmin, metmax);
	TH2D *hmetvsnjetswithb_da_ee = new TH2D("Data_NjetsvsMETwithb_ee", "Data_NjetsvsMETwithb_ee", njetsmax, 0., njetsmax, 100, metmin, metmax);
	TH2D *hmetvsnjetswithb_da_em = new TH2D("Data_NjetsvsMETwithb_em", "Data_NjetsvsMETwithb_em", njetsmax, 0., njetsmax, 100, metmin, metmax);

	TH2D *hmetvsnjetsnob_da_mm = new TH2D("Data_NjetsvsMETnob_mm", "Data_NjetsvsMETnob_mm", njetsmax, 0., njetsmax, 100, metmin, metmax);
	TH2D *hmetvsnjetsnob_da_ee = new TH2D("Data_NjetsvsMETnob_ee", "Data_NjetsvsMETnob_ee", njetsmax, 0., njetsmax, 100, metmin, metmax);
	TH2D *hmetvsnjetsnob_da_em = new TH2D("Data_NjetsvsMETnob_em", "Data_NjetsvsMETnob_em", njetsmax, 0., njetsmax, 100, metmin, metmax);

	TH2D *hmetvsnjetsboth_da_mm = new TH2D("Data_NjetsvsMETboth_mm", "Data_NjetsvsMETboth_mm", njetsmax, 0., njetsmax, 100, metmin, metmax);
	TH2D *hmetvsnjetsboth_da_ee = new TH2D("Data_NjetsvsMETboth_ee", "Data_NjetsvsMETboth_ee", njetsmax, 0., njetsmax, 100, metmin, metmax);
	TH2D *hmetvsnjetsboth_da_em = new TH2D("Data_NjetsvsMETboth_em", "Data_NjetsvsMETboth_em", njetsmax, 0., njetsmax, 100, metmin, metmax);

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

	hmetvsht_da_mm->SetXTitle(KinPlots::axis_label[0]);
	hmetvsht_da_mm->SetYTitle(KinPlots::axis_label[1]);
	hmetvsht_da_ee->SetXTitle(KinPlots::axis_label[0]);
	hmetvsht_da_ee->SetYTitle(KinPlots::axis_label[1]);
	hmetvsht_da_em->SetXTitle(KinPlots::axis_label[0]);
	hmetvsht_da_em->SetYTitle(KinPlots::axis_label[1]);

	//////////////////////////////////////////////////////////
	// Make MET vs njets plot:
	hmetvsnjets_da_mm->SetMarkerStyle(8);
	hmetvsnjets_da_mm->SetMarkerColor(col_mm);
	hmetvsnjets_da_mm->SetMarkerSize(1.5);
	hmetvsnjets_da_mm->GetYaxis()->SetTitleOffset(1.4);
	hmetvsnjets_da_mm->GetXaxis()->SetNdivisions(606);

	hmetvsnjets_da_ee->SetMarkerStyle(21);
	hmetvsnjets_da_ee->SetMarkerColor(col_ee);
	hmetvsnjets_da_ee->SetMarkerSize(1.4);
	hmetvsnjets_da_ee->GetYaxis()->SetTitleOffset(1.4);
	hmetvsnjets_da_ee->GetXaxis()->SetNdivisions(500);

	hmetvsnjets_da_em->SetMarkerStyle(23);
	hmetvsnjets_da_em->SetMarkerColor(col_em);
	hmetvsnjets_da_em->SetMarkerSize(1.7  );
	hmetvsnjets_da_em->GetYaxis()->SetTitleOffset(1.4);
	hmetvsnjets_da_em->GetXaxis()->SetNdivisions(500);

	hmetvsnjets_da_mm->SetXTitle(KinPlots::axis_label[2]);
	hmetvsnjets_da_mm->SetYTitle(KinPlots::axis_label[1]);
	hmetvsnjets_da_ee->SetXTitle(KinPlots::axis_label[2]);
	hmetvsnjets_da_ee->SetYTitle(KinPlots::axis_label[1]);
	hmetvsnjets_da_em->SetXTitle(KinPlots::axis_label[2]);
	hmetvsnjets_da_em->SetYTitle(KinPlots::axis_label[1]);

	//////////////////////////////////////////////////////////
	// Make MET vs njets plot with b-tags
	hmetvsnjetswithb_da_mm->SetMarkerStyle(8);
	hmetvsnjetswithb_da_mm->SetMarkerColor(col_mm);
	hmetvsnjetswithb_da_mm->SetMarkerSize(1.5);
	hmetvsnjetswithb_da_mm->GetYaxis()->SetTitleOffset(1.4);
	hmetvsnjetswithb_da_mm->GetXaxis()->SetNdivisions(606);

	hmetvsnjetswithb_da_ee->SetMarkerStyle(21);
	hmetvsnjetswithb_da_ee->SetMarkerColor(col_ee);
	hmetvsnjetswithb_da_ee->SetMarkerSize(1.4);
	hmetvsnjetswithb_da_ee->GetYaxis()->SetTitleOffset(1.4);
	hmetvsnjetswithb_da_ee->GetXaxis()->SetNdivisions(500);

	hmetvsnjetswithb_da_em->SetMarkerStyle(23);
	hmetvsnjetswithb_da_em->SetMarkerColor(col_em);
	hmetvsnjetswithb_da_em->SetMarkerSize(1.7  );
	hmetvsnjetswithb_da_em->GetYaxis()->SetTitleOffset(1.4);
	hmetvsnjetswithb_da_em->GetXaxis()->SetNdivisions(500);

	hmetvsnjetswithb_da_mm->SetXTitle(KinPlots::axis_label[2]);
	hmetvsnjetswithb_da_mm->SetYTitle(KinPlots::axis_label[1]);
	hmetvsnjetswithb_da_ee->SetXTitle(KinPlots::axis_label[2]);
	hmetvsnjetswithb_da_ee->SetYTitle(KinPlots::axis_label[1]);
	hmetvsnjetswithb_da_em->SetXTitle(KinPlots::axis_label[2]);
	hmetvsnjetswithb_da_em->SetYTitle(KinPlots::axis_label[1]);


	//////////////////////////////////////////////////////////
	// Make MET vs njets plot no b-tags
	hmetvsnjetsnob_da_mm->SetMarkerStyle(8);
	hmetvsnjetsnob_da_mm->SetMarkerColor(col_mm);
	hmetvsnjetsnob_da_mm->SetMarkerSize(1.5);
	hmetvsnjetsnob_da_mm->GetYaxis()->SetTitleOffset(1.4);
	hmetvsnjetsnob_da_mm->GetXaxis()->SetNdivisions(606);

	hmetvsnjetsnob_da_ee->SetMarkerStyle(21);
	hmetvsnjetsnob_da_ee->SetMarkerColor(col_ee);
	hmetvsnjetsnob_da_ee->SetMarkerSize(1.4);
	hmetvsnjetsnob_da_ee->GetYaxis()->SetTitleOffset(1.4);
	hmetvsnjetsnob_da_ee->GetXaxis()->SetNdivisions(500);

	hmetvsnjetsnob_da_em->SetMarkerStyle(23);
	hmetvsnjetsnob_da_em->SetMarkerColor(col_em);
	hmetvsnjetsnob_da_em->SetMarkerSize(1.7  );
	hmetvsnjetsnob_da_em->GetYaxis()->SetTitleOffset(1.4);
	hmetvsnjetsnob_da_em->GetXaxis()->SetNdivisions(500);

	hmetvsnjetsnob_da_mm->SetXTitle(KinPlots::axis_label[2]);
	hmetvsnjetsnob_da_mm->SetYTitle(KinPlots::axis_label[1]);
	hmetvsnjetsnob_da_ee->SetXTitle(KinPlots::axis_label[2]);
	hmetvsnjetsnob_da_ee->SetYTitle(KinPlots::axis_label[1]);
	hmetvsnjetsnob_da_em->SetXTitle(KinPlots::axis_label[2]);
	hmetvsnjetsnob_da_em->SetYTitle(KinPlots::axis_label[1]);

	//////////////////////////////////////////////////////////
	// Make MET vs njets plot both
	hmetvsnjetsboth_da_mm->SetMarkerStyle(8);
	hmetvsnjetsboth_da_mm->SetMarkerColor(col_mm);
	hmetvsnjetsboth_da_mm->SetMarkerSize(1.5);
	hmetvsnjetsboth_da_mm->GetYaxis()->SetTitleOffset(1.4);
	hmetvsnjetsboth_da_mm->GetXaxis()->SetNdivisions(606);

	hmetvsnjetsboth_da_ee->SetMarkerStyle(21);
	hmetvsnjetsboth_da_ee->SetMarkerColor(col_ee);
	hmetvsnjetsboth_da_ee->SetMarkerSize(1.4);
	hmetvsnjetsboth_da_ee->GetYaxis()->SetTitleOffset(1.4);
	hmetvsnjetsboth_da_ee->GetXaxis()->SetNdivisions(500);

	hmetvsnjetsboth_da_em->SetMarkerStyle(23);
	hmetvsnjetsboth_da_em->SetMarkerColor(col_em);
	hmetvsnjetsboth_da_em->SetMarkerSize(1.7  );
	hmetvsnjetsboth_da_em->GetYaxis()->SetTitleOffset(1.4);
	hmetvsnjetsboth_da_em->GetXaxis()->SetNdivisions(500);

	hmetvsnjetsboth_da_mm->SetXTitle(KinPlots::axis_label[2]);
	hmetvsnjetsboth_da_mm->SetYTitle(KinPlots::axis_label[1]);
	hmetvsnjetsboth_da_ee->SetXTitle(KinPlots::axis_label[2]);
	hmetvsnjetsboth_da_ee->SetYTitle(KinPlots::axis_label[1]);
	hmetvsnjetsboth_da_em->SetXTitle(KinPlots::axis_label[2]);
	hmetvsnjetsboth_da_em->SetYTitle(KinPlots::axis_label[1]);

	TLegend *leg = new TLegend(0.80,0.70,0.95,0.88);
	// TLegend *leg = new TLegend(0.67,0.70,0.82,0.88);
	leg->AddEntry(hmetvsht_da_mm, "#mu#mu","p");
	leg->AddEntry(hmetvsht_da_ee, "ee","p");
	leg->AddEntry(hmetvsht_da_em, "e#mu","p");
	leg->SetFillStyle(0);
	leg->SetTextFont(42);
	leg->SetTextSize(0.05);
	leg->SetBorderSize(0);

///////////// FROM SIGEVENTS TREE

	// TGraphs:
	TGraph *gmetvsht_da_mm = getSigEventGraph(Muon, 0., 7000., 120., 7000., 0);
	TGraph *gmetvsht_da_ee = getSigEventGraph(Elec, 0., 7000., 120., 7000., 0);
	TGraph *gmetvsht_da_em = getSigEventGraph(ElMu, 0., 7000., 120., 7000., 0);
	// TGraph *gmetvsht_da_mm = getSigEventGraph(Muon, HT0MET120);
	// TGraph *gmetvsht_da_ee = getSigEventGraph(Elec, HT0MET120);
	// TGraph *gmetvsht_da_em = getSigEventGraph(ElMu, HT0MET120);
	TGraph *gmetvsnjets_da_mm = getSigEventGraph(Muon, 0., 7000., 120., 7000., 1);
	TGraph *gmetvsnjets_da_ee = getSigEventGraph(Elec, 0., 7000., 120., 7000., 1);
	TGraph *gmetvsnjets_da_em = getSigEventGraph(ElMu, 0., 7000., 120., 7000., 1);
	TGraph *gmetvsnjetswithb_da_mm = getSigEventGraph(Muon, 0., 7000., 120., 7000., 2);
	TGraph *gmetvsnjetswithb_da_ee = getSigEventGraph(Elec, 0., 7000., 120., 7000., 2);
	TGraph *gmetvsnjetswithb_da_em = getSigEventGraph(ElMu, 0., 7000., 120., 7000., 2);
	TGraph *gmetvsnjetsnob_da_mm = getSigEventGraph(Muon, 0., 7000., 120., 7000., 3);
	TGraph *gmetvsnjetsnob_da_ee = getSigEventGraph(Elec, 0., 7000., 120., 7000., 3);
	TGraph *gmetvsnjetsnob_da_em = getSigEventGraph(ElMu, 0., 7000., 120., 7000., 3);
	TGraph *gmetvsnjetsboth_da_mm = getSigEventGraph(Muon, 0., 7000., 120., 7000., 4);
	TGraph *gmetvsnjetsboth_da_ee = getSigEventGraph(Elec, 0., 7000., 120., 7000., 4);
	TGraph *gmetvsnjetsboth_da_em = getSigEventGraph(ElMu, 0., 7000., 120., 7000., 4);


	TGraph *gmetvsht_da_mm_ex = getSigEventGraph(Muon, 0., 7000., 0., 120., 0);
	TGraph *gmetvsht_da_ee_ex = getSigEventGraph(Elec, 0., 7000., 0., 120., 0);
	TGraph *gmetvsht_da_em_ex = getSigEventGraph(ElMu, 0., 7000., 0., 120., 0);
	TGraph *gmetvsnjets_da_mm_ex = getSigEventGraph(Muon, 0., 7000., 0., 120., 1);
	TGraph *gmetvsnjets_da_ee_ex = getSigEventGraph(Elec, 0., 7000., 0., 120., 1);
	TGraph *gmetvsnjets_da_em_ex = getSigEventGraph(ElMu, 0., 7000., 0., 120., 1);
	TGraph *gmetvsnjetswithb_da_mm_ex = getSigEventGraph(Muon, 0., 7000., 0., 120., 2);
	TGraph *gmetvsnjetswithb_da_ee_ex = getSigEventGraph(Elec, 0., 7000., 0., 120., 2);
	TGraph *gmetvsnjetswithb_da_em_ex = getSigEventGraph(ElMu, 0., 7000., 0., 120., 2);
	TGraph *gmetvsnjetsnob_da_mm_ex = getSigEventGraph(Muon, 0., 7000., 0., 120., 3);
	TGraph *gmetvsnjetsnob_da_ee_ex = getSigEventGraph(Elec, 0., 7000., 0., 120., 3);
	TGraph *gmetvsnjetsnob_da_em_ex = getSigEventGraph(ElMu, 0., 7000., 0., 120., 3);
	TGraph *gmetvsnjetsboth_da_mm_ex = getSigEventGraph(Muon, 0., 7000., 0., 120., 4);
	TGraph *gmetvsnjetsboth_da_ee_ex = getSigEventGraph(Elec, 0., 7000., 0., 120., 4);
	TGraph *gmetvsnjetsboth_da_em_ex = getSigEventGraph(ElMu, 0., 7000., 0., 120., 4);

	gmetvsht_da_mm_ex->SetMarkerColor(16);
	gmetvsht_da_ee_ex->SetMarkerColor(16);
	gmetvsht_da_em_ex->SetMarkerColor(16);
	gmetvsht_da_mm_ex->SetMarkerSize(1.2);
	gmetvsht_da_ee_ex->SetMarkerSize(1.2);
	gmetvsht_da_em_ex->SetMarkerSize(1.2);
	gmetvsnjets_da_mm_ex->SetMarkerColor(16);
	gmetvsnjets_da_ee_ex->SetMarkerColor(16);
	gmetvsnjets_da_em_ex->SetMarkerColor(16);
	gmetvsnjets_da_mm_ex->SetMarkerSize(1.2);
	gmetvsnjets_da_ee_ex->SetMarkerSize(1.2);
	gmetvsnjets_da_em_ex->SetMarkerSize(1.2);
	gmetvsnjetswithb_da_mm_ex->SetMarkerColor(16);
	gmetvsnjetswithb_da_ee_ex->SetMarkerColor(16);
	gmetvsnjetswithb_da_em_ex->SetMarkerColor(16);
	gmetvsnjetswithb_da_mm_ex->SetMarkerSize(1.2);
	gmetvsnjetswithb_da_ee_ex->SetMarkerSize(1.2);
	gmetvsnjetswithb_da_em_ex->SetMarkerSize(1.2);
	gmetvsnjetsnob_da_mm_ex->SetMarkerColor(16);
	gmetvsnjetsnob_da_ee_ex->SetMarkerColor(16);
	gmetvsnjetsnob_da_em_ex->SetMarkerColor(16);
	gmetvsnjetsnob_da_mm_ex->SetMarkerSize(1.2);
	gmetvsnjetsnob_da_ee_ex->SetMarkerSize(1.2);
	gmetvsnjetsnob_da_em_ex->SetMarkerSize(1.2);
	gmetvsnjetsboth_da_mm_ex->SetMarkerColor(16);
	gmetvsnjetsboth_da_ee_ex->SetMarkerColor(16);
	gmetvsnjetsboth_da_em_ex->SetMarkerColor(16);
	gmetvsnjetsboth_da_mm_ex->SetMarkerSize(1.2);
	gmetvsnjetsboth_da_ee_ex->SetMarkerSize(1.2);
	gmetvsnjetsboth_da_em_ex->SetMarkerSize(1.2);

///////////////////////


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


	// 0 HT tau numbers
	const int nmtev = 19;
	float a_mt_ht [nmtev] = {338.843, 454.149, 96.3014, 549.702, 197.151, 0, 487.826, 228.581, 0, 144.391, 149.425, 218.886, 122.39, 57.3992, 398.942, 334.117, 864.134, 233.842, 215.262};
	float a_mt_met[nmtev] = {151.176, 144.284, 124.364, 143.927, 127.371, 204.468, 161.936, 145.025, 126.596, 128.989, 199.845, 162.472, 142.881, 149.945, 142.304, 163.266, 134.839, 127.221, 122.126};

	const int netev = 9;
	float a_et_ht [netev] = {0, 137.963, 250.429, 0, 312.128, 252.608, 219.162, 48.0821, 134.521};
	float a_et_met[netev] = {129.379, 131.362, 139.101, 228.448, 133.339, 176.02 , 124.487, 142.518, 135.69};

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
	const float lowerht = 0.;
	// const float lowerht = 350.; // tau only
	const float minmet = 120.;

	TWbox *lowhtbox  = new TWbox(0.,      0., lowerht, metmax, kBlack, 0, 0);
	TWbox *lowmetbox = new TWbox(lowerht, metmin, htmax, minmet, kBlack, 0, 0);
	lowhtbox ->SetFillColor(12);
	lowmetbox->SetFillColor(12);
	lowhtbox ->SetFillStyle(3005);
	lowmetbox->SetFillStyle(3005);
	TLine *boxborder1 = new TLine(lowerht, minmet, lowerht, metmax);
	TLine *boxborder2 = new TLine(lowerht, minmet, htmax,   minmet);
	boxborder1->SetLineWidth(1);
	boxborder2->SetLineWidth(1);
	boxborder1->SetLineColor(14);
	boxborder2->SetLineColor(14);

	TWbox *lowhtbox1  = new TWbox(0.,      0., 0., metmax, kBlack, 0, 0);
	TWbox *lowmetbox1 = new TWbox(0., metmin, njetsmax, minmet, kBlack, 0, 0);
	lowhtbox1 ->SetFillColor(12);
	lowmetbox1->SetFillColor(12);
	lowhtbox1 ->SetFillStyle(3005);
	lowmetbox1->SetFillStyle(3005);
	TLine *boxborder11 = new TLine(0., minmet, 0., metmax);
	TLine *boxborder21 = new TLine(0., minmet, njetsmax,   minmet);
	boxborder11->SetLineWidth(1);
	boxborder21->SetLineWidth(1);
	boxborder11->SetLineColor(14);
	boxborder21->SetLineColor(14);

	TLine *sig1 = new TLine(lowerht, 200., htmax, 200.); // met 120 ht 80

	sig1->SetLineWidth(2);
	sig1->SetLineStyle(1);

	// TLegend *regleg = new TLegend(0.70,0.47,0.88,0.6);
	// TLegend *regleg = new TLegend(0.67,0.51,0.87,0.68);
	TLegend *regleg = new TLegend(0.67,0.26,0.93,0.46);
	regleg->AddEntry(sig1, "Signal Region","");
	regleg->SetFillStyle(0);
	regleg->SetTextFont(42);
	regleg->SetTextSize(0.03);
	regleg->SetBorderSize(0);
	
	
	TLegend *regleg1 = new TLegend(0.67,0.18,0.93,0.46);
	regleg1->AddEntry(sig1, "Control Region","");
	regleg1->SetFillStyle(0);
	regleg1->SetTextFont(42);
	regleg1->SetTextSize(0.03);
	regleg1->SetBorderSize(0);
	
	
	TLegend *regleg2 = new TLegend(0.53,0.15,0.93,0.20);
	regleg2->AddEntry(sig1, "Region Not Considered","");
	regleg2->SetFillStyle(0);
	regleg2->SetTextFont(42);
	regleg2->SetTextSize(0.03);
	regleg2->SetBorderSize(0);
	

	// HT vs. MET plot first
	TCanvas *c_temp = new TCanvas("C_HTvsMET", "HT vs MET in Data vs MC", 0, 0, 600, 600);
	c_temp->cd();
	c_temp->SetRightMargin(0.05);
	c_temp->SetLeftMargin(0.13);

	hmetvsht_da_mm->DrawCopy("axis");

	gmetvsht_da_mm_ex->Draw("P");
	gmetvsht_da_ee_ex->Draw("P");
	gmetvsht_da_em_ex->Draw("P");

	lowhtbox ->Draw();
	//	lowmetbox->Draw();
	boxborder1->Draw();
	boxborder2->Draw();

	sig1->Draw();

	// Graphs
	gmetvsht_da_ee->Draw("P");
	gmetvsht_da_em->Draw("P");
	gmetvsht_da_mm->Draw("P");
	
	// gmetvsht_da_mt->Draw("P");
	// gmetvsht_da_et->Draw("P");
	
	leg->Draw();
	// leg2->Draw();
	regleg->Draw();
	regleg1->Draw();
	regleg2->Draw();

	drawTopLine(0.53, 0.8);
	// TPaveText *pave = new TPaveText(0.16, 0.83, 0.55, 0.88, "NDC");
	// pave->SetFillColor(0);
	// pave->SetFillStyle(1001);
	// pave->SetBorderSize(0);
	// pave->SetMargin(0.05);
	// pave->SetTextFont(42);
	// pave->SetTextSize(0.04);
	// pave->SetTextAlign(12);
	// pave->Draw();
	gPad->RedrawAxis();

	// Util::PrintNoEPS(c_temp, "HTvsMET_" + gHiLoLabel[hilo], fOutputDir + fOutputSubDir, NULL);
	Util::PrintPDF(c_temp, "HTvsMET_NoHT", fOutputDir + fOutputSubDir);

	TLatex *lat = new TLatex();
	lat->SetNDC(kTRUE);
	lat->SetTextColor(kBlack);
	lat->SetTextSize(0.04);

	// make lines for every jet multiplicity
	gStyle->SetLineStyleString(11, "30 60");
	gStyle->SetLineStyleString(12, "20 60");
	TLine *njets1 = new TLine(0.5, metmin, 0.5, metmax);
	njets1->SetLineWidth(2);
	njets1->SetLineStyle(12);
	TLine *njets2 = new TLine(1.5, metmin, 1.5, 400.);
	njets2->SetLineWidth(2);
	njets2->SetLineStyle(12);
	TLine *njets3 = new TLine(2.5, metmin, 2.5, metmax);
	njets3->SetLineWidth(2);
	njets3->SetLineStyle(12);
	TLine *njets4 = new TLine(3.5, metmin, 3.5, metmax);
	njets4->SetLineWidth(2);
	njets4->SetLineStyle(12);
	TLine *njets5 = new TLine(4.5, metmin, 4.5, 400.);
	njets5->SetLineWidth(2);
	njets5->SetLineStyle(12);

	// NJETS vs. MET plot second
	TCanvas *c_temp1 = new TCanvas("C_NjetsvsMET", "Njets vs MET in Data vs MC", 0, 0, 600, 600);
	c_temp1->cd();
	c_temp1->SetRightMargin(0.05);
	c_temp1->SetLeftMargin(0.13);

	hmetvsnjets_da_mm->DrawCopy("axis");

	gmetvsnjets_da_mm_ex->Draw("P");
	gmetvsnjets_da_ee_ex->Draw("P");
	gmetvsnjets_da_em_ex->Draw("P");

	lowhtbox1 ->Draw();
	lowmetbox1->Draw();
	boxborder11->Draw();
	boxborder21->Draw();

	// sig1->Draw();
	njets1->Draw();
	njets2->Draw();
	njets3->Draw();
	njets4->Draw();
	njets5->Draw();

	lat->DrawLatex(0.27,0.78, "#splitline{N_{b-jets}^{loose} = 0 (left)}{N_{b-jets}^{loose} #geq 1 (right)}" );
	// Graphs
	gmetvsnjets_da_ee->Draw("P");
	gmetvsnjets_da_em->Draw("P");
	gmetvsnjets_da_mm->Draw("P");
	
	
	leg->Draw();
	// leg2->Draw();
	// regleg->Draw();

	drawTopLine(0.53, 0.8);
	gPad->RedrawAxis();

	// Util::PrintNoEPS(c_temp, "HTvsMET_" + gHiLoLabel[hilo], fOutputDir + fOutputSubDir, NULL);
	Util::PrintPDF(c_temp1, "NjetsvsMET_NoHT", fOutputDir + fOutputSubDir);

	// NJETS vs. MET plot with btags third
	TCanvas *c_temp2 = new TCanvas("C_NjetsvsMET", "Njets vs MET in Data vs MC", 0, 0, 600, 600);
	c_temp2->cd();
	c_temp2->SetRightMargin(0.05);
	c_temp2->SetLeftMargin(0.13);

	hmetvsnjetswithb_da_mm->DrawCopy("axis");

	gmetvsnjetswithb_da_mm_ex->Draw("P");
	gmetvsnjetswithb_da_ee_ex->Draw("P");
	gmetvsnjetswithb_da_em_ex->Draw("P");

	lowhtbox1 ->Draw();
	lowmetbox1->Draw();
	boxborder11->Draw();
	boxborder21->Draw();

	// sig1->Draw();
	njets1->Draw();
	njets2->Draw();
	njets3->Draw();
	njets4->Draw();
	njets5->Draw();

	lat->DrawLatex(0.30,0.78, "N_{b-jets}^{loose} #geq 1" );
	// Graphs
	gmetvsnjetswithb_da_ee->Draw("P");
	gmetvsnjetswithb_da_em->Draw("P");
	gmetvsnjetswithb_da_mm->Draw("P");
	
	
	leg->Draw();
	// leg2->Draw();
	// regleg->Draw();

	drawTopLine(0.53, 0.8);
	gPad->RedrawAxis();

	// Util::PrintNoEPS(c_temp, "HTvsMET_" + gHiLoLabel[hilo], fOutputDir + fOutputSubDir, NULL);
	Util::PrintPDF(c_temp2, "NjetsvsMET_NoHT_withBTags", fOutputDir + fOutputSubDir);
	// Util::SaveAsMacro(c_temp, "HTvsMET_" + gHiLoLabel[hilo], fOutputDir + fOutputSubDir);

	// NJETS vs. MET plot no btags fourth
	TCanvas *c_temp3 = new TCanvas("C_NjetsvsMET", "Njets vs MET in Data vs MC", 0, 0, 600, 600);
	c_temp3->cd();
	c_temp3->SetRightMargin(0.05);
	c_temp3->SetLeftMargin(0.13);

	hmetvsnjetsnob_da_mm->DrawCopy("axis");

	gmetvsnjetsnob_da_mm_ex->Draw("P");
	gmetvsnjetsnob_da_ee_ex->Draw("P");
	gmetvsnjetsnob_da_em_ex->Draw("P");

	lowhtbox1 ->Draw();
	lowmetbox1->Draw();
	boxborder11->Draw();
	boxborder21->Draw();

	// sig1->Draw();
	njets1->Draw();
	njets2->Draw();
	njets3->Draw();
	njets4->Draw();
	njets5->Draw();
	lat->DrawLatex(0.30,0.78, "N_{b-jets}^{loose} = 0" );

	// Graphs
	gmetvsnjetsnob_da_ee->Draw("P");
	gmetvsnjetsnob_da_em->Draw("P");
	gmetvsnjetsnob_da_mm->Draw("P");
	
	
	leg->Draw();
	// leg2->Draw();
	// regleg->Draw();

	drawTopLine(0.53, 0.8);
	gPad->RedrawAxis();

	// Util::PrintNoEPS(c_temp, "HTvsMET_" + gHiLoLabel[hilo], fOutputDir + fOutputSubDir, NULL);
	Util::PrintPDF(c_temp3, "NjetsvsMET_NoHT_noBTags", fOutputDir + fOutputSubDir);
	// Util::SaveAsMacro(c_temp, "HTvsMET_" + gHiLoLabel[hilo], fOutputDir + fOutputSubDir);

	// NJETS vs. MET plot with both fifth
	TCanvas *c_temp4 = new TCanvas("C_NjetsvsMET", "Njets vs MET in Data vs MC", 0, 0, 600, 600);
	c_temp4->cd();
	c_temp4->SetRightMargin(0.05);
	c_temp4->SetLeftMargin(0.13);

	hmetvsnjetsboth_da_mm->DrawCopy("axis");

	gmetvsnjetsboth_da_mm_ex->Draw("P");
	gmetvsnjetsboth_da_ee_ex->Draw("P");
	gmetvsnjetsboth_da_em_ex->Draw("P");

	lowhtbox1 ->Draw();
	lowmetbox1->Draw();
	boxborder11->Draw();
	boxborder21->Draw();

	// sig1->Draw();
	njets1->Draw();
	njets2->Draw();
	njets3->Draw();
	njets4->Draw();
	njets5->Draw();
	lat->DrawLatex(0.30,0.78, "N_{b-jets}^{loose} #geq 0" );

	// Graphs
	gmetvsnjetsboth_da_ee->Draw("P");
	gmetvsnjetsboth_da_em->Draw("P");
	gmetvsnjetsboth_da_mm->Draw("P");
	
	
	leg->Draw();
	// leg2->Draw();
	// regleg->Draw();

	drawTopLine(0.53, 0.8);
	gPad->RedrawAxis();

	// Util::PrintNoEPS(c_temp, "HTvsMET_" + gHiLoLabel[hilo], fOutputDir + fOutputSubDir, NULL);
	Util::PrintPDF(c_temp4, "NjetsvsMET_NoHT_both", fOutputDir + fOutputSubDir);
	// Util::SaveAsMacro(c_temp, "HTvsMET_" + gHiLoLabel[hilo], fOutputDir + fOutputSubDir);

	delete c_temp, c_temp1, c_temp2, c_temp3;
	delete leg, regleg;
	delete hmetvsht_da_mm, hmetvsht_da_ee, hmetvsht_da_em;//, hmetvsht_mc;
	delete gmetvsht_da_mm, gmetvsht_da_ee, gmetvsht_da_em;//, hmetvsht_mc;
	delete hmetvsnjets_da_mm, hmetvsnjets_da_ee, hmetvsnjets_da_em;//, hmetvsnjets_mc;
	delete gmetvsnjets_da_mm, gmetvsnjets_da_ee, gmetvsnjets_da_em;//, hmetvsnjets_mc;
	delete hmetvsnjetswithb_da_mm, hmetvsnjetswithb_da_ee, hmetvsnjetswithb_da_em;
	delete gmetvsnjetswithb_da_mm, gmetvsnjetswithb_da_ee, gmetvsnjetswithb_da_em;
	delete hmetvsnjetsnob_da_mm, hmetvsnjetsnob_da_ee, hmetvsnjetsnob_da_em;
	delete gmetvsnjetsnob_da_mm, gmetvsnjetsnob_da_ee, gmetvsnjetsnob_da_em;
	delete hmetvsnjetsboth_da_mm, hmetvsnjetsboth_da_ee, hmetvsnjetsboth_da_em;
	delete gmetvsnjetsboth_da_mm, gmetvsnjetsboth_da_ee, gmetvsnjetsboth_da_em;
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
void SSDLPlotter::makeChMidvsPtPlots(){
	fOutputSubDir = "ChMisId/";
	char cmd[100];
	sprintf(cmd,"mkdir -p %s%s", fOutputDir.Data(), fOutputSubDir.Data());
	system(cmd);

	TString outputdir = Util::MakeOutputDir(fOutputDir + fOutputSubDir);
	TString outputname = outputdir + "SystematicError.txt";
	TString name = "Electrons";

	TH1D *h_ratio_BB, *h_ratio_EE, *h_ratio_BE;
	h_ratio_BB = new TH1D("Ratio_BB", "P_{cmm} in barrel-barrel", getNPPtBins(Elec), getPPtBins(Elec));
	h_ratio_EE = new TH1D("Ratio_EE", "P_{cmm} in endcap-endcap", getNPPtBins(Elec), getPPtBins(Elec));		
	h_ratio_BE = new TH1D("Ratio_BE", "P_{cmm} in barrel-endcap", getNPPtBins(Elec), getPPtBins(Elec));		

	vector<int> samples = fEGData;
	
	//////////////
	TEfficiency *eff_BB = new TEfficiency("ChMid_BB_pt",  "ChMid_BB_pt",  getNPPtBins(Elec), getPPtBins(Elec));
	TEfficiency *eff_BE = new TEfficiency("ChMid_BE_pt",  "ChMid_BE_pt",  getNPPtBins(Elec), getPPtBins(Elec));
	TEfficiency *eff_EE = new TEfficiency("ChMid_EE_pt",  "ChMid_EE_pt",  getNPPtBins(Elec), getPPtBins(Elec));

	for(size_t i = 0; i < samples.size(); ++i){
		Sample *S  = fSamples[samples[i]];
		Channel *C = &S->region[gRegion[gBaseRegion]][HighPt].ee;

		TEfficiency *tempeff_BB = new TEfficiency(*C->chmid_BB_pt);
		TEfficiency *tempeff_BE = new TEfficiency(*C->chmid_BE_pt);
		TEfficiency *tempeff_EE = new TEfficiency(*C->chmid_EE_pt);
		
		eff_BB->Add(*tempeff_BB);
		eff_BE->Add(*tempeff_BE);
		eff_EE->Add(*tempeff_EE);
		delete tempeff_BB,tempeff_BE,tempeff_EE;
	}
       
	eff_BB->SetName("eff_BB_pt");
	eff_BB->SetLineWidth(2);
	eff_BB->SetMarkerStyle(20);
	eff_BB->SetMarkerColor(kBlack);
	eff_BB->SetMarkerSize(1.5);

	eff_EE->SetName("eff_EE_pt");
	eff_EE->SetLineWidth(2);
	eff_EE->SetLineColor(kRed);
	eff_EE->SetMarkerStyle(20);
	eff_EE->SetMarkerColor(kRed);
	eff_EE->SetMarkerSize(1.5);

	eff_BE->SetName("eff_BE_pt");
	eff_BE->SetLineWidth(2);
	eff_BE->SetLineColor(kBlue);
	eff_BE->SetMarkerStyle(20);
	eff_BE->SetMarkerColor(kBlue);
	eff_BE->SetMarkerSize(1.5);

	//////////////

	float maximum = 0.01;
	h_ratio_BB->SetMaximum(maximum);
	h_ratio_BB->SetMinimum(0.0);
	h_ratio_BB->SetXTitle(convertVarName("ElPt"));
	h_ratio_BB->GetYaxis()->SetTitleOffset(1.2);
	h_ratio_BB->SetYTitle("P_{cmm}");

	h_ratio_BB->SetMarkerColor(kBlack);
	h_ratio_BB->SetMarkerStyle(20);
	h_ratio_BB->SetMarkerSize(1.5);
	h_ratio_BB->SetLineWidth(2);
	h_ratio_BB->SetLineColor(kBlack);
	h_ratio_BB->SetFillColor(kBlack);
	
	h_ratio_EE->SetMarkerColor(kRed);
	h_ratio_EE->SetMarkerStyle(23);
	h_ratio_EE->SetMarkerSize(1.5);
	h_ratio_EE->SetLineWidth(2);
	h_ratio_EE->SetLineColor(kRed);
	h_ratio_EE->SetFillColor(kRed);

	h_ratio_BE->SetMarkerColor(kBlue);
	h_ratio_BE->SetMarkerStyle(23);
	h_ratio_BE->SetMarkerSize(1.5);
	h_ratio_BE->SetLineWidth(2);
	h_ratio_BE->SetLineColor(kBlue);
	h_ratio_BE->SetFillColor(kBlue);

	TLatex *lat = new TLatex();
	lat->SetNDC(kTRUE);
	lat->SetTextColor(kBlack);
	lat->SetTextSize(0.04);
	
	TLegend *leg;
	leg = new TLegend(0.15,0.75,0.35,0.88);
	leg->AddEntry(h_ratio_BB, "Barrel-Barrel", "p");
	leg->AddEntry(h_ratio_BE, "Barrel-EndCap", "p");
	leg->AddEntry(h_ratio_EE, "EndCap-EndCap", "p");
	leg->SetTextSize(0.04);
	leg->SetFillStyle(0);
	leg->SetTextFont(42);
	leg->SetBorderSize(0);

	TCanvas *c_temp = new TCanvas("C_PtRatioPlot", "Ratio vs Pt", 0, 0, 800, 600);
	c_temp->cd();
	h_ratio_BB->DrawCopy("PE 0");
	eff_BB->Draw("PZ 0 same");
	eff_BE->Draw("PZ 0 same");
	eff_EE->Draw("PZ 0 same");
	leg->Draw();
	lat->SetTextSize(0.04);
	lat->DrawLatex(0.62,0.85, name);

	drawTopLine();
	lat->SetTextSize(0.03);
	
	Util::SaveAll (c_temp, fOutputSubDir, fOutputFile);
	Util::PrintPDF(c_temp, "ChMissID_Prob_Pt", fOutputDir + fOutputSubDir);


	float fbb(0.),fee(0.),feb(0.);
	float fbbE(0.),feeE(0.),febE(0.);

	calculateChMisIdProb(fEGData, BB, fbb, fbbE);
	calculateChMisIdProb(fEGData, EB, feb, febE);
	calculateChMisIdProb(fEGData, EE, fee, feeE);

	ofstream OUT(outputname.Data(), ios::trunc);
	OUT << "/////////////////////////////////////////////////////////////////////////////" << endl;
	OUT << " Producing P_{cmm} Syst. error calculation                                   " << endl;
	OUT << "/////////////////////////////////////////////////////////////////////////////" << endl;
	OUT << " " << endl;
 	OUT << "-------------------------------------------------------------------------------------------" << endl;
 	OUT << "       E-ChMisID  ||     Barrel-Barrel    |    Barrel - EndCap   |    Endcap - EndCap    ||" << endl;
 	OUT << "-------------------------------------------------------------------------------------------" << endl;
 	OUT << "  Data            ||";
 	OUT << setw(8) << setprecision(2) << fbb  << " +/- " << setw(8) << setprecision(2) << fbbE  << " |";
 	OUT << setw(8) << setprecision(2) << feb  << " +/- " << setw(8) << setprecision(2) << febE  << " |";
 	OUT << setw(8) << setprecision(2) << fee  << " +/- " << setw(8) << setprecision(2) << feeE  << "  ||" << endl;
 	OUT << "-------------------------------------------------------------------------------------------" << endl << endl;

	int  nbin = getNPPtBins(Elec)-2;
	float  diff_BB[nbin], diff_BE[nbin], diff_EE[nbin];
	float   err_BB[nbin],  err_BE[nbin],  err_EE[nbin];
	float wdiff_BB[nbin],wdiff_BE[nbin],wdiff_EE[nbin];
	
	float syst_BB_max_weighted(0),syst_BB_max_abs(0),syst_BB_mean_weighted(0),syst_BB_mean_abs(0);
	float syst_BE_max_weighted(0),syst_BE_max_abs(0),syst_BE_mean_weighted(0),syst_BE_mean_abs(0);
	float syst_EE_max_weighted(0),syst_EE_max_abs(0),syst_EE_mean_weighted(0),syst_EE_mean_abs(0);
	float syst_max_weighted(0),syst_max_abs(0),syst_mean_weighted(0),syst_mean_abs(0);

	
	OUT << "  ElPt |  P_{cmm}   (barrel-barrel) |   (barrel-endcap)   |  (endcap-endcap) |" << endl;
	OUT << " ---------------------------------------------------------------------------- "  << endl;
	for (int i=0; i<getNPPtBins(Elec)-4; i++){
	  diff_BB [i] =  eff_BB->GetEfficiency(i+3)/2.;
	  diff_BE [i] =  eff_BE->GetEfficiency(i+3)/2.;
	  diff_EE [i] =  eff_EE->GetEfficiency(i+3)/2.;
	  err_BB  [i] = (eff_BB->GetEfficiencyErrorLow(i+3) + eff_BB->GetEfficiencyErrorUp(i+3))/(2.*diff_BB[i]);
	  err_BE  [i] = (eff_BE->GetEfficiencyErrorLow(i+3) + eff_BE->GetEfficiencyErrorUp(i+3))/(2.*diff_BE[i]);
	  err_EE  [i] = (eff_EE->GetEfficiencyErrorLow(i+3) + eff_EE->GetEfficiencyErrorUp(i+3))/(2.*diff_EE[i]);

	  diff_BB [i] = fabs(diff_BB[i] - fbb);
	  diff_BE [i] = fabs(diff_BE[i] - feb);
	  diff_EE [i] = fabs(diff_EE[i] - fee);
	  wdiff_BB[i] = diff_BB[i]*(1/err_BB[i]);
	  wdiff_BE[i] = diff_BE[i]*(1/err_BE[i]);
	  wdiff_EE[i] = diff_EE[i]*(1/err_EE[i]);
	  
	  OUT << Form("  %4.1f |         %6.5f ± %6.5f  |  %6.5f ± %6.5f  |  %5.4f ± %5.4f |", 
		      h_ratio_BB->GetBinCenter(i+3),
		      eff_BB->GetEfficiency(i+3)/2.,err_BB[i]*eff_BB->GetEfficiency(i+3)/2.,
		      eff_BE->GetEfficiency(i+3)/2.,err_BE[i]*eff_BE->GetEfficiency(i+3)/2.,
		      eff_EE->GetEfficiency(i+3)/2.,err_EE[i]*eff_EE->GetEfficiency(i+3)/2.) 
	      << endl;
	  
	  if (wdiff_BB[i] > syst_BB_max_weighted) syst_BB_max_weighted = wdiff_BB[i];
	  if ( diff_BB[i] > syst_BB_max_abs     ) syst_BB_max_abs      = diff_BB[i];
	  if (wdiff_BE[i] > syst_BE_max_weighted) syst_BE_max_weighted = wdiff_BE[i];
	  if ( diff_BE[i] > syst_BE_max_abs     ) syst_BE_max_abs      = diff_BE[i];
	  if (wdiff_EE[i] > syst_EE_max_weighted) syst_EE_max_weighted = wdiff_EE[i];
	  if ( diff_EE[i] > syst_EE_max_abs     ) syst_EE_max_abs      = diff_EE[i];
	  
	  syst_BB_mean_weighted += wdiff_BB[i]/nbin;  
	  syst_BB_mean_abs      +=  diff_BB[i]/nbin;  
	  syst_BE_mean_weighted += wdiff_BE[i]/nbin;  
	  syst_BE_mean_abs      +=  diff_BE[i]/nbin;  
	  syst_EE_mean_weighted += wdiff_EE[i]/nbin;
	  syst_EE_mean_abs      +=  diff_EE[i]/nbin;
	}

	/// all the values are in absolut terms... convert it to % 
	syst_BB_max_weighted  = syst_BB_max_weighted  / fbb * 100;	
	syst_BB_max_abs       = syst_BB_max_abs       / fbb * 100;
	syst_BB_mean_weighted = syst_BB_mean_weighted / fbb * 100;
	syst_BB_mean_abs      = syst_BB_mean_abs      / fbb * 100;

	syst_BE_max_weighted  = syst_BE_max_weighted  / feb * 100;	
	syst_BE_max_abs       = syst_BE_max_abs       / feb * 100;
	syst_BE_mean_weighted = syst_BE_mean_weighted / feb * 100;
	syst_BE_mean_abs      = syst_BE_mean_abs      / feb * 100;

	syst_EE_max_weighted  = syst_EE_max_weighted  / fee * 100;
	syst_EE_max_abs       = syst_EE_max_abs       / fee * 100;
	syst_EE_mean_weighted = syst_EE_mean_weighted / fee * 100;
	syst_EE_mean_abs      = syst_EE_mean_abs      / fee * 100;

	syst_max_weighted   = (syst_BB_max_weighted  + syst_BE_max_weighted  + syst_EE_max_weighted )/3.;
	syst_max_abs        = (syst_BB_max_abs       + syst_BE_max_abs       + syst_EE_max_abs      )/3.;
	syst_mean_weighted  = (syst_BB_mean_weighted + syst_BE_mean_weighted + syst_EE_mean_weighted)/3.;
	syst_mean_abs       = (syst_BB_mean_abs      + syst_BE_mean_abs      + syst_EE_mean_abs     )/3.;
	
	OUT << endl;
	OUT << " ============================================================================== "  << endl;
	OUT << Form(" The max. abs     syst. error is: %4.2f (BB: %4.2f | BE: %4.2f | EE: %4.2f)", syst_max_abs,       syst_BB_max_abs,       syst_BE_max_abs,       syst_EE_max_abs      ) << endl;
	OUT << Form(" The max. weighed syst. error is: %4.2f (BB: %4.2f | BE: %4.2f | EE: %4.2f)", syst_max_weighted,  syst_BB_max_weighted,  syst_BE_max_weighted,  syst_EE_max_weighted ) << endl;
	OUT << Form(" The mean abs     syst. error is: %4.2f (BB: %4.2f | BE: %4.2f | EE: %4.2f)", syst_mean_abs,      syst_BB_mean_abs,      syst_BE_mean_abs,      syst_EE_mean_abs     ) << endl;
	OUT << Form(" The mean weighed syst. error is: %4.2f (BB: %4.2f | BE: %4.2f | EE: %4.2f)", syst_mean_weighted, syst_BB_mean_weighted, syst_BE_mean_weighted, syst_EE_mean_weighted) << endl;
	OUT << " ============================================================================== " <<endl;
	OUT.close();

	delete h_ratio_EE, h_ratio_BB,h_ratio_BE;
	delete eff_BB,eff_EE,eff_BE;
	delete c_temp, lat, leg;
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

	TH1D *h_dummy1, *h_ptratio_data, *h_ptratio_mc, *h_ptratio_qcd, *h_ptratio_wjets, *h_ptratio_ttbar;
	TH2D *h_dummy2;

	TH2D *h2d_ntight, *h2d_nloose;

	if(fp == SigSup){
		h_ptratio_data = new TH1D("Ratio_data", "Tight/Loose Ratio in data", getNFPtBins(chan), getFPtBins(chan));
		h_ptratio_mc   = new TH1D("Ratio_mc",   "Tight/Loose Ratio in MC",   getNFPtBins(chan), getFPtBins(chan));
		h_ptratio_qcd  = new TH1D("Ratio_qcd",  "Tight/Loose Ratio in QCD",  getNFPtBins(chan), getFPtBins(chan));
		h_ptratio_wjets= new TH1D("Ratio_wjets","Tight/Loose Ratio in WJets",getNFPtBins(chan), getFPtBins(chan));
		h_ptratio_ttbar= new TH1D("Ratio_ttbar","Tight/Loose Ratio in ttbar",getNFPtBins(chan), getFPtBins(chan));
		h_dummy1       = new TH1D("dummy1", "dummy1", getNEtaBins(chan), getEtaBins(chan));
		h_dummy2       = new TH2D("dummy2", "dummy2", getNFPtBins(chan), getFPtBins(chan), getNEtaBins(chan), getEtaBins(chan));
	}
	if(fp == ZDecay){
		h_ptratio_data = new TH1D("Ratio_data", "Tight/Loose Ratio in data", getNPPtBins(chan), getPPtBins(chan));
		h_ptratio_mc   = new TH1D("Ratio_mc",   "Tight/Loose Ratio in MC",   getNPPtBins(chan), getPPtBins(chan));
		h_ptratio_qcd  = new TH1D("Ratio_qcd",  "Tight/Loose Ratio in QCD",  getNPPtBins(chan), getPPtBins(chan));
		h_ptratio_wjets= new TH1D("Ratio_wjets","Tight/Loose Ratio in WJets",getNPPtBins(chan), getPPtBins(chan));
		h_ptratio_ttbar= new TH1D("Ratio_ttbar","Tight/Loose Ratio in ttbar",getNPPtBins(chan), getPPtBins(chan));
		h_dummy1       = new TH1D("dummy1", "dummy1", getNEtaBins(chan), getEtaBins(chan));
		h_dummy2       = new TH2D("dummy2", "dummy2", getNPPtBins(chan), getPPtBins(chan), getNEtaBins(chan), getEtaBins(chan));
	}

	vector<int> datasamples;
	vector<int> mcsamples;
	vector<int> qcdsamples;
	vector<int> wjetssamples;
	vector<int> ttbarsamples;

	if(chan == Muon){
		datasamples = fMuData;
//		mcsamples	= fMCBG;
		mcsamples   = fMCBGMuEnr;
//		mcsamples	= fTTJets;
		qcdsamples	= fMuEnr;
		wjetssamples.push_back(WJets);
		ttbarsamples = fTTJets;
	}
	if(chan == Elec){
		datasamples = fEGData;
//		mcsamples	= fMCBG;
		mcsamples   = fMCBGEMEnr;
//		mcsamples	= fTTJets;
		qcdsamples	= fEMEnr;
		wjetssamples.push_back(WJets);
		ttbarsamples = fTTJets;
	}

	calculateRatio(datasamples, chan, fp, h_dummy2, h_ptratio_data, h_dummy1);
	calculateRatio(mcsamples,   chan, fp, h_dummy2, h_ptratio_mc,   h_dummy1);
	calculateRatio(qcdsamples,  chan, fp, h_dummy2, h_ptratio_qcd,   h_dummy1);
	calculateRatio(wjetssamples,  chan, fp, h_dummy2, h_ptratio_wjets,   h_dummy1);
	calculateRatio(ttbarsamples,  chan, fp, h_dummy2, h_ptratio_ttbar,   h_dummy1, false, true);

	//////////////
	float linewidth = 1.;
	TEfficiency *eff_data = getMergedEfficiency(datasamples, chan, fp, 0);
	eff_data->SetName("eff_data_pt");
	eff_data->SetLineWidth(linewidth);
	eff_data->SetMarkerStyle(20);
	eff_data->SetMarkerSize(1.5);

	// TGraphAsymmErrors *eff_mc = getCombEfficiency(mcsamples, chan, fp, 0);
	// eff_mc->SetName("eff_mc_pt");
	// eff_mc->SetMarkerColor(kRed);
	// eff_mc->SetMarkerStyle(23);
	// eff_mc->SetMarkerSize(1.5);
	// eff_mc->SetLineWidth(linewidth);
	// eff_mc->SetLineColor(kRed);
	// eff_mc->SetFillColor(kRed);

	//////////////

	float maximum = 0.2; //0.4;
	if(fp == ZDecay) maximum = 1.1;
	h_ptratio_data->SetMaximum(maximum);
	h_ptratio_mc  ->SetMaximum(maximum);
	h_ptratio_qcd  ->SetMaximum(maximum);
	h_ptratio_wjets  ->SetMaximum(maximum);
	h_ptratio_ttbar  ->SetMaximum(maximum);
	h_ptratio_data->SetMinimum(0.0);
	h_ptratio_mc  ->SetMinimum(0.0);
	h_ptratio_qcd ->SetMinimum(0.0);
	h_ptratio_wjets ->SetMinimum(0.0);
	h_ptratio_ttbar ->SetMinimum(0.0);

	if(chan == Muon) h_ptratio_mc->SetXTitle(convertVarName("MuPt[0]"));
	if(chan == Elec) h_ptratio_mc->SetXTitle(convertVarName("ElPt[0]"));
	h_ptratio_mc->GetYaxis()->SetTitleOffset(1.2);
	h_ptratio_mc->SetYTitle("N_{Tight}/N_{Loose}");

	h_ptratio_data->SetMarkerColor(kBlack);
	h_ptratio_data->SetMarkerStyle(20);
	h_ptratio_data->SetMarkerSize(1.5);
	h_ptratio_data->SetLineWidth(linewidth);
	h_ptratio_data->SetLineColor(kBlack);
	h_ptratio_data->SetFillColor(kBlack);
	
	h_ptratio_mc  ->SetMarkerColor(kRed);
	h_ptratio_mc  ->SetMarkerStyle(23);
	h_ptratio_mc  ->SetMarkerSize(1.5);
	h_ptratio_mc  ->SetLineWidth(linewidth);
	h_ptratio_mc  ->SetLineColor(kRed);
	h_ptratio_mc  ->SetFillColor(kRed);
	h_ptratio_mc  ->GetXaxis()->SetRangeUser(0.,50.);
	
	h_ptratio_qcd  ->SetMarkerColor(kGreen);
	h_ptratio_qcd  ->SetMarkerStyle(23);
	h_ptratio_qcd  ->SetMarkerSize(1.5);
	h_ptratio_qcd  ->SetLineWidth(linewidth);
	h_ptratio_qcd  ->SetLineColor(kGreen);
	h_ptratio_qcd  ->SetFillColor(kGreen);
	h_ptratio_qcd  ->GetXaxis()->SetRangeUser(0.,50.);
	
	h_ptratio_wjets  ->SetMarkerColor(kOrange);
	h_ptratio_wjets  ->SetMarkerStyle(23);
	h_ptratio_wjets  ->SetMarkerSize(1.5);
	h_ptratio_wjets  ->SetLineWidth(2);
	h_ptratio_wjets  ->SetLineColor(kOrange);
	h_ptratio_wjets  ->SetFillColor(kOrange);
	h_ptratio_wjets  ->GetXaxis()->SetRangeUser(0.,50.);
	
	h_ptratio_ttbar  ->SetMarkerColor(kBlue);
	h_ptratio_ttbar  ->SetMarkerStyle(23);
	h_ptratio_ttbar  ->SetMarkerSize(1.5);
	h_ptratio_ttbar  ->SetLineWidth(2);
	h_ptratio_ttbar  ->SetLineColor(kBlue);
	h_ptratio_ttbar  ->SetFillColor(kBlue);
	h_ptratio_ttbar  ->GetXaxis()->SetRangeUser(0.,50.);

	TLatex *lat = new TLatex();
	lat->SetNDC(kTRUE);
	lat->SetTextColor(kBlack);
	lat->SetTextSize(0.04);
	
	TLegend *leg;
	if(fp == SigSup) leg = new TLegend(0.15,0.66,0.35,0.88);
	if(fp == ZDecay) leg = new TLegend(0.15,0.15,0.35,0.37);
	leg->AddEntry(h_ptratio_data, "Data",            "p");
	leg->AddEntry(h_ptratio_mc,   "Simulation",      "p");
	leg->AddEntry(h_ptratio_qcd,  "QCD",             "p");
	leg->AddEntry(h_ptratio_wjets,"WJets",           "p");
	leg->AddEntry(h_ptratio_ttbar,"ttbar (matched)", "p");
	leg->SetTextSize(0.04);
	leg->SetFillStyle(0);
	leg->SetTextFont(42);
	leg->SetBorderSize(0);

	TCanvas *c_temp = new TCanvas("C_PtRatioPlot", "fRatio vs Pt in Data vs MC", 0, 0, 800, 600);
	c_temp->cd();
	// MARC h_ptratio_mc->DrawCopy("axis");
	h_ptratio_mc->DrawCopy("PE 0");
	h_ptratio_qcd->Draw("PE 0 same");
	h_ptratio_ttbar->Draw("PE 0 same");
	// MARC eff_mc->Draw("P same");
	// MARC h_ptratio_data->Draw("PE X0 same");
	eff_data->Draw("PZ 0 same");
	// MARC eff_data->Draw("P same");
	leg->Draw();
	//	lat->DrawLatex(0.70,0.92, Form("L_{int.} = %2.1f fb^{-1}", fLumiNorm/1000.));
	lat->SetTextSize(0.04);
	if(fp == SigSup) lat->DrawLatex(0.62,0.85, pfname + name);
	if(fp == ZDecay) lat->DrawLatex(0.67,0.15, pfname + name);
	double ymean(0.), yrms(0.);
	getWeightedYMeanRMS(h_ptratio_data, ymean, yrms);
	drawTopLine();
	lat->SetTextSize(0.03);
	//	lat->DrawLatex(0.25,0.92, Form("Mean ratio: %4.2f #pm %4.2f", ymean, yrms));

	TString fpname = "F";
	if(fp == ZDecay) fpname = "P";
	
	// Util::PrintNoEPS( c_temp, fpname + "Ratio_" + name + "_Pt", fOutputDir + fOutputSubDir, NULL);
	Util::PrintPDF(c_temp, fpname + "Ratio_" + name + "_Pt", fOutputDir + fOutputSubDir);
	delete h_ptratio_mc, h_ptratio_data, h_ptratio_qcd, h_ptratio_wjets, h_ptratio_ttbar;
	delete c_temp, lat, leg;
	fOutputSubDir = "";
}
void SSDLPlotter::makeFRvsNVPlots(gChannel chan, gFPSwitch fp){
	fOutputSubDir = "Ratios/";
	char cmd[100];
    sprintf(cmd,"mkdir -p %s%s", fOutputDir.Data(), fOutputSubDir.Data());
    system(cmd);

	TString pfname = "Non-prompt ";
	if(fp == ZDecay) pfname = "Prompt ";
	
	TString name;
	if(chan == Muon) name = "Muons";
	if(chan == Elec) name = "Electrons";

	TH1D *h_dummy1, *h_dummy3, *h_nvratio_data, *h_nvratio_mc;
	TH2D *h_dummy2;
	
	TH1D *ratio_dtmc = new TH1D("Ratio_DtMC", "Ratio of Data/MC", 18, 0., 36.); ratio_dtmc->Sumw2();

	if(fp == SigSup){
		h_nvratio_data = new TH1D("Ratio_data", "Tight/Loose Ratio in Data", 18, 0., 36.); h_nvratio_data->Sumw2();
		h_nvratio_mc   = new TH1D("Ratio_mc",   "Tight/Loose Ratio in MC",   18, 0., 36.); h_nvratio_mc->Sumw2();
		h_dummy3       = new TH1D("dummy3", "dummy3", getNFPtBins(chan), getFPtBins(chan));
		h_dummy1       = new TH1D("dummy1", "dummy1", getNEtaBins(chan), getEtaBins(chan));
		h_dummy2       = new TH2D("dummy2", "dummy2", getNFPtBins(chan), getFPtBins(chan), getNEtaBins(chan), getEtaBins(chan));
	}
	if(fp == ZDecay){
		h_nvratio_data = new TH1D("Ratio_data", "Tight/Loose Ratio in Data", 18, 0., 36.); h_nvratio_data->Sumw2();
		h_nvratio_mc   = new TH1D("Ratio_mc",   "Tight/Loose Ratio in MC",   18, 0., 36.); h_nvratio_mc->Sumw2();
		h_dummy3       = new TH1D("dummy3", "dummy3", getNPPtBins(chan), getPPtBins(chan));
		h_dummy1       = new TH1D("dummy1", "dummy1", getNEtaBins(chan), getEtaBins(chan));
		h_dummy2       = new TH2D("dummy2", "dummy2", getNPPtBins(chan), getPPtBins(chan), getNEtaBins(chan), getEtaBins(chan));
	}

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

	calculateRatio(datasamples, chan, fp, h_dummy2, h_dummy3, h_dummy1, h_nvratio_data);
	calculateRatio(mcsamples,   chan, fp, h_dummy2, h_dummy3, h_dummy1, h_nvratio_mc);

	//////////////
	TEfficiency *eff_data = getMergedEfficiency(datasamples, chan, fp, 2);
	eff_data->SetName("eff_data_nv");
	eff_data->SetLineWidth(2);
	eff_data->SetMarkerStyle(20);
	eff_data->SetMarkerSize(1.5);

	// TGraphAsymmErrors *eff_mc = getCombEfficiency(mcsamples, chan, fp, 0);
	// eff_mc->SetName("eff_mc_pt");
	// eff_mc->SetMarkerColor(kRed);
	// eff_mc->SetMarkerStyle(23);
	// eff_mc->SetMarkerSize(1.5);
	// eff_mc->SetLineWidth(2);
	// eff_mc->SetLineColor(kRed);
	// eff_mc->SetFillColor(kRed);

	//////////////

	// Ratio plot:
	float border = 0.3;
	float scale = (1-border)/border;
				
	ratio_dtmc->SetXTitle("N_{Vertices}");
	ratio_dtmc->SetYTitle("");
	ratio_dtmc->GetXaxis()->SetTitleSize(scale * 0.04);
	ratio_dtmc->GetXaxis()->SetLabelSize(scale * ratio_dtmc->GetXaxis()->GetLabelSize());
	ratio_dtmc->GetYaxis()->SetLabelSize(scale * ratio_dtmc->GetYaxis()->GetLabelSize());
	ratio_dtmc->GetXaxis()->SetTickLength(scale * ratio_dtmc->GetXaxis()->GetTickLength());
	ratio_dtmc->GetYaxis()->SetTickLength(ratio_dtmc->GetYaxis()->GetTickLength());

	// ratio_dtmc->GetXaxis()->SetLabelSize(scale*0.06);
	// ratio_dtmc->GetXaxis()->SetLabelOffset(0.02);
	// ratio_dtmc->GetXaxis()->SetTitleOffset(1.20);
	// for(size_t j = 1; j <= ratio_dtmc->GetXaxis()->GetNbins(); ++j) ratio_dtmc->GetXaxis()->SetBinLabel(j, Form("%d", j-1));


	ratio_dtmc->SetFillStyle(1001);
	ratio_dtmc->SetLineWidth(1);
	ratio_dtmc->SetFillColor(  kGray+1);
	ratio_dtmc->SetLineColor(  kGray+1);
	ratio_dtmc->SetMarkerColor(kGray+1);

	ratio_dtmc->Divide(h_nvratio_data, h_nvratio_mc);

	// Plotting range
	float maximum = 0.4;
	if(chan == Muon) maximum = 0.3;
	if(fp == ZDecay) maximum = 1.0;
	h_nvratio_data->SetMaximum(maximum);
	h_nvratio_mc  ->SetMaximum(maximum);
	h_nvratio_data->SetMinimum(0.0);
	h_nvratio_mc  ->SetMinimum(0.0);

	h_nvratio_data->SetMarkerColor(kBlack);
	h_nvratio_data->SetMarkerStyle(20);
	h_nvratio_data->SetMarkerSize(1.5);
	h_nvratio_data->SetLineWidth(2);
	h_nvratio_data->SetLineColor(kBlack);
	h_nvratio_data->SetFillColor(kBlack);
	
	h_nvratio_mc->SetMarkerColor(kRed);
	h_nvratio_mc->SetMarkerStyle(23);
	h_nvratio_mc->SetMarkerSize(1.5);
	h_nvratio_mc->SetLineWidth(2);
	h_nvratio_mc->SetLineColor(kRed);
	h_nvratio_mc->SetFillColor(kRed);
	h_nvratio_mc->SetYTitle("N_{Tight}/N_{Loose}");

	// Remove bin labels for plot pad
	for(size_t j = 1; j <= h_nvratio_data->GetNbinsX(); ++j) h_nvratio_data->GetXaxis()->SetBinLabel(j, "");
	for(size_t j = 1; j <= h_nvratio_mc  ->GetNbinsX(); ++j) h_nvratio_mc  ->GetXaxis()->SetBinLabel(j, "");

	TLatex *lat = new TLatex();
	lat->SetNDC(kTRUE);
	lat->SetTextColor(kBlack);
	lat->SetTextSize(0.04);
	
	TLegend *leg;
	if(fp == SigSup) leg = new TLegend(0.15,0.75,0.35,0.88);
	if(fp == ZDecay) leg = new TLegend(0.15,0.05,0.35,0.23);
	leg->AddEntry(h_nvratio_data, "Data",       "p");
	leg->AddEntry(h_nvratio_mc,   "Simulation", "p");
	leg->SetTextSize(0.04);
	leg->SetFillStyle(0);
	leg->SetTextFont(42);
	leg->SetBorderSize(0);


	TCanvas *c_temp = new TCanvas("C_PtRatioPlot", "fRatio vs NVertices in Data vs MC", 0, 0, 800, 800);
	c_temp->cd();

	TPad *p_plot  = new TPad("plotpad",  "Pad containing the plot", 0.00, border, 1.00, 1.00, 0, 0);
	p_plot->SetBottomMargin(0.015);
	p_plot->Draw();
	TPad *p_ratio = new TPad("ratiopad", "Pad containing the ratio", 0.00, 0.00, 1.00, border, 0, 0);
	p_ratio->SetTopMargin(0.025);
	p_ratio->SetBottomMargin(0.35);
	p_ratio->Draw();

	p_ratio->cd();
	ratio_dtmc->GetYaxis()->SetNdivisions(505);
	ratio_dtmc->SetMaximum(1.5);
	ratio_dtmc->SetMinimum(0.8);
	if(fp == ZDecay) ratio_dtmc->SetMaximum(1.1);
	if(fp == ZDecay) ratio_dtmc->SetMinimum(0.9);
	ratio_dtmc->DrawCopy("E2 ");
	TLine *l3 = new TLine(ratio_dtmc->GetXaxis()->GetXmin(), 1.00, ratio_dtmc->GetXaxis()->GetXmax(), 1.00);
	l3->SetLineWidth(2);
	l3->SetLineStyle(7);
	l3->Draw();
	gPad->RedrawAxis();
	p_ratio->Draw();


	p_plot->cd();
	h_nvratio_mc->DrawCopy("PE 0");
	// eff_mc->Draw("P same");
	// h_nvratio_data->Draw("PE X0 same");
	eff_data->Draw("PZ 0 same");
	leg->Draw();
	// lat->DrawLatex(0.70,0.92, Form("L_{int.} = %2.1f fb^{-1}", fLumiNorm/1000.));
	lat->SetTextSize(0.04);
	if(fp == SigSup) lat->DrawLatex(0.62,0.85, pfname + name);
	if(fp == ZDecay) lat->DrawLatex(0.67,0.07, pfname + name);
	drawTopLine();

	c_temp->Update();
	
	TString fpname = "F";
	if(fp == ZDecay) fpname = "P";
	
	Util::PrintPDF(c_temp, fpname + "Ratio_" + name + "_NV", fOutputDir + fOutputSubDir);
	delete h_nvratio_mc, h_nvratio_data;
	delete c_temp, lat, leg;
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

	TH1D *h_dummy1 = new TH1D("dummy1", "dummy1", getNFPtBins(chan), getFPtBins(chan));
	TH2D *h_dummy2 = new TH2D("dummy2", "dummy2", getNFPtBins(chan), getFPtBins(chan), getNEtaBins(chan), getEtaBins(chan));

	TH1D *h_etaratio_data = new TH1D("Ratio_data", "Tight/Loose Ratio in data", getNEtaBins(chan), getEtaBins(chan));
	TH1D *h_etaratio_mc   = new TH1D("Ratio_mc",   "Tight/Loose Ratio in MC",   getNEtaBins(chan), getEtaBins(chan));
	TH1D *h_etaratio_qcd  = new TH1D("Ratio_qcd",  "Tight/Loose Ratio in QCD",  getNEtaBins(chan), getEtaBins(chan));
	TH1D *h_etaratio_wjets= new TH1D("Ratio_wjets","Tight/Loose Ratio in WJets",getNEtaBins(chan), getEtaBins(chan));

	vector<int> datasamples;
	vector<int> mcsamples;
	vector<int> qcdsamples;
	vector<int> wjetssamples;

	if(chan == Muon){
		datasamples = fMuData;
//		mcsamples	= fMCBG;
		mcsamples   = fMCBGMuEnr;
//		mcsamples	= fTTJets;
		qcdsamples   = fMuEnr;
		wjetssamples.push_back(WJets);
	}
	if(chan == Elec){
		datasamples = fEGData;
//		mcsamples	= fMCBG;
		mcsamples   = fMCBGEMEnr;
//		mcsamples	= fTTJets;
		qcdsamples   = fEMEnr;
		wjetssamples.push_back(WJets);
	}

	calculateRatio(datasamples, chan, SigSup, h_dummy2, h_dummy1, h_etaratio_data);
	calculateRatio(mcsamples,   chan, SigSup, h_dummy2, h_dummy1, h_etaratio_mc);
	calculateRatio(qcdsamples,  chan, SigSup, h_dummy2, h_dummy1, h_etaratio_qcd);
	calculateRatio(wjetssamples,chan, SigSup, h_dummy2, h_dummy1, h_etaratio_wjets);

	//////////////
	TEfficiency *eff_data = getMergedEfficiency(datasamples, chan, SigSup, 1);
	eff_data->SetName("eff_data_eta");
	eff_data->SetLineWidth(2);
	eff_data->SetMarkerStyle(20);
	eff_data->SetMarkerSize(1.5);

	float max = 0.2; //0.4;
	if(chan==Elec) max = 0.8;
	h_etaratio_data->SetMaximum(max);
	h_etaratio_mc  ->SetMaximum(max);
	h_etaratio_qcd  ->SetMaximum(max);
	h_etaratio_wjets  ->SetMaximum(max);
	h_etaratio_data->SetMinimum(0.0);
	h_etaratio_mc  ->SetMinimum(0.0);
	h_etaratio_qcd  ->SetMinimum(0.0);
	h_etaratio_wjets  ->SetMinimum(0.0);

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
	
	h_etaratio_qcd  ->SetMarkerColor(kBlue);
	h_etaratio_qcd  ->SetMarkerStyle(23);
	h_etaratio_qcd  ->SetMarkerSize(1.5);
	h_etaratio_qcd  ->SetLineWidth(2);
	h_etaratio_qcd  ->SetLineColor(kBlue);
	h_etaratio_qcd  ->SetFillColor(kBlue);
	
	h_etaratio_wjets  ->SetMarkerColor(kYellow);
	h_etaratio_wjets  ->SetMarkerStyle(23);
	h_etaratio_wjets  ->SetMarkerSize(1.5);
	h_etaratio_wjets  ->SetLineWidth(2);
	h_etaratio_wjets  ->SetLineColor(kYellow);
	h_etaratio_wjets  ->SetFillColor(kYellow);

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
	leg->AddEntry(h_etaratio_qcd,   "QCD", "p");
	leg->AddEntry(h_etaratio_wjets,   "WJets", "p");
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
	//     lat.DrawLatex(0.23, 0.79, "#int L dt = XXX pb^{-1},   #sqrt{s} = 8 TeV");
	//     lat.DrawLatex(0.83, 0.88, name);
	
	h_etaratio_mc->DrawCopy("PE 0");
	// MARC h_etaratio_mc->DrawCopy("axis");
	// MARC eff_data->Draw("P same");
	h_etaratio_qcd->Draw("PE 0 same");
	h_etaratio_wjets->Draw("PE 0 same");
	eff_data->Draw("PZ 0 same");
	leg->Draw();
	drawTopLine();
	lat->SetTextSize(0.04);
	lat->DrawLatex(0.62,0.85, pfname + name);
	double ymean(0.), yrms(0.);
	getWeightedYMeanRMS(h_etaratio_data, ymean, yrms);
	lat->SetTextSize(0.03);
	//	lat->DrawLatex(0.25,0.92, Form("Mean ratio: %4.2f #pm %4.2f", ymean, yrms));
	
	// Util::PrintNoEPS( c_temp, "FRatio_" + name + "_Eta", fOutputDir + fOutputSubDir, NULL);
	Util::PrintPDF(   c_temp, "FRatio_" + name + "_Eta", fOutputDir + fOutputSubDir);
	delete h_etaratio_mc, h_etaratio_data, h_etaratio_qcd, h_etaratio_wjets;
	// delete c_temp;
	delete c_temp, lat, leg;
	fOutputSubDir = "";
}
void SSDLPlotter::makeRatioPlots(gChannel chan){
	TString name;
	if(chan == Muon) name = "Muons";
	if(chan == Elec) name = "Electrons";

	fOutputSubDir = "Ratios/" + name + "/";
	char cmd[100];
    sprintf(cmd,"mkdir -p %s%s", fOutputDir.Data(), fOutputSubDir.Data());
    system(cmd);

	vector<int> datasamples;
	vector<int> mcsamples;

	if(chan == Muon){
		datasamples = fMuData;
//		mcsamples	= fMCBG;
		mcsamples   = fMCBGMuEnr;
//		mcsamples   = fTTJets;
	}
	if(chan == Elec){
		datasamples = fEGData;
//		mcsamples	= fMCBG;
		mcsamples   = fMCBGEMEnr;
//		mcsamples   = fTTJets;
	}

	// Customization
	TString axis_name[gNRatioVars] = {"N_{Jets}",  "H_{T} (GeV)", "P_{T}(Hardest Jet) (GeV)", "N_{Vertices}", "p_{T}(Closest Jet) (GeV)", "p_{T}(Away Jet) (GeV)", "N_{BJets}", "E_{T}^{miss} (GeV)", "m_{T} (GeV)"};

	for(size_t i = 0; i < gNRatioVars; ++i){
		TH1D *h_ratio_data = getFRatio(datasamples, chan, i);
		TH1D *h_ratio_mc   = getFRatio(mcsamples,   chan, i);
		h_ratio_data->SetName(Form("FRatio_%s_data", FRatioPlots::var_name[i].Data()));
		h_ratio_mc  ->SetName(Form("FRatio_%s_mc",   FRatioPlots::var_name[i].Data()));

		float max = 0.2;
		//		if(i==8) max = 1.0;
		if(chan==Elec) max = 0.5;
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
		h_ratio_mc  ->DrawCopy("AXIS");
		h_ratio_mc  ->DrawCopy("PE X0 SAME");
		h_ratio_data->DrawCopy("PE SAME");
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
void SSDLPlotter::make2DRatioPlots(gChannel chan){
	TH2D* histo = fH2D_MufRatio;
	if(chan == Elec) histo = fH2D_ElfRatio;
	if(!histo){
		cerr << "SSDLPlotter::make2DRatioPlots ==> Warning: ratio histo not filled, exiting" << endl;
		exit(-1);
	}

	TString name;
	if(chan == Muon) name = "Muons";
	if(chan == Elec) name = "Electrons";

	fOutputSubDir = "Ratios/";
	char cmd[100];
	sprintf(cmd,"mkdir -p %s%s", fOutputDir.Data(), fOutputSubDir.Data());
	system(cmd);

	histo->SetMinimum(0.0);
	histo->SetMaximum(1.0);

	TCanvas *c_temp = new TCanvas("C_Temp", "fRatio", 0, 0, 800, 600);
	c_temp->cd();

	// gPad->SetLogy();
	histo->DrawCopy("colz text");
	drawTopLine();

	Util::PrintPDF(c_temp, "FRatio2D_" + name, fOutputDir + fOutputSubDir);
	delete c_temp;
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
		mcsamples   = fMCBGEMEnr;
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
			float scale = fLumiNorm / S->getLumi();
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

void SSDLPlotter::makeFakeGenIDTables(){
	TString outputdir = Util::MakeOutputDir(fOutputDir + "FakeGenIDTables");
// 	for(size_t i = 0; i < gNREGIONS; ++i){
	TString outputnameMuQCD = outputdir + "MuFakeGenIDTable_" + gRegions[gRegion[gBaseRegion]]->sname + "_QCD.txt";
	TString outputnameElQCD = outputdir + "ElFakeGenIDTable_" + gRegions[gRegion[gBaseRegion]]->sname + "_QCD.txt";
	
	makeFakeGenIDTable(fMuEnr, Muon, SigSup, outputnameMuQCD);
	makeFakeGenIDTable(fEMEnr, Elec, SigSup, outputnameElQCD);
	
	TString outputnameMuTTJets = outputdir + "MuFakeGenIDTable_" + gRegions[gRegion[gBaseRegion]]->sname + "_TTJets.txt";
	TString outputnameElTTJets = outputdir + "ElFakeGenIDTable_" + gRegions[gRegion[gBaseRegion]]->sname + "_TTJets.txt";
	
	makeFakeGenIDTable(fTTJets, Muon, SigSup, outputnameMuTTJets);
	makeFakeGenIDTable(fTTJets, Elec, SigSup, outputnameElTTJets);
	
	TString outputnameMuMCBG = outputdir + "MuFakeGenIDTable_" + gRegions[gRegion[gBaseRegion]]->sname + "_MCBG.txt";
	TString outputnameElMCBG = outputdir + "ElFakeGenIDTable_" + gRegions[gRegion[gBaseRegion]]->sname + "_MCBG.txt";
	
	makeFakeGenIDTable(fMCBG, Muon, SigSup, outputnameMuMCBG);
	makeFakeGenIDTable(fMCBG, Elec, SigSup, outputnameElMCBG);
}

void SSDLPlotter::makeFakeGenIDTable(vector<int> samples, gChannel chan, gFPSwitch fp, TString filename){
	ofstream OUT(filename.Data(), ios::trunc);
	TH1D* h_temp_loose         = new TH1D("h_temp_loose"       ,  "h_temp_loose"       ,  1001, -0.5, 1000.5);
	TH1D* h_temp_tight         = new TH1D("h_temp_tight"       ,  "h_temp_tight"       ,  1001, -0.5, 1000.5);
	TH1D* h_temp_loose_sig     = new TH1D("h_temp_loose_sig"   ,  "h_temp_loose_sig"   ,  1001, -0.5, 1000.5);
	TH1D* h_temp_tight_sig     = new TH1D("h_temp_tight_sig"   ,  "h_temp_tight_sig"   ,  1001, -0.5, 1000.5);
	TH1D* h_temp_loose_sigSup  = new TH1D("h_temp_loose_sigSup",  "h_temp_loose_sigSup",  1001, -0.5, 1000.5);
	TH1D* h_temp_tight_sigSup  = new TH1D("h_temp_tight_sigSup",  "h_temp_tight_sigSup",  1001, -0.5, 1000.5);
	OUT << "/////////////////////////////////////////////////////////////////////////////" << endl;
	OUT << "\n\n\n";
	if(chan == Muon) OUT << "muon fake ratio contributions in" << endl;
	if(chan == Elec) OUT << "electron fake ratio contributions in" << endl;
	for(size_t i = 0; i < samples.size(); ++i){
		int index = samples[i];
		Sample *S = fSamples[index];
		Channel *C;
		if(chan == Muon) C = &S->region[gRegion[gBaseRegion]][HighPt].mm;
		if(chan == Elec) C = &S->region[gRegion[gBaseRegion]][HighPt].ee;
		if(fp == SigSup){
			OUT << S->sname << "\t(" << setw(10) << S->ngen << " gen events)" << endl;
//			OUT << endl;
//			OUT << "=============================================" << endl;
//			OUT << "all truth matched              &         same sign pair        &         signal suppressed    \\\\" << endl;
//			OUT << "---------------------------------------------" << endl;
//			OUT << setw(6) << "PDG ID" << " & " << setw(9) << "ntight" << " & " << setw(9) << "nloose" << " & " << setw(9) << "ratio";
//			OUT << "    & " << setw(9) << "ntight" << " & " << setw(9) << "nloose" << " & " << setw(9) << "ratio";
//			OUT << "    & " << setw(9) << "ntight" << " & " << setw(9) << "nloose" << " & " << setw(9) << "ratio" << " \\\\" << endl;
//			OUT << "---------------------------------------------" << endl;
//			for (int i = 0; i < C->fnloose_sigSup_genID->GetNbinsX(); i++){
//				if (C->fnloose_genID->GetBinContent(i) > 0 || C->fnloose_sig_genID->GetBinContent(i) > 0 || C->fnloose_sigSup_genID->GetBinContent(i) > 0) {
//					
//					OUT << setw(6) << i-1 << " & " << setw(9) << C->fntight_genID->GetBinContent(i) << " & " << setw(9) << C->fnloose_genID->GetBinContent(i) << " & " << setw(9) << setprecision(4) << C->fntight_genID->GetBinContent(i)/C->fnloose_genID->GetBinContent(i);
//					OUT << "    & " << setw(9) << C->fntight_sig_genID->GetBinContent(i) << " & " << setw(9) << C->fnloose_sig_genID->GetBinContent(i) << " & " << setw(9) << setprecision(4) << C->fntight_sig_genID->GetBinContent(i)/C->fnloose_sig_genID->GetBinContent(i);
//					OUT << "    & " << setw(9) << C->fntight_sigSup_genID->GetBinContent(i) << " & " << setw(9) << C->fnloose_sigSup_genID->GetBinContent(i) << " & " << setw(9) << setprecision(4) << C->fntight_sigSup_genID->GetBinContent(i)/C->fnloose_sigSup_genID->GetBinContent(i) << " \\\\" << endl;
//					
//					
//					
//				}
//			}
//			OUT << "---------------------------------------------" << endl;
//			OUT << setw(6) << "Total " << " & " << setw(9) << C->fntight_sigSup_genID->Integral() << " & " << setw(9) << C->fnloose_sigSup_genID->Integral() << " & " << setw(9) << setprecision(4) << C->fntight_sigSup_genID->Integral()/C->fnloose_sigSup_genID->Integral() << " \\\\" << endl;
//			//		OUT << setw(6) << "Total " << " & " << setw(9) << C->fntight_sigSup_genID->GetEntries() << " & " << setw(9) << C->fnloose_sigSup_genID->GetEntries() << " & " << setw(9) << setprecision(4) << C->fntight_sigSup_genID->GetEntries()/C->fnloose_sigSup_genID->GetEntries() << " \\\\" << endl;
//			OUT << "=============================================" << endl;
			h_temp_loose       ->Add(C->fnloose_genID);
			h_temp_tight       ->Add(C->fntight_genID);
			h_temp_loose_sig   ->Add(C->fnloose_sig_genID);
			h_temp_tight_sig   ->Add(C->fntight_sig_genID);
			h_temp_loose_sigSup->Add(C->fnloose_sigSup_genID);
			h_temp_tight_sigSup->Add(C->fntight_sigSup_genID);
		}
		
		
	}
	OUT << "\n\n\n";
	OUT << "=============================================" << endl;
	OUT << "=============================================" << endl;
	if(fp == SigSup){
		OUT << "\n\n\n";
		if(chan == Muon) OUT << "muon fake ratio contributions in total" << endl;
		if(chan == Elec) OUT << "electron fake ratio contributions in total" << endl;
		OUT << endl;
		OUT << "\\hline %===================================================================================================================" << endl;
		OUT << "&   \\multicolumn{3}{c|}{all truth matched}   & \\multicolumn{3}{c|}{same sign pair} & \\multicolumn{3}{c}{signal suppressed} \\\\" << endl;
//		OUT << "all truth matched                             &         same sign pair               &         signal suppressed         \\\\" << endl;
		OUT << "\\hline %-------------------------------------------------------------------------------------------------------------------" << endl;
		OUT << setw(6) << "PDG ID" << " & " << setw(9) << "ntight" << " & " << setw(9) << "nloose" << " & " << setw(9) << "ratio";
		OUT << "    & " << setw(9) << "ntight" << " & " << setw(9) << "nloose" << " & " << setw(9) << "ratio";
		OUT << "    & " << setw(9) << "ntight" << " & " << setw(9) << "nloose" << " & " << setw(9) << "ratio" << " \\\\" << endl;
		OUT << "\\hline %-------------------------------------------------------------------------------------------------------------------" << endl;
		for (int i = 0; i < h_temp_loose->GetNbinsX(); i++){
			if (h_temp_loose->GetBinContent(i) > 0 || h_temp_loose_sig->GetBinContent(i) > 0 || h_temp_loose_sigSup->GetBinContent(i) > 0) {
				
				OUT << setw(6) << i-1 << " & " << setw(9) << h_temp_tight       ->GetBinContent(i) << " & " << setw(9) << h_temp_loose       ->GetBinContent(i) << " & " << setw(9) << setprecision(4) << h_temp_tight       ->GetBinContent(i)/h_temp_loose       ->GetBinContent(i);
				OUT << "    & " << setw(9) <<                h_temp_tight_sig   ->GetBinContent(i) << " & " << setw(9) << h_temp_loose_sig   ->GetBinContent(i) << " & " << setw(9) << setprecision(4) << h_temp_tight_sig   ->GetBinContent(i)/h_temp_loose_sig   ->GetBinContent(i);
				OUT << "    & " << setw(9) <<                h_temp_tight_sigSup->GetBinContent(i) << " & " << setw(9) << h_temp_loose_sigSup->GetBinContent(i) << " & " << setw(9) << setprecision(4) << h_temp_tight_sigSup->GetBinContent(i)/h_temp_loose_sigSup->GetBinContent(i) << " \\\\" << endl;
				
				
				
			}
		}
		OUT << "\\hline %-------------------------------------------------------------------------------------------------------------------" << endl;
		OUT << setw(6) << "Total " << " & " << setw(9) << h_temp_tight       ->Integral() << " & " << setw(9) << h_temp_loose       ->Integral() << " & " << setw(9) << setprecision(4) << h_temp_tight       ->Integral()/h_temp_loose       ->Integral();
		OUT << "    & " << setw(9) <<                     h_temp_tight_sig   ->Integral() << " & " << setw(9) << h_temp_loose_sig   ->Integral() << " & " << setw(9) << setprecision(4) << h_temp_tight_sig   ->Integral()/h_temp_loose_sig   ->Integral();
		OUT << "    & " << setw(9) <<                     h_temp_tight_sigSup->Integral() << " & " << setw(9) << h_temp_loose_sigSup->Integral() << " & " << setw(9) << setprecision(4) << h_temp_tight_sigSup->Integral()/h_temp_loose_sigSup->Integral() << " \\\\" << endl;
		//		OUT << setw(6) << "Total " << " & " << setw(9) << h_temp_tight->GetEntries() << " & " << setw(9) << h_temp_loose->GetEntries() << " & " << setw(9) << setprecision(4) << h_temp_tight->GetEntries()/h_temp_loose->GetEntries() << " \\\\" << endl;
		OUT << "\\hline %===================================================================================================================" << endl;
	}
	
	OUT.close();
}
void SSDLPlotter::makeMIDIsolationPlots(gChannel chan, gFPSwitch fp){
	fOutputSubDir = "IsolationPlots/";
	char cmd[100];
	sprintf(cmd,"mkdir -p %s%s", fOutputDir.Data(), fOutputSubDir.Data());
	system(cmd);
	
//	TString pfname = "Non-prompt ";
//	if(fp == ZDecay) pfname = "Prompt ";
	
	TString name;
	if(chan == Muon) name = "MuIsolation_GenID";
	if(chan == Elec) name = "ElIsolation_GenID";
	if(fp == SigSup) name = name + "_sigSup";
	else name = name + "_sig";
	
	vector<int> ttbarsamples;
	vector<int> qcdsamples;
	
	ttbarsamples                 = fTTJets;
	if(chan == Muon) qcdsamples  = fMuEnr;
	if(chan == Elec) qcdsamples  = fEMEnr;
	
	TH1D *hMID24_Iso_ttbar, *hMID500_Iso_ttbar, *hMID400_Iso_ttbar, *hMID15_Iso_ttbar;
	TH1D *hMID24_Iso_qcd, *hMID500_Iso_qcd, *hMID400_Iso_qcd, *hMID15_Iso_qcd;
	
	hMID24_Iso_ttbar	= new TH1D("_sigSup_MID24_Iso_ttbar" , "sigSup_MID24_Iso_ttbar" , 200, 0., 1.);	hMID24_Iso_ttbar ->Sumw2();
	hMID500_Iso_ttbar	= new TH1D("_sigSup_MID500_Iso_ttbar", "sigSup_MID500_Iso_ttbar", 200, 0., 1.);	hMID500_Iso_ttbar->Sumw2();
	hMID400_Iso_ttbar	= new TH1D("_sigSup_MID400_Iso_ttbar", "sigSup_MID400_Iso_ttbar", 200, 0., 1.);	hMID400_Iso_ttbar->Sumw2();
	hMID15_Iso_ttbar	= new TH1D("_sigSup_MID15_Iso_ttbar" , "sigSup_MID15_Iso_ttbar" , 200, 0., 1.);	hMID15_Iso_ttbar ->Sumw2();
	hMID24_Iso_qcd	    = new TH1D("_sigSup_MID24_Iso_qcd"   , "sigSup_MID24_Iso_qcd"   , 200, 0., 1.);	hMID24_Iso_qcd   ->Sumw2();
	hMID500_Iso_qcd	    = new TH1D("_sigSup_MID500_Iso_qcd"  , "sigSup_MID500_Iso_qcd"  , 200, 0., 1.);	hMID500_Iso_qcd  ->Sumw2();
	hMID400_Iso_qcd	    = new TH1D("_sigSup_MID400_Iso_qcd"  , "sigSup_MID400_Iso_qcd"  , 200, 0., 1.);	hMID400_Iso_qcd  ->Sumw2();
	hMID15_Iso_qcd    	= new TH1D("_sigSup_MID15_Iso_qcd"   , "sigSup_MID15_Iso_qcd"   , 200, 0., 1.);	hMID15_Iso_qcd   ->Sumw2();
	
	
	getMIDIsolationPlots(ttbarsamples, chan, fp, hMID24_Iso_ttbar, hMID500_Iso_ttbar, hMID400_Iso_ttbar, hMID15_Iso_ttbar);
	getMIDIsolationPlots(qcdsamples, chan, fp, hMID24_Iso_qcd, hMID500_Iso_qcd, hMID400_Iso_qcd, hMID15_Iso_qcd);
	drawMIDIsolationPlots(name + "_qcd", Muon, fp, hMID24_Iso_qcd, hMID500_Iso_qcd, hMID400_Iso_qcd, hMID15_Iso_qcd);
	drawMIDIsolationPlots(name + "_ttbar", Muon, fp, hMID24_Iso_ttbar, hMID500_Iso_ttbar, hMID400_Iso_ttbar, hMID15_Iso_ttbar);

	delete hMID24_Iso_ttbar, hMID500_Iso_ttbar, hMID400_Iso_ttbar, hMID15_Iso_ttbar;
	delete hMID24_Iso_qcd, hMID500_Iso_qcd, hMID400_Iso_qcd, hMID15_Iso_qcd;
	
	fOutputSubDir = "";
}
void SSDLPlotter::getMIDIsolationPlots(vector<int> samples, gChannel chan, gFPSwitch fp, TH1D*& hMID24iso, TH1D*& hMID500iso, TH1D*& hMID400iso, TH1D*& hMID15iso){
	cout << "SSDLPlotter::getMIDIsolationPlots" << endl;
	for (size_t i = 0; i < samples.size(); ++i) {
		int index = samples[i];
		Sample *S = fSamples[index];
		Channel *C;
		if(chan == Muon) C = &S->region[gRegion[gBaseRegion]][HighPt].mm;
		if(chan == Elec) C = &S->region[gRegion[gBaseRegion]][HighPt].ee;
		float scale = fLumiNorm / S->getLumi();
		scale = 1.;
		cout << "adding histo" << endl;
		if (fp == SigSup) {
			hMID24iso ->Add(C->sigSup_MID24_Iso , scale);
			hMID500iso->Add(C->sigSup_MID500_Iso, scale);
			hMID400iso->Add(C->sigSup_MID400_Iso, scale);
			hMID15iso ->Add(C->sigSup_MID15_Iso , scale);
		}
		else {
			hMID24iso ->Add(C->sig_MID24_Iso , scale);
			hMID500iso->Add(C->sig_MID500_Iso, scale);
			hMID400iso->Add(C->sig_MID400_Iso, scale);
			hMID15iso ->Add(C->sig_MID15_Iso , scale);
		}
	}
}
void SSDLPlotter::drawMIDIsolationPlots(TString title, gChannel chan, gFPSwitch fp, TH1D* hMID24iso, TH1D* hMID500iso, TH1D* hMID400iso, TH1D* hMID15iso){
	
	cout << "SSDLPlotter::drawMIDIsolationPlots" << endl;
//	
//	TH1D *hMID24iso  = (TH1D*)hMID24iso2 ->Clone();
//	TH1D *hMID500iso = (TH1D*)hMID500iso2->Clone();
//	TH1D *hMID400iso = (TH1D*)hMID400iso2->Clone();
//	TH1D *hMID15iso  = (TH1D*)hMID15iso2 ->Clone();
	
	
	hMID24iso ->Rebin(4);
	hMID500iso->Rebin(4);
	hMID400iso->Rebin(4);
	hMID15iso ->Rebin(4);
	
	
	hMID24iso ->Scale(1./hMID24iso ->Integral());
	hMID500iso->Scale(1./hMID500iso->Integral());
	hMID400iso->Scale(1./hMID400iso->Integral());
	hMID15iso ->Scale(1./hMID15iso ->Integral());
	
	float linewidth = 1.;
	float maximum = 0.1;
	hMID24iso ->SetMaximum(maximum);
	hMID500iso->SetMaximum(maximum);
	hMID400iso->SetMaximum(maximum);
	hMID15iso ->SetMaximum(maximum);
	hMID24iso ->SetMinimum(0.);
	hMID500iso->SetMinimum(0.);
	hMID400iso->SetMinimum(0.);
	hMID15iso ->SetMinimum(0.);
	
	
	
//	if(chan == Muon) hMID24iso->SetXTitle(convertVarName("MuPt[0]"));
//	if(chan == Elec) hMID24iso->SetXTitle(convertVarName("ElPt[0]"));
//	hMID500iso ->GetYaxis()->SetTitleOffset(1.2);
	hMID500iso ->SetXTitle("PF Isolation");
	
	hMID24iso ->SetMarkerColor(kBlack);
	hMID24iso ->SetMarkerStyle(20);
	hMID24iso ->SetMarkerSize(1.);
	hMID24iso ->SetLineWidth(linewidth);
	hMID24iso ->SetLineColor(kBlack);
	hMID24iso ->SetFillColor(kBlack);
	
	hMID500iso->SetMarkerColor(kRed);
	hMID500iso->SetMarkerStyle(20);
	hMID500iso->SetMarkerSize(1.);
	hMID500iso->SetLineWidth(linewidth);
	hMID500iso->SetLineColor(kYellow);
	hMID500iso->SetFillColor(kYellow);
//	hMID500iso->GetXaxis()->SetRangeUser(0.,50.);
	
	hMID400iso->SetMarkerColor(kBlue);
	hMID400iso->SetMarkerStyle(20);
	hMID400iso->SetMarkerSize(1.);
	hMID400iso->SetLineWidth(linewidth);
	hMID400iso->SetLineColor(kBlue);
	hMID400iso->SetFillColor(kBlue);
//	hMID400iso->GetXaxis()->SetRangeUser(0.,50.);
	
	hMID15iso ->SetMarkerColor(kGreen);
	hMID15iso ->SetMarkerStyle(20);
	hMID15iso ->SetMarkerSize(1.);
	hMID15iso ->SetLineWidth(linewidth);
	hMID15iso ->SetLineColor(kGreen);
	hMID15iso ->SetFillColor(kGreen);
//	hMID15iso ->GetXaxis()->SetRangeUser(0.,50.);
	
	TLatex *lat = new TLatex();
	lat->SetNDC(kTRUE);
	lat->SetTextColor(kBlack);
	lat->SetTextSize(0.04);
	
	TLegend *leg;
	leg = new TLegend(0.15,0.75,0.35,0.88);
//	if(fp == SigSup) leg = new TLegend(0.15,0.75,0.35,0.88);
//	if(fp == ZDecay) leg = new TLegend(0.15,0.15,0.35,0.28);
	leg->AddEntry(hMID24iso,  "W Mom",             "p");
	leg->AddEntry(hMID500iso, "Bottom Meson Mom",  "p");
	leg->AddEntry(hMID400iso, "Charmed Meson Mom", "p");
	leg->AddEntry(hMID15iso,  "tau Mom",           "p");
	leg->SetTextSize(0.04);
	leg->SetFillStyle(0);
	leg->SetTextFont(42);
	leg->SetBorderSize(0);
	
	TCanvas *c_temp = new TCanvas("C_MIDIsoPlot", "Isolation for different lepton gen moms in " + title, 0, 0, 800, 600);
	c_temp->cd();
	hMID500iso->Draw("PE 0");
	hMID400iso->Draw("PE 0 same");
	hMID15iso ->Draw("PE 0 same");
	hMID24iso ->Draw("PE 0 same");
	leg->Draw();
	//	lat->DrawLatex(0.70,0.92, Form("L_{int.} = %2.1f fb^{-1}", fLumiNorm/1000.));
//	lat->SetTextSize(0.04);
//	lat->DrawLatex(0.62,0.85, title);
//	if(fp == SigSup) lat->DrawLatex(0.62,0.85, pfname + name);
//	if(fp == ZDecay) lat->DrawLatex(0.67,0.15, pfname + name);
//	double ymean(0.), yrms(0.);
//	getWeightedYMeanRMS(h_ptratio_data, ymean, yrms);
	drawTopLine();
	lat->SetTextSize(0.03);
	//	lat->DrawLatex(0.25,0.92, Form("Mean ratio: %4.2f #pm %4.2f", ymean, yrms));
	
//	TString fpname = "F";
//	if(fp == ZDecay) fpname = "P";
	
	// Util::PrintNoEPS( c_temp, fpname + "Ratio_" + name + "_Pt", fOutputDir + fOutputSubDir, NULL);
//	Util::PrintPDF(c_temp, fpname + "Ratio_" + name + "_Pt", fOutputDir + fOutputSubDir);
	Util::PrintPDF(c_temp, title, fOutputDir + fOutputSubDir);
//	delete h_ptratio_mc, h_ptratio_data, h_ptratio_qcd, h_ptratio_wjets;
	delete c_temp, lat, leg;
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

	TGraphAsymmErrors* gr_obs = FR->getGraphPoissonErrors( h_obs );
	gr_obs->SetMarkerColor(kBlack);
	gr_obs->SetMarkerStyle(20);
	gr_obs->SetMarkerSize(2.5);
	gr_obs->SetLineWidth(2);
	gr_obs->SetLineColor(kBlack);
	gr_obs->SetFillColor(kBlack);
	
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
	gr_obs->Draw("P same");
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
	delete gr_obs;
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
	Init(tree);
	fSample = fSamples[sample];
	for (Long64_t jentry=0; jentry<tree->GetEntriesFast();jentry++) {
		tree->GetEntry(jentry);
		printProgress(jentry, tree->GetEntriesFast(), fSamples[sample]->name);
		if(singleMuTrigger() == false) continue;
		int mu1(-1), mu2(-1);
		if(hasLooseMuons(mu1, mu2) < 1) continue;
		setHypLepton1(mu1, Muon);
		if(!passesJet50Cut())           continue;
		if(getNJets() < 1)              continue;
		if(MuMT[mu1] > 20.)             continue;
		if(NMus > 1)                    continue;
		for(size_t j = 0; j < nmetbins; ++j) if(pfMET > metbins[j] && pfMET < metbins[j+1]) hiso[j]->Fill(MuPFIso[mu1], singleMuPrescale()); // MARC
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
		fH2D_MufRatio = fillRatio(  Muon, musamples, SigSup, false);
		fH2D_MupRatio = fillRatio(  Muon, musamples, ZDecay, false);
		fH2D_ElfRatio = fillRatio(  Elec, elsamples, SigSup, false);
		fH2D_ElpRatio = fillRatio(  Elec, elsamples, ZDecay, false);
	}
	if(datamc == 1){
		fH1D_MufRatio_MC = fillRatioPt(Muon, musamples, SigSup, false);
		fH1D_MupRatio_MC = fillRatioPt(Muon, musamples, ZDecay, false);
		fH1D_ElfRatio_MC = fillRatioPt(Elec, elsamples, SigSup, false);
		fH1D_ElpRatio_MC = fillRatioPt(Elec, elsamples, ZDecay, false);
		fH2D_MufRatio_MC = fillRatio(  Muon, musamples, SigSup, false);
		fH2D_MupRatio_MC = fillRatio(  Muon, musamples, ZDecay, false);
		fH2D_ElfRatio_MC = fillRatio(  Elec, elsamples, SigSup, false);
		fH2D_ElpRatio_MC = fillRatio(  Elec, elsamples, ZDecay, false);
	}
}
void SSDLPlotter::fillRatios(vector<int> frmusamples, vector<int> frelsamples, vector<int> prmusamples, vector<int> prelsamples, int datamc){
	if(datamc == 0){
		fH1D_MufRatio = fillRatioPt(Muon, frmusamples, SigSup, false);
		fH1D_MupRatio = fillRatioPt(Muon, prmusamples, ZDecay, false);
		fH1D_ElfRatio = fillRatioPt(Elec, frelsamples, SigSup, false);
		fH1D_ElpRatio = fillRatioPt(Elec, prelsamples, ZDecay, false);
		fH2D_MufRatio = fillRatio(  Muon, frmusamples, SigSup, false);
		fH2D_MupRatio = fillRatio(  Muon, prmusamples, ZDecay, false);
		fH2D_ElfRatio = fillRatio(  Elec, frelsamples, SigSup, false);
		fH2D_ElpRatio = fillRatio(  Elec, prelsamples, ZDecay, false);
	}
	if(datamc == 1){
		fH1D_MufRatio_MC = fillRatioPt(Muon, frmusamples, SigSup, false);
		fH1D_MupRatio_MC = fillRatioPt(Muon, prmusamples, ZDecay, false);
		fH1D_ElfRatio_MC = fillRatioPt(Elec, frelsamples, SigSup, false);
		fH1D_ElpRatio_MC = fillRatioPt(Elec, prelsamples, ZDecay, false);
		fH2D_MufRatio_MC = fillRatio(  Muon, frmusamples, SigSup, false);
		fH2D_MupRatio_MC = fillRatio(  Muon, prmusamples, ZDecay, false);
		fH2D_ElfRatio_MC = fillRatio(  Elec, frelsamples, SigSup, false);
		fH2D_ElpRatio_MC = fillRatio(  Elec, prelsamples, ZDecay, false);
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
	
	TH2D *h_2d;
	TH1D *h_pt, *h_eta;
	if(fp == SigSup){
		h_2d  = new TH2D(Form("%sRatio",   shortname[muelswitch].Data()), Form("Ratio of tight to loose %s vs Pt vs Eta", longname[muelswitch].Data()), getNFPtBins(chan), getFPtBins(chan), getNEtaBins(chan), getEtaBins(chan));
		h_pt  = new TH1D(Form("%sRatioPt" ,shortname[muelswitch].Data()), Form("Ratio of tight to loose %s vs Pt"       , longname[muelswitch].Data()), getNFPtBins(chan), getFPtBins(chan));
		h_eta = new TH1D(Form("%sRatioEta",shortname[muelswitch].Data()), Form("Ratio of tight to loose %s vs Eta"      , longname[muelswitch].Data()), getNEtaBins(chan), getEtaBins(chan));
	}
	if(fp == ZDecay){
		h_2d  = new TH2D(Form("%sRatio",   shortname[muelswitch].Data()), Form("Ratio of tight to loose %s vs Pt vs Eta", longname[muelswitch].Data()), getNPPtBins(chan), getPPtBins(chan), getNEtaBins(chan), getEtaBins(chan));
		h_pt  = new TH1D(Form("%sRatioPt" ,shortname[muelswitch].Data()), Form("Ratio of tight to loose %s vs Pt"       , longname[muelswitch].Data()), getNPPtBins(chan), getPPtBins(chan));
		h_eta = new TH1D(Form("%sRatioEta",shortname[muelswitch].Data()), Form("Ratio of tight to loose %s vs Eta"      , longname[muelswitch].Data()), getNEtaBins(chan), getEtaBins(chan));
	}

	h_pt->SetXTitle("p_{#perp} (GeV)");
	h_pt->SetYTitle("# Tight / # Loose");
	h_pt->GetYaxis()->SetTitleOffset(1.2);

	calculateRatio(samples, chan, fp, h_2d, h_pt, h_eta, output);
	delete h_2d, h_eta;
	return h_pt;
}
TH2D* SSDLPlotter::fillRatio(gChannel chan, int sample, gFPSwitch fp, bool output){
	vector<int> samples; samples.push_back(sample);
	return fillRatio(chan, samples, fp);
}
TH2D* SSDLPlotter::fillRatio(gChannel chan, vector<int> samples, gFPSwitch fp, bool output){
	gStyle->SetOptStat(0);
	TString shortname[2] = {"Mu", "El"};
	TString longname[2] = {"Muons", "Electrons"};
	int muelswitch = 0;
	if(chan == Elec) muelswitch = 1;

	TH2D *h_2d;
	TH1D *h_pt, *h_eta;
	if(fp == SigSup){
		h_2d  = new TH2D(Form("%sRatio",   shortname[muelswitch].Data()), Form("Ratio of tight to loose %s vs Pt vs Eta", longname[muelswitch].Data()), getNFPtBins(chan), getFPtBins(chan), getNEtaBins(chan), getEtaBins(chan));
		h_pt  = new TH1D(Form("%sRatioPt" ,shortname[muelswitch].Data()), Form("Ratio of tight to loose %s vs Pt"       , longname[muelswitch].Data()), getNFPtBins(chan), getFPtBins(chan));
		h_eta = new TH1D(Form("%sRatioEta",shortname[muelswitch].Data()), Form("Ratio of tight to loose %s vs Eta"      , longname[muelswitch].Data()), getNEtaBins(chan), getEtaBins(chan));
	}
	if(fp == ZDecay){
		h_2d  = new TH2D(Form("%sRatio",   shortname[muelswitch].Data()), Form("Ratio of tight to loose %s vs Pt vs Eta", longname[muelswitch].Data()), getNPPtBins(chan), getPPtBins(chan), getNEtaBins(chan), getEtaBins(chan));
		h_pt  = new TH1D(Form("%sRatioPt" ,shortname[muelswitch].Data()), Form("Ratio of tight to loose %s vs Pt"       , longname[muelswitch].Data()), getNPPtBins(chan), getPPtBins(chan));
		h_eta = new TH1D(Form("%sRatioEta",shortname[muelswitch].Data()), Form("Ratio of tight to loose %s vs Eta"      , longname[muelswitch].Data()), getNEtaBins(chan), getEtaBins(chan));
	}

	h_2d->SetXTitle("p_{#perp} (GeV)");
	h_2d->SetYTitle("#eta");
	h_2d->SetZTitle("# Tight / # Loose");

	calculateRatio(samples, chan, fp, h_2d, h_pt, h_eta, output);
	delete h_pt, h_eta;
	return h_2d;
}

//____________________________________________________________________________
void SSDLPlotter::calculateRatio(vector<int> samples, gChannel chan, gFPSwitch fp, TH2D*& h_2d, bool output, bool ttbarMatched){
	TH1D *h_dummy1 = new TH1D("dummy1", "dummy1", 1, 0.,1.);
	TH1D *h_dummy2 = new TH1D("dummy2", "dummy2", 1, 0.,1.);
	calculateRatio(samples, chan, fp, h_2d, h_dummy1, h_dummy2, output, ttbarMatched);
	delete h_dummy1, h_dummy2;
}
void SSDLPlotter::calculateRatio(vector<int> samples, gChannel chan, gFPSwitch fp, TH2D*& h_2d, TH1D*& h_pt, TH1D*& h_eta, bool output, bool ttbarMatched){
	TH1D *h_dummy1 = new TH1D("rat_dummy1", "rat_dummy1", 18, 0.,36.);
	calculateRatio(samples, chan, fp, h_2d, h_pt, h_eta, h_dummy1, output, ttbarMatched);
	delete h_dummy1;
}
void SSDLPlotter::calculateRatio(vector<int> samples, gChannel chan, gFPSwitch fp, TH2D*& h_2d, TH1D*& h_pt, TH1D*& h_eta, TH1D*& h_nv, bool output, bool ttbarMatched){
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

	TH1D *H_ntight_nv = new TH1D("NTight_NV", "NTight Muons NV", 18, 0., 36.); H_ntight_nv->Sumw2();
	TH1D *H_nloose_nv = new TH1D("NLoose_NV", "NLoose Muons NV", 18, 0., 36.); H_nloose_nv->Sumw2();
	if (gRatiosFromTTbar || ttbarMatched) getPassedTotalTTbar(samples, chan, fp, H_ntight, H_nloose, output);
	// marc jan 27 else getPassedTotal(samples, chan, fp, H_ntight, H_nloose, output);
	else getPassedTotal(samples, chan, fp, H_ntight, H_nloose, H_ntight_nv, H_nloose_nv, output);
	h_2d->Divide(H_ntight,    H_nloose,    1., 1., "B");
	h_nv->Divide(H_ntight_nv, H_nloose_nv, 1., 1., "B");

	TH1D *hmuloosept  = H_nloose->ProjectionX();
	TH1D *hmulooseeta = H_nloose->ProjectionY();
	TH1D *hmutightpt  = H_ntight->ProjectionX();
	TH1D *hmutighteta = H_ntight->ProjectionY();

	h_pt ->Divide(hmutightpt,  hmuloosept,  1., 1., "B"); // binomial
	h_eta->Divide(hmutighteta, hmulooseeta, 1., 1., "B"); // weights are ignored
	delete H_ntight, H_nloose, H_ntight_nv, H_nloose_nv, hmuloosept, hmulooseeta, hmutightpt, hmutighteta;
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

		float scale = fLumiNorm/S->getLumi(); // Normalize all
		if(S->datamc == 0) scale = 1;
		if (gRatiosFromTTbar) {
			Channel *C;
			if(chan == Muon) C = &S->region[gRegion[gBaseRegion]][HighPt].mm;
			if(chan == Elec) C = &S->region[gRegion[gBaseRegion]][HighPt].ee;
			if(fp == SigSup){
				ntight = C->fntight_ttbar->GetEntries();
				nloose = C->fnloose_ttbar->GetEntries();
			}
			if(fp == ZDecay){
				ntight = C->pntight_ttbar->GetEntries();
				nloose = C->pnloose_ttbar->GetEntries();
			}
		}
		else {
			if(fp == SigSup){
				// old: pay attention here. does this make sense to change to the baseregion?? it was Baseline before...
				ntight += scale * S->numbers[gRegion[gBaseRegion]][chan].nsst;
				nloose += scale * S->numbers[gRegion[gBaseRegion]][chan].nssl;
				
				ntight_sam += S->numbers[gRegion[gBaseRegion]][chan].nsst;
				nloose_sam += S->numbers[gRegion[gBaseRegion]][chan].nssl;
			}
			if(fp == ZDecay){
				ntight += scale * S->numbers[gRegion[gBaseRegion]][chan].nzt;
				nloose += scale * S->numbers[gRegion[gBaseRegion]][chan].nzl;
				
				ntight_sam += S->numbers[gRegion[gBaseRegion]][chan].nzt;
				nloose_sam += S->numbers[gRegion[gBaseRegion]][chan].nzl;
			}
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

		float scale = fLumiNorm/S->getLumi(); // Normalize all
		if(S->datamc == 0) scale = 1;
		if(fp == SigSup){
			// old: pay attention here. does this make sense to change to the baseregion?? it was Baseline before...
			ntight += scale * S->numbers[gRegion[gBaseRegion]][chan].nsst;
			nloose += scale * S->numbers[gRegion[gBaseRegion]][chan].nssl;

			ntight_sam += S->numbers[gRegion[gBaseRegion]][chan].nsst;
			nloose_sam += S->numbers[gRegion[gBaseRegion]][chan].nssl;
		}
		if(fp == ZDecay){
			ntight += S->numbers[gRegion[gBaseRegion]][chan].nzt;
			nloose += S->numbers[gRegion[gBaseRegion]][chan].nzl;

			ntight_sam += S->numbers[gRegion[gBaseRegion]][chan].nzt;
			nloose_sam += S->numbers[gRegion[gBaseRegion]][chan].nzl;
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
TEfficiency *SSDLPlotter::getMergedEfficiency(vector<int> samples, gChannel chan, gFPSwitch fp, int pteta){
	// Only call for data! For MC call the combination one
	// pteta switch: 0: pt (default), 1: eta
	TEfficiency *eff;
	if(fp == SigSup && pteta == 0) eff = new TEfficiency("fRatio_pt",  "fRatio_pt",  getNFPtBins(chan), getFPtBins(chan));
	if(fp == ZDecay && pteta == 0) eff = new TEfficiency("fRatio_pt",  "fRatio_pt",  getNPPtBins(chan), getPPtBins(chan));
	if(fp == SigSup && pteta == 1) eff = new TEfficiency("fRatio_eta", "fRatio_eta", getNEtaBins(chan), getEtaBins(chan));
	if(fp == ZDecay && pteta == 1) eff = new TEfficiency("fRatio_eta", "fRatio_eta", getNEtaBins(chan), getEtaBins(chan));
	if(fp == SigSup && pteta == 2) eff = new TEfficiency("fRatio_nv",  "fRatio_nv",  18, 0., 36.);
	if(fp == ZDecay && pteta == 2) eff = new TEfficiency("fRatio_nv",  "fRatio_nv",  18, 0., 36.);
	for(size_t i = 0; i < samples.size(); ++i){
		Sample *S = fSamples[samples[i]];

		if(S->datamc > 0) { cout << "Calling the wrong method -> call the combined one" << endl; exit(-1);}
		Channel *C;
			// old: pay attention here. does this make sense to change to the baseregion?? it was Baseline before...
		if(chan == Muon) C = &S->region[gRegion[gBaseRegion]][HighPt].mm;
		if(chan == Elec) C = &S->region[gRegion[gBaseRegion]][HighPt].ee;
		TEfficiency *tempeff;
		if(fp == SigSup && pteta == 0) tempeff = new TEfficiency(*C->fratio_pt);
		if(fp == SigSup && pteta == 1) tempeff = new TEfficiency(*C->fratio_eta);
		if(fp == SigSup && pteta == 2) tempeff = new TEfficiency(*C->fratio_nv);
		if(fp == ZDecay && pteta == 0) tempeff = new TEfficiency(*C->pratio_pt);
		if(fp == ZDecay && pteta == 1) tempeff = new TEfficiency(*C->pratio_eta);
		if(fp == ZDecay && pteta == 2) tempeff = new TEfficiency(*C->pratio_nv);

		eff->Add(*tempeff);
		delete tempeff;
	}
	return eff;
}
TGraphAsymmErrors *SSDLPlotter::getCombEfficiency(vector<int> samples, gChannel chan, gFPSwitch fp, int pteta){
	// pteta switch: 0: pt (default), 1: eta
	TList *list = new TList();
	vector<TEfficiency*> veff;
	for(size_t i = 0; i < samples.size(); ++i){
		Sample *S = fSamples[samples[i]];
		float scale = fLumiNorm / S->getLumi();

		Channel *C;
			// old: pay attention here. does this make sense to change to the baseregion?? it was Baseline before...
		if(chan == Muon) C = &S->region[gRegion[gBaseRegion]][HighPt].mm;
		if(chan == Elec) C = &S->region[gRegion[gBaseRegion]][HighPt].ee;
		if(fp == SigSup && pteta == 0) veff.push_back(new TEfficiency(*C->fratio_pt));
		if(fp == SigSup && pteta == 1) veff.push_back(new TEfficiency(*C->fratio_eta));
		if(fp == ZDecay && pteta == 0) veff.push_back(new TEfficiency(*C->pratio_pt));
		if(fp == ZDecay && pteta == 1) veff.push_back(new TEfficiency(*C->pratio_eta));
		cout << veff[i]->GetName() << endl;
		veff[i]->SetWeight(scale);
		list->Add(veff[i]);
	}
	TEfficiency *eff = new TEfficiency();
	TGraphAsymmErrors *asym = eff->Combine(list);
	for(size_t i = 0; i < veff.size(); ++i) delete veff[i];
	return asym;
}

void SSDLPlotter::getPassedTotal(vector<int> samples, gChannel chan, gFPSwitch fp, TH2D*& h_passed, TH2D*& h_total, TH1D*& h_passed_nv, TH1D*& h_total_nv, bool output){
	if(fVerbose>2) cout << "---------------" << endl;
	for(size_t i = 0; i < samples.size(); ++i){
		Sample *S = fSamples[samples[i]];

		float scale = fLumiNorm / S->getLumi();
		if(S->datamc == 0) scale = 1;

		Channel *C;
		// old: pay attention here. does this make sense to change to the baseregion?? it was Baseline before...
		if(chan == Muon) C = &S->region[gRegion[gBaseRegion]][HighPt].mm;
		if(chan == Elec) C = &S->region[gRegion[gBaseRegion]][HighPt].ee;
		TH2D *ntight, *nloose;
		TH1D *ntight_nv, *nloose_nv;
		if(fp == SigSup){
			ntight = C->fntight;
			nloose = C->fnloose;
			ntight_nv = C->fntight_nv;
			nloose_nv = C->fnloose_nv;
		} else if(fp == ZDecay){
			ntight = C->pntight;
			nloose = C->pnloose;
			ntight_nv = C->pntight_nv;
			nloose_nv = C->pnloose_nv;
		}
		h_passed->Add(ntight, scale);
		h_total ->Add(nloose, scale);
		h_passed_nv->Add(ntight_nv, scale);
		h_total_nv ->Add(nloose_nv, scale);
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
void SSDLPlotter::getPassedTotalTTbar(vector<int> samples, gChannel chan, gFPSwitch fp, TH2D*& h_passed, TH2D*& h_total, bool output){
	if(fVerbose>2) cout << "---------------" << endl;
	for(size_t i = 0; i < samples.size(); ++i){
		Sample *S = fSamples[samples[i]];

		float scale = fLumiNorm / S->getLumi();
		if(S->datamc == 0) scale = 1; // this function doesn't really make sense for data, since the histograms aren't filled.

		Channel *C;
		// old: pay attention here. does this make sense to change to the baseregion?? it was Baseline before...
		if(chan == Muon) C = &S->region[gRegion[gBaseRegion]][HighPt].mm;
		if(chan == Elec) C = &S->region[gRegion[gBaseRegion]][HighPt].ee;
		TH2D *ntight, *nloose;
		if(fp == SigSup){
			ntight = C->fntight_ttbar;
			nloose = C->fnloose_ttbar;
		} else if(fp == ZDecay){
			ntight = C->pntight_ttbar;
			nloose = C->pnloose_ttbar;
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
		printObject(h_passed, TString("Passed") + name + TString("_onlyTTbar"), "colz");
		printObject(h_total,  TString("Total")  + name + TString("_onlyTTbar"), "colz");
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

		float scale = fLumiNorm / S->getLumi();
		if(S->datamc == 0) scale = 1;

		FRatioPlots *RP;
		if(chan == Muon) RP = &S->ratioplots[0];
		if(chan == Elec) RP = &S->ratioplots[1];
		ntight->Add(RP->ntight[ratiovar], scale);
		nloose->Add(RP->nloose[ratiovar], scale);
	}
	
	ratio->Divide(ntight, nloose, 1., 1., "B");

	delete ntight, nloose;
	return ratio;
}
//____________________________________________________________________________
void SSDLPlotter::calculateChMisIdProb(vector<int>  samples, gChMisIdReg chmid_reg, float &chmid, float &chmide){
  
  float ospair(0.),   sspair(0.);
  float ospair_e(0.), sspair_e(0.);
  for (size_t i=0; i<samples.size(); ++i){
    Sample *S = fSamples[samples[i]];

    float scale = fLumiNorm/S->getLumi(); // Normalize all
    if(S->datamc == 0) scale = 1;    

	// old: pay attention here. does this make sense to change to the baseregion?? it was Baseline before...
    Channel *C = &S->region[gRegion[gBaseRegion]][HighPt].ee;
    TH2D *ospairstmp = C->ospairs;
    TH2D *sspairstmp = C->sspairs;
    
    if (chmid_reg == BB){
      ospair   += ospairstmp->GetBinContent(1,1) * scale;
      sspair   += sspairstmp->GetBinContent(1,1) * scale;
      ospair_e += ospairstmp->GetBinError(1,1)   * scale;
      sspair_e += sspairstmp->GetBinError(1,1)   * scale;
    }
    if (chmid_reg == EB){
      ospair   += (ospairstmp->GetBinContent(1,2) + ospairstmp->GetBinContent(2,1)) * scale;
      sspair   += (sspairstmp->GetBinContent(1,2) + sspairstmp->GetBinContent(2,1)) * scale;
      ospair_e += (ospairstmp->GetBinError(1,2)   + ospairstmp->GetBinError(2,1)  ) * scale;
      sspair_e += (sspairstmp->GetBinError(1,2)   + sspairstmp->GetBinError(2,1)  ) * scale;     
    }
    if (chmid_reg == EE){
      ospair   += ospairstmp->GetBinContent(2,2) * scale;
      sspair   += sspairstmp->GetBinContent(2,2) * scale;
      ospair_e += ospairstmp->GetBinError(2,2)   * scale;
      sspair_e += sspairstmp->GetBinError(2,2)   * scale;
    }
  }
  
  ratioWithBinomErrors(sspair, ospair, chmid, chmide);
  
  // Divide to get the per-electron probability...
  chmid  = chmid  / 2.;
  chmide = chmide / 2.;
  
  return;
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
	TString notetable     = outputdir + "NoteTable.tex";
	TString combotool     = outputdir + "eth_signalRegions_highpt.txt";

	fOUTSTREAM.open(tablefilename.Data(), ios::trunc);
	fOUTSTREAM << "==========================================================================================================" << endl;
	fOUTSTREAM << " Table 2 inputs from ETH Analysis" << endl;
	fOUTSTREAM << Form(" Generated on: %s ", asctime(timeinfo)) << endl;
	fOUTSTREAM << endl;
	
	fOUTSTREAM3.open(notetable.Data(), ios::trunc);
	fOUTSTREAM3 << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
	fOUTSTREAM3 << Form("%%%% Generated on: %s ", asctime(timeinfo)) << endl;
	fOUTSTREAM3 << "%% Format is tot, (ee, mm, em)" << endl;
	fOUTSTREAM3 << endl;
	
	fOUTSTREAM4.open(combotool.Data(), ios::trunc);
	fOUTSTREAM4 << "{" << endl;

	vector<SSPrediction> sspreds;
 	for(size_t i = 0; i < gNREGIONS; ++i){
	  TString outputname = outputdir + "DataPred_" + gRegions[i]->sname + ".txt";
	  sspreds.push_back(makeIntPrediction(outputname, i));
 	}

	fOUTSTREAM.close();
	fOUTSTREAM2.close();
	fOUTSTREAM3.close();
	fOUTSTREAM4.open(combotool.Data(), ios::trunc);
	fOUTSTREAM4 << "{" << endl;

	
	if (gDoWZValidation) makeWZValidation(sspreds.at(gRegion["WZEnriched"]));
	
	// Produce Datacards only for certain regions.
	vector<int> ewkregions;
	ewkregions.push_back(gRegion["HT0MET120NJ2bV"]);
	ewkregions.push_back(gRegion["HT0MET200"]);
	ewkregions.push_back(gRegion["HT0MET120NJ2bVlV"]);
	ewkregions.push_back(gRegion["HT0MET200lV"]);
	
	for (size_t i = 0; i<ewkregions.size(); i++) {
	  int reg = ewkregions.at(i);
	  // old gRegion reg = gRegion(ewkregions.at(i));
	  TString datacard = outputdir + "datacard_"+ gRegions[reg]->sname + ".txt";
	  fOUTSTREAM.open(datacard.Data(), ios::trunc);
	  fOUTSTREAM <<      "#=========================================================================================" << endl;
	  fOUTSTREAM <<      "# Systematics table for EWKino analysis, same-sign channel" << endl;
	  fOUTSTREAM << Form("# Generated on: %s ", asctime(timeinfo)) << endl;
	  fOUTSTREAM << Form("# Generated for region: %s ", gRegions[reg]->sname.Data()) << endl;
	  fOUTSTREAM <<      "# Copy between the dashed lines for datacard" << endl;
	  fOUTSTREAM <<      "#-----------------------------------------------------------------------------------------" << endl;
	  fOUTSTREAM <<      "imax 1" << endl;
	  fOUTSTREAM <<      "jmax *" << endl;
	  fOUTSTREAM <<      "kmax *" << endl;
	  fOUTSTREAM << endl << endl;
	  fOUTSTREAM <<      "bin\t\t1" << endl;
	  fOUTSTREAM << Form("observation\t%d", sspreds[reg].obs) << endl;
	  fOUTSTREAM << endl << endl;
	  fOUTSTREAM <<      "bin\t\t1\t\t1\t\t1\t\t1\t\t1" << endl;
	  fOUTSTREAM <<      "process\t\tsig\t\tfake\t\tcmid\t\twz\t\trare" << endl;
	  fOUTSTREAM <<      "process\t\t0\t\t1\t\t2\t\t3\t\t4" << endl;
	  fOUTSTREAM << Form("rate\t\t%5.3f\t\t%5.3f\t\t%5.3f\t\t%5.3f\t\t%5.3f",
			     sspreds[reg].s, sspreds[reg].fake, sspreds[reg].cmid, sspreds[reg].wz, sspreds[reg].rare) << endl;
	  fOUTSTREAM << endl << endl;
	  fOUTSTREAM <<      "lumi     lnN\t1.05\t\t1.05\t\t1.05\t\t1.05\t\t1.05" << endl;
	  //	  fOUTSTREAM << Form("bgUncsf  lnN\t-\t\t%5.3f\t\t-\t\t-\t\t-\t\t-", 1.0+sspreds[reg].df_err/fabs(sspreds[reg].df))     << endl;
	  fOUTSTREAM << Form("sigUnc   lnN\t%5.3f\t\t-\t\t-\t\t-\t\t-", 1.0+0.3/sspreds[reg].s)                       << endl;
	  fOUTSTREAM << Form("bgUncdf  lnN\t-\t\t%5.3f\t\t-\t\t-\t\t-", 1.0+sspreds[reg].fake_err/sspreds[reg].fake) << endl;
	  fOUTSTREAM << Form("bgUnccmi lnN\t-\t\t-\t\t%5.3f\t\t-\t\t-", 1.0+sspreds[reg].cmid_err/sspreds[reg].cmid) << endl;
	  fOUTSTREAM << Form("bgUncwz  lnN\t-\t\t-\t\t-\t\t%5.3f\t\t-", 1.0+sspreds[reg].wz_err  /sspreds[reg].wz)   << endl;
	  fOUTSTREAM << Form("bgUncrar lnN\t-\t\t-\t\t-\t\t-\t\t%5.3f", 1.0+sspreds[reg].rare_err/sspreds[reg].rare) << endl;
	  fOUTSTREAM << endl;
	  fOUTSTREAM.close();
	}

 	TString datacard = outputdir + "EWKdatacard.txt";
 	fOUTSTREAM.open(datacard.Data(), ios::trunc);
 	fOUTSTREAM <<      "#=========================================================================================" << endl;
 	fOUTSTREAM <<      "# Systematics table for EWKino analysis, same-sign channel" << endl;
 	fOUTSTREAM << Form("# Generated on: %s ", asctime(timeinfo)) << endl;
 	fOUTSTREAM <<      "#-----------------------------------------------------------------------------------------" << endl;
	fOUTSTREAM <<      "#<channel label> <observed> <total background> <total bkg uncertainty> <stat uncertainty> <wz> <wz uncertainty> <ttbar+fakes> <ttbar+fakes uncertainty> <zgamma> <zgamma uncertainty> <zz> <z uncertainty><rare> <rare uncertainty>" << endl;
	for (size_t i = 0; i<ewkregions.size(); i++) {
 	  int reg = ewkregions.at(i);
 	  // old gRegion reg = gRegion(ewkregions.at(i));
 	  fOUTSTREAM << Form("# Generated for region: %s ", gRegions[reg]->sname.Data()) << endl;
	  fOUTSTREAM << Form("data %d %d %6.4f %6.4f %6.4f ",i, sspreds[reg].obs, 
			     sspreds[reg].bkg, sspreds[reg].bkg_err, sspreds[reg].bkg_stat);
	  fOUTSTREAM << Form("%6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f",
			     sspreds[reg].wz, sspreds[reg].wz_err, sspreds[reg].fake, sspreds[reg].fake_err,
			     sspreds[reg].cmid,sspreds[reg].cmid_err, 0., 0., sspreds[reg].rare, sspreds[reg].rare_err) << endl;
	  fOUTSTREAM << "#=====================================================================================================" << endl;
 	}

	// PRODUCE SUMMARY PLOTS ONLY FOR CERTAIN REGIONS
	TH1D    *h_obs        = new TH1D("h_observed",   "Observed number of events",  4, 0., 4.);
	TH1D    *h_pred_sfake = new TH1D("h_pred_sfake", "Predicted single fakes", 4, 0., 4.);
	TH1D    *h_pred_dfake = new TH1D("h_pred_dfake", "Predicted double fakes", 4, 0., 4.);
	TH1D    *h_pred_chmid = new TH1D("h_pred_chmid", "Predicted charge mis id", 4, 0., 4.);
	TH1D    *h_pred_mc    = new TH1D("h_pred_mc",    "Predicted Rare SM", 4, 0., 4.);
	TH1D    *h_pred_ttw   = new TH1D("h_pred_ttw",   "Predicted WZ", 4, 0., 4.);
	TH1D    *h_pred_tot   = new TH1D("h_pred_tot",   "Total Prediction", 4, 0., 4.);
	THStack *hs_pred      = new THStack("hs_predicted", "Predicted number of events");
	
	h_pred_tot->Sumw2();
	for (size_t i = 0; i<ewkregions.size(); i++) {
	  int reg = ewkregions.at(i);
 
	  h_obs       ->SetBinContent(i+1, sspreds[reg].obs );
	  h_pred_sfake->SetBinContent(i+1, sspreds[reg].sf  );
	  h_pred_dfake->SetBinContent(i+1, sspreds[reg].df  );
	  h_pred_chmid->SetBinContent(i+1, sspreds[reg].cmid);
	  h_pred_mc   ->SetBinContent(i+1, sspreds[reg].rare);
	  h_pred_ttw  ->SetBinContent(i+1, sspreds[reg].wz  );
	  
	  h_pred_sfake->SetBinError(i+1, sspreds[reg].sf_err  );
	  h_pred_dfake->SetBinError(i+1, sspreds[reg].df_err  );
	  h_pred_chmid->SetBinError(i+1, sspreds[reg].cmid_err);
	  h_pred_mc   ->SetBinError(i+1, sspreds[reg].rare_err);
	  h_pred_ttw  ->SetBinError(i+1, sspreds[reg].wz_err  );
	}
	FakeRatios *FR = new FakeRatios();
	FR->setNToyMCs(100);
	FR->setAddESyst(0.5);

	TGraphAsymmErrors* gr_obs = FR->getGraphPoissonErrors( h_obs );
	gr_obs->SetMarkerColor(kBlack);
	gr_obs->SetMarkerStyle(20);
	gr_obs->SetMarkerSize(2);
	gr_obs->SetLineWidth(2);
	gr_obs->SetLineColor(kBlack);
	gr_obs->SetFillColor(kBlack);
	
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
 	h_pred_ttw  ->SetLineWidth(1);
 	h_pred_sfake->SetLineColor(50);
 	h_pred_sfake->SetFillColor(50);
 	h_pred_dfake->SetLineColor(38);
 	h_pred_dfake->SetFillColor(38);
 	h_pred_chmid->SetLineColor(42);
 	h_pred_chmid->SetFillColor(42);
 	h_pred_mc   ->SetLineColor(31);
 	h_pred_mc   ->SetFillColor(31);
 	h_pred_ttw  ->SetLineColor(29);
 	h_pred_ttw  ->SetFillColor(29);

 	h_pred_tot  ->SetLineWidth(1);
 	h_pred_tot  ->SetFillColor(12);
 	h_pred_tot  ->SetFillStyle(3005);
		
 	double max = std::max(h_obs->GetBinContent(1), h_obs->GetBinContent(2));
 	max = 1.7*std::max(max, h_obs->GetBinContent(3));
	max = 1.7*std::max(max, h_obs->GetBinContent(4));
	
 	h_obs       ->SetMaximum(max>1?max+1:1.);
 	h_pred_sfake->SetMaximum(max>1?max+1:1.);
 	h_pred_dfake->SetMaximum(max>1?max+1:1.);
 	h_pred_chmid->SetMaximum(max>1?max+1:1.);
 	h_pred_mc   ->SetMaximum(max>1?max+1:1.);
 	h_pred_ttw  ->SetMaximum(max>1?max+1:1.);
 	h_pred_tot  ->SetMaximum(max>1?max+1:1.);
 	hs_pred     ->SetMaximum(max>1?max+1:1.); 	

 	h_pred_tot->Add(h_pred_sfake);
 	h_pred_tot->Add(h_pred_dfake);
 	h_pred_tot->Add(h_pred_chmid);
 	h_pred_tot->Add(h_pred_mc);
 	h_pred_tot->Add(h_pred_ttw);
 	
 	hs_pred->Add(h_pred_sfake);
 	hs_pred->Add(h_pred_dfake);
 	hs_pred->Add(h_pred_chmid);
 	hs_pred->Add(h_pred_mc);
 	hs_pred->Add(h_pred_ttw);
 	
 	hs_pred->Draw("goff");
	
	TString hadString;
	for (size_t i = 0; i<ewkregions.size(); i++) {
	  int reg = ewkregions.at(i);
// 	  hadString = "#splitline{";
// 	  if (gRegions[reg]->maxMet < 7900) hadString += Form("%.0f < E_{T}^{miss} < %.0f GeV", gRegions[reg]->minMet, gRegions[reg]->maxMet);
// 	  else                              hadString += Form("E_{T}^{miss} > %.0f GeV", gRegions[reg]->minMet);
// 	  hadString += "}{#splitline{";
	  
// 	  if (gRegions[reg]->maxNjets < 99) hadString += Form("N_{jets} #leq %1d,  N_{bjets} = %1d", gRegions[reg]->maxNjets, gRegions[reg]->maxNbjets);
// 	  hadString += "}{";
	  
// 	  if (gRegions[reg]->app3rdVet != 0)hadString += Form("3^{rd} lepton Veto");
// 	  hadString += "}}";
	  
	  // first line
	  hadString = "#splitline{";
	  if (gRegions[reg]->maxMet < 7900) hadString += Form("%.0f < E_{T}^{miss} < %.0f GeV", gRegions[reg]->minMet, gRegions[reg]->maxMet);
	  else                              hadString += Form("E_{T}^{miss} > %.0f GeV", gRegions[reg]->minMet);
	  hadString += "}{#splitline{";
	  
	  // second line
	  if (gRegions[reg]->maxNjets < 99 && gRegions[reg]->app3rdVet !=0) {
	    hadString += Form("N_{jets} #leq %1d,  N_{bjets} = %1d #color[0]{3^{rd}}", gRegions[reg]->maxNjets, gRegions[reg]->maxNbjets);
	    hadString += Form("}{3^{rd} lepton Veto}}");
	  }
	  else if (gRegions[reg]->maxNjets < 99 && gRegions[reg]->app3rdVet == 0) {
	    hadString += Form("N_{jets} #leq %1d,  N_{bjets} = %1d #color[0]{3^{rd}}", gRegions[reg]->maxNjets, gRegions[reg]->maxNbjets);
	    hadString += Form("}{#color[0]{3^{rd} lepton Veto N_{j}}}}");
	  }
	  else if (gRegions[reg]->maxNjets > 90 && gRegions[reg]->app3rdVet !=0){
	    hadString += Form("3^{rd} lepton Veto #color[0]{N_{j}}");
	    hadString += Form("}{#color[0]{3^{rd} lepton Veto N_{j}}}}");
	  }
	  else if (gRegions[reg]->maxNjets > 90 && gRegions[reg]->app3rdVet ==0){
	    hadString += Form("#color[0]{3^{rd} lepton Veto N_{j}}");
	    hadString += Form("}{#color[0]{3^{rd} lepton Veto N_{j}}}}");
	  }
	  hs_pred->GetXaxis()->SetBinLabel(i+1,hadString);
	}
 	hs_pred->GetXaxis()->SetLabelOffset(0.005);
 	hs_pred->GetXaxis()->SetLabelFont(42);
 	hs_pred->GetXaxis()->SetLabelSize(0.03);
 	
 	TLegend *leg = new TLegend(0.60,0.65,0.95,0.88);
 	leg->AddEntry(h_obs,        "Observed","p");
 	leg->AddEntry(h_pred_sfake, "Single Fakes","f");
 	leg->AddEntry(h_pred_dfake, "Double Fakes","f");
 	leg->AddEntry(h_pred_chmid, "Charge MisID","f");
 	leg->AddEntry(h_pred_mc,    "Irreducible (MC)","f");
 	leg->AddEntry(h_pred_ttw,   "WZ Production","f");
 	leg->AddEntry(h_pred_tot,   "Total Uncertainty","f");
 	leg->SetFillStyle(0);
 	leg->SetTextFont(42);
 	// leg->SetTextSize(0.05);
 	leg->SetBorderSize(0);
	
 	TCanvas *c_temp = new TCanvas("C_ObsPred", "Observed vs Predicted", 0, 0, 600, 700);
 	c_temp->cd();
 
 	hs_pred->Draw("hist");
 	h_pred_tot->DrawCopy("0 E2 same");
	gr_obs->Draw("P same");
 	leg->Draw();
 	
	// 	drawRegionSel(reg);
 	drawTopLine(0.53, 0.8);
	
 	gPad->RedrawAxis();
 	Util::PrintPDF(c_temp,   "SummaryPred", fOutputDir + fOutputSubDir);
 	delete c_temp;	
 	delete h_obs, h_pred_sfake, h_pred_dfake, h_pred_chmid, h_pred_mc, h_pred_ttw, h_pred_tot, hs_pred;
	// 	delete gr_obs;
	// 	delete FR;
	
}

void SSDLPlotter::makeTTWIntPredictions(){
	TString outputdir = Util::MakeOutputDir(fOutputDir + "IntPredictionsTTWZ");
	fOutputSubDir = "IntPredictionsTTWZ/";
	time_t rawtime;
	struct tm* timeinfo;
	time(&rawtime);
	timeinfo = localtime(&rawtime);
	// access a chararray containing the date with asctime(timeinfo)
	
	TString tablefilename = outputdir + "Table2.tex";
	TString notetable     = outputdir + "NoteTable.tex";
	fOUTSTREAM.open(tablefilename.Data(), ios::trunc);
	fOUTSTREAM << "==========================================================================================================" << endl;
	fOUTSTREAM << " Table 2 inputs from ETH Analysis" << endl;
	fOUTSTREAM << Form(" Generated on: %s ", asctime(timeinfo)) << endl;
	fOUTSTREAM << endl;
	
	fOUTSTREAM3.open(notetable.Data(), ios::trunc);
	fOUTSTREAM3 << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
	fOUTSTREAM3 << Form("%%%% Generated on: %s ", asctime(timeinfo)) << endl;
	fOUTSTREAM3 << "%% Format is tot, (ee, mm, em)" << endl;
	fOUTSTREAM3 << endl;
	
	vector<int> ttwregions;
	ttwregions.push_back(gRegion["TTbarWPresel"]);
	ttwregions.push_back(gRegion["TTbarWSel"]);
	ttwregions.push_back(gRegion["TTbarWSelJU"]);
	ttwregions.push_back(gRegion["TTbarWSelJD"]);
	ttwregions.push_back(gRegion["TTbarWSelJS"]);
	ttwregions.push_back(gRegion["TTbarWSelBU"]);
	ttwregions.push_back(gRegion["TTbarWSelBD"]);
	ttwregions.push_back(gRegion["TTbarWSelLU"]);
	ttwregions.push_back(gRegion["TTbarWSelLD"]);
	// MARC ttwregions.push_back(gRegion["TTbarWSelIncl"]);
	
	// ttwregions.push_back(gRegion["TTbarWPreselpp"]);
	// ttwregions.push_back(gRegion["TTbarWSelpp"]);
	vector<TTWZPrediction> ttwzpreds;
	for(size_t i = 0; i < ttwregions.size(); ++i){
		int reg = ttwregions[i];
		TString outputname = outputdir + "DataPred_" + gRegions[reg]->sname + ".txt";
		ttwzpreds.push_back(makeIntPredictionTTW(outputname, reg));
	}
	
	const int inm = 1;
	const int iju = 2;
	const int ijd = 3;
	const int ijs = 4;
	const int ibu = 5;
	const int ibd = 6;
	const int ilu = 7;
	const int ild = 8;
	
	fOUTSTREAM.close();
	fOUTSTREAM2.close();
	fOUTSTREAM3.close();
	
	////////////////////////////////////////////////////
	// Datacards and systematics
	TString datacard;
	const float pu_syst = 1.03;
	const float match_syst_up = 1.01485;
	const float match_syst_dn = 0.997558;
	const float scale_syst_up = 1.02302;
	const float scale_syst_dn = 0.965962;
	
	datacard = outputdir + "datacard_TTWZ.txt";
	fOUTSTREAM.open(datacard.Data(), ios::trunc);
	fOUTSTREAM <<      "#=========================================================================================" << endl;
	fOUTSTREAM <<      "# Systematics table for ttW/Z analysis, same-sign channel" << endl;
	fOUTSTREAM << Form("# Generated on: %s ", asctime(timeinfo)) << endl;
	fOUTSTREAM <<      "# Copy between the dashed lines for datacard" << endl;
	fOUTSTREAM <<      "#-----------------------------------------------------------------------------------------" << endl;
	fOUTSTREAM <<      "imax 1" << endl;
	fOUTSTREAM <<      "jmax 4" << endl;
	fOUTSTREAM <<      "kmax *" << endl;
	fOUTSTREAM << endl << endl;
	fOUTSTREAM <<      "bin\t\t1" << endl;
	fOUTSTREAM << Form("observation\t%d", ttwzpreds[inm].obs) << endl;
	fOUTSTREAM << endl << endl;
	fOUTSTREAM <<      "bin\t\t1\t\t1\t\t1\t\t1\t\t1" << endl;
	fOUTSTREAM <<      "process\t\tttWZ\t\tfake\t\tcmid\t\twz\t\trare" << endl;
	fOUTSTREAM <<      "process\t\t0\t\t1\t\t2\t\t3\t\t4" << endl;
	fOUTSTREAM << Form("rate\t\t%5.3f\t\t%5.3f\t\t%5.3f\t\t%5.3f\t\t%5.3f",
					   ttwzpreds[inm].ttwz, ttwzpreds[inm].fake, ttwzpreds[inm].cmid, ttwzpreds[inm].wz, ttwzpreds[inm].rare) << endl;
	fOUTSTREAM << Form("# separate signal yields: %5.3f (ttW), %5.3f (ttZ)", ttwzpreds[inm].ttw, ttwzpreds[inm].ttz) << endl;
	fOUTSTREAM << endl << endl;
	fOUTSTREAM <<      "#syst" << endl;
	fOUTSTREAM <<      "lumi     lnN\t1.022\t\t1.022\t\t1.022\t\t1.022\t\t1.022" << endl;
	fOUTSTREAM << Form("bgUncfak lnN\t-\t\t%5.3f\t\t-\t\t-\t\t-", 1.0+ttwzpreds[inm].fake_err/ttwzpreds[inm].fake) << endl;
	fOUTSTREAM << Form("bgUnccmi lnN\t-\t\t-\t\t%5.3f\t\t-\t\t-", 1.0+ttwzpreds[inm].cmid_err/ttwzpreds[inm].cmid) << endl;
	fOUTSTREAM << Form("bgUncwz  lnN\t-\t\t-\t\t-\t\t%5.3f\t\t-", 1.0+ttwzpreds[inm].wz_err  /ttwzpreds[inm].wz)   << endl;
	fOUTSTREAM << Form("bgUncrar lnN\t-\t\t-\t\t-\t\t-\t\t%5.3f", 1.0+ttwzpreds[inm].rare_err/ttwzpreds[inm].rare) << endl;
	fOUTSTREAM << Form("lept     lnN\t%5.3f/%5.3f\t-\t\t-\t\t%5.3f/%5.3f\t%5.3f/%5.3f",
					   1.0+(ttwzpreds[ild].ttwz-ttwzpreds[inm].ttwz)/ttwzpreds[inm].ttwz,
					   1.0+(ttwzpreds[ilu].ttwz-ttwzpreds[inm].ttwz)/ttwzpreds[inm].ttwz,
					   1.0+(ttwzpreds[ild].wz  -ttwzpreds[inm].wz  )/ttwzpreds[inm].wz,
					   1.0+(ttwzpreds[ilu].wz  -ttwzpreds[inm].wz  )/ttwzpreds[inm].wz,
					   1.0+(ttwzpreds[ild].rare-ttwzpreds[inm].rare)/ttwzpreds[inm].rare,
					   1.0+(ttwzpreds[ilu].rare-ttwzpreds[inm].rare)/ttwzpreds[inm].rare) << endl;
	fOUTSTREAM << Form("btag     lnN\t%5.3f/%5.3f\t-\t\t-\t\t%5.3f/%5.3f\t%5.3f/%5.3f",
					   1.0+(ttwzpreds[ibd].ttwz-ttwzpreds[inm].ttwz)/ttwzpreds[inm].ttwz,
					   1.0+(ttwzpreds[ibu].ttwz-ttwzpreds[inm].ttwz)/ttwzpreds[inm].ttwz,
					   1.0+(ttwzpreds[ibd].wz  -ttwzpreds[inm].wz  )/ttwzpreds[inm].wz,
					   1.0+(ttwzpreds[ibu].wz  -ttwzpreds[inm].wz  )/ttwzpreds[inm].wz,
					   1.0+(ttwzpreds[ibd].rare-ttwzpreds[inm].rare)/ttwzpreds[inm].rare,
					   1.0+(ttwzpreds[ibu].rare-ttwzpreds[inm].rare)/ttwzpreds[inm].rare) << endl;
	fOUTSTREAM << Form("jes      lnN\t%5.3f/%5.3f\t-\t\t-\t\t%5.3f/%5.3f\t%5.3f/%5.3f",
					   1.0+(ttwzpreds[ijd].ttwz-ttwzpreds[inm].ttwz)/ttwzpreds[inm].ttwz,
					   1.0+(ttwzpreds[iju].ttwz-ttwzpreds[inm].ttwz)/ttwzpreds[inm].ttwz,
					   1.0+(ttwzpreds[ijd].wz  -ttwzpreds[inm].wz  )/ttwzpreds[inm].wz,
					   1.0+(ttwzpreds[iju].wz  -ttwzpreds[inm].wz  )/ttwzpreds[inm].wz,
					   1.0+(ttwzpreds[ijd].rare-ttwzpreds[inm].rare)/ttwzpreds[inm].rare,
					   1.0+(ttwzpreds[iju].rare-ttwzpreds[inm].rare)/ttwzpreds[inm].rare) << endl;
	fOUTSTREAM << Form("jer      lnN\t%5.3f\t\t-\t\t-\t\t%5.3f\t\t%5.3f",
					   1.0+(ttwzpreds[ijs].ttwz-ttwzpreds[inm].ttwz)/ttwzpreds[inm].ttwz,
					   1.0+(ttwzpreds[ijs].wz  -ttwzpreds[inm].wz  )/ttwzpreds[inm].wz,
					   1.0+(ttwzpreds[ijs].rare-ttwzpreds[inm].rare)/ttwzpreds[inm].rare) << endl;
	fOUTSTREAM << Form("pu       lnN\t%5.3f\t\t-\t\t-\t\t%5.3f\t\t%5.3f",
					   pu_syst, pu_syst, pu_syst) << endl;
	fOUTSTREAM << Form("matching lnN\t%5.3f/%5.3f\t-\t\t-\t\t%5.3f/%5.3f\t%5.3f/%5.3f",
					   match_syst_up, match_syst_dn, match_syst_up, match_syst_dn, match_syst_up, match_syst_dn) << endl;
	fOUTSTREAM << Form("scale    lnN\t%5.3f/%5.3f\t-\t\t-\t\t%5.3f/%5.3f\t%5.3f/%5.3f",
					   scale_syst_up, scale_syst_dn, scale_syst_up, scale_syst_dn, scale_syst_up, scale_syst_dn) << endl;
	fOUTSTREAM << endl;
	fOUTSTREAM.close();
	
	datacard = outputdir + "datacard_ssdl_ttW.txt";
	fOUTSTREAM.open(datacard.Data(), ios::trunc);
	fOUTSTREAM <<      "#=========================================================================================" << endl;
	fOUTSTREAM <<      "# Systematics table for ttW/Z analysis, same-sign channel, ttW and ttZ separate" << endl;
	fOUTSTREAM << Form("# Generated on: %s ", asctime(timeinfo)) << endl;
	fOUTSTREAM <<      "# Copy between the dashed lines for datacard" << endl;
	fOUTSTREAM <<      "#-----------------------------------------------------------------------------------------" << endl;
	fOUTSTREAM <<      "imax 1" << endl;
	fOUTSTREAM <<      "jmax 5" << endl;
	fOUTSTREAM <<      "kmax *" << endl;
	fOUTSTREAM << endl << endl;
	fOUTSTREAM <<      "bin\t\t1" << endl;
	fOUTSTREAM << Form("observation\t%d", ttwzpreds[inm].obs) << endl;
	fOUTSTREAM << endl << endl;
	fOUTSTREAM <<      "bin\t\t1\t\t1\t\t1\t\t1\t\t1\t\t1" << endl;
	fOUTSTREAM <<      "process\t\tttW\t\tttZ\t\tfake\t\tcmid\t\twz\t\trare" << endl;
	fOUTSTREAM <<      "process\t\t0\t\t1\t\t2\t\t3\t\t4\t\t5" << endl;
	fOUTSTREAM << Form("rate\t\t%5.3f\t\t%5.3f\t\t%5.3f\t\t%5.3f\t\t%5.3f\t\t%5.3f",
					   ttwzpreds[inm].ttw, ttwzpreds[inm].ttz, ttwzpreds[inm].fake, ttwzpreds[inm].cmid, ttwzpreds[inm].wz, ttwzpreds[inm].rare) << endl;
	fOUTSTREAM << endl << endl;
	fOUTSTREAM <<      "#syst" << endl;
	fOUTSTREAM <<      "lumi     lnN\t1.022\t\t1.022\t\t1.022\t\t1.022\t\t1.022\t\t1.022" << endl;
	fOUTSTREAM << Form("bgUncfak lnN\t-\t\t-\t\t%5.3f\t\t-\t\t-\t\t-", 1.0+ttwzpreds[inm].fake_err/ttwzpreds[inm].fake) << endl;
	fOUTSTREAM << Form("bgUnccmi lnN\t-\t\t-\t\t-\t\t%5.3f\t\t-\t\t-", 1.0+ttwzpreds[inm].cmid_err/ttwzpreds[inm].cmid) << endl;
	fOUTSTREAM << Form("bgUncwz  lnN\t-\t\t-\t\t-\t\t-\t\t%5.3f\t\t-", 1.0+ttwzpreds[inm].wz_err  /ttwzpreds[inm].wz)   << endl;
	fOUTSTREAM << Form("bgUncrar lnN\t-\t\t-\t\t-\t\t-\t\t-\t\t%5.3f", 1.0+ttwzpreds[inm].rare_err/ttwzpreds[inm].rare) << endl;
	fOUTSTREAM << Form("lept     lnN\t%5.3f/%5.3f\t%5.3f/%5.3f\t-\t\t-\t\t%5.3f/%5.3f\t%5.3f/%5.3f",
					   1.0+(ttwzpreds[ild].ttwz-ttwzpreds[inm].ttwz)/ttwzpreds[inm].ttwz,
					   1.0+(ttwzpreds[ilu].ttwz-ttwzpreds[inm].ttwz)/ttwzpreds[inm].ttwz,
					   1.0+(ttwzpreds[ild].ttwz-ttwzpreds[inm].ttwz)/ttwzpreds[inm].ttwz,
					   1.0+(ttwzpreds[ilu].ttwz-ttwzpreds[inm].ttwz)/ttwzpreds[inm].ttwz,
					   1.0+(ttwzpreds[ild].wz  -ttwzpreds[inm].wz  )/ttwzpreds[inm].wz,
					   1.0+(ttwzpreds[ilu].wz  -ttwzpreds[inm].wz  )/ttwzpreds[inm].wz,
					   1.0+(ttwzpreds[ild].rare-ttwzpreds[inm].rare)/ttwzpreds[inm].rare,
					   1.0+(ttwzpreds[ilu].rare-ttwzpreds[inm].rare)/ttwzpreds[inm].rare) << endl;
	fOUTSTREAM << Form("btag     lnN\t%5.3f/%5.3f\t%5.3f/%5.3f\t-\t\t-\t\t%5.3f/%5.3f\t%5.3f/%5.3f",
					   1.0+(ttwzpreds[ibd].ttwz-ttwzpreds[inm].ttwz)/ttwzpreds[inm].ttwz,
					   1.0+(ttwzpreds[ibu].ttwz-ttwzpreds[inm].ttwz)/ttwzpreds[inm].ttwz,
					   1.0+(ttwzpreds[ibd].ttwz-ttwzpreds[inm].ttwz)/ttwzpreds[inm].ttwz,
					   1.0+(ttwzpreds[ibu].ttwz-ttwzpreds[inm].ttwz)/ttwzpreds[inm].ttwz,
					   1.0+(ttwzpreds[ibd].wz  -ttwzpreds[inm].wz  )/ttwzpreds[inm].wz,
					   1.0+(ttwzpreds[ibu].wz  -ttwzpreds[inm].wz  )/ttwzpreds[inm].wz,
					   1.0+(ttwzpreds[ibd].rare-ttwzpreds[inm].rare)/ttwzpreds[inm].rare,
					   1.0+(ttwzpreds[ibu].rare-ttwzpreds[inm].rare)/ttwzpreds[inm].rare) << endl;
	fOUTSTREAM << Form("jes      lnN\t%5.3f/%5.3f\t%5.3f/%5.3f\t-\t\t-\t\t%5.3f/%5.3f\t%5.3f/%5.3f",
					   1.0+(ttwzpreds[ijd].ttwz-ttwzpreds[inm].ttwz)/ttwzpreds[inm].ttwz,
					   1.0+(ttwzpreds[iju].ttwz-ttwzpreds[inm].ttwz)/ttwzpreds[inm].ttwz,
					   1.0+(ttwzpreds[ijd].ttwz-ttwzpreds[inm].ttwz)/ttwzpreds[inm].ttwz,
					   1.0+(ttwzpreds[iju].ttwz-ttwzpreds[inm].ttwz)/ttwzpreds[inm].ttwz,
					   1.0+(ttwzpreds[ijd].wz  -ttwzpreds[inm].wz  )/ttwzpreds[inm].wz,
					   1.0+(ttwzpreds[iju].wz  -ttwzpreds[inm].wz  )/ttwzpreds[inm].wz,
					   1.0+(ttwzpreds[ijd].rare-ttwzpreds[inm].rare)/ttwzpreds[inm].rare,
					   1.0+(ttwzpreds[iju].rare-ttwzpreds[inm].rare)/ttwzpreds[inm].rare) << endl;
	fOUTSTREAM << Form("jer      lnN\t%5.3f\t\t%5.3f\t\t-\t\t-\t\t%5.3f\t\t%5.3f",
					   1.0+(ttwzpreds[ijs].ttwz-ttwzpreds[inm].ttwz)/ttwzpreds[inm].ttwz,
					   1.0+(ttwzpreds[ijs].ttwz-ttwzpreds[inm].ttwz)/ttwzpreds[inm].ttwz,
					   1.0+(ttwzpreds[ijs].wz  -ttwzpreds[inm].wz  )/ttwzpreds[inm].wz,
					   1.0+(ttwzpreds[ijs].rare-ttwzpreds[inm].rare)/ttwzpreds[inm].rare) << endl;
	fOUTSTREAM << Form("pu       lnN\t%5.3f\t\t%5.3f\t\t-\t\t-\t\t%5.3f\t\t%5.3f",
					   pu_syst, pu_syst, pu_syst, pu_syst) << endl;
	fOUTSTREAM << Form("matching lnN\t%5.3f/%5.3f\t%5.3f/%5.3f\t-\t\t-\t\t%5.3f/%5.3f\t%5.3f/%5.3f",
					   match_syst_up, match_syst_dn, match_syst_up, match_syst_dn, match_syst_up, match_syst_dn, match_syst_up, match_syst_dn) << endl;
	fOUTSTREAM << Form("scale    lnN\t%5.3f/%5.3f\t%5.3f/%5.3f\t-\t\t-\t\t%5.3f/%5.3f\t%5.3f/%5.3f",
					   scale_syst_up, scale_syst_dn, scale_syst_up, scale_syst_dn, scale_syst_up, scale_syst_dn, scale_syst_up, scale_syst_dn) << endl;
	fOUTSTREAM << endl;
	fOUTSTREAM.close();
	
	datacard = outputdir + "datacard_ssdl_mm.txt";
	fOUTSTREAM.open(datacard.Data(), ios::trunc);
	fOUTSTREAM <<      "#=========================================================================================" << endl;
	fOUTSTREAM <<      "# Systematics table for ttW/Z analysis, same-sign channel, mumu only" << endl;
	fOUTSTREAM << Form("# Generated on: %s ", asctime(timeinfo)) << endl;
	fOUTSTREAM <<      "# Copy between the dashed lines for datacard" << endl;
	fOUTSTREAM <<      "#-----------------------------------------------------------------------------------------" << endl;
	fOUTSTREAM <<      "imax 1" << endl;
	fOUTSTREAM <<      "jmax 4" << endl;
	fOUTSTREAM <<      "kmax *" << endl;
	fOUTSTREAM << endl << endl;
	fOUTSTREAM <<      "bin\t\t1" << endl;
	fOUTSTREAM << Form("observation\t%d", ttwzpreds[inm].obs_mm) << endl;
	fOUTSTREAM << endl << endl;
	fOUTSTREAM <<      "bin\t\t1\t\t1\t\t1\t\t1\t\t1" << endl;
	fOUTSTREAM <<      "process\t\tttWZ\t\tfake\t\tcmid\t\twz\t\trare" << endl;
	fOUTSTREAM <<      "process\t\t0\t\t1\t\t2\t\t3\t\t4" << endl;
	fOUTSTREAM << Form("rate\t\t%5.3f\t\t%5.3f\t\t%5.3f\t\t%5.3f\t\t%5.3f",
					   ttwzpreds[inm].ttwz_mm, ttwzpreds[inm].fake_mm, 0.0, ttwzpreds[inm].wz_mm, ttwzpreds[inm].rare_mm) << endl;
	fOUTSTREAM << Form("# separate signal yields: %5.3f (ttW), %5.3f (ttZ)", ttwzpreds[inm].ttw_mm, ttwzpreds[inm].ttz_mm) << endl;
	fOUTSTREAM << endl << endl;
	fOUTSTREAM <<      "#syst" << endl;
	fOUTSTREAM <<      "lumi     lnN\t1.022\t\t1.022\t\t1.022\t\t1.022\t\t1.022" << endl;
	fOUTSTREAM << Form("bgUncfak lnN\t-\t\t%5.3f\t\t-\t\t-\t\t-", 1.0+ttwzpreds[inm].fake_err_mm/ttwzpreds[inm].fake_mm) << endl;
	fOUTSTREAM <<      "bgUnccmi lnN\t-\t\t-\t\t-\t\t-\t\t-" << endl;
	fOUTSTREAM << Form("bgUncwz  lnN\t-\t\t-\t\t-\t\t%5.3f\t\t-", 1.0+ttwzpreds[inm].wz_err_mm  /ttwzpreds[inm].wz_mm)   << endl;
	fOUTSTREAM << Form("bgUncrar lnN\t-\t\t-\t\t-\t\t-\t\t%5.3f", 1.0+ttwzpreds[inm].rare_err_mm/ttwzpreds[inm].rare_mm) << endl;
	fOUTSTREAM << Form("lept     lnN\t%5.3f/%5.3f\t-\t\t-\t\t%5.3f/%5.3f\t%5.3f/%5.3f",
					   1.0+(ttwzpreds[ild].ttwz_mm-ttwzpreds[inm].ttwz_mm)/ttwzpreds[inm].ttwz_mm,
					   1.0+(ttwzpreds[ilu].ttwz_mm-ttwzpreds[inm].ttwz_mm)/ttwzpreds[inm].ttwz_mm,
					   1.0+(ttwzpreds[ild].wz_mm  -ttwzpreds[inm].wz_mm  )/ttwzpreds[inm].wz_mm,
					   1.0+(ttwzpreds[ilu].wz_mm  -ttwzpreds[inm].wz_mm  )/ttwzpreds[inm].wz_mm,
					   1.0+(ttwzpreds[ild].rare_mm-ttwzpreds[inm].rare_mm)/ttwzpreds[inm].rare_mm,
					   1.0+(ttwzpreds[ilu].rare_mm-ttwzpreds[inm].rare_mm)/ttwzpreds[inm].rare_mm) << endl;
	fOUTSTREAM << Form("btag     lnN\t%5.3f/%5.3f\t-\t\t-\t\t%5.3f/%5.3f\t%5.3f/%5.3f",
					   1.0+(ttwzpreds[ibd].ttwz_mm-ttwzpreds[inm].ttwz_mm)/ttwzpreds[inm].ttwz_mm,
					   1.0+(ttwzpreds[ibu].ttwz_mm-ttwzpreds[inm].ttwz_mm)/ttwzpreds[inm].ttwz_mm,
					   1.0+(ttwzpreds[ibd].wz_mm  -ttwzpreds[inm].wz_mm  )/ttwzpreds[inm].wz_mm,
					   1.0+(ttwzpreds[ibu].wz_mm  -ttwzpreds[inm].wz_mm  )/ttwzpreds[inm].wz_mm,
					   1.0+(ttwzpreds[ibd].rare_mm-ttwzpreds[inm].rare_mm)/ttwzpreds[inm].rare_mm,
					   1.0+(ttwzpreds[ibu].rare_mm-ttwzpreds[inm].rare_mm)/ttwzpreds[inm].rare_mm) << endl;
	fOUTSTREAM << Form("jes      lnN\t%5.3f/%5.3f\t-\t\t-\t\t%5.3f/%5.3f\t%5.3f/%5.3f",
					   1.0+(ttwzpreds[ijd].ttwz_mm-ttwzpreds[inm].ttwz_mm)/ttwzpreds[inm].ttwz_mm,
					   1.0+(ttwzpreds[iju].ttwz_mm-ttwzpreds[inm].ttwz_mm)/ttwzpreds[inm].ttwz_mm,
					   1.0+(ttwzpreds[ijd].wz_mm  -ttwzpreds[inm].wz_mm  )/ttwzpreds[inm].wz_mm,
					   1.0+(ttwzpreds[iju].wz_mm  -ttwzpreds[inm].wz_mm  )/ttwzpreds[inm].wz_mm,
					   1.0+(ttwzpreds[ijd].rare_mm-ttwzpreds[inm].rare_mm)/ttwzpreds[inm].rare_mm,
					   1.0+(ttwzpreds[iju].rare_mm-ttwzpreds[inm].rare_mm)/ttwzpreds[inm].rare_mm) << endl;
	fOUTSTREAM << Form("jer      lnN\t%5.3f\t\t-\t\t-\t\t%5.3f\t\t%5.3f",
					   1.0+(ttwzpreds[ijs].ttwz_mm-ttwzpreds[inm].ttwz_mm)/ttwzpreds[inm].ttwz_mm,
					   1.0+(ttwzpreds[ijs].wz_mm  -ttwzpreds[inm].wz_mm  )/ttwzpreds[inm].wz_mm,
					   1.0+(ttwzpreds[ijs].rare_mm-ttwzpreds[inm].rare_mm)/ttwzpreds[inm].rare_mm) << endl;
	fOUTSTREAM << Form("pu       lnN\t%5.3f\t\t-\t\t-\t\t%5.3f\t\t%5.3f",
					   pu_syst, pu_syst, pu_syst) << endl;
	fOUTSTREAM << Form("matching lnN\t%5.3f/%5.3f\t-\t\t-\t\t%5.3f/%5.3f\t%5.3f/%5.3f",
					   match_syst_up, match_syst_dn, match_syst_up, match_syst_dn, match_syst_up, match_syst_dn) << endl;
	fOUTSTREAM << Form("scale    lnN\t%5.3f/%5.3f\t-\t\t-\t\t%5.3f/%5.3f\t%5.3f/%5.3f",
					   scale_syst_up, scale_syst_dn, scale_syst_up, scale_syst_dn, scale_syst_up, scale_syst_dn) << endl;
	fOUTSTREAM << endl;
	fOUTSTREAM.close();
	
	datacard = outputdir + "datacard_ssdl_em.txt";
	fOUTSTREAM.open(datacard.Data(), ios::trunc);
	fOUTSTREAM <<      "#=========================================================================================" << endl;
	fOUTSTREAM <<      "# Systematics table for ttW/Z analysis, same-sign channel, emu only" << endl;
	fOUTSTREAM << Form("# Generated on: %s ", asctime(timeinfo)) << endl;
	fOUTSTREAM <<      "# Copy between the dashed lines for datacard" << endl;
	fOUTSTREAM <<      "#-----------------------------------------------------------------------------------------" << endl;
	fOUTSTREAM <<      "imax 1" << endl;
	fOUTSTREAM <<      "jmax 4" << endl;
	fOUTSTREAM <<      "kmax *" << endl;
	fOUTSTREAM << endl << endl;
	fOUTSTREAM <<      "bin\t\t1" << endl;
	fOUTSTREAM << Form("observation\t%d", ttwzpreds[inm].obs_em) << endl;
	fOUTSTREAM << endl << endl;
	fOUTSTREAM <<      "bin\t\t1\t\t1\t\t1\t\t1\t\t1" << endl;
	fOUTSTREAM <<      "process\t\tttWZ\t\tfake\t\tcmid\t\twz\t\trare" << endl;
	fOUTSTREAM <<      "process\t\t0\t\t1\t\t2\t\t3\t\t4" << endl;
	fOUTSTREAM << Form("rate\t\t%5.3f\t\t%5.3f\t\t%5.3f\t\t%5.3f\t\t%5.3f",
					   ttwzpreds[inm].ttwz_em, ttwzpreds[inm].fake_em, ttwzpreds[inm].cmid_em, ttwzpreds[inm].wz_em, ttwzpreds[inm].rare_em) << endl;
	fOUTSTREAM << Form("# separate signal yields: %5.3f (ttW), %5.3f (ttZ)", ttwzpreds[inm].ttw_em, ttwzpreds[inm].ttz_em) << endl;
	fOUTSTREAM << endl << endl;
	fOUTSTREAM <<      "#syst" << endl;
	fOUTSTREAM <<      "lumi     lnN\t1.022\t\t1.022\t\t1.022\t\t1.022\t\t1.022" << endl;
	fOUTSTREAM << Form("bgUncfak lnN\t-\t\t%5.3f\t\t-\t\t-\t\t-", 1.0+ttwzpreds[inm].fake_err_em/ttwzpreds[inm].fake_em) << endl;
	fOUTSTREAM << Form("bgUnccmi lnN\t-\t\t-\t\t%5.3f\t\t-\t\t-", 1.0+ttwzpreds[inm].cmid_err_em/ttwzpreds[inm].cmid_em) << endl;
	fOUTSTREAM << Form("bgUncwz  lnN\t-\t\t-\t\t-\t\t%5.3f\t\t-", 1.0+ttwzpreds[inm].wz_err_em  /ttwzpreds[inm].wz_em)   << endl;
	fOUTSTREAM << Form("bgUncrar lnN\t-\t\t-\t\t-\t\t-\t\t%5.3f", 1.0+ttwzpreds[inm].rare_err_em/ttwzpreds[inm].rare_em) << endl;
	fOUTSTREAM << Form("lept     lnN\t%5.3f/%5.3f\t-\t\t-\t\t%5.3f/%5.3f\t%5.3f/%5.3f",
					   1.0+(ttwzpreds[ild].ttwz_em-ttwzpreds[inm].ttwz_em)/ttwzpreds[inm].ttwz_em,
					   1.0+(ttwzpreds[ilu].ttwz_em-ttwzpreds[inm].ttwz_em)/ttwzpreds[inm].ttwz_em,
					   1.0+(ttwzpreds[ild].wz_em  -ttwzpreds[inm].wz_em  )/ttwzpreds[inm].wz_em,
					   1.0+(ttwzpreds[ilu].wz_em  -ttwzpreds[inm].wz_em  )/ttwzpreds[inm].wz_em,
					   1.0+(ttwzpreds[ild].rare_em-ttwzpreds[inm].rare_em)/ttwzpreds[inm].rare_em,
					   1.0+(ttwzpreds[ilu].rare_em-ttwzpreds[inm].rare_em)/ttwzpreds[inm].rare_em) << endl;
	fOUTSTREAM << Form("btag     lnN\t%5.3f/%5.3f\t-\t\t-\t\t%5.3f/%5.3f\t%5.3f/%5.3f",
					   1.0+(ttwzpreds[ibd].ttwz_em-ttwzpreds[inm].ttwz_em)/ttwzpreds[inm].ttwz_em,
					   1.0+(ttwzpreds[ibu].ttwz_em-ttwzpreds[inm].ttwz_em)/ttwzpreds[inm].ttwz_em,
					   1.0+(ttwzpreds[ibd].wz_em  -ttwzpreds[inm].wz_em  )/ttwzpreds[inm].wz_em,
					   1.0+(ttwzpreds[ibu].wz_em  -ttwzpreds[inm].wz_em  )/ttwzpreds[inm].wz_em,
					   1.0+(ttwzpreds[ibd].rare_em-ttwzpreds[inm].rare_em)/ttwzpreds[inm].rare_em,
					   1.0+(ttwzpreds[ibu].rare_em-ttwzpreds[inm].rare_em)/ttwzpreds[inm].rare_em) << endl;
	fOUTSTREAM << Form("jes      lnN\t%5.3f/%5.3f\t-\t\t-\t\t%5.3f/%5.3f\t%5.3f/%5.3f",
					   1.0+(ttwzpreds[ijd].ttwz_em-ttwzpreds[inm].ttwz_em)/ttwzpreds[inm].ttwz_em,
					   1.0+(ttwzpreds[iju].ttwz_em-ttwzpreds[inm].ttwz_em)/ttwzpreds[inm].ttwz_em,
					   1.0+(ttwzpreds[ijd].wz_em  -ttwzpreds[inm].wz_em  )/ttwzpreds[inm].wz_em,
					   1.0+(ttwzpreds[iju].wz_em  -ttwzpreds[inm].wz_em  )/ttwzpreds[inm].wz_em,
					   1.0+(ttwzpreds[ijd].rare_em-ttwzpreds[inm].rare_em)/ttwzpreds[inm].rare_em,
					   1.0+(ttwzpreds[iju].rare_em-ttwzpreds[inm].rare_em)/ttwzpreds[inm].rare_em) << endl;
	fOUTSTREAM << Form("jer      lnN\t%5.3f\t\t-\t\t-\t\t%5.3f\t\t%5.3f",
					   1.0+(ttwzpreds[ijs].ttwz_em-ttwzpreds[inm].ttwz_em)/ttwzpreds[inm].ttwz_em,
					   1.0+(ttwzpreds[ijs].wz_em  -ttwzpreds[inm].wz_em  )/ttwzpreds[inm].wz_em,
					   1.0+(ttwzpreds[ijs].rare_em-ttwzpreds[inm].rare_em)/ttwzpreds[inm].rare_em) << endl;
	fOUTSTREAM << Form("pu       lnN\t%5.3f\t\t-\t\t-\t\t%5.3f\t\t%5.3f",
					   pu_syst, pu_syst, pu_syst) << endl;
	fOUTSTREAM << Form("matching lnN\t%5.3f/%5.3f\t-\t\t-\t\t%5.3f/%5.3f\t%5.3f/%5.3f",
					   match_syst_up, match_syst_dn, match_syst_up, match_syst_dn, match_syst_up, match_syst_dn) << endl;
	fOUTSTREAM << Form("scale    lnN\t%5.3f/%5.3f\t-\t\t-\t\t%5.3f/%5.3f\t%5.3f/%5.3f",
					   scale_syst_up, scale_syst_dn, scale_syst_up, scale_syst_dn, scale_syst_up, scale_syst_dn) << endl;
	fOUTSTREAM << endl;
	fOUTSTREAM.close();
	
	datacard = outputdir + "datacard_ssdl_ee.txt";
	fOUTSTREAM.open(datacard.Data(), ios::trunc);
	fOUTSTREAM <<      "#=========================================================================================" << endl;
	fOUTSTREAM <<      "# Systematics table for ttW/Z analysis, same-sign channel, ee only" << endl;
	fOUTSTREAM << Form("# Generated on: %s ", asctime(timeinfo)) << endl;
	fOUTSTREAM <<      "# Copy between the dashed lines for datacard" << endl;
	fOUTSTREAM <<      "#-----------------------------------------------------------------------------------------" << endl;
	fOUTSTREAM <<      "imax 1" << endl;
	fOUTSTREAM <<      "jmax 4" << endl;
	fOUTSTREAM <<      "kmax *" << endl;
	fOUTSTREAM << endl << endl;
	fOUTSTREAM <<      "bin\t\t1" << endl;
	fOUTSTREAM << Form("observation\t%d", ttwzpreds[inm].obs_ee) << endl;
	fOUTSTREAM << endl << endl;
	fOUTSTREAM <<      "bin\t\t1\t\t1\t\t1\t\t1\t\t1" << endl;
	fOUTSTREAM <<      "process\t\tttWZ\t\tfake\t\tcmid\t\twz\t\trare" << endl;
	fOUTSTREAM <<      "process\t\t0\t\t1\t\t2\t\t3\t\t4" << endl;
	fOUTSTREAM << Form("rate\t\t%5.3f\t\t%5.3f\t\t%5.3f\t\t%5.3f\t\t%5.3f",
					   ttwzpreds[inm].ttwz_ee, ttwzpreds[inm].fake_ee, ttwzpreds[inm].cmid_ee, ttwzpreds[inm].wz_ee, ttwzpreds[inm].rare_ee) << endl;
	fOUTSTREAM << Form("# separate signal yields: %5.3f (ttW), %5.3f (ttZ)", ttwzpreds[inm].ttw_ee, ttwzpreds[inm].ttz_ee) << endl;
	fOUTSTREAM << endl << endl;
	fOUTSTREAM <<      "#syst" << endl;
	fOUTSTREAM <<      "lumi     lnN\t1.022\t\t1.022\t\t1.022\t\t1.022\t\t1.022" << endl;
	fOUTSTREAM << Form("bgUncfak lnN\t-\t\t%5.3f\t\t-\t\t-\t\t-", 1.0+ttwzpreds[inm].fake_err_ee/ttwzpreds[inm].fake_ee) << endl;
	fOUTSTREAM << Form("bgUnccmi lnN\t-\t\t-\t\t%5.3f\t\t-\t\t-", 1.0+ttwzpreds[inm].cmid_err_ee/ttwzpreds[inm].cmid_ee) << endl;
	fOUTSTREAM << Form("bgUncwz  lnN\t-\t\t-\t\t-\t\t%5.3f\t\t-", 1.0+ttwzpreds[inm].wz_err_ee  /ttwzpreds[inm].wz_ee)   << endl;
	fOUTSTREAM << Form("bgUncrar lnN\t-\t\t-\t\t-\t\t-\t\t%5.3f", 1.0+ttwzpreds[inm].rare_err_ee/ttwzpreds[inm].rare_ee) << endl;
	fOUTSTREAM << Form("lept     lnN\t%5.3f/%5.3f\t-\t\t-\t\t%5.3f/%5.3f\t%5.3f/%5.3f",
					   1.0+(ttwzpreds[ild].ttwz_ee-ttwzpreds[inm].ttwz_ee)/ttwzpreds[inm].ttwz_ee,
					   1.0+(ttwzpreds[ilu].ttwz_ee-ttwzpreds[inm].ttwz_ee)/ttwzpreds[inm].ttwz_ee,
					   1.0+(ttwzpreds[ild].wz_ee  -ttwzpreds[inm].wz_ee  )/ttwzpreds[inm].wz_ee,
					   1.0+(ttwzpreds[ilu].wz_ee  -ttwzpreds[inm].wz_ee  )/ttwzpreds[inm].wz_ee,
					   1.0+(ttwzpreds[ild].rare_ee-ttwzpreds[inm].rare_ee)/ttwzpreds[inm].rare_ee,
					   1.0+(ttwzpreds[ilu].rare_ee-ttwzpreds[inm].rare_ee)/ttwzpreds[inm].rare_ee) << endl;
	fOUTSTREAM << Form("btag     lnN\t%5.3f/%5.3f\t-\t\t-\t\t%5.3f/%5.3f\t%5.3f/%5.3f",
					   1.0+(ttwzpreds[ibd].ttwz_ee-ttwzpreds[inm].ttwz_ee)/ttwzpreds[inm].ttwz_ee,
					   1.0+(ttwzpreds[ibu].ttwz_ee-ttwzpreds[inm].ttwz_ee)/ttwzpreds[inm].ttwz_ee,
					   1.0+(ttwzpreds[ibd].wz_ee  -ttwzpreds[inm].wz_ee  )/ttwzpreds[inm].wz_ee,
					   1.0+(ttwzpreds[ibu].wz_ee  -ttwzpreds[inm].wz_ee  )/ttwzpreds[inm].wz_ee,
					   1.0+(ttwzpreds[ibd].rare_ee-ttwzpreds[inm].rare_ee)/ttwzpreds[inm].rare_ee,
					   1.0+(ttwzpreds[ibu].rare_ee-ttwzpreds[inm].rare_ee)/ttwzpreds[inm].rare_ee) << endl;
	fOUTSTREAM << Form("jes      lnN\t%5.3f/%5.3f\t-\t\t-\t\t%5.3f/%5.3f\t%5.3f/%5.3f",
					   1.0+(ttwzpreds[ijd].ttwz_ee-ttwzpreds[inm].ttwz_ee)/ttwzpreds[inm].ttwz_ee,
					   1.0+(ttwzpreds[iju].ttwz_ee-ttwzpreds[inm].ttwz_ee)/ttwzpreds[inm].ttwz_ee,
					   1.0+(ttwzpreds[ijd].wz_ee  -ttwzpreds[inm].wz_ee  )/ttwzpreds[inm].wz_ee,
					   1.0+(ttwzpreds[iju].wz_ee  -ttwzpreds[inm].wz_ee  )/ttwzpreds[inm].wz_ee,
					   1.0+(ttwzpreds[ijd].rare_ee-ttwzpreds[inm].rare_ee)/ttwzpreds[inm].rare_ee,
					   1.0+(ttwzpreds[iju].rare_ee-ttwzpreds[inm].rare_ee)/ttwzpreds[inm].rare_ee) << endl;
	fOUTSTREAM << Form("jer      lnN\t%5.3f\t\t-\t\t-\t\t%5.3f\t\t%5.3f",
					   1.0+(ttwzpreds[ijs].ttwz_ee-ttwzpreds[inm].ttwz_ee)/ttwzpreds[inm].ttwz_ee,
					   1.0+(ttwzpreds[ijs].wz_ee  -ttwzpreds[inm].wz_ee  )/ttwzpreds[inm].wz_ee,
					   1.0+(ttwzpreds[ijs].rare_ee-ttwzpreds[inm].rare_ee)/ttwzpreds[inm].rare_ee) << endl;
	fOUTSTREAM << Form("pu       lnN\t%5.3f\t\t-\t\t-\t\t%5.3f\t\t%5.3f",
					   pu_syst, pu_syst, pu_syst) << endl;
	fOUTSTREAM << Form("matching lnN\t%5.3f/%5.3f\t-\t\t-\t\t%5.3f/%5.3f\t%5.3f/%5.3f",
					   match_syst_up, match_syst_dn, match_syst_up, match_syst_dn, match_syst_up, match_syst_dn) << endl;
	fOUTSTREAM << Form("scale    lnN\t%5.3f/%5.3f\t-\t\t-\t\t%5.3f/%5.3f\t%5.3f/%5.3f",
					   scale_syst_up, scale_syst_dn, scale_syst_up, scale_syst_dn, scale_syst_up, scale_syst_dn) << endl;
	fOUTSTREAM << endl;
	fOUTSTREAM.close();
	
	datacard = outputdir + "datacard_ssdl_3channels.txt";
	fOUTSTREAM.open(datacard.Data(), ios::trunc);
	fOUTSTREAM <<      "#=========================================================================================" << endl;
	fOUTSTREAM <<      "# Systematics table for ttW/Z analysis, same-sign channel, subchannels" << endl;
	fOUTSTREAM << Form("# Generated on: %s ", asctime(timeinfo)) << endl;
	fOUTSTREAM <<      "# Copy between the dashed lines for datacard" << endl;
	fOUTSTREAM <<      "#-----------------------------------------------------------------------------------------" << endl;
	fOUTSTREAM <<      "imax 3" << endl;
	fOUTSTREAM <<      "jmax 4" << endl;
	fOUTSTREAM <<      "kmax *" << endl;
	fOUTSTREAM << endl << endl;
	fOUTSTREAM <<      "bin\t\t1\t2\t3" << endl;
	fOUTSTREAM << Form("observation\t%d\t%d\t%d", ttwzpreds[inm].obs_mm, ttwzpreds[inm].obs_em, ttwzpreds[inm].obs_ee) << endl;
	fOUTSTREAM << endl << endl;
	fOUTSTREAM <<      "bin\t\t1\t\t1\t\t1\t\t1\t\t1\t\t2\t\t2\t\t2\t\t2\t\t2\t\t3\t\t3\t\t3\t\t3\t\t3" << endl;
	fOUTSTREAM <<      "process\t\tttWZ\t\tfake\t\tcmid\t\twz\t\trare\t\tttWZ\t\tfake\t\tcmid\t\twz\t\trare\t\tttWZ\t\tfake\t\tcmid\t\twz\t\trare" << endl;
	fOUTSTREAM <<      "process\t\t0\t\t1\t\t2\t\t3\t\t4\t\t0\t\t1\t\t2\t\t3\t\t4\t\t0\t\t1\t\t2\t\t3\t\t4" << endl;
	fOUTSTREAM << Form("rate\t\t%5.3f\t\t%5.3f\t\t%5.3f\t\t%5.3f\t\t%5.3f\t\t%5.3f\t\t%5.3f\t\t%5.3f\t\t%5.3f\t\t%5.3f\t\t%5.3f\t\t%5.3f\t\t%5.3f\t\t%5.3f\t\t%5.3f",
					   ttwzpreds[inm].ttwz_mm, ttwzpreds[inm].fake_mm, 0.0                   , ttwzpreds[inm].wz_mm, ttwzpreds[inm].rare_mm,
					   ttwzpreds[inm].ttwz_em, ttwzpreds[inm].fake_em, ttwzpreds[inm].cmid_em, ttwzpreds[inm].wz_em, ttwzpreds[inm].rare_em,
					   ttwzpreds[inm].ttwz_ee, ttwzpreds[inm].fake_ee, ttwzpreds[inm].cmid_ee, ttwzpreds[inm].wz_ee, ttwzpreds[inm].rare_ee) << endl;
	fOUTSTREAM << Form("# separate signal yields: [mm] %5.3f (ttW), %5.3f (ttZ); [em] %5.3f (ttW), %5.3f (ttZ); [ee] %5.3f (ttW), %5.3f (ttZ);",
					   ttwzpreds[inm].ttw_mm, ttwzpreds[inm].ttz_mm,
					   ttwzpreds[inm].ttw_em, ttwzpreds[inm].ttz_em,
					   ttwzpreds[inm].ttw_ee, ttwzpreds[inm].ttz_ee) << endl;
	fOUTSTREAM << endl << endl;
	fOUTSTREAM <<      "#syst" << endl;
	fOUTSTREAM <<      "lumi     lnN\t1.022\t\t1.022\t\t1.022\t\t1.022\t\t1.022\t\t1.022\t\t1.022\t\t1.022\t\t1.022\t\t1.022\t\t1.022\t\t1.022\t\t1.022\t\t1.022\t\t1.022" << endl;
	fOUTSTREAM << Form("bgUncfak lnN\t-\t\t%5.3f\t\t-\t\t-\t\t-\t\t-\t\t%5.3f\t\t-\t\t-\t\t-\t\t-\t\t%5.3f\t\t-\t\t-\t\t-",
					   1.0+ttwzpreds[inm].fake_err_mm/ttwzpreds[inm].fake_mm,
					   1.0+ttwzpreds[inm].fake_err_em/ttwzpreds[inm].fake_em,
					   1.0+ttwzpreds[inm].fake_err_ee/ttwzpreds[inm].fake_ee) << endl;
	fOUTSTREAM << Form("bgUnccmi lnN\t-\t\t-\t\t-\t\t-\t\t-\t\t-\t\t-\t\t%5.3f\t\t-\t\t-\t\t-\t\t-\t\t%5.3f\t\t-\t\t-",
					   1.0+ttwzpreds[inm].cmid_err_em/ttwzpreds[inm].cmid_em,
					   1.0+ttwzpreds[inm].cmid_err_ee/ttwzpreds[inm].cmid_ee) << endl;
	fOUTSTREAM << Form("bgUncwz  lnN\t-\t\t-\t\t-\t\t%5.3f\t\t-\t\t-\t\t-\t\t-\t\t%5.3f\t\t-\t\t-\t\t-\t\t-\t\t%5.3f\t\t-",
					   1.0+ttwzpreds[inm].wz_err_mm  /ttwzpreds[inm].wz_mm,
					   1.0+ttwzpreds[inm].wz_err_em  /ttwzpreds[inm].wz_em,
					   1.0+ttwzpreds[inm].wz_err_ee  /ttwzpreds[inm].wz_ee)   << endl;
	fOUTSTREAM << Form("bgUncrar lnN\t-\t\t-\t\t-\t\t-\t\t%5.3f\t\t-\t\t-\t\t-\t\t-\t\t%5.3f\t\t-\t\t-\t\t-\t\t-\t\t%5.3f",
					   1.0+ttwzpreds[inm].rare_err_mm/ttwzpreds[inm].rare_mm,
					   1.0+ttwzpreds[inm].rare_err_em/ttwzpreds[inm].rare_em,
					   1.0+ttwzpreds[inm].rare_err_ee/ttwzpreds[inm].rare_ee) << endl;
	fOUTSTREAM << Form("lept     lnN\t%5.3f/%5.3f\t-\t\t-\t\t%5.3f/%5.3f\t%5.3f/%5.3f\t%5.3f/%5.3f\t-\t\t-\t\t%5.3f/%5.3f\t%5.3f/%5.3f\t%5.3f/%5.3f\t-\t\t-\t\t%5.3f/%5.3f\t%5.3f/%5.3f",
					   1.0+(ttwzpreds[ild].ttwz_mm-ttwzpreds[inm].ttwz_mm)/ttwzpreds[inm].ttwz_mm,
					   1.0+(ttwzpreds[ilu].ttwz_mm-ttwzpreds[inm].ttwz_mm)/ttwzpreds[inm].ttwz_mm,
					   1.0+(ttwzpreds[ild].wz_mm  -ttwzpreds[inm].wz_mm  )/ttwzpreds[inm].wz_mm,
					   1.0+(ttwzpreds[ilu].wz_mm  -ttwzpreds[inm].wz_mm  )/ttwzpreds[inm].wz_mm,
					   1.0+(ttwzpreds[ild].rare_mm-ttwzpreds[inm].rare_mm)/ttwzpreds[inm].rare_mm,
					   1.0+(ttwzpreds[ilu].rare_mm-ttwzpreds[inm].rare_mm)/ttwzpreds[inm].rare_mm,
					   1.0+(ttwzpreds[ild].ttwz_em-ttwzpreds[inm].ttwz_em)/ttwzpreds[inm].ttwz_em,
					   1.0+(ttwzpreds[ilu].ttwz_em-ttwzpreds[inm].ttwz_em)/ttwzpreds[inm].ttwz_em,
					   1.0+(ttwzpreds[ild].wz_em  -ttwzpreds[inm].wz_em  )/ttwzpreds[inm].wz_em,
					   1.0+(ttwzpreds[ilu].wz_em  -ttwzpreds[inm].wz_em  )/ttwzpreds[inm].wz_em,
					   1.0+(ttwzpreds[ild].rare_em-ttwzpreds[inm].rare_em)/ttwzpreds[inm].rare_em,
					   1.0+(ttwzpreds[ilu].rare_em-ttwzpreds[inm].rare_em)/ttwzpreds[inm].rare_em,
					   1.0+(ttwzpreds[ild].ttwz_ee-ttwzpreds[inm].ttwz_ee)/ttwzpreds[inm].ttwz_ee,
					   1.0+(ttwzpreds[ilu].ttwz_ee-ttwzpreds[inm].ttwz_ee)/ttwzpreds[inm].ttwz_ee,
					   1.0+(ttwzpreds[ild].wz_ee  -ttwzpreds[inm].wz_ee  )/ttwzpreds[inm].wz_ee,
					   1.0+(ttwzpreds[ilu].wz_ee  -ttwzpreds[inm].wz_ee  )/ttwzpreds[inm].wz_ee,
					   1.0+(ttwzpreds[ild].rare_ee-ttwzpreds[inm].rare_ee)/ttwzpreds[inm].rare_ee,
					   1.0+(ttwzpreds[ilu].rare_ee-ttwzpreds[inm].rare_ee)/ttwzpreds[inm].rare_ee) << endl;
	fOUTSTREAM << Form("btag     lnN\t%5.3f/%5.3f\t-\t\t-\t\t%5.3f/%5.3f\t%5.3f/%5.3f\t%5.3f/%5.3f\t-\t\t-\t\t%5.3f/%5.3f\t%5.3f/%5.3f\t%5.3f/%5.3f\t-\t\t-\t\t%5.3f/%5.3f\t%5.3f/%5.3f",
					   1.0+(ttwzpreds[ibd].ttwz_mm-ttwzpreds[inm].ttwz_mm)/ttwzpreds[inm].ttwz_mm,
					   1.0+(ttwzpreds[ibu].ttwz_mm-ttwzpreds[inm].ttwz_mm)/ttwzpreds[inm].ttwz_mm,
					   1.0+(ttwzpreds[ibd].wz_mm  -ttwzpreds[inm].wz_mm  )/ttwzpreds[inm].wz_mm,
					   1.0+(ttwzpreds[ibu].wz_mm  -ttwzpreds[inm].wz_mm  )/ttwzpreds[inm].wz_mm,
					   1.0+(ttwzpreds[ibd].rare_mm-ttwzpreds[inm].rare_mm)/ttwzpreds[inm].rare_mm,
					   1.0+(ttwzpreds[ibu].rare_mm-ttwzpreds[inm].rare_mm)/ttwzpreds[inm].rare_mm,
					   1.0+(ttwzpreds[ibd].ttwz_em-ttwzpreds[inm].ttwz_em)/ttwzpreds[inm].ttwz_em,
					   1.0+(ttwzpreds[ibu].ttwz_em-ttwzpreds[inm].ttwz_em)/ttwzpreds[inm].ttwz_em,
					   1.0+(ttwzpreds[ibd].wz_em  -ttwzpreds[inm].wz_em  )/ttwzpreds[inm].wz_em,
					   1.0+(ttwzpreds[ibu].wz_em  -ttwzpreds[inm].wz_em  )/ttwzpreds[inm].wz_em,
					   1.0+(ttwzpreds[ibd].rare_em-ttwzpreds[inm].rare_em)/ttwzpreds[inm].rare_em,
					   1.0+(ttwzpreds[ibu].rare_em-ttwzpreds[inm].rare_em)/ttwzpreds[inm].rare_em,
					   1.0+(ttwzpreds[ibd].ttwz_ee-ttwzpreds[inm].ttwz_ee)/ttwzpreds[inm].ttwz_ee,
					   1.0+(ttwzpreds[ibu].ttwz_ee-ttwzpreds[inm].ttwz_ee)/ttwzpreds[inm].ttwz_ee,
					   1.0+(ttwzpreds[ibd].wz_ee  -ttwzpreds[inm].wz_ee  )/ttwzpreds[inm].wz_ee,
					   1.0+(ttwzpreds[ibu].wz_ee  -ttwzpreds[inm].wz_ee  )/ttwzpreds[inm].wz_ee,
					   1.0+(ttwzpreds[ibd].rare_ee-ttwzpreds[inm].rare_ee)/ttwzpreds[inm].rare_ee,
					   1.0+(ttwzpreds[ibu].rare_ee-ttwzpreds[inm].rare_ee)/ttwzpreds[inm].rare_ee) << endl;
	fOUTSTREAM << Form("jes      lnN\t%5.3f/%5.3f\t-\t\t-\t\t%5.3f/%5.3f\t%5.3f/%5.3f\t%5.3f/%5.3f\t-\t\t-\t\t%5.3f/%5.3f\t%5.3f/%5.3f\t%5.3f/%5.3f\t-\t\t-\t\t%5.3f/%5.3f\t%5.3f/%5.3f",
					   1.0+(ttwzpreds[ijd].ttwz_mm-ttwzpreds[inm].ttwz_mm)/ttwzpreds[inm].ttwz_mm,
					   1.0+(ttwzpreds[iju].ttwz_mm-ttwzpreds[inm].ttwz_mm)/ttwzpreds[inm].ttwz_mm,
					   1.0+(ttwzpreds[ijd].wz_mm  -ttwzpreds[inm].wz_mm  )/ttwzpreds[inm].wz_mm,
					   1.0+(ttwzpreds[iju].wz_mm  -ttwzpreds[inm].wz_mm  )/ttwzpreds[inm].wz_mm,
					   1.0+(ttwzpreds[ijd].rare_mm-ttwzpreds[inm].rare_mm)/ttwzpreds[inm].rare_mm,
					   1.0+(ttwzpreds[iju].rare_mm-ttwzpreds[inm].rare_mm)/ttwzpreds[inm].rare_mm,
					   1.0+(ttwzpreds[ijd].ttwz_em-ttwzpreds[inm].ttwz_em)/ttwzpreds[inm].ttwz_em,
					   1.0+(ttwzpreds[iju].ttwz_em-ttwzpreds[inm].ttwz_em)/ttwzpreds[inm].ttwz_em,
					   1.0+(ttwzpreds[ijd].wz_em  -ttwzpreds[inm].wz_em  )/ttwzpreds[inm].wz_em,
					   1.0+(ttwzpreds[iju].wz_em  -ttwzpreds[inm].wz_em  )/ttwzpreds[inm].wz_em,
					   1.0+(ttwzpreds[ijd].rare_em-ttwzpreds[inm].rare_em)/ttwzpreds[inm].rare_em,
					   1.0+(ttwzpreds[iju].rare_em-ttwzpreds[inm].rare_em)/ttwzpreds[inm].rare_em,
					   1.0+(ttwzpreds[ijd].ttwz_ee-ttwzpreds[inm].ttwz_ee)/ttwzpreds[inm].ttwz_ee,
					   1.0+(ttwzpreds[iju].ttwz_ee-ttwzpreds[inm].ttwz_ee)/ttwzpreds[inm].ttwz_ee,
					   1.0+(ttwzpreds[ijd].wz_ee  -ttwzpreds[inm].wz_ee  )/ttwzpreds[inm].wz_ee,
					   1.0+(ttwzpreds[iju].wz_ee  -ttwzpreds[inm].wz_ee  )/ttwzpreds[inm].wz_ee,
					   1.0+(ttwzpreds[ijd].rare_ee-ttwzpreds[inm].rare_ee)/ttwzpreds[inm].rare_ee,
					   1.0+(ttwzpreds[iju].rare_ee-ttwzpreds[inm].rare_ee)/ttwzpreds[inm].rare_ee) << endl;
	fOUTSTREAM << Form("jer      lnN\t%5.3f\t\t-\t\t-\t\t%5.3f\t\t%5.3f\t\t%5.3f\t\t-\t\t-\t\t%5.3f\t\t%5.3f\t\t%5.3f\t\t-\t\t-\t\t%5.3f\t\t%5.3f",
					   1.0+(ttwzpreds[ijs].ttwz_mm-ttwzpreds[inm].ttwz_mm)/ttwzpreds[inm].ttwz_mm,
					   1.0+(ttwzpreds[ijs].wz_mm  -ttwzpreds[inm].wz_mm  )/ttwzpreds[inm].wz_mm,
					   1.0+(ttwzpreds[ijs].rare_mm-ttwzpreds[inm].rare_mm)/ttwzpreds[inm].rare_mm,
					   1.0+(ttwzpreds[ijs].ttwz_em-ttwzpreds[inm].ttwz_em)/ttwzpreds[inm].ttwz_em,
					   1.0+(ttwzpreds[ijs].wz_em  -ttwzpreds[inm].wz_em  )/ttwzpreds[inm].wz_em,
					   1.0+(ttwzpreds[ijs].rare_em-ttwzpreds[inm].rare_em)/ttwzpreds[inm].rare_em,
					   1.0+(ttwzpreds[ijs].ttwz_ee-ttwzpreds[inm].ttwz_ee)/ttwzpreds[inm].ttwz_ee,
					   1.0+(ttwzpreds[ijs].wz_ee  -ttwzpreds[inm].wz_ee  )/ttwzpreds[inm].wz_ee,
					   1.0+(ttwzpreds[ijs].rare_ee-ttwzpreds[inm].rare_ee)/ttwzpreds[inm].rare_ee) << endl;
	fOUTSTREAM << Form("pu       lnN\t%5.3f\t\t-\t\t-\t\t%5.3f\t\t%5.3f\t\t%5.3f\t\t-\t\t-\t\t%5.3f\t\t%5.3f\t\t%5.3f\t\t-\t\t-\t\t%5.3f\t\t%5.3f",
					   pu_syst, pu_syst, pu_syst, pu_syst, pu_syst, pu_syst, pu_syst, pu_syst, pu_syst) << endl;
	fOUTSTREAM << Form("matching lnN\t%5.3f/%5.3f\t-\t\t-\t\t%5.3f/%5.3f\t%5.3f/%5.3f\t%5.3f/%5.3f\t-\t\t-\t\t%5.3f/%5.3f\t%5.3f/%5.3f\t%5.3f/%5.3f\t-\t\t-\t\t%5.3f/%5.3f\t%5.3f/%5.3f",
					   match_syst_up, match_syst_dn, match_syst_up, match_syst_dn, match_syst_up, match_syst_dn, match_syst_up, match_syst_dn, match_syst_up, match_syst_dn,
					   match_syst_up, match_syst_dn, match_syst_up, match_syst_dn, match_syst_up, match_syst_dn, match_syst_up, match_syst_dn) << endl;
	fOUTSTREAM << Form("scale    lnN\t%5.3f/%5.3f\t-\t\t-\t\t%5.3f/%5.3f\t%5.3f/%5.3f\t%5.3f/%5.3f\t-\t\t-\t\t%5.3f/%5.3f\t%5.3f/%5.3f\t%5.3f/%5.3f\t-\t\t-\t\t%5.3f/%5.3f\t%5.3f/%5.3f",
					   scale_syst_up, scale_syst_dn, scale_syst_up, scale_syst_dn, scale_syst_up, scale_syst_dn, scale_syst_up, scale_syst_dn, scale_syst_up, scale_syst_dn,
					   scale_syst_up, scale_syst_dn, scale_syst_up, scale_syst_dn, scale_syst_up, scale_syst_dn, scale_syst_up, scale_syst_dn) << endl;
	fOUTSTREAM << endl;
	fOUTSTREAM.close();
	
	
	datacard = outputdir + "datacard_ssdl_3channels_ttW.txt";
	fOUTSTREAM.open(datacard.Data(), ios::trunc);
	fOUTSTREAM <<      "#=========================================================================================" << endl;
	fOUTSTREAM <<      "# Systematics table for ttW/Z analysis, same-sign channel, subchannels" << endl;
	fOUTSTREAM << Form("# Generated on: %s ", asctime(timeinfo)) << endl;
	fOUTSTREAM <<      "# Copy between the dashed lines for datacard" << endl;
	fOUTSTREAM <<      "#-----------------------------------------------------------------------------------------" << endl;
	fOUTSTREAM <<      "imax 3" << endl;
	fOUTSTREAM <<      "jmax 5" << endl;
	fOUTSTREAM <<      "kmax *" << endl;
	fOUTSTREAM << endl << endl;
	fOUTSTREAM <<      "bin\t\t1\t2\t3" << endl;
	fOUTSTREAM << Form("observation\t%d\t%d\t%d", ttwzpreds[inm].obs_mm, ttwzpreds[inm].obs_em, ttwzpreds[inm].obs_ee) << endl;
	fOUTSTREAM << endl << endl;
	fOUTSTREAM <<      "bin\t\t1\t\t1\t\t1\t\t1\t\t1\t\t1\t\t2\t\t2\t\t2\t\t2\t\t2\t\t2\t\t3\t\t3\t\t3\t\t3\t\t3\t\t3" << endl;
	fOUTSTREAM <<      "process\t\tttW\t\tttZ\t\tfake\t\tcmid\t\twz\t\trare\t\tttW\t\tttZ\t\tfake\t\tcmid\t\twz\t\trare\t\tttW\t\tttZ\t\tfake\t\tcmid\t\twz\t\trare" << endl;
	fOUTSTREAM <<      "process\t\t0\t\t1\t\t2\t\t3\t\t4\t\t5\t\t0\t\t1\t\t2\t\t3\t\t4\t\t5\t\t0\t\t1\t\t2\t\t3\t\t4\t\t5" << endl;
	fOUTSTREAM << Form("rate\t\t%5.3f\t\t%5.3f\t\t%5.3f\t\t%5.3f\t\t%5.3f\t\t%5.3f\t\t%5.3f\t\t%5.3f\t\t%5.3f\t\t%5.3f\t\t%5.3f\t\t%5.3f\t\t%5.3f\t\t%5.3f\t\t%5.3f\t\t%5.3f\t\t%5.3f\t\t%5.3f",
					   ttwzpreds[inm].ttw_mm, ttwzpreds[inm].ttz_mm, ttwzpreds[inm].fake_mm, 0.0                   , ttwzpreds[inm].wz_mm, ttwzpreds[inm].rare_mm,
					   ttwzpreds[inm].ttw_em, ttwzpreds[inm].ttz_em, ttwzpreds[inm].fake_em, ttwzpreds[inm].cmid_em, ttwzpreds[inm].wz_em, ttwzpreds[inm].rare_em,
					   ttwzpreds[inm].ttw_ee, ttwzpreds[inm].ttz_ee, ttwzpreds[inm].fake_ee, ttwzpreds[inm].cmid_ee, ttwzpreds[inm].wz_ee, ttwzpreds[inm].rare_ee) << endl;
	fOUTSTREAM << endl << endl;
	fOUTSTREAM <<      "#syst" << endl;
	fOUTSTREAM <<      "lumi     lnN\t1.022\t\t1.022\t\t1.022\t\t1.022\t\t1.022\t\t1.022\t\t1.022\t\t1.022\t\t1.022\t\t1.022\t\t1.022\t\t1.022\t\t1.022\t\t1.022\t\t1.022\t\t1.022\t\t1.022\t\t1.022" << endl;
	fOUTSTREAM << Form("bgUncfak lnN\t-\t\t-\t\t%5.3f\t\t-\t\t-\t\t-\t\t-\t\t-\t\t%5.3f\t\t-\t\t-\t\t-\t\t-\t\t-\t\t%5.3f\t\t-\t\t-\t\t-",
					   1.0+ttwzpreds[inm].fake_err_mm/ttwzpreds[inm].fake_mm,
					   1.0+ttwzpreds[inm].fake_err_em/ttwzpreds[inm].fake_em,
					   1.0+ttwzpreds[inm].fake_err_ee/ttwzpreds[inm].fake_ee) << endl;
	fOUTSTREAM << Form("bgUnccmi lnN\t-\t\t-\t\t-\t\t-\t\t-\t\t-\t\t-\t\t-\t\t-\t\t%5.3f\t\t-\t\t-\t\t-\t\t-\t\t-\t\t%5.3f\t\t-\t\t-",
					   1.0+ttwzpreds[inm].cmid_err_em/ttwzpreds[inm].cmid_em,
					   1.0+ttwzpreds[inm].cmid_err_ee/ttwzpreds[inm].cmid_ee) << endl;
	fOUTSTREAM << Form("bgUncwz  lnN\t-\t\t-\t\t-\t\t-\t\t%5.3f\t\t-\t\t-\t\t-\t\t-\t\t-\t\t%5.3f\t\t-\t\t-\t\t-\t\t-\t\t-\t\t%5.3f\t\t-",
					   1.0+ttwzpreds[inm].wz_err_mm  /ttwzpreds[inm].wz_mm,
					   1.0+ttwzpreds[inm].wz_err_em  /ttwzpreds[inm].wz_em,
					   1.0+ttwzpreds[inm].wz_err_ee  /ttwzpreds[inm].wz_ee)   << endl;
	fOUTSTREAM << Form("bgUncrar lnN\t-\t\t-\t\t-\t\t-\t\t-\t\t%5.3f\t\t-\t\t-\t\t-\t\t-\t\t-\t\t%5.3f\t\t-\t\t-\t\t-\t\t-\t\t-\t\t%5.3f",
					   1.0+ttwzpreds[inm].rare_err_mm/ttwzpreds[inm].rare_mm,
					   1.0+ttwzpreds[inm].rare_err_em/ttwzpreds[inm].rare_em,
					   1.0+ttwzpreds[inm].rare_err_ee/ttwzpreds[inm].rare_ee) << endl;
	fOUTSTREAM << Form("lept     lnN\t%5.3f/%5.3f\t%5.3f/%5.3f\t-\t\t-\t\t%5.3f/%5.3f\t%5.3f/%5.3f\t%5.3f/%5.3f\t%5.3f/%5.3f\t-\t\t-\t\t%5.3f/%5.3f\t%5.3f/%5.3f\t%5.3f/%5.3f\t%5.3f/%5.3f\t-\t\t-\t\t%5.3f/%5.3f\t%5.3f/%5.3f",
					   1.0+(ttwzpreds[ild].ttwz_mm-ttwzpreds[inm].ttwz_mm)/ttwzpreds[inm].ttwz_mm,
					   1.0+(ttwzpreds[ilu].ttwz_mm-ttwzpreds[inm].ttwz_mm)/ttwzpreds[inm].ttwz_mm,
					   1.0+(ttwzpreds[ild].ttwz_mm-ttwzpreds[inm].ttwz_mm)/ttwzpreds[inm].ttwz_mm,
					   1.0+(ttwzpreds[ilu].ttwz_mm-ttwzpreds[inm].ttwz_mm)/ttwzpreds[inm].ttwz_mm,
					   1.0+(ttwzpreds[ild].wz_mm  -ttwzpreds[inm].wz_mm  )/ttwzpreds[inm].wz_mm,
					   1.0+(ttwzpreds[ilu].wz_mm  -ttwzpreds[inm].wz_mm  )/ttwzpreds[inm].wz_mm,
					   1.0+(ttwzpreds[ild].rare_mm-ttwzpreds[inm].rare_mm)/ttwzpreds[inm].rare_mm,
					   1.0+(ttwzpreds[ilu].rare_mm-ttwzpreds[inm].rare_mm)/ttwzpreds[inm].rare_mm,
					   1.0+(ttwzpreds[ild].ttwz_em-ttwzpreds[inm].ttwz_em)/ttwzpreds[inm].ttwz_em,
					   1.0+(ttwzpreds[ilu].ttwz_em-ttwzpreds[inm].ttwz_em)/ttwzpreds[inm].ttwz_em,
					   1.0+(ttwzpreds[ild].ttwz_em-ttwzpreds[inm].ttwz_em)/ttwzpreds[inm].ttwz_em,
					   1.0+(ttwzpreds[ilu].ttwz_em-ttwzpreds[inm].ttwz_em)/ttwzpreds[inm].ttwz_em,
					   1.0+(ttwzpreds[ild].wz_em  -ttwzpreds[inm].wz_em  )/ttwzpreds[inm].wz_em,
					   1.0+(ttwzpreds[ilu].wz_em  -ttwzpreds[inm].wz_em  )/ttwzpreds[inm].wz_em,
					   1.0+(ttwzpreds[ild].rare_em-ttwzpreds[inm].rare_em)/ttwzpreds[inm].rare_em,
					   1.0+(ttwzpreds[ilu].rare_em-ttwzpreds[inm].rare_em)/ttwzpreds[inm].rare_em,
					   1.0+(ttwzpreds[ild].ttwz_ee-ttwzpreds[inm].ttwz_ee)/ttwzpreds[inm].ttwz_ee,
					   1.0+(ttwzpreds[ilu].ttwz_ee-ttwzpreds[inm].ttwz_ee)/ttwzpreds[inm].ttwz_ee,
					   1.0+(ttwzpreds[ild].ttwz_ee-ttwzpreds[inm].ttwz_ee)/ttwzpreds[inm].ttwz_ee,
					   1.0+(ttwzpreds[ilu].ttwz_ee-ttwzpreds[inm].ttwz_ee)/ttwzpreds[inm].ttwz_ee,
					   1.0+(ttwzpreds[ild].wz_ee  -ttwzpreds[inm].wz_ee  )/ttwzpreds[inm].wz_ee,
					   1.0+(ttwzpreds[ilu].wz_ee  -ttwzpreds[inm].wz_ee  )/ttwzpreds[inm].wz_ee,
					   1.0+(ttwzpreds[ild].rare_ee-ttwzpreds[inm].rare_ee)/ttwzpreds[inm].rare_ee,
					   1.0+(ttwzpreds[ilu].rare_ee-ttwzpreds[inm].rare_ee)/ttwzpreds[inm].rare_ee) << endl;
	fOUTSTREAM << Form("btag     lnN\t%5.3f/%5.3f\t%5.3f/%5.3f\t-\t\t-\t\t%5.3f/%5.3f\t%5.3f/%5.3f\t%5.3f/%5.3f\t%5.3f/%5.3f\t-\t\t-\t\t%5.3f/%5.3f\t%5.3f/%5.3f\t%5.3f/%5.3f\t%5.3f/%5.3f\t-\t\t-\t\t%5.3f/%5.3f\t%5.3f/%5.3f",
					   1.0+(ttwzpreds[ibd].ttwz_mm-ttwzpreds[inm].ttwz_mm)/ttwzpreds[inm].ttwz_mm,
					   1.0+(ttwzpreds[ibu].ttwz_mm-ttwzpreds[inm].ttwz_mm)/ttwzpreds[inm].ttwz_mm,
					   1.0+(ttwzpreds[ibd].ttwz_mm-ttwzpreds[inm].ttwz_mm)/ttwzpreds[inm].ttwz_mm,
					   1.0+(ttwzpreds[ibu].ttwz_mm-ttwzpreds[inm].ttwz_mm)/ttwzpreds[inm].ttwz_mm,
					   1.0+(ttwzpreds[ibd].wz_mm  -ttwzpreds[inm].wz_mm  )/ttwzpreds[inm].wz_mm,
					   1.0+(ttwzpreds[ibu].wz_mm  -ttwzpreds[inm].wz_mm  )/ttwzpreds[inm].wz_mm,
					   1.0+(ttwzpreds[ibd].rare_mm-ttwzpreds[inm].rare_mm)/ttwzpreds[inm].rare_mm,
					   1.0+(ttwzpreds[ibu].rare_mm-ttwzpreds[inm].rare_mm)/ttwzpreds[inm].rare_mm,
					   1.0+(ttwzpreds[ibd].ttwz_em-ttwzpreds[inm].ttwz_em)/ttwzpreds[inm].ttwz_em,
					   1.0+(ttwzpreds[ibu].ttwz_em-ttwzpreds[inm].ttwz_em)/ttwzpreds[inm].ttwz_em,
					   1.0+(ttwzpreds[ibd].ttwz_em-ttwzpreds[inm].ttwz_em)/ttwzpreds[inm].ttwz_em,
					   1.0+(ttwzpreds[ibu].ttwz_em-ttwzpreds[inm].ttwz_em)/ttwzpreds[inm].ttwz_em,
					   1.0+(ttwzpreds[ibd].wz_em  -ttwzpreds[inm].wz_em  )/ttwzpreds[inm].wz_em,
					   1.0+(ttwzpreds[ibu].wz_em  -ttwzpreds[inm].wz_em  )/ttwzpreds[inm].wz_em,
					   1.0+(ttwzpreds[ibd].rare_em-ttwzpreds[inm].rare_em)/ttwzpreds[inm].rare_em,
					   1.0+(ttwzpreds[ibu].rare_em-ttwzpreds[inm].rare_em)/ttwzpreds[inm].rare_em,
					   1.0+(ttwzpreds[ibd].ttwz_ee-ttwzpreds[inm].ttwz_ee)/ttwzpreds[inm].ttwz_ee,
					   1.0+(ttwzpreds[ibu].ttwz_ee-ttwzpreds[inm].ttwz_ee)/ttwzpreds[inm].ttwz_ee,
					   1.0+(ttwzpreds[ibd].ttwz_ee-ttwzpreds[inm].ttwz_ee)/ttwzpreds[inm].ttwz_ee,
					   1.0+(ttwzpreds[ibu].ttwz_ee-ttwzpreds[inm].ttwz_ee)/ttwzpreds[inm].ttwz_ee,
					   1.0+(ttwzpreds[ibd].wz_ee  -ttwzpreds[inm].wz_ee  )/ttwzpreds[inm].wz_ee,
					   1.0+(ttwzpreds[ibu].wz_ee  -ttwzpreds[inm].wz_ee  )/ttwzpreds[inm].wz_ee,
					   1.0+(ttwzpreds[ibd].rare_ee-ttwzpreds[inm].rare_ee)/ttwzpreds[inm].rare_ee,
					   1.0+(ttwzpreds[ibu].rare_ee-ttwzpreds[inm].rare_ee)/ttwzpreds[inm].rare_ee) << endl;
	fOUTSTREAM << Form("jes      lnN\t%5.3f/%5.3f\t%5.3f/%5.3f\t-\t\t-\t\t%5.3f/%5.3f\t%5.3f/%5.3f\t%5.3f/%5.3f\t%5.3f/%5.3f\t-\t\t-\t\t%5.3f/%5.3f\t%5.3f/%5.3f\t%5.3f/%5.3f\t%5.3f/%5.3f\t-\t\t-\t\t%5.3f/%5.3f\t%5.3f/%5.3f",
					   1.0+(ttwzpreds[ijd].ttwz_mm-ttwzpreds[inm].ttwz_mm)/ttwzpreds[inm].ttwz_mm,
					   1.0+(ttwzpreds[iju].ttwz_mm-ttwzpreds[inm].ttwz_mm)/ttwzpreds[inm].ttwz_mm,
					   1.0+(ttwzpreds[ijd].ttwz_mm-ttwzpreds[inm].ttwz_mm)/ttwzpreds[inm].ttwz_mm,
					   1.0+(ttwzpreds[iju].ttwz_mm-ttwzpreds[inm].ttwz_mm)/ttwzpreds[inm].ttwz_mm,
					   1.0+(ttwzpreds[ijd].wz_mm  -ttwzpreds[inm].wz_mm  )/ttwzpreds[inm].wz_mm,
					   1.0+(ttwzpreds[iju].wz_mm  -ttwzpreds[inm].wz_mm  )/ttwzpreds[inm].wz_mm,
					   1.0+(ttwzpreds[ijd].rare_mm-ttwzpreds[inm].rare_mm)/ttwzpreds[inm].rare_mm,
					   1.0+(ttwzpreds[iju].rare_mm-ttwzpreds[inm].rare_mm)/ttwzpreds[inm].rare_mm,
					   1.0+(ttwzpreds[ijd].ttwz_em-ttwzpreds[inm].ttwz_em)/ttwzpreds[inm].ttwz_em,
					   1.0+(ttwzpreds[iju].ttwz_em-ttwzpreds[inm].ttwz_em)/ttwzpreds[inm].ttwz_em,
					   1.0+(ttwzpreds[ijd].ttwz_em-ttwzpreds[inm].ttwz_em)/ttwzpreds[inm].ttwz_em,
					   1.0+(ttwzpreds[iju].ttwz_em-ttwzpreds[inm].ttwz_em)/ttwzpreds[inm].ttwz_em,
					   1.0+(ttwzpreds[ijd].wz_em  -ttwzpreds[inm].wz_em  )/ttwzpreds[inm].wz_em,
					   1.0+(ttwzpreds[iju].wz_em  -ttwzpreds[inm].wz_em  )/ttwzpreds[inm].wz_em,
					   1.0+(ttwzpreds[ijd].rare_em-ttwzpreds[inm].rare_em)/ttwzpreds[inm].rare_em,
					   1.0+(ttwzpreds[iju].rare_em-ttwzpreds[inm].rare_em)/ttwzpreds[inm].rare_em,
					   1.0+(ttwzpreds[ijd].ttwz_ee-ttwzpreds[inm].ttwz_ee)/ttwzpreds[inm].ttwz_ee,
					   1.0+(ttwzpreds[iju].ttwz_ee-ttwzpreds[inm].ttwz_ee)/ttwzpreds[inm].ttwz_ee,
					   1.0+(ttwzpreds[ijd].ttwz_ee-ttwzpreds[inm].ttwz_ee)/ttwzpreds[inm].ttwz_ee,
					   1.0+(ttwzpreds[iju].ttwz_ee-ttwzpreds[inm].ttwz_ee)/ttwzpreds[inm].ttwz_ee,
					   1.0+(ttwzpreds[ijd].wz_ee  -ttwzpreds[inm].wz_ee  )/ttwzpreds[inm].wz_ee,
					   1.0+(ttwzpreds[iju].wz_ee  -ttwzpreds[inm].wz_ee  )/ttwzpreds[inm].wz_ee,
					   1.0+(ttwzpreds[ijd].rare_ee-ttwzpreds[inm].rare_ee)/ttwzpreds[inm].rare_ee,
					   1.0+(ttwzpreds[iju].rare_ee-ttwzpreds[inm].rare_ee)/ttwzpreds[inm].rare_ee) << endl;
	fOUTSTREAM << Form("jer      lnN\t%5.3f\t\t%5.3f\t\t-\t\t-\t\t%5.3f\t\t%5.3f\t\t%5.3f\t\t%5.3f\t\t-\t\t-\t\t%5.3f\t\t%5.3f\t\t%5.3f\t\t%5.3f\t\t-\t\t-\t\t%5.3f\t\t%5.3f",
					   1.0+(ttwzpreds[ijs].ttwz_mm-ttwzpreds[inm].ttwz_mm)/ttwzpreds[inm].ttwz_mm,
					   1.0+(ttwzpreds[ijs].ttwz_mm-ttwzpreds[inm].ttwz_mm)/ttwzpreds[inm].ttwz_mm,
					   1.0+(ttwzpreds[ijs].wz_mm  -ttwzpreds[inm].wz_mm  )/ttwzpreds[inm].wz_mm,
					   1.0+(ttwzpreds[ijs].rare_mm-ttwzpreds[inm].rare_mm)/ttwzpreds[inm].rare_mm,
					   1.0+(ttwzpreds[ijs].ttwz_em-ttwzpreds[inm].ttwz_em)/ttwzpreds[inm].ttwz_em,
					   1.0+(ttwzpreds[ijs].ttwz_em-ttwzpreds[inm].ttwz_em)/ttwzpreds[inm].ttwz_em,
					   1.0+(ttwzpreds[ijs].wz_em  -ttwzpreds[inm].wz_em  )/ttwzpreds[inm].wz_em,
					   1.0+(ttwzpreds[ijs].rare_em-ttwzpreds[inm].rare_em)/ttwzpreds[inm].rare_em,
					   1.0+(ttwzpreds[ijs].ttwz_ee-ttwzpreds[inm].ttwz_ee)/ttwzpreds[inm].ttwz_ee,
					   1.0+(ttwzpreds[ijs].ttwz_ee-ttwzpreds[inm].ttwz_ee)/ttwzpreds[inm].ttwz_ee,
					   1.0+(ttwzpreds[ijs].wz_ee  -ttwzpreds[inm].wz_ee  )/ttwzpreds[inm].wz_ee,
					   1.0+(ttwzpreds[ijs].rare_ee-ttwzpreds[inm].rare_ee)/ttwzpreds[inm].rare_ee) << endl;
	fOUTSTREAM << Form("pu       lnN\t%5.3f\t\t%5.3f\t\t-\t\t-\t\t%5.3f\t\t%5.3f\t\t%5.3f\t\t%5.3f\t\t-\t\t-\t\t%5.3f\t\t%5.3f\t\t%5.3f\t\t%5.3f\t\t-\t\t-\t\t%5.3f\t\t%5.3f",
					   pu_syst, pu_syst, pu_syst, pu_syst, pu_syst, pu_syst, pu_syst, pu_syst, pu_syst, pu_syst, pu_syst, pu_syst) << endl;
	fOUTSTREAM << Form("matching lnN\t%5.3f/%5.3f\t%5.3f/%5.3f\t-\t\t-\t\t%5.3f/%5.3f\t%5.3f/%5.3f\t%5.3f/%5.3f\t%5.3f/%5.3f\t-\t\t-\t\t%5.3f/%5.3f\t%5.3f/%5.3f\t%5.3f/%5.3f\t%5.3f/%5.3f\t-\t\t-\t\t%5.3f/%5.3f\t%5.3f/%5.3f",
					   match_syst_up, match_syst_dn, match_syst_up, match_syst_dn, match_syst_up, match_syst_dn, match_syst_up, match_syst_dn, match_syst_up, match_syst_dn, match_syst_up, match_syst_dn,
					   match_syst_up, match_syst_dn, match_syst_up, match_syst_dn, match_syst_up, match_syst_dn, match_syst_up, match_syst_dn, match_syst_up, match_syst_dn, match_syst_up, match_syst_dn) << endl;
	fOUTSTREAM << Form("scale    lnN\t%5.3f/%5.3f\t%5.3f/%5.3f\t-\t\t-\t\t%5.3f/%5.3f\t%5.3f/%5.3f\t%5.3f/%5.3f\t%5.3f/%5.3f\t-\t\t-\t\t%5.3f/%5.3f\t%5.3f/%5.3f\t%5.3f/%5.3f\t%5.3f/%5.3f\t-\t\t-\t\t%5.3f/%5.3f\t%5.3f/%5.3f",
					   scale_syst_up, scale_syst_dn, scale_syst_up, scale_syst_dn, scale_syst_up, scale_syst_dn, scale_syst_up, scale_syst_dn, scale_syst_up, scale_syst_dn, scale_syst_up, scale_syst_dn,
					   scale_syst_up, scale_syst_dn, scale_syst_up, scale_syst_dn, scale_syst_up, scale_syst_dn, scale_syst_up, scale_syst_dn, scale_syst_up, scale_syst_dn, scale_syst_up, scale_syst_dn) << endl;
	fOUTSTREAM << endl;
	fOUTSTREAM.close();
	
	// Produce systematics plots
	TH1D *h_ttwz_nom = new TH1D("h_ttwz",    "Nominal ttWZ",         4, 0., 4.);
	TH1D *h_ttwz_ju  = new TH1D("h_ttwz_ju", "ttWZ jets up",         4, 0., 4.);
	TH1D *h_ttwz_jd  = new TH1D("h_ttwz_jd", "ttWZ jets dn",         4, 0., 4.);
	TH1D *h_ttwz_js  = new TH1D("h_ttwz_js", "ttWZ jets sm",         4, 0., 4.);
	TH1D *h_ttwz_bu  = new TH1D("h_ttwz_bu", "ttWZ b-tags up",       4, 0., 4.);
	TH1D *h_ttwz_bd  = new TH1D("h_ttwz_bd", "ttWZ b-tags dn",       4, 0., 4.);
	TH1D *h_ttwz_lu  = new TH1D("h_ttwz_lu", "ttWZ lepton up",       4, 0., 4.);
	TH1D *h_ttwz_ld  = new TH1D("h_ttwz_ld", "ttWZ lepton dn",       4, 0., 4.);
	TH1D *h_bg_nom   = new TH1D("h_bg",      "Nominal background",   4, 0., 4.);
	TH1D *h_bg_ju    = new TH1D("h_bg_ju",   "background jets up",   4, 0., 4.);
	TH1D *h_bg_jd    = new TH1D("h_bg_jd",   "background jets dn",   4, 0., 4.);
	TH1D *h_bg_js    = new TH1D("h_bg_js",   "background jets sm",   4, 0., 4.);
	TH1D *h_bg_bu    = new TH1D("h_bg_bu",   "background b-tags up", 4, 0., 4.);
	TH1D *h_bg_bd    = new TH1D("h_bg_bd",   "background b-tags dn", 4, 0., 4.);
	TH1D *h_bg_lu    = new TH1D("h_bg_lu",   "background lepton up", 4, 0., 4.);
	TH1D *h_bg_ld    = new TH1D("h_bg_ld",   "background lepton dn", 4, 0., 4.);
	vector<TH1D*> histos;
	histos.push_back(h_ttwz_nom);
	histos.push_back(h_ttwz_ju);
	histos.push_back(h_ttwz_jd);
	histos.push_back(h_ttwz_js);
	histos.push_back(h_ttwz_bu);
	histos.push_back(h_ttwz_bd);
	histos.push_back(h_ttwz_lu);
	histos.push_back(h_ttwz_ld);
	histos.push_back(h_bg_nom);
	histos.push_back(h_bg_ju);
	histos.push_back(h_bg_jd);
	histos.push_back(h_bg_js);
	histos.push_back(h_bg_bu);
	histos.push_back(h_bg_bd);
	histos.push_back(h_bg_lu);
	histos.push_back(h_bg_ld);
	
	h_ttwz_nom->SetBinContent(1, ttwzpreds[inm].ttwz_ee);
	h_ttwz_nom->SetBinContent(2, ttwzpreds[inm].ttwz_mm);
	h_ttwz_nom->SetBinContent(3, ttwzpreds[inm].ttwz_em);
	h_ttwz_nom->SetBinContent(4, ttwzpreds[inm].ttwz);
	
	h_ttwz_ju ->SetBinContent(1, ttwzpreds[iju].ttwz_ee);
	h_ttwz_ju ->SetBinContent(2, ttwzpreds[iju].ttwz_mm);
	h_ttwz_ju ->SetBinContent(3, ttwzpreds[iju].ttwz_em);
	h_ttwz_ju ->SetBinContent(4, ttwzpreds[iju].ttwz);
	h_ttwz_jd ->SetBinContent(1, ttwzpreds[ijd].ttwz_ee);
	h_ttwz_jd ->SetBinContent(2, ttwzpreds[ijd].ttwz_mm);
	h_ttwz_jd ->SetBinContent(3, ttwzpreds[ijd].ttwz_em);
	h_ttwz_jd ->SetBinContent(4, ttwzpreds[ijd].ttwz);
	
	h_ttwz_js ->SetBinContent(1, ttwzpreds[ijs].ttwz_ee);
	h_ttwz_js ->SetBinContent(2, ttwzpreds[ijs].ttwz_mm);
	h_ttwz_js ->SetBinContent(3, ttwzpreds[ijs].ttwz_em);
	h_ttwz_js ->SetBinContent(4, ttwzpreds[ijs].ttwz);
	
	h_ttwz_bu ->SetBinContent(1, ttwzpreds[ibu].ttwz_ee);
	h_ttwz_bu ->SetBinContent(2, ttwzpreds[ibu].ttwz_mm);
	h_ttwz_bu ->SetBinContent(3, ttwzpreds[ibu].ttwz_em);
	h_ttwz_bu ->SetBinContent(4, ttwzpreds[ibu].ttwz);
	h_ttwz_bd ->SetBinContent(1, ttwzpreds[ibd].ttwz_ee);
	h_ttwz_bd ->SetBinContent(2, ttwzpreds[ibd].ttwz_mm);
	h_ttwz_bd ->SetBinContent(3, ttwzpreds[ibd].ttwz_em);
	h_ttwz_bd ->SetBinContent(4, ttwzpreds[ibd].ttwz);
	
	h_ttwz_lu ->SetBinContent(1, ttwzpreds[ilu].ttwz_ee);
	h_ttwz_lu ->SetBinContent(2, ttwzpreds[ilu].ttwz_mm);
	h_ttwz_lu ->SetBinContent(3, ttwzpreds[ilu].ttwz_em);
	h_ttwz_lu ->SetBinContent(4, ttwzpreds[ilu].ttwz);
	h_ttwz_ld ->SetBinContent(1, ttwzpreds[ild].ttwz_ee);
	h_ttwz_ld ->SetBinContent(2, ttwzpreds[ild].ttwz_mm);
	h_ttwz_ld ->SetBinContent(3, ttwzpreds[ild].ttwz_em);
	h_ttwz_ld ->SetBinContent(4, ttwzpreds[ild].ttwz);
	
	h_bg_nom->SetBinContent(1, ttwzpreds[inm].wz_ee+ttwzpreds[inm].rare_ee);
	h_bg_nom->SetBinContent(2, ttwzpreds[inm].wz_mm+ttwzpreds[inm].rare_mm);
	h_bg_nom->SetBinContent(3, ttwzpreds[inm].wz_em+ttwzpreds[inm].rare_em);
	h_bg_nom->SetBinContent(4, ttwzpreds[inm].wz   +ttwzpreds[inm].rare   );
	
	h_bg_ju ->SetBinContent(1, ttwzpreds[iju].wz_ee+ttwzpreds[iju].rare_ee);
	h_bg_ju ->SetBinContent(2, ttwzpreds[iju].wz_mm+ttwzpreds[iju].rare_mm);
	h_bg_ju ->SetBinContent(3, ttwzpreds[iju].wz_em+ttwzpreds[iju].rare_em);
	h_bg_ju ->SetBinContent(4, ttwzpreds[iju].wz   +ttwzpreds[iju].rare   );
	h_bg_jd ->SetBinContent(1, ttwzpreds[ijd].wz_ee+ttwzpreds[ijd].rare_ee);
	h_bg_jd ->SetBinContent(2, ttwzpreds[ijd].wz_mm+ttwzpreds[ijd].rare_mm);
	h_bg_jd ->SetBinContent(3, ttwzpreds[ijd].wz_em+ttwzpreds[ijd].rare_em);
	h_bg_jd ->SetBinContent(4, ttwzpreds[ijd].wz   +ttwzpreds[ijd].rare   );
	
	h_bg_js ->SetBinContent(1, ttwzpreds[ijs].wz_ee+ttwzpreds[ijs].rare_ee);
	h_bg_js ->SetBinContent(2, ttwzpreds[ijs].wz_mm+ttwzpreds[ijs].rare_mm);
	h_bg_js ->SetBinContent(3, ttwzpreds[ijs].wz_em+ttwzpreds[ijs].rare_em);
	h_bg_js ->SetBinContent(4, ttwzpreds[ijs].wz   +ttwzpreds[ijs].rare   );
	
	h_bg_bu ->SetBinContent(1, ttwzpreds[ibu].wz_ee+ttwzpreds[ibu].rare_ee);
	h_bg_bu ->SetBinContent(2, ttwzpreds[ibu].wz_mm+ttwzpreds[ibu].rare_mm);
	h_bg_bu ->SetBinContent(3, ttwzpreds[ibu].wz_em+ttwzpreds[ibu].rare_em);
	h_bg_bu ->SetBinContent(4, ttwzpreds[ibu].wz   +ttwzpreds[ibu].rare   );
	h_bg_bd ->SetBinContent(1, ttwzpreds[ibd].wz_ee+ttwzpreds[ibd].rare_ee);
	h_bg_bd ->SetBinContent(2, ttwzpreds[ibd].wz_mm+ttwzpreds[ibd].rare_mm);
	h_bg_bd ->SetBinContent(3, ttwzpreds[ibd].wz_em+ttwzpreds[ibd].rare_em);
	h_bg_bd ->SetBinContent(4, ttwzpreds[ibd].wz   +ttwzpreds[ibd].rare   );
	
	h_bg_lu ->SetBinContent(1, ttwzpreds[ilu].wz_ee+ttwzpreds[ilu].rare_ee);
	h_bg_lu ->SetBinContent(2, ttwzpreds[ilu].wz_mm+ttwzpreds[ilu].rare_mm);
	h_bg_lu ->SetBinContent(3, ttwzpreds[ilu].wz_em+ttwzpreds[ilu].rare_em);
	h_bg_lu ->SetBinContent(4, ttwzpreds[ilu].wz   +ttwzpreds[ilu].rare   );
	h_bg_ld ->SetBinContent(1, ttwzpreds[ild].wz_ee+ttwzpreds[ild].rare_ee);
	h_bg_ld ->SetBinContent(2, ttwzpreds[ild].wz_mm+ttwzpreds[ild].rare_mm);
	h_bg_ld ->SetBinContent(3, ttwzpreds[ild].wz_em+ttwzpreds[ild].rare_em);
	h_bg_ld ->SetBinContent(4, ttwzpreds[ild].wz   +ttwzpreds[ild].rare   );
	
	h_ttwz_nom->SetLineColor(kBlack);
	h_bg_nom  ->SetLineColor(kBlack);
	
	h_ttwz_js->SetLineColor(kGreen+1);
	h_ttwz_ju->SetLineColor(kGreen+1);
	h_ttwz_bu->SetLineColor(kGreen+1);
	h_ttwz_lu->SetLineColor(kGreen+1);
	h_bg_js  ->SetLineColor(kGreen+1);
	h_bg_ju  ->SetLineColor(kGreen+1);
	h_bg_bu  ->SetLineColor(kGreen+1);
	h_bg_lu  ->SetLineColor(kGreen+1);
	
	h_ttwz_jd->SetLineColor(kRed-3);
	h_ttwz_bd->SetLineColor(kRed-3);
	h_ttwz_ld->SetLineColor(kRed-3);
	h_bg_jd  ->SetLineColor(kRed-3);
	h_bg_bd  ->SetLineColor(kRed-3);
	h_bg_ld  ->SetLineColor(kRed-3);
	
	int cnt = 0;
	for(vector<TH1D*>::iterator it = histos.begin(); it < histos.end(); ++it){
		(*it)->GetXaxis()->SetBinLabel(1, "ee");
		(*it)->GetXaxis()->SetBinLabel(2, "#mu#mu");
		(*it)->GetXaxis()->SetBinLabel(3, "e#mu");
		(*it)->GetXaxis()->SetBinLabel(4, "Total");
		(*it)->GetXaxis()->SetLabelSize(0.06);
		if(cnt<8)(*it)->GetYaxis()->SetTitle("Signal Events");
		if(cnt>7)(*it)->GetYaxis()->SetTitle("Background Events");
		(*it)->GetYaxis()->SetTitleOffset(1.15);
		(*it)->SetMinimum(0);
		if(cnt<8)(*it)->SetMaximum(10);
		if(cnt>7)(*it)->SetMaximum(4);
		(*it)->SetLineWidth(2);
		(*it)->SetFillStyle(0);
		cnt++;
	}
	
	makeSystPlot("Syst_Sig_JES",  "JES",          h_ttwz_nom, h_ttwz_ju, h_ttwz_jd);
	makeSystPlot("Syst_Sig_bTag", "b Tag",        h_ttwz_nom, h_ttwz_bu, h_ttwz_bd);
	makeSystPlot("Syst_Sig_Lept", "Lepton Scale", h_ttwz_nom, h_ttwz_lu, h_ttwz_ld);
	makeSystPlot("Syst_Sig_JER",  "JER",          h_ttwz_nom, h_ttwz_js);
	makeSystPlot("Syst_Bg_JES",   "JES",          h_bg_nom, h_bg_ju, h_bg_jd);
	makeSystPlot("Syst_Bg_bTag",  "b Tag",        h_bg_nom, h_bg_bu, h_bg_bd);
	makeSystPlot("Syst_Bg_Lept",  "Lepton Scale", h_bg_nom, h_bg_lu, h_bg_ld);
	makeSystPlot("Syst_Bg_JER",   "JER",          h_bg_nom, h_bg_js);
	
}
void SSDLPlotter::makeSystPlot(TString outputname, TString label, TH1D *nom, TH1D *plus, TH1D *minus){
	TCanvas *c_temp = new TCanvas("c_temp", "C", 800, 800);
	c_temp->SetRightMargin(0.05);
	c_temp->SetLeftMargin(0.1);
	TLegend *leg = new TLegend(0.15,0.70,0.50,0.88);
	leg->AddEntry(nom,   "Mean",              "l");
	leg->AddEntry(plus,  label+" + 1 #sigma", "l");
	if(minus!=NULL) leg->AddEntry(minus, label+" + 1 #sigma", "l");
	leg->SetFillStyle(0);
	leg->SetTextFont(42);
	leg->SetTextSize(0.03);
	leg->SetBorderSize(0);
	
	plus ->DrawCopy("hist");
	if(minus!=NULL) minus->DrawCopy("hist same");
	nom  ->DrawCopy("hist same");
	leg  ->Draw();
	
	drawTopLineSim(0.50, 1.0, 0.11);
	gPad->RedrawAxis();
	Util::PrintPDF(c_temp, outputname, fOutputDir + fOutputSubDir);
}
SSPrediction SSDLPlotter::makeIntPrediction(TString filename, int reg){
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

	const float WZESyst  = 0.2;
	const float WZESyst2 = WZESyst*WZESyst;

	musamples  = fMuData;
	elsamples  = fEGData;
	emusamples = fMuEGData;

	OUT << "/////////////////////////////////////////////////////////////////////////////" << endl;
	OUT << " Producing integrated predictions for region " << gRegions[reg]->sname << endl;
	OUT << "  scaling MC to " << fLumiNorm << " /pb" << endl << endl;

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
		
		npp_mm += S->numbers[reg][Muon].npp;
		npf_mm += S->numbers[reg][Muon].npf + S->numbers[reg][Muon].nfp;
		nff_mm += S->numbers[reg][Muon].nff;			
	}
	for(size_t i = 0; i < emusamples.size(); ++i){
		Sample *S = fSamples[emusamples[i]];
		nt2_em  += S->numbers[reg][ElMu].nt2;
		nt10_em += S->numbers[reg][ElMu].nt10;
		nt01_em += S->numbers[reg][ElMu].nt01;
		nt0_em  += S->numbers[reg][ElMu].nt0;

		npp_em += S->numbers[reg][ElMu].npp;
		npf_em += S->numbers[reg][ElMu].npf;
		nfp_em += S->numbers[reg][ElMu].nfp;
		nff_em += S->numbers[reg][ElMu].nff;

		nt2_em_BB_os += S->region[reg][HighPt].em.nt20_OS_BB_pt->GetEntries(); // ele in barrel
		nt2_em_EE_os += S->region[reg][HighPt].em.nt20_OS_EE_pt->GetEntries(); // ele in endcal
	}
	for(size_t i = 0; i < elsamples.size(); ++i){
		Sample *S = fSamples[elsamples[i]];
		nt2_ee  += S->numbers[reg][Elec].nt2;
		nt10_ee += S->numbers[reg][Elec].nt10;
		nt0_ee  += S->numbers[reg][Elec].nt0;

		npp_ee += S->numbers[reg][Elec].npp;
		npf_ee += S->numbers[reg][Elec].npf + S->numbers[reg][Elec].nfp;
		nff_ee += S->numbers[reg][Elec].nff;
		
		nt2_ee_BB_os += S->region[reg][HighPt].ee.nt20_OS_BB_pt->GetEntries(); // both in barrel
		nt2_ee_EE_os += S->region[reg][HighPt].ee.nt20_OS_EE_pt->GetEntries(); // both in endcal
		nt2_ee_EB_os += S->region[reg][HighPt].ee.nt20_OS_EB_pt->GetEntries(); // one barrel, one endcap
	}

	///////////////////////////////////////////////////////////////////////////////////
	// PRINTOUT ///////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////
	OUT << "---------------------------------------------------------------------------------------------------------" << endl;
	OUT << "         RATIOS  ||     Mu-fRatio      |     Mu-pRatio      ||     El-fRatio      |     El-pRatio      ||" << endl;
	OUT << "                 ||"<< setw(7)  << setprecision(3) << mufratio_data  << " +/- " << setw(7) << setprecision(3) << mufratio_data_e  << " |";
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


	// Background MC
	for(size_t i = 0; i < fMCBG.size(); ++i){
		Sample *S = fSamples[fMCBG[i]];
		float scale = fLumiNorm / S->getLumi();

		float temp_nt2_mm  = gMMTrigScale*scale*S->region[reg][HighPt].mm.nt20_pt->Integral(0, getNFPtBins(Muon)+1); nt2sum_mm  += temp_nt2_mm ;
		float temp_nt10_mm = gMMTrigScale*scale*S->region[reg][HighPt].mm.nt10_pt->Integral(0, getNFPtBins(Muon)+1); nt10sum_mm += temp_nt10_mm;
		float temp_nt0_mm  = gMMTrigScale*scale*S->region[reg][HighPt].mm.nt00_pt->Integral(0, getNFPtBins(Muon)+1); nt0sum_mm  += temp_nt0_mm ;
		float temp_nt2_em  = gEMTrigScale*scale*S->region[reg][HighPt].em.nt20_pt->Integral(0, getNFPtBins(ElMu)+1); nt2sum_em  += temp_nt2_em ;
		float temp_nt10_em = gEMTrigScale*scale*S->region[reg][HighPt].em.nt10_pt->Integral(0, getNFPtBins(ElMu)+1); nt10sum_em += temp_nt10_em;
		float temp_nt01_em = gEMTrigScale*scale*S->region[reg][HighPt].em.nt01_pt->Integral(0, getNFPtBins(ElMu)+1); nt01sum_em += temp_nt01_em;
		float temp_nt0_em  = gEMTrigScale*scale*S->region[reg][HighPt].em.nt00_pt->Integral(0, getNFPtBins(ElMu)+1); nt0sum_em  += temp_nt0_em ;
		float temp_nt2_ee  = gEETrigScale*scale*S->region[reg][HighPt].ee.nt20_pt->Integral(0, getNFPtBins(Elec)+1); nt2sum_ee  += temp_nt2_ee ;
		float temp_nt10_ee = gEETrigScale*scale*S->region[reg][HighPt].ee.nt10_pt->Integral(0, getNFPtBins(Elec)+1); nt10sum_ee += temp_nt10_ee;
		float temp_nt0_ee  = gEETrigScale*scale*S->region[reg][HighPt].ee.nt00_pt->Integral(0, getNFPtBins(Elec)+1); nt0sum_ee  += temp_nt0_ee ;

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

 	// Signal MC
//  	for(gSample i = sample_begin; i < gNSAMPLES; i=gSample(i+1)){
//  		Sample *S = fSamples[i];
//  		if(S->datamc != 2) continue;
//  		float scale = fLumiNorm / S->getLumi();
 
//  		float temp_nt2_mm  = gMMTrigScale*scale*S->region[reg][HighPt].mm.nt20_pt->Integral(0, getNFPtBins(Muon)+1);
//  		float temp_nt10_mm = gMMTrigScale*scale*S->region[reg][HighPt].mm.nt10_pt->Integral(0, getNFPtBins(Muon)+1);
//  		float temp_nt0_mm  = gMMTrigScale*scale*S->region[reg][HighPt].mm.nt01_pt->Integral(0, getNFPtBins(Muon)+1);
//  		float temp_nt2_em  = gEMTrigScale*scale*S->region[reg][HighPt].em.nt20_pt->Integral(0, getNFPtBins(ElMu)+1);
//  		float temp_nt10_em = gEMTrigScale*scale*S->region[reg][HighPt].em.nt10_pt->Integral(0, getNFPtBins(ElMu)+1);
//  		float temp_nt01_em = gEMTrigScale*scale*S->region[reg][HighPt].em.nt01_pt->Integral(0, getNFPtBins(ElMu)+1);
//  		float temp_nt0_em  = gEMTrigScale*scale*S->region[reg][HighPt].em.nt00_pt->Integral(0, getNFPtBins(ElMu)+1);
//  		float temp_nt2_ee  = gEETrigScale*scale*S->region[reg][HighPt].ee.nt20_pt->Integral(0, getNFPtBins(Elec)+1);
//  		float temp_nt10_ee = gEETrigScale*scale*S->region[reg][HighPt].ee.nt10_pt->Integral(0, getNFPtBins(Elec)+1);
//  		float temp_nt0_ee  = gEETrigScale*scale*S->region[reg][HighPt].ee.nt01_pt->Integral(0, getNFPtBins(Elec)+1);
 
//  		TString tempname = S->sname;
//  		OUT << Form("%16s & %6.2f & %6.2f & %6.2f & %6.2f & %6.2f & %6.2f & %6.2f & %6.2f & %6.2f & %6.2f \\\\\n", (tempname.ReplaceAll("_","\\_")).Data(),
//  		temp_nt2_mm , temp_nt10_mm, temp_nt0_mm,
//  		temp_nt2_em , temp_nt10_em, temp_nt01_em, temp_nt0_em,
//  		temp_nt2_ee , temp_nt10_ee, temp_nt0_ee);
//  	}	
//  	OUT << "\\hline" << endl;
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
 
 	FR->setMFRatio(mufratio_data, mufratio_data_e); // set error to pure statistical of ratio
 	FR->setEFRatio(elfratio_data, elfratio_data_e);
 	FR->setMPRatio(mupratio_data, mupratio_data_e);
 	FR->setEPRatio(elpratio_data, elpratio_data_e);
 
 	FR->setMMNtl(nt2_mm, nt10_mm, nt0_mm);
 	FR->setEENtl(nt2_ee, nt10_ee, nt0_ee);
 	FR->setEMNtl(nt2_em, nt10_em, nt01_em, nt0_em);
 
//  	float nF_mm = npf_mm + nff_mm;
//  	float nF_em = npf_em+nfp_em+nff_em;
//  	float nF_ee = npf_ee+nff_ee;
//  	float nSF   = npf_mm + npf_em + nfp_em + npf_ee;
//  	float nDF   = nff_mm + nff_em + nff_ee;
//  	float nF    = nF_mm + nF_em + nF_ee;
	
// 	npp_mm = FR->getMMNpp();   npf_mm = FR->getMMNpf();	nff_mm = FR->getMMNff();
// 	npp_em = FR->getEMNpp();   npf_em = FR->getEMNpf();	nff_em = FR->getEMNff();
// 	npp_ee = FR->getEENpp();   npf_ee = FR->getEENpf();	nff_ee = FR->getEENff();

 	float nF_mm = npf_mm + nff_mm;
 	float nF_em = npf_em+nfp_em+nff_em;
 	float nF_ee = npf_ee+nff_ee;
 	float nSF   = npf_mm + npf_em + nfp_em + npf_ee;
 	float nDF   = nff_mm + nff_em + nff_ee;
 	float nF    = nF_mm + nF_em + nF_ee;

	
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
	float fbb(0.),fee(0.),feb(0.);
	float fbbE(0.),feeE(0.),febE(0.);
	float fbb_mc(0.),fee_mc(0.),feb_mc(0.);
	float fbbE_mc(0.),feeE_mc(0.),febE_mc(0.);

	calculateChMisIdProb(fEGData, BB, fbb, fbbE);
	calculateChMisIdProb(fEGData, EB, feb, febE);
	calculateChMisIdProb(fEGData, EE, fee, feeE);

	calculateChMisIdProb(fMCBG, BB, fbb_mc, fbbE_mc);
	calculateChMisIdProb(fMCBG, EB, feb_mc, febE_mc);
	calculateChMisIdProb(fMCBG, EE, fee_mc, feeE_mc);
	
	cout << " this is the number of OS events in EE, BB, and EB from the old intpredictions function: " << nt2_ee_EE_os << " " << nt2_ee_BB_os << " " << nt2_ee_EB_os << endl;
	cout << "this is the charge mis ID probabilities: fee, fbb, feb: " << fee << " " << fbb << " " << feb << endl;

  	// Simple error propagation assuming error on number of events is sqrt(N)
 	nt2_ee_chmid    = 2*fbb*nt2_ee_BB_os                           + 2*fee*nt2_ee_EE_os                      + 2*feb*nt2_ee_EB_os;
 	nt2_ee_chmid_e1 = sqrt( 4*fbb*fbb*FR->getEStat2(nt2_ee_BB_os)  + 4*fee*fee*FR->getEStat2(nt2_ee_EE_os)   + 4*feb*feb*FR->getEStat2(nt2_ee_EB_os) ); // stat only
	nt2_ee_chmid_e2 = sqrt( nt2_ee_BB_os*nt2_ee_BB_os*4*fbbE*fbbE  + 4*feeE*feeE*nt2_ee_EE_os*nt2_ee_EE_os   + 4*febE*febE*nt2_ee_EB_os*nt2_ee_EB_os ); // syst only
	
        nt2_em_chmid    = fbb*nt2_em_BB_os + fee*nt2_em_EE_os;
        nt2_em_chmid_e1 = sqrt( fbb*fbb*FR->getEStat2(nt2_em_BB_os) + fee*fee*FR->getEStat2(nt2_em_EE_os) );
        nt2_em_chmid_e2 = sqrt( nt2_em_BB_os*nt2_em_BB_os*fbbE*fbbE + nt2_em_EE_os*nt2_em_EE_os*feeE*feeE );
 
 	///////////////////////////////////////////////////////////////////////////////////
 	// PRINTOUT ///////////////////////////////////////////////////////////////////////
 	///////////////////////////////////////////////////////////////////////////////////
 	OUT << "-------------------------------------------------------------------------------------------" << endl;
 	OUT << "       E-ChMisID  ||     Barrel-Barrel    |    Barrel - EndCap   |    Endcap - EndCap    ||" << endl;
 	OUT << "-------------------------------------------------------------------------------------------" << endl;
 	OUT << "  Data            ||";
 	OUT << setw(8) << setprecision(2) << fbb  << " +/- " << setw(8) << setprecision(2) << fbbE  << " |";
 	OUT << setw(8) << setprecision(2) << feb  << " +/- " << setw(8) << setprecision(2) << febE  << " |";
 	OUT << setw(8) << setprecision(2) << fee  << " +/- " << setw(8) << setprecision(2) << feeE  << "  ||" << endl;
 	OUT << "  MC              ||";
 	OUT << setw(8) << setprecision(2) << fbb_mc  << " +/- " << setw(8) << setprecision(2) << fbbE_mc  << " |";
 	OUT << setw(8) << setprecision(2) << feb_mc  << " +/- " << setw(8) << setprecision(2) << febE_mc  << " |";
 	OUT << setw(8) << setprecision(2) << fee_mc  << " +/- " << setw(8) << setprecision(2) << feeE_mc  << "  ||";
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
 		float scale = fLumiNorm / S->getLumi();
 		
 		mc_os_em_bb_sum += gEMTrigScale*scale*S->region[reg][HighPt].em.nt20_OS_BB_pt->GetEntries();
 		mc_os_em_ee_sum += gEMTrigScale*scale*S->region[reg][HighPt].em.nt20_OS_EE_pt->GetEntries();
 		mc_os_ee_bb_sum += gEETrigScale*scale*S->region[reg][HighPt].ee.nt20_OS_BB_pt->GetEntries();
 		mc_os_ee_eb_sum += gEETrigScale*scale*S->region[reg][HighPt].ee.nt20_OS_EB_pt->GetEntries();
 		mc_os_ee_ee_sum += gEETrigScale*scale*S->region[reg][HighPt].ee.nt20_OS_EE_pt->GetEntries();
 
 		OUT << setw(16) << S->sname << " || ";
 		OUT << setw(7)  << setprecision(2) << gEMTrigScale*scale*S->region[reg][HighPt].em.nt20_OS_BB_pt->GetEntries() << " | ";
 		OUT << setw(7)  << setprecision(2) << gEMTrigScale*scale*S->region[reg][HighPt].em.nt20_OS_EE_pt->GetEntries() << " || ";
 		OUT << setw(7)  << setprecision(2) << gEETrigScale*scale*S->region[reg][HighPt].ee.nt20_OS_BB_pt->GetEntries() << " | ";
 		OUT << setw(7)  << setprecision(2) << gEETrigScale*scale*S->region[reg][HighPt].ee.nt20_OS_EB_pt->GetEntries() << " | ";
 		OUT << setw(7)  << setprecision(2) << gEETrigScale*scale*S->region[reg][HighPt].ee.nt20_OS_EE_pt->GetEntries() << " || ";
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
 	OUT << setw(7) << Form("%6.4f", 0.5*(feb+fbb) * nt2_em_BB_os ) << " | ";
 	OUT << setw(7) << Form("%6.4f", 0.5*(feb+fee) * nt2_em_EE_os ) << " || ";
 	OUT << setw(7) << Form("%6.4f",         fbb   * nt2_ee_BB_os ) << " | ";
 	OUT << setw(7) << Form("%6.4f",         feb   * nt2_ee_EB_os ) << " | ";
 	OUT << setw(7) << Form("%6.4f",         fee   * nt2_ee_EE_os ) << " || ";
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

 	vector<int> mcbkg;
	mcbkg.push_back(ZZ);
// MARC 	mcbkg.push_back(GVJets);
 	mcbkg.push_back(DPSWW);
	mcbkg.push_back(TTbarH);
	mcbkg.push_back(TTbarW);
	mcbkg.push_back(TTbarZ);
	mcbkg.push_back(TTbarG);
	mcbkg.push_back(TbZ);
 	mcbkg.push_back(WpWp);
 	mcbkg.push_back(WmWm);
 	mcbkg.push_back(WWZ);
       	mcbkg.push_back(WZZ);
 	mcbkg.push_back(WWG);
 	mcbkg.push_back(WWW);
 	mcbkg.push_back(ZZZ);
 	for(size_t i = 0; i < mcbkg.size(); ++i){
 		Sample *S = fSamples[mcbkg[i]];
 		float scale = fLumiNorm/S->getLumi();
  		
 		float temp_nt2_mm = gMMTrigScale*scale*S->region[reg][HighPt].mm.nt20_pt->Integral(0, getNFPtBins(Muon)+1);
 		float temp_nt2_em = gEMTrigScale*scale*S->region[reg][HighPt].em.nt20_pt->Integral(0, getNFPtBins(ElMu)+1);
 		float temp_nt2_ee = gEETrigScale*scale*S->region[reg][HighPt].ee.nt20_pt->Integral(0, getNFPtBins(Elec)+1);
 
 		nt2_rare_mc_mm += temp_nt2_mm;
 		nt2_rare_mc_em += temp_nt2_em;
 		nt2_rare_mc_ee += temp_nt2_ee;
 
 		nt2_rare_mc_mm_e1 += gMMTrigScale*gMMTrigScale*scale*scale*S->numbers[reg][Muon].tt_avweight*S->numbers[reg][Muon].tt_avweight * S->getError2(S->region[reg][HighPt].mm.nt20_pt->GetEntries()); // for stat error take actual entries, not pileup weighted integral...
 		nt2_rare_mc_em_e1 += gEMTrigScale*gEMTrigScale*scale*scale*S->numbers[reg][ElMu].tt_avweight*S->numbers[reg][ElMu].tt_avweight * S->getError2(S->region[reg][HighPt].em.nt20_pt->GetEntries());
 		nt2_rare_mc_ee_e1 += gEETrigScale*gEETrigScale*scale*scale*S->numbers[reg][Elec].tt_avweight*S->numbers[reg][Elec].tt_avweight * S->getError2(S->region[reg][HighPt].ee.nt20_pt->GetEntries());
 
 		OUT << Form("%16s || %5.2f ± %5.2f         || %5.2f ± %5.2f         || %5.2f ± %5.2f         ||\n", S->sname.Data(),
 		temp_nt2_mm, gMMTrigScale*S->numbers[reg][Muon].tt_avweight * scale*S->getError(S->region[reg][HighPt].mm.nt20_pt->GetEntries()),
 		temp_nt2_em, gEMTrigScale*S->numbers[reg][ElMu].tt_avweight * scale*S->getError(S->region[reg][HighPt].em.nt20_pt->GetEntries()),
 		temp_nt2_ee, gEETrigScale*S->numbers[reg][Elec].tt_avweight * scale*S->getError(S->region[reg][HighPt].ee.nt20_pt->GetEntries()));
 	}
 	OUT << "----------------------------------------------------------------------------------------------" << endl;
 	OUT << Form("%16s || %5.2f ± %5.2f ± %5.2f || %5.2f ± %5.2f ± %5.2f || %5.2f ± %5.2f ± %5.2f ||\n", "Rare SM (Sum)",
 	nt2_rare_mc_mm, sqrt(nt2_rare_mc_mm_e1), RareESyst*nt2_rare_mc_mm,
 	nt2_rare_mc_em, sqrt(nt2_rare_mc_em_e1), RareESyst*nt2_rare_mc_em,
 	nt2_rare_mc_ee, sqrt(nt2_rare_mc_ee_e1), RareESyst*nt2_rare_mc_ee);
 	OUT << "----------------------------------------------------------------------------------------------" << endl;
 	
 	///////////////////////////////////////////
 	// WZ production
	float wzscale(0.);
	float wz_nt2_mm(0.),wz_nt2_em(0.),wz_nt2_ee(0.);
	float wz_nt2_mm_e1(0.),wz_nt2_em_e1(0.),wz_nt2_ee_e1(0.);
	
	wzscale = fLumiNorm/fSamples[WZ]->getLumi();
	wz_nt2_mm = gMMTrigScale*wzscale*fSamples[WZ]->region[reg][HighPt].mm.nt20_pt->Integral(0, getNFPtBins(Muon)+1);
	wz_nt2_em = gEMTrigScale*wzscale*fSamples[WZ]->region[reg][HighPt].em.nt20_pt->Integral(0, getNFPtBins(ElMu)+1);
	wz_nt2_ee = gEETrigScale*wzscale*fSamples[WZ]->region[reg][HighPt].ee.nt20_pt->Integral(0, getNFPtBins(Elec)+1);
 
	wz_nt2_mm_e1 = gMMTrigScale*gMMTrigScale*wzscale*wzscale * fSamples[WZ]->numbers[reg][Muon].tt_avweight*fSamples[WZ]->numbers[reg][Muon].tt_avweight * fSamples[WZ]->getError2(fSamples[WZ]->region[reg][HighPt].mm.nt20_pt->GetEntries()); // for stat error take actual entries, not pileup weighted integral...
 	wz_nt2_em_e1 = gEMTrigScale*gEMTrigScale*wzscale*wzscale * fSamples[WZ]->numbers[reg][ElMu].tt_avweight*fSamples[WZ]->numbers[reg][ElMu].tt_avweight * fSamples[WZ]->getError2(fSamples[WZ]->region[reg][HighPt].em.nt20_pt->GetEntries());
 	wz_nt2_ee_e1 = gEETrigScale*gEETrigScale*wzscale*wzscale * fSamples[WZ]->numbers[reg][Elec].tt_avweight*fSamples[WZ]->numbers[reg][Elec].tt_avweight * fSamples[WZ]->getError2(fSamples[WZ]->region[reg][HighPt].ee.nt20_pt->GetEntries());
 
 	OUT << Form("%16s || %5.2f ± %5.2f         || %5.2f ± %5.2f         || %5.2f ± %5.2f         ||\n", "WZ Prod.",
	wz_nt2_mm, sqrt(wz_nt2_mm_e1),
 	wz_nt2_em, sqrt(wz_nt2_em_e1),
 	wz_nt2_ee, sqrt(wz_nt2_ee_e1));
 	OUT << "----------------------------------------------------------------------------------------------" << endl;
 	// Just add different errors in quadrature (they are independent)
 	float mm_tot_sqerr1 = FR->getMMTotEStat()*FR->getMMTotEStat()                                   + nt2_rare_mc_mm_e1 + wz_nt2_mm_e1;
 	float em_tot_sqerr1 = FR->getEMTotEStat()*FR->getEMTotEStat() + nt2_em_chmid_e1*nt2_em_chmid_e1 + nt2_rare_mc_em_e1 + wz_nt2_em_e1;
 	float ee_tot_sqerr1 = FR->getEETotEStat()*FR->getEETotEStat() + nt2_ee_chmid_e1*nt2_ee_chmid_e1 + nt2_rare_mc_ee_e1 + wz_nt2_ee_e1;
 	float mm_tot_sqerr2 = nF_mm*nF_mm*FakeESyst2                                   + RareESyst2*nt2_rare_mc_mm*nt2_rare_mc_mm + WZESyst2*wz_nt2_mm*wz_nt2_mm;
 	float em_tot_sqerr2 = nF_em*nF_em*FakeESyst2 + nt2_em_chmid_e2*nt2_em_chmid_e2 + RareESyst2*nt2_rare_mc_em*nt2_rare_mc_em + WZESyst2*wz_nt2_em*wz_nt2_em;
 	float ee_tot_sqerr2 = nF_ee*nF_ee*FakeESyst2 + nt2_ee_chmid_e2*nt2_ee_chmid_e2 + RareESyst2*nt2_rare_mc_ee*nt2_rare_mc_ee + WZESyst2*wz_nt2_ee*wz_nt2_ee;
 	OUT << Form("%16s || %5.2f ± %5.2f ± %5.2f || %5.2f ± %5.2f ± %5.2f || %5.2f ± %5.2f ± %5.2f ||\n", "tot. backgr.",
 	nF_mm                + nt2_rare_mc_mm + wz_nt2_mm, sqrt(mm_tot_sqerr1), sqrt(mm_tot_sqerr2),
 	nF_em + nt2_em_chmid + nt2_rare_mc_em + wz_nt2_em, sqrt(em_tot_sqerr1), sqrt(em_tot_sqerr2),
 	nF_ee + nt2_ee_chmid + nt2_rare_mc_ee + wz_nt2_ee, sqrt(ee_tot_sqerr1), sqrt(ee_tot_sqerr2));
 	OUT << Form("%16s || %5.2f ± %5.2f         || %5.2f ± %5.2f         || %5.2f ± %5.2f         ||\n", "",
 	nF_mm                + nt2_rare_mc_mm + wz_nt2_mm, sqrt(mm_tot_sqerr1 + mm_tot_sqerr2),
 	nF_em + nt2_em_chmid + nt2_rare_mc_em + wz_nt2_em, sqrt(em_tot_sqerr1 + em_tot_sqerr2),
 	nF_ee + nt2_ee_chmid + nt2_rare_mc_ee + wz_nt2_ee, sqrt(ee_tot_sqerr1 + ee_tot_sqerr2));
 	OUT << "----------------------------------------------------------------------------------------------" << endl;
 	OUT << Form("%16s || %5.2f                 || %5.2f                 || %5.2f                 ||\n", "tot. MC", nt2sum_mm, nt2sum_em, nt2sum_ee);
 	OUT << "==============================================================================================" << endl;
 	OUT << Form("%16s || %2.0f                    || %2.0f                    || %2.0f                    ||\n", "observed", nt2_mm, nt2_em, nt2_ee);
 	OUT << "==============================================================================================" << endl;
 	OUT << "        predicted: ";
 	float tot_pred        = nF + nt2_rare_mc_mm + nt2_em_chmid + nt2_rare_mc_em + nt2_ee_chmid + nt2_rare_mc_ee + wz_nt2_mm + wz_nt2_em + wz_nt2_ee;
 	float comb_tot_sqerr1 = FR->getTotEStat()*FR->getTotEStat() + nt2_rare_mc_mm_e1 + nt2_em_chmid_e1*nt2_em_chmid_e1 + nt2_rare_mc_em_e1 + nt2_ee_chmid_e1*nt2_ee_chmid_e1 + nt2_rare_mc_ee_e1 + wz_nt2_mm_e1 + wz_nt2_em_e1 + wz_nt2_ee_e1;
 	float comb_tot_sqerr2 = nF*nF*FakeESyst2 + RareESyst2*(nt2_rare_mc_mm + nt2_rare_mc_em + nt2_rare_mc_ee)*(nt2_rare_mc_mm + nt2_rare_mc_em + nt2_rare_mc_ee) + nt2_em_chmid_e2*nt2_em_chmid_e2 + nt2_ee_chmid_e2*nt2_ee_chmid_e2 + WZESyst2*(wz_nt2_mm+wz_nt2_em+wz_nt2_ee)*(wz_nt2_mm+wz_nt2_em+wz_nt2_ee);
 	// FIXME: Why take 50% on Rare yields on SUM and not on individual channels? because they're 100% correlated
 	OUT << setw(5) << left << Form("%5.2f", tot_pred ) << " ± ";
 	OUT << setw(5) << Form("%5.2f", sqrt(comb_tot_sqerr1)) << " ± ";
 	OUT << setw(5) << Form("%5.2f", sqrt(comb_tot_sqerr2)) << endl;
 	OUT << "      combined MC: ";
 	OUT << setw(5) << left << Form("%5.2f", nt2sum_mm+nt2sum_em+nt2sum_ee ) << endl;
 	OUT << "combined observed: ";
 	OUT << setw(5) << left << Form("%2.0f", nt2_mm+nt2_em+nt2_ee ) << endl;
 	OUT << "==============================================================================================" << endl;
 	OUT.close();
 		
	///////////////////////////////////////////////////////////////////////////////////
 	//  OUTPUT FOR COMBINATION TOOL  //////////////////////////////////////////////////
 	///////////////////////////////////////////////////////////////////////////////////
 	fOUTSTREAM4 << "'"+gRegions[reg]->sname +"':{"<< endl;
 	fOUTSTREAM4 << Form("'f': %6.2f, 'fstat': %6.2f, 'fsyst':  %6.2f,",
 	nff_em + nff_mm + nff_ee + npf_em + nfp_em + npf_mm + npf_ee ,// sum of double and single fakes
	sqrt(FR->getTotDoubleEStat()*FR->getTotDoubleEStat() + FR->getTotSingleEStat()*FR->getTotSingleEStat()),// stat error on the fakes
	0.5*(nff_em + nff_mm + nff_ee + npf_em + nfp_em + npf_mm + npf_ee)) << endl;
 	fOUTSTREAM4 << Form("'c': %6.2f, 'cstat': %6.2f, 'csyst':  %6.2f,",
 	nt2_ee_chmid + nt2_em_chmid, // flips
	sqrt(nt2_ee_chmid_e1*nt2_ee_chmid_e1 + nt2_em_chmid_e1*nt2_em_chmid_e1), //stat error on flips
	sqrt(nt2_ee_chmid_e2*nt2_ee_chmid_e2 + nt2_em_chmid_e2*nt2_em_chmid_e2)) << endl;//syst error on flips
 	fOUTSTREAM4 << Form("'r': %6.2f, 'rstat': %6.2f, 'rsyst':  %6.2f,",
 	nt2_rare_mc_ee + nt2_rare_mc_mm + nt2_rare_mc_em + wz_nt2_mm + wz_nt2_em + wz_nt2_ee,                          // rares
	sqrt(nt2_rare_mc_ee_e1 + nt2_rare_mc_mm_e1 + nt2_rare_mc_em_e1 + wz_nt2_mm_e1 + wz_nt2_em_e1 + wz_nt2_em_e1), // stat error on rare
	sqrt(RareESyst2*(nt2_rare_mc_ee + nt2_rare_mc_mm + nt2_rare_mc_em)*(nt2_rare_mc_ee + nt2_rare_mc_mm + nt2_rare_mc_em) + WZESyst2*(wz_nt2_mm + wz_nt2_em + wz_nt2_ee)*(wz_nt2_mm + wz_nt2_em + wz_nt2_ee))) << endl; //syst error on rares
 	fOUTSTREAM4 << Form("'o': %4.0f", nt2_mm+nt2_em+nt2_ee) << endl; // observed
 	fOUTSTREAM4 << "}," << endl;

 	///////////////////////////////////////////////////////////////////////////////////
 	//  OUTPUT FOR ANALYSIS NOTE  /////////////////////////////////////////////////////
 	///////////////////////////////////////////////////////////////////////////////////
 	fOUTSTREAM3 << "%% " + gRegions[reg]->sname << endl;
 	fOUTSTREAM3 << "-----------------------------------------------------------------" << endl;
 	fOUTSTREAM3 << Form("DF:  %6.2f ± %6.2f  ( %5.2f±%5.2f | %5.2f±%5.2f | %5.2f±%5.2f )\n",
 	nff_em + nff_mm + nff_ee, sqrt(FR->getTotDoubleEStat()*FR->getTotDoubleEStat() + nDF*nDF*FakeESyst2),
 	nff_ee, sqrt(FR->getEENffEStat()*FR->getEENffEStat()+nff_ee*nff_ee*FakeESyst2),
 	nff_mm, sqrt(FR->getMMNffEStat()*FR->getMMNffEStat()+nff_mm*nff_mm*FakeESyst2),
 	nff_em, sqrt(FR->getEMNffEStat()*FR->getEMNffEStat()+nff_em*nff_em*FakeESyst2));
 	fOUTSTREAM3 << Form("SF:  %6.2f ± %6.2f  ( %5.2f±%5.2f | %5.2f±%5.2f | %5.2f±%5.2f )\n",
 	npf_em + nfp_em + npf_mm + npf_ee, sqrt(FR->getTotSingleEStat()*FR->getTotSingleEStat() + nSF*nSF*FakeESyst2),
 	npf_ee,          sqrt(FR->getEENpfEStat()   *FR->getEENpfEStat()    +  npf_ee*npf_ee*FakeESyst2),
 	npf_mm,          sqrt(FR->getMMNpfEStat()   *FR->getMMNpfEStat()    +  npf_mm*npf_mm*FakeESyst2),
 	npf_em + nfp_em, sqrt(FR->getEMSingleEStat()*FR->getEMSingleEStat() + (npf_em+nfp_em)*(npf_em+nfp_em)*FakeESyst2));
 	fOUTSTREAM3 << Form("CM:  %6.2f ± %6.2f  ( %5.2f±%5.2f |   -         | %5.2f±%5.2f )\n",
 	nt2_ee_chmid + nt2_em_chmid, sqrt(nt2_ee_chmid_e1*nt2_ee_chmid_e1 + nt2_ee_chmid_e2*nt2_ee_chmid_e2 + nt2_em_chmid_e1*nt2_em_chmid_e1 + nt2_em_chmid_e2*nt2_em_chmid_e2),
 	nt2_ee_chmid, sqrt(nt2_ee_chmid_e1*nt2_ee_chmid_e1 + nt2_ee_chmid_e2*nt2_ee_chmid_e2),
 	nt2_em_chmid, sqrt(nt2_em_chmid_e1*nt2_em_chmid_e1 + nt2_em_chmid_e2*nt2_em_chmid_e2));
 	fOUTSTREAM3 << Form("MC:  %6.2f ± %6.2f  ( %5.2f±%5.2f | %5.2f±%5.2f | %5.2f±%5.2f )\n",
 	nt2_rare_mc_ee + nt2_rare_mc_mm + nt2_rare_mc_em, sqrt(nt2_rare_mc_ee_e1 + nt2_rare_mc_mm_e1 + nt2_rare_mc_em_e1 + RareESyst2*(nt2_rare_mc_ee + nt2_rare_mc_mm + nt2_rare_mc_em)*(nt2_rare_mc_ee + nt2_rare_mc_mm + nt2_rare_mc_em)),
 	nt2_rare_mc_ee, sqrt(nt2_rare_mc_ee_e1 + RareESyst2*nt2_rare_mc_ee*nt2_rare_mc_ee),
 	nt2_rare_mc_mm, sqrt(nt2_rare_mc_mm_e1 + RareESyst2*nt2_rare_mc_mm*nt2_rare_mc_mm),
 	nt2_rare_mc_em, sqrt(nt2_rare_mc_em_e1 + RareESyst2*nt2_rare_mc_em*nt2_rare_mc_em));
 	fOUTSTREAM3 << Form("WZ:  %6.2f ± %6.2f  ( %5.2f±%5.2f | %5.2f±%5.2f | %5.2f±%5.2f )\n",
 	wz_nt2_mm + wz_nt2_em + wz_nt2_ee,
 	sqrt(wz_nt2_mm_e1 + wz_nt2_em_e1 + wz_nt2_em_e1 + WZESyst2*(wz_nt2_mm + wz_nt2_em + wz_nt2_ee)*(wz_nt2_mm + wz_nt2_em + wz_nt2_ee)),
 	wz_nt2_ee, sqrt(wz_nt2_ee_e1 + WZESyst2*wz_nt2_ee*wz_nt2_ee),
 	wz_nt2_mm, sqrt(wz_nt2_mm_e1 + WZESyst2*wz_nt2_mm*wz_nt2_mm),
 	wz_nt2_em, sqrt(wz_nt2_em_e1 + WZESyst2*wz_nt2_em*wz_nt2_em));
 	fOUTSTREAM3 << Form("Tot: %6.2f ± %6.2f  ( %5.2f±%5.2f | %5.2f±%5.2f | %5.2f±%5.2f )\n",
 	tot_pred, sqrt(comb_tot_sqerr1 + comb_tot_sqerr2),
 	nF_ee + nt2_rare_mc_ee + nt2_ee_chmid + wz_nt2_ee, sqrt(ee_tot_sqerr1 + ee_tot_sqerr2),
 	nF_mm + nt2_rare_mc_mm                + wz_nt2_mm, sqrt(mm_tot_sqerr1 + mm_tot_sqerr2),
 	nF_em + nt2_rare_mc_em + nt2_em_chmid + wz_nt2_em, sqrt(em_tot_sqerr1 + em_tot_sqerr2));
 	fOUTSTREAM3 << "-----------------------------------------------------------------" << endl;
 	fOUTSTREAM3 << Form("Obs: %4.0f             ( %3.0f         | %3.0f         | %3.0f         )\n", nt2_mm+nt2_em+nt2_ee, nt2_ee, nt2_mm, nt2_em);
 	fOUTSTREAM3 << "-----------------------------------------------------------------" << endl;
 	fOUTSTREAM3 << endl;
 	
	
 	///////////////////////////////////////////////////////////////////////////////////
 	//  OUTPUT FOR AN TABLE  //////////////////////////////////////////////////////////
 	///////////////////////////////////////////////////////////////////////////////////
 	fOUTSTREAM << gRegions[reg]->sname << endl;
 	fOUTSTREAM << Form("Double Fakes   & %5.2f $\\pm$ %5.2f & %5.2f $\\pm$ %5.2f & %5.2f $\\pm$ %5.2f & %5.2f $\\pm$ %5.2f \\\\ \n",
 	nff_mm, sqrt(FR->getMMNffEStat()*FR->getMMNffEStat()+nff_mm*nff_mm*FakeESyst2),
 	nff_em, sqrt(FR->getEMNffEStat()*FR->getEMNffEStat()+nff_em*nff_em*FakeESyst2),
 	nff_ee, sqrt(FR->getEENffEStat()*FR->getEENffEStat()+nff_ee*nff_ee*FakeESyst2),
 	nff_em + nff_mm + nff_ee, sqrt(FR->getTotDoubleEStat()*FR->getTotDoubleEStat() + nDF*nDF*FakeESyst2));
 	fOUTSTREAM << Form("Single Fakes   & %5.2f $\\pm$ %5.2f & %5.2f $\\pm$ %5.2f & %5.2f $\\pm$ %5.2f & %5.2f $\\pm$ %5.2f \\\\ \n",
 	npf_mm,          sqrt(FR->getMMNpfEStat()   *FR->getMMNpfEStat()    +  npf_mm*npf_mm*FakeESyst2),
 	npf_em + nfp_em, sqrt(FR->getEMSingleEStat()*FR->getEMSingleEStat() + (npf_em+nfp_em)*(npf_em+nfp_em)*FakeESyst2),
 	npf_ee,          sqrt(FR->getEENpfEStat()   *FR->getEENpfEStat()    +  npf_ee*npf_ee*FakeESyst2),
 	npf_em + nfp_em + npf_mm + npf_ee, sqrt(FR->getTotSingleEStat()*FR->getTotSingleEStat() + nSF*nSF*FakeESyst2));
 	fOUTSTREAM << Form("Charge MisID   &        -          & %5.2f $\\pm$ %5.2f & %5.2f $\\pm$ %5.2f & %5.2f $\\pm$ %5.2f \\\\ \n",
 	nt2_em_chmid, sqrt(nt2_em_chmid_e1*nt2_em_chmid_e1 + nt2_em_chmid_e2*nt2_em_chmid_e2),
 	nt2_ee_chmid, sqrt(nt2_ee_chmid_e1*nt2_ee_chmid_e1 + nt2_ee_chmid_e2*nt2_ee_chmid_e2),
 	nt2_ee_chmid + nt2_em_chmid, sqrt(nt2_ee_chmid_e1*nt2_ee_chmid_e1 + nt2_ee_chmid_e2*nt2_ee_chmid_e2 + nt2_em_chmid_e1*nt2_em_chmid_e1 + nt2_em_chmid_e2*nt2_em_chmid_e2));
 	fOUTSTREAM << Form("Rare SM        & %5.2f $\\pm$ %5.2f & %5.2f $\\pm$ %5.2f & %5.2f $\\pm$ %5.2f & %5.2f $\\pm$ %5.2f \\\\ \n",
 	nt2_rare_mc_mm, sqrt(nt2_rare_mc_mm_e1 + RareESyst2*nt2_rare_mc_mm*nt2_rare_mc_mm),
 	nt2_rare_mc_em, sqrt(nt2_rare_mc_em_e1 + RareESyst2*nt2_rare_mc_em*nt2_rare_mc_em),
 	nt2_rare_mc_ee, sqrt(nt2_rare_mc_ee_e1 + RareESyst2*nt2_rare_mc_ee*nt2_rare_mc_ee),
 	nt2_rare_mc_ee + nt2_rare_mc_mm + nt2_rare_mc_em, sqrt(nt2_rare_mc_ee_e1 + nt2_rare_mc_mm_e1 + nt2_rare_mc_em_e1 + RareESyst2*(nt2_rare_mc_ee + nt2_rare_mc_mm + nt2_rare_mc_em)*(nt2_rare_mc_ee + nt2_rare_mc_mm + nt2_rare_mc_em)));
 	fOUTSTREAM << Form("WZ Prod.       & %5.2f $\\pm$ %5.2f & %5.2f $\\pm$ %5.2f & %5.2f $\\pm$ %5.2f & %5.2f $\\pm$ %5.2f \\\\ \\hline \n",
 	wz_nt2_mm, sqrt(wz_nt2_mm_e1 + WZESyst2*wz_nt2_mm*wz_nt2_mm),
 	wz_nt2_em, sqrt(wz_nt2_em_e1 + WZESyst2*wz_nt2_em*wz_nt2_em),
 	wz_nt2_ee, sqrt(wz_nt2_ee_e1 + WZESyst2*wz_nt2_ee*wz_nt2_ee),
 	wz_nt2_ee + wz_nt2_mm + wz_nt2_em, sqrt(wz_nt2_mm_e1 + wz_nt2_ee_e1 + wz_nt2_em_e1 + WZESyst2*(wz_nt2_ee + wz_nt2_mm + wz_nt2_em)*(wz_nt2_ee + wz_nt2_mm + wz_nt2_em)));
 	fOUTSTREAM << Form("Total Bkg      & %5.2f $\\pm$ %5.2f & %5.2f $\\pm$ %5.2f & %5.2f $\\pm$ %5.2f & %5.2f $\\pm$ %5.2f \\\\ \\hline \n",
 	nF_mm + nt2_rare_mc_mm                + wz_nt2_mm, sqrt(mm_tot_sqerr1 + mm_tot_sqerr2),
 	nF_em + nt2_rare_mc_em + nt2_em_chmid + wz_nt2_em, sqrt(em_tot_sqerr1 + em_tot_sqerr2),
 	nF_ee + nt2_rare_mc_ee + nt2_ee_chmid + wz_nt2_ee, sqrt(ee_tot_sqerr1 + ee_tot_sqerr2),
 	tot_pred, sqrt(comb_tot_sqerr1 + comb_tot_sqerr2));
 	fOUTSTREAM << Form("Observed       & %3.0f               & %3.0f               & %3.0f               & %3.0f               \\\\ \n",
 	nt2_mm, nt2_em, nt2_ee, nt2_mm+nt2_em+nt2_ee);
 	fOUTSTREAM << endl;
 	
 	///////////////////////////////////////////////////////////////////////////////////
 	//  OUTPUT AS DATACARD  ///////////////////////////////////////////////////////////
 	///////////////////////////////////////////////////////////////////////////////////
	SSPrediction pred;
	pred.obs      = nt2_mm+nt2_em+nt2_ee;
 	pred.obs_mm   = nt2_mm;
 	pred.obs_ee   = nt2_ee;
 	pred.obs_em   = nt2_em;
	
	pred.s        = 1.;
	pred.s_mm     = 1.;
	pred.s_ee     = 1.;
	pred.s_em     = 1.;
	
	pred.bkg      = tot_pred;
	pred.bkg_err  = sqrt(comb_tot_sqerr2);
	pred.bkg_stat = sqrt(comb_tot_sqerr1);
	
	pred.fake         = nff_em + nff_mm + nff_ee + npf_em + nfp_em + npf_mm + npf_ee;
	pred.fake_mm      = nff_mm + npf_mm;
	pred.fake_em      = nff_em + npf_em + nfp_em;
	pred.fake_ee      = nff_ee + npf_ee;
	pred.fake_err     = sqrt(nDF*nDF*FakeESyst2 + nSF*nSF*FakeESyst2);
	pred.fake_err_mm  = FakeESyst*nF_mm;
	pred.fake_err_em  = FakeESyst*nF_em;
	pred.fake_err_ee  = FakeESyst*nF_ee;
	pred.fake_stat    = sqrt(FR->getMMTotEStat()*FR->getMMTotEStat()+FR->getEMTotEStat()*FR->getEMTotEStat()+FR->getEETotEStat()*FR->getEETotEStat());
	pred.fake_stat_mm = FR->getMMTotEStat();
	pred.fake_stat_em = FR->getEMTotEStat();
	pred.fake_stat_ee = FR->getEETotEStat();
	
	pred.sf           = npf_em + nfp_em + npf_mm + npf_ee;
	pred.sf_mm        = npf_mm;
	pred.sf_em        = npf_em + nfp_em;
	pred.sf_ee        = npf_ee;
	pred.sf_err       = FakeESyst*(npf_em + nfp_em + npf_mm + npf_ee);
	pred.sf_err_mm    = FakeESyst*npf_mm;
	pred.sf_err_em    = FakeESyst*(npf_em+nfp_em);
	pred.sf_err_ee    = FakeESyst*npf_ee;
	pred.sf_stat      = FR->getTotSingleEStat();
	pred.sf_stat_mm   = FR->getMMNpfEStat();
	pred.sf_stat_em   = FR->getEMSingleEStat();
	pred.sf_stat_ee   = FR->getEENpfEStat();

	pred.df           = nff_em + nff_mm + nff_ee;
	pred.df_mm        = nff_mm;
	pred.df_em        = nff_em;
	pred.df_ee        = nff_ee;
	pred.df_err       = FakeESyst*(nff_em + nff_mm + nff_ee);
	pred.df_err_mm    = FakeESyst*nff_mm;
	pred.df_err_em    = FakeESyst*nff_em;
	pred.df_err_ee    = FakeESyst*nff_ee;
	pred.df_stat      = FR->getTotDoubleEStat();
	pred.df_stat_mm   = FR->getMMNffEStat();
	pred.df_stat_em   = FR->getEMNffEStat();
	pred.df_stat_ee   = FR->getEENffEStat();
	
	pred.cmid         = nt2_ee_chmid + nt2_em_chmid;
	pred.cmid_mm      = 0;
	pred.cmid_em      = nt2_em_chmid;
	pred.cmid_ee      = nt2_ee_chmid;
	pred.cmid_err     = sqrt(nt2_ee_chmid_e2*nt2_ee_chmid_e2 + nt2_em_chmid_e2*nt2_em_chmid_e2);
	pred.cmid_err_mm  = 0;
	pred.cmid_err_em  = nt2_em_chmid_e2;
	pred.cmid_err_ee  = nt2_ee_chmid_e2;
	pred.cmid_stat    = sqrt(nt2_ee_chmid_e1*nt2_ee_chmid_e1);
	pred.cmid_stat_mm = 0;
	pred.cmid_stat_em = nt2_em_chmid_e1;
	pred.cmid_stat_ee = nt2_ee_chmid_e1;

	pred.wz           = wz_nt2_ee + wz_nt2_mm + wz_nt2_em;
	pred.wz_mm        = wz_nt2_mm;
	pred.wz_em        = wz_nt2_em;
	pred.wz_ee        = wz_nt2_ee;
	pred.wz_err       = sqrt(WZESyst2*(wz_nt2_ee + wz_nt2_mm + wz_nt2_em)*(wz_nt2_ee + wz_nt2_mm + wz_nt2_em));
	pred.wz_err_mm    = sqrt(WZESyst2*wz_nt2_mm*wz_nt2_mm);
	pred.wz_err_em    = sqrt(WZESyst2*wz_nt2_em*wz_nt2_em);
	pred.wz_err_ee    = sqrt(WZESyst2*wz_nt2_ee*wz_nt2_ee);
	pred.wz_stat      = sqrt(wz_nt2_mm_e1 + wz_nt2_ee_e1 + wz_nt2_em_e1);
	pred.wz_stat_mm   = sqrt(wz_nt2_mm_e1);
	pred.wz_stat_em   = sqrt(wz_nt2_em_e1);
	pred.wz_stat_ee   = sqrt(wz_nt2_ee_e1);

	pred.rare         = nt2_rare_mc_ee + nt2_rare_mc_mm + nt2_rare_mc_em;
	pred.rare_mm      = nt2_rare_mc_mm;
	pred.rare_em      = nt2_rare_mc_em;
	pred.rare_ee      = nt2_rare_mc_ee;
	pred.rare_err     = sqrt(RareESyst2*(nt2_rare_mc_ee + nt2_rare_mc_mm + nt2_rare_mc_em)*(nt2_rare_mc_ee + nt2_rare_mc_mm + nt2_rare_mc_em));
	pred.rare_err_mm  = sqrt(RareESyst2*nt2_rare_mc_mm*nt2_rare_mc_mm);
	pred.rare_err_em  = sqrt(RareESyst2*nt2_rare_mc_em*nt2_rare_mc_em);
	pred.rare_err_ee  = sqrt(RareESyst2*nt2_rare_mc_ee*nt2_rare_mc_ee);
	pred.rare_stat    = sqrt(nt2_rare_mc_ee_e1 + nt2_rare_mc_mm_e1 + nt2_rare_mc_em_e1);
	pred.rare_stat_mm = sqrt(nt2_rare_mc_mm_e1);
	pred.rare_stat_em = sqrt(nt2_rare_mc_em_e1);
	pred.rare_stat_ee = sqrt(nt2_rare_mc_ee_e1);
	
	
 	///////////////////////////////////////////////////////////////////////////////////
 	//  OUTPUT AS PLOT  ///////////////////////////////////////////////////////////////
 	///////////////////////////////////////////////////////////////////////////////////
 	TH1D    *h_obs        = new TH1D("h_observed",   "Observed number of events",  3, 0., 3.);
 	TH1D    *h_pred_sfake = new TH1D("h_pred_sfake", "Predicted single fakes", 3, 0., 3.);
 	TH1D    *h_pred_dfake = new TH1D("h_pred_dfake", "Predicted double fakes", 3, 0., 3.);
 	TH1D    *h_pred_chmid = new TH1D("h_pred_chmid", "Predicted charge mis id", 3, 0., 3.);
 	TH1D    *h_pred_mc    = new TH1D("h_pred_mc",    "Predicted Rare SM", 3, 0., 3.);
 	TH1D    *h_pred_ttw   = new TH1D("h_pred_ttw",   "Predicted WZ", 3, 0., 3.);
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
 	h_pred_ttw  ->SetLineWidth(1);
 	h_pred_sfake->SetLineColor(50);
 	h_pred_sfake->SetFillColor(50);
 	h_pred_dfake->SetLineColor(38);
 	h_pred_dfake->SetFillColor(38);
 	h_pred_chmid->SetLineColor(42);
 	h_pred_chmid->SetFillColor(42);
 	h_pred_mc   ->SetLineColor(31);
 	h_pred_mc   ->SetFillColor(31);
 	h_pred_ttw  ->SetLineColor(29);
 	h_pred_ttw  ->SetFillColor(29);
 
 	h_pred_tot  ->SetLineWidth(1);
 	// h_pred_tot  ->SetFillColor(kBlack);
 	// h_pred_tot  ->SetFillStyle(3013);
 	h_pred_tot  ->SetFillColor(12);
 	h_pred_tot  ->SetFillStyle(3005);
 	
 	// Add numbers:
 	h_obs->SetBinContent(1, nt2_ee);
 	h_obs->SetBinContent(2, nt2_mm);
 	h_obs->SetBinContent(3, nt2_em);
 	//h_obs->SetBinError(1, FR->getEStat(nt2_ee)); // FIXME
 	//h_obs->SetBinError(2, FR->getEStat(nt2_mm)); // FIXME
 	//h_obs->SetBinError(3, FR->getEStat(nt2_em)); // FIXME
 
 	TGraphAsymmErrors* gr_obs = FR->getGraphPoissonErrors( h_obs );
 	gr_obs->SetMarkerColor(kBlack);
 	gr_obs->SetMarkerStyle(20);
 	gr_obs->SetMarkerSize(2.5);
 	gr_obs->SetLineWidth(2);
 	gr_obs->SetLineColor(kBlack);
 	gr_obs->SetFillColor(kBlack);
 	
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
 	
 	h_pred_ttw->SetBinContent(1, wz_nt2_ee);
 	h_pred_ttw->SetBinContent(2, wz_nt2_mm);
 	h_pred_ttw->SetBinContent(3, wz_nt2_em);
 
 	h_pred_tot->Add(h_pred_sfake);
 	h_pred_tot->Add(h_pred_dfake);
 	h_pred_tot->Add(h_pred_chmid);
 	h_pred_tot->Add(h_pred_mc);
 	h_pred_tot->Add(h_pred_ttw);
 	h_pred_tot->SetBinError(1, sqrt(ee_tot_sqerr1 + ee_tot_sqerr2));
 	h_pred_tot->SetBinError(2, sqrt(mm_tot_sqerr1 + mm_tot_sqerr2));
 	h_pred_tot->SetBinError(3, sqrt(em_tot_sqerr1 + em_tot_sqerr2));
 	
 	hs_pred->Add(h_pred_sfake);
 	hs_pred->Add(h_pred_dfake);
 	hs_pred->Add(h_pred_chmid);
 	hs_pred->Add(h_pred_mc);
 	hs_pred->Add(h_pred_ttw);
 	
 	// double max = h_obs->Integral();
  	// double max = std::max(h_pred_tot->GetBinContent(1), h_pred_tot->GetBinContent(2));
  	// max = 1.7*std::max(max, h_pred_tot->GetBinContent(3));

 	double max = std::max(h_obs->GetBinContent(1), h_obs->GetBinContent(2));
 	max = 1.7*std::max(max, h_obs->GetBinContent(3));
 	
	if (reg == gRegion["HT80MET30bpp"]) max = 8.5; 
 	// if(reg == Baseline)    max = 125.;
 	// if(reg == HT80MET120)  max = 15.;
 	// if(reg == HT200MET120) max = 12.;
 	// if(reg == HT450MET0)   max = 11.;
 	// if(reg == HT450MET50)  max = 7.;
 	// if(reg == HT450MET120) max = 3.;
 	// if(reg == TTbarWSel3) max = 8.;
 	// if(reg == TTbarWSel1) max = 12.;
 	
 	h_obs       ->SetMaximum(max>1?max+1:1.);
 	h_pred_sfake->SetMaximum(max>1?max+1:1.);
 	h_pred_dfake->SetMaximum(max>1?max+1:1.);
 	h_pred_chmid->SetMaximum(max>1?max+1:1.);
 	h_pred_mc   ->SetMaximum(max>1?max+1:1.);
 	h_pred_ttw  ->SetMaximum(max>1?max+1:1.);
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
 	leg->AddEntry(h_pred_ttw,   "WZ Production","f");
 	leg->AddEntry(h_pred_tot,   "Total Uncertainty","f");
 	leg->SetFillStyle(0);
 	leg->SetTextFont(42);
 	// leg->SetTextSize(0.05);
 	leg->SetBorderSize(0);
 	
 	TCanvas *c_temp = new TCanvas("C_ObsPred", "Observed vs Predicted", 0, 0, 600, 600);
 	c_temp->cd();
 
 	hs_pred->Draw("hist");
 	h_pred_tot->DrawCopy("0 E2 same");
 	gr_obs->Draw("P same");
 	leg->Draw();
 	
 	drawRegionSel(reg);
 	drawTopLine(0.53, 0.8);
 	
 	gPad->RedrawAxis();
 	Util::PrintPDF(c_temp,   "ObsPred_" + gRegions[reg]->sname, fOutputDir + fOutputSubDir);
 	delete c_temp;	
 	delete h_obs, h_pred_sfake, h_pred_dfake, h_pred_chmid, h_pred_mc, h_pred_ttw, h_pred_tot, hs_pred;
 	delete gr_obs;
 	delete FR;

	return pred;
}
TTWZPrediction SSDLPlotter::makeIntPredictionTTW(TString filename, int reg){
	ofstream OUT(filename.Data(), ios::trunc);
	
	bool separateTTH = false;
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
	
	const float WZESyst  = 0.2;
	const float WZESyst2 = WZESyst*WZESyst;
	
	musamples = fMuData;
	elsamples = fEGData;
	emusamples = fMuEGData;
	
	OUT << "/////////////////////////////////////////////////////////////////////////////" << endl;
	OUT << " Producing integrated predictions for region " << gRegions[reg]->sname << endl;
	OUT << "  scaling MC to " << fLumiNorm << " /pb" << endl << endl;
	
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
		
		npp_mm += S->numbers[reg][Muon].npp;
		npf_mm += S->numbers[reg][Muon].npf + S->numbers[reg][Muon].nfp;
		nff_mm += S->numbers[reg][Muon].nff;
	}
	for(size_t i = 0; i < emusamples.size(); ++i){
		Sample *S = fSamples[emusamples[i]];
		nt2_em  += S->numbers[reg][ElMu].nt2;
		nt10_em += S->numbers[reg][ElMu].nt10;
		nt01_em += S->numbers[reg][ElMu].nt01;
		nt0_em  += S->numbers[reg][ElMu].nt0;
		
		npp_em += S->numbers[reg][ElMu].npp;
		npf_em += S->numbers[reg][ElMu].npf;
		nfp_em += S->numbers[reg][ElMu].nfp;
		nff_em += S->numbers[reg][ElMu].nff;
		
		nt2_em_BB_os += S->region[reg][HighPt].em.nt20_OS_BB_pt->GetEntries(); // ele in barrel
		nt2_em_EE_os += S->region[reg][HighPt].em.nt20_OS_EE_pt->GetEntries(); // ele in endcal
	}
	for(size_t i = 0; i < elsamples.size(); ++i){
		Sample *S = fSamples[elsamples[i]];
		nt2_ee  += S->numbers[reg][Elec].nt2;
		nt10_ee += S->numbers[reg][Elec].nt10;
		nt0_ee  += S->numbers[reg][Elec].nt0;
		
		npp_ee += S->numbers[reg][Elec].npp;
		npf_ee += S->numbers[reg][Elec].npf + S->numbers[reg][Elec].nfp;
		nff_ee += S->numbers[reg][Elec].nff;
		
		nt2_ee_BB_os += S->region[reg][HighPt].ee.nt20_OS_BB_pt->GetEntries(); // both in barrel
		nt2_ee_EE_os += S->region[reg][HighPt].ee.nt20_OS_EE_pt->GetEntries(); // both in endcal
		nt2_ee_EB_os += S->region[reg][HighPt].ee.nt20_OS_EB_pt->GetEntries(); // one barrel, one endcap
	}
	
	///////////////////////////////////////////////////////////////////////////////////
	// PRINTOUT ///////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////
	OUT << "---------------------------------------------------------------------------------------------------------" << endl;
	OUT << "         RATIOS  ||     Mu-fRatio      |     Mu-pRatio      ||     El-fRatio      |     El-pRatio      ||" << endl;
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
	
	// Background MC
	for(size_t i = 0; i < fMCBG.size(); ++i){
		Sample *S = fSamples[fMCBG[i]];
		float scale = fLumiNorm / S->getLumi();
		
		float temp_nt2_mm  = scale*S->region[reg][HighPt].mm.nt20_pt->Integral(0, getNFPtBins(Muon)+1); nt2sum_mm  += temp_nt2_mm ;
		float temp_nt10_mm = scale*S->region[reg][HighPt].mm.nt10_pt->Integral(0, getNFPtBins(Muon)+1); nt10sum_mm += temp_nt10_mm;
		float temp_nt0_mm  = scale*S->region[reg][HighPt].mm.nt00_pt->Integral(0, getNFPtBins(Muon)+1); nt0sum_mm  += temp_nt0_mm ;
		float temp_nt2_em  = scale*S->region[reg][HighPt].em.nt20_pt->Integral(0, getNFPtBins(ElMu)+1); nt2sum_em  += temp_nt2_em ;
		float temp_nt10_em = scale*S->region[reg][HighPt].em.nt10_pt->Integral(0, getNFPtBins(ElMu)+1); nt10sum_em += temp_nt10_em;
		float temp_nt01_em = scale*S->region[reg][HighPt].em.nt01_pt->Integral(0, getNFPtBins(ElMu)+1); nt01sum_em += temp_nt01_em;
		float temp_nt0_em  = scale*S->region[reg][HighPt].em.nt00_pt->Integral(0, getNFPtBins(ElMu)+1); nt0sum_em  += temp_nt0_em ;
		float temp_nt2_ee  = scale*S->region[reg][HighPt].ee.nt20_pt->Integral(0, getNFPtBins(Elec)+1); nt2sum_ee  += temp_nt2_ee ;
		float temp_nt10_ee = scale*S->region[reg][HighPt].ee.nt10_pt->Integral(0, getNFPtBins(Elec)+1); nt10sum_ee += temp_nt10_ee;
		float temp_nt0_ee  = scale*S->region[reg][HighPt].ee.nt00_pt->Integral(0, getNFPtBins(Elec)+1); nt0sum_ee  += temp_nt0_ee ;
		
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
	
	// Signal MC
	for(gSample i = sample_begin; i < gNSAMPLES; i=gSample(i+1)){
		Sample *S = fSamples[i];
		if(S->datamc != 2) continue;
		float scale = fLumiNorm / S->getLumi();
		
		float temp_nt2_mm  = scale*S->region[reg][HighPt].mm.nt20_pt->Integral(0, getNFPtBins(Muon)+1);
		float temp_nt10_mm = scale*S->region[reg][HighPt].mm.nt10_pt->Integral(0, getNFPtBins(Muon)+1);
		float temp_nt0_mm  = scale*S->region[reg][HighPt].mm.nt01_pt->Integral(0, getNFPtBins(Muon)+1);
		float temp_nt2_em  = scale*S->region[reg][HighPt].em.nt20_pt->Integral(0, getNFPtBins(ElMu)+1);
		float temp_nt10_em = scale*S->region[reg][HighPt].em.nt10_pt->Integral(0, getNFPtBins(ElMu)+1);
		float temp_nt01_em = scale*S->region[reg][HighPt].em.nt01_pt->Integral(0, getNFPtBins(ElMu)+1);
		float temp_nt0_em  = scale*S->region[reg][HighPt].em.nt00_pt->Integral(0, getNFPtBins(ElMu)+1);
		float temp_nt2_ee  = scale*S->region[reg][HighPt].ee.nt20_pt->Integral(0, getNFPtBins(Elec)+1);
		float temp_nt10_ee = scale*S->region[reg][HighPt].ee.nt10_pt->Integral(0, getNFPtBins(Elec)+1);
		float temp_nt0_ee  = scale*S->region[reg][HighPt].ee.nt01_pt->Integral(0, getNFPtBins(Elec)+1);
		
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
	
	nF_mm = npf_mm + nff_mm;
	nF_em = npf_em+nfp_em+nff_em;
	nF_ee = npf_ee+nff_ee;
	nSF   = npf_mm + npf_em + nfp_em + npf_ee;
	nDF   = nff_mm + nff_em + nff_ee;
	nF    = nF_mm + nF_em + nF_ee;
	
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
	float fbb(0.),fee(0.),feb(0.);
	float fbbE(0.),feeE(0.),febE(0.);
	float fbb_mc(0.),fee_mc(0.),feb_mc(0.);
	float fbbE_mc(0.),feeE_mc(0.),febE_mc(0.);
	
	calculateChMisIdProb(fEGData, BB, fbb, fbbE);
	calculateChMisIdProb(fEGData, EB, feb, febE);
	calculateChMisIdProb(fEGData, EE, fee, feeE);
	
	calculateChMisIdProb(fMCBG, BB, fbb_mc, fbbE_mc);
	calculateChMisIdProb(fMCBG, EB, feb_mc, febE_mc);
	calculateChMisIdProb(fMCBG, EE, fee_mc, feeE_mc);
	
	// Simple error propagation assuming error on number of events is sqrt(N)
	nt2_ee_chmid    = 2*fbb*nt2_ee_BB_os                           + 2*fee*nt2_ee_EE_os                      + 2*feb*nt2_ee_EB_os;
	nt2_ee_chmid_e1 = sqrt( 4*fbb*fbb*FR->getEStat2(nt2_ee_BB_os)  + 4*fee*fee*FR->getEStat2(nt2_ee_EE_os)   + 4*feb*feb*FR->getEStat2(nt2_ee_EB_os) ); // stat only
	nt2_ee_chmid_e2 = sqrt( nt2_ee_BB_os*nt2_ee_BB_os*4*fbbE*fbbE  + 4*feeE*feeE*nt2_ee_EE_os*nt2_ee_EE_os   + 4*febE*febE*nt2_ee_EB_os*nt2_ee_EB_os ); // syst only
	
	nt2_em_chmid    = fbb*nt2_em_BB_os + fee*nt2_em_EE_os;
	nt2_em_chmid_e1 = sqrt( fbb*fbb*FR->getEStat2(nt2_em_BB_os) + fee*fee*FR->getEStat2(nt2_em_EE_os) );
	nt2_em_chmid_e2 = sqrt( nt2_em_BB_os*nt2_em_BB_os*fbbE*fbbE + nt2_em_EE_os*nt2_em_EE_os*feeE*feeE );

 	///////////////////////////////////////////////////////////////////////////////////
 	// PRINTOUT ///////////////////////////////////////////////////////////////////////
 	///////////////////////////////////////////////////////////////////////////////////
 	OUT << "-------------------------------------------------------------------------------------------" << endl;
 	OUT << "       E-ChMisID  ||     Barrel-Barrel    |    Barrel - EndCap   |    Endcap - EndCap    ||" << endl;
 	OUT << "-------------------------------------------------------------------------------------------" << endl;
 	OUT << "  Data            ||";
 	OUT << setw(8) << setprecision(2) << fbb  << " +/- " << setw(8) << setprecision(2) << fbbE  << " |";
 	OUT << setw(8) << setprecision(2) << feb  << " +/- " << setw(8) << setprecision(2) << febE  << " |";
 	OUT << setw(8) << setprecision(2) << fee  << " +/- " << setw(8) << setprecision(2) << feeE  << "  ||" << endl;
 	OUT << "  MC              ||";
 	OUT << setw(8) << setprecision(2) << fbb_mc  << " +/- " << setw(8) << setprecision(2) << fbbE_mc  << " |";
 	OUT << setw(8) << setprecision(2) << feb_mc  << " +/- " << setw(8) << setprecision(2) << febE_mc  << " |";
 	OUT << setw(8) << setprecision(2) << fee_mc  << " +/- " << setw(8) << setprecision(2) << feeE_mc  << "  ||";
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
		float scale = fLumiNorm / S->getLumi();
		
		mc_os_em_bb_sum += scale*S->region[reg][HighPt].em.nt20_OS_BB_pt->Integral(0, getNFPtBins(Elec)+1);
		mc_os_em_ee_sum += scale*S->region[reg][HighPt].em.nt20_OS_EE_pt->Integral(0, getNFPtBins(Elec)+1);
		mc_os_ee_bb_sum += scale*S->region[reg][HighPt].ee.nt20_OS_BB_pt->Integral(0, getNFPtBins(Elec)+1);
		mc_os_ee_eb_sum += scale*S->region[reg][HighPt].ee.nt20_OS_EB_pt->Integral(0, getNFPtBins(Elec)+1);
		mc_os_ee_ee_sum += scale*S->region[reg][HighPt].ee.nt20_OS_EE_pt->Integral(0, getNFPtBins(Elec)+1);
		
		OUT << setw(16) << S->sname << " || ";
		OUT << setw(7)  << setprecision(2) << scale*S->region[reg][HighPt].em.nt20_OS_BB_pt->Integral(0, getNFPtBins(Elec)+1) << " | ";
		OUT << setw(7)  << setprecision(2) << scale*S->region[reg][HighPt].em.nt20_OS_EE_pt->Integral(0, getNFPtBins(Elec)+1) << " || ";
		OUT << setw(7)  << setprecision(2) << scale*S->region[reg][HighPt].ee.nt20_OS_BB_pt->Integral(0, getNFPtBins(Elec)+1) << " | ";
		OUT << setw(7)  << setprecision(2) << scale*S->region[reg][HighPt].ee.nt20_OS_EB_pt->Integral(0, getNFPtBins(Elec)+1) << " | ";
		OUT << setw(7)  << setprecision(2) << scale*S->region[reg][HighPt].ee.nt20_OS_EE_pt->Integral(0, getNFPtBins(Elec)+1) << " || ";
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
	OUT << setw(7) << Form("%6.4f", 0.5*(feb+fbb) * nt2_em_BB_os ) << " | ";
	OUT << setw(7) << Form("%6.4f", 0.5*(feb+fee) * nt2_em_EE_os ) << " || ";
	OUT << setw(7) << Form("%6.4f",     fbb       * nt2_ee_BB_os ) << " | ";
	OUT << setw(7) << Form("%6.4f",     feb       * nt2_ee_EB_os ) << " | ";
	OUT << setw(7) << Form("%6.4f",     fee       * nt2_ee_EE_os ) << " || ";
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
	
	vector<int> mcbkg;
	//mcbkg.push_back(WZ);
	mcbkg.push_back(ZZ);
// MARC	mcbkg.push_back(GVJets);
	mcbkg.push_back(DPSWW);
	if (!separateTTH) mcbkg.push_back(TTbarH);
	mcbkg.push_back(TTbarG);
	mcbkg.push_back(TbZ);
	mcbkg.push_back(WpWp);
	mcbkg.push_back(WmWm);
	mcbkg.push_back(WWZ);
	mcbkg.push_back(WZZ);
	mcbkg.push_back(WWG);
	mcbkg.push_back(WWW);
	mcbkg.push_back(ZZZ);

	for(size_t i = 0; i < mcbkg.size(); ++i){
		Sample *S = fSamples[mcbkg[i]];
		float scale = fLumiNorm/S->getLumi();
		
		float temp_nt2_mm = gMMTrigScale*scale*S->region[reg][HighPt].mm.nt20_pt->Integral(0, getNFPtBins(Muon)+1);
		float temp_nt2_em = gEMTrigScale*scale*S->region[reg][HighPt].em.nt20_pt->Integral(0, getNFPtBins(ElMu)+1);
		float temp_nt2_ee = gEETrigScale*scale*S->region[reg][HighPt].ee.nt20_pt->Integral(0, getNFPtBins(Elec)+1);
		
		nt2_rare_mc_mm += temp_nt2_mm;
		nt2_rare_mc_em += temp_nt2_em;
		nt2_rare_mc_ee += temp_nt2_ee;
		
		nt2_rare_mc_mm_e1 += gMMTrigScale*gMMTrigScale*scale*scale * S->numbers[reg][Muon].tt_avweight*S->numbers[reg][Muon].tt_avweight * S->getError2(S->region[reg][HighPt].mm.nt20_pt->GetEntries()); // for stat error take actual entries, not pileup weighted integral...
		nt2_rare_mc_em_e1 += gEMTrigScale*gEMTrigScale*scale*scale * S->numbers[reg][ElMu].tt_avweight*S->numbers[reg][ElMu].tt_avweight * S->getError2(S->region[reg][HighPt].em.nt20_pt->GetEntries());
		nt2_rare_mc_ee_e1 += gEETrigScale*gEETrigScale*scale*scale * S->numbers[reg][Elec].tt_avweight*S->numbers[reg][Elec].tt_avweight * S->getError2(S->region[reg][HighPt].ee.nt20_pt->GetEntries());
		
		OUT << Form("%16s || %5.2f ± %5.2f         || %5.2f ± %5.2f         || %5.2f ± %5.2f         ||\n", S->sname.Data(),
					temp_nt2_mm, gMMTrigScale*scale*S->numbers[reg][Muon].tt_avweight * S->getError(S->region[reg][HighPt].mm.nt20_pt->GetEntries()),
					temp_nt2_em, gEMTrigScale*scale*S->numbers[reg][ElMu].tt_avweight * S->getError(S->region[reg][HighPt].em.nt20_pt->GetEntries()),
					temp_nt2_ee, gEETrigScale*scale*S->numbers[reg][Elec].tt_avweight * S->getError(S->region[reg][HighPt].ee.nt20_pt->GetEntries()));
	}
	OUT << "----------------------------------------------------------------------------------------------" << endl;
	OUT << Form("%16s || %5.2f ± %5.2f ± %5.2f || %5.2f ± %5.2f ± %5.2f || %5.2f ± %5.2f ± %5.2f ||\n", "Rare SM (Sum)",
				nt2_rare_mc_mm, sqrt(nt2_rare_mc_mm_e1), RareESyst*nt2_rare_mc_mm,
				nt2_rare_mc_em, sqrt(nt2_rare_mc_em_e1), RareESyst*nt2_rare_mc_em,
				nt2_rare_mc_ee, sqrt(nt2_rare_mc_ee_e1), RareESyst*nt2_rare_mc_ee);
	OUT << "----------------------------------------------------------------------------------------------" << endl;
	
	///////////////////////////////////////////
	// WZ production
	float wzscale = fLumiNorm/fSamples[WZ]->getLumi();
	float wz_nt2_mm = gMMTrigScale*wzscale*fSamples[WZ]->region[reg][HighPt].mm.nt20_pt->Integral(0, getNFPtBins(Muon)+1);
	float wz_nt2_em = gEMTrigScale*wzscale*fSamples[WZ]->region[reg][HighPt].em.nt20_pt->Integral(0, getNFPtBins(ElMu)+1);
	float wz_nt2_ee = gEETrigScale*wzscale*fSamples[WZ]->region[reg][HighPt].ee.nt20_pt->Integral(0, getNFPtBins(Elec)+1);
	
	float wz_nt2_mm_e1 = gMMTrigScale*gMMTrigScale*wzscale*wzscale * fSamples[WZ]->numbers[reg][Muon].tt_avweight*fSamples[WZ]->numbers[reg][Muon].tt_avweight * fSamples[WZ]->getError2(fSamples[WZ]->region[reg][HighPt].mm.nt20_pt->GetEntries()); // for stat error take actual entries, not pileup weighted integral...
	float wz_nt2_em_e1 = gEMTrigScale*gEMTrigScale*wzscale*wzscale * fSamples[WZ]->numbers[reg][ElMu].tt_avweight*fSamples[WZ]->numbers[reg][ElMu].tt_avweight * fSamples[WZ]->getError2(fSamples[WZ]->region[reg][HighPt].em.nt20_pt->GetEntries());
	float wz_nt2_ee_e1 = gEETrigScale*gEETrigScale*wzscale*wzscale * fSamples[WZ]->numbers[reg][Elec].tt_avweight*fSamples[WZ]->numbers[reg][Elec].tt_avweight * fSamples[WZ]->getError2(fSamples[WZ]->region[reg][HighPt].ee.nt20_pt->GetEntries());
	
	OUT << Form("%16s || %5.2f ± %5.2f         || %5.2f ± %5.2f         || %5.2f ± %5.2f         ||\n", "WZ Prod.",
				wz_nt2_mm, sqrt(wz_nt2_mm_e1),
				wz_nt2_em, sqrt(wz_nt2_em_e1),
				wz_nt2_ee, sqrt(wz_nt2_ee_e1));
	OUT << "----------------------------------------------------------------------------------------------" << endl;
	
	
	///////////////////////////////////////////
	// ttX production
	float ttwscale = fLumiNorm/fSamples[TTbarW]->getLumi();
	float ttw_nt2_mm = gMMTrigScale*ttwscale*fSamples[TTbarW]->region[reg][HighPt].mm.nt20_pt->Integral(0, getNFPtBins(Muon)+1);
	float ttw_nt2_em = gEMTrigScale*ttwscale*fSamples[TTbarW]->region[reg][HighPt].em.nt20_pt->Integral(0, getNFPtBins(ElMu)+1);
	float ttw_nt2_ee = gEETrigScale*ttwscale*fSamples[TTbarW]->region[reg][HighPt].ee.nt20_pt->Integral(0, getNFPtBins(Elec)+1);
	
	float ttw_nt2_mm_e1 = gMMTrigScale*gMMTrigScale*ttwscale*ttwscale * fSamples[TTbarW]->numbers[reg][Muon].tt_avweight*fSamples[TTbarW]->numbers[reg][Muon].tt_avweight * fSamples[TTbarW]->getError2(fSamples[TTbarW]->region[reg][HighPt].mm.nt20_pt->GetEntries()); // for stat error take actual entries, not pileup weighted integral...
	float ttw_nt2_em_e1 = gEMTrigScale*gEMTrigScale*ttwscale*ttwscale * fSamples[TTbarW]->numbers[reg][ElMu].tt_avweight*fSamples[TTbarW]->numbers[reg][ElMu].tt_avweight * fSamples[TTbarW]->getError2(fSamples[TTbarW]->region[reg][HighPt].em.nt20_pt->GetEntries());
	float ttw_nt2_ee_e1 = gEETrigScale*gEETrigScale*ttwscale*ttwscale * fSamples[TTbarW]->numbers[reg][Elec].tt_avweight*fSamples[TTbarW]->numbers[reg][Elec].tt_avweight * fSamples[TTbarW]->getError2(fSamples[TTbarW]->region[reg][HighPt].ee.nt20_pt->GetEntries());
	
	float ttzscale = fLumiNorm/fSamples[TTbarZ]->getLumi();
	float ttz_nt2_mm = gMMTrigScale*ttzscale*fSamples[TTbarZ]->region[reg][HighPt].mm.nt20_pt->Integral(0, getNFPtBins(Muon)+1);
	float ttz_nt2_em = gEMTrigScale*ttzscale*fSamples[TTbarZ]->region[reg][HighPt].em.nt20_pt->Integral(0, getNFPtBins(ElMu)+1);
	float ttz_nt2_ee = gEETrigScale*ttzscale*fSamples[TTbarZ]->region[reg][HighPt].ee.nt20_pt->Integral(0, getNFPtBins(Elec)+1);
	
	float ttz_nt2_mm_e1 = gMMTrigScale*gMMTrigScale*ttzscale*ttzscale * fSamples[TTbarZ]->numbers[reg][Muon].tt_avweight*fSamples[TTbarZ]->numbers[reg][Muon].tt_avweight * fSamples[TTbarZ]->getError2(fSamples[TTbarZ]->region[reg][HighPt].mm.nt20_pt->GetEntries()); // for stat error take actual entries, not pileup weighted integral...
	float ttz_nt2_em_e1 = gEMTrigScale*gEMTrigScale*ttzscale*ttzscale * fSamples[TTbarZ]->numbers[reg][ElMu].tt_avweight*fSamples[TTbarZ]->numbers[reg][ElMu].tt_avweight * fSamples[TTbarZ]->getError2(fSamples[TTbarZ]->region[reg][HighPt].em.nt20_pt->GetEntries());
	float ttz_nt2_ee_e1 = gEETrigScale*gEETrigScale*ttzscale*ttzscale * fSamples[TTbarZ]->numbers[reg][Elec].tt_avweight*fSamples[TTbarZ]->numbers[reg][Elec].tt_avweight * fSamples[TTbarZ]->getError2(fSamples[TTbarZ]->region[reg][HighPt].ee.nt20_pt->GetEntries());

	float tthscale = fLumiNorm/fSamples[TTbarH]->getLumi();
	float tth_nt2_mm = gMMTrigScale*tthscale*fSamples[TTbarH]->region[reg][HighPt].mm.nt20_pt->Integral(0, getNFPtBins(Muon)+1);
	float tth_nt2_em = gEMTrigScale*tthscale*fSamples[TTbarH]->region[reg][HighPt].em.nt20_pt->Integral(0, getNFPtBins(ElMu)+1);
	float tth_nt2_ee = gEETrigScale*tthscale*fSamples[TTbarH]->region[reg][HighPt].ee.nt20_pt->Integral(0, getNFPtBins(Elec)+1);

	float tth_nt2_mm_e1 = gMMTrigScale*gMMTrigScale*tthscale*tthscale * fSamples[TTbarH]->numbers[reg][Muon].tt_avweight*fSamples[TTbarH]->numbers[reg][Muon].tt_avweight * fSamples[TTbarH]->getError2(fSamples[TTbarH]->region[reg][HighPt].mm.nt20_pt->GetEntries()); // for stat error take actual entries, not pileup weighted integral...
	float tth_nt2_em_e1 = gEMTrigScale*gEMTrigScale*tthscale*tthscale * fSamples[TTbarH]->numbers[reg][ElMu].tt_avweight*fSamples[TTbarH]->numbers[reg][ElMu].tt_avweight * fSamples[TTbarH]->getError2(fSamples[TTbarH]->region[reg][HighPt].em.nt20_pt->GetEntries());
	float tth_nt2_ee_e1 = gEETrigScale*gEETrigScale*tthscale*tthscale * fSamples[TTbarH]->numbers[reg][Elec].tt_avweight*fSamples[TTbarH]->numbers[reg][Elec].tt_avweight * fSamples[TTbarH]->getError2(fSamples[TTbarH]->region[reg][HighPt].ee.nt20_pt->GetEntries());

	OUT << Form("%16s || %5.2f ± %5.2f         || %5.2f ± %5.2f         || %5.2f ± %5.2f         ||\n", "ttW Prod.",
				ttw_nt2_mm, sqrt(ttw_nt2_mm_e1),
				ttw_nt2_em, sqrt(ttw_nt2_em_e1),
				ttw_nt2_ee, sqrt(ttw_nt2_ee_e1));
	
	OUT << Form("%16s || %5.2f ± %5.2f         || %5.2f ± %5.2f         || %5.2f ± %5.2f         ||\n", "ttZ Prod.",
	ttz_nt2_mm, sqrt(ttz_nt2_mm_e1),
	ttz_nt2_em, sqrt(ttz_nt2_em_e1),
	ttz_nt2_ee, sqrt(ttz_nt2_ee_e1));

	if (separateTTH) {
		OUT << Form("%16s || %5.2f ± %5.2f         || %5.2f ± %5.2f         || %5.2f ± %5.2f         ||\n", "ttH Prod.",
					tth_nt2_mm, sqrt(tth_nt2_mm_e1),
					tth_nt2_em, sqrt(tth_nt2_em_e1),
					tth_nt2_ee, sqrt(tth_nt2_ee_e1));
	}
	OUT << "----------------------------------------------------------------------------------------------" << endl;
	// Just add different errors in quadrature (they are independent)
	float mm_tot_sqerr1 = FR->getMMTotEStat()*FR->getMMTotEStat()                                   + nt2_rare_mc_mm_e1 + wz_nt2_mm_e1;
	float em_tot_sqerr1 = FR->getEMTotEStat()*FR->getEMTotEStat() + nt2_em_chmid_e1*nt2_em_chmid_e1 + nt2_rare_mc_em_e1 + wz_nt2_em_e1;
	float ee_tot_sqerr1 = FR->getEETotEStat()*FR->getEETotEStat() + nt2_ee_chmid_e1*nt2_ee_chmid_e1 + nt2_rare_mc_ee_e1 + wz_nt2_ee_e1;
	float mm_tot_sqerr2 = nF_mm*nF_mm*FakeESyst2                                   + RareESyst2*nt2_rare_mc_mm*nt2_rare_mc_mm + WZESyst2*wz_nt2_mm*wz_nt2_mm;
	float em_tot_sqerr2 = nF_em*nF_em*FakeESyst2 + nt2_em_chmid_e2*nt2_em_chmid_e2 + RareESyst2*nt2_rare_mc_em*nt2_rare_mc_em + WZESyst2*wz_nt2_em*wz_nt2_em;
	float ee_tot_sqerr2 = nF_ee*nF_ee*FakeESyst2 + nt2_ee_chmid_e2*nt2_ee_chmid_e2 + RareESyst2*nt2_rare_mc_ee*nt2_rare_mc_ee + WZESyst2*wz_nt2_ee*wz_nt2_ee;
	OUT << Form("%16s || %5.2f ± %5.2f ± %5.2f || %5.2f ± %5.2f ± %5.2f || %5.2f ± %5.2f ± %5.2f ||\n", "tot. backgr.",
				nF_mm                + nt2_rare_mc_mm + wz_nt2_mm, sqrt(mm_tot_sqerr1), sqrt(mm_tot_sqerr2),
				nF_em + nt2_em_chmid + nt2_rare_mc_em + wz_nt2_em, sqrt(em_tot_sqerr1), sqrt(em_tot_sqerr2),
				nF_ee + nt2_ee_chmid + nt2_rare_mc_ee + wz_nt2_ee, sqrt(ee_tot_sqerr1), sqrt(ee_tot_sqerr2));
	OUT << Form("%16s || %5.2f ± %5.2f         || %5.2f ± %5.2f         || %5.2f ± %5.2f         ||\n", "",
				nF_mm                + nt2_rare_mc_mm + wz_nt2_mm, sqrt(mm_tot_sqerr1 + mm_tot_sqerr2),
				nF_em + nt2_em_chmid + nt2_rare_mc_em + wz_nt2_em, sqrt(em_tot_sqerr1 + em_tot_sqerr2),
				nF_ee + nt2_ee_chmid + nt2_rare_mc_ee + wz_nt2_ee, sqrt(ee_tot_sqerr1 + ee_tot_sqerr2));
	OUT << "----------------------------------------------------------------------------------------------" << endl;
	OUT << Form("%16s || %5.2f                 || %5.2f                 || %5.2f                 ||\n", "tot. MC", nt2sum_mm, nt2sum_em, nt2sum_ee);
	OUT << "==============================================================================================" << endl;
	OUT << Form("%16s || %2.0f                    || %2.0f                    || %2.0f                    ||\n", "observed", nt2_mm, nt2_em, nt2_ee);
	OUT << "==============================================================================================" << endl;
	OUT << "        predicted: ";
	float tot_pred        = nF + nt2_rare_mc_mm + wz_nt2_mm + nt2_em_chmid + nt2_rare_mc_em + wz_nt2_em + nt2_ee_chmid + nt2_rare_mc_ee + wz_nt2_ee;
	float comb_tot_sqerr1 = FR->getTotEStat()*FR->getTotEStat() + nt2_rare_mc_mm_e1 + nt2_em_chmid_e1*nt2_em_chmid_e1 + nt2_rare_mc_em_e1 + nt2_ee_chmid_e1*nt2_ee_chmid_e1 + nt2_rare_mc_ee_e1 + wz_nt2_mm_e1 + wz_nt2_em_e1 + wz_nt2_ee_e1;
	float comb_tot_sqerr2 = nF*nF*FakeESyst2 + RareESyst2*(nt2_rare_mc_mm + nt2_rare_mc_em + nt2_rare_mc_ee)*(nt2_rare_mc_mm + nt2_rare_mc_em + nt2_rare_mc_ee) + nt2_em_chmid_e2*nt2_em_chmid_e2 + nt2_ee_chmid_e2*nt2_ee_chmid_e2 + WZESyst2*(wz_nt2_mm + wz_nt2_em + wz_nt2_ee)*(wz_nt2_mm + wz_nt2_em + wz_nt2_ee);
	// FIXME: Why take 50% on Rare yields on SUM and not on individual channels?
	OUT << setw(5) << left << Form("%5.2f", tot_pred ) << " ± ";
	OUT << setw(5) << Form("%5.2f", sqrt(comb_tot_sqerr1)) << " ± ";
	OUT << setw(5) << Form("%5.2f", sqrt(comb_tot_sqerr2)) << endl;
	OUT << "      combined MC: ";
	OUT << setw(5) << left << Form("%5.2f", nt2sum_mm+nt2sum_em+nt2sum_ee ) << endl;
	OUT << "combined observed: ";
	OUT << setw(5) << left << Form("%2.0f", nt2_mm+nt2_em+nt2_ee ) << endl;
	OUT << "==============================================================================================" << endl;
	OUT.close();
	
	///////////////////////////////////////////////////////////////////////////////////
	//  OUTPUT FOR ANALYSIS NOTE  /////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////
	fOUTSTREAM3 << "%% " + gRegions[reg]->sname << endl;
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
	fOUTSTREAM3 << Form("WZ:  %6.1f ± %6.1f  ( %5.1f±%5.1f | %5.1f±%5.1f | %5.1f±%5.1f )\n",
						wz_nt2_mm + wz_nt2_em + wz_nt2_ee,
						sqrt(wz_nt2_mm_e1 + wz_nt2_em_e1 + wz_nt2_em_e1 + WZESyst2*(wz_nt2_mm + wz_nt2_em + wz_nt2_ee)*(wz_nt2_mm + wz_nt2_em + wz_nt2_ee)),
						wz_nt2_ee, sqrt(wz_nt2_ee_e1 + WZESyst2*wz_nt2_ee*wz_nt2_ee),
						wz_nt2_mm, sqrt(wz_nt2_mm_e1 + WZESyst2*wz_nt2_mm*wz_nt2_mm),
						wz_nt2_em, sqrt(wz_nt2_em_e1 + WZESyst2*wz_nt2_em*wz_nt2_em));
	fOUTSTREAM3 << Form("Tot: %6.1f ± %6.1f  ( %5.1f±%5.1f | %5.1f±%5.1f | %5.1f±%5.1f )\n",
						tot_pred, sqrt(comb_tot_sqerr1 + comb_tot_sqerr2),
						nF_ee + nt2_rare_mc_ee + nt2_ee_chmid, sqrt(ee_tot_sqerr1 + ee_tot_sqerr2),
						nF_mm + nt2_rare_mc_mm               , sqrt(mm_tot_sqerr1 + mm_tot_sqerr2),
						nF_em + nt2_rare_mc_em + nt2_em_chmid, sqrt(em_tot_sqerr1 + em_tot_sqerr2));
	fOUTSTREAM3 << "-----------------------------------------------------------------" << endl;
	fOUTSTREAM3 << Form("Obs: %4.0f             ( %3.0f         | %3.0f         | %3.0f         )\n", nt2_mm+nt2_em+nt2_ee, nt2_ee, nt2_mm, nt2_em);
	fOUTSTREAM3 << "-----------------------------------------------------------------" << endl;
	fOUTSTREAM3 << Form("ttW: %6.1f ± %6.1f  ( %5.1f±%5.1f | %5.1f±%5.1f | %5.1f±%5.1f )\n",
						ttw_nt2_mm + ttw_nt2_em + ttw_nt2_ee,
						sqrt(ttw_nt2_mm_e1 + ttw_nt2_em_e1 + ttw_nt2_em_e1 + RareESyst2*(ttw_nt2_mm + ttw_nt2_em + ttw_nt2_ee)*(ttw_nt2_mm + ttw_nt2_em + ttw_nt2_ee)),
						ttw_nt2_ee, sqrt(ttw_nt2_ee_e1 + RareESyst2*ttw_nt2_ee*ttw_nt2_ee),
						ttw_nt2_mm, sqrt(ttw_nt2_mm_e1 + RareESyst2*ttw_nt2_mm*ttw_nt2_mm),
						ttw_nt2_em, sqrt(ttw_nt2_em_e1 + RareESyst2*ttw_nt2_em*ttw_nt2_em));
	fOUTSTREAM3 << Form("ttZ: %6.1f ± %6.1f  ( %5.1f±%5.1f | %5.1f±%5.1f | %5.1f±%5.1f )\n",
	ttz_nt2_mm + ttz_nt2_em + ttz_nt2_ee,
	sqrt(ttz_nt2_mm_e1 + ttz_nt2_em_e1 + ttz_nt2_em_e1 + RareESyst2*(ttz_nt2_mm + ttz_nt2_em + ttz_nt2_ee)*(ttz_nt2_mm + ttz_nt2_em + ttz_nt2_ee)),
	ttz_nt2_ee, sqrt(ttz_nt2_ee_e1 + RareESyst2*ttz_nt2_ee*ttz_nt2_ee),
	ttz_nt2_mm, sqrt(ttz_nt2_mm_e1 + RareESyst2*ttz_nt2_mm*ttz_nt2_mm),
	ttz_nt2_em, sqrt(ttz_nt2_em_e1 + RareESyst2*ttz_nt2_em*ttz_nt2_em));
	if (separateTTH) {
		fOUTSTREAM3 << Form("tth: %6.1f ± %6.1f  ( %5.1f±%5.1f | %5.1f±%5.1f | %5.1f±%5.1f )\n",
							tth_nt2_mm + tth_nt2_em + tth_nt2_ee,
							sqrt(tth_nt2_mm_e1 + tth_nt2_em_e1 + tth_nt2_em_e1 + RareESyst2*(tth_nt2_mm + tth_nt2_em + tth_nt2_ee)*(tth_nt2_mm + tth_nt2_em + tth_nt2_ee)),
							tth_nt2_ee, sqrt(tth_nt2_ee_e1 + RareESyst2*tth_nt2_ee*tth_nt2_ee),
							tth_nt2_mm, sqrt(tth_nt2_mm_e1 + RareESyst2*tth_nt2_mm*tth_nt2_mm),
							tth_nt2_em, sqrt(tth_nt2_em_e1 + RareESyst2*tth_nt2_em*tth_nt2_em));
		fOUTSTREAM3 << "-----------------------------------------------------------------" << endl;
	}
	fOUTSTREAM3 << endl;
	
	///////////////////////////////////////////////////////////////////////////////////
	//  OUTPUT FOR AN TABLE  /////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////
	fOUTSTREAM << gRegions[reg]->sname << endl;
	fOUTSTREAM << Form("Double Fakes   & %5.1f &$\\pm$ %5.1f & %5.1f &$\\pm$ %5.1f & %5.1f &$\\pm$ %5.1f & %5.1f &$\\pm$ %5.1f \\\\ \n",
	nff_mm, sqrt(FR->getMMNffEStat()*FR->getMMNffEStat()+nff_mm*nff_mm*FakeESyst2),
	nff_em, sqrt(FR->getEMNffEStat()*FR->getEMNffEStat()+nff_em*nff_em*FakeESyst2),
	nff_ee, sqrt(FR->getEENffEStat()*FR->getEENffEStat()+nff_ee*nff_ee*FakeESyst2),
	nff_em + nff_mm + nff_ee, sqrt(FR->getTotDoubleEStat()*FR->getTotDoubleEStat() + nDF*nDF*FakeESyst2));
	fOUTSTREAM << Form("Single Fakes   & %5.1f &$\\pm$ %5.1f & %5.1f &$\\pm$ %5.1f & %5.1f &$\\pm$ %5.1f & %5.1f &$\\pm$ %5.1f \\\\ \n",
	npf_mm,          sqrt(FR->getMMNpfEStat()   *FR->getMMNpfEStat()    +  npf_mm*npf_mm*FakeESyst2),
	npf_em + nfp_em, sqrt(FR->getEMSingleEStat()*FR->getEMSingleEStat() + (npf_em+nfp_em)*(npf_em+nfp_em)*FakeESyst2),
	npf_ee,          sqrt(FR->getEENpfEStat()   *FR->getEENpfEStat()    +  npf_ee*npf_ee*FakeESyst2),
	npf_em + nfp_em + npf_mm + npf_ee, sqrt(FR->getTotSingleEStat()*FR->getTotSingleEStat() + nSF*nSF*FakeESyst2));
	fOUTSTREAM << Form("Charge MisID   & \\multicolumn{2}{c|}{-} & %5.1f &$\\pm$ %5.1f & %5.1f &$\\pm$ %5.1f & %5.1f &$\\pm$ %5.1f \\\\ \n",
	nt2_em_chmid, sqrt(nt2_em_chmid_e1*nt2_em_chmid_e1 + nt2_em_chmid_e2*nt2_em_chmid_e2),
	nt2_ee_chmid, sqrt(nt2_ee_chmid_e1*nt2_ee_chmid_e1 + nt2_ee_chmid_e2*nt2_ee_chmid_e2),
	nt2_ee_chmid + nt2_em_chmid, sqrt(nt2_ee_chmid_e1*nt2_ee_chmid_e1 + nt2_ee_chmid_e2*nt2_ee_chmid_e2 + nt2_em_chmid_e1*nt2_em_chmid_e1 + nt2_em_chmid_e2*nt2_em_chmid_e2));
	fOUTSTREAM << Form("Rare SM        & %5.1f &$\\pm$ %5.1f & %5.1f &$\\pm$ %5.1f & %5.1f &$\\pm$ %5.1f & %5.1f &$\\pm$ %5.1f \\\\ \n",
	nt2_rare_mc_mm, sqrt(nt2_rare_mc_mm_e1 + RareESyst2*nt2_rare_mc_mm*nt2_rare_mc_mm),
	nt2_rare_mc_em, sqrt(nt2_rare_mc_em_e1 + RareESyst2*nt2_rare_mc_em*nt2_rare_mc_em),
	nt2_rare_mc_ee, sqrt(nt2_rare_mc_ee_e1 + RareESyst2*nt2_rare_mc_ee*nt2_rare_mc_ee),
	nt2_rare_mc_ee + nt2_rare_mc_mm + nt2_rare_mc_em, sqrt(nt2_rare_mc_ee_e1 + nt2_rare_mc_mm_e1 + nt2_rare_mc_em_e1 + RareESyst2*(nt2_rare_mc_ee + nt2_rare_mc_mm + nt2_rare_mc_em)*(nt2_rare_mc_ee + nt2_rare_mc_mm + nt2_rare_mc_em)));
	fOUTSTREAM << Form("WZ Prod.       & %5.1f &$\\pm$ %5.1f & %5.1f &$\\pm$ %5.1f & %5.1f &$\\pm$ %5.1f & %5.1f &$\\pm$ %5.1f \\\\ \\hline \n",
	wz_nt2_mm, sqrt(wz_nt2_mm_e1 + WZESyst2*wz_nt2_mm*wz_nt2_mm),
	wz_nt2_em, sqrt(wz_nt2_em_e1 + WZESyst2*wz_nt2_em*wz_nt2_em),
	wz_nt2_ee, sqrt(wz_nt2_ee_e1 + WZESyst2*wz_nt2_ee*wz_nt2_ee),
	wz_nt2_ee + wz_nt2_mm + wz_nt2_em, sqrt(wz_nt2_mm_e1 + wz_nt2_ee_e1 + wz_nt2_em_e1 + WZESyst2*(wz_nt2_ee + wz_nt2_mm + wz_nt2_em)*(wz_nt2_ee + wz_nt2_mm + wz_nt2_em)));
	fOUTSTREAM << Form("Total Bkg      & %5.1f &$\\pm$ %5.1f & %5.1f &$\\pm$ %5.1f & %5.1f &$\\pm$ %5.1f & %5.1f &$\\pm$ %5.1f \\\\ \\hline \n",
	nF_mm + nt2_rare_mc_mm                + wz_nt2_mm, sqrt(mm_tot_sqerr1 + mm_tot_sqerr2),
	nF_em + nt2_rare_mc_em + nt2_em_chmid + wz_nt2_em, sqrt(em_tot_sqerr1 + em_tot_sqerr2),
	nF_ee + nt2_rare_mc_ee + nt2_ee_chmid + wz_nt2_ee, sqrt(ee_tot_sqerr1 + ee_tot_sqerr2),
	tot_pred, sqrt(comb_tot_sqerr1 + comb_tot_sqerr2));
	fOUTSTREAM << Form("\\bf{Observed}       & \\multicolumn{2}{c|}{\\bf{%3.0f}} & \\multicolumn{2}{c|}{\\bf{%3.0f}}  & \\multicolumn{2}{c|}{\\bf{%3.0f}}  & \\multicolumn{2}{c}{\\bf{%3.0f}}  \\\\ \n",
	nt2_mm, nt2_em, nt2_ee, nt2_mm+nt2_em+nt2_ee);
	fOUTSTREAM << Form("ttW Prod.      & %5.1f &$\\pm$ %5.1f & %5.1f &$\\pm$ %5.1f & %5.1f &$\\pm$ %5.1f & %5.1f &$\\pm$ %5.1f \\\\ \\hline \n",
	ttw_nt2_mm, sqrt(ttw_nt2_mm_e1 + RareESyst2*ttw_nt2_mm*ttw_nt2_mm),
	ttw_nt2_em, sqrt(ttw_nt2_em_e1 + RareESyst2*ttw_nt2_em*ttw_nt2_em),
	ttw_nt2_ee, sqrt(ttw_nt2_ee_e1 + RareESyst2*ttw_nt2_ee*ttw_nt2_ee),
	ttw_nt2_ee + ttw_nt2_mm + ttw_nt2_em, sqrt(ttw_nt2_mm_e1 + ttw_nt2_ee_e1 + ttw_nt2_em_e1 + RareESyst2*(ttw_nt2_ee + ttw_nt2_mm + ttw_nt2_em)*(ttw_nt2_ee + ttw_nt2_mm + ttw_nt2_em)));
	fOUTSTREAM << Form("ttZ Prod.      & %5.1f &$\\pm$ %5.1f & %5.1f &$\\pm$ %5.1f & %5.1f &$\\pm$ %5.1f & %5.1f &$\\pm$ %5.1f \\\\ \\hline \n",
	ttz_nt2_mm, sqrt(ttz_nt2_mm_e1 + RareESyst2*ttz_nt2_mm*ttz_nt2_mm),
	ttz_nt2_em, sqrt(ttz_nt2_em_e1 + RareESyst2*ttz_nt2_em*ttz_nt2_em),
	ttz_nt2_ee, sqrt(ttz_nt2_ee_e1 + RareESyst2*ttz_nt2_ee*ttz_nt2_ee),
	ttz_nt2_ee + ttz_nt2_mm + ttz_nt2_em, sqrt(ttz_nt2_mm_e1 + ttz_nt2_ee_e1 + ttz_nt2_em_e1 + RareESyst2*(ttz_nt2_ee + ttz_nt2_mm + ttz_nt2_em)*(ttz_nt2_ee + ttz_nt2_mm + ttz_nt2_em)));
	if (separateTTH) {
		fOUTSTREAM << Form("ttH Prod.      & %5.1f &$\\pm$ %5.1f & %5.1f &$\\pm$ %5.1f & %5.1f &$\\pm$ %5.1f & %5.1f &$\\pm$ %5.1f \\\\ \\hline \n",
						   tth_nt2_mm, sqrt(tth_nt2_mm_e1 + RareESyst2*tth_nt2_mm*tth_nt2_mm),
						   tth_nt2_em, sqrt(tth_nt2_em_e1 + RareESyst2*tth_nt2_em*tth_nt2_em),
						   tth_nt2_ee, sqrt(tth_nt2_ee_e1 + RareESyst2*tth_nt2_ee*tth_nt2_ee),
						   tth_nt2_ee + tth_nt2_mm + tth_nt2_em, sqrt(tth_nt2_mm_e1 + tth_nt2_ee_e1 + tth_nt2_em_e1 + RareESyst2*(tth_nt2_ee + tth_nt2_mm + tth_nt2_em)*(tth_nt2_ee + tth_nt2_mm + tth_nt2_em)));
	}
	fOUTSTREAM << endl;
	
	
	///////////////////////////////////////////////////////////////////////////////////
	//  OUTPUT FOR DATACARD  //////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////	
	TTWZPrediction pred;
	pred.obs      = nt2_mm+nt2_em+nt2_ee;
	pred.obs_mm   = nt2_mm;
	pred.obs_ee   = nt2_ee;
	pred.obs_em   = nt2_em;
	
	pred.ttw      = ttw_nt2_ee + ttw_nt2_mm + ttw_nt2_em;
	pred.ttw_mm   = ttw_nt2_mm;
	pred.ttw_ee   = ttw_nt2_ee;
	pred.ttw_em   = ttw_nt2_em;
	
	pred.ttz      = ttz_nt2_ee + ttz_nt2_mm + ttz_nt2_em;
	pred.ttz_mm   = ttz_nt2_mm;
	pred.ttz_ee   = ttz_nt2_ee;
	pred.ttz_em   = ttz_nt2_em;
	
	pred.ttwz     = pred.ttw + pred.ttz;
	pred.ttwz_mm  = ttw_nt2_mm + ttz_nt2_mm;
	pred.ttwz_ee  = ttw_nt2_ee + ttz_nt2_ee;
	pred.ttwz_em  = ttw_nt2_em + ttz_nt2_em;
	
	pred.fake     = nF;
	pred.fake_mm  = nF_mm;
	pred.fake_ee  = nF_ee;
	pred.fake_em  = nF_em;
	
	pred.fake_err    = sqrt(FR->getTotEStat()*FR->getTotEStat() + FakeESyst2*nF*nF);
	pred.fake_err_mm = sqrt(FR->getMMTotEStat()*FR->getMMTotEStat() + FakeESyst2*nF_mm*nF_mm);
	pred.fake_err_ee = sqrt(FR->getEETotEStat()*FR->getEETotEStat() + FakeESyst2*nF_ee*nF_ee);
	pred.fake_err_em = sqrt(FR->getEMTotEStat()*FR->getEMTotEStat() + FakeESyst2*nF_em*nF_em);
	
	pred.cmid        = nt2_ee_chmid + nt2_em_chmid;
	pred.cmid_ee     = nt2_ee_chmid;
	pred.cmid_em     = nt2_em_chmid;
	pred.cmid_err    = sqrt(nt2_ee_chmid_e1*nt2_ee_chmid_e1 + nt2_ee_chmid_e2*nt2_ee_chmid_e2 + nt2_em_chmid_e1*nt2_em_chmid_e1 + nt2_em_chmid_e2*nt2_em_chmid_e2);
	pred.cmid_err_ee = sqrt(nt2_ee_chmid_e1*nt2_ee_chmid_e1 + nt2_ee_chmid_e2*nt2_ee_chmid_e2);
	pred.cmid_err_em = sqrt(nt2_em_chmid_e1*nt2_em_chmid_e1 + nt2_em_chmid_e2*nt2_em_chmid_e2);
	
	pred.wz        = wz_nt2_ee + wz_nt2_mm + wz_nt2_em;
	pred.wz_mm     = wz_nt2_mm;
	pred.wz_ee     = wz_nt2_ee;
	pred.wz_em     = wz_nt2_em;
	pred.wz_err    = sqrt(wz_nt2_mm_e1 + wz_nt2_ee_e1 + wz_nt2_em_e1 + WZESyst2*(wz_nt2_ee + wz_nt2_mm + wz_nt2_em)*(wz_nt2_ee + wz_nt2_mm + wz_nt2_em));
	pred.wz_err_mm = sqrt(wz_nt2_mm_e1 + WZESyst2*wz_nt2_mm*wz_nt2_mm);
	pred.wz_err_ee = sqrt(wz_nt2_ee_e1 + WZESyst2*wz_nt2_ee*wz_nt2_ee);
	pred.wz_err_em = sqrt(wz_nt2_em_e1 + WZESyst2*wz_nt2_em*wz_nt2_em);
	
	pred.rare        = nt2_rare_mc_ee + nt2_rare_mc_mm + nt2_rare_mc_em;
	pred.rare_mm     = nt2_rare_mc_mm;
	pred.rare_ee     = nt2_rare_mc_ee;
	pred.rare_em     = nt2_rare_mc_em;
	pred.rare_err    = sqrt(nt2_rare_mc_ee_e1 + nt2_rare_mc_mm_e1 + nt2_rare_mc_em_e1 + RareESyst2*(nt2_rare_mc_ee + nt2_rare_mc_mm + nt2_rare_mc_em)*(nt2_rare_mc_ee + nt2_rare_mc_mm + nt2_rare_mc_em));
	pred.rare_err_mm = sqrt(nt2_rare_mc_mm_e1 + RareESyst2*nt2_rare_mc_mm*nt2_rare_mc_mm);
	pred.rare_err_ee = sqrt(nt2_rare_mc_ee_e1 + RareESyst2*nt2_rare_mc_ee*nt2_rare_mc_ee);
	pred.rare_err_em = sqrt(nt2_rare_mc_em_e1 + RareESyst2*nt2_rare_mc_em*nt2_rare_mc_em);

	///////////////////////////////////////////////////////////////////////////////////
	//  OUTPUT AS PLOT  ///////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////
	TH1D    *h_obs        = new TH1D("h_observed",      "Observed number of events",  4, 0., 4.);
	TH1D    *h_pred_sfake = new TH1D("h_pred_sfake",    "Predicted fakes",            4, 0., 4.);
	TH1D    *h_pred_chmid = new TH1D("h_pred_chmid",    "Predicted charge mis id",    4, 0., 4.);
	TH1D    *h_pred_wz    = new TH1D("h_pred_wz",       "Predicted WZ",               4, 0., 4.);
	TH1D    *h_pred_mc    = new TH1D("h_pred_mc",       "Predicted Rare SM",          4, 0., 4.);
	TH1D    *h_pred_ttw   = new TH1D("h_pred_ttw",      "Predicted ttW",              4, 0., 4.);
	TH1D    *h_pred_ttz   = new TH1D("h_pred_ttz",      "Predicted ttZ",              4, 0., 4.);
	TH1D    *h_pred_tth   = new TH1D("h_pred_tth",      "Predicted tth",              4, 0., 4.);
	TH1D    *h_pred_tot   = new TH1D("h_pred_tot",      "Total Prediction",           4, 0., 4.);
	THStack *hs_pred      = new THStack("hs_predicted", "Predicted number of events");
	
	h_obs->SetMarkerColor(kBlack);
	h_obs->SetMarkerStyle(20);
	h_obs->SetMarkerSize(2);
	h_obs->SetLineWidth(2);
	h_obs->SetLineColor(kBlack);
	h_obs->SetFillColor(kBlack);
	
	h_pred_sfake->SetLineWidth(1);
	h_pred_chmid->SetLineWidth(1);
	h_pred_mc   ->SetLineWidth(1);
	h_pred_wz   ->SetLineWidth(1);
	h_pred_ttw  ->SetLineWidth(1);
	h_pred_ttz  ->SetLineWidth(1);
	h_pred_tth  ->SetLineWidth(1);
	h_pred_sfake->SetFillColor(46);
	h_pred_chmid->SetFillColor(49);
	h_pred_wz   ->SetFillColor(39);
	h_pred_mc   ->SetFillColor(38);
	h_pred_ttw  ->SetFillColor(44);
	h_pred_ttz  ->SetFillColor(42);
	h_pred_tth  ->SetFillColor(45);
	//h_pred_sfake->SetLineColor(46);
	//h_pred_chmid->SetLineColor(42);
	//h_pred_mc   ->SetLineColor(44);
	//h_pred_ttw  ->SetLineColor(29);
	//h_pred_ttz  ->SetLineColor(30);
	//h_pred_tth  ->SetLineColor(31);

	h_pred_tot  ->SetLineWidth(1);
	// h_pred_tot  ->SetFillColor(kBlack);
	// h_pred_tot  ->SetFillStyle(3013);
	h_pred_tot  ->SetFillColor(12);
	h_pred_tot  ->SetFillStyle(3005);
	
	// Add numbers:
	h_obs->SetBinContent(1, nt2_ee);
	h_obs->SetBinContent(2, nt2_mm);
	h_obs->SetBinContent(3, nt2_em);
	h_obs->SetBinContent(4, nt2_ee+nt2_mm+nt2_em);
	
	TGraphAsymmErrors* gr_obs = FR->getGraphPoissonErrors( h_obs );
	gr_obs->SetMarkerColor(kBlack);
	gr_obs->SetMarkerStyle(20);
	gr_obs->SetMarkerSize(2.0);
	gr_obs->SetLineWidth(2);
	gr_obs->SetLineColor(kBlack);
	gr_obs->SetFillColor(kBlack);
	
	h_pred_sfake->SetBinContent(1,nff_ee + npf_ee);
	h_pred_sfake->SetBinContent(2,nff_mm + npf_mm);
	h_pred_sfake->SetBinContent(3,nff_em + npf_em+nfp_em);
	h_pred_sfake->SetBinContent(4,nff_ee+nff_mm+nff_em + npf_ee+npf_mm+npf_em+nfp_em);
	h_pred_sfake->GetXaxis()->SetBinLabel(1, "ee");
	h_pred_sfake->GetXaxis()->SetBinLabel(2, "#mu#mu");
	h_pred_sfake->GetXaxis()->SetBinLabel(3, "e#mu");
	h_pred_sfake->GetXaxis()->SetBinLabel(4, "Total");
	
	h_pred_chmid->SetBinContent(1, nt2_ee_chmid);
	h_pred_chmid->SetBinContent(2, 0.);
	h_pred_chmid->SetBinContent(3, nt2_em_chmid);
	h_pred_chmid->SetBinContent(4, nt2_ee_chmid+nt2_em_chmid);

	h_pred_wz   ->SetBinContent(1, wz_nt2_ee);
	h_pred_wz   ->SetBinContent(2, wz_nt2_mm);
	h_pred_wz   ->SetBinContent(3, wz_nt2_em);
	h_pred_wz   ->SetBinContent(4, wz_nt2_ee + wz_nt2_mm + wz_nt2_em);
	
	h_pred_mc   ->SetBinContent(1, nt2_rare_mc_ee);
	h_pred_mc   ->SetBinContent(2, nt2_rare_mc_mm);
	h_pred_mc   ->SetBinContent(3, nt2_rare_mc_em);
	h_pred_mc   ->SetBinContent(4, nt2_rare_mc_ee + nt2_rare_mc_mm + nt2_rare_mc_em);
	
	h_pred_ttw  ->SetBinContent(1, ttw_nt2_ee);
	h_pred_ttw  ->SetBinContent(2, ttw_nt2_mm);
	h_pred_ttw  ->SetBinContent(3, ttw_nt2_em);
	h_pred_ttw  ->SetBinContent(4, ttw_nt2_ee + ttw_nt2_mm + ttw_nt2_em);
	
	h_pred_ttz  ->SetBinContent(1, ttz_nt2_ee);
	h_pred_ttz  ->SetBinContent(2, ttz_nt2_mm);
	h_pred_ttz  ->SetBinContent(3, ttz_nt2_em);
	h_pred_ttz  ->SetBinContent(4, ttz_nt2_ee + ttz_nt2_mm + ttz_nt2_em);
	
	h_pred_tth->SetBinContent(1, tth_nt2_ee);
	h_pred_tth->SetBinContent(2, tth_nt2_mm);
	h_pred_tth->SetBinContent(3, tth_nt2_em);
	h_pred_tth->SetBinContent(4, tth_nt2_ee + tth_nt2_mm + tth_nt2_em);
	
	h_pred_tot->Add(h_pred_sfake);
	h_pred_tot->Add(h_pred_chmid);
	h_pred_tot->Add(h_pred_wz);
	h_pred_tot->Add(h_pred_mc);
	// h_pred_tot->Add(h_pred_ttw);
	// h_pred_tot->Add(h_pred_ttz);
	h_pred_tot->SetBinError(1, sqrt(ee_tot_sqerr1 + ee_tot_sqerr2));
	h_pred_tot->SetBinError(2, sqrt(mm_tot_sqerr1 + mm_tot_sqerr2));
	h_pred_tot->SetBinError(3, sqrt(em_tot_sqerr1 + em_tot_sqerr2));
	h_pred_tot->SetBinError(4, sqrt(comb_tot_sqerr1 + comb_tot_sqerr2));
	
	hs_pred->Add(h_pred_mc);
	hs_pred->Add(h_pred_wz);
	hs_pred->Add(h_pred_chmid);
	hs_pred->Add(h_pred_sfake);
	hs_pred->Add(h_pred_ttw);
	hs_pred->Add(h_pred_ttz);
	if (separateTTH) hs_pred->Add(h_pred_tth);
	
	float max = 2.0*h_pred_tot->GetBinContent(4);
	
//	if(reg != TTbarWPresel) max = 35.;
	if(gRegions[reg]->sname == "TTbarWPresel") max = 250.;
	if(gRegions[reg]->sname == "TTbarWSel") max = 15.;
	if(gRegions[reg]->sname == "HT200MET50") max = 10.;
	if(gRegions[reg]->sname == "HT200MET50b") max = 2.;
	
	h_obs       ->SetMaximum(max>1?max+1:1.);
	h_pred_sfake->SetMaximum(max>1?max+1:1.);
	h_pred_chmid->SetMaximum(max>1?max+1:1.);
	h_pred_mc   ->SetMaximum(max>1?max+1:1.);
	h_pred_wz   ->SetMaximum(max>1?max+1:1.);
	h_pred_ttw  ->SetMaximum(max>1?max+1:1.);
	h_pred_ttz  ->SetMaximum(max>1?max+1:1.);
	if (separateTTH) h_pred_tth  ->SetMaximum(max>1?max+1:1.);
	h_pred_tot  ->SetMaximum(max>1?max+1:1.);
	hs_pred     ->SetMaximum(max>1?max+1:1.);
	
	hs_pred->Draw("goff");
	hs_pred->GetXaxis()->SetBinLabel(1, "ee");
	hs_pred->GetXaxis()->SetBinLabel(2, "#mu#mu");
	hs_pred->GetXaxis()->SetBinLabel(3, "e#mu");
	hs_pred->GetXaxis()->SetBinLabel(4, "Total");
	hs_pred->GetXaxis()->SetLabelOffset(0.01);
	hs_pred->GetXaxis()->SetLabelFont(42);
	hs_pred->GetXaxis()->SetLabelSize(0.08);
	
	hs_pred->GetYaxis()->SetLabelSize(0.045);
	hs_pred->GetYaxis()->SetTitle("Events");
	hs_pred->GetYaxis()->SetTitleOffset(1.25);
	hs_pred->GetYaxis()->SetTitleSize(0.045);
	
	TLegend *leg = new TLegend(0.16,0.43,0.48,0.85);
	leg->AddEntry(h_obs,        "Data","p");
	leg->AddEntry(h_pred_ttw,   "t#bar{t} + W","f");
	leg->AddEntry(h_pred_ttz,   "t#bar{t} + Z","f");
	if (separateTTH) leg->AddEntry(h_pred_tth,   "t#bar{t} + H","f");
	leg->AddEntry(h_pred_sfake, "Fakes","f");
	leg->AddEntry(h_pred_chmid, "Charge MisID","f");
	leg->AddEntry(h_pred_wz,    "WZ","f");
	leg->AddEntry(h_pred_mc,    "Rare SM","f");
	//leg->AddEntry(h_pred_tot,   "Total Uncertainty","f");
	leg->SetFillStyle(0);
	leg->SetTextFont(42);
	leg->SetTextSize(0.038);
	leg->SetBorderSize(0);
	
	TCanvas *c_temp = new TCanvas("C_ObsPred", "Observed vs Predicted", 0, 0, 600, 600);
	c_temp->SetLeftMargin(0.12);
	c_temp->SetRightMargin(0.04);
	c_temp->cd();
	
	hs_pred->Draw("hist");
	h_pred_tot->DrawCopy("0 E2 same");
	h_pred_tot->SetFillStyle(0);
	h_pred_tot->DrawCopy("hist same");
	gr_obs->Draw("P same");
	leg->Draw();
	
	// drawRegionSel(reg);
	drawTopLine(0.49, 1.0, 0.13);
	
	gPad->RedrawAxis();
	
	// Util::PrintNoEPS(c_temp, "ObsPred_" + Region::sname[reg], fOutputDir + fOutputSubDir, NULL);
	Util::PrintPDF(c_temp,   "ObsPred_" + gRegions[reg]->sname, fOutputDir + fOutputSubDir);
		
	delete c_temp;
	delete h_obs, h_pred_sfake, h_pred_chmid, h_pred_mc, h_pred_wz, h_pred_ttw, h_pred_ttz, h_pred_tth, h_pred_tot, hs_pred;
	delete gr_obs;
	delete FR;
	
	return pred;
}

void SSDLPlotter::makeTTWIntPredictionsSigEvent(float minHT, float maxHT, float minMET, float maxMET, int minNjets, int minNbjetsL, int minNbjetsM, float pT1, float pT2, int chVeto, bool ttw , int flag){
	TString outputdir = Util::MakeOutputDir(fOutputDir + "IntPredictionsSigEventTree");
	fOutputSubDir = "IntPredictionsSigEventTree/";
	time_t rawtime;
	struct tm* timeinfo;
	time(&rawtime);
	timeinfo = localtime(&rawtime);
	// access a chararray containing the date with asctime(timeinfo)
	
	TString chargeString = "";
	if (chVeto == -1) chargeString = "_mimi";
	if (chVeto == +1) chargeString = "_plpl";

	TString tablefilename = outputdir + "Table2.tex";
	TString notetable     = outputdir + "NoteTable.tex";
	fOUTSTREAM.open(tablefilename.Data(), ios::trunc);
	fOUTSTREAM << "==========================================================================================================" << endl;
	fOUTSTREAM << " Table 2 inputs from ETH Analysis" << endl;
	fOUTSTREAM << Form(" Generated on: %s ", asctime(timeinfo)) << endl;
	fOUTSTREAM << endl;
	
	fOUTSTREAM3.open(notetable.Data(), ios::trunc);
	fOUTSTREAM3 << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
	fOUTSTREAM3 << Form("%%%% Generated on: %s ", asctime(timeinfo)) << endl;
	fOUTSTREAM3 << "%% Format is tot, (ee, mm, em)" << endl;
	fOUTSTREAM3 << endl;
	
	// ttwregions.push_back(gRegion["TTbarWPreselpp"]);
	// ttwregions.push_back(gRegion["TTbarWSelpp"]);
	// vector<TTWZPrediction> ttwzpreds;
	std::map< TString, TTWZPrediction > ttwzpreds;
	// for(size_t i = 0; i < ttwregions.size(); ++i){
	gsystIt = gSystematics.begin();
	for( ; gsystIt != gSystematics.end() ; ++gsystIt){
		TString outputname = outputdir + "DataPred_" + gsystIt->first + chargeString + ".txt";
		ttwzpreds[gsystIt->first] = makePredictionSignalEvents(minHT, maxHT, minMET, maxMET, minNjets, minNbjetsL, minNbjetsM, pT1, pT2, chVeto, ttw, gsystIt->second);
	// here we are. fix this stuff
	}
	cout << "------------------------------------------------" << endl;
	cout << "size of gRegion map:        " << gRegion.size() << endl;
	cout << "size of gSystematics map:   " << gSystematics.size() << endl;
	cout << "size of my systematics map: " << ttwzpreds.size() << endl;
	cout << " example: observed in the normal mode:   " << ttwzpreds["Normal"].obs << endl;
	cout << " example: ttwz of the b-scale down mode: " << ttwzpreds["BDown"].ttwz << endl;
	
	fOUTSTREAM.close();
	fOUTSTREAM2.close();
	fOUTSTREAM3.close();
	
	////////////////////////////////////////////////////
	// Datacards and systematics
	TString datacard;
	const float pu_syst = 1.03;
	const float match_syst_up = 1.01485;
	const float match_syst_dn = 0.997558;
	const float scale_syst_up = 1.02302;
	const float scale_syst_dn = 0.965962;
	TString lumiError = "1.045";
	
	datacard = outputdir + "datacard_TTWZ" + chargeString + ".txt";
	fOUTSTREAM.open(datacard.Data(), ios::trunc);
	fOUTSTREAM <<      "#=========================================================================================" << endl;
	fOUTSTREAM <<      "# Systematics table for ttW/Z analysis, same-sign channel" << endl;
	fOUTSTREAM << Form("# Generated on: %s ", asctime(timeinfo)) << endl;
	fOUTSTREAM <<      "# Copy between the dashed lines for datacard" << endl;
	fOUTSTREAM <<      "#-----------------------------------------------------------------------------------------" << endl;
	fOUTSTREAM <<      "imax 1" << endl;
	fOUTSTREAM <<      "jmax 4" << endl;
	fOUTSTREAM <<      "kmax *" << endl;
	fOUTSTREAM << endl << endl;
	fOUTSTREAM <<      "bin\t\t1" << endl;
	if (gFullDataBlind)
		fOUTSTREAM << Form("observation\t%d", 999) << endl;
	else
		fOUTSTREAM << Form("observation\t%d", ttwzpreds["Normal"].obs) << endl;
	fOUTSTREAM << endl << endl;
	fOUTSTREAM <<      "bin\t\t1\t\t1\t\t1\t\t1\t\t1" << endl;
	fOUTSTREAM <<      "process\t\tttWZ\t\tfake\t\tcmid\t\twz\t\trare" << endl;
	fOUTSTREAM <<      "process\t\t0\t\t1\t\t2\t\t3\t\t4" << endl;
	fOUTSTREAM << Form("rate\t\t%5.3f\t\t%5.3f\t\t%5.3f\t\t%5.3f\t\t%5.3f",
					   ttwzpreds["Normal"].ttwz, ttwzpreds["Normal"].fake, ttwzpreds["Normal"].cmid, ttwzpreds["Normal"].wz, ttwzpreds["Normal"].rare) << endl;
	fOUTSTREAM << Form("# separate signal yields: %5.3f (ttW), %5.3f (ttZ)", ttwzpreds["Normal"].ttw, ttwzpreds["Normal"].ttz) << endl;
	fOUTSTREAM << endl << endl;
	fOUTSTREAM <<      "#syst" << endl;
	fOUTSTREAM <<      "lumi     lnN\t"+lumiError+"\t\t"+lumiError+"\t\t"+lumiError+"\t\t"+lumiError+"\t\t"+lumiError+"" << endl;
	fOUTSTREAM << Form("bgUncfak lnN\t-\t\t%5.3f\t\t-\t\t-\t\t-", 1.0+ttwzpreds["Normal"].fake_err/ttwzpreds["Normal"].fake) << endl;
	fOUTSTREAM << Form("bgUnccmi lnN\t-\t\t-\t\t%5.3f\t\t-\t\t-", 1.0+ttwzpreds["Normal"].cmid_err/ttwzpreds["Normal"].cmid) << endl;
	fOUTSTREAM << Form("bgUncwz  lnN\t-\t\t-\t\t-\t\t%5.3f\t\t-", 1.0+ttwzpreds["Normal"].wz_err  /ttwzpreds["Normal"].wz)   << endl;
	fOUTSTREAM << Form("bgUncrar lnN\t-\t\t-\t\t-\t\t-\t\t%5.3f", 1.0+ttwzpreds["Normal"].rare_err/ttwzpreds["Normal"].rare) << endl;
	fOUTSTREAM << Form("lept     lnN\t%5.3f/%5.3f\t-\t\t-\t\t%5.3f/%5.3f\t%5.3f/%5.3f",
					   1.0+(ttwzpreds["LepDown"].ttwz-ttwzpreds["Normal"].ttwz)/ttwzpreds["Normal"].ttwz,
					   1.0+(ttwzpreds["LepUp"]  .ttwz-ttwzpreds["Normal"].ttwz)/ttwzpreds["Normal"].ttwz,
					   1.0+(ttwzpreds["LepDown"].wz  -ttwzpreds["Normal"].wz  )/ttwzpreds["Normal"].wz,
					   1.0+(ttwzpreds["LepUp"]  .wz  -ttwzpreds["Normal"].wz  )/ttwzpreds["Normal"].wz,
					   1.0+(ttwzpreds["LepDown"].rare-ttwzpreds["Normal"].rare)/ttwzpreds["Normal"].rare,
					   1.0+(ttwzpreds["LepUp"]  .rare-ttwzpreds["Normal"].rare)/ttwzpreds["Normal"].rare) << endl;
	fOUTSTREAM << Form("btag     lnN\t%5.3f/%5.3f\t-\t\t-\t\t%5.3f/%5.3f\t%5.3f/%5.3f",
					   1.0+(ttwzpreds["BDown"].ttwz-ttwzpreds["Normal"].ttwz)/ttwzpreds["Normal"].ttwz,
					   1.0+(ttwzpreds["BUp"].ttwz-ttwzpreds["Normal"].ttwz)/ttwzpreds["Normal"].ttwz,
					   1.0+(ttwzpreds["BDown"].wz  -ttwzpreds["Normal"].wz  )/ttwzpreds["Normal"].wz,
					   1.0+(ttwzpreds["BUp"].wz  -ttwzpreds["Normal"].wz  )/ttwzpreds["Normal"].wz,
					   1.0+(ttwzpreds["BDown"].rare-ttwzpreds["Normal"].rare)/ttwzpreds["Normal"].rare,
					   1.0+(ttwzpreds["BUp"].rare-ttwzpreds["Normal"].rare)/ttwzpreds["Normal"].rare) << endl;
	fOUTSTREAM << Form("jes      lnN\t%5.3f/%5.3f\t-\t\t-\t\t%5.3f/%5.3f\t%5.3f/%5.3f",
					   1.0+(ttwzpreds["JetDown"].ttwz-ttwzpreds["Normal"].ttwz)/ttwzpreds["Normal"].ttwz,
					   1.0+(ttwzpreds["JetUp"].ttwz-ttwzpreds["Normal"].ttwz)/ttwzpreds["Normal"].ttwz,
					   1.0+(ttwzpreds["JetDown"].wz  -ttwzpreds["Normal"].wz  )/ttwzpreds["Normal"].wz,
					   1.0+(ttwzpreds["JetUp"].wz  -ttwzpreds["Normal"].wz  )/ttwzpreds["Normal"].wz,
					   1.0+(ttwzpreds["JetDown"].rare-ttwzpreds["Normal"].rare)/ttwzpreds["Normal"].rare,
					   1.0+(ttwzpreds["JetUp"].rare-ttwzpreds["Normal"].rare)/ttwzpreds["Normal"].rare) << endl;
	fOUTSTREAM << Form("jer      lnN\t%5.3f\t\t-\t\t-\t\t%5.3f\t\t%5.3f",
					   1.0+(ttwzpreds["JetSmear"].ttwz-ttwzpreds["Normal"].ttwz)/ttwzpreds["Normal"].ttwz,
					   1.0+(ttwzpreds["JetSmear"].wz  -ttwzpreds["Normal"].wz  )/ttwzpreds["Normal"].wz,
					   1.0+(ttwzpreds["JetSmear"].rare-ttwzpreds["Normal"].rare)/ttwzpreds["Normal"].rare) << endl;
	fOUTSTREAM << Form("pu       lnN\t%5.3f\t\t-\t\t-\t\t%5.3f\t\t%5.3f",
					   pu_syst, pu_syst, pu_syst) << endl;
	fOUTSTREAM << Form("matching lnN\t%5.3f/%5.3f\t-\t\t-\t\t%5.3f/%5.3f\t%5.3f/%5.3f",
					   match_syst_up, match_syst_dn, match_syst_up, match_syst_dn, match_syst_up, match_syst_dn) << endl;
	fOUTSTREAM << Form("scale    lnN\t%5.3f/%5.3f\t-\t\t-\t\t%5.3f/%5.3f\t%5.3f/%5.3f",
					   scale_syst_up, scale_syst_dn, scale_syst_up, scale_syst_dn, scale_syst_up, scale_syst_dn) << endl;
	fOUTSTREAM << endl;
	fOUTSTREAM.close();
	
	datacard = outputdir + "datacard_ssdl_ttW" + chargeString + ".txt";
	fOUTSTREAM.open(datacard.Data(), ios::trunc);
	fOUTSTREAM <<      "#=========================================================================================" << endl;
	fOUTSTREAM <<      "# Systematics table for ttW/Z analysis, same-sign channel, ttW and ttZ separate" << endl;
	fOUTSTREAM << Form("# Generated on: %s ", asctime(timeinfo)) << endl;
	fOUTSTREAM <<      "# Copy between the dashed lines for datacard" << endl;
	fOUTSTREAM <<      "#-----------------------------------------------------------------------------------------" << endl;
	fOUTSTREAM <<      "imax 1" << endl;
	fOUTSTREAM <<      "jmax 5" << endl;
	fOUTSTREAM <<      "kmax *" << endl;
	fOUTSTREAM << endl << endl;
	fOUTSTREAM <<      "bin\t\t1" << endl;
	if (gFullDataBlind)
		fOUTSTREAM << Form("observation\t%d", 999) << endl;
	else
		fOUTSTREAM << Form("observation\t%d", ttwzpreds["Normal"].obs) << endl;
	fOUTSTREAM << endl << endl;
	fOUTSTREAM <<      "bin\t\t1\t\t1\t\t1\t\t1\t\t1\t\t1" << endl;
	fOUTSTREAM <<      "process\t\tttW\t\tttZ\t\tfake\t\tcmid\t\twz\t\trare" << endl;
	fOUTSTREAM <<      "process\t\t0\t\t1\t\t2\t\t3\t\t4\t\t5" << endl;
	fOUTSTREAM << Form("rate\t\t%5.3f\t\t%5.3f\t\t%5.3f\t\t%5.3f\t\t%5.3f\t\t%5.3f",
					   ttwzpreds["Normal"].ttw, ttwzpreds["Normal"].ttz, ttwzpreds["Normal"].fake, ttwzpreds["Normal"].cmid, ttwzpreds["Normal"].wz, ttwzpreds["Normal"].rare) << endl;
	fOUTSTREAM << endl << endl;
	fOUTSTREAM <<      "#syst" << endl;
	fOUTSTREAM <<      "lumi     lnN\t"+lumiError+"\t\t"+lumiError+"\t\t"+lumiError+"\t\t"+lumiError+"\t\t"+lumiError+"\t\t"+lumiError+"" << endl;
	fOUTSTREAM << Form("bgUncttz lnN\t-\t\t%5.3f\t\t-\t\t-\t\t-\t\t-", 1.0+ttwzpreds["Normal"].ttz_err/ttwzpreds["Normal"].ttz)   << endl;
	fOUTSTREAM << Form("bgUncfak lnN\t-\t\t-\t\t%5.3f\t\t-\t\t-\t\t-", 1.0+ttwzpreds["Normal"].fake_err/ttwzpreds["Normal"].fake) << endl;
	fOUTSTREAM << Form("bgUnccmi lnN\t-\t\t-\t\t-\t\t%5.3f\t\t-\t\t-", 1.0+ttwzpreds["Normal"].cmid_err/ttwzpreds["Normal"].cmid) << endl;
	fOUTSTREAM << Form("bgUncwz  lnN\t-\t\t-\t\t-\t\t-\t\t%5.3f\t\t-", 1.0+ttwzpreds["Normal"].wz_err  /ttwzpreds["Normal"].wz)   << endl;
	fOUTSTREAM << Form("bgUncrar lnN\t-\t\t-\t\t-\t\t-\t\t-\t\t%5.3f", 1.0+ttwzpreds["Normal"].rare_err/ttwzpreds["Normal"].rare) << endl;
	fOUTSTREAM << Form("lept     lnN\t%5.3f/%5.3f\t%5.3f/%5.3f\t-\t\t-\t\t%5.3f/%5.3f\t%5.3f/%5.3f",
					   1.0+(ttwzpreds["LepDown"].ttwz-ttwzpreds["Normal"].ttwz)/ttwzpreds["Normal"].ttwz,
					   1.0+(ttwzpreds["LepUp"].ttwz-ttwzpreds["Normal"].ttwz)/ttwzpreds["Normal"].ttwz,
					   1.0+(ttwzpreds["LepDown"].ttwz-ttwzpreds["Normal"].ttwz)/ttwzpreds["Normal"].ttwz,
					   1.0+(ttwzpreds["LepUp"].ttwz-ttwzpreds["Normal"].ttwz)/ttwzpreds["Normal"].ttwz,
					   1.0+(ttwzpreds["LepDown"].wz  -ttwzpreds["Normal"].wz  )/ttwzpreds["Normal"].wz,
					   1.0+(ttwzpreds["LepUp"].wz  -ttwzpreds["Normal"].wz  )/ttwzpreds["Normal"].wz,
					   1.0+(ttwzpreds["LepDown"].rare-ttwzpreds["Normal"].rare)/ttwzpreds["Normal"].rare,
					   1.0+(ttwzpreds["LepUp"].rare-ttwzpreds["Normal"].rare)/ttwzpreds["Normal"].rare) << endl;
	fOUTSTREAM << Form("btag     lnN\t%5.3f/%5.3f\t%5.3f/%5.3f\t-\t\t-\t\t%5.3f/%5.3f\t%5.3f/%5.3f",
					   1.0+(ttwzpreds["BDown"].ttwz-ttwzpreds["Normal"].ttwz)/ttwzpreds["Normal"].ttwz,
					   1.0+(ttwzpreds["BUp"].ttwz-ttwzpreds["Normal"].ttwz)/ttwzpreds["Normal"].ttwz,
					   1.0+(ttwzpreds["BDown"].ttwz-ttwzpreds["Normal"].ttwz)/ttwzpreds["Normal"].ttwz,
					   1.0+(ttwzpreds["BUp"].ttwz-ttwzpreds["Normal"].ttwz)/ttwzpreds["Normal"].ttwz,
					   1.0+(ttwzpreds["BDown"].wz  -ttwzpreds["Normal"].wz  )/ttwzpreds["Normal"].wz,
					   1.0+(ttwzpreds["BUp"].wz  -ttwzpreds["Normal"].wz  )/ttwzpreds["Normal"].wz,
					   1.0+(ttwzpreds["BDown"].rare-ttwzpreds["Normal"].rare)/ttwzpreds["Normal"].rare,
					   1.0+(ttwzpreds["BUp"].rare-ttwzpreds["Normal"].rare)/ttwzpreds["Normal"].rare) << endl;
	fOUTSTREAM << Form("jes      lnN\t%5.3f/%5.3f\t%5.3f/%5.3f\t-\t\t-\t\t%5.3f/%5.3f\t%5.3f/%5.3f",
					   1.0+(ttwzpreds["JetDown"].ttwz-ttwzpreds["Normal"].ttwz)/ttwzpreds["Normal"].ttwz,
					   1.0+(ttwzpreds["JetUp"].ttwz-ttwzpreds["Normal"].ttwz)/ttwzpreds["Normal"].ttwz,
					   1.0+(ttwzpreds["JetDown"].ttwz-ttwzpreds["Normal"].ttwz)/ttwzpreds["Normal"].ttwz,
					   1.0+(ttwzpreds["JetUp"].ttwz-ttwzpreds["Normal"].ttwz)/ttwzpreds["Normal"].ttwz,
					   1.0+(ttwzpreds["JetDown"].wz  -ttwzpreds["Normal"].wz  )/ttwzpreds["Normal"].wz,
					   1.0+(ttwzpreds["JetUp"].wz  -ttwzpreds["Normal"].wz  )/ttwzpreds["Normal"].wz,
					   1.0+(ttwzpreds["JetDown"].rare-ttwzpreds["Normal"].rare)/ttwzpreds["Normal"].rare,
					   1.0+(ttwzpreds["JetUp"].rare-ttwzpreds["Normal"].rare)/ttwzpreds["Normal"].rare) << endl;
	fOUTSTREAM << Form("jer      lnN\t%5.3f\t\t%5.3f\t\t-\t\t-\t\t%5.3f\t\t%5.3f",
					   1.0+(ttwzpreds["JetSmear"].ttwz-ttwzpreds["Normal"].ttwz)/ttwzpreds["Normal"].ttwz,
					   1.0+(ttwzpreds["JetSmear"].ttwz-ttwzpreds["Normal"].ttwz)/ttwzpreds["Normal"].ttwz,
					   1.0+(ttwzpreds["JetSmear"].wz  -ttwzpreds["Normal"].wz  )/ttwzpreds["Normal"].wz,
					   1.0+(ttwzpreds["JetSmear"].rare-ttwzpreds["Normal"].rare)/ttwzpreds["Normal"].rare) << endl;
	fOUTSTREAM << Form("pu       lnN\t%5.3f\t\t%5.3f\t\t-\t\t-\t\t%5.3f\t\t%5.3f",
					   pu_syst, pu_syst, pu_syst, pu_syst) << endl;
	fOUTSTREAM << Form("matching lnN\t%5.3f/%5.3f\t%5.3f/%5.3f\t-\t\t-\t\t%5.3f/%5.3f\t%5.3f/%5.3f",
					   match_syst_up, match_syst_dn, match_syst_up, match_syst_dn, match_syst_up, match_syst_dn, match_syst_up, match_syst_dn) << endl;
	fOUTSTREAM << Form("scale    lnN\t%5.3f/%5.3f\t%5.3f/%5.3f\t-\t\t-\t\t%5.3f/%5.3f\t%5.3f/%5.3f",
					   scale_syst_up, scale_syst_dn, scale_syst_up, scale_syst_dn, scale_syst_up, scale_syst_dn, scale_syst_up, scale_syst_dn) << endl;
	fOUTSTREAM << Form("NLO      lnN\t%5.3f\t\t-\t\t-\t\t-\t\t-\t\t-",
					   1.130) << endl; // old numbers for the moment
	fOUTSTREAM << endl;
	fOUTSTREAM.close();
	
	datacard = outputdir + "datacard_ssdl_mm" + chargeString + ".txt";
	fOUTSTREAM.open(datacard.Data(), ios::trunc);
	fOUTSTREAM <<      "#=========================================================================================" << endl;
	fOUTSTREAM <<      "# Systematics table for ttW/Z analysis, same-sign channel, mumu only" << endl;
	fOUTSTREAM << Form("# Generated on: %s ", asctime(timeinfo)) << endl;
	fOUTSTREAM <<      "# Copy between the dashed lines for datacard" << endl;
	fOUTSTREAM <<      "#-----------------------------------------------------------------------------------------" << endl;
	fOUTSTREAM <<      "imax 1" << endl;
	fOUTSTREAM <<      "jmax 4" << endl;
	fOUTSTREAM <<      "kmax *" << endl;
	fOUTSTREAM << endl << endl;
	fOUTSTREAM <<      "bin\t\t1" << endl;
	if (gFullDataBlind)
		fOUTSTREAM << Form("observation\t%d", 999) << endl;
	else
		fOUTSTREAM << Form("observation\t%d", ttwzpreds["Normal"].obs_mm) << endl;
	fOUTSTREAM << endl << endl;
	fOUTSTREAM <<      "bin\t\t1\t\t1\t\t1\t\t1\t\t1" << endl;
	fOUTSTREAM <<      "process\t\tttWZ\t\tfake\t\tcmid\t\twz\t\trare" << endl;
	fOUTSTREAM <<      "process\t\t0\t\t1\t\t2\t\t3\t\t4" << endl;
	fOUTSTREAM << Form("rate\t\t%5.3f\t\t%5.3f\t\t%5.3f\t\t%5.3f\t\t%5.3f",
					   ttwzpreds["Normal"].ttwz_mm, ttwzpreds["Normal"].fake_mm, 0.0, ttwzpreds["Normal"].wz_mm, ttwzpreds["Normal"].rare_mm) << endl;
	fOUTSTREAM << Form("# separate signal yields: %5.3f (ttW), %5.3f (ttZ)", ttwzpreds["Normal"].ttw_mm, ttwzpreds["Normal"].ttz_mm) << endl;
	fOUTSTREAM << endl << endl;
	fOUTSTREAM <<      "#syst" << endl;
	fOUTSTREAM <<      "lumi     lnN\t"+lumiError+"\t\t"+lumiError+"\t\t"+lumiError+"\t\t"+lumiError+"\t\t"+lumiError+"" << endl;
	fOUTSTREAM << Form("bgUncfak lnN\t-\t\t%5.3f\t\t-\t\t-\t\t-", 1.0+ttwzpreds["Normal"].fake_err_mm/ttwzpreds["Normal"].fake_mm) << endl;
	fOUTSTREAM <<      "bgUnccmi lnN\t-\t\t-\t\t-\t\t-\t\t-" << endl;
	fOUTSTREAM << Form("bgUncwz  lnN\t-\t\t-\t\t-\t\t%5.3f\t\t-", 1.0+ttwzpreds["Normal"].wz_err_mm  /ttwzpreds["Normal"].wz_mm)   << endl;
	fOUTSTREAM << Form("bgUncrar lnN\t-\t\t-\t\t-\t\t-\t\t%5.3f", 1.0+ttwzpreds["Normal"].rare_err_mm/ttwzpreds["Normal"].rare_mm) << endl;
	fOUTSTREAM << Form("lept     lnN\t%5.3f/%5.3f\t-\t\t-\t\t%5.3f/%5.3f\t%5.3f/%5.3f",
					   1.0+(ttwzpreds["LepDown"].ttwz_mm-ttwzpreds["Normal"].ttwz_mm)/ttwzpreds["Normal"].ttwz_mm,
					   1.0+(ttwzpreds["LepUp"].ttwz_mm-ttwzpreds["Normal"].ttwz_mm)/ttwzpreds["Normal"].ttwz_mm,
					   1.0+(ttwzpreds["LepDown"].wz_mm  -ttwzpreds["Normal"].wz_mm  )/ttwzpreds["Normal"].wz_mm,
					   1.0+(ttwzpreds["LepUp"].wz_mm  -ttwzpreds["Normal"].wz_mm  )/ttwzpreds["Normal"].wz_mm,
					   1.0+(ttwzpreds["LepDown"].rare_mm-ttwzpreds["Normal"].rare_mm)/ttwzpreds["Normal"].rare_mm,
					   1.0+(ttwzpreds["LepUp"].rare_mm-ttwzpreds["Normal"].rare_mm)/ttwzpreds["Normal"].rare_mm) << endl;
	fOUTSTREAM << Form("btag     lnN\t%5.3f/%5.3f\t-\t\t-\t\t%5.3f/%5.3f\t%5.3f/%5.3f",
					   1.0+(ttwzpreds["BDown"].ttwz_mm-ttwzpreds["Normal"].ttwz_mm)/ttwzpreds["Normal"].ttwz_mm,
					   1.0+(ttwzpreds["BUp"].ttwz_mm-ttwzpreds["Normal"].ttwz_mm)/ttwzpreds["Normal"].ttwz_mm,
					   1.0+(ttwzpreds["BDown"].wz_mm  -ttwzpreds["Normal"].wz_mm  )/ttwzpreds["Normal"].wz_mm,
					   1.0+(ttwzpreds["BUp"].wz_mm  -ttwzpreds["Normal"].wz_mm  )/ttwzpreds["Normal"].wz_mm,
					   1.0+(ttwzpreds["BDown"].rare_mm-ttwzpreds["Normal"].rare_mm)/ttwzpreds["Normal"].rare_mm,
					   1.0+(ttwzpreds["BUp"].rare_mm-ttwzpreds["Normal"].rare_mm)/ttwzpreds["Normal"].rare_mm) << endl;
	fOUTSTREAM << Form("jes      lnN\t%5.3f/%5.3f\t-\t\t-\t\t%5.3f/%5.3f\t%5.3f/%5.3f",
					   1.0+(ttwzpreds["JetDown"].ttwz_mm-ttwzpreds["Normal"].ttwz_mm)/ttwzpreds["Normal"].ttwz_mm,
					   1.0+(ttwzpreds["JetUp"].ttwz_mm-ttwzpreds["Normal"].ttwz_mm)/ttwzpreds["Normal"].ttwz_mm,
					   1.0+(ttwzpreds["JetDown"].wz_mm  -ttwzpreds["Normal"].wz_mm  )/ttwzpreds["Normal"].wz_mm,
					   1.0+(ttwzpreds["JetUp"].wz_mm  -ttwzpreds["Normal"].wz_mm  )/ttwzpreds["Normal"].wz_mm,
					   1.0+(ttwzpreds["JetDown"].rare_mm-ttwzpreds["Normal"].rare_mm)/ttwzpreds["Normal"].rare_mm,
					   1.0+(ttwzpreds["JetUp"].rare_mm-ttwzpreds["Normal"].rare_mm)/ttwzpreds["Normal"].rare_mm) << endl;
	fOUTSTREAM << Form("jer      lnN\t%5.3f\t\t-\t\t-\t\t%5.3f\t\t%5.3f",
					   1.0+(ttwzpreds["JetSmear"].ttwz_mm-ttwzpreds["Normal"].ttwz_mm)/ttwzpreds["Normal"].ttwz_mm,
					   1.0+(ttwzpreds["JetSmear"].wz_mm  -ttwzpreds["Normal"].wz_mm  )/ttwzpreds["Normal"].wz_mm,
					   1.0+(ttwzpreds["JetSmear"].rare_mm-ttwzpreds["Normal"].rare_mm)/ttwzpreds["Normal"].rare_mm) << endl;
	fOUTSTREAM << Form("pu       lnN\t%5.3f\t\t-\t\t-\t\t%5.3f\t\t%5.3f",
					   pu_syst, pu_syst, pu_syst) << endl;
	fOUTSTREAM << Form("matching lnN\t%5.3f/%5.3f\t-\t\t-\t\t%5.3f/%5.3f\t%5.3f/%5.3f",
					   match_syst_up, match_syst_dn, match_syst_up, match_syst_dn, match_syst_up, match_syst_dn) << endl;
	fOUTSTREAM << Form("scale    lnN\t%5.3f/%5.3f\t-\t\t-\t\t%5.3f/%5.3f\t%5.3f/%5.3f",
					   scale_syst_up, scale_syst_dn, scale_syst_up, scale_syst_dn, scale_syst_up, scale_syst_dn) << endl;
	fOUTSTREAM << endl;
	fOUTSTREAM.close();
	
	datacard = outputdir + "datacard_ssdl_em" + chargeString + ".txt";
	fOUTSTREAM.open(datacard.Data(), ios::trunc);
	fOUTSTREAM <<      "#=========================================================================================" << endl;
	fOUTSTREAM <<      "# Systematics table for ttW/Z analysis, same-sign channel, emu only" << endl;
	fOUTSTREAM << Form("# Generated on: %s ", asctime(timeinfo)) << endl;
	fOUTSTREAM <<      "# Copy between the dashed lines for datacard" << endl;
	fOUTSTREAM <<      "#-----------------------------------------------------------------------------------------" << endl;
	fOUTSTREAM <<      "imax 1" << endl;
	fOUTSTREAM <<      "jmax 4" << endl;
	fOUTSTREAM <<      "kmax *" << endl;
	fOUTSTREAM << endl << endl;
	fOUTSTREAM <<      "bin\t\t1" << endl;
	if (gFullDataBlind)
		fOUTSTREAM << Form("observation\t%d", 999) << endl;
	else
		fOUTSTREAM << Form("observation\t%d", ttwzpreds["Normal"].obs_em) << endl;
	fOUTSTREAM << endl << endl;
	fOUTSTREAM <<      "bin\t\t1\t\t1\t\t1\t\t1\t\t1" << endl;
	fOUTSTREAM <<      "process\t\tttWZ\t\tfake\t\tcmid\t\twz\t\trare" << endl;
	fOUTSTREAM <<      "process\t\t0\t\t1\t\t2\t\t3\t\t4" << endl;
	fOUTSTREAM << Form("rate\t\t%5.3f\t\t%5.3f\t\t%5.3f\t\t%5.3f\t\t%5.3f",
					   ttwzpreds["Normal"].ttwz_em, ttwzpreds["Normal"].fake_em, ttwzpreds["Normal"].cmid_em, ttwzpreds["Normal"].wz_em, ttwzpreds["Normal"].rare_em) << endl;
	fOUTSTREAM << Form("# separate signal yields: %5.3f (ttW), %5.3f (ttZ)", ttwzpreds["Normal"].ttw_em, ttwzpreds["Normal"].ttz_em) << endl;
	fOUTSTREAM << endl << endl;
	fOUTSTREAM <<      "#syst" << endl;
	fOUTSTREAM <<      "lumi     lnN\t"+lumiError+"\t\t"+lumiError+"\t\t"+lumiError+"\t\t"+lumiError+"\t\t"+lumiError+"" << endl;
	fOUTSTREAM << Form("bgUncfak lnN\t-\t\t%5.3f\t\t-\t\t-\t\t-", 1.0+ttwzpreds["Normal"].fake_err_em/ttwzpreds["Normal"].fake_em) << endl;
	fOUTSTREAM << Form("bgUnccmi lnN\t-\t\t-\t\t%5.3f\t\t-\t\t-", 1.0+ttwzpreds["Normal"].cmid_err_em/ttwzpreds["Normal"].cmid_em) << endl;
	fOUTSTREAM << Form("bgUncwz  lnN\t-\t\t-\t\t-\t\t%5.3f\t\t-", 1.0+ttwzpreds["Normal"].wz_err_em  /ttwzpreds["Normal"].wz_em)   << endl;
	fOUTSTREAM << Form("bgUncrar lnN\t-\t\t-\t\t-\t\t-\t\t%5.3f", 1.0+ttwzpreds["Normal"].rare_err_em/ttwzpreds["Normal"].rare_em) << endl;
	fOUTSTREAM << Form("lept     lnN\t%5.3f/%5.3f\t-\t\t-\t\t%5.3f/%5.3f\t%5.3f/%5.3f",
					   1.0+(ttwzpreds["LepDown"].ttwz_em-ttwzpreds["Normal"].ttwz_em)/ttwzpreds["Normal"].ttwz_em,
					   1.0+(ttwzpreds["LepUp"].ttwz_em-ttwzpreds["Normal"].ttwz_em)/ttwzpreds["Normal"].ttwz_em,
					   1.0+(ttwzpreds["LepDown"].wz_em  -ttwzpreds["Normal"].wz_em  )/ttwzpreds["Normal"].wz_em,
					   1.0+(ttwzpreds["LepUp"].wz_em  -ttwzpreds["Normal"].wz_em  )/ttwzpreds["Normal"].wz_em,
					   1.0+(ttwzpreds["LepDown"].rare_em-ttwzpreds["Normal"].rare_em)/ttwzpreds["Normal"].rare_em,
					   1.0+(ttwzpreds["LepUp"].rare_em-ttwzpreds["Normal"].rare_em)/ttwzpreds["Normal"].rare_em) << endl;
	fOUTSTREAM << Form("btag     lnN\t%5.3f/%5.3f\t-\t\t-\t\t%5.3f/%5.3f\t%5.3f/%5.3f",
					   1.0+(ttwzpreds["BDown"].ttwz_em-ttwzpreds["Normal"].ttwz_em)/ttwzpreds["Normal"].ttwz_em,
					   1.0+(ttwzpreds["BUp"].ttwz_em-ttwzpreds["Normal"].ttwz_em)/ttwzpreds["Normal"].ttwz_em,
					   1.0+(ttwzpreds["BDown"].wz_em  -ttwzpreds["Normal"].wz_em  )/ttwzpreds["Normal"].wz_em,
					   1.0+(ttwzpreds["BUp"].wz_em  -ttwzpreds["Normal"].wz_em  )/ttwzpreds["Normal"].wz_em,
					   1.0+(ttwzpreds["BDown"].rare_em-ttwzpreds["Normal"].rare_em)/ttwzpreds["Normal"].rare_em,
					   1.0+(ttwzpreds["BUp"].rare_em-ttwzpreds["Normal"].rare_em)/ttwzpreds["Normal"].rare_em) << endl;
	fOUTSTREAM << Form("jes      lnN\t%5.3f/%5.3f\t-\t\t-\t\t%5.3f/%5.3f\t%5.3f/%5.3f",
					   1.0+(ttwzpreds["JetDown"].ttwz_em-ttwzpreds["Normal"].ttwz_em)/ttwzpreds["Normal"].ttwz_em,
					   1.0+(ttwzpreds["JetUp"].ttwz_em-ttwzpreds["Normal"].ttwz_em)/ttwzpreds["Normal"].ttwz_em,
					   1.0+(ttwzpreds["JetDown"].wz_em  -ttwzpreds["Normal"].wz_em  )/ttwzpreds["Normal"].wz_em,
					   1.0+(ttwzpreds["JetUp"].wz_em  -ttwzpreds["Normal"].wz_em  )/ttwzpreds["Normal"].wz_em,
					   1.0+(ttwzpreds["JetDown"].rare_em-ttwzpreds["Normal"].rare_em)/ttwzpreds["Normal"].rare_em,
					   1.0+(ttwzpreds["JetUp"].rare_em-ttwzpreds["Normal"].rare_em)/ttwzpreds["Normal"].rare_em) << endl;
	fOUTSTREAM << Form("jer      lnN\t%5.3f\t\t-\t\t-\t\t%5.3f\t\t%5.3f",
					   1.0+(ttwzpreds["JetSmear"].ttwz_em-ttwzpreds["Normal"].ttwz_em)/ttwzpreds["Normal"].ttwz_em,
					   1.0+(ttwzpreds["JetSmear"].wz_em  -ttwzpreds["Normal"].wz_em  )/ttwzpreds["Normal"].wz_em,
					   1.0+(ttwzpreds["JetSmear"].rare_em-ttwzpreds["Normal"].rare_em)/ttwzpreds["Normal"].rare_em) << endl;
	fOUTSTREAM << Form("pu       lnN\t%5.3f\t\t-\t\t-\t\t%5.3f\t\t%5.3f",
					   pu_syst, pu_syst, pu_syst) << endl;
	fOUTSTREAM << Form("matching lnN\t%5.3f/%5.3f\t-\t\t-\t\t%5.3f/%5.3f\t%5.3f/%5.3f",
					   match_syst_up, match_syst_dn, match_syst_up, match_syst_dn, match_syst_up, match_syst_dn) << endl;
	fOUTSTREAM << Form("scale    lnN\t%5.3f/%5.3f\t-\t\t-\t\t%5.3f/%5.3f\t%5.3f/%5.3f",
					   scale_syst_up, scale_syst_dn, scale_syst_up, scale_syst_dn, scale_syst_up, scale_syst_dn) << endl;
	fOUTSTREAM << endl;
	fOUTSTREAM.close();
	
	datacard = outputdir + "datacard_ssdl_ee" + chargeString + ".txt";
	fOUTSTREAM.open(datacard.Data(), ios::trunc);
	fOUTSTREAM <<      "#=========================================================================================" << endl;
	fOUTSTREAM <<      "# Systematics table for ttW/Z analysis, same-sign channel, ee only" << endl;
	fOUTSTREAM << Form("# Generated on: %s ", asctime(timeinfo)) << endl;
	fOUTSTREAM <<      "# Copy between the dashed lines for datacard" << endl;
	fOUTSTREAM <<      "#-----------------------------------------------------------------------------------------" << endl;
	fOUTSTREAM <<      "imax 1" << endl;
	fOUTSTREAM <<      "jmax 4" << endl;
	fOUTSTREAM <<      "kmax *" << endl;
	fOUTSTREAM << endl << endl;
	fOUTSTREAM <<      "bin\t\t1" << endl;
	if (gFullDataBlind)
		fOUTSTREAM << Form("observation\t%d", 999) << endl;
	else
		fOUTSTREAM << Form("observation\t%d", ttwzpreds["Normal"].obs_ee) << endl;
	fOUTSTREAM << endl << endl;
	fOUTSTREAM <<      "bin\t\t1\t\t1\t\t1\t\t1\t\t1" << endl;
	fOUTSTREAM <<      "process\t\tttWZ\t\tfake\t\tcmid\t\twz\t\trare" << endl;
	fOUTSTREAM <<      "process\t\t0\t\t1\t\t2\t\t3\t\t4" << endl;
	fOUTSTREAM << Form("rate\t\t%5.3f\t\t%5.3f\t\t%5.3f\t\t%5.3f\t\t%5.3f",
					   ttwzpreds["Normal"].ttwz_ee, ttwzpreds["Normal"].fake_ee, ttwzpreds["Normal"].cmid_ee, ttwzpreds["Normal"].wz_ee, ttwzpreds["Normal"].rare_ee) << endl;
	fOUTSTREAM << Form("# separate signal yields: %5.3f (ttW), %5.3f (ttZ)", ttwzpreds["Normal"].ttw_ee, ttwzpreds["Normal"].ttz_ee) << endl;
	fOUTSTREAM << endl << endl;
	fOUTSTREAM <<      "#syst" << endl;
	fOUTSTREAM <<      "lumi     lnN\t"+lumiError+"\t\t"+lumiError+"\t\t"+lumiError+"\t\t"+lumiError+"\t\t"+lumiError+"" << endl;
	fOUTSTREAM << Form("bgUncfak lnN\t-\t\t%5.3f\t\t-\t\t-\t\t-", 1.0+ttwzpreds["Normal"].fake_err_ee/ttwzpreds["Normal"].fake_ee) << endl;
	fOUTSTREAM << Form("bgUnccmi lnN\t-\t\t-\t\t%5.3f\t\t-\t\t-", 1.0+ttwzpreds["Normal"].cmid_err_ee/ttwzpreds["Normal"].cmid_ee) << endl;
	fOUTSTREAM << Form("bgUncwz  lnN\t-\t\t-\t\t-\t\t%5.3f\t\t-", 1.0+ttwzpreds["Normal"].wz_err_ee  /ttwzpreds["Normal"].wz_ee)   << endl;
	fOUTSTREAM << Form("bgUncrar lnN\t-\t\t-\t\t-\t\t-\t\t%5.3f", 1.0+ttwzpreds["Normal"].rare_err_ee/ttwzpreds["Normal"].rare_ee) << endl;
	fOUTSTREAM << Form("lept     lnN\t%5.3f/%5.3f\t-\t\t-\t\t%5.3f/%5.3f\t%5.3f/%5.3f",
					   1.0+(ttwzpreds["LepDown"].ttwz_ee-ttwzpreds["Normal"].ttwz_ee)/ttwzpreds["Normal"].ttwz_ee,
					   1.0+(ttwzpreds["LepUp"].ttwz_ee-ttwzpreds["Normal"].ttwz_ee)/ttwzpreds["Normal"].ttwz_ee,
					   1.0+(ttwzpreds["LepDown"].wz_ee  -ttwzpreds["Normal"].wz_ee  )/ttwzpreds["Normal"].wz_ee,
					   1.0+(ttwzpreds["LepUp"].wz_ee  -ttwzpreds["Normal"].wz_ee  )/ttwzpreds["Normal"].wz_ee,
					   1.0+(ttwzpreds["LepDown"].rare_ee-ttwzpreds["Normal"].rare_ee)/ttwzpreds["Normal"].rare_ee,
					   1.0+(ttwzpreds["LepUp"].rare_ee-ttwzpreds["Normal"].rare_ee)/ttwzpreds["Normal"].rare_ee) << endl;
	fOUTSTREAM << Form("btag     lnN\t%5.3f/%5.3f\t-\t\t-\t\t%5.3f/%5.3f\t%5.3f/%5.3f",
					   1.0+(ttwzpreds["BDown"].ttwz_ee-ttwzpreds["Normal"].ttwz_ee)/ttwzpreds["Normal"].ttwz_ee,
					   1.0+(ttwzpreds["BUp"].ttwz_ee-ttwzpreds["Normal"].ttwz_ee)/ttwzpreds["Normal"].ttwz_ee,
					   1.0+(ttwzpreds["BDown"].wz_ee  -ttwzpreds["Normal"].wz_ee  )/ttwzpreds["Normal"].wz_ee,
					   1.0+(ttwzpreds["BUp"].wz_ee  -ttwzpreds["Normal"].wz_ee  )/ttwzpreds["Normal"].wz_ee,
					   1.0+(ttwzpreds["BDown"].rare_ee-ttwzpreds["Normal"].rare_ee)/ttwzpreds["Normal"].rare_ee,
					   1.0+(ttwzpreds["BUp"].rare_ee-ttwzpreds["Normal"].rare_ee)/ttwzpreds["Normal"].rare_ee) << endl;
	fOUTSTREAM << Form("jes      lnN\t%5.3f/%5.3f\t-\t\t-\t\t%5.3f/%5.3f\t%5.3f/%5.3f",
					   1.0+(ttwzpreds["JetDown"].ttwz_ee-ttwzpreds["Normal"].ttwz_ee)/ttwzpreds["Normal"].ttwz_ee,
					   1.0+(ttwzpreds["JetUp"].ttwz_ee-ttwzpreds["Normal"].ttwz_ee)/ttwzpreds["Normal"].ttwz_ee,
					   1.0+(ttwzpreds["JetDown"].wz_ee  -ttwzpreds["Normal"].wz_ee  )/ttwzpreds["Normal"].wz_ee,
					   1.0+(ttwzpreds["JetUp"].wz_ee  -ttwzpreds["Normal"].wz_ee  )/ttwzpreds["Normal"].wz_ee,
					   1.0+(ttwzpreds["JetDown"].rare_ee-ttwzpreds["Normal"].rare_ee)/ttwzpreds["Normal"].rare_ee,
					   1.0+(ttwzpreds["JetUp"].rare_ee-ttwzpreds["Normal"].rare_ee)/ttwzpreds["Normal"].rare_ee) << endl;
	fOUTSTREAM << Form("jer      lnN\t%5.3f\t\t-\t\t-\t\t%5.3f\t\t%5.3f",
					   1.0+(ttwzpreds["JetSmear"].ttwz_ee-ttwzpreds["Normal"].ttwz_ee)/ttwzpreds["Normal"].ttwz_ee,
					   1.0+(ttwzpreds["JetSmear"].wz_ee  -ttwzpreds["Normal"].wz_ee  )/ttwzpreds["Normal"].wz_ee,
					   1.0+(ttwzpreds["JetSmear"].rare_ee-ttwzpreds["Normal"].rare_ee)/ttwzpreds["Normal"].rare_ee) << endl;
	fOUTSTREAM << Form("pu       lnN\t%5.3f\t\t-\t\t-\t\t%5.3f\t\t%5.3f",
					   pu_syst, pu_syst, pu_syst) << endl;
	fOUTSTREAM << Form("matching lnN\t%5.3f/%5.3f\t-\t\t-\t\t%5.3f/%5.3f\t%5.3f/%5.3f",
					   match_syst_up, match_syst_dn, match_syst_up, match_syst_dn, match_syst_up, match_syst_dn) << endl;
	fOUTSTREAM << Form("scale    lnN\t%5.3f/%5.3f\t-\t\t-\t\t%5.3f/%5.3f\t%5.3f/%5.3f",
					   scale_syst_up, scale_syst_dn, scale_syst_up, scale_syst_dn, scale_syst_up, scale_syst_dn) << endl;
	fOUTSTREAM << endl;
	fOUTSTREAM.close();
	
	datacard = outputdir + "datacard_ssdl_3channels" + chargeString + ".txt";
	fOUTSTREAM.open(datacard.Data(), ios::trunc);
	fOUTSTREAM <<      "#=========================================================================================" << endl;
	fOUTSTREAM <<      "# Systematics table for ttW/Z analysis, same-sign channel, subchannels" << endl;
	fOUTSTREAM << Form("# Generated on: %s ", asctime(timeinfo)) << endl;
	fOUTSTREAM <<      "# Copy between the dashed lines for datacard" << endl;
	fOUTSTREAM <<      "#-----------------------------------------------------------------------------------------" << endl;
	fOUTSTREAM <<      "imax 3" << endl;
	fOUTSTREAM <<      "jmax 4" << endl;
	fOUTSTREAM <<      "kmax *" << endl;
	fOUTSTREAM << endl << endl;
	fOUTSTREAM <<      "bin\t\t1\t2\t3" << endl;
	if (gFullDataBlind)
		fOUTSTREAM << Form("observation\t%d\t%d\t%d", 999, 999, 999) << endl;
	else
		fOUTSTREAM << Form("observation\t%d\t%d\t%d", ttwzpreds["Normal"].obs_mm, ttwzpreds["Normal"].obs_em, ttwzpreds["Normal"].obs_ee) << endl;
	fOUTSTREAM << endl << endl;
	fOUTSTREAM <<      "bin\t\t1\t\t1\t\t1\t\t1\t\t1\t\t2\t\t2\t\t2\t\t2\t\t2\t\t3\t\t3\t\t3\t\t3\t\t3" << endl;
	fOUTSTREAM <<      "process\t\tttWZ\t\tfake\t\tcmid\t\twz\t\trare\t\tttWZ\t\tfake\t\tcmid\t\twz\t\trare\t\tttWZ\t\tfake\t\tcmid\t\twz\t\trare" << endl;
	fOUTSTREAM <<      "process\t\t0\t\t1\t\t2\t\t3\t\t4\t\t0\t\t1\t\t2\t\t3\t\t4\t\t0\t\t1\t\t2\t\t3\t\t4" << endl;
	fOUTSTREAM << Form("rate\t\t%5.3f\t\t%5.3f\t\t%5.3f\t\t%5.3f\t\t%5.3f\t\t%5.3f\t\t%5.3f\t\t%5.3f\t\t%5.3f\t\t%5.3f\t\t%5.3f\t\t%5.3f\t\t%5.3f\t\t%5.3f\t\t%5.3f",
					   ttwzpreds["Normal"].ttwz_mm, ttwzpreds["Normal"].fake_mm, 0.0                   , ttwzpreds["Normal"].wz_mm, ttwzpreds["Normal"].rare_mm,
					   ttwzpreds["Normal"].ttwz_em, ttwzpreds["Normal"].fake_em, ttwzpreds["Normal"].cmid_em, ttwzpreds["Normal"].wz_em, ttwzpreds["Normal"].rare_em,
					   ttwzpreds["Normal"].ttwz_ee, ttwzpreds["Normal"].fake_ee, ttwzpreds["Normal"].cmid_ee, ttwzpreds["Normal"].wz_ee, ttwzpreds["Normal"].rare_ee) << endl;
	fOUTSTREAM << Form("# separate signal yields: [mm] %5.3f (ttW), %5.3f (ttZ); [em] %5.3f (ttW), %5.3f (ttZ); [ee] %5.3f (ttW), %5.3f (ttZ);",
					   ttwzpreds["Normal"].ttw_mm, ttwzpreds["Normal"].ttz_mm,
					   ttwzpreds["Normal"].ttw_em, ttwzpreds["Normal"].ttz_em,
					   ttwzpreds["Normal"].ttw_ee, ttwzpreds["Normal"].ttz_ee) << endl;
	fOUTSTREAM << endl << endl;
	fOUTSTREAM <<      "#syst" << endl;
	fOUTSTREAM <<      "lumi     lnN\t"+lumiError+"\t\t"+lumiError+"\t\t"+lumiError+"\t\t"+lumiError+"\t\t"+lumiError+"\t\t"+lumiError+"\t\t"+lumiError+"\t\t"+lumiError+"\t\t"+lumiError+"\t\t"+lumiError+"\t\t"+lumiError+"\t\t"+lumiError+"\t\t"+lumiError+"\t\t"+lumiError+"\t\t"+lumiError+"" << endl;
	fOUTSTREAM << Form("bgUncfak lnN\t-\t\t%5.3f\t\t-\t\t-\t\t-\t\t-\t\t%5.3f\t\t-\t\t-\t\t-\t\t-\t\t%5.3f\t\t-\t\t-\t\t-",
					   1.0+ttwzpreds["Normal"].fake_err_mm/ttwzpreds["Normal"].fake_mm,
					   1.0+ttwzpreds["Normal"].fake_err_em/ttwzpreds["Normal"].fake_em,
					   1.0+ttwzpreds["Normal"].fake_err_ee/ttwzpreds["Normal"].fake_ee) << endl;
	fOUTSTREAM << Form("bgUnccmi lnN\t-\t\t-\t\t-\t\t-\t\t-\t\t-\t\t-\t\t%5.3f\t\t-\t\t-\t\t-\t\t-\t\t%5.3f\t\t-\t\t-",
					   1.0+ttwzpreds["Normal"].cmid_err_em/ttwzpreds["Normal"].cmid_em,
					   1.0+ttwzpreds["Normal"].cmid_err_ee/ttwzpreds["Normal"].cmid_ee) << endl;
	fOUTSTREAM << Form("bgUncwz  lnN\t-\t\t-\t\t-\t\t%5.3f\t\t-\t\t-\t\t-\t\t-\t\t%5.3f\t\t-\t\t-\t\t-\t\t-\t\t%5.3f\t\t-",
					   1.0+ttwzpreds["Normal"].wz_err_mm  /ttwzpreds["Normal"].wz_mm,
					   1.0+ttwzpreds["Normal"].wz_err_em  /ttwzpreds["Normal"].wz_em,
					   1.0+ttwzpreds["Normal"].wz_err_ee  /ttwzpreds["Normal"].wz_ee)   << endl;
	fOUTSTREAM << Form("bgUncrar lnN\t-\t\t-\t\t-\t\t-\t\t%5.3f\t\t-\t\t-\t\t-\t\t-\t\t%5.3f\t\t-\t\t-\t\t-\t\t-\t\t%5.3f",
					   1.0+ttwzpreds["Normal"].rare_err_mm/ttwzpreds["Normal"].rare_mm,
					   1.0+ttwzpreds["Normal"].rare_err_em/ttwzpreds["Normal"].rare_em,
					   1.0+ttwzpreds["Normal"].rare_err_ee/ttwzpreds["Normal"].rare_ee) << endl;
	fOUTSTREAM << Form("lept     lnN\t%5.3f/%5.3f\t-\t\t-\t\t%5.3f/%5.3f\t%5.3f/%5.3f\t%5.3f/%5.3f\t-\t\t-\t\t%5.3f/%5.3f\t%5.3f/%5.3f\t%5.3f/%5.3f\t-\t\t-\t\t%5.3f/%5.3f\t%5.3f/%5.3f",
					   1.0+(ttwzpreds["LepDown"].ttwz_mm-ttwzpreds["Normal"].ttwz_mm)/ttwzpreds["Normal"].ttwz_mm,
					   1.0+(ttwzpreds["LepUp"].ttwz_mm-ttwzpreds["Normal"].ttwz_mm)/ttwzpreds["Normal"].ttwz_mm,
					   1.0+(ttwzpreds["LepDown"].wz_mm  -ttwzpreds["Normal"].wz_mm  )/ttwzpreds["Normal"].wz_mm,
					   1.0+(ttwzpreds["LepUp"].wz_mm  -ttwzpreds["Normal"].wz_mm  )/ttwzpreds["Normal"].wz_mm,
					   1.0+(ttwzpreds["LepDown"].rare_mm-ttwzpreds["Normal"].rare_mm)/ttwzpreds["Normal"].rare_mm,
					   1.0+(ttwzpreds["LepUp"].rare_mm-ttwzpreds["Normal"].rare_mm)/ttwzpreds["Normal"].rare_mm,
					   1.0+(ttwzpreds["LepDown"].ttwz_em-ttwzpreds["Normal"].ttwz_em)/ttwzpreds["Normal"].ttwz_em,
					   1.0+(ttwzpreds["LepUp"].ttwz_em-ttwzpreds["Normal"].ttwz_em)/ttwzpreds["Normal"].ttwz_em,
					   1.0+(ttwzpreds["LepDown"].wz_em  -ttwzpreds["Normal"].wz_em  )/ttwzpreds["Normal"].wz_em,
					   1.0+(ttwzpreds["LepUp"].wz_em  -ttwzpreds["Normal"].wz_em  )/ttwzpreds["Normal"].wz_em,
					   1.0+(ttwzpreds["LepDown"].rare_em-ttwzpreds["Normal"].rare_em)/ttwzpreds["Normal"].rare_em,
					   1.0+(ttwzpreds["LepUp"].rare_em-ttwzpreds["Normal"].rare_em)/ttwzpreds["Normal"].rare_em,
					   1.0+(ttwzpreds["LepDown"].ttwz_ee-ttwzpreds["Normal"].ttwz_ee)/ttwzpreds["Normal"].ttwz_ee,
					   1.0+(ttwzpreds["LepUp"].ttwz_ee-ttwzpreds["Normal"].ttwz_ee)/ttwzpreds["Normal"].ttwz_ee,
					   1.0+(ttwzpreds["LepDown"].wz_ee  -ttwzpreds["Normal"].wz_ee  )/ttwzpreds["Normal"].wz_ee,
					   1.0+(ttwzpreds["LepUp"].wz_ee  -ttwzpreds["Normal"].wz_ee  )/ttwzpreds["Normal"].wz_ee,
					   1.0+(ttwzpreds["LepDown"].rare_ee-ttwzpreds["Normal"].rare_ee)/ttwzpreds["Normal"].rare_ee,
					   1.0+(ttwzpreds["LepUp"].rare_ee-ttwzpreds["Normal"].rare_ee)/ttwzpreds["Normal"].rare_ee) << endl;
	fOUTSTREAM << Form("btag     lnN\t%5.3f/%5.3f\t-\t\t-\t\t%5.3f/%5.3f\t%5.3f/%5.3f\t%5.3f/%5.3f\t-\t\t-\t\t%5.3f/%5.3f\t%5.3f/%5.3f\t%5.3f/%5.3f\t-\t\t-\t\t%5.3f/%5.3f\t%5.3f/%5.3f",
					   1.0+(ttwzpreds["BDown"].ttwz_mm-ttwzpreds["Normal"].ttwz_mm)/ttwzpreds["Normal"].ttwz_mm,
					   1.0+(ttwzpreds["BUp"].ttwz_mm-ttwzpreds["Normal"].ttwz_mm)/ttwzpreds["Normal"].ttwz_mm,
					   1.0+(ttwzpreds["BDown"].wz_mm  -ttwzpreds["Normal"].wz_mm  )/ttwzpreds["Normal"].wz_mm,
					   1.0+(ttwzpreds["BUp"].wz_mm  -ttwzpreds["Normal"].wz_mm  )/ttwzpreds["Normal"].wz_mm,
					   1.0+(ttwzpreds["BDown"].rare_mm-ttwzpreds["Normal"].rare_mm)/ttwzpreds["Normal"].rare_mm,
					   1.0+(ttwzpreds["BUp"].rare_mm-ttwzpreds["Normal"].rare_mm)/ttwzpreds["Normal"].rare_mm,
					   1.0+(ttwzpreds["BDown"].ttwz_em-ttwzpreds["Normal"].ttwz_em)/ttwzpreds["Normal"].ttwz_em,
					   1.0+(ttwzpreds["BUp"].ttwz_em-ttwzpreds["Normal"].ttwz_em)/ttwzpreds["Normal"].ttwz_em,
					   1.0+(ttwzpreds["BDown"].wz_em  -ttwzpreds["Normal"].wz_em  )/ttwzpreds["Normal"].wz_em,
					   1.0+(ttwzpreds["BUp"].wz_em  -ttwzpreds["Normal"].wz_em  )/ttwzpreds["Normal"].wz_em,
					   1.0+(ttwzpreds["BDown"].rare_em-ttwzpreds["Normal"].rare_em)/ttwzpreds["Normal"].rare_em,
					   1.0+(ttwzpreds["BUp"].rare_em-ttwzpreds["Normal"].rare_em)/ttwzpreds["Normal"].rare_em,
					   1.0+(ttwzpreds["BDown"].ttwz_ee-ttwzpreds["Normal"].ttwz_ee)/ttwzpreds["Normal"].ttwz_ee,
					   1.0+(ttwzpreds["BUp"].ttwz_ee-ttwzpreds["Normal"].ttwz_ee)/ttwzpreds["Normal"].ttwz_ee,
					   1.0+(ttwzpreds["BDown"].wz_ee  -ttwzpreds["Normal"].wz_ee  )/ttwzpreds["Normal"].wz_ee,
					   1.0+(ttwzpreds["BUp"].wz_ee  -ttwzpreds["Normal"].wz_ee  )/ttwzpreds["Normal"].wz_ee,
					   1.0+(ttwzpreds["BDown"].rare_ee-ttwzpreds["Normal"].rare_ee)/ttwzpreds["Normal"].rare_ee,
					   1.0+(ttwzpreds["BUp"].rare_ee-ttwzpreds["Normal"].rare_ee)/ttwzpreds["Normal"].rare_ee) << endl;
	fOUTSTREAM << Form("jes      lnN\t%5.3f/%5.3f\t-\t\t-\t\t%5.3f/%5.3f\t%5.3f/%5.3f\t%5.3f/%5.3f\t-\t\t-\t\t%5.3f/%5.3f\t%5.3f/%5.3f\t%5.3f/%5.3f\t-\t\t-\t\t%5.3f/%5.3f\t%5.3f/%5.3f",
					   1.0+(ttwzpreds["JetDown"].ttwz_mm-ttwzpreds["Normal"].ttwz_mm)/ttwzpreds["Normal"].ttwz_mm,
					   1.0+(ttwzpreds["JetUp"].ttwz_mm-ttwzpreds["Normal"].ttwz_mm)/ttwzpreds["Normal"].ttwz_mm,
					   1.0+(ttwzpreds["JetDown"].wz_mm  -ttwzpreds["Normal"].wz_mm  )/ttwzpreds["Normal"].wz_mm,
					   1.0+(ttwzpreds["JetUp"].wz_mm  -ttwzpreds["Normal"].wz_mm  )/ttwzpreds["Normal"].wz_mm,
					   1.0+(ttwzpreds["JetDown"].rare_mm-ttwzpreds["Normal"].rare_mm)/ttwzpreds["Normal"].rare_mm,
					   1.0+(ttwzpreds["JetUp"].rare_mm-ttwzpreds["Normal"].rare_mm)/ttwzpreds["Normal"].rare_mm,
					   1.0+(ttwzpreds["JetDown"].ttwz_em-ttwzpreds["Normal"].ttwz_em)/ttwzpreds["Normal"].ttwz_em,
					   1.0+(ttwzpreds["JetUp"].ttwz_em-ttwzpreds["Normal"].ttwz_em)/ttwzpreds["Normal"].ttwz_em,
					   1.0+(ttwzpreds["JetDown"].wz_em  -ttwzpreds["Normal"].wz_em  )/ttwzpreds["Normal"].wz_em,
					   1.0+(ttwzpreds["JetUp"].wz_em  -ttwzpreds["Normal"].wz_em  )/ttwzpreds["Normal"].wz_em,
					   1.0+(ttwzpreds["JetDown"].rare_em-ttwzpreds["Normal"].rare_em)/ttwzpreds["Normal"].rare_em,
					   1.0+(ttwzpreds["JetUp"].rare_em-ttwzpreds["Normal"].rare_em)/ttwzpreds["Normal"].rare_em,
					   1.0+(ttwzpreds["JetDown"].ttwz_ee-ttwzpreds["Normal"].ttwz_ee)/ttwzpreds["Normal"].ttwz_ee,
					   1.0+(ttwzpreds["JetUp"].ttwz_ee-ttwzpreds["Normal"].ttwz_ee)/ttwzpreds["Normal"].ttwz_ee,
					   1.0+(ttwzpreds["JetDown"].wz_ee  -ttwzpreds["Normal"].wz_ee  )/ttwzpreds["Normal"].wz_ee,
					   1.0+(ttwzpreds["JetUp"].wz_ee  -ttwzpreds["Normal"].wz_ee  )/ttwzpreds["Normal"].wz_ee,
					   1.0+(ttwzpreds["JetDown"].rare_ee-ttwzpreds["Normal"].rare_ee)/ttwzpreds["Normal"].rare_ee,
					   1.0+(ttwzpreds["JetUp"].rare_ee-ttwzpreds["Normal"].rare_ee)/ttwzpreds["Normal"].rare_ee) << endl;
	fOUTSTREAM << Form("jer      lnN\t%5.3f\t\t-\t\t-\t\t%5.3f\t\t%5.3f\t\t%5.3f\t\t-\t\t-\t\t%5.3f\t\t%5.3f\t\t%5.3f\t\t-\t\t-\t\t%5.3f\t\t%5.3f",
					   1.0+(ttwzpreds["JetSmear"].ttwz_mm-ttwzpreds["Normal"].ttwz_mm)/ttwzpreds["Normal"].ttwz_mm,
					   1.0+(ttwzpreds["JetSmear"].wz_mm  -ttwzpreds["Normal"].wz_mm  )/ttwzpreds["Normal"].wz_mm,
					   1.0+(ttwzpreds["JetSmear"].rare_mm-ttwzpreds["Normal"].rare_mm)/ttwzpreds["Normal"].rare_mm,
					   1.0+(ttwzpreds["JetSmear"].ttwz_em-ttwzpreds["Normal"].ttwz_em)/ttwzpreds["Normal"].ttwz_em,
					   1.0+(ttwzpreds["JetSmear"].wz_em  -ttwzpreds["Normal"].wz_em  )/ttwzpreds["Normal"].wz_em,
					   1.0+(ttwzpreds["JetSmear"].rare_em-ttwzpreds["Normal"].rare_em)/ttwzpreds["Normal"].rare_em,
					   1.0+(ttwzpreds["JetSmear"].ttwz_ee-ttwzpreds["Normal"].ttwz_ee)/ttwzpreds["Normal"].ttwz_ee,
					   1.0+(ttwzpreds["JetSmear"].wz_ee  -ttwzpreds["Normal"].wz_ee  )/ttwzpreds["Normal"].wz_ee,
					   1.0+(ttwzpreds["JetSmear"].rare_ee-ttwzpreds["Normal"].rare_ee)/ttwzpreds["Normal"].rare_ee) << endl;
	fOUTSTREAM << Form("pu       lnN\t%5.3f\t\t-\t\t-\t\t%5.3f\t\t%5.3f\t\t%5.3f\t\t-\t\t-\t\t%5.3f\t\t%5.3f\t\t%5.3f\t\t-\t\t-\t\t%5.3f\t\t%5.3f",
					   pu_syst, pu_syst, pu_syst, pu_syst, pu_syst, pu_syst, pu_syst, pu_syst, pu_syst) << endl;
	fOUTSTREAM << Form("matching lnN\t%5.3f/%5.3f\t-\t\t-\t\t%5.3f/%5.3f\t%5.3f/%5.3f\t%5.3f/%5.3f\t-\t\t-\t\t%5.3f/%5.3f\t%5.3f/%5.3f\t%5.3f/%5.3f\t-\t\t-\t\t%5.3f/%5.3f\t%5.3f/%5.3f",
					   match_syst_up, match_syst_dn, match_syst_up, match_syst_dn, match_syst_up, match_syst_dn, match_syst_up, match_syst_dn, match_syst_up, match_syst_dn,
					   match_syst_up, match_syst_dn, match_syst_up, match_syst_dn, match_syst_up, match_syst_dn, match_syst_up, match_syst_dn) << endl;
	fOUTSTREAM << Form("scale    lnN\t%5.3f/%5.3f\t-\t\t-\t\t%5.3f/%5.3f\t%5.3f/%5.3f\t%5.3f/%5.3f\t-\t\t-\t\t%5.3f/%5.3f\t%5.3f/%5.3f\t%5.3f/%5.3f\t-\t\t-\t\t%5.3f/%5.3f\t%5.3f/%5.3f",
					   scale_syst_up, scale_syst_dn, scale_syst_up, scale_syst_dn, scale_syst_up, scale_syst_dn, scale_syst_up, scale_syst_dn, scale_syst_up, scale_syst_dn,
					   scale_syst_up, scale_syst_dn, scale_syst_up, scale_syst_dn, scale_syst_up, scale_syst_dn, scale_syst_up, scale_syst_dn) << endl;
	fOUTSTREAM << endl;
	fOUTSTREAM.close();
	
	
	datacard = outputdir + "datacard_ssdl_3channels_ttW" + chargeString + ".txt";
	fOUTSTREAM.open(datacard.Data(), ios::trunc);
	fOUTSTREAM <<      "#=========================================================================================" << endl;
	fOUTSTREAM <<      "# Systematics table for ttW/Z analysis, same-sign channel, subchannels" << endl;
	fOUTSTREAM << Form("# Generated on: %s ", asctime(timeinfo)) << endl;
	fOUTSTREAM <<      "# Copy between the dashed lines for datacard" << endl;
	fOUTSTREAM <<      "#-----------------------------------------------------------------------------------------" << endl;
	fOUTSTREAM <<      "imax 3" << endl;
	fOUTSTREAM <<      "jmax 5" << endl;
	fOUTSTREAM <<      "kmax *" << endl;
	fOUTSTREAM << endl << endl;
	fOUTSTREAM <<      "bin\t\t1\t2\t3" << endl;
	if (gFullDataBlind)
		fOUTSTREAM << Form("observation\t%d\t%d\t%d", 999, 999, 999) << endl;
	else
		fOUTSTREAM << Form("observation\t%d\t%d\t%d", ttwzpreds["Normal"].obs_mm, ttwzpreds["Normal"].obs_em, ttwzpreds["Normal"].obs_ee) << endl;
	fOUTSTREAM << endl << endl;
	fOUTSTREAM <<      "bin\t\t1\t\t1\t\t1\t\t1\t\t1\t\t1\t\t2\t\t2\t\t2\t\t2\t\t2\t\t2\t\t3\t\t3\t\t3\t\t3\t\t3\t\t3" << endl;
	fOUTSTREAM <<      "process\t\tttW\t\tttZ\t\tfake\t\tcmid\t\twz\t\trare\t\tttW\t\tttZ\t\tfake\t\tcmid\t\twz\t\trare\t\tttW\t\tttZ\t\tfake\t\tcmid\t\twz\t\trare" << endl;
	fOUTSTREAM <<      "process\t\t0\t\t1\t\t2\t\t3\t\t4\t\t5\t\t0\t\t1\t\t2\t\t3\t\t4\t\t5\t\t0\t\t1\t\t2\t\t3\t\t4\t\t5" << endl;
	fOUTSTREAM << Form("rate\t\t%5.3f\t\t%5.3f\t\t%5.3f\t\t%5.3f\t\t%5.3f\t\t%5.3f\t\t%5.3f\t\t%5.3f\t\t%5.3f\t\t%5.3f\t\t%5.3f\t\t%5.3f\t\t%5.3f\t\t%5.3f\t\t%5.3f\t\t%5.3f\t\t%5.3f\t\t%5.3f",
					   ttwzpreds["Normal"].ttw_mm, ttwzpreds["Normal"].ttz_mm, ttwzpreds["Normal"].fake_mm, 0.0                   , ttwzpreds["Normal"].wz_mm, ttwzpreds["Normal"].rare_mm,
					   ttwzpreds["Normal"].ttw_em, ttwzpreds["Normal"].ttz_em, ttwzpreds["Normal"].fake_em, ttwzpreds["Normal"].cmid_em, ttwzpreds["Normal"].wz_em, ttwzpreds["Normal"].rare_em,
					   ttwzpreds["Normal"].ttw_ee, ttwzpreds["Normal"].ttz_ee, ttwzpreds["Normal"].fake_ee, ttwzpreds["Normal"].cmid_ee, ttwzpreds["Normal"].wz_ee, ttwzpreds["Normal"].rare_ee) << endl;
	fOUTSTREAM << endl << endl;
	fOUTSTREAM <<      "#syst" << endl;
	fOUTSTREAM <<      "lumi     lnN\t"+lumiError+"\t\t"+lumiError+"\t\t"+lumiError+"\t\t"+lumiError+"\t\t"+lumiError+"\t\t"+lumiError+"\t\t"+lumiError+"\t\t"+lumiError+"\t\t"+lumiError+"\t\t"+lumiError+"\t\t"+lumiError+"\t\t"+lumiError+"\t\t"+lumiError+"\t\t"+lumiError+"\t\t"+lumiError+"\t\t"+lumiError+"\t\t"+lumiError+"\t\t"+lumiError+"" << endl;
	fOUTSTREAM << Form("bgUncttz lnN\t-\t\t%5.3f\t\t-\t\t-\t\t-\t\t-\t\t-\t\t-\t\t%5.3f\t\t-\t\t-\t\t-\t\t-\t\t-\t\t%5.3f\t\t-\t\t-\t\t-",
					   1.0+ttwzpreds["Normal"].ttz_err_mm/ttwzpreds["Normal"].ttz_mm,
					   1.0+ttwzpreds["Normal"].ttz_err_em/ttwzpreds["Normal"].ttz_em,
					   1.0+ttwzpreds["Normal"].ttz_err_ee/ttwzpreds["Normal"].ttz_ee)   << endl;
	fOUTSTREAM << Form("bgUncfak lnN\t-\t\t-\t\t%5.3f\t\t-\t\t-\t\t-\t\t-\t\t-\t\t%5.3f\t\t-\t\t-\t\t-\t\t-\t\t-\t\t%5.3f\t\t-\t\t-\t\t-",
					   1.0+ttwzpreds["Normal"].fake_err_mm/ttwzpreds["Normal"].fake_mm,
					   1.0+ttwzpreds["Normal"].fake_err_em/ttwzpreds["Normal"].fake_em,
					   1.0+ttwzpreds["Normal"].fake_err_ee/ttwzpreds["Normal"].fake_ee) << endl;
	fOUTSTREAM << Form("bgUnccmi lnN\t-\t\t-\t\t-\t\t-\t\t-\t\t-\t\t-\t\t-\t\t-\t\t%5.3f\t\t-\t\t-\t\t-\t\t-\t\t-\t\t%5.3f\t\t-\t\t-",
					   1.0+ttwzpreds["Normal"].cmid_err_em/ttwzpreds["Normal"].cmid_em,
					   1.0+ttwzpreds["Normal"].cmid_err_ee/ttwzpreds["Normal"].cmid_ee) << endl;
	fOUTSTREAM << Form("bgUncwz  lnN\t-\t\t-\t\t-\t\t-\t\t%5.3f\t\t-\t\t-\t\t-\t\t-\t\t-\t\t%5.3f\t\t-\t\t-\t\t-\t\t-\t\t-\t\t%5.3f\t\t-",
					   1.0+ttwzpreds["Normal"].wz_err_mm  /ttwzpreds["Normal"].wz_mm,
					   1.0+ttwzpreds["Normal"].wz_err_em  /ttwzpreds["Normal"].wz_em,
					   1.0+ttwzpreds["Normal"].wz_err_ee  /ttwzpreds["Normal"].wz_ee)   << endl;
	fOUTSTREAM << Form("bgUncrar lnN\t-\t\t-\t\t-\t\t-\t\t-\t\t%5.3f\t\t-\t\t-\t\t-\t\t-\t\t-\t\t%5.3f\t\t-\t\t-\t\t-\t\t-\t\t-\t\t%5.3f",
					   1.0+ttwzpreds["Normal"].rare_err_mm/ttwzpreds["Normal"].rare_mm,
					   1.0+ttwzpreds["Normal"].rare_err_em/ttwzpreds["Normal"].rare_em,
					   1.0+ttwzpreds["Normal"].rare_err_ee/ttwzpreds["Normal"].rare_ee) << endl;
	fOUTSTREAM << Form("lept     lnN\t%5.3f/%5.3f\t%5.3f/%5.3f\t-\t\t-\t\t%5.3f/%5.3f\t%5.3f/%5.3f\t%5.3f/%5.3f\t%5.3f/%5.3f\t-\t\t-\t\t%5.3f/%5.3f\t%5.3f/%5.3f\t%5.3f/%5.3f\t%5.3f/%5.3f\t-\t\t-\t\t%5.3f/%5.3f\t%5.3f/%5.3f",
					   1.0+(ttwzpreds["LepDown"].ttwz_mm-ttwzpreds["Normal"].ttwz_mm)/ttwzpreds["Normal"].ttwz_mm,
					   1.0+(ttwzpreds["LepUp"].ttwz_mm-ttwzpreds["Normal"].ttwz_mm)/ttwzpreds["Normal"].ttwz_mm,
					   1.0+(ttwzpreds["LepDown"].ttwz_mm-ttwzpreds["Normal"].ttwz_mm)/ttwzpreds["Normal"].ttwz_mm,
					   1.0+(ttwzpreds["LepUp"].ttwz_mm-ttwzpreds["Normal"].ttwz_mm)/ttwzpreds["Normal"].ttwz_mm,
					   1.0+(ttwzpreds["LepDown"].wz_mm  -ttwzpreds["Normal"].wz_mm  )/ttwzpreds["Normal"].wz_mm,
					   1.0+(ttwzpreds["LepUp"].wz_mm  -ttwzpreds["Normal"].wz_mm  )/ttwzpreds["Normal"].wz_mm,
					   1.0+(ttwzpreds["LepDown"].rare_mm-ttwzpreds["Normal"].rare_mm)/ttwzpreds["Normal"].rare_mm,
					   1.0+(ttwzpreds["LepUp"].rare_mm-ttwzpreds["Normal"].rare_mm)/ttwzpreds["Normal"].rare_mm,
					   1.0+(ttwzpreds["LepDown"].ttwz_em-ttwzpreds["Normal"].ttwz_em)/ttwzpreds["Normal"].ttwz_em,
					   1.0+(ttwzpreds["LepUp"].ttwz_em-ttwzpreds["Normal"].ttwz_em)/ttwzpreds["Normal"].ttwz_em,
					   1.0+(ttwzpreds["LepDown"].ttwz_em-ttwzpreds["Normal"].ttwz_em)/ttwzpreds["Normal"].ttwz_em,
					   1.0+(ttwzpreds["LepUp"].ttwz_em-ttwzpreds["Normal"].ttwz_em)/ttwzpreds["Normal"].ttwz_em,
					   1.0+(ttwzpreds["LepDown"].wz_em  -ttwzpreds["Normal"].wz_em  )/ttwzpreds["Normal"].wz_em,
					   1.0+(ttwzpreds["LepUp"].wz_em  -ttwzpreds["Normal"].wz_em  )/ttwzpreds["Normal"].wz_em,
					   1.0+(ttwzpreds["LepDown"].rare_em-ttwzpreds["Normal"].rare_em)/ttwzpreds["Normal"].rare_em,
					   1.0+(ttwzpreds["LepUp"].rare_em-ttwzpreds["Normal"].rare_em)/ttwzpreds["Normal"].rare_em,
					   1.0+(ttwzpreds["LepDown"].ttwz_ee-ttwzpreds["Normal"].ttwz_ee)/ttwzpreds["Normal"].ttwz_ee,
					   1.0+(ttwzpreds["LepUp"].ttwz_ee-ttwzpreds["Normal"].ttwz_ee)/ttwzpreds["Normal"].ttwz_ee,
					   1.0+(ttwzpreds["LepDown"].ttwz_ee-ttwzpreds["Normal"].ttwz_ee)/ttwzpreds["Normal"].ttwz_ee,
					   1.0+(ttwzpreds["LepUp"].ttwz_ee-ttwzpreds["Normal"].ttwz_ee)/ttwzpreds["Normal"].ttwz_ee,
					   1.0+(ttwzpreds["LepDown"].wz_ee  -ttwzpreds["Normal"].wz_ee  )/ttwzpreds["Normal"].wz_ee,
					   1.0+(ttwzpreds["LepUp"].wz_ee  -ttwzpreds["Normal"].wz_ee  )/ttwzpreds["Normal"].wz_ee,
					   1.0+(ttwzpreds["LepDown"].rare_ee-ttwzpreds["Normal"].rare_ee)/ttwzpreds["Normal"].rare_ee,
					   1.0+(ttwzpreds["LepUp"].rare_ee-ttwzpreds["Normal"].rare_ee)/ttwzpreds["Normal"].rare_ee) << endl;
	fOUTSTREAM << Form("btag     lnN\t%5.3f/%5.3f\t%5.3f/%5.3f\t-\t\t-\t\t%5.3f/%5.3f\t%5.3f/%5.3f\t%5.3f/%5.3f\t%5.3f/%5.3f\t-\t\t-\t\t%5.3f/%5.3f\t%5.3f/%5.3f\t%5.3f/%5.3f\t%5.3f/%5.3f\t-\t\t-\t\t%5.3f/%5.3f\t%5.3f/%5.3f",
					   1.0+(ttwzpreds["BDown"].ttwz_mm-ttwzpreds["Normal"].ttwz_mm)/ttwzpreds["Normal"].ttwz_mm,
					   1.0+(ttwzpreds["BUp"].ttwz_mm-ttwzpreds["Normal"].ttwz_mm)/ttwzpreds["Normal"].ttwz_mm,
					   1.0+(ttwzpreds["BDown"].ttwz_mm-ttwzpreds["Normal"].ttwz_mm)/ttwzpreds["Normal"].ttwz_mm,
					   1.0+(ttwzpreds["BUp"].ttwz_mm-ttwzpreds["Normal"].ttwz_mm)/ttwzpreds["Normal"].ttwz_mm,
					   1.0+(ttwzpreds["BDown"].wz_mm  -ttwzpreds["Normal"].wz_mm  )/ttwzpreds["Normal"].wz_mm,
					   1.0+(ttwzpreds["BUp"].wz_mm  -ttwzpreds["Normal"].wz_mm  )/ttwzpreds["Normal"].wz_mm,
					   1.0+(ttwzpreds["BDown"].rare_mm-ttwzpreds["Normal"].rare_mm)/ttwzpreds["Normal"].rare_mm,
					   1.0+(ttwzpreds["BUp"].rare_mm-ttwzpreds["Normal"].rare_mm)/ttwzpreds["Normal"].rare_mm,
					   1.0+(ttwzpreds["BDown"].ttwz_em-ttwzpreds["Normal"].ttwz_em)/ttwzpreds["Normal"].ttwz_em,
					   1.0+(ttwzpreds["BUp"].ttwz_em-ttwzpreds["Normal"].ttwz_em)/ttwzpreds["Normal"].ttwz_em,
					   1.0+(ttwzpreds["BDown"].ttwz_em-ttwzpreds["Normal"].ttwz_em)/ttwzpreds["Normal"].ttwz_em,
					   1.0+(ttwzpreds["BUp"].ttwz_em-ttwzpreds["Normal"].ttwz_em)/ttwzpreds["Normal"].ttwz_em,
					   1.0+(ttwzpreds["BDown"].wz_em  -ttwzpreds["Normal"].wz_em  )/ttwzpreds["Normal"].wz_em,
					   1.0+(ttwzpreds["BUp"].wz_em  -ttwzpreds["Normal"].wz_em  )/ttwzpreds["Normal"].wz_em,
					   1.0+(ttwzpreds["BDown"].rare_em-ttwzpreds["Normal"].rare_em)/ttwzpreds["Normal"].rare_em,
					   1.0+(ttwzpreds["BUp"].rare_em-ttwzpreds["Normal"].rare_em)/ttwzpreds["Normal"].rare_em,
					   1.0+(ttwzpreds["BDown"].ttwz_ee-ttwzpreds["Normal"].ttwz_ee)/ttwzpreds["Normal"].ttwz_ee,
					   1.0+(ttwzpreds["BUp"].ttwz_ee-ttwzpreds["Normal"].ttwz_ee)/ttwzpreds["Normal"].ttwz_ee,
					   1.0+(ttwzpreds["BDown"].ttwz_ee-ttwzpreds["Normal"].ttwz_ee)/ttwzpreds["Normal"].ttwz_ee,
					   1.0+(ttwzpreds["BUp"].ttwz_ee-ttwzpreds["Normal"].ttwz_ee)/ttwzpreds["Normal"].ttwz_ee,
					   1.0+(ttwzpreds["BDown"].wz_ee  -ttwzpreds["Normal"].wz_ee  )/ttwzpreds["Normal"].wz_ee,
					   1.0+(ttwzpreds["BUp"].wz_ee  -ttwzpreds["Normal"].wz_ee  )/ttwzpreds["Normal"].wz_ee,
					   1.0+(ttwzpreds["BDown"].rare_ee-ttwzpreds["Normal"].rare_ee)/ttwzpreds["Normal"].rare_ee,
					   1.0+(ttwzpreds["BUp"].rare_ee-ttwzpreds["Normal"].rare_ee)/ttwzpreds["Normal"].rare_ee) << endl;
	fOUTSTREAM << Form("jes      lnN\t%5.3f/%5.3f\t%5.3f/%5.3f\t-\t\t-\t\t%5.3f/%5.3f\t%5.3f/%5.3f\t%5.3f/%5.3f\t%5.3f/%5.3f\t-\t\t-\t\t%5.3f/%5.3f\t%5.3f/%5.3f\t%5.3f/%5.3f\t%5.3f/%5.3f\t-\t\t-\t\t%5.3f/%5.3f\t%5.3f/%5.3f",
					   1.0+(ttwzpreds["JetDown"].ttwz_mm-ttwzpreds["Normal"].ttwz_mm)/ttwzpreds["Normal"].ttwz_mm,
					   1.0+(ttwzpreds["JetUp"].ttwz_mm-ttwzpreds["Normal"].ttwz_mm)/ttwzpreds["Normal"].ttwz_mm,
					   1.0+(ttwzpreds["JetDown"].ttwz_mm-ttwzpreds["Normal"].ttwz_mm)/ttwzpreds["Normal"].ttwz_mm,
					   1.0+(ttwzpreds["JetUp"].ttwz_mm-ttwzpreds["Normal"].ttwz_mm)/ttwzpreds["Normal"].ttwz_mm,
					   1.0+(ttwzpreds["JetDown"].wz_mm  -ttwzpreds["Normal"].wz_mm  )/ttwzpreds["Normal"].wz_mm,
					   1.0+(ttwzpreds["JetUp"].wz_mm  -ttwzpreds["Normal"].wz_mm  )/ttwzpreds["Normal"].wz_mm,
					   1.0+(ttwzpreds["JetDown"].rare_mm-ttwzpreds["Normal"].rare_mm)/ttwzpreds["Normal"].rare_mm,
					   1.0+(ttwzpreds["JetUp"].rare_mm-ttwzpreds["Normal"].rare_mm)/ttwzpreds["Normal"].rare_mm,
					   1.0+(ttwzpreds["JetDown"].ttwz_em-ttwzpreds["Normal"].ttwz_em)/ttwzpreds["Normal"].ttwz_em,
					   1.0+(ttwzpreds["JetUp"].ttwz_em-ttwzpreds["Normal"].ttwz_em)/ttwzpreds["Normal"].ttwz_em,
					   1.0+(ttwzpreds["JetDown"].ttwz_em-ttwzpreds["Normal"].ttwz_em)/ttwzpreds["Normal"].ttwz_em,
					   1.0+(ttwzpreds["JetUp"].ttwz_em-ttwzpreds["Normal"].ttwz_em)/ttwzpreds["Normal"].ttwz_em,
					   1.0+(ttwzpreds["JetDown"].wz_em  -ttwzpreds["Normal"].wz_em  )/ttwzpreds["Normal"].wz_em,
					   1.0+(ttwzpreds["JetUp"].wz_em  -ttwzpreds["Normal"].wz_em  )/ttwzpreds["Normal"].wz_em,
					   1.0+(ttwzpreds["JetDown"].rare_em-ttwzpreds["Normal"].rare_em)/ttwzpreds["Normal"].rare_em,
					   1.0+(ttwzpreds["JetUp"].rare_em-ttwzpreds["Normal"].rare_em)/ttwzpreds["Normal"].rare_em,
					   1.0+(ttwzpreds["JetDown"].ttwz_ee-ttwzpreds["Normal"].ttwz_ee)/ttwzpreds["Normal"].ttwz_ee,
					   1.0+(ttwzpreds["JetUp"].ttwz_ee-ttwzpreds["Normal"].ttwz_ee)/ttwzpreds["Normal"].ttwz_ee,
					   1.0+(ttwzpreds["JetDown"].ttwz_ee-ttwzpreds["Normal"].ttwz_ee)/ttwzpreds["Normal"].ttwz_ee,
					   1.0+(ttwzpreds["JetUp"].ttwz_ee-ttwzpreds["Normal"].ttwz_ee)/ttwzpreds["Normal"].ttwz_ee,
					   1.0+(ttwzpreds["JetDown"].wz_ee  -ttwzpreds["Normal"].wz_ee  )/ttwzpreds["Normal"].wz_ee,
					   1.0+(ttwzpreds["JetUp"].wz_ee  -ttwzpreds["Normal"].wz_ee  )/ttwzpreds["Normal"].wz_ee,
					   1.0+(ttwzpreds["JetDown"].rare_ee-ttwzpreds["Normal"].rare_ee)/ttwzpreds["Normal"].rare_ee,
					   1.0+(ttwzpreds["JetUp"].rare_ee-ttwzpreds["Normal"].rare_ee)/ttwzpreds["Normal"].rare_ee) << endl;
	fOUTSTREAM << Form("jer      lnN\t%5.3f\t\t%5.3f\t\t-\t\t-\t\t%5.3f\t\t%5.3f\t\t%5.3f\t\t%5.3f\t\t-\t\t-\t\t%5.3f\t\t%5.3f\t\t%5.3f\t\t%5.3f\t\t-\t\t-\t\t%5.3f\t\t%5.3f",
					   1.0+(ttwzpreds["JetSmear"].ttwz_mm-ttwzpreds["Normal"].ttwz_mm)/ttwzpreds["Normal"].ttwz_mm,
					   1.0+(ttwzpreds["JetSmear"].ttwz_mm-ttwzpreds["Normal"].ttwz_mm)/ttwzpreds["Normal"].ttwz_mm,
					   1.0+(ttwzpreds["JetSmear"].wz_mm  -ttwzpreds["Normal"].wz_mm  )/ttwzpreds["Normal"].wz_mm,
					   1.0+(ttwzpreds["JetSmear"].rare_mm-ttwzpreds["Normal"].rare_mm)/ttwzpreds["Normal"].rare_mm,
					   1.0+(ttwzpreds["JetSmear"].ttwz_em-ttwzpreds["Normal"].ttwz_em)/ttwzpreds["Normal"].ttwz_em,
					   1.0+(ttwzpreds["JetSmear"].ttwz_em-ttwzpreds["Normal"].ttwz_em)/ttwzpreds["Normal"].ttwz_em,
					   1.0+(ttwzpreds["JetSmear"].wz_em  -ttwzpreds["Normal"].wz_em  )/ttwzpreds["Normal"].wz_em,
					   1.0+(ttwzpreds["JetSmear"].rare_em-ttwzpreds["Normal"].rare_em)/ttwzpreds["Normal"].rare_em,
					   1.0+(ttwzpreds["JetSmear"].ttwz_ee-ttwzpreds["Normal"].ttwz_ee)/ttwzpreds["Normal"].ttwz_ee,
					   1.0+(ttwzpreds["JetSmear"].ttwz_ee-ttwzpreds["Normal"].ttwz_ee)/ttwzpreds["Normal"].ttwz_ee,
					   1.0+(ttwzpreds["JetSmear"].wz_ee  -ttwzpreds["Normal"].wz_ee  )/ttwzpreds["Normal"].wz_ee,
					   1.0+(ttwzpreds["JetSmear"].rare_ee-ttwzpreds["Normal"].rare_ee)/ttwzpreds["Normal"].rare_ee) << endl;
	fOUTSTREAM << Form("pu       lnN\t%5.3f\t\t%5.3f\t\t-\t\t-\t\t%5.3f\t\t%5.3f\t\t%5.3f\t\t%5.3f\t\t-\t\t-\t\t%5.3f\t\t%5.3f\t\t%5.3f\t\t%5.3f\t\t-\t\t-\t\t%5.3f\t\t%5.3f",
					   pu_syst, pu_syst, pu_syst, pu_syst, pu_syst, pu_syst, pu_syst, pu_syst, pu_syst, pu_syst, pu_syst, pu_syst) << endl;
	fOUTSTREAM << Form("matching lnN\t%5.3f/%5.3f\t%5.3f/%5.3f\t-\t\t-\t\t%5.3f/%5.3f\t%5.3f/%5.3f\t%5.3f/%5.3f\t%5.3f/%5.3f\t-\t\t-\t\t%5.3f/%5.3f\t%5.3f/%5.3f\t%5.3f/%5.3f\t%5.3f/%5.3f\t-\t\t-\t\t%5.3f/%5.3f\t%5.3f/%5.3f",
					   match_syst_up, match_syst_dn, match_syst_up, match_syst_dn, match_syst_up, match_syst_dn, match_syst_up, match_syst_dn, match_syst_up, match_syst_dn, match_syst_up, match_syst_dn,
					   match_syst_up, match_syst_dn, match_syst_up, match_syst_dn, match_syst_up, match_syst_dn, match_syst_up, match_syst_dn, match_syst_up, match_syst_dn, match_syst_up, match_syst_dn) << endl;
	fOUTSTREAM << Form("scale    lnN\t%5.3f/%5.3f\t%5.3f/%5.3f\t-\t\t-\t\t%5.3f/%5.3f\t%5.3f/%5.3f\t%5.3f/%5.3f\t%5.3f/%5.3f\t-\t\t-\t\t%5.3f/%5.3f\t%5.3f/%5.3f\t%5.3f/%5.3f\t%5.3f/%5.3f\t-\t\t-\t\t%5.3f/%5.3f\t%5.3f/%5.3f",
					   scale_syst_up, scale_syst_dn, scale_syst_up, scale_syst_dn, scale_syst_up, scale_syst_dn, scale_syst_up, scale_syst_dn, scale_syst_up, scale_syst_dn, scale_syst_up, scale_syst_dn,
					   scale_syst_up, scale_syst_dn, scale_syst_up, scale_syst_dn, scale_syst_up, scale_syst_dn, scale_syst_up, scale_syst_dn, scale_syst_up, scale_syst_dn, scale_syst_up, scale_syst_dn) << endl;
	fOUTSTREAM << Form("NLO      lnN\t%5.3f\t\t-\t\t-\t\t-\t\t-\t\t-\t\t%5.3f\t\t-\t\t-\t\t-\t\t-\t\t-\t\t%5.3f\t\t-\t\t-\t\t-\t\t-\t\t-",
					   1.130, 1.130, 1.130) << endl; // old numbers for the moment
	fOUTSTREAM << endl;
	fOUTSTREAM.close();
	
	// Produce systematics plots
	TH1D *h_ttwz_nom = new TH1D("h_ttwz",    "Nominal ttWZ",         4, 0., 4.);
	TH1D *h_ttwz_ju  = new TH1D("h_ttwz_ju", "ttWZ jets up",         4, 0., 4.);
	TH1D *h_ttwz_jd  = new TH1D("h_ttwz_jd", "ttWZ jets dn",         4, 0., 4.);
	TH1D *h_ttwz_js  = new TH1D("h_ttwz_js", "ttWZ jets sm",         4, 0., 4.);
	TH1D *h_ttwz_bu  = new TH1D("h_ttwz_bu", "ttWZ b-tags up",       4, 0., 4.);
	TH1D *h_ttwz_bd  = new TH1D("h_ttwz_bd", "ttWZ b-tags dn",       4, 0., 4.);
	TH1D *h_ttwz_lu  = new TH1D("h_ttwz_lu", "ttWZ lepton up",       4, 0., 4.);
	TH1D *h_ttwz_ld  = new TH1D("h_ttwz_ld", "ttWZ lepton dn",       4, 0., 4.);
	TH1D *h_bg_nom   = new TH1D("h_bg",      "Nominal background",   4, 0., 4.);
	TH1D *h_bg_ju    = new TH1D("h_bg_ju",   "background jets up",   4, 0., 4.);
	TH1D *h_bg_jd    = new TH1D("h_bg_jd",   "background jets dn",   4, 0., 4.);
	TH1D *h_bg_js    = new TH1D("h_bg_js",   "background jets sm",   4, 0., 4.);
	TH1D *h_bg_bu    = new TH1D("h_bg_bu",   "background b-tags up", 4, 0., 4.);
	TH1D *h_bg_bd    = new TH1D("h_bg_bd",   "background b-tags dn", 4, 0., 4.);
	TH1D *h_bg_lu    = new TH1D("h_bg_lu",   "background lepton up", 4, 0., 4.);
	TH1D *h_bg_ld    = new TH1D("h_bg_ld",   "background lepton dn", 4, 0., 4.);
	vector<TH1D*> histos;
	histos.push_back(h_ttwz_nom);
	histos.push_back(h_ttwz_ju);
	histos.push_back(h_ttwz_jd);
	histos.push_back(h_ttwz_js);
	histos.push_back(h_ttwz_bu);
	histos.push_back(h_ttwz_bd);
	histos.push_back(h_ttwz_lu);
	histos.push_back(h_ttwz_ld);
	histos.push_back(h_bg_nom);
	histos.push_back(h_bg_ju);
	histos.push_back(h_bg_jd);
	histos.push_back(h_bg_js);
	histos.push_back(h_bg_bu);
	histos.push_back(h_bg_bd);
	histos.push_back(h_bg_lu);
	histos.push_back(h_bg_ld);
	
	h_ttwz_nom->SetBinContent(1, ttwzpreds["Normal"].ttwz_ee);
	h_ttwz_nom->SetBinContent(2, ttwzpreds["Normal"].ttwz_mm);
	h_ttwz_nom->SetBinContent(3, ttwzpreds["Normal"].ttwz_em);
	h_ttwz_nom->SetBinContent(4, ttwzpreds["Normal"].ttwz);
	
	h_ttwz_ju ->SetBinContent(1, ttwzpreds["JetUp"].ttwz_ee);
	h_ttwz_ju ->SetBinContent(2, ttwzpreds["JetUp"].ttwz_mm);
	h_ttwz_ju ->SetBinContent(3, ttwzpreds["JetUp"].ttwz_em);
	h_ttwz_ju ->SetBinContent(4, ttwzpreds["JetUp"].ttwz);
	h_ttwz_jd ->SetBinContent(1, ttwzpreds["JetDown"].ttwz_ee);
	h_ttwz_jd ->SetBinContent(2, ttwzpreds["JetDown"].ttwz_mm);
	h_ttwz_jd ->SetBinContent(3, ttwzpreds["JetDown"].ttwz_em);
	h_ttwz_jd ->SetBinContent(4, ttwzpreds["JetDown"].ttwz);
	
	h_ttwz_js ->SetBinContent(1, ttwzpreds["JetSmear"].ttwz_ee);
	h_ttwz_js ->SetBinContent(2, ttwzpreds["JetSmear"].ttwz_mm);
	h_ttwz_js ->SetBinContent(3, ttwzpreds["JetSmear"].ttwz_em);
	h_ttwz_js ->SetBinContent(4, ttwzpreds["JetSmear"].ttwz);
	
	h_ttwz_bu ->SetBinContent(1, ttwzpreds["BUp"].ttwz_ee);
	h_ttwz_bu ->SetBinContent(2, ttwzpreds["BUp"].ttwz_mm);
	h_ttwz_bu ->SetBinContent(3, ttwzpreds["BUp"].ttwz_em);
	h_ttwz_bu ->SetBinContent(4, ttwzpreds["BUp"].ttwz);
	h_ttwz_bd ->SetBinContent(1, ttwzpreds["BDown"].ttwz_ee);
	h_ttwz_bd ->SetBinContent(2, ttwzpreds["BDown"].ttwz_mm);
	h_ttwz_bd ->SetBinContent(3, ttwzpreds["BDown"].ttwz_em);
	h_ttwz_bd ->SetBinContent(4, ttwzpreds["BDown"].ttwz);
	
	h_ttwz_lu ->SetBinContent(1, ttwzpreds["LepUp"].ttwz_ee);
	h_ttwz_lu ->SetBinContent(2, ttwzpreds["LepUp"].ttwz_mm);
	h_ttwz_lu ->SetBinContent(3, ttwzpreds["LepUp"].ttwz_em);
	h_ttwz_lu ->SetBinContent(4, ttwzpreds["LepUp"].ttwz);
	h_ttwz_ld ->SetBinContent(1, ttwzpreds["LepDown"].ttwz_ee);
	h_ttwz_ld ->SetBinContent(2, ttwzpreds["LepDown"].ttwz_mm);
	h_ttwz_ld ->SetBinContent(3, ttwzpreds["LepDown"].ttwz_em);
	h_ttwz_ld ->SetBinContent(4, ttwzpreds["LepDown"].ttwz);
	
	h_bg_nom->SetBinContent(1, ttwzpreds["Normal"].wz_ee+ttwzpreds["Normal"].rare_ee);
	h_bg_nom->SetBinContent(2, ttwzpreds["Normal"].wz_mm+ttwzpreds["Normal"].rare_mm);
	h_bg_nom->SetBinContent(3, ttwzpreds["Normal"].wz_em+ttwzpreds["Normal"].rare_em);
	h_bg_nom->SetBinContent(4, ttwzpreds["Normal"].wz   +ttwzpreds["Normal"].rare   );
	
	h_bg_ju ->SetBinContent(1, ttwzpreds["JetUp"].wz_ee+ttwzpreds["JetUp"].rare_ee);
	h_bg_ju ->SetBinContent(2, ttwzpreds["JetUp"].wz_mm+ttwzpreds["JetUp"].rare_mm);
	h_bg_ju ->SetBinContent(3, ttwzpreds["JetUp"].wz_em+ttwzpreds["JetUp"].rare_em);
	h_bg_ju ->SetBinContent(4, ttwzpreds["JetUp"].wz   +ttwzpreds["JetUp"].rare   );
	h_bg_jd ->SetBinContent(1, ttwzpreds["JetDown"].wz_ee+ttwzpreds["JetDown"].rare_ee);
	h_bg_jd ->SetBinContent(2, ttwzpreds["JetDown"].wz_mm+ttwzpreds["JetDown"].rare_mm);
	h_bg_jd ->SetBinContent(3, ttwzpreds["JetDown"].wz_em+ttwzpreds["JetDown"].rare_em);
	h_bg_jd ->SetBinContent(4, ttwzpreds["JetDown"].wz   +ttwzpreds["JetDown"].rare   );
	
	h_bg_js ->SetBinContent(1, ttwzpreds["JetSmear"].wz_ee+ttwzpreds["JetSmear"].rare_ee);
	h_bg_js ->SetBinContent(2, ttwzpreds["JetSmear"].wz_mm+ttwzpreds["JetSmear"].rare_mm);
	h_bg_js ->SetBinContent(3, ttwzpreds["JetSmear"].wz_em+ttwzpreds["JetSmear"].rare_em);
	h_bg_js ->SetBinContent(4, ttwzpreds["JetSmear"].wz   +ttwzpreds["JetSmear"].rare   );
	
	h_bg_bu ->SetBinContent(1, ttwzpreds["BUp"].wz_ee+ttwzpreds["BUp"].rare_ee);
	h_bg_bu ->SetBinContent(2, ttwzpreds["BUp"].wz_mm+ttwzpreds["BUp"].rare_mm);
	h_bg_bu ->SetBinContent(3, ttwzpreds["BUp"].wz_em+ttwzpreds["BUp"].rare_em);
	h_bg_bu ->SetBinContent(4, ttwzpreds["BUp"].wz   +ttwzpreds["BUp"].rare   );
	h_bg_bd ->SetBinContent(1, ttwzpreds["BDown"].wz_ee+ttwzpreds["BDown"].rare_ee);
	h_bg_bd ->SetBinContent(2, ttwzpreds["BDown"].wz_mm+ttwzpreds["BDown"].rare_mm);
	h_bg_bd ->SetBinContent(3, ttwzpreds["BDown"].wz_em+ttwzpreds["BDown"].rare_em);
	h_bg_bd ->SetBinContent(4, ttwzpreds["BDown"].wz   +ttwzpreds["BDown"].rare   );
	
	h_bg_lu ->SetBinContent(1, ttwzpreds["LepUp"].wz_ee+ttwzpreds["LepUp"].rare_ee);
	h_bg_lu ->SetBinContent(2, ttwzpreds["LepUp"].wz_mm+ttwzpreds["LepUp"].rare_mm);
	h_bg_lu ->SetBinContent(3, ttwzpreds["LepUp"].wz_em+ttwzpreds["LepUp"].rare_em);
	h_bg_lu ->SetBinContent(4, ttwzpreds["LepUp"].wz   +ttwzpreds["LepUp"].rare   );
	h_bg_ld ->SetBinContent(1, ttwzpreds["LepDown"].wz_ee+ttwzpreds["LepDown"].rare_ee);
	h_bg_ld ->SetBinContent(2, ttwzpreds["LepDown"].wz_mm+ttwzpreds["LepDown"].rare_mm);
	h_bg_ld ->SetBinContent(3, ttwzpreds["LepDown"].wz_em+ttwzpreds["LepDown"].rare_em);
	h_bg_ld ->SetBinContent(4, ttwzpreds["LepDown"].wz   +ttwzpreds["LepDown"].rare   );
	
	h_ttwz_nom->SetLineColor(kBlack);
	h_bg_nom  ->SetLineColor(kBlack);
	
	h_ttwz_js->SetLineColor(kGreen+1);
	h_ttwz_ju->SetLineColor(kGreen+1);
	h_ttwz_bu->SetLineColor(kGreen+1);
	h_ttwz_lu->SetLineColor(kGreen+1);
	h_bg_js  ->SetLineColor(kGreen+1);
	h_bg_ju  ->SetLineColor(kGreen+1);
	h_bg_bu  ->SetLineColor(kGreen+1);
	h_bg_lu  ->SetLineColor(kGreen+1);
	
	h_ttwz_jd->SetLineColor(kRed-3);
	h_ttwz_bd->SetLineColor(kRed-3);
	h_ttwz_ld->SetLineColor(kRed-3);
	h_bg_jd  ->SetLineColor(kRed-3);
	h_bg_bd  ->SetLineColor(kRed-3);
	h_bg_ld  ->SetLineColor(kRed-3);
	
	int cnt = 0;
	for(vector<TH1D*>::iterator it = histos.begin(); it < histos.end(); ++it){
		(*it)->GetXaxis()->SetBinLabel(1, "ee");
		(*it)->GetXaxis()->SetBinLabel(2, "#mu#mu");
		(*it)->GetXaxis()->SetBinLabel(3, "e#mu");
		(*it)->GetXaxis()->SetBinLabel(4, "Total");
		(*it)->GetXaxis()->SetLabelSize(0.06);
		if(cnt<8)(*it)->GetYaxis()->SetTitle("Signal Events");
		if(cnt>7)(*it)->GetYaxis()->SetTitle("Background Events");
		(*it)->GetYaxis()->SetTitleOffset(1.15);
		(*it)->SetMinimum(0);
		if(cnt<8)(*it)->SetMaximum(10);
		if(cnt>7)(*it)->SetMaximum(4);
		(*it)->SetLineWidth(2);
		(*it)->SetFillStyle(0);
		cnt++;
	}
	
	makeSystPlot("Syst_Sig_JES"  + chargeString, "JES",          h_ttwz_nom, h_ttwz_ju, h_ttwz_jd);
	makeSystPlot("Syst_Sig_bTag" + chargeString, "b Tag",        h_ttwz_nom, h_ttwz_bu, h_ttwz_bd);
	makeSystPlot("Syst_Sig_Lept" + chargeString, "Lepton Scale", h_ttwz_nom, h_ttwz_lu, h_ttwz_ld);
	makeSystPlot("Syst_Sig_JER"  + chargeString, "JER",          h_ttwz_nom, h_ttwz_js);
	makeSystPlot("Syst_Bg_JES"   + chargeString, "JES",          h_bg_nom, h_bg_ju, h_bg_jd);
	makeSystPlot("Syst_Bg_bTag"  + chargeString, "b Tag",        h_bg_nom, h_bg_bu, h_bg_bd);
	makeSystPlot("Syst_Bg_Lept"  + chargeString, "Lepton Scale", h_bg_nom, h_bg_lu, h_bg_ld);
	makeSystPlot("Syst_Bg_JER"   + chargeString, "JER",          h_bg_nom, h_bg_js);
	
}
TTWZPrediction SSDLPlotter::makePredictionSignalEvents(float minHT, float maxHT, float minMET, float maxMET, int minNjets, int minNbjetsL, int minNbjetsM, float minPt1, float minPt2, int chVeto, bool ttw, int systflag){
	fOutputSubDir = "IntPredictionsSigEventTree/";
	TString jvString = "";
	TString sysString = "";
	TString chargeString = "";

	for( gsystIt = gSystematics.begin(); gsystIt != gSystematics.end(); ++gsystIt) {
		if (gsystIt->second == systflag) { sysString = gsystIt->first; break; }
	}
	if (chVeto == -1) chargeString = "_mm";
	if (chVeto == +1) chargeString = "_pp";

	if (maxHT < 20.) jvString = "JV";
	ofstream OUT(fOutputDir+fOutputSubDir+Form("DataPred_customRegion_HT%.0f"+jvString+"MET%.0fNJ%.0iNbjL%.0iNbjM%.0iPT1%.0fPT2%.0f_"+sysString+chargeString+".txt", minHT, minMET, minNjets, minNbjetsL, minNbjetsM, minPt1, minPt2), ios::trunc);

	TLatex *lat = new TLatex();
	lat->SetNDC(kTRUE);
	lat->SetTextColor(kBlack);
	lat->SetTextSize(0.04);

	const float RareESyst  = 0.5;
	const float RareESyst2 = RareESyst*RareESyst;
	
	const float FakeESyst  = 0.5;
	const float FakeESyst2 = FakeESyst*FakeESyst;

	const float WZESyst  = 0.15;
	const float WZESyst2 = WZESyst*WZESyst;
	
	const float TTZESyst  = 0.5;
	const float TTZESyst2 = TTZESyst*TTZESyst;

	OUT << "/////////////////////////////////////////////////////////////////////////////" << endl;
	OUT << " Producing predictions " ;
	OUT << "  scaling MC to " << fLumiNorm << " /pb" << endl << endl;
	OUT << "-----------------------------------------------------------------------------" << endl;
	OUT << " These are the cuts: " << endl;
	OUT << Form("  minHT:    %4.0f  || maxHT:     %4.0f", minHT   , maxHT    ) << endl;
	OUT << Form("  minMET:   %4.0f  || maxMET:    %4.0f", minMET  , maxMET   ) << endl;
	OUT << Form("  minNjets:   %2i  ", minNjets) << endl;
	OUT << Form("  minNbjetsL: %2i  ", minNbjetsL) << endl;
	OUT << Form("  minNbjetsM: %2i  ", minNbjetsM) << endl;
	OUT << Form("  minpT1:     %2.0f  ", minPt1) << endl;
	OUT << Form("  minpT2:     %2.0f  ", minPt2) << endl;
	OUT << "-----------------------------------------------------------------------------" << endl;

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

	// only take half the events for ++/--
	float chargeFactor = chVeto ? 0.5:1.;

	// rare SM yields
	float nt2_rare_mc_mm(0.),    nt2_rare_mc_em(0.),    nt2_rare_mc_ee(0.);
	float nt2_rare_mc_mm_e2(0.), nt2_rare_mc_em_e2(0.), nt2_rare_mc_ee_e2(0.);

	float nt2_wz_mc_mm(0.),    nt2_wz_mc_em(0.),    nt2_wz_mc_ee(0.);
	float nt2_wz_mc_mm_e2(0.), nt2_wz_mc_em_e2(0.), nt2_wz_mc_ee_e2(0.);

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	TFile *pFile = TFile::Open(fOutputFileName);
	TTree *sigtree; getObjectSafe(pFile, "SigEvents", sigtree);
	
	string *sname = 0;
	int flag;
	int   SType, Flavor, TLCat, NJ, NbJ, NbJmed;
	float puweight, pT1, pT2, HT, MET, MT2, SLumi;
	float eta1, eta2, mll;
	int   event, run;
	int charge;
	int passZVeto;

	sigtree->SetBranchAddress("SystFlag", &flag);
	sigtree->SetBranchAddress("Event",    &event);
	sigtree->SetBranchAddress("Run",      &run);
	sigtree->SetBranchAddress("SName",    &sname);
	sigtree->SetBranchAddress("SType",    &SType);
	sigtree->SetBranchAddress("PUWeight", &puweight);
	sigtree->SetBranchAddress("SLumi",    &SLumi);
	sigtree->SetBranchAddress("Flavor",   &Flavor);
	sigtree->SetBranchAddress("Charge",   &charge);
	sigtree->SetBranchAddress("pT1",      &pT1);
	sigtree->SetBranchAddress("pT2",      &pT2);
	sigtree->SetBranchAddress("eta1",     &eta1);
	sigtree->SetBranchAddress("eta2",     &eta2);
	sigtree->SetBranchAddress("TLCat",    &TLCat);
	sigtree->SetBranchAddress("HT",       &HT);
	sigtree->SetBranchAddress("MET",      &MET);
	sigtree->SetBranchAddress("MT2",      &MT2);
	sigtree->SetBranchAddress("NJ",       &NJ);
	sigtree->SetBranchAddress("NbJ",      &NbJ);
	sigtree->SetBranchAddress("NbJmed",   &NbJmed);
	sigtree->SetBranchAddress("Mll",      &mll);
	sigtree->SetBranchAddress("PassZVeto",&passZVeto);

	FakeRatios *FR = new FakeRatios();

	std::map< std::string, float > rareMapMM;
	std::map< std::string, float > rareMapEM;
	std::map< std::string, float > rareMapEE;

	std::map< std::string, int > rareMapMM_npass;
	std::map< std::string, int > rareMapEM_npass;
	std::map< std::string, int > rareMapEE_npass;

	for(size_t i = 0; i < fMCRareSM.size(); ++i){
		Sample *S = fSamples[fMCRareSM[i]];
		std::string name = (string) S->sname;
		rareMapMM[name] = 0.; rareMapMM_npass[name] = 0;
		rareMapEM[name] = 0.; rareMapEM_npass[name] = 0;
		rareMapEE[name] = 0.; rareMapEE_npass[name] = 0;
	}

	float trigScale[3] = {gMMTrigScale, gEMTrigScale, gEETrigScale};

	for( int i = 0; i < sigtree->GetEntries(); i++ ){
		sigtree->GetEntry(i);
		
		if( flag != systflag ) continue;
		
		if ( mll < 8.) continue;
		if ( HT  < minHT  || HT  > maxHT)  continue;
		if ( MET < minMET || MET > maxMET) continue;
		if ( NJ  < minNjets)      continue;
		if ( NbJ < minNbjetsL)    continue;
		if ( NbJmed < minNbjetsM) continue;

		gChannel chan = gChannel(Flavor);
		if(chan == ElMu || Flavor == 4){
			if(pT1 > pT2){
				if(pT1 < minPt1) continue;
				if(pT2 < minPt2) continue;
			}
			if(pT1 < pT2){
				if(pT1 < minPt2) continue;
				if(pT2 < minPt1) continue;
			}
		}
		else{
			if(pT1 < minPt1) continue;
			if(pT2 < minPt2) continue;
		}

		// GET ALL DATA EVENTS
		if(SType < 3) {             // 0,1,2 are DoubleMu, DoubleEle, MuEG
			if (Flavor < 3) {
				if (gApplyZVeto && passZVeto == 0)  continue;
				if (chVeto && charge != chVeto ) continue;
				Sample *S = fSampleMap[TString(*sname)];

				float npp(0.) , npf(0.) , nfp(0.) , nff(0.);
				float f1(0.)  , f2(0.)  , p1(0.)  , p2(0.);
				f1 = getFRatio(chan, pT1, eta1, S->datamc);
				f2 = getFRatio(chan, pT2, eta2, S->datamc);
				p1 = getPRatio(chan, pT1, S->datamc);
				p2 = getPRatio(chan, pT2, S->datamc);
				if(chan == ElMu){
					f1 = getFRatio(Muon, pT1, eta1, S->datamc);
					f2 = getFRatio(Elec, pT2, eta2, S->datamc);
					p1 = getPRatio(Muon, pT1, S->datamc);
					p2 = getPRatio(Elec, pT2, S->datamc);
				}
				// Get the weights (don't depend on event selection)
				npp = FR->getWpp(FakeRatios::gTLCat(TLCat), f1, f2, p1, p2);
				npf = FR->getWpf(FakeRatios::gTLCat(TLCat), f1, f2, p1, p2);
				nfp = FR->getWfp(FakeRatios::gTLCat(TLCat), f1, f2, p1, p2);
				nff = FR->getWff(FakeRatios::gTLCat(TLCat), f1, f2, p1, p2);			

				if (Flavor == 0) {      // MUMU
					npp_mm += npp;
					npf_mm += (npf+nfp);
					nff_mm += nff;
					if (TLCat == 0)               nt2_mm++;
				 	if (TLCat == 1 || TLCat == 2) nt10_mm++;
				 	if (TLCat == 3)               nt0_mm++;
				}
				if (Flavor == 1) {       // E-MU
					npp_em += npp;
					npf_em += npf;
					nfp_em += nfp;
					nff_em += nff;
					if (TLCat == 0) nt2_em++;
				 	if (TLCat == 1) nt10_em++;
				 	if (TLCat == 2) nt01_em++;
				 	if (TLCat == 3) nt0_em++;
				}
				if (Flavor == 2) {       // E-E
					npp_ee += npp;
					npf_ee += (nfp+npf);
					nff_ee += nff;
					if (TLCat == 0)               nt2_ee++;
				 	if (TLCat == 1 || TLCat == 2) nt10_ee++;
				 	if (TLCat == 3)               nt0_ee++;
				}
			}

			if(Flavor == 4) {       // E-MU OS
				if (TLCat == 0) nt2_em_BB_os+= chargeFactor;
				if (TLCat == 1) nt2_em_EE_os+= chargeFactor;
			}
			if(Flavor == 5) {       // E-E OS
				if (TLCat == 0)               nt2_ee_BB_os+= chargeFactor;
				if (TLCat == 1 || TLCat == 2) nt2_ee_EB_os+= chargeFactor;
				if (TLCat == 3)               nt2_ee_EE_os+= chargeFactor;
			}
		} // end data events


		// GET RARE MC EVENTS
		if (SType == 15 && TLCat == 0) { // tight-tight rare MC events
			if (*sname == "WWTo2L2Nu") continue;
			if (Flavor > 2) continue;
			if (chVeto && charge != chVeto ) continue;
			// make sure here to get the lumi from the sample and not from the Slumi variable. Doesn't work otherwise
			Sample *S = fSampleMap[*sname];
			float scale = fLumiNorm / S->getLumi();
			float tmp_nt2_rare_mc    = puweight*trigScale[Flavor]*scale;
			// float tmp_nt2_rare_mc_e2 = tmp_nt2_rare_mc*tmp_nt2_rare_mc;
			if (Flavor == 0) {
				rareMapMM[*sname] += tmp_nt2_rare_mc;
				rareMapMM_npass[*sname] ++;
			}
			if (Flavor == 1) {
				rareMapEM[*sname] += tmp_nt2_rare_mc;
				rareMapEM_npass[*sname] ++;
			}
			if (Flavor == 2) {
				rareMapEE[*sname] += tmp_nt2_rare_mc;
				rareMapEE_npass[*sname] ++;
			}
		} // end rare mc events
		
	}
	// float nt2_rare_mc_mm_e1 = sqrt(nt2_rare_mc_mm_e2);
	// float nt2_rare_mc_em_e1 = sqrt(nt2_rare_mc_em_e2);
	// float nt2_rare_mc_ee_e1 = sqrt(nt2_rare_mc_ee_e2);

	// float nt2_wz_mc_mm_e1 = sqrt(nt2_wz_mc_mm_e2);
	// float nt2_wz_mc_em_e1 = sqrt(nt2_wz_mc_em_e2);
	// float nt2_wz_mc_ee_e1 = sqrt(nt2_wz_mc_ee_e2);

	///////////////////////////////////////////////////////////////////////////////////
	// PRINTOUT ///////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////
	OUT << "---------------------------------------------------------------------------------------------------------" << endl;
	OUT << "         RATIOS  ||     Mu-fRatio      |     Mu-pRatio      ||     El-fRatio      |     El-pRatio      ||" << endl;
	OUT << "---------------------------------------------------------------------------------------------------------" << endl;
	OUT << setw(16) << "  data stat only ||";
	OUT << setw(7)  << setprecision(3) << mufratio_data  << " ± " << setw(7) << setprecision(3) << mufratio_data_e  << " |";
	OUT << setw(7)  << setprecision(3) << mupratio_data  << " ± " << setw(7) << setprecision(3) << mupratio_data_e  << " ||";
	OUT << setw(7)  << setprecision(3) << elfratio_data  << " ± " << setw(7) << setprecision(3) << elfratio_data_e  << " |";
	OUT << setw(7)  << setprecision(3) << elpratio_data  << " ± " << setw(7) << setprecision(3) << elpratio_data_e  << " ||";
	OUT << endl;
	OUT << "---------------------------------------------------------------------------------------------------------" << endl;

	OUT << "---------------------------------------------------------------------------------------------------------" << endl << endl;
	OUT << "-------------------------------------------------------------------------------------------------------------" << endl;
	OUT << "                 |           Mu/Mu          |                E/Mu               |           E/E            ||" << endl;
	OUT << "         YIELDS  |   Ntt  |   Nt1  |   Nll  |   Ntt  |   Ntl  |   Nlt  |   Nll  |   Ntt  |   Nt1  |   Nll  ||" << endl;
	OUT << "-------------------------------------------------------------------------------------------------------------" << endl;
	OUT << Form("%16s & %6.0f & %6.0f & %6.0f & %6.0f & %6.0f & %6.0f & %6.0f & %6.0f & %6.0f & %6.0f \\\\\n", "Data",
	nt2_mm, nt10_mm, nt0_mm, nt2_em, nt10_em, nt01_em, nt0_em, nt2_ee, nt10_ee, nt0_ee);

	///////////////////////////////////////////////////////////////////////////////////
	// PREDICTIONS ////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////
	// FakeRatios *FR = new FakeRatios();
	FR->setNToyMCs(100);
	FR->setAddESyst(0.5);

	FR->setMFRatio(mufratio_data, mufratio_data_e); // set error to pure statistical of ratio
	FR->setEFRatio(elfratio_data, elfratio_data_e);
	FR->setMPRatio(mupratio_data, mupratio_data_e);
	FR->setEPRatio(elpratio_data, elpratio_data_e);

	FR->setMMNtl(nt2_mm, nt10_mm, nt0_mm);
	FR->setEENtl(nt2_ee, nt10_ee, nt0_ee);
	FR->setEMNtl(nt2_em, nt10_em, nt01_em, nt0_em);

	// Event-by-event differential ratios:
	float nF_mm = npf_mm + nff_mm;
	float nF_em = npf_em + nfp_em + nff_em;
	float nF_ee = npf_ee + nff_ee;
	float nSF   = npf_mm + npf_em + nfp_em + npf_ee;
	float nDF   = nff_mm + nff_em + nff_ee;
	float nF    = nF_mm + nF_em + nF_ee;

	OUT << endl << endl << endl;
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
	float fbb(0.),fee(0.),feb(0.);
	float fbbE(0.),feeE(0.),febE(0.);
	float fbb_mc(0.),fee_mc(0.),feb_mc(0.);
	float fbbE_mc(0.),feeE_mc(0.),febE_mc(0.);
	
	calculateChMisIdProb(fEGData, BB, fbb, fbbE);
	calculateChMisIdProb(fEGData, EB, feb, febE);
	calculateChMisIdProb(fEGData, EE, fee, feeE);
	
	calculateChMisIdProb(fMCBG, BB, fbb_mc, fbbE_mc);
	calculateChMisIdProb(fMCBG, EB, feb_mc, febE_mc);
	calculateChMisIdProb(fMCBG, EE, fee_mc, feeE_mc);
	
	cout << "this is the number of OS events in EE, BB, and EB from the SigEventsTree: " << nt2_ee_EE_os << " " << nt2_ee_BB_os << " " << nt2_ee_EB_os << endl;
	cout << "this is the charge mis ID probabilities: fee, fbb, feb: " << fee << " " << fbb << " " << feb << endl;
	// Simple error propagation assuming error on number of events is sqrt(N)
	nt2_ee_chmid    = 2*fbb*nt2_ee_BB_os                           + 2*fee*nt2_ee_EE_os                      + 2*feb*nt2_ee_EB_os;
	nt2_ee_chmid_e1 = sqrt( 4*fbb*fbb*FR->getEStat2(nt2_ee_BB_os)  + 4*fee*fee*FR->getEStat2(nt2_ee_EE_os)   + 4*feb*feb*FR->getEStat2(nt2_ee_EB_os) ); // stat only
	nt2_ee_chmid_e2 = sqrt( 4*fbbE*fbbE*nt2_ee_BB_os*nt2_ee_BB_os  + 4*feeE*feeE*nt2_ee_EE_os*nt2_ee_EE_os   + 4*febE*febE*nt2_ee_EB_os*nt2_ee_EB_os ); // syst only
	
	nt2_em_chmid    = fbb*nt2_em_BB_os + fee*nt2_em_EE_os;
	nt2_em_chmid_e1 = sqrt( fbb*fbb*FR->getEStat2(nt2_em_BB_os) + fee*fee*FR->getEStat2(nt2_em_EE_os) );
	nt2_em_chmid_e2 = fbbE*fbbE*sqrt( nt2_em_BB_os*nt2_em_BB_os + feeE*feeE*nt2_em_EE_os*nt2_em_EE_os );


	///////////////////////////////////////////////////////////////////////////////////
	// PRINTOUT ///////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////
	OUT << "----------------------------------------------------------------------------------" << endl;
	OUT << "       E-ChMisID  ||       BB           |       BE          |       EE          ||" << endl;
	OUT << "----------------------------------------------------------------------------------" << endl;
	OUT << "                  ||";
	OUT << setw(7)  << setprecision(2) << fbb << " ± " << setw(7) << setprecision(3) << feb  << setw(7) << setprecision(3) << fee  << " ||";
	OUT << endl;
	OUT << "--------------------------------------------------------------" << endl << endl;

	OUT << endl;
	OUT << "---------------------------------------------------------------------------------------------------------" << endl << endl;
	OUT << "-------------------------------------------------------------------------------------------------------------" << endl;
	OUT << "                 |           Mu/Mu          |                E/Mu               |           E/E            ||" << endl;
	OUT << "         YIELDS  |   Ntt  |   Nt1  |   Nll  |   Ntt  |   Ntl  |   Nlt  |   Nll  |   Ntt  |   Nt1  |   Nll  ||" << endl;
	OUT << "-------------------------------------------------------------------------------------------------------------" << endl;
	OUT << Form("%16s & %6.0f & %6.0f & %6.0f & %6.0f & %6.0f & %6.0f & %6.0f & %6.0f & %6.0f & %6.0f \\\\\n", "Data",
	nt2_mm, nt10_mm, nt0_mm, nt2_em, nt10_em, nt01_em, nt0_em, nt2_ee, nt10_ee, nt0_ee);

	OUT << "/////////////////////////////////////////////////////////////////////////////" << endl;
	OUT << endl;
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
	OUT << "----------------------------------------------------------------------------------------------" << endl;

	float nt2_ttw_mc_mm (0.); float nt2_ttw_mc_mm_e2(0.);
	float nt2_ttw_mc_em (0.); float nt2_ttw_mc_em_e2(0.);
	float nt2_ttw_mc_ee (0.); float nt2_ttw_mc_ee_e2(0.);

	int nt2_ttw_Nmc_mm (0);
	int nt2_ttw_Nmc_em (0);
	int nt2_ttw_Nmc_ee (0);

	float nt2_ttz_mc_mm (0.); float nt2_ttz_mc_mm_e2(0.);
	float nt2_ttz_mc_em (0.); float nt2_ttz_mc_em_e2(0.);
	float nt2_ttz_mc_ee (0.); float nt2_ttz_mc_ee_e2(0.);

	int nt2_ttz_Nmc_mm (0);
	int nt2_ttz_Nmc_em (0);
	int nt2_ttz_Nmc_ee (0);

	// DUMP ALL NUMBERS FOR RARE SAMPLES
	std::map<std::string , float >::const_iterator it = rareMapMM.begin();
	for ( ; it != rareMapMM.end() ; it++){
		Sample *S = fSampleMap[it->first];
		float weight = fLumiNorm / S->getLumi();

		float MM_yiel = rareMapMM[it->first];
		float EM_yiel = rareMapEM[it->first];
		float EE_yiel = rareMapEE[it->first];

		float MM_stat = weight*trigScale[0]*(S->getError(rareMapMM_npass[it->first]));
		float EM_stat = weight*trigScale[1]*(S->getError(rareMapEM_npass[it->first]));
		float EE_stat = weight*trigScale[2]*(S->getError(rareMapEE_npass[it->first]));
		// float MM_yiel = rareMapMM.find(it->first) != rareMapMM.end()? rareMapMM[it->first]:0.;
		// float EM_yiel = rareMapEM.find(it->first) != rareMapEM.end()? rareMapEM[it->first]:0.;
		// float EE_yiel = rareMapEE.find(it->first) != rareMapEE.end()? rareMapEE[it->first]:0.;
		// float MM_stat = rareMapMM_stat.find(it->first) != rareMapMM_stat.end() ? weight*S->getError(rareMapMM_npass[it->first]):0.;
		// float EM_stat = rareMapEM_stat.find(it->first) != rareMapEM_stat.end() ? weight*S->getError(rareMapEM_npass[it->first]):0.;
		// float EE_stat = rareMapEE_stat.find(it->first) != rareMapEE_stat.end() ? weight*S->getError(rareMapEE_npass[it->first]):0.;


		if (it->first == "WZTo3LNu") {
			nt2_wz_mc_mm = MM_yiel; nt2_wz_mc_mm_e2 = MM_stat*MM_stat;
			nt2_wz_mc_em = EM_yiel; nt2_wz_mc_em_e2 = EM_stat*EM_stat;
			nt2_wz_mc_ee = EE_yiel; nt2_wz_mc_ee_e2 = EE_stat*EE_stat;
			continue;
		}
		else {
			if (ttw && (it->first == "TTbarW") ){
				nt2_ttw_mc_mm += MM_yiel; nt2_ttw_mc_mm_e2 += MM_stat*MM_stat;
				nt2_ttw_mc_em += EM_yiel; nt2_ttw_mc_em_e2 += EM_stat*EM_stat;
				nt2_ttw_mc_ee += EE_yiel; nt2_ttw_mc_ee_e2 += EE_stat*EE_stat;
				continue;
			}
			else if (ttw && (it->first == "TTbarZ") ){
				nt2_ttz_mc_mm += MM_yiel; nt2_ttz_mc_mm_e2 += MM_stat*MM_stat;
				nt2_ttz_mc_em += EM_yiel; nt2_ttz_mc_em_e2 += EM_stat*EM_stat;
				nt2_ttz_mc_ee += EE_yiel; nt2_ttz_mc_ee_e2 += EE_stat*EE_stat;
				continue;
			}
			else {
				nt2_rare_mc_mm += MM_yiel; nt2_rare_mc_mm_e2 += MM_stat*MM_stat;
				nt2_rare_mc_em += EM_yiel; nt2_rare_mc_em_e2 += EM_stat*EM_stat;
				nt2_rare_mc_ee += EE_yiel; nt2_rare_mc_ee_e2 += EE_stat*EE_stat;
			}
		}


		OUT << Form("%16s || %5.2f ± %5.2f ± %5.2f || %5.2f ± %5.2f ± %5.2f || %5.2f ± %5.2f ± %5.2f ||", (it->first).c_str(), 
		       MM_yiel, sqrt(MM_stat*MM_stat), RareESyst*(MM_yiel),
		       EM_yiel, sqrt(EM_stat*MM_stat), RareESyst*(EM_yiel),
		       EE_yiel, sqrt(EE_stat*MM_stat), RareESyst*(EE_yiel)) << endl;
	}
	OUT << Form("%16s || %5.2f ? %5.2f ? %5.2f || %5.2f ? %5.2f ? %5.2f || %5.2f ? %5.2f ? %5.2f ||", "TTbarZ",
       nt2_ttz_mc_mm, sqrt(nt2_ttz_mc_mm_e2), RareESyst*(nt2_ttz_mc_mm),
       nt2_ttz_mc_em, sqrt(nt2_ttz_mc_em_e2), RareESyst*(nt2_ttz_mc_em),
       nt2_ttz_mc_ee, sqrt(nt2_ttz_mc_ee_e2), RareESyst*(nt2_ttz_mc_ee)) << endl;


	// do this later, we have the total yield anyway...
	// // adding ttz to the total rare MC background
	// nt2_rare_mc_mm += nt2_ttz_mc_mm; nt2_rare_mc_mm_e2 += nt2_ttz_mc_mm_e2;
	// nt2_rare_mc_em += nt2_ttz_mc_em; nt2_rare_mc_em_e2 += nt2_ttz_mc_em_e2;
	// nt2_rare_mc_ee += nt2_ttz_mc_ee; nt2_rare_mc_ee_e2 += nt2_ttz_mc_ee_e2;

	// adding only tt+W to the signal yield
	float nt2_sig_mc_mm = nt2_ttw_mc_mm; /*+ nt2_ttz_mc_mm;*/ float nt2_sig_mc_mm_e2 = nt2_ttw_mc_mm_e2; /*+ nt2_ttz_mc_mm_e2;*/
	float nt2_sig_mc_em = nt2_ttw_mc_em; /*+ nt2_ttz_mc_em;*/ float nt2_sig_mc_em_e2 = nt2_ttw_mc_em_e2; /*+ nt2_ttz_mc_em_e2;*/
	float nt2_sig_mc_ee = nt2_ttw_mc_ee; /*+ nt2_ttz_mc_ee;*/ float nt2_sig_mc_ee_e2 = nt2_ttw_mc_ee_e2; /*+ nt2_ttz_mc_ee_e2;*/

	OUT << "----------------------------------------------------------------------------------------------" << endl;
	// RARE SM BACKGROUND  in case of ttw == true, without TTW
	OUT << Form("%16s || %5.2f ± %5.2f ± %5.2f || %5.2f ± %5.2f ± %5.2f || %5.2f ± %5.2f ± %5.2f ||\n", "Rare SM (Sum)",
	nt2_rare_mc_mm + nt2_ttz_mc_mm, sqrt(nt2_rare_mc_mm_e2 + nt2_ttz_mc_mm_e2), RareESyst*(nt2_rare_mc_mm + nt2_ttz_mc_mm),
	nt2_rare_mc_em + nt2_ttz_mc_em, sqrt(nt2_rare_mc_em_e2 + nt2_ttz_mc_em_e2), RareESyst*(nt2_rare_mc_em + nt2_ttz_mc_em),
	nt2_rare_mc_ee + nt2_ttz_mc_ee, sqrt(nt2_rare_mc_ee_e2 + nt2_ttz_mc_ee_e2), RareESyst*(nt2_rare_mc_ee + nt2_ttz_mc_ee));
	OUT << "----------------------------------------------------------------------------------------------" << endl;

	// LINE FOR WZ BACKGROUND
	OUT << Form("%16s || %5.2f ± %5.2f ± %5.2f || %5.2f ± %5.2f ± %5.2f || %5.2f ± %5.2f ± %5.2f ||\n", "WZ prod.",
	nt2_wz_mc_mm, sqrt(nt2_wz_mc_mm_e2), WZESyst*nt2_wz_mc_mm,
	nt2_wz_mc_em, sqrt(nt2_wz_mc_em_e2), WZESyst*nt2_wz_mc_em,
	nt2_wz_mc_ee, sqrt(nt2_wz_mc_ee_e2), WZESyst*nt2_wz_mc_ee);
	OUT << "----------------------------------------------------------------------------------------------" << endl;
	// Just add different errors in quadrature (they are independent)
	float mm_tot_stat2 = FR->getMMTotEStat()*FR->getMMTotEStat() +                                   nt2_rare_mc_mm_e2 + nt2_ttz_mc_mm_e2 + nt2_wz_mc_mm_e2;
	float em_tot_stat2 = FR->getEMTotEStat()*FR->getEMTotEStat() + nt2_em_chmid_e1*nt2_em_chmid_e1 + nt2_rare_mc_em_e2 + nt2_ttz_mc_em_e2 + nt2_wz_mc_em_e2;
	float ee_tot_stat2 = FR->getEETotEStat()*FR->getEETotEStat() + nt2_ee_chmid_e1*nt2_ee_chmid_e1 + nt2_rare_mc_ee_e2 + nt2_ttz_mc_ee_e2 + nt2_wz_mc_ee_e2;

	//float mm_tot_syst2 = nF_mm*nF_mm*FakeESyst2 + RareESyst2*nt2_rare_mc_mm*nt2_rare_mc_mm + WZESyst2*nt2_wz_mc_mm_e2;
	//float em_tot_syst2 = nF_em*nF_em*FakeESyst2 + nt2_em_chmid_e2*nt2_em_chmid_e2 + RareESyst2*nt2_rare_mc_em*nt2_rare_mc_em + WZESyst2*nt2_wz_mc_em_e2;
	//float ee_tot_syst2 = nF_ee*nF_ee*FakeESyst2 + nt2_ee_chmid_e2*nt2_ee_chmid_e2 + RareESyst2*nt2_rare_mc_ee*nt2_rare_mc_ee + WZESyst2*nt2_wz_mc_ee_e2;

	// include ttZ in the total systematic error
	float mm_tot_syst2 = nF_mm*nF_mm*FakeESyst2 +                                   RareESyst2*((nt2_rare_mc_mm+nt2_ttz_mc_mm)*(nt2_rare_mc_mm+nt2_ttz_mc_mm)) + WZESyst2*nt2_wz_mc_mm_e2;
	float em_tot_syst2 = nF_em*nF_em*FakeESyst2 + nt2_em_chmid_e2*nt2_em_chmid_e2 + RareESyst2*((nt2_rare_mc_em+nt2_ttz_mc_em)*(nt2_rare_mc_em+nt2_ttz_mc_em)) + WZESyst2*nt2_wz_mc_em_e2;
	float ee_tot_syst2 = nF_ee*nF_ee*FakeESyst2 + nt2_ee_chmid_e2*nt2_ee_chmid_e2 + RareESyst2*((nt2_rare_mc_ee+nt2_ttz_mc_ee)*(nt2_rare_mc_ee+nt2_ttz_mc_ee)) + WZESyst2*nt2_wz_mc_ee_e2;

	// benjamins error calculation
	//float mm_tot_sqerr1 = FR->getMMTotEStat()*FR->getMMTotEStat() + nt2_rare_mc_mm_e1 + nt2_wz_mc_mm_e1;
	//float em_tot_sqerr1 = FR->getEMTotEStat()*FR->getEMTotEStat() + nt2_em_chmid_e1*nt2_em_chmid_e1 + nt2_rare_mc_em_e1 + nt2_wz_mc_em_e1;
	//float ee_tot_sqerr1 = FR->getEETotEStat()*FR->getEETotEStat() + nt2_ee_chmid_e1*nt2_ee_chmid_e1 + nt2_rare_mc_ee_e1 + nt2_wz_mc_ee_e1;
	//float mm_tot_sqerr2 = nF_mm*nF_mm*FakeESyst2 + RareESyst2*nt2_rare_mc_mm*nt2_rare_mc_mm + WZESyst*nt2_wz_mc_mm_e1*nt2_wz_mc_mm_e1;
	//float em_tot_sqerr2 = nF_em*nF_em*FakeESyst2 + nt2_em_chmid_e2*nt2_em_chmid_e2 + RareESyst2*nt2_rare_mc_em*nt2_rare_mc_em + WZESyst*nt2_wz_mc_em_e1*nt2_wz_mc_em_e1;
	//float ee_tot_sqerr2 = nF_ee*nF_ee*FakeESyst2 + nt2_ee_chmid_e2*nt2_ee_chmid_e2 + RareESyst2*nt2_rare_mc_ee*nt2_rare_mc_ee + WZESyst*nt2_wz_mc_ee_e1*nt2_wz_mc_ee_e1;

	OUT << Form("%16s || %5.2f ± %5.2f ± %5.2f || %5.2f ± %5.2f ± %5.2f || %5.2f ± %5.2f ± %5.2f ||\n", "tot. backgr.",
	nF_mm + nt2_wz_mc_mm +                nt2_rare_mc_mm + nt2_ttz_mc_mm, sqrt(mm_tot_stat2), sqrt(mm_tot_syst2),
	nF_em + nt2_wz_mc_em + nt2_em_chmid + nt2_rare_mc_em + nt2_ttz_mc_em, sqrt(em_tot_stat2), sqrt(em_tot_syst2),
	nF_ee + nt2_wz_mc_ee + nt2_ee_chmid + nt2_rare_mc_ee + nt2_ttz_mc_ee, sqrt(ee_tot_stat2), sqrt(ee_tot_syst2));
	// OUT << Form("%16s || %5.2f ± %5.2f         || %5.2f ± %5.2f         || %5.2f ± %5.2f         ||\n", "",
	// nF_mm + nt2_wz_mc_mm + nt2_rare_mc_mm,                sqrt(mm_tot_stat2 + mm_tot_syst2),
	// nF_em + nt2_wz_mc_em + nt2_em_chmid + nt2_rare_mc_em, sqrt(em_tot_stat2 + em_tot_syst2),
	// nF_ee + nt2_wz_mc_ee + nt2_ee_chmid + nt2_rare_mc_ee, sqrt(ee_tot_stat2 + ee_tot_syst2));
	OUT << "----------------------------------------------------------------------------------------------" << endl;
	// OUT << Form("%16s || %5.2f                 || %5.2f                 || %5.2f                 ||\n", "tot. MC", nt2sum_mm, nt2sum_em, nt2sum_ee);
	if (ttw) {
		OUT << "==============================================================================================" << endl;
		OUT << Form("%16s || %5.2f   %5.2f         || %5.2f   %5.2f         || %5.2f   %5.2f         ||\n", "ttW",
		nt2_ttw_mc_mm, sqrt(nt2_ttw_mc_mm_e2),
		nt2_ttw_mc_em, sqrt(nt2_ttw_mc_em_e2),
		nt2_ttw_mc_ee, sqrt(nt2_ttw_mc_ee_e2));
		OUT << Form("%16s || %5.2f   %5.2f         || %5.2f   %5.2f         || %5.2f   %5.2f         ||\n", "ttZ",
		nt2_ttz_mc_mm, sqrt(nt2_ttz_mc_mm_e2),
		nt2_ttz_mc_em, sqrt(nt2_ttz_mc_em_e2),
		nt2_ttz_mc_ee, sqrt(nt2_ttz_mc_ee_e2));
		OUT << "==============================================================================================" << endl;
	}
	OUT << "==============================================================================================" << endl;
	OUT << Form("%16s || %2.0f                    || %2.0f                    || %2.0f                    ||\n", "observed", nt2_mm, nt2_em, nt2_ee);
	OUT << "==============================================================================================" << endl;
	OUT << "        predicted: ";
	float tot_pred        = nF +           nt2_rare_mc_mm + nt2_ttz_mc_mm + nt2_wz_mc_mm +
	                        nt2_em_chmid + nt2_rare_mc_em + nt2_ttz_mc_em + nt2_wz_mc_em +
	                        nt2_ee_chmid + nt2_rare_mc_ee + nt2_ttz_mc_ee + nt2_wz_mc_ee ;
	float comb_tot_stat2  = FR->getTotEStat()*FR->getTotEStat() + nt2_rare_mc_mm_e2 + nt2_ttz_mc_mm_e2 + nt2_wz_mc_mm_e2 + 
	                        nt2_em_chmid_e1*nt2_em_chmid_e1     + nt2_rare_mc_em_e2 + nt2_ttz_mc_em_e2 + nt2_wz_mc_em_e2 + 
	                        nt2_ee_chmid_e1*nt2_ee_chmid_e1     + nt2_rare_mc_ee_e2 + nt2_ttz_mc_ee_e2 + nt2_wz_mc_ee_e2;
	float comb_tot_syst2  = nF*nF*FakeESyst2 + 
	                        RareESyst2*(nt2_rare_mc_mm + nt2_rare_mc_em + nt2_rare_mc_ee + nt2_ttz_mc_mm + nt2_ttz_mc_em + nt2_ttz_mc_ee)*
	                                   (nt2_rare_mc_mm + nt2_rare_mc_em + nt2_rare_mc_ee + nt2_ttz_mc_mm + nt2_ttz_mc_em + nt2_ttz_mc_ee)  +
	                        WZESyst2  *(nt2_wz_mc_mm   + nt2_wz_mc_em   + nt2_wz_mc_ee  )*(nt2_wz_mc_mm   + nt2_wz_mc_em   + nt2_wz_mc_ee) + 
	                        nt2_em_chmid_e2*nt2_em_chmid_e2 + nt2_ee_chmid_e2*nt2_ee_chmid_e2;
	OUT << setw(5) << left << Form("%5.2f", tot_pred ) << " ± ";
	OUT << setw(5) << Form("%5.2f", sqrt(comb_tot_stat2)) << " ± ";
	OUT << setw(5) << Form("%5.2f", sqrt(comb_tot_syst2)) << endl;
	// OUT << "      combined MC: ";
	// OUT << setw(5) << left << Form("%5.2f", nt2sum_mm+nt2sum_em+nt2sum_ee ) << endl;
	if (ttw) {
		// OUT << "ttW + ttZ signal : ";
		OUT << "ttW signal : ";
		OUT << setw(5) << left << Form("%5.2f", (nt2_sig_mc_mm + nt2_sig_mc_em + nt2_sig_mc_ee )) << " ±";
		OUT << setw(5) << left << Form("%5.2f", sqrt(nt2_sig_mc_mm_e2 + nt2_sig_mc_em_e2 + nt2_sig_mc_ee_e2)) << endl;
	}
	OUT << "combined observed: ";
	OUT << setw(5) << left << Form("%2.0f", nt2_mm+nt2_em+nt2_ee ) << endl;
	OUT << "==============================================================================================" << endl;
	OUT.close();
	
	///////////////////////////////////////////////////////////////////////////////////
	//  OUTPUT AS PLOT  ///////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////
	TH1D    *h_obs        = new TH1D("h_observed"   , "Observed number of events" , 4 , 0. , 4. );
	TH1D    *h_pred_fake  = new TH1D("h_pred_sfake" , "Predicted fakes"           , 4 , 0. , 4. );
	TH1D    *h_pred_chmid = new TH1D("h_pred_chmid" , "Predicted charge mis id"   , 4 , 0. , 4. );
	TH1D    *h_pred_mc    = new TH1D("h_pred_mc"    , "Predicted Rare SM"         , 4 , 0. , 4. );
	TH1D    *h_pred_wz    = new TH1D("h_pred_wz"    , "Predicted WZ"              , 4 , 0. , 4. );
	TH1D    *h_pred_ttw   = new TH1D("h_pred_ttw"   , "Predicted ttW"             , 4 , 0. , 4. );
	TH1D    *h_pred_ttz   = new TH1D("h_pred_ttz"   , "Predicted ttZ"             , 4 , 0. , 4. );
	TH1D    *h_pred_tot   = new TH1D("h_pred_tot"   , "Total Prediction"          , 4 , 0. , 4. );
	THStack *hs_pred      = new THStack("hs_predicted", "Predicted number of events");

	h_obs->SetMarkerColor(kBlack);
	h_obs->SetMarkerStyle(20);
	h_obs->SetMarkerSize(2.5);
	h_obs->SetLineWidth(2);
	h_obs->SetLineColor(kBlack);
	h_obs->SetFillColor(kBlack);

	h_pred_fake ->SetLineWidth(1);
	h_pred_chmid->SetLineWidth(1);
	h_pred_mc   ->SetLineWidth(1);
	h_pred_wz   ->SetLineWidth(1);
	h_pred_ttw  ->SetLineWidth(1);
	h_pred_ttz  ->SetLineWidth(1);
	h_pred_fake ->SetFillColor(46);
	h_pred_chmid->SetFillColor(49);
	h_pred_wz   ->SetFillColor(39);
	h_pred_mc   ->SetFillColor(38);
	h_pred_ttw  ->SetFillColor(44);
	h_pred_ttz  ->SetFillColor(42);
	// h_pred_ttz  ->SetLineColor(30);
	// h_pred_fake ->SetLineColor(46);
	// h_pred_chmid->SetLineColor(42);
	// h_pred_mc   ->SetLineColor(44);
	// h_pred_wz   ->SetLineColor(39);
	// h_pred_ttw  ->SetLineColor(29);

	h_pred_tot  ->SetLineWidth(3);
	// h_pred_tot  ->SetFillColor(kBlack);
	// h_pred_tot  ->SetFillStyle(3013);
	h_pred_tot  ->SetFillColor(12);
	h_pred_tot  ->SetFillStyle(3005);

	// Add numbers:
	h_obs->SetBinContent(1, nt2_ee);
	h_obs->SetBinContent(2, nt2_mm);
	h_obs->SetBinContent(3, nt2_em);
	h_obs->SetBinContent(4, nt2_ee+nt2_mm+nt2_em);

	TGraphAsymmErrors* gr_obs = FR->getGraphPoissonErrors( h_obs );
	gr_obs->SetMarkerColor(kBlack);
	gr_obs->SetMarkerStyle(20);
	gr_obs->SetMarkerSize(2.5);
	gr_obs->SetLineWidth(2);
	gr_obs->SetLineColor(kBlack);
	gr_obs->SetFillColor(kBlack);
	
	TString eeLabel  = "ee";
	TString emLabel  = "e#mu";
	TString mmLabel  = "#mu#mu";
	TString allLabel = "Total";

	if (chVeto == -1){
		eeLabel  = "e^{-}e^{-}";
		emLabel  = "e^{-}#mu^{-}";
		mmLabel  = "#mu^{-}#mu^{-}";
		// allLabel = "all (--)";
	}
	if (chVeto == +1){
		eeLabel  = "e^{+}e^{+}";
		emLabel  = "e^{+}#mu^{+}";
		mmLabel  = "#mu^{+}#mu^{+}";
		// allLabel = "all (++)";
	}
	
	h_pred_fake->SetBinContent(1, npf_ee+         nff_ee);
	h_pred_fake->SetBinContent(2, npf_mm+         nff_mm);
	h_pred_fake->SetBinContent(3, npf_em+nfp_em + nff_em);
	h_pred_fake->SetBinContent(4, npf_ee+npf_mm+npf_em+nfp_em + nff_ee+nff_mm+nff_em);
	h_pred_fake->GetXaxis()->SetBinLabel(1, eeLabel);
	h_pred_fake->GetXaxis()->SetBinLabel(2, mmLabel);
	h_pred_fake->GetXaxis()->SetBinLabel(3, emLabel);
	h_pred_fake->GetXaxis()->SetBinLabel(4, allLabel);

	h_pred_chmid->SetBinContent(1, nt2_ee_chmid);
	h_pred_chmid->SetBinContent(2, 0.);
	h_pred_chmid->SetBinContent(3, nt2_em_chmid);
	h_pred_chmid->SetBinContent(4, nt2_ee_chmid+nt2_em_chmid);

	h_pred_mc->SetBinContent(1, nt2_rare_mc_ee);
	h_pred_mc->SetBinContent(2, nt2_rare_mc_mm);
	h_pred_mc->SetBinContent(3, nt2_rare_mc_em);
	h_pred_mc->SetBinContent(4, nt2_rare_mc_ee+nt2_rare_mc_mm+nt2_rare_mc_em);

	h_pred_wz->SetBinContent(1, nt2_wz_mc_ee);
	h_pred_wz->SetBinContent(2, nt2_wz_mc_mm);
	h_pred_wz->SetBinContent(3, nt2_wz_mc_em);
	h_pred_wz->SetBinContent(4, nt2_wz_mc_ee+nt2_wz_mc_mm+nt2_wz_mc_em);

	h_pred_ttw->SetBinContent(1, nt2_ttw_mc_ee);
	h_pred_ttw->SetBinContent(2, nt2_ttw_mc_mm);
	h_pred_ttw->SetBinContent(3, nt2_ttw_mc_em);
	h_pred_ttw->SetBinContent(4, nt2_ttw_mc_ee+nt2_ttw_mc_mm+nt2_ttw_mc_em);

	h_pred_ttz->SetBinContent(1, nt2_ttz_mc_ee);
	h_pred_ttz->SetBinContent(1, nt2_ttz_mc_ee);
	h_pred_ttz->SetBinContent(2, nt2_ttz_mc_mm);
	h_pred_ttz->SetBinContent(3, nt2_ttz_mc_em);
	h_pred_ttz->SetBinContent(4, nt2_ttz_mc_ee+nt2_ttz_mc_mm+nt2_ttz_mc_em);

	h_pred_tot->Add(h_pred_fake);
	h_pred_tot->Add(h_pred_chmid);
	h_pred_tot->Add(h_pred_mc);
	h_pred_tot->Add(h_pred_wz);
	h_pred_tot->Add(h_pred_ttz);
	h_pred_tot->SetBinError(1, sqrt(ee_tot_stat2 + ee_tot_syst2));
	h_pred_tot->SetBinError(2, sqrt(mm_tot_stat2 + mm_tot_syst2));
	h_pred_tot->SetBinError(3, sqrt(em_tot_stat2 + em_tot_syst2));
	h_pred_tot->SetBinError(4, sqrt(comb_tot_stat2 + comb_tot_syst2));

	hs_pred->Add(h_pred_fake);
	hs_pred->Add(h_pred_chmid);
	hs_pred->Add(h_pred_mc);
	hs_pred->Add(h_pred_wz);
	hs_pred->Add(h_pred_ttz);
	hs_pred->Add(h_pred_ttw);

	
	// double max = h_obs->Integral();
	double max = 1.7*h_obs->GetBinContent(4);
	if (chVeto) max = 25.;

	// h_obs       ->SetMaximum(max>1?max+1:1.);
	// h_pred_sfake->SetMaximum(max>1?max+1:1.);
	// h_pred_dfake->SetMaximum(max>1?max+1:1.);
	// h_pred_chmid->SetMaximum(max>1?max+1:1.);
	// h_pred_mc   ->SetMaximum(max>1?max+1:1.);
	// h_pred_wz   ->SetMaximum(max>1?max+1:1.);
	// h_pred_ttw  ->SetMaximum(max>1?max+1:1.);
	// h_pred_ttz  ->SetMaximum(max>1?max+1:1.);
	// h_pred_tot  ->SetMaximum(max>1?max+1:1.);
	// hs_pred     ->SetMaximum(max>1?max+1:1.);
	h_obs       ->SetMaximum(max);
	h_pred_fake ->SetMaximum(max);
	h_pred_chmid->SetMaximum(max);
	h_pred_mc   ->SetMaximum(max);
	h_pred_wz   ->SetMaximum(max);
	h_pred_ttw  ->SetMaximum(max);
	h_pred_ttz  ->SetMaximum(max);
	h_pred_tot  ->SetMaximum(max);
	hs_pred     ->SetMaximum(max);
	hs_pred->Draw("goff");
	hs_pred->GetXaxis()->SetBinLabel(1, "ee");
	hs_pred->GetXaxis()->SetBinLabel(2, "#mu#mu");
	hs_pred->GetXaxis()->SetBinLabel(3, "e#mu");
	hs_pred->GetXaxis()->SetBinLabel(4, "all");
	hs_pred->GetXaxis()->SetLabelOffset(0.01);
	hs_pred->GetXaxis()->SetLabelFont(42);
	hs_pred->GetXaxis()->SetLabelSize(0.1);
	hs_pred->GetYaxis()->SetTitle("Events");

	TLegend *leg = new TLegend(0.15,0.62,0.50,0.88);
	leg->AddEntry(h_obs,        "Observed","p");
	leg->AddEntry(h_pred_fake , "Fakes","f");
	leg->AddEntry(h_pred_chmid, "Charge MisID","f");
	leg->AddEntry(h_pred_mc,    "Irreducible (MC)","f");
	leg->AddEntry(h_pred_wz,    "WZ Production","f");
	leg->AddEntry(h_pred_ttz,   "ttZ Production","f");
	leg->AddEntry(h_pred_ttw,   "ttW Production","f");
	//leg->AddEntry(h_pred_tot,   "Total Uncertainty","f");
	leg->SetFillStyle(0);
	leg->SetTextFont(42);
	// leg->SetTextSize(0.05);
	leg->SetBorderSize(0);

	TCanvas *c_temp = new TCanvas("C_ObsPred", "Observed vs Predicted", 0, 0, 600, 600);
	c_temp->cd();
	
	hs_pred->Draw("hist");
	leg->Draw();
	h_pred_tot->DrawCopy("0 E2 same");
	h_pred_tot->SetFillStyle(0);
	h_pred_tot->DrawCopy("hist same");
	gr_obs->Draw("P same");
	
	lat->SetTextSize(0.03);

	if (maxHT < 19.) lat->DrawLatex(0.45,0.85, "N_{Jets} = 0");
	else if (maxHT > 1000.) lat->DrawLatex(0.45,0.85, Form("H_{T} > %.0f GeV, N_{Jets} #geq %1d", minHT, minNjets));
	else lat->DrawLatex(0.45,0.85, Form("%.0f GeV < H_{T} < %.0f GeV, N_{Jets} #geq %1d", minHT, maxHT, minNjets));

	if (minMET > 0. ) {
		if (maxMET > 1000.) lat->DrawLatex(0.45,0.80, Form("E_{T}^{miss} > %.0f GeV", minMET));
		else lat->DrawLatex(0.45,0.80, Form("%.0f GeV < E_{T}^{miss} < %.0f GeV", minMET, maxMET));
	}

	float yPos = (minMET > 0. ? 0.75:0.8);
	if (minNbjetsL != 0) lat->DrawLatex(0.45, yPos, Form("N_{b-Jets} #geq %.0d", minNbjetsL));

	if (chVeto != 0) {
		lat->SetTextFont(62);
		if (chVeto == -1) lat->DrawLatex(0.45, yPos-0.05, "-- only");
		if (chVeto == +1) lat->DrawLatex(0.45, yPos-0.05, "++ only");
	}

	drawTopLine(0.56, 0.8);
	
	gPad->RedrawAxis();
	// Util::PrintNoEPS(c_temp, "ObsPred_" + Region::sname[reg], fOutputDir + fOutputSubDir, NULL);
	Util::PrintPDF (c_temp,   Form("ObsPred_customRegion_HT%.0f"+jvString+"MET%.0fNJ%.0iNbjL%.0iNbjM%.0iPT1%.0fPT2%.0f_"+sysString+chargeString, minHT, minMET, minNjets, minNbjetsL, minNbjetsM, minPt1, minPt2) , fOutputDir + fOutputSubDir);
	Util::PrintROOT(c_temp,   Form("ObsPred_customRegion_HT%.0f"+jvString+"MET%.0fNJ%.0iNbjL%.0iNbjM%.0iPT1%.0fPT2%.0f_"+sysString+chargeString, minHT, minMET, minNjets, minNbjetsL, minNbjetsM, minPt1, minPt2) , fOutputDir + fOutputSubDir);

	delete c_temp;	
	delete h_obs, h_pred_fake, h_pred_chmid, h_pred_mc, h_pred_tot, hs_pred;
	delete gr_obs;
	delete FR;


	TTWZPrediction pred;
	pred.obs      = nt2_mm+nt2_em+nt2_ee;
	pred.obs_mm   = nt2_mm;
	pred.obs_em   = nt2_em;
	pred.obs_ee   = nt2_ee;
	
	pred.ttw      = nt2_ttw_mc_mm + nt2_ttw_mc_em + nt2_ttw_mc_ee;
	pred.ttw_mm   = nt2_ttw_mc_mm;
	pred.ttw_em   = nt2_ttw_mc_em;
	pred.ttw_ee   = nt2_ttw_mc_ee;
	
	pred.ttz      = nt2_ttz_mc_mm + nt2_ttz_mc_em + nt2_ttz_mc_ee;
	pred.ttz_mm   = nt2_ttz_mc_mm;
	pred.ttz_em   = nt2_ttz_mc_em;
	pred.ttz_ee   = nt2_ttz_mc_ee;
	
	pred.ttz_err    = sqrt(nt2_ttz_mc_ee_e2 + nt2_ttz_mc_mm_e2 + nt2_ttz_mc_em_e2 + TTZESyst2*(nt2_ttz_mc_ee + nt2_ttz_mc_mm + nt2_ttz_mc_em)*(nt2_ttz_mc_ee + nt2_ttz_mc_mm + nt2_ttz_mc_em));
	pred.ttz_err_mm = sqrt(nt2_ttz_mc_mm_e2 + TTZESyst2*nt2_ttz_mc_mm*nt2_ttz_mc_mm);
	pred.ttz_err_ee = sqrt(nt2_ttz_mc_em_e2 + TTZESyst2*nt2_ttz_mc_em*nt2_ttz_mc_em);
	pred.ttz_err_em = sqrt(nt2_ttz_mc_ee_e2 + TTZESyst2*nt2_ttz_mc_ee*nt2_ttz_mc_ee);
	
	pred.ttwz     = pred.ttw + pred.ttz;
	pred.ttwz_mm  = pred.ttw_mm + pred.ttz_mm;
	pred.ttwz_em  = pred.ttw_em + pred.ttz_em;
	pred.ttwz_ee  = pred.ttw_ee + pred.ttz_ee;
	
	pred.fake     = nF;
	pred.fake_mm  = nF_mm;
	pred.fake_em  = nF_em;
	pred.fake_ee  = nF_ee;
	
	pred.fake_err    = sqrt(FR->getTotEStat()*FR->getTotEStat() + FakeESyst2*nF*nF);
	pred.fake_err_mm = sqrt(FR->getMMTotEStat()*FR->getMMTotEStat() + FakeESyst2*nF_mm*nF_mm);
	pred.fake_err_em = sqrt(FR->getEMTotEStat()*FR->getEMTotEStat() + FakeESyst2*nF_em*nF_em);
	pred.fake_err_ee = sqrt(FR->getEETotEStat()*FR->getEETotEStat() + FakeESyst2*nF_ee*nF_ee);
	
	pred.cmid        = nt2_ee_chmid + nt2_em_chmid;
	pred.cmid_em     = nt2_em_chmid;
	pred.cmid_ee     = nt2_ee_chmid;
	pred.cmid_err    = sqrt(nt2_ee_chmid_e1*nt2_ee_chmid_e1 + nt2_ee_chmid_e2*nt2_ee_chmid_e2 + nt2_em_chmid_e1*nt2_em_chmid_e1 + nt2_em_chmid_e2*nt2_em_chmid_e2);
	pred.cmid_err_ee = sqrt(nt2_ee_chmid_e1*nt2_ee_chmid_e1 + nt2_ee_chmid_e2*nt2_ee_chmid_e2);
	pred.cmid_err_em = sqrt(nt2_em_chmid_e1*nt2_em_chmid_e1 + nt2_em_chmid_e2*nt2_em_chmid_e2);
	
	pred.wz        = nt2_wz_mc_mm + nt2_wz_mc_em + nt2_wz_mc_ee;
	pred.wz_mm     = nt2_wz_mc_mm;
	pred.wz_em     = nt2_wz_mc_em;
	pred.wz_ee     = nt2_wz_mc_ee;
	pred.wz_err    = sqrt(nt2_wz_mc_mm_e2 + nt2_wz_mc_em_e2 + nt2_wz_mc_ee_e2 + WZESyst2*(nt2_wz_mc_mm + nt2_wz_mc_em + nt2_wz_mc_ee)*(nt2_wz_mc_mm + nt2_wz_mc_em + nt2_wz_mc_ee));
	pred.wz_err_mm = sqrt(nt2_wz_mc_mm_e2  + WZESyst2*nt2_wz_mc_mm*nt2_wz_mc_mm);
	pred.wz_err_em = sqrt(nt2_wz_mc_em_e2  + WZESyst2*nt2_wz_mc_em*nt2_wz_mc_em);
	pred.wz_err_ee = sqrt(nt2_wz_mc_ee_e2  + WZESyst2*nt2_wz_mc_ee*nt2_wz_mc_ee);
	
	pred.rare        = nt2_rare_mc_ee + nt2_rare_mc_mm + nt2_rare_mc_em;
	pred.rare_mm     = nt2_rare_mc_mm;
	pred.rare_em     = nt2_rare_mc_em;
	pred.rare_ee     = nt2_rare_mc_ee;
	pred.rare_err    = sqrt(nt2_rare_mc_ee_e2 + nt2_rare_mc_mm_e2 + nt2_rare_mc_em_e2 + RareESyst2*(nt2_rare_mc_ee + nt2_rare_mc_mm + nt2_rare_mc_em)*(nt2_rare_mc_ee + nt2_rare_mc_mm + nt2_rare_mc_em));
	pred.rare_err_mm = sqrt(nt2_rare_mc_mm_e2 + RareESyst2*nt2_rare_mc_mm*nt2_rare_mc_mm);
	pred.rare_err_em = sqrt(nt2_rare_mc_em_e2 + RareESyst2*nt2_rare_mc_em*nt2_rare_mc_em);
	pred.rare_err_ee = sqrt(nt2_rare_mc_ee_e2 + RareESyst2*nt2_rare_mc_ee*nt2_rare_mc_ee);

	return pred;

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
	
	bool plotRatios = true;
	
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

	// {  0 ,   1  ,    2   ,   3  ,   4  ,   5  ,    6    ,   7   ,      8     ,      9      }
	// {"HT", "MET", "NJets", "MT2", "PT1", "PT2", "NBJets", "MET3", "NBJetsMed", "NBJetsMed2"}
	float binwidthscale[gNDiffVars] = {100., 20., 1., 25., 20., 10., 1., 10., 1., 1.};

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
		// 	float scale = fLumiNorm / S->getLumi();
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
			Sample *S = fSamples[fMCRareSM[i]];
			if (fMCRareSM[i] == WZ) continue;
			float scale = fLumiNorm / S->getLumi();
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

				float estat2_mm = scale*scale*S->getError2(ss_mm);
				float estat2_ee = scale*scale*S->getError2(ss_ee);
				float estat2_em = scale*scale*S->getError2(ss_em);
				// float estat2_mm = scale*scale*S->numbers[reg][Muon].tt_avweight*S->numbers[reg][Muon].tt_avweight*S->getError2(ss_mm); // which region to take for weighing of mc errors?
				// float estat2_ee = scale*scale*S->numbers[reg][Elec].tt_avweight*S->numbers[reg][Elec].tt_avweight*S->getError2(ss_ee);
				// float estat2_em = scale*scale*S->numbers[reg][ElMu].tt_avweight*S->numbers[reg][ElMu].tt_avweight*S->getError2(ss_em);

				float prev    = totbg   ->GetBinError(b+1);
				float prev_mm = totbg_mm->GetBinError(b+1);
				float prev_em = totbg_em->GetBinError(b+1);
				float prev_ee = totbg_ee->GetBinError(b+1);

				totbg   ->SetBinError(b+1, prev    + esyst2_mm + esyst2_ee + esyst2_em + estat2_mm + estat2_ee + estat2_em);
				totbg_mm->SetBinError(b+1, prev_mm + esyst2_mm + estat2_mm);
				totbg_em->SetBinError(b+1, prev_em + esyst2_em + estat2_em);
				totbg_ee->SetBinError(b+1, prev_ee + esyst2_ee + estat2_ee);
			}
		}

		nt11_ss->Add(nt11_mm_ss);
		nt11_ss->Add(nt11_ee_ss);
		nt11_ss->Add(nt11_em_ss);

		totbg   ->Add(nt11_ss);
		totbg_mm->Add(nt11_mm_ss);
		totbg_em->Add(nt11_em_ss);
		totbg_ee->Add(nt11_ee_ss);

		///////////////////////////////////////////////////////////////////////////////////
		// WZ PRODUCTION //////////////////////////////////////////////////////////////////
		///////////////////////////////////////////////////////////////////////////////////
		TH1D *nt11_wz = new TH1D(Form("NT11_WZ_%s", varname.Data()), varname, nbins, bins); nt11_wz->Sumw2();

		TH1D *nt11_mm_wz = new TH1D(Form("NT11_MM_WZ_%s", varname.Data()), varname, nbins, bins); nt11_mm_wz->Sumw2();
		TH1D *nt11_ee_wz = new TH1D(Form("NT11_EE_WZ_%s", varname.Data()), varname, nbins, bins); nt11_ee_wz->Sumw2();
		TH1D *nt11_em_wz = new TH1D(Form("NT11_EM_WZ_%s", varname.Data()), varname, nbins, bins); nt11_em_wz->Sumw2();

		float wzscale = fLumiNorm / fSamples[WZ]->getLumi();
		nt11_mm_wz->Add(fSamples[WZ]->diffyields[Muon].hnt11[j], wzscale);
		nt11_ee_wz->Add(fSamples[WZ]->diffyields[Elec].hnt11[j], wzscale);
		nt11_em_wz->Add(fSamples[WZ]->diffyields[ElMu].hnt11[j], wzscale);

		// Errors
		for(size_t b = 0; b < nbins; ++b){
			float ss_mm = fSamples[WZ]->diffyields[Muon].hnt11[j]->GetBinContent(b+1);
			float ss_ee = fSamples[WZ]->diffyields[Elec].hnt11[j]->GetBinContent(b+1);
			float ss_em = fSamples[WZ]->diffyields[ElMu].hnt11[j]->GetBinContent(b+1);

			float esyst2_mm = 0.25 * ss_mm*ss_mm*wzscale*wzscale;
			float esyst2_ee = 0.25 * ss_ee*ss_ee*wzscale*wzscale;
			float esyst2_em = 0.25 * ss_em*ss_em*wzscale*wzscale;

			float estat2_mm = wzscale*wzscale*fSamples[WZ]->getError2(ss_mm);
			float estat2_ee = wzscale*wzscale*fSamples[WZ]->getError2(ss_ee);
			float estat2_em = wzscale*wzscale*fSamples[WZ]->getError2(ss_em);
			// float estat2_mm = wzscale*wzscale*fSamples[WZ]->numbers[reg][Muon].tt_avweight*fSamples[WZ]->numbers[reg][Muon].tt_avweight*fSamples[WZ]->getError2(ss_mm);
			// float estat2_ee = wzscale*wzscale*fSamples[WZ]->numbers[reg][Elec].tt_avweight*fSamples[WZ]->numbers[reg][Elec].tt_avweight*fSamples[WZ]->getError2(ss_ee);
			// float estat2_em = wzscale*wzscale*fSamples[WZ]->numbers[reg][ElMu].tt_avweight*fSamples[WZ]->numbers[reg][ElMu].tt_avweight*fSamples[WZ]->getError2(ss_em);

			float prev    = totbg   ->GetBinError(b+1);
			float prev_mm = totbg_mm->GetBinError(b+1);
			float prev_em = totbg_em->GetBinError(b+1);
			float prev_ee = totbg_ee->GetBinError(b+1);

			totbg   ->SetBinError(b+1, prev    + esyst2_mm + esyst2_ee + esyst2_em + estat2_mm + estat2_ee + estat2_em);
			totbg_mm->SetBinError(b+1, prev_mm + esyst2_mm + estat2_mm);
			totbg_em->SetBinError(b+1, prev_em + esyst2_em + estat2_em);
			totbg_ee->SetBinError(b+1, prev_ee + esyst2_ee + estat2_ee);
		}

		nt11_wz->Add(nt11_mm_wz);
		nt11_wz->Add(nt11_ee_wz);
		nt11_wz->Add(nt11_em_wz);

		totbg   ->Add(nt11_wz);
		totbg_mm->Add(nt11_mm_wz);
		totbg_em->Add(nt11_em_wz);
		totbg_ee->Add(nt11_ee_wz);

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

		/////////////////////////////////////////////////////
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
		/////////////////////////////////////////////////////

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
			
			///////////////////////////////////////////////////////
			// Flat ratios:
			// nt11_mm_sf->SetBinContent(i+1, FR->getMMNpf());
			// nt11_ee_sf->SetBinContent(i+1, FR->getEENpf());
			// nt11_em_sf->SetBinContent(i+1, FR->getEMNpf() + FR->getEMNfp());
			// nt11_mm_df->SetBinContent(i+1, FR->getMMNff());
			// nt11_ee_df->SetBinContent(i+1, FR->getEENff());
			// nt11_em_df->SetBinContent(i+1, FR->getEMNff());
			///////////////////////////////////////////////////////
			
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
		
		if (gBaseRegion != "WZEnriched") {
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
		}
		///////////////////////////////////////////////////////////////////////////////////
		// E-CHARGE MISID /////////////////////////////////////////////////////////////////
		///////////////////////////////////////////////////////////////////////////////////
		TH1D *nt11_cm = new TH1D(Form("NT11_CM_%s", varname.Data()), varname, nbins, bins); nt11_cm->Sumw2();

		TH1D *nt11_ee_cm = new TH1D(Form("NT11_EE_CM_%s", varname.Data()), varname, nbins, bins); nt11_ee_cm->Sumw2();
		TH1D *nt11_em_cm = new TH1D(Form("NT11_EM_CM_%s", varname.Data()), varname, nbins, bins); nt11_em_cm->Sumw2();

		// Abbreviations
		float fbb(0.),fee(0.),feb(0.);
		float fbbE(0.),feeE(0.),febE(0.);
		
		calculateChMisIdProb(fEGData, BB, fbb, fbbE);
		calculateChMisIdProb(fEGData, EB, feb, febE);
		calculateChMisIdProb(fEGData, EE, fee, feeE);

		for(size_t i = 0; i < nbins; ++i){
			float nt2_ee_BB_os = nt2_os_ee_bb->GetBinContent(i+1);
			float nt2_ee_EB_os = nt2_os_ee_eb->GetBinContent(i+1);
			float nt2_ee_EE_os = nt2_os_ee_ee->GetBinContent(i+1);
			float nt2_em_BB_os = nt2_os_em_bb->GetBinContent(i+1);
			float nt2_em_EE_os = nt2_os_em_ee->GetBinContent(i+1);
			
			// Errors
			FakeRatios *FR = new FakeRatios();

			// Simple error propagation assuming error on number of events is FR->getEStat2()
			nt11_ee_cm->SetBinContent(i+1, 2*fbb*nt2_ee_BB_os + 2*fee*nt2_ee_EE_os + 2*feb*nt2_ee_EB_os);
			float nt11_ee_cm_e1 = sqrt( (4*fbb*fbb*FR->getEStat2(nt2_ee_BB_os)) + (4*fee*fee*FR->getEStat2(nt2_ee_EE_os)) + 4*feb*feb*FR->getEStat2(nt2_ee_EB_os) ); // stat only
			float nt11_ee_cm_e2 = sqrt( (4*nt2_ee_BB_os*nt2_ee_BB_os*fbbE*fbbE) + (4*nt2_ee_EE_os*nt2_ee_EE_os*feeE*feeE) + 4*febE*febE*nt2_ee_EB_os*nt2_ee_EB_os ); // syst only

			nt11_em_cm->SetBinContent(i+i, fbb*nt2_em_BB_os + fee*nt2_em_EE_os);
			float nt11_em_cm_e1 = sqrt( fbb*fbb*FR->getEStat2(nt2_em_BB_os) + fee*fee*FR->getEStat2(nt2_em_EE_os) );
			float nt11_em_cm_e2 = sqrt( nt2_em_BB_os*nt2_em_BB_os * fbbE*fbbE + nt2_em_EE_os*nt2_em_EE_os * feeE*feeE );
			
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
// 		gSample sigsam = LM11;
// 		TH1D *nt11_sig    = new TH1D(Form("NT11_Sig%s",    varname.Data()), varname, nbins, bins); nt11_sig   ->Sumw2();
// 		TH1D *nt11_mm_sig = new TH1D(Form("NT11_mm_Sig%s", varname.Data()), varname, nbins, bins); nt11_mm_sig->Sumw2();
// 		TH1D *nt11_em_sig = new TH1D(Form("NT11_em_Sig%s", varname.Data()), varname, nbins, bins); nt11_em_sig->Sumw2();
// 		TH1D *nt11_ee_sig = new TH1D(Form("NT11_ee_Sig%s", varname.Data()), varname, nbins, bins); nt11_ee_sig->Sumw2();
// 		nt11_mm_sig->Add(fSamples[sigsam]->diffyields[Muon].hnt11[j], fLumiNorm / fSamples[LM4]->getLumi());
// 		nt11_ee_sig->Add(fSamples[sigsam]->diffyields[Elec].hnt11[j], fLumiNorm / fSamples[LM4]->getLumi());
// 		nt11_em_sig->Add(fSamples[sigsam]->diffyields[ElMu].hnt11[j], fLumiNorm / fSamples[LM4]->getLumi());
// 		nt11_sig   ->Add(fSamples[sigsam]->diffyields[Muon].hnt11[j], fLumiNorm / fSamples[LM4]->getLumi());
// 		nt11_sig   ->Add(fSamples[sigsam]->diffyields[Elec].hnt11[j], fLumiNorm / fSamples[LM4]->getLumi());
// 		nt11_sig   ->Add(fSamples[sigsam]->diffyields[ElMu].hnt11[j], fLumiNorm / fSamples[LM4]->getLumi());

		/////////////////
		// RATIO ////////
		/////////////////
		TH1D *rat    = new TH1D(Form("RAT_%s",varname.Data()), varname, nbins, bins);    rat->Sumw2();
		TH1D *rat_mm = new TH1D(Form("RAT_MM_%s",varname.Data()), varname, nbins, bins); rat_mm->Sumw2();
		TH1D *rat_ee = new TH1D(Form("RAT_EE_%s",varname.Data()), varname, nbins, bins); rat_ee->Sumw2();
		TH1D *rat_em = new TH1D(Form("RAT_EM_%s",varname.Data()), varname, nbins, bins); rat_em->Sumw2();
		
		// Ratio plot:
		float border = 0.3;
		float scale = (1-border)/border;
		
		rat->SetFillStyle(1001);
		rat->SetLineWidth(1);
		rat->SetLineWidth(1);
		rat->SetFillColor(  kGray+1);
		rat->SetLineColor(  kGray+1);
		rat->SetMarkerColor(kGray+1);
		rat_mm->SetFillStyle(1001);
		rat_mm->SetLineWidth(1);
		rat_mm->SetLineWidth(1);
		rat_mm->SetFillColor(  kGray+1);
		rat_mm->SetLineColor(  kGray+1);
		rat_mm->SetMarkerColor(kGray+1);
		rat_em->SetFillStyle(1001);
		rat_em->SetLineWidth(1);
		rat_em->SetLineWidth(1);
		rat_em->SetFillColor(  kGray+1);
		rat_em->SetLineColor(  kGray+1);
		rat_em->SetMarkerColor(kGray+1);
		rat_ee->SetFillStyle(1001);
		rat_ee->SetLineWidth(1);
		rat_ee->SetLineWidth(1);
		rat_ee->SetFillColor(  kGray+1);
		rat_ee->SetLineColor(  kGray+1);
		rat_ee->SetMarkerColor(kGray+1);
		
		
		rat->SetXTitle(DiffPredYields::axis_label[j]);
		rat->SetYTitle("");
		rat->GetXaxis()->SetTitleSize(scale * 0.04);
		rat->GetXaxis()->SetLabelSize(scale * nt11->GetXaxis()->GetLabelSize());
		rat->GetYaxis()->SetLabelSize(scale * nt11->GetYaxis()->GetLabelSize());
		rat->GetXaxis()->SetTickLength(scale * nt11->GetXaxis()->GetTickLength());
		rat->GetYaxis()->SetTickLength(nt11->GetYaxis()->GetTickLength());

		rat_mm->SetXTitle(DiffPredYields::axis_label[j]);
		rat_mm->SetYTitle("");
		rat_mm->GetXaxis()->SetTitleSize(scale * 0.04);
		rat_mm->GetXaxis()->SetLabelSize(scale * nt11_mm->GetXaxis()->GetLabelSize());
		rat_mm->GetYaxis()->SetLabelSize(scale * nt11_mm->GetYaxis()->GetLabelSize());
		rat_mm->GetXaxis()->SetTickLength(scale * nt11_mm->GetXaxis()->GetTickLength());
		rat_mm->GetYaxis()->SetTickLength(nt11_mm->GetYaxis()->GetTickLength());

		rat_em->SetXTitle(DiffPredYields::axis_label[j]);
		rat_em->SetYTitle("");
		rat_em->GetXaxis()->SetTitleSize(scale * 0.04);
		rat_em->GetXaxis()->SetLabelSize(scale * nt11_em->GetXaxis()->GetLabelSize());
		rat_em->GetYaxis()->SetLabelSize(scale * nt11_em->GetYaxis()->GetLabelSize());
		rat_em->GetXaxis()->SetTickLength(scale * nt11_em->GetXaxis()->GetTickLength());
		rat_em->GetYaxis()->SetTickLength(nt11_em->GetYaxis()->GetTickLength());

		rat_ee->SetXTitle(DiffPredYields::axis_label[j]);
		rat_ee->SetYTitle("");
		rat_ee->GetXaxis()->SetTitleSize(scale * 0.04);
		rat_ee->GetXaxis()->SetLabelSize(scale * nt11_ee->GetXaxis()->GetLabelSize());
		rat_ee->GetYaxis()->SetLabelSize(scale * nt11_ee->GetYaxis()->GetLabelSize());
		rat_ee->GetXaxis()->SetTickLength(scale * nt11_ee->GetXaxis()->GetTickLength());
		rat_ee->GetYaxis()->SetTickLength(nt11_ee->GetYaxis()->GetTickLength());

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
		nt11_wz->SetLineWidth(1);
		nt11_sf->SetLineColor(50);
		nt11_sf->SetFillColor(50);
		nt11_df->SetLineColor(38);
		nt11_df->SetFillColor(38);
		nt11_cm->SetLineColor(42);
		nt11_cm->SetFillColor(42);
		nt11_ss->SetLineColor(31);
		nt11_ss->SetFillColor(31);
		nt11_wz->SetLineColor(29);
		nt11_wz->SetFillColor(29);

		nt11_mm_sf->SetLineWidth(1);
		nt11_mm_df->SetLineWidth(1);
		nt11_mm_ss->SetLineWidth(1);
		nt11_mm_wz->SetLineWidth(1);
		nt11_mm_sf->SetLineColor(50);
		nt11_mm_sf->SetFillColor(50);
		nt11_mm_df->SetLineColor(38);
		nt11_mm_df->SetFillColor(38);
		nt11_mm_ss->SetLineColor(31);
		nt11_mm_ss->SetFillColor(31);
		nt11_mm_wz->SetLineColor(29);
		nt11_mm_wz->SetFillColor(29);

		nt11_ee_sf->SetLineWidth(1);
		nt11_ee_df->SetLineWidth(1);
		nt11_ee_ss->SetLineWidth(1);
		nt11_ee_wz->SetLineWidth(1);
		nt11_ee_sf->SetLineColor(50);
		nt11_ee_sf->SetFillColor(50);
		nt11_ee_df->SetLineColor(38);
		nt11_ee_df->SetFillColor(38);
		nt11_ee_cm->SetLineColor(42);
		nt11_ee_cm->SetFillColor(42);
		nt11_ee_ss->SetLineColor(31);
		nt11_ee_ss->SetFillColor(31);
		nt11_ee_wz->SetLineColor(29);
		nt11_ee_wz->SetFillColor(29);

		nt11_em_sf->SetLineWidth(1);
		nt11_em_df->SetLineWidth(1);
		nt11_em_ss->SetLineWidth(1);
		nt11_em_wz->SetLineWidth(1);
		nt11_em_sf->SetLineColor(50);
		nt11_em_sf->SetFillColor(50);
		nt11_em_df->SetLineColor(38);
		nt11_em_df->SetFillColor(38);
		nt11_em_cm->SetLineColor(42);
		nt11_em_cm->SetFillColor(42);
		nt11_em_ss->SetLineColor(31);
		nt11_em_ss->SetFillColor(31);
		nt11_em_wz->SetLineColor(29);
		nt11_em_wz->SetFillColor(29);
		
// 		nt11_sig   ->SetLineWidth(2);
// 		nt11_mm_sig->SetLineWidth(2);
// 		nt11_em_sig->SetLineWidth(2);
// 		nt11_ee_sig->SetLineWidth(2);

// 		nt11_sig   ->SetLineColor(kBlue);
// 		nt11_mm_sig->SetLineColor(kBlue);
// 		nt11_em_sig->SetLineColor(kBlue);
// 		nt11_ee_sig->SetLineColor(kBlue);
		
// 		nt11_sig   ->SetFillStyle(0);
// 		nt11_mm_sig->SetFillStyle(0);
// 		nt11_em_sig->SetFillStyle(0);
// 		nt11_ee_sig->SetFillStyle(0);
		
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
		nt11_wz    = normHistBW(nt11_wz, binwidthscale[j]);
		nt11_cm    = normHistBW(nt11_cm, binwidthscale[j]);

		nt11_mm_sf = normHistBW(nt11_mm_sf, binwidthscale[j]);
		nt11_mm_df = normHistBW(nt11_mm_df, binwidthscale[j]);
		nt11_mm_ss = normHistBW(nt11_mm_ss, binwidthscale[j]);
		nt11_mm_wz = normHistBW(nt11_mm_wz, binwidthscale[j]);

		nt11_ee_sf = normHistBW(nt11_ee_sf, binwidthscale[j]);
		nt11_ee_df = normHistBW(nt11_ee_df, binwidthscale[j]);
		nt11_ee_ss = normHistBW(nt11_ee_ss, binwidthscale[j]);
		nt11_ee_wz = normHistBW(nt11_ee_wz, binwidthscale[j]);
		nt11_ee_cm = normHistBW(nt11_ee_cm, binwidthscale[j]);
		
		nt11_em_sf = normHistBW(nt11_em_sf, binwidthscale[j]);
		nt11_em_df = normHistBW(nt11_em_df, binwidthscale[j]);
		nt11_em_ss = normHistBW(nt11_em_ss, binwidthscale[j]);
		nt11_em_wz = normHistBW(nt11_em_wz, binwidthscale[j]);
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
		nt11_tot->Add(nt11_wz);
		nt11_tot->Add(nt11_cm);

		nt11_mm_tot->Add(nt11_mm_sf);
		nt11_mm_tot->Add(nt11_mm_df);
		nt11_mm_tot->Add(nt11_mm_ss);
		nt11_mm_tot->Add(nt11_mm_wz);

		nt11_ee_tot->Add(nt11_ee_sf);
		nt11_ee_tot->Add(nt11_ee_df);
		nt11_ee_tot->Add(nt11_ee_ss);
		nt11_ee_tot->Add(nt11_ee_wz);
		nt11_ee_tot->Add(nt11_ee_cm);

		nt11_em_tot->Add(nt11_em_sf);
		nt11_em_tot->Add(nt11_em_df);
		nt11_em_tot->Add(nt11_em_ss);
		nt11_em_tot->Add(nt11_em_wz);
		nt11_em_tot->Add(nt11_em_cm);
		
		TH1D *h_nt11    = (TH1D*) nt11_sf->Clone("h_nt11");       h_nt11   ->Sumw2();
		TH1D *h_nt11_mm = (TH1D*) nt11_mm_sf->Clone("h_nt11_mm"); h_nt11_mm->Sumw2();
		TH1D *h_nt11_em = (TH1D*) nt11_em_sf->Clone("h_nt11_em"); h_nt11_em->Sumw2();
		TH1D *h_nt11_ee = (TH1D*) nt11_ee_sf->Clone("h_nt11_ee"); h_nt11_ee->Sumw2();
		
		h_nt11->Add(nt11_df);
		h_nt11->Add(nt11_ss);
		h_nt11->Add(nt11_wz);
		h_nt11->Add(nt11_cm);

		h_nt11_mm->Add(nt11_mm_df);
		h_nt11_mm->Add(nt11_mm_ss);
		h_nt11_mm->Add(nt11_mm_wz);
		
		h_nt11_ee->Add(nt11_ee_df);
		h_nt11_ee->Add(nt11_ee_ss);
		h_nt11_ee->Add(nt11_ee_wz);
		h_nt11_ee->Add(nt11_ee_cm);

		h_nt11_em->Add(nt11_em_df);
		h_nt11_em->Add(nt11_em_ss);
		h_nt11_em->Add(nt11_em_wz);
		h_nt11_em->Add(nt11_em_cm);
				
		// Signal
		// nt11_tot->Add(nt11_sig);
		// nt11_mm_tot->Add(nt11_mm_sig);
		// nt11_ee_tot->Add(nt11_ee_sig);
		// nt11_em_tot->Add(nt11_em_sig);
		
		// RATIOS
		rat   ->Divide(nt11   ,h_nt11   );
		rat_mm->Divide(nt11_mm,h_nt11_mm);
		rat_em->Divide(nt11_em,h_nt11_em);
		rat_ee->Divide(nt11_ee,h_nt11_ee);

		TString ytitle = Form("Events / %3.0f GeV", binwidthscale[j]);
		if(j==4 || j==8 || j==10) ytitle = "Events";
		nt11_tot->Draw("goff");
		if (plotRatios) for(size_t i = 1; i <= nt11->GetNbinsX(); ++i)     nt11_tot->GetXaxis()->SetBinLabel(i, "");
		nt11_tot->GetXaxis()->SetTitle(DiffPredYields::axis_label[j]);
		nt11_tot->GetYaxis()->SetTitle(ytitle);
		nt11_mm_tot->Draw("goff");
		if (plotRatios) for(size_t i = 1; i <= nt11_mm->GetNbinsX(); ++i)  nt11_mm_tot->GetXaxis()->SetBinLabel(i, "");
		nt11_mm_tot->GetXaxis()->SetTitle(DiffPredYields::axis_label[j]);
		nt11_mm_tot->GetYaxis()->SetTitle(ytitle);
		nt11_ee_tot->Draw("goff");
		if (plotRatios) for(size_t i = 1; i <= nt11_ee->GetNbinsX(); ++i)  nt11_ee_tot->GetXaxis()->SetBinLabel(i, "");
		nt11_ee_tot->GetXaxis()->SetTitle(DiffPredYields::axis_label[j]);
		nt11_ee_tot->GetYaxis()->SetTitle(ytitle);
		nt11_em_tot->Draw("goff");
		if (plotRatios) for(size_t i = 1; i <= nt11_em->GetNbinsX(); ++i)  nt11_em_tot->GetXaxis()->SetBinLabel(i, "");
		nt11_em_tot->GetXaxis()->SetTitle(DiffPredYields::axis_label[j]);
		nt11_em_tot->GetYaxis()->SetTitle(ytitle);

		nt11_tot   ->SetMinimum(0.5*nt11   ->GetMinimum()+0.001);
		nt11_mm_tot->SetMinimum(0.5*nt11_mm->GetMinimum()+0.001);
		nt11_ee_tot->SetMinimum(0.5*nt11_ee->GetMinimum()+0.001);
		nt11_em_tot->SetMinimum(0.5*nt11_em->GetMinimum()+0.001);

		double max = nt11->Integral();
		nt11    ->SetMaximum(max>1?max+1:1.);
		nt11_sf ->SetMaximum(max>1?max+1:1.);
		nt11_df ->SetMaximum(max>1?max+1:1.);
		nt11_cm ->SetMaximum(max>1?max+1:1.);
		nt11_ss ->SetMaximum(max>1?max+1:1.);
		nt11_wz ->SetMaximum(max>1?max+1:1.);
		nt11_tot->SetMaximum(max>1?max+1:1.);

		double max_mm = nt11_mm->Integral();
		nt11_mm    ->SetMaximum(max_mm>1?max_mm+1:1.);
		nt11_mm_sf ->SetMaximum(max_mm>1?max_mm+1:1.);
		nt11_mm_df ->SetMaximum(max_mm>1?max_mm+1:1.);
		nt11_mm_ss ->SetMaximum(max_mm>1?max_mm+1:1.);
		nt11_mm_wz ->SetMaximum(max_mm>1?max_mm+1:1.);
		nt11_mm_tot->SetMaximum(max_mm>1?max_mm+1:1.);

		double max_ee = nt11_ee->Integral();
		nt11_ee    ->SetMaximum(max_ee>1?max_ee+1:1.);
		nt11_ee_sf ->SetMaximum(max_ee>1?max_ee+1:1.);
		nt11_ee_df ->SetMaximum(max_ee>1?max_ee+1:1.);
		nt11_ee_cm ->SetMaximum(max_ee>1?max_ee+1:1.);
		nt11_ee_ss ->SetMaximum(max_ee>1?max_ee+1:1.);
		nt11_ee_wz ->SetMaximum(max_ee>1?max_ee+1:1.);
		nt11_ee_tot->SetMaximum(max_ee>1?max_ee+1:1.);

		double max_em = nt11_em->Integral();
		nt11_em    ->SetMaximum(max_em>1?max_em+1:1.);
		nt11_em_sf ->SetMaximum(max_em>1?max_em+1:1.);
		nt11_em_df ->SetMaximum(max_em>1?max_em+1:1.);
		nt11_em_cm ->SetMaximum(max_em>1?max_em+1:1.);
		nt11_em_ss ->SetMaximum(max_em>1?max_em+1:1.);
		nt11_em_wz ->SetMaximum(max_em>1?max_em+1:1.);
		nt11_em_tot->SetMaximum(max_em>1?max_em+1:1.);
		
		fOutputSubDir = "DiffPredictionPlots/";
		/////////////////////////////////////////////////////////////////
		TLegend *leg = new TLegend(0.60,0.65,0.90,0.88);
		leg->AddEntry(nt11,    "Observed","p");
		leg->AddEntry(nt11_sf, "Single Fakes","f");
		leg->AddEntry(nt11_df, "Double Fakes","f");
		leg->AddEntry(nt11_ss, "Irreducible (MC)","f");
		leg->AddEntry(nt11_wz, "WZ Production (MC)","f");
		leg->AddEntry(nt11_cm, "Charge MisID","f");
		leg->AddEntry(totbg,   "Total Uncertainty","f");
		// leg->AddEntry(nt11_sig,fSamples[sigsam]->sname,"l");
		leg->SetFillStyle(0);
		leg->SetTextFont(42);
		leg->SetBorderSize(0);
		
		TCanvas *c_temp = new TCanvas("C_ObsPred_" + varname, "Observed vs Predicted", 0, 0, 800, 800);
		c_temp->cd();
		
		TPad *p_plot  = new TPad("plotpad",  "Pad containing the plot", 0.00, border, 1.00, 1.00, 0, 0);
		p_plot->SetBottomMargin(0.015);
		if (plotRatios)	  p_plot->Draw();
		TPad *p_ratio = new TPad("ratiopad", "Pad containing the ratio", 0.00, 0.00, 1.00, border, 0, 0);
		p_ratio->SetTopMargin(0.025);
		p_ratio->SetBottomMargin(0.35);
		TLine *l3 = new TLine(nt11->GetXaxis()->GetXmin(), 1.00, nt11->GetXaxis()->GetXmax(), 1.00);
		if (plotRatios)	{
		  p_ratio->Draw();
		  p_ratio->cd();
		  rat->SetMaximum(1.8);
		  rat->SetMinimum(0.2);
		  rat->DrawCopy("E2");
		  l3->SetLineWidth(2);
		  l3->SetLineStyle(7);
		  l3->Draw();
		  gPad->RedrawAxis();
		  p_ratio->Draw();
		  p_plot->cd();
		}
		gPad->SetLogy();
		
		nt11_tot->Draw("hist");
		// nt11_error->DrawCopy("X0 E1 same");
		nt11->DrawCopy("PE X0 same");
		totbg->DrawCopy("0 E2 same");
		// nt11_sig->DrawCopy("hist same");
		leg->Draw();
		lat->SetTextSize(0.04);
		lat->DrawLatex(0.45,0.92, "#mu#mu/ee/e#mu");
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
		
		c_temp = new TCanvas("C_ObsPred_MM_" + varname, "Observed vs Predicted", 0, 0, 800, 800);
		c_temp->cd();
		
		if (plotRatios){
		  p_plot->Draw();
		  p_ratio->Draw();
		  p_ratio->cd();
		  rat_mm->SetMaximum(1.8);
		  rat_mm->SetMinimum(0.2);
		  rat_mm->DrawCopy("E2");
		  l3 = new TLine(nt11_mm->GetXaxis()->GetXmin(), 1.00, nt11_mm->GetXaxis()->GetXmax(), 1.00);
		  l3->SetLineWidth(2);
		  l3->SetLineStyle(7);
		  l3->Draw();
		  gPad->RedrawAxis();
		  p_ratio->Draw();
		  p_plot->cd();
		}
		gPad->SetLogy();
		
		nt11_mm_tot->Draw("hist");
		// nt11_error->DrawCopy("X0 E1 same");
		nt11_mm->DrawCopy("PE X0 same");
		totbg_mm->DrawCopy("0 E2 same");
		// nt11_mm_sig->DrawCopy("hist same");
		leg_mm->Draw();
		lat->SetTextSize(0.04);
		lat->DrawLatex(0.45,0.92, "#mu#mu");
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
		
		c_temp = new TCanvas("C_ObsPred_EE_" + varname, "Observed vs Predicted", 0, 0, 800, 800);
		c_temp->cd();

		if (plotRatios){
		  p_plot->Draw();
		  p_ratio->Draw();
		  p_ratio->cd();
		  rat_ee->SetMaximum(1.8);
		  rat_ee->SetMinimum(0.2);
		  rat_ee->DrawCopy("E2");
		  l3 = new TLine(nt11_ee->GetXaxis()->GetXmin(), 1.00, nt11_ee->GetXaxis()->GetXmax(), 1.00);
		  l3->SetLineWidth(2);
		  l3->SetLineStyle(7);
		  l3->Draw();
		  gPad->RedrawAxis();
		  p_ratio->Draw();
		  p_plot->cd();
		}
		gPad->SetLogy();
		
		nt11_ee_tot->Draw("hist");
		// nt11_error->DrawCopy("X0 E1 same");
		nt11_ee->DrawCopy("PE X0 same");
		totbg_ee->DrawCopy("0 E2 same");
		// nt11_ee_sig->DrawCopy("hist same");
		leg_ee->Draw();
		lat->SetTextSize(0.04);
		lat->DrawLatex(0.45,0.92, "ee");
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
		
		c_temp = new TCanvas("C_ObsPred_EM_" + varname, "Observed vs Predicted", 0, 0, 800, 800);
		c_temp->cd();
		if (plotRatios) {
		  p_plot->Draw();
		  p_ratio->Draw();
		  p_ratio->cd();
		  rat_em->SetMaximum(1.8);
		  rat_em->SetMinimum(0.2);
		  rat_em->DrawCopy("E2");
		  l3 = new TLine(nt11_em->GetXaxis()->GetXmin(), 1.00, nt11_em->GetXaxis()->GetXmax(), 1.00);
		  l3->SetLineWidth(2);
		  l3->SetLineStyle(7);
		  l3->Draw();
		  gPad->RedrawAxis();
		  p_ratio->Draw();
		  p_plot->cd();
		}
		gPad->SetLogy();
		  
		nt11_em_tot->Draw("hist");
		totbg_em->DrawCopy("0 E2 same");
		nt11_em->DrawCopy("PE X0 same");
		// nt11_em_sig->DrawCopy("hist same");
		leg_em->Draw();
		lat->SetTextSize(0.04);
		lat->DrawLatex(0.45,0.92, "e#mu");
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
		delete p_plot, p_ratio, c_temp, leg, leg_mm, leg_em, leg_ee;
		delete nt11, nt11_mm, nt11_ee, nt11_em;
		//		delete nt11_sig, nt11_mm_sig, nt11_ee_sig, nt11_em_sig;
		delete nt10_mm, nt10_em, nt10_ee, nt01_mm, nt01_em, nt01_ee, nt00_mm, nt00_em, nt00_ee;
		delete nt2_os_ee_bb, nt2_os_ee_eb, nt2_os_ee_ee, nt2_os_em_bb, nt2_os_em_ee;
		delete nt11_ss, nt11_mm_ss, nt11_em_ss, nt11_ee_ss;
		delete nt11_wz, nt11_mm_wz, nt11_em_wz, nt11_ee_wz;
		delete nt11_sf, nt11_mm_sf, nt11_em_sf, nt11_ee_sf;
		delete nt11_df, nt11_mm_df, nt11_em_df, nt11_ee_df;
		delete nt11_cm, nt11_em_cm, nt11_ee_cm;
		// delete nt11_mc, nt11_mm_mc, nt11_ee_mc, nt11_em_mc;
		delete nt11_tot, nt11_mm_tot, nt11_ee_tot, nt11_em_tot;
		delete totbg, totbg_mm, totbg_em, totbg_ee;
	}
	fOutputSubDir = "";
}

void SSDLPlotter::makeTTWDiffPredictions(){
	for(size_t i = 0; i < gNDiffVars; ++i){
		if(i == 7) continue;
		makeDiffPredictionTTW(i);
	}
}
void SSDLPlotter::makeDiffPredictionTTW(int varbin){
	fOutputSubDir = "DiffPredictionPlotsTTWZ/";
	TLatex *lat = new TLatex();
	lat->SetNDC(kTRUE);
	lat->SetTextColor(kBlack);
	lat->SetTextSize(0.04);

	bool separateTTH = false;

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

	// {  0 ,   1  ,    2   ,   3  ,   4  ,   5  ,    6    ,   7   ,      8     ,      9      }
	// {"HT", "MET", "NJets", "MT2", "PT1", "PT2", "NBJets", "MET3", "NBJetsMed", "NBJetsMed2"}
	float binwidthscale[gNDiffVars] = {100., 20., 1., 25., 20., 10., 1., 10., 1., 1.};

	// Loop on the different variables
	TString varname    = DiffPredYields::var_name[varbin];
	const int nbins    = DiffPredYields::nbins[varbin];
	const double *bins = DiffPredYields::bins[varbin];

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
		nt11_mm->Add(S->diffyields[Muon].hnt11[varbin]);
		nt10_mm->Add(S->diffyields[Muon].hnt10[varbin]);
		nt01_mm->Add(S->diffyields[Muon].hnt01[varbin]);
		nt00_mm->Add(S->diffyields[Muon].hnt00[varbin]);
	}
	for(size_t i = 0; i < elsamples.size(); ++i){
		Sample *S = fSamples[elsamples[i]];
		nt11_ee->Add(S->diffyields[Elec].hnt11[varbin]);
		nt10_ee->Add(S->diffyields[Elec].hnt10[varbin]);
		nt01_ee->Add(S->diffyields[Elec].hnt01[varbin]);
		nt00_ee->Add(S->diffyields[Elec].hnt00[varbin]);

		nt2_os_ee_bb->Add(S->diffyields[Elec].hnt2_os_BB[varbin]);
		nt2_os_ee_eb->Add(S->diffyields[Elec].hnt2_os_EB[varbin]);
		nt2_os_ee_ee->Add(S->diffyields[Elec].hnt2_os_EE[varbin]);
	}
	for(size_t i = 0; i < emusamples.size(); ++i){
		Sample *S = fSamples[emusamples[i]];
		nt11_em->Add(S->diffyields[ElMu].hnt11[varbin]);
		nt10_em->Add(S->diffyields[ElMu].hnt10[varbin]);
		nt01_em->Add(S->diffyields[ElMu].hnt01[varbin]);
		nt00_em->Add(S->diffyields[ElMu].hnt00[varbin]);

		nt2_os_em_bb->Add(S->diffyields[ElMu].hnt2_os_BB[varbin]);
		nt2_os_em_ee->Add(S->diffyields[ElMu].hnt2_os_EE[varbin]);
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
	// RARE SM MC /////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////
	TH1D *nt11_ss = new TH1D(Form("NT11_SS_%s", varname.Data()), varname, nbins, bins); nt11_ss->Sumw2();

	TH1D *nt11_mm_ss = new TH1D(Form("NT11_MM_SS_%s", varname.Data()), varname, nbins, bins); nt11_mm_ss->Sumw2();
	TH1D *nt11_ee_ss = new TH1D(Form("NT11_EE_SS_%s", varname.Data()), varname, nbins, bins); nt11_ee_ss->Sumw2();
	TH1D *nt11_em_ss = new TH1D(Form("NT11_EM_SS_%s", varname.Data()), varname, nbins, bins); nt11_em_ss->Sumw2();

	vector<int> raremc;
//	raremc.push_back(WZ);
	raremc.push_back(ZZ);
	// MARC raremc.push_back(GVJets);
	raremc.push_back(DPSWW);
	raremc.push_back(TTbarG);
	if (!separateTTH) raremc.push_back(TTbarH);
	raremc.push_back(WpWp);
	raremc.push_back(WmWm);
	raremc.push_back(WWZ);
	raremc.push_back(WZZ);
	raremc.push_back(WWG);
	raremc.push_back(WWW);
	raremc.push_back(ZZZ);

	for(size_t i = 0; i < raremc.size(); ++i){
		Sample *S = fSamples[raremc[i]];
		float scale = fLumiNorm / S->getLumi();
		nt11_mm_ss->Add(S->diffyields[Muon].hnt11[varbin], scale);
		nt11_ee_ss->Add(S->diffyields[Elec].hnt11[varbin], scale);
		nt11_em_ss->Add(S->diffyields[ElMu].hnt11[varbin], scale);

		// Errors
		for(size_t b = 0; b < nbins; ++b){
			float ss_mm = S->diffyields[Muon].hnt11[varbin]->GetBinContent(b+1);
			float ss_ee = S->diffyields[Elec].hnt11[varbin]->GetBinContent(b+1);
			float ss_em = S->diffyields[ElMu].hnt11[varbin]->GetBinContent(b+1);

			float esyst2_mm = 0.25 * ss_mm*ss_mm*scale*scale;
			float esyst2_ee = 0.25 * ss_ee*ss_ee*scale*scale;
			float esyst2_em = 0.25 * ss_em*ss_em*scale*scale;

			float estat2_mm = scale*scale*S->getError2(ss_mm);
			float estat2_ee = scale*scale*S->getError2(ss_ee);
			float estat2_em = scale*scale*S->getError2(ss_em);
			// float estat2_mm = scale*scale*S->numbers[reg][Muon].tt_avweight*S->numbers[reg][Muon].tt_avweight*S->getError2(ss_mm);
			// float estat2_ee = scale*scale*S->numbers[reg][Elec].tt_avweight*S->numbers[reg][Elec].tt_avweight*S->getError2(ss_ee);
			// float estat2_em = scale*scale*S->numbers[reg][ElMu].tt_avweight*S->numbers[reg][ElMu].tt_avweight*S->getError2(ss_em);

			float prev    = totbg   ->GetBinError(b+1);
			float prev_mm = totbg_mm->GetBinError(b+1);
			float prev_em = totbg_em->GetBinError(b+1);
			float prev_ee = totbg_ee->GetBinError(b+1);

			totbg   ->SetBinError(b+1, prev    + esyst2_mm + esyst2_ee + esyst2_em + estat2_mm + estat2_ee + estat2_em);
			totbg_mm->SetBinError(b+1, prev_mm + esyst2_mm + estat2_mm);
			totbg_em->SetBinError(b+1, prev_em + esyst2_em + estat2_em);
			totbg_ee->SetBinError(b+1, prev_ee + esyst2_ee + estat2_ee);
		}
	}

	nt11_ss->Add(nt11_mm_ss);
	nt11_ss->Add(nt11_ee_ss);
	nt11_ss->Add(nt11_em_ss);

	totbg   ->Add(nt11_ss);
	totbg_mm->Add(nt11_mm_ss);
	totbg_em->Add(nt11_em_ss);
	totbg_ee->Add(nt11_ee_ss);

	///////////////////////////////////////////////////////////////////////////////////
	// WZ PRODUCTION //////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////
	TH1D *nt11_wz = new TH1D(Form("NT11_WZ_%s", varname.Data()), varname, nbins, bins); nt11_wz->Sumw2();

	TH1D *nt11_mm_wz = new TH1D(Form("NT11_MM_WZ_%s", varname.Data()), varname, nbins, bins); nt11_mm_wz->Sumw2();
	TH1D *nt11_ee_wz = new TH1D(Form("NT11_EE_WZ_%s", varname.Data()), varname, nbins, bins); nt11_ee_wz->Sumw2();
	TH1D *nt11_em_wz = new TH1D(Form("NT11_EM_WZ_%s", varname.Data()), varname, nbins, bins); nt11_em_wz->Sumw2();

	float wzscale = fLumiNorm / fSamples[WZ]->getLumi();
	nt11_mm_wz->Add(fSamples[WZ]->diffyields[Muon].hnt11[varbin], wzscale);
	nt11_ee_wz->Add(fSamples[WZ]->diffyields[Elec].hnt11[varbin], wzscale);
	nt11_em_wz->Add(fSamples[WZ]->diffyields[ElMu].hnt11[varbin], wzscale);

	// Errors
	for(size_t b = 0; b < nbins; ++b){
		float ss_mm = fSamples[WZ]->diffyields[Muon].hnt11[varbin]->GetBinContent(b+1);
		float ss_ee = fSamples[WZ]->diffyields[Elec].hnt11[varbin]->GetBinContent(b+1);
		float ss_em = fSamples[WZ]->diffyields[ElMu].hnt11[varbin]->GetBinContent(b+1);

		float esyst2_mm = 0.25 * ss_mm*ss_mm*wzscale*wzscale;
		float esyst2_ee = 0.25 * ss_ee*ss_ee*wzscale*wzscale;
		float esyst2_em = 0.25 * ss_em*ss_em*wzscale*wzscale;

		float estat2_mm = wzscale*wzscale*fSamples[WZ]->getError2(ss_mm);
		float estat2_ee = wzscale*wzscale*fSamples[WZ]->getError2(ss_ee);
		float estat2_em = wzscale*wzscale*fSamples[WZ]->getError2(ss_em);
		// float estat2_mm = wzscale*wzscale*fSamples[WZ]->numbers[reg][Muon].tt_avweight*fSamples[WZ]->numbers[reg][Muon].tt_avweight*fSamples[WZ]->getError2(ss_mm);
		// float estat2_ee = wzscale*wzscale*fSamples[WZ]->numbers[reg][Elec].tt_avweight*fSamples[WZ]->numbers[reg][Elec].tt_avweight*fSamples[WZ]->getError2(ss_ee);
		// float estat2_em = wzscale*wzscale*fSamples[WZ]->numbers[reg][ElMu].tt_avweight*fSamples[WZ]->numbers[reg][ElMu].tt_avweight*fSamples[WZ]->getError2(ss_em);

		float prev    = totbg   ->GetBinError(b+1);
		float prev_mm = totbg_mm->GetBinError(b+1);
		float prev_em = totbg_em->GetBinError(b+1);
		float prev_ee = totbg_ee->GetBinError(b+1);

		totbg   ->SetBinError(b+1, prev    + esyst2_mm + esyst2_ee + esyst2_em + estat2_mm + estat2_ee + estat2_em);
		totbg_mm->SetBinError(b+1, prev_mm + esyst2_mm + estat2_mm);
		totbg_em->SetBinError(b+1, prev_em + esyst2_em + estat2_em);
		totbg_ee->SetBinError(b+1, prev_ee + esyst2_ee + estat2_ee);
	}

	nt11_wz->Add(nt11_mm_wz);
	nt11_wz->Add(nt11_ee_wz);
	nt11_wz->Add(nt11_em_wz);

	totbg   ->Add(nt11_wz);
	totbg_mm->Add(nt11_mm_wz);
	totbg_em->Add(nt11_em_wz);
	totbg_ee->Add(nt11_ee_wz);

	///////////////////////////////////////////////////////////////////////////////////
	// TTW PRODUCTION /////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////
	TH1D *nt11_ttw = new TH1D(Form("NT11_TTW_%s", varname.Data()), varname, nbins, bins); nt11_ttw->Sumw2();

	TH1D *nt11_mm_ttw = new TH1D(Form("NT11_MM_TTW_%s", varname.Data()), varname, nbins, bins); nt11_mm_ttw->Sumw2();
	TH1D *nt11_ee_ttw = new TH1D(Form("NT11_EE_TTW_%s", varname.Data()), varname, nbins, bins); nt11_ee_ttw->Sumw2();
	TH1D *nt11_em_ttw = new TH1D(Form("NT11_EM_TTW_%s", varname.Data()), varname, nbins, bins); nt11_em_ttw->Sumw2();

	float ttwscale = fLumiNorm / fSamples[TTbarW]->getLumi();
	nt11_mm_ttw->Add(fSamples[TTbarW]->diffyields[Muon].hnt11[varbin], ttwscale);
	nt11_ee_ttw->Add(fSamples[TTbarW]->diffyields[Elec].hnt11[varbin], ttwscale);
	nt11_em_ttw->Add(fSamples[TTbarW]->diffyields[ElMu].hnt11[varbin], ttwscale);

	// Errors
	for(size_t b = 0; b < nbins; ++b){
		float ss_mm = fSamples[TTbarW]->diffyields[Muon].hnt11[varbin]->GetBinContent(b+1);
		float ss_ee = fSamples[TTbarW]->diffyields[Elec].hnt11[varbin]->GetBinContent(b+1);
		float ss_em = fSamples[TTbarW]->diffyields[ElMu].hnt11[varbin]->GetBinContent(b+1);

		float esyst2_mm = 0.25 * ss_mm*ss_mm*ttwscale*ttwscale;
		float esyst2_ee = 0.25 * ss_ee*ss_ee*ttwscale*ttwscale;
		float esyst2_em = 0.25 * ss_em*ss_em*ttwscale*ttwscale;

		float estat2_mm = ttwscale*ttwscale*fSamples[TTbarW]->getError2(ss_mm);
		float estat2_ee = ttwscale*ttwscale*fSamples[TTbarW]->getError2(ss_ee);
		float estat2_em = ttwscale*ttwscale*fSamples[TTbarW]->getError2(ss_em);
		// float estat2_mm = ttwscale*ttwscale*fSamples[TTbarW]->numbers[reg][Muon].tt_avweight*fSamples[TTbarW]->numbers[reg][Muon].tt_avweight*fSamples[TTbarW]->getError2(ss_mm);
		// float estat2_ee = ttwscale*ttwscale*fSamples[TTbarW]->numbers[reg][Elec].tt_avweight*fSamples[TTbarW]->numbers[reg][Elec].tt_avweight*fSamples[TTbarW]->getError2(ss_ee);
		// float estat2_em = ttwscale*ttwscale*fSamples[TTbarW]->numbers[reg][ElMu].tt_avweight*fSamples[TTbarW]->numbers[reg][ElMu].tt_avweight*fSamples[TTbarW]->getError2(ss_em);

		float prev    = totbg   ->GetBinError(b+1);
		float prev_mm = totbg_mm->GetBinError(b+1);
		float prev_em = totbg_em->GetBinError(b+1);
		float prev_ee = totbg_ee->GetBinError(b+1);

		totbg   ->SetBinError(b+1, prev    + esyst2_mm + esyst2_ee + esyst2_em + estat2_mm + estat2_ee + estat2_em);
		totbg_mm->SetBinError(b+1, prev_mm + esyst2_mm + estat2_mm);
		totbg_em->SetBinError(b+1, prev_em + esyst2_em + estat2_em);
		totbg_ee->SetBinError(b+1, prev_ee + esyst2_ee + estat2_ee);
	}

	nt11_ttw->Add(nt11_mm_ttw);
	nt11_ttw->Add(nt11_ee_ttw);
	nt11_ttw->Add(nt11_em_ttw);

	// totbg   ->Add(nt11_ttw);
	// totbg_mm->Add(nt11_mm_ttw);
	// totbg_em->Add(nt11_em_ttw);
	// totbg_ee->Add(nt11_ee_ttw);


	///////////////////////////////////////////////////////////////////////////////////
	// TTZ PRODUCTION /////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////
	TH1D *nt11_ttz = new TH1D(Form("NT11_TTZ_%s", varname.Data()), varname, nbins, bins); nt11_ttz->Sumw2();

	TH1D *nt11_mm_ttz = new TH1D(Form("NT11_MM_TTZ_%s", varname.Data()), varname, nbins, bins); nt11_mm_ttz->Sumw2();
	TH1D *nt11_ee_ttz = new TH1D(Form("NT11_EE_TTZ_%s", varname.Data()), varname, nbins, bins); nt11_ee_ttz->Sumw2();
	TH1D *nt11_em_ttz = new TH1D(Form("NT11_EM_TTZ_%s", varname.Data()), varname, nbins, bins); nt11_em_ttz->Sumw2();

	float ttzscale = fLumiNorm / fSamples[TTbarZ]->getLumi();
	nt11_mm_ttz->Add(fSamples[TTbarZ]->diffyields[Muon].hnt11[varbin], ttzscale);
	nt11_ee_ttz->Add(fSamples[TTbarZ]->diffyields[Elec].hnt11[varbin], ttzscale);
	nt11_em_ttz->Add(fSamples[TTbarZ]->diffyields[ElMu].hnt11[varbin], ttzscale);

	// Errors
	for(size_t b = 0; b < nbins; ++b){
		float ss_mm = fSamples[TTbarZ]->diffyields[Muon].hnt11[varbin]->GetBinContent(b+1);
		float ss_ee = fSamples[TTbarZ]->diffyields[Elec].hnt11[varbin]->GetBinContent(b+1);
		float ss_em = fSamples[TTbarZ]->diffyields[ElMu].hnt11[varbin]->GetBinContent(b+1);

		float esyst2_mm = 0.25 * ss_mm*ss_mm*ttzscale*ttzscale;
		float esyst2_ee = 0.25 * ss_ee*ss_ee*ttzscale*ttzscale;
		float esyst2_em = 0.25 * ss_em*ss_em*ttzscale*ttzscale;

		float estat2_mm = ttzscale*ttzscale*fSamples[TTbarZ]->getError2(ss_mm);
		float estat2_ee = ttzscale*ttzscale*fSamples[TTbarZ]->getError2(ss_ee);
		float estat2_em = ttzscale*ttzscale*fSamples[TTbarZ]->getError2(ss_em);
		// float estat2_mm = ttzscale*ttzscale*fSamples[TTbarZ]->numbers[reg][Muon].tt_avweight*fSamples[TTbarZ]->numbers[reg][Muon].tt_avweight*fSamples[TTbarZ]->getError2(ss_mm);
		// float estat2_ee = ttzscale*ttzscale*fSamples[TTbarZ]->numbers[reg][Elec].tt_avweight*fSamples[TTbarZ]->numbers[reg][Elec].tt_avweight*fSamples[TTbarZ]->getError2(ss_ee);
		// float estat2_em = ttzscale*ttzscale*fSamples[TTbarZ]->numbers[reg][ElMu].tt_avweight*fSamples[TTbarZ]->numbers[reg][ElMu].tt_avweight*fSamples[TTbarZ]->getError2(ss_em);

		float prev    = totbg   ->GetBinError(b+1);
		float prev_mm = totbg_mm->GetBinError(b+1);
		float prev_em = totbg_em->GetBinError(b+1);
		float prev_ee = totbg_ee->GetBinError(b+1);

		totbg   ->SetBinError(b+1, prev    + esyst2_mm + esyst2_ee + esyst2_em + estat2_mm + estat2_ee + estat2_em);
		totbg_mm->SetBinError(b+1, prev_mm + esyst2_mm + estat2_mm);
		totbg_em->SetBinError(b+1, prev_em + esyst2_em + estat2_em);
		totbg_ee->SetBinError(b+1, prev_ee + esyst2_ee + estat2_ee);
	}

	nt11_ttz->Add(nt11_mm_ttz);
	nt11_ttz->Add(nt11_ee_ttz);
	nt11_ttz->Add(nt11_em_ttz);

	// totbg   ->Add(nt11_ttz);
	// totbg_mm->Add(nt11_mm_ttz);
	// totbg_em->Add(nt11_em_ttz);
	// totbg_ee->Add(nt11_ee_ttz);

	///////////////////////////////////////////////////////////////////////////////////
	// TTH PRODUCTION /////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////
	TH1D *nt11_tth = new TH1D(Form("NT11_tth_%s", varname.Data()), varname, nbins, bins); nt11_tth->Sumw2();
	
	TH1D *nt11_mm_tth = new TH1D(Form("NT11_MM_tth_%s", varname.Data()), varname, nbins, bins); nt11_mm_tth->Sumw2();
	TH1D *nt11_ee_tth = new TH1D(Form("NT11_EE_tth_%s", varname.Data()), varname, nbins, bins); nt11_ee_tth->Sumw2();
	TH1D *nt11_em_tth = new TH1D(Form("NT11_EM_tth_%s", varname.Data()), varname, nbins, bins); nt11_em_tth->Sumw2();
	
	if (separateTTH) {
		float tthscale = fLumiNorm / fSamples[TTbarH]->getLumi();
		nt11_mm_tth->Add(fSamples[TTbarH]->diffyields[Muon].hnt11[varbin], tthscale);
		nt11_ee_tth->Add(fSamples[TTbarH]->diffyields[Elec].hnt11[varbin], tthscale);
		nt11_em_tth->Add(fSamples[TTbarH]->diffyields[ElMu].hnt11[varbin], tthscale);
	
		// Errors
		for(size_t b = 0; b < nbins; ++b){
			float ss_mm = fSamples[TTbarH]->diffyields[Muon].hnt11[varbin]->GetBinContent(b+1);
			float ss_ee = fSamples[TTbarH]->diffyields[Elec].hnt11[varbin]->GetBinContent(b+1);
			float ss_em = fSamples[TTbarH]->diffyields[ElMu].hnt11[varbin]->GetBinContent(b+1);
	
			float esyst2_mm = 0.25 * ss_mm*ss_mm*tthscale*tthscale;
			float esyst2_ee = 0.25 * ss_ee*ss_ee*tthscale*tthscale;
			float esyst2_em = 0.25 * ss_em*ss_em*tthscale*tthscale;
	
			float estat2_mm = tthscale*tthscale*fSamples[TTbarH]->getError2(ss_mm);
			float estat2_ee = tthscale*tthscale*fSamples[TTbarH]->getError2(ss_ee);
			float estat2_em = tthscale*tthscale*fSamples[TTbarH]->getError2(ss_em);
			// float estat2_mm = tthscale*tthscale*fSamples[TTbarH]->numbers[reg][Muon].tt_avweight*fSamples[TTbarH]->numbers[reg][Muon].tt_avweight*fSamples[TTbarH]->getError2(ss_mm);
			// float estat2_ee = tthscale*tthscale*fSamples[TTbarH]->numbers[reg][Elec].tt_avweight*fSamples[TTbarH]->numbers[reg][Elec].tt_avweight*fSamples[TTbarH]->getError2(ss_ee);
			// float estat2_em = tthscale*tthscale*fSamples[TTbarH]->numbers[reg][ElMu].tt_avweight*fSamples[TTbarH]->numbers[reg][ElMu].tt_avweight*fSamples[TTbarH]->getError2(ss_em);
	
			float prev    = totbg   ->GetBinError(b+1);
			float prev_mm = totbg_mm->GetBinError(b+1);
			float prev_em = totbg_em->GetBinError(b+1);
			float prev_ee = totbg_ee->GetBinError(b+1);
	
			totbg   ->SetBinError(b+1, prev    + esyst2_mm + esyst2_ee + esyst2_em + estat2_mm + estat2_ee + estat2_em);
			totbg_mm->SetBinError(b+1, prev_mm + esyst2_mm + estat2_mm);
			totbg_em->SetBinError(b+1, prev_em + esyst2_em + estat2_em);
			totbg_ee->SetBinError(b+1, prev_ee + esyst2_ee + estat2_ee);
		}
	
		nt11_tth->Add(nt11_mm_tth);
		nt11_tth->Add(nt11_ee_tth);
		nt11_tth->Add(nt11_em_tth);

	} // end if seprarate TTH


	// totbg   ->Add(nt11_tth);
	// totbg_mm->Add(nt11_mm_tth);
	// totbg_em->Add(nt11_em_tth);
	// totbg_ee->Add(nt11_ee_tth);

	///////////////////////////////////////////////////////////////////////////////////
	// FAKE PREDICTIONS ///////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////
	TH1D *nt11_sf = new TH1D(Form("NT11_SF_%s", varname.Data()), varname, nbins, bins); nt11_sf->Sumw2();
	TH1D *nt11_mm_sf = new TH1D(Form("NT11_MM_SF_%s", varname.Data()), varname, nbins, bins); nt11_mm_sf->Sumw2();
	TH1D *nt11_ee_sf = new TH1D(Form("NT11_EE_SF_%s", varname.Data()), varname, nbins, bins); nt11_ee_sf->Sumw2();
	TH1D *nt11_em_sf = new TH1D(Form("NT11_EM_SF_%s", varname.Data()), varname, nbins, bins); nt11_em_sf->Sumw2();

	/////////////////////////////////////////////////////
	// Differential ratios
	for(size_t i = 0; i < musamples.size(); ++i){
		Sample *S = fSamples[musamples[i]];
		nt11_mm_sf->Add(S->diffyields[Muon].hnpf[varbin]);
		nt11_mm_sf->Add(S->diffyields[Muon].hnfp[varbin]);
		nt11_mm_sf->Add(S->diffyields[Muon].hnff[varbin]);
	}
	for(size_t i = 0; i < elsamples.size(); ++i){
		Sample *S = fSamples[elsamples[i]];
		nt11_ee_sf->Add(S->diffyields[Elec].hnpf[varbin]);
		nt11_ee_sf->Add(S->diffyields[Elec].hnfp[varbin]);
		nt11_ee_sf->Add(S->diffyields[Elec].hnff[varbin]);
	}
	for(size_t i = 0; i < emusamples.size(); ++i){
		Sample *S = fSamples[emusamples[i]];
		nt11_em_sf->Add(S->diffyields[ElMu].hnpf[varbin]);
		nt11_em_sf->Add(S->diffyields[ElMu].hnfp[varbin]);
		nt11_em_sf->Add(S->diffyields[ElMu].hnff[varbin]);
	}
	/////////////////////////////////////////////////////

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
		
		float mm_tot_fakes = nt11_mm_sf->GetBinContent(i+1);
		float ee_tot_fakes = nt11_ee_sf->GetBinContent(i+1);
		float em_tot_fakes = nt11_em_sf->GetBinContent(i+1);
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

	totbg   ->Add(nt11_sf);
	totbg_mm->Add(nt11_mm_sf);
	totbg_em->Add(nt11_em_sf);
	totbg_ee->Add(nt11_ee_sf);

	///////////////////////////////////////////////////////////////////////////////////
	// E-CHARGE MISID /////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////
	TH1D *nt11_cm = new TH1D(Form("NT11_CM_%s", varname.Data()), varname, nbins, bins); nt11_cm->Sumw2();

	TH1D *nt11_ee_cm = new TH1D(Form("NT11_EE_CM_%s", varname.Data()), varname, nbins, bins); nt11_ee_cm->Sumw2();
	TH1D *nt11_em_cm = new TH1D(Form("NT11_EM_CM_%s", varname.Data()), varname, nbins, bins); nt11_em_cm->Sumw2();

	// Abbreviations
	float fbb(0.),fee(0.),feb(0.);
	float fbbE(0.),feeE(0.),febE(0.);
	float fbb_mc(0.),fee_mc(0.),feb_mc(0.);
	float fbbE_mc(0.),feeE_mc(0.),febE_mc(0.);
	
	calculateChMisIdProb(fEGData, BB, fbb, fbbE);
	calculateChMisIdProb(fEGData, EB, feb, febE);
	calculateChMisIdProb(fEGData, EE, fee, feeE);
	
	calculateChMisIdProb(fMCBG, BB, fbb_mc, fbbE_mc);
	calculateChMisIdProb(fMCBG, EB, feb_mc, febE_mc);
	calculateChMisIdProb(fMCBG, EE, fee_mc, feeE_mc);

	for(size_t i = 0; i < nbins; ++i){
		float nt2_ee_BB_os = nt2_os_ee_bb->GetBinContent(i+1);
		float nt2_ee_EB_os = nt2_os_ee_eb->GetBinContent(i+1);
		float nt2_ee_EE_os = nt2_os_ee_ee->GetBinContent(i+1);
		float nt2_em_BB_os = nt2_os_em_bb->GetBinContent(i+1);
		float nt2_em_EE_os = nt2_os_em_ee->GetBinContent(i+1);
		
		// Errors
		FakeRatios *FR = new FakeRatios();

		// Simple error propagation assuming error on number of events is FR->getEStat2()
		nt11_ee_cm->SetBinContent(i+1, 2*fbb*nt2_ee_BB_os + 2*fee*nt2_ee_EE_os + 2*feb*nt2_ee_EB_os);
		float nt11_ee_cm_e1 = sqrt( (4*fbb*fbb*FR->getEStat2(nt2_ee_BB_os)) + (4*fee*fee*FR->getEStat2(nt2_ee_EE_os)) + 4*(feb)*(feb)*FR->getEStat2(nt2_ee_EB_os) ); // stat only
		float nt11_ee_cm_e2 = sqrt( (4*nt2_ee_BB_os*nt2_ee_BB_os*fbbE*fbbE) + (4*nt2_ee_EE_os*nt2_ee_EE_os*feeE*feeE) + 4*(febE*febE)*nt2_ee_EB_os*nt2_ee_EB_os ); // syst only

		nt11_em_cm->SetBinContent(i+i, fbb*nt2_em_BB_os + fee*nt2_em_EE_os);
		float nt11_em_cm_e1 = sqrt( fbb*fbb*FR->getEStat2(nt2_em_BB_os) + fee*fee*FR->getEStat2(nt2_em_EE_os) );
		float nt11_em_cm_e2 = sqrt( nt2_em_BB_os*nt2_em_BB_os * fbbE*fbbE + nt2_em_EE_os*nt2_em_EE_os * feeE*feeE );
		
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
	// gSample sigsam = LM11;
	// TH1D *nt11_sig    = new TH1D(Form("NT11_Sig%s",    varname.Data()), varname, nbins, bins); nt11_sig   ->Sumw2();
	// TH1D *nt11_mm_sig = new TH1D(Form("NT11_mm_Sig%s", varname.Data()), varname, nbins, bins); nt11_mm_sig->Sumw2();
	// TH1D *nt11_em_sig = new TH1D(Form("NT11_em_Sig%s", varname.Data()), varname, nbins, bins); nt11_em_sig->Sumw2();
	// TH1D *nt11_ee_sig = new TH1D(Form("NT11_ee_Sig%s", varname.Data()), varname, nbins, bins); nt11_ee_sig->Sumw2();
	// nt11_mm_sig->Add(fSamples[sigsam]->diffyields[Muon].hnt11[varbin], fLumiNorm / fSamples[LM4]->getLumi());
	// nt11_ee_sig->Add(fSamples[sigsam]->diffyields[Elec].hnt11[varbin], fLumiNorm / fSamples[LM4]->getLumi());
	// nt11_em_sig->Add(fSamples[sigsam]->diffyields[ElMu].hnt11[varbin], fLumiNorm / fSamples[LM4]->getLumi());
	// nt11_sig   ->Add(fSamples[sigsam]->diffyields[Muon].hnt11[varbin], fLumiNorm / fSamples[LM4]->getLumi());
	// nt11_sig   ->Add(fSamples[sigsam]->diffyields[Elec].hnt11[varbin], fLumiNorm / fSamples[LM4]->getLumi());
	// nt11_sig   ->Add(fSamples[sigsam]->diffyields[ElMu].hnt11[varbin], fLumiNorm / fSamples[LM4]->getLumi());

	///////////////////////////////////////////////////////////////////////////////////
	// OUTPUT /////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////
	Color_t col_fake = 46;
	Color_t col_cmid = 49;
	Color_t col_wz   = 39;
	Color_t col_rare = 38;
	Color_t col_ttw  = 44;
	Color_t col_ttz  = 42;
	Color_t col_tth  = 45;
	
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
	nt11_cm->SetLineWidth(1);
	nt11_ss->SetLineWidth(1);
	nt11_wz->SetLineWidth(1);
	nt11_ttw->SetLineWidth(1);
	nt11_ttz->SetLineWidth(1);
	nt11_sf ->SetFillColor(col_fake);
	nt11_cm ->SetFillColor(col_cmid);
	nt11_ss ->SetFillColor(col_rare);
	nt11_wz ->SetFillColor(col_wz);
	nt11_ttw->SetFillColor(col_ttw);
	nt11_ttz->SetFillColor(col_ttz);
	if (separateTTH) nt11_tth->SetLineWidth(1);
	if (separateTTH) nt11_tth->SetFillColor(col_tth);

	nt11_mm_sf->SetLineWidth(1);
	nt11_mm_ss->SetLineWidth(1);
	nt11_mm_wz->SetLineWidth(1);
	nt11_mm_ttw->SetLineWidth(1);
	nt11_mm_ttz->SetLineWidth(1);
	if (separateTTH) nt11_mm_tth->SetLineWidth(1);
	nt11_mm_sf ->SetFillColor(col_fake);
	nt11_mm_ss ->SetFillColor(col_rare);
	nt11_mm_wz ->SetFillColor(col_wz);
	nt11_mm_ttw->SetFillColor(col_ttw);
	nt11_mm_ttz->SetFillColor(col_ttz);
	if (separateTTH) nt11_mm_tth->SetFillColor(col_tth);

	nt11_ee_sf->SetLineWidth(1);
	nt11_ee_cm->SetLineWidth(1);
	nt11_ee_ss->SetLineWidth(1);
	nt11_ee_wz->SetLineWidth(1);
	nt11_ee_ttw->SetLineWidth(1);
	nt11_ee_ttz->SetLineWidth(1);
	if (separateTTH) nt11_ee_tth->SetLineWidth(1);
	nt11_ee_sf ->SetFillColor(col_fake);
	nt11_ee_cm ->SetFillColor(col_cmid);
	nt11_ee_ss ->SetFillColor(col_rare);
	nt11_ee_wz ->SetFillColor(col_wz);
	nt11_ee_ttw->SetFillColor(col_ttw);
	nt11_ee_ttz->SetFillColor(col_ttz);
	if (separateTTH) nt11_ee_tth->SetFillColor(col_tth);

	nt11_em_sf->SetLineWidth(1);
	nt11_em_cm->SetLineWidth(1);
	nt11_em_ss->SetLineWidth(1);
	nt11_em_wz->SetLineWidth(1);
	nt11_em_ttw->SetLineWidth(1);
	nt11_em_ttz->SetLineWidth(1);
	if (separateTTH) nt11_em_tth->SetLineWidth(1);
	nt11_em_sf ->SetFillColor(col_fake);
	nt11_em_cm ->SetFillColor(col_cmid);
	nt11_em_ss ->SetFillColor(col_rare);
	nt11_em_wz ->SetFillColor(col_wz);
	nt11_em_ttw->SetFillColor(col_ttw);
	nt11_em_ttz->SetFillColor(col_ttz);
	if (separateTTH) nt11_em_tth->SetFillColor(col_tth);
	
	totbg   ->SetLineWidth(3);
	totbg_mm->SetLineWidth(3);
	totbg_em->SetLineWidth(3);
	totbg_ee->SetLineWidth(3);

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
	nt11_sf    = normHistBW(nt11_sf,  binwidthscale[varbin]);
	nt11_ss    = normHistBW(nt11_ss,  binwidthscale[varbin]);
	nt11_wz    = normHistBW(nt11_wz,  binwidthscale[varbin]);
	nt11_ttw   = normHistBW(nt11_ttw, binwidthscale[varbin]);
	nt11_ttz   = normHistBW(nt11_ttz, binwidthscale[varbin]);
	if (separateTTH) nt11_tth   = normHistBW(nt11_tth, binwidthscale[varbin]);
	nt11_cm    = normHistBW(nt11_cm,  binwidthscale[varbin]);

	nt11_mm_sf  = normHistBW(nt11_mm_sf, binwidthscale[varbin]);
	nt11_mm_ss  = normHistBW(nt11_mm_ss, binwidthscale[varbin]);
	nt11_mm_wz  = normHistBW(nt11_mm_wz, binwidthscale[varbin]);
	nt11_mm_ttw = normHistBW(nt11_mm_ttw, binwidthscale[varbin]);
	nt11_mm_ttz = normHistBW(nt11_mm_ttz, binwidthscale[varbin]);
	if (separateTTH) nt11_mm_tth = normHistBW(nt11_mm_tth, binwidthscale[varbin]);

	nt11_ee_sf  = normHistBW(nt11_ee_sf, binwidthscale[varbin]);
	nt11_ee_ss  = normHistBW(nt11_ee_ss, binwidthscale[varbin]);
	nt11_ee_wz  = normHistBW(nt11_ee_wz, binwidthscale[varbin]);
	nt11_ee_ttw = normHistBW(nt11_ee_ttw, binwidthscale[varbin]);
	nt11_ee_ttz = normHistBW(nt11_ee_ttz, binwidthscale[varbin]);
	if (separateTTH) nt11_ee_tth = normHistBW(nt11_ee_tth, binwidthscale[varbin]);
	nt11_ee_cm  = normHistBW(nt11_ee_cm, binwidthscale[varbin]);
	
	nt11_em_sf  = normHistBW(nt11_em_sf, binwidthscale[varbin]);
	nt11_em_ss  = normHistBW(nt11_em_ss, binwidthscale[varbin]);
	nt11_em_ttw = normHistBW(nt11_em_ttw, binwidthscale[varbin]);
	nt11_em_ttz = normHistBW(nt11_em_ttz, binwidthscale[varbin]);
	if (separateTTH) nt11_em_tth = normHistBW(nt11_em_tth, binwidthscale[varbin]);
	nt11_em_wz  = normHistBW(nt11_em_wz, binwidthscale[varbin]);
	nt11_em_cm  = normHistBW(nt11_em_cm, binwidthscale[varbin]);
	
	totbg      = normHistBW(totbg,    binwidthscale[varbin]);
	totbg_mm   = normHistBW(totbg_mm, binwidthscale[varbin]);
	totbg_em   = normHistBW(totbg_em, binwidthscale[varbin]);
	totbg_ee   = normHistBW(totbg_ee, binwidthscale[varbin]);

	nt11       = normHistBW(nt11,    binwidthscale[varbin]);
	nt11_mm    = normHistBW(nt11_mm, binwidthscale[varbin]);
	nt11_em    = normHistBW(nt11_em, binwidthscale[varbin]);
	nt11_ee    = normHistBW(nt11_ee, binwidthscale[varbin]);

	// Fill stacks
	THStack *nt11_tot    = new THStack("NT11_TotalBG", "NT11_TotalBG");
	THStack *nt11_mm_tot = new THStack("NT11_MM_TotalBG", "NT11_MM_TotalBG");
	THStack *nt11_ee_tot = new THStack("NT11_EE_TotalBG", "NT11_EE_TotalBG");
	THStack *nt11_em_tot = new THStack("NT11_EM_TotalBG", "NT11_EM_TotalBG");
	
	nt11_tot->Add(nt11_ss);
	nt11_tot->Add(nt11_wz);
	nt11_tot->Add(nt11_cm);
	nt11_tot->Add(nt11_sf);
	nt11_tot->Add(nt11_ttw);
	nt11_tot->Add(nt11_ttz);
	if (separateTTH) nt11_tot->Add(nt11_tth);

	nt11_mm_tot->Add(nt11_mm_ss);
	nt11_mm_tot->Add(nt11_mm_wz);
	nt11_mm_tot->Add(nt11_mm_sf);
	nt11_mm_tot->Add(nt11_mm_ttw);
	nt11_mm_tot->Add(nt11_mm_ttz);
	if (separateTTH) nt11_mm_tot->Add(nt11_mm_tth);

	nt11_ee_tot->Add(nt11_ee_ss);
	nt11_ee_tot->Add(nt11_ee_wz);
	nt11_ee_tot->Add(nt11_ee_cm);
	nt11_ee_tot->Add(nt11_ee_sf);
	nt11_ee_tot->Add(nt11_ee_ttw);
	nt11_ee_tot->Add(nt11_ee_ttz);
	if (separateTTH) nt11_ee_tot->Add(nt11_ee_tth);

	nt11_em_tot->Add(nt11_em_ss);
	nt11_em_tot->Add(nt11_em_wz);
	nt11_em_tot->Add(nt11_em_cm);
	nt11_em_tot->Add(nt11_em_sf);
	nt11_em_tot->Add(nt11_em_ttw);
	nt11_em_tot->Add(nt11_em_ttz);
	if (separateTTH) nt11_em_tot->Add(nt11_em_tth);

	// Axis labels
	bool intlabel = false;
	if(varbin == 2 || varbin == 6 || varbin == 8 || varbin == 9) intlabel = true;
	TString ytitle = Form("Events / %3.0f GeV", binwidthscale[varbin]);
	if(intlabel) ytitle = "Events";

	nt11_tot->Draw("goff");
	nt11_tot->GetXaxis()->SetTitle(DiffPredYields::axis_label[varbin]);
	nt11_tot->GetXaxis()->SetTitleOffset(1.07);
	nt11_tot->GetYaxis()->SetTitle(ytitle);
	nt11_tot->GetYaxis()->SetTitleSize(0.045);
	nt11_tot->GetXaxis()->SetTitleSize(0.045);
	nt11_tot->GetYaxis()->SetLabelSize(0.045);
	nt11_tot->GetXaxis()->SetLabelSize(0.045);
	nt11_tot->GetYaxis()->SetTitleOffset(1.25);
	nt11_tot->GetXaxis()->SetTitleOffset(1.065);
	if(varbin==2) nt11_tot->GetYaxis()->SetTitleOffset(1.4);

	nt11_mm_tot->Draw("goff");
	nt11_mm_tot->GetXaxis()->SetTitle(DiffPredYields::axis_label[varbin]);
	nt11_mm_tot->GetYaxis()->SetTitle(ytitle);
	nt11_mm_tot->GetXaxis()->SetTitleOffset(1.07);
	nt11_mm_tot->GetYaxis()->SetTitleSize(0.045);
	nt11_mm_tot->GetXaxis()->SetTitleSize(0.045);
	nt11_mm_tot->GetYaxis()->SetLabelSize(0.045);
	nt11_mm_tot->GetXaxis()->SetLabelSize(0.045);
	nt11_mm_tot->GetYaxis()->SetTitleOffset(1.25);
	nt11_mm_tot->GetXaxis()->SetTitleOffset(1.065);
	if(varbin==2) nt11_mm_tot->GetYaxis()->SetTitleOffset(1.4);

	nt11_ee_tot->Draw("goff");
	nt11_ee_tot->GetXaxis()->SetTitle(DiffPredYields::axis_label[varbin]);
	nt11_ee_tot->GetYaxis()->SetTitle(ytitle);
	nt11_ee_tot->GetXaxis()->SetTitleOffset(1.07);
	nt11_ee_tot->GetYaxis()->SetTitleSize(0.045);
	nt11_ee_tot->GetXaxis()->SetTitleSize(0.045);
	nt11_ee_tot->GetYaxis()->SetLabelSize(0.045);
	nt11_ee_tot->GetXaxis()->SetLabelSize(0.045);
	nt11_ee_tot->GetYaxis()->SetTitleOffset(1.25);
	nt11_ee_tot->GetXaxis()->SetTitleOffset(1.065);
	if(varbin==2) nt11_ee_tot->GetYaxis()->SetTitleOffset(1.4);

	nt11_em_tot->Draw("goff");
	nt11_em_tot->GetXaxis()->SetTitle(DiffPredYields::axis_label[varbin]);
	nt11_em_tot->GetYaxis()->SetTitle(ytitle);
	nt11_em_tot->GetXaxis()->SetTitleOffset(1.07);
	nt11_em_tot->GetYaxis()->SetTitleSize(0.045);
	nt11_em_tot->GetXaxis()->SetTitleSize(0.045);
	nt11_em_tot->GetYaxis()->SetLabelSize(0.045);
	nt11_em_tot->GetXaxis()->SetLabelSize(0.045);
	nt11_em_tot->GetYaxis()->SetTitleOffset(1.25);
	nt11_em_tot->GetXaxis()->SetTitleOffset(1.065);
	if(varbin==2) nt11_em_tot->GetYaxis()->SetTitleOffset(1.4);

	if(intlabel){
		for(size_t i = 1; i <= nt11_sf->GetNbinsX(); ++i) nt11_tot   ->GetXaxis()->SetBinLabel(i, Form("%d", i-1));
		for(size_t i = 1; i <= nt11_sf->GetNbinsX(); ++i) nt11_mm_tot->GetXaxis()->SetBinLabel(i, Form("%d", i-1));
		for(size_t i = 1; i <= nt11_sf->GetNbinsX(); ++i) nt11_ee_tot->GetXaxis()->SetBinLabel(i, Form("%d", i-1));
		for(size_t i = 1; i <= nt11_sf->GetNbinsX(); ++i) nt11_em_tot->GetXaxis()->SetBinLabel(i, Form("%d", i-1));
		nt11_tot   ->GetXaxis()->SetLabelSize(0.07);
		nt11_mm_tot->GetXaxis()->SetLabelSize(0.07);
		nt11_ee_tot->GetXaxis()->SetLabelSize(0.07);
		nt11_em_tot->GetXaxis()->SetLabelSize(0.07);
		nt11_tot   ->GetXaxis()->SetTitleSize(0.045);
		nt11_mm_tot->GetXaxis()->SetTitleSize(0.045);
		nt11_ee_tot->GetXaxis()->SetTitleSize(0.045);
		nt11_em_tot->GetXaxis()->SetTitleSize(0.045);
		nt11_tot   ->GetXaxis()->SetTitleOffset(1.07);
		nt11_mm_tot->GetXaxis()->SetTitleOffset(1.07);
		nt11_ee_tot->GetXaxis()->SetTitleOffset(1.07);
		nt11_em_tot->GetXaxis()->SetTitleOffset(1.07);
	}

	nt11_tot   ->SetMinimum(0.5*nt11   ->GetMinimum());
	nt11_mm_tot->SetMinimum(0.5*nt11_mm->GetMinimum());
	nt11_ee_tot->SetMinimum(0.5*nt11_ee->GetMinimum());
	nt11_em_tot->SetMinimum(0.5*nt11_em->GetMinimum());

	double max = nt11->Integral();
	nt11    ->SetMaximum(max>1?max+1:1.);
	nt11_sf ->SetMaximum(max>1?max+1:1.);
	nt11_cm ->SetMaximum(max>1?max+1:1.);
	nt11_ss ->SetMaximum(max>1?max+1:1.);
	nt11_wz ->SetMaximum(max>1?max+1:1.);
	nt11_ttw->SetMaximum(max>1?max+1:1.);
	nt11_ttz->SetMaximum(max>1?max+1:1.);
	if (separateTTH) nt11_tth->SetMaximum(max>1?max+1:1.);
	nt11_tot->SetMaximum(max>1?max+1:1.);

	double max_mm = nt11_mm->Integral();
	nt11_mm    ->SetMaximum(max_mm>1?max_mm+1:1.);
	nt11_mm_sf ->SetMaximum(max_mm>1?max_mm+1:1.);
	nt11_mm_ss ->SetMaximum(max_mm>1?max_mm+1:1.);
	nt11_mm_wz ->SetMaximum(max_mm>1?max_mm+1:1.);
	nt11_mm_ttw->SetMaximum(max_mm>1?max_mm+1:1.);
	nt11_mm_ttz->SetMaximum(max_mm>1?max_mm+1:1.);
	if (separateTTH) nt11_mm_tth->SetMaximum(max_mm>1?max_mm+1:1.);
	nt11_mm_tot->SetMaximum(max_mm>1?max_mm+1:1.);

	double max_ee = nt11_ee->Integral();
	nt11_ee    ->SetMaximum(max_ee>1?max_ee+1:1.);
	nt11_ee_sf ->SetMaximum(max_ee>1?max_ee+1:1.);
	nt11_ee_cm ->SetMaximum(max_ee>1?max_ee+1:1.);
	nt11_ee_ss ->SetMaximum(max_ee>1?max_ee+1:1.);
	nt11_ee_wz ->SetMaximum(max_ee>1?max_ee+1:1.);
	nt11_ee_ttw->SetMaximum(max_ee>1?max_ee+1:1.);
	nt11_ee_ttz->SetMaximum(max_ee>1?max_ee+1:1.);
	if (separateTTH) nt11_ee_tth->SetMaximum(max_ee>1?max_ee+1:1.);
	nt11_ee_tot->SetMaximum(max_ee>1?max_ee+1:1.);

	double max_em = nt11_em->Integral();
	nt11_em    ->SetMaximum(max_em>1?max_em+1:1.);
	nt11_em_sf ->SetMaximum(max_em>1?max_em+1:1.);
	nt11_em_cm ->SetMaximum(max_em>1?max_em+1:1.);
	nt11_em_ss ->SetMaximum(max_em>1?max_em+1:1.);
	nt11_em_wz ->SetMaximum(max_em>1?max_em+1:1.);
	nt11_em_ttw->SetMaximum(max_em>1?max_em+1:1.);
	nt11_em_ttz->SetMaximum(max_em>1?max_em+1:1.);
	if (separateTTH) nt11_em_tth->SetMaximum(max_em>1?max_em+1:1.);
	nt11_em_tot->SetMaximum(max_em>1?max_em+1:1.);
	
	fOutputSubDir = "DiffPredictionPlotsTTW/";
	/////////////////////////////////////////////////////////////////
	TLegend *leg = new TLegend(0.65,0.47,0.90,0.90);
	leg->AddEntry(nt11,     "Data",        "p");
	leg->AddEntry(nt11_ttz, "t#bar{t} + Z","f");
	if (separateTTH) leg->AddEntry(nt11_tth, "t#bar{t} + H","f");
	leg->AddEntry(nt11_ttw, "t#bar{t} + W","f");
	leg->AddEntry(nt11_sf,  "Fakes",       "f");
	leg->AddEntry(nt11_cm,  "Charge MisID","f");
	leg->AddEntry(nt11_wz,  "WZ",          "f");
	leg->AddEntry(nt11_ss,  "Rare SM",     "f");
	// leg->AddEntry(nt11_sig,fSamples[sigsam]->sname,"l");
	leg->SetFillStyle(0);
	leg->SetTextFont(42);
	leg->SetTextSize(0.038);
	leg->SetBorderSize(0);

	FakeRatios *FR = new FakeRatios();
	TGraphAsymmErrors* gr_obs = FR->getGraphPoissonErrors(nt11);
	gr_obs->SetMarkerColor(kBlack);
	gr_obs->SetMarkerStyle(20);
	gr_obs->SetMarkerSize(2.0);
	gr_obs->SetLineWidth(2);
	gr_obs->SetLineColor(kBlack);
	gr_obs->SetFillColor(kBlack);

	TCanvas *c_temp = new TCanvas("C_ObsPred_" + varname, "Observed vs Predicted", 0, 0, 600, 600);
	c_temp->SetLeftMargin(0.12);
	c_temp->SetRightMargin(0.04);
	c_temp->cd();
	gPad->SetLogy();
	
	nt11_tot->Draw("hist");
	totbg   ->DrawCopy("0 E2 same");
	totbg   ->SetFillStyle(0);
	totbg   ->DrawCopy("hist same");
	gr_obs  ->Draw("P same");
	leg->Draw();
	drawDiffCuts(varbin);
	drawTopLine(0.49, 1.0, 0.13);
	
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

	fOutputSubDir = "DiffPredictionPlotsTTW/IndividualChannels/";
	/////////////////////////////////////////////////////////////////
	
	TGraphAsymmErrors* gr_obs_mm = FR->getGraphPoissonErrors(nt11_mm);
	gr_obs_mm->SetMarkerColor(kBlack);
	gr_obs_mm->SetMarkerStyle(20);
	gr_obs_mm->SetMarkerSize(2.0);
	gr_obs_mm->SetLineWidth(2);
	gr_obs_mm->SetLineColor(kBlack);
	gr_obs_mm->SetFillColor(kBlack);
	
	c_temp = new TCanvas("C_ObsPred_MM_" + varname, "Observed vs Predicted", 0, 0, 600, 600);
	c_temp->SetLeftMargin(0.12);
	c_temp->SetRightMargin(0.04);
	c_temp->cd();
	gPad->SetLogy();

	nt11_mm_tot->Draw("hist");
	totbg_mm   ->DrawCopy("0 E2 same");
	totbg_mm   ->SetFillStyle(0);
	totbg_mm   ->DrawCopy("hist same");
	gr_obs_mm  ->Draw("P same");
	leg->Draw();
	lat->SetTextSize(0.04);
	lat->DrawLatex(0.55,0.92, "#mu#mu");
	drawDiffCuts(varbin);
	drawTopLine(0.39, 1.0, 0.13);
	
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
	
	TGraphAsymmErrors* gr_obs_ee = FR->getGraphPoissonErrors(nt11_ee);
	gr_obs_ee->SetMarkerColor(kBlack);
	gr_obs_ee->SetMarkerStyle(20);
	gr_obs_ee->SetMarkerSize(2.0);
	gr_obs_ee->SetLineWidth(2);
	gr_obs_ee->SetLineColor(kBlack);
	gr_obs_ee->SetFillColor(kBlack);
	
	c_temp = new TCanvas("C_ObsPred_EE_" + varname, "Observed vs Predicted", 0, 0, 600, 600);
	c_temp->SetLeftMargin(0.12);
	c_temp->SetRightMargin(0.04);
	c_temp->cd();
	gPad->SetLogy();
	
	nt11_ee_tot->Draw("hist");
	totbg_ee   ->DrawCopy("0 E2 same");
	totbg_ee   ->SetFillStyle(0);
	totbg_ee   ->DrawCopy("hist same");
	gr_obs_ee  ->Draw("P same");
	leg->Draw();
	lat->SetTextSize(0.04);
	lat->DrawLatex(0.55,0.92, "ee");
	drawDiffCuts(varbin);
	drawTopLine(0.49, 1.0, 0.13);
	
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
	
	TGraphAsymmErrors* gr_obs_em = FR->getGraphPoissonErrors(nt11_em);
	gr_obs_em->SetMarkerColor(kBlack);
	gr_obs_em->SetMarkerStyle(20);
	gr_obs_em->SetMarkerSize(2.0);
	gr_obs_em->SetLineWidth(2);
	gr_obs_em->SetLineColor(kBlack);
	gr_obs_em->SetFillColor(kBlack);
	
	c_temp = new TCanvas("C_ObsPred_EM_" + varname, "Observed vs Predicted", 0, 0, 600, 600);
	c_temp->SetLeftMargin(0.12);
	c_temp->SetRightMargin(0.04);
	c_temp->cd();
	gPad->SetLogy();
	
	nt11_em_tot->Draw("hist");
	totbg_em   ->DrawCopy("0 E2 same");
	totbg_em   ->SetFillStyle(0);
	totbg_em   ->DrawCopy("hist same");
	gr_obs_em  ->Draw("P same");
	leg->Draw();
	lat->SetTextSize(0.04);
	lat->DrawLatex(0.55,0.92, "e#mu");
	drawDiffCuts(varbin);
	drawTopLine(0.49, 1.0, 0.13);
	
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
	delete c_temp, leg;
	delete nt11, nt11_mm, nt11_ee, nt11_em;
	// delete nt11_sig, nt11_mm_sig, nt11_ee_sig, nt11_em_sig;
	delete nt10_mm, nt10_em, nt10_ee, nt01_mm, nt01_em, nt01_ee, nt00_mm, nt00_em, nt00_ee;
	delete nt2_os_ee_bb, nt2_os_ee_eb, nt2_os_ee_ee, nt2_os_em_bb, nt2_os_em_ee;
	delete nt11_ss, nt11_mm_ss, nt11_em_ss, nt11_ee_ss;
	delete nt11_wz, nt11_mm_wz, nt11_em_wz, nt11_ee_wz;
	delete nt11_ttw, nt11_mm_ttw, nt11_em_ttw, nt11_ee_ttw;
	delete nt11_ttz, nt11_mm_ttz, nt11_em_ttz, nt11_ee_ttz;
	if (separateTTH) delete nt11_tth, nt11_mm_tth, nt11_em_tth, nt11_ee_tth;
	delete nt11_sf, nt11_mm_sf, nt11_em_sf, nt11_ee_sf;
	delete nt11_cm, nt11_em_cm, nt11_ee_cm;
	// delete nt11_mc, nt11_mm_mc, nt11_ee_mc, nt11_em_mc;
	delete nt11_tot, nt11_mm_tot, nt11_ee_tot, nt11_em_tot;
	delete totbg, totbg_mm, totbg_em, totbg_ee;
	delete FR, gr_obs, gr_obs_mm, gr_obs_ee, gr_obs_em;

	fOutputSubDir = "";
}

//____________________________________________________________________________
void SSDLPlotter::makeAllClosureTests(){
 	TString outputdir = Util::MakeOutputDir(fOutputDir + "MCClosure");
 	for(size_t i = 0; i < gNREGIONS; ++i){
 		TString outputname = outputdir + "MCClosure_" + gRegions[i]->sname + ".txt";
 		makeIntMCClosure(fMCBGNoQCDNoGJets, outputname, i);
 	}
 	for(size_t i = 0; i < gNREGIONS; ++i){
 		TString outputname = outputdir + "MCClosure_Sig_" + gRegions[i]->sname + ".txt";
 		// TString outputname = outputdir + "MCClosure_Sig_" + Region::sname[i] + ".txt";
 		makeIntMCClosure(fMCBGNoQCDNoGJetsSig, outputname, i);
 	}
	//	makeTTbarClosure();
 	fOutputSubDir = "";
}
void SSDLPlotter::makeIntMCClosure(vector<int> samples, TString filename, int reg){
 	ofstream OUT(filename.Data(), ios::trunc);
 	const int nsamples = samples.size();
 
 	OUT << "/////////////////////////////////////////////////////////////////////////////" << endl;
 	OUT << " Producing integrated predictions" << endl;
 	OUT << "  scaling MC to " << fLumiNorm << " /pb" << endl << endl;
 
 	///////////////////////////////////////////////////////////////////////////////////
 	// RATIOS /////////////////////////////////////////////////////////////////////////
 	///////////////////////////////////////////////////////////////////////////////////
 	float mf(0.), mf_e(0.), mp(0.), mp_e(0.), ef(0.), ef_e(0.), ep(0.), ep_e(0.);
 
 	calculateRatio(fMCBGMuEnr, Muon, SigSup, mf, mf_e);
 	calculateRatio(fMCBGMuEnr, Muon, ZDecay, mp, mp_e);
 	calculateRatio(fMCBGEMEnr, Elec, SigSup, ef, ef_e);
 	calculateRatio(fMCBGEMEnr, Elec, ZDecay, ep, ep_e);
 
 	///////////////////////////////////////////////////////////////////////////////////
 	// OBSERVATIONS ///////////////////////////////////////////////////////////////////
 	///////////////////////////////////////////////////////////////////////////////////
 	vector<float> ntt_mm, ntl_mm, nll_mm;
 	vector<float> ntt_ee, ntl_ee, nll_ee;
 	vector<float> ntt_em, ntl_em, nlt_em, nll_em;
 	
 	vector<float> ntt_mm_e1, ntt_ee_e1, ntt_em_e1; // squared stat. errors
 
 	vector<float> npp_mm, npf_mm, nfp_mm, nff_mm;
 	vector<float> npp_ee, npf_ee, nfp_ee, nff_ee;
 	vector<float> npp_em, npf_em, nfp_em, nff_em;
 
 	vector<float> npp_tt_mm, npf_tt_mm, nfp_tt_mm, nff_tt_mm;
 	vector<float> npp_tt_ee, npf_tt_ee, nfp_tt_ee, nff_tt_ee;
 	vector<float> npp_tt_em, npf_tt_em, nfp_tt_em, nff_tt_em;
 
 	// OS yields
 	vector<float> ntt_os_BB_ee, ntt_os_EE_ee, ntt_os_EB_ee;
 	vector<float> ntt_os_BB_em, ntt_os_EE_em;
 	// Squared stat. errors
 	vector<float> ntt_os_BB_ee_e1, ntt_os_EE_ee_e1, ntt_os_EB_ee_e1, ntt_os_BB_em_e1, ntt_os_EE_em_e1;
 
 	// Charge misid
	vector<float> npp_tt_cm_ee, npp_cm_ee;
 	vector<float> npp_tt_cm_em, npp_cm_em;
 
 	vector<float> scales;
 	vector<TString> names;
 	for(size_t i = 0; i < samples.size(); ++i){
 		Sample *S = fSamples[samples[i]];
 		float scale = fLumiNorm / S->getLumi();
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
 
 		ntt_mm_e1.push_back(S->getError(S->region[reg][HighPt].mm.nt20_pt->GetEntries())); // take unweighted entries
 		ntt_ee_e1.push_back(S->getError(S->region[reg][HighPt].ee.nt20_pt->GetEntries()));
 		ntt_em_e1.push_back(S->getError(S->region[reg][HighPt].em.nt20_pt->GetEntries()));
 
 		npp_mm.push_back(S->region[reg][HighPt].mm.npp_pt->GetEntries());
 		npf_mm.push_back(S->region[reg][HighPt].mm.npf_pt->GetEntries());
 		nfp_mm.push_back(S->region[reg][HighPt].mm.nfp_pt->GetEntries());
 		nff_mm.push_back(S->region[reg][HighPt].mm.nff_pt->GetEntries());
 
 		npp_em.push_back(S->region[reg][HighPt].em.npp_pt->GetEntries());
 		npf_em.push_back(S->region[reg][HighPt].em.npf_pt->GetEntries());
 		nfp_em.push_back(S->region[reg][HighPt].em.nfp_pt->GetEntries());
 		nff_em.push_back(S->region[reg][HighPt].em.nff_pt->GetEntries());
 
 		npp_ee.push_back(S->region[reg][HighPt].ee.npp_pt->GetEntries());
 		npf_ee.push_back(S->region[reg][HighPt].ee.npf_pt->GetEntries());
 		nfp_ee.push_back(S->region[reg][HighPt].ee.nfp_pt->GetEntries());
 		nff_ee.push_back(S->region[reg][HighPt].ee.nff_pt->GetEntries());
 
 		npp_tt_mm.push_back(S->region[reg][HighPt].mm.nt2pp_pt->GetEntries());
 		npf_tt_mm.push_back(S->region[reg][HighPt].mm.nt2pf_pt->GetEntries());
 		nfp_tt_mm.push_back(S->region[reg][HighPt].mm.nt2fp_pt->GetEntries());
 		nff_tt_mm.push_back(S->region[reg][HighPt].mm.nt2ff_pt->GetEntries());
 
 		npp_tt_em.push_back(S->region[reg][HighPt].em.nt2pp_pt->GetEntries());
 		npf_tt_em.push_back(S->region[reg][HighPt].em.nt2pf_pt->GetEntries());
 		nfp_tt_em.push_back(S->region[reg][HighPt].em.nt2fp_pt->GetEntries());
 		nff_tt_em.push_back(S->region[reg][HighPt].em.nt2ff_pt->GetEntries());
 
 		npp_tt_ee.push_back(S->region[reg][HighPt].ee.nt2pp_pt->GetEntries());
 		npf_tt_ee.push_back(S->region[reg][HighPt].ee.nt2pf_pt->GetEntries());
 		nfp_tt_ee.push_back(S->region[reg][HighPt].ee.nt2fp_pt->GetEntries());
 		nff_tt_ee.push_back(S->region[reg][HighPt].ee.nt2ff_pt->GetEntries());
 		
 		ntt_os_BB_em.push_back(S->region[reg][HighPt].em.nt20_OS_BB_pt->GetEntries()); // ele in barrel
 		ntt_os_EE_em.push_back(S->region[reg][HighPt].em.nt20_OS_EE_pt->GetEntries()); // ele in endcal
 		ntt_os_BB_ee.push_back(S->region[reg][HighPt].ee.nt20_OS_BB_pt->GetEntries()); // both in barrel
 		ntt_os_EE_ee.push_back(S->region[reg][HighPt].ee.nt20_OS_EE_pt->GetEntries()); // both in endcal
 		ntt_os_EB_ee.push_back(S->region[reg][HighPt].ee.nt20_OS_EB_pt->GetEntries()); // one barrel, one endcap
 		
 		ntt_os_BB_em_e1.push_back(S->getError(S->region[reg][HighPt].em.nt20_OS_BB_pt->GetEntries()));
 		ntt_os_EE_em_e1.push_back(S->getError(S->region[reg][HighPt].em.nt20_OS_EE_pt->GetEntries()));
 		ntt_os_BB_ee_e1.push_back(S->getError(S->region[reg][HighPt].ee.nt20_OS_BB_pt->GetEntries()));
 		ntt_os_EE_ee_e1.push_back(S->getError(S->region[reg][HighPt].ee.nt20_OS_EE_pt->GetEntries()));
 		ntt_os_EB_ee_e1.push_back(S->getError(S->region[reg][HighPt].ee.nt20_OS_EB_pt->GetEntries()));
 		
 		npp_tt_cm_ee.push_back(scale*S->region[reg][HighPt].ee.nt2pp_cm_pt->GetEntries());
 		npp_cm_ee   .push_back(scale*S->region[reg][HighPt].ee.npp_cm_pt->GetEntries());
 		npp_tt_cm_em.push_back(scale*S->region[reg][HighPt].em.nt2pp_cm_pt->GetEntries());
 		npp_cm_em   .push_back(scale*S->region[reg][HighPt].em.npp_cm_pt->GetEntries());
 	}
 
 	///////////////////////////////////////////////////////////////////////////////////
 	// PREDICTIONS ////////////////////////////////////////////////////////////////////
 	///////////////////////////////////////////////////////////////////////////////////
 	FakeRatios *FR = new FakeRatios();
 	FR->setIsMC(true);
 	FR->setNToyMCs(100);
 	FR->setAddESyst(0.0);
 
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
 		FR->setNGen(S->ngen);
 
 		//// Differential Ratios
 		npp_pred_mm.push_back(S->numbers[reg][Muon].npp);
 		npf_pred_mm.push_back(S->numbers[reg][Muon].npf + S->numbers[reg][Muon].nfp);
 		nff_pred_mm.push_back(S->numbers[reg][Muon].nff);
 		npp_pred_ee.push_back(S->numbers[reg][Elec].npp);
 		npf_pred_ee.push_back(S->numbers[reg][Elec].npf + S->numbers[reg][Elec].nfp);
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
 
 		//// Flat Ratios
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
	float fbb(0.),fee(0.),feb(0.);
	float fbbE(0.),feeE(0.),febE(0.);
 	
	calculateChMisIdProb(fMCBGNoQCDNoGJets, BB, fbb, fbbE);
	calculateChMisIdProb(fMCBGNoQCDNoGJets, EB, feb, febE);
	calculateChMisIdProb(fMCBGNoQCDNoGJets, EE, fee, feeE);
	
 	vector<float> ntt_cm_ee, ntt_cm_em;
 	vector<float> ntt_cm_ee_e1, ntt_cm_em_e1; // squared errors, includes stat errors on yields and errors on probabilities
 	for(size_t i = 0; i < nsamples; ++i){
 		Sample *S = fSamples[samples[i]];
 		ntt_cm_ee.push_back(2*fbb*ntt_os_BB_ee[i] + 2*fee*ntt_os_EE_ee[i] + 2*feb*ntt_os_EB_ee[i]);
 		ntt_cm_ee_e1.push_back((4*fbb*fbb* S->getError(ntt_os_BB_ee[i])) + (4*fee*fee*S->getError(ntt_os_EE_ee[i])) + 4*feb*feb*S->getError(ntt_os_EB_ee[i])
 		    + (4*ntt_os_BB_ee[i]*ntt_os_BB_ee[i]*fbbE*fbbE) + (4*ntt_os_EE_ee[i]*ntt_os_EE_ee[i]*feeE*feeE) + 4*febE*febE*ntt_os_EB_ee[i]*ntt_os_EB_ee[i]);
 
 		ntt_cm_em.push_back(  fbb*ntt_os_BB_em[i] +   fee*ntt_os_EE_em[i]);
 		ntt_cm_em_e1.push_back(fbb*fbb*S->getError(ntt_os_BB_em[i]) + fee*fee*S->getError(ntt_os_EE_em[i]) 
 		    + ntt_os_BB_em[i]*ntt_os_BB_em[i] * fbbE*fbbE + ntt_os_EE_em[i]*ntt_os_EE_em[i] * feeE*feeE);
 	}
 
 	// Sums
 	float ntt_sum_mm(0.), ntl_sum_mm(0.), nll_sum_mm(0.);
 	float ntt_sum_em(0.), ntl_sum_em(0.), nlt_sum_em(0.), nll_sum_em(0.);
 	float ntt_sum_ee(0.), ntl_sum_ee(0.), nll_sum_ee(0.);
 
 	float ntt_sum_mm_e1(0.), ntt_sum_em_e1(0.), ntt_sum_ee_e1(0.);
 
 	float npp_sum_mm(0.), npf_sum_mm(0.), nff_sum_mm(0.);
 	float npp_sum_em(0.), npf_sum_em(0.), nfp_sum_em(0.), nff_sum_em(0.);
 	float npp_sum_ee(0.), npf_sum_ee(0.), nff_sum_ee(0.);
 
 	float npp_pred_sum_mm(0.), npf_pred_sum_mm(0.), nff_pred_sum_mm(0.);
 	float npp_pred_sum_em(0.), npf_pred_sum_em(0.), nfp_pred_sum_em(0.), nff_pred_sum_em(0.);
 	float npp_pred_sum_ee(0.), npf_pred_sum_ee(0.), nff_pred_sum_ee(0.);
 
 	float ntt_cm_sum_ee(0.), ntt_cm_sum_em(0.);
 	float ntt_cm_sum_ee_e1(0.), ntt_cm_sum_em_e1(0.);
 
 	float npp_cm_sum_ee(0.),    npp_cm_sum_em(0.);
 	float npp_tt_cm_sum_ee(0.), npp_tt_cm_sum_em(0.);
 
 	///////////////////////////////////////////
 	// Rare SM
 	float ntt_rare_mm(0.), ntt_rare_em(0.), ntt_rare_ee(0.);
 	float ntt_rare_mm_e1(0.), ntt_rare_em_e1(0.), ntt_rare_ee_e1(0.); // squared stat errors
 	for(size_t i = 0; i < fMCRareSM.size(); ++i){
 		Sample *S = fSamples[fMCRareSM[i]];
 		float scale = fLumiNorm/S->getLumi();
 		ntt_rare_mm += scale*S->numbers[reg][Muon].nt2;
 		ntt_rare_em += scale*S->numbers[reg][ElMu].nt2;
 		ntt_rare_ee += scale*S->numbers[reg][Elec].nt2;
 		ntt_rare_mm_e1 += scale*scale*pow(S->getError(S->region[reg][HighPt].mm.nt20_pt->GetEntries()),2);
 		ntt_rare_em_e1 += scale*scale*pow(S->getError(S->region[reg][HighPt].em.nt20_pt->GetEntries()),2);
 		ntt_rare_ee_e1 += scale*scale*pow(S->getError(S->region[reg][HighPt].ee.nt20_pt->GetEntries()),2);
 	}
 
 	///////////////////////////////////////////
 	// WZ production
 	float wzscale = fLumiNorm/fSamples[WZ]->getLumi();
 	float ntt_wz_mm = wzscale*fSamples[WZ]->numbers[reg][Muon].nt2;
 	float ntt_wz_em = wzscale*fSamples[WZ]->numbers[reg][ElMu].nt2;
 	float ntt_wz_ee = wzscale*fSamples[WZ]->numbers[reg][Elec].nt2;
 
 	float ntt_wz_mm_e1 = wzscale*wzscale*pow(fSamples[WZ]->getError(fSamples[WZ]->region[reg][HighPt].mm.nt20_pt->GetEntries()),2); // for stat error take actual entries, not pileup weighted integral...
 	float ntt_wz_em_e1 = wzscale*wzscale*pow(fSamples[WZ]->getError(fSamples[WZ]->region[reg][HighPt].em.nt20_pt->GetEntries()),2);
 	float ntt_wz_ee_e1 = wzscale*wzscale*pow(fSamples[WZ]->getError(fSamples[WZ]->region[reg][HighPt].ee.nt20_pt->GetEntries()),2);
 
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
 		
 		ntt_cm_sum_ee    += scales[i] * ntt_cm_ee[i];
 		ntt_cm_sum_em    += scales[i] * ntt_cm_em[i];
 
 		ntt_cm_sum_ee_e1 += scales[i]*scales[i] * ntt_cm_ee_e1[i];
 		ntt_cm_sum_em_e1 += scales[i]*scales[i] * ntt_cm_em_e1[i];
 
 		npp_cm_sum_ee    += scales[i] * npp_cm_ee[i];
 		npp_cm_sum_em    += scales[i] * npp_cm_em[i];
 		npp_tt_cm_sum_ee += scales[i] * npp_tt_cm_ee[i];
 		npp_tt_cm_sum_em += scales[i] * npp_tt_cm_em[i];
 		
 		ntt_sum_mm_e1    += scales[i]*scales[i] * ntt_mm_e1[i];
 		ntt_sum_em_e1    += scales[i]*scales[i] * ntt_em_e1[i];
 		ntt_sum_ee_e1    += scales[i]*scales[i] * ntt_ee_e1[i];
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
 
 	OUT << "===================================================================================================" << endl;
 	OUT << "                 ||       Mu/Mu       |       E/Mu        |        E/E        |        Sum        |" << endl;
 	OUT << "---------------------------------------------------------------------------------------------------" << endl;
 	float fakesum = npf_pred_sum_mm+npf_pred_sum_em+nfp_pred_sum_em+npf_pred_sum_ee+nff_pred_sum_mm+nff_pred_sum_em+nff_pred_sum_ee;
 	float fakesum_e1 = nF_pred_sum_mm_e1+nF_pred_sum_em_e1+nF_pred_sum_ee_e1;
 	OUT << Form(" Fakes           || %7.2f\\pm%7.2f & %7.2f\\pm%7.2f & %7.2f\\pm%7.2f & %7.2f\\pm%7.2f |",
 	npf_pred_sum_mm+nff_pred_sum_mm,                 sqrt(nF_pred_sum_mm_e1),
 	npf_pred_sum_em+nfp_pred_sum_em+nff_pred_sum_em, sqrt(nF_pred_sum_em_e1),
 	npf_pred_sum_ee+nff_pred_sum_ee,                 sqrt(nF_pred_sum_ee_e1),
 	fakesum, sqrt(fakesum_e1)) << endl;
 	OUT << Form(" Charge MisID    ||                   & %7.2f\\pm%7.2f & %7.2f\\pm%7.2f & %7.2f\\pm%7.2f |",
 	ntt_cm_sum_em,               sqrt(ntt_cm_sum_em_e1),
 	ntt_cm_sum_ee,               sqrt(ntt_cm_sum_ee_e1),
 	ntt_cm_sum_ee+ntt_cm_sum_em, sqrt(ntt_cm_sum_ee_e1+ntt_cm_sum_em_e1)) << endl;
 	OUT << Form(" Irreducible     || %7.2f\\pm%7.2f & %7.2f\\pm%7.2f & %7.2f\\pm%7.2f & %7.2f\\pm%7.2f |",
 	ntt_rare_mm,                         sqrt(ntt_rare_mm_e1),
 	ntt_rare_em,                         sqrt(ntt_rare_em_e1),
 	ntt_rare_ee,                         sqrt(ntt_rare_ee_e1),
 	ntt_rare_mm+ntt_rare_em+ntt_rare_ee, sqrt(ntt_rare_mm_e1+ntt_rare_em_e1+ntt_rare_ee_e1)) << endl;
 	OUT << Form(" WZ Production   || %7.2f\\pm%7.2f & %7.2f\\pm%7.2f & %7.2f\\pm%7.2f & %7.2f\\pm%7.2f |",
 	ntt_wz_mm,                     sqrt(ntt_wz_mm_e1),
 	ntt_wz_em,                     sqrt(ntt_wz_em_e1),
 	ntt_wz_ee,                     sqrt(ntt_wz_ee_e1),
 	ntt_wz_mm+ntt_wz_em+ntt_wz_ee, sqrt(ntt_wz_mm_e1+ntt_wz_em_e1+ntt_wz_ee_e1)) << endl;
 
 	float tot_bg_mm = npf_pred_sum_mm                +nff_pred_sum_mm              +ntt_rare_mm + ntt_wz_mm;
 	float tot_bg_em = npf_pred_sum_em+nfp_pred_sum_em+nff_pred_sum_em+ntt_cm_sum_em+ntt_rare_em + ntt_wz_em;
 	float tot_bg_ee = npf_pred_sum_ee                +nff_pred_sum_ee+ntt_cm_sum_ee+ntt_rare_ee + ntt_wz_ee;
 	float tot_bg = tot_bg_mm + tot_bg_em + tot_bg_ee;
 	float tot_bg_mm_e1 = nF_pred_sum_mm_e1                    + ntt_rare_mm_e1 + ntt_wz_mm_e1;
 	float tot_bg_em_e1 = nF_pred_sum_em_e1 + ntt_cm_sum_em_e1 + ntt_rare_em_e1 + ntt_wz_em_e1;
 	float tot_bg_ee_e1 = nF_pred_sum_ee_e1 + ntt_cm_sum_ee_e1 + ntt_rare_ee_e1 + ntt_wz_ee_e1;
 	float tot_bg_e1    = tot_bg_mm_e1 + tot_bg_em_e1 + tot_bg_ee_e1;
 	float ntt_sum    = ntt_sum_mm+ntt_sum_em+ntt_sum_ee;
 	float ntt_sum_e1 = ntt_sum_mm_e1+ntt_sum_em_e1+ntt_sum_ee_e1;
 	OUT << Form(" Total Pred.     || %7.2f\\pm%7.2f & %7.2f\\pm%7.2f & %7.2f\\pm%7.2f & %7.2f\\pm%7.2f |",
 	tot_bg_mm, sqrt(tot_bg_mm_e1),
 	tot_bg_em, sqrt(tot_bg_em_e1),
 	tot_bg_ee, sqrt(tot_bg_ee_e1),
 	tot_bg,    sqrt(tot_bg_e1) ) << endl;
 	OUT << Form(" Observed        || %7.2f\\pm%7.2f & %7.2f\\pm%7.2f & %7.2f\\pm%7.2f & %7.2f\\pm%7.2f |",
 	ntt_sum_mm, sqrt(ntt_sum_mm_e1),
 	ntt_sum_em, sqrt(ntt_sum_em_e1),
 	ntt_sum_ee, sqrt(ntt_sum_ee_e1),
 	ntt_sum,    sqrt(ntt_sum_e1)) << endl;
 	OUT << Form(" Pred./Obs.      || %7.2f\\pm%7.2f & %7.2f\\pm%7.2f & %7.2f\\pm%7.2f & %7.2f\\pm%7.2f |",
 	tot_bg_mm/ntt_sum_mm,
 	sqrt(nF_pred_sum_mm_e1                  +  ntt_rare_mm_e1 * pow(((ntt_sum_mm-tot_bg_mm)/tot_bg_mm),2) + (tot_bg_mm_e1 - ntt_rare_mm_e1) * pow(tot_bg_mm/ntt_sum_mm,2)  )/ntt_sum_mm,
 	tot_bg_em/ntt_sum_em,
 	sqrt(nF_pred_sum_em_e1+ntt_cm_sum_em_e1 +  ntt_rare_em_e1 * pow(((ntt_sum_em-tot_bg_em)/tot_bg_em),2) + (tot_bg_em_e1 - ntt_rare_em_e1) * pow(tot_bg_em/ntt_sum_em,2)  )/ntt_sum_em,
 	tot_bg_ee/ntt_sum_ee,
 	sqrt(nF_pred_sum_ee_e1+ntt_cm_sum_ee_e1 +  ntt_rare_ee_e1 * pow(((ntt_sum_ee-tot_bg_ee)/tot_bg_ee),2) + (tot_bg_ee_e1 - ntt_rare_ee_e1) * pow(tot_bg_ee/ntt_sum_ee,2)  )/ntt_sum_ee,
 	tot_bg/ntt_sum,
 	sqrt(fakesum_e1 + ntt_cm_sum_ee_e1+ntt_cm_sum_em_e1 +  (ntt_rare_mm_e1+ntt_rare_em_e1+ntt_rare_ee_e1) * pow(((ntt_sum-tot_bg)/tot_bg),2) + (tot_bg_e1 - (ntt_rare_mm_e1+ntt_rare_em_e1+ntt_rare_ee_e1)) * pow(tot_bg/ntt_sum,2)  )/ntt_sum
 	) << endl;
 
 	OUT << Form(" Pred.-Obs./Pred || %7.2f\\pm%7.2f & %7.2f\\pm%7.2f & %7.2f\\pm%7.2f & %7.2f\\pm%7.2f |",
 	(tot_bg_mm-ntt_sum_mm)/tot_bg_mm,
 	sqrt(nF_pred_sum_mm_e1*pow(ntt_sum_mm/tot_bg_mm,2) + ntt_rare_mm_e1*(pow((ntt_sum_mm-tot_bg_mm)/tot_bg_mm,2)) + ntt_sum_mm_e1 - ntt_rare_mm_e1)/tot_bg_mm,
 
 	(tot_bg_em-ntt_sum_em)/tot_bg_em,
 	sqrt((nF_pred_sum_em_e1+ntt_cm_sum_em_e1)*pow(ntt_sum_em/tot_bg_em,2) + ntt_rare_em_e1*(pow((ntt_sum_em-tot_bg_em)/tot_bg_em,2)) + ntt_sum_em_e1 - ntt_rare_em_e1)/tot_bg_em,
 
 	(tot_bg_ee-ntt_sum_ee)/tot_bg_ee,
 	sqrt((nF_pred_sum_ee_e1+ntt_cm_sum_ee_e1)*pow(ntt_sum_ee/tot_bg_ee,2) + ntt_rare_ee_e1*(pow((ntt_sum_ee-tot_bg_ee)/tot_bg_ee,2)) + ntt_sum_ee_e1 - ntt_rare_ee_e1)/tot_bg_ee,
 
 	(tot_bg-ntt_sum)/tot_bg,
 	sqrt((fakesum_e1+ntt_cm_sum_ee_e1+ntt_cm_sum_em_e1)*pow(ntt_sum/tot_bg,2) + (ntt_rare_mm_e1+ntt_rare_em_e1+ntt_rare_ee_e1)*(pow((ntt_sum-tot_bg)/tot_bg,2)) + ntt_sum_e1 - (ntt_rare_mm_e1+ntt_rare_em_e1+ntt_rare_ee_e1))/tot_bg	
 	) << endl;
	OUT << "===================================================================================================" << endl;
 	OUT << endl;
 
 	// OUT << "==========================================================================================================================" << endl;
 	// OUT << "                 ||            Mu/Mu            ||                   E/Mu                ||             E/E             ||" << endl;
 	// OUT << "     PRED. IN TT ||   Npp   |   Nfp   |   Nff   ||   Npp   |   Npf   |   Nfp   |   Nff   ||   Npp   |   Npf   |   Nff   ||" << endl;
 	// OUT << "--------------------------------------------------------------------------------------------------------------------------" << endl;
 	// OUT << setw(16) << "Npf Truth"  << " || ";
 	// OUT << setw(7) << Form("%6.3f", mp*mp*npp_sum_mm) << " | ";
 	// OUT << setw(7) << Form("%6.3f", mp*mf*npf_sum_mm) << " | ";
 	// OUT << setw(7) << Form("%6.3f", mf*mf*nff_sum_mm) << " || ";
 	// OUT << setw(7) << Form("%6.3f", mp*ep*npp_sum_em) << " | ";
 	// OUT << setw(7) << Form("%6.3f", mp*ef*npf_sum_em) << " | ";
 	// OUT << setw(7) << Form("%6.3f", mf*ep*nfp_sum_em) << " | ";
 	// OUT << setw(7) << Form("%6.3f", mf*ef*nff_sum_em) << " || ";
 	// OUT << setw(7) << Form("%6.3f", ep*ep*npp_sum_ee) << " | ";
 	// OUT << setw(7) << Form("%6.3f", ep*ef*npf_sum_ee) << " | ";
 	// OUT << setw(7) << Form("%6.3f", ef*ef*nff_sum_ee) << " || ";
 	// OUT << endl;
 	// OUT << "--------------------------------------------------------------------------------------------------------------------------" << endl;
 	// OUT << setw(16) << "FR Prediction"  << " || ";
 	// OUT << setw(7) << Form("%6.3f", npp_pred_sum_mm) << " | ";
 	// OUT << setw(7) << Form("%6.3f", npf_pred_sum_mm) << " | ";
 	// OUT << setw(7) << Form("%6.3f", nff_pred_sum_mm) << " || ";
 	// OUT << setw(7) << Form("%6.3f", npp_pred_sum_em) << " | ";
 	// OUT << setw(7) << Form("%6.3f", npf_pred_sum_em) << " | ";
 	// OUT << setw(7) << Form("%6.3f", nfp_pred_sum_em) << " | ";
 	// OUT << setw(7) << Form("%6.3f", nff_pred_sum_em) << " || ";
 	// OUT << setw(7) << Form("%6.3f", npp_pred_sum_ee) << " | ";
 	// OUT << setw(7) << Form("%6.3f", npf_pred_sum_ee) << " | ";
 	// OUT << setw(7) << Form("%6.3f", nff_pred_sum_ee) << " || ";
 	// OUT << endl;
 	// OUT << "==========================================================================================================================" << endl;
 	// OUT << setw(16) << "Pred. Fakes"  << " || ";
 	// OUT << setw(7) << Form("%6.3f", npf_pred_sum_mm+nff_pred_sum_mm) << " |                   || ";
 	// OUT << setw(7) << Form("%6.3f", npf_pred_sum_em+nfp_pred_sum_em+nff_pred_sum_em) << " |                             || ";
 	// OUT << setw(7) << Form("%6.3f", npf_pred_sum_ee+nff_pred_sum_ee) << " |                   || ";
 	// OUT << endl;
 	// OUT << setw(16) << "Pred. Ch-MID"  << " ||         |                   || ";
 	// OUT << setw(7) << Form("%6.3f", ntt_cm_sum_em) << " |                             || ";
 	// OUT << setw(7) << Form("%6.3f", ntt_cm_sum_ee) << " |                   || ";
 	// OUT << endl;
 	// OUT << setw(16) << "Pred. Rare SM"  << " || ";
 	// OUT << setw(7) << Form("%6.3f", ntt_rare_mm) << " |                   || ";
 	// OUT << setw(7) << Form("%6.3f", ntt_rare_em) << " |                             || ";
 	// OUT << setw(7) << Form("%6.3f", ntt_rare_ee) << " |                   || ";
 	// OUT << endl;
 	// OUT << "--------------------------------------------------------------------------------------------------------------------------" << endl;
 	// OUT << setw(16) << "Total BG Pred."  << " || ";
 	// OUT << setw(7) << Form("%6.3f", tot_bg_mm) << " | P/O  | P-O/P |    || ";
 	// OUT << setw(7) << Form("%6.3f", tot_bg_em) << " | P/O  | P-O/P |              || ";
 	// OUT << setw(7) << Form("%6.3f", tot_bg_ee) << " | P/O  | P-O/P |    || ";
 	// OUT << endl;
 	// OUT << "--------------------------------------------------------------------------------------------------------------------------" << endl;
 	// OUT << setw(16) << "Observed"  << " || ";
 	// OUT << setw(7) << Form(" %6.3f | %4.2f | %5.2f |", ntt_sum_mm, tot_bg_mm/ntt_sum_mm, (tot_bg_mm-ntt_sum_mm)/tot_bg_mm) << "    || ";
 	// OUT << setw(7) << Form(" %6.3f | %4.2f | %5.2f |", ntt_sum_em, tot_bg_em/ntt_sum_em, (tot_bg_em-ntt_sum_em)/tot_bg_em) << "              || ";
 	// OUT << setw(7) << Form(" %6.3f | %4.2f | %5.2f |", ntt_sum_ee, tot_bg_ee/ntt_sum_ee, (tot_bg_ee-ntt_sum_ee)/tot_bg_ee) << "    || ";
 	// OUT << endl;
 	// OUT << "==========================================================================================================================" << endl;
 	// OUT << Form("     Comb. Fakes || %6.3f |              ||",
 	// npf_pred_sum_mm+nff_pred_sum_mm + 
 	// npf_pred_sum_em+nfp_pred_sum_em+nff_pred_sum_em +
 	// npf_pred_sum_ee+nff_pred_sum_ee) << endl;
 	// OUT << Form("    Comb. Ch-MID || %6.3f |              ||", ntt_cm_sum_em+ntt_cm_sum_ee) << endl;
 	// OUT << Form("   Comb. Rare SM || %6.3f |              ||", ntt_rare_mm+ntt_rare_em+ntt_rare_ee) << endl;
 	// OUT << "---------------------------------------------" << endl;
 	// OUT << Form("  Total BG Comb. || %6.3f | P/O  | P-O/P ||", tot_bg) << endl;
 	// OUT << "---------------------------------------------" << endl;
 	// OUT << Form("        Observed || %6.3f | %4.2f | %5.2f ||", ntt_sum, tot_bg/ntt_sum, (tot_bg-ntt_sum)/tot_bg) << endl;
 	// OUT << "=============================================" << endl;
 	// OUT << endl;
 
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
 	// 
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
 	// OUT << "  Total fakes:   " <<  Form("%6.3f", npf_pred_sum_ee+nff_pred_sum_ee) << " ±   " << Form("%6.3f", sqrt(nF_pred_sum_ee_e1)) << " (stat)" << endl;
 	// OUT << "==========================================================================================================================" << endl;
 	// OUT << "/////////////////////////////////////////////////////////////////////////////" << endl;
 	
 	OUT.close();
 	delete FR;
}
void SSDLPlotter::makeWZValidation(SSPrediction pred){
        fOutputSubDir = "WZValidation/";
	TString outputdir = Util::MakeOutputDir(fOutputDir + fOutputSubDir);
	TString outputname = outputdir + "WZValidation.txt";
	
	
	////////////////////////////////////////
	//// MAKE OUTPUT AS PLOTS
	///////////////////////////////////////
	TH1D    *hvar_data[gNKinVars];
	TH1D    *hvar_qcd [gNKinVars];
	TH1D    *hvar_ttj [gNKinVars];
	TH1D    *hvar_ewk [gNKinVars];
	TH1D    *hvar_rare[gNKinVars];
	TH1D    *hvar_wz  [gNKinVars];
	TH1D    *hvar_tot [gNKinVars];
	
	TH1D    *hvar_rat [gNKinVars]; // ratio
	THStack *hvar_mc_s[gNKinVars];

	TLatex *lat = new TLatex();
	lat->SetNDC(kTRUE);
	lat->SetTextColor(kBlack);
	lat->SetTextSize(0.04);

	// Create histograms
	for(size_t i = 0; i < gNKinVars; ++i){
		hvar_data[i] = new TH1D("Data_"    + KinPlots::var_name[i], KinPlots::var_name[i] + " in Data", KinPlots::nbins[i], KinPlots::xmin[i], KinPlots::xmax[i]);
		hvar_qcd [i] = new TH1D("QCD_"     + KinPlots::var_name[i], KinPlots::var_name[i] + " in MC",   KinPlots::nbins[i], KinPlots::xmin[i], KinPlots::xmax[i]);
		hvar_ttj [i] = new TH1D("TTjets_"  + KinPlots::var_name[i], KinPlots::var_name[i] + " in MC",   KinPlots::nbins[i], KinPlots::xmin[i], KinPlots::xmax[i]);
		hvar_ewk [i] = new TH1D("EWK_"     + KinPlots::var_name[i], KinPlots::var_name[i] + " in MC",   KinPlots::nbins[i], KinPlots::xmin[i], KinPlots::xmax[i]);
		hvar_rare[i] = new TH1D("Rare_"    + KinPlots::var_name[i], KinPlots::var_name[i] + " in MC",   KinPlots::nbins[i], KinPlots::xmin[i], KinPlots::xmax[i]);
		hvar_wz[i]   = new TH1D("WZ_"      + KinPlots::var_name[i], KinPlots::var_name[i] + " in MC",   KinPlots::nbins[i], KinPlots::xmin[i], KinPlots::xmax[i]);
		hvar_tot[i]  = new TH1D("Tot_"     + KinPlots::var_name[i], KinPlots::var_name[i] + " in MC",   KinPlots::nbins[i], KinPlots::xmin[i], KinPlots::xmax[i]);
		hvar_rat[i]  = new TH1D("Ratio_"   + KinPlots::var_name[i], KinPlots::var_name[i] + " Data/MC", KinPlots::nbins[i], KinPlots::xmin[i], KinPlots::xmax[i]);

		hvar_mc_s[i] = new THStack("MC_stacked_" + KinPlots::var_name[i], KinPlots::var_name[i] + " in MC");
	}
	
	vector<int> mcsamples   = fMCBG; //NoQCDNoGJets;
	vector<int> datasamples = fHighPtData;
		
	//////////////////////////////////////////////////////////
	// Make kin plots
	for(size_t i = 0; i < gNKinVars; ++i){
	  // Create plots
	  bool intlabel = false;
	  if(i == 2 || i == 10 || i == 11) intlabel = true;
	    
	  hvar_data[i]->SetXTitle(KinPlots::axis_label[i]);
	  hvar_data[i]->SetLineWidth(3);
	  hvar_data[i]->SetLineColor(kBlack);
	  hvar_data[i]->SetMarkerStyle(8);
	  hvar_data[i]->SetMarkerColor(kBlack);
	  hvar_data[i]->SetMarkerSize(1.3);
	  
	  // Scale by luminosity
	  for(size_t j = 0; j < gNSAMPLES; ++j){
	    float lumiscale = fLumiNorm / fSamples[j]->getLumi();
	    if(fSamples[j]->datamc == 0) continue;
	    fSamples[j]->kinplots_wz[2].hvar[i]->Scale(lumiscale);
	  }
	  
	  // Fill data histo
	  for(size_t j = 0; j < datasamples.size(); ++j){
	    Sample *S = fSamples[datasamples[j]];
	    hvar_data[i]->Add(S->kinplots_wz[2].hvar[i]);
	    hvar_data[i]->SetXTitle(KinPlots::axis_label[i]);
	  }
	  
	  
	  hvar_qcd [i]->SetFillColor(kYellow-4);
	  hvar_wz  [i]->SetFillColor(kSpring-9);
	  hvar_ewk [i]->SetFillColor(kAzure+8);
	  hvar_ttj [i]->SetFillColor(kAzure-5);
	  hvar_rare[i]->SetFillColor(kGreen+1);


	  // Fill MC stacks
	  for(size_t j = 0; j < mcsamples.size();   ++j){
	    Sample *S = fSamples[mcsamples[j]];
	    TString s_name = S->sname;
	    //			  if (s_name == "TTJets2") continue;
	    
	    if ( S->getType() == 1)                        hvar_qcd [i]->Add( S->kinplots_wz[2].hvar[i] ); // QCD
	    if ( S->getType() == 5 && S->getProc() == 7 )  hvar_wz  [i]->Add( S->kinplots_wz[2].hvar[i] ); // WZ 
	    if ( S->getType() == 2)                        hvar_ttj [i]->Add( S->kinplots_wz[2].hvar[i] ); // TOP      
	    if ( S->getType() == 3 || 
		 (S->getType() == 5 && S->getProc() != 7)) hvar_ewk [i]->Add( S->kinplots_wz[2].hvar[i] ); // EWK+DB      
	    if ( S->getType() == 4)                        hvar_rare[i]->Add( S->kinplots_wz[2].hvar[i] ); // RARE    
	  }
	  
	  hvar_mc_s[i]->Add(hvar_ttj[i]);
	  hvar_mc_s[i]->Add(hvar_ewk[i]);
	  hvar_mc_s[i]->Add(hvar_rare[i]);
	  hvar_mc_s[i]->Add(hvar_wz[i]);
	  hvar_mc_s[i]->Add(hvar_qcd[i]);
	  hvar_mc_s[i]->Draw("goff");
	  hvar_mc_s[i]->GetXaxis()->SetTitle(KinPlots::axis_label[i]);
	  for(size_t j = 1; j <= hvar_data[i]->GetNbinsX(); ++j)            hvar_data[i]->GetXaxis()->SetBinLabel(j, "");
	  for(size_t j = 1; j <= hvar_mc_s[i]->GetXaxis()->GetNbins(); ++j) hvar_mc_s[i]->GetXaxis()->SetBinLabel(j, "");
	  
	  hvar_tot[i]->Add(hvar_qcd[i]);
	  hvar_tot[i]->Add(hvar_wz[i]);
	  hvar_tot[i]->Add(hvar_ewk[i]);
	  hvar_tot[i]->Add(hvar_rare[i]);
	  hvar_tot[i]->Add(hvar_ttj[i]);

	  float binwidth = hvar_data[i]->GetBinWidth(1);
	  TString ytitle = Form("Events / %3.0f GeV", binwidth);
	  if(intlabel) ytitle = "Events"; // Njets, Nbjets
	  hvar_mc_s[i]->GetYaxis()->SetTitle(ytitle);
	  if(i==2) hvar_mc_s[i]->GetYaxis()->SetTitleOffset(1.25);
	  hvar_mc_s[i]->GetYaxis()->SetTitleSize(0.04);
	  
	  double max1 = hvar_mc_s[i]->GetMaximum();
	  double max2 = hvar_data[i]->GetMaximum();
	  double max = max1>max2?max1:max2;
	  hvar_mc_s[i]->SetMaximum(1.5*max);
	  hvar_data[i]->SetMaximum(1.5*max);
	  hvar_mc_s[i]->SetMinimum(0.);
	  hvar_data[i]->SetMinimum(0.);				
	  			
	  // Ratio plot:
	  float border = 0.3;
	  float scale = (1-border)/border;
						
	  hvar_rat[i]->SetXTitle(hvar_data[i]->GetXaxis()->GetTitle());
	  hvar_rat[i]->SetYTitle("");
	  hvar_rat[i]->GetXaxis()->SetTitleSize(scale * 0.04);
	  hvar_rat[i]->GetXaxis()->SetLabelSize(scale * hvar_data[i]->GetXaxis()->GetLabelSize());
	  hvar_rat[i]->GetYaxis()->SetLabelSize(scale * hvar_data[i]->GetYaxis()->GetLabelSize());
	  hvar_rat[i]->GetXaxis()->SetTickLength(scale * hvar_data[i]->GetXaxis()->GetTickLength());
	  hvar_rat[i]->GetYaxis()->SetTickLength(hvar_data[i]->GetYaxis()->GetTickLength());
	  if(intlabel){
	    hvar_rat[i]->GetXaxis()->SetLabelSize(scale*0.06);
	    hvar_rat[i]->GetXaxis()->SetLabelOffset(0.02);
	    hvar_rat[i]->GetXaxis()->SetTitleOffset(1.20);
	    for(size_t j = 1; j <= hvar_rat[i]->GetXaxis()->GetNbins(); ++j) hvar_rat[i]->GetXaxis()->SetBinLabel(j, Form("%d", j-1));
	  }
	  
	  hvar_rat[i]->SetFillStyle(1001);
	  hvar_rat[i]->SetLineWidth(1);
	  hvar_rat[i]->SetFillColor(  kGray+1);
	  hvar_rat[i]->SetLineColor(  kGray+1);
	  hvar_rat[i]->SetMarkerColor(kGray+1);
	  
	  hvar_rat[i]->Divide(hvar_data[i], hvar_tot[i]);
			
	  // Canvas
	  TCanvas *c_temp = new TCanvas("C_" + KinPlots::var_name[i], KinPlots::var_name[i] + " in Data vs MC", 0, 0, 600, 600);
	  c_temp->cd();
	  
	  TPad *p_plot  = new TPad("plotpad",  "Pad containing the plot", 0.00, border, 1.00, 1.00, 0, 0);
	  p_plot->SetBottomMargin(0.015);
	  p_plot->Draw();
	  TPad *p_ratio = new TPad("ratiopad", "Pad containing the ratio", 0.00, 0.00, 1.00, border, 0, 0);
	  p_ratio->SetTopMargin(0.025);
	  p_ratio->SetBottomMargin(0.35);
	  p_ratio->Draw();

	  p_ratio->cd();
	  hvar_rat[i]->GetYaxis()->SetNdivisions(505);
	  // setPlottingRange(hvar_rat[i], 0.3);
	  hvar_rat[i]->SetMaximum(1.99);
	  hvar_rat[i]->SetMinimum(0.0);
	  hvar_rat[i]->DrawCopy("E2 ");
	  TLine *l3 = new TLine(hvar_data[i]->GetXaxis()->GetXmin(), 1.00, hvar_data[i]->GetXaxis()->GetXmax(), 1.00);
	  l3->SetLineWidth(2);
	  l3->SetLineStyle(7);
	  l3->Draw();
	  gPad->RedrawAxis();
	  p_ratio->Draw();
	  
	  FakeRatios *FR = new FakeRatios();
	  TGraphAsymmErrors* gr_obs = FR->getGraphPoissonErrors(hvar_data[i]);
	  gr_obs->SetMarkerColor(kBlack);
	  gr_obs->SetMarkerStyle(8);
	  gr_obs->SetMarkerSize(1.2);
	  gr_obs->SetLineWidth(2);
	  gr_obs->SetLineColor(kBlack);
	  gr_obs->SetFillColor(kBlack);
	  
	  
	  TLegend *leg = new TLegend(0.70,0.62,0.89,0.88);
	  leg->AddEntry(hvar_data[i], "Data",      "p");
	  leg->AddEntry(hvar_qcd[i],  "QCD",       "f");
	  leg->AddEntry(hvar_wz[i],   "WZ",       "f");
	  leg->AddEntry(hvar_rare[i], "Rare SM",   "f");
	  leg->AddEntry(hvar_ewk[i],  "W/Z/WW/ZZ", "f");
	  leg->AddEntry(hvar_ttj[i],  "Top",       "f");
	  leg->AddEntry(hvar_rat[i],  "Ratio",     "f");
	  
	  leg->SetFillStyle(0);
	  leg->SetTextFont(42);
	  leg->SetBorderSize(0);

	  // LINEAR
	  p_plot->cd();
	  p_plot->SetLogy(0);
	  hvar_mc_s[i]->Draw("hist");
	  gr_obs->Draw("P same");
	  leg->Draw();
	   
	  drawTopLine();

	  if(i < 5)  lat->DrawLatex(0.14,0.85, "ee/e#mu/#mu#mu");
	  if(i == 5) lat->DrawLatex(0.14,0.85, "ee/#mu#mu");
	  if(i == 6) lat->DrawLatex(0.14,0.85, "#mu#mu");
	  if(i == 7) lat->DrawLatex(0.14,0.85, "ee");
	  if(i == 8) lat->DrawLatex(0.14,0.85, "e#mu");
	  if(i > 8)  lat->DrawLatex(0.14,0.85, "ee/e#mu/#mu#mu");
	  
	  p_plot->Draw();
	  gPad->RedrawAxis();
	  c_temp->Update();
	  
	  Util::PrintPDF(c_temp, KinPlots::var_name[i], fOutputDir + fOutputSubDir);
	  Util::PrintPNG(c_temp, KinPlots::var_name[i], fOutputDir + fOutputSubDir);
	  
	  // LOG
	  hvar_mc_s[i]->SetMaximum(5.*max);
	  hvar_data[i]->SetMaximum(5.*max);
	  hvar_mc_s[i]->SetMinimum(0.5);
	  hvar_data[i]->SetMinimum(0.5);
	  	 
	  p_plot->cd();
	  p_plot->SetLogy(1);
	  hvar_mc_s[i]->Draw("hist");
	  gr_obs->Draw("P same");
	  leg->Draw();
	  
	  drawTopLine();
	  if(i < 5)  lat->DrawLatex(0.14,0.85, "ee/e#mu/#mu#mu");
	  if(i == 5) lat->DrawLatex(0.14,0.85, "ee/#mu#mu");
	  if(i == 6) lat->DrawLatex(0.14,0.85, "#mu#mu");
	  if(i == 7) lat->DrawLatex(0.14,0.85, "ee");
	  if(i == 8) lat->DrawLatex(0.14,0.85, "e#mu");
	  if(i > 8)  lat->DrawLatex(0.14,0.85, "ee/e#mu/#mu#mu");
	  
	  p_plot->Draw();
	  gPad->RedrawAxis();
	  c_temp->Update();

	  Util::PrintPDF(c_temp, KinPlots::var_name[i]+"_log", fOutputDir + fOutputSubDir);  
	  Util::PrintPNG(c_temp, KinPlots::var_name[i]+"_log", fOutputDir + fOutputSubDir);  
	  
	  delete gr_obs;
	  delete c_temp;
	  delete leg;
	}
	

	ofstream OUT(outputname.Data(), ios::trunc);
	OUT << "/////////////////////////////////////////////////////////////////////////////" << endl;
	OUT << " Producing WZ Validation in a WZ-enriched sample                             " << endl;
	OUT << "  scaling MC to " << fLumiNorm << " /pb                                      " << endl;
	OUT << "/////////////////////////////////////////////////////////////////////////////" << endl;
	OUT << " " << endl;
	
 	float finalWZ  = hvar_wz[1] ->GetEntries();
	float fakes    = hvar_ttj[1]->GetSumOfWeights();
	float fakes_mm = hvar_ttj[6]->GetSumOfWeights();
	float fakes_ee = hvar_ttj[7]->GetSumOfWeights();
	float fakes_em = hvar_ttj[8]->GetSumOfWeights();
	
	float obsbkg_mm(0) , obsbkg_em(0) , obsbkg_ee(0) , obsbkg(0) ;
	float obspred_mm(0), obspred_em(0), obspred_ee(0), obspred(0);
	float obsbkg_err_mm(0) , obsbkg_err_em(0) , obsbkg_err_ee(0) , obsbkg_err(0) ;
	
	obsbkg_mm  = (pred.obs_mm - (pred.cmid_mm+pred.rare_mm)) / pred.wz_mm;
	obsbkg_em  = (pred.obs_em - (pred.cmid_em+pred.rare_em)) / pred.wz_em;
	obsbkg_ee  = (pred.obs_ee - (pred.cmid_ee+pred.rare_ee)) / pred.wz_ee;
	obsbkg     = (pred.obs    - (pred.cmid   +pred.rare   )) / pred.wz   ;
	obspred_mm = (pred.obs_mm - (pred.cmid_mm+pred.rare_mm+pred.wz_mm)) / pred.wz_mm;
	obspred_em = (pred.obs_em - (pred.cmid_em+pred.rare_em+pred.wz_em)) / pred.wz_em;
	obspred_ee = (pred.obs_ee - (pred.cmid_ee+pred.rare_ee+pred.wz_ee)) / pred.wz_ee;
	obspred    = (pred.obs    - (pred.cmid   +pred.rare   +pred.wz   )) / pred.wz   ;

	obsbkg_err_mm  = sqrt((pred.cmid_stat_mm/pred.wz_mm)*(pred.cmid_stat_mm/pred.wz_mm)+(pred.rare_stat_mm/pred.wz_mm)*(pred.rare_stat_mm/pred.wz_mm)+(obsbkg_mm*pred.wz_stat_mm/pred.wz_mm)*(obsbkg_mm*pred.wz_stat_mm/pred.wz_mm));
	obsbkg_err_em  = sqrt((pred.cmid_stat_em/pred.wz_em)*(pred.cmid_stat_em/pred.wz_em)+(pred.rare_stat_em/pred.wz_em)*(pred.rare_stat_em/pred.wz_em)+(obsbkg_em*pred.wz_stat_em/pred.wz_em)*(obsbkg_em*pred.wz_stat_em/pred.wz_em));
	obsbkg_err_ee  = sqrt((pred.cmid_stat_ee/pred.wz_ee)*(pred.cmid_stat_ee/pred.wz_ee)+(pred.rare_stat_ee/pred.wz_ee)*(pred.rare_stat_ee/pred.wz_ee)+(obsbkg_ee*pred.wz_stat_ee/pred.wz_ee)*(obsbkg_ee*pred.wz_stat_ee/pred.wz_ee));
	obsbkg_err     = sqrt((pred.cmid_stat/pred.wz)      *(pred.cmid_stat/pred.wz)      +(pred.rare_stat/pred.wz)      *(pred.rare_stat/pred.wz)      +(obsbkg*pred.wz_stat/pred.wz)         *(obsbkg*pred.wz_stat/pred.wz));

	OUT << endl;
	OUT << "-----------------------------------------------------------------------------------------------------------------------" << endl;
 	OUT << "       SUMMARY   ||         Mu/Mu         ||         E/Mu          ||          E/E          ||         TOTAL         ||" << endl;
 	OUT << "=======================================================================================================================" << endl;
	OUT << Form("%16s ||                       || %5.2f ± %5.2f ± %5.2f || %5.2f ± %5.2f ± %5.2f || %5.2f ± %5.2f ± %5.2f ||\n", "Pred. chmisid",
		    pred.cmid_em, pred.cmid_stat_em, pred.cmid_err_em, 
		    pred.cmid_ee, pred.cmid_stat_ee, pred.cmid_err_ee,
		    pred.cmid   , pred.cmid_stat   , pred.cmid_err   );
 	OUT << Form("%16s || %5.2f ± %5.2f ± %5.2f || %5.2f ± %5.2f ± %5.2f || %5.2f ± %5.2f ± %5.2f || %5.2f ± %5.2f ± %5.2f ||\n", "Rare SM (Sum)",
		    pred.rare_mm, pred.rare_stat_mm, pred.rare_err_mm, 
		    pred.rare_em, pred.rare_stat_em, pred.rare_err_em, 
		    pred.rare_ee, pred.rare_stat_ee, pred.rare_err_ee,
		    pred.rare   , pred.rare_stat   , pred.rare_err   );
	OUT << Form("%16s || %5.2f ± %5.2f ± %5.2f || %5.2f ± %5.2f ± %5.2f || %5.2f ± %5.2f ± %5.2f || %5.2f ± %5.2f ± %5.2f ||\n", "WZ Prod",
		    pred.wz_mm, pred.wz_stat_mm, pred.wz_err_mm, 
		    pred.wz_em, pred.wz_stat_em, pred.wz_err_em, 
		    pred.wz_ee, pred.wz_stat_ee, pred.wz_err_ee,
		    pred.wz   , pred.wz_stat   , pred.wz_err   );
	OUT << "-----------------------------------------------------------------------------------------------------------------------" << endl;
	OUT << Form("%16s || %d                    || %d                    || %d                    || %d                    ||\n", "Observed",
		    pred.obs_mm, pred.obs_em, pred.obs_ee,  pred.obs    );
	OUT << Form("%16s || %5.2f ± %5.2f ± %5.2f || %5.2f ± %5.2f ± %5.2f || %5.2f ± %5.2f ± %5.2f || %5.2f ± %5.2f ± %5.2f ||\n", "Background",
		    pred.cmid_mm+pred.rare_mm+pred.wz_mm+fakes_mm, sqrt(pred.cmid_stat_mm*pred.cmid_stat_mm+pred.rare_stat_mm*pred.rare_stat_mm+pred.wz_stat_mm*pred.wz_stat_mm), sqrt(pred.cmid_err_mm*pred.cmid_err_mm+pred.rare_err_mm*pred.rare_err_mm+pred.wz_err_mm*pred.wz_err_mm), 
		    pred.cmid_em+pred.rare_em+pred.wz_em+fakes_em, sqrt(pred.cmid_stat_em*pred.cmid_stat_em+pred.rare_stat_em*pred.rare_stat_em+pred.wz_stat_em*pred.wz_stat_em), sqrt(pred.cmid_err_em*pred.cmid_err_em+pred.rare_err_em*pred.rare_err_em+pred.wz_err_em*pred.wz_err_em), 
		    pred.cmid_ee+pred.rare_ee+pred.wz_ee+fakes_ee, sqrt(pred.cmid_stat_ee*pred.cmid_stat_ee+pred.rare_stat_ee*pred.rare_stat_ee+pred.wz_stat_ee*pred.wz_stat_ee), sqrt(pred.cmid_err_ee*pred.cmid_err_ee+pred.rare_err_ee*pred.rare_err_ee+pred.wz_err_ee*pred.wz_err_ee), 
		    pred.cmid   +pred.rare   +pred.wz   +fakes   , sqrt(pred.cmid_stat   *pred.cmid_stat   +pred.rare_stat   *pred.rare_stat   +pred.wz_stat   *pred.wz_stat   ), sqrt(pred.cmid_err   *pred.cmid_err   +pred.rare_err   *pred.rare_err   +pred.wz_err   *pred.wz_err   )); 
	OUT << "-----------------------------------------------------------------------------------------------------------------------" << endl;
	OUT << Form("%16s || %5.2f ± %5.2f         || %5.2f ± %5.2f         || %5.2f ± %5.2f         || %5.2f ± %5.2f         ||\n", "Obs. - Bkg./ WZ",
		    obsbkg_mm, obsbkg_err_mm, obsbkg_em, obsbkg_err_em, obsbkg_ee, obsbkg_err_ee, obsbkg,    obsbkg_err     );
	OUT << Form("%16s || %5.2f ± %5.2f         || %5.2f ± %5.2f         || %5.2f ± %5.2f         || %5.2f ± %5.2f         ||\n", "Obs. - Pred./ WZ",
		    obspred_mm, obsbkg_err_mm, obspred_em, obsbkg_err_em, obspred_ee, obsbkg_err_ee, obspred,    obsbkg_err    );
	OUT << "------------------------------------------------------------------------------------------------------------------------" << endl;
	OUT << endl;
	OUT << "Including fakes taken from ttbar MC sample... " << endl;

	float syst1(0), syst2(0);
	syst1 = sqrt( obspred_mm*obspred_mm*((float)pred.obs_mm/(float)pred.obs) + obspred_em*obspred_em*((float)pred.obs_em/(float)pred.obs) + obspred_ee*obspred_ee*((float)pred.obs_ee/(float)pred.obs));
	syst2 = sqrt(obspred_mm*obspred_mm)*((float)pred.obs_mm/(float)pred.obs) + sqrt(obspred_em*obspred_em)*((float)pred.obs_em/(float)pred.obs) + sqrt(obspred_ee*obspred_ee)*((float)pred.obs_ee/(float)pred.obs);

	OUT << " Taking as systematic the mean syst. error among all channels... " << endl;
	OUT << Form("Syst = %4.3f", syst2) << endl;
	OUT << endl;
	OUT << endl;

	/// CALCULATE XSECTION
	float eff    = finalWZ / 2011253;
	float W2e    = 0.1075;
	float W2m    = 0.1057;
	float W2tau  = 0.1125;
	float Z2ll   = 0.033658;
	float WZ3lnu = 3*Z2ll*(W2e + W2m + W2tau);	//	float WZ3lnu = 0.04036309056;
	float xsec   = (pred.obs - (pred.rare + pred.cmid + fakes))/(eff * fLumiNorm * WZ3lnu);
	OUT << endl;
	OUT << endl;
	OUT << "The total WZ cross-section is: " << xsec << " ± " << syst2*xsec <<" (syst) ± " << obsbkg_err*xsec << " (stat) ± " << 0.044*xsec << " (lumi)" << endl;
	xsec = (pred.obs - (pred.rare + pred.cmid + fakes))/(eff * fLumiNorm);
	OUT << "The WZ->3l3N cross-section is: " << xsec << " ± " << syst2*xsec <<" (syst) ± " << obsbkg_err*xsec << " (stat) ± " << 0.044*xsec << " (lumi)" << endl;
	OUT << "-------------------------------------------------------------" << endl;

	fOutputSubDir = "IntPredictions/";
}
void SSDLPlotter::makeAllClosureTestsTTW(){
 	TString outputdir = Util::MakeOutputDir(fOutputDir + "MCClosureTTWZ");
 	for(size_t i = 0; i < gNREGIONS; ++i){
 		TString outputname = outputdir + "MCClosure_" + gRegions[i]->sname + ".txt";
 		makeIntMCClosureTTW(fClosureSamples, outputname, i);
// 		makeIntMCClosureTTW(fTTJets, outputname, i);
 	}
// LUKAS 	for(size_t i = 0; i < gNREGIONS; ++i){
// LUKAS 		TString outputname = outputdir + "MCClosure_Sig_" + gRegions[i]->sname + ".txt";
// LUKAS 		// TString outputname = outputdir + "MCClosure_Sig_" + Region::sname[i] + ".txt";
// LUKAS 		makeIntMCClosureTTW(fMCBGNoQCDNoGJetsSig, outputname, i);
// LUKAS 	}
 	fOutputSubDir = "";
}
void SSDLPlotter::makeIntMCClosureTTW(vector<int> samples, TString filename, int reg){
 	ofstream OUT(filename.Data(), ios::trunc);
 	const int nsamples = samples.size();
	
 	OUT << "/////////////////////////////////////////////////////////////////////////////" << endl;
 	OUT << " Producing integrated predictions" << endl;
 	OUT << "  scaling MC to " << fLumiNorm << " /pb" << endl << endl;
	
 	///////////////////////////////////////////////////////////////////////////////////
 	// RATIOS /////////////////////////////////////////////////////////////////////////
 	///////////////////////////////////////////////////////////////////////////////////
 	float mf(0.), mf_e(0.), mp(0.), mp_e(0.), ef(0.), ef_e(0.), ep(0.), ep_e(0.);
	
// 	calculateRatio(fTTJets, Muon, SigSup, mf, mf_e);	// LUKAS: fMCBG before
// 	calculateRatio(fTTJets, Muon, ZDecay, mp, mp_e);
// 	calculateRatio(fTTJets, Elec, SigSup, ef, ef_e);
// 	calculateRatio(fTTJets, Elec, ZDecay, ep, ep_e);
	
	calculateRatio(fMCBGMuEnr, Muon, SigSup, mf, mf_e);	// LUKAS: fMCBG before
 	calculateRatio(fMCBGMuEnr, Muon, ZDecay, mp, mp_e);
 	calculateRatio(fMCBGEMEnr, Elec, SigSup, ef, ef_e);
 	calculateRatio(fMCBGEMEnr, Elec, ZDecay, ep, ep_e);
	
//	calculateRatio(fMuEnr, Muon, SigSup, mf, mf_e);	// LUKAS: fMCBG before
// 	calculateRatio(fMCBGMuEnr, Muon, ZDecay, mp, mp_e);
// 	calculateRatio(fEMEnr, Elec, SigSup, ef, ef_e);
// 	calculateRatio(fMCBGEMEnr, Elec, ZDecay, ep, ep_e);
	
//	calculateRatio(fMuEnr, Muon, SigSup, mf, mf_e);	// LUKAS: fMCBG before
// 	calculateRatio(fDYJets, Muon, ZDecay, mp, mp_e);
// 	calculateRatio(fEMEnr, Elec, SigSup, ef, ef_e);
// 	calculateRatio(fDYJets, Elec, ZDecay, ep, ep_e);
	
//	calculateRatio(fMCBG, Muon, SigSup, mf, mf_e);	// LUKAS: fMCBG before
// 	calculateRatio(fMCBG, Muon, ZDecay, mp, mp_e);
// 	calculateRatio(fMCBG, Elec, SigSup, ef, ef_e);
// 	calculateRatio(fMCBG, Elec, ZDecay, ep, ep_e);
	
//	calculateRatio(fMuData, Muon, SigSup, mf, mf_e);	// LUKAS: fMCBG before
// 	calculateRatio(fMuData, Muon, ZDecay, mp, mp_e);
// 	calculateRatio(fEGData, Elec, SigSup, ef, ef_e);
// 	calculateRatio(fEGData, Elec, ZDecay, ep, ep_e);
	
 	///////////////////////////////////////////////////////////////////////////////////
 	// OBSERVATIONS ///////////////////////////////////////////////////////////////////
 	///////////////////////////////////////////////////////////////////////////////////
 	vector<float> ntt_mm, ntl_mm, nll_mm;
 	vector<float> ntt_ee, ntl_ee, nll_ee;
 	vector<float> ntt_em, ntl_em, nlt_em, nll_em;
 	
 	vector<float> ntt_mm_e1, ntt_ee_e1, ntt_em_e1; // squared stat. errors
	
 	vector<float> npp_mm, npf_mm, nfp_mm, nff_mm;
 	vector<float> npp_ee, npf_ee, nfp_ee, nff_ee;
 	vector<float> npp_em, npf_em, nfp_em, nff_em;
	
 	vector<float> npp_tt_mm, npf_tt_mm, nfp_tt_mm, nff_tt_mm;
 	vector<float> npp_tt_ee, npf_tt_ee, nfp_tt_ee, nff_tt_ee;
 	vector<float> npp_tt_em, npf_tt_em, nfp_tt_em, nff_tt_em;
	
	// LUKAS test
	vector<float> nt20_mm, nt10_mm, nt00_mm;
	vector<float> nt20_ee, nt10_ee, nt00_ee;
	vector<float> nt20_em, nt10_em, nt00_em;
	
 	// OS yields
 	vector<float> ntt_os_BB_ee, ntt_os_EE_ee, ntt_os_EB_ee;
 	vector<float> ntl_os_BB_ee, ntl_os_EE_ee, ntl_os_EB_ee;
 	vector<float> nlt_os_BB_ee, nlt_os_EE_ee, nlt_os_EB_ee;
	vector<float> nll_os_BB_ee, nll_os_EE_ee, nll_os_EB_ee;
 	vector<float> ntt_os_BB_em, ntt_os_EE_em;
 	vector<float> ntl_os_BB_em, ntl_os_EE_em;
 	vector<float> nlt_os_BB_em, nlt_os_EE_em;
 	vector<float> nll_os_BB_em, nll_os_EE_em;
 	// Squared stat. errors
 	vector<float> ntt_os_BB_ee_e1, ntt_os_EE_ee_e1, ntt_os_EB_ee_e1, ntt_os_BB_em_e1, ntt_os_EE_em_e1;
 	vector<float> ntl_os_BB_ee_e1, ntl_os_EE_ee_e1, ntl_os_EB_ee_e1, ntl_os_BB_em_e1, ntl_os_EE_em_e1;
 	vector<float> nlt_os_BB_ee_e1, nlt_os_EE_ee_e1, nlt_os_EB_ee_e1, nlt_os_BB_em_e1, nlt_os_EE_em_e1;
 	vector<float> nll_os_BB_ee_e1, nll_os_EE_ee_e1, nll_os_EB_ee_e1, nll_os_BB_em_e1, nll_os_EE_em_e1;
	
 	// Charge misid
	vector<float> npp_tt_cm_ee, npp_cm_ee;
 	vector<float> npp_tt_cm_em, npp_cm_em;
	
 	vector<float> scales;
 	vector<TString> names;
 	for(size_t i = 0; i < samples.size(); ++i){
 		Sample *S = fSamples[samples[i]];
 		float scale = fLumiNorm / S->getLumi();
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
		
 		ntt_mm_e1.push_back(S->getError(S->region[reg][HighPt].mm.nt20_pt->GetEntries())); // take unweighted entries
 		ntt_ee_e1.push_back(S->getError(S->region[reg][HighPt].ee.nt20_pt->GetEntries()));
 		ntt_em_e1.push_back(S->getError(S->region[reg][HighPt].em.nt20_pt->GetEntries()));
		
 		npp_mm.push_back(S->region[reg][HighPt].mm.npp_pt->GetEntries());
 		npf_mm.push_back(S->region[reg][HighPt].mm.npf_pt->GetEntries());
 		nfp_mm.push_back(S->region[reg][HighPt].mm.nfp_pt->GetEntries());
 		nff_mm.push_back(S->region[reg][HighPt].mm.nff_pt->GetEntries());
		
 		npp_em.push_back(S->region[reg][HighPt].em.npp_pt->GetEntries());
 		npf_em.push_back(S->region[reg][HighPt].em.npf_pt->GetEntries());
 		nfp_em.push_back(S->region[reg][HighPt].em.nfp_pt->GetEntries());
 		nff_em.push_back(S->region[reg][HighPt].em.nff_pt->GetEntries());
		
 		npp_ee.push_back(S->region[reg][HighPt].ee.npp_pt->GetEntries());
 		npf_ee.push_back(S->region[reg][HighPt].ee.npf_pt->GetEntries());
 		nfp_ee.push_back(S->region[reg][HighPt].ee.nfp_pt->GetEntries());
 		nff_ee.push_back(S->region[reg][HighPt].ee.nff_pt->GetEntries());
		
 		npp_tt_mm.push_back(S->region[reg][HighPt].mm.nt2pp_pt->GetEntries());
 		npf_tt_mm.push_back(S->region[reg][HighPt].mm.nt2pf_pt->GetEntries());
 		nfp_tt_mm.push_back(S->region[reg][HighPt].mm.nt2fp_pt->GetEntries());
 		nff_tt_mm.push_back(S->region[reg][HighPt].mm.nt2ff_pt->GetEntries());
		
 		npp_tt_em.push_back(S->region[reg][HighPt].em.nt2pp_pt->GetEntries());
 		npf_tt_em.push_back(S->region[reg][HighPt].em.nt2pf_pt->GetEntries());
 		nfp_tt_em.push_back(S->region[reg][HighPt].em.nt2fp_pt->GetEntries());
 		nff_tt_em.push_back(S->region[reg][HighPt].em.nt2ff_pt->GetEntries());
		
 		npp_tt_ee.push_back(S->region[reg][HighPt].ee.nt2pp_pt->GetEntries());
 		npf_tt_ee.push_back(S->region[reg][HighPt].ee.nt2pf_pt->GetEntries());
 		nfp_tt_ee.push_back(S->region[reg][HighPt].ee.nt2fp_pt->GetEntries());
 		nff_tt_ee.push_back(S->region[reg][HighPt].ee.nt2ff_pt->GetEntries());
 		
 		ntt_os_BB_em.push_back(S->region[reg][HighPt].em.nt20_OS_BB_pt->GetEntries()); // ele in barrel
 		ntt_os_EE_em.push_back(S->region[reg][HighPt].em.nt20_OS_EE_pt->GetEntries()); // ele in endcal
 		ntt_os_BB_ee.push_back(S->region[reg][HighPt].ee.nt20_OS_BB_pt->GetEntries()); // both in barrel
 		ntt_os_EE_ee.push_back(S->region[reg][HighPt].ee.nt20_OS_EE_pt->GetEntries()); // both in endcal
 		ntt_os_EB_ee.push_back(S->region[reg][HighPt].ee.nt20_OS_EB_pt->GetEntries()); // one barrel, one endcap
		
		ntl_os_BB_em.push_back(S->region[reg][HighPt].em.nt10_OS_BB_pt->GetEntries());
		ntl_os_EE_em.push_back(S->region[reg][HighPt].em.nt10_OS_EE_pt->GetEntries());
		ntl_os_BB_ee.push_back(S->region[reg][HighPt].ee.nt10_OS_BB_pt->GetEntries());
		ntl_os_EE_ee.push_back(S->region[reg][HighPt].ee.nt10_OS_EE_pt->GetEntries());
		ntl_os_EB_ee.push_back(S->region[reg][HighPt].ee.nt10_OS_EB_pt->GetEntries());
		
		nlt_os_BB_em.push_back(S->region[reg][HighPt].em.nt01_OS_BB_pt->GetEntries());
		nlt_os_EE_em.push_back(S->region[reg][HighPt].em.nt01_OS_EE_pt->GetEntries());
		nlt_os_BB_ee.push_back(S->region[reg][HighPt].ee.nt01_OS_BB_pt->GetEntries());
		nlt_os_EE_ee.push_back(S->region[reg][HighPt].ee.nt01_OS_EE_pt->GetEntries());
		nlt_os_EB_ee.push_back(S->region[reg][HighPt].ee.nt01_OS_EB_pt->GetEntries());
		
		nll_os_BB_em.push_back(S->region[reg][HighPt].em.nt00_OS_BB_pt->GetEntries());
		nll_os_EE_em.push_back(S->region[reg][HighPt].em.nt00_OS_EE_pt->GetEntries());
		nll_os_BB_ee.push_back(S->region[reg][HighPt].ee.nt00_OS_BB_pt->GetEntries());
		nll_os_EE_ee.push_back(S->region[reg][HighPt].ee.nt00_OS_EE_pt->GetEntries());
		nll_os_EB_ee.push_back(S->region[reg][HighPt].ee.nt00_OS_EB_pt->GetEntries());
 		
 		ntt_os_BB_em_e1.push_back(S->getError(S->region[reg][HighPt].em.nt20_OS_BB_pt->GetEntries()));
 		ntt_os_EE_em_e1.push_back(S->getError(S->region[reg][HighPt].em.nt20_OS_EE_pt->GetEntries()));
 		ntt_os_BB_ee_e1.push_back(S->getError(S->region[reg][HighPt].ee.nt20_OS_BB_pt->GetEntries()));
 		ntt_os_EE_ee_e1.push_back(S->getError(S->region[reg][HighPt].ee.nt20_OS_EE_pt->GetEntries()));
 		ntt_os_EB_ee_e1.push_back(S->getError(S->region[reg][HighPt].ee.nt20_OS_EB_pt->GetEntries()));
 		
 		ntl_os_BB_em_e1.push_back(S->getError(S->region[reg][HighPt].em.nt10_OS_BB_pt->GetEntries()));
 		ntl_os_EE_em_e1.push_back(S->getError(S->region[reg][HighPt].em.nt10_OS_EE_pt->GetEntries()));
 		ntl_os_BB_ee_e1.push_back(S->getError(S->region[reg][HighPt].ee.nt10_OS_BB_pt->GetEntries()));
 		ntl_os_EE_ee_e1.push_back(S->getError(S->region[reg][HighPt].ee.nt10_OS_EE_pt->GetEntries()));
 		ntl_os_EB_ee_e1.push_back(S->getError(S->region[reg][HighPt].ee.nt10_OS_EB_pt->GetEntries()));
 		
 		nlt_os_BB_em_e1.push_back(S->getError(S->region[reg][HighPt].em.nt01_OS_BB_pt->GetEntries()));
 		nlt_os_EE_em_e1.push_back(S->getError(S->region[reg][HighPt].em.nt01_OS_EE_pt->GetEntries()));
 		nlt_os_BB_ee_e1.push_back(S->getError(S->region[reg][HighPt].ee.nt01_OS_BB_pt->GetEntries()));
 		nlt_os_EE_ee_e1.push_back(S->getError(S->region[reg][HighPt].ee.nt01_OS_EE_pt->GetEntries()));
 		nlt_os_EB_ee_e1.push_back(S->getError(S->region[reg][HighPt].ee.nt01_OS_EB_pt->GetEntries()));
 		
 		nll_os_BB_em_e1.push_back(S->getError(S->region[reg][HighPt].em.nt00_OS_BB_pt->GetEntries()));
 		nll_os_EE_em_e1.push_back(S->getError(S->region[reg][HighPt].em.nt00_OS_EE_pt->GetEntries()));
 		nll_os_BB_ee_e1.push_back(S->getError(S->region[reg][HighPt].ee.nt00_OS_BB_pt->GetEntries()));
 		nll_os_EE_ee_e1.push_back(S->getError(S->region[reg][HighPt].ee.nt00_OS_EE_pt->GetEntries()));
 		nll_os_EB_ee_e1.push_back(S->getError(S->region[reg][HighPt].ee.nt00_OS_EB_pt->GetEntries()));
 		
 		npp_tt_cm_ee.push_back(scale*S->region[reg][HighPt].ee.nt2pp_cm_pt->GetEntries());
 		npp_cm_ee   .push_back(scale*S->region[reg][HighPt].ee.npp_cm_pt->GetEntries());
 		npp_tt_cm_em.push_back(scale*S->region[reg][HighPt].em.nt2pp_cm_pt->GetEntries());
 		npp_cm_em   .push_back(scale*S->region[reg][HighPt].em.npp_cm_pt->GetEntries());
		
		// LUKAS test
//		nt20_mm.push_back(S->region[reg][HighPt].mm.nt20_pt->GetEntries());
//		nt10_mm.push_back(S->region[reg][HighPt].mm.nt10_pt->GetEntries());
//		nt00_mm.push_back(S->region[reg][HighPt].mm.nt00_pt->GetEntries());
//		nt20_ee.push_back(S->region[reg][HighPt].ee.nt20_pt->GetEntries());
//		nt10_ee.push_back(S->region[reg][HighPt].ee.nt10_pt->GetEntries());
//		nt00_ee.push_back(S->region[reg][HighPt].ee.nt00_pt->GetEntries());
//		nt20_em.push_back(S->region[reg][HighPt].em.nt20_pt->GetEntries());
//		nt10_em.push_back(S->region[reg][HighPt].em.nt10_pt->GetEntries());
//		nt00_em.push_back(S->region[reg][HighPt].em.nt00_pt->GetEntries());
 	}
	
 	///////////////////////////////////////////////////////////////////////////////////
 	// PREDICTIONS ////////////////////////////////////////////////////////////////////
 	///////////////////////////////////////////////////////////////////////////////////
 	
 	// Charge MisID Predictions
 	// Abbreviations
	float fbb(0.),fee(0.),feb(0.);
	float fbbE(0.),feeE(0.),febE(0.);
 	
	calculateChMisIdProb(fMCBG, BB, fbb, fbbE);
	calculateChMisIdProb(fMCBG, EB, feb, febE);
	calculateChMisIdProb(fMCBG, EE, fee, feeE);
	
 	vector<float> ntt_cm_ee, ntt_cm_em;
 	vector<float> ntl_cm_ee, ntl_cm_em;
 	vector<float> nlt_cm_ee, nlt_cm_em;
 	vector<float> nll_cm_ee, nll_cm_em;
 	vector<float> ntt_cm_ee_e1, ntt_cm_em_e1; // squared errors, includes stat errors on yields and errors on probabilities
 	vector<float> ntl_cm_ee_e1, ntl_cm_em_e1;
 	vector<float> nlt_cm_ee_e1, nlt_cm_em_e1;
 	vector<float> nll_cm_ee_e1, nll_cm_em_e1;
 	for(size_t i = 0; i < nsamples; ++i){
 		Sample *S = fSamples[samples[i]];
		//tight-tight
 		ntt_cm_ee.push_back(2*fbb*ntt_os_BB_ee[i] + 2*fee*ntt_os_EE_ee[i] + 2*feb*ntt_os_EB_ee[i]);
 		ntt_cm_ee_e1.push_back((4*fbb*fbb* S->getError(ntt_os_BB_ee[i])) + (4*fee*fee*S->getError(ntt_os_EE_ee[i])) + 4*feb*feb*S->getError(ntt_os_EB_ee[i])
							   + (4*ntt_os_BB_ee[i]*ntt_os_BB_ee[i]*fbbE*fbbE) + (4*ntt_os_EE_ee[i]*ntt_os_EE_ee[i]*feeE*feeE) + 4*febE*febE*ntt_os_EB_ee[i]*ntt_os_EB_ee[i]);
		
 		ntt_cm_em.push_back(  fbb*ntt_os_BB_em[i] +   fee*ntt_os_EE_em[i]);
 		ntt_cm_em_e1.push_back(fbb*fbb*S->getError(ntt_os_BB_em[i]) + fee*fee*S->getError(ntt_os_EE_em[i])
							   + ntt_os_BB_em[i]*ntt_os_BB_em[i] * fbbE*fbbE + ntt_os_EE_em[i]*ntt_os_EE_em[i] * feeE*feeE);
		//tight-loose
 		ntl_cm_ee.push_back(2*fbb*ntl_os_BB_ee[i] + 2*fee*ntl_os_EE_ee[i] + 2*feb*ntl_os_EB_ee[i]);
 		ntl_cm_ee_e1.push_back((4*fbb*fbb* S->getError(ntl_os_BB_ee[i])) + (4*fee*fee*S->getError(ntl_os_EE_ee[i])) + 4*feb*feb*S->getError(ntl_os_EB_ee[i])
							   + (4*ntl_os_BB_ee[i]*ntl_os_BB_ee[i]*fbbE*fbbE) + (4*ntl_os_EE_ee[i]*ntl_os_EE_ee[i]*feeE*feeE) + 4*febE*febE*ntl_os_EB_ee[i]*ntl_os_EB_ee[i]);
		
 		ntl_cm_em.push_back(  fbb*ntl_os_BB_em[i] +   fee*ntl_os_EE_em[i]);
 		ntl_cm_em_e1.push_back(fbb*fbb*S->getError(ntl_os_BB_em[i]) + fee*fee*S->getError(ntl_os_EE_em[i])
							   + ntl_os_BB_em[i]*ntl_os_BB_em[i] * fbbE*fbbE + ntl_os_EE_em[i]*ntl_os_EE_em[i] * feeE*feeE);
		//loose-tight
 		nlt_cm_ee.push_back(2*fbb*nlt_os_BB_ee[i] + 2*fee*nlt_os_EE_ee[i] + 2*feb*nlt_os_EB_ee[i]);
 		nlt_cm_ee_e1.push_back((4*fbb*fbb* S->getError(nlt_os_BB_ee[i])) + (4*fee*fee*S->getError(nlt_os_EE_ee[i])) + 4*feb*feb*S->getError(nlt_os_EB_ee[i])
							   + (4*nlt_os_BB_ee[i]*nlt_os_BB_ee[i]*fbbE*fbbE) + (4*nlt_os_EE_ee[i]*nlt_os_EE_ee[i]*feeE*feeE) + 4*febE*febE*nlt_os_EB_ee[i]*nlt_os_EB_ee[i]);
		
 		nlt_cm_em.push_back(  fbb*nlt_os_BB_em[i] +   fee*nlt_os_EE_em[i]);
 		nlt_cm_em_e1.push_back(fbb*fbb*S->getError(nlt_os_BB_em[i]) + fee*fee*S->getError(nlt_os_EE_em[i])
							   + nlt_os_BB_em[i]*nlt_os_BB_em[i] * fbbE*fbbE + nlt_os_EE_em[i]*nlt_os_EE_em[i] * feeE*feeE);
		//loose-loose
 		nll_cm_ee.push_back(2*fbb*nll_os_BB_ee[i] + 2*fee*nll_os_EE_ee[i] + 2*feb*nll_os_EB_ee[i]);
 		nll_cm_ee_e1.push_back((4*fbb*fbb* S->getError(nll_os_BB_ee[i])) + (4*fee*fee*S->getError(nll_os_EE_ee[i])) + 4*feb*feb*S->getError(nll_os_EB_ee[i])
							   + (4*nll_os_BB_ee[i]*nll_os_BB_ee[i]*fbbE*fbbE) + (4*nll_os_EE_ee[i]*nll_os_EE_ee[i]*feeE*feeE) + 4*febE*febE*nll_os_EB_ee[i]*nll_os_EB_ee[i]);
		
 		nll_cm_em.push_back(  fbb*nll_os_BB_em[i] +   fee*nll_os_EE_em[i]);
 		nll_cm_em_e1.push_back(fbb*fbb*S->getError(nll_os_BB_em[i]) + fee*fee*S->getError(nll_os_EE_em[i])
							   + nll_os_BB_em[i]*nll_os_BB_em[i] * fbbE*fbbE + nll_os_EE_em[i]*nll_os_EE_em[i] * feeE*feeE);
 	}
	
	// Fake Predictions
 	FakeRatios *FR = new FakeRatios();
 	FR->setIsMC(true);
 	FR->setNToyMCs(100);
 	FR->setAddESyst(0.0);
	
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
		
		// subtract charge misID before predicting fakes (works only for flat ratios)
// 		FR->setMMNtl(ntt_mm[i]                , ntl_mm[i]                                           , nll_mm[i]                );
// 		FR->setEENtl(ntt_ee[i] - ntt_cm_ee[i] , ntl_ee[i] - ntl_cm_ee[i] - nlt_cm_ee[i]             , nll_ee[i] - nll_cm_ee[i] );
// 		FR->setEMNtl(ntt_em[i] - ntt_cm_em[i] , ntl_em[i] - ntl_cm_em[i] , nlt_em[i] - nlt_cm_em[i] , nll_em[i] - nll_cm_em[i] );
		
 		Sample *S = fSamples[samples[i]];
 		FR->setNGen(S->ngen);
		
 		//// Differential Ratios
 		npp_pred_mm.push_back(S->numbers[reg][Muon].npp);
 		npf_pred_mm.push_back(S->numbers[reg][Muon].npf + S->numbers[reg][Muon].nfp);
 		nff_pred_mm.push_back(S->numbers[reg][Muon].nff);
 		npp_pred_ee.push_back(S->numbers[reg][Elec].npp);
 		npf_pred_ee.push_back(S->numbers[reg][Elec].npf + S->numbers[reg][Elec].nfp);
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
		
 		//// Flat Ratios
// 		npp_pred_mm   .push_back(FR->getMMNpp());
// 		npp_pred_mm_e1.push_back(FR->getMMNppEStat());
// 		npf_pred_mm   .push_back(FR->getMMNpf());
// 		npf_pred_mm_e1.push_back(FR->getMMNpfEStat());
// 		nff_pred_mm   .push_back(FR->getMMNff());
// 		nff_pred_mm_e1.push_back(FR->getMMNffEStat());
// 		nF_pred_mm_e1 .push_back(FR->getMMTotEStat());
// 		
// 		npp_pred_ee   .push_back(FR->getEENpp());
// 		npp_pred_ee_e1.push_back(FR->getEENppEStat());
// 		npf_pred_ee   .push_back(FR->getEENpf());
// 		npf_pred_ee_e1.push_back(FR->getEENpfEStat());
// 		nff_pred_ee   .push_back(FR->getEENff());
// 		nff_pred_ee_e1.push_back(FR->getEENffEStat());
// 		nF_pred_ee_e1 .push_back(FR->getEETotEStat());
// 		
// 		npp_pred_em   .push_back(FR->getEMNpp());
// 		npp_pred_em_e1.push_back(FR->getEMNppEStat());
// 		npf_pred_em   .push_back(FR->getEMNpf());
// 		npf_pred_em_e1.push_back(FR->getEMNpfEStat());
// 		nfp_pred_em   .push_back(FR->getEMNfp());
// 		nfp_pred_em_e1.push_back(FR->getEMNfpEStat());
// 		nff_pred_em   .push_back(FR->getEMNff());
// 		nff_pred_em_e1.push_back(FR->getEMNffEStat());
// 		nF_pred_em_e1 .push_back(FR->getEMTotEStat());
 	}
	
 	// Sums
 	float ntt_sum_mm(0.), ntl_sum_mm(0.), nll_sum_mm(0.);
 	float ntt_sum_em(0.), ntl_sum_em(0.), nlt_sum_em(0.), nll_sum_em(0.);
 	float ntt_sum_ee(0.), ntl_sum_ee(0.), nll_sum_ee(0.);
	
 	float ntt_sum_mm_e1(0.), ntt_sum_em_e1(0.), ntt_sum_ee_e1(0.);
	
 	float npp_sum_mm(0.), npf_sum_mm(0.), nff_sum_mm(0.);
 	float npp_sum_em(0.), npf_sum_em(0.), nfp_sum_em(0.), nff_sum_em(0.);
 	float npp_sum_ee(0.), npf_sum_ee(0.), nff_sum_ee(0.);
	
 	float npp_pred_sum_mm(0.), npf_pred_sum_mm(0.), nff_pred_sum_mm(0.);
 	float npp_pred_sum_em(0.), npf_pred_sum_em(0.), nfp_pred_sum_em(0.), nff_pred_sum_em(0.);
 	float npp_pred_sum_ee(0.), npf_pred_sum_ee(0.), nff_pred_sum_ee(0.);
	
 	float ntt_cm_sum_ee(0.), ntt_cm_sum_em(0.);
 	float ntt_cm_sum_ee_e1(0.), ntt_cm_sum_em_e1(0.);
	float ntl_cm_sum_ee(0.), ntl_cm_sum_em(0.);
 	float ntl_cm_sum_ee_e1(0.), ntl_cm_sum_em_e1(0.);
	float nlt_cm_sum_ee(0.), nlt_cm_sum_em(0.);
 	float nlt_cm_sum_ee_e1(0.), nlt_cm_sum_em_e1(0.);
	float nll_cm_sum_ee(0.), nll_cm_sum_em(0.);
 	float nll_cm_sum_ee_e1(0.), nll_cm_sum_em_e1(0.);
	
 	float npp_cm_sum_ee(0.),    npp_cm_sum_em(0.);
 	float npp_tt_cm_sum_ee(0.), npp_tt_cm_sum_em(0.);
	
	float npp_tt_sum_mm(0.), npf_tt_sum_mm(0.), nff_tt_sum_mm(0.);
	float npp_tt_sum_em(0.), npf_tt_sum_em(0.), nfp_tt_sum_em(0.), nff_tt_sum_em(0.);
	float npp_tt_sum_ee(0.), npf_tt_sum_ee(0.), nff_tt_sum_ee(0.);
	
 	///////////////////////////////////////////
 	// Rare SM
// LUKAS 	float ntt_rare_mm(0.), ntt_rare_em(0.), ntt_rare_ee(0.);
// LUKAS 	float ntt_rare_mm_e1(0.), ntt_rare_em_e1(0.), ntt_rare_ee_e1(0.); // squared stat errors
// LUKAS 	for(size_t i = 0; i < fMCRareSM.size(); ++i){
// LUKAS 		Sample *S = fSamples[fMCRareSM[i]];
// LUKAS 		float scale = fLumiNorm/S->getLumi();
// LUKAS 		ntt_rare_mm += scale*S->numbers[reg][Muon].nt2;
// LUKAS 		ntt_rare_em += scale*S->numbers[reg][ElMu].nt2;
// LUKAS 		ntt_rare_ee += scale*S->numbers[reg][Elec].nt2;
// LUKAS 		ntt_rare_mm_e1 += scale*scale*pow(S->getError(S->region[reg][HighPt].mm.nt20_pt->GetEntries()),2);
// LUKAS 		ntt_rare_em_e1 += scale*scale*pow(S->getError(S->region[reg][HighPt].em.nt20_pt->GetEntries()),2);
// LUKAS 		ntt_rare_ee_e1 += scale*scale*pow(S->getError(S->region[reg][HighPt].ee.nt20_pt->GetEntries()),2);
// LUKAS 	}
	
 	///////////////////////////////////////////
 	// WZ production
// LUKAS 	float wzscale = fLumiNorm/fSamples[WZ]->getLumi();
// LUKAS 	float ntt_wz_mm = wzscale*fSamples[WZ]->numbers[reg][Muon].nt2;
// LUKAS 	float ntt_wz_em = wzscale*fSamples[WZ]->numbers[reg][ElMu].nt2;
// LUKAS 	float ntt_wz_ee = wzscale*fSamples[WZ]->numbers[reg][Elec].nt2;
// LUKAS
// LUKAS 	float ntt_wz_mm_e1 = wzscale*wzscale*pow(fSamples[WZ]->getError(fSamples[WZ]->region[reg][HighPt].mm.nt20_pt->GetEntries()),2); // for stat error take actual entries, not pileup weighted integral...
// LUKAS 	float ntt_wz_em_e1 = wzscale*wzscale*pow(fSamples[WZ]->getError(fSamples[WZ]->region[reg][HighPt].em.nt20_pt->GetEntries()),2);
// LUKAS 	float ntt_wz_ee_e1 = wzscale*wzscale*pow(fSamples[WZ]->getError(fSamples[WZ]->region[reg][HighPt].ee.nt20_pt->GetEntries()),2);
	
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
 		
 		ntt_cm_sum_ee    += scales[i] * ntt_cm_ee[i];
 		ntt_cm_sum_em    += scales[i] * ntt_cm_em[i];
		ntl_cm_sum_ee    += scales[i] * ntl_cm_ee[i];
 		ntl_cm_sum_em    += scales[i] * ntl_cm_em[i];
		nlt_cm_sum_ee    += scales[i] * nlt_cm_ee[i];
 		nlt_cm_sum_em    += scales[i] * nlt_cm_em[i];
		nll_cm_sum_ee    += scales[i] * nll_cm_ee[i];
 		nll_cm_sum_em    += scales[i] * nll_cm_em[i];
		
 		ntt_cm_sum_ee_e1 += scales[i]*scales[i] * ntt_cm_ee_e1[i];
 		ntt_cm_sum_em_e1 += scales[i]*scales[i] * ntt_cm_em_e1[i];
		
 		ntl_cm_sum_ee_e1 += scales[i]*scales[i] * ntl_cm_ee_e1[i];
 		ntl_cm_sum_em_e1 += scales[i]*scales[i] * ntl_cm_em_e1[i];
		
 		nlt_cm_sum_ee_e1 += scales[i]*scales[i] * nlt_cm_ee_e1[i];
 		nlt_cm_sum_em_e1 += scales[i]*scales[i] * nlt_cm_em_e1[i];
		
 		nll_cm_sum_ee_e1 += scales[i]*scales[i] * nll_cm_ee_e1[i];
 		nll_cm_sum_em_e1 += scales[i]*scales[i] * nll_cm_em_e1[i];
		
 		npp_cm_sum_ee    += scales[i] * npp_cm_ee[i];
 		npp_cm_sum_em    += scales[i] * npp_cm_em[i];
 		npp_tt_cm_sum_ee += scales[i] * npp_tt_cm_ee[i];
 		npp_tt_cm_sum_em += scales[i] * npp_tt_cm_em[i];
 		
 		ntt_sum_mm_e1    += scales[i]*scales[i] * ntt_mm_e1[i];
 		ntt_sum_em_e1    += scales[i]*scales[i] * ntt_em_e1[i];
 		ntt_sum_ee_e1    += scales[i]*scales[i] * ntt_ee_e1[i];
		
		npp_tt_sum_mm	+= scales[i] * npp_tt_mm[i];
		npf_tt_sum_mm	+= scales[i] * npf_tt_mm[i] + scales[i] * nfp_tt_mm[i];
		nff_tt_sum_mm	+= scales[i] * nff_tt_mm[i];
		npp_tt_sum_em	+= scales[i] * npp_tt_em[i];
		npf_tt_sum_em	+= scales[i] * npf_tt_em[i];
		nfp_tt_sum_em	+= scales[i] * nfp_tt_em[i];
		nff_tt_sum_em	+= scales[i] * nff_tt_em[i];
		npp_tt_sum_ee	+= scales[i] * npp_tt_ee[i];
		npf_tt_sum_ee	+= scales[i] * npf_tt_ee[i] + scales[i] * nfp_tt_ee[i];
		nff_tt_sum_ee	+= scales[i] * nff_tt_ee[i];
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
 		OUT << setw(7)  << Form("%6.3f", scales[i]*npp_tt_mm[i]) << " | ";
		OUT << setw(7)  << Form("%6.3f", scales[i]*(npf_tt_mm[i]+nfp_tt_mm[i])) << " | ";
		OUT << setw(7)  << Form("%6.3f", scales[i]*nff_tt_mm[i]) << " || ";
		OUT << setw(7)  << Form("%6.3f", scales[i]*npp_tt_em[i]) << " | ";
		OUT << setw(7)  << Form("%6.3f", scales[i]*npf_tt_em[i]) << " | ";
		OUT << setw(7)  << Form("%6.3f", scales[i]*nfp_tt_em[i]) << " | ";
		OUT << setw(7)  << Form("%6.3f", scales[i]*nff_tt_em[i]) << " || ";
		OUT << setw(7)  << Form("%6.3f", scales[i]*npp_tt_ee[i]) << " | ";
		OUT << setw(7)  << Form("%6.3f", scales[i]*(npf_tt_ee[i]+nfp_tt_ee[i])) << " | ";
		OUT << setw(7)  << Form("%6.3f", scales[i]*nff_tt_ee[i]) << " || ";
 		OUT << endl;
 	}
 	OUT << "--------------------------------------------------------------------------------------------------------------------------" << endl;
 	OUT << setw(16) << "Sum"  << " || ";
	OUT << setw(7) << Form("%6.3f", npp_tt_sum_mm) << " | ";
 	OUT << setw(7) << Form("%6.3f", npf_tt_sum_mm) << " | ";
 	OUT << setw(7) << Form("%6.3f", nff_tt_sum_mm) << " || ";
 	OUT << setw(7) << Form("%6.3f", npp_tt_sum_em) << " | ";
 	OUT << setw(7) << Form("%6.3f", npf_tt_sum_em) << " | ";
 	OUT << setw(7) << Form("%6.3f", nfp_tt_sum_em) << " | ";
 	OUT << setw(7) << Form("%6.3f", nff_tt_sum_em) << " || ";
 	OUT << setw(7) << Form("%6.3f", npp_tt_sum_ee) << " | ";
 	OUT << setw(7) << Form("%6.3f", npf_tt_sum_ee) << " | ";
 	OUT << setw(7) << Form("%6.3f", nff_tt_sum_ee) << " || ";
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
	
	OUT << "===========================================================================================" << endl;
 	OUT << "                 ||                   E/Mu                ||             E/E             ||" << endl;
 	OUT << "    CHARGE MISID ||   Ntt   |   Ntl   |   Nlt   |   Nll   ||   Ntt   |   Ntl   |   Nll   ||" << endl;
 	OUT << "-------------------------------------------------------------------------------------------" << endl;
 	for(size_t i = 0; i < nsamples; ++i){
 		OUT << setw(16) << names[i] << " || ";
 		OUT << setw(7)  << Form("%6.3f", scales[i]*ntt_cm_em[i]) << " | ";
 		OUT << setw(7)  << Form("%6.3f", scales[i]*ntl_cm_em[i]) << " | ";
 		OUT << setw(7)  << Form("%6.3f", scales[i]*nlt_cm_em[i]) << " | ";
 		OUT << setw(7)  << Form("%6.3f", scales[i]*nll_cm_em[i]) << " || ";
 		OUT << setw(7)  << Form("%6.3f", scales[i]*ntt_cm_ee[i]) << " | ";
 		OUT << setw(7)  << Form("%6.3f", scales[i]*(ntl_cm_ee[i]+nlt_cm_ee[i])) << " | ";
 		OUT << setw(7)  << Form("%6.3f", scales[i]*nll_cm_ee[i]) << " || ";
 		OUT << endl;
 	}
 	OUT << "-------------------------------------------------------------------------------------------" << endl;
 	OUT << setw(16) << "Sum"  << " || ";
 	OUT << setw(7) << Form("%6.3f", ntt_cm_sum_em) << " | ";
 	OUT << setw(7) << Form("%6.3f", ntl_cm_sum_em) << " | ";
 	OUT << setw(7) << Form("%6.3f", nlt_cm_sum_em) << " | ";
 	OUT << setw(7) << Form("%6.3f", nll_cm_sum_em) << " || ";
 	OUT << setw(7) << Form("%6.3f", ntt_cm_sum_ee) << " | ";
 	OUT << setw(7) << Form("%6.3f", ntl_cm_sum_ee+nlt_cm_sum_ee) << " | ";
 	OUT << setw(7) << Form("%6.3f", nll_cm_sum_ee) << " || ";
 	OUT << endl;
// 	OUT << setw(16) << "Channels sum"  << " || ";
// 	OUT << Form("                     %6.3f ||                               %6.3f ||                      %6.3f || ",
//				ntt_sum_mm+ntl_sum_mm+nll_sum_mm, ntt_sum_em+ntl_sum_em+nlt_sum_em+nll_sum_em, ntt_sum_ee+ntl_sum_ee+nll_sum_ee) << endl;
 	OUT << "===========================================================================================" << endl;
 	OUT << endl;
	
 	OUT << "===================================================================================================" << endl;
 	OUT << "                 ||       Mu/Mu       |       E/Mu        |        E/E        |        Sum        |" << endl;
 	OUT << "---------------------------------------------------------------------------------------------------" << endl;
 	float fakesum = npf_pred_sum_mm+npf_pred_sum_em+nfp_pred_sum_em+npf_pred_sum_ee+nff_pred_sum_mm+nff_pred_sum_em+nff_pred_sum_ee;
 	float fakesum_e1 = nF_pred_sum_mm_e1+nF_pred_sum_em_e1+nF_pred_sum_ee_e1;
 	OUT << Form(" Fakes           || %7.2f\\pm%7.2f & %7.2f\\pm%7.2f & %7.2f\\pm%7.2f & %7.2f\\pm%7.2f |",
				npf_pred_sum_mm+nff_pred_sum_mm,                 sqrt(nF_pred_sum_mm_e1),
				npf_pred_sum_em+nfp_pred_sum_em+nff_pred_sum_em, sqrt(nF_pred_sum_em_e1),
				npf_pred_sum_ee+nff_pred_sum_ee,                 sqrt(nF_pred_sum_ee_e1),
				fakesum, sqrt(fakesum_e1)) << endl;
 	OUT << Form(" Charge MisID    ||                   & %7.2f\\pm%7.2f & %7.2f\\pm%7.2f & %7.2f\\pm%7.2f |",
				ntt_cm_sum_em,               sqrt(ntt_cm_sum_em_e1),
				ntt_cm_sum_ee,               sqrt(ntt_cm_sum_ee_e1),
				ntt_cm_sum_ee+ntt_cm_sum_em, sqrt(ntt_cm_sum_ee_e1+ntt_cm_sum_em_e1)) << endl;
// LUKAS 	OUT << Form(" Irreducible     || %7.2f\\pm%7.2f & %7.2f\\pm%7.2f & %7.2f\\pm%7.2f & %7.2f\\pm%7.2f |",
// LUKAS 	ntt_rare_mm,                         sqrt(ntt_rare_mm_e1),
// LUKAS 	ntt_rare_em,                         sqrt(ntt_rare_em_e1),
// LUKAS 	ntt_rare_ee,                         sqrt(ntt_rare_ee_e1),
// LUKAS 	ntt_rare_mm+ntt_rare_em+ntt_rare_ee, sqrt(ntt_rare_mm_e1+ntt_rare_em_e1+ntt_rare_ee_e1)) << endl;
// LUKAS 	OUT << Form(" WZ Production   || %7.2f\\pm%7.2f & %7.2f\\pm%7.2f & %7.2f\\pm%7.2f & %7.2f\\pm%7.2f |",
// LUKAS 	ntt_wz_mm,                     sqrt(ntt_wz_mm_e1),
// LUKAS 	ntt_wz_em,                     sqrt(ntt_wz_em_e1),
// LUKAS 	ntt_wz_ee,                     sqrt(ntt_wz_ee_e1),
// LUKAS 	ntt_wz_mm+ntt_wz_em+ntt_wz_ee, sqrt(ntt_wz_mm_e1+ntt_wz_em_e1+ntt_wz_ee_e1)) << endl;
	
	float tot_bg_mm = npf_pred_sum_mm                +nff_pred_sum_mm              ;
 	float tot_bg_em = npf_pred_sum_em+nfp_pred_sum_em+nff_pred_sum_em+ntt_cm_sum_em;
 	float tot_bg_ee = npf_pred_sum_ee                +nff_pred_sum_ee+ntt_cm_sum_ee;
 	float tot_bg = tot_bg_mm + tot_bg_em + tot_bg_ee;
 	float tot_bg_mm_e1 = nF_pred_sum_mm_e1                   ;
 	float tot_bg_em_e1 = nF_pred_sum_em_e1 + ntt_cm_sum_em_e1;
 	float tot_bg_ee_e1 = nF_pred_sum_ee_e1 + ntt_cm_sum_ee_e1;
 	float tot_bg_e1    = tot_bg_mm_e1 + tot_bg_em_e1 + tot_bg_ee_e1;
// LUKAS 	float tot_bg_mm = npf_pred_sum_mm                +nff_pred_sum_mm              +ntt_rare_mm + ntt_wz_mm;
// LUKAS 	float tot_bg_em = npf_pred_sum_em+nfp_pred_sum_em+nff_pred_sum_em+ntt_cm_sum_em+ntt_rare_em + ntt_wz_em;
// LUKAS 	float tot_bg_ee = npf_pred_sum_ee                +nff_pred_sum_ee+ntt_cm_sum_ee+ntt_rare_ee + ntt_wz_ee;
// LUKAS 	float tot_bg = tot_bg_mm + tot_bg_em + tot_bg_ee;
// LUKAS 	float tot_bg_mm_e1 = nF_pred_sum_mm_e1                    + ntt_rare_mm_e1 + ntt_wz_mm_e1;
// LUKAS 	float tot_bg_em_e1 = nF_pred_sum_em_e1 + ntt_cm_sum_em_e1 + ntt_rare_em_e1 + ntt_wz_em_e1;
// LUKAS 	float tot_bg_ee_e1 = nF_pred_sum_ee_e1 + ntt_cm_sum_ee_e1 + ntt_rare_ee_e1 + ntt_wz_ee_e1;
// LUKAS 	float tot_bg_e1    = tot_bg_mm_e1 + tot_bg_em_e1 + tot_bg_ee_e1;
 	float ntt_sum    = ntt_sum_mm+ntt_sum_em+ntt_sum_ee;
 	float ntt_sum_e1 = ntt_sum_mm_e1+ntt_sum_em_e1+ntt_sum_ee_e1;
 	OUT << Form(" Total Pred.     || %7.2f\\pm%7.2f & %7.2f\\pm%7.2f & %7.2f\\pm%7.2f & %7.2f\\pm%7.2f |",
				tot_bg_mm, sqrt(tot_bg_mm_e1),
				tot_bg_em, sqrt(tot_bg_em_e1),
				tot_bg_ee, sqrt(tot_bg_ee_e1),
				tot_bg,    sqrt(tot_bg_e1) ) << endl;
 	OUT << Form(" Observed        || %7.2f\\pm%7.2f & %7.2f\\pm%7.2f & %7.2f\\pm%7.2f & %7.2f\\pm%7.2f |",
				ntt_sum_mm, sqrt(ntt_sum_mm_e1),
				ntt_sum_em, sqrt(ntt_sum_em_e1),
				ntt_sum_ee, sqrt(ntt_sum_ee_e1),
				ntt_sum,    sqrt(ntt_sum_e1)) << endl;
	OUT << Form(" Pred./Obs.      || %7.2f\\pm%7.2f & %7.2f\\pm%7.2f & %7.2f\\pm%7.2f & %7.2f\\pm%7.2f |",
				tot_bg_mm/ntt_sum_mm,
//				0.,
				sqrt(nF_pred_sum_mm_e1                  /*+  ntt_rare_mm_e1 * pow(((ntt_sum_mm-tot_bg_mm)/tot_bg_mm),2)*/ + (tot_bg_mm_e1 /*- ntt_rare_mm_e1*/) * pow(tot_bg_mm/ntt_sum_mm,2)  )/ntt_sum_mm,
				
				tot_bg_em/ntt_sum_em,
//				0.,
				sqrt(nF_pred_sum_em_e1+ntt_cm_sum_em_e1 /*+  ntt_rare_em_e1 * pow(((ntt_sum_em-tot_bg_em)/tot_bg_em),2)*/ + (tot_bg_em_e1 /*- ntt_rare_em_e1*/) * pow(tot_bg_em/ntt_sum_em,2)  )/ntt_sum_em,
				
				tot_bg_ee/ntt_sum_ee,
//				0.,
				sqrt(nF_pred_sum_ee_e1+ntt_cm_sum_ee_e1 +  /*ntt_rare_ee_e1 * pow(((ntt_sum_ee-tot_bg_ee)/tot_bg_ee),2)*/ + (tot_bg_ee_e1 /*- ntt_rare_ee_e1*/) * pow(tot_bg_ee/ntt_sum_ee,2)  )/ntt_sum_ee,
				
				tot_bg/ntt_sum,
//				0.
				sqrt(fakesum_e1 + ntt_cm_sum_ee_e1+ntt_cm_sum_em_e1 +  /*(ntt_rare_mm_e1+ntt_rare_em_e1+ntt_rare_ee_e1) * pow(((ntt_sum-tot_bg)/tot_bg),2)*/ + (tot_bg_e1 /*- (ntt_rare_mm_e1+ntt_rare_em_e1+ntt_rare_ee_e1)*/) * pow(tot_bg/ntt_sum,2)  )/ntt_sum
				) << endl;
// 	OUT << Form(" Pred./Obs.      || %7.2f\\pm%7.2f & %7.2f\\pm%7.2f & %7.2f\\pm%7.2f & %7.2f\\pm%7.2f |",
//				tot_bg_mm/ntt_sum_mm,
//				sqrt(nF_pred_sum_mm_e1                  +  ntt_rare_mm_e1 * pow(((ntt_sum_mm-tot_bg_mm)/tot_bg_mm),2) + (tot_bg_mm_e1 - ntt_rare_mm_e1) * pow(tot_bg_mm/ntt_sum_mm,2)  )/ntt_sum_mm,
//				tot_bg_em/ntt_sum_em,
//				sqrt(nF_pred_sum_em_e1+ntt_cm_sum_em_e1 +  ntt_rare_em_e1 * pow(((ntt_sum_em-tot_bg_em)/tot_bg_em),2) + (tot_bg_em_e1 - ntt_rare_em_e1) * pow(tot_bg_em/ntt_sum_em,2)  )/ntt_sum_em,
//				tot_bg_ee/ntt_sum_ee,
//				sqrt(nF_pred_sum_ee_e1+ntt_cm_sum_ee_e1 +  ntt_rare_ee_e1 * pow(((ntt_sum_ee-tot_bg_ee)/tot_bg_ee),2) + (tot_bg_ee_e1 - ntt_rare_ee_e1) * pow(tot_bg_ee/ntt_sum_ee,2)  )/ntt_sum_ee,
//				tot_bg/ntt_sum,
//				sqrt(fakesum_e1 + ntt_cm_sum_ee_e1+ntt_cm_sum_em_e1 +  (ntt_rare_mm_e1+ntt_rare_em_e1+ntt_rare_ee_e1) * pow(((ntt_sum-tot_bg)/tot_bg),2) + (tot_bg_e1 - (ntt_rare_mm_e1+ntt_rare_em_e1+ntt_rare_ee_e1)) * pow(tot_bg/ntt_sum,2)  )/ntt_sum
//				) << endl;
	OUT << Form(" Pred.-Obs./Pred || %7.2f\\pm%7.2f & %7.2f\\pm%7.2f & %7.2f\\pm%7.2f & %7.2f\\pm%7.2f |",
				(tot_bg_mm-ntt_sum_mm)/tot_bg_mm,
//				0.,
				sqrt(nF_pred_sum_mm_e1*pow(ntt_sum_mm/tot_bg_mm,2) /*+ ntt_rare_mm_e1*(pow((ntt_sum_mm-tot_bg_mm)/tot_bg_mm,2))*/ + ntt_sum_mm_e1 /*- ntt_rare_mm_e1*/)/tot_bg_mm,
				
				(tot_bg_em-ntt_sum_em)/tot_bg_em,
//				0.,
				sqrt((nF_pred_sum_em_e1+ntt_cm_sum_em_e1)*pow(ntt_sum_em/tot_bg_em,2) + /*ntt_rare_em_e1*(pow((ntt_sum_em-tot_bg_em)/tot_bg_em,2))*/ + ntt_sum_em_e1 /*- ntt_rare_em_e1*/)/tot_bg_em,
				
				(tot_bg_ee-ntt_sum_ee)/tot_bg_ee,
//				0.,
				sqrt((nF_pred_sum_ee_e1+ntt_cm_sum_ee_e1)*pow(ntt_sum_ee/tot_bg_ee,2) + /*ntt_rare_ee_e1*(pow((ntt_sum_ee-tot_bg_ee)/tot_bg_ee,2))*/ + ntt_sum_ee_e1 /*- ntt_rare_ee_e1*/)/tot_bg_ee,
				
				(tot_bg-ntt_sum)/tot_bg,
//				0.
				sqrt((fakesum_e1+ntt_cm_sum_ee_e1+ntt_cm_sum_em_e1)*pow(ntt_sum/tot_bg,2) /*+ (ntt_rare_mm_e1+ntt_rare_em_e1+ntt_rare_ee_e1)*(pow((ntt_sum-tot_bg)/tot_bg,2))*/ + ntt_sum_e1 /*- (ntt_rare_mm_e1+ntt_rare_em_e1+ntt_rare_ee_e1)*/)/tot_bg
				) << endl;
// 	OUT << Form(" Pred.-Obs./Pred || %7.2f\\pm%7.2f & %7.2f\\pm%7.2f & %7.2f\\pm%7.2f & %7.2f\\pm%7.2f |",
//				(tot_bg_mm-ntt_sum_mm)/tot_bg_mm,
//				sqrt(nF_pred_sum_mm_e1*pow(ntt_sum_mm/tot_bg_mm,2) + ntt_rare_mm_e1*(pow((ntt_sum_mm-tot_bg_mm)/tot_bg_mm,2)) + ntt_sum_mm_e1 - ntt_rare_mm_e1)/tot_bg_mm,
//				
//				(tot_bg_em-ntt_sum_em)/tot_bg_em,
//				sqrt((nF_pred_sum_em_e1+ntt_cm_sum_em_e1)*pow(ntt_sum_em/tot_bg_em,2) + ntt_rare_em_e1*(pow((ntt_sum_em-tot_bg_em)/tot_bg_em,2)) + ntt_sum_em_e1 - ntt_rare_em_e1)/tot_bg_em,
//				
//				(tot_bg_ee-ntt_sum_ee)/tot_bg_ee,
//				sqrt((nF_pred_sum_ee_e1+ntt_cm_sum_ee_e1)*pow(ntt_sum_ee/tot_bg_ee,2) + ntt_rare_ee_e1*(pow((ntt_sum_ee-tot_bg_ee)/tot_bg_ee,2)) + ntt_sum_ee_e1 - ntt_rare_ee_e1)/tot_bg_ee,
//				
//				(tot_bg-ntt_sum)/tot_bg,
//				sqrt((fakesum_e1+ntt_cm_sum_ee_e1+ntt_cm_sum_em_e1)*pow(ntt_sum/tot_bg,2) + (ntt_rare_mm_e1+ntt_rare_em_e1+ntt_rare_ee_e1)*(pow((ntt_sum-tot_bg)/tot_bg,2)) + ntt_sum_e1 - (ntt_rare_mm_e1+ntt_rare_em_e1+ntt_rare_ee_e1))/tot_bg
//				) << endl;
	OUT << "===================================================================================================" << endl;
 	OUT << endl;
	
 	// OUT << "==========================================================================================================================" << endl;
 	// OUT << "                 ||            Mu/Mu            ||                   E/Mu                ||             E/E             ||" << endl;
 	// OUT << "     PRED. IN TT ||   Npp   |   Nfp   |   Nff   ||   Npp   |   Npf   |   Nfp   |   Nff   ||   Npp   |   Npf   |   Nff   ||" << endl;
 	// OUT << "--------------------------------------------------------------------------------------------------------------------------" << endl;
 	// OUT << setw(16) << "Npf Truth"  << " || ";
 	// OUT << setw(7) << Form("%6.3f", mp*mp*npp_sum_mm) << " | ";
 	// OUT << setw(7) << Form("%6.3f", mp*mf*npf_sum_mm) << " | ";
 	// OUT << setw(7) << Form("%6.3f", mf*mf*nff_sum_mm) << " || ";
 	// OUT << setw(7) << Form("%6.3f", mp*ep*npp_sum_em) << " | ";
 	// OUT << setw(7) << Form("%6.3f", mp*ef*npf_sum_em) << " | ";
 	// OUT << setw(7) << Form("%6.3f", mf*ep*nfp_sum_em) << " | ";
 	// OUT << setw(7) << Form("%6.3f", mf*ef*nff_sum_em) << " || ";
 	// OUT << setw(7) << Form("%6.3f", ep*ep*npp_sum_ee) << " | ";
 	// OUT << setw(7) << Form("%6.3f", ep*ef*npf_sum_ee) << " | ";
 	// OUT << setw(7) << Form("%6.3f", ef*ef*nff_sum_ee) << " || ";
 	// OUT << endl;
 	// OUT << "--------------------------------------------------------------------------------------------------------------------------" << endl;
 	// OUT << setw(16) << "FR Prediction"  << " || ";
 	// OUT << setw(7) << Form("%6.3f", npp_pred_sum_mm) << " | ";
 	// OUT << setw(7) << Form("%6.3f", npf_pred_sum_mm) << " | ";
 	// OUT << setw(7) << Form("%6.3f", nff_pred_sum_mm) << " || ";
 	// OUT << setw(7) << Form("%6.3f", npp_pred_sum_em) << " | ";
 	// OUT << setw(7) << Form("%6.3f", npf_pred_sum_em) << " | ";
 	// OUT << setw(7) << Form("%6.3f", nfp_pred_sum_em) << " | ";
 	// OUT << setw(7) << Form("%6.3f", nff_pred_sum_em) << " || ";
 	// OUT << setw(7) << Form("%6.3f", npp_pred_sum_ee) << " | ";
 	// OUT << setw(7) << Form("%6.3f", npf_pred_sum_ee) << " | ";
 	// OUT << setw(7) << Form("%6.3f", nff_pred_sum_ee) << " || ";
 	// OUT << endl;
 	// OUT << "==========================================================================================================================" << endl;
 	// OUT << setw(16) << "Pred. Fakes"  << " || ";
 	// OUT << setw(7) << Form("%6.3f", npf_pred_sum_mm+nff_pred_sum_mm) << " |                   || ";
 	// OUT << setw(7) << Form("%6.3f", npf_pred_sum_em+nfp_pred_sum_em+nff_pred_sum_em) << " |                             || ";
 	// OUT << setw(7) << Form("%6.3f", npf_pred_sum_ee+nff_pred_sum_ee) << " |                   || ";
 	// OUT << endl;
 	// OUT << setw(16) << "Pred. Ch-MID"  << " ||         |                   || ";
 	// OUT << setw(7) << Form("%6.3f", ntt_cm_sum_em) << " |                             || ";
 	// OUT << setw(7) << Form("%6.3f", ntt_cm_sum_ee) << " |                   || ";
 	// OUT << endl;
 	// OUT << setw(16) << "Pred. Rare SM"  << " || ";
 	// OUT << setw(7) << Form("%6.3f", ntt_rare_mm) << " |                   || ";
 	// OUT << setw(7) << Form("%6.3f", ntt_rare_em) << " |                             || ";
 	// OUT << setw(7) << Form("%6.3f", ntt_rare_ee) << " |                   || ";
 	// OUT << endl;
 	// OUT << "--------------------------------------------------------------------------------------------------------------------------" << endl;
 	// OUT << setw(16) << "Total BG Pred."  << " || ";
 	// OUT << setw(7) << Form("%6.3f", tot_bg_mm) << " | P/O  | P-O/P |    || ";
 	// OUT << setw(7) << Form("%6.3f", tot_bg_em) << " | P/O  | P-O/P |              || ";
 	// OUT << setw(7) << Form("%6.3f", tot_bg_ee) << " | P/O  | P-O/P |    || ";
 	// OUT << endl;
 	// OUT << "--------------------------------------------------------------------------------------------------------------------------" << endl;
 	// OUT << setw(16) << "Observed"  << " || ";
 	// OUT << setw(7) << Form(" %6.3f | %4.2f | %5.2f |", ntt_sum_mm, tot_bg_mm/ntt_sum_mm, (tot_bg_mm-ntt_sum_mm)/tot_bg_mm) << "    || ";
 	// OUT << setw(7) << Form(" %6.3f | %4.2f | %5.2f |", ntt_sum_em, tot_bg_em/ntt_sum_em, (tot_bg_em-ntt_sum_em)/tot_bg_em) << "              || ";
 	// OUT << setw(7) << Form(" %6.3f | %4.2f | %5.2f |", ntt_sum_ee, tot_bg_ee/ntt_sum_ee, (tot_bg_ee-ntt_sum_ee)/tot_bg_ee) << "    || ";
 	// OUT << endl;
 	// OUT << "==========================================================================================================================" << endl;
 	// OUT << Form("     Comb. Fakes || %6.3f |              ||",
 	// npf_pred_sum_mm+nff_pred_sum_mm +
 	// npf_pred_sum_em+nfp_pred_sum_em+nff_pred_sum_em +
 	// npf_pred_sum_ee+nff_pred_sum_ee) << endl;
 	// OUT << Form("    Comb. Ch-MID || %6.3f |              ||", ntt_cm_sum_em+ntt_cm_sum_ee) << endl;
 	// OUT << Form("   Comb. Rare SM || %6.3f |              ||", ntt_rare_mm+ntt_rare_em+ntt_rare_ee) << endl;
 	// OUT << "---------------------------------------------" << endl;
 	// OUT << Form("  Total BG Comb. || %6.3f | P/O  | P-O/P ||", tot_bg) << endl;
 	// OUT << "---------------------------------------------" << endl;
 	// OUT << Form("        Observed || %6.3f | %4.2f | %5.2f ||", ntt_sum, tot_bg/ntt_sum, (tot_bg-ntt_sum)/tot_bg) << endl;
 	// OUT << "=============================================" << endl;
 	// OUT << endl;
	
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
 	//
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
 	// OUT << "  Total fakes:   " <<  Form("%6.3f", npf_pred_sum_ee+nff_pred_sum_ee) << " ±   " << Form("%6.3f", sqrt(nF_pred_sum_ee_e1)) << " (stat)" << endl;
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

	vector<int> ttjets;  ttjets.push_back(TTJets);
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
	Init(tree);

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
	float scale = fLumiNorm / fSamples[TTJets]->getLumi();
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
	OUT << setw(7)  << setprecision(2) << mufratio_allmc << " ± " << setw(7) << setprecision(2) << mufratio_allmc_e << " |";
	OUT << setw(7)  << setprecision(2) << mupratio_allmc << " ± " << setw(7) << setprecision(2) << mupratio_allmc_e << " ||";
	OUT << setw(7)  << setprecision(2) << elfratio_allmc << " ± " << setw(7) << setprecision(2) << elfratio_allmc_e << " |";
	OUT << setw(7)  << setprecision(2) << elpratio_allmc << " ± " << setw(7) << setprecision(2) << elpratio_allmc_e << " ||";
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

void SSDLPlotter::storeWeightedPred(int baseRegion){
	TFile *pFile = TFile::Open(fOutputFileName);
	TTree *sigtree; getObjectSafe(pFile, "SigEvents", sigtree);

	string *sname = 0;
	int flag;
	int   stype, flav, cat, njets, nbjets, nbjetsmed;
	int   passZVeto, pass3rdVeto;
	float puweight, slumi, pT1, pT2, HT, MET, MT2;
	float eta1, eta2;
	int event, run;
	signed int charge;

	
	sigtree->SetBranchAddress("SystFlag",    &flag);
	sigtree->SetBranchAddress("Event",       &event);
	sigtree->SetBranchAddress("Run",         &run);
	sigtree->SetBranchAddress("SName",       &sname);
	sigtree->SetBranchAddress("SType",       &stype);
	sigtree->SetBranchAddress("SLumi",       &slumi);
	sigtree->SetBranchAddress("PUWeight",    &puweight);
	sigtree->SetBranchAddress("Flavor",      &flav);
	sigtree->SetBranchAddress("pT1",         &pT1);
	sigtree->SetBranchAddress("pT2",         &pT2);
	sigtree->SetBranchAddress("eta1",        &eta1);
	sigtree->SetBranchAddress("eta2",        &eta2);
	sigtree->SetBranchAddress("TLCat",       &cat);
	sigtree->SetBranchAddress("HT",          &HT);
	sigtree->SetBranchAddress("MET",         &MET);
	sigtree->SetBranchAddress("MT2",         &MT2);
	sigtree->SetBranchAddress("NJ",          &njets);
	sigtree->SetBranchAddress("NbJ",         &nbjets);
	sigtree->SetBranchAddress("NbJmed",      &nbjetsmed);
	sigtree->SetBranchAddress("Charge",      &charge);
	sigtree->SetBranchAddress("Pass3rdVeto", &pass3rdVeto);
	FakeRatios *FR = new FakeRatios();

	float npp(0.), npf(0.), nfp(0.), nff(0.);
	float f1(0.), f2(0.), p1(0.), p2(0.);
	
	TFile* file_opt;
	TTree* tree_opt;
	float eventWeight;

	if( fDO_OPT ) {
   		system( "mkdir -p OPT_ttW" );
		std::string outputdir_str(fOutputDir);
		outputdir_str.erase(outputdir_str.end()-1); //get rid of last slash
		TString file_opt_name = "OPT_ttW/opt_ttW_" + outputdir_str + ".root";
		file_opt = TFile::Open(file_opt_name, "recreate");
		tree_opt = new TTree("tree_opt", "");

		tree_opt->Branch( "SName",       &*sname );
		tree_opt->Branch( "Flavor",      &flav,        "&flav/I" );
		tree_opt->Branch( "eventWeight", &eventWeight, "eventWeight/F" );
		tree_opt->Branch( "NJ",          &njets,       "njets/I" );
		tree_opt->Branch( "NbJ",         &nbjets,      "nbjets/I" );
		tree_opt->Branch( "NbJmed",      &nbjetsmed,   "nbjetsmed/I" );
		tree_opt->Branch( "pT1",         &pT1,         "pT1/F" );
		tree_opt->Branch( "pT2",         &pT2,         "pT2/F" );
		tree_opt->Branch( "HT",          &HT,          "HT/F" );
		tree_opt->Branch( "MET",         &MET,         "MET/F" );
	}

	
	for( int i = 0; i < sigtree->GetEntries(); i++ ){
		sigtree->GetEntry(i);
		if( flav > 2 ) continue; // OS events
		Sample *S = fSampleMap[TString(*sname)];

		int datamc = S->datamc;
		
		gChannel chan = gChannel(flav);
		f1 = getFRatio(chan, pT1, eta1, S->datamc);
		f2 = getFRatio(chan, pT2, eta2, S->datamc);
		p1 = getPRatio(chan, pT1, S->datamc);
		p2 = getPRatio(chan, pT2, S->datamc);
		if(chan == ElMu){
			f1 = getFRatio(Muon, pT1, eta1, S->datamc);
			f2 = getFRatio(Elec, pT2, eta2, S->datamc);
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
		for(int r = 0; r < gNREGIONS; r++){
			// Select correct incarnation for each region
// 			if(r  < TTbarWSelJU && flag != 0) continue;
// 			if(r == TTbarWSelJU && flag != 1) continue;
// 			if(r == TTbarWSelJD && flag != 2) continue;
// 			if(r == TTbarWSelJS && flag != 3) continue;
// 			if(r == TTbarWSelBU && flag != 4) continue;
// 			if(r == TTbarWSelBD && flag != 5) continue;
// 			if(r == TTbarWSelLU && flag != 6) continue;
// 			if(r == TTbarWSelLD && flag != 7) continue;
			
			// Event Selection:
			if(HT        < gRegions[r]->minHT     || HT        > gRegions[r]->maxHT    )   continue;
			if(MET       < gRegions[r]->minMet    || MET       > gRegions[r]->maxMet   )   continue;
			if(njets     < gRegions[r]->minNjets  || njets     > gRegions[r]->maxNjets )   continue;
			if(nbjets    < gRegions[r]->minNbjets || nbjets    > gRegions[r]->maxNbjets)   continue;
			if(nbjetsmed < gRegions[r]->minNbjmed || nbjetsmed > gRegions[r]->maxNbjmed)   continue;
			if(gRegions[r]->app3rdVet  != 0       && pass3rdVeto == 1)                     continue;
			if(gRegions[r]->chargeVeto != 0       && (charge == gRegions[r]->chargeVeto) ) continue;
			   
			if(passesPtCuts(pT1, pT2, r, chan) == false) continue;
			
			S->numbers[r][chan].npp += puweight * npp;
			S->numbers[r][chan].npf += puweight * npf;
			S->numbers[r][chan].nfp += puweight * nfp;
			S->numbers[r][chan].nff += puweight * nff;
		}
		
		// Differential predictions
		if(flag != 0) continue;
		float maxpt = TMath::Max(pT1, pT2);
 		float minpt = TMath::Min(pT1, pT2);
		
 		if(MET > 30. && maxpt > 20. && minpt > 20.){
		  // not in benjamins code if(HT        >= gRegions[baseRegion]->minHT     &&
		  // not in benjamins code    HT        <  gRegions[baseRegion]->maxHT     &&
		  // not in benjamins code    njets     >= gRegions[baseRegion]->minNjets  &&
		  // not in benjamins code    njets     <= gRegions[baseRegion]->maxNjets  && 
		  // not in benjamins code    nbjets    >= gRegions[baseRegion]->minNbjets &&
		  // not in benjamins code    nbjets    <= gRegions[baseRegion]->maxNbjets &&
		  // not in benjamins code    nbjetsmed >= gRegions[baseRegion]->minNbjmed &&
		  // not in benjamins code    nbjetsmed <= gRegions[baseRegion]->maxNbjmed){
			fillWithoutOF(S->diffyields[chan].hnpp[7], MET, puweight * npp);
			fillWithoutOF(S->diffyields[chan].hnpf[7], MET, puweight * npf);
			fillWithoutOF(S->diffyields[chan].hnfp[7], MET, puweight * nfp);
			fillWithoutOF(S->diffyields[chan].hnff[7], MET, puweight * nff);
		  // not in benjamins code    }
		}

		// Check pt cuts of baseRegion:
 		bool passespt = passesPtCuts(pT1, pT2, baseRegion, chan);
		
		if(HT        >= gRegions[baseRegion]->minHT     &&
		   HT        <  gRegions[baseRegion]->maxHT     &&
 		   MET       >  gRegions[baseRegion]->minMet    &&
 		   MET       <  gRegions[baseRegion]->maxMet    &&
 		   nbjets    >= gRegions[baseRegion]->minNbjets &&
 		   nbjets    <= gRegions[baseRegion]->maxNbjets &&
 		   nbjetsmed >= gRegions[baseRegion]->minNbjmed &&
		   nbjetsmed <= gRegions[baseRegion]->maxNbjmed &&
 		   passespt) 		  
		{
			fillWithoutOF(S->diffyields[chan].hnpp[2], njets+0.5, puweight * npp);
			fillWithoutOF(S->diffyields[chan].hnpf[2], njets+0.5, puweight * npf);
			fillWithoutOF(S->diffyields[chan].hnfp[2], njets+0.5, puweight * nfp);
 			fillWithoutOF(S->diffyields[chan].hnff[2], njets+0.5, puweight * nff);
  	
 		if(njets     >= gRegions[baseRegion]->minNjets  &&
		   njets     <= gRegions[baseRegion]->maxNjets  )
		   {

			fillWithoutOF(S->diffyields[chan].hnpp[0], HT, puweight * npp);
			fillWithoutOF(S->diffyields[chan].hnpf[0], HT, puweight * npf);
			fillWithoutOF(S->diffyields[chan].hnfp[0], HT, puweight * nfp);
			fillWithoutOF(S->diffyields[chan].hnff[0], HT, puweight * nff);
			
			fillWithoutOF(S->diffyields[chan].hnpp[1], MET, puweight * npp);
			fillWithoutOF(S->diffyields[chan].hnpf[1], MET, puweight * npf);
			fillWithoutOF(S->diffyields[chan].hnfp[1], MET, puweight * nfp);
			fillWithoutOF(S->diffyields[chan].hnff[1], MET, puweight * nff);
			
			fillWithoutOF(S->diffyields[chan].hnpp[3], MT2, puweight * npp);
			fillWithoutOF(S->diffyields[chan].hnpf[3], MT2, puweight * npf);
			fillWithoutOF(S->diffyields[chan].hnfp[3], MT2, puweight * nfp);
			fillWithoutOF(S->diffyields[chan].hnff[3], MT2, puweight * nff);
			
			fillWithoutOF(S->diffyields[chan].hnpp[4], std::max(pT1,pT2), puweight * npp);
			fillWithoutOF(S->diffyields[chan].hnpf[4], std::max(pT1,pT2), puweight * npf);
			fillWithoutOF(S->diffyields[chan].hnfp[4], std::max(pT1,pT2), puweight * nfp);
			fillWithoutOF(S->diffyields[chan].hnff[4], std::max(pT1,pT2), puweight * nff);
			
			fillWithoutOF(S->diffyields[chan].hnpp[5], std::min(pT1,pT2), puweight * npp);
			fillWithoutOF(S->diffyields[chan].hnpf[5], std::min(pT1,pT2), puweight * npf);
			fillWithoutOF(S->diffyields[chan].hnfp[5], std::min(pT1,pT2), puweight * nfp);
			fillWithoutOF(S->diffyields[chan].hnff[5], std::min(pT1,pT2), puweight * nff);
			
			fillWithoutOF(S->diffyields[chan].hnpp[6], nbjets+0.5, puweight * npp);
			fillWithoutOF(S->diffyields[chan].hnpf[6], nbjets+0.5, puweight * npf);
			fillWithoutOF(S->diffyields[chan].hnfp[6], nbjets+0.5, puweight * nfp);
			fillWithoutOF(S->diffyields[chan].hnff[6], nbjets+0.5, puweight * nff);
			
			fillWithoutOF(S->diffyields[chan].hnpp[8], nbjetsmed+0.5, puweight * npp);
			fillWithoutOF(S->diffyields[chan].hnpf[8], nbjetsmed+0.5, puweight * npf);
			fillWithoutOF(S->diffyields[chan].hnfp[8], nbjetsmed+0.5, puweight * nfp);
			fillWithoutOF(S->diffyields[chan].hnff[8], nbjetsmed+0.5, puweight * nff);
			
			fillWithoutOF(S->diffyields[chan].hnpp[9], nbjetsmed+0.5, puweight * npp);
			fillWithoutOF(S->diffyields[chan].hnpf[9], nbjetsmed+0.5, puweight * npf);
			fillWithoutOF(S->diffyields[chan].hnfp[9], nbjetsmed+0.5, puweight * nfp);
			fillWithoutOF(S->diffyields[chan].hnff[9], nbjetsmed+0.5, puweight * nff);
		   }
		}

		if( fDO_OPT && flav<3 && flag == 0){
			float lumi_pb = 9000.;
			eventWeight = lumi_pb*puweight/slumi;
			tree_opt->Fill();
		}
	}

	if( fDO_OPT ){
		file_opt->cd();
		tree_opt->Write();
		file_opt->Close();
	}
	delete FR;
}

float SSDLPlotter::getFRatio(gChannel chan, float pt, int datamc){
	// if(chan == Muon) return 0.0672; // flat ratios
	// if(chan == Elec) return 0.224;

	const float mu_flatout = 30.;
	const float el_flatout = 30.;
	
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
float SSDLPlotter::getFRatio(gChannel chan, float pt, float eta, int datamc){
	// if(chan == Muon) return 0.0672; // flat ratios
	// if(chan == Elec) return 0.224;
	const float mu_flatout = 30.;
	const float el_flatout = 30.;
	
	eta = fabs(eta); // make sure we take the absolute value

	if(chan == Muon){
		TH2D *histo          = fH2D_MufRatio;
		if(datamc > 0) histo = fH2D_MufRatio_MC;
		if(!histo){
			cerr << "SSDLPlotter::getPRatio ==> Warning: ratio histo not filled, exiting" << endl;
			exit(-1);
		}
		if(pt >= mu_flatout){
			Int_t asym_val_bin = histo->FindBin(mu_flatout-1., eta);
			return histo->GetBinContent(asym_val_bin); // findbin(35) returns bin from 35-45, want the previous one
		}
		Int_t binnumb = histo->FindBin(pt, eta);
		return histo->GetBinContent(binnumb);
	}
	if(chan == Elec){
		TH2D *histo          = fH2D_ElfRatio;
		if(datamc > 0) histo = fH2D_ElfRatio_MC;
		if(!histo){
			cerr << "SSDLPlotter::getPRatio ==> Warning: ratio histo not filled, exiting" << endl;
			exit(-1);
		}
		if(pt >= el_flatout){
			Int_t asym_val_bin = histo->FindBin(el_flatout-1., eta);
			return histo->GetBinContent(asym_val_bin);
		}
		Int_t binnumb = histo->FindBin(pt, eta);
		return histo->GetBinContent(binnumb);
	}
}
float SSDLPlotter::getPRatio(gChannel chan, float pt, int datamc){
	// if(chan == Muon) return 0.941; // flat ratios
	// if(chan == Elec) return 0.922;
	
//	return 1.;
	
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
		NumberSet numbers = S->numbers[gRegion[gBaseRegion]][chan];
		float scale = lumiscale / S->getLumi();
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
		nt20[Muon]    += S->numbers[gRegion[gBaseRegion]][Muon].nt2;
		nt10[Muon]    += S->numbers[gRegion[gBaseRegion]][Muon].nt10;
		nt00[Muon]    += S->numbers[gRegion[gBaseRegion]][Muon].nt0;
	}

	for(size_t i = 0; i < emusamples.size(); ++i){
		Sample *S = fSamples[emusamples[i]];
		nt20[ElMu]    += S->numbers[gRegion[gBaseRegion]][ElMu].nt2;
		nt10[ElMu]    += S->numbers[gRegion[gBaseRegion]][ElMu].nt10;
		nt01[ElMu]    += S->numbers[gRegion[gBaseRegion]][ElMu].nt01;
		nt00[ElMu]    += S->numbers[gRegion[gBaseRegion]][ElMu].nt0;
	}		

	for(size_t i = 0; i < elsamples.size(); ++i){
		Sample *S = fSamples[elsamples[i]];
		nt20[Elec]    += S->numbers[gRegion[gBaseRegion]][Elec].nt2;
		nt10[Elec]    += S->numbers[gRegion[gBaseRegion]][Elec].nt10;
		nt00[Elec]    += S->numbers[gRegion[gBaseRegion]][Elec].nt0;
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
		float scale = luminorm / S->getLumi();
		if(luminorm < 0) scale = 1;
		nt2sum_mumu  += scale*S->numbers[gRegion[gBaseRegion]][Muon]    .nt2;
		nt10sum_mumu += scale*S->numbers[gRegion[gBaseRegion]][Muon]    .nt10;
		nt0sum_mumu  += scale*S->numbers[gRegion[gBaseRegion]][Muon]    .nt0;
		nt2sum_emu   += scale*S->numbers[gRegion[gBaseRegion]][ElMu]     .nt2;
		nt10sum_emu  += scale*S->numbers[gRegion[gBaseRegion]][ElMu]     .nt10;
		nt01sum_emu  += scale*S->numbers[gRegion[gBaseRegion]][ElMu]     .nt01;
		nt0sum_emu   += scale*S->numbers[gRegion[gBaseRegion]][ElMu]     .nt0;
		nt2sum_ee    += scale*S->numbers[gRegion[gBaseRegion]][Elec].nt2;
		nt10sum_ee   += scale*S->numbers[gRegion[gBaseRegion]][Elec].nt10;
		nt0sum_ee    += scale*S->numbers[gRegion[gBaseRegion]][Elec].nt0;

		cout << setw(9) << S->sname << " || ";
		cout << setw(7) << scale*S->numbers[gRegion[gBaseRegion]][Muon]    .nt2  << " | ";
		cout << setw(7) << scale*S->numbers[gRegion[gBaseRegion]][Muon]    .nt10 << " | ";
		cout << setw(7) << scale*S->numbers[gRegion[gBaseRegion]][Muon]    .nt0  << " || ";
		cout << setw(7) << scale*S->numbers[gRegion[gBaseRegion]][ElMu]     .nt2  << " | ";
		cout << setw(7) << scale*S->numbers[gRegion[gBaseRegion]][ElMu]     .nt10 << " | ";
		cout << setw(7) << scale*S->numbers[gRegion[gBaseRegion]][ElMu]     .nt01 << " | ";
		cout << setw(7) << scale*S->numbers[gRegion[gBaseRegion]][ElMu]     .nt0  << " || ";
		cout << setw(7) << scale*S->numbers[gRegion[gBaseRegion]][Elec].nt2  << " | ";
		cout << setw(7) << scale*S->numbers[gRegion[gBaseRegion]][Elec].nt10 << " | ";
		cout << setw(7) << scale*S->numbers[gRegion[gBaseRegion]][Elec].nt0  << " || ";
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
	// MARC cout << setw(9) << fSamples[LM0]->sname << " || ";
	// MARC float scale = luminorm / fSamples[LM0]->getLumi();
	// MARC if(luminorm < 0) scale = 1;
	// MARC cout << setw(7) << scale * fSamples[LM0]->numbers[Baseline][Muon]    .nt2  << " | ";
	// MARC cout << setw(7) << scale * fSamples[LM0]->numbers[Baseline][Muon]    .nt10 << " | ";
	// MARC cout << setw(7) << scale * fSamples[LM0]->numbers[Baseline][Muon]    .nt0  << " || ";
	// MARC cout << setw(7) << scale * fSamples[LM0]->numbers[Baseline][ElMu]     .nt2  << " | ";
	// MARC cout << setw(7) << scale * fSamples[LM0]->numbers[Baseline][ElMu]     .nt10 << " | ";
	// MARC cout << setw(7) << scale * fSamples[LM0]->numbers[Baseline][ElMu]     .nt01 << " | ";
	// MARC cout << setw(7) << scale * fSamples[LM0]->numbers[Baseline][ElMu]     .nt0  << " || ";
	// MARC cout << setw(7) << scale * fSamples[LM0]->numbers[Baseline][Elec].nt2  << " | ";
	// MARC cout << setw(7) << scale * fSamples[LM0]->numbers[Baseline][Elec].nt10 << " | ";
	// MARC cout << setw(7) << scale * fSamples[LM0]->numbers[Baseline][Elec].nt0  << " || ";
	// MARC cout << endl;
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
void SSDLPlotter::printAllYieldTables(){
	TString outputdir = Util::MakeOutputDir(fOutputDir + "YieldTables");
	for(size_t i = 0; i < gNREGIONS; ++i){
		TString outputname = outputdir + "YieldTable_" + gRegions[i]->sname + ".txt";
		printMCYieldTable(outputname, i);
	}
	
}
void SSDLPlotter::printMCYieldTable(TString filename, int reg){
	ofstream OUT(filename.Data(), ios::trunc);
	OUT << "/////////////////////////////////////////////////////////////////////////////" << endl;
	OUT << " Printing yields for region " << gRegions[reg]->sname << endl;
	OUT << "  scaling to " << fLumiNorm << " /pb" << endl << endl;

	vector<int> musamples, elsamples, emusamples, mcsamples;

	musamples  = fMuData;
	elsamples  = fEGData;
	emusamples = fMuEGData;
	mcsamples  = fMCBG;

	int nt20[gNCHANNELS], nt10[gNCHANNELS], nt01[gNCHANNELS], nt00[gNCHANNELS]; // data yields
	const int nprocs = fSamples[DoubleMu1]->getNProcs();
	float ntt_mm[nprocs], ntl_mm[nprocs], nll_mm[nprocs]; // yields per process
	float ntt_em[nprocs], ntl_em[nprocs], nlt_em[nprocs], nll_em[nprocs];
	float ntt_ee[nprocs], ntl_ee[nprocs], nll_ee[nprocs];
	float ntt_mm_e2[nprocs], ntt_em_e2[nprocs], ntt_ee_e2[nprocs]; // squared errors
	float ntt_sum_mm(0.), ntl_sum_mm(0.), nll_sum_mm(0.); // total sum
	float ntt_sum_em(0.), ntl_sum_em(0.), nlt_sum_em(0.), nll_sum_em(0.);
	float ntt_sum_ee(0.), ntl_sum_ee(0.), nll_sum_ee(0.);
	float ntt_sum_mm_e2(0.), ntt_sum_em_e2(0.), ntt_sum_ee_e2(0.);
	float ntt_sum(0.), ntl_sum(0.), nlt_sum(0.), nll_sum(0.);
	
	for(size_t i = 0; i < nprocs; ++i){ // reset everything
		ntt_mm[i] = 0.; ntl_mm[i] = 0.; nll_mm[i] = 0.;
		ntt_em[i] = 0.; ntl_em[i] = 0.; nlt_em[i] = 0.; nll_em[i] = 0.;
		ntt_ee[i] = 0.; ntl_ee[i] = 0.; nll_ee[i] = 0.;
		ntt_mm_e2[i] = 0.; ntt_em_e2[i] = 0.; ntt_ee_e2[i] = 0.; // squared errors
	}
	for(size_t i = 0; i < gNCHANNELS; ++i){ // reset everything
		nt20[i] = 0; nt10[i] = 0; nt01[i] = 0; nt00[i] = 0;
	}

	for(size_t i = 0; i < mcsamples.size(); ++i){
		Sample *S = fSamples[mcsamples[i]];
		float scale = fLumiNorm / S->getLumi();
		int proc = S->getProc();
		if(proc == 0 || proc >= nprocs) continue; // safety
		
		float ntt_mm_temp = gMMTrigScale*scale*S->region[reg][HighPt].mm.nt20_pt->Integral(0, getNFPtBins(Muon)+1);
		float ntl_mm_temp = gMMTrigScale*scale*S->region[reg][HighPt].mm.nt10_pt->Integral(0, getNFPtBins(Muon)+1);
		float nll_mm_temp = gMMTrigScale*scale*S->region[reg][HighPt].mm.nt00_pt->Integral(0, getNFPtBins(Muon)+1);
		float ntt_em_temp = gEMTrigScale*scale*S->region[reg][HighPt].em.nt20_pt->Integral(0, getNFPtBins(ElMu)+1);
		float ntl_em_temp = gEMTrigScale*scale*S->region[reg][HighPt].em.nt10_pt->Integral(0, getNFPtBins(ElMu)+1);
		float nlt_em_temp = gEMTrigScale*scale*S->region[reg][HighPt].em.nt01_pt->Integral(0, getNFPtBins(ElMu)+1);
		float nll_em_temp = gEMTrigScale*scale*S->region[reg][HighPt].em.nt00_pt->Integral(0, getNFPtBins(ElMu)+1);
		float ntt_ee_temp = gEETrigScale*scale*S->region[reg][HighPt].ee.nt20_pt->Integral(0, getNFPtBins(Elec)+1);
		float ntl_ee_temp = gEETrigScale*scale*S->region[reg][HighPt].ee.nt10_pt->Integral(0, getNFPtBins(Elec)+1);
		float nll_ee_temp = gEETrigScale*scale*S->region[reg][HighPt].ee.nt00_pt->Integral(0, getNFPtBins(Elec)+1);

		ntt_mm[proc] += ntt_mm_temp;
		ntl_mm[proc] += ntl_mm_temp;
		nll_mm[proc] += nll_mm_temp;
		ntt_em[proc] += ntt_em_temp;
		ntl_em[proc] += ntl_em_temp;
		nlt_em[proc] += nlt_em_temp;
		nll_em[proc] += nll_em_temp;
		ntt_ee[proc] += ntt_ee_temp;
		ntl_ee[proc] += ntl_ee_temp;
		nll_ee[proc] += nll_ee_temp;

		ntt_sum_mm += ntt_mm_temp;
		ntl_sum_mm += ntl_mm_temp;
		nll_sum_mm += nll_mm_temp;
		ntt_sum_em += ntt_em_temp;
		ntl_sum_em += ntl_em_temp;
		nlt_sum_em += nlt_em_temp;
		nll_sum_em += nll_em_temp;
		ntt_sum_ee += ntt_ee_temp;
		ntl_sum_ee += ntl_ee_temp;
		nll_sum_ee += nll_ee_temp;
		
		if(proc == 5 || proc == 17) continue; // fuck QCD and g+jets for the errors
		// Careful, I assume now that they will always have 0 events in tight-tight,
		// otherwise my yield table won't be self consistent anymore

		// Errors
		float ntt_mm_e2_temp = pow(gMMTrigScale*scale*S->numbers[reg][Muon].tt_avweight*S->getError(S->region[reg][HighPt].mm.nt20_pt->GetEntries()),2);
		float ntt_em_e2_temp = pow(gEMTrigScale*scale*S->numbers[reg][ElMu].tt_avweight*S->getError(S->region[reg][HighPt].em.nt20_pt->GetEntries()),2);
		float ntt_ee_e2_temp = pow(gEETrigScale*scale*S->numbers[reg][Elec].tt_avweight*S->getError(S->region[reg][HighPt].ee.nt20_pt->GetEntries()),2);
		ntt_mm_e2[proc] += ntt_mm_e2_temp;
		ntt_em_e2[proc] += ntt_em_e2_temp;
		ntt_ee_e2[proc] += ntt_ee_e2_temp;
		ntt_sum_mm_e2 += ntt_mm_e2_temp;
		ntt_sum_em_e2 += ntt_em_e2_temp;
		ntt_sum_ee_e2 += ntt_ee_e2_temp;
	}

	for(size_t i = 0; i < musamples.size(); ++i){
		Sample *S = fSamples[musamples[i]];
		nt20[Muon] += S->numbers[reg][Muon].nt2;
		nt10[Muon] += S->numbers[reg][Muon].nt10;
		nt00[Muon] += S->numbers[reg][Muon].nt0;
	}
	for(size_t i = 0; i < emusamples.size(); ++i){
		Sample *S = fSamples[emusamples[i]];
		nt20[ElMu] += S->numbers[reg][ElMu].nt2;
		nt10[ElMu] += S->numbers[reg][ElMu].nt10;
		nt01[ElMu] += S->numbers[reg][ElMu].nt01;
		nt00[ElMu] += S->numbers[reg][ElMu].nt0;
	}		
	for(size_t i = 0; i < elsamples.size(); ++i){
		Sample *S = fSamples[elsamples[i]];
		nt20[Elec] += S->numbers[reg][Elec].nt2;
		nt10[Elec] += S->numbers[reg][Elec].nt10;
		nt00[Elec] += S->numbers[reg][Elec].nt0;
	}

	OUT << "--------------------------------------------------------------------------------------------------------------------" << endl;
	OUT << "                         |          Mu/Mu           |               E/Mu                |           E/E            |" << endl;
	OUT << "  Process                |  Ntt   |  Ntl   |  Nll   |  Ntt   |  Ntl   |  Nlt   |  Nll   |  Ntt   |  Ntl   |  Nll   |" << endl;
	OUT << "--------------------------------------------------------------------------------------------------------------------" << endl;

	for(size_t i = 1; i < nprocs; ++i){ // skip data
		OUT << Form("%24s & %6.2f & %6.2f & %6.2f & %6.2f & %6.2f & %6.2f & %6.2f & %6.2f & %6.2f & %6.2f \\\\",
		fSamples[DoubleMu1]->getProcName(i).Data(),
		ntt_mm[i], ntl_mm[i], nll_mm[i],
		ntt_em[i], ntl_em[i], nlt_em[i], nll_em[i],
		ntt_ee[i], ntl_ee[i], nll_ee[i]) << endl;
	}
	OUT << "\\hline" << endl;
	OUT << Form("%24s & %6.2f & %6.2f & %6.2f & %6.2f & %6.2f & %6.2f & %6.2f & %6.2f & %6.2f & %6.2f \\\\",
	"Sum", ntt_sum_mm, ntl_sum_mm, nll_sum_mm, ntt_sum_em, ntl_sum_em, nlt_sum_em, nll_sum_em, ntt_sum_ee, ntl_sum_ee, nll_sum_ee) << endl;
	OUT << "\\hline" << endl;
	for(gSample i = sample_begin; i < gNSAMPLES; i=gSample(i+1)){
		Sample *S = fSamples[i];
		if(S->datamc != 2) continue;
		float scale = fLumiNorm / S->getLumi();

		float temp_nt2_mm  = gMMTrigScale*scale*S->region[reg][HighPt].mm.nt20_pt->Integral(0, getNFPtBins(Muon)+1);
		float temp_nt10_mm = gMMTrigScale*scale*S->region[reg][HighPt].mm.nt10_pt->Integral(0, getNFPtBins(Muon)+1);
		float temp_nt0_mm  = gMMTrigScale*scale*S->region[reg][HighPt].mm.nt00_pt->Integral(0, getNFPtBins(Muon)+1);
		float temp_nt2_em  = gEMTrigScale*scale*S->region[reg][HighPt].em.nt20_pt->Integral(0, getNFPtBins(ElMu)+1);
		float temp_nt10_em = gEMTrigScale*scale*S->region[reg][HighPt].em.nt10_pt->Integral(0, getNFPtBins(ElMu)+1);
		float temp_nt01_em = gEMTrigScale*scale*S->region[reg][HighPt].em.nt01_pt->Integral(0, getNFPtBins(ElMu)+1);
		float temp_nt0_em  = gEMTrigScale*scale*S->region[reg][HighPt].em.nt00_pt->Integral(0, getNFPtBins(ElMu)+1);
		float temp_nt2_ee  = gEETrigScale*scale*S->region[reg][HighPt].ee.nt20_pt->Integral(0, getNFPtBins(Elec)+1);
		float temp_nt10_ee = gEETrigScale*scale*S->region[reg][HighPt].ee.nt10_pt->Integral(0, getNFPtBins(Elec)+1);
		float temp_nt0_ee  = gEETrigScale*scale*S->region[reg][HighPt].ee.nt00_pt->Integral(0, getNFPtBins(Elec)+1);

		TString tempname = S->sname;
		OUT << Form("%24s & %6.2f & %6.2f & %6.2f & %6.2f & %6.2f & %6.2f & %6.2f & %6.2f & %6.2f & %6.2f \\\\\n", (tempname.ReplaceAll("_","\\_")).Data(),
		temp_nt2_mm , temp_nt10_mm, temp_nt0_mm,
		temp_nt2_em , temp_nt10_em, temp_nt01_em, temp_nt0_em,
		temp_nt2_ee , temp_nt10_ee, temp_nt0_ee);
	}	

	OUT << "\\hline" << endl;
	OUT << Form("%24s & %6d & %6d & %6d & %6d & %6d & %6d & %6d & %6d & %6d & %6d \\\\",
	"Data", nt20[Muon], nt10[Muon], nt00[Muon], nt20[ElMu], nt10[ElMu], nt01[ElMu], nt00[ElMu], nt20[Elec], nt10[Elec], nt00[Elec]) << endl;
	OUT << "--------------------------------------------------------------------------------------------------------------------" << endl;

	OUT << endl;
	OUT << "REMINDER: CHECK THAT QCD AND G+JETS REALLY HAVE ZERO YIELDS IN TIGHT-TIGHT EVERYWHERE" << endl;
	OUT << "          OTHERWISE THE SUMS IN THE TABLE BELOW ARE NOT CONSISTENT!" << endl;
	OUT << endl;


	OUT << "------------------------------------------------------------------------------------------------------------------" << endl;
	OUT << "  Process                |  Mu/Mu              |  E/Mu               |  E/E                |  Tot                |" << endl;
	OUT << "------------------------------------------------------------------------------------------------------------------" << endl;
	
	for(size_t i = 1; i < nprocs; ++i){ // skip data
		if(i == 5 || i == 16) continue; // fuck QCD and g+jets
		OUT << Form("%24s & %6.2f $\\pm$ %6.2f & %6.2f $\\pm$ %6.2f & %6.2f $\\pm$ %6.2f & %6.2f $\\pm$ %6.2f \\\\",
		fSamples[DoubleMu1]->getProcName(i).Data(),
		ntt_mm[i], sqrt(ntt_mm_e2[i]), ntt_em[i], sqrt(ntt_em_e2[i]), ntt_ee[i], sqrt(ntt_ee_e2[i]),
		ntt_mm[i]+ntt_em[i]+ntt_ee[i], sqrt(ntt_mm_e2[i]+ntt_em_e2[i]+ntt_ee_e2[i])) << endl;
	}
	OUT << "\\hline \\hline" << endl;
	OUT << Form("%24s & %6.2f $\\pm$ %6.2f & %6.2f $\\pm$ %6.2f & %6.2f $\\pm$ %6.2f & %6.2f $\\pm$ %6.2f \\\\",
	"Sum", ntt_sum_mm, sqrt(ntt_sum_mm_e2), ntt_sum_em, sqrt(ntt_sum_em_e2), ntt_sum_ee, sqrt(ntt_sum_ee_e2),
	       ntt_sum_mm+ntt_sum_em+ntt_sum_ee, sqrt(ntt_sum_mm_e2+ntt_sum_em_e2+ntt_sum_ee_e2)) << endl;
	
	OUT << "\\hline \\hline" << endl;
	OUT << Form("%24s & %6d          & %6d          & %6d          & %6d          \\\\",
	"Data", nt20[Muon], nt20[ElMu], nt20[Elec], nt20[Muon]+nt20[ElMu]+nt20[Elec]) << endl;
OUT << "------------------------------------------------------------------------------------------------------------------" << endl;
	OUT << endl;
	OUT.close();
}

//____________________________________________________________________________
TGraph* SSDLPlotter::getSigEventGraph(gChannel chan, int reg){
	TString channame = "MM";
	if(chan == Elec) channame = "EE";
	if(chan == ElMu)      channame = "EM";
	vector<float> ht;
	vector<float> met;

	TFile *pFile = TFile::Open(fOutputFileName);
	TTree *sigtree; getObjectSafe(pFile, "SigEvents", sigtree);

	string *sname = 0;
	int flag;
	int   stype, flav, cat, njets, nbjets;
	float puweight, pT1, pT2, HT, MET, MT2;

	sigtree->SetBranchAddress("SystFlag",    &flag);
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
	
	for( int i = 0; i < sigtree->GetEntries(); i++ ){
		sigtree->GetEntry(i);
		Sample *S = fSampleMap[TString(*sname)];
		int datamc = S->datamc;
		
		if(flag != 0)              continue; // Only choose nominal
		if(stype > 2)              continue; // 0,1,2 are DoubleMu, DoubleEle, MuEG
		if(cat != 0)               continue; // tight-tight selection
		if(gChannel(flav) != chan) continue; // channel selection
		
		// Region selections
		if(HT     < gRegions[reg]->minHT     || HT  >    gRegions[reg]->maxHT ) continue;
		if(MET    < gRegions[reg]->minMet    || MET >    gRegions[reg]->maxMet) continue;
		if(njets  < gRegions[reg]->minNjets  || njets  > gRegions[reg]->maxNjets ) continue;
		if(nbjets < gRegions[reg]->minNbjets || nbjets > gRegions[reg]->maxNbjets) continue;
		
		ht.push_back(HT);
		met.push_back(MET);
	}
	
	const int nsig = ht.size();
	float ht_a [nsig];
	float met_a[nsig];
	for(size_t i = 0; i < ht.size(); ++i){
		ht_a[i] = ht[i];
		met_a[i] = met[i];
	}
	
	Color_t color[3] = {kBlack, kBlue, kRed};
	Size_t size = 1.5;
	Style_t style[3] = {8, 23, 21};
	
	TGraph *sigevents = new TGraph(nsig, ht_a, met_a);
	sigevents->SetName(Form("%s_%s_SigEvents", gRegions[reg]->sname.Data(), channame.Data()));
	
	sigevents->SetMarkerColor(color[chan]);
	sigevents->SetMarkerStyle(style[chan]);
	sigevents->SetMarkerSize(size);
	
	return sigevents;
}
TGraph* SSDLPlotter::getSigEventGraph(gChannel chan, float HTmin, float HTmax, float METmin, float METmax, int which){
	TString channame = "MM";
	if(chan == Elec) channame = "EE";
	if(chan == ElMu) channame = "EM";
	vector<float> ht;
	vector<float> met;
	vector<float> nj;
	vector<int> nb;

	TFile *pFile = TFile::Open(fOutputFileName);
	TTree *sigtree; getObjectSafe(pFile, "SigEvents", sigtree);

	string *sname = 0;
	int flag;
	int   stype, flav, cat;
	float HT, MET;
	int njets, nbjets;

	sigtree->SetBranchAddress("SystFlag", &flag);
	sigtree->SetBranchAddress("SName",    &sname);
	sigtree->SetBranchAddress("SType",    &stype);
	sigtree->SetBranchAddress("Flavor",   &flav);
	sigtree->SetBranchAddress("TLCat",    &cat);
	sigtree->SetBranchAddress("HT",       &HT);
	sigtree->SetBranchAddress("MET",      &MET);
	sigtree->SetBranchAddress("NbJ",      &nbjets);
	sigtree->SetBranchAddress("NJ",       &njets);
	
	for( int i = 0; i < sigtree->GetEntries(); i++ ){
		sigtree->GetEntry(i);
		Sample *S = fSampleMap[TString(*sname)];
		int datamc = S->datamc;
		
		if(flag != 0)              continue; // Only choose nominal
		if(stype > 2)              continue; // 0,1,2 are DoubleMu, DoubleEle, MuEG
		if(cat != 0)               continue; // tight-tight selection
		if(gChannel(flav) != chan) continue; // channel selection
		
		// Region selections
		if(HT  < HTmin  || HT  > HTmax ) continue;
		if(MET < METmin || MET > METmax) continue;
		
		// ht.push_back(HT);
		// met.push_back(MET);
		if ( which == 0 ) { nj.push_back(njets); ht.push_back(HT); met.push_back(MET);}
		else if(which == 1) {
			if (njets > 0 && nbjets >  0) {
				if      (njets <  5) {nj.push_back(njets + 0.25); ht.push_back(HT); met.push_back(MET);}
				else if (njets >= 5) {nj.push_back(5.25); ht.push_back(HT); met.push_back(MET);}
			}
			else if (njets > 0 && nbjets == 0) {
				if      (njets <  5) {nj.push_back(njets - 0.25); ht.push_back(HT); met.push_back(MET);}
				else if (njets >= 5) {nj.push_back(4.75); ht.push_back(HT); met.push_back(MET);}
			}
			else if (njets == 0) {nj.push_back(njets); ht.push_back(HT); met.push_back(MET);}
		}
		else if ( which == 2 && nbjets > 0){
			if      (njets <= 5) {nj.push_back(njets); ht.push_back(HT); met.push_back(MET);}
			else if (njets >  5) {nj.push_back(5); ht.push_back(HT); met.push_back(MET);}
		}
		else if ( which == 3 && nbjets == 0){
			if      (njets <= 5) {nj.push_back(njets); ht.push_back(HT); met.push_back(MET);}
			else if (njets >  5) {nj.push_back(5); ht.push_back(HT); met.push_back(MET);}
		}
		else if ( which == 4){
			if      (njets <= 5) {nj.push_back(njets); ht.push_back(HT); met.push_back(MET);}
			else if (njets >  5) {nj.push_back(5); ht.push_back(HT); met.push_back(MET);}
		}
		// nj.push_back(njets);
		// nb.push_back(nbjets);
	}
	
	const int nsig  = ht.size();
	float ht_a [nsig];
	float met_a[nsig];
	float nj_a[nj.size()];
	for(size_t i = 0; i < ht.size(); ++i){
		ht_a[i] = ht[i];
		met_a[i] = met[i];
	}
	for(size_t i = 0; i < nj.size(); ++i){
		nj_a[i] = nj[i];
	}
	
	Color_t color[3] = {kBlack, kBlue, kRed};
	Size_t size = 1.5;
	Style_t style[3] = {8, 23, 21};
	
	TGraph *sigevents;// = new TGraph(nsig, ht_a, met_a);
	if (which == 0) {
		sigevents = new TGraph(nsig, ht_a, met_a);
		sigevents->SetName(Form("HT%4.0f-%4.0f_MET%4.0f-%4.0f_%s_SigEvents", HTmin, HTmax, METmin, METmax, channame.Data()));
		
		sigevents->SetMarkerColor(color[chan]);
		sigevents->SetMarkerStyle(style[chan]);
		sigevents->SetMarkerSize(size);
	}
	if (which > 0) {
		sigevents = new TGraph(nj.size(), nj_a, met_a);
		sigevents->SetName(Form("HT%4.0f-%4.0f_MET%4.0f-%4.0f_%s_SigEvents", HTmin, HTmax, METmin, METmax, channame.Data()));
		
		sigevents->SetMarkerColor(color[chan]);
		sigevents->SetMarkerStyle(style[chan]);
		sigevents->SetMarkerSize(size);
	}
	
	return sigevents;
}

//////////////////////////////////////////////////////////////////////////////
// Geninfo stuff
//____________________________________________________________________________
void SSDLPlotter::makeOriginPlots(int reg){
	gStyle->SetPaintTextFormat("5.2f");
	// useNiceColorPalette();
	bool hasBjets = (gRegions[reg]->minNbjets > 0);

	fOutputSubDir = "Origins/";
	// make the histograms first. one for ttjets and one for all mc (without signal). this in each channel
	TH2D    *horigin_tt [gNCHANNELS];
	TH2D    *horigin_mc [gNCHANNELS];
	int nbins(12);

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
			float scale = fLumiNorm / fSamples[*sampleInd]->getLumi();
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
		
		TCanvas *c_temp = new TCanvas("Origin" + SSDLDumper::gChanLabel[i], "Origin plot in region " + SSDLDumper::gRegions[reg]->sname, 0, 0, 600, 600);
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
		if (gRegions[reg]->maxHT < 39.) lat->DrawLatex(latX, latY, Form("#splitline{N_{Jets} = 0}{E_{T}^{miss} > %.0f GeV}", gRegions[reg]->minMet));
		else lat->DrawLatex(latX, latY, Form("#splitline{H_{T} > %.0f GeV}{E_{T}^{miss} > %.0f GeV}", gRegions[reg]->minHT, gRegions[reg]->minMet));
		if (hasBjets) lat->DrawLatex(latX, latY-0.09, Form("N_{b-jets} #geq %1d", gRegions[reg]->minNbjets));
		lat->SetTextSize(0.03);
		lat->DrawLatex(0.75, 0.93, Form("Exp. # ev.: %5.2f", total_mc));
		Util::PrintPDF(c_temp, "Origin_" + SSDLDumper::gChanLabel[i] + "_" + SSDLDumper::gRegions[reg]->sname, fOutputDir + fOutputSubDir);

		// TTJets only histos
		float total_tt = horigin_tt[i]->Integral();
		horigin_tt[i]->Scale(100./total_tt);
		horigin_tt[i]->Draw("col text");
		lat->SetTextSize(0.04);
		lat->DrawLatex(0.05, 0.93, "Origins in "+SSDLDumper::gChanLabel[i]+" channel (ttbar only)");
		if (gRegions[reg]->maxHT < 39.) lat->DrawLatex(latX, latY, Form("#splitline{N_{Jets} = 0}{E_{T}^{miss} > %.0f GeV}", gRegions[reg]->minMet));
		else lat->DrawLatex(latX, latY, Form("#splitline{H_{T} > %.0f GeV}{E_{T}^{miss} > %.0f GeV}", gRegions[reg]->minHT, gRegions[reg]->minMet));
		if (hasBjets) lat->DrawLatex(latX, latY-0.09, Form("N_{b-jets} #geq %1d", gRegions[reg]->minNbjets));
		lat->SetTextSize(0.03);
		lat->DrawLatex(0.75, 0.93, Form("Exp. # ev.: %5.2f", total_tt));
		Util::PrintPDF(c_temp, "Origin_TTJets_" + SSDLDumper::gChanLabel[i] + "_" + SSDLDumper::gRegions[reg]->sname, fOutputDir + fOutputSubDir);

		delete c_temp;
	} // end channel loop
	fOutputSubDir = "";
}
void SSDLPlotter::printOrigins(int reg){
	TString filename = fOutputDir + "Origins_"+SSDLDumper::gRegions[reg]->sname+".txt";
	fOUTSTREAM.open(filename.Data(), ios::trunc);
	printMuOriginTable(reg);
	fOUTSTREAM << endl << endl;
	printElOriginTable(reg);
	fOUTSTREAM << endl << endl;
	printEMuOriginTable(reg);
	fOUTSTREAM << endl;
	fOUTSTREAM.close();
}
void SSDLPlotter::printMuOriginTable(int reg){
	fOUTSTREAM << "-------------------------------------------" << endl;
	fOUTSTREAM << " Printing origins for Mu/Mu channel..." << endl;

	printMuOriginHeader("NT20");
	for(gSample i = sample_begin; i < gNSAMPLES; i=gSample(i+1)){
		if(i > QCDMuEnr15) continue;
		print2MuOriginsFromSample(fSamples[i], 2, reg);
	}
	fOUTSTREAM << "-------------------------------------------------------------------------------------------------------------------------" << endl;
	printOriginSummary2L(fMCBGMuEnr, 2, Muon, reg);
	fOUTSTREAM << "=========================================================================================================================" << endl << endl;
	printMuOriginHeader("NT10");
	for(gSample i = sample_begin; i < gNSAMPLES; i=gSample(i+1)){
		if(i > QCDMuEnr15) continue;
		print2MuOriginsFromSample(fSamples[i], 1, reg);
	}
	fOUTSTREAM << "-------------------------------------------------------------------------------------------------------------------------" << endl;
	printOriginSummary2L(fMCBGMuEnr, 1, Muon, reg);
	fOUTSTREAM << "=========================================================================================================================" << endl << endl;
	printMuOriginHeader("NT00");
	for(gSample i = sample_begin; i < gNSAMPLES; i=gSample(i+1)){
		if(i > QCDMuEnr15) continue;
		print2MuOriginsFromSample(fSamples[i], 0, reg);
	}
	fOUTSTREAM << "-------------------------------------------------------------------------------------------------------------------------" << endl;
	printOriginSummary2L(fMCBGMuEnr, 0, Muon, reg);
	fOUTSTREAM << "=========================================================================================================================" << endl << endl;

	printMuOriginHeader("SSTi");
	for(gSample i = sample_begin; i < gNSAMPLES; i=gSample(i+1)){
		if(i > QCDMuEnr15) continue;
		printMuOriginFromSample(fSamples[i], 1, reg);
	}
	fOUTSTREAM << "-------------------------------------------------------------------------------------------------------------------------" << endl;
	printOriginSummary(fMCBGMuEnr, 1, Muon, reg);
	fOUTSTREAM << "=========================================================================================================================" << endl << endl;
	printMuOriginHeader("SSLo");
	for(gSample i = sample_begin; i < gNSAMPLES; i=gSample(i+1)){
		if(i > QCDMuEnr15) continue;
		printMuOriginFromSample(fSamples[i], 2, reg);
	}
	fOUTSTREAM << "-------------------------------------------------------------------------------------------------------------------------" << endl;
	printOriginSummary(fMCBGMuEnr, 2, Muon, reg);
	fOUTSTREAM << "=========================================================================================================================" << endl << endl;
	printMuOriginHeader("Z Ti");
	for(gSample i = sample_begin; i < gNSAMPLES; i=gSample(i+1)){
		if(i > QCDMuEnr15) continue;
		printMuOriginFromSample(fSamples[i], 3, reg);
	}
	fOUTSTREAM << "-------------------------------------------------------------------------------------------------------------------------" << endl;
	printOriginSummary(fMCBGMuEnr, 3, Muon, reg);
	fOUTSTREAM << "=========================================================================================================================" << endl << endl;
	printMuOriginHeader("Z Lo");
	for(gSample i = sample_begin; i < gNSAMPLES; i=gSample(i+1)){
		if(i > QCDMuEnr15) continue;
		printMuOriginFromSample(fSamples[i], 4, reg);
	}
	fOUTSTREAM << "-------------------------------------------------------------------------------------------------------------------------" << endl;
	printOriginSummary(fMCBGMuEnr, 4, Muon, reg);
	fOUTSTREAM << "=========================================================================================================================" << endl << endl;
}
void SSDLPlotter::printElOriginTable(int reg){
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
void SSDLPlotter::printEMuOriginTable(int reg){
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

void SSDLPlotter::printMuOriginFromSample(Sample *S, int toggle, int reg, gHiLoSwitch hilo){
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
void SSDLPlotter::printElOriginFromSample(Sample *S, int toggle, int reg, gHiLoSwitch hilo){
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

void SSDLPlotter::print2MuOriginsFromSample(Sample *S, int toggle, int reg, gHiLoSwitch hilo){
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
void SSDLPlotter::print2ElOriginsFromSample(Sample *S, int toggle, int reg, gHiLoSwitch hilo){
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
void SSDLPlotter::printEMuOriginsFromSample(Sample *S, int toggle, int reg, gHiLoSwitch hilo){
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

void SSDLPlotter::printOriginSummary(vector<int> samples, int toggle, gChannel chan, int reg, gHiLoSwitch hilo){
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

		float scale = fLumiNorm / S->getLumi();
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
void SSDLPlotter::printOriginSummary2L(vector<int> samples, int toggle, gChannel chan, int reg, gHiLoSwitch hilo){
	if(toggle != 0 && toggle != 1 && toggle != 10 && toggle != 2) return;
	TH1D *histosum1 = new TH1D("SST_Origin_Sum1", "SSTOrigin",  12, 0, 12);
	TH1D *histosum2 = new TH1D("SST_Origin_Sum2", "SSTOrigin",  12, 0, 12);
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

		float scale = fLumiNorm / S->getLumi();
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
void SSDLPlotter::drawTopLine(float rightedge, float scale, float leftedge){
	fLatex->SetTextFont(62);
	fLatex->SetTextSize(scale*0.05);
	fLatex->DrawLatex(leftedge,0.92, "CMS Preliminary");
	fLatex->SetTextFont(42);
	fLatex->SetTextSize(scale*0.04);
	fLatex->DrawLatex(rightedge,0.92, Form("L = %4.2f fb^{-1} at #sqrt{s} = 8 TeV", fLumiNorm/1000.));
	// fLatex->DrawLatex(0.70,0.92, Form("L_{int.} = %4.0f pb^{-1}", fLumiNorm));
	return;
}
void SSDLPlotter::drawTopLineSim(float rightedge, float scale, float leftedge){
	fLatex->SetTextFont(62);
	fLatex->SetTextSize(scale*0.05);
	fLatex->DrawLatex(leftedge,0.92, "CMS Simulation");
	fLatex->SetTextFont(42);
	fLatex->SetTextSize(scale*0.04);
	fLatex->DrawLatex(rightedge,0.92, Form("L = %4.2f fb^{-1} at #sqrt{s} = 8 TeV", fLumiNorm/1000.));
	// fLatex->DrawLatex(0.70,0.92, Form("L_{int.} = %4.0f pb^{-1}", fLumiNorm));
	return;
}
void SSDLPlotter::drawDiffCuts(int j){
	// fLatex->SetTextFont(42);
	// fLatex->SetTextSize(0.03);
	// if(j==8)       fLatex->DrawLatex(0.13,0.85, "E_{T}^{miss} > 30 GeV, H_{T} > 200 GeV, N_{Jets} #geq 2");
	// if(j>4&&j<8)   fLatex->DrawLatex(0.13,0.85, "E_{T}^{miss} > 30 GeV, H_{T} > 80 GeV, N_{Jets} #geq 2");
	// if(j==0)       fLatex->DrawLatex(0.13,0.85, "E_{T}^{miss} > 50 GeV");
	// if(j==1||j==4) fLatex->DrawLatex(0.13,0.85, "E_{T}^{miss} > 120 GeV");
	// if(j==2)       fLatex->DrawLatex(0.13,0.85, "H_{T} > 200 GeV, N_{Jets} #geq 2");
	// if(j==3)       fLatex->DrawLatex(0.13,0.85, "H_{T} > 450 GeV, N_{Jets} #geq 2");
	return;
}
void SSDLPlotter::drawRegionSel(int reg){
	TLatex *lat = new TLatex();
	lat->SetNDC(kTRUE);
	lat->SetTextColor(kBlack);
	lat->SetTextSize(0.03);
	TString hadString;
	if (gRegions[reg]->minNbjmed > 0 ) {
		// cout << Region::sname[reg] << ": minNbjets: "<<Region::minNbjets[reg] << " maxnbjmed: " << Region::maxNbjmed[reg] << " maxnbjets: " << Region::maxNbjets[reg] << endl;
		// lat->DrawLatex(0.52, 0.85, Form("H_{T} > %.0f GeV, N_{bTags}^{med} #geq %1d",Region::minHT[reg], Region::minNbjmed[reg]));
		hadString = Form("H_{T} > %.0f GeV, N_{bTags}^{med} #geq %1d",gRegions[reg]->minHT, gRegions[reg]->minNbjmed);
	}
	else if (gRegions[reg]->maxNbjets == gRegions[reg]->minNbjets ) {
		// cout << Region::sname[reg] << ": minNbjets: "<<Region::minNbjets[reg] << " maxnbjmed: " << Region::maxNbjmed[reg] << " maxnbjets: " << Region::maxNbjets[reg] << endl;
		// lat->DrawLatex(0.52, 0.85, Form("H_{T} > %.0f GeV, N_{bTags}^{med} = 0",Region::minHT[reg]));
		hadString = Form("H_{T} > %.0f GeV, N_{bTags} = %1d",gRegions[reg]->minHT, gRegions[reg]->maxNbjets);
	}
	else {
		hadString = Form("H_{T} > %.0f GeV, N_{jets} #geq %1d",gRegions[reg]->minHT, gRegions[reg]->minNjets);
	}
	lat->DrawLatex(0.52, 0.85, hadString);
	if (gRegions[reg]->maxNjets < 10) {
		lat->DrawLatex(0.52, 0.75, Form("N_{jets} #leq %1d ",gRegions[reg]->maxNjets));
	}

	// lat->DrawLatex(0.52, 0.85, Form("H_{T} > %.0f GeV, N_{bTags}(med) #geq %1d",Region::minHT[reg], Region::minNbjmed[reg]));
	lat->DrawLatex(0.52, 0.80, Form("E_{T}^{miss} > %.0f GeV",    gRegions[reg]->minMet));
	if (reg == gRegion["HT80MET30bpp"]) 
	  lat->DrawLatex(0.52,0.75, Form("Only ++"));
// 	if(reg == TTbarWPresel){
// 		lat->DrawLatex(0.55,0.85, Form("H_{T} > %.0f GeV, N_{Jets} #geq %1d",        Region::minHT[reg], Region::minNjets[reg]));
// 		lat->DrawLatex(0.55,0.80, Form("N_{bTags} (medium) #geq %1d",                Region::minNbjmed[reg]));

// 	}
	// else if(Region::maxHT[reg] < 19.)    lat->DrawLatex(0.55,0.85, "N_{Jets} = 0");
	// else if(Region::minHT[reg] == 0.)    lat->DrawLatex(0.55,0.85, Form("H_{T} #geq %.0f GeV, N_{Jets} #geq %1d", Region::minHT[reg], Region::minNjets[reg]));
	// else                                 lat->DrawLatex(0.55,0.85, Form("H_{T} > %.0f GeV, N_{Jets} #geq %1d",    Region::minHT[reg], Region::minNjets[reg]));
	// // else                               lat->DrawLatex(0.55,0.85, Form("H_{T} > %.0f GeV, N_{Jets} #geq %1d", Region::minHT[reg], Region::minNjets[reg]));
	// // if(     Region::minMet[reg] == 0.)   lat->DrawLatex(0.55,0.80, Form("E_{T}^{miss} #geq %.0f GeV", Region::minMet[reg]));
	// // if(     Region::minMet[reg] == 0.);
	// // else if(Region::minMet[reg] > 0. )   
	// // if(     Region::minNbjets[reg] > 0 ) lat->DrawLatex(0.55,0.75, Form("N_{bTags} #geq %1d",         Region::minNbjets[reg]));
	// if(     Region::minNbjmed[reg] > 0 ) lat->DrawLatex(0.55,0.80, Form("N_{bTags} (medium) #geq %1d",         Region::minNbjets[reg]));
}

void SSDLPlotter::msugraKfacs(TFile * results){
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
void SSDLPlotter::msugraLOxsecs(TFile * results){
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
void SSDLPlotter::msugraNLOxsecs(TFile * results) {
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
void SSDLPlotter::scanMSUGRA( const char * filestring){

	bool isLeptonSkim = true;
	TString specialString = "_iso0p1_newZveto";

	TString lepString = "";
	if (isLeptonSkim) lepString = "_leptonSkim";

	fC_minMu1pt = 20.;
	fC_minMu2pt = 10.;
	fC_minEl1pt = 20.;
	fC_minEl2pt = 10.;
	fC_minHT    = 0.;
	fC_maxHT    = 29.;
	fC_minMet   = 120.;
	fC_maxMet   = 7000.;
	fC_minNjets = 0;
	
	TString htString;
	TString htTitleString;

	if (fC_maxHT < 40.) {
		htString = "HT0JV";
		htTitleString = "N_{Jets} = 0";
	}
	else if (fC_maxHT == 7000.) {
		htString = "HT0";
		htTitleString = "H_{T} > 0 GeV";
	}
	else {
		htString = Form("HT%3.0f", fC_maxHT);
		htTitleString = Form("H_{T} > %3.0f GeV", fC_maxHT);
	}
	
	TFile * res_ = new TFile(Form("msugraSSDL/mSugraresults_"+htString+"_MET%3.0f_PT%2.0f_%2.0f"+lepString+specialString+".root", fC_minMet, fC_minMu1pt, fC_minMu2pt), "RECREATE", "res_");

	bool verbose = false;
	if (verbose) fOUTSTREAM.open("msugraSSDL/output"+lepString+specialString+".txt", ios::trunc);
	// comment the next two lines in case you already have the necessary files at hand.
	// you can also just let them in, it doesn't really affect the performance much
	SSDLPlotter::msugraKfacs(res_);
	SSDLPlotter::msugraLOxsecs(res_);
	SSDLPlotter::msugraNLOxsecs(res_);

	TH2D  * nPass_     = new TH2D("msugra_nPass"  , "msugra_nPass"  , gM0bins , gM0min+10 , gM0max+10 , gM12bins , gM12min+10 , gM12max+10);
	TH2D  * wPass_     = new TH2D("msugra_wPass"  , "msugra_wPass"  , gM0bins , gM0min+10 , gM0max+10 , gM12bins , gM12min+10 , gM12max+10);
	TH2D  * wCount_    = new TH2D("msugra_wcount" , "msugra_wcount" , gM0bins , gM0min+10 , gM0max+10 , gM12bins , gM12min+10 , gM12max+10);
	TH2D  * yield_     = new TH2D("msugra_yield"  , "msugra_yield"  , gM0bins , gM0min+10 , gM0max+10 , gM12bins , gM12min+10 , gM12max+10);

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

	TH2D  * filterEff_;
	if (isLeptonSkim) {
		TFile * filterEffFile_ = new TFile("msugraSSDL/FilterEfficiencyv3.root", "READ", "filterEffFile_");
		filterEff_ = (TH2D* ) filterEffFile_->Get("FilterEfficiency");
	}
//////--------------------------------------------------------

	TTree * tree_= (TTree *) file_->Get("Analysis");
	tree_->ResetBranchAddresses();

	Init(tree_);
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
				if (verbose) fOUTSTREAM << Form("MuMu - ev %11.0d - m0 %4.0f - m12 %4.0f - process %2.0i - HT(#J/#bJ) %6.2f(%1d/%1d) MET %6.2f Pt1 %6.2f Pt2 %6.2f Charge %2d", Event, m0, m12, process, getHT(), getNJets(), getNBTags(), pfMET, MuPt[mu1], MuPt[mu2], MuCharge[mu1]) << endl ;
				n_tot++;
				int bin = nlo_[process-1]->FindBin(m0, m12);
				float weight = fLumiNorm * nlo_[process-1]->GetBinContent(bin)/ kCounts[process-1]->GetBinContent(bin);
				if (isLeptonSkim){
					int newbin = filterEff_->FindBin(m0, m12);
					//if (bin != newbin) cout << "Binning is different in the leptonic filter efficiency skim!" << endl;
					//float pointFilterEff = filterEff_->GetBinContent(newbin);
					//cout << " m0 , m12 : " << m0 << " , " << m12 << " filterEff: " << pointFilterEff << endl;
					weight *= filterEff_->GetBinContent(newbin);
				}
				if (verbose) fOUTSTREAM << " m0: " << m0 << " m12: " << m12 << " process " << process << " weight: " << weight << endl;
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
				if (verbose) fOUTSTREAM << Form("ElMu - ev %11.0d - m0 %4.0f - m12 %4.0f - process %2.0i - HT(#J/#bJ) %6.2f(%1d/%1d) MET %6.2f Pt1 %6.2f Pt2 %6.2f Charge %2d", Event, m0, m12, process, getHT(), getNJets(), getNBTags(), pfMET, MuPt[mu], ElPt[el], MuCharge[mu]) << endl ;
				n_tot++;
				int bin = nlo_[process-1]->FindBin(m0, m12);
				float weight = fLumiNorm * nlo_[process-1]->GetBinContent(bin)/ kCounts[process-1]->GetBinContent(bin);
				if (isLeptonSkim){
					int newbin = filterEff_->FindBin(m0, m12);
					weight *= filterEff_->GetBinContent(newbin);
				}
				if (verbose) fOUTSTREAM << " m0 " << m0 << " m12: " << m12 << " process " << process << " weight: " << weight << endl;
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
				if (verbose) fOUTSTREAM << Form("ElEl - ev %11.0d - m0 %4.0f - m12 %4.0f - process %2.0i - HT(#J/#bJ) %6.2f(%1d/%1d) MET %6.2f Pt1 %6.2f Pt2 %6.2f Charge %2d", Event, m0, m12, process, getHT(), getNJets(), getNBTags(), pfMET, ElPt[el1], ElPt[el2], ElCharge[el1]) << endl ;
				n_tot++;
				int bin = nlo_[process-1]->FindBin(m0, m12);
				float weight = fLumiNorm * nlo_[process-1]->GetBinContent(bin)/ kCounts[process-1]->GetBinContent(bin);
				if (isLeptonSkim){
					int newbin = filterEff_->FindBin(m0, m12);
					weight *= filterEff_->GetBinContent(newbin);
				}
				if (verbose) fOUTSTREAM << " m0 " << m0 << " m12: " << m12 << " process " << process << " weight: " << weight << endl;
				nPass_->Fill(m0, m12);
				yield_->Fill(m0, m12, weight);
				nEE++;
			}
		}
	}
	cout << Form("Total Number of SS events (with MET > %3.0f GeV and HT > %3.0f GeV and HT < %3.0f): ", fC_minMet, fC_minHT, fC_maxHT) << n_tot << endl;
	cout << "nEE: " << nEE << " nEM: " << nEM << " nMM: " << nMM << endl;
	if (verbose) fOUTSTREAM << "Total Number of SS events (with MET > 120 GeV and HT > 450 GeV): " << n_tot << endl;
	if (verbose) fOUTSTREAM << "Total number of tight pairs: " << tightTot << " total number of signal pairs: " << signalTot << " resulting efficiency: " << signalTot/tightTot << endl;
	if (verbose) fOUTSTREAM << "nEE: " << nEE << " nEM: " << nEM << " nMM: " << nMM << endl;

	TH2D *eff_   = new TH2D("msugra_eff"   , "msugra_eff"  , gM0bins, gM0min+10, gM0max+10, gM12bins, gM12min+10, gM12max+10);
	eff_->Divide(nPass_, count_, 1., 1.);

	res_->cd();

	for (int i = 0 ; i< 10; i++){
	    nlo_[i]->Write();
	}
	yield_->Write();
	count_->Write();
	nPass_->Write();
	eff_->Write();

}
void SSDLPlotter::scanSMS( const char * filestring, int reg){

        bool verbose = false;
	if (verbose) fOUTSTREAM.open( "SMSoutput_"+gRegions[reg]->sname+".txt" );

	// put the x-values that are in the scan here. those are percentage numbers. i.e. if x == 0.05 put 5
	int nx(3);
	float xvals[nx];
	xvals[0] = 0.05;
	xvals[1] = 0.50;
	xvals[2] = 0.95;

	TH2D * TChi_yield_      [nx];
	TH2D * TChiRight_yield_ [nx];

	for (int i = 0; i<nx; i++) {
		TChi_yield_              [i] = new TH2D(Form("TChi_yield%.0f", 100*xvals[i])      , Form("TChi_yield%.0f", 100*xvals[i])      , 101, -5, 1005    , 101 , -5 , 1005);
		TChiRight_yield_         [i] = new TH2D(Form("TChiRight_yield%.0f", 100* xvals[i]), Form("TChiRight_yield%.0f", 100*xvals[i]) , 101, -5, 1005    , 101 , -5 , 1005);
		TChi_yield_              [i]->Sumw2();
		TChiRight_yield_         [i]->Sumw2();
	}

	int nSyst(9);
	TString systs[nSyst];
	systs[0] = "norm";
	systs[1] = "metsmear";
	systs[2] = "lepup";
	systs[3] = "lepdown";
	systs[4] = "JESup";
	systs[5] = "JESdown";
	systs[6] = "JER";
	systs[7] = "METup";
	systs[8] = "METdown";

	TH2D  * TChi_nPass_     [nx][nSyst];
	TH2D  * TChiRight_nPass_[nx][nSyst];
	for (int i = 0; i<nx; i++) {
		for (int j = 0; j<nSyst; j++) {
			TChi_nPass_      [i][j] = new TH2D(Form("TChi_nPass_x%.0f_"+systs[j]     , 100*xvals[i]), Form("TChi_nPass_x%.0f_"+systs[j]     , 100*xvals[i]) , 101, -5, 1005, 101, -5, 1005);
			TChiRight_nPass_ [i][j] = new TH2D(Form("TChiRight_nPass_x%.0f_"+systs[j], 100*xvals[i]), Form("TChiRight_nPass_x%.0f_"+systs[j], 100*xvals[i]) , 101, -5, 1005, 101, -5, 1005);
			TChi_nPass_      [i][j]->Sumw2();
			TChiRight_nPass_ [i][j]->Sumw2();
		}
	}

	// get the histo with the x-secs
	TFile * xsecFile_ = new TFile("msugraSSDL/C1N2_8TeV_xsecs.root", "READ", "xsecFile_");
	TH1D  * xsecs     = (TH1D *) xsecFile_->Get("xsecs");

	// get the histo with the count for each point
	TFile * file_ = TFile::Open(filestring); //, "READ", "file_");
	TH2D  * TChi_nTot_ [nx];
	TH2D  * TChiRight_nTot_[nx];
	for (int i = 0; i<nx; i++) {
		TChi_nTot_[i]      = (TH2D  *) file_->Get(Form("ModelCount%.0f", 100*xvals[i]));
		TChiRight_nTot_[i] = (TH2D  *) file_->Get(Form("RightHandedSlepCount%.0f", 100*xvals[i]));
		TChi_nTot_[i]     ->Sumw2();
		TChiRight_nTot_[i]->Sumw2();
	}
	// get the Analysis tree from the file, initialize it
	TTree * tree_ = (TTree *) file_->Get("Analysis");
	tree_->ResetBranchAddresses();
	Init(tree_);

	// just for counting, really
	Long64_t tot_events = tree_->GetEntriesFast();
	cout << "Total Number of entries: " << tot_events << endl;
	int n_tot = 0;
	int signalTot(0);
	int nEE(0), nEM(0), nMM(0);
	int nSlepSnu(0), nSlepSlep(0);

	bool doSystematic = true;
	// make a dummy sample for the systematic functions
	Sample *S = new Sample();
	S->datamc = 1;
	// need to initialize the TRandom object for the systematic studies
	SSDLDumper::fRand3 = new TRandom3(0);

	for (Long64_t jentry=0; jentry<tree_->GetEntriesFast();jentry++) {
		setRegionCuts(reg);
		printProgress(jentry, tot_events, "SMS Scan "+gRegions[reg]->sname);
		for (int i = 0; i<nSyst; i++) {
			tree_->GetEntry(jentry); // have to reload the entry for each systematic
			int x=0;                   // get the value of x. yes, in the treee it's m0
			if (m0 == xvals[0]) x = 0;
			if (m0 == xvals[1]) x = 1;
			if (m0 == xvals[2]) x = 2;
			 
			if (!doSystematic && i!=0) continue;
			if (i == 1) smearMET(S);
			if (i == 2) scaleLeptons(S, 1);
			if (i == 3) scaleLeptons(S, 2);
			if (i == 4) smearJetPts(S, 1);
			if (i == 5) smearJetPts(S, 2);
			if (i == 6) smearJetPts(S, 3);
			if (i == 7) scaleMET(S, 0);
			if (i == 8) scaleMET(S, 1);
			// fancier now if (i == 9) {pfMET = pfMET + 0.05 * pfMET; }
			// fancier now if (i ==10) {pfMET = pfMET - 0.05 * pfMET; }
			int   xsecBin      = xsecs->FindBin(mGlu);
			float nloXsec      = xsecs->GetBinContent(xsecBin);
			int   nGenBin      = TChi_nTot_[x]->FindBin(mGlu, mLSP);
			float nGen         = TChi_nTot_[x]->GetBinContent(nGenBin);
			int   nGenBinRight = TChiRight_nTot_[x]->FindBin(mGlu, mLSP);
			float nGenRight    = TChiRight_nTot_[x]->GetBinContent(nGenBin);
			float weight       = fLumiNorm * nloXsec / nGen;
			float weightRight  = fLumiNorm * nloXsec / nGenRight;

			int mu1(-1), mu2(-1);
			if( isSSLLMuEvent(mu1, mu2) ){ // Same-sign loose-loose di muon event
				// if(pass3rdVeto && !passes3rdLepVeto()) continue;
				if(isTightMuon(mu1) &&  isTightMuon(mu2) ){ // Tight-tight
					n_tot++;
					if ( IsSignalMuon[mu1] != 1 || IsSignalMuon[mu2] != 1 ) continue;
					signalTot++;
					if (verbose) fOUTSTREAM << Form("MM - mGlu %4.0f - mLSP %4.0f - HT %4.2f - MET %6.2f Pt1 %6.2f Pt2 %6.2f | %2d | SlepSnu: %2i", mGlu, mLSP, getHT(), pfMET, MuPt[mu1], MuPt[mu2], MuCharge[mu1], isTChiSlepSnu) << endl ;
					TChi_nPass_[x] [i]-> Fill(mGlu, mLSP);
					if (i==0) {
						TChi_yield_ [x]-> Fill(mGlu, mLSP, weight * gMMTrigScale);
					}
					if (isRightHanded == 1 && isTChiSlepSnu == 0){
						TChiRight_nPass_[x] [i]-> Fill(mGlu, mLSP);
						if (i==0) {
							TChiRight_yield_ [x]-> Fill(mGlu, mLSP, weightRight * gMMTrigScale);
						}
					}
					nMM++;
					continue;
				}
				resetHypLeptons();
			}
			int mu(-1), el(-1);
			if( isSSLLElMuEvent(mu, el) ){
				// if(pass3rdVeto && !passes3rdLepVeto()) continue;
				if(  isTightElectron(el) &&  isTightMuon(mu) ){ // Tight-tight
					n_tot++;
					if ( IsSignalMuon[mu] != 1 || IsSignalElectron[el] != 1 ) continue;
					signalTot++;
					if (verbose) fOUTSTREAM << Form("EM - mGlu %4.0f - mLSP %4.0f - HT %4.2f - MET %6.2f Pt1 %6.2f Pt2 %6.2f | %2d | SlepSnu: %2i", mGlu, mLSP, getHT(), pfMET, MuPt[mu], ElPt[el], MuCharge[mu1], isTChiSlepSnu) << endl ;
					TChi_nPass_[x] [i]-> Fill(mGlu, mLSP);
					if (i==0) {
						TChi_yield_ [x]-> Fill(mGlu, mLSP, weight * gEMTrigScale);
					}
					if (isRightHanded == 1 && isTChiSlepSnu == 0){
						TChiRight_nPass_[x] [i]-> Fill(mGlu, mLSP);
						if (i==0) {
							TChiRight_yield_ [x]-> Fill(mGlu, mLSP, weightRight * gEMTrigScale);
						}
					}
					nEM++;
					continue;
				}
				resetHypLeptons();
			}
			int el1(-1), el2(-1);
			if( isSSLLElEvent(el1, el2) ){
				// if(pass3rdVeto && !passes3rdLepVeto()) continue;
				if(  isTightElectron(el1) &&  isTightElectron(el2) ){ // Tight-tight
					n_tot++;
					if ( IsSignalElectron[el1] != 1 || IsSignalElectron[el2] != 1 ) continue;
					signalTot++;
					if (verbose) fOUTSTREAM << Form("EE - mGlu %4.0f - mLSP %4.0f - HT %4.2f - MET %6.2f Pt1 %6.2f Pt2 %6.2f | %2d | SlepSnu: %2i", mGlu, mLSP, getHT(), pfMET, ElPt[el1], ElPt[el2], ElCharge[el1], isTChiSlepSnu) << endl ;
					TChi_nPass_[x] [i]-> Fill(mGlu, mLSP);
					if (i==0) {
						TChi_yield_ [x]-> Fill(mGlu, mLSP, weight * gEETrigScale);
					}
					if (isRightHanded == 1 && isTChiSlepSnu == 0){
						TChiRight_nPass_[x] [i]-> Fill(mGlu, mLSP);
						if (i==0) {
							TChiRight_yield_ [x]-> Fill(mGlu, mLSP, weightRight* gEETrigScale);
						}
					}
					nEE++;
				}
			}
		}
	}
	if (verbose) cout << "Total Number of SS events: " << n_tot << endl;
	if (verbose) cout << "nEE: " << nEE << " nEM: " << nEM << " nMM: " << nMM << endl;
	if (verbose) fOUTSTREAM << "Total number of tight pairs: " << n_tot << " total number of signal pairs: " << signalTot << " resulting efficiency: " << signalTot/n_tot << endl;
	if (verbose) fOUTSTREAM << "nEE: " << nEE << " nEM: " << nEM << " nMM: " << nMM << endl;
	if (verbose) fOUTSTREAM << "nSlepSnu: " << nSlepSnu << " nSlepSlep: " << nSlepSlep << endl;

	TH2D * TChi_eff_      [nx] [nSyst];
	TH2D * TChiRight_eff_ [nx] [nSyst];
	for (int x=0; x<nx;x++) {
		for (int i=0; i<nSyst;i++) {
			TChi_eff_      [x][i] = new TH2D(Form("TChi_eff_x%.0f_"+systs[i]      , 100*xvals[x]) , Form("TChi_eff_x%.0f_"+systs[i]     , 100*xvals[x]) , 101, -5, 1005, 101, -5, 1005);
			TChiRight_eff_ [x][i] = new TH2D(Form("TChiRight_eff_x%.0f_"+systs[i] , 100*xvals[x]) , Form("TChiRight_eff_x%.0f_"+systs[i], 100*xvals[x]) , 101, -5, 1005, 101, -5, 1005);
			TChi_eff_      [x][i]->Sumw2();
			TChiRight_eff_ [x][i]->Sumw2();
			TChi_eff_      [x][i]->Divide(TChi_nPass_[x][i]      , TChi_nTot_[x]      , 1. , 1.);
			TChiRight_eff_ [x][i]->Divide(TChiRight_nPass_[x][i] , TChiRight_nTot_[x] , 1. , 1.);
		}
	}

	TFile * res_;
	res_ = new TFile("SMSresults_"+gRegions[reg]->sname+".root", "RECREATE", "res_");
	res_   -> cd();

	for (int x=0; x<nx;x++) {
		TChi_yield_     [x]->Write();
		TChi_nTot_      [x]->Write();
		TChiRight_yield_[x]->Write();
		TChiRight_nTot_ [x]->Write();

		for (int i=0; i<nSyst; i++) {
			TChi_eff_       [x][i]-> Write();
			TChi_nPass_     [x][i]->Write();
			TChiRight_eff_  [x][i]-> Write();
			TChiRight_nPass_[x][i]->Write();
		}
	}
	file_->Close();
	res_->Close();

	//	delete S;
	delete res_;
}
void SSDLPlotter::scanModelGeneric( const char * filestring, int reg){

	TString name = "T6ttWW";
	bool verbose = false;
	if (verbose) fOUTSTREAM.open( name+"_"+gRegions[reg]->sname+".txt" );


	TH2D * Model_yield_ ;
	Model_yield_  = new TH2D(name+"_yield"      , name+"_yield"  , 101, -5, 1005    , 101 , -5 , 1005);
	Model_yield_ ->Sumw2();

	int nSyst(9);
	TString systs[nSyst];
	systs[0] = "norm";
	systs[1] = "metsmear";
	systs[2] = "lepup";
	systs[3] = "lepdown";
	systs[4] = "JESup";
	systs[5] = "JESdown";
	systs[6] = "JER";
	systs[7] = "METup";
	systs[8] = "METdown";

	TH2D  * Model_nPass_    [nSyst];
	for (int j = 0; j<nSyst; j++) {
		Model_nPass_      [j] = new TH2D(name+"_nPass_"+systs[j], name+"_nPass_"+systs[j] , 101, -5, 1005, 101, -5, 1005);
		Model_nPass_      [j]->Sumw2();
	}

	// get the histo with the x-secs


	TFile * xsecFile_ = new TFile("/shome/mdunser/xsecfiles/"+name+"_xsecs.root", "READ", "xsecFile_");
	TH1D  * xsecs     = (TH1D *) xsecFile_->Get("xsecs");
	TF1   * xsecfit =  xsecs->GetFunction("xsec_fit1");

	// get the histo with the count for each point
	TFile * file_ = TFile::Open(filestring); //, "READ", "file_");
	TH2D  * Model_nTot_ ;
	Model_nTot_ = (TH2D  *) file_->Get("ModelCountAll");
	Model_nTot_ ->Sumw2();
	// get the Analysis tree from the file, initialize it
	TTree * tree_ = (TTree *) file_->Get("Analysis");
	tree_->ResetBranchAddresses();
	Init(tree_);

	// just for counting, really
	Long64_t tot_events = tree_->GetEntriesFast();
	cout << "Total Number of entries: " << tot_events << endl;
	int n_tot = 0;
	int signalTot(0);
	int nEE(0), nEM(0), nMM(0);
	int nSlepSnu(0), nSlepSlep(0);

	bool doSystematic = true;
	// make a dummy sample for the systematic functions
	Sample *S = new Sample();
	S->datamc = 1;
	// need to initialize the TRandom object for the systematic studies
	SSDLDumper::fRand3 = new TRandom3(0);

	for (Long64_t jentry=0; jentry<tree_->GetEntriesFast();jentry++) {
		setRegionCuts(reg);
		printProgress(jentry, tot_events, name+" Scan "+gRegions[reg]->sname);
		for (int i = 0; i<nSyst; i++) {
			tree_->GetEntry(jentry); // have to reload the entry for each systematic
			 
			if (!doSystematic && i!=0) continue;
			if (i == 1) smearMET(S);
			if (i == 2) scaleLeptons(S, 1);
			if (i == 3) scaleLeptons(S, 2);
			if (i == 4) smearJetPts(S, 1);
			if (i == 5) smearJetPts(S, 2);
			if (i == 6) smearJetPts(S, 3);
			if (i == 7) scaleMET(S, 0);
			if (i == 8) scaleMET(S, 1);
			float nloXsec      = xsecfit->Eval(mGlu);
			int   nGenBin      = Model_nTot_->FindBin(mGlu, mChi);
			float nGen         = Model_nTot_->GetBinContent(nGenBin);
			float weight       = fLumiNorm * nloXsec / nGen;

			int mu1(-1), mu2(-1);
			if( isSSLLMuEvent(mu1, mu2) ){ // Same-sign loose-loose di muon event
				// if(pass3rdVeto && !passes3rdLepVeto()) continue;
				if(isTightMuon(mu1) &&  isTightMuon(mu2) ){ // Tight-tight
					n_tot++;
					if ( IsSignalMuon[mu1] != 1 || IsSignalMuon[mu2] != 1 ) continue;
					signalTot++;
					if (verbose) fOUTSTREAM << Form("MM - mGlu %4.0f - mLSP %4.0f - HT %4.2f - MET %6.2f Pt1 %6.2f Pt2 %6.2f | %2d | SlepSnu: %2i", mGlu, mLSP, getHT(), pfMET, MuPt[mu1], MuPt[mu2], MuCharge[mu1], isTChiSlepSnu) << endl ;
					Model_nPass_[i]-> Fill(mGlu, mChi);
					if (i==0) {
						Model_yield_-> Fill(mGlu, mChi, weight * gMMTrigScale);
					}
					nMM++;
					continue;
				}
				resetHypLeptons();
			} // end loop on mumu channel
			int mu(-1), el(-1);
			if( isSSLLElMuEvent(mu, el) ){
				// if(pass3rdVeto && !passes3rdLepVeto()) continue;
				if(  isTightElectron(el) &&  isTightMuon(mu) ){ // Tight-tight
					n_tot++;
					if ( IsSignalMuon[mu] != 1 || IsSignalElectron[el] != 1 ) continue;
					signalTot++;
					if (verbose) fOUTSTREAM << Form("EM - mGlu %4.0f - mLSP %4.0f - HT %4.2f - MET %6.2f Pt1 %6.2f Pt2 %6.2f | %2d | SlepSnu: %2i", mGlu, mLSP, getHT(), pfMET, MuPt[mu], ElPt[el], MuCharge[mu1], isTChiSlepSnu) << endl ;
					Model_nPass_[i]-> Fill(mGlu, mChi);
					if (i==0) {
						Model_yield_-> Fill(mGlu, mChi, weight * gEMTrigScale);
					}
					nEM++;
					continue;
				}
				resetHypLeptons();
			} // end loop on emu channel
			int el1(-1), el2(-1);
			if( isSSLLElEvent(el1, el2) ){
				// if(pass3rdVeto && !passes3rdLepVeto()) continue;
				if(  isTightElectron(el1) &&  isTightElectron(el2) ){ // Tight-tight
					n_tot++;
					if ( IsSignalElectron[el1] != 1 || IsSignalElectron[el2] != 1 ) continue;
					signalTot++;
					if (verbose) fOUTSTREAM << Form("EE - mGlu %4.0f - mLSP %4.0f - HT %4.2f - MET %6.2f Pt1 %6.2f Pt2 %6.2f | %2d | SlepSnu: %2i", mGlu, mLSP, getHT(), pfMET, ElPt[el1], ElPt[el2], ElCharge[el1], isTChiSlepSnu) << endl ;
					Model_nPass_[i]-> Fill(mGlu, mChi);
					if (i==0) {
						Model_yield_-> Fill(mGlu, mChi, weight * gEETrigScale);
					}
					nEE++;
				}
			} // end loop on ee channel
		} // end loop on all systematics
	} // end loop on all events
	if (verbose) cout << "Total Number of SS events: " << n_tot << endl;
	if (verbose) cout << "nEE: " << nEE << " nEM: " << nEM << " nMM: " << nMM << endl;
	if (verbose) fOUTSTREAM << "Total number of tight pairs: " << n_tot << " total number of signal pairs: " << signalTot << " resulting efficiency: " << signalTot/n_tot << endl;
	if (verbose) fOUTSTREAM << "nEE: " << nEE << " nEM: " << nEM << " nMM: " << nMM << endl;

	TH2D * Model_eff_     [nSyst];
	for (int i=0; i<nSyst;i++) {
		Model_eff_ [i] = new TH2D(name+"_eff_"+systs[i], name+"_eff_"+systs[i], 101, -5, 1005, 101, -5, 1005);
		Model_eff_ [i]->Sumw2();
		Model_eff_ [i]->Divide(Model_nPass_[i], Model_nTot_ , 1. , 1.);
	}

	TFile * res_;
	res_ = new TFile(fOutputDir+name+"_results_"+gRegions[reg]->sname+".root", "RECREATE", "res_");
	res_   -> cd();

	Model_yield_ ->Write();
	Model_nTot_  ->Write();

	for (int i=0; i<nSyst; i++) {
		Model_eff_   [i]->Write();
		Model_nPass_ [i]->Write();
	}
	file_->Close();
	res_->Close();

	//	delete S;
	delete res_;
}
