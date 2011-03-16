/*****************************************************************************
*   Collection of tools for producing plots for Muon Fake Rate Analysis      *
*                                                                            *
*                                                  (c) 2010 Benjamin Stieger *
*****************************************************************************/
#include "MuonPlotter.hh"
#include "helper/AnaClass.hh"
#include "helper/Utilities.hh"
#include "helper/FPRatios.hh"
#include "helper/Monitor.hh"

#include "TLorentzVector.h"
#include "TGraphAsymmErrors.h"
#include "TEfficiency.h"

#include <iomanip>

using namespace std;

//////////////////////////////////////////////////////////////////////////////////
// Global parameters:

static const float gMMU = 0.1057;
static const float gMEL = 0.0005;
static const float gMZ  = 91.2;

static const TEfficiency::EStatOption gStatOpt = TEfficiency::kBBayesian; // kBBayesian (alpha, beta), kBUniform (1,1), kBJeffrey (0.5, 0.5)
static const double gStatBetaAlpha = 1.;
static const double gStatBetaBeta  = 1.;

// Muon Binning //////////////////////////////////////////////////////////////////
static const int gNMuPtbins_0  = 6;
static const double gMuPtbins_0 [gNMuPtbins_0+1]  = {20., 30., 40., 50., 60., 70., 100.};

static const int gNMuPtbins_1  = 7;
static const double gMuPtbins_1 [gNMuPtbins_1+1]  = { 5., 10., 20., 30., 40., 50., 60., 70.}; // ., 50., 65., 80.};

static const int gNMuPt2bins_0  = 7;
static const double gMuPt2bins_0 [gNMuPt2bins_0+1]  = {10., 20., 30., 40., 50., 60., 70., 100.};

static const int gNMuPt2bins_1 = 7;
static const double gMuPt2bins_1[gNMuPt2bins_1+1] = { 5., 10., 20., 30., 40., 50., 60., 70.}; // ., 50., 65., 80.};

static const int    gNMuEtabins = 1;
static const double gMuEtabins[gNMuEtabins+1] = {-2.4, 2.4};

// Electron Binning //////////////////////////////////////////////////////////////
static const int gNElPtbins_0  = 6;
static const double gElPtbins_0 [gNElPtbins_0+1]  = {20., 30., 40., 50., 60., 70., 100.};

static const int gNElPtbins_1  = 6;
static const double gElPtbins_1 [gNElPtbins_1+1]  = {10., 20., 30., 40., 50., 60., 70.}; // ., 50., 65., 80.};

static const int gNElPt2bins_0 = 7;
static const double gElPt2bins_0[gNElPt2bins_0+1] = {10., 20., 30., 40., 50., 60., 70., 100.};

static const int gNElPt2bins_1 = 6;
static const double gElPt2bins_1[gNElPt2bins_1+1] = {10., 20., 30., 40., 50., 60., 70.}; // ., 50., 65., 80.};

static const int    gNElEtabins = 1;
static const double gElEtabins[gNElEtabins+1] = {-2.4, 2.4};
//////////////////////////////////////////////////////////////////////////////////


//____________________________________________________________________________
MuonPlotter::MuonPlotter(){
// Default constructor, no samples are set
}

//____________________________________________________________________________
MuonPlotter::MuonPlotter(TString outputdir){
// Explicit constructor with output directory
	setOutputDir(outputdir);
}

//____________________________________________________________________________
MuonPlotter::MuonPlotter(TString outputdir, TString outputfile){
// Explicit constructor with output directory and output file
	setOutputDir(outputdir);
	setOutputFile(outputfile);
}

//____________________________________________________________________________
MuonPlotter::~MuonPlotter(){
	if(fOutputFile->IsOpen()) fOutputFile->Close();
}

//____________________________________________________________________________
void MuonPlotter::init(TString filename){
	if(fVerbose > 0) cout << "------------------------------------" << endl;
	if(fVerbose > 0) cout << "Initializing MuonPlotter ... " << endl;
	if(fVerbose > 0){
		cout << " ... using ";
		if(fSelectionSwitch == 0)      cout << "UCSD/UCSB/FNAL selection." << endl;
		else if(fSelectionSwitch == 1) cout << "UFlorida selection." << endl;
		else cout << "???" << endl;
		if(fChargeSwitch == 1) cout << " ... using opposite-sign charge selection" << endl;
	}
	Util::SetStyle();
	loadSamples(filename);
	readVarNames("anavarnames.dat");
	fOutputFileName = fOutputDir + "Yields.root";

	// fLumiNorm = fSamples[0].lumi; // Normalize everything to this lumi in /pb
	fLumiNorm = 1000; // Normalize everything to this lumi in /pb
	fBinWidthScale = 10.; // Normalize Y axis to this binwidth
	fDoCounting = false; // Disable counters by default
	fMinPt1 = 20.;
	fMinPt2 = 10.;

	// Prevent root from adding histograms to current file
	TH1::AddDirectory(kFALSE);

	fMCBG.push_back(TTbar);
	fMCBG.push_back(WJets);
	fMCBG.push_back(ZJets);
	fMCBG.push_back(AstarJets);
	fMCBG.push_back(VVJets);
	fMCBG.push_back(QCD15);
	fMCBG.push_back(QCD30);
	fMCBG.push_back(QCD80);
	fMCBG.push_back(QCD170);
	fMCBG.push_back(SSWWDPS);
	fMCBG.push_back(SSWWSPSPos);
	fMCBG.push_back(SSWWSPSNeg);

	fMCBGMuEnr.push_back(TTbar);
	fMCBGMuEnr.push_back(WJets);
	fMCBGMuEnr.push_back(ZJets);
	fMCBGMuEnr.push_back(AstarJets);
	fMCBGMuEnr.push_back(VVJets);
	fMCBGMuEnr.push_back(InclMu);
	fMCBGMuEnr.push_back(SSWWDPS);
	fMCBGMuEnr.push_back(SSWWSPSPos);
	fMCBGMuEnr.push_back(SSWWSPSNeg);

	fMCBGSig = fMCBG;
	fMCBGSig.push_back(LM0);

	fMCBGMuEnrSig = fMCBGMuEnr;
	fMCBGMuEnrSig.push_back(LM0);

	fMuData.push_back(MuA);
	fMuData.push_back(MuB);
	fEGData.push_back(EGA);
	fEGData.push_back(EGB);
	fJMData.push_back(JMA);
	fJMData.push_back(JMB);
	fJMData.push_back(MultiJet);

	fAllSamples.push_back(MuA);
	fAllSamples.push_back(MuB);
	fAllSamples.push_back(EGA);
	fAllSamples.push_back(EGB);
	fAllSamples.push_back(JMA);
	fAllSamples.push_back(JMB);
	fAllSamples.push_back(MultiJet);
	fAllSamples.push_back(TTbar);
	fAllSamples.push_back(WJets);
	fAllSamples.push_back(ZJets);
	fAllSamples.push_back(AstarJets);
	fAllSamples.push_back(VVJets);
	fAllSamples.push_back(QCD15);
	fAllSamples.push_back(QCD30);
	fAllSamples.push_back(QCD80);
	fAllSamples.push_back(QCD170);
	fAllSamples.push_back(SSWWDPS);
	fAllSamples.push_back(SSWWSPSPos);
	fAllSamples.push_back(SSWWSPSNeg);
	fAllSamples.push_back(LM0);
	fAllSamples.push_back(InclMu);

	bookHistos();
	bookRatioHistos();
}

//____________________________________________________________________________
void MuonPlotter::loadSamples(const char* filename){
	char buffer[200];
	ifstream IN(filename);

	char ParName[100], StringValue[1000];
	float ParValue;

	if(fVerbose > 2) cout << "------------------------------------" << endl;
	if(fVerbose > 2) cout << "Sample File  " << filename << endl;
	int counter(0);

	while( IN.getline(buffer, 200, '\n') ){
		// ok = false;
		if (buffer[0] == '#') continue; // Skip lines commented with '#'
		if( !strcmp(buffer, "SAMPLE")){
			Sample s;

			IN.getline(buffer, 200, '\n');
			sscanf(buffer, "Name\t%s", StringValue);
			s.name = TString(StringValue);

			IN.getline(buffer, 200, '\n');
			sscanf(buffer, "SName\t%s", StringValue);
			s.sname = TString(StringValue);

			IN.getline(buffer, 200, '\n');
			sscanf(buffer, "File\t%s", StringValue);
			TFile *f = TFile::Open(StringValue);
			s.file = f;
			s.tree = (TTree*)f->Get("Analysis");
			if(s.tree == NULL){ cout << " Tree not found in file!" << endl; break; }

			IN.getline(buffer, 200, '\n');
			sscanf(buffer, "Lumi\t%f", &ParValue);
			s.lumi = ParValue;

			IN.getline(buffer, 200, '\n');
			sscanf(buffer, "IsData\t%f", &ParValue);
			s.isdata = (bool)ParValue;

			IN.getline(buffer, 200, '\n');
			sscanf(buffer, "Color\t%f", &ParValue);
			s.color = ParValue;

			if(fVerbose > 2){
				cout << " ---- " << endl;
				cout << "  New sample added: " << s.name << endl;
				cout << "   Sample no.  " << counter << endl;
				cout << "   Short name: " << s.sname << endl;
				cout << "   File:       " << (s.file)->GetName() << endl;
				cout << "   Events:     " << s.tree->GetEntries() << endl;
				cout << "   Lumi:       " << s.lumi << endl;
				cout << "   Color:      " << s.color << endl;
				cout << "   IsData:     " << s.isdata << endl;
			}
			
			for(gRegion r = region_begin; r < gNREGIONS; r=gRegion(r+1)){
				Region *R = &(s.region[r]);
				if(r == Signal){
					R->name  = "Signal";
					R->sname = "Sig";
				}
				if(r == SigSup){
					R->name  = "Signal Suppressed";
					R->sname = "SigSup";
				}
				if(r == ZDecay){
					R->name  = "Z Decay";
					R->sname = "ZDecay";
				}
				for(gChannel c = channels_begin; c < gNCHANNELS; c=gChannel(c+1)){
					Channel *C;
					if(c == Muon){
						C = &R->mm;
						C->name  = "Mu/Mu";
						C->sname = "MM";
					}
					if(c == Electron){
						C = &R->ee;
						C->name  = "El/El";
						C->sname = "EE";
					}
					if(c == EMu){
						C = &R->em;
						C->name  = "El/Mu";
						C->sname = "EM";
					}
				}
			}
			
			fSampleMap[s.sname] = counter;
			fSamples.push_back(s);
			counter++;
		}
	}
	if(fVerbose > 2) cout << "------------------------------------" << endl;
}

//____________________________________________________________________________
const int     MuonPlotter::getNPtBins (gChannel chan){
	if(chan == Muon || chan == EMu){
		if(fSelectionSwitch == 0) return gNMuPtbins_0;
		if(fSelectionSwitch == 1) return gNMuPtbins_1;		
	}
	if(chan == Electron){
		if(fSelectionSwitch == 0) return gNElPtbins_0;
		if(fSelectionSwitch == 1) return gNElPtbins_1;		
	}
}
const double *MuonPlotter::getPtBins  (gChannel chan){
	if(chan == Muon || chan == EMu){
		if(fSelectionSwitch == 0) return gMuPtbins_0;
		if(fSelectionSwitch == 1) return gMuPtbins_1;
	}
	if(chan == Electron){
		if(fSelectionSwitch == 0) return gElPtbins_0;
		if(fSelectionSwitch == 1) return gElPtbins_1;
	}
}
const int     MuonPlotter::getNPt2Bins(gChannel chan){
	if(chan == Muon || chan == EMu){
		if(fSelectionSwitch == 0) return gNMuPt2bins_0;
		if(fSelectionSwitch == 1) return gNMuPt2bins_1;
	}
	if(chan == Electron){
		if(fSelectionSwitch == 0) return gNElPt2bins_0;
		if(fSelectionSwitch == 1) return gNElPt2bins_1;
	}
}
const double *MuonPlotter::getPt2Bins (gChannel chan){
	if(chan == Muon || chan == EMu){
		if(fSelectionSwitch == 0) return gMuPt2bins_0;
		if(fSelectionSwitch == 1) return gMuPt2bins_1;
	}
	if(chan == Electron){
		if(fSelectionSwitch == 0) return gElPt2bins_0;
		if(fSelectionSwitch == 1) return gElPt2bins_1;
	}	
}
const int     MuonPlotter::getNEtaBins(gChannel chan){
	if(chan == Muon || chan == EMu) return gNMuEtabins;
	if(chan == Electron)            return gNElEtabins;
}
const double *MuonPlotter::getEtaBins (gChannel chan){
	if(chan == Muon || chan == EMu) return gMuEtabins;
	if(chan == Electron)            return gElEtabins;
}

//____________________________________________________________________________
void MuonPlotter::doAnalysis(){
	// if(readHistos(fOutputFileName) != 0) return;

	// sandBox();
	
	fLumiNorm = 35.;
	// printYields(Muon);
	// printYields(Electron);
	// printYields(EMu);

	// printYields(Muon,     fLumiNorm);
	// printYields(Electron, fLumiNorm);
	// printYields(EMu,      fLumiNorm);

	// printYieldsShort();

	// makeIntPrediction(fOutputDir + "Yields.txt");
	makeMuIsolationPlot();
	
	// makeMufEffPlots(false);
	// makeElfEffPlots(false);
	// makeMupEffPlots(false);
	// makeElpEffPlots(false);

	// makeMufRatioPlots(false);
	// makeMupRatioPlots(false);
	// makeElfRatioPlots(false);
	// makeElpRatioPlots(false);
	// makeMCClosurePlots(fMCBG);
	
	// makeDataClosurePlots();
	// makeNT012Plots(EMu, fJMData, &MuonPlotter::isOSLLElMuEventTRG, "JMData_OS_");
	// makeNT012Plots(EMu, fMCBG,   &MuonPlotter::isOSLLElMuEventTRG, "MC_OS_");
	// makeNT012Plots(Muon, fJMData,    &MuonPlotter::isSSLLMuEventHTControlTRG, "JMData");
	// makeNT012Plots(Muon, fMCBGMuEnr, &MuonPlotter::isSSLLMuEventHTControlTRG, "MC");

	// makeMuIsoVsPtPlot(fMCBG, 0, &MuonPlotter::isSigSupMuEvent, &MuonPlotter::isLooseMuon, fMuData, 0, &MuonPlotter::isSigSupMuEventTRG, &MuonPlotter::isLooseMuon, "IsoVsPt_SigSuppressed", false);
}

//____________________________________________________________________________
void MuonPlotter::sandBox(){
	fOutputSubDir = "sandbox";
	// TEfficiency *eff = getEfficiency(&fSamples[LM0], Muon, SigSup, 1);
	// printObject(eff, "Test", "Test", "");
	
	vector<int> samples = fMCBG;
	
	fLumiNorm = 10.;
	TGraphAsymmErrors *f_mc = combineMCEfficiencies(samples, Muon, SigSup, true);
	// TEfficiency *f_mc = getEfficiency(samples, Muon, SigSup, true);
	printObject(f_mc, "TestEff", "TestEff", "AP");
	TH1D *h_fallmc = fillMuRatioPt(samples, SigSup, true);
	printObject(h_fallmc, "TestRat", "TestRat", "PE1");

	fOutputSubDir = "";
}

//____________________________________________________________________________
void MuonPlotter::doLoop(){
	vector<int> samples = fAllSamples;
	int step = 5000;
	fDoCounting = true;
	TString eventfilename = fOutputDir + "InterestingEvents.txt";
	fOUTSTREAM.open(eventfilename.Data(), ios::trunc);

	// Sample loop
	for(size_t i = 0; i < samples.size(); ++i){
		int index = samples[i];
		Sample *S = &fSamples[index];
		fCurrentSample = gSample(index);

		TTree *tree = S->tree;
		const bool isdata = S->isdata;


		// Stuff to execute for each sample BEFORE looping on the events
		initCounters(index);

		// Event loop
		tree->ResetBranchAddresses();
		Init(tree);
		if (fChain == 0) return;
		Long64_t nentries = fChain->GetEntriesFast();
		Long64_t nbytes = 0, nb = 0;
		for (Long64_t jentry=0; jentry<nentries;jentry++) {
			step = nentries/20;
			if( step < 200 ) step = 200;
			if( step > 10000 ) step = 10000;
			printProgress(jentry, nentries, S->name, step);

			Long64_t ientry = LoadTree(jentry);
			if (ientry < 0) break;
			nb = fChain->GetEntry(jentry);   nbytes += nb;

			fCounters[fCurrentSample][Muon]    .fill("All events");
			fCounters[fCurrentSample][EMu]     .fill("All events");
			fCounters[fCurrentSample][Electron].fill("All events");
		
			// Select mutually exclusive runs for Jet and MultiJet datasets
			if(!isGoodRun(index)) continue;
		
			fCounters[fCurrentSample][Muon]    .fill(" ... is good run");
			fCounters[fCurrentSample][EMu]     .fill(" ... is good run");
			fCounters[fCurrentSample][Electron].fill(" ... is good run");
		
			fillYields(S);

		}
		cout << endl;
		
		// Stuff to execute for each sample AFTER looping on the events
		storeNumbers(S, Muon);
		storeNumbers(S, Electron);
		storeNumbers(S, EMu);

	}
	fOUTSTREAM.close();
	writeHistos();
	printCutFlows(fOutputDir + "CutFlow.txt");
	fDoCounting = false;
}

//____________________________________________________________________________
void MuonPlotter::makeMCClosurePlots(vector<int> samples){
	// Fill the ratios
	fLumiNorm = 1000.;
	
	fillMuElRatios(samples);
	fixPRatios(); // puts p ratios to 1 whereever they are 0 (no entries)

	makeSSMuMuPredictionPlots(samples, false);
	makeSSElElPredictionPlots(samples, false);	
	makeSSElMuPredictionPlots(samples, false);

	makeNT012Plots(samples, Muon,     Signal);
	makeNT012Plots(samples, Electron, Signal);
	makeNT012Plots(samples, EMu,      Signal);
}

//____________________________________________________________________________
void MuonPlotter::makeDataClosurePlots(){
	// Fill the ratios
	// fOutputSubDir = "DataClosurePlots";
	fLumiNorm = 35.;

	vector<int> samples;
	if(fSelectionSwitch == 0) samples = fMuData;
	if(fSelectionSwitch == 1) samples = fJMData;

	const int    munptbins             = 1;
	const double muptbins[munptbins+1]   = {5., 1000.};
	const int    munetabins            = 1;
	const double muetabins[munetabins+1] = {-2.4, 2.4};
	const int    elnptbins             = 1;
	const double elptbins[elnptbins+1]   = {10., 1000.};
	const int    elnetabins            = 1;
	const double eletabins[elnetabins+1] = {-2.4, 2.4};

	TH2D *MufRatio2D  = new TH2D("TempMufRatio",    "Ratio of tight to loose Muons vs Pt vs Eta", munptbins, muptbins, munetabins, muetabins);
	TH1D *MufRatiopt  = new TH1D("TempMufRatioPt",  "Ratio of tight to loose Muons vs Pt",        munptbins, muptbins);
	TH1D *MufRatioeta = new TH1D("TempMufRatioEta", "Ratio of tight to loose Muons vs Eta",       munetabins, muetabins);
	TH2D *MupRatio2D  = new TH2D("TempMupRatio",    "Ratio of tight to loose Muons vs Pt vs Eta", munptbins, muptbins, munetabins, muetabins);
	TH1D *MupRatiopt  = new TH1D("TempMupRatioPt",  "Ratio of tight to loose Muons vs Pt",        munptbins, muptbins);
	TH1D *MupRatioeta = new TH1D("TempMupRatioEta", "Ratio of tight to loose Muons vs Eta",       munetabins, muetabins);
	TH2D *ElfRatio2D  = new TH2D("TempElfRatio",    "Ratio of tight to loose Electrons vs Pt vs Eta", elnptbins,  elptbins, elnetabins, eletabins);
	TH1D *ElfRatiopt  = new TH1D("TempElfRatioPt",  "Ratio of tight to loose Electrons vs Pt",        elnptbins,  elptbins);
	TH1D *ElfRatioeta = new TH1D("TempElfRatioEta", "Ratio of tight to loose Electrons vs Eta",       elnetabins, eletabins);
	TH2D *ElpRatio2D  = new TH2D("TempElpRatio",    "Ratio of tight to loose Electrons vs Pt vs Eta", elnptbins,  elptbins, elnetabins, eletabins);
	TH1D *ElpRatiopt  = new TH1D("TempElpRatioPt",  "Ratio of tight to loose Electrons vs Pt",        elnptbins,  elptbins);
	TH1D *ElpRatioeta = new TH1D("TempElpRatioEta", "Ratio of tight to loose Electrons vs Eta",       elnetabins, eletabins);

	produceRatio(Muon,     fJMData, 0, &MuonPlotter::isSigSupMuEventTRG, &MuonPlotter::isLooseMuon,     MufRatio2D, MufRatiopt, MufRatioeta); // f Ratio
	produceRatio(Muon,     fMuData, 0, &MuonPlotter::isZMuMuEventTRG,    &MuonPlotter::isLooseMuon,     MupRatio2D, MupRatiopt, MupRatioeta); // p Ratio
	produceRatio(Electron, fJMData, 0, &MuonPlotter::isSigSupElEventTRG, &MuonPlotter::isLooseElectron, ElfRatio2D, ElfRatiopt, ElfRatioeta); // f Ratio
	produceRatio(Electron, fEGData, 0, &MuonPlotter::isZElElEventTRG,    &MuonPlotter::isLooseElectron, ElpRatio2D, ElpRatiopt, ElpRatioeta); // p Ratio
	delete MufRatiopt, MufRatioeta, MupRatiopt, MupRatioeta; // don't need those
	delete ElfRatiopt, ElfRatioeta, ElpRatiopt, ElpRatioeta; //

	// Prediction: /////////////////////////////////////////////////////////////////////
	// Observation: ////////////////////////////////////////////////////////////////////
	TH1D *H_npppred = new TH1D("Npppred", "Predicted N_sig in Pt1 bins",       munptbins,  muptbins);
	TH1D *H_npfpred = new TH1D("Npfpred", "Predicted N_pf in Pt1 bins",        munptbins,  muptbins);
	TH1D *H_nfppred = new TH1D("Nfppred", "Predicted N_fp in Pt1 bins",        munptbins,  muptbins);
	TH1D *H_nffpred = new TH1D("Nffpred", "Predicted N_ff in Pt1 bins",        munptbins,  muptbins);
	TH1D *H_nFpred  = new TH1D("NFpred",  "Total predicted fakes in Pt1 bins", munptbins,  muptbins);
	TH1D *H_nt2obs  = new TH1D("Nt2obs",  "Observed Nt2 in Pt1 bins",          munptbins,  muptbins);

	H_npppred->Sumw2();
	H_npfpred->Sumw2();
	H_nfppred->Sumw2();
    H_nffpred->Sumw2();
    H_nFpred ->Sumw2();
    H_nt2obs ->Sumw2();

	bool output = true;

	for(size_t i = 0; i < samples.size(); ++i){
		int index = samples[i];
		fCurrentSample = gSample(index);
		// float scale = fLumiNorm/fSamples[samples[i]].lumi;

		TH2D *h_nt20_temp = new TH2D(Form("Nt20_temp_%s", fSamples[index].sname.Data()), "Nt20_temp", munptbins, muptbins, elnptbins, elptbins);
		TH2D *h_nt10_temp = new TH2D(Form("Nt10_temp_%s", fSamples[index].sname.Data()), "Nt10_temp", munptbins, muptbins, elnptbins, elptbins);
		TH2D *h_nt01_temp = new TH2D(Form("Nt01_temp_%s", fSamples[index].sname.Data()), "Nt01_temp", munptbins, muptbins, elnptbins, elptbins);
		TH2D *h_nt00_temp = new TH2D(Form("Nt00_temp_%s", fSamples[index].sname.Data()), "Nt00_temp", munptbins, muptbins, elnptbins, elptbins);
		// TH2D *h_nt20_temp = new TH2D(Form("Nt20_temp_%s", fSamples[index].sname.Data()), "Nt20_temp", munptbins, muptbins, munptbins, muptbins);
		// TH2D *h_nt10_temp = new TH2D(Form("Nt10_temp_%s", fSamples[index].sname.Data()), "Nt10_temp", munptbins, muptbins, munptbins, muptbins);
		// TH2D *h_nt01_temp = new TH2D(Form("Nt01_temp_%s", fSamples[index].sname.Data()), "Nt01_temp", munptbins, muptbins, munptbins, muptbins);
		// TH2D *h_nt00_temp = new TH2D(Form("Nt00_temp_%s", fSamples[index].sname.Data()), "Nt00_temp", munptbins, muptbins, munptbins, muptbins);
		
		// Event loop
		TTree *tree = fSamples[index].tree;
		tree->ResetBranchAddresses();
		Init(tree);
		if (fChain == 0) return;
		Long64_t nentries = fChain->GetEntriesFast();
		Long64_t nbytes = 0, nb = 0;
		for (Long64_t jentry=0; jentry<nentries;jentry++) {
			printProgress(jentry, nentries, fSamples[index].name);

			Long64_t ientry = LoadTree(jentry);
			if (ientry < 0) break;
			nb = fChain->GetEntry(jentry);   nbytes += nb;

			int mu(-1), el(-1);
			if(!isOSLLElMuEventTRG(mu, el)) continue; // OS channel
			// if(!isSSLLElMuEventInvMETTRG(mu, el)) continue; // signal window with inverted MET cut (MET < 30 GeV)
			
			if( isTightMuon(mu) &&  isTightElectron(el)) H_nt2obs->Fill(MuPt[mu]);
			
			if( isTightMuon(mu) &&  isTightElectron(el)) h_nt20_temp->Fill(MuPt[mu], ElPt[el]);
			if( isTightMuon(mu) && !isTightElectron(el)) h_nt10_temp->Fill(MuPt[mu], ElPt[el]);
			if(!isTightMuon(mu) &&  isTightElectron(el)) h_nt01_temp->Fill(MuPt[mu], ElPt[el]);
			if(!isTightMuon(mu) && !isTightElectron(el)) h_nt00_temp->Fill(MuPt[mu], ElPt[el]);

			// int mu1(-1), mu2(-1);
			// // if(!isSSLLMuEventInvMETTRG(mu1, mu2)) continue; // signal window with inverted MET cut (MET < 30 GeV)
			// if(!isSSLLMuEventHTControlTRG(mu1, mu2)) continue; // signal window in HT control region (250 GeV < HT < 350 GeV)
			// // if(!isSigSupSSMuMuEvent(mu1, mu2)) continue;
			// 
			// if( isTightMuon(mu1) &&  isTightMuon(mu2)) H_nt2obs->Fill(MuPt[mu1]);
			// 
			// if( isTightMuon(mu1) &&  isTightMuon(mu2)) h_nt20_temp->Fill(MuPt[mu1], MuPt[mu2]);
			// if( isTightMuon(mu1) && !isTightMuon(mu2)) h_nt10_temp->Fill(MuPt[mu1], MuPt[mu2]);
			// if(!isTightMuon(mu1) &&  isTightMuon(mu2)) h_nt10_temp->Fill(MuPt[mu2], MuPt[mu1]);
			// if(!isTightMuon(mu1) && !isTightMuon(mu2)) h_nt00_temp->Fill(MuPt[mu1], MuPt[mu2]);
		}
		vector<TH1D*> prediction = ElMuFPPrediction(MufRatio2D, MupRatio2D, ElfRatio2D, ElpRatio2D, h_nt20_temp, h_nt10_temp, h_nt01_temp, h_nt00_temp, output);
		// vector<TH1D*> prediction = MuMuFPPrediction(MufRatio2D, MupRatio2D, h_nt20_temp, h_nt10_temp, h_nt00_temp, output);
		H_npppred->Add(prediction[0]);
		H_npfpred->Add(prediction[1]);
		H_nfppred->Add(prediction[2]);
		H_nffpred->Add(prediction[3]);
		
		delete h_nt20_temp;
		delete h_nt10_temp;
		delete h_nt01_temp;
		delete h_nt00_temp;
	}

	// Output
	H_nFpred->Add(H_nfppred);
	H_nFpred->Add(H_npfpred);
	H_nFpred->Add(H_nffpred);
	H_nFpred->SetXTitle(convertVarName("MuPt[0]"));
	H_nFpred->SetYTitle(Form("Events / %2.0f GeV", 95.));
	// H_nFpred->SetYTitle(Form("Events / %2.0f GeV", fBinWidthScale));
	H_nt2obs->SetXTitle(convertVarName("MuPt[0]"));
	H_nt2obs->SetYTitle(Form("Events / %2.0f GeV", 95.));
	// H_nt2obs->SetYTitle(Form("Events / %2.0f GeV", fBinWidthScale));

	// Normalize to binwidth
	// H_npppred = normHistBW(H_npppred, fBinWidthScale);
	// H_nfppred = normHistBW(H_nfppred, fBinWidthScale);
	// H_nffpred = normHistBW(H_nffpred, fBinWidthScale);
	// H_nFpred  = normHistBW(H_nFpred,  fBinWidthScale);
	// H_nt2obs  = normHistBW(H_nt2obs,  fBinWidthScale);

	H_nt2obs->SetFillColor(kBlue);
	H_nt2obs->SetLineColor(kBlue);
	H_nt2obs->SetFillStyle(3004);
	H_nt2obs->SetLineWidth(2);


	int ppcolor = 8;
	H_npppred->SetFillColor(  ppcolor);
	H_npppred->SetMarkerColor(ppcolor);
	H_npppred->SetLineColor(  ppcolor);
	H_npppred->SetMarkerStyle(20);
	H_npppred->SetLineWidth(2);

	int fpcolor = 2;
	H_nfppred->SetFillColor(  fpcolor);
	H_nfppred->SetMarkerColor(fpcolor);
	H_nfppred->SetLineColor(  fpcolor);
	H_nfppred->SetMarkerStyle(20);
	H_nfppred->SetLineWidth(2);

	int pfcolor = 51;
	H_npfpred->SetFillColor(  pfcolor);
	H_npfpred->SetMarkerColor(pfcolor);
	H_npfpred->SetLineColor(  pfcolor);
	H_npfpred->SetMarkerStyle(20);
	H_npfpred->SetLineWidth(2);

	int ffcolor = 13;
	H_nffpred->SetFillColor(  ffcolor);
	H_nffpred->SetMarkerColor(ffcolor);
	H_nffpred->SetLineColor(  ffcolor);
	H_nffpred->SetMarkerStyle(20);
	H_nffpred->SetLineWidth(2);

	plotOverlay5H(H_nt2obs, "N_{ t2}", H_npppred, "N_{ pp}" , H_nfppred, "N_{ f p}" , H_npfpred, "N_{ p f}", H_nffpred, "N_{ f f}");

	H_nFpred->SetMinimum(0.);
	H_nt2obs->SetMinimum(0.);
	H_npppred->SetMinimum(0.);

	H_npppred->SetMarkerColor(kRed);
	H_nFpred->SetMarkerColor(kRed);
	H_nFpred->SetMarkerStyle(20);
	// H_npppred->SetMaximum(14.);
	// H_npppred->SetMinimum(0.);

	plotPredOverlay2HWithRatio(H_nFpred, "Predicted Fakes", H_nt2obs,  "Obs. N_{t2}", false, false);
	// fOutputSubDir = "";
}

//____________________________________________________________________________
void MuonPlotter::makeNT012Plots(vector<int> mcsamples, gChannel chan, gRegion reg){
	TString name;
	if(chan == Muon)     name = "MuMu";
	if(chan == Electron) name = "ElEl";
	if(chan == EMu)      name = "ElMu";
	
	fOutputSubDir = name + "Predictions";
	
	THStack *hnt2_stack = new THStack(Form("%s_nt2_stack", name.Data()), "Observed Nt2");
	THStack *hnt1_stack = new THStack(Form("%s_nt1_stack", name.Data()), "Observed Nt1");
	THStack *hnt0_stack = new THStack(Form("%s_nt0_stack", name.Data()), "Observed Nt0");
	const unsigned int nmcsamples = mcsamples.size();
	TH1D* hnt2[nmcsamples];
	TH1D* hnt1[nmcsamples];
	TH1D* hnt0[nmcsamples];

	for(size_t i = 0; i < mcsamples.size(); ++i){
		Sample S = fSamples[mcsamples[i]];
		float scale = fLumiNorm / S.lumi;
		Channel *cha;
		if(chan == Muon)     cha = &S.region[reg].mm;
		if(chan == Electron) cha = &S.region[reg].ee;
		if(chan == EMu)      cha = &S.region[reg].em;
		hnt2[i] = (TH1D*)(cha->nt20_pt->ProjectionX())->Clone();
		hnt1[i] = (TH1D*)(cha->nt10_pt->ProjectionX())->Clone();
		hnt0[i] = (TH1D*)(cha->nt00_pt->ProjectionX())->Clone();
		hnt2[i]->Scale(scale);
		hnt1[i]->Scale(scale);
		hnt0[i]->Scale(scale);
		
		hnt2[i]->SetFillColor(S.color);
		hnt1[i]->SetFillColor(S.color);
		hnt0[i]->SetFillColor(S.color);

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
		leg->AddEntry(hnt2[i], fSamples[index].sname.Data(), "f");
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

//____________________________________________________________________________
void MuonPlotter::makeNT012Plots(gChannel chan, vector<int> mcsamples, bool(MuonPlotter::*eventSelector)(int&, int&), TString tag){
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
			hnt20[i]->SetFillColor(fSamples[index].color);
			hnt10[i] = (TH1D*)list10->At(i);
			hnt10[i]->SetFillColor(fSamples[index].color);
			hnt01[i] = (TH1D*)list01->At(i);
			hnt01[i]->SetFillColor(fSamples[index].color);
			hnt00[i] = (TH1D*)list00->At(i);
			hnt00[i]->SetFillColor(fSamples[index].color);
		}
	}
	
	if(!read){
		TTree *tree = NULL;
		for(size_t i = 0; i < mcsamples.size(); ++i){
			int index = mcsamples[i];
			tree = fSamples[index].tree;
			hnt20[i] = new TH1D(Form("nt20_%s", fSamples[index].sname.Data()), "Observed Nt20", getNPt2Bins(Muon), getPt2Bins(Muon));
			hnt10[i] = new TH1D(Form("nt10_%s", fSamples[index].sname.Data()), "Observed Nt10", getNPt2Bins(Muon), getPt2Bins(Muon));
			hnt01[i] = new TH1D(Form("nt01_%s", fSamples[index].sname.Data()), "Observed Nt01", getNPt2Bins(Muon), getPt2Bins(Muon));
			hnt00[i] = new TH1D(Form("nt00_%s", fSamples[index].sname.Data()), "Observed Nt00", getNPt2Bins(Muon), getPt2Bins(Muon));
			hnt20[i]->SetFillColor(fSamples[index].color);
			hnt10[i]->SetFillColor(fSamples[index].color);
			hnt01[i]->SetFillColor(fSamples[index].color);
			hnt00[i]->SetFillColor(fSamples[index].color);
			hnt20[i]->Sumw2();
			hnt10[i]->Sumw2();
			hnt01[i]->Sumw2();
			hnt00[i]->Sumw2();
			float scale = fLumiNorm / fSamples[index].lumi;
			if(fSamples[index].isdata) scale = 1;
			tree->ResetBranchAddresses();
			Init(tree);
			if (fChain == 0) return;
			Long64_t nentries = fChain->GetEntriesFast();
			Long64_t nbytes = 0, nb = 0;
			for (Long64_t jentry=0; jentry<nentries;jentry++) {
				Long64_t ientry = LoadTree(jentry);
				if (ientry < 0) break;
				if(fVerbose > 1) printProgress(jentry, nentries, fSamples[index].name);
				nb = fChain->GetEntry(jentry);   nbytes += nb;

				int ind1(-1), ind2(-1);
				if((*this.*eventSelector)(ind1, ind2) == false) continue;

				if(chan == Muon){
					if( isTightMuon(ind1) &&  isTightMuon(ind2)) hnt20[i]->Fill(MuPt[ind1], scale);
					if( isTightMuon(ind1) && !isTightMuon(ind2)) hnt10[i]->Fill(MuPt[ind1], scale);
					if(!isTightMuon(ind1) &&  isTightMuon(ind2)) hnt10[i]->Fill(MuPt[ind2], scale);
					if(!isTightMuon(ind1) && !isTightMuon(ind2)) hnt00[i]->Fill(MuPt[ind1], scale);
				}
				if(chan == Electron){
					if( isTightElectron(ind1) &&  isTightElectron(ind2)) hnt20[i]->Fill(ElPt[ind1], scale);
					if( isTightElectron(ind1) && !isTightElectron(ind2)) hnt10[i]->Fill(ElPt[ind1], scale);
					if(!isTightElectron(ind1) &&  isTightElectron(ind2)) hnt10[i]->Fill(ElPt[ind2], scale);
					if(!isTightElectron(ind1) && !isTightElectron(ind2)) hnt00[i]->Fill(ElPt[ind1], scale);
				}
				if(chan == EMu){
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
		leg->AddEntry(hnt20[i], fSamples[index].sname.Data(), "f");
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
void MuonPlotter::makeMufRatioPlots(bool output){
	fOutputSubDir = "Ratios/Muon";
	fLumiNorm = 35;
	// TH1D *h_fdata  = fillRatioPt(fMuData, 0, &MuonPlotter::isSigSupMuEventTRG, &MuonPlotter::isLooseMuon);      // JetMET Dataset (Single Muon Selection)
	// TH1D *h_fttbar = fillRatioPt(TTbar,   0, &MuonPlotter::isGoodMuEvent,                  &MuonPlotter::isFakeTTbarMuon); // TTbarJets MC
	// TH1D *h_fallmc = fillRatioPt(fMCBG,   0, &MuonPlotter::isSigSupMuEvent,    &MuonPlotter::isLooseMuon);      // QCD MC


	TH1D *h_fdata1 = fillMuRatioPt(fJMData, SigSup, output);      // JetMET Dataset (Single Muon Selection)
	TH1D *h_fdata2 = fillMuRatioPt(fMuData, SigSup, output);
	TH1D *h_fallmc = fillMuRatioPt(fMCBG,   SigSup, output);      // QCD MC
	TH1D *h_fttbar = fillMuRatioPt(TTbar,   0, &MuonPlotter::isSigSupMuEvent, &MuonPlotter::isFakeTTbarMuon); // TTbarJets MC
	h_fdata1->SetName("MufRatioData");
	h_fdata2->SetName("MufRatioDataMu");
	h_fttbar->SetName("MufRatioTTbar");
	h_fallmc->SetName("MufRatioAllMC");

	setPlottingRange(h_fdata1, h_fttbar, h_fallmc);

	// h_fdata1->SetMinimum(0.);
	// h_fdata2->SetMinimum(0.);
	// h_fttbar->SetMinimum(0.);
	// h_fallmc->SetMinimum(0.);

	h_fdata1->SetMarkerColor(kBlack);
	h_fdata2->SetMarkerColor(kBlue);
	h_fttbar->SetMarkerColor(kBlue);
	h_fallmc->SetMarkerColor(kRed);

	h_fdata1->SetMarkerStyle(20);
	h_fdata2->SetMarkerStyle(20);
	h_fttbar->SetMarkerStyle(20);
	h_fallmc->SetMarkerStyle(20);

	// h_fdata1->SetMarkerSize(2);
	// h_fdata2->SetMarkerSize(2);
	// h_fttbar->SetMarkerSize(2);
	// h_fallmc->SetMarkerSize(2);

	h_fdata1->SetLineWidth(2);
	h_fdata2->SetLineWidth(2);
	h_fttbar->SetLineWidth(2);
	h_fallmc->SetLineWidth(2);

	h_fdata1->SetLineColor(kBlack);
	h_fdata2->SetLineColor(kBlue);
	h_fttbar->SetLineColor(kBlue);
	h_fallmc->SetLineColor(kRed);

	h_fdata1->SetFillColor(kBlack);
	h_fdata2->SetFillColor(kBlue);
	h_fttbar->SetFillColor(kBlue);
	h_fallmc->SetFillColor(kRed);

	plotRatioOverlay3H(h_fdata1, "Data (Jet, L = 21.7 pb^{-1})", h_fttbar, "t#bar{t} Fake GenMatch", h_fallmc, "QCD, t#bar{t}+jets, V+jets");
	// setPlottingRange(h_fdata1, h_fdata2);
	plotRatioOverlay2H(h_fdata1, "Data (Jet, L = 21.7 pb^{-1})", h_fdata2, "Data (Muon, L = 21.3 pb^{-1})");
	// plotRatioOverlay3H(h_fdata2, "Data (Mu, L = 21.3 pb^{-1})", h_fttbar, "t#bar{t} Fake GenMatch", h_fallmc, "QCD, t#bar{t}+jets, V+jets");

	vector<int> mcsamples = fMCBGSig;
	THStack *hntight_stack = new THStack("MufTightStack", "Stack of tight muons in sig. supp. selection");
	THStack *hnloose_stack = new THStack("MufLooseStack", "Stack of loose muons in sig. supp. selection");
	const unsigned int nmcsamples = mcsamples.size();
	TH1D *hntight[nmcsamples];
	TH1D *hnloose[nmcsamples];

	for(size_t i = 0; i < mcsamples.size(); ++i){
		Sample *S = &fSamples[mcsamples[i]];
		Channel *cha = &S->region[Signal].mm;
		hntight[i] = cha->fntight->ProjectionX();
		hnloose[i] = cha->fnloose->ProjectionX();
		hntight[i]->SetFillColor(S->color);
		hnloose[i]->SetFillColor(S->color);
		float scale = fLumiNorm / S->lumi;
		hntight[i]->Scale(scale);
		hnloose[i]->Scale(scale);
		hntight_stack->Add(hntight[i]);
		hnloose_stack->Add(hnloose[i]);
	}
	hntight_stack->Draw();
	hntight_stack->GetXaxis()->SetTitle(convertVarName("MuPt[0]"));
	hnloose_stack->Draw();
	hnloose_stack->GetXaxis()->SetTitle(convertVarName("MuPt[0]"));
	
	TCanvas *c_tight = new TCanvas("MufStackTight", "Tight Stack", 0, 0, 800, 600);
	TCanvas *c_loose = new TCanvas("MufStackLoose", "Loose Stack", 0, 0, 800, 600);

	// TLegend *leg = new TLegend(0.13,0.60,0.38,0.88);
	TLegend *leg = new TLegend(0.75,0.60,0.89,0.88);
	for(size_t i = 0; i < nmcsamples; ++i){
		int index = mcsamples[i];
		leg->AddEntry(hntight[i], fSamples[index].sname.Data(), "f");
	}
	leg->SetFillStyle(0);
	leg->SetTextFont(42);
	leg->SetBorderSize(0);

	c_tight->cd();
	gPad->SetLogy();
	hntight_stack->Draw("hist");
	leg->Draw();

	c_loose->cd();
	gPad->SetLogy();
	hnloose_stack->Draw("hist");
	leg->Draw();

	Util::PrintNoEPS(c_tight, "MufRatioTightStack", fOutputDir + fOutputSubDir, fOutputFile);
	Util::PrintNoEPS(c_loose, "MufRatioLooseStack", fOutputDir + fOutputSubDir, fOutputFile);	
	fOutputSubDir = "";
}
void MuonPlotter::makeElfRatioPlots(bool output){
	fOutputSubDir = "Ratios/Electron";
	fLumiNorm = 35.;
	TH1D *h_fdata1 = fillElRatioPt(fJMData, SigSup, output);      // JetMET Dataset (Single Muon Selection)
	TH1D *h_fdata2 = fillElRatioPt(fEGData, SigSup, output);
	TH1D *h_fallmc = fillElRatioPt(fMCBG,   SigSup, output);      // QCD MC
	h_fdata1->SetName("ElfRatioData");
	h_fdata2->SetName("ElfRatioDataEl");
	h_fallmc->SetName("ElfRatioAllMC");

	setPlottingRange(h_fdata1, h_fdata2, h_fallmc);

	// h_fdata1->SetMinimum(0.);
	// h_fdata2->SetMinimum(0.);
	// h_fallmc->SetMinimum(0.);

	h_fdata1->SetMarkerColor(kBlack);
	h_fdata2->SetMarkerColor(kBlue);
	h_fallmc->SetMarkerColor(kRed);

	h_fdata1->SetMarkerStyle(20);
	h_fdata2->SetMarkerStyle(20);
	h_fallmc->SetMarkerStyle(20);

	// h_fdata1->SetMarkerSize(2);
	// h_fdata2->SetMarkerSize(2);
	// h_fallmc->SetMarkerSize(2);

	h_fdata1->SetLineWidth(2);
	h_fdata2->SetLineWidth(2);
	h_fallmc->SetLineWidth(2);

	h_fdata1->SetLineColor(kBlack);
	h_fdata2->SetLineColor(kBlue);
	h_fallmc->SetLineColor(kRed);

	h_fdata1->SetFillColor(kBlack);
	h_fdata2->SetFillColor(kBlue);
	h_fallmc->SetFillColor(kRed);

	plotRatioOverlay2H(h_fdata1, "Data (Jet, L = 35 pb^{-1})", h_fallmc, "QCD, t#bar{t}+jets, V+jets");
	// setPlottingRange(h_fdata1, h_fdata2);
	plotRatioOverlay2H(h_fdata1, "Data (Jet, L = 35 pb^{-1})", h_fdata2, "Data (EGamma, L = 35 pb^{-1})");

	vector<int> mcsamples = fMCBGSig;
	THStack *hntight_stack = new THStack("ElfTightStack", "Stack of tight electrons in sig. supp. selection");
	THStack *hnloose_stack = new THStack("ElfLooseStack", "Stack of loose electrons in sig. supp. selection");
	const unsigned int nmcsamples = mcsamples.size();
	TH1D *hntight[nmcsamples];
	TH1D *hnloose[nmcsamples];

	for(size_t i = 0; i < mcsamples.size(); ++i){
		Sample S = fSamples[mcsamples[i]];
		Channel *cha = &S.region[Signal].ee;
		hntight[i] = cha->fntight->ProjectionX();
		hnloose[i] = cha->fnloose->ProjectionX();
		hntight[i]->SetFillColor(S.color);
		hnloose[i]->SetFillColor(S.color);
		float scale = fLumiNorm / S.lumi;
		hntight[i]->Scale(scale);
		hnloose[i]->Scale(scale);
		hntight_stack->Add(hntight[i]);
		hnloose_stack->Add(hnloose[i]);
	}
	hntight_stack->Draw();
	hntight_stack->GetXaxis()->SetTitle(convertVarName("ElPt[0]"));
	hnloose_stack->Draw();
	hnloose_stack->GetXaxis()->SetTitle(convertVarName("ElPt[0]"));
	
	TCanvas *c_tight = new TCanvas("ElfStackTight", "Tight Stack", 0, 0, 800, 600);
	TCanvas *c_loose = new TCanvas("ElfStackLoose", "Loose Stack", 0, 0, 800, 600);

	// TLegend *leg = new TLegend(0.13,0.60,0.38,0.88);
	TLegend *leg = new TLegend(0.75,0.60,0.89,0.88);
	for(size_t i = 0; i < nmcsamples; ++i){
		Sample S = fSamples[mcsamples[i]];
		leg->AddEntry(hntight[i], S.sname.Data(), "f");
	}
	leg->SetFillStyle(0);
	leg->SetTextFont(42);
	leg->SetBorderSize(0);

	c_tight->cd();
	// gPad->SetLogy();
	hntight_stack->Draw("hist");
	leg->Draw();

	c_loose->cd();
	// gPad->SetLogy();
	hnloose_stack->Draw("hist");
	leg->Draw();

	Util::PrintNoEPS(c_tight, "ElfRatioTightStack", fOutputDir + fOutputSubDir, fOutputFile);
	Util::PrintNoEPS(c_loose, "ElfRatioLooseStack", fOutputDir + fOutputSubDir, fOutputFile);	
	fOutputSubDir = "";
}
void MuonPlotter::makeMupRatioPlots(bool output){
	fOutputSubDir = "Ratios/Muon";
	fLumiNorm = 35.;
	// TH1D *h_pdata  = fillRatioPt(fMuData, 0, &MuonPlotter::isZMuMuEventTRG, &MuonPlotter::isLooseMuon); // Mu Dataset (Di Muon Selection)
	// TH1D *h_pttbar = fillRatioPt(TTbar,   0, &MuonPlotter::isGoodMuEvent, &MuonPlotter::isPromptTTbarMuon); // TTbar
	// TH1D *h_pallmc = fillRatioPt(fMCBG,   0, &MuonPlotter::isZMuMuEvent,    &MuonPlotter::isLooseMuon); // all MC
	TH1D *h_pdata  = fillMuRatioPt(fMuData, ZDecay, output);
	TH1D *h_pallmc = fillMuRatioPt(fMCBG,   ZDecay, output);
	TH1D *h_pttbar = fillMuRatioPt(TTbar,   0, &MuonPlotter::isGoodMuEvent, &MuonPlotter::isPromptTTbarMuon); // TTbar
	h_pdata ->SetName("MupRatioData");
	h_pttbar->SetName("MupRatioTTbar");
	h_pallmc->SetName("MupRatioAllMC");

	setPlottingRange(h_pdata, h_pttbar, h_pallmc);

	h_pdata ->Draw("goff");
	h_pttbar->Draw("goff");
	h_pallmc->Draw("goff");

	// h_pdata ->SetMinimum(0.);
	// h_pttbar->SetMinimum(0.);
	// h_pallmc->SetMinimum(0.);
	// h_pttbar->SetMaximum(1.3);
	// h_pallmc->SetMaximum(1.3);

	h_pdata ->SetLineWidth(2);
	h_pttbar->SetLineWidth(2);
	h_pallmc->SetLineWidth(2);

	h_pdata ->SetMarkerColor(kBlack);
	h_pttbar->SetMarkerColor(kBlue);
	h_pallmc->SetMarkerColor(kRed);

	h_pdata ->SetMarkerStyle(20);
	h_pttbar->SetMarkerStyle(20);
	h_pallmc->SetMarkerStyle(20);

	// h_pdata ->SetMarkerSize(2);
	// h_pttbar->SetMarkerSize(2);
	// h_pallmc->SetMarkerSize(2);

	h_pdata ->SetLineColor(kBlack);
	h_pttbar->SetLineColor(kBlue);
	h_pallmc->SetLineColor(kRed);

	h_pdata ->SetFillColor(kBlack);
	h_pttbar->SetFillColor(kBlue);
	h_pallmc->SetFillColor(kRed);

	h_pdata ->SetDrawOption("E1");
	h_pttbar->SetDrawOption("E1");
	h_pallmc->SetDrawOption("E1");

	// plotRatioOverlay3H(h_pdata, "Data (Jet, L = 21.7 pb^{-1})", h_pttbar, "t#bar{t} Prompt GenMatch", h_pallmc, "QCD, t#bar{t}+jets, V+jets");
	plotRatioOverlay3H(h_pdata, "Data (Mu, L = 21.3 pb^{-1})", h_pttbar, "t#bar{t} Prompt GenMatch", h_pallmc, "QCD, t#bar{t}+jets, V+jets");

	vector<int> mcsamples = fMCBGSig;
	THStack *hntight_stack = new THStack("MupTightStack", "Stack of tight muons in Z decay selection");
	THStack *hnloose_stack = new THStack("MupLooseStack", "Stack of loose muons in Z decay selection");
	const unsigned int nmcsamples = mcsamples.size();
	TH1D *hntight[nmcsamples];
	TH1D *hnloose[nmcsamples];

	for(size_t i = 0; i < mcsamples.size(); ++i){
		Sample S = fSamples[mcsamples[i]];
		Channel *cha = &S.region[Signal].mm;
		hntight[i] = cha->pntight->ProjectionX();
		hnloose[i] = cha->pnloose->ProjectionX();
		hntight[i]->SetFillColor(S.color);
		hnloose[i]->SetFillColor(S.color);
		float scale = fLumiNorm / S.lumi;
		hntight[i]->Scale(scale);
		hnloose[i]->Scale(scale);
		hntight_stack->Add(hntight[i]);
		hnloose_stack->Add(hnloose[i]);
	}
	hntight_stack->Draw();
	hntight_stack->GetXaxis()->SetTitle(convertVarName("MuPt[0]"));
	hnloose_stack->Draw();
	hnloose_stack->GetXaxis()->SetTitle(convertVarName("MuPt[0]"));
	
	TCanvas *c_tight = new TCanvas("MupStackTight", "Tight Stack", 0, 0, 800, 600);
	TCanvas *c_loose = new TCanvas("MupStackLoose", "Loose Stack", 0, 0, 800, 600);

	// TLegend *leg = new TLegend(0.13,0.60,0.38,0.88);
	TLegend *leg = new TLegend(0.75,0.60,0.89,0.88);
	for(size_t i = 0; i < nmcsamples; ++i){
		Sample S = fSamples[mcsamples[i]];
		leg->AddEntry(hntight[i], S.sname.Data(), "f");
	}
	leg->SetFillStyle(0);
	leg->SetTextFont(42);
	leg->SetBorderSize(0);

	c_tight->cd();
	// gPad->SetLogy();
	hntight_stack->Draw("hist");
	leg->Draw();

	c_loose->cd();
	// gPad->SetLogy();
	hnloose_stack->Draw("hist");
	leg->Draw();

	Util::PrintNoEPS(c_tight, "MupRatioTightStack", fOutputDir + fOutputSubDir, fOutputFile);
	Util::PrintNoEPS(c_loose, "MupRatioLooseStack", fOutputDir + fOutputSubDir, fOutputFile);	
	fOutputSubDir = "";
}
void MuonPlotter::makeElpRatioPlots(bool output){
	fOutputSubDir = "Ratios/Electron";
	fLumiNorm = 35.;
	TH1D *h_pdata  = fillElRatioPt(fEGData, ZDecay, output);
	TH1D *h_pallmc = fillElRatioPt(fMCBG,   ZDecay, output);
	h_pdata ->SetName("ElpRatioData");
	h_pallmc->SetName("ElpRatioAllMC");

	setPlottingRange(h_pdata, h_pallmc);

	h_pdata ->Draw("goff");
	h_pallmc->Draw("goff");

	// h_pdata ->SetMinimum(0.);
	// h_pallmc->SetMinimum(0.);

	h_pdata ->SetLineWidth(2);
	h_pallmc->SetLineWidth(2);

	h_pdata ->SetMarkerColor(kBlack);
	h_pallmc->SetMarkerColor(kRed);

	h_pdata ->SetMarkerStyle(20);
	h_pallmc->SetMarkerStyle(20);

	h_pdata ->SetLineColor(kBlack);
	h_pallmc->SetLineColor(kRed);

	h_pdata ->SetFillColor(kBlack);
	h_pallmc->SetFillColor(kRed);

	h_pdata ->SetDrawOption("E1");
	h_pallmc->SetDrawOption("E1");

	setPlottingRange(h_pdata, h_pallmc);
	plotRatioOverlay2H(h_pdata, "Data (EG, L = 35 pb^{-1})", h_pallmc, "QCD, t#bar{t}+jets, V+jets");

	vector<int> mcsamples = fMCBGSig;
	THStack *hntight_stack = new THStack("ElpTightStack", "Stack of tight electrons in Z decay selection");
	THStack *hnloose_stack = new THStack("ElpLooseStack", "Stack of loose electrons in Z decay selection");
	const unsigned int nmcsamples = mcsamples.size();
	TH1D *hntight[nmcsamples];
	TH1D *hnloose[nmcsamples];

	for(size_t i = 0; i < mcsamples.size(); ++i){
		Sample S = fSamples[mcsamples[i]];
		Channel *cha = &S.region[Signal].ee;
		hntight[i] = cha->pntight->ProjectionX();
		hnloose[i] = cha->pnloose->ProjectionX();
		hntight[i]->SetFillColor(S.color);
		hnloose[i]->SetFillColor(S.color);
		float scale = fLumiNorm / S.lumi;
		hntight[i]->Scale(scale);
		hnloose[i]->Scale(scale);
		hntight_stack->Add(hntight[i]);
		hnloose_stack->Add(hnloose[i]);
	}
	hntight_stack->Draw();
	hntight_stack->GetXaxis()->SetTitle(convertVarName("ElPt[0]"));
	hnloose_stack->Draw();
	hnloose_stack->GetXaxis()->SetTitle(convertVarName("ElPt[0]"));

	TCanvas *c_tight = new TCanvas("ElpStackTight", "Tight Stack", 0, 0, 800, 600);
	TCanvas *c_loose = new TCanvas("ElpStackLoose", "Loose Stack", 0, 0, 800, 600);

	// TLegend *leg = new TLegend(0.13,0.60,0.38,0.88);
	TLegend *leg = new TLegend(0.75,0.60,0.89,0.88);
	for(size_t i = 0; i < nmcsamples; ++i){
		Sample S = fSamples[mcsamples[i]];
		leg->AddEntry(hntight[i], S.sname.Data(), "f");
	}
	leg->SetFillStyle(0);
	leg->SetTextFont(42);
	leg->SetBorderSize(0);

	c_tight->cd();
	// gPad->SetLogy();
	hntight_stack->Draw("hist");
	leg->Draw();

	c_loose->cd();
	// gPad->SetLogy();
	hnloose_stack->Draw("hist");
	leg->Draw();

	Util::PrintNoEPS(c_tight, "ElpRatioTightStack", fOutputDir + fOutputSubDir, fOutputFile);
	Util::PrintNoEPS(c_loose, "ElpRatioLooseStack", fOutputDir + fOutputSubDir, fOutputFile);	
	fOutputSubDir = "";
}

//____________________________________________________________________________
void MuonPlotter::makeMufEffPlots(bool output){
	fOutputSubDir = "Ratios/Muon/";
	bool stacks = true;
	fLumiNorm = 35;
	// vector<int> mcsamples = fMCBGMuEnr;
	vector<int> mcsamples = fMCBG;
	TEfficiency       *e_fdata1 = mergeDataEfficiencies(fJMData,   Muon, SigSup, output, gStatOpt, gStatBetaAlpha, gStatBetaBeta);
	TEfficiency       *e_fdata2 = mergeDataEfficiencies(fMuData,   Muon, SigSup, output, gStatOpt, gStatBetaAlpha, gStatBetaBeta);
	TGraphAsymmErrors *g_fallmc = combineMCEfficiencies(mcsamples, Muon, SigSup, output, gStatOpt, gStatBetaAlpha, gStatBetaBeta);
	e_fdata1->SetName("MufEffData");
	e_fdata2->SetName("MufEffDataMu");
	g_fallmc->SetName("MufEffAllMC");

	e_fdata1->SetMarkerColor(kBlack);
	e_fdata2->SetMarkerColor(kBlue);
	g_fallmc->SetMarkerColor(kRed);

	e_fdata1->SetMarkerStyle(20);
	e_fdata2->SetMarkerStyle(20);
	g_fallmc->SetMarkerStyle(20);

	e_fdata1->SetLineWidth(2);
	e_fdata2->SetLineWidth(2);
	g_fallmc->SetLineWidth(2);

	e_fdata1->SetLineColor(kBlack);
	e_fdata2->SetLineColor(kBlue);
	g_fallmc->SetLineColor(kRed);

	e_fdata1->SetFillColor(kBlack);
	e_fdata2->SetFillColor(kBlue);
	g_fallmc->SetFillColor(kRed);

	plotEffOverlayEG(e_fdata1, "Data (Jet, L = 21.7 pb^{-1})", g_fallmc, "QCD, t#bar{t}+jets, V+jets");
	plotEffOverlayEE(e_fdata1, "Data (Jet, L = 21.7 pb^{-1})", e_fdata2, "Data (Muon, L = 21.3 pb^{-1})");
	// plotRatioOverlay3H(h_fdata1, "Data (Jet, L = 21.7 pb^{-1})", h_fttbar, "t#bar{t} Fake GenMatch", h_fallmc, "QCD, t#bar{t}+jets, V+jets");

	////// STACKS ///////////////////////////////////////////////////////////////////////////////////////////////////
	if(stacks){
		THStack *hntight_stack = new THStack("MufTightStack", "Stack of tight muons in sig. supp. selection");
		THStack *hnloose_stack = new THStack("MufLooseStack", "Stack of loose muons in sig. supp. selection");
		const unsigned int nmcsamples = mcsamples.size();
		TH1D *hntight[nmcsamples];
		TH1D *hnloose[nmcsamples];
	
		for(size_t i = 0; i < mcsamples.size(); ++i){
			Sample *S = &fSamples[mcsamples[i]];
			Channel *cha = &S->region[Signal].mm;
			hntight[i] = cha->fntight->ProjectionX();
			hnloose[i] = cha->fnloose->ProjectionX();
			hntight[i]->SetFillColor(S->color);
			hnloose[i]->SetFillColor(S->color);
			float scale = fLumiNorm / S->lumi;
			hntight[i]->Scale(scale);
			hnloose[i]->Scale(scale);
			hntight_stack->Add(hntight[i]);
			hnloose_stack->Add(hnloose[i]);
		}
		hntight_stack->Draw();
		hntight_stack->GetXaxis()->SetTitle(convertVarName("MuPt[0]"));
		hnloose_stack->Draw();
		hnloose_stack->GetXaxis()->SetTitle(convertVarName("MuPt[0]"));
	
		TCanvas *c_tight = new TCanvas("MufStackTight", "Tight Stack", 0, 0, 800, 600);
		TCanvas *c_loose = new TCanvas("MufStackLoose", "Loose Stack", 0, 0, 800, 600);
	
		// TLegend *leg = new TLegend(0.13,0.60,0.38,0.88);
		TLegend *leg = new TLegend(0.75,0.60,0.89,0.88);
		for(size_t i = 0; i < nmcsamples; ++i){
			int index = mcsamples[i];
			leg->AddEntry(hntight[i], fSamples[index].sname.Data(), "f");
		}
		leg->SetFillStyle(0);
		leg->SetTextFont(42);
		leg->SetBorderSize(0);
	
		c_tight->cd();
		gPad->SetLogy();
		hntight_stack->Draw("hist");
		leg->Draw();
	
		c_loose->cd();
		gPad->SetLogy();
		hnloose_stack->Draw("hist");
		leg->Draw();
	
		Util::PrintNoEPS(c_tight, "MufRatioTightStack", fOutputDir + fOutputSubDir, fOutputFile);
		Util::PrintNoEPS(c_loose, "MufRatioLooseStack", fOutputDir + fOutputSubDir, fOutputFile);		
	}
	fOutputSubDir = "";
}
void MuonPlotter::makeElfEffPlots(bool output){
	fOutputSubDir = "Ratios/Electron/";
	bool stacks = true;
	fLumiNorm = 35;
	// vector<int> mcsamples = fMCBGMuEnr;
	vector<int> mcsamples = fMCBG;
	TEfficiency       *e_fdata1 = mergeDataEfficiencies(fJMData,   Electron, SigSup, output, gStatOpt, gStatBetaAlpha, gStatBetaBeta);
	TEfficiency       *e_fdata2 = mergeDataEfficiencies(fEGData,   Electron, SigSup, output, gStatOpt, gStatBetaAlpha, gStatBetaBeta);
	TGraphAsymmErrors *g_fallmc = combineMCEfficiencies(mcsamples, Electron, SigSup, output, gStatOpt, gStatBetaAlpha, gStatBetaBeta);
	e_fdata1->SetName("ElfEffData");
	e_fdata2->SetName("ElfEffDataMu");
	g_fallmc->SetName("ElfEffAllMC");

	e_fdata1->SetMarkerColor(kBlack);
	e_fdata2->SetMarkerColor(kBlue);
	g_fallmc->SetMarkerColor(kRed);

	e_fdata1->SetMarkerStyle(20);
	e_fdata2->SetMarkerStyle(20);
	g_fallmc->SetMarkerStyle(20);

	e_fdata1->SetLineWidth(2);
	e_fdata2->SetLineWidth(2);
	g_fallmc->SetLineWidth(2);

	e_fdata1->SetLineColor(kBlack);
	e_fdata2->SetLineColor(kBlue);
	g_fallmc->SetLineColor(kRed);

	e_fdata1->SetFillColor(kBlack);
	e_fdata2->SetFillColor(kBlue);
	g_fallmc->SetFillColor(kRed);

	plotEffOverlayEG(e_fdata1, "Data (Jet, L = 21.7 pb^{-1})", g_fallmc, "QCD, t#bar{t}+jets, V+jets");
	plotEffOverlayEE(e_fdata1, "Data (Jet, L = 21.7 pb^{-1})", e_fdata2, "Data (EG, L = 21.3 pb^{-1})");
	// plotRatioOverlay3H(h_fdata1, "Data (Jet, L = 21.7 pb^{-1})", h_fttbar, "t#bar{t} Fake GenMatch", h_fallmc, "QCD, t#bar{t}+jets, V+jets");

	////// STACKS ///////////////////////////////////////////////////////////////////////////////////////////////////
	if(stacks){
		THStack *hntight_stack = new THStack("ElfTightStack", "Stack of tight electrons in sig. supp. selection");
		THStack *hnloose_stack = new THStack("ElfLooseStack", "Stack of loose electrons in sig. supp. selection");
		const unsigned int nmcsamples = mcsamples.size();
		TH1D *hntight[nmcsamples];
		TH1D *hnloose[nmcsamples];

		for(size_t i = 0; i < mcsamples.size(); ++i){
			Sample S = fSamples[mcsamples[i]];
			Channel *cha = &S.region[Signal].ee;
			hntight[i] = cha->fntight->ProjectionX();
			hnloose[i] = cha->fnloose->ProjectionX();
			hntight[i]->SetFillColor(S.color);
			hnloose[i]->SetFillColor(S.color);
			float scale = fLumiNorm / S.lumi;
			hntight[i]->Scale(scale);
			hnloose[i]->Scale(scale);
			hntight_stack->Add(hntight[i]);
			hnloose_stack->Add(hnloose[i]);
		}
		hntight_stack->Draw();
		hntight_stack->GetXaxis()->SetTitle(convertVarName("ElPt[0]"));
		hnloose_stack->Draw();
		hnloose_stack->GetXaxis()->SetTitle(convertVarName("ElPt[0]"));

		TCanvas *c_tight = new TCanvas("ElfStackTight", "Tight Stack", 0, 0, 800, 600);
		TCanvas *c_loose = new TCanvas("ElfStackLoose", "Loose Stack", 0, 0, 800, 600);

		// TLegend *leg = new TLegend(0.13,0.60,0.38,0.88);
		TLegend *leg = new TLegend(0.75,0.60,0.89,0.88);
		for(size_t i = 0; i < nmcsamples; ++i){
			Sample S = fSamples[mcsamples[i]];
			leg->AddEntry(hntight[i], S.sname.Data(), "f");
		}
		leg->SetFillStyle(0);
		leg->SetTextFont(42);
		leg->SetBorderSize(0);

		c_tight->cd();
		gPad->SetLogy();
		hntight_stack->Draw("hist");
		leg->Draw();

		c_loose->cd();
		gPad->SetLogy();
		hnloose_stack->Draw("hist");
		leg->Draw();

		Util::PrintNoEPS(c_tight, "ElfRatioTightStack", fOutputDir + fOutputSubDir, fOutputFile);
		Util::PrintNoEPS(c_loose, "ElfRatioLooseStack", fOutputDir + fOutputSubDir, fOutputFile);		
	}
	fOutputSubDir = "";
}
void MuonPlotter::makeMupEffPlots(bool output){
	fOutputSubDir = "Ratios/Muon/";
	bool stacks = true;
	fLumiNorm = 35;
	// vector<int> mcsamples = fMCBGMuEnr;
	vector<int> mcsamples = fMCBG;
	TEfficiency       *e_pdata1 = mergeDataEfficiencies(fJMData,   Muon, ZDecay, output, gStatOpt, gStatBetaAlpha, gStatBetaBeta);
	TEfficiency       *e_pdata2 = mergeDataEfficiencies(fMuData,   Muon, ZDecay, output, gStatOpt, gStatBetaAlpha, gStatBetaBeta);
	TGraphAsymmErrors *g_pallmc = combineMCEfficiencies(mcsamples, Muon, ZDecay, output, gStatOpt, gStatBetaAlpha, gStatBetaBeta);
	e_pdata1->SetName("MupEffData");
	e_pdata2->SetName("MupEffDataMu");
	g_pallmc->SetName("MupEffAllMC");

	e_pdata1->SetMarkerColor(kBlack);
	e_pdata2->SetMarkerColor(kBlue);
	g_pallmc->SetMarkerColor(kRed);

	e_pdata1->SetMarkerStyle(20);
	e_pdata2->SetMarkerStyle(20);
	g_pallmc->SetMarkerStyle(20);

	e_pdata1->SetLineWidth(2);
	e_pdata2->SetLineWidth(2);
	g_pallmc->SetLineWidth(2);

	e_pdata1->SetLineColor(kBlack);
	e_pdata2->SetLineColor(kBlue);
	g_pallmc->SetLineColor(kRed);

	e_pdata1->SetFillColor(kBlack);
	e_pdata2->SetFillColor(kBlue);
	g_pallmc->SetFillColor(kRed);

	plotEffOverlayEG(e_pdata1, "Data (Jet, L = 21.7 pb^{-1})", g_pallmc, "QCD, t#bar{t}+jets, V+jets");
	plotEffOverlayEE(e_pdata1, "Data (Jet, L = 21.7 pb^{-1})", e_pdata2, "Data (Muon, L = 21.3 pb^{-1})");
	// plotRatioOverlay3H(h_fdata1, "Data (Jet, L = 21.7 pb^{-1})", h_fttbar, "t#bar{t} Fake GenMatch", h_fallmc, "QCD, t#bar{t}+jets, V+jets");
	
	if(stacks){
		vector<int> mcsamples = fMCBGSig;
		THStack *hntight_stack = new THStack("MupTightStack", "Stack of tight muons in Z decay selection");
		THStack *hnloose_stack = new THStack("MupLooseStack", "Stack of loose muons in Z decay selection");
		const unsigned int nmcsamples = mcsamples.size();
		TH1D *hntight[nmcsamples];
		TH1D *hnloose[nmcsamples];

		for(size_t i = 0; i < mcsamples.size(); ++i){
			Sample S = fSamples[mcsamples[i]];
			Channel *cha = &S.region[Signal].mm;
			hntight[i] = cha->pntight->ProjectionX();
			hnloose[i] = cha->pnloose->ProjectionX();
			hntight[i]->SetFillColor(S.color);
			hnloose[i]->SetFillColor(S.color);
			float scale = fLumiNorm / S.lumi;
			hntight[i]->Scale(scale);
			hnloose[i]->Scale(scale);
			hntight_stack->Add(hntight[i]);
			hnloose_stack->Add(hnloose[i]);
		}
		hntight_stack->Draw();
		hntight_stack->GetXaxis()->SetTitle(convertVarName("MuPt[0]"));
		hnloose_stack->Draw();
		hnloose_stack->GetXaxis()->SetTitle(convertVarName("MuPt[0]"));

		TCanvas *c_tight = new TCanvas("MupStackTight", "Tight Stack", 0, 0, 800, 600);
		TCanvas *c_loose = new TCanvas("MupStackLoose", "Loose Stack", 0, 0, 800, 600);

		// TLegend *leg = new TLegend(0.13,0.60,0.38,0.88);
		TLegend *leg = new TLegend(0.75,0.60,0.89,0.88);
		for(size_t i = 0; i < nmcsamples; ++i){
			Sample S = fSamples[mcsamples[i]];
			leg->AddEntry(hntight[i], S.sname.Data(), "f");
		}
		leg->SetFillStyle(0);
		leg->SetTextFont(42);
		leg->SetBorderSize(0);

		c_tight->cd();
		// gPad->SetLogy();
		hntight_stack->Draw("hist");
		leg->Draw();

		c_loose->cd();
		// gPad->SetLogy();
		hnloose_stack->Draw("hist");
		leg->Draw();

		Util::PrintNoEPS(c_tight, "MupRatioTightStack", fOutputDir + fOutputSubDir, fOutputFile);
		Util::PrintNoEPS(c_loose, "MupRatioLooseStack", fOutputDir + fOutputSubDir, fOutputFile);	
		
	}
	fOutputSubDir = "";
}
void MuonPlotter::makeElpEffPlots(bool output){
	fOutputSubDir = "Ratios/Electron/";
	bool stacks = true;
	fLumiNorm = 35;
	// vector<int> mcsamples = fMCBGMuEnr;
	vector<int> mcsamples = fMCBG;
	TEfficiency       *e_pdata1 = mergeDataEfficiencies(fJMData,   Electron, ZDecay, output, gStatOpt, gStatBetaAlpha, gStatBetaBeta);
	TEfficiency       *e_pdata2 = mergeDataEfficiencies(fEGData,   Electron, ZDecay, output, gStatOpt, gStatBetaAlpha, gStatBetaBeta);
	TGraphAsymmErrors *g_pallmc = combineMCEfficiencies(mcsamples, Electron, ZDecay, output, gStatOpt, gStatBetaAlpha, gStatBetaBeta);
	e_pdata1->SetName("ElpEffData");
	e_pdata2->SetName("ElpEffDataMu");
	g_pallmc->SetName("ElpEffAllMC");

	e_pdata1->SetMarkerColor(kBlack);
	e_pdata2->SetMarkerColor(kBlue);
	g_pallmc->SetMarkerColor(kRed);

	e_pdata1->SetMarkerStyle(20);
	e_pdata2->SetMarkerStyle(20);
	g_pallmc->SetMarkerStyle(20);

	e_pdata1->SetLineWidth(2);
	e_pdata2->SetLineWidth(2);
	g_pallmc->SetLineWidth(2);

	e_pdata1->SetLineColor(kBlack);
	e_pdata2->SetLineColor(kBlue);
	g_pallmc->SetLineColor(kRed);

	e_pdata1->SetFillColor(kBlack);
	e_pdata2->SetFillColor(kBlue);
	g_pallmc->SetFillColor(kRed);

	plotEffOverlayEG(e_pdata1, "Data (Jet, L = 21.7 pb^{-1})", g_pallmc, "QCD, t#bar{t}+jets, V+jets");
	plotEffOverlayEE(e_pdata1, "Data (Jet, L = 21.7 pb^{-1})", e_pdata2, "Data (EG, L = 21.3 pb^{-1})");
	// plotRatioOverlay3H(h_fdata1, "Data (Jet, L = 21.7 pb^{-1})", h_fttbar, "t#bar{t} Fake GenMatch", h_fallmc, "QCD, t#bar{t}+jets, V+jets");

	////// STACKS ///////////////////////////////////////////////////////////////////////////////////////////////////
	if(stacks){
		THStack *hntight_stack = new THStack("ElpTightStack", "Stack of tight electrons in Z decay selection");
		THStack *hnloose_stack = new THStack("ElpLooseStack", "Stack of loose electrons in Z decay selection");
		const unsigned int nmcsamples = mcsamples.size();
		TH1D *hntight[nmcsamples];
		TH1D *hnloose[nmcsamples];

		for(size_t i = 0; i < mcsamples.size(); ++i){
			Sample S = fSamples[mcsamples[i]];
			Channel *cha = &S.region[Signal].ee;
			hntight[i] = cha->pntight->ProjectionX();
			hnloose[i] = cha->pnloose->ProjectionX();
			hntight[i]->SetFillColor(S.color);
			hnloose[i]->SetFillColor(S.color);
			float scale = fLumiNorm / S.lumi;
			hntight[i]->Scale(scale);
			hnloose[i]->Scale(scale);
			hntight_stack->Add(hntight[i]);
			hnloose_stack->Add(hnloose[i]);
		}
		hntight_stack->Draw();
		hntight_stack->GetXaxis()->SetTitle(convertVarName("ElPt[0]"));
		hnloose_stack->Draw();
		hnloose_stack->GetXaxis()->SetTitle(convertVarName("ElPt[0]"));

		TCanvas *c_tight = new TCanvas("ElpStackTight", "Tight Stack", 0, 0, 800, 600);
		TCanvas *c_loose = new TCanvas("ElpStackLoose", "Loose Stack", 0, 0, 800, 600);

		// TLegend *leg = new TLegend(0.13,0.60,0.38,0.88);
		TLegend *leg = new TLegend(0.75,0.60,0.89,0.88);
		for(size_t i = 0; i < nmcsamples; ++i){
			Sample S = fSamples[mcsamples[i]];
			leg->AddEntry(hntight[i], S.sname.Data(), "f");
		}
		leg->SetFillStyle(0);
		leg->SetTextFont(42);
		leg->SetBorderSize(0);

		c_tight->cd();
		// gPad->SetLogy();
		hntight_stack->Draw("hist");
		leg->Draw();

		c_loose->cd();
		// gPad->SetLogy();
		hnloose_stack->Draw("hist");
		leg->Draw();

		Util::PrintNoEPS(c_tight, "ElpRatioTightStack", fOutputDir + fOutputSubDir, fOutputFile);
		Util::PrintNoEPS(c_loose, "ElpRatioLooseStack", fOutputDir + fOutputSubDir, fOutputFile);	
	}
	fOutputSubDir = "";
}

//____________________________________________________________________________
void MuonPlotter::makeMuIsolationPlot(){
	const bool read = true;
	TString filename = "IsoHistos.root";
	const unsigned nselections = 2;
	const unsigned nbins = 20.;
	const float ptcut = 20.;

	TString sel_name[nselections];
	sel_name[0] = "Base";
	sel_name[1] = "SigSup";

	TH1D    *hiso     [nselections][gNSAMPLES];
	TH1D    *hiso_data[nselections];
	THStack *hiso_mc  [nselections];
	for(size_t i = 0; i < nselections; ++i){
		for(gSample j = sample_begin; j < gNSAMPLES; j=gSample(j+1)){
			hiso[i][j] = new TH1D(Form("hiso_%s_%s", sel_name[i].Data(), fSamples[j].sname.Data()), fSamples[j].sname.Data(), nbins, 0., 1.);
			hiso[i][j]->SetFillColor(fSamples[j].color);
			hiso[i][j]->Sumw2();
		}
		hiso_data[i] = new TH1D("IsoData_"  + sel_name[i], "Muon Isolation in Data for " + sel_name[i], nbins, 0., 1.);
		hiso_mc[i]   = new THStack("IsoMC_" + sel_name[i], "Muon Isolation in MC for "   + sel_name[i]);
	}
	// hiso_mc[i] = new TH1D(Form("iso_mc_%s", fSamples[index].sname.Data()), "Muon Isolation in MC", nbins, 0., 1.);

	TFile *file;
	if(!read){
		file = new TFile(fOutputDir + filename, "RECREATE");
		file->cd();
		fDoCounting = false;
		// Sample loop
		for(gSample i = sample_begin; i < gNSAMPLES; i=gSample(i+1)){
			Sample *S = &fSamples[i];
			fCurrentSample = gSample(i);

			TTree *tree = S->tree;

			// Event loop
			tree->ResetBranchAddresses();
			Init(tree);
			if (fChain == 0) return;
			Long64_t nentries = fChain->GetEntriesFast();
			Long64_t nbytes = 0, nb = 0;
			for (Long64_t jentry=0; jentry<nentries;jentry++) {
				printProgress(jentry, nentries, S->name);

				Long64_t ientry = LoadTree(jentry);
				if (ientry < 0) break;
				nb = fChain->GetEntry(jentry);   nbytes += nb;

				// Select mutually exclusive runs for Jet and MultiJet datasets
				if(!isGoodRun(i)) continue;

				////////////////////////////////////////////////////
				// MOST LOOSE SELECTION
				if(isGoodMuEvent() && isMuTriggeredEvent()){
					if(isLooseMuon(0) == false) continue;
					if(MuPt[0] < ptcut) continue;

					hiso[0][i]->Fill(MuIso[0]);
				}
				////////////////////////////////////////////////////
				// SIGNAL SUPPRESSED SELECTION
				if(isSigSupMuEventTRG()){
					if(isLooseMuon(0) == false) continue;
					if(MuPt[0] < ptcut) continue;

					hiso[1][i]->Fill(MuIso[0]);
				}
				////////////////////////////////////////////////////
			}
			cout << endl;

			// Write histos
			for(size_t j = 0; j < nselections; ++j){
				TDirectory* dir = Util::FindOrCreate(sel_name[j], file);
				dir->cd();
				hiso[j][i]->Write(hiso[j][i]->GetName(), TObject::kWriteDelete);
			}
		}
		file->Close();
	}
	else{
		file = TFile::Open(fOutputDir + filename, "READ");
		if(file == NULL){
			cout << "File " << filename << " does not exist!" << endl;
			return;
		}
		file->cd();
		
		for(size_t i = 0; i < nselections; ++i){
			for(gSample j = sample_begin; j < gNSAMPLES; j=gSample(j+1)){
				TString getname = sel_name[i] + "/" + hiso[0][j]->GetName(); // had bug when storing file, change this after rerunning!
				hiso[i][j] = (TH1D*)file->Get(getname);
				hiso[i][j]->SetFillColor(fSamples[j].color);
			}
		}
	}
	
	// vector<int> mcsamples = fMCBG;
	vector<int> mcsamples = fMCBGMuEnr;
	// vector<int> datasamples = fJMData;
	vector<int> datasamples = fMuData;
	
	for(size_t i = 0; i < nselections; ++i){
		hiso_data[i]->SetXTitle(convertVarName("MuIso[0]"));
		hiso_data[i]->SetLineWidth(3);
		hiso_data[i]->SetLineColor(kBlack);
		hiso_data[i]->SetMarkerStyle(8);
		hiso_data[i]->SetMarkerColor(kBlack);
		hiso_data[i]->SetMarkerSize(1.2);
	
		// Apply weights
		for(size_t j = 0; j < gNSAMPLES; ++j){
			float scale = fLumiNorm / fSamples[j].lumi;
			if(fSamples[j].isdata) continue;
			else hiso[i][j]->Scale(scale);
		}
	
		// Fill stacks
		for(size_t j = 0; j < datasamples.size(); ++j) hiso_data[i]->Add(hiso[i][datasamples[j]]);
		for(size_t j = 0; j < mcsamples.size();   ++j) hiso_mc[i]  ->Add(hiso[i][mcsamples  [j]]);
	
		double max1 = hiso_mc  [i]->GetMaximum();
		double max2 = hiso_data[i]->GetMaximum();
		double max = max1>max2?max1:max2;
	
		// hiso_mc  [i]->SetMinimum(10);
		// hiso_data[i]->SetMinimum(10);
		hiso_mc  [i]->SetMaximum(1.2*max);
	
		TCanvas *c_temp = new TCanvas("Iso" + sel_name[i], "Isolation in Data vs MC", 0, 0, 800, 600);
		c_temp->cd();
	
		TLegend *leg = new TLegend(0.13,0.60,0.38,0.88);
		// TLegend *leg = new TLegend(0.75,0.60,0.89,0.88);
		leg->AddEntry(hiso_data[i], "Data","p");
		for(size_t j = 0; j < mcsamples.size(); ++j) leg->AddEntry(hiso[i][mcsamples[j]], fSamples[mcsamples[j]].sname.Data(), "f");
		leg->SetFillStyle(0);
		leg->SetTextFont(42);
		leg->SetBorderSize(0);
		
		// gPad->SetLogy();
		hiso_mc[i]->Draw("hist");
		hiso_data[i]->DrawCopy("PE X0 same");
		leg->Draw();
	
		Util::PrintNoEPS(c_temp, "Iso" + sel_name[i], fOutputDir, NULL);	
	}	
}

//____________________________________________________________________________
void MuonPlotter::makeMuIsolationPlots(){
	// void MuonPlotter::makeIsoVsPtPlot(vector<int> samples1, int muon1, bool(MuonPlotter::*eventSelector1)(), bool(MuonPlotter::*muonSelector1)(int), vector<int> samples2, int muon2, bool(MuonPlotter::*eventSelector2)(), bool(MuonPlotter::*muonSelector2)(int), TString outputname, bool logy){
	const bool logy = false;
	const int nbins = 30;
	TH2D *h2_data = new TH2D("h2_data", "Isolation vs Pt for muons in data", nbins, 0., 1., getNPt2Bins(Muon), getPt2Bins(Muon));
	h2_data->SetXTitle(convertVarName("MuIso[0]"));
	h2_data->SetYTitle(convertVarName("MuPt[0]"));
	const unsigned int nmcsamples = fMCBG.size();
	TH2D *h2_mc[nmcsamples];
	for(size_t i = 0; i < nmcsamples; ++i){
		TString name = fSamples[i].sname;
		h2_mc[i] = new TH2D(Form("h2_mc_%s", name.Data()), Form("Isolation vs Pt for muons in %s", name.Data()), nbins, 0., 1., getNPt2Bins(Muon), getPt2Bins(Muon));
		h2_mc[i]  ->SetXTitle(convertVarName("MuIso[0]"));
		h2_mc[i]  ->SetYTitle(convertVarName("MuPt[0]"));
	}

	TTree *tree = NULL;
	for(size_t i = 0; i < fJMData.size(); ++i){
		int index = fJMData[i];

		tree = fSamples[index].tree;
		tree->ResetBranchAddresses();
		Init(tree);
		if (fChain == 0) return;
		Long64_t nentries = fChain->GetEntriesFast();
		Long64_t nbytes = 0, nb = 0;
		for (Long64_t jentry=0; jentry<nentries;jentry++) {
			Long64_t ientry = LoadTree(jentry);
			if (ientry < 0) break;
			nb = fChain->GetEntry(jentry);   nbytes += nb;

			if(isSigSupMuEventTRG() == false) continue;
			if(isLooseMuon(0)       == false) continue;

			h2_data->Fill(MuIso[0], MuPt[0]);
		}
	}

	for(size_t i = 0; i < nmcsamples; ++i){
		int index = fMCBG[i];

		tree = fSamples[index].tree;
		float scale = fLumiNorm / fSamples[index].lumi;
		tree->ResetBranchAddresses();
		Init(tree);
		if (fChain == 0) return;
		Long64_t nentries = fChain->GetEntriesFast();
		Long64_t nbytes = 0, nb = 0;
		for (Long64_t jentry=0; jentry<nentries;jentry++) {
			Long64_t ientry = LoadTree(jentry);
			if (ientry < 0) break;
			nb = fChain->GetEntry(jentry);   nbytes += nb;

			if(isSigSupMuEvent() == false) continue;
			if(isLooseMuon(0)    == false) continue;
			h2_mc[i]->Fill(MuIso[0], MuPt[0], scale);
		}
	}

	TLatex *lat = new TLatex();
	lat->SetNDC(kTRUE);
	lat->SetTextColor(kBlack);
	lat->SetTextSize(0.06);

	TCanvas *c_temp = new TCanvas("IsoVsPt", "Isolating in Pt bins", 0, 0, 1200, 800);
	c_temp->Divide(3,2);
	// TCanvas *c_temp = new TCanvas("IsoVsPt", "Isolating in Pt bins", 0, 0, 1200, 1200);
	// c_temp->Divide(3,3);
	// c_temp->cd(9);
	// if(logy) gPad->SetLogz(1);
	// h2_data->DrawCopy("colz");
	// lat->DrawLatex(0.11,0.92, fSamples[TTbar].sname);

	for(unsigned i = 1; i <= getNPt2Bins(Muon); ++i){
		c_temp->cd(i);
		gStyle->SetOptStat(1111);
		TH1D *h1 = h2_data->ProjectionX(Form("h2_datax_%d", i), i, i);
		THStack *hs2 = new THStack(Form("h2_mcx_stack_%d", i), "Stack of MC Histos");
		for(size_t j = 0; j < nmcsamples; ++j){
			int index = fMCBG[i];
			TH1D *h2_tmp = h2_mc[j]->ProjectionX(Form("h2_bgx_%d_%s", i, fSamples[index].sname.Data()), i, i);
			h2_tmp->SetFillColor(fSamples[index].color);
			hs2->Add(h2_tmp);
		}
		h1->SetXTitle(convertVarName("MuIso[0]"));
		h1->SetLineWidth(2);
		h1->SetLineColor(kBlack);
		h1->SetMarkerStyle(8);
		h1->SetMarkerColor(kBlack);
		// h1->SetFillColor(15);
		// h1->SetFillStyle(1001);

		// hs2->SetLineWidth(2);
		// hs2->SetLineColor(kBlue);
		// hs2->SetFillColor(kBlue);
		// hs2->SetFillStyle(3004);

		if(logy) gPad->SetLogy(1);
		gPad->SetFillStyle(0);
		h1->Sumw2();
		// hs2->Sumw2();

		// Scaling
		// if(h1->GetEntries() > 0 ) h1->Scale(1.0/h1->Integral());
		// if(hs2->GetEntries() > 0 ) hs2->Scale(1.0/h2->Integral());

		// setPlottingRange(h1, h2);

		// Determine plotting range
		double max1 = h1->GetMaximum();
		double max2 = hs2->GetMaximum();
		double max  = (max1>max2)?max1:max2;
		if(logy) max = 5*max;
		else max = 1.05*max;
		h1->SetMaximum(max);
		hs2->SetMaximum(max);
		if(!logy){
			h1->SetMinimum(0.0);
			hs2->SetMinimum(0.0);
		}

		TLegend *leg = new TLegend(0.45,0.75,0.65,0.88);
		if(i == 1){
			leg->AddEntry(h1,  "JM Data","f");
			leg->AddEntry(hs2, "All MC","f");
			leg->SetFillStyle(0);
			leg->SetTextFont(42);
			leg->SetBorderSize(0);
		}

		TH1D *h1_temp = (TH1D*)h1->DrawCopy("PE");
		hs2->Draw("same");
		h1_temp->SetName(Form("h1_%d",i));
		// hs2_temp->SetName(Form("hs2_%d",i));
		gPad->Update();
		TPaveStats *s1 = (TPaveStats*)h1_temp->GetListOfFunctions()->FindObject("stats");
		// TPaveStats *s2 = (TPaveStats*)hs2->GetListOfFunctions()->FindObject("stats");
		// s2->SetTextColor(kBlue); s2->SetLineColor(kBlue);
		// s2->SetY1NDC(s1->GetY1NDC() - (s1->GetY2NDC() - s1->GetY1NDC()));
		// s2->SetY2NDC(s1->GetY1NDC());

		if(i==1) leg->Draw();

		// double max1 = h1->GetYaxis()->GetXmax();
		// double max2 = h2->GetYaxis()->GetXmax();
		// double max  = (max1>max2)?max1:max2;
		double min1 = h1->GetYaxis()->GetXmin();
		double min2 = hs2->GetYaxis()->GetXmin();
		double min  = (min1<min2)?min1:min2;

		TLine *l1 = new TLine(0.15, min, 0.15, max);
		l1->SetLineColor(kRed);
		l1->SetLineWidth(2);
		l1->Draw();

		lat->SetTextColor(kBlack);
		lat->SetTextSize(0.05);
		lat->DrawLatex(0.11,0.92, Form("p_{T}(#mu) %3.0f - %3.0f GeV", getPt2Bins(Muon)[i-1], getPt2Bins(Muon)[i]));

		// int bin0 = h1->FindBin(0.00);
		// int bin15 = h1->FindBin(0.15);
		// int bin1 = h1->FindBin(1.00);
		// float f1 = h1->Integral(bin0, bin15) / h1->Integral(bin0, bin1);
		// float f2 = hs2->Integral(bin0, bin15) / hs2->Integral(bin0, bin1);
		// lat->SetTextSize(0.04);
		// lat->DrawLatex(0.55,0.905, Form("ratio = %4.2f", f1));
		// lat->SetTextColor(kBlue);
		// lat->DrawLatex(0.55,0.945, Form("ratio = %4.2f", f2));

		gPad->RedrawAxis();
	}
	Util::PrintNoEPS(c_temp, "IsoPlots", fOutputDir, fOutputFile);
}

//____________________________________________________________________________
void MuonPlotter::makeMuPtPlots(){
	const int nbins = 40;
	TH1D *h_prompt   = new TH1D("h_prompt",  "Pt for prompt Muons in WJets",                   nbins, 0., 100.);
	TH1D *h_fakew    = new TH1D("h_fakew",   "Pt for fake Muons in WJets events",              nbins, 0., 100.);
	TH1D *h_fakeqcd  = new TH1D("h_fakeqcd", "Pt for fake Muons in QCD events",                nbins, 0., 100.);
	TH1D *h_wtau     = new TH1D("h_wtau",    "Pt for fake Muons from tau in WJets events",     nbins, 0., 100.);
	TH1D *h_wnotau   = new TH1D("h_wnotau",  "Pt for fake Muons not from tau in WJets events", nbins, 0., 100.);
	TH1D *h_ttp      = new TH1D("h_ttp",     "Pt for prompt Muons in ttbar",                   nbins, 0., 100.);
	TH1D *h_ttp2     = new TH1D("h_ttp2",    "Pt for prompt Muons in ttbar",                   nbins, 0., 200.);
	TH1D *h_ttf      = new TH1D("h_ttf",     "Pt for non prompt Muons in ttbar",               nbins, 0., 100.);
	TH1D *h_ttftau   = new TH1D("h_ttftau",  "Pt for Muons from tau in ttbar",                 nbins, 0., 100.);
	TH1D *h_ttfnotau = new TH1D("h_ttfnotau","Pt for Muons not from tau in ttbar",             nbins, 0., 100.);
	TH1D *h_qcdb     = new TH1D("h_qcdb",   "Pt for muons from bottom hadrons in QCD",         nbins, 0., 100.);
	TH1D *h_qcdpik   = new TH1D("h_qcdpik", "Pt for muons from pions/kaons in QCD",            nbins, 0., 100.);
	TH1D *h_ttb      = new TH1D("h_ttb",    "Pt for muons from bottom hadrons in ttbar",       nbins, 0., 100.);
	TH1D *h_z        = new TH1D("h_z",      "Pt for muons from Z boson decays",                nbins, 0., 100.);

	h_fakeqcd->SetXTitle(convertVarName("MuPt[0]"));
	h_qcdb->SetXTitle(convertVarName("MuPt[0]"));
	h_prompt->SetXTitle(convertVarName("MuPt[0]"));
	h_ttp2->SetXTitle(convertVarName("MuPt[0]"));

	fSamples[6].tree->Project("h_prompt",  "MuPt[0]", "abs(MuGenMoID[0])==24");
	fSamples[5].tree->Project("h_z",       "MuPt[0]", "abs(MuGenMoID[0])==23");
	fSamples[1].tree->Project("h_fakew",   "MuPt[1]", "abs(MuGenMoID[1])!=24&&abs(MuGenGMoID[1])!=24");
	fSamples[0].tree->Project("h_fakeqcd", "MuPt[0]", "");
	fSamples[1].tree->Project("h_wtau",    "MuPt[1]", "abs(MuGenMoID[1])==15");
	fSamples[1].tree->Project("h_wnotau",  "MuPt[1]", "abs(MuGenMoID[1])!=24&&abs(MuGenMoID[1])!=15");
	fSamples[2].tree->Project("h_ttp",     "MuPt[0]", "abs(MuGenMoID[0])==24");
	fSamples[2].tree->Project("h_ttp2",    "MuPt[0]", "abs(MuGenMoID[0])==24");
	fSamples[2].tree->Project("h_ttf",     "MuPt[0]", "abs(MuGenMoID[0])!=24");
	fSamples[2].tree->Project("h_ttftau",  "MuPt[0]", "abs(MuGenMoID[0])==15");
	fSamples[2].tree->Project("h_ttfnotau","MuPt[0]", "abs(MuGenMoID[0])!=24&&abs(MuGenMoID[0])!=15");
	fSamples[2].tree->Project("h_ttb",     "MuPt[1]", "MuGenMoType[1]==15||MuGenMoType[1]==17||MuGenMoType[1]==21");
	fSamples[4].tree->Project("h_qcdb",    "MuPt[0]", "MuGenMoType[0]==15||MuGenMoType[0]==17||MuGenMoType[0]==21");
	fSamples[4].tree->Project("h_qcdpik",  "MuPt[0]", "MuGenMoType[0]==11||MuGenMoType[0]==12||MuGenMoType[0]==13");

	cout << h_fakeqcd->GetEntries() << " " << h_prompt->GetEntries() << " " << h_fakew->GetEntries() << endl;

	// h_ttp2->SetMinimum(0);
	printObject(h_ttp2, "TTbarPt", "Pt of prompt ttbar Muons", "hist");

	plotOverlay3H(h_fakeqcd, "QCD",         h_prompt, "W-jets",       h_ttp,      "ttbar",         true);
	plotOverlay3H(h_z,       "Z-jets",      h_prompt, "W-jets",       h_ttp,      "ttbar",         true);
	plotOverlay3H(h_fakeqcd, "QCD",         h_wnotau, "W-jets",       h_ttfnotau, "ttbar",         true);

	// plotOverlay3H(h_fakeqcd, "Fake in QCD", h_prompt, "Prompt",       h_fakew,    "Fake in WJets", true);
	// plotOverlay3H(h_fakeqcd, "QCD",         h_wtau,   "W: tau",       h_wnotau,   "W: No tau",     true);
	// plotOverlay3H(h_prompt,  "Prompt W",    h_ttp,    "Prompt ttbar", h_ttftau,   "ttbar tau",     true);
	// plotOverlay3H(h_fakeqcd, "QCD",         h_wnotau, "W: no tau",    h_ttfnotau, "ttbar: no tau", true);
	// plotOverlay3H(h_qcdb,    "QCD: b",      h_qcdpik, "QCD: #pi/K",   h_ttb,      "ttbar: b",      true);

	// plotOverlay2H(h_fakeqcd, "QCD", h_fakew, "WJets", false, 0.15);
}

//____________________________________________________________________________
void MuonPlotter::makeMuIsoVsPtPlot(int sample1, int muon1, TCut c1, int sample2, int muon2, TCut c2, TString outputname, bool logy){
	const int nbins = 20;
	TH2D *h2_sig = new TH2D("h2_sig", "Isolation vs Pt for fake muons in signal", nbins, 0., 1., getNPtBins(Muon), getPtBins(Muon));
	TH2D *h2_bg  = new TH2D("h2_bg", "Isolation vs Pt for muons in background",   nbins, 0., 1., getNPtBins(Muon), getPtBins(Muon));
	h2_sig->SetXTitle(convertVarName("MuIso[0]"));
	h2_sig->SetYTitle(convertVarName("MuPt[0]"));
	h2_bg->SetXTitle(convertVarName("MuIso[0]"));
	h2_bg->SetYTitle(convertVarName("MuPt[0]"));

	fSamples[sample1].tree->Project("h2_bg", Form("MuPt[%d]:MuIso[%d]", muon1, muon1), c1);
	fSamples[sample2].tree->Project("h2_sig", Form("MuPt[%d]:MuIso[%d]", muon2, muon2), c2);

	TLatex *lat = new TLatex();
	lat->SetNDC(kTRUE);
	lat->SetTextColor(kBlack);
	lat->SetTextSize(0.06);

	TCanvas *c_temp = new TCanvas("IsoVsPt", "Isolating in Pt bins", 0, 0, 1200, 800);
	c_temp->Divide(3,2);
	c_temp->cd(6);
	if(logy) gPad->SetLogz(1);
	h2_bg->DrawCopy("colz");
	lat->DrawLatex(0.11,0.92, fSamples[sample1].sname);

	for(size_t i = 1; i <= 5; ++i){
		c_temp->cd(i);
		gStyle->SetOptStat(1111);
		TH1D *h1 = h2_bg->ProjectionX("h2_bgx",i, i);	
		TH1D *h2 = h2_sig->ProjectionX("h2_sigx",i, i);
		h1->SetXTitle(convertVarName("MuIso[0]"));
		h1->SetLineWidth(2);
		h1->SetFillColor(15);
		h1->SetFillStyle(1001);
		h2->SetLineWidth(2);
		h2->SetLineColor(kBlue);
		h2->SetFillColor(kBlue);
		h2->SetFillStyle(3004);

		if(logy) gPad->SetLogy(1);
		gPad->SetFillStyle(0);
		h1->Sumw2();
		h2->Sumw2();

		// Scaling
		if(h1->GetEntries() > 0 ) h1->Scale(1.0/h1->Integral());
		if(h2->GetEntries() > 0 ) h2->Scale(1.0/h2->Integral());

		// Determine plotting range
		double max1 = h1->GetMaximum();
		double max2 = h2->GetMaximum();
		double max  = (max1>max2)?max1:max2;
		if(logy) max = 5*max;
		else max = 1.05*max;
		h1->SetMaximum(max);
		h2->SetMaximum(max);
		if(!logy){
			h1->SetMinimum(0.0);
			h2->SetMinimum(0.0);
		}

		TLegend *leg = new TLegend(0.55,0.75,0.75,0.88);
		if(i == 1){
			leg->AddEntry(h1, fSamples[sample1].sname,"f");
			leg->AddEntry(h2, fSamples[sample2].sname,"f");
			leg->SetFillStyle(0);
			leg->SetTextFont(42);
			leg->SetBorderSize(0);
		}

		TH1D *h1_temp = (TH1D*)h1->DrawCopy("hist");
		TH1D *h2_temp = (TH1D*)h2->DrawCopy("histsames");
		gPad->Update();
		TPaveStats *s1 = (TPaveStats*)h1_temp->GetListOfFunctions()->FindObject("stats");
		TPaveStats *s2 = (TPaveStats*)h2_temp->GetListOfFunctions()->FindObject("stats");
		s2->SetTextColor(kBlue); s2->SetLineColor(kBlue);
		s2->SetY1NDC(s1->GetY1NDC() - (s1->GetY2NDC() - s1->GetY1NDC()));
		s2->SetY2NDC(s1->GetY1NDC());

		if(i==1) leg->Draw();

		double min1 = h1->GetYaxis()->GetXmin();
		double min2 = h2->GetYaxis()->GetXmin();
		double min  = (min1<min2)?min1:min2;

		TLine *l1 = new TLine(0.15, min, 0.15, max);
		l1->SetLineColor(kRed);
		l1->SetLineWidth(2);
		l1->Draw();

		lat->SetTextColor(kBlack);
		lat->SetTextSize(0.05);
		lat->DrawLatex(0.11,0.92, Form("p_{T}(#mu) %3.0f - %3.0f GeV", getPtBins(Muon)[i-1], getPtBins(Muon)[i]));

		int bin0 = h1->FindBin(0.00);
		int bin15 = h1->FindBin(0.15);
		int bin1 = h1->FindBin(1.00);
		float f1 = h1->Integral(bin0, bin15) / h1->Integral(bin0, bin1);
		float f2 = h2->Integral(bin0, bin15) / h2->Integral(bin0, bin1);
		lat->SetTextSize(0.04);
		lat->DrawLatex(0.55,0.905, Form("ratio = %4.2f", f1));
		lat->SetTextColor(kBlue);
		lat->DrawLatex(0.55,0.945, Form("ratio = %4.2f", f2));

		gPad->RedrawAxis();
	}
	Util::PrintNoEPS(c_temp, outputname, fOutputDir, fOutputFile);
}

//____________________________________________________________________________
void MuonPlotter::makeMuIsoVsPtPlot(int sample1, int muon1, bool(MuonPlotter::*eventSelector1)(), bool(MuonPlotter::*muonSelector1)(int), int sample2, int muon2, bool(MuonPlotter::*eventSelector2)(), bool(MuonPlotter::*muonSelector2)(int), TString outputname, bool logy){
	vector<int> samples1; samples1.push_back(sample1);
	vector<int> samples2; samples2.push_back(sample2);
	makeMuIsoVsPtPlot(samples1, muon1, eventSelector1, muonSelector1, samples2, muon2, eventSelector2, muonSelector2, outputname, logy);
}
void MuonPlotter::makeMuIsoVsPtPlot(vector<int> samples1, int muon1, bool(MuonPlotter::*eventSelector1)(), bool(MuonPlotter::*muonSelector1)(int), vector<int> samples2, int muon2, bool(MuonPlotter::*eventSelector2)(), bool(MuonPlotter::*muonSelector2)(int), TString outputname, bool logy){
	const int nbins = 30;
	TH2D *h2_sig = new TH2D("h2_sig", "Isolation vs Pt for muons in data",       nbins, 0., 1., getNPt2Bins(Muon), getPt2Bins(Muon));
	TH2D *h2_bg  = new TH2D("h2_bg",  "Isolation vs Pt for fake muons in ttbar", nbins, 0., 1., getNPt2Bins(Muon), getPt2Bins(Muon));
	h2_sig->SetXTitle(convertVarName("MuIso[0]"));
	h2_bg ->SetXTitle(convertVarName("MuIso[0]"));
	h2_sig->SetYTitle(convertVarName("MuPt[0]"));
	h2_bg ->SetYTitle(convertVarName("MuPt[0]"));
	
	TTree *tree = NULL;
	for(size_t i = 0; i < samples1.size(); ++i){
		int index = samples1[i];
	
		tree = fSamples[index].tree;
		float scale = fLumiNorm / fSamples[index].lumi;
		if(fSamples[index].isdata) scale = 1.;
		tree->ResetBranchAddresses();
		Init(tree);
		if (fChain == 0) return;
		Long64_t nentries = fChain->GetEntriesFast();
		Long64_t nbytes = 0, nb = 0;
		for (Long64_t jentry=0; jentry<nentries;jentry++) {
			Long64_t ientry = LoadTree(jentry);
			if (ientry < 0) break;
			nb = fChain->GetEntry(jentry);   nbytes += nb;

			if((*this.*eventSelector1)() == false) continue;
			if((*this.*muonSelector1)(muon1) == false) continue;
			h2_sig->Fill(MuIso[muon1], MuPt[muon1], scale);
		}
	}

	for(size_t i = 0; i < samples1.size(); ++i){
		int index = samples2[i];

		tree = fSamples[index].tree;
		float scale = fLumiNorm / fSamples[index].lumi;
		if(fSamples[index].isdata) scale = 1.;
		tree->ResetBranchAddresses();
		Init(tree);
		if (fChain == 0) return;
		Long64_t nentries = fChain->GetEntriesFast();
		Long64_t nbytes = 0, nb = 0;
		for (Long64_t jentry=0; jentry<nentries;jentry++) {
			Long64_t ientry = LoadTree(jentry);
			if (ientry < 0) break;
			nb = fChain->GetEntry(jentry);   nbytes += nb;

			if((*this.*eventSelector2)() == false) continue;
			if((*this.*muonSelector2)(muon2) == false) continue;
			h2_bg->Fill(MuIso[muon2], MuPt[muon2], scale);
		}
	}

	// fSamples[sample1].tree->Project("h2_bg", Form("MuPt[%d]:MuIso[%d]", muon1, muon1), c1);
	// fSamples[sample2].tree->Project("h2_sig", Form("MuPt[%d]:MuIso[%d]", muon2, muon2), c2);

	TLatex *lat = new TLatex();
	lat->SetNDC(kTRUE);
	lat->SetTextColor(kBlack);
	lat->SetTextSize(0.06);

	TCanvas *c_temp = new TCanvas("IsoVsPt", "Isolating in Pt bins", 0, 0, 1200, 800);
	c_temp->Divide(3,2);
	// TCanvas *c_temp = new TCanvas("IsoVsPt", "Isolating in Pt bins", 0, 0, 1200, 1200);
	// c_temp->Divide(3,3);
	// c_temp->cd(9);
	// if(logy) gPad->SetLogz(1);
	// h2_sig->DrawCopy("colz");
	// lat->DrawLatex(0.11,0.92, fSamples[TTbar].sname);

	for(unsigned i = 1; i <= getNPt2Bins(Muon); ++i){
		c_temp->cd(i);
		gStyle->SetOptStat(1111);
		TH1D *h1 = h2_sig->ProjectionX(Form("h2_sigx_%d", i), i, i);	
		TH1D *h2 = h2_bg ->ProjectionX(Form("h2_bgx_%d",  i), i, i);
		h1->SetXTitle(convertVarName("MuIso[0]"));
		h1->SetLineWidth(2);
		h1->SetFillColor(15);
		h1->SetFillStyle(1001);
		h2->SetLineWidth(2);
		h2->SetLineColor(kBlue);
		h2->SetFillColor(kBlue);
		h2->SetFillStyle(3004);

		if(logy) gPad->SetLogy(1);
		gPad->SetFillStyle(0);
		h1->Sumw2();
		h2->Sumw2();

		// Scaling
		if(h1->GetEntries() > 0 ) h1->Scale(1.0/h1->Integral());
		if(h2->GetEntries() > 0 ) h2->Scale(1.0/h2->Integral());

		// setPlottingRange(h1, h2);

		// Determine plotting range
		double max1 = h1->GetMaximum();
		double max2 = h2->GetMaximum();
		double max  = (max1>max2)?max1:max2;
		if(logy) max = 5*max;
		else max = 1.05*max;
		h1->SetMaximum(max);
		h2->SetMaximum(max);
		if(!logy){
			h1->SetMinimum(0.0);
			h2->SetMinimum(0.0);
		}

		TLegend *leg = new TLegend(0.45,0.75,0.65,0.88);
		if(i == 1){
			leg->AddEntry(h1, fSamples[TTbar].sname,"f");
			leg->AddEntry(h2, fSamples[MuB].sname,"f");
			leg->SetFillStyle(0);
			leg->SetTextFont(42);
			leg->SetBorderSize(0);
		}

		TH1D *h1_temp = (TH1D*)h1->DrawCopy("hist");
		TH1D *h2_temp = (TH1D*)h2->DrawCopy("histsames");
		h1_temp->SetName(Form("h1_%d",i));
		h2_temp->SetName(Form("h2_%d",i));
		gPad->Update();
		TPaveStats *s1 = (TPaveStats*)h1_temp->GetListOfFunctions()->FindObject("stats");
		TPaveStats *s2 = (TPaveStats*)h2_temp->GetListOfFunctions()->FindObject("stats");
		s2->SetTextColor(kBlue); s2->SetLineColor(kBlue);
		s2->SetY1NDC(s1->GetY1NDC() - (s1->GetY2NDC() - s1->GetY1NDC()));
		s2->SetY2NDC(s1->GetY1NDC());

		if(i==1) leg->Draw();

		// double max1 = h1->GetYaxis()->GetXmax();
		// double max2 = h2->GetYaxis()->GetXmax();
		// double max  = (max1>max2)?max1:max2;
		double min1 = h1->GetYaxis()->GetXmin();
		double min2 = h2->GetYaxis()->GetXmin();
		double min  = (min1<min2)?min1:min2;

		TLine *l1 = new TLine(0.15, min, 0.15, max);
		l1->SetLineColor(kRed);
		l1->SetLineWidth(2);
		l1->Draw();

		lat->SetTextColor(kBlack);
		lat->SetTextSize(0.05);
		lat->DrawLatex(0.11,0.92, Form("p_{T}(#mu) %3.0f - %3.0f GeV", getPt2Bins(Muon)[i-1], getPt2Bins(Muon)[i]));

		int bin0 = h1->FindBin(0.00);
		int bin15 = h1->FindBin(0.15);
		int bin1 = h1->FindBin(1.00);
		float f1 = h1->Integral(bin0, bin15) / h1->Integral(bin0, bin1);
		float f2 = h2->Integral(bin0, bin15) / h2->Integral(bin0, bin1);
		lat->SetTextSize(0.04);
		lat->DrawLatex(0.55,0.905, Form("ratio = %4.2f", f1));
		lat->SetTextColor(kBlue);
		lat->DrawLatex(0.55,0.945, Form("ratio = %4.2f", f2));

		gPad->RedrawAxis();
	}
	Util::PrintNoEPS(c_temp, outputname, fOutputDir, fOutputFile);
}

//____________________________________________________________________________
void MuonPlotter::makeMuIsoVsNJetsPlot(int sample1, int muon1, TCut c1, int sample2, int muon2, TCut c2, TString outputname, bool logy){
	const int nbins = 20;
	const int Nnjetbins = 4;
	const double njetbins[Nnjetbins+1] = {0.,1.,2.,3.,4.};
	TH2D *h2_sig = new TH2D("h2_sig", "Isolation vs NJets for fake muons in signal", nbins, 0., 1., Nnjetbins, njetbins);
	TH2D *h2_bg  = new TH2D("h2_bg", "Isolation vs NJets for muons in background",   nbins, 0., 1., Nnjetbins, njetbins);
	h2_sig->SetXTitle(convertVarName("MuIso[0]"));
	h2_sig->SetYTitle("NJets");
	h2_bg->SetXTitle(convertVarName("MuIso[0]"));
	h2_bg->SetYTitle("NJets");

	fSamples[sample1].tree->Project("h2_bg",  Form("NJets:MuIso[%d]", muon1), c1);
	fSamples[sample2].tree->Project("h2_sig", Form("NJets:MuIso[%d]", muon2), c2);


	TLatex *lat = new TLatex();
	lat->SetNDC(kTRUE);
	lat->SetTextSize(0.06);
	lat->SetTextColor(kBlack);

	TCanvas *c_temp = new TCanvas("IsoVsNJets", "Isolating in bins of jet multiplicity", 0, 0, 1200, 800);
	c_temp->Divide(3,2);
	c_temp->cd(5);
	if(logy) gPad->SetLogz(1);
	h2_sig->DrawCopy("colz");
	lat->DrawLatex(0.11,0.92, fSamples[sample2].sname);
	c_temp->cd(6);
	if(logy) gPad->SetLogz(1);
	h2_bg->DrawCopy("colz");
	lat->DrawLatex(0.11,0.92, fSamples[sample1].sname);

	for(size_t i = 1; i <= Nnjetbins; ++i){
		c_temp->cd(i);
		gStyle->SetOptStat(1111);
		int lastbin = i;
		if(i == Nnjetbins) lastbin = i+1;
		TH1D *h1 = h2_bg->ProjectionX("h2_bgx",i, lastbin);
		TH1D *h2 = h2_sig->ProjectionX("h2_sigx",i, lastbin);
		h1->SetXTitle(convertVarName("MuIso[0]"));
		h1->SetLineWidth(2);
		h1->SetFillColor(15);
		h1->SetFillStyle(1001);
		h2->SetLineWidth(2);
		h2->SetLineColor(kBlue);
		h2->SetFillColor(kBlue);
		h2->SetFillStyle(3004);

		if(logy) gPad->SetLogy(1);
		gPad->SetFillStyle(0);
		h1->Sumw2();
		h2->Sumw2();
		if(h1->GetEntries() > 0 ) h1->Scale(1.0/h1->Integral());
		if(h2->GetEntries() > 0 ) h2->Scale(1.0/h2->Integral());

		// Determine plotting range
		double max1 = h1->GetMaximum();
		double max2 = h2->GetMaximum();
		double max  = (max1>max2)?max1:max2;
		if(logy) max = 5*max;
		else max = 1.05*max;
		h1->SetMaximum(max);
		h2->SetMaximum(max);
		if(!logy){
			h1->SetMinimum(0.0);
			h2->SetMinimum(0.0);
		}

		TLegend *leg = new TLegend(0.55,0.75,0.75,0.88);
		if(i == 1){
			leg->AddEntry(h1, fSamples[sample1].sname,"f");
			leg->AddEntry(h2, fSamples[sample2].sname,"f");
			leg->SetFillStyle(0);
			leg->SetTextFont(42);
			leg->SetBorderSize(0);
		}

		TH1D *h1_temp = (TH1D*)h1->DrawCopy("hist");
		TH1D *h2_temp = (TH1D*)h2->DrawCopy("histsames");
		gPad->Update();
		TPaveStats *s1 = (TPaveStats*)h1_temp->GetListOfFunctions()->FindObject("stats");
		TPaveStats *s2 = (TPaveStats*)h2_temp->GetListOfFunctions()->FindObject("stats");
		s2->SetTextColor(kBlue); s2->SetLineColor(kBlue);
		s2->SetY1NDC(s1->GetY1NDC() - (s1->GetY2NDC() - s1->GetY1NDC()));
		s2->SetY2NDC(s1->GetY1NDC());

		if(i==1) leg->Draw();

		double min1 = h1->GetYaxis()->GetXmin();
		double min2 = h2->GetYaxis()->GetXmin();
		double min  = (min1<min2)?min1:min2;

		TLine *l1 = new TLine(0.15, min, 0.15, max);
		l1->SetLineColor(kRed);
		l1->SetLineWidth(2);
		l1->Draw();

		lat->SetTextColor(kBlack);
		lat->SetTextSize(0.05);
		if(i < 4) lat->DrawLatex(0.11,0.92, Form("NJets = %1.0f",  njetbins[i-1]));
		else      lat->DrawLatex(0.11,0.92, Form("NJets >= %1.0f", njetbins[i-1]));

		int bin0 = h1->FindBin(0.00);
		int bin15 = h1->FindBin(0.15);
		int bin1 = h1->FindBin(1.00);
		float f1 = h1->Integral(bin0, bin15) / h1->Integral(bin0, bin1);
		float f2 = h2->Integral(bin0, bin15) / h2->Integral(bin0, bin1);
		lat->SetTextSize(0.04);
		lat->DrawLatex(0.55,0.905, Form("ratio = %4.2f", f1));
		lat->SetTextColor(kBlue);
		lat->DrawLatex(0.55,0.945, Form("ratio = %4.2f", f2));

		gPad->RedrawAxis();
	}
	Util::PrintNoEPS(c_temp, outputname, fOutputDir, fOutputFile);
}

//____________________________________________________________________________
void MuonPlotter::produceRatio(gChannel chan, int sample, int index, bool(MuonPlotter::*eventSelector)(), bool(MuonPlotter::*objSelector)(int), TH2D *&h_2d, TH1D *&h_pt, TH1D *&h_eta, bool output){
	vector<int> samples; samples.push_back(sample);
	produceRatio(chan, samples, index, eventSelector, objSelector, h_2d, h_pt, h_eta, output);
}
void MuonPlotter::produceRatio(gChannel chan, vector<int> samples, int index, bool(MuonPlotter::*eventSelector)(), bool(MuonPlotter::*objSelector)(int), TH2D *&h_2d, TH1D *&h_pt, TH1D *&h_eta, bool output){
// Base function for production of all ratios
/*
	TODO Fix treatment of statistical errors and luminosity scaling here!
*/
	gStyle->SetOptStat(0);
	h_2d->Sumw2();
	h_pt->Sumw2();
	h_eta->Sumw2();
	TString sname = "Mu";
	TString fname = "Muon";
	if(chan == Electron){
		sname = "El";
		fname = "Electron";
	}
	TH2D *H_ntight = new TH2D(Form("%sNTight", sname.Data()), Form("NTight %ss", fname.Data()), h_2d->GetNbinsX(), h_2d->GetXaxis()->GetXbins()->GetArray(), h_2d->GetNbinsY(),  h_2d->GetYaxis()->GetXbins()->GetArray());
	TH2D *H_nloose = new TH2D(Form("%sNLoose", sname.Data()), Form("NLoose %ss", fname.Data()), h_2d->GetNbinsX(), h_2d->GetXaxis()->GetXbins()->GetArray(), h_2d->GetNbinsY(),  h_2d->GetYaxis()->GetXbins()->GetArray());
	H_ntight->Sumw2();
	H_nloose->Sumw2();

	if(fVerbose>2) cout << "---------------" << endl;
	for(size_t i = 0; i < samples.size(); ++i){
		int sample = samples[i];
		TTree *tree = fSamples[sample].tree;
		if(fVerbose>2) cout << "Producing ratios for " << fSamples[sample].sname << endl;
		tree->ResetBranchAddresses();
		Init(tree);
		if (fChain == 0) return;
		Long64_t nentries = fChain->GetEntriesFast();

		float scale = fLumiNorm / fSamples[sample].lumi;
		if(fSamples[sample].isdata) scale = 1;

		Long64_t nbytes = 0, nb = 0;
		for (Long64_t jentry=0; jentry<nentries;jentry++) {
			Long64_t ientry = LoadTree(jentry);
			if (ientry < 0) break;
			nb = fChain->GetEntry(jentry);   nbytes += nb;
			printProgress(jentry, nentries, fSamples[sample].name);

			if((*this.*eventSelector)() == false) continue;
			if((*this.*objSelector)(index) == false) continue;

			if(chan == Muon){
				if(isLooseMuon(index)) H_nloose->Fill(MuPt[index], MuEta[index], scale); // Tight or loose
				if(isTightMuon(index)) H_ntight->Fill(MuPt[index], MuEta[index], scale); // Tight
			}
			if(chan == Electron){
				if(isLooseElectron(index)) H_nloose->Fill(ElPt[index], ElEta[index], scale); // Tight or loose
				if(isTightElectron(index)) H_ntight->Fill(ElPt[index], ElEta[index], scale); // Tight
			}

		}
		cout << endl;

		if(fVerbose>2) cout << " Tight entries so far: " << H_ntight->GetEntries() << " / " << H_ntight->Integral() << endl;
		if(fVerbose>2) cout << " Loose entries so far: " << H_nloose->GetEntries() << " / " << H_nloose->Integral() << endl;
		if(fVerbose>2) cout << "  Ratio: " << (double)H_ntight->GetEntries()/(double)H_nloose->GetEntries() << endl;
	}
	h_2d->Divide(H_ntight, H_nloose, 1., 1., "B"); // binomial, weights are ignored

	TH1D *hloosept  = H_nloose->ProjectionX();
	TH1D *hlooseeta = H_nloose->ProjectionY();
	TH1D *htightpt  = H_ntight->ProjectionX();
	TH1D *htighteta = H_ntight->ProjectionY();
	h_pt ->Divide(htightpt,  hloosept,  1., 1., "B");
	h_eta->Divide(htighteta, hlooseeta, 1., 1., "B");

	h_pt ->SetXTitle(convertVarName(sname + "Pt[0]"));
	h_eta->SetXTitle(convertVarName(sname + "Eta[0]"));
	h_pt ->SetYTitle("# Tight / # Loose");
	h_eta->SetYTitle("# Tight / # Loose");
	h_2d->SetXTitle(convertVarName(sname + "Pt[0]"));
	h_2d->SetYTitle(convertVarName(sname + "Eta[0]"));
	TString name = "";
	for(size_t i = 0; i < samples.size(); ++i){
		int sample = samples[i];
		name += h_2d->GetName();
		name += "_";
		name += fSamples[sample].sname;
	}
	if(output){
		printObject(h_2d,  sname + "Ratio"    + name, "Muon Fake Ratio vs pt/eta", "colz");
		printObject(h_pt,  sname + "RatioPt"  + name, "Muon Fake Ratio vs pt",     "PE1");
		printObject(h_eta, sname + "RatioEta" + name, "Muon Fake Ratio vs eta",    "PE1");
	}
	delete H_ntight, H_nloose, hloosept, hlooseeta, htightpt, htighteta;
}

//____________________________________________________________________________
void MuonPlotter::plotMuRatio(int sample, int muon, bool(MuonPlotter::*eventSelector)(), bool(MuonPlotter::*muonSelector)(int), TString tag){
	vector<int> samples; samples.push_back(sample);
	plotMuRatio(samples, muon, eventSelector, muonSelector, tag);
}
void MuonPlotter::plotMuRatio(vector<int> samples, int muon, bool(MuonPlotter::*eventSelector)(), bool(MuonPlotter::*muonSelector)(int), TString tag){
	gStyle->SetOptStat(0);
	TH2D *h_2d  = new TH2D("MuRatio",    Form("Ratio of tight to loose Muons vs Pt vs Eta : %s", tag.Data()), getNPtBins(Muon), getPtBins(Muon), getNEtaBins(Muon), getEtaBins(Muon));
	TH1D *h_pt  = new TH1D("MuRatioPt",  Form("Ratio of tight to loose Muons vs Pt : %s",        tag.Data()), getNPtBins(Muon), getPtBins(Muon));
	TH1D *h_eta = new TH1D("MuRatioEta", Form("Ratio of tight to loose Muons vs Eta : %s",       tag.Data()), getNEtaBins(Muon), getEtaBins(Muon));
	h_pt->SetXTitle(convertVarName("MuPt[0]"));
	h_eta->SetXTitle(convertVarName("MuEta[0]"));
	h_2d->SetXTitle(convertVarName("MuPt[0]"));
	h_2d->SetYTitle(convertVarName("MuEta[0]"));
	h_pt ->SetYTitle("# Tight / # Loose");
	h_eta->SetYTitle("# Tight / # Loose");
	h_pt->GetYaxis()->SetTitleOffset(1.2);
	h_eta->GetYaxis()->SetTitleOffset(1.2);

	produceRatio(Muon, samples, muon, eventSelector, muonSelector, h_2d, h_pt, h_eta, true);
}

//____________________________________________________________________________
TH1D* MuonPlotter::fillMuRatioPt(int sample, int muon, bool(MuonPlotter::*eventSelector)(), bool(MuonPlotter::*muonSelector)(int), bool output){
	vector<int> samples; samples.push_back(sample);
	return fillMuRatioPt(samples, muon, eventSelector, muonSelector, output);
}
TH1D* MuonPlotter::fillMuRatioPt(vector<int> samples, int muon, bool(MuonPlotter::*eventSelector)(), bool(MuonPlotter::*muonSelector)(int), bool output){
	gStyle->SetOptStat(0);
	TH2D *h_2d  = new TH2D("MuRatio",    "Ratio of tight to loose Muons vs Pt vs Eta", getNPt2Bins(Muon), getPt2Bins(Muon), getNEtaBins(Muon), getEtaBins(Muon));
	TH1D *h_pt  = new TH1D("MuRatioPt",  "Ratio of tight to loose Muons vs Pt",        getNPt2Bins(Muon), getPt2Bins(Muon));
	TH1D *h_eta = new TH1D("MuRatioEta", "Ratio of tight to loose Muons vs Eta",       getNEtaBins(Muon), getEtaBins(Muon));

	h_pt->SetXTitle(convertVarName("MuPt[0]"));
	h_pt ->SetYTitle("# Tight / # Loose");
	h_pt->GetYaxis()->SetTitleOffset(1.2);

	produceRatio(Muon, samples, muon, eventSelector, muonSelector, h_2d, h_pt, h_eta, output);
	return h_pt;
}
TH1D* MuonPlotter::fillMuRatioPt(vector<int> samples, int muon, bool(MuonPlotter::*eventSelector)(), bool(MuonPlotter::*muonSelector)(int), const int nptbins, const double* ptbins, const int netabins, const double* etabins, bool output){
	gStyle->SetOptStat(0);
	TH2D *h_2d  = new TH2D("MuRatio",    "Ratio of tight to loose Muons vs Pt vs Eta", nptbins, ptbins, netabins, etabins);
	TH1D *h_pt  = new TH1D("MuRatioPt",  "Ratio of tight to loose Muons vs Pt",        nptbins, ptbins);
	TH1D *h_eta = new TH1D("MuRatioEta", "Ratio of tight to loose Muons vs Eta",       netabins, etabins);

	h_pt->SetXTitle(convertVarName("MuPt[0]"));
	h_pt ->SetYTitle("# Tight / # Loose");
	h_pt->GetYaxis()->SetTitleOffset(1.2);

	produceRatio(Muon, samples, muon, eventSelector, muonSelector, h_2d, h_pt, h_eta, output);
	return h_pt;
}

//____________________________________________________________________________
TH2D* MuonPlotter::fillMuRatio(int sample, int muon, bool(MuonPlotter::*eventSelector)(), bool(MuonPlotter::*muonSelector)(int), const int nptbins, const double* ptbins, const int netabins, const double* etabins){
	vector<int> samples; samples.push_back(sample);
	return fillMuRatio(samples, muon, eventSelector, muonSelector, nptbins, ptbins, netabins, etabins);
}
TH2D* MuonPlotter::fillMuRatio(vector<int> samples, int muon, bool(MuonPlotter::*eventSelector)(), bool(MuonPlotter::*muonSelector)(int), const int nptbins, const double* ptbins, const int netabins, const double* etabins){
	gStyle->SetOptStat(0);
	TH2D *h_2d  = new TH2D("MuRatio",    "Ratio of tight to loose Muons vs Pt vs Eta", nptbins, ptbins, netabins, etabins);
	TH1D *h_pt  = new TH1D("MuRatioPt",  "Ratio of tight to loose Muons vs Pt",        nptbins, ptbins);
	TH1D *h_eta = new TH1D("MuRatioEta", "Ratio of tight to loose Muons vs Eta",       netabins, etabins);

	h_2d->SetXTitle(convertVarName("MuPt[0]"));
	h_2d->SetYTitle(convertVarName("MuEta[0]"));

	produceRatio(Muon, samples, muon, eventSelector, muonSelector, h_2d, h_pt, h_eta, false);
	return h_2d;
}

//____________________________________________________________________________
void MuonPlotter::fillMuElRatios(vector<int> samples){
	calculateRatio(samples, Muon,     SigSup, fH2D_MufRatio, fH1D_MufRatioPt, fH1D_MufRatioEta);
	calculateRatio(samples, Muon,     ZDecay, fH2D_MupRatio, fH1D_MupRatioPt, fH1D_MupRatioEta);
	calculateRatio(samples, Electron, SigSup, fH2D_ElfRatio, fH1D_ElfRatioPt, fH1D_ElfRatioEta);
	calculateRatio(samples, Electron, ZDecay, fH2D_ElpRatio, fH1D_ElpRatioPt, fH1D_ElpRatioEta);
}

//____________________________________________________________________________
TH1D* MuonPlotter::fillMuRatioPt(int sample, gRegion reg, bool output){
	vector<int> samples; samples.push_back(sample);
	return fillMuRatioPt(samples, reg);
}
TH1D* MuonPlotter::fillMuRatioPt(vector<int> samples, gRegion reg, bool output){
	gStyle->SetOptStat(0);
	TH2D *h_2d  = new TH2D("MuRatio",    "Ratio of tight to loose Muons vs Pt vs Eta", getNPt2Bins(Muon), getPt2Bins(Muon), getNEtaBins(Muon), getEtaBins(Muon));
	TH1D *h_pt  = new TH1D("MuRatioPt",  "Ratio of tight to loose Muons vs Pt",        getNPt2Bins(Muon), getPt2Bins(Muon));
	TH1D *h_eta = new TH1D("MuRatioEta", "Ratio of tight to loose Muons vs Eta",       getNEtaBins(Muon), getEtaBins(Muon));

	h_pt->SetXTitle(convertVarName("MuPt[0]"));
	h_pt->SetYTitle("# Tight / # Loose");
	h_pt->GetYaxis()->SetTitleOffset(1.2);

	calculateRatio(samples, Muon, reg, h_2d, h_pt, h_eta, output);
	delete h_2d, h_eta;
	return h_pt;
};

//____________________________________________________________________________
TH1D* MuonPlotter::fillElRatioPt(int sample, gRegion reg, bool output){
	vector<int> samples; samples.push_back(sample);
	return fillElRatioPt(samples, reg, output);
}
TH1D* MuonPlotter::fillElRatioPt(vector<int> samples, gRegion reg, bool output){
	gStyle->SetOptStat(0);
	TH2D *h_2d  = new TH2D("ElRatio",    "Ratio of tight to loose Electrons vs Pt vs Eta", getNPt2Bins(Electron), getPt2Bins(Electron), getNEtaBins(Electron), getEtaBins(Electron));
	TH1D *h_pt  = new TH1D("ElRatioPt",  "Ratio of tight to loose Electrons vs Pt",        getNPt2Bins(Electron), getPt2Bins(Electron));
	TH1D *h_eta = new TH1D("ElRatioEta", "Ratio of tight to loose Electrons vs Eta",       getNEtaBins(Electron), getEtaBins(Electron));

	h_pt->SetXTitle(convertVarName("ElPt[0]"));
	h_pt->SetYTitle("# Tight / # Loose");
	h_pt->GetYaxis()->SetTitleOffset(1.2);

	calculateRatio(samples, Electron, reg, h_2d, h_pt, h_eta, output);
	return h_pt;
};

//____________________________________________________________________________
void MuonPlotter::calculateRatio(vector<int> samples, gChannel chan, gRegion reg, TH2D*& h_2d, bool output){
	TH1D *h_dummy1 = new TH1D("dummy1", "dummy1", 1, 0.,1.);
	TH1D *h_dummy2 = new TH1D("dummy2", "dummy2", 1, 0.,1.);
	calculateRatio(samples, chan, reg, h_2d, h_dummy1, h_dummy2, output);
	delete h_dummy1, h_dummy2;
}
void MuonPlotter::calculateRatio(vector<int> samples, gChannel chan, gRegion reg, TH2D*& h_2d, TH1D*& h_pt, TH1D*&h_eta, bool output){
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

	getPassedTotal(samples, chan, reg, H_ntight, H_nloose, output);
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
		name += fSamples[sample].sname;
	}
	if(output){
		printObject(h_2d,  TString("Ratio")    + name, "Fake Ratio vs pt/eta", "colz");
		printObject(h_pt,  TString("RatioPt")  + name, "Fake Ratio vs pt",     "PE1");
		printObject(h_eta, TString("RatioEta") + name, "Fake Ratio vs eta",    "PE1");
	}
}
void MuonPlotter::calculateRatio(vector<int> samples, gChannel chan, gRegion reg, float &ratio, float &ratioe){
	double ntight(0.), nloose(0.);
	double ntighte2(0.), nloosee2(0.);
	vector<int> v_ntight, v_nloose;
	vector<float> v_scale;
	vector<TString> v_name;
	for(size_t i = 0; i < samples.size(); ++i){
		int index = samples[i];
		Sample S = fSamples[index];
		
		int ntight_sam(0), nloose_sam(0);
		v_name.push_back(S.sname);
		
		float scale = fLumiNorm/S.lumi; // Normalize all
		if(S.isdata) scale = 1;
		if(reg == SigSup){
			ntight += scale * S.numbers[chan].nsst;
			nloose += scale * S.numbers[chan].nssl;
			
			ntight_sam += S.numbers[chan].nsst;
			nloose_sam += S.numbers[chan].nssl;
		}
		if(reg == ZDecay){
			ntight += scale * S.numbers[chan].nzt;
			nloose += scale * S.numbers[chan].nzl;

			ntight_sam += S.numbers[chan].nzt;
			nloose_sam += S.numbers[chan].nzl;
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
		if(reg == SigSup) s_forp = "f Ratio";
		if(reg == ZDecay) s_forp = "p Ratio";
		TString s_channel;
		if(chan == Muon)     s_channel = "Muon";
		if(chan == Electron) s_channel = "Electron";
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
void MuonPlotter::calculateRatio(vector<int> samples, gChannel chan, gRegion reg, float &ratio, float &ratioeup, float &ratioelow){
	// Careful, this method only takes integer numbers for passed/total events, therefore
	// only makes sense for application on data right now.
	int ntight(0), nloose(0);
	float ntighte2(0.), nloosee2(0.);
	vector<int> v_ntight, v_nloose;
	vector<TString> v_name;
	for(size_t i = 0; i < samples.size(); ++i){
		Sample S = fSamples[samples[i]];
		
		int ntight_sam(0), nloose_sam(0);
		v_name.push_back(S.sname);
		
		float scale = fLumiNorm/S.lumi; // Normalize all
		if(S.isdata) scale = 1;
		// Channel *cha;
		// if(chan == Muon)     cha = &S.region[reg].mm;
		// if(chan == Electron) cha = &S.region[reg].ee;
		if(reg == SigSup){
			ntight += scale * S.numbers[chan].nsst;
			nloose += scale * S.numbers[chan].nssl;
			
			ntight_sam += S.numbers[chan].nsst;
			nloose_sam += S.numbers[chan].nssl;
		}
		if(reg == ZDecay){
			ntight += S.numbers[chan].nzt;
			nloose += S.numbers[chan].nzl;

			ntight_sam += S.numbers[chan].nzt;
			nloose_sam += S.numbers[chan].nzl;
		}
		v_ntight.push_back(ntight_sam);
		v_nloose.push_back(nloose_sam);
	}
	ratioWithAsymmCPErrors(ntight, nloose, ratio, ratioeup, ratioelow);
	if(fVerbose > 2){
		cout << "--------------------------------------------------------" << endl;
		TString s_forp;
		if(reg == SigSup) s_forp = "f Ratio";
		if(reg == ZDecay) s_forp = "p Ratio";
		TString s_channel;
		if(chan == Muon)     s_channel = "Muon";
		if(chan == Electron) s_channel = "Electron";
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
TGraphAsymmErrors* MuonPlotter::combineMCEfficiencies(vector<int> samples, gChannel chan, gRegion reg, bool output, TEfficiency::EStatOption statopt, double betaalpha, double betabeta){
	bool data = fSamples[samples[0]].isdata;
	for(size_t i = 0; i < samples.size(); ++i){ // check if samples contains only data or only MC
		if(fSamples[samples[i]].isdata == data) continue;
		else{
			TGraphAsymmErrors *null = NULL;
			cout << "MuonPlotter::calculateEfficiency ==> sample is not pure data or pure MC, aborting..." << endl;
			return null;
		}
	}

	TString regchan = "";
	if(chan == Muon)     regchan += "Mu_";
	if(chan == Electron) regchan += "El_";
	if(reg == SigSup) regchan += "f_";
	if(reg == ZDecay) regchan += "p_";
	TString name1 = regchan + "TLEff";
	TString name2 = regchan + "TLGraph";
	for(size_t i = 0; i < samples.size(); ++i){
		name1 += "_";
		name2 += "_";
		name1 += fSamples[samples[i]].sname;
		name2 += fSamples[samples[i]].sname;
	}
	
	TEfficiency *total_eff = new TEfficiency(name1, "TotalEfficiency", getNPt2Bins(chan), getPt2Bins(chan));
	total_eff->SetName(name1);
	total_eff->SetStatisticOption(statopt); // Needs to be bayesian in order to combine different effs
	total_eff->SetBetaAlpha(betaalpha);
	total_eff->SetBetaBeta(betabeta);
	total_eff->SetConfidenceLevel(0.683); // 1-sigma = default
	TList *pList = new TList();
	for(size_t i = 0; i < samples.size(); ++i){ // get the individual efficiencies
		Sample *S = &fSamples[samples[i]];
		TEfficiency *eff = getEfficiency(S, chan, reg, 1, output);
		eff->SetStatisticOption(total_eff->GetStatisticOption());
		eff->SetBetaAlpha(total_eff->GetBetaAlpha());
		eff->SetBetaBeta(total_eff->GetBetaBeta());

		// Weighting
		eff->SetWeight(fLumiNorm / S->lumi);
		// if(eff->GetTotalHistogram()->GetEntries() == 0) eff->SetWeight(0.);
		// else eff->SetWeight(fLumiNorm / S->lumi);
		pList->Add(eff);
		// if(eff->GetTotalHistogram()->GetEntries() > 0) pList->Add(eff);
		if(output){
			if(!fOutputSubDir.EndsWith("/")) fOutputSubDir += "/";
			TString temp = fOutputSubDir;
			fOutputSubDir = temp + "Effs/";
			printObject(eff, regchan + "Eff_" + S->sname, "Test", "");
			fOutputSubDir = temp;
		}
	}
	
	TGraphAsymmErrors *result = new TGraphAsymmErrors();
	result = total_eff->Combine(pList, "mode central");
	result->SetName(name2);
	return result;
}
TEfficiency* MuonPlotter::mergeDataEfficiencies(vector<int> samples, gChannel chan, gRegion reg, bool output, TEfficiency::EStatOption statopt, double betaalpha, double betabeta){
	bool data = fSamples[samples[0]].isdata;
	for(size_t i = 0; i < samples.size(); ++i){ // check if samples contains only data or only MC
		if(fSamples[samples[i]].isdata == data) continue;
		else{
			TEfficiency *null = NULL;
			cout << "MuonPlotter::calculateEfficiency ==> sample is not pure data or pure MC, aborting..." << endl;
			return null;
		}
	}

	TString regchan = "";
	if(chan == Muon)     regchan += "Mu_";
	if(chan == Electron) regchan += "El_";
	if(reg == SigSup) regchan += "f_";
	if(reg == ZDecay) regchan += "p_";
	TString name1 = regchan + "TLEff";
	TString name2 = regchan + "TLGraph";
	for(size_t i = 0; i < samples.size(); ++i){
		name1 += "_";
		name2 += "_";
		name1 += fSamples[samples[i]].sname;
		name2 += fSamples[samples[i]].sname;
	}

	TEfficiency *total_eff = new TEfficiency(name1, "TotalEfficiency", getNPt2Bins(chan), getPt2Bins(chan));
	total_eff->SetName(name1);
	total_eff->SetStatisticOption(statopt);
	total_eff->SetBetaAlpha(betaalpha);
	total_eff->SetBetaBeta(betabeta);
	total_eff->SetConfidenceLevel(0.683); // 1-sigma = default
	TList *pList = new TList();
	for(size_t i = 0; i < samples.size(); ++i){ // get the individual efficiencies
		Sample *S = &fSamples[samples[i]];
		double scale = fLumiNorm / S->lumi;
		TEfficiency *eff = getEfficiency(S, chan, reg, 1, output);
		pList->Add(eff);
		if(output){
			if(!fOutputSubDir.EndsWith("/")) fOutputSubDir += "/";
			TString temp = fOutputSubDir;
			fOutputSubDir = temp + "Effs/";
			printObject(eff, regchan + "Eff_" + S->sname, "Test", "");
			fOutputSubDir = temp;
		}
	}
	total_eff->Merge(pList);
	return total_eff;
}

//____________________________________________________________________________
TEfficiency* MuonPlotter::getEfficiency(Sample *S, gChannel chan, gRegion reg, int pteta2d, bool output){
// Basic method for retrieving efficiency of a sample, used for both data and mc
// pteta2d switches between different binnings: 0 = 2d (pt vs eta), 1 = pt, 2 = eta
	Region *R = &S->region[Signal];
	Channel *C;
	if(chan == Muon)     C = &R->mm;
	if(chan == Electron) C = &R->ee;

	TH2D *ntight, *nloose;
	if(reg == SigSup){
		ntight = C->fntight;
		nloose = C->fnloose;
	} else if(reg == ZDecay){
		ntight = C->pntight;
		nloose = C->pnloose;		
	}


	if(output){
		if(!fOutputSubDir.EndsWith("/")) fOutputSubDir += "/";
		TString temp = fOutputSubDir;
		fOutputSubDir = temp + "NTight/";
		printObject(ntight->ProjectionX(), S->region[reg].sname + "_" + C->sname + "_Tight_" + S->name, "Number of loose leptons vs pt", "PE1");
		fOutputSubDir = temp + "NLoose/";
		printObject(nloose->ProjectionX(), S->region[reg].sname + "_" + C->sname + "_Loose_" + S->name, "Number of tight leptons vs pt", "PE1");
		TH1D *rat = new TH1D(S->region[reg].sname + "_" + C->sname + "_Ratio_" + S->sname, "Title", getNPt2Bins(chan), getPt2Bins(chan));
		rat->Divide(ntight->ProjectionX(), nloose->ProjectionX(), 1., 1., "B"); // binomial
		fOutputSubDir = temp + "Ratio/";
		printObject(rat, rat->GetName(), "Title", "PE1");
		fOutputSubDir = temp;
	}

	TEfficiency *eff;
	if(pteta2d == 0){
		eff = new TEfficiency(*ntight, *nloose);
		eff->SetName(S->sname + "_TLEff_2d");
		eff->SetTitle("Tight/Loose Efficiency for " + S->sname + " vs pt vs eta");
	}
	if(pteta2d == 1){
		TH1D *h_passed = ntight->ProjectionX();
		TH1D *h_total  = nloose->ProjectionX();
		eff = new TEfficiency(*h_passed, *h_total);	
		eff->SetName(S->sname + "_TLEff_pt");
		eff->SetTitle("Tight/Loose Efficiency for " + S->sname + " vs pt");
	} 
	if(pteta2d == 2){
		TH1D *h_passed = ntight->ProjectionY();
		TH1D *h_total  = nloose->ProjectionY();
		eff = new TEfficiency(*h_passed, *h_total);	
		eff->SetName(S->sname + "_TLEff_eta");
		eff->SetTitle("Tight/Loose Efficiency for " + S->sname + " vs eta");
	}
	return eff;
}

//____________________________________________________________________________
void MuonPlotter::getPassedTotal(vector<int> samples, gChannel chan, gRegion reg, TH2D*& h_passed, TH2D*& h_total, bool output){
	if(fVerbose>2) cout << "---------------" << endl;
	for(size_t i = 0; i < samples.size(); ++i){
		Sample *S = &fSamples[samples[i]];

		float scale = fLumiNorm / S->lumi;
		if(S->isdata) scale = 1;

		Channel *C;
		if(chan == Muon)     C = &S->region[Signal].mm;
		if(chan == Electron) C = &S->region[Signal].ee;
		TH2D *ntight, *nloose;
		if(reg == SigSup){
			ntight = C->fntight;
			nloose = C->fnloose;
		} else if(reg == ZDecay){
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
		name += fSamples[sample].sname;
	}
	if(output){
		printObject(h_passed, TString("Passed") + name, "Passed vs pt vs eta", "colz");
		printObject(h_total,  TString("Total")  + name, "Total vs pt vs eta",  "colz");
	}	
}

//____________________________________________________________________________
void MuonPlotter::ratioWithBinomErrors(float ntight, float nloose, float &ratio, float &error){
	ratio = ntight/nloose;
	error = TMath::Sqrt( ntight*(1.0-ntight/nloose) ) / nloose;                  // Binomial
}
void MuonPlotter::ratioWithPoissErrors(float ntight, float nloose, float &ratio, float &error){
	ratio = ntight/nloose;
	error = TMath::Sqrt( ntight*ntight*(nloose+ntight)/(nloose*nloose*nloose) ); // Poissonian	
}
void MuonPlotter::ratioWithAsymmCPErrors(int passed, int total, float &ratio, float &upper, float &lower){
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
void MuonPlotter::makeSSMuMuPredictionPlots(vector<int> samples, bool output){
	fOutputSubDir = "MuMuPredictions";
	// Need filled muon ratios before calling this function!

	// Prediction: /////////////////////////////////////////////////////////////////////
	TH1D *H_nsigpred = new TH1D("MuMuNsigpred", "Predicted N_sig in Pt1 bins",       getNPt2Bins(Muon),  getPt2Bins(Muon));
	TH1D *H_nfppred  = new TH1D("MuMuNfppred",  "Predicted N_fp in Pt1 bins",        getNPt2Bins(Muon),  getPt2Bins(Muon));
	TH1D *H_nffpred  = new TH1D("MuMuNffpred",  "Predicted N_ff in Pt1 bins",        getNPt2Bins(Muon),  getPt2Bins(Muon));
	TH1D *H_nFpred   = new TH1D("MuMuNFpred",   "Total predicted fakes in Pt1 bins", getNPt2Bins(Muon),  getPt2Bins(Muon));
	for(size_t i = 0; i < samples.size(); ++i){
		Sample *S = &fSamples[samples[i]];
		float scale = fLumiNorm/S->lumi;
		Channel *C = &S->region[Signal].mm;

		vector<TH1D*> prediction = MuMuFPPrediction(fH2D_MufRatio, fH2D_MupRatio, C->nt20_pt, C->nt10_pt, C->nt00_pt, output);
		H_nsigpred->Add(prediction[0], scale);
		H_nfppred ->Add(prediction[1], scale);
		H_nffpred ->Add(prediction[2], scale);
	}

	// Observation: ////////////////////////////////////////////////////////////////////
	TH1D *H_nsigobs  = new TH1D("Nsigobs", "Observed N_sig in Pt1 bins",  getNPt2Bins(Muon),  getPt2Bins(Muon));
	vector<int> lm0sample; lm0sample.push_back(LM0);
	NObs(Muon, H_nsigobs, lm0sample, &MuonPlotter::isGenMatchedSUSYDiLepEvent);

	TH1D *H_nt2obs  = new TH1D("Nt2obs", "Observed N_t2 in Pt1 bins",  getNPt2Bins(Muon),  getPt2Bins(Muon));
	NObs(Muon, H_nt2obs, samples);
	// NObs(Muon, H_nt2obs, samples, &MuonPlotter::isSSTTEvent);

	// THStack *HS_nt2obs = new THStack("Nt2obs_stacked", "Observed N_t2 in Pt1 bins");
	// NObs(Muon, HS_nt2obs, samples);

	TH1D *H_nt2obsSM = new TH1D("Nt2obsSM", "Observed N_t2 in Pt1 bins, ttbar only",  getNPt2Bins(Muon),  getPt2Bins(Muon));
	NObs(Muon, H_nt2obsSM, fMCBG);
	// vector<int> ttbarsample; ttbarsample.push_back(TTbar);
	// NObs(Muon, H_nt2obsttbar, ttbarsample, &MuonPlotter::isSSTTEvent);	

	// Output
	H_nsigobs->SetXTitle(convertVarName("MuPt[0]"));
	H_nsigobs->SetYTitle(Form("Events / %2.0f GeV", fBinWidthScale));

	H_nFpred->Add(H_nfppred);
	H_nFpred->Add(H_nffpred);
	H_nFpred->SetXTitle(convertVarName("MuPt[0]"));
	H_nFpred->SetYTitle(Form("Events / %2.0f GeV", fBinWidthScale));
	H_nt2obs->SetXTitle(convertVarName("MuPt[0]"));
	H_nt2obs->SetYTitle(Form("Events / %2.0f GeV", fBinWidthScale));
	H_nt2obsSM->SetXTitle(convertVarName("MuPt[0]"));
	H_nt2obsSM->SetYTitle(Form("Events / %2.0f GeV", fBinWidthScale));

	// Normalize to binwidth
	H_nsigpred = normHistBW(H_nsigpred, fBinWidthScale);
	H_nsigobs  = normHistBW(H_nsigobs,  fBinWidthScale);
	H_nfppred  = normHistBW(H_nfppred,  fBinWidthScale);
	H_nffpred  = normHistBW(H_nffpred,  fBinWidthScale);
	H_nFpred   = normHistBW(H_nFpred,   fBinWidthScale);
	H_nt2obs   = normHistBW(H_nt2obs,   fBinWidthScale);
	H_nt2obsSM = normHistBW(H_nt2obsSM, fBinWidthScale);

	H_nt2obs->SetFillColor(kBlue);
	H_nt2obs->SetLineColor(kBlue);
	H_nt2obs->SetFillStyle(3004);
	H_nt2obs->SetLineWidth(2);

	H_nsigpred->SetFillColor(8);
	H_nsigpred->SetMarkerColor(8);
	H_nsigpred->SetMarkerStyle(20);
	H_nsigpred->SetLineColor(8);
	H_nsigpred->SetLineWidth(2);

	H_nfppred->SetFillColor(kRed);
	H_nfppred->SetMarkerColor(kRed);
	H_nfppred->SetMarkerStyle(20);
	H_nfppred->SetLineColor(kRed);
	H_nfppred->SetLineWidth(2);

	H_nffpred->SetFillColor(13);
	H_nffpred->SetMarkerColor(13);
	H_nffpred->SetLineColor(13);
	H_nffpred->SetMarkerStyle(20);
	H_nffpred->SetLineWidth(2);

	plotOverlay4H(H_nt2obs, "N_{ t2}", H_nsigpred, "N_{ pp}" , H_nfppred, "N_{ f p}", H_nffpred, "N_{ f f}");

	H_nsigpred->SetMarkerColor(kRed);
	H_nFpred->SetMarkerColor(kRed);
	H_nFpred->SetMarkerStyle(20);
	
	plotPredOverlay2HWithRatio(H_nsigobs, "Observed N_{sig}", H_nsigpred, "Predicted N_{sig}");
	plotPredOverlay2HWithRatio(H_nt2obs,  "Observed N_{t2}",  H_nFpred,   "Predicted Fakes");
	plotPredOverlay3HWithRatio(H_nFpred,  "Predicted Fakes",  H_nt2obs,   "Obs. N_{t2} (QCD, t#bar{t}+jets, V+jets, LM0)", H_nt2obsSM, "Obs. N_{t2} (QCD, t#bar{t}+jets, V+jets)", false, false);
	fOutputSubDir = "";
}
void MuonPlotter::makeSSElElPredictionPlots(vector<int> samples, bool output){
	fOutputSubDir = "ElElPredictions";
	// Need filled electron ratios before calling this function!

	// Prediction: /////////////////////////////////////////////////////////////////////
	TH1D *H_nsigpred = new TH1D("ElElNsigpred", "Predicted N_sig in Pt1 bins",       getNPt2Bins(Electron),  getPt2Bins(Electron));
	TH1D *H_nfppred  = new TH1D("ElElNfppred",  "Predicted N_fp in Pt1 bins",        getNPt2Bins(Electron),  getPt2Bins(Electron));
	TH1D *H_nffpred  = new TH1D("ElElNffpred",  "Predicted N_ff in Pt1 bins",        getNPt2Bins(Electron),  getPt2Bins(Electron));
	TH1D *H_nFpred   = new TH1D("ElElNFpred",   "Total predicted fakes in Pt1 bins", getNPt2Bins(Electron),  getPt2Bins(Electron));
	for(size_t i = 0; i < samples.size(); ++i){
		Sample *S = &fSamples[samples[i]];
		float scale = fLumiNorm/S->lumi;

		Channel *C = &S->region[Signal].ee;
		vector<TH1D*> prediction = ElElFPPrediction(fH2D_ElfRatio, fH2D_ElpRatio, C->nt20_pt, C->nt10_pt, C->nt00_pt, output);
		H_nsigpred->Add(prediction[0], scale);
		H_nfppred ->Add(prediction[1], scale);
		H_nffpred ->Add(prediction[2], scale);
	}

	// Observation: ////////////////////////////////////////////////////////////////////
	TH1D *H_nsigobs  = new TH1D("Nsigobs", "Observed N_sig in Pt1 bins",  getNPt2Bins(Electron),  getPt2Bins(Electron));
	vector<int> lm0sample; lm0sample.push_back(LM0);
	NObs(Electron, H_nsigobs, lm0sample, &MuonPlotter::isGenMatchedSUSYEEEvent);

	TH1D *H_nt2obs  = new TH1D("Nt2obs", "Observed N_t2 in Pt1 bins",  getNPt2Bins(Electron),  getPt2Bins(Electron));
	NObs(Electron, H_nt2obs, samples);
	// NObs(H_nt2obs, samples, &MuonPlotter::isSSTTEvent);

	TH1D *H_nt2obsSM = new TH1D("Nt2obsSM", "Observed N_t2 in Pt1 bins, ttbar only",  getNPt2Bins(Electron),  getPt2Bins(Electron));
	// vector<int> ttbarsample; ttbarsample.push_back(TTbar);
	NObs(Electron, H_nt2obsSM, fMCBG);
	// NObs(Electron, H_nt2obsttbar, ttbarsample, &MuonPlotter::isSSTTEvent);	

	// Output
	H_nsigobs->SetXTitle(convertVarName("ElPt[0]"));
	H_nsigobs->SetYTitle(Form("Events / %2.0f GeV", fBinWidthScale));

	H_nFpred->Add(H_nfppred);
	H_nFpred->Add(H_nffpred);
	H_nFpred->SetXTitle(convertVarName("ElPt[0]"));
	H_nFpred->SetYTitle(Form("Events / %2.0f GeV", fBinWidthScale));
	H_nt2obs->SetXTitle(convertVarName("ElPt[0]"));
	H_nt2obs->SetYTitle(Form("Events / %2.0f GeV", fBinWidthScale));
	H_nt2obsSM->SetXTitle(convertVarName("ElPt[0]"));
	H_nt2obsSM->SetYTitle(Form("Events / %2.0f GeV", fBinWidthScale));

	// Normalize to binwidth
	H_nsigpred = normHistBW(H_nsigpred, fBinWidthScale);
	H_nsigobs  = normHistBW(H_nsigobs,  fBinWidthScale);
	H_nfppred  = normHistBW(H_nfppred,  fBinWidthScale);
	H_nffpred  = normHistBW(H_nffpred,  fBinWidthScale);
	H_nFpred   = normHistBW(H_nFpred,   fBinWidthScale);
	H_nt2obs   = normHistBW(H_nt2obs,   fBinWidthScale);
	H_nt2obsSM = normHistBW(H_nt2obsSM, fBinWidthScale);

	H_nt2obs->SetFillColor(kBlue);
	H_nt2obs->SetLineColor(kBlue);
	H_nt2obs->SetFillStyle(3004);
	H_nt2obs->SetLineWidth(2);

	H_nsigpred->SetFillColor(8);
	H_nsigpred->SetMarkerColor(8);
	H_nsigpred->SetMarkerStyle(20);
	H_nsigpred->SetLineColor(8);
	H_nsigpred->SetLineWidth(2);

	H_nfppred->SetFillColor(kRed);
	H_nfppred->SetMarkerColor(kRed);
	H_nfppred->SetMarkerStyle(20);
	H_nfppred->SetLineColor(kRed);
	H_nfppred->SetLineWidth(2);

	H_nffpred->SetFillColor(13);
	H_nffpred->SetMarkerColor(13);
	H_nffpred->SetLineColor(13);
	H_nffpred->SetMarkerStyle(20);
	H_nffpred->SetLineWidth(2);

	plotOverlay4H(H_nt2obs, "N_{ t2}", H_nsigpred, "N_{ pp}" , H_nfppred, "N_{ f p}", H_nffpred, "N_{ f f}");

	H_nsigpred->SetMarkerColor(kRed);
	H_nFpred->SetMarkerColor(kRed);
	H_nFpred->SetMarkerStyle(20);
		
	plotPredOverlay2HWithRatio(H_nsigobs, "Observed N_{sig}", H_nsigpred, "Predicted N_{sig}");
	plotPredOverlay2HWithRatio(H_nt2obs, "Observed N_{t2}", H_nFpred, "Predicted Fakes");
	plotPredOverlay3HWithRatio(H_nFpred, "Predicted Fakes", H_nt2obs,  "Obs. N_{t2} (QCD, t#bar{t}+jets, V+jets, LM0)", H_nt2obsSM, "Obs. N_{t2} (QCD, t#bar{t}+jets, V+jets)", false, false);
	fOutputSubDir = "";
}
void MuonPlotter::makeSSElMuPredictionPlots(vector<int> samples, bool output){
	fOutputSubDir = "ElMuPredictions";
	// Need filled electron and muon ratios before calling this function!

	// Prediction: /////////////////////////////////////////////////////////////////////
	TH1D *H_npppred = new TH1D("ElMuNpppred", "Predicted N_pp in MuPt bins",        getNPt2Bins(Muon),  getPt2Bins(Muon));
	TH1D *H_nfppred = new TH1D("ElMuNfppred", "Predicted N_fp in MuPt bins",        getNPt2Bins(Muon),  getPt2Bins(Muon));
	TH1D *H_npfpred = new TH1D("ElMuNpfpred", "Predicted N_pf in MuPt bins",        getNPt2Bins(Muon),  getPt2Bins(Muon));
	TH1D *H_nffpred = new TH1D("ElMuNffpred", "Predicted N_ff in MuPt bins",        getNPt2Bins(Muon),  getPt2Bins(Muon));
	TH1D *H_nFpred  = new TH1D("ElMuNFpred",  "Total predicted fakes in MuPt bins", getNPt2Bins(Muon),  getPt2Bins(Muon));
	for(size_t i = 0; i < samples.size(); ++i){
		Sample *S = &fSamples[samples[i]];
		float scale = fLumiNorm/S->lumi;

		Channel *C = &S->region[Signal].em;
		if(fVerbose > 2) cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
		if(fVerbose > 2) cout << "Calling FPRatios for " << S->sname << endl;
		vector<TH1D*> prediction = ElMuFPPrediction(fH2D_MufRatio, fH2D_MupRatio, fH2D_ElfRatio, fH2D_ElpRatio, C->nt20_pt, C->nt10_pt, C->nt01_pt, C->nt00_pt, output);
		H_npppred->Add(prediction[0], scale);
		H_nfppred->Add(prediction[1], scale);
		H_npfpred->Add(prediction[2], scale);
		H_nffpred->Add(prediction[3], scale);
	}

	// Observation: ////////////////////////////////////////////////////////////////////
	TH1D *H_nsigobs  = new TH1D("Nsigobs", "Observed N_sig in MuPt bins",  getNPt2Bins(Muon),  getPt2Bins(Muon));
	vector<int> lm0sample; lm0sample.push_back(LM0);
	NObs(EMu, H_nsigobs, lm0sample, &MuonPlotter::isGenMatchedSUSYEMuEvent);

	TH1D *H_nt2obs  = new TH1D("Nt2obs", "Observed N_t2 in MuPt bins",  getNPt2Bins(Muon),  getPt2Bins(Muon));
	NObs(EMu, H_nt2obs, samples);

	TH1D *H_nt2obsSM = new TH1D("Nt2obsSM", "Observed N_t2 in MuPt bins, ttbar only",  getNPt2Bins(Muon),  getPt2Bins(Muon));
	NObs(EMu, H_nt2obsSM, fMCBG);

	// Output
	H_nsigobs->SetXTitle(convertVarName("MuPt[0]"));
	H_nsigobs->SetYTitle(Form("Events / %2.0f GeV", fBinWidthScale));

	H_nFpred->Add(H_nfppred);
	H_nFpred->Add(H_npfpred);
	H_nFpred->Add(H_nffpred);
	H_nFpred->SetXTitle(convertVarName("MuPt[0]"));
	H_nFpred->SetYTitle(Form("Events / %2.0f GeV", fBinWidthScale));
	H_nt2obs->SetXTitle(convertVarName("MuPt[0]"));
	H_nt2obs->SetYTitle(Form("Events / %2.0f GeV", fBinWidthScale));
	H_nt2obsSM->SetXTitle(convertVarName("MuPt[0]"));
	H_nt2obsSM->SetYTitle(Form("Events / %2.0f GeV", fBinWidthScale));

	// Normalize to binwidth
	H_npppred  = normHistBW(H_npppred,  fBinWidthScale);
	H_nsigobs  = normHistBW(H_nsigobs,  fBinWidthScale);
	H_nfppred  = normHistBW(H_nfppred,  fBinWidthScale);
	H_npfpred  = normHistBW(H_npfpred,  fBinWidthScale);
	H_nffpred  = normHistBW(H_nffpred,  fBinWidthScale);
	H_nFpred   = normHistBW(H_nFpred,   fBinWidthScale);
	H_nt2obs   = normHistBW(H_nt2obs,   fBinWidthScale);
	H_nt2obsSM = normHistBW(H_nt2obsSM, fBinWidthScale);

	H_nt2obs->SetFillColor(kBlue);
	H_nt2obs->SetLineColor(kBlue);
	H_nt2obs->SetFillStyle(3004);
	H_nt2obs->SetLineWidth(2);

	int ppcolor = 8;
	H_npppred->SetFillColor(  ppcolor);
	H_npppred->SetMarkerColor(ppcolor);
	H_npppred->SetLineColor(  ppcolor);
	H_npppred->SetMarkerStyle(20);
	H_npppred->SetLineWidth(2);

	int fpcolor = 2;
	H_nfppred->SetFillColor(  fpcolor);
	H_nfppred->SetMarkerColor(fpcolor);
	H_nfppred->SetLineColor(  fpcolor);
	H_nfppred->SetMarkerStyle(20);
	H_nfppred->SetLineWidth(2);

	int pfcolor = 51;
	H_npfpred->SetFillColor(  pfcolor);
	H_npfpred->SetMarkerColor(pfcolor);
	H_npfpred->SetLineColor(  pfcolor);
	H_npfpred->SetMarkerStyle(20);
	H_npfpred->SetLineWidth(2);

	int ffcolor = 13;
	H_nffpred->SetFillColor(  ffcolor);
	H_nffpred->SetMarkerColor(ffcolor);
	H_nffpred->SetLineColor(  ffcolor);
	H_nffpred->SetMarkerStyle(20);
	H_nffpred->SetLineWidth(2);

	plotOverlay5H(H_nt2obs, "N_{ t2}", H_npppred, "N_{ pp}" , H_nfppred, "N_{ f p}", H_npfpred, "N_{ p f}", H_nffpred, "N_{ f f}");

	H_npppred->SetMarkerColor(kRed);
	H_nFpred->SetMarkerColor(kRed);
	H_nFpred->SetMarkerStyle(20);
	plotPredOverlay2HWithRatio(H_nsigobs, "Observed N_{sig}", H_npppred, "Predicted N_{sig}");
	plotPredOverlay2HWithRatio(H_nt2obs,  "Observed N_{t2}",  H_nFpred,  "Predicted Fakes");
	plotPredOverlay3HWithRatio(H_nFpred,  "Predicted Fakes",  H_nt2obs,  "Obs. N_{t2} (QCD, t#bar{t}+jets, V+jets, LM0)", H_nt2obsSM, "Obs. N_{t2} (QCD, t#bar{t}+jets, V+jets)", false, false);
	fOutputSubDir = "";
}

//____________________________________________________________________________
void MuonPlotter::NObs(gChannel chan, TH1D *&hist, vector<int> samples, bool(MuonPlotter::*eventSelector)()){
/* This fills a histogram with the pt of the first muon for a given selection */
	hist->Sumw2();
	for(size_t i = 0; i < samples.size(); ++i){
		Sample *S = &fSamples[samples[i]];
		float scale = fLumiNorm / S->lumi;

		TTree *tree = S->tree;
		tree->ResetBranchAddresses();
		Init(tree);
		if (fChain == 0) return;
		Long64_t nentries = fChain->GetEntriesFast();
		Long64_t nbytes = 0, nb = 0;
		for (Long64_t jentry=0; jentry<nentries;jentry++) {
			Long64_t ientry = LoadTree(jentry);
			if (ientry < 0) break;
			nb = fChain->GetEntry(jentry);   nbytes += nb;
			printProgress(jentry, nentries, S->name);

			if((*this.*eventSelector)() == false) continue;

			if((chan == Muon || chan == EMu) && NMus > 0) hist->Fill(MuPt[0], scale);
			if( chan == Electron             && NEls > 0) hist->Fill(ElPt[0], scale);
		}
	}	
}
void MuonPlotter::NObs(gChannel chan, TH1D *&hist, vector<int> samples, gRegion reg){
	hist->Sumw2();
	for(size_t i = 0; i < samples.size(); ++i){
		Sample *S = &fSamples[samples[i]];
		float scale = fLumiNorm / S->lumi;
		Channel *C;
		if(chan == Muon)     C = &S->region[reg].mm;
		if(chan == EMu)      C = &S->region[reg].em;
		if(chan == Electron) C = &S->region[reg].ee;
		hist->Add(C->nt20_pt->ProjectionX(), scale);
	}	
}
void MuonPlotter::NObs(gChannel chan, THStack *&stack, vector<int> samples, gRegion reg){
	// hist->Sumw2();
	const int nsamples = samples.size();
	TH1D *hnt2[nsamples];
	for(size_t i = 0; i < nsamples; ++i){
		Sample *S = &fSamples[samples[i]];
		float scale = fLumiNorm / S->lumi;
		Channel *C;
		if(chan == Muon)     C = &S->region[reg].mm;
		if(chan == EMu)      C = &S->region[reg].em;
		if(chan == Electron) C = &S->region[reg].ee;
		hnt2[i] = (TH1D*)C->nt20_pt->ProjectionX()->Clone(Form("Nt2_%s", S->sname.Data())); // not sure if i need the clone here...
		hnt2[i]->SetFillColor(S->color);
		hnt2[i]->Sumw2();
		hnt2[i]->Scale(scale);
		stack->Add(hnt2[i]);
	}	
}

//____________________________________________________________________________
void MuonPlotter::makeIntPrediction(TString filename){
	ofstream OUT(filename.Data(), ios::trunc);

	vector<int> musamples;
	vector<int> elsamples;
	vector<int> emusamples;

	if(fSelectionSwitch == 0){
		musamples = fMuData;
		elsamples = fEGData;
		emusamples.push_back(MuA);
		emusamples.push_back(MuB);
		emusamples.push_back(EGA);
		emusamples.push_back(EGB);
	}
	if(fSelectionSwitch == 1){
		musamples  = fJMData;
		elsamples  = fJMData;
		emusamples = fJMData;
	}

	OUT << "/////////////////////////////////////////////////////////////////////////////" << endl;
	OUT << " Producing integrated predictions" << endl << endl;

	///////////////////////////////////////////////////////////////////////////////////
	// RATIOS /////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////
	float mufratio_data(0.),  mufratio_data_elow(0.), mufratio_data_eup(0.);
	float mupratio_data(0.),  mupratio_data_elow(0.), mupratio_data_eup(0.);
	float elfratio_data(0.),  elfratio_data_elow(0.), elfratio_data_eup(0.);
	float elpratio_data(0.),  elpratio_data_elow(0.), elpratio_data_eup(0.);
	float mufratio_allmc(0.), mufratio_allmc_e(0.);
	float mupratio_allmc(0.), mupratio_allmc_e(0.);
	float elfratio_allmc(0.), elfratio_allmc_e(0.);
	float elpratio_allmc(0.), elpratio_allmc_e(0.);

	calculateRatio(fJMData, Muon, SigSup, mufratio_data, mufratio_data_elow);
	calculateRatio(fMuData, Muon, ZDecay, mupratio_data, mupratio_data_elow);
	// calculateRatio(fJMData, Muon, SigSup, mufratio_data, mufratio_data_elow , mufratio_data_eup);
	// calculateRatio(fMuData, Muon, ZDecay, mupratio_data, mupratio_data_elow , mupratio_data_eup);

	calculateRatio(elsamples, Electron, SigSup, elfratio_data, elfratio_data_elow);
	calculateRatio(fEGData,   Electron, ZDecay, elpratio_data, elpratio_data_elow);
	// calculateRatio(elsamples, Electron, 1, elfratio_data, elfratio_data_elow, elfratio_data_eup);
	// calculateRatio(fEGData,  Electron, 2, elpratio_data, elpratio_data_elow, elpratio_data_eup);

	calculateRatio(fMCBG, Muon,     SigSup, mufratio_allmc, mufratio_allmc_e);
	calculateRatio(fMCBG, Muon,     ZDecay, mupratio_allmc, mupratio_allmc_e);
	calculateRatio(fMCBG, Electron, SigSup, elfratio_allmc, elfratio_allmc_e);
	calculateRatio(fMCBG, Electron, ZDecay, elpratio_allmc, elpratio_allmc_e);

	///////////////////////////////////////////////////////////////////////////////////
	// SYSTEMATICS ////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////
	// // TTbar genmatch ratios for muons:
	// const int nptbins = 1; // Dummy binning
	// const double ptbins[nptbins+1] = {5., 1000.};
	// const int netabins = 1;
	// const double etabins[netabins+1] = {-2.5, 2.5};
	// 
	// TH2D *H_fratio_ttbar = fillRatio(TTbar, 0, &MuonPlotter::isGoodMuEvent, &MuonPlotter::isFakeTTbarMuon,   nptbins, ptbins, netabins, etabins);
	// TH2D *H_pratio_ttbar = fillRatio(TTbar, 0, &MuonPlotter::isGoodMuEvent, &MuonPlotter::isPromptTTbarMuon, nptbins, ptbins, netabins, etabins);
	// H_fratio_ttbar->SetName("MufRatioTTbar");
	// H_pratio_ttbar->SetName("MupRatioTTbar");
	// 
	// float mufratio_ttbar   = H_fratio_ttbar->GetBinContent(1,1);
	// float mufratio_ttbar_e = H_fratio_ttbar->GetBinError(1,1);
	// float mupratio_ttbar   = H_pratio_ttbar->GetBinContent(1,1);
	// float mupratio_ttbar_e = H_pratio_ttbar->GetBinError(1,1);

	// // Add deviation as systematic (only muons)
	// float deviation(0.), olderror(0.), newerror(0.);
	// deviation = fabs(mufratio_allmc - mufratio_ttbar);
	// // Add to mc ratios
	// olderror = mufratio_allmc_e;
	// newerror = sqrt(olderror*olderror + deviation*deviation);
	// mufratio_allmc_e = newerror;
	// // Add to data ratios
	// olderror = mufratio_data_e;
	// newerror = sqrt(olderror*olderror + deviation*deviation);
	// mufratio_data_e = newerror;
	// 
	// deviation  = fabs(mupratio_allmc - mupratio_ttbar);
	// // Add to mc ratios
	// olderror = mupratio_allmc_e;
	// newerror = sqrt(olderror*olderror + deviation*deviation);
	// mupratio_allmc_e = newerror;
	// // Add to data ratios
	// olderror = mupratio_data_e;
	// newerror = sqrt(olderror*olderror + deviation*deviation);
	// mupratio_data_e = newerror;
	// 
	// // Add hardcoded systematic error on electron ratios:
	// float adderror = 0.05; // add this to e f ratio
	// olderror = elfratio_data_e;
	// newerror = sqrt(olderror*olderror + adderror*adderror);
	// elfratio_data_e = newerror;
	// 
	// adderror = 0.01; // add this to e p ratio
	// olderror = elpratio_data_e;
	// newerror = sqrt(olderror*olderror + adderror*adderror);
	// elpratio_data_e = newerror;

	///////////////////////////////////////////////////////////////////////////////////
	// OBSERVATIONS ///////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////
	float nt2_mumu(0.),    nt10_mumu(0.),    nt0_mumu(0.);
	float nt2_mumu_e2(0.), nt10_mumu_e2(0.), nt0_mumu_e2(0.);
	float nt2_emu(0.),     nt10_emu(0.),     nt01_emu(0.),    nt0_emu(0.);
	float nt2_emu_e2(0.),  nt10_emu_e2(0.),  nt01_emu_e2(0.), nt0_emu_e2(0.);	
	float nt2_ee(0.),      nt10_ee(0.),      nt0_ee(0.);
	float nt2_ee_e2(0.),   nt10_ee_e2(0.),   nt0_ee_e2(0.);
	

	for(size_t i = 0; i < musamples.size(); ++i){
		Sample S = fSamples[musamples[i]];
		nt2_mumu     += S.numbers[Muon].nt2;
		nt10_mumu    += S.numbers[Muon].nt10;
		nt0_mumu     += S.numbers[Muon].nt0;
		nt2_mumu_e2  += S.numbers[Muon].nt2;
		nt10_mumu_e2 += S.numbers[Muon].nt10;
		nt0_mumu_e2  += S.numbers[Muon].nt0;
	}		
	for(size_t i = 0; i < emusamples.size(); ++i){
		Sample S = fSamples[emusamples[i]];
		nt2_emu     += S.numbers[EMu].nt2;
		nt10_emu    += S.numbers[EMu].nt10;
		nt01_emu    += S.numbers[EMu].nt01;
		nt0_emu     += S.numbers[EMu].nt0;
		nt2_emu_e2  += S.numbers[EMu].nt2;
		nt10_emu_e2 += S.numbers[EMu].nt10;
		nt01_emu_e2 += S.numbers[EMu].nt01;
		nt0_emu_e2  += S.numbers[EMu].nt0;
	}		
	for(size_t i = 0; i < elsamples.size(); ++i){
		Sample S = fSamples[elsamples[i]];
		nt2_ee     += S.numbers[Electron].nt2;
		nt10_ee    += S.numbers[Electron].nt10;
		nt0_ee     += S.numbers[Electron].nt0;
		nt2_ee_e2  += S.numbers[Electron].nt2;
		nt10_ee_e2 += S.numbers[Electron].nt10;
		nt0_ee_e2  += S.numbers[Electron].nt0;
	}		
	
	///////////////////////////////////////////////////////////////////////////////////
	// PREDICTIONS ////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////
	FPRatios *fMuMuFPRatios_lo = new FPRatios();
	// FPRatios *fMuMuFPRatios_up = new FPRatios();
	FPRatios *fEMuFPRatios     = new FPRatios();
	FPRatios *fEEFPRatios      = new FPRatios();

	fMuMuFPRatios_lo->SetVerbose(fVerbose);
	// fMuMuFPRatios_up->SetVerbose(fVerbose);
	fEMuFPRatios ->SetVerbose(fVerbose);
	fEEFPRatios  ->SetVerbose(fVerbose);

	fEEFPRatios  ->SetElFratios(elfratio_data, elfratio_data_elow);
	fEEFPRatios  ->SetElPratios(elpratio_data, elpratio_data_elow);

	fEMuFPRatios ->SetFratios(elfratio_data, elfratio_data_elow, mufratio_data, mufratio_data_elow);
	fEMuFPRatios ->SetPratios(elpratio_data, elpratio_data_elow, mupratio_data, mupratio_data_elow);


	// MuMu Channel
	fMuMuFPRatios_lo->SetMuFratios(mufratio_data, mufratio_data_elow);
	fMuMuFPRatios_lo->SetMuPratios(mupratio_data, mupratio_data_elow);
	// fMuMuFPRatios_up->SetMuFratios(mufratio_data, mufratio_data_eup);
	// fMuMuFPRatios_up->SetMuPratios(mupratio_data, mupratio_data_eup);

	vector<double> nt_mumu;
	nt_mumu.push_back(nt0_mumu);
	nt_mumu.push_back(nt10_mumu); // mu passes
	nt_mumu.push_back(nt2_mumu);
	
	fMuMuFPRatios_lo->NevtTopol(0, 2, nt_mumu);
	// fMuMuFPRatios_up->NevtTopol(0, 2, nt_mumu);

	vector<double> vpt, veta;
	vpt.push_back(30.); vpt.push_back(30.); // Fake pts and etas (first electron then muon)
	veta.push_back(0.); veta.push_back(0.);

	vector<double> MuMu_Nev      = fMuMuFPRatios_lo->NevtPass(vpt, veta);
	vector<double> MuMu_Estat_lo = fMuMuFPRatios_lo->NevtPassErrStat();
	vector<double> MuMu_Esyst_lo = fMuMuFPRatios_lo->NevtPassErrSyst();

	// vector<double> MuMu_Nev_up   = fMuMuFPRatios_up->NevtPass(vpt, veta);
	// vector<double> MuMu_Estat_up = fMuMuFPRatios_up->NevtPassErrStat();
	// vector<double> MuMu_Esyst_up = fMuMuFPRatios_up->NevtPassErrSyst();

	
	// EMu Channel
	vector<double> nt_emu;
	nt_emu.push_back(nt0_emu);
	nt_emu.push_back(nt01_emu); // e passes
	nt_emu.push_back(nt10_emu); // mu passes
	nt_emu.push_back(nt2_emu);
	
	fEMuFPRatios->NevtTopol(1, 1, nt_emu);

	vector<double> EMu_nevFP      = fEMuFPRatios->NevtPass(vpt, veta);
	vector<double> EMu_nevFPEstat = fEMuFPRatios->NevtPassErrStat();
	vector<double> EMu_nevFPEsyst = fEMuFPRatios->NevtPassErrSyst();
	
	// EE Channel
	vector<double> nt_ee;
	nt_ee.push_back(nt0_ee);
	nt_ee.push_back(nt10_ee); // el passes
	nt_ee.push_back(nt2_ee);
	
	fEEFPRatios->NevtTopol(2, 0, nt_ee);

	vector<double> EE_nevFP      = fEEFPRatios->NevtPass(vpt, veta);
	vector<double> EE_nevFPEstat = fEEFPRatios->NevtPassErrStat();
	vector<double> EE_nevFPEsyst = fEEFPRatios->NevtPassErrSyst();
	
	///////////////////////////////////////////////////////////////////////////////////
	// PRINTOUT ///////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////
	OUT << "--------------------------------------------------------------------------------------------------" << endl;
	OUT << "  RATIOS  ||     Mu-fRatio      |     Mu-pRatio      ||     El-fRatio      |     El-pRatio      ||" << endl;
	OUT << "--------------------------------------------------------------------------------------------------" << endl;
	OUT << setw(10) << " allMC    ||";
	OUT << setw(7)  << setprecision(2) << mufratio_allmc << " +/- " << setw(7) << setprecision(2) << mufratio_allmc_e << " |";
	OUT << setw(7)  << setprecision(2) << mupratio_allmc << " +/- " << setw(7) << setprecision(2) << mupratio_allmc_e << " ||";
	OUT << setw(7)  << setprecision(2) << elfratio_allmc << " +/- " << setw(7) << setprecision(2) << elfratio_allmc_e << " |";
	OUT << setw(7)  << setprecision(2) << elpratio_allmc << " +/- " << setw(7) << setprecision(2) << elpratio_allmc_e << " ||";
	OUT << endl;
	OUT << setw(10) << " data     ||";
	OUT << setw(7)  << setprecision(2) << mufratio_data  << " +/- " << setw(7) << setprecision(2) << mufratio_data_elow  << " |";
	OUT << setw(7)  << setprecision(2) << mupratio_data  << " +/- " << setw(7) << setprecision(2) << mupratio_data_elow  << " ||";
	OUT << setw(7)  << setprecision(2) << elfratio_data  << " +/- " << setw(7) << setprecision(2) << elfratio_data_elow  << " |";
	OUT << setw(7)  << setprecision(2) << elpratio_data  << " +/- " << setw(7) << setprecision(2) << elpratio_data_elow  << " ||";
	OUT << endl;
	OUT << "--------------------------------------------------------------------------------------------------" << endl << endl;
	OUT << "-------------------------------------------------------------------------------------------------------------------" << endl;
	OUT << "          ||            Mu/Mu            ||                   E/Mu                ||             E/E             ||" << endl;
	OUT << "  YIELDS  ||   Nt2   |   Nt1   |   Nt0   ||   Nt2   |   Nt10  |   Nt01  |   Nt0   ||   Nt2   |   Nt1   |   Nt0   ||" << endl;
	OUT << "-------------------------------------------------------------------------------------------------------------------" << endl;
	
	float nt2sum_mumu(0.), nt10sum_mumu(0.), nt0sum_mumu(0.);
	float nt2sum_emu(0.), nt10sum_emu(0.), nt01sum_emu(0.), nt0sum_emu(0.);
	float nt2sum_ee(0.), nt10sum_ee(0.), nt0sum_ee(0.);
	for(size_t i = 0; i < fMCBG.size(); ++i){
		int index = fMCBG[i];
		float scale = fLumiNorm / fSamples[index].lumi;
		nt2sum_mumu  += scale*fSamples[index].numbers[Muon]    .nt2;
		nt10sum_mumu += scale*fSamples[index].numbers[Muon]    .nt10;
		nt0sum_mumu  += scale*fSamples[index].numbers[Muon]    .nt0;
		nt2sum_emu   += scale*fSamples[index].numbers[EMu]     .nt2;
		nt10sum_emu  += scale*fSamples[index].numbers[EMu]     .nt10;
		nt01sum_emu  += scale*fSamples[index].numbers[EMu]     .nt01;
		nt0sum_emu   += scale*fSamples[index].numbers[EMu]     .nt0;
		nt2sum_ee    += scale*fSamples[index].numbers[Electron].nt2;
		nt10sum_ee   += scale*fSamples[index].numbers[Electron].nt10;
		nt0sum_ee    += scale*fSamples[index].numbers[Electron].nt0;

		OUT << setw(9) << fSamples[index].sname << " || ";
		OUT << setw(7)  << setprecision(2) << scale*fSamples[index].numbers[Muon]    .nt2  << " | ";
		OUT << setw(7)  << setprecision(2) << scale*fSamples[index].numbers[Muon]    .nt10 << " | ";
		OUT << setw(7)  << setprecision(2) << scale*fSamples[index].numbers[Muon]    .nt0  << " || ";
		OUT << setw(7)  << setprecision(2) << scale*fSamples[index].numbers[EMu]     .nt2  << " | ";
		OUT << setw(7)  << setprecision(2) << scale*fSamples[index].numbers[EMu]     .nt10 << " | ";
		OUT << setw(7)  << setprecision(2) << scale*fSamples[index].numbers[EMu]     .nt01 << " | ";
		OUT << setw(7)  << setprecision(2) << scale*fSamples[index].numbers[EMu]     .nt0  << " || ";
		OUT << setw(7)  << setprecision(2) << scale*fSamples[index].numbers[Electron].nt2  << " | ";
		OUT << setw(7)  << setprecision(2) << scale*fSamples[index].numbers[Electron].nt10 << " | ";
		OUT << setw(7)  << setprecision(2) << scale*fSamples[index].numbers[Electron].nt0  << " || ";
		OUT << endl;
	}	
	OUT << "-------------------------------------------------------------------------------------------------------------------" << endl;
	OUT << setw(9) << "MC sum" << " || ";
	OUT << setw(7) << setprecision(2) << nt2sum_mumu  << " | ";
	OUT << setw(7) << setprecision(2) << nt10sum_mumu << " | ";
	OUT << setw(7) << setprecision(2) << nt0sum_mumu  << " || ";
	OUT << setw(7) << setprecision(2) << nt2sum_emu   << " | ";
	OUT << setw(7) << setprecision(2) << nt10sum_emu  << " | ";
	OUT << setw(7) << setprecision(2) << nt01sum_emu  << " | ";
	OUT << setw(7) << setprecision(2) << nt0sum_emu   << " || ";
	OUT << setw(7) << setprecision(2) << nt2sum_ee    << " | ";
	OUT << setw(7) << setprecision(2) << nt10sum_ee   << " | ";
	OUT << setw(7) << setprecision(2) << nt0sum_ee    << " || ";
	OUT << endl;
	OUT << "-------------------------------------------------------------------------------------------------------------------" << endl;
	OUT << setw(9) << fSamples[LM0].sname << " || ";
	OUT << setw(7)  << setprecision(2) << fLumiNorm / fSamples[LM0].lumi * fSamples[LM0].numbers[Muon]    .nt2  << " | ";
	OUT << setw(7)  << setprecision(2) << fLumiNorm / fSamples[LM0].lumi * fSamples[LM0].numbers[Muon]    .nt10 << " | ";
	OUT << setw(7)  << setprecision(2) << fLumiNorm / fSamples[LM0].lumi * fSamples[LM0].numbers[Muon]    .nt0  << " || ";
	OUT << setw(7)  << setprecision(2) << fLumiNorm / fSamples[LM0].lumi * fSamples[LM0].numbers[EMu]     .nt2  << " | ";
	OUT << setw(7)  << setprecision(2) << fLumiNorm / fSamples[LM0].lumi * fSamples[LM0].numbers[EMu]     .nt10 << " | ";
	OUT << setw(7)  << setprecision(2) << fLumiNorm / fSamples[LM0].lumi * fSamples[LM0].numbers[EMu]     .nt01 << " | ";
	OUT << setw(7)  << setprecision(2) << fLumiNorm / fSamples[LM0].lumi * fSamples[LM0].numbers[EMu]     .nt0  << " || ";
	OUT << setw(7)  << setprecision(2) << fLumiNorm / fSamples[LM0].lumi * fSamples[LM0].numbers[Electron].nt2  << " | ";
	OUT << setw(7)  << setprecision(2) << fLumiNorm / fSamples[LM0].lumi * fSamples[LM0].numbers[Electron].nt10 << " | ";
	OUT << setw(7)  << setprecision(2) << fLumiNorm / fSamples[LM0].lumi * fSamples[LM0].numbers[Electron].nt0  << " || ";
	OUT << endl;
	OUT << "-------------------------------------------------------------------------------------------------------------------" << endl;
	OUT << setw(9) << "data"  << " || ";
	OUT << setw(7) << setprecision(2) << nt2_mumu  << " | ";
	OUT << setw(7) << setprecision(2) << nt10_mumu << " | ";
	OUT << setw(7) << setprecision(2) << nt0_mumu  << " || ";
	OUT << setw(7) << setprecision(2) << nt2_emu   << " | ";
	OUT << setw(7) << setprecision(2) << nt10_emu  << " | ";
	OUT << setw(7) << setprecision(2) << nt01_emu  << " | ";
	OUT << setw(7) << setprecision(2) << nt0_emu   << " || ";
	OUT << setw(7) << setprecision(2) << nt2_ee    << " | ";
	OUT << setw(7) << setprecision(2) << nt10_ee   << " | ";
	OUT << setw(7) << setprecision(2) << nt0_ee    << " || ";
	OUT << endl;
	OUT << "-------------------------------------------------------------------------------------------------------------------" << endl;
	OUT << endl;
	OUT << "  PREDICTIONS" << endl;
	OUT << "--------------------------------------------------------------" << endl;
	OUT << " Mu/Mu Channel:" << endl;
	// OUT << "  Npp:           " << setw(7) << setprecision(3) << MuMu_Nev[0];
	// OUT << " +/- "             << setw(7) << setprecision(3) << MuMu_Estat_lo[0] << " (stat)";
	// OUT << " + "               << setw(7) << setprecision(3) << MuMu_Esyst_up[0];
	// OUT << " - "               << setw(7) << setprecision(3) << MuMu_Esyst_lo[0] << " (syst)" << endl;
	// OUT << "  Nfp:           " << setw(7) << setprecision(3) << MuMu_Nev[1];
	// OUT << " +/- "             << setw(7) << setprecision(3) << MuMu_Estat_lo[1] << " (stat)";
	// OUT << " + "               << setw(7) << setprecision(3) << MuMu_Esyst_up[1];
	// OUT << " - "               << setw(7) << setprecision(3) << MuMu_Esyst_lo[1] << " (syst)" << endl;
	// OUT << "  Nff:           " << setw(7) << setprecision(3) << MuMu_Nev[2];
	// OUT << " +/- "             << setw(7) << setprecision(3) << MuMu_Estat_lo[2] << " (stat)";
	// OUT << " + "               << setw(7) << setprecision(3) << MuMu_Esyst_up[2];
	// OUT << " - "               << setw(7) << setprecision(3) << MuMu_Esyst_lo[2] << " (syst)" << endl;
	// OUT << "  Total fakes:   " << setw(7) << setprecision(3) << MuMu_Nev[1]+MuMu_Nev[2];
	// OUT << " +/- "             << setw(7) << setprecision(3) << TMath::Sqrt(MuMu_Estat_lo[1]*MuMu_Estat_lo[1] + MuMu_Estat_lo[2]*MuMu_Estat_lo[2]) << " (stat)";
	// OUT << " + "               << setw(7) << setprecision(3) << TMath::Sqrt(MuMu_Esyst_up[1]*MuMu_Esyst_up[1] + MuMu_Esyst_up[2]*MuMu_Esyst_up[2]);
	// OUT << " - "               << setw(7) << setprecision(3) << TMath::Sqrt(MuMu_Esyst_lo[1]*MuMu_Esyst_lo[1] + MuMu_Esyst_lo[2]*MuMu_Esyst_lo[2]) << " (syst)" << endl;
	OUT << "  Npp:           " << setw(7) << setprecision(3) << MuMu_Nev[0];
	OUT << " +/- "             << setw(7) << setprecision(3) << MuMu_Estat_lo[0];
	OUT << " (stat) +/- "      << setw(7) << setprecision(3) << MuMu_Esyst_lo[0] << " (syst)" << endl;
	OUT << "  Nfp:           " << setw(7) << setprecision(3) << MuMu_Nev[1];
	OUT << " +/- "             << setw(7) << setprecision(3) << MuMu_Estat_lo[1];
	OUT << " (stat) +/- "      << setw(7) << setprecision(3) << MuMu_Esyst_lo[1] << " (syst)" << endl;
	OUT << "  Nff:           " << setw(7) << setprecision(3) << MuMu_Nev[2];
	OUT << " +/- "             << setw(7) << setprecision(3) << MuMu_Estat_lo[2];
	OUT << " (stat) +/- "      << setw(7) << setprecision(3) << MuMu_Esyst_lo[2] << " (syst)" << endl;
	OUT << "  Total fakes:   " << setw(7) << setprecision(3) << MuMu_Nev[1]+MuMu_Nev[2];
	OUT << " +/- "             << setw(7) << setprecision(3) << TMath::Sqrt(MuMu_Estat_lo[1]*MuMu_Estat_lo[1] + MuMu_Estat_lo[2]*MuMu_Estat_lo[2]);
	OUT << " (stat) +/- "      << setw(7) << setprecision(3) << TMath::Sqrt(MuMu_Esyst_lo[1]*MuMu_Esyst_lo[1] + MuMu_Esyst_lo[2]*MuMu_Esyst_lo[2]) << " (syst)" << endl;
	OUT << "--------------------------------------------------------------" << endl;
	OUT << " E/Mu Channel:" << endl;
	OUT << "  Npp:           " << setw(7) << setprecision(3) << EMu_nevFP[0];
	OUT << " +/- "             << setw(7) << setprecision(3) << EMu_nevFPEstat[0];
	OUT << " (stat) +/- "      << setw(7) << setprecision(3) << EMu_nevFPEsyst[0] << " (syst)" << endl;
	OUT << "  Nfp:           " << setw(7) << setprecision(3) << EMu_nevFP[1];
	OUT << " +/- "             << setw(7) << setprecision(3) << EMu_nevFPEstat[1];
	OUT << " (stat) +/- "      << setw(7) << setprecision(3) << EMu_nevFPEsyst[1] << " (syst)" << endl;
	OUT << "  Npf:           " << setw(7) << setprecision(3) << EMu_nevFP[2];
	OUT << " +/- "             << setw(7) << setprecision(3) << EMu_nevFPEstat[2];
	OUT << " (stat) +/- "      << setw(7) << setprecision(3) << EMu_nevFPEsyst[2] << " (syst)" << endl;
	OUT << "  Nff:           " << setw(7) << setprecision(3) << EMu_nevFP[3];
	OUT << " +/- "             << setw(7) << setprecision(3) << EMu_nevFPEstat[3];
	OUT << " (stat) +/- "      << setw(7) << setprecision(3) << EMu_nevFPEsyst[3] << " (syst)" << endl;
	OUT << "  Total fakes:   " << setw(7) << setprecision(3) << EMu_nevFP[1]+EMu_nevFP[2]+EMu_nevFP[3];
	OUT << " +/- "             << setw(7) << setprecision(3) << TMath::Sqrt(EMu_nevFPEstat[1]*EMu_nevFPEstat[1] + EMu_nevFPEstat[2]*EMu_nevFPEstat[2] + EMu_nevFPEstat[3]*EMu_nevFPEstat[3]);
	OUT << " (stat) +/- "      << setw(7) << setprecision(3) << TMath::Sqrt(EMu_nevFPEsyst[1]*EMu_nevFPEsyst[1] + EMu_nevFPEsyst[2]*EMu_nevFPEsyst[2] + EMu_nevFPEsyst[3]*EMu_nevFPEsyst[3]) << " (syst)" << endl;
	OUT << "--------------------------------------------------------------" << endl;
	OUT << " E/E Channel:" << endl;
	OUT << "  Npp:           " << setw(7) << setprecision(3) << EE_nevFP[0];
	OUT << " +/- "             << setw(7) << setprecision(3) << EE_nevFPEstat[0];
	OUT << " (stat) +/- "      << setw(7) << setprecision(3) << EE_nevFPEsyst[0] << " (syst)" << endl;
	OUT << "  Nfp:           " << setw(7) << setprecision(3) << EE_nevFP[1];
	OUT << " +/- "             << setw(7) << setprecision(3) << EE_nevFPEstat[1];
	OUT << " (stat) +/- "      << setw(7) << setprecision(3) << EE_nevFPEsyst[1] << " (syst)" << endl;
	OUT << "  Nff:           " << setw(7) << setprecision(3) << EE_nevFP[2];
	OUT << " +/- "             << setw(7) << setprecision(3) << EE_nevFPEstat[2];
	OUT << " (stat) +/- "      << setw(7) << setprecision(3) << EE_nevFPEsyst[2] << " (syst)" << endl;
	OUT << "  Total fakes:   " << setw(7) << setprecision(3) << EE_nevFP[1]+EE_nevFP[2];
	OUT << " +/- "             << setw(7) << setprecision(3) << TMath::Sqrt(EE_nevFPEstat[1]*EE_nevFPEstat[1] + EE_nevFPEstat[2]*EE_nevFPEstat[2]);
	OUT << " (stat) +/- "      << setw(7) << setprecision(3) << TMath::Sqrt(EE_nevFPEsyst[1]*EE_nevFPEsyst[1] + EE_nevFPEsyst[2]*EE_nevFPEsyst[2]) << " (syst)" << endl;
	OUT << "--------------------------------------------------------------" << endl;
	OUT << "/////////////////////////////////////////////////////////////////////////////" << endl;

	OUT.close();
	delete fMuMuFPRatios_lo;
	// delete fMuMuFPRatios_up;
	delete fEMuFPRatios;
	delete fEEFPRatios;
}

//____________________________________________________________________________
vector<TH1D*> MuonPlotter::MuMuFPPrediction(TH2D* H_fratio, TH2D* H_pratio, TH2D* H_nt2mes, TH2D* H_nt1mes, TH2D* H_nt0mes, bool output){
	///////////////////////////////////////////////////////////////////////////
	// Note: Careful, this is only a workaround at the moment, and is only
	//       really valid for one single eta bin!
	//       In the future we should rewrite the interface to FPRatios and
	//       pass it the ratios directly to have full control
	///////////////////////////////////////////////////////////////////////////
	vector<TH1D*> res;
	const int nptbins = getNPt2Bins(Muon);

	TH1D *H_dummy    = new TH1D("MuDummyhist","Dummy",                       H_nt2mes->GetNbinsX(), H_nt2mes->GetXaxis()->GetXbins()->GetArray());
	TH1D *H_nsigpred = new TH1D("MuNsigpred", "Predicted N_sig in Pt1 bins", H_nt2mes->GetNbinsX(), H_nt2mes->GetXaxis()->GetXbins()->GetArray());
	TH1D *H_nfppred  = new TH1D("MuNfppred",  "Predicted N_fp in Pt1 bins",  H_nt2mes->GetNbinsX(), H_nt2mes->GetXaxis()->GetXbins()->GetArray());
	TH1D *H_nffpred  = new TH1D("MuNffpred",  "Predicted N_ff in Pt1 bins",  H_nt2mes->GetNbinsX(), H_nt2mes->GetXaxis()->GetXbins()->GetArray());	

	if(fVerbose > 2){
		cout << " Nt2 = " << H_nt2mes->GetEntries();
		cout << " Nt1 = " << H_nt1mes->GetEntries();
		cout << " Nt0 = " << H_nt0mes->GetEntries() << endl;
	}

	for(size_t i = 1; i <= nptbins; ++i){
		double pt1 = H_dummy->GetBinCenter(i);
		if(fVerbose > 2){
			cout << "=======================================" << endl;
			cout << "MuPt1 = " << pt1 << "  MufRatio = " << H_fratio->GetBinContent(i, 1) << "  MupRatio = " << H_pratio->GetBinContent(i, 1) << endl;
		}
		double eta1 = 0.0;
		double npppred(0.0), npppredEstat2(0.0), npppredEsyst2(0.0);
		double nfppred(0.0), nfppredEstat2(0.0), nfppredEsyst2(0.0);
		double nffpred(0.0), nffpredEstat2(0.0), nffpredEsyst2(0.0);
		for(size_t j = 1; j <= nptbins; ++j){
			if(fVerbose > 2) cout << " --------" << endl;
			double pt2 = H_dummy->GetBinCenter(j);

			double eta2 = 0.0;
			double nt2 = H_nt2mes->GetBinContent(i,j);
			double nt1 = H_nt1mes->GetBinContent(i,j);
			double nt0 = H_nt0mes->GetBinContent(i,j);

			if(fVerbose > 2) cout << "MuPt2 = " << pt2 << "  MufRatio = " << H_fratio->GetBinContent(j, 1) << "  MupRatio = " << H_pratio->GetBinContent(j, 1) << endl;
			if(fVerbose > 2) cout << "   nt2: " << nt2 << "  nt1: " << nt1 << "  nt0: " << nt0 << endl;

			FPRatios *fFPRatios = new FPRatios();
			fFPRatios->SetVerbose(fVerbose);
			fFPRatios->SetMuFratios(H_fratio);
			fFPRatios->SetMuPratios(H_pratio);
			vector<double> nt;
			nt.push_back(nt0); nt.push_back(nt1); nt.push_back(nt2);
			fFPRatios->NevtTopol(0, 2, nt);

			vector<double> vpt, veta;
			vpt.push_back(pt1); vpt.push_back(pt2);
			veta.push_back(eta1); veta.push_back(eta2);

			vector<double> nevFP = fFPRatios->NevtPass(vpt, veta);
			vector<double> nevFPEstat = fFPRatios->NevtPassErrStat();
			vector<double> nevFPEsyst = fFPRatios->NevtPassErrSyst();

			if(fVerbose > 2) cout << "   npp: " << nevFP[0] << "  nfp: " << nevFP[1] << "  nff: " << nevFP[2] << endl;
			npppred += nevFP[0];
			nfppred += nevFP[1];
			nffpred += nevFP[2];
			npppredEstat2 += nevFPEstat[0]*nevFPEstat[0];
			npppredEsyst2 += nevFPEsyst[0]*nevFPEsyst[0];
			nfppredEstat2 += nevFPEstat[1]*nevFPEstat[1];
			nfppredEsyst2 += nevFPEsyst[1]*nevFPEsyst[1];
			nffpredEstat2 += nevFPEstat[2]*nevFPEstat[2];
			nffpredEsyst2 += nevFPEsyst[2]*nevFPEsyst[2];
			delete fFPRatios;
		}
		Int_t bin = H_nsigpred->FindBin(pt1);
		H_nsigpred->SetBinContent(bin, npppred);
		H_nfppred ->SetBinContent(bin, nfppred);
		H_nffpred ->SetBinContent(bin, nffpred);
		H_nsigpred->SetBinError(bin, sqrt(npppredEstat2 + npppredEsyst2));
		H_nfppred ->SetBinError(bin, sqrt(nfppredEstat2 + nfppredEsyst2));
		H_nffpred ->SetBinError(bin, sqrt(nffpredEstat2 + nffpredEsyst2));
	}	

	if(fVerbose > 2) cout << " Predict " << H_nsigpred->Integral() << " signal events (Nsig = p^2*Npp) from this sample" << endl;
	if(fVerbose > 2) cout << " Predict " << H_nfppred->Integral() <<  " fake-prompt events (f*p*Nfp) from this sample" << endl;
	if(fVerbose > 2) cout << " Predict " << H_nffpred->Integral() <<  " fake-fake events (f*f*Nff) from this sample" << endl;
	// Output
	H_nsigpred->SetXTitle(convertVarName("MuPt[0]"));
	if(output){
		H_nt2mes->SetXTitle(convertVarName("MuPt[0]"));
		H_nt2mes->SetYTitle(convertVarName("MuPt[0]"));
		H_nt1mes->SetXTitle(convertVarName("MuPt[0]"));
		H_nt1mes->SetYTitle(convertVarName("MuPt[0]"));
		H_nt0mes->SetXTitle(convertVarName("MuPt[0]"));
		H_nt0mes->SetYTitle(convertVarName("MuPt[0]"));
		// if(H_nsigpred->GetMinimum() < 0) H_nsigpred->SetMaximum(0);
		// if(H_nsigpred->GetMinimum() > 0) H_nsigpred->SetMinimum(0);
		printObject(H_nsigpred, H_nsigpred->GetName(), H_nsigpred->GetTitle());
		printObject(H_nt2mes, H_nt2mes->GetName(), H_nt2mes->GetTitle(), "colz");
		printObject(H_nt1mes, H_nt1mes->GetName(), H_nt1mes->GetTitle(), "colz");
		printObject(H_nt0mes, H_nt0mes->GetName(), H_nt0mes->GetTitle(), "colz");
	}
	res.push_back(H_nsigpred);
	res.push_back(H_nfppred);
	res.push_back(H_nffpred);
	return res;
}
vector<TH1D*> MuonPlotter::ElElFPPrediction(TH2D* H_fratio, TH2D* H_pratio, TH2D* H_nt2mes, TH2D* H_nt1mes, TH2D* H_nt0mes, bool output){
	///////////////////////////////////////////////////////////////////////////
	// Note: Careful, this is only a workaround at the moment, and is only
	//       really valid for one single eta bin!
	//       In the future we should rewrite the interface to FPRatios and
	//       pass it the ratios directly to have full control
	///////////////////////////////////////////////////////////////////////////
	vector<TH1D*> res;
	const int nptbins = getNPt2Bins(Electron);

	TH1D *H_dummy    = new TH1D("ElDummyhist", "Dummy",                      H_nt2mes->GetNbinsX(), H_nt2mes->GetXaxis()->GetXbins()->GetArray());
	TH1D *H_nsigpred = new TH1D("ElNsigpred", "Predicted N_sig in Pt1 bins", H_nt2mes->GetNbinsX(), H_nt2mes->GetXaxis()->GetXbins()->GetArray());
	TH1D *H_nfppred  = new TH1D("ElNfppred",  "Predicted N_fp in Pt1 bins",  H_nt2mes->GetNbinsX(), H_nt2mes->GetXaxis()->GetXbins()->GetArray());
	TH1D *H_nffpred  = new TH1D("ElNffpred",  "Predicted N_ff in Pt1 bins",  H_nt2mes->GetNbinsX(), H_nt2mes->GetXaxis()->GetXbins()->GetArray());	

	if(fVerbose > 2){
		cout << " Nt2 = " << H_nt2mes->GetEntries();
		cout << " Nt1 = " << H_nt1mes->GetEntries();
		cout << " Nt0 = " << H_nt0mes->GetEntries() << endl;
	}

	for(size_t i = 1; i <= nptbins; ++i){
		double pt1 = H_dummy->GetBinCenter(i);
		if(fVerbose > 2){
			cout << "=======================================" << endl;
			cout << "Pt1 = " << pt1 << endl;
		}
		double eta1 = 0.0;
		double npppred(0.0), npppredEstat2(0.0), npppredEsyst2(0.0);
		double nfppred(0.0), nfppredEstat2(0.0), nfppredEsyst2(0.0);
		double nffpred(0.0), nffpredEstat2(0.0), nffpredEsyst2(0.0);
		for(size_t j = 1; j <= nptbins; ++j){
			if(fVerbose > 2) cout << " --------" << endl;
			double pt2 = H_dummy->GetBinCenter(j);

			double eta2 = 0.0;
			double nt2 = H_nt2mes->GetBinContent(i,j);
			double nt1 = H_nt1mes->GetBinContent(i,j);
			double nt0 = H_nt0mes->GetBinContent(i,j);

			if(fVerbose > 2) cout << "   Pt2 = " << pt2 << endl;
			if(fVerbose > 2) cout << "   nt2: " << nt2 << "  nt1: " << nt1 << "  nt0: " << nt0 << endl;

			FPRatios *fFPRatios = new FPRatios();
			fFPRatios->SetVerbose(fVerbose);
			fFPRatios->SetElFratios(H_fratio);
			fFPRatios->SetElPratios(H_pratio);
			vector<double> nt;
			nt.push_back(nt0); nt.push_back(nt1); nt.push_back(nt2);
			fFPRatios->NevtTopol(2, 0, nt);

			vector<double> vpt, veta;
			vpt.push_back(pt1); vpt.push_back(pt2);
			veta.push_back(eta1); veta.push_back(eta2);

			vector<double> nevFP = fFPRatios->NevtPass(vpt, veta);
			vector<double> nevFPEstat = fFPRatios->NevtPassErrStat();
			vector<double> nevFPEsyst = fFPRatios->NevtPassErrSyst();
			npppred += nevFP[0];
			nfppred += nevFP[1];
			nffpred += nevFP[2];
			npppredEstat2 += nevFPEstat[0]*nevFPEstat[0];
			npppredEsyst2 += nevFPEsyst[0]*nevFPEsyst[0];
			nfppredEstat2 += nevFPEstat[1]*nevFPEstat[1];
			nfppredEsyst2 += nevFPEsyst[1]*nevFPEsyst[1];
			nffpredEstat2 += nevFPEstat[2]*nevFPEstat[2];
			nffpredEsyst2 += nevFPEsyst[2]*nevFPEsyst[2];
			delete fFPRatios;
		}
		Int_t bin = H_nsigpred->FindBin(pt1);
		H_nsigpred->SetBinContent(bin, npppred);
		H_nfppred ->SetBinContent(bin, nfppred);
		H_nffpred ->SetBinContent(bin, nffpred);
		H_nsigpred->SetBinError(bin, sqrt(npppredEstat2 + npppredEsyst2));
		H_nfppred ->SetBinError(bin, sqrt(nfppredEstat2 + nfppredEsyst2));
		H_nffpred ->SetBinError(bin, sqrt(nffpredEstat2 + nffpredEsyst2));
	}	

	if(fVerbose > 2) cout << " Predict " << H_nsigpred->Integral() << " signal events (Nsig = p^2*Npp) from this sample" << endl;
	if(fVerbose > 2) cout << " Predict " << H_nfppred->Integral() <<  " fake-prompt events (f*p*Nfp) from this sample" << endl;
	if(fVerbose > 2) cout << " Predict " << H_nffpred->Integral() <<  " fake-fake events (f*f*Nff) from this sample" << endl;
	// Output
	H_nsigpred->SetXTitle(convertVarName("ElPt[0]"));
	if(output){
		H_nt2mes->SetXTitle(convertVarName("ElPt[0]"));
		H_nt2mes->SetYTitle(convertVarName("ElPt[0]"));
		H_nt1mes->SetXTitle(convertVarName("ElPt[0]"));
		H_nt1mes->SetYTitle(convertVarName("ElPt[0]"));
		H_nt0mes->SetXTitle(convertVarName("ElPt[0]"));
		H_nt0mes->SetYTitle(convertVarName("ElPt[0]"));
		// if(H_nsigpred->GetMinimum() < 0) H_nsigpred->SetMaximum(0);
		// if(H_nsigpred->GetMinimum() > 0) H_nsigpred->SetMinimum(0);
		printObject(H_nsigpred, H_nsigpred->GetName(), H_nsigpred->GetTitle());
		printObject(H_nt2mes, H_nt2mes->GetName(), H_nt2mes->GetTitle(), "colz");
		printObject(H_nt1mes, H_nt1mes->GetName(), H_nt1mes->GetTitle(), "colz");
		printObject(H_nt0mes, H_nt0mes->GetName(), H_nt0mes->GetTitle(), "colz");
	}
	res.push_back(H_nsigpred);
	res.push_back(H_nfppred);
	res.push_back(H_nffpred);
	return res;
}
vector<TH1D*> MuonPlotter::ElMuFPPrediction(TH2D* H_mufratio, TH2D* H_mupratio, TH2D* H_elfratio, TH2D* H_elpratio, TH2D* H_nt2mes, TH2D* H_nt10mes, TH2D* H_nt01mes, TH2D* H_nt0mes, bool output){
	///////////////////////////////////////////////////////////////////////////
	// Note: Careful, this is only a workaround at the moment, and is only
	//       really valid for one single eta bin!
	//       In the future we should rewrite the interface to FPRatios and
	//       pass it the ratios directly to have full control
	///////////////////////////////////////////////////////////////////////////
	vector<TH1D*> res;

	TH1D *H_mudummy = new TH1D("ElMuDummy",   "Dummy1",                      H_nt2mes->GetNbinsX(), H_nt2mes->GetXaxis()->GetXbins()->GetArray()); // x axis is muon
	TH1D *H_eldummy = new TH1D("MuElDummy",   "Dummy2",                      H_nt2mes->GetNbinsY(), H_nt2mes->GetYaxis()->GetXbins()->GetArray()); // y axis is electron
	TH1D *H_npppred = new TH1D("ElMuNpppred", "Predicted N_pp in MuPt bins", H_nt2mes->GetNbinsX(), H_nt2mes->GetXaxis()->GetXbins()->GetArray());
	TH1D *H_nfppred = new TH1D("ElMuNfppred", "Predicted N_fp in MuPt bins", H_nt2mes->GetNbinsX(), H_nt2mes->GetXaxis()->GetXbins()->GetArray());
	TH1D *H_npfpred = new TH1D("ElMuNpfpred", "Predicted N_pf in MuPt bins", H_nt2mes->GetNbinsX(), H_nt2mes->GetXaxis()->GetXbins()->GetArray());
	TH1D *H_nffpred = new TH1D("ElMuNffpred", "Predicted N_ff in MuPt bins", H_nt2mes->GetNbinsX(), H_nt2mes->GetXaxis()->GetXbins()->GetArray());	

	if(fVerbose > 2){
		cout << " Nt2  = " << H_nt2mes->GetEntries();
		cout << " Nt10 = " << H_nt10mes->GetEntries();
		cout << " Nt01 = " << H_nt01mes->GetEntries();
		cout << " Nt0  = " << H_nt0mes->GetEntries() << endl;
	}

	const double eta = 0.0;
	for(size_t i = 1; i <= H_mudummy->GetNbinsX(); ++i){ // loop on mu bins
		double mupt = H_mudummy->GetBinCenter(i);
		if(fVerbose > 2){
			cout << "=======================================" << endl;
			cout << "MuPt = " << mupt << "  MufRatio = " << H_mufratio->GetBinContent(i, 1) << "  MupRatio = " << H_mupratio->GetBinContent(i, 1) << endl;
		}
		double npppred(0.0), npppredEstat2(0.0), npppredEsyst2(0.0);
		double nfppred(0.0), nfppredEstat2(0.0), nfppredEsyst2(0.0);
		double npfpred(0.0), npfpredEstat2(0.0), npfpredEsyst2(0.0);
		double nffpred(0.0), nffpredEstat2(0.0), nffpredEsyst2(0.0);
		for(size_t j = 1; j <= H_eldummy->GetNbinsX(); ++j){ // loop on el bins
			if(fVerbose > 2) cout << " --------" << endl;
			double elpt = H_eldummy->GetBinCenter(j);

			double nt2  = H_nt2mes->GetBinContent(i,j);
			double nt10 = H_nt10mes->GetBinContent(i,j);
			double nt01 = H_nt01mes->GetBinContent(i,j);
			double nt0  = H_nt0mes->GetBinContent(i,j);

			if(fVerbose > 2) cout << "ElPt = " << elpt << "  ElfRatio = " << H_elfratio->GetBinContent(j, 1) << "  ElpRatio = " << H_elpratio->GetBinContent(j, 1) << endl;
			if(fVerbose > 2) cout << "   nt2: " << nt2 << "  nt10: " << nt10  << "  nt01: " << nt01 << "  nt0: " << nt0 << endl;

			FPRatios *fFPRatios = new FPRatios();
			fFPRatios->SetVerbose(fVerbose);
			fFPRatios->SetMuFratios(H_mufratio);
			fFPRatios->SetMuPratios(H_mupratio);
			fFPRatios->SetElFratios(H_elfratio);
			fFPRatios->SetElPratios(H_elpratio);
			vector<double> nt;
			nt.push_back(nt0);  // none pass
			nt.push_back(nt01); // e passes   >> Here I turn them around because FPRatios treats the first index as the electron,
			nt.push_back(nt10); // mu passes  >> and I have the opposite convention since muons rule electrons as everyone knows
			nt.push_back(nt2);  // both pass
			fFPRatios->NevtTopol(1, 1, nt);

			vector<double> vpt, veta;
			vpt.push_back(elpt); vpt.push_back(mupt);
			veta.push_back(eta); veta.push_back(eta);

			vector<double> nevFP      = fFPRatios->NevtPass(vpt, veta);
			vector<double> nevFPEstat = fFPRatios->NevtPassErrStat();
			vector<double> nevFPEsyst = fFPRatios->NevtPassErrSyst();
			npppred += nevFP[0];
			nfppred += nevFP[1];
			npfpred += nevFP[2];
			nffpred += nevFP[3];
			if(fVerbose > 2) cout << "   npp: " << nevFP[0] << "  nfp: " << nevFP[1]  << "  npf: " << nevFP[2] << "  nff: " << nevFP[3] << endl;
			npppredEstat2 += nevFPEstat[0]*nevFPEstat[0];
			npppredEsyst2 += nevFPEsyst[0]*nevFPEsyst[0];
			nfppredEstat2 += nevFPEstat[1]*nevFPEstat[1];
			nfppredEsyst2 += nevFPEsyst[1]*nevFPEsyst[1];
			npfpredEstat2 += nevFPEstat[2]*nevFPEstat[2];
			npfpredEsyst2 += nevFPEsyst[2]*nevFPEsyst[2];
			nffpredEstat2 += nevFPEstat[3]*nevFPEstat[3];
			nffpredEsyst2 += nevFPEsyst[3]*nevFPEsyst[3];
			delete fFPRatios;
		}
		Int_t bin = H_npppred->FindBin(mupt);
		H_npppred->SetBinContent(bin, npppred);
		H_nfppred->SetBinContent(bin, nfppred);
		H_npfpred->SetBinContent(bin, npfpred);
		H_nffpred->SetBinContent(bin, nffpred);
		H_npppred->SetBinError(bin, sqrt(npppredEstat2 + npppredEsyst2));
		H_nfppred->SetBinError(bin, sqrt(nfppredEstat2 + nfppredEsyst2));
		H_npfpred->SetBinError(bin, sqrt(npfpredEstat2 + npfpredEsyst2));
		H_nffpred->SetBinError(bin, sqrt(nffpredEstat2 + nffpredEsyst2));
	}	

	if(fVerbose > 2) cout << " Total: ";
	if(fVerbose > 2) cout << "   npp: " << H_npppred->Integral() << "  nfp: " << H_nfppred->Integral()  << "  npf: " << H_npfpred->Integral() << "  nff: " << H_nffpred->Integral() << endl;

	// Output
	H_npppred->SetXTitle(convertVarName("MuPt[0]"));
	if(output){
		H_nt2mes->SetXTitle(convertVarName("MuPt[0]"));
		H_nt2mes->SetYTitle(convertVarName("MuPt[0]"));
		H_nt10mes->SetXTitle(convertVarName("MuPt[0]"));
		H_nt10mes->SetYTitle(convertVarName("MuPt[0]"));
		H_nt01mes->SetXTitle(convertVarName("MuPt[0]"));
		H_nt01mes->SetYTitle(convertVarName("MuPt[0]"));
		H_nt0mes->SetXTitle(convertVarName("MuPt[0]"));
		H_nt0mes->SetYTitle(convertVarName("MuPt[0]"));
		printObject(H_npppred, H_npppred->GetName(), H_npppred->GetTitle());
		printObject(H_nt2mes,  H_nt2mes->GetName(),  H_nt2mes->GetTitle(),  "colz");
		printObject(H_nt10mes, H_nt10mes->GetName(), H_nt10mes->GetTitle(), "colz");
		printObject(H_nt01mes, H_nt01mes->GetName(), H_nt01mes->GetTitle(), "colz");
		printObject(H_nt0mes,  H_nt0mes->GetName(),  H_nt0mes->GetTitle(),  "colz");
	}
	res.push_back(H_npppred);
	res.push_back(H_nfppred);
	res.push_back(H_npfpred);
	res.push_back(H_nffpred);
	return res;
}

//____________________________________________________________________________
void MuonPlotter::fillYields(Sample *S){
	// MuMu Channel
	fCurrentChannel = Muon;
	int mu1(-1), mu2(-1);
	bool isdata = S->isdata;
	if( (S->isdata && isSSLLMuEventTRG(mu1, mu2)) || (!S->isdata && isSSLLMuEvent(mu1, mu2)) ){
		if(  isTightMuon(mu1) &&  isTightMuon(mu2) ){ // Tight-tight
			fCounters[fCurrentSample][Muon].fill(" ... first muon passes tight cut");
			fCounters[fCurrentSample][Muon].fill(" ... second muon passes tight cut");
			fCounters[fCurrentSample][Muon].fill(" ... both muons pass tight cut");
			S->region[Signal].mm.nt20_pt ->Fill(MuPt [mu1], MuPt [mu2]);
			S->region[Signal].mm.nt20_eta->Fill(MuEta[mu1], MuEta[mu2]);
			if(S->isdata) fOUTSTREAM << " Mu/Mu Tight-Tight event in run " << setw(7) << Run << " event " << setw(13) << Event << " lumisection " << setw(5) << LumiSec << " in dataset " << setw(9) << S->sname << endl;
		}
		if(  isTightMuon(mu1) && !isTightMuon(mu2) ){ // Tight-loose
			fCounters[fCurrentSample][Muon].fill(" ... first muon passes tight cut");
			S->region[Signal].mm.nt10_pt ->Fill(MuPt [mu1], MuPt [mu2]);
			S->region[Signal].mm.nt10_eta->Fill(MuEta[mu1], MuEta[mu2]);
		}
		if( !isTightMuon(mu1) &&  isTightMuon(mu2) ){ // Loose-tight
			fCounters[fCurrentSample][Muon].fill(" ... second muon passes tight cut");
			S->region[Signal].mm.nt10_pt ->Fill(MuPt [mu2], MuPt [mu1]); // tight one always in x axis; fill same again
			S->region[Signal].mm.nt10_eta->Fill(MuEta[mu2], MuEta[mu1]);
		}
		if( !isTightMuon(mu1) && !isTightMuon(mu2) ){ // Loose-loose
			S->region[Signal].mm.nt00_pt ->Fill(MuPt [mu1], MuPt [mu2]);
			S->region[Signal].mm.nt00_eta->Fill(MuEta[mu1], MuEta[mu2]);
		}
	}
	if( (S->isdata && isSigSupMuEventTRG()) || (!S->isdata && isSigSupMuEvent()) ){ // f Ratio
		if( isTightMuon(0) ) S->region[Signal].mm.fntight->Fill(MuPt[0], MuEta[0]);
		if( isLooseMuon(0) ) S->region[Signal].mm.fnloose->Fill(MuPt[0], MuEta[0]);
	}
	if( (S->isdata && isZMuMuEventTRG())    || (!S->isdata && isZMuMuEvent()) ){ // p Ratio
		if( isTightMuon(0) ) S->region[Signal].mm.pntight->Fill(MuPt[0], MuEta[0]);
		if( isLooseMuon(0) ) S->region[Signal].mm.pnloose->Fill(MuPt[0], MuEta[0]);
	}				
	
	// EE Channel
	fCurrentChannel = Electron;
	int el1(-1), el2(-1);
	if( (S->isdata && isSSLLElEventTRG(el1, el2)) || (!S->isdata && isSSLLElEvent(el1, el2)) ){
		if(  isTightElectron(el1) &&  isTightElectron(el2) ){ // Tight-tight
			fCounters[fCurrentSample][Electron].fill(" ... first electron passes tight cut");
			fCounters[fCurrentSample][Electron].fill(" ... second electron passes tight cut");
			fCounters[fCurrentSample][Electron].fill(" ... both electrons pass tight cut");
			S->region[Signal].ee.nt20_pt ->Fill(ElPt [el1], ElPt [el2]);
			S->region[Signal].ee.nt20_eta->Fill(ElEta[el1], ElEta[el2]);
			if(S->isdata) fOUTSTREAM << " E/E Tight-Tight event in run   " << setw(7) << Run << " event " << setw(13) << Event << " lumisection " << setw(5) << LumiSec << " in dataset " << setw(9) << S->sname << endl;
		}
		if(  isTightElectron(el1) && !isTightElectron(el2) ){ // Tight-loose
			fCounters[fCurrentSample][Electron].fill(" ... first electron passes tight cut");
			S->region[Signal].ee.nt10_pt ->Fill(ElPt [el1], ElPt [el2]);
			S->region[Signal].ee.nt10_eta->Fill(ElEta[el1], ElEta[el2]);
		}
		if( !isTightElectron(el1) &&  isTightElectron(el2) ){ // Loose-tight
			fCounters[fCurrentSample][Electron].fill(" ... second electron passes tight cut");
			S->region[Signal].ee.nt10_pt ->Fill(ElPt [el2], ElPt [el1]); // tight one always in x axis; fill same again
			S->region[Signal].ee.nt10_eta->Fill(ElEta[el2], ElEta[el1]);
		}
		if( !isTightElectron(el1) && !isTightElectron(el2) ){ // Loose-loose
			S->region[Signal].ee.nt00_pt ->Fill(ElPt [el1], ElPt [el2]);
			S->region[Signal].ee.nt00_eta->Fill(ElEta[el1], ElEta[el2]);
		}
	}
	if( (S->isdata && isSigSupElEventTRG())   || (!S->isdata && isSigSupElEvent()) ){ // f Ratio
		if( isTightElectron(0) ) S->region[Signal].ee.fntight->Fill(ElPt[0], ElEta[0]);
		if( isLooseElectron(0) ) S->region[Signal].ee.fnloose->Fill(ElPt[0], ElEta[0]);
	}
	int elind;
	if( (S->isdata && isZElElEventTRG(elind)) || (!S->isdata && isZElElEvent(elind)) ){ // p Ratio
		if( isTightElectron(elind) ) S->region[Signal].ee.pntight->Fill(ElPt[elind], ElEta[elind]);
		if( isLooseElectron(elind) ) S->region[Signal].ee.pnloose->Fill(ElPt[elind], ElEta[elind]);
	}
				
	// EMu Channel
	fCurrentChannel = EMu;
	int mu(-1), el(-1);
	if( (S->isdata && isSSLLElMuEventTRG(mu, el)) || (!S->isdata && isSSLLElMuEvent(mu, el)) ){
		if(  isTightElectron(el) &&  isTightMuon(mu) ){ // Tight-tight
			fCounters[fCurrentSample][EMu].fill(" ... muon passes tight cut");
			fCounters[fCurrentSample][EMu].fill(" ... electron passes tight cut");
			fCounters[fCurrentSample][EMu].fill(" ... both e and mu pass tight cuts");
			S->region[Signal].em.nt20_pt ->Fill(MuPt [mu], ElPt [el]);
			S->region[Signal].em.nt20_eta->Fill(MuEta[mu], ElEta[el]);
			if(S->isdata) fOUTSTREAM << " E/Mu Tight-Tight event in run  " << setw(7) << Run << " event " << setw(13) << Event << " lumisection " << setw(5) << LumiSec << " in dataset " << setw(9) << S->sname << endl;
		}
		if( !isTightElectron(el) &&  isTightMuon(mu) ){ // Tight-loose
			fCounters[fCurrentSample][EMu].fill(" ... muon passes tight cut");
			S->region[Signal].em.nt10_pt ->Fill(MuPt [mu], ElPt [el]);
			S->region[Signal].em.nt10_eta->Fill(MuEta[mu], ElEta[el]);
		}
		if(  isTightElectron(el) && !isTightMuon(mu) ){ // Loose-tight
			fCounters[fCurrentSample][EMu].fill(" ... electron passes tight cut");
			S->region[Signal].em.nt01_pt ->Fill(MuPt [mu], ElPt [el]); // muon always in x axis for e/mu
			S->region[Signal].em.nt01_eta->Fill(MuEta[mu], ElEta[el]);
		}
		if( !isTightElectron(0) && !isTightMuon(0) ){ // Loose-loose
			S->region[Signal].em.nt00_pt ->Fill(MuPt [mu], ElPt [el]);
			S->region[Signal].em.nt00_eta->Fill(MuEta[mu], ElEta[el]);
		}
	}
}

//____________________________________________________________________________
void MuonPlotter::storeNumbers(Sample *S, gChannel chan){
	Channel *C;
	if(chan == Muon)     C = &S->region[Signal].mm;
	if(chan == Electron) C = &S->region[Signal].ee;
	if(chan == EMu)      C = &S->region[Signal].em;
	S->numbers[chan].nt2  = C->nt20_pt->GetEntries();
	S->numbers[chan].nt10 = C->nt10_pt->GetEntries();
	S->numbers[chan].nt01 = C->nt01_pt->GetEntries();
	S->numbers[chan].nt0  = C->nt00_pt->GetEntries();
	if(chan != EMu){
		S->numbers[chan].nsst = C->fntight->GetEntries();
		S->numbers[chan].nssl = C->fnloose->GetEntries();
		S->numbers[chan].nzt  = C->pntight->GetEntries();
		S->numbers[chan].nzl  = C->pnloose->GetEntries();
	}
}

//____________________________________________________________________________
void MuonPlotter::initCounters(int sample){
	if(sample < 0) sample = fCurrentSample;
	fCounters[sample][Muon]    .fill("All events",                             0.);
	fCounters[sample][Muon]    .fill(" ... is good run",                       0.);
	fCounters[sample][Muon]    .fill(" ... passes triggers",                   0.);
	fCounters[sample][Muon]    .fill(" ... has at least 2 jets, 1 loose muon", 0.);
	fCounters[sample][Muon]    .fill(" ... has 2 loose muons",                 0.);
	fCounters[sample][Muon]    .fill(" ... passes Z veto",                     0.);
	fCounters[sample][Muon]    .fill(" ... passes Minv veto",                  0.);
	fCounters[sample][Muon]    .fill(" ... passes HT cut",                     0.);
	fCounters[sample][Muon]    .fill(" ... passes MET cut",                    0.);
	if(fChargeSwitch == 0) fCounters[sample][Muon]    .fill(" ... has same-sign muons",     0.);
	if(fChargeSwitch == 1) fCounters[sample][Muon]    .fill(" ... has opposite-sign muons", 0.);
	fCounters[sample][Muon]    .fill(" ... second muon passes pt cut",         0.);
	fCounters[sample][Muon]    .fill(" ... first muon passes pt cut",          0.);
	fCounters[sample][Muon]    .fill(" ... first muon passes tight cut",       0.);
	fCounters[sample][Muon]    .fill(" ... second muon passes tight cut",      0.);
	fCounters[sample][Muon]    .fill(" ... both muons pass tight cut",         0.);

	fCounters[sample][EMu]     .fill("All events",                             0.);
	fCounters[sample][EMu]     .fill(" ... is good run",                       0.);
	fCounters[sample][EMu]     .fill(" ... passes triggers",                   0.);
	fCounters[sample][EMu]     .fill(" ... has at least 1 j, loose e/mu pair", 0.);
	fCounters[sample][EMu]     .fill(" ... passes Z veto",                     0.);
	fCounters[sample][EMu]     .fill(" ... passes Minv veto",                  0.);
	fCounters[sample][EMu]     .fill(" ... passes HT cut",                     0.);
	fCounters[sample][EMu]     .fill(" ... passes MET cut",                    0.);
	if(fChargeSwitch == 0) fCounters[sample][EMu]     .fill(" ... has same-sign electron muon pair", 0.);
	if(fChargeSwitch == 1) fCounters[sample][EMu]     .fill(" ... has opp.-sign electron muon pair", 0.);
	fCounters[sample][EMu]     .fill(" ... muon passes pt cut",                0.);
	fCounters[sample][EMu]     .fill(" ... electron passes pt cut",            0.);
	fCounters[sample][EMu]     .fill(" ... muon passes tight cut",             0.);
	fCounters[sample][EMu]     .fill(" ... electron passes tight cut",         0.);
	fCounters[sample][EMu]     .fill(" ... both e and mu pass tight cuts",     0.);

	fCounters[sample][Electron].fill("All events",                             0.);
	fCounters[sample][Electron].fill(" ... is good run",                       0.);
	fCounters[sample][Electron].fill(" ... passes triggers",                   0.);
	fCounters[sample][Electron].fill(" ... has at least 2 jets, 1 loose el",   0.);
	fCounters[sample][Electron].fill(" ... has 2 loose electrons",             0.);
	fCounters[sample][Electron].fill(" ... passes Z veto",                     0.);
	fCounters[sample][Electron].fill(" ... passes Minv veto",                  0.);
	fCounters[sample][Electron].fill(" ... passes HT cut",                     0.);
	fCounters[sample][Electron].fill(" ... passes MET cut",                    0.);
	if(fChargeSwitch == 0) fCounters[sample][Electron].fill(" ... has same-sign electrons",     0.);
	if(fChargeSwitch == 1) fCounters[sample][Electron].fill(" ... has opposite-sign electrons", 0.);
	fCounters[sample][Electron].fill(" ... second electron passes pt cut",     0.);
	fCounters[sample][Electron].fill(" ... first electron passes pt cut",      0.);
	fCounters[sample][Electron].fill(" ... first electron passes tight cut",   0.);
	fCounters[sample][Electron].fill(" ... second electron passes tight cut",  0.);
	fCounters[sample][Electron].fill(" ... both electrons pass tight cut",     0.);
}
void MuonPlotter::printCutFlows(TString filename){
	// Remove existing cutflow file
	// char cmd[100];
	// sprintf(cmd,"rm -f %s", filename.Data());
	// system(cmd);

	ofstream OUT(filename.Data(), ios::trunc);
	vector<string>::const_iterator countit;
	
	OUT << " Printing Cutflow for Mu/Mu channel..." << endl;
	OUT << "------------------------------------------------------------------------------------------------------------------------------" << endl;
	OUT << " Cutname                                 | ";
	for(gSample i = sample_begin; i < gNSAMPLES; i=gSample(i+1)){
		if(fSamples[i].isdata == 0) continue;
		OUT << setw(9) << fSamples[i].sname;
		OUT << " | ";
	}
	OUT << endl;
	OUT << "------------------------------------------------------------------------------------------------------------------------------" << endl;

    for( countit = fCounters[0][Muon].countNames.begin(); countit != fCounters[0][Muon].countNames.end(); ++countit ){
		OUT << setw(40) << *countit;
		OUT << " | ";
		
		for(gSample i = sample_begin; i < gNSAMPLES; i=gSample(i+1)){
			if(fSamples[i].isdata == 0) continue;
			OUT << setw(9) << setprecision(3) << fCounters[i][Muon].counts(*countit) << " | ";
		}
		OUT << endl;
	}
	OUT << "------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" << endl;
	OUT << " Cutname                                 | ";
	for(gSample i = sample_begin; i < gNSAMPLES; i=gSample(i+1)){
		if(fSamples[i].isdata == 1) continue;
		OUT << setw(9) << fSamples[i].sname;
		OUT << " | ";
	}
	OUT << endl;
	OUT << "------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" << endl;

    for( countit = fCounters[0][Muon].countNames.begin(); countit != fCounters[0][Muon].countNames.end(); ++countit ){
		OUT << setw(40) << *countit;
		OUT << " | ";
		
		for(gSample i = sample_begin; i < gNSAMPLES; i=gSample(i+1)){
			if(fSamples[i].isdata == 1) continue;
			OUT << setw(9) << setprecision(3) << fCounters[i][Muon].counts(*countit) << " | ";
		}
		OUT << endl;
	}
	OUT << "------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" << endl;
	OUT << endl << endl;
	OUT << " Printing Cutflow for E/Mu channel..." << endl;
	OUT << "------------------------------------------------------------------------------------------------------------------------------" << endl;
	OUT << " Cutname                                 | ";
	for(gSample i = sample_begin; i < gNSAMPLES; i=gSample(i+1)){
		if(fSamples[i].isdata == 0) continue;
		OUT << setw(9) << fSamples[i].sname;
		OUT << " | ";
	}
	OUT << endl;
	OUT << "------------------------------------------------------------------------------------------------------------------------------" << endl;

    for( countit = fCounters[0][EMu].countNames.begin(); countit != fCounters[0][EMu].countNames.end(); ++countit ){
		OUT << setw(40) << *countit;
		OUT << " | ";
		
		for(gSample i = sample_begin; i < gNSAMPLES; i=gSample(i+1)){
			if(fSamples[i].isdata == 0) continue;
			OUT << setw(9) << setprecision(3) << fCounters[i][EMu].counts(*countit) << " | ";
		}
		OUT << endl;
	}
	OUT << "------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" << endl;
	OUT << " Cutname                                 | ";
	for(gSample i = sample_begin; i < gNSAMPLES; i=gSample(i+1)){
		if(fSamples[i].isdata == 1) continue;
		OUT << setw(9) << fSamples[i].sname;
		OUT << " | ";
	}
	OUT << endl;
	OUT << "------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" << endl;

    for( countit = fCounters[0][EMu].countNames.begin(); countit != fCounters[0][EMu].countNames.end(); ++countit ){
		OUT << setw(40) << *countit;
		OUT << " | ";
		
		for(gSample i = sample_begin; i < gNSAMPLES; i=gSample(i+1)){
			if(fSamples[i].isdata == 1) continue;
			OUT << setw(9) << setprecision(3) << fCounters[i][EMu].counts(*countit) << " | ";
		}
		OUT << endl;
	}
	OUT << "------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" << endl;
	OUT << endl << endl;
	OUT << " Printing Cutflow for E/E channel..." << endl;
	OUT << "------------------------------------------------------------------------------------------------------------------------------" << endl;
	OUT << " Cutname                                 | ";
	for(gSample i = sample_begin; i < gNSAMPLES; i=gSample(i+1)){
		if(fSamples[i].isdata == 0) continue;
		OUT << setw(9) << fSamples[i].sname;
		OUT << " | ";
	}
	OUT << endl;
	OUT << "------------------------------------------------------------------------------------------------------------------------------" << endl;

    for( countit = fCounters[0][Electron].countNames.begin(); countit != fCounters[0][Electron].countNames.end(); ++countit ){
		OUT << setw(40) << *countit;
		OUT << " | ";
		
		for(gSample i = sample_begin; i < gNSAMPLES; i=gSample(i+1)){
			if(fSamples[i].isdata == 0) continue;
			OUT << setw(9) << setprecision(3) << fCounters[i][Electron].counts(*countit) << " | ";
		}
		OUT << endl;
	}
	OUT << "------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" << endl;
	OUT << " Cutname                                 | ";
	for(gSample i = sample_begin; i < gNSAMPLES; i=gSample(i+1)){
		if(fSamples[i].isdata == 1) continue;
		OUT << setw(9) << fSamples[i].sname;
		OUT << " | ";
	}
	OUT << endl;
	OUT << "------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" << endl;

    for( countit = fCounters[0][Electron].countNames.begin(); countit != fCounters[0][Electron].countNames.end(); ++countit ){
		OUT << setw(40) << *countit;
		OUT << " | ";
		
		for(gSample i = sample_begin; i < gNSAMPLES; i=gSample(i+1)){
			if(fSamples[i].isdata == 1) continue;
			OUT << setw(9) << setprecision(3) << fCounters[i][Electron].counts(*countit) << " | ";
		}
		OUT << endl;
	}
	OUT << "------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" << endl;
	OUT.close();
}

//____________________________________________________________________________
void MuonPlotter::printYields(gChannel chan){ printYields(chan, fAllSamples); }
void MuonPlotter::printYields(gChannel chan, float scale){ printYields(chan, fAllSamples, scale); }
void MuonPlotter::printYields(gChannel chan, int sample, float lumiscale){ vector<int> samples; samples.push_back(sample); printYields(chan, samples, lumiscale); }
void MuonPlotter::printYields(gChannel chan, vector<int> samples, float lumiscale){
	cout << setfill('-') << setw(94) << "-" << endl;
	TString name = "Mu/Mu";
	if(chan == Electron) name = "E/E";
	if(chan == EMu) name = "E/Mu";
	cout << " Printing yields for " << name << " channel..." << endl;
	if(lumiscale > -1.0) cout << " Numbers scaled to " << lumiscale << " /pb" << endl;
	cout << "        Name |   Nt2   |   Nt10  |   Nt01  |   Nt0   |   Nsst  |   Nssl  |   NZ t  |   NZ l  |" << endl;
	cout << setfill('-') << setw(94) << "-" << endl;
	cout << setfill(' ');
	for(size_t i = 0; i < samples.size(); ++i){
		Sample S = fSamples[samples[i]];
		NumberSet numbers = S.numbers[chan];
		float scale = lumiscale / S.lumi;
		if(S.isdata || scale < 0) scale = 1;
		cout << setw(12) << S.sname << " |";
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

	cout << setfill('-') << setw(94) << "-" << endl;
}
void MuonPlotter::printYieldsShort(float luminorm){
	vector<int> musamples;
	vector<int> elsamples;
	vector<int> emusamples;

	if(fSelectionSwitch == 0){
		cout << "---------------------------------------------------" << endl;
		cout << "Printing yields for UCSD/UCSB/FNAL selection" << endl;
		musamples = fMuData;
		elsamples = fEGData;
		emusamples.push_back(MuA);
		emusamples.push_back(MuB);
		emusamples.push_back(EGA);
		emusamples.push_back(EGB);
	}
	if(fSelectionSwitch == 1){
		cout << "---------------------------------------------------" << endl;
		cout << "Printing yields for Florida selection" << endl;
		musamples  = fJMData;
		elsamples  = fJMData;
		emusamples = fJMData;
	}

	float nt20[gNCHANNELS],    nt10[gNCHANNELS],    nt01[gNCHANNELS],    nt00[gNCHANNELS];

	// float nt2_mumu(0.),    nt10_mumu(0.),    nt0_mumu(0.);
	// float nt2_emu(0.),    nt10_emu(0.),    nt01_emu(0.),    nt0_emu(0.);
	// float nt2_ee(0.),    nt10_ee(0.),    nt0_ee(0.);

	for(size_t i = 0; i < musamples.size(); ++i){
		Sample S = fSamples[musamples[i]];
		nt20[Muon]    += S.numbers[Muon].nt2;
		nt10[Muon]    += S.numbers[Muon].nt10;
		nt00[Muon]    += S.numbers[Muon].nt0;
	}		

	for(size_t i = 0; i < emusamples.size(); ++i){
		Sample S = fSamples[emusamples[i]];
		nt20[EMu]    += S.numbers[EMu].nt2;
		nt10[EMu]    += S.numbers[EMu].nt10;
		nt01[EMu]    += S.numbers[EMu].nt01;
		nt00[EMu]    += S.numbers[EMu].nt0;
	}		

	for(size_t i = 0; i < elsamples.size(); ++i){
		Sample S = fSamples[elsamples[i]];
		nt20[Electron]    += S.numbers[Electron].nt2;
		nt10[Electron]    += S.numbers[Electron].nt10;
		nt00[Electron]    += S.numbers[Electron].nt0;
	}		

	// for(size_t i = 0; i < musamples.size(); ++i){
	// 	int index = musamples[i];
	// 	Channel *mumu = &fSamples[index].mm;
	// 	nt2_mumu  += mumu->numbers.nt2;
	// 	nt10_mumu += mumu->numbers.nt10;
	// 	nt0_mumu  += mumu->numbers.nt0;
	// 	nt2_mumu_e2  += mumu->numbers.nt2;
	// 	nt10_mumu_e2 += mumu->numbers.nt10;
	// 	nt0_mumu_e2  += mumu->numbers.nt0;
	// }		
	// for(size_t i = 0; i < emusamples.size(); ++i){
	// 	int index = emusamples[i];
	// 	Channel *emu = &fSamples[index].em;
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
	// 	Channel *ee = &fSamples[index].ee;
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
		Sample S = fSamples[index];
		float scale = luminorm / S.lumi;
		if(luminorm < 0) scale = 1;
		nt2sum_mumu  += scale*S.numbers[Muon]    .nt2;
		nt10sum_mumu += scale*S.numbers[Muon]    .nt10;
		nt0sum_mumu  += scale*S.numbers[Muon]    .nt0;
		nt2sum_emu   += scale*S.numbers[EMu]     .nt2;
		nt10sum_emu  += scale*S.numbers[EMu]     .nt10;
		nt01sum_emu  += scale*S.numbers[EMu]     .nt01;
		nt0sum_emu   += scale*S.numbers[EMu]     .nt0;
		nt2sum_ee    += scale*S.numbers[Electron].nt2;
		nt10sum_ee   += scale*S.numbers[Electron].nt10;
		nt0sum_ee    += scale*S.numbers[Electron].nt0;

		cout << setw(9) << S.sname << " || ";
		cout << setw(7) << scale*S.numbers[Muon]    .nt2  << " | ";
		cout << setw(7) << scale*S.numbers[Muon]    .nt10 << " | ";
		cout << setw(7) << scale*S.numbers[Muon]    .nt0  << " || ";
		cout << setw(7) << scale*S.numbers[EMu]     .nt2  << " | ";
		cout << setw(7) << scale*S.numbers[EMu]     .nt10 << " | ";
		cout << setw(7) << scale*S.numbers[EMu]     .nt01 << " | ";
		cout << setw(7) << scale*S.numbers[EMu]     .nt0  << " || ";
		cout << setw(7) << scale*S.numbers[Electron].nt2  << " | ";
		cout << setw(7) << scale*S.numbers[Electron].nt10 << " | ";
		cout << setw(7) << scale*S.numbers[Electron].nt0  << " || ";
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
	cout << setw(9) << fSamples[LM0].sname << " || ";
	float scale = luminorm / fSamples[LM0].lumi;
	if(luminorm < 0) scale = 1;
	cout << setw(7) << scale * fSamples[LM0].numbers[Muon]    .nt2  << " | ";
	cout << setw(7) << scale * fSamples[LM0].numbers[Muon]    .nt10 << " | ";
	cout << setw(7) << scale * fSamples[LM0].numbers[Muon]    .nt0  << " || ";
	cout << setw(7) << scale * fSamples[LM0].numbers[EMu]     .nt2  << " | ";
	cout << setw(7) << scale * fSamples[LM0].numbers[EMu]     .nt10 << " | ";
	cout << setw(7) << scale * fSamples[LM0].numbers[EMu]     .nt01 << " | ";
	cout << setw(7) << scale * fSamples[LM0].numbers[EMu]     .nt0  << " || ";
	cout << setw(7) << scale * fSamples[LM0].numbers[Electron].nt2  << " | ";
	cout << setw(7) << scale * fSamples[LM0].numbers[Electron].nt10 << " | ";
	cout << setw(7) << scale * fSamples[LM0].numbers[Electron].nt0  << " || ";
	cout << endl;
	cout << "-------------------------------------------------------------------------------------------------------------------" << endl;
	cout << setw(9) << "data"  << " || ";
	cout << setw(7) << nt20[Muon]     << " | ";
	cout << setw(7) << nt10[Muon]     << " | ";
	cout << setw(7) << nt00[Muon]     << " || ";
	cout << setw(7) << nt20[EMu]      << " | ";
	cout << setw(7) << nt10[EMu]      << " | ";
	cout << setw(7) << nt01[EMu]      << " | ";
	cout << setw(7) << nt00[EMu]      << " || ";
	cout << setw(7) << nt20[Electron] << " | ";
	cout << setw(7) << nt10[Electron] << " | ";
	cout << setw(7) << nt00[Electron] << " || ";
	cout << endl;
	cout << "-------------------------------------------------------------------------------------------------------------------" << endl;
	
}

//____________________________________________________________________________
void MuonPlotter::bookHistos(){
	for(gSample i = sample_begin; i < gNSAMPLES; i=gSample(i+1)){
		Sample *S = &fSamples[i];
		for(gRegion r = region_begin; r < gNREGIONS; r=gRegion(r+1)){
			Region *R = &S->region[r];
			for(gChannel c = channels_begin; c < gNCHANNELS; c=gChannel(c+1)){
				Channel *C;
				if(c == Muon || c == Electron){
					if(c == Muon)     C = &R->mm;
					if(c == Electron) C = &R->ee;
					C->nt20_pt  = new TH2D(S->sname + "_" + R->sname + "_" + C->sname + "_NT20_pt",  "NT20_pt",  getNPt2Bins(c), getPt2Bins(c), getNPt2Bins(c), getPt2Bins(c)); C->nt20_pt ->Sumw2();
					C->nt10_pt  = new TH2D(S->sname + "_" + R->sname + "_" + C->sname + "_NT10_pt",  "NT10_pt",  getNPt2Bins(c), getPt2Bins(c), getNPt2Bins(c), getPt2Bins(c)); C->nt10_pt ->Sumw2();
					C->nt01_pt  = new TH2D(S->sname + "_" + R->sname + "_" + C->sname + "_NT01_pt",  "NT01_pt",  getNPt2Bins(c), getPt2Bins(c), getNPt2Bins(c), getPt2Bins(c)); C->nt01_pt ->Sumw2();
					C->nt00_pt  = new TH2D(S->sname + "_" + R->sname + "_" + C->sname + "_NT00_pt",  "NT00_pt",  getNPt2Bins(c), getPt2Bins(c), getNPt2Bins(c), getPt2Bins(c)); C->nt00_pt ->Sumw2();
					C->nt20_eta = new TH2D(S->sname + "_" + R->sname + "_" + C->sname + "_NT20_eta", "NT20_eta", getNEtaBins(c), getEtaBins(c), getNEtaBins(c), getEtaBins(c)); C->nt20_eta->Sumw2();
					C->nt10_eta = new TH2D(S->sname + "_" + R->sname + "_" + C->sname + "_NT10_eta", "NT10_eta", getNEtaBins(c), getEtaBins(c), getNEtaBins(c), getEtaBins(c)); C->nt10_eta->Sumw2();
					C->nt01_eta = new TH2D(S->sname + "_" + R->sname + "_" + C->sname + "_NT01_eta", "NT01_eta", getNEtaBins(c), getEtaBins(c), getNEtaBins(c), getEtaBins(c)); C->nt01_eta->Sumw2();
					C->nt00_eta = new TH2D(S->sname + "_" + R->sname + "_" + C->sname + "_NT00_eta", "NT00_eta", getNEtaBins(c), getEtaBins(c), getNEtaBins(c), getEtaBins(c)); C->nt00_eta->Sumw2();
					C->fntight  = new TH2D(S->sname + "_" + R->sname + "_" + C->sname + "_fNTight",  "fNTight",  getNPt2Bins(c), getPt2Bins(c), getNEtaBins(c), getEtaBins(c)); C->fntight  ->Sumw2();
					C->fnloose  = new TH2D(S->sname + "_" + R->sname + "_" + C->sname + "_fNLoose",  "fNLoose",  getNPt2Bins(c), getPt2Bins(c), getNEtaBins(c), getEtaBins(c)); C->fnloose  ->Sumw2();					
					C->pntight  = new TH2D(S->sname + "_" + R->sname + "_" + C->sname + "_pNTight",  "pNTight",  getNPt2Bins(c), getPt2Bins(c), getNEtaBins(c), getEtaBins(c)); C->pntight  ->Sumw2();
					C->pnloose  = new TH2D(S->sname + "_" + R->sname + "_" + C->sname + "_pNLoose",  "pNLoose",  getNPt2Bins(c), getPt2Bins(c), getNEtaBins(c), getEtaBins(c)); C->pnloose  ->Sumw2();					
				}
				else if(c == EMu){
					C = &R->em;
					C->nt20_pt  = new TH2D(S->sname + "_" + R->sname + "_" + C->sname + "_NT20_pt",  "NT20_pt",  getNPt2Bins(Muon), getPt2Bins(Muon), getNPt2Bins(Electron), getPt2Bins(Electron)); C->nt20_pt ->Sumw2();
					C->nt10_pt  = new TH2D(S->sname + "_" + R->sname + "_" + C->sname + "_NT10_pt",  "NT10_pt",  getNPt2Bins(Muon), getPt2Bins(Muon), getNPt2Bins(Electron), getPt2Bins(Electron)); C->nt10_pt ->Sumw2();
					C->nt01_pt  = new TH2D(S->sname + "_" + R->sname + "_" + C->sname + "_NT01_pt",  "NT01_pt",  getNPt2Bins(Muon), getPt2Bins(Muon), getNPt2Bins(Electron), getPt2Bins(Electron)); C->nt01_pt ->Sumw2();
					C->nt00_pt  = new TH2D(S->sname + "_" + R->sname + "_" + C->sname + "_NT00_pt",  "NT00_pt",  getNPt2Bins(Muon), getPt2Bins(Muon), getNPt2Bins(Electron), getPt2Bins(Electron)); C->nt00_pt ->Sumw2();
					C->nt20_eta = new TH2D(S->sname + "_" + R->sname + "_" + C->sname + "_NT20_eta", "NT20_eta", getNEtaBins(Muon), getEtaBins(Muon), getNEtaBins(Electron), getEtaBins(Electron)); C->nt20_eta->Sumw2();
					C->nt10_eta = new TH2D(S->sname + "_" + R->sname + "_" + C->sname + "_NT10_eta", "NT10_eta", getNEtaBins(Muon), getEtaBins(Muon), getNEtaBins(Electron), getEtaBins(Electron)); C->nt10_eta->Sumw2();
					C->nt01_eta = new TH2D(S->sname + "_" + R->sname + "_" + C->sname + "_NT01_eta", "NT01_eta", getNEtaBins(Muon), getEtaBins(Muon), getNEtaBins(Electron), getEtaBins(Electron)); C->nt01_eta->Sumw2();
					C->nt00_eta = new TH2D(S->sname + "_" + R->sname + "_" + C->sname + "_NT00_eta", "NT00_eta", getNEtaBins(Muon), getEtaBins(Muon), getNEtaBins(Electron), getEtaBins(Electron)); C->nt00_eta->Sumw2();
					// C->ntight   = new TH2D(S->sname + "_" + R->sname + "_" + C->sname + "_NTight",   "NTight",   getNPt2Bins(Muon), getPt2Bins(Muon), getNEtaBins(Muon), getEtaBins(Muon)); C->ntight  ->Sumw2(); // Don't really matter
					// C->nloose   = new TH2D(S->sname + "_" + R->sname + "_" + C->sname + "_NLoose",   "NLoose",   getNPt2Bins(Muon), getPt2Bins(Muon), getNEtaBins(Muon), getEtaBins(Muon)); C->nloose  ->Sumw2();					
				}
			}
		}
	}
}
void MuonPlotter::writeHistos(){
	TFile *pFile = new TFile(fOutputFileName, "RECREATE");
	pFile->cd();

	for(gSample i = sample_begin; i < gNSAMPLES; i=gSample(i+1)){
		Sample *S = &fSamples[i];
		TDirectory* sdir = Util::FindOrCreate(S->sname, pFile);
		sdir->cd();

		for(gRegion r = region_begin; r < gNREGIONS; r=gRegion(r+1)){
			Region *R = &S->region[r];
			TString temp = S->sname + "/" + R->sname;
			TDirectory* rdir = Util::FindOrCreate(temp, pFile);
			rdir->cd();
			
			for(gChannel ch = channels_begin; ch < gNCHANNELS; ch=gChannel(ch+1)){ // Loop over channels, mumu, emu, ee
				Channel *C;
				if(ch == Muon)     C = &R->mm;
				if(ch == Electron) C = &R->ee;
				if(ch == EMu)      C = &R->em;
				C->nt20_pt ->Write(C->nt20_pt ->GetName(), TObject::kWriteDelete);
				C->nt10_pt ->Write(C->nt10_pt ->GetName(), TObject::kWriteDelete);
				C->nt01_pt ->Write(C->nt01_pt ->GetName(), TObject::kWriteDelete);
				C->nt00_pt ->Write(C->nt00_pt ->GetName(), TObject::kWriteDelete);
				C->nt20_eta->Write(C->nt20_eta->GetName(), TObject::kWriteDelete);
				C->nt10_eta->Write(C->nt10_eta->GetName(), TObject::kWriteDelete);
				C->nt01_eta->Write(C->nt01_eta->GetName(), TObject::kWriteDelete);
				C->nt00_eta->Write(C->nt00_eta->GetName(), TObject::kWriteDelete);
				if(ch != EMu){
					C->fntight->Write(C->fntight->GetName(), TObject::kWriteDelete);
					C->fnloose->Write(C->fnloose->GetName(), TObject::kWriteDelete);
					C->pntight->Write(C->pntight->GetName(), TObject::kWriteDelete);
					C->pnloose->Write(C->pnloose->GetName(), TObject::kWriteDelete);
				}
			}
		}
	}
	pFile->Write();
	pFile->Close();
}
int  MuonPlotter::readHistos(TString filename){
	TFile *pFile = TFile::Open(filename, "READ");
	if(pFile == NULL){
		cout << "File " << filename << " does not exist!" << endl;
		return 1;
	}

	pFile->cd();
	if(gNSAMPLES != fSamples.size()){
		cout << "Mismatch in number of samples! Help!" << endl;
		return 1;
	}

	for(gSample i = sample_begin; i < gNSAMPLES; i=gSample(i+1)){
		Sample *S = &fSamples[i];
		for(gRegion r = region_begin; r < gNREGIONS; r=gRegion(r+1)){ // Loop over regions
			Region *R = &S->region[r];
			for(gChannel ch = channels_begin; ch < gNCHANNELS; ch=gChannel(ch+1)){ // Loop over channels, mumu, emu, ee
				Channel *C;
				if(ch == Muon)     C = &R->mm;
				if(ch == Electron) C = &R->ee;
				if(ch == EMu)      C = &R->em;
				TString rootdir = S->sname + "/" + R->sname + "/";
				C->nt20_pt  = (TH2D*)pFile->Get(rootdir + C->nt20_pt ->GetName() );
				C->nt10_pt  = (TH2D*)pFile->Get(rootdir + C->nt10_pt ->GetName() );
				C->nt01_pt  = (TH2D*)pFile->Get(rootdir + C->nt01_pt ->GetName() );
				C->nt00_pt  = (TH2D*)pFile->Get(rootdir + C->nt00_pt ->GetName() );
				C->nt20_eta = (TH2D*)pFile->Get(rootdir + C->nt20_eta->GetName() );
				C->nt10_eta = (TH2D*)pFile->Get(rootdir + C->nt10_eta->GetName() );
				C->nt01_eta = (TH2D*)pFile->Get(rootdir + C->nt01_eta->GetName() );
				C->nt00_eta = (TH2D*)pFile->Get(rootdir + C->nt00_eta->GetName() );
				if(ch != EMu){
					C->fntight = (TH2D*)pFile->Get(rootdir + C->fntight->GetName() );
					C->fnloose = (TH2D*)pFile->Get(rootdir + C->fnloose->GetName() );
					C->pntight = (TH2D*)pFile->Get(rootdir + C->pntight->GetName() );
					C->pnloose = (TH2D*)pFile->Get(rootdir + C->pnloose->GetName() );					
				}
			}
		}
		storeNumbers(S, Muon);
		storeNumbers(S, Electron);
		storeNumbers(S, EMu);
	}
	return 0;
}

//____________________________________________________________________________
void MuonPlotter::bookRatioHistos(){
	gStyle->SetOptStat(0);

	fH2D_MufRatio    = new TH2D("MufRatio",    "Ratio of tight to loose Muons vs Pt vs Eta", getNPt2Bins(Muon), getPt2Bins(Muon), getNEtaBins(Muon), getEtaBins(Muon));
	fH1D_MufRatioPt  = new TH1D("MufRatioPt",  "Ratio of tight to loose Muons vs Pt",        getNPt2Bins(Muon), getPt2Bins(Muon));
	fH1D_MufRatioEta = new TH1D("MufRatioEta", "Ratio of tight to loose Muons vs Eta",       getNEtaBins(Muon), getEtaBins(Muon));
	fH1D_MufRatioPt->SetXTitle(convertVarName("MuPt[0]"));
	fH1D_MufRatioEta->SetXTitle(convertVarName("MuEta[0]"));
	fH2D_MufRatio->SetXTitle(convertVarName("MuPt[0]"));
	fH2D_MufRatio->SetYTitle(convertVarName("MuEta[0]"));
	fH1D_MufRatioPt ->SetYTitle("# Tight / # Loose");
	fH1D_MufRatioEta->SetYTitle("# Tight / # Loose");
	fH1D_MufRatioPt->GetYaxis()->SetTitleOffset(1.2);
	fH1D_MufRatioEta->GetYaxis()->SetTitleOffset(1.2);
	
	fH2D_ElfRatio    = new TH2D("ElfRatio",    "Ratio of tight to loose Electrons vs Pt vs Eta", getNPt2Bins(Electron), getPt2Bins(Electron), getNEtaBins(Electron), getEtaBins(Electron));
	fH1D_ElfRatioPt  = new TH1D("ElfRatioPt",  "Ratio of tight to loose Electrons vs Pt",        getNPt2Bins(Electron), getPt2Bins(Electron));
	fH1D_ElfRatioEta = new TH1D("ElfRatioEta", "Ratio of tight to loose Electrons vs Eta",       getNEtaBins(Electron), getEtaBins(Electron));
	fH1D_ElfRatioPt->SetXTitle(convertVarName("ElPt[0]"));
	fH1D_ElfRatioEta->SetXTitle(convertVarName("ElEta[0]"));
	fH2D_ElfRatio->SetXTitle(convertVarName("ElPt[0]"));
	fH2D_ElfRatio->SetYTitle(convertVarName("ElEta[0]"));
	fH1D_ElfRatioPt ->SetYTitle("# Tight / # Loose");
	fH1D_ElfRatioEta->SetYTitle("# Tight / # Loose");
	fH1D_ElfRatioPt->GetYaxis()->SetTitleOffset(1.2);
	fH1D_ElfRatioEta->GetYaxis()->SetTitleOffset(1.2);
	
	fH2D_MupRatio    = new TH2D("MupRatio",    "Ratio of tight to loose Muons vs Pt vs Eta", getNPt2Bins(Muon), getPt2Bins(Muon), getNEtaBins(Muon), getEtaBins(Muon));
	fH1D_MupRatioPt  = new TH1D("MupRatioPt",  "Ratio of tight to loose Muons vs Pt",        getNPt2Bins(Muon), getPt2Bins(Muon));
	fH1D_MupRatioEta = new TH1D("MupRatioEta", "Ratio of tight to loose Muons vs Eta",       getNEtaBins(Muon), getEtaBins(Muon));
	fH1D_MupRatioPt->SetXTitle(convertVarName("MuPt[0]"));
	fH1D_MupRatioEta->SetXTitle(convertVarName("MuEta[0]"));
	fH2D_MupRatio->SetXTitle(convertVarName("MuPt[0]"));
	fH2D_MupRatio->SetYTitle(convertVarName("MuEta[0]"));
	fH1D_MupRatioPt ->SetYTitle("# Tight / # Loose");
	fH1D_MupRatioEta->SetYTitle("# Tight / # Loose");
	fH1D_MupRatioPt->GetYaxis()->SetTitleOffset(1.2);
	fH1D_MupRatioEta->GetYaxis()->SetTitleOffset(1.2);
	
	fH2D_ElpRatio    = new TH2D("ElpRatio",    "Ratio of tight to loose Electrons vs Pt vs Eta", getNPt2Bins(Electron), getPt2Bins(Electron), getNEtaBins(Electron), getEtaBins(Electron));
	fH1D_ElpRatioPt  = new TH1D("ElpRatioPt",  "Ratio of tight to loose Electrons vs Pt",        getNPt2Bins(Electron), getPt2Bins(Electron));
	fH1D_ElpRatioEta = new TH1D("ElpRatioEta", "Ratio of tight to loose Electrons vs Eta",       getNEtaBins(Electron), getEtaBins(Electron));
	fH1D_ElpRatioPt->SetXTitle(convertVarName("ElPt[0]"));
	fH1D_ElpRatioEta->SetXTitle(convertVarName("ElEta[0]"));
	fH2D_ElpRatio->SetXTitle(convertVarName("ElPt[0]"));
	fH2D_ElpRatio->SetYTitle(convertVarName("ElEta[0]"));
	fH1D_ElpRatioPt ->SetYTitle("# Tight / # Loose");
	fH1D_ElpRatioEta->SetYTitle("# Tight / # Loose");
	fH1D_ElpRatioPt->GetYaxis()->SetTitleOffset(1.2);
	fH1D_ElpRatioEta->GetYaxis()->SetTitleOffset(1.2);	
}
void MuonPlotter::fixPRatios(){
	// Checks if any bin of p ratio histo is empty and puts it to 1
	for(size_t i = 1; i <= fH2D_MupRatio->GetNbinsX(); ++i){
		for(size_t j = 1; j <= fH2D_MupRatio->GetNbinsY(); ++j){
			if(fH2D_MupRatio->GetBinContent(i,j) == 0.0){
				fH2D_MupRatio->SetBinContent(i,j, 1.0);
				fH2D_MupRatio->SetBinError(i,j, 0.5);
			}
		}
	}
	for(size_t i = 1; i <= fH2D_ElpRatio->GetNbinsX(); ++i){
		for(size_t j = 1; j <= fH2D_ElpRatio->GetNbinsY(); ++j){
			if(fH2D_ElpRatio->GetBinContent(i,j) == 0.0){
				fH2D_ElpRatio->SetBinContent(i,j, 1.0);
				fH2D_ElpRatio->SetBinError(i,j, 0.5);
			}
		}
	}
}

//////////////////////////////////////////////////////////////////////////////
// Event Selections:
//____________________________________________________________________________
bool MuonPlotter::isGoodEvent(){
	// Some global cuts, select events with >1 jets
	if(!passesNJetCut(2)) return false;
	return true;
}

//____________________________________________________________________________
bool MuonPlotter::isGoodMuEvent(){
	// Ask for >0 loose muons, if 2 muons ask for second to be loose too
	if(!isGoodEvent()) return false;
	if(NMus < 1) return false;
	if(isLooseMuon(0) == false) return false;
	if(NMus > 1) if(isLooseMuon(1) == false) return false;
	return true;
}

//____________________________________________________________________________
bool MuonPlotter::isGoodElEvent(){
	// Ask for >0 loose electrons, if 2 electrons ask for second to be loose too
	if(!isGoodEvent()) return false;
	if(NEls < 1) return false;
	if(isLooseElectron(0) == false) return false;
	if(NEls > 1) if(isLooseElectron(1) == false) return false;
	return true;
}

//____________________________________________________________________________
bool MuonPlotter::isGoodElMuEvent(){
	// Ask for >0 loose electrons and muons
	if(!isGoodEvent()) return false;
	if(NEls < 1 || NMus < 1) return false;
	for(size_t i = 0; i < NEls; ++i) if(!isLooseElectron(i)) return false;
	for(size_t i = 0; i < NMus; ++i) if(!isLooseMuon(i)) return false;
	// if(isLooseElectron(0) == false) return false;
	// if(isLooseMuon(0) == false) return false;
	return true;
}

//____________________________________________________________________________
int MuonPlotter::isSSLLEvent(int &ind1, int &ind2){
	// Looks for a SS loose/loose pair of leptons
	// Return the channel: 0 = none found
	//                     1 / -1 = mu+mu+ / mu-mu- pair
	//                     2 / -2 = e+e+   / e-e-   pair
	//                     3 / -3 = mu+e+  / mu-e-  pair
	// The indices in the argument given are sorted by pt unless
	// it's a e/mu event when they are sorted such that the muon
	// is ind1
	// The event selected is the one with hardest pt1 + pt2
	const float MMU = 0.1057;
	const float MEL = 0.0005;
	vector<lepton> tmp_looseLeps_p;
	vector<lepton> tmp_looseLeps_m;

	// First store all loose leptons in two vectors according to their charges
	for(size_t i = 0; i < NMus; ++i){
		if(!isLooseMuon(i)) continue;
		if(MuCharge[i] == 1 ){
			lepton tmpLepton;
			TLorentzVector pmu;
			pmu.SetPtEtaPhiM(MuPt[i], MuEta[i], MuPhi[i], gMMU);
			tmpLepton.p      = pmu;
			tmpLepton.charge = 1;
			tmpLepton.type   = 0;
			tmpLepton.index  = i;
			tmp_looseLeps_p.push_back(tmpLepton);
		}
		if(MuCharge[i] == -1){
			lepton tmpLepton;
			TLorentzVector p;
			p.SetPtEtaPhiM(MuPt[i], MuEta[i], MuPhi[i], gMMU);
			tmpLepton.p      = p;
			tmpLepton.charge = -1;
			tmpLepton.type   = 0;
			tmpLepton.index  = i;
			tmp_looseLeps_m.push_back(tmpLepton);
		}
	}
	for(size_t i = 0; i < NEls; ++i){
		if(!isLooseElectron(i)) continue;
		if(ElCh[i] == 1 ){
			lepton tmpLepton;
			TLorentzVector p;
			p.SetPtEtaPhiM(ElPt[i], ElEta[i], ElPhi[i], gMEL);
			tmpLepton.p      = p;
			tmpLepton.charge = 1;
			tmpLepton.type   = 1;
			tmpLepton.index  = i;
			tmp_looseLeps_p.push_back(tmpLepton);
		}
		if(ElCh[i] == -1){
			lepton tmpLepton;
			TLorentzVector p;
			p.SetPtEtaPhiM(ElPt[i], ElEta[i], ElPhi[i], gMEL);
			tmpLepton.p      = p;
			tmpLepton.charge = -1;
			tmpLepton.type   = 1;
			tmpLepton.index  = i;
			tmp_looseLeps_m.push_back(tmpLepton);
		}
	}

	// Sort these vectors by pt
	vector<lepton> looseLeps_p;
	vector<lepton> looseLeps_m;
	looseLeps_p = sortLeptonsByPt(tmp_looseLeps_p);
	looseLeps_m = sortLeptonsByPt(tmp_looseLeps_m);

	// Proceed to select one ss pair
	double ptmax(0.);
	int select(0);
	if(looseLeps_p.size() > 1){
		ptmax = looseLeps_p[0].p.Pt() + looseLeps_p[1].p.Pt();
		select = 1;
	}
	if(looseLeps_m.size() > 1){
		double ptsum = looseLeps_m[0].p.Pt() + looseLeps_m[1].p.Pt();		
		if(ptsum > ptmax){
			ptmax = ptsum;
			select = -1;
		}
		// if(looseLeps_p.size() > 1) cout << " Event with TWO SS pairs: r" << setw(7) << Run << "/e" << setw(13) << Event << "/l" << setw(5) << LumiSec << " in dataset " << setw(9) << fSamples[fCurrentSample].sname << endl;		
	}
	if(select == 0) return 0; // this ensures we have at least one pair
	
	vector<lepton> selectedPair;
	if(select == 1){ // positive
		selectedPair.push_back(looseLeps_p[0]);
		selectedPair.push_back(looseLeps_p[1]);
	}
	if(select == -1){ // negative
		selectedPair.push_back(looseLeps_m[0]);
		selectedPair.push_back(looseLeps_m[1]);
	}
	// Switch if el/mu combination (want ind1 to be mu, ind2 to be el)
	if(selectedPair[0].type == 1 && selectedPair[1].type == 0){
		lepton el = selectedPair[0];
		lepton mu = selectedPair[1];
		selectedPair[0] = mu;
		selectedPair[1] = el;
	}

	int result = 0;
	if(selectedPair[0].type == 0 && selectedPair[1].type == 0) result = 1; // mu/mu
	if(selectedPair[0].type == 1 && selectedPair[1].type == 1) result = 2; // el/el
	if(selectedPair[0].type == 0 && selectedPair[1].type == 1) result = 3; // mu/el
	result *= select; // Add charge to result
	
	// Result
	ind1 = selectedPair[0].index;
	ind2 = selectedPair[1].index;
	return result;
}

//____________________________________________________________________________
int MuonPlotter::isOSLLEvent(int &ind1, int &ind2){
	// Looks for a OS loose/loose pair of leptons
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
	const float MMU = 0.1057;
	const float MEL = 0.0005;
	vector<lepton> tmp_looseLeps;

	// First store all loose leptons in a vector
	for(size_t i = 0; i < NMus; ++i){
		if(!isLooseMuon(i)) continue;
			lepton tmpLepton;
			TLorentzVector pmu;
			pmu.SetPtEtaPhiM(MuPt[i], MuEta[i], MuPhi[i], gMMU);
			tmpLepton.p      = pmu;
			tmpLepton.charge = MuCharge[i];
			tmpLepton.type   = 0;
			tmpLepton.index  = i;
			tmp_looseLeps.push_back(tmpLepton);
	}
	for(size_t i = 0; i < NEls; ++i){
		if(!isLooseElectron(i)) continue;
			lepton tmpLepton;
			TLorentzVector p;
			p.SetPtEtaPhiM(ElPt[i], ElEta[i], ElPhi[i], gMEL);
			tmpLepton.p      = p;
			tmpLepton.charge = ElCh[i];
			tmpLepton.type   = 1;
			tmpLepton.index  = i;
			tmp_looseLeps.push_back(tmpLepton);
	}

	// Sort these vector by pt
	vector<lepton> looseLeps;
	looseLeps = sortLeptonsByPt(tmp_looseLeps);

	// Proceed to select one os pair
	if(looseLeps.size() < 2) return 0;
	
	vector<lepton> selectedPair;
	selectedPair.push_back(looseLeps[0]);
	for(size_t i = 1; i < looseLeps.size(); ++i){
		if(selectedPair[0].charge == looseLeps[i].charge) continue;
		selectedPair.push_back(looseLeps[i]);
		break;
	}
	if(selectedPair.size() < 2) return 0;

	// Switch if el/mu combination (want ind1 to be mu, ind2 to be el)
	if(selectedPair[0].type == 1 && selectedPair[1].type == 0){
		lepton el = selectedPair[0];
		lepton mu = selectedPair[1];
		selectedPair[0] = mu;
		selectedPair[1] = el;
	}

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
bool momentumComparator(MuonPlotter::lepton i, MuonPlotter::lepton j){ return (i.p.Pt()>j.p.Pt()); }
vector<MuonPlotter::lepton> MuonPlotter::sortLeptonsByPt(vector<lepton>& leptons){
	vector<lepton> theLep = leptons;
	sort (theLep.begin(), theLep.end(), momentumComparator);
	return theLep;
}

//____________________________________________________________________________
bool MuonPlotter::passesNJetCut(int cut){
	int njets(0);
	for(size_t i = 0; i < NJets; ++i) if(isGoodJet(i)) njets++;
	return njets>=cut;
}

//____________________________________________________________________________
bool MuonPlotter::passesNJetCut_LooseLep(int cut){
	// This is TIGHTER than passesNJetCut, so can be called on top of the other
	int njets(0);
	for(size_t i = 0; i < NJets; ++i) if(isGoodJet_LooseLep(i)) njets++;
	return njets>=cut;
}

//____________________________________________________________________________
bool MuonPlotter::passesHTCut(float cut){
	float ht(0.);
	for(size_t i = 0; i < NJets; ++i) if(isGoodJet(i)) ht += JetPt[i];
	return ht >= cut;
}

//____________________________________________________________________________
bool MuonPlotter::passesMETCut(float cut){
	if(pfMET < cut) return false;
	return true;
}

//____________________________________________________________________________
bool MuonPlotter::passesZVeto(float dm){
// Checks if any combination of opposite sign, same flavor tight leptons (e or mu)
// has invariant mass closer than dm to the Z mass, returns true if none found
// Default for dm is 15 GeV

	if(NMus > 1){
		unsigned start = 0;
		if(fCurrentChannel == Muon) start = 1; // For mumu, ignore first mu
		// First muon
		for(size_t i = start; i < NMus-1; ++i){
			if(isTightMuon(i)){
				TLorentzVector pmu1, pmu2;
				pmu1.SetPtEtaPhiM(MuPt[i], MuEta[i], MuPhi[i], gMMU);

				// Second muon
				for(size_t j = i+1; j < NMus; ++j){ 
					if(isTightMuon(j) && (MuCharge[i] != MuCharge[j]) ){
						pmu2.SetPtEtaPhiM(MuPt[j], MuEta[j], MuPhi[j], gMMU);
						if(fabs((pmu1+pmu2).M() - gMZ) < dm) return false;
					}
				}
			}
		}
	}
	
	if(NEls > 1){
		unsigned start = 0;
		if(fCurrentChannel == Electron) start = 1; // For ee, ignore first e
		// First electron
		for(size_t i = start; i < NEls-1; ++i){
			if(isTightElectron(i)){
				TLorentzVector pel1, pel2;
				pel1.SetPtEtaPhiM(ElPt[i], ElEta[i], ElPhi[i], gMEL);

				// Second electron
				for(size_t j = i+1; j < NEls; ++j){
					if(isTightElectron(j) && (ElCh[i] != ElCh[j]) ){
						pel2.SetPtEtaPhiM(ElPt[j], ElEta[j], ElPhi[j], gMEL);
						if(fabs((pel1+pel2).M() - gMZ) < dm) return false;
					}
				}
			}
		}		
	}
	return true;
}

//____________________________________________________________________________
bool MuonPlotter::passesMllEventVeto(float cut){
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
					if(isTightMuon(j)){/*
						TODO Check if they really use OS or also SS
					*/
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

//____________________________________________________________________________
bool MuonPlotter::isMuTriggeredEvent(){
	if(HLT_Mu9 == 0 &&
	   HLT_Mu11 == 0 &&
	   HLT_Mu13_v1 == 0 &&
	   HLT_Mu15 == 0 &&
	   HLT_Mu15_v1 == 0 &&
	   HLT_DoubleMu0 == 0 &&
	   HLT_DoubleMu3 == 0 &&
	   HLT_DoubleMu3_v2 == 0 &&
	   HLT_DoubleMu5_v2 == 0
	   ) return false;
	return true;
}

//____________________________________________________________________________
bool MuonPlotter::isElTriggeredEvent(){
	// Leptonic triggers from UCSD/UCSB/FNAL list
	if(Run == 1)                       return  HLT_Ele10_LW_L1R;
	if(Run >  1      && Run <  138000) return( HLT_Ele10_LW_L1R ||
	                                           HLT_Ele10_SW_L1R ||
	                                           HLT_Ele15_LW_L1R ||
	                                           HLT_DoubleEle5_SW_L1R );
	if(Run >= 138000 && Run <= 141900) return( HLT_Ele15_LW_L1R ||
	                                           HLT_Ele15_SW_L1R ||
	                                           HLT_Ele10_LW_EleId_L1R ||
	                                           HLT_DoubleEle5_SW_L1R);
	if(Run >  141900)                  return( HLT_Ele10_SW_EleId_L1R ||
	                                           HLT_Ele15_SW_CaloEleId_L1R ||
	                                           HLT_Ele15_SW_EleId_L1R ||
	                                           HLT_Ele17_SW_LooseEleId_L1R ||
	                                           HLT_Ele17_SW_CaloEleId_L1R ||
	                                           HLT_Ele17_SW_EleId_L1R ||
	                                           HLT_Ele17_SW_TightEleId_L1R ||
	                                           HLT_Ele17_SW_TighterEleId_L1R_v1 ||
	                                           HLT_Ele20_SW_L1R ||
	                                           HLT_Ele22_SW_TighterEleId_L1R_v2 ||
	                                           HLT_Ele22_SW_TighterEleId_L1R_v3 ||
	                                           HLT_Ele27_SW_TightCaloEleIdTrack_L1R_v1 ||
	                                           HLT_Ele32_SW_TightCaloEleIdTrack_L1R_v1 ||
	                                           HLT_Ele32_SW_TighterEleId_L1R_v2 ||
	                                           HLT_DoubleEle10_SW_L1R ||
	                                           HLT_DoubleEle15_SW_L1R_v1 ||
	                                           HLT_DoubleEle17_SW_L1R_v1 );
	return false;
}

//____________________________________________________________________________
bool MuonPlotter::isJetTriggeredEvent(){
	if(HLT_Jet15U == 0  && 
	   HLT_Jet30U == 0  && 
	   HLT_Jet50U == 0  &&
	   HLT_Jet70U == 0  &&
	   HLT_Jet100U == 0 &&
	   HLT_Jet100U_v2 == 0 &&
	   HLT_Jet100U_v3 == 0
	) return false;
	return true;
}

//____________________________________________________________________________
bool MuonPlotter::isHTTriggeredEvent(){
	// Selects run ranges and HT triggers
	if(Run == 1) if(HLT_HT100U == 1) return true;
	if(Run >= 140160 && Run <= 147116 ) if(HLT_HT100U    == 1) return true; // RunA
	if(Run >= 147196 && Run <= 148058 ) if(HLT_HT140U    == 1) return true; // Jet dataset
	if(Run >= 148822 && Run <= 149294 ) if(HLT_HT150U_v3 == 1) return true; // Multijet dataset
	if(!passesHTCut(300.)) return false;
	return false;
}

//____________________________________________________________________________
bool MuonPlotter::isGoodRun(int sample){
	// Select runs such that JetB and MultiJet datasets are mutually exclusive
	// if(gSample(sample) == JMB)      if(Run > 147195) return false;
	// if(gSample(sample) == MultiJet) if(Run < 147196) return false;
	if(gSample(sample) == JMB)      if(Run > 148058) return false;
	if(gSample(sample) == MultiJet) if(Run < 148822) return false;
	return true;
}

//____________________________________________________________________________
bool MuonPlotter::isSigSupMuEvent(){
	if(isGoodMuEvent() == false) return false;
	if(fSelectionSwitch == 1) if(!passesHTCut(300.)) return false;
	if(MuMT > 20.) return false;
	if(pfMET > 20.) return false;
	if(NMus > 1) return false;
	return true;
}

//____________________________________________________________________________
bool MuonPlotter::isSigSupMuEventTRG(){
	if(!isHTTriggeredEvent()) return false;
	// if(!isMuTriggeredEvent()) return false;
	if(!isSigSupMuEvent()) return false;
	return true;
}

//____________________________________________________________________________
bool MuonPlotter::isSigSupOSMuMuEvent(int& mu1, int& mu2){
	// This should include all the cuts for the final selection
	if(NMus < 1) return false;
	if(isLooseMuon(0) == false) return false;

	// if(!isGoodMuEvent()) return false; // >1 jets, >0 loose muons

	// if(MuMT  > 20.) return false;
	if(pfMET > 20.) return false;
	if(NMus < 2)    return false;    // >1 muons

	if(abs(isOSLLEvent(mu1, mu2)) != 1) return false;

	TLorentzVector pmu1, pmu2;
	pmu1.SetPtEtaPhiM(MuPt[mu1], MuEta[mu1], MuPhi[mu1], gMMU);
	pmu2.SetPtEtaPhiM(MuPt[mu2], MuEta[mu2], MuPhi[mu2], gMMU);
	if(fabs((pmu1+pmu2).M() - gMZ) < 25.) return false;

	return true;
}

//____________________________________________________________________________
bool MuonPlotter::isSigSupSSMuMuEvent(int& mu1, int& mu2){
	// This should include all the cuts for the final selection
	if(NMus < 1) return false;
	if(isLooseMuon(0) == false) return false;

	// if(!isGoodMuEvent()) return false; // >1 jets, >0 loose muons

	// if(MuMT  > 20.) return false;
	if(pfMET > 20.) return false;
	if(NMus < 2)    return false;    // >1 muons

	if(abs(isSSLLEvent(mu1, mu2)) != 1) return false;

	TLorentzVector pmu1, pmu2;
	pmu1.SetPtEtaPhiM(MuPt[mu1], MuEta[mu1], MuPhi[mu1], gMMU);
	pmu2.SetPtEtaPhiM(MuPt[mu2], MuEta[mu2], MuPhi[mu2], gMMU);
	if(fabs((pmu1+pmu2).M() - gMZ) < 25.) return false;

	return true;
}

//____________________________________________________________________________
bool MuonPlotter::isZMuMuEvent(){
	if(!passesNJetCut_LooseLep(2)) return false;
	if(isGoodMuEvent() == false) return false;
	if(NMus != 2) return false;
	if(!isLooseMuon(0) || !isLooseMuon(1)) return false; // both loose
	// if(!isTightMuon(0) && !isTightMuon(1)) return false; // at least one tight

	if(MuCharge[0] == MuCharge[1]) return false; // os

	// Z mass window cut
	TLorentzVector p1, p2;
	p1.SetPtEtaPhiM(MuPt[0], MuEta[0], MuPhi[0], 0.1057);
	p2.SetPtEtaPhiM(MuPt[1], MuEta[1], MuPhi[1], 0.1057);
	double m = (p1+p2).M();
	if(fabs(gMZ - m) > 15.) return false;

	if(pfMET > 20.) return false;

	return true;
}

//____________________________________________________________________________
bool MuonPlotter::isZMuMuEventTRG(){
	if(isMuTriggeredEvent() == false) return false;
	if(isZMuMuEvent() == false) return false;
	return true;
}

//____________________________________________________________________________
bool MuonPlotter::isZElElEvent(int &elind){
	if(!passesNJetCut_LooseLep(2)) return false;
	if(isGoodElEvent() == false) return false;
	if(NEls != 2) return false;
	if(!isLooseElectron(0) || !isLooseElectron(1)) return false; // both loose
	if(!isTightElectron(0) && !isTightElectron(1)) return false; // at least one tight

	if(ElCh[0] == ElCh[1]) return false; // os

	// Z mass window cut
	TLorentzVector p1, p2;
	p1.SetPtEtaPhiM(ElPt[0], ElEta[0], ElPhi[0], 0.0005);
	p2.SetPtEtaPhiM(ElPt[1], ElEta[1], ElPhi[1], 0.0005);
	double m = (p1+p2).M();
	if(fabs(gMZ - m) > 15.) return false;

	if(pfMET > 20.) return false;

	// If only the harder one tight or both tight, return the softer one
	// If only the softer one tight, return the harder one
	elind = 1;
	if(isTightElectron(1) && !isTightElectron(0)) elind = 0;
	return true;
}

//____________________________________________________________________________
bool MuonPlotter::isZElElEventTRG(){
	// Overload to be able to pass it as an argument
	int temp = 0;
	return isZElElEventTRG(temp);
}
bool MuonPlotter::isZElElEventTRG(int &elind){
	if(isElTriggeredEvent() == false) return false;
	if(isZElElEvent(elind) == false) return false;
	return true;
}

//____________________________________________________________________________
bool MuonPlotter::isSigSupElEvent(){
	if(!isGoodElEvent()) return false;
	if(fSelectionSwitch == 1) if(!passesHTCut(300.)) return false;
	if(ElMT[0] > 20.) return false;
	if(pfMET > 20.)   return false;
	if(NEls > 1)      return false;
	return true;
}

//____________________________________________________________________________
bool MuonPlotter::isSigSupElEventTRG(){
	// if(!isElTriggeredEvent()) return false;
	if(fSelectionSwitch == 0) if(!isElTriggeredEvent()) return false;
	if(fSelectionSwitch == 1) if(!isHTTriggeredEvent()) return false;
	if(!isSigSupElEvent()) return false;
	return true;
}

//____________________________________________________________________________
bool MuonPlotter::isGenMatchedSUSYDiLepEvent(){
	int ind1(-1), ind2(-1);
	return isGenMatchedSUSYDiLepEvent(ind1, ind2);
}
bool MuonPlotter::isGenMatchedSUSYDiLepEvent(int &mu1, int &mu2){
	if(!isGoodMuEvent()) return false;
	// if(isMuTriggeredEvent() == false) return false;
	if(!isSSTTMuEvent(mu1, mu2)) return false;
	if(isPromptSUSYMuon(mu1) && isPromptSUSYMuon(mu2)) return true;
	return false;
}

//____________________________________________________________________________
bool MuonPlotter::isGenMatchedSUSYEEEvent(){
	int ind1(-1), ind2(-1);
	if(!isSSTTElEvent(ind1, ind2)) return false;
	if(isPromptSUSYElectron(ind1) && isPromptSUSYElectron(ind2)) return true;
	return false;
}

//____________________________________________________________________________
bool MuonPlotter::isGenMatchedSUSYEMuEvent(){
	int muind(-1), elind(-1);
	if(!isSSTTElMuEvent(muind, elind)) return false;
	if(isPromptSUSYMuon(muind) && isPromptSUSYElectron(elind)) return true;
	return false;
}

//____________________________________________________________________________
bool MuonPlotter::isSSLLMuEvent(int& mu1, int& mu2){
	// This should include all the cuts for the final selection
	if(fDoCounting) fCounters[fCurrentSample][Muon].fill(" ... passes triggers");

	if(!isGoodMuEvent()) return false; // >1 jets, >0 loose muons
	if(fDoCounting) fCounters[fCurrentSample][Muon].fill(" ... has at least 2 jets, 1 loose muon");
	if(NMus < 2) return false;         // >1 muons
	if(fDoCounting) fCounters[fCurrentSample][Muon].fill(" ... has 2 loose muons");

	if(!passesZVeto()) return false; // no Zs in event
	if(fDoCounting) fCounters[fCurrentSample][Muon].fill(" ... passes Z veto");
	
	if(fSelectionSwitch == 0) if(!passesMllEventVeto(12.)) return false; // no low mass OSSF pairs
	if(fSelectionSwitch == 1) if(!passesMllEventVeto(5.) ) return false;
	if(fDoCounting) fCounters[fCurrentSample][Muon].fill(" ... passes Minv veto");
	
	if(fSelectionSwitch == 0) if(!passesHTCut(60.) )  return false;    // ht cut
	if(fSelectionSwitch == 1) if(!passesHTCut(300.))  return false;
	if(fDoCounting) fCounters[fCurrentSample][Muon].fill(" ... passes HT cut");
	
	if(fSelectionSwitch == 0) if(!passesMETCut(30.) ) return false;    // met cut
	if(fSelectionSwitch == 1) if(!passesMETCut(30.))  return false;
	if(fDoCounting) fCounters[fCurrentSample][Muon].fill(" ... passes MET cut");

	if(fChargeSwitch == 0){
		if(abs(isSSLLEvent(mu1, mu2)) != 1) return false;
		if(fDoCounting) fCounters[fCurrentSample][Muon].fill(" ... has same-sign muons");
	}
	if(fChargeSwitch == 1){
		if(abs(isOSLLEvent(mu1, mu2)) != 1) return false;
		if(fDoCounting) fCounters[fCurrentSample][Muon].fill(" ... has opposite-sign muons");		
	}

	if(!isGoodSecMuon(mu2)) return false; // pt cuts
	if(fDoCounting) fCounters[fCurrentSample][Muon].fill(" ... second muon passes pt cut");

	if(!isGoodPrimMuon(mu1)) return false;
	if(fDoCounting) fCounters[fCurrentSample][Muon].fill(" ... first muon passes pt cut");		

	return true;
}
bool MuonPlotter::isSSLLMuEventInvMETTRG(int& mu1, int& mu2){
	if(fSelectionSwitch == 0) if(!isMuTriggeredEvent()) return false;
	if(fSelectionSwitch == 1) if(!isHTTriggeredEvent()) return false;

	if(!isGoodMuEvent()) return false; // >1 jets, >0 loose muons
	if(NMus < 2) return false;         // >1 muons

	if(!passesZVeto()) return false; // no Zs in event
	
	if(fSelectionSwitch == 0) if(!passesMllEventVeto(12.)) return false; // no low mass OSSF pairs
	if(fSelectionSwitch == 1) if(!passesMllEventVeto(5.) ) return false;
	
	if(fSelectionSwitch == 0) if(!passesHTCut(60.) )  return false;    // ht cut
	if(fSelectionSwitch == 1) if(!passesHTCut(300.))  return false;
	
	if(fSelectionSwitch == 0) if(passesMETCut(30.) ) return false;    // INVERTED met cut
	if(fSelectionSwitch == 1) if(passesMETCut(30.))  return false;

	if(fChargeSwitch == 0) if(abs(isSSLLEvent(mu1, mu2)) != 1) return false;
	if(fChargeSwitch == 1) if(abs(isOSLLEvent(mu1, mu2)) != 1) return false;

	if(!isGoodSecMuon(mu2)) return false; // pt cuts
	if(!isGoodPrimMuon(mu1)) return false;
	return true;
}
bool MuonPlotter::isSSLLMuEventHTControlTRG(int& mu1, int& mu2){
	// only use for florida selection right now! there is no HT control region in UCSB selection
	if(fSelectionSwitch == 0) if(!isMuTriggeredEvent()) return false;
	if(fSelectionSwitch == 1) if(!isHTTriggeredEvent()) return false;

	if(!isGoodMuEvent()) return false; // >1 jets, >0 loose muons
	if(NMus < 2) return false;         // >1 muons

	if(!passesZVeto()) return false; // no Zs in event
	
	if(fSelectionSwitch == 0) if(!passesMllEventVeto(12.)) return false; // no low mass OSSF pairs
	if(fSelectionSwitch == 1) if(!passesMllEventVeto(5.) ) return false;
	
	if(fSelectionSwitch == 0) if(!passesHTCut(60.) )  return false;	// ht cut
	if(fSelectionSwitch == 1) if(!passesHTCut(250.))  return false;	// HT control region between 250 and 350
	if(fSelectionSwitch == 1) if( passesHTCut(350.))  return false;
	
	if(fSelectionSwitch == 0) if(!passesMETCut(30.) ) return false;
	if(fSelectionSwitch == 1) if(!passesMETCut(30.))  return false;

	if(fChargeSwitch == 0) if(abs(isSSLLEvent(mu1, mu2)) != 1) return false;
	if(fChargeSwitch == 1) if(abs(isOSLLEvent(mu1, mu2)) != 1) return false;

	if(!isGoodSecMuon(mu2)) return false; // pt cuts
	if(!isGoodPrimMuon(mu1)) return false;
	return true;
}
bool MuonPlotter::isSSLLMuEventTRG(int& mu1, int& mu2){
	if(fSelectionSwitch == 0) if(!isMuTriggeredEvent()) return false;
	if(fSelectionSwitch == 1) if(!isHTTriggeredEvent()) return false;
	if(!isSSLLMuEvent(mu1, mu2)) return false;
	return true;
}

//____________________________________________________________________________
bool MuonPlotter::isSSTTMuEvent(int& mu1, int& mu2){
	if(!isSSLLMuEvent(mu1, mu2)) return false;
	if(!isTightMuon(mu1) || !isTightMuon(mu2)) return false;
	return true;
}
bool MuonPlotter::isSSTTMuEventTRG(int& mu1, int& mu2){
	if(!isSSLLMuEventTRG(mu1, mu2)) return false;
	if(!isTightMuon(mu1) || !isTightMuon(mu2)) return false;
	return true;
}

//____________________________________________________________________________
bool MuonPlotter::isSSLLElEvent(int& el1, int& el2){
	// This should include all the cuts for the final selection
	if(fDoCounting) fCounters[fCurrentSample][Electron].fill(" ... passes triggers");

	if(!isGoodElEvent()) return false; // >1 jets, >0 loose eles
	if(fDoCounting) fCounters[fCurrentSample][Electron].fill(" ... has at least 2 jets, 1 loose el");
	if(NEls < 2) return false;         // >1 eles
	if(fDoCounting) fCounters[fCurrentSample][Electron].fill(" ... has 2 loose electrons");

	if(!passesZVeto()) return false; // no Zs in event
	if(fDoCounting) fCounters[fCurrentSample][Electron].fill(" ... passes Z veto");

	if(fSelectionSwitch == 0) if(!passesMllEventVeto(12.)) return false; // no low mass OSSF pairs
	if(fSelectionSwitch == 1) if(!passesMllEventVeto(5.) ) return false;
	if(fDoCounting) fCounters[fCurrentSample][Electron].fill(" ... passes Minv veto");

	if(fSelectionSwitch == 0) if(!passesHTCut(60.) )  return false;    // ht cut
	if(fSelectionSwitch == 1) if(!passesHTCut(300.))  return false;
	if(fDoCounting) fCounters[fCurrentSample][Electron].fill(" ... passes HT cut");

	if(fSelectionSwitch == 0) if(!passesMETCut(30.) ) return false;    // met cut
	if(fSelectionSwitch == 1) if(!passesMETCut(30.))  return false;
	if(fDoCounting) fCounters[fCurrentSample][Electron].fill(" ... passes MET cut");

	if(fChargeSwitch == 0){
		if(abs(isSSLLEvent(el1, el2)) != 2) return false;
		if(fDoCounting) fCounters[fCurrentSample][Electron].fill(" ... has same-sign electrons");
	}
	if(fChargeSwitch == 1){
		if(abs(isOSLLEvent(el1, el2)) != 2) return false;
		if(fDoCounting) fCounters[fCurrentSample][Electron].fill(" ... has opposite-sign electrons");
	}

	if(!isGoodSecElectron(el2)) return false; // pt cuts
	if(fDoCounting) fCounters[fCurrentSample][Electron].fill(" ... second electron passes pt cut");
	
	if(!isGoodPrimElectron(el1)) return false;
	if(fDoCounting) fCounters[fCurrentSample][Electron].fill(" ... first electron passes pt cut");
	
	return true;
}
bool MuonPlotter::isSSLLElEventTRG(int& el1, int& el2){
	if(fSelectionSwitch == 0) if(!isElTriggeredEvent()) return false;
	if(fSelectionSwitch == 1) if(!isHTTriggeredEvent()) return false;
	if(!isSSLLElEvent(el1, el2)) return false;
	return true;
}

//____________________________________________________________________________
bool MuonPlotter::isSSTTElEvent(int& el1, int& el2){
	if(!isSSLLElEvent(el1, el2)) return false;
	if(!isTightElectron(el1) || !isTightElectron(el2)) return false;
	return true;
}
bool MuonPlotter::isSSTTElEventTRG(int& el1, int& el2){
	if(!isSSLLElEventTRG(el1, el2)) return false;
	if(!isTightElectron(el1) || !isTightElectron(el2)) return false;
	return true;
}

//____________________________________________________________________________
bool MuonPlotter::isSSLLElMuEvent(int& mu, int& el){
	// This should include all the cuts for the final selection
	if(fDoCounting) fCounters[fCurrentSample][EMu].fill(" ... passes triggers");

	if(!isGoodElMuEvent()) return false; // >1 jets, >0 loose eles or muons
	if(fDoCounting) fCounters[fCurrentSample][EMu].fill(" ... has at least 1 j, loose e/mu pair");

	if(!passesZVeto())       return false;
	if(fDoCounting) fCounters[fCurrentSample][EMu].fill(" ... passes Z veto");

	if(fSelectionSwitch == 0) if(!passesMllEventVeto(12.)) return false; // no low mass OSSF pairs
	if(fSelectionSwitch == 1) if(!passesMllEventVeto(5.) ) return false;
	if(fDoCounting) fCounters[fCurrentSample][EMu].fill(" ... passes Minv veto");

	if(fSelectionSwitch == 0) if(!passesHTCut(60.) )  return false;    // ht cut
	if(fSelectionSwitch == 1) if(!passesHTCut(300.))  return false;
	if(fDoCounting) fCounters[fCurrentSample][EMu].fill(" ... passes HT cut");

	if(fSelectionSwitch == 0) if(!passesMETCut(20.) ) return false;    // met cut
	if(fSelectionSwitch == 1) if(!passesMETCut(30.))  return false;
	if(fDoCounting) fCounters[fCurrentSample][EMu].fill(" ... passes MET cut");

	if(fChargeSwitch == 0){
		if(abs(isSSLLEvent(mu, el)) != 3) return false;
		if(fDoCounting) fCounters[fCurrentSample][EMu].fill(" ... has same-sign electron muon pair");
	}
	if(fChargeSwitch == 1){
		if(abs(isOSLLEvent(mu, el)) != 3) return false;
		if(fDoCounting) fCounters[fCurrentSample][EMu].fill(" ... has opp.-sign electron muon pair");
	}

	if(MuPt[mu] > ElPt[el]){
		if(!isGoodPrimMuon(mu))    return false;
		if(fDoCounting) fCounters[fCurrentSample][EMu].fill(" ... muon passes pt cut");
		if(!isGoodSecElectron(el)) return false;
		if(fDoCounting) fCounters[fCurrentSample][EMu].fill(" ... electron passes pt cut");
	}
	else if(MuPt[mu] < ElPt[el]){
		if(!isGoodPrimElectron(el)) return false;
		if(fDoCounting) fCounters[fCurrentSample][EMu].fill(" ... electron passes pt cut");
		if(!isGoodSecMuon(mu))      return false;
		if(fDoCounting) fCounters[fCurrentSample][EMu].fill(" ... muon passes pt cut");
	}
	return true;
}
bool MuonPlotter::isSSLLElMuEventInvMETTRG(int& mu, int& el){
	if(fSelectionSwitch == 0){
		// Take all muon triggered events from muon datasets
		if(fCurrentSample == MuA || fCurrentSample == MuB) if(!isMuTriggeredEvent()) return false;
		// Take only electron but NOT muon triggered events from electron datasets
		if(fCurrentSample == EGA || fCurrentSample == EGB) if(!isElTriggeredEvent() ||  isMuTriggeredEvent()) return false;
		if(!isElTriggeredEvent() && !isMuTriggeredEvent()) return false;
	}
	if(fSelectionSwitch == 1) if(!isHTTriggeredEvent()) return false;

	if(!isGoodElMuEvent()) return false; // >1 jets, >0 loose eles or muons

	if(!passesZVeto())       return false;

	if(fSelectionSwitch == 0) if(!passesMllEventVeto(12.)) return false; // no low mass OSSF pairs
	if(fSelectionSwitch == 1) if(!passesMllEventVeto(5.) ) return false;

	if(fSelectionSwitch == 0) if(!passesHTCut(60.) )  return false;    // ht cut
	if(fSelectionSwitch == 1) if(!passesHTCut(300.))  return false;

	if(fSelectionSwitch == 0) if( passesMETCut(20.) ) return false;    // INVERTED MET cut
	if(fSelectionSwitch == 1) if( passesMETCut(30.))  return false;

	if(abs(isSSLLEvent(mu, el)) != 3) return false;

	if(MuPt[mu] > ElPt[el]){
		if(!isGoodPrimMuon(mu))    return false;
		if(!isGoodSecElectron(el)) return false;
	}
	else if(MuPt[mu] < ElPt[el]){
		if(!isGoodPrimElectron(el)) return false;
		if(!isGoodSecMuon(mu))      return false;
	}
	return true;
}
bool MuonPlotter::isOSLLElMuEventTRG(int& mu, int& el){
	if(fSelectionSwitch == 0){
		// Take all muon triggered events from muon datasets
		if(fCurrentSample == MuA || fCurrentSample == MuB) if(!isMuTriggeredEvent()) return false;
		// Take only electron but NOT muon triggered events from electron datasets
		if(fCurrentSample == EGA || fCurrentSample == EGB) if(!isElTriggeredEvent() ||  isMuTriggeredEvent()) return false;
		if(!isElTriggeredEvent() && !isMuTriggeredEvent()) return false;
	}
	if(fSelectionSwitch == 1) if(!isHTTriggeredEvent()) return false;

	if(!isGoodElMuEvent()) return false; // >1 jets, >0 loose eles or muons

	if(!passesZVeto())       return false;

	if(fSelectionSwitch == 0) if(!passesMllEventVeto(12.)) return false; // no low mass OSSF pairs
	if(fSelectionSwitch == 1) if(!passesMllEventVeto(5.) ) return false;

	if(fSelectionSwitch == 0) if(!passesHTCut(60.) )  return false;    // ht cut
	if(fSelectionSwitch == 1) if(!passesHTCut(300.))  return false;

	if(fSelectionSwitch == 0) if(!passesMETCut(20.) ) return false;
	if(fSelectionSwitch == 1) if(!passesMETCut(30.))  return false;

	if(abs(isOSLLEvent(mu, el)) != 3) return false;

	if(MuPt[mu] > ElPt[el]){
		if(!isGoodPrimMuon(mu))    return false;
		if(!isGoodSecElectron(el)) return false;
	}
	else if(MuPt[mu] < ElPt[el]){
		if(!isGoodPrimElectron(el)) return false;
		if(!isGoodSecMuon(mu))      return false;
	}
	return true;
}
bool MuonPlotter::isSSLLElMuEventTRG(int& mu, int& el){
	// For UCSD/SB/FNAL just use OR of all lepton triggers
	if(fSelectionSwitch == 0){
		// Take all muon triggered events from muon datasets
		if(fCurrentSample == MuA || fCurrentSample == MuB) if(!isMuTriggeredEvent()) return false;
		// Take only electron but NOT muon triggered events from electron datasets
		if(fCurrentSample == EGA || fCurrentSample == EGB) if(!isElTriggeredEvent() ||  isMuTriggeredEvent()) return false;
		if(!isElTriggeredEvent() && !isMuTriggeredEvent()) return false;
	}
	if(fSelectionSwitch == 1) if(!isHTTriggeredEvent()) return false;
	if(isSSLLElMuEvent(mu, el)) return true;
	return false;
}

//____________________________________________________________________________
bool MuonPlotter::isSSTTElMuEvent(int& mu, int& el){
	if(!isSSLLElMuEvent(mu, el)) return false;
	if(!isTightElectron(el) || !isTightMuon(mu)) return false;
	return true;
}
bool MuonPlotter::isSSTTElMuEventTRG(int& mu, int& el){
	/*
		TODO fixme: If this method is ever used, check trigger and dataset selections!
	*/
	if(!isSSLLElMuEventTRG(mu, el)) return false;
	if(!isTightElectron(el) || !isTightMuon(mu)) return false;
	return true;
}

//////////////////////////////////////////////////////////////////////////////
// Object selections:
//////////////////////////////////////////////////////////////////////////////
// Muons
//____________________________________________________________________________
bool MuonPlotter::isGoodMuon(int muon){
	if(muon >= NMus) return false; // Sanity check
	float ptcut(5.);
	if(fSelectionSwitch == 0) ptcut = 10.;
	if(fSelectionSwitch == 1) ptcut = 5.;
	if(MuPt[muon] < ptcut) return false;
	return true;
}

//____________________________________________________________________________
bool MuonPlotter::isLooseMuon(int muon){
	if(isGoodMuon(muon) == false)  return false;
	if(fSelectionSwitch == 0) if(MuIsoHybrid[muon] > 1.00) return false;
	if(fSelectionSwitch == 1) if(MuIso[muon]       > 1.00) return false;
	return true;
}

//____________________________________________________________________________
bool MuonPlotter::isTightMuon(int muon){
	if(isGoodMuon(muon) == false)  return false;
	if(isLooseMuon(muon) == false) return false;
	if(fSelectionSwitch == 0) if(MuIsoHybrid[muon] > 0.10) return false;
	if(fSelectionSwitch == 1) if(MuIso[muon]       > 0.15) return false;
	return true;
}

//____________________________________________________________________________
bool MuonPlotter::isGoodPrimMuon(int muon){
	if(isLooseMuon(muon) == false) return false;
	if(fSelectionSwitch == 0) if(MuPt[muon] < 20.) return false;
	if(fSelectionSwitch == 1) if(MuPt[muon] < 5. ) return false;
	return true;
}

//____________________________________________________________________________
bool MuonPlotter::isGoodSecMuon(int muon){
	if(isLooseMuon(muon) == false) return false;
	if(fSelectionSwitch == 0) if(MuPt[muon] < 10.) return false;
	if(fSelectionSwitch == 1) if(MuPt[muon] < 5. ) return false;
	return true;
}

//____________________________________________________________________________
bool MuonPlotter::isFakeTTbarMuon(int muon){
	if(isLooseMuon(muon) == false) return false;
	if(abs(MuGenMoID[muon]) == 24 || abs(MuGenMoID[muon]) == 15) return false;
	return true;
}

//____________________________________________________________________________
bool MuonPlotter::isPromptTTbarMuon(int muon){
	if(isLooseMuon(muon) == false) return false;
	if(abs(MuGenMoID[muon] == 24 && abs(MuGenGMoID[muon]) == 6))  return true;
	if(abs(MuGenMoID[muon] == 15 && abs(MuGenGMoID[muon]) == 24)) return true;
	return false;
}

//____________________________________________________________________________
bool MuonPlotter::isPromptSUSYMuon(int muon){
	if(isLooseMuon(muon) == false) return false;
	if( abs(MuGenMoType[muon]) == 9 || abs(MuGenMoType[muon]) == 4  || abs(MuGenMoType[muon]) == 2 ) return true;
	return false;
}

//////////////////////////////////////////////////////////////////////////////
// Electrons
//____________________________________________________________________________
bool MuonPlotter::isGoodElectron(int ele){
	if(ele >= NEls) return false; // Sanity check
	return true;
}

//____________________________________________________________________________
bool MuonPlotter::isLooseElectron(int ele){
	// All electrons are already loose in the high-pt selection (hybiso)
	if(isGoodElectron(ele) == false) return false;
	if(fSelectionSwitch == 1){
		if( fabs(ElEta[ele]) < 1.479 ) if(ElRelIso[ele] > 1.00) return false;
		else                           if(ElRelIso[ele] > 0.60) return false;		
	}
	if(ElChIsCons[ele] != 1) return false;
	return true;
}

//____________________________________________________________________________
bool MuonPlotter::isTightElectron(int ele){
	if(!isLooseElectron(ele))       return false;
	if(ElIsGoodElId_WP80[ele] != 1) return false;

	if(fSelectionSwitch == 0) if(ElHybRelIso[ele] > 0.10) return false;
	if(fSelectionSwitch == 1) if(ElRelIso[ele]    > 0.15) return false;
	
	return true;
}

//____________________________________________________________________________
bool MuonPlotter::isGoodPrimElectron(int ele){
	if(isLooseElectron(ele) == false) return false;
	if(fSelectionSwitch == 0) if(ElPt[ele] < 20.) return false;
	if(fSelectionSwitch == 1) if(ElPt[ele] < 10. ) return false;
	return true;
}

//____________________________________________________________________________
bool MuonPlotter::isGoodSecElectron(int ele){
	if(isLooseElectron(ele) == false) return false;
	if(fSelectionSwitch == 0) if(ElPt[ele] < 10.) return false;
	if(fSelectionSwitch == 1) if(ElPt[ele] < 10. ) return false;
	return true;
}

//____________________________________________________________________________
bool MuonPlotter::isPromptSUSYElectron(int ele){
	if(isLooseElectron(ele) == false) return false;
	if( abs(ElGenMType[ele]) == 9 || abs(ElGenMType[ele]) == 4  || abs(ElGenMType[ele]) == 2 ) return true;
	return false;
}

//____________________________________________________________________________
bool MuonPlotter::isGoodJet(int jet){
	if(jet >= NJets) return false; // Sanity check
	float minDR = 0.4;
	for(size_t imu = 0; imu < NMus; ++imu){
		if(!isTightMuon(imu)) continue;
		if(!isGoodSecMuon(imu)) continue;
		if(Util::GetDeltaR(MuEta[imu], JetEta[jet], MuPhi[imu], JetPhi[jet]) > minDR ) continue;
		return false;
	}
	for(size_t iel = 0; iel < NEls; ++iel){
		if(!isTightElectron(iel)) continue;
		if(!isGoodSecElectron(iel)) continue;
		if(Util::GetDeltaR(ElEta[iel], JetEta[jet], ElPhi[iel], JetPhi[jet]) > minDR ) continue;
		return false;
	}
	return true;
}

//____________________________________________________________________________
bool MuonPlotter::isGoodJet_LooseLep(int jet){
	if(jet >= NJets) return false; // Sanity check
	float minDR = 0.4;
	for(size_t imu = 0; imu < NMus; ++imu){
		if(!isLooseMuon(imu)) continue;
		if(!isGoodSecMuon(imu)) continue;
		if(Util::GetDeltaR(MuEta[imu], JetEta[jet], MuPhi[imu], JetPhi[jet]) > minDR ) continue;
		return false;
	}
	for(size_t iel = 0; iel < NEls; ++iel){
		if(!isLooseElectron(iel)) continue;
		if(!isGoodSecElectron(iel)) continue;
		if(Util::GetDeltaR(ElEta[iel], JetEta[jet], ElPhi[iel], JetPhi[jet]) > minDR ) continue;
		return false;
	}
	return true;
}
