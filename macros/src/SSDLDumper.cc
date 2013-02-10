/******************************************************************************
*   Collection of tools for producing plots for same-sign dilepton analysis  *
*                                                                            *
*                                                  (c) 2010 Benjamin Stieger *
******************************************************************************
* This runs on SSDLTrees, applies all object and event selections for signal *
* and control regions of the SSDL analysis, and stores the relevant numbers  *
* in a set of histograms in a ROOT file.                                     *
* It also writes a small tree contraining only same-sign events for quick    *
* browsing of signal candidate events                                        *
*****************************************************************************/
#include "SSDLDumper.hh"

#include "helper/AnaClass.hh"
#include "helper/Utilities.hh"
#include "helper/FPRatios.hh"
#include "helper/Davismt2.h"
#include "helper/FakeRatios.hh"
#include "helper/Monitor.hh"

//#include "mcbtagSFuncert.h" // claudios btag scales
//
//#include "helper/BTagSFUtil/BTagSFUtil.h"
// #include "helper/BTagSF.hh"
// #include "helper/GoodRunList.h"

#include "TLorentzVector.h"
#include "TGraphAsymmErrors.h"
#include "TEfficiency.h"
#include "TRandom3.h"

#include "TWbox.h"
#include "TMultiGraph.h"
#include "TGaxis.h"

#include <iostream>
#include <iomanip>
#include <time.h> // access to date/time

#include <typeinfo>

int gDEBUG_EVENTNUMBER_ = -1;  
int gDEBUG_RUNNUMBER_ = -1;

using namespace std;

//////////////////////////////////////////////////////////////////////////////////
// Global parameters:

static const bool  gApplyTauVeto = true;
static       bool  gSmearMET     = false;

bool gDoSystStudies;
static const bool gDoSyncExercise = false;
float gMuMaxIso     ;
float gElMaxIso     ;
float gMinJetPt     ;
float gMaxJetEta    ;
bool  gTTWZ         ;
//bool  gApplyZVeto   ;
bool  gInvertZVeto  ;
bool  tmp_gApplyZVeto;
bool  gApplyGStarVeto    ;
TString tmp_gBaseRegion  ;
TString gJSONfile        ;
bool  tmp_gDoWZValidation;

// std::vector< SSDLDumper::Region* > gRegions;
// std::vector< SSDLDumper::Region* >::iterator regIt;
// std::map<TString , int> gRegion;
// int gNREGIONS;


//////////////////////////////////////////////////////////////////////////////////
static const float gMMU = 0.1057;
static const float gMEL = 0.0005;
static const float gMZ  = 91.;

// Muon Binning //////////////////////////////////////////////////////////////////
double SSDLDumper::gMuFPtBins[gNMuFPtBins+1] = {20., 25., 30., 35., 40., 50., 60.}; // fake ratios
double SSDLDumper::gMuPPtbins[gNMuPPtbins+1] = {20., 25., 30., 35., 40., 50., 60., 70., 80., 90., 100.}; // prompt ratios
double SSDLDumper::gMuEtabins[gNMuEtabins+1] = {0., 0.5, 1.0, 1.479, 2.0, 2.5};

// Electron Binning //////////////////////////////////////////////////////////////
double SSDLDumper::gElFPtBins[gNElFPtBins+1]   = {20., 25., 30., 40., 50., 60., 70., 80., 100.}; // fake ratios
double SSDLDumper::gElPPtbins[gNElPPtbins+1]   = {20., 25., 30., 35., 40., 50., 60., 70., 80., 90., 100.}; // prompt ratios
double SSDLDumper::gElEtabins[gNElEtabins+1]   = {0., 0.5, 1.0, 1.479, 2.0, 2.5};
double SSDLDumper::gElCMIdbins[gNElCMIdbins+1] = {0., 1.479, 2.5};
//////////////////////////////////////////////////////////////////////////////////

// NVrtx Binning /////////////////////////////////////////////////////////////////
double SSDLDumper::gNVrtxBins[gNNVrtxBins+1]  = {5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25};
//OLD double SSDLDumper::gNVrtxBins[gNNVrtxBins+1]  = {0, 2, 4, 6, 8, 10, 12, 14, 16, 18};
//////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////
TString SSDLDumper::gKinSelNames[gNKinSels] = {"LL", "TT", "Sig"};
TString SSDLDumper::KinPlots::var_name[SSDLDumper::gNKinVars] = {"HT", "MET", "NJets", "Pt1", "Pt2", "InvMassSF", "InvMassMM", "InvMassEE", "InvMassEM", "MT2", "NbJets", "NbJetsMed"};
int     SSDLDumper::KinPlots::nbins[SSDLDumper::gNKinVars]    = {  20 ,  20 ,      8 ,   14 ,   14 ,        14  ,        14  ,        14  ,        14  ,   20 ,       5 ,          5 };
float   SSDLDumper::KinPlots::xmin[SSDLDumper::gNKinVars]     = {   0.,  30.,      0.,   20.,   20.,        20. ,        20. ,        20. ,        20. ,    0.,       0.,          0.};
float   SSDLDumper::KinPlots::xmax[SSDLDumper::gNKinVars]     = {1000., 330.,      8.,  300.,  160.,       300. ,       300. ,       300. ,       300. ,  100.,       5.,          5.};
TString SSDLDumper::KinPlots::axis_label[SSDLDumper::gNKinVars] = {"H_{T} [GeV]",
                                                                   "Particle Flow ME_{T} [GeV]",
                                                                   "Jet Multiplicity",
                                                                   "Leading Lepton p_{T} [GeV]",
                                                                   "Subleading Lepton p_{T} [GeV]",
                                                                   "m_{ll} (SF) [GeV]",
                                                                   "m_{#mu#mu} [GeV]",
                                                                   "m_{ee} [GeV]",
                                                                   "m_{ll} (OF) [GeV]",
                                                                   "M_{T2} [GeV]",
                                                                   "b-Jet Multiplicity (loose)",
                                                                   "b-Jet Multiplicity (medium)"};

//////////////////////////////////////////////////////////////////////////////////
double SSDLDumper::gDiffHTBins  [gNDiffHTBins+1]   = { 80., 120., 200., 320., 400., 500., 600.};
double SSDLDumper::gDiffMETBins [gNDiffMETBins+1]  = { 0.,   30.,  50.,  70.,  80., 100., 120.};
double SSDLDumper::gDiffNJBins  [gNDiffNJBins+1]   = { 0.,    1.,   2.,   3.,   4.,   5.,   6.}; // fill NJets + 0.5 to hit the right bin
double SSDLDumper::gDiffMT2Bins [gNDiffMT2Bins+1]  = { 0.,   25.,  50., 100.                  };
double SSDLDumper::gDiffPT1Bins [gNDiffPT1Bins+1]  = { 20., 40., 60., 80., 100., 120., 140., 160., 180., 200.};
double SSDLDumper::gDiffPT2Bins [gNDiffPT2Bins+1]  = { 20., 30., 40., 50.,  60.,  70.,  80.,  90., 100.      };
double SSDLDumper::gDiffNBJBins [gNDiffNBJBins+1]  = { 0., 1., 2., 3., 4.};
double SSDLDumper::gDiffNBJMBins[gNDiffNBJMBins+1] = { 0., 1., 2., 3.    };
double SSDLDumper::gDiffMET3Bins[gNDiffMET3Bins+1] = {30., 40., 50., 60., 80., 100., 120., 160., 200., 250.};

//////////////////////////////////////////////////////////////////////////////////
TString SSDLDumper::DiffPredYields::var_name  [SSDLDumper::gNDiffVars] = {"HT", "MET", "NJets", "MT2", "PT1", "PT2", "NBJets", "MET3", "NBJetsMed", "NBJetsMed2"};
int     SSDLDumper::DiffPredYields::nbins     [SSDLDumper::gNDiffVars] = {gNDiffHTBins, gNDiffMETBins, gNDiffNJBins, gNDiffMT2Bins, gNDiffPT1Bins, gNDiffPT2Bins, gNDiffNBJBins, gNDiffMET3Bins, gNDiffNBJMBins, gNDiffNBJMBins};
double* SSDLDumper::DiffPredYields::bins      [SSDLDumper::gNDiffVars] = { gDiffHTBins,  gDiffMETBins,  gDiffNJBins,  gDiffMT2Bins,  gDiffPT1Bins,  gDiffPT2Bins,  gDiffNBJBins,  gDiffMET3Bins,  gDiffNBJMBins,  gDiffNBJMBins};
TString SSDLDumper::DiffPredYields::axis_label[SSDLDumper::gNDiffVars] = {"H_{T} [GeV]",
                                                                          "Particle Flow ME_{T} [GeV]",
                                                                          "Jet Multiplicity",
                                                                          "M_{T2} [GeV]",
                                                                          "Leading Lepton p_{T} [GeV]",
                                                                          "Subleading Lepton p_{T} [GeV]",
                                                                          "b-Jet Multiplicity (loose)",
                                                                          "Particle Flow ME_{T} [GeV]",
                                                                          "b-Jet Multiplicity (medium)",
                                                                          "b-Jet Multiplicity (medium)"};

//////////////////////////////////////////////////////////////////////////////////
TString SSDLDumper::FRatioPlots::var_name[SSDLDumper::gNRatioVars] = {"NJets",  "HT", "MaxJPt", "NVertices", "ClosJetPt", "AwayJetPt", "NBJets", "MET",  "MT"};
int     SSDLDumper::FRatioPlots::nbins[SSDLDumper::gNRatioVars]    = {     7 ,   10 ,      10 ,        10  ,        10  ,        10  ,       3 ,   15 ,   10 };
float   SSDLDumper::FRatioPlots::xmin[SSDLDumper::gNRatioVars]     = {     1.,   50.,      30.,         5. ,        30. ,        30. ,       0.,    0.,    0.};
float   SSDLDumper::FRatioPlots::xmax[SSDLDumper::gNRatioVars]     = {     8.,  500.,     300.,        25. ,       150. ,       300. ,       3.,  150.,  100.};

//////////////////////////////////////////////////////////////////////////////////
TString SSDLDumper::IsoPlots::sel_name[SSDLDumper::gNSels] = {"Base", "SigSup"};
int     SSDLDumper::IsoPlots::nbins[SSDLDumper::gNSels]    = {12, 12};
TString SSDLDumper::IdPlots::sel_name[SSDLDumper::gNSels] = {"Base", "SigSup"};
int     SSDLDumper::IdPlots::nbins[SSDLDumper::gNSels]    = {20, 20};

TString SSDLDumper::gEMULabel[2] = {"mu", "el"};
TString SSDLDumper::gChanLabel[3] = {"MM", "EM", "EE"}; // make SURE this is the same order as gChannel enum!
TString SSDLDumper::gHiLoLabel[3] = {"HighPt", "LowPt", "TauChan"};

void setVariables(char buffer[1000]){
	char va[1], t[100], n[100], val[100];
	TString v, type, name, value;
	if( sscanf(buffer, "%s\t%s\t%s\t%s", va, t, n, val) > 3){
		v = va; type = t; name = n; value = val;
		if (v != "v") {cout << "ERROR in reading variables!!" << endl; exit(1); }
		if      (type == "bool"    && name =="gDoSystStudies" ) gDoSystStudies   = ((value == "1" || value == "true") ? true:false);
		else if (type == "bool"    && name =="gTTWZ"          ) gTTWZ            = ((value == "1" || value == "true") ? true:false);
		else if (type == "bool"    && name =="gInvertZVeto"   ) gInvertZVeto        = ((value == "1" || value == "true") ? true:false);
		else if (type == "bool"    && name =="gApplyGStarVeto") gApplyGStarVeto     = ((value == "1" || value == "true") ? true:false);
		else if (type == "TString" && name =="gBaseRegion"    ) tmp_gBaseRegion     = value; // this and the next are the only ones used in the plotter
		else if (type == "TString" && name =="gJSONfile"      ) gJSONfile           = value;
		else if (type == "bool"    && name =="gApplyZVeto"   ) tmp_gApplyZVeto  = ((value == "1" || value == "true") ? true:false);
		else if (type == "float"   && name =="gMuMaxIso"     ) gMuMaxIso        = value.Atof();
		else if (type == "float"   && name =="gElMaxIso"     ) gElMaxIso        = value.Atof();
		else if (type == "float"   && name =="gMaxJetEta"    ) gMaxJetEta       = value.Atof();
		else if (type == "float"   && name =="gMinJetPt"     ) gMinJetPt        = value.Atof();
		else {cout << "ERROR in reading variables!!" << endl; exit(1); }
	}
	else{
		cout << " SSDLDumper::setVariables ==> Wrong configfile format for special variables. Have a look in the dumperconfig.cfg for more info. Aborting..." << endl;
		exit(1);
	}

}
//____________________________________________________________________________
SSDLDumper::SSDLDumper(TString configfile){
	char buffer[1000];
	ifstream IN(configfile);
	if ( !(IN.is_open()) ) {
		cout << "ERROR detected: dumper configuration file " << configfile << " could not be opened! Are you sure it exists??" << endl;
		cout << " ... exiting ..." << endl;
		exit(1);
	}

	// TString name;
	float miMu1, miMu2, miEl1, miEl2;
	float miHT, maHT, miMET, maMET, miJetPt  ;
	int   miNj, maNj, miNb, maNb, miNbm, maNbm, ve3rd, veTTZ, veCha, veGStar;

	char name[100];
	cout << "====================================" << endl;
	cout << "Reading Configuration file " << configfile << endl;
	int counter(0);

	while( IN.getline(buffer, 700, '\n') ){
		if (buffer[0] == '#') continue; // Skip lines commented with '#'
		if (strlen(buffer) == 0)  continue; // Skip empty lines
		if (buffer[0] == 'v') { //set the special variables from the config file
			setVariables(buffer);
			continue;
		}
		if( sscanf(buffer, "%s\t%f\t%f\t%f\t%f\t%d\t%d\t%d\t%d\t%d\t%d\t%f\t%f\t%f\t%f\t%d\t%d\t%d", 
			   name, &miHT, &maHT, &miMET, &maMET, &miNj, &maNj, &miNb, &maNb, &miNbm, &maNbm, &miMu1, &miMu2, &miEl1, &miEl2, &ve3rd, &veTTZ, &veCha) > 17){
               //             name, &miHT, &maHT, &miMET, &maMET, &miNj, &maNj, &miNb, &maNb, &miNbm, &maNbm, &miMu1, &miMu2, &miEl1, &miEl2, &ve3rd, &veTTZ, &veCha) > 17){
			SSDLDumper::Region * tmp_region  = new SSDLDumper::Region();
			tmp_region->sname      = (TString) name;
			tmp_region->minHT      = miHT;
			tmp_region->maxHT      = maHT;
			tmp_region->minMet     = miMET;
			tmp_region->maxMet     = maMET;
			tmp_region->minNjets   = miNj;
			tmp_region->maxNjets   = maNj;
			tmp_region->minNbjets  = miNb;
			tmp_region->maxNbjets  = maNb;
			tmp_region->minNbjmed  = miNbm;
			tmp_region->maxNbjmed  = maNbm;
			tmp_region->minMu1pt   = miMu1;
			tmp_region->minMu2pt   = miMu2;
			tmp_region->minEl1pt   = miEl1;
			tmp_region->minEl2pt   = miEl2;
			tmp_region->app3rdVet  = ve3rd;
			tmp_region->vetoTTZSel = veTTZ;
			tmp_region->chargeVeto = veCha;
			//			tmp_region->GStarVeto  = veGStar;
			SSDLDumper::gRegion[name] = counter;
			counter++;


			SSDLDumper::gRegions.push_back(tmp_region);
			if(fVerbose > 0) {
				cout << "----------------  LOADING REGION  ------------------------------" << endl;
				cout << Form("%13s/%3.0f/%4.0f/%3.0f/%3.0f/%3.0f/%3.0f",
				tmp_region->sname.Data()       ,
				tmp_region->minHT       ,
				tmp_region->maxHT       ,
				tmp_region->minMet      ,
				tmp_region->maxMet      ,
				tmp_region->minMu1pt    ,
				tmp_region->minMu2pt  ) << endl;
				cout << "--------------------------------------------------------" << endl;
			}
			//delete tmp_region;
		}
		else{
			cout << " SSDLDumper::readConfig ==> Wrong configfile format in file " << configfile << " detected, aborting..." << endl;
			exit(1);
		}


	}
	SSDLDumper::gNREGIONS = SSDLDumper::gRegions.size();
	IN.close();

	cout << "================  GLOBAL PARAMETERS  ===================" << endl;
	cout << Form("gTTWZ = %d\ngApplyZVeto = %d\ngInvertZVeto = %d\ngMuMaxIso = %3.2f\ngElMaxIso = %3.2f\ngMinJetPt = %3.2f\ngMaxJetEta = %3.2f\ngBaseRegion = %s",
		gTTWZ       ,
		gApplyZVeto ,
		gInvertZVeto,
		gMuMaxIso   ,
		gElMaxIso   ,
		gMinJetPt   ,
		gMaxJetEta   ,
		tmp_gBaseRegion.Data()  ) << endl;
	cout << "========================================================" << endl;
	SSDLDumper::gBaseRegion = tmp_gBaseRegion;
	SSDLDumper::gApplyZVeto = tmp_gApplyZVeto;
	SSDLDumper::gDoWZValidation = tmp_gDoWZValidation;

	// initializing all the systematics here out of a lack of other places
	gSystematics["Normal"]   = 0;
	gSystematics["JetUp"]    = 1;
	gSystematics["JetDown"]  = 2;
	gSystematics["JetSmear"] = 3;
	gSystematics["BUp"]      = 4;
	gSystematics["BDown"]    = 5;
	gSystematics["LepUp"]    = 6;
	gSystematics["LepDown"]  = 7;

}

//____________________________________________________________________________
SSDLDumper::~SSDLDumper(){
	if(fOutputFile != NULL && fOutputFile->IsOpen()) fOutputFile->Close();
	fChain = 0;
}

//____________________________________________________________________________
void SSDLDumper::init(TString inputfile, TString sname, int datamc, float xsec, int chan){
	fSamples.push_back(new Sample(inputfile, sname, datamc, xsec, gNREGIONS, chan));

	if(fVerbose > 0) cout << "------------------------------------" << endl;
	if(fVerbose > 0) cout << " Initializing SSDLDumper ... " << endl;
	if(fVerbose > 0) cout << "   Running on:      " << fSamples.back()->location << endl;
	if(fVerbose > 0) cout << "   Naming it:       " << fSamples.back()->sname << endl;
	if(fVerbose > 0) cout << "   Is data/mc:      " << fSamples.back()->datamc << endl;
	if(fVerbose > 0) cout << "   x-section is:    " << fSamples.back()->xsec << endl;
	if(fVerbose > 0) cout << "   Int. Lumi is:    " << fSamples.back()->getLumi() << endl;
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

	// fBTagSFUtil = new BTagSFUtil("CSV", 28);
	fBTagSF = new BTagSF();
	fRand3 = new TRandom3(0);

	resetHypLeptons();
	initCutNames();
	
	setRegionCuts(gRegion[gBaseRegion]); // no argument = reset to gBaseRegion

 	fC_maxMet_Control = 20.;
	fC_maxMt_Control  = 20.;

	// Prevent root from adding histograms to current file
	TH1::AddDirectory(kFALSE);
}


void SSDLDumper::readDatacard(TString cardfile){
	char buffer[1000];
	ifstream IN(cardfile);

	char inputfile[500], sname[100];
	int datamc, chan;
	float lumi;
	float xsec;
	if(fVerbose > 2) cout << "------------------------------------" << endl;
	if(fVerbose > 2) cout << "Reading datacard  " << cardfile << endl;
	int counter(0);

	while( IN.getline(buffer, 600, '\n') ){
		lumi = 1.0;
		xsec = 1.0;
		if (buffer[0] == '#') continue; // Skip lines commented with '#'
		if( sscanf(buffer, "%s\t%s\t%d\t%d\t%f", sname, inputfile, &datamc, &chan, &xsec) > 3){

			fSamples.push_back(new Sample(inputfile, sname, datamc, xsec, gNREGIONS, chan));
			fSampleMap[sname] = fSamples.back(); // S;
			if(fVerbose > 1 && fVerbose < 3){
				cout << Form("%13s/%2d/%2d/%12.2f - %s",
				fSamples.back()->sname.Data()   ,
				fSamples.back()->datamc         ,
				fSamples.back()->chansel        ,
				fSamples.back()->xsec           ,
				fSamples.back()->location.Data()) << endl;
			}
			if(fVerbose > 2){
				cout << " ---- " << endl;
				cout << "  New sample added: " << fSamples.back()->sname     << endl;
				cout << "   Input:      "      << fSamples.back()->location  << endl;
				cout << "   DataMC:     "      << fSamples.back()->datamc    << endl;
				cout << "   Channel:    "      << fSamples.back()->chansel   << endl;
				cout << "   XSec:       "      << fSamples.back()->xsec      << endl;
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
const int     SSDLDumper::getNFPtBins(gChannel chan){ // fake ratios
	if(chan == Muon || chan == ElMu) return gNMuFPtBins;
	if(chan == Elec) return gNElFPtBins;
}
const double *SSDLDumper::getFPtBins (gChannel chan){
	if(chan == Muon || chan == ElMu) return gMuFPtBins;
	if(chan == Elec) return gElFPtBins;
}
const int     SSDLDumper::getNPPtBins(gChannel chan){ // prompt ratios
	if(chan == Muon || chan == ElMu) return gNMuPPtbins;
	if(chan == Elec) return gNElPPtbins;
}
const double *SSDLDumper::getPPtBins (gChannel chan){
	if(chan == Muon || chan == ElMu) return gMuPPtbins;
	if(chan == Elec) return gElPPtbins;
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
        /////////////////////// 
        //READ JSON FILE
        /////////////////////// 
        fGoodRunList = new GoodRunList(gJSONfile);
	//////////////////////
	for(size_t i = 0; i < fSamples.size(); ++i){
		fSample = fSamples[i]; // TODO: Clean this up, call the triggers with an argument
		fOutputFileName = fOutputDir + fSample->sname + "_Yields.root";
		
		fCurLumi = -1;
		fCurRun  = -1;
		skipRun  = false;
		skipLumi = false;
		
		loopEvents(fSample);
	}

	delete fGoodRunList;
}
void SSDLDumper::loopEvents(Sample *S){
	fDoCounting = true;
	if(S->datamc == 0){
		TString eventfilename  = fOutputDir + S->sname + "_SignalEvents.txt";
		fOUTSTREAM.open(eventfilename.Data(), ios::trunc);		
	}

	TFile *pFile = new TFile(fOutputFileName, "RECREATE");

	bookHistos(S);
	bookSigEvTree();
	
	TTree *tree = S->getTree();
	S->evcount->Add(S->getEvCount());	
	//cout << "MARC: evcount: " << S->getEvCount()->GetEntries() << endl;
	// MARC nov 29 S->ngen = S->getEvCount()->GetEntries();
	
	// Stuff to execute for each sample BEFORE looping on the events
	initCounters();

	// Event loop
	tree->ResetBranchAddresses();
	// Init(tree);
	Init(tree);

	if (fChain == 0) return;
	Long64_t nentries = fChain->GetEntriesFast();
	Long64_t nbytes = 0, nb = 0;
	for (Long64_t jentry=0; jentry<nentries;jentry++) {
		printProgress(jentry, nentries, S->sname);

		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		nb = fChain->GetEntry(jentry);   nbytes += nb;

		
		/////////////////////////////////////////////
		//   APPLY JSON
		/////////////////////////////////////////////
		if (S->datamc == 0 && !IsInJSON()) continue;
		
		// if (S->datamc == 0 && Run > 200991) continue; // ATTENTION HERE, THIS IS JUST FOR TTWZ AND TEMPORARY
		/////////////////////////////////////////////
		// DEBUG
		//		if (Run!=gDEBUG_RUNNUMBER_) continue;
		// if(!(Event==gDEBUG_EVENTNUMBER_ && Run==gDEBUG_RUNNUMBER_)) continue;
		// if(jentry > 10000) break;
		/////////////////////////////////////////////
		
		/////////////////////////////////////////////
		// REJECT DOUBLE COUNTED EVENTS
		/////////////////////////////////////////////
		// double counting of fakes if (S->datamc!=0 && AvoidDoubleCountingOfFakes(S)) continue;

		/////////////////////////////////////////////
		// Event modifications
		scaleBTags(S, 0); // this applies the bTagSF which comes from pandolfis function
		saveBTags(); // this just saves the btag values for each jet.
		/////////////////////////////////////////////
		
		fCounter[Muon].fill(fMMCutNames[0]);
		fCounter[ElMu].fill(fEMCutNames[0]);
		fCounter[Elec].fill(fEECutNames[0]);

		// Select mutually exclusive runs for Jet and MultiJet datasets
		if(!isGoodRun(S)) continue;

		fCounter[Muon].fill(fMMCutNames[1]);
		fCounter[ElMu].fill(fEMCutNames[1]);
		fCounter[Elec].fill(fEECutNames[1]);


		// Compute event-by-event weights:
		gHLTSF = 1.;
		gBtagSF  = 1.;
		gBtagSF1 = 1.;
		gBtagSF2 = 1.;

		fDoCounting = false;
		// disable for now if( S->datamc!=0 ) {
		// disable for now 	int ind1(-1), ind2(-1);
		// disable for now 	fDoCounting = false;
		// disable for now 	if(isSSLLMuEvent(       ind1, ind2)) gHLTSF = diMuonHLTSF2012();
		// disable for now 	else if(isSSLLElMuEvent(ind1, ind2)) gHLTSF = muEleHLTSF2012();
		// disable for now 	else if(isSSLLElEvent(  ind1, ind2)) gHLTSF = diEleHLTSF2012();
		// disable for now 	fDoCounting = true;
		// disable for now }

		// disable for now // marcs try on the b-tag scale factors:
		// disable for now if( S->datamc!=0 ) {
		// disable for now 	bool isFastsim = false;
		// disable for now 	int ind1(-1), ind2(-1);
		// disable for now 	fDoCounting = false;
		// disable for now 	float myweight1(1.), myweight2(1.);
		// disable for now 	std::vector< int > btags = getNBTagsMedIndices();
		// disable for now 	int nbtags = btags.size();
		// disable for now 	if(isSSLLMuEvent(ind1, ind2) || isSSLLElMuEvent(ind1, ind2) || isSSLLElEvent(ind1, ind2) ) {
		// disable for now 		if (fVerbose > 2 && nbtags > 1) {
		// disable for now 			cout << "There are " << nbtags << " btags in this event!" << endl;
		// disable for now 			cout << "jet pt of 0th btagged-jet: " << getJetPt(btags[0]) << endl;
		// disable for now 			cout << "      index and tagger value: " << btags[0] << "  " << JetCSVBTag[btags[0]] << endl;
		// disable for now 			cout << "jet pt of 1th btagged-jet: " << getJetPt(btags[1]) << endl;
		// disable for now 			cout << "      index and tagger value: " << btags[1] << "  " << JetCSVBTag[btags[1]] << endl;
		// disable for now 			cout << endl;
		// disable for now 		}
		// disable for now 		if (nbtags == 2) {
		// disable for now 			myweight1 = btagEventWeight (2, getJetPt(btags[0]), getJetPt(btags[1]), isFastsim);
		// disable for now 		}
		// disable for now 		if (nbtags == 3) {
		// disable for now 			myweight1 = btagEventWeight (3, getJetPt(btags[0]), getJetPt(btags[1]), getJetPt(btags[2]), isFastsim);
		// disable for now 			myweight2 = btagEventWeight3(3, getJetPt(btags[0]), getJetPt(btags[1]), getJetPt(btags[2]), isFastsim);
		// disable for now 		}
		// disable for now 		if (nbtags  > 3) {
		// disable for now 			myweight1 = btagEventWeight (4, getJetPt(btags[0]), getJetPt(btags[1]), getJetPt(btags[2]), getJetPt(btags[3]), isFastsim );
		// disable for now 			myweight2 = btagEventWeight3(4, getJetPt(btags[0]), getJetPt(btags[1]), getJetPt(btags[2]), getJetPt(btags[3]), isFastsim );
		// disable for now 		}
		// disable for now 	}
		// disable for now 	gBtagSF1 = myweight1;
		// disable for now 	gBtagSF2 = myweight2;
		// disable for now 	fDoCounting = true;
		// disable for now }
		
		
		gEventWeight  = gHLTSF * gBtagSF;
		if( S->datamc!=4 ) gEventWeight *= PUWeight; // no pu weight for mc with no pu, that really doesn't exist anymore in 2012, so it's fine


		fillKinPlots(S,gRegion[gBaseRegion]);
		if (gDoWZValidation) fillKinPlots(S,gRegion["WZEnriched"]);
		for(regIt = gRegions.begin(); regIt != gRegions.end(); regIt++) {
			if ( (*regIt)->sname == "TTbarWSel") fDoCounting = true;
			fillYields(S, gRegion[(*regIt)->sname]);
			if ( (*regIt)->sname == "TTbarWSel") fDoCounting = false;
		}

		
		// disable for now // reset the event weight to what it should be without a region cut applied
		// disable for now gEventWeight  = gHLTSF * 1.;
		// disable for now if( S->datamc!=4 ) { 			//no pu weight for mc with no pu
		// disable for now 	gEventWeight *= PUWeight; // PU weight is always 1. for data, so no problem here
		// disable for now }

		// fillYields(S, gRegion["TTbarWPresel"]);

		// fDoCounting = true;
		// fillYields(S, gRegion["TTbarWSel"]);
		// fDoCounting = false;

		fillSigEventTree(S, 0);
		fillDiffYields(S);
		fillRatioPlots(S);
		fillMuIsoPlots(S);
		fillElIsoPlots(S);
		fillElIdPlots(S);
		
		//		fillSyncCounters(S);
		/////////////////////////////////////////////
		// Systematic studies
		if(!gDoSystStudies) continue;
		fChain->GetEntry(jentry); // reset tree vars
		resetBTags();

		// Jet pts scaled down
		smearJetPts(S, 1);
 		// fillYields(S, gRegion["TTbarWSelJU"]);
		fillSigEventTree(S, gSystematics["JetUp"]);

		// Jet pts scaled down
		fChain->GetEntry(jentry); // reset tree vars
		resetBTags();
		smearJetPts(S, 2);
 		// fillYields(S, gRegion["TTbarWSelJD"]);
		fillSigEventTree(S, gSystematics["JetDown"]);

		// Jet pts smeared
		fChain->GetEntry(jentry); // reset tree vars
		resetBTags();
		smearJetPts(S, 3);
 		// fillYields(S, gRegion["TTbarWSelJS"]);
		fillSigEventTree(S, gSystematics["JetSmear"]);

		// Btags scaled up
		fChain->GetEntry(jentry); // reset tree vars
		scaleBTags(S, 1);
 		// fillYields(S, gRegion["TTbarWSelBU"]);
		fillSigEventTree(S, gSystematics["BUp"]);

		// Btags scaled down
		fChain->GetEntry(jentry); // reset tree vars
		scaleBTags(S, 2);
 		// fillYields(S, gRegion["TTbarWSelBD"]);
		fillSigEventTree(S, gSystematics["BDown"]);

		// Lepton pts scaled up
		fChain->GetEntry(jentry); // reset tree vars
		resetBTags();
		scaleLeptons(S, 1);
 		// fillYields(S, gRegion["TTbarWSelLU"]);
		fillSigEventTree(S, gSystematics["LepUp"]);

		// Lepton pts scaled down
		fChain->GetEntry(jentry); // reset tree vars
		resetBTags();
		scaleLeptons(S, 2);
 		// fillYields(S, gRegion["TTbarWSelLD"]);
		fillSigEventTree(S, gSystematics["LepDown"]);
	}

	// Stuff to execute for each sample AFTER looping on the events
	fillCutFlowHistos(S);
	// marc printCutFlow(Muon);
	// marc printCutFlow(Elec);
	// marc printCutFlow(ElMu);
	
	writeHistos(S, pFile);
	writeSigGraphs(S, Muon, pFile);
	writeSigGraphs(S, Elec, pFile);
	writeSigGraphs(S, ElMu, pFile);

	writeSigEvTree(pFile);

	deleteHistos(S);
	S->cleanUp();

	pFile->Write();
	pFile->Close();

	if(S->datamc == 0) fOUTSTREAM.close();
	fDoCounting = false;
}
//____________________________________________________________________________
bool SSDLDumper::IsInJSON(){
        if (!fGoodRunList->IsInitialized()) {
	        cout << "[ERROR]: GoodRunList is not initialized " << endl;
		return false;
	}
	
	// Run processing
	if( fCurRun != Run ) { // new run
 	        fCurRun = Run;
		skipRun = false;
		if ( fGoodRunList->CheckRun(Run) == false ) skipRun = true;
	}
	
	// Check if new lumi is in JSON file
	if( fCurLumi != LumiSec ) { // new lumisection
	        fCurLumi = LumiSec;
		skipLumi = false;
		if ( fGoodRunList->CheckRunLumi(Run,LumiSec) == false ) skipLumi = true;
	}
	if(skipRun || skipLumi) return false;

	return true;
}

//____________________________________________________________________________
void SSDLDumper::fillYields(Sample *S, int reg){
	///////////////////////////////////////////////////
	// Set custom event selections here:
	setRegionCuts(reg);
	
	///////////////////////////////////////////////////
	// SS YIELDS
	// MuMu Channel
	resetHypLeptons();
	fDoCounting = false;
	
	if(reg == gRegion[gBaseRegion])   fDoCounting  = true;

	if(reg == gRegion["WZEnriched"])  gInvertZVeto = true;
	else                              gInvertZVeto = false;
	
	fCurrentChannel = Muon;
	int mu1(-1), mu2(-1);
	if(mumuSignalTrigger()){ // Trigger selection
		if(fDoCounting) fCounter[Muon].fill(fMMCutNames[2]);
		if(isSSLLMuEvent(mu1, mu2)){ // Same-sign loose-loose di muon event
			if(  isTightMuon(mu1) &&  isTightMuon(mu2) ){ // Tight-tight
				if(fDoCounting) fCounter[Muon].fill(fMMCutNames[15]); // ... first muon passes tight cut
				if(fDoCounting) fCounter[Muon].fill(fMMCutNames[16]); // ... second muon passes tight cut
				if(fDoCounting) fCounter[Muon].fill(fMMCutNames[17]); // ... both muons pass tight cut
				S->region[reg][HighPt].mm.nt20_pt ->Fill(MuPt [mu1], MuPt [mu2], gEventWeight);
				S->region[reg][HighPt].mm.nt20_eta->Fill(fabs(MuEta[mu1]), fabs(MuEta[mu2]), gEventWeight);
				if(S->datamc == 0 && reg == gRegion[gBaseRegion]){
				  fOUTSTREAM << Form("%12s: MuMu - run %6.0d / ls %5.0d / ev %11.0d - HT(#J/#bJ) %6.2f(%1d/%1d) MET %6.2f MT2 %6.2f Pt1 %6.2f Pt2 %6.2f Charge %2d", S->sname.Data(), Run, LumiSec, Event, getHT(), getNJets(), getNBTags(), getMET(), getMT2(mu1,mu2,Muon), MuPt[mu1], MuPt[mu2], MuCharge[mu1]) << endl ;
				}
				if(S->datamc > 0 ){
					S->region[reg][HighPt].mm.nt11_origin->Fill(muIndexToBin(mu1)-0.5, muIndexToBin(mu2)-0.5, gEventWeight);
					if(isPromptMuon(mu1) && isPromptMuon(mu2)) S->region[reg][HighPt].mm.nt2pp_pt->Fill(MuPt[mu1], MuPt[mu2], gEventWeight);
					if(isPromptMuon(mu1) &&   isFakeMuon(mu2)) S->region[reg][HighPt].mm.nt2pf_pt->Fill(MuPt[mu1], MuPt[mu2], gEventWeight);
					if(  isFakeMuon(mu1) && isPromptMuon(mu2)) S->region[reg][HighPt].mm.nt2fp_pt->Fill(MuPt[mu1], MuPt[mu2], gEventWeight);
					if(  isFakeMuon(mu1) &&   isFakeMuon(mu2)) S->region[reg][HighPt].mm.nt2ff_pt->Fill(MuPt[mu1], MuPt[mu2], gEventWeight);
				}
			}
			if(  isTightMuon(mu1) && !isTightMuon(mu2) ){ // Tight-loose
				if(fDoCounting) fCounter[Muon].fill(fMMCutNames[15]); // ... first muon passes tight cut
				S->region[reg][HighPt].mm.nt10_pt ->Fill(MuPt [mu1], MuPt [mu2], gEventWeight);
				S->region[reg][HighPt].mm.nt10_eta->Fill(fabs(MuEta[mu1]), fabs(MuEta[mu2]), gEventWeight);
				if(S->datamc > 0) S->region[reg][HighPt].mm.nt10_origin->Fill(muIndexToBin(mu1)-0.5, muIndexToBin(mu2)-0.5, gEventWeight);
			}
			if( !isTightMuon(mu1) &&  isTightMuon(mu2) ){ // Loose-tight
				if(fDoCounting) fCounter[Muon].fill(fMMCutNames[16]); // ... second muon passes tight cut
				S->region[reg][HighPt].mm.nt10_pt ->Fill(MuPt [mu2], MuPt [mu1], gEventWeight); // tight one always in x axis; fill same again
				S->region[reg][HighPt].mm.nt10_eta->Fill(fabs(MuEta[mu2]), fabs(MuEta[mu1]), gEventWeight);
				if(S->datamc > 0) S->region[reg][HighPt].mm.nt10_origin->Fill(muIndexToBin(mu2)-0.5, muIndexToBin(mu1)-0.5, gEventWeight);
			}
			if( !isTightMuon(mu1) && !isTightMuon(mu2) ){ // Loose-loose
				S->region[reg][HighPt].mm.nt00_pt ->Fill(MuPt [mu1], MuPt [mu2], gEventWeight);
				S->region[reg][HighPt].mm.nt00_eta->Fill(fabs(MuEta[mu1]), fabs(MuEta[mu2]), gEventWeight);
				if(S->datamc > 0) S->region[reg][HighPt].mm.nt00_origin->Fill(muIndexToBin(mu1)-0.5, muIndexToBin(mu2)-0.5, gEventWeight);
			}
			if(S->datamc > 0){
				if(isPromptMuon(mu1) && isPromptMuon(mu2)) S->region[reg][HighPt].mm.npp_pt->Fill(MuPt[mu1], MuPt[mu2], gEventWeight);
				if(isPromptMuon(mu1) &&   isFakeMuon(mu2)) S->region[reg][HighPt].mm.npf_pt->Fill(MuPt[mu1], MuPt[mu2], gEventWeight);
				if(  isFakeMuon(mu1) && isPromptMuon(mu2)) S->region[reg][HighPt].mm.nfp_pt->Fill(MuPt[mu1], MuPt[mu2], gEventWeight);
				if(  isFakeMuon(mu1) &&   isFakeMuon(mu2)) S->region[reg][HighPt].mm.nff_pt->Fill(MuPt[mu1], MuPt[mu2], gEventWeight);			
			}
		}
		resetHypLeptons();
	}
	
	// filling histos for ttbar only ratios. this is a fairly random place to do this, but ok...
	if(S->sname == "TTJets"){
		for (int i =0 ; i < NMus; ++i) {
			if( isLooseMuon(i) && !IsSignalMuon[i] ) S->region[reg][HighPt].mm.fnloose_ttbar->Fill(MuPt[i], fabs(MuEta[i]), gEventWeight);
			if( isTightMuon(i) && !IsSignalMuon[i] ) S->region[reg][HighPt].mm.fntight_ttbar->Fill(MuPt[i], fabs(MuEta[i]), gEventWeight);
			if( isLooseMuon(i) &&  IsSignalMuon[i] ) S->region[reg][HighPt].mm.pnloose_ttbar->Fill(MuPt[i], fabs(MuEta[i]), gEventWeight);
			if( isTightMuon(i) &&  IsSignalMuon[i] ) S->region[reg][HighPt].mm.pntight_ttbar->Fill(MuPt[i], fabs(MuEta[i]), gEventWeight);
		}
		for (int i =0 ; i < NEls; ++i) {
			if ( abs(ElGenID[i])  == 13 ) continue;
			if ( abs(ElGenID[i]) == 0  ) continue;
			if ( abs(ElGenMID[i]) == 13 ) continue;
			if( isLooseElectron(i) && !IsSignalElectron[i] ) S->region[reg][HighPt].ee.fnloose_ttbar->Fill(ElPt[i], fabs(ElEta[i]), gEventWeight);
			if( isTightElectron(i) && !IsSignalElectron[i] ) S->region[reg][HighPt].ee.fntight_ttbar->Fill(ElPt[i], fabs(ElEta[i]), gEventWeight);
			if( isLooseElectron(i) &&  IsSignalElectron[i] ) S->region[reg][HighPt].ee.pnloose_ttbar->Fill(ElPt[i], fabs(ElEta[i]), gEventWeight);
			if( isTightElectron(i) &&  IsSignalElectron[i] ) S->region[reg][HighPt].ee.pntight_ttbar->Fill(ElPt[i], fabs(ElEta[i]), gEventWeight);
		}
	}
	// end filling the ttbar ratio histograms
	
	// filling histos for with gen ID of fake leptons
	if (isSigSupMuEvent() > -1 && getNBTagsMed() > 0) {
		for (int i = 0 ; i < NMus; ++i) {
			if (getClosestJetDPhi(i, Muon) < 2.0) continue;
//			if (getClosestJetDR(i, Muon) < 0.5) continue;
			if( isLooseMuon(i) && !IsSignalMuon[i] ) S->region[reg][HighPt].mm.fnloose_genID->Fill(fabs(MuGenID[i])/*, gEventWeight*/);
			if( isTightMuon(i) && !IsSignalMuon[i] ) S->region[reg][HighPt].mm.fntight_genID->Fill(fabs(MuGenID[i])/*, gEventWeight*/);
		}
	}
	int elec1(-1), elec2(-1);
	if (isSigSupElEvent() > -1) {
		for (int i = 0 ; i < NEls; ++i) {
			if ( abs(ElGenID[i])  == 13 ) continue;
			if ( abs(ElGenMID[i]) == 13 ) continue;
//			if (getClosestJetDR(i, Elec) < 0.5) continue;
			if( isLooseElectron(i) && !IsSignalElectron[i] ) S->region[reg][HighPt].ee.fnloose_genID->Fill(fabs(ElGenID[i])/*, gEventWeight*/);
			if( isTightElectron(i) && !IsSignalElectron[i] ) S->region[reg][HighPt].ee.fntight_genID->Fill(fabs(ElGenID[i])/*, gEventWeight*/);
		}
	}
	// end filling the ttbar ratio histograms

	int looseMuInd = isSigSupMuEvent();
	if(singleMuTrigger() && isSigSupMuEvent() > -1){
		if( isTightMuon(looseMuInd) ){
			S->region[reg][HighPt].mm.fntight->Fill(MuPt[looseMuInd], fabs(MuEta[looseMuInd]), gEventWeight);
			// marc S->region[reg][HighPt].mm.fntight->Fill(MuPt[0], fabs(MuEta[0]), singleMuPrescale() * gEventWeight);
			S->region[reg][HighPt].mm.fntight_nv->Fill(NVrtx,                   gEventWeight);
			if(S->datamc > 0) S->region[reg][HighPt].mm.sst_origin->Fill(muIndexToBin(0)-0.5, gEventWeight);
		}
		if( isLooseMuon(looseMuInd) ){
			S->region[reg][HighPt].mm.fratio_pt ->Fill(isTightMuon(looseMuInd), MuPt[looseMuInd]);
			S->region[reg][HighPt].mm.fratio_eta->Fill(isTightMuon(looseMuInd), fabs(MuEta[looseMuInd]));
			S->region[reg][HighPt].mm.fratio_nv ->Fill(isTightMuon(looseMuInd), NVrtx);
			// marc S->region[reg][HighPt].mm.fratio_pt ->FillWeighted(isTightMuon(0), singleMuPrescale() * gEventWeight, MuPt[0]);
			// marc S->region[reg][HighPt].mm.fratio_eta->FillWeighted(isTightMuon(0), singleMuPrescale() * gEventWeight, fabs(MuEta[0]));

			S->region[reg][HighPt].mm.fnloose->Fill(MuPt[looseMuInd], fabs(MuEta[looseMuInd]), gEventWeight);
			S->region[reg][HighPt].mm.fnloose_nv->Fill(NVrtx,                   gEventWeight);
			// marc S->region[reg][HighPt].mm.fnloose->Fill(MuPt[0], fabs(MuEta[0]), singleMuPrescale() * gEventWeight);
			if(S->datamc > 0) S->region[reg][HighPt].mm.ssl_origin->Fill(muIndexToBin(0)-0.5, gEventWeight);
		}
	}
	if(doubleMuTrigger() && isZMuMuEvent(mu1, mu2)){
		if( isTightMuon(mu2) ){
			S->region[reg][HighPt].mm.pntight->Fill(MuPt[mu2], fabs(MuEta[mu2]), gEventWeight);
			S->region[reg][HighPt].mm.pntight_nv->Fill(NVrtx,                    gEventWeight);
			if(S->datamc > 0) S->region[reg][HighPt].mm.zt_origin->Fill(muIndexToBin(mu2)-0.5, gEventWeight);
		}
		if( isLooseMuon(mu2) ){
			S->region[reg][HighPt].mm.pratio_pt ->Fill(isTightMuon(mu2), MuPt[mu2]);
			S->region[reg][HighPt].mm.pratio_eta->Fill(isTightMuon(mu2), fabs(MuEta[mu2]));
			S->region[reg][HighPt].mm.pratio_nv ->Fill(isTightMuon(mu2), NVrtx);
			// S->region[reg][HighPt].mm.pratio_pt ->FillWeighted(isTightMuon(mu2), gEventWeight, MuPt[mu2]);
			// S->region[reg][HighPt].mm.pratio_eta->FillWeighted(isTightMuon(mu2), gEventWeight, fabs(MuEta[mu2]));

			S->region[reg][HighPt].mm.pnloose->Fill(MuPt[mu2], fabs(MuEta[mu2]), gEventWeight);
			S->region[reg][HighPt].mm.pnloose_nv->Fill(NVrtx,                       gEventWeight);
			if(S->datamc > 0) S->region[reg][HighPt].mm.zl_origin->Fill(muIndexToBin(mu2)-0.5, gEventWeight);
		}
	}
	resetHypLeptons();


	// EE Channel
	fCurrentChannel = Elec;
	int el1(-1), el2(-1);
	if(elelSignalTrigger()){
		if(fDoCounting) fCounter[Elec].fill(fEECutNames[2]);
		if( isSSLLElEvent(el1, el2) ){
			if(  isTightElectron(el1) &&  isTightElectron(el2) ){ // Tight-tight
				if(fDoCounting) fCounter[Elec].fill(fEECutNames[15]); // " ... first electron passes tight cut
				if(fDoCounting) fCounter[Elec].fill(fEECutNames[16]); // " ... second electron passes tight cut
				if(fDoCounting) fCounter[Elec].fill(fEECutNames[17]); // " ... both electrons pass tight cut
				S->region[reg][HighPt].ee.nt20_pt ->Fill(ElPt [el1], ElPt [el2], gEventWeight);
				S->region[reg][HighPt].ee.nt20_eta->Fill(fabs(ElEta[el1]), fabs(ElEta[el2]), gEventWeight);
				if(S->datamc == 0 && reg == gRegion[gBaseRegion] && HighPt == HighPt){
					fOUTSTREAM << Form("%12s: ElEl - run %6.0d / ls %5.0d / ev %11.0d - HT(#J/#bJ) %6.2f(%1d/%1d) MET %6.2f MT2 %6.2f Pt1 %6.2f Pt2 %6.2f Charge %2d", S->sname.Data(), Run, LumiSec, Event, getHT(), getNJets(), getNBTags(), getMET(), getMT2(el1,el2,Elec), ElPt[el1], ElPt[el2], ElCharge[el1]) << endl ;
				}
				if(S->datamc > 0 ){
					S->region[reg][HighPt].ee.nt11_origin->Fill(elIndexToBin(el1)-0.5, elIndexToBin(el2)-0.5, gEventWeight);
					if(isPromptElectron(el1) && isPromptElectron(el2)){
						S->region[reg][HighPt].ee.nt2pp_pt->Fill(ElPt[el1], ElPt[el2], gEventWeight);
						if(!isChargeMatchedElectron(el1) || !isChargeMatchedElectron(el2)){
							S->region[reg][HighPt].ee.nt2pp_cm_pt->Fill(ElPt[el1], ElPt[el2], gEventWeight);							
						}
					}
					if(isPromptElectron(el1) &&   isFakeElectron(el2)) S->region[reg][HighPt].ee.nt2pf_pt->Fill(ElPt[el1], ElPt[el2], gEventWeight);
					if(  isFakeElectron(el1) && isPromptElectron(el2)) S->region[reg][HighPt].ee.nt2fp_pt->Fill(ElPt[el1], ElPt[el2], gEventWeight);
					if(  isFakeElectron(el1) &&   isFakeElectron(el2)) S->region[reg][HighPt].ee.nt2ff_pt->Fill(ElPt[el1], ElPt[el2], gEventWeight);
				}
			}
			if(  isTightElectron(el1) && !isTightElectron(el2) ){ // Tight-loose
				if(fDoCounting) fCounter[Elec].fill(fEECutNames[15]);
				S->region[reg][HighPt].ee.nt10_pt ->Fill(ElPt [el1], ElPt [el2], gEventWeight);
				S->region[reg][HighPt].ee.nt10_eta->Fill(fabs(ElEta[el1]), fabs(ElEta[el2]), gEventWeight);
				if(S->datamc > 0) S->region[reg][HighPt].ee.nt10_origin->Fill(elIndexToBin(el1)-0.5, elIndexToBin(el2)-0.5, gEventWeight);
			}
			if( !isTightElectron(el1) &&  isTightElectron(el2) ){ // Loose-tight
				if(fDoCounting) fCounter[Elec].fill(fEECutNames[16]);
				S->region[reg][HighPt].ee.nt10_pt ->Fill(ElPt [el2], ElPt [el1], gEventWeight); // tight one always in x axis; fill same again
				S->region[reg][HighPt].ee.nt10_eta->Fill(fabs(ElEta[el2]), fabs(ElEta[el2]), gEventWeight);
				if(S->datamc > 0) S->region[reg][HighPt].ee.nt10_origin->Fill(elIndexToBin(el2)-0.5, elIndexToBin(el1)-0.5, gEventWeight);
			}
			if( !isTightElectron(el1) && !isTightElectron(el2) ){ // Loose-loose
				S->region[reg][HighPt].ee.nt00_pt ->Fill(ElPt [el1], ElPt [el2], gEventWeight);
				S->region[reg][HighPt].ee.nt00_eta->Fill(fabs(ElEta[el1]), fabs(ElEta[el2]), gEventWeight);
				if(S->datamc > 0) S->region[reg][HighPt].ee.nt00_origin->Fill(elIndexToBin(el1)-0.5, elIndexToBin(el2)-0.5, gEventWeight);
			}
			if(S->datamc > 0 ){
				if(isPromptElectron(el1) && isPromptElectron(el2)){
					S->region[reg][HighPt].ee.npp_pt->Fill(ElPt[el1], ElPt[el2], gEventWeight);
					if(!isChargeMatchedElectron(el1) || !isChargeMatchedElectron(el2)){
						S->region[reg][HighPt].ee.npp_cm_pt->Fill(ElPt[el1], ElPt[el2], gEventWeight);							
					}
				}
				if(isPromptElectron(el1) &&   isFakeElectron(el2)) S->region[reg][HighPt].ee.npf_pt->Fill(ElPt[el1], ElPt[el2], gEventWeight);
				if(  isFakeElectron(el1) && isPromptElectron(el2)) S->region[reg][HighPt].ee.nfp_pt->Fill(ElPt[el1], ElPt[el2], gEventWeight);
				if(  isFakeElectron(el1) &&   isFakeElectron(el2)) S->region[reg][HighPt].ee.nff_pt->Fill(ElPt[el1], ElPt[el2], gEventWeight);
			}
		}
		resetHypLeptons();
	}
	int looseElInd = isSigSupElEvent();
	if(singleElTrigger() && isSigSupElEvent() > -1){
		if( isTightElectron(looseElInd) ){
			S->region[reg][HighPt].ee.fntight->Fill(ElPt[looseElInd], fabs(ElEta[looseElInd]), gEventWeight);
			S->region[reg][HighPt].ee.fntight_nv->Fill(NVrtx,                   gEventWeight);
			// marc S->region[reg][HighPt].ee.fntight->Fill(ElPt[0], fabs(ElEta[0]), singleElPrescale() * gEventWeight);
			if(S->datamc > 0) S->region[reg][HighPt].ee.sst_origin->Fill(elIndexToBin(0)-0.5, gEventWeight);
		}
		if( isLooseElectron(looseElInd) ){
			S->region[reg][HighPt].ee.fratio_pt ->Fill(isTightElectron(looseElInd), ElPt[looseElInd]);
			S->region[reg][HighPt].ee.fratio_eta->Fill(isTightElectron(looseElInd), fabs(ElEta[looseElInd]));
			S->region[reg][HighPt].ee.fratio_nv ->Fill(isTightElectron(looseElInd), NVrtx);
			// marc S->region[reg][HighPt].ee.fratio_pt ->FillWeighted(isTightElectron(0), singleElPrescale() * gEventWeight, ElPt[0]);
			// marc S->region[reg][HighPt].ee.fratio_eta->FillWeighted(isTightElectron(0), singleElPrescale() * gEventWeight, fabs(ElEta[0]));

			S->region[reg][HighPt].ee.fnloose->Fill(ElPt[looseElInd], fabs(ElEta[looseElInd]), gEventWeight);
			S->region[reg][HighPt].ee.fnloose_nv->Fill(NVrtx,                   gEventWeight);
			// marc S->region[reg][HighPt].ee.fnloose->Fill(ElPt[0], fabs(ElEta[0]), singleElPrescale() * gEventWeight);
			if(S->datamc > 0) S->region[reg][HighPt].ee.ssl_origin->Fill(elIndexToBin(0)-0.5, gEventWeight);
		}
	}
	int elind;
	if(doubleElTrigger() && isZElElEvent(el1, el2)){
		if( isTightElectron(el2) ){
			S->region[reg][HighPt].ee.pntight->Fill(ElPt[el2], fabs(ElEta[el2]), gEventWeight);
			S->region[reg][HighPt].ee.pntight_nv->Fill(NVrtx,                       gEventWeight);
			if(S->datamc > 0) S->region[reg][HighPt].ee.zt_origin->Fill(elIndexToBin(el2)-0.5, gEventWeight);
		}
		if( isLooseElectron(el2) ){
			S->region[reg][HighPt].ee.pratio_pt ->Fill(isTightElectron(el2), ElPt[el2]);
			S->region[reg][HighPt].ee.pratio_eta->Fill(isTightElectron(el2), fabs(ElEta[el2]));
			S->region[reg][HighPt].ee.pratio_nv ->Fill(isTightElectron(el2), NVrtx);
			// S->region[reg][HighPt].ee.pratio_pt ->FillWeighted(isTightElectron(el2), gEventWeight, ElPt[el2]);
			// S->region[reg][HighPt].ee.pratio_eta->FillWeighted(isTightElectron(el2), gEventWeight, fabs(ElEta[el2]));

			S->region[reg][HighPt].ee.pnloose->Fill(ElPt[el2], fabs(ElEta[el2]), gEventWeight);
			S->region[reg][HighPt].ee.pnloose_nv->Fill(NVrtx,                       gEventWeight);
			if(S->datamc > 0) S->region[reg][HighPt].ee.zl_origin->Fill(elIndexToBin(el2)-0.5, gEventWeight);
		}
	}
	if (doubleElTrigger() && isZElElChMisIdEvent(el1, el2)){
	  if ( ElCharge[el1] != ElCharge[el2]) {  //OS pair
	    S->region[reg][HighPt].ee.ospairs->Fill(fabs(ElEta[el1]), fabs(ElEta[el2]), gEventWeight);
	  }
	  if ( ElCharge[el1] == ElCharge[el2]) {  //SS pair
	    S->region[reg][HighPt].ee.sspairs->Fill(fabs(ElEta[el1]), fabs(ElEta[el2]), gEventWeight);
	  }
	}
	
	resetHypLeptons();
	
	// EMu Channel
	fCurrentChannel = ElMu;
	int mu(-1), el(-1);
	if(elmuSignalTrigger()){
		if(fDoCounting) fCounter[ElMu].fill(fEMCutNames[2]);
		if( isSSLLElMuEvent(mu, el) ){
			if(  isTightElectron(el) &&  isTightMuon(mu) ){ // Tight-tight
				if(fDoCounting) fCounter[ElMu].fill(fEMCutNames[15]);
				if(fDoCounting) fCounter[ElMu].fill(fEMCutNames[16]);
				if(fDoCounting) fCounter[ElMu].fill(fEMCutNames[17]);
				S->region[reg][HighPt].em.nt20_pt ->Fill(MuPt [mu], ElPt [el], gEventWeight);
				S->region[reg][HighPt].em.nt20_eta->Fill(fabs(MuEta[mu]), fabs(ElEta[el]), gEventWeight);
				if(S->datamc == 0 && reg == gRegion[gBaseRegion] && HighPt == HighPt){
					fOUTSTREAM << Form("%12s: ElMu - run %6.0d / ls %5.0d / ev %11.0d - HT(#J/#bJ) %6.2f(%1d/%1d) MET %6.2f MT2 %6.2f Pt1 %6.2f Pt2 %6.2f Charge %2d", S->sname.Data(), Run, LumiSec, Event, getHT(), getNJets(), getNBTags(), getMET(), getMT2(mu,el,ElMu), MuPt[mu], ElPt[el], ElCharge[el]) << endl;
				}
				
				if(S->datamc > 0){
					S->region[reg][HighPt].em.nt11_origin->Fill(muIndexToBin(mu)-0.5, elIndexToBin(el)-0.5, gEventWeight);
					if(isPromptMuon(mu) && isPromptElectron(el)) S->region[reg][HighPt].em.nt2pp_pt->Fill(MuPt[mu], ElPt[el], gEventWeight);
					if(isPromptMuon(mu) &&   isFakeElectron(el)) S->region[reg][HighPt].em.nt2pf_pt->Fill(MuPt[mu], ElPt[el], gEventWeight);
					if(  isFakeMuon(mu) && isPromptElectron(el)) S->region[reg][HighPt].em.nt2fp_pt->Fill(MuPt[mu], ElPt[el], gEventWeight);
					if(  isFakeMuon(mu) &&   isFakeElectron(el)) S->region[reg][HighPt].em.nt2ff_pt->Fill(MuPt[mu], ElPt[el], gEventWeight);
				}
			}
			if( !isTightElectron(el) &&  isTightMuon(mu) ){ // Tight-loose
				if(fDoCounting) fCounter[ElMu].fill(fEMCutNames[15]);
				S->region[reg][HighPt].em.nt10_pt ->Fill(MuPt [mu], ElPt [el], gEventWeight);
				S->region[reg][HighPt].em.nt10_eta->Fill(fabs(MuEta[mu]), fabs(ElEta[el]), gEventWeight);
				if(S->datamc > 0) S->region[reg][HighPt].em.nt10_origin->Fill(muIndexToBin(mu)-0.5, elIndexToBin(el)-0.5, gEventWeight);
			}
			if(  isTightElectron(el) && !isTightMuon(mu) ){ // Loose-tight
				if(fDoCounting) fCounter[ElMu].fill(fEMCutNames[16]);
				S->region[reg][HighPt].em.nt01_pt ->Fill(MuPt [mu], ElPt [el], gEventWeight); // muon always in x axis for e/mu
				S->region[reg][HighPt].em.nt01_eta->Fill(fabs(MuEta[mu]), fabs(ElEta[el]), gEventWeight);
				if(S->datamc > 0) S->region[reg][HighPt].em.nt01_origin->Fill(muIndexToBin(mu)-0.5, elIndexToBin(el)-0.5, gEventWeight);
			}
			if( !isTightElectron(el) && !isTightMuon(mu) ){ // Loose-loose
				S->region[reg][HighPt].em.nt00_pt ->Fill(MuPt [mu], ElPt [el], gEventWeight);
				S->region[reg][HighPt].em.nt00_eta->Fill(fabs(MuEta[mu]), fabs(ElEta[el]), gEventWeight);
				if(S->datamc > 0) S->region[reg][HighPt].em.nt00_origin->Fill(muIndexToBin(mu)-0.5, elIndexToBin(el)-0.5, gEventWeight);
			}
			if(S->datamc > 0){
				if(isPromptMuon(mu) && isPromptElectron(el)) S->region[reg][HighPt].em.npp_pt->Fill(MuPt[mu], ElPt[el], gEventWeight);
				if(isPromptMuon(mu) &&   isFakeElectron(el)) S->region[reg][HighPt].em.npf_pt->Fill(MuPt[mu], ElPt[el], gEventWeight);
				if(  isFakeMuon(mu) && isPromptElectron(el)) S->region[reg][HighPt].em.nfp_pt->Fill(MuPt[mu], ElPt[el], gEventWeight);
				if(  isFakeMuon(mu) &&   isFakeElectron(el)) S->region[reg][HighPt].em.nff_pt->Fill(MuPt[mu], ElPt[el], gEventWeight);
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
				if( isBarrelElectron(el1) &&  isBarrelElectron(el2)) S->region[reg][HighPt].ee.nt20_OS_BB_pt->Fill(ElPt[el1], ElPt[el2], gEventWeight);
				if(!isBarrelElectron(el1) && !isBarrelElectron(el2)) S->region[reg][HighPt].ee.nt20_OS_EE_pt->Fill(ElPt[el1], ElPt[el2], gEventWeight);
				if( isBarrelElectron(el1) && !isBarrelElectron(el2)) S->region[reg][HighPt].ee.nt20_OS_EB_pt->Fill(ElPt[el1], ElPt[el2], gEventWeight);
				if(!isBarrelElectron(el1) &&  isBarrelElectron(el2)) S->region[reg][HighPt].ee.nt20_OS_EB_pt->Fill(ElPt[el1], ElPt[el2], gEventWeight);
			}
			if(  isTightElectron(el1) && !isTightElectron(el2) ){ // Tight-loose
				if( isBarrelElectron(el1) &&  isBarrelElectron(el2)) S->region[reg][HighPt].ee.nt10_OS_BB_pt->Fill(ElPt[el1], ElPt[el2], gEventWeight);
				if(!isBarrelElectron(el1) && !isBarrelElectron(el2)) S->region[reg][HighPt].ee.nt10_OS_EE_pt->Fill(ElPt[el1], ElPt[el2], gEventWeight);
				if( isBarrelElectron(el1) && !isBarrelElectron(el2)) S->region[reg][HighPt].ee.nt10_OS_EB_pt->Fill(ElPt[el1], ElPt[el2], gEventWeight);
				if(!isBarrelElectron(el1) &&  isBarrelElectron(el2)) S->region[reg][HighPt].ee.nt10_OS_EB_pt->Fill(ElPt[el1], ElPt[el2], gEventWeight);
			}
			if( !isTightElectron(el1) &&  isTightElectron(el2) ){ // Loose-tight
				if( isBarrelElectron(el1) &&  isBarrelElectron(el2)) S->region[reg][HighPt].ee.nt01_OS_BB_pt->Fill(ElPt[el1], ElPt[el2], gEventWeight);
				if(!isBarrelElectron(el1) && !isBarrelElectron(el2)) S->region[reg][HighPt].ee.nt01_OS_EE_pt->Fill(ElPt[el1], ElPt[el2], gEventWeight);
				if( isBarrelElectron(el1) && !isBarrelElectron(el2)) S->region[reg][HighPt].ee.nt01_OS_EB_pt->Fill(ElPt[el1], ElPt[el2], gEventWeight);
				if(!isBarrelElectron(el1) &&  isBarrelElectron(el2)) S->region[reg][HighPt].ee.nt01_OS_EB_pt->Fill(ElPt[el1], ElPt[el2], gEventWeight);
			}
			if( !isTightElectron(el1) && !isTightElectron(el2) ){ // Loose-loose
				if( isBarrelElectron(el1) &&  isBarrelElectron(el2)) S->region[reg][HighPt].ee.nt00_OS_BB_pt->Fill(ElPt[el1], ElPt[el2], gEventWeight);
				if(!isBarrelElectron(el1) && !isBarrelElectron(el2)) S->region[reg][HighPt].ee.nt00_OS_EE_pt->Fill(ElPt[el1], ElPt[el2], gEventWeight);
				if( isBarrelElectron(el1) && !isBarrelElectron(el2)) S->region[reg][HighPt].ee.nt00_OS_EB_pt->Fill(ElPt[el1], ElPt[el2], gEventWeight);
				if(!isBarrelElectron(el1) &&  isBarrelElectron(el2)) S->region[reg][HighPt].ee.nt00_OS_EB_pt->Fill(ElPt[el1], ElPt[el2], gEventWeight);
			}
		}
		resetHypLeptons();
	}

	// EMu Channel
	fCurrentChannel = ElMu;
	if(elmuSignalTrigger()){
		if( isSSLLElMuEvent(mu, el) ){ // this selects now OS events with the exact same cuts
			if(  isTightElectron(el) &&  isTightMuon(mu) ){ // Tight-tight
				if( isBarrelElectron(el)) S->region[reg][HighPt].em.nt20_OS_BB_pt->Fill(MuPt[mu], ElPt[el], gEventWeight);
				if(!isBarrelElectron(el)) S->region[reg][HighPt].em.nt20_OS_EE_pt->Fill(MuPt[mu], ElPt[el], gEventWeight);
			}
			if( !isTightElectron(el) &&  isTightMuon(mu) ){ // Tight-loose
				if( isBarrelElectron(el)) S->region[reg][HighPt].em.nt10_OS_BB_pt->Fill(MuPt[mu], ElPt[el], gEventWeight);
				if(!isBarrelElectron(el)) S->region[reg][HighPt].em.nt10_OS_EE_pt->Fill(MuPt[mu], ElPt[el], gEventWeight);
			}
			if(  isTightElectron(el) && !isTightMuon(mu) ){ // Loose-tight
				if( isBarrelElectron(el)) S->region[reg][HighPt].em.nt01_OS_BB_pt->Fill(MuPt[mu], ElPt[el], gEventWeight);
				if(!isBarrelElectron(el)) S->region[reg][HighPt].em.nt01_OS_EE_pt->Fill(MuPt[mu], ElPt[el], gEventWeight);
			}
			if( !isTightElectron(el) && !isTightMuon(mu) ){ // Loose-loose
				if( isBarrelElectron(el)) S->region[reg][HighPt].em.nt00_OS_BB_pt->Fill(MuPt[mu], ElPt[el], gEventWeight);
				if(!isBarrelElectron(el)) S->region[reg][HighPt].em.nt00_OS_EE_pt->Fill(MuPt[mu], ElPt[el], gEventWeight);
			}
		}
		resetHypLeptons();
	}
	fChargeSwitch = 0;
	fDoCounting = false;
	resetHypLeptons();
}
void SSDLDumper::fillDiffYields(Sample *S){
	// {  0 ,   1  ,    2   ,   3  ,   4  ,   5  ,    6    ,   7   ,      8     ,      9      }
	// {"HT", "MET", "NJets", "MT2", "PT1", "PT2", "NBJets", "MET3", "NBJetsMed", "NBJetsMed2"}
	///////////////////////////////////////////////////
	setRegionCuts(gRegion[gBaseRegion]);
	
	////////////////////////////////////////////////////////////////////////////////
	// SS YIELDS
	fDoCounting = false;
	resetHypLeptons();

	////////////////////////////////////////////////////////////////////////////////
	// MuMu Channel
	fCurrentChannel = Muon;
	int mu1(-1), mu2(-1);
	if(mumuSignalTrigger()){ // Trigger selection
		////////////////////////////
		// njets binning, ttbarW preselection
		setRegionCuts(gRegion[gBaseRegion]);
		fC_minNjets = 0;
		if(isSSLLMuEvent(mu1, mu2)) fillDiffVar(S, mu1, mu2, getNJets()+0.5, 2, Muon);
		resetHypLeptons();

		////////////////////////////
		// Other binnings, ttbarW preselection
		setRegionCuts(gRegion[gBaseRegion]);
		if(isSSLLMuEvent(mu1, mu2)){
			fillDiffVar(S, mu1, mu2, getHT(),                0, Muon);
			fillDiffVar(S, mu1, mu2, getMET(),               1, Muon);
			fillDiffVar(S, mu1, mu2, getMT2(mu1, mu2, Muon), 3, Muon);
			fillDiffVar(S, mu1, mu2, MuPt[mu1],              4, Muon);
			fillDiffVar(S, mu1, mu2, MuPt[mu2],              5, Muon);
			fillDiffVar(S, mu1, mu2, getNBTags()+0.5,        6, Muon);
			fillDiffVar(S, mu1, mu2, getNBTagsMed()+0.5,     8, Muon);
		}
		resetHypLeptons();

		////////////////////////////
		// MET BINNING, HT > 0
		setRegionCuts(gRegion[gBaseRegion]);
		fC_minMet = 30.;
// cout << __LINE__ << "min HT: " << fC_minHT << " minMET: " << fC_minMet << " minNjets: " << fC_minNjets << " minNbjets: " << fC_minNbjets << endl;
		fC_minHT =  0.;
		fC_minNjets = 0;
		if(isSSLLMuEvent(mu1, mu2)) fillDiffVar(S, mu1, mu2, getMET(), 7, Muon);
		resetHypLeptons();
		
		////////////////////////////
		// NBtag binning, ttbarW Selection
		setRegionCuts(gRegion[gBaseRegion]);
		fC_minNbjets = 0;
		fC_minNbjmed = 0;
		if(isSSLLMuEvent(mu1, mu2)) fillDiffVar(S, mu1, mu2, getNBTagsMed()+0.5, 9, Muon);
		resetHypLeptons();
	}
	setRegionCuts(gRegion[gBaseRegion]);

	////////////////////////////////////////////////////////////////////////////////
	// EE Channel
	fCurrentChannel = Elec;
	int el1(-1), el2(-1);
	if(elelSignalTrigger()){
		////////////////////////////
		// njets binning, ttbarW preselection
		setRegionCuts(gRegion[gBaseRegion]);
		fC_minNjets = 0;
		if(isSSLLElEvent(el1, el2)) fillDiffVar(S, el1, el2, getNJets()+0.5, 2, Elec);
		resetHypLeptons();

		////////////////////////////
		// Other binnings, ttbarW preselection
		setRegionCuts(gRegion[gBaseRegion]);
		if(isSSLLElEvent(el1, el2)){
			fillDiffVar(S, el1, el2, getHT(),                0, Elec);
			fillDiffVar(S, el1, el2, getMET(),               1, Elec);
			fillDiffVar(S, el1, el2, getMT2(el1, el2, Elec), 3, Elec);
			fillDiffVar(S, el1, el2, ElPt[el1],              4, Elec);
			fillDiffVar(S, el1, el2, ElPt[el2],              5, Elec);
			fillDiffVar(S, el1, el2, getNBTags()+0.5,        6, Elec);
			fillDiffVar(S, el1, el2, getNBTagsMed()+0.5,     8, Elec);
		}
		resetHypLeptons();

		////////////////////////////
		// MET BINNING, HT > 0
		setRegionCuts(gRegion[gBaseRegion]);
		fC_minMet = 30.;
		fC_minHT =  0.;
		fC_minNjets = 0.;
		if(isSSLLElEvent(el1, el2)) fillDiffVar(S, el1, el2, getMET(), 7, Elec);
		resetHypLeptons();
		
		////////////////////////////
		// NBtag binning, ttbarW Selection
		setRegionCuts(gRegion[gBaseRegion]);
		fC_minNbjets = 0;
		fC_minNbjmed = 0;
		if(isSSLLElEvent(el1, el2)) fillDiffVar(S, el1, el2, getNBTagsMed()+0.5, 9, Elec);
		resetHypLeptons();
	}
	setRegionCuts(gRegion[gBaseRegion]);

	////////////////////////////////////////////////////////////////////////////////
	// EMu Channel
	fCurrentChannel = ElMu;
	int mu(-1), el(-1);
	if(elmuSignalTrigger()){
		////////////////////////////
		// njets binning, ttbarW preselection
		setRegionCuts(gRegion[gBaseRegion]);
		fC_minNjets = 0;
		if(isSSLLElMuEvent(mu, el)) fillDiffVar(S, mu, el, getNJets()+0.5, 2, ElMu);
		resetHypLeptons();

		////////////////////////////
		// Other binnings, ttbarW preselection
		setRegionCuts(gRegion[gBaseRegion]);
		if(isSSLLElMuEvent(mu, el)){
			fillDiffVar(S, mu, el, getHT(),              0, ElMu);
			fillDiffVar(S, mu, el, getMET(),             1, ElMu);
			fillDiffVar(S, mu, el, getMT2(mu, el, ElMu), 3, ElMu);
			fillDiffVar(S, mu, el, getNBTags()+0.5,      6, ElMu);
			fillDiffVar(S, mu, el, getNBTagsMed()+0.5,   8, ElMu);
			float ptmax = MuPt[mu];
			float ptmin = ElPt[el];
			if(ptmin > ptmax){
				ptmin = MuPt[mu];
				ptmax = ElPt[el];
			}
			fillDiffVar(S, mu, el, ptmax,          4, ElMu);
			fillDiffVar(S, mu, el, ptmin,          5, ElMu);
		}
		resetHypLeptons();

		////////////////////////////
		// MET BINNING, HT > 0
		setRegionCuts(gRegion[gBaseRegion]);
		fC_minMet = 30.;
		fC_minHT =  0.;
		fC_minNjets = 0;
		if(isSSLLElMuEvent(mu, el)) fillDiffVar(S, mu, el, getMET(), 7, ElMu);
		resetHypLeptons();
		
		////////////////////////////
		// NBtag binning, ttbarW Selection
		setRegionCuts(gRegion[gBaseRegion]);
		fC_minNbjets = 0;
		fC_minNbjmed = 0;
		if(isSSLLElMuEvent(mu, el)) fillDiffVar(S, mu, el, getNBTagsMed()+0.5, 9, ElMu);
		resetHypLeptons();
	}
	setRegionCuts(gRegion[gBaseRegion]);

	///////////////////////////////////////////////////
	// OS YIELDS
	fChargeSwitch = 1;

	////////////////////////////////////////////////////////////////////////////////
	// EE Channel
	fCurrentChannel = Elec;
	if(elelSignalTrigger()){

		////////////////////////////
		// njets binning, ttbarW preselection
		setRegionCuts(gRegion[gBaseRegion]);
		fC_minNjets = 0;
		if(isSSLLElEvent(el1, el2)) fillDiffVarOS(S, el1, el2, getNJets()+0.5, 2, Elec);
		resetHypLeptons();

		////////////////////////////
		// Other binnings, ttbarW preselection
		setRegionCuts(gRegion[gBaseRegion]);
		if(isSSLLElEvent(el1, el2)){
			fillDiffVarOS(S, el1, el2, getHT(),                0, Elec);
			fillDiffVarOS(S, el1, el2, getMET(),               1, Elec);
			fillDiffVarOS(S, el1, el2, getMT2(el1, el2, Elec), 3, Elec);
			fillDiffVarOS(S, el1, el2, ElPt[el1],              4, Elec);
			fillDiffVarOS(S, el1, el2, ElPt[el2],              5, Elec);
			fillDiffVarOS(S, el1, el2, getNBTags()+0.5,        6, Elec);
			fillDiffVarOS(S, el1, el2, getNBTagsMed()+0.5,     8, Elec);
		}
		resetHypLeptons();

		////////////////////////////
		// MET BINNING, HT > 0
		setRegionCuts(gRegion[gBaseRegion]);
		fC_minMet = 30.;
		fC_minHT =  0.;
		fC_minNjets = 0;
		if(isSSLLElEvent(el1, el2)) fillDiffVarOS(S, el1, el2, getMET(), 7, Elec);
		resetHypLeptons();
		
		////////////////////////////
		// NBtag binning, ttbarW Selection
		setRegionCuts(gRegion[gBaseRegion]);
		fC_minNbjets = 0;
		fC_minNbjmed = 0;
		if(isSSLLElEvent(el1, el2)) fillDiffVarOS(S, el1, el2, getNBTagsMed()+0.5, 9, Elec);
		resetHypLeptons();
	}
	setRegionCuts(gRegion[gBaseRegion]);
	
	////////////////////////////////////////////////////////////////////////////////
	// EMu Channel
	fCurrentChannel = ElMu;
	if(elmuSignalTrigger()){
		////////////////////////////
		// njets binning, ttbarW preselection
		setRegionCuts(gRegion[gBaseRegion]);
		fC_minNjets = 0;
		if(isSSLLElMuEvent(mu, el)) fillDiffVarOS(S, mu, el, getNJets()+0.5, 2, ElMu);
		resetHypLeptons();

		////////////////////////////
		// Other binnings, ttbarW preselection
		setRegionCuts(gRegion[gBaseRegion]);
		if(isSSLLElMuEvent(mu, el)){
			fillDiffVarOS(S, mu, el, getHT(),                        0, ElMu);
			fillDiffVarOS(S, mu, el, getMET(),                       1, ElMu);
			fillDiffVarOS(S, mu, el, getMT2(mu, el, ElMu),           3, ElMu);
			fillDiffVarOS(S, mu, el, TMath::Max(MuPt[mu], ElPt[el]), 4, ElMu);
			fillDiffVarOS(S, mu, el, TMath::Min(MuPt[mu], ElPt[el]), 5, ElMu);
			fillDiffVarOS(S, mu, el, getNBTags()+0.5,                6, ElMu);
			fillDiffVarOS(S, mu, el, getNBTagsMed()+0.5,             8, ElMu);
		}
		resetHypLeptons();

		////////////////////////////
		// MET BINNING, HT > 0
		setRegionCuts(gRegion[gBaseRegion]);
		fC_minMet = 30.;
		fC_minHT =  0.;
		fC_minNjets = 0;
		if(isSSLLElMuEvent(mu, el)) fillDiffVarOS(S, mu, el, getMET(), 7, ElMu);
		resetHypLeptons();
		
		////////////////////////////
		// NBtag binning, ttbarW Selection
		setRegionCuts(gRegion[gBaseRegion]);
		fC_minNbjets = 0;
		fC_minNbjmed = 0;
		if(isSSLLElMuEvent(mu, el)) fillDiffVarOS(S, mu, el, getNBTagsMed()+0.5, 9, ElMu);
		resetHypLeptons();
	}

	fChargeSwitch = 0;
	resetHypLeptons();
	setRegionCuts(gRegion[gBaseRegion]);
}
void SSDLDumper::fillDiffVar(  Sample *S, int lep1, int lep2, float val, int bin, gChannel chan){
	if(bin > gNDiffVars-1){
		cerr << "SSDLDumper::fillDiffVar ==> Trying to fill diff var out of bounds! Check number of vars!" << endl;
		exit(-1);
	}
	int tlcat = -1;
	if(chan == Muon){
		if(  isTightMuon(lep1) &&  isTightMuon(lep2) )         tlcat = 0; 
		if(  isTightMuon(lep1) && !isTightMuon(lep2) )         tlcat = 1;
		if( !isTightMuon(lep1) &&  isTightMuon(lep2) )         tlcat = 2;
		if( !isTightMuon(lep1) && !isTightMuon(lep2) )         tlcat = 3;
	}
	if(chan == Elec){
		if(  isTightElectron(lep1) &&  isTightElectron(lep2) ) tlcat = 0; 
		if(  isTightElectron(lep1) && !isTightElectron(lep2) ) tlcat = 1;
		if( !isTightElectron(lep1) &&  isTightElectron(lep2) ) tlcat = 2;
		if( !isTightElectron(lep1) && !isTightElectron(lep2) ) tlcat = 3;
	}
	if(chan == ElMu){
		if(  isTightMuon(lep1) &&  isTightElectron(lep2) )     tlcat = 0; 
		if(  isTightMuon(lep1) && !isTightElectron(lep2) )     tlcat = 1;
		if( !isTightMuon(lep1) &&  isTightElectron(lep2) )     tlcat = 2;
		if( !isTightMuon(lep1) && !isTightElectron(lep2) )     tlcat = 3;
	}

	if(tlcat == 0) fillWithoutOF(S->diffyields[chan].hnt11[bin], val, gEventWeight);
	if(tlcat == 1) fillWithoutOF(S->diffyields[chan].hnt10[bin], val, gEventWeight);
	if(tlcat == 2) fillWithoutOF(S->diffyields[chan].hnt01[bin], val, gEventWeight);
	if(tlcat == 3) fillWithoutOF(S->diffyields[chan].hnt00[bin], val, gEventWeight);
}
void SSDLDumper::fillDiffVarOS(Sample *S, int lep1, int lep2, float val, int bin, gChannel chan){
	if(bin > gNDiffVars-1){
		cerr << "SSDLDumper::fillDiffVar ==> Trying to fill diff var out of bounds! Check number of vars!" << endl;
		exit(-1);
	}
	int ebcat = -1;
	if(chan == Muon) return;
	if(chan == Elec){
		if( !isTightElectron(lep1) || !isTightElectron(lep2) ) return; // tight-tight

		if(  isBarrelElectron(lep1) &&  isBarrelElectron(lep2) ) ebcat = 0; 
		if( !isBarrelElectron(lep1) && !isBarrelElectron(lep2) ) ebcat = 1;
		if(  isBarrelElectron(lep1) && !isBarrelElectron(lep2) ) ebcat = 2;
		if( !isBarrelElectron(lep1) &&  isBarrelElectron(lep2) ) ebcat = 3;
	}


	if(chan == ElMu){
		if( !isTightMuon(lep1) || !isTightElectron(lep2) ) return; // tight-tight

		if(  isBarrelElectron(lep2) ) ebcat = 0; 
		if( !isBarrelElectron(lep2) ) ebcat = 1;
	}

	if(ebcat == 0) fillWithoutOF(S->diffyields[chan].hnt2_os_BB[bin], val, gEventWeight);
	if(ebcat == 1) fillWithoutOF(S->diffyields[chan].hnt2_os_EE[bin], val, gEventWeight);
	if(ebcat == 2) fillWithoutOF(S->diffyields[chan].hnt2_os_EB[bin], val, gEventWeight);
	if(ebcat == 3) fillWithoutOF(S->diffyields[chan].hnt2_os_EB[bin], val, gEventWeight);
}

void SSDLDumper::fillRatioPlots(Sample *S){
	resetHypLeptons();
	fDoCounting = false;
	fChargeSwitch = 0;

	// Reset event selections to baseline:
	setRegionCuts(gRegion[gBaseRegion]);

	FRatioPlots *RP0 = &S->ratioplots[0];
	FRatioPlots *RP1 = &S->ratioplots[1];

	fCurrentChannel = Muon;
	if(singleMuTrigger()){
		int looseMuInd = isSigSupMuEvent();
		if(isSigSupMuEvent() > -1){
			if( isTightMuon(looseMuInd) ){
				RP0->ntight[0]->Fill(getNJets(),               gEventWeight);
				RP0->ntight[1]->Fill(getHT(),                  gEventWeight);
				RP0->ntight[2]->Fill(getMaxJPt(),              gEventWeight);
				RP0->ntight[3]->Fill(NVrtx,                    gEventWeight);
				RP0->ntight[4]->Fill(getClosestJetPt(0, Muon), gEventWeight);
				RP0->ntight[5]->Fill(getAwayJetPt(0, Muon),    gEventWeight);
				RP0->ntight[6]->Fill(getNBTags(),              gEventWeight);
			}
			if( isLooseMuon(looseMuInd) ){
				RP0->nloose[0]->Fill(getNJets(),               gEventWeight);
				RP0->nloose[1]->Fill(getHT(),                  gEventWeight);
				RP0->nloose[2]->Fill(getMaxJPt(),              gEventWeight);
				RP0->nloose[3]->Fill(NVrtx,                    gEventWeight);
				RP0->nloose[4]->Fill(getClosestJetPt(0, Muon), gEventWeight);
				RP0->nloose[5]->Fill(getAwayJetPt(0, Muon),    gEventWeight);
				RP0->nloose[6]->Fill(getNBTags(),              gEventWeight);
			}
		}
		fC_maxMet_Control = 1000.;
		looseMuInd = isSigSupMuEvent();
		if(isSigSupMuEvent() > -1){
			if( isTightMuon(looseMuInd) ){
				RP0->ntight[7]->Fill(getMET(),                    gEventWeight);
			}
			if( isLooseMuon(looseMuInd) ){
				RP0->nloose[7]->Fill(getMET(),                    gEventWeight);
			}
		}		
		fC_maxMet_Control = 20.;
		fC_maxMt_Control = 1000.;
		looseMuInd = isSigSupMuEvent();
		if(isSigSupMuEvent() > -1){
			if( isTightMuon(looseMuInd) ){
				RP0->ntight[8]->Fill(MuMT[0],                  gEventWeight);
			}
			if( isLooseMuon(looseMuInd) ){
				RP0->nloose[8]->Fill(MuMT[0],                  gEventWeight);
			}
		}		
		fC_maxMt_Control = 20.;
	}
	resetHypLeptons();
	setRegionCuts(gRegion[gBaseRegion]);
	if(singleElTrigger()){
		int looseElInd = isSigSupElEvent();
		if(isSigSupElEvent() > -1){
			if( isTightElectron(looseElInd) ){
				RP1->ntight[0]->Fill(getNJets(),                   gEventWeight);
				RP1->ntight[1]->Fill(getHT(),                      gEventWeight);
				RP1->ntight[2]->Fill(getMaxJPt(),                  gEventWeight);
				RP1->ntight[3]->Fill(NVrtx,                        gEventWeight);
				RP1->ntight[4]->Fill(getClosestJetPt(0, Elec), gEventWeight);
				RP1->ntight[5]->Fill(getAwayJetPt(0, Elec),    gEventWeight);
				RP1->ntight[6]->Fill(getNBTags(),                  gEventWeight);
			}
			if( isLooseElectron(looseElInd) ){
				RP1->nloose[0]->Fill(getNJets(),                   gEventWeight);
				RP1->nloose[1]->Fill(getHT(),                      gEventWeight);
				RP1->nloose[2]->Fill(getMaxJPt(),                  gEventWeight);
				RP1->nloose[3]->Fill(NVrtx,                        gEventWeight);
				RP1->nloose[4]->Fill(getClosestJetPt(0, Elec), gEventWeight);
				RP1->nloose[5]->Fill(getAwayJetPt(0, Elec),    gEventWeight);
				RP1->nloose[6]->Fill(getNBTags(),                  gEventWeight);
			}
		}
		fC_maxMet_Control = 1000.;
		looseElInd = isSigSupElEvent();
		if(isSigSupElEvent() > -1){
			if( isTightElectron(looseElInd) ){
				RP1->ntight[7]->Fill(getMET(),                    gEventWeight);
			}
			if( isLooseElectron(looseElInd) ){
				RP1->nloose[7]->Fill(getMET(),                    gEventWeight);
			}
		}		
		fC_maxMet_Control = 20.;
		fC_maxMt_Control = 1000.;
		looseElInd = isSigSupElEvent();
		if(isSigSupElEvent() > -1){
			if( isTightElectron(looseElInd) ){
				RP1->ntight[8]->Fill(ElMT[0],                  gEventWeight);
			}
			if( isLooseElectron(looseElInd) ){
				RP1->nloose[8]->Fill(ElMT[0],                  gEventWeight);
			}
		}		
		fC_maxMt_Control = 20.;
	}
	setRegionCuts(gRegion[gBaseRegion]);
	resetHypLeptons();
}
void SSDLDumper::fillSigEventTree(Sample *S, int flag=0){
	resetSigEventTree();
	resetHypLeptons();
	fDoCounting = false;

	///////////////////////////////////////////////////
	// Set custom event selections here:
	setRegionCuts(gRegion[gBaseRegion]);
	fC_minMu1pt  = 10.; // lower pt cuts for sig tree
	fC_minEl1pt  = 10.;
	fC_minNjets  = 0;
	fC_minMet    = 0.;
	fC_app3rdVet  = 0;
	fC_chargeVeto = 0;
	// gApplyZVeto   = false;

	fSETree_SystFlag = flag;
	fSETree_PUWeight = PUWeight;
	fSETree_HLTSF    = gHLTSF;
	fSETree_BtagSF1  = gBtagSF1;
	fSETree_BtagSF2  = gBtagSF2;
	fSETree_SLumi    = S->getLumi();
	fSETree_SName    = S->sname.Data();
	fSETree_SType    = getSampleType(S);
	fSETree_Run      = Run;
	fSETree_LS       = LumiSec;
	fSETree_Event    = Event;
	fSETree_MET      = getMET();

	fSETree_NM = getNTightMuons();
	fSETree_NE = getNTightElectrons();

//	fChargeSwitch = 1;
	
	if( Event==gDEBUG_EVENTNUMBER_ && Run==gDEBUG_RUNNUMBER_ ) {

		std::cout << std::endl << std::endl << "------------------------------------------------------------" << std::endl;
		std::cout << "  Debug log for run: " << Run << "  LS: " << LumiSec << "  Event: " << Event << std::endl;

		std::cout << std::endl << "Here are the jets: " << std::endl;
		for(size_t i = 0; i < NJets; ++i) 
			std::cout << "Pt: " << getJetPt(i) << "  Eta: " << JetEta[i] << " TCHE: " << JetCSVBTag[i] << std::endl;

	}


	int ind1(-1), ind2(-1);
	
	////////////////////////////////////////////////////////////////////////////////////////////////////////
	// MUMU CHANNEL:  //////////////////////////////////////////////////////////////////////////////////////
	if(mumuSignalTrigger() && isSSLLMuEvent(ind1, ind2)){ // trigger && select loose mu/mu pair
		fSETree_M3      = getM3();
		fSETree_MT2     = getMT2(ind1, ind2, Muon);
		fSETree_Mll     = getMll(ind1, ind2, Muon);
		fSETree_HT      = getHT();
		fSETree_NJ      = getNJets();
		fSETree_NbJ     = getNBTags();
		fSETree_NbJmed  = getNBTagsMed();
		fSETree_Flavor  = 0;
		fSETree_Charge  = MuCharge[ind1];
		fSETree_pT1     = MuPt[ind1];
		fSETree_pT2     = MuPt[ind2];
		fSETree_eta1    = MuEta[ind1];
		fSETree_eta2    = MuEta[ind2];
		fSETree_ZVeto   = passesZVeto()?1:0;
		fSETree_3rdVeto = passes3rdLepVeto()?1:0;
		fSETree_ttZSel  = passesTTZSel()?1:0;
		fSETree_PFIso1  = MuPFIso[ind1];
		fSETree_PFIso2  = MuPFIso[ind2];
		if( isTightMuon(ind1)&& isTightMuon(ind2)) fSETree_TLCat = 0;
		if( isTightMuon(ind1)&&!isTightMuon(ind2)) fSETree_TLCat = 1;
		if(!isTightMuon(ind1)&& isTightMuon(ind2)) fSETree_TLCat = 2;
		if(!isTightMuon(ind1)&&!isTightMuon(ind2)) fSETree_TLCat = 3;
		fSigEv_Tree->Fill();

		if( Event==gDEBUG_EVENTNUMBER_ && Run==gDEBUG_RUNNUMBER_ ) {
			std::cout << " -> This is a mu-mu event." << std::endl;
			std::cout << "SystFlag :" << fSETree_SystFlag << std::endl;
			std::cout << "MT2      :" << fSETree_MT2    << std::endl;
			std::cout << "Mll      :" << fSETree_Mll    << std::endl;
			std::cout << "HT       :" << fSETree_HT     << std::endl;
			std::cout << "NJ       :" << fSETree_NJ     << std::endl;
			std::cout << "NbJ      :" << fSETree_NbJ    << std::endl;
			std::cout << "Flavor   :" << fSETree_Flavor << std::endl;
			std::cout << "Charge   :" << fSETree_Charge << std::endl;
			std::cout << "pT1      :" << fSETree_pT1    << std::endl;
			std::cout << "pT2      :" << fSETree_pT2    << std::endl;
			std::cout << "eta1     :" << fSETree_eta1   << std::endl;
			std::cout << "eta2     :" << fSETree_eta2   << std::endl;
			std::cout << "TLCat    :" << fSETree_TLCat  << std::endl;
		}

		resetHypLeptons();
	}

	////////////////////////////////////////////////////////////////////////////////////////////////////////
	// EMU CHANNEL:  ///////////////////////////////////////////////////////////////////////////////////////
	else if(elmuSignalTrigger() && isSSLLElMuEvent(ind1, ind2)){ // trigger && select loose e/mu pair
		fSETree_M3     = getM3();
		fSETree_MT2    = getMT2(ind1, ind2, ElMu);
		fSETree_Mll    = getMll(ind1, ind2, ElMu);
		fSETree_HT     = getHT();
		fSETree_NJ     = getNJets();
		fSETree_NbJ    = getNBTags();
		fSETree_NbJmed = getNBTagsMed();
		fSETree_Flavor = 1;
		fSETree_Charge = MuCharge[ind1];
		fSETree_pT1    = MuPt[ind1];
		fSETree_pT2    = ElPt[ind2];
		fSETree_eta1   = MuEta[ind1];
		fSETree_eta2   = ElEta[ind2];
		fSETree_ZVeto   = passesZVeto()?1:0;
		fSETree_3rdVeto = passes3rdLepVeto()?1:0;
		fSETree_ttZSel  = passesTTZSel()?1:0;
		fSETree_PFIso1  = MuPFIso[ind1];
		fSETree_PFIso2  = ElPFIso[ind2];
		// fSETree_MVAID1  = ElMVAIDTrig[ind1];
		fSETree_MVAID2  = ElMVAIDTrig[ind2];
		// fSETree_medWP1  = ElIsGoodElId_MediumWP[ind1];
		fSETree_medWP2  = ElIsGoodElId_MediumWP[ind2];
		if( isTightMuon(ind1)&& isTightElectron(ind2)) fSETree_TLCat = 0;
		if( isTightMuon(ind1)&&!isTightElectron(ind2)) fSETree_TLCat = 1;
		if(!isTightMuon(ind1)&& isTightElectron(ind2)) fSETree_TLCat = 2;
		if(!isTightMuon(ind1)&&!isTightElectron(ind2)) fSETree_TLCat = 3;
		fSigEv_Tree->Fill();

		if( Event==gDEBUG_EVENTNUMBER_ && Run==gDEBUG_RUNNUMBER_ ) {
			std::cout << " -> This is a mu-ele event." << std::endl;
			std::cout << "SystFlag :" << fSETree_SystFlag << std::endl;
			std::cout << "MT2      :" << fSETree_MT2    << std::endl;
			std::cout << "Mll      :" << fSETree_Mll    << std::endl;
			std::cout << "HT       :" << fSETree_HT     << std::endl;
			std::cout << "NJ       :" << fSETree_NJ     << std::endl;
			std::cout << "NbJ      :" << fSETree_NbJ    << std::endl;
			std::cout << "Flavor   :" << fSETree_Flavor << std::endl;
			std::cout << "Charge   :" << fSETree_Charge << std::endl;
			std::cout << "pT1      :" << fSETree_pT1    << std::endl;
			std::cout << "pT2      :" << fSETree_pT2    << std::endl;
			std::cout << "eta1     :" << fSETree_eta1   << std::endl;
			std::cout << "eta2     :" << fSETree_eta2   << std::endl;
			std::cout << "TLCat    :" << fSETree_TLCat  << std::endl;
		}

		resetHypLeptons();
	}

	////////////////////////////////////////////////////////////////////////////////////////////////////////
	// EE CHANNEL:  ////////////////////////////////////////////////////////////////////////////////////////
	else if(elelSignalTrigger() && isSSLLElEvent(ind1, ind2)){ // trigger && select loose e/e pair
		fSETree_M3     = getM3();
		fSETree_MT2    = getMT2(ind1, ind2, Elec);
		fSETree_Mll    = getMll(ind1, ind2, Elec);
		fSETree_HT     = getHT();
		fSETree_NJ     = getNJets();
		fSETree_NbJ    = getNBTags();
		fSETree_NbJmed = getNBTagsMed();
		fSETree_Flavor = 2;
		fSETree_Charge = ElCharge[ind1];
		fSETree_pT1    = ElPt[ind1];
		fSETree_pT2    = ElPt[ind2];
		fSETree_eta1   = ElEta[ind1];
		fSETree_eta2   = ElEta[ind2];

		fSETree_ZVeto   = passesZVeto()?1:0;
		fSETree_3rdVeto = passes3rdLepVeto()?1:0;
		fSETree_ttZSel  = passesTTZSel()?1:0;
		fSETree_PFIso1  = ElPFIso[ind1];
		fSETree_PFIso2  = ElPFIso[ind2];
		fSETree_MVAID1  = ElMVAIDTrig[ind1];
		fSETree_MVAID2  = ElMVAIDTrig[ind2];
		fSETree_medWP1  = ElIsGoodElId_MediumWP[ind1];
		fSETree_medWP2  = ElIsGoodElId_MediumWP[ind2];
		if( isTightElectron(ind1)&& isTightElectron(ind2)) fSETree_TLCat = 0;
		if( isTightElectron(ind1)&&!isTightElectron(ind2)) fSETree_TLCat = 1;
		if(!isTightElectron(ind1)&& isTightElectron(ind2)) fSETree_TLCat = 2;
		if(!isTightElectron(ind1)&&!isTightElectron(ind2)) fSETree_TLCat = 3;
		fSigEv_Tree->Fill();

		if( Event==gDEBUG_EVENTNUMBER_ && Run==gDEBUG_RUNNUMBER_ ) {
			std::cout << " -> This is an ele-ele event." << std::endl;
			std::cout << "SystFlag :" << fSETree_SystFlag << std::endl;
			std::cout << "MT2      :" << fSETree_MT2    << std::endl;
			std::cout << "Mll      :" << fSETree_Mll    << std::endl;
			std::cout << "HT       :" << fSETree_HT     << std::endl;
			std::cout << "NJ       :" << fSETree_NJ     << std::endl;
			std::cout << "NbJ      :" << fSETree_NbJ    << std::endl;
			std::cout << "Flavor   :" << fSETree_Flavor << std::endl;
			std::cout << "Charge   :" << fSETree_Charge << std::endl;
			std::cout << "pT1      :" << fSETree_pT1    << std::endl;
			std::cout << "pT2      :" << fSETree_pT2    << std::endl;
			std::cout << "eta1     :" << fSETree_eta1   << std::endl;
			std::cout << "eta2     :" << fSETree_eta2   << std::endl;
			std::cout << "TLCat    :" << fSETree_TLCat  << std::endl;
		}

		resetHypLeptons();
	}

	gApplyZVeto   = true;
	/// OS YIELDS only for data:
	if (S->datamc == 0) {
		fChargeSwitch = 1;

		////////////////////////////////////////////////////////////////////////////////////////////////////////
		// EM CHANNEL:  OS  ////////////////////////////////////////////////////////////////////////////////////////
		if (elmuSignalTrigger() && isSSLLElMuEvent(ind1, ind2)){ // trigger && select loose e/mu pair
			if ( isTightMuon(ind1) && isTightElectron(ind2)) {
				fSETree_M3      = getM3();
				fSETree_MT2     = getMT2(ind1, ind2, ElMu);
				fSETree_Mll     = getMll(ind1, ind2, ElMu);
				fSETree_HT      = getHT();
				fSETree_NJ      = getNJets();
				fSETree_NbJ     = getNBTags();
				fSETree_NbJmed  = getNBTagsMed();
				fSETree_Flavor  = 4;
				fSETree_Charge  = MuCharge[ind1];
				fSETree_pT1     = MuPt[ind1];
				fSETree_pT2     = ElPt[ind2];
				fSETree_eta1    = MuEta[ind1];
				fSETree_eta2    = ElEta[ind2];
				fSETree_ZVeto   = passesZVeto()?1:0;
				//fSETree_ZVeto   = passesZVetoNew(ind1, ind2, 1)?1:0;
				fSETree_3rdVeto = passes3rdLepVeto()?1:0;
				fSETree_ttZSel  = passesTTZSel()?1:0;
				fSETree_PFIso1  = MuPFIso[ind1];
				fSETree_PFIso2  = ElPFIso[ind2];
				if( isBarrelElectron(ind2)) fSETree_TLCat = 0; // TLCat == 0 if Barrel-Electron
				if(!isBarrelElectron(ind2)) fSETree_TLCat = 1; // TLCat == 1 if Endcap-Electron
				fSigEv_Tree->Fill();
			}
			resetHypLeptons();
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////////
		// EE CHANNEL:  OS  ////////////////////////////////////////////////////////////////////////////////////////
		if (elelSignalTrigger() && isSSLLElEvent(ind1, ind2)){ // trigger && select loose e/e pair
			if ( isTightElectron(ind1) && isTightElectron(ind2)) {
				fSETree_M3      = getM3();
				fSETree_MT2     = getMT2(ind1, ind2, Elec);
				fSETree_Mll     = getMll(ind1, ind2, Elec);
				fSETree_HT      = getHT();
				fSETree_NJ      = getNJets();
				fSETree_NbJ     = getNBTags();
				fSETree_NbJmed  = getNBTagsMed();
				fSETree_Flavor  = 5;
				fSETree_Charge  = ElCharge[ind1];
				fSETree_pT1     = ElPt[ind1];
				fSETree_pT2     = ElPt[ind2];
				fSETree_eta1    = ElEta[ind1];
				fSETree_eta2    = ElEta[ind2];
				fSETree_ZVeto   = passesZVeto()?1:0;
				//fSETree_ZVeto   = passesZVetoNew(ind1, ind2, 2)?1:0;
				fSETree_3rdVeto = passes3rdLepVeto()?1:0;
				fSETree_ttZSel  = passesTTZSel()?1:0;
				fSETree_PFIso1  = ElPFIso[ind1];
				fSETree_PFIso2  = ElPFIso[ind2];
				if( isBarrelElectron(ind1)&& isBarrelElectron(ind2)) fSETree_TLCat = 0; // TLCat == 0 if Barrel/Barrel
				if( isBarrelElectron(ind1)&&!isBarrelElectron(ind2)) fSETree_TLCat = 1; // TLCat == 1 if Barrel/Endcap
				if(!isBarrelElectron(ind1)&& isBarrelElectron(ind2)) fSETree_TLCat = 2; // TLCat == 2 if Endcap/Barrel
				if(!isBarrelElectron(ind1)&&!isBarrelElectron(ind2)) fSETree_TLCat = 3; // TLCat == 3 if Endcap/Endcap
				fSigEv_Tree->Fill();
			}
			resetHypLeptons();
		}
		fChargeSwitch = 0;
	}
//	fChargeSwitch = 0;
	setRegionCuts(gRegion[gBaseRegion]);
	return;
}
// OS version void SSDLDumper::fillSigEventTree(Sample *S, int flag=0){
// OS version 	resetSigEventTree();
// OS version 	resetHypLeptons();
// OS version 	fDoCounting = false;
// OS version 
// OS version 	///////////////////////////////////////////////////
// OS version 	// Set custom event selections here:
// OS version 	setRegionCuts();
// OS version 	fC_minMu1pt  = 10.; // lower pt cuts for sig tree
// OS version 	fC_minEl1pt  = 10.;
// OS version 
// OS version 	fSETree_SystFlag = flag;
// OS version 	fSETree_PUWeight = PUWeight;
// OS version 	fSETree_HLTSF    = gHLTSF;
// OS version 	fSETree_BtagSF1   = gBtagSF1;
// OS version 	fSETree_BtagSF2   = gBtagSF2;
// OS version 	fSETree_SLumi    = S->getLumi();
// OS version 	fSETree_SName    = S->sname.Data();
// OS version 	fSETree_SType    = getSampleType(S);
// OS version 	fSETree_Run      = Run;
// OS version 	fSETree_LS       = LumiSec;
// OS version 	fSETree_Event    = Event;
// OS version 	fSETree_MET      = getMET();
// OS version 
// OS version 	fSETree_NM = getNTightMuons();
// OS version 	fSETree_NE = getNTightElectrons();
// OS version 	
// OS version 	if( Event==gDEBUG_EVENTNUMBER_ && Run==gDEBUG_RUNNUMBER_ ) {
// OS version 
// OS version 		std::cout << std::endl << std::endl << "------------------------------------------------------------" << std::endl;
// OS version 		std::cout << "  Debug log for run: " << Run << "  LS: " << LumiSec << "  Event: " << Event << std::endl;
// OS version 
// OS version 		std::cout << std::endl << "Here are the jets: " << std::endl;
// OS version 		for(size_t i = 0; i < NJets; ++i) 
// OS version 			std::cout << "Pt: " << getJetPt(i) << "  Eta: " << JetEta[i] << " TCHE: " << JetCSVBTag[i] << std::endl;
// OS version 
// OS version 	}
// OS version 
// OS version 
// OS version 	int ind1(-1), ind2(-1);
// OS version 	
// OS version 	////////////////////////////////////////////////////////////////////////////////////////////////////////
// OS version 	// MUMU CHANNEL:  //////////////////////////////////////////////////////////////////////////////////////
// OS version 	if(mumuSignalTrigger() && isSSLLMuEvent(ind1, ind2)){ // trigger && select loose mu/mu pair
// OS version 		fSETree_MT2     = getMT2(ind1, ind2, Muon);
// OS version 		fSETree_Mll     = getMll(ind1, ind2, Muon);
// OS version 		fSETree_HT      = getHT();
// OS version 		fSETree_NJ      = getNJets();
// OS version 		fSETree_NbJ     = getNBTags();
// OS version 		fSETree_NbJmed  = getNBTagsMed();
// OS version 		fSETree_Flavor  = 0;
// OS version 		fSETree_Charge  = MuCharge[ind1];
// OS version 		fSETree_pT1     = MuPt[ind1];
// OS version 		fSETree_pT2     = MuPt[ind2];
// OS version 		fSETree_eta1    = MuEta[ind1];
// OS version 		fSETree_eta2    = MuEta[ind2];
// OS version 		fSETree_ZVeto   = passesZVeto()?1:0;
// OS version 		fSETree_3rdVeto = passes3rdLepVeto()?1:0;
// OS version 		fSETree_ttZSel  = passesTTZSel()?1:0;
// OS version 		fSETree_PFIso1  = MuPFIso[ind1];
// OS version 		fSETree_PFIso2  = MuPFIso[ind2];
// OS version 		if( isTightMuon(ind1)&& isTightMuon(ind2)) fSETree_TLCat = 0;
// OS version 		if( isTightMuon(ind1)&&!isTightMuon(ind2)) fSETree_TLCat = 1;
// OS version 		if(!isTightMuon(ind1)&& isTightMuon(ind2)) fSETree_TLCat = 2;
// OS version 		if(!isTightMuon(ind1)&&!isTightMuon(ind2)) fSETree_TLCat = 3;
// OS version 		fSigEv_Tree->Fill();
// OS version 
// OS version 		if( Event==gDEBUG_EVENTNUMBER_ && Run==gDEBUG_RUNNUMBER_ ) {
// OS version 			std::cout << " -> This is a mu-mu event." << std::endl;
// OS version 			std::cout << "SystFlag :" << fSETree_SystFlag << std::endl;
// OS version 			std::cout << "MT2      :" << fSETree_MT2    << std::endl;
// OS version 			std::cout << "Mll      :" << fSETree_Mll    << std::endl;
// OS version 			std::cout << "HT       :" << fSETree_HT     << std::endl;
// OS version 			std::cout << "NJ       :" << fSETree_NJ     << std::endl;
// OS version 			std::cout << "NbJ      :" << fSETree_NbJ    << std::endl;
// OS version 			std::cout << "Flavor   :" << fSETree_Flavor << std::endl;
// OS version 			std::cout << "Charge   :" << fSETree_Charge << std::endl;
// OS version 			std::cout << "pT1      :" << fSETree_pT1    << std::endl;
// OS version 			std::cout << "pT2      :" << fSETree_pT2    << std::endl;
// OS version 			std::cout << "eta1     :" << fSETree_eta1   << std::endl;
// OS version 			std::cout << "eta2     :" << fSETree_eta2   << std::endl;
// OS version 			std::cout << "TLCat    :" << fSETree_TLCat  << std::endl;
// OS version 		}
// OS version 
// OS version 		resetHypLeptons();
// OS version 	}
// OS version 
// OS version 	////////////////////////////////////////////////////////////////////////////////////////////////////////
// OS version 	// EMU CHANNEL:  ///////////////////////////////////////////////////////////////////////////////////////
// OS version 	else if(elmuSignalTrigger() && isSSLLElMuEvent(ind1, ind2)){ // trigger && select loose e/mu pair
// OS version 		fSETree_MT2    = getMT2(ind1, ind2, ElMu);
// OS version 		fSETree_Mll    = getMll(ind1, ind2, ElMu);
// OS version 		fSETree_HT     = getHT();
// OS version 		fSETree_NJ     = getNJets();
// OS version 		fSETree_NbJ    = getNBTags();
// OS version 		fSETree_NbJmed = getNBTagsMed();
// OS version 		fSETree_Flavor = 1;
// OS version 		fSETree_Charge = MuCharge[ind1];
// OS version 		fSETree_pT1    = MuPt[ind1];
// OS version 		fSETree_pT2    = ElPt[ind2];
// OS version 		fSETree_eta1   = MuEta[ind1];
// OS version 		fSETree_eta2   = ElEta[ind2];
// OS version 		fSETree_ZVeto   = passesZVeto()?1:0;
// OS version 		fSETree_3rdVeto = passes3rdLepVeto()?1:0;
// OS version 		fSETree_ttZSel  = passesTTZSel()?1:0;
// OS version 		fSETree_PFIso1  = MuPFIso[ind1];
// OS version 		fSETree_PFIso2  = ElPFIso[ind2];
// OS version 		if( isTightMuon(ind1)&& isTightElectron(ind2)) fSETree_TLCat = 0;
// OS version 		if( isTightMuon(ind1)&&!isTightElectron(ind2)) fSETree_TLCat = 1;
// OS version 		if(!isTightMuon(ind1)&& isTightElectron(ind2)) fSETree_TLCat = 2;
// OS version 		if(!isTightMuon(ind1)&&!isTightElectron(ind2)) fSETree_TLCat = 3;
// OS version 		fSigEv_Tree->Fill();
// OS version 
// OS version 		if( Event==gDEBUG_EVENTNUMBER_ && Run==gDEBUG_RUNNUMBER_ ) {
// OS version 			std::cout << " -> This is a mu-ele event." << std::endl;
// OS version 			std::cout << "SystFlag :" << fSETree_SystFlag << std::endl;
// OS version 			std::cout << "MT2      :" << fSETree_MT2    << std::endl;
// OS version 			std::cout << "Mll      :" << fSETree_Mll    << std::endl;
// OS version 			std::cout << "HT       :" << fSETree_HT     << std::endl;
// OS version 			std::cout << "NJ       :" << fSETree_NJ     << std::endl;
// OS version 			std::cout << "NbJ      :" << fSETree_NbJ    << std::endl;
// OS version 			std::cout << "Flavor   :" << fSETree_Flavor << std::endl;
// OS version 			std::cout << "Charge   :" << fSETree_Charge << std::endl;
// OS version 			std::cout << "pT1      :" << fSETree_pT1    << std::endl;
// OS version 			std::cout << "pT2      :" << fSETree_pT2    << std::endl;
// OS version 			std::cout << "eta1     :" << fSETree_eta1   << std::endl;
// OS version 			std::cout << "eta2     :" << fSETree_eta2   << std::endl;
// OS version 			std::cout << "TLCat    :" << fSETree_TLCat  << std::endl;
// OS version 		}
// OS version 
// OS version 		resetHypLeptons();
// OS version 	}
// OS version 
// OS version 	////////////////////////////////////////////////////////////////////////////////////////////////////////
// OS version 	// EE CHANNEL:  ////////////////////////////////////////////////////////////////////////////////////////
// OS version 	else if(elelSignalTrigger() && isSSLLElEvent(ind1, ind2)){ // trigger && select loose e/e pair
// OS version 		fSETree_MT2    = getMT2(ind1, ind2, Elec);
// OS version 		fSETree_Mll    = getMll(ind1, ind2, Elec);
// OS version 		fSETree_HT     = getHT();
// OS version 		fSETree_NJ     = getNJets();
// OS version 		fSETree_NbJ    = getNBTags();
// OS version 		fSETree_NbJmed = getNBTagsMed();
// OS version 		fSETree_Flavor = 2;
// OS version 		fSETree_Charge = ElCharge[ind1];
// OS version 		fSETree_pT1    = ElPt[ind1];
// OS version 		fSETree_pT2    = ElPt[ind2];
// OS version 		fSETree_eta1   = ElEta[ind1];
// OS version 		fSETree_eta2   = ElEta[ind2];
// OS version 		fSETree_ZVeto   = passesZVeto()?1:0;
// OS version 		fSETree_3rdVeto = passes3rdLepVeto()?1:0;
// OS version 		fSETree_ttZSel  = passesTTZSel()?1:0;
// OS version 		fSETree_PFIso1  = ElPFIso[ind1];
// OS version 		fSETree_PFIso2  = ElPFIso[ind2];
// OS version 		if( isTightElectron(ind1)&& isTightElectron(ind2)) fSETree_TLCat = 0;
// OS version 		if( isTightElectron(ind1)&&!isTightElectron(ind2)) fSETree_TLCat = 1;
// OS version 		if(!isTightElectron(ind1)&& isTightElectron(ind2)) fSETree_TLCat = 2;
// OS version 		if(!isTightElectron(ind1)&&!isTightElectron(ind2)) fSETree_TLCat = 3;
// OS version 		fSigEv_Tree->Fill();
// OS version 
// OS version 		if( Event==gDEBUG_EVENTNUMBER_ && Run==gDEBUG_RUNNUMBER_ ) {
// OS version 			std::cout << " -> This is an ele-ele event." << std::endl;
// OS version 			std::cout << "SystFlag :" << fSETree_SystFlag << std::endl;
// OS version 			std::cout << "MT2      :" << fSETree_MT2    << std::endl;
// OS version 			std::cout << "Mll      :" << fSETree_Mll    << std::endl;
// OS version 			std::cout << "HT       :" << fSETree_HT     << std::endl;
// OS version 			std::cout << "NJ       :" << fSETree_NJ     << std::endl;
// OS version 			std::cout << "NbJ      :" << fSETree_NbJ    << std::endl;
// OS version 			std::cout << "Flavor   :" << fSETree_Flavor << std::endl;
// OS version 			std::cout << "Charge   :" << fSETree_Charge << std::endl;
// OS version 			std::cout << "pT1      :" << fSETree_pT1    << std::endl;
// OS version 			std::cout << "pT2      :" << fSETree_pT2    << std::endl;
// OS version 			std::cout << "eta1     :" << fSETree_eta1   << std::endl;
// OS version 			std::cout << "eta2     :" << fSETree_eta2   << std::endl;
// OS version 			std::cout << "TLCat    :" << fSETree_TLCat  << std::endl;
// OS version 		}
// OS version 
// OS version 		resetHypLeptons();
// OS version 	}
// OS version 	/// OS YIELDS only for data:
// OS version 	if (S->datamc == 0) {
// OS version 		fChargeSwitch = 1;
// OS version 
// OS version 		////////////////////////////////////////////////////////////////////////////////////////////////////////
// OS version 		// MM CHANNEL:  OS  this is only for the OS people ////////////////////////////////////////////////////////////////////////////////////////
// OS version 		if (mumuSignalTrigger() && isSSLLMuEvent(ind1, ind2)){ // trigger && select loose mu/mu pair
// OS version 				fSETree_MT2     = getMT2(ind1, ind2, Muon);
// OS version 				fSETree_Mll     = getMll(ind1, ind2, Muon);
// OS version 				fSETree_HT      = getHT();
// OS version 				fSETree_NJ      = getNJets();
// OS version 				fSETree_NbJ     = getNBTags();
// OS version 				fSETree_NbJmed  = getNBTagsMed();
// OS version 				fSETree_Flavor  = 3;
// OS version 				fSETree_Charge  = MuCharge[ind1];
// OS version 				fSETree_pT1     = MuPt[ind1];
// OS version 				fSETree_pT2     = MuPt[ind2];
// OS version 				fSETree_eta1    = MuEta[ind1];
// OS version 				fSETree_eta2    = MuEta[ind2];
// OS version 				fSETree_ZVeto   = passesZVeto()?1:0;
// OS version 				fSETree_3rdVeto = passes3rdLepVeto()?1:0;
// OS version 				fSETree_ttZSel  = passesTTZSel()?1:0;
// OS version 				fSETree_PFIso1  = MuPFIso[ind1];
// OS version 				fSETree_PFIso2  = MuPFIso[ind2];
// OS version 				if( isTightMuon(ind1)&& isTightMuon(ind2)) fSETree_TLCat = 0;
// OS version 				if( isTightMuon(ind1)&&!isTightMuon(ind2)) fSETree_TLCat = 1;
// OS version 				if(!isTightMuon(ind1)&& isTightMuon(ind2)) fSETree_TLCat = 2;
// OS version 				if(!isTightMuon(ind1)&&!isTightMuon(ind2)) fSETree_TLCat = 3;
// OS version 				fSigEv_Tree->Fill();
// OS version 			resetHypLeptons();
// OS version 			}
// OS version 
// OS version 		////////////////////////////////////////////////////////////////////////////////////////////////////////
// OS version 		// EM CHANNEL:  OS  ////////////////////////////////////////////////////////////////////////////////////////
// OS version 		else if (elmuSignalTrigger() && isSSLLElMuEvent(ind1, ind2)){ // trigger && select loose e/mu pair
// OS version 				fSETree_MT2     = getMT2(ind1, ind2, ElMu);
// OS version 				fSETree_Mll     = getMll(ind1, ind2, ElMu);
// OS version 				fSETree_HT      = getHT();
// OS version 				fSETree_NJ      = getNJets();
// OS version 				fSETree_NbJ     = getNBTags();
// OS version 				fSETree_NbJmed  = getNBTagsMed();
// OS version 				fSETree_Flavor  = 4;
// OS version 				fSETree_Charge  = MuCharge[ind1];
// OS version 				fSETree_pT1     = MuPt[ind1];
// OS version 				fSETree_pT2     = ElPt[ind2];
// OS version 				fSETree_eta1    = MuEta[ind1];
// OS version 				fSETree_eta2    = ElEta[ind2];
// OS version 				fSETree_ZVeto   = passesZVeto()?1:0;
// OS version 				fSETree_3rdVeto = passes3rdLepVeto()?1:0;
// OS version 				fSETree_ttZSel  = passesTTZSel()?1:0;
// OS version 				fSETree_PFIso1  = MuPFIso[ind1];
// OS version 				fSETree_PFIso2  = ElPFIso[ind2];
// OS version 				if( isTightMuon(ind1)&& isTightElectron(ind2)) fSETree_TLCat = 0;
// OS version 				if( isTightMuon(ind1)&&!isTightElectron(ind2)) fSETree_TLCat = 1;
// OS version 				if(!isTightMuon(ind1)&& isTightElectron(ind2)) fSETree_TLCat = 2;
// OS version 				if(!isTightMuon(ind1)&&!isTightElectron(ind2)) fSETree_TLCat = 3;
// OS version 				fSigEv_Tree->Fill();
// OS version 			resetHypLeptons();
// OS version 			}
// OS version 
// OS version 		////////////////////////////////////////////////////////////////////////////////////////////////////////
// OS version 		// EE CHANNEL:  OS  ////////////////////////////////////////////////////////////////////////////////////////
// OS version 		else if (elelSignalTrigger() && isSSLLElEvent(ind1, ind2)){ // trigger && select loose e/e pair
// OS version 				fSETree_MT2     = getMT2(ind1, ind2, Elec);
// OS version 				fSETree_Mll     = getMll(ind1, ind2, Elec);
// OS version 				fSETree_HT      = getHT();
// OS version 				fSETree_NJ      = getNJets();
// OS version 				fSETree_NbJ     = getNBTags();
// OS version 				fSETree_NbJmed  = getNBTagsMed();
// OS version 				fSETree_Flavor  = 5;
// OS version 				fSETree_Charge  = ElCharge[ind1];
// OS version 				fSETree_pT1     = ElPt[ind1];
// OS version 				fSETree_pT2     = ElPt[ind2];
// OS version 				fSETree_eta1    = ElEta[ind1];
// OS version 				fSETree_eta2    = ElEta[ind2];
// OS version 				fSETree_ZVeto   = passesZVeto()?1:0;
// OS version 				fSETree_3rdVeto = passes3rdLepVeto()?1:0;
// OS version 				fSETree_ttZSel  = passesTTZSel()?1:0;
// OS version 				fSETree_PFIso1  = ElPFIso[ind1];
// OS version 				fSETree_PFIso2  = ElPFIso[ind2];
// OS version 				if( isTightElectron(ind1)&& isTightElectron(ind2)) fSETree_TLCat = 0;
// OS version 				if( isTightElectron(ind1)&&!isTightElectron(ind2)) fSETree_TLCat = 1;
// OS version 				if(!isTightElectron(ind1)&& isTightElectron(ind2)) fSETree_TLCat = 2;
// OS version 				if(!isTightElectron(ind1)&&!isTightElectron(ind2)) fSETree_TLCat = 3;
// OS version 				fSigEv_Tree->Fill();
// OS version 		resetHypLeptons();
// OS version 		}
// OS version 		fChargeSwitch = 0;
// OS version 	}
// OS version 
// OS version 	setRegionCuts();
// OS version 	return;
// OS version }
void SSDLDumper::fillKinPlots(Sample *S, int reg){
	resetHypLeptons();
	fDoCounting = false;
	KinPlots *KP0(0); 
	KinPlots *KP1(0); 
	KinPlots *KP2(0); 
	
	if (reg == gRegion["WZEnriched"] && gDoWZValidation){
	  KP0 = &S->kinplots_wz[0];
	  KP1 = &S->kinplots_wz[1];
	  KP2 = &S->kinplots_wz[2];
	}
	else {
	  KP0 = &S->kinplots[0][HighPt];
	  KP1 = &S->kinplots[1][HighPt];
	  KP2 = &S->kinplots[2][HighPt];	  
	}
	// cout << "[DEBUG] Initialized KP0 " << KP0 << endl;

	int ind1(-1), ind2(-1);
	int mu(-1), el(-1); // for e/mu channel, easier readability
	
	///////////////////////////////////////////////////
	// Set custom event selections here:
	setRegionCuts(reg); 
	
	////////////////////////////////////////////////////////////////////////////////////////////////////////
	// MUMU CHANNEL:  //////////////////////////////////////////////////////////////////////////////////////
	if(mumuSignalTrigger() && abs(isSSLLEvent(ind1, ind2)) == 1){ // trigger && select loose mu/mu pair
		if(MuPt[ind1] > fC_minMu1pt && MuPt[ind2] > fC_minMu2pt){ // pt cuts
			// Fill histos
			fillWithoutOF(KP0->hvar[0],  getHT(),                  gEventWeight);
			fillWithoutOF(KP0->hvar[1],  getMET(),                 gEventWeight);
			fillWithoutOF(KP0->hvar[2],  getNJets(),               gEventWeight);
			fillWithoutOF(KP0->hvar[3],  MuPt[ind1],               gEventWeight);
			fillWithoutOF(KP0->hvar[4],  MuPt[ind2],               gEventWeight);
			fillWithoutOF(KP0->hvar[5],  getMll(ind1, ind2, Muon), gEventWeight); // SF
			fillWithoutOF(KP0->hvar[6],  getMll(ind1, ind2, Muon), gEventWeight); // MM
			fillWithoutOF(KP0->hvar[9],  getMT2(ind1, ind2, Muon), gEventWeight);
			fillWithoutOF(KP0->hvar[10], getNBTags(),              gEventWeight);
			fillWithoutOF(KP0->hvar[11], getNBTagsMed(),           gEventWeight);

			if(isTightMuon(ind1) && isTightMuon(ind2)){ // tight-tight
				fillWithoutOF(KP1->hvar[0],  getHT(),                  gEventWeight);
				fillWithoutOF(KP1->hvar[1],  getMET(),                 gEventWeight);
				fillWithoutOF(KP1->hvar[2],  getNJets(),               gEventWeight);
				fillWithoutOF(KP1->hvar[3],  MuPt[ind1],               gEventWeight);
				fillWithoutOF(KP1->hvar[4],  MuPt[ind2],               gEventWeight);
				fillWithoutOF(KP1->hvar[5],  getMll(ind1, ind2, Muon), gEventWeight); // SF
				fillWithoutOF(KP1->hvar[6],  getMll(ind1, ind2, Muon), gEventWeight); // MM
				fillWithoutOF(KP1->hvar[9],  getMT2(ind1, ind2, Muon), gEventWeight);
				fillWithoutOF(KP1->hvar[10], getNBTags(),              gEventWeight);
				fillWithoutOF(KP1->hvar[11], getNBTagsMed(),           gEventWeight);

				if(passesMllEventVeto(ind1, ind2, 1, 8.)){
					// Want to fill njets with all but njets cut
					fillWithoutOF(KP2->hvar[2],  getNJets(),               gEventWeight);
					if(isSSLLMuEvent(ind1, ind2)){ // signal region
						fSigEv_HI_MM_HT .push_back(getHT());
						fSigEv_HI_MM_MET.push_back(getMET());
					
						fillWithoutOF(KP2->hvar[0],  getHT(),                  gEventWeight);
						fillWithoutOF(KP2->hvar[1],  getMET(),                 gEventWeight);
						fillWithoutOF(KP2->hvar[3],  MuPt[ind1],               gEventWeight);
						fillWithoutOF(KP2->hvar[4],  MuPt[ind2],               gEventWeight);
						fillWithoutOF(KP2->hvar[5],  getMll(ind1, ind2, Muon), gEventWeight); // SF
						fillWithoutOF(KP2->hvar[6],  getMll(ind1, ind2, Muon), gEventWeight); // MM
						fillWithoutOF(KP2->hvar[9],  getMT2(ind1, ind2, Muon), gEventWeight);
						fillWithoutOF(KP2->hvar[10], getNBTags(),              gEventWeight);
						fillWithoutOF(KP2->hvar[11], getNBTagsMed(),           gEventWeight);
					}
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
			fillWithoutOF(KP0->hvar[0],  getHT(),                  gEventWeight);
			fillWithoutOF(KP0->hvar[1],  getMET(),                 gEventWeight);
			fillWithoutOF(KP0->hvar[2],  getNJets(),               gEventWeight);
			fillWithoutOF(KP0->hvar[3],  ElPt[ind1],               gEventWeight);
			fillWithoutOF(KP0->hvar[4],  ElPt[ind2],               gEventWeight);
			fillWithoutOF(KP0->hvar[5],  getMll(ind1, ind2, Elec), gEventWeight); // SF
			fillWithoutOF(KP0->hvar[7],  getMll(ind1, ind2, Elec), gEventWeight); // EE
			fillWithoutOF(KP0->hvar[9],  getMT2(ind1, ind2, Elec), gEventWeight);
			fillWithoutOF(KP0->hvar[10], getNBTags(),              gEventWeight);
			fillWithoutOF(KP0->hvar[11], getNBTagsMed(),           gEventWeight);
			
			if(isTightElectron(ind1) && isTightElectron(ind2)){ // tight-tight
				fillWithoutOF(KP1->hvar[0],  getHT(),                  gEventWeight);
				fillWithoutOF(KP1->hvar[1],  getMET(),                 gEventWeight);
				fillWithoutOF(KP1->hvar[2],  getNJets(),               gEventWeight);
				fillWithoutOF(KP1->hvar[3],  ElPt[ind1],               gEventWeight);
				fillWithoutOF(KP1->hvar[4],  ElPt[ind2],               gEventWeight);
				fillWithoutOF(KP1->hvar[5],  getMll(ind1, ind2, Elec), gEventWeight); // SF
				fillWithoutOF(KP1->hvar[7],  getMll(ind1, ind2, Elec), gEventWeight); // EE
				fillWithoutOF(KP1->hvar[9],  getMT2(ind1, ind2, Elec), gEventWeight);
				fillWithoutOF(KP1->hvar[10], getNBTags(),              gEventWeight);
				fillWithoutOF(KP1->hvar[11], getNBTagsMed(),           gEventWeight);

				if(passesMllEventVeto(ind1, ind2, 2, 8.)){
					// Want to fill njets with all but njets cut
					fillWithoutOF(KP2->hvar[2],  getNJets(),               gEventWeight);

					if(isSSLLElEvent(ind1, ind2)){
						fSigEv_HI_EE_HT .push_back(getHT());
						fSigEv_HI_EE_MET.push_back(getMET());
										
						fillWithoutOF(KP2->hvar[0],  getHT(),                  gEventWeight);
						fillWithoutOF(KP2->hvar[1],  getMET(),                 gEventWeight);
						fillWithoutOF(KP2->hvar[3],  ElPt[ind1],               gEventWeight);
						fillWithoutOF(KP2->hvar[4],  ElPt[ind2],               gEventWeight);
						fillWithoutOF(KP2->hvar[5],  getMll(ind1, ind2, Elec), gEventWeight); // SF
						fillWithoutOF(KP2->hvar[7],  getMll(ind1, ind2, Elec), gEventWeight); // EE
						fillWithoutOF(KP2->hvar[9],  getMT2(ind1, ind2, Elec), gEventWeight);
						fillWithoutOF(KP2->hvar[10], getNBTags(),              gEventWeight);
						fillWithoutOF(KP2->hvar[11], getNBTagsMed(),           gEventWeight);
					}
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
			fillWithoutOF(KP0->hvar[0],  getHT(),              gEventWeight);
			fillWithoutOF(KP0->hvar[1],  getMET(),             gEventWeight);
			fillWithoutOF(KP0->hvar[2],  getNJets(),           gEventWeight);
			float ptmax = MuPt[mu];
			float ptmin = ElPt[el];
			if(ptmin > ptmax){
				ptmin = MuPt[mu];
				ptmax = ElPt[el];
			}
			fillWithoutOF(KP0->hvar[3],  ptmax,                gEventWeight);
			fillWithoutOF(KP0->hvar[4],  ptmin,                gEventWeight);
			fillWithoutOF(KP0->hvar[8],  getMll(mu, el, ElMu), gEventWeight); // SF
			fillWithoutOF(KP0->hvar[9],  getMT2(mu, el, ElMu), gEventWeight);
			fillWithoutOF(KP0->hvar[10], getNBTags(),          gEventWeight);			
			fillWithoutOF(KP0->hvar[11], getNBTagsMed(),       gEventWeight);

			if(isTightMuon(mu) && isTightElectron(el)){ // tight-tight
				fillWithoutOF(KP1->hvar[0],  getHT(),              gEventWeight);
				fillWithoutOF(KP1->hvar[1],  getMET(),             gEventWeight);
				fillWithoutOF(KP1->hvar[2],  getNJets(),           gEventWeight);
				fillWithoutOF(KP1->hvar[3],  ptmax,                gEventWeight);
				fillWithoutOF(KP1->hvar[4],  ptmin,                gEventWeight);
				fillWithoutOF(KP1->hvar[8],  getMll(mu, el, ElMu), gEventWeight); // SF
				fillWithoutOF(KP1->hvar[9],  getMT2(mu, el, ElMu), gEventWeight);
				fillWithoutOF(KP1->hvar[10], getNBTags(),          gEventWeight);			
				fillWithoutOF(KP1->hvar[11], getNBTagsMed(),       gEventWeight);

				fillWithoutOF(KP2->hvar[2],  getNJets(),           gEventWeight);
				if( isSSLLElMuEvent(mu, el) ){
					fSigEv_HI_EM_HT .push_back(getHT());
					fSigEv_HI_EM_MET.push_back(getMET());
					fillWithoutOF(KP2->hvar[0],  getHT(),              gEventWeight);
					fillWithoutOF(KP2->hvar[1],  getMET(),             gEventWeight);
					fillWithoutOF(KP2->hvar[3],  ptmax,                gEventWeight);
					fillWithoutOF(KP2->hvar[4],  ptmin,                gEventWeight);
					fillWithoutOF(KP2->hvar[8],  getMll(mu, el, ElMu), gEventWeight); // SF
					fillWithoutOF(KP2->hvar[9],  getMT2(mu, el, ElMu), gEventWeight);
					fillWithoutOF(KP2->hvar[10], getNBTags(),          gEventWeight);			
					fillWithoutOF(KP2->hvar[11], getNBTagsMed(),       gEventWeight);
				}
			}
		}
	}
	resetHypLeptons();
	setRegionCuts(gRegion[gBaseRegion]);
	return;
}
void SSDLDumper::fillMuIsoPlots(Sample *S){
	resetHypLeptons();
	fDoCounting = false;
	IsoPlots *IP0 = &S->isoplots[0]; // mu

	// Reset event selections to baseline:
	setRegionCuts(gRegion[gBaseRegion]);

	int muind1(-1), muind2(-1);
	if(hasLooseMuons(muind1, muind2) > 0){
		setHypLepton1(muind1, Muon);
		// Common trigger selection
		if(!singleMuTrigger()) return;
		float prescale = singleMuPrescale();
		float scale = prescale * gEventWeight;

		// Common event selections
		if(!passesJet50Cut()) return; // make trigger 100% efficient

		// Common object selections
		if(!isLooseMuon(muind1)) return;
		if(MuPt[muind1] < fC_minMu2pt) return;
		if(MuPt[muind1] > gMuFPtBins[gNMuFPtBins]) return;

		////////////////////////////////////////////////////
		// MOST LOOSE SELECTION
		IP0->hiso[0]->Fill(MuPFIso[muind1], scale);
		for(size_t k = 0; k < gNMuFPtBins; ++k){
			if(MuPt[muind1] < gMuFPtBins[k]) continue;
			if(MuPt[muind1] > gMuFPtBins[k+1]) continue;
			IP0->hiso_pt[0][k]->Fill(MuPFIso[muind1], scale);
		}
		for(size_t k = 0; k < gNNVrtxBins; ++k){
			if(NVrtx <  gNVrtxBins[k]) continue;
			if(NVrtx >= gNVrtxBins[k+1]) continue;
			IP0->hiso_nv[0][k]->Fill(MuPFIso[muind1], scale);
		}

		////////////////////////////////////////////////////
		// SIGNAL SUPPRESSED SELECTION
		if(isSigSupMuEvent() > -1){
			IP0->hiso[1]->Fill(MuPFIso[muind1], scale);
			for(size_t k = 0; k < gNMuFPtBins; ++k){
				if(MuPt[muind1] < gMuFPtBins[k]) continue;
				if(MuPt[muind1] > gMuFPtBins[k+1]) continue;
				IP0->hiso_pt[1][k]->Fill(MuPFIso[muind1], scale);
			}
			for(size_t k = 0; k < gNNVrtxBins; ++k){
				if(NVrtx <  gNVrtxBins[k]) continue;
				if(NVrtx >= gNVrtxBins[k+1]) continue;
				IP0->hiso_nv[1][k]->Fill(MuPFIso[muind1], scale);
			}
		}
		resetHypLeptons();
		////////////////////////////////////////////////////
	}
	setRegionCuts(gRegion[gBaseRegion]);
	return;
}
void SSDLDumper::fillElIdPlots(Sample *S){
	resetHypLeptons();
	fDoCounting = false;
	IdPlots *IdP = &S->idplots;

	// Reset event selections to baseline:
	setRegionCuts(gRegion[gBaseRegion]);

	int elind1(-1), elind2(-1);
	if(hasLooseElectrons(elind1, elind2) > 0){
		setHypLepton1(elind1, Elec);
		// Common trigger selection
		if(!singleElTrigger()) return;
		float prescale = singleElPrescale();
		float scale = prescale * gEventWeight;

		// Common event selections
		if(!passesJet50Cut()) return; // make trigger 100% efficient

		// Common object selections
		if(!isLooseElectron(elind1)) return;

		if(ElPt[elind1] < fC_minEl2pt) return;
		if(ElPt[elind1] > gElFPtBins[gNElFPtBins]) return;

		////////////////////////////////////////////////////
		// MOST LOOSE SELECTION
		IdP->hhoe   [0]->Fill(ElHoverE[elind1], scale);
		IdP->hsiesie[0]->Fill(ElSigmaIetaIeta[elind1], scale);
		IdP->hdeta  [0]->Fill(ElDEta[elind1], scale);
		IdP->hdphi  [0]->Fill(ElDPhi[elind1], scale);
		IdP->hmvaid [0]->Fill(ElMVAIDTrig[elind1], scale);
		IdP->hmedwp [0]->Fill(ElIsGoodElId_MediumWP[elind1], scale);
		//for(size_t k = 0; k < gNElFPtBins; ++k){
		//	if(ElPt[elind1] < gElFPtBins[k]) continue;
		//	if(ElPt[elind1] > gElFPtBins[k+1]) continue;
		//	IdP->hiso_pt[0][k]->Fill(ElHoverE[elind1], scale);
		//}
		//for(size_t k = 0; k < gNNVrtxBins; ++k){
		//	if(NVrtx <  gNVrtxBins[k]) continue;
		//	if(NVrtx >= gNVrtxBins[k+1]) continue;
		//	IdP->hiso_nv[0][k]->Fill(ElHoverE[elind1], scale);
		//}

		////////////////////////////////////////////////////
		// SIGNAL SUPPRESSED SELECTION
		if(isSigSupElEvent() > -1){
			IdP->hhoe   [1]->Fill(ElHoverE[elind1], scale);
			IdP->hsiesie[1]->Fill(ElSigmaIetaIeta[elind1], scale);
			IdP->hdeta  [1]->Fill(ElDEta[elind1], scale);
			IdP->hdphi  [1]->Fill(ElDPhi[elind1], scale);
			IdP->hmvaid [1]->Fill(ElMVAIDTrig[elind1], scale);
			IdP->hmedwp [1]->Fill(ElIsGoodElId_MediumWP[elind1], scale);
			//IP->hiso[1]->Fill(ElPFIso[elind1], scale);
			//for(size_t k = 0; k < gNElFPtBins; ++k){
			//	if(ElPt[elind1] < gElFPtBins[k]) continue;
			//	if(ElPt[elind1] > gElFPtBins[k+1]) continue;
			//	IP->hiso_pt[1][k]->Fill(ElPFIso[elind1], scale);
			//}
			//for(size_t k = 0; k < gNNVrtxBins; ++k){
			//	if(NVrtx <  gNVrtxBins[k]) continue;
			//	if(NVrtx >= gNVrtxBins[k+1]) continue;
			//	IP->hiso_nv[1][k]->Fill(ElPFIso[elind1], scale);
			//}
		}
		//////////////////////////////////////////////////
		resetHypLeptons();
	}
	setRegionCuts(gRegion[gBaseRegion]);
	return;
}
void SSDLDumper::fillElIsoPlots(Sample *S){
	resetHypLeptons();
	fDoCounting = false;
	IsoPlots *IP = &S->isoplots[1]; // el

	// Reset event selections to baseline:
	setRegionCuts(gRegion[gBaseRegion]);

	int elind1(-1), elind2(-1);
	if(hasLooseElectrons(elind1, elind2) > 0){
		setHypLepton1(elind1, Elec);
		// Common trigger selection
		if(!singleElTrigger()) return;
		float prescale = singleElPrescale();
		float scale = prescale * gEventWeight;

		// Common event selections
		if(!passesJet50Cut()) return; // make trigger 100% efficient

		// Common object selections
		if(!isLooseElectron(elind1)) return;

		if(ElPt[elind1] < fC_minEl2pt) return;
		if(ElPt[elind1] > gElFPtBins[gNElFPtBins]) return;

		////////////////////////////////////////////////////
		// MOST LOOSE SELECTION
		IP->hiso[0]->Fill(ElPFIso[elind1], scale);
		for(size_t k = 0; k < gNElFPtBins; ++k){
			if(ElPt[elind1] < gElFPtBins[k]) continue;
			if(ElPt[elind1] > gElFPtBins[k+1]) continue;
			IP->hiso_pt[0][k]->Fill(ElPFIso[elind1], scale);
		}
		for(size_t k = 0; k < gNNVrtxBins; ++k){
			if(NVrtx <  gNVrtxBins[k]) continue;
			if(NVrtx >= gNVrtxBins[k+1]) continue;
			IP->hiso_nv[0][k]->Fill(ElPFIso[elind1], scale);
		}

		////////////////////////////////////////////////////
		// SIGNAL SUPPRESSED SELECTION
		if(isSigSupElEvent() > -1){
			IP->hiso[1]->Fill(ElPFIso[elind1], scale);
			for(size_t k = 0; k < gNElFPtBins; ++k){
				if(ElPt[elind1] < gElFPtBins[k]) continue;
				if(ElPt[elind1] > gElFPtBins[k+1]) continue;
				IP->hiso_pt[1][k]->Fill(ElPFIso[elind1], scale);
			}
			for(size_t k = 0; k < gNNVrtxBins; ++k){
				if(NVrtx <  gNVrtxBins[k]) continue;
				if(NVrtx >= gNVrtxBins[k+1]) continue;
				IP->hiso_nv[1][k]->Fill(ElPFIso[elind1], scale);
			}
		}
		////////////////////////////////////////////////////
		resetHypLeptons();
	}
	setRegionCuts(gRegion[gBaseRegion]);
	return;
}
void SSDLDumper::fillSyncCounters(Sample *S){
  resetHypLeptons();
  if (!gDoSyncExercise) return;
  setRegionCuts(gRegion[gBaseRegion]);
  //////////////////////////////////////////////////////////////////
  fCurrentChannel = Muon;
  int mu1(-1), mu2(-1);
  if(abs(isOSLLEvent(mu1, mu2)) == 1) {
    fCounterSync[Muon].fill(fSyncCutNames[0]);
    
    if ((S->sname).Contains("DoubleMu"))
      fOUTSTREAM << Run <<" "<< LumiSec <<" "<< Event <<" "<< getNJets() <<" "
		 << getNBTagsMed() <<" "<< (pfMET>10) <<" "<<mumuSignalTrigger() <<" "<< elelSignalTrigger() 
		 <<" "<< elmuSignalTrigger() << endl;
    
    if (mumuSignalTrigger()) fCounterSync[Muon].fill(fSyncCutNames[1]);
    
    setHypLepton1(mu1, Muon);
    setHypLepton2(mu2, Muon);
    
    if (getNBTagsMed() == 0) fCounterSync[Muon].fill(fSyncCutNames[2]);
    if (getNJets() == 0)     fCounterSync[Muon].fill(fSyncCutNames[3]);
    if (getNJets() == 1)     fCounterSync[Muon].fill(fSyncCutNames[4]);
    if (getNJets() >= 2)     fCounterSync[Muon].fill(fSyncCutNames[5]);
    
    if (pfMET      > 10)     fCounterSync[Muon].fill(fSyncCutNames[6]);
    if (pfMETType1 > 10)     fCounterSync[Muon].fill(fSyncCutNames[7]);
  }
  //////////////////////////////////////////////////////////////////
  
  //////////////////////////////////////////////////////////////////
  fCurrentChannel = Elec;
  resetHypLeptons();
  int el1(-1), el2(-1);
  if(abs(isOSLLEvent(el1, el2)) == 2) {
    fCounterSync[Elec].fill(fSyncCutNames[0]);
    if ((S->sname).Contains("DoubleEl"))
      fOUTSTREAM << Run <<" "<< LumiSec <<" "<< Event <<" "<< getNJets() <<" "
		 << getNBTagsMed() <<" "<< (pfMET>10) <<" "<<mumuSignalTrigger() <<" "<< elelSignalTrigger() 
		 <<" "<< elmuSignalTrigger() << endl;
    
    if (elelSignalTrigger()) fCounterSync[Elec].fill(fSyncCutNames[1]);
    setHypLepton1(el1, Elec);
    setHypLepton2(el2, Elec);
    
    if (getNBTagsMed() == 0) fCounterSync[Elec].fill(fSyncCutNames[2]);
    if (getNJets() == 0)     fCounterSync[Elec].fill(fSyncCutNames[3]);
    if (getNJets() == 1)     fCounterSync[Elec].fill(fSyncCutNames[4]);
    if (getNJets() >= 2)     fCounterSync[Elec].fill(fSyncCutNames[5]);
    
    if (pfMET      > 10)     fCounterSync[Elec].fill(fSyncCutNames[6]);
    if (pfMETType1 > 10)     fCounterSync[Elec].fill(fSyncCutNames[7]);
  }
  //////////////////////////////////////////////////////////////////
  
  //////////////////////////////////////////////////////////////////
  fCurrentChannel = ElMu;
  int mu(-1), el(-1);
  resetHypLeptons();
  if(abs(isOSLLEvent(mu, el)) == 3) {
    fCounterSync[ElMu].fill(fSyncCutNames[0]);
    if ((S->sname).Contains("MuEG"))
      fOUTSTREAM << Run <<" "<< LumiSec <<" "<< Event <<" "<< getNJets() <<" "
		 << getNBTagsMed() <<" "<< (pfMET>10) <<" "<<mumuSignalTrigger() <<" "<< elelSignalTrigger() 
		 <<" "<< elmuSignalTrigger() << endl;
    if (elmuSignalTrigger()) fCounterSync[ElMu].fill(fSyncCutNames[1]);
    
    setHypLepton1(mu, Muon);
    setHypLepton2(el, Elec);
    
    if (getNBTagsMed() == 0) fCounterSync[ElMu].fill(fSyncCutNames[2]);
    if (getNJets() == 0)     fCounterSync[ElMu].fill(fSyncCutNames[3]);
    if (getNJets() == 1)     fCounterSync[ElMu].fill(fSyncCutNames[4]);
    if (getNJets() >= 2)     fCounterSync[ElMu].fill(fSyncCutNames[5]);
    
    if (pfMET      > 10)     fCounterSync[ElMu].fill(fSyncCutNames[6]);
    if (pfMETType1 > 10)     fCounterSync[ElMu].fill(fSyncCutNames[7]);
  }
  //////////////////////////////////////////////////////////////////

}
//____________________________________________________________________________
void SSDLDumper::storeNumbers(Sample *S, gChannel chan, int reg){
	Channel *C;
	if(chan == Muon) C = &S->region[reg][HighPt].mm;
	if(chan == Elec) C = &S->region[reg][HighPt].ee;
	if(chan == ElMu) C = &S->region[reg][HighPt].em;
	S->numbers[reg][chan].npp = 0.;
	S->numbers[reg][chan].npf = 0.;
	S->numbers[reg][chan].nfp = 0.;
	S->numbers[reg][chan].nff = 0.;
	S->numbers[reg][chan].nt2  = C->nt20_pt->Integral(0, getNFPtBins(chan)+1);
	S->numbers[reg][chan].nt10 = C->nt10_pt->Integral(0, getNFPtBins(chan)+1);
	S->numbers[reg][chan].nt01 = C->nt01_pt->Integral(0, getNFPtBins(chan)+1);
	S->numbers[reg][chan].nt0  = C->nt00_pt->Integral(0, getNFPtBins(chan)+1);
	
	S->numbers[reg][chan].tt_avweight = C->nt20_pt->GetEntries()>0?C->nt20_pt->Integral(0, getNFPtBins(chan)+1)/C->nt20_pt->GetEntries():1.;
	S->numbers[reg][chan].tl_avweight = C->nt10_pt->GetEntries()>0?C->nt10_pt->Integral(0, getNFPtBins(chan)+1)/C->nt10_pt->GetEntries():1.;
	S->numbers[reg][chan].lt_avweight = C->nt01_pt->GetEntries()>0?C->nt01_pt->Integral(0, getNFPtBins(chan)+1)/C->nt01_pt->GetEntries():1.;
	S->numbers[reg][chan].ll_avweight = C->nt00_pt->GetEntries()>0?C->nt00_pt->Integral(0, getNFPtBins(chan)+1)/C->nt00_pt->GetEntries():1.;
	
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
	fMMCutNames.push_back(" ... passes 3rd lep veto"); //              = 8
	fMMCutNames.push_back(" ... passes NJets cut"); //                 = 9
	fMMCutNames.push_back(" ... passes NbJets cut"); //                = 10
	fMMCutNames.push_back(" ... passes HT cut"); //                    = 11
	fMMCutNames.push_back(" ... passes MET cut"); //                   = 12
	fMMCutNames.push_back(" ... second muon passes pt cut"); //        = 13
	fMMCutNames.push_back(" ... first muon passes pt cut"); //         = 14
	fMMCutNames.push_back(" ... first muon passes tight cut"); //      = 15
	fMMCutNames.push_back(" ... second muon passes tight cut"); //     = 16
	fMMCutNames.push_back(" ... both muons pass tight cut"); //        = 17

	// Electron channel
	fEECutNames.push_back("All events"); //                            = 0
	fEECutNames.push_back(" ... is good run"); //                      = 1
	fEECutNames.push_back(" ... passes triggers"); //                  = 2
	fEECutNames.push_back(" ... has 1 loose electron"); //             = 3
	fEECutNames.push_back(" ... has 2 loose electrons"); //            = 4
	fEECutNames.push_back(" ... has same-sign electrons"); //          = 5
	fEECutNames.push_back(" ... passes Z veto"); //                    = 6
	fEECutNames.push_back(" ... passes Minv veto"); //                 = 7
	fEECutNames.push_back(" ... passes 3rd lep veto"); //              = 8
	fEECutNames.push_back(" ... passes NJets cut"); //                 = 9
	fEECutNames.push_back(" ... passes NbJets cut"); //                = 10
	fEECutNames.push_back(" ... passes HT cut"); //                    = 11
	fEECutNames.push_back(" ... passes MET cut"); //                   = 12
	fEECutNames.push_back(" ... second electron passes pt cut"); //    = 13
	fEECutNames.push_back(" ... first electron passes pt cut"); //     = 14
	fEECutNames.push_back(" ... first electron passes tight cut"); //  = 15
	fEECutNames.push_back(" ... second electron passes tight cut"); // = 16
	fEECutNames.push_back(" ... both electrons pass tight cut"); //    = 17
	
	// E-Mu channel
	fEMCutNames.push_back("All events"); //                            = 0
	fEMCutNames.push_back(" ... is good run"); //                      = 1
	fEMCutNames.push_back(" ... passes triggers"); //                  = 2
	fEMCutNames.push_back(" ... has a loose muon"); //                 = 3
	fEMCutNames.push_back(" ... has a loose electron"); //             = 4
	fEMCutNames.push_back(" ... has both"); //                         = 5
	fEMCutNames.push_back(" ... has same-sign electron muon pair"); // = 6
	fEMCutNames.push_back(" ... passes Z veto"); //                    = 7
	fEMCutNames.push_back(" ... passes 3rd lep veto"); //              = 8
	fEMCutNames.push_back(" ... passes NJets cut"); //                 = 9
	fEMCutNames.push_back(" ... passes NbJets cut"); //                = 10
	fEMCutNames.push_back(" ... passes HT cut"); //                    = 11
	fEMCutNames.push_back(" ... passes MET cut"); //                   = 12
	fEMCutNames.push_back(" ... muon passes pt cut"); //               = 13
	fEMCutNames.push_back(" ... electron passes pt cut"); //           = 14
	fEMCutNames.push_back(" ... muon passes tight cut"); //            = 15
	fEMCutNames.push_back(" ... electron passes tight cut"); //        = 16
	fEMCutNames.push_back(" ... both e and mu pass tight cuts"); //    = 17

	fSyncCutNames.push_back("             "); //0
	fSyncCutNames.push_back(" + triggers  "); //1
	fSyncCutNames.push_back(" + nbjets = 0"); //2
	fSyncCutNames.push_back(" + njets  = 0"); //3
	fSyncCutNames.push_back(" + njets  = 1"); //4
	fSyncCutNames.push_back(" + njets >= 2"); //5
	fSyncCutNames.push_back(" + pfmet > 10"); //6
	fSyncCutNames.push_back(" + t1met > 10"); //7

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
	fCounter[Muon].fill(fMMCutNames[17], 0.);

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
	fCounter[Elec].fill(fEECutNames[17], 0.);

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
	fCounter[ElMu].fill(fEMCutNames[17], 0.);
 
	fCounterSync[Muon].fill(fSyncCutNames[0], 0.);
	fCounterSync[Muon].fill(fSyncCutNames[1], 0.);
	fCounterSync[Muon].fill(fSyncCutNames[2], 0.);
	fCounterSync[Muon].fill(fSyncCutNames[3], 0.);
	fCounterSync[Muon].fill(fSyncCutNames[4], 0.);
	fCounterSync[Muon].fill(fSyncCutNames[5], 0.);
	fCounterSync[Muon].fill(fSyncCutNames[6], 0.);
	fCounterSync[Muon].fill(fSyncCutNames[7], 0.);

	fCounterSync[Elec].fill(fSyncCutNames[0], 0.);
	fCounterSync[Elec].fill(fSyncCutNames[1], 0.);
	fCounterSync[Elec].fill(fSyncCutNames[2], 0.);
	fCounterSync[Elec].fill(fSyncCutNames[3], 0.);
	fCounterSync[Elec].fill(fSyncCutNames[4], 0.);
	fCounterSync[Elec].fill(fSyncCutNames[5], 0.);
	fCounterSync[Elec].fill(fSyncCutNames[6], 0.);
	fCounterSync[Elec].fill(fSyncCutNames[7], 0.);

	fCounterSync[ElMu].fill(fSyncCutNames[0], 0.);
	fCounterSync[ElMu].fill(fSyncCutNames[1], 0.);
	fCounterSync[ElMu].fill(fSyncCutNames[2], 0.);
	fCounterSync[ElMu].fill(fSyncCutNames[3], 0.);
	fCounterSync[ElMu].fill(fSyncCutNames[4], 0.);
	fCounterSync[ElMu].fill(fSyncCutNames[5], 0.);
	fCounterSync[ElMu].fill(fSyncCutNames[6], 0.);
	fCounterSync[ElMu].fill(fSyncCutNames[7], 0.);
	
}
void SSDLDumper::fillCutFlowHistos(Sample *S){
	for(int i=0; i<fMMCutNames.size(); i++) S->cutFlowHisto[Muon]->SetBinContent(i+1, fCounter[Muon].counts(fMMCutNames[i]));
	for(int i=0; i<fEECutNames.size(); i++) S->cutFlowHisto[Elec]->SetBinContent(i+1, fCounter[Elec].counts(fEECutNames[i]));
	for(int i=0; i<fEMCutNames.size(); i++) S->cutFlowHisto[ElMu]->SetBinContent(i+1, fCounter[ElMu].counts(fEMCutNames[i]));	
}
void SSDLDumper::printCutFlow(gChannel chan){
  TString ch = "";
  if(chan == Muon) ch = "MM";
  if(chan == Elec) ch = "EE";
  if(chan == ElMu) ch = "EM";

  fOUTSTREAM << "--------------------------------------------------------" << endl;
  fOUTSTREAM << " Cutname                                 | " << endl;  
  for( int c = 0; c < fSyncCutNames.size(); c++ ){
    fOUTSTREAM << setw(40) << ch + fSyncCutNames[c] << " | ";
    fOUTSTREAM << setw(11) << setprecision(11) << fCounterSync[chan].counts(fSyncCutNames[c]) << endl;
  }
  
}

void SSDLDumper::printCutFlow(gChannel chan, gSample indmin, gSample indmax){
	vector<string> names;
	if(chan == Muon) names = fMMCutNames;
	if(chan == Elec) names = fEECutNames;
	if(chan == ElMu) names = fEMCutNames;

	fOUTSTREAM << "------------------------------------------";
	for(int i = 0; i <= indmax-indmin; i++) fOUTSTREAM << "--------------";
	fOUTSTREAM << endl;


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
	printCutFlow(Muon, TTJets, GJets400);
	printCutFlow(Muon, WW, WmWm);
	// printCutFlow(Muon, LM0, LM6);
	// printCutFlow(Muon, LM7, LM13);
	printCutFlow(Muon, QCDMuEnr15, QCD800);
	// printCutFlow(Muon, QCD600, gNSAMPLES);
	fOUTSTREAM << endl << endl;
	
	fOUTSTREAM << " Printing Cutflow for E/Mu channel..." << endl;
	printCutFlow(ElMu, MuEG1, MuEG4);
	printCutFlow(ElMu, TTJets, GJets400);
	printCutFlow(ElMu, WW, WmWm);
	// printCutFlow(ElMu, LM0, LM6);
	// printCutFlow(ElMu, LM7, LM13);
	printCutFlow(ElMu, QCDMuEnr15, QCD800);
	// printCutFlow(ElMu, QCD600, gNSAMPLES);
	fOUTSTREAM << endl << endl;
	
	fOUTSTREAM << " Printing Cutflow for E/E channel..." << endl;
	printCutFlow(Elec, DoubleEle1, DoubleEle4);
	printCutFlow(Elec, TTJets, GJets400);
	printCutFlow(Elec, WW, WmWm);
	// printCutFlow(Elec, LM0, LM6);
	// printCutFlow(Elec, LM7, LM13);
	printCutFlow(Elec, QCDEM20, QCDEM250);
	// printCutFlow(Elec, QCD600, gNSAMPLES);
	fOUTSTREAM << endl << endl;

	fOUTSTREAM.close();
}

//____________________________________________________________________________
void SSDLDumper::bookSigEvTree(){
	fSigEv_Tree = new TTree("SigEvents", "SigEventTree");
	fSigEv_Tree->Branch("PUWeight",    &fSETree_PUWeight, "PUWeight/F");
	fSigEv_Tree->Branch("SystFlag",    &fSETree_SystFlag, "SystFlag/I");
	fSigEv_Tree->Branch("HLTSF",       &fSETree_HLTSF   , "HLTSF/F");
	fSigEv_Tree->Branch("BtagSF1",      &fSETree_BtagSF1  , "BtagSF1/F");
	fSigEv_Tree->Branch("BtagSF2",      &fSETree_BtagSF2  , "BtagSF2/F");
	fSigEv_Tree->Branch("SLumi",       &fSETree_SLumi   , "SLumi/F");
	fSigEv_Tree->Branch("SName",       &fSETree_SName);
	fSigEv_Tree->Branch("SType",       &fSETree_SType   , "SType/I");
	fSigEv_Tree->Branch("Run",         &fSETree_Run     , "Run/I");
	fSigEv_Tree->Branch("LS",          &fSETree_LS      , "LS/I");
	fSigEv_Tree->Branch("Event",       &fSETree_Event   , "Event/I");
	fSigEv_Tree->Branch("Flavor",      &fSETree_Flavor  , "Flavor/I");
	fSigEv_Tree->Branch("Charge",      &fSETree_Charge  , "Charge/I");
	fSigEv_Tree->Branch("TLCat",       &fSETree_TLCat   , "TLCat/I");
	fSigEv_Tree->Branch("Pass3rdVeto", &fSETree_3rdVeto , "Pass3rdVeto/I");
	fSigEv_Tree->Branch("PassTTZSel",  &fSETree_ttZSel  , "PassTTZSel/I");
	fSigEv_Tree->Branch("PassZVeto",   &fSETree_ZVeto   , "PassZVeto/I");
	fSigEv_Tree->Branch("HT",          &fSETree_HT      , "HT/F");
	fSigEv_Tree->Branch("MET",         &fSETree_MET     , "MET/F");
	fSigEv_Tree->Branch("NM",          &fSETree_NM      , "NM/I");
	fSigEv_Tree->Branch("NE",          &fSETree_NE      , "NE/I");
	fSigEv_Tree->Branch("NJ",          &fSETree_NJ      , "NJ/I");
	fSigEv_Tree->Branch("NbJ",         &fSETree_NbJ     , "NbJ/I");
	fSigEv_Tree->Branch("NbJmed",      &fSETree_NbJmed  , "NbJmed/I");
	fSigEv_Tree->Branch("M3",          &fSETree_M3      , "M3/F");
	fSigEv_Tree->Branch("MT2",         &fSETree_MT2     , "MT2/F");
	fSigEv_Tree->Branch("Mll",         &fSETree_Mll     , "Mll/F");
	fSigEv_Tree->Branch("pT1",         &fSETree_pT1     , "pT1/F");
	fSigEv_Tree->Branch("pT2",         &fSETree_pT2     , "pT2/F");
	fSigEv_Tree->Branch("eta1",        &fSETree_eta1    , "eta1/F");
	fSigEv_Tree->Branch("eta2",        &fSETree_eta2    , "eta2/F");
	fSigEv_Tree->Branch("PFIso1",      &fSETree_PFIso1  , "PFIso1/F");
	fSigEv_Tree->Branch("PFIso2",      &fSETree_PFIso2  , "PFIso2/F");
	fSigEv_Tree->Branch("MVAID1",      &fSETree_MVAID1  , "MVAID1/F");
	fSigEv_Tree->Branch("MVAID2",      &fSETree_MVAID2  , "MVAID2/F");
	fSigEv_Tree->Branch("medWP1",      &fSETree_medWP1  , "medWP1/F");
	fSigEv_Tree->Branch("medWP2",      &fSETree_medWP2  , "medWP2/F");
}
void SSDLDumper::resetSigEventTree(){
	fSETree_SystFlag = -1;
	fSETree_PUWeight = -1.;
	fSETree_HLTSF    = -1.;
	fSETree_BtagSF1   = -1.;
	fSETree_BtagSF2   = -1.;
	fSETree_SLumi    = -1.;
	fSETree_SName    = "?";
	fSETree_SType    = -1;
	fSETree_Run      = -1;
	fSETree_LS       = -1;
	fSETree_Event    = -1;
	fSETree_Flavor   = -1;
	fSETree_Charge   = -99;
	fSETree_TLCat    = -1;
	fSETree_ZVeto    = -1;
	fSETree_3rdVeto  = -1;
	fSETree_ttZSel   = -1;
	fSETree_HT       = -1.;
	fSETree_MET      = -1.;
	fSETree_NM       = -1;
	fSETree_NE       = -1;
	fSETree_NJ       = -1;
	fSETree_NbJ      = -1;
	fSETree_NbJmed   = -1;
	fSETree_M3       = -1.;
	fSETree_MT2      = -1.;
	fSETree_Mll      = -1.;
	fSETree_pT1      = -1.;
	fSETree_pT2      = -1.;
	fSETree_eta1     = -1.;
	fSETree_eta2     = -1.;
	fSETree_PFIso1   = -1.;
	fSETree_PFIso2   = -1.;
	fSETree_MVAID1   = -999.;
	fSETree_MVAID2   = -999.;
	fSETree_medWP1   = -1.;
	fSETree_medWP2   = -1.;
}
void SSDLDumper::writeSigEvTree(TFile *pFile){
	pFile->cd();
	fSigEv_Tree->Write("SigEvents", TObject::kWriteDelete);	
}
int SSDLDumper::getSampleType(Sample *S){
	if(S->datamc == 2) return 20; // signal
	if(S->datamc > 0 && (S->getType() != 4 && S->getType() != 5) ) return 10;  // SM MC without RARE and di-BOSON
	if(S->getType() == 4 || S->getType() == 5  ) return 15;  // RARE MC + di-BOSON
	if(S->sname.Contains("DoubleMu"))  return 0; // Data
	if(S->sname.Contains("DoubleEle")) return 1;
	if(S->sname.Contains("MuEG"))      return 2;
	if(S->sname.Contains("MuHad"))     return 3;
	if(S->sname.Contains("EleHad"))    return 4;
	else return -1;
}

//____________________________________________________________________________
void SSDLDumper::bookHistos(Sample *S){
	// Event counter
	S->evcount = new TH1F(S->sname + "_EventCount", "Event Counter", 1, 0., 1.);	
	
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

			name = Form("%s_%s_NPP_%s", S->sname.Data(), DiffPredYields::var_name[j].Data(), gChanLabel[k].Data());
			S->diffyields[k].hnpp[j] = new TH1D(name, DiffPredYields::var_name[j].Data(), DiffPredYields::nbins[j], DiffPredYields::bins[j]);
			S->diffyields[k].hnpp[j]->SetFillColor(S->color);
			S->diffyields[k].hnpp[j]->SetXTitle(DiffPredYields::axis_label[j]);
			S->diffyields[k].hnpp[j]->Sumw2();
			name = Form("%s_%s_NPF_%s", S->sname.Data(), DiffPredYields::var_name[j].Data(), gChanLabel[k].Data());
			S->diffyields[k].hnpf[j] = new TH1D(name, DiffPredYields::var_name[j].Data(), DiffPredYields::nbins[j], DiffPredYields::bins[j]);
			S->diffyields[k].hnpf[j]->SetFillColor(S->color);
			S->diffyields[k].hnpf[j]->SetXTitle(DiffPredYields::axis_label[j]);
			S->diffyields[k].hnpf[j]->Sumw2();
			name = Form("%s_%s_NFP_%s", S->sname.Data(), DiffPredYields::var_name[j].Data(), gChanLabel[k].Data());
			S->diffyields[k].hnfp[j] = new TH1D(name, DiffPredYields::var_name[j].Data(), DiffPredYields::nbins[j], DiffPredYields::bins[j]);
			S->diffyields[k].hnfp[j]->SetFillColor(S->color);
			S->diffyields[k].hnfp[j]->SetXTitle(DiffPredYields::axis_label[j]);
			S->diffyields[k].hnfp[j]->Sumw2();
			name = Form("%s_%s_NFF_%s", S->sname.Data(), DiffPredYields::var_name[j].Data(), gChanLabel[k].Data());
			S->diffyields[k].hnff[j] = new TH1D(name, DiffPredYields::var_name[j].Data(), DiffPredYields::nbins[j], DiffPredYields::bins[j]);
			S->diffyields[k].hnff[j]->SetFillColor(S->color);
			S->diffyields[k].hnff[j]->SetXTitle(DiffPredYields::axis_label[j]);
			S->diffyields[k].hnff[j]->Sumw2();

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
	for(size_t k = 0; k < gNKinSels; ++k){
		TString name = Form("%s_%s_HTvsMET", S->sname.Data(), gKinSelNames[k].Data());
		for(size_t j = 0; j < gNKinVars; ++j){
			name = Form("%s_%s_%s", S->sname.Data(), gKinSelNames[k].Data(), KinPlots::var_name[j].Data());
			S->kinplots[k][HighPt].hvar[j] = new TH1D(name, KinPlots::var_name[j].Data(), KinPlots::nbins[j], KinPlots::xmin[j], KinPlots::xmax[j]);
			S->kinplots[k][HighPt].hvar[j]->SetFillColor(S->color);
			S->kinplots[k][HighPt].hvar[j]->SetXTitle(KinPlots::axis_label[j]);
			S->kinplots[k][HighPt].hvar[j]->Sumw2();
		}
	}
	
	// WZ Kinematical histos
	for (size_t k = 0; k < gNKinSels; ++k){
		for(size_t j = 0; j < gNKinVars; ++j){
			TString name = Form("%s_%s_%s", S->sname.Data(), gKinSelNames[k].Data(), KinPlots::var_name[j].Data());
			S->kinplots_wz[k].hvar[j] = new TH1D(name, KinPlots::var_name[j].Data(), KinPlots::nbins[j], KinPlots::xmin[j], KinPlots::xmax[j]);
			S->kinplots_wz[k].hvar[j]->SetFillColor(S->color);
			S->kinplots_wz[k].hvar[j]->SetXTitle(KinPlots::axis_label[j]);
			S->kinplots_wz[k].hvar[j]->Sumw2();
		}
	}
	
	// id histos for electrons
	for(size_t j = 0; j < gNSels; ++j){
		TString hoename    = Form("%s_%s_%shoe"   , S->sname.Data(), IdPlots::sel_name[j].Data(), gEMULabel[1].Data());
		TString siesiename = Form("%s_%s_%ssiesie", S->sname.Data(), IdPlots::sel_name[j].Data(), gEMULabel[1].Data());
		TString detaname   = Form("%s_%s_%sdeta"  , S->sname.Data(), IdPlots::sel_name[j].Data(), gEMULabel[1].Data());
		TString dphiname   = Form("%s_%s_%sdphi"  , S->sname.Data(), IdPlots::sel_name[j].Data(), gEMULabel[1].Data());
		TString mvaidname  = Form("%s_%s_%smvaid" , S->sname.Data(), IdPlots::sel_name[j].Data(), gEMULabel[1].Data());
		TString medwpname  = Form("%s_%s_%smedwp" , S->sname.Data(), IdPlots::sel_name[j].Data(), gEMULabel[1].Data());
		S->idplots.hhoe[j] = new TH1D(hoename, Form("%shoe", gEMULabel[1].Data()), IdPlots::nbins[j], 0., 0.15);
		S->idplots.hhoe[j]->SetFillColor(S->color);
		S->idplots.hhoe[j]->SetXTitle("HoE");
		S->idplots.hhoe[j]->Sumw2();
		S->idplots.hsiesie[j] = new TH1D(siesiename, Form("%ssiesie", gEMULabel[1].Data()), IdPlots::nbins[j], 0., 0.035);
		S->idplots.hsiesie[j]->SetFillColor(S->color);
		S->idplots.hsiesie[j]->SetXTitle("#sigma_{i#eta , i#eta}");
		S->idplots.hsiesie[j]->Sumw2();
		S->idplots.hdeta[j] = new TH1D(detaname, Form("%sdeta", gEMULabel[1].Data()), IdPlots::nbins[j], -0.01, 0.01);
		S->idplots.hdeta[j]->SetFillColor(S->color);
		S->idplots.hdeta[j]->SetXTitle("#Delta#eta");
		S->idplots.hdeta[j]->Sumw2();
		S->idplots.hdphi[j] = new TH1D(dphiname, Form("%sdphi", gEMULabel[1].Data()), IdPlots::nbins[j], -0.15, 0.15);
		S->idplots.hdphi[j]->SetFillColor(S->color);
		S->idplots.hdphi[j]->SetXTitle("#Delta#phi");
		S->idplots.hdphi[j]->Sumw2();
		S->idplots.hmvaid[j] = new TH1D(mvaidname, Form("%smvaID", gEMULabel[1].Data()), IdPlots::nbins[j], -1., 1.);
		S->idplots.hmvaid[j]->SetFillColor(S->color);
		S->idplots.hmvaid[j]->SetXTitle("MVA-ID disc.");
		S->idplots.hmvaid[j]->Sumw2();
		S->idplots.hmedwp[j] = new TH1D(medwpname, Form("%smedwp", gEMULabel[1].Data()), 3, 0, 2);
		S->idplots.hmedwp[j]->SetFillColor(S->color);
		S->idplots.hmedwp[j]->SetXTitle("medium WP");
		S->idplots.hmedwp[j]->Sumw2();
	}
	// isolation histos for electrons
	for(size_t j = 0; j < gNSels; ++j){
		TString name = Form("%s_%s_%siso", S->sname.Data(), IsoPlots::sel_name[j].Data(), gEMULabel[1].Data());
		S->isoplots[1].hiso[j] = new TH1D(name, Form("%siso", gEMULabel[1].Data()), IsoPlots::nbins[j], 0., 0.6);
		S->isoplots[1].hiso[j]->SetFillColor(S->color);
		S->isoplots[1].hiso[j]->SetXTitle("RelIso");
		S->isoplots[1].hiso[j]->Sumw2();
		for(int k = 0; k < gNMuFPtBins; ++k){
			name = Form("%s_%s_%siso_pt%d", S->sname.Data(), IsoPlots::sel_name[j].Data(), gEMULabel[1].Data(), k);
			S->isoplots[1].hiso_pt[j][k] = new TH1D(name, Form("%siso_pt%d", gEMULabel[1].Data(), k), IsoPlots::nbins[j], 0., 0.6);
			S->isoplots[1].hiso_pt[j][k]->SetFillColor(S->color);
			S->isoplots[1].hiso_pt[j][k]->SetXTitle("RelIso");
			S->isoplots[1].hiso_pt[j][k]->Sumw2();
		}
		for(int k = 0; k < gNNVrtxBins; ++k){
			name = Form("%s_%s_%siso_nv%d", S->sname.Data(), IsoPlots::sel_name[j].Data(), gEMULabel[1].Data(), k);
			S->isoplots[1].hiso_nv[j][k] = new TH1D(name, Form("%siso_nv%d", gEMULabel[1].Data(), k), IsoPlots::nbins[j], 0., 0.6);
			S->isoplots[1].hiso_nv[j][k]->SetFillColor(S->color);
			S->isoplots[1].hiso_nv[j][k]->SetXTitle("RelIso");
			S->isoplots[1].hiso_nv[j][k]->Sumw2();
		}
	}
	// isolation histos for muons
	for(size_t j = 0; j < gNSels; ++j){
		TString name = Form("%s_%s_%siso", S->sname.Data(), IsoPlots::sel_name[j].Data(), gEMULabel[0].Data());
		S->isoplots[0].hiso[j] = new TH1D(name, Form("%siso", gEMULabel[0].Data()), IsoPlots::nbins[j], 0., 1.0);
		S->isoplots[0].hiso[j]->SetFillColor(S->color);
		S->isoplots[0].hiso[j]->SetXTitle("RelIso");
		S->isoplots[0].hiso[j]->Sumw2();
		for(int k = 0; k < gNMuFPtBins; ++k){
			name = Form("%s_%s_%siso_pt%d", S->sname.Data(), IsoPlots::sel_name[j].Data(), gEMULabel[0].Data(), k);
			S->isoplots[0].hiso_pt[j][k] = new TH1D(name, Form("%siso_pt%d", gEMULabel[0].Data(), k), IsoPlots::nbins[j], 0., 1.0);
			S->isoplots[0].hiso_pt[j][k]->SetFillColor(S->color);
			S->isoplots[0].hiso_pt[j][k]->SetXTitle("RelIso");
			S->isoplots[0].hiso_pt[j][k]->Sumw2();
		}
		for(int k = 0; k < gNNVrtxBins; ++k){
			name = Form("%s_%s_%siso_nv%d", S->sname.Data(), IsoPlots::sel_name[j].Data(), gEMULabel[0].Data(), k);
			S->isoplots[0].hiso_nv[j][k] = new TH1D(name, Form("%siso_nv%d", gEMULabel[0].Data(), k), IsoPlots::nbins[j], 0., 1.0);
			S->isoplots[0].hiso_nv[j][k]->SetFillColor(S->color);
			S->isoplots[0].hiso_nv[j][k]->SetXTitle("RelIso");
			S->isoplots[0].hiso_nv[j][k]->Sumw2();
		}
	}


	for(size_t l = 0; l < 2; ++l){
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

	for(regIt = gRegions.begin(); regIt != gRegions.end(); regIt++){
		int r = gRegion[(*regIt)->sname];
		Region *R = &S->region[r][HighPt];
		for(gChannel c = channels_begin; c < gNCHANNELS; c=gChannel(c+1)){
			Channel *C;
			if(c == Muon) C = &R->mm;
			if(c == Elec) C = &R->ee;
			if(c == ElMu) C = &R->em;
			TString rootname = S->sname + "_" + (*regIt)->sname + "_" + gChanLabel[c];
			// Yields common for all channels and data-mc:
			C->nt20_pt  = new TH2D(rootname + "_NT20_pt",  "NT20_pt",  getNFPtBins(c), getFPtBins(c), getNFPtBins(c), getFPtBins(c)); C->nt20_pt ->Sumw2();
			C->nt10_pt  = new TH2D(rootname + "_NT10_pt",  "NT10_pt",  getNFPtBins(c), getFPtBins(c), getNFPtBins(c), getFPtBins(c)); C->nt10_pt ->Sumw2();
			C->nt01_pt  = new TH2D(rootname + "_NT01_pt",  "NT01_pt",  getNFPtBins(c), getFPtBins(c), getNFPtBins(c), getFPtBins(c)); C->nt01_pt ->Sumw2();
			C->nt00_pt  = new TH2D(rootname + "_NT00_pt",  "NT00_pt",  getNFPtBins(c), getFPtBins(c), getNFPtBins(c), getFPtBins(c)); C->nt00_pt ->Sumw2();
			C->nt20_eta = new TH2D(rootname + "_NT20_eta", "NT20_eta", getNEtaBins(c), getEtaBins(c), getNEtaBins(c), getEtaBins(c)); C->nt20_eta->Sumw2();
			C->nt10_eta = new TH2D(rootname + "_NT10_eta", "NT10_eta", getNEtaBins(c), getEtaBins(c), getNEtaBins(c), getEtaBins(c)); C->nt10_eta->Sumw2();
			C->nt01_eta = new TH2D(rootname + "_NT01_eta", "NT01_eta", getNEtaBins(c), getEtaBins(c), getNEtaBins(c), getEtaBins(c)); C->nt01_eta->Sumw2();
			C->nt00_eta = new TH2D(rootname + "_NT00_eta", "NT00_eta", getNEtaBins(c), getEtaBins(c), getNEtaBins(c), getEtaBins(c)); C->nt00_eta->Sumw2();

			// MC truth info
			if(S->datamc > 0){
				C->npp_pt   = new TH2D(rootname + "_NPP_pt",   "NPP_pt",   getNFPtBins(c), getFPtBins(c), getNFPtBins(c), getFPtBins(c)); C->npp_pt->Sumw2();
				C->nfp_pt   = new TH2D(rootname + "_NFP_pt",   "NFP_pt",   getNFPtBins(c), getFPtBins(c), getNFPtBins(c), getFPtBins(c)); C->nfp_pt->Sumw2();
				C->npf_pt   = new TH2D(rootname + "_NPF_pt",   "NPF_pt",   getNFPtBins(c), getFPtBins(c), getNFPtBins(c), getFPtBins(c)); C->npf_pt->Sumw2();
				C->nff_pt   = new TH2D(rootname + "_NFF_pt",   "NFF_pt",   getNFPtBins(c), getFPtBins(c), getNFPtBins(c), getFPtBins(c)); C->nff_pt->Sumw2();
				C->nt2pp_pt = new TH2D(rootname + "_NT2PP_pt", "NT2PP_pt", getNFPtBins(c), getFPtBins(c), getNFPtBins(c), getFPtBins(c)); C->nt2pp_pt->Sumw2();
				C->nt2fp_pt = new TH2D(rootname + "_NT2FP_pt", "NT2FP_pt", getNFPtBins(c), getFPtBins(c), getNFPtBins(c), getFPtBins(c)); C->nt2fp_pt->Sumw2();
				C->nt2pf_pt = new TH2D(rootname + "_NT2PF_pt", "NT2PF_pt", getNFPtBins(c), getFPtBins(c), getNFPtBins(c), getFPtBins(c)); C->nt2pf_pt->Sumw2();
				C->nt2ff_pt = new TH2D(rootname + "_NT2FF_pt", "NT2FF_pt", getNFPtBins(c), getFPtBins(c), getNFPtBins(c), getFPtBins(c)); C->nt2ff_pt->Sumw2();

				C->nt11_origin = new TH2D(rootname + "_NT20_Origin",  "NT2Origin",  12, 0, 12, 12, 0, 12);
				C->nt10_origin = new TH2D(rootname + "_NT10_Origin",  "NT1Origin",  12, 0, 12, 12, 0, 12);
				C->nt01_origin = new TH2D(rootname + "_NT01_Origin",  "NT01Origin", 12, 0, 12, 12, 0, 12);
				C->nt00_origin = new TH2D(rootname + "_NT00_Origin",  "NT0Origin",  12, 0, 12, 12, 0, 12);
				label2OriginAxes(C->nt11_origin->GetXaxis(), C->nt11_origin->GetYaxis(), c);
				label2OriginAxes(C->nt10_origin->GetXaxis(), C->nt10_origin->GetYaxis(), c);
				label2OriginAxes(C->nt01_origin->GetXaxis(), C->nt01_origin->GetYaxis(), c);
				label2OriginAxes(C->nt00_origin->GetXaxis(), C->nt00_origin->GetYaxis(), c);					
			}

			// Charge misid truth
			if(c != Muon){
				C->npp_cm_pt   = new TH2D(rootname + "_NPP_CM_pt",   "NPP_CM_pt",   getNFPtBins(c), getFPtBins(c), getNFPtBins(c), getFPtBins(c)); C->npp_cm_pt->Sumw2();
				C->nt2pp_cm_pt = new TH2D(rootname + "_NT2PP_CM_pt", "NT2PP_CM_pt", getNFPtBins(c), getFPtBins(c), getNFPtBins(c), getFPtBins(c)); C->nt2pp_cm_pt->Sumw2();					
			}

			// OS Yields
			if(c == Elec){
				C->nt20_OS_BB_pt = new TH2D(rootname + "_NT20_OS_BB_pt",  "NT20_OS_BB_pt",  getNFPtBins(c), getFPtBins(c), getNFPtBins(c), getFPtBins(c)); C->nt20_OS_BB_pt ->Sumw2();
				C->nt20_OS_EE_pt = new TH2D(rootname + "_NT20_OS_EE_pt",  "NT20_OS_EE_pt",  getNFPtBins(c), getFPtBins(c), getNFPtBins(c), getFPtBins(c)); C->nt20_OS_EE_pt ->Sumw2();
				C->nt20_OS_EB_pt = new TH2D(rootname + "_NT20_OS_EB_pt",  "NT20_OS_EB_pt",  getNFPtBins(c), getFPtBins(c), getNFPtBins(c), getFPtBins(c)); C->nt20_OS_EB_pt ->Sumw2();
				
				C->nt10_OS_BB_pt = new TH2D(rootname + "_NT10_OS_BB_pt",  "NT10_OS_BB_pt",  getNFPtBins(c), getFPtBins(c), getNFPtBins(c), getFPtBins(c)); C->nt10_OS_BB_pt ->Sumw2();
				C->nt10_OS_EE_pt = new TH2D(rootname + "_NT10_OS_EE_pt",  "NT10_OS_EE_pt",  getNFPtBins(c), getFPtBins(c), getNFPtBins(c), getFPtBins(c)); C->nt10_OS_EE_pt ->Sumw2();
				C->nt10_OS_EB_pt = new TH2D(rootname + "_NT10_OS_EB_pt",  "NT10_OS_EB_pt",  getNFPtBins(c), getFPtBins(c), getNFPtBins(c), getFPtBins(c)); C->nt10_OS_EB_pt ->Sumw2();
				
				C->nt01_OS_BB_pt = new TH2D(rootname + "_NT01_OS_BB_pt",  "NT01_OS_BB_pt",  getNFPtBins(c), getFPtBins(c), getNFPtBins(c), getFPtBins(c)); C->nt01_OS_BB_pt ->Sumw2();
				C->nt01_OS_EE_pt = new TH2D(rootname + "_NT01_OS_EE_pt",  "NT01_OS_EE_pt",  getNFPtBins(c), getFPtBins(c), getNFPtBins(c), getFPtBins(c)); C->nt01_OS_EE_pt ->Sumw2();
				C->nt01_OS_EB_pt = new TH2D(rootname + "_NT01_OS_EB_pt",  "NT01_OS_EB_pt",  getNFPtBins(c), getFPtBins(c), getNFPtBins(c), getFPtBins(c)); C->nt01_OS_EB_pt ->Sumw2();
				
				C->nt00_OS_BB_pt = new TH2D(rootname + "_NT00_OS_BB_pt",  "NT00_OS_BB_pt",  getNFPtBins(c), getFPtBins(c), getNFPtBins(c), getFPtBins(c)); C->nt00_OS_BB_pt ->Sumw2();
				C->nt00_OS_EE_pt = new TH2D(rootname + "_NT00_OS_EE_pt",  "NT00_OS_EE_pt",  getNFPtBins(c), getFPtBins(c), getNFPtBins(c), getFPtBins(c)); C->nt00_OS_EE_pt ->Sumw2();
				C->nt00_OS_EB_pt = new TH2D(rootname + "_NT00_OS_EB_pt",  "NT00_OS_EB_pt",  getNFPtBins(c), getFPtBins(c), getNFPtBins(c), getFPtBins(c)); C->nt00_OS_EB_pt ->Sumw2();
			}
			if(c == ElMu){
				C->nt20_OS_BB_pt = new TH2D(rootname + "_NT20_OS_BB_pt",  "NT20_OS_BB_pt",  getNFPtBins(Muon), getFPtBins(Muon), getNFPtBins(Elec), getFPtBins(Elec)); C->nt20_OS_BB_pt ->Sumw2();
				C->nt20_OS_EE_pt = new TH2D(rootname + "_NT20_OS_EE_pt",  "NT20_OS_EE_pt",  getNFPtBins(Muon), getFPtBins(Muon), getNFPtBins(Elec), getFPtBins(Elec)); C->nt20_OS_EE_pt ->Sumw2();
				
				C->nt10_OS_BB_pt = new TH2D(rootname + "_NT10_OS_BB_pt",  "NT10_OS_BB_pt",  getNFPtBins(Muon), getFPtBins(Muon), getNFPtBins(Elec), getFPtBins(Elec)); C->nt10_OS_BB_pt ->Sumw2();
				C->nt10_OS_EE_pt = new TH2D(rootname + "_NT10_OS_EE_pt",  "NT10_OS_EE_pt",  getNFPtBins(Muon), getFPtBins(Muon), getNFPtBins(Elec), getFPtBins(Elec)); C->nt10_OS_EE_pt ->Sumw2();
				
				C->nt01_OS_BB_pt = new TH2D(rootname + "_NT01_OS_BB_pt",  "NT01_OS_BB_pt",  getNFPtBins(Muon), getFPtBins(Muon), getNFPtBins(Elec), getFPtBins(Elec)); C->nt01_OS_BB_pt ->Sumw2();
				C->nt01_OS_EE_pt = new TH2D(rootname + "_NT01_OS_EE_pt",  "NT01_OS_EE_pt",  getNFPtBins(Muon), getFPtBins(Muon), getNFPtBins(Elec), getFPtBins(Elec)); C->nt01_OS_EE_pt ->Sumw2();
				
				C->nt00_OS_BB_pt = new TH2D(rootname + "_NT00_OS_BB_pt",  "NT00_OS_BB_pt",  getNFPtBins(Muon), getFPtBins(Muon), getNFPtBins(Elec), getFPtBins(Elec)); C->nt00_OS_BB_pt ->Sumw2();
				C->nt00_OS_EE_pt = new TH2D(rootname + "_NT00_OS_EE_pt",  "NT00_OS_EE_pt",  getNFPtBins(Muon), getFPtBins(Muon), getNFPtBins(Elec), getFPtBins(Elec)); C->nt00_OS_EE_pt ->Sumw2();
			}

			// Ratios
			if(c != ElMu){
				C->fntight  = new TH2D(rootname + "_fNTight",  "fNTight",  getNFPtBins(c), getFPtBins(c), getNEtaBins(c), getEtaBins(c)); C->fntight ->Sumw2();
				C->fnloose  = new TH2D(rootname + "_fNLoose",  "fNLoose",  getNFPtBins(c), getFPtBins(c), getNEtaBins(c), getEtaBins(c)); C->fnloose ->Sumw2();
				C->pntight  = new TH2D(rootname + "_pNTight",  "pNTight",  getNPPtBins(c), getPPtBins(c), getNEtaBins(c), getEtaBins(c)); C->pntight ->Sumw2();
				C->pnloose  = new TH2D(rootname + "_pNLoose",  "pNLoose",  getNPPtBins(c), getPPtBins(c), getNEtaBins(c), getEtaBins(c)); C->pnloose ->Sumw2();

				C->fntight_nv  = new TH1D(rootname + "_fNTight_nv",  "fNTight_nv", 18, 0., 36.); C->fntight_nv ->Sumw2();
				C->fnloose_nv  = new TH1D(rootname + "_fNLoose_nv",  "fNLoose_nv", 18, 0., 36.); C->fnloose_nv ->Sumw2();
				C->pntight_nv  = new TH1D(rootname + "_pNTight_nv",  "pNTight_nv", 18, 0., 36.); C->pntight_nv ->Sumw2();
				C->pnloose_nv  = new TH1D(rootname + "_pNLoose_nv",  "pNLoose_nv", 18, 0., 36.); C->pnloose_nv ->Sumw2();
				// duplicate for only ttbar ratios
				C->fntight_ttbar  = new TH2D(rootname + "_fNTight_ttbar",  "fNTight_ttbar",  getNFPtBins(c), getFPtBins(c), getNEtaBins(c), getEtaBins(c)); C->fntight_ttbar ->Sumw2();
				C->fnloose_ttbar  = new TH2D(rootname + "_fNLoose_ttbar",  "fNLoose_ttbar",  getNFPtBins(c), getFPtBins(c), getNEtaBins(c), getEtaBins(c)); C->fnloose_ttbar ->Sumw2();
				C->pntight_ttbar  = new TH2D(rootname + "_pNTight_ttbar",  "pNTight_ttbar",  getNPPtBins(c), getPPtBins(c), getNEtaBins(c), getEtaBins(c)); C->pntight_ttbar ->Sumw2();
				C->pnloose_ttbar  = new TH2D(rootname + "_pNLoose_ttbar",  "pNLoose_ttbar",  getNPPtBins(c), getPPtBins(c), getNEtaBins(c), getEtaBins(c)); C->pnloose_ttbar ->Sumw2();
				
				// gen ID
				C->fntight_genID  = new TH1D(rootname + "_fNTight_genID",  "fNTight_genID",  1001, -0.5, 1000.5); C->fntight_genID ->Sumw2();
				C->fnloose_genID  = new TH1D(rootname + "_fNLoose_genID",  "fNLoose_genID",  1001, -0.5, 1000.5); C->fnloose_genID ->Sumw2();
				C->pntight_genID  = new TH1D(rootname + "_pNTight_genID",  "pNTight_genID",  1001, -0.5, 1000.5); C->pntight_genID ->Sumw2();
				C->pnloose_genID  = new TH1D(rootname + "_pNLoose_genID",  "pNLoose_genID",  1001, -0.5, 1000.5); C->pnloose_genID ->Sumw2();
				
				C->fratio_pt  = new TEfficiency(rootname + "_fRatio_pt",  "fRatio_pt",  getNFPtBins(c), getFPtBins(c));
				C->fratio_eta = new TEfficiency(rootname + "_fRatio_eta", "fRatio_eta", getNEtaBins(c), getEtaBins(c));
				C->pratio_pt  = new TEfficiency(rootname + "_pRatio_pt",  "pRatio_pt",  getNPPtBins(c), getPPtBins(c));
				C->pratio_eta = new TEfficiency(rootname + "_pRatio_eta", "pRatio_eta", getNEtaBins(c), getEtaBins(c));
				C->fratio_nv  = new TEfficiency(rootname + "_fRatio_nv",  "fRatio_nv",  18, 0., 36.);
				C->pratio_nv  = new TEfficiency(rootname + "_pRatio_nv",  "pRatio_nv",  18, 0., 36.);

				// C->fratio_pt ->SetUseWeightedEvents();
				// C->fratio_eta->SetUseWeightedEvents();
				// C->pratio_pt ->SetUseWeightedEvents();
				// C->pratio_eta->SetUseWeightedEvents();

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
			
			// Charge Mis id
			if (c == Elec){
			        C->ospairs = new TH2D(rootname + "_ospairs", "ospairs", getNCMidbins(), getCMIdbins(), getNCMidbins(), getCMIdbins()); C->ospairs->Sumw2();
			        C->sspairs = new TH2D(rootname + "_sspairs", "sspairs", getNCMidbins(), getCMIdbins(), getNCMidbins(), getCMIdbins()); C->sspairs->Sumw2();

			}

		}
	}
	// }
}
void SSDLDumper::deleteHistos(Sample *S){
	delete S->cutFlowHisto[Muon];
	delete S->cutFlowHisto[Elec];
	delete S->cutFlowHisto[ElMu];
	
	// Kinematical histos
	for(size_t k = 0; k < gNKinSels; ++k){
		for(size_t j = 0; j < gNKinVars; ++j) delete S->kinplots[k][HighPt].hvar[j];
	}

	// WZ Kinematical histos
	for(size_t k = 0; k < gNKinSels; ++k){
		for(size_t j = 0; j < gNKinVars; ++j) delete S->kinplots_wz[k].hvar[j];
	}

	// Histos for differential yields
	for(size_t k = 0; k < gNCHANNELS; ++k){
		for(size_t j = 0; j < gNDiffVars; ++j){
			delete S->diffyields[k].hnt11[j];
			delete S->diffyields[k].hnt10[j];
			delete S->diffyields[k].hnt01[j];
			delete S->diffyields[k].hnt00[j];
			delete S->diffyields[k].hnpp[j];
			delete S->diffyields[k].hnpf[j];
			delete S->diffyields[k].hnfp[j];
			delete S->diffyields[k].hnff[j];
			if(k == Muon) continue;
			delete S->diffyields[k].hnt2_os_BB[j];
			delete S->diffyields[k].hnt2_os_EE[j];
			if(k == ElMu) continue;
			delete S->diffyields[k].hnt2_os_EB[j];
		}
	}

	// id histos for electrons
	for(size_t j = 0; j < gNSels; ++j){
		delete S->idplots.hhoe   [j];
		delete S->idplots.hsiesie[j];
		delete S->idplots.hdeta  [j];
		delete S->idplots.hdphi  [j];
		delete S->idplots.hmvaid [j];
		delete S->idplots.hmedwp [j];
	}

	for(size_t l = 0; l < 2; ++l){
		// Isolation histos
		for(size_t j = 0; j < gNSels; ++j){
			delete S->isoplots[l].hiso[j];
			for(int k = 0; k < gNMuFPtBins; ++k){
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

	for(regIt = gRegions.begin(); regIt != gRegions.end() ; regIt++){
		int r = gRegion[(*regIt)->sname];
		Region *R = &S->region[r][HighPt];
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
				
				delete C->nt10_OS_BB_pt;
				delete C->nt10_OS_EE_pt;
				delete C->nt10_OS_EB_pt;
				
				delete C->nt01_OS_BB_pt;
				delete C->nt01_OS_EE_pt;
				delete C->nt01_OS_EB_pt;
				
				delete C->nt00_OS_BB_pt;
				delete C->nt00_OS_EE_pt;
				delete C->nt00_OS_EB_pt;
			}
			if(c == ElMu){
				delete C->nt20_OS_BB_pt;
				delete C->nt20_OS_EE_pt;
				
				delete C->nt10_OS_BB_pt;
				delete C->nt10_OS_EE_pt;
				
				delete C->nt01_OS_BB_pt;
				delete C->nt01_OS_EE_pt;
				
				delete C->nt00_OS_BB_pt;
				delete C->nt00_OS_EE_pt;
			}

			// Ratios
			if(c != ElMu){
				delete C->fntight;
				delete C->fnloose;
				delete C->pntight;
				delete C->pnloose;
				delete C->fntight_nv;
				delete C->fnloose_nv;
				delete C->pntight_nv;
				delete C->pnloose_nv;
				
				delete C->fratio_nv;
				delete C->pratio_nv;
				// duplicate for ttbar only ratios
				delete C->fntight_ttbar;
				delete C->fnloose_ttbar;
				delete C->pntight_ttbar;
				delete C->pnloose_ttbar;
				
				// gen ID
				delete C->fntight_genID;
				delete C->fnloose_genID;
				delete C->pntight_genID;
				delete C->pnloose_genID;
				
				delete C->fratio_pt;
				delete C->pratio_pt;
				delete C->fratio_eta;
				delete C->pratio_eta;

				if(S->datamc > 0){
					delete C->sst_origin;
					delete C->ssl_origin;
					delete C->zt_origin;
					delete C->zl_origin;
				}
			}
			//  Charge Miss-ID
			if (c == Elec){
			       delete C->ospairs;
			       delete C->sspairs;
			}
			
		}
	}
	// }
}
void SSDLDumper::writeHistos(Sample *S, TFile *pFile){
	pFile->cd();
	TDirectory* sdir = Util::FindOrCreate(S->sname, pFile);
	sdir->cd();

	TString temp;
	TDirectory *rdir;

	// Event counter
	S->evcount->Write(S->evcount->GetName(), TObject::kWriteDelete);

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
			S->diffyields[k].hnpp[j] ->Write(S->diffyields[k].hnpp[j] ->GetName(), TObject::kWriteDelete);
			S->diffyields[k].hnpf[j] ->Write(S->diffyields[k].hnpf[j] ->GetName(), TObject::kWriteDelete);
			S->diffyields[k].hnfp[j] ->Write(S->diffyields[k].hnfp[j] ->GetName(), TObject::kWriteDelete);
			S->diffyields[k].hnff[j] ->Write(S->diffyields[k].hnff[j] ->GetName(), TObject::kWriteDelete);
			if(k == Muon) continue;
			S->diffyields[k].hnt2_os_BB[j]->Write(S->diffyields[k].hnt2_os_BB[j]->GetName(), TObject::kWriteDelete);
			S->diffyields[k].hnt2_os_EE[j]->Write(S->diffyields[k].hnt2_os_EE[j]->GetName(), TObject::kWriteDelete);
			if(k == ElMu) continue;
			S->diffyields[k].hnt2_os_EB[j]->Write(S->diffyields[k].hnt2_os_EB[j]->GetName(), TObject::kWriteDelete);
		}
	}	

	// Kinematic histos
	temp = S->sname + "/KinPlots/";
	rdir = Util::FindOrCreate(temp, pFile);
	rdir->cd();
	for(size_t k = 0; k < gNKinSels; ++k){
		KinPlots *kp = &S->kinplots[k][HighPt];
		for(size_t j = 0; j < gNKinVars; ++j) kp->hvar[j]->Write(kp->hvar[j]->GetName(), TObject::kWriteDelete);
	}


	// WZ Kinematic histos
	temp = S->sname + "/WZValidation/";
	rdir = Util::FindOrCreate(temp, pFile);
	rdir->cd();
	for(size_t k = 0; k < gNKinSels; ++k){
		KinPlots *kp = &S->kinplots_wz[k];
		for(size_t j = 0; j < gNKinVars; ++j) kp->hvar[j]->Write(kp->hvar[j]->GetName(), TObject::kWriteDelete);
	}


	// Id histos for electrons
	temp = S->sname + "/IdPlots/";
	rdir = Util::FindOrCreate(temp, pFile);
	rdir->cd();
	IdPlots *idp = &S->idplots;
	for(size_t j = 0; j < gNSels; ++j){
		idp->hhoe   [j]->Write(idp->hhoe   [j]->GetName(), TObject::kWriteDelete);
		idp->hsiesie[j]->Write(idp->hsiesie[j]->GetName(), TObject::kWriteDelete);
		idp->hdeta  [j]->Write(idp->hdeta  [j]->GetName(), TObject::kWriteDelete);
		idp->hdphi  [j]->Write(idp->hdphi  [j]->GetName(), TObject::kWriteDelete);
		idp->hmvaid [j]->Write(idp->hmvaid [j]->GetName(), TObject::kWriteDelete);
		idp->hmedwp [j]->Write(idp->hmedwp [j]->GetName(), TObject::kWriteDelete);
	}

	// Isolation histos
	temp = S->sname + "/IsoPlots/";
	rdir = Util::FindOrCreate(temp, pFile);
	rdir->cd();
	for(size_t l = 0; l < 2; ++l){
		IsoPlots *ip = &S->isoplots[l];
		for(size_t j = 0; j < gNSels; ++j){
			ip->hiso[j]->Write(ip->hiso[j]->GetName(), TObject::kWriteDelete);
			for(int k = 0; k < gNMuFPtBins; ++k) ip->hiso_pt[j][k]->Write(ip->hiso_pt[j][k]->GetName(), TObject::kWriteDelete);
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
	for(regIt = gRegions.begin(); regIt != gRegions.end(); regIt++){
		int r = gRegion[(*regIt)->sname];
		Region *R = &S->region[r][HighPt];
		TString temp = S->sname + "/" + gRegions[r]->sname;
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
				C->nt10_OS_BB_pt->Write(C->nt10_OS_BB_pt->GetName(), TObject::kWriteDelete);
				C->nt10_OS_EE_pt->Write(C->nt10_OS_EE_pt->GetName(), TObject::kWriteDelete);
				C->nt01_OS_BB_pt->Write(C->nt01_OS_BB_pt->GetName(), TObject::kWriteDelete);
				C->nt01_OS_EE_pt->Write(C->nt01_OS_EE_pt->GetName(), TObject::kWriteDelete);
				C->nt00_OS_BB_pt->Write(C->nt00_OS_BB_pt->GetName(), TObject::kWriteDelete);
				C->nt00_OS_EE_pt->Write(C->nt00_OS_EE_pt->GetName(), TObject::kWriteDelete);
				if(ch == Elec) {
					C->nt20_OS_EB_pt->Write(C->nt20_OS_EB_pt->GetName(), TObject::kWriteDelete);
					C->nt10_OS_EB_pt->Write(C->nt10_OS_EB_pt->GetName(), TObject::kWriteDelete);
					C->nt01_OS_EB_pt->Write(C->nt01_OS_EB_pt->GetName(), TObject::kWriteDelete);
					C->nt00_OS_EB_pt->Write(C->nt00_OS_EB_pt->GetName(), TObject::kWriteDelete);
				}
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
				C->fntight_nv ->Write(C->fntight_nv ->GetName(), TObject::kWriteDelete);
				C->fnloose_nv ->Write(C->fnloose_nv ->GetName(), TObject::kWriteDelete);
				C->pntight_nv ->Write(C->pntight_nv ->GetName(), TObject::kWriteDelete);
				C->pnloose_nv ->Write(C->pnloose_nv ->GetName(), TObject::kWriteDelete);

				// duplicate for ttbar only ratios
				C->fntight_ttbar    ->Write(C->fntight_ttbar    ->GetName(), TObject::kWriteDelete);
				C->fnloose_ttbar    ->Write(C->fnloose_ttbar    ->GetName(), TObject::kWriteDelete);
				C->pntight_ttbar    ->Write(C->pntight_ttbar    ->GetName(), TObject::kWriteDelete);
				C->pnloose_ttbar    ->Write(C->pnloose_ttbar    ->GetName(), TObject::kWriteDelete);
				
				// gen ID
				C->fntight_genID    ->Write(C->fntight_genID    ->GetName(), TObject::kWriteDelete);
				C->fnloose_genID    ->Write(C->fnloose_genID    ->GetName(), TObject::kWriteDelete);
				C->pntight_genID    ->Write(C->pntight_genID    ->GetName(), TObject::kWriteDelete);
				C->pnloose_genID    ->Write(C->pnloose_genID    ->GetName(), TObject::kWriteDelete);

				C->fratio_pt ->Write(C->fratio_pt ->GetName(), TObject::kWriteDelete);
				C->pratio_pt ->Write(C->pratio_pt ->GetName(), TObject::kWriteDelete);
				C->fratio_eta->Write(C->fratio_eta->GetName(), TObject::kWriteDelete);
				C->pratio_eta->Write(C->pratio_eta->GetName(), TObject::kWriteDelete);
				C->fratio_nv ->Write(C->fratio_nv ->GetName(), TObject::kWriteDelete);
				C->pratio_nv ->Write(C->pratio_nv ->GetName(), TObject::kWriteDelete);

				if(S->datamc > 0){
					C->sst_origin ->Write(C->sst_origin ->GetName(), TObject::kWriteDelete);
					C->ssl_origin ->Write(C->ssl_origin ->GetName(), TObject::kWriteDelete);
					C->zt_origin  ->Write(C->zt_origin  ->GetName(), TObject::kWriteDelete);
					C->zl_origin  ->Write(C->zl_origin  ->GetName(), TObject::kWriteDelete);						
				}
			}
			if (ch == Elec){
                          	C->ospairs->Write(C->ospairs->GetName(), TObject::kWriteDelete);
                          	C->sspairs->Write(C->sspairs->GetName(), TObject::kWriteDelete);
			}
		}
	}
}
void SSDLDumper::writeSigGraphs(Sample *S, gChannel chan, TFile *pFile){
	TString channame = "MM";
	if(chan == Elec) channame = "EE";
	if(chan == ElMu)      channame = "EM";
	vector<float> ht;
	vector<float> met;
	if(chan == Muon){
		ht  = fSigEv_HI_MM_HT;
		met = fSigEv_HI_MM_MET;
	}
	if(chan == Elec){
		ht  = fSigEv_HI_EE_HT;
		met = fSigEv_HI_EE_MET;
	}
	if(chan == ElMu){
		ht  = fSigEv_HI_EM_HT;
		met = fSigEv_HI_EM_MET;
	}
	
	const int nsig = ht.size();
	float ht_a [nsig];
	float met_a[nsig];
	for(size_t i = 0; i < ht.size(); ++i){
		ht_a[i] = ht[i];
		met_a[i] = met[i];
	}
	
	TGraph *sigevents = new TGraph(nsig, ht_a, met_a);
	sigevents->SetName(Form("%s_%s_SigEvents", S->sname.Data(), channame.Data()));

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
		if(fVerbose > 2) cout << "Reading histos for " << S->sname << endl;

		// Event count histo
		getObjectSafe(pFile, S->sname + "/" + S->sname + "_EventCount", S->evcount);
		S->ngen = S->evcount->GetEntries();

		// Cut flow histos
		getObjectSafe(pFile, S->sname + "/MMCutFlow", S->cutFlowHisto[Muon]);
		getObjectSafe(pFile, S->sname + "/EECutFlow", S->cutFlowHisto[Elec]);
		getObjectSafe(pFile, S->sname + "/EMCutFlow", S->cutFlowHisto[ElMu]);

		// Histos for differential yields
		for(size_t k = 0; k < gNCHANNELS; ++k){
			TString name;
			for(size_t j = 0; j < gNDiffVars; ++j){
				getname = Form("%s_%s_NT11_%s", S->sname.Data(), DiffPredYields::var_name[j].Data(), gChanLabel[k].Data());
				getObjectSafe(pFile, S->sname + "/DiffYields/" + getname, S->diffyields[k].hnt11[j]);
				getname = Form("%s_%s_NT10_%s", S->sname.Data(), DiffPredYields::var_name[j].Data(), gChanLabel[k].Data());
				getObjectSafe(pFile, S->sname + "/DiffYields/" + getname, S->diffyields[k].hnt10[j]);
				getname = Form("%s_%s_NT01_%s", S->sname.Data(), DiffPredYields::var_name[j].Data(), gChanLabel[k].Data());
				getObjectSafe(pFile, S->sname + "/DiffYields/" + getname, S->diffyields[k].hnt01[j]);
				getname = Form("%s_%s_NT00_%s", S->sname.Data(), DiffPredYields::var_name[j].Data(), gChanLabel[k].Data());
				getObjectSafe(pFile, S->sname + "/DiffYields/" + getname, S->diffyields[k].hnt00[j]);
				getname = Form("%s_%s_NPP_%s", S->sname.Data(), DiffPredYields::var_name[j].Data(), gChanLabel[k].Data());
				getObjectSafe(pFile, S->sname + "/DiffYields/" + getname, S->diffyields[k].hnpp[j]);
				getname = Form("%s_%s_NPF_%s", S->sname.Data(), DiffPredYields::var_name[j].Data(), gChanLabel[k].Data());
				getObjectSafe(pFile, S->sname + "/DiffYields/" + getname, S->diffyields[k].hnpf[j]);
				getname = Form("%s_%s_NFP_%s", S->sname.Data(), DiffPredYields::var_name[j].Data(), gChanLabel[k].Data());
				getObjectSafe(pFile, S->sname + "/DiffYields/" + getname, S->diffyields[k].hnfp[j]);
				getname = Form("%s_%s_NFF_%s", S->sname.Data(), DiffPredYields::var_name[j].Data(), gChanLabel[k].Data());
				getObjectSafe(pFile, S->sname + "/DiffYields/" + getname, S->diffyields[k].hnff[j]);
				if(k == Muon) continue;
				getname = Form("%s_%s_NT11_OS_BB_%s", S->sname.Data(), DiffPredYields::var_name[j].Data(), gChanLabel[k].Data());
				getObjectSafe(pFile, S->sname + "/DiffYields/" + getname, S->diffyields[k].hnt2_os_BB[j]);
				getname = Form("%s_%s_NT11_OS_EE_%s", S->sname.Data(), DiffPredYields::var_name[j].Data(), gChanLabel[k].Data());
				getObjectSafe(pFile, S->sname + "/DiffYields/" + getname, S->diffyields[k].hnt2_os_EE[j]);
				if(k == ElMu) continue;
				getname = Form("%s_%s_NT11_OS_EB_%s", S->sname.Data(), DiffPredYields::var_name[j].Data(), gChanLabel[k].Data());
				getObjectSafe(pFile, S->sname + "/DiffYields/" + getname, S->diffyields[k].hnt2_os_EB[j]);
			}
		}	


		// Kinematic histos
		for(size_t k = 0; k < gNKinSels; ++k){
			KinPlots *kp = &S->kinplots[k][HighPt];
			getname = Form("%s_%s_HTvsMET", S->sname.Data(), gKinSelNames[k].Data());
			for(size_t j = 0; j < gNKinVars; ++j){
				getname = Form("%s_%s_%s", S->sname.Data(), gKinSelNames[k].Data(), KinPlots::var_name[j].Data());
				getObjectSafe(pFile, S->sname + "/KinPlots/" + getname, kp->hvar[j]);
				kp->hvar[j]->SetFillColor(S->color);
			}
		}


		// WZ Kinematic histos
		for(size_t k = 0; k < gNKinSels; ++k){
			KinPlots *kp = &S->kinplots_wz[k];
			for(size_t j = 0; j < gNKinVars; ++j){
				getname = Form("%s_%s_%s", S->sname.Data(), gKinSelNames[k].Data(), KinPlots::var_name[j].Data());
				getObjectSafe(pFile, S->sname + "/WZValidation/" + getname, kp->hvar[j]);
				kp->hvar[j]->SetFillColor(S->color);
			}
		}


		// Id histos for electrons only
		IdPlots *idp = &S->idplots;
		for(size_t j = 0; j < gNSels; ++j){
			// hoe
			getname = Form("%s_%s_%shoe", S->sname.Data(), IdPlots::sel_name[j].Data(), gEMULabel[1].Data());
			getObjectSafe(pFile, S->sname + "/IdPlots/" + getname, idp->hhoe[j]);
			idp->hhoe[j]->SetFillColor(S->color);
			// sigma ieta ieta
			getname = Form("%s_%s_%ssiesie", S->sname.Data(), IdPlots::sel_name[j].Data(), gEMULabel[1].Data());
			getObjectSafe(pFile, S->sname + "/IdPlots/" + getname, idp->hsiesie[j]);
			idp->hsiesie[j]->SetFillColor(S->color);
			// delta eta
			getname = Form("%s_%s_%sdeta", S->sname.Data(), IdPlots::sel_name[j].Data(), gEMULabel[1].Data());
			getObjectSafe(pFile, S->sname + "/IdPlots/" + getname, idp->hdeta[j]);
			idp->hdeta[j]->SetFillColor(S->color);
			// delta phi
			getname = Form("%s_%s_%sdphi", S->sname.Data(), IdPlots::sel_name[j].Data(), gEMULabel[1].Data());
			getObjectSafe(pFile, S->sname + "/IdPlots/" + getname, idp->hdphi[j]);
			idp->hdphi[j]->SetFillColor(S->color);
			// mva id
			getname = Form("%s_%s_%smvaid", S->sname.Data(), IdPlots::sel_name[j].Data(), gEMULabel[1].Data());
			getObjectSafe(pFile, S->sname + "/IdPlots/" + getname, idp->hmvaid[j]);
			idp->hmvaid[j]->SetFillColor(S->color);
			// medium ID WP pass
			getname = Form("%s_%s_%smedwp", S->sname.Data(), IdPlots::sel_name[j].Data(), gEMULabel[1].Data());
			getObjectSafe(pFile, S->sname + "/IdPlots/" + getname, idp->hmedwp[j]);
			idp->hmedwp[j]->SetFillColor(S->color);
		}

		for(size_t lep = 0; lep < 2; ++lep){ // e-mu loop
			// Isolation histos
			IsoPlots *ip = &S->isoplots[lep];
			for(size_t j = 0; j < gNSels; ++j){
				getname = Form("%s_%s_%siso", S->sname.Data(), IsoPlots::sel_name[j].Data(), gEMULabel[lep].Data());
				getObjectSafe(pFile, S->sname + "/IsoPlots/" + getname, ip->hiso[j]);
				ip->hiso[j]->SetFillColor(S->color);
				for(int k = 0; k < gNMuFPtBins; ++k){
					getname = Form("%s_%s_%siso_pt%d", S->sname.Data(), IsoPlots::sel_name[j].Data(), gEMULabel[lep].Data(), k);
					getObjectSafe(pFile, S->sname + "/IsoPlots/" + getname, ip->hiso_pt[j][k]);
					ip->hiso_pt[j][k]->SetFillColor(S->color);
				}
				for(int k = 0; k < gNNVrtxBins; ++k){
					getname = Form("%s_%s_%siso_nv%d", S->sname.Data(), IsoPlots::sel_name[j].Data(), gEMULabel[lep].Data(), k);
					getObjectSafe(pFile, S->sname + "/IsoPlots/" + getname, ip->hiso_nv[j][k]);
					ip->hiso_nv[j][k]->SetFillColor(S->color);
				}
			}

			// Ratio histos
			FRatioPlots *rp = &S->ratioplots[lep];
			for(size_t j = 0; j < gNRatioVars; ++j){
				getname = Form("%s_%s_ntight_%s", S->sname.Data(), gEMULabel[lep].Data(), FRatioPlots::var_name[j].Data());
				getObjectSafe(pFile, S->sname + "/FRatioPlots/" + getname, rp->ntight[j]);
				getname = Form("%s_%s_nloose_%s", S->sname.Data(), gEMULabel[lep].Data(), FRatioPlots::var_name[j].Data());
				getObjectSafe(pFile, S->sname + "/FRatioPlots/" + getname, rp->nloose[j]);
			}
		}

		// Yields
		for(size_t HighPt = 0; HighPt < 2; ++HighPt){
			for(regIt = gRegions.begin(); regIt != gRegions.end(); regIt++){ // Loop over regions
				int r = gRegion[(*regIt)->sname];
				Region *R = &S->region[r][HighPt];
				for(gChannel ch = channels_begin; ch < gNCHANNELS; ch=gChannel(ch+1)){ // Loop over channels, mumu, emu, ee
					Channel *C;
					if(ch == Muon) C = &R->mm;
					if(ch == Elec) C = &R->ee;
					if(ch == ElMu) C = &R->em;
					TString root = S->sname +"/"+ gRegions[r]->sname +"/"+ S->sname +"_"+ gRegions[r]->sname +"_"+ gChanLabel[ch];
					getObjectSafe(pFile, root + "_NT20_pt",  C->nt20_pt );
					getObjectSafe(pFile, root + "_NT10_pt" , C->nt10_pt );
					getObjectSafe(pFile, root + "_NT01_pt" , C->nt01_pt );
					getObjectSafe(pFile, root + "_NT00_pt" , C->nt00_pt );
					getObjectSafe(pFile, root + "_NT20_eta", C->nt20_eta);
					getObjectSafe(pFile, root + "_NT10_eta", C->nt10_eta);
					getObjectSafe(pFile, root + "_NT01_eta", C->nt01_eta);
					getObjectSafe(pFile, root + "_NT00_eta", C->nt00_eta);
					if(S->datamc > 0){
						getObjectSafe(pFile, root + "_NPP_pt"     , C->npp_pt     );
						getObjectSafe(pFile, root + "_NFP_pt"     , C->nfp_pt     );
						getObjectSafe(pFile, root + "_NPF_pt"     , C->npf_pt     );
						getObjectSafe(pFile, root + "_NFF_pt"     , C->nff_pt     );
						getObjectSafe(pFile, root + "_NT2PP_pt"   , C->nt2pp_pt   );
						getObjectSafe(pFile, root + "_NT2FP_pt"   , C->nt2fp_pt   );
						getObjectSafe(pFile, root + "_NT2PF_pt"   , C->nt2pf_pt   );
						getObjectSafe(pFile, root + "_NT2FF_pt"   , C->nt2ff_pt   );
						getObjectSafe(pFile, root + "_NT20_Origin", C->nt11_origin);
						getObjectSafe(pFile, root + "_NT10_Origin", C->nt10_origin);
						getObjectSafe(pFile, root + "_NT01_Origin", C->nt01_origin);
						getObjectSafe(pFile, root + "_NT00_Origin", C->nt00_origin);
						
						if(ch != Muon){
							getObjectSafe(pFile, root + "_NPP_CM_pt"  , C->npp_cm_pt  );
							getObjectSafe(pFile, root + "_NT2PP_CM_pt", C->nt2pp_cm_pt);
						}
					}
					if(ch == Elec || ch == ElMu){
						getObjectSafe(pFile, root + "_NT20_OS_BB_pt", C->nt20_OS_BB_pt);
						getObjectSafe(pFile, root + "_NT20_OS_EE_pt", C->nt20_OS_EE_pt);
						getObjectSafe(pFile, root + "_NT10_OS_BB_pt", C->nt10_OS_BB_pt);
						getObjectSafe(pFile, root + "_NT10_OS_EE_pt", C->nt10_OS_EE_pt);
						getObjectSafe(pFile, root + "_NT01_OS_BB_pt", C->nt01_OS_BB_pt);
						getObjectSafe(pFile, root + "_NT01_OS_EE_pt", C->nt01_OS_EE_pt);
						getObjectSafe(pFile, root + "_NT00_OS_BB_pt", C->nt00_OS_BB_pt);
						getObjectSafe(pFile, root + "_NT00_OS_EE_pt", C->nt00_OS_EE_pt);
						if(ch == Elec) {
							getObjectSafe(pFile, root + "_NT20_OS_EB_pt", C->nt20_OS_EB_pt);
							getObjectSafe(pFile, root + "_NT10_OS_EB_pt", C->nt10_OS_EB_pt);
							getObjectSafe(pFile, root + "_NT01_OS_EB_pt", C->nt01_OS_EB_pt);
							getObjectSafe(pFile, root + "_NT00_OS_EB_pt", C->nt00_OS_EB_pt);
						}
					}

					if(ch != ElMu){
						getObjectSafe(pFile, root + "_fNTight", C->fntight);
						getObjectSafe(pFile, root + "_fNLoose", C->fnloose);
						getObjectSafe(pFile, root + "_pNTight", C->pntight);
						getObjectSafe(pFile, root + "_pNLoose", C->pnloose);

						getObjectSafe(pFile, root + "_fNTight_nv", C->fntight_nv);
						getObjectSafe(pFile, root + "_fNLoose_nv", C->fnloose_nv);
						getObjectSafe(pFile, root + "_pNTight_nv", C->pntight_nv);
						getObjectSafe(pFile, root + "_pNLoose_nv", C->pnloose_nv);
						// duplicate for ttbar only ratios
						getObjectSafe(pFile, root + "_fNTight_ttbar", C->fntight_ttbar);
						getObjectSafe(pFile, root + "_fNLoose_ttbar", C->fnloose_ttbar);
						getObjectSafe(pFile, root + "_pNTight_ttbar", C->pntight_ttbar);
						getObjectSafe(pFile, root + "_pNLoose_ttbar", C->pnloose_ttbar);
						
						// gen ID
						getObjectSafe(pFile, root + "_fNTight_genID", C->fntight_genID);
						getObjectSafe(pFile, root + "_fNLoose_genID", C->fnloose_genID);
						getObjectSafe(pFile, root + "_pNTight_genID", C->pntight_genID);
						getObjectSafe(pFile, root + "_pNLoose_genID", C->pnloose_genID);

						getObjectSafe(pFile, root + "_fRatio_pt",  C->fratio_pt);
						getObjectSafe(pFile, root + "_pRatio_pt",  C->pratio_pt);
						getObjectSafe(pFile, root + "_fRatio_eta", C->fratio_eta);
						getObjectSafe(pFile, root + "_pRatio_eta", C->pratio_eta);
						getObjectSafe(pFile, root + "_fRatio_nv", C->fratio_nv);
						getObjectSafe(pFile, root + "_pRatio_nv", C->pratio_nv);

						if(S->datamc > 0){
							getObjectSafe(pFile, root + "_fTOrigin", C->sst_origin);
							getObjectSafe(pFile, root + "_fLOrigin", C->ssl_origin);
							getObjectSafe(pFile, root + "_pTOrigin", C->zt_origin );
							getObjectSafe(pFile, root + "_pLOrigin", C->zl_origin );
						}
					}
					if (ch == Elec){
					       getObjectSafe(pFile, root + "_ospairs", C->ospairs);
					       getObjectSafe(pFile, root + "_sspairs", C->sspairs);

					}
				}
				if(HighPt == HighPt){
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
		for(size_t HighPt = 0; HighPt < 2; ++HighPt){
			for(gChannel ch = channels_begin; ch < gNCHANNELS; ch=gChannel(ch+1)){ // Loop over channels, mumu, emu, ee
				if(ch == Muon)     channame = "MM";
				if(ch == Elec) channame = "EE";
				if(ch == ElMu)      channame = "EM";
				getname = Form("%s/SigGraphs/%s_%s_SigEvents", S->sname.Data(), S->sname.Data(), channame.Data());
				getObjectSafe(pFile, getname, S->sigevents[ch][HighPt]);
				S->sigevents[ch][HighPt]->SetMarkerColor(color[ch]);
				S->sigevents[ch][HighPt]->SetMarkerStyle(style[ch]);
				S->sigevents[ch][HighPt]->SetMarkerSize(size);
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
	if(mtype == 1)                                return 2; // W/Z skipped in madgraph event (WGstar/TTV samples)
	if(mid == 24)                                 return 2; // W
	if(mid == 23)                                 return 3; // Z
	if(mtype == 2)                                return 4; // tau
	if(mtype == 11 || mtype == 12 || mtype == 18) return 5; // light hadrons
	if(mtype == 13 || mtype == 19)                return 6; // strange hadrons
	if(mtype == 14 || mtype == 16 || mtype == 20) return 7; // charmed hadrons
	if(mtype == 15 || mtype == 17 || mtype == 21) return 8; // bottom hadrons
	if(mtype == 91 || mtype == 92)                return 9; // pythia strings
	return 12;                                              // uid
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
		return 12;                                            // uid
	}
	if(mtype == 1)                                 return 4;  // W/Z skipped in madgraph event (WGstar/TTV samples)
	if(mid == 24)                                  return 4;  // W
	if(mid == 23)                                  return 5;  // Z
	if(mtype == 2)                                 return 6;  // tau

	if(mtype == 11 || mtype == 12 || mtype == 18)  return 7;  // light hadrons
	if(mtype == 13 || mtype == 19)                 return 8;  // strange hadrons
	if(mtype == 14 || mtype == 16 || mtype == 20)  return 9;  // charmed hadrons
	if(mtype == 15 || mtype == 17 || mtype == 21)  return 10; // bottom hadrons
	if(mtype == 91 || mtype == 92)                 return 11; // pythia strings
	return 12;                                                // uid
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
		case 12: return "Unidentified";
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
		case 12: return "Unidentified";
		default: return "?";
	}
}
void SSDLDumper::labelOriginAxis(TAxis *axis, gChannel chan){
	if(chan == Muon){
		axis->SetTitle("#mu Origin");
		for(size_t i = 1; i <= 12; ++i){
			axis->SetBinLabel(i, muBinToLabel(i));
		}		
	}
	if(chan == Elec){
		axis->SetTitle("e Origin");
		for(size_t i = 1; i <= 12; ++i){
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
	return ( doubleMuTrigger() );
}
bool SSDLDumper::elelSignalTrigger(){
	return ( doubleElTrigger() );
}
bool SSDLDumper::elmuSignalTrigger(){
	return ( eMuTrigger() );
}

//____________________________________________________________________________
// The following triggers will always return true for MC samples, and they
// will always return false if they are called in the 'wrong' dataset
bool SSDLDumper::singleMuTrigger(){
	// Pretend MC samples always fire trigger
	if(fSample->datamc > 0) return true;
	return ( HLT_MU17 > 0  );
	// marc return ( (HLT_MU8 > 0 || HLT_MU17 > 0 ) );
}
float SSDLDumper::singleMuPrescale(){
	// Pretend MC samples have prescale 1.
	if(fSample->datamc > 0) return 1.;
	// Get the prescale factor for whichever of these triggers fired
	// Only correct if they are mutually exclusive!
	// marc if(HLT_MU8_PS > 0) return HLT_MU8_PS;
	if(HLT_MU17_PS > 0) return HLT_MU17_PS;
	return 1;
}
bool SSDLDumper::singleElTrigger(){
  // Pretend MC samples always fire trigger
	if(fSample->datamc > 0) return true;
	return (HLT_ELE17_JET30_TIGHT > 0);
	// marc return ((HLT_ELE8_TIGHT > 0) || (HLT_ELE8_JET30_TIGHT > 0) || (HLT_ELE17_JET30_TIGHT > 0) || (HLT_ELE17_TIGHT > 0));
}
float SSDLDumper::singleElPrescale(){
	// Pretend MC samples have prescale 1.
  if(fSample->datamc > 0) return 1.;
	// Get the prescale factor for whichever of these triggers fired
	// Only correct if they are mutually exclusive!
	if( HLT_ELE17_JET30_TIGHT_PS > 0 ) return HLT_ELE17_JET30_TIGHT_PS;
	// marc if( HLT_ELE8_JET30_TIGHT_PS > 0 )  return HLT_ELE8_JET30_TIGHT_PS;
	// marc if( HLT_ELE17_TIGHT_PS > 0 )       return HLT_ELE17_TIGHT_PS;
	// marc if( HLT_ELE8_TIGHT_PS > 0 ) return HLT_ELE8_TIGHT_PS;
	return 1.;
}

bool SSDLDumper::doubleMuTrigger(){
	// Pretend MC samples always fire trigger
	if(fSample->datamc > 0) return true;
	// Only look for mumu events
	if(fSample->chansel != -1 && fSample->chansel != 0) return false;
	return ( (HLT_MU17_MU8  > 0) );
}
bool SSDLDumper::doubleElTrigger(){
	// Pretend MC samples always fire trigger
	if(fSample->datamc > 0) return true;
	// Only look for elel events
	if(fSample->chansel != -1 && fSample->chansel != 1) return false;

	return ( (HLT_ELE17_ELE8_TIGHT > 0) );
}
bool SSDLDumper::eMuTrigger(){
	// Pretend MC samples always fire trigger
	if(fSample->datamc > 0) return true;
	// Only look for emu events
	if(fSample->chansel != -1 && fSample->chansel != 2) return false;
	return ( (HLT_MU17_ELE8_TIGHT > 0) ||
	         (HLT_MU8_ELE17_TIGHT > 0) );
}
//////////////////////////////////////////////////////////////////////////////
// Helper functions:
//____________________________________________________________________________
void SSDLDumper::scaleBTags(Sample *S, int flag){
	// for now supports only CSVM b-tagger. can be extended if need be
	if(S->datamc == 0) return; // don't smear data
	for(size_t i = 0; i < NJets; ++i){
		if(isGoodJet(i) == false) continue;
		bool is_tagged_lse = JetCSVBTag[i] > 0.244;
		bool is_tagged_med = JetCSVBTag[i] > 0.679;
		float random = fRand3->Uniform(0,1);
		string meanminmax = "mean";
		if(flag == 1) meanminmax = "max";
		if(flag == 2) meanminmax = "min";

		bool newTag = fBTagSF->modifyBTagsWithSF(is_tagged_med, JetPt[i], JetEta[i], JetPartonID[i], meanminmax, random); // WARNING: change this to partonflavor once the minitrees are ready!!!
		if(!newTag) JetCSVBTag[i] = 0.1; // not tagged
		if( newTag) JetCSVBTag[i] = 1.0; // tagged
	}
}

// void SSDLDumper::scaleBTags(Sample *S, int flag){
// 	if(S->datamc == 0) return; // don't smear data
// 	for(size_t i = 0; i < NJets; ++i){
// 		// if(isGoodJet(i) == false) continue;
// 		bool is_tagged_lse = JetCSVBTag[i] > 0.244;
// 		bool is_tagged_med = JetCSVBTag[i] > 0.679;
// 		int pdgid = 0;
// 		if(JetPartonID[i] == 0 || JetPartonID[i] == -2) pdgid = -999;
// 		else if(JetPartonID[i] == -1)                   pdgid = 0;
// 		else pdgid = JetPartonID[i];
// 		
// 		string meanminmax = "mean";
// 		if(flag == 1) meanminmax = "max";
// 		if(flag == 2) meanminmax = "min";
// 		fBTagSFUtil->modifyBTagsWithSF_fast(is_tagged_lse, is_tagged_med, JetPt[i], JetEta[i], pdgid, meanminmax);
// 		if(!is_tagged_lse && !is_tagged_med) JetCSVBTag[i] = 0.1; // not tagged
// 		if( is_tagged_lse && !is_tagged_med) JetCSVBTag[i] = 0.5; // loose tagged
// 		if( is_tagged_lse &&  is_tagged_med) JetCSVBTag[i] = 1.0; // medium tagged
// 	}
// }
void SSDLDumper::saveBTags(){
	// Saves the current tagger values for all jets in a vector
	fSaved_Tags.clear();
	for(size_t i = 0; i < NJets; ++i) fSaved_Tags.push_back(JetCSVBTag[i]);
}
void SSDLDumper::resetBTags(){
	// Resets the b tags for all jets to saved values
	if(fSaved_Tags.size() < 1) return;
	for(size_t i = 0; i < NJets; ++i) JetCSVBTag[i] = fSaved_Tags[i];
}
void SSDLDumper::scaleMET(Sample *S, int flag){
	// first try on MET uncertainty
	if(S->datamc == 0) return; // don't scale data
	TLorentzVector umet, jets, leps, tmp;
	tmp.SetPtEtaPhiM(getMET(), 0., getMETPhi(), 0.);
	umet += tmp; // add met
	tmp.SetPtEtaPhiM(0., 0., 0., 0.); // reset
	// subtract jets
	float tmp_minJetPt = gMinJetPt;
	gMinJetPt = 15.;
	for (int i=0; i<NJets; i++) {
		if (!isGoodJet(i)) continue;
		tmp.SetPtEtaPhiE(JetPt[i], JetEta[i], JetPhi[i], JetEnergy[i]);
		umet += tmp;
		jets += tmp;
		tmp.SetPtEtaPhiE(0., 0., 0., 0.);
	}
	// subtract muons
	for (int i=0; i<NMus; i++) {
		if (!isGoodMuon(i, 5.)) continue;
		tmp.SetPtEtaPhiM(MuPt[i], MuEta[i], MuPhi[i], gMMU);
		umet += tmp;
		leps += tmp;
		tmp.SetPtEtaPhiM(0., 0., 0., 0.);
	}
	// subtract electrons
	for (int i=0; i<NEls; i++) {
		if (!isGoodElectron(i, 5.)) continue;
		tmp.SetPtEtaPhiM(ElPt[i], ElEta[i], ElPhi[i], gMEL);
		umet += tmp;
		leps += tmp;
		tmp.SetPtEtaPhiM(0., 0., 0., 0.);
	}
	// scale the unclustered energy by 10%
	if (flag == 0) tmp.SetPtEtaPhiE(1.1 * umet.Pt(), umet.Eta(), umet.Phi(), umet.E());
	if (flag == 1) tmp.SetPtEtaPhiE(0.9 * umet.Pt(), umet.Eta(), umet.Phi(), umet.E());
	// subtract the leptons and jets again
	tmp -= leps;
	tmp -= jets;
	// reset the minPt cut for the jets . this is important!
	gMinJetPt = tmp_minJetPt;
	// set the new MET value -- ATTENTION: this is only for pfMET, not Type1 corrected MET
	pfMET = tmp.Pt();
	
}
void SSDLDumper::propagateMET(TLorentzVector nVec, TLorentzVector oVec){
	TLorentzVector met;
	met.SetPtEtaPhiM(getMET(), 0., getMETPhi(), 0.);
	// set the pfMET to the old MET minus original vector plus new vector
	pfMET = (met-oVec+nVec).Pt();
}
void SSDLDumper::smearJetPts(Sample *S, int flag){
	// Modify the jet pt for systematics studies
	// Either shifted or smeared
	// propagate to the MET!!
	if(S->datamc == 0) return; // don't smear data
	if(flag == 0) return;
	TLorentzVector ojets, jets, tmp;
	for(size_t i = 0; i < NJets; ++i){
		tmp.SetPtEtaPhiE(JetPt[i], JetEta[i], JetPhi[i], JetEnergy[i]);
		ojets += tmp;
		if(flag == 1){ JetPt[i] += JetJEC[i]*JetPt[i]; continue; }
		if(flag == 2){ JetPt[i] -= JetJEC[i]*JetPt[i]; continue; }
		if(flag == 3){
			if(JetGenPt[i] > -100.)	JetPt[i] = TMath::Max(float(0.), JetGenPt[i] + getJERScale(i)*(JetPt[i] - JetGenPt[i]));
			else                    JetPt[i] = JetPt[i] * fRand3->Gaus(1., fabs(1. - getJERScale(i)));
		}
		tmp.SetPtEtaPhiE(JetPt[i], JetEta[i], JetPhi[i], JetEnergy[i]);
		jets += tmp;
	}
	propagateMET(jets, ojets);
}
void SSDLDumper::scaleLeptons(Sample *S, int flag){
	// Shift the lepton pts for systematics studies
	if(S->datamc == 0) return; // don't smear data
	if(flag == 0) return;
	float scale = 0.02;
	TLorentzVector oleps, leps, tmp;
	for(size_t i = 0; i < NMus; ++i){
		tmp.SetPtEtaPhiM(MuPt[i], MuEta[i], MuPhi[i], gMMU);
		oleps += tmp;
		// marc scale = getMuScale(MuPt[i], MuEta[i]);
		if(flag == 1) MuPt[i] += scale*MuPt[i];
		if(flag == 2) MuPt[i] -= scale*MuPt[i];
		tmp.SetPtEtaPhiM(MuPt[i], MuEta[i], MuPhi[i], gMMU);
		leps += tmp;
	}
	for(size_t i = 0; i < NEls; ++i){
		tmp.SetPtEtaPhiM(ElPt[i], ElEta[i], ElPhi[i], gMEL);
		oleps += tmp;
		// marc scale = getElScale(ElPt[i], ElEta[i]);
		if(flag == 1) ElPt[i] += scale*ElPt[i];
		if(flag == 2) ElPt[i] -= scale*ElPt[i];
		tmp.SetPtEtaPhiM(ElPt[i], ElEta[i], ElPhi[i], gMEL);
		leps += tmp;
	}
	propagateMET(leps, oleps);
}
void SSDLDumper::smearMET(Sample *S){
	if(S->datamc == 0) return; // don't smear data
	//TRandom3 * myRand = new TRandom3(0);
	//float sm_met = pfMET + myRand->Gaus(0, 0.05) * pfMET;
	float sm_met = pfMET + fRand3->Gaus(0, 0.05) * pfMET;
	pfMET = sm_met;
}
float SSDLDumper::getJetPt(int i){
	return JetPt[i];
}
float SSDLDumper::getMET(){
  	return pfMETType1;
	//return pfMET;
}
float SSDLDumper::getMETPhi(){
	// Return the METPhi, either the true one or the one corrected for applied
	// JES/JER smearing/scaling
  	// float phi = pfMETType1Phi;
	float phi = pfMETPhi;
  	return phi;
}
float SSDLDumper::getM3(){
	// Return M3
	int njets = getNJets();
	if(njets < 3) return -1.;

	float triJetPtMax = 0.;
	float m3 = 0.;
	for(size_t i = 0; i < NJets; ++i){
		if(!isGoodJet(i)) continue;
		for(size_t j = i+1; j < NJets; ++j){
			if(!isGoodJet(j)) continue;
			for(size_t k = j+1; k < NJets; ++k){
				if(!isGoodJet(k)) continue;
				TLorentzVector jet1; jet1.SetPtEtaPhiE(JetPt[i], JetEta[i], JetPhi[i], JetEnergy[i]);
				TLorentzVector jet2; jet2.SetPtEtaPhiE(JetPt[j], JetEta[j], JetPhi[j], JetEnergy[j]);
				TLorentzVector jet3; jet3.SetPtEtaPhiE(JetPt[k], JetEta[k], JetPhi[k], JetEnergy[k]);
				TLorentzVector triJet = jet1 + jet2 + jet3;
				if( triJet.Pt() > triJetPtMax ){
					triJetPtMax = triJet.Pt();
					m3 = triJet.M();
				}
			}
		}
	}	
	return m3;
}
int SSDLDumper::getNJets(){
	int njets(0);
	for(size_t i = 0; i < NJets; ++i) if(isGoodJet(i)) njets++;
	return njets;
}
int SSDLDumper::getNBTags(){
	int ntags(0);
	for(size_t i = 0; i < NJets; ++i) if(isGoodJet(i) && JetCSVBTag[i] > 0.244) ntags++;
	return ntags;
}
int SSDLDumper::getNBTagsMed(){
	int ntags(0);
	for(size_t i = 0; i < NJets; ++i) if(isGoodJet(i) && JetCSVBTag[i] > 0.679) ntags++;
	return ntags;
}
std::vector< int > SSDLDumper::getNBTagsMedIndices(){
	std::vector< int > tagIndices;
	for(size_t i = 0; i < NJets; ++i) {
		if(isGoodJet(i) && JetCSVBTag[i] > 0.679) {
			tagIndices.push_back(i);
		}
	}
	return tagIndices;
}
float SSDLDumper::getHT(){
	float ht(0.);
	for(size_t i = 0; i < NJets; ++i) if(isGoodJet(i)) ht += getJetPt(i);
	return ht;
}
float SSDLDumper::getWeightedHT(){
	float ht(0.);
	for(size_t i = 0; i < NJets; ++i) if(isGoodJet(i)) ht += getJetPt(i)*exp(-fabs(JetEta[i]));
	return ht;
}
float SSDLDumper::getMT(int ind, gChannel chan){
  // Calculates MT
  
  TLorentzVector pmet, plep;
  if (chan == Muon) plep.SetPtEtaPhiM(MuPt[ind], MuEta[ind], MuPhi[ind], gMMU);
  if (chan == Elec) plep.SetPtEtaPhiM(ElPt[ind], ElEta[ind], ElPhi[ind], gMEL);

  pmet.SetPtEtaPhiM(getMET(), 0., getMETPhi(), 0.);
  double ETlept = sqrt(plep.M2() + plep.Perp2());

  double MT = sqrt( 2*(ETlept*getMET() - plep.Px()*pmet.Px() - plep.Py()*pmet.Py()));
  return MT;
}

float SSDLDumper::getMT2(int ind1, int ind2, gChannel chan){
	// Calculate MT2 variable for two leptons and missing energy,
	// assuming zero testmass
	double pa[3];
	double pb[3];
	double pmiss[3];

	TLorentzVector pmet, pl1, pl2;

	if(chan == Muon){ // mumu
		pl1.SetPtEtaPhiM(MuPt[ind1], MuEta[ind1], MuPhi[ind1], gMMU);
		pl2.SetPtEtaPhiM(MuPt[ind2], MuEta[ind2], MuPhi[ind2], gMMU);			
	}
	if(chan == Elec){ // ee
		pl1.SetPtEtaPhiM(ElPt[ind1], ElEta[ind1], ElPhi[ind1], gMEL);
		pl2.SetPtEtaPhiM(ElPt[ind2], ElEta[ind2], ElPhi[ind2], gMEL);			
	}
	if(chan == ElMu){ // emu
		pl1.SetPtEtaPhiM(MuPt[ind1], MuEta[ind1], MuPhi[ind1], gMMU);
		pl2.SetPtEtaPhiM(ElPt[ind2], ElEta[ind2], ElPhi[ind2], gMEL);			
	}
	
	pmet.SetPtEtaPhiM(getMET(), 0., getMETPhi(), 0.);
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
float SSDLDumper::getMll(int ind1, int ind2, gChannel chan){
	// Calculate inv mass for two leptons
	TLorentzVector pl1, pl2;

	if(chan == Muon){ // mumu
		pl1.SetPtEtaPhiM(MuPt[ind1], MuEta[ind1], MuPhi[ind1], gMMU);
		pl2.SetPtEtaPhiM(MuPt[ind2], MuEta[ind2], MuPhi[ind2], gMMU);			
	}
	if(chan == Elec){ // ee
		pl1.SetPtEtaPhiM(ElPt[ind1], ElEta[ind1], ElPhi[ind1], gMEL);
		pl2.SetPtEtaPhiM(ElPt[ind2], ElEta[ind2], ElPhi[ind2], gMEL);			
	}
	if(chan == ElMu){ // emu
		pl1.SetPtEtaPhiM(MuPt[ind1], MuEta[ind1], MuPhi[ind1], gMMU);
		pl2.SetPtEtaPhiM(ElPt[ind2], ElEta[ind2], ElPhi[ind2], gMEL);			
	}
	
	float mass = (pl1+pl2).M();
	return mass;
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
	if(jind > -1) return getJetPt(jind);
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
float SSDLDumper::getClosestJetDPhi(int ind, gChannel chan){
	// Get delta phi to the closest jet
	float lepphi = (chan==Muon)?MuPhi[ind]:ElPhi[ind];
	int jind = getClosestJet(ind, chan);
	if(jind > -1) return fabs(lepphi - JetPhi[jind]);
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
		return getJetPt(i); // assume sorted by pt, so this will return the hardest one
	}
	return 0.;
}
float SSDLDumper::getMaxJPt(){
	float maxpt(0.);
	for(size_t i = 0; i < NJets; ++i){
		if(!isGoodJet(i)) continue;
		if(getJetPt(i) < maxpt) continue;
		maxpt = getJetPt(i);
	}
	return maxpt;
}

int SSDLDumper::getNTightMuons(){
	int nmus = 0;
	for(size_t i = 0; i < NMus; ++i) if(isTightMuon(i)) nmus++;
	return nmus;
}
int SSDLDumper::getNTightElectrons(){
	int nels = 0;
	for(size_t i = 0; i < NEls; ++i) if(isTightElectron(i)) nels++;
	return nels;
}

void SSDLDumper::resetHypLeptons(){
	TLorentzVector vec(0., 0., 0., 0.);
	fHypLepton1 = lepton(vec, 0, -1, -1);
	fHypLepton2 = lepton(vec, 0, -1, -1);
	fHypLepton3 = lepton(vec, 0, -1, -1);
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
void SSDLDumper::setHypLepton3(int index, gChannel chan){
	TLorentzVector vec(0., 0., 0., 0.);
	if(chan == Muon){
		vec.SetPtEtaPhiM(MuPt[index], MuEta[index], MuPhi[index], gMMU);
		fHypLepton3 = lepton(vec, MuCharge[index], 0, index);
	}
	else if(chan == Elec){
		vec.SetPtEtaPhiM(ElPt[index], ElEta[index], ElPhi[index], gMEL);
		fHypLepton3 = lepton(vec, ElCharge[index], 1, index);
	}
	else exit(-1);
}

//////////////////////////////////////////////////////////////////////////////
// Event Selections:
void SSDLDumper::setRegionCuts(int reg){
	fC_minMu1pt   = gRegions[reg]->minMu1pt   ;
	fC_minMu2pt   = gRegions[reg]->minMu2pt   ;
	fC_minEl1pt   = gRegions[reg]->minEl1pt   ;
	fC_minEl2pt   = gRegions[reg]->minEl2pt   ;
	fC_minHT      = gRegions[reg]->minHT      ;
	fC_maxHT      = gRegions[reg]->maxHT      ;
	fC_minMet     = gRegions[reg]->minMet     ;
	fC_maxMet     = gRegions[reg]->maxMet     ;
	fC_minNjets   = gRegions[reg]->minNjets   ;
	fC_maxNjets   = gRegions[reg]->maxNjets   ;
	fC_minNbjets  = gRegions[reg]->minNbjets  ;
	fC_maxNbjets  = gRegions[reg]->maxNbjets  ;
	fC_minNbjmed  = gRegions[reg]->minNbjmed  ;
	fC_maxNbjmed  = gRegions[reg]->maxNbjmed  ;
	fC_app3rdVet  = gRegions[reg]->app3rdVet  ;
	fC_vetoTTZSel = gRegions[reg]->vetoTTZSel ;
	fC_chargeVeto = gRegions[reg]->chargeVeto ;

	if(reg == gRegion["WZEnriched"])  gInvertZVeto = true;
	else                              gInvertZVeto = false;
}

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

	// Sort these vectors by their flavor and their pt
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
	return (getMET() >= min && getMET() < max);
}
bool SSDLDumper::passesZVetoNew(int l1, int l2, int toggle, float dm){
	if (toggle == 0 && NMus > 2){
		TLorentzVector pmu1, pmu2, pmu3;
		pmu1.SetPtEtaPhiM(MuPt[l1], MuEta[l1], MuPhi[l1], gMMU);
		pmu2.SetPtEtaPhiM(MuPt[l2], MuEta[l2], MuPhi[l2], gMMU);
		// third muon
		for(size_t i = 0; i < NMus; ++i){
			if (i == l1 || i == l2) continue;
			if(isGoodMuonForZVeto(i) == false) continue;
			if(MuCharge[i] == MuCharge[l1]) continue;
			pmu3.SetPtEtaPhiM(MuPt[i], MuEta[i], MuPhi[i], gMMU);
			if( ( fabs((pmu1+pmu3).M()) - gMZ ) < dm) return false;
			if( ( fabs((pmu2+pmu3).M()) - gMZ ) < dm) return false;
		  }
	}

	if (toggle == 1 && (NMus > 1 || NEls > 1) ){
		TLorentzVector pmu1, pel1, pl3;
		pmu1.SetPtEtaPhiM(MuPt[l1], MuEta[l1], MuPhi[l1], gMMU);
		pel1.SetPtEtaPhiM(ElPt[l2], ElEta[l2], ElPhi[l2], gMEL);
		// third lepton
		for(size_t i = 0; i < NMus; ++i){
			if (i == l1) continue;
			if(isGoodMuonForZVeto(i) == false) continue;
			if(MuCharge[i] == MuCharge[l1]) continue;
			pl3.SetPtEtaPhiM(MuPt[i], MuEta[i], MuPhi[i], gMMU);
			if( ( fabs((pmu1+pl3).M()) - gMZ ) < dm) return false;
		  }
		for(size_t i = 0; i < NEls; ++i){
			if (i == l2) continue;
			if(isGoodEleForZVeto(i) == false) continue;
			if(ElCharge[i] == ElCharge[l2]) continue;
			pl3.SetPtEtaPhiM(ElPt[i], ElEta[i], ElPhi[i], gMEL);
			if( ( fabs((pel1+pl3).M()) - gMZ )  < dm) return false;
		  }
	}

	if (toggle == 2 && NEls > 2){
		TLorentzVector pel1, pel2, pel3;
		pel1.SetPtEtaPhiM(ElPt[l1], ElEta[l1], ElPhi[l1], gMEL);
		pel2.SetPtEtaPhiM(ElPt[l2], ElEta[l2], ElPhi[l2], gMEL);
		// third electron
		for(size_t i = 0; i < NEls; ++i){
			if (i == l1 || i == l2) continue;
			if(isGoodEleForZVeto(i) == false) continue;
			if(ElCharge[i] == ElCharge[l1]) continue;
			pel3.SetPtEtaPhiM(ElPt[i], ElEta[i], ElPhi[i], gMEL);
			if( ( fabs((pel1+pel3).M()) - gMZ ) < dm) return false;
			if( ( fabs((pel2+pel3).M()) - gMZ ) < dm) return false;
		  }
	}
	if (gInvertZVeto) return false;
	return true;
}
bool SSDLDumper::passesGammaStarVeto(int l1, int l2, int toggle, float mass){
	if (toggle == 0 && NMus > 2){
		TLorentzVector pmu1, pmu2, pmu3;
		pmu1.SetPtEtaPhiM(MuPt[l1], MuEta[l1], MuPhi[l1], gMMU);
		pmu2.SetPtEtaPhiM(MuPt[l2], MuEta[l2], MuPhi[l2], gMMU);
		// third muon
		for(size_t i = 0; i < NMus; ++i){
			if (i == l1 || i == l2) continue;
			if(isGoodMuonForGammaStarVeto(i) == false) continue;
			if(MuCharge[i] == MuCharge[l1]) continue;
			pmu3.SetPtEtaPhiM(MuPt[i], MuEta[i], MuPhi[i], gMMU);
			if( (pmu1+pmu3).M() < mass) return false;
			if( (pmu2+pmu3).M() < mass) return false;
		  }
	}

	if (toggle == 1 && (NMus > 1 || NEls > 1) ){
		TLorentzVector pmu1, pel1, pl3;
		pmu1.SetPtEtaPhiM(MuPt[l1], MuEta[l1], MuPhi[l1], gMMU);
		pel1.SetPtEtaPhiM(ElPt[l2], ElEta[l2], ElPhi[l2], gMEL);
		// third lepton
		for(size_t i = 0; i < NMus; ++i){
			if (i == l1) continue;
			if(isGoodMuonForGammaStarVeto(i) == false) continue;
			if(MuCharge[i] == MuCharge[l1]) continue;
			pl3.SetPtEtaPhiM(MuPt[i], MuEta[i], MuPhi[i], gMMU);
			if((pmu1+pl3).M()  < mass) return false;
		  }
		for(size_t i = 0; i < NEls; ++i){
			if (i == l2) continue;
			if(isGoodEleForGammaStarVeto(i) == false) continue;
			if(ElCharge[i] == ElCharge[l2]) continue;
			pl3.SetPtEtaPhiM(ElPt[i], ElEta[i], ElPhi[i], gMEL);
			if((pel1+pl3).M()  < mass) return false;
		  }
	}

	if (toggle == 2 && NEls > 2){
		TLorentzVector pel1, pel2, pel3;
		pel1.SetPtEtaPhiM(ElPt[l1], ElEta[l1], ElPhi[l1], gMEL);
		pel2.SetPtEtaPhiM(ElPt[l2], ElEta[l2], ElPhi[l2], gMEL);
		// third electron
		for(size_t i = 0; i < NEls; ++i){
			if (i == l1 || i == l2) continue;
			if(isGoodEleForGammaStarVeto(i) == false) continue;
			if(ElCharge[i] == ElCharge[l1]) continue;
			pel3.SetPtEtaPhiM(ElPt[i], ElEta[i], ElPhi[i], gMEL);
			if( (pel1+pel3).M() < mass) return false;
			if( (pel2+pel3).M() < mass) return false;
		  }
	}
	return true;
}
bool SSDLDumper::passesZVeto(bool(SSDLDumper::*muonSelector)(int), bool(SSDLDumper::*eleSelector)(int), float dm){
// Checks if any combination of opposite sign, same flavor leptons (e or mu)
// has invariant mass closer than dm to the Z mass, returns true if none found
// Default for dm is 15 GeV
  if(NMus > 1){
    for(size_t i = 0; i < NMus-1; ++i){
      if((*this.*muonSelector)(i) == false) continue;
      TLorentzVector pmu1, pmu2;
      pmu1.SetPtEtaPhiM(MuPt[i], MuEta[i], MuPhi[i], gMMU);
      
      // Second muon
      for(size_t j = i+1; j < NMus; ++j){ 
	if((*this.*muonSelector)(j) == false) continue;
	pmu2.SetPtEtaPhiM(MuPt[j], MuEta[j], MuPhi[j], gMMU);
	if(MuCharge[i] == MuCharge[j]) continue;
	if(fabs((pmu1+pmu2).M() - gMZ) < dm) return false;
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
	if(!gApplyZVeto) return true;
	
	bool passZVeto = passesZVeto(&SSDLDumper::isGoodMuonForZVeto, &SSDLDumper::isGoodEleForZVeto, dm);
	// return passesZVeto(&SSDLDumper::isTightMuon, &SSDLDumper::isTightElectron, dm);
	// return passesZVeto(&SSDLDumper::isLooseMuon, &SSDLDumper::isLooseElectron, dm);
	// return !passesZVeto(&SSDLDumper::isGoodMuonForZVeto, &SSDLDumper::isGoodEleForZVeto, dm); //inverted Zveto
	if (gInvertZVeto)  return !passZVeto;
	else               return  passZVeto;
}
bool SSDLDumper::passesChVeto(int ch){
  if (ch == 0) return true;  //DON'T apply the veto
  
  int ind1(-1),ind2(-1);
  if (ch > 0 && isSSLLEvent(ind1, ind2) > 0) return true;
  if (ch < 0 && isSSLLEvent(ind1, ind2) < 0) return true;
  
  return false;
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
bool SSDLDumper::passes3rdLepVeto(){
	// Checks if there is a 3rd lepton (other than hyp leptons) in the event
	// Obviously requires for hyp leptons to be set...
	for(size_t i = 0; i < NMus; ++i){
		if(fHypLepton1.type == 0 && fHypLepton1.index == i) continue;
		if(fHypLepton2.type == 0 && fHypLepton2.index == i) continue;
		if(isGoodMuonFor3rdLepVeto(i)) return false;
	}
	for(size_t i = 0; i < NEls; ++i){
		if(fHypLepton1.type == 1 && fHypLepton1.index == i) continue;
		if(fHypLepton2.type == 1 && fHypLepton2.index == i) continue;
		if(isGoodEleFor3rdLepVeto(i)) return false;
	}
	if (gApplyTauVeto && !passesTauVeto()) return false;

	return true;
}
bool SSDLDumper::passesTauVeto(){
	// Checks if there is a good tau in the event
	for(size_t i = 0; i < NTaus; ++i){
		if(isGoodTau(i)) return false;
	}
	return true;
}
bool SSDLDumper::passesTTZSel(){
	int ind1(-1), ind2(-1);
	gChannel chan = Muon;
	bool pass = true;

	// OS Dilepton:
	// Store leptons in vectors
	vector<lepton> tmp_Leptons_p;
	vector<lepton> tmp_Leptons_m;
	TLorentzVector plep;
	for(size_t i = 0; i < NMus; ++i){
		if(isGoodMuonForTTZ(i, 20.) == false) continue;
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
		if(isGoodEleForTTZ(i, 20.) == false) continue;
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

	// Select OSSF pair with best Z mass
	float bestMZ = 999999999.;
	for(size_t iPlus = 0; iPlus < tmp_Leptons_p.size(); ++iPlus){
		for(size_t iMinus = 0; iMinus < tmp_Leptons_m.size(); ++iMinus){
			// OS pair
			lepton l1(tmp_Leptons_p[iPlus]);
			lepton l2(tmp_Leptons_m[iMinus]);

			// Same flavor
			if(l2.type != l1.type) continue;

			// Minimize difference to Z mass
			TLorentzVector dilepton = l1.p+l2.p;
			if(fabs(dilepton.M() - gMZ) > fabs(bestMZ - gMZ)) continue;
			ind1 = l1.index;
			ind2 = l2.index;
			bestMZ = dilepton.M();
			if(l1.type == 1) chan = Elec;
		}
	}
	if(bestMZ > 1000) pass = false;
	tmp_Leptons_p.clear();
	tmp_Leptons_m.clear();
	
	// Want to use hyp leptons only temporary! Have to set them to original later
	lepton orig_hyplep1 = lepton(fHypLepton1.p, fHypLepton1.charge, fHypLepton1.type, fHypLepton1.index);
	lepton orig_hyplep2 = lepton(fHypLepton2.p, fHypLepton2.charge, fHypLepton2.type, fHypLepton2.index);
	
	setHypLepton1(ind1, chan);
	setHypLepton2(ind2, chan);
	
	TLorentzVector pz = fHypLepton1.p + fHypLepton2.p;
	
	// Cuts on inv. mass and pt:
	if(pz.M() < 81. || pz.M() > 101.) return false;
	if(pz.Pt() < 35.) return false;
	
	// Require 3rd lepton:
	bool found3rd = false;
	for(size_t i = 0; i < NMus; ++i){
		if(fHypLepton1.type == 0 && fHypLepton1.index == i) continue;
		if(fHypLepton2.type == 0 && fHypLepton2.index == i) continue;
		if(!isGoodMuonForTTZ(i, 10.)) continue;
		found3rd = true;
		setHypLepton3(i, Muon);
		break;
	}
	if(found3rd == false){ // only look for 3rd ele if no 3rd mu found
		for(size_t i = 0; i < NEls; ++i){
			if(fHypLepton1.type == 1 && fHypLepton1.index == i) continue;
			if(fHypLepton2.type == 1 && fHypLepton2.index == i) continue;
			if(!isGoodEleForTTZ(i, 10.)) continue;
			found3rd = true;
			setHypLepton3(i, Elec);
			break;
		}
	}
	if(found3rd == false) pass = false;
	
	// Cut on jets and btags
	if(getNJets() < 3)     pass = false;
	if(getHT() < 120.)     pass = false;
	if(getNBTags() < 2)    pass = false;
	if(getNBTagsMed() < 1) pass = false;

	resetHypLeptons();
	lepton fHypLepton1 = lepton(orig_hyplep1.p, orig_hyplep1.charge, orig_hyplep1.type, orig_hyplep1.index);
	lepton fHypLepton2 = lepton(orig_hyplep2.p, orig_hyplep2.charge, orig_hyplep2.type, orig_hyplep2.index);
	return pass;
}
// old gstar bool SSDLDumper::passesGammaStarVeto(){
// old gstar   // Return false if there is an additional loose muon within DR=0.2 of an OS hypothesis muon if M(mm) < 10 GeV.
// old gstar   //
// old gstar   if(NMus > 1){
// old gstar     for(size_t i = 0; i < NMus; ++i){
// old gstar       if(fHypLepton1.type == 0 && fHypLepton1.index == i) continue;
// old gstar       if(fHypLepton2.type == 0 && fHypLepton2.index == i) continue;
// old gstar       if(isGoodMuonForGammaStarVeto(i) == false)          continue;
// old gstar       
// old gstar       if(MuCharge[i] != fHypLepton1.charge && fHypLepton1.type == 0) { 
// old gstar 	TLorentzVector pmu;
// old gstar 	pmu.SetPtEtaPhiM(MuPt[i], MuEta[i], MuPhi[i], gMMU);
// old gstar 	float mll = (pmu+fHypLepton1.p).M();
// old gstar 	float dr  = Util::GetDeltaR(fHypLepton1.p.Eta(), pmu.Eta(), fHypLepton1.p.Phi(), pmu.Phi());
// old gstar 	if(mll < 10. && dr < 0.2) return false;
// old gstar       }
// old gstar       if(MuCharge[i] != fHypLepton2.charge && fHypLepton2.type == 0) { 
// old gstar 	TLorentzVector pmu;
// old gstar 	pmu.SetPtEtaPhiM(MuPt[i], MuEta[i], MuPhi[i], gMMU);
// old gstar 	float mll = (pmu+fHypLepton2.p).M();
// old gstar 	float dr  = Util::GetDeltaR(fHypLepton2.p.Eta(), pmu.Eta(), fHypLepton2.p.Phi(), pmu.Phi());
// old gstar 	if(mll < 10. && dr < 0.2) return false;
// old gstar       }
// old gstar     }
// old gstar   }
// old gstar        
// old gstar   if(NEls > 1){
// old gstar     for(size_t i = 0; i < NEls; ++i){
// old gstar       if(fHypLepton1.type == 1 && fHypLepton1.index == i) continue;
// old gstar       if(fHypLepton2.type == 1 && fHypLepton2.index == i) continue;
// old gstar       if(isGoodEleForGammaStarVeto(i) == false)           continue;
// old gstar       
// old gstar       if(ElCharge[i] != fHypLepton1.charge && fHypLepton1.type == 1) { 
// old gstar 	TLorentzVector pel;
// old gstar 	pel.SetPtEtaPhiM(ElPt[i], ElEta[i], ElPhi[i], gMEL);
// old gstar 	float mll = (pel+fHypLepton1.p).M();
// old gstar 	float dr  = Util::GetDeltaR(fHypLepton1.p.Eta(), pel.Eta(), fHypLepton1.p.Phi(), pel.Phi());
// old gstar 	if(mll < 10. && dr < 0.2) return false;
// old gstar       }
// old gstar       if(ElCharge[i] != fHypLepton2.charge && fHypLepton2.type == 1) { 
// old gstar 	TLorentzVector pel;
// old gstar 	pel.SetPtEtaPhiM(ElPt[i], ElEta[i], ElPhi[i], gMEL);
// old gstar 	float mll = (pel+fHypLepton2.p).M();
// old gstar 	float dr  = Util::GetDeltaR(fHypLepton2.p.Eta(), pel.Eta(), fHypLepton2.p.Phi(), pel.Phi());
// old gstar 	if(mll < 10. && dr < 0.2) return false;
// old gstar       }
// old gstar     }
// old gstar   }
// old gstar   
// old gstar   return true;
// old gstar }
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
int SSDLDumper::isSigSupMuEvent(){
	int mu1(-1), mu2(-1);
	if(hasLooseMuons(mu1, mu2) < 1) return -1;
	setHypLepton1(mu1, Muon);
	if(!passesJet50Cut())  return -1;
	if(getNJets() < 1)     return -1;
	if(getMT(mu1,Muon) > fC_maxMt_Control)  return -1;
	if(getMET() > fC_maxMet_Control)   return -1;
	// int nmus(0);
	// for (int i=0; i< NMus; ++i){
	// 	if (isLooseMuon(i)) nmus++;
	// }
	// if (nmus>1) return -1;
	if(NMus > 1)                    return -1;
	return mu1;
}
bool SSDLDumper::isZMuMuEvent(int &mu1, int &mu2){
	if(hasLooseMuons(mu1, mu2) < 2)  return false;
	// if(!isTightMuon(0) && !isTightMuon(1)) return false; // at least one tight

	if(MuCharge[mu1] == MuCharge[mu2]) return false; // os

	// Z mass window cut
	TLorentzVector p1, p2;
	p1.SetPtEtaPhiM(MuPt[mu1], MuEta[mu1], MuPhi[mu1], gMMU);
	p2.SetPtEtaPhiM(MuPt[mu2], MuEta[mu2], MuPhi[mu2], gMMU);
	double m = (p1+p2).M();
	if(fabs(gMZ - m) > 15.) return false;

	setHypLepton1(mu1, Muon);
	setHypLepton2(mu2, Muon);

	if(getMET() > 20.) return false;
	if(getNJets() < 2) return false;
	return true;
}
int SSDLDumper::isSigSupElEvent(){
	int el1(-1), el2(-1);
	if(hasLooseElectrons(el1, el2) < 1) return -1;
	setHypLepton1(el1, Elec);
	if(!passesJet50Cut())          return -1;
	if(getNJets() < 1)             return -1;
	//	if(ElMT[0] > fC_maxMt_Control) return false;
	if(getMT(el1,Elec) > fC_maxMt_Control) return -1;
	if(getMET() > fC_maxMet_Control)  return -1;
	// int nels(0);
	// for (int i = 0; i < NEls; i++) {
	// 	if (isLooseElectron(i)) nels++;
	// }
	// if (nels > 1) return -1;
	if(NEls > 1)                   return -1;
	return el1;
}
bool SSDLDumper::isZElElEvent(int &el1, int &el2){
	if(hasLooseElectrons(el1, el2) < 2)  return false;

	if(ElCharge[el1] == ElCharge[el2]) return false; // os

	// Z mass window cut
	TLorentzVector p1, p2;
	p1.SetPtEtaPhiM(ElPt[el1], ElEta[el1], ElPhi[el1], gMEL);
	p2.SetPtEtaPhiM(ElPt[el2], ElEta[el2], ElPhi[el2], gMEL);
	double m = (p1+p2).M();
	if(fabs(gMZ - m) > 15.) return false;

	setHypLepton1(el1, Elec);
	setHypLepton2(el2, Elec);

	if(getMET() > 20.) return false;
	if(getNJets() < 2) return false;
	return true;
}
bool SSDLDumper::isZElElChMisIdEvent(int &el1, int &el2){
	if (hasLooseElectrons(el1,el2) < 2) return false;
	if (!isTightElectron(el1))          return false;
	if (!isTightElectron(el2))          return false;
	if (getMET() > 30.)                 return false;
	if (ElMT[0] > 25.)                  return false;
	
	TLorentzVector p1, p2;
	p1.SetPtEtaPhiM(ElPt[el1], ElEta[el1], ElPhi[el1], gMEL);
	p2.SetPtEtaPhiM(ElPt[el2], ElEta[el2], ElPhi[el2], gMEL);
	double m = (p1+p2).M();
	if(fabs(gMZ - m) > 15.)             return false;
	
	setHypLepton1(el1, Elec);
	setHypLepton2(el2, Elec);
	
	return true;
}
//____________________________________________________________________________
bool SSDLDumper::AvoidDoubleCountingOfFakes(Sample *S){
	if (S->getType() == 4) { //Check for non-prompt leptons only if the sample is rare MC
		// If the event is not matched to a pair of prompt leptons in GEN the event is rejected.
		if (!isGenMatchedSUSYMuMuEvent() && 
		    !isGenMatchedSUSYEEEvent()   && 
		    !isGenMatchedSUSYEMuEvent()) return true;
	}
  
  return false;
}

bool SSDLDumper::isGenMatchedSUSYMuMuEvent(){
        int ind1(-1), ind2(-1);
	if(!isSSLLMuEvent(ind1, ind2)) return false;
	if(!isTightMuon(ind1) || !isTightMuon(ind2)) return false;
	if(isPromptMuon(ind1) && isPromptMuon(ind2)) return true;
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
		if(!passesChVeto(fC_chargeVeto))    return false;
		if(fDoCounting) fCounter[Muon].fill(fMMCutNames[5]);
	}
	if(fChargeSwitch == 1){
		if(abs(isOSLLEvent(mu1, mu2)) != 1) return false;
	}
	
	
	// Define hypothesis leptons
	setHypLepton1(mu1, Muon);
	setHypLepton2(mu2, Muon);
	
	if(fChargeSwitch == 0){
		// Need to remove the Z veto in OS case, since these OS events are exactly what 
		// contributes to the SS yield. They would NOT fire the Z veto, since they are
		// misidd as SS events
		if(!passesZVeto()) return false; // no Zs in event
		//if(!passesZVetoNew(mu1, mu2, 0)) return false; // no Zs in event
		if(fDoCounting) fCounter[Muon].fill(fMMCutNames[6]);
	}
	
	if(!passesMllEventVeto(mu1, mu2, 1, 8.)) return false; // no low mass OSSF pairs
	if(!passesGammaStarVeto(mu1, mu2, 0)) return false; // reject GStar
	if(fDoCounting) fCounter[Muon].fill(fMMCutNames[7]);


	if(fC_app3rdVet && !passes3rdLepVeto()) return false; // 3rd lepton veto
	if(fDoCounting) fCounter[Muon].fill(fMMCutNames[8]);

	if(fC_vetoTTZSel && passesTTZSel()) return false; // ttZ veto
	// if(fDoCounting) fCounter[Muon].fill(fMMCutNames[8]);

	if(getNJets() < fC_minNjets || getNJets() > fC_maxNjets) return false;    // njets cut
	if(fDoCounting) fCounter[Muon].fill(fMMCutNames[9]);

	if(getNBTags() < fC_minNbjets || getNBTagsMed() < fC_minNbjmed) return false;    // nbjets cut
	if(getNBTags() > fC_maxNbjets || getNBTagsMed() > fC_maxNbjmed) return false;    // nbjets cut
	if(fDoCounting) fCounter[Muon].fill(fMMCutNames[10]); // FIXME

	if(!passesHTCut(fC_minHT, fC_maxHT) )  return false;    // ht cut
	if(fDoCounting) fCounter[Muon].fill(fMMCutNames[11]);

	if(!passesMETCut(fC_minMet, fC_maxMet) ) return false;    // met cut
	if(fDoCounting) fCounter[Muon].fill(fMMCutNames[12]);

	if(!isGoodSecMuon(mu2)) return false; // pt cuts
	if(fDoCounting) fCounter[Muon].fill(fMMCutNames[13]);

	if(!isGoodPrimMuon(mu1)) return false;
	if(fDoCounting) fCounter[Muon].fill(fMMCutNames[14]);

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
		if(!passesChVeto(fC_chargeVeto))    return false;
		if(fDoCounting) fCounter[Elec].fill(fEECutNames[5]);
	}
	if(fChargeSwitch == 1){
		if(abs(isOSLLEvent(el1, el2)) != 2) return false;
	}

	// Define hypothesis leptons
	setHypLepton1(el1, Elec);
	setHypLepton2(el2, Elec);

	
	if(fChargeSwitch == 0){
		// Need to remove the Z veto in OS case, since these OS events are exactly what 
		// contributes to the SS yield. They would NOT fire the Z veto, since they are
		// misidd as SS events
		if(!passesZVeto()) return false; // no Zs in event
		//if(!passesZVetoNew(el1, el2, 2)) return false; // no Zs in event
		if(fDoCounting) fCounter[Elec].fill(fEECutNames[6]);
	}

	if(!passesMllEventVeto(el1, el2, 2, 8.))      return false; // no low mass OSSF pairs
	if(!passesGammaStarVeto(el1, el2, 2)) return false; // reject GStar
	if(fDoCounting) fCounter[Elec].fill(fEECutNames[7]);


	if(fC_app3rdVet && !passes3rdLepVeto()) return false; // 3rd lepton veto
	if(fDoCounting) fCounter[Elec].fill(fEECutNames[8]);

	if(fC_vetoTTZSel && passesTTZSel()) return false; // ttZ veto
	// if(fDoCounting) fCounter[Elec].fill(fEECutNames[8]);

	if(getNJets() < fC_minNjets || getNJets() > fC_maxNjets) return false;    // njets cut
	if(fDoCounting) fCounter[Elec].fill(fEECutNames[9]);

	if(getNBTags() < fC_minNbjets || getNBTagsMed() < fC_minNbjmed) return false;    // nbjets cut
	if(getNBTags() > fC_maxNbjets || getNBTagsMed() > fC_maxNbjmed) return false;    // nbjets cut
	if(fDoCounting) fCounter[Elec].fill(fEECutNames[10]); // FIXME

	if(!passesHTCut(fC_minHT, fC_maxHT) )  return false;    // ht cut
	if(fDoCounting) fCounter[Elec].fill(fEECutNames[11]);

	if(!passesMETCut(fC_minMet, fC_maxMet) ) return false;    // met cut
	if(fDoCounting) fCounter[Elec].fill(fEECutNames[12]);

	if(!isGoodSecElectron(el2)) return false; // pt cuts
	if(fDoCounting) fCounter[Elec].fill(fEECutNames[13]);

	if(!isGoodPrimElectron(el1)) return false;
	if(fDoCounting) fCounter[Elec].fill(fEECutNames[14]);

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
		if(!passesChVeto(fC_chargeVeto))  return false;
		if(fDoCounting) fCounter[ElMu].fill(fEMCutNames[6]);
	}
	if(fChargeSwitch == 1){
		if(abs(isOSLLEvent(mu, el)) != 3) return false;
		if(fDoCounting) fCounter[ElMu].fill(fEMCutNames[6]);
	}

	// Define hypothesis leptons
	setHypLepton1(mu, Muon);
	setHypLepton2(el, Elec);


	if(fChargeSwitch == 0){
		// Need to remove the Z veto in OS case, since these OS events are exactly what 
		// contributes to the SS yield. They would NOT fire the Z veto, since they are
		// misidd as SS events
		if(!passesZVeto()) return false;
		//if(!passesZVetoNew(mu, el, 1)) return false;
		if(fDoCounting) fCounter[ElMu].fill(fEMCutNames[7]);
	}

	if(!passesGammaStarVeto(mu, el, 1)) return false; // reject GStar

	if(fC_app3rdVet && !passes3rdLepVeto()) return false; // 3rd lepton veto
	if(fDoCounting) fCounter[ElMu].fill(fEMCutNames[8]);

	if(fC_vetoTTZSel && passesTTZSel()) return false; // ttZ veto
	// if(fDoCounting) fCounter[ElMu].fill(fEMCutNames[8]);

	if(getNJets() < fC_minNjets || getNJets() > fC_maxNjets) return false;    // njets cut
	if(fDoCounting) fCounter[ElMu].fill(fEMCutNames[9]);

	if(getNBTags() < fC_minNbjets || getNBTagsMed() < fC_minNbjmed) return false;    // nbjets cut
	if(getNBTags() > fC_maxNbjets || getNBTagsMed() > fC_maxNbjmed) return false;    // nbjets cut	if(getNBTags() > fC_maxNbjets) return false;
	if(fDoCounting) fCounter[ElMu].fill(fEMCutNames[10]); // FIXME

	if(!passesHTCut(fC_minHT, fC_maxHT) )  return false;    // ht cut
	if(fDoCounting) fCounter[ElMu].fill(fEMCutNames[11]);

	if(!passesMETCut(fC_minMet, fC_maxMet) ) return false;    // met cut
	if(fDoCounting) fCounter[ElMu].fill(fEMCutNames[12]);

	if(MuPt[mu] > ElPt[el]){
		if(!isGoodPrimMuon(mu))    return false;
		if(fDoCounting) fCounter[ElMu].fill(fEMCutNames[13]);
		if(!isGoodSecElectron(el)) return false;
		if(fDoCounting) fCounter[ElMu].fill(fEMCutNames[14]);
	}
	else if(MuPt[mu] < ElPt[el]){
		if(!isGoodPrimElectron(el)) return false;
		if(fDoCounting) fCounter[ElMu].fill(fEMCutNames[14]);
		if(!isGoodSecMuon(mu))      return false;
		if(fDoCounting) fCounter[ElMu].fill(fEMCutNames[13]);
	}
	return true;
}

//////////////////////////////////////////////////////////////////////////////
// Object selections:
//////////////////////////////////////////////////////////////////////////////
bool SSDLDumper::passesPtCuts(float pT1, float pT2, int reg, gChannel chan){
	if(chan == ElMu){
		if(pT1 > pT2){
			if(pT1 < gRegions[reg]->minMu1pt) return false;
			if(pT2 < gRegions[reg]->minEl2pt) return false;
		}
		if(pT1 < pT2){
			if(pT1 < gRegions[reg]->minMu2pt) return false;
			if(pT2 < gRegions[reg]->minEl1pt) return false;
		}
	}
	else{
		if(chan == Muon && pT1 < gRegions[reg]->minMu1pt) return false;
		if(chan == Muon && pT2 < gRegions[reg]->minMu2pt) return false;
		if(chan == Elec && pT1 < gRegions[reg]->minEl1pt) return false;
		if(chan == Elec && pT2 < gRegions[reg]->minEl2pt) return false;
	}
	return true;
}

// Muons
//____________________________________________________________________________
bool SSDLDumper::isGoodMuon(int muon, float ptcut){
	if(muon >= NMus) return false; // Sanity check
	if(ptcut < 0.) ptcut = fC_minMu2pt;
	if(MuPt[muon] < ptcut) return false;

	return true;
}
bool SSDLDumper::isGoodMuonForZVeto(int muon){
	// Remove stupid PtE/Pt cut
	if(muon >= NMus) return false; // Sanity check
	
	if(MuPt[muon] < 10.) return false;
       	
	if (MuPFIso[muon] > 0.2) return false;
	return true;
}
bool SSDLDumper::isGoodMuonForGammaStarVeto(int muon){
        // Method for selecting very loose muons to reject gamma* bkg.
        if(muon >= NMus) return false; // Sanity check
	
	if(MuPt[muon] <  5.) return false;
       	
	if (MuPFIso[muon] > 0.2) return false;
	return true;
}
bool SSDLDumper::isGoodMuonFor3rdLepVeto(int muon){
	if(isGoodMuon(muon, 10.) == false)  return false;
	if(MuPassesTightID[muon] != 1) return false;
	if(MuPFIso[muon] > 0.15) return false;

	return true;
}
bool SSDLDumper::isGoodMuonForTTZ(int muon, float pt){
	if(isGoodMuon(muon, pt) == false)  return false;
	if(MuPFIso[muon] > 0.15) return false;

	return true;
}
bool SSDLDumper::isLooseMuon(int muon){
	if(isGoodMuon(muon) == false)  return false; // pT cut

	// Veto dep cuts previously in SSDLAnalysis presel
	if(MuEMVetoEt[muon]  > 4.0)      return false;
	if(MuHadVetoEt[muon] > 6.0)      return false;
	// no longer in 2013 if(MuPtE[muon]/MuPt[muon] > 0.1) return false;
	
	// require to pass tight ID:
	if(MuPassesTightID[muon] != 1) return false;

	// passes ISOLATION
	if (gTTWZ) {
		if(MuPFIso[muon] > 1.00) return false; // using detector isolation for ttWZ as requested by f.p.
	}
	else {
		if(MuPFIso[muon] > 1.00) return false;
	}
	// diable for normal running if (MuD0[muon] > 0.005) return false; // this is for testing only!!!
	return true;
}
bool SSDLDumper::isTightMuon(int muon){
	if(isGoodMuon(muon) == false)  return false; // again
	if(isLooseMuon(muon) == false) return false;
	if (gTTWZ) {
		if(MuPFIso[muon] > gMuMaxIso) return false; // using detector isolation for ttWZ as requested by f.p.
	}
	else {
		if(MuPFIso[muon] > gMuMaxIso) return false;
	}

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
	if(isPromptMuon(muon)) return false;
	return true;
}
bool SSDLDumper::isPromptMuon(int muon){
	if(isLooseMuon(muon) == false) return false;

	// Mother or Grandmother is a SM hadron:
	if(MuGenMType[muon] > 10 || MuGenGMType[muon] > 10) return false;

	// Mother or Grandmother is a SUSY particle:
	if(MuGenMType[muon] == 9 || MuGenGMType[muon] == 9) return true;
	// sparticle -> SM hadron -> lepton is taken care of by previous

	// Mother is Gaugeboson but not gluon or gamma:
	if(MuGenMType[muon] == 4 && abs(MuGenMID[muon]) != 21 && abs(MuGenMID[muon]) != 22) return true;

	// Mother is a tau from W or Z
	if(abs(MuGenMID[muon]) == 15 && MuGenGMType[muon] == 4) return true;
	
	// Mother is a top
	if(abs(MuGenMID[muon]) == 6) return true;
	
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

	// 	if(ElEcalRecHitSumEt[ele]/ElPt[ele] > 0.2) return false; // when using "CaloIsoVL" triggers

	// Reject electrons closer than 0.1 in DR to tight muons
	for(size_t i = 0; i < NMus; ++i){
		if(!isLooseMuon(i)) continue;
		if(!isTightMuon(i)) continue;
		if(Util::GetDeltaR(MuEta[i], ElEta[ele], MuPhi[i], ElPhi[ele]) > 0.1 ) continue;
		return false;
	}
	
	return true;
}
bool SSDLDumper::isGoodEleForZVeto(int ele){
	// Don't care about charge consistency or trigger efficiency
	if(ele >= NEls) return false; // Sanity check
	if(ElPt[ele] < 10.) return false;
	
	// Apply IsoCUT
	if(ElPFIso[ele] > 0.2) return false;
	return true;	
}
bool SSDLDumper::isGoodEleForGammaStarVeto(int ele){
        // Method for selecting very loose muons to reject gamma* bkg.
        if(ele >= NEls) return false; // Sanity check
	if(ElPt[ele] <  5.) return false;
       	
	if(ElPFIso[ele] > 0.2) return false;
	return true;
}
bool SSDLDumper::isGoodEleFor3rdLepVeto(int ele){
	// Don't care about charge consistency or trigger efficiency
	if(ele >= NEls) return false; // Sanity check
	if(ElPt[ele] < 10.) return false;

	if(ElPFIso[ele] > 0.15) return false;

	// Reject electrons closer than 0.1 in DR to tight muons
	for(size_t i = 0; i < NMus; ++i){
		if(!isTightMuon(i)) continue;
		if(Util::GetDeltaR(MuEta[i], ElEta[ele], MuPhi[i], ElPhi[ele]) > 0.1 ) continue;
		return false;
	}

	if(ElIsGoodElId_LooseWP[ele] != 1) return false;
	  
// 	// Additional cuts from AN-12-059_v5:
// 	if( fabs(ElEta[ele]) < 1.479 ){ // Barrel
// 		if(ElHoverE[ele] > 0.10)     return false;
// 		if(fabs(ElDPhi[ele]) > 0.15) return false;
// 	}
// 	if( fabs(ElEta[ele]) >= 1.479 ){ // Endcap
// 		if(ElHoverE[ele] > 0.075)    return false;
// 		if(fabs(ElDPhi[ele]) > 0.10) return false;
// 	}
	
	return true;	
}
bool SSDLDumper::isGoodEleForChMId(int ele, float ptcut){
	// Don't care about charge consistency or trigger efficiency
	if(ele >= NEls) return false; // Sanity check
	if(ElPt[ele] < ptcut) return false;
	
	if (!isLooseElectron(ele)) return false;
	float RelIso = (ElEcalRecHitSumEt[ele] + ElHcalTowerSumEt[ele] + ElTkSumPt[ele])/ElPt[ele];
  
	if (RelIso > 0.15) return false;
	
	return true;
}
bool SSDLDumper::isGoodEleForTTZ(int ele, float pt){
	// Don't care about charge consistency or trigger efficiency
	if(ele >= NEls) return false; // Sanity check
	if(ElPt[ele] < pt) return false;

	// Additional cuts from AN-12-059_v5:
	if( fabs(ElEta[ele]) < 1.479 ){ // Barrel
		if(ElPFIso[ele] > 0.15)     return false;

		if(ElSigmaIetaIeta[ele] > 0.01)  return false;
		if(fabs(ElDPhi[ele])    > 0.15)  return false;
		if(fabs(ElDEta[ele])    > 0.007) return false;
		if(ElHoverE[ele]        > 0.15)  return false;
	}
	if( fabs(ElEta[ele]) >= 1.479 ){ // Endcap
		if(ElPFIso[ele]  > 0.10)     return false;

		if(ElSigmaIetaIeta[ele] > 0.03)  return false;
		if(fabs(ElDPhi[ele])    > 0.10)  return false;
		if(fabs(ElDEta[ele])    > 0.01)  return false;
		if(ElHoverE[ele]        > 0.10)  return false;
	}

	// Reject electrons closer than 0.1 in DR to tight muons
	for(size_t i = 0; i < NMus; ++i){
		if(!isTightMuon(i)) continue;
		if(Util::GetDeltaR(MuEta[i], ElEta[ele], MuPhi[i], ElPhi[ele]) > 0.1 ) continue;
		return false;
	}
	
	return true;	
}
bool SSDLDumper::isLooseElectron(int ele){
	if(isGoodElectron(ele) == false) return false; // pT cut and dR cut with muons (0.1)
	if(ElChIsCons[ele] != 1) return false;
	
	// Additional cuts for CaloIsoVL and TrkIsoVL
	// temporary if(ElEcalRecHitSumEt[ele]/ElPt[ele] > 0.2) return false; // CaloIsoVL
	// temporary if(ElHcalTowerSumEt [ele]/ElPt[ele] > 0.2) return false; // CaloIsoVL
	// temporary if(ElTkSumPt        [ele]/ElPt[ele] > 0.2) return false; // TrkIsoVL
	
	// ISO AND ID FOR TT+W
	if (gTTWZ) {
		if(ElPFIso[ele] > 0.6) return false; // using detector isolation for ttWZ as requested by f.p.
		if(ElIsGoodTriggerEl[ele] != 1) return false; // trigger ID cuts
	}
	// ISO AND ID FOR SUSY
	else {
		if(ElPFIso[ele] > 0.6) return false;
		if(ElIsGoodTriggerEl[ele] != 1) return false; // trigger ID cuts
	}
	// disable for normal running if (ElD0[ele] > 0.01) return false; // this is for testing!!!!

	// JUST FOR FUTURE REFERENCE, THOSE ARE THE VALUES FOR MVA-ID WE USED
	// loose MVA-ID values 	if(fabs(ElEta[ele]) < 0.8   &&                             ElMVAIDTrig[ele] < 0.949) return false;
	// loose MVA-ID values 	if(fabs(ElEta[ele]) > 0.8   && fabs(ElEta[ele]) < 1.479 && ElMVAIDTrig[ele] < 0.921) return false;
	// loose MVA-ID values 	if(fabs(ElEta[ele]) > 1.479 && fabs(ElEta[ele]) < 2.5   && ElMVAIDTrig[ele] < 0.956) return false;

	return true;
}
bool SSDLDumper::isTightElectron(int ele){
	if(!isLooseElectron(ele))       return false;

	// ISO AND ID FOR TT+W
	if (gTTWZ) {
		if(ElPFIso[ele] > gElMaxIso) return false; // using detector isolation for ttWZ as requested by f.p.
		if(ElIsGoodElId_MediumWP[ele] != 1) return false; // using cut based ID again
	}
	// ISO AND ID FOR SUSY
	else {
		if(ElPFIso[ele] > gElMaxIso) return false;
		if(ElIsGoodElId_MediumWP[ele] != 1) return false;
	}

	// JUST FOR FUTURE REFERENCE, THOSE ARE THE VALUES FOR MVA-ID WE USED
	// tight MVA-ID values if(fabs(ElEta[ele]) < 0.8   &&                             ElMVAIDTrig[ele] < 0.956) return false;
	// tight MVA-ID values if(fabs(ElEta[ele]) > 0.8   && fabs(ElEta[ele]) < 1.479 && ElMVAIDTrig[ele] < 0.949) return false;
	// tight MVA-ID values if(fabs(ElEta[ele]) > 1.479 && fabs(ElEta[ele]) < 2.5   && ElMVAIDTrig[ele] < 0.968) return false;

	return true;
}
bool SSDLDumper::isGoodPrimElectron(int ele, float ptcut){
	if(ptcut < 0.) ptcut = fC_minEl1pt;
	if(isLooseElectron(ele) == false) return false;
	if(ElPt[ele] < ptcut)             return false;
	return true;
}
bool SSDLDumper::isGoodSecElectron(int ele, float ptcut){
	if(ptcut < 0.) ptcut = fC_minEl2pt;
	if(isLooseElectron(ele) == false) return false;
	if(ElPt[ele] < ptcut)             return false;
	return true;
}

//____________________________________________________________________________
bool SSDLDumper::isFakeElectron(int ele){
	if(isPromptElectron(ele)) return false;
	return true;
}
bool SSDLDumper::isPromptElectron(int ele){
	if(isLooseElectron(ele) == false) return false;

	// Matched to electron
	if(abs(ElGenID[ele]) != 11) return false;

	// Mother or Grandmother is a SM hadron:
	if(ElGenMType[ele] > 10 || ElGenGMType[ele] > 10) return false;

	// Mother or Grandmother is a SUSY particle:
	if(ElGenMType[ele] == 9 || ElGenGMType[ele] == 9) return true;
	// sparticle -> SM hadron -> lepton is taken care of by previous

	// Mother is Gaugeboson but not gluon or gamma:
	if(ElGenMType[ele] == 4 && abs(ElGenMID[ele]) != 21 && abs(ElGenMID[ele]) != 22) return true;

	// Mother is a tau from W or Z
	if(abs(ElGenMID[ele]) == 15 && ElGenGMType[ele] == 4) return true;
	
	// Mother is a top
	if(abs(ElGenMID[ele]) == 6) return true;
	
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
// Taus
//____________________________________________________________________________
bool SSDLDumper::isGoodTau(int tau){
// Good taus according to EWino people
	if(tau >= NTaus) return false; // Sanity check
	if(TauPt[tau]  < 20) return false;
	if(fabs(TauEta[tau]) > 2.3) return false;
	if( TauMVAElRej[tau] < 0.5 ) return false;
	if( TauTightMuRej[tau]   < 0.5 ) return false;
	if( TauLCombIsoDB[tau] < 0.5 ) return false;

	// Reject taus closer than 0.1 in DR to tight muons and tight electrons
	for(size_t i = 0; i < NMus; ++i){
		if(!isGoodMuonFor3rdLepVeto(i)) continue;
		if(Util::GetDeltaR(MuEta[i], TauEta[tau], MuPhi[i], TauPhi[tau]) > 0.1 ) continue;
		return false;
	}
	for(size_t i = 0; i < NEls; ++i){
		if(!isGoodEleFor3rdLepVeto(i)) continue;
		if(Util::GetDeltaR(ElEta[i], TauEta[tau], ElPhi[i], TauPhi[tau]) > 0.1 ) continue;
		return false;
	}
	
	return true;
}

//////////////////////////////////////////////////////////////////////////////
// Jets
//____________________________________________________________________________
bool SSDLDumper::isGoodJet(int jet, float pt){
	if(jet >= NJets) return false; // Sanity check
	// JET - LEPTON CLEANING PARAMETER!! should switch to 0.5!!!
	float minDR = 0.4;

	if(getJetPt(jet) < gMinJetPt) return false;
	if(getJetPt(jet) < pt) return false;
	
	if(fabs(JetEta[jet]) > gMaxJetEta) return false; // btagging only up to 2.4
	
	// Remove jets close to hypothesis leptons
	if(fHypLepton1.index > -1) if(Util::GetDeltaR(fHypLepton1.p.Eta(), JetEta[jet], fHypLepton1.p.Phi(), JetPhi[jet]) < minDR) return false;
	if(fHypLepton2.index > -1) if(Util::GetDeltaR(fHypLepton2.p.Eta(), JetEta[jet], fHypLepton2.p.Phi(), JetPhi[jet]) < minDR) return false;
	if(fHypLepton3.index > -1) if(Util::GetDeltaR(fHypLepton3.p.Eta(), JetEta[jet], fHypLepton3.p.Phi(), JetPhi[jet]) < minDR) return false;

	// Remove jets close to all tight leptons
	for(size_t imu = 0; imu < NMus; ++imu){
		if(!isTightMuon(imu)) continue;
		if(!isGoodSecMuon(imu)) continue; // pt  > 10
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
float SSDLDumper::getJERScale(int jet){
	float eta = JetEta[jet];
	if(     fabs(eta) < 0.5) return 1.052;
	else if(fabs(eta) < 1.1) return 1.057;
	else if(fabs(eta) < 1.7) return 1.096;
	else if(fabs(eta) < 2.3) return 1.134;
	else                     return 1.288;
}

//////////////////////////////////////////////////////////////////////////////
// HLT Scale Factors
//____________________________________________________________________________
float SSDLDumper::getHLTSF_DoubleElectron( float pt, float eta, const std::string& runPeriod ) {
	return 1.;
}
float SSDLDumper::getHLTSF_MuEG( float pt, float eta, const std::string& runPeriod ) {
	return sqrt(0.95);
}
float SSDLDumper::diMuonHLTSF2012(){
	return 0.872;
}
float SSDLDumper::muEleHLTSF2012(){
	return 0.925;
}
float SSDLDumper::diEleHLTSF2012(){
	return 0.954;
}
SSDLDumper::ValueAndError SSDLDumper::getMuMuHLTSF( float pt, float eta, const std::string& runPeriod ) {
	// These numbers taken from AN2011-399-v4
	float hltsf = 0.;
	float hltsf_err = 0.;

	if( runPeriod=="Run2011A" ) {
		if( fabs(eta)<0.8 ) {
			hltsf = 0.975;
			hltsf_err = 0.004;
		} else if( fabs(eta)<2.1 ) {
			hltsf = 0.950;
			hltsf_err = 0.005;
		} else  {
			hltsf = 0.910;
			hltsf_err = 0.01;
		}
	} else if( runPeriod=="Run2011B" ) {
		if( fabs(eta)<0.8 ) {
			if( pt<40. ) {
				hltsf = 0.977;
				hltsf_err = 0.001;
			} else {
				hltsf = 0.975;
				hltsf_err = 0.001;
			}
		} else if( fabs(eta)<2.1 ){
			if( pt<40. ) {
				hltsf = 0.955;
				hltsf_err = 0.002;
			} else {
				hltsf = 0.955;
				hltsf_err = 0.001;
			}
		} else { // eta 2.1 -> 2.4
			if( pt<40. ) {
				hltsf = 0.89;
				hltsf_err = 0.007;
			} else {
				hltsf = 0.90;
				hltsf_err = 0.006;
			}
		}
	} else if( runPeriod=="Run2011") {
		hltsf     = ( 2.3*this->getMuMuHLTSF(pt, eta, "Run2011A").val + 2.7*this->getMuMuHLTSF(pt, eta, "Run2011B").val ) / 5.;
		hltsf_err = ( 2.3*this->getMuMuHLTSF(pt, eta, "Run2011A").err + 2.7*this->getMuMuHLTSF(pt, eta, "Run2011B").err ) / 5.;
	} else {
		std::cout << "WARNING! Unknown run period: " << runPeriod << "! Returning HLTSF=0." << std::endl;
	}

	ValueAndError ve_hlt;
	ve_hlt.val = hltsf;
	ve_hlt.err = hltsf_err;

	return ve_hlt;
}

float SSDLDumper::getMuScale(float pt, float eta){
	// Returns the combined error of hlt/id/iso/reco scale factors
	ValueAndError HLTSF_Run2011A = getMuMuHLTSF(pt, eta, "Run2011A");
	ValueAndError HLTSF_Run2011B = getMuMuHLTSF(pt, eta, "Run2011B");

	ValueAndError hltscale;
	hltscale.val = ((217.+920.+1000.) * HLTSF_Run2011A.val + 2100. * HLTSF_Run2011B.val)/(217.+920.+1000.+2100.);
	hltscale.err = ((217.+920.+1000.) * HLTSF_Run2011A.err + 2100. * HLTSF_Run2011B.err)/(217.+920.+1000.+2100.);
	ValueAndError recoscale = getMuonRecoSF(pt, eta);
	ValueAndError isoscale  = getMuonIsoSF(pt, eta);
	
	float sqscale = hltscale.err*hltscale.err * recoscale.val*recoscale.val * isoscale.val*isoscale.val
	              + hltscale.val*hltscale.val * recoscale.err*recoscale.err * isoscale.val*isoscale.val
	              + hltscale.val*hltscale.val * recoscale.val*recoscale.val * isoscale.err*isoscale.err;
	return sqrt(sqscale);
}
float SSDLDumper::getElScale(float pt, float eta){
	// Returns the combined error of hlt/id/iso/reco scale factors
	ValueAndError hltscale; hltscale.val = 1.0; hltscale.err = 0.01;
	ValueAndError recoscale = getElectronRecoSF(pt, eta);
	ValueAndError isoscale  = getElectronIsoSF(pt, eta);
	
	float sqscale = hltscale.err*hltscale.err * recoscale.val*recoscale.val * isoscale.val*isoscale.val
	              + hltscale.val*hltscale.val * recoscale.err*recoscale.err * isoscale.val*isoscale.val
	              + hltscale.val*hltscale.val * recoscale.val*recoscale.val * isoscale.err*isoscale.err;
	return sqrt(sqscale);
}

SSDLDumper::ValueAndError SSDLDumper::getMuonRecoSF( float pt, float eta ) {

	float recoSF;
	float recoSF_err;

	if( fabs(eta)<1.2 ) {
		recoSF = 0.996;
		recoSF_err = 0.001;
	} else {
		recoSF = 0.986;
		recoSF_err = 0.001;
	}

	ValueAndError ve_reco;
	ve_reco.val = recoSF;
	ve_reco.err = recoSF_err;

	return ve_reco;

}
SSDLDumper::ValueAndError SSDLDumper::getElectronRecoSF( float pt, float eta ) {

	float recoSF;
	float recoSF_err;

	if( fabs(eta)<0.8 ) {
		recoSF = 0.999;
		recoSF_err = 0.005;
	} else if( fabs(eta)<1.44 ) {
		recoSF = 0.964;
		recoSF_err = 0.003;
	} else if( fabs(eta)<1.57 ) {
		recoSF = 0.99;
		recoSF_err = 0.04;
	} else if( fabs(eta)<2.0 ) {
		recoSF = 0.992;
		recoSF_err = 0.006;
	} else {
		recoSF = 1.001;
		recoSF_err = 0.006;
	}

	ValueAndError ve_reco;
	ve_reco.val = recoSF;
	ve_reco.err = recoSF_err;

	return ve_reco;

}
SSDLDumper::ValueAndError SSDLDumper::getMuonIsoSF( float pt, float eta ) {

	float recoSF;
	float recoSF_err;

	if( pt<40. ) {

		if( fabs(eta)<0.9 ) {
			recoSF = 0.987;
			recoSF_err = 0.006;
		} else {
			recoSF = 0.995;
			recoSF_err = 0.005;
		}

	} else {

		if( fabs(eta)<0.9 ) {
			recoSF = 0.994;
			recoSF_err = 0.002;
		} else {
			recoSF = 0.996;
			recoSF_err = 0.002;
		}

	}


	ValueAndError ve_reco;
	ve_reco.val = recoSF;
	ve_reco.err = recoSF_err;

	return ve_reco;

}
SSDLDumper::ValueAndError SSDLDumper::getElectronIsoSF( float pt, float eta ) {

	float recoSF;
	float recoSF_err;

	if( pt<40. ) {

		if( fabs(eta)<1.5 ) {
			recoSF = 0.988;
			recoSF_err = 0.006;
		} else {
			recoSF = 0.998;
			recoSF_err = 0.011;
		}

	} else {

		if( fabs(eta)<1.5 ) {
			recoSF = 0.988;
			recoSF_err = 0.003;
		} else {
			recoSF = 1.016;
			recoSF_err = 0.064;
		}

	}


	ValueAndError ve_reco;
	ve_reco.val = recoSF;
	ve_reco.err = recoSF_err;

	return ve_reco;

}
