/*****************************************************************************
*     Collection of tools for producing plots for October Exercise           *
*                                                                            *
*                                                  (c) 2009 Benjamin Stieger *
*****************************************************************************/
#ifndef ANACLASS_HH
#define ANACLASS_HH

#include "TLatex.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TLine.h"
#include "TArrow.h"
#include "TBox.h"
#include "TMath.h"

#include "TString.h"
#include "TObject.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TH2.h"
#include "TH3.h"
#include "TF1.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TCut.h"
#include "TStyle.h"
#include "TRandom.h"
#include "TROOT.h"
#include "TVirtualPad.h"
#include "TClonesArray.h"

#include "TGraph2D.h"
#include "TMatrixD.h"
#include "TMatrixT.h"
#include "TVectorD.h"
#include "TPaveStats.h"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <vector>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <map>
#include <time.h> // access to date/time

#include "Utilities.hh"

// class AnaClass: public TObject {
class AnaClass{

public:
	AnaClass();
	AnaClass(const char*, bool = false);
	virtual ~AnaClass();

/*****************************************************************************
##################| Initialization and Setup |################################
*****************************************************************************/
	void init(bool = false);
	void loadSamples(TString parfile = "parfile.dat", bool verbose = false);
	void readParms(TString filename, bool verbose = false);
	void readVarNames(const char* = "varnames.dat");
/*****************************************************************************
##################| Produce Plots |###########################################
*****************************************************************************/
	void plotPlotList(const char* = "plotlist.dat", TTree *tree = NULL, TString tag = "tag");
	void plotPlotList2D(const char* = "plotlist2d.dat", TTree *tree = NULL, TString tag = "tag");
	void makePlotPanel(const TCut req = "", const int nbins = 100);
	void plotAllBranches(TTree*, TString);
	
	void plotAllBranches(int);
/*****************************************************************************
###################| Main Methods |###########################################
*****************************************************************************/
	void plotEID(TCut, TTree*, TString);

/*****************************************************************************
##################| Utilities |###############################################
*****************************************************************************/
	TTree* getTree(TString treename, TString filename, TString subdir = "");
	TH1D* drawTree1D(const char* arg, const TCut reqs, const char* histn, const int nbins, const double xmin, const double xmax, TTree* t, bool draw, const char* drawopt = "");
	TH2D* drawTree2D(const char* arg1, const char* arg2, const TCut reqs, const char* histn, const int nbinsx, const double xmin, const double xmax, const int nbinsy, const double ymin, const double ymax, TTree *tree, bool draw, const char* drawopt = "");
	void plotVar(const char* var, const TCut reqs, TTree *tree, TString tag, int nbins, double xmin, double xmax, TString ofilename="ofilename", bool logy = false, double line1x = -999., double line2x = -999.);
	void refValues(const char* var, TH1D* h);
	double tailFraction(TH1D* h, double frac);
	void printCheckList(const char* var, TH1D* h, const char* filename);
	TString printTailFraction(const char* var, TH1D* h, double frac);
	TString printAverage(const char* var, TH1D* h);
	TString printRatio(const char* var, TH1D* h, double x1, double x2, double y1, double y2);
	void plotVar2D(const char* var1, const char* var2, const TCut reqs, TTree *tree, TString tag, int nbinsx, double xmin, double xmax, int nbinsy, double ymin, double ymax, Option_t *topt ="", int markstyle = 0, bool logx = false, bool logy = false, bool logz = false, double line1x = -999., double line2x = -999., double line1y = -999., double line2y = -999.);
	void plotOverlay2T(const char* var, const TCut reqs, int index1, int index2, int nbins, double xmin, double xmax, bool logy = false, double line1x = -999., double line2x = -999.);
	void plotOverlay1T2V(const char* var1, const char* var2, const TCut reqs, int sampleindex, int nbins, double xmin, double xmax, bool logy = false, double line1x = -999., double line2x = -999.);
	void plotOverlay2C(const char* var, const TCut req1, const TCut req2, int sampleindex, TString tag1, TString tag2, int nbins, double xmin, double xmax, bool logy = false);
	void plotOverlay3T(const char* var, const TCut reqs, int index1, int index2, int index3, int nbins, double xmin, double xmax, bool logy = false, double line1x = -999., double line2x = -999.);
	void plotOverlay3C(const char* var, const TCut req1, TString tag1, const TCut req2, TString tag2, const TCut req3, TString tag3, int sampleindex, int nbins, double xmin, double xmax, bool logy = false);
	void plotOverlay4T(const char* var, const TCut reqs, int index1, int index2, int index3, int index4, int nbins, double xmin, double xmax, bool logy = false, double line1x = -999., double line2x = -999.);
	TString convertVarName(const char* var);
	TString convertVarName2(const char* var);
	int OptNBins(int);
	TH1D* normHist(const TH1D *ihist);
	TH2D* normHist(const TH2D *ihist);
	TCanvas* makeCanvas(const char*);
	TString numbForm(double);
	int getExp(double e);
	
/*****************************************************************************
##################| Variables |###############################################
*****************************************************************************/
	inline void setOutputDir(TString dir){ fOutputDir = Util::MakeOutputDir(dir); };
	inline void setOutputSubDir(TString dir){ fOutputSubDir = Util::MakeOutputDir(dir); };
	inline void setCheckListFile(TString file){ fChecklistFile = file; };
	inline void setGlobalTag(TString tag){ fGlobalTag = tag; };
	TFile *fFile[20];
	TString fFileName[20];
	TString fTreeSubDirName[20]; // Name of subdirectory inside rootfile where tree is
	TTree *fTree[20];
	TString fTag[20];
	double fNorm[20];
	std::map<TString, TString> fVarNameMap;	// Mapping of axis names

	int fFont;				// Global font value for labels, titles, legends
	int fBGColor;			// Choose color for background (temporary)
	
	TString fGlobalTag;	// Global Tag for a certain parameter set
	TLatex *fTlat;			// TLatex object for generic use
	TString fParFile;		// FileName of Parameter File

	TCut fCuts;
	TCut fL1Cuts;			// General cuts
	TCut fL2Cuts;			// Fake suppression cuts
	TCut fL1L2Cuts;		// AND of L1 and L2 cuts

	TString fOutputDir;	  // Output directory for plots
	TString fOutputSubDir; // Output subdirectory for plots, changes for each use
	TString fChecklistFile; // Name of checklist file

private:
	// ClassDef(AnaClass,2)
};

#endif
