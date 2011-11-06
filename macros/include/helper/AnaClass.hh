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
#include "TTreeFormula.h"
#include "TStyle.h"
#include "TRandom.h"
#include "TROOT.h"
#include "TVirtualPad.h"
#include "TClonesArray.h"

#include "TEfficiency.h"
#include "TGraphAsymmErrors.h"

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
#include "MetaTreeClassBase.h"

// class AnaClass: public TObject {
class AnaClass: public MetaTreeClassBase{

public:
	AnaClass();
	AnaClass(const char*, bool = false);
	virtual ~AnaClass();

/*****************************************************************************
##################| Initialization and Setup |################################
*****************************************************************************/
	virtual void init(bool = false); // Careful, MakeClass produces Init with capital I!
	virtual void loadSamples(TString parfile = "parfile.dat", bool verbose = false);
	virtual void readParms(TString filename, bool verbose = false);
	virtual void readVarNames(const char* = "varnames.dat");
/*****************************************************************************
##################| Produce Plots |###########################################
*****************************************************************************/
	virtual void plotPlotList(const char* = "plotlist.dat", TTree *tree = NULL, TString tag = "tag", TCut cut ="", TFile* file = 0 );
	virtual void plotPlotList2D(const char* = "plotlist2d.dat", TTree *tree = NULL, TString tag = "tag", TFile* file = 0 );
	virtual void plotAllBranches(TTree*, TString);
	
	virtual void plotAllBranches(int);
	virtual void plotEID(TCut, TTree*, TString, TFile* file = 0 );

	virtual void plotVar(const char* var, const TCut reqs, TTree *tree, TString tag, int nbins, double xmin, double xmax, TString ofilename="ofilename", bool logy = false, double line1x = -999., double line2x = -999., TFile* file = 0 );
	virtual void plotVar2D(const char* var1, const char* var2, const TCut reqs, TTree *tree, TString tag, int nbinsx, double xmin, double xmax, int nbinsy, double ymin, double ymax, Option_t *topt ="", int markstyle = 0, bool logx = false, bool logy = false, bool logz = false, double line1x = -999., double line2x = -999., double line1y = -999., double line2y = -999., TFile* file = 0 );
	virtual void plotOverlay2H(TH1D *h1, TString tag1, TH1D *h2, TString tag2, bool logy = false, double line1x = -999.99, double line2x = -999.99);
	virtual void plotOverlay2HNorm(TH1D *h1, TString tag1, TH1D *h2, TString tag2, bool logy = false, double line1x = -999.99, double line2x = -999.99);
	virtual void plotOverlay3H(TH1D *h1, TString tag1, TH1D *h2, TString tag2, TH1D *h3, TString tag3, bool logy = false, double line1x = -999.99, double line2x = -999.99);
	virtual void plotOverlay2T(const char* var, const TCut reqs, int index1, int index2, int nbins, double xmin, double xmax, bool logy = false, double line1x = -999., double line2x = -999.);
	virtual void plotOverlay1T2V(const char* var1, const char* var2, const TCut reqs, int sampleindex, int nbins, double xmin, double xmax, bool logy = false, double line1x = -999., double line2x = -999.);
	virtual void plotOverlay2C(const char* var, const TCut req1, const TCut req2, int sampleindex, TString tag1, TString tag2, int nbins, double xmin, double xmax, bool logy = false);
	virtual void plotOverlay3T(const char* var, const TCut reqs, int index1, int index2, int index3, int nbins, double xmin, double xmax, bool logy = false, double line1x = -999., double line2x = -999.);
	virtual void plotOverlay3C(const char* var, const TCut req1, TString tag1, const TCut req2, TString tag2, const TCut req3, TString tag3, int sampleindex, int nbins, double xmin, double xmax, bool logy = false);
	virtual void plotOverlay4T(const char* var, const TCut reqs, int index1, int index2, int index3, int index4, int nbins, double xmin, double xmax, bool logy = false, double line1x = -999., double line2x = -999.);

	virtual void plotPredOverlay2H(TH1D *h1, TString tag1, TH1D *h2, TString tag2, bool logy = false, double line1x = -999.99, double line2x = -999.99);	
	virtual void plotPredOverlay2HWithRatio(TH1D *h1, TString tag1, TH1D *h2, TString tag2, bool logy = false, bool ratio = true, double line1x = -999.99, double line2x = -999.99);	
	virtual void plotPredOverlay2HWithRatio(THStack *h1s, TH1D *h1, TString tag1, TH1D *h2, TString tag2, bool logy = false, bool ratio = true, double line1x = -999.99, double line2x = -999.99);	
	virtual void plotPredOverlay3HWithRatio(TH1D *h1, TString tag1, TH1D *h2, TString tag2, TH1D *h3, TString tag3, bool logy = false, bool ratio = true, double line1x = -999.99, double line2x = -999.99);
	virtual void plotRatioOverlay2H(TH1D *h1, TString tag1, TH1D *h2, TString tag2, bool logy = false, double line1x = -999.99, double line2x = -999.99);
	virtual void plotEffOverlayEG(TEfficiency *h1, TString tag1, TGraphAsymmErrors *h2, TString tag2, bool logy = false);
	virtual void plotEffOverlayEE(TEfficiency *h1, TString tag1, TEfficiency *h2, TString tag2, bool logy = false);
	virtual void plotRatioOverlay3H(TH1D *h1, TString tag1, TH1D *h2, TString tag2, TH1D *h3, TString tag3, bool logy = false, double line1x = -999.99, double line2x = -999.99);
	virtual void plotOverlay3HData(TH1F *h1, TString tag1, TH1F *h2, TString tag2, TH1F *h3, TString tag3, bool logy = false, double line1x = -999.99, double line2x = -999.99);

	virtual void plotOverlay4H(TH1D *h1, TString tag1, TH1D *h2, TString tag2, TH1D *h3, TString tag3, TH1D *h4, TString tag4, bool logy = false, double line1x = -999.99, double line2x = -999.99);
	virtual void plotOverlay5H(TH1D *h1, TString tag1, TH1D *h2, TString tag2, TH1D *h3, TString tag3, TH1D *h4, TString tag4, TH1D *h5, TString tag5, bool logy = false, double line1x = -999.99, double line2x = -999.99);

/*****************************************************************************
##################| Utilities |###############################################
*****************************************************************************/
	template <class T> inline void getObjectSafe(TFile* pFile, TString name, T*& object){
		pFile->GetObject(name, object);
		if(!object){
			std::cout << name + " not found!" << std::endl;
			exit(-1);
		}
		return;
	};
	virtual TTree* getTree(TString treename, TString filename, TString subdir = "");
	virtual TH1D* drawTree1D(const char* arg, const TCut reqs, const char* histn, const int nbins, const double xmin, const double xmax, TTree* t, bool draw = false, const char* drawopt = "");
	virtual TH1D* drawTree1D(const char* arg, const TCut reqs, const char* histn, const int nbins, const double *xbins, TTree* t, bool draw = false, const char* drawopt = "");
	virtual TH2D* drawTree2D(const char* arg1, const char* arg2, const TCut reqs, const char* histn, const int nbinsx, const double xmin, const double xmax, const int nbinsy, const double ymin, const double ymax, TTree *tree, bool draw, const char* drawopt = "");
	virtual TString convertVarName(const char* var);
	virtual TString convertVarName2(const char* var);
	virtual int OptNBins(int);
	virtual TH1D* normHist(const TH1D *ihist);
	virtual TH2D* normHist(const TH2D *ihist);
	virtual TH1D* normHistBW(const TH1D *ihist, float scale = 1);
	virtual const Double_t* getBinning(const TH1D*);
	virtual void setPlottingRange(TH1D *&, float = 0.05, bool = false);
	virtual void setPlottingRange(TH1D *&, TH1D *&, float = 0.05, bool = false);
	virtual void setPlottingRange(TH1D *&, TH1D *&, TH1D *&, float = 0.05, bool = false);
	virtual void setPlottingRange(std::vector<TH1D*>&, float = 0.05, bool = false);
	virtual void getPlottingRange(float&, float&, std::vector<TH1D*> , float = 0.05, bool = false);

	virtual float getMaxYExtension(TH1*);
	virtual float getMinYExtension(TH1*);
	virtual void setZeroBinError(TH1D*);
	virtual void fillWithoutOF(TH1D *&, double, double=1.);
	
	virtual TCanvas* makeCanvas(const char*);
	virtual void printObject(TObject* o, TString name, Option_t *drawopt = "", bool logy = false);
	virtual TH1D* bookTH1D(const char*, const char*, int, double, double);
	virtual void printProgress(int, const int, TString, const int = -1);
	virtual TString numbForm(double);
	virtual int getExp(double e);
	virtual void refValues(const char* var, TH1D* h);

	virtual void getWeightedYMeanRMS(TH1D*, double&, double&);

	virtual double tailFraction(TH1D* h, double frac);
	virtual void printCheckList(const char* var, TH1D* h, const char* filename);
	virtual TString printTailFraction(const char* var, TH1D* h, double frac);
	virtual TString printAverage(const char* var, TH1D* h);
	virtual TString printRatio(const char* var, TH1D* h, double x1, double x2, double y1, double y2);
/*****************************************************************************
##################| Variables |###############################################
*****************************************************************************/
	inline virtual void setOutputDir(TString dir){ fOutputDir = Util::MakeOutputDir(dir); };
	inline virtual void setOutputSubDir(TString dir){ fOutputSubDir = Util::MakeOutputDir(dir); };
	inline virtual void setOutputFile(TString filename){ fOutputFile = Util::MakeOutputFile(fOutputDir + filename); };
	inline virtual void setCheckListFile(TString file){ fChecklistFile = file; };
	inline virtual void setGlobalTag(TString tag){ fGlobalTag = tag; };
	inline virtual void setVerbose(int v){ fVerbose = v;};

	int fVerbose;
	
	std::vector<TFile*>  fFile;
	std::vector<TString> fFileName;
	std::vector<TString> fTreeSubDirName; // Name of subdirectory inside rootfile where tree is
	std::vector<TTree*>  fTree;
	std::vector<TString> fTag;
	std::vector<TCut>    fCut;
	std::vector<double>  fNorm;

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
	TFile *fOutputFile;  // File containing output histograms

private:
	// ClassDef(AnaClass,2)
};

#endif
