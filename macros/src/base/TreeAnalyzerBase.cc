#include <stdlib.h>
#include "base/TreeAnalyzerBase.hh"
#include <TTree.h>
#include <TString.h>

using namespace std;

TreeAnalyzerBase::TreeAnalyzerBase(TTree *tree) {
	fTR = new TreeReader(tree);
	fVerbose = false;
}

TreeAnalyzerBase::~TreeAnalyzerBase(){
	if(!fTR->fChain) cout << "TreeAnalyzerBase ==> No chain!" << endl;
}

void TreeAnalyzerBase::SetOutputDir(TString dir){
	if(!dir.EndsWith("/")) dir += "/";
	fOutputDir = dir;
	// Create directory if needed
	//  >> NOTE: This function needs to be called before the booking functions!
	char cmd[100];
	sprintf(cmd,"mkdir -p %s", fOutputDir.Data());
	system(cmd);
};

// Method for looping over the tree
void TreeAnalyzerBase::Loop(){
	Long64_t nentries = fTR->GetEntries();
	cout << " total events in ntuples = " << fTR->GetEntries() << endl;
	// nentries = 10;
	for( Long64_t jentry = 0; jentry < nentries; jentry++ ){
		if( jentry%200 == 0 ) cout << ">>> Processing event # " << jentry << endl;
		fTR->GetEntry(jentry);

	}
}

// Method called before starting the event loop
void TreeAnalyzerBase::BeginJob(){}

// Method called after finishing the event loop
void TreeAnalyzerBase::EndJob(){}

void TreeAnalyzerBase::SetStyle(){
	fStyle = new TStyle("ETHStyle", "Standard Plain");
	fStyle->SetCanvasColor(0);
	fStyle->SetFrameFillColor(0);
	fStyle->SetFrameBorderMode(0);
	fStyle->SetFrameBorderSize(0);
	fStyle->SetPalette(1,0);
	fStyle->SetOptTitle(0);
	fStyle->SetOptStat(111111);
	fStyle->SetStatColor(0);
	fStyle->SetStatStyle(3001);
	fStyle->SetStatBorderSize(1);

	// Fonts
	Int_t font = 42;
	fStyle->SetStatFont(font);
	fStyle->SetTextFont(font);
	fStyle->SetLabelFont(font, "xyz");
	fStyle->SetTitleFont(font, "xyz");

	// Histograms
	fStyle->SetHistFillColor(15);
	fStyle->SetHistFillStyle(1001);
	fStyle->SetHistLineWidth(2);
	gROOT->SetStyle("ETHStyle");
	gROOT->ForceStyle();
}
