/*****************************************************************************
*     Collection of tools for producing plots for October Exercise           *
*                                                                            *
*                                                  (c) 2009 Benjamin Stieger *
*****************************************************************************/
#include "../include/AnaClass.hh"


// ClassImp(AnaClass);
using namespace std;

/****************************************************************************/
AnaClass::AnaClass(){
// Default constructor, no samples are set
	init();
}

/****************************************************************************/
AnaClass::AnaClass(const char* parfile, bool verbose){
// Explicit constructor, reading from a parameter file
	init(verbose);
	readVarNames("varnames.dat");
	loadSamples(parfile, verbose);
}

/****************************************************************************/
AnaClass::~AnaClass(){
}

/****************************************************************************/
void AnaClass::init(bool verbose){
	if(verbose) cout << "------------------------------------" << endl;
	if(verbose) cout << "Initializing AnaClass ... " << endl;
	fOutputSubDir = "";
	fFont = 42;
	fTlat = new TLatex();
	gROOT->SetStyle("Plain");
	gStyle->SetOptTitle(0);
	gStyle->SetOptStat("");
	gStyle->SetStatFont(fFont);
	gStyle->SetStatColor(0);
	gStyle->SetStatStyle(3001);
	gStyle->SetTextFont(fFont);
	fBGColor = 10; // White, use 16 for dark gray
	gStyle->SetLabelFont(fFont, "xyz");
	gStyle->SetTitleFont(fFont, "xyz");
}

/****************************************************************************/
void AnaClass::loadSamples(TString parfile, bool verbose){
	fParFile = parfile;
	for(size_t i = 0; i < 20; ++i){ // Reset subdirnames before reading parfile
		fTreeSubDirName[i] = "";
	}
	readParms(fParFile, verbose);
	cout << " Samples:" << endl;
	for(size_t i = 0; i < 20; ++i){
		if(fFileName[i] == "") continue;
		fTree[i] = getTree("Analysis", fFileName[i], fTreeSubDirName[i]);
		cout << "  Sample Index " << i << " > " << fTag[i] << endl;
	}
}

/****************************************************************************/
void AnaClass::readParms(TString filename, bool verbose){
/*		-	Reads in parameters for analysis from filename							*/
	char buffer[200];
	fParFile = filename;
	sprintf(buffer, "%s", fParFile.Data());
	ifstream IN(buffer);
	char ParName[100], StringValue[100];
	float ParValue;

	bool ok(false);
	// TString fn(fParFile.Data());

	if(verbose) cout << "------------------------------------" << endl;
	if(verbose) cout << "Parameter File  " << fParFile.Data() << endl;
	if(verbose) cout << "------------------------------------" << endl;

	while( IN.getline(buffer, 200, '\n') ){
		ok = false;
		if (buffer[0] == '#') {continue;} // Skip lines commented with '#'
		sscanf(buffer, "%s %f", ParName, &ParValue);

		// ----------------------------------------------------------
		// Numbers
		if( !strcmp(ParName, "fNorm1") ){
			fNorm[0] = double(ParValue); ok = true;
			if(verbose) cout << "fNorm1:\t\t\t" << fNorm[0] << endl;
		}
		if( !strcmp(ParName, "fNorm2") ){
			fNorm[1] = double(ParValue); ok = true;
			if(verbose) cout << "fNorm2:\t\t\t" << fNorm[1] << endl;
		}
		if( !strcmp(ParName, "fNorm3") ){
			fNorm[2] = double(ParValue); ok = true;
			if(verbose) cout << "fNorm3:\t\t\t" << fNorm[2] << endl;
		}
		if( !strcmp(ParName, "fNorm4") ){
			fNorm[3] = double(ParValue); ok = true;
			if(verbose) cout << "fNorm4:\t\t\t" << fNorm[3] << endl;
		}
		if( !strcmp(ParName, "fNorm5") ){
			fNorm[4] = double(ParValue); ok = true;
			if(verbose) cout << "fNorm5:\t\t\t" << fNorm[4] << endl;
		}
		if( !strcmp(ParName, "fNorm6") ){
			fNorm[5] = double(ParValue); ok = true;
			if(verbose) cout << "fNorm6:\t\t\t" << fNorm[5] << endl;
		}
		if( !strcmp(ParName, "fNorm7") ){
			fNorm[6] = double(ParValue); ok = true;
			if(verbose) cout << "fNorm7:\t\t\t" << fNorm[6] << endl;
		}
		if( !strcmp(ParName, "fNorm8") ){
			fNorm[7] = double(ParValue); ok = true;
			if(verbose) cout << "fNorm8:\t\t\t" << fNorm[7] << endl;
		}
		if( !strcmp(ParName, "fNorm9") ){
			fNorm[8] = double(ParValue); ok = true;
			if(verbose) cout << "fNorm9:\t\t\t" << fNorm[8] << endl;
		}
		if( !strcmp(ParName, "fNorm10") ){
			fNorm[9] = double(ParValue); ok = true;
			if(verbose) cout << "fNorm10:\t\t\t" << fNorm[9] << endl;
		}

		// ----------------------------------------------------------
		// Strings
		sscanf(buffer, "%s %s", ParName, StringValue);
		if( !strcmp(ParName, "fTag") ){
			fGlobalTag = TString(StringValue); ok = true;
			if(verbose) cout << "fTag:\t\t\t" << fGlobalTag << endl;
		}
		if( !strcmp(ParName, "fOutputDir") ){
			fOutputDir = TString(StringValue); ok = true;
			if(verbose) cout << "fOutputDir:\t\t" << fOutputDir << endl;
		}
		if( !strcmp(ParName, "fFile1Name") ){
			fFileName[0] = TString(StringValue); ok = true;
			if(verbose) cout << "fFile1Name:\t\t" << fFileName[0] << endl;
		}
		if( !strcmp(ParName, "fFile2Name") ){
			fFileName[1] = TString(StringValue); ok = true;
			if(verbose) cout << "fFile2Name:\t\t" << fFileName[1] << endl;
		}
		if( !strcmp(ParName, "fFile3Name") ){
			fFileName[2] = TString(StringValue); ok = true;
			if(verbose) cout << "fFile3Name:\t\t" << fFileName[2] << endl;
		}
		if( !strcmp(ParName, "fFile4Name") ){
			fFileName[3] = TString(StringValue); ok = true;
			if(verbose) cout << "fFile4Name:\t\t" << fFileName[3] << endl;
		}
		if( !strcmp(ParName, "fFile5Name") ){
			fFileName[4] = TString(StringValue); ok = true;
			if(verbose) cout << "fFile5Name:\t\t" << fFileName[4] << endl;
		}
		if( !strcmp(ParName, "fFile6Name") ){
			fFileName[5] = TString(StringValue); ok = true;
			if(verbose) cout << "fFile6Name:\t\t" << fFileName[5] << endl;
		}
		if( !strcmp(ParName, "fFile7Name") ){
			fFileName[6] = TString(StringValue); ok = true;
			if(verbose) cout << "fFile7Name:\t\t" << fFileName[6] << endl;
		}
		if( !strcmp(ParName, "fFile8Name") ){
			fFileName[7] = TString(StringValue); ok = true;
			if(verbose) cout << "fFile8Name:\t\t" << fFileName[7] << endl;
		}
		if( !strcmp(ParName, "fFile9Name") ){
			fFileName[8] = TString(StringValue); ok = true;
			if(verbose) cout << "fFile9Name:\t\t" << fFileName[8] << endl;
		}
		if( !strcmp(ParName, "fFile10Name") ){
			fFileName[9] = TString(StringValue); ok = true;
			if(verbose) cout << "fFile10Name:\t\t" << fFileName[9] << endl;
		}
		if( !strcmp(ParName, "fTreeSubDirName1") ){
			fTreeSubDirName[0] = TString(StringValue); ok = true;
			if(verbose) cout << "fTreeSubDirName1:\t\t" << fTreeSubDirName[0] << endl;
		}
		if( !strcmp(ParName, "fTreeSubDirName2") ){
			fTreeSubDirName[1] = TString(StringValue); ok = true;
			if(verbose) cout << "fTreeSubDirName2:\t\t" << fTreeSubDirName[1] << endl;
		}
		if( !strcmp(ParName, "fTreeSubDirName3") ){
			fTreeSubDirName[2] = TString(StringValue); ok = true;
			if(verbose) cout << "fTreeSubDirName3:\t\t" << fTreeSubDirName[2] << endl;
		}
		if( !strcmp(ParName, "fTreeSubDirName4") ){
			fTreeSubDirName[3] = TString(StringValue); ok = true;
			if(verbose) cout << "fTreeSubDirName4:\t\t" << fTreeSubDirName[3] << endl;
		}
		if( !strcmp(ParName, "fTreeSubDirName5") ){
			fTreeSubDirName[4] = TString(StringValue); ok = true;
			if(verbose) cout << "fTreeSubDirName5:\t\t" << fTreeSubDirName[4] << endl;
		}
		if( !strcmp(ParName, "fTreeSubDirName6") ){
			fTreeSubDirName[5] = TString(StringValue); ok = true;
			if(verbose) cout << "fTreeSubDirName6:\t\t" << fTreeSubDirName[5] << endl;
		}
		if( !strcmp(ParName, "fTreeSubDirName7") ){
			fTreeSubDirName[6] = TString(StringValue); ok = true;
			if(verbose) cout << "fTreeSubDirName7:\t\t" << fTreeSubDirName[6] << endl;
		}
		if( !strcmp(ParName, "fTreeSubDirName8") ){
			fTreeSubDirName[7] = TString(StringValue); ok = true;
			if(verbose) cout << "fTreeSubDirName8:\t\t" << fTreeSubDirName[7] << endl;
		}
		if( !strcmp(ParName, "fTreeSubDirName9") ){
			fTreeSubDirName[8] = TString(StringValue); ok = true;
			if(verbose) cout << "fTreeSubDirName9:\t\t" << fTreeSubDirName[8] << endl;
		}
		if( !strcmp(ParName, "fTreeSubDirName10") ){
			fTreeSubDirName[9] = TString(StringValue); ok = true;
			if(verbose) cout << "fTreeSubDirName10:\t\t" << fTreeSubDirName[9] << endl;
		}
		if( !strcmp(ParName, "fTag1") ){
			fTag[0] = TString(StringValue); ok = true;
			if(verbose) cout << "fTag1:\t\t\t" << fTag[0] << endl;
		}
		if( !strcmp(ParName, "fTag2") ){
			fTag[1] = TString(StringValue); ok = true;
			if(verbose) cout << "fTag2:\t\t\t" << fTag[1] << endl;
		}
		if( !strcmp(ParName, "fTag3") ){
			fTag[2] = TString(StringValue); ok = true;
			if(verbose) cout << "fTag3:\t\t\t" << fTag[2] << endl;
		}
		if( !strcmp(ParName, "fTag4") ){
			fTag[3] = TString(StringValue); ok = true;
			if(verbose) cout << "fTag4:\t\t\t" << fTag[3] << endl;
		}
		if( !strcmp(ParName, "fTag5") ){
			fTag[4] = TString(StringValue); ok = true;
			if(verbose) cout << "fTag5:\t\t\t" << fTag[4] << endl;
		}
		if( !strcmp(ParName, "fTag6") ){
			fTag[5] = TString(StringValue); ok = true;
			if(verbose) cout << "fTag6:\t\t\t" << fTag[5] << endl;
		}
		if( !strcmp(ParName, "fTag7") ){
			fTag[6] = TString(StringValue); ok = true;
			if(verbose) cout << "fTag7:\t\t\t" << fTag[6] << endl;
		}
		if( !strcmp(ParName, "fTag8") ){
			fTag[7] = TString(StringValue); ok = true;
			if(verbose) cout << "fTag8:\t\t\t" << fTag[7] << endl;
		}
		if( !strcmp(ParName, "fTag9") ){
			fTag[8] = TString(StringValue); ok = true;
			if(verbose) cout << "fTag9:\t\t\t" << fTag[8] << endl;
		}
		if( !strcmp(ParName, "fTag10") ){
			fTag[9] = TString(StringValue); ok = true;
			if(verbose) cout << "fTag10:\t\t\t" << fTag[9] << endl;
		}
		if( !strcmp(ParName, "fL1Cuts") ){
			fL1Cuts = TCut(StringValue); ok = true;
			if(verbose) cout << "fL1Cuts:\t\t" << fL1Cuts << endl;
		}
		if( !strcmp(ParName, "fL2Cuts") ){
			fL2Cuts = TCut(StringValue); ok = true;
			if(verbose) cout << "fL2Cuts:\t\t" << fL2Cuts << endl;
		}
		if(!ok) cout << "%% ERROR: Unknown variable " << ParName << endl;
	}
	fL1L2Cuts = fL1Cuts && fL2Cuts;
	fL1Cuts.SetName("L1Cuts");
	fL2Cuts.SetName("L2Cuts");
	fL1L2Cuts.SetName("L1L2Cuts");
	if(verbose){
		cout << "------------------------------------" << endl;
	}
}

/****************************************************************************/
void AnaClass::readVarNames(const char* filename){
// Fills the VarNameMap map from the varnames.dat file
	TString branchname, texname;
	ifstream IN(filename);
	char buff1[200], buff2[200];
	// Loop over lines of datafile
	while( IN.getline(buff1, 200, '\t') && IN.getline(buff2, 200, '\n') ){
		// Convert chararrays to TStrings
		branchname = TString(buff1); 
		texname = TString(buff2);
		texname.ReplaceAll("\t",""); // Remove tabs
		// Fill map
		fVarNameMap[branchname] = texname;
	}
}

/*****************************************************************************
###################| Produce Plots |##########################################
*****************************************************************************/

/*****************************************************************************
###################| Main Methods |###########################################
*****************************************************************************/

/****************************************************************************/
void AnaClass::plotPlotList(const char* filename, TTree *tree, TString tag){
	int sampleindex(0), nbins(0), logy(0);
	float xmin(0.), xmax(0.), x1(-999.), x2(-999.);
	ifstream IN(filename);
	// char buff1[200], buff2[200], buff3[200];
	char readbuff[200];
	char path[200], varname[200], reqbuff[200];
	TCut req;
	// Loop over lines of datafile
	while( IN.getline(readbuff, 200, '\n') ){
		if (readbuff[0] == '#') {continue;} // Skip lines commented with '#'
		int nargs = sscanf(readbuff, "%s %s %d %d %f %f %d %f %f %s", path, varname, &sampleindex, &nbins, &xmin, &xmax, &logy, &x1, &x2, reqbuff);
		if(nargs < 8){ x1 = -999; x2 = -999; }
		else if(nargs < 9) x2 = -999;
		else if(nargs < 10) req = "";
		else req = TCut(reqbuff);
		// cout << "varname=" << varname << " req=" << req << " sampleindex=" << sampleindex << " nbins=" << nbins << " xmin=" << xmin << " xmax=" << xmax << " logy=" << logy << " x1=" << x1 << " x2=" << x2 << endl;
		if(tag=="tag") tag = fTag[sampleindex];
		if(tree==NULL) tree = fTree[sampleindex];
		fOutputSubDir = tag + "/" + TString(path);
		if(!strcmp(varname, "ElID")) plotEID(req, tree, tag);
		else plotVar(varname, req, tree, tag, nbins, xmin, xmax, logy, x1, x2);
	}
}

/****************************************************************************/
void AnaClass::plotPlotList2D(const char* filename, TTree *tree, TString tag){
	int sampleindex(0), nbinsx(0), nbinsy(0), logx(0), logy(0), logz(0), mrkstl(0);
	float xmin(0.), xmax(0.), x1(-999.), x2(-999.);
	float ymin(0.), ymax(0.), y1(-999.), y2(-999.);
	ifstream IN(filename);
	// char buff1[200], buff2[200], buff3[200];
	char readbuff[200];
	char path[200], varname1[200], varname2[200], reqbuff[200], topt[200];
	TCut req;
	// Loop over lines of datafile
	while( IN.getline(readbuff, 200, '\n') ){
		if (readbuff[0] == '#') {continue;} // Skip lines commented with '#'
		int nargs = sscanf(readbuff, "%s %s %s %d %d %f %f %d %f %f %s %d %d %d %d %f %f %f %f %s", path, varname1, varname2, &sampleindex, &nbinsx, &xmin, &xmax, &nbinsy, &ymin, &ymax, topt, &mrkstl, &logx, &logy, &logz, &x1, &x2, &y1, &y2, reqbuff);
		if(nargs < 16){ x1 = -999.; x2 = -999.; y1 = -999.; y2 = -999.;}
		else if(nargs < 17){ x2 = -999.; y1 = -999.; y2 = -999.;}
		else if(nargs < 18){ y1 = -999.; y2 = -999.;}
		else if(nargs < 19) y2 = -999.;
		else if(nargs < 20) req = "";
		else req = TCut(reqbuff);
		// cout << "varname=" << varname << " req=" << req << " sampleindex=" << sampleindex << " nbins=" << nbins << " xmin=" << xmin << " xmax=" << xmax << " logy=" << logy << " x1=" << x1 << " x2=" << x2 << endl;
		if(tag=="tag") tag = fTag[sampleindex];
		if(tree==NULL) tree = fTree[sampleindex];
		fOutputSubDir = TString(path);
		plotVar2D(varname1, varname2, req, tree, tag, nbinsx, xmin, xmax, nbinsy, ymin, ymax, topt, mrkstl, logx, logy, logz, x1, x2, y1, y2);
	}
}

/****************************************************************************/
void AnaClass::plotAllBranches(TTree *tree, TString tag){
//  Plots all branches of a Tree and saves the histograms in a rootfile
	gStyle->SetOptStat(1111);
	fOutputSubDir = tag + "/" + "All/";
	TObjArray *array = tree->GetListOfBranches();
	TH1 *hists[array->GetEntries()];
	for(size_t i = 0; i < array->GetEntries(); ++i){
		const char* branchname = array->At(i)->GetName();
		TCanvas *c = makeCanvas(Form("c_%s", branchname));
		tree->Draw(branchname);
		hists[i] = (TH1*)gROOT->FindObject("htemp");
		hists[i]->SetName(Form("h_%s", branchname));
		hists[i]->SetLineWidth(2);
		hists[i]->SetFillColor(15);
		hists[i]->SetFillStyle(1001);
		TString outputname = tag + "_" + branchname;
		printPNG(c, outputname, fOutputSubDir);
		printEPS(c, outputname, fOutputSubDir);
	}
}
void AnaClass::plotAllBranches(int sampleindex){
//  Plots all branches of a Tree and saves the histograms in a rootfile
//  -- Overloaded to enable calling by the index in the parameter file
	gStyle->SetOptStat(1111);
	fOutputSubDir = fTag[sampleindex] + "/" + "AllBranches/";
	TTree *tree = fTree[sampleindex];
	TObjArray *array = tree->GetListOfBranches();
	TH1 *hists[array->GetEntries()];
	for(size_t i = 0; i < array->GetEntries(); ++i){
		const char* branchname = array->At(i)->GetName();
		TCanvas *c = makeCanvas(Form("c_%s", branchname));
		tree->Draw(branchname);
		hists[i] = (TH1*)gROOT->FindObject("htemp");
		hists[i]->SetName(Form("h_%s", branchname));
		hists[i]->SetLineWidth(2);
		hists[i]->SetFillColor(15);
		hists[i]->SetFillStyle(1001);
		TString outputname = fTag[sampleindex] + "_" + branchname;
		printPNG(c, outputname, fOutputSubDir);
		printEPS(c, outputname, fOutputSubDir);
	}
}

/****************************************************************************/
void AnaClass::plotEID(TCut req, TTree *t, TString tag){
	gStyle->SetOptStat(1111111);
	if( TH1D *h = (TH1D*)gROOT->FindObject("hfir")) h->Delete();
	TH1D *hfir = new TH1D("ElID","eID",5,0,5);
	if(!hfir){ cout << "AnaClass::plotVar() ==> Error, missing input histogram ..." << endl; return;}

	TCut tempreq = req && "ElID[0][0]>0";
	long n_tight  = t->Draw("ElID[0][0]", tempreq, "goff");
	tempreq = req && "ElID[0][1]>0";
	long n_loose  = t->Draw("ElID[0][1]", tempreq, "goff");
	tempreq = req && "ElID[0][2]>0";
	long n_rtight = t->Draw("ElID[0][2]", tempreq, "goff");
	tempreq = req && "ElID[0][3]>0";
	long n_rloose = t->Draw("ElID[0][3]", tempreq, "goff");
	tempreq = req && "ElID[0][0]==0&&ElID[0][1]==0&&ElID[0][2]==0&&ElID[0][3]==0";
	long n_none = t->Draw("ElID[0][0]", tempreq, "goff");

	hfir->SetBinContent(1, n_none);
	hfir->GetXaxis()->SetBinLabel(1, "None");
	hfir->SetBinContent(2, n_tight);
	hfir->GetXaxis()->SetBinLabel(2, "Tight");
	hfir->SetBinContent(3, n_loose);
	hfir->GetXaxis()->SetBinLabel(3, "Loose");
	hfir->SetBinContent(4, n_rtight);
	hfir->GetXaxis()->SetBinLabel(4, "RobustTight");
	hfir->SetBinContent(5, n_rloose);
	hfir->GetXaxis()->SetBinLabel(5, "RobustLoose");

	hfir->SetLineWidth(2);
	hfir->SetFillColor(15);
	hfir->SetFillStyle(1001);
	char ctitle[1000];
	sprintf(ctitle,"Plot of ElID: %s", tag.Data());
	TCanvas *col = makeCanvas(ctitle);
	col->SetFillStyle(0);
	col->SetFrameFillStyle(0);
	col->cd();
	gPad->SetFillStyle(0);

	// Determine plotting range
	double max = hfir->GetMaximum();
	max = 1.05*max;
	hfir->SetMaximum(max);
	hfir->SetMinimum(0.0);
	hfir->DrawCopy("hist");

	TString outputname = tag + "_ElID";
	printPNG(col, outputname, fOutputSubDir);
	printEPS(col, outputname, fOutputSubDir);
}

/*****************************************************************************
###################| Utilities |##############################################
*****************************************************************************/

/****************************************************************************/
void AnaClass::setOutputDir(TString dir){
	if(!dir.EndsWith("/")) dir += "/";
	fOutputDir = dir;
	// Create directory if needed
	char cmd[100];
	sprintf(cmd,"mkdir -p %s", fOutputDir.Data());
	system(cmd);
}

/****************************************************************************/
TTree* AnaClass::getTree(TString treename, TString filename, TString subdir){
	TFile *file = NULL;
	TTree *tree = NULL;
	file = TFile::Open(filename);
	if(file == NULL){
		cout << "AnaClass::getTree ==> Input file '" << filename << "' not found, breaking!" << endl;
		return NULL;
	}
	if(subdir != ""){
		TString treepath = TString(subdir) + "/" + TString(treename);
		tree = (TTree*)file->Get(treepath);
	}
	else tree = (TTree*)file->Get(treename);
	if(tree == NULL){
		if(subdir == "") cout << "AnaClass::getTree ==> Tree '" << treename << "' not found, breaking!" << endl;
		if(subdir != "") cout << "AnaClass::getTree ==> Tree '" << treename << "' in subdir '"<< subdir <<"' not found, breaking!" << endl;
		return NULL;
	}
	return tree;
}

/****************************************************************************/
TH1D* AnaClass::drawTree1D(const char* arg, const TCut reqs, const char* histn, const int nbins, const double xmin, const double xmax, TTree* tree, bool draw, const char* drawopt){
	int nbins_auto = nbins;
	if(nbins == 0)	nbins_auto = OptNBins(tree->Draw(arg, reqs, "goff"));
	TH1D* h1 = new TH1D(histn,histn,nbins_auto,xmin,xmax);
	tree->Project(histn,arg,reqs);
	if(draw) h1->Draw(drawopt);
	return h1;
}

/****************************************************************************/
TH2D* AnaClass::drawTree2D(const char* arg1, const char* arg2, const TCut reqs, const char* histn, const int nbinsx, const double xmin, const double xmax, const int nbinsy, const double ymin, const double ymax, TTree* tree, bool draw, const char* drawopt){
	char out[1000];
	TH2D* h1 = new TH2D(histn,histn,nbinsx,xmin,xmax,nbinsy,ymin,ymax);
	sprintf(out,"%s:%s",arg2,arg1);
	tree->Project(histn,out,reqs);
	if(draw) h1->Draw(drawopt);
	return h1;
}

/****************************************************************************/
void AnaClass::plotVar(const char* var, const TCut reqs, TTree *tree, TString tag, int nbins, double xmin, double xmax, bool logy, double line1x, double line2x){
	gStyle->SetOptStat(1111111);
	if( TH1D *h = (TH1D*)gROOT->FindObject("hfir")) h->Delete();
	TH1D *hfir = drawTree1D(var,reqs,convertVarName2(var),nbins,xmin,xmax,tree,false);
	if(!hfir){ cout << "AnaClass::plotVar() ==> Error, missing input histogram ..." << endl; return;}

	hfir->SetXTitle(convertVarName(var));
	hfir->SetLineWidth(2);
	hfir->SetFillColor(15);
	hfir->SetFillStyle(1001);
	char ctitle[1000];
	sprintf(ctitle,"Plot of %s: %s", var, tag.Data());
	TCanvas *col = makeCanvas(ctitle);
	col->SetFillStyle(0);
	col->SetFrameFillStyle(0);
	col->cd();
	gPad->SetFillStyle(0);
	if(logy) col->SetLogy(1);

	// Determine plotting range
	double max = hfir->GetMaximum();
	if(logy) max = 5*max;
	else max = 1.05*max;
	hfir->SetMaximum(max);

	if(!logy) hfir->SetMinimum(0.0);

	hfir->DrawCopy("hist");

	double min = hfir->GetYaxis()->GetXmin();

	TLine *l1, *l2;
	if(line1x != -999.){
		l1 = new TLine(line1x,min,line1x,max);
		l1->SetLineColor(kRed);
		l1->SetLineWidth(2);
		l1->Draw();
	}

	if(line2x != -999.){
		l2 = new TLine(line2x,min,line2x,max);
		l2->SetLineColor(kRed);
		l2->SetLineWidth(2);
		l2->Draw();
	}
	TString outputname = tag + "_" + convertVarName2(var);
	printPNG(col, outputname, fOutputSubDir);
	printEPS(col, outputname, fOutputSubDir);
}

/****************************************************************************/
void AnaClass::plotVar2D(const char* var1, const char* var2, const TCut reqs, TTree *tree, TString tag, int nbinsx, double xmin, double xmax, int nbinsy, double ymin, double ymax, Option_t *topt, int markstyle, bool logx, bool logy, bool logz, double line1x, double line2x, double line1y, double line2y){
	gStyle->SetOptStat(1111111);
	TString histn = Form("%svs%s", convertVarName2(var1).Data(), convertVarName2(var2).Data());
	TH2D *h = drawTree2D(var1, var2, reqs, histn , nbinsx, xmin, xmax, nbinsy, ymin, ymax, tree, false);
	if(!h){ cout << "AnaClass::plotVar2D() ==> Error, missing input histogram ..." << endl; return;}

	h->SetXTitle(convertVarName(var1));
	h->SetYTitle(convertVarName(var2));
	h->SetMarkerStyle(markstyle);

	char ctitle[1000];
	sprintf(ctitle,"Plot of %svs%s: %s", var1, var2, tag.Data());
	TCanvas *col = makeCanvas(ctitle);
	col->SetFillStyle(0);
	col->SetFrameFillStyle(0);
	col->cd();
	gPad->SetFillStyle(0);
	if(logx) col->SetLogx(1);
	if(logy) col->SetLogy(1);
	if(logz) col->SetLogz(1);

	// // Determine plotting range
	// double max = h->GetMaximum();
	// if(logz) max = 5*max;
	// else max = 1.05*max;
	// h->SetMaximum(max);
	// 
	// if(!logz) h->SetMinimum(0.0);

	h->DrawCopy(topt);

	double minx = h->GetXaxis()->GetXmin();
	double maxx = h->GetXaxis()->GetXmax();
	double miny = h->GetYaxis()->GetXmin();
	double maxy = h->GetYaxis()->GetXmax();

	TLine *l1, *l2, *l3, *l4;
	if(line1x != -999.){
		l1 = new TLine(line1x,miny,line1x,maxy);
		l1->SetLineColor(kRed);
		l1->SetLineWidth(2);
		l1->Draw();
	}
	
	if(line2x != -999.){
		l2 = new TLine(line2x,miny,line2x,maxy);
		l2->SetLineColor(kRed);
		l2->SetLineWidth(2);
		l2->Draw();
	}
	if(line1y != -999.){
		l3 = new TLine(minx, line1y, maxx, line1y);
		l3->SetLineColor(kRed);
		l3->SetLineWidth(2);
		l3->Draw();
	}
	
	if(line2y != -999.){
		l4 = new TLine(minx, line2y, maxx, line2y);
		l4->SetLineColor(kRed);
		l4->SetLineWidth(2);
		l4->Draw();
	}

	TString outputname = tag + "_" + convertVarName2(var1) + "-" + convertVarName2(var2);
	printPNG(col, outputname, fOutputSubDir);
	printEPS(col, outputname, fOutputSubDir);
}

/****************************************************************************/
void AnaClass::plotOverlay2T(const char* var, const TCut reqs, int index1, int index2, int nbins, double xmin, double xmax, bool logy, double line1x, double line2x){
	gStyle->SetOptStat("");
	if( TH1D *h = (TH1D*)gROOT->FindObject("hfir")) h->Delete();
	if( TH1D *h = (TH1D*)gROOT->FindObject("hsec")) h->Delete();
	TH1D *hfir = drawTree1D(var,reqs,"hfir",nbins,xmin,xmax,fTree[index1],false);
	TH1D *hsec = drawTree1D(var,reqs,"hsec",nbins,xmin,xmax,fTree[index2],false);

	if(!hfir){ cout << "AnaClass::plotOverlay2T() ==> Error missing input histogram ..." << endl; return;}
	if(!hsec){ cout << "AnaClass::plotOverlay2T() ==> Error missing input histogram ..." << endl; return;}

	hfir->SetXTitle(convertVarName(var));
	hfir->SetLineWidth(2);
	hfir->SetFillColor(15);
	hfir->SetFillStyle(1001);
	hsec->SetLineWidth(2);
	hsec->SetLineColor(kBlue);
	hsec->SetFillColor(kBlue);
	hsec->SetFillStyle(3005);
	char ctitle[1000];
	sprintf(ctitle,"Overlay of %s: %s vs %s", var, fTag[index1].Data(), fTag[index2].Data());
	TCanvas *col = makeCanvas(ctitle);
	col->SetFillStyle(0);
	col->SetFrameFillStyle(0);
	col->cd();
	gPad->SetFillStyle(0);
	if(logy) col->SetLogy(1);
	hfir->Scale(1.0/hfir->Integral());
	hfir = normHist(hfir);
	hsec = normHist(hsec);

	// Determine plotting range
	double max1 = hfir->GetMaximum();
	double max2 = hsec->GetMaximum();
	double max = (max1>max2)?max1:max2;
	if(logy) max = 5*max;
	else max = 1.05*max;
	hfir->SetMaximum(max);
	hsec->SetMaximum(max);

	if(!logy){
		hfir->SetMinimum(0.0);
		hsec->SetMinimum(0.0);
	}

	TLegend *leg = new TLegend(0.65,0.75,0.886,0.88);
	leg->AddEntry(hfir,fTag[index1],"f");
	leg->AddEntry(hsec,fTag[index2],"f");
	leg->SetFillStyle(0);
	leg->SetTextFont(fFont);

	hfir->DrawCopy("hist");
	hsec->DrawCopy("histsame");
	leg->Draw();

	double min1 = hfir->GetYaxis()->GetXmin();
	double min2 = hsec->GetYaxis()->GetXmin();
	double min = (min1<min2)?min1:min2;

	TLine *l1, *l2;
	if(line1x != -999.){
		l1 = new TLine(line1x,min,line1x,max);
		l1->SetLineColor(kRed);
		l1->SetLineWidth(2);
		l1->Draw();
	}

	if(line2x != -999.){
		l2 = new TLine(line2x,min,line2x,max);
		l2->SetLineColor(kRed);
		l2->SetLineWidth(2);
		l2->Draw();
	}

	TString outputname = fTag[index1] + "_" + fTag[index2] + "_" + convertVarName2(var);
	printPNG(col, outputname, fOutputSubDir);
	printEPS(col, outputname, fOutputSubDir);
}

/****************************************************************************/
void AnaClass::plotOverlay1T2V(const char* var1, const char* var2, const TCut reqs, int sampleindex, int nbins, double xmin, double xmax, bool logy, double line1x, double line2x){
	gStyle->SetOptStat("");
	if( TH1D *h = (TH1D*)gROOT->FindObject("hfir")) h->Delete();
	if( TH1D *h = (TH1D*)gROOT->FindObject("hsec")) h->Delete();
	TH1D *hfir = drawTree1D(var1,reqs,"hfir",nbins,xmin,xmax,fTree[sampleindex],false);
	TH1D *hsec = drawTree1D(var2,reqs,"hsec",nbins,xmin,xmax,fTree[sampleindex],false);

	if(!hfir){ cout << "AnaClass::plotOverlay1() ==> Error missing input histogram ..." << endl; return;}
	if(!hsec){ cout << "AnaClass::plotOverlay1() ==> Error missing input histogram ..." << endl; return;}

	hfir->SetXTitle(convertVarName(var1));
	hfir->SetLineWidth(2);
	hfir->SetFillColor(15);
	hfir->SetFillStyle(1001);
	hsec->SetLineWidth(2);
	hsec->SetLineColor(kRed);
	hsec->SetFillColor(kRed);
	hsec->SetFillStyle(3005);
	char ctitle[1000];
	sprintf(ctitle,"Overlay of %s vs %s", var1, var2);
	TCanvas *col = makeCanvas(ctitle);
	col->SetFillStyle(0);
	col->SetFrameFillStyle(0);
	col->cd();
	gPad->SetFillStyle(0);
	if(logy) col->SetLogy(1);
	hfir = normHist(hfir);
	hsec = normHist(hsec);

	// Determine plotting range
	double max1 = hfir->GetMaximum();
	double max2 = hsec->GetMaximum();
	double max = (max1>max2)?max1:max2;
	if(logy) max = 5*max;
	else max = 1.05*max;
	hfir->SetMaximum(max);
	hsec->SetMaximum(max);

	TLegend *leg = new TLegend(0.65,0.75,0.886,0.88);
	leg->AddEntry(hfir,convertVarName(var1),"f");
	leg->AddEntry(hsec,convertVarName(var2),"f");
	leg->SetFillStyle(0);
	leg->SetTextFont(fFont);

	hfir->DrawCopy("hist");
	hsec->DrawCopy("histsame");
	leg->Draw();

	TLine *l1, *l2;
	if(line1x != -999.){
		l1 = new TLine(line1x,0,line1x,max);
		l1->SetLineColor(kRed);
		l1->SetLineWidth(2);
		l1->Draw();
	}

	if(line2x != -999.){
		l2 = new TLine(line2x,0,line2x,max);
		l2->SetLineColor(kRed);
		l2->SetLineWidth(2);
		l2->Draw();
	}
	char out[100];
	sprintf(out, "%s_%s-%s", fTag[sampleindex].Data(), convertVarName2(var1).Data(), convertVarName2(var2).Data());
	printPNG(col, out, fOutputSubDir);
	printEPS(col, out, fOutputSubDir);
}

/****************************************************************************/
void AnaClass::plotOverlay2C(const char* var, const TCut req1, const TCut req2, int sampleindex, TString tag1, TString tag2, int nbins, double xmin, double xmax, bool logy){
/*		-	Creates an normalized overlay from a tree variable with two conditions
		-	Arguments:
			 var: tree. var to be drawn
			 req1, req2: arguments to be used
			 file: file to be used
			 nbins, xmin, xmax: specification for the histogram
			 logy toggle logarithmiy plot                                       */
	gStyle->SetOptStat("");
	if( TH1D *h = (TH1D*)gROOT->FindObject("hfir")) h->Delete();
	if( TH1D *h = (TH1D*)gROOT->FindObject("hsec")) h->Delete();
	TH1D *hfir = drawTree1D(var,req1,"hfir",nbins,xmin,xmax,fTree[sampleindex],false);
	TH1D *hsec = drawTree1D(var,req2,"hsec",nbins,xmin,xmax,fTree[sampleindex],false);

	if(!hfir){ cout << "AnaClass::plotOverlay2C() ==> Error missing input histogram ..." << endl; return;}
	if(!hsec){ cout << "AnaClass::plotOverlay2C() ==> Error missing input histogram ..." << endl; return;}

	hfir->SetXTitle(convertVarName(var));
	hfir->SetLineWidth(2);
	hfir->SetFillColor(15);
	hfir->SetFillStyle(1001);
	hsec->SetLineWidth(2);
	hsec->SetLineColor(kBlue);
	hsec->SetFillColor(kBlue);
	hsec->SetFillStyle(3005);
	TCanvas *col = makeCanvas(Form("Overlay of %s", var));
	col->cd();
	if(logy) col->SetLogy(1);
	hfir = normHist(hfir);
	hsec = normHist(hsec);
	double max1 = hfir->GetMaximum();
	double max2 = hsec->GetMaximum();
	double max = (max1>max2)?max1:max2;
	if(logy) max = 5*max;
	else max = 1.05*max;
	hfir->SetMaximum(max);
	hsec->SetMaximum(max);

	TLegend *leg = new TLegend(0.6,0.73,0.917,0.88);
	leg->AddEntry(hfir,tag1,"f");
	leg->AddEntry(hsec,tag2,"f");
	leg->SetFillColor(0);
	leg->SetTextFont(fFont);

	hfir->DrawCopy("hist");
	hsec->DrawCopy("histsame");
	leg->Draw();
	TString outputname = fTag[sampleindex] + "_" + convertVarName2(var) + "_" + tag1 + "-" + tag2;
	printPNG(col, outputname, fOutputSubDir);
	printEPS(col, outputname, fOutputSubDir);
}

/****************************************************************************/
void AnaClass::plotOverlay3T(const char* var, const TCut reqs, int index1, int index2, int index3, int nbins, double xmin, double xmax, bool logy, double line1x, double line2x){
	gStyle->SetOptStat("");
	if( TH1D *h = (TH1D*)gROOT->FindObject("h1")) h->Delete();
	if( TH1D *h = (TH1D*)gROOT->FindObject("h2")) h->Delete();
	if( TH1D *h = (TH1D*)gROOT->FindObject("h3")) h->Delete();
	TH1D *h1 = drawTree1D(var,reqs,"h1",nbins,xmin,xmax,fTree[index1],false);
	TH1D *h2 = drawTree1D(var,reqs,"h2",nbins,xmin,xmax,fTree[index2],false);
	TH1D *h3 = drawTree1D(var,reqs,"h3",nbins,xmin,xmax,fTree[index3],false);

	if(!h1 || !h2 || !h3){
		cout << "AnaClass::plotOverlay3T() ==> Error missing input histogram ..." << endl;
		return;
	}

	h1->SetXTitle(convertVarName(var));
	h1->SetLineWidth(2);
	h1->SetFillColor(15);
	h1->SetFillStyle(1001);
	h2->SetLineWidth(2);
	h2->SetLineColor(kBlue);
	h2->SetFillColor(kBlue);
	h2->SetFillStyle(3005);
	h3->SetLineWidth(2);
	h3->SetLineColor(kRed);
	h3->SetFillColor(kRed);
	h3->SetFillStyle(3004);
	char ctitle[1000];
	sprintf(ctitle,"Overlay of %s: %s vs %s vs %s", var, fTag[index1].Data(), fTag[index2].Data(), fTag[index3].Data());
	TCanvas *col = makeCanvas(ctitle);
	col->SetFillStyle(0);
	col->SetFrameFillStyle(0);
	col->cd();
	gPad->SetFillStyle(0);
	if(logy) col->SetLogy(1);
	h1 = normHist(h1);
	h2 = normHist(h2);
	h3 = normHist(h3);

	// Determine plotting range
	double max1 = h1->GetMaximum();
	double max2 = h2->GetMaximum();
	double max3 = h3->GetMaximum();
	double max12 = (max1>max2)?max1:max2;
	double max = (max12>max3)?max12:max3;
	if(logy) max = 5*max;
	else max = 1.05*max;
	h1->SetMaximum(max);
	h2->SetMaximum(max);
	h3->SetMaximum(max);
	if(!logy){
		h1->SetMinimum(0.0);
		h2->SetMinimum(0.0);
		h3->SetMinimum(0.0);
	}

	TLegend *leg = new TLegend(0.65,0.75,0.886,0.88);
	leg->AddEntry(h1,fTag[index1],"f");
	leg->AddEntry(h2,fTag[index2],"f");
	leg->AddEntry(h3,fTag[index3],"f");
	leg->SetFillStyle(0);
	leg->SetTextFont(fFont);

	h1->DrawCopy("hist");
	h2->DrawCopy("histsame");
	h3->DrawCopy("histsame");
	leg->Draw();

	double min1 = h1->GetYaxis()->GetXmin();
	double min2 = h2->GetYaxis()->GetXmin();
	double min3 = h3->GetYaxis()->GetXmin();
	double min12 = (min1<min2)?min1:min2;
	double min = (min12<min3)?min12:min3;

	TLine *l1, *l2;
	if(line1x != -999.){
		l1 = new TLine(line1x,min,line1x,max);
		l1->SetLineColor(kRed);
		l1->SetLineWidth(2);
		l1->Draw();
	}

	if(line2x != -999.){
		l2 = new TLine(line2x,min,line2x,max);
		l2->SetLineColor(kRed);
		l2->SetLineWidth(2);
		l2->Draw();
	}
	TString outputname = convertVarName2(var) + "_" + fTag[index1] + "_" + fTag[index2] + "_" + fTag[index3];
	printPNG(col, outputname, fOutputSubDir);
	printEPS(col, outputname, fOutputSubDir);
}

/****************************************************************************/
void AnaClass::plotOverlay3C(const char* var, const TCut req1, TString tag1, const TCut req2, TString tag2, const TCut req3, TString tag3, int sampleindex, int nbins, double xmin, double xmax, bool logy){
/*		-	Creates an normalized overlay from a tree variable with three
			conditions
		-	Arguments:
			 var: tree. var to be drawn
			 req1, req2, req3: arguments to be used
			 file: file to be used
			 nbins, xmin, xmax: specification for the histogram
			 logy toggle logarithmiy plot                                       */
	gStyle->SetOptStat("");
	if( TH1D *h = (TH1D*)gROOT->FindObject("hfir")) h->Delete();
	if( TH1D *h = (TH1D*)gROOT->FindObject("hsec")) h->Delete();
	if( TH1D *h = (TH1D*)gROOT->FindObject("hthr")) h->Delete();
	TH1D *hfir = drawTree1D(var,req1,"hfir",nbins,xmin,xmax,fTree[sampleindex],false);
	TH1D *hsec = drawTree1D(var,req2,"hsec",nbins,xmin,xmax,fTree[sampleindex],false);
	TH1D *hthr = drawTree1D(var,req3,"hthr",nbins,xmin,xmax,fTree[sampleindex],false);

	if(!hfir){ cout << "AnaClass::plotOverlay3C() ==> Error missing input histogram ..." << endl; return;}
	if(!hsec){ cout << "AnaClass::plotOverlay3C() ==> Error missing input histogram ..." << endl; return;}
	if(!hthr){ cout << "AnaClass::plotOverlay3C() ==> Error missing input histogram ..." << endl; return;}

	hfir->SetXTitle(convertVarName(var));
	hfir->SetLineWidth(2);
	hfir->SetFillColor(15);
	hfir->SetFillStyle(1001);
	hsec->SetLineWidth(2);
	hsec->SetLineColor(kBlue);
	hsec->SetFillColor(kBlue);
	hsec->SetFillStyle(3005);

	hthr->SetLineWidth(2);
	hthr->SetLineColor(kRed);
	hthr->SetFillColor(kRed);
	hthr->SetFillStyle(3003);

	TCanvas *col = makeCanvas(Form("Overlay of %s", var));
	col->cd();
	if(logy) col->SetLogy(1);
	hfir = normHist(hfir);
	hsec = normHist(hsec);
	hthr = normHist(hthr);
	double max1 = hfir->GetMaximum();
	double max2 = hsec->GetMaximum();
	double max3 = hthr->GetMaximum();
	double tempmax = (max1>max2)?max1:max2;
	double max = (tempmax>max3)?tempmax:max3;
	if(logy) max = 5*max;
	else max = 1.05*max;
	hfir->SetMaximum(max);
	hsec->SetMaximum(max);
	hthr->SetMaximum(max);

	TLegend *leg = new TLegend(0.7,0.73,0.917,0.88);
	leg->AddEntry(hfir,tag1,"f");
	leg->AddEntry(hsec,tag2,"f");
	leg->AddEntry(hthr,tag3,"f");
	leg->SetFillColor(0);
	leg->SetTextFont(fFont);

	hfir->DrawCopy("hist");
	hsec->DrawCopy("histsame");
	hthr->DrawCopy("histsame");
	leg->Draw();
	TString outputname = fTag[sampleindex] + "_" + convertVarName2(var) + "_" + tag1 + "-" + tag2 + "-" + tag3;
	printPNG(col, outputname, fOutputSubDir);
	printEPS(col, outputname, fOutputSubDir);
}

/****************************************************************************/
void AnaClass::plotOverlay4T(const char* var, const TCut reqs, int index1, int index2, int index3, int index4, int nbins, double xmin, double xmax, bool logy, double line1x, double line2x){
	gStyle->SetOptStat("");
	if( TH1D *h = (TH1D*)gROOT->FindObject("h1")) h->Delete();
	if( TH1D *h = (TH1D*)gROOT->FindObject("h2")) h->Delete();
	if( TH1D *h = (TH1D*)gROOT->FindObject("h3")) h->Delete();
	if( TH1D *h = (TH1D*)gROOT->FindObject("h4")) h->Delete();
	TH1D *h1 = drawTree1D(var,reqs,"h1",nbins,xmin,xmax,fTree[index1],false);
	TH1D *h2 = drawTree1D(var,reqs,"h2",nbins,xmin,xmax,fTree[index2],false);
	TH1D *h3 = drawTree1D(var,reqs,"h3",nbins,xmin,xmax,fTree[index3],false);
	TH1D *h4 = drawTree1D(var,reqs,"h4",nbins,xmin,xmax,fTree[index4],false);

	if(!h1 || !h2 || !h3 || !h4){
		cout << "AnaClass::plotOverlay3T() ==> Error missing input histogram ..." << endl;
		return;
	}

	h1->SetXTitle(convertVarName(var));
	h1->SetLineWidth(2);
	h1->SetFillColor(15);
	h1->SetFillStyle(1001);
	h2->SetLineWidth(2);
	h2->SetLineColor(kBlue);
	h2->SetFillColor(kBlue);
	h2->SetFillStyle(3005);
	h3->SetLineWidth(2);
	h3->SetLineColor(kRed);
	h3->SetFillColor(kRed);
	h3->SetFillStyle(3004);
	h4->SetLineWidth(2);
	h4->SetLineColor(kGreen);
	h4->SetFillColor(kGreen);
	h4->SetFillStyle(3003);
	char ctitle[1000];
	sprintf(ctitle,"Overlay of %s: %s vs %s vs %s vs %s", var, fTag[index1].Data(), fTag[index2].Data(), fTag[index3].Data(), fTag[index4].Data());
	TCanvas *col = makeCanvas(ctitle);
	col->SetFillStyle(0);
	col->SetFrameFillStyle(0);
	col->cd();
	gPad->SetFillStyle(0);
	if(logy) col->SetLogy(1);
	h1 = normHist(h1);
	h2 = normHist(h2);
	h3 = normHist(h3);
	h4 = normHist(h4);

	// Determine plotting range
	double max1 = h1->GetMaximum();
	double max2 = h2->GetMaximum();
	double max3 = h3->GetMaximum();
	double max4 = h4->GetMaximum();
	double max12 = (max1>max2)?max1:max2;
	double max123 = (max12>max3)?max12:max3;
	double max = (max123>max4)?max123:max4;
	if(logy) max = 5*max;
	else max = 1.05*max;
	h1->SetMaximum(max);
	h2->SetMaximum(max);
	h3->SetMaximum(max);
	h4->SetMaximum(max);

	if(!logy){
		h1->SetMinimum(0.0);
		h2->SetMinimum(0.0);
		h3->SetMinimum(0.0);
		h4->SetMinimum(0.0);
	}

	TLegend *leg = new TLegend(0.65,0.72,0.886,0.88);
	leg->AddEntry(h1,fTag[index1],"f");
	leg->AddEntry(h2,fTag[index2],"f");
	leg->AddEntry(h3,fTag[index3],"f");
	leg->AddEntry(h4,fTag[index4],"f");
	leg->SetFillStyle(0);
	leg->SetTextFont(fFont);

	h1->DrawCopy("hist");
	h2->DrawCopy("histsame");
	h3->DrawCopy("histsame");
	h4->DrawCopy("histsame");
	leg->Draw();

	double min1 = h1->GetYaxis()->GetXmin();
	double min2 = h2->GetYaxis()->GetXmin();
	double min3 = h3->GetYaxis()->GetXmin();
	double min4 = h4->GetYaxis()->GetXmin();
	double min12 = (min1<min2)?min1:min2;
	double min123 = (min12<min3)?min12:min3;
	double min = (min123<min4)?min123:min4;

	TLine *l1, *l2;
	if(line1x != -999.){
		l1 = new TLine(line1x,min,line1x,max);
		l1->SetLineColor(kRed);
		l1->SetLineWidth(2);
		l1->Draw();
	}

	if(line2x != -999.){
		l2 = new TLine(line2x,min,line2x,max);
		l2->SetLineColor(kRed);
		l2->SetLineWidth(2);
		l2->Draw();
	}
	TString outputname = convertVarName2(var) + "_" + fTag[index1] + "-" + fTag[index2] + "-" + fTag[index3] + "-" + fTag[index4];
	printPNG(col, outputname, fOutputSubDir);
	printEPS(col, outputname, fOutputSubDir);
}

/****************************************************************************/
TString AnaClass::convertVarName(const char* var){
/*  - Converts the name of a tree variable into a label for the x-axis      */
	TString temp = TString(var);
	TString outp = fVarNameMap[temp];
	if(outp == ""){
		cout << "AnaClass::convertVarName() ==> Variable " << temp << " unknown..." << endl;
		return temp;
	}
	else return outp;
}

/****************************************************************************/
TString AnaClass::convertVarName2(const char* var){
/*  - Removes bracket signs from filenames to enable saving as .eps files  */
	TString out_str = TString(var);
	out_str.ReplaceAll("(","");
	out_str.ReplaceAll(")","");
	out_str.ReplaceAll("[","");
	out_str.ReplaceAll("]","");
	out_str.ReplaceAll("/","");
	return out_str;
}

/****************************************************************************/
int AnaClass::OptNBins(int nentries){
// Returns an 'optimal' number of bins for some number of entries
	if(nentries < 100)  return 10;
	if(nentries < 2000) return 60;
	if(nentries < 5000) return 80;
	if(nentries < 10000) return 100;
	else return 200;
}

/****************************************************************************/
void AnaClass::printPNG(TCanvas *cin, const char* name, const char* dir){
/*		-	Prints a ROOT TCanvas Object to a .png file
			name is the bare output filename, e.g. "fit_4_8",
			dir is the output directory (inside the overall output dir.)
			including the last "/", e.g. "fits/SG/"									 */
	// Create sub directories if needed
	char cmd[100];
	sprintf(cmd,"mkdir -p %s%s",fOutputDir.Data(),dir);
	system(cmd);

	TString filen = TString(name);
	TString direc = TString(dir);
	TString file = fOutputDir;
	file += direc;
	file += filen;
	file += ".png";
	cin->Print(file,"png");
}

/****************************************************************************/
void AnaClass::printEPS(TCanvas *cin, const char* name, const char* dir){
/*		-	Prints a ROOT TCanvas Object to a .eps file
			name is the bare output filename, e.g. "fit_4_8",
			dir is the output directory (inside the overall output dir.)
			including the last "/", e.g. "fits/SG/"									 */
	// Create sub directories if needed
	char cmd[100];
	sprintf(cmd,"mkdir -p %s%seps/",fOutputDir.Data(),dir);
	system(cmd);

	TString filen = TString(name);
	TString direc = TString(dir);
	TString file = fOutputDir;
	file += direc;
	file += "eps/";
	file += filen;
	file += ".eps";
	cin->SaveAs(file);
}

/****************************************************************************/
TH1D* AnaClass::normHist(const TH1D *ihist){
/*		-	Normalizes a histogram (incl. errors) to unit integral
			I.e. divides each entry by the integral									*/
	TH1D *ohist = new TH1D(*ihist);
	double scale = ihist->Integral();
	for( int i = 0; i < ihist->GetNbinsX()+2; i++ ){
		ohist->SetBinContent(i,ihist->GetBinContent(i)/scale);
		ohist->SetBinError(i,ihist->GetBinError(i)/scale);
	}
	return ohist;
}

/****************************************************************************/
TH2D* AnaClass::normHist(const TH2D *ihist){
/*		-	Normalizes a histogram (incl. errors) to unit integral
			I.e. divides each entry by the integral									*/
	TH2D *ohist = new TH2D(*ihist);
	double scale = ihist->Integral();
	for( int i = 0; i < ihist->GetNbinsX()+2; i++ ){
		for( int j = 0; j < ihist->GetNbinsY()+2; j++ ){
			ohist->SetBinContent(i,j,ihist->GetBinContent(i,j)/scale);
			ohist->SetBinError(i,j,ihist->GetBinError(i,j)/scale);
		}
	}
	return ohist;
}

/****************************************************************************/
TCanvas* AnaClass::makeCanvas(const char* name){
/*		-	Creates a TCanvas at a somewhat random position							 */
	int x = (int)gRandom->Uniform(10,200);
	int y = (int)gRandom->Uniform(0,70);
	TCanvas* c1 = new TCanvas(name, name, x, y, 900, 700);
	return c1;
}

/****************************************************************************/
TString AnaClass::numbForm(double n){
/*		-	Converts a number into a TString 											 */
	int expo = getExp(n);
	TString s(Form("1.e%i", expo));
	double scale = s.Atof();
	TString out("-");
	if(isnan(n))       out.Form("NaN");
	else if(expo < -2) out.Form("%1.2f E%i", n/scale, expo);
	else if(expo < -1) out.Form("%4.4f", n);
	else if(expo < 0 ) out.Form("%4.3f", n);
	else if(expo < 1 ) out.Form("%4.2f", n);
	else if(expo < 2 ) out.Form("%4.1f", n);
	else if(expo < 3 ) out.Form("%4.0f", n);
	else if(expo < 4 ) out.Form("%4.0f", n);
	else if(expo < 5 ) out.Form("%4.0f", n);
	else if(expo < 6 ) out.Form("%4.0f", n);
	// else if(expo < 7 ) out.Form("%4.0f", n);
	else               out.Form("%1.2f E%i", n/scale, expo);
	return out;
}

/****************************************************************************/
int AnaClass::getExp(double e){
/*		-	Returns exponent of input double												 */
	int expo(0);
	TString s(Form("%1.8e", e));
	if(e<0) s.Replace(0,11,"", 0);
	if(e>=0) s.Replace(0,10,"", 0);
	s.ReplaceAll("e", "");
	expo = s.Atoi();
	return expo;
}
