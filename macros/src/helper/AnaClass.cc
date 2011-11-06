/*****************************************************************************
*     Collection of tools for producing plots for October Exercise           *
*                                                                            *
*                                                  (c) 2009 Benjamin Stieger *
*****************************************************************************/
#include "helper/AnaClass.hh"
#include "helper/Utilities.hh"

// ClassImp(AnaClass);
using namespace std;

//____________________________________________________________________________
AnaClass::AnaClass(){
// Default constructor, no samples are set
	init();
}

//____________________________________________________________________________
AnaClass::AnaClass(const char* parfile, bool verbose){
// Explicit constructor, reading from a parameter file
	init(verbose);
	readVarNames("varnames.dat");
	loadSamples(parfile, verbose);
}

//____________________________________________________________________________
AnaClass::~AnaClass(){}

//____________________________________________________________________________
void AnaClass::init(bool verbose){
	if(verbose) cout << "------------------------------------" << endl;
	if(verbose) cout << "Initializing AnaClass ... " << endl;
	fOutputSubDir = "";
	fChecklistFile = "checklist.txt";
	Util::SetStyle();
}

//____________________________________________________________________________
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

//____________________________________________________________________________
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

//____________________________________________________________________________
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

//____________________________________________________________________________
void AnaClass::plotPlotList(const char* filename, TTree *tree, TString tag, TCut cut, TFile* file){
	int twod(0), sampleindex(0), nbinsx(0), nbinsy(0), logx(0), logy(0), logz(0), mrkstl(0);
	float xmin(0.), xmax(0.), x1(-999.), x2(-999.);
	float ymin(0.), ymax(0.), y1(-999.), y2(-999.);
	ifstream IN(filename);
	char readbuff[200];
	char path[200], var1name[200], var2name[200], reqbuff[200], topt[200];
	TCut req;
	TString temp = fOutputSubDir;
	fChecklistFile = TString(filename);
	fChecklistFile.ReplaceAll(".dat", "");
	fChecklistFile = fChecklistFile + "_checklist.txt";
	// Remove existing checklist file:
	char cmd[100];
	TString checklistfile = fOutputDir + fChecklistFile;
	sprintf(cmd,"rm -f %s", checklistfile.Data());
	system(cmd);

	// Loop over lines of datafile
	while( IN.getline(readbuff, 200, '\n') ){
		if( readbuff[0] == '#' ) continue; // Skip lines commented with '#'
		if( readbuff[0] != '2' ){ // 1D plotting
			int nargs = sscanf(readbuff, "%s %s %d %d %f %f %d %f %f %s", path, var1name, &sampleindex, &nbinsx, &xmin, &xmax, &logy, &x1, &x2, reqbuff);
			if(nargs < 7){ logy = 0; x1 = -999; x2 = -999; req=""; }
			else if(nargs < 8){ x1 = -999; x2 = -999; req=""; }
			else if(nargs < 9){ x2 = -999; req=""; }
			else if(nargs < 10) req = "";
			else req = TCut(reqbuff);
			if(tag=="tag") tag = fTag[sampleindex];
			if(tree==NULL) tree = fTree[sampleindex];
			fOutputSubDir = TString(path);
			if(!strcmp(var1name, "ElID")) plotEID(req, tree, tag);
			else plotVar(var1name, req&&cut, tree, tag, nbinsx, xmin, xmax, "ofilename", logy, x1, x2, file );
		}
		if( readbuff[0] == '2' ){ // 2D plotting
			int nargs = sscanf(readbuff, "%d %s %s %s %d %d %f %f %d %f %f %s %d %d %d %d %f %f %f %f %s", &twod, path, var1name, var2name, &sampleindex, &nbinsx, &xmin, &xmax, &nbinsy, &ymin, &ymax, topt, &mrkstl, &logx, &logy, &logz, &x1, &x2, &y1, &y2, reqbuff);
			if(     nargs < 13){ mrkstl = 1; logx = 0; logy = 0; logz = 0; x1 = -999.; x2 = -999.; y1 = -999.; y2 = -999.;}
			else if(nargs < 14){ logx = 0; logy = 0; logz = 0; x1 = -999.; x2 = -999.; y1 = -999.; y2 = -999.;}
			else if(nargs < 15){ logy = 0; logz = 0; x1 = -999.; x2 = -999.; y1 = -999.; y2 = -999.;}
			else if(nargs < 16){ logz = 0; x1 = -999.; x2 = -999.; y1 = -999.; y2 = -999.;}
			else if(nargs < 17){ x1 = -999.; x2 = -999.; y1 = -999.; y2 = -999.;}
			else if(nargs < 18){ x2 = -999.; y1 = -999.; y2 = -999.;}
			else if(nargs < 19){ y1 = -999.; y2 = -999.;}
			else if(nargs < 20) y2 = -999.;
			else if(nargs < 21) req = "";
			else req = TCut(reqbuff);
			if(tag=="tag") tag = fTag[sampleindex];
			if(tree==NULL) tree = fTree[sampleindex];
			fOutputSubDir = TString(path);
			plotVar2D(var1name, var2name, req&&cut, tree, tag, nbinsx, xmin, xmax, nbinsy, ymin, ymax, topt, mrkstl, logx, logy, logz, x1, x2, y1, y2, file );
		}
	}
	fOutputSubDir = temp;
}

//____________________________________________________________________________
void AnaClass::plotPlotList2D(const char* filename, TTree *tree, TString tag, TFile* file){
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
		if(tag=="tag") tag = fTag[sampleindex];
		if(tree==NULL) tree = fTree[sampleindex];
		fOutputSubDir = TString(path);
		plotVar2D(varname1, varname2, req, tree, tag, nbinsx, xmin, xmax, nbinsy, ymin, ymax, topt, mrkstl, logx, logy, logz, x1, x2, y1, y2, file );
	}
}

//____________________________________________________________________________
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
		Util::Print(c, outputname, fOutputSubDir);
	}
}

//____________________________________________________________________________
void AnaClass::plotAllBranches(int sampleindex){
//  Plots all branches of a Tree and saves the histograms in a rootfile
//  -- Overloaded to enable calling by the index in the parameter file
	gStyle->SetOptStat(1111);
	fOutputSubDir = fTag[sampleindex] + "/" + "AllBranches/";
	TTree *tree = fTree[sampleindex];
	TString tag = fTag[sampleindex];

	plotAllBranches(tree, tag);
}

//____________________________________________________________________________
void AnaClass::plotEID(TCut req, TTree *t, TString tag, TFile* file){
	gStyle->SetOptStat(1111111);
	if( TH1D *h = (TH1D*)gROOT->FindObject("hfir")) h->Delete();


		// The electron IDs to plot (add any and increments nIDs)
	string IDs[] = { "Loose", "Tight", "RobustTight", "RobustLoose" };
	int nIDs = 4;
	gROOT->cd(); // Store this in root directory
	TH1D* hfir = new TH1D("ElID","eID",nIDs+1,0.,nIDs+1.0);      

		// Use TTreeFormulae and loop over tree entries
	TTreeFormula* formulae[nIDs+1];
	TCut none(req);
	for ( int i=0; i<nIDs; ++i ) {
		char cut[256]; sprintf(cut,"ElID%s[0]>0",IDs[i].c_str());
		TCut tempreq = req&&cut;
		none *= !TCut(cut);
		formulae[i] = new TTreeFormula(IDs[i].c_str(),tempreq,t);
	}
	formulae[nIDs] = new TTreeFormula("None",none,t);
	cout << "Plotting eID: this can take a while... " << flush;
		// Have to do some gymnastics with the TChain
	Long_t treenumber = t->GetTreeNumber(); 
	for( int ientry=0; ientry<t->GetEntriesFast(); ++ientry ) 
	{ 
		if ( t->LoadTree(ientry) < 0 ) break; 
		if (t->GetTreeNumber() != treenumber) {
			for ( int i=0; i<=nIDs; ++i ) formulae[i]->UpdateFormulaLeaves(); 
			treenumber = t->GetTreeNumber(); 
		} 
		for ( int i=0; i<=nIDs; ++i )
			if ( formulae[i]->GetNdata()>0 && formulae[i]->EvalInstance() )
			hfir->Fill(i);
	}
	cout << "done." << endl;

		// Set bin labels to ID names
	for ( int i=0; i<nIDs; ++i )
		hfir->GetXaxis()->SetBinLabel(i+1, IDs[i].c_str());
	hfir->GetXaxis()->SetBinLabel(nIDs+1, "None");

		// Styling
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

	TString outputname = tag + "ElID";
	fOutputSubDir = "ElID";
	Util::Print(col, outputname, fOutputDir + fOutputSubDir, file);

	delete col;
	delete hfir;
}

//____________________________________________________________________________
void AnaClass::plotVar(const char* var, const TCut reqs, TTree *tree, TString tag, int nbins, double xmin, double xmax, TString ofilename, bool logy, double line1x, double line2x, TFile* file){
	gStyle->SetOptStat(1111111);
	if( TH1D *h = (TH1D*)gROOT->FindObject("hfir")) h->Delete();
	TH1D *hfir = drawTree1D(var,reqs,convertVarName2(var),nbins,xmin,xmax,tree,false);
	if(!hfir){ cout << "AnaClass::plotVar() ==> Error, missing input histogram ..." << endl; return;}
	if(ofilename == "ofilename") ofilename = convertVarName2(var);
	char ctitle[1000];
	sprintf(ctitle,"Plot of %s: %s", var, tag.Data());
	TCanvas *col = makeCanvas(ctitle);
	col->cd();
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
		l1 = new TLine(line1x, min, line1x, max);
		l1->SetLineColor(kRed);
		l1->SetLineWidth(2);
		l1->Draw();
	}

	if(line2x != -999.){
		l2 = new TLine(line2x, min, line2x, max);
		l2->SetLineColor(kRed);
		l2->SetLineWidth(2);
		l2->Draw();
	}

	refValues(var, hfir);

	TString outputname = tag + ofilename;
	Util::Print(col, outputname, fOutputDir + fOutputSubDir, file );

	printCheckList(var, hfir, fOutputDir + fChecklistFile);

	delete col;
	delete hfir;
}

//____________________________________________________________________________
void AnaClass::plotVar2D(const char* var1, const char* var2, const TCut reqs, TTree *tree, TString tag, int nbinsx, double xmin, double xmax, int nbinsy, double ymin, double ymax, Option_t *topt, int markstyle, bool logx, bool logy, bool logz, double line1x, double line2x, double line1y, double line2y, TFile* file ){
	gStyle->SetOptStat(1111111);
	TString histn = Form("%svs%s", convertVarName2(var1).Data(), convertVarName2(var2).Data());
	TH2D *h = drawTree2D(var1, var2, reqs, histn , nbinsx, xmin, xmax, nbinsy, ymin, ymax, tree, false);
	if(!h){ cout << "AnaClass::plotVar2D() ==> Error, missing input histogram ..." << endl; return;}

	h->GetYaxis()->SetTitleOffset(1.25);
	h->SetMarkerStyle(markstyle);

	char ctitle[1000];
	sprintf(ctitle,"Plot of %svs%s: %s", var1, var2, tag.Data());
	TCanvas *col = makeCanvas(ctitle);
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

	TString outputname = tag + convertVarName2(var1) + "-" + convertVarName2(var2);
	Util::Print(col, outputname, fOutputDir + fOutputSubDir, file );

	delete col;
	delete h;
}

//____________________________________________________________________________
void AnaClass::plotOverlay2H(TH1D *h1, TString tag1, TH1D *h2, TString tag2, bool logy, double line1x, double line2x){
	gStyle->SetOptStat("");
	h1->SetLineWidth(2);
	h1->SetFillColor(15);
	h1->SetFillStyle(1001);
	h2->SetLineWidth(2);
	h2->SetLineColor(kBlue);
	h2->SetFillColor(kBlue);
	h2->SetFillStyle(3005);

	char canvtitle[100], canvname[100];
	sprintf(canvtitle,"%s vs %s", h1->GetName(), h2->GetName());
	sprintf(canvname,"%s:%s", h1->GetName(), h2->GetName());
	TCanvas *col = new TCanvas(canvname, canvtitle, 0, 0, 900, 700);
	col->SetFillStyle(0);
	col->SetFrameFillStyle(0);
	col->cd();
	gPad->SetFillStyle(0);
	if(logy) col->SetLogy(1);
	h1->Sumw2();
	h2->Sumw2();

	h1->Scale(1.0/h1->Integral());
	h2->Scale(1.0/h2->Integral());

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

	TLegend *leg = new TLegend(0.65,0.75,0.886,0.88);
	leg->AddEntry(h1, tag1,"f");
	leg->AddEntry(h2, tag2,"f");
	leg->SetFillStyle(0);
	leg->SetTextFont(42);
	leg->SetBorderSize(1);

	h1->DrawCopy("hist");
	h2->DrawCopy("histsame");
	leg->Draw();

	double min1 = h1->GetYaxis()->GetXmin();
	double min2 = h2->GetYaxis()->GetXmin();
	double min  = (min1<min2)?min1:min2;

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
	gPad->RedrawAxis();
	TString outputname = TString(h1->GetName()) + "_" + TString(h2->GetName());
	Util::PrintNoEPS(col, outputname, fOutputDir, fOutputFile);
	// Util::Print(col, outputname, fOutputDir, fOutputFile);
}
void AnaClass::plotOverlay2HNorm(TH1D *h1, TString tag1, TH1D *h2, TString tag2, bool logy, double line1x, double line2x){
	h1->Sumw2();
	h2->Sumw2();

	h1->Scale(1.0/h1->Integral());
	h2->Scale(1.0/h2->Integral());
	plotOverlay2H(h1, tag1, h2, tag2, logy, line1x, line2x);
}
void AnaClass::plotOverlay3H(TH1D *h1, TString tag1, TH1D *h2, TString tag2, TH1D *h3, TString tag3, bool logy, double line1x, double line2x){
	gStyle->SetOptStat("");
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

	char canvtitle[100], canvname[100];
	sprintf(canvtitle,"%s vs %s vs %s", h1->GetName(), h2->GetName(), h3->GetName());
	sprintf(canvname,"%s:%s:%s", h1->GetName(), h2->GetName(), h3->GetName());
	TCanvas *col = new TCanvas(canvname, canvtitle, 0, 0, 900, 700);
	col->SetFillStyle(0);
	col->SetFrameFillStyle(0);
	col->cd();
	gPad->SetFillStyle(0);
	if(logy) col->SetLogy(1);
	h1->Sumw2();
	h2->Sumw2();
	h3->Sumw2();

	h1->Scale(1.0/h1->Integral());
	h2->Scale(1.0/h2->Integral());
	h3->Scale(1.0/h3->Integral());

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

	TLegend *leg = new TLegend(0.65,0.69,0.886,0.88);
	leg->AddEntry(h1, tag1,"f");
	leg->AddEntry(h2, tag2,"f");
	leg->AddEntry(h3, tag3,"f");
	leg->SetFillStyle(0);
	leg->SetTextFont(42);
	leg->SetBorderSize(1);

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
	gPad->RedrawAxis();
	TString outputname = TString(h1->GetName()) + "_" + TString(h2->GetName()) + "_" + TString(h3->GetName());
	// Util::Print(col, outputname, fOutputDir, fOutputFile);
	Util::PrintNoEPS(col, outputname, fOutputDir, fOutputFile);
}
void AnaClass::plotOverlay2T(const char* var, const TCut reqs, int index1, int index2, int nbins, double xmin, double xmax, bool logy, double line1x, double line2x){
	gStyle->SetOptStat("");
	if( TH1D *h = (TH1D*)gROOT->FindObject("hfir")) h->Delete();
	if( TH1D *h = (TH1D*)gROOT->FindObject("hsec")) h->Delete();
	TH1D *hfir = drawTree1D(var,reqs,"hfir",nbins,xmin,xmax,fTree[index1],false);
	TH1D *hsec = drawTree1D(var,reqs,"hsec",nbins,xmin,xmax,fTree[index2],false);

	if(!hfir){ cout << "AnaClass::plotOverlay2T() ==> Error missing input histogram ..." << endl; return;}
	if(!hsec){ cout << "AnaClass::plotOverlay2T() ==> Error missing input histogram ..." << endl; return;}

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
	Util::Print(col, outputname, fOutputSubDir);
}
void AnaClass::plotOverlay1T2V(const char* var1, const char* var2, const TCut reqs, int sampleindex, int nbins, double xmin, double xmax, bool logy, double line1x, double line2x){
	gStyle->SetOptStat("");
	if( TH1D *h = (TH1D*)gROOT->FindObject("hfir")) h->Delete();
	if( TH1D *h = (TH1D*)gROOT->FindObject("hsec")) h->Delete();
	TH1D *hfir = drawTree1D(var1,reqs,"hfir",nbins,xmin,xmax,fTree[sampleindex],false);
	TH1D *hsec = drawTree1D(var2,reqs,"hsec",nbins,xmin,xmax,fTree[sampleindex],false);

	if(!hfir){ cout << "AnaClass::plotOverlay1() ==> Error missing input histogram ..." << endl; return;}
	if(!hsec){ cout << "AnaClass::plotOverlay1() ==> Error missing input histogram ..." << endl; return;}

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
	Util::Print(col, out, fOutputSubDir);
}
void AnaClass::plotOverlay2C(const char* var, const TCut req1, const TCut req2, int sampleindex, TString tag1, TString tag2, int nbins, double xmin, double xmax, bool logy){
// Creates a normalized overlay from a tree variable with two conditions
// -   Arguments:
//     var: tree. var to be drawn
//     req1, req2: arguments to be used
//     file: file to be used
//     nbins, xmin, xmax: specification for the histogram
//     logy toggle logarithmiy plot                                      

	gStyle->SetOptStat("");
	if( TH1D *h = (TH1D*)gROOT->FindObject("hfir")) h->Delete();
	if( TH1D *h = (TH1D*)gROOT->FindObject("hsec")) h->Delete();
	TH1D *hfir = drawTree1D(var,req1,"hfir",nbins,xmin,xmax,fTree[sampleindex],false);
	TH1D *hsec = drawTree1D(var,req2,"hsec",nbins,xmin,xmax,fTree[sampleindex],false);

	if(!hfir){ cout << "AnaClass::plotOverlay2C() ==> Error missing input histogram ..." << endl; return;}
	if(!hsec){ cout << "AnaClass::plotOverlay2C() ==> Error missing input histogram ..." << endl; return;}

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
	Util::Print(col, outputname, fOutputSubDir);
}
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
	Util::Print(col, outputname, fOutputSubDir);
}
void AnaClass::plotOverlay3C(const char* var, const TCut req1, TString tag1, const TCut req2, TString tag2, const TCut req3, TString tag3, int sampleindex, int nbins, double xmin, double xmax, bool logy){
// - Creates a normalized overlay from a tree variable with three conditions
// -	Arguments:
//    var: tree. var to be drawn
//    req1, req2, req3: arguments to be used
//    file: file to be used
//    nbins, xmin, xmax: specification for the histogram
//   logy toggle logarithmiy plotgStyle->SetOptStat("");
	if( TH1D *h = (TH1D*)gROOT->FindObject("hfir")) h->Delete();
	if( TH1D *h = (TH1D*)gROOT->FindObject("hsec")) h->Delete();
	if( TH1D *h = (TH1D*)gROOT->FindObject("hthr")) h->Delete();
	TH1D *hfir = drawTree1D(var,req1,"hfir",nbins,xmin,xmax,fTree[sampleindex],false);
	TH1D *hsec = drawTree1D(var,req2,"hsec",nbins,xmin,xmax,fTree[sampleindex],false);
	TH1D *hthr = drawTree1D(var,req3,"hthr",nbins,xmin,xmax,fTree[sampleindex],false);

	if(!hfir){ cout << "AnaClass::plotOverlay3C() ==> Error missing input histogram ..." << endl; return;}
	if(!hsec){ cout << "AnaClass::plotOverlay3C() ==> Error missing input histogram ..." << endl; return;}
	if(!hthr){ cout << "AnaClass::plotOverlay3C() ==> Error missing input histogram ..." << endl; return;}

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
	Util::Print(col, outputname, fOutputSubDir);
}
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
		cout << "AnaClass::plotOverlay4T() ==> Error missing input histogram ..." << endl;
		return;
	}

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
	Util::Print(col, outputname, fOutputSubDir);
}

//____________________________________________________________________________
void AnaClass::plotPredOverlay2H(TH1D *h1, TString tag1, TH1D *h2, TString tag2, bool logy, double line1x, double line2x){
	gStyle->SetOptStat("");
	h1->SetFillColor(kBlue);
	h1->SetLineColor(kBlue);
	h1->SetLineWidth(2);
	h1->SetLineStyle(0);
	h1->SetFillStyle(3004);

	h2->SetLineWidth(2);
	h2->SetLineColor(kRed);
	h2->SetFillColor(kRed);
	// h2->SetFillStyle(3005);

	char canvtitle[100], canvname[100];
	sprintf(canvtitle,"%s vs %s", h1->GetName(), h2->GetName());
	sprintf(canvname,"%s:%s", h1->GetName(), h2->GetName());
	TCanvas *col = new TCanvas(canvname, canvtitle, 0, 0, 900, 700);
	col->SetFillStyle(0);
	col->SetFrameFillStyle(0);
	col->cd();
	gPad->SetFillStyle(0);
	if(logy) col->SetLogy(1);

	setPlottingRange(h1, h2, 0.05, logy);

	// TLegend *leg = new TLegend(0.65,0.15,0.886,0.28); // Lower right
	TLegend *leg = new TLegend(0.65,0.75,0.886,0.88); // Upper right
	leg->AddEntry(h1, tag1,"f");
	leg->AddEntry(h2, tag2,"f");
	leg->SetFillStyle(0);
	leg->SetTextFont(42);
	leg->SetBorderSize(0);

	h1->DrawCopy("E2");
	h2->DrawCopy("PE1same");
	leg->Draw();

	double min1 = h1->GetYaxis()->GetXmin();
	double min2 = h2->GetYaxis()->GetXmin();
	double min  = (min1<min2)?min1:min2;

	double max1 = h1->GetYaxis()->GetXmax();
	double max2 = h2->GetYaxis()->GetXmax();
	double max  = (max1>max2)?max1:max2;

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
	gPad->RedrawAxis();
	TString outputname = TString(h1->GetName()) + "_" + TString(h2->GetName());
	if(logy) outputname += "_log";
	Util::PrintNoEPS(col, outputname, fOutputDir, fOutputFile);
}
void AnaClass::plotPredOverlay2HWithRatio(TH1D *hist1, TString tag1, TH1D *hist2, TString tag2, bool logy, bool ratio, double line1x, double line2x){
	TH1D *h1 = new TH1D(*hist1);
	TH1D *h2 = new TH1D(*hist2);
	
	gStyle->SetOptStat("");
	TH1D *h_ratio = new TH1D("h_ratio", "Ratio histogram", h1->GetNbinsX(), getBinning(h1));
	Int_t color = 4;
	h1->SetFillColor(color);
	h1->SetLineColor(color);
	h1->SetLineWidth(2);
	h1->SetLineStyle(0);
	h1->SetFillStyle(3004);

	h2->SetLineWidth(2);
	h2->SetFillColor(kRed);
	h2->SetLineColor(kRed);
	h2->SetFillColor(kRed);
	// h2->SetFillStyle(3005);

	float border = 0.3;
	float scale = (1-border)/border;
	h_ratio->SetXTitle(h1->GetXaxis()->GetTitle());
	h_ratio->SetYTitle("Obs./Pred.");
	h_ratio->GetXaxis()->SetTitleSize(scale * h1->GetXaxis()->GetTitleSize());
	h_ratio->GetYaxis()->SetTitleSize(scale * h1->GetYaxis()->GetTitleSize());
	h_ratio->GetXaxis()->SetLabelSize(scale * h1->GetXaxis()->GetLabelSize());
	h_ratio->GetYaxis()->SetLabelSize(scale * h1->GetYaxis()->GetLabelSize());
	h_ratio->GetXaxis()->SetTickLength(scale * h1->GetXaxis()->GetTickLength());
	h_ratio->GetYaxis()->SetTickLength(h1->GetYaxis()->GetTickLength());
	
	h_ratio->SetFillStyle(3004);
	h_ratio->SetLineWidth(2);
	h_ratio->SetFillColor(color);
	h_ratio->SetLineColor(color);

	char canvtitle[100], canvname[100];
	sprintf(canvtitle,"%s vs %s", h1->GetName(), h2->GetName());
	sprintf(canvname,"%s:%s", h1->GetName(), h2->GetName());
	TCanvas *col = new TCanvas(canvname, canvtitle, 0, 0, 900, 900);
	col->cd();
	// col->SetFillStyle(0);
	// col->SetFrameFillStyle(0);
	// gPad->SetFillStyle(0);

	TPad *p_plot  = new TPad("plotpad",  "Pad containing the overlay plot", 0.00, border, 1.00, 1.00, 0, 0);
	p_plot->SetBottomMargin(0);
	p_plot->Draw();
	TPad *p_ratio = new TPad("ratiopad", "Pad containing the ratio",        0.00, 0.00, 1.00, border, 0, 0);
	p_ratio->SetTopMargin(0);
	p_ratio->SetBottomMargin(0.35);
	p_ratio->Draw();

	p_plot->cd();
	if(logy) p_plot->SetLogy(1);

	setPlottingRange(h1, h2, 0.05, logy);
	
	// TLegend *leg = new TLegend(0.65,0.15,0.886,0.28); // Lower right
	TLegend *leg = new TLegend(0.65,0.75,0.886,0.88); // Upper right
	leg->AddEntry(h1, tag1,"f");
	leg->AddEntry(h2, tag2,"f");
	leg->SetFillStyle(0);
	leg->SetTextFont(42);
	leg->SetBorderSize(0);

	h1->DrawCopy("E2");
	h2->DrawCopy("PE1same");
	leg->Draw();

	float axmin1 = h1->GetYaxis()->GetXmin();
	float axmin2 = h2->GetYaxis()->GetXmin();
	float axmin  = (axmin1<axmin2)?axmin1:axmin2;
	
	float axmax1 = h1->GetYaxis()->GetXmax();
	float axmax2 = h2->GetYaxis()->GetXmax();
	float axmax  = (axmax1>axmax2)?axmax1:axmax2;
	
	TLine *l1, *l2;
	if(line1x != -999.){
		l1 = new TLine(line1x,axmin,line1x,axmax);
		l1->SetLineColor(kRed);
		l1->SetLineWidth(2);
		l1->Draw();
	}

	if(line2x != -999.){
		l2 = new TLine(line2x,axmin,line2x,axmax);
		l2->SetLineColor(kRed);
		l2->SetLineWidth(2);
		l2->Draw();
	}
	gPad->RedrawAxis();
	p_plot->Draw();
	// p_plot->Update();

	p_ratio->cd();
	if(ratio) h_ratio->Divide(h2,h1);
	else      h_ratio->Divide(h1,h2);
	h_ratio->GetYaxis()->SetTitleOffset(h1->GetYaxis()->GetTitleOffset());
	h_ratio->GetYaxis()->SetTitle("Ratio");
	setPlottingRange(h_ratio, 0.3);
	h_ratio->DrawCopy("E2 ");
	TLine *l3 = new TLine(h1->GetXaxis()->GetXmin(), 1.00, h1->GetXaxis()->GetXmax(), 1.00);
	l3->SetLineWidth(2);
	l3->SetLineStyle(7);
	l3->Draw();
	gPad->RedrawAxis();
	p_ratio->Draw();

	col->Update();

	TString outputname = TString(h1->GetName()) + "_" + TString(h2->GetName());
	if(logy) outputname += "_log";
	Util::PrintNoEPS(col, outputname, fOutputDir + fOutputSubDir, fOutputFile);
	// delete h1;
	// delete h2;
}
void AnaClass::plotPredOverlay2HWithRatio(THStack *stack, TH1D *hist1, TString tag1, TH1D *hist2, TString tag2, bool logy, bool ratio, double line1x, double line2x){
	THStack *h1s = new THStack(*stack);
	TH1D *h1 = new TH1D(*hist1);
	TH1D *h2 = new TH1D(*hist2);
	
	gStyle->SetOptStat("");
	TH1D *h_ratio = new TH1D("h_ratio", "Ratio histogram", h2->GetNbinsX(), getBinning(h2));
	for(size_t i = 0; i < h1s->GetHists()->GetSize(); ++i){
		TH1D *hist = (TH1D*)h1s->GetHists()->At(i);
		hist->SetLineWidth(1);
		hist->SetFillStyle(1);
		// hist->SetFillStyle(3002);
	}

	Int_t color1 = 1;
	h1->SetLineWidth(2);
	h1->SetLineColor(color1);
	h1->SetFillColor(color1);
	h1->SetFillStyle(3001);

	Int_t color2 = 2;
	h2->SetLineWidth(2.5);
	h2->SetFillColor(color2);
	h2->SetLineColor(color2);
	h2->SetMarkerColor(color2);
	h2->SetMarkerSize(1.5);
	// h2->SetFillStyle(3005);

	float border = 0.3;
	float scale = (1-border)/border;
	h_ratio->SetXTitle(h2->GetXaxis()->GetTitle());
	h_ratio->SetYTitle("Obs./Pred.");
	h_ratio->GetXaxis()->SetTitleSize(scale * h2->GetXaxis()->GetTitleSize());
	h_ratio->GetYaxis()->SetTitleSize(scale * h2->GetYaxis()->GetTitleSize());
	h_ratio->GetXaxis()->SetLabelSize(scale * h2->GetXaxis()->GetLabelSize());
	h_ratio->GetYaxis()->SetLabelSize(scale * h2->GetYaxis()->GetLabelSize());
	h_ratio->GetXaxis()->SetTickLength(scale * h2->GetXaxis()->GetTickLength());
	h_ratio->GetYaxis()->SetTickLength(h2->GetYaxis()->GetTickLength());
	
	Int_t color = 1;
	// h_ratio->SetFillStyle(3004);
	h_ratio->SetLineWidth(2);
	h_ratio->SetFillColor(color);
	h_ratio->SetLineColor(color);
	h_ratio->SetMarkerColor(color);
	// h_ratio->SetMarkerSize(1.5);
	h_ratio->SetMarkerStyle(20);

	char canvtitle[100], canvname[100];
	sprintf(canvtitle,"%s vs %s", h1s->GetName(), h2->GetName());
	sprintf(canvname,"%s:%s", h1s->GetName(), h2->GetName());
	TCanvas *col = new TCanvas(canvname, canvtitle, 0, 0, 900, 900);
	col->cd();
	// col->SetFillStyle(0);
	// col->SetFrameFillStyle(0);
	// gPad->SetFillStyle(0);

	TPad *p_plot  = new TPad("plotpad",  "Pad containing the overlay plot", 0.00, border, 1.00, 1.00, 0, 0);
	p_plot->SetBottomMargin(0);
	p_plot->Draw();
	TPad *p_ratio = new TPad("ratiopad", "Pad containing the ratio",        0.00, 0.00, 1.00, border, 0, 0);
	p_ratio->SetTopMargin(0);
	p_ratio->SetBottomMargin(0.35);
	p_ratio->Draw();

	p_plot->cd();
	if(logy) p_plot->SetLogy(1);
	
	// TLegend *leg = new TLegend(0.65,0.15,0.886,0.28); // Lower right
	TLegend *leg = new TLegend(0.62,0.50,0.99,0.88); // Upper right
	leg->AddEntry(h2, tag2,"f");
	for(size_t i = 0; i < h1s->GetHists()->GetSize(); ++i)leg->AddEntry((TH1D*)h1s->GetHists()->At(i), ((TH1D*)h1s->GetHists()->At(i))->GetName(), "f");

	leg->SetFillStyle(0);
	leg->SetTextFont(42);
	leg->SetBorderSize(0);

	setPlottingRange(h1, h2, 0.05, logy);

	h1s->Draw("hist");
	// h1->DrawCopy("E2 same");
	h2->DrawCopy("PE X0 same");
	leg->Draw();

	float axmin1 = h1->GetYaxis()->GetXmin();
	float axmin2 = h2->GetYaxis()->GetXmin();
	float axmin  = (axmin1<axmin2)?axmin1:axmin2;
	
	float axmax1 = h1->GetYaxis()->GetXmax();
	float axmax2 = h2->GetYaxis()->GetXmax();
	float axmax  = (axmax1>axmax2)?axmax1:axmax2;
	
	TLine *l1, *l2;
	if(line1x != -999.){
		l1 = new TLine(line1x,axmin,line1x,axmax);
		l1->SetLineColor(kRed);
		l1->SetLineWidth(2);
		l1->Draw();
	}

	if(line2x != -999.){
		l2 = new TLine(line2x,axmin,line2x,axmax);
		l2->SetLineColor(kRed);
		l2->SetLineWidth(2);
		l2->Draw();
	}
	gPad->RedrawAxis();
	p_plot->Draw();

	p_ratio->cd();
	if(ratio) h_ratio->Divide(h2,h1);
	else      h_ratio->Divide(h1,h2);
	h_ratio->GetYaxis()->SetTitleOffset(h1->GetYaxis()->GetTitleOffset());
	h_ratio->GetYaxis()->SetTitle("Ratio");
	setPlottingRange(h_ratio, 0.2);
	h_ratio->DrawCopy("PE X0");
	TLine *l3 = new TLine(h1->GetXaxis()->GetXmin(), 1.00, h1->GetXaxis()->GetXmax(), 1.00);
	l3->SetLineWidth(2);
	l3->SetLineStyle(7);
	l3->Draw();
	gPad->RedrawAxis();	
	p_ratio->Draw();

	col->Update();

	TString outputname = TString(h1s->GetName()) + "_" + TString(h2->GetName());
	if(logy) outputname += "_log";
	Util::PrintNoEPS(col, outputname, fOutputDir + fOutputSubDir, fOutputFile);
}
void AnaClass::plotPredOverlay3HWithRatio(TH1D *hist1, TString tag1, TH1D *hist2, TString tag2, TH1D *hist3, TString tag3, bool logy, bool ratio, double line1x, double line2x){
	TH1D *h1 = new TH1D(*hist1);
	TH1D *h2 = new TH1D(*hist2);
	TH1D *h3 = new TH1D(*hist3);
	gStyle->SetOptStat("");
	TH1D *h_ratio1 = new TH1D("h_ratio1", "Ratio histogram 1", h1->GetNbinsX(), getBinning(h1));
	TH1D *h_ratio2 = new TH1D("h_ratio2", "Ratio histogram 2", h1->GetNbinsX(), getBinning(h1));
	Int_t color1 = 4;
	Int_t fillstyle1 = 3004;
	Int_t color2 = 8;
	Int_t fillstyle2 = 3005;
	h1->SetLineWidth(2);
	h1->SetFillColor(kRed);
	h1->SetLineColor(kRed);
	h1->SetFillColor(kRed);
	// h1->SetFillStyle(3005);

	h2->SetFillColor(color1);
	h2->SetLineColor(color1);
	h2->SetLineWidth(2);
	h2->SetLineStyle(0);
	h2->SetFillStyle(fillstyle1);

	h3->SetFillColor(color2);
	h3->SetLineColor(color2);
	h3->SetLineWidth(2);
	h3->SetLineStyle(0);
	h3->SetFillStyle(fillstyle2);

	float border = 0.3;
	float scale = (1-border)/border;
	h_ratio1->SetXTitle(h1->GetXaxis()->GetTitle());
	h_ratio1->SetYTitle("Obs./Pred.");
	h_ratio1->GetXaxis()->SetTitleSize(scale * h1->GetXaxis()->GetTitleSize());
	h_ratio1->GetYaxis()->SetTitleSize(scale * h1->GetYaxis()->GetTitleSize());
	h_ratio1->GetXaxis()->SetLabelSize(scale * h1->GetXaxis()->GetLabelSize());
	h_ratio1->GetYaxis()->SetLabelSize(scale * h1->GetYaxis()->GetLabelSize());
	h_ratio1->GetXaxis()->SetTickLength(scale * h1->GetXaxis()->GetTickLength());
	h_ratio1->GetYaxis()->SetTickLength(h1->GetYaxis()->GetTickLength());
	
	h_ratio2->SetXTitle(h1->GetXaxis()->GetTitle());
	h_ratio2->SetYTitle("Obs./Pred.");
	h_ratio2->GetXaxis()->SetTitleSize(scale * h1->GetXaxis()->GetTitleSize());
	h_ratio2->GetYaxis()->SetTitleSize(scale * h1->GetYaxis()->GetTitleSize());
	h_ratio2->GetXaxis()->SetLabelSize(scale * h1->GetXaxis()->GetLabelSize());
	h_ratio2->GetYaxis()->SetLabelSize(scale * h1->GetYaxis()->GetLabelSize());
	h_ratio2->GetXaxis()->SetTickLength(scale * h1->GetXaxis()->GetTickLength());
	h_ratio2->GetYaxis()->SetTickLength(h1->GetYaxis()->GetTickLength());
	
	h_ratio1->SetFillStyle(fillstyle1);
	h_ratio1->SetLineWidth(2);
	h_ratio1->SetFillColor(color1);
	h_ratio1->SetLineColor(color1);

	h_ratio2->SetFillStyle(fillstyle2);
	h_ratio2->SetLineWidth(2);
	h_ratio2->SetFillColor(color2);
	h_ratio2->SetLineColor(color2);

	char canvtitle[100], canvname[100];
	sprintf(canvtitle,"%s vs %s vs %s", h1->GetName(), h2->GetName(), h3->GetName());
	sprintf(canvname,"%s:%s:%s", h1->GetName(), h2->GetName(), h3->GetName());
	TCanvas *col = new TCanvas(canvname, canvtitle, 0, 0, 900, 900);
	col->cd();
	// col->SetFillStyle(0);
	// col->SetFrameFillStyle(0);
	// gPad->SetFillStyle(0);

	TPad *p_plot  = new TPad("plotpad",  "Pad containing the overlay plot", 0.00, border, 1.00, 1.00, 0, 0);
	p_plot->SetBottomMargin(0);
	p_plot->Draw();
	TPad *p_ratio = new TPad("ratiopad", "Pad containing the ratio",        0.00, 0.00, 1.00, border, 0, 0);
	p_ratio->SetTopMargin(0);
	p_ratio->SetBottomMargin(0.35);
	p_ratio->Draw();

	p_plot->cd();
	if(logy) p_plot->SetLogy(1);

	setPlottingRange(h1, h2, h3, 0.05, logy);
	
	// TLegend *leg = new TLegend(0.65,0.15,0.886,0.28); // Lower right
	// TLegend *leg = new TLegend(0.55,0.75,0.886,0.88); // Upper right
	TLegend *leg = new TLegend(0.45,0.75,0.886,0.88); // Upper right
	leg->AddEntry(h1, tag1,"f");
	leg->AddEntry(h2, tag2,"f");
	leg->AddEntry(h3, tag3,"f");
	leg->SetFillStyle(0);
	leg->SetTextFont(42);
	leg->SetBorderSize(0);

	h2->DrawCopy("E2");
	h3->DrawCopy("E2, same");
	h1->DrawCopy("PE1, same");
	leg->Draw();

	float axmin1  = h1->GetYaxis()->GetXmin();
	float axmin2  = h2->GetYaxis()->GetXmin();
	float axmin3  = h3->GetYaxis()->GetXmin();
	float axmin12 = (axmin1<axmin2)?axmin1:axmin2;
	float axmin   = (axmin12<axmin3)?axmin12:axmin3;
	
	float axmax1  = h1->GetYaxis()->GetXmax();
	float axmax2  = h2->GetYaxis()->GetXmax();
	float axmax3  = h3->GetYaxis()->GetXmax();
	float axmax12 = (axmax1>axmax2)?axmax1:axmax2;
	float axmax   = (axmax12>axmax3)?axmax12:axmax3;
	
	TLine *l1, *l2;
	if(line1x != -999.){
		l1 = new TLine(line1x,axmin,line1x,axmax);
		l1->SetLineColor(kRed);
		l1->SetLineWidth(2);
		l1->Draw();
	}

	if(line2x != -999.){
		l2 = new TLine(line2x,axmin,line2x,axmax);
		l2->SetLineColor(kRed);
		l2->SetLineWidth(2);
		l2->Draw();
	}
	gPad->RedrawAxis();
	p_plot->Draw();
	// p_plot->Update();

	p_ratio->cd();
	if(ratio){ h_ratio1->Divide(h1,h2); h_ratio2->Divide(h1,h3); }
	else{      h_ratio1->Divide(h2,h1); h_ratio2->Divide(h3,h1); }
	setPlottingRange(h_ratio1, h_ratio2, 0.2);
	h_ratio1->GetYaxis()->SetTitleOffset(h1->GetYaxis()->GetTitleOffset());
	h_ratio1->GetYaxis()->SetTitle("Ratio");
	h_ratio1->DrawCopy("E2 ");
	h_ratio2->GetYaxis()->SetTitleOffset(h1->GetYaxis()->GetTitleOffset());
	h_ratio2->GetYaxis()->SetTitle("Ratio");
	h_ratio2->DrawCopy("E2 same");
	TLine *l3 = new TLine(h1->GetXaxis()->GetXmin(), 1.00, h1->GetXaxis()->GetXmax(), 1.00);
	l3->SetLineWidth(2);
	l3->SetLineStyle(7);
	l3->Draw();
	gPad->RedrawAxis();
	p_ratio->Draw();

	col->Update();

	TString outputname = TString(h1->GetName()) + "_" + TString(h2->GetName()) + "_" + TString(h3->GetName());
	if(logy) outputname += "_log";
	Util::PrintNoEPS(col, outputname, fOutputDir + fOutputSubDir, fOutputFile);
	delete h1;
	delete h2;
	delete h3;
}

//____________________________________________________________________________
void AnaClass::plotOverlay3HData(TH1F *h1, TString tag1, TH1F *h2, TString tag2, TH1F *h3, TString tag3, bool logy, double line1x, double line2x){
	gStyle->SetOptStat("");

	h1->SetLineWidth(2);
	h1->SetLineColor(kBlack);
	h1->SetFillColor(kBlack);
	h1->SetMarkerColor(kBlack);
	h1->SetMarkerStyle(8);
	h1->SetMarkerSize(1.2);

	h2->SetLineWidth(2);
	h2->SetLineColor(kBlue);
	h2->SetFillColor(kBlue);
	h2->SetMarkerColor(kBlue);
	h2->SetMarkerStyle(8);
	h2->SetMarkerSize(1.2);

	h3->SetLineWidth(2);
	h3->SetLineColor(kRed);
	h3->SetFillColor(kRed);
	h3->SetMarkerColor(kRed);
	h3->SetMarkerStyle(8);
	h3->SetMarkerSize(1.2);

	char canvtitle[100], canvname[100];
	sprintf(canvtitle,"%s vs %s vs %s", h1->GetName(), h2->GetName(), h3->GetName());
	sprintf(canvname,"%s:%s:%s", h1->GetName(), h2->GetName(), h3->GetName());
	TCanvas *col = new TCanvas(canvname, canvtitle, 0, 0, 900, 700);
	col->SetFillStyle(0);
	col->SetFrameFillStyle(0);
	col->cd();
	gPad->SetFillStyle(0);
	if(logy) col->SetLogy(1);
	h1->Sumw2();
	h2->Sumw2();
	h3->Sumw2();

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

	TLegend *leg = new TLegend(0.15,0.70,0.35,0.88);
	// TLegend *leg = new TLegend(0.65,0.69,0.886,0.88);
	leg->AddEntry(h1, tag1,"f");
	leg->AddEntry(h2, tag2,"f");
	leg->AddEntry(h3, tag3,"f");
	leg->SetFillStyle(0);
	leg->SetTextFont(42);
	leg->SetBorderSize(0);

	h1->DrawCopy("PE1");
	h2->DrawCopy("PE1, same");
	h3->DrawCopy("PE1, same");
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

	TLatex *lat = new TLatex();
	lat->SetNDC(kTRUE);
	lat->SetTextColor(kBlack);
	lat->SetTextSize(0.04);
	lat->DrawLatex(0.58,0.85, "2.67 pb^{ -1} at  #sqrt{s} = 7 TeV");

	gPad->RedrawAxis();
	TString outputname = TString(h1->GetName()) + "_" + TString(h2->GetName()) + "_" + TString(h3->GetName());
	// Util::Print(col, outputname, fOutputDir, fOutputFile);
	Util::PrintNoEPS(col, outputname, fOutputDir, fOutputFile);
}
void AnaClass::plotOverlay4H(TH1D *h1in, TString tag1, TH1D *h2in, TString tag2, TH1D *h3in, TString tag3, TH1D *h4in, TString tag4, bool logy, double line1x, double line2x){
	gStyle->SetOptStat("");

	TH1D *h1 = new TH1D(*h1in);
	TH1D *h2 = new TH1D(*h2in);
	TH1D *h3 = new TH1D(*h3in);
	TH1D *h4 = new TH1D(*h4in);

	char canvtitle[100], canvname[100];
	sprintf(canvtitle,"%s vs %s vs %s vs %s", h1->GetName(), h2->GetName(), h3->GetName(), h4->GetName());
	sprintf(canvname,"%s:%s:%s:%s", h1->GetName(), h2->GetName(), h3->GetName(), h4->GetName());
	TCanvas *col = new TCanvas(canvname, canvtitle, 0, 0, 900, 700);
	col->SetFillStyle(0);
	col->SetFrameFillStyle(0);
	col->cd();
	gPad->SetFillStyle(0);
	if(logy) col->SetLogy(1);

	vector<TH1D*> hists;
	hists.push_back(h1);
	hists.push_back(h2);
	hists.push_back(h3);
	hists.push_back(h4);
	setPlottingRange(hists, 0.05, logy);
	
	// TLegend *leg = new TLegend(0.15,0.70,0.35,0.88);
	// TLegend *leg = new TLegend(0.65,0.69,0.886,0.88);
	TLegend *leg = new TLegend(0.65,0.66,0.886,0.88);
	leg->AddEntry(h1, tag1,"f");
	leg->AddEntry(h2, tag2,"f");
	leg->AddEntry(h3, tag3,"f");
	leg->AddEntry(h4, tag4,"f");
	leg->SetFillStyle(0);
	leg->SetTextFont(42);
	leg->SetBorderSize(0);

	h1->DrawCopy("PE2");
	h2->DrawCopy("PE1, same");
	h3->DrawCopy("PE1, same");
	h4->DrawCopy("PE1, same");
	leg->Draw();

	double max1 = h1->GetYaxis()->GetXmax();
	double max2 = h2->GetYaxis()->GetXmax();
	double max3 = h3->GetYaxis()->GetXmax();
	double max4 = h4->GetYaxis()->GetXmax();
	double max12 = (max1>max2)?max1:max2;
	double max34 = (max3>max4)?max3:max4;
	double max = (max12>max34)?max12:max34;

	double min1 = h1->GetYaxis()->GetXmin();
	double min2 = h2->GetYaxis()->GetXmin();
	double min3 = h3->GetYaxis()->GetXmin();
	double min4 = h4->GetYaxis()->GetXmin();
	double min12 = (min1<min2)?min1:min2;
	double min34 = (min3<min4)?min3:min4;
	double min = (min12<min34)?min12:min34;

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

	TLine *l3 = new TLine(h1->GetXaxis()->GetXmin(), 0.00, h1->GetXaxis()->GetXmax(), 0.00);
	l3->SetLineWidth(2);
	l3->SetLineStyle(7);
	l3->Draw();

	// TLatex *lat = new TLatex();
	// lat->SetNDC(kTRUE);
	// lat->SetTextColor(kBlack);
	// lat->SetTextSize(0.04);
	// lat->DrawLatex(0.58,0.85, "2.67 pb^{ -1} at  #sqrt{s} = 7 TeV");

	gPad->RedrawAxis();
	TString outputname = TString(h1->GetName()) + "_" + TString(h2->GetName()) + "_" + TString(h3->GetName()) + "_" + TString(h4->GetName());
	// Util::Print(col, outputname, fOutputDir, fOutputFile);
	Util::PrintNoEPS(col, outputname, fOutputDir + fOutputSubDir, fOutputFile);
}
void AnaClass::plotOverlay5H(TH1D *h1in, TString tag1, TH1D *h2in, TString tag2, TH1D *h3in, TString tag3, TH1D *h4in, TString tag4, TH1D *h5in, TString tag5, bool logy, double line1x, double line2x){
	gStyle->SetOptStat("");

	TH1D *h1 = new TH1D(*h1in);
	TH1D *h2 = new TH1D(*h2in);
	TH1D *h3 = new TH1D(*h3in);
	TH1D *h4 = new TH1D(*h4in);
	TH1D *h5 = new TH1D(*h5in);

	char canvtitle[100], canvname[100];
	sprintf(canvtitle,"%s vs %s vs %s vs %s vs %s", h1->GetName(), h2->GetName(), h3->GetName(), h4->GetName(), h5->GetName());
	sprintf(canvname,"%s:%s:%s:%s:%s", h1->GetName(), h2->GetName(), h3->GetName(), h4->GetName(), h5->GetName());
	TCanvas *col = new TCanvas(canvname, canvtitle, 0, 0, 900, 700);
	col->SetFillStyle(0);
	col->SetFrameFillStyle(0);
	col->cd();
	gPad->SetFillStyle(0);
	if(logy) col->SetLogy(1);

	vector<TH1D*> hists;
	hists.push_back(h1);
	hists.push_back(h2);
	hists.push_back(h3);
	hists.push_back(h4);
	hists.push_back(h5);
	setPlottingRange(hists, 0.05, logy);
	
	// TLegend *leg = new TLegend(0.15,0.70,0.35,0.88);
	// TLegend *leg = new TLegend(0.65,0.69,0.886,0.88);
	TLegend *leg = new TLegend(0.65,0.66,0.886,0.88);
	leg->AddEntry(h1, tag1,"f");
	leg->AddEntry(h2, tag2,"f");
	leg->AddEntry(h3, tag3,"f");
	leg->AddEntry(h4, tag4,"f");
	leg->AddEntry(h5, tag5,"f");
	leg->SetFillStyle(0);
	leg->SetTextFont(42);
	leg->SetBorderSize(0);

	h1->DrawCopy("PE2");
	h2->DrawCopy("PE1, same");
	h3->DrawCopy("PE1, same");
	h4->DrawCopy("PE1, same");
	h5->DrawCopy("PE1, same");
	leg->Draw();

	double max1 = h1->GetYaxis()->GetXmax();
	double max2 = h2->GetYaxis()->GetXmax();
	double max3 = h3->GetYaxis()->GetXmax();
	double max4 = h4->GetYaxis()->GetXmax();
	double max5 = h5->GetYaxis()->GetXmax();
	double max12   = std::max(max1, max2);
	double max34   = std::max(max3, max4);
	double max1234 = std::max(max12, max34);
	double max     = std::max(max1234, max5);

	double min1 = h1->GetYaxis()->GetXmin();
	double min2 = h2->GetYaxis()->GetXmin();
	double min3 = h3->GetYaxis()->GetXmin();
	double min4 = h4->GetYaxis()->GetXmin();
	double min5 = h5->GetYaxis()->GetXmin();
	double min12   = std::min(min1, min2);
	double min34   = std::min(min3, min4);
	double min1234 = std::min(min12, min34);
	double min     = std::min(min1234, min5);

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

	TLine *l3 = new TLine(h1->GetXaxis()->GetXmin(), 0.00, h1->GetXaxis()->GetXmax(), 0.00);
	l3->SetLineWidth(2);
	l3->SetLineStyle(7);
	l3->Draw();

	// TLatex *lat = new TLatex();
	// lat->SetNDC(kTRUE);
	// lat->SetTextColor(kBlack);
	// lat->SetTextSize(0.04);
	// lat->DrawLatex(0.58,0.85, "2.67 pb^{ -1} at  #sqrt{s} = 7 TeV");

	gPad->RedrawAxis();
	TString outputname = TString(h1->GetName()) + "_" + TString(h2->GetName()) + "_" + TString(h3->GetName()) + "_" + TString(h4->GetName()) + "_" + TString(h5->GetName());
	// Util::Print(col, outputname, fOutputDir, fOutputFile);
	Util::PrintNoEPS(col, outputname, fOutputDir + fOutputSubDir, fOutputFile);
}

//____________________________________________________________________________
void AnaClass::plotRatioOverlay2H(TH1D *h1in, TString tag1, TH1D *h2in, TString tag2, bool logy, double line1x, double line2x){
	gStyle->SetOptStat("");
	TH1D *h1 = new TH1D(*h1in);
	TH1D *h2 = new TH1D(*h2in);

	char canvtitle[100], canvname[100];
	sprintf(canvtitle,"%s vs %s", h1->GetName(), h2->GetName());
	sprintf(canvname,"%s:%s", h1->GetName(), h2->GetName());
	TCanvas *col = new TCanvas(canvname, canvtitle, 0, 0, 900, 700);
	col->SetFillStyle(0);
	col->SetFrameFillStyle(0);
	col->cd();
	gPad->SetFillStyle(0);
	gPad->SetGridy(1);
	if(logy) col->SetLogy(1);

	// setPlottingRange(h1, h2, 0.05, logy);

	// TLegend *leg = new TLegend(0.65,0.45,0.886,0.58);
	TLegend *leg = new TLegend(0.60,0.15,0.886,0.28); // Lower right
	// TLegend *leg = new TLegend(0.60,0.75,0.886,0.88);
	leg->AddEntry(h1, tag1,"f");
	leg->AddEntry(h2, tag2,"f");
	leg->SetFillStyle(0);
	leg->SetTextFont(42);
	leg->SetBorderSize(0);

	h1->DrawCopy("PE1");
	h2->DrawCopy("PE1,same");
	leg->Draw();

	double min1 = h1->GetYaxis()->GetXmin();
	double min2 = h2->GetYaxis()->GetXmin();
	double min  = (min1<min2)?min1:min2;
	double max1 = h1->GetYaxis()->GetXmax();
	double max2 = h2->GetYaxis()->GetXmax();
	double max  = (max1>max2)?max1:max2;

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
	gPad->RedrawAxis();
	TString outputname = TString(h1->GetName()) + "_" + TString(h2->GetName());
	Util::PrintNoEPS(col, outputname, fOutputDir + fOutputSubDir, fOutputFile);
	// Util::Print(col, outputname, fOutputDir + fOutputSubDir, fOutputFile);
}
void AnaClass::plotRatioOverlay3H(TH1D *h1in, TString tag1, TH1D *h2in, TString tag2, TH1D *h3in, TString tag3, bool logy, double line1x, double line2x){
	gStyle->SetOptStat("");
	TH1D *h1 = new TH1D(*h1in);
	TH1D *h2 = new TH1D(*h2in);
	TH1D *h3 = new TH1D(*h3in);

	char canvtitle[100], canvname[100];
	sprintf(canvtitle,"%s vs %s vs %s", h1->GetName(), h2->GetName(), h3->GetName());
	sprintf(canvname,"%s:%s:%s", h1->GetName(), h2->GetName(), h3->GetName());
	TCanvas *col = new TCanvas(canvname, canvtitle, 0, 0, 900, 700);
	col->SetFillStyle(0);
	col->SetFrameFillStyle(0);
	col->cd();
	gPad->SetFillStyle(0);
	gPad->SetGridy(1);
	if(logy) col->SetLogy(1);

	setPlottingRange(h1, h2, h3, 0.05, logy);

	// TLegend *leg = new TLegend(0.65,0.45,0.886,0.58);
	TLegend *leg = new TLegend(0.15,0.70,0.55,0.88);
	leg->AddEntry(h1, tag1,"f");
	leg->AddEntry(h2, tag2,"f");
	leg->AddEntry(h3, tag3,"f");
	leg->SetFillStyle(0);
	leg->SetTextFont(42);
	leg->SetBorderSize(0);

	h1->DrawCopy("PE1");
	h2->DrawCopy("PE1,same");
	h3->DrawCopy("PE1,same");
	leg->Draw();

	double min1  = h1->GetYaxis()->GetXmin();
	double min2  = h2->GetYaxis()->GetXmin();
	double min3  = h3->GetYaxis()->GetXmin();
	double min12 = (min1<min2)?min1:min2;
	double min   = (min12<min3)?min12:min3;
	double max1  = h1->GetYaxis()->GetXmax();
	double max2  = h2->GetYaxis()->GetXmax();
	double max3  = h3->GetYaxis()->GetXmax();
	double max12 = (max1>max2)?max1:max2;
	double max   = (max12>max3)?max12:max3;

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
	gPad->RedrawAxis();
	TString outputname = TString(h1->GetName()) + "_" + TString(h2->GetName()) + "_" + TString(h3->GetName());
	Util::PrintNoEPS(col, outputname, fOutputDir + fOutputSubDir, fOutputFile);
	// Util::Print(col, outputname, fOutputDir + fOutputSubDir, fOutputFile);
}

//____________________________________________________________________________
void AnaClass::plotEffOverlayEG(TEfficiency *h1in, TString tag1, TGraphAsymmErrors *h2in, TString tag2, bool logy){
	gStyle->SetOptStat("");
	TEfficiency       *h1 = new TEfficiency(*h1in);
	TGraphAsymmErrors *h2 = new TGraphAsymmErrors(*h2in);

	char canvtitle[100], canvname[100];
	sprintf(canvtitle,"%s vs %s", h1->GetName(), h2->GetName());
	sprintf(canvname,"%s:%s", h1->GetName(), h2->GetName());
	TCanvas *col = new TCanvas(canvname, canvtitle, 0, 0, 900, 700);
	col->SetFillStyle(0);
	col->SetFrameFillStyle(0);
	col->cd();
	gPad->SetFillStyle(0);
	gPad->SetGridy(1);
	if(logy) col->SetLogy(1);

	// TLegend *leg = new TLegend(0.60,0.15,0.886,0.28); // Lower right
	TLegend *leg = new TLegend(0.60,0.75,0.886,0.88); // Upper right
	leg->AddEntry(h1, tag1,"P");
	leg->AddEntry(h2, tag2,"f");
	leg->SetFillStyle(0);
	leg->SetTextFont(42);
	leg->SetBorderSize(0);

	h1->Draw("AP");
	h2->Draw("Psame");
	leg->Draw();

	gPad->RedrawAxis();
	TString outputname = TString(h1->GetName()) + "_" + TString(h2->GetName());
	Util::PrintNoEPS(col, outputname, fOutputDir + fOutputSubDir, fOutputFile);
	// Util::Print(col, outputname, fOutputDir + fOutputSubDir, fOutputFile);
}
void AnaClass::plotEffOverlayEE(TEfficiency *h1in, TString tag1, TEfficiency *h2in, TString tag2, bool logy){
	gStyle->SetOptStat("");
	TEfficiency *h1 = new TEfficiency(*h1in);
	TEfficiency *h2 = new TEfficiency(*h2in);

	char canvtitle[100], canvname[100];
	sprintf(canvtitle,"%s vs %s", h1->GetName(), h2->GetName());
	sprintf(canvname,"%s:%s", h1->GetName(), h2->GetName());
	TCanvas *col = new TCanvas(canvname, canvtitle, 0, 0, 900, 700);
	col->SetFillStyle(0);
	col->SetFrameFillStyle(0);
	col->cd();
	gPad->SetFillStyle(0);
	gPad->SetGridy(1);
	if(logy) col->SetLogy(1);

	// TLegend *leg = new TLegend(0.60,0.15,0.886,0.28); // Lower right
	TLegend *leg = new TLegend(0.60,0.75,0.886,0.88); // Upper right
	leg->AddEntry(h1, tag1,"P");
	leg->AddEntry(h2, tag2,"f");
	leg->SetFillStyle(0);
	leg->SetTextFont(42);
	leg->SetBorderSize(0);

	h1->Draw("AP");
	h2->Draw("Psame");
	leg->Draw();

	gPad->RedrawAxis();
	TString outputname = TString(h1->GetName()) + "_" + TString(h2->GetName());
	Util::PrintNoEPS(col, outputname, fOutputDir + fOutputSubDir, fOutputFile);
	// Util::Print(col, outputname, fOutputDir + fOutputSubDir, fOutputFile);
}


/*****************************************************************************
###################| Utilities |##############################################
*****************************************************************************/
//____________________________________________________________________________
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

//____________________________________________________________________________
TH1D* AnaClass::drawTree1D(const char* arg, const TCut reqs, const char* histn, const int nbins, const double xmin, const double xmax, TTree* tree, bool draw, const char* drawopt){
	int nbins_auto = nbins;
	if(nbins == 0)	nbins_auto = OptNBins(tree->Draw(arg, reqs, "goff"));
	TH1D* h1 = new TH1D(histn,histn,nbins_auto,xmin,xmax);
	h1->SetXTitle(convertVarName(arg));
	tree->Project(histn,arg,reqs);
	if(draw) h1->Draw(drawopt);
	return h1;
}
TH1D* AnaClass::drawTree1D(const char* arg, const TCut reqs, const char* histn, const int nbins, const double *xbins, TTree* tree, bool draw, const char* drawopt){
	int nbins_auto = nbins;
	if(nbins == 0)	nbins_auto = OptNBins(tree->Draw(arg, reqs, "goff"));
	TH1D* h1 = new TH1D(histn,histn,nbins_auto,xbins);
	h1->SetXTitle(convertVarName(arg));
	tree->Project(histn,arg,reqs);
	if(draw) h1->Draw(drawopt);
	return h1;
}
TH2D* AnaClass::drawTree2D(const char* arg1, const char* arg2, const TCut reqs, const char* histn, const int nbinsx, const double xmin, const double xmax, const int nbinsy, const double ymin, const double ymax, TTree* tree, bool draw, const char* drawopt){
	char out[1000];
	int nbinsx_auto(nbinsx), nbinsy_auto(nbinsy);
	if(nbinsx == 0) nbinsx_auto = OptNBins(tree->Draw(arg1, reqs, "goff"));
	if(nbinsy == 0) nbinsy_auto = OptNBins(tree->Draw(arg2, reqs, "goff"));
	TH2D* h1 = new TH2D(histn,histn,nbinsx_auto,xmin,xmax,nbinsy_auto,ymin,ymax);
	h1->SetXTitle(convertVarName(arg1));
	h1->SetYTitle(convertVarName(arg2));
	sprintf(out,"%s:%s",arg2,arg1);
	tree->Project(histn,out,reqs);
	if(draw) h1->Draw(drawopt);
	return h1;
}

//____________________________________________________________________________
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

//____________________________________________________________________________
int AnaClass::OptNBins(int nentries){
// Returns an 'optimal' number of bins for some number of entries
	if(nentries < 100)  return 10;
	if(nentries < 2000) return 60;
	if(nentries < 5000) return 80;
	if(nentries < 10000) return 100;
	else return 200;
}

//____________________________________________________________________________
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
TH1D* AnaClass::normHistBW(const TH1D *ihist, float scale){
/*		-	Normalizes a histogram (incl. errors) with variable binwidth.
			I.e. divides each entry by the binwidth									*/
	TH1D *ohist = new TH1D(*ihist);
	for( int i = 0; i < ihist->GetNbinsX()+2; i++ ){
		ohist->SetBinContent(i, scale * ihist->GetBinContent(i)/ihist->GetBinWidth(i));
		ohist->SetBinError(i,   scale * ihist->GetBinError(i)/ihist->GetBinWidth(i));
	}
	return ohist;
}

//____________________________________________________________________________
const Double_t* AnaClass::getBinning(const TH1D *hist){
	return hist->GetXaxis()->GetXbins()->GetArray();
}

//____________________________________________________________________________
void AnaClass::setPlottingRange(TH1D *&h1, float margin, bool logy){
	// Determine plotting range
	// Default margin is 0.05	
	float max = getMaxYExtension(h1);
	float min = getMinYExtension(h1);
	
	float range = max-min;

	if(logy){
		max *= 1.5;
		if(min <= 0.){
			float minval = h1->GetBinContent(h1->GetMinimumBin());
			min = minval/1.5;
		}
		else min /= 1.5;
	}
	else{
		max = max + margin * range;
		min = min - margin * range;
	}

	h1->SetMinimum(min);
	h1->SetMaximum(max);
}
void AnaClass::setPlottingRange(TH1D *&h1, TH1D *&h2, float margin, bool logy){
	// Determine plotting range
	// Default margin is 0.05	
	float max1 = getMaxYExtension(h1);
	float max2 = getMaxYExtension(h2);
	
	float max  = (max1>max2)?max1:max2;

	float min1 = getMinYExtension(h1);
	float min2 = getMinYExtension(h2);
	float min  = (min1<min2)?min1:min2;
	
	float range = max-min;

	if(logy){
		max *= 1.5;
		if(min <= 0.){
			float minval1 = h1->GetBinContent(h1->GetMinimumBin());
			float minval2 = h2->GetBinContent(h2->GetMinimumBin());
			float minval = minval1<minval2?minval1:minval2;
			min = minval/1.5;
		}
		else min /= 1.5;
	}
	else{
		max = max + margin * range;
		min = min - margin * range;
	}

	h1->SetMinimum(min);
	h1->SetMaximum(max);
	h2->SetMinimum(min);
	h2->SetMaximum(max);
}
void AnaClass::setPlottingRange(TH1D *&h1, TH1D *&h2, TH1D *&h3, float margin, bool logy){
	// Determine plotting range
	// Default margin is 0.05
	float max1 = getMaxYExtension(h1);
	float max2 = getMaxYExtension(h2);
	float max3 = getMaxYExtension(h3);
	float max12 = (max1>max2)?max1:max2;
	float max   = (max12>max3)?max12:max3;

	float min1 = getMinYExtension(h1);
	float min2 = getMinYExtension(h2);
	float min3 = getMinYExtension(h3);
	float min12 = (min1<min2)?min1:min2;
	float min   = (min12<min3)?min12:min3;
	
	float range = max-min;

	if(logy){
		max *= 1.5;
		if(min <= 0.){
			float minval1 = h1->GetBinContent(h1->GetMinimumBin());
			float minval2 = h2->GetBinContent(h2->GetMinimumBin());
			float minval3 = h3->GetBinContent(h3->GetMinimumBin());
			float minval12 = minval1<minval2?minval1:minval2;
			float minval = minval12<minval3?minval12:minval3;
			min = minval/1.5;
		}
		else min /= 1.5;
	}
	else{
		max = max + margin * range;
		min = min - margin * range;
	}

	h1->SetMinimum(min);
	h1->SetMaximum(max);
	h2->SetMinimum(min);
	h2->SetMaximum(max);
	h3->SetMinimum(min);
	h3->SetMaximum(max);
}
void AnaClass::setPlottingRange(std::vector<TH1D*> &hists, float margin, bool logy){
	// Determine plotting range
	// Default margin is 0.05
	float max = hists[0]->GetMaximum();
	float min = hists[0]->GetMinimum();
	for(size_t i = 0; i < hists.size(); ++i){
		float tempmax = getMaxYExtension(hists[i]);
		float tempmin = getMinYExtension(hists[i]);
		if(tempmax > max) max = tempmax;
		if(tempmin < min) min = tempmin;
	}
	float range = max-min;

	if(logy){
		max *= 1.5;
		if(min <= 0.){
			float minval = hists[0]->GetBinContent(hists[0]->GetMinimumBin());
			for(size_t i = 0; i < hists.size(); ++i){
				float tempval = hists[i]->GetBinContent(hists[i]->GetMinimumBin());
				if(tempval < minval) minval = tempval;
			}
			min = minval/1.5;
		}
		else min /= 1.5;
	}
	else{
		max = max + margin * range;
		min = min - margin * range;
	}

	for(size_t i = 0; i < hists.size(); ++i){
		hists[i]->SetMinimum(min);
		hists[i]->SetMaximum(max);
	}
}
void AnaClass::getPlottingRange(float &minopt, float &maxopt, std::vector<TH1D*> hists, float margin, bool logy){
	// Determine plotting range
	// Default margin is 0.05
	float max = hists[0]->GetBinContent(1);
	float min = hists[0]->GetBinContent(1);
	for(size_t i = 0; i < hists.size(); ++i){
		float tempmax = getMaxYExtension(hists[i]);
		float tempmin = getMinYExtension(hists[i]);
		if(tempmax > max) max = tempmax;
		if(tempmin < min) min = tempmin;
	}
	float range = max-min;

	if(logy){
		max *= 1.5;
		if(min <= 0.){
			float minval = hists[0]->GetBinContent(hists[0]->GetMinimumBin());
			for(size_t i = 0; i < hists.size(); ++i){
				float tempval = hists[i]->GetBinContent(hists[i]->GetMinimumBin());
				if(tempval < minval) minval = tempval;
			}
			min = minval/1.5;
		}
		else min /= 1.5;
	}
	else{
		max = max + margin * range;
		min = min - margin * range;
	}

	minopt = min;
	maxopt = max;
	return;
}

//____________________________________________________________________________
float AnaClass::getMaxYExtension(TH1 *h){
	float max = h->GetBinContent(1);
	// float max = h->GetMaximum();
	for(size_t i = 1; i <= h->GetNbinsX(); ++i){
		float temp = h->GetBinContent(i) + h->GetBinError(i);
		if(temp > max) max = temp;
	}
	return max;
}
float AnaClass::getMinYExtension(TH1 *h){
	float min = h->GetBinContent(1);
	// float min = h->GetMinimum();
	for(size_t i = 1; i <= h->GetNbinsX(); ++i){
		float temp = h->GetBinContent(i) - h->GetBinError(i);
		if(temp < min) min = temp;
	}
	return min;
}

//____________________________________________________________________________
void AnaClass::setZeroBinError(TH1D *ihist){
//    - Sets all bin errors to zero
	for( int i = 0; i < ihist->GetNbinsX()+2; i++ ) ihist->SetBinError(i,0.);
}

//____________________________________________________________________________
void AnaClass::fillWithoutOF(TH1D *&ihist, double x, double w){
	double xmax = ihist->GetBinLowEdge(ihist->GetNbinsX());
	// double xmax = ihist->GetBinLowEdge(ihist->GetMaximumBin());
	double bw = ihist->GetBinWidth(ihist->GetMaximumBin());
	if(x > xmax) ihist->Fill(xmax + bw*0.5, w); // always increment last bin (i.e. never the overflow)
	else ihist->Fill(x, w);
}

//____________________________________________________________________________
TCanvas* AnaClass::makeCanvas(const char* name){
/*		-	Creates a TCanvas at a somewhat random position							 */
	int x = (int)gRandom->Uniform(10,200);
	int y = (int)gRandom->Uniform(0,70);
	TCanvas* c1 = new TCanvas(name, name, x, y, 900, 700);
	return c1;
}

//____________________________________________________________________________
void AnaClass::printObject(TObject* o, TString name, Option_t *drawopt, bool logy){
	TCanvas *col = new TCanvas(o->GetName(), o->GetTitle(), 0, 0, 900, 700);
	col->SetFillStyle(0);
	col->SetFrameFillStyle(0);
	col->cd();
	gPad->SetFillStyle(0);
	o->Draw(drawopt);
	gPad->RedrawAxis();
	if(logy) gPad->SetLogy(1);
	Util::PrintNoEPS(col, name, fOutputDir + fOutputSubDir, fOutputFile);
}

//____________________________________________________________________________
TH1D* AnaClass::bookTH1D(const char* name, const char* title, int nbins, double xlow, double xup){
	TH1D *h = new TH1D(name, title, nbins, xlow, xup);
	h->SetXTitle(convertVarName(name));
	return h;
}

//____________________________________________________________________________
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

//____________________________________________________________________________
void AnaClass::printProgress(int entry, const int nentries, TString header, int step){
	if(step == -1){
		step = nentries/20;
		if( step < 200 )   step = 200;
		if( step > 10000 ) step = 10000;
	}
	if(entry%step != 0 && (entry+1 != nentries) ) return;
	
	float progress_f = (float)(entry+1)/(float)(nentries)*100.;
	char progress[10];
	sprintf(progress, "%5.1f", progress_f);
	cout << " Processing " << setw(50) << left << header << setw(6) << right << progress << " %      \r" << flush;
	if(entry+1 == nentries) cout << endl;
}

//____________________________________________________________________________
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

//____________________________________________________________________________
void AnaClass::refValues(const char* var, TH1D* h){
	double percent1 = 0.05;
	double percent2 = 0.01;
	if (!strcmp(var, "MuPt")         || !strcmp(var, "ElPt")      || !strcmp(var, "JPt") ||
		!strcmp(var, "MuJESCorrMET") || !strcmp(var, "TCMET")     ||
		!strcmp(var, "PFMET")        || !strcmp(var, "SumEt")     ||
		!strcmp(var, "ECALSumEt")    || !strcmp(var, "HCALSumEt") ||
	!strcmp(var, "PrimVtxPtSum") || !strcmp(var, "TrkPtSum") ) {
		tailFraction(h, percent1);
		tailFraction(h, percent2);
	}
}

//____________________________________________________________________________
void AnaClass::getWeightedYMeanRMS(TH1D *h, double &mean, double &rms){
	vector<double> values;
	vector<double> weights;
	for(size_t i = 1; i <= h->GetNbinsX(); ++i){
		if(h->GetBinError(i) == 0.) continue;
		values.push_back(h->GetBinContent(i));
		weights.push_back(1./h->GetBinError(i));
	}
	mean = TMath::Mean(values.begin(), values.end(), weights.begin());
	rms = TMath::RMS(values.begin(), values.end());
	return;
}

//____________________________________________________________________________
double AnaClass::tailFraction(TH1D* h, double frac){
	double binValue = -999.;
	int nbins = h->GetNbinsX();
	double tail = frac * h->Integral();
	double tailSum = 0.;

	int loc = 0;
	for (int i = nbins+1; i >= 0; --i) {
		tailSum += h->GetBinContent(i);
		if (tailSum < tail) loc = i;
		if (tailSum >= tail) break;
	}
	if (loc > 0 && loc < nbins) {
		binValue = h->GetBinCenter(loc);
		double maxY = 0.5 * h->GetMaximum();
		double minY = 0.05 * h->GetMaximum();
		TArrow* arr = new TArrow(binValue, maxY, binValue, minY, 0.02, "|>");
		arr->SetLineColor(4);
		arr->SetFillColor(4);
		arr->Draw();
		TLatex l;
		l.SetTextColor(4);
		l.SetTextSize(0.04);
		char tit[50];
		l.SetTextAlign(23);
		sprintf(tit, "%6.2f ", (float)frac);
		l.DrawLatex(binValue, 1.1*maxY, tit);
		l.SetTextAlign(12);
		sprintf(tit, "%8.3f", (float)binValue);
		l.DrawLatex(binValue, 0.8*maxY, tit);
	}

	return binValue;
}
void AnaClass::printCheckList(const char* var, TH1D* h, const char* filename){
// Prints the CheckList to file filename
	ofstream file;
	file.open(filename, ios::app);
	double percent1 = 0.05;
	double percent2 = 0.01;
	double etaLow = 1.44, etaHigh = 3.;
	double d0PVLow = 0.02, d0PVHigh = 0.5;
	double isoLow = 0.3, isoHigh = 1.;

	if( !strcmp(var, "MuPt")         || !strcmp(var, "ElPt")         || !strcmp(var, "JPt") ||
		!strcmp(var, "TCMET")        || !strcmp(var, "MuJESCorrMET") || !strcmp(var, "PFMET") ||
		!strcmp(var, "SumEt")        || !strcmp(var, "ECALSumEt")    || !strcmp(var, "HCALSumEt") ||
	!strcmp(var, "PrimVtxPtSum") || !strcmp(var, "TrkPtSum") ){
		file << "* " << var << endl;
		file << printAverage(var, h) << endl;
		file << printTailFraction(var, h, percent1) << endl;
		file << printTailFraction(var, h, percent2) << endl;
	}
	if( !strcmp(var, "MuDzPV")       || !strcmp(var, "MuNChi2")      ||
		!strcmp(var, "MuNTkHits")      || !strcmp(var, "ElDzPV")       ||
		!strcmp(var, "Elfbrem")        || !strcmp(var, "JEMfrac")      || !strcmp(var, "JCHfrac")      ||
		!strcmp(var, "JNConstituents") || !strcmp(var, "MuIso03SumPt") ||
		!strcmp(var, "MuIso03EmEt")    || !strcmp(var, "MuIso03HadEt") || 
	!strcmp(var, "ElPtSum")        || !strcmp(var, "ElEmEtSum")    || !strcmp(var, "ElHadEtSum") ){
		file << "* " << var << endl;
		file << printAverage(var, h) << endl;
	}
	if( !strcmp(var, "MuEta") || !strcmp(var, "ElEta") || !strcmp(var, "JEta") ){
		file << "* " << var << endl;
		file << printAverage(var, h) << endl;
		file << printRatio(var, h, etaLow, etaHigh, -etaHigh, etaHigh) << endl;
		file << printRatio(var, h, -etaHigh, -etaLow, -etaHigh, etaHigh) << endl;
		file << printRatio(var, h, 0., etaLow, -etaLow, 0.) << endl;
	}
	if( !strcmp(var, "PrimVtxx")       || !strcmp(var, "PrimVtxy")     || !strcmp(var, "PrimVtxz") ||
		!strcmp(var, "PrimVtxNdof") || !strcmp(var, "PrimVtxNChi2") ||!strcmp(var, "NTracks")   ||
	!strcmp(var, "MuEem")          || !strcmp(var, "MuEHad") ){
		file << "* " << var << endl;
		file << printAverage(var, h) << endl;
	}
	if( !strcmp(var, "MuD0PV") || !strcmp(var, "ElD0PV") ){
		file << "* " << var << endl;
		file << printAverage(var, h) << endl;
		file << printRatio(var, h, -d0PVLow, d0PVLow, -d0PVHigh, d0PVHigh) << endl;
	}
	if( !strcmp(var, "MuRelIso03") || !strcmp(var, "ElRelIso04") ){
		file << "* " << var << endl;
		file << printAverage(var, h) << endl;;
		file << printRatio(var, h, 0., isoLow, 0., isoHigh) << endl;
	}
}
TString AnaClass::printTailFraction(const char* var, TH1D* h, double frac){
// Prints the value of the variable for which frac remains in the tail

	double binValue = -999.;
	double dbinValue = -999.;
	int nbins = h->GetNbinsX();
	double nevts = h->Integral();
	double tail = frac * nevts;
	double dtail = frac * sqrt(nevts);
	double tailp = tail + dtail;
	double tailm = tail - dtail;
	double tailSum = 0.;

	int loc = 0, locp = 0, locm = 0;
	for (int i = nbins+1; i >= 0; --i) {
		tailSum += h->GetBinContent(i);
		if (tailSum < tail) loc = i;
		if (tailSum < tailp) locp = i;
		if (tailSum < tailm) locm = i;
		if (tailSum >= tailp) break;
	}
	if (loc > 0 && loc < nbins){
		binValue = h->GetBinCenter(loc);
		if (locm > 0 && locp < nbins) {
			dbinValue = 0.5*(h->GetBinCenter(locm)-h->GetBinCenter(locp));
		} else if (locp <= 0 && locm < nbins) {
			dbinValue = binValue - h->GetBinCenter(locp);
		} else if (locp > 0 && locm >= nbins) {
			dbinValue = h->GetBinCenter(locp) - binValue;
		}
	}
	if (nevts > 0) {
		double binSize = h->GetBinWidth(loc);
		if (binSize > dbinValue) dbinValue = binSize;
	} else {
		binValue = 0.;
		dbinValue = 0.;
	}
	TString result = Form("  For %f of tail %s = %f +- %f", frac, var, binValue, dbinValue);
	return result;
}
TString AnaClass::printAverage(const char* var, TH1D* h) {
// Prints the average of the histogram

	double nevts = h->GetEntries();
	double aver, daver;
	if (nevts > 0) {
		aver = h->GetMean(1);
		double rms = h->GetRMS(1);
		daver = rms / sqrt(nevts);
	} else {
		aver = 0.;
		daver = 0.;
	}
	TString result = Form("  Mean value of %s = %f +- %f", var, aver, daver);
	return result;
}
TString AnaClass::printRatio(const char* var, TH1D* h, double x1, double x2, double y1, double y2){
// Prints the ratio of entries for which (x1<var<x2) / (y1<var<y2)

	int nbins = h->GetNbinsX();
	double sumx1y1 = 0.;
	double sumx2y2 = 0.;
	double sumx1y2 = 0.;
	double sumx2y1 = 0.;
	double xmin = h->GetXaxis()->GetXmin();
	double xmax = h->GetXaxis()->GetXmax();
	if (x1 < xmin) x1 = xmin;
	if (x2 > xmax) x2 = xmax;
	if (y1 < xmin) y1 = xmin;
	if (y2 > xmax) y2 = xmax;

	for (int i = 1; i < nbins; ++i) {
		double xy = h->GetBinCenter(i);
		double content = h->GetBinContent(i);
		if ( (xy - x1)*(xy - y1) < 0.) sumx1y1 += content;
		if ( (xy - x2)*(xy - y2) < 0.) sumx2y2 += content;
		if ( (xy - x1)*(xy - y2) < 0.) sumx1y2 += content;
		if ( (xy - x2)*(xy - y1) < 0.) sumx2y1 += content;
	}
	double xy1 = (x1 - y1)*(x1 - y2);
	double xy2 = (x2 - y1)*(x2 - y2);

	double xunc = 0.;
	double yunc = 0.;
	double xycor = 0.;
	if (xy1 <= 0. && xy2 <= 0.) {
		xunc  = 0.;
		yunc  = sumx1y1 + sumx2y2;
		xycor = sumx1y2 - sumx2y2;
	} else if (xy1 > 0. && xy2 <= 0.) {
		xunc  = sumx1y1;
		yunc  = sumx2y2;
		xycor = sumx2y1;
	} else if (xy1 > 0. && xy2 > 0.) {
		xunc  = sumx1y1 - sumx2y1;
		yunc  = sumx2y2 - sumx2y1;
		xycor = 0.;
	} else if (xy1 <= 0. && xy2 > 0.) {
		xunc  = sumx2y2;
		yunc  = sumx1y1;
		xycor = sumx1y2;
	}

	double rat = -999., drat = -999.;
	if (yunc+xycor != 0.) {
		rat = (xunc+xycor) / (yunc+xycor);
		drat = sqrt( (yunc+xycor)*(yunc+xycor)*xunc
			+ (xunc+xycor)*(xunc+xycor)*yunc
			+ (yunc-xunc)*(yunc-xunc)*xycor )
			/ ( (yunc+xycor)*(yunc+xycor) );
	}
	TString result = Form("  Ratio (%f<%s<%f) / (%f<%s<%f) = %f +- %f", x1, var, x2, y1, var, y2, rat, drat);
	return result;
}

