#include "MultiplicityAnalysis.hh"
#include "base/TreeReader.hh"
#include "helper/LeptJetStat.h"
#include "helper/Utilities.hh"

using namespace std;

MultiplicityAnalysis::MultiplicityAnalysis(TreeReader *tr) : UserAnalysisBase(tr){
}

MultiplicityAnalysis::~MultiplicityAnalysis(){
}

void MultiplicityAnalysis::Begin(const char* filename){
	fMPHistFile = new TFile(fOutputDir + TString(filename), "RECREATE");
	// Temp objects:
	TString hname, htit;

	fMyLeptJetStat = new LeptJetStat();
	fHljMult  = new TH2D("ljMult", "Lepton / Jets multiplicity", 13, 0, 13, 7, 0, 7);
	fHemuMult = new TH2D("emuMult", "e/mu multiplicity", 18, 0, 18, 7, 0, 7);
	fHemuEff  = new TH1F("emuEffic", "e/mu Efficiency", 13, 0, 13);

	fHljMult->SetStats(false);
	fHemuMult->SetStats(false);
	fHemuEff->SetStats(false);
}

void MultiplicityAnalysis::Analyze(){
	//define an array for leptons (es+mus) with
	//information on pt, eta and category
	//category 1: e+
	//category 2: e-
	//category 3: mu+
	//category 4: mu-
	const int NLepts = 40;
	int LeptCat[NLepts];
	double LeptPt[NLepts];
	double LeptEta[NLepts];

	for(int tmp=0; tmp<20; ++tmp){
		LeptCat[tmp]=0;
		LeptPt[tmp]=-999.;
		LeptEta[tmp]=-999.;
	}
	//loop over es and mus and fill the leptons array
	int ILept=-1;
	for(int i=0; i< fTR->NEles; ++i){
		//////////////////////////////////
		// Your electron selection here //
		//////////////////////////////////
		ILept++;
		LeptPt[ILept]  = fTR->ElPt[i];
		LeptEta[ILept] = fTR->ElEta[i];
		if(fTR->ElCharge[i]>0) {
			LeptCat[ILept]=1;
		} else {
			LeptCat[ILept]=2;
		}
	}
	for(int i=0; i< fTR->NMus; ++i){
		//////////////////////////////////
		// Your muon selection here     //
		//////////////////////////////////
		ILept++;
		LeptPt[ILept]  = fTR->MuPt[i];
		LeptEta[ILept] = fTR->MuEta[i];
		if(fTR->MuCharge[i]>0) {
			LeptCat[ILept]=3;
		} else {
			LeptCat[ILept]=4;
		}
	}
	//total number of leptons in the event
	int NLeptons= ILept+1;

// reorder here ...
	bool changed = false;
	do {
		changed = false;
		for (int i = 0; i < NLeptons; i++){
			for (int j = i+1; j < NLeptons; j++){
				if (LeptPt[i] < LeptPt[j] ){
					int cat=LeptCat[i];
					double pt=LeptPt[i];
					double eta=LeptEta[i];
					LeptCat[i]=LeptCat[j];
					LeptPt[i]=LeptPt[j];
					LeptEta[i]=LeptEta[j];
					LeptCat[j]=cat;
					LeptPt[j]=pt;
					LeptEta[j]=eta;
					changed = true;
				}
			}
		}
	} while (changed);

// if too many leptons, remove the softest ones
	if (NLeptons > 4){
		LeptCat[0]=5;
	}

// count number of good jets (beni)
	unsigned int NQJets = 0;
	for(int i=0; i < fTR->NJets; ++i){
		//////////////////////////////////
		// Your jet selection here      //
		//////////////////////////////////		
		NQJets++;
	}
// convert the lepton config into the index and count
	fMyLeptJetStat->FillLeptJetStat(LeptCat, NQJets, 0);
}

void MultiplicityAnalysis::End(){
	fMyLeptJetStat->FillShortTables();
	if(fVerbose) fMyLeptJetStat->PrintConfigs();
	PlotMPSummary();
	PlotMPEffic();

	fTlat->SetTextColor(kBlack);
	fTlat->SetNDC(kTRUE);
	fTlat->SetTextSize(0.04);

	// const int ncol = 11;
	// int colors[ncol] = { 10, 16, 5, 28, 29, 8, 4, 9, 45, 46, 2};

	TString subdir = "MultiplicityPlots";
	TCanvas *canv;
	TString canvtitle = "Lepton / Jets multiplicity";
	canv = new TCanvas("ljMult", canvtitle , 0, 0, 900, 700);
	canv->SetRightMargin(0.15);
	// gStyle->SetPalette(ncol, colors);
	gPad->SetTheta(50);
	gPad->SetPhi(240);
	fHljMult->SetMinimum(0);
	fHljMult->DrawCopy("colz");
	// fHljMult->DrawCopy("lego2 Z");
	fTlat->DrawLatex(0.11,0.92, canvtitle);
	Util::PrintPNG(canv, fTag + "_ljMult", fOutputDir + subdir);
	Util::PrintEPS(canv, fTag + "_ljMult", fOutputDir + subdir);

	canvtitle = "e/mu multiplicity";
	canv = new TCanvas("emuMult", canvtitle , 0, 0, 900, 700);
	canv->SetRightMargin(0.15);
	gPad->SetTheta(50);
	gPad->SetPhi(240);
	fHemuMult->SetMinimum(0);
	fHemuMult->DrawCopy("colz");
	// fHemuMult->DrawCopy("lego2 Z");
	fTlat->DrawLatex(0.11,0.92, canvtitle);
	Util::PrintPNG(canv, fTag + "_emuMult", fOutputDir + subdir);
	Util::PrintEPS(canv, fTag + "_emuMult", fOutputDir + subdir);

	canvtitle = "e/mu Efficiency";
	canv = new TCanvas("emuEffic", canvtitle , 0, 0, 900, 700);
	canv->SetRightMargin(0.15);
	gPad->SetTheta(50);
	gPad->SetPhi(240);
	fHemuEff->DrawCopy();
	fTlat->DrawLatex(0.11,0.92, canvtitle);
	Util::PrintPNG(canv, fTag + "_emuEffic", fOutputDir+subdir);
	Util::PrintEPS(canv, fTag + "_emuEffic", fOutputDir+subdir);


	fMPHistFile->cd();
	fHljMult->Write();
	fHemuMult->Write();
	fHemuEff->Write();
	fMPHistFile->Close();
}

void MultiplicityAnalysis::PlotMPSummary(){
// Makes a 2D plot of lepton configurations versus jet multiplicity
//  in the lepton configurations, e and mu are summed

//Collect the entries for the plot
	int index;
	int nperJets[9];
	const int nlept = 13;
	const int njets = 7;
	int multable[nlept][njets];
	const char * lablx[nlept] = {"0l    ", "1l+   ", "1l-   ", "OSll  ", "OSem  ",
		"SSll  ", "SSem  ", "OS2ll+", "OS2ll-", "OSeml+",
		"OSeml-", "SS3l  ", "4lincl"};
	const char * lably[njets] = {" 0j"," 1j", " 2j", " 3j", " 4j", " 5j", ">5j"};
	int indTab2l[10] = { 5,  3,  6,  4,  5,  4,  6,  5,  3,  5};
	int indTab3l[20] = {11,  7, 11,  9,  8,  7,  8, 11,  7, 10,
		11, 10, 11,  9,  8, 11, 11,  7,  8, 11};

// initialize the multiplicity table
	for (int i = 0; i < nlept; ++i) {
		for (int j = 0; j < njets; ++j) {
			multable[i][j] = 0;
		}
	}
	int imax = njets + 1;
// nber of leptons = 0
	index = fMyLeptJetStat->GetConfigfrOrder(0);
	fMyLeptJetStat->NEntriesPerJetMult(index, nperJets);
	for (int i = 1; i < imax; ++i) {multable[0][i-1] = nperJets[i];}

// nber of leptons = 1
	index = fMyLeptJetStat->GetConfigfrOrder(1);
	fMyLeptJetStat->NEntriesPerJetMult(index, nperJets);
	for (int i = 1; i < imax; ++i) {multable[1][i-1] = nperJets[i];}

	index = fMyLeptJetStat->GetConfigfrOrder(3);
	fMyLeptJetStat->NEntriesPerJetMult(index, nperJets);
	for (int i = 1; i < imax; ++i) {multable[1][i-1] += nperJets[i];}

	index = fMyLeptJetStat->GetConfigfrOrder(2);
	fMyLeptJetStat->NEntriesPerJetMult(index, nperJets);
	for (int i = 1; i < imax; ++i) {multable[2][i-1] = nperJets[i];}

	index = fMyLeptJetStat->GetConfigfrOrder(4);
	fMyLeptJetStat->NEntriesPerJetMult(index, nperJets);
	for (int i = 1; i < imax; ++i) {multable[2][i-1] += nperJets[i];}

// nber of leptons = 2
	for (int i = 0; i < 10; ++i) {
		int ii = i + 5;
		index = fMyLeptJetStat->GetConfigfrOrder(ii);
		fMyLeptJetStat->NEntriesPerJetMult(index, nperJets);
//    cout << " ii = << ii << "
		for (int j = 1; j < 8; ++j) {multable[indTab2l[i]][j-1] += nperJets[j];}
	}
//  cout << " 2l done " << endl;

// nber of leptons = 3
	for (int i = 0; i < 20; ++i) {
		int ii = i + 15;
		index = fMyLeptJetStat->GetConfigfrOrder(ii);
		fMyLeptJetStat->NEntriesPerJetMult(index, nperJets);
		for (int j = 1; j < imax; ++j) {multable[indTab3l[i]][j-1] += nperJets[j];}
	}
//  cout << " 3l done " << endl;

// nber of leptons >= 4
	for (int m = 35; m < 71; ++m){
		index = fMyLeptJetStat->GetConfigfrOrder(m);
		fMyLeptJetStat->NEntriesPerJetMult(index, nperJets);
		for (int i = 1; i < imax; ++i) {multable[12][i-1] += nperJets[i];}
	}

// now fill the 2D plot
	if(fVerbose) cout << endl;
	if(fVerbose) 	cout << " Contents of the lepton/jets multiplicity 2D plot " << endl;
	for (int i = 0; i < nlept; ++i) {
		if(fVerbose) cout << "  " << lablx[i];
		for (int j = 0; j < njets; ++j) {
			fHljMult->Fill(lablx[i], lably[j], multable[i][j]);
			if(fVerbose) cout << "  " << multable[i][j];
		}
		if(fVerbose) cout << endl;
	}
	return;
}

void MultiplicityAnalysis::PlotMPEffic(){
// Makes a 2D plot of lepton configurations versus jet multiplicity
//  in the lepton configurations, + and - charges are summed
// Makes a profile histogram of e/mu efficiency ratios

//Collect the entries for the plot
	int index;
	int nperJets[9];
	const int nlept = 18;
	const int njets = 7;
	int multable[nlept][njets];
	const char* lablx[nlept] = {"1m    ", "1e    ", "OSmm  ", "OSem  ", "OSee  ",
		"SSmm  ", "SSem  ", "SSee  ", "OSmm1m", "OSmm1e",
		"OSee1m", "OSee1e", "SSmm1e", "SSee1m", "3SSmmm",
		"3SSmme", "3SSeem", "3SSeee"};
	const char* lably[njets] = {" 0j"," 1j", " 2j", " 3j", " 4j", " 5j", ">5j"};

// initialize the multiplicity table
	for (int i = 0; i < nlept; ++i) {
		for (int j = 0; j < njets; ++j) {
			multable[i][j] = 0;
		}
	}
	int imax = njets + 1;
// nber of leptons = 1m
	index = fMyLeptJetStat->GetConfigfrOrder(3);
	fMyLeptJetStat->NEntriesPerJetMult(index, nperJets);
	for (int i = 1; i < imax; ++i) {multable[0][i-1] = nperJets[i];}
	index = fMyLeptJetStat->GetConfigfrOrder(4);
	fMyLeptJetStat->NEntriesPerJetMult(index, nperJets);
	for (int i = 1; i < imax; ++i) {multable[0][i-1] += nperJets[i];}

// nber of leptons = 1e
	index = fMyLeptJetStat->GetConfigfrOrder(1);
	fMyLeptJetStat->NEntriesPerJetMult(index, nperJets);
	for (int i = 1; i < imax; ++i) {multable[1][i-1] = nperJets[i];}
	index = fMyLeptJetStat->GetConfigfrOrder(2);
	fMyLeptJetStat->NEntriesPerJetMult(index, nperJets);
	for (int i = 1; i < imax; ++i) {multable[1][i-1] += nperJets[i];}

// nber of leptons = OSmm
	index = fMyLeptJetStat->GetConfigfrOrder(13);
	fMyLeptJetStat->NEntriesPerJetMult(index, nperJets);
	for (int i = 1; i < imax; ++i) {multable[2][i-1] = nperJets[i];}

// nber of leptons = OSem
	index = fMyLeptJetStat->GetConfigfrOrder(8);
	fMyLeptJetStat->NEntriesPerJetMult(index, nperJets);
	for (int i = 1; i < imax; ++i) {multable[3][i-1] = nperJets[i];}
	index = fMyLeptJetStat->GetConfigfrOrder(10);
	fMyLeptJetStat->NEntriesPerJetMult(index, nperJets);
	for (int i = 1; i < imax; ++i) {multable[3][i-1] += nperJets[i];}

// nber of leptons = OSee
	index = fMyLeptJetStat->GetConfigfrOrder(6);
	fMyLeptJetStat->NEntriesPerJetMult(index, nperJets);
	for (int i = 1; i < imax; ++i) {multable[4][i-1] = nperJets[i];}

// nber of leptons = SSmm
	index = fMyLeptJetStat->GetConfigfrOrder(12);
	fMyLeptJetStat->NEntriesPerJetMult(index, nperJets);
	for (int i = 1; i < imax; ++i) {multable[5][i-1] = nperJets[i];}
	index = fMyLeptJetStat->GetConfigfrOrder(14);
	fMyLeptJetStat->NEntriesPerJetMult(index, nperJets);
	for (int i = 1; i < imax; ++i) {multable[5][i-1] += nperJets[i];}

// nber of leptons = SSem
	index = fMyLeptJetStat->GetConfigfrOrder(7);
	fMyLeptJetStat->NEntriesPerJetMult(index, nperJets);
	for (int i = 1; i < imax; ++i) {multable[6][i-1] = nperJets[i];}
	index = fMyLeptJetStat->GetConfigfrOrder(11);
	fMyLeptJetStat->NEntriesPerJetMult(index, nperJets);
	for (int i = 1; i < imax; ++i) {multable[6][i-1] += nperJets[i];}

// nber of leptons = SSee
	index = fMyLeptJetStat->GetConfigfrOrder(5);
	fMyLeptJetStat->NEntriesPerJetMult(index, nperJets);
	for (int i = 1; i < imax; ++i) {multable[7][i-1] = nperJets[i];}
	index = fMyLeptJetStat->GetConfigfrOrder(9);
	fMyLeptJetStat->NEntriesPerJetMult(index, nperJets);
	for (int i = 1; i < imax; ++i) {multable[7][i-1] += nperJets[i];}

// nber of leptons = OSmm1m
	index = fMyLeptJetStat->GetConfigfrOrder(32);
	fMyLeptJetStat->NEntriesPerJetMult(index, nperJets);
	for (int i = 1; i < imax; ++i) {multable[8][i-1] = nperJets[i];}
	index = fMyLeptJetStat->GetConfigfrOrder(33);
	fMyLeptJetStat->NEntriesPerJetMult(index, nperJets);
	for (int i = 1; i < imax; ++i) {multable[8][i-1] += nperJets[i];}

// nber of leptons = OSmm1e
	index = fMyLeptJetStat->GetConfigfrOrder(23);
	fMyLeptJetStat->NEntriesPerJetMult(index, nperJets);
	for (int i = 1; i < imax; ++i) {multable[9][i-1] = nperJets[i];}
	index = fMyLeptJetStat->GetConfigfrOrder(29);
	fMyLeptJetStat->NEntriesPerJetMult(index, nperJets);
	for (int i = 1; i < imax; ++i) {multable[9][i-1] += nperJets[i];}

// nber of leptons = OSee1m
	index = fMyLeptJetStat->GetConfigfrOrder(20);
	fMyLeptJetStat->NEntriesPerJetMult(index, nperJets);
	for (int i = 1; i < imax; ++i) {multable[10][i-1] = nperJets[i];}
	index = fMyLeptJetStat->GetConfigfrOrder(21);
	fMyLeptJetStat->NEntriesPerJetMult(index, nperJets);
	for (int i = 1; i < imax; ++i) {multable[10][i-1] += nperJets[i];}

// nber of leptons = OSee1e
	index = fMyLeptJetStat->GetConfigfrOrder(16);
	fMyLeptJetStat->NEntriesPerJetMult(index, nperJets);
	for (int i = 1; i < imax; ++i) {multable[11][i-1] = nperJets[i];}
	index = fMyLeptJetStat->GetConfigfrOrder(19);
	fMyLeptJetStat->NEntriesPerJetMult(index, nperJets);
	for (int i = 1; i < imax; ++i) {multable[11][i-1] += nperJets[i];}

// nber of leptons = SSmm1e
	index = fMyLeptJetStat->GetConfigfrOrder(28);
	fMyLeptJetStat->NEntriesPerJetMult(index, nperJets);
	for (int i = 1; i < imax; ++i) {multable[12][i-1] = nperJets[i];}
	index = fMyLeptJetStat->GetConfigfrOrder(24);
	fMyLeptJetStat->NEntriesPerJetMult(index, nperJets);
	for (int i = 1; i < imax; ++i) {multable[12][i-1] += nperJets[i];}

// nber of leptons = SSee1m
	index = fMyLeptJetStat->GetConfigfrOrder(18);
	fMyLeptJetStat->NEntriesPerJetMult(index, nperJets);
	for (int i = 1; i < imax; ++i) {multable[13][i-1] = nperJets[i];}
	index = fMyLeptJetStat->GetConfigfrOrder(26);
	fMyLeptJetStat->NEntriesPerJetMult(index, nperJets);
	for (int i = 1; i < imax; ++i) {multable[13][i-1] += nperJets[i];}

// nber of leptons = 3SSmmm
	index = fMyLeptJetStat->GetConfigfrOrder(31);
	fMyLeptJetStat->NEntriesPerJetMult(index, nperJets);
	for (int i = 1; i < imax; ++i) {multable[14][i-1] = nperJets[i];}
	index = fMyLeptJetStat->GetConfigfrOrder(34);
	fMyLeptJetStat->NEntriesPerJetMult(index, nperJets);
	for (int i = 1; i < imax; ++i) {multable[14][i-1] += nperJets[i];}

// nber of leptons = 3SSmme
	index = fMyLeptJetStat->GetConfigfrOrder(22);
	fMyLeptJetStat->NEntriesPerJetMult(index, nperJets);
	for (int i = 1; i < imax; ++i) {multable[15][i-1] = nperJets[i];}
	index = fMyLeptJetStat->GetConfigfrOrder(30);
	fMyLeptJetStat->NEntriesPerJetMult(index, nperJets);
	for (int i = 1; i < imax; ++i) {multable[15][i-1] += nperJets[i];}

// nber of leptons = 3SSmee
	index = fMyLeptJetStat->GetConfigfrOrder(17);
	fMyLeptJetStat->NEntriesPerJetMult(index, nperJets);
	for (int i = 1; i < imax; ++i) {multable[16][i-1] = nperJets[i];}
	index = fMyLeptJetStat->GetConfigfrOrder(27);
	fMyLeptJetStat->NEntriesPerJetMult(index, nperJets);
	for (int i = 1; i < imax; ++i) {multable[16][i-1] += nperJets[i];}

// nber of leptons = 3SSeee
	index = fMyLeptJetStat->GetConfigfrOrder(15);
	fMyLeptJetStat->NEntriesPerJetMult(index, nperJets);
	for (int i = 1; i < imax; ++i) {multable[17][i-1] = nperJets[i];}
	index = fMyLeptJetStat->GetConfigfrOrder(25);
	fMyLeptJetStat->NEntriesPerJetMult(index, nperJets);
	for (int i = 1; i < imax; ++i) {multable[17][i-1] += nperJets[i];}

// now fill the 2D plot
	if(fVerbose) cout << endl;
	if(fVerbose) cout << " Contents of the e/mu multiplicity 2D plot " << endl;
	for (int i = 0; i < nlept; ++i) {
		if(fVerbose) cout << "  " << lablx[i];
		for (int j = 0; j < njets; ++j) {
			fHemuMult->Fill(lablx[i], lably[j], multable[i][j]);
			if(fVerbose) cout << "  " << multable[i][j];
		}
		if(fVerbose) cout << endl;
	}

// Compute e/mu efficiency ratios
	const int neffRat = 13;
	double effRat[neffRat], deffRat[neffRat];
	const char * lablRat[neffRat] =
		{"R1l    ", "ROSll  ", "RSSll  ", "RSSel  ", "RSSml  ",
		"ROSmm1l", "ROSee1l", "ROSll1m", "ROSll1e", "RSSll1l",
		"R3SSllm", "R3SSlle", "RSS3l  "};
	float n1m = 0., n1e = 0.;
	float nOSmm = 0., nOSem = 0., nOSee = 0.;
	float nSSmm = 0., nSSem = 0., nSSee = 0.;
	float nOSmm1m = 0., nOSmm1e = 0., nOSee1m = 0., nOSee1e = 0.;
	float nSSmm1e = 0., nSSee1m = 0., nSS3m = 0.,
		nSSmme = 0., nSSeem = 0., nSS3e = 0.;
	for (int i = 1; i < imax; ++i) {
		n1m += multable[0][i-1];
		n1e += multable[1][i-1];
		nOSmm += multable[2][i-1];
		nOSem += multable[3][i-1];
		nOSee += multable[4][i-1];
		nSSmm += multable[5][i-1];
		nSSem += multable[6][i-1];
		nSSee += multable[7][i-1];
		nOSmm1m += multable[8][i-1];
		nOSmm1e += multable[9][i-1];
		nOSee1m += multable[10][i-1];
		nOSee1e += multable[11][i-1];
		nSSmm1e += multable[12][i-1];
		nSSee1m += multable[13][i-1];
		nSS3m += multable[14][i-1];
		nSSmme += multable[15][i-1];
		nSSeem += multable[16][i-1];
		nSS3e += multable[17][i-1];
	}

	for (int i = 0; i < neffRat; ++i) {
		effRat[i] = 0.;
		deffRat[i] = 0.;
	}
	if (n1m != 0.) {
		effRat[0] = n1e / n1m;
		if (n1e != 0.) {
			deffRat[0] = sqrt(1./n1e + 1./n1m) * effRat[0];
		}
	}
	float nOSZee = nOSee-0.5*nOSem;
	float nOSZmm = nOSmm-0.5*nOSem;
	if (nOSZmm > 0. && nOSZmm*nOSZee > 0.) {
		effRat[1] = sqrt( nOSZee / nOSZmm );
		if (nOSZee > 0.) {
			deffRat[1] = sqrt((nOSee+0.25*nOSem)/(nOSZee*nOSZee)
				+ (nOSmm+0.25*nOSem)/(nOSZmm*nOSZmm) ) * effRat[1];
		}
	}
	if (nSSmm != 0. && nSSee*nSSmm > 0.) {
		effRat[2] = sqrt( nSSee / nSSmm );
		if (nSSee != 0.) {
			deffRat[2] = 0.5 * sqrt(1./nSSee + 1./nSSmm) * effRat[2];
		}
	}
	if (nSSem != 0.) {
		effRat[3] = 2. * nSSee / nSSem;
		if (nSSee != 0.) {
			deffRat[3] = 2. * sqrt(1./nSSee + 1./nSSem) * effRat[3];
		}
	}
	if (nSSmm != 0.) {
		effRat[4] = 0.5 * nSSem / nSSmm;
		if (nSSem != 0.) {
			deffRat[4] = 0.5 * sqrt(1./nSSem + 1./nSSmm) * effRat[4];
		}
	}
	if (nOSmm1m != 0.) {
		effRat[5] = nOSmm1e / nOSmm1m;
		if (nOSmm1e != 0.) {
			deffRat[5] = sqrt(1./nOSmm1e + 1./nOSmm1m) * effRat[5];
		}
	}
	if (nOSee1m != 0.) {
		effRat[6] = nOSee1e / nOSee1m;
		if (nOSee1e != 0.) {
			deffRat[6] = sqrt(1./nOSee1e + 1./nOSee1m) * effRat[6];
		}
	}
	if (nOSmm1m != 0. && nOSee1m*nOSmm1m > 0.) {
		effRat[7] = sqrt( nOSee1m / nOSmm1m );
		if (nOSee1m != 0.) {
			deffRat[7] = 0.5 * sqrt(1./nOSee1m + 1./nOSmm1m) * effRat[7];
		}
	}
	if (nOSmm1e != 0. && nOSee1e*nOSmm1e > 0.) {
		effRat[8] = sqrt( nOSee1e / nOSmm1e );
		if (nOSee1e != 0.) {
			deffRat[8] = 0.5 * sqrt(1./nOSee1e + 1./nOSmm1e) * effRat[8];
		}
	}
	if (nSSmm1e != 0.) {
		effRat[9] = nSSee1m / nSSmm1e;
		if (nSSee1m != 0.) {
			deffRat[9] = sqrt(1./nSSee1m + 1./nSSmm1e) * effRat[9];
		}
	}
	if (nSS3m != 0. && nSSeem*nSS3m > 0.) {
		effRat[10] = sqrt( nSSeem / nSS3m );
		if (nSSeem != 0.) {
			deffRat[10] = 0.5 * sqrt(1./nSSeem + 1./nSS3m) * effRat[10];
		}
	}
	if (nSSmme != 0. && nSS3e*nSSmme > 0.) {
		effRat[11] = sqrt( nSS3e / nSSmme );
		if (nSS3e != 0.) {
			deffRat[11] = 0.5 * sqrt(1./nSS3e + 1./nSSmme) * effRat[11];
		}
	}
	if (nSS3m != 0. && nSS3e*nSS3m > 0.) {
		effRat[12] = pow((nSS3e / nSS3m), (float)(1./3.) );
		if (nSS3e != 0.) {
			deffRat[12] = sqrt(1./nSS3e + 1./nSS3m) / 3. * effRat[12];
		}
	}

// Now fill the Profile histogram
	if(fVerbose){
		cout << endl;
		cout << " Efficiency ratios e/mu" << endl;
		cout << "  for lepton = l, the ratio e/mu is implied" << endl;
		cout << "  for ratios " << lablRat[1] << ", " << lablRat[2] << ", "
		<< lablRat[7] << ", " << lablRat[8] << ", "
		<< lablRat[10] << ", " << lablRat[11]
		<< " a sqrt is implied" << endl;
		cout << "  for ratio " << lablRat[10] << " a 3rd root is implied" << endl;
	}
	for (int i = 0; i < neffRat; ++i) {
		fHemuEff-> GetXaxis()->SetBinLabel(i+1, lablRat[i]);
		fHemuEff->SetBinContent(i+1, effRat[i]);
		fHemuEff->SetBinError(i+1, deffRat[i]);
		if(fVerbose) cout << "  " << lablRat[i] << " " << effRat[i]
			<< " +- " << deffRat[i] << endl;
	}
	return;
}

