#include "TreeReader.hh"
using namespace std;

TreeReader::TreeReader(TTree *tree, int flag) : TreeClassBase(tree){
	fDiLep = false;
	fMPHist = false;
	fSignHist = false;

	if((flag/100)%10) fDiLep    = true;
	if((flag/10)%10)  fMPHist   = true;
	if((flag/1)%10)   fSignHist = true;

	fTlat = new TLatex();
	DefStyle();
	gROOT->SetStyle("ETHStyle");
	gROOT->ForceStyle();
}

TreeReader::~TreeReader(){
	if(!fChain) cout << "no chain!" << endl;
}

void TreeReader::DefStyle(){
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
}

// Method called before starting the event loop
void TreeReader::BeginJob(){
	if(fDiLep)    InitDiLepTree();
	if(fMPHist)   BookMPHistos();
	if(fSignHist) BookSignHists();
}

// Method called after finishing the event loop
void TreeReader::EndJob(){
	if(fDiLep)    WriteDiLepTree();
	if(fMPHist)   PrintMPOutput();
	if(fSignHist) WriteSignHists();
}

// Method for looping over the tree
void TreeReader::Loop(){
	Long64_t nbytes = 0, nb = 0;
	Long64_t nentries = fChain->GetEntries();

	for(Long64_t jentry=0; jentry<nentries;jentry++){
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		nb = fChain->GetEntry(jentry);   nbytes += nb;
		if (jentry%10000 == 0) cout << ">>> Processing event # " << jentry << endl;

		// Put here any method that needs to be called once every event
		if(fDiLep)    FillDiLepTree();
		if(fMPHist)   FillMPHistos();
		if(fSignHist) for(size_t i = 0; i < 5; ++i) FillSignHists(i);
	}
}

///////////////////////////////////////////////////////////////////////////////////////////////
void TreeReader::SetOutputDir(TString dir){
	if(!dir.EndsWith("/")) dir += "/";
	fOutputDir = dir;
	// Create directory if needed
	//  >> NOTE: This function needs to be called before the booking functions!
	char cmd[100];
	sprintf(cmd,"mkdir -p %s", fOutputDir.Data());
	system(cmd);
};

///////////////////////////////////////////////////////////////////////////////////////////////
// Significance Plots Stuff ///////////////////////////////////////////////////////////////////
void TreeReader::BookSignHists(const char* filename){
	fSignHistsFile = new TFile(fOutputDir + TString(filename), "RECREATE");
	fSignHistsFile->cd();
	fNBinsEta[0] = 10;
	fNBinsEta[1] = 10;
	fNBinsEta[2] = 10;
	fNBinsEta[3] = 20;
	fNBinsEta[4] = 20;
	fNBinsPhi = 20;
	const float etamin[5] = {-3.0, -2.5, -5.0, -10.0, -10.0};
	const float etamax[5] = { 3.0,  2.5,  5.0,  10.0,  10.0};
	const float phimin = -3.1416, phimax = 3.1416;
	TString partname[5];
	partname[0] = "electrons";
	partname[1] = "muons";
	partname[2] = "jets";
	partname[3] = "ecalMET";
	partname[4] = "hcalMET";

	for(int i=0; i<5; i++){
		fH_ptdev[i]    = new TH2D(Form("h_ptdev%d", i),    Form("pT deviation in eta vs phi for %s", partname[i].Data()),  fNBinsEta[i], etamin[i], etamax[i], fNBinsPhi, phimin, phimax);
		fH_ptsum[i]    = new TH2D(Form("h_ptsum%d", i),    Form("pT sum in eta vs phi for %s", partname[i].Data()),        fNBinsEta[i], etamin[i], etamax[i], fNBinsPhi, phimin, phimax);
		fH_pt2sum[i]   = new TH2D(Form("h_pt2sum%d", i),   Form("pT square sum in eta vs phi for %s", partname[i].Data()), fNBinsEta[i], etamin[i], etamax[i], fNBinsPhi, phimin, phimax);
		fH_ptevt[i]    = new TH2I(Form("h_ptevt%d", i),    Form("Number of ev. in eta vs phi for %s", partname[i].Data()), fNBinsEta[i], etamin[i], etamax[i], fNBinsPhi, phimin, phimax);
		fH_ptavg[i]    = new TH2D(Form("h_ptavg%d", i),    Form("pT average in eta vs phi for %s", partname[i].Data()),    fNBinsEta[i], etamin[i], etamax[i], fNBinsPhi, phimin, phimax);
		fH_ptsumeta[i] = new TH1D(Form("h_ptsumeta%d", i), Form("pT sum in eta for %s", partname[i].Data()),               fNBinsEta[i], etamin[i], etamax[i]);
		fH_ptevteta[i] = new TH1I(Form("h_ptevteta%d", i), Form("Number of ev. in eta for %s", partname[i].Data()),        fNBinsEta[i], etamin[i], etamax[i]);
		fH_ptdev[i]->SetXTitle("#eta");
		fH_ptsum[i]->SetXTitle("#eta");
		fH_ptevt[i]->SetXTitle("#eta");
		fH_ptavg[i]->SetXTitle("#eta");
		fH_ptdev[i]->SetYTitle("#phi");
		fH_ptsum[i]->SetYTitle("#phi");
		fH_ptevt[i]->SetYTitle("#phi");
		fH_ptavg[i]->SetYTitle("#phi");

		fH_ptdev[i]->SetStats(false);
		fH_ptavg[i]->SetStats(false);
		fH_ptsum[i]->SetStats(false);
		fH_ptevt[i]->SetStats(false);
	}
}

void TreeReader::FillSignHists(Int_t part){
// Makes 2D plots of pT 
// in eta vs phi for
//  part = 0: electrons
//       = 1: muons
//       = 2: jets (default)
//       = 3: ECAL recoil
//       = 4: HCAL recoil

	fSignHistsFile->cd();
	double EcalEta(0.);
	double HcalEta(0.);

	// electrons
	if (part == 0) {
		for (int ip = 0; ip < NEles; ++ ip) {
			if(IsGoodObj(ip, &fElCuts) == false) continue;
			fH_ptsum[part]->Fill(ElEta[ip], ElPhi[ip], ElPt[ip]);
			fH_pt2sum[part]->Fill(ElEta[ip], ElPhi[ip], ElPt[ip]*ElPt[ip]);
			fH_ptevt[part]->Fill(ElEta[ip], ElPhi[ip]);
			fH_ptsumeta[part]->Fill(ElEta[ip], ElPt[ip]);
			fH_ptevteta[part]->Fill(ElEta[ip]);
		}
	}
	// muons
	else if (part == 1) {
		for (int ip = 0; ip < NMus; ++ ip) {
			if(IsGoodObj(ip, &fMuCuts) == false) continue;
			fH_ptsum[part]->Fill(MuEta[ip], MuPhi[ip], MuPt[ip]);
			fH_pt2sum[part]->Fill(MuEta[ip], MuPhi[ip], MuPt[ip]*MuPt[ip]);
			fH_ptevt[part]->Fill(MuEta[ip], MuPhi[ip]);
			fH_ptsumeta[part]->Fill(MuEta[ip], MuPt[ip]);
			fH_ptevteta[part]->Fill(MuEta[ip]);
		}
	}
	// jets
	else if (part == 2) {
		for (int ip = 0; ip < NJets; ++ ip) {
			if(IsGoodObj(ip, &fJetCuts) == false) continue;
			fH_ptsum[part]->Fill(JEta[ip], JPhi[ip], JPt[ip]);
			fH_pt2sum[part]->Fill(JEta[ip], JPhi[ip], JPt[ip]*JPt[ip]);
			fH_ptevt[part]->Fill(JEta[ip], JPhi[ip]);
			fH_ptsumeta[part]->Fill(JEta[ip], JPt[ip]);
			fH_ptevteta[part]->Fill(JEta[ip]);
		}
	}
	// ECAL MET
	else if (part == 3) {
		// this is to protect against NAN
		if (ECALEsumx != ECALEsumx) {return;}
		EcalEta = getEta(ECALEsumx, ECALEsumy, ECALEsumz);
		fH_ptsum[part]->Fill(EcalEta, ECALMETPhi, ECALMET);
		fH_pt2sum[part]->Fill(EcalEta, ECALMETPhi, ECALMET*ECALMET);
		fH_ptevt[part]->Fill(EcalEta, ECALMETPhi);
		fH_ptsumeta[part]->Fill(EcalEta, ECALMET);
		fH_ptevteta[part]->Fill(EcalEta);
	}
	// HCAL MET
	else if (part == 4) {
	// this is to protect against NAN
		if (HCALEsumx != HCALEsumx) {return;}
		HcalEta = getEta(HCALEsumx, HCALEsumy, HCALEsumz);
		fH_ptsum[part]->Fill(HcalEta, HCALMETPhi, HCALMET);
		fH_pt2sum[part]->Fill(HcalEta, HCALMETPhi, HCALMET*HCALMET);
		fH_ptevt[part]->Fill(HcalEta, HCALMETPhi);
		fH_ptsumeta[part]->Fill(HcalEta, HCALMET);
		fH_ptevteta[part]->Fill(HcalEta);
	}
}

void TreeReader::WriteSignHists(){
	fSignHistsFile->cd();
	float zminDev = -4., zmaxDev = 4.;

	// Calculate averages and deviations
	for(size_t part = 0; part < 5; ++part){
		for (int i = 0; i < fNBinsEta[part]; ++i) {
			float ptaverEta = 0.;
			if (fH_ptevteta[part]->GetBinContent(i+1) > 0) {
				ptaverEta = fH_ptsumeta[part]->GetBinContent(i+1) / (float)fH_ptevteta[part]->GetBinContent(i+1);
			}
			for (int j = 0; j < fNBinsPhi; ++j) {
				float ptaver = 0.;
				float ptdev = 0.;
				int nptij = (int)fH_ptevt[part]->GetBinContent(i+1, j+1);
				if (nptij > 0) {
					float ptsumij  = fH_ptsum[part]->GetBinContent(i+1, j+1);
					float pt2sumij = fH_pt2sum[part]->GetBinContent(i+1, j+1);
					ptaver = ptsumij / (float)nptij;
					float pterr = sqrt(pt2sumij - (float)nptij*ptaver*ptaver);
					if (pterr <= 0.) {pterr = 0.1;}
					ptdev = (ptsumij - (float)nptij*ptaverEta) / pterr;
				}
				if (ptdev > zmaxDev) {ptdev = zmaxDev;}
				if (ptdev < zminDev) {ptdev = zminDev;}
				fH_ptdev[part]->SetBinContent(i+1, j+1, ptdev);
				fH_ptavg[part]->SetBinContent(i+1, j+1, ptaver);
			}
		}
	}
		
	TString pnames[5];
	pnames[0] = "el";
	pnames[1] = "mu";
	pnames[2] = "jet";
	pnames[3] = "eMET";
	pnames[4] = "hMET";
	
	TString pnamel[5];
	pnamel[0] = "Electrons";
	pnamel[1] = "Muons";
	pnamel[2] = "Jets";
	pnamel[3] = "ECALMET";
	pnamel[4] = "HCALMET";
	
	fTlat->SetTextColor(kBlack);
	fTlat->SetNDC(kTRUE);
	fTlat->SetTextSize(0.04);

	// TStyle *style = gROOT->GetStyle("ETHStyle");
	// const int ncol1 = 11;
	// int colors1[ncol1]= {10, 16, 5, 28, 29, 8, 4, 9, 45, 46, 2};
	// const int ncol2 = 8;
	// // int colors2[ncol2]= {6, 28, 5, 29, 8, 7, 4, 2};
	// int colors2[ncol2]= {4, 28, 5, 29, 8, 7, 47, 2}; // new colors from luc: 21/10/09
	// const int ncont2 = 9;
	// double contours2[ncont2] = { -4.,-3.,-2.,-1., 0., 1., 2., 3., 4.};

	// Plot the histograms
	TString subdir = "SignificancePlots";
	TCanvas *canv;
	for(size_t i = 0; i < 5; ++i){
		// gStyle->SetPalette(ncol2, colors2);
		TString canvname = "PTDev_" + pnames[i];
		TString canvtitle = "pT Deviation for " + pnamel[i];
		canv = new TCanvas(canvname, canvtitle, 0, 0, 900, 700);
		canv->SetRightMargin(0.15);
		// fH_ptdev[i]->SetContour(ncont2, contours2);
		fH_ptdev[i]->SetMinimum(-4);
		fH_ptdev[i]->SetMaximum(4);
		fH_ptdev[i]->DrawCopy("colz");
		fTlat->DrawLatex(0.11,0.92, canvtitle);
		printPNG(canv, fTag + "_" + canvname, fOutputDir+subdir);
		printEPS(canv, fTag + "_" + canvname, fOutputDir+subdir);

		// gStyle->SetPalette(ncol1, colors1);

		canvname = "PTSum_" + pnames[i];
		canvtitle = "pT sum for " + pnamel[i];
		canv = new TCanvas(canvname, canvtitle, 0, 0, 900, 700);
		canv->SetRightMargin(0.15);
		fH_ptsum[i]->SetMinimum(0);
		fH_ptsum[i]->DrawCopy("lego2 z");
		fTlat->DrawLatex(0.11,0.92, canvtitle);
		printPNG(canv, fTag + "_" + canvname, fOutputDir+subdir);
		printEPS(canv, fTag + "_" + canvname, fOutputDir+subdir);

		canvname = "PTEvt_" + pnames[i];
		canvtitle = "Event multiplicity for " + pnamel[i];
		canv = new TCanvas(canvname, canvtitle, 0, 0, 900, 700);
		canv->SetRightMargin(0.15);
		fH_ptevt[i]->SetMinimum(0);
		fH_ptevt[i]->DrawCopy("lego2 z");
		fTlat->DrawLatex(0.11,0.92, canvtitle);
		printPNG(canv, fTag + "_" + canvname, fOutputDir+subdir);
		printEPS(canv, fTag + "_" + canvname, fOutputDir+subdir);

		canvname = "PTAvg_" + pnames[i];
		canvtitle = "pT Average for " + pnamel[i];
		canv = new TCanvas(canvname, canvtitle, 0, 0, 900, 700);
		canv->SetRightMargin(0.15);
		fH_ptavg[i]->SetMinimum(0);
		fH_ptavg[i]->DrawCopy("lego2 z");
		fTlat->DrawLatex(0.11,0.92, canvtitle);
		printPNG(canv, fTag + "_" + canvname, fOutputDir+subdir);
		printEPS(canv, fTag + "_" + canvname, fOutputDir+subdir);
	}

	// Write the histograms
	for(size_t i = 0; i < 5; ++i){
		fH_ptdev[i]->Write();
		fH_ptsum[i]->Write();
		fH_ptevt[i]->Write();
		fH_ptavg[i]->Write();
	}
	fSignHistsFile->Close();
}

///////////////////////////////////////////////////////////////////////////////////////////////
// Multiplicity Plots Stuff ///////////////////////////////////////////////////////////////////
void TreeReader::BookMPHistos(const char* filename){
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

void TreeReader::FillMPHistos(){
	// for (int i=0 ; i<NJets; i++){
	// cout<< "jet with pt = "<< JPt[i] << endl;
	// }

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
	for(int i=0; i< NEles; ++i){
		if(IsGoodObj(i, &fElCuts) == false) continue;
		ILept++;
		LeptPt[ILept]=ElPt[i];
		LeptEta[ILept]=ElEta[i];
		if(ElCharge[i]>0) {
			LeptCat[ILept]=1;
		} else {
			LeptCat[ILept]=2;
		}
	}
	for(int i=0; i< NMus; ++i){
		if(IsGoodObj(i, &fMuCuts) == false) continue;
		ILept++;
		LeptPt[ILept]=MuPt[i];
		LeptEta[ILept]=MuEta[i];
		if(MuCharge[i]>0) {
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
	for(int i=0; i< NJets; ++i){
		if(IsGoodObj(i, &fJetCuts)) NQJets++;
	}

// convert the lepton config into the index and count
	fMyLeptJetStat->FillLeptJetStat(LeptCat, NQJets, 0);
}

void TreeReader::PrintMPOutput(){
	fMyLeptJetStat->FillShortTables();
	fMyLeptJetStat->PrintConfigs ();
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
	printPNG(canv, fTag + "_ljMult", fOutputDir + subdir);
	printEPS(canv, fTag + "_ljMult", fOutputDir + subdir);

	canvtitle = "e/mu multiplicity";
	canv = new TCanvas("emuMult", canvtitle , 0, 0, 900, 700);
	canv->SetRightMargin(0.15);
   gPad->SetTheta(50); 
   gPad->SetPhi(240);
	fHemuMult->SetMinimum(0);
	fHemuMult->DrawCopy("colz");
	// fHemuMult->DrawCopy("lego2 Z");
	fTlat->DrawLatex(0.11,0.92, canvtitle);
	printPNG(canv, fTag + "_emuMult", fOutputDir + subdir);
	printEPS(canv, fTag + "_emuMult", fOutputDir + subdir);

	canvtitle = "e/mu Efficiency";
	canv = new TCanvas("emuEffic", canvtitle , 0, 0, 900, 700);
	canv->SetRightMargin(0.15);
   gPad->SetTheta(50); 
   gPad->SetPhi(240);
	fHemuEff->DrawCopy();
	fTlat->DrawLatex(0.11,0.92, canvtitle);
	printPNG(canv, fTag + "_emuEffic", fOutputDir+subdir);
	printEPS(canv, fTag + "_emuEffic", fOutputDir+subdir);


	fMPHistFile->cd();
	fHljMult->Write();
	fHemuMult->Write();
	fHemuEff->Write();
	fMPHistFile->Close();	
}

void TreeReader::PlotMPSummary(){
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
	cout << endl;
	cout << " Contents of the lepton/jets multiplicity 2D plot " << endl;
	for (int i = 0; i < nlept; ++i) {
		cout << "  " << lablx[i];
		for (int j = 0; j < njets; ++j) {
			fHljMult->Fill(lablx[i], lably[j], multable[i][j]);
			cout << "  " << multable[i][j];
		}
		cout << endl;
	}
	return;
}

void TreeReader::PlotMPEffic(){
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
	cout << endl;
	cout << " Contents of the e/mu multiplicity 2D plot " << endl;
	for (int i = 0; i < nlept; ++i) {
		cout << "  " << lablx[i];
		for (int j = 0; j < njets; ++j) {
			fHemuMult->Fill(lablx[i], lably[j], multable[i][j]);
			cout << "  " << multable[i][j];
		}
		cout << endl;
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
	cout << endl;
	cout << " Efficiency ratios e/mu" << endl;
	cout << "  for lepton = l, the ratio e/mu is implied" << endl;
	cout << "  for ratios " << lablRat[1] << ", " << lablRat[2] << ", "
		<< lablRat[7] << ", " << lablRat[8] << ", "
		<< lablRat[10] << ", " << lablRat[11]
		<< " a sqrt is implied" << endl;
	cout << "  for ratio " << lablRat[10] 
		<< " a 3rd root is implied" << endl;
	for (int i = 0; i < neffRat; ++i) {
		fHemuEff-> GetXaxis()->SetBinLabel(i+1, lablRat[i]);
		fHemuEff->SetBinContent(i+1, effRat[i]);
		fHemuEff->SetBinError(i+1, deffRat[i]);
		cout << "  " << lablRat[i] << " " << effRat[i]
			<< " +- " << deffRat[i] << endl;
	}
	return;
}

///////////////////////////////////////////////////////////////////////////////////////////////
// Di Lepton Tree Stuff ///////////////////////////////////////////////////////////////////////
void TreeReader::InitDiLepTree(const char* filename){
	fDiLepTreeFile = new TFile(fOutputDir + TString(filename), "RECREATE");
	fDiLepTree = new TTree("Analysis", "NanoAnalysisTree");
	fDiLepTree->Branch("Run",           &fTRunNumber,      "Run/I");
	fDiLepTree->Branch("Event",         &fTEventNumber,    "Event/I");
	fDiLepTree->Branch("LumiSec",       &fTLumiSection,    "LumiSec/I");
	fDiLepTree->Branch("NQMus",         &fTNqualmu,        "NQMus/I");
	fDiLepTree->Branch("Mu1Ch",         &fTMu1charge,      "Mu1Ch/I");
	fDiLepTree->Branch("Mu2Ch",         &fTMu2charge,      "Mu2Ch/I");
	fDiLepTree->Branch("Mu1Pt",         &fTMu1pt,          "Mu1Pt/D");
	fDiLepTree->Branch("Mu2Pt",         &fTMu2pt,          "Mu2Pt/D");
	fDiLepTree->Branch("Mu1Eta",        &fTMu1eta,         "Mu1Eta/D");
	fDiLepTree->Branch("Mu2Eta",        &fTMu2eta,         "Mu2Eta/D");
	fDiLepTree->Branch("Mu1Iso",        &fTMu1iso,         "Mu1Iso/D");
	fDiLepTree->Branch("Mu2Iso",        &fTMu2iso,         "Mu2Iso/D");
	fDiLepTree->Branch("Mu1D0",         &fTMu1d0,          "Mu1D0/D");
	fDiLepTree->Branch("Mu2D0",         &fTMu2d0,          "Mu2D0/D");
	fDiLepTree->Branch("Mu1NTkHits",    &fTMu1ntkhits,     "Mu1NTkHits/D");
	fDiLepTree->Branch("Mu2NTkHits",    &fTMu2ntkhits,     "Mu2NTkHits/D");
	fDiLepTree->Branch("MuMInv",        &fTMuminv,         "MuMInv/D");
	fDiLepTree->Branch("MuMT2_50",      &fTMumt2_50,       "MuMT2_50/D");
	fDiLepTree->Branch("MuMT2_100",     &fTMumt2_100,      "MuMT2_100/D");
	fDiLepTree->Branch("NQEls",         &fTNqualel,        "NQEls/I");
	fDiLepTree->Branch("El1Ch",         &fTEl1charge,      "El1Ch/I");
	fDiLepTree->Branch("El2Ch",         &fTEl2charge,      "El2Ch/I");
	fDiLepTree->Branch("El1Pt",         &fTEl1pt,          "El1Pt/D");
	fDiLepTree->Branch("El2Pt",         &fTEl2pt,          "El2Pt/D");
	fDiLepTree->Branch("El1Eta",        &fTEl1eta,         "El1Eta/D");
	fDiLepTree->Branch("El2Eta",        &fTEl2eta,         "El2Eta/D");
	fDiLepTree->Branch("El1Iso",        &fTEl1iso,         "El1Iso/D");
	fDiLepTree->Branch("El2Iso",        &fTEl2iso,         "El2Iso/D");
	fDiLepTree->Branch("El1D0",         &fTEl1d0,          "El1D0/D");
	fDiLepTree->Branch("El2D0",         &fTEl2d0,          "El2D0/D");
	fDiLepTree->Branch("ElMInv",        &fTElminv,         "ElMInv/D");
	fDiLepTree->Branch("ElMT2_50",      &fTElmt2_50,       "ElMT2_50/D");
	fDiLepTree->Branch("ElMT2_100",     &fTElmt2_100,      "ElMT2_100/D");
}

void TreeReader::FillDiLepTree(){
	ResetDiLepTree();
	// Select events with either 2 muons or 2 electrons

	if(!IsGoodEvt(&fEvtSelCuts)) return;

	if( NMus < 2 && NEles < 2 ) return;
	vector<int> qualMuInd;
	for(size_t imu = 0; imu < NMus; ++imu){
		// Muon selection
		if(IsGoodObj(imu, &fMuCuts) == false) continue;
		qualMuInd.push_back(imu);
	}
	vector<int> qualElInd;
	for(size_t iel = 0; iel < NEles; ++iel){
		// Electron selection
		if(IsGoodObj(iel, &fElCuts) == false) continue;
		qualElInd.push_back(iel);
	}

	int nqmus = qualMuInd.size();
	int nqels = qualElInd.size();

	// Select events with either 2 qualified muons or electrons
	if(nqmus < 2 && nqels < 2) return;

	fTNqualmu = nqmus;
	fTNqualel = nqels;

	fMT2 = new Davismt2();
	double pa[3], pb[3], pmiss[3];

	int lep1index(-1), lep2index(-1);
	// DiMuons:
	if(nqmus > 1){
		// Find the two hardest muons
		double maxmupt = 0.;
		for(size_t i = 0; i < nqmus; ++i){
			int index = qualMuInd[i];
			if(MuPt[index] < maxmupt) continue;
			maxmupt = MuPt[index];
			lep1index = index;
		}
		maxmupt = 0.;
		for(size_t i = 0; i < nqmus; ++i){
			int index = qualMuInd[i];
			if(index == lep1index) continue;
			if(MuPt[index] < maxmupt) continue;
			maxmupt = MuPt[index];
			lep2index = index;
		}

		fTMu1pt       = MuPt[lep1index];
		fTMu2pt       = MuPt[lep2index];
		fTMu1charge   = MuCharge[lep1index];
		fTMu2charge   = MuCharge[lep2index];
		fTMu1eta      = MuEta[lep1index];
		fTMu2eta      = MuEta[lep2index];
		fTMu1iso      = MuIso[lep1index];
		fTMu2iso      = MuIso[lep2index];
		fTMu1d0       = MuD0BS[lep1index];
		fTMu2d0       = MuD0BS[lep2index];
		fTMu1ntkhits  = MuNTkHits[lep1index];
		fTMu2ntkhits  = MuNTkHits[lep2index];

		// Calculate invariant mass of the pair			
		TLorentzVector p1(MuPx[lep1index], MuPy[lep1index], MuPz[lep1index], MuE[lep1index]);
		TLorentzVector p2(MuPx[lep2index], MuPy[lep2index], MuPz[lep2index], MuE[lep2index]);
		fTMuminv = (p1+p2).Mag();

		pa[0] = 0.; pb[0] = 0.; pmiss[0] = 0.;
		pa[1] = MuPx[lep1index];
		pa[2] = MuPz[lep1index];
		pb[1] = MuPx[lep2index];
		pb[2] = MuPy[lep2index];
		pmiss[1] = TCMETpx;
		pmiss[2] = TCMETpy;
		fMT2->set_momenta(pa, pb, pmiss);
		fMT2->set_mn(50.);
		fTMumt2_50 = fMT2->get_mt2();
		fMT2->set_mn(100.);
		fTMumt2_100 = fMT2->get_mt2();			
	}else if(nqmus > 0){ // Meaning there is only one muon
		fTMu1pt       = MuPt[qualMuInd[0]];
		fTMu1charge   = MuCharge[qualMuInd[0]];
		fTMu1eta      = MuEta[qualMuInd[0]];
		fTMu1iso      = MuIso[qualMuInd[0]];
		fTMu1d0       = MuD0BS[qualMuInd[0]];
		fTMu1ntkhits  = MuNTkHits[qualMuInd[0]];
	}

	lep1index = -1;
	lep2index = -1;

	// DiElectrons
	if(nqels > 1){
		// Find the two hardest muons
		double maxelpt = 0.;
		for(size_t i = 0; i < nqels; ++i){
			int index = qualElInd[i];
			if(ElPt[index] < maxelpt) continue;
			maxelpt = ElPt[index];
			lep1index = index;
		}
		maxelpt = 0.;
		for(size_t i = 0; i < nqels; ++i){
			int index = qualElInd[i];
			if(index == lep1index) continue;
			if(ElPt[index] < maxelpt) continue;
			maxelpt = ElPt[index];
			lep2index = index;
		}

		fTEl1charge   = ElCharge[lep1index];
		fTEl2charge   = ElCharge[lep2index];
		fTEl1pt       = ElPt[lep1index];
		fTEl2pt       = ElPt[lep2index];
		fTEl1eta      = ElEta[lep1index];
		fTEl2eta      = ElEta[lep2index];
		fTEl1iso      = ElIso[lep1index];
		fTEl2iso      = ElIso[lep2index];
		fTEl1d0       = ElD0BS[lep1index];
		fTEl2d0       = ElD0BS[lep2index];

		// Calculate invariant mass of the pair
		TLorentzVector p1(ElPx[lep1index], ElPy[lep1index], ElPz[lep1index], ElE[lep1index]);
		TLorentzVector p2(ElPx[lep2index], ElPy[lep2index], ElPz[lep2index], ElE[lep2index]);
		fTElminv = (p1+p2).Mag();

		pa[0] = 0.; pb[0] = 0.; pmiss[0] = 0.;
		pa[1] = ElPx[lep1index];
		pa[2] = ElPz[lep1index];
		pb[1] = ElPx[lep2index];
		pb[2] = ElPy[lep2index];
		pmiss[1] = TCMETpx;
		pmiss[2] = TCMETpy;
		fMT2->set_momenta(pa, pb, pmiss);
		fMT2->set_mn(50.);
		fTElmt2_50 = fMT2->get_mt2();
		fMT2->set_mn(100.);
		fTElmt2_100 = fMT2->get_mt2();
	}else if (nqels>0){ // Meaning there is only one electron
		fTEl1pt       = ElPt[qualElInd[0]];
		fTEl1charge   = ElCharge[qualElInd[0]];
		fTEl1eta      = ElEta[qualElInd[0]];
		fTEl1iso      = ElIso[qualElInd[0]];
		fTEl1d0       = ElD0BS[qualElInd[0]];
	}
	fDiLepTree->Fill();
	delete fMT2;
}

void TreeReader::ResetDiLepTree(){
	fTRunNumber   = -999;
	fTEventNumber = -999;
	fTLumiSection = -999;
	fTNqualmu     = -999;
	fTMu1charge   = -999;
	fTMu2charge   = -999;
	fTMu1pt       = -999.99;
	fTMu2pt       = -999.99;
	fTMu1eta      = -999.99;
	fTMu2eta      = -999.99;
	fTMu1iso      = -999.99;
	fTMu2iso      = -999.99;
	fTMu1d0       = -999.99;
	fTMu2d0       = -999.99;
	fTMu1ntkhits  = -999.99;
	fTMu2ntkhits  = -999.99;
	fTMuminv      = -999.99;
	fTMumt2_50    = -999.99;
	fTMumt2_100   = -999.99;
	fTNqualel     = -999;
	fTEl1charge   = -999;
	fTEl2charge   = -999;
	fTEl1pt       = -999.99;
	fTEl2pt       = -999.99;
	fTEl1eta      = -999.99;
	fTEl2eta      = -999.99;
	fTEl1iso      = -999.99;
	fTEl2iso      = -999.99;
	fTEl1d0       = -999.99;
	fTEl2d0       = -999.99;
	fTElminv      = -999.99;
	fTElmt2_50    = -999.99;
	fTElmt2_100   = -999.99;
}

void TreeReader::WriteDiLepTree(){
	fDiLepTreeFile->cd();
	fDiLepTree->Write();
	fDiLepTreeFile->Close();
}

///////////////////////////////////////////////////////////////////////////////////////////////
// Cut Stuff //////////////////////////////////////////////////////////////////////////////////
void TreeReader::ReadObjCuts(const char* filename){
// Fills the vectors containing object quality cuts for muons, electrons and jets
	ifstream IN(filename);
	float inf(0.), sup(0.);
	char readbuff[200], branch[200];
	cout << "Reading object selection from " << filename << endl;
	// Loop over lines of datafile
	while( IN.getline(readbuff, 200, '\n') ){
		if (readbuff[0] == '#') {continue;} // Skip lines commented with '#'
		sscanf(readbuff, "%s %f %f", branch, &inf, &sup);
		TBranch *tempbranch = NULL;
		if((tempbranch = fChain->FindBranch(branch)) == NULL){
			cout << "TreeReader::ReadCuts ==> Branch \"" << branch << "\" not found, continuing... " << endl;
			continue;
		}
		Cut cut;
		if (branch[0] == 'M'){
			cut.branch = tempbranch;
			cut.upperbound = sup;
			cut.lowerbound = inf;
			cout << "  Cutting  " << (cut.branch)->GetName() << "\t\tbetween " << cut.lowerbound << " and " << cut.upperbound << endl;
			fMuCuts.push_back(cut);
		}
		else if (branch[0] == 'E'){
			cut.branch = tempbranch;
			cut.upperbound = sup;
			cut.lowerbound = inf;
			cout << "  Cutting  " << (cut.branch)->GetName() << "\t\tbetween " << cut.lowerbound << " and " << cut.upperbound << endl;
			fElCuts.push_back(cut);
		}
		else if (branch[0] == 'J'){
			cut.branch = tempbranch;
			cut.upperbound = sup;
			cut.lowerbound = inf;
			cout << "  Cutting  " << (cut.branch)->GetName() << "\t\tbetween " << cut.lowerbound << " and " << cut.upperbound << endl;
			fJetCuts.push_back(cut);
		}
	}
	cout << " ----------- " << endl;
}

void TreeReader::ReadEvtSel(const char* filename){
// Fills the vectors containing event selection cuts
	ifstream IN(filename);
	float inf(0.), sup(0.);
	char readbuff[200], branch[200];
	cout << "Reading event selection from " << filename << endl;
	// Loop over lines of datafile
	while( IN.getline(readbuff, 200, '\n') ){
		if (readbuff[0] == '#') {continue;} // Skip lines commented with '#'
		sscanf(readbuff, "%s %f %f", branch, &inf, &sup);
		TBranch *tempbranch = NULL;
		if((tempbranch = fChain->FindBranch(branch)) != NULL){
			Cut cut;
			cut.branch = tempbranch;
			cut.upperbound = sup;
			cut.lowerbound = inf;
			cout << "  Cutting  " << (cut.branch)->GetName() << "\t\tbetween " << cut.lowerbound << " and " << cut.upperbound << endl;
			fEvtSelCuts.push_back(cut);
		}
		// else if();
		// Maybe implement here additional, more complicated cuts, to be read from the same file
	}
	cout << " ----------- " << endl;
}

bool TreeReader::IsGoodObj(int i, vector<Cut> *cutVec){
// This loops on the vector of cuts as filled in ReadObjCuts() and 
// determines if an object is good
	double *val;
	for(vector<Cut>::const_iterator it = cutVec->begin(); it != cutVec->end(); it++){
		val = (double*)it->branch->GetAddress();
		// cout << " Object " << i << " has " << it->branch->GetName() << " " << val[i] << endl;
		if(val[i] < it->lowerbound) return false;
		if(val[i] > it->upperbound) return false;
	}
	return true;
}

bool TreeReader::IsGoodEvt(vector<Cut> *cutVec){
// This loops on the vector of cuts as filled in ReadEvtSel()
// and determines if an event is selected
	double *val;
	for(vector<Cut>::const_iterator it = cutVec->begin(); it != cutVec->end(); it++){
		val = (double*)it->branch->GetAddress();
		// cout << " Event has " << it->branch->GetName() << " with " << *val << endl;
		if(*val < it->lowerbound) return false;
		if(*val > it->upperbound) return false;
	}
	return true;
}

///////////////////////////////////////////////////////////////////////////////////////////////
// Utilities //////////////////////////////////////////////////////////////////////////////////
void TreeReader::printPNG(TCanvas *cin, TString name, TString dir){
/*		-	Prints a ROOT TCanvas Object to a .png file
			name is the bare output filename, e.g. "fit_4_8",
			dir is the output directory (inside the overall output dir.)       */
	// Create sub directories if needed
	if(!dir.EndsWith("/")) dir += "/";
	char cmd[100];
	sprintf(cmd,"mkdir -p %s", dir.Data());
	system(cmd);

	dir += name;
	dir += ".png";
	cin->Print(dir,"png");
}

void TreeReader::printEPS(TCanvas *cin, TString name, TString dir){
/*		-	Prints a ROOT TCanvas Object to a .eps file
			name is the bare output filename, e.g. "fit_4_8",
			dir is the output directory (inside the overall output dir.)       */
	// Create sub directories if needed
	if(!dir.EndsWith("/")) dir += "/";
	char cmd[100];
	sprintf(cmd,"mkdir -p %seps/", dir.Data());
	system(cmd);

	dir += "eps/";
	dir += name;
	dir += ".eps";
	cin->SaveAs(dir);
}

double TreeReader::getEta(double x, double y, double z){
	if(fabs(z) <1.0e-5) return 0;
	double theta = atan(sqrt(x*x+y*y)/z);
	if(theta < 0.) theta = theta + 3.141592654;
	return -log(tan(0.5*theta));
}
