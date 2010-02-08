#include "TreeReader.hh"
using namespace std;

TreeReader::TreeReader(TTree *tree, int flag) : TreeClassBase(tree){
	fClean = true;
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
	if(fClean){
		ReadCleaningParameters();
		StatInit();
	}
	if(fDiLep)    InitDiLepTree();
	if(fMPHist)   BookMPHistos();
	if(fSignHist) BookSignHists();

	fCleanTreeFile = new TFile(fOutputDir + "CleanTree.root", "RECREATE");
	fCleanTreeFile->mkdir("analyze", "analyze");
	fCleanTreeFile->cd("analyze");
	fCleanTree = fChain->CloneTree(0);
	fCleanTree->CopyAddresses(fChain);
}

// Method called after finishing the event loop
void TreeReader::EndJob(){
	if(fDiLep)    WriteDiLepTree();
	if(fMPHist)   PrintMPOutput();
	if(fSignHist) WriteSignHists();
	if(fClean){
		StatPrint();
		StatHistos();
	}
	
	fCleanTreeFile->cd("analyze");
	fCleanTree->Write();
	fCleanTreeFile->Close();
	fHstatFile->cd();
	fHstatHistos->Write();
	fHstatFile->Close();
}

// Method for looping over the tree
void TreeReader::Loop(){
	Long64_t nbytes = 0, nb = 0;
	Long64_t nentries = fChain->GetEntries();
	cout << " total events in ntuples = " << fChain->GetEntries() << endl;
	// nentries = 10;

	for( Long64_t jentry = 0; jentry < nentries; jentry++ ){
		Long64_t ientry = LoadTree(jentry);
		if( ientry < 0 ) break;
		nb = fChain->GetEntry(jentry);
		nbytes += nb;
		if( jentry%200 == 0 ) cout << ">>> Processing event # " << jentry << endl;
		// if( jentry%10000 == 0 ) cout << ">>> Processing event # " << jentry << endl;

		// Event Selection
		if(!IsGoodEvt(&fEvtSelCuts)) continue;

		// Put here any method that needs to be called once every event
		// Do the object cleaning on the NTuples

		InitCleaning();
		TagCleanObjects();
		if( fClean ){
			DecideIso();
			DoCleanObjects();

			StatFill();

			// Write out a cleaned and skimmed tree
			int iclean = 0;
			for( int i = 0; i < NMus; i++ ){
				if( MuGood[i] != 0 ) continue;
				PutMuon(iclean, i);
				iclean++;
			}
			NMus = iclean;
			
			iclean = 0;
			for( int i = 0; i < NEles; i++ ){
				if( ElGood[i] != 0 ) continue;
				PutElectron(iclean, i);
				iclean++;
			}
			NEles = iclean;
			
			iclean = 0;
			for( int i = 0; i < NPhotons; i++ ){
				if( PhoGood[i] != 0 ) continue;
				PutPhoton(iclean, i);
				iclean++;
			}
			NPhotons = iclean;
			
			iclean = 0;
			for( int i = 0; i < NJets; i++ ){
				if( JGood[i] != 0 ) continue;
				PutJet(iclean, i);
				iclean++;
			}
			NJets = iclean;
			
			// At this point the arrays only contain clean objects
			// The non-clean object are overwritten
			fCleanTree->Fill();
		}

		if( fDiLep )    FillDiLepTree();
		if( fMPHist )   FillMPHistos();
		if( fSignHist ) for(size_t i = 0; i < 5; ++i) FillSignHists(i);
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
	for(int i=0; i < NJets; ++i){
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
		fTMu1iso      = MuRelIso03[lep1index];
		fTMu2iso      = MuRelIso03[lep2index];
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
		fTMu1iso      = MuRelIso03[qualMuInd[0]];
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
	TBranch *b;
	TLeaf *l;
	double *fval;
	int *ival;
	for(vector<Cut>::const_iterator it = cutVec->begin(); it != cutVec->end(); it++){
		b = it->branch;
		l = b->GetLeaf(b->GetName());
		if(!strcmp(l->GetTypeName(), "Double_t")){
			fval = (Double_t*)l->GetValuePointer();
			// cout << " Object " << i << " has " << it->branch->GetName() << " " << fval[i] << endl;
			if(fval[i] < it->lowerbound) return false;
			if(fval[i] > it->upperbound) return false;
		}
		if(!strcmp(l->GetTypeName(), "Int_t")){
			ival = (Int_t*)l->GetValuePointer();
			// cout << " Object " << i << " has " << it->branch->GetName() << " " << ival[i] << endl;
			if(ival[i] < it->lowerbound) return false;
			if(ival[i] > it->upperbound) return false;
		}
	}
	return true;
}

bool TreeReader::IsGoodEvt(vector<Cut> *cutVec){
// This loops on the vector of cuts as filled in ReadEvtSel()
// and determines if an event is selected
	TBranch *b;
	TLeaf *l;
	double *fval;
	int *ival;
		// double *val;
	for(vector<Cut>::const_iterator it = cutVec->begin(); it != cutVec->end(); it++){
		b = it->branch;
		l = b->GetLeaf(b->GetName());
		if(!strcmp(l->GetTypeName(), "Double_t")){
			fval = (double*)l->GetValuePointer();
			// cout << " Object " << i << " has " << it->branch->GetName() << " " << fval[i] << endl;
			if(*fval < it->lowerbound) return false;
			if(*fval > it->upperbound) return false;
		}
		if(!strcmp(l->GetTypeName(), "Int_t")){
			ival = (int*)l->GetValuePointer();
			// cout << " Object " << i << " has " << it->branch->GetName() << " " << ival[i] << endl;
			if(*ival < it->lowerbound) return false;
			if(*ival > it->upperbound) return false;
		}
	}
	return true;
}

///////////////////////////////////////////////////////////////////////////////////////////////
// Cleaning Tools /////////////////////////////////////////////////////////////////////////////
void TreeReader::ReadCleaningParameters(const char* filename){
	ifstream IN(filename);
	char buffer[200];
	char ParName[100];
	float ParValue;

	bool verbose(true);
	bool ok(false);

	while( IN.getline(buffer, 200, '\n') ){
		ok = false;
		if (buffer[0] == '#') {continue;} // Skip lines commented with '#'
		sscanf(buffer, "%s %f", ParName, &ParValue);

		// -- Primary vertex:
		if( !strcmp(ParName, "chisqVxmax") ){
			fClean_chisqVxmax = float(ParValue); ok = true;
		}
		if( !strcmp(ParName, "dRVxmax") ){
			fClean_dRVxmax = float(ParValue); ok = true;
		}
		if( !strcmp(ParName, "dzVxmax") ){
			fClean_dzVxmax = float(ParValue); ok = true;
		}
		if( !strcmp(ParName, "sumPtTkfromVxmin") ){
			fClean_sumPtTkfromVxmin = float(ParValue); ok = true;
		}
		if( !strcmp(ParName, "distVxmax") ){
			fClean_distVxmax = float(ParValue); ok = true;
		}

		// -- Muons:
		if( !strcmp(ParName, "MuonDPbyPmax") ){
			fClean_MuonDPbyPmax = float(ParValue); ok = true;
		}
		if( !strcmp(ParName, "MuonChi2max") ){
			fClean_MuonChi2max = float(ParValue); ok = true;
		}
		if( !strcmp(ParName, "MuonNHitsmin") ){
			fClean_MuonNHitsmin = float(ParValue); ok = true;
		}
		if( !strcmp(ParName, "dRSSmuonmax") ){
			fClean_dRSSmuonmax = float(ParValue); ok = true;
		}

		// -- Electrons:
		if( !strcmp(ParName, "ElecHoverEBarmax") ){
			fClean_ElecHoverEBarmax = float(ParValue); ok = true;
		}
		if( !strcmp(ParName, "ElecHoverEEndmax") ){
			fClean_ElecHoverEEndmax = float(ParValue); ok = true;
		}
		if( !strcmp(ParName, "ElecSigmaEtaEtaBarmax") ){
			fClean_ElecSigmaEtaEtaBarmax = float(ParValue); ok = true;
		}
		if( !strcmp(ParName, "ElecSigmaEtaEtaEndmax") ){
			fClean_ElecSigmaEtaEtaEndmax = float(ParValue); ok = true;
		}
		if( !strcmp(ParName, "ElecEoverPInBarmin") ){
			fClean_ElecEoverPInBarmin = float(ParValue); ok = true;
		}
		if( !strcmp(ParName, "ElecEoverPInEndmin") ){
			fClean_ElecEoverPInEndmin = float(ParValue); ok = true;
		}
		if( !strcmp(ParName, "ElecDeltaEtaInBarmax") ){
			fClean_ElecDeltaEtaInBarmax = float(ParValue); ok = true;
		}
		if( !strcmp(ParName, "ElecDeltaEtaInEndmax") ){
			fClean_ElecDeltaEtaInEndmax = float(ParValue); ok = true;
		}
		if( !strcmp(ParName, "ElecDeltaPhiInBarmax") ){
			fClean_ElecDeltaPhiInBarmax = float(ParValue); ok = true;
		}
		if( !strcmp(ParName, "ElecDeltaPhiInEndmax") ){
			fClean_ElecDeltaPhiInEndmax = float(ParValue); ok = true;
		}
		if( !strcmp(ParName, "ElecDeltaPhiOutBarmax") ){
			fClean_ElecDeltaPhiOutBarmax = float(ParValue); ok = true;
		}
		if( !strcmp(ParName, "ElecDeltaPhiOutEndmax") ){
			fClean_ElecDeltaPhiOutEndmax = float(ParValue); ok = true;
		}
		if( !strcmp(ParName, "dRSSelecmax") ){
			fClean_dRSSelecmax = float(ParValue); ok = true;
		}

		// -- Photons:
		if( !strcmp(ParName, "PhotHoverEBarmax") ){
			fClean_PhotHoverEBarmax = float(ParValue); ok = true;
		}
		if( !strcmp(ParName, "PhotHoverEEndmax") ){
			fClean_PhotHoverEEndmax = float(ParValue); ok = true;
		}
		if( !strcmp(ParName, "PhotSigmaEtaEtaBarmax") ){
			fClean_PhotSigmaEtaEtaBarmax = float(ParValue); ok = true;
		}
		if( !strcmp(ParName, "PhotSigmaEtaEtaEndmax") ){
			fClean_PhotSigmaEtaEtaEndmax = float(ParValue); ok = true;
		}

		// -- Jets:
		if( !strcmp(ParName, "FracEmmaxJet") ){
			fClean_FracEmmaxJet = float(ParValue); ok = true;
		}
		if( !strcmp(ParName, "FracEmminJet") ){
			fClean_FracEmminJet = float(ParValue); ok = true;
		}
		if( !strcmp(ParName, "FracChminJet") ){
			fClean_FracChminJet = float(ParValue); ok = true;
		}
		if( !strcmp(ParName, "deltaRElecJetmax") ){
			fClean_deltaRElecJetmax = float(ParValue); ok = true;
		}
		if( !strcmp(ParName, "elecbyJetEratio") ){
			fClean_elecbyJetEratio = float(ParValue); ok = true;
		}

		// -- Isolation:
		if( !strcmp(ParName, "MuonIsomax") ){
			fClean_MuonIsomax = float(ParValue); ok = true;
		}
		if( !strcmp(ParName, "ElecIsomax") ){
			fClean_ElecIsomax = float(ParValue); ok = true;
		}
		if( !strcmp(ParName, "PhotIsomax") ){
			fClean_PhotIsomax = float(ParValue); ok = true;
		}

		// -- Event Cleaning:
		if( !strcmp(ParName, "FracChmin") ){
			fClean_FracChmin = float(ParValue); ok = true;
		}
		if( !strcmp(ParName, "FracEmmin") ){
			fClean_FracEmmin = float(ParValue); ok = true;
		}

		// -- MET:
		if( !strcmp(ParName, "METmin") ){
			fClean_METmin = float(ParValue); ok = true;
		}
		if( !strcmp(ParName, "dPhiJetMETmin") ){
			fClean_dPhiJetMETmin = float(ParValue); ok = true;
		}
		if( !strcmp(ParName, "dR12min") ){
			fClean_dR12min = float(ParValue); ok = true;
		}
		if( !strcmp(ParName, "dR21min") ){
			fClean_dR21min = float(ParValue); ok = true;
		}

		if(!ok) cout << "%% TreeReader::ReadCleaningParameters ==> ERROR: Unknown variable " << ParName << endl;
	}
}

void TreeReader::TagCleanObjects(void){
// Steering for tagging clean/bad objects

	// Primary vertex
	PrimVtxGood = CleanPrimaryVertex();

	// Muons
	for( int ichk = 0; ichk < NMus; ++ichk ){
		MuGood[ichk] = 10*IsFromPrimaryVx(1, ichk);
		MuGood[ichk] += CleanMuon(ichk);;
	}

	// Electrons
	for( int ichk = 0; ichk < NEles; ++ichk ){
		ElGood[ichk] = 10*IsFromPrimaryVx(2, ichk);
		ElGood[ichk] += CleanElectron(ichk);
	}

	// Photons
	for( int ichk = 0; ichk < NPhotons; ++ichk ){
		PhoGood[ichk] = CleanPhoton(ichk);
	}

	// Jets
	for( int ichk = 0; ichk < NJets; ++ichk ){
		JGood[ichk] = 10*IsFromPrimaryVx(4, ichk);
		JGood[ichk] += CleanJet(ichk);
	}

	// Duplication (only after cleanness has been checked)
	for( int ichk = 0; ichk < NMus; ++ichk )  if( DuplicateMuon(ichk)     ) MuGood[ichk] += 100;
	for( int ichk = 0; ichk < NEles; ++ichk ) if( DuplicateElectron(ichk) ) ElGood[ichk] += 100;
	for( int ichk = 0; ichk < NJets; ++ichk ) if( ElectronJet(ichk)       ) JGood[ichk]  += 100;

	// Event and MET can only be done later (after bad objects are removed)
	GoodEvent = 0;
	return;
}

int TreeReader::CleanPrimaryVertex(void){
// Verifies the primary vertex quality
// returns iBad = 1 for no charged tracks
//              = 2 for bad normalized chi squared AND 1
//              = 3 for incompatible with beamspot AND 1 AND 2
//              = 4 for insufficient track pT AND 1 AND 2 AND 3
	int iBad = 0;
// Check that there are tracks at the Primary Vertex
	if( PrimVtxNTracks <= 0 ){
		iBad = 1;
		return iBad;
	}

// Check the chi2/ndof
	if( PrimVtxNChi2 > fClean_chisqVxmax || PrimVtxNChi2 < 0. ){
		iBad = 2;
		return iBad;
	}

// Check compatibility of vertex with beam spot
	double xVx = PrimVtxx - Beamspotx;
	double yVx = PrimVtxy - Beamspoty;
	double zVx = PrimVtxz - Beamspotz;
	double rVx = sqrt(xVx*xVx + yVx*yVx);
	if (rVx > fClean_dRVxmax || fabs(zVx) > fClean_dzVxmax) {
		iBad = 3;
		return iBad;
	}

// Check that there is sufficient Et in the tracks
	if( PrimVtxPtSum < fClean_sumPtTkfromVxmin ){
		iBad = 4;
		return iBad;
	}
	return iBad;
}

int TreeReader::IsFromPrimaryVx(int ipart, int ichk){
// Checks whether the object is compatible with the primary vertex
// ipart = 1 for muon
//       = 2 for electron
//       = 3 for photon
//       = 4 for jet
// returns iBad = 1 for incompatible with primary vertex
	int iBad = 0;
// take error from vertex and from track extrapolation into account
	double drVxsq = PrimVtxxE*PrimVtxxE + PrimVtxyE*PrimVtxyE;
	double d0, dd0, dz, ddz;
	if( ipart <= 0 || ichk < 0 ){
		return iBad;
	} else if( ipart == 1 ){ // Muons
		d0  = MuD0PV[ichk];
		dd0 = sqrt(MuD0E[ichk]*MuD0E[ichk] + drVxsq);
		dz  = MuDzPV[ichk];
		ddz = sqrt(MuDzE[ichk]*MuDzE[ichk] + PrimVtxzE*PrimVtxzE);
	} else if( ipart == 2 ){ // Electrons
		d0  = ElD0PV[ichk];
		dd0 = sqrt(ElD0E[ichk]*ElD0E[ichk] + drVxsq);
		dz  = ElDzPV[ichk];
		ddz = sqrt(ElDzE[ichk]*ElDzE[ichk] + PrimVtxzE*PrimVtxzE);
	} else if( ipart == 3 ){ // Photons
		return 1; // Photons not implemented yet
	} else if( ipart == 4 ){ // Jets
		d0  = sqrt((JVtxx[ichk]-PrimVtxx)*(JVtxx[ichk]-PrimVtxx)
			+ (JVtxy[ichk]-PrimVtxy)*(JVtxy[ichk]-PrimVtxy) );
		dd0 = sqrt(JVtxExx[ichk] + JVtxEyy[ichk] + drVxsq);
		dz  = JVtxz[ichk] - PrimVtxz;
		ddz = sqrt(JVtxEzz[ichk] + PrimVtxzE*PrimVtxzE);
	}

// test that the distance is not too large
	if( fabs(d0) > fClean_distVxmax * dd0 || fabs(dz) > fClean_distVxmax * ddz ){
		iBad = 1;
	}
	return iBad;
}

int TreeReader::CleanMuon(int ichk){
// Verifies the muon identification quality
// returns iBad = 0 for good muons
//              = 1 for bad Delta pT / pT
//              = 2 for bad normalized chi squared
//              = 3 for too few valid hits in tracker
	int iBad = 0;
	if( ichk < 0 ) return -1;

	if( MuPtE[ichk] >= fClean_MuonDPbyPmax * MuPt[ichk] ) iBad = 1; // Maximum Delta p / p
	else if( MuNChi2[ichk] > fClean_MuonChi2max )         iBad = 2; // Maximum Chisquared
	else if( MuNTkHits[ichk] < fClean_MuonNHitsmin )      iBad = 3; // Minimum number of valid hits
	return iBad;
}

bool TreeReader::DuplicateMuon(int ichk){
// Checks for duplicate muons
	// fClean_dRSSmuonmax = 0.1;

	if( ichk < 0 ) return false;

	for( int j = 0; j < NMus; ++j ){
		if( j == ichk ) continue;
		if( MuCharge[ichk] != MuCharge[j] ) continue;

		double deltaR = GetDeltaR(MuEta[ichk], MuEta[j], MuPhi[ichk], MuPhi[j]);
		if( deltaR > fClean_dRSSmuonmax ) continue;

		// Both are bad or both are good -> compare them
		if( (MuGood[ichk] == 0 && MuGood[j] == 0) || (MuGood[ichk] != 0 && MuGood[j] != 0) ){
			if( MuPtE[ichk]/MuPt[ichk] >= MuPtE[j]/MuPt[j] ) return true;
		}
		// One is good, one is bad -> take the good one
		else if( MuGood[j] == 0 && MuGood[ichk] != 0 ) return true;
	}
	return false;
}

int TreeReader::CleanElectron(int ichk){
// Verifies the electron identification quality
// returns iBad = 0 for good muons
//              = 1 for bad H/E
//              = 2 for bad shower shape
//              = 3 for bad matching of Ecal and track
	bool useHoverE      = true;
	bool useSigmaEtaEta = true;
	bool useEoverPIn    = true;
	bool useDeltaEtaIn  = true;
	bool useDeltaPhiIn  = true;
	bool useDeltaPhiOut = true;
	if( ichk < 0 ) return -1;

	double hOverE      = ElHcalOverEcal[ichk];
	double sigmaee     = ElSigmaIetaIeta[ichk];
	double eOverPin    = ElESuperClusterOverP[ichk];
	double deltaEtaIn  = ElDeltaEtaSuperClusterAtVtx[ichk];
	double deltaPhiIn  = ElDeltaPhiSuperClusterAtVtx[ichk];
	double deltaPhiOut = ElDeltaPhiSeedClusterAtCalo[ichk];

	// Distinguish between barrel and endcap electrons
	if( fabs(ElEta[ichk]) < 1.479 ){ // Barrel
		if(useHoverE)      if( hOverE            > fClean_ElecHoverEBarmax      ) return 1;
		if(useSigmaEtaEta) if( sigmaee           > fClean_ElecSigmaEtaEtaBarmax ) return 2;
		if(useEoverPIn)    if( eOverPin          < fClean_ElecEoverPInBarmin    ) return 3;
		if(useDeltaEtaIn)  if( fabs(deltaEtaIn)  > fClean_ElecDeltaEtaInBarmax  ) return 3;
		if(useDeltaPhiIn)  if( fabs(deltaPhiIn)  > fClean_ElecDeltaPhiInBarmax  ) return 3;
		if(useDeltaPhiOut) if( fabs(deltaPhiOut) > fClean_ElecDeltaPhiOutBarmax ) return 3;
	}
	else{ // EndCap
		if(useHoverE)      if( hOverE            > fClean_ElecHoverEEndmax      ) return 1;
		if(useSigmaEtaEta) if( sigmaee           > fClean_ElecSigmaEtaEtaEndmax ) return 2;
		if(useEoverPIn)    if( eOverPin          < fClean_ElecEoverPInEndmin    ) return 3;
		if(useDeltaEtaIn)  if( fabs(deltaEtaIn)  > fClean_ElecDeltaEtaInEndmax  ) return 3;
		if(useDeltaPhiIn)  if( fabs(deltaPhiIn)  > fClean_ElecDeltaPhiInEndmax  ) return 3;
		if(useDeltaPhiOut) if( fabs(deltaPhiOut) > fClean_ElecDeltaPhiOutEndmax ) return 3;
	}
	return 0;
}

bool TreeReader::DuplicateElectron(int ichk){
// Checks for duplicate electrons
	if( ichk < 0 ) return false;

	int j = ElDuplicateEl[ichk];
	if( j < 0 ) return false;

	if( GetDeltaR(ElEta[ichk], ElEta[j], ElPhi[ichk], ElPhi[j]) < fClean_dRSSelecmax ){
		// Both are bad or both are good -> compare them
		if( (ElGood[ichk] == 0 && ElGood[j] == 0) || (ElGood[ichk] != 0 && ElGood[j] != 0) ){
			double elecEoP = ElESuperClusterOverP[ichk];
			double newEoP  = ElESuperClusterOverP[j];
			if( fabs(ElESuperClusterOverP[ichk]-1.) > fabs(ElESuperClusterOverP[j]-1.) ) return true;
		}
		// One is good, one is bad -> take good one
		else if( ElGood[j] == 0 && ElGood[ichk] != 0 ) return true;
	}
	return false;
}

int TreeReader::CleanPhoton(int ichk){
// Verifies the photon identification quality
// returns iBad = 1 for bad H/E
// Still to be completed ****
	if( ichk < 0 ) return 0;
	bool useHoverE      = true;
	bool useSigmaEtaEta = true;

	// Distinguish between barrel and endcap photons
	if( fabs(PhoEta[ichk]) < 1.479 ){ // Barrel
		if(useHoverE)      if( PhoHoverE[ichk]        > fClean_PhotHoverEBarmax      ) return 1;
		if(useSigmaEtaEta) if( PhoSigmaIetaIeta[ichk] > fClean_PhotSigmaEtaEtaBarmax ) return 2;
	}
	else{ // EndCap
		if(useHoverE)      if( PhoHoverE[ichk]        > fClean_PhotHoverEEndmax      ) return 1;
		if(useSigmaEtaEta) if( PhoSigmaIetaIeta[ichk] > fClean_PhotSigmaEtaEtaEndmax ) return 2;
	}

	return 0;
}

int TreeReader::CleanJet(int ichk){
// Verifies the jet reconstruction quality
// and flags jets made from electrons
// (electrons should be filled before cleaning the jets)
// returns iBad = 0 for good jets
//              = 1 Et < Pt
//              = 2 for too large EM fraction
//              = 3 for too small EM fraction
//              = 4 for bad Trk pT fraction
	int iBad = 0;
	if( ichk < 0 ) return 0;

// veto jets with E<p
	if( JEt[ichk] - JPt[ichk] < -0.0001 ) return 1;

// check EM and track pT fractions
	if( JEMfrac[ichk] > fClean_FracEmmaxJet ) return 2;
	if( JEMfrac[ichk] < fClean_FracEmminJet ) return 3;
	if( JChfrac[ichk] < fClean_FracChminJet && fabs(JEta[ichk]) < 2.1 ) return 4;
	return 0;
}

bool TreeReader::ElectronJet(int ichk){
// checks for jets made from electrons
// ichk = index of the jet
	if( ichk < 0 ) return false;

	bool isDuplicate = false;

// veto jets made of electrons
	for( int j = 0; j < NEles; ++j ){
		if( ElIsInJet[j] < 0 ) continue;
		if( ElIsInJet[j] != ichk ) continue;
		if( GetDeltaR(JEta[ichk], ElEta[j], JPhi[ichk], ElPhi[j]) > fClean_deltaRElecJetmax ) continue;
		if( ElSharedEnergy[j] > fClean_elecbyJetEratio * JE[ichk] ){
			isDuplicate = true;
			break;
		}
	}
	return isDuplicate;
}

void TreeReader::DecideIso(void){
// decide whether objects are isolated or not
	// Muons
	for( int ichk = 0; ichk < NMus; ++ichk ){
		if( MuRelIso03[ichk] < fClean_MuonIsomax ) MuIsIso[ichk] = 1;
		else MuIsIso[ichk] = 0;
	}

	// Electrons
	for( int ichk = 0; ichk < NEles; ++ichk ){
		if( ElIso[ichk] < fClean_ElecIsomax ) ElIsIso[ichk] = 1;
		else ElIsIso[ichk] = 0;
	}

	// Photons
	for ( int ichk = 0; ichk < NPhotons; ++ichk){
		if ( PhoIso03[ichk] < fClean_PhotIsomax ) PhoIsIso[ichk] = 1;
		else PhoIsIso[ichk] = 0;
	}

	return;
}

void TreeReader::InitCleaning(){
// Initializes the object cleaning
// Has to be called once each event!
	if(fClean){
		fNMuClean = 0;
		fNElClean = 0;
		fNJClean  = 0;
	} else {
		fNMuClean = NMus;
		fNElClean = NEles;
		fNJClean  = NJets;
	}
	return;
}

void TreeReader::DoCleanObjects(void){
// Steering for removal of bad objects
// Should be initialized by calling InitCleaning
// Should only be called after deciding on isolation

	// Check the primary vertex
	if( PrimVtxGood != 0 ){
		fNMuClean = 0;
		fNElClean = 0;
		fNJClean = 0;
	}

	// Muons
	fNMuClean = 0;
	for( int ichk = 0; ichk < NMus; ++ichk ){
		if( MuGood[ichk] != 0 || MuIsIso[ichk] != 1 ) continue;
		fNMuClean++;
	}

	// Electrons
	fNElClean = 0;
	for( int ichk = 0; ichk < NEles; ++ichk ){
		if( ElGood[ichk] != 0 || ElIsIso[ichk] != 1 ) continue;
		fNElClean++;
	}

	// Photons
	fNPhClean = 0;
	for( int ichk = 0; ichk < NPhotons; ++ichk ){
		if( PhoGood[ichk] != 0 || PhoIsIso[ichk] != 1 ) continue;
		fNPhClean++;
	}

	// Clean the jets
	fNJClean = 0;
	for( int ichk = 0; ichk < NJets; ++ichk ){
		if( JGood[ichk] != 0 ) continue;
		fNJClean++;
		for( int iel = 0; iel < NEles; ++iel ){
			if( ElIsInJet[iel] != ichk ) continue;
			if( ElIsIso[iel] ) SubtrFromJet(2, iel, ichk);
			else AddToJet(2, iel, ichk);
		}
	}

	// Check the event
	// (Should do CleanEvent only on remaining objects)
	GoodEvent = CleanEvent();
	if( fNMuClean + fNElClean + fNJClean <= 0 || GoodEvent%10 != 0 ){
		fNMuClean = 0;
		fNElClean = 0;
		fNJClean = 0;
		return;
	}

// Check the MET
// (Should do CleanMET only on remaining objects)
// Could we choose one of the MET types?
// Or do we let the user choose???
	GoodEvent += 10   * CleanMET(TCMET, TCMETphi);
	GoodEvent += 100  * CleanMET(MuJESCorrMET, MuJESCorrMETphi);
	GoodEvent += 1000 * CleanMET(PFMET, PFMETphi);
	if( GoodEvent/10 == 111 ){
		fNMuClean = 0;
		fNElClean = 0;
		fNJClean = 0;
	}
	return;
}

void TreeReader::AddToJet(int ipart, int ichk, int iJet){
// adds an object (e or mu) to its nearest jet
// ipart = 1 for muon
//       = 2 for electron
//       = 3 for photon
//       = 4 for jet
// ichk = index of the object to be added
// iJet = index of the jet

	if( ichk >= 0 && iJet >= 0 ){
		double pxadd, pyadd, pzadd, eadd, ptadd, emadd;
		if( ipart == 1 ){
			pxadd = MuPx[ichk];
			pyadd = MuPy[ichk];
			pzadd = MuPz[ichk];
			eadd  = MuE[ichk];
			ptadd = MuPt[ichk];
			emadd = 0.;
		} else if( ipart == 2 ){
			pxadd = ElPx[ichk] - ElSharedPx[ichk];
			pyadd = ElPy[ichk] - ElSharedPy[ichk];
			pzadd = ElPz[ichk] - ElSharedPz[ichk];
			eadd  = ElE[ichk] - ElSharedEnergy[ichk];
			ptadd = ElPt[ichk];
			double sharedpt = sqrt(ElSharedPx[ichk]*ElSharedPx[ichk]
				+ ElSharedPy[ichk]*ElSharedPy[ichk]);
			emadd = ElPt[ichk] - sharedpt;
		} else if( ipart == 3 ){
			// cout << " *** Problem: adding photons is not foreseen" << endl;
			return;
		} else if( ipart == 4 ){
			pxadd = JPx[ichk];
			pyadd = JPy[ichk];
			pzadd = JPz[ichk];
			eadd  = JE[ichk];
			ptadd = JChfrac[ichk] * JPt[ichk];
			emadd = JEMfrac[ichk] * JPt[ichk];
		}

		JPx[iJet] += pxadd;
		JPy[iJet] += pyadd;
		JPz[iJet] += pzadd;
		JE[iJet] += eadd;
	// ??? or do we want to keep the jets massless?
		JPt[iJet] = sqrt(JPx[iJet]*JPx[iJet] + JPy[iJet]*JPy[iJet]);
		JEt[iJet] = JPt[iJet];
		if( fabs(JPz[iJet]) < 1.0e-5 ) JEta[iJet] = 0.;
		else {
			double theta = atan(JPt[iJet]/JPz[iJet]);
			if( theta < 0. ) theta += 3.141592654;
			JEta[iJet] =  -log(tan(0.5*theta));
		}
		JPhi[iJet] =  atan2(JPy[iJet],JPx[iJet]);
		JChfrac[iJet] += ptadd / JPt[iJet];
		JEMfrac[iJet] += emadd / JPt[iJet];
	}
	return;
}

void TreeReader::SubtrFromJet(int ipart, int ichk, int iJet){
// subtracts an object from its nearest jet
// ipart = 1 for muon
//       = 2 for electron
//       = 3 for photon
//       = 4 for jet
// ichk = index of the object to be subtracted
// iJet = index of the jet

	if( ichk >= 0 && iJet >= 0 ){
		double pxadd, pyadd, pzadd, eadd, ptadd, emadd;
		if( ipart == 1 ){
			pxadd = MuPx[ichk];
			pyadd = MuPy[ichk];
			pzadd = MuPz[ichk];
			eadd  = MuE[ichk];
			ptadd = MuPt[ichk];
			emadd = 0.;
		} else if( ipart == 2 ){
			pxadd = ElSharedPx[ichk];
			pyadd = ElSharedPy[ichk];
			pzadd = ElSharedPz[ichk];
			eadd  = ElSharedEnergy[ichk];
			ptadd = ElPt[ichk];
			double emshpt = sqrt(ElSharedPx[ichk]*ElSharedPx[ichk]
				+ ElSharedPy[ichk]*ElSharedPy[ichk]);
			emadd = ElPt[ichk] - emshpt;
		} else if( ipart == 3 ){
			// cout << " *** Problem: adding photons is not foreseen" << endl;
			return;
		} else if( ipart == 4 ){
			pxadd = JPx[ichk];
			pyadd = JPy[ichk];
			pzadd = JPz[ichk];
			eadd  = JE[ichk];
			ptadd = JChfrac[ichk] * JPt[ichk];
			emadd = JEMfrac[ichk] * JPt[ichk];
		}

		JPx[iJet] -= pxadd;
		JPy[iJet] -= pyadd;
		JPz[iJet] -= pzadd;
		JE[iJet] -= eadd;
	// ??? or do we want to keep the jets massless?
		double jmom = sqrt (JPx[iJet]*JPx[iJet] + JPy[iJet]*JPy[iJet] + JPz[iJet]*JPz[iJet]);
		if( JE[iJet] < jmom ){
			double egy = JE[iJet];
			if (egy <= 0.){ egy = 0.001;}
			double scale = egy / jmom;
			JPx[iJet] *= scale;
			JPy[iJet] *= scale;
			JPz[iJet] *= scale;
		}
		JPt[iJet] = sqrt(JPx[iJet]*JPx[iJet] + JPy[iJet]*JPy[iJet]);
		JEt[iJet] = JPt[iJet];
		if( fabs(JPz[iJet]) <1.0e-5 ) JEta[iJet] = 0.;
		else {
			double theta = atan(JPt[iJet]/JPz[iJet]);
			if( theta < 0. ) theta += 3.141592654;
			JEta[iJet] = -log(tan(0.5*theta));
		}
		JPhi[iJet] = atan2(JPy[iJet],JPx[iJet]);
		JChfrac[iJet] = ptadd / JPt[iJet];
		JEMfrac[iJet] = emadd / JPt[iJet];
	}
	return;
}

int TreeReader::CleanEvent(void){
// To veto events from beam halo, cosmics or noise
// will also need the primary vertex
// tests on Fem and Ftrk (careful for mu-mu and e-e)
// needs to be modified when photons become available
// returns iBad = 0 event is good
//              = 1 event is empty
//              = 2 for too small EM fraction
//              = 3 for bad Trk pT fraction
// test that data still exist
	if( fNMuClean + fNElClean + fNJClean <= 0 ) return 1;

// test total Fem and Ftrk in event
	int nPhot = 0;
	int nChObj = 0;
	double pt_track = 0.;
	double et_em = 0.;
	double et_had = 0.;
	for( int i = 0; i < NMus; ++i ){
		if(MuGood[i] != 0) continue;
		pt_track += MuPt[i];
		nChObj++;
	}
	for( int i = 0; i < NEles; ++i ){
		if(ElGood[i] != 0) continue;
		pt_track += ElPt[i];
		et_em += ElEt[i];
		nChObj++;
	}
	for( int i = 0; i < NJets; ++i ){
		if(JGood[i] != 0) continue;
		pt_track += JChfrac[i] * JPt[i];
		et_em    += JEMfrac[i] * JEt[i];
		et_had   += (1.-JEMfrac[i]) * JEt[i];
		if( JChfrac[i] > 0. ) nChObj++;
	}

	double fracCh = 0.;
	double fracEm = 0.;
	if( et_em + et_had <= 0. ){
		if( fNMuClean < 1 ) return 1; // bad prim vtx will trigger this...
		fracCh = 1.;
		fracEm = 1.;
	} else {
		fracCh = pt_track / (et_em + et_had);
		fracEm = et_em / (et_em + et_had);
	}
	if( fracEm < fClean_FracEmmin ) return 2;
	if( fracCh < fClean_FracChmin && (nPhot < 1 || nChObj > 0) ) return 3;
	return 0;
}

int TreeReader::CleanMET(double met, double metphi){
// The MET should not be aligned with any jet
// and should not be along one and opposite to the other of the 2 leading jets
// to veto QCD events with jet "mismeasurements"
// (see Jet + MET analysis in ptdr2)
// returns iBad = 0 for good MET
//              = 1 for MET aligned with a jet
//              = 2 for MET within Rij limit
	int iBad = 0;
	if( met == 0 ) return 0;
	// Care only if MET is large enough
	if (met < fClean_METmin) return 0;

	double etmax1 = 0.;
	double etmax2 = 0.;
	int imax1 = -1;
	int imax2 = -1;

	// Loop over all jets
	for( int i = 0; i < NJets; ++i ){
		// Reject if the MET is along the jet
		if(JGood[i] != 0) continue;
		double dPhi =  DeltaPhi(JPhi[i], metphi);
		if( dPhi < fClean_dPhiJetMETmin ) return 1;
		// Else, pick up the 2 leading jets
		if( JPt[i] > etmax1 ){
			etmax2 = etmax1;
			imax2  = imax1;
			etmax1 = JPt[i];
			imax1  = i;
		} else if( JPt[i] > etmax2 ){
			etmax2 = JPt[i];
			imax2  = i;
		}
	}

	// Check dR12 and dR21
	if( imax2 >= 0 ){
		double dPhi1 = DeltaPhi(JPhi[imax1], metphi );
		double dPhi2 = DeltaPhi(JPhi[imax2], metphi );
		double pi = 3.141592654;
		double r12 = sqrt(dPhi1*dPhi1 + (pi-dPhi2)*(pi-dPhi2) );
		double r21 = sqrt(dPhi2*dPhi2 + (pi-dPhi1)*(pi-dPhi1) );
		if( r12 < fClean_dR12min || r21 < fClean_dR21min ) return 2;
	}
	return 0;
}

int TreeReader::FindNearestJet(double eta, double phi){
// Looks for the nearest jet in deltaR to a given object (e or mu or photon)
// and returns its index
// returns -1 if no nearest jet

	int iJetMin = -1;

	double deltaRmin = 999.;
	for(int i = 0; i < NJets; ++i){
		double deltaR = GetDeltaR(eta, JEta[i], phi, JPhi[i]);
		if (deltaR < deltaRmin){
			deltaRmin = deltaR;
			iJetMin = i;
		}
	}
	return iJetMin;
}

double TreeReader::DeltaPhi(double v1, double v2){
// Computes the correctly normalized phi difference
// v1, v2 = phi of object 1 and 2
	const double pi    = 3.141592654;
	double twopi = 6.283185307;

	double diff = fabs(v2 - v1);
	double corr = twopi - diff;
	if (diff < pi){ return diff;} else { return corr;}
}

double TreeReader::GetDeltaR(double eta1, double eta2, double phi1, double phi2){
// Computes the DeltaR of two objects from their eta and phi values
	return sqrt( (eta1-eta2)*(eta1-eta2)
		+ DeltaPhi(phi1, phi2)*DeltaPhi(phi1, phi2) );
}

void TreeReader::StatInit(const char* filename){
// initializes the cleaning statistics
// to be called once at the beginning of the job
	fNumTotEvt               = 0;
	fNumTotEvtReject         = 0;
	fNumTotEvtEmpty          = 0;
	fNumTotEvtCleanEmpty     = 0;
	fNumTotEvtLtFem          = 0;
	fNumTotEvtLtFch          = 0;
	fNumTotEvtPfMETJet       = 0;
	fNumTotEvtPfMETRij       = 0;
	fNumTotEvtCaMETJet       = 0;
	fNumTotEvtCaMETRij       = 0;
	fNumTotEvtTcMETJet       = 0;
	fNumTotEvtTcMETRij       = 0;
	fNumTotEvtBadHardJet     = 0;

	fNumTotMuons             = 0;  
	fNumTotMuonGoodIso       = 0;  
	fNumTotMuonGoodNonIso    = 0;  
	fNumTotMuonBadIso        = 0;  
	fNumTotMuonBadNonIso     = 0;  
	fNumTotMuonDupl          = 0;
	fNumTotMuonNotPrimaryTrk = 0;
	fNumTotMuonNotClean      = 0;
	fNumTotMuonBadDpop       = 0;
	fNumTotMuonBadChi2       = 0;  
	fNumTotMuonBadNhit       = 0;  

	fNumTotElectrons         = 0;
	fNumTotElecGoodIso       = 0;
	fNumTotElecGoodNonIso    = 0;
	fNumTotElecBadIso        = 0;
	fNumTotElecBadNonIso     = 0;
	fNumTotElecDupl          = 0;
	fNumTotElecNotPrimaryTrk = 0;
	fNumTotElecNotClean      = 0;
	fNumTotElecBadHoE        = 0;
	fNumTotElecBadShsh       = 0;
	fNumTotElecBadTmat       = 0;

	fNumTotPhotons           = 0;
	fNumTotPhotGoodIso       = 0;
	fNumTotPhotGoodNonIso    = 0;
	fNumTotPhotBadIso        = 0;
	fNumTotPhotBadNonIso     = 0;
	fNumTotPhotDupl          = 0;
	fNumTotPhotNotClean      = 0;
	fNumTotPhotBadHoE        = 0;
	fNumTotPhotBadShsh       = 0;

	fNumTotJets              = 0;  
	fNumTotJetGood           = 0;  
	fNumTotJetBad            = 0;  
	fNumTotJetDuplElJet      = 0;
	fNumTotJetNotPrimaryTrk  = 0;
	fNumTotJetNotClean       = 0;
	fNumTotJetPgtE           = 0;
	fNumTotJetGtFem          = 0;
	fNumTotJetLtFem          = 0;
	fNumTotJetLtFch          = 0;
	fNumTotBJets             = 0;  

	// initialize cleaning statistics histogram
	fHstatFile = new TFile(fOutputDir + TString(filename), "RECREATE");
	fHstatHistos  = new TH1D("CleanStats", "cleaning statistics", 100, 0, 100);

	return;
}

void TreeReader::StatFill(void){
// accumulates the cleaning statistics
// to be called once per event, after cleaning is done

	fNumTotEvt++;
	int nObjects = 0;

	fNumTotMuons += NMus;
	for( int i=0; i < NMus; ++i ){
		if(      MuGood[i] == 0 && MuIsIso[i] == 1 ){ fNumTotMuonGoodIso++; nObjects++; }
		else if( MuGood[i] == 0 && MuIsIso[i] == 0 ) fNumTotMuonGoodNonIso++;
		else if( MuGood[i] != 0 && MuIsIso[i] == 1 ) fNumTotMuonBadIso++;
		else if( MuGood[i] != 0 && MuIsIso[i] == 0 ) fNumTotMuonBadNonIso++;
		if( MuGood[i] !=0 ){
			if( MuGood[i] / 100 != 0 ) fNumTotMuonDupl++;
			if( MuGood[i] % 100 / 10 != 0 ) fNumTotMuonNotPrimaryTrk++;
			int muClean = MuGood[i] % 10;
			if( muClean != 0 ){
				fNumTotMuonNotClean++;
				if( muClean == 1 ) fNumTotMuonBadDpop++;
				if( muClean == 2 ) fNumTotMuonBadChi2++;
				if( muClean == 3 ) fNumTotMuonBadNhit++;
			}
		}
	}

	fNumTotElectrons += NEles;
	for( int i=0; i < NEles; ++i ){
		if(      ElGood[i] == 0 && ElIsIso[i] == 1 ){ fNumTotElecGoodIso++; nObjects++; }
		else if( ElGood[i] == 0 && ElIsIso[i] == 0 ) fNumTotElecGoodNonIso++;
		else if( ElGood[i] != 0 && ElIsIso[i] == 1 ) fNumTotElecBadIso++;
		else if( ElGood[i] != 0 && ElIsIso[i] == 0 ) fNumTotElecBadNonIso++;
		if( ElGood[i] != 0 ){
			if( ElGood[i] / 100 != 0 ) fNumTotElecDupl++;
			if( ElGood[i] % 100 / 10 != 0 ) fNumTotElecNotPrimaryTrk++;
			int elClean = ElGood[i] % 10;
			if( elClean != 0 ){
				fNumTotElecNotClean++;
				if( elClean == 1 ) fNumTotElecBadHoE++;
				if( elClean == 2 ) fNumTotElecBadShsh++;
				if( elClean == 3 ) fNumTotElecBadTmat++;
			}
		}
	}

	fNumTotPhotons += NPhotons;
	for( int i=0; i < NPhotons; ++i ){
		if(      PhoGood[i] == 0 && PhoIsIso[i] == 1 ){ fNumTotPhotGoodIso++; nObjects++; }
		else if( PhoGood[i] == 0 && PhoIsIso[i] == 0 ) fNumTotPhotGoodNonIso++;
		else if( PhoGood[i] != 0 && PhoIsIso[i] == 1 ) fNumTotPhotBadIso++;
		else if( PhoGood[i] != 0 && PhoIsIso[i] == 0 ) fNumTotPhotBadNonIso++;
		if( PhoGood[i] != 0 ){
			int phClean = PhoGood[i] % 10;
			if( phClean != 0 ){
				fNumTotPhotNotClean++;
				if( phClean == 1 ) fNumTotPhotBadHoE++;
				if( phClean == 2 ) fNumTotPhotBadShsh++;
			}
		}
	}

	fNumTotJets += NJets;
	for( int i=0; i < NJets; ++i ){
		if( JGood[i] == 0 ){ fNumTotJetGood++; nObjects++; }
		else {
			fNumTotJetBad++;
			if( JGood[i] / 100 != 0 ) fNumTotJetDuplElJet++;
			if( JGood[i] % 100 / 10 != 0 ) fNumTotJetNotPrimaryTrk++;
			int jetClean = JGood[i] % 10;
			if( jetClean != 0 ){
				fNumTotJetNotClean++;
				if( jetClean == 1 ) fNumTotJetPgtE++;
				if( jetClean == 2 ) fNumTotJetGtFem++;
				if( jetClean == 3 ) fNumTotJetLtFem++;
				if( jetClean == 4 ) fNumTotJetLtFch++;
			}
		}
	}

	if( nObjects <= 0) fNumTotEvtEmpty++;
	int evClean = GoodEvent % 10;
	if( evClean != 0 || nObjects <= 0 ) fNumTotEvtReject++;
	if( evClean == 1 ) fNumTotEvtCleanEmpty++;
	if( evClean == 2 ) fNumTotEvtLtFem++;
	if( evClean == 3 ) fNumTotEvtLtFch++;

	evClean = GoodEvent / 1000;
	if( evClean == 1 ) fNumTotEvtPfMETJet++;
	if( evClean == 2 ) fNumTotEvtPfMETRij++;
	evClean = GoodEvent % 1000 / 100;
	if( evClean == 1 ) fNumTotEvtCaMETJet++;
	if( evClean == 2 ) fNumTotEvtCaMETRij++;
	evClean = GoodEvent % 100 / 10;
	if( evClean == 1 ) fNumTotEvtTcMETJet++;
	if( evClean == 2 ) fNumTotEvtTcMETRij++;
	return;
}

void TreeReader::StatPrint(void){
// prints the cleaning statistics
// to be called once at the end of the job


	cout << endl;
	cout << "Cleaning statistics from TreeReader " << endl;
	// cout << "  " << fNumTotElecGoodIso << "  " << fNumTotElecGoodNonIso
	//         << "  " << fNumTotElecBadIso  << "  " << fNumTotElecBadNonIso << endl;

	cout << endl;
	cout << " Total number of events processed = " << fNumTotEvt << endl;

// Statistics for events
	if (fNumTotEvt > 0) {
		cout << endl;
		cout << "   events accepted                  = " 
			<< fNumTotEvt-fNumTotEvtReject << endl;
		cout << "   events rejected (total)          = " << fNumTotEvtReject
			<< "  = " << 100.*(float)fNumTotEvtReject / (float)fNumTotEvt << " %" << endl;
		cout << endl;
		cout << "    empty after cleaning+Iso        = " << fNumTotEvtEmpty
			<< "  = " << 100.*(float)fNumTotEvtEmpty / (float)fNumTotEvt << " %" << endl;
		cout << "    empty after cleaning            = " << fNumTotEvtCleanEmpty
			<< "  = " << 100.*(float)fNumTotEvtCleanEmpty / (float)fNumTotEvt << " %" << endl;
		cout << "    too small em fraction           = " << fNumTotEvtLtFem
			<< "  = " << 100.*(float)fNumTotEvtLtFem / (float)fNumTotEvt << " %" << endl;
		cout << "    too small track pT fraction     = " << fNumTotEvtLtFch
			<< "  = " << 100.*(float)fNumTotEvtLtFch / (float)fNumTotEvt << " %" << endl;
		cout << "    jet with pT > " <<  "30" <<  " and bad        = " << fNumTotEvtBadHardJet
			<< "  = " << 100.*(float)fNumTotEvtBadHardJet / (float)fNumTotEvt << " %" << endl;
		cout << endl;
		cout << "    caloMET aligned with jet        = " << fNumTotEvtCaMETJet
			<< "  = " << 100.*(float)fNumTotEvtCaMETJet  / (float)fNumTotEvt << " %" << endl;
		cout << "    caloMET bad Rij                 = " << fNumTotEvtCaMETRij
			<< "  = " << 100.*(float)fNumTotEvtCaMETRij  / (float)fNumTotEvt << " %" << endl;
		cout << "    tcMET aligned with jet          = " << fNumTotEvtCaMETJet
			<< "  = " << 100.*(float)fNumTotEvtCaMETJet  / (float)fNumTotEvt << " %" << endl;
		cout << "    tcMET bad Rij                   = " << fNumTotEvtCaMETRij
			<< "  = " << 100.*(float)fNumTotEvtCaMETRij  / (float)fNumTotEvt << " %" << endl;
		cout << "    pfMET aligned with jet          = " << fNumTotEvtCaMETJet
			<< "  = " << 100.*(float)fNumTotEvtCaMETJet  / (float)fNumTotEvt << " %" << endl;
		cout << "    pfMET bad Rij                   = " << fNumTotEvtCaMETRij
			<< "  = " << 100.*(float)fNumTotEvtCaMETRij  / (float)fNumTotEvt << " %" << endl ;
	}

// Statistics for muons
	cout << endl;
	cout << " Total number of muons              = " << fNumTotMuons << endl;
	if (fNumTotMuons > 0) {
		cout << "  muons good+Iso                    = " << fNumTotMuonGoodIso
			<< "  = " << 100.*(float)fNumTotMuonGoodIso / (float)fNumTotMuons << " %" << endl;
		cout << "  muons good+non-Iso                = " << fNumTotMuonGoodNonIso
			<< "  = " << 100.*(float)fNumTotMuonGoodNonIso / (float)fNumTotMuons << " %" << endl;
		cout << "  muons bad+Iso                     = " << fNumTotMuonBadIso
			<< "  = " << 100.*(float)fNumTotMuonBadIso / (float)fNumTotMuons << " %" << endl;
		cout << "  muons bad+non-Iso                 = " << fNumTotMuonBadNonIso
			<< "  = " << 100.*(float)fNumTotMuonBadNonIso / (float)fNumTotMuons << " %" << endl;
		cout << endl;
		int mubad = fNumTotMuonBadIso + fNumTotMuonBadNonIso;
		cout << "  muons bad total                   = " << mubad
			<< "  = " << 100.*(float)mubad / (float)fNumTotMuons << " %" << endl;
		cout << "   muons duplicated                 = " << fNumTotMuonDupl
			<< "  = " << 100.*(float)fNumTotMuonDupl / (float)fNumTotMuons << " %" << endl;
		cout << "   muons not from Primary Vertex    = " << fNumTotMuonNotPrimaryTrk
			<< "  = " << 100.*(float)fNumTotMuonNotPrimaryTrk / (float)fNumTotMuons << " %" << endl;
		cout << "   muons not clean                  = " << fNumTotMuonNotClean
			<< "  = " << 100.*(float)fNumTotMuonNotClean / (float)fNumTotMuons << " %" << endl;
		cout << "    muons with bad dP/P             = " << fNumTotMuonBadDpop
			<< "  = " << 100.*(float)fNumTotMuonBadDpop / (float)fNumTotMuons << " %" << endl;
		cout << "    muons with bad chisquared       = " << fNumTotMuonBadChi2
			<< "  = " << 100.*(float)fNumTotMuonBadChi2 / (float)fNumTotMuons << " %" << endl;
		cout << "    muons with bad #tracker hits    = " << fNumTotMuonBadNhit
			<< "  = " << 100.*(float)fNumTotMuonBadNhit / (float)fNumTotMuons << " %" << endl;
	}

// Statistics for electrons
	cout << endl;
	cout << " Total number of electrons          = " << fNumTotElectrons << endl;
	if (fNumTotElectrons > 0) {
		cout << "  elecs good+Iso                    = " << fNumTotElecGoodIso
			<< "  = " << 100.*(float)fNumTotElecGoodIso / (float)fNumTotElectrons << " %" << endl;
		cout << "  elecs good+non-Iso                = " << fNumTotElecGoodNonIso
			<< "  = " << 100.*(float)fNumTotElecGoodNonIso / (float)fNumTotElectrons << " %" << endl;
		cout << "  elecs bad+Iso                     = " << fNumTotElecBadIso
			<< "  = " << 100.*(float)fNumTotElecBadIso / (float)fNumTotElectrons << " %" << endl;
		cout << "  elecs bad+non-Iso                 = " << fNumTotElecBadNonIso
			<< "  = " << 100.*(float)fNumTotElecBadNonIso / (float)fNumTotElectrons << " %" << endl;
		cout << endl;
		int elebad = fNumTotElecBadIso + fNumTotElecBadNonIso;
		cout << "  elecs bad total                   = " << elebad
			<< "  = " << 100.*(float)elebad / (float)fNumTotElectrons << " %" << endl;
		cout << "   elecs duplicated                 = " << fNumTotElecDupl
			<< "  = " << 100.*(float)fNumTotElecDupl / (float)fNumTotElectrons << " %" << endl;
		cout << "   elecs not from Primary Vertex    = " << fNumTotElecNotPrimaryTrk
			<< "  = " << 100.*(float)fNumTotElecNotPrimaryTrk / (float)fNumTotElectrons << " %" << endl;
		cout << "   elecs not clean                  = " << fNumTotElecNotClean
			<< "  = " << 100.*(float)fNumTotElecNotClean / (float)fNumTotElectrons << " %" << endl;
		cout << "    elecs with bad H/E              = " << fNumTotElecBadHoE
			<< "  = " << 100.*(float)fNumTotElecBadHoE / (float)fNumTotElectrons << " %" << endl;
		cout << "    elecs with bad shower shape     = " << fNumTotElecBadShsh
			<< "  = " << 100.*(float)fNumTotElecBadShsh / (float)fNumTotElectrons << " %" << endl;
		cout << "    elecs with bad track matching   = " << fNumTotElecBadTmat
			<< "  = " << 100.*(float)fNumTotElecBadTmat / (float)fNumTotElectrons << " %" << endl;
	}

// Statistics for photons
	cout << endl;
	cout << " Total number of photons            = " << fNumTotPhotons << endl;
	if (fNumTotPhotons > 0) {
		cout << "  photons good+Iso                  = " << fNumTotPhotGoodIso
			<< "  = " << 100.*(float)fNumTotPhotGoodIso / (float)fNumTotPhotons << " %" << endl;
		cout << "  photons good+non-Iso              = " << fNumTotPhotGoodNonIso
			<< "  = " << 100.*(float)fNumTotPhotGoodNonIso / (float)fNumTotPhotons << " %" << endl;
		cout << "  photons bad+Iso                   = " << fNumTotPhotBadIso
			<< "  = " << 100.*(float)fNumTotPhotBadIso / (float)fNumTotPhotons << " %" << endl;
		cout << "  photons bad+non-Iso               = " << fNumTotPhotBadNonIso
			<< "  = " << 100.*(float)fNumTotPhotBadNonIso / (float)fNumTotPhotons << " %" << endl;
		cout << endl;
		int photbad = fNumTotPhotBadIso + fNumTotPhotBadNonIso;
		cout << "  photons bad total                 = " << photbad
			<< "  = " << 100.*(float)photbad / (float)fNumTotPhotons << " %" << endl;
		cout << "   photons duplicated               = " << fNumTotPhotDupl
			<< "  = " << 100.*(float)fNumTotPhotDupl / (float)fNumTotPhotons << " %" << endl;
		cout << "   photons not clean                = " << fNumTotPhotNotClean
			<< "  = " << 100.*(float)fNumTotPhotNotClean / (float)fNumTotPhotons << " %" << endl;
		cout << "    photons with bad H/E            = " << fNumTotPhotBadHoE
			<< "  = " << 100.*(float)fNumTotPhotBadHoE / (float)fNumTotPhotons << " %" << endl;
		cout << "    photons with bad shower shape   = " << fNumTotPhotBadShsh
			<< "  = " << 100.*(float)fNumTotPhotBadShsh / (float)fNumTotPhotons << " %" << endl;
	}

// Statistics for jets
	cout << endl;
	cout << " Total number of jets               = " << fNumTotJets << endl;
	if (fNumTotJets > 0) {
		cout << "  jets good                         = " << fNumTotJetGood
			<< "  = " << 100.*(float)fNumTotJetGood / (float)fNumTotJets << " %" << endl;
		cout << "  jets bad                          = " << fNumTotJetBad
			<< "  = " << 100.*(float)fNumTotJetBad / (float)fNumTotJets << " %" << endl;
		cout << endl;
		cout << "   jets duplicated with elec        = " << fNumTotJetDuplElJet
			<< "  = " << 100.*(float)fNumTotJetDuplElJet / (float)fNumTotJets << " %" << endl;
		cout << "   jets not from Primary Vertex     = " << fNumTotJetNotPrimaryTrk
			<< "  = " << 100.*(float)fNumTotJetNotPrimaryTrk / (float)fNumTotJets << " %" << endl;
		cout << "   jets not clean                   = " << fNumTotJetNotClean
			<< "  = " << 100.*(float)fNumTotJetNotClean / (float)fNumTotJets << " %" << endl;
		cout << "    jets with P > E                 = " << fNumTotJetPgtE
			<< "  = " << 100.*(float)fNumTotJetPgtE / (float)fNumTotJets << " %" << endl;
		cout << "    jets with too large Fem         = " << fNumTotJetGtFem
			<< "  = " << 100.*(float)fNumTotJetGtFem / (float)fNumTotJets << " %" << endl;
		cout << "    jets with too small Fem         = " << fNumTotJetLtFem
			<< "  = " << 100.*(float)fNumTotJetLtFem / (float)fNumTotJets << " %" << endl;
		cout << "    jets with too small Fch         = " << fNumTotJetLtFch
			<< "  = " << 100.*(float)fNumTotJetLtFch / (float)fNumTotJets << " %" << endl;
		cout << "  Total number of b-jets            = " << fNumTotBJets
			<< "  = " << 100.*(float)fNumTotBJets / (float)fNumTotJets << " %" << endl;
	}
	return;  
}

void TreeReader::StatHistos(void){
// saves the cleaning statistics into histograms (from Kostas)
// to be called once at the end of the job

	fHstatHistos->SetFillColor(kRed);
	fHstatHistos->SetDrawOption("hbar");
	fHstatHistos->GetXaxis()->SetBinLabel(1,  "TotEvts");           fHstatHistos->SetBinContent(1,  fNumTotEvt);
	fHstatHistos->GetXaxis()->SetBinLabel(2,  "EvtReject");         fHstatHistos->SetBinContent(2,  fNumTotEvtReject);
	fHstatHistos->GetXaxis()->SetBinLabel(3,  "EvtEmpty");          fHstatHistos->SetBinContent(3,  fNumTotEvtEmpty);
	fHstatHistos->GetXaxis()->SetBinLabel(4,  "EvtCleanEmpty");     fHstatHistos->SetBinContent(4,  fNumTotEvtCleanEmpty);
	fHstatHistos->GetXaxis()->SetBinLabel(5,  "EvtLtFem");          fHstatHistos->SetBinContent(5,  fNumTotEvtLtFem);
	fHstatHistos->GetXaxis()->SetBinLabel(6,  "EvtLtFch");          fHstatHistos->SetBinContent(6,  fNumTotEvtLtFch);
	fHstatHistos->GetXaxis()->SetBinLabel(7,  "EvtPfMETJet");       fHstatHistos->SetBinContent(7,  fNumTotEvtPfMETJet);
	fHstatHistos->GetXaxis()->SetBinLabel(8,  "EvtPfMETRij");       fHstatHistos->SetBinContent(8,  fNumTotEvtPfMETRij);
	fHstatHistos->GetXaxis()->SetBinLabel(9,  "EvtCaMETJet");       fHstatHistos->SetBinContent(9,  fNumTotEvtCaMETJet);
	fHstatHistos->GetXaxis()->SetBinLabel(10, "EvtCaMETRij");       fHstatHistos->SetBinContent(10, fNumTotEvtCaMETRij);
	fHstatHistos->GetXaxis()->SetBinLabel(11, "EvtTcMETJet");       fHstatHistos->SetBinContent(11, fNumTotEvtTcMETJet);
	fHstatHistos->GetXaxis()->SetBinLabel(12, "EvtTcMETRij");       fHstatHistos->SetBinContent(12, fNumTotEvtTcMETRij);
	fHstatHistos->GetXaxis()->SetBinLabel(13, "EvtBadHardJet");     fHstatHistos->SetBinContent(13, fNumTotEvtBadHardJet);

	fHstatHistos->GetXaxis()->SetBinLabel(21, "TotMuons");          fHstatHistos->SetBinContent(21, fNumTotMuons);
	fHstatHistos->GetXaxis()->SetBinLabel(22, "MuonGoodIso");       fHstatHistos->SetBinContent(22, fNumTotMuonGoodIso);
	fHstatHistos->GetXaxis()->SetBinLabel(23, "MuonGoodNonIso");    fHstatHistos->SetBinContent(23, fNumTotMuonGoodNonIso);
	fHstatHistos->GetXaxis()->SetBinLabel(24, "MuonBadIso");        fHstatHistos->SetBinContent(24, fNumTotMuonBadIso);
	fHstatHistos->GetXaxis()->SetBinLabel(25, "MuonBadNonIso");     fHstatHistos->SetBinContent(25, fNumTotMuonBadNonIso);
	fHstatHistos->GetXaxis()->SetBinLabel(26, "MuonDupl");          fHstatHistos->SetBinContent(26, fNumTotMuonDupl);
	fHstatHistos->GetXaxis()->SetBinLabel(27, "MuonNotPrimaryTrk"); fHstatHistos->SetBinContent(27, fNumTotMuonNotPrimaryTrk);
	fHstatHistos->GetXaxis()->SetBinLabel(28, "MuonNotClean");      fHstatHistos->SetBinContent(28, fNumTotMuonNotClean);
	fHstatHistos->GetXaxis()->SetBinLabel(29, "MuonBadDpop");       fHstatHistos->SetBinContent(29, fNumTotMuonBadDpop);
	fHstatHistos->GetXaxis()->SetBinLabel(30, "MuonBadChi2");       fHstatHistos->SetBinContent(30, fNumTotMuonBadChi2);
	fHstatHistos->GetXaxis()->SetBinLabel(31, "MuonBadNhit");       fHstatHistos->SetBinContent(31, fNumTotMuonBadNhit);

	fHstatHistos->GetXaxis()->SetBinLabel(41, "TotElectrons");      fHstatHistos->SetBinContent(41, fNumTotElectrons);
	fHstatHistos->GetXaxis()->SetBinLabel(42, "ElecGoodIso");       fHstatHistos->SetBinContent(42, fNumTotElecGoodIso);
	fHstatHistos->GetXaxis()->SetBinLabel(43, "ElecGoodNonIso");    fHstatHistos->SetBinContent(43, fNumTotElecGoodNonIso);
	fHstatHistos->GetXaxis()->SetBinLabel(44, "ElecBadIso");        fHstatHistos->SetBinContent(44, fNumTotElecBadIso);
	fHstatHistos->GetXaxis()->SetBinLabel(45, "ElecBadNonIso");     fHstatHistos->SetBinContent(45, fNumTotElecBadNonIso);
	fHstatHistos->GetXaxis()->SetBinLabel(46, "ElecDupl");          fHstatHistos->SetBinContent(46, fNumTotElecDupl);
	fHstatHistos->GetXaxis()->SetBinLabel(47, "ElecNotPrimaryTrk"); fHstatHistos->SetBinContent(47, fNumTotElecNotPrimaryTrk);
	fHstatHistos->GetXaxis()->SetBinLabel(48, "ElecNotClean");      fHstatHistos->SetBinContent(48, fNumTotElecNotClean);
	fHstatHistos->GetXaxis()->SetBinLabel(49, "ElecBadHoE");        fHstatHistos->SetBinContent(49, fNumTotElecBadHoE);
	fHstatHistos->GetXaxis()->SetBinLabel(50, "ElecBadShsh");       fHstatHistos->SetBinContent(50, fNumTotElecBadShsh);
	fHstatHistos->GetXaxis()->SetBinLabel(51, "ElecBadTmat");       fHstatHistos->SetBinContent(51, fNumTotElecBadTmat);

	fHstatHistos->GetXaxis()->SetBinLabel(61, "TotPhotons");        fHstatHistos->SetBinContent(61, fNumTotPhotons);
	fHstatHistos->GetXaxis()->SetBinLabel(62, "PhotGoodIso");       fHstatHistos->SetBinContent(62, fNumTotPhotGoodIso);
	fHstatHistos->GetXaxis()->SetBinLabel(63, "PhotGoodNonIso");    fHstatHistos->SetBinContent(63, fNumTotPhotGoodNonIso);
	fHstatHistos->GetXaxis()->SetBinLabel(64, "PhotBadIso");        fHstatHistos->SetBinContent(64, fNumTotPhotBadIso);
	fHstatHistos->GetXaxis()->SetBinLabel(65, "PhotBadNonIso");     fHstatHistos->SetBinContent(65, fNumTotPhotBadNonIso);
	fHstatHistos->GetXaxis()->SetBinLabel(66, "PhotDupl");          fHstatHistos->SetBinContent(66, fNumTotPhotDupl);
	fHstatHistos->GetXaxis()->SetBinLabel(67, "PhotNotClean");      fHstatHistos->SetBinContent(67, fNumTotPhotNotClean);
	fHstatHistos->GetXaxis()->SetBinLabel(68, "PhotBadHoE");        fHstatHistos->SetBinContent(68, fNumTotPhotBadHoE);
	fHstatHistos->GetXaxis()->SetBinLabel(69, "PhotBadShsh");       fHstatHistos->SetBinContent(69, fNumTotPhotBadShsh);

	fHstatHistos->GetXaxis()->SetBinLabel(81, "TotJets");           fHstatHistos->SetBinContent(81, fNumTotJets);
	fHstatHistos->GetXaxis()->SetBinLabel(82, "JetGood");           fHstatHistos->SetBinContent(82, fNumTotJetGood);
	fHstatHistos->GetXaxis()->SetBinLabel(83, "JetBad");            fHstatHistos->SetBinContent(83, fNumTotJetBad);
	fHstatHistos->GetXaxis()->SetBinLabel(84, "JetDuplElJet");      fHstatHistos->SetBinContent(84, fNumTotJetDuplElJet);
	fHstatHistos->GetXaxis()->SetBinLabel(85, "JetNotPrimaryTrk");  fHstatHistos->SetBinContent(85, fNumTotJetNotPrimaryTrk);
	fHstatHistos->GetXaxis()->SetBinLabel(86, "JetNotClean");       fHstatHistos->SetBinContent(86, fNumTotJetNotClean);
	fHstatHistos->GetXaxis()->SetBinLabel(87, "JetPgtE");           fHstatHistos->SetBinContent(87, fNumTotJetPgtE);
	fHstatHistos->GetXaxis()->SetBinLabel(88, "JetGtFem");          fHstatHistos->SetBinContent(88, fNumTotJetGtFem);
	fHstatHistos->GetXaxis()->SetBinLabel(89, "JetLtFem");          fHstatHistos->SetBinContent(89, fNumTotJetLtFem);
	fHstatHistos->GetXaxis()->SetBinLabel(90, "JetLtFch");          fHstatHistos->SetBinContent(90, fNumTotJetLtFch);
	fHstatHistos->GetXaxis()->SetBinLabel(91, "TotBJets");          fHstatHistos->SetBinContent(91, fNumTotBJets);

	return;
}

///////////////////////////////////////////////////////////////////////////////////////////////
// Utilities //////////////////////////////////////////////////////////////////////////////////
void TreeReader::PutMuon(int inew, int iold){
	// This needs to be UPDATED every time the tree content changes!
	MuGood        [inew] = MuGood        [iold];
	MuIsIso       [inew] = MuIsIso       [iold];
	MuPx          [inew] = MuPx          [iold];
	MuPy          [inew] = MuPy          [iold];
	MuPz          [inew] = MuPz          [iold];
	MuPt          [inew] = MuPt          [iold];
	MuPtE         [inew] = MuPtE         [iold];
	MuE           [inew] = MuE           [iold];
	MuEt          [inew] = MuEt          [iold];
	MuEta         [inew] = MuEta         [iold];
	MuPhi         [inew] = MuPhi         [iold];
	MuCharge      [inew] = MuCharge      [iold];
	MuRelIso03    [inew] = MuRelIso03    [iold];
	MuIso03SumPt  [inew] = MuIso03SumPt  [iold];
	MuIso03EmEt   [inew] = MuIso03EmEt   [iold];
	MuIso03HadEt  [inew] = MuIso03HadEt  [iold];
	MuIso05SumPt  [inew] = MuIso05SumPt  [iold];
	MuIso05EmEt   [inew] = MuIso05EmEt   [iold];
	MuIso05HadEt  [inew] = MuIso05HadEt  [iold];
	MuEem         [inew] = MuEem         [iold];
	MuEhad        [inew] = MuEhad        [iold];
	MuD0BS        [inew] = MuD0BS        [iold];
	MuD0PV        [inew] = MuD0PV        [iold];
	MuD0E         [inew] = MuD0E         [iold];
	MuDzBS        [inew] = MuDzBS        [iold];
	MuDzPV        [inew] = MuDzPV        [iold];
	MuDzE         [inew] = MuDzE         [iold];
	MuNChi2       [inew] = MuNChi2       [iold];
	MuNGlHits     [inew] = MuNGlHits     [iold];
	MuNMuHits     [inew] = MuNMuHits     [iold];
	MuNTkHits     [inew] = MuNTkHits     [iold];
	MuNMatches    [inew] = MuNMatches    [iold];
	MuNChambers   [inew] = MuNChambers   [iold];
	MuCaloComp    [inew] = MuCaloComp    [iold];
	MuSegmComp    [inew] = MuSegmComp    [iold];
	MuTrackerMu   [inew] = MuTrackerMu   [iold];
	MuGMPT        [inew] = MuGMPT        [iold];
	MuID          [inew] = MuID          [iold];
	MuMID         [inew] = MuMID         [iold];
}

void TreeReader::PutElectron(int inew, int iold){
	// This needs to be UPDATED every time the tree content changes!
	ElGood                      [inew] = ElGood                      [iold];
	ElIsIso                     [inew] = ElIsIso                     [iold];
	ElPx                        [inew] = ElPx                        [iold];
	ElPy                        [inew] = ElPy                        [iold];
	ElPz                        [inew] = ElPz                        [iold];
	ElPt                        [inew] = ElPt                        [iold];
	ElPtE                       [inew] = ElPtE                       [iold];
	ElE                         [inew] = ElE                         [iold];
	ElEt                        [inew] = ElEt                        [iold];
	ElEta                       [inew] = ElEta                       [iold];
	ElTheta                     [inew] = ElTheta                       [iold];
	ElPhi                       [inew] = ElPhi                       [iold];
	ElD0BS                      [inew] = ElD0BS                      [iold];
	ElD0PV                      [inew] = ElD0PV                      [iold];
	ElD0E                       [inew] = ElD0E                       [iold];
	ElDzBS                      [inew] = ElDzBS                      [iold];
	ElDzPV                      [inew] = ElDzPV                      [iold];
	ElDzE                       [inew] = ElDzE                       [iold];
	ElIso                       [inew] = ElIso                       [iold];
	ElPtSum                     [inew] = ElPtSum                     [iold];
	ElEmEtSum                   [inew] = ElEmEtSum                   [iold];
	ElHadEtSum                  [inew] = ElHadEtSum                  [iold];
	ElNChi2                     [inew] = ElNChi2                     [iold];
	ElCharge                    [inew] = ElCharge                    [iold];
	ElIDTight                   [inew] = ElIDTight                   [iold];
	ElIDLoose                   [inew] = ElIDLoose                   [iold];
	ElIDRobustTight             [inew] = ElIDRobustTight             [iold];
	ElIDRobustLoose             [inew] = ElIDRobustLoose             [iold];
	ElInGap                     [inew] = ElInGap                     [iold];
	ElEcalDriven                [inew] = ElEcalDriven                [iold];
	ElTrackerDriven             [inew] = ElTrackerDriven             [iold];
	ElBasicClustersSize         [inew] = ElBasicClustersSize         [iold];
	Elfbrem                     [inew] = Elfbrem                     [iold];
	ElHcalOverEcal              [inew] = ElHcalOverEcal              [iold];
	ElE5x5                      [inew] = ElE5x5                      [iold];
	ElE2x5Max                   [inew] = ElE2x5Max                   [iold];
	ElSigmaIetaIeta             [inew] = ElSigmaIetaIeta             [iold];
	ElDeltaPhiSeedClusterAtCalo [inew] = ElDeltaPhiSeedClusterAtCalo [iold];
	ElDeltaEtaSeedClusterAtCalo [inew] = ElDeltaEtaSeedClusterAtCalo [iold];
	ElDeltaPhiSuperClusterAtVtx [inew] = ElDeltaPhiSuperClusterAtVtx [iold];
	ElDeltaEtaSuperClusterAtVtx [inew] = ElDeltaEtaSuperClusterAtVtx [iold];
	ElCaloEnergy                [inew] = ElCaloEnergy                [iold];
	ElTrkMomAtVtx               [inew] = ElTrkMomAtVtx               [iold];
	ElESuperClusterOverP        [inew] = ElESuperClusterOverP        [iold];
	ElIsInJet                   [inew] = ElIsInJet                   [iold];
	ElSharedPx                  [inew] = ElSharedPx                  [iold];
	ElSharedPy                  [inew] = ElSharedPy                  [iold];
	ElSharedPz                  [inew] = ElSharedPz                  [iold];
	ElSharedEnergy              [inew] = ElSharedEnergy              [iold];
	ElDuplicateEl               [inew] = ElDuplicateEl               [iold];
	ElDR03TkSumPt               [inew] = ElDR03TkSumPt               [iold];
	ElDR04TkSumPt               [inew] = ElDR04TkSumPt               [iold];
	ElDR03EcalRecHitSumEt       [inew] = ElDR03EcalRecHitSumEt       [iold];
	ElDR04EcalRecHitSumEt       [inew] = ElDR04EcalRecHitSumEt       [iold];
	ElDR03HcalTowerSumEt        [inew] = ElDR03HcalTowerSumEt        [iold];
	ElDR04HcalTowerSumEt        [inew] = ElDR04HcalTowerSumEt        [iold];
	ElID                        [inew] = ElID                        [iold];
	ElMID                       [inew] = ElMID                       [iold];
}

void TreeReader::PutPhoton(int inew, int iold){
	// This needs to be UPDATED every time the tree content changes!
	PhoGood           [inew] = PhoGood           [iold];
	PhoIsIso          [inew] = PhoIsIso          [iold];
	PhoPt             [inew] = PhoPt             [iold];
	PhoPx             [inew] = PhoPx             [iold];
	PhoPy             [inew] = PhoPy             [iold];
	PhoPz             [inew] = PhoPz             [iold];
	PhoEta            [inew] = PhoEta            [iold];
	PhoPhi            [inew] = PhoPhi            [iold];
	PhoEnergy         [inew] = PhoEnergy         [iold];
	PhoIso03Ecal      [inew] = PhoIso03Ecal      [iold];
	PhoIso03Hcal      [inew] = PhoIso03Hcal      [iold];
	PhoIso03TrkSolid  [inew] = PhoIso03TrkSolid  [iold];
	PhoIso03TrkHollow [inew] = PhoIso03TrkHollow [iold];
	PhoIso03          [inew] = PhoIso03          [iold];
	PhoCaloPositionX  [inew] = PhoCaloPositionX  [iold];
	PhoCaloPositionY  [inew] = PhoCaloPositionY  [iold];
	PhoCaloPositionZ  [inew] = PhoCaloPositionZ  [iold];
	PhoHoverE         [inew] = PhoHoverE         [iold];
	PhoH1overE        [inew] = PhoH1overE        [iold];
	PhoH2overE        [inew] = PhoH2overE        [iold];
	PhoSigmaIetaIeta  [inew] = PhoSigmaIetaIeta  [iold];
	PhoHasPixSeed     [inew] = PhoHasPixSeed     [iold];
	PhoHasConvTrks    [inew] = PhoHasConvTrks    [iold];
	PhoIsInJet        [inew] = PhoIsInJet        [iold];
	PhoSharedPx       [inew] = PhoSharedPx       [iold];
	PhoSharedPy       [inew] = PhoSharedPy       [iold];
	PhoSharedPz       [inew] = PhoSharedPz       [iold];
	PhoSharedEnergy   [inew] = PhoSharedEnergy   [iold];
}

void TreeReader::PutJet(int inew, int iold){
	// This needs to be UPDATED every time the tree content changes!
	JGood          [inew] = JGood          [iold];
	JPx            [inew] = JPx            [iold];
	JPy            [inew] = JPy            [iold];
	JPz            [inew] = JPz            [iold];
	JPt            [inew] = JPt            [iold];
	JE             [inew] = JE             [iold];
	JEt            [inew] = JEt            [iold];
	JEta           [inew] = JEta           [iold];
	JPhi           [inew] = JPhi           [iold];
	JEMfrac        [inew] = JEMfrac        [iold];
	JNConstituents [inew] = JNConstituents [iold];
	JID_HPD        [inew] = JID_HPD        [iold];
	JID_RBX        [inew] = JID_RBX        [iold];
	JID_n90Hits    [inew] = JID_n90Hits    [iold];
	JID_SubDet1    [inew] = JID_SubDet1    [iold];
	JID_SubDet2    [inew] = JID_SubDet2    [iold];
	JID_SubDet3    [inew] = JID_SubDet3    [iold];
	JID_SubDet4    [inew] = JID_SubDet4    [iold];
	JID_resEMF     [inew] = JID_resEMF     [iold];
	JID_HCALTow    [inew] = JID_HCALTow    [iold];
	JID_ECALTow    [inew] = JID_ECALTow    [iold];
	JEtaEMrms      [inew] = JEtaEMrms      [iold];
	JEtaHADrms     [inew] = JEtaHADrms     [iold];
	JPhiEMrms      [inew] = JPhiEMrms      [iold];
	JPhiHADrms     [inew] = JPhiHADrms     [iold];
	JbTagProb      [inew] = JbTagProb      [iold];
	JChfrac        [inew] = JChfrac        [iold];
	JMass          [inew] = JMass          [iold];
	JNAssoTracks   [inew] = JNAssoTracks   [iold];
	Jtrk1px        [inew] = Jtrk1px        [iold];
	Jtrk1py        [inew] = Jtrk1py        [iold];
	Jtrk1pz        [inew] = Jtrk1pz        [iold];
	Jtrk2px        [inew] = Jtrk2px        [iold];
	Jtrk2py        [inew] = Jtrk2py        [iold];
	Jtrk2pz        [inew] = Jtrk2pz        [iold];
	Jtrk3px        [inew] = Jtrk3px        [iold];
	Jtrk3py        [inew] = Jtrk3py        [iold];
	Jtrk3pz        [inew] = Jtrk3pz        [iold];
	JEcorr         [inew] = JEcorr         [iold];
	JeMinDR        [inew] = JeMinDR        [iold];
	JVtxx          [inew] = JVtxx          [iold];
	JVtxy          [inew] = JVtxy          [iold];
	JVtxz          [inew] = JVtxz          [iold];
	JVtxExx        [inew] = JVtxExx        [iold];
	JVtxEyx        [inew] = JVtxEyx        [iold];
	JVtxEyy        [inew] = JVtxEyy        [iold];
	JVtxEzy        [inew] = JVtxEzy        [iold];
	JVtxEzz        [inew] = JVtxEzz        [iold];
	JVtxEzx        [inew] = JVtxEzx        [iold];
	JVtxNChi2      [inew] = JVtxNChi2      [iold];
}

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
