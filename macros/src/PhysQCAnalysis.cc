#include <TChainElement.h>

#include "base/TreeReader.hh"
#include "helper/AnaClass.hh"
#include "base/TreeAnalyzerBase.hh"
#include "TreeCleaner.hh"
#include "helper/Utilities.hh"

#include "PhysQCAnalysis.hh"

using namespace std;

PhysQCAnalysis::PhysQCAnalysis(TreeReader *tr, TreeCleaner *tc) : UserAnalysisBase(tr){
	fTC = tc;
	fTlat = new TLatex();
	fAC = new AnaClass();
	Util::SetStyle();

	// Histos for uncleaned objects
	const unsigned int nmubins = 100;
	fMuHistos[0] = new TH1D("MuDeltaPOverP", "Mu DeltaP/P;Mu DeltaP/P",    nmubins, 0., 0.25 );
	fMuHistos[1] = new TH1D("Mud0signif", "Mu d0 significance;Mu d0 significance",  nmubins, 0., 15. );
	fMuHistos[2] = new TH1D("Mudzsignif", "Mu dz significance;Mu dz significance",  nmubins, 0., 15. );
	fMuHistos[3] = new TH1D("MuDRSS", "Mu DR Same Sign;Mu DR Same Sign",  nmubins, 0., 3.2 );
	
	const unsigned int nelbins = 100;
	fElHistos[0]  = new TH1D("ElHcalOverEcalBar", "El H/E Bar",  nelbins, 0., 0.05);
	fElHistos[1]  = new TH1D("ElSigmaIetaIetaBar", "El sigmaIetaIeta Bar",  nelbins, 0., 0.05);
	fElHistos[2]  = new TH1D("ElDeltaPhiSeedClusterAtCaloBar", "El DeltaPhiSeedClusterAtCalo Bar",  nelbins, -0.1, 0.1);
	fElHistos[3]  = new TH1D("ElDeltaEtaSeedClusterAtCaloBar", "El DeltaEtaSeedClusterAtCalo Bar",  nelbins, -0.05, 0.05);
	fElHistos[4]  = new TH1D("ElDeltaPhiSuperClusterAtVtxBar", "El DeltaPhiSuperClusterAtVtx Bar",  nelbins, -0.1, 0.1);
	fElHistos[5]  = new TH1D("ElDeltaEtaSuperClusterAtVtxBar", "El DeltaEtaSuperClusterAtVtx Bar",  nelbins, -0.05, 0.05);
	fElHistos[6]  = new TH1D("ElESuperClusterOverPBar", "El ESuperClusterOverP Bar",  nelbins, 0., 5.);
	fElHistos[7]  = new TH1D("ElHcalOverEcalEnd", "El H/E End",  nelbins, 0., 0.05);
	fElHistos[8]  = new TH1D("ElSigmaIetaIetaEnd", "El sigmaIetaIeta End",  nelbins, 0., 0.05);
	fElHistos[9]  = new TH1D("ElDeltaPhiSeedClusterAtCaloEnd", "El DeltaPhiSeedClusterAtCalo End",  nelbins, -0.1, 0.1);
	fElHistos[10] = new TH1D("ElDeltaEtaSeedClusterAtCaloEnd", "El DeltaEtaSeedClusterAtCalo End",  nelbins, -0.05, 0.05);
	fElHistos[11] = new TH1D("ElDeltaPhiSuperClusterAtVtxEnd", "El DeltaPhiSuperClusterAtVtx End",  nelbins, -0.1, 0.1);
	fElHistos[12] = new TH1D("ElDeltaEtaSuperClusterAtVtxEnd", "El DeltaEtaSuperClusterAtVtx End",  nelbins, -0.05, 0.05);
	fElHistos[13] = new TH1D("ElESuperClusterOverPEnd", "El ESuperClusterOverP End",  nelbins, 0., 5.);
	fElHistos[14] = new TH1D("Eld0signif", "El d0 significance;El d0 significance",  nelbins, 0., 15. );
	fElHistos[15] = new TH1D("Eldzsignif", "El dz significance;El dz significance",  nelbins, 0., 15. );
	fElHistos[16] = new TH1D("ElDRSS", "El DR Same Sign;El DR Same Sign",  nelbins, 0., 3.2 );
	fElHistos[17] = new TH1D("ElDROS", "El DR Opp Sign;El DR Opp Sign",  nelbins, 0., 3.2 );
	fElHistos[18] = new TH1D("ElInvMSS", "El InvM Same Sign;El InvM Same Sign", nelbins, 0., 10. );
	fElHistos[19] = new TH1D("ElInvMOS", "El InvM Opp Sign;El InvM Opp Sign",  nelbins, 0., 10. );

	const unsigned int nphbins = 100;
	fPhHistos[0]  = new TH1D("PhoHoverEBar", "Phot H/E Bar",  nphbins, 0., 0.05);
	fPhHistos[1]  = new TH1D("PhoSigmaIetaIetaBar", "Phot sigmaIetaIeta Bar",  nphbins, 0., 0.05);
	fPhHistos[2]  = new TH1D("PhoHoverEEnd", "Phot H/E End",  nphbins, 0., 0.05);
	fPhHistos[3]  = new TH1D("PhoSigmaIetaIetaEnd", "Phot sigmaIetaIeta End",  nphbins, 0., 0.05);
	fPhHistos[4] = new TH1D("PhDR", "Phot DR;Phot DR",  nphbins, 0., 3.2 );
	fPhHistos[5] = new TH1D("PhInvM", "Phot InvM;Phot InvM", nphbins, 0., 15. );

	const unsigned int njetbins = 100;
	fJHistos[0] = new TH1D("Jetd0PV", "Jet d0 PV;Jet d0 to PV",  njetbins, -0.5, 0.5 );
	fJHistos[1] = new TH1D("JetdzPV", "Jet dz PV;Jet dz to PV",  njetbins, -2., 2. );
	fJHistos[2] = new TH1D("Jetd0signif", "Jet d0 significance;Jet d0 significance",  njetbins, 0., 15. );
	fJHistos[3] = new TH1D("Jetdzsignif", "Jet dz significance;Jet dz significance",  njetbins, 0., 15. );
	
	const unsigned int nevbins = 100;
	fMETHistos[0] = new TH1D("DphiMETJet1", "Dphi(MET,j1)", nevbins, 0, 3.1416);
	fMETHistos[1] = new TH1D("DphiMETJet2", "Dphi(MET,j2)", nevbins, 0, 3.1416);
	fMETHistos[2] = new TH1D("METR12", "MET R12;MET R12", nevbins, 0, 6.283);
	fMETHistos[3] = new TH1D("METR21", "MET R21;MET R21", nevbins, 0, 6.283);
	fMETHistos[4] = new TH1D("METRsum", "MET Rsum;MET Rsum", nevbins, 0, 6.283);
	fMETHistos[5] = new TH1D("EvtEmFrac", "Evt Em Frac;Evt Em Frac", nevbins, 0, 1.05);
	fMETHistos[6] = new TH1D("EvtChFrac", "Evt Ch Frac;Evt Ch Frac", nevbins, 0, 1.05);
	fMETDphi12    = new TH2D("METDphi12", "MET Dphi 12;Dphi(MET,j1);Dphi(MET,j2)", nevbins, 0, 3.142, nevbins, 0, 3.142);
	fMETDphi12->SetMarkerSize(7);

	// Histos for cleaned objects
}

PhysQCAnalysis::~PhysQCAnalysis(){
}

void PhysQCAnalysis::Begin(){
	fAC->readVarNames("varnames.dat");
	fAC->setOutputDir(fOutputDir);
	fAC->setGlobalTag(fTag);
}

void PhysQCAnalysis::Analyze1(){
// Analysis before cleaning
	
	// skip empty events
//	if (fTR->GoodEvent > 0) return;
	if(fTR->NMus + fTR->NEles + fTR->NPhotons + fTR->NJets <= 0) return;
	
	double p1[4], p2[4], minv;
	
	double drVxsq = fTR->PrimVtxxE*fTR->PrimVtxxE + fTR->PrimVtxyE*fTR->PrimVtxyE;

	for(size_t i = 0; i < fTR->NMus; ++i){
		fMuHistos[0]->Fill(fTR->MuPtE[i]/fTR->MuPt[i]);
		double d0  = fTR->MuD0PV[i];
		double dd0 = sqrt(fTR->MuD0E[i]*fTR->MuD0E[i] + drVxsq);
		if (dd0 <= 0.) dd0 = 0.001;
		double dz  = fTR->MuDzPV[i];
		double ddz = sqrt(fTR->MuDzE[i]*fTR->MuDzE[i] + fTR->PrimVtxzE*fTR->PrimVtxzE);
		if (ddz <= 0.) ddz = 0.001;
		fMuHistos[1]->Fill(fabs(d0) / dd0);
		fMuHistos[2]->Fill(fabs(dz) / ddz);
	}

	for(size_t i = 0; i < fTR->NMus; ++i){
		for(size_t j = i+1; j < fTR->NMus; ++j){
			if (fTR->MuCharge[i] == fTR->MuCharge[j]) {
				double deltaR = Util::GetDeltaR(fTR->MuEta[i], fTR->MuEta[j], fTR->MuPhi[i], fTR->MuPhi[j]);
				fMuHistos[3]->Fill(deltaR);
			}
		}
	}

	for(size_t i = 0; i < fTR->NEles; ++i){
		if (fTR->ElEta[i] < 1.479) {
			fElHistos[0]->Fill(fTR->ElHcalOverEcal[i]);
			fElHistos[1]->Fill(fTR->ElSigmaIetaIeta[i]);
			fElHistos[2]->Fill(fTR->ElDeltaPhiSeedClusterAtCalo[i]);
			fElHistos[3]->Fill(fTR->ElDeltaEtaSeedClusterAtCalo[i]);
			fElHistos[4]->Fill(fTR->ElDeltaPhiSuperClusterAtVtx[i]);
			fElHistos[5]->Fill(fTR->ElDeltaEtaSuperClusterAtVtx[i]);
			fElHistos[6]->Fill(fTR->ElESuperClusterOverP[i]);
		} else {
			fElHistos[7]->Fill(fTR->ElHcalOverEcal[i]);
			fElHistos[8]->Fill(fTR->ElSigmaIetaIeta[i]);
			fElHistos[9]->Fill(fTR->ElDeltaPhiSeedClusterAtCalo[i]);
			fElHistos[10]->Fill(fTR->ElDeltaEtaSeedClusterAtCalo[i]);
			fElHistos[11]->Fill(fTR->ElDeltaPhiSuperClusterAtVtx[i]);
			fElHistos[12]->Fill(fTR->ElDeltaEtaSuperClusterAtVtx[i]);
			fElHistos[13]->Fill(fTR->ElESuperClusterOverP[i]);
		}
		double d0  = fTR->ElD0PV[i];
		double dd0 = sqrt(fTR->ElD0E[i]*fTR->ElD0E[i] + drVxsq);
		if (dd0 <= 0.) dd0 = 0.001;
		double dz  = fTR->ElDzPV[i];
		double ddz = sqrt(fTR->ElDzE[i]*fTR->ElDzE[i] + fTR->PrimVtxzE*fTR->PrimVtxzE);
		if (ddz <= 0.) ddz = 0.001;
		fElHistos[14]->Fill(fabs(d0) / dd0);
		fElHistos[15]->Fill(fabs(dz) / ddz);
	}

	for(size_t i = 0; i < fTR->NEles; ++i){
		for(size_t j = i+1; j < fTR->NEles; ++j){
			double deltaR = Util::GetDeltaR(fTR->ElEta[i], fTR->ElEta[j], fTR->ElPhi[i], fTR->ElPhi[j]);
			if (fTR->ElCharge[i] == fTR->ElCharge[j]) {
				fElHistos[16]->Fill(deltaR);
				if (deltaR < 0.5) {
					p1[0] = fTR->ElPx[i];
					p1[1] = fTR->ElPy[i];
					p1[2] = fTR->ElPz[i];
					p1[3] = fTR->ElE[i];
					p2[0] = fTR->ElPx[j];
					p2[1] = fTR->ElPy[j];
					p2[2] = fTR->ElPz[j];
					p2[3] = fTR->ElE[j];
					minv = invMass(p1, p2);
					fElHistos[18]->Fill(minv);
				}
			} else {
				fElHistos[17]->Fill(deltaR);
				if (deltaR < 0.5) {
					p1[0] = fTR->ElPx[i];
					p1[1] = fTR->ElPy[i];
					p1[2] = fTR->ElPz[i];
					p1[3] = fTR->ElE[i];
					p2[0] = fTR->ElPx[j];
					p2[1] = fTR->ElPy[j];
					p2[2] = fTR->ElPz[j];
					p2[3] = fTR->ElE[j];
					minv = invMass(p1, p2);
					fElHistos[19]->Fill(minv);
				}
			}
		}
	}
	
	for(size_t i = 0; i < fTR->NPhotons; ++i){
		if (fTR->PhoEta[i] < 1.479) {
			fPhHistos[0]->Fill(fTR->PhoHoverE[i]);
			fPhHistos[1]->Fill(fTR->PhoSigmaIetaIeta[i]);
		} else {
			fPhHistos[2]->Fill(fTR->PhoHoverE[i]);
			fPhHistos[3]->Fill(fTR->PhoSigmaIetaIeta[i]);
		}
	}
	
	for(size_t i = 0; i < fTR->NPhotons; ++i){
		for(size_t j = i+1; j < fTR->NPhotons; ++j){
			double deltaR = Util::GetDeltaR(fTR->PhoEta[i], fTR->PhoEta[j], fTR->PhoPhi[i], fTR->PhoPhi[j]);
			fPhHistos[4]->Fill(deltaR);
			if (deltaR < 0.5) {
				p1[0] = fTR->PhoPx[i];
				p1[1] = fTR->PhoPy[i];
				p1[2] = fTR->PhoPz[i];
				p1[3] = fTR->PhoEnergy[i];
				p2[0] = fTR->PhoPx[j];
				p2[1] = fTR->PhoPy[j];
				p2[2] = fTR->PhoPz[j];
				p2[3] = fTR->PhoEnergy[j];
				minv = invMass(p1, p2);
				fPhHistos[5]->Fill(minv);
			}
		}
	}
	
	for(size_t i = 0; i < fTR->NJets; ++i){
		double d0  = 0.;
		double dd0 = 0.001;
		double dz  = 0.;
		double ddz = 0.001;
//		if (fTR->JVtxNChi2 > 0) { // not defined in this file
		double dztry = fTR->JVtxz[i] - fTR->PrimVtxz;
		if (dztry > -100.) {
			d0  = sqrt((fTR->JVtxx[i]-fTR->PrimVtxx)*(fTR->JVtxx[i]-fTR->PrimVtxx)
			+ (fTR->JVtxy[i]-fTR->PrimVtxy)*(fTR->JVtxy[i]-fTR->PrimVtxy) );
			dd0 = sqrt(fTR->JVtxExx[i] + fTR->JVtxEyy[i] + drVxsq);
			if (dd0 <= 0.) dd0 = 0.001;
			dz  = fTR->JVtxz[i] - fTR->PrimVtxz;
			ddz = sqrt(fTR->JVtxEzz[i] + fTR->PrimVtxzE*fTR->PrimVtxzE);
			if (ddz <= 0.) ddz = 0.001;
		}
		fJHistos[0]->Fill(d0);
		fJHistos[1]->Fill(dz);
		fJHistos[2]->Fill(fabs(d0) / dd0);
		fJHistos[3]->Fill(fabs(dz) / ddz);
//		if (dz < -2.) {
//		cout << " dz = " << dz << ", JVtxNChi2 = " << fTR->JVtxNChi2
//		     << ", JNAssoTracks = " << fTR->JNAssoTracks << endl;
//		}
	}

	double METPhi = fTR->MuCorrMETphi;
	double MET = fTR->MuCorrMET;
	double METBadJetmin = 20.;
	if (fTR->NJets > 1 && MET > METBadJetmin) {
		double dPhiMJ1 = Util::DeltaPhi(fTR->JPhi[0], METPhi);
		double dPhiMJ2 = Util::DeltaPhi(fTR->JPhi[1], METPhi);
		fMETHistos[0]->Fill(dPhiMJ1);
		fMETHistos[1]->Fill(dPhiMJ2);
		fMETDphi12->Fill(dPhiMJ1, dPhiMJ2);
		double pi = 3.141592654;
		double fR12 = sqrt(dPhiMJ1*dPhiMJ1 + (pi-dPhiMJ2)*(pi-dPhiMJ2) );
		double fR21 = sqrt(dPhiMJ2*dPhiMJ2 + (pi-dPhiMJ1)*(pi-dPhiMJ1) );
		fMETHistos[2]->Fill(fR12);
		fMETHistos[3]->Fill(fR21);
		fMETHistos[4]->Fill(fR12+fR21);
	}
	
	double fracEm, fracCh;
	GetEvtEmChFrac(fracEm, fracCh);
	fMETHistos[5]->Fill(fracEm);
	fMETHistos[6]->Fill(fracCh);
}

void PhysQCAnalysis::Analyze2(){
// Analysis after cleaning
	
	// skip empty events
	if (fTR->NMus + fTR->NEles + fTR->NPhotons + fTR->NJets <= 0) return;
	
}

void PhysQCAnalysis::End(){

	double distVxmax = 10.;
	double dRSSmuonmax = 0.1;
	double ElecHoverEBarmax      = 0.045;
	double ElecHoverEEndmax      = 0.05;
	double ElecSigmaEtaEtaBarmax = 0.011;
	double ElecSigmaEtaEtaEndmax = 0.025;
	double ElecEoverPInBarmin    = 0.3;
	double ElecEoverPInEndmin    = 0.4;
	double ElecDeltaEtaInBarmax  = 0.007;
	double ElecDeltaEtaInEndmax  = 0.007;
	double ElecDeltaPhiInBarmax  = 0.06;
	double ElecDeltaPhiInEndmax  = 0.06;
	double ElecDeltaPhiOutBarmax = 999.0;
	double ElecDeltaPhiOutEndmax = 999.0;
	double PhotHoverEBarmax      = 0.2;
	double PhotHoverEEndmax      = 0.2;
	double PhotSigmaEtaEtaBarmax = 0.011;
	double PhotSigmaEtaEtaEndmax = 0.025;
	double dR12min = 0.5;
	double dR21min = 0.5;
	double FracChmin = 0.1;
	double FracEmmin = 0.175;

	TString fChecklistFile = fOutputDir + "plots_Cleaning_checklist.txt";
	ofstream file;
	file.open(fChecklistFile, ios::out);
	TLine *line, *l1, *l2;
	double maxy;
	
	TCanvas *canv;
	TString tempstring1, tempstring2, outputdir;
	outputdir = fOutputDir + "Cleaning/Mu/";
	tempstring1 = "Mu DeltaP/P";
	tempstring2 = "MuDeltaPOverP";
	canv = new TCanvas(tempstring2, tempstring1 , 0, 0, 900, 700);
	fMuHistos[0]->DrawCopy();
	Util::PrintBoth(canv, tempstring2, outputdir);
	file << "* " << tempstring2 << endl;
	file << fAC->printAverage(tempstring2, fMuHistos[0]) << endl;

	tempstring1 = "Mu d0signif";
	tempstring2 = "Mud0signif";
	canv = new TCanvas(tempstring2, tempstring1 , 0, 0, 900, 700);
	maxy = fMuHistos[1]->GetMaximum();
	maxy = 1.05*maxy;
	fMuHistos[1]->SetMaximum(maxy);
	fMuHistos[1]->DrawCopy();
	line = new TLine(distVxmax,0,distVxmax,maxy);
	line->SetLineColor(kRed);
	line->SetLineWidth(2);
	line->Draw();
	Util::PrintBoth(canv, tempstring2, outputdir);
	file << "* " << tempstring2 << endl;
	file << fAC->printAverage(tempstring2, fMuHistos[1]) << endl;
	file << fAC->printRatio(tempstring2, fMuHistos[1], distVxmax, 100., 0., 100.) << endl;

	tempstring1 = "Mu dzsignif";
	tempstring2 = "Mudzsignif";
	canv = new TCanvas(tempstring2, tempstring1 , 0, 0, 900, 700);
	maxy = fMuHistos[2]->GetMaximum();
	maxy = 1.05*maxy;
	fMuHistos[2]->SetMaximum(maxy);
	fMuHistos[2]->DrawCopy();
	line = new TLine(distVxmax,0,distVxmax,maxy);
	line->SetLineColor(kRed);
	line->SetLineWidth(2);
	line->Draw();
	Util::PrintBoth(canv, tempstring2, outputdir);
	file << "* " << tempstring2 << endl;
	file << fAC->printAverage(tempstring2, fMuHistos[2]) << endl;
	file << fAC->printRatio(tempstring2, fMuHistos[2], distVxmax, 100., 0., 100.) << endl;

	tempstring1 = "Mu DR SS";
	tempstring2 = "MuDRSS";
	canv = new TCanvas(tempstring2, tempstring1 , 0, 0, 900, 700);
	maxy = fMuHistos[3]->GetMaximum();
	maxy = 1.05*maxy;
	fMuHistos[3]->SetMaximum(maxy);
	fMuHistos[3]->DrawCopy();
	line = new TLine(dRSSmuonmax,0,dRSSmuonmax,maxy);
	line->SetLineColor(kRed);
	line->SetLineWidth(2);
	line->Draw();
	Util::PrintBoth(canv, tempstring2, outputdir);
	file << "* " << tempstring2 << endl;
	file << fAC->printAverage(tempstring2, fMuHistos[3]) << endl;

	outputdir = fOutputDir + "Cleaning/El/";
	tempstring1 = "El H/E Bar";
	tempstring2 = "ElHcalOverEcalBar";
	canv = new TCanvas(tempstring2, tempstring1 , 0, 0, 900, 700);
	maxy = fElHistos[0]->GetMaximum();
	maxy = 1.05*maxy;
	fElHistos[0]->SetMaximum(maxy);
	fElHistos[0]->DrawCopy();
	line = new TLine(ElecHoverEBarmax,0,ElecHoverEBarmax,maxy);
	line->SetLineColor(kRed);
	line->SetLineWidth(2);
	line->Draw();
	Util::PrintBoth(canv, tempstring2, outputdir);
	file << "* " << tempstring2 << endl;
	file << fAC->printAverage(tempstring2, fElHistos[0]) << endl;
	file << fAC->printRatio(tempstring2, fElHistos[0], ElecHoverEBarmax, 100., 0., 100.) << endl;

	tempstring1 = "El SigmaIetaIeta Bar";
	tempstring2 = "ElSigmaIetaIetaBar";
	canv = new TCanvas(tempstring2, tempstring1 , 0, 0, 900, 700);
	maxy = fElHistos[1]->GetMaximum();
	maxy = 1.05*maxy;
	fElHistos[1]->SetMaximum(maxy);
	fElHistos[1]->DrawCopy();
	line = new TLine(ElecSigmaEtaEtaBarmax,0,ElecSigmaEtaEtaBarmax,maxy);
	line->SetLineColor(kRed);
	line->SetLineWidth(2);
	line->Draw();
	Util::PrintBoth(canv, tempstring2, outputdir);
	file << "* " << tempstring2 << endl;
	file << fAC->printAverage(tempstring2, fElHistos[1]) << endl;
	file << fAC->printRatio(tempstring2, fElHistos[1], ElecSigmaEtaEtaBarmax, 100., 0., 100.) << endl;

	tempstring1 = "El DeltaPhiSeedClusterAtCalo Bar";
	tempstring2 = "ElDeltaPhiSeedClusterAtCaloBar";
	canv = new TCanvas(tempstring2, tempstring1 , 0, 0, 900, 700);
	maxy = fElHistos[2]->GetMaximum();
	maxy = 1.05*maxy;
	fElHistos[2]->SetMaximum(maxy);
	fElHistos[2]->DrawCopy();
	l1 = new TLine(ElecDeltaPhiOutBarmax,0,ElecDeltaPhiOutBarmax,maxy);
	l1->SetLineColor(kRed);
	l1->SetLineWidth(2);
	l1->Draw();
	l2 = new TLine(-ElecDeltaPhiOutBarmax,0,-ElecDeltaPhiOutBarmax,maxy);
	l2->SetLineColor(kRed);
	l2->SetLineWidth(2);
	l2->Draw();
	Util::PrintBoth(canv, tempstring2, outputdir);
	file << "* " << tempstring2 << endl;
	file << fAC->printAverage(tempstring2, fElHistos[2]) << endl;
	file << fAC->printRatio(tempstring2, fElHistos[2], -ElecDeltaPhiOutBarmax, ElecDeltaPhiOutBarmax, -100., 100.) << endl;

	tempstring1 = "El DeltaEtaSeedClusterAtCalo Bar";
	tempstring2 = "ElDeltaEtaSeedClusterAtCaloBar";
	canv = new TCanvas(tempstring2, tempstring1 , 0, 0, 900, 700);
	fElHistos[3]->DrawCopy();
	Util::PrintBoth(canv, tempstring2, outputdir);
	file << "* " << tempstring2 << endl;
	file << fAC->printAverage(tempstring2, fElHistos[3]) << endl;

	tempstring1 = "El DeltaPhiSuperClusterAtVtx Bar";
	tempstring2 = "ElDeltaPhiSuperClusterAtVtxBar";
	canv = new TCanvas(tempstring2, tempstring1 , 0, 0, 900, 700);
	maxy = fElHistos[4]->GetMaximum();
	maxy = 1.05*maxy;
	fElHistos[4]->SetMaximum(maxy);
	fElHistos[4]->DrawCopy();
	l1 = new TLine(ElecDeltaPhiInBarmax,0,ElecDeltaPhiInBarmax,maxy);
	l1->SetLineColor(kRed);
	l1->SetLineWidth(2);
	l1->Draw();
	l2 = new TLine(-ElecDeltaPhiInBarmax,0,-ElecDeltaPhiInBarmax,maxy);
	l2->SetLineColor(kRed);
	l2->SetLineWidth(2);
	l2->Draw();
	Util::PrintBoth(canv, tempstring2, outputdir);
	file << "* " << tempstring2 << endl;
	file << fAC->printAverage(tempstring2, fElHistos[4]) << endl;
	file << fAC->printRatio(tempstring2, fElHistos[4], -ElecDeltaPhiInBarmax, ElecDeltaPhiInBarmax, -100., 100.) << endl;

	tempstring1 = "El DeltaEtaSuperClusterAtVtx Bar";
	tempstring2 = "ElDeltaEtaSuperClusterAtVtxBar";
	canv = new TCanvas(tempstring2, tempstring1 , 0, 0, 900, 700);
	maxy = fElHistos[5]->GetMaximum();
	maxy = 1.05*maxy;
	fElHistos[5]->SetMaximum(maxy);
	fElHistos[5]->DrawCopy();
	l1 = new TLine(ElecDeltaEtaInBarmax,0,ElecDeltaEtaInBarmax,maxy);
	l1->SetLineColor(kRed);
	l1->SetLineWidth(2);
	l1->Draw();
	l2 = new TLine(-ElecDeltaEtaInBarmax,0,-ElecDeltaEtaInBarmax,maxy);
	l2->SetLineColor(kRed);
	l2->SetLineWidth(2);
	l2->Draw();
	Util::PrintBoth(canv, tempstring2, outputdir);
	file << "* " << tempstring2 << endl;
	file << fAC->printAverage(tempstring2, fElHistos[5]) << endl;
	file << fAC->printRatio(tempstring2, fElHistos[5], -ElecDeltaEtaInBarmax, ElecDeltaEtaInBarmax, -100., 100.) << endl;

	tempstring1 = "El ESuperClusterOverP Bar";
	tempstring2 = "ElESuperClusterOverPBar";
	canv = new TCanvas(tempstring2, tempstring1 , 0, 0, 900, 700);
	maxy = fElHistos[6]->GetMaximum();
	maxy = 1.05*maxy;
	fElHistos[6]->SetMaximum(maxy);
	fElHistos[6]->DrawCopy();
	line = new TLine(ElecEoverPInBarmin,0,ElecEoverPInBarmin,maxy);
	line->SetLineColor(kRed);
	line->SetLineWidth(2);
	line->Draw();
	Util::PrintBoth(canv, tempstring2, outputdir);
	file << "* " << tempstring2 << endl;
	file << fAC->printAverage(tempstring2, fElHistos[6]) << endl;
	file << fAC->printRatio(tempstring2, fElHistos[6], 0., ElecEoverPInBarmin, 0., 100.) << endl;

	tempstring1 = "El H/E End";
	tempstring2 = "ElHcalOverEcalEnd";
	canv = new TCanvas(tempstring2, tempstring1 , 0, 0, 900, 700);
	maxy = fElHistos[7]->GetMaximum();
	maxy = 1.05*maxy;
	fElHistos[7]->SetMaximum(maxy);
	fElHistos[7]->DrawCopy();
	line = new TLine(ElecHoverEEndmax,0,ElecHoverEEndmax,maxy);
	line->SetLineColor(kRed);
	line->SetLineWidth(2);
	line->Draw();
	Util::PrintBoth(canv, tempstring2, outputdir);
	file << "* " << tempstring2 << endl;
	file << fAC->printAverage(tempstring2, fElHistos[7]) << endl;
	file << fAC->printRatio(tempstring2, fElHistos[7], ElecHoverEEndmax, 100., 0., 100.) << endl;

	tempstring1 = "El SigmaIetaIeta End";
	tempstring2 = "ElSigmaIetaIetaEnd";
	canv = new TCanvas(tempstring2, tempstring1 , 0, 0, 900, 700);
	maxy = fElHistos[8]->GetMaximum();
	maxy = 1.05*maxy;
	fElHistos[8]->SetMaximum(maxy);
	fElHistos[8]->DrawCopy();
	line = new TLine(ElecSigmaEtaEtaEndmax,0,ElecSigmaEtaEtaEndmax,maxy);
	line->SetLineColor(kRed);
	line->SetLineWidth(2);
	line->Draw();
	Util::PrintBoth(canv, tempstring2, outputdir);
	file << "* " << tempstring2 << endl;
	file << fAC->printAverage(tempstring2, fElHistos[8]) << endl;
	file << fAC->printRatio(tempstring2, fElHistos[8], ElecSigmaEtaEtaEndmax, 100., 0., 100.) << endl;

	tempstring1 = "El DeltaPhiSeedClusterAtCalo End";
	tempstring2 = "ElDeltaPhiSeedClusterAtCaloEnd";
	canv = new TCanvas(tempstring2, tempstring1 , 0, 0, 900, 700);
	maxy = fElHistos[9]->GetMaximum();
	maxy = 1.05*maxy;
	fElHistos[9]->SetMaximum(maxy);
	fElHistos[9]->DrawCopy();
	l1 = new TLine(ElecDeltaPhiOutEndmax,0,ElecDeltaPhiOutEndmax,maxy);
	l1->SetLineColor(kRed);
	l1->SetLineWidth(2);
	l1->Draw();
	l2 = new TLine(-ElecDeltaPhiOutEndmax,0,-ElecDeltaPhiOutEndmax,maxy);
	l2->SetLineColor(kRed);
	l2->SetLineWidth(2);
	l2->Draw();
	Util::PrintBoth(canv, tempstring2, outputdir);
	file << "* " << tempstring2 << endl;
	file << fAC->printAverage(tempstring2, fElHistos[9]) << endl;
	file << fAC->printRatio(tempstring2, fElHistos[9], -ElecDeltaPhiOutEndmax, ElecDeltaPhiOutEndmax, -100., 100.) << endl;

	tempstring1 = "El DeltaEtaSeedClusterAtCalo End";
	tempstring2 = "ElDeltaEtaSeedClusterAtCaloEnd";
	canv = new TCanvas(tempstring2, tempstring1 , 0, 0, 900, 700);
	fElHistos[10]->DrawCopy();
	Util::PrintBoth(canv, tempstring2, outputdir);
	file << "* " << tempstring2 << endl;
	file << fAC->printAverage(tempstring2, fElHistos[10]) << endl;

	tempstring1 = "El DeltaPhiSuperClusterAtVtx End";
	tempstring2 = "ElDeltaPhiSuperClusterAtVtxEnd";
	canv = new TCanvas(tempstring2, tempstring1 , 0, 0, 900, 700);
	maxy = fElHistos[11]->GetMaximum();
	maxy = 1.05*maxy;
	fElHistos[11]->SetMaximum(maxy);
	fElHistos[11]->DrawCopy();
	l1 = new TLine(ElecDeltaPhiInEndmax,0,ElecDeltaPhiInEndmax,maxy);
	l1->SetLineColor(kRed);
	l1->SetLineWidth(2);
	l1->Draw();
	l2 = new TLine(-ElecDeltaPhiInEndmax,0,-ElecDeltaPhiInEndmax,maxy);
	l2->SetLineColor(kRed);
	l2->SetLineWidth(2);
	l2->Draw();
	Util::PrintBoth(canv, tempstring2, outputdir);
	file << "* " << tempstring2 << endl;
	file << fAC->printAverage(tempstring2, fElHistos[11]) << endl;
	file << fAC->printRatio(tempstring2, fElHistos[11], -ElecDeltaPhiInEndmax, ElecDeltaPhiInEndmax, -100., 100.) << endl;

	tempstring1 = "El DeltaEtaSuperClusterAtVtx End";
	tempstring2 = "ElDeltaEtaSuperClusterAtVtxEnd";
	canv = new TCanvas(tempstring2, tempstring1 , 0, 0, 900, 700);
	maxy = fElHistos[12]->GetMaximum();
	maxy = 1.05*maxy;
	fElHistos[12]->SetMaximum(maxy);
	fElHistos[12]->DrawCopy();
	l1 = new TLine(ElecDeltaEtaInEndmax,0,ElecDeltaEtaInEndmax,maxy);
	l1->SetLineColor(kRed);
	l1->SetLineWidth(2);
	l1->Draw();
	l2 = new TLine(-ElecDeltaEtaInEndmax,0,-ElecDeltaEtaInEndmax,maxy);
	l2->SetLineColor(kRed);
	l2->SetLineWidth(2);
	l2->Draw();
	Util::PrintBoth(canv, tempstring2, outputdir);
	file << "* " << tempstring2 << endl;
	file << fAC->printAverage(tempstring2, fElHistos[12]) << endl;
	file << fAC->printRatio(tempstring2, fElHistos[12], -ElecDeltaEtaInEndmax, ElecDeltaEtaInEndmax, -100., 100.) << endl;

	tempstring1 = "El ESuperClusterOverP End";
	tempstring2 = "ElESuperClusterOverPEnd";
	canv = new TCanvas(tempstring2, tempstring1 , 0, 0, 900, 700);
	maxy = fElHistos[13]->GetMaximum();
	maxy = 1.05*maxy;
	fElHistos[13]->SetMaximum(maxy);
	fElHistos[13]->DrawCopy();
	line = new TLine(ElecEoverPInEndmin,0,ElecEoverPInEndmin,maxy);
	line->SetLineColor(kRed);
	line->SetLineWidth(2);
	line->Draw();
	Util::PrintBoth(canv, tempstring2, outputdir);
	file << "* " << tempstring2 << endl;
	file << fAC->printAverage(tempstring2, fElHistos[13]) << endl;
	file << fAC->printRatio(tempstring2, fElHistos[13], 0., ElecEoverPInEndmin, 0., 100.) << endl;
	
	tempstring1 = "El d0signif";
	tempstring2 = "Eld0signif";
	canv = new TCanvas(tempstring2, tempstring1 , 0, 0, 900, 700);
	maxy = fElHistos[14]->GetMaximum();
	maxy = 1.05*maxy;
	fElHistos[14]->SetMaximum(maxy);
	fElHistos[14]->DrawCopy();
	line = new TLine(distVxmax,0,distVxmax,maxy);
	line->SetLineColor(kRed);
	line->SetLineWidth(2);
	line->Draw();
	Util::PrintBoth(canv, tempstring2, outputdir);
	file << "* " << tempstring2 << endl;
	file << fAC->printAverage(tempstring2, fElHistos[14]) << endl;
	file << fAC->printRatio(tempstring2, fElHistos[14], distVxmax, 100., 0., 100.) << endl;
	
	tempstring1 = "El dzsignif";
	tempstring2 = "Eldzsignif";
	canv = new TCanvas(tempstring2, tempstring1 , 0, 0, 900, 700);
	maxy = fElHistos[15]->GetMaximum();
	maxy = 1.05*maxy;
	fElHistos[15]->SetMaximum(maxy);
	fElHistos[15]->DrawCopy();
	line = new TLine(distVxmax,0,distVxmax,maxy);
	line->SetLineColor(kRed);
	line->SetLineWidth(2);
	line->Draw();
	Util::PrintBoth(canv, tempstring2, outputdir);
	file << "* " << tempstring2 << endl;
	file << fAC->printAverage(tempstring2, fElHistos[15]) << endl;
	file << fAC->printRatio(tempstring2, fElHistos[15], distVxmax, 100., 0., 100.) << endl;
	
	tempstring1 = "El DR SS";
	tempstring2 = "ElDRSS";
	canv = new TCanvas(tempstring2, tempstring1 , 0, 0, 900, 700);
	fElHistos[16]->DrawCopy();
	Util::PrintBoth(canv, tempstring2, outputdir);
	file << "* " << tempstring2 << endl;
	file << fAC->printAverage(tempstring2, fElHistos[16]) << endl;
	
	tempstring1 = "El DR OS";
	tempstring2 = "ElDROS";
	canv = new TCanvas(tempstring2, tempstring1 , 0, 0, 900, 700);
	fElHistos[17]->DrawCopy();
	Util::PrintBoth(canv, tempstring2, outputdir);
	file << "* " << tempstring2 << endl;
	file << fAC->printAverage(tempstring2, fElHistos[17]) << endl;

	tempstring1 = "El InvM Same Sign";
	tempstring2 = "ElInvMSS";
	canv = new TCanvas(tempstring2, tempstring1 , 0, 0, 900, 700);
	fElHistos[18]->DrawCopy();
	Util::PrintBoth(canv, tempstring2, outputdir);
	file << "* " << tempstring2 << endl;
	file << fAC->printAverage(tempstring2, fElHistos[18]) << endl;

	tempstring1 = "El InvM Opp Sign";
	tempstring2 = "ElInvMOS";
	canv = new TCanvas(tempstring2, tempstring1 , 0, 0, 900, 700);
	fElHistos[19]->DrawCopy();
	Util::PrintBoth(canv, tempstring2, outputdir);
	file << "* " << tempstring2 << endl;
	file << fAC->printAverage(tempstring2, fElHistos[19]) << endl;

	outputdir = fOutputDir + "Cleaning/Ph/";
	tempstring1 = "Phot H/E Bar";
	tempstring2 = "PhoHoverEBar";
	canv = new TCanvas(tempstring2, tempstring1 , 0, 0, 900, 700);
	maxy = fPhHistos[0]->GetMaximum();
	maxy = 1.05*maxy;
	fPhHistos[0]->SetMaximum(maxy);
	fPhHistos[0]->DrawCopy();
	line = new TLine(PhotHoverEBarmax,0,PhotHoverEBarmax,maxy);
	line->SetLineColor(kRed);
	line->SetLineWidth(2);
	line->Draw();
	Util::PrintBoth(canv, tempstring2, outputdir);
	file << "* " << tempstring2 << endl;
	file << fAC->printAverage(tempstring2, fPhHistos[0]) << endl;
	file << fAC->printRatio(tempstring2, fPhHistos[0], PhotHoverEBarmax, 100., 0., 100.) << endl;

	tempstring1 = "Phot sigmaIetaIeta Bar";
	tempstring2 = "PhoSigmaIetaIetaBar";
	canv = new TCanvas(tempstring2, tempstring1 , 0, 0, 900, 700);
	maxy = fPhHistos[1]->GetMaximum();
	maxy = 1.05*maxy;
	fPhHistos[1]->SetMaximum(maxy);
	fPhHistos[1]->DrawCopy();
	line = new TLine(PhotSigmaEtaEtaBarmax,0,PhotSigmaEtaEtaBarmax,maxy);
	line->SetLineColor(kRed);
	line->SetLineWidth(2);
	line->Draw();
	Util::PrintBoth(canv, tempstring2, outputdir);
	file << "* " << tempstring2 << endl;
	file << fAC->printAverage(tempstring2, fPhHistos[1]) << endl;
	file << fAC->printRatio(tempstring2, fPhHistos[1], PhotSigmaEtaEtaBarmax, 100., 0., 100.) << endl;

	tempstring1 = "Phot H/E End";
	tempstring2 = "PhoHoverEEnd";
	canv = new TCanvas(tempstring2, tempstring1 , 0, 0, 900, 700);
	maxy = fPhHistos[2]->GetMaximum();
	maxy = 1.05*maxy;
	fPhHistos[2]->SetMaximum(maxy);
	fPhHistos[2]->DrawCopy();
	line = new TLine(PhotHoverEEndmax,0,PhotHoverEEndmax,maxy);
	line->SetLineColor(kRed);
	line->SetLineWidth(2);
	line->Draw();
	Util::PrintBoth(canv, tempstring2, outputdir);
	file << "* " << tempstring2 << endl;
	file << fAC->printAverage(tempstring2, fPhHistos[2]) << endl;
	file << fAC->printRatio(tempstring2, fPhHistos[2], PhotHoverEEndmax, 100., 0., 100.) << endl;

	tempstring1 = "Phot sigmaIetaIeta End";
	tempstring2 = "PhoSigmaIetaIetaEnd";
	canv = new TCanvas(tempstring2, tempstring1 , 0, 0, 900, 700);
	maxy = fPhHistos[3]->GetMaximum();
	maxy = 1.05*maxy;
	fPhHistos[3]->SetMaximum(maxy);
	fPhHistos[3]->DrawCopy();
	line = new TLine(PhotSigmaEtaEtaEndmax,0,PhotSigmaEtaEtaEndmax,maxy);
	line->SetLineColor(kRed);
	line->SetLineWidth(2);
	line->Draw();
	Util::PrintBoth(canv, tempstring2, outputdir);
	file << "* " << tempstring2 << endl;
	file << fAC->printAverage(tempstring2, fPhHistos[3]) << endl;
	file << fAC->printRatio(tempstring2, fPhHistos[3], PhotSigmaEtaEtaEndmax, 100., 0., 100.) << endl;

	tempstring1 = "Phot DR";
	tempstring2 = "PhDR";
	canv = new TCanvas(tempstring2, tempstring1 , 0, 0, 900, 700);
	fPhHistos[4]->DrawCopy();
	Util::PrintBoth(canv, tempstring2, outputdir);
	file << "* " << tempstring2 << endl;
	file << fAC->printAverage(tempstring2, fPhHistos[4]) << endl;

	tempstring1 = "Phot InvM";
	tempstring2 = "PhInvM";
	canv = new TCanvas(tempstring2, tempstring1 , 0, 0, 900, 700);
	fPhHistos[5]->DrawCopy();
	Util::PrintBoth(canv, tempstring2, outputdir);
	file << "* " << tempstring2 << endl;
	file << fAC->printAverage(tempstring2, fPhHistos[5]) << endl;
	
	outputdir = fOutputDir + "Cleaning/Jet/";
	tempstring1 = "Jet d0 PV";
	tempstring2 = "Jetd0PV";
	canv = new TCanvas(tempstring2, tempstring1 , 0, 0, 900, 700);
	fJHistos[0]->DrawCopy();
	Util::PrintBoth(canv, tempstring2, outputdir);
	file << "* " << tempstring2 << endl;
	file << fAC->printAverage(tempstring2, fJHistos[0]) << endl;
	file << fAC->printRatio(tempstring2, fJHistos[0], distVxmax, 100., 0., 100.) << endl;
	
	tempstring1 = "Jet dz PV";
	tempstring2 = "JetdzPV";
	canv = new TCanvas(tempstring2, tempstring1 , 0, 0, 900, 700);
	fJHistos[1]->DrawCopy();
	Util::PrintBoth(canv, tempstring2, outputdir);
	file << "* " << tempstring2 << endl;
	file << fAC->printAverage(tempstring2, fJHistos[1]) << endl;
	file << fAC->printRatio(tempstring2, fJHistos[1], distVxmax, 100., 0., 100.) << endl;
	
	tempstring1 = "Jet d0signif";
	tempstring2 = "Jetd0signif";
	canv = new TCanvas(tempstring2, tempstring1 , 0, 0, 900, 700);
	maxy = fJHistos[2]->GetMaximum();
	maxy = 1.05*maxy;
	fJHistos[2]->SetMaximum(maxy);
	fJHistos[2]->DrawCopy();
	line = new TLine(distVxmax,0,distVxmax,maxy);
	line->SetLineColor(kRed);
	line->SetLineWidth(2);
	line->Draw();
	Util::PrintBoth(canv, tempstring2, outputdir);
	file << "* " << tempstring2 << endl;
	file << fAC->printAverage(tempstring2, fJHistos[2]) << endl;
	
	tempstring1 = "Jet dzsignif";
	tempstring2 = "Jetdzsignif";
	canv = new TCanvas(tempstring2, tempstring1 , 0, 0, 900, 700);
	maxy = fJHistos[3]->GetMaximum();
	maxy = 1.05*maxy;
	fJHistos[3]->SetMaximum(maxy);
	fJHistos[3]->DrawCopy();
	line = new TLine(distVxmax,0,distVxmax,maxy);
	line->SetLineColor(kRed);
	line->SetLineWidth(2);
	line->Draw();
	Util::PrintBoth(canv, tempstring2, outputdir);
	file << "* " << tempstring2 << endl;
	file << fAC->printAverage(tempstring2, fJHistos[3]) << endl;

	outputdir = fOutputDir + "Cleaning/MET/";
	tempstring1 = "Dphi(MET,j1)";
	tempstring2 = "DphiMETJet1";
	canv = new TCanvas(tempstring2, tempstring1 , 0, 0, 900, 700);
	fMETHistos[0]->DrawCopy();
	Util::PrintBoth(canv, tempstring2, outputdir);
	file << "* " << tempstring2 << endl;
	file << fAC->printAverage(tempstring2, fMETHistos[0]) << endl;

	tempstring1 = "Dphi(MET,j2)";
	tempstring2 = "DphiMETJet2";
	canv = new TCanvas(tempstring2, tempstring1 , 0, 0, 900, 700);
	fMETHistos[1]->DrawCopy();
	Util::PrintBoth(canv, tempstring2, outputdir);
	file << "* " << tempstring2 << endl;
	file << fAC->printAverage(tempstring2, fMETHistos[1]) << endl;

	tempstring1 = "MET R12";
	tempstring2 = "METR12";
	canv = new TCanvas(tempstring2, tempstring1 , 0, 0, 900, 700);
	maxy = fMETHistos[2]->GetMaximum();
	maxy = 1.05*maxy;
	fMETHistos[2]->SetMaximum(maxy);
	fMETHistos[2]->DrawCopy();
	line = new TLine(dR12min,0,dR12min,maxy);
	line->SetLineColor(kRed);
	line->SetLineWidth(2);
	line->Draw();
	Util::PrintBoth(canv, tempstring2, outputdir);
	file << "* " << tempstring2 << endl;
	file << fAC->printAverage(tempstring2, fMETHistos[2]) << endl;

	tempstring1 = "MET R21";
	tempstring2 = "METR21";
	canv = new TCanvas(tempstring2, tempstring1 , 0, 0, 900, 700);
	maxy = fMETHistos[3]->GetMaximum();
	maxy = 1.05*maxy;
	fMETHistos[3]->SetMaximum(maxy);
	fMETHistos[3]->DrawCopy();
	line = new TLine(dR21min,0,dR21min,maxy);
	line->SetLineColor(kRed);
	line->SetLineWidth(2);
	line->Draw();
	Util::PrintBoth(canv, tempstring2, outputdir);
	file << "* " << tempstring2 << endl;
	file << fAC->printAverage(tempstring2, fMETHistos[3]) << endl;

	tempstring1 = "MET Rsum";
	tempstring2 = "METRsum";
	canv = new TCanvas(tempstring2, tempstring1 , 0, 0, 900, 700);
	fMETHistos[4]->DrawCopy();
	Util::PrintBoth(canv, tempstring2, outputdir);
	file << "* " << tempstring2 << endl;
	file << fAC->printAverage(tempstring2, fMETHistos[4]) << endl;

	tempstring1 = "MET Dphi 12";
	tempstring2 = "METDphi12";
	canv = new TCanvas(tempstring2, tempstring1 , 0, 0, 900, 700);
	if(fTR->GetEntries() < 100000) fMETDphi12->SetMarkerStyle(6);
	fMETDphi12->DrawCopy();
	Util::PrintBoth(canv, tempstring2, outputdir);

	tempstring1 = "Evt Em Frac";
	tempstring2 = "EvtEmFrac";
	canv = new TCanvas(tempstring2, tempstring1 , 0, 0, 900, 700);
	maxy = fMETHistos[5]->GetMaximum();
	maxy = 1.05*maxy;
	fMETHistos[5]->SetMaximum(maxy);
	fMETHistos[5]->DrawCopy();
	line = new TLine(FracEmmin,0,FracEmmin,maxy);
	line->SetLineColor(kRed);
	line->SetLineWidth(2);
	line->Draw();
	Util::PrintBoth(canv, tempstring2, outputdir);
	file << "* " << tempstring2 << endl;
	file << fAC->printAverage(tempstring2, fMETHistos[5]) << endl;

	tempstring1 = "Evt Ch Frac";
	tempstring2 = "EvtChFrac";
	canv = new TCanvas(tempstring2, tempstring1 , 0, 0, 900, 700);
	maxy = fMETHistos[6]->GetMaximum();
	maxy = 1.05*maxy;
	fMETHistos[6]->SetMaximum(maxy);
	fMETHistos[6]->DrawCopy();
	line = new TLine(FracChmin,0,FracChmin,maxy);
	line->SetLineColor(kRed);
	line->SetLineWidth(2);
	line->Draw();
	Util::PrintBoth(canv, tempstring2, outputdir);
	file << "* " << tempstring2 << endl;
	file << fAC->printAverage(tempstring2, fMETHistos[6]) << endl;
	
	file.close();
	
	// write out the cleaning statistics
	TString cleanerstatfile = fOutputDir + "cleanerStats.txt";
	fTC->StatWrite(cleanerstatfile);

	// write the info file
	const int nentries = fTR->fChain->GetEntries();
	PrintInfoStart(nentries);

}

void PhysQCAnalysis::PrintInfoStart(int nEntries){
// output the text file with run info

	ofstream file;
	file.open(fOutputDir + "/info.txt");
	time_t tim;
	time(&tim);
//	cout << ctime(&tim) << endl;
	file << ctime(&tim) << endl;
	file << nEntries << endl;

	file.close();
	cout << " Info file made in " << fOutputDir << endl;

}

void PhysQCAnalysis::GetEvtEmChFrac(double & fracEm, double & fracCh){
// Computes the event EM and Charged fractions

	double pt_track = 0.;
	double et_em = 0.;
	double et_had = 0.;
	for( int i = 0; i < fTR->NMus; ++i ){
		if(fTR->MuGood[i] != 0) continue;
		pt_track += fTR->MuPt[i];
	}
	for( int i = 0; i < fTR->NEles; ++i ){
		if(fTR->ElGood[i] != 0) continue;
		pt_track += fTR->ElPt[i];
		et_em += fTR->ElEt[i];
	}
	for( int i = 0; i < fTR->NPhotons; ++i ){
		if(fTR->PhoGood[i] != 0) continue;
		et_em += fTR->PhoPt[i];
	}
	for( int i = 0; i < fTR->NJets; ++i ){
		if(fTR->JGood[i] != 0) continue;
		pt_track += fTR->JChfrac[i] * fTR->JPt[i];
		et_em    += fTR->JEMfrac[i] * fTR->JEt[i];
		et_had   += (1.-fTR->JEMfrac[i]) * fTR->JEt[i];
	}

	fracCh = 0.;
	fracEm = 0.;
	if( et_em + et_had <= 0. ){
		if( fTR->NMus < 1 ) return;
		fracCh = 1.;
		fracEm = 1.;
	} else {
		fracCh = pt_track / (et_em + et_had);
		fracEm = et_em / (et_em + et_had);
	}
	return;

}

void PhysQCAnalysis::MakePlots(TString plotlist, TTree *tree, TCut select){
	fAC->plotPlotList(plotlist, tree, "", select);
}

void PhysQCAnalysis::MakeElIDPlots(TTree *tree, TCut select){
	fAC->plotEID(select, tree, "");
}

void PhysQCAnalysis::PlotTriggerStats(){
	const int nentries = fTR->fChain->GetEntries();

        TH1I* hlt_stats,* l1p_stats,* l1t_stats;
        TString hltPath("analyze/HLTTriggerStats");
        TString l1pPath("analyze/L1PhysTriggerStats");
        TString l1tPath("analyze/L1TechTriggerStats");
        TString subdir("TriggerStats");

        // Loop over all files (if chain) and add histograms
        if ( !fTR->isChain() ) {
          TFile *f = fTR->fChain->GetCurrentFile();
          hlt_stats = (TH1I*)f->Get(hltPath);
          l1p_stats = (TH1I*)f->Get(l1pPath);
          l1t_stats = (TH1I*)f->Get(l1tPath);
        } else {
          TObjArray *fileElements = ((TChain*)fTR->fChain)->GetListOfFiles();
          TIter next(fileElements);
          TChainElement *chEl=0;
          bool firstFile = true;
          while (( chEl=(TChainElement*)next() )) {
            TFile* f = TFile::Open(chEl->GetTitle());
            if ( f == NULL ) continue;
            if ( firstFile ) {
              firstFile = false;
              hlt_stats = (TH1I*)f->Get(hltPath)->Clone();
              l1p_stats = (TH1I*)f->Get(l1pPath)->Clone();
              l1t_stats = (TH1I*)f->Get(l1tPath)->Clone();
              hlt_stats->SetDirectory(0);
              l1p_stats->SetDirectory(0);
              l1t_stats->SetDirectory(0);
            } else {
              hlt_stats->Add( (TH1I*)f->Get(hltPath) );
              l1p_stats->Add( (TH1I*)f->Get(l1pPath) );
              l1t_stats->Add( (TH1I*)f->Get(l1tPath) );
            }
          }
        }

        // Set style
	hlt_stats->GetXaxis()->LabelsOption("v");
	hlt_stats->GetXaxis()->SetLabelSize(0.035);
	hlt_stats->SetMaximum(nentries + 0.05*nentries);

	l1p_stats->GetXaxis()->LabelsOption("v");
	l1p_stats->GetXaxis()->SetLabelSize(0.033);
	l1p_stats->SetMaximum(nentries + 0.05*nentries);

	l1t_stats->SetMaximum(nentries + 0.05*nentries);

	fTlat->SetTextColor(kBlack);
	fTlat->SetNDC(kTRUE);
	fTlat->SetTextSize(0.04);
	gStyle->SetOptStat(0);

	TString entries = Form("# Entries = %d", nentries);
	TLine *l1;

	TCanvas *canv;
	TString tempstring;
	tempstring = "HLT Trigger Bits (0-50)";
	canv = new TCanvas("HLTStats1", tempstring , 0, 0, 900, 700);
	canv->SetBottomMargin(0.50);
	gPad->SetGridy();
	hlt_stats->GetXaxis()->SetRange(0,50);
	hlt_stats->DrawCopy();
	l1 = new TLine(0, nentries, 50, nentries);
	l1->SetLineColor(kRed);
	l1->SetLineWidth(2);
	l1->SetLineStyle(2);
	l1->Draw();
	fTlat->DrawLatex(0.17,0.92, tempstring);
	fTlat->DrawLatex(0.60,0.92, entries);
	Util::PrintBoth(canv, "HLTStats1", fOutputDir+subdir);

	tempstring = "HLT Trigger Bits (51-100)";
	canv = new TCanvas("HLTStats2", tempstring , 0, 0, 900, 700);
	canv->SetBottomMargin(0.50);
	gPad->SetGridy();
	hlt_stats->GetXaxis()->SetRange(51,100);
	hlt_stats->DrawCopy();
	l1 = new TLine(50, nentries, 100, nentries);
	l1->SetLineColor(kRed);
	l1->SetLineWidth(2);
	l1->SetLineStyle(2);
	l1->Draw();
	fTlat->DrawLatex(0.17,0.92, tempstring);
	fTlat->DrawLatex(0.60,0.92, entries);	
	Util::PrintBoth(canv, "HLTStats2", fOutputDir+subdir);

	tempstring = "L1 Phys Bits (0-63)";
	canv = new TCanvas("L1PStats1", tempstring , 0, 0, 900, 700);
	canv->SetBottomMargin(0.42);
	gPad->SetGridy();
	l1p_stats->GetXaxis()->SetRange(0,63);
	l1p_stats->DrawCopy();
	l1 = new TLine(0, nentries, 63, nentries);
	l1->SetLineColor(kRed);
	l1->SetLineWidth(2);
	l1->SetLineStyle(2);
	l1->Draw();
	fTlat->DrawLatex(0.17,0.92, tempstring);
	fTlat->DrawLatex(0.60,0.92, entries);
	Util::PrintBoth(canv, "L1PStats1", fOutputDir+subdir);

	tempstring = "L1 Phys Bits (64-128)";
	canv = new TCanvas("L1PStats2", tempstring , 0, 0, 900, 700);
	canv->SetBottomMargin(0.42);
	gPad->SetGridy();
	l1p_stats->GetXaxis()->SetRange(64,128);
	l1p_stats->DrawCopy();
	l1 = new TLine(63, nentries, 128, nentries);
	l1->SetLineColor(kRed);
	l1->SetLineWidth(2);
	l1->SetLineStyle(2);
	l1->Draw();
	fTlat->DrawLatex(0.17,0.92, tempstring);
	fTlat->DrawLatex(0.60,0.92, entries);
	Util::PrintBoth(canv, "L1PStats2", fOutputDir+subdir);

	tempstring = "L1 Tech Bits";
	canv = new TCanvas("L1TStats", tempstring , 0, 0, 900, 700);
	gPad->SetGridy();
	l1t_stats->GetXaxis()->SetRange(0,64);
	l1t_stats->DrawCopy();
	l1 = new TLine(0, nentries, 64, nentries);
	l1->SetLineColor(kRed);
	l1->SetLineWidth(2);
	l1->SetLineStyle(2);
	l1->Draw();
	fTlat->DrawLatex(0.17,0.92, tempstring);
	fTlat->DrawLatex(0.60,0.92, entries);
	Util::PrintBoth(canv, "L1TStats", fOutputDir+subdir);

        gROOT->cd(); // Leave local file
}

double PhysQCAnalysis::invMass(double p1[], double p2[]){
	double msq = (p1[3]+p2[3])*(p1[3]+p2[3]) - (p1[0]+p2[0])*(p1[0]+p2[0])
	           - (p1[1]+p2[1])*(p1[1]+p2[1]) - (p1[2]+p2[2])*(p1[2]+p2[2]);
	if (msq < 0.) msq = 0.;
	return sqrt(msq);
}

