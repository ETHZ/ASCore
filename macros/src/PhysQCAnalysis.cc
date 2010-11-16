#include <TChainElement.h>
#include <TLorentzVector.h>

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
}
PhysQCAnalysis::~PhysQCAnalysis(){
	delete fTlat;
	delete fAC;
}

void PhysQCAnalysis::Begin(const char* filename){
	fHistFile = new TFile(fOutputDir + TString(filename), "RECREATE");
	fAC->readVarNames("varnames.dat");
	fAC->setOutputDir(fOutputDir);
	fAC->setGlobalTag(fTag);
	BookHistos();
}


void PhysQCAnalysis::Analyze1(){
// Analysis before cleaning

	// Skip empty events
	if((fTR->NMus + fTR->NEles + fTR->NJets) <= 0 && fTR->NPhotons < 1) return;

	// Primary vertex plots
	fPvxHistos[0]->Fill(fTR->PrimVtxNChi2);
	double xVx = fTR->PrimVtxx - fTR->Beamspotx;
	double yVx = fTR->PrimVtxy - fTR->Beamspoty;
	double zVx = fTR->PrimVtxz - fTR->Beamspotz;
	double rVx = sqrt(xVx*xVx + yVx*yVx);
	fPvxHistos[1]->Fill(rVx);
	fPvxHistos[2]->Fill(zVx);
	fPvxHistos[3]->Fill(fTR->PrimVtxPtSum);

	TLorentzVector p1(0.,0.,0.,0.), p2(0.,0.,0.,0.), psum(0.,0.,0.,0.);
	double minv = -999.99;

	double drVxsq = fTR->PrimVtxxE*fTR->PrimVtxxE + fTR->PrimVtxyE*fTR->PrimVtxyE;

	// Muon plots
	for(size_t i = 0; i < fTR->NMus; ++i){
		if(fTR->MuIsGlobalMuon[i] == 0) continue;
		fMuHistos[0]->Fill(fTR->MuPtE[i]/fTR->MuPt[i]);
		double d0  = fTR->MuD0PV[i];
		double dd0 = sqrt(fTR->MuD0E[i]*fTR->MuD0E[i] + drVxsq);
		if (dd0 <= 0.) dd0 = 0.001;
		double dz  = fTR->MuDzPV[i];
		double ddz = sqrt(fTR->MuDzE[i]*fTR->MuDzE[i] + fTR->PrimVtxzE*fTR->PrimVtxzE);
		if (ddz <= 0.) ddz = 0.001;
		fMuHistos[1]->Fill(fabs(d0) / dd0);
		fMuHistos[2]->Fill(fabs(dz) / ddz);
		fMuHistos[3]->Fill(fTR->MuNChi2[i]);
		fMuHistos[4]->Fill(fTR->MuNTkHits[i]);
	}

	for(size_t i = 0; i < fTR->NMus; ++i){
		if(fTR->MuIsGlobalMuon[i] == 0) continue;
		for(size_t j = i+1; j < fTR->NMus; ++j){
			if(fTR->MuIsGlobalMuon[j] == 0) continue;
			double deltaR = Util::GetDeltaR(fTR->MuEta[i], fTR->MuEta[j], fTR->MuPhi[i], fTR->MuPhi[j]);
			p1.SetPxPyPzE(fTR->MuPx[i],fTR->MuPy[i],fTR->MuPz[i],fTR->MuE[i]);
			p2.SetPxPyPzE(fTR->MuPx[j],fTR->MuPy[j],fTR->MuPz[j],fTR->MuE[j]);
			psum = p1 + p2;
			minv = psum.M();
			if (fTR->MuCharge[i] == fTR->MuCharge[j]) {
				fMuHistos[5]->Fill(deltaR);
				if (deltaR < 0.5) {
					fMuHistos[7]->Fill(minv);
				}
			} else {
				fMuHistos[6]->Fill(deltaR);
				if (deltaR < 0.5) {
					fMuHistos[8]->Fill(minv);
				}
			}
		}
	}

	// Electron plots
	for(size_t i = 0; i < fTR->NEles; ++i){
		if (fabs(fTR->ElEta[i]) < 1.479) {
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
			p1.SetPxPyPzE(fTR->ElPx[i],fTR->ElPy[i],fTR->ElPz[i],fTR->ElE[i]);
			p2.SetPxPyPzE(fTR->ElPx[j],fTR->ElPy[j],fTR->ElPz[j],fTR->ElE[j]);
			psum = p1 + p2;
			minv = psum.M();
			if (fTR->ElCharge[i] == fTR->ElCharge[j]) {
				fElHistos[16]->Fill(deltaR);
				if (deltaR < 0.5) {
					fElHistos[18]->Fill(minv);
				}
			} else {
				fElHistos[17]->Fill(deltaR);
				if (deltaR < 0.5) {
					fElHistos[19]->Fill(minv);
				}
			}
		}
	}

	// Photon plots
	for(size_t i = 0; i < fTR->NPhotons; ++i){
		if (fabs(fTR->PhoEta[i]) < 1.479) {
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
				p1.SetPxPyPzE(fTR->PhoPx[i],fTR->PhoPy[i],fTR->PhoPz[i],fTR->PhoEnergy[i]);
				p2.SetPxPyPzE(fTR->PhoPx[j],fTR->PhoPy[j],fTR->PhoPz[j],fTR->PhoEnergy[j]);
				psum = p1 + p2;
				minv = psum.M();
				fPhHistos[5]->Fill(minv);
				fPhHistos[6]->Fill(minv);
			}
		}
	}

	// Jet plots
	for(size_t i = 0; i < fTR->NJets; ++i){
		double d0  = 0.;
		double dd0 = 0.001;
		double dz  = 0.;
		double ddz = 0.001;
		double dztry = fTR->JVtxz[i] - fTR->PrimVtxz;
		// Ignore jets with badly reconstructed vertices
		if (fTR->JVtxz[i] != -888.88 && fTR->JVtxz[i] != -777.77) {
			d0  = sqrt((fTR->JVtxx[i]-fTR->PrimVtxx)*(fTR->JVtxx[i]-fTR->PrimVtxx)
				+ (fTR->JVtxy[i]-fTR->PrimVtxy)*(fTR->JVtxy[i]-fTR->PrimVtxy) );
			dd0 = sqrt(fTR->JVtxExx[i] + fTR->JVtxEyy[i] + drVxsq);
			if (dd0 <= 0.) dd0 = 0.001;
			dz  = fTR->JVtxz[i] - fTR->PrimVtxz;
			ddz = sqrt(fTR->JVtxEzz[i] + fTR->PrimVtxzE*fTR->PrimVtxzE);
			if (ddz <= 0.) ddz = 0.001;
			fJHistos[0]->Fill(d0);
			fJHistos[1]->Fill(dz);
			fJHistos[2]->Fill(fabs(d0) / dd0);
			fJHistos[3]->Fill(fabs(dz) / ddz);
		}
		fJHistos[4]->Fill(fTR->JEMfrac[i]);
		fJHistos[5]->Fill(fTR->JChfrac[i]);
		fJHistos[6]->Fill(fTR->JID_n90Hits[i]);
		fJHistos[7]->Fill(fTR->JID_HPD[i]);
		fJChEMfrac->Fill(fTR->JEMfrac[i], fTR->JChfrac[i]);
		fJEMfracEta->Fill(fTR->JEta[i], fTR->JEMfrac[i]);
		fJChfracEta->Fill(fTR->JEta[i], fTR->JChfrac[i]);
	}

	// MET plots
	double METPhi = fTR->MuJESCorrMETphi;
	double MET = fTR->MuJESCorrMET;
	double METBadJetmin = 20.;
	if (fTR->NJets > 1 && MET > METBadJetmin) {
		double dPhiMJ1 = Util::DeltaPhi(fTR->JPhi[0], METPhi);
		double dPhiMJ2 = Util::DeltaPhi(fTR->JPhi[1], METPhi);
		double R12 = sqrt(dPhiMJ1*dPhiMJ1 + (TMath::Pi()-dPhiMJ2)*(TMath::Pi()-dPhiMJ2) );
		double R21 = sqrt(dPhiMJ2*dPhiMJ2 + (TMath::Pi()-dPhiMJ1)*(TMath::Pi()-dPhiMJ1) );
		double dPhij12 = Util::DeltaPhi(fTR->JPhi[0], fTR->JPhi[1]);
		double dRj12 = Util::GetDeltaR(fTR->JEta[0], fTR->JEta[1], fTR->JPhi[0], fTR->JPhi[1]);
		fMETHistos[0]->Fill(dPhiMJ1);
		fMETHistos[1]->Fill(dPhiMJ2);
		fMETHistos[2]->Fill(R12);
		fMETHistos[3]->Fill(R21);
		fMETHistos[4]->Fill(R12+R21);
		fMETDphi12->Fill(dPhiMJ1, dPhiMJ2);
		fMETR12R21->Fill(R12, R21);
		fMETR12Dphij12->Fill(R12, dPhij12);
		fMETR12Etaj12->Fill(R12, fTR->JEta[0]);
		fMETR12Etaj12->Fill(R12, fTR->JEta[1]);
		fMETR12j12PtRat->Fill(R12, fTR->JPt[1]/fTR->JPt[0]);
		fMETR12dRj12->Fill(R12, dRj12);
	}

	// Event plots
	double fracEm(0.), fracCh(0.);
	GetEvtEmChFrac(fracEm, fracCh);
	fEvtHistos[0]->Fill(fracEm);
	fEvtHistos[1]->Fill(fracCh);
}

void PhysQCAnalysis::Analyze2(){
// Analysis after cleaning

	// skip empty and bad events
	if ((fTR->NMus + fTR->NEles + fTR->NJets) <= 0 && fTR->NPhotons < 1) return;
	if (fTR->GoodEvent > 0) return;

	TLorentzVector p1(0.,0.,0.,0.), p2(0.,0.,0.,0.), psum(0.,0.,0.,0.);
	double minv = -999.99;

	fMuCIHistos[0]->Fill(fTR->NGMus);
	for(size_t i = 0; i < fTR->NMus; ++i){
		if(fTR->MuIsGlobalMuon[i] == 0) continue;
		if (fTR->MuIsIso[i] == 0) continue;
		fMuCIHistos[1]->Fill(fTR->MuPt[i]);
		fMuCIHistos[2]->Fill(fTR->MuEta[i]);
		fMuCIHistos[3]->Fill(fTR->MuPhi[i]);
		for(size_t j = i+1; j < fTR->NMus; ++j){
			if(fTR->MuIsGlobalMuon[j] == 0) continue;
			p1.SetPxPyPzE(fTR->MuPx[i],fTR->MuPy[i],fTR->MuPz[i],fTR->MuE[i]);
			p2.SetPxPyPzE(fTR->MuPx[j],fTR->MuPy[j],fTR->MuPz[j],fTR->MuE[j]);
			psum = p1 + p2;
			minv = psum.M();
			if (fTR->MuCharge[i] == fTR->MuCharge[j]) {
				fInvMHistos[0]->Fill(minv);
				fInvMHistos[1]->Fill(minv);
			} else {
				fInvMHistos[2]->Fill(minv);
				fInvMHistos[3]->Fill(minv);
			}
		}
	}

	fElCIHistos[0]->Fill(fTR->NEles);
	for(size_t i = 0; i < fTR->NEles; ++i){
		if (fTR->ElIsIso[i] == 0) continue;
		fElCIHistos[1]->Fill(fTR->ElPt[i]);
		fElCIHistos[2]->Fill(fTR->ElEta[i]);
		fElCIHistos[3]->Fill(fTR->ElPhi[i]);
		fElCIHistos[4]->Fill(fTR->ElSigmaIetaIeta[i]);
		for(size_t j = i+1; j < fTR->NEles; ++j){
			p1.SetPxPyPzE(fTR->ElPx[i],fTR->ElPy[i],fTR->ElPz[i],fTR->ElE[i]);
			p2.SetPxPyPzE(fTR->ElPx[j],fTR->ElPy[j],fTR->ElPz[j],fTR->ElE[j]);
			psum = p1 + p2;
			minv = psum.M();
			if (fTR->ElCharge[i] == fTR->ElCharge[j]) {
				fInvMHistos[4]->Fill(minv);
				fInvMHistos[5]->Fill(minv);
			} else {
				fInvMHistos[6]->Fill(minv);
				fInvMHistos[7]->Fill(minv);
			}
		}
	}

	fPhCIHistos[0]->Fill(fTR->NPhotons);
	for(size_t i = 0; i < fTR->NPhotons; ++i){
		if (fTR->PhoIsIso[i] == 0) continue;
		fPhCIHistos[1]->Fill(fTR->PhoPt[i]);
		fPhCIHistos[2]->Fill(fTR->PhoEta[i]);
		fPhCIHistos[3]->Fill(fTR->PhoPhi[i]);
		for(size_t j = i+1; j < fTR->NPhotons; ++j){
			p1.SetPxPyPzE(fTR->PhoPx[i],fTR->PhoPy[i],fTR->PhoPz[i],fTR->PhoEnergy[i]);
			p2.SetPxPyPzE(fTR->PhoPx[j],fTR->PhoPy[j],fTR->PhoPz[j],fTR->PhoEnergy[j]);
			psum = p1 + p2;
			minv = psum.M();
			fInvMHistos[8]->Fill(minv);
			fInvMHistos[9]->Fill(minv);
		}
	}

	fJCIHistos[0]->Fill(fTR->NJets);
	for(size_t i = 0; i < fTR->NJets; ++i){
		fJCIHistos[1]->Fill(fTR->JPt[i]);
		fJCIHistos[2]->Fill(fTR->JEta[i]);
		fJCIHistos[3]->Fill(fTR->JPhi[i]);
		fJCIHistos[4]->Fill(fTR->JEMfrac[i]);
		fJCIHistos[5]->Fill(fTR->JChfrac[i]);
		for(size_t j = i+1; j < fTR->NJets; ++j){
			p1.SetPxPyPzE(fTR->JPx[i],fTR->JPy[i],fTR->JPz[i],fTR->JE[i]);
			p2.SetPxPyPzE(fTR->JPx[j],fTR->JPy[j],fTR->JPz[j],fTR->JE[j]);
			psum = p1 + p2;
			minv = psum.M();
			fInvMHistos[10]->Fill(minv);
			if (fTR->NJets == 2) continue;
			fInvMHistos[11]->Fill(minv);
		}
		if (fTR->NJets == 1) {
			fJCIHistos[10]->Fill(fTR->JPt[i]);
			fJCIHistos[11]->Fill(fTR->JEta[i]);
			fJCIHistos[12]->Fill(fTR->JPhi[i]);
		}
	}

	double MET = fTR->MuCorrMET;
	double METBadJetmin = 20.;
	if (fTR->NJets > 1 && MET > METBadJetmin) {
		double METPhi = fTR->MuCorrMETphi;
		double dPhiMJ1 = Util::DeltaPhi(fTR->JPhi[0], METPhi);
		double dPhiMJ2 = Util::DeltaPhi(fTR->JPhi[1], METPhi);
		double R12 = sqrt(dPhiMJ1*dPhiMJ1 + (TMath::Pi()-dPhiMJ2)*(TMath::Pi()-dPhiMJ2) );
		double R21 = sqrt(dPhiMJ2*dPhiMJ2 + (TMath::Pi()-dPhiMJ1)*(TMath::Pi()-dPhiMJ1) );
		fJCIHistos[6]->Fill(dPhiMJ1);
		fJCIHistos[7]->Fill(dPhiMJ2);
		fJCIHistos[8]->Fill(R12);
		fJCIHistos[9]->Fill(R21);
	}
	if (fTR->NJets == 1 && MET > METBadJetmin) {
		double METPhi = fTR->MuCorrMETphi;
		double dPhiMJ1 = Util::DeltaPhi(fTR->JPhi[0], METPhi);
		fJCIHistos[13]->Fill(dPhiMJ1);
	}

}

void PhysQCAnalysis::End(){

	double chisqVxmax = fTC->fClean_chisqVxmax;
	double dRVxmax = fTC->fClean_dRVxmax;
	double dzVxmax = fTC->fClean_dzVxmax;
	double distVxmax = fTC->fClean_distVxmax;
	double PrimVtxPtSum = fTR->PrimVtxPtSum;
	double MuonDPbyPmax = fTC->fClean_MuonDPbyPmax;
	double MuonChi2max = fTC->fClean_MuonChi2max;
	double MuonNHitsmin = fTC->fClean_MuonNHitsmin;
	double dRSSmuonmax = fTC->fClean_dRSSmuonmax;
	double ElecHoverEBarmax      = fTC->fClean_ElecHoverEBarmax;
	double ElecHoverEEndmax      = fTC->fClean_ElecHoverEEndmax;
	double ElecSigmaEtaEtaBarmax = fTC->fClean_ElecSigmaEtaEtaBarmax;
	double ElecSigmaEtaEtaEndmax = fTC->fClean_ElecSigmaEtaEtaEndmax;
	double ElecEoverPInBarmin    = fTC->fClean_ElecEoverPInBarmin;
	double ElecEoverPInEndmin    = fTC->fClean_ElecEoverPInEndmin;
	double ElecDeltaEtaInBarmax  = fTC->fClean_ElecDeltaEtaInBarmax;
	double ElecDeltaEtaInEndmax  = fTC->fClean_ElecDeltaEtaInEndmax;
	double ElecDeltaPhiInBarmax  = fTC->fClean_ElecDeltaPhiInBarmax;
	double ElecDeltaPhiInEndmax  = fTC->fClean_ElecDeltaPhiInEndmax;
	double ElecDeltaPhiOutBarmax = fTC->fClean_ElecDeltaPhiOutBarmax;
	double ElecDeltaPhiOutEndmax = fTC->fClean_ElecDeltaPhiOutEndmax;
	double PhotHoverEBarmax      = fTC->fClean_PhotHoverEBarmax;
	double PhotHoverEEndmax      = fTC->fClean_PhotHoverEEndmax;
	double PhotSigmaEtaEtaBarmax = fTC->fClean_PhotSigmaEtaEtaBarmax;
	double PhotSigmaEtaEtaEndmax = fTC->fClean_PhotSigmaEtaEtaEndmax;
	double FracEmminJet      = fTC->fClean_FracEmminJet;
	double FracEmmaxJet      = fTC->fClean_FracEmmaxJet;
	double FracChminJet      = fTC->fClean_FracChminJet;
	double JID_n90Hitsmin    = fTC->fClean_JID_n90Hitsmin;
	double JID_HPDmax        = fTC->fClean_JID_HPDmax;
	double FracChmin = fTC->fClean_FracChmin;
	double FracEmmin = fTC->fClean_FracEmmin;
	double dPhiJetMETmin = fTC->fClean_dPhiJetMETmin;
	double dR12min = fTC->fClean_dR12min;
	double dR21min = fTC->fClean_dR21min;
	double MuPtmin = 10.;
	double MuEtamax = 2.4;
	double ElPtmin = 10.;
	double ElEtamax = 2.4;
	double PhPtmin = 10.;
	double PhEtamax = 2.4;
	double JPtmin = 30.;
	double JEtamax = 3.;

	TString fChecklistFile = fOutputDir + "plots_Cleaning_checklist.txt";
	ofstream file;
	file.open(fChecklistFile, ios::out);
	TLine *line, *l1, *l2;
	double maxy;
	TString tempstring1, tempstring2, outputdir;
	TCanvas * canv;

	// Primary vertices
	outputdir = fOutputDir + "Cleaning/PVx/";
	tempstring1 = "Prim Vtx NChi2";
	tempstring2 = "PrimVtxNChi2";
	PrintHisto(fPvxHistos[0], tempstring1, tempstring2, outputdir, chisqVxmax);
	file << "* " << tempstring2 << endl;
	file << fAC->printAverage(tempstring2, fPvxHistos[0]) << endl;

	tempstring1 = "Prim Vtx d0 BS";
	tempstring2 = "PrimVtxd0BS";
	PrintHisto(fPvxHistos[1], tempstring1, tempstring2, outputdir, dRVxmax);
	file << "* " << tempstring2 << endl;
	file << fAC->printAverage(tempstring2, fPvxHistos[1]) << endl;

	tempstring1 = "Prim Vtx dz BS";
	tempstring2 = "PrimVtxdzBS";
	PrintHisto(fPvxHistos[2], tempstring1, tempstring2, outputdir, dzVxmax);
	file << "* " << tempstring2 << endl;
	file << fAC->printAverage(tempstring2, fPvxHistos[2]) << endl;

	tempstring1 = "Prim Vtx PtSum";
	tempstring2 = "PrimVtxPtSum";
	PrintHisto(fPvxHistos[3], tempstring1, tempstring2, outputdir, PrimVtxPtSum, 999., 1, 1);
	file << "* " << tempstring2 << endl;
	file << fAC->printAverage(tempstring2, fPvxHistos[3]) << endl;
	file << fAC->printTailFraction(tempstring2, fPvxHistos[3], 0.05) << endl;
	file << fAC->printTailFraction(tempstring2, fPvxHistos[3], 0.01) << endl;

	// uncleaned muons
	outputdir = fOutputDir + "Cleaning/Mu/";
	tempstring1 = "Mu DeltaP/P";
	tempstring2 = "MuDeltaPOverP";
	PrintHisto(fMuHistos[0], tempstring1, tempstring2, outputdir, MuonDPbyPmax);
	file << "* " << tempstring2 << endl;
	file << fAC->printAverage(tempstring2, fMuHistos[0]) << endl;

	tempstring1 = "Mu d0signif";
	tempstring2 = "Mud0signif";
	PrintHisto(fMuHistos[1], tempstring1, tempstring2, outputdir, distVxmax);
	file << "* " << tempstring2 << endl;
	file << fAC->printAverage(tempstring2, fMuHistos[1]) << endl;
	file << fAC->printRatio(tempstring2, fMuHistos[1], distVxmax, 100., 0., 100.) << endl;

	tempstring1 = "Mu dzsignif";
	tempstring2 = "Mudzsignif";
	PrintHisto(fMuHistos[2], tempstring1, tempstring2, outputdir, distVxmax);
	file << "* " << tempstring2 << endl;
	file << fAC->printAverage(tempstring2, fMuHistos[2]) << endl;
	file << fAC->printRatio(tempstring2, fMuHistos[2], distVxmax, 100., 0., 100.) << endl;

	tempstring1 = "Mu NChi2";
	tempstring2 = "MuNChi2";
	PrintHisto(fMuHistos[3], tempstring1, tempstring2, outputdir, MuonChi2max);
	file << "* " << tempstring2 << endl;
	file << fAC->printAverage(tempstring2, fMuHistos[3]) << endl;

	tempstring1 = "Mu NTkHits";
	tempstring2 = "MuNTkHits";
	PrintHisto(fMuHistos[4], tempstring1, tempstring2, outputdir, MuonNHitsmin);
	file << "* " << tempstring2 << endl;
	file << fAC->printAverage(tempstring2, fMuHistos[4]) << endl;

	tempstring1 = "Mu DR SS";
	tempstring2 = "MuDRSS";
	PrintHisto(fMuHistos[5], tempstring1, tempstring2, outputdir, dRSSmuonmax);
	file << "* " << tempstring2 << endl;
	file << fAC->printAverage(tempstring2, fMuHistos[5]) << endl;

	tempstring1 = "Mu DR OS";
	tempstring2 = "MuDROS";
	PrintHisto(fMuHistos[6], tempstring1, tempstring2, outputdir, dRSSmuonmax);
	file << "* " << tempstring2 << endl;
	file << fAC->printAverage(tempstring2, fMuHistos[6]) << endl;

	tempstring1 = "Mu InvM SS";
	tempstring2 = "MuInvMSS";
	PrintHisto(fMuHistos[7], tempstring1, tempstring2, outputdir);
	file << "* " << tempstring2 << endl;
	file << fAC->printAverage(tempstring2, fMuHistos[7]) << endl;

	tempstring1 = "Mu InvM OS";
	tempstring2 = "MuInvMOS";
	PrintHisto(fMuHistos[8], tempstring1, tempstring2, outputdir);
	file << "* " << tempstring2 << endl;
	file << fAC->printAverage(tempstring2, fMuHistos[8]) << endl;

	// uncleaned electrons
	outputdir = fOutputDir + "Cleaning/El/";
	tempstring1 = "El H/E Bar";
	tempstring2 = "ElHcalOverEcalBar";
	PrintHisto(fElHistos[0], tempstring1, tempstring2, outputdir, ElecHoverEBarmax, 999., 1);
	file << "* " << tempstring2 << endl;
	file << fAC->printAverage(tempstring2, fElHistos[0]) << endl;
	file << fAC->printRatio(tempstring2, fElHistos[0], ElecHoverEBarmax, 100., 0., 100.) << endl;

	tempstring1 = "El SigmaIetaIeta Bar";
	tempstring2 = "ElSigmaIetaIetaBar";
	PrintHisto(fElHistos[1], tempstring1, tempstring2, outputdir, ElecSigmaEtaEtaBarmax);
	file << "* " << tempstring2 << endl;
	file << fAC->printAverage(tempstring2, fElHistos[1]) << endl;
	file << fAC->printRatio(tempstring2, fElHistos[1], ElecSigmaEtaEtaBarmax, 100., 0., 100.) << endl;

	tempstring1 = "El DeltaPhiSeedClusterAtCalo Bar";
	tempstring2 = "ElDeltaPhiSeedClusterAtCaloBar";
	PrintHisto(fElHistos[2], tempstring1, tempstring2, outputdir, ElecDeltaPhiOutBarmax, -ElecDeltaPhiOutBarmax);
	file << "* " << tempstring2 << endl;
	file << fAC->printAverage(tempstring2, fElHistos[2]) << endl;
	file << fAC->printRatio(tempstring2, fElHistos[2], -ElecDeltaPhiOutBarmax, ElecDeltaPhiOutBarmax, -100., 100.) << endl;

	tempstring1 = "El DeltaEtaSeedClusterAtCalo Bar";
	tempstring2 = "ElDeltaEtaSeedClusterAtCaloBar";
	PrintHisto(fElHistos[3], tempstring1, tempstring2, outputdir);
	file << "* " << tempstring2 << endl;
	file << fAC->printAverage(tempstring2, fElHistos[3]) << endl;

	tempstring1 = "El DeltaPhiSuperClusterAtVtx Bar";
	tempstring2 = "ElDeltaPhiSuperClusterAtVtxBar";
	PrintHisto(fElHistos[4], tempstring1, tempstring2, outputdir, ElecDeltaPhiInBarmax, -ElecDeltaPhiInBarmax);
	file << "* " << tempstring2 << endl;
	file << fAC->printAverage(tempstring2, fElHistos[4]) << endl;
	file << fAC->printRatio(tempstring2, fElHistos[4], -ElecDeltaPhiInBarmax, ElecDeltaPhiInBarmax, -100., 100.) << endl;

	tempstring1 = "El DeltaEtaSuperClusterAtVtx Bar";
	tempstring2 = "ElDeltaEtaSuperClusterAtVtxBar";
	PrintHisto(fElHistos[5], tempstring1, tempstring2, outputdir, ElecDeltaEtaInBarmax, -ElecDeltaEtaInBarmax);
	file << "* " << tempstring2 << endl;
	file << fAC->printAverage(tempstring2, fElHistos[5]) << endl;
	file << fAC->printRatio(tempstring2, fElHistos[5], -ElecDeltaEtaInBarmax, ElecDeltaEtaInBarmax, -100., 100.) << endl;

	tempstring1 = "El ESuperClusterOverP Bar";
	tempstring2 = "ElESuperClusterOverPBar";
	PrintHisto(fElHistos[6], tempstring1, tempstring2, outputdir, ElecEoverPInBarmin);
	file << "* " << tempstring2 << endl;
	file << fAC->printAverage(tempstring2, fElHistos[6]) << endl;
	file << fAC->printRatio(tempstring2, fElHistos[6], 0., ElecEoverPInBarmin, 0., 100.) << endl;

	tempstring1 = "El H/E End";
	tempstring2 = "ElHcalOverEcalEnd";
	PrintHisto(fElHistos[7], tempstring1, tempstring2, outputdir, ElecHoverEEndmax, 999., 1);
	file << "* " << tempstring2 << endl;
	file << fAC->printAverage(tempstring2, fElHistos[7]) << endl;
	file << fAC->printRatio(tempstring2, fElHistos[7], ElecHoverEEndmax, 100., 0., 100.) << endl;

	tempstring1 = "El SigmaIetaIeta End";
	tempstring2 = "ElSigmaIetaIetaEnd";
	PrintHisto(fElHistos[8], tempstring1, tempstring2, outputdir, ElecSigmaEtaEtaEndmax);
	file << "* " << tempstring2 << endl;
	file << fAC->printAverage(tempstring2, fElHistos[8]) << endl;
	file << fAC->printRatio(tempstring2, fElHistos[8], ElecSigmaEtaEtaEndmax, 100., 0., 100.) << endl;

	tempstring1 = "El DeltaPhiSeedClusterAtCalo End";
	tempstring2 = "ElDeltaPhiSeedClusterAtCaloEnd";
	PrintHisto(fElHistos[9], tempstring1, tempstring2, outputdir, ElecDeltaPhiOutEndmax, -ElecDeltaPhiOutEndmax);
	file << "* " << tempstring2 << endl;
	file << fAC->printAverage(tempstring2, fElHistos[9]) << endl;
	file << fAC->printRatio(tempstring2, fElHistos[9], -ElecDeltaPhiOutEndmax, ElecDeltaPhiOutEndmax, -100., 100.) << endl;

	tempstring1 = "El DeltaEtaSeedClusterAtCalo End";
	tempstring2 = "ElDeltaEtaSeedClusterAtCaloEnd";
	PrintHisto(fElHistos[10], tempstring1, tempstring2, outputdir);
	file << "* " << tempstring2 << endl;
	file << fAC->printAverage(tempstring2, fElHistos[10]) << endl;

	tempstring1 = "El DeltaPhiSuperClusterAtVtx End";
	tempstring2 = "ElDeltaPhiSuperClusterAtVtxEnd";
	PrintHisto(fElHistos[11], tempstring1, tempstring2, outputdir, ElecDeltaPhiInEndmax, -ElecDeltaPhiInEndmax);
	file << "* " << tempstring2 << endl;
	file << fAC->printAverage(tempstring2, fElHistos[11]) << endl;
	file << fAC->printRatio(tempstring2, fElHistos[11], -ElecDeltaPhiInEndmax, ElecDeltaPhiInEndmax, -100., 100.) << endl;

	tempstring1 = "El DeltaEtaSuperClusterAtVtx End";
	tempstring2 = "ElDeltaEtaSuperClusterAtVtxEnd";
	PrintHisto(fElHistos[12], tempstring1, tempstring2, outputdir, ElecDeltaEtaInEndmax, -ElecDeltaEtaInEndmax);
	file << "* " << tempstring2 << endl;
	file << fAC->printAverage(tempstring2, fElHistos[12]) << endl;
	file << fAC->printRatio(tempstring2, fElHistos[12], -ElecDeltaEtaInEndmax, ElecDeltaEtaInEndmax, -100., 100.) << endl;

	tempstring1 = "El ESuperClusterOverP End";
	tempstring2 = "ElESuperClusterOverPEnd";
	PrintHisto(fElHistos[13], tempstring1, tempstring2, outputdir, ElecEoverPInEndmin);
	file << "* " << tempstring2 << endl;
	file << fAC->printAverage(tempstring2, fElHistos[13]) << endl;
	file << fAC->printRatio(tempstring2, fElHistos[13], 0., ElecEoverPInEndmin, 0., 100.) << endl;

	tempstring1 = "El d0signif";
	tempstring2 = "Eld0signif";
	PrintHisto(fElHistos[14], tempstring1, tempstring2, outputdir, distVxmax);
	file << "* " << tempstring2 << endl;
	file << fAC->printAverage(tempstring2, fElHistos[14]) << endl;
	file << fAC->printRatio(tempstring2, fElHistos[14], distVxmax, 100., 0., 100.) << endl;

	tempstring1 = "El dzsignif";
	tempstring2 = "Eldzsignif";
	PrintHisto(fElHistos[15], tempstring1, tempstring2, outputdir, distVxmax);
	file << "* " << tempstring2 << endl;
	file << fAC->printAverage(tempstring2, fElHistos[15]) << endl;
	file << fAC->printRatio(tempstring2, fElHistos[15], distVxmax, 100., 0., 100.) << endl;

	tempstring1 = "El DR SS";
	tempstring2 = "ElDRSS";
	PrintHisto(fElHistos[16], tempstring1, tempstring2, outputdir);
	file << "* " << tempstring2 << endl;
	file << fAC->printAverage(tempstring2, fElHistos[16]) << endl;

	tempstring1 = "El DR OS";
	tempstring2 = "ElDROS";
	PrintHisto(fElHistos[17], tempstring1, tempstring2, outputdir);
	file << "* " << tempstring2 << endl;
	file << fAC->printAverage(tempstring2, fElHistos[17]) << endl;

	tempstring1 = "El InvM Same Sign";
	tempstring2 = "ElInvMSS";
	PrintHisto(fElHistos[18], tempstring1, tempstring2, outputdir);
	file << "* " << tempstring2 << endl;
	file << fAC->printAverage(tempstring2, fElHistos[18]) << endl;

	tempstring1 = "El InvM Opp Sign";
	tempstring2 = "ElInvMOS";
	PrintHisto(fElHistos[19], tempstring1, tempstring2, outputdir);
	file << "* " << tempstring2 << endl;
	file << fAC->printAverage(tempstring2, fElHistos[19]) << endl;

	// uncleaned photons
	outputdir = fOutputDir + "Cleaning/Ph/";
	tempstring1 = "Phot H/E Bar";
	tempstring2 = "PhoHoverEBar";
	PrintHisto(fPhHistos[0], tempstring1, tempstring2, outputdir, PhotHoverEBarmax, 999., 1);
	file << "* " << tempstring2 << endl;
	file << fAC->printAverage(tempstring2, fPhHistos[0]) << endl;
	file << fAC->printRatio(tempstring2, fPhHistos[0], PhotHoverEBarmax, 100., 0., 100.) << endl;

	tempstring1 = "Phot sigmaIetaIeta Bar";
	tempstring2 = "PhoSigmaIetaIetaBar";
	PrintHisto(fPhHistos[1], tempstring1, tempstring2, outputdir, PhotSigmaEtaEtaBarmax);
	file << "* " << tempstring2 << endl;
	file << fAC->printAverage(tempstring2, fPhHistos[1]) << endl;
	file << fAC->printRatio(tempstring2, fPhHistos[1], PhotSigmaEtaEtaBarmax, 100., 0., 100.) << endl;

	tempstring1 = "Phot H/E End";
	tempstring2 = "PhoHoverEEnd";
	PrintHisto(fPhHistos[2], tempstring1, tempstring2, outputdir, PhotHoverEEndmax, 999., 1);
	file << "* " << tempstring2 << endl;
	file << fAC->printAverage(tempstring2, fPhHistos[2]) << endl;
	file << fAC->printRatio(tempstring2, fPhHistos[2], PhotHoverEEndmax, 100., 0., 100.) << endl;

	tempstring1 = "Phot sigmaIetaIeta End";
	tempstring2 = "PhoSigmaIetaIetaEnd";
	PrintHisto(fPhHistos[3], tempstring1, tempstring2, outputdir, PhotSigmaEtaEtaEndmax);
	file << "* " << tempstring2 << endl;
	file << fAC->printAverage(tempstring2, fPhHistos[3]) << endl;
	file << fAC->printRatio(tempstring2, fPhHistos[3], PhotSigmaEtaEtaEndmax, 100., 0., 100.) << endl;

	tempstring1 = "Phot DR";
	tempstring2 = "PhDR";
	PrintHisto(fPhHistos[4], tempstring1, tempstring2, outputdir);
	file << "* " << tempstring2 << endl;
	file << fAC->printAverage(tempstring2, fPhHistos[4]) << endl;

	tempstring1 = "Phot InvM";
	tempstring2 = "PhInvM";
	PrintHisto(fPhHistos[5], tempstring1, tempstring2, outputdir);
	file << "* " << tempstring2 << endl;
	file << fAC->printAverage(tempstring2, fPhHistos[5]) << endl;

	tempstring1 = "Phot InvM Low";
	tempstring2 = "PhInvMLow";
	PrintHisto(fPhHistos[6], tempstring1, tempstring2, outputdir);
	file << "* " << tempstring2 << endl;
	file << fAC->printAverage(tempstring2, fPhHistos[6]) << endl;

	// uncleaned jets
	outputdir = fOutputDir + "Cleaning/Jet/";
	tempstring1 = "Jet d0 PV";
	tempstring2 = "Jetd0PV";
	PrintHisto(fJHistos[0], tempstring1, tempstring2, outputdir);
	file << "* " << tempstring2 << endl;
	file << fAC->printAverage(tempstring2, fJHistos[0]) << endl;
	file << fAC->printRatio(tempstring2, fJHistos[0], distVxmax, 100., 0., 100.) << endl;

	tempstring1 = "Jet dz PV";
	tempstring2 = "JetdzPV";
	PrintHisto(fJHistos[1], tempstring1, tempstring2, outputdir);
	file << "* " << tempstring2 << endl;
	file << fAC->printAverage(tempstring2, fJHistos[1]) << endl;
	file << fAC->printRatio(tempstring2, fJHistos[1], distVxmax, 100., 0., 100.) << endl;

	tempstring1 = "Jet d0signif";
	tempstring2 = "Jetd0signif";
	PrintHisto(fJHistos[2], tempstring1, tempstring2, outputdir, distVxmax);
	file << "* " << tempstring2 << endl;
	file << fAC->printAverage(tempstring2, fJHistos[2]) << endl;

	tempstring1 = "Jet dzsignif";
	tempstring2 = "Jetdzsignif";
	PrintHisto(fJHistos[3], tempstring1, tempstring2, outputdir, distVxmax);
	file << "* " << tempstring2 << endl;
	file << fAC->printAverage(tempstring2, fJHistos[3]) << endl;

	tempstring1 = "Jet EMfrac";
	tempstring2 = "JetEMfrac";
	PrintHisto(fJHistos[4], tempstring1, tempstring2, outputdir, FracEmminJet, FracEmmaxJet);
	file << "* " << tempstring2 << endl;
	file << fAC->printAverage(tempstring2, fJHistos[4]) << endl;

	tempstring1 = "Jet Chfrac";
	tempstring2 = "JetChfrac";
	PrintHisto(fJHistos[5], tempstring1, tempstring2, outputdir, FracChminJet);
	file << "* " << tempstring2 << endl;
	file << fAC->printAverage(tempstring2, fJHistos[5]) << endl;

	tempstring1 = "Jet ID n90Hits";
	tempstring2 = "JID_n90Hits";
	PrintHisto(fJHistos[6], tempstring1, tempstring2, outputdir, JID_n90Hitsmin);
	file << "* " << tempstring2 << endl;
	file << fAC->printAverage(tempstring2, fJHistos[6]) << endl;

	tempstring1 = "Jet ID HPD";
	tempstring2 = "JID_HPD";
	PrintHisto(fJHistos[7], tempstring1, tempstring2, outputdir, JID_HPDmax, 999., 1);
	file << "* " << tempstring2 << endl;
	file << fAC->printAverage(tempstring2, fJHistos[7]) << endl;

	tempstring1 = "Jet Ch vs EMfrac";
	tempstring2 = "JChEMfrac";
	canv = new TCanvas(tempstring2, tempstring1 , 0, 0, 900, 700);
	if(fJChEMfrac->GetEntries() < 100000) fJChEMfrac->SetMarkerStyle(6);
	fJChEMfrac->DrawCopy();
	Util::Print(canv, tempstring2, outputdir, fHistFile);

	tempstring1 = "Jet EMfrac vs Eta";
	tempstring2 = "JEMfracEta";
	canv = new TCanvas(tempstring2, tempstring1 , 0, 0, 900, 700);
	if(fJEMfracEta->GetEntries() < 100000) fJEMfracEta->SetMarkerStyle(6);
	fJEMfracEta->DrawCopy();
	Util::Print(canv, tempstring2, outputdir, fHistFile);

	tempstring1 = "Jet Chfrac vs Eta";
	tempstring2 = "JChfracEta";
	canv = new TCanvas(tempstring2, tempstring1 , 0, 0, 900, 700);
	if(fJChfracEta->GetEntries() < 100000) fJChfracEta->SetMarkerStyle(6);
	fJChfracEta->DrawCopy();
	Util::Print(canv, tempstring2, outputdir, fHistFile);

	// uncleaned Event and MET
	outputdir = fOutputDir + "Cleaning/MET/";
	tempstring1 = "Dphi(MET,j1)";
	tempstring2 = "DphiMETJet1";
	PrintHisto(fMETHistos[0], tempstring1, tempstring2, outputdir, dPhiJetMETmin);
	file << "* " << tempstring2 << endl;
	file << fAC->printAverage(tempstring2, fMETHistos[0]) << endl;

	tempstring1 = "Dphi(MET,j2)";
	tempstring2 = "DphiMETJet2";
	PrintHisto(fMETHistos[1], tempstring1, tempstring2, outputdir, dPhiJetMETmin);
	file << "* " << tempstring2 << endl;
	file << fAC->printAverage(tempstring2, fMETHistos[1]) << endl;

	tempstring1 = "MET R12";
	tempstring2 = "METR12";
	PrintHisto(fMETHistos[2], tempstring1, tempstring2, outputdir, dR12min);
	file << "* " << tempstring2 << endl;
	file << fAC->printAverage(tempstring2, fMETHistos[2]) << endl;

	tempstring1 = "MET R21";
	tempstring2 = "METR21";
	PrintHisto(fMETHistos[3], tempstring1, tempstring2, outputdir, dR21min);
	file << "* " << tempstring2 << endl;
	file << fAC->printAverage(tempstring2, fMETHistos[3]) << endl;

	tempstring1 = "MET Rsum";
	tempstring2 = "METRsum";
	PrintHisto(fMETHistos[4], tempstring1, tempstring2, outputdir);
	file << "* " << tempstring2 << endl;
	file << fAC->printAverage(tempstring2, fMETHistos[4]) << endl;

	tempstring1 = "MET Dphi 12";
	tempstring2 = "METDphi12";
	canv = new TCanvas(tempstring2, tempstring1 , 0, 0, 900, 700);
	if(fMETDphi12->GetEntries() < 100000) fMETDphi12->SetMarkerStyle(6);
	fMETDphi12->DrawCopy();
	Util::Print(canv, tempstring2, outputdir, fHistFile);

	tempstring1 = "MET R12 R21";
	tempstring2 = "METR12R21";
	canv = new TCanvas(tempstring2, tempstring1 , 0, 0, 900, 700);
	if(fMETR12R21->GetEntries() < 100000) fMETR12R21->SetMarkerStyle(6);
	fMETR12R21->DrawCopy();
	Util::Print(canv, tempstring2, outputdir, fHistFile);

	tempstring1 = "MET R12 Dphij12";
	tempstring2 = "METR12Dphij12";
	canv = new TCanvas(tempstring2, tempstring1 , 0, 0, 900, 700);
	if(fMETR12Dphij12->GetEntries() < 100000) fMETR12Dphij12->SetMarkerStyle(6);
	fMETR12Dphij12->DrawCopy();
	Util::Print(canv, tempstring2, outputdir, fHistFile);

	tempstring1 = "MET R12 Etaj1,2";
	tempstring2 = "METR12Etaj12";
	canv = new TCanvas(tempstring2, tempstring1 , 0, 0, 900, 700);
	if(fMETR12Etaj12->GetEntries() < 100000) fMETR12Etaj12->SetMarkerStyle(6);
	fMETR12Etaj12->DrawCopy();
	Util::Print(canv, tempstring2, outputdir, fHistFile);

	tempstring1 = "MET R12 Ptj1/Ptj2";
	tempstring2 = "METR12j12PtRat";
	canv = new TCanvas(tempstring2, tempstring1 , 0, 0, 900, 700);
	if(fMETR12j12PtRat->GetEntries() < 100000) fMETR12j12PtRat->SetMarkerStyle(6);
	fMETR12j12PtRat->DrawCopy();
	Util::Print(canv, tempstring2, outputdir, fHistFile);

	tempstring1 = "MET R12 dR(j1,j2)";
	tempstring2 = "METR12dRj12";
	canv = new TCanvas(tempstring2, tempstring1 , 0, 0, 900, 700);
	if(fMETR12dRj12->GetEntries() < 100000) fMETR12dRj12->SetMarkerStyle(6);
	fMETR12dRj12->DrawCopy();
	Util::Print(canv, tempstring2, outputdir, fHistFile);

	tempstring1 = "Evt Em Frac";
	tempstring2 = "EvtEmFrac";
	PrintHisto(fEvtHistos[0], tempstring1, tempstring2, outputdir, FracEmmin, 999., 1);
	file << "* " << tempstring2 << endl;
	file << fAC->printAverage(tempstring2, fEvtHistos[0]) << endl;

	tempstring1 = "Evt Ch Frac";
	tempstring2 = "EvtChFrac";
	PrintHisto(fEvtHistos[1], tempstring1, tempstring2, outputdir, FracChmin, 999., 1);
	file << "* " << tempstring2 << endl;
	file << fAC->printAverage(tempstring2, fEvtHistos[1]) << endl;

	// Clean+iso objects distributions
	outputdir = fOutputDir + "Cleaning/Obj/";
	tempstring1 = "Number of Mus";
	tempstring2 = "NMus";
	PrintHisto(fMuCIHistos[0], tempstring1, tempstring2, outputdir, 999., 999., 1);
	file << "* " << tempstring2 << endl;
	file << fAC->printAverage(tempstring2, fMuCIHistos[0]) << endl;

	tempstring1 = "Muon Pt";
	tempstring2 = "MuPt";
	PrintHisto(fMuCIHistos[1], tempstring1, tempstring2, outputdir, MuPtmin, 999., 1, 1);
	file << "* " << tempstring2 << endl;
	file << fAC->printAverage(tempstring2, fMuCIHistos[1]) << endl;

	tempstring1 = "Muon Eta";
	tempstring2 = "MuEta";
	PrintHisto(fMuCIHistos[2], tempstring1, tempstring2, outputdir, -MuEtamax, MuEtamax);
	file << "* " << tempstring2 << endl;
	file << fAC->printAverage(tempstring2, fMuCIHistos[2]) << endl;

	tempstring1 = "Muon Phi";
	tempstring2 = "MuPhi";
	PrintHisto(fMuCIHistos[3], tempstring1, tempstring2, outputdir);
	file << "* " << tempstring2 << endl;
	file << fAC->printAverage(tempstring2, fMuCIHistos[3]) << endl;

	tempstring1 = "Number of Eles";
	tempstring2 = "NEles";
	PrintHisto(fElCIHistos[0], tempstring1, tempstring2, outputdir, 999., 999., 1);
	file << "* " << tempstring2 << endl;
	file << fAC->printAverage(tempstring2, fElCIHistos[0]) << endl;

	tempstring1 = "Elec Pt";
	tempstring2 = "ElPt";
	PrintHisto(fElCIHistos[1], tempstring1, tempstring2, outputdir, ElPtmin, 999., 1, 1);
	file << "* " << tempstring2 << endl;
	file << fAC->printAverage(tempstring2, fElCIHistos[1]) << endl;

	tempstring1 = "Elec Eta";
	tempstring2 = "ElEta";
	PrintHisto(fElCIHistos[2], tempstring1, tempstring2, outputdir, -ElEtamax, ElEtamax);
	file << "* " << tempstring2 << endl;
	file << fAC->printAverage(tempstring2, fElCIHistos[2]) << endl;

	tempstring1 = "Elec Phi";
	tempstring2 = "ElPhi";
	PrintHisto(fElCIHistos[3], tempstring1, tempstring2, outputdir);
	file << "* " << tempstring2 << endl;
	file << fAC->printAverage(tempstring2, fElCIHistos[3]) << endl;

	tempstring1 = "Elec SigmaIetaIeta";
	tempstring2 = "ElSigmaIetaIetaClean";
	PrintHisto(fElCIHistos[4], tempstring1, tempstring2, outputdir);
	file << "* " << tempstring2 << endl;
	file << fAC->printAverage(tempstring2, fElCIHistos[4]) << endl;

	tempstring1 = "Number of Photons";
	tempstring2 = "NPhotons";
	PrintHisto(fPhCIHistos[0], tempstring1, tempstring2, outputdir, 999., 999., 1);
	file << "* " << tempstring2 << endl;
	file << fAC->printAverage(tempstring2, fPhCIHistos[0]) << endl;

	tempstring1 = "Pho Pt";
	tempstring2 = "PhoPt";
	PrintHisto(fPhCIHistos[1], tempstring1, tempstring2, outputdir, PhPtmin, 999., 1, 1);
	file << "* " << tempstring2 << endl;
	file << fAC->printAverage(tempstring2, fPhCIHistos[1]) << endl;

	tempstring1 = "Pho Eta";
	tempstring2 = "PhoEta";
	PrintHisto(fPhCIHistos[2], tempstring1, tempstring2, outputdir, -PhEtamax, PhEtamax);
	file << "* " << tempstring2 << endl;
	file << fAC->printAverage(tempstring2, fPhCIHistos[2]) << endl;

	tempstring1 = "Pho Phi";
	tempstring2 = "PhoPhi";
	PrintHisto(fPhCIHistos[3], tempstring1, tempstring2, outputdir);
	file << "* " << tempstring2 << endl;
	file << fAC->printAverage(tempstring2, fPhCIHistos[3]) << endl;

	tempstring1 = "Number of Jets";
	tempstring2 = "NJets";
	PrintHisto(fJCIHistos[0], tempstring1, tempstring2, outputdir, 999., 999., 1);
	file << "* " << tempstring2 << endl;
	file << fAC->printAverage(tempstring2, fJCIHistos[0]) << endl;

	tempstring1 = "Jet Pt";
	tempstring2 = "JPt";
	PrintHisto(fJCIHistos[1], tempstring1, tempstring2, outputdir, JPtmin, 999., 1, 1);
	file << "* " << tempstring2 << endl;
	file << fAC->printAverage(tempstring2, fJCIHistos[1]) << endl;

	tempstring1 = "Jet Eta";
	tempstring2 = "JEta";
	PrintHisto(fJCIHistos[2], tempstring1, tempstring2, outputdir, -JEtamax, JEtamax);
	file << "* " << tempstring2 << endl;
	file << fAC->printAverage(tempstring2, fJCIHistos[2]) << endl;

	tempstring1 = "Jet Phi";
	tempstring2 = "JPhi";
	PrintHisto(fJCIHistos[3], tempstring1, tempstring2, outputdir);
	file << "* " << tempstring2 << endl;
	file << fAC->printAverage(tempstring2, fJCIHistos[3]) << endl;

	tempstring1 = "Jet EMfrac";
	tempstring2 = "JetEMfracClean";
	PrintHisto(fJCIHistos[4], tempstring1, tempstring2, outputdir, FracEmminJet, FracEmmaxJet);
	file << "* " << tempstring2 << endl;
	file << fAC->printAverage(tempstring2, fJCIHistos[4]) << endl;

	tempstring1 = "Jet Chfrac";
	tempstring2 = "JetChfracClean";
	PrintHisto(fJCIHistos[5], tempstring1, tempstring2, outputdir, FracChminJet);
	file << "* " << tempstring2 << endl;
	file << fAC->printAverage(tempstring2, fJCIHistos[5]) << endl;

	tempstring1 = "Dphi(MET,j1)";
	tempstring2 = "DphiMETJet1Clean";
	PrintHisto(fJCIHistos[6], tempstring1, tempstring2, outputdir, dPhiJetMETmin);
	file << "* " << tempstring2 << endl;
	file << fAC->printAverage(tempstring2, fJCIHistos[6]) << endl;

	tempstring1 = "Dphi(MET,j2)";
	tempstring2 = "DphiMETJet2Clean";
	PrintHisto(fJCIHistos[7], tempstring1, tempstring2, outputdir, dPhiJetMETmin);
	file << "* " << tempstring2 << endl;
	file << fAC->printAverage(tempstring2, fJCIHistos[7]) << endl;

	tempstring1 = "MET R12";
	tempstring2 = "METR12Clean";
	PrintHisto(fJCIHistos[8], tempstring1, tempstring2, outputdir, dR12min);
	file << "* " << tempstring2 << endl;
	file << fAC->printAverage(tempstring2, fJCIHistos[8]) << endl;

	tempstring1 = "MET R21";
	tempstring2 = "METR21Clean";
	PrintHisto(fJCIHistos[9], tempstring1, tempstring2, outputdir, dR12min);
	file << "* " << tempstring2 << endl;
	file << fAC->printAverage(tempstring2, fJCIHistos[9]) << endl;

	tempstring1 = "mono-Jet Pt";
	tempstring2 = "JPt1j";
	PrintHisto(fJCIHistos[10], tempstring1, tempstring2, outputdir, JPtmin, 999., 1, 1);
	file << "* " << tempstring2 << endl;
	file << fAC->printAverage(tempstring2, fJCIHistos[10]) << endl;

	tempstring1 = "mono-Jet Eta";
	tempstring2 = "JEta1j";
	PrintHisto(fJCIHistos[11], tempstring1, tempstring2, outputdir);
	file << "* " << tempstring2 << endl;
	file << fAC->printAverage(tempstring2, fJCIHistos[11]) << endl;

	tempstring1 = "mono-Jet Phi";
	tempstring2 = "JPhi1j";
	PrintHisto(fJCIHistos[12], tempstring1, tempstring2, outputdir);
	file << "* " << tempstring2 << endl;
	file << fAC->printAverage(tempstring2, fJCIHistos[12]) << endl;

	tempstring1 = "Dphi(MET,mono-jet)";
	tempstring2 = "DphiMETmonoJetClean";
	PrintHisto(fJCIHistos[13], tempstring1, tempstring2, outputdir, dPhiJetMETmin);
	file << "* " << tempstring2 << endl;
	file << fAC->printAverage(tempstring2, fJCIHistos[13]) << endl;

	// Clean+iso objects, Inv Mass distributions
	tempstring1 = "Mu Inv Mass SS Clean";
	tempstring2 = "MuInvMSSClean";
	PrintHisto(fInvMHistos[0], tempstring1, tempstring2, outputdir);
	file << "* " << tempstring2 << endl;
	file << fAC->printAverage(tempstring2, fInvMHistos[0]) << endl;

	tempstring1 = "Mu Low Inv Mass SS Clean";
	tempstring2 = "MuInvMSSCleanLow";
	PrintHisto(fInvMHistos[1], tempstring1, tempstring2, outputdir);
	file << "* " << tempstring2 << endl;
	file << fAC->printAverage(tempstring2, fInvMHistos[1]) << endl;

	tempstring1 = "Mu Inv Mass OS Clean";
	tempstring2 = "MuInvMOSClean";
	PrintHisto(fInvMHistos[2], tempstring1, tempstring2, outputdir);
	file << "* " << tempstring2 << endl;
	file << fAC->printAverage(tempstring2, fInvMHistos[2]) << endl;

	tempstring1 = "Mu Low Inv Mass OS Clean";
	tempstring2 = "MuInvMOSCleanLow";
	PrintHisto(fInvMHistos[3], tempstring1, tempstring2, outputdir);
	file << "* " << tempstring2 << endl;
	file << fAC->printAverage(tempstring2, fInvMHistos[3]) << endl;

	tempstring1 = "El Inv Mass SS Clean";
	tempstring2 = "ElInvMSSClean";
	PrintHisto(fInvMHistos[4], tempstring1, tempstring2, outputdir);
	file << "* " << tempstring2 << endl;
	file << fAC->printAverage(tempstring2, fInvMHistos[4]) << endl;

	tempstring1 = "El Low Inv Mass SS Clean";
	tempstring2 = "ElInvMSSCleanLow";
	PrintHisto(fInvMHistos[5], tempstring1, tempstring2, outputdir);
	file << "* " << tempstring2 << endl;
	file << fAC->printAverage(tempstring2, fInvMHistos[5]) << endl;

	tempstring1 = "El Inv Mass OS Clean";
	tempstring2 = "ElInvMOSClean";
	PrintHisto(fInvMHistos[6], tempstring1, tempstring2, outputdir);
	file << "* " << tempstring2 << endl;
	file << fAC->printAverage(tempstring2, fInvMHistos[6]) << endl;

	tempstring1 = "El Low Inv Mass OS Clean";
	tempstring2 = "ElInvMOSCleanLow";
	PrintHisto(fInvMHistos[7], tempstring1, tempstring2, outputdir);
	file << "* " << tempstring2 << endl;
	file << fAC->printAverage(tempstring2, fInvMHistos[7]) << endl;

	tempstring1 = "Ph Inv Mass Clean";
	tempstring2 = "PhInvMClean";
	PrintHisto(fInvMHistos[8], tempstring1, tempstring2, outputdir);
	file << "* " << tempstring2 << endl;
	file << fAC->printAverage(tempstring2, fInvMHistos[8]) << endl;

	tempstring1 = "Ph Low Inv Mass Clean";
	tempstring2 = "PhInvMCleanLow";
	PrintHisto(fInvMHistos[9], tempstring1, tempstring2, outputdir);
	file << "* " << tempstring2 << endl;
	file << fAC->printAverage(tempstring2, fInvMHistos[9]) << endl;

	tempstring1 = "diJet Inv Mass Clean";
	tempstring2 = "JInvM2jClean";
	PrintHisto(fInvMHistos[10], tempstring1, tempstring2, outputdir);
	file << "* " << tempstring2 << endl;
	file << fAC->printAverage(tempstring2, fInvMHistos[10]) << endl;

	tempstring1 = "GT diJet Inv Mass Clean";
	tempstring2 = "JInvMGT2jClean";
	PrintHisto(fInvMHistos[11], tempstring1, tempstring2, outputdir);
	file << "* " << tempstring2 << endl;
	file << fAC->printAverage(tempstring2, fInvMHistos[11]) << endl;

	file.close();

	// write out the cleaning statistics
	TString cleanerstatfile = fOutputDir + "cleanerStats.txt";
	fTC->StatWrite(cleanerstatfile);

	// write the info file
	const int nentries = fTR->fChain->GetEntries();
	PrintInfoStart(nentries);

		// Close files
	fHistFile->Close();

}

void PhysQCAnalysis::BookHistos(){
	// Histos for uncleaned objects
	const unsigned int nvxbins = 100;
	fPvxHistos[0] = fAC->bookTH1D("PrimVtxNChi2", "Prim Vtx NChi2", nvxbins, 0., 10.);	
	fPvxHistos[1] = fAC->bookTH1D("PrimVtxd0BS", "Prim Vtx d0 to BS", nvxbins, 0., 0.3 );
	fPvxHistos[2] = fAC->bookTH1D("PrimVtxdzBS", "Prim Vtx dz to BS", nvxbins, 0., 20.);
	fPvxHistos[3] = fAC->bookTH1D("PrimVtxPtSum", "Prim Vtx PtSum", nvxbins, 0., 300.);

	const unsigned int nmubins = 100;
	fMuHistos[0] = fAC->bookTH1D("MuDeltaPOverP", "Mu DeltaP/P;Mu DeltaP/P",    nmubins, 0., 0.25 );
	fMuHistos[1] = fAC->bookTH1D("Mud0signif", "Mu d0 significance;Mu d0 significance",  nmubins, 0., 15. );
	fMuHistos[2] = fAC->bookTH1D("Mudzsignif", "Mu dz significance;Mu dz significance",  nmubins, 0., 15. );
	fMuHistos[3] = fAC->bookTH1D("MuNChi2", "Mu NChi2;Mu NChi2",  nmubins, 0., 20. );
	fMuHistos[4] = fAC->bookTH1D("MuNTkHits", "Mu NTkHits;Mu NTkHits",  30, 0., 30 );
	fMuHistos[5] = fAC->bookTH1D("MuDRSS", "Mu DR Same Sign;Mu DR Same Sign",  nmubins, 0., 3.2 );
	fMuHistos[6] = fAC->bookTH1D("MuDROS", "Mu DR Opp Sign;Mu DR Opp Sign",  nmubins, 0., 3.2 );
	fMuHistos[7] = fAC->bookTH1D("MuInvMSS", "Mu InvM Same Sign;Mu InvM Same Sign", nmubins, 0., 100. );
	fMuHistos[8] = fAC->bookTH1D("MuInvMOS", "Mu InvM Opp Sign;Mu InvM Opp Sign",  nmubins, 0., 100. );

	const unsigned int nelbins = 100;
	fElHistos[0]  = fAC->bookTH1D("ElHcalOverEcalBar", "El H/E Bar",  nelbins, 0., 0.2);
	fElHistos[1]  = fAC->bookTH1D("ElSigmaIetaIetaBar", "El sigmaIetaIeta Bar",  nelbins, 0., 0.05);
	fElHistos[2]  = fAC->bookTH1D("ElDeltaPhiSeedClusterAtCaloBar", "El DeltaPhiSeedClusterAtCalo Bar",  nelbins, -0.1, 0.1);
	fElHistos[3]  = fAC->bookTH1D("ElDeltaEtaSeedClusterAtCaloBar", "El DeltaEtaSeedClusterAtCalo Bar",  nelbins, -0.05, 0.05);
	fElHistos[4]  = fAC->bookTH1D("ElDeltaPhiSuperClusterAtVtxBar", "El DeltaPhiSuperClusterAtVtx Bar",  nelbins, -0.1, 0.1);
	fElHistos[5]  = fAC->bookTH1D("ElDeltaEtaSuperClusterAtVtxBar", "El DeltaEtaSuperClusterAtVtx Bar",  nelbins, -0.05, 0.05);
	fElHistos[6]  = fAC->bookTH1D("ElESuperClusterOverPBar", "El ESuperClusterOverP Bar",  nelbins, 0., 5.);
	fElHistos[7]  = fAC->bookTH1D("ElHcalOverEcalEnd", "El H/E End",  nelbins, 0., 0.2);
	fElHistos[8]  = fAC->bookTH1D("ElSigmaIetaIetaEnd", "El sigmaIetaIeta End",  nelbins, 0., 0.05);
	fElHistos[9]  = fAC->bookTH1D("ElDeltaPhiSeedClusterAtCaloEnd", "El DeltaPhiSeedClusterAtCalo End",  nelbins, -0.1, 0.1);
	fElHistos[10] = fAC->bookTH1D("ElDeltaEtaSeedClusterAtCaloEnd", "El DeltaEtaSeedClusterAtCalo End",  nelbins, -0.05, 0.05);
	fElHistos[11] = fAC->bookTH1D("ElDeltaPhiSuperClusterAtVtxEnd", "El DeltaPhiSuperClusterAtVtx End",  nelbins, -0.1, 0.1);
	fElHistos[12] = fAC->bookTH1D("ElDeltaEtaSuperClusterAtVtxEnd", "El DeltaEtaSuperClusterAtVtx End",  nelbins, -0.05, 0.05);
	fElHistos[13] = fAC->bookTH1D("ElESuperClusterOverPEnd", "El ESuperClusterOverP End",  nelbins, 0., 5.);
	fElHistos[14] = fAC->bookTH1D("Eld0signif", "El d0 significance;El d0 significance",  nelbins, 0., 15. );
	fElHistos[15] = fAC->bookTH1D("Eldzsignif", "El dz significance;El dz significance",  nelbins, 0., 15. );
	fElHistos[16] = fAC->bookTH1D("ElDRSS", "El DR Same Sign;El DR Same Sign",  nelbins, 0., 3.2 );
	fElHistos[17] = fAC->bookTH1D("ElDROS", "El DR Opp Sign;El DR Opp Sign",  nelbins, 0., 3.2 );
	fElHistos[18] = fAC->bookTH1D("ElInvMSS", "El InvM Same Sign;El InvM Same Sign", nelbins, 0., 100. );
	fElHistos[19] = fAC->bookTH1D("ElInvMOS", "El InvM Opp Sign;El InvM Opp Sign",  nelbins, 0., 100. );

	const unsigned int nphbins = 100;
	fPhHistos[0] = fAC->bookTH1D("PhoHoverEBar", "Phot H/E Bar",  nphbins, 0., 0.2);
	fPhHistos[1] = fAC->bookTH1D("PhoSigmaIetaIetaBar", "Phot sigmaIetaIeta Bar",  nphbins, 0., 0.05);
	fPhHistos[2] = fAC->bookTH1D("PhoHoverEEnd", "Phot H/E End",  nphbins, 0., 0.05);
	fPhHistos[3] = fAC->bookTH1D("PhoSigmaIetaIetaEnd", "Phot sigmaIetaIeta End",  nphbins, 0., 0.05);
	fPhHistos[4] = fAC->bookTH1D("PhDR", "Phot DR;Phot DR",  nphbins, 0., 3.2 );
	fPhHistos[5] = fAC->bookTH1D("PhInvM", "Phot InvM;Phot InvM", nphbins, 0., 50. );
	fPhHistos[6] = fAC->bookTH1D("PhInvMLow", "Phot InvM;Phot InvM", nphbins, 0., 5. );

	const unsigned int njetbins = 100;
	fJHistos[0] = fAC->bookTH1D("Jetd0PV", "Jet d0 PV;Jet d0 to PV",  njetbins, 0., 0.5 );
	fJHistos[1] = fAC->bookTH1D("JetdzPV", "Jet dz PV;Jet dz to PV",  njetbins, -2., 2. );
	fJHistos[2] = fAC->bookTH1D("Jetd0signif", "Jet d0 significance;Jet d0 significance",  njetbins, 0., 15. );
	fJHistos[3] = fAC->bookTH1D("Jetdzsignif", "Jet dz significance;Jet dz significance",  njetbins, 0., 15. );
	fJHistos[4] = fAC->bookTH1D("JEMfrac", "Jet EM fraction;Jet EM fraction",  njetbins, -0.01, 1.01 );
	fJHistos[5] = fAC->bookTH1D("JChfrac", "Jet Charged fraction;Jet Charged fraction",  njetbins, -1.01, 5. );
	fJHistos[6] = fAC->bookTH1D("JID_n90Hits", "Jet ID n90Hits;Jet ID n90Hits",  30, 0, 30 );
	fJHistos[7] = fAC->bookTH1D("JID_HPD", "Jet ID HPD;Jet ID HPD",  njetbins, 0, 1. );
	fJChEMfrac  = new TH2D("JChEMfrac", "Jet Ch vs EM frac;Jet EM frac;Jet Ch frac", njetbins, -0.01, 1.01, njetbins, -1.01, 2.);
	fJChEMfrac->SetXTitle(fAC->convertVarName("JEMfrac"));
	fJChEMfrac->SetYTitle(fAC->convertVarName("JChfrac"));
	fJEMfracEta = new TH2D("JEMfracEta", "Jet EM frac vs Eta;Jet Eta;Jet EM frac", njetbins, -5., 5., njetbins, -1.01, 1.01);
	fJEMfracEta->SetXTitle(fAC->convertVarName("JEta"));
	fJEMfracEta->SetYTitle(fAC->convertVarName("JEMfrac"));
	fJChfracEta = new TH2D("JChfracEta", "Jet Ch frac vs Eta;Jet Eta;Jet Ch frac", njetbins, -5., 5., njetbins, -1.01, 2.01);
	fJChfracEta->SetXTitle(fAC->convertVarName("JEta"));
	fJChfracEta->SetYTitle(fAC->convertVarName("JChfrac"));

	const unsigned int nevbins = 100;
	fMETHistos[0]   = fAC->bookTH1D("DphiMETJet1", "Dphi(MET,j1)", nevbins, 0, 3.1416);
	fMETHistos[1]   = fAC->bookTH1D("DphiMETJet2", "Dphi(MET,j2)", nevbins, 0, 3.1416);
	fMETHistos[2]   = fAC->bookTH1D("METR12", "MET R12;MET R12", nevbins, 0, 6.283);
	fMETHistos[3]   = fAC->bookTH1D("METR21", "MET R21;MET R21", nevbins, 0, 6.283);
	fMETHistos[4]   = fAC->bookTH1D("METRsum", "MET Rsum;MET Rsum", nevbins, 3., 6.283);

	fMETDphi12      = new TH2D("METDphi12", "MET Dphi 12;Dphi(MET,j1);Dphi(MET,j2)", nevbins, 0, 3.142, nevbins, 0, 3.142);
	fMETDphi12->SetXTitle(fAC->convertVarName("DphiMETJet1"));
	fMETDphi12->SetYTitle(fAC->convertVarName("DphiMETJet2"));
	fMETR12R21      = new TH2D("METR12R21", "MET R21 vs R12;R12;R21", nevbins, 0, 6.283, nevbins, 0, 6.283);
	fMETR12R21->SetXTitle(fAC->convertVarName("METR12"));
	fMETR12R21->SetYTitle(fAC->convertVarName("METR21"));
	fMETR12Dphij12  = new TH2D("METR12Dphij12", "MET R21 vs Dphijet12;R12;Dphij12", nevbins, 0, 6.283, nevbins, 0, 3.142);
	fMETR12Dphij12->SetXTitle(fAC->convertVarName("METR12"));
	fMETR12Dphij12->SetYTitle(fAC->convertVarName("DphiJ1J21"));
	fMETR12Etaj12   = new TH2D("METR12Etaj12", "MET R21 vs Eta jet12;R12;Eta jet 1,2", nevbins, 0, 6.283, nevbins, -5., 5.);
	fMETR12Etaj12->SetXTitle(fAC->convertVarName("METR12"));
	fMETR12Etaj12->SetYTitle(fAC->convertVarName("DetaJ1J21"));
	fMETR12j12PtRat = new TH2D("METR12j12PtRat", "MET R21 vs Eta jet12;R12;Pt2/Pt1", nevbins, 0, 6.283, nevbins, 0., 1.);
	fMETR12j12PtRat->SetXTitle(fAC->convertVarName("METR12"));
	fMETR12j12PtRat->SetYTitle(fAC->convertVarName("J1J2PtRatio"));
	fMETR12dRj12    = new TH2D("METR12dRj12", "MET R21 vs dR(j1,j2);R12;dR(j1,j2)", nevbins, 0, 6.283, nevbins, 0., 6.283);
	fMETR12dRj12->SetXTitle(fAC->convertVarName("METR12"));
	fMETR12dRj12->SetYTitle(fAC->convertVarName("DRJ1J2"));

	fEvtHistos[0]   = fAC->bookTH1D("EvtEmFrac", "Evt Em Frac;Evt Em Frac", nevbins, 0, 1.05);
	fEvtHistos[1]   = fAC->bookTH1D("EvtChFrac", "Evt Ch Frac;Evt Ch Frac", nevbins, 0, 1.05);

	// Histos for cleaned objects, CI = Clean+Isolated
	const unsigned int nCIbins = 100;
	fMuCIHistos[0] = fAC->bookTH1D("NMus", "Nber of Mus clean+iso", 10, 0, 10);
	fMuCIHistos[1] = fAC->bookTH1D("MuPt", "Pt of Mus clean+iso", nCIbins, 0, 100.);
	fMuCIHistos[2] = fAC->bookTH1D("MuEta", "Eta of Mus clean+iso", nCIbins, -2.5, 2.5);
	fMuCIHistos[3] = fAC->bookTH1D("MuPhi", "Phi of Mus clean+iso", nCIbins, -3.2, 3.2);
	fElCIHistos[0] = fAC->bookTH1D("NEles", "Nber of Eles clean+iso", 10, 0, 10);
	fElCIHistos[1] = fAC->bookTH1D("ElPt", "Pt of Eles clean+iso", nCIbins, 0, 100.);
	fElCIHistos[2] = fAC->bookTH1D("ElEta", "Eta of Eles clean+iso", nCIbins, -2.5, 2.5);
	fElCIHistos[3] = fAC->bookTH1D("ElPhi", "Phi of Eles clean+iso", nCIbins, -3.2, 3.2);
	fElCIHistos[4] = fAC->bookTH1D("ElSigmaIetaIetaClean", "ElSigmaIetaIeta of Eles clean+iso", nCIbins, 0., 0.015);
	fPhCIHistos[0] = fAC->bookTH1D("NPhotons", "Nber of Photons clean+iso", 20, 0, 20);
	fPhCIHistos[1] = fAC->bookTH1D("PhoPt", "Pt of Photons clean+iso", nCIbins, 0, 100.);
	fPhCIHistos[2] = fAC->bookTH1D("PhoEta", "Eta of Photons clean+iso", nCIbins, -2.5, 2.5);
	fPhCIHistos[3] = fAC->bookTH1D("PhoPhi", "Phi of Photons clean+iso", nCIbins, -3.2, 3.2);
	fJCIHistos[0]  = fAC->bookTH1D("NJets", "Nber of Jets clean+iso", 20, 0, 20);
	fJCIHistos[1]  = fAC->bookTH1D("JPt", "Pt of Jets clean+iso", nCIbins, 0, 250.);
	fJCIHistos[2]  = fAC->bookTH1D("JEta", "Eta of Jets clean+iso", nCIbins, -5., 5.);
	fJCIHistos[3]  = fAC->bookTH1D("JPhi", "Phi of Jets clean+iso", nCIbins, -3.2, 3.2);
	fJCIHistos[4]  = fAC->bookTH1D("JEMfracClean", "Jet EM fraction;Jet EM fraction",  njetbins, -0.01, 1.01 );
	fJCIHistos[5]  = fAC->bookTH1D("JChfracClean", "Jet Charged fraction;Jet Charged fraction",  njetbins, -1.01, 5. );
	fJCIHistos[6]  = fAC->bookTH1D("DphiMETJet1Clean", "Dphi(MET,j1)", nevbins, 0, 3.1416);
	fJCIHistos[7]  = fAC->bookTH1D("DphiMETJet2Clean", "Dphi(MET,j2)", nevbins, 0, 3.1416);
	fJCIHistos[8]  = fAC->bookTH1D("METR12Clean", "MET R12;MET R12", nevbins, 0, 6.283);
	fJCIHistos[9]  = fAC->bookTH1D("METR21Clean", "MET R21;MET R21", nevbins, 0, 6.283);
	fJCIHistos[10] = fAC->bookTH1D("JPt1j", "Pt of mono-Jets clean+iso", nCIbins, 0, 250.);
	fJCIHistos[11] = fAC->bookTH1D("JEta1j", "Eta of mono-Jets clean+iso", nCIbins, -5., 5.);
	fJCIHistos[12] = fAC->bookTH1D("JPhi1j", "Phi of mono-Jets clean+iso", nCIbins, -3.2, 3.2);
	fJCIHistos[13] = fAC->bookTH1D("DphiMETmonoJetClean", "Dphi(MET,monoJet)", nevbins, 0, 3.1416);

	// Inv mass Histos for cleaned objects, CI = Clean+Isolated
	const unsigned int ninvmbins = 100;
	fInvMHistos[0]  = fAC->bookTH1D("MuInvMSSClean", "Mu Clean InvM Same Sign;Mu InvM Same Sign", ninvmbins, 0., 100. );
	fInvMHistos[1]  = fAC->bookTH1D("MuInvMSSCleanLow", "Mu Clean InvM Low Same Sign;Mu InvM Same Sign", ninvmbins, 0., 5. );
	fInvMHistos[2]  = fAC->bookTH1D("MuInvMOSClean", "Mu Clean InvM Opp Sign;Mu InvM Opp Sign",  ninvmbins, 0., 100. );
	fInvMHistos[3]  = fAC->bookTH1D("MuInvMOSCleanLow", "Mu Clean InvM Low Opp Sign;Mu InvM Opp Sign",  ninvmbins, 0., 5. );
	fInvMHistos[4]  = fAC->bookTH1D("ElInvMSSClean", "El Clean InvM Same Sign;El InvM Same Sign", ninvmbins, 0., 100. );
	fInvMHistos[5]  = fAC->bookTH1D("ElInvMSSCleanLow", "El Clean InvM Low Same Sign;El InvM Same Sign", ninvmbins, 0., 5. );
	fInvMHistos[6]  = fAC->bookTH1D("ElInvMOSClean", "El Clean InvM Opp Sign;El InvM Opp Sign",  ninvmbins, 0., 100. );
	fInvMHistos[7]  = fAC->bookTH1D("ElInvMOSCleanLow", "El Clean InvM Low Opp Sign;El InvM Opp Sign",  ninvmbins, 0., 5. );
	fInvMHistos[8]  = fAC->bookTH1D("PhInvMClean", "Phot Clean InvM;Phot InvM", ninvmbins, 0., 100. );
	fInvMHistos[9]  = fAC->bookTH1D("PhInvMCleanLow", "Phot Clean InvM Low;Phot InvM", ninvmbins, 0., 5. );
	fInvMHistos[10] = fAC->bookTH1D("JInvM2jClean", "Jets Clean InvM dijets;dijet InvM", ninvmbins, 0., 300. );
	fInvMHistos[11] = fAC->bookTH1D("JInvMGT2jClean", "Jets Clean InvM >dijets;GT dijet InvM", ninvmbins, 0., 300. );
}

void PhysQCAnalysis::PrintHisto(TH1D* hist, TString tempstring1, TString tempstring2, TString outputdir, double x1, double x2, int logy, int fract){
	double percent1 = 0.05;
	double percent2 = 0.01;

	TCanvas * canv = new TCanvas(tempstring2, tempstring1 , 0, 0, 900, 700);
	if (logy == 1) gPad->SetLogy();
	hist->DrawCopy();

	double maxy = hist->GetMaximum();
	maxy = 1.05*maxy;
	hist->SetMaximum(maxy);
	if (x1 != 999.) {
		hist->DrawCopy();
		TLine * line = new TLine(x1, 0, x1, maxy);
		line->SetLineColor(kRed);
		line->SetLineWidth(2);
		line->Draw();
	}
	if (x2 != 999.) {
		TLine * line = new TLine(x2, 0, x2, maxy);
		line->SetLineColor(kRed);
		line->SetLineWidth(2);
		line->Draw();
	}
	if (fract == 1) {
		fAC->tailFraction(hist, percent1);
		fAC->tailFraction(hist, percent2);
	}

	Util::Print(canv, tempstring2, outputdir, fHistFile);
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

void PhysQCAnalysis::MakePlots(TString plotlist, TCut select, TTree *tree){
	fAC->plotPlotList(plotlist, tree, "", select, fHistFile);
}

void PhysQCAnalysis::MakeElIDPlots(TCut select, TTree *tree){
	fAC->plotEID(select, tree, "", fHistFile);
}

void PhysQCAnalysis::PlotTriggerStats(){
	const int nentries = fTR->fChain->GetEntries();

	TH1I *hlt_stats, *l1p_stats, *l1t_stats;
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
			TFile *f = TFile::Open(chEl->GetTitle());
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
	gPad->SetLogy();
	hlt_stats->SetName("HLTStats1");
	hlt_stats->GetXaxis()->SetRange(0,50);
	hlt_stats->DrawCopy();
	l1 = new TLine(0, nentries, 50, nentries);
	l1->SetLineColor(kRed);
	l1->SetLineWidth(2);
	l1->SetLineStyle(2);
	l1->Draw();
	fTlat->DrawLatex(0.17,0.92, tempstring);
	fTlat->DrawLatex(0.60,0.92, entries);
	Util::Print(canv, "HLTStats1", fOutputDir+subdir, fHistFile);

	tempstring = "HLT Trigger Bits (51-100)";
	canv = new TCanvas("HLTStats2", tempstring , 0, 0, 900, 700);
	canv->SetBottomMargin(0.50);
	gPad->SetGridy();
	gPad->SetLogy();
	hlt_stats->SetName("HLTStats2");
	hlt_stats->GetXaxis()->SetRange(51,100);
	hlt_stats->DrawCopy();
	l1 = new TLine(50, nentries, 100, nentries);
	l1->SetLineColor(kRed);
	l1->SetLineWidth(2);
	l1->SetLineStyle(2);
	l1->Draw();
	fTlat->DrawLatex(0.17,0.92, tempstring);
	fTlat->DrawLatex(0.60,0.92, entries);	
	Util::Print(canv, "HLTStats2", fOutputDir+subdir, fHistFile);

	tempstring = "L1 Phys Bits (0-63)";
	canv = new TCanvas("L1PStats1", tempstring , 0, 0, 900, 700);
	canv->SetBottomMargin(0.42);
	gPad->SetGridy();
	gPad->SetLogy();
	l1p_stats->SetName("L1PStats1");
	l1p_stats->GetXaxis()->SetRange(0,63);
	l1p_stats->DrawCopy();
	l1 = new TLine(0, nentries, 63, nentries);
	l1->SetLineColor(kRed);
	l1->SetLineWidth(2);
	l1->SetLineStyle(2);
	l1->Draw();
	fTlat->DrawLatex(0.17,0.92, tempstring);
	fTlat->DrawLatex(0.60,0.92, entries);
	Util::Print(canv, "L1PStats1", fOutputDir+subdir, fHistFile);

	tempstring = "L1 Phys Bits (64-128)";
	canv = new TCanvas("L1PStats2", tempstring , 0, 0, 900, 700);
	canv->SetBottomMargin(0.42);
	gPad->SetGridy();
	gPad->SetLogy();
	l1p_stats->SetName("L1PStats2");
	l1p_stats->GetXaxis()->SetRange(64,128);
	l1p_stats->DrawCopy();
	l1 = new TLine(63, nentries, 128, nentries);
	l1->SetLineColor(kRed);
	l1->SetLineWidth(2);
	l1->SetLineStyle(2);
	l1->Draw();
	fTlat->DrawLatex(0.17,0.92, tempstring);
	fTlat->DrawLatex(0.60,0.92, entries);
	Util::Print(canv, "L1PStats2", fOutputDir+subdir, fHistFile);

	tempstring = "L1 Tech Bits";
	canv = new TCanvas("L1TStats", tempstring , 0, 0, 900, 700);
	gPad->SetGridy();
	gPad->SetLogy();
	l1t_stats->SetName("L1TStats");
	l1t_stats->GetXaxis()->SetRange(0,64);
	l1t_stats->DrawCopy();
	l1 = new TLine(0, nentries, 64, nentries);
	l1->SetLineColor(kRed);
	l1->SetLineWidth(2);
	l1->SetLineStyle(2);
	l1->Draw();
	fTlat->DrawLatex(0.17,0.92, tempstring);
	fTlat->DrawLatex(0.60,0.92, entries);
	Util::Print(canv, "L1TStats", fOutputDir+subdir, fHistFile);

	gROOT->cd(); // Leave local file
}
