#include "helper/Utilities.hh"
#include "ZeeAnalysis.hh"

using namespace std;

ZeeAnalysis::ZeeAnalysis(TreeReader *tr) : UserAnalysisBase(tr){
	Util::SetStyle();	
        elecorr = new EnergyCorrection("electrons");
}

ZeeAnalysis::~ZeeAnalysis(){
  delete elecorr;
}

void ZeeAnalysis::Begin(){
	// Define the output file of histograms
	const char* filename = "histos_ee.root";
	fHistFile = new TFile(fOutputDir + TString(filename), "RECREATE");
	
	// Define the histograms

	fHElPt   = new TH1D("ElPt", "Pt of electrons", 100, 0., 500.);
	fHElPtCorr   = new TH1D("ElPtCorr", "Corrected Pt of electrons", 100, 0., 500.);
	fHInvMass0 = new TH1D("eeInvMass0","Invariant ee mass", 100, 80,100);
	fHInvMass15 = new TH1D("eeInvMass15","Invariant ee mass corr15", 100, 80,100);
	fHInvMass16 = new TH1D("eeInvMassCorr16","Invariant ee mass corr16", 100, 80,100);
	fHInvMass20 = new TH1D("eeInvMassCorr20","Invariant ee mass corr20", 100, 80,100);
}

void ZeeAnalysis::Analyze(){

  for (int i=0; i<fTR->NEles; i++){

    fHElPt -> Fill(fTR->ElPt[i]);
    float correnergy = elecorr->get_correctedenergy(fTR,i,15);
    fHElPtCorr -> Fill(fTR->ElPt[i]*correnergy/fTR->ElE[i]);
 
  }

  if (fTR->NEles<2) return;
  if (fTR->ElPt[1]<20) return;

  bool id=true;
  for (int i=0; i<2; i++){
    if (fTR->ElIDsimpleWP80relIso[i]<7) id=false;
  }
  if (!id) return;

  TLorentzVector elec[2];
  for (int i=0; i<2; i++){
    elec[i].SetPtEtaPhiE(fTR->ElPt[i],fTR->ElEta[i],fTR->ElPhi[i],fTR->ElE[i]);
  }

  bool masswindow=false;
  if (fabs((elec[0]+elec[1]).M()-91.2)<10) masswindow=true;
  if (!masswindow) return;

  fHInvMass0->Fill((elec[0]+elec[1]).M());

  fHInvMass15 ->Fill((CorrElectron(fTR,0,15)+CorrElectron(fTR,1,15)).M());
  fHInvMass16 ->Fill((CorrElectron(fTR,0,16)+CorrElectron(fTR,1,16)).M());
  fHInvMass20 ->Fill((CorrElectron(fTR,0,20)+CorrElectron(fTR,1,20)).M());

}

void ZeeAnalysis::End(){
	fHistFile->cd();	

	fHElPt   ->Write();
	fHElPtCorr -> Write();

	fHInvMass0 -> Write();
	fHInvMass15 -> Write();
	fHInvMass16 -> Write();
	fHInvMass20 -> Write();

	fHistFile->Close();
}

TLorentzVector ZeeAnalysis::CorrElectron(TreeReader *fTR, int i, int mode){
  TLorentzVector corr;
  corr.SetPtEtaPhiE(fTR->ElPt[i],fTR->ElEta[i],fTR->ElPhi[i],fTR->ElE[i]);
  corr.SetE(elecorr->get_correctedenergy(fTR,i,mode));
  corr.SetRho(elecorr->get_correctedenergy(fTR,i,mode));
  return corr;
}
