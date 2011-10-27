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
	fHInvMass = new TH1D("eeInvMass","Invariant ee mass", 100, 80,100);
	fHInvMassCorr = new TH1D("eeInvMassCorr","Invariant ee mass corr", 100, 80,100);
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
    if (fTR->ElIDsimpleWP90relIso[i]<7) id=false;
  }
  if (!id) return;

  TLorentzVector elec[2];
  for (int i=0; i<2; i++){
    elec[i].SetPtEtaPhiE(fTR->ElPt[i],fTR->ElEta[i],fTR->ElPhi[i],fTR->ElE[i]);
  }

  bool masswindow=false;
  if (fabs((elec[0]+elec[1]).M()-91.2)<10) masswindow=true;
  if (!masswindow) return;

  TLorentzVector correlec[2];
  for (int i=0; i<2; i++){
    correlec[i]=elec[i];
    correlec[i].SetE(elecorr->get_correctedenergy(fTR,i,15));
    correlec[i].SetRho(elecorr->get_correctedenergy(fTR,i,15));
  }

  fHInvMass->Fill((elec[0]+elec[1]).M());
  fHInvMassCorr ->Fill((correlec[0]+correlec[1]).M());


}

void ZeeAnalysis::End(){
	fHistFile->cd();	
	fHElPt   ->Write();
	fHElPtCorr -> Write();
	fHInvMass -> Write();
	fHInvMassCorr -> Write();
	fHistFile->Close();
}
