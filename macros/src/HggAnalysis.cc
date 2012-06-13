#include "helper/Utilities.hh"
#include "HggAnalysis.hh"

using namespace std;

HggAnalysis::HggAnalysis(TreeReader *tr) : UserAnalysisBase(tr){
	Util::SetStyle();	
        phocorr = new EnergyCorrection("photons");
}

HggAnalysis::~HggAnalysis(){
  delete phocorr;
}

void HggAnalysis::Begin(){
	// Define the output file of histograms
	const char* filename = "histos_gg.root";
	fHistFile = new TFile(fOutputDir + TString(filename), "RECREATE");
	
	// Define the histograms

	fHPhoPt   = new TH1D("PhoPt", "Pt of photons", 100, 0., 500.);
	fHPhoPtCorr   = new TH1D("PhoPtCorr", "Corrected Pt of photons", 100, 0., 500.);
	fHInvMass = new TH1D("ggInvMass","Invariant gg mass", 100, 100,140);
	fHInvMassCorr = new TH1D("ggInvMassCorr","Invariant gg mass corr", 100, 100,140);
}

void HggAnalysis::Analyze(){

  for (int i=0; i<fTR->NPhotons; i++){
    fHPhoPt -> Fill(fTR->PhoPt[i]);
    float correnergy = phocorr->get_correctedenergy(fTR,i,5);
    fHPhoPtCorr -> Fill(fTR->PhoPt[i]*correnergy/fTR->PhoEnergy[i]);
   }


  if (fTR->NPhotons<2) return;
  if (fTR->PhoPt[1]<20) return;

  TLorentzVector phot[2];
  for (int i=0; i<2; i++){
    phot[i].SetPtEtaPhiE(fTR->PhoPt[i],fTR->PhoEta[i],fTR->PhoPhi[i],fTR->PhoEnergy[i]);
  }

  TLorentzVector corrphot[2];
  for (int i=0; i<2; i++){
    corrphot[i]=phot[i];
    corrphot[i].SetE(phocorr->get_correctedenergy(fTR,i,5));
    corrphot[i].SetRho(phocorr->get_correctedenergy(fTR,i,5));
  }

  fHInvMass->Fill((phot[0]+phot[1]).M());
  fHInvMassCorr ->Fill((corrphot[0]+corrphot[1]).M());


}

void HggAnalysis::End(){
	fHistFile->cd();	
	fHPhoPt   ->Write();
	fHPhoPtCorr -> Write();
	fHInvMass -> Write();
	fHInvMassCorr -> Write();
	fHistFile->Close();
}
