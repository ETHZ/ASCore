#include "helper/OnTheFlyCorrections.hh"
#include <math.h>

using namespace std;

OnTheFlyCorrections::OnTheFlyCorrections(string gt, bool isdata){
	JetCorrectorParameters *ResJetPar = new JetCorrectorParameters("src/helper/JetCorrectionFiles/"+gt+"_L2L3Residual_AK5PF.txt");
	JetCorrectorParameters *L3JetPar  = new JetCorrectorParameters("src/helper/JetCorrectionFiles/"+gt+"_L3Absolute_AK5PF.txt");
	JetCorrectorParameters *L2JetPar  = new JetCorrectorParameters("src/helper/JetCorrectionFiles/"+gt+"_L2Relative_AK5PF.txt");
	JetCorrectorParameters *L1JetPar  = new JetCorrectorParameters("src/helper/JetCorrectionFiles/"+gt+"_L1FastJet_AK5PF.txt");
	fJetCorPar.push_back(*L1JetPar);
	fJetCorPar.push_back(*L2JetPar);
	fJetCorPar.push_back(*L3JetPar);
	if (isdata) fJetCorPar.push_back(*ResJetPar);
	fJetCorrector = new FactorizedJetCorrector(fJetCorPar);
	
	delete L1JetPar;
	delete L2JetPar;
	delete L3JetPar; 
	delete ResJetPar; 
}

float OnTheFlyCorrections::getJetCorrection(float pt, float corr, float eta, float rho, float area, string level = "L3Absolute"){
	// sets the pT back to raw and returns the raw-pT correction factor
	return getJetCorrectionRawPt(pt/corr, eta, rho, area, level);
}

float OnTheFlyCorrections::getJetCorrectionRawPt(float rawpt, float eta, float rho, float area, string level = "L3Absolute"){
	// slighly redundant considering we have what we have below, but I think that's what frederic was thinking about
	fJetCorrector->setJetEta(eta);
	fJetCorrector->setRho(rho);
	fJetCorrector->setJetA(area);
	fJetCorrector->setJetPt(rawpt); // new raw-pT here...! this is called with fTR->JPt[jetindex]/fTR->JEcorr[jetindex]; in the useranalysisbase.
	std::vector< float > corrections = fJetCorrector->getSubCorrections();

	if (level == "L1FastJet")    return corrections.front();
	if (level == "L2Relative")   return corrections[1];
	if (level == "L2L3Residual") return corrections.back();
	return corrections[2];
}

float OnTheFlyCorrections::getJetCorrectedPt(float pt, float corr, float eta, float rho, float area){
	// slighly redundant considering we have what we have below, but I think that's what frederic was thinking about
	fJetCorrector->setJetEta(eta);
	fJetCorrector->setRho(rho);
	fJetCorrector->setJetA(area);
	fJetCorrector->setJetPt(pt/corr);
	std::vector< float > corrections = fJetCorrector->getSubCorrections();

	return corrections.back() * pt/corr;

}

std::pair< float, float > 
OnTheFlyCorrections::getCorrections(float rawpt, float raweta, float rawnomupt, float phi, float emf, float rho, float area, string level) {
  
  std::pair< float, float > corr(0., 0.); // pair with zeroes that gets returned if jet fails a cut
  if (emf > 0.9)       return corr;       // skip jets with EMF > 0.9

  fJetCorrector->setJetEta(raweta);
  fJetCorrector->setRho(rho);
  fJetCorrector->setJetA(area);
  fJetCorrector->setJetPt(rawpt);      // the jets we have in the collection are raw already
  
  std::vector< float > corrections = fJetCorrector->getSubCorrections();
  

  float l1corrpt   = rawnomupt*corrections.front(); // l1fastjet corrections were pushed pack first
  float fullcorrpt = rawnomupt*corrections.back();  // full corrections are the last in the vector

  if (fullcorrpt < 10.)     return corr;        // skip all jets that have corrected pt below 10 GeV

  // the corretions for the MET are the difference between l1fastjet and the full corrections on the jet!
  
  corr.first  = getPx(l1corrpt - fullcorrpt, phi); // fill the px of the correction in the pair.first
  corr.second = getPy(l1corrpt - fullcorrpt, phi); // and the py in the pair.second
  
  return corr;
}
