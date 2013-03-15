#include "helper/OnTheFlyCorrections.hh"
#include <math.h>

using namespace std;

OnTheFlyCorrections::OnTheFlyCorrections(std::string gt, bool isdata){
	std::string path="/shome/mdunser/jetfiles/";
	JetCorrectorParameters *ResJetPar = new JetCorrectorParameters(path+gt+"_L2L3Residual_AK5PF.txt");
	JetCorrectorParameters *L3JetPar  = new JetCorrectorParameters(path+gt+"_L3Absolute_AK5PF.txt");
	JetCorrectorParameters *L2JetPar  = new JetCorrectorParameters(path+gt+"_L2Relative_AK5PF.txt");
	JetCorrectorParameters *L1JetPar  = new JetCorrectorParameters(path+gt+"_L1FastJet_AK5PF.txt");
	fJetCorPar.push_back(*L1JetPar);
	fJetCorPar.push_back(*L2JetPar);
	fJetCorPar.push_back(*L3JetPar);
	if (isdata) fJetCorPar.push_back(*ResJetPar);
	fJetCorrector = new FactorizedJetCorrector(fJetCorPar);
	fIsData = isdata;
	
	delete L1JetPar;
	delete L2JetPar;
	delete L3JetPar; 
	delete ResJetPar; 

	fJetUncertainty = new JetCorrectionUncertainty(path+gt+"_Uncertainty_AK5PF.txt");
}

float OnTheFlyCorrections::getJECUncertainty(float pt, float eta){
	if      (eta> 5.0) eta = 5.0;
	else if (eta<-5.0) eta =-5.0;
	fJetUncertainty->setJetPt(pt);
	fJetUncertainty->setJetEta(eta);
	float uncert= fJetUncertainty->getUncertainty(true);
	return uncert;

}

std::vector< float > OnTheFlyCorrections::getCorrPtECorr(float jetpt, float jeteta, float jetenergy, float jecorr, float jetarea, float rho){
// give this function a corrected jet and it will return a vector with the 
// corrected jet-pt, jet-energy and the new correction according to the global tag that is set

	float rawpt = jetpt/jecorr;
	float rawe  = jetenergy/jecorr;
	fJetCorrector->setJetPt(rawpt); // raw-pT here
	fJetCorrector->setJetEta(jeteta);
	fJetCorrector->setRho(rho);
	fJetCorrector->setJetA(jetarea);

	float corr = fJetCorrector->getCorrection();
	std::vector< float > corrected;

	// new jetpt = old jet-pt * new correction / old correction
	corrected.push_back(corr*rawpt); // new corrected pt as 0th item
	corrected.push_back(corr*rawe);  // new corrected energy as 1st item
	corrected.push_back(corr);       // new correction as 2nd item
	return corrected;

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

float OnTheFlyCorrections::getJetPtNoResidual(float pt, float eta, float ecorr, float area, float rho){
	if (!fIsData) return pt; // no residual corrections for MC
	float rawpt = pt/ecorr;
	fJetCorrector->setJetEta(eta);
	fJetCorrector->setRho(rho);
	fJetCorrector->setJetA(area);
	fJetCorrector->setJetPt(rawpt);      // the jets we have in the collection are raw already
	
	vector<float> factors = fJetCorrector->getSubCorrections();
	if(fabs(factors[3] - ecorr) > 0.000001) cout << "OnTheFlyCorrections::getJetPtNoResidual ==> WARNING: Your JECs don't seem to be consistent!" << endl;
	float l1l2l3scale = factors[2];
	return rawpt*l1l2l3scale;
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
