#include "helper/OnTheFlyMet.hh"
#include <math.h>

using namespace std;

OnTheFlyCorrections::OnTheFlyCorrections(string gt, bool isdata){
	JetCorrectorParameters *L1JetPar  = new JetCorrectorParameters("JetCorrectionFiles/"+gt+"_L1FastJet_AK5PF.txt");
	JetCorrectorParameters *L2JetPar  = new JetCorrectorParameters("JetCorrectionFiles/"+gt+"_L2Relative_AK5PF.txt");
	JetCorrectorParameters *L3JetPar  = new JetCorrectorParameters("JetCorrectionFiles/"+gt+"_L3Absolute_AK5PF.txt");
	JetCorrectorParameters *ResJetPar = new JetCorrectorParameters("JetCorrectionFiles/"+gt+"_L2L3Residual_AK5PF.txt");
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
