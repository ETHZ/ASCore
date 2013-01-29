#include <vector>
#include <iostream>
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "../base/TreeReader.hh"

class OnTheFlyCorrections {
	public:
		OnTheFlyCorrections (std::string globaltag, bool isdata);
		~OnTheFlyCorrections(){};

		std::vector<JetCorrectorParameters> fJetCorPar;
		FactorizedJetCorrector *fJetCorrector;
		// FactorizedJetCorrector(OnTheFlyCorrections::fJetCorPar) *fJetCorrector;

		std::pair< float, float > getCorrections( float rawpt, float raweta, float rawnomupt, 
                                                  float phi, float emf, float rho, float area, std::string level = "" ); // for on the fly corrections
		float getJetCorrection     (float pt, float corr, float eta, float rho, float area, std::string level );     // this function returns, for a given jet the correction factor
		float getJetCorrectionRawPt(float pt,             float eta, float rho, float area, std::string level );     // same as above, for people who want to call it with the raw-pt already
		float getJetCorrectedPt    (float pt, float corr, float eta, float rho, float area                    );     // returns really the new jet-pT with the maximum corrections applied

		float getPx(float pt, float phi){ return pt*cos(phi); };
		float getPy(float pt, float phi){ return pt*sin(phi); };
};

