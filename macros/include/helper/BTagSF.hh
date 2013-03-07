#include <vector>
#include <iostream>
#include <map>
#include "TF1.h"
#include "TString.h"


class BTagSF {
	public:
		BTagSF ();
		~BTagSF(){};
		
		struct SFandErr {
			float sf;
			float err;
		};

		double efficiency (float pt, float eta, int flavor, TString meanminmax, int lightUnc = 0);
		double scalefactor(float pt, float eta, int flavor, TString meanminmax, bool isFastsim=false, TString model="");
		float getSFLight(float jetpt, float jeteta, TString meanminmax);
		bool   applySF(bool& isBTagged, float SF, float eff, float random);
		bool   modifyBTagsWithSF(bool&, float, float, int, TString, float, bool = false, TString = ""); // just to have it in a similar format as before
		SFandErr getSMSSFandError(float, float, int, TString);

		TF1 * fLightSF;
		std::map< TString, float> fMeanminmax;

		TF1 * fLightSFeta0mean ;
		TF1 * fLightSFeta0min  ;
		TF1 * fLightSFeta0max  ;
		TF1 * fLightSFeta1mean ;
		TF1 * fLightSFeta1min  ;
		TF1 * fLightSFeta1max  ;
		TF1 * fLightSFeta2mean ;
		TF1 * fLightSFeta2min  ;
		TF1 * fLightSFeta2max  ;
		
};

