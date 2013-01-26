#include <vector>
#include <iostream>
#include <map>
#include "TF1.h"
#include "TString.h"


class BTagSF {
	public:
		BTagSF ();
		~BTagSF(){};
		
		double efficiency (float pt, float eta, int flavor, TString meanminmax, int lightUnc = 0);
		double scalefactor(float pt, float eta, int flavor, TString meanminmax);
		void getSFLightFunction(float jeteta, TString meanminmax);
		bool   applySF(bool& isBTagged, float SF, float eff, float random);
		bool   modifyBTagsWithSF(bool& is_tagged, float pt, float eta, int flavor, TString meanminmax, float random); // just to have it in a similar format as before

		TF1 * fLightSF;
		std::map< TString, float> fMeanminmax;
		
};

