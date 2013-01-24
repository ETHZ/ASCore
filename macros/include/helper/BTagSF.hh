#include <vector>
#include <iostream>

class BTagSF {
	public:
		BTagSF ();
		~BTagSF(){};
		
		double efficiency (float flavor, double pt, int unc, int lightUnc);
		double scalefactor(float flavor, double pt, int unc, int lightUnc);
		bool   applySF(bool& isBTagged, float SF, float eff, float random);
		bool   modifyBTagsWithSF(bool& is_tagged, int flavor, int pt, float random, int flag); // just to have it in a similar format as before
		
};

