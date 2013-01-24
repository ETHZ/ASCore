#include "helper/BTagSF.hh"
#include <math.h>
#include <iostream>

BTagSF::BTagSF() {};

double BTagSF::efficiency(float flavor, double jetpt, int uncertainty, int uncertaintyLight){
	if      (fabs(flavor)==5) return 0.7195+uncertainty*(0.7665-0.7195);
	else if (fabs(flavor)==4) return 0.19249+uncertainty*(0.7665 - 0.7195)/5;
	else                      return (0.0113428+(5.18983e-05*jetpt))+(-2.59881e-08*(jetpt*jetpt))*(1+uncertaintyLight*0.5);
}

double BTagSF::scalefactor(float flavor, double jetpt, int uncertainty, int uncertaintyLight){

	float ptmin[] = {30, 40, 50, 60, 70, 80, 100, 120, 160, 210, 260, 320, 400, 500, 600};
	int binnumber=0;
	for (int i=0;i<sizeof(ptmin)/sizeof(int)-1;i++){
		if(jetpt>ptmin[i] && jetpt<ptmin[i+1]) {
			binnumber=i;
			break;
		}
	}

	// B-JETS
	if (fabs(flavor)==5) {
		double SFb = 0.6981*((1.+(0.414063*jetpt))/(1.+(0.300155*jetpt)));
		double SFb_error[] = {
			0.0295675,
			0.0295095,
			0.0210867,
			0.0219349,
			0.0227033,
			0.0204062,
			0.0185857,
			0.0256242,
			0.0383341,
			0.0409675,
			0.0420284,
			0.0541299,
			0.0578761,
			0.0655432 };
		return SFb+uncertainty*1.5*SFb_error[binnumber];
	}
	// C-JETS
	if (fabs(flavor)==4) {
		double SFb = 0.6981*((1.+(0.414063*jetpt))/(1.+(0.300155*jetpt)));
		double SFb_error[] = {
			0.0295675,
			0.0295095,
			0.0210867,
			0.0219349,
			0.0227033,
			0.0204062,
			0.0185857,
			0.0256242,
			0.0383341,
			0.0409675,
			0.0420284,
			0.0541299,
			0.0578761,
			0.0655432 };
		return SFb+uncertainty*1.5*2*SFb_error[binnumber];   ///---? Where this come from 
	}
	// LIGHT JETS
	else {
		double SF2012=1.10422 + -0.000523856*jetpt + 1.14251e-06*jetpt*jetpt;
		if(uncertaintyLight==-1){
			double SFmin = ((0.962627+(0.000448344*jetpt))+(-1.25579e-06*(jetpt*jetpt)))+(4.82283e-10*(jetpt*(jetpt*jetpt)));
			return SFmin*SF2012;
		}
		if(uncertaintyLight==1){
			double SFmax = ((1.12368+(0.00124806*jetpt))+(-3.9032e-06*(jetpt*jetpt)))+(2.80083e-09*(jetpt*(jetpt*jetpt)));
			return SFmax*SF2012;
		}
		else{
			double SFmean = ((1.04318+(0.000848162*jetpt))+(-2.5795e-06*(jetpt*jetpt)))+(1.64156e-09*(jetpt*(jetpt*jetpt)));
			return SFmean*SF2012;
		}
	}

}

bool BTagSF::applySF(bool& isBTagged, float Btag_SF, float Btag_eff, float random){
	bool newBTag = isBTagged;
	if (Btag_SF == 1) return newBTag; //no correction needed 
	if (Btag_SF > 1){                                                            // use this if SF>1
		if( !isBTagged ) {
			float mistagPercent = (1.0 - Btag_SF) / (1.0 - (Btag_SF/Btag_eff) );  //fraction of jets that need to be upgraded      
			if( random < mistagPercent ) newBTag = true;                          //upgrade to tagged
		}
	}
	else {                                                                    // use this if SF<1
	      if( isBTagged && random > Btag_SF ) newBTag = false;                   //downgrade tagged to untagged
	
	}
	return newBTag;
}

bool BTagSF::modifyBTagsWithSF(bool& is_tagged, int flavor, int pt, float random, int flag) {
	int factor(0);
	if      (flag == 1) factor =  1;
	else if (flag == 2) factor = -1;
	double btageff = efficiency (flavor, pt, factor, 0);
	double btagSF  = scalefactor(flavor, pt, factor, 0);
	return applySF(is_tagged, btagSF, btageff, random);                       ///--->> Apply scale factor
}
