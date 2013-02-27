#include "helper/BTagSF.hh"
#include <math.h>
#include <iostream>
#include <TString.h>

using namespace std;

BTagSF::BTagSF() {
	fMeanminmax["mean"] =  0.;
	fMeanminmax["min" ] = -1.;
	fMeanminmax["max" ] =  1.;

	float ptmax = 800;
	// eta between 0 and 0.8
	fLightSFeta0mean  = new TF1("SFlight"   , "((1.06238+(0.00198635*x))+(-4.89082e-06*(x*x)))+(3.29312e-09*(x*(x*x)))" , 20., ptmax);
	fLightSFeta0min   = new TF1("SFlightMin", "((0.972746+(0.00104424*x))+(-2.36081e-06*(x*x)))+(1.53438e-09*(x*(x*x)))", 20., ptmax);
	fLightSFeta0max   = new TF1("SFlightMax", "((1.15201+(0.00292575*x))+(-7.41497e-06*(x*x)))+(5.0512e-09*(x*(x*x)))"  , 20., ptmax);
	// eta between 0.8 and 1.6
	fLightSFeta1mean  = new TF1("SFlight"   , "((1.08048+(0.00110831*x))+(-2.96189e-06*(x*x)))+(2.16266e-09*(x*(x*x)))" , 20., ptmax);
	fLightSFeta1min   = new TF1("SFlightMin", "((0.9836+(0.000649761*x))+(-1.59773e-06*(x*x)))+(1.14324e-09*(x*(x*x)))" , 20., ptmax);
	fLightSFeta1max   = new TF1("SFlightMax", "((1.17735+(0.00156533*x))+(-4.32257e-06*(x*x)))+(3.18197e-09*(x*(x*x)))" , 20., ptmax);
	// eta between 1.6 and 2.4
	fLightSFeta2mean  = new TF1("SFlight"   , "((1.09145+(0.000687171*x))+(-2.45054e-06*(x*x)))+(1.7844e-09*(x*(x*x)))" , 20., ptmax);
	fLightSFeta2min   = new TF1("SFlightMin", "((1.00616+(0.000358884*x))+(-1.23768e-06*(x*x)))+(6.86678e-10*(x*(x*x)))", 20., ptmax);
	fLightSFeta2max   = new TF1("SFlightMax", "((1.17671+(0.0010147*x))+(-3.66269e-06*(x*x)))+(2.88425e-09*(x*(x*x)))"  , 20., ptmax);
};

double BTagSF::efficiency(float jetpt, float jeteta, int flavor, TString meanminmax, int uncertaintyLight){
	if      (fabs(flavor)==5) return 0.7195 +fMeanminmax[meanminmax]*(0.7665-0.7195);
	else if (fabs(flavor)==4) return 0.19249+fMeanminmax[meanminmax]*(0.7665 - 0.7195)/5;
	else                      return (0.0113428+(5.18983e-05*jetpt))+(-2.59881e-08*(jetpt*jetpt))*(1+uncertaintyLight*0.5);
}

float BTagSF::getSFLight(float jetpt, float jeteta, TString meanminmax) {
	// assuming tagger is always CSV medium. might want to extend that in the future
	// float ptmax = 800.;
	float eta = fabs(jeteta); // making sure we're taking the absolute eta
	float sf(-999.);
	if ( eta >= 0.  && eta <= 0.8 ) {
		if (meanminmax == "mean" ) sf = fLightSFeta0mean ->Eval(jetpt);
		if (meanminmax == "min"  ) sf = fLightSFeta0min  ->Eval(jetpt);
		if (meanminmax == "max"  ) sf = fLightSFeta0max  ->Eval(jetpt);
	}
	else if ( eta  > 0.8 && eta <= 1.6 ) {
		if (meanminmax == "mean" ) sf = fLightSFeta1mean ->Eval(jetpt);
		if (meanminmax == "min"  ) sf = fLightSFeta1min  ->Eval(jetpt);
		if (meanminmax == "max"  ) sf = fLightSFeta1max  ->Eval(jetpt);
	}
	else if ( eta  > 1.6 && eta <= 2.4 ) {
		if (meanminmax == "mean" ) sf = fLightSFeta2mean ->Eval(jetpt);
		if (meanminmax == "min"  ) sf = fLightSFeta2min  ->Eval(jetpt);
		if (meanminmax == "max"  ) sf = fLightSFeta2max  ->Eval(jetpt);
	}
	return sf;
}

double BTagSF::scalefactor(float jetpt, float jeteta, int flavor, TString meanminmax){

	float ptmin[] = {20, 30, 40, 50, 60, 70, 80, 100, 120, 160, 210, 260, 320, 400, 500, 600};
	float ptmax[] = {30, 40, 50, 60, 70, 80,100, 120, 160, 210, 260, 320, 400, 500, 600, 800};
	int binnumber=0;
	int nbins = sizeof(ptmin)/sizeof(int)-1;
	bool tooHigh = false;
	if ( jetpt>ptmax[nbins] ) tooHigh = true;
	else {
		for (int i=0; i<nbins; ++i) {
			if(jetpt>ptmin[i] && jetpt<ptmax[i]) {
				binnumber=i;
				break;
			}
		}
	}

	// B-JETS
	if (fabs(flavor)==5) {
		double SFb = 0.726981*((1.+(0.253238*jetpt))/(1.+(0.188389*jetpt)));
		double SFb_error[] = {
		 0.0554504,
		 0.0209663,
		 0.0207019,
		 0.0230073,
		 0.0208719,
		 0.0200453,
		 0.0264232,
		 0.0240102,
		 0.0229375,
		 0.0184615,
		 0.0216242,
		 0.0248119,
		 0.0465748,
		 0.0474666,
		 0.0718173,
		 0.0717567 };
		if ( tooHigh ) return SFb+2*fMeanminmax[meanminmax]*SFb_error[nbins];
		return SFb+fMeanminmax[meanminmax]*SFb_error[binnumber];
	}
	// C-JETS
	if (fabs(flavor)==4) {
		double SFb = 0.726981*((1.+(0.253238*jetpt))/(1.+(0.188389*jetpt)));
		double SFb_error[] = {
		 0.0554504,
		 0.0209663,
		 0.0207019,
		 0.0230073,
		 0.0208719,
		 0.0200453,
		 0.0264232,
		 0.0240102,
		 0.0229375,
		 0.0184615,
		 0.0216242,
		 0.0248119,
		 0.0465748,
		 0.0474666,
		 0.0718173,
		 0.0717567 };
		if ( tooHigh ) return SFb+2*2*fMeanminmax[meanminmax]*SFb_error[nbins];
		return SFb+fMeanminmax[meanminmax]*2*SFb_error[binnumber];
	}
	// LIGHT JETS
	else {
		// if (fLightSF == NULL ) cout << Form("jetpt: %5.2f jeteta: %3.2f jetflavor: %d ", jetpt, jeteta, flavor) << endl;
		float SFlight = getSFLight(jetpt, jeteta, meanminmax);
		// delete fLightSF;
		return SFlight; // return only the central value for now. should adapt for jets > 800 GeV etc.
	}
}

bool BTagSF::applySF(bool& isBTagged, float Btag_SF, float Btag_eff, float random){
	bool newBTag = isBTagged;
	if (Btag_SF == 1) return newBTag; //no correction needed 
	if (Btag_SF > 1){                                                             // use this if SF>1
		if( !isBTagged ) {
			float mistagPercent = (1.0 - Btag_SF) / (1.0 - (Btag_SF/Btag_eff) );  //fraction of jets that need to be upgraded      
			if( random < mistagPercent ) newBTag = true;                          //upgrade to tagged
		}
	}
	else {                                                                        // use this if SF<1
	      if( isBTagged && random > Btag_SF ) newBTag = false;                    //downgrade tagged to untagged
	
	}
	return newBTag;
}

bool BTagSF::modifyBTagsWithSF(bool& is_tagged, float pt, float eta, int flavor, TString meanminmax, float random) {
	if (eta < -2.399 ) eta = -2.399;
	if (eta >  2.399 ) eta =  2.399;
	double btageff = efficiency (pt, eta, flavor, meanminmax);
	double btagSF  = scalefactor(pt, eta, flavor, meanminmax);
	return applySF(is_tagged, btagSF, btageff, random);                       ///--->> Apply scale factor
}
