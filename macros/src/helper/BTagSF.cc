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
	if      (fabs(flavor)==5) return 0.7195 +fMeanminmax[meanminmax]*(0.7665 - 0.7195);
	else if (fabs(flavor)==4) return 0.19249+fMeanminmax[meanminmax]*(0.7665 - 0.7195)/5;
	else                      return (0.0113428+(5.18983e-05*jetpt))+(-2.59881e-08*(jetpt*jetpt))*(1+uncertaintyLight*0.5);
}

float BTagSF::getSFLight(float jetpt, float jeteta, TString meanminmax) {
	// assuming tagger is always CSV medium. might want to extend that in the future
	// float ptmax = 800.;
	if (jetpt > 800.) jetpt = 800.;
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

double BTagSF::scalefactor(float jetpt, float jeteta, int flavor, TString meanminmax, bool isFastsim, TString model){

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
	BTagSF::SFandErr fsFactors;
	if (isFastsim){
		fsFactors = getSMSSFandError(jetpt, jeteta, flavor, model);
	}

	// B-JETS
	if (fabs(flavor)==5) {
		double SFb = 0.726981*((1.+(0.253238*jetpt))/(1.+(0.188389*jetpt)));
		double SFb_error[] = { 0.0554504, 0.0209663, 0.0207019, 0.0230073, 0.0208719, 0.0200453, 0.0264232, 0.0240102,
		                       0.0229375, 0.0184615, 0.0216242, 0.0248119, 0.0465748, 0.0474666, 0.0718173, 0.0717567 };
		if (isFastsim){
			if (tooHigh) return SFb*fsFactors.sf + fMeanminmax[meanminmax]*sqrt(4*SFb_error[nbins]*SFb_error[nbins] + fsFactors.err*fsFactors.err);
			return SFb*fsFactors.sf + fMeanminmax[meanminmax]*sqrt(SFb_error[nbins]*SFb_error[nbins] + fsFactors.err*fsFactors.err);
		}
		if ( tooHigh ) return SFb+2*fMeanminmax[meanminmax]*SFb_error[nbins]; // for jets with pT > 800, use twice the uncertainty on the last SF
		return SFb+fMeanminmax[meanminmax]*SFb_error[binnumber];
	}
	// C-JETS
	if (fabs(flavor)==4) {
		double SFb = 0.726981*((1.+(0.253238*jetpt))/(1.+(0.188389*jetpt)));
		double SFb_error[] = { 0.0554504, 0.0209663, 0.0207019, 0.0230073, 0.0208719, 0.0200453, 0.0264232, 0.0240102, 
                               0.0229375, 0.0184615, 0.0216242, 0.0248119, 0.0465748, 0.0474666, 0.0718173, 0.0717567 };
		if (isFastsim){
			if (tooHigh) return SFb*fsFactors.sf + fMeanminmax[meanminmax]*sqrt(4*SFb_error[nbins]*4*SFb_error[nbins] + fsFactors.err*fsFactors.err);
			return SFb*fsFactors.sf + fMeanminmax[meanminmax]*sqrt(2*SFb_error[nbins]*2*SFb_error[nbins] + fsFactors.err*fsFactors.err);
		}
		if ( tooHigh ) return SFb+2*2*fMeanminmax[meanminmax]*SFb_error[nbins];
		return SFb+fMeanminmax[meanminmax]*2*SFb_error[binnumber]; // use scale factors from bjets with twice the uncertainty
	}
	// LIGHT JETS
	else {
		// if (fLightSF == NULL ) cout << Form("jetpt: %5.2f jeteta: %3.2f jetflavor: %d ", jetpt, jeteta, flavor) << endl;
		float SFlight = getSFLight(jetpt, jeteta, meanminmax);
		if (isFastsim){
			return SFlight*fsFactors.sf + fMeanminmax[meanminmax]*sqrt(fsFactors.err*fsFactors.err); // not 100% correct, but doesn't make any difference
		}
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

bool BTagSF::modifyBTagsWithSF(bool& is_tagged, float pt, float eta, int flavor, TString meanminmax, float random, bool isFastsim, TString model) {
	if (eta < -2.39999 ) eta = -2.39999; // protection against one buggy event in WW sample
	if (eta >  2.39999 ) eta =  2.39999;
	double btageff = efficiency (pt, eta, flavor, meanminmax);
	double btagSF  = scalefactor(pt, eta, flavor, meanminmax, isFastsim, model);
	return applySF(is_tagged, btagSF, btageff, random);                       ///--->> Apply scale factor
}

//==========================================================
// ALL THE FASTSIM SCALE FACTORS AND ERRORS WILL FOLLOW HERE
//==========================================================

BTagSF::SFandErr BTagSF::getSMSSFandError(float jetpt, float jeteta, int jetflavor, TString model){

	// i hereby apologize to anyone who ever has to read this code...
	// but i was somewhat constrained by the format people provide these
	// b tagging scale factors in... why not use TF1s?? what's so hard?

	// make the vector for the scale factors
	std::vector<float> vecSF;
	// make the vector for the PU and stat uncertainties
	std::vector<float> vecErr;
	// make the vector for the model uncertainties
	std::vector<float> vec;

	if (fabs(jetflavor) == 5){ // FOR B-JETS
		float SF[]         = {0.982194,0.980998,0.992014,0.994472,0.996825,0.999822,1.00105,1.00023,0.991994,0.979123,0.947207,0.928006,0.874260,0.839610};
		float SF_error[]   = {0.00253112,0.00296453,0.00113963,0.00128363,0.00232566,0.00232353,0.00219086,0.00156856,0.00322279,0.00400414,0.00737465,0.0105033,0.0171706,0.0344172}; // stat + PU
		vecSF .assign(SF, SF + sizeof(SF) / sizeof(SF[0]));
		vecErr.assign(SF_error, SF_error + sizeof(SF_error) / sizeof(SF_error[0]));
		if (model == "T1") {
			float err[] = {0.0305974,0.0251517,0.0205015,0.0187029,0.0138344,0.0155380,0.0153906,0.0210581,0.0175900,-0.00234255,0.000241935,-0.0287645,-0.0472476,-0.0841584};
			vec.assign(err, err + sizeof(err) / sizeof(err[0]));
		}
		if (model == "T1bbbb") {
			float err[] = {0.00897405,0.00984249,0.00694051,0.00454724,0.00505632,0.00173861,0.00184828,0.00124377,-0.00265479,-0.0100402,-0.0112412,-0.0261436,-0.0221387,-0.0377308};
			vec.assign(err, err + sizeof(err) / sizeof(err[0]));
		}
		if (model == "T1tttt") {
			float err[] = {0.0112096,0.0127103,0.0107696,0.0105987,0.0102283,0.00953639,0.0107003,0.0118546,0.00837368,0.000790179,-0.00111371,-0.0146178,-0.00818416,-0.0197257};
			vec.assign(err, err + sizeof(err) / sizeof(err[0]));
		}
		if (model == "T2")  { 
			float err[] = {0.0197094,0.0218538,0.00671038,0.00481349,0.00234514,-0.00960910,-0.00872135,-0.0109075,-0.0185559,-0.0352550,-0.0374648,-0.0604555,-0.0752100,-0.0999645};
			vec.assign(err, err + sizeof(err) / sizeof(err[0]));
		}
		if (model == "T2bb")   {
			float err[] = {0.0125569,0.0119411,0.0100657,0.0106521,0.00982046,0.00745928,0.00802320,0.00942034,0.00741357,0.00160137,0.00219074,-0.00892913,0.00172952,-0.000213087};
			vec.assign(err, err + sizeof(err) / sizeof(err[0]));
		}
		if (model == "T2bw")   {
			float err[] = {0.0111744,0.0112791,0.00760594,0.00597137,0.00484192,0.00301468,0.00359970,0.00540084,0.00215334,-0.00427964,-0.00468144,-0.0184798,-0.0110016,-0.0187086};
			vec.assign(err, err + sizeof(err) / sizeof(err[0]));
		}
		if (model == "T2tt")   {
			float err[] = {0.00574604,0.00677246,0.00509557,0.00374240,0.00314873,0.000637591,-0.000242591,-4.16636e-05,-0.00292352,-0.00581479,-0.000461876,-0.00676391,0.00488830,3.05474e-05};
			vec.assign(err, err + sizeof(err) / sizeof(err[0]));
		}
		if (model == "T3w")    {
			float err[] = {0.0197131,0.0171196,0.0159192,0.0127636,0.0132435,0.00963777,0.00937313,0.00896174,0.00418186,-0.00353286,-0.00389037,-0.0171415,-0.0120094,-0.0215860};
			vec.assign(err, err + sizeof(err) / sizeof(err[0]));
		}
	} // end bjets

	else if (fabs(jetflavor) == 4){ // FOR CHARM JETS
		float SF[]         = {0.988545,0.981714,1.00946,1.01591,1.02810,1.02195,1.02590,1.01936,0.991228,0.955343,0.944433,0.917282,0.935018,1.06375};
		float SF_error[]   = {0.00746259,0.00661831,0.00968682,0.00751322,0.00675507,0.00562821,0.00862890,0.00768003,0.0188981,0.0261163,0.0450601,0.0448453,0.148805,0.177157}; // stat + PU
		vecSF .assign(SF, SF + sizeof(SF) / sizeof(SF[0]));
		vecErr.assign(SF_error, SF_error + sizeof(SF_error) / sizeof(SF_error[0]));
		if (model == "T1")      {
			float err[] = {-0.00104785,-0.00216204,0.00839430,0.00884008,0.00204188,0.00832790,0.0169768,0.0229713,0.0157189,-0.00730190,-0.0692086,-0.108517,-0.137035,-0.181932};
			vec.assign(err, err + sizeof(err) / sizeof(err[0]));
		}
		if (model == "T2bbbb")  { 
			float err[] = {-0.0169490,-0.0109324,-0.0173578,-0.0226300,-0.0354031,-0.0380664,-0.0406916,-0.0448566,-0.0634652,-0.0916214,-0.142743,-0.168372,-0.179460,-0.223442};
			vec.assign(err, err + sizeof(err) / sizeof(err[0]));
		}
		if (model == "T1tttt")  {
			float err[] = {-0.00769350,0.00246567,0.00672805,0.00625175,0.0121922,0.0183616,0.0224260,0.0350031,0.0361672,0.0372230,0.0116431,0.0207569,0.0382855,0.0252644};
			vec.assign(err, err + sizeof(err) / sizeof(err[0]));
		}
		if (model == "T2")      {
			float err[] = {-0.0147888,-0.00520468,-0.00901467,-0.0194454,-0.00635600,-0.00759417,-0.00953454,-0.0174082,-0.0184701,-0.0257653,-0.0740010,-0.0899951,-0.0860117,-0.0738075};
			vec.assign(err, err + sizeof(err) / sizeof(err[0]));
		}
		if (model == "T2bb")    {
			float err[] = {-0.0126963,0.00240847,-0.0237588,-0.0202803,-0.0362858,-0.0296324,-0.0417801,-0.0566426,-0.0675621,-0.0768022,-0.141505,-0.160204,-0.199828,-0.237504};
			vec.assign(err, err + sizeof(err) / sizeof(err[0]));
		}
		if (model == "T2bw")    {
			float err[] = {-0.0177944,-0.0168491,-0.0145971,-0.0171311,-0.0170042,-0.0143744,-0.0160470,-0.0149559,-0.0172561,-0.0137762,-0.0203696,0.00322482,0.0229054,0.0400957};
			vec.assign(err, err + sizeof(err) / sizeof(err[0]));
		}
		if (model == "T2tt")    {
			float err[] = {-0.0183669,-0.0125071,-0.0174156,-0.0164738,-0.0167200,-0.0149260,-0.0180894,-0.0154648,-0.0141536,-0.0119079,-0.0206974,0.000753522,0.0221000,0.0209901};
			vec.assign(err, err + sizeof(err) / sizeof(err[0]));
		}
		if (model == "T3bw")    {
			float err[] = {-0.0126896,-0.00609615,-0.00314039,-0.00273418,-0.00209022,0.000352204,0.000533044,0.00463945,-0.000409096,-0.00550145,-0.0442329,-0.0519994,-0.0384817,-0.0126860};
			vec.assign(err, err + sizeof(err) / sizeof(err[0]));
		}
	} // end c-jets

	else { // FOR LIGHT JETS
		if (fabs(jeteta) < 1.2){ // CENTRAL
			float SF[]       = {1.21878,1.28615,1.37535,1.38966,1.40320,1.49835,1.44308,1.58198,1.55687,1.65790,1.90233,1.92259,2.66174,3.08688};
			float SF_error[] = {0.0182686,0.0373732,0.0461870,0.0288973,0.0333528,0.0513836,0.0420353,0.106627,0.0658359,0.117285,0.185533,0.214071,0.487274,0.871502}; // stat + PU
			vecSF .assign(SF, SF + sizeof(SF) / sizeof(SF[0]));
			vecErr.assign(SF_error, SF_error + sizeof(SF_error) / sizeof(SF_error[0]));
			if (model == "T1")      {
				float err[]   = {0.386488,0.524838,0.679631,0.682134,0.731629,0.757618,0.695844,0.724127,0.623855,0.660598,0.829350,0.905624,1.40528,1.87998};
				vec.assign(err, err + sizeof(err) / sizeof(err[0]));
			}
			if (model == "T1bbbb")  {
				float err[]   = {0.269423,0.377897,0.460441,0.456512,0.470195,0.481113,0.438144,0.464775,0.347294,0.411220,0.550301,0.623299,1.14485,1.53694};
				vec.assign(err, err + sizeof(err) / sizeof(err[0]));
			}
			if (model == "T1tttt")  {
				float err[]   = {0.118477,0.162964,0.223318,0.220063,0.222306,0.267305,0.222287,0.283804,0.252221,0.324747,0.527015,0.659528,1.19317,1.50547};
				vec.assign(err, err + sizeof(err) / sizeof(err[0]));
			}
			if (model == "T2")      {
				float err[]   = {0.411044,0.558016,0.670788,0.776938,0.802923,0.895418,0.806768,0.812508,0.671626,0.715303,0.865253,0.889535,1.36529,1.83958};
				vec.assign(err, err + sizeof(err) / sizeof(err[0]));
			}
			if (model == "T2bb")    {
				float err[]   = {0.419325,0.559732,0.664588,0.701051,0.760269,0.876007,0.749539,0.827054,0.657627,0.702294,0.858618,0.837998,1.36137,1.75727};
				vec.assign(err, err + sizeof(err) / sizeof(err[0]));
			}
			if (model == "T2bw")    {
				float err[]   = {0.249750,0.331114,0.370544,0.380683,0.375024,0.433907,0.370687,0.430421,0.397210,0.479439,0.675053,0.815746,1.36142,1.86164};
				vec.assign(err, err + sizeof(err) / sizeof(err[0]));
			}
			if (model == "T2tt")    {
				float err[]   = {0.241447,0.297617,0.365921,0.372697,0.378869,0.434225,0.385061,0.452832,0.412124,0.498940,0.675028,0.813003,1.31961,1.57929};
				vec.assign(err, err + sizeof(err) / sizeof(err[0]));
			}
			if (model == "T3w")     {
				float err[]   = {0.287246,0.388381,0.480550,0.504640,0.531340,0.572774,0.532622,0.586227,0.529575,0.600174,0.792410,0.882505,1.42788,1.91256};
				vec.assign(err, err + sizeof(err) / sizeof(err[0]));
			}
		}
		else if (fabs(jeteta) >= 1.2 && fabs(jeteta) <2.4){ // FORWARD
			float SF[]       = {1.46970,1.48732,1.69024,1.64494,1.79297,1.90760,1.99867,2.21659,2.20103,2.42645,2.67594,4.24735,3.98979,15.0457};
			float SF_error[] = {0.104716,0.0392025,0.106315,0.115751,0.106807,0.0642086,0.138742,0.182345,0.169922,0.297889,0.320088,0.927736,1.24666,15.1860}; // stat + PU
			vecSF .assign(SF, SF + sizeof(SF) / sizeof(SF[0]));
			vecErr.assign(SF_error, SF_error + sizeof(SF_error) / sizeof(SF_error[0]));
			if (model == "T1")      {
				float err[]   = {0.686015,0.720420,0.991786,1.02970,1.22030,1.34286,1.64405,2.09951,2.92008,4.41435,6.25081,9.86965,11.9982,25.3907};
				vec.assign(err, err + sizeof(err) / sizeof(err[0]));
			}
			if (model == "T1bbbb")  {
				float err[]   = {0.618506,0.635309,0.794568,0.803646,0.886742,0.988190,1.10888,1.31924,1.49423,1.96107,2.36883,3.67770,4.80525,11.2853};
				vec.assign(err, err + sizeof(err) / sizeof(err[0]));
			}
			if (model == "T1tttt")  {
				float err[]   = {0.772327,0.874528,1.19814,1.24806,1.49608,1.73841,2.00430,2.54257,3.27898,4.35726,5.31846,7.44186,9.19039,15.6896};
				vec.assign(err, err + sizeof(err) / sizeof(err[0]));
			}
			if (model == "T2")      {
				float err[]   = {0.549720,0.580865,0.765356,0.788296,0.913463,1.03193,1.19510,1.40819,1.77895,2.69320,3.44912,5.91765,9.20944,18.0392};
				vec.assign(err, err + sizeof(err) / sizeof(err[0]));
			}
			if (model == "T2bb")    {
				float err[]   = {0.569517,0.548840,0.763820,0.726966,0.927079,0.959964,1.09951,1.29263,1.39602,1.97896,2.41141,3.67147,4.17557,11.7192};
				vec.assign(err, err + sizeof(err) / sizeof(err[0]));
			}
			if (model == "T2bw")    {
				float err[]   = {0.697431,0.759470,1.03429,1.05697,1.23279,1.38067,1.52550,1.89634,2.29738,2.87713,3.64427,5.54452,6.93274,13.9094};
				vec.assign(err, err + sizeof(err) / sizeof(err[0]));
			}
			if (model == "T2tt")    {
				float err[]   = {0.694594,0.753930,1.01105,1.02488,1.18455,1.36025,1.57676,1.87545,2.24691,2.81635,3.46050,5.75946,6.89900,15.8855};
				vec.assign(err, err + sizeof(err) / sizeof(err[0]));
			}
			if (model == "T3w")     {
				float err[]   = {0.773925,0.839802,1.08844,1.16056,1.35051,1.53349,1.76781,2.23936,2.95149,4.32639,5.98244,9.13821,12.3203,24.2016};
				vec.assign(err, err + sizeof(err) / sizeof(err[0]));
			}
		}
	} // end light jets

	// ---------------------------------------------------------------------------------
	// the sms specific scale factor is to be added in quadrature to the pu + stat scale 
	// factor which is due to fastsim vs. fullsim and is denoted SF_error!
	// ---------------------------------------------------------------------------------
	float ptmin[] = {30, 40, 50, 60, 70, 80, 100, 120, 160, 210, 260, 320, 400, 500};
	float ptmax[] = {40, 50, 60, 70, 80,100, 120, 160, 210, 260, 320, 400, 500, 670};
	int binnumber=0;
	int nbins = sizeof(ptmin)/sizeof(int)-1;
	for (int i=0; i<nbins; ++i) {
		if(jetpt>ptmin[i] && jetpt<ptmax[i]) {
			binnumber=i;
			break;
		}
	}
	if (jetpt > ptmax[nbins]) binnumber = nbins; //take the last bin for very energetic jets... no increase of uncertainty for now.

	SFandErr valerr;
	// cout << Form("jetflavor: %d jetpt: %f jeteta: %f binnumber: %d size of vecSF: %d", jetflavor, jetpt, jeteta, binnumber, vecSF.size()) << endl;
	valerr.sf  = vecSF[binnumber];
	valerr.err = sqrt(vecErr[binnumber]*vecErr[binnumber] + vec[binnumber]*vec[binnumber]);

	return valerr;
}
