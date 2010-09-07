#include <cmath>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <cstdlib>

#include "helper/FPRatios.hh"

using namespace std;

FPRatios::FPRatios() : verbose(1), initTopol(0), initVar(0){}

FPRatios::~FPRatios(){}

// *************************************************************************
int FPRatios::NevtTopol(int nemaxi, int nmmaxi, vector<double> npassi){
// Defines the topologies
// nemaxi = number of electrons
// nmmaxi = number of muons
// npassi = vector with number of events in each topology (Nti in the note)

	ne = nemaxi;
	nm = nmmaxi;
	ntopol = npassi.size();
	if (ntopol != (ne+1)*(nm+1) ){ 
		if(verbose > 5) cout << " ntopol = " << ntopol <<  " ne = " << ne << " nm = " << nm << " is wrong!" << endl;
		return 0; 
	}

// define the number of events passing and failing
// i.e. the lepton topologies, N_tx in the note
// remember: #topologies = #of fake/prompt combinations (called FPcnfigs)
	initTopol = 1;
	for (int i = 0; i < ntopol; ++i) {
		npass.push_back(npassi[i]);
	}

	if(verbose > 5) cout << endl;
	if(verbose > 5) cout << " Calculation of prompts and fakes" << endl;
	if(verbose > 5) cout << " ================================" << endl;
	if(verbose > 5) cout << endl;
	for (int i = 0; i < ntopol; ++i) {
		if(verbose > 5) cout << " Nt" << i << " = " << npass[i];
		if (i < ntopol-1) if(verbose > 5) cout << ", ";
	}
	if(verbose > 5) cout << endl;

// define the combinations of fakes and prompts for each topology
	int ife = -1;
	for (int i = 0; i < ntopol; ++i) {
		int ifm = i % (nm+1);
		if (ifm  == 0 ) {ife++;}
	//    if(verbose > 5) cout << " ife = " << ife << ", ifm = " << ifm << endl;
		nfetop.push_back(ife);
		npetop.push_back(ne-ife);
		nfmtop.push_back(ifm);
		npmtop.push_back(nm-ifm);
	}

	return 1;

}

// *************************************************************************
void FPRatios::SetFratios(double pRatioei, double pRatiomi, double pRatioErrei, double pRatioErrmi){
	fFRatioe1 =  pRatioei;
	fFRatiom1 =  pRatiomi;
	fFRatioErre1 = pRatioErrei;
	fFRatioErrm1 = pRatioErrmi;
	if(verbose > 5) cout << " fFRatioe1 = " << fFRatioe1 << " fFRatiom1 = " << fFRatiom1 
	     << " fFRatioErre1 =  " << fFRatioErre1 << " fFRatioErrm1 = " << fFRatioErrm1 << endl;
}

// *************************************************************************
void FPRatios::SetPratios(double pRatioei, double pRatiomi, double pRatioErrei, double pRatioErrmi){
	fPRatioe1 =  pRatioei;
	fPRatiom1 =  pRatiomi;
	fPRatioErre1 = pRatioErrei;
	fPRatioErrm1 = pRatioErrmi;
	if(verbose > 5) cout << " fPRatioe1 = " << fPRatioe1 << " fPRatiom1 = " << fPRatiom1 
	     << " fPRatioErre1 =  " << fPRatioErre1 << " fPRatioErrm1 = " << fPRatioErrm1 << endl;
}

// *************************************************************************
void FPRatios::SetFratios2(double pRatioei, double pRatiomi, double pRatioErrei, double pRatioErrmi){
	fFRatioe2 =  pRatioei;
	fFRatiom2 =  pRatiomi;
	fFRatioErre2 = pRatioErrei;
	fFRatioErrm2 = pRatioErrmi;
	if(verbose > 5) cout << " fFRatioe2 = " << fFRatioe2 << " fFRatiom2 = " << fFRatiom2 
	     << " fFRatioErre2 =  " << fFRatioErre2 << " fFRatioErrm2 = " << fFRatioErrm2 << endl;
}

// *************************************************************************
void FPRatios::SetPratios2(double pRatioei, double pRatiomi, double pRatioErrei, double pRatioErrmi){
	fPRatioe2 =  pRatioei;
	fPRatiom2 =  pRatiomi;
	fPRatioErre2 = pRatioErrei;
	fPRatioErrm2 = pRatioErrmi;
	if(verbose > 5) cout << " fPRatioe2 = " << fPRatioe2 << " fPRatiom2 = " << fPRatiom2 
	     << " fPRatioErre2 =  " << fPRatioErre2 << " fPRatioErrm2 = " << fPRatioErrm2 << endl;
}

// *************************************************************************
vector<double> FPRatios::NevtPass(){
	vector<double> vpt, veta;
	vpt.push_back(0.);
	veta.push_back(0.);
	return NevtPass(vpt, veta);
}

// *************************************************************************
vector<double> FPRatios::NevtPass(vector<double> vpt, vector<double> veta){
// returns a vector with the number of events for all FPconfigs
	vector<double> result;
// loop over all FPconfigs
	int isav = 0;
	for (int it = 0; it < ntopol; ++it) {

	// define the number of fakes and prompts for this FPconfigs
		int nfe = nfetop[it];
		int npe = npetop[it];
		int nfm = nfmtop[it];
		int npm = npmtop[it];
	// start with the topology with all leptons failing
		vector<int> ipass;
		for (int i = 0; i < ne; ++i) ipass.push_back(0);
		for (int i = 0; i < nm; ++i) ipass.push_back(0);
		if (verbose > 5) {
			if(verbose > 5) cout << endl;
			if(verbose > 5) cout << " nfe = " << nfe << ", npe = " << npe
				<< ", nfm = " << nfm << ", npm = " << npm << endl;
		}

	// loop over all topologies of passing and failing leptons
		int ie = 0;
		int im = 0;
		double wgttot = 0., dntotdpe = 0., dntotdfe = 0., dntotdpm = 0., dntotdfm= 0.;
		double ntot = 0., dntotstat = 0., dntotsyst = 0.;
		int itop = 0;
	// loop: increase the number of muons passing
		while (1) {
		// loop: increase the number of electrons passing
			while (1) {
				if (verbose > 5) {
					if(verbose > 5) cout << " pass = ";
					for (int ipr = 0; ipr < ne+nm; ++ipr) if(verbose > 5) cout << " " << ipass[ipr];
					if(verbose > 5) cout << endl;
				}

	// compute the weight for this combination of prompt and fakes
	//   and each topology in turn
	// (FPWeight also computes the derivatives for the error calculation)
				double wgt = FPWeight(npe, nfe, npm, nfm, vpt, veta, ipass);
	//	int isgn = abs(nfe + nfm - ie - im);

				int isgn = abs(npe + npm + ie + im);
	//      if(verbose > 5) cout << " isgn = " << isgn << endl;
				double sgn = 1.;
				if (isgn%2 != 0) sgn = -1.;
				wgt *= sgn;

	// save the signed weights into the vector wgtitot (w_i in the note)
				wgtitot.push_back(wgt);
	// also save the weight contributions from e and mu
				wgtitote.push_back(wgtie);
				wgtitotm.push_back(wgtim);

	// compute the number of events for this topology and its errors
				double ntop = wgt * npass[itop];
				double dntopstat = wgt * wgt * npass[itop];
				double dntopdpe = sgn * FPWeightError(1) * npass[itop];
				double dntopdfe = sgn * FPWeightError(2) * npass[itop];
				double dntopdpm = sgn * FPWeightError(3) * npass[itop];
				double dntopdfm = sgn * FPWeightError(4) * npass[itop];

	// add the weight of this topology to the total one and also the errors (for e and mu)
				wgttot += wgt;
				dntotdpe += dntopdpe;
				dntotdfe += dntopdfe;
				dntotdpm += dntopdpm;
				dntotdfm += dntopdfm;
				if (verbose > 5) {
					if(verbose > 5) cout << " weight = " << wgt << ", dndpe = " << dntopdpe << ", dndfe = " << dntopdfe
						<< ", dndpm = " << dntopdpm<< ", dndfm = " << dntopdfm << endl;
					if(verbose > 5) cout << " ntop = " << ntop << " +- " << sqrt(dntopstat) << "_stat" << endl;
				}

	// save errors on the weights for the calculation of correlations
				vector<double> temp = GetWgtError(1);
	//      for (int ite = 0; ite < temp.size(); ++ite) if(verbose > 5) cout << " temp(1) = " << temp[ite] << endl;
				dwgtdptote.push_back(temp);
				temp.clear();
				temp = GetWgtError(2);
	//      for (int ite = 0; ite < temp.size(); ++ite) if(verbose > 5) cout << " temp(2) = " << temp[ite] << endl;
				dwgtdftote.push_back(temp);
				temp.clear();
				temp = GetWgtError(3);
	//      for (int ite = 0; ite < temp.size(); ++ite) if(verbose > 5) cout << " temp(1) = " << temp[ite] << endl;
				dwgtdptotm.push_back(temp);
				temp.clear();
				temp = GetWgtError(4);
	//      for (int ite = 0; ite < temp.size(); ++ite) if(verbose > 5) cout << " temp(2) = " << temp[ite] << endl;
				dwgtdftotm.push_back(temp);
				if (verbose > 5) {
					for (int iv = 0; iv < temp.size(); ++iv) {
						if(verbose > 5) cout << " dwgtdptote = " << dwgtdptote[isav][iv] 
							<< ", dwgtdftote = " << dwgtdftote[isav][iv]
							<< ", dwgtdptotm = " << dwgtdptotm[isav][iv] 
							<< ", dwgtdftotm = " << dwgtdftotm[isav][iv] << endl;
					}
				}
				temp.clear();

				isav++;
				itop++;
				ntot += ntop;
				dntotstat += dntopstat;

	// declare next muons as passing
				ipass[ne+im] = 1;
				im++;
				if (im > nm) break;
			}
		// declare next electron as passing
			ipass[ie] = 1;
			for (int j = 0; j < nm; ++j) {
				ipass[ne+j] = 0;
			}
			im = 0;
			ie++;
			if (ie > ne) break;
		}
		if (verbose > 6) {
			if(verbose > 5) cout << " Total weight = " << wgttot 
				<< ", dntotdpe = " << dntotdpe << ", dntotdfe = " << dntotdfe
				<< ", dntotdpm = " << dntotdpm << ", dntotdfm = " << dntotdfm << endl;
		}
	// compute for this FPconfig the final statistical and systematic errors
		dntotstat = sqrt(dntotstat);
		dntotsyst = sqrt(dntotdpe*dntotdpe + dntotdfe*dntotdfe + dntotdpm*dntotdpm + dntotdfm*dntotdfm);
		if (verbose > 5) {
			if(verbose > 5) cout << " N_nfe,npe,nfm,npm = N_" << nfe << "," << npe << "," << nfm << "," << npm << " = " 
				<< ntot << " +- " << dntotstat << "_stat +- " 
				<< dntotsyst << "_syst" << endl;
		}
		result.push_back(ntot);
		ntotfp.push_back(ntot);
		dNtotstat.push_back(dntotstat);
		dNtotsyst.push_back(dntotsyst);
	}
// end loop over FPconfigs

	return result;
	// return ntotfp;
}

// *************************************************************************
vector<double> FPRatios::NevtPass(vector<double> f, vector<double> ferr, vector<double> p, vector<double> perr){
// returns a vector with the number of events for all FPconfigs
	vector<double> result;
// loop over all FPconfigs
	int isav = 0;
	for (int it = 0; it < ntopol; ++it) {

	// define the number of fakes and prompts for this FPconfigs
		int nfe = nfetop[it];
		int npe = npetop[it];
		int nfm = nfmtop[it];
		int npm = npmtop[it];
	// start with the topology with all leptons failing
		vector<int> ipass;
		for (int i = 0; i < ne; ++i) ipass.push_back(0);
		for (int i = 0; i < nm; ++i) ipass.push_back(0);
		if (verbose > 5) {
			if(verbose > 5) cout << endl;
			if(verbose > 5) cout << " nfe = " << nfe << ", npe = " << npe
				<< ", nfm = " << nfm << ", npm = " << npm << endl;
		}

	// loop over all topologies of passing and failing leptons
		int ie = 0;
		int im = 0;
		double wgttot = 0., dntotdpe = 0., dntotdfe = 0., dntotdpm = 0., dntotdfm= 0.;
		double ntot = 0., dntotstat = 0., dntotsyst = 0.;
		int itop = 0;
	// loop: increase the number of muons passing
		while (1) {
		// loop: increase the number of electrons passing
			while (1) {
				if (verbose > 5) {
					if(verbose > 5) cout << " pass = ";
					for (int ipr = 0; ipr < ne+nm; ++ipr) if(verbose > 5) cout << " " << ipass[ipr];
					if(verbose > 5) cout << endl;
				}

	// compute the weight for this combination of prompt and fakes
	//   and each topology in turn
	// (FPWeight also computes the derivatives for the error calculation)
				double wgt = FPWeight(npe, nfe, npm, nfm, f, ferr, p, perr, ipass);
	//	int isgn = abs(nfe + nfm - ie - im);

				int isgn = abs(npe + npm + ie + im);
	//      if(verbose > 5) cout << " isgn = " << isgn << endl;
				double sgn = 1.;
				if (isgn%2 != 0) sgn = -1.;
				wgt *= sgn;

	// save the signed weights into the vector wgtitot (w_i in the note)
				wgtitot.push_back(wgt);
	// also save the weight contributions from e and mu
				wgtitote.push_back(wgtie);
				wgtitotm.push_back(wgtim);

	// compute the number of events for this topology and its errors
				double ntop = wgt * npass[itop];
				double dntopstat = wgt * wgt * npass[itop];
				double dntopdpe = sgn * FPWeightError(1) * npass[itop];
				double dntopdfe = sgn * FPWeightError(2) * npass[itop];
				double dntopdpm = sgn * FPWeightError(3) * npass[itop];
				double dntopdfm = sgn * FPWeightError(4) * npass[itop];

	// add the weight of this topology to the total one and also the errors (for e and mu)
				wgttot += wgt;
				dntotdpe += dntopdpe;
				dntotdfe += dntopdfe;
				dntotdpm += dntopdpm;
				dntotdfm += dntopdfm;
				if (verbose > 5) {
					if(verbose > 5) cout << " weight = " << wgt << ", dndpe = " << dntopdpe << ", dndfe = " << dntopdfe
						<< ", dndpm = " << dntopdpm<< ", dndfm = " << dntopdfm << endl;
					if(verbose > 5) cout << " ntop = " << ntop << " +- " << sqrt(dntopstat) << "_stat" << endl;
				}

	// save errors on the weights for the calculation of correlations
				vector<double> temp = GetWgtError(1);
	//      for (int ite = 0; ite < temp.size(); ++ite) if(verbose > 5) cout << " temp(1) = " << temp[ite] << endl;
				dwgtdptote.push_back(temp);
				temp.clear();
				temp = GetWgtError(2);
	//      for (int ite = 0; ite < temp.size(); ++ite) if(verbose > 5) cout << " temp(2) = " << temp[ite] << endl;
				dwgtdftote.push_back(temp);
				temp.clear();
				temp = GetWgtError(3);
	//      for (int ite = 0; ite < temp.size(); ++ite) if(verbose > 5) cout << " temp(1) = " << temp[ite] << endl;
				dwgtdptotm.push_back(temp);
				temp.clear();
				temp = GetWgtError(4);
	//      for (int ite = 0; ite < temp.size(); ++ite) if(verbose > 5) cout << " temp(2) = " << temp[ite] << endl;
				dwgtdftotm.push_back(temp);
				if (verbose > 5) {
					for (int iv = 0; iv < temp.size(); ++iv) {
						if(verbose > 5) cout << " dwgtdptote = " << dwgtdptote[isav][iv] 
							<< ", dwgtdftote = " << dwgtdftote[isav][iv]
							<< ", dwgtdptotm = " << dwgtdptotm[isav][iv] 
							<< ", dwgtdftotm = " << dwgtdftotm[isav][iv] << endl;
					}
				}
				temp.clear();

				isav++;
				itop++;
				ntot += ntop;
				dntotstat += dntopstat;

	// declare next muons as passing
				ipass[ne+im] = 1;
				im++;
				if (im > nm) break;
			}
		// declare next electron as passing
			ipass[ie] = 1;
			for (int j = 0; j < nm; ++j) {
				ipass[ne+j] = 0;
			}
			im = 0;
			ie++;
			if (ie > ne) break;
		}
		if (verbose > 6) {
			if(verbose > 5) cout << " Total weight = " << wgttot 
				<< ", dntotdpe = " << dntotdpe << ", dntotdfe = " << dntotdfe
				<< ", dntotdpm = " << dntotdpm << ", dntotdfm = " << dntotdfm << endl;
		}
	// compute for this FPconfig the final statistical and systematic errors
		dntotstat = sqrt(dntotstat);
		dntotsyst = sqrt(dntotdpe*dntotdpe + dntotdfe*dntotdfe + dntotdpm*dntotdpm + dntotdfm*dntotdfm);
		if (verbose > 5) {
			if(verbose > 5) cout << " N_nfe,npe,nfm,npm = N_" << nfe << "," << npe << "," << nfm << "," << npm << " = " 
				<< ntot << " +- " << dntotstat << "_stat +- " 
				<< dntotsyst << "_syst" << endl;
		}
		result.push_back(ntot);
		ntotfp.push_back(ntot);
		dNtotstat.push_back(dntotstat);
		dNtotsyst.push_back(dntotsyst);
	}
// end loop over FPconfigs

	return result;
	// return ntotfp;
}

// *************************************************************************
double FPRatios::FPWeightError (int ityp){
// returns the error components on the various weights
// for ityp = 1: error from dp for electron
// for ityp = 2: error from df for electron
// for ityp = 3: error from dp for muon
// for ityp = 4: error from df for muon

	double dwgt = 0.;
// for e and mu systematic errors
	double we = wgtie;
	double wm = wgtim;
	if (we == 0.) we = 1.;
	if (wm == 0.) wm = 1.;
	if (ityp == 1)      dwgt = dwgtdpe*wm;
	else if (ityp == 2) dwgt = dwgtdfe*wm;
	else if (ityp == 3) dwgt = dwgtdpm*we;
	else if (ityp == 4) dwgt = dwgtdfm*we;

	return dwgt;
}

// *************************************************************************
vector<double> FPRatios::GetWgtError(int ityp){
// for ityp = 1: error from dp for electron
// for ityp = 2: error from df for electron
// for ityp = 3: error from dp for muon
// for ityp = 4: error from df for muon

	vector<double> dwgt;
	if (ityp == 1) {
		for (int i = 0; i < dwgtdptope.size(); ++i) {
			dwgt.push_back(dwgtdptope[i]);
		}
	} else if (ityp == 2) {
		for (int i = 0; i < dwgtdftope.size(); ++i) {
			dwgt.push_back(dwgtdftope[i]);
		}
	} else if (ityp == 3) {
		for (int i = 0; i < dwgtdptopm.size(); ++i) {
			dwgt.push_back(dwgtdptopm[i]);
		}
	} else {
		for (int i = 0; i < dwgtdftopm.size(); ++i) {
			dwgt.push_back(dwgtdftopm[i]);
		}
	}

	return dwgt;
}

// *************************************************************************
vector<double> FPRatios::NevtPassErrStat(){
// return the vector of statistical errors for all FPconfigs

	return dNtotstat;
}

// *************************************************************************
vector<double> FPRatios::NevtPassErrSyst(){
// return the vector of systematic errors for all FPconfigs

	return dNtotsyst;
}

// *************************************************************************
void FPRatios::Varnf1nf2(){
// Computes the variances for all the cross combinations of N_nf

//  if(verbose > 5) cout << endl;
//  if(verbose > 5) cout << " Computing correlations " << endl;
	if (verbose > 5) {
		for (int iv = 0; iv < dwgtdptote.size(); ++iv) {
			for (int i = 0; i < dwgtdptote[iv].size(); ++i) {
				if(verbose > 5) cout << " " << iv << " " << i << " dwgtdptote = " << dwgtdptote[iv][i]
					<< ", dwgtdftote = " << dwgtdftote[iv][i]
					<< ", dwgtdptotm = " << dwgtdptotm[iv][i]
					<< ", dwgtdftotm = " << dwgtdftotm[iv][i] << endl;
			}
		}
	}
// loop over the first FPconfig
	for (int nf1 = 0; nf1 < ntopol; ++nf1) {
	// loop over the second FPconfig
		for (int nf2 = nf1 + 1; nf2 < ntopol; ++nf2) {
		//    for (int nf2 = 0; nf2 < ntopol; ++nf2) {
			if (verbose > 5) {
				if(verbose > 5) cout << endl;
				if(verbose > 5) cout << " nf1 = " << nf1 << ", nf2 = " << nf2 << endl;
			}
		// get the variances on the weights into vector varwi for these two FPconfigs
			vector<double> varwi;
			varwi = GetVarWeights(ntopol, nf1, nf2);
			double vstat = 0.;
			double vsyst = 0;
			int iwi = 0;
		// loop over the topologies for the first FPconfig
			for (int i = 0; i < ntopol; ++i) {
				int ivar1 = nf1 * ntopol + i;
				int ivar2 = nf2 * ntopol + i;
	// the statistical variance
				vstat += wgtitot[ivar1] * wgtitot[ivar2] * npass[i];

				int nfefake1 = nfetop[nf1];
				int nfmfake1 = nfmtop[nf1];
				int nepass1 = i / (nm+1);
				int nmpass1 = i % (nm+1);
				int isgn;
				isgn = nfefake1 - ne + nepass1 + nfmfake1 - nm + nmpass1;
				double sgn1 = 1.;
				if (isgn%2 != 0) sgn1 = -1.;
	//	if(verbose > 5) cout << "  nfake1 = " << nfefake1 << " " << nfmfake1 << ", npass = " << nepass1 << " " << nmpass1 
	//	     << "sgn = " << sgn1 << endl;

	// loop over the topologies for the second FPconfig
				for (int k = 0; k < ntopol; ++k) {
					int nfefake2 = nfetop[nf2];
					int nfmfake2 = nfmtop[nf2];
					int nepass2 = k / (nm+1);
					int nmpass2 = k % (nm+1);
					isgn = nfefake2 - ne + nepass2 + nfmfake2 - nm + nmpass2;
					double sgn2 = 1.;
					if (isgn%2 != 0) sgn2 = -1.;
	//	  if(verbose > 5) cout << "  nfake2 = " << nfefake2 << " " << nfmfake2 << ", npass = " << nepass2 << " " << nmpass2 
	//	       << "sgn = " << sgn2 << endl;
	// multiply them with the nbers of events in the topologies
					double varnew = sgn1 * npass[i] * sgn2 * npass[k] * varwi[iwi];
					if (verbose > 5) {
						if(verbose > 5) cout << " sgn1 = " << sgn1 << " npass = " << npass[i]
							<< ", sgn2 = " << sgn2 << " npass = " << npass[k] << ", varwi = " << varwi[iwi] 
							<< ", varnew = " << varnew << endl;
					}
					vsyst += varnew;
					iwi++;
				}
			}
			if (verbose > 5) {
				if(verbose > 5) cout << " Variance stat = " << vstat << ", syst = " << vsyst << endl;
			}

		// save the stat and syst errors in global vectors
			varstat.push_back(vstat);
			varsyst.push_back(vsyst);
		}
	}


	return;
}

// *************************************************************************
vector<double> FPRatios::GetVarWeights(int ntopol, int itop1, int itop2){
// computes the variances V(wi, wk) for a chosen FPconfig

	int i1 = itop1 * ntopol;
	int i2 = itop2 * ntopol;
	vector<double> varwi;
	if (verbose > 6) {
		if(verbose > 5) cout << " topol = " << itop1 << " " << itop2 << ", index = " << i1 << " " << i2 << endl;
	}

// store the V(wi,wk) for all topologies in a vector
	for (int j1 = i1; j1 < i1+ntopol; ++j1) {
		double wgte1 = wgtitote[j1];
		double wgtm1 = wgtitotm[j1];
		if (wgte1 == 0.) wgte1 = 1.;
		if (wgtm1 == 0.) wgtm1 = 1.;
		for (int j2 = i2; j2 < i2+ntopol; ++j2) {
			double wgte2 = wgtitote[j2];
			double wgtm2 = wgtitotm[j2];
			if (wgte2 == 0.) wgte2 = 1.;
			if (wgtm2 == 0.) wgtm2 = 1.;
			double sumwij = 0.;
			double sumdwijdpe1 = 0.;
			double sumdwijdfe1 = 0.;
			double sumdwijdpm1 = 0.;
			double sumdwijdfm1 = 0.;
			double sumdwijdpe2 = 0.;
			double sumdwijdfe2 = 0.;
			double sumdwijdpm2 = 0.;
			double sumdwijdfm2 = 0.;
			for (int k1 = 0; k1 < dwgtdptote[j1].size(); ++k1) {
				sumdwijdpe1 += dwgtdptote[j1][k1];
				sumdwijdfe1 += dwgtdftote[j1][k1];
				sumdwijdpm1 += dwgtdptotm[j1][k1];
				sumdwijdfm1 += dwgtdftotm[j1][k1];
			}
			for (int k2 = 0; k2 < dwgtdptote[j2].size(); ++k2) {
				sumdwijdpe2 += dwgtdptote[j2][k2];
				sumdwijdfe2 += dwgtdftote[j2][k2];
				sumdwijdpm2 += dwgtdptotm[j2][k2];
				sumdwijdfm2 += dwgtdftotm[j2][k2];
			}
			sumwij = sumdwijdpe1 * wgtm1 * sumdwijdpe2 * wgtm2
				+ sumdwijdfe1 * wgtm1 * sumdwijdfe2 * wgtm2
				+ sumdwijdpm1 * wgte1 * sumdwijdpm2 * wgte2
				+ sumdwijdfm1 * wgte1 * sumdwijdfm2 * wgte2;
	//	  if(verbose > 5) cout << " ( " << j1 << ", " << k1 << ") * ( " << j2 << ", " << k2 << ")" << endl;
			varwi.push_back(sumwij);
			if (verbose > 6) {
				if(verbose > 5) cout << " " << j1 << " " << j2 << " varwi = " << sumwij << endl;
			}
		} 
	}

	return varwi;
}

// *************************************************************************
vector<double> FPRatios::VarStat(){
// returns the statistical variances for all FPratio combinations into a vector
	if (!initVar) {Varnf1nf2();}
	initVar = 1;
	return varstat;
}

// *************************************************************************
vector<double> FPRatios::VarSyst(){
// returns the systematic variances for all FPratio combinations into a vector
	if (!initVar) {Varnf1nf2();}
	initVar = 1;
	return varsyst;
}

// *************************************************************************
double FPRatios::FPWeight (int npe, int nfe, int npm, int nfm, vector<double> vpt, vector<double> veta, vector<int> ipass){
// computes the weight per event for a given combination of prompt and fakes
// and fixed numbers of passing and failing leptons (e or mu)
// Note: it is the weight for a given configuration N(npe, nfe, npm, nfm) 
//       for events passing the tight (analysis) criteria
//       i.e. the wi in the note
// npe = number of prompt electrons
// nfe = number of fake electrons
// npm = number of prompt muons
// nfm = number of fake muons
// the vectors contain the lepton variables in the order:
//    npe prompt electrons, nfe fake electrons,
//    npm prompts muons, nfm fake muons
// vpt   = lepton pT
// veta  = lepton eta
// ipass = 1 if the lepton passes tight cuts, 0 if it fails
//
// it uses 4 functions to be provided by the user (here dummies are included)
// FakeRatio (int emu, double pt, double eta) for the fake ratio
// PromptRatio (int emu, double pt, double eta) for the propmt ratio
// FakeRatioError (int emu, double pt, double eta) for the fake ratio error
// PromptRatioError (int emu, double pt, double eta) for the propmt ratio error
// with emu = 1 for electrons, 2 for muons
//
// At the same time, it computes the errors and correlations, stored in global variables

// if = 0 is in loose cuts, = 1 in tight cuts
	int inTight = 1;

	vector<int> ind;
	vector<double> fakrat, prorat;
	vector<double> fakerr, proerr;
	double wgt = 1.;
	dwgtdpe = 0.;
	dwgtdfe = 0.;
	dwgtdpm = 0.;
	dwgtdfm = 0.;
	dwgtdptope.clear(); 
	dwgtdftope.clear();
	dwgtdptopm.clear(); 
	dwgtdftopm.clear();

	int ne = npe + nfe;
	int ipeinit = npe - 1;
	int ifeinit = ne - 1;
	int nm = npm + nfm;
	int ipminit = ne + npm - 1;
	int ifminit = ne + nm - 1;

	for (int i = 0; i < ne+nm; ++i) {
		if (i < npe) {
			ind.push_back(1);
			fakrat.push_back(FakeRatio(1, vpt[i], veta[i]) );
			prorat.push_back(PromptRatio(1, vpt[i], veta[i]) );
			fakerr.push_back(FakeRatioError(1, vpt[i], veta[i]) );
			proerr.push_back(PromptRatioError(1, vpt[i], veta[i]) );
		} else if (i < ne) {
			ind.push_back(0); 
			fakrat.push_back(FakeRatio(1, vpt[i], veta[i]) );
			prorat.push_back(PromptRatio(1, vpt[i], veta[i]) );
			fakerr.push_back(FakeRatioError(1, vpt[i], veta[i]) );
			proerr.push_back(PromptRatioError(1, vpt[i], veta[i]) );
		} else if (i < ne+npm) {
			ind.push_back(1);
			fakrat.push_back(FakeRatio(2, vpt[i], veta[i]) );
			prorat.push_back(PromptRatio(2, vpt[i], veta[i]) );
			fakerr.push_back(FakeRatioError(2, vpt[i], veta[i]) );
			proerr.push_back(PromptRatioError(2, vpt[i], veta[i]) );
		} else {
			ind.push_back(0);
			fakrat.push_back(FakeRatio(2, vpt[i], veta[i]) );
			prorat.push_back(PromptRatio(2, vpt[i], veta[i]) );
			fakerr.push_back(FakeRatioError(2, vpt[i], veta[i]) );
			proerr.push_back(PromptRatioError(2, vpt[i], veta[i]) );
		}
	//     if(verbose > 5) cout << " i = " << i << ", ind = " << ind[i] 
	//	  << ", pT = " << vpt[i] << ", eta = " << veta[i] 
	//	  << ", fakrat = " << fakrat[i] << ", prorat = " << prorat[i] << endl;
	}

	int ipe = ipeinit;
	int ife = ifeinit;
	int ipemv = ipeinit;
	int ifemv = ifeinit;
	int ipm = ipminit;
	int ifm = ifminit;
	int ipmmv = ipminit;
	int ifmmv = ifminit;

	int icyce = 0;
	int icycm = 3;
	int nowmu = 0;
	if (ne == 0) nowmu = 1;
	wgtie = 0.;
	wgtim = 0.;
	wgt = 0.;

// sum over all the weights wij of the various FPconfigs
	while (1){

		// for (int i = 0; i < ne; ++i) {
		// 	if(verbose > 5) cout << " " << ind[i];
		// }
		// if(verbose > 5) cout << " | ";
		// for (int i = ne; i < ne+nm; ++i) {
		// 	if(verbose > 5) cout << " " << ind[i];
		// }
		// if(verbose > 5) cout << endl;
		// 

		// compute the weight for this combination, 
		//   multiplying all leptons of a given FPconfig
		double wgtnewe = 0., wgtnewm = 0.;
		double dwgtnewdpe = 0., dwgtnewdfe = 0.;
		double dwgtnewdpm = 0., dwgtnewdfm = 0.;
		double cpdpe = 0., cfdfe = 0.;
		double cpdpm = 0., cfdfm = 0.;
		for (int i = 0; i < ne+nm; ++i) {
			double p = prorat[i];
			double f = fakrat[i];
			double dp = proerr[i];
			double df = fakerr[i];

			if (i < ne) { // Electron
				if (wgtnewe != 0.) wgtnewe /= (p - f);
				else wgtnewe = 1. / (p - f);
				if (ind[i] && !ipass[i]) {
					wgtnewe *= f;
					if (inTight) wgtnewe *= p;
					if (p > 0.) cpdpe += (-1./(p-f) + 1./p);
					else cpdpe = 0.;
					if (f > 0.) cfdfe += (1./(p-f) + 1./f);
					else cfdfe = 0.;
				}
				if (ind[i] && ipass[i]) {
					wgtnewe *= (1. - f);
					if (inTight) wgtnewe *= p;
					if (p > 0.) cpdpe += (-1./(p-f) + 1./p);
					else cpdpe = 0.;
					cfdfe += (1./(p-f) - 1./(1.-f));
				}
				if (!ind[i] && !ipass[i]) {
					wgtnewe *= p;
					if (inTight) wgtnewe *= f;
					if (p > 0.) cpdpe += (-1./(p-f) + 1./p);
					else cpdpe = 0.;
					if (f > 0.) cfdfe += (1./(p-f) + 1./f);
					else cfdfe = 0.;
				}
				if (!ind[i] && ipass[i]) {
					wgtnewe *= (1. - p);
					if (inTight) wgtnewe *= f;
					if (p < 1.) cpdpe += (-1./(p-f) - 1./(1.-p));
					else cpdpe = 0.;
					if (f > 0.) cfdfe += (1./(p-f) + 1./f);
					else cfdfe = 0.;
				}
			} else { // Muons
				if (wgtnewm != 0.) wgtnewm /= (p - f);
				else wgtnewm = 1. / (p - f);
				if (ind[i] && !ipass[i]) {
					wgtnewm *= f;
					if (inTight) wgtnewm *= p;
					if (p > 0.) cpdpm += (-1./(p-f) + 1./p);
					else cpdpm = 0.;
					if (f > 0.) cfdfm += (1./(p-f) + 1./f);
					else cfdfm = 0.;
				}
				if (ind[i] && ipass[i]) {
					wgtnewm *= (1. - f);
					if (inTight) wgtnewm *= p;
					if (p > 0.) cpdpm += (-1./(p-f) + 1./p);
					else cpdpm = 0.;
					cfdfm += (1./(p-f) - 1./(1.-f));
				}
				if (!ind[i] && !ipass[i]) {
					wgtnewm *= p;
					if (inTight) wgtnewm *= f;
					if (p > 0.) cpdpm += (-1./(p-f) + 1./p);
					else cpdpm = 0.;
					if (f > 0.) cfdfm += (1./(p-f) + 1./f);
					else cfdfm = 0.;
				}
				if (!ind[i] && ipass[i]) {
					wgtnewm *= (1. - p);
					if (inTight) wgtnewm *= f;
					if (p < 1.) cpdpm += (-1./(p-f) - 1./(1.-p));
					else cpdpm = 0.;
					if (f > 0.) cfdfm += (1./(p-f) + 1./f);
					else cfdfm = 0.;
				}
			}
			dwgtnewdpe = wgtnewe * cpdpe * dp;
			dwgtnewdfe = wgtnewe * cfdfe * df;
			dwgtnewdpm = wgtnewm * cpdpm * dp;
			dwgtnewdfm = wgtnewm * cfdfm * df;
			if (verbose > 6) {
				if (i < ne) {
					if(verbose > 5) cout << " " << i << " wgtnewe = " << wgtnewe  << ", cpdpe = " << cpdpe << ", cfdfe = " << cfdfe
						<< ", dwgtnewdpe = " << dwgtnewdpe << ", dwgtnewdfe = " << dwgtnewdfe << endl;
				} else {
					if(verbose > 5) cout << " " << i << " wgtnewm = " << wgtnewm  << ", cpdpm = " << cpdpm << ", cfdfm = " << cfdfm
						<< ", dwgtnewdpm = " << dwgtnewdpm << ", dwgtnewdfm = " << dwgtnewdfm << endl;
				}
			}
		}
		wgtie += wgtnewe;
		wgtim += wgtnewm;
		double wgtnew;
		if (wgtnewe == 0.) wgtnew = wgtnewm;
		else if (wgtnewm == 0.) wgtnew = wgtnewe;
		else wgtnew = wgtnewe * wgtnewm;
		wgt += wgtnew;
		dwgtdpe += dwgtnewdpe;
		dwgtdfe += dwgtnewdfe;
		dwgtdpm += dwgtnewdpm;
		dwgtdfm += dwgtnewdfm;
		dwgtdptope.push_back(dwgtnewdpe); 
		dwgtdftope.push_back(dwgtnewdfe);
		dwgtdptopm.push_back(dwgtnewdpm); 
		dwgtdftopm.push_back(dwgtnewdfm);

	// generate next combination for electrons
		if (icyce == 1 && !nowmu) {
		//       if(verbose > 5) cout << " cycle over the fakes, ife = " << ife << endl;
			int temp = ind[ipeinit+1];
			for (int i = ipeinit+1; i < ifeinit+1; ++i) {
				ind[i] = ind[i+1];
			}
			ind[ifeinit] = temp;
			if (ife > ipeinit + 1 && ifemv > ipeinit) {
				ife--;
			} else {
				ife = ifeinit;
				icyce = 2;
			}
		}

		if (icyce == 2 && !nowmu) {
		//       if(verbose > 5) cout << " cycle over the prompts, ipe = " << ipe << endl;
			int temp = ind[0];
			for (int i = 0; i < ipeinit+1; ++i) {
				ind[i] = ind[i+1];
			}
			ind[ipeinit] = temp;
			if (ipe > 0 && ipemv >= 0) {
				ipe--;
				icyce = 1;
			} else {
				ipe = ipeinit;
				icyce = 0;
			}
		}

		if (icyce == 0 && !nowmu) {
			if (ipemv >= 0 && ifemv > ipeinit) {
	//          if(verbose > 5) cout << " move " << ipemv << " " << ifemv << endl;
				ind[ifemv] = 1;
				ind[ipemv] = 0;
				ipemv--;
				ifemv--;
				icyce = 1;
			} else {
				for (int i = 0; i < ne; ++i) {
					if (i < npe) ind[i] = 1;
					else ind[i] = 0; 
				}
				ipe = ipeinit;
				ife = ifeinit;
				ipemv = ipeinit;
				ifemv = ifeinit;
	//          if(verbose > 5) cout << " end cycling over electrons " << endl;
				nowmu = 1;
			}
		}

	// generate next combination for muons
		if (icycm == 4 && nowmu) {
		//       if(verbose > 5) cout << " cycle over the fakes, ifm = " << ifm << endl;
			int temp = ind[ipminit+1];
			for (int i = ipminit+1; i < ifminit+1; ++i) {
				ind[i] = ind[i+1];
			}
			ind[ifminit] = temp;
			if (ifm > ipminit + 1 && ifmmv > ipminit) {
				ifm--;
				if (ne > 0) nowmu = 0;
			} else {
				ifm = ifminit;
				icycm = 5;
			}
		}

		if (icycm == 5 && nowmu) {
		//       if(verbose > 5) cout << " cycle over the prompts, ipm = " << ipm << endl;
			int temp = ind[ne];
			for (int i = ne; i < ipminit+1; ++i) {
				ind[i] = ind[i+1];
			}
			ind[ipminit] = temp;
			if (ipm > ne && ipmmv >= ne) {
				ipm--;
				icycm = 4;
				if (ne > 0) nowmu = 0;
			} else {
				ipm = ipminit;
				icycm = 3;
			}
		}

		if (icycm == 3 && nowmu) {
			if (ipmmv >= ne && ifmmv > ipminit) {
	//          if(verbose > 5) cout << " move " << ipmmv << " " << ifmmv << endl;
				ind[ifmmv] = 1;
				ind[ipmmv] = 0;
				ipmmv--;
				ifmmv--;
				if (ne > 0) {
					icyce = 0;
					nowmu = 0;
				}
				icycm = 4;
			} else {
	//          if(verbose > 5) cout << " end cycling over muons " << endl;
				break;
			}
		}

	}

	if (verbose > 6) {
		if(verbose > 5) cout << " wgt = " << wgt << ", dwgtdpe = " << dwgtdpe << ", dwgtdfe = " << dwgtdfe
			<< ", dwgtdpm = " << dwgtdpm << ", dwgtdfm = " << dwgtdfm << endl;
	}

	ind.clear();
	fakrat.clear(); 
	prorat.clear();
	fakerr.clear(); 
	proerr.clear();

	return wgt;
}

// *************************************************************************
double FPRatios::FPWeight (int npe, int nfe, int npm, int nfm, vector<double> f, vector<double> ferr, vector<double> p, vector<double> perr, vector<int> ipass){
// computes the weight per event for a given combination of prompt and fakes
// and fixed numbers of passing and failing leptons (e or mu)
// Note: it is the weight for a given configuration N(npe, nfe, npm, nfm) 
//       for events passing the tight (analysis) criteria
//       i.e. the wi in the note
// npe = number of prompt electrons
// nfe = number of fake electrons
// npm = number of prompt muons
// nfm = number of fake muons
// the vectors contain the lepton variables in the order:
//    npe prompt electrons, nfe fake electrons,
//    npm prompts muons, nfm fake muons
// vpt   = lepton pT
// veta  = lepton eta
// ipass = 1 if the lepton passes tight cuts, 0 if it fails
//
// it uses 4 functions to be provided by the user (here dummies are included)
// FakeRatio (int emu, double pt, double eta) for the fake ratio
// PromptRatio (int emu, double pt, double eta) for the propmt ratio
// FakeRatioError (int emu, double pt, double eta) for the fake ratio error
// PromptRatioError (int emu, double pt, double eta) for the propmt ratio error
// with emu = 1 for electrons, 2 for muons
//
// At the same time, it computes the errors and correlations, stored in global variables

// if = 0 is in loose cuts, = 1 in tight cuts
	int inTight = 1;

	vector<int> ind;
	vector<double> fakrat, prorat;
	vector<double> fakerr, proerr;
	double wgt = 1.;
	dwgtdpe = 0.;
	dwgtdfe = 0.;
	dwgtdpm = 0.;
	dwgtdfm = 0.;
	dwgtdptope.clear(); 
	dwgtdftope.clear();
	dwgtdptopm.clear(); 
	dwgtdftopm.clear();

	int ne = npe + nfe;
	int ipeinit = npe - 1;
	int ifeinit = ne - 1;
	int nm = npm + nfm;
	int ipminit = ne + npm - 1;
	int ifminit = ne + nm - 1;

	for (int i = 0; i < ne+nm; ++i) {
		if (i < npe) {
			ind.push_back(1);
			fakrat.push_back(f[i]);
			prorat.push_back(p[i]);
			fakerr.push_back(ferr[i]);
			proerr.push_back(perr[i]);
		} else if (i < ne) {
			ind.push_back(0); 
			fakrat.push_back(f[i]);
			prorat.push_back(p[i]);
			fakerr.push_back(ferr[i]);
			proerr.push_back(perr[i]);
		} else if (i < ne+npm) {
			ind.push_back(1);
			fakrat.push_back(f[i]);
			prorat.push_back(p[i]);
			fakerr.push_back(ferr[i]);
			proerr.push_back(perr[i]);
		} else {
			ind.push_back(0);
			fakrat.push_back(f[i]);
			prorat.push_back(p[i]);
			fakerr.push_back(ferr[i]);
			proerr.push_back(perr[i]);
		}
	}

	// for (int i = 0; i < ne+nm; ++i) {
	// 	if (i < npe) {
	// 		ind.push_back(1);
	// 		fakrat.push_back(FakeRatio(1, vpt[i], veta[i]) );
	// 		prorat.push_back(PromptRatio(1, vpt[i], veta[i]) );
	// 		fakerr.push_back(FakeRatioError(1, vpt[i], veta[i]) );
	// 		proerr.push_back(PromptRatioError(1, vpt[i], veta[i]) );
	// 	} else if (i < ne) {
	// 		ind.push_back(0); 
	// 		fakrat.push_back(FakeRatio(1, vpt[i], veta[i]) );
	// 		prorat.push_back(PromptRatio(1, vpt[i], veta[i]) );
	// 		fakerr.push_back(FakeRatioError(1, vpt[i], veta[i]) );
	// 		proerr.push_back(PromptRatioError(1, vpt[i], veta[i]) );
	// 	} else if (i < ne+npm) {
	// 		ind.push_back(1);
	// 		fakrat.push_back(FakeRatio(2, vpt[i], veta[i]) );
	// 		prorat.push_back(PromptRatio(2, vpt[i], veta[i]) );
	// 		fakerr.push_back(FakeRatioError(2, vpt[i], veta[i]) );
	// 		proerr.push_back(PromptRatioError(2, vpt[i], veta[i]) );
	// 	} else {
	// 		ind.push_back(0);
	// 		fakrat.push_back(FakeRatio(2, vpt[i], veta[i]) );
	// 		prorat.push_back(PromptRatio(2, vpt[i], veta[i]) );
	// 		fakerr.push_back(FakeRatioError(2, vpt[i], veta[i]) );
	// 		proerr.push_back(PromptRatioError(2, vpt[i], veta[i]) );
	// 	}
	// //     if(verbose > 5) cout << " i = " << i << ", ind = " << ind[i] 
	// //	  << ", pT = " << vpt[i] << ", eta = " << veta[i] 
	// //	  << ", fakrat = " << fakrat[i] << ", prorat = " << prorat[i] << endl;
	// }

	int ipe = ipeinit;
	int ife = ifeinit;
	int ipemv = ipeinit;
	int ifemv = ifeinit;
	int ipm = ipminit;
	int ifm = ifminit;
	int ipmmv = ipminit;
	int ifmmv = ifminit;

	int icyce = 0;
	int icycm = 3;
	int nowmu = 0;
	if (ne == 0) nowmu = 1;
	wgtie = 0.;
	wgtim = 0.;
	wgt = 0.;

// sum over all the weights wij of the various FPconfigs
	while (1){

		// for (int i = 0; i < ne; ++i) {
		// 	if(verbose > 5) cout << " " << ind[i];
		// }
		// if(verbose > 5) cout << " | ";
		// for (int i = ne; i < ne+nm; ++i) {
		// 	if(verbose > 5) cout << " " << ind[i];
		// }
		// if(verbose > 5) cout << endl;
		// 

		// compute the weight for this combination, 
		//   multiplying all leptons of a given FPconfig
		double wgtnewe = 0., wgtnewm = 0.;
		double dwgtnewdpe = 0., dwgtnewdfe = 0.;
		double dwgtnewdpm = 0., dwgtnewdfm = 0.;
		double cpdpe = 0., cfdfe = 0.;
		double cpdpm = 0., cfdfm = 0.;
		for (int i = 0; i < ne+nm; ++i) {
			double p = prorat[i];
			double f = fakrat[i];
			double dp = proerr[i];
			double df = fakerr[i];

			if (i < ne) { // Electron
				if (wgtnewe != 0.) wgtnewe /= (p - f);
				else wgtnewe = 1. / (p - f);
				if (ind[i] && !ipass[i]) {
					wgtnewe *= f;
					if (inTight) wgtnewe *= p;
					if (p > 0.) cpdpe += (-1./(p-f) + 1./p);
					else cpdpe = 0.;
					if (f > 0.) cfdfe += (1./(p-f) + 1./f);
					else cfdfe = 0.;
				}
				if (ind[i] && ipass[i]) {
					wgtnewe *= (1. - f);
					if (inTight) wgtnewe *= p;
					if (p > 0.) cpdpe += (-1./(p-f) + 1./p);
					else cpdpe = 0.;
					cfdfe += (1./(p-f) - 1./(1.-f));
				}
				if (!ind[i] && !ipass[i]) {
					wgtnewe *= p;
					if (inTight) wgtnewe *= f;
					if (p > 0.) cpdpe += (-1./(p-f) + 1./p);
					else cpdpe = 0.;
					if (f > 0.) cfdfe += (1./(p-f) + 1./f);
					else cfdfe = 0.;
				}
				if (!ind[i] && ipass[i]) {
					wgtnewe *= (1. - p);
					if (inTight) wgtnewe *= f;
					if (p < 1.) cpdpe += (-1./(p-f) - 1./(1.-p));
					else cpdpe = 0.;
					if (f > 0.) cfdfe += (1./(p-f) + 1./f);
					else cfdfe = 0.;
				}
			} else { // Muons
				if (wgtnewm != 0.) wgtnewm /= (p - f);
				else wgtnewm = 1. / (p - f);
				if (ind[i] && !ipass[i]) {
					wgtnewm *= f;
					if (inTight) wgtnewm *= p;
					if (p > 0.) cpdpm += (-1./(p-f) + 1./p);
					else cpdpm = 0.;
					if (f > 0.) cfdfm += (1./(p-f) + 1./f);
					else cfdfm = 0.;
				}
				if (ind[i] && ipass[i]) {
					wgtnewm *= (1. - f);
					if (inTight) wgtnewm *= p;
					if (p > 0.) cpdpm += (-1./(p-f) + 1./p);
					else cpdpm = 0.;
					cfdfm += (1./(p-f) - 1./(1.-f));
				}
				if (!ind[i] && !ipass[i]) {
					wgtnewm *= p;
					if (inTight) wgtnewm *= f;
					if (p > 0.) cpdpm += (-1./(p-f) + 1./p);
					else cpdpm = 0.;
					if (f > 0.) cfdfm += (1./(p-f) + 1./f);
					else cfdfm = 0.;
				}
				if (!ind[i] && ipass[i]) {
					wgtnewm *= (1. - p);
					if (inTight) wgtnewm *= f;
					if (p < 1.) cpdpm += (-1./(p-f) - 1./(1.-p));
					else cpdpm = 0.;
					if (f > 0.) cfdfm += (1./(p-f) + 1./f);
					else cfdfm = 0.;
				}
			}
			dwgtnewdpe = wgtnewe * cpdpe * dp;
			dwgtnewdfe = wgtnewe * cfdfe * df;
			dwgtnewdpm = wgtnewm * cpdpm * dp;
			dwgtnewdfm = wgtnewm * cfdfm * df;
			if (verbose > 6) {
				if (i < ne) {
					if(verbose > 5) cout << " " << i << " wgtnewe = " << wgtnewe  << ", cpdpe = " << cpdpe << ", cfdfe = " << cfdfe
						<< ", dwgtnewdpe = " << dwgtnewdpe << ", dwgtnewdfe = " << dwgtnewdfe << endl;
				} else {
					if(verbose > 5) cout << " " << i << " wgtnewm = " << wgtnewm  << ", cpdpm = " << cpdpm << ", cfdfm = " << cfdfm
						<< ", dwgtnewdpm = " << dwgtnewdpm << ", dwgtnewdfm = " << dwgtnewdfm << endl;
				}
			}
		}
		wgtie += wgtnewe;
		wgtim += wgtnewm;
		double wgtnew;
		if (wgtnewe == 0.) wgtnew = wgtnewm;
		else if (wgtnewm == 0.) wgtnew = wgtnewe;
		else wgtnew = wgtnewe * wgtnewm;
		wgt += wgtnew;
		dwgtdpe += dwgtnewdpe;
		dwgtdfe += dwgtnewdfe;
		dwgtdpm += dwgtnewdpm;
		dwgtdfm += dwgtnewdfm;
		dwgtdptope.push_back(dwgtnewdpe); 
		dwgtdftope.push_back(dwgtnewdfe);
		dwgtdptopm.push_back(dwgtnewdpm); 
		dwgtdftopm.push_back(dwgtnewdfm);

	// generate next combination for electrons
		if (icyce == 1 && !nowmu) {
		//       if(verbose > 5) cout << " cycle over the fakes, ife = " << ife << endl;
			int temp = ind[ipeinit+1];
			for (int i = ipeinit+1; i < ifeinit+1; ++i) {
				ind[i] = ind[i+1];
			}
			ind[ifeinit] = temp;
			if (ife > ipeinit + 1 && ifemv > ipeinit) {
				ife--;
			} else {
				ife = ifeinit;
				icyce = 2;
			}
		}

		if (icyce == 2 && !nowmu) {
		//       if(verbose > 5) cout << " cycle over the prompts, ipe = " << ipe << endl;
			int temp = ind[0];
			for (int i = 0; i < ipeinit+1; ++i) {
				ind[i] = ind[i+1];
			}
			ind[ipeinit] = temp;
			if (ipe > 0 && ipemv >= 0) {
				ipe--;
				icyce = 1;
			} else {
				ipe = ipeinit;
				icyce = 0;
			}
		}

		if (icyce == 0 && !nowmu) {
			if (ipemv >= 0 && ifemv > ipeinit) {
	//          if(verbose > 5) cout << " move " << ipemv << " " << ifemv << endl;
				ind[ifemv] = 1;
				ind[ipemv] = 0;
				ipemv--;
				ifemv--;
				icyce = 1;
			} else {
				for (int i = 0; i < ne; ++i) {
					if (i < npe) ind[i] = 1;
					else ind[i] = 0; 
				}
				ipe = ipeinit;
				ife = ifeinit;
				ipemv = ipeinit;
				ifemv = ifeinit;
	//          if(verbose > 5) cout << " end cycling over electrons " << endl;
				nowmu = 1;
			}
		}

	// generate next combination for muons
		if (icycm == 4 && nowmu) {
		//       if(verbose > 5) cout << " cycle over the fakes, ifm = " << ifm << endl;
			int temp = ind[ipminit+1];
			for (int i = ipminit+1; i < ifminit+1; ++i) {
				ind[i] = ind[i+1];
			}
			ind[ifminit] = temp;
			if (ifm > ipminit + 1 && ifmmv > ipminit) {
				ifm--;
				if (ne > 0) nowmu = 0;
			} else {
				ifm = ifminit;
				icycm = 5;
			}
		}

		if (icycm == 5 && nowmu) {
		//       if(verbose > 5) cout << " cycle over the prompts, ipm = " << ipm << endl;
			int temp = ind[ne];
			for (int i = ne; i < ipminit+1; ++i) {
				ind[i] = ind[i+1];
			}
			ind[ipminit] = temp;
			if (ipm > ne && ipmmv >= ne) {
				ipm--;
				icycm = 4;
				if (ne > 0) nowmu = 0;
			} else {
				ipm = ipminit;
				icycm = 3;
			}
		}

		if (icycm == 3 && nowmu) {
			if (ipmmv >= ne && ifmmv > ipminit) {
	//          if(verbose > 5) cout << " move " << ipmmv << " " << ifmmv << endl;
				ind[ifmmv] = 1;
				ind[ipmmv] = 0;
				ipmmv--;
				ifmmv--;
				if (ne > 0) {
					icyce = 0;
					nowmu = 0;
				}
				icycm = 4;
			} else {
	//          if(verbose > 5) cout << " end cycling over muons " << endl;
				break;
			}
		}

	}

	if (verbose > 6) {
		if(verbose > 5) cout << " wgt = " << wgt << ", dwgtdpe = " << dwgtdpe << ", dwgtdfe = " << dwgtdfe
			<< ", dwgtdpm = " << dwgtdpm << ", dwgtdfm = " << dwgtdfm << endl;
	}

	ind.clear();
	fakrat.clear(); 
	prorat.clear();
	fakerr.clear(); 
	proerr.clear();

	return wgt;
}

// *************************************************************************
double FPRatios::FakeRatio (int emu, double pt, double eta){
// Defines the fake ratio
// emu = 1 for electrons
//     = 2 for muons
	if(emu == 1){ // Electrons
		if(fElefRatio == NULL || fElefRatio->GetEntries() == 0){
			// Return fixed value if histogram is not set
			return fFRatioe1;
		}
		return fElefRatio->GetBinContent(fElefRatio->FindBin(pt, eta));
	}

	if(emu == 2){ // Muons
		if(fMufRatio == NULL || fMufRatio->GetEntries() == 0){
			// Return fixed value if histogram is not set
			return fFRatiom1;
		}
		return fMufRatio->GetBinContent(fMufRatio->FindBin(pt, eta));
	}
}

// *************************************************************************
double FPRatios::FakeRatioError (int emu, double pt, double eta){
// Defines the fake ratio error
// emu = 1 for electrons
//     = 2 for muons
	if(emu == 1){ // Electrons
		if(fElefRatio == NULL || fElefRatio->GetEntries() == 0){
			// Return fixed value if histogram is not set
			return fFRatioErre1;
		}
		return fElefRatio->GetBinError(fElefRatio->FindBin(pt, eta));
	}

	if(emu == 2){ // Muons
		if(fMufRatio == NULL || fMufRatio->GetEntries() == 0){
			// Return fixed value if histogram is not set
			return fFRatioErrm1;
		}
		return fMufRatio->GetBinError(fMufRatio->FindBin(pt, eta));
	}
}

// *************************************************************************
double FPRatios::PromptRatio (int emu, double pt, double eta){
// Defines the prompt ratio
// emu = 1 for electrons
//     = 2 for muons
// Dummy function
	if(emu == 1){ // Electrons
		if(fElepRatio == NULL || fElepRatio->GetEntries() == 0){
			// Return fixed value if histogram is not set
			return fPRatioe1;
		}
		return fElepRatio->GetBinContent(fElepRatio->FindBin(pt, eta));
	}

	if(emu == 2){ // Muons
		if(fMupRatio == NULL || fMupRatio->GetEntries() == 0){
			// Return fixed value if histogram is not set
			return fPRatiom1;
		}
		return fMupRatio->GetBinContent(fMupRatio->FindBin(pt, eta));
	}
}

// *************************************************************************
double FPRatios::PromptRatioError (int emu, double pt, double eta){
// Defines the prompt ratio error
// emu = 1 for electrons
//     = 2 for muons
	if(emu == 1){ // Electrons
		if(fElepRatio == NULL || fElepRatio->GetEntries() == 0){
			// Return fixed value if histogram is not set
			return fPRatioErre1;
		}
		return fElepRatio->GetBinError(fElepRatio->FindBin(pt, eta));
	}

	if(emu == 2){ // Muons
		if(fMupRatio == NULL || fMupRatio->GetEntries() == 0){
			// Return fixed value if histogram is not set
			return fPRatioErrm1;
		}
		return fMupRatio->GetBinError(fMupRatio->FindBin(pt, eta));
	}
}

