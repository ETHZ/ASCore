/**************************************************************************************
 * Class to calculate di-lepton fake rate predictions                                 *
 *                                                                                    *
 * Detailed description of underlying methods in AN-10-261                            *
 * (c) Jun 2011, Benjamin Stieger, stiegerb@cern.ch                                   *
 *************************************************************************************/

#include <cmath>
#include <vector>
#include <iostream>
#include <iomanip>

#include "helper/FakeRatios.hh"
#include "TMath.h"
#include "TRandom3.h"

using namespace std;

FakeRatios::FakeRatios() : fVerbose(0), fNToyMCs(100), fAddESyst(0.0) {
	// Initialize arrays to avoid segfaults
	fMMNtl[0] = -1.; fMMNtl[1] = -1.; fMMNtl[2] = -1.; 
	fEENtl[0] = -1.; fEENtl[1] = -1.; fEENtl[2] = -1.; 
	fEMNtl[0] = -1.; fEMNtl[1] = -1.; fEMNtl[2] = -1.; fEMNtl[3] = -1.; 
	
	fMFRatio[0] = 0.; fMFRatio[1] = 0.; 
	fMPRatio[0] = 1.; fMPRatio[1] = 0.; 
	fEFRatio[0] = 0.; fEFRatio[1] = 0.; 
	fEPRatio[0] = 1.; fEPRatio[1] = 0.; 
}

FakeRatios::~FakeRatios(){}

//____________________________________________________________________________________
// Input
void FakeRatios::setMMNtl(float Ntt, float Ntl, float Nll){
	fMMNtl[0] = Ntt;
	fMMNtl[1] = Ntl;
	fMMNtl[2] = Nll;
}
void FakeRatios::setEENtl(float Ntt, float Ntl, float Nll){
	fEENtl[0] = Ntt;
	fEENtl[1] = Ntl;
	fEENtl[2] = Nll;
}
void FakeRatios::setEMNtl(float Ntt, float Ntl, float Nlt, float Nll){
	fEMNtl[0] = Ntt;
	fEMNtl[1] = Ntl;
	fEMNtl[2] = Nlt;
	fEMNtl[3] = Nll;
}

void FakeRatios::setMFRatio(float ratio, float error){
	fMFRatio[0] = ratio;
	fMFRatio[1] = error;
}
void FakeRatios::setMPRatio(float ratio, float error){
	fMPRatio[0] = ratio;
	fMPRatio[1] = error;	
}
void FakeRatios::setEFRatio(float ratio, float error){
	fEFRatio[0] = ratio;
	fEFRatio[1] = error;
}
void FakeRatios::setEPRatio(float ratio, float error){
	fEPRatio[0] = ratio;
	fEPRatio[1] = error;
}

//____________________________________________________________________________________
// Output and combinations
float FakeRatios::getMMNpp(){
	return getNpp(           fMMNtl[0], 0.5*fMMNtl[1], 0.5*fMMNtl[1], fMMNtl[2], fMFRatio[0], fMFRatio[0], fMPRatio[0], fMPRatio[0]);
}
float FakeRatios::getMMNppEStat(){
	return getNppEStat(      fMMNtl[0], 0.5*fMMNtl[1], 0.5*fMMNtl[1], fMMNtl[2], fMFRatio[0], fMFRatio[0], fMPRatio[0], fMPRatio[0]);
}
float FakeRatios::getMMNppESyst(){
	float fromtoys = getESystFromToys2(fMMNtl[0], 0.5*fMMNtl[1], 0.5*fMMNtl[1], fMMNtl[2], fMFRatio[0], fMFRatio[0], fMPRatio[0], fMPRatio[0], fMFRatio[1], fMFRatio[1], fMPRatio[1], fMPRatio[1], &FakeRatios::getNpp);
	float addsyst = fAddESyst * getMMNpp();
	return sqrt(fromtoys + addsyst*addsyst);
}
float FakeRatios::getMMNppETot(){
	return sqrt(getMMNppEStat()*getMMNppEStat() + getMMNppESyst()*getMMNppESyst());
}

float FakeRatios::getMMNpf(){
	return 2.*getNpf(           fMMNtl[0], 0.5*fMMNtl[1], 0.5*fMMNtl[1], fMMNtl[2], fMFRatio[0], fMFRatio[0], fMPRatio[0], fMPRatio[0]);
}
float FakeRatios::getMMNpfEStat(){
	return 2.*getNpfEStat(      fMMNtl[0], 0.5*fMMNtl[1], 0.5*fMMNtl[1], fMMNtl[2], fMFRatio[0], fMFRatio[0], fMPRatio[0], fMPRatio[0]);
}
float FakeRatios::getMMNpfESyst(){
	float fromtoys = 4.*getESystFromToys2(fMMNtl[0], 0.5*fMMNtl[1], 0.5*fMMNtl[1], fMMNtl[2], fMFRatio[0], fMFRatio[0], fMPRatio[0], fMPRatio[0], fMFRatio[1], fMFRatio[1], fMPRatio[1], fMPRatio[1], &FakeRatios::getNpf);
	float addsyst = fAddESyst * getMMNpf();
	return sqrt(fromtoys + addsyst*addsyst);	
}
float FakeRatios::getMMNpfETot(){
	return sqrt(getMMNpfEStat()*getMMNpfEStat() + getMMNpfESyst()*getMMNpfESyst());
}

float FakeRatios::getMMNff(){
	return getNff(           fMMNtl[0], 0.5*fMMNtl[1], 0.5*fMMNtl[1], fMMNtl[2], fMFRatio[0], fMFRatio[0], fMPRatio[0], fMPRatio[0]);
}
float FakeRatios::getMMNffEStat(){
	return getNffEStat(      fMMNtl[0], 0.5*fMMNtl[1], 0.5*fMMNtl[1], fMMNtl[2], fMFRatio[0], fMFRatio[0], fMPRatio[0], fMPRatio[0]);
}
float FakeRatios::getMMNffESyst(){
	float fromtoys = getESystFromToys2(fMMNtl[0], 0.5*fMMNtl[1], 0.5*fMMNtl[1], fMMNtl[2], fMFRatio[0], fMFRatio[0], fMPRatio[0], fMPRatio[0], fMFRatio[1], fMFRatio[1], fMPRatio[1], fMPRatio[1], &FakeRatios::getNff);
	float addsyst = fAddESyst * getMMNff();
	return sqrt(fromtoys + addsyst*addsyst);	
}
float FakeRatios::getMMNffETot(){
	return sqrt(getMMNffEStat()*getMMNffEStat() + getMMNffESyst()*getMMNffESyst());
}

//____________________________________________________________________________________
float FakeRatios::getEENpp(){
	return getNpp(   fEENtl[0], 0.5*fEENtl[1], 0.5*fEENtl[1], fEENtl[2], fEFRatio[0], fEFRatio[0], fEPRatio[0], fEPRatio[0]);
}
float FakeRatios::getEENppEStat(){
	return getNppEStat(      fEENtl[0], 0.5*fEENtl[1], 0.5*fEENtl[1], fEENtl[2], fEFRatio[0], fEFRatio[0], fEPRatio[0], fEPRatio[0]);
}
float FakeRatios::getEENppESyst(){
	float fromtoys = getESystFromToys2(fEENtl[0], 0.5*fEENtl[1], 0.5*fEENtl[1], fEENtl[2], fEFRatio[0], fEFRatio[0], fEPRatio[0], fEPRatio[0], fEFRatio[1], fEFRatio[1], fEPRatio[1], fEPRatio[1], &FakeRatios::getNpp);
	float addsyst = fAddESyst * getEENpp();
	return sqrt(fromtoys + addsyst*addsyst);	
}
float FakeRatios::getEENppETot(){
	return sqrt(getEENppEStat()*getEENppEStat() + getEENppESyst()*getEENppESyst());
}

float FakeRatios::getEENpf(){
	return 2.*getNfp(fEENtl[0], 0.5*fEENtl[1], 0.5*fEENtl[1], fEENtl[2], fEFRatio[0], fEFRatio[0], fEPRatio[0], fEPRatio[0]);
}
float FakeRatios::getEENpfEStat(){
	return 2.*getNfpEStat(      fEENtl[0], 0.5*fEENtl[1], 0.5*fEENtl[1], fEENtl[2], fEFRatio[0], fEFRatio[0], fEPRatio[0], fEPRatio[0]);
}
float FakeRatios::getEENpfESyst(){
	float fromtoys =  4.*getESystFromToys2(fEENtl[0], 0.5*fEENtl[1], 0.5*fEENtl[1], fEENtl[2], fEFRatio[0], fEFRatio[0], fEPRatio[0], fEPRatio[0], fEFRatio[1], fEFRatio[1], fEPRatio[1], fEPRatio[1], &FakeRatios::getNfp);
	float addsyst = fAddESyst * getEENpf();
	return sqrt(fromtoys + addsyst*addsyst);	
}
float FakeRatios::getEENpfETot(){
	return sqrt(getEENpfEStat()*getEENpfEStat() + getEENpfESyst()*getEENpfESyst());
}

float FakeRatios::getEENff(){
	return getNff(   fEENtl[0], 0.5*fEENtl[1], 0.5*fEENtl[1], fEENtl[2], fEFRatio[0], fEFRatio[0], fEPRatio[0], fEPRatio[0]);
}
float FakeRatios::getEENffEStat(){
	return getNffEStat(      fEENtl[0], 0.5*fEENtl[1], 0.5*fEENtl[1], fEENtl[2], fEFRatio[0], fEFRatio[0], fEPRatio[0], fEPRatio[0]);
}
float FakeRatios::getEENffESyst(){
	float fromtoys = getESystFromToys2(fEENtl[0], 0.5*fEENtl[1], 0.5*fEENtl[1], fEENtl[2], fEFRatio[0], fEFRatio[0], fEPRatio[0], fEPRatio[0], fEFRatio[1], fEFRatio[1], fEPRatio[1], fEPRatio[1], &FakeRatios::getNff);
	float addsyst = fAddESyst * getEENff();
	return sqrt(fromtoys + addsyst*addsyst);	
}
float FakeRatios::getEENffETot(){
	return sqrt(getEENffEStat()*getEENffEStat() + getEENffESyst()*getEENffESyst());
}

//____________________________________________________________________________________
float FakeRatios::getEMNpp(){
	return getNpp(           fEMNtl[0], fEMNtl[1], fEMNtl[2], fEMNtl[3], fMFRatio[0], fEFRatio[0], fMPRatio[0], fEPRatio[0]);
}
float FakeRatios::getEMNppEStat(){
	return getNppEStat(      fEMNtl[0], fEMNtl[1], fEMNtl[2], fEMNtl[3], fMFRatio[0], fEFRatio[0], fMPRatio[0], fEPRatio[0]);
}
float FakeRatios::getEMNppESyst(){
	float fromtoys = getESystFromToys2(fEMNtl[0], fEMNtl[1], fEMNtl[2], fEMNtl[3], fMFRatio[0], fEFRatio[0], fMPRatio[0], fEPRatio[0], fMFRatio[1], fEFRatio[1], fMPRatio[1], fEPRatio[1], &FakeRatios::getNpp);	
	float addsyst = fAddESyst * getEMNpp();
	return sqrt(fromtoys + addsyst*addsyst);	
}
float FakeRatios::getEMNppETot(){
	return sqrt(getEMNppEStat()*getEMNppEStat() + getEMNppESyst()*getEMNppESyst());
}

float FakeRatios::getEMNpf(){
	return getNpf(           fEMNtl[0], fEMNtl[1], fEMNtl[2], fEMNtl[3], fMFRatio[0], fEFRatio[0], fMPRatio[0], fEPRatio[0]);
}
float FakeRatios::getEMNpfEStat(){
	return getNpfEStat(      fEMNtl[0], fEMNtl[1], fEMNtl[2], fEMNtl[3], fMFRatio[0], fEFRatio[0], fMPRatio[0], fEPRatio[0]);
}
float FakeRatios::getEMNpfESyst(){
	float fromtoys = getESystFromToys2(fEMNtl[0], fEMNtl[1], fEMNtl[2], fEMNtl[3], fMFRatio[0], fEFRatio[0], fMPRatio[0], fEPRatio[0], fMFRatio[1], fEFRatio[1], fMPRatio[1], fEPRatio[1], &FakeRatios::getNpf);
	float addsyst = fAddESyst * getEMNpf();
	return sqrt(fromtoys + addsyst*addsyst);	
}
float FakeRatios::getEMNpfETot(){
	return sqrt(getEMNpfEStat()*getEMNpfEStat() + getEMNpfESyst()*getEMNpfESyst());
}

float FakeRatios::getEMNfp(){
	return getNfp(           fEMNtl[0], fEMNtl[1], fEMNtl[2], fEMNtl[3], fMFRatio[0], fEFRatio[0], fMPRatio[0], fEPRatio[0]);
}
float FakeRatios::getEMNfpEStat(){
	return getNfpEStat(      fEMNtl[0], fEMNtl[1], fEMNtl[2], fEMNtl[3], fMFRatio[0], fEFRatio[0], fMPRatio[0], fEPRatio[0]);
}
float FakeRatios::getEMNfpESyst(){
	float fromtoys = getESystFromToys2(fEMNtl[0], fEMNtl[1], fEMNtl[2], fEMNtl[3], fMFRatio[0], fEFRatio[0], fMPRatio[0], fEPRatio[0], fMFRatio[1], fEFRatio[1], fMPRatio[1], fEPRatio[1], &FakeRatios::getNfp);
	float addsyst = fAddESyst * getEMNfp();
	return sqrt(fromtoys + addsyst*addsyst);	
}
float FakeRatios::getEMNfpETot(){
	return sqrt(getEMNfpEStat()*getEMNfpEStat() + getEMNfpESyst()*getEMNfpESyst());
}

float FakeRatios::getEMSingleFakes(){
	return getEMNpf() + getEMNfp();
}
float FakeRatios::getEMSingleEStat(){
	return getNfpNpfSumEStat(fEMNtl[0], fEMNtl[1], fEMNtl[2], fEMNtl[3], fMFRatio[0], fEFRatio[0], fMPRatio[0], fEPRatio[0]);
}
float FakeRatios::getEMSingleESyst(){
	float fromtoys = getESystFromToys2(fEMNtl[0], fEMNtl[1], fEMNtl[2], fEMNtl[3], fMFRatio[0], fEFRatio[0], fMPRatio[0], fEPRatio[0], fMFRatio[1], fEFRatio[1], fMPRatio[1], fEPRatio[1], &FakeRatios::getNfpNpfSum);
	float addsyst = fAddESyst * (getEMNfp() + getEMNpf());
	return sqrt(fromtoys + addsyst*addsyst);
}
float FakeRatios::getEMSingleETot(){
	return sqrt(getEMSingleEStat()*getEMSingleEStat() + getEMSingleESyst()*getEMSingleESyst());	
}

float FakeRatios::getEMNff(){
	return getNff(           fEMNtl[0], fEMNtl[1], fEMNtl[2], fEMNtl[3], fMFRatio[0], fEFRatio[0], fMPRatio[0], fEPRatio[0]);
}
float FakeRatios::getEMNffEStat(){
	return getNffEStat(      fEMNtl[0], fEMNtl[1], fEMNtl[2], fEMNtl[3], fMFRatio[0], fEFRatio[0], fMPRatio[0], fEPRatio[0]);
}
float FakeRatios::getEMNffESyst(){
	float fromtoys = getESystFromToys2(fEMNtl[0], fEMNtl[1], fEMNtl[2], fEMNtl[3], fMFRatio[0], fEFRatio[0], fMPRatio[0], fEPRatio[0], fMFRatio[1], fEFRatio[1], fMPRatio[1], fEPRatio[1], &FakeRatios::getNff);
	float addsyst = fAddESyst * getEMNff();
	return sqrt(fromtoys + addsyst*addsyst);	
}
float FakeRatios::getEMNffETot(){
	return sqrt(getEMNffEStat()*getEMNffEStat() + getEMNffESyst()*getEMNffESyst());
}

//____________________________________________________________________________________
// Single channel combinations
float FakeRatios::getMMTotFakes(){
	return getNfpNpfNffSum(fMMNtl[0], 0.5*fMMNtl[1], 0.5*fMMNtl[1], fMMNtl[2], fMFRatio[0], fMFRatio[0], fMPRatio[0], fMPRatio[0]);
}
float FakeRatios::getMMTotEStat(){
	return getNfpNpfNffSumEStat(fMMNtl[0], 0.5*fMMNtl[1], 0.5*fMMNtl[1], fMMNtl[2], fMFRatio[0], fMFRatio[0], fMPRatio[0], fMPRatio[0]);
}
float FakeRatios::getMMTotESyst(){
	float fromtoys = sqrt(getESystFromToys2(fMMNtl[0], 0.5*fMMNtl[1], 0.5*fMMNtl[1], fMMNtl[2], fMFRatio[0], fMFRatio[0], fMPRatio[0], fMPRatio[0], fMFRatio[1], fMFRatio[1], fMPRatio[1], fMPRatio[1], &FakeRatios::getNfpNpfNffSum));
	float addsyst = fAddESyst * getMMTotFakes();
	return sqrt(fromtoys + addsyst*addsyst);
}
float FakeRatios::getEETotFakes(){
	return getNfpNpfNffSum(fEENtl[0], 0.5*fEENtl[1], 0.5*fEENtl[1], fEENtl[2], fEFRatio[0], fEFRatio[0], fEPRatio[0], fEPRatio[0]);
}
float FakeRatios::getEETotEStat(){
	return getNfpNpfNffSumEStat(fEENtl[0], 0.5*fEENtl[1], 0.5*fEENtl[1], fEENtl[2], fEFRatio[0], fEFRatio[0], fEPRatio[0], fEPRatio[0]);
}
float FakeRatios::getEETotESyst(){
	float fromtoys = sqrt(getESystFromToys2(fEENtl[0], 0.5*fEENtl[1], 0.5*fEENtl[1], fEENtl[2], fEFRatio[0], fEFRatio[0], fEPRatio[0], fEPRatio[0], fEFRatio[1], fEFRatio[1], fEPRatio[1], fEPRatio[1], &FakeRatios::getNfpNpfNffSum));
	float addsyst = fAddESyst * getEETotFakes();
	return sqrt(fromtoys + addsyst*addsyst);
}
float FakeRatios::getEMTotFakes(){
	return getNfpNpfNffSum(fEMNtl[0], fEMNtl[1], fEMNtl[2], fEMNtl[3], fMFRatio[0], fEFRatio[0], fMPRatio[0], fEPRatio[0]);
}
float FakeRatios::getEMTotEStat(){
	return getNfpNpfNffSumEStat(fEMNtl[0], fEMNtl[1], fEMNtl[2], fEMNtl[3], fMFRatio[0], fEFRatio[0], fMPRatio[0], fEPRatio[0]);
}
float FakeRatios::getEMTotESyst(){
	float fromtoys = sqrt(getESystFromToys2(fEMNtl[0], fEMNtl[1], fEMNtl[2], fEMNtl[3], fMFRatio[0], fEFRatio[0], fMPRatio[0], fEPRatio[0], fMFRatio[1], fEFRatio[1], fMPRatio[1], fEPRatio[1], &FakeRatios::getNfpNpfNffSum));
	float addsyst = fAddESyst * getEMTotFakes();
	return sqrt(fromtoys + addsyst*addsyst);
}

//____________________________________________________________________________________
// Multi channel combination
float FakeRatios::getTotFakes(){
	// Simple sum
	return getMMTotFakes() + getEETotFakes() + getEMTotFakes();
}
float FakeRatios::getTotSingleFakes(){
	return getMMNpf() + getEENpf() + getEMNpf() + getEMNfp();
}
float FakeRatios::getTotDoubleFakes(){
	return getMMNff() + getEENff() + getEMNff();
}
float FakeRatios::getTotEStat(){
	// Control yields in different channels are independent!
	float mm = getMMTotEStat();
	float ee = getEETotEStat();
	float em = getEMTotEStat();
	return sqrt( mm*mm + ee*ee + em*em );
}
float FakeRatios::getTotSingleEStat(){
	// Control yields in different channels are independent!
	float mm = getMMNpfEStat();
	float ee = getEENpfEStat();
	float em = getEMSingleEStat();
	return sqrt( mm*mm + ee*ee + em*em );
}
float FakeRatios::getTotDoubleEStat(){
	// Control yields in different channels are independent!
	float mm = getMMNffEStat();
	float ee = getEENffEStat();
	float em = getEMNffEStat();
	return sqrt( mm*mm + ee*ee + em*em );
}
float FakeRatios::getTotESyst(){
	// Assume errors of p and f are uncorrelated
	// Throw toys in a gaussian around f and p with df and dp as their sigmas
	// Distributions for f and p are cut off at 0 and 1
	float f1  = fMFRatio[0];
	float df1 = fMFRatio[1];
	float f2  = fEFRatio[0];
	float df2 = fEFRatio[1];
	float p1  = fMPRatio[0];
	float dp1 = fMPRatio[1];
	float p2  = fEPRatio[0];
	float dp2 = fEPRatio[1];	

	if(fVerbose > 2) cout << "FakeRatios::getTotESyst ..." << endl;
	TRandom3 *rand = new TRandom3();
	rand->SetSeed();
	vector<float> f_results;
	for(size_t i = 0; i < fNToyMCs; ++i){ // vary f1
		float f1_v = rand->Gaus(f1, df1);
		if(f1_v > 1. || f1_v < 0.){--i; continue;} // throw again if f<0 or f>1
		for(size_t j = 0; j < fNToyMCs; ++j){ // vary f2
			float f2_v = rand->Gaus(f2, df2);
			if(f2_v > 1. || f2_v < 0.){--j; continue;} // throw again if f<0 or f>1
			float result = getTotFakes(f1_v, f2_v, p1, p2);
			f_results.push_back(result);
			if(fVerbose > 2) cout << result << endl;
		}
	}
	float rms_f = TMath::RMS(f_results.begin(), f_results.end());
	if(fVerbose > 2) cout << " RMS = " << rms_f << endl;
	if(fVerbose > 2) cout << endl;
	vector<float> p_results;
	for(size_t i = 0; i < fNToyMCs; ++i){ // vary p1
		float p1_v = rand->Gaus(p1, dp1);
		if(p1_v > 1. || p1_v < 0.){--i; continue;} // throw again if p<0 or p>1
		for(size_t j = 0; j < fNToyMCs; ++j){ // vary p2
			float p2_v = rand->Gaus(p2, dp2);
			if(p2_v > 1. || p2_v < 0.){--j; continue;} // throw again if p<0 or p>1
			float result = getTotFakes(f1, f2, p1_v, p2_v);
			p_results.push_back(result);
			if(fVerbose > 2) cout << result << endl;
		}
	}
	float rms_p = TMath::RMS(p_results.begin(), p_results.end());
	if(fVerbose > 2) cout << " RMS = " << rms_p << endl;
	
	float fromtoys = rms_f*rms_f + rms_p*rms_p;
	float addsyst = fAddESyst * getTotFakes();
	return sqrt(fromtoys + addsyst*addsyst);
}
float FakeRatios::getTotSingleESyst(){
	// Assume errors of p and f are uncorrelated
	// Throw toys in a gaussian around f and p with df and dp as their sigmas
	// Distributions for f and p are cut off at 0 and 1
	float f1  = fMFRatio[0];
	float df1 = fMFRatio[1];
	float f2  = fEFRatio[0];
	float df2 = fEFRatio[1];
	float p1  = fMPRatio[0];
	float dp1 = fMPRatio[1];
	float p2  = fEPRatio[0];
	float dp2 = fEPRatio[1];	

	if(fVerbose > 2) cout << "FakeRatios::getTotESyst ..." << endl;
	TRandom3 *rand = new TRandom3();
	rand->SetSeed();
	vector<float> f_results;
	for(size_t i = 0; i < fNToyMCs; ++i){ // vary f1
		float f1_v = rand->Gaus(f1, df1);
		if(f1_v > 1. || f1_v < 0.){--i; continue;} // throw again if f<0 or f>1
		for(size_t j = 0; j < fNToyMCs; ++j){ // vary f2
			float f2_v = rand->Gaus(f2, df2);
			if(f2_v > 1. || f2_v < 0.){--j; continue;} // throw again if f<0 or f>1
			float result = getTotSingleFakes(f1_v, f2_v, p1, p2);
			f_results.push_back(result);
			if(fVerbose > 2) cout << result << endl;
		}
	}
	float rms_f = TMath::RMS(f_results.begin(), f_results.end());
	if(fVerbose > 2) cout << " RMS = " << rms_f << endl;
	if(fVerbose > 2) cout << endl;
	vector<float> p_results;
	for(size_t i = 0; i < fNToyMCs; ++i){ // vary p1
		float p1_v = rand->Gaus(p1, dp1);
		if(p1_v > 1. || p1_v < 0.){--i; continue;} // throw again if p<0 or p>1
		for(size_t j = 0; j < fNToyMCs; ++j){ // vary p2
			float p2_v = rand->Gaus(p2, dp2);
			if(p2_v > 1. || p2_v < 0.){--j; continue;} // throw again if p<0 or p>1
			float result = getTotSingleFakes(f1, f2, p1_v, p2_v);
			p_results.push_back(result);
			if(fVerbose > 2) cout << result << endl;
		}
	}
	float rms_p = TMath::RMS(p_results.begin(), p_results.end());
	if(fVerbose > 2) cout << " RMS = " << rms_p << endl;
	
	float fromtoys = rms_f*rms_f + rms_p*rms_p;
	float addsyst = fAddESyst * getTotSingleFakes();
	return sqrt(fromtoys + addsyst*addsyst);
}
float FakeRatios::getTotDoubleESyst(){
	// Assume errors of p and f are uncorrelated
	// Throw toys in a gaussian around f and p with df and dp as their sigmas
	// Distributions for f and p are cut off at 0 and 1
	float f1  = fMFRatio[0];
	float df1 = fMFRatio[1];
	float f2  = fEFRatio[0];
	float df2 = fEFRatio[1];
	float p1  = fMPRatio[0];
	float dp1 = fMPRatio[1];
	float p2  = fEPRatio[0];
	float dp2 = fEPRatio[1];	

	if(fVerbose > 2) cout << "FakeRatios::getTotESyst ..." << endl;
	TRandom3 *rand = new TRandom3();
	rand->SetSeed();
	vector<float> f_results;
	for(size_t i = 0; i < fNToyMCs; ++i){ // vary f1
		float f1_v = rand->Gaus(f1, df1);
		if(f1_v > 1. || f1_v < 0.){--i; continue;} // throw again if f<0 or f>1
		for(size_t j = 0; j < fNToyMCs; ++j){ // vary f2
			float f2_v = rand->Gaus(f2, df2);
			if(f2_v > 1. || f2_v < 0.){--j; continue;} // throw again if f<0 or f>1
			float result = getTotDoubleFakes(f1_v, f2_v, p1, p2);
			f_results.push_back(result);
			if(fVerbose > 2) cout << result << endl;
		}
	}
	float rms_f = TMath::RMS(f_results.begin(), f_results.end());
	if(fVerbose > 2) cout << " RMS = " << rms_f << endl;
	if(fVerbose > 2) cout << endl;
	vector<float> p_results;
	for(size_t i = 0; i < fNToyMCs; ++i){ // vary p1
		float p1_v = rand->Gaus(p1, dp1);
		if(p1_v > 1. || p1_v < 0.){--i; continue;} // throw again if p<0 or p>1
		for(size_t j = 0; j < fNToyMCs; ++j){ // vary p2
			float p2_v = rand->Gaus(p2, dp2);
			if(p2_v > 1. || p2_v < 0.){--j; continue;} // throw again if p<0 or p>1
			float result = getTotDoubleFakes(f1, f2, p1_v, p2_v);
			p_results.push_back(result);
			if(fVerbose > 2) cout << result << endl;
		}
	}
	float rms_p = TMath::RMS(p_results.begin(), p_results.end());
	if(fVerbose > 2) cout << " RMS = " << rms_p << endl;
	
	float fromtoys = rms_f*rms_f + rms_p*rms_p;
	float addsyst = fAddESyst * getTotDoubleFakes();
	return sqrt(fromtoys + addsyst*addsyst);
}
float FakeRatios::getTotETot(){
	return sqrt(getTotEStat()*getTotEStat() + getTotESyst()*getTotESyst());
}
float FakeRatios::getTotSingleETot(){
	return sqrt(getTotSingleEStat()*getTotSingleEStat() + getTotSingleESyst()*getTotSingleESyst());
}
float FakeRatios::getTotDoubleETot(){
	return sqrt(getTotDoubleEStat()*getTotDoubleEStat() + getTotDoubleESyst()*getTotDoubleESyst());
}

//____________________________________________________________________________________
// Helper methods for inside the toy loops
float FakeRatios::getTotFakes(float f1, float f2, float p1, float p2){
	// Need this with custom ratios inside the ESyst toy MC loops
	float mm = getNfpNpfNffSum(fMMNtl[0], 0.5*fMMNtl[1], 0.5*fMMNtl[1], fMMNtl[2], f1, f2, p1, p2);	
	float ee = getNfpNpfNffSum(fEENtl[0], 0.5*fEENtl[1], 0.5*fEENtl[1], fEENtl[2], f1, f2, p1, p2);
	float em = getNfpNpfNffSum(fEMNtl[0],     fEMNtl[1],     fEMNtl[2], fEMNtl[3], f1, f2, p1, p2);
	return mm + ee + em;
}
float FakeRatios::getTotSingleFakes(float f1, float f2, float p1, float p2){
	// Need this with custom ratios inside the ESyst toy MC loops
	float mm = getNfpNpfSum(fMMNtl[0], 0.5*fMMNtl[1], 0.5*fMMNtl[1], fMMNtl[2], f1, f2, p1, p2);	
	float ee = getNfpNpfSum(fEENtl[0], 0.5*fEENtl[1], 0.5*fEENtl[1], fEENtl[2], f1, f2, p1, p2);
	float em = getNfpNpfSum(fEMNtl[0],     fEMNtl[1],     fEMNtl[2], fEMNtl[3], f1, f2, p1, p2);
	return mm + ee + em;
}
float FakeRatios::getTotDoubleFakes(float f1, float f2, float p1, float p2){
	// Need this with custom ratios inside the ESyst toy MC loops
	float mm = getNff(fMMNtl[0], 0.5*fMMNtl[1], 0.5*fMMNtl[1], fMMNtl[2], f1, f2, p1, p2);	
	float ee = getNff(fEENtl[0], 0.5*fEENtl[1], 0.5*fEENtl[1], fEENtl[2], f1, f2, p1, p2);
	float em = getNff(fEMNtl[0],     fEMNtl[1],     fEMNtl[2], fEMNtl[3], f1, f2, p1, p2);
	return mm + ee + em;
}

//____________________________________________________________________________________
// Print
void FakeRatios::printOutput(){
	cout << "--------------------------------------------------------------------------------------------------------------------------" << endl;
	cout << " Output of FakeRatios class" << endl;
	cout << "--------------------------------------------------------------------------------------------------------------------------" << endl;
	cout << "                 ||            Mu/Mu            ||                   E/Mu                ||             E/E             ||" << endl;
	cout << " INPUT YIELDS    ||   Nt2   |   Nt1   |   Nt0   ||   Nt2   |   Nt10  |   Nt01  |   Nt0   ||   Nt2   |   Nt1   |   Nt0   ||" << endl;
	cout << "--------------------------------------------------------------------------------------------------------------------------" << endl;
	cout << setw(16) << " || ";
	cout << setw(7) << Form("%5.1f", fMMNtl[0] ) << " | ";
	cout << setw(7) << Form("%5.1f", fMMNtl[1] ) << " | ";
	cout << setw(7) << Form("%5.1f", fMMNtl[2] ) << " || ";
	cout << setw(7) << Form("%5.1f", fEMNtl[0] ) << " | ";
	cout << setw(7) << Form("%5.1f", fEMNtl[1] ) << " | ";
	cout << setw(7) << Form("%5.1f", fEMNtl[2] ) << " | ";
	cout << setw(7) << Form("%5.1f", fEMNtl[3] ) << " || ";
	cout << setw(7) << Form("%5.1f", fEENtl[0] ) << " | ";
	cout << setw(7) << Form("%5.1f", fEENtl[1] ) << " | ";
	cout << setw(7) << Form("%5.1f", fEENtl[2] ) << " || ";
	cout << endl;
	cout << "--------------------------------------------------------------------------------------------------------------------------" << endl;
	cout << "                 ||            Mu/Mu            ||                   E/Mu                ||             E/E             ||" << endl;
	cout << " PREDICTIONS     ||   Npp   |   Npf   |   Nff   ||   Npp   |   Npf   |   Nfp   |   Nff   ||   Npp   |   Npf   |   Nff   ||" << endl;
	cout << "--------------------------------------------------------------------------------------------------------------------------" << endl;
	cout << setw(16) << " || ";
	cout << setw(7) << Form("%5.1f", getMMNpp() ) << " | ";
	cout << setw(7) << Form("%5.1f", getMMNpf() ) << " | ";
	cout << setw(7) << Form("%5.1f", getMMNff() ) << " || ";
	cout << setw(7) << Form("%5.1f", getEMNpp() ) << " | ";
	cout << setw(7) << Form("%5.1f", getEMNpf() ) << " | ";
	cout << setw(7) << Form("%5.1f", getEMNfp() ) << " | ";
	cout << setw(7) << Form("%5.1f", getEMNff() ) << " || ";
	cout << setw(7) << Form("%5.1f", getEENpp() ) << " | ";
	cout << setw(7) << Form("%5.1f", getEENpf() ) << " | ";
	cout << setw(7) << Form("%5.1f", getEENff() ) << " || ";
	cout << endl;
	cout << "--------------------------------------------------------------------------------------------------------------------------" << endl;
}

/*************************************************************************************
###  ENGINE  #########################################################################
**************************************************************************************/
//____________________________________________________________________________________
// Dilepton formulas from AN-10-261
//  -- Note: these methods return yields in tight-tight region
float FakeRatios::getNpp(float Ntt, float Ntl, float Nlt, float Nll, float f1, float f2, float p1, float p2){
	// (f1-p1)(f2-p2) Npp   =  (f1-1)(f2-1) Ntt  +  (f1-1)f2 Ntl  +  f1(f2-1) Nlt  +  f1f2 Nll
	return p1*p2/((f1-p1)*(f2-p2)) * ( (f1-1.)*(f2-1.)*Ntt
	                                 + (f1-1.)* f2    *Ntl
	                                 +  f1    *(f2-1.)*Nlt
	                                 +  f1    * f2    *Nll );
}
float FakeRatios::getNpf(float Ntt, float Ntl, float Nlt, float Nll, float f1, float f2, float p1, float p2){
	// (f1-p1)(f2-p2) Npf   =  (f1-1)(1-p2) Ntt  -  (f1-1)p2 Ntl  +  f1(1-p2) Nlt  -  f1p2 Nll
	return p1*f2/((f1-p1)*(f2-p2)) * ( (f1-1.)*(1.-p2)*Ntt
	                                 - (f1-1.)*    p2 *Ntl
	                                 +  f1    *(1.-p2)*Nlt
	                                 -  f1    *    p2 *Nll );
}
float FakeRatios::getNfp(float Ntt, float Ntl, float Nlt, float Nll, float f1, float f2, float p1, float p2){
	// (f1-p1)(f2-p2) Nfp   =  (1-p1)(f2-1) Ntt  +  (1-p1)f2 Ntl  -  p1(f2-1) Nlt  -  p1f2 Nll
	return f1*p2/((f1-p1)*(f2-p2)) * ( (1.-p1)*(f2-1.)*Ntt
	                                 + (1.-p1)* f2    *Ntl
	                                 -     p1 *(f2-1.)*Nlt
	                                 -     p1 * f2    *Nll );
}
float FakeRatios::getNff(float Ntt, float Ntl, float Nlt, float Nll, float f1, float f2, float p1, float p2){
	// (f1-p1)(f2-p2) Nff   =  (1-p1)(1-p2) Ntt  -  (1-p1)p2 Ntl  -  p1(1-p2) Nlt  +  p1p2 Nll
	return f1*f2/((f1-p1)*(f2-p2)) * ( (1.-p1)*(1.-p2)*Ntt
	                                 - (1.-p1)*    p2 *Ntl
	                                 -     p1 *(1.-p2)*Nlt
	                                 +     p1 *    p2 *Nll );
}
float FakeRatios::getNfpNpfNffSum(float Ntt, float Ntl, float Nlt, float Nll, float f1, float f2, float p1, float p2){
	// This is equivalent to p1*f2*Npf + f1*p2*Nfp + f1*f2*Nff
	return 1./((f1-p1)*(f2-p2)) * ( (    p1*f2*(f1-1.)*(1.-p2) + f1*p2*(1.-p1)*(f2-1.) + f1*f2*(1.-p1)*(1.-p2))*Ntt
	                             +  (-1.*p1*f2*(f1-1.)*   p2   + f1*p2*(1.-p1)* f2     - f1*f2*(1.-p1)*    p2 )*Ntl
	                             +  (    p1*f2* f1    *(1.-p2) - f1*p2*    p1 *(f2-1.) - f1*f2*    p1 *(1.-p2))*Nlt
	                             +  (-1.*p1*f2* f1    *   p2   - f1*p2*    p1 * f2     + f1*f2*    p1 *    p2 )*Nll ); 
}
float FakeRatios::getNfpNpfSum(float Ntt, float Ntl, float Nlt, float Nll, float f1, float f2, float p1, float p2){
	// This is equivalent to p1*f2*Npf + f1*p2*Nfp
	return 1./((f1-p1)*(f2-p2)) * ( (    p1*f2*(f1-1.)*(1.-p2) + f1*p2*(1.-p1)*(f2-1.))*Ntt
	                             +  (-1.*p1*f2*(f1-1.)*   p2   + f1*p2*(1.-p1)* f2    )*Ntl
	                             +  (    p1*f2* f1    *(1.-p2) - f1*p2*    p1 *(f2-1.))*Nlt
	                             +  (-1.*p1*f2* f1    *   p2   - f1*p2*    p1 * f2    )*Nll ); 
}

//____________________________________________________________________________________
// Simple error propagation on above formulas for stat errors
float FakeRatios::getNppEStat(float Ntt, float Ntl, float Nlt, float Nll, float f1, float f2, float p1, float p2){
	return p1*p2/((f1-p1)*(f2-p2)) * sqrt( (f1-1.)*(f2-1.)*(f1-1.)*(f2-1.)*getEStat2(Ntt)
	                                     + (f1-1.)* f2    *(f1-1.)* f2    *getEStat2(Ntl)
	                                     +  f1    *(f2-1.)* f1    *(f2-1.)*getEStat2(Nlt)
	                                     +  f1    * f2    * f1    * f2    *getEStat2(Nll) );
}
float FakeRatios::getNpfEStat(float Ntt, float Ntl, float Nlt, float Nll, float f1, float f2, float p1, float p2){
	return p1*f2/((f1-p1)*(f2-p2)) * sqrt( (f1-1.)*(1.-p2)*(f1-1.)*(1.-p2)*getEStat2(Ntt)
	                                     + (f1-1.)*    p2 *(f1-1.)*    p2 *getEStat2(Ntl)
	                                     +  f1    *(1.-p2)* f1    *(1.-p2)*getEStat2(Nlt)
	                                     +  f1    *    p2 * f1    *    p2 *getEStat2(Nll) );
}
float FakeRatios::getNfpEStat(float Ntt, float Ntl, float Nlt, float Nll, float f1, float f2, float p1, float p2){
	return f1*p2/((f1-p1)*(f2-p2)) * sqrt( (1.-p1)*(f2-1.)*(1.-p1)*(f2-1.)*getEStat2(Ntt)
	                                     + (1.-p1)* f2    *(1.-p1)* f2    *getEStat2(Ntl)
	                                     +     p1 *(f2-1.)*    p1 *(f2-1.)*getEStat2(Nlt)
	                                     +     p1 * f2    *    p1 * f2    *getEStat2(Nll) );
}
float FakeRatios::getNffEStat(float Ntt, float Ntl, float Nlt, float Nll, float f1, float f2, float p1, float p2){
	return f1*f2/((f1-p1)*(f2-p2)) * sqrt( (1.-p1)*(1.-p2)*(1.-p1)*(1.-p2)*getEStat2(Ntt)
	                                     + (1.-p1)*    p2 *(1.-p1)*    p2 *getEStat2(Ntl)
	                                     +     p1 *(1.-p2)*    p1 *(1.-p2)*getEStat2(Nlt)
	                                     +     p1 *    p2 *    p1 *    p2 *getEStat2(Nll) );
}
float FakeRatios::getNfpNpfNffSumEStat(float Ntt, float Ntl, float Nlt, float Nll, float f1, float f2, float p1, float p2){
	return 1./((f1-p1)*(f2-p2)) * sqrt( (    p1*f2*(f1-1)*(1-p2) + f1*p2*(1-p1)*(f2-1) + f1*f2*(1-p1)*(1-p2))*(    p1*f2*(f1-1)*(1-p2) + f1*p2*(1-p1)*(f2-1) + f1*f2*(1-p1)*(1-p2))*getEStat2(Ntt)
	                                 +  (-1.*p1*f2*(f1-1)*   p2  + f1*p2*(1-p1)* f2    - f1*f2*(1-p1)*   p2 )*(-1.*p1*f2*(f1-1)*   p2  + f1*p2*(1-p1)* f2    - f1*f2*(1-p1)*   p2 )*getEStat2(Ntl)
	                                 +  (    p1*f2* f1   *(1-p2) - f1*p2*   p1 *(f2-1) - f1*f2*   p1 *(1-p2))*(    p1*f2* f1   *(1-p2) - f1*p2*   p1 *(f2-1) - f1*f2*   p1 *(1-p2))*getEStat2(Nlt)
	                                 +  (-1.*p1*f2* f1   *   p2  - f1*p2*   p1 * f2    + f1*f2*   p1 *   p2 )*(-1.*p1*f2* f1   *   p2  - f1*p2*   p1 * f2    + f1*f2*   p1 *   p2 )*getEStat2(Nll) ); 
}
float FakeRatios::getNfpNpfSumEStat(float Ntt, float Ntl, float Nlt, float Nll, float f1, float f2, float p1, float p2){
	return 1./((f1-p1)*(f2-p2)) * sqrt( (    p1*f2*(f1-1)*(1-p2) + f1*p2*(1-p1)*(f2-1))*(    p1*f2*(f1-1)*(1-p2) + f1*p2*(1-p1)*(f2-1))*getEStat2(Ntt)
	                                 +  (-1.*p1*f2*(f1-1)*   p2  + f1*p2*(1-p1)* f2   )*(-1.*p1*f2*(f1-1)*   p2  + f1*p2*(1-p1)* f2   )*getEStat2(Ntl)
	                                 +  (    p1*f2* f1   *(1-p2) - f1*p2*   p1 *(f2-1))*(    p1*f2* f1   *(1-p2) - f1*p2*   p1 *(f2-1))*getEStat2(Nlt)
	                                 +  (-1.*p1*f2* f1   *   p2  - f1*p2*   p1 * f2   )*(-1.*p1*f2* f1   *   p2  - f1*p2*   p1 * f2   )*getEStat2(Nll) ); 
}

//____________________________________________________________________________________
// Throw toy MC to determine systematics from errors on ratios
float FakeRatios::getESystFromToys2(float Ntt, float Ntl, float Nlt, float Nll, float f1, float f2, float p1, float p2, float df1, float df2, float dp1, float dp2, float(FakeRatios::*func)(float, float, float, float, float, float, float, float)){
	// Assume errors of p and f are uncorrelated
	// Throw toys in a gaussian around f and p with df and dp as their sigmas
	// Distributions for f and p are cut off at 0 and 1
	if(fVerbose > 2) cout << "FakeRatios::getESystFromToys2 ..." << endl;
	TRandom3 *rand = new TRandom3();
	rand->SetSeed();
	vector<float> f_results;
	for(size_t i = 0; i < fNToyMCs; ++i){ // vary f1
		float f1_v = rand->Gaus(f1, df1);
		if(f1_v > 1. || f1_v < 0.){--i; continue;} // throw again if f<0 or f>1
		for(size_t j = 0; j < fNToyMCs; ++j){ // vary f2
			float f2_v = rand->Gaus(f2, df2);
			if(f2_v > 1. || f2_v < 0.){--j; continue;} // throw again if f<0 or f>1
			float result = (*this.*func)(Ntt, Ntl, Nlt, Nll, f1_v, f2_v, p1, p2);
			f_results.push_back(result);
			if(fVerbose > 2) cout << result << endl;
		}
	}
	float rms_f = TMath::RMS(f_results.begin(), f_results.end());
	if(fVerbose > 2) cout << " RMS = " << rms_f << endl;
	if(fVerbose > 2) cout << endl;
	vector<float> p_results;
	for(size_t i = 0; i < fNToyMCs; ++i){ // vary p1
		float p1_v = rand->Gaus(p1, dp1);
		if(p1_v > 1. || p1_v < 0.){--i; continue;} // throw again if p<0 or p>1
		for(size_t j = 0; j < fNToyMCs; ++j){ // vary p2
			float p2_v = rand->Gaus(p2, dp2);
			if(p2_v > 1. || p2_v < 0.){--j; continue;} // throw again if p<0 or p>1
			float result = (*this.*func)(Ntt, Ntl, Nlt, Nll, f1, f2, p1_v, p2_v);
			p_results.push_back(result);
			if(fVerbose > 2) cout << result << endl;
		}
	}
	float rms_p = TMath::RMS(p_results.begin(), p_results.end());
	if(fVerbose > 2) cout << " RMS = " << rms_p << endl;
	
	return rms_f*rms_f + rms_p*rms_p;
}

//____________________________________________________________________________________
// This is the central place to fix how to calculate statistical errors
float FakeRatios::getEStat2(float N){
	if(N > 3.) return N;
	// Get these limits from TMath::ChisquaredQuantile() according to formulas
	// 32.51a and 32.51b in the PDG:
	// nulo(N) = 1/2 TMath::ChisquaredQuantile(1-0.6827, 2N)
	// nuup(N) = 1/2 TMath::ChisquaredQuantile(0.6827, 2(N+1))
	float up = 0.5 * TMath::ChisquareQuantile(0.6827, 2*(N+1));
	// float lo = 0.5 * TMath::ChisquareQuantile(1.-0.6827, 2*N);
	
	// Always just return the upper one (will always be larger)
	return ((up-N)*(up-N));
}

