#ifndef FakeRatios_hh
#define FakeRatios_hh


#include <cmath>
#include <vector>
#include <iostream>

using namespace std;

class FakeRatios {
	
public:
	FakeRatios();
	virtual ~FakeRatios();

	inline void setNToyMCs(int n){fNToyMCs = n;};
	inline void setVerbose(int n){fVerbose = n;};

	//_____________________________________________________________________________
	// Input
	// Syntax is Ntt, Ntl, (Nlt), Nll
	void setMMNtl(float, float, float);
	void setEENtl(float, float, float);
	void setEMNtl(float, float, float, float);
	// Syntax is ratio, error
	void setMFRatio(float, float);
	void setMPRatio(float, float);
	void setEFRatio(float, float);
	void setEPRatio(float, float);

	//_____________________________________________________________________________
	// Output
	// Syntax is value, estat, esyst
	float getMMNpp();
	float getMMNppEStat();
	float getMMNppESyst();
	float getMMNpf();
	float getMMNpfEStat();
	float getMMNpfESyst();
	float getMMNff();
	float getMMNffEStat();
	float getMMNffESyst();

	float getEENpp();
	float getEENppEStat();
	float getEENppESyst();
	float getEENpf();
	float getEENpfEStat();
	float getEENpfESyst();
	float getEENff();
	float getEENffEStat();
	float getEENffESyst();
	
	float getEMNpp();
	float getEMNppEStat();
	float getEMNppESyst();
	float getEMNpf();
	float getEMNpfEStat();
	float getEMNpfESyst();
	float getEMNfp();
	float getEMNfpEStat();
	float getEMNfpESyst();
	float getEMNff();
	float getEMNffEStat();
	float getEMNffESyst();
	
	float getMMTotFakes();
	float getMMTotEStat();
	float getMMTotESyst();
	float getEETotFakes();
	float getEETotEStat();
	float getEETotESyst();
	float getEMTotFakes();
	float getEMTotEStat();
	float getEMTotESyst();

	float getTotFakes();
	float getTotEStat();
	float getTotESyst();
	float getTotFakes(float, float, float, float);

	void printOutput();

	//_____________________________________________________________________________
	// Engine
	// Input syntax is: Ntt, Ntl, Nlt, Nll, f1, f2, p1, p2
	// In case of 1=2, use the 0.5*Nt1 for Ntl and Nlt and use 2*Npf
	float getNpp(float, float, float, float, float, float, float, float);
	float getNpf(float, float, float, float, float, float, float, float);
	float getNfp(float, float, float, float, float, float, float, float);
	float getNff(float, float, float, float, float, float, float, float);
	float getNfpNpfNffSum(float, float, float, float, float, float, float, float);
	float getNppEStat(float, float, float, float, float, float, float, float);
	float getNpfEStat(float, float, float, float, float, float, float, float);
	float getNfpEStat(float, float, float, float, float, float, float, float);
	float getNffEStat(float, float, float, float, float, float, float, float);
	float getNfpNpfNffSumEStat(float, float, float, float, float, float, float, float);

	float getESystFromToys2(float, float, float, float, float, float, float, float, float, float, float, float, float(FakeRatios::*)(float, float, float, float, float, float, float, float));

	// Central place to fix how to handle statistical errors
	float getEStat2(float);
	inline float getEStat(float N){return sqrt(getEStat2(N));};

private:
	int fVerbose; // default 0
	int fNToyMCs; // default 100
	
	float fMMNtl[3]; // tt, tl, ll
	float fEENtl[3]; // tt, tl, ll
	float fEMNtl[4]; // tt, tl, lt, ll
	
	float fMFRatio[2]; // ratio, error
	float fMPRatio[2]; // ratio, error
	float fEFRatio[2]; // ratio, error
	float fEPRatio[2]; // ratio, error

};

#endif
