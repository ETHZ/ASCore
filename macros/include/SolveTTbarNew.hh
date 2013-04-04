#ifndef SolveTTbarNew_h
#define SolveTTbarNew_h


//*************************************************************
//
// \class SolveTTbar: Computes a discriminant to distinguish between ttbar and stop direct production
//
//*************************************************************
#include <cmath>
#include <vector>
#include <string>
#include <iostream>
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TRandom2.h"
#include "TError.h"
#include "TH1.h"
#include "TF1.h"
#include "TF2.h"
#include "TStopwatch.h"
#include "TSystem.h"
#include "TRandom3.h"
#include "TVirtualFitter.h"

#include "TH2.h"


using namespace std;


// ************************************************************
// class SolveTTbarNew
//
// The class SolveTTbarNew contains one function attempting to find the minimum of the discriminant using root minimizers
//  -NumericalMinimization
//
// ************************************************************


class SolveTTbarNew {

private:

// private functions:
// some kinematics functions
virtual double Minv2 (double p1[], double p2[]);
virtual int InterParabol (double dz1, double dz2, double D1, double D2, double D3, double dz[]);

//put parameters into one array p[]
virtual bool GetEquationParameters(double pb1[], double pb2[], double pl1[], double pl2[], double pTm[], double p[]);

//same as equation put parameters p and free variables x in a global variable xglob
double EquationGlob(const double *xglob);

// Solves the constraint equations and returns a discriminant of type D2
virtual double Equation(double *x, double *p);



// internal constants
int debug;
double mt, mW;
double qz1lo, qz1hi, qz2lo, qz2hi;
int nscan; 
double qz1strt, qz2strt;
double discconverge;
double qz1min, qz2min;
int status;

public:

// Constructor
// constFileName = name of the file of constants (character string)
SolveTTbarNew ();

// Destructor
~SolveTTbarNew();

// Compute the discriminator *************************************
// by some in root defined minimizers compute the discriminator
double NumericalMinimization(int &statuss, double &edm, double &q1z, double&q2z, double pb1[], double pb2[], double pl1[], double pl2[], double pTm[], const char * minName = "Minuit2", const char *algoName = "" , int randomSeed = -1);

// Setting constants *************************************

// Set the debug flag
void SetDebug (int dbg) {
	debug = dbg;
}
// Set the W mass (default = 80)
void SetmW (double mWin) {
	mW = mWin;
}

// Set the top mass (default = 173)
void Setmt (double mtin) {
	mt = mtin;
}

// Set the lower and higher limts on the qz1 window for scan with ScanEqn
void SetQz1Window (double qz1l, double qz1h) {
	qz1lo = qz1l;
	qz1hi = qz1h;
}

// Set the lower and higher limts on the qz2 window for scan with ScanEqn
void SetQz2Window (double qz2l, double qz2h) {
	qz2lo = qz2l;
	qz2hi = qz2h;
}

// Set the number of bins (or steps) for scan in qz1 and qz2 with ScanEqn or IterDeriv
void SetNscan (int nsc) {
	nscan = nsc;
}

// Set the precison on the discriminant for convergence of the scan in qz1 and qz2 with ScanEqn or IterDeriv
void SetDiscrConv (double dconv) {
	discconverge = dconv;
}

// Set the starting value of qz1 for the search with IterDeriv or IterEqn
void SetStartqz1 (double qz) {
	qz1strt = qz;
}

// Set the starting value of qz2 for the search with IterDeriv or IterEqn
void SetStartqz2 (double qz) {
	qz2strt = qz;
}

// Getting more information *************************************

// getting qz1 at minimum
double GetQz1Min (void) {
	return qz1min;
}

// getting qz2 at minimum
double GetQz2Min (void) {
	return qz2min;
}

// getting the status of ending of IterDeriv
// = 0 if converged
// = 1 if ended in looping
// = 2 if ended for too many steps
int GetStatus (void) {
	return status;
}

// printing the event contents
virtual void EvtPrint (double pb1[], double pb2[], double pl1[], double pl2[], double pn1[], double pn2[]);

 




};

#endif
