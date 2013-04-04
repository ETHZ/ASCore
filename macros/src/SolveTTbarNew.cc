// Compute a discriminant bewteen ttbar and stop direct production

#include <cmath>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>

#include "helper/Davismt2.h"
#include "SolveTTbarNew.hh"


//****************************************************************************
// Constructor
//****************************************************************************

SolveTTbarNew::SolveTTbarNew () :
debug(0),
mt(173.), mW(80.4),
qz1lo(-3000.), qz1hi(3000.), qz2lo(-3000.), qz2hi(3000.), nscan(1000), qz1strt(0.), qz2strt(0.), discconverge(0.000001),
qz1min(0.), qz2min(0.),
status(-1)
{

 return;
}

//****************************************************************************
// Destructor
//****************************************************************************

SolveTTbarNew::~SolveTTbarNew () {

}

//NEW
bool SolveTTbarNew::GetEquationParameters(double pb1[], double pb2[], double pl1[], double pl2[], double pTm[], double p[]){
//set all fixed parameters into one array

	for(int i = 0; i<=3; ++i){ p[i] = pb1[i]; p[i+4] = pb2[i]; p[i+8] = pl1[i]; p[i+12] = pl2[i]; if(i<2) p[i+16] = pTm[i]; }
	return true;
}

double SolveTTbarNew::EquationGlob(const double *xglob){
// same as equation but putting free parameters x and fixed parameters p into one array.
	double x[2]; double p[18];
	for(int i = 0; i<2; ++i)  x[i] = xglob[i];
	for(int i = 0; i<18; ++i) p[i] = xglob[i+2];
	return Equation(x,p);

}

double SolveTTbarNew::Equation(double *x, double *p){
// Solves the constraint equations and returns a discriminant of type D2
// expects qz in q[2] to be filled for q1 and q2
// also completes the arrays q1 and q2
//x are free parameters, i.e. qzs,[2 free parameters]
//p are fixed parameters, i.e. met,b1,b2,l1,l2;[18 fixed parameters]

	double pb1[4],pb2[4],pl1[4],pl2[4],pTm[2];
	//p[i] is effectively xglob[i+2]!!
	for(int i = 0; i<=3; ++i){ pb1[i] = p[i]; pb2[i] = p[i+4]; pl1[i] = p[i+8]; pl2[i] = p[i+12]; if(i<2) pTm[i] = p[i+16]; }
	double q1[4], q2[4];

	double pb1mod = pb1[3];
	double pl1mod = pl1[3];
	double pb2mod = pb2[3];
	double pl2mod = pl2[3];
    
	//double qz1 = q1[2];
	//double qz2 = q2[2];
	//NEW
	double qz1 = x[0];
	double qz2 = x[1];

	// compute the constants
	double A1 = 0.5 * (mt*mt/pb1mod - mW*mW/pl1mod);
	double A2 = 0.5 * (mt*mt/pb2mod - mW*mW/pl2mod);
	double v1[3], v2[3];
	v1[0] = pb1[0] / pb1[3] - pl1[0] / pl1[3];
	v1[1] = pb1[1] / pb1[3] - pl1[1] / pl1[3];
	v1[2] = pb1[2] / pb1[3] - pl1[2] / pl1[3];
	v2[0] = pb2[0] / pb2[3] - pl2[0] / pl2[3];
	v2[1] = pb2[1] / pb2[3] - pl2[1] / pl2[3];
	v2[2] = pb2[2] / pb2[3] - pl2[2] / pl2[3];
	double Ml1b1 = Minv2 (pl1, pb1);
	double Ml2b2 = Minv2 (pl2, pb2);
	double c1 = (Ml1b1*Ml1b1 + mW*mW) / (2.*pb1mod) - A1;
	double c2 = (Ml2b2*Ml2b2 + mW*mW) / (2.*pb2mod) - A2 - v2[0]*pTm[0] - v2[1]*pTm[1];

	
	// Now compute qx1 and qy1
	double num1 = v1[2]*qz1 - c1;
	double num2 = v2[2]*qz2 - c2;
	double denom = -v1[0]*v2[1] + v1[1]*v2[0];
	q1[0] =  (v2[1]*num1 + v1[1]*num2) / denom;
	q1[1] = -(v2[0]*num1 + v1[0]*num2) / denom;
	q1[2] = qz1;
	q2[0] = pTm[0] - q1[0];
	q2[1] = pTm[1] - q1[1];
	q2[2] = qz2;
//	cout << " q1 px, py, pz = " << q1[0] << ", " << q1[1] << ", " << q1[2] << endl;
//	cout << " q2 px, py, pz = " << q2[0] << ", " << q2[1] << ", " << q2[2] << endl;
	
	double q1mod = sqrt(q1[0]*q1[0] + q1[1]*q1[1] + q1[2]*q1[2]);
	double q2mod = sqrt(q2[0]*q2[0] + q2[1]*q2[1] + q2[2]*q2[2]);
	q1[3] = q1mod;
	q2[3] = q2mod;
		
	double mW1 = sqrt(2.*(pl1[3]*q1[3] - pl1[0]*q1[0] - pl1[1]*q1[1] - pl1[2]*q1[2]));
	double mt1 = sqrt(2.*(pb1[3]*q1[3] - pb1[0]*q1[0] - pb1[1]*q1[1] - pb1[2]*q1[2]) + Ml1b1*Ml1b1 + mW*mW);
	double mW2 = sqrt(2.*(pl2[3]*q2[3] - pl2[0]*q2[0] - pl2[1]*q2[1] - pl2[2]*q2[2]));
	double mt2 = sqrt(2.*(pb2[3]*q2[3] - pb2[0]*q2[0] - pb2[1]*q2[1] - pb2[2]*q2[2]) + Ml2b2*Ml2b2 + mW*mW);

	double D = 0.5 * sqrt( ((mW1-mW)*(mW1-mW) + (mW2-mW)*(mW2-mW)) + ((mt1-mt)*(mt1-mt) + (mt2-mt)*(mt2-mt)) );

	if (debug > 3) {
		cout << " qz1 = " << qz1 << ", qz2 = " << qz2 << endl;
		cout << " mW1 = " << mW1 << ", mW2 = " << mW2 << ", mt1 = " << mt1 << ", mt2 = " << mt2 << endl;
		cout << " D = " << D << endl;
	}

	return D;

}

// does Minimization of equation using some minimizers defined in root
// output is discriminant of type D2
double SolveTTbarNew::NumericalMinimization(int &statuss, double &edm, double &q1z, double&q2z, double pb1[], double pb2[], double pl1[], double pl2[], double pTm[], const char *minName, const char *algoName, int randomSeed)
{
   // create minimizer giving a name and a name (optionally) for the specific
   // algorithm
   // possible choices are: 
   //     minName                  algoName
   // Minuit /Minuit2             Migrad, Simplex,Combined,Scan  (default is Migrad)
   //  Minuit2                     Fumili2
   //  Fumili
   //  GSLMultiMin                ConjugateFR, ConjugatePR, BFGS, 
   //                              BFGS2, SteepestDescent
   //  GSLMultiFit
   //   GSLSimAn
   //   Genetic
   ROOT::Math::Minimizer* min = 
      ROOT::Math::Factory::CreateMinimizer(minName, algoName);

   // set tolerance , etc...
   min->SetMaxFunctionCalls(1000000); // for Minuit/Minuit2 
   min->SetMaxIterations(10000);  // for GSL 
   min->SetTolerance(0.001);
   min->SetPrintLevel(0);

   // create funciton wrapper for minmizer
   // a IMultiGenFunction type 
   SolveTTbarNew dummy;
   ROOT::Math::Functor f(&dummy, &SolveTTbarNew::EquationGlob, 20); 
   double step[2] = {0.001,0.001};
   // starting point
    
   double variable[2] = { q1z,q2z};
   //double variable[2] = { 0.,0.};
   if (randomSeed >= 0) { 
      TRandom2 r(randomSeed);
      variable[0] = r.Uniform(-20,20);
      variable[1] = r.Uniform(-20,20);
   }
 
   min->SetFunction(f);
 
   // Set the free variables to be minimized!
//   min->SetVariable(0,"x",variable[0], step[0]);
//   min->SetVariable(1,"y",variable[1], step[1]);
   min->SetLimitedVariable(0,"x",variable[0], step[0],-8050.,8050.);
   min->SetLimitedVariable(1,"y",variable[1], step[1],-8050.,8050.);

   // Set fixed variables
   min->SetFixedVariable(2,  "pb1x",    pb1[0]);
   min->SetFixedVariable(3,  "pb1y",    pb1[1]);
   min->SetFixedVariable(4,  "pb1z",    pb1[2]);
   min->SetFixedVariable(5,  "pb1mag",  pb1[3]);
   min->SetFixedVariable(6,  "pb2x",    pb2[0]);
   min->SetFixedVariable(7,  "pb2y",    pb2[1]);
   min->SetFixedVariable(8,  "pb2z",    pb2[2]);
   min->SetFixedVariable(9,  "pb2mag",  pb2[3]);
   min->SetFixedVariable(10, "pl1x",    pl1[0]);
   min->SetFixedVariable(11, "pl1y",    pl1[1]);
   min->SetFixedVariable(12, "pl1z",    pl1[2]);
   min->SetFixedVariable(13, "pl1mag",  pl1[3]);
   min->SetFixedVariable(14, "pl2x",    pl2[0]);
   min->SetFixedVariable(15, "pl2y",    pl2[1]);
   min->SetFixedVariable(16, "pl2z",    pl2[2]);
   min->SetFixedVariable(17, "pl2mag",  pl2[3]);
   min->SetFixedVariable(18, "pTmx",    pTm[0]);
   min->SetFixedVariable(19, "pTmy",    pTm[1]);
   // Set fixed but smeared Variables --> very slow ~ 5 h
/*   double st = 0.01;
   min->SetLimitedVariable(2,  "pb1x",    pb1[0], st, 0.8*pb1[0], 1.2*pb1[0]);
   min->SetLimitedVariable(3,  "pb1y",    pb1[1], st, 0.8*pb1[1], 1.2*pb1[1]);
   min->SetLimitedVariable(4,  "pb1z",    pb1[2], st, 0.8*pb1[2], 1.2*pb1[2]);
   min->SetLimitedVariable(5,  "pb1mag",  pb1[3], st, 0.8*pb1[3], 1.2*pb1[3]);
   min->SetLimitedVariable(6,  "pb2x",    pb2[0], st, 0.8*pb2[0], 1.2*pb2[0]);
   min->SetLimitedVariable(7,  "pb2y",    pb2[1], st, 0.8*pb2[1], 1.2*pb2[1]);
   min->SetLimitedVariable(8,  "pb2z",    pb2[2], st, 0.8*pb2[2], 1.2*pb2[2]);
   min->SetLimitedVariable(9,  "pb2mag",  pb2[3], st, 0.8*pb2[3], 1.2*pb2[3]);
   min->SetLimitedVariable(10, "pl1x",    pl1[0], st, 0.99*pl1[0],1.01*pl1[0]);
   min->SetLimitedVariable(11, "pl1y",    pl1[1], st, 0.99*pl1[1],1.01*pl1[1]);
   min->SetLimitedVariable(12, "pl1z",    pl1[2], st, 0.99*pl1[2],1.01*pl1[2]);
   min->SetLimitedVariable(13, "pl1mag",  pl1[3], st, 0.99*pl1[3],1.01*pl1[3]);
   min->SetLimitedVariable(14, "pl2x",    pl2[0], st, 0.99*pl2[0],1.01*pl2[0]);
   min->SetLimitedVariable(15, "pl2y",    pl2[1], st, 0.99*pl2[1],1.01*pl2[1]);
   min->SetLimitedVariable(16, "pl2z",    pl2[2], st, 0.99*pl2[2],1.01*pl2[2]);
   min->SetLimitedVariable(17, "pl2mag",  pl2[3], st, 0.99*pl2[3],1.01*pl2[3]);
   min->SetLimitedVariable(18, "pTmx",    pTm[0], st, 0.8*pTm[0], 1.2*pTm[0]);
   min->SetLimitedVariable(19, "pTmy",    pTm[1], st, 0.8*pTm[1], 1.2*pTm[1]);
*/
   // do the minimization
   min->Minimize(); 
   statuss = min->Status();
   const double *xs = min->X();
   double returnValue = 9999.99;
   returnValue = min->MinValue();
   edm = min->Edm();
   q1z = xs[0], q2z = xs[1];

   if(debug>0 /*|| status>0 || returnValue>300.*/ /*|| (fabs(q1z)==1000.||fabs(q2z)==1000.||q1z==0.||q2z==0.||edm>0.01||returnValue>800.) || algoName=="Combined" || algoName=="Fumili2" || algoName=="Migrad" || (fabs(returnValue-edm)<0.1 && edm>0.01) *//*|| returnValue==edm || (algoName=="Migrad"&&statuss==0) */){
   std::cout << "Minimum: f(" << xs[0] << ", " << xs[1] << "): " 
             << "D " <<  min->MinValue()  << std::endl;
   std::cout << "Additional Info for method " << minName << " with algorithm " << algoName << ": EDM " << min->Edm() << " NCalls " << min->NCalls() << " Nfree " << min->NFree() << " errorscale " << min->ErrorDef() << " status " << statuss << " strategy " << min->Strategy() << std::endl;
//   if(debug>1){
	min->Hesse();
	std::cout << " now Provide error: " << std::endl;
	min->ProvidesError();
   	const double *gs = min->MinGradient();
   	std::cout << " gradient at minimum: ";
	std::cout << gs;// << ",  " << gs[1];// << std::endl;
   	const double *es = min->Errors();
   	std::cout << " errors at minimum: " << es[0] << ",  " << es[1] << std::endl;
//	if(debug>2){
	   std::cout << "Print results: " << std::endl;
	   min->PrintResults();
//	}
//    }
	std::cout << std::endl;
   }

//   delete xs;
   delete min;
   return returnValue;

}

void SolveTTbarNew::EvtPrint (double pb1[], double pb2[], double pl1[], double pl2[], double pn1[], double pn2[]) {


	double pTm[2];
	pTm[0] = -(pb1[0] + pb2[0] + pl1[0] + pl2[0]);
	pTm[1] = -(pb1[1] + pb2[1] + pl1[1] + pl2[1]);

	cout << " Event print: " << endl;
	cout << " b1 px, py, pz = " << pb1[0] << ", " << pb1[1] << ", " << pb1[2] << endl;
	cout << " l1 px, py, pz = " << pl1[0] << ", " << pl1[1] << ", " << pl1[2] << endl;
	cout << " b2 px, py, pz = " << pb2[0] << ", " << pb2[1] << ", " << pb2[2] << endl;
	cout << " l2 px, py, pz = " << pl2[0] << ", " << pl2[1] << ", " << pl2[2] << endl;
	cout << " pTmiss x, y =   " << pTm[0] << ", " << pTm[1] << endl;
	cout << " Generated neutrinos: " << endl;
	cout << " n1 px, py, pz = " << pn1[0] << ", " << pn1[1] << ", " << pn1[2] << endl;
	cout << " n2 px, py, pz = " << pn2[0] << ", " << pn2[1] << ", " << pn2[2] << endl;
	
	double pTsum[] = {0., 0.};
	for (int i = 0; i < 2; ++i) {
		pTsum[i] = pb1[i] + pb2[i] + pl1[i] + pl2[i] + pn1[i] + pn2[i];
	}
	cout << " pT sum x = " << pTsum[0] << ", y = " << pTsum[1] << endl;

	return;

}
/*
void SolveTTbarNew::FCNformula{int &npar, double *gin, double &f, double *par, int flag){


}
*/

// double Minv2 (vector<double> p1, vector<double> p2) {
double SolveTTbarNew::Minv2 (double p1[], double p2[]) {
	
	double minvsq = ((p1[3]+p2[3])*(p1[3]+p2[3])  - (p1[0]+p2[0])*(p1[0]+p2[0])
					 - (p1[1]+p2[1])*(p1[1]+p2[1]) - (p1[2]+p2[2])*(p1[2]+p2[2]) );
	if (minvsq < 0.) minvsq = 0.;
	return sqrt(minvsq);
	
}


// double Minv2 (vector<double> p1, vector<double> p2) {
int SolveTTbarNew::InterParabol (double dz1, double dz2, double D1, double D2, double D3, double dz[]) {
// finds the lowest point of the parabola a(x-x0)^2 + c = D, where x is the length along the joining line
// returns 0 if no solution
	
	if (D2 > D1 && D2 > D3) return 0;
	
	double x3 = sqrt(dz1*dz1 + dz2*dz2);
	double x2 = 0.5 * x3;
	
	double a = (D1 * (x2-x3) + D2*x3 - D3*x2) / (x2*x3*(x2-x3));
	double x0 = ((D1-D2)*x3*x3 - (D1-D3)*x2*x2) / (2.*a*x2*x3*(x3-x2));
	dz[0] = dz1 * (1. - x0 / x3);
	dz[1] = dz2 * x0 / x3;
	
//	double c = D1 - a * x0*x0;
//	cout << " Interp dz1 = " << dz1 << ", dz2 = " << dz2 << ", D1 = " << D1 << ", D2 = " << D2 << ", D3 = " << D3 << endl;
//	cout << "  a = " << a << ", x0 = " << x0 << ", discr min = " << c
//		<< ", dz[0] = " << dz[0] << ", dz[1] = " << dz[1] << endl;

	return 1;

}

 
