#include "helper/Utilities.hh"
#include "RunEfficiency.hh"
#include "TF1.h"
#include <time.h>
#include <TRandom.h>
#include "TF1.h"


using namespace std;


#define jMax 30  // do not touch this
#define metMax 30
#define rMax 30
#define nAnalysis 3



//********************************* FUNCTIONS THAT ARE STILL MISSING !!!!

const bool RunEfficiency::ElPassingIsoTag(int a, int b) {
cout << "ElPassingIsoTag NOT DEFINED!" << endl;
}

const bool RunEfficiency::ElPassingIsoProbe(int a, int b) {
cout << "ElPassingIsoProbe NOT DEFINED!" << endl;
}

const bool RunEfficiency::ElPassingIDTag(int a, int b) {
cout << "ElPassingIDTag NOT DEFINED!" << endl;
}

const bool RunEfficiency::ElPassingIDPProbe(int a,int b,int c) {
cout << "ElPassingIDPProbe NOT DEFINED!" << endl;
}

const bool RunEfficiency::ElPassingRecoTag(int a, int b) {
cout << "ElPassingRecoTag NOT DEFINED!" << endl;
}
const bool RunEfficiency::ElPassingRecoPProbe(int a, int b) {
cout << "ElPassingRecoPProbe NOT DEFINED!" << endl;
}
const bool RunEfficiency::ElPassingIsoPProbe(int a, int b) {
cout << "ElPassingIsoPProbe NOT DEFINED!" << endl;
}