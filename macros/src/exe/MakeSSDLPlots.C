// C++ includes
#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>

// ROOT includes
#include <TROOT.h>
#include <TTree.h>
#include <TChain.h>

#include "SSDLPlotter.hh"

using namespace std;

//_____________________________________________________________________________________
// Print out usage
void usage( int status = 0 ) {
	cout << "Usage: RunSSDLPlotter [-d dir] [-v verbose] [-s region] [-m mode] [-c datacard]" << endl;
	cout << "  where:" << endl;
	cout << "     dir        is the output directory               " << endl;
	cout << "                 default is TempOutput/               " << endl;
	cout << "     verbose    sets the verbose level                " << endl;
	cout << "                 default is 0 (quiet mode)            " << endl;
	cout << "     datacard   is the datacard to be used            " << endl;
	cout << "                 default is 'datacard.dat'            " << endl;
	cout << "     region     is the search region you want to apply to the SMS scans  " << endl;
	cout << "                 default is none                      " << endl;
	cout << endl;
	exit(status);
}

//_____________________________________________________________________________________
int main(int argc, char* argv[]) {
// Default options
	TString outputdir = "SSDLPlots/";
	TString datacard  = "DataCard_SSDL.dat";
	int verbose = 0;
	TString region = "";

// Parse options
	char ch;
	while ((ch = getopt(argc, argv, "d:c:v:s:m:lh?")) != -1 ) {
		switch (ch) {
			case 'd': outputdir = TString(optarg); break;
			case 'v': verbose = atoi(optarg); break;
			case 'c': datacard = TString(optarg); break;
 		        case 's': region = TString(optarg); break; 
			case '?':
			case 'h': usage(0); break;
			default:
			cerr << "*** Error: unknown option " << optarg << std::endl;
			usage(-1);
		}
	}

	cout << "--------------" << endl;
	cout << "OutputDir is:      " << outputdir << endl;
	cout << "Verbose level is:  " << verbose << endl;
	cout << "--------------" << endl;

	SSDLPlotter *tA = new SSDLPlotter();
	tA->setOutputDir(outputdir);
	tA->setVerbose(verbose);
	tA->init(datacard);
	tA->doAnalysis();
	if (region != "")	
	  tA->doSMSscans(region);

	delete tA;
	return 0;
}

