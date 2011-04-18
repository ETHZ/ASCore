// C++ includes
#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>

// ROOT includes
#include <TROOT.h>
#include <TTree.h>
#include <TChain.h>

#include "MuonPlotter.hh"

using namespace std;

//_____________________________________________________________________________________
// Print out usage
void usage( int status = 0 ) {
	cout << "Usage: RunUserAnalyzer [-d dir] [-v verbose] [-c charge] [-m mode]" << endl;
	cout << "  where:" << endl;
	cout << "     dir        is the output directory               " << endl;
	cout << "                 default is TempOutput/               " << endl;
	cout << "     verbose    sets the verbose level                " << endl;
	cout << "                 default is 0 (quiet mode)            " << endl;
	// cout << "     sample   is the file with the list of samples" << endl;
	cout << "     charge     switches between SS and OS selections" << endl;
	cout << "                 default is 0 (SS), 1 is OS    " << endl;
	cout << "     mode       switches between different paths" << endl;
	cout << endl;
	exit(status);
}

//_____________________________________________________________________________________
int main(int argc, char* argv[]) {
// Default options
	TString outputdir = "MuonPlots/";
	TString samples = "samples.dat";
	int verbose = 0;
	int mode = 0;
	int charge = 0; // SS default

// Parse options
	char ch;
	while ((ch = getopt(argc, argv, "d:c:v:m:lh?")) != -1 ) {
		switch (ch) {
			case 'd': outputdir = TString(optarg); break;
			case 'c': charge = atoi(optarg); break;
			// case 's': samples = TString(optarg); break;
			case 'v': verbose = atoi(optarg); break;
			case 'm': mode = atoi(optarg); break;
			case '?':
			case 'h': usage(0); break;
			default:
			cerr << "*** Error: unknown option " << optarg << std::endl;
			usage(-1);
		}
	}

	// argc -= optind;
	// argv += optind;
	
	// Check arguments
	// if( argc<1 ) {
	// 	usage(-1);
	// }

	cout << "--------------" << endl;
	cout << "OutputDir is:      " << outputdir << endl;
	cout << "Verbose level is:  " << verbose << endl;
	cout << "Using sample file: " << samples << endl;
	cout << "Using charge sel.: " << charge << endl;
	string mode_str = "read";
	if(mode == 2) mode_str = "write";
	cout << "Mode:              " << mode_str << endl;
	cout << "--------------" << endl;

	MuonPlotter *tA = new MuonPlotter();
	tA->setOutputDir(outputdir);
	tA->setOutputFile("MuonPlots.root");
	tA->setCharge(charge);
	tA->setVerbose(verbose);
	tA->init(samples);
	if(mode == 2) tA->doLoop();
	else tA->doAnalysis();
	delete tA;
	return 0;
}

