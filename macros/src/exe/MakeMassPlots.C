// C++ includes
#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>

// ROOT includes
#include <TROOT.h>
#include <TTree.h>
#include <TChain.h>

#include "MassPlotter.hh"

using namespace std;

//_____________________________________________________________________________________
// Print out usage
void usage( int status = 0 ) {
	cout << "Usage: RunUserAnalyzer [-d dir] [-v verbose] [-s sample]" << endl;
	cout << "  where:" << endl;
	cout << "     dir      is the output directory               " << endl;
	cout << "               default is TempOutput/               " << endl;
	cout << "     verbose  sets the verbose level                " << endl;
	cout << "               default is 0 (quiet mode)            " << endl;
	cout << "     sample   is the file with the list of samples" << endl;
	cout << endl;
	exit(status);
}

//_____________________________________________________________________________________
int main(int argc, char* argv[]) {
// Default options
	TString outputdir = "MassPlots/";
	TString samples = "samples.dat";
	int verbose = 0;

// Parse options
	char ch;
	while ((ch = getopt(argc, argv, "d:v:s:lh?")) != -1 ) {
		switch (ch) {
			case 'd': outputdir = TString(optarg); break;
			case 's': samples = TString(optarg); break;
			case 'v': verbose = atoi(optarg); break;
			case '?':
			case 'h': usage(0); break;
			default:
			cerr << "*** Error: unknown option " << optarg << std::endl;
			usage(-1);
		}
	}
	
	cout << "outputdir " << outputdir << " samples " << samples << endl;	
	argc -= optind;
	argv += optind;
	
	// Check arguments
	if( argc<1 ) {
//	 	usage(-1);
	}

	cout << "--------------" << endl;
	cout << "OutputDir is:     " << outputdir << endl;
	cout << "Verbose level is: " << verbose << endl;
	cout << "Using sample file: " << samples << endl;
	cout << "--------------" << endl;

	MassPlotter *tA = new MassPlotter(outputdir, "MassPlots.root");
	tA->setVerbose(verbose);
	tA->init(samples);
	tA->makePlots();
//	tA->makeZnunu();
//	tA->makeSmallCopy(100000, 0);
	delete tA;
	return 0;
}

