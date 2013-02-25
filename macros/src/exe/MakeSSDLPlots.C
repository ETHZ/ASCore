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
	cout << "Usage: RunSSDLPlotter [-d dir] [-v verbose] [-m mode] [-c datacard] [-p configfile] [-s region] [-i input_file]" << endl;
	cout << "  where:" << endl;
	cout << "     dir        is the output directory               " << endl;
	cout << "                 default is TempOutput/               " << endl;
	cout << "     verbose    sets the verbose level                " << endl;
	cout << "                 default is 0 (quiet mode)            " << endl;
	cout << "     datacard   is the datacard to be used            " << endl;
	cout << "                 default is 'datacard.dat'            " << endl;
	cout << "     configfile is the configuration file with the    " << endl;
	cout << "                 regions and global variables.        " << endl;
	cout << "                 default is <dir>/dumperconfig.cfg    " << endl;
	cout << endl;
	cout << "     region     is the search region you want to apply to the SMS scans  " << endl;
	cout << "                 default is none                      " << endl;
	cout << "     file       is the file on which the plotter runs " << endl;
	cout << "                the scanSMS function                  " << endl;
	cout << endl;
	exit(status);
}

//_____________________________________________________________________________________
int main(int argc, char* argv[]) {
// Default options
	TString outputdir = "SSDLPlots/";
	TString datacard  = "DataCard_SSDL.dat";
	int verbose = 0;
	TString configfile;
	bool cfg = false;
	TString region = "";
	TString file = "";

// Parse options
	char ch;
	while ((ch = getopt(argc, argv, "d:c:v:m:s:i:p:lh?")) != -1 ) {
		switch (ch) {
			case 'd': outputdir = TString(optarg); break;
			case 'v': verbose = atoi(optarg); break;
			case 'c': datacard = TString(optarg); break;
			case 'p': configfile = TString(optarg); break;
			case 's': region = TString(optarg); break; 
			case 'i': file   = TString(optarg); break; 
			case '?':
			case 'h': usage(0); break;
			default:
			cerr << "*** Error: unknown option " << optarg << std::endl;
			usage(-1);
		}
	}

	if(configfile.Length() > 0) {
		cfg = true;
		cout << "Using specified config file: " << configfile << endl;
	}
	else cout << "Using default configuration file: " << outputdir << "/dumperconfig.cfg"<< endl;

	cout << "--------------" << endl;
	cout << "OutputDir is:      " << outputdir << endl;
	cout << "Verbose level is:  " << verbose << endl;
	cout << "--------------" << endl;

	SSDLPlotter *tA;
	if (cfg) tA = new SSDLPlotter(configfile);
	else tA = new SSDLPlotter(outputdir+"/dumperconfig.cfg");
	tA->setOutputDir(outputdir);
	tA->setOutputFile("SSDLHistos.root");
	tA->setVerbose(verbose);
	tA->init(datacard);
	if (region != "")	
	  tA->doSMSscans(region, file);
	else 
	  tA->doAnalysis();

	delete tA;
	return 0;
}

