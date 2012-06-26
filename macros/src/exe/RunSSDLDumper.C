// C++ includes
#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>

// ROOT includes
#include <TROOT.h>
#include <TTree.h>
#include <TChain.h>

#include "SSDLDumper.hh"

using namespace std;

//_____________________________________________________________________________________
// Print out usage
void usage( int status = 0 ) {
	cout << "Usage: RunSSDLDumper [-v verbose] [-i input]/[-l datacard] [-n name] [-m datamc] [-o output] [-c channel] [-x cross-section]" << endl;
	cout << "  where:" << endl;
	cout << "     verbose         sets the verbose level                  " << endl;
	cout << "                        default is 0 (quiet mode)            " << endl;
	cout << "     input           is the file containing an SSDLTree      " << endl;
	cout << "     name            is the short tag included in the histos " << endl;
	cout << "     datamc          switches between                        " << endl;
	cout << "                        (0) data, i.e. no pileup weights     " << endl;
	cout << "                        (1) SM MC                            " << endl;
	cout << "                        (2) Signal MC                        " << endl;
	cout << "                        (3) Rare SM MC                       " << endl;
	cout << "                        (4) Rare SM MC with no pile up       " << endl;
	cout << "     channel         Determines which type (mumu/elel/elmu)  " << endl;
	cout << "                     of events are considered from this      " << endl;
	cout << "                     sample. Used to avoid double counting   " << endl;
	cout << "                     of events between different channels.   " << endl;
	cout << "                       (-1) Ignore (use for MC)              " << endl;
	cout << "                        (0) MuMu                             " << endl;
	cout << "                        (1) ElEl                             " << endl;
	cout << "                        (2) ElMu                             " << endl;
	cout << "     cross-section   Is the cross-section for each MC sample " << endl;
	cout << "     datacard        switches to input of a datacard         " << endl;
	cout << "                     containing several input files          " << endl;
	cout << "     output          is the output directory                 " << endl;
	cout << endl;
	exit(status);
}

//_____________________________________________________________________________________
int main(int argc, char* argv[]) {
// Default options
	TString datacard = "";
	TString inputfile = "";
	TString outputdir = "";
	TString name = "";
	int channel = -1; // ignore(-1), mumu(0), elel(1), elmu(2)
	int verbose = 0;
	int datamc = 0;
	double xsec = 1.;
	
	bool card = false; // toggle between running on single file or datacard

	// Parse options
	char ch;
	while ((ch = getopt(argc, argv, "v:l:i:n:m:o:c:g:x:h?")) != -1 ) {
		switch (ch) {
			case 'v': verbose    = atoi(optarg);         break;
			case 'm': datamc     = atoi(optarg);         break;
			case 'c': channel    = atoi(optarg);         break;
			case 'x': xsec       = strtod(optarg, NULL); break;
			case 'l': datacard   = TString(optarg);      break;
			case 'i': inputfile  = TString(optarg);      break;
			case 'n': name       = TString(optarg);      break;
			case 'o': outputdir  = TString(optarg);      break;
			case '?':
			case 'h': usage(0); break;
			default:
			cerr << "*** Error: unknown option " << optarg << std::endl;
			usage(-1);
		}
	}

	// Check arguments
	if( argc<1 ) {
		usage(-1);
	}

	TString type = "";
	switch(datamc){
		case 0: type  = "Data";                   break;
		case 1: type  = "SM MC";                  break;
		case 2: type  = "Signal MC";              break;
		case 3: type  = "Rare SM MC";             break;
		case 4: type  = "Rare SM MC (no pileup)"; break;
		default: type = "Unknown";
	}

	TString chan = "";
	switch(channel){
		case -1: chan = "Ignored";   break;
		case  0: chan = "MuMu";      break;
		case  1: chan = "ElEl";      break;
		case  2: chan = "ElMu";      break;
		default: chan = "Ignored";
	}
	
	if(datacard.Length() > 0) card = true;

	if(verbose > 0)          cout << "------------------------------------" << endl;
	if(verbose > 0)          cout << " Verbose level is:  " << verbose << endl;
	if(verbose > 0 && !card) cout << " Inputfile is:      " << inputfile << endl;
	if(verbose > 0 && !card) cout << " Name is:           " << name << endl;
	if(verbose > 0 && !card) cout << " Type is:           " << type << endl;
	if(verbose > 0 && !card) cout << " Channel is:        " << chan << endl;
	if(verbose > 0 && !card) cout << " Cross-section is:  " << xsec << endl;
	if(verbose > 0 &&  card) cout << " Datacard is:       " << datacard << endl;
	if(verbose > 0)          cout << " Outputdir is:      " << outputdir << endl;

	SSDLDumper *tA = new SSDLDumper();
	tA->setVerbose(verbose);
	tA->setOutputDir(outputdir);
	if(!card) tA->init(inputfile, name, datamc, xsec, channel);
	if( card) tA->init(datacard);
	tA->loop();
	delete tA;
	return 0;
}

