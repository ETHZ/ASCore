// C++ includes
#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>

// ROOT includes
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>

#include "SSDLAnalyzer.hh"

using namespace std;

//________________________________________________________________________________________
// Print out usage
void usage( int status = 0 ) {
	cout << "Usage: RunSSDLAnalyzer [-d dir] [-v verbose] [-m maxevents] [-p pthat] [-l] file1 [... filen]" << endl;
	cout << "  where:" << endl;
	cout << "     dir       is the output directory               " << endl;
	cout << "                default is TempOutput/               " << endl;
	cout << "     verbose   sets the verbose level                " << endl;
	cout << "                default is 0 (quiet mode)            " << endl;
	cout << "     maxevents are the number of events to run over  " << endl;
	cout << "               default is -1 (all)                   " << endl;
	cout << "     pthat     sets the upper cut on PtHat           " << endl;
	cout << "                default is -1.0 (no cut)                " << endl;
	cout << "     filen     are the input files (by default: ROOT files)" << endl;
	cout << "               with option -l, these are read as text files" << endl;
	cout << "               with one ROOT file name per line      " << endl;
	cout << endl;
	exit(status);
}

//________________________________________________________________________________________
int main(int argc, char* argv[]) {
// Default options
	bool isList = false;
	TString outputdir = "TempOutput/";
	int verbose = 0;
	Long64_t maxevents = -1;
	float pthatcut = -1.0;

// Parse options
	char ch;
	while ((ch = getopt(argc, argv, "d:v:p:m:lh?")) != -1 ) {
		switch (ch) {
			case 'd': outputdir = TString(optarg); break;
			case 'l': isList = true; break;
			case 'v': verbose = atoi(optarg); break;
			case 'p': pthatcut = atof(optarg); break;
			case 'm': maxevents = atoi(optarg); break;
			case '?':
			case 'h': usage(0); break;
			default:
				cerr << "*** Error: unknown option " << optarg << std::endl;
				usage(-1);
		}
	}
	
	argc -= optind;
	argv += optind;

// Check arguments
	if( argc<1 ) {
		usage(-1);
	}

	TChain *theChain = new TChain("analyze/Analysis");
	for(int i = 0; i < argc; i++){
		if( !isList ){
			theChain->Add(argv[i]);
			printf(" Adding file: %s\n",argv[i]);
		} else {
			TString rootFile;
			ifstream is(argv[i]);
			while(rootFile.ReadLine(is) && (!rootFile.IsNull())){
				if(rootFile[0] == '#') continue;
				theChain->Add(rootFile);
				printf(" Adding file: %s\n", rootFile.Data());
			}
		}
	}

	cout << "--------------" << endl;
	cout << "OutputDir is:     " << outputdir << endl;
	cout << "Verbose level is: " << verbose << endl;
	cout << "Number of events: " << theChain->GetEntries() << endl;
	if(pthatcut > -1.) cout << "Lower pthat cut: " << pthatcut << endl;
	cout << "--------------" << endl;

	SSDLAnalyzer *tA = new SSDLAnalyzer(theChain);
	tA->SetOutputDir(outputdir);
	tA->SetVerbose(verbose);
	tA->SetMaxEvents(maxevents);
	tA->SetPtHatCut(pthatcut);
	tA->BeginJob();
	tA->Loop();
	tA->EndJob();
	delete tA;
	return 0;	
}

