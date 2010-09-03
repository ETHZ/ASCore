// C++ includes
#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>

// ROOT includes
#include <TROOT.h>
#include <TTree.h>
#include <TChain.h>

#include "PhysQCAnalyzer.hh"

using namespace std;

//________________________________________________________________________________________
// Print out usage
void usage( int status = 0 ) {
	cout << "Usage: RunPhysQCAnalyzer [-d dir] [-v verbose] [-n max] [-p prescale] [-l] file1 [... filen]" << endl;
	cout << "  where:" << endl;
	cout << "     dir      is the output directory                       " << endl;
	cout << "               default is TempOutput/                       " << endl;
	cout << "     verbose  sets the verbose level                        " << endl;
	cout << "               default is 0 (quiet mode)                    " << endl;
        cout << "     max      is the maximum number of events to process    " << endl;
        cout << "     prescale is the (inverse) fraction of events to process" << endl;
        cout << "              (e.g., '-p 10' will process every 10 events   " << endl;
	cout << "     filen    are the input files (by default: ROOT files)  " << endl;
	cout << "              with option -l, these are read as text files  " << endl;
	cout << "              with one ROOT file name per line              " << endl;
	cout << endl;
	exit(status);
}

//________________________________________________________________________________________
int main(int argc, char* argv[]) {
// Default options
	bool isList = false;
	TString outputdir = "TempOutput/";
	int verbose = 0;
        Long64_t maxEvents = -1;
        Int_t prescale = 0;

// Parse options
	char ch;
	while ((ch = getopt(argc, argv, "d:v:n:p:lh?")) != -1 ) {
		switch (ch) {
			case 'd': outputdir = TString(optarg); break;
			case 'v': verbose = atoi(optarg); break;
                        case 'n': maxEvents = atoi(optarg); break;
                        case 'p': prescale  = atoi(optarg); break;
			case 'l': isList = true; break;
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

        Long_t nentries = theChain->GetEntries();
        if ( nentries<1 ) {
          cerr << "ERROR: Couldn't find any entry in files." << endl;
          cerr << "       Check that files exist and URL is correct." << endl;
          cerr << "       e.g.: dcap://t3se01.psi.ch:22125/pnfs/psi.ch/cms/trivcat/store/user/foo/foo.root" << endl;
          exit(-1);
        }

	cout << "--------------" << endl;
	cout << "OutputDir is:     " << outputdir << endl;
	cout << "Verbose level is: " << verbose << endl;
	cout << "Number of events: " << nentries << endl;
	cout << "--------------" << endl;

	PhysQCAnalyzer *tA = new PhysQCAnalyzer(theChain);
	tA->SetOutputDir(outputdir);
	tA->SetVerbose(verbose);
	tA->BeginJob();
	tA->Loop(maxEvents,prescale);
	tA->EndJob();
	delete tA;
	return 0;
}

