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
	cout << "Usage: RunSSDLAnalyzer [-o outfile] [-s] [-v verbose] [-m maxevents] [-j JSON] [-p pthat] [-l] file1 [... filen]" << endl;
	cout << "  where:" << endl;
	cout << "     outfile   is the output file                    " << endl;
	cout << "                default is ssdlfile.root             " << endl;
	cout << "     -s        toggles between data/mc               " << endl;
	cout << "                default is data                      " << endl;
	cout << "     verbose   sets the verbose level                " << endl;
	cout << "                default is 0 (quiet mode)            " << endl;
	cout << "     maxevents are the number of events to run over  " << endl;
	cout << "               default is -1 (all)                   " << endl;
	cout << "     JSON      path of a JSON file to use            " << endl;
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
	bool isdata = true;
	// TString outputdir = "TempOutput/";
	TString outputfile = "ssdlfile.root";
	string jsonfile = "";
	int verbose = 0;
	Long64_t maxevents = -1;
	float pthatcut = -1.0;

// Parse options
	char ch;
	while ((ch = getopt(argc, argv, "o:sv:p:m:j:lh?")) != -1 ) {
		switch (ch) {
			case 'o': outputfile = TString(optarg); break;
			case 'l': isList     = true; break;
			case 's': isdata     = false; break;
			case 'v': verbose    = atoi(optarg); break;
			case 'p': pthatcut   = atof(optarg); break;
			case 'm': maxevents  = atoi(optarg); break;
			case 'j': jsonfile   = string(optarg); break;
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
	cout << "OutputFile is:    " << outputfile << endl;
	cout << "Verbose level is: " << verbose << endl;
	cout << "JSON file is:     " << (jsonfile.length()>0?jsonfile:"empty") << endl;
	cout << "Number of events: " << theChain->GetEntries() << endl;
	if(pthatcut > -1.) cout << "Lower pthat cut: " << pthatcut << endl;
	cout << "Running on " << (isdata?"data":"MC") << endl;
	cout << "--------------" << endl;

	SSDLAnalyzer *tA = new SSDLAnalyzer(theChain);
	tA->SetOutputFile(outputfile);
	tA->SetData(isdata);
	tA->SetVerbose(verbose);
	tA->SetMaxEvents(maxevents);
	if ( jsonfile.length() ) tA->ReadJSON(jsonfile.c_str());
	tA->SetPtHatCut(pthatcut);
	tA->BeginJob();
	tA->Loop();
	tA->EndJob();
	delete tA;
	return 0;	
}

