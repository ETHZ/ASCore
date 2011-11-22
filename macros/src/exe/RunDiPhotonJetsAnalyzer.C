// C++ includes
#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>

// ROOT includes
#include <TROOT.h>
#include <TTree.h>
#include <TChain.h>

#include "DiPhotonJetsAnalyzer.hh"

using namespace std;

//________________________________________________________________________________________
// Print out usage
void usage( int status = 0 ) {
	cout << "Usage: RunDiPhotonJetsAnalyzer [-f datamc] [-d dir] [-v verbose] [-p datapileup] [-P MCpileup] [-n MaxEvents] [-j jsonfile] [-x xsec(pb)] [-L nlumi(/fb)] [-l] file1 [... filen]" << endl;
	cout << "  where:" << endl;
	cout << "     dir      is the output directory               " << endl;
	cout << "               default is TempOutput/               " << endl;
	cout << "     verbose  sets the verbose level                " << endl;
	cout << "               default is 0 (quiet mode)            " << endl;
	cout << "     filen    are the input files (by default: ROOT files)" << endl;
	cout << "              with option -l, these are read as text files" << endl;
	cout << "              with one ROOT file name per line      " << endl;
	cout << endl;
	exit(status);
}

//________________________________________________________________________________________
int main(int argc, char* argv[]) {
// Default options
	bool isList = false;
	TString outputdir = "TempOutput/";
	int verbose = 0;
	
	int maxEvents=-1;
	string jsonFileName = "";
	string  data_PileUp = "";
	string  mc_PileUp = "";
	string dataType = "";
	double xsec=-1;
	double nlumi=-1;
	
	// Parse options
	char ch;
	while ((ch = getopt(argc, argv, "f:d:v:j:p:P:n:x:L:lh?")) != -1 ) {
	  switch (ch) {
	  case 'd': outputdir = TString(optarg); break;
	  case 'f': dataType = TString(optarg); break;
	  case 'v': verbose = atoi(optarg); break;
	  case 'l': isList = true; break;
	  case '?':
	  case 'h': usage(0); break;
	  case 'p': data_PileUp     = string(optarg); break;
	  case 'P': mc_PileUp       = string(optarg); break;
	  case 'n': maxEvents = atoi(optarg); break;
	  case 'j': jsonFileName = string(optarg); break;
	  case 'x': xsec = atof(optarg); break;
	  case 'L': nlumi = atof(optarg); break;
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
	cout << "Events to process: " << maxEvents << endl;
	cout << "Input cross section: " << xsec << " pb" << endl;
	cout << "Lumi to normalize to: " << nlumi << " /fb" << endl;
	cout << "JSON file is: " << (jsonFileName.length()>0?jsonFileName:"empty") << endl;
	cout << "Data/MC flag: " << dataType << endl;
	cout << "MC_PileUp file:                 " << (mc_PileUp.length()>0?mc_PileUp:"empty") << endl;
	cout << "Data_PileUp file:               " << (data_PileUp.length()>0?data_PileUp:"empty") << endl;
	cout << "--------------" << endl;

	double AddWeight;
	if (nlumi==-1 || xsec==-1) AddWeight=1;
	else {
	  double effentries;
	  effentries = maxEvents==-1 ? theChain->GetEntries() : maxEvents;
	  AddWeight=nlumi/(effentries/(1000*xsec));
	}

	if (verbose) cout << "Reweighting factor for luminosity rescaling: " << AddWeight << endl;

	DiPhotonJetsAnalyzer *tA = new DiPhotonJetsAnalyzer(theChain,dataType,AddWeight);
	tA->SetOutputDir(outputdir);
	tA->SetVerbose(verbose);
	tA->SetMaxEvents(maxEvents);
	if ( jsonFileName.length() ) tA->ReadJSON(jsonFileName.c_str());
	tA->BeginJob(data_PileUp, mc_PileUp);
	tA->Loop();
	tA->EndJob();
	delete tA;
	return 0;
}

