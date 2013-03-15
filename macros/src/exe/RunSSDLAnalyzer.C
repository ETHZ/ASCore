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

#include "FWCore/FWLite/interface/AutoLibraryLoader.h"

#include "SSDLAnalyzer.hh"

using namespace std;

//________________________________________________________________________________________
// Print out usage
void usage( int status = 0 ) {
	cout << "Usage: RunSSDLAnalyzer [-o outfile] [-s] [-v verbose] [-m maxevents] [-e filleff] [-j JSON] [-p pthat] [-l] [-x data PU file] [-y MC PU file] [-g GlobalTag] file1 [... filen]" << endl;
	cout << "  where:" << endl;
	cout << "     outfile   is the output file                    " << endl;
	cout << "                default is ssdlfile.root             " << endl;
	cout << "     -s        toggles between data/mc               " << endl;
	cout << "                default is data                      " << endl;
	cout << "     verbose   sets the verbose level                " << endl;
	cout << "                default is 0 (quiet mode)            " << endl;
	cout << "     maxevents are the number of events to run over  " << endl;
	cout << "               default is -1 (all)                   " << endl;
	cout << "     filleff   switches on the efficiency tree       " << endl;
	cout << "               default is 0 (all)                   " << endl;
	cout << "     JSON      path of a JSON file to use            " << endl;
	cout << "     pthat     sets the upper cut on PtHat           " << endl;
	cout << "                default is -1.0 (no cut)                "    << endl ;
	cout << "     filen     are the input files (by default: ROOT files)" << endl ;
	cout << "               with option -l, these are read as text files" << endl ;
	cout << "               with one ROOT file name per line      "       << endl ;
	cout << "     -x        path of PU distribution file of data"         << endl ;
	cout << "     -y        path of PU distribution file of MC"           << endl ;
	cout << "     -g        Global Tag to do the jet corrections"         << endl ;
	cout << "               default is \"\", which uses the corrections"  << endl ;
	cout << "               in the original ntuple."                      << endl ;
	cout << "               Make sure that the correction files are there"<< endl ;
	cout << endl;
	exit(status);
}

//________________________________________________________________________________________
int main(int argc, char* argv[]) {
	AutoLibraryLoader::enable();
// Default options
	bool isList = false;
	bool isdata = true;
	// TString outputdir = "TempOutput/";
	TString outputfile = "ssdlfile.root";
	std::string datapufile;
	std::string mcpufile;
	std::string globaltag = "";
	string jsonfile = "";
	int verbose = 0;
	Long64_t maxevents = -1;
	bool doeff = 0;
	float pthatcut = -1.0;

// Parse options
	char ch;
	while ((ch = getopt(argc, argv, "o:sv:p:m:j:e:x:y:g:lh?")) != -1 ) {
		switch (ch) {
			case 'o': outputfile = TString(optarg); break;
			case 'l': isList     = true; break;
			case 's': isdata     = false; break;
			case 'v': verbose    = atoi(optarg); break;
			case 'p': pthatcut   = atof(optarg); break;
			case 'm': maxevents  = atoi(optarg); break;
			case 'j': jsonfile   = string(optarg); break;
			case 'e': doeff      = atoi(optarg); break;
			case 'x': datapufile = string(optarg); break;
			case 'y': mcpufile   = string(optarg); break;
			case 'g': globaltag  = string(optarg); break;
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

	// MARC TChain *theChain = new TChain("analyze/Analysis");
	// MARC for(int i = 0; i < argc; i++){
	// MARC 	if( !isList ){
	// MARC 		theChain->Add(argv[i]);
	// MARC 		printf(" Adding file: %s\n",argv[i]);
	// MARC 	} else {
	// MARC 		TString rootFile;
	// MARC 		ifstream is(argv[i]);
	// MARC 		while(rootFile.ReadLine(is) && (!rootFile.IsNull())){
	// MARC 			if(rootFile[0] == '#') continue;
	// MARC 			theChain->Add(rootFile);
	// MARC 			printf(" Adding file: %s\n", rootFile.Data());
	// MARC 		}
	// MARC 	}
	// MARC }
	std::vector<std::string> fileList;
		for(int i = 0; i < argc; i++){
			if( !isList ){
				fileList.push_back(argv[i]);
				printf(" Adding file: %s\n",argv[i]);
			} else {
				TString rootFile;
				ifstream is(argv[i]);
				while(rootFile.ReadLine(is) && (!rootFile.IsNull())){
					if(rootFile[0] == '#') continue;
					fileList.push_back(rootFile.Data());
					printf(" Adding file: %s\n", rootFile.Data());
				}
			}
		}

	cout << "--------------" << endl;
	cout << "OutputFile is:    " << outputfile << endl;
	cout << "Verbose level is: " << verbose << endl;
	cout << "JSON file is:     " << (jsonfile.length()>0?jsonfile:"empty") << endl;
	cout << "PU distribution file of data is: " << datapufile << endl;
	cout << "PU distribution file of MC is: " << mcpufile << endl;
	if (globaltag != "") cout << "GlobalTag for corrections is: " << globaltag << endl;
	cout << "Number of files: " << fileList.size() << endl;
	// MARC cout << "Number of events: " << theChain->GetEntries() << endl;
	if(pthatcut > -1.) cout << "Lower pthat cut: " << pthatcut << endl;
	cout << "Running on " << (isdata?"data":"MC") << endl;
	cout << "--------------" << endl;

	// test mar 14 SSDLAnalyzer *tA = new SSDLAnalyzer(fileList);
	SSDLAnalyzer *tA = new SSDLAnalyzer(fileList, isdata, globaltag);
	// MARC SSDLAnalyzer *tA = new SSDLAnalyzer(theChain);
	tA->SetOutputFile(outputfile);
	//tA->SetPuFiles(datapufile, mcpufile);
	tA->SetData(isdata);
	tA->SetVerbose(verbose);
	tA->SetMaxEvents(maxevents);
	if ( jsonfile.length() ) tA->ReadJSON(jsonfile.c_str());
	tA->SetPtHatCut(pthatcut);
	tA->DoFillEffTree(doeff);
	if (isdata) tA->BeginJob();
	else tA->BeginJob(datapufile, mcpufile);
	tA->Loop();
	tA->EndJob();
	delete tA;
	return 0;	
}

