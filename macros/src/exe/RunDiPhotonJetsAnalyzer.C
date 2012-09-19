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
	cout << "Usage: RunDiPhotonJetsAnalyzer [-o outfile] [-f datamc] [-d dir] [-v verbose] [-p datapileup] [-P MCpileup] [-n MaxEvents] [-j jsonfile] [-x xsec(pb)] [-L nlumi(/fb)] [-N events_in_dset] [-G gg k factor] [-g gj k factor] [-J jj k factor] [-Y minthrpfphotoncandEB] [-y minthrpfphotoncandEE] [-l] file1 [... filen]" << endl;
	cout << "  where:" << endl;
	cout << "     dir      is the output directory               " << endl;
	cout << "               default is current directory               " << endl;
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
	TString outputdir = "./";
	TString outputfile = "outfile.root";
	int verbose = 0;
	
	int maxEvents=-1;
	string jsonFileName = "";
	string  data_PileUp = "";
	string  mc_PileUp = "";
	string dataType = "";
	Float_t xsec=-1;
	Float_t nlumi=-1;
	Int_t nevtsindset=-1;
	Float_t minthrpfphotoncandEB = 0;
	Float_t minthrpfphotoncandEE = 0;

	Float_t kfactors[3]={1,1,1};

	// Parse options
	char ch;
	while ((ch = getopt(argc, argv, "N:G:g:J:o:f:d:v:j:p:P:n:x:L:Y:y:lh?")) != -1 ) {
	  switch (ch) {
	  case 'N': nevtsindset = atoi(optarg); break;
	  case 'G': kfactors[0] = atof(optarg); break;
	  case 'g': kfactors[1] = atof(optarg); break;
	  case 'J': kfactors[2] = atof(optarg); break;
	  case 'd': outputdir = TString(optarg); break;
	  case 'o': outputfile = TString(optarg); break;
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
	  case 'Y': minthrpfphotoncandEB = atof(optarg); break;
	  case 'y': minthrpfphotoncandEE = atof(optarg); break;
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
	cout << "OutputFile is:    " << outputfile << endl;
	cout << "Verbose level is: " << verbose << endl;
	cout << "Number of events: " << theChain->GetEntries() << endl;
	cout << "Events to process: " << maxEvents << endl;
	cout << "Input cross section: " << xsec << " pb" << endl;
	cout << "Lumi to normalize to: " << nlumi << " /fb" << endl;
	cout << "Number of events in dset for normalization: " << nevtsindset << endl;
	cout << "JSON file is: " << (jsonFileName.length()>0?jsonFileName:"empty") << endl;
	cout << "Data/MC flag: " << dataType << endl;
	cout << "MC_PileUp file:                 " << (mc_PileUp.length()>0?mc_PileUp:"empty") << endl;
	cout << "Data_PileUp file:               " << (data_PileUp.length()>0?data_PileUp:"empty") << endl;
	cout << "gg,gj,jj k-factors:  " << kfactors[0] << " " << kfactors[1] << " " << kfactors[2] << endl;
	cout << "minthrpfphotoncandEB = " << minthrpfphotoncandEB << endl;
	cout << "minthrpfphotoncandEE = " << minthrpfphotoncandEE << endl;
	cout << "--------------" << endl;

	Float_t AddWeight;
	if (nlumi==-1 || xsec==-1 || nevtsindset==-1) AddWeight=1;
	else AddWeight=nlumi*xsec*1000.0/nevtsindset;
	

	if (verbose) cout << "Reweighting factor for luminosity rescaling: " << AddWeight << endl;

	DiPhotonJetsAnalyzer *tA = new DiPhotonJetsAnalyzer(theChain,dataType,AddWeight,kfactors,minthrpfphotoncandEB,minthrpfphotoncandEE);
	tA->SetOutputDir(outputdir);
	tA->SetOutputFile(outputfile);
	tA->SetVerbose(verbose);
	tA->SetMaxEvents(maxEvents);
	if ( jsonFileName.length() ) tA->ReadJSON(jsonFileName.c_str());
	tA->BeginJob(data_PileUp, mc_PileUp);
	tA->Loop();
	tA->EndJob();
	delete tA;
	return 0;
}

