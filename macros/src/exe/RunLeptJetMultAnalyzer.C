// C++ includes
#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>

// ROOT includes
#include <TROOT.h>
#include <TTree.h>
#include <TChain.h>

#include "LeptJetMultAnalyzer.hh"

using namespace std;

//________________________________________________________________________________________
// Print out usage
void usage( int status = 0 ) {
	cout << "Usage: RunLeptJetMultAnalyzer [-d dir] [-o filename] [-v verbose] [-j json] [-m set_of_cuts] [-n maxEvents] [-t type] [-x lumi] [-l] file1 [... filen]" << endl;
	cout << "  where:" << endl;
	cout << "     dir           is the output directory                                   " << endl;
	cout << "                   default is TempOutput/                                    " << endl;
	cout << "     filename      is the output filename for the MassAnalysis               " << endl;
	cout << "     verbose       sets the verbose level                                    " << endl;
	cout << "                   default is 0 (quiet mode)                                 " << endl;
	cout << "     json          json file to be read                                      " << endl;
	cout << "     set_of_cuts   optional cuts for MultiplicityAnalysis                    " << endl;
	cout << "                   default is cuts from cleaning                             " << endl;
	cout << "     lumi          integrated lumi (Monte Carlo only!!)                      " << endl;
	cout << "                   scale multiplicity plots with xsection to lumi            " << endl;
	cout << "                   this affects only the MultiplicityAnalysis                " << endl;
	cout << "     type          data or mc=default                                        " << endl;
	cout << "     filen         are the input files (by default: ROOT files)              " << endl;
	cout << "                   with option -l, these are read as text files              " << endl;
	cout << "                   with one ROOT file name per line                          " << endl;
	cout << endl;
	exit(status);
}

//________________________________________________________________________________________
int main(int argc, char* argv[]) {
// Default options
	bool isList = false;
	TString outputdir = "TempOutput/";
	TString filename  = "MassTree.root";
	TString setofcuts = "default";
  	string  jsonFileName = " ";
	string type = "mc";
	bool isData = false;
	int verbose  = 0;
	int maxEvents=-1;
	float lumi   = -999.99;
	

// Parse options
	char ch;
	while ((ch = getopt(argc, argv, "d:o:v:j:m:n:t:x:lh?")) != -1 ) {
		switch (ch) {
			case 'd': outputdir       = TString(optarg); break;
			case 'o': filename        = TString(optarg); break;
			case 'v': verbose         = atoi(optarg); break;
			case 'j': jsonFileName    = string(optarg); break;
			case 'm': setofcuts       = TString(optarg); break;
			case 'n': maxEvents       = atoi(optarg); break;
			case 't': type            = string(optarg); break;
			case 'x': lumi     	  = atof(optarg); break;
			case 'l': isList          = true; break;
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
	if      (type=="data") isData =true;
	else if (type=="mc"  ) isData =false;
	else    usage(-1);

	// setofcuts="/shome/pnef/SUSY/SUSY_macros/LeptJetMult/20101129_macros/multiplicity_cuts/"+setofcuts+".dat";
	setofcuts="/shome/pnef/SUSY/SUSY_macros/LeptJetMult/multiplicity_cuts/"+setofcuts+".dat";
	
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
	cout << "OutputDir is:                   " << outputdir << endl;
	cout << "Type is:                        " << type << endl;
	cout << "Verbose level is:               " << verbose << endl;
  	cout << "JSON file is:                   " << (jsonFileName.length()>0?jsonFileName:"empty") << endl;
	if(setofcuts!="multiplicity_cuts/default"){
		cout << "Set of Cuts is:                 " << setofcuts << endl;
	}
	if(lumi > 0){
		cout << "scaling multiplicity plots with xsection for int lumi= " << lumi << endl;
	}

	cout << "Number of events:               " << theChain->GetEntries() << endl;
	cout << "--------------" << endl;

	LeptJetMultAnalyzer *tA = new LeptJetMultAnalyzer(theChain);
	tA->SetOutputDir(outputdir);
	tA->SetVerbose(verbose);
	tA->SetMaxEvents(maxEvents);
  	if (jsonFileName!=" ") tA->ReadJSON(jsonFileName.c_str());
	tA->BeginJob(filename, setofcuts, lumi, isData);
	tA->Loop();
	tA->EndJob();
	delete tA;
	return 0;
}

