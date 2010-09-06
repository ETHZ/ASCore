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
	cout << "Usage: RunLeptJetMultAnalyzer [-d dir] [-v verbose] [-m set_of_cuts] [-x lumi] [-t HLT] [-l] file1 [... filen]" << endl;
	cout << "  where:" << endl;
	cout << "     dir           is the output directory                                   " << endl;
	cout << "                   default is TempOutput/                                    " << endl;
	cout << "     verbose       sets the verbose level                                    " << endl;
	cout << "                   default is 0 (quiet mode)                                 " << endl;
	cout << "     set_of_cuts   optional cuts for MultiplicityAnalysis                    " << endl;
	cout << "                   default is cuts from cleaning                             " << endl;
	cout << "     lumi          integrated lumi (Monte Carlo only!!)                      " << endl;
	cout << "                   scale multiplicity plots with xsection to lumi            " << endl;
	cout << "     HLT           set to 1 if HLT bits required/vetoed (default 0)          " << endl;
	cout << "                   HLT triggers specified in HLTfile.dat                     " << endl;
	cout << "                   this affects only the MultiplicityAnalysis                " << endl;
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
	TString setofcuts = "default";
	int HLT      = 0;
	int verbose  = 0;
	float lumi   = -999.99;
	

// Parse options
	char ch;
	while ((ch = getopt(argc, argv, "d:v:m:x:t:lh?")) != -1 ) {
		switch (ch) {
			case 'd': outputdir       = TString(optarg); break;
			case 'v': verbose         = atoi(optarg); break;
			case 'm': setofcuts       = TString(optarg); break;
			case 'x': lumi     	      = atof(optarg); break; 	
			case 't': HLT             = atoi(optarg); break; 				
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
	TString  filename="multiplicity.root";
	setofcuts="multiplicity_cuts/"+setofcuts+".dat";
	
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
	
	
	vector<std::string> requiredHLT;
	vector<std::string> vetoedHLT;
	if(HLT == 1){
		ifstream IN("HLTfile.dat");
		char buffer[200];
		char StringValue[100];
		char ParName[100];	
		bool ok(false);
	
		while( IN.getline(buffer, 200, '\n') ){
			ok = false;
			if (buffer[0] == '#') {continue;} // Skip lines commented with '#'
			sscanf(buffer, "%s %s", ParName, StringValue);
			if( !strcmp(ParName, "required") ){
				requiredHLT.push_back(StringValue); ok = true;
			}else if ( !strcmp(ParName, "vetoed") ){
				vetoedHLT.push_back(StringValue); ok = true;
			}
			if(!ok) cout << "%% RunLeptJetMultAnalyzer ==> ERROR: " << ParName << " should be either required or vetoed." << endl;
		}
	}
	
	cout << "--------------" << endl;
	cout << "OutputDir is:                   " << outputdir << endl;
	cout << "Verbose level is:               " << verbose << endl;
	if(setofcuts!="multiplicity_cuts/default"){
		cout << "Set of Cuts is:                 " << setofcuts << endl;
	}
	if(lumi > 0){
		cout << "scaling multiplicity plots with xsection for int lumi= " << lumi << endl;
	}
	if(requiredHLT.size()!=0){
		cout << "require HLT trigger (logic OR): ";
		for(int i=0; i<(requiredHLT.size()); ++i){
			cout << requiredHLT[i] << "  ";
		}
		cout << endl;
	}
	if(vetoedHLT.size()!=0){
		cout << "veto HLT trigger (logic OR):    ";
		for(int i=0; i<vetoedHLT.size(); ++i){
			cout << vetoedHLT[i] << " ";
		}
		cout << endl;
	}
	cout << "Number of events:               " << theChain->GetEntries() << endl;
	cout << "--------------" << endl;

	LeptJetMultAnalyzer *tA = new LeptJetMultAnalyzer(theChain);
	tA->SetOutputDir(outputdir);
	tA->SetVerbose(verbose);
	tA->BeginJob(filename, setofcuts, lumi, &requiredHLT, &vetoedHLT);
	tA->Loop();
	tA->EndJob();
	delete tA;
	return 0;
}

