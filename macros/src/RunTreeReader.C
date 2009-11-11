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

#include "TreeReader.hh"
#include "AnaClass.hh"

using namespace std;

//________________________________________________________________________________________
// Print out usage
void usage( int status = 0 ) {
	cout << "Usage: RunTreeReader -f flag -t tag [-d dir] [-o objcuts] [-e evtcuts] [-l] file1 [... filen]" << endl;
	cout << "  where:" << endl;
	cout << "     flag     is a 5-digits bitmap that controls    " << endl;
	cout << "              the following components:             " << endl;
	cout << "               1. plotting all branches             " << endl;
	cout << "               2. plotting the plotlist             " << endl;
	cout << "               3. producing the dilepton tree       " << endl;
	cout << "               4. plotting the multiplicity plots   " << endl;
	cout << "               5. plotting the significance plots   " << endl;
	cout << "              E.g., 10110 will produce all branches," << endl;
	cout << "              the dilep. tree and the sig. plots    " << endl;
	cout << "     tag      is the name of the output directory,  " << endl;
	cout << "              and will appear in the file names     " << endl;
	cout << "     dir      is the output directory               " << endl;
	cout << "     objcuts  is the object selection filename      " << endl;
	cout << "     evtcuts  is the event selection filename       " << endl;
	cout << "     filen    are the input files (by default: ROOT files)" << endl;
	cout << "              with option -l, these are read as text files" << endl;
	cout << "              with one ROOT file name per line      " << endl;
	cout << endl;
	exit(status);
}

//________________________________________________________________________________________
int main(int argc, char* argv[]) {
// Default options
	unsigned int flag=0;
	TString tag;
	bool isList = false;
	TString outputdir = "/data/wwwhome/susy/ETHPromptAnalysis/";
	TString objselfile = "objsel.dat";
	TString evtselfile = "evtsel.dat";

// Parse options
	char ch;
	while ((ch = getopt(argc, argv, "d:c:l:f:t:o:e:h?")) != -1 ) {
		switch (ch) {
			case 'd': outputdir = TString(optarg); break;
			case 'l': isList = true; break;
			case 'f': flag = atoi(optarg); break;
			case 't': tag = TString(optarg); break;
			case 'o': objselfile = TString(optarg); break;
			case 'e': evtselfile = TString(optarg); break;
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
	if ( !flag || !(tag.Length()>0) || argc<1 ) {
		usage(-1);
	}

	TChain *theChain = new TChain("analyze/Analysis");
	for(int i = 0; i < argc; i++){
		if ( !isList ) {
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
	cout << "Using flag        " << flag << "..." << endl;
	cout << "Using tag         " << tag << "..." << endl;
	cout << "OutputDir is:     " << outputdir << endl;
	cout << "ObjSelFile is:    " << objselfile << endl;
	cout << "EvtSelFile is:    " << evtselfile << endl;
	cout << "Number of events: " << theChain->GetEntries() << endl;
	cout << "--------------" << endl;

// Check which functions have to be called
	bool allbranches(false), plotlist(false), treeread(false);
	if((flag/10000)%10) allbranches = true;
	if((flag/1000)%10)  plotlist = true;
	if((flag/100%10) || (flag/10)%10 || flag%10) treeread = true;
	cout << "Will " << (allbranches?"":"not ") << "plot all branches" << endl;
	cout << "Will " << (plotlist?"":"not ")    << "plot plotlist"     << endl;
	cout << "Will " << (treeread?"":"not ")    << "read tree"    << endl;
	cout << "--------------" << endl;

	AnaClass *ana;
	if (allbranches || plotlist) {
		ana = new AnaClass();
		ana->readVarNames("varnames.dat");
		ana->setOutputDir(outputdir);
		ana->setGlobalTag(tag);
		if(allbranches) ana->plotAllBranches(theChain, tag);		
		if(plotlist)    ana->plotPlotList("plotlist.dat", theChain, tag);
		delete ana;
	}

	TreeReader *tR;
	if(treeread){
		tR = new TreeReader(theChain, flag%1000);
		tR->ReadObjCuts(objselfile);
		tR->ReadEvtSel(evtselfile);
		tR->SetOutputDir(outputdir+tag);
		tR->SetTag(tag);
		tR->BeginJob();
		tR->Loop();
		tR->EndJob();
		delete tR;
	}
	return 0;
}

