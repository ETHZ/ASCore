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
	cout << "Usage: RunLeptJetMultAnalyzer [-d dir] [-o filename] [-v verbose] [-j json] " << endl;
	cout << "                              [-m set_of_cuts] [-n maxEvents] [-t type] [-x lumi] " << endl;
	cout << "                              [-p data_PileUp] [-P mc_PileUP]                " << endl; 
        cout << "                              [-s S3,noPU]                                        " << endl;
	//cout << "                              [-noPU]                                        " << endl;
	cout << "                              [-l] file1 [... filen]"                          << endl;
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
	cout << "     data_PileUp   root file from which the expected # pile-up               " << endl;
	cout << "                   interactions is read                                      " << endl;
	cout << "     mc_PileUP     root file from which the generated # pile up              " << endl;
	cout << "                   interactions is read                                      " << endl;
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
	string puScenario = "";
  	string  jsonFileName = "";
	string  data_PileUp = "";
	string  mc_PileUp = "";
	string type = "mc";
	bool isData = false;
	int verbose  = 0;
	int maxEvents=-1;
	float lumi   = -999.99;
	bool isS3 = false;
	bool noPU = false;

// Parse options
	char ch;
	while ((ch = getopt(argc, argv, "s:d:o:v:j:m:n:p:P:t:x:lh?")) != -1 ) {
	  switch (ch) {
	  case 'd': outputdir       = TString(optarg); break;
	  case 'o': filename        = TString(optarg); break;
	  case 'v': verbose         = atoi(optarg); break;
	  case 'j': jsonFileName    = string(optarg); break;
	  case 'm': setofcuts       = TString(optarg); break;
	  case 'n': maxEvents       = atoi(optarg); break;
	  case 'p': data_PileUp     = string(optarg); break;
	  case 'P': mc_PileUp       = string(optarg); break;
	  case 't': type            = string(optarg); break;
	  case 'x': lumi     	  = atof(optarg); break;
	  case 's': puScenario      = string(optarg); break;
	  case 'l': isList          = true; break;
	    //case 'noPU': noPU = true; break;  
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
	if      (!isData && data_PileUp.length()==0  ) {
		cout << "WARNING: need data_PileUp to run on MC " << endl;
	}
	if      ( isData && (data_PileUp.length() >0 || mc_PileUp.length() >0)  ) {
		cout << "ERROR: you are running on data, no reweighting needed... " << endl; exit(-1);
	}

//	setofcuts   ="/shome/leo/Analysis/cuts/"+setofcuts+".dat";
	setofcuts   ="/shome/leo/Analysis/runManagerT3_MT2Trees/cuts/"+setofcuts+".dat";
	//if(data_PileUp.length()!=0){data_PileUp ="/shome/pnef/SUSY/SUSY_macros/LeptJetMult/certification/pileUp_data/"+data_PileUp;}
	//if(mc_PileUp.length()  !=0){mc_PileUp   ="/shome/pnef/SUSY/SUSY_macros/LeptJetMult/certification/pileUp_mc/"  +mc_PileUp;}
	if(data_PileUp.length()!=0){data_PileUp ="/shome/leo/Analysis/runManagerT3_MT2Trees/"            + data_PileUp;}
	if(mc_PileUp.length()  !=0){mc_PileUp   ="/shome/leo/Analysis/runManagerT3_MT2Trees/pileup_mc/"  + mc_PileUp;}

	if(jsonFileName.length() !=0){jsonFileName="/shome/pnef/SUSY/SUSY_macros/LeptJetMult/certification/"           +jsonFileName;}

	if(puScenario=="S3") isS3=true;
	else if(puScenario=="noPU") noPU=true;

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
  	cout << "MC_PileUp file:                 " << (mc_PileUp.length()>0?mc_PileUp:"empty") << endl;
  	cout << "Data_PileUp file:               " << (data_PileUp.length()>0?data_PileUp:"empty") << endl;
	if(noPU) cout << "WARNING: NoPU option set, all the PU weights will be set to 1" << endl;
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
	tA->isS3 = isS3;
	tA->noPU = noPU;
  	if (jsonFileName!="") tA->ReadJSON(jsonFileName.c_str());
	tA->BeginJob(filename, setofcuts, lumi, isData, data_PileUp, mc_PileUp);
	tA->Loop();
	tA->EndJob();

	delete tA;
	return 0;
}

