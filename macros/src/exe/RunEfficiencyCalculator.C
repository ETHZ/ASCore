// C++ includes
#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>

// ROOT includes
#include <TROOT.h>
#include <TTree.h>
#include <TChain.h>

#include "EfficiencyCalculator.hh"

using namespace std;

//________________________________________________________________________________________
// Print out usage
void usage( int status = 0 ) {
  cout << "Usage: RunEfficiencyCalculator [-o filename] [-v verbose] [-n maxEvents] [-j JSON] [-t type] [-c] [-l] [-p data_PileUp] [-P mc_PileUP] file1 [... filen]" << endl;
  cout << "  where:" << endl;
  cout << "     -c       runs full lepton cleaning                                       " << endl;
  cout << "     filename    is the output filename                                       " << endl;
  cout << "               default is /tmp/delete.root                                    " << endl;
  cout << "     verbose  sets the verbose level                                          " << endl;
  cout << "               default is 0 (quiet mode)                                      " << endl;
  cout << "     maxEvents number of events to process                                    " << endl;
  cout << "     JSON     path of a JSON file to use                                      " << endl;
  cout << "     data_PileUp   root file from which the expected # pile-up                " << endl;
  cout << "                   interactions is read                                       " << endl;
  cout << "     mc_PileUP     root file from which the generated # pile up               " << endl;
  cout << "     type     is 'el', 'mu' or 'mc' (default)                                 " << endl;
  cout << "     filen    are the input files (by default: ROOT files)                    " << endl;
  cout << "              with option -l, these are read as text files                    " << endl;
  cout << "              with one ROOT file name per line                                " << endl;
  cout << endl;
  exit(status);
}

//________________________________________________________________________________________
int main(int argc, char* argv[]) {
  // Default options
  bool isList = false;
  bool fullCleaning = false;
  //	TString outputfile = "/tmp/delete.root";
  string outputFileName = "/tmp/delete.root";
  string jsonFileName = "";
  string  data_PileUp = "";
  string  mc_PileUp = "";

  int verbose = 0;
  int maxEvents=-1;
  string type = "data";

  // Parse options
  char ch;
  while ((ch = getopt(argc, argv, "o:v:n:j:t:lh?cp:P")) != -1 ) {
    switch (ch) {
    case 'o': outputFileName = string(optarg); break;
    case 'v': verbose = atoi(optarg); break;
    case 'l': isList = true; break;
    case '?':
    case 'h': usage(0); break;
    case 'n': maxEvents = atoi(optarg); break;
    case 'j': jsonFileName = string(optarg); break;
    case 't': type = string(optarg); break;
    case 'c': fullCleaning = true; break;
    case 'p': data_PileUp     = string(optarg); break;
    case 'P': mc_PileUp       = string(optarg); break;

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
  cout << "outputFileName is:     " << outputFileName << endl;
  cout << "Verbose level is: " << verbose << endl;
  cout << "Number of events: " << theChain->GetEntries() << endl;
  cout << "Events to process: " << maxEvents << endl;
  cout << "JSON file is: " << (jsonFileName.length()>0?jsonFileName:"empty") << endl;
  cout << "Type is: " << type << endl;
  cout << "Full cleaning is " << (fullCleaning?"ON":"OFF") << endl;
  cout << "MC_PileUp file:                 " << (mc_PileUp.length()>0?mc_PileUp:"empty") << endl;
  cout << "Data_PileUp file:               " << (data_PileUp.length()>0?data_PileUp:"empty") << endl;

  cout << "--------------" << endl;

  EfficiencyCalculator *tA = new EfficiencyCalculator(theChain,type,fullCleaning);
  //	tA->SetOutputFile(outputfile);
  tA->SetOutputFileName(outputFileName);
  tA->SetVerbose(verbose);
  tA->SetMaxEvents(maxEvents);
  if ( jsonFileName.length() ) tA->ReadJSON(jsonFileName.c_str());
  tA->BeginJob(data_PileUp, mc_PileUp);
  tA->Loop();
  tA->EndJob();
  delete tA;
  return 0;
}

