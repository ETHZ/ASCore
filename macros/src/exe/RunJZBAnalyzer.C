// C++ includes
#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>

// ROOT includes
#include <TROOT.h>
#include <TTree.h>
#include <TChain.h>

#include "JZBAnalyzer.hh"

using namespace std;

//________________________________________________________________________________________
// Print out usage
void usage( int status = 0 ) {
  cout << "Usage: RunJZBAnalyzer [-o filename] [-v verbose] [-n maxEvents] [-t type] [-c] [-l] file1 [... filen]" << endl;
  cout << "  where:" << endl;
  cout << "     -c       runs full lepton cleaning             " << endl;
  cout << "     filename    is the output filename             " << endl;
  cout << "               default is /tmp/delete.root          " << endl;
  cout << "     verbose  sets the verbose level                " << endl;
  cout << "               default is 0 (quiet mode)            " << endl;
  cout << "     maxEvents number of events to process          " << endl;
  cout << "     type     is 'el', 'mu' or 'mc' (default)       " << endl;
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
  bool fullCleaning = false;
  //	TString outputfile = "/tmp/delete.root";
  string outputFileName = "/tmp/delete.root";
  int verbose = 0;
  int maxEvents=-1;
  string type = "data";

  // Parse options
  char ch;
  while ((ch = getopt(argc, argv, "o:v:n:t:lh?c")) != -1 ) {
    switch (ch) {
    case 'o': outputFileName = string(optarg); break;
    case 'v': verbose = atoi(optarg); break;
    case 'l': isList = true; break;
    case '?':
    case 'h': usage(0); break;
    case 'n': maxEvents = atoi(optarg); break;
    case 't': type = string(optarg); break;
    case 'c': fullCleaning = true; break;
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
  cout << "Type is: " << type << endl;
  cout << "Full cleaning is " << (fullCleaning?"ON":"OFF") << endl;
  cout << "--------------" << endl;

  JZBAnalyzer *tA = new JZBAnalyzer(theChain,type,fullCleaning);
  //	tA->SetOutputFile(outputfile);
  tA->SetOutputFileName(outputFileName);
  tA->SetVerbose(verbose);
  tA->SetMaxEvents(maxEvents);
  tA->BeginJob();
  tA->Loop();
  tA->EndJob();
  delete tA;
  return 0;
}

