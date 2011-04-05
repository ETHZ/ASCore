//
// Script to remove duplicates based on run, lumi and event.
// Result is a merged file, with duplicates removed.
//
// To compile: g++ -o removeDuplicates removeDuplicates.C `root-config --cflags --libs`
//
// Usage: removeDuplicates -o output <file to clean> <reference file>
//    example: removeDuplicates -o All_nodup.root EG.root Mu.root
//             this will merge EG.root and Mu.root into All_nodup.root, where
//             events from EG.root that are also in Mu.root are removed.
//
#include <iostream>
#include <iomanip>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <TTree.h>
#include <TFile.h>
#include <TChain.h>

bool duplicate( TTree* ref, Int_t& run, Int_t& lumi, Int_t& event ) {
    
    // Check if this event already exists in the tree
    char cut[256]; sprintf(cut,"runNum==%d&&lumi==%d&&eventNum==%d",run,lumi,event);
    return (ref->GetEntries(cut)>0);
    
}

//________________________________________________________________________________________
// Print out usage
void usage( int status = 0 ) {
  std::cout << "Usage: removeDuplicates [-o filename] [-v verbose] [file to clean] [reference file]" << std::endl;
  std::cout << "  where:" << std::endl;
  std::cout << "     filename is the merged output filename         " << std::endl;
  std::cout << "     verbose  sets the verbose level                " << std::endl;
  std::cout << "               default is 0 (quiet mode)            " << std::endl;
  std::cout << "     file to clean   is the file to clean up        " << std::endl;
  std::cout << "     reference file  is the file to check against   " << std::endl;
  std::cout << std::endl;
  exit(status);
}



//int removeDuplicates() {
int main(int argc, char** argv) {

    TString outputFileName;
    int verbose = 0;
    
    // Parse options
    char ch;
    while ((ch = getopt(argc, argv, "o:v:lh?c")) != -1 ) {
      switch (ch) {
      case 'o': outputFileName = std::string(optarg); break;
      case 'v': verbose = atoi(optarg); break;
      case '?':
      case 'h': usage(0); break;
      default:
        std::cerr << "*** Error: unknown option " << optarg << std::endl;
        usage(-1);
      }
    }

    argc -= optind;
    argv += optind;

    // Check arguments
    if( argc<2 || outputFileName.Length()==0 ) { usage(-1); }

    TString inFile(argv[0]);
    TString refFile(argv[1]);

    TFile* f1 = new TFile(inFile);
    TFile* f2 = new TFile(refFile);
    
    TTree* t1 = (TTree*)f1->Get("events");
    TTree* t2 = (TTree*)f2->Get("events");
    
    // Just copy the structure
    TFile* newFile = new TFile("tmp.root","RECREATE");
    TTree* newTree = t1->CloneTree(0);
    
    Int_t egRun, egLumi, egEvt;
    t1->SetBranchAddress("eventNum",&egEvt);
    t1->SetBranchAddress("runNum",  &egRun);
    t1->SetBranchAddress("lumi",    &egLumi);
    
    Int_t nentries = (Int_t)t1->GetEntries();
    Int_t duplicates = 0;
    int freq = nentries/100;
    std::cout << "Processing ...   0% " << std::flush;
    for (Int_t i=0;i<nentries; i++) {
        if ( freq>0 && !(i%freq) ) { // Counter
              std::cout << "\b\b\b\b\b" << std::setprecision(0) << std::setw(3) << std::fixed 
                        << i/static_cast<double>(nentries)*100 << "% " << std::flush;
              std::cout << std::setprecision(4);
        }
        t1->GetEntry(i);
        if ( !duplicate(t2, egRun, egLumi, egEvt) ) newTree->Fill();
        else 
            ++duplicates;
    }
    std::cout << std::endl;
    newTree->AutoSave();

    std::cout << "Found " << duplicates << " duplicate events in " << nentries;
    std::cout << " (" << duplicates/static_cast<float>(nentries)*100 << "%)" << std::endl;
    
    delete f1;
    delete f2;
    delete newFile;

    TString cmd("hadd "+outputFileName+" tmp.root "+refFile);

    std::cout << "Merging files: " << cmd << std::endl;
    gSystem->Exec(cmd);

    return 0;
    
}
