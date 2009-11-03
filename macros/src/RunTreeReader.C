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

/**********************************************************************
>	Usage:
	single file: ./RunTreeReader 0 flag tag file1path file2path
	filelist:    ./RunTreeReader 1 flag tag filelistpath
>  Tags: Output is created in a subdir with this name,
         plots have that name in their filename
>  Flags: 5 digits correspond to 5 components:
     1. plotting all branches
     2. plotting the plotlist
     3. produce the dilepton tree
     4. produce the multiplicity plots
     5. produce the siginificance plots
   E.g. 10110 will produce all branches, the dilep tree
        and the sign plots
**********************************************************************/

int main(int argc, char* argv[]) {
	if(argc<5){ // Invalid use:
		printf("Usage: \n \t Single files  > ./RunTreeReader 0 10101 tag file1path file2path ... \n");
		printf("\t List of files > ./RunTreeReader 1 01001 tag listpath \n");
		printf("\t\t The second number flags which plots to produce. \n");
		printf("\t\t First digit is plot all branches, second is plot list \n");
		printf("\t\t third is dilepton tree, fourth is multiplicity plots \n");
		printf("\t\t and fifth is significance plots. \n");
		return -1;
	}
	TChain *theChain = new TChain("analyze/Analysis");
	int flag = atoi(argv[2]);
	TString tag = TString(argv[3]);
	cout << "Using flag " << flag << "..." << endl;
	cout << "Using tag  " << tag << " ..." << endl;
	if(atoi(argv[1]) == 0){ // case 0 (single files)
		for(int i = 4; i < argc; i++){
			theChain->Add(argv[i]);
			printf(" Adding file: %s\n",argv[i]);
		} 
	}else{ // case 1 (list of files)
		TString rootFile;
		ifstream is(argv[4]);
		while(rootFile.ReadLine(is) && (!rootFile.IsNull())){
			if(rootFile[0] == '#') continue;
			printf(" Adding file: %s\n", rootFile.Data());
			theChain->Add(rootFile);
		}
	}

	TString outputdir = "/data/wwwhome/susy/ETHPromptAnalysis/";
	TString cutfilename = "cutfile.dat";
	
	// Read parameters
	ifstream IN("treepars.dat");
	char buffer[200];
	char ParName[100], StringValue[100];
	while( IN.getline(buffer, 200, '\n') ){
		if (buffer[0] == '#') {continue;} // Skip lines commented with '#'
		sscanf(buffer, "%s %s", ParName, StringValue);
		if( !strcmp(ParName, "OutputDir") ){
			outputdir = TString(StringValue);
		}
		if( !strcmp(ParName, "CutFile") ){
			cutfilename = TString(StringValue);
		}
	}
	if(!outputdir.EndsWith("/")) outputdir += "/";

	// Print parameters
	cout << "OutputDir is:     " << outputdir << endl;
	cout << "CutFile is:       " << cutfilename << endl;
	cout << "Number of events: " << theChain->GetEntries() << endl;


	// Check which functions have to be called
	bool allbranches(false), plotlist(false), treeread(false);
	if((flag/10000)%10) allbranches = true;
	if((flag/1000)%10)  plotlist = true;
	if((flag/100%10) || (flag/10)%10 || flag%10) treeread = true;

	AnaClass *ana;
	if(allbranches || plotlist){
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
		tR->ReadCuts("cutfile.dat");
		tR->SetOutputDir(outputdir+tag);
		tR->SetTag(tag);
		tR->BeginJob();
		tR->Loop();
		tR->EndJob();
		delete tR;
	}
	return 0;
}

