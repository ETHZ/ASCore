#ifndef TreeReader_hh
#define TreeReader_hh


#include <vector>
#include <cmath>
#include <string>
#include <iostream>
#include <fstream>

#include <TCanvas.h>
#include <TTree.h>
#include <TStyle.h>
#include <TLatex.h>
#include <TLeaf.h>
#include <TBranch.h>

#include "TreeClassBase.h"

class TreeReader : public TreeClassBase{
public:
	TreeReader(TTree *tree = NULL);
	virtual ~TreeReader();
	inline virtual Long64_t GetEntries(){return fChain->GetEntries();};
  inline virtual bool isChain() { return fIsChain; }
private:
  bool fIsChain;
	
};
#endif
