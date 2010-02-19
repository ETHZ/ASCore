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
	TreeReader(TTree *tree=0, int flag = 111);
	virtual ~TreeReader();
	
private:
	
};
#endif
