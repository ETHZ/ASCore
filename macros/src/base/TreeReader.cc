#include "base/TreeReader.hh"
using namespace std;

TreeReader::TreeReader(TTree *tree, int flag){
	if( tree == 0 ) cout << "TreeReader ==> No tree!" << endl;
	
	// Set all branch addresses, fChain = tree:
	Init(tree);
}

TreeReader::~TreeReader(){
	if(!fChain) cout << "TreeReader ==> No chain!" << endl;
}
