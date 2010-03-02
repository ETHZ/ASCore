#include "base/TreeReader.hh"
using namespace std;

TreeReader::TreeReader(TTree *tree) : TreeClassBase(tree){
	if( tree == 0 ) cout << "TreeReader ==> No tree!" << endl;
}

TreeReader::~TreeReader(){
	if(!fChain) cout << "TreeReader ==> No chain!" << endl;
}
