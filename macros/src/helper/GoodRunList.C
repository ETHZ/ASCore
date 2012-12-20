///////////////////////////////////////////////////////////////////////
//
//    FILE: GoodRunList.C
//   CLASS: GoodRunList
// AUTHORS: Santiago Folgueras, Benjamin Stiegerb
//    DATE: 23/03/2011
//
// CONTENT: An utility class to check whether the current Run/lumi
//          is on a JSON file
//            
//
//////////////////////////////////////////////////////////////////////
#include "helper/GoodRunList.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>

using namespace std;

//#define DEBUG 1

GoodRunList::GoodRunList(const char* JSONpath){
  ReadJSON(JSONpath);
}

void GoodRunList::ReadJSON(const char* JSONpath){
  string line;
  ifstream JSON(JSONpath);
  if (JSON.is_open()){
    getline (JSON,line);
  }
  JSON.close();

  line.erase(0,1);
  line.erase(line.size()-1);
  
  // split JSON into runs with associated lumis
  string c =", \"";
  vector<string> splitted;
  while(line.find(c) != string::npos){
    int loc = line.find(c);
    string tmp=line.substr(0, (loc));
    splitted.push_back(tmp);
    line=line.erase(0, tmp.size()+2);
  }       
  splitted.push_back(line); // for the last entry
  
  // loop over runs (with lumis)
  for(unsigned int i=0; i<splitted.size(); ++i){
    fCurRunLumi.run = 0;
    fCurRunLumi.lumi_min.clear();
    fCurRunLumi.lumi_max.clear();
    
    // first seperate runs and save them
    if(splitted[i].find(":") != string::npos){
      int loc = splitted[i].find(":");
      string tmp = splitted[i].substr(1,loc-2);
      stringstream ss ( tmp);
      ss >> fCurRunLumi.run;
      splitted[i].erase(0, loc+1);
    }else {cout << "ERROR in ReadJSON: seperate runs" << endl; return;}

    // now we're left with the lumis, remove the surrounding "[]"
    if(splitted[i].find("[[") == string::npos || splitted[i].find("]]") == string::npos){ cout << "ERROR in ReadJSON: lumi parsing" << endl; return;}
    else{
      int loc  = splitted[i].find("[[");
      int loc2 = splitted[i].find("]]");
      splitted[i].erase(loc   , 1);
      splitted[i].erase(loc2  , 1);
    }
    // left with lumis a la [a, b], [c, d], ...
    while(splitted[i].find("[")!= string::npos){
      int loc1 = splitted[i].find("[");
      int loc2 = splitted[i].find("]");
      int loc3 = splitted[i].find(",", loc1);
      string lumi1= splitted[i].substr(loc1+1, loc3-1);
      string lumi2= splitted[i].substr(loc3+1, loc2-1);
      stringstream ss (lumi1); int lumimin;
      ss >> lumimin;
      fCurRunLumi.lumi_min.push_back(lumimin);
      stringstream ss2 (lumi2); int lumimax;
      ss2 >> lumimax;
      fCurRunLumi.lumi_max.push_back(lumimax);
      splitted[i].erase(loc1, loc2-loc1+3);
    }
    // push_back run and lumis
    fRunLumis.push_back(fCurRunLumi);
  }
  
#ifdef DEBUG
  cout << "-----------------------------------------------" << endl;
  cout << "Reading the following runs and lumi sections from JSON file:" << endl; 
  for(unsigned int i=0; i<fRunLumis.size(); ++i){
    cout << "  Run: " << fRunLumis[i].run << " Lumis: " ;
    for(unsigned int k=0; k<fRunLumis[i].lumi_min.size(); ++k){
      cout << " [" << fRunLumis[i].lumi_min[k] << ", " << fRunLumis[i].lumi_max[k] << "] ";
    }
    cout << endl;
  }
  cout << "-----------------------------------------------" << endl;
#endif 
}

const bool GoodRunList::CheckRunLumi(int RunNumber, int LumiSection) const {
  bool good(false);
  if(fRunLumis.size()==0) {
#ifdef DEBUG
    cout << "CheckRunLumi: fRunLumis is empty. Assuming no JSON." << endl; 
#endif    
    return true;
  }
  for(unsigned int r=0; r<fRunLumis.size(); ++r){
    if(RunNumber != fRunLumis[r].run) continue;
    for(unsigned int l=0; l<fRunLumis[r].lumi_min.size(); ++l){
      if(LumiSection < fRunLumis[r].lumi_min[l]) continue;
      if(LumiSection > fRunLumis[r].lumi_max[l]) continue;
      good = true;
    }
  }
#ifdef DEBUG
  if(good) cout <<  "JSON: accepted ";
  else     cout <<  "JSON: rejected ";
  cout <<  RunNumber << ":" << LumiSection << endl;
#endif
  
  return good;
}

//________________________________________________________________________________________
const bool GoodRunList::CheckRun(int RunNumber) const {
  // Check if current run is in JSON file
  bool runInJson(false);
  
  // Not initialized: assume no JSON
  if(fRunLumis.size()==0) {
#ifdef DEBUG
    cout << "CheckRunLumi: fRunLumis is empty. Assuming no JSON." << endl; 
#endif
    return true;
  }
  for(unsigned int r=0; r<fRunLumis.size(); ++r){
    if(RunNumber == fRunLumis[r].run)
      runInJson = true;
  }

#ifdef DEBUG
  if (!runInJson)
    cout << "Run " <<  RunNumber << " not in JSON file" << endl;
#endif
  
  return runInJson;
}
void GoodRunList::Dump(ostream& os) const {
  os << "{";
  
  for (unsigned int i = 0; i < fRunLumis.size(); i++) {
    os << "\"" << fRunLumis[i].run << "\": [";
    
    for (unsigned int j = 0; j < fRunLumis[i].lumi_min.size(); j++) {
      os << "[" << fRunLumis[i].lumi_min[j]
	 << ", " << fRunLumis[i].lumi_max[j] << "]";
      if (j != fRunLumis[i].lumi_min.size()-1)
	os << ", ";
    }
    os << "]";
    if (i != fRunLumis.size()-1)
      os << ", ";
  }
  
  os << "}";
}
// const bool IsInitialized(){
//   return (!(fRunLumis.size()==0));
// };
  
