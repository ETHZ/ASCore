///////////////////////////////////////////////////////////////////////
//
//    FILE: GoodRunList.h
//   CLASS: GoodRunList
// AUTHORS: Santiago Folgueras, Benjamin Stiegerb
//    DATE: 20/06/2011
//
// CONTENT: An utility class to check whether current Run/Lumi is or  not
//          in any certified JSON file. 
//
///////////////////////////////////////////////////////////////////////
#ifndef GOODRUNLIST_H
#define GOODRUNLIST_H 1

#include <vector>
#include <iostream>

class GoodRunList {
 public:
  GoodRunList(const char* JSONpath);
  ~GoodRunList(){};
  
  const bool IsInitialized() { 
    return (!(fRunLumis.size()==0));
  }
  void ReadJSON(const char* JSONpath);
  const bool CheckRunLumi(int Run, int LumiSection) const;
  const bool CheckRun(int RunNumber) const;
  
  void Dump(std::ostream& os) const;

 protected:
  struct RunLumi{
    int run;
    std::vector<int> lumi_min;
    std::vector<int> lumi_max;
  } fCurRunLumi; 

  std::vector<RunLumi> fRunLumis;
  
};

#endif

  

