#ifndef RunEfficiencyCalculator_hh
#define RunEfficiencyCalculator_hh

#include <TTree.h>
#include <TStyle.h>

#include "base/TreeAnalyzerBase.hh"
#include "base/TreeReader.hh"
#include "RunEfficiency.hh"
#include <string>

class EfficiencyCalculator : public TreeAnalyzerBase {
public:
  EfficiencyCalculator(TTree *tree = 0, std::string dataType="mc", bool fullCleaning=false );
  virtual ~EfficiencyCalculator();
  void BeginJob(string data_PileUp, string mc_PileUp);
  void EndJob();
  void Loop();
  void SetMaxEvents(int a){fMaxEvents=a;}
  void SetOutputFile(TString a){fOutputFile=a;}
  void SetOutputFileName(string a){outputFileName_=a;}

private:
  RunEfficiency *runEfficiency;
  int fMaxEvents;
  TString fOutputFile;
  string outputFileName_;

};
#endif
