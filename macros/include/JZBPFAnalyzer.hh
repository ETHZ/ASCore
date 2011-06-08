#ifndef JZBAnalyzer_hh
#define JZBAnalyzer_hh

#include <TTree.h>
#include <TStyle.h>

#include "base/TreeAnalyzerBase.hh"
#include "base/TreeReader.hh"
#include "JZBPFAnalysis.hh"
#include <string>

class JZBPFAnalyzer : public TreeAnalyzerBase {
public:
  JZBPFAnalyzer(TTree *tree = 0, std::string dataType="mc", bool fullCleaning=false );
  virtual ~JZBPFAnalyzer();
  void BeginJob(string data_PileUp, string mc_PileUp);
  void EndJob();
  void Loop();
  void SetMaxEvents(int a){fMaxEvents=a;}
  void SetOutputFile(TString a){fOutputFile=a;}
  void SetOutputFileName(string a){outputFileName_=a;}

private:
  JZBPFAnalysis *fJZBPFAnalysis;
  int fMaxEvents;
  TString fOutputFile;
  string outputFileName_;

};
#endif
