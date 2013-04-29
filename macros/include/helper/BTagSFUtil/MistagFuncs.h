#ifndef MistagFuncs_h
#define MistagFuncs_h

#include "TF1.h"
#include "TLegend.h"
#include <TCanvas.h>


class MistagFuncs{

public:

  MistagFuncs();

  TF1* GetMistag(TString meanminmax, TString tagger, TString TaggerStrength, Float_t Etamin, Float_t Etamax);
  
  float GetMistag_fast(const TString& tagger, const TString& TaggerStrength, float pt, float eta );
  
  TF1* GetMistagmean(TString tagger, TString TaggerStrength, float Etamin, float Etamax);
  
  TF1* GetMistagmin(TString tagger, TString TaggerStrength, float Etamin, float Etamax);
  
  TF1* GetMistagmax(TString tagger, TString TaggerStrength, float Etamin, float Etamax);
  
  TF1* plotmean(TString tagger, TString TaggerStrength, float Etamin, float Etamax, TString opt = "" , int col = 1, float lineWidth = 1, int lineStyle = 1);
  
  TF1* plotmin(TString tagger, TString TaggerStrength, float Etamin, float Etamax, TString opt = "" , int col = 1, float lineWidth = 1, int lineStyle = 1);
  
  TF1* plotmax(TString tagger, TString TaggerStrength, float Etamin, float Etamax, TString opt = "" , int col = 1, float lineWidth = 1, int lineStyle = 1);
  
  void plotmean(TCanvas *yourC, int yourzone, TString tagger, TString TaggerStrength);
  
  TCanvas *plotmean(TString tagger, TString TaggerStrength);
  
  TCanvas *plotmean(TString selecter);
  
  TCanvas *plotmean();
  

  //TF1* GetMistag(TString meanminmax, TString tagger, TString TaggerStrength, Float_t Etamin, Float_t Etamax);

};
#endif
