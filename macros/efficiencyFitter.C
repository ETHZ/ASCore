/*______________________________________________________________________
________________________________________________________________________
__________________  Efficiency Fitter  _________________________________
________________________________________________________________________
______________________________________________________________________*/


//ROOT Files
#include "TH1F.h"
#include "TF1.h"
#include "TTree.h"
#include "TFile.h"
#include "TCut.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TGraphErrors.h"
#include "TGraph2DErrors.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TH1D.h"
#include "TRandom.h"

//RooFit
#ifndef __CINT__
#include "/swshare/cms/slc5_amd64_gcc434/lcg/roofit/5.28.00a-cms5/include/RooGlobalFunc.h"
#endif
#include "/swshare/cms/slc5_amd64_gcc434/lcg/roofit/5.28.00a-cms5/include/RooRealVar.h"
#include "/swshare/cms/slc5_amd64_gcc434/lcg/roofit/5.28.00a-cms5/include/RooDataSet.h"
#include "/swshare/cms/slc5_amd64_gcc434/lcg/roofit/5.28.00a-cms5/include/RooDataHist.h"
#include "/swshare/cms/slc5_amd64_gcc434/lcg/roofit/5.28.00a-cms5/include/RooGaussian.h"
#include "/swshare/cms/slc5_amd64_gcc434/lcg/roofit/5.28.00a-cms5/include/RooPlot.h"
#include "/swshare/cms/slc5_amd64_gcc434/lcg/roofit/5.28.00a-cms5/include/RooExtendPdf.h"
#include "/swshare/cms/slc5_amd64_gcc434/lcg/roofit/5.28.00a-cms5/include/RooGenericPdf.h"
#include "/swshare/cms/slc5_amd64_gcc434/lcg/roofit/5.28.00a-cms5/include/RooBreitWigner.h"
#include "/swshare/cms/slc5_amd64_gcc434/lcg/roofit/5.28.00a-cms5/include/RooChebychev.h"
#include "/swshare/cms/slc5_amd64_gcc434/lcg/roofit/5.28.00a-cms5/include/RooAddPdf.h"
#include "/swshare/cms/slc5_amd64_gcc434/lcg/roofit/5.28.00a-cms5/include/RooSimultaneous.h"
#include "/swshare/cms/slc5_amd64_gcc434/lcg/roofit/5.28.00a-cms5/include/RooFitResult.h"
#include "/swshare/cms/slc5_amd64_gcc434/lcg/roofit/5.28.00a-cms5/include/RooCategory.h"
#include "/swshare/cms/slc5_amd64_gcc434/lcg/roofit/5.28.00a-cms5/include/RooVoigtian.h"
#include "/swshare/cms/slc5_amd64_gcc434/lcg/roofit/5.28.00a-cms5/include/RooExponential.h"

//C++
#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>

using namespace RooFit;

float kefficiency,  ke_efficiency1, ke_efficiency2;
float kpt, keta, kdr, knjets;
float kptwidth, ketawidth, kdrwidth, knjetswidth;


//________________________EfficiencyFitter class________________________________________
class EfficiencyFitter {

public:
  
  EfficiencyFitter(std::string, std::string, std::string);
  
  ~EfficiencyFitter();
  
  void doPlotAgainstVar(std::string, std::vector<float>, std::string, std::string, std::string);

  void doPlotAgainstVarMC(std::string, std::vector<float>, std::string, std::string, std::string);

  void doTriggerPlotAgainstVar(std::string, std::vector<float>, std::string, std::string, std::string);

  void doPlotAgainstVarAndVar(std::string, std::string, std::vector<float>, std::vector<float>, std::string, std::string, std::string, std::string);
  
  void doPlotAgainstVarAndVarMC(std::string, std::string, std::vector<float>, std::vector<float>, std::string, std::string, std::string, std::string);



private:

  TGraphErrors * produceGraph(std::vector<float>, std::vector<float>, std::vector<float>, std::string, std::string, int);

  void calculateTriggerEfficiency(TCut, float &, float &, std::string);

  TGraph2DErrors * produce2DGraph(std::vector<std::vector<float> >, std::vector<std::vector<float> >, std::vector<std::vector<float> >, std::vector<std::vector<float> >, std::string, std::string, std::string, int );

  void calculateMCNEfficiency(TCut, float &, float &);

  void calculateFitEfficiency(std::string, float &, float &, std::string);

  TCut convertIntoCut(std::string, float, float);

  void calculateNEfficiency(TCut, float &, float &);

  //data_______________________________________________
  TFile *f, *foutput;
  TTree *myTree, *s;
  template<class T> std::string any2string(T i);
  TCut tpp, tp, standard, tppMC, ntpp;

};


//__________________________________________________________________________
EfficiencyFitter::EfficiencyFitter(std::string fileName, std::string outputFileName, std::string theCut) {
  
  
  f = new TFile(fileName.c_str());
  myTree = (TTree *) f->Get("events");

  ntpp = "((ch1>0 &&tag1 == 1 && pprobe2 == 0) || (ch2>0 && tag2 == 1 && pprobe1 == 0))";
  tpp = "((ch1>0 &&tag1==1 &&pprobe2==1) || (ch2>0 && tag2==1 && pprobe1 == 1))";
  tp = "((ch1>0 &&tag1==1 &&probe2==1) || (ch2>0 && tag2==1 && probe1==1))";
  tppMC = "((ch1>0 && pprobe2==1) || (ch2>0 && pprobe1==1))";
  
  
  standard = theCut.c_str();

  foutput = new TFile(outputFileName.c_str(), "RECREATE");
  foutput->cd();

}


//__________________________________________________________________________
EfficiencyFitter::~EfficiencyFitter() {

  foutput->Write();
  foutput->Close();
  f->Close();

};


//__________________________________________________________________________
template<class T>
std::string EfficiencyFitter::any2string(T i)
{
  std::ostringstream buffer;
  buffer << i;
  return buffer.str();
}




//__________________________________________________________________________
TCut EfficiencyFitter::convertIntoCut(std::string var, float minv, float maxv) {

  std::string a = "(" + var + " > " + any2string(minv) + " && " +  var + " < " + any2string(maxv) + ")";
  TCut *theCut = new TCut(a.c_str());
  return *theCut;

}


//___________________________________________________________________________
void EfficiencyFitter::calculateNEfficiency(TCut acceptance, float &efficiency, float &efficiency_error) {

  float nTagProbe = myTree->Draw("mll", acceptance+standard+tp, "N");
  float nTagPProbe = myTree->Draw("mll", acceptance+standard+tpp, "N");
  
  efficiency = nTagPProbe/nTagProbe;
  efficiency_error = sqrt(efficiency*(1-efficiency)/nTagProbe);

}



//___________________________________________________________________________
void EfficiencyFitter::calculateMCNEfficiency(TCut acceptance, float &efficiency, float &efficiency_error) {


  cout << acceptance << endl;
  float nTagProbe = myTree->Draw("mll", acceptance+standard+"genMll>60&&genMll<120", "N");
  float nTagPProbe = myTree->Draw("mll", acceptance+standard+tppMC+"genMll>60&&genMll<120&&pt1!=0", "N");
  
  efficiency = nTagPProbe/nTagProbe;
  efficiency_error = sqrt(efficiency*(1-efficiency)/nTagProbe);

}




//____________________________________________________________________________
void EfficiencyFitter::calculateTriggerEfficiency(TCut acceptance, float &efficiency, float &efficiency_error, std::string type) {

  TCut trig = tpp;
  if(type == "em") {
    trig += "passingem_trigger";
  } else if(type == "ee") {
    trig += "passingee_trigger";
  } else {
    trig += "passingmm_trigger";
  }

  float nTagProbe = myTree->Draw("mll", acceptance+standard+tpp, "N");
  float nTagPProbe = myTree->Draw("mll", acceptance+standard+trig, "N");

  efficiency = nTagPProbe/nTagProbe;
  efficiency_error = sqrt(efficiency*(1-efficiency)/nTagProbe);

}



//___________________________________________________________________________
void EfficiencyFitter::calculateFitEfficiency(std::string acceptance, float &efficiency, float &efficiency_error, std::string id) {

  int nbins = 50;
  float nmin = 60;
  float nmax = 120;

  std::string nameCanvas = id + "_canvas";
  
  std::string theAcceptance = "&& " + acceptance;

  // Create observables
  RooRealVar mll("mll","mll", nmin, nmax);
  RooRealVar tag1("tag1","tag1", -1, 2);
  RooRealVar tag2("tag2","tag2", -1, 2);
  RooRealVar pprobe1("pprobe1","pprobe1", -1, 2);
  RooRealVar pprobe2("pprobe2","pprobe2", -1, 2);
  RooRealVar probe1("probe1","probe1", -1, 2);
  RooRealVar probe2("probe2","probe2", -1, 2);
  RooRealVar ch1("ch1","ch1", -1, 2);
  RooRealVar ch2("ch2","ch2", -1, 2);
  RooRealVar ptn("ptn","ptn", 0, 1000);
  RooRealVar etan("etan","etan", -3, 3);
  RooRealVar drn("drn","drn", 0, 1009);
  RooRealVar pfJetGoodNum("pfJetGoodNum","pfJetGoodNum", 0, 20);
  

  
  RooArgSet mySet(mll, tag1, tag2, pprobe1, pprobe2, probe1, probe2, ch1, ch2);
  mySet.add(ptn); mySet.add(etan); mySet.add(drn); mySet.add(pfJetGoodNum);
  RooArgList myList(mll, tag1, tag2, pprobe1, pprobe2, probe1, probe2, ch1, ch2);
  myList.add(ptn); myList.add(etan); myList.add(drn); myList.add(pfJetGoodNum);
 
  std::string formpp = "mll>60&&mll<120&&((ch1>0 &&tag1==1 &&pprobe2==1) || (ch2>0 && tag2==1 && pprobe1==1))"+theAcceptance;
  std::string formnp = "mll>60&&mll<120&&((ch1>0 &&tag1==1 &&(pprobe2==0 && probe2==1)) || (ch2>0 && tag2==1 && (pprobe1==0 && probe1==1)))"+theAcceptance;
  std::string formall = "mll>60&&mll<120&&((ch1>0 &&tag1==1 &&( probe2==1)) || (ch2>0 && tag2==1 && (probe1==1)))"+theAcceptance;
 
  
  RooFormulaVar formulapp("cutpp", formpp.c_str(), myList);
  
  RooFormulaVar formulanp("cutnp", formnp.c_str(), myList);
  RooFormulaVar formulaall("cutnp", formall.c_str(), myList);

  std::string idpp = id + "_pproodataset";
  std::string idnp = id + "_nproodataset";
  std::string idall = id + "_allroodataset";
  RooDataSet passingpp(idpp.c_str(), "pp", myTree, mySet, formulapp);
  //passingpp.write("datos.txt");
  //return;

  RooDataSet passingnp(idnp.c_str(), "np", myTree, mySet, formulanp);
  RooDataSet passingall(idall.c_str(), "all", myTree, mySet, formulaall);

  RooRealVar meansignalpp("meansignalpp", "meansignalpp", 91, 89, 93);
  RooRealVar sigmasignalpp("sigmasignalpp", "sigmasignalpp", 3, 0, 20);
  RooRealVar widthsignalpp("widthsignalpp", "widthsignalpp", 2.4, 0, 20);
  RooRealVar var1backpp("var1backpp", "var1backpp", -0.02138, -2.0, 2.0);
  RooRealVar var2backpp("var2backpp", "var2backpp", 5, -1000000, 10000000);
  
  RooRealVar meansignalnp("meansignalnp", "meansignalnp", 91, 89, 93);
  RooRealVar sigmasignalnp("sigmasignalnp", "sigmasignalnp", 3, 1, 20);
  RooRealVar widthsignalnp("widthsignalnp", "widthsignalnp", 2.4, 1, 20);
  RooRealVar var1backnp("var1backnp", "var1backnp", -0.0521, -2, 2);
  RooRealVar var2backnp("var2backnp", "var2backnp", 5, 0, 1000);
  RooRealVar normpp("normpp", "normpp", 10, 0.0, 1000000);
  RooRealVar normnp("normnp", "normnp", 10, 0.0, 1000000);
  RooRealVar a("a", "a", 10, 0.0, 1000000);
  RooRealVar fpp("fpp", "fpp", 0.5, 0.0, 1);
  RooRealVar fnp("fnp", "fnp", 0.5, 0.0, 1);
  RooRealVar eff("eff", "eff", 0.99, 0.4, 1);


  //Equations
  RooVoigtian signalpp("signalpp", "signalpp", mll, meansignalpp, widthsignalpp, sigmasignalpp);
  RooVoigtian signalnp("signalnp", "signalnp", mll, meansignalnp, widthsignalnp, sigmasignalnp);
  RooExponential backpp("backpp","backpp", mll , var1backpp);
  //RooChebychev backpp("backpp","backpp", mll , RooArgSet(var1backpp,var2backpp));
  RooExponential backnp("backnp","backnp", mll , var1backnp);


  //RooGenericPdf funcpp_("funcpp_", "FunctionPP", "eff*signalpp+backpp", RooArgList(mll, eff, signalpp, backpp));
  //RooGenericPdf funcpp_("funcpp_", "FunctionPP", "eff*signalpp", RooArgList(mll, eff, signalpp));
  //RooGenericPdf funcnp_("funcnp_", "FunctionNP", "(1-eff)*signalpp", RooArgList(mll, eff, signalpp));
  RooGenericPdf funcpp_("funcpp_", "FunctionPP", "eff*signalpp", RooArgList(mll, eff, signalpp));
  RooGenericPdf funcnp_("funcnp_", "FunctionNP", "a*backnp", RooArgList(mll, eff, a, backnp));

  RooAddPdf funcpp("funcpp", "funcpp", RooArgList(signalpp, backpp), RooArgList(fpp));
  //RooAddPdf funcnp("funcnp", "funcnp", RooArgList(backnp, signalnp), RooArgList(fpp));
  RooAddPdf funcnp("funcnp", "funcnp", RooArgList(backnp, signalnp), RooArgList(fpp));
  
  //RooExtendPdf funcpp("funcpp", "funcpp", funcpp_, normpp);
  //RooExtendPdf funcnp("funcnp", "funcnp", funcnp_, normnp);

  RooCategory sample("pagging","passing") ;
  sample.defineType("passing") ;
  sample.defineType("nopassing") ;

  RooDataSet combData("combData","combined data",mll,Index(sample),Import("passing",passingpp),Import("nopassing",passingnp));

  RooSimultaneous final("simPdf","simultaneous pdf",sample) ;

  final.addPdf(funcpp,"passing") ;
  final.addPdf(funcnp,"nopassing") ;


  RooFitResult *result = final.fitTo(combData, Save()); 

  RooRealVar* e = (RooRealVar*) result->floatParsFinal().find("fpp");
  efficiency = e->getVal();
  efficiency_error = e->getError();

  result->Print();
  
  //RooPlot* frame1 = mll.frame(Bins(30),Title("Physics sample")) ;
  //passingall.plotOn(frame1);
  //final.plotOn(frame1) ;
 
  RooPlot* frame2 = mll.frame(Bins(30),Title("Physics sample")) ;
  passingpp.plotOn(frame2);
  funcpp.plotOn(frame2) ;
  
  RooPlot* frame3 = mll.frame(Bins(30),Title("Physics sample")) ;
  passingnp.plotOn(frame3);
  funcnp.plotOn(frame3) ;

  TCanvas* ca = new TCanvas(nameCanvas.c_str(), "Fit",800,400) ;
  ca->Divide(2,1);
  //ca->cd(1) ; gPad->SetLeftMargin(0.15) ; frame1->GetYaxis()->SetTitleOffset(1.4) ; frame1->Draw();
  ca->cd(1) ; gPad->SetLeftMargin(0.15) ; frame2->GetYaxis()->SetTitleOffset(1.4) ; frame2->Draw();
  ca->cd(2) ; gPad->SetLeftMargin(0.15) ; frame3->GetYaxis()->SetTitleOffset(1.4) ; frame3->Draw();
  ca->Write();    

  
  delete e;
  //delete frame1;
  delete frame2;
  delete frame3;
  //delete resultpp;
  //delete resultnp;
  //delete result;
 

 
}



//________________________________________________________
TGraphErrors * EfficiencyFitter::produceGraph(std::vector<float> x, std::vector<float> y, std::vector<float> e, std::string xlabel, std::string ylabel, int col) {

  float *theX = new float [x.size()];
  float *theY = new float [y.size()];
  float *theE = new float [e.size()];

  for(int c = 0; c < x.size(); ++c) {

    theX[c] = x.at(c);
    theY[c] = y.at(c);
    theE[c] = e.at(c);

  }


  TGraphErrors *theGraph = new TGraphErrors(x.size(), theX, theY, NULL, theE);
  theGraph->GetXaxis()->SetTitle(xlabel.c_str());
  theGraph->GetYaxis()->SetTitle(ylabel.c_str());
  theGraph->SetMarkerColor(col);
  theGraph->SetTitle("");
  theGraph->GetYaxis()->SetRangeUser(0, 1);


  return theGraph;

}



//________________________________________________________
TGraph2DErrors * EfficiencyFitter::produce2DGraph(std::vector<std::vector<float> > x, std::vector<std::vector<float> > y, std::vector<std::vector<float> > e, std::vector<std::vector<float> > er, std::string xlabel, std::string ylabel, std::string zlabel, int col) {

  int rsize = x.size();
  int csize = x[0].size();

  double *theX = new double [rsize*csize];
  double *theY = new double [rsize*csize];
  double *theE = new double [rsize*csize];
  double *theEE = new double [rsize*csize];

  for(int c = 0; c < rsize; ++c) {
    for(int d = 0; d < csize; ++d) {

      theX[csize*c+d] = x[c].at(d);
      theY[csize*c+d] = y[c].at(d);
      theE[csize*c+d] = e[c].at(d);
      theEE[csize*c+d] = er[c].at(d);
  
    }
  }


  TGraph2DErrors *theGraph = new TGraph2DErrors(rsize*csize, theX, theY, theE, NULL, NULL, theEE);
  theGraph->GetXaxis()->SetTitle(xlabel.c_str());
  theGraph->GetYaxis()->SetTitle(ylabel.c_str());
  theGraph->GetZaxis()->SetTitle(zlabel.c_str());
  theGraph->SetMarkerColor(col);
  theGraph->SetTitle("");
  theGraph->GetZaxis()->SetRangeUser(0, 1);


  return theGraph;

}



//___----____----____----____----___----____----____----___----____----____----____----____----___
//                          *****Public methods***** 
//___----____----____----____----___----____----____----___----____----____----____----____----___



//_________________________________________________________
void EfficiencyFitter::doPlotAgainstVar(std::string var, std::vector<float> range, std::string xlabel, std::string graphid, std::string tC) {

  std::vector<float> eff;
  std::vector<float> eeff;
  std::vector<float> x;
  
  for(int c = 0; c < range.size()-1; c++) {
 
    std::string tCut = "(" + var + " > " + any2string(range.at(c)) + " && " + var + " < " + any2string(range.at(c+1)) + " && " + tC + ")";
    std::string id = graphid + "_" + var + "_" + any2string(c) + "_" + any2string(c+1);
    TCut provCut(tCut.c_str());
    float theEff, theEEff, theX;
    theX = (range.at(c)+range.at(c+1))/2.0;
    //theEff = 0.5;
    //theEEff = 0.0;
    calculateFitEfficiency(tCut, theEff, theEEff, id); 
    eff.push_back(theEff);
    eeff.push_back(theEEff);
    x.push_back(theX);

  }

  TGraphErrors * myGraph = produceGraph(x, eff, eeff, xlabel, "Efficiency", kBlack);
  
  TCanvas *theC = new TCanvas(graphid.c_str());
  theC->cd();
  myGraph->Draw("AP");
  theC->Write();

}

//_________________________________________________________
void EfficiencyFitter::doPlotAgainstVarMC(std::string var, std::vector<float> range, std::string xlabel, std::string graphid, std::string tC) {
  
  std::vector<float> eff;
  std::vector<float> eeff;
  std::vector<float> x; 

  for(int c = 0; c < range.size()-1; c++) { 

    std::string tCut = "(" + var + " > " + any2string(range.at(c)) + " && " + var + " < " + any2string(range.at(c+1)) + " && " + tC + ")";
    std::string id = graphid + "_" + var + "_" + any2string(c) + "_" + any2string(c+1);
    TCut provCut(tCut.c_str());
    float theEff, theEEff, theX;
    theX = (range.at(c)+range.at(c+1))/2.0;
    calculateMCNEfficiency(provCut, theEff, theEEff);
    eff.push_back(theEff);
    eeff.push_back(theEEff);
    x.push_back(theX);

  }

  TGraphErrors * myGraph = produceGraph(x, eff, eeff, xlabel, "Efficiency", kBlack);
    
  TCanvas *theC = new TCanvas(graphid.c_str());
  theC->cd();
  myGraph->Draw("AP");
  theC->Write();

}

//_________________________________________________________
void EfficiencyFitter::doPlotAgainstVarAndVar(std::string var1, std::string var2, std::vector<float> range1, std::vector<float> range2, 
                                              std::string xlabel, std::string ylabel, std::string graphid, std::string tC) {


  std::vector<std::vector<float> > eff;
  std::vector<std::vector<float> > eeff;
  std::vector<std::vector<float> > x;
  std::vector<std::vector<float> > y;

  for(int c = 0; c < range1.size()-1; c++) {
    std::vector<float> efi, efier, xi, yi;
    for(int d = 0; d < range2.size()-1; d++) {
      std::string tCut1 = "(" + var1 + " > " + any2string(range1.at(c)) + " && " + var1 + " < " + any2string(range1.at(c+1)) + " && " + tC + ")";
      std::string tCut2 = "(" + var2 + " > " + any2string(range2.at(d)) + " && " + var2 + " < " + any2string(range2.at(d+1)) + " && " + tC + ")";
      std::string tCut = "(" + tCut1 + "&&" + tCut2 + ")";
      std::string id = graphid + "_" + var1 + "_" + any2string(c) + "_" + any2string(c+1) +  "_" + var2 + "_" + any2string(d) + "_" + any2string(d+1);

      TCut provCut(tCut.c_str());
      float theEff, theEEff, theX, theY;
      theX = (range1.at(c)+range1.at(c+1))/2.0;
      theY = (range2.at(d)+range2.at(d+1))/2.0;
      //theEff = 0.7;
      //theEEff = 0.1;
      calculateFitEfficiency(tCut, theEff, theEEff, id);
      efi.push_back(theEff);
      efier.push_back(theEEff);
      xi.push_back(theX);
      yi.push_back(theY);
    }
    eff.push_back(efi);
    eeff.push_back(efier);
    x.push_back(xi);
    y.push_back(yi);
  }

  /*

  TGraph2DErrors * myGraph = produce2DGraph(x, y, eff, eeff, xlabel, ylabel, "Efficiency", kBlack);

  TCanvas *theC = new TCanvas(graphid.c_str());
  theC->cd();
  myGraph->GetXaxis()->SetTitle(xlabel.c_str());
  myGraph->GetYaxis()->SetTitle(ylabel.c_str());
  myGraph->GetZaxis()->SetTitle("Efficiency");
  myGraph->SetTitle("");
  myGraph->GetZaxis()->SetRangeUser(0, 1);
  myGraph->SetFillColor(29);
  myGraph->SetMarkerSize(0.8);
  myGraph->SetMarkerStyle(20);
  myGraph->SetMarkerColor(kRed);
  myGraph->SetLineColor(kBlue-3);
  myGraph->SetLineWidth(2);
  myGraph->Draw("COLZ");

  theC->Write();

  */


  std::string idtree = "tree_" + graphid;
  s = new TTree(idtree.c_str(), idtree.c_str());

  s->Branch("efficiency", &kefficiency, "efficiency/F");
  s->Branch("e_efficiency1", &ke_efficiency1, "e_efficiency1/F");
  s->Branch("e_efficiency2", &ke_efficiency2, "e_efficiency2/F");
  s->Branch("eta", &keta, "eta/F");
  s->Branch("pt", &kpt, "pt/F");
  s->Branch("dr", &kdr, "dr/F");
  s->Branch("njets", &knjets, "njets/F");
  s->Branch("ptwidth", &kptwidth, "ptwidth/F");
  s->Branch("etawidth", &ketawidth, "etawidth/F");
  s->Branch("drwidth", &kdrwidth, "drwidth/F");
  s->Branch("njetswidth", &knjetswidth, "njetswidth/F");


  for(int c = 0; c < range1.size()-1; c++) {
    for(int d = 0; d < range2.size()-1; d++) {
      std::vector<float> theeff = eff.at(c);
      std::vector<float> theeeff = eeff.at(c);
      kefficiency = theeff.at(d);
      ke_efficiency1 = theeeff.at(d);
      ke_efficiency2 = theeeff.at(d);
      kpt = (range1.at(c)+range1.at(c+1))/2.0;
      kdr = (range2.at(d)+range2.at(d+1))/2.0;
      knjets = 2;
      keta = 0;
      kptwidth = range1.at(c)-range1.at(c+1);
      kdrwidth = range2.at(d)-range2.at(d+1);
      ketawidth = 1;
      knjetswidth = 1;
      s->Fill();
    }
  }

  s->Write();


  
}



//_________________________________________________________
void EfficiencyFitter::doPlotAgainstVarAndVarMC(std::string var1, std::string var2, std::vector<float> range1, std::vector<float> range2,
                                              std::string xlabel, std::string ylabel, std::string graphid, std::string tC) {


  std::vector<std::vector<float> > eff;
  std::vector<std::vector<float> > eeff;
  std::vector<std::vector<float> > x;
  std::vector<std::vector<float> > y;

  for(int c = 0; c < range1.size()-1; c++) {
    std::vector<float> efi, efier, xi, yi;
    for(int d = 0; d < range2.size()-1; d++) {
      std::string tCut1 = "(" + var1 + " > " + any2string(range1.at(c)) + " && " + var1 + " < " + any2string(range1.at(c+1)) + " && " + tC + ")";
      std::string tCut2 = "(" + var2 + " > " + any2string(range2.at(d)) + " && " + var2 + " < " + any2string(range2.at(d+1)) + " && " + tC + ")";
      std::string tCut = "(" + tCut1 + "&&" + tCut2 + ")";
      std::string id = graphid + "_" + var1 + "_" + any2string(c) + "_" + any2string(c+1) +  "_" + var2 + "_" + any2string(d) + "_" + any2string(d+1);

      TCut provCut(tCut.c_str());
      float theEff, theEEff, theX, theY;
      theX = (range1.at(c)+range1.at(c+1))/2.0;
      theY = (range2.at(d)+range2.at(d+1))/2.0;
      calculateMCNEfficiency(provCut, theEff, theEEff);
      efi.push_back(theEff);
      efier.push_back(theEEff);
      xi.push_back(theX);
      yi.push_back(theY);
    }
    eff.push_back(efi);
    eeff.push_back(efier);
    x.push_back(xi);
    y.push_back(yi);
  }


  std::string idtree = "tree_" + graphid;
  s = new TTree(idtree.c_str(), idtree.c_str());

  s->Branch("efficiency", &kefficiency, "efficiency/F");
  s->Branch("e_efficiency1", &ke_efficiency1, "e_efficiency1/F");
  s->Branch("e_efficiency2", &ke_efficiency2, "e_efficiency2/F");
  s->Branch("eta", &keta, "eta/F");
  s->Branch("pt", &kpt, "pt/F");
  s->Branch("dr", &kdr, "dr/F");
  s->Branch("njets", &knjets, "njets/F");
  s->Branch("ptwidth", &kptwidth, "ptwidth/F");
  s->Branch("etawidth", &ketawidth, "etawidth/F");
  s->Branch("drwidth", &kdrwidth, "drwidth/F");
  s->Branch("njetswidth", &knjetswidth, "njetswidth/F");


  for(int c = 0; c < range1.size()-1; c++) {
    for(int d = 0; d < range2.size()-1; d++) {
      std::vector<float> theeff = eff.at(c);
      std::vector<float> theeeff = eeff.at(c);
      kefficiency = theeff.at(d);
      ke_efficiency1 = theeeff.at(d);
      ke_efficiency2 = theeeff.at(d);
      kpt = (range1.at(c)+range1.at(c+1))/2.0;
      kdr = (range2.at(d)+range2.at(d+1))/2.0;
      knjets = 2;
      keta = 0;
      kptwidth = range1.at(c)-range1.at(c+1);
      kdrwidth = range2.at(d)-range2.at(d+1);
      ketawidth = 1;
      knjetswidth = 1;
      s->Fill();
    }
  }

  s->Write();

}


//_________________________________________________________
void EfficiencyFitter::doTriggerPlotAgainstVar(std::string var, std::vector<float> range, std::string xlabel, std::string graphid, std::string tC) {

  std::vector<float> eff;
  std::vector<float> eeff;
  std::vector<float> x;

  for(int c = 0; c < range.size()-1; c++) {

    std::string tCut = "(" + var + " > " + any2string(range.at(c)) + " && " + var + " < " + any2string(range.at(c+1)) + " && " + tC + ")";
    std::string id = graphid + "_" + var + "_" + any2string(c) + "_" + any2string(c+1);
    TCut provCut(tCut.c_str());
    float theEff, theEEff, theX;
    theX = (range.at(c)+range.at(c+1))/2.0;
    calculateTriggerEfficiency(provCut, theEff, theEEff, id);
    eff.push_back(theEff);
    eeff.push_back(theEEff);
    x.push_back(theX);

  }

  TGraphErrors * myGraph = produceGraph(x, eff, eeff, xlabel, "Efficiency", kBlack);

  std::string thegraphid = "trigger_" + graphid;
  TCanvas *theC = new TCanvas(thegraphid.c_str());
  theC->cd();
  myGraph->Draw("AP");
  theC->Write();


}


