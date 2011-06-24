#ifndef Utilities_hh
#define Utilities_hh

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <vector>
#include <string>
#include <cmath>

#include "TROOT.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TString.h"
#include "TMath.h"
#include "TFile.h"
#include "TH1.h"


namespace Util {
  //__________________________________________________________________________
  inline TString MakeOutputDir(TString dir){
    if(!dir.EndsWith("/")) dir += "/";
    // Create directory if needed
    //  >> NOTE: This function needs to be called before the booking functions!
    char cmd[100];
    sprintf(cmd,"mkdir -p %s", dir.Data());
    system(cmd);
    return dir;
  }

  //__________________________________________________________________________
  inline TFile* MakeOutputFile(TString filename){
    if(!filename.EndsWith(".root")) filename += ".root";
    TFile *file = new TFile(filename, "RECREATE");
    return file;
  }

  //__________________________________________________________________________
  inline void SetStyle(){
    TStyle *style = new TStyle("ETHStyle", "Standard Plain");
    style->SetCanvasColor(0);
    style->SetFrameFillColor(0);
    style->SetFrameBorderMode(0);
    style->SetFrameBorderSize(0);
    style->SetPalette(1,0);
    style->SetOptTitle(0);
    style->SetOptStat(111111);
    style->SetStatColor(0);
    style->SetStatStyle(3001);
    style->SetStatBorderSize(1);

    // Fonts
    Int_t font = 42;
    style->SetStatFont(font);
    style->SetTextFont(font);
    style->SetLabelFont(font, "xyz");
    style->SetTitleFont(font, "xyz");

    // Histograms
    style->SetHistFillColor(15);
    style->SetHistFillStyle(1001);
    style->SetHistLineWidth(2);
    gROOT->SetStyle("ETHStyle");
    gROOT->ForceStyle();
  }

 //__________________________________________________________________________
  inline void PrintPNG(TCanvas *cin, TString name, TString dir){
    // Prints a ROOT TCanvas Object to a .png file
    //  name is the bare output filename, e.g. "fit_4_8",
    //  dir is the output directory (inside the overall output dir.)
    // Create sub directories if needed
	dir = MakeOutputDir(dir);

    dir += name;
    dir += ".png";
    cin->Print(dir,"png");
  }

  //__________________________________________________________________________
  inline void PrintPDF(TCanvas *cin, TString name, TString dir){
    // Prints a ROOT TCanvas Object to a .png file
    //  name is the bare output filename, e.g. "fit_4_8",
    //  dir is the output directory (inside the overall output dir.)
    // Create sub directories if needed
	dir = MakeOutputDir(dir);

    dir += name;
    dir += ".pdf";
    cin->Print(dir,"pdf");
  }

  //__________________________________________________________________________
  inline void PrintEPS(TCanvas *cin, TString name, TString dir){
    // Prints a ROOT TCanvas Object to a .eps file
    //  name is the bare output filename, e.g. "fit_4_8",
    //  dir is the output directory (inside the overall output dir.)
    // Create sub directories if needed
	dir = MakeOutputDir(dir);

    dir += name;
    dir += ".eps";
    cin->Print(dir,"eps");
  }

  //__________________________________________________________________________
  inline void SaveAsMacro(TCanvas *cin, TString name, TString dir){
    // Prints a ROOT TCanvas Object to a .eps file
    //  name is the bare output filename, e.g. "fit_4_8",
    //  dir is the output directory (inside the overall output dir.)
    // Create sub directories if needed
	dir = MakeOutputDir(dir);

    dir += name;
    dir += ".C";
    cin->SaveAs(dir);
  }

  //__________________________________________________________________________
  inline TDirectory* FindOrCreate( TString& dir, TFile* file ) {
    // Look for a directory and create it if it does not exist

    // Start from the root of the file
    file->cd();
    // Remove deadly '/'s
    while ( dir.BeginsWith("/") ) dir = dir.Strip(TString::kLeading,'/');
    dir.ReplaceAll("//","/");

    // Loop over sub-directories to create (ROOT's mkdir has no -p option...)
    TString cdir(dir);
    while ( cdir.First('/')>0 || cdir.Length()>0 ) {
       // Create new subdirectory
       Size_t index = (cdir.First('/')>0 ? cdir.First('/') : cdir.Length());
       TString subdir = cdir(0,index);
       if ( !gDirectory->GetDirectory(subdir) ) {
          std::cout << "Creating directory " << subdir.Data() << std::endl;
          gDirectory->mkdir( subdir.Data() );
       }
       gDirectory->cd(subdir);
       cdir = cdir(index+1,cdir.Length());
    }
    return file->GetDirectory(dir);
    
  }
  
  //__________________________________________________________________________
  inline void SaveAll(TCanvas *cin, TString dir, TFile* file) {
    // Save all objects in a canvas to a file
    //   dir is a sub-directory in the file
    //   file is the file object (need to be already open)
    
    // A few checks
    if ( !file || !file->IsOpen() ) {
      std::cerr << "*** Util::SaveAll: file " << (file?file->GetName():"") << " does not exist" << std::endl;
      exit(-1);
    } 

    // Go to directory (create it if needed)
    TDirectory* cdir = Util::FindOrCreate(dir,file);
    if ( !cdir) {
      std::cerr << "Couldn't create directory " << dir << std::endl;
      exit(-1);
    }
    cdir->cd();
    
    // Loop over canvas object and save some of them
    TIter next(cin->GetListOfPrimitives());
    while (TObject *obj = next()) {
      if ( !strcmp(obj->ClassName(),"TFrame") ) continue;
      if ( !strcmp(obj->ClassName(),"TLine") ) continue;
      if ( !strcmp(obj->ClassName(),"TArrow") ) continue;
      if ( !strcmp(obj->ClassName(),"TLatex") ) continue;
      obj->Write(obj->GetName(),TObject::kOverwrite);
    }
  }

  //__________________________________________________________________________
  inline void Print(TCanvas *cin, TString name, TString dir, TFile* file=0) {
    // Print plot (PNG, EPS and to file)
    Util::PrintPNG(cin,name,dir);
    Util::PrintEPS(cin,name,dir);
    if ( file ) Util::SaveAll(cin,dir,file); 
  }

  //__________________________________________________________________________
  inline void PrintNoEPS(TCanvas *cin, TString name, TString dir, TFile* file=0) {
    // Print plot (PNG and to file)
    Util::PrintPNG(cin,name,dir);
    if ( file ) Util::SaveAll(cin,dir,file);
  }

  //__________________________________________________________________________
  inline double DeltaPhi(double phi1, double phi2){
    // From cmssw reco::deltaPhi()
    double result = phi1 - phi2;
    while( result >   TMath::Pi() ) result -= TMath::TwoPi();
    while( result <= -TMath::Pi() ) result += TMath::TwoPi();
    return TMath::Abs(result);
  }

  //__________________________________________________________________________
  inline double GetDeltaR(double eta1, double eta2, double phi1, double phi2){
    double deta = eta1 - eta2;
    double dphi = Util::DeltaPhi(phi1, phi2);
    return sqrt( deta*deta + dphi*dphi );
  }

  //__________________________________________________________________________
  template<class T> inline std::vector<int> VSort(std::vector<T> vec, bool asc = false){
    // Sort a vector and return the vector of sorted indices
    // Simple bubble sort algorithm, don't use for more than a few entries!
    std::vector<int> ind;
    if(vec.size() == 0) return ind; // Return original empty vector

    std::vector<T> vecClone = vec; // clone orignal vector
    sort(vecClone.begin(), vecClone.end()); // ascending order
    if(!asc)reverse(vecClone.begin(), vecClone.end());  // sort in descending order  

    int collectionSize = vec.size();
//    ind.reserve(collectionSize); // better to initialize to a code value
    for(int i=0;i<collectionSize;i++)ind.push_back(-999);
   
    
    for(int i =0;i<collectionSize;i++) // loop in the sorted collection
    {
	for(int j =0;j<collectionSize;j++) // loop in the original unsorted collection
	{
	    if(vecClone[i]==vec[j]){ind[i]=j;}
	}
    }

    bool success = true; // failsafe test 
    for(int i=0;i<collectionSize;i++)
    {
       	if(ind[i]==-999)success=false;
    }

    if(!success)
    {
	std::cout<<"problem with the sorting"<<std::endl;
	std::vector<int> dummy; 
	return dummy;
    }
    else
    {
	return ind;
    }

  }

  //__________________________________________________________________________
  template<class T> inline std::vector<int> VSort(std::vector<int> ind, std::vector<T> vec, bool asc = false){
    // Sort a vector of ints (ind) according to the values in a second vector (vec)
    // of the same length
    // Simple bubble sort algorithm, don't use for more than a few entries!
    if(ind.size()!=vec.size()){
      std::cout << "Util::VSort ==> Vectors don't match in size! Returning unsorted vector..." << std::endl;
      return ind;
    }
    if(ind.size() == 0) return ind; // Return original empty vector
    std::vector<int> ind2 = VSort(vec, asc);
    std::vector<int> ind3;
    for(size_t i = 0; i < vec.size(); ++i) ind3.push_back(ind[ind2[i]]);
    return ind3;
  }

  //__________________________________________________________________________
  inline std::string removeFunnyChar(const std::string& input){
    const std::string excluded = " ()[]/\\~;:{}*&$<>`!@#%^+|\'\",?";
    std::string answer(input);
    std::string::size_type pos = answer.find_first_of(excluded);
    while (pos != std::string::npos) {
      answer.erase(pos, 1);
      pos = answer.find_first_of(excluded);
    }    
    return answer;
  }
  //__________________________________________________________________________
  inline double IntegralAndError(TH1 *hist, int bin1, int bin2, double &err){  // not implemented before ROOT v2.26
    double_t integral = 0;
    double_t igerr2 = 0;
    for (Int_t bin = bin1; bin <= bin2; ++bin) {
      integral += hist->GetBinContent(bin);
      igerr2   += hist->GetBinError(bin)*hist->GetBinError(bin);
      
    }
    err = TMath::Sqrt(igerr2);
    return integral;
  }
}

#endif
