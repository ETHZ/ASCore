#ifndef Lumi3DReWeighting_standalone_hh
#define Lumi3DReWeighting_standalone_hh

#include "TH1.h"
#include "TFile.h"
#include <cmath>
#include <string>
#include <vector>
#include <iostream>


class Lumi3DReWeighting{

public:
  Lumi3DReWeighting( std::string generatedFile,
		     std::string dataFile,
		     std::string GenHistName,
		     std::string DataHistName);
  
  Lumi3DReWeighting( std::vector< float > MC_distr, std::vector< float > Lumi_distr);
  
  Lumi3DReWeighting ( ) { } ;
  
  double weight3D( int, int, int );

  void weight3D_set( std::string generatedFile, std::string dataFile, std::string GenHistName, std::string DataHistName);

  void weight3D_init( float Scale );

  void weight3D_init( std::string WeightFileName );  // initialize from root file

  void weight3D_init( std::string MCFileName, std::string DataFileName );  // initialize from root files


protected:

  std::string generatedFileName_;
  std::string dataFileName_;
  std::string GenHistName_;
  std::string DataHistName_;
  TFile*     generatedFile_;
  TFile*     dataFile_;
  TH1*      weights_;

  //keep copies of normalized distributions:

  TH1*      MC_distr_;
  TH1*      Data_distr_;


  double Weight3D_[50][50][50];


};


#endif
