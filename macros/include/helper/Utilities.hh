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
	inline void SetTDRStyle(){
		TStyle *tdrStyle = new TStyle("tdrStyle","Style for P-TDR");

		// For the canvas:
		tdrStyle->SetCanvasBorderMode(0);
		tdrStyle->SetCanvasColor(kWhite);
		tdrStyle->SetCanvasDefH(600); //Height of canvas
		tdrStyle->SetCanvasDefW(600); //Width of canvas
		tdrStyle->SetCanvasDefX(0);   //POsition on screen
		tdrStyle->SetCanvasDefY(0);

		// For the Pad:
		tdrStyle->SetPadBorderMode(0);
		// tdrStyle->SetPadBorderSize(Width_t size = 1);
		tdrStyle->SetPadColor(kWhite);
		tdrStyle->SetPadGridX(false);
		tdrStyle->SetPadGridY(false);
		tdrStyle->SetGridColor(0);
		tdrStyle->SetGridStyle(3);
		tdrStyle->SetGridWidth(1);

		// For the frame:
		tdrStyle->SetFrameBorderMode(0);
		tdrStyle->SetFrameBorderSize(1);
		tdrStyle->SetFrameFillColor(0);
		tdrStyle->SetFrameFillStyle(0);
		tdrStyle->SetFrameLineColor(1);
		tdrStyle->SetFrameLineStyle(1);
		tdrStyle->SetFrameLineWidth(1);

		// For the histo:
		// tdrStyle->SetHistFillColor(1);
		// tdrStyle->SetHistFillStyle(0);
		tdrStyle->SetHistLineColor(1);
		tdrStyle->SetHistLineStyle(0);
		tdrStyle->SetHistLineWidth(1);
		// tdrStyle->SetLegoInnerR(Float_t rad = 0.5);
		// tdrStyle->SetNumberContours(Int_t number = 20);

		tdrStyle->SetEndErrorSize(2);
		// tdrStyle->SetErrorMarker(20);
		tdrStyle->SetErrorX(0.);

		tdrStyle->SetMarkerStyle(20);

		//For the fit/function:
		tdrStyle->SetOptFit(1);
		tdrStyle->SetFitFormat("5.4g");
		tdrStyle->SetFuncColor(2);
		tdrStyle->SetFuncStyle(1);
		tdrStyle->SetFuncWidth(1);

		//For the date:
		tdrStyle->SetOptDate(0);
		// tdrStyle->SetDateX(Float_t x = 0.01);
		// tdrStyle->SetDateY(Float_t y = 0.01);

		// For the statistics box:
		tdrStyle->SetOptFile(0);
		tdrStyle->SetOptStat(0); // To display the mean and RMS:   SetOptStat("mr");
		tdrStyle->SetStatColor(kWhite);
		tdrStyle->SetStatFont(42);
		tdrStyle->SetStatFontSize(0.025);
		tdrStyle->SetStatTextColor(1);
		tdrStyle->SetStatFormat("6.4g");
		tdrStyle->SetStatBorderSize(1);
		tdrStyle->SetStatH(0.1);
		tdrStyle->SetStatW(0.15);
		// tdrStyle->SetStatStyle(Style_t style = 1001);
		// tdrStyle->SetStatX(Float_t x = 0);
		// tdrStyle->SetStatY(Float_t y = 0);

		// Margins:
		tdrStyle->SetPadTopMargin(0.05);
		tdrStyle->SetPadBottomMargin(0.13);
		tdrStyle->SetPadLeftMargin(0.16);
		tdrStyle->SetPadRightMargin(0.02);

		// For the Global title:
		tdrStyle->SetOptTitle(0);
		tdrStyle->SetTitleFont(42);
		tdrStyle->SetTitleColor(1);
		tdrStyle->SetTitleTextColor(1);
		tdrStyle->SetTitleFillColor(10);
		tdrStyle->SetTitleFontSize(0.05);
		// tdrStyle->SetTitleH(0); // Set the height of the title box
		// tdrStyle->SetTitleW(0); // Set the width of the title box
		// tdrStyle->SetTitleX(0); // Set the position of the title box
		// tdrStyle->SetTitleY(0.985); // Set the position of the title box
		// tdrStyle->SetTitleStyle(Style_t style = 1001);
		// tdrStyle->SetTitleBorderSize(2);

		// For the axis titles:
		tdrStyle->SetTitleColor(1, "XYZ");
		tdrStyle->SetTitleFont(42, "XYZ");
		tdrStyle->SetTitleSize(0.06, "XYZ");
		// tdrStyle->SetTitleXSize(Float_t size = 0.02); // Another way to set the size?
		// tdrStyle->SetTitleYSize(Float_t size = 0.02);
		tdrStyle->SetTitleXOffset(0.9);
		tdrStyle->SetTitleYOffset(1.25);
		// tdrStyle->SetTitleOffset(1.1, "Y"); // Another way to set the Offset

		// For the axis labels:
		tdrStyle->SetLabelColor(1, "XYZ");
		tdrStyle->SetLabelFont(42, "XYZ");
		tdrStyle->SetLabelOffset(0.007, "XYZ");
		tdrStyle->SetLabelSize(0.05, "XYZ");

		// For the axis:
		tdrStyle->SetAxisColor(1, "XYZ");
		tdrStyle->SetStripDecimals(kTRUE);
		tdrStyle->SetTickLength(0.03, "XYZ");
		tdrStyle->SetNdivisions(510, "XYZ");
		tdrStyle->SetPadTickX(1);  // To get tick marks on the opposite side of the frame
		tdrStyle->SetPadTickY(1);

		// Change for log plots:
		tdrStyle->SetOptLogx(0);
		tdrStyle->SetOptLogy(0);
		tdrStyle->SetOptLogz(0);

		// Postscript options:
		tdrStyle->SetPaperSize(20.,20.);
		// tdrStyle->SetLineScalePS(Float_t scale = 3);
		// tdrStyle->SetLineStyleString(Int_t i, const char* text);
		// tdrStyle->SetHeaderPS(const char* header);
		// tdrStyle->SetTitlePS(const char* pstitle);

		// tdrStyle->SetBarOffset(Float_t baroff = 0.5);
		// tdrStyle->SetBarWidth(Float_t barwidth = 0.5);
		// tdrStyle->SetPaintTextFormat(const char* format = "g");
		// tdrStyle->SetPalette(Int_t ncolors = 0, Int_t* colors = 0);
		// tdrStyle->SetTimeOffset(Double_t toffset);
		// tdrStyle->SetHistMinimumZero(kTRUE);

		tdrStyle->cd();

		gROOT->SetStyle("tdrStyle");
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
		// Saves a ROOT TCanvas Object to a .C file
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
		// ind.reserve(collectionSize); // better to initialize to a code value
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
	inline double IntegralAndError(TH1 *hist, int bin1, int bin2, double &err){
		// not implemented before ROOT v2.26
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
