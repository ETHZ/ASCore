#include "EnergyCorrection.hh"
#include "base/TreeReader.hh"

EnergyCorrection::EnergyCorrection(TString tuning)
{
  
  {
    Double_t       leftEtatemp  [nBinsEta]   = { 0.02, 0.25, 0.46, 0.81, 0.91, 1.01, 1.16,           etaCrackMax,  1.653,  1.8, 2.0, 2.2, 2.3, 2.4 };
    Double_t       rightEtatemp [nBinsEta]   = { 0.25, 0.42, 0.77, 0.91, 1.01, 1.13, etaCrackMin,    1.653,        1.8  ,  2.0, 2.2, 2.3, 2.4, 2.5 };
    
    for (int i=0; i<nBinsEta; i++) leftEta[i]=leftEtatemp[i];
    for (int i=0; i<nBinsEta; i++) rightEta[i]=rightEtatemp[i];
  }
  
  if (tuning=="photons") forphotons=true;
  else if (tuning=="electrons") forphotons=false;
  else {
    std::cout << "Wrong initialization of EnergyCorrection!" << std::endl;
  }

  if (forphotons){


    xcorr[0]=1.00506;
    xcorr[1]=1.00697;
    xcorr[2]=1.00595;
    xcorr[3]=1.00595;
    xcorr[4]=1.00595;
    xcorr[5]=1.00595;
    xcorr[6]=1.00595;
    xcorr[7]=0.966651;
    xcorr[8]=0.97381;
    xcorr[9]=0.976516;
    xcorr[10]=0.983254;
    xcorr[11]=0.98502;
    xcorr[12]=0.98502;
    xcorr[13]=0.978472;

    par0[0] = 0.00132382 ;
    par1[0] = 2.17664 ;
    par2[0] = -0.00467206 ;
    par3[0] = 0.988994 ;
    par4[0] = 17.5858 ;

    par0[1] = -0.00590257 ;
    par1[1] = 1.90733 ;
    par2[1] = 0.000684327 ;
    par3[1] = 0.986431 ;
    par4[1] = 16.6698 ;

    par0[2] = 0.00265109 ;
    par1[2] = 1.73272 ;
    par2[2] = -0.00107022 ;
    par3[2] = 0.989322 ;
    par4[2] = 15.4911 ;

    par0[3] = 0.00231631 ;
    par1[3] = 1.3463 ;
    par2[3] = -0.00369555 ;
    par3[3] = 0.987133 ;
    par4[3] = 10.9233 ;

    par0[4] = 0.00984253 ;
    par1[4] = 1.33889 ;
    par2[4] = -0.00392593 ;
    par3[4] = 0.979191 ;
    par4[4] = 9.35276 ;

    par0[5] = 0.023683 ;
    par1[5] = 1.31198 ;
    par2[5] = -0.00947317 ;
    par3[5] = 0.963352 ;
    par4[5] = 7.5597 ;

    par0[6] = 0.0851133 ;
    par1[6] = 1.38097 ;
    par2[6] = -0.0340201 ;
    par3[6] = 0.969502 ;
    par4[6] = 4.17983 ;

    par0[7] = 6.71705 ;
    par1[7] = 5034.26 ;
    par2[7] = -2.68669 ;
    par3[7] = 0.970174 ;
    par4[7] = 1.00288 ;

    par0[8] = 1306.82 ;
    par1[8] = 472004 ;
    par2[8] = -1.86145 ;
    par3[8] = 0.981714 ;
    par4[8] = -0.25644 ;

    par0[9] = 0.317121 ;
    par1[9] = 3.22717 ;
    par2[9] = -0.126848 ;
    par3[9] = 0.957792 ;
    par4[9] = 2.01028 ;

    par0[10] = 0.275225 ;
    par1[10] = 2.20686 ;
    par2[10] = -0.11009 ;
    par3[10] = 0.93922 ;
    par4[10] = 2.69958 ;

    par0[11] = 0.0639875 ;
    par1[11] = 1.40045 ;
    par2[11] = -0.0255853 ;
    par3[11] = 0.821566 ;
    par4[11] = 7.3297 ;

    par0[12] = 0.030488 ;
    par1[12] = 1.37842 ;
    par2[12] = -0.0121879 ;
    par3[12] = 0.8173 ;
    par4[12] = 9.29944 ;

    par0[13] = 0.213906 ;
    par1[13] = 1.67471 ;
    par2[13] = -0.0860589 ;
    par3[13] = 0.893636 ;
    par4[13] = 3.78218 ;

  } // end for photons

  else { // for electrons

    xcorr[0]=1.00227;
    xcorr[1]=1.00252;
    xcorr[2]=1.00225;
    xcorr[3]=1.00159;
    xcorr[4]=0.999475;
    xcorr[5]=0.997203;
    xcorr[6]=0.993886;
    xcorr[7]=0.971262;
    xcorr[8]=0.975922;
    xcorr[9]=0.979087;
    xcorr[10]=0.98495;
    xcorr[11]=0.98781;
    xcorr[12]=0.989546;
    xcorr[13]=0.989638;

    par0[0] = 1.00718;
    par1[0] = -0.00187886;
    par2[0] = 0 ;

    par0[1] = 1.00713;
    par1[1] = -0.00227574;
    par2[1] = 0 ;

    par0[2] = 1.00641;
    par1[2] = -0.00259935;
    par2[2] = 0 ;

    par0[3] = 1.00761;
    par1[3] = -0.00433692;
    par2[3] = 0 ;

    par0[4] = 1.00682;
    par1[4] = -0.00551324;
    par2[4] = 0 ;

    par0[5] = 1.0073;
    par1[5] = -0.00799669;
    par2[5] = 0 ;

    par0[6] = 1.00462;
    par1[6] = -0.00870057;
    par2[6] = 0 ;

    par0[7] = 0.972798;
    par1[7] = -0.000771577;
    par2[7] = -0.00276696;

    par0[8] = 0.981672;
    par1[8] = -0.00202028;
    par2[8] = -0.00471028;

    par0[9] = 0.98251;
    par1[9] = 0.00441308;
    par2[9] = -0.00809139;

    par0[10] = 0.986123;
    par1[10] = 0.00832913;
    par2[10] = -0.00944584;

    par0[11] = 0.990124;
    par1[11] = 0.00742879;
    par2[11] = -0.00960462;

    par0[12] = 0.990187;
    par1[12] = 0.0094608;
    par2[12] = -0.010172;

    par0[13] = 0.99372;
    par1[13] = 0.00560406;
    par2[13] = -0.00943169;

  }

};

 

EnergyCorrection::~EnergyCorrection() 
{
};

Double_t EnergyCorrection::applyScCorrectionsBrEta_photons(Double_t eta, Double_t sigmaPhiSigmaEta){  

  // extra protections																					   
  // fix sigmaPhiSigmaEta boundaries 
  if (sigmaPhiSigmaEta < 0.8)  sigmaPhiSigmaEta = 0.8; 
  if (sigmaPhiSigmaEta > 5  )  sigmaPhiSigmaEta = 5; 

  // eta = 0																						   
  if (TMath::Abs(eta)  <  leftEta[0]            ) { eta = 0.02 ; }																   
  // outside acceptance																					   
  if (TMath::Abs(eta)  >=  rightEta[nBinsEta-1] ) { eta = 2.49; if (DBG) std::cout << " WARNING [applyScCorrections]: TMath::Abs(eta)  >=  rightEta[nBinsEta-1] " << std::endl;}  
  																								   
  Int_t tmpEta = -1;                                                                                                                                                                         
  for (Int_t iEta = 0; iEta < nBinsEta; ++iEta){								              								      	   
    if ( leftEta[iEta] <= TMath::Abs(eta) && TMath::Abs(eta) <rightEta[iEta] ){				       									      	   
      tmpEta = iEta;											       										   
    }													       										   
  }													       										           

  // Interpolation													         
  Double_t tmpInter = 1;													         
  // In eta cracks/gaps 													         
  if (tmpEta == -1 ) { // need to interpolate    
    for (Int_t iEta = 0; iEta < nBinsEta-1; ++iEta){									         
      if (rightEta[iEta] <= TMath::Abs(eta) && TMath::Abs(eta) <leftEta[iEta+1] ){					         
	if (sigmaPhiSigmaEta >= 1)  tmpInter =   (par0[iEta  ]*(1.-exp(-(sigmaPhiSigmaEta-par4[iEta  ])/par1[iEta  ]))*par2[iEta  ]*sigmaPhiSigmaEta + par3[iEta  ]+
						  par0[iEta+1]*(1.-exp(-(sigmaPhiSigmaEta-par4[iEta+1])/par1[iEta+1]))*par2[iEta+1]*sigmaPhiSigmaEta + par3[iEta+1] ) /2.;
	// ( fcorr[iEta].Eval(sigmaPhiSigmaEta) +fcorr[iEta+1].Eval(sigmaPhiSigmaEta) ) / 2. ;
	else tmpInter = (xcorr[iEta] + xcorr[iEta+1])/2.; 
      }															         
    }															         
    return tmpInter;													         
  }  															         
  if (sigmaPhiSigmaEta >= 1) return par0[tmpEta  ]*(1.-exp(-(sigmaPhiSigmaEta-par4[tmpEta  ])/par1[tmpEta  ]))*par2[tmpEta  ]*sigmaPhiSigmaEta + par3[tmpEta  ]; // fcorr[tmpEta].Eval(sigmaPhiSigmaEta);  
  else return xcorr[tmpEta]; 

};

Double_t EnergyCorrection::applyScCorrectionsET_EB_photons(Double_t ET){      

  Double_t par0 =  1; 
  Double_t par1 =  1.00348; 
  Double_t par2 =  1.001;	   
  Double_t par3 =  -9.17302e-06;	   
  Double_t par4 =  0.999688;     

  if (             ET <   5 ) return         1.;  
  if (  5 <= ET && ET <  10 ) return         par0 ;  
  if ( 10 <= ET && ET <  20 ) return         par1 ;  
  if ( 20 <= ET && ET < 140 ) return         par2 + par3*ET ;  
  if (140 <= ET             ) return         par4;  
 						  
};                               

Double_t EnergyCorrection::applyScCorrectionsET_EE_photons(Double_t ET){      
   					  
  Double_t par0 =  1; 
  Double_t par1 =  0.996931; 
  Double_t par2 =  0.999497;	   
  Double_t par3 =  0.992617;	   
  Double_t par4 =  7.52128e-05;     
  Double_t par5 =  -1.2845e-07;     
  Double_t par6 =  1.00231;     

  if (             ET <   5 ) return         1.;  
  if (  5 <= ET && ET <  10 ) return          par0 ;  
  if ( 10 <= ET && ET <  20 ) return          par1 ;  
  if ( 20 <= ET && ET <  30 ) return          par2 ;  
  if ( 30 <= ET && ET < 200 ) return          par3 + par4 *ET + par5 *ET*ET ;  
  if ( 200 <= ET            ) return          par6 ;
 						  
};

Double_t EnergyCorrection::applyScCorrectionsE_EE_photons(Double_t E){      
  				 	  
  Double_t par0 = 850;               
  Double_t par1 = 0.994169 ;	  
  Double_t par2 = 1.28629e-05 ;     
  				 	  
  if (E  > par0 ) E = par0 ;   		  
  if (            E <   0     ) return      1.;  
  if (  0 <= E && E <=  par0  ) return      par1 + E*par2; 
						   
};  

Double_t EnergyCorrection::applyScCorrectionsBrEta_electrons(Double_t eta, Double_t sigmaPhiSigmaEta){

  // extra protections																					   
  // fix sigmaPhiSigmaEta boundaries 
  if (sigmaPhiSigmaEta < 0.8)  sigmaPhiSigmaEta = 0.8; 
  if (sigmaPhiSigmaEta > 5  )  sigmaPhiSigmaEta = 5; 

  // eta = 0																						   
  if (TMath::Abs(eta)  <  leftEta[0]            ) { eta = 0.02 ; }																   
  // outside acceptance																					   
  if (TMath::Abs(eta)  >=  rightEta[nBinsEta-1] ) { eta = 2.49; if (DBG) std::cout << " WARNING [applyScCorrections]: TMath::Abs(eta)  >=  rightEta[nBinsEta-1] " << std::endl;}  
  																								   
  Int_t tmpEta = -1;                                                                                                                                                                         
  for (Int_t iEta = 0; iEta < nBinsEta; ++iEta){								              								      	   
    if ( leftEta[iEta] <= TMath::Abs(eta) && TMath::Abs(eta) <rightEta[iEta] ){				       									      	   
      tmpEta = iEta;											       										   
    }													       										   
  }													       										           

  // Interpolation																					         
  Double_t tmpInter = 1;																				         
  // In eta cracks/gaps 																				         
  if (tmpEta == -1 ) { // need to interpolate    
    for (Int_t iEta = 0; iEta < nBinsEta-1; ++iEta){								       								         
      if (rightEta[iEta] <= TMath::Abs(eta) && TMath::Abs(eta) <leftEta[iEta+1] ){													         
	if (sigmaPhiSigmaEta >= 1.2)  tmpInter = ( par0[iEta] + sigmaPhiSigmaEta*par1[iEta] + sigmaPhiSigmaEta*sigmaPhiSigmaEta*par2[iEta] +  
						   par0[iEta+1] + sigmaPhiSigmaEta*par1[iEta+1] + sigmaPhiSigmaEta*sigmaPhiSigmaEta*par2[iEta+1]) / 2. ; 
	else tmpInter = (xcorr[iEta] + xcorr[iEta+1])/2.; 
      }																						         
    }																						         
    return tmpInter;																					         
  }  																							         
  if (sigmaPhiSigmaEta >= 1.2) return par0[tmpEta] + sigmaPhiSigmaEta*par1[tmpEta] + sigmaPhiSigmaEta*sigmaPhiSigmaEta*par2[tmpEta]; 
  else return xcorr[tmpEta]; 

};

Double_t EnergyCorrection::applyScCorrectionsET_EB_electrons(Double_t ET){      
  
  Double_t par0 = 0.97213; 
  Double_t par1 = 0.999528; 
  Double_t par2 = 5.61192e-06; 
  Double_t par3 = 0.0143269; 
  Double_t par4 = -17.1776; 

  if (ET > 200) ET =200;   		  
  if (             ET <    5 ) return         1.;  
  if (  5 <= ET && ET <   10 ) return         par0 ;  
  if ( 10 <= ET && ET <= 200 ) return         (par1  + ET*par2)*(1- par3*exp(ET/ par4));
 						  
}    

Double_t EnergyCorrection::applyScCorrectionsET_EE_electrons(Double_t ET){      
  
  Double_t par0 = 0.930081; 
  Double_t par1 = 0.996683; 
  Double_t par2 = 3.54079e-05; 
  Double_t par3 = 0.0460187; 
  Double_t par4 = -23.2461; 

  if (ET > 200) ET =200;   		  
  if (             ET <    5 ) return         1.;  
  if (  5 <= ET && ET <   10 ) return         par0;  
  if ( 10 <= ET && ET <= 200 ) return         ( par1  + ET*par2)*(1-par3*exp(ET/par4));
 						  
}    


Double_t EnergyCorrection::applyScCorrectionsE_EE_electrons(Double_t E){      
   				 	  
  Double_t par0 = 400;               
  Double_t par1 = 0.982475; 
  Double_t par2 = 4.95413e-05; 
  Double_t par3 = 0.16886; 
  Double_t par4 = -30.1517; 
   				 	  
  if (E > par0) E =par0;   		  
  if (             E <   0 ) return         1.;  
  if (  0 <= E && E <= par0 ) return         (par1 + E*par2 )*(1- par3*exp(E/par4 ));
 						  
}         




double EnergyCorrection::f5x5( double iEta ) {
  if ( iEta < 40.2198 ) return 1;
  return 1 - 3.03103e-6*(iEta - 40.2198)*(iEta - 40.2198);
};

float EnergyCorrection::getEtaCorrectionBarrel(float eta){
  return 1.0/f5x5((int)(TMath::Abs(eta)*(5/0.087)));
};

bool EnergyCorrection::isInPhiCracks(double phi, double eta){


  // tranform radiants [-pi,pi] in degrees [0,360]
  phi = (phi+TMath::Pi()) *180/TMath::Pi();

  // each supermodule is 20 degrees wide
  Double_t moduleWidth = 20;

  // the first module is centered at phi=0, so the first cracks are +10 and -10
  Double_t phi0 = 10.;

  // set a fiducial cut around the crack of +2-2 degrees
  Double_t fiducialCut = 2.;

  bool OK = false;
  if (fabs(eta)<1.44){
    for (Int_t i = 0 ; i < 18; ++i){
      if ((phi0 + moduleWidth*i -fiducialCut) <= phi && phi <= (phi0 + moduleWidth*i + fiducialCut)) OK = true;
      //        cout << " PHI " << (phi0 + moduleWidth*i -fiducialCut) << " " << phi << " " <<  (phi0 + moduleWidth*i + fiducialCut)  << " " << OK << endl ;
    }
  }

  //  cout << "is in phi crack ? " << OK << endl;
  return OK;
};

double EnergyCorrection::get_correctedenergy(TreeReader *fTR, int ind,  int mode){
  if (forphotons) return getPho_correctedenergy(fTR, ind, mode);
  else return getEl_correctedenergy(fTR, ind, mode);
};

double EnergyCorrection::getPho_correctedenergy(TreeReader *fTR, int pi, int mode){

  if (mode!=0 && mode!=5 && mode!=6) {std::cout << "wrong call mode" << std::endl; return -999;}

  if (mode==0) return fTR->PhoEnergy[pi];

  int si=fTR->PhotSCindex[pi];

  if (si<0) {
    std::cout << "Calling photon energy correction without SC matching! Returning default photon energy." << std::endl;
    return fTR->PhoEnergy[pi];
  }

  float sc_eta = fTR->SCEta[si];
  float sc_phi = fTR->SCPhi[si];
  float sc_brem = fTR->SCBrem[si];

  bool isbarrel = fTR->PhoisEB[pi];
  bool isendcap = !isbarrel;

  bool highr9;
  double minr9;
  if (isbarrel) minr9=0.94; else minr9=0.95;
  if (fTR->PhoR9[pi]>minr9) highr9=true; else highr9=false;
  bool lowr9=!highr9;

  float energy=0;

  if (lowr9 && isbarrel) energy = fTR->SCRaw[si]*getEtaCorrectionBarrel(sc_eta);
  if (lowr9 && isendcap) energy = fTR->SCRaw[si]+fTR->SCPre[si];
  if (highr9 && isbarrel) energy = fTR->PhoE5x5[pi]; // e5x5 is already ceta corrected
  if (highr9 && isendcap) energy = fTR->PhoE5x5[pi]+fTR->SCPre[si];

  if (( TMath::Abs(sc_eta) < 1.44 ) && lowr9 && (energy/cosh(sc_eta)>10)){
    float et = energy/cosh(sc_eta);
    if (mode==5 || mode==6){
      et /= applyScCorrectionsBrEta_photons(sc_eta,sc_brem);
      et /= applyScCorrectionsET_EB_photons(et);
    }
    energy = et*cosh(sc_eta);
  }

  if ( ( TMath::Abs(sc_eta) < 2.5 ) && ( TMath::Abs(sc_eta) > 1.56 ) && lowr9 && (energy/cosh(sc_eta)>10)){
    float et = energy/cosh(sc_eta);
    if (mode==5 || mode==6){
      et /= applyScCorrectionsBrEta_photons(sc_eta,sc_brem);
      //      et /= applyScCorrectionsET_EE_photons(et);
      et /= applyScCorrectionsE_EE_photons(et*cosh(sc_eta));
    }
    energy = et*cosh(sc_eta);
  }
  
  if (mode==5 || mode==6){
    if (lowr9 && isbarrel) energy *= fTR->SCcrackcorr[si];
    if (highr9 && isbarrel) energy *= fTR->SCcrackcorrseedfactor[si];
  }

  if (mode==5){
    if (lowr9 && isbarrel) energy *= fTR->SClocalcorr[si];
    if (highr9 && isbarrel) energy *= fTR->SClocalcorrseedfactor[si];
  }

  return energy;
};

double EnergyCorrection::getEl_correctedenergy(TreeReader *fTR, int ei, int mode){

  if (mode!=0 && mode!=15 && mode!=16 && mode!=17 && mode!=20) {std::cout << "wrong call mode" << std::endl; return -999;}

  if (mode==0) return fTR->ElE[ei];

  int si=fTR->ElSCindex[ei];

  if (si<0) {
    std::cout << "Calling electron energy correction without SC matching! Returning default electron energy." << std::endl;
    return fTR->ElE[ei];
  }

  float sc_eta = fTR->SCEta[si];
  float sc_phi = fTR->SCPhi[si];
  float sc_brem = fTR->SCBrem[si];

  bool isbarrel = (fabs(sc_eta)<1.4442);
  bool isendcap = !isbarrel;

  float energy=0;

  if (isbarrel) energy = fTR->SCRaw[si]*getEtaCorrectionBarrel(sc_eta);
  if (isendcap) energy = fTR->SCRaw[si]+fTR->SCPre[si];

  if (mode==20) return fTR->SCEnergy[si];

  if (( TMath::Abs(sc_eta) < 1.4442 ) && (energy/cosh(sc_eta)>5)){
    float et = energy/cosh(sc_eta);

      et /= applyScCorrectionsBrEta_electrons(sc_eta,sc_brem);
      et /= applyScCorrectionsET_EB_electrons(et);

    energy = et*cosh(sc_eta);
  }

  if ( ( TMath::Abs(sc_eta) < 2.5) && ( TMath::Abs(sc_eta) > 1.56 ) && (energy/cosh(sc_eta)>5)){
    float et = energy/cosh(sc_eta);

      et /= applyScCorrectionsBrEta_electrons(sc_eta,sc_brem);
      // et /= applyScCorrectionsET_EE_electrons(et);
      et /= applyScCorrectionsE_EE_electrons(et*cosh(sc_eta));

    energy = et*cosh(sc_eta);
  }
  


  if (mode==15) if (isbarrel) energy *= fTR->SCcrackcorr[si];
  if (mode==16) if (isbarrel) energy *= fTR->SCcrackcorr[si]*fTR->SClocalcorr[si];
  if (mode==17) if (isbarrel) energy *= fTR->SCcrackcorrseed[si];

  return energy;
};



