#include "EnergyCorrection.hh"
#include "base/TreeReader.hh"
#include <assert.h>

EnergyCorrection::EnergyCorrection(TString tuning)
{
  
  {
    Double_t       leftEtatemp  [nBinsEta]   = { 0.02, 0.25, 0.46, 0.81, 0.91, 1.01, 1.16, 1.22, 1.33,           etaCrackMax,  1.653,  1.8, 2.0, 2.2, 2.3, 2.4 };
    Double_t       rightEtatemp [nBinsEta]   = { 0.25, 0.42, 0.77, 0.91, 1.01, 1.13, 1.22, 1.33, etaCrackMin,    1.653,        1.8  ,  2.0, 2.2, 2.3, 2.4, 2.5 };
    
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
  xcorr[7]=1.00595;
  xcorr[8]=1.00595;
  xcorr[9]=0.966651;
  xcorr[10]=0.97381;
  xcorr[11]=0.976516;
  xcorr[12]=0.983254;
  xcorr[13]=0.98502;
  xcorr[14]=0.98502;
  xcorr[15]=0.978472;

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

  par0[6] = 0.134234 ;
  par1[6] = 1.78462 ;
  par2[6] = -0.0536956 ;
  par3[6] = 0.987307 ;
  par4[6] = 2.54916 ;

  par0[7] = 0.111337 ;
  par1[7] = 1.33644 ;
  par2[7] = -0.0445354 ;
  par3[7] = 0.98227 ;
  par4[7] = 2.86281 ;

  par0[8] = 0.0265652 ;
  par1[8] = 1.48372 ;
  par2[8] = -0.0921363 ;
  par3[8] = 0.963707 ;
  par4[8] = 4.81597 ;

  par0[9] = 6.596 ;
  par1[9] = 4851.13 ;
  par2[9] = -2.63843 ;
  par3[9] = 0.970162 ;
  par4[9] = 1.00567 ;

  par0[10] = 2406.05 ;
  par1[10] = 855573 ;
  par2[10] = -1.91938 ;
  par3[10] = 0.980931 ;
  par4[10] = -0.0740673 ;

  par0[11] = 0.829199 ;
  par1[11] = 3.22717 ;
  par2[11] = -0.0485123 ;
  par3[11] = 0.957792 ;
  par4[11] = 2.01028 ;

  par0[12] = 0.275226 ;
  par1[12] = 2.20687 ;
  par2[12] = -0.110091 ;
  par3[12] = 0.93922 ;
  par4[12] = 2.69956 ;

  par0[13] = 0.169408 ;
  par1[13] = 1.53546 ;
  par2[13] = -0.0677257 ;
  par3[13] = 0.860703 ;
  par4[13] = 4.63723 ;

  par0[14] = 0.0780066 ;
  par1[14] = 1.39628 ;
  par2[14] = -0.0312043 ;
  par3[14] = 0.826576 ;
  par4[14] = 6.71971 ;

  par0[15] = 0.214497 ;
  par1[15] = 1.6747 ;
  par2[15] = -0.0858208 ;
  par3[15] = 0.893635 ;
  par4[15] = 3.7822 ;


  
  } // end for photons

  else { // for electrons

  xcorr[0]=1.00227;
  xcorr[1]=1.00252;
  xcorr[2]=1.00225;
  xcorr[3]=1.00159;
  xcorr[4]=0.999475;
  xcorr[5]=0.997203;
  xcorr[6]=0.99628;
  xcorr[7]=0.993651;
  xcorr[8]=0.993058;
  xcorr[9]=0.971262;
  xcorr[10]=0.975922;
  xcorr[11]=0.979087;
  xcorr[12]=0.98495;
  xcorr[13]=0.98781;
  xcorr[14]=0.989546;
  xcorr[15]=0.989638;


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

  par0[6] = 1.00565;
  par1[6] = -0.00856895;
  par2[6] = 0 ;

  par0[7] = 1.0052;
  par1[7] = -0.00886579;
  par2[7] = 0 ;

  par0[8] = 1.00325;
  par1[8] = -0.00849542;
  par2[8] = 0 ;

  par0[9] = 0.972798;
  par1[9] = -0.000771577;
  par2[9] = -0.00276696;

  par0[10] = 0.981672;
  par1[10] = -0.00202028;
  par2[10] = -0.00471028;

  par0[11] = 0.98251;
  par1[11] = 0.00441308;
  par2[11] = -0.00809139;

  par0[12] = 0.986123;
  par1[12] = 0.00832913;
  par2[12] = -0.00944584;

  par0[13] = 0.990124;
  par1[13] = 0.00742879;
  par2[13] = -0.00960462;

  par0[14] = 0.990187;
  par1[14] = 0.0094608;
  par2[14] = -0.010172;

  par0[15] = 0.99372;
  par1[15] = 0.00560406;
  par2[15] = -0.00943169;


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

Double_t EnergyCorrection::applyScCorrectionsETETA_working_photons(Double_t ET, Double_t eta){  
    							   
  // protect against high ET				   
								   
  if (ET > 200 ) ET =200;    				   
  if (fabs(eta) <1.33  || ET < 10 ) return 1;              

  else if (1.33 <= fabs(eta) && fabs(eta) < 1.44 ) { 
    if ( 10 <= ET && ET <  50 ) return     0.984352 + 0.000316509 *ET ;
    if ( 50 <= ET             ) return         1; 
  }                                               
  else if (1.56 <= fabs(eta) && fabs(eta) < 1.653 )        return     0.992191 + 5.85881e-05 *ET ;
  else if (1.653 <= fabs(eta) && fabs(eta) < 1.8 )        return     0.993616 + 4.58595e-05 *ET ;
  else if (1.8 <= fabs(eta) && fabs(eta) < 2 )        return     0.995683 + 3.29924e-05 *ET ;
  else if (2 <= fabs(eta) && fabs(eta) < 2.2 )        return     0.993543 + 5.77925e-05 *ET ;
  else if (2.2 <= fabs(eta) && fabs(eta) < 2.3 )        return     0.995799 + 2.5475e-05 *ET ;
  else if (2.3 <= fabs(eta) && fabs(eta) < 2.4 ) return 1;
  else if (2.4 <= fabs(eta) && fabs(eta) < 2.5 ) return 1;
  return 1;
}                                           

Double_t EnergyCorrection::applyScCorrectionsETETA_photons(Double_t ET, Double_t eta){ 

  // std::cout << "true eta: " << eta << std::endl;

  // eta = 0																						   
  if (TMath::Abs(eta)  <  leftEta[0]            ) { eta = 0.02 ; }	

  // outside acceptance

  if (TMath::Abs(eta)  >=  rightEta[nBinsEta-1] ) { eta = 2.49; if (DBG) std::cout << " WARNING [applyScCorrections]: TMath::Abs(eta)  >=  rightEta[nBinsEta-1] " << std::endl;}  
  																								   
 
  bool needinterpolation=true;
                               
  for (Int_t iEta = 0; iEta < nBinsEta; ++iEta){
    if ( leftEta[iEta] <= TMath::Abs(eta) && TMath::Abs(eta) <rightEta[iEta] ){
       needinterpolation=false;
    }	
  }													       										
  if (isInEBEtaCracks(eta) && !needinterpolation) {
    std::cout << "ERROR: not interpolating in EB eta crack" << std::endl;
    std::cout << eta << std::endl;
  };

  if (!needinterpolation) { // no need to interpolate
    // std::cout << "calling not interpolated with " << tmpEta << std::endl;
    return applyScCorrectionsETETA_working_photons(ET,eta);
  } 
  else { // need to interpolate    
   
    // std::cout << "calling interpolation" << std::endl;

    Double_t tmpLeftEta, tmpRightEta;
    
    for (Int_t iEta = 0; iEta < nBinsEta-1; ++iEta){								       								         if (rightEta[iEta] <= TMath::Abs(eta) && TMath::Abs(eta) <leftEta[iEta+1]) {							      
	tmpLeftEta=leftEta[iEta];         
	tmpRightEta=leftEta[iEta+1];
      }																						      }																						         
    // std::cout << "calling interpolation with "<< eta << " " << tmpLeftEta << " " << tmpRightEta << std::endl;
    return (applyScCorrectionsETETA_working_photons(ET,tmpLeftEta)+applyScCorrectionsETETA_working_photons(ET,tmpRightEta))/2;
      
  }
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



 

double EnergyCorrection::applyScCorrectionsETETA_working_electrons(Double_t ET, Double_t eta){ 

  Double_t par0 = 0;
  Double_t par1 = 0;
  Double_t par2 = 0;
  Double_t par3 = 0;
  Double_t par4 = 0;
  
  // std::cout << ET << std::endl;
  // std::cout << eta << std::endl;

  

  if (0.02 <= fabs(eta) && fabs(eta) <0.25) {
    par0 = 0.977764; 
    par1 = 1.00082; 
    par2 = -5.57895e-06; 
    par3 = 0.0189922; 
    par4 = -13.7688; 
  } 

  if (0.25 <= fabs(eta) && fabs(eta) <0.42) {
    par0 = 0.971133; 
    par1 = 1.00034; 
    par2 = -7.22258e-07; 
    par3 = 0.0149721; 
    par4 = -19.9583; 
  }                                             

  if (0.46 <= fabs(eta) && fabs(eta) <0.77) {
    par0 = 0.975339; 
    par1 = 1.0003; 
    par2 = -8.28272e-07; 
    par3 = 0.0112622; 
    par4 = -19.8705; 
  }                                           
   				 	  
  if (0.81 <= fabs(eta) && fabs(eta) <0.91) {
    par0 = 0.958131; 
    par1 = 0.999076; 
    par2 = 7.84516e-06; 
    par3 = 0.0549452; 
    par4 = -11.6424; 
  }                                           

   				 	  
  if (0.91 <= fabs(eta) && fabs(eta) <1.01) {
    par0 = 0.961929; 
    par1 = 0.996587; 
    par2 = 2.68649e-05; 
    par3 = 0.033021; 
    par4 = -16.3223; 
  }                                           

   				 	  
  if (1.01 <= fabs(eta) && fabs(eta) <1.13) {
    par0 = 0.917478; 
    par1 = 0.995432; 
    par2 = 3.42621e-05; 
    par3 = 0.0794779; 
    par4 = -16.9804; 
  }                                           

   				 	  
  if (1.16 <= fabs(eta) && fabs(eta) <1.22) {
    par0 = 0.912609; 
    par1 = 0.993364; 
    par2 = 5.12739e-05; 
    par3 = 0.0819801; 
    par4 = -23.3512; 
  }                                           

   				 	  
  if (1.22 <= fabs(eta) && fabs(eta) <1.33) {
    par0 = 0.864972; 
    par1 = 0.994976; 
    par2 = 4.13714e-05; 
    par3 = 0.122349; 
    par4 = -20.8449; 
  }                                           

   				 	  
  if (1.33 <= fabs(eta) && fabs(eta) <etaCrackMin) {
    par0 = 0.864351; 
    par1 = 0.991656; 
    par2 = 6.2759e-05; 
    par3 = 0.172499; 
    par4 = -16.6057; 
  }                                           

   				 	  
  if (etaCrackMax <= fabs(eta) && fabs(eta) <1.653) {
    par0 = 0.857992; 
    par1 = 0.98803; 
    par2 = 9.8499e-05; 
    par3 = 0.142157; 
    par4 = -17.6396; 

  }                                           
                					        
   				 	  
  if (1.653 <= fabs(eta) && fabs(eta) <1.8) {
    par0 = 0.877181; 
    par1 = 0.9858; 
    par2 = 0.000107183; 
    par3 = 0.182178; 
    par4 = -13.1828; 

  }                                           
                					        
   				 	  
  if (1.8 <= fabs(eta) && fabs(eta) <2) {
    par0 = 0.902391; 
    par1 = 0.993132; 
    par2 = 6.34629e-05; 
    par3 = 0.0868843; 
    par4 = -22.0789; 

  }                                           
                					        
   				 	  
  if (2 <= fabs(eta) && fabs(eta) <2.2) {
    par0 = 0.926399; 
    par1 = 0.996113; 
    par2 = 3.72639e-05; 
    par3 = 0.0512621; 
    par4 = -21.0456; 

  }                                           
                					        
   				 	  
  if (2.2 <= fabs(eta) && fabs(eta) <2.3) {
    par0 = 0.9419; 
    par1 = 0.997974; 
    par2 = 2.11951e-05; 
    par3 = 0.0436488; 
    par4 = -16.9869; 

  }                                           
                					        
   				 	  
  if (2.3 <= fabs(eta) && fabs(eta) <2.4) {
    par0 = 0.973133; 
    par1 = 1.00242; 
    par2 = -7.26518e-06; 
    par3 = 0.0241854; 
    par4 = -29.3233; 

  }                                           
                					        
   				 	  
  if (2.4 <= fabs(eta) && fabs(eta) <2.5) {
    par0 = 0.984542; 
    par1 = 1.00329; 
    par2 = -8.84497e-06; 
    par3 = 0.019583; 
    par4 = -38.8689; 

  }                                           


  if (par4 == 0) {
    std::cout << "WARNING YOU'RE ASKING FOR CORRECTIONS OUTSIDE ETA BOUNDARIES " << std::endl;
    std::cout << "input eta: " << eta << std::endl;
    std::cout << "input ET: " << ET << std::endl;
    return -999;
  }     
           					        
  if (ET > 200) ET =200;   		  
  if (             ET <    5 ) {
    std::cout << "WARNING YOU'RE ASKING FOR CORRECTIONS OUTSIDE ET LOW BOUNDARY " << std::endl;
    std::cout << "input eta: " << eta << std::endl;
    std::cout << "input ET: " << ET << std::endl;
    return         -999;  
  }
  if (  5 <= ET && ET <   10 ) return         par0 ;  
  if ( 10 <= ET && ET <= 200 ) return         (par1  + ET*par2)*(1- par3*exp(ET/ par4));



  return -999;

};


double EnergyCorrection::applyScCorrectionsETETA_electrons(Double_t ET, Double_t eta){ 

  // std::cout << "true eta: " << eta << std::endl;

  // eta = 0																						   
  if (TMath::Abs(eta)  <  leftEta[0]            ) { eta = 0.02 ; }	

  // outside acceptance

  if (TMath::Abs(eta)  >=  rightEta[nBinsEta-1] ) { eta = 2.49; if (DBG) std::cout << " WARNING [applyScCorrections]: TMath::Abs(eta)  >=  rightEta[nBinsEta-1] " << std::endl;}  
  																								   
 
  bool needinterpolation=true;
                               
  for (Int_t iEta = 0; iEta < nBinsEta; ++iEta){
    if ( leftEta[iEta] <= TMath::Abs(eta) && TMath::Abs(eta) <rightEta[iEta] ){
       needinterpolation=false;
    }	
  }													       										
  if (isInEBEtaCracks(eta) && !needinterpolation) {
    std::cout << "ERROR: not interpolating in EB eta crack" << std::endl;
    std::cout << eta << std::endl;
  };

  if (!needinterpolation) { // no need to interpolate
    // std::cout << "calling not interpolated with " << tmpEta << std::endl;
    return applyScCorrectionsETETA_working_electrons(ET,eta);
  } 
  else { // need to interpolate    
   
    // std::cout << "calling interpolation" << std::endl;

    Double_t tmpLeftEta, tmpRightEta;
    
    for (Int_t iEta = 0; iEta < nBinsEta-1; ++iEta){								       								         if (rightEta[iEta] <= TMath::Abs(eta) && TMath::Abs(eta) <leftEta[iEta+1]) {							      
	tmpLeftEta=leftEta[iEta];         
	tmpRightEta=leftEta[iEta+1];
      }																						      }																						         
    // std::cout << "calling interpolation with "<< eta << " " << tmpLeftEta << " " << tmpRightEta << std::endl;
    return (applyScCorrectionsETETA_working_electrons(ET,tmpLeftEta)+applyScCorrectionsETETA_working_electrons(ET,tmpRightEta))/2;
      
  }
};


double EnergyCorrection::f5x5( double iEta ) {
  if ( iEta < 40.2198 ) return 1;
  return 1 - 3.03103e-6*(iEta - 40.2198)*(iEta - 40.2198);
};

float EnergyCorrection::getEtaCorrectionBarrel(float eta){
  return 1.0/f5x5((int)(TMath::Abs(eta)*(5/0.087)));
};

bool EnergyCorrection::isInEBEtaCracks(double eta){

  if (TMath::Abs(eta)>=1.4442) return false;

  if (0.02<= TMath::Abs(eta) && TMath::Abs(eta) < 0.25) return false;
  if (0.25<= TMath::Abs(eta) && TMath::Abs(eta) < 0.42) return false; 
  if (0.46<= TMath::Abs(eta) && TMath::Abs(eta) < 0.77) return false; 
  if (0.81<= TMath::Abs(eta) && TMath::Abs(eta) < 0.91) return false; 
  if (0.91<= TMath::Abs(eta) && TMath::Abs(eta) < 1.01) return false; 
  if (1.01<= TMath::Abs(eta) && TMath::Abs(eta) < 1.13) return false; 
  if (1.16<= TMath::Abs(eta) && TMath::Abs(eta) < 1.22) return false; 
  if (1.22<= TMath::Abs(eta) && TMath::Abs(eta) < 1.33) return false; 
  if (1.33<= TMath::Abs(eta) && TMath::Abs(eta) < 1.4442) return false; 

  return true;

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
  if (fabs(eta)<1.4442){
    for (Int_t i = 0 ; i < 18; ++i){
      if ((phi0 + moduleWidth*i -fiducialCut) <= phi && phi <= (phi0 + moduleWidth*i + fiducialCut)) OK = true;
      //	std::cout << (int)OK << std::endl;
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

  if (( TMath::Abs(sc_eta) < 1.4442 ) && lowr9 && (energy/cosh(sc_eta)>10)){
    float et = energy/cosh(sc_eta);
    if (mode==5 || mode==6){
      et /= applyScCorrectionsBrEta_photons(sc_eta,sc_brem);
      //et /= applyScCorrectionsET_EB_photons(et);
      et /= applyScCorrectionsETETA_photons(et,sc_eta);
    }
    energy = et*cosh(sc_eta);
  }

  if ( ( TMath::Abs(sc_eta) < 2.5 ) && ( TMath::Abs(sc_eta) > 1.56 ) && lowr9 && (energy/cosh(sc_eta)>10)){
    float et = energy/cosh(sc_eta);
    if (mode==5 || mode==6){
      et /= applyScCorrectionsBrEta_photons(sc_eta,sc_brem);
      //      et /= applyScCorrectionsET_EE_photons(et);
      // et /= applyScCorrectionsE_EE_photons(et*cosh(sc_eta));
      et /= applyScCorrectionsETETA_photons(et,sc_eta);
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

  if (mode != 999 && mode!=0 && mode!=15 && mode!=16 && mode!=17 && mode!=18 && mode!=20 && mode!=21) {std::cout << "wrong call mode" << std::endl; return -999;}

  if (mode==0) return fTR->ElE[ei];

  if (mode==999) return fTR->ElGenE[ei];

  int si=fTR->ElSCindex[ei];

  if (si<0) {
    std::cout << "Calling electron energy correction without SC matching! Returning default electron energy." << std::endl;
    return fTR->ElE[ei];
  }

  float sc_eta = fTR->SCEta[si];
  float sc_phi = fTR->SCPhi[si];
  float sc_brem = fTR->SCBrem[si];

  bool isbarrel = (fabs(sc_eta)<1.4442);


  float energy=0;

  if (isbarrel) energy = fTR->SCRaw[si]*getEtaCorrectionBarrel(sc_eta);
  if (fabs(sc_eta)>1.56) energy = fTR->SCRaw[si]+fTR->SCPre[si];

  if (mode==21) return energy;

  if (mode==20) return fTR->SCEnergy[si];

  if (( TMath::Abs(sc_eta) < 1.4442 ) && (energy/cosh(sc_eta)>5)){
    float et = energy/cosh(sc_eta);

      et /= applyScCorrectionsBrEta_electrons(sc_eta,sc_brem);
      //  et /= applyScCorrectionsET_EB_electrons(et);
      et /= applyScCorrectionsETETA_electrons(et,sc_eta);

    energy = et*cosh(sc_eta);
  }

  if ( ( TMath::Abs(sc_eta) < 2.5) && ( TMath::Abs(sc_eta) > 1.56 ) && (energy/cosh(sc_eta)>5)){
    float et = energy/cosh(sc_eta);

      et /= applyScCorrectionsBrEta_electrons(sc_eta,sc_brem);
      // et /= applyScCorrectionsET_EE_electrons(et);
      // et /= applyScCorrectionsE_EE_electrons(et*cosh(sc_eta));
      et /= applyScCorrectionsETETA_electrons(et,sc_eta);

    energy = et*cosh(sc_eta);
  }
  



  if (mode==15) if (isbarrel) energy *= fTR->SCcrackcorr[si];
  if (mode==16) if (isbarrel) energy *= fTR->SCcrackcorr[si]*fTR->SClocalcorr[si];
  if (mode==17) if (isbarrel) energy *= fTR->SCcrackcorrseed[si];
  if (mode==18) ;

  /*
  if (mode==17){
    std::cout << "eta " << sc_eta << std::endl;
    std::cout << "gen " << fTR->ElGenE[ei] << std::endl;
    std::cout << "17 " << energy << std::endl;
    std::cout << "20 " << fTR->SCEnergy[si]  << std::endl;

  }
  */

  return energy;
};



