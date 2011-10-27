#include "EnergyCorrection.hh"
#include "base/TreeReader.hh"

EnergyCorrection::EnergyCorrection(TString tuning)
{

  if (tuning=="photons") forphotons=true;
  else if (tuning=="electrons") forphotons=false;
  else {
    std::cout << "Wrong initialization of EnergyCorrection!" << std::endl;
  }

  leftEta[0] =  0.02 ; 
  leftEta[1] =  0.25 ; 
  leftEta[2] =  0.46 ; 
  leftEta[3] =  0.81 ; 
  leftEta[4] =  0.91 ; 
  leftEta[5] =  1.01 ; 
  leftEta[6] =  1.16 ; 
  leftEta[7] =  1.56 ; 
  leftEta[8] =  1.653 ; 
  leftEta[9] =  1.8 ; 
  leftEta[10] =  2 ; 
  leftEta[11] =  2.2 ; 
  leftEta[12] =  2.3 ; 
  leftEta[13] =  2.4 ; 
  rightEta[0] =  0.25 ; 
  rightEta[1] =  0.42 ; 
  rightEta[2] =  0.77 ; 
  rightEta[3] =  0.91 ; 
  rightEta[4] =  1.01 ; 
  rightEta[5] =  1.13 ; 
  rightEta[6] =  1.44 ; 
  rightEta[7] =  1.653 ; 
  rightEta[8] =  1.8 ; 
  rightEta[9] =  2 ; 
  rightEta[10] =  2.2 ; 
  rightEta[11] =  2.3 ; 
  rightEta[12] =  2.4 ; 
  rightEta[13] =  2.5 ; 
  leftBr[0] =  0.8 ; 
  leftBr[1] =  1 ; 
  leftBr[2] =  1.2 ; 
  leftBr[3] =  1.4 ; 
  leftBr[4] =  1.6 ; 
  leftBr[5] =  1.8 ; 
  leftBr[6] =  2 ; 
  leftBr[7] =  2.2 ; 
  leftBr[8] =  2.4 ; 
  leftBr[9] =  2.6 ; 
  leftBr[10] =  2.8 ; 
  leftBr[11] =  3 ; 
  leftBr[12] =  3.2 ; 
  leftBr[13] =  3.4 ; 
  leftBr[14] =  3.6 ; 
  leftBr[15] =  3.8 ; 
  leftBr[16] =  4 ; 
  leftBr[17] =  5 ; 
  rightBr[0] =  1 ; 
  rightBr[1] =  1.2 ; 
  rightBr[2] =  1.4 ; 
  rightBr[3] =  1.6 ; 
  rightBr[4] =  1.8 ; 
  rightBr[5] =  2 ; 
  rightBr[6] =  2.2 ; 
  rightBr[7] =  2.4 ; 
  rightBr[8] =  2.6 ; 
  rightBr[9] =  2.8 ; 
  rightBr[10] =  3 ; 
  rightBr[11] =  3.2 ; 
  rightBr[12] =  3.4 ; 
  rightBr[13] =  3.6 ; 
  rightBr[14] =  3.8 ; 
  rightBr[15] =  4 ; 
  rightBr[16] =  5 ; 
  rightBr[17] =  10 ; 
  brbins[0] =  0.8 ; 
  brbins[1] =  1 ; 
  brbins[2] =  1 ; 
  brbins[3] =  1.2 ; 
  brbins[4] =  1.2 ; 
  brbins[5] =  1.4 ; 
  brbins[6] =  1.4 ; 
  brbins[7] =  1.6 ; 
  brbins[8] =  1.6 ; 
  brbins[9] =  1.8 ; 
  brbins[10] =  1.8 ; 
  brbins[11] =  2 ; 
  brbins[12] =  2 ; 
  brbins[13] =  2.2 ; 
  brbins[14] =  2.2 ; 
  brbins[15] =  2.4 ; 
  brbins[16] =  2.4 ; 
  brbins[17] =  2.6 ; 
  brbins[18] =  2.6 ; 
  brbins[19] =  2.8 ; 
  brbins[20] =  2.8 ; 
  brbins[21] =  3 ; 
  brbins[22] =  3 ; 
  brbins[23] =  3.2 ; 
  brbins[24] =  3.2 ; 
  brbins[25] =  3.4 ; 
  brbins[26] =  3.4 ; 
  brbins[27] =  3.6 ; 
  brbins[28] =  3.6 ; 
  brbins[29] =  3.8 ; 
  brbins[30] =  3.8 ; 
  brbins[31] =  4 ; 
  brbins[32] =  4 ; 
  brbins[33] =  5 ; 
  brbins[34] =  5 ; 
  brbins[35] =  10 ; 
  leftET[0] =  5 ; 
  leftET[1] =  10 ; 
  leftET[2] =  20 ; 
  leftET[3] =  30 ; 
  leftET[4] =  40 ; 
  leftET[5] =  50 ; 
  leftET[6] =  60 ; 
  leftET[7] =  80 ; 
  leftET[8] =  100 ; 
  leftET[9] =  120 ; 
  leftET[10] =  140 ; 
  leftET[11] =  160 ; 
  leftET[12] =  180 ; 
  leftET[13] =  200 ; 
  rightET[0] =  10 ; 
  rightET[1] =  20 ; 
  rightET[2] =  30 ; 
  rightET[3] =  40 ; 
  rightET[4] =  50 ; 
  rightET[5] =  60 ; 
  rightET[6] =  80 ; 
  rightET[7] =  100 ; 
  rightET[8] =  120 ; 
  rightET[9] =  140 ; 
  rightET[10] =  160 ; 
  rightET[11] =  180 ; 
  rightET[12] =  200 ; 
  rightET[13] =  250 ; 
  ETBins[0] =  5 ; 
  ETBins[1] =  10 ; 
  ETBins[2] =  10 ; 
  ETBins[3] =  20 ; 
  ETBins[4] =  20 ; 
  ETBins[5] =  30 ; 
  ETBins[6] =  30 ; 
  ETBins[7] =  40 ; 
  ETBins[8] =  40 ; 
  ETBins[9] =  50 ; 
  ETBins[10] =  50 ; 
  ETBins[11] =  60 ; 
  ETBins[12] =  60 ; 
  ETBins[13] =  80 ; 
  ETBins[14] =  80 ; 
  ETBins[15] =  100 ; 
  ETBins[16] =  100 ; 
  ETBins[17] =  120 ; 
  ETBins[18] =  120 ; 
  ETBins[19] =  140 ; 
  ETBins[20] =  140 ; 
  ETBins[21] =  160 ; 
  ETBins[22] =  160 ; 
  ETBins[23] =  180 ; 
  ETBins[24] =  180 ; 
  ETBins[25] =  200 ; 
  ETBins[26] =  200 ; 
  ETBins[27] =  250 ; 

  for (Int_t i = 0; i<nBinsEta; ++i){                                                   
    h_corr[i] = new TH1F(Form("h_corr_%d",i),Form("h_corr_%d",i),nBinsBr*2-1, brbins);  
  }                                                                                     

  h_CBET_EB    = new TH1F("h_CBET_EB"    ,"h_CBET_EB"    ,nBinsET*2-1, ETBins); 

  h_CBET_EE    = new TH1F("h_CBET_EE"    ,"h_CBET_EE"    ,nBinsET*2-1, ETBins); 

  if (forphotons){

    /* "legacy" binned corrections
  h_corr[0]->SetBinContent(1,1.00477 );
  h_corr[0]->SetBinContent(3,1.0005 );
  h_corr[0]->SetBinContent(5,1.00293 );
  h_corr[0]->SetBinContent(7,1.0046 );
  h_corr[0]->SetBinContent(9,1.00419 );
  h_corr[0]->SetBinContent(11,1.00477 );
  h_corr[0]->SetBinContent(13,1.00384 );
  h_corr[0]->SetBinContent(15,1.00448 );
  h_corr[0]->SetBinContent(17,1.00402 );
  h_corr[0]->SetBinContent(19,1.0051 );
  h_corr[0]->SetBinContent(21,1.00306 );
  h_corr[0]->SetBinContent(23,1.00324 );
  h_corr[0]->SetBinContent(25,1.00344 );
  h_corr[0]->SetBinContent(27,1.00536 );
  h_corr[0]->SetBinContent(29,1.0025 );
  h_corr[0]->SetBinContent(31,1.0021 );
  h_corr[0]->SetBinContent(33,1.00103 );
  h_corr[0]->SetBinContent(35,0.99688 );

  h_corr[1]->SetBinContent(1,1.00684 );
  h_corr[1]->SetBinContent(3,1.00073 );
  h_corr[1]->SetBinContent(5,1.00374 );
  h_corr[1]->SetBinContent(7,1.00425 );
  h_corr[1]->SetBinContent(9,1.0046 );
  h_corr[1]->SetBinContent(11,1.0038 );
  h_corr[1]->SetBinContent(13,1.00412 );
  h_corr[1]->SetBinContent(15,1.00263 );
  h_corr[1]->SetBinContent(17,1.00319 );
  h_corr[1]->SetBinContent(19,1.00138 );
  h_corr[1]->SetBinContent(21,1.00185 );
  h_corr[1]->SetBinContent(23,1.0011 );
  h_corr[1]->SetBinContent(25,0.999708 );
  h_corr[1]->SetBinContent(27,1.00156 );
  h_corr[1]->SetBinContent(29,1.00507 );
  h_corr[1]->SetBinContent(31,1.00024 );
  h_corr[1]->SetBinContent(33,0.999653 );
  h_corr[1]->SetBinContent(35,0.996048 );

  h_corr[2]->SetBinContent(1,1.00365 );
  h_corr[2]->SetBinContent(3,1.00016 );
  h_corr[2]->SetBinContent(5,1.00223 );
  h_corr[2]->SetBinContent(7,1.00357 );
  h_corr[2]->SetBinContent(9,1.00321 );
  h_corr[2]->SetBinContent(11,1.0029 );
  h_corr[2]->SetBinContent(13,1.0028 );
  h_corr[2]->SetBinContent(15,1.00168 );
  h_corr[2]->SetBinContent(17,1.00192 );
  h_corr[2]->SetBinContent(19,1.00164 );
  h_corr[2]->SetBinContent(21,1.00057 );
  h_corr[2]->SetBinContent(23,1.00231 );
  h_corr[2]->SetBinContent(25,0.99915 );
  h_corr[2]->SetBinContent(27,1.00004 );
  h_corr[2]->SetBinContent(29,0.999525 );
  h_corr[2]->SetBinContent(31,0.998512 );
  h_corr[2]->SetBinContent(33,0.997322 );
  h_corr[2]->SetBinContent(35,0.994791 );

  h_corr[3]->SetBinContent(1,1.00365 );
  h_corr[3]->SetBinContent(3,1.00091 );
  h_corr[3]->SetBinContent(5,1.00165 );
  h_corr[3]->SetBinContent(7,1.00144 );
  h_corr[3]->SetBinContent(9,1.00161 );
  h_corr[3]->SetBinContent(11,1.00033 );
  h_corr[3]->SetBinContent(13,0.998908 );
  h_corr[3]->SetBinContent(15,0.998933 );
  h_corr[3]->SetBinContent(17,0.998971 );
  h_corr[3]->SetBinContent(19,0.996114 );
  h_corr[3]->SetBinContent(21,0.998552 );
  h_corr[3]->SetBinContent(23,0.997014 );
  h_corr[3]->SetBinContent(25,0.996017 );
  h_corr[3]->SetBinContent(27,0.994346 );
  h_corr[3]->SetBinContent(29,0.997474 );
  h_corr[3]->SetBinContent(31,0.994367 );
  h_corr[3]->SetBinContent(33,0.989633 );
  h_corr[3]->SetBinContent(35,0.983556 );

  h_corr[4]->SetBinContent(1,1.00365 );
  h_corr[4]->SetBinContent(3,0.999573 );
  h_corr[4]->SetBinContent(5,0.999298 );
  h_corr[4]->SetBinContent(7,0.999856 );
  h_corr[4]->SetBinContent(9,0.999745 );
  h_corr[4]->SetBinContent(11,0.998014 );
  h_corr[4]->SetBinContent(13,0.998564 );
  h_corr[4]->SetBinContent(15,0.996768 );
  h_corr[4]->SetBinContent(17,0.995129 );
  h_corr[4]->SetBinContent(19,0.992972 );
  h_corr[4]->SetBinContent(21,0.99556 );
  h_corr[4]->SetBinContent(23,0.992829 );
  h_corr[4]->SetBinContent(25,0.990584 );
  h_corr[4]->SetBinContent(27,0.988897 );
  h_corr[4]->SetBinContent(29,0.987413 );
  h_corr[4]->SetBinContent(31,0.986436 );
  h_corr[4]->SetBinContent(33,0.987777 );
  h_corr[4]->SetBinContent(35,0.968688 );

  h_corr[5]->SetBinContent(1,1.00365 );
  h_corr[5]->SetBinContent(3,0.997488 );
  h_corr[5]->SetBinContent(5,0.997056 );
  h_corr[5]->SetBinContent(7,0.996287 );
  h_corr[5]->SetBinContent(9,0.997043 );
  h_corr[5]->SetBinContent(11,0.995985 );
  h_corr[5]->SetBinContent(13,0.993078 );
  h_corr[5]->SetBinContent(15,0.992365 );
  h_corr[5]->SetBinContent(17,0.98878 );
  h_corr[5]->SetBinContent(19,0.985867 );
  h_corr[5]->SetBinContent(21,0.988083 );
  h_corr[5]->SetBinContent(23,0.982983 );
  h_corr[5]->SetBinContent(25,0.984689 );
  h_corr[5]->SetBinContent(27,0.986268 );
  h_corr[5]->SetBinContent(29,0.977702 );
  h_corr[5]->SetBinContent(31,0.976403 );
  h_corr[5]->SetBinContent(33,0.970385 );
  h_corr[5]->SetBinContent(35,0.95356 );

  h_corr[6]->SetBinContent(1,1.00365 );
  h_corr[6]->SetBinContent(3,0.996276 );
  h_corr[6]->SetBinContent(5,0.995901 );
  h_corr[6]->SetBinContent(7,0.995306 );
  h_corr[6]->SetBinContent(9,0.994644 );
  h_corr[6]->SetBinContent(11,0.993144 );
  h_corr[6]->SetBinContent(13,0.991161 );
  h_corr[6]->SetBinContent(15,0.988541 );
  h_corr[6]->SetBinContent(17,0.986377 );
  h_corr[6]->SetBinContent(19,0.985698 );
  h_corr[6]->SetBinContent(21,0.983213 );
  h_corr[6]->SetBinContent(23,0.982471 );
  h_corr[6]->SetBinContent(25,0.978945 );
  h_corr[6]->SetBinContent(27,0.979868 );
  h_corr[6]->SetBinContent(29,0.973695 );
  h_corr[6]->SetBinContent(31,0.971409 );
  h_corr[6]->SetBinContent(33,0.967811 );
  h_corr[6]->SetBinContent(35,0.95356 );

  h_corr[7]->SetBinContent(1,0.966757 );
  h_corr[7]->SetBinContent(3,0.973499 );
  h_corr[7]->SetBinContent(5,0.969414 );
  h_corr[7]->SetBinContent(7,0.969698 );
  h_corr[7]->SetBinContent(9,0.96632 );
  h_corr[7]->SetBinContent(11,0.965057 );
  h_corr[7]->SetBinContent(13,0.962489 );
  h_corr[7]->SetBinContent(15,0.962289 );
  h_corr[7]->SetBinContent(17,0.95726 );
  h_corr[7]->SetBinContent(19,0.95717 );
  h_corr[7]->SetBinContent(21,0.950804 );
  h_corr[7]->SetBinContent(23,0.946758 );
  h_corr[7]->SetBinContent(25,0.943674 );
  h_corr[7]->SetBinContent(27,0.942974 );
  h_corr[7]->SetBinContent(29,0.939583 );
  h_corr[7]->SetBinContent(31,0.936795 );
  h_corr[7]->SetBinContent(33,0.91318 );
  h_corr[7]->SetBinContent(35,0.864583 );

  h_corr[8]->SetBinContent(1,0.973319 );
  h_corr[8]->SetBinContent(3,0.974327 );
  h_corr[8]->SetBinContent(5,0.97135 );
  h_corr[8]->SetBinContent(7,0.969466 );
  h_corr[8]->SetBinContent(9,0.96644 );
  h_corr[8]->SetBinContent(11,0.96251 );
  h_corr[8]->SetBinContent(13,0.960385 );
  h_corr[8]->SetBinContent(15,0.954011 );
  h_corr[8]->SetBinContent(17,0.945243 );
  h_corr[8]->SetBinContent(19,0.95272 );
  h_corr[8]->SetBinContent(21,0.942462 );
  h_corr[8]->SetBinContent(23,0.933115 );
  h_corr[8]->SetBinContent(25,0.91591 );
  h_corr[8]->SetBinContent(27,0.919095 );
  h_corr[8]->SetBinContent(29,0.8951 );
  h_corr[8]->SetBinContent(31,0.91149 );
  h_corr[8]->SetBinContent(33,0.872187 );
  h_corr[8]->SetBinContent(35,0.831379 );

  h_corr[9]->SetBinContent(1,0.97587 );
  h_corr[9]->SetBinContent(3,0.972743 );
  h_corr[9]->SetBinContent(5,0.971705 );
  h_corr[9]->SetBinContent(7,0.968167 );
  h_corr[9]->SetBinContent(9,0.968744 );
  h_corr[9]->SetBinContent(11,0.96537 );
  h_corr[9]->SetBinContent(13,0.958599 );
  h_corr[9]->SetBinContent(15,0.954942 );
  h_corr[9]->SetBinContent(17,0.946286 );
  h_corr[9]->SetBinContent(19,0.931168 );
  h_corr[9]->SetBinContent(21,0.929851 );
  h_corr[9]->SetBinContent(23,0.904241 );
  h_corr[9]->SetBinContent(25,0.920334 );
  h_corr[9]->SetBinContent(27,1.01417 );
  h_corr[9]->SetBinContent(29,0.894637 );
  h_corr[9]->SetBinContent(31,0.894975 );
  h_corr[9]->SetBinContent(33,0.856302 );
  h_corr[9]->SetBinContent(35,0.831379 );

  h_corr[10]->SetBinContent(1,0.980946 );
  h_corr[10]->SetBinContent(3,0.975896 );
  h_corr[10]->SetBinContent(5,0.97371 );
  h_corr[10]->SetBinContent(7,0.973454 );
  h_corr[10]->SetBinContent(9,0.971114 );
  h_corr[10]->SetBinContent(11,0.967695 );
  h_corr[10]->SetBinContent(13,0.965971 );
  h_corr[10]->SetBinContent(15,0.953794 );
  h_corr[10]->SetBinContent(17,0.95019 );
  h_corr[10]->SetBinContent(19,0.927728 );
  h_corr[10]->SetBinContent(21,0.928867 );
  h_corr[10]->SetBinContent(23,0.920211 );
  h_corr[10]->SetBinContent(25,0.920334 );
  h_corr[10]->SetBinContent(27,1.01417 );
  h_corr[10]->SetBinContent(29,0.894637 );
  h_corr[10]->SetBinContent(31,0.894975 );
  h_corr[10]->SetBinContent(33,0.856302 );
  h_corr[10]->SetBinContent(35,0.831379 );

  h_corr[11]->SetBinContent(1,0.985509 );
  h_corr[11]->SetBinContent(3,0.975446 );
  h_corr[11]->SetBinContent(5,0.976873 );
  h_corr[11]->SetBinContent(7,0.975401 );
  h_corr[11]->SetBinContent(9,0.974248 );
  h_corr[11]->SetBinContent(11,0.972033 );
  h_corr[11]->SetBinContent(13,0.968997 );
  h_corr[11]->SetBinContent(15,0.950584 );
  h_corr[11]->SetBinContent(17,0.95019 );
  h_corr[11]->SetBinContent(19,0.927728 );
  h_corr[11]->SetBinContent(21,0.928867 );
  h_corr[11]->SetBinContent(23,0.920211 );
  h_corr[11]->SetBinContent(25,0.920334 );
  h_corr[11]->SetBinContent(27,1.01417 );
  h_corr[11]->SetBinContent(29,0.894637 );
  h_corr[11]->SetBinContent(31,0.894975 );
  h_corr[11]->SetBinContent(33,0.856302 );
  h_corr[11]->SetBinContent(35,0.831379 );

  h_corr[12]->SetBinContent(1,0.985509 );
  h_corr[12]->SetBinContent(3,0.974643 );
  h_corr[12]->SetBinContent(5,0.975387 );
  h_corr[12]->SetBinContent(7,0.977001 );
  h_corr[12]->SetBinContent(9,0.974089 );
  h_corr[12]->SetBinContent(11,0.967086 );
  h_corr[12]->SetBinContent(13,0.964069 );
  h_corr[12]->SetBinContent(15,0.950584 );
  h_corr[12]->SetBinContent(17,0.95019 );
  h_corr[12]->SetBinContent(19,0.927728 );
  h_corr[12]->SetBinContent(21,0.928867 );
  h_corr[12]->SetBinContent(23,0.920211 );
  h_corr[12]->SetBinContent(25,0.920334 );
  h_corr[12]->SetBinContent(27,1.01417 );
  h_corr[12]->SetBinContent(29,0.894637 );
  h_corr[12]->SetBinContent(31,0.894975 );
  h_corr[12]->SetBinContent(33,0.856302 );
  h_corr[12]->SetBinContent(35,0.831379 );

  h_corr[13]->SetBinContent(1,0.979224 );
  h_corr[13]->SetBinContent(3,0.975708 );
  h_corr[13]->SetBinContent(5,0.973844 );
  h_corr[13]->SetBinContent(7,0.971448 );
  h_corr[13]->SetBinContent(9,0.975281 );
  h_corr[13]->SetBinContent(11,0.976356 );
  h_corr[13]->SetBinContent(13,0.964069 );
  h_corr[13]->SetBinContent(15,0.950584 );
  h_corr[13]->SetBinContent(17,0.95019 );
  h_corr[13]->SetBinContent(19,0.927728 );
  h_corr[13]->SetBinContent(21,0.928867 );
  h_corr[13]->SetBinContent(23,0.920211 );
  h_corr[13]->SetBinContent(25,0.920334 );
  h_corr[13]->SetBinContent(27,1.01417 );
  h_corr[13]->SetBinContent(29,0.894637 );
  h_corr[13]->SetBinContent(31,0.894975 );
  h_corr[13]->SetBinContent(33,0.856302 );
  h_corr[13]->SetBinContent(35,0.831379 );

  h_CBET_EB->SetBinContent(1, 0.996172); 
  h_CBET_EB->SetBinContent(3, 1.00346); 
  h_CBET_EB->SetBinContent(5, 1.00096); 
  h_CBET_EB->SetBinContent(7, 1.00067); 
  h_CBET_EB->SetBinContent(9, 1.00105); 
  h_CBET_EB->SetBinContent(11, 1.00053); 
  h_CBET_EB->SetBinContent(13, 1.00052); 
  h_CBET_EB->SetBinContent(15, 1.00023); 
  h_CBET_EB->SetBinContent(17, 0.999966); 
  h_CBET_EB->SetBinContent(19, 1.00009); 
  h_CBET_EB->SetBinContent(21, 0.999689); 
  h_CBET_EB->SetBinContent(23, 0.999704); 
  h_CBET_EB->SetBinContent(25, 0.999898); 
  h_CBET_EB->SetBinContent(27, 0.999869); 

  h_CBET_EE->SetBinContent(1, 0.962984); 
  h_CBET_EE->SetBinContent(3, 0.998347); 
  h_CBET_EE->SetBinContent(5, 0.995259); 
  h_CBET_EE->SetBinContent(7, 0.994463); 
  h_CBET_EE->SetBinContent(9, 0.994789); 
  h_CBET_EE->SetBinContent(11, 0.995445); 
  h_CBET_EE->SetBinContent(13, 0.996237); 
  h_CBET_EE->SetBinContent(15, 0.997266); 
  h_CBET_EE->SetBinContent(17, 0.998239); 
  h_CBET_EE->SetBinContent(19, 0.999248); 
  h_CBET_EE->SetBinContent(21, 0.999833); 
  h_CBET_EE->SetBinContent(23, 1.0003); 
  h_CBET_EE->SetBinContent(25, 1.00502); 
  h_CBET_EE->SetBinContent(27, 1.00295); 
    */

  xcorr[0]=1.00477;
  xcorr[1]=1.00684;
  xcorr[2]=1.00365;
  xcorr[3]=1.00365;
  xcorr[4]=1.00365;
  xcorr[5]=1.00365;
  xcorr[6]=1.00365;
  xcorr[7]=0.966757;
  xcorr[8]=0.973319;
  xcorr[9]=0.97587;
  xcorr[10]=0.980946;
  xcorr[11]=0.985509;
  xcorr[12]=0.985509;
  xcorr[13]=0.979224;

  fcorr[0] = new TF1("ftest_0","[0]*(1.-exp(-(x-[4])/[1]))*[2]*x + [3]",1,10); 
  fcorr[0]->SetParameters(-0.00319973, 2.2396, 0.00197239, 0.986858, 18.2498 );
  fcorr[1] = new TF1("ftest_1","[0]*(1.-exp(-(x-[4])/[1]))*[2]*x + [3]",1,10); 
  fcorr[1]->SetParameters(-0.00291167, 1.9014, 0.00135785, 0.985441, 16.7636 );
  fcorr[2] = new TF1("ftest_2","[0]*(1.-exp(-(x-[4])/[1]))*[2]*x + [3]",1,10); 
  fcorr[2]->SetParameters(0.000247919, 1.876, -0.0217462, 0.986375, 15.7761 );
  fcorr[3] = new TF1("ftest_3","[0]*(1.-exp(-(x-[4])/[1]))*[2]*x + [3]",1,10); 
  fcorr[3]->SetParameters(0.00375001, 1.35455, -0.00239691, 0.986979, 10.9435 );
  fcorr[4] = new TF1("ftest_4","[0]*(1.-exp(-(x-[4])/[1]))*[2]*x + [3]",1,10); 
  fcorr[4]->SetParameters(0.0120053, 1.35895, -0.00464863, 0.979589, 8.95397 );
  fcorr[5] = new TF1("ftest_5","[0]*(1.-exp(-(x-[4])/[1]))*[2]*x + [3]",1,10); 
  fcorr[5]->SetParameters(0.0830868, 1.51088, -0.0339125, 0.970732, 4.4903 );
  fcorr[6] = new TF1("ftest_6","[0]*(1.-exp(-(x-[4])/[1]))*[2]*x + [3]",1,10); 
  fcorr[6]->SetParameters(0.0939344, 1.40516, -0.0375749, 0.975119, 3.71635 );
  fcorr[7] = new TF1("ftest_7","[0]*(1.-exp(-(x-[4])/[1]))*[2]*x + [3]",1,10); 
  fcorr[7]->SetParameters(0.276827, 319.43, -3.60775, 0.973674, 0.463167 );
  fcorr[8] = new TF1("ftest_8","[0]*(1.-exp(-(x-[4])/[1]))*[2]*x + [3]",1,10); 
  fcorr[8]->SetParameters(1337.84, 485627, -2.19017, 0.976675, 0.689757 );
  fcorr[9] = new TF1("ftest_9","[0]*(1.-exp(-(x-[4])/[1]))*[2]*x + [3]",1,10); 
  fcorr[9]->SetParameters(0.311055, 2.72391, -0.124422, 0.942146, 2.56506 );
  fcorr[10] = new TF1("ftest_10","[0]*(1.-exp(-(x-[4])/[1]))*[2]*x + [3]",1,10); 
  fcorr[10]->SetParameters(0.310272, 2.69723, -0.124112, 0.940468, 2.70685 );
  fcorr[11] = new TF1("ftest_11","[0]*(1.-exp(-(x-[4])/[1]))*[2]*x + [3]",1,10); 
  fcorr[11]->SetParameters(0.231957, 1.7835, -0.0927829, 0.886792, 3.87139 );
  fcorr[12] = new TF1("ftest_12","[0]*(1.-exp(-(x-[4])/[1]))*[2]*x + [3]",1,10); 
  fcorr[12]->SetParameters(0.203507, 1.63421, -0.0810196, 0.873407, 4.17771 );
  fcorr[13] = new TF1("ftest_13","[0]*(1.-exp(-(x-[4])/[1]))*[2]*x + [3]",1,10); 
  fcorr[13]->SetParameters(2.33719, 2.17042, -0.0131261, 0.923385, 3.09554 );

  } // end for photons

  else { // for electrons

  xcorr[0]=1.00081;
  xcorr[1]=1.00271;
  xcorr[2]=1.0022;
  xcorr[3]=1.00199;
  xcorr[4]=0.999823;
  xcorr[5]=0.997862;
  xcorr[6]=0.995133;
  xcorr[7]=0.97329;
  xcorr[8]=0.97803;
  xcorr[9]=0.98094;
  xcorr[10]=0.985973;
  xcorr[11]=0.988175;
  xcorr[12]=0.99007;
  xcorr[13]=0.990097;

  fcorr[0] = new TF1("ftest_0","pol1",1,10); 
  fcorr[0]->SetParameters(1.00689, -0.00176103);
  fcorr[1] = new TF1("ftest_1","pol1",1,10); 
  fcorr[1]->SetParameters(1.00695, -0.00216774);
  fcorr[2] = new TF1("ftest_2","pol1",1,10); 
  fcorr[2]->SetParameters(1.00633, -0.00248825);
  fcorr[3] = new TF1("ftest_3","pol1",1,10); 
  fcorr[3]->SetParameters(1.00757, -0.00421121);
  fcorr[4] = new TF1("ftest_4","pol1",1,10); 
  fcorr[4]->SetParameters(1.0066, -0.00508913);
  fcorr[5] = new TF1("ftest_5","pol1",1,10); 
  fcorr[5]->SetParameters(1.00805, -0.00799981);
  fcorr[6] = new TF1("ftest_6","pol1",1,10); 
  fcorr[6]->SetParameters(1.00574, -0.00853869);
  fcorr[7] = new TF1("ftest_7","pol2",1,10); 
  fcorr[7]->SetParameters(0.974138, 0.000681559, -0.003205 );
  fcorr[8] = new TF1("ftest_8","pol2",1,10); 
  fcorr[8]->SetParameters(0.982062, 0.000505085, -0.00540531 );
  fcorr[9] = new TF1("ftest_9","pol2",1,10); 
  fcorr[9]->SetParameters(0.98278, 0.00707044, -0.00886468 );
  fcorr[10] = new TF1("ftest_10","pol2",1,10); 
  fcorr[10]->SetParameters(0.984711, 0.0110689, -0.0100418 );
  fcorr[11] = new TF1("ftest_11","pol2",1,10); 
  fcorr[11]->SetParameters(0.988272, 0.0102602, -0.0102262 );
  fcorr[12] = new TF1("ftest_12","pol2",1,10); 
  fcorr[12]->SetParameters(0.99132, 0.00962309, -0.0102865 );
  fcorr[13] = new TF1("ftest_13","pol2",1,10); 
  fcorr[13]->SetParameters(0.993064, 0.00748717, -0.0098858 );

  }

};

 Double_t EnergyCorrection::applyScCorrectionsET_EB_FIT(Double_t ET){      
   					  
   if (             ET <   5 ) return         1.;  
   if (  5 <= ET && ET <  10 ) return         0.996172;  
   if ( 10 <= ET && ET <  20 ) return         1.00346;  
   if ( 20 <= ET && ET < 140 ) return         1.00107 + -8.59261e-06*ET ;  
   if (140 <= ET             ) return         0.999828;  
 						  
 };
                					        

 Double_t EnergyCorrection::applyScCorrectionsET_EE_FIT(Double_t ET){      
   					  
   if (             ET <   5 ) return         1.;  
   if (  5 <= ET && ET <  10 ) return         0.962984;  
   if ( 10 <= ET && ET <  20 ) return         0.998347;  
   if ( 20 <= ET && ET <  30 ) return         0.995259;  
   if ( 30 <= ET && ET < 200 ) return         0.991964 + (6.48159e-05)*ET    + (-6.30917e-08)*ET*ET ;  
   if ( 200 <= ET            ) return         1.00295;
 						  
 };

Double_t EnergyCorrection::applyScCorrectionsET_EB_FIT_electrons(Double_t ET){      
                                          
   if (ET > 250) ET =250;                 
   if (             ET <    5 ) return         1.;  
   if (  5 <= ET && ET <   10 ) return         0.972175;  
   if ( 10 <= ET && ET <= 250 ) return         (0.999631 + ET*4.64281e-06)*(1-0.0163644*exp(ET/-16.2768));
                                                  
};

 Double_t EnergyCorrection::applyScCorrectionsET_EE_FIT_electrons(Double_t ET){      
                                          
   if (ET > 250) ET =250;                 
   if (             ET <    5 ) return         1.;  
   if (  5 <= ET && ET <   10 ) return         0.928851;  
   if ( 10 <= ET && ET <= 250 ) return         (0.993435 + ET*4.75233e-05)*(1-0.0521102*exp(ET/-18.176));
                                                  
 }                                   


EnergyCorrection::~EnergyCorrection() 
{
  for (Int_t i = 0; i<nBinsEta; ++i) delete h_corr[i];
  delete h_CBET_EB;
  delete h_CBET_EE;
  for (Int_t i = 0; i<nBinsEta; ++i) delete fcorr[i];
}

Double_t EnergyCorrection::applyScCorrectionsBrEta_FIT(Double_t eta, Double_t sigmaPhiSigmaEta){  
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
 	        if (sigmaPhiSigmaEta >= 1)  tmpInter = ( fcorr[iEta]->Eval(sigmaPhiSigmaEta) + 															          
 		         	    		         fcorr[iEta+1]->Eval(sigmaPhiSigmaEta) ) / 2. ;														          
 	        else tmpInter = (xcorr[iEta] + xcorr[iEta+1])/2.; 
        }																						         
      }																						         
      return tmpInter;																					         
    }  																							         
    if (sigmaPhiSigmaEta >= 1) return fcorr[tmpEta]->Eval(sigmaPhiSigmaEta);						       						       
    else return xcorr[tmpEta]; 
};

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

  if (mode!=0 && mode!=1 && mode!=2 && mode!=5 && mode!=6) {std::cout << "wrong call mode" << std::endl; return -999;}

  if (mode==0) return fTR->PhoEnergy[pi];

  if (mode==1 || mode==2) {
    std::cout << "Legacy call mode (binned corrections, not fit), not fully implemented here. Returning error flag." << std::endl;
    return -999;
  }

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

  if (( TMath::Abs(sc_eta) < 1.44 ) && lowr9 && (energy/cosh(sc_eta)>5)){
    float et = energy/cosh(sc_eta);
    if (mode==1 || mode==2){
      //    et /= applyScCorrectionsBrEta(sc_eta,sc_brem);
      //    et /= applyScCorrectionsET_EB(et);
    }
    else if (mode==5 || mode==6){
      et /= applyScCorrectionsBrEta_FIT(sc_eta,sc_brem);
      et /= applyScCorrectionsET_EB_FIT(et);
    }
    energy = et*cosh(sc_eta);
  }

  if ( ( TMath::Abs(sc_eta) < 2.5 ) && ( TMath::Abs(sc_eta) > 1.56 ) && lowr9 && (energy/cosh(sc_eta)>5)){
    float et = energy/cosh(sc_eta);
    if (mode==1 || mode==2){
      //      et /= applyScCorrectionsBrEta(sc_eta,sc_brem);
      //      et /= applyScCorrectionsET_EE(et);
    }
    else if (mode==5 || mode==6){
      et /= applyScCorrectionsBrEta_FIT(sc_eta,sc_brem);
      et /= applyScCorrectionsET_EE_FIT(et);
    }
    energy = et*cosh(sc_eta);
  }
  
  if (lowr9 && isbarrel) energy *= fTR->SCcrackcorr[si];
  if (highr9 && isbarrel) energy *= fTR->SCcrackcorrseedfactor[si];


  if (mode==1 || mode==5){
  if (lowr9 && isbarrel) energy *= fTR->SClocalcorr[si];
  if (highr9 && isbarrel) energy *= fTR->SClocalcorrseedfactor[si];
  }


  return energy;
};

double EnergyCorrection::getEl_correctedenergy(TreeReader *fTR, int ei, int mode){

  if (mode!=0 && mode!=15 && mode!=16 && mode!=20) {std::cout << "wrong call mode" << std::endl; return -999;}

  if (mode==0) return fTR->ElE[ei];

  int si=fTR->ElSCindex[ei];

  if (si<0) {
    std::cout << "Calling electron energy correction without SC matching! Returning default electron energy." << std::endl;
    return fTR->ElE[ei];
  }

  float sc_eta = fTR->SCEta[si];
  float sc_phi = fTR->SCPhi[si];
  float sc_brem = fTR->SCBrem[si];

  bool isbarrel = (fabs(sc_eta)<1.44);
  bool isendcap = !isbarrel;

  float energy=0;

  if (isbarrel) energy = fTR->SCRaw[si]*getEtaCorrectionBarrel(sc_eta);
  if (isendcap) energy = fTR->SCRaw[si]+fTR->SCPre[si];

  if (mode==20) return energy;

  if (( TMath::Abs(sc_eta) < 1.44 ) && (energy/cosh(sc_eta)>5)){
    float et = energy/cosh(sc_eta);
    if (mode==15 || mode==16){
      et /= applyScCorrectionsBrEta_FIT(sc_eta,sc_brem);
      et /= applyScCorrectionsET_EB_FIT_electrons(et);
    }
    energy = et*cosh(sc_eta);
  }

  if ( ( TMath::Abs(sc_eta) < 2.5) && ( TMath::Abs(sc_eta) > 1.56 ) && (energy/cosh(sc_eta)>5)){
    float et = energy/cosh(sc_eta);
    if (mode==15 || mode==16){
      et /= applyScCorrectionsBrEta_FIT(sc_eta,sc_brem);
      et /= applyScCorrectionsET_EE_FIT_electrons(et);
    }
    energy = et*cosh(sc_eta);
  }
  
  if (mode==15 || mode==16) if (isbarrel) energy *= fTR->SCcrackcorr[si];

  if (mode==15) if (isbarrel) energy *= fTR->SClocalcorr[si];


  return energy;
};



/*
// legacy code for binned corrections


Double_t PhoEnergyCorrection::applyScCorrectionsBrEta(Double_t eta, Double_t sigmaPhiSigmaEta){  
  
  Int_t tmpEta = -1;                                                                                             
  Int_t tmpBr  = -1;	

  // Extra protections										       
   
  if (TMath::Abs(eta)  <   leftEta[0]           ) { tmpEta = 0;          if (DBG) std::cout << " WARNING [applyScCorrections]: (TMath::Abs(eta)  <   leftEta[0]          " << std::endl;}
  if (TMath::Abs(eta)  >=  rightEta[nBinsEta-1] ) { tmpEta = nBinsEta-1; if (DBG) std::cout << " WARNING [applyScCorrections]: TMath::Abs(eta)  >=  rightEta[nBinsEta-1] " << std::endl;}

  if (sigmaPhiSigmaEta <  leftBr [0]            ) {tmpBr = 0;            if (DBG) std::cout << " WARNING [applyScCorrections]: sigmaPhiSigmaEta <  leftBr [0]            " << std::endl;}
  if (sigmaPhiSigmaEta >= rightBr[nBinsBr]      ) {tmpBr = nBinsBr  -1;  if (DBG) std::cout << " WARNING [applyScCorrections]: sigmaPhiSigmaEta >= rightBr[nBinsBr]      " << std::endl;}
  
  for (Int_t iEta = 0; iEta < nBinsEta; ++iEta){								       
    if ( leftEta[iEta] <= TMath::Abs(eta) && TMath::Abs(eta) <rightEta[iEta] ){				       
      tmpEta = iEta;											       
    }													       
  }													         

  for (Int_t iSigmaPhiSigmaEta = 0; iSigmaPhiSigmaEta < nBinsBr; ++iSigmaPhiSigmaEta){			       
    if (leftBr [iSigmaPhiSigmaEta]  <= sigmaPhiSigmaEta && sigmaPhiSigmaEta <rightBr [iSigmaPhiSigmaEta] ){      
      tmpBr = iSigmaPhiSigmaEta;										       
    }													       
  }													       
  
  // Interpolation
  Double_t tmpInter = 1;
  // In eta cracks/gaps 
  if (tmpEta == -1 && tmpBr != -1 ) { // need to interpolate only eta, if br is out of range skip this
    
    if (TMath::Abs(eta) >=rightEta[nBinsEta-1] ) { // out of ECAL boundary
      //      for (Int_t i = 0; i < nBinsEta; ++i) delete h_corr[i];    
      return 1; // don't correct
    }
    
    // central bin at eta = 0
    if (0 <= TMath::Abs(eta) && TMath::Abs(eta) <leftEta[0] ) {
      
      tmpInter = h_corr[0]->GetBinContent(2*tmpBr+1);
      
    }
    else { // all other crack-bins
      
      for (Int_t iEta = 0; iEta < nBinsEta-1; ++iEta){								       
	if (rightEta[iEta] <= TMath::Abs(eta) && TMath::Abs(eta) <leftEta[iEta+1] ){
	  tmpInter = ( h_corr[iEta]  ->GetBinContent(2*tmpBr+1) + 
		       h_corr[iEta+1]->GetBinContent(2*tmpBr+1) ) / 2. ;
	}
      }
    }
    //    for (Int_t i = 0; i < nBinsEta; ++i) delete h_corr[i];    
    return tmpInter;
  }  
  // end interpolation
  
  if (tmpEta == -1 || tmpBr == -1){									       
    //    for (Int_t i = 0; i < nBinsEta; ++i) delete h_corr[i];                                  		       
    return  1; // don't correct												       
  }													       
  
  Double_t tmp = h_corr[tmpEta]->GetBinContent(2*tmpBr+1);				   		       
  //  for (Int_t i = 0; i < nBinsEta; ++i) delete h_corr[i];                                  		       
   return  tmp;	                                                               
  
};													   
                                                                         
Double_t PhoEnergyCorrection::applyScCorrectionsET_EB(Double_t ET){  							   
 
  Int_t iET      = -1;
  for (Int_t iiET  = 0; iiET  < nBinsET;  ++iiET ){ 
    if ( leftET [iiET]  <= (ET)      &&       (ET) < rightET[iiET] ) {
      iET  = iiET;  
    }
  }	
  if (iET == -1) {return 1;}

  Int_t binET =  2*iET+1 ; // h_CBET_EB->FindBin(ET);                                           
  Double_t tmp = h_CBET_EB->GetBinContent(binET);                                   
  if ( 0< binET && binET < 2*nBinsET+1) return tmp;                          
  else return 1.;                                                               
};                                                                               
                                                                                  
                                                                         
Double_t PhoEnergyCorrection::applyScCorrectionsET_EE(Double_t ET){  							   
 

  Int_t iET      = -1;
  for (Int_t iiET  = 0; iiET  < nBinsET;  ++iiET ){ 
    if ( leftET [iiET]  <= (ET)      &&       (ET) < rightET[iiET] ) {
      iET  = iiET;  
      //      cout << leftET [iiET]  << " " << ET  << " " << rightET[iiET] << " " << iET << endl;
    }
  }	
  if (iET == -1) {return 1;}


  Int_t binET =  2*iET+1; // h_CBET_EE->FindBin(ET);                                           
  Double_t tmp = h_CBET_EE->GetBinContent(binET);                                   
  //  cout << tmp  << endl;
  if ( 0< binET && binET < 2*nBinsET+1) return tmp;                          
  else return 1.;                                                               
};                                                                              
*/
